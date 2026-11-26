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


#include "SimdThreeCenterElectronRepulsionVrrRecSLH.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_slh_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t skh0,
                                                          const size_t skg, const size_t skh1,
                                                          const size_t slf0, const size_t slf1,
                                                          const size_t slg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
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
    const auto f_12 = 3.5 / q;
    const auto f_13 = 3.0 / q;

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

    const auto *skh0_0 = buffer.data(skh0 + 0);
    const auto *skh0_3 = buffer.data(skh0 + 3);
    const auto *skh0_5 = buffer.data(skh0 + 5);
    const auto *skh0_6 = buffer.data(skh0 + 6);
    const auto *skh0_9 = buffer.data(skh0 + 9);
    const auto *skh0_15 = buffer.data(skh0 + 15);
    const auto *skh0_20 = buffer.data(skh0 + 20);
    const auto *skh0_24 = buffer.data(skh0 + 24);
    const auto *skh0_27 = buffer.data(skh0 + 27);
    const auto *skh0_36 = buffer.data(skh0 + 36);
    const auto *skh0_42 = buffer.data(skh0 + 42);
    const auto *skh0_47 = buffer.data(skh0 + 47);
    const auto *skh0_51 = buffer.data(skh0 + 51);
    const auto *skh0_62 = buffer.data(skh0 + 62);

    const auto *skg_0 = buffer.data(skg + 0);
    const auto *skg_1 = buffer.data(skg + 1);
    const auto *skg_2 = buffer.data(skg + 2);
    const auto *skg_3 = buffer.data(skg + 3);
    const auto *skg_5 = buffer.data(skg + 5);
    const auto *skg_6 = buffer.data(skg + 6);
    const auto *skg_9 = buffer.data(skg + 9);
    const auto *skg_10 = buffer.data(skg + 10);
    const auto *skg_11 = buffer.data(skg + 11);
    const auto *skg_12 = buffer.data(skg + 12);
    const auto *skg_13 = buffer.data(skg + 13);
    const auto *skg_14 = buffer.data(skg + 14);
    const auto *skg_15 = buffer.data(skg + 15);
    const auto *skg_17 = buffer.data(skg + 17);
    const auto *skg_18 = buffer.data(skg + 18);
    const auto *skg_20 = buffer.data(skg + 20);
    const auto *skg_25 = buffer.data(skg + 25);
    const auto *skg_26 = buffer.data(skg + 26);
    const auto *skg_27 = buffer.data(skg + 27);
    const auto *skg_28 = buffer.data(skg + 28);
    const auto *skg_29 = buffer.data(skg + 29);
    const auto *skg_30 = buffer.data(skg + 30);
    const auto *skg_32 = buffer.data(skg + 32);
    const auto *skg_33 = buffer.data(skg + 33);
    const auto *skg_35 = buffer.data(skg + 35);
    const auto *skg_40 = buffer.data(skg + 40);
    const auto *skg_41 = buffer.data(skg + 41);
    const auto *skg_42 = buffer.data(skg + 42);
    const auto *skg_43 = buffer.data(skg + 43);
    const auto *skg_44 = buffer.data(skg + 44);
    const auto *skg_45 = buffer.data(skg + 45);
    const auto *skg_48 = buffer.data(skg + 48);
    const auto *skg_50 = buffer.data(skg + 50);
    const auto *skg_51 = buffer.data(skg + 51);
    const auto *skg_54 = buffer.data(skg + 54);
    const auto *skg_55 = buffer.data(skg + 55);
    const auto *skg_56 = buffer.data(skg + 56);
    const auto *skg_57 = buffer.data(skg + 57);
    const auto *skg_58 = buffer.data(skg + 58);
    const auto *skg_59 = buffer.data(skg + 59);
    const auto *skg_70 = buffer.data(skg + 70);
    const auto *skg_71 = buffer.data(skg + 71);
    const auto *skg_72 = buffer.data(skg + 72);
    const auto *skg_73 = buffer.data(skg + 73);
    const auto *skg_74 = buffer.data(skg + 74);
    const auto *skg_75 = buffer.data(skg + 75);
    const auto *skg_78 = buffer.data(skg + 78);
    const auto *skg_80 = buffer.data(skg + 80);
    const auto *skg_81 = buffer.data(skg + 81);
    const auto *skg_84 = buffer.data(skg + 84);
    const auto *skg_85 = buffer.data(skg + 85);
    const auto *skg_86 = buffer.data(skg + 86);
    const auto *skg_87 = buffer.data(skg + 87);
    const auto *skg_88 = buffer.data(skg + 88);
    const auto *skg_89 = buffer.data(skg + 89);

    const auto *skh1_0 = buffer.data(skh1 + 0);
    const auto *skh1_3 = buffer.data(skh1 + 3);
    const auto *skh1_5 = buffer.data(skh1 + 5);
    const auto *skh1_6 = buffer.data(skh1 + 6);
    const auto *skh1_9 = buffer.data(skh1 + 9);
    const auto *skh1_15 = buffer.data(skh1 + 15);
    const auto *skh1_20 = buffer.data(skh1 + 20);
    const auto *skh1_24 = buffer.data(skh1 + 24);
    const auto *skh1_27 = buffer.data(skh1 + 27);
    const auto *skh1_36 = buffer.data(skh1 + 36);
    const auto *skh1_42 = buffer.data(skh1 + 42);
    const auto *skh1_47 = buffer.data(skh1 + 47);
    const auto *skh1_51 = buffer.data(skh1 + 51);
    const auto *skh1_62 = buffer.data(skh1 + 62);

    const auto *slf0_0 = buffer.data(slf0 + 0);
    const auto *slf0_3 = buffer.data(slf0 + 3);
    const auto *slf0_5 = buffer.data(slf0 + 5);
    const auto *slf0_6 = buffer.data(slf0 + 6);
    const auto *slf0_8 = buffer.data(slf0 + 8);
    const auto *slf0_9 = buffer.data(slf0 + 9);
    const auto *slf0_16 = buffer.data(slf0 + 16);
    const auto *slf0_18 = buffer.data(slf0 + 18);
    const auto *slf0_19 = buffer.data(slf0 + 19);
    const auto *slf0_28 = buffer.data(slf0 + 28);
    const auto *slf0_29 = buffer.data(slf0 + 29);
    const auto *slf0_30 = buffer.data(slf0 + 30);
    const auto *slf0_33 = buffer.data(slf0 + 33);
    const auto *slf0_35 = buffer.data(slf0 + 35);
    const auto *slf0_36 = buffer.data(slf0 + 36);
    const auto *slf0_38 = buffer.data(slf0 + 38);
    const auto *slf0_39 = buffer.data(slf0 + 39);
    const auto *slf0_48 = buffer.data(slf0 + 48);
    const auto *slf0_49 = buffer.data(slf0 + 49);
    const auto *slf0_50 = buffer.data(slf0 + 50);
    const auto *slf0_53 = buffer.data(slf0 + 53);
    const auto *slf0_55 = buffer.data(slf0 + 55);
    const auto *slf0_56 = buffer.data(slf0 + 56);
    const auto *slf0_58 = buffer.data(slf0 + 58);
    const auto *slf0_59 = buffer.data(slf0 + 59);

    const auto *slf1_0 = buffer.data(slf1 + 0);
    const auto *slf1_3 = buffer.data(slf1 + 3);
    const auto *slf1_5 = buffer.data(slf1 + 5);
    const auto *slf1_6 = buffer.data(slf1 + 6);
    const auto *slf1_8 = buffer.data(slf1 + 8);
    const auto *slf1_9 = buffer.data(slf1 + 9);
    const auto *slf1_16 = buffer.data(slf1 + 16);
    const auto *slf1_18 = buffer.data(slf1 + 18);
    const auto *slf1_19 = buffer.data(slf1 + 19);
    const auto *slf1_28 = buffer.data(slf1 + 28);
    const auto *slf1_29 = buffer.data(slf1 + 29);
    const auto *slf1_30 = buffer.data(slf1 + 30);
    const auto *slf1_33 = buffer.data(slf1 + 33);
    const auto *slf1_35 = buffer.data(slf1 + 35);
    const auto *slf1_36 = buffer.data(slf1 + 36);
    const auto *slf1_38 = buffer.data(slf1 + 38);
    const auto *slf1_39 = buffer.data(slf1 + 39);
    const auto *slf1_48 = buffer.data(slf1 + 48);
    const auto *slf1_49 = buffer.data(slf1 + 49);
    const auto *slf1_50 = buffer.data(slf1 + 50);
    const auto *slf1_53 = buffer.data(slf1 + 53);
    const auto *slf1_55 = buffer.data(slf1 + 55);
    const auto *slf1_56 = buffer.data(slf1 + 56);
    const auto *slf1_58 = buffer.data(slf1 + 58);
    const auto *slf1_59 = buffer.data(slf1 + 59);

    const auto *slg_0 = buffer.data(slg + 0);
    const auto *slg_2 = buffer.data(slg + 2);
    const auto *slg_3 = buffer.data(slg + 3);
    const auto *slg_5 = buffer.data(slg + 5);
    const auto *slg_6 = buffer.data(slg + 6);
    const auto *slg_9 = buffer.data(slg + 9);
    const auto *slg_10 = buffer.data(slg + 10);
    const auto *slg_11 = buffer.data(slg + 11);
    const auto *slg_12 = buffer.data(slg + 12);
    const auto *slg_13 = buffer.data(slg + 13);
    const auto *slg_14 = buffer.data(slg + 14);
    const auto *slg_15 = buffer.data(slg + 15);
    const auto *slg_17 = buffer.data(slg + 17);
    const auto *slg_18 = buffer.data(slg + 18);
    const auto *slg_20 = buffer.data(slg + 20);
    const auto *slg_25 = buffer.data(slg + 25);
    const auto *slg_26 = buffer.data(slg + 26);
    const auto *slg_27 = buffer.data(slg + 27);
    const auto *slg_28 = buffer.data(slg + 28);
    const auto *slg_29 = buffer.data(slg + 29);
    const auto *slg_30 = buffer.data(slg + 30);
    const auto *slg_32 = buffer.data(slg + 32);
    const auto *slg_33 = buffer.data(slg + 33);
    const auto *slg_35 = buffer.data(slg + 35);
    const auto *slg_40 = buffer.data(slg + 40);
    const auto *slg_41 = buffer.data(slg + 41);
    const auto *slg_42 = buffer.data(slg + 42);
    const auto *slg_43 = buffer.data(slg + 43);
    const auto *slg_44 = buffer.data(slg + 44);
    const auto *slg_45 = buffer.data(slg + 45);
    const auto *slg_47 = buffer.data(slg + 47);
    const auto *slg_48 = buffer.data(slg + 48);
    const auto *slg_50 = buffer.data(slg + 50);
    const auto *slg_51 = buffer.data(slg + 51);
    const auto *slg_54 = buffer.data(slg + 54);
    const auto *slg_55 = buffer.data(slg + 55);
    const auto *slg_56 = buffer.data(slg + 56);
    const auto *slg_57 = buffer.data(slg + 57);
    const auto *slg_58 = buffer.data(slg + 58);
    const auto *slg_59 = buffer.data(slg + 59);
    const auto *slg_60 = buffer.data(slg + 60);
    const auto *slg_62 = buffer.data(slg + 62);
    const auto *slg_63 = buffer.data(slg + 63);
    const auto *slg_65 = buffer.data(slg + 65);
    const auto *slg_70 = buffer.data(slg + 70);
    const auto *slg_71 = buffer.data(slg + 71);
    const auto *slg_72 = buffer.data(slg + 72);
    const auto *slg_73 = buffer.data(slg + 73);
    const auto *slg_74 = buffer.data(slg + 74);
    const auto *slg_75 = buffer.data(slg + 75);
    const auto *slg_77 = buffer.data(slg + 77);
    const auto *slg_78 = buffer.data(slg + 78);
    const auto *slg_80 = buffer.data(slg + 80);
    const auto *slg_81 = buffer.data(slg + 81);
    const auto *slg_84 = buffer.data(slg + 84);
    const auto *slg_85 = buffer.data(slg + 85);
    const auto *slg_86 = buffer.data(slg + 86);
    const auto *slg_87 = buffer.data(slg + 87);
    const auto *slg_88 = buffer.data(slg + 88);
    const auto *slg_89 = buffer.data(slg + 89);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, skg_0, skg_3, slf0_0, slf0_3, \
                         slf1_0, slf1_3, slg_0, slg_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * skg_0[k]
                 + f_1 * slf0_0[k]
                 - f_2 * slf1_0[k]
                 + f_3 * pc_x[k] * slg_0[k];

        t_1[k] = f_3 * pc_y[k] * slg_0[k];

        t_2[k] = f_3 * pc_z[k] * slg_0[k];

        t_3[k] = f_0 * skg_3[k]
                 + f_4 * slf0_3[k]
                 - f_5 * slf1_3[k]
                 + f_3 * pc_x[k] * slg_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pc_x, pc_y, skg_5, skg_6, slf0_5, slf0_6, slf1_5, \
                         slf1_6, slg_2, slg_5, slg_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * pc_y[k] * slg_2[k];

        t_5[k] = f_0 * skg_5[k]
                 + f_4 * slf0_5[k]
                 - f_5 * slf1_5[k]
                 + f_3 * pc_x[k] * slg_5[k];

        t_6[k] = f_0 * skg_6[k]
                 + f_6 * slf0_6[k]
                 - f_7 * slf1_6[k]
                 + f_3 * pc_x[k] * slg_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, pc_x, pc_y, pc_z, skg_9, skg_10, slf0_9, slf1_9, \
                         slg_3, slg_5, slg_9, slg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_3 * pc_z[k] * slg_3[k];

        t_8[k] = f_3 * pc_y[k] * slg_5[k];

        t_9[k] = f_0 * skg_9[k]
                 + f_6 * slf0_9[k]
                 - f_7 * slf1_9[k]
                 + f_3 * pc_x[k] * slg_9[k];

        t_10[k] = f_0 * skg_10[k]
                  + f_3 * pc_x[k] * slg_10[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, pc_x, skg_11, skg_12, skg_13, skg_14, slg_11, \
                         slg_12, slg_13, slg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_0 * skg_11[k]
                  + f_3 * pc_x[k] * slg_11[k];

        t_12[k] = f_0 * skg_12[k]
                  + f_3 * pc_x[k] * slg_12[k];

        t_13[k] = f_0 * skg_13[k]
                  + f_3 * pc_x[k] * slg_13[k];

        t_14[k] = f_0 * skg_14[k]
                  + f_3 * pc_x[k] * slg_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pc_y, pc_z, slf0_6, slf0_8, slf0_9, slf1_6, \
                         slf1_8, slf1_9, slg_10, slg_12, slg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_1 * slf0_6[k]
                  - f_2 * slf1_6[k]
                  + f_3 * pc_y[k] * slg_10[k];

        t_16[k] = f_3 * pc_z[k] * slg_10[k];

        t_17[k] = f_4 * slf0_8[k]
                  - f_5 * slf1_8[k]
                  + f_3 * pc_y[k] * slg_12[k];

        t_18[k] = f_6 * slf0_9[k]
                  - f_7 * slf1_9[k]
                  + f_3 * pc_y[k] * slg_13[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pb_y, pc_y, pc_z, skh0_0, skg_0, \
                         skh1_0, slf0_9, slf1_9, slg_14, slg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * pc_y[k] * slg_14[k];

        t_20[k] = f_1 * slf0_9[k]
                  - f_2 * slf1_9[k]
                  + f_3 * pc_z[k] * slg_14[k];

        t_21[k] = pb_y[k] * skh0_0[k]
                  - f_8 * pc_y[k] * skh1_0[k];

        t_22[k] = f_9 * skg_0[k]
                  + f_3 * pc_y[k] * slg_15[k];

        t_23[k] = f_3 * pc_z[k] * slg_15[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pb_y, pc_y, skh0_3, skh0_5, skh0_6, skg_1, \
                         skg_2, skg_3, skh1_3, skh1_5, skh1_6, slg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pb_y[k] * skh0_3[k]
                  + f_10 * skg_1[k]
                  - f_8 * pc_y[k] * skh1_3[k];

        t_25[k] = f_9 * skg_2[k]
                  + f_3 * pc_y[k] * slg_17[k];

        t_26[k] = pb_y[k] * skh0_5[k]
                  - f_8 * pc_y[k] * skh1_5[k];

        t_27[k] = pb_y[k] * skh0_6[k]
                  + f_11 * skg_3[k]
                  - f_8 * pc_y[k] * skh1_6[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pb_y, pc_x, pc_y, pc_z, skh0_9, skg_5, \
                         skg_25, skh1_9, slg_18, slg_20, slg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_3 * pc_z[k] * slg_18[k];

        t_29[k] = f_9 * skg_5[k]
                  + f_3 * pc_y[k] * slg_20[k];

        t_30[k] = pb_y[k] * skh0_9[k]
                  - f_8 * pc_y[k] * skh1_9[k];

        t_31[k] = f_12 * skg_25[k]
                  + f_3 * pc_x[k] * slg_25[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pc_x, skg_26, skg_27, skg_28, skg_29, slg_26, \
                         slg_27, slg_28, slg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_12 * skg_26[k]
                  + f_3 * pc_x[k] * slg_26[k];

        t_33[k] = f_12 * skg_27[k]
                  + f_3 * pc_x[k] * slg_27[k];

        t_34[k] = f_12 * skg_28[k]
                  + f_3 * pc_x[k] * slg_28[k];

        t_35[k] = f_12 * skg_29[k]
                  + f_3 * pc_x[k] * slg_29[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pc_y, pc_z, skg_10, skg_12, slf0_16, slf0_18, \
                         slf1_16, slf1_18, slg_25, slg_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_9 * skg_10[k]
                  + f_1 * slf0_16[k]
                  - f_2 * slf1_16[k]
                  + f_3 * pc_y[k] * slg_25[k];

        t_37[k] = f_3 * pc_z[k] * slg_25[k];

        t_38[k] = f_9 * skg_12[k]
                  + f_4 * slf0_18[k]
                  - f_5 * slf1_18[k]
                  + f_3 * pc_y[k] * slg_27[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pb_y, pc_y, skh0_20, skg_13, skg_14, skh1_20, \
                         slf0_19, slf1_19, slg_28, slg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_9 * skg_13[k]
                  + f_6 * slf0_19[k]
                  - f_7 * slf1_19[k]
                  + f_3 * pc_y[k] * slg_28[k];

        t_40[k] = f_9 * skg_14[k]
                  + f_3 * pc_y[k] * slg_29[k];

        t_41[k] = pb_y[k] * skh0_20[k]
                  - f_8 * pc_y[k] * skh1_20[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, pb_z, pc_y, pc_z, skh0_0, skh0_3, \
                         skg_0, skh1_0, skh1_3, slg_30, slg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_z[k] * skh0_0[k]
                  - f_8 * pc_z[k] * skh1_0[k];

        t_43[k] = f_3 * pc_y[k] * slg_30[k];

        t_44[k] = f_9 * skg_0[k]
                  + f_3 * pc_z[k] * slg_30[k];

        t_45[k] = pb_z[k] * skh0_3[k]
                  - f_8 * pc_z[k] * skh1_3[k];

        t_46[k] = f_3 * pc_y[k] * slg_32[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pb_z, pc_y, pc_z, skh0_5, skh0_6, skg_2, \
                         skg_3, skh1_5, skh1_6, slg_33, slg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = pb_z[k] * skh0_5[k]
                  + f_10 * skg_2[k]
                  - f_8 * pc_z[k] * skh1_5[k];

        t_48[k] = pb_z[k] * skh0_6[k]
                  - f_8 * pc_z[k] * skh1_6[k];

        t_49[k] = f_9 * skg_3[k]
                  + f_3 * pc_z[k] * slg_33[k];

        t_50[k] = f_3 * pc_y[k] * slg_35[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pb_z, pc_x, pc_z, skh0_9, skg_5, skg_40, \
                         skg_41, skg_42, skh1_9, slg_40, slg_41, \
                         slg_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = pb_z[k] * skh0_9[k]
                  + f_11 * skg_5[k]
                  - f_8 * pc_z[k] * skh1_9[k];

        t_52[k] = f_12 * skg_40[k]
                  + f_3 * pc_x[k] * slg_40[k];

        t_53[k] = f_12 * skg_41[k]
                  + f_3 * pc_x[k] * slg_41[k];

        t_54[k] = f_12 * skg_42[k]
                  + f_3 * pc_x[k] * slg_42[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pb_z, pc_x, pc_z, skh0_15, skg_10, skg_43, \
                         skg_44, skh1_15, slg_40, slg_43, slg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_12 * skg_43[k]
                  + f_3 * pc_x[k] * slg_43[k];

        t_56[k] = f_12 * skg_44[k]
                  + f_3 * pc_x[k] * slg_44[k];

        t_57[k] = pb_z[k] * skh0_15[k]
                  - f_8 * pc_z[k] * skh1_15[k];

        t_58[k] = f_9 * skg_10[k]
                  + f_3 * pc_z[k] * slg_40[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pc_y, pc_z, skg_14, slf0_28, slf0_29, \
                         slf1_28, slf1_29, slg_42, slg_43, slg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_4 * slf0_28[k]
                  - f_5 * slf1_28[k]
                  + f_3 * pc_y[k] * slg_42[k];

        t_60[k] = f_6 * slf0_29[k]
                  - f_7 * slf1_29[k]
                  + f_3 * pc_y[k] * slg_43[k];

        t_61[k] = f_3 * pc_y[k] * slg_44[k];

        t_62[k] = f_9 * skg_14[k]
                  + f_1 * slf0_29[k]
                  - f_2 * slf1_29[k]
                  + f_3 * pc_z[k] * slg_44[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pc_x, pc_y, pc_z, skg_15, skg_45, skg_48, \
                         slf0_30, slf0_33, slf1_30, slf1_33, slg_45, \
                         slg_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_13 * skg_45[k]
                  + f_1 * slf0_30[k]
                  - f_2 * slf1_30[k]
                  + f_3 * pc_x[k] * slg_45[k];

        t_64[k] = f_10 * skg_15[k]
                  + f_3 * pc_y[k] * slg_45[k];

        t_65[k] = f_3 * pc_z[k] * slg_45[k];

        t_66[k] = f_13 * skg_48[k]
                  + f_4 * slf0_33[k]
                  - f_5 * slf1_33[k]
                  + f_3 * pc_x[k] * slg_48[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, pc_x, pc_y, skg_17, skg_50, skg_51, slf0_35, \
                         slf0_36, slf1_35, slf1_36, slg_47, slg_50, \
                         slg_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_10 * skg_17[k]
                  + f_3 * pc_y[k] * slg_47[k];

        t_68[k] = f_13 * skg_50[k]
                  + f_4 * slf0_35[k]
                  - f_5 * slf1_35[k]
                  + f_3 * pc_x[k] * slg_50[k];

        t_69[k] = f_13 * skg_51[k]
                  + f_6 * slf0_36[k]
                  - f_7 * slf1_36[k]
                  + f_3 * pc_x[k] * slg_51[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pc_x, pc_y, pc_z, skg_20, skg_54, skg_55, \
                         slf0_39, slf1_39, slg_48, slg_50, slg_54, \
                         slg_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_3 * pc_z[k] * slg_48[k];

        t_71[k] = f_10 * skg_20[k]
                  + f_3 * pc_y[k] * slg_50[k];

        t_72[k] = f_13 * skg_54[k]
                  + f_6 * slf0_39[k]
                  - f_7 * slf1_39[k]
                  + f_3 * pc_x[k] * slg_54[k];

        t_73[k] = f_13 * skg_55[k]
                  + f_3 * pc_x[k] * slg_55[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pc_x, skg_56, skg_57, skg_58, skg_59, slg_56, \
                         slg_57, slg_58, slg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_13 * skg_56[k]
                  + f_3 * pc_x[k] * slg_56[k];

        t_75[k] = f_13 * skg_57[k]
                  + f_3 * pc_x[k] * slg_57[k];

        t_76[k] = f_13 * skg_58[k]
                  + f_3 * pc_x[k] * slg_58[k];

        t_77[k] = f_13 * skg_59[k]
                  + f_3 * pc_x[k] * slg_59[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pc_y, pc_z, skg_25, skg_27, slf0_36, slf0_38, \
                         slf1_36, slf1_38, slg_55, slg_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_10 * skg_25[k]
                  + f_1 * slf0_36[k]
                  - f_2 * slf1_36[k]
                  + f_3 * pc_y[k] * slg_55[k];

        t_79[k] = f_3 * pc_z[k] * slg_55[k];

        t_80[k] = f_10 * skg_27[k]
                  + f_4 * slf0_38[k]
                  - f_5 * slf1_38[k]
                  + f_3 * pc_y[k] * slg_57[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pb_y, pc_y, pc_z, skh0_42, skg_28, skg_29, \
                         skh1_42, slf0_39, slf1_39, slg_58, slg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_10 * skg_28[k]
                  + f_6 * slf0_39[k]
                  - f_7 * slf1_39[k]
                  + f_3 * pc_y[k] * slg_58[k];

        t_82[k] = f_10 * skg_29[k]
                  + f_3 * pc_y[k] * slg_59[k];

        t_83[k] = f_1 * slf0_39[k]
                  - f_2 * slf1_39[k]
                  + f_3 * pc_z[k] * slg_59[k];

        t_84[k] = pb_y[k] * skh0_42[k]
                  - f_8 * pc_y[k] * skh1_42[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pb_z, pc_y, pc_z, skh0_24, skg_15, skg_30, \
                         skg_32, skh1_24, slg_60, slg_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_9 * skg_30[k]
                  + f_3 * pc_y[k] * slg_60[k];

        t_86[k] = f_9 * skg_15[k]
                  + f_3 * pc_z[k] * slg_60[k];

        t_87[k] = pb_z[k] * skh0_24[k]
                  - f_8 * pc_z[k] * skh1_24[k];

        t_88[k] = f_9 * skg_32[k]
                  + f_3 * pc_y[k] * slg_62[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, pb_y, pb_z, pc_y, pc_z, skh0_27, skh0_47, \
                         skg_18, skg_35, skh1_27, skh1_47, slg_63, \
                         slg_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = pb_y[k] * skh0_47[k]
                  - f_8 * pc_y[k] * skh1_47[k];

        t_90[k] = pb_z[k] * skh0_27[k]
                  - f_8 * pc_z[k] * skh1_27[k];

        t_91[k] = f_9 * skg_18[k]
                  + f_3 * pc_z[k] * slg_63[k];

        t_92[k] = f_9 * skg_35[k]
                  + f_3 * pc_y[k] * slg_65[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, pb_y, pc_x, pc_y, skh0_51, skg_70, skg_71, \
                         skg_72, skh1_51, slg_70, slg_71, slg_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = pb_y[k] * skh0_51[k]
                  - f_8 * pc_y[k] * skh1_51[k];

        t_94[k] = f_13 * skg_70[k]
                  + f_3 * pc_x[k] * slg_70[k];

        t_95[k] = f_13 * skg_71[k]
                  + f_3 * pc_x[k] * slg_71[k];

        t_96[k] = f_13 * skg_72[k]
                  + f_3 * pc_x[k] * slg_72[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pb_z, pc_x, pc_z, skh0_36, skg_25, skg_73, \
                         skg_74, skh1_36, slg_70, slg_73, slg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_13 * skg_73[k]
                  + f_3 * pc_x[k] * slg_73[k];

        t_98[k] = f_13 * skg_74[k]
                  + f_3 * pc_x[k] * slg_74[k];

        t_99[k] = pb_z[k] * skh0_36[k]
                  - f_8 * pc_z[k] * skh1_36[k];

        t_100[k] = f_9 * skg_25[k]
                   + f_3 * pc_z[k] * slg_70[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, pc_y, skg_42, skg_43, skg_44, slf0_48, slf0_49, \
                         slf1_48, slf1_49, slg_72, slg_73, slg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_9 * skg_42[k]
                   + f_4 * slf0_48[k]
                   - f_5 * slf1_48[k]
                   + f_3 * pc_y[k] * slg_72[k];

        t_102[k] = f_9 * skg_43[k]
                   + f_6 * slf0_49[k]
                   - f_7 * slf1_49[k]
                   + f_3 * pc_y[k] * slg_73[k];

        t_103[k] = f_9 * skg_44[k]
                   + f_3 * pc_y[k] * slg_74[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pb_y, pc_x, pc_y, pc_z, skh0_62, skg_30, \
                         skg_75, skh1_62, slf0_50, slf1_50, slg_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pb_y[k] * skh0_62[k]
                   - f_8 * pc_y[k] * skh1_62[k];

        t_105[k] = f_13 * skg_75[k]
                   + f_1 * slf0_50[k]
                   - f_2 * slf1_50[k]
                   + f_3 * pc_x[k] * slg_75[k];

        t_106[k] = f_3 * pc_y[k] * slg_75[k];

        t_107[k] = f_10 * skg_30[k]
                   + f_3 * pc_z[k] * slg_75[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, pc_x, pc_y, skg_78, skg_80, slf0_53, slf0_55, \
                         slf1_53, slf1_55, slg_77, slg_78, slg_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_13 * skg_78[k]
                   + f_4 * slf0_53[k]
                   - f_5 * slf1_53[k]
                   + f_3 * pc_x[k] * slg_78[k];

        t_109[k] = f_3 * pc_y[k] * slg_77[k];

        t_110[k] = f_13 * skg_80[k]
                   + f_4 * slf0_55[k]
                   - f_5 * slf1_55[k]
                   + f_3 * pc_x[k] * slg_80[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, pc_x, pc_y, pc_z, skg_33, skg_81, slf0_56, \
                         slf1_56, slg_78, slg_80, slg_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_13 * skg_81[k]
                   + f_6 * slf0_56[k]
                   - f_7 * slf1_56[k]
                   + f_3 * pc_x[k] * slg_81[k];

        t_112[k] = f_10 * skg_33[k]
                   + f_3 * pc_z[k] * slg_78[k];

        t_113[k] = f_3 * pc_y[k] * slg_80[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pc_x, skg_84, skg_85, skg_86, skg_87, \
                         slf0_59, slf1_59, slg_84, slg_85, slg_86, \
                         slg_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_13 * skg_84[k]
                   + f_6 * slf0_59[k]
                   - f_7 * slf1_59[k]
                   + f_3 * pc_x[k] * slg_84[k];

        t_115[k] = f_13 * skg_85[k]
                   + f_3 * pc_x[k] * slg_85[k];

        t_116[k] = f_13 * skg_86[k]
                   + f_3 * pc_x[k] * slg_86[k];

        t_117[k] = f_13 * skg_87[k]
                   + f_3 * pc_x[k] * slg_87[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, pc_x, pc_y, pc_z, skg_40, skg_88, skg_89, \
                         slf0_56, slf1_56, slg_85, slg_88, slg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_13 * skg_88[k]
                   + f_3 * pc_x[k] * slg_88[k];

        t_119[k] = f_13 * skg_89[k]
                   + f_3 * pc_x[k] * slg_89[k];

        t_120[k] = f_1 * slf0_56[k]
                   - f_2 * slf1_56[k]
                   + f_3 * pc_y[k] * slg_85[k];

        t_121[k] = f_10 * skg_40[k]
                   + f_3 * pc_z[k] * slg_85[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, pc_y, pc_z, skg_44, slf0_58, slf0_59, \
                         slf1_58, slf1_59, slg_87, slg_88, slg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_4 * slf0_58[k]
                   - f_5 * slf1_58[k]
                   + f_3 * pc_y[k] * slg_87[k];

        t_123[k] = f_6 * slf0_59[k]
                   - f_7 * slf1_59[k]
                   + f_3 * pc_y[k] * slg_88[k];

        t_124[k] = f_3 * pc_y[k] * slg_89[k];

        t_125[k] = f_10 * skg_44[k]
                   + f_1 * slf0_59[k]
                   - f_2 * slf1_59[k]
                   + f_3 * pc_z[k] * slg_89[k];
    }
}

static auto
compute_prim_slh_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t skh0,
                                                          const size_t skg, const size_t skh1,
                                                          const size_t slf0, const size_t slf1,
                                                          const size_t slg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_14 = 2.5 / q;
    const auto f_15 = 2.0 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *skh0_63 = buffer.data(skh0 + 63);
    const auto *skh0_66 = buffer.data(skh0 + 66);
    const auto *skh0_69 = buffer.data(skh0 + 69);
    const auto *skh0_78 = buffer.data(skh0 + 78);
    const auto *skh0_105 = buffer.data(skh0 + 105);
    const auto *skh0_108 = buffer.data(skh0 + 108);
    const auto *skh0_110 = buffer.data(skh0 + 110);
    const auto *skh0_111 = buffer.data(skh0 + 111);
    const auto *skh0_114 = buffer.data(skh0 + 114);
    const auto *skh0_125 = buffer.data(skh0 + 125);
    const auto *skh0_126 = buffer.data(skh0 + 126);
    const auto *skh0_129 = buffer.data(skh0 + 129);
    const auto *skh0_132 = buffer.data(skh0 + 132);

    const auto *skg_45 = buffer.data(skg + 45);
    const auto *skg_47 = buffer.data(skg + 47);
    const auto *skg_48 = buffer.data(skg + 48);
    const auto *skg_50 = buffer.data(skg + 50);
    const auto *skg_55 = buffer.data(skg + 55);
    const auto *skg_57 = buffer.data(skg + 57);
    const auto *skg_58 = buffer.data(skg + 58);
    const auto *skg_59 = buffer.data(skg + 59);
    const auto *skg_60 = buffer.data(skg + 60);
    const auto *skg_62 = buffer.data(skg + 62);
    const auto *skg_63 = buffer.data(skg + 63);
    const auto *skg_65 = buffer.data(skg + 65);
    const auto *skg_70 = buffer.data(skg + 70);
    const auto *skg_72 = buffer.data(skg + 72);
    const auto *skg_73 = buffer.data(skg + 73);
    const auto *skg_74 = buffer.data(skg + 74);
    const auto *skg_75 = buffer.data(skg + 75);
    const auto *skg_76 = buffer.data(skg + 76);
    const auto *skg_77 = buffer.data(skg + 77);
    const auto *skg_78 = buffer.data(skg + 78);
    const auto *skg_80 = buffer.data(skg + 80);
    const auto *skg_85 = buffer.data(skg + 85);
    const auto *skg_87 = buffer.data(skg + 87);
    const auto *skg_88 = buffer.data(skg + 88);
    const auto *skg_89 = buffer.data(skg + 89);
    const auto *skg_90 = buffer.data(skg + 90);
    const auto *skg_92 = buffer.data(skg + 92);
    const auto *skg_93 = buffer.data(skg + 93);
    const auto *skg_95 = buffer.data(skg + 95);
    const auto *skg_96 = buffer.data(skg + 96);
    const auto *skg_99 = buffer.data(skg + 99);
    const auto *skg_100 = buffer.data(skg + 100);
    const auto *skg_101 = buffer.data(skg + 101);
    const auto *skg_102 = buffer.data(skg + 102);
    const auto *skg_103 = buffer.data(skg + 103);
    const auto *skg_104 = buffer.data(skg + 104);
    const auto *skg_105 = buffer.data(skg + 105);
    const auto *skg_107 = buffer.data(skg + 107);
    const auto *skg_110 = buffer.data(skg + 110);
    const auto *skg_114 = buffer.data(skg + 114);
    const auto *skg_115 = buffer.data(skg + 115);
    const auto *skg_116 = buffer.data(skg + 116);
    const auto *skg_117 = buffer.data(skg + 117);
    const auto *skg_118 = buffer.data(skg + 118);
    const auto *skg_119 = buffer.data(skg + 119);
    const auto *skg_130 = buffer.data(skg + 130);
    const auto *skg_131 = buffer.data(skg + 131);
    const auto *skg_132 = buffer.data(skg + 132);
    const auto *skg_133 = buffer.data(skg + 133);
    const auto *skg_134 = buffer.data(skg + 134);
    const auto *skg_135 = buffer.data(skg + 135);
    const auto *skg_138 = buffer.data(skg + 138);
    const auto *skg_140 = buffer.data(skg + 140);
    const auto *skg_141 = buffer.data(skg + 141);
    const auto *skg_144 = buffer.data(skg + 144);
    const auto *skg_145 = buffer.data(skg + 145);
    const auto *skg_146 = buffer.data(skg + 146);
    const auto *skg_147 = buffer.data(skg + 147);
    const auto *skg_148 = buffer.data(skg + 148);
    const auto *skg_149 = buffer.data(skg + 149);
    const auto *skg_150 = buffer.data(skg + 150);
    const auto *skg_153 = buffer.data(skg + 153);
    const auto *skg_155 = buffer.data(skg + 155);
    const auto *skg_156 = buffer.data(skg + 156);
    const auto *skg_159 = buffer.data(skg + 159);
    const auto *skg_160 = buffer.data(skg + 160);
    const auto *skg_161 = buffer.data(skg + 161);
    const auto *skg_162 = buffer.data(skg + 162);
    const auto *skg_163 = buffer.data(skg + 163);
    const auto *skg_164 = buffer.data(skg + 164);
    const auto *skg_170 = buffer.data(skg + 170);
    const auto *skg_174 = buffer.data(skg + 174);
    const auto *skg_175 = buffer.data(skg + 175);
    const auto *skg_176 = buffer.data(skg + 176);

    const auto *skh1_63 = buffer.data(skh1 + 63);
    const auto *skh1_66 = buffer.data(skh1 + 66);
    const auto *skh1_69 = buffer.data(skh1 + 69);
    const auto *skh1_78 = buffer.data(skh1 + 78);
    const auto *skh1_105 = buffer.data(skh1 + 105);
    const auto *skh1_108 = buffer.data(skh1 + 108);
    const auto *skh1_110 = buffer.data(skh1 + 110);
    const auto *skh1_111 = buffer.data(skh1 + 111);
    const auto *skh1_114 = buffer.data(skh1 + 114);
    const auto *skh1_125 = buffer.data(skh1 + 125);
    const auto *skh1_126 = buffer.data(skh1 + 126);
    const auto *skh1_129 = buffer.data(skh1 + 129);
    const auto *skh1_132 = buffer.data(skh1 + 132);

    const auto *slf0_60 = buffer.data(slf0 + 60);
    const auto *slf0_63 = buffer.data(slf0 + 63);
    const auto *slf0_65 = buffer.data(slf0 + 65);
    const auto *slf0_66 = buffer.data(slf0 + 66);
    const auto *slf0_68 = buffer.data(slf0 + 68);
    const auto *slf0_69 = buffer.data(slf0 + 69);
    const auto *slf0_75 = buffer.data(slf0 + 75);
    const auto *slf0_78 = buffer.data(slf0 + 78);
    const auto *slf0_79 = buffer.data(slf0 + 79);
    const auto *slf0_86 = buffer.data(slf0 + 86);
    const auto *slf0_88 = buffer.data(slf0 + 88);
    const auto *slf0_89 = buffer.data(slf0 + 89);
    const auto *slf0_90 = buffer.data(slf0 + 90);
    const auto *slf0_93 = buffer.data(slf0 + 93);
    const auto *slf0_95 = buffer.data(slf0 + 95);
    const auto *slf0_96 = buffer.data(slf0 + 96);
    const auto *slf0_98 = buffer.data(slf0 + 98);
    const auto *slf0_99 = buffer.data(slf0 + 99);
    const auto *slf0_100 = buffer.data(slf0 + 100);
    const auto *slf0_103 = buffer.data(slf0 + 103);
    const auto *slf0_105 = buffer.data(slf0 + 105);
    const auto *slf0_106 = buffer.data(slf0 + 106);
    const auto *slf0_108 = buffer.data(slf0 + 108);
    const auto *slf0_109 = buffer.data(slf0 + 109);
    const auto *slf0_115 = buffer.data(slf0 + 115);
    const auto *slf0_119 = buffer.data(slf0 + 119);

    const auto *slf1_60 = buffer.data(slf1 + 60);
    const auto *slf1_63 = buffer.data(slf1 + 63);
    const auto *slf1_65 = buffer.data(slf1 + 65);
    const auto *slf1_66 = buffer.data(slf1 + 66);
    const auto *slf1_68 = buffer.data(slf1 + 68);
    const auto *slf1_69 = buffer.data(slf1 + 69);
    const auto *slf1_75 = buffer.data(slf1 + 75);
    const auto *slf1_78 = buffer.data(slf1 + 78);
    const auto *slf1_79 = buffer.data(slf1 + 79);
    const auto *slf1_86 = buffer.data(slf1 + 86);
    const auto *slf1_88 = buffer.data(slf1 + 88);
    const auto *slf1_89 = buffer.data(slf1 + 89);
    const auto *slf1_90 = buffer.data(slf1 + 90);
    const auto *slf1_93 = buffer.data(slf1 + 93);
    const auto *slf1_95 = buffer.data(slf1 + 95);
    const auto *slf1_96 = buffer.data(slf1 + 96);
    const auto *slf1_98 = buffer.data(slf1 + 98);
    const auto *slf1_99 = buffer.data(slf1 + 99);
    const auto *slf1_100 = buffer.data(slf1 + 100);
    const auto *slf1_103 = buffer.data(slf1 + 103);
    const auto *slf1_105 = buffer.data(slf1 + 105);
    const auto *slf1_106 = buffer.data(slf1 + 106);
    const auto *slf1_108 = buffer.data(slf1 + 108);
    const auto *slf1_109 = buffer.data(slf1 + 109);
    const auto *slf1_115 = buffer.data(slf1 + 115);
    const auto *slf1_119 = buffer.data(slf1 + 119);

    const auto *slg_90 = buffer.data(slg + 90);
    const auto *slg_92 = buffer.data(slg + 92);
    const auto *slg_93 = buffer.data(slg + 93);
    const auto *slg_95 = buffer.data(slg + 95);
    const auto *slg_96 = buffer.data(slg + 96);
    const auto *slg_99 = buffer.data(slg + 99);
    const auto *slg_100 = buffer.data(slg + 100);
    const auto *slg_101 = buffer.data(slg + 101);
    const auto *slg_102 = buffer.data(slg + 102);
    const auto *slg_103 = buffer.data(slg + 103);
    const auto *slg_104 = buffer.data(slg + 104);
    const auto *slg_105 = buffer.data(slg + 105);
    const auto *slg_107 = buffer.data(slg + 107);
    const auto *slg_108 = buffer.data(slg + 108);
    const auto *slg_110 = buffer.data(slg + 110);
    const auto *slg_114 = buffer.data(slg + 114);
    const auto *slg_115 = buffer.data(slg + 115);
    const auto *slg_116 = buffer.data(slg + 116);
    const auto *slg_117 = buffer.data(slg + 117);
    const auto *slg_118 = buffer.data(slg + 118);
    const auto *slg_119 = buffer.data(slg + 119);
    const auto *slg_120 = buffer.data(slg + 120);
    const auto *slg_122 = buffer.data(slg + 122);
    const auto *slg_123 = buffer.data(slg + 123);
    const auto *slg_125 = buffer.data(slg + 125);
    const auto *slg_130 = buffer.data(slg + 130);
    const auto *slg_131 = buffer.data(slg + 131);
    const auto *slg_132 = buffer.data(slg + 132);
    const auto *slg_133 = buffer.data(slg + 133);
    const auto *slg_134 = buffer.data(slg + 134);
    const auto *slg_135 = buffer.data(slg + 135);
    const auto *slg_137 = buffer.data(slg + 137);
    const auto *slg_138 = buffer.data(slg + 138);
    const auto *slg_140 = buffer.data(slg + 140);
    const auto *slg_141 = buffer.data(slg + 141);
    const auto *slg_144 = buffer.data(slg + 144);
    const auto *slg_145 = buffer.data(slg + 145);
    const auto *slg_146 = buffer.data(slg + 146);
    const auto *slg_147 = buffer.data(slg + 147);
    const auto *slg_148 = buffer.data(slg + 148);
    const auto *slg_149 = buffer.data(slg + 149);
    const auto *slg_150 = buffer.data(slg + 150);
    const auto *slg_152 = buffer.data(slg + 152);
    const auto *slg_153 = buffer.data(slg + 153);
    const auto *slg_155 = buffer.data(slg + 155);
    const auto *slg_156 = buffer.data(slg + 156);
    const auto *slg_159 = buffer.data(slg + 159);
    const auto *slg_160 = buffer.data(slg + 160);
    const auto *slg_161 = buffer.data(slg + 161);
    const auto *slg_162 = buffer.data(slg + 162);
    const auto *slg_163 = buffer.data(slg + 163);
    const auto *slg_164 = buffer.data(slg + 164);
    const auto *slg_165 = buffer.data(slg + 165);
    const auto *slg_167 = buffer.data(slg + 167);
    const auto *slg_168 = buffer.data(slg + 168);
    const auto *slg_170 = buffer.data(slg + 170);
    const auto *slg_174 = buffer.data(slg + 174);
    const auto *slg_175 = buffer.data(slg + 175);
    const auto *slg_176 = buffer.data(slg + 176);

#pragma omp simd aligned(t_126, t_127, t_128, t_129, pc_x, pc_y, pc_z, skg_45, skg_90, skg_93, \
                         slf0_60, slf0_63, slf1_60, slf1_63, slg_90, \
                         slg_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_14 * skg_90[k]
                   + f_1 * slf0_60[k]
                   - f_2 * slf1_60[k]
                   + f_3 * pc_x[k] * slg_90[k];

        t_127[k] = f_11 * skg_45[k]
                   + f_3 * pc_y[k] * slg_90[k];

        t_128[k] = f_3 * pc_z[k] * slg_90[k];

        t_129[k] = f_14 * skg_93[k]
                   + f_4 * slf0_63[k]
                   - f_5 * slf1_63[k]
                   + f_3 * pc_x[k] * slg_93[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, pc_x, pc_y, skg_47, skg_95, skg_96, slf0_65, \
                         slf0_66, slf1_65, slf1_66, slg_92, slg_95, \
                         slg_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_11 * skg_47[k]
                   + f_3 * pc_y[k] * slg_92[k];

        t_131[k] = f_14 * skg_95[k]
                   + f_4 * slf0_65[k]
                   - f_5 * slf1_65[k]
                   + f_3 * pc_x[k] * slg_95[k];

        t_132[k] = f_14 * skg_96[k]
                   + f_6 * slf0_66[k]
                   - f_7 * slf1_66[k]
                   + f_3 * pc_x[k] * slg_96[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, pc_x, pc_y, pc_z, skg_50, skg_99, \
                         skg_100, slf0_69, slf1_69, slg_93, slg_95, slg_99, \
                         slg_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_3 * pc_z[k] * slg_93[k];

        t_134[k] = f_11 * skg_50[k]
                   + f_3 * pc_y[k] * slg_95[k];

        t_135[k] = f_14 * skg_99[k]
                   + f_6 * slf0_69[k]
                   - f_7 * slf1_69[k]
                   + f_3 * pc_x[k] * slg_99[k];

        t_136[k] = f_14 * skg_100[k]
                   + f_3 * pc_x[k] * slg_100[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, t_140, pc_x, skg_101, skg_102, skg_103, skg_104, \
                         slg_101, slg_102, slg_103, slg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_14 * skg_101[k]
                   + f_3 * pc_x[k] * slg_101[k];

        t_138[k] = f_14 * skg_102[k]
                   + f_3 * pc_x[k] * slg_102[k];

        t_139[k] = f_14 * skg_103[k]
                   + f_3 * pc_x[k] * slg_103[k];

        t_140[k] = f_14 * skg_104[k]
                   + f_3 * pc_x[k] * slg_104[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, pc_y, pc_z, skg_55, skg_57, slf0_66, slf0_68, \
                         slf1_66, slf1_68, slg_100, slg_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_11 * skg_55[k]
                   + f_1 * slf0_66[k]
                   - f_2 * slf1_66[k]
                   + f_3 * pc_y[k] * slg_100[k];

        t_142[k] = f_3 * pc_z[k] * slg_100[k];

        t_143[k] = f_11 * skg_57[k]
                   + f_4 * slf0_68[k]
                   - f_5 * slf1_68[k]
                   + f_3 * pc_y[k] * slg_102[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, pb_z, pc_y, pc_z, skh0_63, skg_58, \
                         skg_59, skh1_63, slf0_69, slf1_69, slg_103, \
                         slg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_11 * skg_58[k]
                   + f_6 * slf0_69[k]
                   - f_7 * slf1_69[k]
                   + f_3 * pc_y[k] * slg_103[k];

        t_145[k] = f_11 * skg_59[k]
                   + f_3 * pc_y[k] * slg_104[k];

        t_146[k] = f_1 * slf0_69[k]
                   - f_2 * slf1_69[k]
                   + f_3 * pc_z[k] * slg_104[k];

        t_147[k] = pb_z[k] * skh0_63[k]
                   - f_8 * pc_z[k] * skh1_63[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, t_151, pb_z, pc_y, pc_z, skh0_66, skg_45, \
                         skg_60, skg_62, skh1_66, slg_105, slg_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_10 * skg_60[k]
                   + f_3 * pc_y[k] * slg_105[k];

        t_149[k] = f_9 * skg_45[k]
                   + f_3 * pc_z[k] * slg_105[k];

        t_150[k] = pb_z[k] * skh0_66[k]
                   - f_8 * pc_z[k] * skh1_66[k];

        t_151[k] = f_10 * skg_62[k]
                   + f_3 * pc_y[k] * slg_107[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, pb_z, pc_x, pc_z, skh0_69, skg_48, skg_110, \
                         skh1_69, slf0_75, slf1_75, slg_108, slg_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_14 * skg_110[k]
                   + f_4 * slf0_75[k]
                   - f_5 * slf1_75[k]
                   + f_3 * pc_x[k] * slg_110[k];

        t_153[k] = pb_z[k] * skh0_69[k]
                   - f_8 * pc_z[k] * skh1_69[k];

        t_154[k] = f_9 * skg_48[k]
                   + f_3 * pc_z[k] * slg_108[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, pc_x, pc_y, skg_65, skg_114, skg_115, \
                         skg_116, slf0_79, slf1_79, slg_110, slg_114, slg_115, \
                         slg_116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_10 * skg_65[k]
                   + f_3 * pc_y[k] * slg_110[k];

        t_156[k] = f_14 * skg_114[k]
                   + f_6 * slf0_79[k]
                   - f_7 * slf1_79[k]
                   + f_3 * pc_x[k] * slg_114[k];

        t_157[k] = f_14 * skg_115[k]
                   + f_3 * pc_x[k] * slg_115[k];

        t_158[k] = f_14 * skg_116[k]
                   + f_3 * pc_x[k] * slg_116[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, pb_z, pc_x, pc_z, skh0_78, skg_117, \
                         skg_118, skg_119, skh1_78, slg_117, slg_118, \
                         slg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_14 * skg_117[k]
                   + f_3 * pc_x[k] * slg_117[k];

        t_160[k] = f_14 * skg_118[k]
                   + f_3 * pc_x[k] * slg_118[k];

        t_161[k] = f_14 * skg_119[k]
                   + f_3 * pc_x[k] * slg_119[k];

        t_162[k] = pb_z[k] * skh0_78[k]
                   - f_8 * pc_z[k] * skh1_78[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, pc_y, pc_z, skg_55, skg_72, skg_73, slf0_78, \
                         slf0_79, slf1_78, slf1_79, slg_115, slg_117, \
                         slg_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_9 * skg_55[k]
                   + f_3 * pc_z[k] * slg_115[k];

        t_164[k] = f_10 * skg_72[k]
                   + f_4 * slf0_78[k]
                   - f_5 * slf1_78[k]
                   + f_3 * pc_y[k] * slg_117[k];

        t_165[k] = f_10 * skg_73[k]
                   + f_6 * slf0_79[k]
                   - f_7 * slf1_79[k]
                   + f_3 * pc_y[k] * slg_118[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, pb_y, pc_y, pc_z, skh0_105, skg_59, \
                         skg_74, skg_75, skh1_105, slf0_79, slf1_79, slg_119, \
                         slg_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_10 * skg_74[k]
                   + f_3 * pc_y[k] * slg_119[k];

        t_167[k] = f_9 * skg_59[k]
                   + f_1 * slf0_79[k]
                   - f_2 * slf1_79[k]
                   + f_3 * pc_z[k] * slg_119[k];

        t_168[k] = pb_y[k] * skh0_105[k]
                   - f_8 * pc_y[k] * skh1_105[k];

        t_169[k] = f_9 * skg_75[k]
                   + f_3 * pc_y[k] * slg_120[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, pb_y, pc_y, pc_z, skh0_108, skh0_110, \
                         skg_60, skg_76, skg_77, skh1_108, skh1_110, slg_120, \
                         slg_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_10 * skg_60[k]
                   + f_3 * pc_z[k] * slg_120[k];

        t_171[k] = pb_y[k] * skh0_108[k]
                   + f_10 * skg_76[k]
                   - f_8 * pc_y[k] * skh1_108[k];

        t_172[k] = f_9 * skg_77[k]
                   + f_3 * pc_y[k] * slg_122[k];

        t_173[k] = pb_y[k] * skh0_110[k]
                   - f_8 * pc_y[k] * skh1_110[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, pb_y, pc_y, pc_z, skh0_111, skh0_114, \
                         skg_63, skg_78, skg_80, skh1_111, skh1_114, slg_123, \
                         slg_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = pb_y[k] * skh0_111[k]
                   + f_11 * skg_78[k]
                   - f_8 * pc_y[k] * skh1_111[k];

        t_175[k] = f_10 * skg_63[k]
                   + f_3 * pc_z[k] * slg_123[k];

        t_176[k] = f_9 * skg_80[k]
                   + f_3 * pc_y[k] * slg_125[k];

        t_177[k] = pb_y[k] * skh0_114[k]
                   - f_8 * pc_y[k] * skh1_114[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, t_181, t_182, pc_x, skg_130, skg_131, skg_132, \
                         skg_133, skg_134, slg_130, slg_131, slg_132, slg_133, \
                         slg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_14 * skg_130[k]
                   + f_3 * pc_x[k] * slg_130[k];

        t_179[k] = f_14 * skg_131[k]
                   + f_3 * pc_x[k] * slg_131[k];

        t_180[k] = f_14 * skg_132[k]
                   + f_3 * pc_x[k] * slg_132[k];

        t_181[k] = f_14 * skg_133[k]
                   + f_3 * pc_x[k] * slg_133[k];

        t_182[k] = f_14 * skg_134[k]
                   + f_3 * pc_x[k] * slg_134[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, pc_y, pc_z, skg_70, skg_85, skg_87, slf0_86, \
                         slf0_88, slf1_86, slf1_88, slg_130, slg_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_9 * skg_85[k]
                   + f_1 * slf0_86[k]
                   - f_2 * slf1_86[k]
                   + f_3 * pc_y[k] * slg_130[k];

        t_184[k] = f_10 * skg_70[k]
                   + f_3 * pc_z[k] * slg_130[k];

        t_185[k] = f_9 * skg_87[k]
                   + f_4 * slf0_88[k]
                   - f_5 * slf1_88[k]
                   + f_3 * pc_y[k] * slg_132[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, pb_y, pc_y, skh0_125, skg_88, skg_89, skh1_125, \
                         slf0_89, slf1_89, slg_133, slg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_9 * skg_88[k]
                   + f_6 * slf0_89[k]
                   - f_7 * slf1_89[k]
                   + f_3 * pc_y[k] * slg_133[k];

        t_187[k] = f_9 * skg_89[k]
                   + f_3 * pc_y[k] * slg_134[k];

        t_188[k] = pb_y[k] * skh0_125[k]
                   - f_8 * pc_y[k] * skh1_125[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, pc_x, pc_y, pc_z, skg_75, skg_135, \
                         skg_138, slf0_90, slf0_93, slf1_90, slf1_93, slg_135, \
                         slg_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_14 * skg_135[k]
                   + f_1 * slf0_90[k]
                   - f_2 * slf1_90[k]
                   + f_3 * pc_x[k] * slg_135[k];

        t_190[k] = f_3 * pc_y[k] * slg_135[k];

        t_191[k] = f_11 * skg_75[k]
                   + f_3 * pc_z[k] * slg_135[k];

        t_192[k] = f_14 * skg_138[k]
                   + f_4 * slf0_93[k]
                   - f_5 * slf1_93[k]
                   + f_3 * pc_x[k] * slg_138[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, pc_x, pc_y, skg_140, skg_141, slf0_95, slf0_96, \
                         slf1_95, slf1_96, slg_137, slg_140, slg_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_3 * pc_y[k] * slg_137[k];

        t_194[k] = f_14 * skg_140[k]
                   + f_4 * slf0_95[k]
                   - f_5 * slf1_95[k]
                   + f_3 * pc_x[k] * slg_140[k];

        t_195[k] = f_14 * skg_141[k]
                   + f_6 * slf0_96[k]
                   - f_7 * slf1_96[k]
                   + f_3 * pc_x[k] * slg_141[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, pc_x, pc_y, pc_z, skg_78, skg_144, \
                         skg_145, slf0_99, slf1_99, slg_138, slg_140, slg_144, \
                         slg_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = f_11 * skg_78[k]
                   + f_3 * pc_z[k] * slg_138[k];

        t_197[k] = f_3 * pc_y[k] * slg_140[k];

        t_198[k] = f_14 * skg_144[k]
                   + f_6 * slf0_99[k]
                   - f_7 * slf1_99[k]
                   + f_3 * pc_x[k] * slg_144[k];

        t_199[k] = f_14 * skg_145[k]
                   + f_3 * pc_x[k] * slg_145[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, pc_x, skg_146, skg_147, skg_148, skg_149, \
                         slg_146, slg_147, slg_148, slg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_14 * skg_146[k]
                   + f_3 * pc_x[k] * slg_146[k];

        t_201[k] = f_14 * skg_147[k]
                   + f_3 * pc_x[k] * slg_147[k];

        t_202[k] = f_14 * skg_148[k]
                   + f_3 * pc_x[k] * slg_148[k];

        t_203[k] = f_14 * skg_149[k]
                   + f_3 * pc_x[k] * slg_149[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, pc_y, pc_z, skg_85, slf0_96, slf0_98, \
                         slf0_99, slf1_96, slf1_98, slf1_99, slg_145, slg_147, \
                         slg_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_1 * slf0_96[k]
                   - f_2 * slf1_96[k]
                   + f_3 * pc_y[k] * slg_145[k];

        t_205[k] = f_11 * skg_85[k]
                   + f_3 * pc_z[k] * slg_145[k];

        t_206[k] = f_4 * slf0_98[k]
                   - f_5 * slf1_98[k]
                   + f_3 * pc_y[k] * slg_147[k];

        t_207[k] = f_6 * slf0_99[k]
                   - f_7 * slf1_99[k]
                   + f_3 * pc_y[k] * slg_148[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, pc_x, pc_y, pc_z, skg_89, skg_90, \
                         skg_150, slf0_99, slf0_100, slf1_99, slf1_100, slg_149, \
                         slg_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_3 * pc_y[k] * slg_149[k];

        t_209[k] = f_11 * skg_89[k]
                   + f_1 * slf0_99[k]
                   - f_2 * slf1_99[k]
                   + f_3 * pc_z[k] * slg_149[k];

        t_210[k] = f_15 * skg_150[k]
                   + f_1 * slf0_100[k]
                   - f_2 * slf1_100[k]
                   + f_3 * pc_x[k] * slg_150[k];

        t_211[k] = f_15 * skg_90[k]
                   + f_3 * pc_y[k] * slg_150[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, pc_x, pc_y, pc_z, skg_92, skg_153, slf0_103, \
                         slf1_103, slg_150, slg_152, slg_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_3 * pc_z[k] * slg_150[k];

        t_213[k] = f_15 * skg_153[k]
                   + f_4 * slf0_103[k]
                   - f_5 * slf1_103[k]
                   + f_3 * pc_x[k] * slg_153[k];

        t_214[k] = f_15 * skg_92[k]
                   + f_3 * pc_y[k] * slg_152[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, pc_x, pc_z, skg_155, skg_156, slf0_105, \
                         slf0_106, slf1_105, slf1_106, slg_153, slg_155, \
                         slg_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = f_15 * skg_155[k]
                   + f_4 * slf0_105[k]
                   - f_5 * slf1_105[k]
                   + f_3 * pc_x[k] * slg_155[k];

        t_216[k] = f_15 * skg_156[k]
                   + f_6 * slf0_106[k]
                   - f_7 * slf1_106[k]
                   + f_3 * pc_x[k] * slg_156[k];

        t_217[k] = f_3 * pc_z[k] * slg_153[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, t_221, pc_x, pc_y, skg_95, skg_159, skg_160, \
                         skg_161, slf0_109, slf1_109, slg_155, slg_159, slg_160, \
                         slg_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_15 * skg_95[k]
                   + f_3 * pc_y[k] * slg_155[k];

        t_219[k] = f_15 * skg_159[k]
                   + f_6 * slf0_109[k]
                   - f_7 * slf1_109[k]
                   + f_3 * pc_x[k] * slg_159[k];

        t_220[k] = f_15 * skg_160[k]
                   + f_3 * pc_x[k] * slg_160[k];

        t_221[k] = f_15 * skg_161[k]
                   + f_3 * pc_x[k] * slg_161[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, pc_x, pc_y, skg_100, skg_162, skg_163, \
                         skg_164, slf0_106, slf1_106, slg_160, slg_162, slg_163, \
                         slg_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_15 * skg_162[k]
                   + f_3 * pc_x[k] * slg_162[k];

        t_223[k] = f_15 * skg_163[k]
                   + f_3 * pc_x[k] * slg_163[k];

        t_224[k] = f_15 * skg_164[k]
                   + f_3 * pc_x[k] * slg_164[k];

        t_225[k] = f_15 * skg_100[k]
                   + f_1 * slf0_106[k]
                   - f_2 * slf1_106[k]
                   + f_3 * pc_y[k] * slg_160[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, pc_y, pc_z, skg_102, skg_103, slf0_108, \
                         slf0_109, slf1_108, slf1_109, slg_160, slg_162, \
                         slg_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = f_3 * pc_z[k] * slg_160[k];

        t_227[k] = f_15 * skg_102[k]
                   + f_4 * slf0_108[k]
                   - f_5 * slf1_108[k]
                   + f_3 * pc_y[k] * slg_162[k];

        t_228[k] = f_15 * skg_103[k]
                   + f_6 * slf0_109[k]
                   - f_7 * slf1_109[k]
                   + f_3 * pc_y[k] * slg_163[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, pb_z, pc_y, pc_z, skh0_126, skg_104, \
                         skg_105, skh1_126, slf0_109, slf1_109, slg_164, \
                         slg_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = f_15 * skg_104[k]
                   + f_3 * pc_y[k] * slg_164[k];

        t_230[k] = f_1 * slf0_109[k]
                   - f_2 * slf1_109[k]
                   + f_3 * pc_z[k] * slg_164[k];

        t_231[k] = pb_z[k] * skh0_126[k]
                   - f_8 * pc_z[k] * skh1_126[k];

        t_232[k] = f_11 * skg_105[k]
                   + f_3 * pc_y[k] * slg_165[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, pb_z, pc_y, pc_z, skh0_129, skg_90, skg_107, \
                         skh1_129, slg_165, slg_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = f_9 * skg_90[k]
                   + f_3 * pc_z[k] * slg_165[k];

        t_234[k] = pb_z[k] * skh0_129[k]
                   - f_8 * pc_z[k] * skh1_129[k];

        t_235[k] = f_11 * skg_107[k]
                   + f_3 * pc_y[k] * slg_167[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, pb_z, pc_x, pc_z, skh0_132, skg_93, skg_170, \
                         skh1_132, slf0_115, slf1_115, slg_168, \
                         slg_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = f_15 * skg_170[k]
                   + f_4 * slf0_115[k]
                   - f_5 * slf1_115[k]
                   + f_3 * pc_x[k] * slg_170[k];

        t_237[k] = pb_z[k] * skh0_132[k]
                   - f_8 * pc_z[k] * skh1_132[k];

        t_238[k] = f_9 * skg_93[k]
                   + f_3 * pc_z[k] * slg_168[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, pc_x, pc_y, skg_110, skg_174, skg_175, \
                         skg_176, slf0_119, slf1_119, slg_170, slg_174, slg_175, \
                         slg_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_11 * skg_110[k]
                   + f_3 * pc_y[k] * slg_170[k];

        t_240[k] = f_15 * skg_174[k]
                   + f_6 * slf0_119[k]
                   - f_7 * slf1_119[k]
                   + f_3 * pc_x[k] * slg_174[k];

        t_241[k] = f_15 * skg_175[k]
                   + f_3 * pc_x[k] * slg_175[k];

        t_242[k] = f_15 * skg_176[k]
                   + f_3 * pc_x[k] * slg_176[k];
    }
}

static auto
compute_prim_slh_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t skh0,
                                                          const size_t skg, const size_t skh1,
                                                          const size_t slf0, const size_t slf1,
                                                          const size_t slg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_14 = 2.5 / q;
    const auto f_15 = 2.0 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *skh0_141 = buffer.data(skh0 + 141);
    const auto *skh0_189 = buffer.data(skh0 + 189);
    const auto *skh0_192 = buffer.data(skh0 + 192);
    const auto *skh0_194 = buffer.data(skh0 + 194);
    const auto *skh0_195 = buffer.data(skh0 + 195);
    const auto *skh0_198 = buffer.data(skh0 + 198);
    const auto *skh0_209 = buffer.data(skh0 + 209);
    const auto *skh0_210 = buffer.data(skh0 + 210);
    const auto *skh0_213 = buffer.data(skh0 + 213);
    const auto *skh0_216 = buffer.data(skh0 + 216);
    const auto *skh0_225 = buffer.data(skh0 + 225);

    const auto *skg_100 = buffer.data(skg + 100);
    const auto *skg_104 = buffer.data(skg + 104);
    const auto *skg_105 = buffer.data(skg + 105);
    const auto *skg_108 = buffer.data(skg + 108);
    const auto *skg_115 = buffer.data(skg + 115);
    const auto *skg_117 = buffer.data(skg + 117);
    const auto *skg_118 = buffer.data(skg + 118);
    const auto *skg_119 = buffer.data(skg + 119);
    const auto *skg_120 = buffer.data(skg + 120);
    const auto *skg_122 = buffer.data(skg + 122);
    const auto *skg_123 = buffer.data(skg + 123);
    const auto *skg_125 = buffer.data(skg + 125);
    const auto *skg_130 = buffer.data(skg + 130);
    const auto *skg_132 = buffer.data(skg + 132);
    const auto *skg_133 = buffer.data(skg + 133);
    const auto *skg_134 = buffer.data(skg + 134);
    const auto *skg_135 = buffer.data(skg + 135);
    const auto *skg_136 = buffer.data(skg + 136);
    const auto *skg_137 = buffer.data(skg + 137);
    const auto *skg_138 = buffer.data(skg + 138);
    const auto *skg_140 = buffer.data(skg + 140);
    const auto *skg_145 = buffer.data(skg + 145);
    const auto *skg_147 = buffer.data(skg + 147);
    const auto *skg_148 = buffer.data(skg + 148);
    const auto *skg_149 = buffer.data(skg + 149);
    const auto *skg_150 = buffer.data(skg + 150);
    const auto *skg_152 = buffer.data(skg + 152);
    const auto *skg_153 = buffer.data(skg + 153);
    const auto *skg_155 = buffer.data(skg + 155);
    const auto *skg_160 = buffer.data(skg + 160);
    const auto *skg_162 = buffer.data(skg + 162);
    const auto *skg_163 = buffer.data(skg + 163);
    const auto *skg_164 = buffer.data(skg + 164);
    const auto *skg_165 = buffer.data(skg + 165);
    const auto *skg_167 = buffer.data(skg + 167);
    const auto *skg_170 = buffer.data(skg + 170);
    const auto *skg_177 = buffer.data(skg + 177);
    const auto *skg_178 = buffer.data(skg + 178);
    const auto *skg_179 = buffer.data(skg + 179);
    const auto *skg_180 = buffer.data(skg + 180);
    const auto *skg_183 = buffer.data(skg + 183);
    const auto *skg_185 = buffer.data(skg + 185);
    const auto *skg_186 = buffer.data(skg + 186);
    const auto *skg_189 = buffer.data(skg + 189);
    const auto *skg_190 = buffer.data(skg + 190);
    const auto *skg_191 = buffer.data(skg + 191);
    const auto *skg_192 = buffer.data(skg + 192);
    const auto *skg_193 = buffer.data(skg + 193);
    const auto *skg_194 = buffer.data(skg + 194);
    const auto *skg_205 = buffer.data(skg + 205);
    const auto *skg_206 = buffer.data(skg + 206);
    const auto *skg_207 = buffer.data(skg + 207);
    const auto *skg_208 = buffer.data(skg + 208);
    const auto *skg_209 = buffer.data(skg + 209);
    const auto *skg_210 = buffer.data(skg + 210);
    const auto *skg_213 = buffer.data(skg + 213);
    const auto *skg_215 = buffer.data(skg + 215);
    const auto *skg_216 = buffer.data(skg + 216);
    const auto *skg_219 = buffer.data(skg + 219);
    const auto *skg_220 = buffer.data(skg + 220);
    const auto *skg_221 = buffer.data(skg + 221);
    const auto *skg_222 = buffer.data(skg + 222);
    const auto *skg_223 = buffer.data(skg + 223);
    const auto *skg_224 = buffer.data(skg + 224);
    const auto *skg_225 = buffer.data(skg + 225);
    const auto *skg_228 = buffer.data(skg + 228);
    const auto *skg_230 = buffer.data(skg + 230);
    const auto *skg_231 = buffer.data(skg + 231);
    const auto *skg_234 = buffer.data(skg + 234);
    const auto *skg_235 = buffer.data(skg + 235);
    const auto *skg_236 = buffer.data(skg + 236);
    const auto *skg_237 = buffer.data(skg + 237);
    const auto *skg_238 = buffer.data(skg + 238);
    const auto *skg_239 = buffer.data(skg + 239);
    const auto *skg_245 = buffer.data(skg + 245);
    const auto *skg_249 = buffer.data(skg + 249);
    const auto *skg_250 = buffer.data(skg + 250);
    const auto *skg_251 = buffer.data(skg + 251);
    const auto *skg_252 = buffer.data(skg + 252);
    const auto *skg_253 = buffer.data(skg + 253);
    const auto *skg_254 = buffer.data(skg + 254);
    const auto *skg_255 = buffer.data(skg + 255);

    const auto *skh1_141 = buffer.data(skh1 + 141);
    const auto *skh1_189 = buffer.data(skh1 + 189);
    const auto *skh1_192 = buffer.data(skh1 + 192);
    const auto *skh1_194 = buffer.data(skh1 + 194);
    const auto *skh1_195 = buffer.data(skh1 + 195);
    const auto *skh1_198 = buffer.data(skh1 + 198);
    const auto *skh1_209 = buffer.data(skh1 + 209);
    const auto *skh1_210 = buffer.data(skh1 + 210);
    const auto *skh1_213 = buffer.data(skh1 + 213);
    const auto *skh1_216 = buffer.data(skh1 + 216);
    const auto *skh1_225 = buffer.data(skh1 + 225);

    const auto *slf0_118 = buffer.data(slf0 + 118);
    const auto *slf0_119 = buffer.data(slf0 + 119);
    const auto *slf0_120 = buffer.data(slf0 + 120);
    const auto *slf0_123 = buffer.data(slf0 + 123);
    const auto *slf0_125 = buffer.data(slf0 + 125);
    const auto *slf0_126 = buffer.data(slf0 + 126);
    const auto *slf0_128 = buffer.data(slf0 + 128);
    const auto *slf0_129 = buffer.data(slf0 + 129);
    const auto *slf0_136 = buffer.data(slf0 + 136);
    const auto *slf0_138 = buffer.data(slf0 + 138);
    const auto *slf0_139 = buffer.data(slf0 + 139);
    const auto *slf0_140 = buffer.data(slf0 + 140);
    const auto *slf0_143 = buffer.data(slf0 + 143);
    const auto *slf0_145 = buffer.data(slf0 + 145);
    const auto *slf0_146 = buffer.data(slf0 + 146);
    const auto *slf0_148 = buffer.data(slf0 + 148);
    const auto *slf0_149 = buffer.data(slf0 + 149);
    const auto *slf0_150 = buffer.data(slf0 + 150);
    const auto *slf0_153 = buffer.data(slf0 + 153);
    const auto *slf0_155 = buffer.data(slf0 + 155);
    const auto *slf0_156 = buffer.data(slf0 + 156);
    const auto *slf0_158 = buffer.data(slf0 + 158);
    const auto *slf0_159 = buffer.data(slf0 + 159);
    const auto *slf0_165 = buffer.data(slf0 + 165);
    const auto *slf0_168 = buffer.data(slf0 + 168);
    const auto *slf0_169 = buffer.data(slf0 + 169);
    const auto *slf0_170 = buffer.data(slf0 + 170);

    const auto *slf1_118 = buffer.data(slf1 + 118);
    const auto *slf1_119 = buffer.data(slf1 + 119);
    const auto *slf1_120 = buffer.data(slf1 + 120);
    const auto *slf1_123 = buffer.data(slf1 + 123);
    const auto *slf1_125 = buffer.data(slf1 + 125);
    const auto *slf1_126 = buffer.data(slf1 + 126);
    const auto *slf1_128 = buffer.data(slf1 + 128);
    const auto *slf1_129 = buffer.data(slf1 + 129);
    const auto *slf1_136 = buffer.data(slf1 + 136);
    const auto *slf1_138 = buffer.data(slf1 + 138);
    const auto *slf1_139 = buffer.data(slf1 + 139);
    const auto *slf1_140 = buffer.data(slf1 + 140);
    const auto *slf1_143 = buffer.data(slf1 + 143);
    const auto *slf1_145 = buffer.data(slf1 + 145);
    const auto *slf1_146 = buffer.data(slf1 + 146);
    const auto *slf1_148 = buffer.data(slf1 + 148);
    const auto *slf1_149 = buffer.data(slf1 + 149);
    const auto *slf1_150 = buffer.data(slf1 + 150);
    const auto *slf1_153 = buffer.data(slf1 + 153);
    const auto *slf1_155 = buffer.data(slf1 + 155);
    const auto *slf1_156 = buffer.data(slf1 + 156);
    const auto *slf1_158 = buffer.data(slf1 + 158);
    const auto *slf1_159 = buffer.data(slf1 + 159);
    const auto *slf1_165 = buffer.data(slf1 + 165);
    const auto *slf1_168 = buffer.data(slf1 + 168);
    const auto *slf1_169 = buffer.data(slf1 + 169);
    const auto *slf1_170 = buffer.data(slf1 + 170);

    const auto *slg_175 = buffer.data(slg + 175);
    const auto *slg_177 = buffer.data(slg + 177);
    const auto *slg_178 = buffer.data(slg + 178);
    const auto *slg_179 = buffer.data(slg + 179);
    const auto *slg_180 = buffer.data(slg + 180);
    const auto *slg_182 = buffer.data(slg + 182);
    const auto *slg_183 = buffer.data(slg + 183);
    const auto *slg_185 = buffer.data(slg + 185);
    const auto *slg_186 = buffer.data(slg + 186);
    const auto *slg_189 = buffer.data(slg + 189);
    const auto *slg_190 = buffer.data(slg + 190);
    const auto *slg_191 = buffer.data(slg + 191);
    const auto *slg_192 = buffer.data(slg + 192);
    const auto *slg_193 = buffer.data(slg + 193);
    const auto *slg_194 = buffer.data(slg + 194);
    const auto *slg_195 = buffer.data(slg + 195);
    const auto *slg_197 = buffer.data(slg + 197);
    const auto *slg_198 = buffer.data(slg + 198);
    const auto *slg_200 = buffer.data(slg + 200);
    const auto *slg_205 = buffer.data(slg + 205);
    const auto *slg_206 = buffer.data(slg + 206);
    const auto *slg_207 = buffer.data(slg + 207);
    const auto *slg_208 = buffer.data(slg + 208);
    const auto *slg_209 = buffer.data(slg + 209);
    const auto *slg_210 = buffer.data(slg + 210);
    const auto *slg_212 = buffer.data(slg + 212);
    const auto *slg_213 = buffer.data(slg + 213);
    const auto *slg_215 = buffer.data(slg + 215);
    const auto *slg_216 = buffer.data(slg + 216);
    const auto *slg_219 = buffer.data(slg + 219);
    const auto *slg_220 = buffer.data(slg + 220);
    const auto *slg_221 = buffer.data(slg + 221);
    const auto *slg_222 = buffer.data(slg + 222);
    const auto *slg_223 = buffer.data(slg + 223);
    const auto *slg_224 = buffer.data(slg + 224);
    const auto *slg_225 = buffer.data(slg + 225);
    const auto *slg_227 = buffer.data(slg + 227);
    const auto *slg_228 = buffer.data(slg + 228);
    const auto *slg_230 = buffer.data(slg + 230);
    const auto *slg_231 = buffer.data(slg + 231);
    const auto *slg_234 = buffer.data(slg + 234);
    const auto *slg_235 = buffer.data(slg + 235);
    const auto *slg_236 = buffer.data(slg + 236);
    const auto *slg_237 = buffer.data(slg + 237);
    const auto *slg_238 = buffer.data(slg + 238);
    const auto *slg_239 = buffer.data(slg + 239);
    const auto *slg_240 = buffer.data(slg + 240);
    const auto *slg_242 = buffer.data(slg + 242);
    const auto *slg_243 = buffer.data(slg + 243);
    const auto *slg_245 = buffer.data(slg + 245);
    const auto *slg_249 = buffer.data(slg + 249);
    const auto *slg_250 = buffer.data(slg + 250);
    const auto *slg_251 = buffer.data(slg + 251);
    const auto *slg_252 = buffer.data(slg + 252);
    const auto *slg_253 = buffer.data(slg + 253);
    const auto *slg_254 = buffer.data(slg + 254);
    const auto *slg_255 = buffer.data(slg + 255);

#pragma omp simd aligned(t_243, t_244, t_245, t_246, pb_z, pc_x, pc_z, skh0_141, skg_177, \
                         skg_178, skg_179, skh1_141, slg_177, slg_178, \
                         slg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = f_15 * skg_177[k]
                   + f_3 * pc_x[k] * slg_177[k];

        t_244[k] = f_15 * skg_178[k]
                   + f_3 * pc_x[k] * slg_178[k];

        t_245[k] = f_15 * skg_179[k]
                   + f_3 * pc_x[k] * slg_179[k];

        t_246[k] = pb_z[k] * skh0_141[k]
                   - f_8 * pc_z[k] * skh1_141[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, pc_y, pc_z, skg_100, skg_117, skg_118, slf0_118, \
                         slf0_119, slf1_118, slf1_119, slg_175, slg_177, \
                         slg_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = f_9 * skg_100[k]
                   + f_3 * pc_z[k] * slg_175[k];

        t_248[k] = f_11 * skg_117[k]
                   + f_4 * slf0_118[k]
                   - f_5 * slf1_118[k]
                   + f_3 * pc_y[k] * slg_177[k];

        t_249[k] = f_11 * skg_118[k]
                   + f_6 * slf0_119[k]
                   - f_7 * slf1_119[k]
                   + f_3 * pc_y[k] * slg_178[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, pc_x, pc_y, pc_z, skg_104, skg_119, skg_180, \
                         slf0_119, slf0_120, slf1_119, slf1_120, slg_179, \
                         slg_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = f_11 * skg_119[k]
                   + f_3 * pc_y[k] * slg_179[k];

        t_251[k] = f_9 * skg_104[k]
                   + f_1 * slf0_119[k]
                   - f_2 * slf1_119[k]
                   + f_3 * pc_z[k] * slg_179[k];

        t_252[k] = f_15 * skg_180[k]
                   + f_1 * slf0_120[k]
                   - f_2 * slf1_120[k]
                   + f_3 * pc_x[k] * slg_180[k];
    }

#pragma omp simd aligned(t_253, t_254, t_255, t_256, pc_x, pc_y, pc_z, skg_105, skg_120, \
                         skg_122, skg_183, slf0_123, slf1_123, slg_180, slg_182, \
                         slg_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_253[k] = f_10 * skg_120[k]
                   + f_3 * pc_y[k] * slg_180[k];

        t_254[k] = f_10 * skg_105[k]
                   + f_3 * pc_z[k] * slg_180[k];

        t_255[k] = f_15 * skg_183[k]
                   + f_4 * slf0_123[k]
                   - f_5 * slf1_123[k]
                   + f_3 * pc_x[k] * slg_183[k];

        t_256[k] = f_10 * skg_122[k]
                   + f_3 * pc_y[k] * slg_182[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, pc_x, pc_z, skg_108, skg_185, skg_186, slf0_125, \
                         slf0_126, slf1_125, slf1_126, slg_183, slg_185, \
                         slg_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_15 * skg_185[k]
                   + f_4 * slf0_125[k]
                   - f_5 * slf1_125[k]
                   + f_3 * pc_x[k] * slg_185[k];

        t_258[k] = f_15 * skg_186[k]
                   + f_6 * slf0_126[k]
                   - f_7 * slf1_126[k]
                   + f_3 * pc_x[k] * slg_186[k];

        t_259[k] = f_10 * skg_108[k]
                   + f_3 * pc_z[k] * slg_183[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pc_x, pc_y, skg_125, skg_189, skg_190, \
                         skg_191, slf0_129, slf1_129, slg_185, slg_189, slg_190, \
                         slg_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_10 * skg_125[k]
                   + f_3 * pc_y[k] * slg_185[k];

        t_261[k] = f_15 * skg_189[k]
                   + f_6 * slf0_129[k]
                   - f_7 * slf1_129[k]
                   + f_3 * pc_x[k] * slg_189[k];

        t_262[k] = f_15 * skg_190[k]
                   + f_3 * pc_x[k] * slg_190[k];

        t_263[k] = f_15 * skg_191[k]
                   + f_3 * pc_x[k] * slg_191[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, pc_x, pc_y, skg_130, skg_192, skg_193, \
                         skg_194, slf0_126, slf1_126, slg_190, slg_192, slg_193, \
                         slg_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_15 * skg_192[k]
                   + f_3 * pc_x[k] * slg_192[k];

        t_265[k] = f_15 * skg_193[k]
                   + f_3 * pc_x[k] * slg_193[k];

        t_266[k] = f_15 * skg_194[k]
                   + f_3 * pc_x[k] * slg_194[k];

        t_267[k] = f_10 * skg_130[k]
                   + f_1 * slf0_126[k]
                   - f_2 * slf1_126[k]
                   + f_3 * pc_y[k] * slg_190[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, pc_y, pc_z, skg_115, skg_132, skg_133, slf0_128, \
                         slf0_129, slf1_128, slf1_129, slg_190, slg_192, \
                         slg_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_10 * skg_115[k]
                   + f_3 * pc_z[k] * slg_190[k];

        t_269[k] = f_10 * skg_132[k]
                   + f_4 * slf0_128[k]
                   - f_5 * slf1_128[k]
                   + f_3 * pc_y[k] * slg_192[k];

        t_270[k] = f_10 * skg_133[k]
                   + f_6 * slf0_129[k]
                   - f_7 * slf1_129[k]
                   + f_3 * pc_y[k] * slg_193[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, t_274, pb_y, pc_y, pc_z, skh0_189, skg_119, \
                         skg_134, skg_135, skh1_189, slf0_129, slf1_129, slg_194, \
                         slg_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = f_10 * skg_134[k]
                   + f_3 * pc_y[k] * slg_194[k];

        t_272[k] = f_10 * skg_119[k]
                   + f_1 * slf0_129[k]
                   - f_2 * slf1_129[k]
                   + f_3 * pc_z[k] * slg_194[k];

        t_273[k] = pb_y[k] * skh0_189[k]
                   - f_8 * pc_y[k] * skh1_189[k];

        t_274[k] = f_9 * skg_135[k]
                   + f_3 * pc_y[k] * slg_195[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, pb_y, pc_y, pc_z, skh0_192, skh0_194, \
                         skg_120, skg_136, skg_137, skh1_192, skh1_194, slg_195, \
                         slg_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = f_11 * skg_120[k]
                   + f_3 * pc_z[k] * slg_195[k];

        t_276[k] = pb_y[k] * skh0_192[k]
                   + f_10 * skg_136[k]
                   - f_8 * pc_y[k] * skh1_192[k];

        t_277[k] = f_9 * skg_137[k]
                   + f_3 * pc_y[k] * slg_197[k];

        t_278[k] = pb_y[k] * skh0_194[k]
                   - f_8 * pc_y[k] * skh1_194[k];
    }

#pragma omp simd aligned(t_279, t_280, t_281, t_282, pb_y, pc_y, pc_z, skh0_195, skh0_198, \
                         skg_123, skg_138, skg_140, skh1_195, skh1_198, slg_198, \
                         slg_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_279[k] = pb_y[k] * skh0_195[k]
                   + f_11 * skg_138[k]
                   - f_8 * pc_y[k] * skh1_195[k];

        t_280[k] = f_11 * skg_123[k]
                   + f_3 * pc_z[k] * slg_198[k];

        t_281[k] = f_9 * skg_140[k]
                   + f_3 * pc_y[k] * slg_200[k];

        t_282[k] = pb_y[k] * skh0_198[k]
                   - f_8 * pc_y[k] * skh1_198[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, t_286, t_287, pc_x, skg_205, skg_206, skg_207, \
                         skg_208, skg_209, slg_205, slg_206, slg_207, slg_208, \
                         slg_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_15 * skg_205[k]
                   + f_3 * pc_x[k] * slg_205[k];

        t_284[k] = f_15 * skg_206[k]
                   + f_3 * pc_x[k] * slg_206[k];

        t_285[k] = f_15 * skg_207[k]
                   + f_3 * pc_x[k] * slg_207[k];

        t_286[k] = f_15 * skg_208[k]
                   + f_3 * pc_x[k] * slg_208[k];

        t_287[k] = f_15 * skg_209[k]
                   + f_3 * pc_x[k] * slg_209[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, pc_y, pc_z, skg_130, skg_145, skg_147, slf0_136, \
                         slf0_138, slf1_136, slf1_138, slg_205, \
                         slg_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = f_9 * skg_145[k]
                   + f_1 * slf0_136[k]
                   - f_2 * slf1_136[k]
                   + f_3 * pc_y[k] * slg_205[k];

        t_289[k] = f_11 * skg_130[k]
                   + f_3 * pc_z[k] * slg_205[k];

        t_290[k] = f_9 * skg_147[k]
                   + f_4 * slf0_138[k]
                   - f_5 * slf1_138[k]
                   + f_3 * pc_y[k] * slg_207[k];
    }

#pragma omp simd aligned(t_291, t_292, t_293, pb_y, pc_y, skh0_209, skg_148, skg_149, \
                         skh1_209, slf0_139, slf1_139, slg_208, \
                         slg_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_291[k] = f_9 * skg_148[k]
                   + f_6 * slf0_139[k]
                   - f_7 * slf1_139[k]
                   + f_3 * pc_y[k] * slg_208[k];

        t_292[k] = f_9 * skg_149[k]
                   + f_3 * pc_y[k] * slg_209[k];

        t_293[k] = pb_y[k] * skh0_209[k]
                   - f_8 * pc_y[k] * skh1_209[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, pc_x, pc_y, pc_z, skg_135, skg_210, \
                         skg_213, slf0_140, slf0_143, slf1_140, slf1_143, slg_210, \
                         slg_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_15 * skg_210[k]
                   + f_1 * slf0_140[k]
                   - f_2 * slf1_140[k]
                   + f_3 * pc_x[k] * slg_210[k];

        t_295[k] = f_3 * pc_y[k] * slg_210[k];

        t_296[k] = f_15 * skg_135[k]
                   + f_3 * pc_z[k] * slg_210[k];

        t_297[k] = f_15 * skg_213[k]
                   + f_4 * slf0_143[k]
                   - f_5 * slf1_143[k]
                   + f_3 * pc_x[k] * slg_213[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, pc_x, pc_y, skg_215, skg_216, slf0_145, \
                         slf0_146, slf1_145, slf1_146, slg_212, slg_215, \
                         slg_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = f_3 * pc_y[k] * slg_212[k];

        t_299[k] = f_15 * skg_215[k]
                   + f_4 * slf0_145[k]
                   - f_5 * slf1_145[k]
                   + f_3 * pc_x[k] * slg_215[k];

        t_300[k] = f_15 * skg_216[k]
                   + f_6 * slf0_146[k]
                   - f_7 * slf1_146[k]
                   + f_3 * pc_x[k] * slg_216[k];
    }

#pragma omp simd aligned(t_301, t_302, t_303, t_304, pc_x, pc_y, pc_z, skg_138, skg_219, \
                         skg_220, slf0_149, slf1_149, slg_213, slg_215, slg_219, \
                         slg_220 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_301[k] = f_15 * skg_138[k]
                   + f_3 * pc_z[k] * slg_213[k];

        t_302[k] = f_3 * pc_y[k] * slg_215[k];

        t_303[k] = f_15 * skg_219[k]
                   + f_6 * slf0_149[k]
                   - f_7 * slf1_149[k]
                   + f_3 * pc_x[k] * slg_219[k];

        t_304[k] = f_15 * skg_220[k]
                   + f_3 * pc_x[k] * slg_220[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, pc_x, skg_221, skg_222, skg_223, skg_224, \
                         slg_221, slg_222, slg_223, slg_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = f_15 * skg_221[k]
                   + f_3 * pc_x[k] * slg_221[k];

        t_306[k] = f_15 * skg_222[k]
                   + f_3 * pc_x[k] * slg_222[k];

        t_307[k] = f_15 * skg_223[k]
                   + f_3 * pc_x[k] * slg_223[k];

        t_308[k] = f_15 * skg_224[k]
                   + f_3 * pc_x[k] * slg_224[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, pc_y, pc_z, skg_145, slf0_146, slf0_148, \
                         slf0_149, slf1_146, slf1_148, slf1_149, slg_220, slg_222, \
                         slg_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_1 * slf0_146[k]
                   - f_2 * slf1_146[k]
                   + f_3 * pc_y[k] * slg_220[k];

        t_310[k] = f_15 * skg_145[k]
                   + f_3 * pc_z[k] * slg_220[k];

        t_311[k] = f_4 * slf0_148[k]
                   - f_5 * slf1_148[k]
                   + f_3 * pc_y[k] * slg_222[k];

        t_312[k] = f_6 * slf0_149[k]
                   - f_7 * slf1_149[k]
                   + f_3 * pc_y[k] * slg_223[k];
    }

#pragma omp simd aligned(t_313, t_314, t_315, t_316, pc_x, pc_y, pc_z, skg_149, skg_150, \
                         skg_225, slf0_149, slf0_150, slf1_149, slf1_150, slg_224, \
                         slg_225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_313[k] = f_3 * pc_y[k] * slg_224[k];

        t_314[k] = f_15 * skg_149[k]
                   + f_1 * slf0_149[k]
                   - f_2 * slf1_149[k]
                   + f_3 * pc_z[k] * slg_224[k];

        t_315[k] = f_11 * skg_225[k]
                   + f_1 * slf0_150[k]
                   - f_2 * slf1_150[k]
                   + f_3 * pc_x[k] * slg_225[k];

        t_316[k] = f_14 * skg_150[k]
                   + f_3 * pc_y[k] * slg_225[k];
    }

#pragma omp simd aligned(t_317, t_318, t_319, pc_x, pc_y, pc_z, skg_152, skg_228, slf0_153, \
                         slf1_153, slg_225, slg_227, slg_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_317[k] = f_3 * pc_z[k] * slg_225[k];

        t_318[k] = f_11 * skg_228[k]
                   + f_4 * slf0_153[k]
                   - f_5 * slf1_153[k]
                   + f_3 * pc_x[k] * slg_228[k];

        t_319[k] = f_14 * skg_152[k]
                   + f_3 * pc_y[k] * slg_227[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, pc_x, pc_z, skg_230, skg_231, slf0_155, \
                         slf0_156, slf1_155, slf1_156, slg_228, slg_230, \
                         slg_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = f_11 * skg_230[k]
                   + f_4 * slf0_155[k]
                   - f_5 * slf1_155[k]
                   + f_3 * pc_x[k] * slg_230[k];

        t_321[k] = f_11 * skg_231[k]
                   + f_6 * slf0_156[k]
                   - f_7 * slf1_156[k]
                   + f_3 * pc_x[k] * slg_231[k];

        t_322[k] = f_3 * pc_z[k] * slg_228[k];
    }

#pragma omp simd aligned(t_323, t_324, t_325, t_326, pc_x, pc_y, skg_155, skg_234, skg_235, \
                         skg_236, slf0_159, slf1_159, slg_230, slg_234, slg_235, \
                         slg_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_323[k] = f_14 * skg_155[k]
                   + f_3 * pc_y[k] * slg_230[k];

        t_324[k] = f_11 * skg_234[k]
                   + f_6 * slf0_159[k]
                   - f_7 * slf1_159[k]
                   + f_3 * pc_x[k] * slg_234[k];

        t_325[k] = f_11 * skg_235[k]
                   + f_3 * pc_x[k] * slg_235[k];

        t_326[k] = f_11 * skg_236[k]
                   + f_3 * pc_x[k] * slg_236[k];
    }

#pragma omp simd aligned(t_327, t_328, t_329, t_330, pc_x, pc_y, skg_160, skg_237, skg_238, \
                         skg_239, slf0_156, slf1_156, slg_235, slg_237, slg_238, \
                         slg_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_327[k] = f_11 * skg_237[k]
                   + f_3 * pc_x[k] * slg_237[k];

        t_328[k] = f_11 * skg_238[k]
                   + f_3 * pc_x[k] * slg_238[k];

        t_329[k] = f_11 * skg_239[k]
                   + f_3 * pc_x[k] * slg_239[k];

        t_330[k] = f_14 * skg_160[k]
                   + f_1 * slf0_156[k]
                   - f_2 * slf1_156[k]
                   + f_3 * pc_y[k] * slg_235[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, pc_y, pc_z, skg_162, skg_163, slf0_158, \
                         slf0_159, slf1_158, slf1_159, slg_235, slg_237, \
                         slg_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = f_3 * pc_z[k] * slg_235[k];

        t_332[k] = f_14 * skg_162[k]
                   + f_4 * slf0_158[k]
                   - f_5 * slf1_158[k]
                   + f_3 * pc_y[k] * slg_237[k];

        t_333[k] = f_14 * skg_163[k]
                   + f_6 * slf0_159[k]
                   - f_7 * slf1_159[k]
                   + f_3 * pc_y[k] * slg_238[k];
    }

#pragma omp simd aligned(t_334, t_335, t_336, t_337, pb_z, pc_y, pc_z, skh0_210, skg_164, \
                         skg_165, skh1_210, slf0_159, slf1_159, slg_239, \
                         slg_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_334[k] = f_14 * skg_164[k]
                   + f_3 * pc_y[k] * slg_239[k];

        t_335[k] = f_1 * slf0_159[k]
                   - f_2 * slf1_159[k]
                   + f_3 * pc_z[k] * slg_239[k];

        t_336[k] = pb_z[k] * skh0_210[k]
                   - f_8 * pc_z[k] * skh1_210[k];

        t_337[k] = f_15 * skg_165[k]
                   + f_3 * pc_y[k] * slg_240[k];
    }

#pragma omp simd aligned(t_338, t_339, t_340, pb_z, pc_y, pc_z, skh0_213, skg_150, skg_167, \
                         skh1_213, slg_240, slg_242 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_338[k] = f_9 * skg_150[k]
                   + f_3 * pc_z[k] * slg_240[k];

        t_339[k] = pb_z[k] * skh0_213[k]
                   - f_8 * pc_z[k] * skh1_213[k];

        t_340[k] = f_15 * skg_167[k]
                   + f_3 * pc_y[k] * slg_242[k];
    }

#pragma omp simd aligned(t_341, t_342, t_343, pb_z, pc_x, pc_z, skh0_216, skg_153, skg_245, \
                         skh1_216, slf0_165, slf1_165, slg_243, \
                         slg_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_341[k] = f_11 * skg_245[k]
                   + f_4 * slf0_165[k]
                   - f_5 * slf1_165[k]
                   + f_3 * pc_x[k] * slg_245[k];

        t_342[k] = pb_z[k] * skh0_216[k]
                   - f_8 * pc_z[k] * skh1_216[k];

        t_343[k] = f_9 * skg_153[k]
                   + f_3 * pc_z[k] * slg_243[k];
    }

#pragma omp simd aligned(t_344, t_345, t_346, t_347, pc_x, pc_y, skg_170, skg_249, skg_250, \
                         skg_251, slf0_169, slf1_169, slg_245, slg_249, slg_250, \
                         slg_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_344[k] = f_15 * skg_170[k]
                   + f_3 * pc_y[k] * slg_245[k];

        t_345[k] = f_11 * skg_249[k]
                   + f_6 * slf0_169[k]
                   - f_7 * slf1_169[k]
                   + f_3 * pc_x[k] * slg_249[k];

        t_346[k] = f_11 * skg_250[k]
                   + f_3 * pc_x[k] * slg_250[k];

        t_347[k] = f_11 * skg_251[k]
                   + f_3 * pc_x[k] * slg_251[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, t_351, pb_z, pc_x, pc_z, skh0_225, skg_252, \
                         skg_253, skg_254, skh1_225, slg_252, slg_253, \
                         slg_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = f_11 * skg_252[k]
                   + f_3 * pc_x[k] * slg_252[k];

        t_349[k] = f_11 * skg_253[k]
                   + f_3 * pc_x[k] * slg_253[k];

        t_350[k] = f_11 * skg_254[k]
                   + f_3 * pc_x[k] * slg_254[k];

        t_351[k] = pb_z[k] * skh0_225[k]
                   - f_8 * pc_z[k] * skh1_225[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, pc_y, pc_z, skg_160, skg_177, skg_178, slf0_168, \
                         slf0_169, slf1_168, slf1_169, slg_250, slg_252, \
                         slg_253 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = f_9 * skg_160[k]
                   + f_3 * pc_z[k] * slg_250[k];

        t_353[k] = f_15 * skg_177[k]
                   + f_4 * slf0_168[k]
                   - f_5 * slf1_168[k]
                   + f_3 * pc_y[k] * slg_252[k];

        t_354[k] = f_15 * skg_178[k]
                   + f_6 * slf0_169[k]
                   - f_7 * slf1_169[k]
                   + f_3 * pc_y[k] * slg_253[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, pc_x, pc_y, pc_z, skg_164, skg_179, skg_255, \
                         slf0_169, slf0_170, slf1_169, slf1_170, slg_254, \
                         slg_255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = f_15 * skg_179[k]
                   + f_3 * pc_y[k] * slg_254[k];

        t_356[k] = f_9 * skg_164[k]
                   + f_1 * slf0_169[k]
                   - f_2 * slf1_169[k]
                   + f_3 * pc_z[k] * slg_254[k];

        t_357[k] = f_11 * skg_255[k]
                   + f_1 * slf0_170[k]
                   - f_2 * slf1_170[k]
                   + f_3 * pc_x[k] * slg_255[k];
    }
}

static auto
compute_prim_slh_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t skh0,
                                                          const size_t skg, const size_t skh1,
                                                          const size_t slf0, const size_t slf1,
                                                          const size_t slg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_13 = 3.0 / q;
    const auto f_14 = 2.5 / q;
    const auto f_15 = 2.0 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *skh0_294 = buffer.data(skh0 + 294);
    const auto *skh0_297 = buffer.data(skh0 + 297);
    const auto *skh0_299 = buffer.data(skh0 + 299);
    const auto *skh0_300 = buffer.data(skh0 + 300);
    const auto *skh0_303 = buffer.data(skh0 + 303);
    const auto *skh0_314 = buffer.data(skh0 + 314);
    const auto *skh0_315 = buffer.data(skh0 + 315);
    const auto *skh0_318 = buffer.data(skh0 + 318);
    const auto *skh0_321 = buffer.data(skh0 + 321);

    const auto *skg_165 = buffer.data(skg + 165);
    const auto *skg_168 = buffer.data(skg + 168);
    const auto *skg_175 = buffer.data(skg + 175);
    const auto *skg_179 = buffer.data(skg + 179);
    const auto *skg_180 = buffer.data(skg + 180);
    const auto *skg_182 = buffer.data(skg + 182);
    const auto *skg_183 = buffer.data(skg + 183);
    const auto *skg_185 = buffer.data(skg + 185);
    const auto *skg_190 = buffer.data(skg + 190);
    const auto *skg_192 = buffer.data(skg + 192);
    const auto *skg_193 = buffer.data(skg + 193);
    const auto *skg_194 = buffer.data(skg + 194);
    const auto *skg_195 = buffer.data(skg + 195);
    const auto *skg_197 = buffer.data(skg + 197);
    const auto *skg_198 = buffer.data(skg + 198);
    const auto *skg_200 = buffer.data(skg + 200);
    const auto *skg_205 = buffer.data(skg + 205);
    const auto *skg_207 = buffer.data(skg + 207);
    const auto *skg_208 = buffer.data(skg + 208);
    const auto *skg_209 = buffer.data(skg + 209);
    const auto *skg_210 = buffer.data(skg + 210);
    const auto *skg_211 = buffer.data(skg + 211);
    const auto *skg_212 = buffer.data(skg + 212);
    const auto *skg_213 = buffer.data(skg + 213);
    const auto *skg_215 = buffer.data(skg + 215);
    const auto *skg_220 = buffer.data(skg + 220);
    const auto *skg_222 = buffer.data(skg + 222);
    const auto *skg_223 = buffer.data(skg + 223);
    const auto *skg_224 = buffer.data(skg + 224);
    const auto *skg_225 = buffer.data(skg + 225);
    const auto *skg_227 = buffer.data(skg + 227);
    const auto *skg_228 = buffer.data(skg + 228);
    const auto *skg_230 = buffer.data(skg + 230);
    const auto *skg_235 = buffer.data(skg + 235);
    const auto *skg_237 = buffer.data(skg + 237);
    const auto *skg_238 = buffer.data(skg + 238);
    const auto *skg_239 = buffer.data(skg + 239);
    const auto *skg_240 = buffer.data(skg + 240);
    const auto *skg_242 = buffer.data(skg + 242);
    const auto *skg_245 = buffer.data(skg + 245);
    const auto *skg_258 = buffer.data(skg + 258);
    const auto *skg_260 = buffer.data(skg + 260);
    const auto *skg_261 = buffer.data(skg + 261);
    const auto *skg_264 = buffer.data(skg + 264);
    const auto *skg_265 = buffer.data(skg + 265);
    const auto *skg_266 = buffer.data(skg + 266);
    const auto *skg_267 = buffer.data(skg + 267);
    const auto *skg_268 = buffer.data(skg + 268);
    const auto *skg_269 = buffer.data(skg + 269);
    const auto *skg_270 = buffer.data(skg + 270);
    const auto *skg_273 = buffer.data(skg + 273);
    const auto *skg_275 = buffer.data(skg + 275);
    const auto *skg_276 = buffer.data(skg + 276);
    const auto *skg_279 = buffer.data(skg + 279);
    const auto *skg_280 = buffer.data(skg + 280);
    const auto *skg_281 = buffer.data(skg + 281);
    const auto *skg_282 = buffer.data(skg + 282);
    const auto *skg_283 = buffer.data(skg + 283);
    const auto *skg_284 = buffer.data(skg + 284);
    const auto *skg_295 = buffer.data(skg + 295);
    const auto *skg_296 = buffer.data(skg + 296);
    const auto *skg_297 = buffer.data(skg + 297);
    const auto *skg_298 = buffer.data(skg + 298);
    const auto *skg_299 = buffer.data(skg + 299);
    const auto *skg_300 = buffer.data(skg + 300);
    const auto *skg_303 = buffer.data(skg + 303);
    const auto *skg_305 = buffer.data(skg + 305);
    const auto *skg_306 = buffer.data(skg + 306);
    const auto *skg_309 = buffer.data(skg + 309);
    const auto *skg_310 = buffer.data(skg + 310);
    const auto *skg_311 = buffer.data(skg + 311);
    const auto *skg_312 = buffer.data(skg + 312);
    const auto *skg_313 = buffer.data(skg + 313);
    const auto *skg_314 = buffer.data(skg + 314);
    const auto *skg_315 = buffer.data(skg + 315);
    const auto *skg_318 = buffer.data(skg + 318);
    const auto *skg_320 = buffer.data(skg + 320);
    const auto *skg_321 = buffer.data(skg + 321);
    const auto *skg_324 = buffer.data(skg + 324);
    const auto *skg_325 = buffer.data(skg + 325);
    const auto *skg_326 = buffer.data(skg + 326);
    const auto *skg_327 = buffer.data(skg + 327);
    const auto *skg_328 = buffer.data(skg + 328);
    const auto *skg_329 = buffer.data(skg + 329);
    const auto *skg_335 = buffer.data(skg + 335);
    const auto *skg_339 = buffer.data(skg + 339);
    const auto *skg_340 = buffer.data(skg + 340);
    const auto *skg_341 = buffer.data(skg + 341);

    const auto *skh1_294 = buffer.data(skh1 + 294);
    const auto *skh1_297 = buffer.data(skh1 + 297);
    const auto *skh1_299 = buffer.data(skh1 + 299);
    const auto *skh1_300 = buffer.data(skh1 + 300);
    const auto *skh1_303 = buffer.data(skh1 + 303);
    const auto *skh1_314 = buffer.data(skh1 + 314);
    const auto *skh1_315 = buffer.data(skh1 + 315);
    const auto *skh1_318 = buffer.data(skh1 + 318);
    const auto *skh1_321 = buffer.data(skh1 + 321);

    const auto *slf0_173 = buffer.data(slf0 + 173);
    const auto *slf0_175 = buffer.data(slf0 + 175);
    const auto *slf0_176 = buffer.data(slf0 + 176);
    const auto *slf0_178 = buffer.data(slf0 + 178);
    const auto *slf0_179 = buffer.data(slf0 + 179);
    const auto *slf0_180 = buffer.data(slf0 + 180);
    const auto *slf0_183 = buffer.data(slf0 + 183);
    const auto *slf0_185 = buffer.data(slf0 + 185);
    const auto *slf0_186 = buffer.data(slf0 + 186);
    const auto *slf0_188 = buffer.data(slf0 + 188);
    const auto *slf0_189 = buffer.data(slf0 + 189);
    const auto *slf0_196 = buffer.data(slf0 + 196);
    const auto *slf0_198 = buffer.data(slf0 + 198);
    const auto *slf0_199 = buffer.data(slf0 + 199);
    const auto *slf0_200 = buffer.data(slf0 + 200);
    const auto *slf0_203 = buffer.data(slf0 + 203);
    const auto *slf0_205 = buffer.data(slf0 + 205);
    const auto *slf0_206 = buffer.data(slf0 + 206);
    const auto *slf0_208 = buffer.data(slf0 + 208);
    const auto *slf0_209 = buffer.data(slf0 + 209);
    const auto *slf0_210 = buffer.data(slf0 + 210);
    const auto *slf0_213 = buffer.data(slf0 + 213);
    const auto *slf0_215 = buffer.data(slf0 + 215);
    const auto *slf0_216 = buffer.data(slf0 + 216);
    const auto *slf0_218 = buffer.data(slf0 + 218);
    const auto *slf0_219 = buffer.data(slf0 + 219);
    const auto *slf0_225 = buffer.data(slf0 + 225);
    const auto *slf0_229 = buffer.data(slf0 + 229);

    const auto *slf1_173 = buffer.data(slf1 + 173);
    const auto *slf1_175 = buffer.data(slf1 + 175);
    const auto *slf1_176 = buffer.data(slf1 + 176);
    const auto *slf1_178 = buffer.data(slf1 + 178);
    const auto *slf1_179 = buffer.data(slf1 + 179);
    const auto *slf1_180 = buffer.data(slf1 + 180);
    const auto *slf1_183 = buffer.data(slf1 + 183);
    const auto *slf1_185 = buffer.data(slf1 + 185);
    const auto *slf1_186 = buffer.data(slf1 + 186);
    const auto *slf1_188 = buffer.data(slf1 + 188);
    const auto *slf1_189 = buffer.data(slf1 + 189);
    const auto *slf1_196 = buffer.data(slf1 + 196);
    const auto *slf1_198 = buffer.data(slf1 + 198);
    const auto *slf1_199 = buffer.data(slf1 + 199);
    const auto *slf1_200 = buffer.data(slf1 + 200);
    const auto *slf1_203 = buffer.data(slf1 + 203);
    const auto *slf1_205 = buffer.data(slf1 + 205);
    const auto *slf1_206 = buffer.data(slf1 + 206);
    const auto *slf1_208 = buffer.data(slf1 + 208);
    const auto *slf1_209 = buffer.data(slf1 + 209);
    const auto *slf1_210 = buffer.data(slf1 + 210);
    const auto *slf1_213 = buffer.data(slf1 + 213);
    const auto *slf1_215 = buffer.data(slf1 + 215);
    const auto *slf1_216 = buffer.data(slf1 + 216);
    const auto *slf1_218 = buffer.data(slf1 + 218);
    const auto *slf1_219 = buffer.data(slf1 + 219);
    const auto *slf1_225 = buffer.data(slf1 + 225);
    const auto *slf1_229 = buffer.data(slf1 + 229);

    const auto *slg_255 = buffer.data(slg + 255);
    const auto *slg_257 = buffer.data(slg + 257);
    const auto *slg_258 = buffer.data(slg + 258);
    const auto *slg_260 = buffer.data(slg + 260);
    const auto *slg_261 = buffer.data(slg + 261);
    const auto *slg_264 = buffer.data(slg + 264);
    const auto *slg_265 = buffer.data(slg + 265);
    const auto *slg_266 = buffer.data(slg + 266);
    const auto *slg_267 = buffer.data(slg + 267);
    const auto *slg_268 = buffer.data(slg + 268);
    const auto *slg_269 = buffer.data(slg + 269);
    const auto *slg_270 = buffer.data(slg + 270);
    const auto *slg_272 = buffer.data(slg + 272);
    const auto *slg_273 = buffer.data(slg + 273);
    const auto *slg_275 = buffer.data(slg + 275);
    const auto *slg_276 = buffer.data(slg + 276);
    const auto *slg_279 = buffer.data(slg + 279);
    const auto *slg_280 = buffer.data(slg + 280);
    const auto *slg_281 = buffer.data(slg + 281);
    const auto *slg_282 = buffer.data(slg + 282);
    const auto *slg_283 = buffer.data(slg + 283);
    const auto *slg_284 = buffer.data(slg + 284);
    const auto *slg_285 = buffer.data(slg + 285);
    const auto *slg_287 = buffer.data(slg + 287);
    const auto *slg_288 = buffer.data(slg + 288);
    const auto *slg_290 = buffer.data(slg + 290);
    const auto *slg_295 = buffer.data(slg + 295);
    const auto *slg_296 = buffer.data(slg + 296);
    const auto *slg_297 = buffer.data(slg + 297);
    const auto *slg_298 = buffer.data(slg + 298);
    const auto *slg_299 = buffer.data(slg + 299);
    const auto *slg_300 = buffer.data(slg + 300);
    const auto *slg_302 = buffer.data(slg + 302);
    const auto *slg_303 = buffer.data(slg + 303);
    const auto *slg_305 = buffer.data(slg + 305);
    const auto *slg_306 = buffer.data(slg + 306);
    const auto *slg_309 = buffer.data(slg + 309);
    const auto *slg_310 = buffer.data(slg + 310);
    const auto *slg_311 = buffer.data(slg + 311);
    const auto *slg_312 = buffer.data(slg + 312);
    const auto *slg_313 = buffer.data(slg + 313);
    const auto *slg_314 = buffer.data(slg + 314);
    const auto *slg_315 = buffer.data(slg + 315);
    const auto *slg_317 = buffer.data(slg + 317);
    const auto *slg_318 = buffer.data(slg + 318);
    const auto *slg_320 = buffer.data(slg + 320);
    const auto *slg_321 = buffer.data(slg + 321);
    const auto *slg_324 = buffer.data(slg + 324);
    const auto *slg_325 = buffer.data(slg + 325);
    const auto *slg_326 = buffer.data(slg + 326);
    const auto *slg_327 = buffer.data(slg + 327);
    const auto *slg_328 = buffer.data(slg + 328);
    const auto *slg_329 = buffer.data(slg + 329);
    const auto *slg_330 = buffer.data(slg + 330);
    const auto *slg_332 = buffer.data(slg + 332);
    const auto *slg_333 = buffer.data(slg + 333);
    const auto *slg_335 = buffer.data(slg + 335);
    const auto *slg_339 = buffer.data(slg + 339);
    const auto *slg_340 = buffer.data(slg + 340);
    const auto *slg_341 = buffer.data(slg + 341);

#pragma omp simd aligned(t_358, t_359, t_360, t_361, pc_x, pc_y, pc_z, skg_165, skg_180, \
                         skg_182, skg_258, slf0_173, slf1_173, slg_255, slg_257, \
                         slg_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_358[k] = f_11 * skg_180[k]
                   + f_3 * pc_y[k] * slg_255[k];

        t_359[k] = f_10 * skg_165[k]
                   + f_3 * pc_z[k] * slg_255[k];

        t_360[k] = f_11 * skg_258[k]
                   + f_4 * slf0_173[k]
                   - f_5 * slf1_173[k]
                   + f_3 * pc_x[k] * slg_258[k];

        t_361[k] = f_11 * skg_182[k]
                   + f_3 * pc_y[k] * slg_257[k];
    }

#pragma omp simd aligned(t_362, t_363, t_364, pc_x, pc_z, skg_168, skg_260, skg_261, slf0_175, \
                         slf0_176, slf1_175, slf1_176, slg_258, slg_260, \
                         slg_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = f_11 * skg_260[k]
                   + f_4 * slf0_175[k]
                   - f_5 * slf1_175[k]
                   + f_3 * pc_x[k] * slg_260[k];

        t_363[k] = f_11 * skg_261[k]
                   + f_6 * slf0_176[k]
                   - f_7 * slf1_176[k]
                   + f_3 * pc_x[k] * slg_261[k];

        t_364[k] = f_10 * skg_168[k]
                   + f_3 * pc_z[k] * slg_258[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, pc_x, pc_y, skg_185, skg_264, skg_265, \
                         skg_266, slf0_179, slf1_179, slg_260, slg_264, slg_265, \
                         slg_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = f_11 * skg_185[k]
                   + f_3 * pc_y[k] * slg_260[k];

        t_366[k] = f_11 * skg_264[k]
                   + f_6 * slf0_179[k]
                   - f_7 * slf1_179[k]
                   + f_3 * pc_x[k] * slg_264[k];

        t_367[k] = f_11 * skg_265[k]
                   + f_3 * pc_x[k] * slg_265[k];

        t_368[k] = f_11 * skg_266[k]
                   + f_3 * pc_x[k] * slg_266[k];
    }

#pragma omp simd aligned(t_369, t_370, t_371, t_372, pc_x, pc_y, skg_190, skg_267, skg_268, \
                         skg_269, slf0_176, slf1_176, slg_265, slg_267, slg_268, \
                         slg_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_369[k] = f_11 * skg_267[k]
                   + f_3 * pc_x[k] * slg_267[k];

        t_370[k] = f_11 * skg_268[k]
                   + f_3 * pc_x[k] * slg_268[k];

        t_371[k] = f_11 * skg_269[k]
                   + f_3 * pc_x[k] * slg_269[k];

        t_372[k] = f_11 * skg_190[k]
                   + f_1 * slf0_176[k]
                   - f_2 * slf1_176[k]
                   + f_3 * pc_y[k] * slg_265[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, pc_y, pc_z, skg_175, skg_192, skg_193, slf0_178, \
                         slf0_179, slf1_178, slf1_179, slg_265, slg_267, \
                         slg_268 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = f_10 * skg_175[k]
                   + f_3 * pc_z[k] * slg_265[k];

        t_374[k] = f_11 * skg_192[k]
                   + f_4 * slf0_178[k]
                   - f_5 * slf1_178[k]
                   + f_3 * pc_y[k] * slg_267[k];

        t_375[k] = f_11 * skg_193[k]
                   + f_6 * slf0_179[k]
                   - f_7 * slf1_179[k]
                   + f_3 * pc_y[k] * slg_268[k];
    }

#pragma omp simd aligned(t_376, t_377, t_378, pc_x, pc_y, pc_z, skg_179, skg_194, skg_270, \
                         slf0_179, slf0_180, slf1_179, slf1_180, slg_269, \
                         slg_270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_376[k] = f_11 * skg_194[k]
                   + f_3 * pc_y[k] * slg_269[k];

        t_377[k] = f_10 * skg_179[k]
                   + f_1 * slf0_179[k]
                   - f_2 * slf1_179[k]
                   + f_3 * pc_z[k] * slg_269[k];

        t_378[k] = f_11 * skg_270[k]
                   + f_1 * slf0_180[k]
                   - f_2 * slf1_180[k]
                   + f_3 * pc_x[k] * slg_270[k];
    }

#pragma omp simd aligned(t_379, t_380, t_381, t_382, pc_x, pc_y, pc_z, skg_180, skg_195, \
                         skg_197, skg_273, slf0_183, slf1_183, slg_270, slg_272, \
                         slg_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_379[k] = f_10 * skg_195[k]
                   + f_3 * pc_y[k] * slg_270[k];

        t_380[k] = f_11 * skg_180[k]
                   + f_3 * pc_z[k] * slg_270[k];

        t_381[k] = f_11 * skg_273[k]
                   + f_4 * slf0_183[k]
                   - f_5 * slf1_183[k]
                   + f_3 * pc_x[k] * slg_273[k];

        t_382[k] = f_10 * skg_197[k]
                   + f_3 * pc_y[k] * slg_272[k];
    }

#pragma omp simd aligned(t_383, t_384, t_385, pc_x, pc_z, skg_183, skg_275, skg_276, slf0_185, \
                         slf0_186, slf1_185, slf1_186, slg_273, slg_275, \
                         slg_276 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_383[k] = f_11 * skg_275[k]
                   + f_4 * slf0_185[k]
                   - f_5 * slf1_185[k]
                   + f_3 * pc_x[k] * slg_275[k];

        t_384[k] = f_11 * skg_276[k]
                   + f_6 * slf0_186[k]
                   - f_7 * slf1_186[k]
                   + f_3 * pc_x[k] * slg_276[k];

        t_385[k] = f_11 * skg_183[k]
                   + f_3 * pc_z[k] * slg_273[k];
    }

#pragma omp simd aligned(t_386, t_387, t_388, t_389, pc_x, pc_y, skg_200, skg_279, skg_280, \
                         skg_281, slf0_189, slf1_189, slg_275, slg_279, slg_280, \
                         slg_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_386[k] = f_10 * skg_200[k]
                   + f_3 * pc_y[k] * slg_275[k];

        t_387[k] = f_11 * skg_279[k]
                   + f_6 * slf0_189[k]
                   - f_7 * slf1_189[k]
                   + f_3 * pc_x[k] * slg_279[k];

        t_388[k] = f_11 * skg_280[k]
                   + f_3 * pc_x[k] * slg_280[k];

        t_389[k] = f_11 * skg_281[k]
                   + f_3 * pc_x[k] * slg_281[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, pc_x, pc_y, skg_205, skg_282, skg_283, \
                         skg_284, slf0_186, slf1_186, slg_280, slg_282, slg_283, \
                         slg_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = f_11 * skg_282[k]
                   + f_3 * pc_x[k] * slg_282[k];

        t_391[k] = f_11 * skg_283[k]
                   + f_3 * pc_x[k] * slg_283[k];

        t_392[k] = f_11 * skg_284[k]
                   + f_3 * pc_x[k] * slg_284[k];

        t_393[k] = f_10 * skg_205[k]
                   + f_1 * slf0_186[k]
                   - f_2 * slf1_186[k]
                   + f_3 * pc_y[k] * slg_280[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, pc_y, pc_z, skg_190, skg_207, skg_208, slf0_188, \
                         slf0_189, slf1_188, slf1_189, slg_280, slg_282, \
                         slg_283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = f_11 * skg_190[k]
                   + f_3 * pc_z[k] * slg_280[k];

        t_395[k] = f_10 * skg_207[k]
                   + f_4 * slf0_188[k]
                   - f_5 * slf1_188[k]
                   + f_3 * pc_y[k] * slg_282[k];

        t_396[k] = f_10 * skg_208[k]
                   + f_6 * slf0_189[k]
                   - f_7 * slf1_189[k]
                   + f_3 * pc_y[k] * slg_283[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, t_400, pb_y, pc_y, pc_z, skh0_294, skg_194, \
                         skg_209, skg_210, skh1_294, slf0_189, slf1_189, slg_284, \
                         slg_285 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = f_10 * skg_209[k]
                   + f_3 * pc_y[k] * slg_284[k];

        t_398[k] = f_11 * skg_194[k]
                   + f_1 * slf0_189[k]
                   - f_2 * slf1_189[k]
                   + f_3 * pc_z[k] * slg_284[k];

        t_399[k] = pb_y[k] * skh0_294[k]
                   - f_8 * pc_y[k] * skh1_294[k];

        t_400[k] = f_9 * skg_210[k]
                   + f_3 * pc_y[k] * slg_285[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, t_404, pb_y, pc_y, pc_z, skh0_297, skh0_299, \
                         skg_195, skg_211, skg_212, skh1_297, skh1_299, slg_285, \
                         slg_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_15 * skg_195[k]
                   + f_3 * pc_z[k] * slg_285[k];

        t_402[k] = pb_y[k] * skh0_297[k]
                   + f_10 * skg_211[k]
                   - f_8 * pc_y[k] * skh1_297[k];

        t_403[k] = f_9 * skg_212[k]
                   + f_3 * pc_y[k] * slg_287[k];

        t_404[k] = pb_y[k] * skh0_299[k]
                   - f_8 * pc_y[k] * skh1_299[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, pb_y, pc_y, pc_z, skh0_300, skh0_303, \
                         skg_198, skg_213, skg_215, skh1_300, skh1_303, slg_288, \
                         slg_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = pb_y[k] * skh0_300[k]
                   + f_11 * skg_213[k]
                   - f_8 * pc_y[k] * skh1_300[k];

        t_406[k] = f_15 * skg_198[k]
                   + f_3 * pc_z[k] * slg_288[k];

        t_407[k] = f_9 * skg_215[k]
                   + f_3 * pc_y[k] * slg_290[k];

        t_408[k] = pb_y[k] * skh0_303[k]
                   - f_8 * pc_y[k] * skh1_303[k];
    }

#pragma omp simd aligned(t_409, t_410, t_411, t_412, t_413, pc_x, skg_295, skg_296, skg_297, \
                         skg_298, skg_299, slg_295, slg_296, slg_297, slg_298, \
                         slg_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_409[k] = f_11 * skg_295[k]
                   + f_3 * pc_x[k] * slg_295[k];

        t_410[k] = f_11 * skg_296[k]
                   + f_3 * pc_x[k] * slg_296[k];

        t_411[k] = f_11 * skg_297[k]
                   + f_3 * pc_x[k] * slg_297[k];

        t_412[k] = f_11 * skg_298[k]
                   + f_3 * pc_x[k] * slg_298[k];

        t_413[k] = f_11 * skg_299[k]
                   + f_3 * pc_x[k] * slg_299[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, pc_y, pc_z, skg_205, skg_220, skg_222, slf0_196, \
                         slf0_198, slf1_196, slf1_198, slg_295, \
                         slg_297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_9 * skg_220[k]
                   + f_1 * slf0_196[k]
                   - f_2 * slf1_196[k]
                   + f_3 * pc_y[k] * slg_295[k];

        t_415[k] = f_15 * skg_205[k]
                   + f_3 * pc_z[k] * slg_295[k];

        t_416[k] = f_9 * skg_222[k]
                   + f_4 * slf0_198[k]
                   - f_5 * slf1_198[k]
                   + f_3 * pc_y[k] * slg_297[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, pb_y, pc_y, skh0_314, skg_223, skg_224, \
                         skh1_314, slf0_199, slf1_199, slg_298, \
                         slg_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = f_9 * skg_223[k]
                   + f_6 * slf0_199[k]
                   - f_7 * slf1_199[k]
                   + f_3 * pc_y[k] * slg_298[k];

        t_418[k] = f_9 * skg_224[k]
                   + f_3 * pc_y[k] * slg_299[k];

        t_419[k] = pb_y[k] * skh0_314[k]
                   - f_8 * pc_y[k] * skh1_314[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, pc_x, pc_y, pc_z, skg_210, skg_300, \
                         skg_303, slf0_200, slf0_203, slf1_200, slf1_203, slg_300, \
                         slg_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_11 * skg_300[k]
                   + f_1 * slf0_200[k]
                   - f_2 * slf1_200[k]
                   + f_3 * pc_x[k] * slg_300[k];

        t_421[k] = f_3 * pc_y[k] * slg_300[k];

        t_422[k] = f_14 * skg_210[k]
                   + f_3 * pc_z[k] * slg_300[k];

        t_423[k] = f_11 * skg_303[k]
                   + f_4 * slf0_203[k]
                   - f_5 * slf1_203[k]
                   + f_3 * pc_x[k] * slg_303[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, pc_x, pc_y, skg_305, skg_306, slf0_205, \
                         slf0_206, slf1_205, slf1_206, slg_302, slg_305, \
                         slg_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = f_3 * pc_y[k] * slg_302[k];

        t_425[k] = f_11 * skg_305[k]
                   + f_4 * slf0_205[k]
                   - f_5 * slf1_205[k]
                   + f_3 * pc_x[k] * slg_305[k];

        t_426[k] = f_11 * skg_306[k]
                   + f_6 * slf0_206[k]
                   - f_7 * slf1_206[k]
                   + f_3 * pc_x[k] * slg_306[k];
    }

#pragma omp simd aligned(t_427, t_428, t_429, t_430, pc_x, pc_y, pc_z, skg_213, skg_309, \
                         skg_310, slf0_209, slf1_209, slg_303, slg_305, slg_309, \
                         slg_310 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_427[k] = f_14 * skg_213[k]
                   + f_3 * pc_z[k] * slg_303[k];

        t_428[k] = f_3 * pc_y[k] * slg_305[k];

        t_429[k] = f_11 * skg_309[k]
                   + f_6 * slf0_209[k]
                   - f_7 * slf1_209[k]
                   + f_3 * pc_x[k] * slg_309[k];

        t_430[k] = f_11 * skg_310[k]
                   + f_3 * pc_x[k] * slg_310[k];
    }

#pragma omp simd aligned(t_431, t_432, t_433, t_434, pc_x, skg_311, skg_312, skg_313, skg_314, \
                         slg_311, slg_312, slg_313, slg_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_431[k] = f_11 * skg_311[k]
                   + f_3 * pc_x[k] * slg_311[k];

        t_432[k] = f_11 * skg_312[k]
                   + f_3 * pc_x[k] * slg_312[k];

        t_433[k] = f_11 * skg_313[k]
                   + f_3 * pc_x[k] * slg_313[k];

        t_434[k] = f_11 * skg_314[k]
                   + f_3 * pc_x[k] * slg_314[k];
    }

#pragma omp simd aligned(t_435, t_436, t_437, t_438, pc_y, pc_z, skg_220, slf0_206, slf0_208, \
                         slf0_209, slf1_206, slf1_208, slf1_209, slg_310, slg_312, \
                         slg_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_435[k] = f_1 * slf0_206[k]
                   - f_2 * slf1_206[k]
                   + f_3 * pc_y[k] * slg_310[k];

        t_436[k] = f_14 * skg_220[k]
                   + f_3 * pc_z[k] * slg_310[k];

        t_437[k] = f_4 * slf0_208[k]
                   - f_5 * slf1_208[k]
                   + f_3 * pc_y[k] * slg_312[k];

        t_438[k] = f_6 * slf0_209[k]
                   - f_7 * slf1_209[k]
                   + f_3 * pc_y[k] * slg_313[k];
    }

#pragma omp simd aligned(t_439, t_440, t_441, t_442, pc_x, pc_y, pc_z, skg_224, skg_225, \
                         skg_315, slf0_209, slf0_210, slf1_209, slf1_210, slg_314, \
                         slg_315 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_439[k] = f_3 * pc_y[k] * slg_314[k];

        t_440[k] = f_14 * skg_224[k]
                   + f_1 * slf0_209[k]
                   - f_2 * slf1_209[k]
                   + f_3 * pc_z[k] * slg_314[k];

        t_441[k] = f_10 * skg_315[k]
                   + f_1 * slf0_210[k]
                   - f_2 * slf1_210[k]
                   + f_3 * pc_x[k] * slg_315[k];

        t_442[k] = f_13 * skg_225[k]
                   + f_3 * pc_y[k] * slg_315[k];
    }

#pragma omp simd aligned(t_443, t_444, t_445, pc_x, pc_y, pc_z, skg_227, skg_318, slf0_213, \
                         slf1_213, slg_315, slg_317, slg_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_443[k] = f_3 * pc_z[k] * slg_315[k];

        t_444[k] = f_10 * skg_318[k]
                   + f_4 * slf0_213[k]
                   - f_5 * slf1_213[k]
                   + f_3 * pc_x[k] * slg_318[k];

        t_445[k] = f_13 * skg_227[k]
                   + f_3 * pc_y[k] * slg_317[k];
    }

#pragma omp simd aligned(t_446, t_447, t_448, pc_x, pc_z, skg_320, skg_321, slf0_215, \
                         slf0_216, slf1_215, slf1_216, slg_318, slg_320, \
                         slg_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_446[k] = f_10 * skg_320[k]
                   + f_4 * slf0_215[k]
                   - f_5 * slf1_215[k]
                   + f_3 * pc_x[k] * slg_320[k];

        t_447[k] = f_10 * skg_321[k]
                   + f_6 * slf0_216[k]
                   - f_7 * slf1_216[k]
                   + f_3 * pc_x[k] * slg_321[k];

        t_448[k] = f_3 * pc_z[k] * slg_318[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, t_452, pc_x, pc_y, skg_230, skg_324, skg_325, \
                         skg_326, slf0_219, slf1_219, slg_320, slg_324, slg_325, \
                         slg_326 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_13 * skg_230[k]
                   + f_3 * pc_y[k] * slg_320[k];

        t_450[k] = f_10 * skg_324[k]
                   + f_6 * slf0_219[k]
                   - f_7 * slf1_219[k]
                   + f_3 * pc_x[k] * slg_324[k];

        t_451[k] = f_10 * skg_325[k]
                   + f_3 * pc_x[k] * slg_325[k];

        t_452[k] = f_10 * skg_326[k]
                   + f_3 * pc_x[k] * slg_326[k];
    }

#pragma omp simd aligned(t_453, t_454, t_455, t_456, pc_x, pc_y, skg_235, skg_327, skg_328, \
                         skg_329, slf0_216, slf1_216, slg_325, slg_327, slg_328, \
                         slg_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_453[k] = f_10 * skg_327[k]
                   + f_3 * pc_x[k] * slg_327[k];

        t_454[k] = f_10 * skg_328[k]
                   + f_3 * pc_x[k] * slg_328[k];

        t_455[k] = f_10 * skg_329[k]
                   + f_3 * pc_x[k] * slg_329[k];

        t_456[k] = f_13 * skg_235[k]
                   + f_1 * slf0_216[k]
                   - f_2 * slf1_216[k]
                   + f_3 * pc_y[k] * slg_325[k];
    }

#pragma omp simd aligned(t_457, t_458, t_459, pc_y, pc_z, skg_237, skg_238, slf0_218, \
                         slf0_219, slf1_218, slf1_219, slg_325, slg_327, \
                         slg_328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_457[k] = f_3 * pc_z[k] * slg_325[k];

        t_458[k] = f_13 * skg_237[k]
                   + f_4 * slf0_218[k]
                   - f_5 * slf1_218[k]
                   + f_3 * pc_y[k] * slg_327[k];

        t_459[k] = f_13 * skg_238[k]
                   + f_6 * slf0_219[k]
                   - f_7 * slf1_219[k]
                   + f_3 * pc_y[k] * slg_328[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, t_463, pb_z, pc_y, pc_z, skh0_315, skg_239, \
                         skg_240, skh1_315, slf0_219, slf1_219, slg_329, \
                         slg_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = f_13 * skg_239[k]
                   + f_3 * pc_y[k] * slg_329[k];

        t_461[k] = f_1 * slf0_219[k]
                   - f_2 * slf1_219[k]
                   + f_3 * pc_z[k] * slg_329[k];

        t_462[k] = pb_z[k] * skh0_315[k]
                   - f_8 * pc_z[k] * skh1_315[k];

        t_463[k] = f_14 * skg_240[k]
                   + f_3 * pc_y[k] * slg_330[k];
    }

#pragma omp simd aligned(t_464, t_465, t_466, pb_z, pc_y, pc_z, skh0_318, skg_225, skg_242, \
                         skh1_318, slg_330, slg_332 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_464[k] = f_9 * skg_225[k]
                   + f_3 * pc_z[k] * slg_330[k];

        t_465[k] = pb_z[k] * skh0_318[k]
                   - f_8 * pc_z[k] * skh1_318[k];

        t_466[k] = f_14 * skg_242[k]
                   + f_3 * pc_y[k] * slg_332[k];
    }

#pragma omp simd aligned(t_467, t_468, t_469, pb_z, pc_x, pc_z, skh0_321, skg_228, skg_335, \
                         skh1_321, slf0_225, slf1_225, slg_333, \
                         slg_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_467[k] = f_10 * skg_335[k]
                   + f_4 * slf0_225[k]
                   - f_5 * slf1_225[k]
                   + f_3 * pc_x[k] * slg_335[k];

        t_468[k] = pb_z[k] * skh0_321[k]
                   - f_8 * pc_z[k] * skh1_321[k];

        t_469[k] = f_9 * skg_228[k]
                   + f_3 * pc_z[k] * slg_333[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, pc_x, pc_y, skg_245, skg_339, skg_340, \
                         skg_341, slf0_229, slf1_229, slg_335, slg_339, slg_340, \
                         slg_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = f_14 * skg_245[k]
                   + f_3 * pc_y[k] * slg_335[k];

        t_471[k] = f_10 * skg_339[k]
                   + f_6 * slf0_229[k]
                   - f_7 * slf1_229[k]
                   + f_3 * pc_x[k] * slg_339[k];

        t_472[k] = f_10 * skg_340[k]
                   + f_3 * pc_x[k] * slg_340[k];

        t_473[k] = f_10 * skg_341[k]
                   + f_3 * pc_x[k] * slg_341[k];
    }
}

static auto
compute_prim_slh_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t skh0,
                                                          const size_t skg, const size_t skh1,
                                                          const size_t slf0, const size_t slf1,
                                                          const size_t slg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_13 = 3.0 / q;
    const auto f_14 = 2.5 / q;
    const auto f_15 = 2.0 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *skh0_330 = buffer.data(skh0 + 330);
    const auto *skh0_420 = buffer.data(skh0 + 420);
    const auto *skh0_423 = buffer.data(skh0 + 423);
    const auto *skh0_425 = buffer.data(skh0 + 425);
    const auto *skh0_426 = buffer.data(skh0 + 426);
    const auto *skh0_429 = buffer.data(skh0 + 429);
    const auto *skh0_440 = buffer.data(skh0 + 440);

    const auto *skg_235 = buffer.data(skg + 235);
    const auto *skg_239 = buffer.data(skg + 239);
    const auto *skg_240 = buffer.data(skg + 240);
    const auto *skg_243 = buffer.data(skg + 243);
    const auto *skg_250 = buffer.data(skg + 250);
    const auto *skg_252 = buffer.data(skg + 252);
    const auto *skg_253 = buffer.data(skg + 253);
    const auto *skg_254 = buffer.data(skg + 254);
    const auto *skg_255 = buffer.data(skg + 255);
    const auto *skg_257 = buffer.data(skg + 257);
    const auto *skg_258 = buffer.data(skg + 258);
    const auto *skg_260 = buffer.data(skg + 260);
    const auto *skg_265 = buffer.data(skg + 265);
    const auto *skg_267 = buffer.data(skg + 267);
    const auto *skg_268 = buffer.data(skg + 268);
    const auto *skg_269 = buffer.data(skg + 269);
    const auto *skg_270 = buffer.data(skg + 270);
    const auto *skg_272 = buffer.data(skg + 272);
    const auto *skg_273 = buffer.data(skg + 273);
    const auto *skg_275 = buffer.data(skg + 275);
    const auto *skg_280 = buffer.data(skg + 280);
    const auto *skg_282 = buffer.data(skg + 282);
    const auto *skg_283 = buffer.data(skg + 283);
    const auto *skg_284 = buffer.data(skg + 284);
    const auto *skg_285 = buffer.data(skg + 285);
    const auto *skg_287 = buffer.data(skg + 287);
    const auto *skg_288 = buffer.data(skg + 288);
    const auto *skg_290 = buffer.data(skg + 290);
    const auto *skg_295 = buffer.data(skg + 295);
    const auto *skg_297 = buffer.data(skg + 297);
    const auto *skg_298 = buffer.data(skg + 298);
    const auto *skg_299 = buffer.data(skg + 299);
    const auto *skg_300 = buffer.data(skg + 300);
    const auto *skg_301 = buffer.data(skg + 301);
    const auto *skg_302 = buffer.data(skg + 302);
    const auto *skg_303 = buffer.data(skg + 303);
    const auto *skg_305 = buffer.data(skg + 305);
    const auto *skg_310 = buffer.data(skg + 310);
    const auto *skg_312 = buffer.data(skg + 312);
    const auto *skg_313 = buffer.data(skg + 313);
    const auto *skg_314 = buffer.data(skg + 314);
    const auto *skg_342 = buffer.data(skg + 342);
    const auto *skg_343 = buffer.data(skg + 343);
    const auto *skg_344 = buffer.data(skg + 344);
    const auto *skg_345 = buffer.data(skg + 345);
    const auto *skg_348 = buffer.data(skg + 348);
    const auto *skg_350 = buffer.data(skg + 350);
    const auto *skg_351 = buffer.data(skg + 351);
    const auto *skg_354 = buffer.data(skg + 354);
    const auto *skg_355 = buffer.data(skg + 355);
    const auto *skg_356 = buffer.data(skg + 356);
    const auto *skg_357 = buffer.data(skg + 357);
    const auto *skg_358 = buffer.data(skg + 358);
    const auto *skg_359 = buffer.data(skg + 359);
    const auto *skg_360 = buffer.data(skg + 360);
    const auto *skg_363 = buffer.data(skg + 363);
    const auto *skg_365 = buffer.data(skg + 365);
    const auto *skg_366 = buffer.data(skg + 366);
    const auto *skg_369 = buffer.data(skg + 369);
    const auto *skg_370 = buffer.data(skg + 370);
    const auto *skg_371 = buffer.data(skg + 371);
    const auto *skg_372 = buffer.data(skg + 372);
    const auto *skg_373 = buffer.data(skg + 373);
    const auto *skg_374 = buffer.data(skg + 374);
    const auto *skg_375 = buffer.data(skg + 375);
    const auto *skg_378 = buffer.data(skg + 378);
    const auto *skg_380 = buffer.data(skg + 380);
    const auto *skg_381 = buffer.data(skg + 381);
    const auto *skg_384 = buffer.data(skg + 384);
    const auto *skg_385 = buffer.data(skg + 385);
    const auto *skg_386 = buffer.data(skg + 386);
    const auto *skg_387 = buffer.data(skg + 387);
    const auto *skg_388 = buffer.data(skg + 388);
    const auto *skg_389 = buffer.data(skg + 389);
    const auto *skg_400 = buffer.data(skg + 400);
    const auto *skg_401 = buffer.data(skg + 401);
    const auto *skg_402 = buffer.data(skg + 402);
    const auto *skg_403 = buffer.data(skg + 403);
    const auto *skg_404 = buffer.data(skg + 404);
    const auto *skg_405 = buffer.data(skg + 405);
    const auto *skg_408 = buffer.data(skg + 408);
    const auto *skg_410 = buffer.data(skg + 410);
    const auto *skg_411 = buffer.data(skg + 411);
    const auto *skg_414 = buffer.data(skg + 414);
    const auto *skg_415 = buffer.data(skg + 415);
    const auto *skg_416 = buffer.data(skg + 416);
    const auto *skg_417 = buffer.data(skg + 417);
    const auto *skg_418 = buffer.data(skg + 418);
    const auto *skg_419 = buffer.data(skg + 419);

    const auto *skh1_330 = buffer.data(skh1 + 330);
    const auto *skh1_420 = buffer.data(skh1 + 420);
    const auto *skh1_423 = buffer.data(skh1 + 423);
    const auto *skh1_425 = buffer.data(skh1 + 425);
    const auto *skh1_426 = buffer.data(skh1 + 426);
    const auto *skh1_429 = buffer.data(skh1 + 429);
    const auto *skh1_440 = buffer.data(skh1 + 440);

    const auto *slf0_228 = buffer.data(slf0 + 228);
    const auto *slf0_229 = buffer.data(slf0 + 229);
    const auto *slf0_230 = buffer.data(slf0 + 230);
    const auto *slf0_233 = buffer.data(slf0 + 233);
    const auto *slf0_235 = buffer.data(slf0 + 235);
    const auto *slf0_236 = buffer.data(slf0 + 236);
    const auto *slf0_238 = buffer.data(slf0 + 238);
    const auto *slf0_239 = buffer.data(slf0 + 239);
    const auto *slf0_240 = buffer.data(slf0 + 240);
    const auto *slf0_243 = buffer.data(slf0 + 243);
    const auto *slf0_245 = buffer.data(slf0 + 245);
    const auto *slf0_246 = buffer.data(slf0 + 246);
    const auto *slf0_248 = buffer.data(slf0 + 248);
    const auto *slf0_249 = buffer.data(slf0 + 249);
    const auto *slf0_250 = buffer.data(slf0 + 250);
    const auto *slf0_253 = buffer.data(slf0 + 253);
    const auto *slf0_255 = buffer.data(slf0 + 255);
    const auto *slf0_256 = buffer.data(slf0 + 256);
    const auto *slf0_258 = buffer.data(slf0 + 258);
    const auto *slf0_259 = buffer.data(slf0 + 259);
    const auto *slf0_266 = buffer.data(slf0 + 266);
    const auto *slf0_268 = buffer.data(slf0 + 268);
    const auto *slf0_269 = buffer.data(slf0 + 269);
    const auto *slf0_270 = buffer.data(slf0 + 270);
    const auto *slf0_273 = buffer.data(slf0 + 273);
    const auto *slf0_275 = buffer.data(slf0 + 275);
    const auto *slf0_276 = buffer.data(slf0 + 276);
    const auto *slf0_278 = buffer.data(slf0 + 278);
    const auto *slf0_279 = buffer.data(slf0 + 279);

    const auto *slf1_228 = buffer.data(slf1 + 228);
    const auto *slf1_229 = buffer.data(slf1 + 229);
    const auto *slf1_230 = buffer.data(slf1 + 230);
    const auto *slf1_233 = buffer.data(slf1 + 233);
    const auto *slf1_235 = buffer.data(slf1 + 235);
    const auto *slf1_236 = buffer.data(slf1 + 236);
    const auto *slf1_238 = buffer.data(slf1 + 238);
    const auto *slf1_239 = buffer.data(slf1 + 239);
    const auto *slf1_240 = buffer.data(slf1 + 240);
    const auto *slf1_243 = buffer.data(slf1 + 243);
    const auto *slf1_245 = buffer.data(slf1 + 245);
    const auto *slf1_246 = buffer.data(slf1 + 246);
    const auto *slf1_248 = buffer.data(slf1 + 248);
    const auto *slf1_249 = buffer.data(slf1 + 249);
    const auto *slf1_250 = buffer.data(slf1 + 250);
    const auto *slf1_253 = buffer.data(slf1 + 253);
    const auto *slf1_255 = buffer.data(slf1 + 255);
    const auto *slf1_256 = buffer.data(slf1 + 256);
    const auto *slf1_258 = buffer.data(slf1 + 258);
    const auto *slf1_259 = buffer.data(slf1 + 259);
    const auto *slf1_266 = buffer.data(slf1 + 266);
    const auto *slf1_268 = buffer.data(slf1 + 268);
    const auto *slf1_269 = buffer.data(slf1 + 269);
    const auto *slf1_270 = buffer.data(slf1 + 270);
    const auto *slf1_273 = buffer.data(slf1 + 273);
    const auto *slf1_275 = buffer.data(slf1 + 275);
    const auto *slf1_276 = buffer.data(slf1 + 276);
    const auto *slf1_278 = buffer.data(slf1 + 278);
    const auto *slf1_279 = buffer.data(slf1 + 279);

    const auto *slg_340 = buffer.data(slg + 340);
    const auto *slg_342 = buffer.data(slg + 342);
    const auto *slg_343 = buffer.data(slg + 343);
    const auto *slg_344 = buffer.data(slg + 344);
    const auto *slg_345 = buffer.data(slg + 345);
    const auto *slg_347 = buffer.data(slg + 347);
    const auto *slg_348 = buffer.data(slg + 348);
    const auto *slg_350 = buffer.data(slg + 350);
    const auto *slg_351 = buffer.data(slg + 351);
    const auto *slg_354 = buffer.data(slg + 354);
    const auto *slg_355 = buffer.data(slg + 355);
    const auto *slg_356 = buffer.data(slg + 356);
    const auto *slg_357 = buffer.data(slg + 357);
    const auto *slg_358 = buffer.data(slg + 358);
    const auto *slg_359 = buffer.data(slg + 359);
    const auto *slg_360 = buffer.data(slg + 360);
    const auto *slg_362 = buffer.data(slg + 362);
    const auto *slg_363 = buffer.data(slg + 363);
    const auto *slg_365 = buffer.data(slg + 365);
    const auto *slg_366 = buffer.data(slg + 366);
    const auto *slg_369 = buffer.data(slg + 369);
    const auto *slg_370 = buffer.data(slg + 370);
    const auto *slg_371 = buffer.data(slg + 371);
    const auto *slg_372 = buffer.data(slg + 372);
    const auto *slg_373 = buffer.data(slg + 373);
    const auto *slg_374 = buffer.data(slg + 374);
    const auto *slg_375 = buffer.data(slg + 375);
    const auto *slg_377 = buffer.data(slg + 377);
    const auto *slg_378 = buffer.data(slg + 378);
    const auto *slg_380 = buffer.data(slg + 380);
    const auto *slg_381 = buffer.data(slg + 381);
    const auto *slg_384 = buffer.data(slg + 384);
    const auto *slg_385 = buffer.data(slg + 385);
    const auto *slg_386 = buffer.data(slg + 386);
    const auto *slg_387 = buffer.data(slg + 387);
    const auto *slg_388 = buffer.data(slg + 388);
    const auto *slg_389 = buffer.data(slg + 389);
    const auto *slg_390 = buffer.data(slg + 390);
    const auto *slg_392 = buffer.data(slg + 392);
    const auto *slg_393 = buffer.data(slg + 393);
    const auto *slg_395 = buffer.data(slg + 395);
    const auto *slg_400 = buffer.data(slg + 400);
    const auto *slg_401 = buffer.data(slg + 401);
    const auto *slg_402 = buffer.data(slg + 402);
    const auto *slg_403 = buffer.data(slg + 403);
    const auto *slg_404 = buffer.data(slg + 404);
    const auto *slg_405 = buffer.data(slg + 405);
    const auto *slg_407 = buffer.data(slg + 407);
    const auto *slg_408 = buffer.data(slg + 408);
    const auto *slg_410 = buffer.data(slg + 410);
    const auto *slg_411 = buffer.data(slg + 411);
    const auto *slg_414 = buffer.data(slg + 414);
    const auto *slg_415 = buffer.data(slg + 415);
    const auto *slg_416 = buffer.data(slg + 416);
    const auto *slg_417 = buffer.data(slg + 417);
    const auto *slg_418 = buffer.data(slg + 418);
    const auto *slg_419 = buffer.data(slg + 419);

#pragma omp simd aligned(t_474, t_475, t_476, t_477, pb_z, pc_x, pc_z, skh0_330, skg_342, \
                         skg_343, skg_344, skh1_330, slg_342, slg_343, \
                         slg_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = f_10 * skg_342[k]
                   + f_3 * pc_x[k] * slg_342[k];

        t_475[k] = f_10 * skg_343[k]
                   + f_3 * pc_x[k] * slg_343[k];

        t_476[k] = f_10 * skg_344[k]
                   + f_3 * pc_x[k] * slg_344[k];

        t_477[k] = pb_z[k] * skh0_330[k]
                   - f_8 * pc_z[k] * skh1_330[k];
    }

#pragma omp simd aligned(t_478, t_479, t_480, pc_y, pc_z, skg_235, skg_252, skg_253, slf0_228, \
                         slf0_229, slf1_228, slf1_229, slg_340, slg_342, \
                         slg_343 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_478[k] = f_9 * skg_235[k]
                   + f_3 * pc_z[k] * slg_340[k];

        t_479[k] = f_14 * skg_252[k]
                   + f_4 * slf0_228[k]
                   - f_5 * slf1_228[k]
                   + f_3 * pc_y[k] * slg_342[k];

        t_480[k] = f_14 * skg_253[k]
                   + f_6 * slf0_229[k]
                   - f_7 * slf1_229[k]
                   + f_3 * pc_y[k] * slg_343[k];
    }

#pragma omp simd aligned(t_481, t_482, t_483, pc_x, pc_y, pc_z, skg_239, skg_254, skg_345, \
                         slf0_229, slf0_230, slf1_229, slf1_230, slg_344, \
                         slg_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_481[k] = f_14 * skg_254[k]
                   + f_3 * pc_y[k] * slg_344[k];

        t_482[k] = f_9 * skg_239[k]
                   + f_1 * slf0_229[k]
                   - f_2 * slf1_229[k]
                   + f_3 * pc_z[k] * slg_344[k];

        t_483[k] = f_10 * skg_345[k]
                   + f_1 * slf0_230[k]
                   - f_2 * slf1_230[k]
                   + f_3 * pc_x[k] * slg_345[k];
    }

#pragma omp simd aligned(t_484, t_485, t_486, t_487, pc_x, pc_y, pc_z, skg_240, skg_255, \
                         skg_257, skg_348, slf0_233, slf1_233, slg_345, slg_347, \
                         slg_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_484[k] = f_15 * skg_255[k]
                   + f_3 * pc_y[k] * slg_345[k];

        t_485[k] = f_10 * skg_240[k]
                   + f_3 * pc_z[k] * slg_345[k];

        t_486[k] = f_10 * skg_348[k]
                   + f_4 * slf0_233[k]
                   - f_5 * slf1_233[k]
                   + f_3 * pc_x[k] * slg_348[k];

        t_487[k] = f_15 * skg_257[k]
                   + f_3 * pc_y[k] * slg_347[k];
    }

#pragma omp simd aligned(t_488, t_489, t_490, pc_x, pc_z, skg_243, skg_350, skg_351, slf0_235, \
                         slf0_236, slf1_235, slf1_236, slg_348, slg_350, \
                         slg_351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_488[k] = f_10 * skg_350[k]
                   + f_4 * slf0_235[k]
                   - f_5 * slf1_235[k]
                   + f_3 * pc_x[k] * slg_350[k];

        t_489[k] = f_10 * skg_351[k]
                   + f_6 * slf0_236[k]
                   - f_7 * slf1_236[k]
                   + f_3 * pc_x[k] * slg_351[k];

        t_490[k] = f_10 * skg_243[k]
                   + f_3 * pc_z[k] * slg_348[k];
    }

#pragma omp simd aligned(t_491, t_492, t_493, t_494, pc_x, pc_y, skg_260, skg_354, skg_355, \
                         skg_356, slf0_239, slf1_239, slg_350, slg_354, slg_355, \
                         slg_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_491[k] = f_15 * skg_260[k]
                   + f_3 * pc_y[k] * slg_350[k];

        t_492[k] = f_10 * skg_354[k]
                   + f_6 * slf0_239[k]
                   - f_7 * slf1_239[k]
                   + f_3 * pc_x[k] * slg_354[k];

        t_493[k] = f_10 * skg_355[k]
                   + f_3 * pc_x[k] * slg_355[k];

        t_494[k] = f_10 * skg_356[k]
                   + f_3 * pc_x[k] * slg_356[k];
    }

#pragma omp simd aligned(t_495, t_496, t_497, t_498, pc_x, pc_y, skg_265, skg_357, skg_358, \
                         skg_359, slf0_236, slf1_236, slg_355, slg_357, slg_358, \
                         slg_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_495[k] = f_10 * skg_357[k]
                   + f_3 * pc_x[k] * slg_357[k];

        t_496[k] = f_10 * skg_358[k]
                   + f_3 * pc_x[k] * slg_358[k];

        t_497[k] = f_10 * skg_359[k]
                   + f_3 * pc_x[k] * slg_359[k];

        t_498[k] = f_15 * skg_265[k]
                   + f_1 * slf0_236[k]
                   - f_2 * slf1_236[k]
                   + f_3 * pc_y[k] * slg_355[k];
    }

#pragma omp simd aligned(t_499, t_500, t_501, pc_y, pc_z, skg_250, skg_267, skg_268, slf0_238, \
                         slf0_239, slf1_238, slf1_239, slg_355, slg_357, \
                         slg_358 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_499[k] = f_10 * skg_250[k]
                   + f_3 * pc_z[k] * slg_355[k];

        t_500[k] = f_15 * skg_267[k]
                   + f_4 * slf0_238[k]
                   - f_5 * slf1_238[k]
                   + f_3 * pc_y[k] * slg_357[k];

        t_501[k] = f_15 * skg_268[k]
                   + f_6 * slf0_239[k]
                   - f_7 * slf1_239[k]
                   + f_3 * pc_y[k] * slg_358[k];
    }

#pragma omp simd aligned(t_502, t_503, t_504, pc_x, pc_y, pc_z, skg_254, skg_269, skg_360, \
                         slf0_239, slf0_240, slf1_239, slf1_240, slg_359, \
                         slg_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_502[k] = f_15 * skg_269[k]
                   + f_3 * pc_y[k] * slg_359[k];

        t_503[k] = f_10 * skg_254[k]
                   + f_1 * slf0_239[k]
                   - f_2 * slf1_239[k]
                   + f_3 * pc_z[k] * slg_359[k];

        t_504[k] = f_10 * skg_360[k]
                   + f_1 * slf0_240[k]
                   - f_2 * slf1_240[k]
                   + f_3 * pc_x[k] * slg_360[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, pc_x, pc_y, pc_z, skg_255, skg_270, \
                         skg_272, skg_363, slf0_243, slf1_243, slg_360, slg_362, \
                         slg_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = f_11 * skg_270[k]
                   + f_3 * pc_y[k] * slg_360[k];

        t_506[k] = f_11 * skg_255[k]
                   + f_3 * pc_z[k] * slg_360[k];

        t_507[k] = f_10 * skg_363[k]
                   + f_4 * slf0_243[k]
                   - f_5 * slf1_243[k]
                   + f_3 * pc_x[k] * slg_363[k];

        t_508[k] = f_11 * skg_272[k]
                   + f_3 * pc_y[k] * slg_362[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, pc_x, pc_z, skg_258, skg_365, skg_366, slf0_245, \
                         slf0_246, slf1_245, slf1_246, slg_363, slg_365, \
                         slg_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = f_10 * skg_365[k]
                   + f_4 * slf0_245[k]
                   - f_5 * slf1_245[k]
                   + f_3 * pc_x[k] * slg_365[k];

        t_510[k] = f_10 * skg_366[k]
                   + f_6 * slf0_246[k]
                   - f_7 * slf1_246[k]
                   + f_3 * pc_x[k] * slg_366[k];

        t_511[k] = f_11 * skg_258[k]
                   + f_3 * pc_z[k] * slg_363[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, t_515, pc_x, pc_y, skg_275, skg_369, skg_370, \
                         skg_371, slf0_249, slf1_249, slg_365, slg_369, slg_370, \
                         slg_371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = f_11 * skg_275[k]
                   + f_3 * pc_y[k] * slg_365[k];

        t_513[k] = f_10 * skg_369[k]
                   + f_6 * slf0_249[k]
                   - f_7 * slf1_249[k]
                   + f_3 * pc_x[k] * slg_369[k];

        t_514[k] = f_10 * skg_370[k]
                   + f_3 * pc_x[k] * slg_370[k];

        t_515[k] = f_10 * skg_371[k]
                   + f_3 * pc_x[k] * slg_371[k];
    }

#pragma omp simd aligned(t_516, t_517, t_518, t_519, pc_x, pc_y, skg_280, skg_372, skg_373, \
                         skg_374, slf0_246, slf1_246, slg_370, slg_372, slg_373, \
                         slg_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_516[k] = f_10 * skg_372[k]
                   + f_3 * pc_x[k] * slg_372[k];

        t_517[k] = f_10 * skg_373[k]
                   + f_3 * pc_x[k] * slg_373[k];

        t_518[k] = f_10 * skg_374[k]
                   + f_3 * pc_x[k] * slg_374[k];

        t_519[k] = f_11 * skg_280[k]
                   + f_1 * slf0_246[k]
                   - f_2 * slf1_246[k]
                   + f_3 * pc_y[k] * slg_370[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, pc_y, pc_z, skg_265, skg_282, skg_283, slf0_248, \
                         slf0_249, slf1_248, slf1_249, slg_370, slg_372, \
                         slg_373 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = f_11 * skg_265[k]
                   + f_3 * pc_z[k] * slg_370[k];

        t_521[k] = f_11 * skg_282[k]
                   + f_4 * slf0_248[k]
                   - f_5 * slf1_248[k]
                   + f_3 * pc_y[k] * slg_372[k];

        t_522[k] = f_11 * skg_283[k]
                   + f_6 * slf0_249[k]
                   - f_7 * slf1_249[k]
                   + f_3 * pc_y[k] * slg_373[k];
    }

#pragma omp simd aligned(t_523, t_524, t_525, pc_x, pc_y, pc_z, skg_269, skg_284, skg_375, \
                         slf0_249, slf0_250, slf1_249, slf1_250, slg_374, \
                         slg_375 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_523[k] = f_11 * skg_284[k]
                   + f_3 * pc_y[k] * slg_374[k];

        t_524[k] = f_11 * skg_269[k]
                   + f_1 * slf0_249[k]
                   - f_2 * slf1_249[k]
                   + f_3 * pc_z[k] * slg_374[k];

        t_525[k] = f_10 * skg_375[k]
                   + f_1 * slf0_250[k]
                   - f_2 * slf1_250[k]
                   + f_3 * pc_x[k] * slg_375[k];
    }

#pragma omp simd aligned(t_526, t_527, t_528, t_529, pc_x, pc_y, pc_z, skg_270, skg_285, \
                         skg_287, skg_378, slf0_253, slf1_253, slg_375, slg_377, \
                         slg_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_526[k] = f_10 * skg_285[k]
                   + f_3 * pc_y[k] * slg_375[k];

        t_527[k] = f_15 * skg_270[k]
                   + f_3 * pc_z[k] * slg_375[k];

        t_528[k] = f_10 * skg_378[k]
                   + f_4 * slf0_253[k]
                   - f_5 * slf1_253[k]
                   + f_3 * pc_x[k] * slg_378[k];

        t_529[k] = f_10 * skg_287[k]
                   + f_3 * pc_y[k] * slg_377[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, pc_x, pc_z, skg_273, skg_380, skg_381, slf0_255, \
                         slf0_256, slf1_255, slf1_256, slg_378, slg_380, \
                         slg_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = f_10 * skg_380[k]
                   + f_4 * slf0_255[k]
                   - f_5 * slf1_255[k]
                   + f_3 * pc_x[k] * slg_380[k];

        t_531[k] = f_10 * skg_381[k]
                   + f_6 * slf0_256[k]
                   - f_7 * slf1_256[k]
                   + f_3 * pc_x[k] * slg_381[k];

        t_532[k] = f_15 * skg_273[k]
                   + f_3 * pc_z[k] * slg_378[k];
    }

#pragma omp simd aligned(t_533, t_534, t_535, t_536, pc_x, pc_y, skg_290, skg_384, skg_385, \
                         skg_386, slf0_259, slf1_259, slg_380, slg_384, slg_385, \
                         slg_386 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_533[k] = f_10 * skg_290[k]
                   + f_3 * pc_y[k] * slg_380[k];

        t_534[k] = f_10 * skg_384[k]
                   + f_6 * slf0_259[k]
                   - f_7 * slf1_259[k]
                   + f_3 * pc_x[k] * slg_384[k];

        t_535[k] = f_10 * skg_385[k]
                   + f_3 * pc_x[k] * slg_385[k];

        t_536[k] = f_10 * skg_386[k]
                   + f_3 * pc_x[k] * slg_386[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, t_540, pc_x, pc_y, skg_295, skg_387, skg_388, \
                         skg_389, slf0_256, slf1_256, slg_385, slg_387, slg_388, \
                         slg_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = f_10 * skg_387[k]
                   + f_3 * pc_x[k] * slg_387[k];

        t_538[k] = f_10 * skg_388[k]
                   + f_3 * pc_x[k] * slg_388[k];

        t_539[k] = f_10 * skg_389[k]
                   + f_3 * pc_x[k] * slg_389[k];

        t_540[k] = f_10 * skg_295[k]
                   + f_1 * slf0_256[k]
                   - f_2 * slf1_256[k]
                   + f_3 * pc_y[k] * slg_385[k];
    }

#pragma omp simd aligned(t_541, t_542, t_543, pc_y, pc_z, skg_280, skg_297, skg_298, slf0_258, \
                         slf0_259, slf1_258, slf1_259, slg_385, slg_387, \
                         slg_388 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_541[k] = f_15 * skg_280[k]
                   + f_3 * pc_z[k] * slg_385[k];

        t_542[k] = f_10 * skg_297[k]
                   + f_4 * slf0_258[k]
                   - f_5 * slf1_258[k]
                   + f_3 * pc_y[k] * slg_387[k];

        t_543[k] = f_10 * skg_298[k]
                   + f_6 * slf0_259[k]
                   - f_7 * slf1_259[k]
                   + f_3 * pc_y[k] * slg_388[k];
    }

#pragma omp simd aligned(t_544, t_545, t_546, t_547, pb_y, pc_y, pc_z, skh0_420, skg_284, \
                         skg_299, skg_300, skh1_420, slf0_259, slf1_259, slg_389, \
                         slg_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_544[k] = f_10 * skg_299[k]
                   + f_3 * pc_y[k] * slg_389[k];

        t_545[k] = f_15 * skg_284[k]
                   + f_1 * slf0_259[k]
                   - f_2 * slf1_259[k]
                   + f_3 * pc_z[k] * slg_389[k];

        t_546[k] = pb_y[k] * skh0_420[k]
                   - f_8 * pc_y[k] * skh1_420[k];

        t_547[k] = f_9 * skg_300[k]
                   + f_3 * pc_y[k] * slg_390[k];
    }

#pragma omp simd aligned(t_548, t_549, t_550, t_551, pb_y, pc_y, pc_z, skh0_423, skh0_425, \
                         skg_285, skg_301, skg_302, skh1_423, skh1_425, slg_390, \
                         slg_392 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_548[k] = f_14 * skg_285[k]
                   + f_3 * pc_z[k] * slg_390[k];

        t_549[k] = pb_y[k] * skh0_423[k]
                   + f_10 * skg_301[k]
                   - f_8 * pc_y[k] * skh1_423[k];

        t_550[k] = f_9 * skg_302[k]
                   + f_3 * pc_y[k] * slg_392[k];

        t_551[k] = pb_y[k] * skh0_425[k]
                   - f_8 * pc_y[k] * skh1_425[k];
    }

#pragma omp simd aligned(t_552, t_553, t_554, t_555, pb_y, pc_y, pc_z, skh0_426, skh0_429, \
                         skg_288, skg_303, skg_305, skh1_426, skh1_429, slg_393, \
                         slg_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_552[k] = pb_y[k] * skh0_426[k]
                   + f_11 * skg_303[k]
                   - f_8 * pc_y[k] * skh1_426[k];

        t_553[k] = f_14 * skg_288[k]
                   + f_3 * pc_z[k] * slg_393[k];

        t_554[k] = f_9 * skg_305[k]
                   + f_3 * pc_y[k] * slg_395[k];

        t_555[k] = pb_y[k] * skh0_429[k]
                   - f_8 * pc_y[k] * skh1_429[k];
    }

#pragma omp simd aligned(t_556, t_557, t_558, t_559, t_560, pc_x, skg_400, skg_401, skg_402, \
                         skg_403, skg_404, slg_400, slg_401, slg_402, slg_403, \
                         slg_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_556[k] = f_10 * skg_400[k]
                   + f_3 * pc_x[k] * slg_400[k];

        t_557[k] = f_10 * skg_401[k]
                   + f_3 * pc_x[k] * slg_401[k];

        t_558[k] = f_10 * skg_402[k]
                   + f_3 * pc_x[k] * slg_402[k];

        t_559[k] = f_10 * skg_403[k]
                   + f_3 * pc_x[k] * slg_403[k];

        t_560[k] = f_10 * skg_404[k]
                   + f_3 * pc_x[k] * slg_404[k];
    }

#pragma omp simd aligned(t_561, t_562, t_563, pc_y, pc_z, skg_295, skg_310, skg_312, slf0_266, \
                         slf0_268, slf1_266, slf1_268, slg_400, \
                         slg_402 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_561[k] = f_9 * skg_310[k]
                   + f_1 * slf0_266[k]
                   - f_2 * slf1_266[k]
                   + f_3 * pc_y[k] * slg_400[k];

        t_562[k] = f_14 * skg_295[k]
                   + f_3 * pc_z[k] * slg_400[k];

        t_563[k] = f_9 * skg_312[k]
                   + f_4 * slf0_268[k]
                   - f_5 * slf1_268[k]
                   + f_3 * pc_y[k] * slg_402[k];
    }

#pragma omp simd aligned(t_564, t_565, t_566, pb_y, pc_y, skh0_440, skg_313, skg_314, \
                         skh1_440, slf0_269, slf1_269, slg_403, \
                         slg_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_564[k] = f_9 * skg_313[k]
                   + f_6 * slf0_269[k]
                   - f_7 * slf1_269[k]
                   + f_3 * pc_y[k] * slg_403[k];

        t_565[k] = f_9 * skg_314[k]
                   + f_3 * pc_y[k] * slg_404[k];

        t_566[k] = pb_y[k] * skh0_440[k]
                   - f_8 * pc_y[k] * skh1_440[k];
    }

#pragma omp simd aligned(t_567, t_568, t_569, t_570, pc_x, pc_y, pc_z, skg_300, skg_405, \
                         skg_408, slf0_270, slf0_273, slf1_270, slf1_273, slg_405, \
                         slg_408 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_567[k] = f_10 * skg_405[k]
                   + f_1 * slf0_270[k]
                   - f_2 * slf1_270[k]
                   + f_3 * pc_x[k] * slg_405[k];

        t_568[k] = f_3 * pc_y[k] * slg_405[k];

        t_569[k] = f_13 * skg_300[k]
                   + f_3 * pc_z[k] * slg_405[k];

        t_570[k] = f_10 * skg_408[k]
                   + f_4 * slf0_273[k]
                   - f_5 * slf1_273[k]
                   + f_3 * pc_x[k] * slg_408[k];
    }

#pragma omp simd aligned(t_571, t_572, t_573, pc_x, pc_y, skg_410, skg_411, slf0_275, \
                         slf0_276, slf1_275, slf1_276, slg_407, slg_410, \
                         slg_411 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_571[k] = f_3 * pc_y[k] * slg_407[k];

        t_572[k] = f_10 * skg_410[k]
                   + f_4 * slf0_275[k]
                   - f_5 * slf1_275[k]
                   + f_3 * pc_x[k] * slg_410[k];

        t_573[k] = f_10 * skg_411[k]
                   + f_6 * slf0_276[k]
                   - f_7 * slf1_276[k]
                   + f_3 * pc_x[k] * slg_411[k];
    }

#pragma omp simd aligned(t_574, t_575, t_576, t_577, pc_x, pc_y, pc_z, skg_303, skg_414, \
                         skg_415, slf0_279, slf1_279, slg_408, slg_410, slg_414, \
                         slg_415 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_574[k] = f_13 * skg_303[k]
                   + f_3 * pc_z[k] * slg_408[k];

        t_575[k] = f_3 * pc_y[k] * slg_410[k];

        t_576[k] = f_10 * skg_414[k]
                   + f_6 * slf0_279[k]
                   - f_7 * slf1_279[k]
                   + f_3 * pc_x[k] * slg_414[k];

        t_577[k] = f_10 * skg_415[k]
                   + f_3 * pc_x[k] * slg_415[k];
    }

#pragma omp simd aligned(t_578, t_579, t_580, t_581, pc_x, skg_416, skg_417, skg_418, skg_419, \
                         slg_416, slg_417, slg_418, slg_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_578[k] = f_10 * skg_416[k]
                   + f_3 * pc_x[k] * slg_416[k];

        t_579[k] = f_10 * skg_417[k]
                   + f_3 * pc_x[k] * slg_417[k];

        t_580[k] = f_10 * skg_418[k]
                   + f_3 * pc_x[k] * slg_418[k];

        t_581[k] = f_10 * skg_419[k]
                   + f_3 * pc_x[k] * slg_419[k];
    }

#pragma omp simd aligned(t_582, t_583, t_584, t_585, pc_y, pc_z, skg_310, slf0_276, slf0_278, \
                         slf0_279, slf1_276, slf1_278, slf1_279, slg_415, slg_417, \
                         slg_418 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_582[k] = f_1 * slf0_276[k]
                   - f_2 * slf1_276[k]
                   + f_3 * pc_y[k] * slg_415[k];

        t_583[k] = f_13 * skg_310[k]
                   + f_3 * pc_z[k] * slg_415[k];

        t_584[k] = f_4 * slf0_278[k]
                   - f_5 * slf1_278[k]
                   + f_3 * pc_y[k] * slg_417[k];

        t_585[k] = f_6 * slf0_279[k]
                   - f_7 * slf1_279[k]
                   + f_3 * pc_y[k] * slg_418[k];
    }
}

static auto
compute_prim_slh_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t skh0,
                                                          const size_t skg, const size_t skh1,
                                                          const size_t slf0, const size_t slf1,
                                                          const size_t slg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_12 = 3.5 / q;
    const auto f_13 = 3.0 / q;
    const auto f_14 = 2.5 / q;
    const auto f_15 = 2.0 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *skh0_441 = buffer.data(skh0 + 441);
    const auto *skh0_444 = buffer.data(skh0 + 444);
    const auto *skh0_447 = buffer.data(skh0 + 447);
    const auto *skh0_588 = buffer.data(skh0 + 588);
    const auto *skh0_591 = buffer.data(skh0 + 591);
    const auto *skh0_593 = buffer.data(skh0 + 593);
    const auto *skh0_594 = buffer.data(skh0 + 594);
    const auto *skh0_597 = buffer.data(skh0 + 597);
    const auto *skh0_603 = buffer.data(skh0 + 603);
    const auto *skh0_605 = buffer.data(skh0 + 605);
    const auto *skh0_606 = buffer.data(skh0 + 606);
    const auto *skh0_608 = buffer.data(skh0 + 608);
    const auto *skh0_614 = buffer.data(skh0 + 614);
    const auto *skh0_618 = buffer.data(skh0 + 618);
    const auto *skh0_624 = buffer.data(skh0 + 624);
    const auto *skh0_626 = buffer.data(skh0 + 626);
    const auto *skh0_627 = buffer.data(skh0 + 627);
    const auto *skh0_629 = buffer.data(skh0 + 629);
    const auto *skh0_630 = buffer.data(skh0 + 630);
    const auto *skh0_633 = buffer.data(skh0 + 633);
    const auto *skh0_635 = buffer.data(skh0 + 635);
    const auto *skh0_636 = buffer.data(skh0 + 636);
    const auto *skh0_639 = buffer.data(skh0 + 639);
    const auto *skh0_645 = buffer.data(skh0 + 645);
    const auto *skh0_647 = buffer.data(skh0 + 647);
    const auto *skh0_648 = buffer.data(skh0 + 648);
    const auto *skh0_650 = buffer.data(skh0 + 650);
    const auto *skh0_651 = buffer.data(skh0 + 651);
    const auto *skh0_654 = buffer.data(skh0 + 654);
    const auto *skh0_656 = buffer.data(skh0 + 656);
    const auto *skh0_657 = buffer.data(skh0 + 657);
    const auto *skh0_660 = buffer.data(skh0 + 660);
    const auto *skh0_666 = buffer.data(skh0 + 666);
    const auto *skh0_668 = buffer.data(skh0 + 668);
    const auto *skh0_669 = buffer.data(skh0 + 669);
    const auto *skh0_671 = buffer.data(skh0 + 671);
    const auto *skh0_672 = buffer.data(skh0 + 672);
    const auto *skh0_675 = buffer.data(skh0 + 675);
    const auto *skh0_677 = buffer.data(skh0 + 677);
    const auto *skh0_678 = buffer.data(skh0 + 678);
    const auto *skh0_681 = buffer.data(skh0 + 681);
    const auto *skh0_687 = buffer.data(skh0 + 687);
    const auto *skh0_689 = buffer.data(skh0 + 689);
    const auto *skh0_690 = buffer.data(skh0 + 690);
    const auto *skh0_692 = buffer.data(skh0 + 692);
    const auto *skh0_693 = buffer.data(skh0 + 693);
    const auto *skh0_696 = buffer.data(skh0 + 696);
    const auto *skh0_698 = buffer.data(skh0 + 698);
    const auto *skh0_699 = buffer.data(skh0 + 699);
    const auto *skh0_702 = buffer.data(skh0 + 702);

    const auto *skg_314 = buffer.data(skg + 314);
    const auto *skg_315 = buffer.data(skg + 315);
    const auto *skg_317 = buffer.data(skg + 317);
    const auto *skg_318 = buffer.data(skg + 318);
    const auto *skg_320 = buffer.data(skg + 320);
    const auto *skg_325 = buffer.data(skg + 325);
    const auto *skg_329 = buffer.data(skg + 329);
    const auto *skg_330 = buffer.data(skg + 330);
    const auto *skg_332 = buffer.data(skg + 332);
    const auto *skg_333 = buffer.data(skg + 333);
    const auto *skg_335 = buffer.data(skg + 335);
    const auto *skg_340 = buffer.data(skg + 340);
    const auto *skg_344 = buffer.data(skg + 344);
    const auto *skg_345 = buffer.data(skg + 345);
    const auto *skg_347 = buffer.data(skg + 347);
    const auto *skg_348 = buffer.data(skg + 348);
    const auto *skg_350 = buffer.data(skg + 350);
    const auto *skg_355 = buffer.data(skg + 355);
    const auto *skg_359 = buffer.data(skg + 359);
    const auto *skg_360 = buffer.data(skg + 360);
    const auto *skg_362 = buffer.data(skg + 362);
    const auto *skg_363 = buffer.data(skg + 363);
    const auto *skg_365 = buffer.data(skg + 365);
    const auto *skg_370 = buffer.data(skg + 370);
    const auto *skg_374 = buffer.data(skg + 374);
    const auto *skg_375 = buffer.data(skg + 375);
    const auto *skg_377 = buffer.data(skg + 377);
    const auto *skg_378 = buffer.data(skg + 378);
    const auto *skg_380 = buffer.data(skg + 380);
    const auto *skg_389 = buffer.data(skg + 389);
    const auto *skg_390 = buffer.data(skg + 390);
    const auto *skg_392 = buffer.data(skg + 392);
    const auto *skg_395 = buffer.data(skg + 395);
    const auto *skg_420 = buffer.data(skg + 420);
    const auto *skg_423 = buffer.data(skg + 423);
    const auto *skg_425 = buffer.data(skg + 425);
    const auto *skg_426 = buffer.data(skg + 426);
    const auto *skg_429 = buffer.data(skg + 429);
    const auto *skg_430 = buffer.data(skg + 430);
    const auto *skg_431 = buffer.data(skg + 431);
    const auto *skg_432 = buffer.data(skg + 432);
    const auto *skg_433 = buffer.data(skg + 433);
    const auto *skg_434 = buffer.data(skg + 434);
    const auto *skg_440 = buffer.data(skg + 440);
    const auto *skg_444 = buffer.data(skg + 444);
    const auto *skg_445 = buffer.data(skg + 445);
    const auto *skg_446 = buffer.data(skg + 446);
    const auto *skg_447 = buffer.data(skg + 447);
    const auto *skg_448 = buffer.data(skg + 448);
    const auto *skg_449 = buffer.data(skg + 449);
    const auto *skg_450 = buffer.data(skg + 450);
    const auto *skg_453 = buffer.data(skg + 453);
    const auto *skg_455 = buffer.data(skg + 455);
    const auto *skg_456 = buffer.data(skg + 456);
    const auto *skg_459 = buffer.data(skg + 459);
    const auto *skg_460 = buffer.data(skg + 460);
    const auto *skg_461 = buffer.data(skg + 461);
    const auto *skg_462 = buffer.data(skg + 462);
    const auto *skg_463 = buffer.data(skg + 463);
    const auto *skg_464 = buffer.data(skg + 464);
    const auto *skg_465 = buffer.data(skg + 465);
    const auto *skg_468 = buffer.data(skg + 468);
    const auto *skg_470 = buffer.data(skg + 470);
    const auto *skg_471 = buffer.data(skg + 471);
    const auto *skg_474 = buffer.data(skg + 474);
    const auto *skg_475 = buffer.data(skg + 475);
    const auto *skg_476 = buffer.data(skg + 476);
    const auto *skg_477 = buffer.data(skg + 477);
    const auto *skg_478 = buffer.data(skg + 478);
    const auto *skg_479 = buffer.data(skg + 479);
    const auto *skg_480 = buffer.data(skg + 480);
    const auto *skg_483 = buffer.data(skg + 483);
    const auto *skg_485 = buffer.data(skg + 485);
    const auto *skg_486 = buffer.data(skg + 486);
    const auto *skg_489 = buffer.data(skg + 489);
    const auto *skg_490 = buffer.data(skg + 490);
    const auto *skg_491 = buffer.data(skg + 491);
    const auto *skg_492 = buffer.data(skg + 492);
    const auto *skg_493 = buffer.data(skg + 493);
    const auto *skg_494 = buffer.data(skg + 494);
    const auto *skg_495 = buffer.data(skg + 495);
    const auto *skg_498 = buffer.data(skg + 498);
    const auto *skg_500 = buffer.data(skg + 500);
    const auto *skg_501 = buffer.data(skg + 501);
    const auto *skg_504 = buffer.data(skg + 504);
    const auto *skg_505 = buffer.data(skg + 505);
    const auto *skg_506 = buffer.data(skg + 506);
    const auto *skg_507 = buffer.data(skg + 507);

    const auto *skh1_441 = buffer.data(skh1 + 441);
    const auto *skh1_444 = buffer.data(skh1 + 444);
    const auto *skh1_447 = buffer.data(skh1 + 447);
    const auto *skh1_588 = buffer.data(skh1 + 588);
    const auto *skh1_591 = buffer.data(skh1 + 591);
    const auto *skh1_593 = buffer.data(skh1 + 593);
    const auto *skh1_594 = buffer.data(skh1 + 594);
    const auto *skh1_597 = buffer.data(skh1 + 597);
    const auto *skh1_603 = buffer.data(skh1 + 603);
    const auto *skh1_605 = buffer.data(skh1 + 605);
    const auto *skh1_606 = buffer.data(skh1 + 606);
    const auto *skh1_608 = buffer.data(skh1 + 608);
    const auto *skh1_614 = buffer.data(skh1 + 614);
    const auto *skh1_618 = buffer.data(skh1 + 618);
    const auto *skh1_624 = buffer.data(skh1 + 624);
    const auto *skh1_626 = buffer.data(skh1 + 626);
    const auto *skh1_627 = buffer.data(skh1 + 627);
    const auto *skh1_629 = buffer.data(skh1 + 629);
    const auto *skh1_630 = buffer.data(skh1 + 630);
    const auto *skh1_633 = buffer.data(skh1 + 633);
    const auto *skh1_635 = buffer.data(skh1 + 635);
    const auto *skh1_636 = buffer.data(skh1 + 636);
    const auto *skh1_639 = buffer.data(skh1 + 639);
    const auto *skh1_645 = buffer.data(skh1 + 645);
    const auto *skh1_647 = buffer.data(skh1 + 647);
    const auto *skh1_648 = buffer.data(skh1 + 648);
    const auto *skh1_650 = buffer.data(skh1 + 650);
    const auto *skh1_651 = buffer.data(skh1 + 651);
    const auto *skh1_654 = buffer.data(skh1 + 654);
    const auto *skh1_656 = buffer.data(skh1 + 656);
    const auto *skh1_657 = buffer.data(skh1 + 657);
    const auto *skh1_660 = buffer.data(skh1 + 660);
    const auto *skh1_666 = buffer.data(skh1 + 666);
    const auto *skh1_668 = buffer.data(skh1 + 668);
    const auto *skh1_669 = buffer.data(skh1 + 669);
    const auto *skh1_671 = buffer.data(skh1 + 671);
    const auto *skh1_672 = buffer.data(skh1 + 672);
    const auto *skh1_675 = buffer.data(skh1 + 675);
    const auto *skh1_677 = buffer.data(skh1 + 677);
    const auto *skh1_678 = buffer.data(skh1 + 678);
    const auto *skh1_681 = buffer.data(skh1 + 681);
    const auto *skh1_687 = buffer.data(skh1 + 687);
    const auto *skh1_689 = buffer.data(skh1 + 689);
    const auto *skh1_690 = buffer.data(skh1 + 690);
    const auto *skh1_692 = buffer.data(skh1 + 692);
    const auto *skh1_693 = buffer.data(skh1 + 693);
    const auto *skh1_696 = buffer.data(skh1 + 696);
    const auto *skh1_698 = buffer.data(skh1 + 698);
    const auto *skh1_699 = buffer.data(skh1 + 699);
    const auto *skh1_702 = buffer.data(skh1 + 702);

    const auto *slf0_279 = buffer.data(slf0 + 279);

    const auto *slf1_279 = buffer.data(slf1 + 279);

    const auto *slg_419 = buffer.data(slg + 419);
    const auto *slg_420 = buffer.data(slg + 420);
    const auto *slg_422 = buffer.data(slg + 422);
    const auto *slg_423 = buffer.data(slg + 423);
    const auto *slg_425 = buffer.data(slg + 425);
    const auto *slg_430 = buffer.data(slg + 430);
    const auto *slg_431 = buffer.data(slg + 431);
    const auto *slg_432 = buffer.data(slg + 432);
    const auto *slg_433 = buffer.data(slg + 433);
    const auto *slg_434 = buffer.data(slg + 434);
    const auto *slg_435 = buffer.data(slg + 435);
    const auto *slg_437 = buffer.data(slg + 437);
    const auto *slg_438 = buffer.data(slg + 438);
    const auto *slg_440 = buffer.data(slg + 440);
    const auto *slg_445 = buffer.data(slg + 445);
    const auto *slg_446 = buffer.data(slg + 446);
    const auto *slg_447 = buffer.data(slg + 447);
    const auto *slg_448 = buffer.data(slg + 448);
    const auto *slg_449 = buffer.data(slg + 449);
    const auto *slg_450 = buffer.data(slg + 450);
    const auto *slg_452 = buffer.data(slg + 452);
    const auto *slg_453 = buffer.data(slg + 453);
    const auto *slg_455 = buffer.data(slg + 455);
    const auto *slg_460 = buffer.data(slg + 460);
    const auto *slg_461 = buffer.data(slg + 461);
    const auto *slg_462 = buffer.data(slg + 462);
    const auto *slg_463 = buffer.data(slg + 463);
    const auto *slg_464 = buffer.data(slg + 464);
    const auto *slg_465 = buffer.data(slg + 465);
    const auto *slg_467 = buffer.data(slg + 467);
    const auto *slg_468 = buffer.data(slg + 468);
    const auto *slg_470 = buffer.data(slg + 470);
    const auto *slg_475 = buffer.data(slg + 475);
    const auto *slg_476 = buffer.data(slg + 476);
    const auto *slg_477 = buffer.data(slg + 477);
    const auto *slg_478 = buffer.data(slg + 478);
    const auto *slg_479 = buffer.data(slg + 479);
    const auto *slg_480 = buffer.data(slg + 480);
    const auto *slg_482 = buffer.data(slg + 482);
    const auto *slg_483 = buffer.data(slg + 483);
    const auto *slg_485 = buffer.data(slg + 485);
    const auto *slg_490 = buffer.data(slg + 490);
    const auto *slg_491 = buffer.data(slg + 491);
    const auto *slg_492 = buffer.data(slg + 492);
    const auto *slg_493 = buffer.data(slg + 493);
    const auto *slg_494 = buffer.data(slg + 494);
    const auto *slg_495 = buffer.data(slg + 495);
    const auto *slg_497 = buffer.data(slg + 497);
    const auto *slg_498 = buffer.data(slg + 498);
    const auto *slg_500 = buffer.data(slg + 500);
    const auto *slg_505 = buffer.data(slg + 505);
    const auto *slg_506 = buffer.data(slg + 506);
    const auto *slg_507 = buffer.data(slg + 507);

#pragma omp simd aligned(t_586, t_587, t_588, pb_x, pc_x, pc_y, pc_z, skh0_588, skg_314, \
                         skg_420, skh1_588, slf0_279, slf1_279, \
                         slg_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_586[k] = f_3 * pc_y[k] * slg_419[k];

        t_587[k] = f_13 * skg_314[k]
                   + f_1 * slf0_279[k]
                   - f_2 * slf1_279[k]
                   + f_3 * pc_z[k] * slg_419[k];

        t_588[k] = pb_x[k] * skh0_588[k]
                   + f_14 * skg_420[k]
                   - f_8 * pc_x[k] * skh1_588[k];
    }

#pragma omp simd aligned(t_589, t_590, t_591, t_592, pb_x, pc_x, pc_y, pc_z, skh0_591, \
                         skg_315, skg_317, skg_423, skh1_591, slg_420, \
                         slg_422 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_589[k] = f_12 * skg_315[k]
                   + f_3 * pc_y[k] * slg_420[k];

        t_590[k] = f_3 * pc_z[k] * slg_420[k];

        t_591[k] = pb_x[k] * skh0_591[k]
                   + f_11 * skg_423[k]
                   - f_8 * pc_x[k] * skh1_591[k];

        t_592[k] = f_12 * skg_317[k]
                   + f_3 * pc_y[k] * slg_422[k];
    }

#pragma omp simd aligned(t_593, t_594, t_595, pb_x, pc_x, pc_z, skh0_593, skh0_594, skg_425, \
                         skg_426, skh1_593, skh1_594, slg_423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_593[k] = pb_x[k] * skh0_593[k]
                   + f_11 * skg_425[k]
                   - f_8 * pc_x[k] * skh1_593[k];

        t_594[k] = pb_x[k] * skh0_594[k]
                   + f_10 * skg_426[k]
                   - f_8 * pc_x[k] * skh1_594[k];

        t_595[k] = f_3 * pc_z[k] * slg_423[k];
    }

#pragma omp simd aligned(t_596, t_597, t_598, t_599, pb_x, pc_x, pc_y, skh0_597, skg_320, \
                         skg_429, skg_430, skg_431, skh1_597, slg_425, slg_430, \
                         slg_431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_596[k] = f_12 * skg_320[k]
                   + f_3 * pc_y[k] * slg_425[k];

        t_597[k] = pb_x[k] * skh0_597[k]
                   + f_10 * skg_429[k]
                   - f_8 * pc_x[k] * skh1_597[k];

        t_598[k] = f_9 * skg_430[k]
                   + f_3 * pc_x[k] * slg_430[k];

        t_599[k] = f_9 * skg_431[k]
                   + f_3 * pc_x[k] * slg_431[k];
    }

#pragma omp simd aligned(t_600, t_601, t_602, t_603, pb_x, pc_x, skh0_603, skg_432, skg_433, \
                         skg_434, skh1_603, slg_432, slg_433, slg_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = f_9 * skg_432[k]
                   + f_3 * pc_x[k] * slg_432[k];

        t_601[k] = f_9 * skg_433[k]
                   + f_3 * pc_x[k] * slg_433[k];

        t_602[k] = f_9 * skg_434[k]
                   + f_3 * pc_x[k] * slg_434[k];

        t_603[k] = pb_x[k] * skh0_603[k]
                   - f_8 * pc_x[k] * skh1_603[k];
    }

#pragma omp simd aligned(t_604, t_605, t_606, t_607, pb_x, pc_x, pc_y, pc_z, skh0_605, \
                         skh0_606, skg_329, skh1_605, skh1_606, slg_430, \
                         slg_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_604[k] = f_3 * pc_z[k] * slg_430[k];

        t_605[k] = pb_x[k] * skh0_605[k]
                   - f_8 * pc_x[k] * skh1_605[k];

        t_606[k] = pb_x[k] * skh0_606[k]
                   - f_8 * pc_x[k] * skh1_606[k];

        t_607[k] = f_12 * skg_329[k]
                   + f_3 * pc_y[k] * slg_434[k];
    }

#pragma omp simd aligned(t_608, t_609, t_610, t_611, pb_x, pb_z, pc_x, pc_y, pc_z, skh0_441, \
                         skh0_608, skg_315, skg_330, skh1_441, skh1_608, \
                         slg_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_608[k] = pb_x[k] * skh0_608[k]
                   - f_8 * pc_x[k] * skh1_608[k];

        t_609[k] = pb_z[k] * skh0_441[k]
                   - f_8 * pc_z[k] * skh1_441[k];

        t_610[k] = f_13 * skg_330[k]
                   + f_3 * pc_y[k] * slg_435[k];

        t_611[k] = f_9 * skg_315[k]
                   + f_3 * pc_z[k] * slg_435[k];
    }

#pragma omp simd aligned(t_612, t_613, t_614, pb_x, pb_z, pc_x, pc_y, pc_z, skh0_444, \
                         skh0_614, skg_332, skg_440, skh1_444, skh1_614, \
                         slg_437 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_612[k] = pb_z[k] * skh0_444[k]
                   - f_8 * pc_z[k] * skh1_444[k];

        t_613[k] = f_13 * skg_332[k]
                   + f_3 * pc_y[k] * slg_437[k];

        t_614[k] = pb_x[k] * skh0_614[k]
                   + f_11 * skg_440[k]
                   - f_8 * pc_x[k] * skh1_614[k];
    }

#pragma omp simd aligned(t_615, t_616, t_617, pb_z, pc_y, pc_z, skh0_447, skg_318, skg_335, \
                         skh1_447, slg_438, slg_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_615[k] = pb_z[k] * skh0_447[k]
                   - f_8 * pc_z[k] * skh1_447[k];

        t_616[k] = f_9 * skg_318[k]
                   + f_3 * pc_z[k] * slg_438[k];

        t_617[k] = f_13 * skg_335[k]
                   + f_3 * pc_y[k] * slg_440[k];
    }

#pragma omp simd aligned(t_618, t_619, t_620, t_621, pb_x, pc_x, skh0_618, skg_444, skg_445, \
                         skg_446, skg_447, skh1_618, slg_445, slg_446, \
                         slg_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_618[k] = pb_x[k] * skh0_618[k]
                   + f_10 * skg_444[k]
                   - f_8 * pc_x[k] * skh1_618[k];

        t_619[k] = f_9 * skg_445[k]
                   + f_3 * pc_x[k] * slg_445[k];

        t_620[k] = f_9 * skg_446[k]
                   + f_3 * pc_x[k] * slg_446[k];

        t_621[k] = f_9 * skg_447[k]
                   + f_3 * pc_x[k] * slg_447[k];
    }

#pragma omp simd aligned(t_622, t_623, t_624, t_625, pb_x, pc_x, pc_z, skh0_624, skg_325, \
                         skg_448, skg_449, skh1_624, slg_445, slg_448, \
                         slg_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_622[k] = f_9 * skg_448[k]
                   + f_3 * pc_x[k] * slg_448[k];

        t_623[k] = f_9 * skg_449[k]
                   + f_3 * pc_x[k] * slg_449[k];

        t_624[k] = pb_x[k] * skh0_624[k]
                   - f_8 * pc_x[k] * skh1_624[k];

        t_625[k] = f_9 * skg_325[k]
                   + f_3 * pc_z[k] * slg_445[k];
    }

#pragma omp simd aligned(t_626, t_627, t_628, t_629, pb_x, pc_x, pc_y, skh0_626, skh0_627, \
                         skh0_629, skg_344, skh1_626, skh1_627, skh1_629, \
                         slg_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_626[k] = pb_x[k] * skh0_626[k]
                   - f_8 * pc_x[k] * skh1_626[k];

        t_627[k] = pb_x[k] * skh0_627[k]
                   - f_8 * pc_x[k] * skh1_627[k];

        t_628[k] = f_13 * skg_344[k]
                   + f_3 * pc_y[k] * slg_449[k];

        t_629[k] = pb_x[k] * skh0_629[k]
                   - f_8 * pc_x[k] * skh1_629[k];
    }

#pragma omp simd aligned(t_630, t_631, t_632, pb_x, pc_x, pc_y, pc_z, skh0_630, skg_330, \
                         skg_345, skg_450, skh1_630, slg_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_630[k] = pb_x[k] * skh0_630[k]
                   + f_14 * skg_450[k]
                   - f_8 * pc_x[k] * skh1_630[k];

        t_631[k] = f_14 * skg_345[k]
                   + f_3 * pc_y[k] * slg_450[k];

        t_632[k] = f_10 * skg_330[k]
                   + f_3 * pc_z[k] * slg_450[k];
    }

#pragma omp simd aligned(t_633, t_634, t_635, pb_x, pc_x, pc_y, skh0_633, skh0_635, skg_347, \
                         skg_453, skg_455, skh1_633, skh1_635, \
                         slg_452 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_633[k] = pb_x[k] * skh0_633[k]
                   + f_11 * skg_453[k]
                   - f_8 * pc_x[k] * skh1_633[k];

        t_634[k] = f_14 * skg_347[k]
                   + f_3 * pc_y[k] * slg_452[k];

        t_635[k] = pb_x[k] * skh0_635[k]
                   + f_11 * skg_455[k]
                   - f_8 * pc_x[k] * skh1_635[k];
    }

#pragma omp simd aligned(t_636, t_637, t_638, pb_x, pc_x, pc_y, pc_z, skh0_636, skg_333, \
                         skg_350, skg_456, skh1_636, slg_453, slg_455 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_636[k] = pb_x[k] * skh0_636[k]
                   + f_10 * skg_456[k]
                   - f_8 * pc_x[k] * skh1_636[k];

        t_637[k] = f_10 * skg_333[k]
                   + f_3 * pc_z[k] * slg_453[k];

        t_638[k] = f_14 * skg_350[k]
                   + f_3 * pc_y[k] * slg_455[k];
    }

#pragma omp simd aligned(t_639, t_640, t_641, t_642, pb_x, pc_x, skh0_639, skg_459, skg_460, \
                         skg_461, skg_462, skh1_639, slg_460, slg_461, \
                         slg_462 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_639[k] = pb_x[k] * skh0_639[k]
                   + f_10 * skg_459[k]
                   - f_8 * pc_x[k] * skh1_639[k];

        t_640[k] = f_9 * skg_460[k]
                   + f_3 * pc_x[k] * slg_460[k];

        t_641[k] = f_9 * skg_461[k]
                   + f_3 * pc_x[k] * slg_461[k];

        t_642[k] = f_9 * skg_462[k]
                   + f_3 * pc_x[k] * slg_462[k];
    }

#pragma omp simd aligned(t_643, t_644, t_645, t_646, pb_x, pc_x, pc_z, skh0_645, skg_340, \
                         skg_463, skg_464, skh1_645, slg_460, slg_463, \
                         slg_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_643[k] = f_9 * skg_463[k]
                   + f_3 * pc_x[k] * slg_463[k];

        t_644[k] = f_9 * skg_464[k]
                   + f_3 * pc_x[k] * slg_464[k];

        t_645[k] = pb_x[k] * skh0_645[k]
                   - f_8 * pc_x[k] * skh1_645[k];

        t_646[k] = f_10 * skg_340[k]
                   + f_3 * pc_z[k] * slg_460[k];
    }

#pragma omp simd aligned(t_647, t_648, t_649, t_650, pb_x, pc_x, pc_y, skh0_647, skh0_648, \
                         skh0_650, skg_359, skh1_647, skh1_648, skh1_650, \
                         slg_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_647[k] = pb_x[k] * skh0_647[k]
                   - f_8 * pc_x[k] * skh1_647[k];

        t_648[k] = pb_x[k] * skh0_648[k]
                   - f_8 * pc_x[k] * skh1_648[k];

        t_649[k] = f_14 * skg_359[k]
                   + f_3 * pc_y[k] * slg_464[k];

        t_650[k] = pb_x[k] * skh0_650[k]
                   - f_8 * pc_x[k] * skh1_650[k];
    }

#pragma omp simd aligned(t_651, t_652, t_653, pb_x, pc_x, pc_y, pc_z, skh0_651, skg_345, \
                         skg_360, skg_465, skh1_651, slg_465 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_651[k] = pb_x[k] * skh0_651[k]
                   + f_14 * skg_465[k]
                   - f_8 * pc_x[k] * skh1_651[k];

        t_652[k] = f_15 * skg_360[k]
                   + f_3 * pc_y[k] * slg_465[k];

        t_653[k] = f_11 * skg_345[k]
                   + f_3 * pc_z[k] * slg_465[k];
    }

#pragma omp simd aligned(t_654, t_655, t_656, pb_x, pc_x, pc_y, skh0_654, skh0_656, skg_362, \
                         skg_468, skg_470, skh1_654, skh1_656, \
                         slg_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_654[k] = pb_x[k] * skh0_654[k]
                   + f_11 * skg_468[k]
                   - f_8 * pc_x[k] * skh1_654[k];

        t_655[k] = f_15 * skg_362[k]
                   + f_3 * pc_y[k] * slg_467[k];

        t_656[k] = pb_x[k] * skh0_656[k]
                   + f_11 * skg_470[k]
                   - f_8 * pc_x[k] * skh1_656[k];
    }

#pragma omp simd aligned(t_657, t_658, t_659, pb_x, pc_x, pc_y, pc_z, skh0_657, skg_348, \
                         skg_365, skg_471, skh1_657, slg_468, slg_470 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_657[k] = pb_x[k] * skh0_657[k]
                   + f_10 * skg_471[k]
                   - f_8 * pc_x[k] * skh1_657[k];

        t_658[k] = f_11 * skg_348[k]
                   + f_3 * pc_z[k] * slg_468[k];

        t_659[k] = f_15 * skg_365[k]
                   + f_3 * pc_y[k] * slg_470[k];
    }

#pragma omp simd aligned(t_660, t_661, t_662, t_663, pb_x, pc_x, skh0_660, skg_474, skg_475, \
                         skg_476, skg_477, skh1_660, slg_475, slg_476, \
                         slg_477 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_660[k] = pb_x[k] * skh0_660[k]
                   + f_10 * skg_474[k]
                   - f_8 * pc_x[k] * skh1_660[k];

        t_661[k] = f_9 * skg_475[k]
                   + f_3 * pc_x[k] * slg_475[k];

        t_662[k] = f_9 * skg_476[k]
                   + f_3 * pc_x[k] * slg_476[k];

        t_663[k] = f_9 * skg_477[k]
                   + f_3 * pc_x[k] * slg_477[k];
    }

#pragma omp simd aligned(t_664, t_665, t_666, t_667, pb_x, pc_x, pc_z, skh0_666, skg_355, \
                         skg_478, skg_479, skh1_666, slg_475, slg_478, \
                         slg_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_664[k] = f_9 * skg_478[k]
                   + f_3 * pc_x[k] * slg_478[k];

        t_665[k] = f_9 * skg_479[k]
                   + f_3 * pc_x[k] * slg_479[k];

        t_666[k] = pb_x[k] * skh0_666[k]
                   - f_8 * pc_x[k] * skh1_666[k];

        t_667[k] = f_11 * skg_355[k]
                   + f_3 * pc_z[k] * slg_475[k];
    }

#pragma omp simd aligned(t_668, t_669, t_670, t_671, pb_x, pc_x, pc_y, skh0_668, skh0_669, \
                         skh0_671, skg_374, skh1_668, skh1_669, skh1_671, \
                         slg_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_668[k] = pb_x[k] * skh0_668[k]
                   - f_8 * pc_x[k] * skh1_668[k];

        t_669[k] = pb_x[k] * skh0_669[k]
                   - f_8 * pc_x[k] * skh1_669[k];

        t_670[k] = f_15 * skg_374[k]
                   + f_3 * pc_y[k] * slg_479[k];

        t_671[k] = pb_x[k] * skh0_671[k]
                   - f_8 * pc_x[k] * skh1_671[k];
    }

#pragma omp simd aligned(t_672, t_673, t_674, pb_x, pc_x, pc_y, pc_z, skh0_672, skg_360, \
                         skg_375, skg_480, skh1_672, slg_480 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_672[k] = pb_x[k] * skh0_672[k]
                   + f_14 * skg_480[k]
                   - f_8 * pc_x[k] * skh1_672[k];

        t_673[k] = f_11 * skg_375[k]
                   + f_3 * pc_y[k] * slg_480[k];

        t_674[k] = f_15 * skg_360[k]
                   + f_3 * pc_z[k] * slg_480[k];
    }

#pragma omp simd aligned(t_675, t_676, t_677, pb_x, pc_x, pc_y, skh0_675, skh0_677, skg_377, \
                         skg_483, skg_485, skh1_675, skh1_677, \
                         slg_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_675[k] = pb_x[k] * skh0_675[k]
                   + f_11 * skg_483[k]
                   - f_8 * pc_x[k] * skh1_675[k];

        t_676[k] = f_11 * skg_377[k]
                   + f_3 * pc_y[k] * slg_482[k];

        t_677[k] = pb_x[k] * skh0_677[k]
                   + f_11 * skg_485[k]
                   - f_8 * pc_x[k] * skh1_677[k];
    }

#pragma omp simd aligned(t_678, t_679, t_680, pb_x, pc_x, pc_y, pc_z, skh0_678, skg_363, \
                         skg_380, skg_486, skh1_678, slg_483, slg_485 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_678[k] = pb_x[k] * skh0_678[k]
                   + f_10 * skg_486[k]
                   - f_8 * pc_x[k] * skh1_678[k];

        t_679[k] = f_15 * skg_363[k]
                   + f_3 * pc_z[k] * slg_483[k];

        t_680[k] = f_11 * skg_380[k]
                   + f_3 * pc_y[k] * slg_485[k];
    }

#pragma omp simd aligned(t_681, t_682, t_683, t_684, pb_x, pc_x, skh0_681, skg_489, skg_490, \
                         skg_491, skg_492, skh1_681, slg_490, slg_491, \
                         slg_492 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_681[k] = pb_x[k] * skh0_681[k]
                   + f_10 * skg_489[k]
                   - f_8 * pc_x[k] * skh1_681[k];

        t_682[k] = f_9 * skg_490[k]
                   + f_3 * pc_x[k] * slg_490[k];

        t_683[k] = f_9 * skg_491[k]
                   + f_3 * pc_x[k] * slg_491[k];

        t_684[k] = f_9 * skg_492[k]
                   + f_3 * pc_x[k] * slg_492[k];
    }

#pragma omp simd aligned(t_685, t_686, t_687, t_688, pb_x, pc_x, pc_z, skh0_687, skg_370, \
                         skg_493, skg_494, skh1_687, slg_490, slg_493, \
                         slg_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_685[k] = f_9 * skg_493[k]
                   + f_3 * pc_x[k] * slg_493[k];

        t_686[k] = f_9 * skg_494[k]
                   + f_3 * pc_x[k] * slg_494[k];

        t_687[k] = pb_x[k] * skh0_687[k]
                   - f_8 * pc_x[k] * skh1_687[k];

        t_688[k] = f_15 * skg_370[k]
                   + f_3 * pc_z[k] * slg_490[k];
    }

#pragma omp simd aligned(t_689, t_690, t_691, t_692, pb_x, pc_x, pc_y, skh0_689, skh0_690, \
                         skh0_692, skg_389, skh1_689, skh1_690, skh1_692, \
                         slg_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_689[k] = pb_x[k] * skh0_689[k]
                   - f_8 * pc_x[k] * skh1_689[k];

        t_690[k] = pb_x[k] * skh0_690[k]
                   - f_8 * pc_x[k] * skh1_690[k];

        t_691[k] = f_11 * skg_389[k]
                   + f_3 * pc_y[k] * slg_494[k];

        t_692[k] = pb_x[k] * skh0_692[k]
                   - f_8 * pc_x[k] * skh1_692[k];
    }

#pragma omp simd aligned(t_693, t_694, t_695, pb_x, pc_x, pc_y, pc_z, skh0_693, skg_375, \
                         skg_390, skg_495, skh1_693, slg_495 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_693[k] = pb_x[k] * skh0_693[k]
                   + f_14 * skg_495[k]
                   - f_8 * pc_x[k] * skh1_693[k];

        t_694[k] = f_10 * skg_390[k]
                   + f_3 * pc_y[k] * slg_495[k];

        t_695[k] = f_14 * skg_375[k]
                   + f_3 * pc_z[k] * slg_495[k];
    }

#pragma omp simd aligned(t_696, t_697, t_698, pb_x, pc_x, pc_y, skh0_696, skh0_698, skg_392, \
                         skg_498, skg_500, skh1_696, skh1_698, \
                         slg_497 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_696[k] = pb_x[k] * skh0_696[k]
                   + f_11 * skg_498[k]
                   - f_8 * pc_x[k] * skh1_696[k];

        t_697[k] = f_10 * skg_392[k]
                   + f_3 * pc_y[k] * slg_497[k];

        t_698[k] = pb_x[k] * skh0_698[k]
                   + f_11 * skg_500[k]
                   - f_8 * pc_x[k] * skh1_698[k];
    }

#pragma omp simd aligned(t_699, t_700, t_701, pb_x, pc_x, pc_y, pc_z, skh0_699, skg_378, \
                         skg_395, skg_501, skh1_699, slg_498, slg_500 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_699[k] = pb_x[k] * skh0_699[k]
                   + f_10 * skg_501[k]
                   - f_8 * pc_x[k] * skh1_699[k];

        t_700[k] = f_14 * skg_378[k]
                   + f_3 * pc_z[k] * slg_498[k];

        t_701[k] = f_10 * skg_395[k]
                   + f_3 * pc_y[k] * slg_500[k];
    }

#pragma omp simd aligned(t_702, t_703, t_704, t_705, pb_x, pc_x, skh0_702, skg_504, skg_505, \
                         skg_506, skg_507, skh1_702, slg_505, slg_506, \
                         slg_507 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_702[k] = pb_x[k] * skh0_702[k]
                   + f_10 * skg_504[k]
                   - f_8 * pc_x[k] * skh1_702[k];

        t_703[k] = f_9 * skg_505[k]
                   + f_3 * pc_x[k] * slg_505[k];

        t_704[k] = f_9 * skg_506[k]
                   + f_3 * pc_x[k] * slg_506[k];

        t_705[k] = f_9 * skg_507[k]
                   + f_3 * pc_x[k] * slg_507[k];
    }
}

static auto
compute_prim_slh_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t skh0,
                                                          const size_t skg, const size_t skh1,
                                                          const size_t slf0, const size_t slf1,
                                                          const size_t slg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
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
    const auto f_12 = 3.5 / q;
    const auto f_13 = 3.0 / q;
    const auto f_14 = 2.5 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *skh0_567 = buffer.data(skh0 + 567);
    const auto *skh0_572 = buffer.data(skh0 + 572);
    const auto *skh0_576 = buffer.data(skh0 + 576);
    const auto *skh0_588 = buffer.data(skh0 + 588);
    const auto *skh0_591 = buffer.data(skh0 + 591);
    const auto *skh0_594 = buffer.data(skh0 + 594);
    const auto *skh0_603 = buffer.data(skh0 + 603);
    const auto *skh0_605 = buffer.data(skh0 + 605);
    const auto *skh0_606 = buffer.data(skh0 + 606);
    const auto *skh0_708 = buffer.data(skh0 + 708);
    const auto *skh0_710 = buffer.data(skh0 + 710);
    const auto *skh0_711 = buffer.data(skh0 + 711);
    const auto *skh0_713 = buffer.data(skh0 + 713);
    const auto *skh0_717 = buffer.data(skh0 + 717);
    const auto *skh0_720 = buffer.data(skh0 + 720);
    const auto *skh0_729 = buffer.data(skh0 + 729);
    const auto *skh0_731 = buffer.data(skh0 + 731);
    const auto *skh0_732 = buffer.data(skh0 + 732);
    const auto *skh0_734 = buffer.data(skh0 + 734);
    const auto *skh0_735 = buffer.data(skh0 + 735);
    const auto *skh0_738 = buffer.data(skh0 + 738);
    const auto *skh0_740 = buffer.data(skh0 + 740);
    const auto *skh0_741 = buffer.data(skh0 + 741);
    const auto *skh0_744 = buffer.data(skh0 + 744);
    const auto *skh0_750 = buffer.data(skh0 + 750);
    const auto *skh0_752 = buffer.data(skh0 + 752);
    const auto *skh0_753 = buffer.data(skh0 + 753);
    const auto *skh0_755 = buffer.data(skh0 + 755);

    const auto *skg_385 = buffer.data(skg + 385);
    const auto *skg_390 = buffer.data(skg + 390);
    const auto *skg_393 = buffer.data(skg + 393);
    const auto *skg_400 = buffer.data(skg + 400);
    const auto *skg_404 = buffer.data(skg + 404);
    const auto *skg_405 = buffer.data(skg + 405);
    const auto *skg_407 = buffer.data(skg + 407);
    const auto *skg_408 = buffer.data(skg + 408);
    const auto *skg_410 = buffer.data(skg + 410);
    const auto *skg_415 = buffer.data(skg + 415);
    const auto *skg_419 = buffer.data(skg + 419);
    const auto *skg_420 = buffer.data(skg + 420);
    const auto *skg_422 = buffer.data(skg + 422);
    const auto *skg_423 = buffer.data(skg + 423);
    const auto *skg_425 = buffer.data(skg + 425);
    const auto *skg_430 = buffer.data(skg + 430);
    const auto *skg_431 = buffer.data(skg + 431);
    const auto *skg_432 = buffer.data(skg + 432);
    const auto *skg_433 = buffer.data(skg + 433);
    const auto *skg_434 = buffer.data(skg + 434);
    const auto *skg_435 = buffer.data(skg + 435);
    const auto *skg_437 = buffer.data(skg + 437);
    const auto *skg_438 = buffer.data(skg + 438);
    const auto *skg_440 = buffer.data(skg + 440);
    const auto *skg_445 = buffer.data(skg + 445);
    const auto *skg_449 = buffer.data(skg + 449);
    const auto *skg_450 = buffer.data(skg + 450);
    const auto *skg_452 = buffer.data(skg + 452);
    const auto *skg_453 = buffer.data(skg + 453);
    const auto *skg_455 = buffer.data(skg + 455);
    const auto *skg_460 = buffer.data(skg + 460);
    const auto *skg_462 = buffer.data(skg + 462);
    const auto *skg_463 = buffer.data(skg + 463);
    const auto *skg_464 = buffer.data(skg + 464);
    const auto *skg_465 = buffer.data(skg + 465);
    const auto *skg_467 = buffer.data(skg + 467);
    const auto *skg_470 = buffer.data(skg + 470);
    const auto *skg_475 = buffer.data(skg + 475);
    const auto *skg_508 = buffer.data(skg + 508);
    const auto *skg_509 = buffer.data(skg + 509);
    const auto *skg_513 = buffer.data(skg + 513);
    const auto *skg_516 = buffer.data(skg + 516);
    const auto *skg_520 = buffer.data(skg + 520);
    const auto *skg_521 = buffer.data(skg + 521);
    const auto *skg_522 = buffer.data(skg + 522);
    const auto *skg_523 = buffer.data(skg + 523);
    const auto *skg_524 = buffer.data(skg + 524);
    const auto *skg_525 = buffer.data(skg + 525);
    const auto *skg_528 = buffer.data(skg + 528);
    const auto *skg_530 = buffer.data(skg + 530);
    const auto *skg_531 = buffer.data(skg + 531);
    const auto *skg_534 = buffer.data(skg + 534);
    const auto *skg_535 = buffer.data(skg + 535);
    const auto *skg_536 = buffer.data(skg + 536);
    const auto *skg_537 = buffer.data(skg + 537);
    const auto *skg_538 = buffer.data(skg + 538);
    const auto *skg_539 = buffer.data(skg + 539);

    const auto *skh1_567 = buffer.data(skh1 + 567);
    const auto *skh1_572 = buffer.data(skh1 + 572);
    const auto *skh1_576 = buffer.data(skh1 + 576);
    const auto *skh1_588 = buffer.data(skh1 + 588);
    const auto *skh1_591 = buffer.data(skh1 + 591);
    const auto *skh1_594 = buffer.data(skh1 + 594);
    const auto *skh1_603 = buffer.data(skh1 + 603);
    const auto *skh1_605 = buffer.data(skh1 + 605);
    const auto *skh1_606 = buffer.data(skh1 + 606);
    const auto *skh1_708 = buffer.data(skh1 + 708);
    const auto *skh1_710 = buffer.data(skh1 + 710);
    const auto *skh1_711 = buffer.data(skh1 + 711);
    const auto *skh1_713 = buffer.data(skh1 + 713);
    const auto *skh1_717 = buffer.data(skh1 + 717);
    const auto *skh1_720 = buffer.data(skh1 + 720);
    const auto *skh1_729 = buffer.data(skh1 + 729);
    const auto *skh1_731 = buffer.data(skh1 + 731);
    const auto *skh1_732 = buffer.data(skh1 + 732);
    const auto *skh1_734 = buffer.data(skh1 + 734);
    const auto *skh1_735 = buffer.data(skh1 + 735);
    const auto *skh1_738 = buffer.data(skh1 + 738);
    const auto *skh1_740 = buffer.data(skh1 + 740);
    const auto *skh1_741 = buffer.data(skh1 + 741);
    const auto *skh1_744 = buffer.data(skh1 + 744);
    const auto *skh1_750 = buffer.data(skh1 + 750);
    const auto *skh1_752 = buffer.data(skh1 + 752);
    const auto *skh1_753 = buffer.data(skh1 + 753);
    const auto *skh1_755 = buffer.data(skh1 + 755);

    const auto *slf0_360 = buffer.data(slf0 + 360);
    const auto *slf0_363 = buffer.data(slf0 + 363);
    const auto *slf0_365 = buffer.data(slf0 + 365);
    const auto *slf0_366 = buffer.data(slf0 + 366);
    const auto *slf0_368 = buffer.data(slf0 + 368);
    const auto *slf0_369 = buffer.data(slf0 + 369);
    const auto *slf0_375 = buffer.data(slf0 + 375);
    const auto *slf0_379 = buffer.data(slf0 + 379);
    const auto *slf0_380 = buffer.data(slf0 + 380);
    const auto *slf0_383 = buffer.data(slf0 + 383);
    const auto *slf0_385 = buffer.data(slf0 + 385);
    const auto *slf0_386 = buffer.data(slf0 + 386);
    const auto *slf0_388 = buffer.data(slf0 + 388);
    const auto *slf0_389 = buffer.data(slf0 + 389);
    const auto *slf0_390 = buffer.data(slf0 + 390);
    const auto *slf0_393 = buffer.data(slf0 + 393);
    const auto *slf0_395 = buffer.data(slf0 + 395);
    const auto *slf0_396 = buffer.data(slf0 + 396);
    const auto *slf0_399 = buffer.data(slf0 + 399);

    const auto *slf1_360 = buffer.data(slf1 + 360);
    const auto *slf1_363 = buffer.data(slf1 + 363);
    const auto *slf1_365 = buffer.data(slf1 + 365);
    const auto *slf1_366 = buffer.data(slf1 + 366);
    const auto *slf1_368 = buffer.data(slf1 + 368);
    const auto *slf1_369 = buffer.data(slf1 + 369);
    const auto *slf1_375 = buffer.data(slf1 + 375);
    const auto *slf1_379 = buffer.data(slf1 + 379);
    const auto *slf1_380 = buffer.data(slf1 + 380);
    const auto *slf1_383 = buffer.data(slf1 + 383);
    const auto *slf1_385 = buffer.data(slf1 + 385);
    const auto *slf1_386 = buffer.data(slf1 + 386);
    const auto *slf1_388 = buffer.data(slf1 + 388);
    const auto *slf1_389 = buffer.data(slf1 + 389);
    const auto *slf1_390 = buffer.data(slf1 + 390);
    const auto *slf1_393 = buffer.data(slf1 + 393);
    const auto *slf1_395 = buffer.data(slf1 + 395);
    const auto *slf1_396 = buffer.data(slf1 + 396);
    const auto *slf1_399 = buffer.data(slf1 + 399);

    const auto *slg_505 = buffer.data(slg + 505);
    const auto *slg_508 = buffer.data(slg + 508);
    const auto *slg_509 = buffer.data(slg + 509);
    const auto *slg_510 = buffer.data(slg + 510);
    const auto *slg_512 = buffer.data(slg + 512);
    const auto *slg_513 = buffer.data(slg + 513);
    const auto *slg_515 = buffer.data(slg + 515);
    const auto *slg_520 = buffer.data(slg + 520);
    const auto *slg_521 = buffer.data(slg + 521);
    const auto *slg_522 = buffer.data(slg + 522);
    const auto *slg_523 = buffer.data(slg + 523);
    const auto *slg_524 = buffer.data(slg + 524);
    const auto *slg_525 = buffer.data(slg + 525);
    const auto *slg_527 = buffer.data(slg + 527);
    const auto *slg_528 = buffer.data(slg + 528);
    const auto *slg_530 = buffer.data(slg + 530);
    const auto *slg_535 = buffer.data(slg + 535);
    const auto *slg_536 = buffer.data(slg + 536);
    const auto *slg_537 = buffer.data(slg + 537);
    const auto *slg_538 = buffer.data(slg + 538);
    const auto *slg_539 = buffer.data(slg + 539);
    const auto *slg_540 = buffer.data(slg + 540);
    const auto *slg_542 = buffer.data(slg + 542);
    const auto *slg_543 = buffer.data(slg + 543);
    const auto *slg_545 = buffer.data(slg + 545);
    const auto *slg_546 = buffer.data(slg + 546);
    const auto *slg_549 = buffer.data(slg + 549);
    const auto *slg_550 = buffer.data(slg + 550);
    const auto *slg_551 = buffer.data(slg + 551);
    const auto *slg_552 = buffer.data(slg + 552);
    const auto *slg_553 = buffer.data(slg + 553);
    const auto *slg_554 = buffer.data(slg + 554);
    const auto *slg_555 = buffer.data(slg + 555);
    const auto *slg_557 = buffer.data(slg + 557);
    const auto *slg_558 = buffer.data(slg + 558);
    const auto *slg_560 = buffer.data(slg + 560);
    const auto *slg_564 = buffer.data(slg + 564);
    const auto *slg_565 = buffer.data(slg + 565);
    const auto *slg_566 = buffer.data(slg + 566);
    const auto *slg_567 = buffer.data(slg + 567);
    const auto *slg_568 = buffer.data(slg + 568);
    const auto *slg_569 = buffer.data(slg + 569);
    const auto *slg_570 = buffer.data(slg + 570);
    const auto *slg_572 = buffer.data(slg + 572);
    const auto *slg_573 = buffer.data(slg + 573);
    const auto *slg_575 = buffer.data(slg + 575);
    const auto *slg_576 = buffer.data(slg + 576);
    const auto *slg_579 = buffer.data(slg + 579);
    const auto *slg_580 = buffer.data(slg + 580);
    const auto *slg_581 = buffer.data(slg + 581);
    const auto *slg_582 = buffer.data(slg + 582);
    const auto *slg_583 = buffer.data(slg + 583);
    const auto *slg_584 = buffer.data(slg + 584);
    const auto *slg_585 = buffer.data(slg + 585);
    const auto *slg_587 = buffer.data(slg + 587);
    const auto *slg_588 = buffer.data(slg + 588);
    const auto *slg_590 = buffer.data(slg + 590);
    const auto *slg_591 = buffer.data(slg + 591);
    const auto *slg_594 = buffer.data(slg + 594);
    const auto *slg_595 = buffer.data(slg + 595);
    const auto *slg_596 = buffer.data(slg + 596);
    const auto *slg_597 = buffer.data(slg + 597);
    const auto *slg_598 = buffer.data(slg + 598);
    const auto *slg_599 = buffer.data(slg + 599);

#pragma omp simd aligned(t_706, t_707, t_708, t_709, pb_x, pc_x, pc_z, skh0_708, skg_385, \
                         skg_508, skg_509, skh1_708, slg_505, slg_508, \
                         slg_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_706[k] = f_9 * skg_508[k]
                   + f_3 * pc_x[k] * slg_508[k];

        t_707[k] = f_9 * skg_509[k]
                   + f_3 * pc_x[k] * slg_509[k];

        t_708[k] = pb_x[k] * skh0_708[k]
                   - f_8 * pc_x[k] * skh1_708[k];

        t_709[k] = f_14 * skg_385[k]
                   + f_3 * pc_z[k] * slg_505[k];
    }

#pragma omp simd aligned(t_710, t_711, t_712, t_713, pb_x, pc_x, pc_y, skh0_710, skh0_711, \
                         skh0_713, skg_404, skh1_710, skh1_711, skh1_713, \
                         slg_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_710[k] = pb_x[k] * skh0_710[k]
                   - f_8 * pc_x[k] * skh1_710[k];

        t_711[k] = pb_x[k] * skh0_711[k]
                   - f_8 * pc_x[k] * skh1_711[k];

        t_712[k] = f_10 * skg_404[k]
                   + f_3 * pc_y[k] * slg_509[k];

        t_713[k] = pb_x[k] * skh0_713[k]
                   - f_8 * pc_x[k] * skh1_713[k];
    }

#pragma omp simd aligned(t_714, t_715, t_716, pb_y, pc_y, pc_z, skh0_567, skg_390, skg_405, \
                         skh1_567, slg_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_714[k] = pb_y[k] * skh0_567[k]
                   - f_8 * pc_y[k] * skh1_567[k];

        t_715[k] = f_9 * skg_405[k]
                   + f_3 * pc_y[k] * slg_510[k];

        t_716[k] = f_13 * skg_390[k]
                   + f_3 * pc_z[k] * slg_510[k];
    }

#pragma omp simd aligned(t_717, t_718, t_719, pb_x, pb_y, pc_x, pc_y, skh0_572, skh0_717, \
                         skg_407, skg_513, skh1_572, skh1_717, \
                         slg_512 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_717[k] = pb_x[k] * skh0_717[k]
                   + f_11 * skg_513[k]
                   - f_8 * pc_x[k] * skh1_717[k];

        t_718[k] = f_9 * skg_407[k]
                   + f_3 * pc_y[k] * slg_512[k];

        t_719[k] = pb_y[k] * skh0_572[k]
                   - f_8 * pc_y[k] * skh1_572[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, pb_x, pc_x, pc_y, pc_z, skh0_720, skg_393, \
                         skg_410, skg_516, skh1_720, slg_513, slg_515 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = pb_x[k] * skh0_720[k]
                   + f_10 * skg_516[k]
                   - f_8 * pc_x[k] * skh1_720[k];

        t_721[k] = f_13 * skg_393[k]
                   + f_3 * pc_z[k] * slg_513[k];

        t_722[k] = f_9 * skg_410[k]
                   + f_3 * pc_y[k] * slg_515[k];
    }

#pragma omp simd aligned(t_723, t_724, t_725, t_726, pb_y, pc_x, pc_y, skh0_576, skg_520, \
                         skg_521, skg_522, skh1_576, slg_520, slg_521, \
                         slg_522 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_723[k] = pb_y[k] * skh0_576[k]
                   - f_8 * pc_y[k] * skh1_576[k];

        t_724[k] = f_9 * skg_520[k]
                   + f_3 * pc_x[k] * slg_520[k];

        t_725[k] = f_9 * skg_521[k]
                   + f_3 * pc_x[k] * slg_521[k];

        t_726[k] = f_9 * skg_522[k]
                   + f_3 * pc_x[k] * slg_522[k];
    }

#pragma omp simd aligned(t_727, t_728, t_729, t_730, pb_x, pc_x, pc_z, skh0_729, skg_400, \
                         skg_523, skg_524, skh1_729, slg_520, slg_523, \
                         slg_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_727[k] = f_9 * skg_523[k]
                   + f_3 * pc_x[k] * slg_523[k];

        t_728[k] = f_9 * skg_524[k]
                   + f_3 * pc_x[k] * slg_524[k];

        t_729[k] = pb_x[k] * skh0_729[k]
                   - f_8 * pc_x[k] * skh1_729[k];

        t_730[k] = f_13 * skg_400[k]
                   + f_3 * pc_z[k] * slg_520[k];
    }

#pragma omp simd aligned(t_731, t_732, t_733, t_734, pb_x, pc_x, pc_y, skh0_731, skh0_732, \
                         skh0_734, skg_419, skh1_731, skh1_732, skh1_734, \
                         slg_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_731[k] = pb_x[k] * skh0_731[k]
                   - f_8 * pc_x[k] * skh1_731[k];

        t_732[k] = pb_x[k] * skh0_732[k]
                   - f_8 * pc_x[k] * skh1_732[k];

        t_733[k] = f_9 * skg_419[k]
                   + f_3 * pc_y[k] * slg_524[k];

        t_734[k] = pb_x[k] * skh0_734[k]
                   - f_8 * pc_x[k] * skh1_734[k];
    }

#pragma omp simd aligned(t_735, t_736, t_737, t_738, pb_x, pc_x, pc_y, pc_z, skh0_735, \
                         skh0_738, skg_405, skg_525, skg_528, skh1_735, skh1_738, \
                         slg_525 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_735[k] = pb_x[k] * skh0_735[k]
                   + f_14 * skg_525[k]
                   - f_8 * pc_x[k] * skh1_735[k];

        t_736[k] = f_3 * pc_y[k] * slg_525[k];

        t_737[k] = f_12 * skg_405[k]
                   + f_3 * pc_z[k] * slg_525[k];

        t_738[k] = pb_x[k] * skh0_738[k]
                   + f_11 * skg_528[k]
                   - f_8 * pc_x[k] * skh1_738[k];
    }

#pragma omp simd aligned(t_739, t_740, t_741, pb_x, pc_x, pc_y, skh0_740, skh0_741, skg_530, \
                         skg_531, skh1_740, skh1_741, slg_527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_739[k] = f_3 * pc_y[k] * slg_527[k];

        t_740[k] = pb_x[k] * skh0_740[k]
                   + f_11 * skg_530[k]
                   - f_8 * pc_x[k] * skh1_740[k];

        t_741[k] = pb_x[k] * skh0_741[k]
                   + f_10 * skg_531[k]
                   - f_8 * pc_x[k] * skh1_741[k];
    }

#pragma omp simd aligned(t_742, t_743, t_744, t_745, pb_x, pc_x, pc_y, pc_z, skh0_744, \
                         skg_408, skg_534, skg_535, skh1_744, slg_528, slg_530, \
                         slg_535 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_742[k] = f_12 * skg_408[k]
                   + f_3 * pc_z[k] * slg_528[k];

        t_743[k] = f_3 * pc_y[k] * slg_530[k];

        t_744[k] = pb_x[k] * skh0_744[k]
                   + f_10 * skg_534[k]
                   - f_8 * pc_x[k] * skh1_744[k];

        t_745[k] = f_9 * skg_535[k]
                   + f_3 * pc_x[k] * slg_535[k];
    }

#pragma omp simd aligned(t_746, t_747, t_748, t_749, pc_x, skg_536, skg_537, skg_538, skg_539, \
                         slg_536, slg_537, slg_538, slg_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_746[k] = f_9 * skg_536[k]
                   + f_3 * pc_x[k] * slg_536[k];

        t_747[k] = f_9 * skg_537[k]
                   + f_3 * pc_x[k] * slg_537[k];

        t_748[k] = f_9 * skg_538[k]
                   + f_3 * pc_x[k] * slg_538[k];

        t_749[k] = f_9 * skg_539[k]
                   + f_3 * pc_x[k] * slg_539[k];
    }

#pragma omp simd aligned(t_750, t_751, t_752, t_753, pb_x, pc_x, pc_z, skh0_750, skh0_752, \
                         skh0_753, skg_415, skh1_750, skh1_752, skh1_753, \
                         slg_535 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_750[k] = pb_x[k] * skh0_750[k]
                   - f_8 * pc_x[k] * skh1_750[k];

        t_751[k] = f_12 * skg_415[k]
                   + f_3 * pc_z[k] * slg_535[k];

        t_752[k] = pb_x[k] * skh0_752[k]
                   - f_8 * pc_x[k] * skh1_752[k];

        t_753[k] = pb_x[k] * skh0_753[k]
                   - f_8 * pc_x[k] * skh1_753[k];
    }

#pragma omp simd aligned(t_754, t_755, t_756, t_757, t_758, pb_x, pc_x, pc_y, pc_z, skh0_755, \
                         skg_420, skh1_755, slf0_360, slf1_360, slg_539, \
                         slg_540 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_754[k] = f_3 * pc_y[k] * slg_539[k];

        t_755[k] = pb_x[k] * skh0_755[k]
                   - f_8 * pc_x[k] * skh1_755[k];

        t_756[k] = f_1 * slf0_360[k]
                   - f_2 * slf1_360[k]
                   + f_3 * pc_x[k] * slg_540[k];

        t_757[k] = f_0 * skg_420[k]
                   + f_3 * pc_y[k] * slg_540[k];

        t_758[k] = f_3 * pc_z[k] * slg_540[k];
    }

#pragma omp simd aligned(t_759, t_760, t_761, pc_x, pc_y, skg_422, slf0_363, slf0_365, \
                         slf1_363, slf1_365, slg_542, slg_543, \
                         slg_545 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_759[k] = f_4 * slf0_363[k]
                   - f_5 * slf1_363[k]
                   + f_3 * pc_x[k] * slg_543[k];

        t_760[k] = f_0 * skg_422[k]
                   + f_3 * pc_y[k] * slg_542[k];

        t_761[k] = f_4 * slf0_365[k]
                   - f_5 * slf1_365[k]
                   + f_3 * pc_x[k] * slg_545[k];
    }

#pragma omp simd aligned(t_762, t_763, t_764, t_765, pc_x, pc_y, pc_z, skg_425, slf0_366, \
                         slf0_369, slf1_366, slf1_369, slg_543, slg_545, slg_546, \
                         slg_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_762[k] = f_6 * slf0_366[k]
                   - f_7 * slf1_366[k]
                   + f_3 * pc_x[k] * slg_546[k];

        t_763[k] = f_3 * pc_z[k] * slg_543[k];

        t_764[k] = f_0 * skg_425[k]
                   + f_3 * pc_y[k] * slg_545[k];

        t_765[k] = f_6 * slf0_369[k]
                   - f_7 * slf1_369[k]
                   + f_3 * pc_x[k] * slg_549[k];
    }

#pragma omp simd aligned(t_766, t_767, t_768, t_769, t_770, t_771, pc_x, pc_y, skg_430, \
                         slf0_366, slf1_366, slg_550, slg_551, slg_552, slg_553, \
                         slg_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_766[k] = f_3 * pc_x[k] * slg_550[k];

        t_767[k] = f_3 * pc_x[k] * slg_551[k];

        t_768[k] = f_3 * pc_x[k] * slg_552[k];

        t_769[k] = f_3 * pc_x[k] * slg_553[k];

        t_770[k] = f_3 * pc_x[k] * slg_554[k];

        t_771[k] = f_0 * skg_430[k]
                   + f_1 * slf0_366[k]
                   - f_2 * slf1_366[k]
                   + f_3 * pc_y[k] * slg_550[k];
    }

#pragma omp simd aligned(t_772, t_773, t_774, pc_y, pc_z, skg_432, skg_433, slf0_368, \
                         slf0_369, slf1_368, slf1_369, slg_550, slg_552, \
                         slg_553 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_772[k] = f_3 * pc_z[k] * slg_550[k];

        t_773[k] = f_0 * skg_432[k]
                   + f_4 * slf0_368[k]
                   - f_5 * slf1_368[k]
                   + f_3 * pc_y[k] * slg_552[k];

        t_774[k] = f_0 * skg_433[k]
                   + f_6 * slf0_369[k]
                   - f_7 * slf1_369[k]
                   + f_3 * pc_y[k] * slg_553[k];
    }

#pragma omp simd aligned(t_775, t_776, t_777, t_778, pb_z, pc_y, pc_z, skh0_588, skg_434, \
                         skg_435, skh1_588, slf0_369, slf1_369, slg_554, \
                         slg_555 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_775[k] = f_0 * skg_434[k]
                   + f_3 * pc_y[k] * slg_554[k];

        t_776[k] = f_1 * slf0_369[k]
                   - f_2 * slf1_369[k]
                   + f_3 * pc_z[k] * slg_554[k];

        t_777[k] = pb_z[k] * skh0_588[k]
                   - f_8 * pc_z[k] * skh1_588[k];

        t_778[k] = f_12 * skg_435[k]
                   + f_3 * pc_y[k] * slg_555[k];
    }

#pragma omp simd aligned(t_779, t_780, t_781, pb_z, pc_y, pc_z, skh0_591, skg_420, skg_437, \
                         skh1_591, slg_555, slg_557 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_779[k] = f_9 * skg_420[k]
                   + f_3 * pc_z[k] * slg_555[k];

        t_780[k] = pb_z[k] * skh0_591[k]
                   - f_8 * pc_z[k] * skh1_591[k];

        t_781[k] = f_12 * skg_437[k]
                   + f_3 * pc_y[k] * slg_557[k];
    }

#pragma omp simd aligned(t_782, t_783, t_784, t_785, pb_z, pc_x, pc_y, pc_z, skh0_594, \
                         skg_423, skg_440, skh1_594, slf0_375, slf1_375, slg_558, \
                         slg_560 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_782[k] = f_4 * slf0_375[k]
                   - f_5 * slf1_375[k]
                   + f_3 * pc_x[k] * slg_560[k];

        t_783[k] = pb_z[k] * skh0_594[k]
                   - f_8 * pc_z[k] * skh1_594[k];

        t_784[k] = f_9 * skg_423[k]
                   + f_3 * pc_z[k] * slg_558[k];

        t_785[k] = f_12 * skg_440[k]
                   + f_3 * pc_y[k] * slg_560[k];
    }

#pragma omp simd aligned(t_786, t_787, t_788, t_789, t_790, t_791, pc_x, slf0_379, slf1_379, \
                         slg_564, slg_565, slg_566, slg_567, slg_568, \
                         slg_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_786[k] = f_6 * slf0_379[k]
                   - f_7 * slf1_379[k]
                   + f_3 * pc_x[k] * slg_564[k];

        t_787[k] = f_3 * pc_x[k] * slg_565[k];

        t_788[k] = f_3 * pc_x[k] * slg_566[k];

        t_789[k] = f_3 * pc_x[k] * slg_567[k];

        t_790[k] = f_3 * pc_x[k] * slg_568[k];

        t_791[k] = f_3 * pc_x[k] * slg_569[k];
    }

#pragma omp simd aligned(t_792, t_793, t_794, t_795, pb_z, pc_z, skh0_603, skh0_605, skh0_606, \
                         skg_430, skg_431, skg_432, skh1_603, skh1_605, skh1_606, \
                         slg_565 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_792[k] = pb_z[k] * skh0_603[k]
                   - f_8 * pc_z[k] * skh1_603[k];

        t_793[k] = f_9 * skg_430[k]
                   + f_3 * pc_z[k] * slg_565[k];

        t_794[k] = pb_z[k] * skh0_605[k]
                   + f_10 * skg_431[k]
                   - f_8 * pc_z[k] * skh1_605[k];

        t_795[k] = pb_z[k] * skh0_606[k]
                   + f_11 * skg_432[k]
                   - f_8 * pc_z[k] * skh1_606[k];
    }

#pragma omp simd aligned(t_796, t_797, t_798, t_799, pc_x, pc_y, pc_z, skg_434, skg_449, \
                         skg_450, slf0_379, slf0_380, slf1_379, slf1_380, slg_569, \
                         slg_570 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_796[k] = f_12 * skg_449[k]
                   + f_3 * pc_y[k] * slg_569[k];

        t_797[k] = f_9 * skg_434[k]
                   + f_1 * slf0_379[k]
                   - f_2 * slf1_379[k]
                   + f_3 * pc_z[k] * slg_569[k];

        t_798[k] = f_1 * slf0_380[k]
                   - f_2 * slf1_380[k]
                   + f_3 * pc_x[k] * slg_570[k];

        t_799[k] = f_13 * skg_450[k]
                   + f_3 * pc_y[k] * slg_570[k];
    }

#pragma omp simd aligned(t_800, t_801, t_802, pc_x, pc_y, pc_z, skg_435, skg_452, slf0_383, \
                         slf1_383, slg_570, slg_572, slg_573 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_800[k] = f_10 * skg_435[k]
                   + f_3 * pc_z[k] * slg_570[k];

        t_801[k] = f_4 * slf0_383[k]
                   - f_5 * slf1_383[k]
                   + f_3 * pc_x[k] * slg_573[k];

        t_802[k] = f_13 * skg_452[k]
                   + f_3 * pc_y[k] * slg_572[k];
    }

#pragma omp simd aligned(t_803, t_804, t_805, t_806, pc_x, pc_y, pc_z, skg_438, skg_455, \
                         slf0_385, slf0_386, slf1_385, slf1_386, slg_573, slg_575, \
                         slg_576 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_803[k] = f_4 * slf0_385[k]
                   - f_5 * slf1_385[k]
                   + f_3 * pc_x[k] * slg_575[k];

        t_804[k] = f_6 * slf0_386[k]
                   - f_7 * slf1_386[k]
                   + f_3 * pc_x[k] * slg_576[k];

        t_805[k] = f_10 * skg_438[k]
                   + f_3 * pc_z[k] * slg_573[k];

        t_806[k] = f_13 * skg_455[k]
                   + f_3 * pc_y[k] * slg_575[k];
    }

#pragma omp simd aligned(t_807, t_808, t_809, t_810, t_811, t_812, pc_x, slf0_389, slf1_389, \
                         slg_579, slg_580, slg_581, slg_582, slg_583, \
                         slg_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_807[k] = f_6 * slf0_389[k]
                   - f_7 * slf1_389[k]
                   + f_3 * pc_x[k] * slg_579[k];

        t_808[k] = f_3 * pc_x[k] * slg_580[k];

        t_809[k] = f_3 * pc_x[k] * slg_581[k];

        t_810[k] = f_3 * pc_x[k] * slg_582[k];

        t_811[k] = f_3 * pc_x[k] * slg_583[k];

        t_812[k] = f_3 * pc_x[k] * slg_584[k];
    }

#pragma omp simd aligned(t_813, t_814, t_815, pc_y, pc_z, skg_445, skg_460, skg_462, slf0_386, \
                         slf0_388, slf1_386, slf1_388, slg_580, \
                         slg_582 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_813[k] = f_13 * skg_460[k]
                   + f_1 * slf0_386[k]
                   - f_2 * slf1_386[k]
                   + f_3 * pc_y[k] * slg_580[k];

        t_814[k] = f_10 * skg_445[k]
                   + f_3 * pc_z[k] * slg_580[k];

        t_815[k] = f_13 * skg_462[k]
                   + f_4 * slf0_388[k]
                   - f_5 * slf1_388[k]
                   + f_3 * pc_y[k] * slg_582[k];
    }

#pragma omp simd aligned(t_816, t_817, t_818, pc_y, pc_z, skg_449, skg_463, skg_464, slf0_389, \
                         slf1_389, slg_583, slg_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_816[k] = f_13 * skg_463[k]
                   + f_6 * slf0_389[k]
                   - f_7 * slf1_389[k]
                   + f_3 * pc_y[k] * slg_583[k];

        t_817[k] = f_13 * skg_464[k]
                   + f_3 * pc_y[k] * slg_584[k];

        t_818[k] = f_10 * skg_449[k]
                   + f_1 * slf0_389[k]
                   - f_2 * slf1_389[k]
                   + f_3 * pc_z[k] * slg_584[k];
    }

#pragma omp simd aligned(t_819, t_820, t_821, t_822, pc_x, pc_y, pc_z, skg_450, skg_465, \
                         slf0_390, slf0_393, slf1_390, slf1_393, slg_585, \
                         slg_588 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_819[k] = f_1 * slf0_390[k]
                   - f_2 * slf1_390[k]
                   + f_3 * pc_x[k] * slg_585[k];

        t_820[k] = f_14 * skg_465[k]
                   + f_3 * pc_y[k] * slg_585[k];

        t_821[k] = f_11 * skg_450[k]
                   + f_3 * pc_z[k] * slg_585[k];

        t_822[k] = f_4 * slf0_393[k]
                   - f_5 * slf1_393[k]
                   + f_3 * pc_x[k] * slg_588[k];
    }

#pragma omp simd aligned(t_823, t_824, t_825, pc_x, pc_y, skg_467, slf0_395, slf0_396, \
                         slf1_395, slf1_396, slg_587, slg_590, \
                         slg_591 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_823[k] = f_14 * skg_467[k]
                   + f_3 * pc_y[k] * slg_587[k];

        t_824[k] = f_4 * slf0_395[k]
                   - f_5 * slf1_395[k]
                   + f_3 * pc_x[k] * slg_590[k];

        t_825[k] = f_6 * slf0_396[k]
                   - f_7 * slf1_396[k]
                   + f_3 * pc_x[k] * slg_591[k];
    }

#pragma omp simd aligned(t_826, t_827, t_828, t_829, pc_x, pc_y, pc_z, skg_453, skg_470, \
                         slf0_399, slf1_399, slg_588, slg_590, slg_594, \
                         slg_595 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_826[k] = f_11 * skg_453[k]
                   + f_3 * pc_z[k] * slg_588[k];

        t_827[k] = f_14 * skg_470[k]
                   + f_3 * pc_y[k] * slg_590[k];

        t_828[k] = f_6 * slf0_399[k]
                   - f_7 * slf1_399[k]
                   + f_3 * pc_x[k] * slg_594[k];

        t_829[k] = f_3 * pc_x[k] * slg_595[k];
    }

#pragma omp simd aligned(t_830, t_831, t_832, t_833, t_834, pc_x, pc_y, skg_475, slf0_396, \
                         slf1_396, slg_595, slg_596, slg_597, slg_598, \
                         slg_599 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_830[k] = f_3 * pc_x[k] * slg_596[k];

        t_831[k] = f_3 * pc_x[k] * slg_597[k];

        t_832[k] = f_3 * pc_x[k] * slg_598[k];

        t_833[k] = f_3 * pc_x[k] * slg_599[k];

        t_834[k] = f_14 * skg_475[k]
                   + f_1 * slf0_396[k]
                   - f_2 * slf1_396[k]
                   + f_3 * pc_y[k] * slg_595[k];
    }
}

static auto
compute_prim_slh_three_center_electron_repulsion_0_piece7(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t skh0,
                                                          const size_t skg, const size_t skh1,
                                                          const size_t slf0, const size_t slf1,
                                                          const size_t slg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
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
    const auto f_12 = 3.5 / q;
    const auto f_13 = 3.0 / q;
    const auto f_14 = 2.5 / q;
    const auto f_15 = 2.0 / q;

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

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *skh0_735 = buffer.data(skh0 + 735);
    const auto *skh0_740 = buffer.data(skh0 + 740);
    const auto *skh0_744 = buffer.data(skh0 + 744);
    const auto *skh0_750 = buffer.data(skh0 + 750);
    const auto *skh0_752 = buffer.data(skh0 + 752);
    const auto *skh0_753 = buffer.data(skh0 + 753);
    const auto *skh0_755 = buffer.data(skh0 + 755);

    const auto *skg_460 = buffer.data(skg + 460);
    const auto *skg_464 = buffer.data(skg + 464);
    const auto *skg_465 = buffer.data(skg + 465);
    const auto *skg_468 = buffer.data(skg + 468);
    const auto *skg_475 = buffer.data(skg + 475);
    const auto *skg_477 = buffer.data(skg + 477);
    const auto *skg_478 = buffer.data(skg + 478);
    const auto *skg_479 = buffer.data(skg + 479);
    const auto *skg_480 = buffer.data(skg + 480);
    const auto *skg_482 = buffer.data(skg + 482);
    const auto *skg_483 = buffer.data(skg + 483);
    const auto *skg_485 = buffer.data(skg + 485);
    const auto *skg_490 = buffer.data(skg + 490);
    const auto *skg_492 = buffer.data(skg + 492);
    const auto *skg_493 = buffer.data(skg + 493);
    const auto *skg_494 = buffer.data(skg + 494);
    const auto *skg_495 = buffer.data(skg + 495);
    const auto *skg_497 = buffer.data(skg + 497);
    const auto *skg_498 = buffer.data(skg + 498);
    const auto *skg_500 = buffer.data(skg + 500);
    const auto *skg_505 = buffer.data(skg + 505);
    const auto *skg_507 = buffer.data(skg + 507);
    const auto *skg_508 = buffer.data(skg + 508);
    const auto *skg_509 = buffer.data(skg + 509);
    const auto *skg_510 = buffer.data(skg + 510);
    const auto *skg_512 = buffer.data(skg + 512);
    const auto *skg_513 = buffer.data(skg + 513);
    const auto *skg_515 = buffer.data(skg + 515);
    const auto *skg_520 = buffer.data(skg + 520);
    const auto *skg_522 = buffer.data(skg + 522);
    const auto *skg_523 = buffer.data(skg + 523);
    const auto *skg_524 = buffer.data(skg + 524);
    const auto *skg_525 = buffer.data(skg + 525);
    const auto *skg_527 = buffer.data(skg + 527);
    const auto *skg_528 = buffer.data(skg + 528);
    const auto *skg_530 = buffer.data(skg + 530);
    const auto *skg_535 = buffer.data(skg + 535);
    const auto *skg_537 = buffer.data(skg + 537);
    const auto *skg_538 = buffer.data(skg + 538);
    const auto *skg_539 = buffer.data(skg + 539);

    const auto *skh1_735 = buffer.data(skh1 + 735);
    const auto *skh1_740 = buffer.data(skh1 + 740);
    const auto *skh1_744 = buffer.data(skh1 + 744);
    const auto *skh1_750 = buffer.data(skh1 + 750);
    const auto *skh1_752 = buffer.data(skh1 + 752);
    const auto *skh1_753 = buffer.data(skh1 + 753);
    const auto *skh1_755 = buffer.data(skh1 + 755);

    const auto *slf0_398 = buffer.data(slf0 + 398);
    const auto *slf0_399 = buffer.data(slf0 + 399);
    const auto *slf0_400 = buffer.data(slf0 + 400);
    const auto *slf0_403 = buffer.data(slf0 + 403);
    const auto *slf0_405 = buffer.data(slf0 + 405);
    const auto *slf0_406 = buffer.data(slf0 + 406);
    const auto *slf0_408 = buffer.data(slf0 + 408);
    const auto *slf0_409 = buffer.data(slf0 + 409);
    const auto *slf0_410 = buffer.data(slf0 + 410);
    const auto *slf0_413 = buffer.data(slf0 + 413);
    const auto *slf0_415 = buffer.data(slf0 + 415);
    const auto *slf0_416 = buffer.data(slf0 + 416);
    const auto *slf0_418 = buffer.data(slf0 + 418);
    const auto *slf0_419 = buffer.data(slf0 + 419);
    const auto *slf0_420 = buffer.data(slf0 + 420);
    const auto *slf0_423 = buffer.data(slf0 + 423);
    const auto *slf0_425 = buffer.data(slf0 + 425);
    const auto *slf0_426 = buffer.data(slf0 + 426);
    const auto *slf0_428 = buffer.data(slf0 + 428);
    const auto *slf0_429 = buffer.data(slf0 + 429);
    const auto *slf0_433 = buffer.data(slf0 + 433);
    const auto *slf0_436 = buffer.data(slf0 + 436);
    const auto *slf0_440 = buffer.data(slf0 + 440);
    const auto *slf0_443 = buffer.data(slf0 + 443);
    const auto *slf0_445 = buffer.data(slf0 + 445);
    const auto *slf0_446 = buffer.data(slf0 + 446);
    const auto *slf0_448 = buffer.data(slf0 + 448);
    const auto *slf0_449 = buffer.data(slf0 + 449);

    const auto *slf1_398 = buffer.data(slf1 + 398);
    const auto *slf1_399 = buffer.data(slf1 + 399);
    const auto *slf1_400 = buffer.data(slf1 + 400);
    const auto *slf1_403 = buffer.data(slf1 + 403);
    const auto *slf1_405 = buffer.data(slf1 + 405);
    const auto *slf1_406 = buffer.data(slf1 + 406);
    const auto *slf1_408 = buffer.data(slf1 + 408);
    const auto *slf1_409 = buffer.data(slf1 + 409);
    const auto *slf1_410 = buffer.data(slf1 + 410);
    const auto *slf1_413 = buffer.data(slf1 + 413);
    const auto *slf1_415 = buffer.data(slf1 + 415);
    const auto *slf1_416 = buffer.data(slf1 + 416);
    const auto *slf1_418 = buffer.data(slf1 + 418);
    const auto *slf1_419 = buffer.data(slf1 + 419);
    const auto *slf1_420 = buffer.data(slf1 + 420);
    const auto *slf1_423 = buffer.data(slf1 + 423);
    const auto *slf1_425 = buffer.data(slf1 + 425);
    const auto *slf1_426 = buffer.data(slf1 + 426);
    const auto *slf1_428 = buffer.data(slf1 + 428);
    const auto *slf1_429 = buffer.data(slf1 + 429);
    const auto *slf1_433 = buffer.data(slf1 + 433);
    const auto *slf1_436 = buffer.data(slf1 + 436);
    const auto *slf1_440 = buffer.data(slf1 + 440);
    const auto *slf1_443 = buffer.data(slf1 + 443);
    const auto *slf1_445 = buffer.data(slf1 + 445);
    const auto *slf1_446 = buffer.data(slf1 + 446);
    const auto *slf1_448 = buffer.data(slf1 + 448);
    const auto *slf1_449 = buffer.data(slf1 + 449);

    const auto *slg_595 = buffer.data(slg + 595);
    const auto *slg_597 = buffer.data(slg + 597);
    const auto *slg_598 = buffer.data(slg + 598);
    const auto *slg_599 = buffer.data(slg + 599);
    const auto *slg_600 = buffer.data(slg + 600);
    const auto *slg_602 = buffer.data(slg + 602);
    const auto *slg_603 = buffer.data(slg + 603);
    const auto *slg_605 = buffer.data(slg + 605);
    const auto *slg_606 = buffer.data(slg + 606);
    const auto *slg_609 = buffer.data(slg + 609);
    const auto *slg_610 = buffer.data(slg + 610);
    const auto *slg_611 = buffer.data(slg + 611);
    const auto *slg_612 = buffer.data(slg + 612);
    const auto *slg_613 = buffer.data(slg + 613);
    const auto *slg_614 = buffer.data(slg + 614);
    const auto *slg_615 = buffer.data(slg + 615);
    const auto *slg_617 = buffer.data(slg + 617);
    const auto *slg_618 = buffer.data(slg + 618);
    const auto *slg_620 = buffer.data(slg + 620);
    const auto *slg_621 = buffer.data(slg + 621);
    const auto *slg_624 = buffer.data(slg + 624);
    const auto *slg_625 = buffer.data(slg + 625);
    const auto *slg_626 = buffer.data(slg + 626);
    const auto *slg_627 = buffer.data(slg + 627);
    const auto *slg_628 = buffer.data(slg + 628);
    const auto *slg_629 = buffer.data(slg + 629);
    const auto *slg_630 = buffer.data(slg + 630);
    const auto *slg_632 = buffer.data(slg + 632);
    const auto *slg_633 = buffer.data(slg + 633);
    const auto *slg_635 = buffer.data(slg + 635);
    const auto *slg_636 = buffer.data(slg + 636);
    const auto *slg_639 = buffer.data(slg + 639);
    const auto *slg_640 = buffer.data(slg + 640);
    const auto *slg_641 = buffer.data(slg + 641);
    const auto *slg_642 = buffer.data(slg + 642);
    const auto *slg_643 = buffer.data(slg + 643);
    const auto *slg_644 = buffer.data(slg + 644);
    const auto *slg_645 = buffer.data(slg + 645);
    const auto *slg_647 = buffer.data(slg + 647);
    const auto *slg_648 = buffer.data(slg + 648);
    const auto *slg_650 = buffer.data(slg + 650);
    const auto *slg_651 = buffer.data(slg + 651);
    const auto *slg_655 = buffer.data(slg + 655);
    const auto *slg_656 = buffer.data(slg + 656);
    const auto *slg_657 = buffer.data(slg + 657);
    const auto *slg_658 = buffer.data(slg + 658);
    const auto *slg_659 = buffer.data(slg + 659);
    const auto *slg_660 = buffer.data(slg + 660);
    const auto *slg_662 = buffer.data(slg + 662);
    const auto *slg_663 = buffer.data(slg + 663);
    const auto *slg_665 = buffer.data(slg + 665);
    const auto *slg_666 = buffer.data(slg + 666);
    const auto *slg_669 = buffer.data(slg + 669);
    const auto *slg_670 = buffer.data(slg + 670);
    const auto *slg_671 = buffer.data(slg + 671);
    const auto *slg_672 = buffer.data(slg + 672);
    const auto *slg_673 = buffer.data(slg + 673);
    const auto *slg_674 = buffer.data(slg + 674);

#pragma omp simd aligned(t_835, t_836, t_837, pc_y, pc_z, skg_460, skg_477, skg_478, slf0_398, \
                         slf0_399, slf1_398, slf1_399, slg_595, slg_597, \
                         slg_598 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_835[k] = f_11 * skg_460[k]
                   + f_3 * pc_z[k] * slg_595[k];

        t_836[k] = f_14 * skg_477[k]
                   + f_4 * slf0_398[k]
                   - f_5 * slf1_398[k]
                   + f_3 * pc_y[k] * slg_597[k];

        t_837[k] = f_14 * skg_478[k]
                   + f_6 * slf0_399[k]
                   - f_7 * slf1_399[k]
                   + f_3 * pc_y[k] * slg_598[k];
    }

#pragma omp simd aligned(t_838, t_839, t_840, t_841, pc_x, pc_y, pc_z, skg_464, skg_479, \
                         skg_480, slf0_399, slf0_400, slf1_399, slf1_400, slg_599, \
                         slg_600 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_838[k] = f_14 * skg_479[k]
                   + f_3 * pc_y[k] * slg_599[k];

        t_839[k] = f_11 * skg_464[k]
                   + f_1 * slf0_399[k]
                   - f_2 * slf1_399[k]
                   + f_3 * pc_z[k] * slg_599[k];

        t_840[k] = f_1 * slf0_400[k]
                   - f_2 * slf1_400[k]
                   + f_3 * pc_x[k] * slg_600[k];

        t_841[k] = f_15 * skg_480[k]
                   + f_3 * pc_y[k] * slg_600[k];
    }

#pragma omp simd aligned(t_842, t_843, t_844, pc_x, pc_y, pc_z, skg_465, skg_482, slf0_403, \
                         slf1_403, slg_600, slg_602, slg_603 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_842[k] = f_15 * skg_465[k]
                   + f_3 * pc_z[k] * slg_600[k];

        t_843[k] = f_4 * slf0_403[k]
                   - f_5 * slf1_403[k]
                   + f_3 * pc_x[k] * slg_603[k];

        t_844[k] = f_15 * skg_482[k]
                   + f_3 * pc_y[k] * slg_602[k];
    }

#pragma omp simd aligned(t_845, t_846, t_847, t_848, pc_x, pc_y, pc_z, skg_468, skg_485, \
                         slf0_405, slf0_406, slf1_405, slf1_406, slg_603, slg_605, \
                         slg_606 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_845[k] = f_4 * slf0_405[k]
                   - f_5 * slf1_405[k]
                   + f_3 * pc_x[k] * slg_605[k];

        t_846[k] = f_6 * slf0_406[k]
                   - f_7 * slf1_406[k]
                   + f_3 * pc_x[k] * slg_606[k];

        t_847[k] = f_15 * skg_468[k]
                   + f_3 * pc_z[k] * slg_603[k];

        t_848[k] = f_15 * skg_485[k]
                   + f_3 * pc_y[k] * slg_605[k];
    }

#pragma omp simd aligned(t_849, t_850, t_851, t_852, t_853, t_854, pc_x, slf0_409, slf1_409, \
                         slg_609, slg_610, slg_611, slg_612, slg_613, \
                         slg_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_849[k] = f_6 * slf0_409[k]
                   - f_7 * slf1_409[k]
                   + f_3 * pc_x[k] * slg_609[k];

        t_850[k] = f_3 * pc_x[k] * slg_610[k];

        t_851[k] = f_3 * pc_x[k] * slg_611[k];

        t_852[k] = f_3 * pc_x[k] * slg_612[k];

        t_853[k] = f_3 * pc_x[k] * slg_613[k];

        t_854[k] = f_3 * pc_x[k] * slg_614[k];
    }

#pragma omp simd aligned(t_855, t_856, t_857, pc_y, pc_z, skg_475, skg_490, skg_492, slf0_406, \
                         slf0_408, slf1_406, slf1_408, slg_610, \
                         slg_612 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_855[k] = f_15 * skg_490[k]
                   + f_1 * slf0_406[k]
                   - f_2 * slf1_406[k]
                   + f_3 * pc_y[k] * slg_610[k];

        t_856[k] = f_15 * skg_475[k]
                   + f_3 * pc_z[k] * slg_610[k];

        t_857[k] = f_15 * skg_492[k]
                   + f_4 * slf0_408[k]
                   - f_5 * slf1_408[k]
                   + f_3 * pc_y[k] * slg_612[k];
    }

#pragma omp simd aligned(t_858, t_859, t_860, pc_y, pc_z, skg_479, skg_493, skg_494, slf0_409, \
                         slf1_409, slg_613, slg_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_858[k] = f_15 * skg_493[k]
                   + f_6 * slf0_409[k]
                   - f_7 * slf1_409[k]
                   + f_3 * pc_y[k] * slg_613[k];

        t_859[k] = f_15 * skg_494[k]
                   + f_3 * pc_y[k] * slg_614[k];

        t_860[k] = f_15 * skg_479[k]
                   + f_1 * slf0_409[k]
                   - f_2 * slf1_409[k]
                   + f_3 * pc_z[k] * slg_614[k];
    }

#pragma omp simd aligned(t_861, t_862, t_863, t_864, pc_x, pc_y, pc_z, skg_480, skg_495, \
                         slf0_410, slf0_413, slf1_410, slf1_413, slg_615, \
                         slg_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_861[k] = f_1 * slf0_410[k]
                   - f_2 * slf1_410[k]
                   + f_3 * pc_x[k] * slg_615[k];

        t_862[k] = f_11 * skg_495[k]
                   + f_3 * pc_y[k] * slg_615[k];

        t_863[k] = f_14 * skg_480[k]
                   + f_3 * pc_z[k] * slg_615[k];

        t_864[k] = f_4 * slf0_413[k]
                   - f_5 * slf1_413[k]
                   + f_3 * pc_x[k] * slg_618[k];
    }

#pragma omp simd aligned(t_865, t_866, t_867, pc_x, pc_y, skg_497, slf0_415, slf0_416, \
                         slf1_415, slf1_416, slg_617, slg_620, \
                         slg_621 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_865[k] = f_11 * skg_497[k]
                   + f_3 * pc_y[k] * slg_617[k];

        t_866[k] = f_4 * slf0_415[k]
                   - f_5 * slf1_415[k]
                   + f_3 * pc_x[k] * slg_620[k];

        t_867[k] = f_6 * slf0_416[k]
                   - f_7 * slf1_416[k]
                   + f_3 * pc_x[k] * slg_621[k];
    }

#pragma omp simd aligned(t_868, t_869, t_870, t_871, pc_x, pc_y, pc_z, skg_483, skg_500, \
                         slf0_419, slf1_419, slg_618, slg_620, slg_624, \
                         slg_625 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_868[k] = f_14 * skg_483[k]
                   + f_3 * pc_z[k] * slg_618[k];

        t_869[k] = f_11 * skg_500[k]
                   + f_3 * pc_y[k] * slg_620[k];

        t_870[k] = f_6 * slf0_419[k]
                   - f_7 * slf1_419[k]
                   + f_3 * pc_x[k] * slg_624[k];

        t_871[k] = f_3 * pc_x[k] * slg_625[k];
    }

#pragma omp simd aligned(t_872, t_873, t_874, t_875, t_876, pc_x, pc_y, skg_505, slf0_416, \
                         slf1_416, slg_625, slg_626, slg_627, slg_628, \
                         slg_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_872[k] = f_3 * pc_x[k] * slg_626[k];

        t_873[k] = f_3 * pc_x[k] * slg_627[k];

        t_874[k] = f_3 * pc_x[k] * slg_628[k];

        t_875[k] = f_3 * pc_x[k] * slg_629[k];

        t_876[k] = f_11 * skg_505[k]
                   + f_1 * slf0_416[k]
                   - f_2 * slf1_416[k]
                   + f_3 * pc_y[k] * slg_625[k];
    }

#pragma omp simd aligned(t_877, t_878, t_879, pc_y, pc_z, skg_490, skg_507, skg_508, slf0_418, \
                         slf0_419, slf1_418, slf1_419, slg_625, slg_627, \
                         slg_628 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_877[k] = f_14 * skg_490[k]
                   + f_3 * pc_z[k] * slg_625[k];

        t_878[k] = f_11 * skg_507[k]
                   + f_4 * slf0_418[k]
                   - f_5 * slf1_418[k]
                   + f_3 * pc_y[k] * slg_627[k];

        t_879[k] = f_11 * skg_508[k]
                   + f_6 * slf0_419[k]
                   - f_7 * slf1_419[k]
                   + f_3 * pc_y[k] * slg_628[k];
    }

#pragma omp simd aligned(t_880, t_881, t_882, t_883, pc_x, pc_y, pc_z, skg_494, skg_509, \
                         skg_510, slf0_419, slf0_420, slf1_419, slf1_420, slg_629, \
                         slg_630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_880[k] = f_11 * skg_509[k]
                   + f_3 * pc_y[k] * slg_629[k];

        t_881[k] = f_14 * skg_494[k]
                   + f_1 * slf0_419[k]
                   - f_2 * slf1_419[k]
                   + f_3 * pc_z[k] * slg_629[k];

        t_882[k] = f_1 * slf0_420[k]
                   - f_2 * slf1_420[k]
                   + f_3 * pc_x[k] * slg_630[k];

        t_883[k] = f_10 * skg_510[k]
                   + f_3 * pc_y[k] * slg_630[k];
    }

#pragma omp simd aligned(t_884, t_885, t_886, pc_x, pc_y, pc_z, skg_495, skg_512, slf0_423, \
                         slf1_423, slg_630, slg_632, slg_633 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_884[k] = f_13 * skg_495[k]
                   + f_3 * pc_z[k] * slg_630[k];

        t_885[k] = f_4 * slf0_423[k]
                   - f_5 * slf1_423[k]
                   + f_3 * pc_x[k] * slg_633[k];

        t_886[k] = f_10 * skg_512[k]
                   + f_3 * pc_y[k] * slg_632[k];
    }

#pragma omp simd aligned(t_887, t_888, t_889, t_890, pc_x, pc_y, pc_z, skg_498, skg_515, \
                         slf0_425, slf0_426, slf1_425, slf1_426, slg_633, slg_635, \
                         slg_636 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_887[k] = f_4 * slf0_425[k]
                   - f_5 * slf1_425[k]
                   + f_3 * pc_x[k] * slg_635[k];

        t_888[k] = f_6 * slf0_426[k]
                   - f_7 * slf1_426[k]
                   + f_3 * pc_x[k] * slg_636[k];

        t_889[k] = f_13 * skg_498[k]
                   + f_3 * pc_z[k] * slg_633[k];

        t_890[k] = f_10 * skg_515[k]
                   + f_3 * pc_y[k] * slg_635[k];
    }

#pragma omp simd aligned(t_891, t_892, t_893, t_894, t_895, t_896, pc_x, slf0_429, slf1_429, \
                         slg_639, slg_640, slg_641, slg_642, slg_643, \
                         slg_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_891[k] = f_6 * slf0_429[k]
                   - f_7 * slf1_429[k]
                   + f_3 * pc_x[k] * slg_639[k];

        t_892[k] = f_3 * pc_x[k] * slg_640[k];

        t_893[k] = f_3 * pc_x[k] * slg_641[k];

        t_894[k] = f_3 * pc_x[k] * slg_642[k];

        t_895[k] = f_3 * pc_x[k] * slg_643[k];

        t_896[k] = f_3 * pc_x[k] * slg_644[k];
    }

#pragma omp simd aligned(t_897, t_898, t_899, pc_y, pc_z, skg_505, skg_520, skg_522, slf0_426, \
                         slf0_428, slf1_426, slf1_428, slg_640, \
                         slg_642 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_897[k] = f_10 * skg_520[k]
                   + f_1 * slf0_426[k]
                   - f_2 * slf1_426[k]
                   + f_3 * pc_y[k] * slg_640[k];

        t_898[k] = f_13 * skg_505[k]
                   + f_3 * pc_z[k] * slg_640[k];

        t_899[k] = f_10 * skg_522[k]
                   + f_4 * slf0_428[k]
                   - f_5 * slf1_428[k]
                   + f_3 * pc_y[k] * slg_642[k];
    }

#pragma omp simd aligned(t_900, t_901, t_902, t_903, pb_y, pc_y, pc_z, skh0_735, skg_509, \
                         skg_523, skg_524, skh1_735, slf0_429, slf1_429, slg_643, \
                         slg_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_900[k] = f_10 * skg_523[k]
                   + f_6 * slf0_429[k]
                   - f_7 * slf1_429[k]
                   + f_3 * pc_y[k] * slg_643[k];

        t_901[k] = f_10 * skg_524[k]
                   + f_3 * pc_y[k] * slg_644[k];

        t_902[k] = f_13 * skg_509[k]
                   + f_1 * slf0_429[k]
                   - f_2 * slf1_429[k]
                   + f_3 * pc_z[k] * slg_644[k];

        t_903[k] = pb_y[k] * skh0_735[k]
                   - f_8 * pc_y[k] * skh1_735[k];
    }

#pragma omp simd aligned(t_904, t_905, t_906, t_907, pc_x, pc_y, pc_z, skg_510, skg_525, \
                         skg_527, slf0_433, slf1_433, slg_645, slg_647, \
                         slg_648 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_904[k] = f_9 * skg_525[k]
                   + f_3 * pc_y[k] * slg_645[k];

        t_905[k] = f_12 * skg_510[k]
                   + f_3 * pc_z[k] * slg_645[k];

        t_906[k] = f_4 * slf0_433[k]
                   - f_5 * slf1_433[k]
                   + f_3 * pc_x[k] * slg_648[k];

        t_907[k] = f_9 * skg_527[k]
                   + f_3 * pc_y[k] * slg_647[k];
    }

#pragma omp simd aligned(t_908, t_909, t_910, pb_y, pc_x, pc_y, pc_z, skh0_740, skg_513, \
                         skh1_740, slf0_436, slf1_436, slg_648, \
                         slg_651 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_908[k] = pb_y[k] * skh0_740[k]
                   - f_8 * pc_y[k] * skh1_740[k];

        t_909[k] = f_6 * slf0_436[k]
                   - f_7 * slf1_436[k]
                   + f_3 * pc_x[k] * slg_651[k];

        t_910[k] = f_12 * skg_513[k]
                   + f_3 * pc_z[k] * slg_648[k];
    }

#pragma omp simd aligned(t_911, t_912, t_913, t_914, t_915, pb_y, pc_x, pc_y, skh0_744, \
                         skg_530, skh1_744, slg_650, slg_655, slg_656, \
                         slg_657 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_911[k] = f_9 * skg_530[k]
                   + f_3 * pc_y[k] * slg_650[k];

        t_912[k] = pb_y[k] * skh0_744[k]
                   - f_8 * pc_y[k] * skh1_744[k];

        t_913[k] = f_3 * pc_x[k] * slg_655[k];

        t_914[k] = f_3 * pc_x[k] * slg_656[k];

        t_915[k] = f_3 * pc_x[k] * slg_657[k];
    }

#pragma omp simd aligned(t_916, t_917, t_918, t_919, pb_y, pc_x, pc_y, pc_z, skh0_750, \
                         skg_520, skg_535, skh1_750, slg_655, slg_658, \
                         slg_659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_916[k] = f_3 * pc_x[k] * slg_658[k];

        t_917[k] = f_3 * pc_x[k] * slg_659[k];

        t_918[k] = pb_y[k] * skh0_750[k]
                   + f_14 * skg_535[k]
                   - f_8 * pc_y[k] * skh1_750[k];

        t_919[k] = f_12 * skg_520[k]
                   + f_3 * pc_z[k] * slg_655[k];
    }

#pragma omp simd aligned(t_920, t_921, t_922, t_923, pb_y, pc_y, skh0_752, skh0_753, skh0_755, \
                         skg_537, skg_538, skg_539, skh1_752, skh1_753, skh1_755, \
                         slg_659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_920[k] = pb_y[k] * skh0_752[k]
                   + f_11 * skg_537[k]
                   - f_8 * pc_y[k] * skh1_752[k];

        t_921[k] = pb_y[k] * skh0_753[k]
                   + f_10 * skg_538[k]
                   - f_8 * pc_y[k] * skh1_753[k];

        t_922[k] = f_9 * skg_539[k]
                   + f_3 * pc_y[k] * slg_659[k];

        t_923[k] = pb_y[k] * skh0_755[k]
                   - f_8 * pc_y[k] * skh1_755[k];
    }

#pragma omp simd aligned(t_924, t_925, t_926, t_927, t_928, pc_x, pc_y, pc_z, skg_525, \
                         slf0_440, slf0_443, slf1_440, slf1_443, slg_660, slg_662, \
                         slg_663 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_924[k] = f_1 * slf0_440[k]
                   - f_2 * slf1_440[k]
                   + f_3 * pc_x[k] * slg_660[k];

        t_925[k] = f_3 * pc_y[k] * slg_660[k];

        t_926[k] = f_0 * skg_525[k]
                   + f_3 * pc_z[k] * slg_660[k];

        t_927[k] = f_4 * slf0_443[k]
                   - f_5 * slf1_443[k]
                   + f_3 * pc_x[k] * slg_663[k];

        t_928[k] = f_3 * pc_y[k] * slg_662[k];
    }

#pragma omp simd aligned(t_929, t_930, t_931, t_932, pc_x, pc_y, pc_z, skg_528, slf0_445, \
                         slf0_446, slf1_445, slf1_446, slg_663, slg_665, \
                         slg_666 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_929[k] = f_4 * slf0_445[k]
                   - f_5 * slf1_445[k]
                   + f_3 * pc_x[k] * slg_665[k];

        t_930[k] = f_6 * slf0_446[k]
                   - f_7 * slf1_446[k]
                   + f_3 * pc_x[k] * slg_666[k];

        t_931[k] = f_0 * skg_528[k]
                   + f_3 * pc_z[k] * slg_663[k];

        t_932[k] = f_3 * pc_y[k] * slg_665[k];
    }

#pragma omp simd aligned(t_933, t_934, t_935, t_936, t_937, t_938, pc_x, slf0_449, slf1_449, \
                         slg_669, slg_670, slg_671, slg_672, slg_673, \
                         slg_674 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_933[k] = f_6 * slf0_449[k]
                   - f_7 * slf1_449[k]
                   + f_3 * pc_x[k] * slg_669[k];

        t_934[k] = f_3 * pc_x[k] * slg_670[k];

        t_935[k] = f_3 * pc_x[k] * slg_671[k];

        t_936[k] = f_3 * pc_x[k] * slg_672[k];

        t_937[k] = f_3 * pc_x[k] * slg_673[k];

        t_938[k] = f_3 * pc_x[k] * slg_674[k];
    }

#pragma omp simd aligned(t_939, t_940, t_941, t_942, pc_y, pc_z, skg_535, slf0_446, slf0_448, \
                         slf0_449, slf1_446, slf1_448, slf1_449, slg_670, slg_672, \
                         slg_673 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_939[k] = f_1 * slf0_446[k]
                   - f_2 * slf1_446[k]
                   + f_3 * pc_y[k] * slg_670[k];

        t_940[k] = f_0 * skg_535[k]
                   + f_3 * pc_z[k] * slg_670[k];

        t_941[k] = f_4 * slf0_448[k]
                   - f_5 * slf1_448[k]
                   + f_3 * pc_y[k] * slg_672[k];

        t_942[k] = f_6 * slf0_449[k]
                   - f_7 * slf1_449[k]
                   + f_3 * pc_y[k] * slg_673[k];
    }

#pragma omp simd aligned(t_943, t_944, pc_y, pc_z, skg_539, slf0_449, slf1_449, \
                         slg_674 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_943[k] = f_3 * pc_y[k] * slg_674[k];

        t_944[k] = f_0 * skg_539[k]
                   + f_1 * slf0_449[k]
                   - f_2 * slf1_449[k]
                   + f_3 * pc_z[k] * slg_674[k];
    }
}

auto
compute_prim_slh_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t skh0, const size_t skg,
                                                   const size_t skh1, const size_t slf0,
                                                   const size_t slf1, const size_t slg,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_slh_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, skh0, skg,
                                                              skh1, slf0, slf1, slg, ncols,
                                                              gamma, p, q);

    compute_prim_slh_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, skh0, skg,
                                                              skh1, slf0, slf1, slg, ncols,
                                                              gamma, p, q);

    compute_prim_slh_three_center_electron_repulsion_0_piece2(buffer, target, pb, pc, skh0, skg,
                                                              skh1, slf0, slf1, slg, ncols,
                                                              gamma, p, q);

    compute_prim_slh_three_center_electron_repulsion_0_piece3(buffer, target, pb, pc, skh0, skg,
                                                              skh1, slf0, slf1, slg, ncols,
                                                              gamma, p, q);

    compute_prim_slh_three_center_electron_repulsion_0_piece4(buffer, target, pb, pc, skh0, skg,
                                                              skh1, slf0, slf1, slg, ncols,
                                                              gamma, p, q);

    compute_prim_slh_three_center_electron_repulsion_0_piece5(buffer, target, pb, pc, skh0, skg,
                                                              skh1, slf0, slf1, slg, ncols,
                                                              gamma, p, q);

    compute_prim_slh_three_center_electron_repulsion_0_piece6(buffer, target, pb, pc, skh0, skg,
                                                              skh1, slf0, slf1, slg, ncols,
                                                              gamma, p, q);

    compute_prim_slh_three_center_electron_repulsion_0_piece7(buffer, target, pb, pc, skh0, skg,
                                                              skh1, slf0, slf1, slg, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
