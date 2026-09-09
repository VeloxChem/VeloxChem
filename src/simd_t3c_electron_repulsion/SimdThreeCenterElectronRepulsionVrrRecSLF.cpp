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


#include "SimdThreeCenterElectronRepulsionVrrRecSLF.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_slf_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t skf0,
                                                          const size_t skd, const size_t skf1,
                                                          const size_t slp0, const size_t slp1,
                                                          const size_t sld, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 3.5 / q;
    const auto f_7 = 3.0 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 2.5 / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 2.0 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *skf0_0 = buffer.data(skf0 + 0);
    const auto *skf0_6 = buffer.data(skf0 + 6);
    const auto *skf0_9 = buffer.data(skf0 + 9);
    const auto *skf0_16 = buffer.data(skf0 + 16);
    const auto *skf0_20 = buffer.data(skf0 + 20);
    const auto *skf0_29 = buffer.data(skf0 + 29);
    const auto *skf0_30 = buffer.data(skf0 + 30);
    const auto *skf0_36 = buffer.data(skf0 + 36);
    const auto *skf0_50 = buffer.data(skf0 + 50);
    const auto *skf0_59 = buffer.data(skf0 + 59);
    const auto *skf0_60 = buffer.data(skf0 + 60);
    const auto *skf0_66 = buffer.data(skf0 + 66);

    const auto *skd_0 = buffer.data(skd + 0);
    const auto *skd_3 = buffer.data(skd + 3);
    const auto *skd_4 = buffer.data(skd + 4);
    const auto *skd_5 = buffer.data(skd + 5);
    const auto *skd_6 = buffer.data(skd + 6);
    const auto *skd_9 = buffer.data(skd + 9);
    const auto *skd_10 = buffer.data(skd + 10);
    const auto *skd_11 = buffer.data(skd + 11);
    const auto *skd_12 = buffer.data(skd + 12);
    const auto *skd_15 = buffer.data(skd + 15);
    const auto *skd_16 = buffer.data(skd + 16);
    const auto *skd_17 = buffer.data(skd + 17);
    const auto *skd_18 = buffer.data(skd + 18);
    const auto *skd_21 = buffer.data(skd + 21);
    const auto *skd_22 = buffer.data(skd + 22);
    const auto *skd_23 = buffer.data(skd + 23);
    const auto *skd_24 = buffer.data(skd + 24);
    const auto *skd_27 = buffer.data(skd + 27);
    const auto *skd_28 = buffer.data(skd + 28);
    const auto *skd_29 = buffer.data(skd + 29);
    const auto *skd_30 = buffer.data(skd + 30);
    const auto *skd_33 = buffer.data(skd + 33);
    const auto *skd_34 = buffer.data(skd + 34);
    const auto *skd_35 = buffer.data(skd + 35);
    const auto *skd_36 = buffer.data(skd + 36);
    const auto *skd_39 = buffer.data(skd + 39);
    const auto *skd_40 = buffer.data(skd + 40);
    const auto *skd_41 = buffer.data(skd + 41);
    const auto *skd_42 = buffer.data(skd + 42);
    const auto *skd_45 = buffer.data(skd + 45);
    const auto *skd_46 = buffer.data(skd + 46);
    const auto *skd_47 = buffer.data(skd + 47);
    const auto *skd_48 = buffer.data(skd + 48);
    const auto *skd_51 = buffer.data(skd + 51);
    const auto *skd_52 = buffer.data(skd + 52);
    const auto *skd_53 = buffer.data(skd + 53);
    const auto *skd_54 = buffer.data(skd + 54);
    const auto *skd_57 = buffer.data(skd + 57);
    const auto *skd_58 = buffer.data(skd + 58);
    const auto *skd_59 = buffer.data(skd + 59);
    const auto *skd_60 = buffer.data(skd + 60);
    const auto *skd_63 = buffer.data(skd + 63);
    const auto *skd_64 = buffer.data(skd + 64);
    const auto *skd_65 = buffer.data(skd + 65);
    const auto *skd_69 = buffer.data(skd + 69);
    const auto *skd_70 = buffer.data(skd + 70);
    const auto *skd_71 = buffer.data(skd + 71);
    const auto *skd_72 = buffer.data(skd + 72);
    const auto *skd_75 = buffer.data(skd + 75);
    const auto *skd_76 = buffer.data(skd + 76);
    const auto *skd_77 = buffer.data(skd + 77);

    const auto *skf1_0 = buffer.data(skf1 + 0);
    const auto *skf1_6 = buffer.data(skf1 + 6);
    const auto *skf1_9 = buffer.data(skf1 + 9);
    const auto *skf1_16 = buffer.data(skf1 + 16);
    const auto *skf1_20 = buffer.data(skf1 + 20);
    const auto *skf1_29 = buffer.data(skf1 + 29);
    const auto *skf1_30 = buffer.data(skf1 + 30);
    const auto *skf1_36 = buffer.data(skf1 + 36);
    const auto *skf1_50 = buffer.data(skf1 + 50);
    const auto *skf1_59 = buffer.data(skf1 + 59);
    const auto *skf1_60 = buffer.data(skf1 + 60);
    const auto *skf1_66 = buffer.data(skf1 + 66);

    const auto *slp0_0 = buffer.data(slp0 + 0);
    const auto *slp0_1 = buffer.data(slp0 + 1);
    const auto *slp0_2 = buffer.data(slp0 + 2);
    const auto *slp0_4 = buffer.data(slp0 + 4);
    const auto *slp0_8 = buffer.data(slp0 + 8);
    const auto *slp0_9 = buffer.data(slp0 + 9);
    const auto *slp0_10 = buffer.data(slp0 + 10);
    const auto *slp0_11 = buffer.data(slp0 + 11);
    const auto *slp0_15 = buffer.data(slp0 + 15);
    const auto *slp0_16 = buffer.data(slp0 + 16);
    const auto *slp0_17 = buffer.data(slp0 + 17);
    const auto *slp0_18 = buffer.data(slp0 + 18);
    const auto *slp0_19 = buffer.data(slp0 + 19);
    const auto *slp0_20 = buffer.data(slp0 + 20);
    const auto *slp0_23 = buffer.data(slp0 + 23);
    const auto *slp0_25 = buffer.data(slp0 + 25);
    const auto *slp0_27 = buffer.data(slp0 + 27);
    const auto *slp0_28 = buffer.data(slp0 + 28);
    const auto *slp0_29 = buffer.data(slp0 + 29);
    const auto *slp0_30 = buffer.data(slp0 + 30);
    const auto *slp0_31 = buffer.data(slp0 + 31);
    const auto *slp0_32 = buffer.data(slp0 + 32);
    const auto *slp0_35 = buffer.data(slp0 + 35);
    const auto *slp0_36 = buffer.data(slp0 + 36);
    const auto *slp0_37 = buffer.data(slp0 + 37);
    const auto *slp0_38 = buffer.data(slp0 + 38);

    const auto *slp1_0 = buffer.data(slp1 + 0);
    const auto *slp1_1 = buffer.data(slp1 + 1);
    const auto *slp1_2 = buffer.data(slp1 + 2);
    const auto *slp1_4 = buffer.data(slp1 + 4);
    const auto *slp1_8 = buffer.data(slp1 + 8);
    const auto *slp1_9 = buffer.data(slp1 + 9);
    const auto *slp1_10 = buffer.data(slp1 + 10);
    const auto *slp1_11 = buffer.data(slp1 + 11);
    const auto *slp1_15 = buffer.data(slp1 + 15);
    const auto *slp1_16 = buffer.data(slp1 + 16);
    const auto *slp1_17 = buffer.data(slp1 + 17);
    const auto *slp1_18 = buffer.data(slp1 + 18);
    const auto *slp1_19 = buffer.data(slp1 + 19);
    const auto *slp1_20 = buffer.data(slp1 + 20);
    const auto *slp1_23 = buffer.data(slp1 + 23);
    const auto *slp1_25 = buffer.data(slp1 + 25);
    const auto *slp1_27 = buffer.data(slp1 + 27);
    const auto *slp1_28 = buffer.data(slp1 + 28);
    const auto *slp1_29 = buffer.data(slp1 + 29);
    const auto *slp1_30 = buffer.data(slp1 + 30);
    const auto *slp1_31 = buffer.data(slp1 + 31);
    const auto *slp1_32 = buffer.data(slp1 + 32);
    const auto *slp1_35 = buffer.data(slp1 + 35);
    const auto *slp1_36 = buffer.data(slp1 + 36);
    const auto *slp1_37 = buffer.data(slp1 + 37);
    const auto *slp1_38 = buffer.data(slp1 + 38);

    const auto *sld_0 = buffer.data(sld + 0);
    const auto *sld_3 = buffer.data(sld + 3);
    const auto *sld_4 = buffer.data(sld + 4);
    const auto *sld_5 = buffer.data(sld + 5);
    const auto *sld_6 = buffer.data(sld + 6);
    const auto *sld_9 = buffer.data(sld + 9);
    const auto *sld_10 = buffer.data(sld + 10);
    const auto *sld_11 = buffer.data(sld + 11);
    const auto *sld_12 = buffer.data(sld + 12);
    const auto *sld_15 = buffer.data(sld + 15);
    const auto *sld_16 = buffer.data(sld + 16);
    const auto *sld_17 = buffer.data(sld + 17);
    const auto *sld_18 = buffer.data(sld + 18);
    const auto *sld_21 = buffer.data(sld + 21);
    const auto *sld_22 = buffer.data(sld + 22);
    const auto *sld_23 = buffer.data(sld + 23);
    const auto *sld_24 = buffer.data(sld + 24);
    const auto *sld_27 = buffer.data(sld + 27);
    const auto *sld_28 = buffer.data(sld + 28);
    const auto *sld_29 = buffer.data(sld + 29);
    const auto *sld_30 = buffer.data(sld + 30);
    const auto *sld_33 = buffer.data(sld + 33);
    const auto *sld_34 = buffer.data(sld + 34);
    const auto *sld_35 = buffer.data(sld + 35);
    const auto *sld_36 = buffer.data(sld + 36);
    const auto *sld_39 = buffer.data(sld + 39);
    const auto *sld_40 = buffer.data(sld + 40);
    const auto *sld_41 = buffer.data(sld + 41);
    const auto *sld_42 = buffer.data(sld + 42);
    const auto *sld_45 = buffer.data(sld + 45);
    const auto *sld_46 = buffer.data(sld + 46);
    const auto *sld_47 = buffer.data(sld + 47);
    const auto *sld_48 = buffer.data(sld + 48);
    const auto *sld_51 = buffer.data(sld + 51);
    const auto *sld_52 = buffer.data(sld + 52);
    const auto *sld_53 = buffer.data(sld + 53);
    const auto *sld_54 = buffer.data(sld + 54);
    const auto *sld_57 = buffer.data(sld + 57);
    const auto *sld_58 = buffer.data(sld + 58);
    const auto *sld_59 = buffer.data(sld + 59);
    const auto *sld_60 = buffer.data(sld + 60);
    const auto *sld_63 = buffer.data(sld + 63);
    const auto *sld_64 = buffer.data(sld + 64);
    const auto *sld_65 = buffer.data(sld + 65);
    const auto *sld_66 = buffer.data(sld + 66);
    const auto *sld_69 = buffer.data(sld + 69);
    const auto *sld_70 = buffer.data(sld + 70);
    const auto *sld_71 = buffer.data(sld + 71);
    const auto *sld_72 = buffer.data(sld + 72);
    const auto *sld_75 = buffer.data(sld + 75);
    const auto *sld_76 = buffer.data(sld + 76);
    const auto *sld_77 = buffer.data(sld + 77);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pc_x, pc_y, pc_z, skd_0, skd_3, skd_4, \
                         slp0_0, slp1_0, sld_0, sld_3, sld_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * skd_0[k]
                 + f_1 * slp0_0[k]
                 - f_2 * slp1_0[k]
                 + f_3 * pc_x[k] * sld_0[k];

        t_1[k] = f_3 * pc_y[k] * sld_0[k];

        t_2[k] = f_3 * pc_z[k] * sld_0[k];

        t_3[k] = f_0 * skd_3[k]
                 + f_3 * pc_x[k] * sld_3[k];

        t_4[k] = f_0 * skd_4[k]
                 + f_3 * pc_x[k] * sld_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pc_x, pc_y, pc_z, skd_5, slp0_1, slp0_2, \
                         slp1_1, slp1_2, sld_3, sld_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * skd_5[k]
                 + f_3 * pc_x[k] * sld_5[k];

        t_6[k] = f_1 * slp0_1[k]
                 - f_2 * slp1_1[k]
                 + f_3 * pc_y[k] * sld_3[k];

        t_7[k] = f_3 * pc_z[k] * sld_3[k];

        t_8[k] = f_3 * pc_y[k] * sld_5[k];

        t_9[k] = f_1 * slp0_2[k]
                 - f_2 * slp1_2[k]
                 + f_3 * pc_z[k] * sld_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pb_y, pc_x, pc_y, pc_z, skf0_0, skd_0, skd_9, \
                         skf1_0, sld_6, sld_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pb_y[k] * skf0_0[k]
                  - f_4 * pc_y[k] * skf1_0[k];

        t_11[k] = f_5 * skd_0[k]
                  + f_3 * pc_y[k] * sld_6[k];

        t_12[k] = f_3 * pc_z[k] * sld_6[k];

        t_13[k] = f_6 * skd_9[k]
                  + f_3 * pc_x[k] * sld_9[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pc_x, pc_y, pc_z, skd_3, skd_10, skd_11, \
                         slp0_4, slp1_4, sld_9, sld_10, sld_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_6 * skd_10[k]
                  + f_3 * pc_x[k] * sld_10[k];

        t_15[k] = f_6 * skd_11[k]
                  + f_3 * pc_x[k] * sld_11[k];

        t_16[k] = f_5 * skd_3[k]
                  + f_1 * slp0_4[k]
                  - f_2 * slp1_4[k]
                  + f_3 * pc_y[k] * sld_9[k];

        t_17[k] = f_3 * pc_z[k] * sld_9[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pb_y, pb_z, pc_y, pc_z, skf0_0, skf0_9, \
                         skd_5, skf1_0, skf1_9, sld_11, sld_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * skd_5[k]
                  + f_3 * pc_y[k] * sld_11[k];

        t_19[k] = pb_y[k] * skf0_9[k]
                  - f_4 * pc_y[k] * skf1_9[k];

        t_20[k] = pb_z[k] * skf0_0[k]
                  - f_4 * pc_z[k] * skf1_0[k];

        t_21[k] = f_3 * pc_y[k] * sld_12[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pc_x, pc_z, skd_0, skd_15, skd_16, skd_17, \
                         sld_12, sld_15, sld_16, sld_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_5 * skd_0[k]
                  + f_3 * pc_z[k] * sld_12[k];

        t_23[k] = f_6 * skd_15[k]
                  + f_3 * pc_x[k] * sld_15[k];

        t_24[k] = f_6 * skd_16[k]
                  + f_3 * pc_x[k] * sld_16[k];

        t_25[k] = f_6 * skd_17[k]
                  + f_3 * pc_x[k] * sld_17[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pb_z, pc_y, pc_z, skf0_6, skd_3, skd_5, \
                         skf1_6, slp0_8, slp1_8, sld_15, sld_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pb_z[k] * skf0_6[k]
                  - f_4 * pc_z[k] * skf1_6[k];

        t_27[k] = f_5 * skd_3[k]
                  + f_3 * pc_z[k] * sld_15[k];

        t_28[k] = f_3 * pc_y[k] * sld_17[k];

        t_29[k] = f_5 * skd_5[k]
                  + f_1 * slp0_8[k]
                  - f_2 * slp1_8[k]
                  + f_3 * pc_z[k] * sld_17[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pc_x, pc_y, pc_z, skd_6, skd_18, skd_21, \
                         slp0_9, slp1_9, sld_18, sld_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_7 * skd_18[k]
                  + f_1 * slp0_9[k]
                  - f_2 * slp1_9[k]
                  + f_3 * pc_x[k] * sld_18[k];

        t_31[k] = f_8 * skd_6[k]
                  + f_3 * pc_y[k] * sld_18[k];

        t_32[k] = f_3 * pc_z[k] * sld_18[k];

        t_33[k] = f_7 * skd_21[k]
                  + f_3 * pc_x[k] * sld_21[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pc_x, pc_y, pc_z, skd_9, skd_22, skd_23, \
                         slp0_10, slp1_10, sld_21, sld_22, sld_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_7 * skd_22[k]
                  + f_3 * pc_x[k] * sld_22[k];

        t_35[k] = f_7 * skd_23[k]
                  + f_3 * pc_x[k] * sld_23[k];

        t_36[k] = f_8 * skd_9[k]
                  + f_1 * slp0_10[k]
                  - f_2 * slp1_10[k]
                  + f_3 * pc_y[k] * sld_21[k];

        t_37[k] = f_3 * pc_z[k] * sld_21[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pb_y, pc_y, pc_z, skf0_20, skd_11, skd_12, \
                         skf1_20, slp0_11, slp1_11, sld_23, sld_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_8 * skd_11[k]
                  + f_3 * pc_y[k] * sld_23[k];

        t_39[k] = f_1 * slp0_11[k]
                  - f_2 * slp1_11[k]
                  + f_3 * pc_z[k] * sld_23[k];

        t_40[k] = pb_y[k] * skf0_20[k]
                  - f_4 * pc_y[k] * skf1_20[k];

        t_41[k] = f_5 * skd_12[k]
                  + f_3 * pc_y[k] * sld_24[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pc_x, pc_z, skd_6, skd_27, skd_28, skd_29, \
                         sld_24, sld_27, sld_28, sld_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_5 * skd_6[k]
                  + f_3 * pc_z[k] * sld_24[k];

        t_43[k] = f_7 * skd_27[k]
                  + f_3 * pc_x[k] * sld_27[k];

        t_44[k] = f_7 * skd_28[k]
                  + f_3 * pc_x[k] * sld_28[k];

        t_45[k] = f_7 * skd_29[k]
                  + f_3 * pc_x[k] * sld_29[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pb_y, pb_z, pc_y, pc_z, skf0_16, skf0_29, \
                         skd_9, skd_17, skf1_16, skf1_29, sld_27, \
                         sld_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = pb_z[k] * skf0_16[k]
                  - f_4 * pc_z[k] * skf1_16[k];

        t_47[k] = f_5 * skd_9[k]
                  + f_3 * pc_z[k] * sld_27[k];

        t_48[k] = f_5 * skd_17[k]
                  + f_3 * pc_y[k] * sld_29[k];

        t_49[k] = pb_y[k] * skf0_29[k]
                  - f_4 * pc_y[k] * skf1_29[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pc_x, pc_y, pc_z, skd_12, skd_30, skd_33, \
                         slp0_15, slp1_15, sld_30, sld_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_7 * skd_30[k]
                  + f_1 * slp0_15[k]
                  - f_2 * slp1_15[k]
                  + f_3 * pc_x[k] * sld_30[k];

        t_51[k] = f_3 * pc_y[k] * sld_30[k];

        t_52[k] = f_8 * skd_12[k]
                  + f_3 * pc_z[k] * sld_30[k];

        t_53[k] = f_7 * skd_33[k]
                  + f_3 * pc_x[k] * sld_33[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, pc_x, pc_y, pc_z, skd_15, skd_34, \
                         skd_35, slp0_16, slp1_16, sld_33, sld_34, \
                         sld_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_7 * skd_34[k]
                  + f_3 * pc_x[k] * sld_34[k];

        t_55[k] = f_7 * skd_35[k]
                  + f_3 * pc_x[k] * sld_35[k];

        t_56[k] = f_1 * slp0_16[k]
                  - f_2 * slp1_16[k]
                  + f_3 * pc_y[k] * sld_33[k];

        t_57[k] = f_8 * skd_15[k]
                  + f_3 * pc_z[k] * sld_33[k];

        t_58[k] = f_3 * pc_y[k] * sld_35[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pc_x, pc_y, pc_z, skd_17, skd_18, skd_36, \
                         slp0_17, slp0_18, slp1_17, slp1_18, sld_35, \
                         sld_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_8 * skd_17[k]
                  + f_1 * slp0_17[k]
                  - f_2 * slp1_17[k]
                  + f_3 * pc_z[k] * sld_35[k];

        t_60[k] = f_9 * skd_36[k]
                  + f_1 * slp0_18[k]
                  - f_2 * slp1_18[k]
                  + f_3 * pc_x[k] * sld_36[k];

        t_61[k] = f_10 * skd_18[k]
                  + f_3 * pc_y[k] * sld_36[k];

        t_62[k] = f_3 * pc_z[k] * sld_36[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pc_x, pc_y, skd_21, skd_39, skd_40, skd_41, \
                         slp0_19, slp1_19, sld_39, sld_40, sld_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_9 * skd_39[k]
                  + f_3 * pc_x[k] * sld_39[k];

        t_64[k] = f_9 * skd_40[k]
                  + f_3 * pc_x[k] * sld_40[k];

        t_65[k] = f_9 * skd_41[k]
                  + f_3 * pc_x[k] * sld_41[k];

        t_66[k] = f_10 * skd_21[k]
                  + f_1 * slp0_19[k]
                  - f_2 * slp1_19[k]
                  + f_3 * pc_y[k] * sld_39[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, pb_z, pc_y, pc_z, skf0_30, skd_23, skf1_30, \
                         slp0_20, slp1_20, sld_39, sld_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_3 * pc_z[k] * sld_39[k];

        t_68[k] = f_10 * skd_23[k]
                  + f_3 * pc_y[k] * sld_41[k];

        t_69[k] = f_1 * slp0_20[k]
                  - f_2 * slp1_20[k]
                  + f_3 * pc_z[k] * sld_41[k];

        t_70[k] = pb_z[k] * skf0_30[k]
                  - f_4 * pc_z[k] * skf1_30[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pc_x, pc_y, pc_z, skd_18, skd_24, skd_45, \
                         skd_46, sld_42, sld_45, sld_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_8 * skd_24[k]
                  + f_3 * pc_y[k] * sld_42[k];

        t_72[k] = f_5 * skd_18[k]
                  + f_3 * pc_z[k] * sld_42[k];

        t_73[k] = f_9 * skd_45[k]
                  + f_3 * pc_x[k] * sld_45[k];

        t_74[k] = f_9 * skd_46[k]
                  + f_3 * pc_x[k] * sld_46[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pb_z, pc_x, pc_y, pc_z, skf0_36, skd_21, \
                         skd_29, skd_47, skf1_36, sld_45, sld_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_9 * skd_47[k]
                  + f_3 * pc_x[k] * sld_47[k];

        t_76[k] = pb_z[k] * skf0_36[k]
                  - f_4 * pc_z[k] * skf1_36[k];

        t_77[k] = f_5 * skd_21[k]
                  + f_3 * pc_z[k] * sld_45[k];

        t_78[k] = f_8 * skd_29[k]
                  + f_3 * pc_y[k] * sld_47[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pb_y, pc_y, pc_z, skf0_50, skd_23, skd_24, \
                         skd_30, skf1_50, slp0_23, slp1_23, sld_47, \
                         sld_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_5 * skd_23[k]
                  + f_1 * slp0_23[k]
                  - f_2 * slp1_23[k]
                  + f_3 * pc_z[k] * sld_47[k];

        t_80[k] = pb_y[k] * skf0_50[k]
                  - f_4 * pc_y[k] * skf1_50[k];

        t_81[k] = f_5 * skd_30[k]
                  + f_3 * pc_y[k] * sld_48[k];

        t_82[k] = f_8 * skd_24[k]
                  + f_3 * pc_z[k] * sld_48[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pc_x, pc_y, skd_33, skd_51, skd_52, skd_53, \
                         slp0_25, slp1_25, sld_51, sld_52, sld_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_9 * skd_51[k]
                  + f_3 * pc_x[k] * sld_51[k];

        t_84[k] = f_9 * skd_52[k]
                  + f_3 * pc_x[k] * sld_52[k];

        t_85[k] = f_9 * skd_53[k]
                  + f_3 * pc_x[k] * sld_53[k];

        t_86[k] = f_5 * skd_33[k]
                  + f_1 * slp0_25[k]
                  - f_2 * slp1_25[k]
                  + f_3 * pc_y[k] * sld_51[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, pb_y, pc_y, pc_z, skf0_59, skd_27, skd_35, skf1_59, \
                         sld_51, sld_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_8 * skd_27[k]
                  + f_3 * pc_z[k] * sld_51[k];

        t_88[k] = f_5 * skd_35[k]
                  + f_3 * pc_y[k] * sld_53[k];

        t_89[k] = pb_y[k] * skf0_59[k]
                  - f_4 * pc_y[k] * skf1_59[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, pc_x, pc_y, pc_z, skd_30, skd_54, skd_57, \
                         slp0_27, slp1_27, sld_54, sld_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_9 * skd_54[k]
                  + f_1 * slp0_27[k]
                  - f_2 * slp1_27[k]
                  + f_3 * pc_x[k] * sld_54[k];

        t_91[k] = f_3 * pc_y[k] * sld_54[k];

        t_92[k] = f_10 * skd_30[k]
                  + f_3 * pc_z[k] * sld_54[k];

        t_93[k] = f_9 * skd_57[k]
                  + f_3 * pc_x[k] * sld_57[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, pc_x, pc_y, pc_z, skd_33, skd_58, \
                         skd_59, slp0_28, slp1_28, sld_57, sld_58, \
                         sld_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_9 * skd_58[k]
                  + f_3 * pc_x[k] * sld_58[k];

        t_95[k] = f_9 * skd_59[k]
                  + f_3 * pc_x[k] * sld_59[k];

        t_96[k] = f_1 * slp0_28[k]
                  - f_2 * slp1_28[k]
                  + f_3 * pc_y[k] * sld_57[k];

        t_97[k] = f_10 * skd_33[k]
                  + f_3 * pc_z[k] * sld_57[k];

        t_98[k] = f_3 * pc_y[k] * sld_59[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, pc_x, pc_y, pc_z, skd_35, skd_36, skd_60, \
                         slp0_29, slp0_30, slp1_29, slp1_30, sld_59, \
                         sld_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_10 * skd_35[k]
                  + f_1 * slp0_29[k]
                  - f_2 * slp1_29[k]
                  + f_3 * pc_z[k] * sld_59[k];

        t_100[k] = f_11 * skd_60[k]
                   + f_1 * slp0_30[k]
                   - f_2 * slp1_30[k]
                   + f_3 * pc_x[k] * sld_60[k];

        t_101[k] = f_11 * skd_36[k]
                   + f_3 * pc_y[k] * sld_60[k];

        t_102[k] = f_3 * pc_z[k] * sld_60[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, pc_x, pc_y, skd_39, skd_63, skd_64, \
                         skd_65, slp0_31, slp1_31, sld_63, sld_64, \
                         sld_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_11 * skd_63[k]
                   + f_3 * pc_x[k] * sld_63[k];

        t_104[k] = f_11 * skd_64[k]
                   + f_3 * pc_x[k] * sld_64[k];

        t_105[k] = f_11 * skd_65[k]
                   + f_3 * pc_x[k] * sld_65[k];

        t_106[k] = f_11 * skd_39[k]
                   + f_1 * slp0_31[k]
                   - f_2 * slp1_31[k]
                   + f_3 * pc_y[k] * sld_63[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, pb_z, pc_y, pc_z, skf0_60, skd_41, \
                         skf1_60, slp0_32, slp1_32, sld_63, sld_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_3 * pc_z[k] * sld_63[k];

        t_108[k] = f_11 * skd_41[k]
                   + f_3 * pc_y[k] * sld_65[k];

        t_109[k] = f_1 * slp0_32[k]
                   - f_2 * slp1_32[k]
                   + f_3 * pc_z[k] * sld_65[k];

        t_110[k] = pb_z[k] * skf0_60[k]
                   - f_4 * pc_z[k] * skf1_60[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, pc_x, pc_y, pc_z, skd_36, skd_42, skd_69, \
                         skd_70, sld_66, sld_69, sld_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_10 * skd_42[k]
                   + f_3 * pc_y[k] * sld_66[k];

        t_112[k] = f_5 * skd_36[k]
                   + f_3 * pc_z[k] * sld_66[k];

        t_113[k] = f_11 * skd_69[k]
                   + f_3 * pc_x[k] * sld_69[k];

        t_114[k] = f_11 * skd_70[k]
                   + f_3 * pc_x[k] * sld_70[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, pb_z, pc_x, pc_y, pc_z, skf0_66, skd_39, \
                         skd_47, skd_71, skf1_66, sld_69, sld_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_11 * skd_71[k]
                   + f_3 * pc_x[k] * sld_71[k];

        t_116[k] = pb_z[k] * skf0_66[k]
                   - f_4 * pc_z[k] * skf1_66[k];

        t_117[k] = f_5 * skd_39[k]
                   + f_3 * pc_z[k] * sld_69[k];

        t_118[k] = f_10 * skd_47[k]
                   + f_3 * pc_y[k] * sld_71[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, pc_x, pc_y, pc_z, skd_41, skd_48, skd_72, \
                         slp0_35, slp0_36, slp1_35, slp1_36, sld_71, \
                         sld_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = f_5 * skd_41[k]
                   + f_1 * slp0_35[k]
                   - f_2 * slp1_35[k]
                   + f_3 * pc_z[k] * sld_71[k];

        t_120[k] = f_11 * skd_72[k]
                   + f_1 * slp0_36[k]
                   - f_2 * slp1_36[k]
                   + f_3 * pc_x[k] * sld_72[k];

        t_121[k] = f_8 * skd_48[k]
                   + f_3 * pc_y[k] * sld_72[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, pc_x, pc_z, skd_42, skd_75, skd_76, \
                         skd_77, sld_72, sld_75, sld_76, sld_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_8 * skd_42[k]
                   + f_3 * pc_z[k] * sld_72[k];

        t_123[k] = f_11 * skd_75[k]
                   + f_3 * pc_x[k] * sld_75[k];

        t_124[k] = f_11 * skd_76[k]
                   + f_3 * pc_x[k] * sld_76[k];

        t_125[k] = f_11 * skd_77[k]
                   + f_3 * pc_x[k] * sld_77[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, pc_y, pc_z, skd_45, skd_47, skd_51, \
                         skd_53, slp0_37, slp0_38, slp1_37, slp1_38, sld_75, \
                         sld_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_8 * skd_51[k]
                   + f_1 * slp0_37[k]
                   - f_2 * slp1_37[k]
                   + f_3 * pc_y[k] * sld_75[k];

        t_127[k] = f_8 * skd_45[k]
                   + f_3 * pc_z[k] * sld_75[k];

        t_128[k] = f_8 * skd_53[k]
                   + f_3 * pc_y[k] * sld_77[k];

        t_129[k] = f_8 * skd_47[k]
                   + f_1 * slp0_38[k]
                   - f_2 * slp1_38[k]
                   + f_3 * pc_z[k] * sld_77[k];
    }
}

static auto
compute_prim_slf_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t skf0,
                                                          const size_t skd, const size_t skf1,
                                                          const size_t slp0, const size_t slp1,
                                                          const size_t sld, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_7 = 3.0 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 2.5 / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 2.0 / q;

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
    auto *t_254 = buffer.data(target + 254);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *skf0_90 = buffer.data(skf0 + 90);
    const auto *skf0_99 = buffer.data(skf0 + 99);
    const auto *skf0_100 = buffer.data(skf0 + 100);
    const auto *skf0_106 = buffer.data(skf0 + 106);
    const auto *skf0_140 = buffer.data(skf0 + 140);
    const auto *skf0_149 = buffer.data(skf0 + 149);
    const auto *skf0_150 = buffer.data(skf0 + 150);
    const auto *skf0_156 = buffer.data(skf0 + 156);

    const auto *skd_48 = buffer.data(skd + 48);
    const auto *skd_51 = buffer.data(skd + 51);
    const auto *skd_54 = buffer.data(skd + 54);
    const auto *skd_57 = buffer.data(skd + 57);
    const auto *skd_59 = buffer.data(skd + 59);
    const auto *skd_60 = buffer.data(skd + 60);
    const auto *skd_63 = buffer.data(skd + 63);
    const auto *skd_65 = buffer.data(skd + 65);
    const auto *skd_66 = buffer.data(skd + 66);
    const auto *skd_69 = buffer.data(skd + 69);
    const auto *skd_71 = buffer.data(skd + 71);
    const auto *skd_72 = buffer.data(skd + 72);
    const auto *skd_75 = buffer.data(skd + 75);
    const auto *skd_77 = buffer.data(skd + 77);
    const auto *skd_78 = buffer.data(skd + 78);
    const auto *skd_81 = buffer.data(skd + 81);
    const auto *skd_82 = buffer.data(skd + 82);
    const auto *skd_83 = buffer.data(skd + 83);
    const auto *skd_84 = buffer.data(skd + 84);
    const auto *skd_87 = buffer.data(skd + 87);
    const auto *skd_88 = buffer.data(skd + 88);
    const auto *skd_89 = buffer.data(skd + 89);
    const auto *skd_90 = buffer.data(skd + 90);
    const auto *skd_93 = buffer.data(skd + 93);
    const auto *skd_94 = buffer.data(skd + 94);
    const auto *skd_95 = buffer.data(skd + 95);
    const auto *skd_96 = buffer.data(skd + 96);
    const auto *skd_99 = buffer.data(skd + 99);
    const auto *skd_100 = buffer.data(skd + 100);
    const auto *skd_101 = buffer.data(skd + 101);
    const auto *skd_102 = buffer.data(skd + 102);
    const auto *skd_105 = buffer.data(skd + 105);
    const auto *skd_106 = buffer.data(skd + 106);
    const auto *skd_107 = buffer.data(skd + 107);
    const auto *skd_108 = buffer.data(skd + 108);
    const auto *skd_111 = buffer.data(skd + 111);
    const auto *skd_112 = buffer.data(skd + 112);
    const auto *skd_113 = buffer.data(skd + 113);
    const auto *skd_114 = buffer.data(skd + 114);
    const auto *skd_117 = buffer.data(skd + 117);
    const auto *skd_118 = buffer.data(skd + 118);
    const auto *skd_119 = buffer.data(skd + 119);
    const auto *skd_120 = buffer.data(skd + 120);
    const auto *skd_123 = buffer.data(skd + 123);
    const auto *skd_124 = buffer.data(skd + 124);
    const auto *skd_125 = buffer.data(skd + 125);
    const auto *skd_126 = buffer.data(skd + 126);
    const auto *skd_129 = buffer.data(skd + 129);
    const auto *skd_130 = buffer.data(skd + 130);
    const auto *skd_131 = buffer.data(skd + 131);
    const auto *skd_135 = buffer.data(skd + 135);
    const auto *skd_136 = buffer.data(skd + 136);
    const auto *skd_137 = buffer.data(skd + 137);
    const auto *skd_138 = buffer.data(skd + 138);
    const auto *skd_141 = buffer.data(skd + 141);
    const auto *skd_142 = buffer.data(skd + 142);
    const auto *skd_143 = buffer.data(skd + 143);
    const auto *skd_144 = buffer.data(skd + 144);
    const auto *skd_147 = buffer.data(skd + 147);
    const auto *skd_148 = buffer.data(skd + 148);
    const auto *skd_149 = buffer.data(skd + 149);
    const auto *skd_150 = buffer.data(skd + 150);
    const auto *skd_153 = buffer.data(skd + 153);
    const auto *skd_154 = buffer.data(skd + 154);

    const auto *skf1_90 = buffer.data(skf1 + 90);
    const auto *skf1_99 = buffer.data(skf1 + 99);
    const auto *skf1_100 = buffer.data(skf1 + 100);
    const auto *skf1_106 = buffer.data(skf1 + 106);
    const auto *skf1_140 = buffer.data(skf1 + 140);
    const auto *skf1_149 = buffer.data(skf1 + 149);
    const auto *skf1_150 = buffer.data(skf1 + 150);
    const auto *skf1_156 = buffer.data(skf1 + 156);

    const auto *slp0_40 = buffer.data(slp0 + 40);
    const auto *slp0_42 = buffer.data(slp0 + 42);
    const auto *slp0_43 = buffer.data(slp0 + 43);
    const auto *slp0_44 = buffer.data(slp0 + 44);
    const auto *slp0_45 = buffer.data(slp0 + 45);
    const auto *slp0_46 = buffer.data(slp0 + 46);
    const auto *slp0_47 = buffer.data(slp0 + 47);
    const auto *slp0_50 = buffer.data(slp0 + 50);
    const auto *slp0_51 = buffer.data(slp0 + 51);
    const auto *slp0_52 = buffer.data(slp0 + 52);
    const auto *slp0_53 = buffer.data(slp0 + 53);
    const auto *slp0_54 = buffer.data(slp0 + 54);
    const auto *slp0_55 = buffer.data(slp0 + 55);
    const auto *slp0_56 = buffer.data(slp0 + 56);
    const auto *slp0_58 = buffer.data(slp0 + 58);
    const auto *slp0_60 = buffer.data(slp0 + 60);
    const auto *slp0_61 = buffer.data(slp0 + 61);
    const auto *slp0_62 = buffer.data(slp0 + 62);
    const auto *slp0_63 = buffer.data(slp0 + 63);
    const auto *slp0_64 = buffer.data(slp0 + 64);
    const auto *slp0_65 = buffer.data(slp0 + 65);
    const auto *slp0_68 = buffer.data(slp0 + 68);
    const auto *slp0_69 = buffer.data(slp0 + 69);
    const auto *slp0_70 = buffer.data(slp0 + 70);
    const auto *slp0_71 = buffer.data(slp0 + 71);
    const auto *slp0_72 = buffer.data(slp0 + 72);
    const auto *slp0_73 = buffer.data(slp0 + 73);
    const auto *slp0_74 = buffer.data(slp0 + 74);
    const auto *slp0_75 = buffer.data(slp0 + 75);

    const auto *slp1_40 = buffer.data(slp1 + 40);
    const auto *slp1_42 = buffer.data(slp1 + 42);
    const auto *slp1_43 = buffer.data(slp1 + 43);
    const auto *slp1_44 = buffer.data(slp1 + 44);
    const auto *slp1_45 = buffer.data(slp1 + 45);
    const auto *slp1_46 = buffer.data(slp1 + 46);
    const auto *slp1_47 = buffer.data(slp1 + 47);
    const auto *slp1_50 = buffer.data(slp1 + 50);
    const auto *slp1_51 = buffer.data(slp1 + 51);
    const auto *slp1_52 = buffer.data(slp1 + 52);
    const auto *slp1_53 = buffer.data(slp1 + 53);
    const auto *slp1_54 = buffer.data(slp1 + 54);
    const auto *slp1_55 = buffer.data(slp1 + 55);
    const auto *slp1_56 = buffer.data(slp1 + 56);
    const auto *slp1_58 = buffer.data(slp1 + 58);
    const auto *slp1_60 = buffer.data(slp1 + 60);
    const auto *slp1_61 = buffer.data(slp1 + 61);
    const auto *slp1_62 = buffer.data(slp1 + 62);
    const auto *slp1_63 = buffer.data(slp1 + 63);
    const auto *slp1_64 = buffer.data(slp1 + 64);
    const auto *slp1_65 = buffer.data(slp1 + 65);
    const auto *slp1_68 = buffer.data(slp1 + 68);
    const auto *slp1_69 = buffer.data(slp1 + 69);
    const auto *slp1_70 = buffer.data(slp1 + 70);
    const auto *slp1_71 = buffer.data(slp1 + 71);
    const auto *slp1_72 = buffer.data(slp1 + 72);
    const auto *slp1_73 = buffer.data(slp1 + 73);
    const auto *slp1_74 = buffer.data(slp1 + 74);
    const auto *slp1_75 = buffer.data(slp1 + 75);

    const auto *sld_78 = buffer.data(sld + 78);
    const auto *sld_81 = buffer.data(sld + 81);
    const auto *sld_82 = buffer.data(sld + 82);
    const auto *sld_83 = buffer.data(sld + 83);
    const auto *sld_84 = buffer.data(sld + 84);
    const auto *sld_87 = buffer.data(sld + 87);
    const auto *sld_88 = buffer.data(sld + 88);
    const auto *sld_89 = buffer.data(sld + 89);
    const auto *sld_90 = buffer.data(sld + 90);
    const auto *sld_93 = buffer.data(sld + 93);
    const auto *sld_94 = buffer.data(sld + 94);
    const auto *sld_95 = buffer.data(sld + 95);
    const auto *sld_96 = buffer.data(sld + 96);
    const auto *sld_99 = buffer.data(sld + 99);
    const auto *sld_100 = buffer.data(sld + 100);
    const auto *sld_101 = buffer.data(sld + 101);
    const auto *sld_102 = buffer.data(sld + 102);
    const auto *sld_105 = buffer.data(sld + 105);
    const auto *sld_106 = buffer.data(sld + 106);
    const auto *sld_107 = buffer.data(sld + 107);
    const auto *sld_108 = buffer.data(sld + 108);
    const auto *sld_111 = buffer.data(sld + 111);
    const auto *sld_112 = buffer.data(sld + 112);
    const auto *sld_113 = buffer.data(sld + 113);
    const auto *sld_114 = buffer.data(sld + 114);
    const auto *sld_117 = buffer.data(sld + 117);
    const auto *sld_118 = buffer.data(sld + 118);
    const auto *sld_119 = buffer.data(sld + 119);
    const auto *sld_120 = buffer.data(sld + 120);
    const auto *sld_123 = buffer.data(sld + 123);
    const auto *sld_124 = buffer.data(sld + 124);
    const auto *sld_125 = buffer.data(sld + 125);
    const auto *sld_126 = buffer.data(sld + 126);
    const auto *sld_129 = buffer.data(sld + 129);
    const auto *sld_130 = buffer.data(sld + 130);
    const auto *sld_131 = buffer.data(sld + 131);
    const auto *sld_132 = buffer.data(sld + 132);
    const auto *sld_135 = buffer.data(sld + 135);
    const auto *sld_136 = buffer.data(sld + 136);
    const auto *sld_137 = buffer.data(sld + 137);
    const auto *sld_138 = buffer.data(sld + 138);
    const auto *sld_141 = buffer.data(sld + 141);
    const auto *sld_142 = buffer.data(sld + 142);
    const auto *sld_143 = buffer.data(sld + 143);
    const auto *sld_144 = buffer.data(sld + 144);
    const auto *sld_147 = buffer.data(sld + 147);
    const auto *sld_148 = buffer.data(sld + 148);
    const auto *sld_149 = buffer.data(sld + 149);
    const auto *sld_150 = buffer.data(sld + 150);
    const auto *sld_153 = buffer.data(sld + 153);
    const auto *sld_154 = buffer.data(sld + 154);

#pragma omp simd aligned(t_130, t_131, t_132, t_133, pb_y, pc_x, pc_y, pc_z, skf0_90, skd_48, \
                         skd_54, skd_81, skf1_90, sld_78, sld_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = pb_y[k] * skf0_90[k]
                   - f_4 * pc_y[k] * skf1_90[k];

        t_131[k] = f_5 * skd_54[k]
                   + f_3 * pc_y[k] * sld_78[k];

        t_132[k] = f_10 * skd_48[k]
                   + f_3 * pc_z[k] * sld_78[k];

        t_133[k] = f_11 * skd_81[k]
                   + f_3 * pc_x[k] * sld_81[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, pc_x, pc_y, pc_z, skd_51, skd_57, skd_82, \
                         skd_83, slp0_40, slp1_40, sld_81, sld_82, \
                         sld_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_11 * skd_82[k]
                   + f_3 * pc_x[k] * sld_82[k];

        t_135[k] = f_11 * skd_83[k]
                   + f_3 * pc_x[k] * sld_83[k];

        t_136[k] = f_5 * skd_57[k]
                   + f_1 * slp0_40[k]
                   - f_2 * slp1_40[k]
                   + f_3 * pc_y[k] * sld_81[k];

        t_137[k] = f_10 * skd_51[k]
                   + f_3 * pc_z[k] * sld_81[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, pb_y, pc_x, pc_y, skf0_99, skd_59, \
                         skd_84, skf1_99, slp0_42, slp1_42, sld_83, \
                         sld_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_5 * skd_59[k]
                   + f_3 * pc_y[k] * sld_83[k];

        t_139[k] = pb_y[k] * skf0_99[k]
                   - f_4 * pc_y[k] * skf1_99[k];

        t_140[k] = f_11 * skd_84[k]
                   + f_1 * slp0_42[k]
                   - f_2 * slp1_42[k]
                   + f_3 * pc_x[k] * sld_84[k];

        t_141[k] = f_3 * pc_y[k] * sld_84[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, pc_x, pc_z, skd_54, skd_87, skd_88, \
                         skd_89, sld_84, sld_87, sld_88, sld_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_11 * skd_54[k]
                   + f_3 * pc_z[k] * sld_84[k];

        t_143[k] = f_11 * skd_87[k]
                   + f_3 * pc_x[k] * sld_87[k];

        t_144[k] = f_11 * skd_88[k]
                   + f_3 * pc_x[k] * sld_88[k];

        t_145[k] = f_11 * skd_89[k]
                   + f_3 * pc_x[k] * sld_89[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pc_y, pc_z, skd_57, skd_59, slp0_43, \
                         slp0_44, slp1_43, slp1_44, sld_87, sld_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_1 * slp0_43[k]
                   - f_2 * slp1_43[k]
                   + f_3 * pc_y[k] * sld_87[k];

        t_147[k] = f_11 * skd_57[k]
                   + f_3 * pc_z[k] * sld_87[k];

        t_148[k] = f_3 * pc_y[k] * sld_89[k];

        t_149[k] = f_11 * skd_59[k]
                   + f_1 * slp0_44[k]
                   - f_2 * slp1_44[k]
                   + f_3 * pc_z[k] * sld_89[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, pc_x, pc_y, pc_z, skd_60, skd_90, skd_93, \
                         slp0_45, slp1_45, sld_90, sld_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_10 * skd_90[k]
                   + f_1 * slp0_45[k]
                   - f_2 * slp1_45[k]
                   + f_3 * pc_x[k] * sld_90[k];

        t_151[k] = f_9 * skd_60[k]
                   + f_3 * pc_y[k] * sld_90[k];

        t_152[k] = f_3 * pc_z[k] * sld_90[k];

        t_153[k] = f_10 * skd_93[k]
                   + f_3 * pc_x[k] * sld_93[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, pc_x, pc_y, pc_z, skd_63, skd_94, skd_95, \
                         slp0_46, slp1_46, sld_93, sld_94, sld_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = f_10 * skd_94[k]
                   + f_3 * pc_x[k] * sld_94[k];

        t_155[k] = f_10 * skd_95[k]
                   + f_3 * pc_x[k] * sld_95[k];

        t_156[k] = f_9 * skd_63[k]
                   + f_1 * slp0_46[k]
                   - f_2 * slp1_46[k]
                   + f_3 * pc_y[k] * sld_93[k];

        t_157[k] = f_3 * pc_z[k] * sld_93[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, pb_z, pc_y, pc_z, skf0_100, skd_65, \
                         skd_66, skf1_100, slp0_47, slp1_47, sld_95, \
                         sld_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_9 * skd_65[k]
                   + f_3 * pc_y[k] * sld_95[k];

        t_159[k] = f_1 * slp0_47[k]
                   - f_2 * slp1_47[k]
                   + f_3 * pc_z[k] * sld_95[k];

        t_160[k] = pb_z[k] * skf0_100[k]
                   - f_4 * pc_z[k] * skf1_100[k];

        t_161[k] = f_11 * skd_66[k]
                   + f_3 * pc_y[k] * sld_96[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, pc_x, pc_z, skd_60, skd_99, skd_100, \
                         skd_101, sld_96, sld_99, sld_100, sld_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_5 * skd_60[k]
                   + f_3 * pc_z[k] * sld_96[k];

        t_163[k] = f_10 * skd_99[k]
                   + f_3 * pc_x[k] * sld_99[k];

        t_164[k] = f_10 * skd_100[k]
                   + f_3 * pc_x[k] * sld_100[k];

        t_165[k] = f_10 * skd_101[k]
                   + f_3 * pc_x[k] * sld_101[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, pb_z, pc_y, pc_z, skf0_106, skd_63, \
                         skd_65, skd_71, skf1_106, slp0_50, slp1_50, sld_99, \
                         sld_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = pb_z[k] * skf0_106[k]
                   - f_4 * pc_z[k] * skf1_106[k];

        t_167[k] = f_5 * skd_63[k]
                   + f_3 * pc_z[k] * sld_99[k];

        t_168[k] = f_11 * skd_71[k]
                   + f_3 * pc_y[k] * sld_101[k];

        t_169[k] = f_5 * skd_65[k]
                   + f_1 * slp0_50[k]
                   - f_2 * slp1_50[k]
                   + f_3 * pc_z[k] * sld_101[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, pc_x, pc_y, pc_z, skd_66, skd_72, \
                         skd_102, skd_105, slp0_51, slp1_51, sld_102, \
                         sld_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_10 * skd_102[k]
                   + f_1 * slp0_51[k]
                   - f_2 * slp1_51[k]
                   + f_3 * pc_x[k] * sld_102[k];

        t_171[k] = f_10 * skd_72[k]
                   + f_3 * pc_y[k] * sld_102[k];

        t_172[k] = f_8 * skd_66[k]
                   + f_3 * pc_z[k] * sld_102[k];

        t_173[k] = f_10 * skd_105[k]
                   + f_3 * pc_x[k] * sld_105[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, pc_x, pc_y, pc_z, skd_69, skd_75, \
                         skd_106, skd_107, slp0_52, slp1_52, sld_105, sld_106, \
                         sld_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_10 * skd_106[k]
                   + f_3 * pc_x[k] * sld_106[k];

        t_175[k] = f_10 * skd_107[k]
                   + f_3 * pc_x[k] * sld_107[k];

        t_176[k] = f_10 * skd_75[k]
                   + f_1 * slp0_52[k]
                   - f_2 * slp1_52[k]
                   + f_3 * pc_y[k] * sld_105[k];

        t_177[k] = f_8 * skd_69[k]
                   + f_3 * pc_z[k] * sld_105[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, pc_x, pc_y, pc_z, skd_71, skd_77, skd_108, \
                         slp0_53, slp0_54, slp1_53, slp1_54, sld_107, \
                         sld_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_10 * skd_77[k]
                   + f_3 * pc_y[k] * sld_107[k];

        t_179[k] = f_8 * skd_71[k]
                   + f_1 * slp0_53[k]
                   - f_2 * slp1_53[k]
                   + f_3 * pc_z[k] * sld_107[k];

        t_180[k] = f_10 * skd_108[k]
                   + f_1 * slp0_54[k]
                   - f_2 * slp1_54[k]
                   + f_3 * pc_x[k] * sld_108[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pc_x, pc_y, pc_z, skd_72, skd_78, \
                         skd_111, skd_112, sld_108, sld_111, sld_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_8 * skd_78[k]
                   + f_3 * pc_y[k] * sld_108[k];

        t_182[k] = f_10 * skd_72[k]
                   + f_3 * pc_z[k] * sld_108[k];

        t_183[k] = f_10 * skd_111[k]
                   + f_3 * pc_x[k] * sld_111[k];

        t_184[k] = f_10 * skd_112[k]
                   + f_3 * pc_x[k] * sld_112[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pc_x, pc_y, pc_z, skd_75, skd_81, skd_83, \
                         skd_113, slp0_55, slp1_55, sld_111, sld_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_10 * skd_113[k]
                   + f_3 * pc_x[k] * sld_113[k];

        t_186[k] = f_8 * skd_81[k]
                   + f_1 * slp0_55[k]
                   - f_2 * slp1_55[k]
                   + f_3 * pc_y[k] * sld_111[k];

        t_187[k] = f_10 * skd_75[k]
                   + f_3 * pc_z[k] * sld_111[k];

        t_188[k] = f_8 * skd_83[k]
                   + f_3 * pc_y[k] * sld_113[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, pb_y, pc_y, pc_z, skf0_140, skd_77, \
                         skd_78, skd_84, skf1_140, slp0_56, slp1_56, sld_113, \
                         sld_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_10 * skd_77[k]
                   + f_1 * slp0_56[k]
                   - f_2 * slp1_56[k]
                   + f_3 * pc_z[k] * sld_113[k];

        t_190[k] = pb_y[k] * skf0_140[k]
                   - f_4 * pc_y[k] * skf1_140[k];

        t_191[k] = f_5 * skd_84[k]
                   + f_3 * pc_y[k] * sld_114[k];

        t_192[k] = f_11 * skd_78[k]
                   + f_3 * pc_z[k] * sld_114[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, pc_x, pc_y, skd_87, skd_117, skd_118, \
                         skd_119, slp0_58, slp1_58, sld_117, sld_118, \
                         sld_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_10 * skd_117[k]
                   + f_3 * pc_x[k] * sld_117[k];

        t_194[k] = f_10 * skd_118[k]
                   + f_3 * pc_x[k] * sld_118[k];

        t_195[k] = f_10 * skd_119[k]
                   + f_3 * pc_x[k] * sld_119[k];

        t_196[k] = f_5 * skd_87[k]
                   + f_1 * slp0_58[k]
                   - f_2 * slp1_58[k]
                   + f_3 * pc_y[k] * sld_117[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, pb_y, pc_y, pc_z, skf0_149, skd_81, skd_89, \
                         skf1_149, sld_117, sld_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_11 * skd_81[k]
                   + f_3 * pc_z[k] * sld_117[k];

        t_198[k] = f_5 * skd_89[k]
                   + f_3 * pc_y[k] * sld_119[k];

        t_199[k] = pb_y[k] * skf0_149[k]
                   - f_4 * pc_y[k] * skf1_149[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, pc_x, pc_y, pc_z, skd_84, skd_120, \
                         skd_123, slp0_60, slp1_60, sld_120, sld_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_10 * skd_120[k]
                   + f_1 * slp0_60[k]
                   - f_2 * slp1_60[k]
                   + f_3 * pc_x[k] * sld_120[k];

        t_201[k] = f_3 * pc_y[k] * sld_120[k];

        t_202[k] = f_9 * skd_84[k]
                   + f_3 * pc_z[k] * sld_120[k];

        t_203[k] = f_10 * skd_123[k]
                   + f_3 * pc_x[k] * sld_123[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, t_208, pc_x, pc_y, pc_z, skd_87, skd_124, \
                         skd_125, slp0_61, slp1_61, sld_123, sld_124, \
                         sld_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_10 * skd_124[k]
                   + f_3 * pc_x[k] * sld_124[k];

        t_205[k] = f_10 * skd_125[k]
                   + f_3 * pc_x[k] * sld_125[k];

        t_206[k] = f_1 * slp0_61[k]
                   - f_2 * slp1_61[k]
                   + f_3 * pc_y[k] * sld_123[k];

        t_207[k] = f_9 * skd_87[k]
                   + f_3 * pc_z[k] * sld_123[k];

        t_208[k] = f_3 * pc_y[k] * sld_125[k];
    }

#pragma omp simd aligned(t_209, t_210, t_211, t_212, pc_x, pc_y, pc_z, skd_89, skd_90, \
                         skd_126, slp0_62, slp0_63, slp1_62, slp1_63, sld_125, \
                         sld_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_209[k] = f_9 * skd_89[k]
                   + f_1 * slp0_62[k]
                   - f_2 * slp1_62[k]
                   + f_3 * pc_z[k] * sld_125[k];

        t_210[k] = f_8 * skd_126[k]
                   + f_1 * slp0_63[k]
                   - f_2 * slp1_63[k]
                   + f_3 * pc_x[k] * sld_126[k];

        t_211[k] = f_7 * skd_90[k]
                   + f_3 * pc_y[k] * sld_126[k];

        t_212[k] = f_3 * pc_z[k] * sld_126[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, pc_x, pc_y, skd_93, skd_129, skd_130, \
                         skd_131, slp0_64, slp1_64, sld_129, sld_130, \
                         sld_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_8 * skd_129[k]
                   + f_3 * pc_x[k] * sld_129[k];

        t_214[k] = f_8 * skd_130[k]
                   + f_3 * pc_x[k] * sld_130[k];

        t_215[k] = f_8 * skd_131[k]
                   + f_3 * pc_x[k] * sld_131[k];

        t_216[k] = f_7 * skd_93[k]
                   + f_1 * slp0_64[k]
                   - f_2 * slp1_64[k]
                   + f_3 * pc_y[k] * sld_129[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, pb_z, pc_y, pc_z, skf0_150, skd_95, \
                         skf1_150, slp0_65, slp1_65, sld_129, sld_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = f_3 * pc_z[k] * sld_129[k];

        t_218[k] = f_7 * skd_95[k]
                   + f_3 * pc_y[k] * sld_131[k];

        t_219[k] = f_1 * slp0_65[k]
                   - f_2 * slp1_65[k]
                   + f_3 * pc_z[k] * sld_131[k];

        t_220[k] = pb_z[k] * skf0_150[k]
                   - f_4 * pc_z[k] * skf1_150[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, pc_x, pc_y, pc_z, skd_90, skd_96, \
                         skd_135, skd_136, sld_132, sld_135, sld_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_9 * skd_96[k]
                   + f_3 * pc_y[k] * sld_132[k];

        t_222[k] = f_5 * skd_90[k]
                   + f_3 * pc_z[k] * sld_132[k];

        t_223[k] = f_8 * skd_135[k]
                   + f_3 * pc_x[k] * sld_135[k];

        t_224[k] = f_8 * skd_136[k]
                   + f_3 * pc_x[k] * sld_136[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, pb_z, pc_x, pc_y, pc_z, skf0_156, skd_93, \
                         skd_101, skd_137, skf1_156, sld_135, sld_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_8 * skd_137[k]
                   + f_3 * pc_x[k] * sld_137[k];

        t_226[k] = pb_z[k] * skf0_156[k]
                   - f_4 * pc_z[k] * skf1_156[k];

        t_227[k] = f_5 * skd_93[k]
                   + f_3 * pc_z[k] * sld_135[k];

        t_228[k] = f_9 * skd_101[k]
                   + f_3 * pc_y[k] * sld_137[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, pc_x, pc_y, pc_z, skd_95, skd_102, skd_138, \
                         slp0_68, slp0_69, slp1_68, slp1_69, sld_137, \
                         sld_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = f_5 * skd_95[k]
                   + f_1 * slp0_68[k]
                   - f_2 * slp1_68[k]
                   + f_3 * pc_z[k] * sld_137[k];

        t_230[k] = f_8 * skd_138[k]
                   + f_1 * slp0_69[k]
                   - f_2 * slp1_69[k]
                   + f_3 * pc_x[k] * sld_138[k];

        t_231[k] = f_11 * skd_102[k]
                   + f_3 * pc_y[k] * sld_138[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, t_235, pc_x, pc_z, skd_96, skd_141, skd_142, \
                         skd_143, sld_138, sld_141, sld_142, sld_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = f_8 * skd_96[k]
                   + f_3 * pc_z[k] * sld_138[k];

        t_233[k] = f_8 * skd_141[k]
                   + f_3 * pc_x[k] * sld_141[k];

        t_234[k] = f_8 * skd_142[k]
                   + f_3 * pc_x[k] * sld_142[k];

        t_235[k] = f_8 * skd_143[k]
                   + f_3 * pc_x[k] * sld_143[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, t_239, pc_y, pc_z, skd_99, skd_101, skd_105, \
                         skd_107, slp0_70, slp0_71, slp1_70, slp1_71, sld_141, \
                         sld_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = f_11 * skd_105[k]
                   + f_1 * slp0_70[k]
                   - f_2 * slp1_70[k]
                   + f_3 * pc_y[k] * sld_141[k];

        t_237[k] = f_8 * skd_99[k]
                   + f_3 * pc_z[k] * sld_141[k];

        t_238[k] = f_11 * skd_107[k]
                   + f_3 * pc_y[k] * sld_143[k];

        t_239[k] = f_8 * skd_101[k]
                   + f_1 * slp0_71[k]
                   - f_2 * slp1_71[k]
                   + f_3 * pc_z[k] * sld_143[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, pc_x, pc_y, pc_z, skd_102, skd_108, \
                         skd_144, skd_147, slp0_72, slp1_72, sld_144, \
                         sld_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_8 * skd_144[k]
                   + f_1 * slp0_72[k]
                   - f_2 * slp1_72[k]
                   + f_3 * pc_x[k] * sld_144[k];

        t_241[k] = f_10 * skd_108[k]
                   + f_3 * pc_y[k] * sld_144[k];

        t_242[k] = f_10 * skd_102[k]
                   + f_3 * pc_z[k] * sld_144[k];

        t_243[k] = f_8 * skd_147[k]
                   + f_3 * pc_x[k] * sld_147[k];
    }

#pragma omp simd aligned(t_244, t_245, t_246, t_247, pc_x, pc_y, pc_z, skd_105, skd_111, \
                         skd_148, skd_149, slp0_73, slp1_73, sld_147, sld_148, \
                         sld_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_244[k] = f_8 * skd_148[k]
                   + f_3 * pc_x[k] * sld_148[k];

        t_245[k] = f_8 * skd_149[k]
                   + f_3 * pc_x[k] * sld_149[k];

        t_246[k] = f_10 * skd_111[k]
                   + f_1 * slp0_73[k]
                   - f_2 * slp1_73[k]
                   + f_3 * pc_y[k] * sld_147[k];

        t_247[k] = f_10 * skd_105[k]
                   + f_3 * pc_z[k] * sld_147[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, pc_x, pc_y, pc_z, skd_107, skd_113, skd_150, \
                         slp0_74, slp0_75, slp1_74, slp1_75, sld_149, \
                         sld_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_10 * skd_113[k]
                   + f_3 * pc_y[k] * sld_149[k];

        t_249[k] = f_10 * skd_107[k]
                   + f_1 * slp0_74[k]
                   - f_2 * slp1_74[k]
                   + f_3 * pc_z[k] * sld_149[k];

        t_250[k] = f_8 * skd_150[k]
                   + f_1 * slp0_75[k]
                   - f_2 * slp1_75[k]
                   + f_3 * pc_x[k] * sld_150[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, pc_x, pc_y, pc_z, skd_108, skd_114, \
                         skd_153, skd_154, sld_150, sld_153, sld_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = f_8 * skd_114[k]
                   + f_3 * pc_y[k] * sld_150[k];

        t_252[k] = f_11 * skd_108[k]
                   + f_3 * pc_z[k] * sld_150[k];

        t_253[k] = f_8 * skd_153[k]
                   + f_3 * pc_x[k] * sld_153[k];

        t_254[k] = f_8 * skd_154[k]
                   + f_3 * pc_x[k] * sld_154[k];
    }
}

static auto
compute_prim_slf_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t skf0,
                                                          const size_t skd, const size_t skf1,
                                                          const size_t slp0, const size_t slp1,
                                                          const size_t sld, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 3.5 / q;
    const auto f_7 = 3.0 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 2.5 / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 2.0 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *skf0_200 = buffer.data(skf0 + 200);
    const auto *skf0_209 = buffer.data(skf0 + 209);
    const auto *skf0_210 = buffer.data(skf0 + 210);
    const auto *skf0_270 = buffer.data(skf0 + 270);
    const auto *skf0_280 = buffer.data(skf0 + 280);
    const auto *skf0_286 = buffer.data(skf0 + 286);
    const auto *skf0_289 = buffer.data(skf0 + 289);
    const auto *skf0_296 = buffer.data(skf0 + 296);
    const auto *skf0_299 = buffer.data(skf0 + 299);
    const auto *skf0_300 = buffer.data(skf0 + 300);
    const auto *skf0_306 = buffer.data(skf0 + 306);
    const auto *skf0_309 = buffer.data(skf0 + 309);
    const auto *skf0_310 = buffer.data(skf0 + 310);
    const auto *skf0_316 = buffer.data(skf0 + 316);
    const auto *skf0_319 = buffer.data(skf0 + 319);
    const auto *skf0_320 = buffer.data(skf0 + 320);
    const auto *skf0_326 = buffer.data(skf0 + 326);
    const auto *skf0_329 = buffer.data(skf0 + 329);
    const auto *skf0_330 = buffer.data(skf0 + 330);
    const auto *skf0_336 = buffer.data(skf0 + 336);
    const auto *skf0_339 = buffer.data(skf0 + 339);
    const auto *skf0_346 = buffer.data(skf0 + 346);
    const auto *skf0_349 = buffer.data(skf0 + 349);
    const auto *skf0_350 = buffer.data(skf0 + 350);
    const auto *skf0_356 = buffer.data(skf0 + 356);
    const auto *skf0_359 = buffer.data(skf0 + 359);

    const auto *skd_111 = buffer.data(skd + 111);
    const auto *skd_113 = buffer.data(skd + 113);
    const auto *skd_114 = buffer.data(skd + 114);
    const auto *skd_117 = buffer.data(skd + 117);
    const auto *skd_119 = buffer.data(skd + 119);
    const auto *skd_120 = buffer.data(skd + 120);
    const auto *skd_123 = buffer.data(skd + 123);
    const auto *skd_125 = buffer.data(skd + 125);
    const auto *skd_126 = buffer.data(skd + 126);
    const auto *skd_129 = buffer.data(skd + 129);
    const auto *skd_131 = buffer.data(skd + 131);
    const auto *skd_132 = buffer.data(skd + 132);
    const auto *skd_135 = buffer.data(skd + 135);
    const auto *skd_137 = buffer.data(skd + 137);
    const auto *skd_138 = buffer.data(skd + 138);
    const auto *skd_141 = buffer.data(skd + 141);
    const auto *skd_143 = buffer.data(skd + 143);
    const auto *skd_144 = buffer.data(skd + 144);
    const auto *skd_147 = buffer.data(skd + 147);
    const auto *skd_149 = buffer.data(skd + 149);
    const auto *skd_150 = buffer.data(skd + 150);
    const auto *skd_153 = buffer.data(skd + 153);
    const auto *skd_155 = buffer.data(skd + 155);
    const auto *skd_156 = buffer.data(skd + 156);
    const auto *skd_159 = buffer.data(skd + 159);
    const auto *skd_160 = buffer.data(skd + 160);
    const auto *skd_161 = buffer.data(skd + 161);
    const auto *skd_162 = buffer.data(skd + 162);
    const auto *skd_165 = buffer.data(skd + 165);
    const auto *skd_166 = buffer.data(skd + 166);
    const auto *skd_167 = buffer.data(skd + 167);
    const auto *skd_168 = buffer.data(skd + 168);
    const auto *skd_171 = buffer.data(skd + 171);
    const auto *skd_172 = buffer.data(skd + 172);
    const auto *skd_173 = buffer.data(skd + 173);
    const auto *skd_174 = buffer.data(skd + 174);
    const auto *skd_177 = buffer.data(skd + 177);
    const auto *skd_178 = buffer.data(skd + 178);
    const auto *skd_179 = buffer.data(skd + 179);
    const auto *skd_180 = buffer.data(skd + 180);
    const auto *skd_183 = buffer.data(skd + 183);
    const auto *skd_184 = buffer.data(skd + 184);
    const auto *skd_185 = buffer.data(skd + 185);
    const auto *skd_186 = buffer.data(skd + 186);
    const auto *skd_189 = buffer.data(skd + 189);
    const auto *skd_190 = buffer.data(skd + 190);
    const auto *skd_191 = buffer.data(skd + 191);
    const auto *skd_192 = buffer.data(skd + 192);
    const auto *skd_195 = buffer.data(skd + 195);
    const auto *skd_196 = buffer.data(skd + 196);
    const auto *skd_197 = buffer.data(skd + 197);
    const auto *skd_198 = buffer.data(skd + 198);
    const auto *skd_201 = buffer.data(skd + 201);
    const auto *skd_202 = buffer.data(skd + 202);
    const auto *skd_203 = buffer.data(skd + 203);
    const auto *skd_207 = buffer.data(skd + 207);
    const auto *skd_208 = buffer.data(skd + 208);
    const auto *skd_209 = buffer.data(skd + 209);
    const auto *skd_210 = buffer.data(skd + 210);
    const auto *skd_213 = buffer.data(skd + 213);
    const auto *skd_214 = buffer.data(skd + 214);
    const auto *skd_215 = buffer.data(skd + 215);

    const auto *skf1_200 = buffer.data(skf1 + 200);
    const auto *skf1_209 = buffer.data(skf1 + 209);
    const auto *skf1_210 = buffer.data(skf1 + 210);
    const auto *skf1_270 = buffer.data(skf1 + 270);
    const auto *skf1_280 = buffer.data(skf1 + 280);
    const auto *skf1_286 = buffer.data(skf1 + 286);
    const auto *skf1_289 = buffer.data(skf1 + 289);
    const auto *skf1_296 = buffer.data(skf1 + 296);
    const auto *skf1_299 = buffer.data(skf1 + 299);
    const auto *skf1_300 = buffer.data(skf1 + 300);
    const auto *skf1_306 = buffer.data(skf1 + 306);
    const auto *skf1_309 = buffer.data(skf1 + 309);
    const auto *skf1_310 = buffer.data(skf1 + 310);
    const auto *skf1_316 = buffer.data(skf1 + 316);
    const auto *skf1_319 = buffer.data(skf1 + 319);
    const auto *skf1_320 = buffer.data(skf1 + 320);
    const auto *skf1_326 = buffer.data(skf1 + 326);
    const auto *skf1_329 = buffer.data(skf1 + 329);
    const auto *skf1_330 = buffer.data(skf1 + 330);
    const auto *skf1_336 = buffer.data(skf1 + 336);
    const auto *skf1_339 = buffer.data(skf1 + 339);
    const auto *skf1_346 = buffer.data(skf1 + 346);
    const auto *skf1_349 = buffer.data(skf1 + 349);
    const auto *skf1_350 = buffer.data(skf1 + 350);
    const auto *skf1_356 = buffer.data(skf1 + 356);
    const auto *skf1_359 = buffer.data(skf1 + 359);

    const auto *slp0_76 = buffer.data(slp0 + 76);
    const auto *slp0_77 = buffer.data(slp0 + 77);
    const auto *slp0_79 = buffer.data(slp0 + 79);
    const auto *slp0_81 = buffer.data(slp0 + 81);
    const auto *slp0_82 = buffer.data(slp0 + 82);
    const auto *slp0_83 = buffer.data(slp0 + 83);
    const auto *slp0_108 = buffer.data(slp0 + 108);
    const auto *slp0_109 = buffer.data(slp0 + 109);
    const auto *slp0_110 = buffer.data(slp0 + 110);
    const auto *slp0_113 = buffer.data(slp0 + 113);
    const auto *slp0_114 = buffer.data(slp0 + 114);
    const auto *slp0_115 = buffer.data(slp0 + 115);

    const auto *slp1_76 = buffer.data(slp1 + 76);
    const auto *slp1_77 = buffer.data(slp1 + 77);
    const auto *slp1_79 = buffer.data(slp1 + 79);
    const auto *slp1_81 = buffer.data(slp1 + 81);
    const auto *slp1_82 = buffer.data(slp1 + 82);
    const auto *slp1_83 = buffer.data(slp1 + 83);
    const auto *slp1_108 = buffer.data(slp1 + 108);
    const auto *slp1_109 = buffer.data(slp1 + 109);
    const auto *slp1_110 = buffer.data(slp1 + 110);
    const auto *slp1_113 = buffer.data(slp1 + 113);
    const auto *slp1_114 = buffer.data(slp1 + 114);
    const auto *slp1_115 = buffer.data(slp1 + 115);

    const auto *sld_153 = buffer.data(sld + 153);
    const auto *sld_155 = buffer.data(sld + 155);
    const auto *sld_156 = buffer.data(sld + 156);
    const auto *sld_159 = buffer.data(sld + 159);
    const auto *sld_160 = buffer.data(sld + 160);
    const auto *sld_161 = buffer.data(sld + 161);
    const auto *sld_162 = buffer.data(sld + 162);
    const auto *sld_165 = buffer.data(sld + 165);
    const auto *sld_166 = buffer.data(sld + 166);
    const auto *sld_167 = buffer.data(sld + 167);
    const auto *sld_168 = buffer.data(sld + 168);
    const auto *sld_171 = buffer.data(sld + 171);
    const auto *sld_172 = buffer.data(sld + 172);
    const auto *sld_173 = buffer.data(sld + 173);
    const auto *sld_174 = buffer.data(sld + 174);
    const auto *sld_177 = buffer.data(sld + 177);
    const auto *sld_178 = buffer.data(sld + 178);
    const auto *sld_179 = buffer.data(sld + 179);
    const auto *sld_180 = buffer.data(sld + 180);
    const auto *sld_183 = buffer.data(sld + 183);
    const auto *sld_184 = buffer.data(sld + 184);
    const auto *sld_185 = buffer.data(sld + 185);
    const auto *sld_186 = buffer.data(sld + 186);
    const auto *sld_189 = buffer.data(sld + 189);
    const auto *sld_190 = buffer.data(sld + 190);
    const auto *sld_191 = buffer.data(sld + 191);
    const auto *sld_192 = buffer.data(sld + 192);
    const auto *sld_195 = buffer.data(sld + 195);
    const auto *sld_196 = buffer.data(sld + 196);
    const auto *sld_197 = buffer.data(sld + 197);
    const auto *sld_198 = buffer.data(sld + 198);
    const auto *sld_201 = buffer.data(sld + 201);
    const auto *sld_202 = buffer.data(sld + 202);
    const auto *sld_203 = buffer.data(sld + 203);
    const auto *sld_204 = buffer.data(sld + 204);
    const auto *sld_207 = buffer.data(sld + 207);
    const auto *sld_208 = buffer.data(sld + 208);
    const auto *sld_209 = buffer.data(sld + 209);
    const auto *sld_210 = buffer.data(sld + 210);
    const auto *sld_213 = buffer.data(sld + 213);
    const auto *sld_214 = buffer.data(sld + 214);
    const auto *sld_215 = buffer.data(sld + 215);
    const auto *sld_216 = buffer.data(sld + 216);
    const auto *sld_219 = buffer.data(sld + 219);
    const auto *sld_220 = buffer.data(sld + 220);
    const auto *sld_221 = buffer.data(sld + 221);
    const auto *sld_222 = buffer.data(sld + 222);
    const auto *sld_225 = buffer.data(sld + 225);
    const auto *sld_226 = buffer.data(sld + 226);
    const auto *sld_227 = buffer.data(sld + 227);
    const auto *sld_228 = buffer.data(sld + 228);
    const auto *sld_231 = buffer.data(sld + 231);
    const auto *sld_232 = buffer.data(sld + 232);
    const auto *sld_233 = buffer.data(sld + 233);

#pragma omp simd aligned(t_255, t_256, t_257, t_258, pc_x, pc_y, pc_z, skd_111, skd_117, \
                         skd_119, skd_155, slp0_76, slp1_76, sld_153, \
                         sld_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_8 * skd_155[k]
                   + f_3 * pc_x[k] * sld_155[k];

        t_256[k] = f_8 * skd_117[k]
                   + f_1 * slp0_76[k]
                   - f_2 * slp1_76[k]
                   + f_3 * pc_y[k] * sld_153[k];

        t_257[k] = f_11 * skd_111[k]
                   + f_3 * pc_z[k] * sld_153[k];

        t_258[k] = f_8 * skd_119[k]
                   + f_3 * pc_y[k] * sld_155[k];
    }

#pragma omp simd aligned(t_259, t_260, t_261, t_262, pb_y, pc_y, pc_z, skf0_200, skd_113, \
                         skd_114, skd_120, skf1_200, slp0_77, slp1_77, sld_155, \
                         sld_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_259[k] = f_11 * skd_113[k]
                   + f_1 * slp0_77[k]
                   - f_2 * slp1_77[k]
                   + f_3 * pc_z[k] * sld_155[k];

        t_260[k] = pb_y[k] * skf0_200[k]
                   - f_4 * pc_y[k] * skf1_200[k];

        t_261[k] = f_5 * skd_120[k]
                   + f_3 * pc_y[k] * sld_156[k];

        t_262[k] = f_9 * skd_114[k]
                   + f_3 * pc_z[k] * sld_156[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, t_266, pc_x, pc_y, skd_123, skd_159, skd_160, \
                         skd_161, slp0_79, slp1_79, sld_159, sld_160, \
                         sld_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = f_8 * skd_159[k]
                   + f_3 * pc_x[k] * sld_159[k];

        t_264[k] = f_8 * skd_160[k]
                   + f_3 * pc_x[k] * sld_160[k];

        t_265[k] = f_8 * skd_161[k]
                   + f_3 * pc_x[k] * sld_161[k];

        t_266[k] = f_5 * skd_123[k]
                   + f_1 * slp0_79[k]
                   - f_2 * slp1_79[k]
                   + f_3 * pc_y[k] * sld_159[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, pb_y, pc_y, pc_z, skf0_209, skd_117, skd_125, \
                         skf1_209, sld_159, sld_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_9 * skd_117[k]
                   + f_3 * pc_z[k] * sld_159[k];

        t_268[k] = f_5 * skd_125[k]
                   + f_3 * pc_y[k] * sld_161[k];

        t_269[k] = pb_y[k] * skf0_209[k]
                   - f_4 * pc_y[k] * skf1_209[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, pc_x, pc_y, pc_z, skd_120, skd_162, \
                         skd_165, slp0_81, slp1_81, sld_162, sld_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_8 * skd_162[k]
                   + f_1 * slp0_81[k]
                   - f_2 * slp1_81[k]
                   + f_3 * pc_x[k] * sld_162[k];

        t_271[k] = f_3 * pc_y[k] * sld_162[k];

        t_272[k] = f_7 * skd_120[k]
                   + f_3 * pc_z[k] * sld_162[k];

        t_273[k] = f_8 * skd_165[k]
                   + f_3 * pc_x[k] * sld_165[k];
    }

#pragma omp simd aligned(t_274, t_275, t_276, t_277, t_278, pc_x, pc_y, pc_z, skd_123, \
                         skd_166, skd_167, slp0_82, slp1_82, sld_165, sld_166, \
                         sld_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_274[k] = f_8 * skd_166[k]
                   + f_3 * pc_x[k] * sld_166[k];

        t_275[k] = f_8 * skd_167[k]
                   + f_3 * pc_x[k] * sld_167[k];

        t_276[k] = f_1 * slp0_82[k]
                   - f_2 * slp1_82[k]
                   + f_3 * pc_y[k] * sld_165[k];

        t_277[k] = f_7 * skd_123[k]
                   + f_3 * pc_z[k] * sld_165[k];

        t_278[k] = f_3 * pc_y[k] * sld_167[k];
    }

#pragma omp simd aligned(t_279, t_280, t_281, pb_x, pc_x, pc_y, pc_z, skf0_280, skd_125, \
                         skd_126, skd_168, skf1_280, slp0_83, slp1_83, sld_167, \
                         sld_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_279[k] = f_7 * skd_125[k]
                   + f_1 * slp0_83[k]
                   - f_2 * slp1_83[k]
                   + f_3 * pc_z[k] * sld_167[k];

        t_280[k] = pb_x[k] * skf0_280[k]
                   + f_10 * skd_168[k]
                   - f_4 * pc_x[k] * skf1_280[k];

        t_281[k] = f_6 * skd_126[k]
                   + f_3 * pc_y[k] * sld_168[k];
    }

#pragma omp simd aligned(t_282, t_283, t_284, t_285, pc_x, pc_z, skd_171, skd_172, skd_173, \
                         sld_168, sld_171, sld_172, sld_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_282[k] = f_3 * pc_z[k] * sld_168[k];

        t_283[k] = f_5 * skd_171[k]
                   + f_3 * pc_x[k] * sld_171[k];

        t_284[k] = f_5 * skd_172[k]
                   + f_3 * pc_x[k] * sld_172[k];

        t_285[k] = f_5 * skd_173[k]
                   + f_3 * pc_x[k] * sld_173[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, pb_x, pc_x, pc_y, pc_z, skf0_286, \
                         skf0_289, skd_131, skf1_286, skf1_289, sld_171, \
                         sld_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = pb_x[k] * skf0_286[k]
                   - f_4 * pc_x[k] * skf1_286[k];

        t_287[k] = f_3 * pc_z[k] * sld_171[k];

        t_288[k] = f_6 * skd_131[k]
                   + f_3 * pc_y[k] * sld_173[k];

        t_289[k] = pb_x[k] * skf0_289[k]
                   - f_4 * pc_x[k] * skf1_289[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, pb_z, pc_x, pc_y, pc_z, skf0_210, \
                         skd_126, skd_132, skd_177, skf1_210, sld_174, \
                         sld_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = pb_z[k] * skf0_210[k]
                   - f_4 * pc_z[k] * skf1_210[k];

        t_291[k] = f_7 * skd_132[k]
                   + f_3 * pc_y[k] * sld_174[k];

        t_292[k] = f_5 * skd_126[k]
                   + f_3 * pc_z[k] * sld_174[k];

        t_293[k] = f_5 * skd_177[k]
                   + f_3 * pc_x[k] * sld_177[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, pb_x, pc_x, pc_z, skf0_296, skd_129, \
                         skd_178, skd_179, skf1_296, sld_177, sld_178, \
                         sld_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_5 * skd_178[k]
                   + f_3 * pc_x[k] * sld_178[k];

        t_295[k] = f_5 * skd_179[k]
                   + f_3 * pc_x[k] * sld_179[k];

        t_296[k] = pb_x[k] * skf0_296[k]
                   - f_4 * pc_x[k] * skf1_296[k];

        t_297[k] = f_5 * skd_129[k]
                   + f_3 * pc_z[k] * sld_177[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, t_301, pb_x, pc_x, pc_y, skf0_299, skf0_300, \
                         skd_137, skd_138, skd_180, skf1_299, skf1_300, sld_179, \
                         sld_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = f_7 * skd_137[k]
                   + f_3 * pc_y[k] * sld_179[k];

        t_299[k] = pb_x[k] * skf0_299[k]
                   - f_4 * pc_x[k] * skf1_299[k];

        t_300[k] = pb_x[k] * skf0_300[k]
                   + f_10 * skd_180[k]
                   - f_4 * pc_x[k] * skf1_300[k];

        t_301[k] = f_9 * skd_138[k]
                   + f_3 * pc_y[k] * sld_180[k];
    }

#pragma omp simd aligned(t_302, t_303, t_304, t_305, pc_x, pc_z, skd_132, skd_183, skd_184, \
                         skd_185, sld_180, sld_183, sld_184, sld_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_302[k] = f_8 * skd_132[k]
                   + f_3 * pc_z[k] * sld_180[k];

        t_303[k] = f_5 * skd_183[k]
                   + f_3 * pc_x[k] * sld_183[k];

        t_304[k] = f_5 * skd_184[k]
                   + f_3 * pc_x[k] * sld_184[k];

        t_305[k] = f_5 * skd_185[k]
                   + f_3 * pc_x[k] * sld_185[k];
    }

#pragma omp simd aligned(t_306, t_307, t_308, t_309, pb_x, pc_x, pc_y, pc_z, skf0_306, \
                         skf0_309, skd_135, skd_143, skf1_306, skf1_309, sld_183, \
                         sld_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_306[k] = pb_x[k] * skf0_306[k]
                   - f_4 * pc_x[k] * skf1_306[k];

        t_307[k] = f_8 * skd_135[k]
                   + f_3 * pc_z[k] * sld_183[k];

        t_308[k] = f_9 * skd_143[k]
                   + f_3 * pc_y[k] * sld_185[k];

        t_309[k] = pb_x[k] * skf0_309[k]
                   - f_4 * pc_x[k] * skf1_309[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, pb_x, pc_x, pc_y, pc_z, skf0_310, \
                         skd_138, skd_144, skd_186, skd_189, skf1_310, sld_186, \
                         sld_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = pb_x[k] * skf0_310[k]
                   + f_10 * skd_186[k]
                   - f_4 * pc_x[k] * skf1_310[k];

        t_311[k] = f_11 * skd_144[k]
                   + f_3 * pc_y[k] * sld_186[k];

        t_312[k] = f_10 * skd_138[k]
                   + f_3 * pc_z[k] * sld_186[k];

        t_313[k] = f_5 * skd_189[k]
                   + f_3 * pc_x[k] * sld_189[k];
    }

#pragma omp simd aligned(t_314, t_315, t_316, t_317, pb_x, pc_x, pc_z, skf0_316, skd_141, \
                         skd_190, skd_191, skf1_316, sld_189, sld_190, \
                         sld_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_314[k] = f_5 * skd_190[k]
                   + f_3 * pc_x[k] * sld_190[k];

        t_315[k] = f_5 * skd_191[k]
                   + f_3 * pc_x[k] * sld_191[k];

        t_316[k] = pb_x[k] * skf0_316[k]
                   - f_4 * pc_x[k] * skf1_316[k];

        t_317[k] = f_10 * skd_141[k]
                   + f_3 * pc_z[k] * sld_189[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, t_321, pb_x, pc_x, pc_y, skf0_319, skf0_320, \
                         skd_149, skd_150, skd_192, skf1_319, skf1_320, sld_191, \
                         sld_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = f_11 * skd_149[k]
                   + f_3 * pc_y[k] * sld_191[k];

        t_319[k] = pb_x[k] * skf0_319[k]
                   - f_4 * pc_x[k] * skf1_319[k];

        t_320[k] = pb_x[k] * skf0_320[k]
                   + f_10 * skd_192[k]
                   - f_4 * pc_x[k] * skf1_320[k];

        t_321[k] = f_10 * skd_150[k]
                   + f_3 * pc_y[k] * sld_192[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, t_325, pc_x, pc_z, skd_144, skd_195, skd_196, \
                         skd_197, sld_192, sld_195, sld_196, sld_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = f_11 * skd_144[k]
                   + f_3 * pc_z[k] * sld_192[k];

        t_323[k] = f_5 * skd_195[k]
                   + f_3 * pc_x[k] * sld_195[k];

        t_324[k] = f_5 * skd_196[k]
                   + f_3 * pc_x[k] * sld_196[k];

        t_325[k] = f_5 * skd_197[k]
                   + f_3 * pc_x[k] * sld_197[k];
    }

#pragma omp simd aligned(t_326, t_327, t_328, t_329, pb_x, pc_x, pc_y, pc_z, skf0_326, \
                         skf0_329, skd_147, skd_155, skf1_326, skf1_329, sld_195, \
                         sld_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_326[k] = pb_x[k] * skf0_326[k]
                   - f_4 * pc_x[k] * skf1_326[k];

        t_327[k] = f_11 * skd_147[k]
                   + f_3 * pc_z[k] * sld_195[k];

        t_328[k] = f_10 * skd_155[k]
                   + f_3 * pc_y[k] * sld_197[k];

        t_329[k] = pb_x[k] * skf0_329[k]
                   - f_4 * pc_x[k] * skf1_329[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, pb_x, pc_x, pc_y, pc_z, skf0_330, \
                         skd_150, skd_156, skd_198, skd_201, skf1_330, sld_198, \
                         sld_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = pb_x[k] * skf0_330[k]
                   + f_10 * skd_198[k]
                   - f_4 * pc_x[k] * skf1_330[k];

        t_331[k] = f_8 * skd_156[k]
                   + f_3 * pc_y[k] * sld_198[k];

        t_332[k] = f_9 * skd_150[k]
                   + f_3 * pc_z[k] * sld_198[k];

        t_333[k] = f_5 * skd_201[k]
                   + f_3 * pc_x[k] * sld_201[k];
    }

#pragma omp simd aligned(t_334, t_335, t_336, t_337, pb_x, pc_x, pc_z, skf0_336, skd_153, \
                         skd_202, skd_203, skf1_336, sld_201, sld_202, \
                         sld_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_334[k] = f_5 * skd_202[k]
                   + f_3 * pc_x[k] * sld_202[k];

        t_335[k] = f_5 * skd_203[k]
                   + f_3 * pc_x[k] * sld_203[k];

        t_336[k] = pb_x[k] * skf0_336[k]
                   - f_4 * pc_x[k] * skf1_336[k];

        t_337[k] = f_9 * skd_153[k]
                   + f_3 * pc_z[k] * sld_201[k];
    }

#pragma omp simd aligned(t_338, t_339, t_340, t_341, pb_x, pb_y, pc_x, pc_y, skf0_270, \
                         skf0_339, skd_161, skd_162, skf1_270, skf1_339, sld_203, \
                         sld_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_338[k] = f_8 * skd_161[k]
                   + f_3 * pc_y[k] * sld_203[k];

        t_339[k] = pb_x[k] * skf0_339[k]
                   - f_4 * pc_x[k] * skf1_339[k];

        t_340[k] = pb_y[k] * skf0_270[k]
                   - f_4 * pc_y[k] * skf1_270[k];

        t_341[k] = f_5 * skd_162[k]
                   + f_3 * pc_y[k] * sld_204[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, t_345, pc_x, pc_z, skd_156, skd_207, skd_208, \
                         skd_209, sld_204, sld_207, sld_208, sld_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = f_7 * skd_156[k]
                   + f_3 * pc_z[k] * sld_204[k];

        t_343[k] = f_5 * skd_207[k]
                   + f_3 * pc_x[k] * sld_207[k];

        t_344[k] = f_5 * skd_208[k]
                   + f_3 * pc_x[k] * sld_208[k];

        t_345[k] = f_5 * skd_209[k]
                   + f_3 * pc_x[k] * sld_209[k];
    }

#pragma omp simd aligned(t_346, t_347, t_348, t_349, pb_x, pc_x, pc_y, pc_z, skf0_346, \
                         skf0_349, skd_159, skd_167, skf1_346, skf1_349, sld_207, \
                         sld_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_346[k] = pb_x[k] * skf0_346[k]
                   - f_4 * pc_x[k] * skf1_346[k];

        t_347[k] = f_7 * skd_159[k]
                   + f_3 * pc_z[k] * sld_207[k];

        t_348[k] = f_5 * skd_167[k]
                   + f_3 * pc_y[k] * sld_209[k];

        t_349[k] = pb_x[k] * skf0_349[k]
                   - f_4 * pc_x[k] * skf1_349[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, pb_x, pc_x, pc_y, pc_z, skf0_350, \
                         skd_162, skd_210, skd_213, skf1_350, sld_210, \
                         sld_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = pb_x[k] * skf0_350[k]
                   + f_10 * skd_210[k]
                   - f_4 * pc_x[k] * skf1_350[k];

        t_351[k] = f_3 * pc_y[k] * sld_210[k];

        t_352[k] = f_6 * skd_162[k]
                   + f_3 * pc_z[k] * sld_210[k];

        t_353[k] = f_5 * skd_213[k]
                   + f_3 * pc_x[k] * sld_213[k];
    }

#pragma omp simd aligned(t_354, t_355, t_356, t_357, pb_x, pc_x, pc_z, skf0_356, skd_165, \
                         skd_214, skd_215, skf1_356, sld_213, sld_214, \
                         sld_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_354[k] = f_5 * skd_214[k]
                   + f_3 * pc_x[k] * sld_214[k];

        t_355[k] = f_5 * skd_215[k]
                   + f_3 * pc_x[k] * sld_215[k];

        t_356[k] = pb_x[k] * skf0_356[k]
                   - f_4 * pc_x[k] * skf1_356[k];

        t_357[k] = f_6 * skd_165[k]
                   + f_3 * pc_z[k] * sld_213[k];
    }

#pragma omp simd aligned(t_358, t_359, t_360, t_361, t_362, pb_x, pc_x, pc_y, pc_z, skf0_359, \
                         skd_168, skf1_359, slp0_108, slp1_108, sld_215, \
                         sld_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_358[k] = f_3 * pc_y[k] * sld_215[k];

        t_359[k] = pb_x[k] * skf0_359[k]
                   - f_4 * pc_x[k] * skf1_359[k];

        t_360[k] = f_1 * slp0_108[k]
                   - f_2 * slp1_108[k]
                   + f_3 * pc_x[k] * sld_216[k];

        t_361[k] = f_0 * skd_168[k]
                   + f_3 * pc_y[k] * sld_216[k];

        t_362[k] = f_3 * pc_z[k] * sld_216[k];
    }

#pragma omp simd aligned(t_363, t_364, t_365, t_366, t_367, t_368, pc_x, pc_y, pc_z, skd_171, \
                         skd_173, slp0_109, slp1_109, sld_219, sld_220, \
                         sld_221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = f_3 * pc_x[k] * sld_219[k];

        t_364[k] = f_3 * pc_x[k] * sld_220[k];

        t_365[k] = f_3 * pc_x[k] * sld_221[k];

        t_366[k] = f_0 * skd_171[k]
                   + f_1 * slp0_109[k]
                   - f_2 * slp1_109[k]
                   + f_3 * pc_y[k] * sld_219[k];

        t_367[k] = f_3 * pc_z[k] * sld_219[k];

        t_368[k] = f_0 * skd_173[k]
                   + f_3 * pc_y[k] * sld_221[k];
    }

#pragma omp simd aligned(t_369, t_370, t_371, t_372, pb_z, pc_y, pc_z, skf0_280, skd_168, \
                         skd_174, skf1_280, slp0_110, slp1_110, sld_221, \
                         sld_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_369[k] = f_1 * slp0_110[k]
                   - f_2 * slp1_110[k]
                   + f_3 * pc_z[k] * sld_221[k];

        t_370[k] = pb_z[k] * skf0_280[k]
                   - f_4 * pc_z[k] * skf1_280[k];

        t_371[k] = f_6 * skd_174[k]
                   + f_3 * pc_y[k] * sld_222[k];

        t_372[k] = f_5 * skd_168[k]
                   + f_3 * pc_z[k] * sld_222[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, t_376, t_377, pb_z, pc_x, pc_z, skf0_286, \
                         skd_171, skf1_286, sld_225, sld_226, sld_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = f_3 * pc_x[k] * sld_225[k];

        t_374[k] = f_3 * pc_x[k] * sld_226[k];

        t_375[k] = f_3 * pc_x[k] * sld_227[k];

        t_376[k] = pb_z[k] * skf0_286[k]
                   - f_4 * pc_z[k] * skf1_286[k];

        t_377[k] = f_5 * skd_171[k]
                   + f_3 * pc_z[k] * sld_225[k];
    }

#pragma omp simd aligned(t_378, t_379, t_380, t_381, pc_x, pc_y, pc_z, skd_173, skd_179, \
                         skd_180, slp0_113, slp0_114, slp1_113, slp1_114, sld_227, \
                         sld_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_378[k] = f_6 * skd_179[k]
                   + f_3 * pc_y[k] * sld_227[k];

        t_379[k] = f_5 * skd_173[k]
                   + f_1 * slp0_113[k]
                   - f_2 * slp1_113[k]
                   + f_3 * pc_z[k] * sld_227[k];

        t_380[k] = f_1 * slp0_114[k]
                   - f_2 * slp1_114[k]
                   + f_3 * pc_x[k] * sld_228[k];

        t_381[k] = f_7 * skd_180[k]
                   + f_3 * pc_y[k] * sld_228[k];
    }

#pragma omp simd aligned(t_382, t_383, t_384, t_385, t_386, pc_x, pc_y, pc_z, skd_174, \
                         skd_183, slp0_115, slp1_115, sld_228, sld_231, sld_232, \
                         sld_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_382[k] = f_8 * skd_174[k]
                   + f_3 * pc_z[k] * sld_228[k];

        t_383[k] = f_3 * pc_x[k] * sld_231[k];

        t_384[k] = f_3 * pc_x[k] * sld_232[k];

        t_385[k] = f_3 * pc_x[k] * sld_233[k];

        t_386[k] = f_7 * skd_183[k]
                   + f_1 * slp0_115[k]
                   - f_2 * slp1_115[k]
                   + f_3 * pc_y[k] * sld_231[k];
    }
}

static auto
compute_prim_slf_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t skf0,
                                                          const size_t skd, const size_t skf1,
                                                          const size_t slp0, const size_t slp1,
                                                          const size_t sld, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 3.5 / q;
    const auto f_7 = 3.0 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 2.5 / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 2.0 / q;

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

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *skf0_350 = buffer.data(skf0 + 350);
    const auto *skf0_356 = buffer.data(skf0 + 356);
    const auto *skf0_359 = buffer.data(skf0 + 359);

    const auto *skd_177 = buffer.data(skd + 177);
    const auto *skd_179 = buffer.data(skd + 179);
    const auto *skd_180 = buffer.data(skd + 180);
    const auto *skd_183 = buffer.data(skd + 183);
    const auto *skd_185 = buffer.data(skd + 185);
    const auto *skd_186 = buffer.data(skd + 186);
    const auto *skd_189 = buffer.data(skd + 189);
    const auto *skd_191 = buffer.data(skd + 191);
    const auto *skd_192 = buffer.data(skd + 192);
    const auto *skd_195 = buffer.data(skd + 195);
    const auto *skd_197 = buffer.data(skd + 197);
    const auto *skd_198 = buffer.data(skd + 198);
    const auto *skd_201 = buffer.data(skd + 201);
    const auto *skd_203 = buffer.data(skd + 203);
    const auto *skd_204 = buffer.data(skd + 204);
    const auto *skd_207 = buffer.data(skd + 207);
    const auto *skd_209 = buffer.data(skd + 209);
    const auto *skd_210 = buffer.data(skd + 210);
    const auto *skd_213 = buffer.data(skd + 213);
    const auto *skd_215 = buffer.data(skd + 215);

    const auto *skf1_350 = buffer.data(skf1 + 350);
    const auto *skf1_356 = buffer.data(skf1 + 356);
    const auto *skf1_359 = buffer.data(skf1 + 359);

    const auto *slp0_116 = buffer.data(slp0 + 116);
    const auto *slp0_117 = buffer.data(slp0 + 117);
    const auto *slp0_118 = buffer.data(slp0 + 118);
    const auto *slp0_119 = buffer.data(slp0 + 119);
    const auto *slp0_120 = buffer.data(slp0 + 120);
    const auto *slp0_121 = buffer.data(slp0 + 121);
    const auto *slp0_122 = buffer.data(slp0 + 122);
    const auto *slp0_123 = buffer.data(slp0 + 123);
    const auto *slp0_124 = buffer.data(slp0 + 124);
    const auto *slp0_125 = buffer.data(slp0 + 125);
    const auto *slp0_126 = buffer.data(slp0 + 126);
    const auto *slp0_127 = buffer.data(slp0 + 127);
    const auto *slp0_128 = buffer.data(slp0 + 128);
    const auto *slp0_132 = buffer.data(slp0 + 132);
    const auto *slp0_133 = buffer.data(slp0 + 133);
    const auto *slp0_134 = buffer.data(slp0 + 134);

    const auto *slp1_116 = buffer.data(slp1 + 116);
    const auto *slp1_117 = buffer.data(slp1 + 117);
    const auto *slp1_118 = buffer.data(slp1 + 118);
    const auto *slp1_119 = buffer.data(slp1 + 119);
    const auto *slp1_120 = buffer.data(slp1 + 120);
    const auto *slp1_121 = buffer.data(slp1 + 121);
    const auto *slp1_122 = buffer.data(slp1 + 122);
    const auto *slp1_123 = buffer.data(slp1 + 123);
    const auto *slp1_124 = buffer.data(slp1 + 124);
    const auto *slp1_125 = buffer.data(slp1 + 125);
    const auto *slp1_126 = buffer.data(slp1 + 126);
    const auto *slp1_127 = buffer.data(slp1 + 127);
    const auto *slp1_128 = buffer.data(slp1 + 128);
    const auto *slp1_132 = buffer.data(slp1 + 132);
    const auto *slp1_133 = buffer.data(slp1 + 133);
    const auto *slp1_134 = buffer.data(slp1 + 134);

    const auto *sld_231 = buffer.data(sld + 231);
    const auto *sld_233 = buffer.data(sld + 233);
    const auto *sld_234 = buffer.data(sld + 234);
    const auto *sld_237 = buffer.data(sld + 237);
    const auto *sld_238 = buffer.data(sld + 238);
    const auto *sld_239 = buffer.data(sld + 239);
    const auto *sld_240 = buffer.data(sld + 240);
    const auto *sld_243 = buffer.data(sld + 243);
    const auto *sld_244 = buffer.data(sld + 244);
    const auto *sld_245 = buffer.data(sld + 245);
    const auto *sld_246 = buffer.data(sld + 246);
    const auto *sld_249 = buffer.data(sld + 249);
    const auto *sld_250 = buffer.data(sld + 250);
    const auto *sld_251 = buffer.data(sld + 251);
    const auto *sld_252 = buffer.data(sld + 252);
    const auto *sld_255 = buffer.data(sld + 255);
    const auto *sld_256 = buffer.data(sld + 256);
    const auto *sld_257 = buffer.data(sld + 257);
    const auto *sld_258 = buffer.data(sld + 258);
    const auto *sld_261 = buffer.data(sld + 261);
    const auto *sld_262 = buffer.data(sld + 262);
    const auto *sld_263 = buffer.data(sld + 263);
    const auto *sld_264 = buffer.data(sld + 264);
    const auto *sld_267 = buffer.data(sld + 267);
    const auto *sld_268 = buffer.data(sld + 268);
    const auto *sld_269 = buffer.data(sld + 269);

#pragma omp simd aligned(t_387, t_388, t_389, pc_y, pc_z, skd_177, skd_179, skd_185, slp0_116, \
                         slp1_116, sld_231, sld_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_387[k] = f_8 * skd_177[k]
                   + f_3 * pc_z[k] * sld_231[k];

        t_388[k] = f_7 * skd_185[k]
                   + f_3 * pc_y[k] * sld_233[k];

        t_389[k] = f_8 * skd_179[k]
                   + f_1 * slp0_116[k]
                   - f_2 * slp1_116[k]
                   + f_3 * pc_z[k] * sld_233[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, pc_x, pc_y, pc_z, skd_180, \
                         skd_186, slp0_117, slp1_117, sld_234, sld_237, \
                         sld_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = f_1 * slp0_117[k]
                   - f_2 * slp1_117[k]
                   + f_3 * pc_x[k] * sld_234[k];

        t_391[k] = f_9 * skd_186[k]
                   + f_3 * pc_y[k] * sld_234[k];

        t_392[k] = f_10 * skd_180[k]
                   + f_3 * pc_z[k] * sld_234[k];

        t_393[k] = f_3 * pc_x[k] * sld_237[k];

        t_394[k] = f_3 * pc_x[k] * sld_238[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, pc_x, pc_y, pc_z, skd_183, skd_189, \
                         skd_191, slp0_118, slp1_118, sld_237, \
                         sld_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = f_3 * pc_x[k] * sld_239[k];

        t_396[k] = f_9 * skd_189[k]
                   + f_1 * slp0_118[k]
                   - f_2 * slp1_118[k]
                   + f_3 * pc_y[k] * sld_237[k];

        t_397[k] = f_10 * skd_183[k]
                   + f_3 * pc_z[k] * sld_237[k];

        t_398[k] = f_9 * skd_191[k]
                   + f_3 * pc_y[k] * sld_239[k];
    }

#pragma omp simd aligned(t_399, t_400, t_401, t_402, pc_x, pc_y, pc_z, skd_185, skd_186, \
                         skd_192, slp0_119, slp0_120, slp1_119, slp1_120, sld_239, \
                         sld_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_399[k] = f_10 * skd_185[k]
                   + f_1 * slp0_119[k]
                   - f_2 * slp1_119[k]
                   + f_3 * pc_z[k] * sld_239[k];

        t_400[k] = f_1 * slp0_120[k]
                   - f_2 * slp1_120[k]
                   + f_3 * pc_x[k] * sld_240[k];

        t_401[k] = f_11 * skd_192[k]
                   + f_3 * pc_y[k] * sld_240[k];

        t_402[k] = f_11 * skd_186[k]
                   + f_3 * pc_z[k] * sld_240[k];
    }

#pragma omp simd aligned(t_403, t_404, t_405, t_406, t_407, pc_x, pc_y, pc_z, skd_189, \
                         skd_195, slp0_121, slp1_121, sld_243, sld_244, \
                         sld_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_403[k] = f_3 * pc_x[k] * sld_243[k];

        t_404[k] = f_3 * pc_x[k] * sld_244[k];

        t_405[k] = f_3 * pc_x[k] * sld_245[k];

        t_406[k] = f_11 * skd_195[k]
                   + f_1 * slp0_121[k]
                   - f_2 * slp1_121[k]
                   + f_3 * pc_y[k] * sld_243[k];

        t_407[k] = f_11 * skd_189[k]
                   + f_3 * pc_z[k] * sld_243[k];
    }

#pragma omp simd aligned(t_408, t_409, t_410, t_411, pc_x, pc_y, pc_z, skd_191, skd_197, \
                         skd_198, slp0_122, slp0_123, slp1_122, slp1_123, sld_245, \
                         sld_246 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_408[k] = f_11 * skd_197[k]
                   + f_3 * pc_y[k] * sld_245[k];

        t_409[k] = f_11 * skd_191[k]
                   + f_1 * slp0_122[k]
                   - f_2 * slp1_122[k]
                   + f_3 * pc_z[k] * sld_245[k];

        t_410[k] = f_1 * slp0_123[k]
                   - f_2 * slp1_123[k]
                   + f_3 * pc_x[k] * sld_246[k];

        t_411[k] = f_10 * skd_198[k]
                   + f_3 * pc_y[k] * sld_246[k];
    }

#pragma omp simd aligned(t_412, t_413, t_414, t_415, t_416, pc_x, pc_y, pc_z, skd_192, \
                         skd_201, slp0_124, slp1_124, sld_246, sld_249, sld_250, \
                         sld_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_412[k] = f_9 * skd_192[k]
                   + f_3 * pc_z[k] * sld_246[k];

        t_413[k] = f_3 * pc_x[k] * sld_249[k];

        t_414[k] = f_3 * pc_x[k] * sld_250[k];

        t_415[k] = f_3 * pc_x[k] * sld_251[k];

        t_416[k] = f_10 * skd_201[k]
                   + f_1 * slp0_124[k]
                   - f_2 * slp1_124[k]
                   + f_3 * pc_y[k] * sld_249[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, pc_y, pc_z, skd_195, skd_197, skd_203, slp0_125, \
                         slp1_125, sld_249, sld_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = f_9 * skd_195[k]
                   + f_3 * pc_z[k] * sld_249[k];

        t_418[k] = f_10 * skd_203[k]
                   + f_3 * pc_y[k] * sld_251[k];

        t_419[k] = f_9 * skd_197[k]
                   + f_1 * slp0_125[k]
                   - f_2 * slp1_125[k]
                   + f_3 * pc_z[k] * sld_251[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, pc_x, pc_y, pc_z, skd_198, \
                         skd_204, slp0_126, slp1_126, sld_252, sld_255, \
                         sld_256 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_1 * slp0_126[k]
                   - f_2 * slp1_126[k]
                   + f_3 * pc_x[k] * sld_252[k];

        t_421[k] = f_8 * skd_204[k]
                   + f_3 * pc_y[k] * sld_252[k];

        t_422[k] = f_7 * skd_198[k]
                   + f_3 * pc_z[k] * sld_252[k];

        t_423[k] = f_3 * pc_x[k] * sld_255[k];

        t_424[k] = f_3 * pc_x[k] * sld_256[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, pc_x, pc_y, pc_z, skd_201, skd_207, \
                         skd_209, slp0_127, slp1_127, sld_255, \
                         sld_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = f_3 * pc_x[k] * sld_257[k];

        t_426[k] = f_8 * skd_207[k]
                   + f_1 * slp0_127[k]
                   - f_2 * slp1_127[k]
                   + f_3 * pc_y[k] * sld_255[k];

        t_427[k] = f_7 * skd_201[k]
                   + f_3 * pc_z[k] * sld_255[k];

        t_428[k] = f_8 * skd_209[k]
                   + f_3 * pc_y[k] * sld_257[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, t_432, pb_y, pc_y, pc_z, skf0_350, skd_203, \
                         skd_204, skd_210, skf1_350, slp0_128, slp1_128, sld_257, \
                         sld_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = f_7 * skd_203[k]
                   + f_1 * slp0_128[k]
                   - f_2 * slp1_128[k]
                   + f_3 * pc_z[k] * sld_257[k];

        t_430[k] = pb_y[k] * skf0_350[k]
                   - f_4 * pc_y[k] * skf1_350[k];

        t_431[k] = f_5 * skd_210[k]
                   + f_3 * pc_y[k] * sld_258[k];

        t_432[k] = f_6 * skd_204[k]
                   + f_3 * pc_z[k] * sld_258[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, t_436, t_437, pb_y, pc_x, pc_y, pc_z, skf0_356, \
                         skd_207, skd_213, skf1_356, sld_261, sld_262, \
                         sld_263 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = f_3 * pc_x[k] * sld_261[k];

        t_434[k] = f_3 * pc_x[k] * sld_262[k];

        t_435[k] = f_3 * pc_x[k] * sld_263[k];

        t_436[k] = pb_y[k] * skf0_356[k]
                   + f_10 * skd_213[k]
                   - f_4 * pc_y[k] * skf1_356[k];

        t_437[k] = f_6 * skd_207[k]
                   + f_3 * pc_z[k] * sld_261[k];
    }

#pragma omp simd aligned(t_438, t_439, t_440, t_441, pb_y, pc_x, pc_y, skf0_359, skd_215, \
                         skf1_359, slp0_132, slp1_132, sld_263, \
                         sld_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_438[k] = f_5 * skd_215[k]
                   + f_3 * pc_y[k] * sld_263[k];

        t_439[k] = pb_y[k] * skf0_359[k]
                   - f_4 * pc_y[k] * skf1_359[k];

        t_440[k] = f_1 * slp0_132[k]
                   - f_2 * slp1_132[k]
                   + f_3 * pc_x[k] * sld_264[k];

        t_441[k] = f_3 * pc_y[k] * sld_264[k];
    }

#pragma omp simd aligned(t_442, t_443, t_444, t_445, t_446, pc_x, pc_y, pc_z, skd_210, \
                         slp0_133, slp1_133, sld_264, sld_267, sld_268, \
                         sld_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_442[k] = f_0 * skd_210[k]
                   + f_3 * pc_z[k] * sld_264[k];

        t_443[k] = f_3 * pc_x[k] * sld_267[k];

        t_444[k] = f_3 * pc_x[k] * sld_268[k];

        t_445[k] = f_3 * pc_x[k] * sld_269[k];

        t_446[k] = f_1 * slp0_133[k]
                   - f_2 * slp1_133[k]
                   + f_3 * pc_y[k] * sld_267[k];
    }

#pragma omp simd aligned(t_447, t_448, t_449, pc_y, pc_z, skd_213, skd_215, slp0_134, \
                         slp1_134, sld_267, sld_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_447[k] = f_0 * skd_213[k]
                   + f_3 * pc_z[k] * sld_267[k];

        t_448[k] = f_3 * pc_y[k] * sld_269[k];

        t_449[k] = f_0 * skd_215[k]
                   + f_1 * slp0_134[k]
                   - f_2 * slp1_134[k]
                   + f_3 * pc_z[k] * sld_269[k];
    }
}

auto
compute_prim_slf_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t skf0, const size_t skd,
                                                   const size_t skf1, const size_t slp0,
                                                   const size_t slp1, const size_t sld,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_slf_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, skf0, skd,
                                                              skf1, slp0, slp1, sld, ncols,
                                                              gamma, p, q);

    compute_prim_slf_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, skf0, skd,
                                                              skf1, slp0, slp1, sld, ncols,
                                                              gamma, p, q);

    compute_prim_slf_three_center_electron_repulsion_0_piece2(buffer, target, pb, pc, skf0, skd,
                                                              skf1, slp0, slp1, sld, ncols,
                                                              gamma, p, q);

    compute_prim_slf_three_center_electron_repulsion_0_piece3(buffer, target, pb, pc, skf0, skd,
                                                              skf1, slp0, slp1, sld, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
