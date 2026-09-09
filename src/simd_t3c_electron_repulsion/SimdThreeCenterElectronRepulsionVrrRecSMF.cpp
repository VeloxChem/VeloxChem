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


#include "SimdThreeCenterElectronRepulsionVrrRecSMF.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_smf_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t slf0,
                                                          const size_t sld, const size_t slf1,
                                                          const size_t smp0, const size_t smp1,
                                                          const size_t smd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / q;
    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 4.0 / q;
    const auto f_7 = 3.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 3.0 / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 2.5 / q;
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
    auto *t_126 = buffer.data(target + 126);
    auto *t_127 = buffer.data(target + 127);
    auto *t_128 = buffer.data(target + 128);
    auto *t_129 = buffer.data(target + 129);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *slf0_0 = buffer.data(slf0 + 0);
    const auto *slf0_6 = buffer.data(slf0 + 6);
    const auto *slf0_9 = buffer.data(slf0 + 9);
    const auto *slf0_16 = buffer.data(slf0 + 16);
    const auto *slf0_20 = buffer.data(slf0 + 20);
    const auto *slf0_29 = buffer.data(slf0 + 29);
    const auto *slf0_30 = buffer.data(slf0 + 30);
    const auto *slf0_36 = buffer.data(slf0 + 36);
    const auto *slf0_50 = buffer.data(slf0 + 50);
    const auto *slf0_59 = buffer.data(slf0 + 59);
    const auto *slf0_60 = buffer.data(slf0 + 60);
    const auto *slf0_66 = buffer.data(slf0 + 66);

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
    const auto *sld_69 = buffer.data(sld + 69);
    const auto *sld_70 = buffer.data(sld + 70);
    const auto *sld_71 = buffer.data(sld + 71);
    const auto *sld_72 = buffer.data(sld + 72);
    const auto *sld_75 = buffer.data(sld + 75);
    const auto *sld_76 = buffer.data(sld + 76);
    const auto *sld_77 = buffer.data(sld + 77);

    const auto *slf1_0 = buffer.data(slf1 + 0);
    const auto *slf1_6 = buffer.data(slf1 + 6);
    const auto *slf1_9 = buffer.data(slf1 + 9);
    const auto *slf1_16 = buffer.data(slf1 + 16);
    const auto *slf1_20 = buffer.data(slf1 + 20);
    const auto *slf1_29 = buffer.data(slf1 + 29);
    const auto *slf1_30 = buffer.data(slf1 + 30);
    const auto *slf1_36 = buffer.data(slf1 + 36);
    const auto *slf1_50 = buffer.data(slf1 + 50);
    const auto *slf1_59 = buffer.data(slf1 + 59);
    const auto *slf1_60 = buffer.data(slf1 + 60);
    const auto *slf1_66 = buffer.data(slf1 + 66);

    const auto *smp0_0 = buffer.data(smp0 + 0);
    const auto *smp0_1 = buffer.data(smp0 + 1);
    const auto *smp0_2 = buffer.data(smp0 + 2);
    const auto *smp0_4 = buffer.data(smp0 + 4);
    const auto *smp0_8 = buffer.data(smp0 + 8);
    const auto *smp0_9 = buffer.data(smp0 + 9);
    const auto *smp0_10 = buffer.data(smp0 + 10);
    const auto *smp0_11 = buffer.data(smp0 + 11);
    const auto *smp0_15 = buffer.data(smp0 + 15);
    const auto *smp0_16 = buffer.data(smp0 + 16);
    const auto *smp0_17 = buffer.data(smp0 + 17);
    const auto *smp0_18 = buffer.data(smp0 + 18);
    const auto *smp0_19 = buffer.data(smp0 + 19);
    const auto *smp0_20 = buffer.data(smp0 + 20);
    const auto *smp0_23 = buffer.data(smp0 + 23);
    const auto *smp0_25 = buffer.data(smp0 + 25);
    const auto *smp0_27 = buffer.data(smp0 + 27);
    const auto *smp0_28 = buffer.data(smp0 + 28);
    const auto *smp0_29 = buffer.data(smp0 + 29);
    const auto *smp0_30 = buffer.data(smp0 + 30);
    const auto *smp0_31 = buffer.data(smp0 + 31);
    const auto *smp0_32 = buffer.data(smp0 + 32);
    const auto *smp0_35 = buffer.data(smp0 + 35);
    const auto *smp0_36 = buffer.data(smp0 + 36);
    const auto *smp0_37 = buffer.data(smp0 + 37);
    const auto *smp0_38 = buffer.data(smp0 + 38);

    const auto *smp1_0 = buffer.data(smp1 + 0);
    const auto *smp1_1 = buffer.data(smp1 + 1);
    const auto *smp1_2 = buffer.data(smp1 + 2);
    const auto *smp1_4 = buffer.data(smp1 + 4);
    const auto *smp1_8 = buffer.data(smp1 + 8);
    const auto *smp1_9 = buffer.data(smp1 + 9);
    const auto *smp1_10 = buffer.data(smp1 + 10);
    const auto *smp1_11 = buffer.data(smp1 + 11);
    const auto *smp1_15 = buffer.data(smp1 + 15);
    const auto *smp1_16 = buffer.data(smp1 + 16);
    const auto *smp1_17 = buffer.data(smp1 + 17);
    const auto *smp1_18 = buffer.data(smp1 + 18);
    const auto *smp1_19 = buffer.data(smp1 + 19);
    const auto *smp1_20 = buffer.data(smp1 + 20);
    const auto *smp1_23 = buffer.data(smp1 + 23);
    const auto *smp1_25 = buffer.data(smp1 + 25);
    const auto *smp1_27 = buffer.data(smp1 + 27);
    const auto *smp1_28 = buffer.data(smp1 + 28);
    const auto *smp1_29 = buffer.data(smp1 + 29);
    const auto *smp1_30 = buffer.data(smp1 + 30);
    const auto *smp1_31 = buffer.data(smp1 + 31);
    const auto *smp1_32 = buffer.data(smp1 + 32);
    const auto *smp1_35 = buffer.data(smp1 + 35);
    const auto *smp1_36 = buffer.data(smp1 + 36);
    const auto *smp1_37 = buffer.data(smp1 + 37);
    const auto *smp1_38 = buffer.data(smp1 + 38);

    const auto *smd_0 = buffer.data(smd + 0);
    const auto *smd_3 = buffer.data(smd + 3);
    const auto *smd_4 = buffer.data(smd + 4);
    const auto *smd_5 = buffer.data(smd + 5);
    const auto *smd_6 = buffer.data(smd + 6);
    const auto *smd_9 = buffer.data(smd + 9);
    const auto *smd_10 = buffer.data(smd + 10);
    const auto *smd_11 = buffer.data(smd + 11);
    const auto *smd_12 = buffer.data(smd + 12);
    const auto *smd_15 = buffer.data(smd + 15);
    const auto *smd_16 = buffer.data(smd + 16);
    const auto *smd_17 = buffer.data(smd + 17);
    const auto *smd_18 = buffer.data(smd + 18);
    const auto *smd_21 = buffer.data(smd + 21);
    const auto *smd_22 = buffer.data(smd + 22);
    const auto *smd_23 = buffer.data(smd + 23);
    const auto *smd_24 = buffer.data(smd + 24);
    const auto *smd_27 = buffer.data(smd + 27);
    const auto *smd_28 = buffer.data(smd + 28);
    const auto *smd_29 = buffer.data(smd + 29);
    const auto *smd_30 = buffer.data(smd + 30);
    const auto *smd_33 = buffer.data(smd + 33);
    const auto *smd_34 = buffer.data(smd + 34);
    const auto *smd_35 = buffer.data(smd + 35);
    const auto *smd_36 = buffer.data(smd + 36);
    const auto *smd_39 = buffer.data(smd + 39);
    const auto *smd_40 = buffer.data(smd + 40);
    const auto *smd_41 = buffer.data(smd + 41);
    const auto *smd_42 = buffer.data(smd + 42);
    const auto *smd_45 = buffer.data(smd + 45);
    const auto *smd_46 = buffer.data(smd + 46);
    const auto *smd_47 = buffer.data(smd + 47);
    const auto *smd_48 = buffer.data(smd + 48);
    const auto *smd_51 = buffer.data(smd + 51);
    const auto *smd_52 = buffer.data(smd + 52);
    const auto *smd_53 = buffer.data(smd + 53);
    const auto *smd_54 = buffer.data(smd + 54);
    const auto *smd_57 = buffer.data(smd + 57);
    const auto *smd_58 = buffer.data(smd + 58);
    const auto *smd_59 = buffer.data(smd + 59);
    const auto *smd_60 = buffer.data(smd + 60);
    const auto *smd_63 = buffer.data(smd + 63);
    const auto *smd_64 = buffer.data(smd + 64);
    const auto *smd_65 = buffer.data(smd + 65);
    const auto *smd_66 = buffer.data(smd + 66);
    const auto *smd_69 = buffer.data(smd + 69);
    const auto *smd_70 = buffer.data(smd + 70);
    const auto *smd_71 = buffer.data(smd + 71);
    const auto *smd_72 = buffer.data(smd + 72);
    const auto *smd_75 = buffer.data(smd + 75);
    const auto *smd_76 = buffer.data(smd + 76);
    const auto *smd_77 = buffer.data(smd + 77);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pc_x, pc_y, pc_z, sld_0, sld_3, sld_4, \
                         smp0_0, smp1_0, smd_0, smd_3, smd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sld_0[k]
                 + f_1 * smp0_0[k]
                 - f_2 * smp1_0[k]
                 + f_3 * pc_x[k] * smd_0[k];

        t_1[k] = f_3 * pc_y[k] * smd_0[k];

        t_2[k] = f_3 * pc_z[k] * smd_0[k];

        t_3[k] = f_0 * sld_3[k]
                 + f_3 * pc_x[k] * smd_3[k];

        t_4[k] = f_0 * sld_4[k]
                 + f_3 * pc_x[k] * smd_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pc_x, pc_y, pc_z, sld_5, smp0_1, smp0_2, \
                         smp1_1, smp1_2, smd_3, smd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * sld_5[k]
                 + f_3 * pc_x[k] * smd_5[k];

        t_6[k] = f_1 * smp0_1[k]
                 - f_2 * smp1_1[k]
                 + f_3 * pc_y[k] * smd_3[k];

        t_7[k] = f_3 * pc_z[k] * smd_3[k];

        t_8[k] = f_3 * pc_y[k] * smd_5[k];

        t_9[k] = f_1 * smp0_2[k]
                 - f_2 * smp1_2[k]
                 + f_3 * pc_z[k] * smd_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pb_y, pc_x, pc_y, pc_z, slf0_0, sld_0, sld_9, \
                         slf1_0, smd_6, smd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pb_y[k] * slf0_0[k]
                  - f_4 * pc_y[k] * slf1_0[k];

        t_11[k] = f_5 * sld_0[k]
                  + f_3 * pc_y[k] * smd_6[k];

        t_12[k] = f_3 * pc_z[k] * smd_6[k];

        t_13[k] = f_6 * sld_9[k]
                  + f_3 * pc_x[k] * smd_9[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pc_x, pc_y, pc_z, sld_3, sld_10, sld_11, \
                         smp0_4, smp1_4, smd_9, smd_10, smd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_6 * sld_10[k]
                  + f_3 * pc_x[k] * smd_10[k];

        t_15[k] = f_6 * sld_11[k]
                  + f_3 * pc_x[k] * smd_11[k];

        t_16[k] = f_5 * sld_3[k]
                  + f_1 * smp0_4[k]
                  - f_2 * smp1_4[k]
                  + f_3 * pc_y[k] * smd_9[k];

        t_17[k] = f_3 * pc_z[k] * smd_9[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pb_y, pb_z, pc_y, pc_z, slf0_0, slf0_9, \
                         sld_5, slf1_0, slf1_9, smd_11, smd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * sld_5[k]
                  + f_3 * pc_y[k] * smd_11[k];

        t_19[k] = pb_y[k] * slf0_9[k]
                  - f_4 * pc_y[k] * slf1_9[k];

        t_20[k] = pb_z[k] * slf0_0[k]
                  - f_4 * pc_z[k] * slf1_0[k];

        t_21[k] = f_3 * pc_y[k] * smd_12[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pc_x, pc_z, sld_0, sld_15, sld_16, sld_17, \
                         smd_12, smd_15, smd_16, smd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_5 * sld_0[k]
                  + f_3 * pc_z[k] * smd_12[k];

        t_23[k] = f_6 * sld_15[k]
                  + f_3 * pc_x[k] * smd_15[k];

        t_24[k] = f_6 * sld_16[k]
                  + f_3 * pc_x[k] * smd_16[k];

        t_25[k] = f_6 * sld_17[k]
                  + f_3 * pc_x[k] * smd_17[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pb_z, pc_y, pc_z, slf0_6, sld_3, sld_5, \
                         slf1_6, smp0_8, smp1_8, smd_15, smd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pb_z[k] * slf0_6[k]
                  - f_4 * pc_z[k] * slf1_6[k];

        t_27[k] = f_5 * sld_3[k]
                  + f_3 * pc_z[k] * smd_15[k];

        t_28[k] = f_3 * pc_y[k] * smd_17[k];

        t_29[k] = f_5 * sld_5[k]
                  + f_1 * smp0_8[k]
                  - f_2 * smp1_8[k]
                  + f_3 * pc_z[k] * smd_17[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pc_x, pc_y, pc_z, sld_6, sld_18, sld_21, \
                         smp0_9, smp1_9, smd_18, smd_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_7 * sld_18[k]
                  + f_1 * smp0_9[k]
                  - f_2 * smp1_9[k]
                  + f_3 * pc_x[k] * smd_18[k];

        t_31[k] = f_8 * sld_6[k]
                  + f_3 * pc_y[k] * smd_18[k];

        t_32[k] = f_3 * pc_z[k] * smd_18[k];

        t_33[k] = f_7 * sld_21[k]
                  + f_3 * pc_x[k] * smd_21[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pc_x, pc_y, pc_z, sld_9, sld_22, sld_23, \
                         smp0_10, smp1_10, smd_21, smd_22, smd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_7 * sld_22[k]
                  + f_3 * pc_x[k] * smd_22[k];

        t_35[k] = f_7 * sld_23[k]
                  + f_3 * pc_x[k] * smd_23[k];

        t_36[k] = f_8 * sld_9[k]
                  + f_1 * smp0_10[k]
                  - f_2 * smp1_10[k]
                  + f_3 * pc_y[k] * smd_21[k];

        t_37[k] = f_3 * pc_z[k] * smd_21[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pb_y, pc_y, pc_z, slf0_20, sld_11, sld_12, \
                         slf1_20, smp0_11, smp1_11, smd_23, smd_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_8 * sld_11[k]
                  + f_3 * pc_y[k] * smd_23[k];

        t_39[k] = f_1 * smp0_11[k]
                  - f_2 * smp1_11[k]
                  + f_3 * pc_z[k] * smd_23[k];

        t_40[k] = pb_y[k] * slf0_20[k]
                  - f_4 * pc_y[k] * slf1_20[k];

        t_41[k] = f_5 * sld_12[k]
                  + f_3 * pc_y[k] * smd_24[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pc_x, pc_z, sld_6, sld_27, sld_28, sld_29, \
                         smd_24, smd_27, smd_28, smd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_5 * sld_6[k]
                  + f_3 * pc_z[k] * smd_24[k];

        t_43[k] = f_7 * sld_27[k]
                  + f_3 * pc_x[k] * smd_27[k];

        t_44[k] = f_7 * sld_28[k]
                  + f_3 * pc_x[k] * smd_28[k];

        t_45[k] = f_7 * sld_29[k]
                  + f_3 * pc_x[k] * smd_29[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pb_y, pb_z, pc_y, pc_z, slf0_16, slf0_29, \
                         sld_9, sld_17, slf1_16, slf1_29, smd_27, \
                         smd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = pb_z[k] * slf0_16[k]
                  - f_4 * pc_z[k] * slf1_16[k];

        t_47[k] = f_5 * sld_9[k]
                  + f_3 * pc_z[k] * smd_27[k];

        t_48[k] = f_5 * sld_17[k]
                  + f_3 * pc_y[k] * smd_29[k];

        t_49[k] = pb_y[k] * slf0_29[k]
                  - f_4 * pc_y[k] * slf1_29[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pc_x, pc_y, pc_z, sld_12, sld_30, sld_33, \
                         smp0_15, smp1_15, smd_30, smd_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_7 * sld_30[k]
                  + f_1 * smp0_15[k]
                  - f_2 * smp1_15[k]
                  + f_3 * pc_x[k] * smd_30[k];

        t_51[k] = f_3 * pc_y[k] * smd_30[k];

        t_52[k] = f_8 * sld_12[k]
                  + f_3 * pc_z[k] * smd_30[k];

        t_53[k] = f_7 * sld_33[k]
                  + f_3 * pc_x[k] * smd_33[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, pc_x, pc_y, pc_z, sld_15, sld_34, \
                         sld_35, smp0_16, smp1_16, smd_33, smd_34, \
                         smd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_7 * sld_34[k]
                  + f_3 * pc_x[k] * smd_34[k];

        t_55[k] = f_7 * sld_35[k]
                  + f_3 * pc_x[k] * smd_35[k];

        t_56[k] = f_1 * smp0_16[k]
                  - f_2 * smp1_16[k]
                  + f_3 * pc_y[k] * smd_33[k];

        t_57[k] = f_8 * sld_15[k]
                  + f_3 * pc_z[k] * smd_33[k];

        t_58[k] = f_3 * pc_y[k] * smd_35[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pc_x, pc_y, pc_z, sld_17, sld_18, sld_36, \
                         smp0_17, smp0_18, smp1_17, smp1_18, smd_35, \
                         smd_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_8 * sld_17[k]
                  + f_1 * smp0_17[k]
                  - f_2 * smp1_17[k]
                  + f_3 * pc_z[k] * smd_35[k];

        t_60[k] = f_9 * sld_36[k]
                  + f_1 * smp0_18[k]
                  - f_2 * smp1_18[k]
                  + f_3 * pc_x[k] * smd_36[k];

        t_61[k] = f_10 * sld_18[k]
                  + f_3 * pc_y[k] * smd_36[k];

        t_62[k] = f_3 * pc_z[k] * smd_36[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pc_x, pc_y, sld_21, sld_39, sld_40, sld_41, \
                         smp0_19, smp1_19, smd_39, smd_40, smd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_9 * sld_39[k]
                  + f_3 * pc_x[k] * smd_39[k];

        t_64[k] = f_9 * sld_40[k]
                  + f_3 * pc_x[k] * smd_40[k];

        t_65[k] = f_9 * sld_41[k]
                  + f_3 * pc_x[k] * smd_41[k];

        t_66[k] = f_10 * sld_21[k]
                  + f_1 * smp0_19[k]
                  - f_2 * smp1_19[k]
                  + f_3 * pc_y[k] * smd_39[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, pb_z, pc_y, pc_z, slf0_30, sld_23, slf1_30, \
                         smp0_20, smp1_20, smd_39, smd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_3 * pc_z[k] * smd_39[k];

        t_68[k] = f_10 * sld_23[k]
                  + f_3 * pc_y[k] * smd_41[k];

        t_69[k] = f_1 * smp0_20[k]
                  - f_2 * smp1_20[k]
                  + f_3 * pc_z[k] * smd_41[k];

        t_70[k] = pb_z[k] * slf0_30[k]
                  - f_4 * pc_z[k] * slf1_30[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pc_x, pc_y, pc_z, sld_18, sld_24, sld_45, \
                         sld_46, smd_42, smd_45, smd_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_8 * sld_24[k]
                  + f_3 * pc_y[k] * smd_42[k];

        t_72[k] = f_5 * sld_18[k]
                  + f_3 * pc_z[k] * smd_42[k];

        t_73[k] = f_9 * sld_45[k]
                  + f_3 * pc_x[k] * smd_45[k];

        t_74[k] = f_9 * sld_46[k]
                  + f_3 * pc_x[k] * smd_46[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pb_z, pc_x, pc_y, pc_z, slf0_36, sld_21, \
                         sld_29, sld_47, slf1_36, smd_45, smd_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_9 * sld_47[k]
                  + f_3 * pc_x[k] * smd_47[k];

        t_76[k] = pb_z[k] * slf0_36[k]
                  - f_4 * pc_z[k] * slf1_36[k];

        t_77[k] = f_5 * sld_21[k]
                  + f_3 * pc_z[k] * smd_45[k];

        t_78[k] = f_8 * sld_29[k]
                  + f_3 * pc_y[k] * smd_47[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pb_y, pc_y, pc_z, slf0_50, sld_23, sld_24, \
                         sld_30, slf1_50, smp0_23, smp1_23, smd_47, \
                         smd_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_5 * sld_23[k]
                  + f_1 * smp0_23[k]
                  - f_2 * smp1_23[k]
                  + f_3 * pc_z[k] * smd_47[k];

        t_80[k] = pb_y[k] * slf0_50[k]
                  - f_4 * pc_y[k] * slf1_50[k];

        t_81[k] = f_5 * sld_30[k]
                  + f_3 * pc_y[k] * smd_48[k];

        t_82[k] = f_8 * sld_24[k]
                  + f_3 * pc_z[k] * smd_48[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pc_x, pc_y, sld_33, sld_51, sld_52, sld_53, \
                         smp0_25, smp1_25, smd_51, smd_52, smd_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_9 * sld_51[k]
                  + f_3 * pc_x[k] * smd_51[k];

        t_84[k] = f_9 * sld_52[k]
                  + f_3 * pc_x[k] * smd_52[k];

        t_85[k] = f_9 * sld_53[k]
                  + f_3 * pc_x[k] * smd_53[k];

        t_86[k] = f_5 * sld_33[k]
                  + f_1 * smp0_25[k]
                  - f_2 * smp1_25[k]
                  + f_3 * pc_y[k] * smd_51[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, pb_y, pc_y, pc_z, slf0_59, sld_27, sld_35, slf1_59, \
                         smd_51, smd_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_8 * sld_27[k]
                  + f_3 * pc_z[k] * smd_51[k];

        t_88[k] = f_5 * sld_35[k]
                  + f_3 * pc_y[k] * smd_53[k];

        t_89[k] = pb_y[k] * slf0_59[k]
                  - f_4 * pc_y[k] * slf1_59[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, pc_x, pc_y, pc_z, sld_30, sld_54, sld_57, \
                         smp0_27, smp1_27, smd_54, smd_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_9 * sld_54[k]
                  + f_1 * smp0_27[k]
                  - f_2 * smp1_27[k]
                  + f_3 * pc_x[k] * smd_54[k];

        t_91[k] = f_3 * pc_y[k] * smd_54[k];

        t_92[k] = f_10 * sld_30[k]
                  + f_3 * pc_z[k] * smd_54[k];

        t_93[k] = f_9 * sld_57[k]
                  + f_3 * pc_x[k] * smd_57[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, pc_x, pc_y, pc_z, sld_33, sld_58, \
                         sld_59, smp0_28, smp1_28, smd_57, smd_58, \
                         smd_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_9 * sld_58[k]
                  + f_3 * pc_x[k] * smd_58[k];

        t_95[k] = f_9 * sld_59[k]
                  + f_3 * pc_x[k] * smd_59[k];

        t_96[k] = f_1 * smp0_28[k]
                  - f_2 * smp1_28[k]
                  + f_3 * pc_y[k] * smd_57[k];

        t_97[k] = f_10 * sld_33[k]
                  + f_3 * pc_z[k] * smd_57[k];

        t_98[k] = f_3 * pc_y[k] * smd_59[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, pc_x, pc_y, pc_z, sld_35, sld_36, sld_60, \
                         smp0_29, smp0_30, smp1_29, smp1_30, smd_59, \
                         smd_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_10 * sld_35[k]
                  + f_1 * smp0_29[k]
                  - f_2 * smp1_29[k]
                  + f_3 * pc_z[k] * smd_59[k];

        t_100[k] = f_11 * sld_60[k]
                   + f_1 * smp0_30[k]
                   - f_2 * smp1_30[k]
                   + f_3 * pc_x[k] * smd_60[k];

        t_101[k] = f_12 * sld_36[k]
                   + f_3 * pc_y[k] * smd_60[k];

        t_102[k] = f_3 * pc_z[k] * smd_60[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, pc_x, pc_y, sld_39, sld_63, sld_64, \
                         sld_65, smp0_31, smp1_31, smd_63, smd_64, \
                         smd_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_11 * sld_63[k]
                   + f_3 * pc_x[k] * smd_63[k];

        t_104[k] = f_11 * sld_64[k]
                   + f_3 * pc_x[k] * smd_64[k];

        t_105[k] = f_11 * sld_65[k]
                   + f_3 * pc_x[k] * smd_65[k];

        t_106[k] = f_12 * sld_39[k]
                   + f_1 * smp0_31[k]
                   - f_2 * smp1_31[k]
                   + f_3 * pc_y[k] * smd_63[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, pb_z, pc_y, pc_z, slf0_60, sld_41, \
                         slf1_60, smp0_32, smp1_32, smd_63, smd_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_3 * pc_z[k] * smd_63[k];

        t_108[k] = f_12 * sld_41[k]
                   + f_3 * pc_y[k] * smd_65[k];

        t_109[k] = f_1 * smp0_32[k]
                   - f_2 * smp1_32[k]
                   + f_3 * pc_z[k] * smd_65[k];

        t_110[k] = pb_z[k] * slf0_60[k]
                   - f_4 * pc_z[k] * slf1_60[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, pc_x, pc_y, pc_z, sld_36, sld_42, sld_69, \
                         sld_70, smd_66, smd_69, smd_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_10 * sld_42[k]
                   + f_3 * pc_y[k] * smd_66[k];

        t_112[k] = f_5 * sld_36[k]
                   + f_3 * pc_z[k] * smd_66[k];

        t_113[k] = f_11 * sld_69[k]
                   + f_3 * pc_x[k] * smd_69[k];

        t_114[k] = f_11 * sld_70[k]
                   + f_3 * pc_x[k] * smd_70[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, pb_z, pc_x, pc_y, pc_z, slf0_66, sld_39, \
                         sld_47, sld_71, slf1_66, smd_69, smd_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_11 * sld_71[k]
                   + f_3 * pc_x[k] * smd_71[k];

        t_116[k] = pb_z[k] * slf0_66[k]
                   - f_4 * pc_z[k] * slf1_66[k];

        t_117[k] = f_5 * sld_39[k]
                   + f_3 * pc_z[k] * smd_69[k];

        t_118[k] = f_10 * sld_47[k]
                   + f_3 * pc_y[k] * smd_71[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, pc_x, pc_y, pc_z, sld_41, sld_48, sld_72, \
                         smp0_35, smp0_36, smp1_35, smp1_36, smd_71, \
                         smd_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = f_5 * sld_41[k]
                   + f_1 * smp0_35[k]
                   - f_2 * smp1_35[k]
                   + f_3 * pc_z[k] * smd_71[k];

        t_120[k] = f_11 * sld_72[k]
                   + f_1 * smp0_36[k]
                   - f_2 * smp1_36[k]
                   + f_3 * pc_x[k] * smd_72[k];

        t_121[k] = f_8 * sld_48[k]
                   + f_3 * pc_y[k] * smd_72[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, pc_x, pc_z, sld_42, sld_75, sld_76, \
                         sld_77, smd_72, smd_75, smd_76, smd_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_8 * sld_42[k]
                   + f_3 * pc_z[k] * smd_72[k];

        t_123[k] = f_11 * sld_75[k]
                   + f_3 * pc_x[k] * smd_75[k];

        t_124[k] = f_11 * sld_76[k]
                   + f_3 * pc_x[k] * smd_76[k];

        t_125[k] = f_11 * sld_77[k]
                   + f_3 * pc_x[k] * smd_77[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, pc_y, pc_z, sld_45, sld_47, sld_51, \
                         sld_53, smp0_37, smp0_38, smp1_37, smp1_38, smd_75, \
                         smd_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_8 * sld_51[k]
                   + f_1 * smp0_37[k]
                   - f_2 * smp1_37[k]
                   + f_3 * pc_y[k] * smd_75[k];

        t_127[k] = f_8 * sld_45[k]
                   + f_3 * pc_z[k] * smd_75[k];

        t_128[k] = f_8 * sld_53[k]
                   + f_3 * pc_y[k] * smd_77[k];

        t_129[k] = f_8 * sld_47[k]
                   + f_1 * smp0_38[k]
                   - f_2 * smp1_38[k]
                   + f_3 * pc_z[k] * smd_77[k];
    }
}

static auto
compute_prim_smf_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t slf0,
                                                          const size_t sld, const size_t slf1,
                                                          const size_t smp0, const size_t smp1,
                                                          const size_t smd, const size_t ncols,
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
    const auto f_8 = 1.0 / q;
    const auto f_9 = 3.0 / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 2.5 / q;
    const auto f_12 = 2.0 / q;

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

    const auto *slf0_90 = buffer.data(slf0 + 90);
    const auto *slf0_99 = buffer.data(slf0 + 99);
    const auto *slf0_100 = buffer.data(slf0 + 100);
    const auto *slf0_106 = buffer.data(slf0 + 106);
    const auto *slf0_140 = buffer.data(slf0 + 140);
    const auto *slf0_149 = buffer.data(slf0 + 149);
    const auto *slf0_150 = buffer.data(slf0 + 150);
    const auto *slf0_156 = buffer.data(slf0 + 156);

    const auto *sld_48 = buffer.data(sld + 48);
    const auto *sld_51 = buffer.data(sld + 51);
    const auto *sld_54 = buffer.data(sld + 54);
    const auto *sld_57 = buffer.data(sld + 57);
    const auto *sld_59 = buffer.data(sld + 59);
    const auto *sld_60 = buffer.data(sld + 60);
    const auto *sld_63 = buffer.data(sld + 63);
    const auto *sld_65 = buffer.data(sld + 65);
    const auto *sld_66 = buffer.data(sld + 66);
    const auto *sld_69 = buffer.data(sld + 69);
    const auto *sld_71 = buffer.data(sld + 71);
    const auto *sld_72 = buffer.data(sld + 72);
    const auto *sld_75 = buffer.data(sld + 75);
    const auto *sld_77 = buffer.data(sld + 77);
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

    const auto *slf1_90 = buffer.data(slf1 + 90);
    const auto *slf1_99 = buffer.data(slf1 + 99);
    const auto *slf1_100 = buffer.data(slf1 + 100);
    const auto *slf1_106 = buffer.data(slf1 + 106);
    const auto *slf1_140 = buffer.data(slf1 + 140);
    const auto *slf1_149 = buffer.data(slf1 + 149);
    const auto *slf1_150 = buffer.data(slf1 + 150);
    const auto *slf1_156 = buffer.data(slf1 + 156);

    const auto *smp0_40 = buffer.data(smp0 + 40);
    const auto *smp0_42 = buffer.data(smp0 + 42);
    const auto *smp0_43 = buffer.data(smp0 + 43);
    const auto *smp0_44 = buffer.data(smp0 + 44);
    const auto *smp0_45 = buffer.data(smp0 + 45);
    const auto *smp0_46 = buffer.data(smp0 + 46);
    const auto *smp0_47 = buffer.data(smp0 + 47);
    const auto *smp0_50 = buffer.data(smp0 + 50);
    const auto *smp0_51 = buffer.data(smp0 + 51);
    const auto *smp0_52 = buffer.data(smp0 + 52);
    const auto *smp0_53 = buffer.data(smp0 + 53);
    const auto *smp0_54 = buffer.data(smp0 + 54);
    const auto *smp0_55 = buffer.data(smp0 + 55);
    const auto *smp0_56 = buffer.data(smp0 + 56);
    const auto *smp0_58 = buffer.data(smp0 + 58);
    const auto *smp0_60 = buffer.data(smp0 + 60);
    const auto *smp0_61 = buffer.data(smp0 + 61);
    const auto *smp0_62 = buffer.data(smp0 + 62);
    const auto *smp0_63 = buffer.data(smp0 + 63);
    const auto *smp0_64 = buffer.data(smp0 + 64);
    const auto *smp0_65 = buffer.data(smp0 + 65);
    const auto *smp0_68 = buffer.data(smp0 + 68);
    const auto *smp0_69 = buffer.data(smp0 + 69);
    const auto *smp0_70 = buffer.data(smp0 + 70);
    const auto *smp0_71 = buffer.data(smp0 + 71);
    const auto *smp0_72 = buffer.data(smp0 + 72);
    const auto *smp0_73 = buffer.data(smp0 + 73);
    const auto *smp0_74 = buffer.data(smp0 + 74);
    const auto *smp0_75 = buffer.data(smp0 + 75);

    const auto *smp1_40 = buffer.data(smp1 + 40);
    const auto *smp1_42 = buffer.data(smp1 + 42);
    const auto *smp1_43 = buffer.data(smp1 + 43);
    const auto *smp1_44 = buffer.data(smp1 + 44);
    const auto *smp1_45 = buffer.data(smp1 + 45);
    const auto *smp1_46 = buffer.data(smp1 + 46);
    const auto *smp1_47 = buffer.data(smp1 + 47);
    const auto *smp1_50 = buffer.data(smp1 + 50);
    const auto *smp1_51 = buffer.data(smp1 + 51);
    const auto *smp1_52 = buffer.data(smp1 + 52);
    const auto *smp1_53 = buffer.data(smp1 + 53);
    const auto *smp1_54 = buffer.data(smp1 + 54);
    const auto *smp1_55 = buffer.data(smp1 + 55);
    const auto *smp1_56 = buffer.data(smp1 + 56);
    const auto *smp1_58 = buffer.data(smp1 + 58);
    const auto *smp1_60 = buffer.data(smp1 + 60);
    const auto *smp1_61 = buffer.data(smp1 + 61);
    const auto *smp1_62 = buffer.data(smp1 + 62);
    const auto *smp1_63 = buffer.data(smp1 + 63);
    const auto *smp1_64 = buffer.data(smp1 + 64);
    const auto *smp1_65 = buffer.data(smp1 + 65);
    const auto *smp1_68 = buffer.data(smp1 + 68);
    const auto *smp1_69 = buffer.data(smp1 + 69);
    const auto *smp1_70 = buffer.data(smp1 + 70);
    const auto *smp1_71 = buffer.data(smp1 + 71);
    const auto *smp1_72 = buffer.data(smp1 + 72);
    const auto *smp1_73 = buffer.data(smp1 + 73);
    const auto *smp1_74 = buffer.data(smp1 + 74);
    const auto *smp1_75 = buffer.data(smp1 + 75);

    const auto *smd_78 = buffer.data(smd + 78);
    const auto *smd_81 = buffer.data(smd + 81);
    const auto *smd_82 = buffer.data(smd + 82);
    const auto *smd_83 = buffer.data(smd + 83);
    const auto *smd_84 = buffer.data(smd + 84);
    const auto *smd_87 = buffer.data(smd + 87);
    const auto *smd_88 = buffer.data(smd + 88);
    const auto *smd_89 = buffer.data(smd + 89);
    const auto *smd_90 = buffer.data(smd + 90);
    const auto *smd_93 = buffer.data(smd + 93);
    const auto *smd_94 = buffer.data(smd + 94);
    const auto *smd_95 = buffer.data(smd + 95);
    const auto *smd_96 = buffer.data(smd + 96);
    const auto *smd_99 = buffer.data(smd + 99);
    const auto *smd_100 = buffer.data(smd + 100);
    const auto *smd_101 = buffer.data(smd + 101);
    const auto *smd_102 = buffer.data(smd + 102);
    const auto *smd_105 = buffer.data(smd + 105);
    const auto *smd_106 = buffer.data(smd + 106);
    const auto *smd_107 = buffer.data(smd + 107);
    const auto *smd_108 = buffer.data(smd + 108);
    const auto *smd_111 = buffer.data(smd + 111);
    const auto *smd_112 = buffer.data(smd + 112);
    const auto *smd_113 = buffer.data(smd + 113);
    const auto *smd_114 = buffer.data(smd + 114);
    const auto *smd_117 = buffer.data(smd + 117);
    const auto *smd_118 = buffer.data(smd + 118);
    const auto *smd_119 = buffer.data(smd + 119);
    const auto *smd_120 = buffer.data(smd + 120);
    const auto *smd_123 = buffer.data(smd + 123);
    const auto *smd_124 = buffer.data(smd + 124);
    const auto *smd_125 = buffer.data(smd + 125);
    const auto *smd_126 = buffer.data(smd + 126);
    const auto *smd_129 = buffer.data(smd + 129);
    const auto *smd_130 = buffer.data(smd + 130);
    const auto *smd_131 = buffer.data(smd + 131);
    const auto *smd_132 = buffer.data(smd + 132);
    const auto *smd_135 = buffer.data(smd + 135);
    const auto *smd_136 = buffer.data(smd + 136);
    const auto *smd_137 = buffer.data(smd + 137);
    const auto *smd_138 = buffer.data(smd + 138);
    const auto *smd_141 = buffer.data(smd + 141);
    const auto *smd_142 = buffer.data(smd + 142);
    const auto *smd_143 = buffer.data(smd + 143);
    const auto *smd_144 = buffer.data(smd + 144);
    const auto *smd_147 = buffer.data(smd + 147);
    const auto *smd_148 = buffer.data(smd + 148);
    const auto *smd_149 = buffer.data(smd + 149);
    const auto *smd_150 = buffer.data(smd + 150);
    const auto *smd_153 = buffer.data(smd + 153);
    const auto *smd_154 = buffer.data(smd + 154);

#pragma omp simd aligned(t_130, t_131, t_132, t_133, pb_y, pc_x, pc_y, pc_z, slf0_90, sld_48, \
                         sld_54, sld_81, slf1_90, smd_78, smd_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = pb_y[k] * slf0_90[k]
                   - f_4 * pc_y[k] * slf1_90[k];

        t_131[k] = f_5 * sld_54[k]
                   + f_3 * pc_y[k] * smd_78[k];

        t_132[k] = f_10 * sld_48[k]
                   + f_3 * pc_z[k] * smd_78[k];

        t_133[k] = f_11 * sld_81[k]
                   + f_3 * pc_x[k] * smd_81[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, pc_x, pc_y, pc_z, sld_51, sld_57, sld_82, \
                         sld_83, smp0_40, smp1_40, smd_81, smd_82, \
                         smd_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_11 * sld_82[k]
                   + f_3 * pc_x[k] * smd_82[k];

        t_135[k] = f_11 * sld_83[k]
                   + f_3 * pc_x[k] * smd_83[k];

        t_136[k] = f_5 * sld_57[k]
                   + f_1 * smp0_40[k]
                   - f_2 * smp1_40[k]
                   + f_3 * pc_y[k] * smd_81[k];

        t_137[k] = f_10 * sld_51[k]
                   + f_3 * pc_z[k] * smd_81[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, pb_y, pc_x, pc_y, slf0_99, sld_59, \
                         sld_84, slf1_99, smp0_42, smp1_42, smd_83, \
                         smd_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_5 * sld_59[k]
                   + f_3 * pc_y[k] * smd_83[k];

        t_139[k] = pb_y[k] * slf0_99[k]
                   - f_4 * pc_y[k] * slf1_99[k];

        t_140[k] = f_11 * sld_84[k]
                   + f_1 * smp0_42[k]
                   - f_2 * smp1_42[k]
                   + f_3 * pc_x[k] * smd_84[k];

        t_141[k] = f_3 * pc_y[k] * smd_84[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, pc_x, pc_z, sld_54, sld_87, sld_88, \
                         sld_89, smd_84, smd_87, smd_88, smd_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_12 * sld_54[k]
                   + f_3 * pc_z[k] * smd_84[k];

        t_143[k] = f_11 * sld_87[k]
                   + f_3 * pc_x[k] * smd_87[k];

        t_144[k] = f_11 * sld_88[k]
                   + f_3 * pc_x[k] * smd_88[k];

        t_145[k] = f_11 * sld_89[k]
                   + f_3 * pc_x[k] * smd_89[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pc_y, pc_z, sld_57, sld_59, smp0_43, \
                         smp0_44, smp1_43, smp1_44, smd_87, smd_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_1 * smp0_43[k]
                   - f_2 * smp1_43[k]
                   + f_3 * pc_y[k] * smd_87[k];

        t_147[k] = f_12 * sld_57[k]
                   + f_3 * pc_z[k] * smd_87[k];

        t_148[k] = f_3 * pc_y[k] * smd_89[k];

        t_149[k] = f_12 * sld_59[k]
                   + f_1 * smp0_44[k]
                   - f_2 * smp1_44[k]
                   + f_3 * pc_z[k] * smd_89[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, pc_x, pc_y, pc_z, sld_60, sld_90, sld_93, \
                         smp0_45, smp1_45, smd_90, smd_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_12 * sld_90[k]
                   + f_1 * smp0_45[k]
                   - f_2 * smp1_45[k]
                   + f_3 * pc_x[k] * smd_90[k];

        t_151[k] = f_11 * sld_60[k]
                   + f_3 * pc_y[k] * smd_90[k];

        t_152[k] = f_3 * pc_z[k] * smd_90[k];

        t_153[k] = f_12 * sld_93[k]
                   + f_3 * pc_x[k] * smd_93[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, pc_x, pc_y, pc_z, sld_63, sld_94, sld_95, \
                         smp0_46, smp1_46, smd_93, smd_94, smd_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = f_12 * sld_94[k]
                   + f_3 * pc_x[k] * smd_94[k];

        t_155[k] = f_12 * sld_95[k]
                   + f_3 * pc_x[k] * smd_95[k];

        t_156[k] = f_11 * sld_63[k]
                   + f_1 * smp0_46[k]
                   - f_2 * smp1_46[k]
                   + f_3 * pc_y[k] * smd_93[k];

        t_157[k] = f_3 * pc_z[k] * smd_93[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, pb_z, pc_y, pc_z, slf0_100, sld_65, \
                         sld_66, slf1_100, smp0_47, smp1_47, smd_95, \
                         smd_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_11 * sld_65[k]
                   + f_3 * pc_y[k] * smd_95[k];

        t_159[k] = f_1 * smp0_47[k]
                   - f_2 * smp1_47[k]
                   + f_3 * pc_z[k] * smd_95[k];

        t_160[k] = pb_z[k] * slf0_100[k]
                   - f_4 * pc_z[k] * slf1_100[k];

        t_161[k] = f_12 * sld_66[k]
                   + f_3 * pc_y[k] * smd_96[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, pc_x, pc_z, sld_60, sld_99, sld_100, \
                         sld_101, smd_96, smd_99, smd_100, smd_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_5 * sld_60[k]
                   + f_3 * pc_z[k] * smd_96[k];

        t_163[k] = f_12 * sld_99[k]
                   + f_3 * pc_x[k] * smd_99[k];

        t_164[k] = f_12 * sld_100[k]
                   + f_3 * pc_x[k] * smd_100[k];

        t_165[k] = f_12 * sld_101[k]
                   + f_3 * pc_x[k] * smd_101[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, pb_z, pc_y, pc_z, slf0_106, sld_63, \
                         sld_65, sld_71, slf1_106, smp0_50, smp1_50, smd_99, \
                         smd_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = pb_z[k] * slf0_106[k]
                   - f_4 * pc_z[k] * slf1_106[k];

        t_167[k] = f_5 * sld_63[k]
                   + f_3 * pc_z[k] * smd_99[k];

        t_168[k] = f_12 * sld_71[k]
                   + f_3 * pc_y[k] * smd_101[k];

        t_169[k] = f_5 * sld_65[k]
                   + f_1 * smp0_50[k]
                   - f_2 * smp1_50[k]
                   + f_3 * pc_z[k] * smd_101[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, pc_x, pc_y, pc_z, sld_66, sld_72, \
                         sld_102, sld_105, smp0_51, smp1_51, smd_102, \
                         smd_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_12 * sld_102[k]
                   + f_1 * smp0_51[k]
                   - f_2 * smp1_51[k]
                   + f_3 * pc_x[k] * smd_102[k];

        t_171[k] = f_10 * sld_72[k]
                   + f_3 * pc_y[k] * smd_102[k];

        t_172[k] = f_8 * sld_66[k]
                   + f_3 * pc_z[k] * smd_102[k];

        t_173[k] = f_12 * sld_105[k]
                   + f_3 * pc_x[k] * smd_105[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, pc_x, pc_y, pc_z, sld_69, sld_75, \
                         sld_106, sld_107, smp0_52, smp1_52, smd_105, smd_106, \
                         smd_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_12 * sld_106[k]
                   + f_3 * pc_x[k] * smd_106[k];

        t_175[k] = f_12 * sld_107[k]
                   + f_3 * pc_x[k] * smd_107[k];

        t_176[k] = f_10 * sld_75[k]
                   + f_1 * smp0_52[k]
                   - f_2 * smp1_52[k]
                   + f_3 * pc_y[k] * smd_105[k];

        t_177[k] = f_8 * sld_69[k]
                   + f_3 * pc_z[k] * smd_105[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, pc_x, pc_y, pc_z, sld_71, sld_77, sld_108, \
                         smp0_53, smp0_54, smp1_53, smp1_54, smd_107, \
                         smd_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_10 * sld_77[k]
                   + f_3 * pc_y[k] * smd_107[k];

        t_179[k] = f_8 * sld_71[k]
                   + f_1 * smp0_53[k]
                   - f_2 * smp1_53[k]
                   + f_3 * pc_z[k] * smd_107[k];

        t_180[k] = f_12 * sld_108[k]
                   + f_1 * smp0_54[k]
                   - f_2 * smp1_54[k]
                   + f_3 * pc_x[k] * smd_108[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pc_x, pc_y, pc_z, sld_72, sld_78, \
                         sld_111, sld_112, smd_108, smd_111, smd_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_8 * sld_78[k]
                   + f_3 * pc_y[k] * smd_108[k];

        t_182[k] = f_10 * sld_72[k]
                   + f_3 * pc_z[k] * smd_108[k];

        t_183[k] = f_12 * sld_111[k]
                   + f_3 * pc_x[k] * smd_111[k];

        t_184[k] = f_12 * sld_112[k]
                   + f_3 * pc_x[k] * smd_112[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pc_x, pc_y, pc_z, sld_75, sld_81, sld_83, \
                         sld_113, smp0_55, smp1_55, smd_111, smd_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_12 * sld_113[k]
                   + f_3 * pc_x[k] * smd_113[k];

        t_186[k] = f_8 * sld_81[k]
                   + f_1 * smp0_55[k]
                   - f_2 * smp1_55[k]
                   + f_3 * pc_y[k] * smd_111[k];

        t_187[k] = f_10 * sld_75[k]
                   + f_3 * pc_z[k] * smd_111[k];

        t_188[k] = f_8 * sld_83[k]
                   + f_3 * pc_y[k] * smd_113[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, pb_y, pc_y, pc_z, slf0_140, sld_77, \
                         sld_78, sld_84, slf1_140, smp0_56, smp1_56, smd_113, \
                         smd_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_10 * sld_77[k]
                   + f_1 * smp0_56[k]
                   - f_2 * smp1_56[k]
                   + f_3 * pc_z[k] * smd_113[k];

        t_190[k] = pb_y[k] * slf0_140[k]
                   - f_4 * pc_y[k] * slf1_140[k];

        t_191[k] = f_5 * sld_84[k]
                   + f_3 * pc_y[k] * smd_114[k];

        t_192[k] = f_12 * sld_78[k]
                   + f_3 * pc_z[k] * smd_114[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, pc_x, pc_y, sld_87, sld_117, sld_118, \
                         sld_119, smp0_58, smp1_58, smd_117, smd_118, \
                         smd_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_12 * sld_117[k]
                   + f_3 * pc_x[k] * smd_117[k];

        t_194[k] = f_12 * sld_118[k]
                   + f_3 * pc_x[k] * smd_118[k];

        t_195[k] = f_12 * sld_119[k]
                   + f_3 * pc_x[k] * smd_119[k];

        t_196[k] = f_5 * sld_87[k]
                   + f_1 * smp0_58[k]
                   - f_2 * smp1_58[k]
                   + f_3 * pc_y[k] * smd_117[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, pb_y, pc_y, pc_z, slf0_149, sld_81, sld_89, \
                         slf1_149, smd_117, smd_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_12 * sld_81[k]
                   + f_3 * pc_z[k] * smd_117[k];

        t_198[k] = f_5 * sld_89[k]
                   + f_3 * pc_y[k] * smd_119[k];

        t_199[k] = pb_y[k] * slf0_149[k]
                   - f_4 * pc_y[k] * slf1_149[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, pc_x, pc_y, pc_z, sld_84, sld_120, \
                         sld_123, smp0_60, smp1_60, smd_120, smd_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_12 * sld_120[k]
                   + f_1 * smp0_60[k]
                   - f_2 * smp1_60[k]
                   + f_3 * pc_x[k] * smd_120[k];

        t_201[k] = f_3 * pc_y[k] * smd_120[k];

        t_202[k] = f_11 * sld_84[k]
                   + f_3 * pc_z[k] * smd_120[k];

        t_203[k] = f_12 * sld_123[k]
                   + f_3 * pc_x[k] * smd_123[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, t_208, pc_x, pc_y, pc_z, sld_87, sld_124, \
                         sld_125, smp0_61, smp1_61, smd_123, smd_124, \
                         smd_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_12 * sld_124[k]
                   + f_3 * pc_x[k] * smd_124[k];

        t_205[k] = f_12 * sld_125[k]
                   + f_3 * pc_x[k] * smd_125[k];

        t_206[k] = f_1 * smp0_61[k]
                   - f_2 * smp1_61[k]
                   + f_3 * pc_y[k] * smd_123[k];

        t_207[k] = f_11 * sld_87[k]
                   + f_3 * pc_z[k] * smd_123[k];

        t_208[k] = f_3 * pc_y[k] * smd_125[k];
    }

#pragma omp simd aligned(t_209, t_210, t_211, t_212, pc_x, pc_y, pc_z, sld_89, sld_90, \
                         sld_126, smp0_62, smp0_63, smp1_62, smp1_63, smd_125, \
                         smd_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_209[k] = f_11 * sld_89[k]
                   + f_1 * smp0_62[k]
                   - f_2 * smp1_62[k]
                   + f_3 * pc_z[k] * smd_125[k];

        t_210[k] = f_10 * sld_126[k]
                   + f_1 * smp0_63[k]
                   - f_2 * smp1_63[k]
                   + f_3 * pc_x[k] * smd_126[k];

        t_211[k] = f_9 * sld_90[k]
                   + f_3 * pc_y[k] * smd_126[k];

        t_212[k] = f_3 * pc_z[k] * smd_126[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, pc_x, pc_y, sld_93, sld_129, sld_130, \
                         sld_131, smp0_64, smp1_64, smd_129, smd_130, \
                         smd_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_10 * sld_129[k]
                   + f_3 * pc_x[k] * smd_129[k];

        t_214[k] = f_10 * sld_130[k]
                   + f_3 * pc_x[k] * smd_130[k];

        t_215[k] = f_10 * sld_131[k]
                   + f_3 * pc_x[k] * smd_131[k];

        t_216[k] = f_9 * sld_93[k]
                   + f_1 * smp0_64[k]
                   - f_2 * smp1_64[k]
                   + f_3 * pc_y[k] * smd_129[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, pb_z, pc_y, pc_z, slf0_150, sld_95, \
                         slf1_150, smp0_65, smp1_65, smd_129, smd_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = f_3 * pc_z[k] * smd_129[k];

        t_218[k] = f_9 * sld_95[k]
                   + f_3 * pc_y[k] * smd_131[k];

        t_219[k] = f_1 * smp0_65[k]
                   - f_2 * smp1_65[k]
                   + f_3 * pc_z[k] * smd_131[k];

        t_220[k] = pb_z[k] * slf0_150[k]
                   - f_4 * pc_z[k] * slf1_150[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, pc_x, pc_y, pc_z, sld_90, sld_96, \
                         sld_135, sld_136, smd_132, smd_135, smd_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_11 * sld_96[k]
                   + f_3 * pc_y[k] * smd_132[k];

        t_222[k] = f_5 * sld_90[k]
                   + f_3 * pc_z[k] * smd_132[k];

        t_223[k] = f_10 * sld_135[k]
                   + f_3 * pc_x[k] * smd_135[k];

        t_224[k] = f_10 * sld_136[k]
                   + f_3 * pc_x[k] * smd_136[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, pb_z, pc_x, pc_y, pc_z, slf0_156, sld_93, \
                         sld_101, sld_137, slf1_156, smd_135, smd_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_10 * sld_137[k]
                   + f_3 * pc_x[k] * smd_137[k];

        t_226[k] = pb_z[k] * slf0_156[k]
                   - f_4 * pc_z[k] * slf1_156[k];

        t_227[k] = f_5 * sld_93[k]
                   + f_3 * pc_z[k] * smd_135[k];

        t_228[k] = f_11 * sld_101[k]
                   + f_3 * pc_y[k] * smd_137[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, pc_x, pc_y, pc_z, sld_95, sld_102, sld_138, \
                         smp0_68, smp0_69, smp1_68, smp1_69, smd_137, \
                         smd_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = f_5 * sld_95[k]
                   + f_1 * smp0_68[k]
                   - f_2 * smp1_68[k]
                   + f_3 * pc_z[k] * smd_137[k];

        t_230[k] = f_10 * sld_138[k]
                   + f_1 * smp0_69[k]
                   - f_2 * smp1_69[k]
                   + f_3 * pc_x[k] * smd_138[k];

        t_231[k] = f_12 * sld_102[k]
                   + f_3 * pc_y[k] * smd_138[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, t_235, pc_x, pc_z, sld_96, sld_141, sld_142, \
                         sld_143, smd_138, smd_141, smd_142, smd_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = f_8 * sld_96[k]
                   + f_3 * pc_z[k] * smd_138[k];

        t_233[k] = f_10 * sld_141[k]
                   + f_3 * pc_x[k] * smd_141[k];

        t_234[k] = f_10 * sld_142[k]
                   + f_3 * pc_x[k] * smd_142[k];

        t_235[k] = f_10 * sld_143[k]
                   + f_3 * pc_x[k] * smd_143[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, t_239, pc_y, pc_z, sld_99, sld_101, sld_105, \
                         sld_107, smp0_70, smp0_71, smp1_70, smp1_71, smd_141, \
                         smd_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = f_12 * sld_105[k]
                   + f_1 * smp0_70[k]
                   - f_2 * smp1_70[k]
                   + f_3 * pc_y[k] * smd_141[k];

        t_237[k] = f_8 * sld_99[k]
                   + f_3 * pc_z[k] * smd_141[k];

        t_238[k] = f_12 * sld_107[k]
                   + f_3 * pc_y[k] * smd_143[k];

        t_239[k] = f_8 * sld_101[k]
                   + f_1 * smp0_71[k]
                   - f_2 * smp1_71[k]
                   + f_3 * pc_z[k] * smd_143[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, pc_x, pc_y, pc_z, sld_102, sld_108, \
                         sld_144, sld_147, smp0_72, smp1_72, smd_144, \
                         smd_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_10 * sld_144[k]
                   + f_1 * smp0_72[k]
                   - f_2 * smp1_72[k]
                   + f_3 * pc_x[k] * smd_144[k];

        t_241[k] = f_10 * sld_108[k]
                   + f_3 * pc_y[k] * smd_144[k];

        t_242[k] = f_10 * sld_102[k]
                   + f_3 * pc_z[k] * smd_144[k];

        t_243[k] = f_10 * sld_147[k]
                   + f_3 * pc_x[k] * smd_147[k];
    }

#pragma omp simd aligned(t_244, t_245, t_246, t_247, pc_x, pc_y, pc_z, sld_105, sld_111, \
                         sld_148, sld_149, smp0_73, smp1_73, smd_147, smd_148, \
                         smd_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_244[k] = f_10 * sld_148[k]
                   + f_3 * pc_x[k] * smd_148[k];

        t_245[k] = f_10 * sld_149[k]
                   + f_3 * pc_x[k] * smd_149[k];

        t_246[k] = f_10 * sld_111[k]
                   + f_1 * smp0_73[k]
                   - f_2 * smp1_73[k]
                   + f_3 * pc_y[k] * smd_147[k];

        t_247[k] = f_10 * sld_105[k]
                   + f_3 * pc_z[k] * smd_147[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, pc_x, pc_y, pc_z, sld_107, sld_113, sld_150, \
                         smp0_74, smp0_75, smp1_74, smp1_75, smd_149, \
                         smd_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_10 * sld_113[k]
                   + f_3 * pc_y[k] * smd_149[k];

        t_249[k] = f_10 * sld_107[k]
                   + f_1 * smp0_74[k]
                   - f_2 * smp1_74[k]
                   + f_3 * pc_z[k] * smd_149[k];

        t_250[k] = f_10 * sld_150[k]
                   + f_1 * smp0_75[k]
                   - f_2 * smp1_75[k]
                   + f_3 * pc_x[k] * smd_150[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, pc_x, pc_y, pc_z, sld_108, sld_114, \
                         sld_153, sld_154, smd_150, smd_153, smd_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = f_8 * sld_114[k]
                   + f_3 * pc_y[k] * smd_150[k];

        t_252[k] = f_12 * sld_108[k]
                   + f_3 * pc_z[k] * smd_150[k];

        t_253[k] = f_10 * sld_153[k]
                   + f_3 * pc_x[k] * smd_153[k];

        t_254[k] = f_10 * sld_154[k]
                   + f_3 * pc_x[k] * smd_154[k];
    }
}

static auto
compute_prim_smf_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t slf0,
                                                          const size_t sld, const size_t slf1,
                                                          const size_t smp0, const size_t smp1,
                                                          const size_t smd, const size_t ncols,
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
    const auto f_6 = 4.0 / q;
    const auto f_7 = 3.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 3.0 / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 2.5 / q;
    const auto f_12 = 2.0 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *slf0_200 = buffer.data(slf0 + 200);
    const auto *slf0_209 = buffer.data(slf0 + 209);
    const auto *slf0_210 = buffer.data(slf0 + 210);
    const auto *slf0_216 = buffer.data(slf0 + 216);
    const auto *slf0_270 = buffer.data(slf0 + 270);
    const auto *slf0_279 = buffer.data(slf0 + 279);
    const auto *slf0_280 = buffer.data(slf0 + 280);
    const auto *slf0_360 = buffer.data(slf0 + 360);
    const auto *slf0_366 = buffer.data(slf0 + 366);
    const auto *slf0_369 = buffer.data(slf0 + 369);
    const auto *slf0_376 = buffer.data(slf0 + 376);

    const auto *sld_111 = buffer.data(sld + 111);
    const auto *sld_113 = buffer.data(sld + 113);
    const auto *sld_114 = buffer.data(sld + 114);
    const auto *sld_117 = buffer.data(sld + 117);
    const auto *sld_119 = buffer.data(sld + 119);
    const auto *sld_120 = buffer.data(sld + 120);
    const auto *sld_123 = buffer.data(sld + 123);
    const auto *sld_125 = buffer.data(sld + 125);
    const auto *sld_126 = buffer.data(sld + 126);
    const auto *sld_129 = buffer.data(sld + 129);
    const auto *sld_131 = buffer.data(sld + 131);
    const auto *sld_132 = buffer.data(sld + 132);
    const auto *sld_135 = buffer.data(sld + 135);
    const auto *sld_137 = buffer.data(sld + 137);
    const auto *sld_138 = buffer.data(sld + 138);
    const auto *sld_141 = buffer.data(sld + 141);
    const auto *sld_143 = buffer.data(sld + 143);
    const auto *sld_144 = buffer.data(sld + 144);
    const auto *sld_147 = buffer.data(sld + 147);
    const auto *sld_149 = buffer.data(sld + 149);
    const auto *sld_150 = buffer.data(sld + 150);
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
    const auto *sld_225 = buffer.data(sld + 225);
    const auto *sld_226 = buffer.data(sld + 226);
    const auto *sld_227 = buffer.data(sld + 227);

    const auto *slf1_200 = buffer.data(slf1 + 200);
    const auto *slf1_209 = buffer.data(slf1 + 209);
    const auto *slf1_210 = buffer.data(slf1 + 210);
    const auto *slf1_216 = buffer.data(slf1 + 216);
    const auto *slf1_270 = buffer.data(slf1 + 270);
    const auto *slf1_279 = buffer.data(slf1 + 279);
    const auto *slf1_280 = buffer.data(slf1 + 280);
    const auto *slf1_360 = buffer.data(slf1 + 360);
    const auto *slf1_366 = buffer.data(slf1 + 366);
    const auto *slf1_369 = buffer.data(slf1 + 369);
    const auto *slf1_376 = buffer.data(slf1 + 376);

    const auto *smp0_76 = buffer.data(smp0 + 76);
    const auto *smp0_77 = buffer.data(smp0 + 77);
    const auto *smp0_79 = buffer.data(smp0 + 79);
    const auto *smp0_81 = buffer.data(smp0 + 81);
    const auto *smp0_82 = buffer.data(smp0 + 82);
    const auto *smp0_83 = buffer.data(smp0 + 83);
    const auto *smp0_84 = buffer.data(smp0 + 84);
    const auto *smp0_85 = buffer.data(smp0 + 85);
    const auto *smp0_86 = buffer.data(smp0 + 86);
    const auto *smp0_89 = buffer.data(smp0 + 89);
    const auto *smp0_90 = buffer.data(smp0 + 90);
    const auto *smp0_91 = buffer.data(smp0 + 91);
    const auto *smp0_92 = buffer.data(smp0 + 92);
    const auto *smp0_93 = buffer.data(smp0 + 93);
    const auto *smp0_94 = buffer.data(smp0 + 94);
    const auto *smp0_95 = buffer.data(smp0 + 95);
    const auto *smp0_96 = buffer.data(smp0 + 96);
    const auto *smp0_97 = buffer.data(smp0 + 97);
    const auto *smp0_98 = buffer.data(smp0 + 98);
    const auto *smp0_99 = buffer.data(smp0 + 99);
    const auto *smp0_100 = buffer.data(smp0 + 100);
    const auto *smp0_101 = buffer.data(smp0 + 101);
    const auto *smp0_103 = buffer.data(smp0 + 103);
    const auto *smp0_105 = buffer.data(smp0 + 105);
    const auto *smp0_106 = buffer.data(smp0 + 106);
    const auto *smp0_107 = buffer.data(smp0 + 107);

    const auto *smp1_76 = buffer.data(smp1 + 76);
    const auto *smp1_77 = buffer.data(smp1 + 77);
    const auto *smp1_79 = buffer.data(smp1 + 79);
    const auto *smp1_81 = buffer.data(smp1 + 81);
    const auto *smp1_82 = buffer.data(smp1 + 82);
    const auto *smp1_83 = buffer.data(smp1 + 83);
    const auto *smp1_84 = buffer.data(smp1 + 84);
    const auto *smp1_85 = buffer.data(smp1 + 85);
    const auto *smp1_86 = buffer.data(smp1 + 86);
    const auto *smp1_89 = buffer.data(smp1 + 89);
    const auto *smp1_90 = buffer.data(smp1 + 90);
    const auto *smp1_91 = buffer.data(smp1 + 91);
    const auto *smp1_92 = buffer.data(smp1 + 92);
    const auto *smp1_93 = buffer.data(smp1 + 93);
    const auto *smp1_94 = buffer.data(smp1 + 94);
    const auto *smp1_95 = buffer.data(smp1 + 95);
    const auto *smp1_96 = buffer.data(smp1 + 96);
    const auto *smp1_97 = buffer.data(smp1 + 97);
    const auto *smp1_98 = buffer.data(smp1 + 98);
    const auto *smp1_99 = buffer.data(smp1 + 99);
    const auto *smp1_100 = buffer.data(smp1 + 100);
    const auto *smp1_101 = buffer.data(smp1 + 101);
    const auto *smp1_103 = buffer.data(smp1 + 103);
    const auto *smp1_105 = buffer.data(smp1 + 105);
    const auto *smp1_106 = buffer.data(smp1 + 106);
    const auto *smp1_107 = buffer.data(smp1 + 107);

    const auto *smd_153 = buffer.data(smd + 153);
    const auto *smd_155 = buffer.data(smd + 155);
    const auto *smd_156 = buffer.data(smd + 156);
    const auto *smd_159 = buffer.data(smd + 159);
    const auto *smd_160 = buffer.data(smd + 160);
    const auto *smd_161 = buffer.data(smd + 161);
    const auto *smd_162 = buffer.data(smd + 162);
    const auto *smd_165 = buffer.data(smd + 165);
    const auto *smd_166 = buffer.data(smd + 166);
    const auto *smd_167 = buffer.data(smd + 167);
    const auto *smd_168 = buffer.data(smd + 168);
    const auto *smd_171 = buffer.data(smd + 171);
    const auto *smd_172 = buffer.data(smd + 172);
    const auto *smd_173 = buffer.data(smd + 173);
    const auto *smd_174 = buffer.data(smd + 174);
    const auto *smd_177 = buffer.data(smd + 177);
    const auto *smd_178 = buffer.data(smd + 178);
    const auto *smd_179 = buffer.data(smd + 179);
    const auto *smd_180 = buffer.data(smd + 180);
    const auto *smd_183 = buffer.data(smd + 183);
    const auto *smd_184 = buffer.data(smd + 184);
    const auto *smd_185 = buffer.data(smd + 185);
    const auto *smd_186 = buffer.data(smd + 186);
    const auto *smd_189 = buffer.data(smd + 189);
    const auto *smd_190 = buffer.data(smd + 190);
    const auto *smd_191 = buffer.data(smd + 191);
    const auto *smd_192 = buffer.data(smd + 192);
    const auto *smd_195 = buffer.data(smd + 195);
    const auto *smd_196 = buffer.data(smd + 196);
    const auto *smd_197 = buffer.data(smd + 197);
    const auto *smd_198 = buffer.data(smd + 198);
    const auto *smd_201 = buffer.data(smd + 201);
    const auto *smd_202 = buffer.data(smd + 202);
    const auto *smd_203 = buffer.data(smd + 203);
    const auto *smd_204 = buffer.data(smd + 204);
    const auto *smd_207 = buffer.data(smd + 207);
    const auto *smd_208 = buffer.data(smd + 208);
    const auto *smd_209 = buffer.data(smd + 209);
    const auto *smd_210 = buffer.data(smd + 210);
    const auto *smd_213 = buffer.data(smd + 213);
    const auto *smd_214 = buffer.data(smd + 214);
    const auto *smd_215 = buffer.data(smd + 215);
    const auto *smd_216 = buffer.data(smd + 216);
    const auto *smd_219 = buffer.data(smd + 219);
    const auto *smd_220 = buffer.data(smd + 220);
    const auto *smd_221 = buffer.data(smd + 221);
    const auto *smd_222 = buffer.data(smd + 222);
    const auto *smd_225 = buffer.data(smd + 225);
    const auto *smd_226 = buffer.data(smd + 226);
    const auto *smd_227 = buffer.data(smd + 227);

#pragma omp simd aligned(t_255, t_256, t_257, t_258, pc_x, pc_y, pc_z, sld_111, sld_117, \
                         sld_119, sld_155, smp0_76, smp1_76, smd_153, \
                         smd_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_10 * sld_155[k]
                   + f_3 * pc_x[k] * smd_155[k];

        t_256[k] = f_8 * sld_117[k]
                   + f_1 * smp0_76[k]
                   - f_2 * smp1_76[k]
                   + f_3 * pc_y[k] * smd_153[k];

        t_257[k] = f_12 * sld_111[k]
                   + f_3 * pc_z[k] * smd_153[k];

        t_258[k] = f_8 * sld_119[k]
                   + f_3 * pc_y[k] * smd_155[k];
    }

#pragma omp simd aligned(t_259, t_260, t_261, t_262, pb_y, pc_y, pc_z, slf0_200, sld_113, \
                         sld_114, sld_120, slf1_200, smp0_77, smp1_77, smd_155, \
                         smd_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_259[k] = f_12 * sld_113[k]
                   + f_1 * smp0_77[k]
                   - f_2 * smp1_77[k]
                   + f_3 * pc_z[k] * smd_155[k];

        t_260[k] = pb_y[k] * slf0_200[k]
                   - f_4 * pc_y[k] * slf1_200[k];

        t_261[k] = f_5 * sld_120[k]
                   + f_3 * pc_y[k] * smd_156[k];

        t_262[k] = f_11 * sld_114[k]
                   + f_3 * pc_z[k] * smd_156[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, t_266, pc_x, pc_y, sld_123, sld_159, sld_160, \
                         sld_161, smp0_79, smp1_79, smd_159, smd_160, \
                         smd_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = f_10 * sld_159[k]
                   + f_3 * pc_x[k] * smd_159[k];

        t_264[k] = f_10 * sld_160[k]
                   + f_3 * pc_x[k] * smd_160[k];

        t_265[k] = f_10 * sld_161[k]
                   + f_3 * pc_x[k] * smd_161[k];

        t_266[k] = f_5 * sld_123[k]
                   + f_1 * smp0_79[k]
                   - f_2 * smp1_79[k]
                   + f_3 * pc_y[k] * smd_159[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, pb_y, pc_y, pc_z, slf0_209, sld_117, sld_125, \
                         slf1_209, smd_159, smd_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_11 * sld_117[k]
                   + f_3 * pc_z[k] * smd_159[k];

        t_268[k] = f_5 * sld_125[k]
                   + f_3 * pc_y[k] * smd_161[k];

        t_269[k] = pb_y[k] * slf0_209[k]
                   - f_4 * pc_y[k] * slf1_209[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, pc_x, pc_y, pc_z, sld_120, sld_162, \
                         sld_165, smp0_81, smp1_81, smd_162, smd_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_10 * sld_162[k]
                   + f_1 * smp0_81[k]
                   - f_2 * smp1_81[k]
                   + f_3 * pc_x[k] * smd_162[k];

        t_271[k] = f_3 * pc_y[k] * smd_162[k];

        t_272[k] = f_9 * sld_120[k]
                   + f_3 * pc_z[k] * smd_162[k];

        t_273[k] = f_10 * sld_165[k]
                   + f_3 * pc_x[k] * smd_165[k];
    }

#pragma omp simd aligned(t_274, t_275, t_276, t_277, t_278, pc_x, pc_y, pc_z, sld_123, \
                         sld_166, sld_167, smp0_82, smp1_82, smd_165, smd_166, \
                         smd_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_274[k] = f_10 * sld_166[k]
                   + f_3 * pc_x[k] * smd_166[k];

        t_275[k] = f_10 * sld_167[k]
                   + f_3 * pc_x[k] * smd_167[k];

        t_276[k] = f_1 * smp0_82[k]
                   - f_2 * smp1_82[k]
                   + f_3 * pc_y[k] * smd_165[k];

        t_277[k] = f_9 * sld_123[k]
                   + f_3 * pc_z[k] * smd_165[k];

        t_278[k] = f_3 * pc_y[k] * smd_167[k];
    }

#pragma omp simd aligned(t_279, t_280, t_281, t_282, pc_x, pc_y, pc_z, sld_125, sld_126, \
                         sld_168, smp0_83, smp0_84, smp1_83, smp1_84, smd_167, \
                         smd_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_279[k] = f_9 * sld_125[k]
                   + f_1 * smp0_83[k]
                   - f_2 * smp1_83[k]
                   + f_3 * pc_z[k] * smd_167[k];

        t_280[k] = f_8 * sld_168[k]
                   + f_1 * smp0_84[k]
                   - f_2 * smp1_84[k]
                   + f_3 * pc_x[k] * smd_168[k];

        t_281[k] = f_7 * sld_126[k]
                   + f_3 * pc_y[k] * smd_168[k];

        t_282[k] = f_3 * pc_z[k] * smd_168[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, t_286, pc_x, pc_y, sld_129, sld_171, sld_172, \
                         sld_173, smp0_85, smp1_85, smd_171, smd_172, \
                         smd_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_8 * sld_171[k]
                   + f_3 * pc_x[k] * smd_171[k];

        t_284[k] = f_8 * sld_172[k]
                   + f_3 * pc_x[k] * smd_172[k];

        t_285[k] = f_8 * sld_173[k]
                   + f_3 * pc_x[k] * smd_173[k];

        t_286[k] = f_7 * sld_129[k]
                   + f_1 * smp0_85[k]
                   - f_2 * smp1_85[k]
                   + f_3 * pc_y[k] * smd_171[k];
    }

#pragma omp simd aligned(t_287, t_288, t_289, t_290, pb_z, pc_y, pc_z, slf0_210, sld_131, \
                         slf1_210, smp0_86, smp1_86, smd_171, smd_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_287[k] = f_3 * pc_z[k] * smd_171[k];

        t_288[k] = f_7 * sld_131[k]
                   + f_3 * pc_y[k] * smd_173[k];

        t_289[k] = f_1 * smp0_86[k]
                   - f_2 * smp1_86[k]
                   + f_3 * pc_z[k] * smd_173[k];

        t_290[k] = pb_z[k] * slf0_210[k]
                   - f_4 * pc_z[k] * slf1_210[k];
    }

#pragma omp simd aligned(t_291, t_292, t_293, t_294, pc_x, pc_y, pc_z, sld_126, sld_132, \
                         sld_177, sld_178, smd_174, smd_177, smd_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_291[k] = f_9 * sld_132[k]
                   + f_3 * pc_y[k] * smd_174[k];

        t_292[k] = f_5 * sld_126[k]
                   + f_3 * pc_z[k] * smd_174[k];

        t_293[k] = f_8 * sld_177[k]
                   + f_3 * pc_x[k] * smd_177[k];

        t_294[k] = f_8 * sld_178[k]
                   + f_3 * pc_x[k] * smd_178[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, pb_z, pc_x, pc_y, pc_z, slf0_216, \
                         sld_129, sld_137, sld_179, slf1_216, smd_177, \
                         smd_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = f_8 * sld_179[k]
                   + f_3 * pc_x[k] * smd_179[k];

        t_296[k] = pb_z[k] * slf0_216[k]
                   - f_4 * pc_z[k] * slf1_216[k];

        t_297[k] = f_5 * sld_129[k]
                   + f_3 * pc_z[k] * smd_177[k];

        t_298[k] = f_9 * sld_137[k]
                   + f_3 * pc_y[k] * smd_179[k];
    }

#pragma omp simd aligned(t_299, t_300, t_301, pc_x, pc_y, pc_z, sld_131, sld_138, sld_180, \
                         smp0_89, smp0_90, smp1_89, smp1_90, smd_179, \
                         smd_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_299[k] = f_5 * sld_131[k]
                   + f_1 * smp0_89[k]
                   - f_2 * smp1_89[k]
                   + f_3 * pc_z[k] * smd_179[k];

        t_300[k] = f_8 * sld_180[k]
                   + f_1 * smp0_90[k]
                   - f_2 * smp1_90[k]
                   + f_3 * pc_x[k] * smd_180[k];

        t_301[k] = f_11 * sld_138[k]
                   + f_3 * pc_y[k] * smd_180[k];
    }

#pragma omp simd aligned(t_302, t_303, t_304, t_305, pc_x, pc_z, sld_132, sld_183, sld_184, \
                         sld_185, smd_180, smd_183, smd_184, smd_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_302[k] = f_8 * sld_132[k]
                   + f_3 * pc_z[k] * smd_180[k];

        t_303[k] = f_8 * sld_183[k]
                   + f_3 * pc_x[k] * smd_183[k];

        t_304[k] = f_8 * sld_184[k]
                   + f_3 * pc_x[k] * smd_184[k];

        t_305[k] = f_8 * sld_185[k]
                   + f_3 * pc_x[k] * smd_185[k];
    }

#pragma omp simd aligned(t_306, t_307, t_308, t_309, pc_y, pc_z, sld_135, sld_137, sld_141, \
                         sld_143, smp0_91, smp0_92, smp1_91, smp1_92, smd_183, \
                         smd_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_306[k] = f_11 * sld_141[k]
                   + f_1 * smp0_91[k]
                   - f_2 * smp1_91[k]
                   + f_3 * pc_y[k] * smd_183[k];

        t_307[k] = f_8 * sld_135[k]
                   + f_3 * pc_z[k] * smd_183[k];

        t_308[k] = f_11 * sld_143[k]
                   + f_3 * pc_y[k] * smd_185[k];

        t_309[k] = f_8 * sld_137[k]
                   + f_1 * smp0_92[k]
                   - f_2 * smp1_92[k]
                   + f_3 * pc_z[k] * smd_185[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, pc_x, pc_y, pc_z, sld_138, sld_144, \
                         sld_186, sld_189, smp0_93, smp1_93, smd_186, \
                         smd_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = f_8 * sld_186[k]
                   + f_1 * smp0_93[k]
                   - f_2 * smp1_93[k]
                   + f_3 * pc_x[k] * smd_186[k];

        t_311[k] = f_12 * sld_144[k]
                   + f_3 * pc_y[k] * smd_186[k];

        t_312[k] = f_10 * sld_138[k]
                   + f_3 * pc_z[k] * smd_186[k];

        t_313[k] = f_8 * sld_189[k]
                   + f_3 * pc_x[k] * smd_189[k];
    }

#pragma omp simd aligned(t_314, t_315, t_316, t_317, pc_x, pc_y, pc_z, sld_141, sld_147, \
                         sld_190, sld_191, smp0_94, smp1_94, smd_189, smd_190, \
                         smd_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_314[k] = f_8 * sld_190[k]
                   + f_3 * pc_x[k] * smd_190[k];

        t_315[k] = f_8 * sld_191[k]
                   + f_3 * pc_x[k] * smd_191[k];

        t_316[k] = f_12 * sld_147[k]
                   + f_1 * smp0_94[k]
                   - f_2 * smp1_94[k]
                   + f_3 * pc_y[k] * smd_189[k];

        t_317[k] = f_10 * sld_141[k]
                   + f_3 * pc_z[k] * smd_189[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, pc_x, pc_y, pc_z, sld_143, sld_149, sld_192, \
                         smp0_95, smp0_96, smp1_95, smp1_96, smd_191, \
                         smd_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = f_12 * sld_149[k]
                   + f_3 * pc_y[k] * smd_191[k];

        t_319[k] = f_10 * sld_143[k]
                   + f_1 * smp0_95[k]
                   - f_2 * smp1_95[k]
                   + f_3 * pc_z[k] * smd_191[k];

        t_320[k] = f_8 * sld_192[k]
                   + f_1 * smp0_96[k]
                   - f_2 * smp1_96[k]
                   + f_3 * pc_x[k] * smd_192[k];
    }

#pragma omp simd aligned(t_321, t_322, t_323, t_324, pc_x, pc_y, pc_z, sld_144, sld_150, \
                         sld_195, sld_196, smd_192, smd_195, smd_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_321[k] = f_10 * sld_150[k]
                   + f_3 * pc_y[k] * smd_192[k];

        t_322[k] = f_12 * sld_144[k]
                   + f_3 * pc_z[k] * smd_192[k];

        t_323[k] = f_8 * sld_195[k]
                   + f_3 * pc_x[k] * smd_195[k];

        t_324[k] = f_8 * sld_196[k]
                   + f_3 * pc_x[k] * smd_196[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, pc_x, pc_y, pc_z, sld_147, sld_153, \
                         sld_155, sld_197, smp0_97, smp1_97, smd_195, \
                         smd_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = f_8 * sld_197[k]
                   + f_3 * pc_x[k] * smd_197[k];

        t_326[k] = f_10 * sld_153[k]
                   + f_1 * smp0_97[k]
                   - f_2 * smp1_97[k]
                   + f_3 * pc_y[k] * smd_195[k];

        t_327[k] = f_12 * sld_147[k]
                   + f_3 * pc_z[k] * smd_195[k];

        t_328[k] = f_10 * sld_155[k]
                   + f_3 * pc_y[k] * smd_197[k];
    }

#pragma omp simd aligned(t_329, t_330, t_331, pc_x, pc_y, pc_z, sld_149, sld_156, sld_198, \
                         smp0_98, smp0_99, smp1_98, smp1_99, smd_197, \
                         smd_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_329[k] = f_12 * sld_149[k]
                   + f_1 * smp0_98[k]
                   - f_2 * smp1_98[k]
                   + f_3 * pc_z[k] * smd_197[k];

        t_330[k] = f_8 * sld_198[k]
                   + f_1 * smp0_99[k]
                   - f_2 * smp1_99[k]
                   + f_3 * pc_x[k] * smd_198[k];

        t_331[k] = f_8 * sld_156[k]
                   + f_3 * pc_y[k] * smd_198[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, t_335, pc_x, pc_z, sld_150, sld_201, sld_202, \
                         sld_203, smd_198, smd_201, smd_202, smd_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = f_11 * sld_150[k]
                   + f_3 * pc_z[k] * smd_198[k];

        t_333[k] = f_8 * sld_201[k]
                   + f_3 * pc_x[k] * smd_201[k];

        t_334[k] = f_8 * sld_202[k]
                   + f_3 * pc_x[k] * smd_202[k];

        t_335[k] = f_8 * sld_203[k]
                   + f_3 * pc_x[k] * smd_203[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, t_339, pc_y, pc_z, sld_153, sld_155, sld_159, \
                         sld_161, smp0_100, smp0_101, smp1_100, smp1_101, smd_201, \
                         smd_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = f_8 * sld_159[k]
                   + f_1 * smp0_100[k]
                   - f_2 * smp1_100[k]
                   + f_3 * pc_y[k] * smd_201[k];

        t_337[k] = f_11 * sld_153[k]
                   + f_3 * pc_z[k] * smd_201[k];

        t_338[k] = f_8 * sld_161[k]
                   + f_3 * pc_y[k] * smd_203[k];

        t_339[k] = f_11 * sld_155[k]
                   + f_1 * smp0_101[k]
                   - f_2 * smp1_101[k]
                   + f_3 * pc_z[k] * smd_203[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, pb_y, pc_x, pc_y, pc_z, slf0_270, \
                         sld_156, sld_162, sld_207, slf1_270, smd_204, \
                         smd_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = pb_y[k] * slf0_270[k]
                   - f_4 * pc_y[k] * slf1_270[k];

        t_341[k] = f_5 * sld_162[k]
                   + f_3 * pc_y[k] * smd_204[k];

        t_342[k] = f_9 * sld_156[k]
                   + f_3 * pc_z[k] * smd_204[k];

        t_343[k] = f_8 * sld_207[k]
                   + f_3 * pc_x[k] * smd_207[k];
    }

#pragma omp simd aligned(t_344, t_345, t_346, t_347, pc_x, pc_y, pc_z, sld_159, sld_165, \
                         sld_208, sld_209, smp0_103, smp1_103, smd_207, smd_208, \
                         smd_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_344[k] = f_8 * sld_208[k]
                   + f_3 * pc_x[k] * smd_208[k];

        t_345[k] = f_8 * sld_209[k]
                   + f_3 * pc_x[k] * smd_209[k];

        t_346[k] = f_5 * sld_165[k]
                   + f_1 * smp0_103[k]
                   - f_2 * smp1_103[k]
                   + f_3 * pc_y[k] * smd_207[k];

        t_347[k] = f_9 * sld_159[k]
                   + f_3 * pc_z[k] * smd_207[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, t_351, pb_y, pc_x, pc_y, slf0_279, sld_167, \
                         sld_210, slf1_279, smp0_105, smp1_105, smd_209, \
                         smd_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = f_5 * sld_167[k]
                   + f_3 * pc_y[k] * smd_209[k];

        t_349[k] = pb_y[k] * slf0_279[k]
                   - f_4 * pc_y[k] * slf1_279[k];

        t_350[k] = f_8 * sld_210[k]
                   + f_1 * smp0_105[k]
                   - f_2 * smp1_105[k]
                   + f_3 * pc_x[k] * smd_210[k];

        t_351[k] = f_3 * pc_y[k] * smd_210[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, t_355, pc_x, pc_z, sld_162, sld_213, sld_214, \
                         sld_215, smd_210, smd_213, smd_214, smd_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = f_7 * sld_162[k]
                   + f_3 * pc_z[k] * smd_210[k];

        t_353[k] = f_8 * sld_213[k]
                   + f_3 * pc_x[k] * smd_213[k];

        t_354[k] = f_8 * sld_214[k]
                   + f_3 * pc_x[k] * smd_214[k];

        t_355[k] = f_8 * sld_215[k]
                   + f_3 * pc_x[k] * smd_215[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, t_359, pc_y, pc_z, sld_165, sld_167, smp0_106, \
                         smp0_107, smp1_106, smp1_107, smd_213, \
                         smd_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_1 * smp0_106[k]
                   - f_2 * smp1_106[k]
                   + f_3 * pc_y[k] * smd_213[k];

        t_357[k] = f_7 * sld_165[k]
                   + f_3 * pc_z[k] * smd_213[k];

        t_358[k] = f_3 * pc_y[k] * smd_215[k];

        t_359[k] = f_7 * sld_167[k]
                   + f_1 * smp0_107[k]
                   - f_2 * smp1_107[k]
                   + f_3 * pc_z[k] * smd_215[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, pb_x, pc_x, pc_y, pc_z, slf0_360, \
                         sld_168, sld_216, sld_219, slf1_360, smd_216, \
                         smd_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = pb_x[k] * slf0_360[k]
                   + f_10 * sld_216[k]
                   - f_4 * pc_x[k] * slf1_360[k];

        t_361[k] = f_6 * sld_168[k]
                   + f_3 * pc_y[k] * smd_216[k];

        t_362[k] = f_3 * pc_z[k] * smd_216[k];

        t_363[k] = f_5 * sld_219[k]
                   + f_3 * pc_x[k] * smd_219[k];
    }

#pragma omp simd aligned(t_364, t_365, t_366, t_367, pb_x, pc_x, pc_z, slf0_366, sld_220, \
                         sld_221, slf1_366, smd_219, smd_220, smd_221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_364[k] = f_5 * sld_220[k]
                   + f_3 * pc_x[k] * smd_220[k];

        t_365[k] = f_5 * sld_221[k]
                   + f_3 * pc_x[k] * smd_221[k];

        t_366[k] = pb_x[k] * slf0_366[k]
                   - f_4 * pc_x[k] * slf1_366[k];

        t_367[k] = f_3 * pc_z[k] * smd_219[k];
    }

#pragma omp simd aligned(t_368, t_369, t_370, pb_x, pb_z, pc_x, pc_y, pc_z, slf0_280, \
                         slf0_369, sld_173, slf1_280, slf1_369, \
                         smd_221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_368[k] = f_6 * sld_173[k]
                   + f_3 * pc_y[k] * smd_221[k];

        t_369[k] = pb_x[k] * slf0_369[k]
                   - f_4 * pc_x[k] * slf1_369[k];

        t_370[k] = pb_z[k] * slf0_280[k]
                   - f_4 * pc_z[k] * slf1_280[k];
    }

#pragma omp simd aligned(t_371, t_372, t_373, t_374, pc_x, pc_y, pc_z, sld_168, sld_174, \
                         sld_225, sld_226, smd_222, smd_225, smd_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_371[k] = f_7 * sld_174[k]
                   + f_3 * pc_y[k] * smd_222[k];

        t_372[k] = f_5 * sld_168[k]
                   + f_3 * pc_z[k] * smd_222[k];

        t_373[k] = f_5 * sld_225[k]
                   + f_3 * pc_x[k] * smd_225[k];

        t_374[k] = f_5 * sld_226[k]
                   + f_3 * pc_x[k] * smd_226[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, pb_x, pc_x, pc_y, pc_z, slf0_376, \
                         sld_171, sld_179, sld_227, slf1_376, smd_225, \
                         smd_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = f_5 * sld_227[k]
                   + f_3 * pc_x[k] * smd_227[k];

        t_376[k] = pb_x[k] * slf0_376[k]
                   - f_4 * pc_x[k] * slf1_376[k];

        t_377[k] = f_5 * sld_171[k]
                   + f_3 * pc_z[k] * smd_225[k];

        t_378[k] = f_7 * sld_179[k]
                   + f_3 * pc_y[k] * smd_227[k];
    }
}

static auto
compute_prim_smf_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t slf0,
                                                          const size_t sld, const size_t slf1,
                                                          const size_t smp0, const size_t smp1,
                                                          const size_t smd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / q;
    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 4.0 / q;
    const auto f_7 = 3.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 3.0 / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 2.5 / q;
    const auto f_12 = 2.0 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *slf0_350 = buffer.data(slf0 + 350);
    const auto *slf0_360 = buffer.data(slf0 + 360);
    const auto *slf0_366 = buffer.data(slf0 + 366);
    const auto *slf0_379 = buffer.data(slf0 + 379);
    const auto *slf0_380 = buffer.data(slf0 + 380);
    const auto *slf0_386 = buffer.data(slf0 + 386);
    const auto *slf0_389 = buffer.data(slf0 + 389);
    const auto *slf0_390 = buffer.data(slf0 + 390);
    const auto *slf0_396 = buffer.data(slf0 + 396);
    const auto *slf0_399 = buffer.data(slf0 + 399);
    const auto *slf0_400 = buffer.data(slf0 + 400);
    const auto *slf0_406 = buffer.data(slf0 + 406);
    const auto *slf0_409 = buffer.data(slf0 + 409);
    const auto *slf0_410 = buffer.data(slf0 + 410);
    const auto *slf0_416 = buffer.data(slf0 + 416);
    const auto *slf0_419 = buffer.data(slf0 + 419);
    const auto *slf0_420 = buffer.data(slf0 + 420);
    const auto *slf0_426 = buffer.data(slf0 + 426);
    const auto *slf0_429 = buffer.data(slf0 + 429);
    const auto *slf0_436 = buffer.data(slf0 + 436);
    const auto *slf0_439 = buffer.data(slf0 + 439);
    const auto *slf0_440 = buffer.data(slf0 + 440);
    const auto *slf0_446 = buffer.data(slf0 + 446);
    const auto *slf0_449 = buffer.data(slf0 + 449);

    const auto *sld_174 = buffer.data(sld + 174);
    const auto *sld_177 = buffer.data(sld + 177);
    const auto *sld_180 = buffer.data(sld + 180);
    const auto *sld_183 = buffer.data(sld + 183);
    const auto *sld_185 = buffer.data(sld + 185);
    const auto *sld_186 = buffer.data(sld + 186);
    const auto *sld_189 = buffer.data(sld + 189);
    const auto *sld_191 = buffer.data(sld + 191);
    const auto *sld_192 = buffer.data(sld + 192);
    const auto *sld_195 = buffer.data(sld + 195);
    const auto *sld_197 = buffer.data(sld + 197);
    const auto *sld_198 = buffer.data(sld + 198);
    const auto *sld_201 = buffer.data(sld + 201);
    const auto *sld_203 = buffer.data(sld + 203);
    const auto *sld_204 = buffer.data(sld + 204);
    const auto *sld_207 = buffer.data(sld + 207);
    const auto *sld_209 = buffer.data(sld + 209);
    const auto *sld_210 = buffer.data(sld + 210);
    const auto *sld_213 = buffer.data(sld + 213);
    const auto *sld_215 = buffer.data(sld + 215);
    const auto *sld_216 = buffer.data(sld + 216);
    const auto *sld_219 = buffer.data(sld + 219);
    const auto *sld_221 = buffer.data(sld + 221);
    const auto *sld_222 = buffer.data(sld + 222);
    const auto *sld_225 = buffer.data(sld + 225);
    const auto *sld_227 = buffer.data(sld + 227);
    const auto *sld_228 = buffer.data(sld + 228);
    const auto *sld_231 = buffer.data(sld + 231);
    const auto *sld_232 = buffer.data(sld + 232);
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
    const auto *sld_261 = buffer.data(sld + 261);
    const auto *sld_262 = buffer.data(sld + 262);
    const auto *sld_263 = buffer.data(sld + 263);
    const auto *sld_264 = buffer.data(sld + 264);
    const auto *sld_267 = buffer.data(sld + 267);
    const auto *sld_268 = buffer.data(sld + 268);
    const auto *sld_269 = buffer.data(sld + 269);

    const auto *slf1_350 = buffer.data(slf1 + 350);
    const auto *slf1_360 = buffer.data(slf1 + 360);
    const auto *slf1_366 = buffer.data(slf1 + 366);
    const auto *slf1_379 = buffer.data(slf1 + 379);
    const auto *slf1_380 = buffer.data(slf1 + 380);
    const auto *slf1_386 = buffer.data(slf1 + 386);
    const auto *slf1_389 = buffer.data(slf1 + 389);
    const auto *slf1_390 = buffer.data(slf1 + 390);
    const auto *slf1_396 = buffer.data(slf1 + 396);
    const auto *slf1_399 = buffer.data(slf1 + 399);
    const auto *slf1_400 = buffer.data(slf1 + 400);
    const auto *slf1_406 = buffer.data(slf1 + 406);
    const auto *slf1_409 = buffer.data(slf1 + 409);
    const auto *slf1_410 = buffer.data(slf1 + 410);
    const auto *slf1_416 = buffer.data(slf1 + 416);
    const auto *slf1_419 = buffer.data(slf1 + 419);
    const auto *slf1_420 = buffer.data(slf1 + 420);
    const auto *slf1_426 = buffer.data(slf1 + 426);
    const auto *slf1_429 = buffer.data(slf1 + 429);
    const auto *slf1_436 = buffer.data(slf1 + 436);
    const auto *slf1_439 = buffer.data(slf1 + 439);
    const auto *slf1_440 = buffer.data(slf1 + 440);
    const auto *slf1_446 = buffer.data(slf1 + 446);
    const auto *slf1_449 = buffer.data(slf1 + 449);

    const auto *smp0_135 = buffer.data(smp0 + 135);
    const auto *smp0_136 = buffer.data(smp0 + 136);
    const auto *smp0_137 = buffer.data(smp0 + 137);
    const auto *smp0_140 = buffer.data(smp0 + 140);
    const auto *smp0_141 = buffer.data(smp0 + 141);
    const auto *smp0_142 = buffer.data(smp0 + 142);
    const auto *smp0_143 = buffer.data(smp0 + 143);
    const auto *smp0_144 = buffer.data(smp0 + 144);
    const auto *smp0_145 = buffer.data(smp0 + 145);
    const auto *smp0_146 = buffer.data(smp0 + 146);
    const auto *smp0_147 = buffer.data(smp0 + 147);
    const auto *smp0_148 = buffer.data(smp0 + 148);
    const auto *smp0_149 = buffer.data(smp0 + 149);
    const auto *smp0_150 = buffer.data(smp0 + 150);
    const auto *smp0_151 = buffer.data(smp0 + 151);
    const auto *smp0_152 = buffer.data(smp0 + 152);
    const auto *smp0_153 = buffer.data(smp0 + 153);

    const auto *smp1_135 = buffer.data(smp1 + 135);
    const auto *smp1_136 = buffer.data(smp1 + 136);
    const auto *smp1_137 = buffer.data(smp1 + 137);
    const auto *smp1_140 = buffer.data(smp1 + 140);
    const auto *smp1_141 = buffer.data(smp1 + 141);
    const auto *smp1_142 = buffer.data(smp1 + 142);
    const auto *smp1_143 = buffer.data(smp1 + 143);
    const auto *smp1_144 = buffer.data(smp1 + 144);
    const auto *smp1_145 = buffer.data(smp1 + 145);
    const auto *smp1_146 = buffer.data(smp1 + 146);
    const auto *smp1_147 = buffer.data(smp1 + 147);
    const auto *smp1_148 = buffer.data(smp1 + 148);
    const auto *smp1_149 = buffer.data(smp1 + 149);
    const auto *smp1_150 = buffer.data(smp1 + 150);
    const auto *smp1_151 = buffer.data(smp1 + 151);
    const auto *smp1_152 = buffer.data(smp1 + 152);
    const auto *smp1_153 = buffer.data(smp1 + 153);

    const auto *smd_228 = buffer.data(smd + 228);
    const auto *smd_231 = buffer.data(smd + 231);
    const auto *smd_232 = buffer.data(smd + 232);
    const auto *smd_233 = buffer.data(smd + 233);
    const auto *smd_234 = buffer.data(smd + 234);
    const auto *smd_237 = buffer.data(smd + 237);
    const auto *smd_238 = buffer.data(smd + 238);
    const auto *smd_239 = buffer.data(smd + 239);
    const auto *smd_240 = buffer.data(smd + 240);
    const auto *smd_243 = buffer.data(smd + 243);
    const auto *smd_244 = buffer.data(smd + 244);
    const auto *smd_245 = buffer.data(smd + 245);
    const auto *smd_246 = buffer.data(smd + 246);
    const auto *smd_249 = buffer.data(smd + 249);
    const auto *smd_250 = buffer.data(smd + 250);
    const auto *smd_251 = buffer.data(smd + 251);
    const auto *smd_252 = buffer.data(smd + 252);
    const auto *smd_255 = buffer.data(smd + 255);
    const auto *smd_256 = buffer.data(smd + 256);
    const auto *smd_257 = buffer.data(smd + 257);
    const auto *smd_258 = buffer.data(smd + 258);
    const auto *smd_261 = buffer.data(smd + 261);
    const auto *smd_262 = buffer.data(smd + 262);
    const auto *smd_263 = buffer.data(smd + 263);
    const auto *smd_264 = buffer.data(smd + 264);
    const auto *smd_267 = buffer.data(smd + 267);
    const auto *smd_268 = buffer.data(smd + 268);
    const auto *smd_269 = buffer.data(smd + 269);
    const auto *smd_270 = buffer.data(smd + 270);
    const auto *smd_273 = buffer.data(smd + 273);
    const auto *smd_274 = buffer.data(smd + 274);
    const auto *smd_275 = buffer.data(smd + 275);
    const auto *smd_276 = buffer.data(smd + 276);
    const auto *smd_279 = buffer.data(smd + 279);
    const auto *smd_280 = buffer.data(smd + 280);
    const auto *smd_281 = buffer.data(smd + 281);
    const auto *smd_282 = buffer.data(smd + 282);
    const auto *smd_285 = buffer.data(smd + 285);
    const auto *smd_286 = buffer.data(smd + 286);
    const auto *smd_287 = buffer.data(smd + 287);
    const auto *smd_288 = buffer.data(smd + 288);
    const auto *smd_291 = buffer.data(smd + 291);
    const auto *smd_292 = buffer.data(smd + 292);
    const auto *smd_293 = buffer.data(smd + 293);
    const auto *smd_294 = buffer.data(smd + 294);
    const auto *smd_297 = buffer.data(smd + 297);
    const auto *smd_298 = buffer.data(smd + 298);
    const auto *smd_299 = buffer.data(smd + 299);
    const auto *smd_300 = buffer.data(smd + 300);
    const auto *smd_303 = buffer.data(smd + 303);
    const auto *smd_304 = buffer.data(smd + 304);
    const auto *smd_305 = buffer.data(smd + 305);
    const auto *smd_306 = buffer.data(smd + 306);

#pragma omp simd aligned(t_379, t_380, t_381, t_382, pb_x, pc_x, pc_y, pc_z, slf0_379, \
                         slf0_380, sld_174, sld_180, sld_228, slf1_379, slf1_380, \
                         smd_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_379[k] = pb_x[k] * slf0_379[k]
                   - f_4 * pc_x[k] * slf1_379[k];

        t_380[k] = pb_x[k] * slf0_380[k]
                   + f_10 * sld_228[k]
                   - f_4 * pc_x[k] * slf1_380[k];

        t_381[k] = f_9 * sld_180[k]
                   + f_3 * pc_y[k] * smd_228[k];

        t_382[k] = f_8 * sld_174[k]
                   + f_3 * pc_z[k] * smd_228[k];
    }

#pragma omp simd aligned(t_383, t_384, t_385, t_386, pb_x, pc_x, slf0_386, sld_231, sld_232, \
                         sld_233, slf1_386, smd_231, smd_232, smd_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_383[k] = f_5 * sld_231[k]
                   + f_3 * pc_x[k] * smd_231[k];

        t_384[k] = f_5 * sld_232[k]
                   + f_3 * pc_x[k] * smd_232[k];

        t_385[k] = f_5 * sld_233[k]
                   + f_3 * pc_x[k] * smd_233[k];

        t_386[k] = pb_x[k] * slf0_386[k]
                   - f_4 * pc_x[k] * slf1_386[k];
    }

#pragma omp simd aligned(t_387, t_388, t_389, pb_x, pc_x, pc_y, pc_z, slf0_389, sld_177, \
                         sld_185, slf1_389, smd_231, smd_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_387[k] = f_8 * sld_177[k]
                   + f_3 * pc_z[k] * smd_231[k];

        t_388[k] = f_9 * sld_185[k]
                   + f_3 * pc_y[k] * smd_233[k];

        t_389[k] = pb_x[k] * slf0_389[k]
                   - f_4 * pc_x[k] * slf1_389[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, pb_x, pc_x, pc_y, pc_z, slf0_390, \
                         sld_180, sld_186, sld_234, sld_237, slf1_390, smd_234, \
                         smd_237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = pb_x[k] * slf0_390[k]
                   + f_10 * sld_234[k]
                   - f_4 * pc_x[k] * slf1_390[k];

        t_391[k] = f_11 * sld_186[k]
                   + f_3 * pc_y[k] * smd_234[k];

        t_392[k] = f_10 * sld_180[k]
                   + f_3 * pc_z[k] * smd_234[k];

        t_393[k] = f_5 * sld_237[k]
                   + f_3 * pc_x[k] * smd_237[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, t_397, pb_x, pc_x, pc_z, slf0_396, sld_183, \
                         sld_238, sld_239, slf1_396, smd_237, smd_238, \
                         smd_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = f_5 * sld_238[k]
                   + f_3 * pc_x[k] * smd_238[k];

        t_395[k] = f_5 * sld_239[k]
                   + f_3 * pc_x[k] * smd_239[k];

        t_396[k] = pb_x[k] * slf0_396[k]
                   - f_4 * pc_x[k] * slf1_396[k];

        t_397[k] = f_10 * sld_183[k]
                   + f_3 * pc_z[k] * smd_237[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, t_401, pb_x, pc_x, pc_y, slf0_399, slf0_400, \
                         sld_191, sld_192, sld_240, slf1_399, slf1_400, smd_239, \
                         smd_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = f_11 * sld_191[k]
                   + f_3 * pc_y[k] * smd_239[k];

        t_399[k] = pb_x[k] * slf0_399[k]
                   - f_4 * pc_x[k] * slf1_399[k];

        t_400[k] = pb_x[k] * slf0_400[k]
                   + f_10 * sld_240[k]
                   - f_4 * pc_x[k] * slf1_400[k];

        t_401[k] = f_12 * sld_192[k]
                   + f_3 * pc_y[k] * smd_240[k];
    }

#pragma omp simd aligned(t_402, t_403, t_404, t_405, pc_x, pc_z, sld_186, sld_243, sld_244, \
                         sld_245, smd_240, smd_243, smd_244, smd_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_402[k] = f_12 * sld_186[k]
                   + f_3 * pc_z[k] * smd_240[k];

        t_403[k] = f_5 * sld_243[k]
                   + f_3 * pc_x[k] * smd_243[k];

        t_404[k] = f_5 * sld_244[k]
                   + f_3 * pc_x[k] * smd_244[k];

        t_405[k] = f_5 * sld_245[k]
                   + f_3 * pc_x[k] * smd_245[k];
    }

#pragma omp simd aligned(t_406, t_407, t_408, t_409, pb_x, pc_x, pc_y, pc_z, slf0_406, \
                         slf0_409, sld_189, sld_197, slf1_406, slf1_409, smd_243, \
                         smd_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_406[k] = pb_x[k] * slf0_406[k]
                   - f_4 * pc_x[k] * slf1_406[k];

        t_407[k] = f_12 * sld_189[k]
                   + f_3 * pc_z[k] * smd_243[k];

        t_408[k] = f_12 * sld_197[k]
                   + f_3 * pc_y[k] * smd_245[k];

        t_409[k] = pb_x[k] * slf0_409[k]
                   - f_4 * pc_x[k] * slf1_409[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, pb_x, pc_x, pc_y, pc_z, slf0_410, \
                         sld_192, sld_198, sld_246, sld_249, slf1_410, smd_246, \
                         smd_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = pb_x[k] * slf0_410[k]
                   + f_10 * sld_246[k]
                   - f_4 * pc_x[k] * slf1_410[k];

        t_411[k] = f_10 * sld_198[k]
                   + f_3 * pc_y[k] * smd_246[k];

        t_412[k] = f_11 * sld_192[k]
                   + f_3 * pc_z[k] * smd_246[k];

        t_413[k] = f_5 * sld_249[k]
                   + f_3 * pc_x[k] * smd_249[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, t_417, pb_x, pc_x, pc_z, slf0_416, sld_195, \
                         sld_250, sld_251, slf1_416, smd_249, smd_250, \
                         smd_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_5 * sld_250[k]
                   + f_3 * pc_x[k] * smd_250[k];

        t_415[k] = f_5 * sld_251[k]
                   + f_3 * pc_x[k] * smd_251[k];

        t_416[k] = pb_x[k] * slf0_416[k]
                   - f_4 * pc_x[k] * slf1_416[k];

        t_417[k] = f_11 * sld_195[k]
                   + f_3 * pc_z[k] * smd_249[k];
    }

#pragma omp simd aligned(t_418, t_419, t_420, t_421, pb_x, pc_x, pc_y, slf0_419, slf0_420, \
                         sld_203, sld_204, sld_252, slf1_419, slf1_420, smd_251, \
                         smd_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_418[k] = f_10 * sld_203[k]
                   + f_3 * pc_y[k] * smd_251[k];

        t_419[k] = pb_x[k] * slf0_419[k]
                   - f_4 * pc_x[k] * slf1_419[k];

        t_420[k] = pb_x[k] * slf0_420[k]
                   + f_10 * sld_252[k]
                   - f_4 * pc_x[k] * slf1_420[k];

        t_421[k] = f_8 * sld_204[k]
                   + f_3 * pc_y[k] * smd_252[k];
    }

#pragma omp simd aligned(t_422, t_423, t_424, t_425, pc_x, pc_z, sld_198, sld_255, sld_256, \
                         sld_257, smd_252, smd_255, smd_256, smd_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_422[k] = f_9 * sld_198[k]
                   + f_3 * pc_z[k] * smd_252[k];

        t_423[k] = f_5 * sld_255[k]
                   + f_3 * pc_x[k] * smd_255[k];

        t_424[k] = f_5 * sld_256[k]
                   + f_3 * pc_x[k] * smd_256[k];

        t_425[k] = f_5 * sld_257[k]
                   + f_3 * pc_x[k] * smd_257[k];
    }

#pragma omp simd aligned(t_426, t_427, t_428, t_429, pb_x, pc_x, pc_y, pc_z, slf0_426, \
                         slf0_429, sld_201, sld_209, slf1_426, slf1_429, smd_255, \
                         smd_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_426[k] = pb_x[k] * slf0_426[k]
                   - f_4 * pc_x[k] * slf1_426[k];

        t_427[k] = f_9 * sld_201[k]
                   + f_3 * pc_z[k] * smd_255[k];

        t_428[k] = f_8 * sld_209[k]
                   + f_3 * pc_y[k] * smd_257[k];

        t_429[k] = pb_x[k] * slf0_429[k]
                   - f_4 * pc_x[k] * slf1_429[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, pb_y, pc_x, pc_y, pc_z, slf0_350, \
                         sld_204, sld_210, sld_261, slf1_350, smd_258, \
                         smd_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = pb_y[k] * slf0_350[k]
                   - f_4 * pc_y[k] * slf1_350[k];

        t_431[k] = f_5 * sld_210[k]
                   + f_3 * pc_y[k] * smd_258[k];

        t_432[k] = f_7 * sld_204[k]
                   + f_3 * pc_z[k] * smd_258[k];

        t_433[k] = f_5 * sld_261[k]
                   + f_3 * pc_x[k] * smd_261[k];
    }

#pragma omp simd aligned(t_434, t_435, t_436, t_437, pb_x, pc_x, pc_z, slf0_436, sld_207, \
                         sld_262, sld_263, slf1_436, smd_261, smd_262, \
                         smd_263 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_434[k] = f_5 * sld_262[k]
                   + f_3 * pc_x[k] * smd_262[k];

        t_435[k] = f_5 * sld_263[k]
                   + f_3 * pc_x[k] * smd_263[k];

        t_436[k] = pb_x[k] * slf0_436[k]
                   - f_4 * pc_x[k] * slf1_436[k];

        t_437[k] = f_7 * sld_207[k]
                   + f_3 * pc_z[k] * smd_261[k];
    }

#pragma omp simd aligned(t_438, t_439, t_440, t_441, pb_x, pc_x, pc_y, slf0_439, slf0_440, \
                         sld_215, sld_264, slf1_439, slf1_440, smd_263, \
                         smd_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_438[k] = f_5 * sld_215[k]
                   + f_3 * pc_y[k] * smd_263[k];

        t_439[k] = pb_x[k] * slf0_439[k]
                   - f_4 * pc_x[k] * slf1_439[k];

        t_440[k] = pb_x[k] * slf0_440[k]
                   + f_10 * sld_264[k]
                   - f_4 * pc_x[k] * slf1_440[k];

        t_441[k] = f_3 * pc_y[k] * smd_264[k];
    }

#pragma omp simd aligned(t_442, t_443, t_444, t_445, pc_x, pc_z, sld_210, sld_267, sld_268, \
                         sld_269, smd_264, smd_267, smd_268, smd_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_442[k] = f_6 * sld_210[k]
                   + f_3 * pc_z[k] * smd_264[k];

        t_443[k] = f_5 * sld_267[k]
                   + f_3 * pc_x[k] * smd_267[k];

        t_444[k] = f_5 * sld_268[k]
                   + f_3 * pc_x[k] * smd_268[k];

        t_445[k] = f_5 * sld_269[k]
                   + f_3 * pc_x[k] * smd_269[k];
    }

#pragma omp simd aligned(t_446, t_447, t_448, t_449, pb_x, pc_x, pc_y, pc_z, slf0_446, \
                         slf0_449, sld_213, slf1_446, slf1_449, smd_267, \
                         smd_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_446[k] = pb_x[k] * slf0_446[k]
                   - f_4 * pc_x[k] * slf1_446[k];

        t_447[k] = f_6 * sld_213[k]
                   + f_3 * pc_z[k] * smd_267[k];

        t_448[k] = f_3 * pc_y[k] * smd_269[k];

        t_449[k] = pb_x[k] * slf0_449[k]
                   - f_4 * pc_x[k] * slf1_449[k];
    }

#pragma omp simd aligned(t_450, t_451, t_452, t_453, t_454, t_455, pc_x, pc_y, pc_z, sld_216, \
                         smp0_135, smp1_135, smd_270, smd_273, smd_274, \
                         smd_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_450[k] = f_1 * smp0_135[k]
                   - f_2 * smp1_135[k]
                   + f_3 * pc_x[k] * smd_270[k];

        t_451[k] = f_0 * sld_216[k]
                   + f_3 * pc_y[k] * smd_270[k];

        t_452[k] = f_3 * pc_z[k] * smd_270[k];

        t_453[k] = f_3 * pc_x[k] * smd_273[k];

        t_454[k] = f_3 * pc_x[k] * smd_274[k];

        t_455[k] = f_3 * pc_x[k] * smd_275[k];
    }

#pragma omp simd aligned(t_456, t_457, t_458, t_459, pc_y, pc_z, sld_219, sld_221, smp0_136, \
                         smp0_137, smp1_136, smp1_137, smd_273, \
                         smd_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_456[k] = f_0 * sld_219[k]
                   + f_1 * smp0_136[k]
                   - f_2 * smp1_136[k]
                   + f_3 * pc_y[k] * smd_273[k];

        t_457[k] = f_3 * pc_z[k] * smd_273[k];

        t_458[k] = f_0 * sld_221[k]
                   + f_3 * pc_y[k] * smd_275[k];

        t_459[k] = f_1 * smp0_137[k]
                   - f_2 * smp1_137[k]
                   + f_3 * pc_z[k] * smd_275[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, t_463, t_464, pb_z, pc_x, pc_y, pc_z, slf0_360, \
                         sld_216, sld_222, slf1_360, smd_276, smd_279, \
                         smd_280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = pb_z[k] * slf0_360[k]
                   - f_4 * pc_z[k] * slf1_360[k];

        t_461[k] = f_6 * sld_222[k]
                   + f_3 * pc_y[k] * smd_276[k];

        t_462[k] = f_5 * sld_216[k]
                   + f_3 * pc_z[k] * smd_276[k];

        t_463[k] = f_3 * pc_x[k] * smd_279[k];

        t_464[k] = f_3 * pc_x[k] * smd_280[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, t_468, pb_z, pc_x, pc_y, pc_z, slf0_366, \
                         sld_219, sld_227, slf1_366, smd_279, smd_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = f_3 * pc_x[k] * smd_281[k];

        t_466[k] = pb_z[k] * slf0_366[k]
                   - f_4 * pc_z[k] * slf1_366[k];

        t_467[k] = f_5 * sld_219[k]
                   + f_3 * pc_z[k] * smd_279[k];

        t_468[k] = f_6 * sld_227[k]
                   + f_3 * pc_y[k] * smd_281[k];
    }

#pragma omp simd aligned(t_469, t_470, t_471, t_472, pc_x, pc_y, pc_z, sld_221, sld_222, \
                         sld_228, smp0_140, smp0_141, smp1_140, smp1_141, smd_281, \
                         smd_282 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_469[k] = f_5 * sld_221[k]
                   + f_1 * smp0_140[k]
                   - f_2 * smp1_140[k]
                   + f_3 * pc_z[k] * smd_281[k];

        t_470[k] = f_1 * smp0_141[k]
                   - f_2 * smp1_141[k]
                   + f_3 * pc_x[k] * smd_282[k];

        t_471[k] = f_7 * sld_228[k]
                   + f_3 * pc_y[k] * smd_282[k];

        t_472[k] = f_8 * sld_222[k]
                   + f_3 * pc_z[k] * smd_282[k];
    }

#pragma omp simd aligned(t_473, t_474, t_475, t_476, t_477, pc_x, pc_y, pc_z, sld_225, \
                         sld_231, smp0_142, smp1_142, smd_285, smd_286, \
                         smd_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_473[k] = f_3 * pc_x[k] * smd_285[k];

        t_474[k] = f_3 * pc_x[k] * smd_286[k];

        t_475[k] = f_3 * pc_x[k] * smd_287[k];

        t_476[k] = f_7 * sld_231[k]
                   + f_1 * smp0_142[k]
                   - f_2 * smp1_142[k]
                   + f_3 * pc_y[k] * smd_285[k];

        t_477[k] = f_8 * sld_225[k]
                   + f_3 * pc_z[k] * smd_285[k];
    }

#pragma omp simd aligned(t_478, t_479, t_480, t_481, pc_x, pc_y, pc_z, sld_227, sld_233, \
                         sld_234, smp0_143, smp0_144, smp1_143, smp1_144, smd_287, \
                         smd_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_478[k] = f_7 * sld_233[k]
                   + f_3 * pc_y[k] * smd_287[k];

        t_479[k] = f_8 * sld_227[k]
                   + f_1 * smp0_143[k]
                   - f_2 * smp1_143[k]
                   + f_3 * pc_z[k] * smd_287[k];

        t_480[k] = f_1 * smp0_144[k]
                   - f_2 * smp1_144[k]
                   + f_3 * pc_x[k] * smd_288[k];

        t_481[k] = f_9 * sld_234[k]
                   + f_3 * pc_y[k] * smd_288[k];
    }

#pragma omp simd aligned(t_482, t_483, t_484, t_485, t_486, pc_x, pc_y, pc_z, sld_228, \
                         sld_237, smp0_145, smp1_145, smd_288, smd_291, smd_292, \
                         smd_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_482[k] = f_10 * sld_228[k]
                   + f_3 * pc_z[k] * smd_288[k];

        t_483[k] = f_3 * pc_x[k] * smd_291[k];

        t_484[k] = f_3 * pc_x[k] * smd_292[k];

        t_485[k] = f_3 * pc_x[k] * smd_293[k];

        t_486[k] = f_9 * sld_237[k]
                   + f_1 * smp0_145[k]
                   - f_2 * smp1_145[k]
                   + f_3 * pc_y[k] * smd_291[k];
    }

#pragma omp simd aligned(t_487, t_488, t_489, pc_y, pc_z, sld_231, sld_233, sld_239, smp0_146, \
                         smp1_146, smd_291, smd_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_487[k] = f_10 * sld_231[k]
                   + f_3 * pc_z[k] * smd_291[k];

        t_488[k] = f_9 * sld_239[k]
                   + f_3 * pc_y[k] * smd_293[k];

        t_489[k] = f_10 * sld_233[k]
                   + f_1 * smp0_146[k]
                   - f_2 * smp1_146[k]
                   + f_3 * pc_z[k] * smd_293[k];
    }

#pragma omp simd aligned(t_490, t_491, t_492, t_493, t_494, pc_x, pc_y, pc_z, sld_234, \
                         sld_240, smp0_147, smp1_147, smd_294, smd_297, \
                         smd_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_490[k] = f_1 * smp0_147[k]
                   - f_2 * smp1_147[k]
                   + f_3 * pc_x[k] * smd_294[k];

        t_491[k] = f_11 * sld_240[k]
                   + f_3 * pc_y[k] * smd_294[k];

        t_492[k] = f_12 * sld_234[k]
                   + f_3 * pc_z[k] * smd_294[k];

        t_493[k] = f_3 * pc_x[k] * smd_297[k];

        t_494[k] = f_3 * pc_x[k] * smd_298[k];
    }

#pragma omp simd aligned(t_495, t_496, t_497, t_498, pc_x, pc_y, pc_z, sld_237, sld_243, \
                         sld_245, smp0_148, smp1_148, smd_297, \
                         smd_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_495[k] = f_3 * pc_x[k] * smd_299[k];

        t_496[k] = f_11 * sld_243[k]
                   + f_1 * smp0_148[k]
                   - f_2 * smp1_148[k]
                   + f_3 * pc_y[k] * smd_297[k];

        t_497[k] = f_12 * sld_237[k]
                   + f_3 * pc_z[k] * smd_297[k];

        t_498[k] = f_11 * sld_245[k]
                   + f_3 * pc_y[k] * smd_299[k];
    }

#pragma omp simd aligned(t_499, t_500, t_501, t_502, pc_x, pc_y, pc_z, sld_239, sld_240, \
                         sld_246, smp0_149, smp0_150, smp1_149, smp1_150, smd_299, \
                         smd_300 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_499[k] = f_12 * sld_239[k]
                   + f_1 * smp0_149[k]
                   - f_2 * smp1_149[k]
                   + f_3 * pc_z[k] * smd_299[k];

        t_500[k] = f_1 * smp0_150[k]
                   - f_2 * smp1_150[k]
                   + f_3 * pc_x[k] * smd_300[k];

        t_501[k] = f_12 * sld_246[k]
                   + f_3 * pc_y[k] * smd_300[k];

        t_502[k] = f_11 * sld_240[k]
                   + f_3 * pc_z[k] * smd_300[k];
    }

#pragma omp simd aligned(t_503, t_504, t_505, t_506, t_507, pc_x, pc_y, pc_z, sld_243, \
                         sld_249, smp0_151, smp1_151, smd_303, smd_304, \
                         smd_305 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_503[k] = f_3 * pc_x[k] * smd_303[k];

        t_504[k] = f_3 * pc_x[k] * smd_304[k];

        t_505[k] = f_3 * pc_x[k] * smd_305[k];

        t_506[k] = f_12 * sld_249[k]
                   + f_1 * smp0_151[k]
                   - f_2 * smp1_151[k]
                   + f_3 * pc_y[k] * smd_303[k];

        t_507[k] = f_11 * sld_243[k]
                   + f_3 * pc_z[k] * smd_303[k];
    }

#pragma omp simd aligned(t_508, t_509, t_510, t_511, pc_x, pc_y, pc_z, sld_245, sld_251, \
                         sld_252, smp0_152, smp0_153, smp1_152, smp1_153, smd_305, \
                         smd_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_508[k] = f_12 * sld_251[k]
                   + f_3 * pc_y[k] * smd_305[k];

        t_509[k] = f_11 * sld_245[k]
                   + f_1 * smp0_152[k]
                   - f_2 * smp1_152[k]
                   + f_3 * pc_z[k] * smd_305[k];

        t_510[k] = f_1 * smp0_153[k]
                   - f_2 * smp1_153[k]
                   + f_3 * pc_x[k] * smd_306[k];

        t_511[k] = f_10 * sld_252[k]
                   + f_3 * pc_y[k] * smd_306[k];
    }
}

static auto
compute_prim_smf_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t slf0,
                                                          const size_t sld, const size_t slf1,
                                                          const size_t smp0, const size_t smp1,
                                                          const size_t smd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / q;
    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 4.0 / q;
    const auto f_7 = 3.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 3.0 / q;
    const auto f_10 = 1.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *slf0_440 = buffer.data(slf0 + 440);
    const auto *slf0_446 = buffer.data(slf0 + 446);
    const auto *slf0_449 = buffer.data(slf0 + 449);

    const auto *sld_246 = buffer.data(sld + 246);
    const auto *sld_249 = buffer.data(sld + 249);
    const auto *sld_251 = buffer.data(sld + 251);
    const auto *sld_252 = buffer.data(sld + 252);
    const auto *sld_255 = buffer.data(sld + 255);
    const auto *sld_257 = buffer.data(sld + 257);
    const auto *sld_258 = buffer.data(sld + 258);
    const auto *sld_261 = buffer.data(sld + 261);
    const auto *sld_263 = buffer.data(sld + 263);
    const auto *sld_264 = buffer.data(sld + 264);
    const auto *sld_267 = buffer.data(sld + 267);
    const auto *sld_269 = buffer.data(sld + 269);

    const auto *slf1_440 = buffer.data(slf1 + 440);
    const auto *slf1_446 = buffer.data(slf1 + 446);
    const auto *slf1_449 = buffer.data(slf1 + 449);

    const auto *smp0_154 = buffer.data(smp0 + 154);
    const auto *smp0_155 = buffer.data(smp0 + 155);
    const auto *smp0_156 = buffer.data(smp0 + 156);
    const auto *smp0_157 = buffer.data(smp0 + 157);
    const auto *smp0_158 = buffer.data(smp0 + 158);
    const auto *smp0_162 = buffer.data(smp0 + 162);
    const auto *smp0_163 = buffer.data(smp0 + 163);
    const auto *smp0_164 = buffer.data(smp0 + 164);

    const auto *smp1_154 = buffer.data(smp1 + 154);
    const auto *smp1_155 = buffer.data(smp1 + 155);
    const auto *smp1_156 = buffer.data(smp1 + 156);
    const auto *smp1_157 = buffer.data(smp1 + 157);
    const auto *smp1_158 = buffer.data(smp1 + 158);
    const auto *smp1_162 = buffer.data(smp1 + 162);
    const auto *smp1_163 = buffer.data(smp1 + 163);
    const auto *smp1_164 = buffer.data(smp1 + 164);

    const auto *smd_306 = buffer.data(smd + 306);
    const auto *smd_309 = buffer.data(smd + 309);
    const auto *smd_310 = buffer.data(smd + 310);
    const auto *smd_311 = buffer.data(smd + 311);
    const auto *smd_312 = buffer.data(smd + 312);
    const auto *smd_315 = buffer.data(smd + 315);
    const auto *smd_316 = buffer.data(smd + 316);
    const auto *smd_317 = buffer.data(smd + 317);
    const auto *smd_318 = buffer.data(smd + 318);
    const auto *smd_321 = buffer.data(smd + 321);
    const auto *smd_322 = buffer.data(smd + 322);
    const auto *smd_323 = buffer.data(smd + 323);
    const auto *smd_324 = buffer.data(smd + 324);
    const auto *smd_327 = buffer.data(smd + 327);
    const auto *smd_328 = buffer.data(smd + 328);
    const auto *smd_329 = buffer.data(smd + 329);

#pragma omp simd aligned(t_512, t_513, t_514, t_515, t_516, pc_x, pc_y, pc_z, sld_246, \
                         sld_255, smp0_154, smp1_154, smd_306, smd_309, smd_310, \
                         smd_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = f_9 * sld_246[k]
                   + f_3 * pc_z[k] * smd_306[k];

        t_513[k] = f_3 * pc_x[k] * smd_309[k];

        t_514[k] = f_3 * pc_x[k] * smd_310[k];

        t_515[k] = f_3 * pc_x[k] * smd_311[k];

        t_516[k] = f_10 * sld_255[k]
                   + f_1 * smp0_154[k]
                   - f_2 * smp1_154[k]
                   + f_3 * pc_y[k] * smd_309[k];
    }

#pragma omp simd aligned(t_517, t_518, t_519, pc_y, pc_z, sld_249, sld_251, sld_257, smp0_155, \
                         smp1_155, smd_309, smd_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_517[k] = f_9 * sld_249[k]
                   + f_3 * pc_z[k] * smd_309[k];

        t_518[k] = f_10 * sld_257[k]
                   + f_3 * pc_y[k] * smd_311[k];

        t_519[k] = f_9 * sld_251[k]
                   + f_1 * smp0_155[k]
                   - f_2 * smp1_155[k]
                   + f_3 * pc_z[k] * smd_311[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, t_523, t_524, pc_x, pc_y, pc_z, sld_252, \
                         sld_258, smp0_156, smp1_156, smd_312, smd_315, \
                         smd_316 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = f_1 * smp0_156[k]
                   - f_2 * smp1_156[k]
                   + f_3 * pc_x[k] * smd_312[k];

        t_521[k] = f_8 * sld_258[k]
                   + f_3 * pc_y[k] * smd_312[k];

        t_522[k] = f_7 * sld_252[k]
                   + f_3 * pc_z[k] * smd_312[k];

        t_523[k] = f_3 * pc_x[k] * smd_315[k];

        t_524[k] = f_3 * pc_x[k] * smd_316[k];
    }

#pragma omp simd aligned(t_525, t_526, t_527, t_528, pc_x, pc_y, pc_z, sld_255, sld_261, \
                         sld_263, smp0_157, smp1_157, smd_315, \
                         smd_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_525[k] = f_3 * pc_x[k] * smd_317[k];

        t_526[k] = f_8 * sld_261[k]
                   + f_1 * smp0_157[k]
                   - f_2 * smp1_157[k]
                   + f_3 * pc_y[k] * smd_315[k];

        t_527[k] = f_7 * sld_255[k]
                   + f_3 * pc_z[k] * smd_315[k];

        t_528[k] = f_8 * sld_263[k]
                   + f_3 * pc_y[k] * smd_317[k];
    }

#pragma omp simd aligned(t_529, t_530, t_531, t_532, pb_y, pc_y, pc_z, slf0_440, sld_257, \
                         sld_258, sld_264, slf1_440, smp0_158, smp1_158, smd_317, \
                         smd_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_529[k] = f_7 * sld_257[k]
                   + f_1 * smp0_158[k]
                   - f_2 * smp1_158[k]
                   + f_3 * pc_z[k] * smd_317[k];

        t_530[k] = pb_y[k] * slf0_440[k]
                   - f_4 * pc_y[k] * slf1_440[k];

        t_531[k] = f_5 * sld_264[k]
                   + f_3 * pc_y[k] * smd_318[k];

        t_532[k] = f_6 * sld_258[k]
                   + f_3 * pc_z[k] * smd_318[k];
    }

#pragma omp simd aligned(t_533, t_534, t_535, t_536, t_537, pb_y, pc_x, pc_y, pc_z, slf0_446, \
                         sld_261, sld_267, slf1_446, smd_321, smd_322, \
                         smd_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_533[k] = f_3 * pc_x[k] * smd_321[k];

        t_534[k] = f_3 * pc_x[k] * smd_322[k];

        t_535[k] = f_3 * pc_x[k] * smd_323[k];

        t_536[k] = pb_y[k] * slf0_446[k]
                   + f_10 * sld_267[k]
                   - f_4 * pc_y[k] * slf1_446[k];

        t_537[k] = f_6 * sld_261[k]
                   + f_3 * pc_z[k] * smd_321[k];
    }

#pragma omp simd aligned(t_538, t_539, t_540, t_541, pb_y, pc_x, pc_y, slf0_449, sld_269, \
                         slf1_449, smp0_162, smp1_162, smd_323, \
                         smd_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_538[k] = f_5 * sld_269[k]
                   + f_3 * pc_y[k] * smd_323[k];

        t_539[k] = pb_y[k] * slf0_449[k]
                   - f_4 * pc_y[k] * slf1_449[k];

        t_540[k] = f_1 * smp0_162[k]
                   - f_2 * smp1_162[k]
                   + f_3 * pc_x[k] * smd_324[k];

        t_541[k] = f_3 * pc_y[k] * smd_324[k];
    }

#pragma omp simd aligned(t_542, t_543, t_544, t_545, t_546, pc_x, pc_y, pc_z, sld_264, \
                         smp0_163, smp1_163, smd_324, smd_327, smd_328, \
                         smd_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_542[k] = f_0 * sld_264[k]
                   + f_3 * pc_z[k] * smd_324[k];

        t_543[k] = f_3 * pc_x[k] * smd_327[k];

        t_544[k] = f_3 * pc_x[k] * smd_328[k];

        t_545[k] = f_3 * pc_x[k] * smd_329[k];

        t_546[k] = f_1 * smp0_163[k]
                   - f_2 * smp1_163[k]
                   + f_3 * pc_y[k] * smd_327[k];
    }

#pragma omp simd aligned(t_547, t_548, t_549, pc_y, pc_z, sld_267, sld_269, smp0_164, \
                         smp1_164, smd_327, smd_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_547[k] = f_0 * sld_267[k]
                   + f_3 * pc_z[k] * smd_327[k];

        t_548[k] = f_3 * pc_y[k] * smd_329[k];

        t_549[k] = f_0 * sld_269[k]
                   + f_1 * smp0_164[k]
                   - f_2 * smp1_164[k]
                   + f_3 * pc_z[k] * smd_329[k];
    }
}

auto
compute_prim_smf_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t slf0, const size_t sld,
                                                   const size_t slf1, const size_t smp0,
                                                   const size_t smp1, const size_t smd,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_smf_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, slf0, sld,
                                                              slf1, smp0, smp1, smd, ncols,
                                                              gamma, p, q);

    compute_prim_smf_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, slf0, sld,
                                                              slf1, smp0, smp1, smd, ncols,
                                                              gamma, p, q);

    compute_prim_smf_three_center_electron_repulsion_0_piece2(buffer, target, pb, pc, slf0, sld,
                                                              slf1, smp0, smp1, smd, ncols,
                                                              gamma, p, q);

    compute_prim_smf_three_center_electron_repulsion_0_piece3(buffer, target, pb, pc, slf0, sld,
                                                              slf1, smp0, smp1, smd, ncols,
                                                              gamma, p, q);

    compute_prim_smf_three_center_electron_repulsion_0_piece4(buffer, target, pb, pc, slf0, sld,
                                                              slf1, smp0, smp1, smd, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
