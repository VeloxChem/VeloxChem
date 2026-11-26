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


#include "SimdThreeCenterElectronRepulsionVrrRecSKF.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_skf_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sif0,
                                                          const size_t sid, const size_t sif1,
                                                          const size_t skp0, const size_t skp1,
                                                          const size_t skd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / q;
    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 3.0 / q;
    const auto f_7 = 2.5 / q;
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

    const auto *sif0_0 = buffer.data(sif0 + 0);
    const auto *sif0_6 = buffer.data(sif0 + 6);
    const auto *sif0_9 = buffer.data(sif0 + 9);
    const auto *sif0_16 = buffer.data(sif0 + 16);
    const auto *sif0_20 = buffer.data(sif0 + 20);
    const auto *sif0_29 = buffer.data(sif0 + 29);
    const auto *sif0_30 = buffer.data(sif0 + 30);
    const auto *sif0_36 = buffer.data(sif0 + 36);
    const auto *sif0_50 = buffer.data(sif0 + 50);
    const auto *sif0_59 = buffer.data(sif0 + 59);
    const auto *sif0_60 = buffer.data(sif0 + 60);
    const auto *sif0_66 = buffer.data(sif0 + 66);

    const auto *sid_0 = buffer.data(sid + 0);
    const auto *sid_3 = buffer.data(sid + 3);
    const auto *sid_4 = buffer.data(sid + 4);
    const auto *sid_5 = buffer.data(sid + 5);
    const auto *sid_6 = buffer.data(sid + 6);
    const auto *sid_9 = buffer.data(sid + 9);
    const auto *sid_10 = buffer.data(sid + 10);
    const auto *sid_11 = buffer.data(sid + 11);
    const auto *sid_12 = buffer.data(sid + 12);
    const auto *sid_15 = buffer.data(sid + 15);
    const auto *sid_16 = buffer.data(sid + 16);
    const auto *sid_17 = buffer.data(sid + 17);
    const auto *sid_18 = buffer.data(sid + 18);
    const auto *sid_21 = buffer.data(sid + 21);
    const auto *sid_22 = buffer.data(sid + 22);
    const auto *sid_23 = buffer.data(sid + 23);
    const auto *sid_24 = buffer.data(sid + 24);
    const auto *sid_27 = buffer.data(sid + 27);
    const auto *sid_28 = buffer.data(sid + 28);
    const auto *sid_29 = buffer.data(sid + 29);
    const auto *sid_30 = buffer.data(sid + 30);
    const auto *sid_33 = buffer.data(sid + 33);
    const auto *sid_34 = buffer.data(sid + 34);
    const auto *sid_35 = buffer.data(sid + 35);
    const auto *sid_36 = buffer.data(sid + 36);
    const auto *sid_39 = buffer.data(sid + 39);
    const auto *sid_40 = buffer.data(sid + 40);
    const auto *sid_41 = buffer.data(sid + 41);
    const auto *sid_42 = buffer.data(sid + 42);
    const auto *sid_45 = buffer.data(sid + 45);
    const auto *sid_46 = buffer.data(sid + 46);
    const auto *sid_47 = buffer.data(sid + 47);
    const auto *sid_48 = buffer.data(sid + 48);
    const auto *sid_51 = buffer.data(sid + 51);
    const auto *sid_52 = buffer.data(sid + 52);
    const auto *sid_53 = buffer.data(sid + 53);
    const auto *sid_54 = buffer.data(sid + 54);
    const auto *sid_57 = buffer.data(sid + 57);
    const auto *sid_58 = buffer.data(sid + 58);
    const auto *sid_59 = buffer.data(sid + 59);
    const auto *sid_60 = buffer.data(sid + 60);
    const auto *sid_63 = buffer.data(sid + 63);
    const auto *sid_64 = buffer.data(sid + 64);
    const auto *sid_65 = buffer.data(sid + 65);
    const auto *sid_69 = buffer.data(sid + 69);
    const auto *sid_70 = buffer.data(sid + 70);
    const auto *sid_71 = buffer.data(sid + 71);
    const auto *sid_72 = buffer.data(sid + 72);
    const auto *sid_75 = buffer.data(sid + 75);
    const auto *sid_76 = buffer.data(sid + 76);
    const auto *sid_77 = buffer.data(sid + 77);

    const auto *sif1_0 = buffer.data(sif1 + 0);
    const auto *sif1_6 = buffer.data(sif1 + 6);
    const auto *sif1_9 = buffer.data(sif1 + 9);
    const auto *sif1_16 = buffer.data(sif1 + 16);
    const auto *sif1_20 = buffer.data(sif1 + 20);
    const auto *sif1_29 = buffer.data(sif1 + 29);
    const auto *sif1_30 = buffer.data(sif1 + 30);
    const auto *sif1_36 = buffer.data(sif1 + 36);
    const auto *sif1_50 = buffer.data(sif1 + 50);
    const auto *sif1_59 = buffer.data(sif1 + 59);
    const auto *sif1_60 = buffer.data(sif1 + 60);
    const auto *sif1_66 = buffer.data(sif1 + 66);

    const auto *skp0_0 = buffer.data(skp0 + 0);
    const auto *skp0_1 = buffer.data(skp0 + 1);
    const auto *skp0_2 = buffer.data(skp0 + 2);
    const auto *skp0_4 = buffer.data(skp0 + 4);
    const auto *skp0_8 = buffer.data(skp0 + 8);
    const auto *skp0_9 = buffer.data(skp0 + 9);
    const auto *skp0_10 = buffer.data(skp0 + 10);
    const auto *skp0_11 = buffer.data(skp0 + 11);
    const auto *skp0_15 = buffer.data(skp0 + 15);
    const auto *skp0_16 = buffer.data(skp0 + 16);
    const auto *skp0_17 = buffer.data(skp0 + 17);
    const auto *skp0_18 = buffer.data(skp0 + 18);
    const auto *skp0_19 = buffer.data(skp0 + 19);
    const auto *skp0_20 = buffer.data(skp0 + 20);
    const auto *skp0_23 = buffer.data(skp0 + 23);
    const auto *skp0_25 = buffer.data(skp0 + 25);
    const auto *skp0_27 = buffer.data(skp0 + 27);
    const auto *skp0_28 = buffer.data(skp0 + 28);
    const auto *skp0_29 = buffer.data(skp0 + 29);
    const auto *skp0_30 = buffer.data(skp0 + 30);
    const auto *skp0_31 = buffer.data(skp0 + 31);
    const auto *skp0_32 = buffer.data(skp0 + 32);
    const auto *skp0_35 = buffer.data(skp0 + 35);
    const auto *skp0_36 = buffer.data(skp0 + 36);
    const auto *skp0_37 = buffer.data(skp0 + 37);
    const auto *skp0_38 = buffer.data(skp0 + 38);

    const auto *skp1_0 = buffer.data(skp1 + 0);
    const auto *skp1_1 = buffer.data(skp1 + 1);
    const auto *skp1_2 = buffer.data(skp1 + 2);
    const auto *skp1_4 = buffer.data(skp1 + 4);
    const auto *skp1_8 = buffer.data(skp1 + 8);
    const auto *skp1_9 = buffer.data(skp1 + 9);
    const auto *skp1_10 = buffer.data(skp1 + 10);
    const auto *skp1_11 = buffer.data(skp1 + 11);
    const auto *skp1_15 = buffer.data(skp1 + 15);
    const auto *skp1_16 = buffer.data(skp1 + 16);
    const auto *skp1_17 = buffer.data(skp1 + 17);
    const auto *skp1_18 = buffer.data(skp1 + 18);
    const auto *skp1_19 = buffer.data(skp1 + 19);
    const auto *skp1_20 = buffer.data(skp1 + 20);
    const auto *skp1_23 = buffer.data(skp1 + 23);
    const auto *skp1_25 = buffer.data(skp1 + 25);
    const auto *skp1_27 = buffer.data(skp1 + 27);
    const auto *skp1_28 = buffer.data(skp1 + 28);
    const auto *skp1_29 = buffer.data(skp1 + 29);
    const auto *skp1_30 = buffer.data(skp1 + 30);
    const auto *skp1_31 = buffer.data(skp1 + 31);
    const auto *skp1_32 = buffer.data(skp1 + 32);
    const auto *skp1_35 = buffer.data(skp1 + 35);
    const auto *skp1_36 = buffer.data(skp1 + 36);
    const auto *skp1_37 = buffer.data(skp1 + 37);
    const auto *skp1_38 = buffer.data(skp1 + 38);

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
    const auto *skd_66 = buffer.data(skd + 66);
    const auto *skd_69 = buffer.data(skd + 69);
    const auto *skd_70 = buffer.data(skd + 70);
    const auto *skd_71 = buffer.data(skd + 71);
    const auto *skd_72 = buffer.data(skd + 72);
    const auto *skd_75 = buffer.data(skd + 75);
    const auto *skd_76 = buffer.data(skd + 76);
    const auto *skd_77 = buffer.data(skd + 77);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pc_x, pc_y, pc_z, sid_0, sid_3, sid_4, \
                         skp0_0, skp1_0, skd_0, skd_3, skd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sid_0[k]
                 + f_1 * skp0_0[k]
                 - f_2 * skp1_0[k]
                 + f_3 * pc_x[k] * skd_0[k];

        t_1[k] = f_3 * pc_y[k] * skd_0[k];

        t_2[k] = f_3 * pc_z[k] * skd_0[k];

        t_3[k] = f_0 * sid_3[k]
                 + f_3 * pc_x[k] * skd_3[k];

        t_4[k] = f_0 * sid_4[k]
                 + f_3 * pc_x[k] * skd_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pc_x, pc_y, pc_z, sid_5, skp0_1, skp0_2, \
                         skp1_1, skp1_2, skd_3, skd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * sid_5[k]
                 + f_3 * pc_x[k] * skd_5[k];

        t_6[k] = f_1 * skp0_1[k]
                 - f_2 * skp1_1[k]
                 + f_3 * pc_y[k] * skd_3[k];

        t_7[k] = f_3 * pc_z[k] * skd_3[k];

        t_8[k] = f_3 * pc_y[k] * skd_5[k];

        t_9[k] = f_1 * skp0_2[k]
                 - f_2 * skp1_2[k]
                 + f_3 * pc_z[k] * skd_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pb_y, pc_x, pc_y, pc_z, sif0_0, sid_0, sid_9, \
                         sif1_0, skd_6, skd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pb_y[k] * sif0_0[k]
                  - f_4 * pc_y[k] * sif1_0[k];

        t_11[k] = f_5 * sid_0[k]
                  + f_3 * pc_y[k] * skd_6[k];

        t_12[k] = f_3 * pc_z[k] * skd_6[k];

        t_13[k] = f_6 * sid_9[k]
                  + f_3 * pc_x[k] * skd_9[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pc_x, pc_y, pc_z, sid_3, sid_10, sid_11, \
                         skp0_4, skp1_4, skd_9, skd_10, skd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_6 * sid_10[k]
                  + f_3 * pc_x[k] * skd_10[k];

        t_15[k] = f_6 * sid_11[k]
                  + f_3 * pc_x[k] * skd_11[k];

        t_16[k] = f_5 * sid_3[k]
                  + f_1 * skp0_4[k]
                  - f_2 * skp1_4[k]
                  + f_3 * pc_y[k] * skd_9[k];

        t_17[k] = f_3 * pc_z[k] * skd_9[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pb_y, pb_z, pc_y, pc_z, sif0_0, sif0_9, \
                         sid_5, sif1_0, sif1_9, skd_11, skd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * sid_5[k]
                  + f_3 * pc_y[k] * skd_11[k];

        t_19[k] = pb_y[k] * sif0_9[k]
                  - f_4 * pc_y[k] * sif1_9[k];

        t_20[k] = pb_z[k] * sif0_0[k]
                  - f_4 * pc_z[k] * sif1_0[k];

        t_21[k] = f_3 * pc_y[k] * skd_12[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pc_x, pc_z, sid_0, sid_15, sid_16, sid_17, \
                         skd_12, skd_15, skd_16, skd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_5 * sid_0[k]
                  + f_3 * pc_z[k] * skd_12[k];

        t_23[k] = f_6 * sid_15[k]
                  + f_3 * pc_x[k] * skd_15[k];

        t_24[k] = f_6 * sid_16[k]
                  + f_3 * pc_x[k] * skd_16[k];

        t_25[k] = f_6 * sid_17[k]
                  + f_3 * pc_x[k] * skd_17[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pb_z, pc_y, pc_z, sif0_6, sid_3, sid_5, \
                         sif1_6, skp0_8, skp1_8, skd_15, skd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pb_z[k] * sif0_6[k]
                  - f_4 * pc_z[k] * sif1_6[k];

        t_27[k] = f_5 * sid_3[k]
                  + f_3 * pc_z[k] * skd_15[k];

        t_28[k] = f_3 * pc_y[k] * skd_17[k];

        t_29[k] = f_5 * sid_5[k]
                  + f_1 * skp0_8[k]
                  - f_2 * skp1_8[k]
                  + f_3 * pc_z[k] * skd_17[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pc_x, pc_y, pc_z, sid_6, sid_18, sid_21, \
                         skp0_9, skp1_9, skd_18, skd_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_7 * sid_18[k]
                  + f_1 * skp0_9[k]
                  - f_2 * skp1_9[k]
                  + f_3 * pc_x[k] * skd_18[k];

        t_31[k] = f_8 * sid_6[k]
                  + f_3 * pc_y[k] * skd_18[k];

        t_32[k] = f_3 * pc_z[k] * skd_18[k];

        t_33[k] = f_7 * sid_21[k]
                  + f_3 * pc_x[k] * skd_21[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pc_x, pc_y, pc_z, sid_9, sid_22, sid_23, \
                         skp0_10, skp1_10, skd_21, skd_22, skd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_7 * sid_22[k]
                  + f_3 * pc_x[k] * skd_22[k];

        t_35[k] = f_7 * sid_23[k]
                  + f_3 * pc_x[k] * skd_23[k];

        t_36[k] = f_8 * sid_9[k]
                  + f_1 * skp0_10[k]
                  - f_2 * skp1_10[k]
                  + f_3 * pc_y[k] * skd_21[k];

        t_37[k] = f_3 * pc_z[k] * skd_21[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pb_y, pc_y, pc_z, sif0_20, sid_11, sid_12, \
                         sif1_20, skp0_11, skp1_11, skd_23, skd_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_8 * sid_11[k]
                  + f_3 * pc_y[k] * skd_23[k];

        t_39[k] = f_1 * skp0_11[k]
                  - f_2 * skp1_11[k]
                  + f_3 * pc_z[k] * skd_23[k];

        t_40[k] = pb_y[k] * sif0_20[k]
                  - f_4 * pc_y[k] * sif1_20[k];

        t_41[k] = f_5 * sid_12[k]
                  + f_3 * pc_y[k] * skd_24[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pc_x, pc_z, sid_6, sid_27, sid_28, sid_29, \
                         skd_24, skd_27, skd_28, skd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_5 * sid_6[k]
                  + f_3 * pc_z[k] * skd_24[k];

        t_43[k] = f_7 * sid_27[k]
                  + f_3 * pc_x[k] * skd_27[k];

        t_44[k] = f_7 * sid_28[k]
                  + f_3 * pc_x[k] * skd_28[k];

        t_45[k] = f_7 * sid_29[k]
                  + f_3 * pc_x[k] * skd_29[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pb_y, pb_z, pc_y, pc_z, sif0_16, sif0_29, \
                         sid_9, sid_17, sif1_16, sif1_29, skd_27, \
                         skd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = pb_z[k] * sif0_16[k]
                  - f_4 * pc_z[k] * sif1_16[k];

        t_47[k] = f_5 * sid_9[k]
                  + f_3 * pc_z[k] * skd_27[k];

        t_48[k] = f_5 * sid_17[k]
                  + f_3 * pc_y[k] * skd_29[k];

        t_49[k] = pb_y[k] * sif0_29[k]
                  - f_4 * pc_y[k] * sif1_29[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pc_x, pc_y, pc_z, sid_12, sid_30, sid_33, \
                         skp0_15, skp1_15, skd_30, skd_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_7 * sid_30[k]
                  + f_1 * skp0_15[k]
                  - f_2 * skp1_15[k]
                  + f_3 * pc_x[k] * skd_30[k];

        t_51[k] = f_3 * pc_y[k] * skd_30[k];

        t_52[k] = f_8 * sid_12[k]
                  + f_3 * pc_z[k] * skd_30[k];

        t_53[k] = f_7 * sid_33[k]
                  + f_3 * pc_x[k] * skd_33[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, pc_x, pc_y, pc_z, sid_15, sid_34, \
                         sid_35, skp0_16, skp1_16, skd_33, skd_34, \
                         skd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_7 * sid_34[k]
                  + f_3 * pc_x[k] * skd_34[k];

        t_55[k] = f_7 * sid_35[k]
                  + f_3 * pc_x[k] * skd_35[k];

        t_56[k] = f_1 * skp0_16[k]
                  - f_2 * skp1_16[k]
                  + f_3 * pc_y[k] * skd_33[k];

        t_57[k] = f_8 * sid_15[k]
                  + f_3 * pc_z[k] * skd_33[k];

        t_58[k] = f_3 * pc_y[k] * skd_35[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pc_x, pc_y, pc_z, sid_17, sid_18, sid_36, \
                         skp0_17, skp0_18, skp1_17, skp1_18, skd_35, \
                         skd_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_8 * sid_17[k]
                  + f_1 * skp0_17[k]
                  - f_2 * skp1_17[k]
                  + f_3 * pc_z[k] * skd_35[k];

        t_60[k] = f_9 * sid_36[k]
                  + f_1 * skp0_18[k]
                  - f_2 * skp1_18[k]
                  + f_3 * pc_x[k] * skd_36[k];

        t_61[k] = f_10 * sid_18[k]
                  + f_3 * pc_y[k] * skd_36[k];

        t_62[k] = f_3 * pc_z[k] * skd_36[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pc_x, pc_y, sid_21, sid_39, sid_40, sid_41, \
                         skp0_19, skp1_19, skd_39, skd_40, skd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_9 * sid_39[k]
                  + f_3 * pc_x[k] * skd_39[k];

        t_64[k] = f_9 * sid_40[k]
                  + f_3 * pc_x[k] * skd_40[k];

        t_65[k] = f_9 * sid_41[k]
                  + f_3 * pc_x[k] * skd_41[k];

        t_66[k] = f_10 * sid_21[k]
                  + f_1 * skp0_19[k]
                  - f_2 * skp1_19[k]
                  + f_3 * pc_y[k] * skd_39[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, pb_z, pc_y, pc_z, sif0_30, sid_23, sif1_30, \
                         skp0_20, skp1_20, skd_39, skd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_3 * pc_z[k] * skd_39[k];

        t_68[k] = f_10 * sid_23[k]
                  + f_3 * pc_y[k] * skd_41[k];

        t_69[k] = f_1 * skp0_20[k]
                  - f_2 * skp1_20[k]
                  + f_3 * pc_z[k] * skd_41[k];

        t_70[k] = pb_z[k] * sif0_30[k]
                  - f_4 * pc_z[k] * sif1_30[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pc_x, pc_y, pc_z, sid_18, sid_24, sid_45, \
                         sid_46, skd_42, skd_45, skd_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_8 * sid_24[k]
                  + f_3 * pc_y[k] * skd_42[k];

        t_72[k] = f_5 * sid_18[k]
                  + f_3 * pc_z[k] * skd_42[k];

        t_73[k] = f_9 * sid_45[k]
                  + f_3 * pc_x[k] * skd_45[k];

        t_74[k] = f_9 * sid_46[k]
                  + f_3 * pc_x[k] * skd_46[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pb_z, pc_x, pc_y, pc_z, sif0_36, sid_21, \
                         sid_29, sid_47, sif1_36, skd_45, skd_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_9 * sid_47[k]
                  + f_3 * pc_x[k] * skd_47[k];

        t_76[k] = pb_z[k] * sif0_36[k]
                  - f_4 * pc_z[k] * sif1_36[k];

        t_77[k] = f_5 * sid_21[k]
                  + f_3 * pc_z[k] * skd_45[k];

        t_78[k] = f_8 * sid_29[k]
                  + f_3 * pc_y[k] * skd_47[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pb_y, pc_y, pc_z, sif0_50, sid_23, sid_24, \
                         sid_30, sif1_50, skp0_23, skp1_23, skd_47, \
                         skd_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_5 * sid_23[k]
                  + f_1 * skp0_23[k]
                  - f_2 * skp1_23[k]
                  + f_3 * pc_z[k] * skd_47[k];

        t_80[k] = pb_y[k] * sif0_50[k]
                  - f_4 * pc_y[k] * sif1_50[k];

        t_81[k] = f_5 * sid_30[k]
                  + f_3 * pc_y[k] * skd_48[k];

        t_82[k] = f_8 * sid_24[k]
                  + f_3 * pc_z[k] * skd_48[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pc_x, pc_y, sid_33, sid_51, sid_52, sid_53, \
                         skp0_25, skp1_25, skd_51, skd_52, skd_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_9 * sid_51[k]
                  + f_3 * pc_x[k] * skd_51[k];

        t_84[k] = f_9 * sid_52[k]
                  + f_3 * pc_x[k] * skd_52[k];

        t_85[k] = f_9 * sid_53[k]
                  + f_3 * pc_x[k] * skd_53[k];

        t_86[k] = f_5 * sid_33[k]
                  + f_1 * skp0_25[k]
                  - f_2 * skp1_25[k]
                  + f_3 * pc_y[k] * skd_51[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, pb_y, pc_y, pc_z, sif0_59, sid_27, sid_35, sif1_59, \
                         skd_51, skd_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_8 * sid_27[k]
                  + f_3 * pc_z[k] * skd_51[k];

        t_88[k] = f_5 * sid_35[k]
                  + f_3 * pc_y[k] * skd_53[k];

        t_89[k] = pb_y[k] * sif0_59[k]
                  - f_4 * pc_y[k] * sif1_59[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, pc_x, pc_y, pc_z, sid_30, sid_54, sid_57, \
                         skp0_27, skp1_27, skd_54, skd_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_9 * sid_54[k]
                  + f_1 * skp0_27[k]
                  - f_2 * skp1_27[k]
                  + f_3 * pc_x[k] * skd_54[k];

        t_91[k] = f_3 * pc_y[k] * skd_54[k];

        t_92[k] = f_10 * sid_30[k]
                  + f_3 * pc_z[k] * skd_54[k];

        t_93[k] = f_9 * sid_57[k]
                  + f_3 * pc_x[k] * skd_57[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, pc_x, pc_y, pc_z, sid_33, sid_58, \
                         sid_59, skp0_28, skp1_28, skd_57, skd_58, \
                         skd_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_9 * sid_58[k]
                  + f_3 * pc_x[k] * skd_58[k];

        t_95[k] = f_9 * sid_59[k]
                  + f_3 * pc_x[k] * skd_59[k];

        t_96[k] = f_1 * skp0_28[k]
                  - f_2 * skp1_28[k]
                  + f_3 * pc_y[k] * skd_57[k];

        t_97[k] = f_10 * sid_33[k]
                  + f_3 * pc_z[k] * skd_57[k];

        t_98[k] = f_3 * pc_y[k] * skd_59[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, pc_x, pc_y, pc_z, sid_35, sid_36, sid_60, \
                         skp0_29, skp0_30, skp1_29, skp1_30, skd_59, \
                         skd_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_10 * sid_35[k]
                  + f_1 * skp0_29[k]
                  - f_2 * skp1_29[k]
                  + f_3 * pc_z[k] * skd_59[k];

        t_100[k] = f_10 * sid_60[k]
                   + f_1 * skp0_30[k]
                   - f_2 * skp1_30[k]
                   + f_3 * pc_x[k] * skd_60[k];

        t_101[k] = f_9 * sid_36[k]
                   + f_3 * pc_y[k] * skd_60[k];

        t_102[k] = f_3 * pc_z[k] * skd_60[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, pc_x, pc_y, sid_39, sid_63, sid_64, \
                         sid_65, skp0_31, skp1_31, skd_63, skd_64, \
                         skd_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_10 * sid_63[k]
                   + f_3 * pc_x[k] * skd_63[k];

        t_104[k] = f_10 * sid_64[k]
                   + f_3 * pc_x[k] * skd_64[k];

        t_105[k] = f_10 * sid_65[k]
                   + f_3 * pc_x[k] * skd_65[k];

        t_106[k] = f_9 * sid_39[k]
                   + f_1 * skp0_31[k]
                   - f_2 * skp1_31[k]
                   + f_3 * pc_y[k] * skd_63[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, pb_z, pc_y, pc_z, sif0_60, sid_41, \
                         sif1_60, skp0_32, skp1_32, skd_63, skd_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_3 * pc_z[k] * skd_63[k];

        t_108[k] = f_9 * sid_41[k]
                   + f_3 * pc_y[k] * skd_65[k];

        t_109[k] = f_1 * skp0_32[k]
                   - f_2 * skp1_32[k]
                   + f_3 * pc_z[k] * skd_65[k];

        t_110[k] = pb_z[k] * sif0_60[k]
                   - f_4 * pc_z[k] * sif1_60[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, pc_x, pc_y, pc_z, sid_36, sid_42, sid_69, \
                         sid_70, skd_66, skd_69, skd_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_10 * sid_42[k]
                   + f_3 * pc_y[k] * skd_66[k];

        t_112[k] = f_5 * sid_36[k]
                   + f_3 * pc_z[k] * skd_66[k];

        t_113[k] = f_10 * sid_69[k]
                   + f_3 * pc_x[k] * skd_69[k];

        t_114[k] = f_10 * sid_70[k]
                   + f_3 * pc_x[k] * skd_70[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, pb_z, pc_x, pc_y, pc_z, sif0_66, sid_39, \
                         sid_47, sid_71, sif1_66, skd_69, skd_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_10 * sid_71[k]
                   + f_3 * pc_x[k] * skd_71[k];

        t_116[k] = pb_z[k] * sif0_66[k]
                   - f_4 * pc_z[k] * sif1_66[k];

        t_117[k] = f_5 * sid_39[k]
                   + f_3 * pc_z[k] * skd_69[k];

        t_118[k] = f_10 * sid_47[k]
                   + f_3 * pc_y[k] * skd_71[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, pc_x, pc_y, pc_z, sid_41, sid_48, sid_72, \
                         skp0_35, skp0_36, skp1_35, skp1_36, skd_71, \
                         skd_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = f_5 * sid_41[k]
                   + f_1 * skp0_35[k]
                   - f_2 * skp1_35[k]
                   + f_3 * pc_z[k] * skd_71[k];

        t_120[k] = f_10 * sid_72[k]
                   + f_1 * skp0_36[k]
                   - f_2 * skp1_36[k]
                   + f_3 * pc_x[k] * skd_72[k];

        t_121[k] = f_8 * sid_48[k]
                   + f_3 * pc_y[k] * skd_72[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, pc_x, pc_z, sid_42, sid_75, sid_76, \
                         sid_77, skd_72, skd_75, skd_76, skd_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_8 * sid_42[k]
                   + f_3 * pc_z[k] * skd_72[k];

        t_123[k] = f_10 * sid_75[k]
                   + f_3 * pc_x[k] * skd_75[k];

        t_124[k] = f_10 * sid_76[k]
                   + f_3 * pc_x[k] * skd_76[k];

        t_125[k] = f_10 * sid_77[k]
                   + f_3 * pc_x[k] * skd_77[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, pc_y, pc_z, sid_45, sid_47, sid_51, \
                         sid_53, skp0_37, skp0_38, skp1_37, skp1_38, skd_75, \
                         skd_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_8 * sid_51[k]
                   + f_1 * skp0_37[k]
                   - f_2 * skp1_37[k]
                   + f_3 * pc_y[k] * skd_75[k];

        t_127[k] = f_8 * sid_45[k]
                   + f_3 * pc_z[k] * skd_75[k];

        t_128[k] = f_8 * sid_53[k]
                   + f_3 * pc_y[k] * skd_77[k];

        t_129[k] = f_8 * sid_47[k]
                   + f_1 * skp0_38[k]
                   - f_2 * skp1_38[k]
                   + f_3 * pc_z[k] * skd_77[k];
    }
}

static auto
compute_prim_skf_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sif0,
                                                          const size_t sid, const size_t sif1,
                                                          const size_t skp0, const size_t skp1,
                                                          const size_t skd, const size_t ncols,
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
    const auto f_6 = 3.0 / q;
    const auto f_7 = 2.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 2.0 / q;
    const auto f_10 = 1.5 / q;

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
    auto *t_255 = buffer.data(target + 255);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sif0_90 = buffer.data(sif0 + 90);
    const auto *sif0_99 = buffer.data(sif0 + 99);
    const auto *sif0_100 = buffer.data(sif0 + 100);
    const auto *sif0_106 = buffer.data(sif0 + 106);
    const auto *sif0_140 = buffer.data(sif0 + 140);
    const auto *sif0_149 = buffer.data(sif0 + 149);
    const auto *sif0_150 = buffer.data(sif0 + 150);
    const auto *sif0_210 = buffer.data(sif0 + 210);
    const auto *sif0_216 = buffer.data(sif0 + 216);
    const auto *sif0_219 = buffer.data(sif0 + 219);
    const auto *sif0_226 = buffer.data(sif0 + 226);
    const auto *sif0_229 = buffer.data(sif0 + 229);
    const auto *sif0_230 = buffer.data(sif0 + 230);
    const auto *sif0_236 = buffer.data(sif0 + 236);
    const auto *sif0_239 = buffer.data(sif0 + 239);
    const auto *sif0_240 = buffer.data(sif0 + 240);
    const auto *sif0_246 = buffer.data(sif0 + 246);
    const auto *sif0_249 = buffer.data(sif0 + 249);
    const auto *sif0_250 = buffer.data(sif0 + 250);

    const auto *sid_48 = buffer.data(sid + 48);
    const auto *sid_51 = buffer.data(sid + 51);
    const auto *sid_54 = buffer.data(sid + 54);
    const auto *sid_57 = buffer.data(sid + 57);
    const auto *sid_59 = buffer.data(sid + 59);
    const auto *sid_60 = buffer.data(sid + 60);
    const auto *sid_63 = buffer.data(sid + 63);
    const auto *sid_65 = buffer.data(sid + 65);
    const auto *sid_66 = buffer.data(sid + 66);
    const auto *sid_69 = buffer.data(sid + 69);
    const auto *sid_71 = buffer.data(sid + 71);
    const auto *sid_72 = buffer.data(sid + 72);
    const auto *sid_75 = buffer.data(sid + 75);
    const auto *sid_77 = buffer.data(sid + 77);
    const auto *sid_78 = buffer.data(sid + 78);
    const auto *sid_81 = buffer.data(sid + 81);
    const auto *sid_82 = buffer.data(sid + 82);
    const auto *sid_83 = buffer.data(sid + 83);
    const auto *sid_84 = buffer.data(sid + 84);
    const auto *sid_87 = buffer.data(sid + 87);
    const auto *sid_88 = buffer.data(sid + 88);
    const auto *sid_89 = buffer.data(sid + 89);
    const auto *sid_90 = buffer.data(sid + 90);
    const auto *sid_93 = buffer.data(sid + 93);
    const auto *sid_94 = buffer.data(sid + 94);
    const auto *sid_95 = buffer.data(sid + 95);
    const auto *sid_96 = buffer.data(sid + 96);
    const auto *sid_99 = buffer.data(sid + 99);
    const auto *sid_100 = buffer.data(sid + 100);
    const auto *sid_101 = buffer.data(sid + 101);
    const auto *sid_102 = buffer.data(sid + 102);
    const auto *sid_105 = buffer.data(sid + 105);
    const auto *sid_106 = buffer.data(sid + 106);
    const auto *sid_107 = buffer.data(sid + 107);
    const auto *sid_108 = buffer.data(sid + 108);
    const auto *sid_111 = buffer.data(sid + 111);
    const auto *sid_112 = buffer.data(sid + 112);
    const auto *sid_113 = buffer.data(sid + 113);
    const auto *sid_114 = buffer.data(sid + 114);
    const auto *sid_117 = buffer.data(sid + 117);
    const auto *sid_118 = buffer.data(sid + 118);
    const auto *sid_119 = buffer.data(sid + 119);
    const auto *sid_120 = buffer.data(sid + 120);
    const auto *sid_123 = buffer.data(sid + 123);
    const auto *sid_124 = buffer.data(sid + 124);
    const auto *sid_125 = buffer.data(sid + 125);
    const auto *sid_126 = buffer.data(sid + 126);
    const auto *sid_129 = buffer.data(sid + 129);
    const auto *sid_130 = buffer.data(sid + 130);
    const auto *sid_131 = buffer.data(sid + 131);
    const auto *sid_135 = buffer.data(sid + 135);
    const auto *sid_136 = buffer.data(sid + 136);
    const auto *sid_137 = buffer.data(sid + 137);
    const auto *sid_138 = buffer.data(sid + 138);
    const auto *sid_141 = buffer.data(sid + 141);
    const auto *sid_142 = buffer.data(sid + 142);
    const auto *sid_143 = buffer.data(sid + 143);
    const auto *sid_144 = buffer.data(sid + 144);
    const auto *sid_147 = buffer.data(sid + 147);
    const auto *sid_148 = buffer.data(sid + 148);
    const auto *sid_149 = buffer.data(sid + 149);
    const auto *sid_150 = buffer.data(sid + 150);
    const auto *sid_153 = buffer.data(sid + 153);
    const auto *sid_154 = buffer.data(sid + 154);
    const auto *sid_155 = buffer.data(sid + 155);

    const auto *sif1_90 = buffer.data(sif1 + 90);
    const auto *sif1_99 = buffer.data(sif1 + 99);
    const auto *sif1_100 = buffer.data(sif1 + 100);
    const auto *sif1_106 = buffer.data(sif1 + 106);
    const auto *sif1_140 = buffer.data(sif1 + 140);
    const auto *sif1_149 = buffer.data(sif1 + 149);
    const auto *sif1_150 = buffer.data(sif1 + 150);
    const auto *sif1_210 = buffer.data(sif1 + 210);
    const auto *sif1_216 = buffer.data(sif1 + 216);
    const auto *sif1_219 = buffer.data(sif1 + 219);
    const auto *sif1_226 = buffer.data(sif1 + 226);
    const auto *sif1_229 = buffer.data(sif1 + 229);
    const auto *sif1_230 = buffer.data(sif1 + 230);
    const auto *sif1_236 = buffer.data(sif1 + 236);
    const auto *sif1_239 = buffer.data(sif1 + 239);
    const auto *sif1_240 = buffer.data(sif1 + 240);
    const auto *sif1_246 = buffer.data(sif1 + 246);
    const auto *sif1_249 = buffer.data(sif1 + 249);
    const auto *sif1_250 = buffer.data(sif1 + 250);

    const auto *skp0_40 = buffer.data(skp0 + 40);
    const auto *skp0_42 = buffer.data(skp0 + 42);
    const auto *skp0_43 = buffer.data(skp0 + 43);
    const auto *skp0_44 = buffer.data(skp0 + 44);
    const auto *skp0_45 = buffer.data(skp0 + 45);
    const auto *skp0_46 = buffer.data(skp0 + 46);
    const auto *skp0_47 = buffer.data(skp0 + 47);
    const auto *skp0_50 = buffer.data(skp0 + 50);
    const auto *skp0_51 = buffer.data(skp0 + 51);
    const auto *skp0_52 = buffer.data(skp0 + 52);
    const auto *skp0_53 = buffer.data(skp0 + 53);
    const auto *skp0_54 = buffer.data(skp0 + 54);
    const auto *skp0_55 = buffer.data(skp0 + 55);
    const auto *skp0_56 = buffer.data(skp0 + 56);
    const auto *skp0_58 = buffer.data(skp0 + 58);
    const auto *skp0_60 = buffer.data(skp0 + 60);
    const auto *skp0_61 = buffer.data(skp0 + 61);
    const auto *skp0_62 = buffer.data(skp0 + 62);

    const auto *skp1_40 = buffer.data(skp1 + 40);
    const auto *skp1_42 = buffer.data(skp1 + 42);
    const auto *skp1_43 = buffer.data(skp1 + 43);
    const auto *skp1_44 = buffer.data(skp1 + 44);
    const auto *skp1_45 = buffer.data(skp1 + 45);
    const auto *skp1_46 = buffer.data(skp1 + 46);
    const auto *skp1_47 = buffer.data(skp1 + 47);
    const auto *skp1_50 = buffer.data(skp1 + 50);
    const auto *skp1_51 = buffer.data(skp1 + 51);
    const auto *skp1_52 = buffer.data(skp1 + 52);
    const auto *skp1_53 = buffer.data(skp1 + 53);
    const auto *skp1_54 = buffer.data(skp1 + 54);
    const auto *skp1_55 = buffer.data(skp1 + 55);
    const auto *skp1_56 = buffer.data(skp1 + 56);
    const auto *skp1_58 = buffer.data(skp1 + 58);
    const auto *skp1_60 = buffer.data(skp1 + 60);
    const auto *skp1_61 = buffer.data(skp1 + 61);
    const auto *skp1_62 = buffer.data(skp1 + 62);

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
    const auto *skd_132 = buffer.data(skd + 132);
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
    const auto *skd_155 = buffer.data(skd + 155);

#pragma omp simd aligned(t_130, t_131, t_132, t_133, pb_y, pc_x, pc_y, pc_z, sif0_90, sid_48, \
                         sid_54, sid_81, sif1_90, skd_78, skd_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = pb_y[k] * sif0_90[k]
                   - f_4 * pc_y[k] * sif1_90[k];

        t_131[k] = f_5 * sid_54[k]
                   + f_3 * pc_y[k] * skd_78[k];

        t_132[k] = f_10 * sid_48[k]
                   + f_3 * pc_z[k] * skd_78[k];

        t_133[k] = f_10 * sid_81[k]
                   + f_3 * pc_x[k] * skd_81[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, pc_x, pc_y, pc_z, sid_51, sid_57, sid_82, \
                         sid_83, skp0_40, skp1_40, skd_81, skd_82, \
                         skd_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_10 * sid_82[k]
                   + f_3 * pc_x[k] * skd_82[k];

        t_135[k] = f_10 * sid_83[k]
                   + f_3 * pc_x[k] * skd_83[k];

        t_136[k] = f_5 * sid_57[k]
                   + f_1 * skp0_40[k]
                   - f_2 * skp1_40[k]
                   + f_3 * pc_y[k] * skd_81[k];

        t_137[k] = f_10 * sid_51[k]
                   + f_3 * pc_z[k] * skd_81[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, pb_y, pc_x, pc_y, sif0_99, sid_59, \
                         sid_84, sif1_99, skp0_42, skp1_42, skd_83, \
                         skd_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_5 * sid_59[k]
                   + f_3 * pc_y[k] * skd_83[k];

        t_139[k] = pb_y[k] * sif0_99[k]
                   - f_4 * pc_y[k] * sif1_99[k];

        t_140[k] = f_10 * sid_84[k]
                   + f_1 * skp0_42[k]
                   - f_2 * skp1_42[k]
                   + f_3 * pc_x[k] * skd_84[k];

        t_141[k] = f_3 * pc_y[k] * skd_84[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, pc_x, pc_z, sid_54, sid_87, sid_88, \
                         sid_89, skd_84, skd_87, skd_88, skd_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_9 * sid_54[k]
                   + f_3 * pc_z[k] * skd_84[k];

        t_143[k] = f_10 * sid_87[k]
                   + f_3 * pc_x[k] * skd_87[k];

        t_144[k] = f_10 * sid_88[k]
                   + f_3 * pc_x[k] * skd_88[k];

        t_145[k] = f_10 * sid_89[k]
                   + f_3 * pc_x[k] * skd_89[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pc_y, pc_z, sid_57, sid_59, skp0_43, \
                         skp0_44, skp1_43, skp1_44, skd_87, skd_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_1 * skp0_43[k]
                   - f_2 * skp1_43[k]
                   + f_3 * pc_y[k] * skd_87[k];

        t_147[k] = f_9 * sid_57[k]
                   + f_3 * pc_z[k] * skd_87[k];

        t_148[k] = f_3 * pc_y[k] * skd_89[k];

        t_149[k] = f_9 * sid_59[k]
                   + f_1 * skp0_44[k]
                   - f_2 * skp1_44[k]
                   + f_3 * pc_z[k] * skd_89[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, pc_x, pc_y, pc_z, sid_60, sid_90, sid_93, \
                         skp0_45, skp1_45, skd_90, skd_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_8 * sid_90[k]
                   + f_1 * skp0_45[k]
                   - f_2 * skp1_45[k]
                   + f_3 * pc_x[k] * skd_90[k];

        t_151[k] = f_7 * sid_60[k]
                   + f_3 * pc_y[k] * skd_90[k];

        t_152[k] = f_3 * pc_z[k] * skd_90[k];

        t_153[k] = f_8 * sid_93[k]
                   + f_3 * pc_x[k] * skd_93[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, pc_x, pc_y, pc_z, sid_63, sid_94, sid_95, \
                         skp0_46, skp1_46, skd_93, skd_94, skd_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = f_8 * sid_94[k]
                   + f_3 * pc_x[k] * skd_94[k];

        t_155[k] = f_8 * sid_95[k]
                   + f_3 * pc_x[k] * skd_95[k];

        t_156[k] = f_7 * sid_63[k]
                   + f_1 * skp0_46[k]
                   - f_2 * skp1_46[k]
                   + f_3 * pc_y[k] * skd_93[k];

        t_157[k] = f_3 * pc_z[k] * skd_93[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, pb_z, pc_y, pc_z, sif0_100, sid_65, \
                         sid_66, sif1_100, skp0_47, skp1_47, skd_95, \
                         skd_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_7 * sid_65[k]
                   + f_3 * pc_y[k] * skd_95[k];

        t_159[k] = f_1 * skp0_47[k]
                   - f_2 * skp1_47[k]
                   + f_3 * pc_z[k] * skd_95[k];

        t_160[k] = pb_z[k] * sif0_100[k]
                   - f_4 * pc_z[k] * sif1_100[k];

        t_161[k] = f_9 * sid_66[k]
                   + f_3 * pc_y[k] * skd_96[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, pc_x, pc_z, sid_60, sid_99, sid_100, \
                         sid_101, skd_96, skd_99, skd_100, skd_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_5 * sid_60[k]
                   + f_3 * pc_z[k] * skd_96[k];

        t_163[k] = f_8 * sid_99[k]
                   + f_3 * pc_x[k] * skd_99[k];

        t_164[k] = f_8 * sid_100[k]
                   + f_3 * pc_x[k] * skd_100[k];

        t_165[k] = f_8 * sid_101[k]
                   + f_3 * pc_x[k] * skd_101[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, pb_z, pc_y, pc_z, sif0_106, sid_63, \
                         sid_65, sid_71, sif1_106, skp0_50, skp1_50, skd_99, \
                         skd_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = pb_z[k] * sif0_106[k]
                   - f_4 * pc_z[k] * sif1_106[k];

        t_167[k] = f_5 * sid_63[k]
                   + f_3 * pc_z[k] * skd_99[k];

        t_168[k] = f_9 * sid_71[k]
                   + f_3 * pc_y[k] * skd_101[k];

        t_169[k] = f_5 * sid_65[k]
                   + f_1 * skp0_50[k]
                   - f_2 * skp1_50[k]
                   + f_3 * pc_z[k] * skd_101[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, pc_x, pc_y, pc_z, sid_66, sid_72, \
                         sid_102, sid_105, skp0_51, skp1_51, skd_102, \
                         skd_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_8 * sid_102[k]
                   + f_1 * skp0_51[k]
                   - f_2 * skp1_51[k]
                   + f_3 * pc_x[k] * skd_102[k];

        t_171[k] = f_10 * sid_72[k]
                   + f_3 * pc_y[k] * skd_102[k];

        t_172[k] = f_8 * sid_66[k]
                   + f_3 * pc_z[k] * skd_102[k];

        t_173[k] = f_8 * sid_105[k]
                   + f_3 * pc_x[k] * skd_105[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, pc_x, pc_y, pc_z, sid_69, sid_75, \
                         sid_106, sid_107, skp0_52, skp1_52, skd_105, skd_106, \
                         skd_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_8 * sid_106[k]
                   + f_3 * pc_x[k] * skd_106[k];

        t_175[k] = f_8 * sid_107[k]
                   + f_3 * pc_x[k] * skd_107[k];

        t_176[k] = f_10 * sid_75[k]
                   + f_1 * skp0_52[k]
                   - f_2 * skp1_52[k]
                   + f_3 * pc_y[k] * skd_105[k];

        t_177[k] = f_8 * sid_69[k]
                   + f_3 * pc_z[k] * skd_105[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, pc_x, pc_y, pc_z, sid_71, sid_77, sid_108, \
                         skp0_53, skp0_54, skp1_53, skp1_54, skd_107, \
                         skd_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_10 * sid_77[k]
                   + f_3 * pc_y[k] * skd_107[k];

        t_179[k] = f_8 * sid_71[k]
                   + f_1 * skp0_53[k]
                   - f_2 * skp1_53[k]
                   + f_3 * pc_z[k] * skd_107[k];

        t_180[k] = f_8 * sid_108[k]
                   + f_1 * skp0_54[k]
                   - f_2 * skp1_54[k]
                   + f_3 * pc_x[k] * skd_108[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pc_x, pc_y, pc_z, sid_72, sid_78, \
                         sid_111, sid_112, skd_108, skd_111, skd_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_8 * sid_78[k]
                   + f_3 * pc_y[k] * skd_108[k];

        t_182[k] = f_10 * sid_72[k]
                   + f_3 * pc_z[k] * skd_108[k];

        t_183[k] = f_8 * sid_111[k]
                   + f_3 * pc_x[k] * skd_111[k];

        t_184[k] = f_8 * sid_112[k]
                   + f_3 * pc_x[k] * skd_112[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pc_x, pc_y, pc_z, sid_75, sid_81, sid_83, \
                         sid_113, skp0_55, skp1_55, skd_111, skd_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_8 * sid_113[k]
                   + f_3 * pc_x[k] * skd_113[k];

        t_186[k] = f_8 * sid_81[k]
                   + f_1 * skp0_55[k]
                   - f_2 * skp1_55[k]
                   + f_3 * pc_y[k] * skd_111[k];

        t_187[k] = f_10 * sid_75[k]
                   + f_3 * pc_z[k] * skd_111[k];

        t_188[k] = f_8 * sid_83[k]
                   + f_3 * pc_y[k] * skd_113[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, pb_y, pc_y, pc_z, sif0_140, sid_77, \
                         sid_78, sid_84, sif1_140, skp0_56, skp1_56, skd_113, \
                         skd_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_10 * sid_77[k]
                   + f_1 * skp0_56[k]
                   - f_2 * skp1_56[k]
                   + f_3 * pc_z[k] * skd_113[k];

        t_190[k] = pb_y[k] * sif0_140[k]
                   - f_4 * pc_y[k] * sif1_140[k];

        t_191[k] = f_5 * sid_84[k]
                   + f_3 * pc_y[k] * skd_114[k];

        t_192[k] = f_9 * sid_78[k]
                   + f_3 * pc_z[k] * skd_114[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, pc_x, pc_y, sid_87, sid_117, sid_118, \
                         sid_119, skp0_58, skp1_58, skd_117, skd_118, \
                         skd_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_8 * sid_117[k]
                   + f_3 * pc_x[k] * skd_117[k];

        t_194[k] = f_8 * sid_118[k]
                   + f_3 * pc_x[k] * skd_118[k];

        t_195[k] = f_8 * sid_119[k]
                   + f_3 * pc_x[k] * skd_119[k];

        t_196[k] = f_5 * sid_87[k]
                   + f_1 * skp0_58[k]
                   - f_2 * skp1_58[k]
                   + f_3 * pc_y[k] * skd_117[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, pb_y, pc_y, pc_z, sif0_149, sid_81, sid_89, \
                         sif1_149, skd_117, skd_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_9 * sid_81[k]
                   + f_3 * pc_z[k] * skd_117[k];

        t_198[k] = f_5 * sid_89[k]
                   + f_3 * pc_y[k] * skd_119[k];

        t_199[k] = pb_y[k] * sif0_149[k]
                   - f_4 * pc_y[k] * sif1_149[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, pc_x, pc_y, pc_z, sid_84, sid_120, \
                         sid_123, skp0_60, skp1_60, skd_120, skd_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_8 * sid_120[k]
                   + f_1 * skp0_60[k]
                   - f_2 * skp1_60[k]
                   + f_3 * pc_x[k] * skd_120[k];

        t_201[k] = f_3 * pc_y[k] * skd_120[k];

        t_202[k] = f_7 * sid_84[k]
                   + f_3 * pc_z[k] * skd_120[k];

        t_203[k] = f_8 * sid_123[k]
                   + f_3 * pc_x[k] * skd_123[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, t_208, pc_x, pc_y, pc_z, sid_87, sid_124, \
                         sid_125, skp0_61, skp1_61, skd_123, skd_124, \
                         skd_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_8 * sid_124[k]
                   + f_3 * pc_x[k] * skd_124[k];

        t_205[k] = f_8 * sid_125[k]
                   + f_3 * pc_x[k] * skd_125[k];

        t_206[k] = f_1 * skp0_61[k]
                   - f_2 * skp1_61[k]
                   + f_3 * pc_y[k] * skd_123[k];

        t_207[k] = f_7 * sid_87[k]
                   + f_3 * pc_z[k] * skd_123[k];

        t_208[k] = f_3 * pc_y[k] * skd_125[k];
    }

#pragma omp simd aligned(t_209, t_210, t_211, pb_x, pc_x, pc_y, pc_z, sif0_210, sid_89, \
                         sid_90, sid_126, sif1_210, skp0_62, skp1_62, skd_125, \
                         skd_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_209[k] = f_7 * sid_89[k]
                   + f_1 * skp0_62[k]
                   - f_2 * skp1_62[k]
                   + f_3 * pc_z[k] * skd_125[k];

        t_210[k] = pb_x[k] * sif0_210[k]
                   + f_10 * sid_126[k]
                   - f_4 * pc_x[k] * sif1_210[k];

        t_211[k] = f_6 * sid_90[k]
                   + f_3 * pc_y[k] * skd_126[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, pc_x, pc_z, sid_129, sid_130, sid_131, \
                         skd_126, skd_129, skd_130, skd_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_3 * pc_z[k] * skd_126[k];

        t_213[k] = f_5 * sid_129[k]
                   + f_3 * pc_x[k] * skd_129[k];

        t_214[k] = f_5 * sid_130[k]
                   + f_3 * pc_x[k] * skd_130[k];

        t_215[k] = f_5 * sid_131[k]
                   + f_3 * pc_x[k] * skd_131[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, pb_x, pc_x, pc_y, pc_z, sif0_216, \
                         sif0_219, sid_95, sif1_216, sif1_219, skd_129, \
                         skd_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = pb_x[k] * sif0_216[k]
                   - f_4 * pc_x[k] * sif1_216[k];

        t_217[k] = f_3 * pc_z[k] * skd_129[k];

        t_218[k] = f_6 * sid_95[k]
                   + f_3 * pc_y[k] * skd_131[k];

        t_219[k] = pb_x[k] * sif0_219[k]
                   - f_4 * pc_x[k] * sif1_219[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, pb_z, pc_x, pc_y, pc_z, sif0_150, sid_90, \
                         sid_96, sid_135, sif1_150, skd_132, skd_135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = pb_z[k] * sif0_150[k]
                   - f_4 * pc_z[k] * sif1_150[k];

        t_221[k] = f_7 * sid_96[k]
                   + f_3 * pc_y[k] * skd_132[k];

        t_222[k] = f_5 * sid_90[k]
                   + f_3 * pc_z[k] * skd_132[k];

        t_223[k] = f_5 * sid_135[k]
                   + f_3 * pc_x[k] * skd_135[k];
    }

#pragma omp simd aligned(t_224, t_225, t_226, t_227, pb_x, pc_x, pc_z, sif0_226, sid_93, \
                         sid_136, sid_137, sif1_226, skd_135, skd_136, \
                         skd_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_224[k] = f_5 * sid_136[k]
                   + f_3 * pc_x[k] * skd_136[k];

        t_225[k] = f_5 * sid_137[k]
                   + f_3 * pc_x[k] * skd_137[k];

        t_226[k] = pb_x[k] * sif0_226[k]
                   - f_4 * pc_x[k] * sif1_226[k];

        t_227[k] = f_5 * sid_93[k]
                   + f_3 * pc_z[k] * skd_135[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, t_231, pb_x, pc_x, pc_y, sif0_229, sif0_230, \
                         sid_101, sid_102, sid_138, sif1_229, sif1_230, skd_137, \
                         skd_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_7 * sid_101[k]
                   + f_3 * pc_y[k] * skd_137[k];

        t_229[k] = pb_x[k] * sif0_229[k]
                   - f_4 * pc_x[k] * sif1_229[k];

        t_230[k] = pb_x[k] * sif0_230[k]
                   + f_10 * sid_138[k]
                   - f_4 * pc_x[k] * sif1_230[k];

        t_231[k] = f_9 * sid_102[k]
                   + f_3 * pc_y[k] * skd_138[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, t_235, pc_x, pc_z, sid_96, sid_141, sid_142, \
                         sid_143, skd_138, skd_141, skd_142, skd_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = f_8 * sid_96[k]
                   + f_3 * pc_z[k] * skd_138[k];

        t_233[k] = f_5 * sid_141[k]
                   + f_3 * pc_x[k] * skd_141[k];

        t_234[k] = f_5 * sid_142[k]
                   + f_3 * pc_x[k] * skd_142[k];

        t_235[k] = f_5 * sid_143[k]
                   + f_3 * pc_x[k] * skd_143[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, t_239, pb_x, pc_x, pc_y, pc_z, sif0_236, \
                         sif0_239, sid_99, sid_107, sif1_236, sif1_239, skd_141, \
                         skd_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = pb_x[k] * sif0_236[k]
                   - f_4 * pc_x[k] * sif1_236[k];

        t_237[k] = f_8 * sid_99[k]
                   + f_3 * pc_z[k] * skd_141[k];

        t_238[k] = f_9 * sid_107[k]
                   + f_3 * pc_y[k] * skd_143[k];

        t_239[k] = pb_x[k] * sif0_239[k]
                   - f_4 * pc_x[k] * sif1_239[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, pb_x, pc_x, pc_y, pc_z, sif0_240, \
                         sid_102, sid_108, sid_144, sid_147, sif1_240, skd_144, \
                         skd_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = pb_x[k] * sif0_240[k]
                   + f_10 * sid_144[k]
                   - f_4 * pc_x[k] * sif1_240[k];

        t_241[k] = f_10 * sid_108[k]
                   + f_3 * pc_y[k] * skd_144[k];

        t_242[k] = f_10 * sid_102[k]
                   + f_3 * pc_z[k] * skd_144[k];

        t_243[k] = f_5 * sid_147[k]
                   + f_3 * pc_x[k] * skd_147[k];
    }

#pragma omp simd aligned(t_244, t_245, t_246, t_247, pb_x, pc_x, pc_z, sif0_246, sid_105, \
                         sid_148, sid_149, sif1_246, skd_147, skd_148, \
                         skd_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_244[k] = f_5 * sid_148[k]
                   + f_3 * pc_x[k] * skd_148[k];

        t_245[k] = f_5 * sid_149[k]
                   + f_3 * pc_x[k] * skd_149[k];

        t_246[k] = pb_x[k] * sif0_246[k]
                   - f_4 * pc_x[k] * sif1_246[k];

        t_247[k] = f_10 * sid_105[k]
                   + f_3 * pc_z[k] * skd_147[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, t_251, pb_x, pc_x, pc_y, sif0_249, sif0_250, \
                         sid_113, sid_114, sid_150, sif1_249, sif1_250, skd_149, \
                         skd_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_10 * sid_113[k]
                   + f_3 * pc_y[k] * skd_149[k];

        t_249[k] = pb_x[k] * sif0_249[k]
                   - f_4 * pc_x[k] * sif1_249[k];

        t_250[k] = pb_x[k] * sif0_250[k]
                   + f_10 * sid_150[k]
                   - f_4 * pc_x[k] * sif1_250[k];

        t_251[k] = f_8 * sid_114[k]
                   + f_3 * pc_y[k] * skd_150[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, pc_x, pc_z, sid_108, sid_153, sid_154, \
                         sid_155, skd_150, skd_153, skd_154, skd_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = f_9 * sid_108[k]
                   + f_3 * pc_z[k] * skd_150[k];

        t_253[k] = f_5 * sid_153[k]
                   + f_3 * pc_x[k] * skd_153[k];

        t_254[k] = f_5 * sid_154[k]
                   + f_3 * pc_x[k] * skd_154[k];

        t_255[k] = f_5 * sid_155[k]
                   + f_3 * pc_x[k] * skd_155[k];
    }
}

static auto
compute_prim_skf_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sif0,
                                                          const size_t sid, const size_t sif1,
                                                          const size_t skp0, const size_t skp1,
                                                          const size_t skd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / q;
    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 3.0 / q;
    const auto f_7 = 2.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 2.0 / q;
    const auto f_10 = 1.5 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sif0_200 = buffer.data(sif0 + 200);
    const auto *sif0_210 = buffer.data(sif0 + 210);
    const auto *sif0_216 = buffer.data(sif0 + 216);
    const auto *sif0_256 = buffer.data(sif0 + 256);
    const auto *sif0_259 = buffer.data(sif0 + 259);
    const auto *sif0_266 = buffer.data(sif0 + 266);
    const auto *sif0_269 = buffer.data(sif0 + 269);
    const auto *sif0_270 = buffer.data(sif0 + 270);
    const auto *sif0_276 = buffer.data(sif0 + 276);
    const auto *sif0_279 = buffer.data(sif0 + 279);

    const auto *sid_111 = buffer.data(sid + 111);
    const auto *sid_114 = buffer.data(sid + 114);
    const auto *sid_117 = buffer.data(sid + 117);
    const auto *sid_119 = buffer.data(sid + 119);
    const auto *sid_120 = buffer.data(sid + 120);
    const auto *sid_123 = buffer.data(sid + 123);
    const auto *sid_125 = buffer.data(sid + 125);
    const auto *sid_126 = buffer.data(sid + 126);
    const auto *sid_129 = buffer.data(sid + 129);
    const auto *sid_131 = buffer.data(sid + 131);
    const auto *sid_132 = buffer.data(sid + 132);
    const auto *sid_135 = buffer.data(sid + 135);
    const auto *sid_137 = buffer.data(sid + 137);
    const auto *sid_138 = buffer.data(sid + 138);
    const auto *sid_141 = buffer.data(sid + 141);
    const auto *sid_143 = buffer.data(sid + 143);
    const auto *sid_144 = buffer.data(sid + 144);
    const auto *sid_147 = buffer.data(sid + 147);
    const auto *sid_149 = buffer.data(sid + 149);
    const auto *sid_150 = buffer.data(sid + 150);
    const auto *sid_153 = buffer.data(sid + 153);
    const auto *sid_155 = buffer.data(sid + 155);
    const auto *sid_156 = buffer.data(sid + 156);
    const auto *sid_159 = buffer.data(sid + 159);
    const auto *sid_160 = buffer.data(sid + 160);
    const auto *sid_161 = buffer.data(sid + 161);
    const auto *sid_162 = buffer.data(sid + 162);
    const auto *sid_165 = buffer.data(sid + 165);
    const auto *sid_166 = buffer.data(sid + 166);
    const auto *sid_167 = buffer.data(sid + 167);

    const auto *sif1_200 = buffer.data(sif1 + 200);
    const auto *sif1_210 = buffer.data(sif1 + 210);
    const auto *sif1_216 = buffer.data(sif1 + 216);
    const auto *sif1_256 = buffer.data(sif1 + 256);
    const auto *sif1_259 = buffer.data(sif1 + 259);
    const auto *sif1_266 = buffer.data(sif1 + 266);
    const auto *sif1_269 = buffer.data(sif1 + 269);
    const auto *sif1_270 = buffer.data(sif1 + 270);
    const auto *sif1_276 = buffer.data(sif1 + 276);
    const auto *sif1_279 = buffer.data(sif1 + 279);

    const auto *skp0_84 = buffer.data(skp0 + 84);
    const auto *skp0_85 = buffer.data(skp0 + 85);
    const auto *skp0_86 = buffer.data(skp0 + 86);
    const auto *skp0_89 = buffer.data(skp0 + 89);
    const auto *skp0_90 = buffer.data(skp0 + 90);
    const auto *skp0_91 = buffer.data(skp0 + 91);
    const auto *skp0_92 = buffer.data(skp0 + 92);
    const auto *skp0_93 = buffer.data(skp0 + 93);
    const auto *skp0_94 = buffer.data(skp0 + 94);
    const auto *skp0_95 = buffer.data(skp0 + 95);
    const auto *skp0_96 = buffer.data(skp0 + 96);
    const auto *skp0_97 = buffer.data(skp0 + 97);
    const auto *skp0_98 = buffer.data(skp0 + 98);
    const auto *skp0_99 = buffer.data(skp0 + 99);
    const auto *skp0_100 = buffer.data(skp0 + 100);
    const auto *skp0_101 = buffer.data(skp0 + 101);
    const auto *skp0_105 = buffer.data(skp0 + 105);
    const auto *skp0_106 = buffer.data(skp0 + 106);
    const auto *skp0_107 = buffer.data(skp0 + 107);

    const auto *skp1_84 = buffer.data(skp1 + 84);
    const auto *skp1_85 = buffer.data(skp1 + 85);
    const auto *skp1_86 = buffer.data(skp1 + 86);
    const auto *skp1_89 = buffer.data(skp1 + 89);
    const auto *skp1_90 = buffer.data(skp1 + 90);
    const auto *skp1_91 = buffer.data(skp1 + 91);
    const auto *skp1_92 = buffer.data(skp1 + 92);
    const auto *skp1_93 = buffer.data(skp1 + 93);
    const auto *skp1_94 = buffer.data(skp1 + 94);
    const auto *skp1_95 = buffer.data(skp1 + 95);
    const auto *skp1_96 = buffer.data(skp1 + 96);
    const auto *skp1_97 = buffer.data(skp1 + 97);
    const auto *skp1_98 = buffer.data(skp1 + 98);
    const auto *skp1_99 = buffer.data(skp1 + 99);
    const auto *skp1_100 = buffer.data(skp1 + 100);
    const auto *skp1_101 = buffer.data(skp1 + 101);
    const auto *skp1_105 = buffer.data(skp1 + 105);
    const auto *skp1_106 = buffer.data(skp1 + 106);
    const auto *skp1_107 = buffer.data(skp1 + 107);

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
    const auto *skd_204 = buffer.data(skd + 204);
    const auto *skd_207 = buffer.data(skd + 207);
    const auto *skd_208 = buffer.data(skd + 208);
    const auto *skd_209 = buffer.data(skd + 209);
    const auto *skd_210 = buffer.data(skd + 210);
    const auto *skd_213 = buffer.data(skd + 213);
    const auto *skd_214 = buffer.data(skd + 214);
    const auto *skd_215 = buffer.data(skd + 215);

#pragma omp simd aligned(t_256, t_257, t_258, t_259, pb_x, pc_x, pc_y, pc_z, sif0_256, \
                         sif0_259, sid_111, sid_119, sif1_256, sif1_259, skd_153, \
                         skd_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = pb_x[k] * sif0_256[k]
                   - f_4 * pc_x[k] * sif1_256[k];

        t_257[k] = f_9 * sid_111[k]
                   + f_3 * pc_z[k] * skd_153[k];

        t_258[k] = f_8 * sid_119[k]
                   + f_3 * pc_y[k] * skd_155[k];

        t_259[k] = pb_x[k] * sif0_259[k]
                   - f_4 * pc_x[k] * sif1_259[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pb_y, pc_x, pc_y, pc_z, sif0_200, \
                         sid_114, sid_120, sid_159, sif1_200, skd_156, \
                         skd_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = pb_y[k] * sif0_200[k]
                   - f_4 * pc_y[k] * sif1_200[k];

        t_261[k] = f_5 * sid_120[k]
                   + f_3 * pc_y[k] * skd_156[k];

        t_262[k] = f_7 * sid_114[k]
                   + f_3 * pc_z[k] * skd_156[k];

        t_263[k] = f_5 * sid_159[k]
                   + f_3 * pc_x[k] * skd_159[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, pb_x, pc_x, pc_z, sif0_266, sid_117, \
                         sid_160, sid_161, sif1_266, skd_159, skd_160, \
                         skd_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_5 * sid_160[k]
                   + f_3 * pc_x[k] * skd_160[k];

        t_265[k] = f_5 * sid_161[k]
                   + f_3 * pc_x[k] * skd_161[k];

        t_266[k] = pb_x[k] * sif0_266[k]
                   - f_4 * pc_x[k] * sif1_266[k];

        t_267[k] = f_7 * sid_117[k]
                   + f_3 * pc_z[k] * skd_159[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, t_271, pb_x, pc_x, pc_y, sif0_269, sif0_270, \
                         sid_125, sid_162, sif1_269, sif1_270, skd_161, \
                         skd_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_5 * sid_125[k]
                   + f_3 * pc_y[k] * skd_161[k];

        t_269[k] = pb_x[k] * sif0_269[k]
                   - f_4 * pc_x[k] * sif1_269[k];

        t_270[k] = pb_x[k] * sif0_270[k]
                   + f_10 * sid_162[k]
                   - f_4 * pc_x[k] * sif1_270[k];

        t_271[k] = f_3 * pc_y[k] * skd_162[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, pc_x, pc_z, sid_120, sid_165, sid_166, \
                         sid_167, skd_162, skd_165, skd_166, skd_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_6 * sid_120[k]
                   + f_3 * pc_z[k] * skd_162[k];

        t_273[k] = f_5 * sid_165[k]
                   + f_3 * pc_x[k] * skd_165[k];

        t_274[k] = f_5 * sid_166[k]
                   + f_3 * pc_x[k] * skd_166[k];

        t_275[k] = f_5 * sid_167[k]
                   + f_3 * pc_x[k] * skd_167[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pb_x, pc_x, pc_y, pc_z, sif0_276, \
                         sif0_279, sid_123, sif1_276, sif1_279, skd_165, \
                         skd_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = pb_x[k] * sif0_276[k]
                   - f_4 * pc_x[k] * sif1_276[k];

        t_277[k] = f_6 * sid_123[k]
                   + f_3 * pc_z[k] * skd_165[k];

        t_278[k] = f_3 * pc_y[k] * skd_167[k];

        t_279[k] = pb_x[k] * sif0_279[k]
                   - f_4 * pc_x[k] * sif1_279[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, t_285, pc_x, pc_y, pc_z, sid_126, \
                         skp0_84, skp1_84, skd_168, skd_171, skd_172, \
                         skd_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_1 * skp0_84[k]
                   - f_2 * skp1_84[k]
                   + f_3 * pc_x[k] * skd_168[k];

        t_281[k] = f_0 * sid_126[k]
                   + f_3 * pc_y[k] * skd_168[k];

        t_282[k] = f_3 * pc_z[k] * skd_168[k];

        t_283[k] = f_3 * pc_x[k] * skd_171[k];

        t_284[k] = f_3 * pc_x[k] * skd_172[k];

        t_285[k] = f_3 * pc_x[k] * skd_173[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, pc_y, pc_z, sid_129, sid_131, skp0_85, \
                         skp0_86, skp1_85, skp1_86, skd_171, skd_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = f_0 * sid_129[k]
                   + f_1 * skp0_85[k]
                   - f_2 * skp1_85[k]
                   + f_3 * pc_y[k] * skd_171[k];

        t_287[k] = f_3 * pc_z[k] * skd_171[k];

        t_288[k] = f_0 * sid_131[k]
                   + f_3 * pc_y[k] * skd_173[k];

        t_289[k] = f_1 * skp0_86[k]
                   - f_2 * skp1_86[k]
                   + f_3 * pc_z[k] * skd_173[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, pb_z, pc_x, pc_y, pc_z, sif0_210, \
                         sid_126, sid_132, sif1_210, skd_174, skd_177, \
                         skd_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = pb_z[k] * sif0_210[k]
                   - f_4 * pc_z[k] * sif1_210[k];

        t_291[k] = f_6 * sid_132[k]
                   + f_3 * pc_y[k] * skd_174[k];

        t_292[k] = f_5 * sid_126[k]
                   + f_3 * pc_z[k] * skd_174[k];

        t_293[k] = f_3 * pc_x[k] * skd_177[k];

        t_294[k] = f_3 * pc_x[k] * skd_178[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, pb_z, pc_x, pc_y, pc_z, sif0_216, \
                         sid_129, sid_137, sif1_216, skd_177, skd_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = f_3 * pc_x[k] * skd_179[k];

        t_296[k] = pb_z[k] * sif0_216[k]
                   - f_4 * pc_z[k] * sif1_216[k];

        t_297[k] = f_5 * sid_129[k]
                   + f_3 * pc_z[k] * skd_177[k];

        t_298[k] = f_6 * sid_137[k]
                   + f_3 * pc_y[k] * skd_179[k];
    }

#pragma omp simd aligned(t_299, t_300, t_301, t_302, pc_x, pc_y, pc_z, sid_131, sid_132, \
                         sid_138, skp0_89, skp0_90, skp1_89, skp1_90, skd_179, \
                         skd_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_299[k] = f_5 * sid_131[k]
                   + f_1 * skp0_89[k]
                   - f_2 * skp1_89[k]
                   + f_3 * pc_z[k] * skd_179[k];

        t_300[k] = f_1 * skp0_90[k]
                   - f_2 * skp1_90[k]
                   + f_3 * pc_x[k] * skd_180[k];

        t_301[k] = f_7 * sid_138[k]
                   + f_3 * pc_y[k] * skd_180[k];

        t_302[k] = f_8 * sid_132[k]
                   + f_3 * pc_z[k] * skd_180[k];
    }

#pragma omp simd aligned(t_303, t_304, t_305, t_306, t_307, pc_x, pc_y, pc_z, sid_135, \
                         sid_141, skp0_91, skp1_91, skd_183, skd_184, \
                         skd_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_303[k] = f_3 * pc_x[k] * skd_183[k];

        t_304[k] = f_3 * pc_x[k] * skd_184[k];

        t_305[k] = f_3 * pc_x[k] * skd_185[k];

        t_306[k] = f_7 * sid_141[k]
                   + f_1 * skp0_91[k]
                   - f_2 * skp1_91[k]
                   + f_3 * pc_y[k] * skd_183[k];

        t_307[k] = f_8 * sid_135[k]
                   + f_3 * pc_z[k] * skd_183[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, t_311, pc_x, pc_y, pc_z, sid_137, sid_143, \
                         sid_144, skp0_92, skp0_93, skp1_92, skp1_93, skd_185, \
                         skd_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = f_7 * sid_143[k]
                   + f_3 * pc_y[k] * skd_185[k];

        t_309[k] = f_8 * sid_137[k]
                   + f_1 * skp0_92[k]
                   - f_2 * skp1_92[k]
                   + f_3 * pc_z[k] * skd_185[k];

        t_310[k] = f_1 * skp0_93[k]
                   - f_2 * skp1_93[k]
                   + f_3 * pc_x[k] * skd_186[k];

        t_311[k] = f_9 * sid_144[k]
                   + f_3 * pc_y[k] * skd_186[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, t_315, t_316, pc_x, pc_y, pc_z, sid_138, \
                         sid_147, skp0_94, skp1_94, skd_186, skd_189, skd_190, \
                         skd_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = f_10 * sid_138[k]
                   + f_3 * pc_z[k] * skd_186[k];

        t_313[k] = f_3 * pc_x[k] * skd_189[k];

        t_314[k] = f_3 * pc_x[k] * skd_190[k];

        t_315[k] = f_3 * pc_x[k] * skd_191[k];

        t_316[k] = f_9 * sid_147[k]
                   + f_1 * skp0_94[k]
                   - f_2 * skp1_94[k]
                   + f_3 * pc_y[k] * skd_189[k];
    }

#pragma omp simd aligned(t_317, t_318, t_319, pc_y, pc_z, sid_141, sid_143, sid_149, skp0_95, \
                         skp1_95, skd_189, skd_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_317[k] = f_10 * sid_141[k]
                   + f_3 * pc_z[k] * skd_189[k];

        t_318[k] = f_9 * sid_149[k]
                   + f_3 * pc_y[k] * skd_191[k];

        t_319[k] = f_10 * sid_143[k]
                   + f_1 * skp0_95[k]
                   - f_2 * skp1_95[k]
                   + f_3 * pc_z[k] * skd_191[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, pc_x, pc_y, pc_z, sid_144, \
                         sid_150, skp0_96, skp1_96, skd_192, skd_195, \
                         skd_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = f_1 * skp0_96[k]
                   - f_2 * skp1_96[k]
                   + f_3 * pc_x[k] * skd_192[k];

        t_321[k] = f_10 * sid_150[k]
                   + f_3 * pc_y[k] * skd_192[k];

        t_322[k] = f_9 * sid_144[k]
                   + f_3 * pc_z[k] * skd_192[k];

        t_323[k] = f_3 * pc_x[k] * skd_195[k];

        t_324[k] = f_3 * pc_x[k] * skd_196[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, pc_x, pc_y, pc_z, sid_147, sid_153, \
                         sid_155, skp0_97, skp1_97, skd_195, skd_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = f_3 * pc_x[k] * skd_197[k];

        t_326[k] = f_10 * sid_153[k]
                   + f_1 * skp0_97[k]
                   - f_2 * skp1_97[k]
                   + f_3 * pc_y[k] * skd_195[k];

        t_327[k] = f_9 * sid_147[k]
                   + f_3 * pc_z[k] * skd_195[k];

        t_328[k] = f_10 * sid_155[k]
                   + f_3 * pc_y[k] * skd_197[k];
    }

#pragma omp simd aligned(t_329, t_330, t_331, t_332, pc_x, pc_y, pc_z, sid_149, sid_150, \
                         sid_156, skp0_98, skp0_99, skp1_98, skp1_99, skd_197, \
                         skd_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_329[k] = f_9 * sid_149[k]
                   + f_1 * skp0_98[k]
                   - f_2 * skp1_98[k]
                   + f_3 * pc_z[k] * skd_197[k];

        t_330[k] = f_1 * skp0_99[k]
                   - f_2 * skp1_99[k]
                   + f_3 * pc_x[k] * skd_198[k];

        t_331[k] = f_8 * sid_156[k]
                   + f_3 * pc_y[k] * skd_198[k];

        t_332[k] = f_7 * sid_150[k]
                   + f_3 * pc_z[k] * skd_198[k];
    }

#pragma omp simd aligned(t_333, t_334, t_335, t_336, t_337, pc_x, pc_y, pc_z, sid_153, \
                         sid_159, skp0_100, skp1_100, skd_201, skd_202, \
                         skd_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_333[k] = f_3 * pc_x[k] * skd_201[k];

        t_334[k] = f_3 * pc_x[k] * skd_202[k];

        t_335[k] = f_3 * pc_x[k] * skd_203[k];

        t_336[k] = f_8 * sid_159[k]
                   + f_1 * skp0_100[k]
                   - f_2 * skp1_100[k]
                   + f_3 * pc_y[k] * skd_201[k];

        t_337[k] = f_7 * sid_153[k]
                   + f_3 * pc_z[k] * skd_201[k];
    }

#pragma omp simd aligned(t_338, t_339, t_340, t_341, pb_y, pc_y, pc_z, sif0_270, sid_155, \
                         sid_161, sid_162, sif1_270, skp0_101, skp1_101, skd_203, \
                         skd_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_338[k] = f_8 * sid_161[k]
                   + f_3 * pc_y[k] * skd_203[k];

        t_339[k] = f_7 * sid_155[k]
                   + f_1 * skp0_101[k]
                   - f_2 * skp1_101[k]
                   + f_3 * pc_z[k] * skd_203[k];

        t_340[k] = pb_y[k] * sif0_270[k]
                   - f_4 * pc_y[k] * sif1_270[k];

        t_341[k] = f_5 * sid_162[k]
                   + f_3 * pc_y[k] * skd_204[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, t_345, pc_x, pc_z, sid_156, skd_204, skd_207, \
                         skd_208, skd_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = f_6 * sid_156[k]
                   + f_3 * pc_z[k] * skd_204[k];

        t_343[k] = f_3 * pc_x[k] * skd_207[k];

        t_344[k] = f_3 * pc_x[k] * skd_208[k];

        t_345[k] = f_3 * pc_x[k] * skd_209[k];
    }

#pragma omp simd aligned(t_346, t_347, t_348, t_349, pb_y, pc_y, pc_z, sif0_276, sif0_279, \
                         sid_159, sid_165, sid_167, sif1_276, sif1_279, skd_207, \
                         skd_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_346[k] = pb_y[k] * sif0_276[k]
                   + f_10 * sid_165[k]
                   - f_4 * pc_y[k] * sif1_276[k];

        t_347[k] = f_6 * sid_159[k]
                   + f_3 * pc_z[k] * skd_207[k];

        t_348[k] = f_5 * sid_167[k]
                   + f_3 * pc_y[k] * skd_209[k];

        t_349[k] = pb_y[k] * sif0_279[k]
                   - f_4 * pc_y[k] * sif1_279[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, t_355, pc_x, pc_y, pc_z, sid_162, \
                         skp0_105, skp1_105, skd_210, skd_213, skd_214, \
                         skd_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = f_1 * skp0_105[k]
                   - f_2 * skp1_105[k]
                   + f_3 * pc_x[k] * skd_210[k];

        t_351[k] = f_3 * pc_y[k] * skd_210[k];

        t_352[k] = f_0 * sid_162[k]
                   + f_3 * pc_z[k] * skd_210[k];

        t_353[k] = f_3 * pc_x[k] * skd_213[k];

        t_354[k] = f_3 * pc_x[k] * skd_214[k];

        t_355[k] = f_3 * pc_x[k] * skd_215[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, t_359, pc_y, pc_z, sid_165, sid_167, skp0_106, \
                         skp0_107, skp1_106, skp1_107, skd_213, \
                         skd_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_1 * skp0_106[k]
                   - f_2 * skp1_106[k]
                   + f_3 * pc_y[k] * skd_213[k];

        t_357[k] = f_0 * sid_165[k]
                   + f_3 * pc_z[k] * skd_213[k];

        t_358[k] = f_3 * pc_y[k] * skd_215[k];

        t_359[k] = f_0 * sid_167[k]
                   + f_1 * skp0_107[k]
                   - f_2 * skp1_107[k]
                   + f_3 * pc_z[k] * skd_215[k];
    }
}

auto
compute_prim_skf_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t sif0, const size_t sid,
                                                   const size_t sif1, const size_t skp0,
                                                   const size_t skp1, const size_t skd,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_skf_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, sif0, sid,
                                                              sif1, skp0, skp1, skd, ncols,
                                                              gamma, p, q);

    compute_prim_skf_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, sif0, sid,
                                                              sif1, skp0, skp1, skd, ncols,
                                                              gamma, p, q);

    compute_prim_skf_three_center_electron_repulsion_0_piece2(buffer, target, pb, pc, sif0, sid,
                                                              sif1, skp0, skp1, skd, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
