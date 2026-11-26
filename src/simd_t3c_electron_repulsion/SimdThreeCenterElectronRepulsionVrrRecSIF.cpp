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


#include "SimdThreeCenterElectronRepulsionVrrRecSIF.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_sif_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t shf0,
                                                          const size_t shd, const size_t shf1,
                                                          const size_t sip0, const size_t sip1,
                                                          const size_t sid, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 2.5 / q;
    const auto f_7 = 2.0 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 1.5 / q;

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

    const auto *shf0_0 = buffer.data(shf0 + 0);
    const auto *shf0_6 = buffer.data(shf0 + 6);
    const auto *shf0_9 = buffer.data(shf0 + 9);
    const auto *shf0_16 = buffer.data(shf0 + 16);
    const auto *shf0_20 = buffer.data(shf0 + 20);
    const auto *shf0_29 = buffer.data(shf0 + 29);
    const auto *shf0_30 = buffer.data(shf0 + 30);
    const auto *shf0_36 = buffer.data(shf0 + 36);
    const auto *shf0_50 = buffer.data(shf0 + 50);
    const auto *shf0_59 = buffer.data(shf0 + 59);
    const auto *shf0_60 = buffer.data(shf0 + 60);
    const auto *shf0_66 = buffer.data(shf0 + 66);

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
    const auto *shd_69 = buffer.data(shd + 69);
    const auto *shd_70 = buffer.data(shd + 70);
    const auto *shd_71 = buffer.data(shd + 71);
    const auto *shd_72 = buffer.data(shd + 72);
    const auto *shd_75 = buffer.data(shd + 75);
    const auto *shd_76 = buffer.data(shd + 76);
    const auto *shd_77 = buffer.data(shd + 77);

    const auto *shf1_0 = buffer.data(shf1 + 0);
    const auto *shf1_6 = buffer.data(shf1 + 6);
    const auto *shf1_9 = buffer.data(shf1 + 9);
    const auto *shf1_16 = buffer.data(shf1 + 16);
    const auto *shf1_20 = buffer.data(shf1 + 20);
    const auto *shf1_29 = buffer.data(shf1 + 29);
    const auto *shf1_30 = buffer.data(shf1 + 30);
    const auto *shf1_36 = buffer.data(shf1 + 36);
    const auto *shf1_50 = buffer.data(shf1 + 50);
    const auto *shf1_59 = buffer.data(shf1 + 59);
    const auto *shf1_60 = buffer.data(shf1 + 60);
    const auto *shf1_66 = buffer.data(shf1 + 66);

    const auto *sip0_0 = buffer.data(sip0 + 0);
    const auto *sip0_1 = buffer.data(sip0 + 1);
    const auto *sip0_2 = buffer.data(sip0 + 2);
    const auto *sip0_4 = buffer.data(sip0 + 4);
    const auto *sip0_8 = buffer.data(sip0 + 8);
    const auto *sip0_9 = buffer.data(sip0 + 9);
    const auto *sip0_10 = buffer.data(sip0 + 10);
    const auto *sip0_11 = buffer.data(sip0 + 11);
    const auto *sip0_15 = buffer.data(sip0 + 15);
    const auto *sip0_16 = buffer.data(sip0 + 16);
    const auto *sip0_17 = buffer.data(sip0 + 17);
    const auto *sip0_18 = buffer.data(sip0 + 18);
    const auto *sip0_19 = buffer.data(sip0 + 19);
    const auto *sip0_20 = buffer.data(sip0 + 20);
    const auto *sip0_23 = buffer.data(sip0 + 23);
    const auto *sip0_25 = buffer.data(sip0 + 25);
    const auto *sip0_27 = buffer.data(sip0 + 27);
    const auto *sip0_28 = buffer.data(sip0 + 28);
    const auto *sip0_29 = buffer.data(sip0 + 29);
    const auto *sip0_30 = buffer.data(sip0 + 30);
    const auto *sip0_31 = buffer.data(sip0 + 31);
    const auto *sip0_32 = buffer.data(sip0 + 32);
    const auto *sip0_35 = buffer.data(sip0 + 35);
    const auto *sip0_36 = buffer.data(sip0 + 36);
    const auto *sip0_37 = buffer.data(sip0 + 37);
    const auto *sip0_38 = buffer.data(sip0 + 38);

    const auto *sip1_0 = buffer.data(sip1 + 0);
    const auto *sip1_1 = buffer.data(sip1 + 1);
    const auto *sip1_2 = buffer.data(sip1 + 2);
    const auto *sip1_4 = buffer.data(sip1 + 4);
    const auto *sip1_8 = buffer.data(sip1 + 8);
    const auto *sip1_9 = buffer.data(sip1 + 9);
    const auto *sip1_10 = buffer.data(sip1 + 10);
    const auto *sip1_11 = buffer.data(sip1 + 11);
    const auto *sip1_15 = buffer.data(sip1 + 15);
    const auto *sip1_16 = buffer.data(sip1 + 16);
    const auto *sip1_17 = buffer.data(sip1 + 17);
    const auto *sip1_18 = buffer.data(sip1 + 18);
    const auto *sip1_19 = buffer.data(sip1 + 19);
    const auto *sip1_20 = buffer.data(sip1 + 20);
    const auto *sip1_23 = buffer.data(sip1 + 23);
    const auto *sip1_25 = buffer.data(sip1 + 25);
    const auto *sip1_27 = buffer.data(sip1 + 27);
    const auto *sip1_28 = buffer.data(sip1 + 28);
    const auto *sip1_29 = buffer.data(sip1 + 29);
    const auto *sip1_30 = buffer.data(sip1 + 30);
    const auto *sip1_31 = buffer.data(sip1 + 31);
    const auto *sip1_32 = buffer.data(sip1 + 32);
    const auto *sip1_35 = buffer.data(sip1 + 35);
    const auto *sip1_36 = buffer.data(sip1 + 36);
    const auto *sip1_37 = buffer.data(sip1 + 37);
    const auto *sip1_38 = buffer.data(sip1 + 38);

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
    const auto *sid_66 = buffer.data(sid + 66);
    const auto *sid_69 = buffer.data(sid + 69);
    const auto *sid_70 = buffer.data(sid + 70);
    const auto *sid_71 = buffer.data(sid + 71);
    const auto *sid_72 = buffer.data(sid + 72);
    const auto *sid_75 = buffer.data(sid + 75);
    const auto *sid_76 = buffer.data(sid + 76);
    const auto *sid_77 = buffer.data(sid + 77);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pc_x, pc_y, pc_z, shd_0, shd_3, shd_4, \
                         sip0_0, sip1_0, sid_0, sid_3, sid_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * shd_0[k]
                 + f_1 * sip0_0[k]
                 - f_2 * sip1_0[k]
                 + f_3 * pc_x[k] * sid_0[k];

        t_1[k] = f_3 * pc_y[k] * sid_0[k];

        t_2[k] = f_3 * pc_z[k] * sid_0[k];

        t_3[k] = f_0 * shd_3[k]
                 + f_3 * pc_x[k] * sid_3[k];

        t_4[k] = f_0 * shd_4[k]
                 + f_3 * pc_x[k] * sid_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pc_x, pc_y, pc_z, shd_5, sip0_1, sip0_2, \
                         sip1_1, sip1_2, sid_3, sid_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * shd_5[k]
                 + f_3 * pc_x[k] * sid_5[k];

        t_6[k] = f_1 * sip0_1[k]
                 - f_2 * sip1_1[k]
                 + f_3 * pc_y[k] * sid_3[k];

        t_7[k] = f_3 * pc_z[k] * sid_3[k];

        t_8[k] = f_3 * pc_y[k] * sid_5[k];

        t_9[k] = f_1 * sip0_2[k]
                 - f_2 * sip1_2[k]
                 + f_3 * pc_z[k] * sid_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pb_y, pc_x, pc_y, pc_z, shf0_0, shd_0, shd_9, \
                         shf1_0, sid_6, sid_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pb_y[k] * shf0_0[k]
                  - f_4 * pc_y[k] * shf1_0[k];

        t_11[k] = f_5 * shd_0[k]
                  + f_3 * pc_y[k] * sid_6[k];

        t_12[k] = f_3 * pc_z[k] * sid_6[k];

        t_13[k] = f_6 * shd_9[k]
                  + f_3 * pc_x[k] * sid_9[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pc_x, pc_y, pc_z, shd_3, shd_10, shd_11, \
                         sip0_4, sip1_4, sid_9, sid_10, sid_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_6 * shd_10[k]
                  + f_3 * pc_x[k] * sid_10[k];

        t_15[k] = f_6 * shd_11[k]
                  + f_3 * pc_x[k] * sid_11[k];

        t_16[k] = f_5 * shd_3[k]
                  + f_1 * sip0_4[k]
                  - f_2 * sip1_4[k]
                  + f_3 * pc_y[k] * sid_9[k];

        t_17[k] = f_3 * pc_z[k] * sid_9[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pb_y, pb_z, pc_y, pc_z, shf0_0, shf0_9, \
                         shd_5, shf1_0, shf1_9, sid_11, sid_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * shd_5[k]
                  + f_3 * pc_y[k] * sid_11[k];

        t_19[k] = pb_y[k] * shf0_9[k]
                  - f_4 * pc_y[k] * shf1_9[k];

        t_20[k] = pb_z[k] * shf0_0[k]
                  - f_4 * pc_z[k] * shf1_0[k];

        t_21[k] = f_3 * pc_y[k] * sid_12[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pc_x, pc_z, shd_0, shd_15, shd_16, shd_17, \
                         sid_12, sid_15, sid_16, sid_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_5 * shd_0[k]
                  + f_3 * pc_z[k] * sid_12[k];

        t_23[k] = f_6 * shd_15[k]
                  + f_3 * pc_x[k] * sid_15[k];

        t_24[k] = f_6 * shd_16[k]
                  + f_3 * pc_x[k] * sid_16[k];

        t_25[k] = f_6 * shd_17[k]
                  + f_3 * pc_x[k] * sid_17[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pb_z, pc_y, pc_z, shf0_6, shd_3, shd_5, \
                         shf1_6, sip0_8, sip1_8, sid_15, sid_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pb_z[k] * shf0_6[k]
                  - f_4 * pc_z[k] * shf1_6[k];

        t_27[k] = f_5 * shd_3[k]
                  + f_3 * pc_z[k] * sid_15[k];

        t_28[k] = f_3 * pc_y[k] * sid_17[k];

        t_29[k] = f_5 * shd_5[k]
                  + f_1 * sip0_8[k]
                  - f_2 * sip1_8[k]
                  + f_3 * pc_z[k] * sid_17[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pc_x, pc_y, pc_z, shd_6, shd_18, shd_21, \
                         sip0_9, sip1_9, sid_18, sid_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_7 * shd_18[k]
                  + f_1 * sip0_9[k]
                  - f_2 * sip1_9[k]
                  + f_3 * pc_x[k] * sid_18[k];

        t_31[k] = f_8 * shd_6[k]
                  + f_3 * pc_y[k] * sid_18[k];

        t_32[k] = f_3 * pc_z[k] * sid_18[k];

        t_33[k] = f_7 * shd_21[k]
                  + f_3 * pc_x[k] * sid_21[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pc_x, pc_y, pc_z, shd_9, shd_22, shd_23, \
                         sip0_10, sip1_10, sid_21, sid_22, sid_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_7 * shd_22[k]
                  + f_3 * pc_x[k] * sid_22[k];

        t_35[k] = f_7 * shd_23[k]
                  + f_3 * pc_x[k] * sid_23[k];

        t_36[k] = f_8 * shd_9[k]
                  + f_1 * sip0_10[k]
                  - f_2 * sip1_10[k]
                  + f_3 * pc_y[k] * sid_21[k];

        t_37[k] = f_3 * pc_z[k] * sid_21[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pb_y, pc_y, pc_z, shf0_20, shd_11, shd_12, \
                         shf1_20, sip0_11, sip1_11, sid_23, sid_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_8 * shd_11[k]
                  + f_3 * pc_y[k] * sid_23[k];

        t_39[k] = f_1 * sip0_11[k]
                  - f_2 * sip1_11[k]
                  + f_3 * pc_z[k] * sid_23[k];

        t_40[k] = pb_y[k] * shf0_20[k]
                  - f_4 * pc_y[k] * shf1_20[k];

        t_41[k] = f_5 * shd_12[k]
                  + f_3 * pc_y[k] * sid_24[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pc_x, pc_z, shd_6, shd_27, shd_28, shd_29, \
                         sid_24, sid_27, sid_28, sid_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_5 * shd_6[k]
                  + f_3 * pc_z[k] * sid_24[k];

        t_43[k] = f_7 * shd_27[k]
                  + f_3 * pc_x[k] * sid_27[k];

        t_44[k] = f_7 * shd_28[k]
                  + f_3 * pc_x[k] * sid_28[k];

        t_45[k] = f_7 * shd_29[k]
                  + f_3 * pc_x[k] * sid_29[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pb_y, pb_z, pc_y, pc_z, shf0_16, shf0_29, \
                         shd_9, shd_17, shf1_16, shf1_29, sid_27, \
                         sid_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = pb_z[k] * shf0_16[k]
                  - f_4 * pc_z[k] * shf1_16[k];

        t_47[k] = f_5 * shd_9[k]
                  + f_3 * pc_z[k] * sid_27[k];

        t_48[k] = f_5 * shd_17[k]
                  + f_3 * pc_y[k] * sid_29[k];

        t_49[k] = pb_y[k] * shf0_29[k]
                  - f_4 * pc_y[k] * shf1_29[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pc_x, pc_y, pc_z, shd_12, shd_30, shd_33, \
                         sip0_15, sip1_15, sid_30, sid_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_7 * shd_30[k]
                  + f_1 * sip0_15[k]
                  - f_2 * sip1_15[k]
                  + f_3 * pc_x[k] * sid_30[k];

        t_51[k] = f_3 * pc_y[k] * sid_30[k];

        t_52[k] = f_8 * shd_12[k]
                  + f_3 * pc_z[k] * sid_30[k];

        t_53[k] = f_7 * shd_33[k]
                  + f_3 * pc_x[k] * sid_33[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, pc_x, pc_y, pc_z, shd_15, shd_34, \
                         shd_35, sip0_16, sip1_16, sid_33, sid_34, \
                         sid_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_7 * shd_34[k]
                  + f_3 * pc_x[k] * sid_34[k];

        t_55[k] = f_7 * shd_35[k]
                  + f_3 * pc_x[k] * sid_35[k];

        t_56[k] = f_1 * sip0_16[k]
                  - f_2 * sip1_16[k]
                  + f_3 * pc_y[k] * sid_33[k];

        t_57[k] = f_8 * shd_15[k]
                  + f_3 * pc_z[k] * sid_33[k];

        t_58[k] = f_3 * pc_y[k] * sid_35[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pc_x, pc_y, pc_z, shd_17, shd_18, shd_36, \
                         sip0_17, sip0_18, sip1_17, sip1_18, sid_35, \
                         sid_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_8 * shd_17[k]
                  + f_1 * sip0_17[k]
                  - f_2 * sip1_17[k]
                  + f_3 * pc_z[k] * sid_35[k];

        t_60[k] = f_9 * shd_36[k]
                  + f_1 * sip0_18[k]
                  - f_2 * sip1_18[k]
                  + f_3 * pc_x[k] * sid_36[k];

        t_61[k] = f_9 * shd_18[k]
                  + f_3 * pc_y[k] * sid_36[k];

        t_62[k] = f_3 * pc_z[k] * sid_36[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pc_x, pc_y, shd_21, shd_39, shd_40, shd_41, \
                         sip0_19, sip1_19, sid_39, sid_40, sid_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_9 * shd_39[k]
                  + f_3 * pc_x[k] * sid_39[k];

        t_64[k] = f_9 * shd_40[k]
                  + f_3 * pc_x[k] * sid_40[k];

        t_65[k] = f_9 * shd_41[k]
                  + f_3 * pc_x[k] * sid_41[k];

        t_66[k] = f_9 * shd_21[k]
                  + f_1 * sip0_19[k]
                  - f_2 * sip1_19[k]
                  + f_3 * pc_y[k] * sid_39[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, pb_z, pc_y, pc_z, shf0_30, shd_23, shf1_30, \
                         sip0_20, sip1_20, sid_39, sid_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_3 * pc_z[k] * sid_39[k];

        t_68[k] = f_9 * shd_23[k]
                  + f_3 * pc_y[k] * sid_41[k];

        t_69[k] = f_1 * sip0_20[k]
                  - f_2 * sip1_20[k]
                  + f_3 * pc_z[k] * sid_41[k];

        t_70[k] = pb_z[k] * shf0_30[k]
                  - f_4 * pc_z[k] * shf1_30[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pc_x, pc_y, pc_z, shd_18, shd_24, shd_45, \
                         shd_46, sid_42, sid_45, sid_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_8 * shd_24[k]
                  + f_3 * pc_y[k] * sid_42[k];

        t_72[k] = f_5 * shd_18[k]
                  + f_3 * pc_z[k] * sid_42[k];

        t_73[k] = f_9 * shd_45[k]
                  + f_3 * pc_x[k] * sid_45[k];

        t_74[k] = f_9 * shd_46[k]
                  + f_3 * pc_x[k] * sid_46[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pb_z, pc_x, pc_y, pc_z, shf0_36, shd_21, \
                         shd_29, shd_47, shf1_36, sid_45, sid_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_9 * shd_47[k]
                  + f_3 * pc_x[k] * sid_47[k];

        t_76[k] = pb_z[k] * shf0_36[k]
                  - f_4 * pc_z[k] * shf1_36[k];

        t_77[k] = f_5 * shd_21[k]
                  + f_3 * pc_z[k] * sid_45[k];

        t_78[k] = f_8 * shd_29[k]
                  + f_3 * pc_y[k] * sid_47[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pb_y, pc_y, pc_z, shf0_50, shd_23, shd_24, \
                         shd_30, shf1_50, sip0_23, sip1_23, sid_47, \
                         sid_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_5 * shd_23[k]
                  + f_1 * sip0_23[k]
                  - f_2 * sip1_23[k]
                  + f_3 * pc_z[k] * sid_47[k];

        t_80[k] = pb_y[k] * shf0_50[k]
                  - f_4 * pc_y[k] * shf1_50[k];

        t_81[k] = f_5 * shd_30[k]
                  + f_3 * pc_y[k] * sid_48[k];

        t_82[k] = f_8 * shd_24[k]
                  + f_3 * pc_z[k] * sid_48[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pc_x, pc_y, shd_33, shd_51, shd_52, shd_53, \
                         sip0_25, sip1_25, sid_51, sid_52, sid_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_9 * shd_51[k]
                  + f_3 * pc_x[k] * sid_51[k];

        t_84[k] = f_9 * shd_52[k]
                  + f_3 * pc_x[k] * sid_52[k];

        t_85[k] = f_9 * shd_53[k]
                  + f_3 * pc_x[k] * sid_53[k];

        t_86[k] = f_5 * shd_33[k]
                  + f_1 * sip0_25[k]
                  - f_2 * sip1_25[k]
                  + f_3 * pc_y[k] * sid_51[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, pb_y, pc_y, pc_z, shf0_59, shd_27, shd_35, shf1_59, \
                         sid_51, sid_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_8 * shd_27[k]
                  + f_3 * pc_z[k] * sid_51[k];

        t_88[k] = f_5 * shd_35[k]
                  + f_3 * pc_y[k] * sid_53[k];

        t_89[k] = pb_y[k] * shf0_59[k]
                  - f_4 * pc_y[k] * shf1_59[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, pc_x, pc_y, pc_z, shd_30, shd_54, shd_57, \
                         sip0_27, sip1_27, sid_54, sid_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_9 * shd_54[k]
                  + f_1 * sip0_27[k]
                  - f_2 * sip1_27[k]
                  + f_3 * pc_x[k] * sid_54[k];

        t_91[k] = f_3 * pc_y[k] * sid_54[k];

        t_92[k] = f_9 * shd_30[k]
                  + f_3 * pc_z[k] * sid_54[k];

        t_93[k] = f_9 * shd_57[k]
                  + f_3 * pc_x[k] * sid_57[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, pc_x, pc_y, pc_z, shd_33, shd_58, \
                         shd_59, sip0_28, sip1_28, sid_57, sid_58, \
                         sid_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_9 * shd_58[k]
                  + f_3 * pc_x[k] * sid_58[k];

        t_95[k] = f_9 * shd_59[k]
                  + f_3 * pc_x[k] * sid_59[k];

        t_96[k] = f_1 * sip0_28[k]
                  - f_2 * sip1_28[k]
                  + f_3 * pc_y[k] * sid_57[k];

        t_97[k] = f_9 * shd_33[k]
                  + f_3 * pc_z[k] * sid_57[k];

        t_98[k] = f_3 * pc_y[k] * sid_59[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, pc_x, pc_y, pc_z, shd_35, shd_36, shd_60, \
                         sip0_29, sip0_30, sip1_29, sip1_30, sid_59, \
                         sid_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_9 * shd_35[k]
                  + f_1 * sip0_29[k]
                  - f_2 * sip1_29[k]
                  + f_3 * pc_z[k] * sid_59[k];

        t_100[k] = f_8 * shd_60[k]
                   + f_1 * sip0_30[k]
                   - f_2 * sip1_30[k]
                   + f_3 * pc_x[k] * sid_60[k];

        t_101[k] = f_7 * shd_36[k]
                   + f_3 * pc_y[k] * sid_60[k];

        t_102[k] = f_3 * pc_z[k] * sid_60[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, pc_x, pc_y, shd_39, shd_63, shd_64, \
                         shd_65, sip0_31, sip1_31, sid_63, sid_64, \
                         sid_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_8 * shd_63[k]
                   + f_3 * pc_x[k] * sid_63[k];

        t_104[k] = f_8 * shd_64[k]
                   + f_3 * pc_x[k] * sid_64[k];

        t_105[k] = f_8 * shd_65[k]
                   + f_3 * pc_x[k] * sid_65[k];

        t_106[k] = f_7 * shd_39[k]
                   + f_1 * sip0_31[k]
                   - f_2 * sip1_31[k]
                   + f_3 * pc_y[k] * sid_63[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, pb_z, pc_y, pc_z, shf0_60, shd_41, \
                         shf1_60, sip0_32, sip1_32, sid_63, sid_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_3 * pc_z[k] * sid_63[k];

        t_108[k] = f_7 * shd_41[k]
                   + f_3 * pc_y[k] * sid_65[k];

        t_109[k] = f_1 * sip0_32[k]
                   - f_2 * sip1_32[k]
                   + f_3 * pc_z[k] * sid_65[k];

        t_110[k] = pb_z[k] * shf0_60[k]
                   - f_4 * pc_z[k] * shf1_60[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, pc_x, pc_y, pc_z, shd_36, shd_42, shd_69, \
                         shd_70, sid_66, sid_69, sid_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_9 * shd_42[k]
                   + f_3 * pc_y[k] * sid_66[k];

        t_112[k] = f_5 * shd_36[k]
                   + f_3 * pc_z[k] * sid_66[k];

        t_113[k] = f_8 * shd_69[k]
                   + f_3 * pc_x[k] * sid_69[k];

        t_114[k] = f_8 * shd_70[k]
                   + f_3 * pc_x[k] * sid_70[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, pb_z, pc_x, pc_y, pc_z, shf0_66, shd_39, \
                         shd_47, shd_71, shf1_66, sid_69, sid_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_8 * shd_71[k]
                   + f_3 * pc_x[k] * sid_71[k];

        t_116[k] = pb_z[k] * shf0_66[k]
                   - f_4 * pc_z[k] * shf1_66[k];

        t_117[k] = f_5 * shd_39[k]
                   + f_3 * pc_z[k] * sid_69[k];

        t_118[k] = f_9 * shd_47[k]
                   + f_3 * pc_y[k] * sid_71[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, pc_x, pc_y, pc_z, shd_41, shd_48, shd_72, \
                         sip0_35, sip0_36, sip1_35, sip1_36, sid_71, \
                         sid_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = f_5 * shd_41[k]
                   + f_1 * sip0_35[k]
                   - f_2 * sip1_35[k]
                   + f_3 * pc_z[k] * sid_71[k];

        t_120[k] = f_8 * shd_72[k]
                   + f_1 * sip0_36[k]
                   - f_2 * sip1_36[k]
                   + f_3 * pc_x[k] * sid_72[k];

        t_121[k] = f_8 * shd_48[k]
                   + f_3 * pc_y[k] * sid_72[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, pc_x, pc_z, shd_42, shd_75, shd_76, \
                         shd_77, sid_72, sid_75, sid_76, sid_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_8 * shd_42[k]
                   + f_3 * pc_z[k] * sid_72[k];

        t_123[k] = f_8 * shd_75[k]
                   + f_3 * pc_x[k] * sid_75[k];

        t_124[k] = f_8 * shd_76[k]
                   + f_3 * pc_x[k] * sid_76[k];

        t_125[k] = f_8 * shd_77[k]
                   + f_3 * pc_x[k] * sid_77[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, pc_y, pc_z, shd_45, shd_47, shd_51, \
                         shd_53, sip0_37, sip0_38, sip1_37, sip1_38, sid_75, \
                         sid_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_8 * shd_51[k]
                   + f_1 * sip0_37[k]
                   - f_2 * sip1_37[k]
                   + f_3 * pc_y[k] * sid_75[k];

        t_127[k] = f_8 * shd_45[k]
                   + f_3 * pc_z[k] * sid_75[k];

        t_128[k] = f_8 * shd_53[k]
                   + f_3 * pc_y[k] * sid_77[k];

        t_129[k] = f_8 * shd_47[k]
                   + f_1 * sip0_38[k]
                   - f_2 * sip1_38[k]
                   + f_3 * pc_z[k] * sid_77[k];
    }
}

static auto
compute_prim_sif_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t shf0,
                                                          const size_t shd, const size_t shf1,
                                                          const size_t sip0, const size_t sip1,
                                                          const size_t sid, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 2.5 / q;
    const auto f_7 = 2.0 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 1.5 / q;

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
    auto *t_256 = buffer.data(target + 256);
    auto *t_257 = buffer.data(target + 257);
    auto *t_258 = buffer.data(target + 258);
    auto *t_259 = buffer.data(target + 259);
    auto *t_260 = buffer.data(target + 260);
    auto *t_261 = buffer.data(target + 261);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *shf0_90 = buffer.data(shf0 + 90);
    const auto *shf0_99 = buffer.data(shf0 + 99);
    const auto *shf0_100 = buffer.data(shf0 + 100);
    const auto *shf0_140 = buffer.data(shf0 + 140);
    const auto *shf0_150 = buffer.data(shf0 + 150);
    const auto *shf0_156 = buffer.data(shf0 + 156);
    const auto *shf0_159 = buffer.data(shf0 + 159);
    const auto *shf0_166 = buffer.data(shf0 + 166);
    const auto *shf0_169 = buffer.data(shf0 + 169);
    const auto *shf0_170 = buffer.data(shf0 + 170);
    const auto *shf0_176 = buffer.data(shf0 + 176);
    const auto *shf0_179 = buffer.data(shf0 + 179);
    const auto *shf0_180 = buffer.data(shf0 + 180);
    const auto *shf0_186 = buffer.data(shf0 + 186);
    const auto *shf0_189 = buffer.data(shf0 + 189);
    const auto *shf0_196 = buffer.data(shf0 + 196);
    const auto *shf0_199 = buffer.data(shf0 + 199);
    const auto *shf0_200 = buffer.data(shf0 + 200);
    const auto *shf0_206 = buffer.data(shf0 + 206);
    const auto *shf0_209 = buffer.data(shf0 + 209);

    const auto *shd_48 = buffer.data(shd + 48);
    const auto *shd_51 = buffer.data(shd + 51);
    const auto *shd_54 = buffer.data(shd + 54);
    const auto *shd_57 = buffer.data(shd + 57);
    const auto *shd_59 = buffer.data(shd + 59);
    const auto *shd_60 = buffer.data(shd + 60);
    const auto *shd_63 = buffer.data(shd + 63);
    const auto *shd_65 = buffer.data(shd + 65);
    const auto *shd_66 = buffer.data(shd + 66);
    const auto *shd_69 = buffer.data(shd + 69);
    const auto *shd_71 = buffer.data(shd + 71);
    const auto *shd_72 = buffer.data(shd + 72);
    const auto *shd_75 = buffer.data(shd + 75);
    const auto *shd_77 = buffer.data(shd + 77);
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

    const auto *shf1_90 = buffer.data(shf1 + 90);
    const auto *shf1_99 = buffer.data(shf1 + 99);
    const auto *shf1_100 = buffer.data(shf1 + 100);
    const auto *shf1_140 = buffer.data(shf1 + 140);
    const auto *shf1_150 = buffer.data(shf1 + 150);
    const auto *shf1_156 = buffer.data(shf1 + 156);
    const auto *shf1_159 = buffer.data(shf1 + 159);
    const auto *shf1_166 = buffer.data(shf1 + 166);
    const auto *shf1_169 = buffer.data(shf1 + 169);
    const auto *shf1_170 = buffer.data(shf1 + 170);
    const auto *shf1_176 = buffer.data(shf1 + 176);
    const auto *shf1_179 = buffer.data(shf1 + 179);
    const auto *shf1_180 = buffer.data(shf1 + 180);
    const auto *shf1_186 = buffer.data(shf1 + 186);
    const auto *shf1_189 = buffer.data(shf1 + 189);
    const auto *shf1_196 = buffer.data(shf1 + 196);
    const auto *shf1_199 = buffer.data(shf1 + 199);
    const auto *shf1_200 = buffer.data(shf1 + 200);
    const auto *shf1_206 = buffer.data(shf1 + 206);
    const auto *shf1_209 = buffer.data(shf1 + 209);

    const auto *sip0_40 = buffer.data(sip0 + 40);
    const auto *sip0_42 = buffer.data(sip0 + 42);
    const auto *sip0_43 = buffer.data(sip0 + 43);
    const auto *sip0_44 = buffer.data(sip0 + 44);
    const auto *sip0_63 = buffer.data(sip0 + 63);
    const auto *sip0_64 = buffer.data(sip0 + 64);
    const auto *sip0_65 = buffer.data(sip0 + 65);
    const auto *sip0_68 = buffer.data(sip0 + 68);
    const auto *sip0_69 = buffer.data(sip0 + 69);
    const auto *sip0_70 = buffer.data(sip0 + 70);
    const auto *sip0_71 = buffer.data(sip0 + 71);
    const auto *sip0_72 = buffer.data(sip0 + 72);
    const auto *sip0_73 = buffer.data(sip0 + 73);
    const auto *sip0_74 = buffer.data(sip0 + 74);
    const auto *sip0_75 = buffer.data(sip0 + 75);
    const auto *sip0_76 = buffer.data(sip0 + 76);
    const auto *sip0_77 = buffer.data(sip0 + 77);

    const auto *sip1_40 = buffer.data(sip1 + 40);
    const auto *sip1_42 = buffer.data(sip1 + 42);
    const auto *sip1_43 = buffer.data(sip1 + 43);
    const auto *sip1_44 = buffer.data(sip1 + 44);
    const auto *sip1_63 = buffer.data(sip1 + 63);
    const auto *sip1_64 = buffer.data(sip1 + 64);
    const auto *sip1_65 = buffer.data(sip1 + 65);
    const auto *sip1_68 = buffer.data(sip1 + 68);
    const auto *sip1_69 = buffer.data(sip1 + 69);
    const auto *sip1_70 = buffer.data(sip1 + 70);
    const auto *sip1_71 = buffer.data(sip1 + 71);
    const auto *sip1_72 = buffer.data(sip1 + 72);
    const auto *sip1_73 = buffer.data(sip1 + 73);
    const auto *sip1_74 = buffer.data(sip1 + 74);
    const auto *sip1_75 = buffer.data(sip1 + 75);
    const auto *sip1_76 = buffer.data(sip1 + 76);
    const auto *sip1_77 = buffer.data(sip1 + 77);

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
    const auto *sid_132 = buffer.data(sid + 132);
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
    const auto *sid_156 = buffer.data(sid + 156);

#pragma omp simd aligned(t_130, t_131, t_132, t_133, pb_y, pc_x, pc_y, pc_z, shf0_90, shd_48, \
                         shd_54, shd_81, shf1_90, sid_78, sid_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = pb_y[k] * shf0_90[k]
                   - f_4 * pc_y[k] * shf1_90[k];

        t_131[k] = f_5 * shd_54[k]
                   + f_3 * pc_y[k] * sid_78[k];

        t_132[k] = f_9 * shd_48[k]
                   + f_3 * pc_z[k] * sid_78[k];

        t_133[k] = f_8 * shd_81[k]
                   + f_3 * pc_x[k] * sid_81[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, pc_x, pc_y, pc_z, shd_51, shd_57, shd_82, \
                         shd_83, sip0_40, sip1_40, sid_81, sid_82, \
                         sid_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_8 * shd_82[k]
                   + f_3 * pc_x[k] * sid_82[k];

        t_135[k] = f_8 * shd_83[k]
                   + f_3 * pc_x[k] * sid_83[k];

        t_136[k] = f_5 * shd_57[k]
                   + f_1 * sip0_40[k]
                   - f_2 * sip1_40[k]
                   + f_3 * pc_y[k] * sid_81[k];

        t_137[k] = f_9 * shd_51[k]
                   + f_3 * pc_z[k] * sid_81[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, pb_y, pc_x, pc_y, shf0_99, shd_59, \
                         shd_84, shf1_99, sip0_42, sip1_42, sid_83, \
                         sid_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_5 * shd_59[k]
                   + f_3 * pc_y[k] * sid_83[k];

        t_139[k] = pb_y[k] * shf0_99[k]
                   - f_4 * pc_y[k] * shf1_99[k];

        t_140[k] = f_8 * shd_84[k]
                   + f_1 * sip0_42[k]
                   - f_2 * sip1_42[k]
                   + f_3 * pc_x[k] * sid_84[k];

        t_141[k] = f_3 * pc_y[k] * sid_84[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, pc_x, pc_z, shd_54, shd_87, shd_88, \
                         shd_89, sid_84, sid_87, sid_88, sid_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_7 * shd_54[k]
                   + f_3 * pc_z[k] * sid_84[k];

        t_143[k] = f_8 * shd_87[k]
                   + f_3 * pc_x[k] * sid_87[k];

        t_144[k] = f_8 * shd_88[k]
                   + f_3 * pc_x[k] * sid_88[k];

        t_145[k] = f_8 * shd_89[k]
                   + f_3 * pc_x[k] * sid_89[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pc_y, pc_z, shd_57, shd_59, sip0_43, \
                         sip0_44, sip1_43, sip1_44, sid_87, sid_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_1 * sip0_43[k]
                   - f_2 * sip1_43[k]
                   + f_3 * pc_y[k] * sid_87[k];

        t_147[k] = f_7 * shd_57[k]
                   + f_3 * pc_z[k] * sid_87[k];

        t_148[k] = f_3 * pc_y[k] * sid_89[k];

        t_149[k] = f_7 * shd_59[k]
                   + f_1 * sip0_44[k]
                   - f_2 * sip1_44[k]
                   + f_3 * pc_z[k] * sid_89[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, pb_x, pc_x, pc_y, pc_z, shf0_150, shd_60, \
                         shd_90, shd_93, shf1_150, sid_90, sid_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = pb_x[k] * shf0_150[k]
                   + f_9 * shd_90[k]
                   - f_4 * pc_x[k] * shf1_150[k];

        t_151[k] = f_6 * shd_60[k]
                   + f_3 * pc_y[k] * sid_90[k];

        t_152[k] = f_3 * pc_z[k] * sid_90[k];

        t_153[k] = f_5 * shd_93[k]
                   + f_3 * pc_x[k] * sid_93[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, pb_x, pc_x, pc_z, shf0_156, shd_94, \
                         shd_95, shf1_156, sid_93, sid_94, sid_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = f_5 * shd_94[k]
                   + f_3 * pc_x[k] * sid_94[k];

        t_155[k] = f_5 * shd_95[k]
                   + f_3 * pc_x[k] * sid_95[k];

        t_156[k] = pb_x[k] * shf0_156[k]
                   - f_4 * pc_x[k] * shf1_156[k];

        t_157[k] = f_3 * pc_z[k] * sid_93[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, pb_x, pb_z, pc_x, pc_y, pc_z, shf0_100, \
                         shf0_159, shd_65, shf1_100, shf1_159, sid_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_6 * shd_65[k]
                   + f_3 * pc_y[k] * sid_95[k];

        t_159[k] = pb_x[k] * shf0_159[k]
                   - f_4 * pc_x[k] * shf1_159[k];

        t_160[k] = pb_z[k] * shf0_100[k]
                   - f_4 * pc_z[k] * shf1_100[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, t_164, pc_x, pc_y, pc_z, shd_60, shd_66, shd_99, \
                         shd_100, sid_96, sid_99, sid_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = f_7 * shd_66[k]
                   + f_3 * pc_y[k] * sid_96[k];

        t_162[k] = f_5 * shd_60[k]
                   + f_3 * pc_z[k] * sid_96[k];

        t_163[k] = f_5 * shd_99[k]
                   + f_3 * pc_x[k] * sid_99[k];

        t_164[k] = f_5 * shd_100[k]
                   + f_3 * pc_x[k] * sid_100[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, pb_x, pc_x, pc_y, pc_z, shf0_166, shd_63, \
                         shd_71, shd_101, shf1_166, sid_99, sid_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = f_5 * shd_101[k]
                   + f_3 * pc_x[k] * sid_101[k];

        t_166[k] = pb_x[k] * shf0_166[k]
                   - f_4 * pc_x[k] * shf1_166[k];

        t_167[k] = f_5 * shd_63[k]
                   + f_3 * pc_z[k] * sid_99[k];

        t_168[k] = f_7 * shd_71[k]
                   + f_3 * pc_y[k] * sid_101[k];
    }

#pragma omp simd aligned(t_169, t_170, t_171, t_172, pb_x, pc_x, pc_y, pc_z, shf0_169, \
                         shf0_170, shd_66, shd_72, shd_102, shf1_169, shf1_170, \
                         sid_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = pb_x[k] * shf0_169[k]
                   - f_4 * pc_x[k] * shf1_169[k];

        t_170[k] = pb_x[k] * shf0_170[k]
                   + f_9 * shd_102[k]
                   - f_4 * pc_x[k] * shf1_170[k];

        t_171[k] = f_9 * shd_72[k]
                   + f_3 * pc_y[k] * sid_102[k];

        t_172[k] = f_8 * shd_66[k]
                   + f_3 * pc_z[k] * sid_102[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, pb_x, pc_x, shf0_176, shd_105, shd_106, \
                         shd_107, shf1_176, sid_105, sid_106, sid_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = f_5 * shd_105[k]
                   + f_3 * pc_x[k] * sid_105[k];

        t_174[k] = f_5 * shd_106[k]
                   + f_3 * pc_x[k] * sid_106[k];

        t_175[k] = f_5 * shd_107[k]
                   + f_3 * pc_x[k] * sid_107[k];

        t_176[k] = pb_x[k] * shf0_176[k]
                   - f_4 * pc_x[k] * shf1_176[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, pb_x, pc_x, pc_y, pc_z, shf0_179, shd_69, \
                         shd_77, shf1_179, sid_105, sid_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = f_8 * shd_69[k]
                   + f_3 * pc_z[k] * sid_105[k];

        t_178[k] = f_9 * shd_77[k]
                   + f_3 * pc_y[k] * sid_107[k];

        t_179[k] = pb_x[k] * shf0_179[k]
                   - f_4 * pc_x[k] * shf1_179[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pb_x, pc_x, pc_y, pc_z, shf0_180, shd_72, \
                         shd_78, shd_108, shd_111, shf1_180, sid_108, \
                         sid_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = pb_x[k] * shf0_180[k]
                   + f_9 * shd_108[k]
                   - f_4 * pc_x[k] * shf1_180[k];

        t_181[k] = f_8 * shd_78[k]
                   + f_3 * pc_y[k] * sid_108[k];

        t_182[k] = f_9 * shd_72[k]
                   + f_3 * pc_z[k] * sid_108[k];

        t_183[k] = f_5 * shd_111[k]
                   + f_3 * pc_x[k] * sid_111[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, pb_x, pc_x, pc_z, shf0_186, shd_75, \
                         shd_112, shd_113, shf1_186, sid_111, sid_112, \
                         sid_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_5 * shd_112[k]
                   + f_3 * pc_x[k] * sid_112[k];

        t_185[k] = f_5 * shd_113[k]
                   + f_3 * pc_x[k] * sid_113[k];

        t_186[k] = pb_x[k] * shf0_186[k]
                   - f_4 * pc_x[k] * shf1_186[k];

        t_187[k] = f_9 * shd_75[k]
                   + f_3 * pc_z[k] * sid_111[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, pb_x, pb_y, pc_x, pc_y, shf0_140, \
                         shf0_189, shd_83, shd_84, shf1_140, shf1_189, sid_113, \
                         sid_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = f_8 * shd_83[k]
                   + f_3 * pc_y[k] * sid_113[k];

        t_189[k] = pb_x[k] * shf0_189[k]
                   - f_4 * pc_x[k] * shf1_189[k];

        t_190[k] = pb_y[k] * shf0_140[k]
                   - f_4 * pc_y[k] * shf1_140[k];

        t_191[k] = f_5 * shd_84[k]
                   + f_3 * pc_y[k] * sid_114[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, pc_x, pc_z, shd_78, shd_117, shd_118, \
                         shd_119, sid_114, sid_117, sid_118, sid_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_7 * shd_78[k]
                   + f_3 * pc_z[k] * sid_114[k];

        t_193[k] = f_5 * shd_117[k]
                   + f_3 * pc_x[k] * sid_117[k];

        t_194[k] = f_5 * shd_118[k]
                   + f_3 * pc_x[k] * sid_118[k];

        t_195[k] = f_5 * shd_119[k]
                   + f_3 * pc_x[k] * sid_119[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, pb_x, pc_x, pc_y, pc_z, shf0_196, \
                         shf0_199, shd_81, shd_89, shf1_196, shf1_199, sid_117, \
                         sid_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = pb_x[k] * shf0_196[k]
                   - f_4 * pc_x[k] * shf1_196[k];

        t_197[k] = f_7 * shd_81[k]
                   + f_3 * pc_z[k] * sid_117[k];

        t_198[k] = f_5 * shd_89[k]
                   + f_3 * pc_y[k] * sid_119[k];

        t_199[k] = pb_x[k] * shf0_199[k]
                   - f_4 * pc_x[k] * shf1_199[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, pb_x, pc_x, pc_y, pc_z, shf0_200, shd_84, \
                         shd_120, shd_123, shf1_200, sid_120, sid_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = pb_x[k] * shf0_200[k]
                   + f_9 * shd_120[k]
                   - f_4 * pc_x[k] * shf1_200[k];

        t_201[k] = f_3 * pc_y[k] * sid_120[k];

        t_202[k] = f_6 * shd_84[k]
                   + f_3 * pc_z[k] * sid_120[k];

        t_203[k] = f_5 * shd_123[k]
                   + f_3 * pc_x[k] * sid_123[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, pb_x, pc_x, pc_z, shf0_206, shd_87, \
                         shd_124, shd_125, shf1_206, sid_123, sid_124, \
                         sid_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_5 * shd_124[k]
                   + f_3 * pc_x[k] * sid_124[k];

        t_205[k] = f_5 * shd_125[k]
                   + f_3 * pc_x[k] * sid_125[k];

        t_206[k] = pb_x[k] * shf0_206[k]
                   - f_4 * pc_x[k] * shf1_206[k];

        t_207[k] = f_6 * shd_87[k]
                   + f_3 * pc_z[k] * sid_123[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, t_212, pb_x, pc_x, pc_y, pc_z, shf0_209, \
                         shd_90, shf1_209, sip0_63, sip1_63, sid_125, \
                         sid_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_3 * pc_y[k] * sid_125[k];

        t_209[k] = pb_x[k] * shf0_209[k]
                   - f_4 * pc_x[k] * shf1_209[k];

        t_210[k] = f_1 * sip0_63[k]
                   - f_2 * sip1_63[k]
                   + f_3 * pc_x[k] * sid_126[k];

        t_211[k] = f_0 * shd_90[k]
                   + f_3 * pc_y[k] * sid_126[k];

        t_212[k] = f_3 * pc_z[k] * sid_126[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, t_217, t_218, pc_x, pc_y, pc_z, shd_93, \
                         shd_95, sip0_64, sip1_64, sid_129, sid_130, \
                         sid_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_3 * pc_x[k] * sid_129[k];

        t_214[k] = f_3 * pc_x[k] * sid_130[k];

        t_215[k] = f_3 * pc_x[k] * sid_131[k];

        t_216[k] = f_0 * shd_93[k]
                   + f_1 * sip0_64[k]
                   - f_2 * sip1_64[k]
                   + f_3 * pc_y[k] * sid_129[k];

        t_217[k] = f_3 * pc_z[k] * sid_129[k];

        t_218[k] = f_0 * shd_95[k]
                   + f_3 * pc_y[k] * sid_131[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, t_222, pb_z, pc_y, pc_z, shf0_150, shd_90, \
                         shd_96, shf1_150, sip0_65, sip1_65, sid_131, \
                         sid_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = f_1 * sip0_65[k]
                   - f_2 * sip1_65[k]
                   + f_3 * pc_z[k] * sid_131[k];

        t_220[k] = pb_z[k] * shf0_150[k]
                   - f_4 * pc_z[k] * shf1_150[k];

        t_221[k] = f_6 * shd_96[k]
                   + f_3 * pc_y[k] * sid_132[k];

        t_222[k] = f_5 * shd_90[k]
                   + f_3 * pc_z[k] * sid_132[k];
    }

#pragma omp simd aligned(t_223, t_224, t_225, t_226, t_227, pb_z, pc_x, pc_z, shf0_156, \
                         shd_93, shf1_156, sid_135, sid_136, sid_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = f_3 * pc_x[k] * sid_135[k];

        t_224[k] = f_3 * pc_x[k] * sid_136[k];

        t_225[k] = f_3 * pc_x[k] * sid_137[k];

        t_226[k] = pb_z[k] * shf0_156[k]
                   - f_4 * pc_z[k] * shf1_156[k];

        t_227[k] = f_5 * shd_93[k]
                   + f_3 * pc_z[k] * sid_135[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, t_231, pc_x, pc_y, pc_z, shd_95, shd_101, \
                         shd_102, sip0_68, sip0_69, sip1_68, sip1_69, sid_137, \
                         sid_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_6 * shd_101[k]
                   + f_3 * pc_y[k] * sid_137[k];

        t_229[k] = f_5 * shd_95[k]
                   + f_1 * sip0_68[k]
                   - f_2 * sip1_68[k]
                   + f_3 * pc_z[k] * sid_137[k];

        t_230[k] = f_1 * sip0_69[k]
                   - f_2 * sip1_69[k]
                   + f_3 * pc_x[k] * sid_138[k];

        t_231[k] = f_7 * shd_102[k]
                   + f_3 * pc_y[k] * sid_138[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, t_235, t_236, pc_x, pc_y, pc_z, shd_96, shd_105, \
                         sip0_70, sip1_70, sid_138, sid_141, sid_142, \
                         sid_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = f_8 * shd_96[k]
                   + f_3 * pc_z[k] * sid_138[k];

        t_233[k] = f_3 * pc_x[k] * sid_141[k];

        t_234[k] = f_3 * pc_x[k] * sid_142[k];

        t_235[k] = f_3 * pc_x[k] * sid_143[k];

        t_236[k] = f_7 * shd_105[k]
                   + f_1 * sip0_70[k]
                   - f_2 * sip1_70[k]
                   + f_3 * pc_y[k] * sid_141[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, pc_y, pc_z, shd_99, shd_101, shd_107, sip0_71, \
                         sip1_71, sid_141, sid_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = f_8 * shd_99[k]
                   + f_3 * pc_z[k] * sid_141[k];

        t_238[k] = f_7 * shd_107[k]
                   + f_3 * pc_y[k] * sid_143[k];

        t_239[k] = f_8 * shd_101[k]
                   + f_1 * sip0_71[k]
                   - f_2 * sip1_71[k]
                   + f_3 * pc_z[k] * sid_143[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, pc_x, pc_y, pc_z, shd_102, \
                         shd_108, sip0_72, sip1_72, sid_144, sid_147, \
                         sid_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_1 * sip0_72[k]
                   - f_2 * sip1_72[k]
                   + f_3 * pc_x[k] * sid_144[k];

        t_241[k] = f_9 * shd_108[k]
                   + f_3 * pc_y[k] * sid_144[k];

        t_242[k] = f_9 * shd_102[k]
                   + f_3 * pc_z[k] * sid_144[k];

        t_243[k] = f_3 * pc_x[k] * sid_147[k];

        t_244[k] = f_3 * pc_x[k] * sid_148[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, pc_x, pc_y, pc_z, shd_105, shd_111, \
                         shd_113, sip0_73, sip1_73, sid_147, sid_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_3 * pc_x[k] * sid_149[k];

        t_246[k] = f_9 * shd_111[k]
                   + f_1 * sip0_73[k]
                   - f_2 * sip1_73[k]
                   + f_3 * pc_y[k] * sid_147[k];

        t_247[k] = f_9 * shd_105[k]
                   + f_3 * pc_z[k] * sid_147[k];

        t_248[k] = f_9 * shd_113[k]
                   + f_3 * pc_y[k] * sid_149[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, t_252, pc_x, pc_y, pc_z, shd_107, shd_108, \
                         shd_114, sip0_74, sip0_75, sip1_74, sip1_75, sid_149, \
                         sid_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = f_9 * shd_107[k]
                   + f_1 * sip0_74[k]
                   - f_2 * sip1_74[k]
                   + f_3 * pc_z[k] * sid_149[k];

        t_250[k] = f_1 * sip0_75[k]
                   - f_2 * sip1_75[k]
                   + f_3 * pc_x[k] * sid_150[k];

        t_251[k] = f_8 * shd_114[k]
                   + f_3 * pc_y[k] * sid_150[k];

        t_252[k] = f_7 * shd_108[k]
                   + f_3 * pc_z[k] * sid_150[k];
    }

#pragma omp simd aligned(t_253, t_254, t_255, t_256, t_257, pc_x, pc_y, pc_z, shd_111, \
                         shd_117, sip0_76, sip1_76, sid_153, sid_154, \
                         sid_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_253[k] = f_3 * pc_x[k] * sid_153[k];

        t_254[k] = f_3 * pc_x[k] * sid_154[k];

        t_255[k] = f_3 * pc_x[k] * sid_155[k];

        t_256[k] = f_8 * shd_117[k]
                   + f_1 * sip0_76[k]
                   - f_2 * sip1_76[k]
                   + f_3 * pc_y[k] * sid_153[k];

        t_257[k] = f_7 * shd_111[k]
                   + f_3 * pc_z[k] * sid_153[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, t_261, pb_y, pc_y, pc_z, shf0_200, shd_113, \
                         shd_119, shd_120, shf1_200, sip0_77, sip1_77, sid_155, \
                         sid_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = f_8 * shd_119[k]
                   + f_3 * pc_y[k] * sid_155[k];

        t_259[k] = f_7 * shd_113[k]
                   + f_1 * sip0_77[k]
                   - f_2 * sip1_77[k]
                   + f_3 * pc_z[k] * sid_155[k];

        t_260[k] = pb_y[k] * shf0_200[k]
                   - f_4 * pc_y[k] * shf1_200[k];

        t_261[k] = f_5 * shd_120[k]
                   + f_3 * pc_y[k] * sid_156[k];
    }
}

static auto
compute_prim_sif_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t shf0,
                                                          const size_t shd, const size_t shf1,
                                                          const size_t sip0, const size_t sip1,
                                                          const size_t sid, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 2.5 / q;
    const auto f_9 = 1.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *shf0_206 = buffer.data(shf0 + 206);
    const auto *shf0_209 = buffer.data(shf0 + 209);

    const auto *shd_114 = buffer.data(shd + 114);
    const auto *shd_117 = buffer.data(shd + 117);
    const auto *shd_120 = buffer.data(shd + 120);
    const auto *shd_123 = buffer.data(shd + 123);
    const auto *shd_125 = buffer.data(shd + 125);

    const auto *shf1_206 = buffer.data(shf1 + 206);
    const auto *shf1_209 = buffer.data(shf1 + 209);

    const auto *sip0_81 = buffer.data(sip0 + 81);
    const auto *sip0_82 = buffer.data(sip0 + 82);
    const auto *sip0_83 = buffer.data(sip0 + 83);

    const auto *sip1_81 = buffer.data(sip1 + 81);
    const auto *sip1_82 = buffer.data(sip1 + 82);
    const auto *sip1_83 = buffer.data(sip1 + 83);

    const auto *sid_156 = buffer.data(sid + 156);
    const auto *sid_159 = buffer.data(sid + 159);
    const auto *sid_160 = buffer.data(sid + 160);
    const auto *sid_161 = buffer.data(sid + 161);
    const auto *sid_162 = buffer.data(sid + 162);
    const auto *sid_165 = buffer.data(sid + 165);
    const auto *sid_166 = buffer.data(sid + 166);
    const auto *sid_167 = buffer.data(sid + 167);

#pragma omp simd aligned(t_262, t_263, t_264, t_265, pc_x, pc_z, shd_114, sid_156, sid_159, \
                         sid_160, sid_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_262[k] = f_6 * shd_114[k]
                   + f_3 * pc_z[k] * sid_156[k];

        t_263[k] = f_3 * pc_x[k] * sid_159[k];

        t_264[k] = f_3 * pc_x[k] * sid_160[k];

        t_265[k] = f_3 * pc_x[k] * sid_161[k];
    }

#pragma omp simd aligned(t_266, t_267, t_268, t_269, pb_y, pc_y, pc_z, shf0_206, shf0_209, \
                         shd_117, shd_123, shd_125, shf1_206, shf1_209, sid_159, \
                         sid_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_266[k] = pb_y[k] * shf0_206[k]
                   + f_9 * shd_123[k]
                   - f_4 * pc_y[k] * shf1_206[k];

        t_267[k] = f_6 * shd_117[k]
                   + f_3 * pc_z[k] * sid_159[k];

        t_268[k] = f_5 * shd_125[k]
                   + f_3 * pc_y[k] * sid_161[k];

        t_269[k] = pb_y[k] * shf0_209[k]
                   - f_4 * pc_y[k] * shf1_209[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, t_275, pc_x, pc_y, pc_z, shd_120, \
                         sip0_81, sip1_81, sid_162, sid_165, sid_166, \
                         sid_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_1 * sip0_81[k]
                   - f_2 * sip1_81[k]
                   + f_3 * pc_x[k] * sid_162[k];

        t_271[k] = f_3 * pc_y[k] * sid_162[k];

        t_272[k] = f_0 * shd_120[k]
                   + f_3 * pc_z[k] * sid_162[k];

        t_273[k] = f_3 * pc_x[k] * sid_165[k];

        t_274[k] = f_3 * pc_x[k] * sid_166[k];

        t_275[k] = f_3 * pc_x[k] * sid_167[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pc_y, pc_z, shd_123, shd_125, sip0_82, \
                         sip0_83, sip1_82, sip1_83, sid_165, sid_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_1 * sip0_82[k]
                   - f_2 * sip1_82[k]
                   + f_3 * pc_y[k] * sid_165[k];

        t_277[k] = f_0 * shd_123[k]
                   + f_3 * pc_z[k] * sid_165[k];

        t_278[k] = f_3 * pc_y[k] * sid_167[k];

        t_279[k] = f_0 * shd_125[k]
                   + f_1 * sip0_83[k]
                   - f_2 * sip1_83[k]
                   + f_3 * pc_z[k] * sid_167[k];
    }
}

auto
compute_prim_sif_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t shf0, const size_t shd,
                                                   const size_t shf1, const size_t sip0,
                                                   const size_t sip1, const size_t sid,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_sif_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, shf0, shd,
                                                              shf1, sip0, sip1, sid, ncols,
                                                              gamma, p, q);

    compute_prim_sif_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, shf0, shd,
                                                              shf1, sip0, sip1, sid, ncols,
                                                              gamma, p, q);

    compute_prim_sif_three_center_electron_repulsion_0_piece2(buffer, target, pb, pc, shf0, shd,
                                                              shf1, sip0, sip1, sid, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
