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


#include "SimdThreeCenterElectronRepulsionVrrRecSKD.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_skd_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sid0,
                                                          const size_t sip, const size_t sid1,
                                                          const size_t sks0, const size_t sks1,
                                                          const size_t skp, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / q;
    const auto f_1 = 0.5 / gamma;
    const auto f_2 = 0.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 3.0 / q;
    const auto f_6 = 0.5 / q;
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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sid0_0 = buffer.data(sid0 + 0);
    const auto *sid0_3 = buffer.data(sid0 + 3);
    const auto *sid0_5 = buffer.data(sid0 + 5);
    const auto *sid0_9 = buffer.data(sid0 + 9);
    const auto *sid0_12 = buffer.data(sid0 + 12);
    const auto *sid0_17 = buffer.data(sid0 + 17);
    const auto *sid0_18 = buffer.data(sid0 + 18);
    const auto *sid0_21 = buffer.data(sid0 + 21);
    const auto *sid0_30 = buffer.data(sid0 + 30);
    const auto *sid0_35 = buffer.data(sid0 + 35);
    const auto *sid0_36 = buffer.data(sid0 + 36);
    const auto *sid0_39 = buffer.data(sid0 + 39);
    const auto *sid0_54 = buffer.data(sid0 + 54);
    const auto *sid0_59 = buffer.data(sid0 + 59);
    const auto *sid0_60 = buffer.data(sid0 + 60);
    const auto *sid0_63 = buffer.data(sid0 + 63);
    const auto *sid0_84 = buffer.data(sid0 + 84);
    const auto *sid0_89 = buffer.data(sid0 + 89);

    const auto *sip_0 = buffer.data(sip + 0);
    const auto *sip_1 = buffer.data(sip + 1);
    const auto *sip_2 = buffer.data(sip + 2);
    const auto *sip_4 = buffer.data(sip + 4);
    const auto *sip_5 = buffer.data(sip + 5);
    const auto *sip_7 = buffer.data(sip + 7);
    const auto *sip_8 = buffer.data(sip + 8);
    const auto *sip_9 = buffer.data(sip + 9);
    const auto *sip_10 = buffer.data(sip + 10);
    const auto *sip_11 = buffer.data(sip + 11);
    const auto *sip_13 = buffer.data(sip + 13);
    const auto *sip_14 = buffer.data(sip + 14);
    const auto *sip_15 = buffer.data(sip + 15);
    const auto *sip_16 = buffer.data(sip + 16);
    const auto *sip_17 = buffer.data(sip + 17);
    const auto *sip_18 = buffer.data(sip + 18);
    const auto *sip_19 = buffer.data(sip + 19);
    const auto *sip_20 = buffer.data(sip + 20);
    const auto *sip_22 = buffer.data(sip + 22);
    const auto *sip_23 = buffer.data(sip + 23);
    const auto *sip_25 = buffer.data(sip + 25);
    const auto *sip_26 = buffer.data(sip + 26);
    const auto *sip_27 = buffer.data(sip + 27);
    const auto *sip_28 = buffer.data(sip + 28);
    const auto *sip_29 = buffer.data(sip + 29);
    const auto *sip_30 = buffer.data(sip + 30);
    const auto *sip_31 = buffer.data(sip + 31);
    const auto *sip_32 = buffer.data(sip + 32);
    const auto *sip_34 = buffer.data(sip + 34);
    const auto *sip_35 = buffer.data(sip + 35);
    const auto *sip_36 = buffer.data(sip + 36);
    const auto *sip_37 = buffer.data(sip + 37);
    const auto *sip_38 = buffer.data(sip + 38);
    const auto *sip_40 = buffer.data(sip + 40);
    const auto *sip_41 = buffer.data(sip + 41);
    const auto *sip_42 = buffer.data(sip + 42);
    const auto *sip_43 = buffer.data(sip + 43);
    const auto *sip_44 = buffer.data(sip + 44);
    const auto *sip_45 = buffer.data(sip + 45);
    const auto *sip_46 = buffer.data(sip + 46);
    const auto *sip_47 = buffer.data(sip + 47);
    const auto *sip_49 = buffer.data(sip + 49);
    const auto *sip_50 = buffer.data(sip + 50);
    const auto *sip_51 = buffer.data(sip + 51);
    const auto *sip_52 = buffer.data(sip + 52);
    const auto *sip_53 = buffer.data(sip + 53);
    const auto *sip_54 = buffer.data(sip + 54);
    const auto *sip_55 = buffer.data(sip + 55);
    const auto *sip_56 = buffer.data(sip + 56);
    const auto *sip_58 = buffer.data(sip + 58);
    const auto *sip_59 = buffer.data(sip + 59);

    const auto *sid1_0 = buffer.data(sid1 + 0);
    const auto *sid1_3 = buffer.data(sid1 + 3);
    const auto *sid1_5 = buffer.data(sid1 + 5);
    const auto *sid1_9 = buffer.data(sid1 + 9);
    const auto *sid1_12 = buffer.data(sid1 + 12);
    const auto *sid1_17 = buffer.data(sid1 + 17);
    const auto *sid1_18 = buffer.data(sid1 + 18);
    const auto *sid1_21 = buffer.data(sid1 + 21);
    const auto *sid1_30 = buffer.data(sid1 + 30);
    const auto *sid1_35 = buffer.data(sid1 + 35);
    const auto *sid1_36 = buffer.data(sid1 + 36);
    const auto *sid1_39 = buffer.data(sid1 + 39);
    const auto *sid1_54 = buffer.data(sid1 + 54);
    const auto *sid1_59 = buffer.data(sid1 + 59);
    const auto *sid1_60 = buffer.data(sid1 + 60);
    const auto *sid1_63 = buffer.data(sid1 + 63);
    const auto *sid1_84 = buffer.data(sid1 + 84);
    const auto *sid1_89 = buffer.data(sid1 + 89);

    const auto *sks0_0 = buffer.data(sks0 + 0);
    const auto *sks0_1 = buffer.data(sks0 + 1);
    const auto *sks0_2 = buffer.data(sks0 + 2);
    const auto *sks0_3 = buffer.data(sks0 + 3);
    const auto *sks0_5 = buffer.data(sks0 + 5);
    const auto *sks0_6 = buffer.data(sks0 + 6);
    const auto *sks0_7 = buffer.data(sks0 + 7);
    const auto *sks0_8 = buffer.data(sks0 + 8);
    const auto *sks0_9 = buffer.data(sks0 + 9);
    const auto *sks0_10 = buffer.data(sks0 + 10);
    const auto *sks0_11 = buffer.data(sks0 + 11);
    const auto *sks0_12 = buffer.data(sks0 + 12);
    const auto *sks0_13 = buffer.data(sks0 + 13);
    const auto *sks0_14 = buffer.data(sks0 + 14);
    const auto *sks0_15 = buffer.data(sks0 + 15);
    const auto *sks0_16 = buffer.data(sks0 + 16);
    const auto *sks0_17 = buffer.data(sks0 + 17);
    const auto *sks0_18 = buffer.data(sks0 + 18);
    const auto *sks0_19 = buffer.data(sks0 + 19);

    const auto *sks1_0 = buffer.data(sks1 + 0);
    const auto *sks1_1 = buffer.data(sks1 + 1);
    const auto *sks1_2 = buffer.data(sks1 + 2);
    const auto *sks1_3 = buffer.data(sks1 + 3);
    const auto *sks1_5 = buffer.data(sks1 + 5);
    const auto *sks1_6 = buffer.data(sks1 + 6);
    const auto *sks1_7 = buffer.data(sks1 + 7);
    const auto *sks1_8 = buffer.data(sks1 + 8);
    const auto *sks1_9 = buffer.data(sks1 + 9);
    const auto *sks1_10 = buffer.data(sks1 + 10);
    const auto *sks1_11 = buffer.data(sks1 + 11);
    const auto *sks1_12 = buffer.data(sks1 + 12);
    const auto *sks1_13 = buffer.data(sks1 + 13);
    const auto *sks1_14 = buffer.data(sks1 + 14);
    const auto *sks1_15 = buffer.data(sks1 + 15);
    const auto *sks1_16 = buffer.data(sks1 + 16);
    const auto *sks1_17 = buffer.data(sks1 + 17);
    const auto *sks1_18 = buffer.data(sks1 + 18);
    const auto *sks1_19 = buffer.data(sks1 + 19);

    const auto *skp_0 = buffer.data(skp + 0);
    const auto *skp_1 = buffer.data(skp + 1);
    const auto *skp_2 = buffer.data(skp + 2);
    const auto *skp_4 = buffer.data(skp + 4);
    const auto *skp_5 = buffer.data(skp + 5);
    const auto *skp_7 = buffer.data(skp + 7);
    const auto *skp_8 = buffer.data(skp + 8);
    const auto *skp_9 = buffer.data(skp + 9);
    const auto *skp_10 = buffer.data(skp + 10);
    const auto *skp_11 = buffer.data(skp + 11);
    const auto *skp_13 = buffer.data(skp + 13);
    const auto *skp_14 = buffer.data(skp + 14);
    const auto *skp_15 = buffer.data(skp + 15);
    const auto *skp_16 = buffer.data(skp + 16);
    const auto *skp_17 = buffer.data(skp + 17);
    const auto *skp_18 = buffer.data(skp + 18);
    const auto *skp_19 = buffer.data(skp + 19);
    const auto *skp_20 = buffer.data(skp + 20);
    const auto *skp_22 = buffer.data(skp + 22);
    const auto *skp_23 = buffer.data(skp + 23);
    const auto *skp_25 = buffer.data(skp + 25);
    const auto *skp_26 = buffer.data(skp + 26);
    const auto *skp_27 = buffer.data(skp + 27);
    const auto *skp_28 = buffer.data(skp + 28);
    const auto *skp_29 = buffer.data(skp + 29);
    const auto *skp_30 = buffer.data(skp + 30);
    const auto *skp_31 = buffer.data(skp + 31);
    const auto *skp_32 = buffer.data(skp + 32);
    const auto *skp_34 = buffer.data(skp + 34);
    const auto *skp_35 = buffer.data(skp + 35);
    const auto *skp_36 = buffer.data(skp + 36);
    const auto *skp_37 = buffer.data(skp + 37);
    const auto *skp_38 = buffer.data(skp + 38);
    const auto *skp_40 = buffer.data(skp + 40);
    const auto *skp_41 = buffer.data(skp + 41);
    const auto *skp_42 = buffer.data(skp + 42);
    const auto *skp_43 = buffer.data(skp + 43);
    const auto *skp_44 = buffer.data(skp + 44);
    const auto *skp_45 = buffer.data(skp + 45);
    const auto *skp_46 = buffer.data(skp + 46);
    const auto *skp_47 = buffer.data(skp + 47);
    const auto *skp_49 = buffer.data(skp + 49);
    const auto *skp_50 = buffer.data(skp + 50);
    const auto *skp_51 = buffer.data(skp + 51);
    const auto *skp_52 = buffer.data(skp + 52);
    const auto *skp_53 = buffer.data(skp + 53);
    const auto *skp_54 = buffer.data(skp + 54);
    const auto *skp_55 = buffer.data(skp + 55);
    const auto *skp_56 = buffer.data(skp + 56);
    const auto *skp_58 = buffer.data(skp + 58);
    const auto *skp_59 = buffer.data(skp + 59);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pc_x, pc_y, sip_0, sip_1, sip_2, sks0_0, \
                         sks1_0, skp_0, skp_1, skp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sip_0[k]
                 + f_1 * sks0_0[k]
                 - f_2 * sks1_0[k]
                 + f_3 * pc_x[k] * skp_0[k];

        t_1[k] = f_0 * sip_1[k]
                 + f_3 * pc_x[k] * skp_1[k];

        t_2[k] = f_0 * sip_2[k]
                 + f_3 * pc_x[k] * skp_2[k];

        t_3[k] = f_1 * sks0_0[k]
                 - f_2 * sks1_0[k]
                 + f_3 * pc_y[k] * skp_1[k];

        t_4[k] = f_3 * pc_y[k] * skp_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pb_y, pc_x, pc_y, pc_z, sid0_0, sip_4, sid1_0, sks0_0, \
                         sks1_0, skp_2, skp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * sks0_0[k]
                 - f_2 * sks1_0[k]
                 + f_3 * pc_z[k] * skp_2[k];

        t_6[k] = pb_y[k] * sid0_0[k]
                 - f_4 * pc_y[k] * sid1_0[k];

        t_7[k] = f_5 * sip_4[k]
                 + f_3 * pc_x[k] * skp_4[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pb_y, pc_x, pc_y, sid0_5, sip_1, sip_2, sip_5, \
                         sid1_5, sks0_1, sks1_1, skp_4, skp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_5 * sip_5[k]
                 + f_3 * pc_x[k] * skp_5[k];

        t_9[k] = f_6 * sip_1[k]
                 + f_1 * sks0_1[k]
                 - f_2 * sks1_1[k]
                 + f_3 * pc_y[k] * skp_4[k];

        t_10[k] = f_6 * sip_2[k]
                  + f_3 * pc_y[k] * skp_5[k];

        t_11[k] = pb_y[k] * sid0_5[k]
                  - f_4 * pc_y[k] * sid1_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pb_z, pc_x, pc_z, sid0_0, sid0_3, sip_7, \
                         sip_8, sid1_0, sid1_3, skp_7, skp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pb_z[k] * sid0_0[k]
                  - f_4 * pc_z[k] * sid1_0[k];

        t_13[k] = f_5 * sip_7[k]
                  + f_3 * pc_x[k] * skp_7[k];

        t_14[k] = f_5 * sip_8[k]
                  + f_3 * pc_x[k] * skp_8[k];

        t_15[k] = pb_z[k] * sid0_3[k]
                  - f_4 * pc_z[k] * sid1_3[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pc_x, pc_y, pc_z, sip_2, sip_9, sks0_2, sks0_3, \
                         sks1_2, sks1_3, skp_8, skp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * pc_y[k] * skp_8[k];

        t_17[k] = f_6 * sip_2[k]
                  + f_1 * sks0_2[k]
                  - f_2 * sks1_2[k]
                  + f_3 * pc_z[k] * skp_8[k];

        t_18[k] = f_7 * sip_9[k]
                  + f_1 * sks0_3[k]
                  - f_2 * sks1_3[k]
                  + f_3 * pc_x[k] * skp_9[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pc_x, pc_y, pc_z, sip_4, sip_5, sip_10, \
                         sip_11, sks0_3, sks1_3, skp_10, skp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_7 * sip_10[k]
                  + f_3 * pc_x[k] * skp_10[k];

        t_20[k] = f_7 * sip_11[k]
                  + f_3 * pc_x[k] * skp_11[k];

        t_21[k] = f_8 * sip_4[k]
                  + f_1 * sks0_3[k]
                  - f_2 * sks1_3[k]
                  + f_3 * pc_y[k] * skp_10[k];

        t_22[k] = f_8 * sip_5[k]
                  + f_3 * pc_y[k] * skp_11[k];

        t_23[k] = f_1 * sks0_3[k]
                  - f_2 * sks1_3[k]
                  + f_3 * pc_z[k] * skp_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pb_y, pc_x, pc_y, sid0_12, sip_13, sip_14, sid1_12, \
                         skp_13, skp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pb_y[k] * sid0_12[k]
                  - f_4 * pc_y[k] * sid1_12[k];

        t_25[k] = f_7 * sip_13[k]
                  + f_3 * pc_x[k] * skp_13[k];

        t_26[k] = f_7 * sip_14[k]
                  + f_3 * pc_x[k] * skp_14[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pb_y, pb_z, pc_y, pc_z, sid0_9, sid0_17, sip_8, \
                         sid1_9, sid1_17, skp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = pb_z[k] * sid0_9[k]
                  - f_4 * pc_z[k] * sid1_9[k];

        t_28[k] = f_6 * sip_8[k]
                  + f_3 * pc_y[k] * skp_14[k];

        t_29[k] = pb_y[k] * sid0_17[k]
                  - f_4 * pc_y[k] * sid1_17[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pc_x, pc_y, sip_15, sip_16, sip_17, \
                         sks0_5, sks1_5, skp_15, skp_16, skp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_7 * sip_15[k]
                  + f_1 * sks0_5[k]
                  - f_2 * sks1_5[k]
                  + f_3 * pc_x[k] * skp_15[k];

        t_31[k] = f_7 * sip_16[k]
                  + f_3 * pc_x[k] * skp_16[k];

        t_32[k] = f_7 * sip_17[k]
                  + f_3 * pc_x[k] * skp_17[k];

        t_33[k] = f_1 * sks0_5[k]
                  - f_2 * sks1_5[k]
                  + f_3 * pc_y[k] * skp_16[k];

        t_34[k] = f_3 * pc_y[k] * skp_17[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pc_x, pc_z, sip_8, sip_18, sip_19, sks0_5, sks0_6, \
                         sks1_5, sks1_6, skp_17, skp_18, skp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_8 * sip_8[k]
                  + f_1 * sks0_5[k]
                  - f_2 * sks1_5[k]
                  + f_3 * pc_z[k] * skp_17[k];

        t_36[k] = f_9 * sip_18[k]
                  + f_1 * sks0_6[k]
                  - f_2 * sks1_6[k]
                  + f_3 * pc_x[k] * skp_18[k];

        t_37[k] = f_9 * sip_19[k]
                  + f_3 * pc_x[k] * skp_19[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pc_x, pc_y, pc_z, sip_10, sip_11, sip_20, \
                         sks0_6, sks1_6, skp_19, skp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_9 * sip_20[k]
                  + f_3 * pc_x[k] * skp_20[k];

        t_39[k] = f_10 * sip_10[k]
                  + f_1 * sks0_6[k]
                  - f_2 * sks1_6[k]
                  + f_3 * pc_y[k] * skp_19[k];

        t_40[k] = f_10 * sip_11[k]
                  + f_3 * pc_y[k] * skp_20[k];

        t_41[k] = f_1 * sks0_6[k]
                  - f_2 * sks1_6[k]
                  + f_3 * pc_z[k] * skp_20[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pb_z, pc_x, pc_z, sid0_18, sid0_21, sip_22, \
                         sip_23, sid1_18, sid1_21, skp_22, skp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_z[k] * sid0_18[k]
                  - f_4 * pc_z[k] * sid1_18[k];

        t_43[k] = f_9 * sip_22[k]
                  + f_3 * pc_x[k] * skp_22[k];

        t_44[k] = f_9 * sip_23[k]
                  + f_3 * pc_x[k] * skp_23[k];

        t_45[k] = pb_z[k] * sid0_21[k]
                  - f_4 * pc_z[k] * sid1_21[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, pb_y, pc_y, pc_z, sid0_30, sip_11, sip_14, sid1_30, \
                         sks0_7, sks1_7, skp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_8 * sip_14[k]
                  + f_3 * pc_y[k] * skp_23[k];

        t_47[k] = f_6 * sip_11[k]
                  + f_1 * sks0_7[k]
                  - f_2 * sks1_7[k]
                  + f_3 * pc_z[k] * skp_23[k];

        t_48[k] = pb_y[k] * sid0_30[k]
                  - f_4 * pc_y[k] * sid1_30[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pc_x, pc_y, sip_16, sip_17, sip_25, sip_26, \
                         sks0_8, sks1_8, skp_25, skp_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_9 * sip_25[k]
                  + f_3 * pc_x[k] * skp_25[k];

        t_50[k] = f_9 * sip_26[k]
                  + f_3 * pc_x[k] * skp_26[k];

        t_51[k] = f_6 * sip_16[k]
                  + f_1 * sks0_8[k]
                  - f_2 * sks1_8[k]
                  + f_3 * pc_y[k] * skp_25[k];

        t_52[k] = f_6 * sip_17[k]
                  + f_3 * pc_y[k] * skp_26[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pb_y, pc_x, pc_y, sid0_35, sip_27, sip_28, sid1_35, \
                         sks0_9, sks1_9, skp_27, skp_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = pb_y[k] * sid0_35[k]
                  - f_4 * pc_y[k] * sid1_35[k];

        t_54[k] = f_9 * sip_27[k]
                  + f_1 * sks0_9[k]
                  - f_2 * sks1_9[k]
                  + f_3 * pc_x[k] * skp_27[k];

        t_55[k] = f_9 * sip_28[k]
                  + f_3 * pc_x[k] * skp_28[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pc_x, pc_y, pc_z, sip_17, sip_29, sks0_9, \
                         sks1_9, skp_28, skp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_9 * sip_29[k]
                  + f_3 * pc_x[k] * skp_29[k];

        t_57[k] = f_1 * sks0_9[k]
                  - f_2 * sks1_9[k]
                  + f_3 * pc_y[k] * skp_28[k];

        t_58[k] = f_3 * pc_y[k] * skp_29[k];

        t_59[k] = f_10 * sip_17[k]
                  + f_1 * sks0_9[k]
                  - f_2 * sks1_9[k]
                  + f_3 * pc_z[k] * skp_29[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pc_x, pc_y, sip_19, sip_30, sip_31, sip_32, \
                         sks0_10, sks1_10, skp_30, skp_31, skp_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_10 * sip_30[k]
                  + f_1 * sks0_10[k]
                  - f_2 * sks1_10[k]
                  + f_3 * pc_x[k] * skp_30[k];

        t_61[k] = f_10 * sip_31[k]
                  + f_3 * pc_x[k] * skp_31[k];

        t_62[k] = f_10 * sip_32[k]
                  + f_3 * pc_x[k] * skp_32[k];

        t_63[k] = f_9 * sip_19[k]
                  + f_1 * sks0_10[k]
                  - f_2 * sks1_10[k]
                  + f_3 * pc_y[k] * skp_31[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pb_z, pc_x, pc_y, pc_z, sid0_36, sip_20, \
                         sip_34, sid1_36, sks0_10, sks1_10, skp_32, \
                         skp_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_9 * sip_20[k]
                  + f_3 * pc_y[k] * skp_32[k];

        t_65[k] = f_1 * sks0_10[k]
                  - f_2 * sks1_10[k]
                  + f_3 * pc_z[k] * skp_32[k];

        t_66[k] = pb_z[k] * sid0_36[k]
                  - f_4 * pc_z[k] * sid1_36[k];

        t_67[k] = f_10 * sip_34[k]
                  + f_3 * pc_x[k] * skp_34[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pb_z, pc_x, pc_y, pc_z, sid0_39, sip_20, \
                         sip_23, sip_35, sid1_39, sks0_11, sks1_11, \
                         skp_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_10 * sip_35[k]
                  + f_3 * pc_x[k] * skp_35[k];

        t_69[k] = pb_z[k] * sid0_39[k]
                  - f_4 * pc_z[k] * sid1_39[k];

        t_70[k] = f_10 * sip_23[k]
                  + f_3 * pc_y[k] * skp_35[k];

        t_71[k] = f_6 * sip_20[k]
                  + f_1 * sks0_11[k]
                  - f_2 * sks1_11[k]
                  + f_3 * pc_z[k] * skp_35[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pc_x, pc_y, sip_25, sip_36, sip_37, sip_38, \
                         sks0_12, sks1_12, skp_36, skp_37, skp_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_10 * sip_36[k]
                  + f_1 * sks0_12[k]
                  - f_2 * sks1_12[k]
                  + f_3 * pc_x[k] * skp_36[k];

        t_73[k] = f_10 * sip_37[k]
                  + f_3 * pc_x[k] * skp_37[k];

        t_74[k] = f_10 * sip_38[k]
                  + f_3 * pc_x[k] * skp_38[k];

        t_75[k] = f_8 * sip_25[k]
                  + f_1 * sks0_12[k]
                  - f_2 * sks1_12[k]
                  + f_3 * pc_y[k] * skp_37[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, pb_y, pc_y, pc_z, sid0_54, sip_23, sip_26, sid1_54, \
                         sks0_12, sks1_12, skp_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_8 * sip_26[k]
                  + f_3 * pc_y[k] * skp_38[k];

        t_77[k] = f_8 * sip_23[k]
                  + f_1 * sks0_12[k]
                  - f_2 * sks1_12[k]
                  + f_3 * pc_z[k] * skp_38[k];

        t_78[k] = pb_y[k] * sid0_54[k]
                  - f_4 * pc_y[k] * sid1_54[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pc_x, pc_y, sip_28, sip_29, sip_40, sip_41, \
                         sks0_13, sks1_13, skp_40, skp_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_10 * sip_40[k]
                  + f_3 * pc_x[k] * skp_40[k];

        t_80[k] = f_10 * sip_41[k]
                  + f_3 * pc_x[k] * skp_41[k];

        t_81[k] = f_6 * sip_28[k]
                  + f_1 * sks0_13[k]
                  - f_2 * sks1_13[k]
                  + f_3 * pc_y[k] * skp_40[k];

        t_82[k] = f_6 * sip_29[k]
                  + f_3 * pc_y[k] * skp_41[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, pb_y, pc_x, pc_y, sid0_59, sip_42, sip_43, sid1_59, \
                         sks0_14, sks1_14, skp_42, skp_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = pb_y[k] * sid0_59[k]
                  - f_4 * pc_y[k] * sid1_59[k];

        t_84[k] = f_10 * sip_42[k]
                  + f_1 * sks0_14[k]
                  - f_2 * sks1_14[k]
                  + f_3 * pc_x[k] * skp_42[k];

        t_85[k] = f_10 * sip_43[k]
                  + f_3 * pc_x[k] * skp_43[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pc_x, pc_y, pc_z, sip_29, sip_44, sks0_14, \
                         sks1_14, skp_43, skp_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_10 * sip_44[k]
                  + f_3 * pc_x[k] * skp_44[k];

        t_87[k] = f_1 * sks0_14[k]
                  - f_2 * sks1_14[k]
                  + f_3 * pc_y[k] * skp_43[k];

        t_88[k] = f_3 * pc_y[k] * skp_44[k];

        t_89[k] = f_9 * sip_29[k]
                  + f_1 * sks0_14[k]
                  - f_2 * sks1_14[k]
                  + f_3 * pc_z[k] * skp_44[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, pc_x, pc_y, sip_31, sip_45, sip_46, sip_47, \
                         sks0_15, sks1_15, skp_45, skp_46, skp_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_8 * sip_45[k]
                  + f_1 * sks0_15[k]
                  - f_2 * sks1_15[k]
                  + f_3 * pc_x[k] * skp_45[k];

        t_91[k] = f_8 * sip_46[k]
                  + f_3 * pc_x[k] * skp_46[k];

        t_92[k] = f_8 * sip_47[k]
                  + f_3 * pc_x[k] * skp_47[k];

        t_93[k] = f_7 * sip_31[k]
                  + f_1 * sks0_15[k]
                  - f_2 * sks1_15[k]
                  + f_3 * pc_y[k] * skp_46[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, pb_z, pc_x, pc_y, pc_z, sid0_60, sip_32, \
                         sip_49, sid1_60, sks0_15, sks1_15, skp_47, \
                         skp_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_7 * sip_32[k]
                  + f_3 * pc_y[k] * skp_47[k];

        t_95[k] = f_1 * sks0_15[k]
                  - f_2 * sks1_15[k]
                  + f_3 * pc_z[k] * skp_47[k];

        t_96[k] = pb_z[k] * sid0_60[k]
                  - f_4 * pc_z[k] * sid1_60[k];

        t_97[k] = f_8 * sip_49[k]
                  + f_3 * pc_x[k] * skp_49[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, pb_z, pc_x, pc_y, pc_z, sid0_63, sip_32, \
                         sip_35, sip_50, sid1_63, sks0_16, sks1_16, \
                         skp_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_8 * sip_50[k]
                  + f_3 * pc_x[k] * skp_50[k];

        t_99[k] = pb_z[k] * sid0_63[k]
                  - f_4 * pc_z[k] * sid1_63[k];

        t_100[k] = f_9 * sip_35[k]
                   + f_3 * pc_y[k] * skp_50[k];

        t_101[k] = f_6 * sip_32[k]
                   + f_1 * sks0_16[k]
                   - f_2 * sks1_16[k]
                   + f_3 * pc_z[k] * skp_50[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, pc_x, pc_y, sip_37, sip_51, sip_52, \
                         sip_53, sks0_17, sks1_17, skp_51, skp_52, \
                         skp_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_8 * sip_51[k]
                   + f_1 * sks0_17[k]
                   - f_2 * sks1_17[k]
                   + f_3 * pc_x[k] * skp_51[k];

        t_103[k] = f_8 * sip_52[k]
                   + f_3 * pc_x[k] * skp_52[k];

        t_104[k] = f_8 * sip_53[k]
                   + f_3 * pc_x[k] * skp_53[k];

        t_105[k] = f_10 * sip_37[k]
                   + f_1 * sks0_17[k]
                   - f_2 * sks1_17[k]
                   + f_3 * pc_y[k] * skp_52[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, pc_x, pc_y, pc_z, sip_35, sip_38, sip_54, \
                         sks0_17, sks0_18, sks1_17, sks1_18, skp_53, \
                         skp_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_10 * sip_38[k]
                   + f_3 * pc_y[k] * skp_53[k];

        t_107[k] = f_8 * sip_35[k]
                   + f_1 * sks0_17[k]
                   - f_2 * sks1_17[k]
                   + f_3 * pc_z[k] * skp_53[k];

        t_108[k] = f_8 * sip_54[k]
                   + f_1 * sks0_18[k]
                   - f_2 * sks1_18[k]
                   + f_3 * pc_x[k] * skp_54[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pc_x, pc_y, sip_40, sip_41, sip_55, \
                         sip_56, sks0_18, sks1_18, skp_55, skp_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_8 * sip_55[k]
                   + f_3 * pc_x[k] * skp_55[k];

        t_110[k] = f_8 * sip_56[k]
                   + f_3 * pc_x[k] * skp_56[k];

        t_111[k] = f_8 * sip_40[k]
                   + f_1 * sks0_18[k]
                   - f_2 * sks1_18[k]
                   + f_3 * pc_y[k] * skp_55[k];

        t_112[k] = f_8 * sip_41[k]
                   + f_3 * pc_y[k] * skp_56[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, pb_y, pc_x, pc_y, pc_z, sid0_84, sip_38, sip_58, \
                         sid1_84, sks0_18, sks1_18, skp_56, skp_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_10 * sip_38[k]
                   + f_1 * sks0_18[k]
                   - f_2 * sks1_18[k]
                   + f_3 * pc_z[k] * skp_56[k];

        t_114[k] = pb_y[k] * sid0_84[k]
                   - f_4 * pc_y[k] * sid1_84[k];

        t_115[k] = f_8 * sip_58[k]
                   + f_3 * pc_x[k] * skp_58[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, pb_y, pc_x, pc_y, sid0_89, sip_43, \
                         sip_44, sip_59, sid1_89, sks0_19, sks1_19, skp_58, \
                         skp_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_8 * sip_59[k]
                   + f_3 * pc_x[k] * skp_59[k];

        t_117[k] = f_6 * sip_43[k]
                   + f_1 * sks0_19[k]
                   - f_2 * sks1_19[k]
                   + f_3 * pc_y[k] * skp_58[k];

        t_118[k] = f_6 * sip_44[k]
                   + f_3 * pc_y[k] * skp_59[k];

        t_119[k] = pb_y[k] * sid0_89[k]
                   - f_4 * pc_y[k] * sid1_89[k];
    }
}

static auto
compute_prim_skd_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sid0,
                                                          const size_t sip, const size_t sid1,
                                                          const size_t sks0, const size_t sks1,
                                                          const size_t skp, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / q;
    const auto f_1 = 0.5 / gamma;
    const auto f_2 = 0.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 3.0 / q;
    const auto f_6 = 0.5 / q;
    const auto f_7 = 2.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 2.0 / q;
    const auto f_10 = 1.5 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sid0_90 = buffer.data(sid0 + 90);
    const auto *sid0_120 = buffer.data(sid0 + 120);
    const auto *sid0_126 = buffer.data(sid0 + 126);
    const auto *sid0_129 = buffer.data(sid0 + 129);
    const auto *sid0_131 = buffer.data(sid0 + 131);
    const auto *sid0_135 = buffer.data(sid0 + 135);
    const auto *sid0_137 = buffer.data(sid0 + 137);
    const auto *sid0_138 = buffer.data(sid0 + 138);
    const auto *sid0_141 = buffer.data(sid0 + 141);
    const auto *sid0_143 = buffer.data(sid0 + 143);
    const auto *sid0_144 = buffer.data(sid0 + 144);
    const auto *sid0_147 = buffer.data(sid0 + 147);
    const auto *sid0_149 = buffer.data(sid0 + 149);
    const auto *sid0_150 = buffer.data(sid0 + 150);
    const auto *sid0_153 = buffer.data(sid0 + 153);
    const auto *sid0_155 = buffer.data(sid0 + 155);
    const auto *sid0_159 = buffer.data(sid0 + 159);
    const auto *sid0_161 = buffer.data(sid0 + 161);
    const auto *sid0_162 = buffer.data(sid0 + 162);
    const auto *sid0_165 = buffer.data(sid0 + 165);
    const auto *sid0_167 = buffer.data(sid0 + 167);

    const auto *sip_44 = buffer.data(sip + 44);
    const auto *sip_47 = buffer.data(sip + 47);
    const auto *sip_50 = buffer.data(sip + 50);
    const auto *sip_53 = buffer.data(sip + 53);
    const auto *sip_56 = buffer.data(sip + 56);
    const auto *sip_59 = buffer.data(sip + 59);
    const auto *sip_60 = buffer.data(sip + 60);
    const auto *sip_61 = buffer.data(sip + 61);
    const auto *sip_62 = buffer.data(sip + 62);
    const auto *sip_63 = buffer.data(sip + 63);
    const auto *sip_64 = buffer.data(sip + 64);
    const auto *sip_65 = buffer.data(sip + 65);
    const auto *sip_67 = buffer.data(sip + 67);
    const auto *sip_68 = buffer.data(sip + 68);
    const auto *sip_69 = buffer.data(sip + 69);
    const auto *sip_70 = buffer.data(sip + 70);
    const auto *sip_71 = buffer.data(sip + 71);
    const auto *sip_72 = buffer.data(sip + 72);
    const auto *sip_73 = buffer.data(sip + 73);
    const auto *sip_74 = buffer.data(sip + 74);
    const auto *sip_75 = buffer.data(sip + 75);
    const auto *sip_76 = buffer.data(sip + 76);
    const auto *sip_77 = buffer.data(sip + 77);
    const auto *sip_79 = buffer.data(sip + 79);
    const auto *sip_80 = buffer.data(sip + 80);
    const auto *sip_81 = buffer.data(sip + 81);
    const auto *sip_82 = buffer.data(sip + 82);
    const auto *sip_83 = buffer.data(sip + 83);

    const auto *sid1_90 = buffer.data(sid1 + 90);
    const auto *sid1_120 = buffer.data(sid1 + 120);
    const auto *sid1_126 = buffer.data(sid1 + 126);
    const auto *sid1_129 = buffer.data(sid1 + 129);
    const auto *sid1_131 = buffer.data(sid1 + 131);
    const auto *sid1_135 = buffer.data(sid1 + 135);
    const auto *sid1_137 = buffer.data(sid1 + 137);
    const auto *sid1_138 = buffer.data(sid1 + 138);
    const auto *sid1_141 = buffer.data(sid1 + 141);
    const auto *sid1_143 = buffer.data(sid1 + 143);
    const auto *sid1_144 = buffer.data(sid1 + 144);
    const auto *sid1_147 = buffer.data(sid1 + 147);
    const auto *sid1_149 = buffer.data(sid1 + 149);
    const auto *sid1_150 = buffer.data(sid1 + 150);
    const auto *sid1_153 = buffer.data(sid1 + 153);
    const auto *sid1_155 = buffer.data(sid1 + 155);
    const auto *sid1_159 = buffer.data(sid1 + 159);
    const auto *sid1_161 = buffer.data(sid1 + 161);
    const auto *sid1_162 = buffer.data(sid1 + 162);
    const auto *sid1_165 = buffer.data(sid1 + 165);
    const auto *sid1_167 = buffer.data(sid1 + 167);

    const auto *sks0_20 = buffer.data(sks0 + 20);
    const auto *sks0_28 = buffer.data(sks0 + 28);
    const auto *sks0_29 = buffer.data(sks0 + 29);
    const auto *sks0_30 = buffer.data(sks0 + 30);
    const auto *sks0_31 = buffer.data(sks0 + 31);
    const auto *sks0_32 = buffer.data(sks0 + 32);
    const auto *sks0_33 = buffer.data(sks0 + 33);
    const auto *sks0_35 = buffer.data(sks0 + 35);

    const auto *sks1_20 = buffer.data(sks1 + 20);
    const auto *sks1_28 = buffer.data(sks1 + 28);
    const auto *sks1_29 = buffer.data(sks1 + 29);
    const auto *sks1_30 = buffer.data(sks1 + 30);
    const auto *sks1_31 = buffer.data(sks1 + 31);
    const auto *sks1_32 = buffer.data(sks1 + 32);
    const auto *sks1_33 = buffer.data(sks1 + 33);
    const auto *sks1_35 = buffer.data(sks1 + 35);

    const auto *skp_60 = buffer.data(skp + 60);
    const auto *skp_61 = buffer.data(skp + 61);
    const auto *skp_62 = buffer.data(skp + 62);
    const auto *skp_64 = buffer.data(skp + 64);
    const auto *skp_65 = buffer.data(skp + 65);
    const auto *skp_67 = buffer.data(skp + 67);
    const auto *skp_68 = buffer.data(skp + 68);
    const auto *skp_70 = buffer.data(skp + 70);
    const auto *skp_71 = buffer.data(skp + 71);
    const auto *skp_73 = buffer.data(skp + 73);
    const auto *skp_74 = buffer.data(skp + 74);
    const auto *skp_76 = buffer.data(skp + 76);
    const auto *skp_77 = buffer.data(skp + 77);
    const auto *skp_79 = buffer.data(skp + 79);
    const auto *skp_80 = buffer.data(skp + 80);
    const auto *skp_82 = buffer.data(skp + 82);
    const auto *skp_83 = buffer.data(skp + 83);
    const auto *skp_84 = buffer.data(skp + 84);
    const auto *skp_85 = buffer.data(skp + 85);
    const auto *skp_86 = buffer.data(skp + 86);
    const auto *skp_88 = buffer.data(skp + 88);
    const auto *skp_89 = buffer.data(skp + 89);
    const auto *skp_90 = buffer.data(skp + 90);
    const auto *skp_91 = buffer.data(skp + 91);
    const auto *skp_92 = buffer.data(skp + 92);
    const auto *skp_93 = buffer.data(skp + 93);
    const auto *skp_94 = buffer.data(skp + 94);
    const auto *skp_95 = buffer.data(skp + 95);
    const auto *skp_96 = buffer.data(skp + 96);
    const auto *skp_97 = buffer.data(skp + 97);
    const auto *skp_98 = buffer.data(skp + 98);
    const auto *skp_99 = buffer.data(skp + 99);
    const auto *skp_100 = buffer.data(skp + 100);
    const auto *skp_101 = buffer.data(skp + 101);
    const auto *skp_103 = buffer.data(skp + 103);
    const auto *skp_104 = buffer.data(skp + 104);
    const auto *skp_105 = buffer.data(skp + 105);
    const auto *skp_106 = buffer.data(skp + 106);
    const auto *skp_107 = buffer.data(skp + 107);

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, pc_x, pc_y, sip_60, sip_61, \
                         sip_62, sks0_20, sks1_20, skp_60, skp_61, \
                         skp_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_8 * sip_60[k]
                   + f_1 * sks0_20[k]
                   - f_2 * sks1_20[k]
                   + f_3 * pc_x[k] * skp_60[k];

        t_121[k] = f_8 * sip_61[k]
                   + f_3 * pc_x[k] * skp_61[k];

        t_122[k] = f_8 * sip_62[k]
                   + f_3 * pc_x[k] * skp_62[k];

        t_123[k] = f_1 * sks0_20[k]
                   - f_2 * sks1_20[k]
                   + f_3 * pc_y[k] * skp_61[k];

        t_124[k] = f_3 * pc_y[k] * skp_62[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, pb_x, pc_x, pc_z, sid0_126, sip_44, sip_63, \
                         sip_64, sid1_126, sks0_20, sks1_20, skp_62, \
                         skp_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_7 * sip_44[k]
                   + f_1 * sks0_20[k]
                   - f_2 * sks1_20[k]
                   + f_3 * pc_z[k] * skp_62[k];

        t_126[k] = pb_x[k] * sid0_126[k]
                   + f_8 * sip_63[k]
                   - f_4 * pc_x[k] * sid1_126[k];

        t_127[k] = f_6 * sip_64[k]
                   + f_3 * pc_x[k] * skp_64[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, pb_x, pc_x, pc_y, sid0_129, sid0_131, \
                         sip_47, sip_65, sid1_129, sid1_131, skp_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_6 * sip_65[k]
                   + f_3 * pc_x[k] * skp_65[k];

        t_129[k] = pb_x[k] * sid0_129[k]
                   - f_4 * pc_x[k] * sid1_129[k];

        t_130[k] = f_5 * sip_47[k]
                   + f_3 * pc_y[k] * skp_65[k];

        t_131[k] = pb_x[k] * sid0_131[k]
                   - f_4 * pc_x[k] * sid1_131[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pb_x, pb_z, pc_x, pc_z, sid0_90, \
                         sid0_135, sip_67, sip_68, sid1_90, sid1_135, skp_67, \
                         skp_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = pb_z[k] * sid0_90[k]
                   - f_4 * pc_z[k] * sid1_90[k];

        t_133[k] = f_6 * sip_67[k]
                   + f_3 * pc_x[k] * skp_67[k];

        t_134[k] = f_6 * sip_68[k]
                   + f_3 * pc_x[k] * skp_68[k];

        t_135[k] = pb_x[k] * sid0_135[k]
                   - f_4 * pc_x[k] * sid1_135[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, pb_x, pc_x, pc_y, sid0_137, sid0_138, \
                         sip_50, sip_69, sip_70, sid1_137, sid1_138, skp_68, \
                         skp_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_7 * sip_50[k]
                   + f_3 * pc_y[k] * skp_68[k];

        t_137[k] = pb_x[k] * sid0_137[k]
                   - f_4 * pc_x[k] * sid1_137[k];

        t_138[k] = pb_x[k] * sid0_138[k]
                   + f_8 * sip_69[k]
                   - f_4 * pc_x[k] * sid1_138[k];

        t_139[k] = f_6 * sip_70[k]
                   + f_3 * pc_x[k] * skp_70[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, pb_x, pc_x, pc_y, sid0_141, sid0_143, \
                         sip_53, sip_71, sid1_141, sid1_143, skp_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_6 * sip_71[k]
                   + f_3 * pc_x[k] * skp_71[k];

        t_141[k] = pb_x[k] * sid0_141[k]
                   - f_4 * pc_x[k] * sid1_141[k];

        t_142[k] = f_9 * sip_53[k]
                   + f_3 * pc_y[k] * skp_71[k];

        t_143[k] = pb_x[k] * sid0_143[k]
                   - f_4 * pc_x[k] * sid1_143[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, pb_x, pc_x, sid0_144, sid0_147, sip_72, \
                         sip_73, sip_74, sid1_144, sid1_147, skp_73, \
                         skp_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = pb_x[k] * sid0_144[k]
                   + f_8 * sip_72[k]
                   - f_4 * pc_x[k] * sid1_144[k];

        t_145[k] = f_6 * sip_73[k]
                   + f_3 * pc_x[k] * skp_73[k];

        t_146[k] = f_6 * sip_74[k]
                   + f_3 * pc_x[k] * skp_74[k];

        t_147[k] = pb_x[k] * sid0_147[k]
                   - f_4 * pc_x[k] * sid1_147[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, t_151, pb_x, pc_x, pc_y, sid0_149, sid0_150, \
                         sip_56, sip_75, sip_76, sid1_149, sid1_150, skp_74, \
                         skp_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_10 * sip_56[k]
                   + f_3 * pc_y[k] * skp_74[k];

        t_149[k] = pb_x[k] * sid0_149[k]
                   - f_4 * pc_x[k] * sid1_149[k];

        t_150[k] = pb_x[k] * sid0_150[k]
                   + f_8 * sip_75[k]
                   - f_4 * pc_x[k] * sid1_150[k];

        t_151[k] = f_6 * sip_76[k]
                   + f_3 * pc_x[k] * skp_76[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, pb_x, pc_x, pc_y, sid0_153, sid0_155, \
                         sip_59, sip_77, sid1_153, sid1_155, skp_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_6 * sip_77[k]
                   + f_3 * pc_x[k] * skp_77[k];

        t_153[k] = pb_x[k] * sid0_153[k]
                   - f_4 * pc_x[k] * sid1_153[k];

        t_154[k] = f_8 * sip_59[k]
                   + f_3 * pc_y[k] * skp_77[k];

        t_155[k] = pb_x[k] * sid0_155[k]
                   - f_4 * pc_x[k] * sid1_155[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, pb_x, pb_y, pc_x, pc_y, sid0_120, \
                         sid0_159, sip_79, sip_80, sid1_120, sid1_159, skp_79, \
                         skp_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = pb_y[k] * sid0_120[k]
                   - f_4 * pc_y[k] * sid1_120[k];

        t_157[k] = f_6 * sip_79[k]
                   + f_3 * pc_x[k] * skp_79[k];

        t_158[k] = f_6 * sip_80[k]
                   + f_3 * pc_x[k] * skp_80[k];

        t_159[k] = pb_x[k] * sid0_159[k]
                   - f_4 * pc_x[k] * sid1_159[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, pb_x, pc_x, pc_y, sid0_161, sid0_162, \
                         sip_62, sip_81, sip_82, sid1_161, sid1_162, skp_80, \
                         skp_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = f_6 * sip_62[k]
                   + f_3 * pc_y[k] * skp_80[k];

        t_161[k] = pb_x[k] * sid0_161[k]
                   - f_4 * pc_x[k] * sid1_161[k];

        t_162[k] = pb_x[k] * sid0_162[k]
                   + f_8 * sip_81[k]
                   - f_4 * pc_x[k] * sid1_162[k];

        t_163[k] = f_6 * sip_82[k]
                   + f_3 * pc_x[k] * skp_82[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, pb_x, pc_x, pc_y, sid0_165, sid0_167, \
                         sip_83, sid1_165, sid1_167, skp_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_6 * sip_83[k]
                   + f_3 * pc_x[k] * skp_83[k];

        t_165[k] = pb_x[k] * sid0_165[k]
                   - f_4 * pc_x[k] * sid1_165[k];

        t_166[k] = f_3 * pc_y[k] * skp_83[k];

        t_167[k] = pb_x[k] * sid0_167[k]
                   - f_4 * pc_x[k] * sid1_167[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, t_172, t_173, pc_x, pc_y, pc_z, sip_64, \
                         sip_65, sks0_28, sks1_28, skp_84, skp_85, \
                         skp_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_1 * sks0_28[k]
                   - f_2 * sks1_28[k]
                   + f_3 * pc_x[k] * skp_84[k];

        t_169[k] = f_3 * pc_x[k] * skp_85[k];

        t_170[k] = f_3 * pc_x[k] * skp_86[k];

        t_171[k] = f_0 * sip_64[k]
                   + f_1 * sks0_28[k]
                   - f_2 * sks1_28[k]
                   + f_3 * pc_y[k] * skp_85[k];

        t_172[k] = f_0 * sip_65[k]
                   + f_3 * pc_y[k] * skp_86[k];

        t_173[k] = f_1 * sks0_28[k]
                   - f_2 * sks1_28[k]
                   + f_3 * pc_z[k] * skp_86[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, t_178, pb_z, pc_x, pc_y, pc_z, sid0_126, \
                         sid0_129, sip_68, sid1_126, sid1_129, skp_88, \
                         skp_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = pb_z[k] * sid0_126[k]
                   - f_4 * pc_z[k] * sid1_126[k];

        t_175[k] = f_3 * pc_x[k] * skp_88[k];

        t_176[k] = f_3 * pc_x[k] * skp_89[k];

        t_177[k] = pb_z[k] * sid0_129[k]
                   - f_4 * pc_z[k] * sid1_129[k];

        t_178[k] = f_5 * sip_68[k]
                   + f_3 * pc_y[k] * skp_89[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, t_182, pc_x, pc_z, sip_65, sks0_29, sks0_30, \
                         sks1_29, sks1_30, skp_89, skp_90, skp_91, \
                         skp_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = f_6 * sip_65[k]
                   + f_1 * sks0_29[k]
                   - f_2 * sks1_29[k]
                   + f_3 * pc_z[k] * skp_89[k];

        t_180[k] = f_1 * sks0_30[k]
                   - f_2 * sks1_30[k]
                   + f_3 * pc_x[k] * skp_90[k];

        t_181[k] = f_3 * pc_x[k] * skp_91[k];

        t_182[k] = f_3 * pc_x[k] * skp_92[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, pc_y, pc_z, sip_68, sip_70, sip_71, sks0_30, \
                         sks1_30, skp_91, skp_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_7 * sip_70[k]
                   + f_1 * sks0_30[k]
                   - f_2 * sks1_30[k]
                   + f_3 * pc_y[k] * skp_91[k];

        t_184[k] = f_7 * sip_71[k]
                   + f_3 * pc_y[k] * skp_92[k];

        t_185[k] = f_8 * sip_68[k]
                   + f_1 * sks0_30[k]
                   - f_2 * sks1_30[k]
                   + f_3 * pc_z[k] * skp_92[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, t_189, t_190, pc_x, pc_y, sip_73, sip_74, \
                         sks0_31, sks1_31, skp_93, skp_94, skp_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_1 * sks0_31[k]
                   - f_2 * sks1_31[k]
                   + f_3 * pc_x[k] * skp_93[k];

        t_187[k] = f_3 * pc_x[k] * skp_94[k];

        t_188[k] = f_3 * pc_x[k] * skp_95[k];

        t_189[k] = f_9 * sip_73[k]
                   + f_1 * sks0_31[k]
                   - f_2 * sks1_31[k]
                   + f_3 * pc_y[k] * skp_94[k];

        t_190[k] = f_9 * sip_74[k]
                   + f_3 * pc_y[k] * skp_95[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, t_194, pc_x, pc_z, sip_71, sks0_31, sks0_32, \
                         sks1_31, sks1_32, skp_95, skp_96, skp_97, \
                         skp_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = f_10 * sip_71[k]
                   + f_1 * sks0_31[k]
                   - f_2 * sks1_31[k]
                   + f_3 * pc_z[k] * skp_95[k];

        t_192[k] = f_1 * sks0_32[k]
                   - f_2 * sks1_32[k]
                   + f_3 * pc_x[k] * skp_96[k];

        t_193[k] = f_3 * pc_x[k] * skp_97[k];

        t_194[k] = f_3 * pc_x[k] * skp_98[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, pc_y, pc_z, sip_74, sip_76, sip_77, sks0_32, \
                         sks1_32, skp_97, skp_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = f_10 * sip_76[k]
                   + f_1 * sks0_32[k]
                   - f_2 * sks1_32[k]
                   + f_3 * pc_y[k] * skp_97[k];

        t_196[k] = f_10 * sip_77[k]
                   + f_3 * pc_y[k] * skp_98[k];

        t_197[k] = f_9 * sip_74[k]
                   + f_1 * sks0_32[k]
                   - f_2 * sks1_32[k]
                   + f_3 * pc_z[k] * skp_98[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, t_202, pc_x, pc_y, sip_79, sip_80, \
                         sks0_33, sks1_33, skp_99, skp_100, skp_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_1 * sks0_33[k]
                   - f_2 * sks1_33[k]
                   + f_3 * pc_x[k] * skp_99[k];

        t_199[k] = f_3 * pc_x[k] * skp_100[k];

        t_200[k] = f_3 * pc_x[k] * skp_101[k];

        t_201[k] = f_8 * sip_79[k]
                   + f_1 * sks0_33[k]
                   - f_2 * sks1_33[k]
                   + f_3 * pc_y[k] * skp_100[k];

        t_202[k] = f_8 * sip_80[k]
                   + f_3 * pc_y[k] * skp_101[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, pb_y, pc_x, pc_y, pc_z, sid0_162, sip_77, \
                         sid1_162, sks0_33, sks1_33, skp_101, skp_103, \
                         skp_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_7 * sip_77[k]
                   + f_1 * sks0_33[k]
                   - f_2 * sks1_33[k]
                   + f_3 * pc_z[k] * skp_101[k];

        t_204[k] = pb_y[k] * sid0_162[k]
                   - f_4 * pc_y[k] * sid1_162[k];

        t_205[k] = f_3 * pc_x[k] * skp_103[k];

        t_206[k] = f_3 * pc_x[k] * skp_104[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, pb_y, pc_y, sid0_165, sid0_167, sip_82, sip_83, \
                         sid1_165, sid1_167, skp_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = pb_y[k] * sid0_165[k]
                   + f_8 * sip_82[k]
                   - f_4 * pc_y[k] * sid1_165[k];

        t_208[k] = f_6 * sip_83[k]
                   + f_3 * pc_y[k] * skp_104[k];

        t_209[k] = pb_y[k] * sid0_167[k]
                   - f_4 * pc_y[k] * sid1_167[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, t_215, pc_x, pc_y, pc_z, sip_83, \
                         sks0_35, sks1_35, skp_105, skp_106, skp_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_1 * sks0_35[k]
                   - f_2 * sks1_35[k]
                   + f_3 * pc_x[k] * skp_105[k];

        t_211[k] = f_3 * pc_x[k] * skp_106[k];

        t_212[k] = f_3 * pc_x[k] * skp_107[k];

        t_213[k] = f_1 * sks0_35[k]
                   - f_2 * sks1_35[k]
                   + f_3 * pc_y[k] * skp_106[k];

        t_214[k] = f_3 * pc_y[k] * skp_107[k];

        t_215[k] = f_0 * sip_83[k]
                   + f_1 * sks0_35[k]
                   - f_2 * sks1_35[k]
                   + f_3 * pc_z[k] * skp_107[k];
    }
}

auto
compute_prim_skd_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t sid0, const size_t sip,
                                                   const size_t sid1, const size_t sks0,
                                                   const size_t sks1, const size_t skp,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_skd_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, sid0, sip,
                                                              sid1, sks0, sks1, skp, ncols,
                                                              gamma, p, q);

    compute_prim_skd_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, sid0, sip,
                                                              sid1, sks0, sks1, skp, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
