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


#include "SimdThreeCenterElectronRepulsionVrrRecSOD.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_sod_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t snd0,
                                                          const size_t snp, const size_t snd1,
                                                          const size_t sos0, const size_t sos1,
                                                          const size_t sop, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.5 / q;
    const auto f_1 = 0.5 / gamma;
    const auto f_2 = 0.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 5.0 / q;
    const auto f_6 = 0.5 / q;
    const auto f_7 = 4.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 4.0 / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 3.5 / q;
    const auto f_12 = 2.0 / q;
    const auto f_13 = 3.0 / q;
    const auto f_14 = 2.5 / q;

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

    const auto *snd0_0 = buffer.data(snd0 + 0);
    const auto *snd0_3 = buffer.data(snd0 + 3);
    const auto *snd0_5 = buffer.data(snd0 + 5);
    const auto *snd0_9 = buffer.data(snd0 + 9);
    const auto *snd0_12 = buffer.data(snd0 + 12);
    const auto *snd0_17 = buffer.data(snd0 + 17);
    const auto *snd0_18 = buffer.data(snd0 + 18);
    const auto *snd0_21 = buffer.data(snd0 + 21);
    const auto *snd0_30 = buffer.data(snd0 + 30);
    const auto *snd0_35 = buffer.data(snd0 + 35);
    const auto *snd0_36 = buffer.data(snd0 + 36);
    const auto *snd0_39 = buffer.data(snd0 + 39);
    const auto *snd0_54 = buffer.data(snd0 + 54);
    const auto *snd0_59 = buffer.data(snd0 + 59);
    const auto *snd0_60 = buffer.data(snd0 + 60);
    const auto *snd0_63 = buffer.data(snd0 + 63);
    const auto *snd0_84 = buffer.data(snd0 + 84);
    const auto *snd0_89 = buffer.data(snd0 + 89);

    const auto *snp_0 = buffer.data(snp + 0);
    const auto *snp_1 = buffer.data(snp + 1);
    const auto *snp_2 = buffer.data(snp + 2);
    const auto *snp_4 = buffer.data(snp + 4);
    const auto *snp_5 = buffer.data(snp + 5);
    const auto *snp_7 = buffer.data(snp + 7);
    const auto *snp_8 = buffer.data(snp + 8);
    const auto *snp_9 = buffer.data(snp + 9);
    const auto *snp_10 = buffer.data(snp + 10);
    const auto *snp_11 = buffer.data(snp + 11);
    const auto *snp_13 = buffer.data(snp + 13);
    const auto *snp_14 = buffer.data(snp + 14);
    const auto *snp_15 = buffer.data(snp + 15);
    const auto *snp_16 = buffer.data(snp + 16);
    const auto *snp_17 = buffer.data(snp + 17);
    const auto *snp_18 = buffer.data(snp + 18);
    const auto *snp_19 = buffer.data(snp + 19);
    const auto *snp_20 = buffer.data(snp + 20);
    const auto *snp_22 = buffer.data(snp + 22);
    const auto *snp_23 = buffer.data(snp + 23);
    const auto *snp_25 = buffer.data(snp + 25);
    const auto *snp_26 = buffer.data(snp + 26);
    const auto *snp_27 = buffer.data(snp + 27);
    const auto *snp_28 = buffer.data(snp + 28);
    const auto *snp_29 = buffer.data(snp + 29);
    const auto *snp_30 = buffer.data(snp + 30);
    const auto *snp_31 = buffer.data(snp + 31);
    const auto *snp_32 = buffer.data(snp + 32);
    const auto *snp_34 = buffer.data(snp + 34);
    const auto *snp_35 = buffer.data(snp + 35);
    const auto *snp_36 = buffer.data(snp + 36);
    const auto *snp_37 = buffer.data(snp + 37);
    const auto *snp_38 = buffer.data(snp + 38);
    const auto *snp_40 = buffer.data(snp + 40);
    const auto *snp_41 = buffer.data(snp + 41);
    const auto *snp_42 = buffer.data(snp + 42);
    const auto *snp_43 = buffer.data(snp + 43);
    const auto *snp_44 = buffer.data(snp + 44);
    const auto *snp_45 = buffer.data(snp + 45);
    const auto *snp_46 = buffer.data(snp + 46);
    const auto *snp_47 = buffer.data(snp + 47);
    const auto *snp_49 = buffer.data(snp + 49);
    const auto *snp_50 = buffer.data(snp + 50);
    const auto *snp_51 = buffer.data(snp + 51);
    const auto *snp_52 = buffer.data(snp + 52);
    const auto *snp_53 = buffer.data(snp + 53);
    const auto *snp_54 = buffer.data(snp + 54);
    const auto *snp_55 = buffer.data(snp + 55);
    const auto *snp_56 = buffer.data(snp + 56);
    const auto *snp_58 = buffer.data(snp + 58);
    const auto *snp_59 = buffer.data(snp + 59);

    const auto *snd1_0 = buffer.data(snd1 + 0);
    const auto *snd1_3 = buffer.data(snd1 + 3);
    const auto *snd1_5 = buffer.data(snd1 + 5);
    const auto *snd1_9 = buffer.data(snd1 + 9);
    const auto *snd1_12 = buffer.data(snd1 + 12);
    const auto *snd1_17 = buffer.data(snd1 + 17);
    const auto *snd1_18 = buffer.data(snd1 + 18);
    const auto *snd1_21 = buffer.data(snd1 + 21);
    const auto *snd1_30 = buffer.data(snd1 + 30);
    const auto *snd1_35 = buffer.data(snd1 + 35);
    const auto *snd1_36 = buffer.data(snd1 + 36);
    const auto *snd1_39 = buffer.data(snd1 + 39);
    const auto *snd1_54 = buffer.data(snd1 + 54);
    const auto *snd1_59 = buffer.data(snd1 + 59);
    const auto *snd1_60 = buffer.data(snd1 + 60);
    const auto *snd1_63 = buffer.data(snd1 + 63);
    const auto *snd1_84 = buffer.data(snd1 + 84);
    const auto *snd1_89 = buffer.data(snd1 + 89);

    const auto *sos0_0 = buffer.data(sos0 + 0);
    const auto *sos0_1 = buffer.data(sos0 + 1);
    const auto *sos0_2 = buffer.data(sos0 + 2);
    const auto *sos0_3 = buffer.data(sos0 + 3);
    const auto *sos0_5 = buffer.data(sos0 + 5);
    const auto *sos0_6 = buffer.data(sos0 + 6);
    const auto *sos0_7 = buffer.data(sos0 + 7);
    const auto *sos0_8 = buffer.data(sos0 + 8);
    const auto *sos0_9 = buffer.data(sos0 + 9);
    const auto *sos0_10 = buffer.data(sos0 + 10);
    const auto *sos0_11 = buffer.data(sos0 + 11);
    const auto *sos0_12 = buffer.data(sos0 + 12);
    const auto *sos0_13 = buffer.data(sos0 + 13);
    const auto *sos0_14 = buffer.data(sos0 + 14);
    const auto *sos0_15 = buffer.data(sos0 + 15);
    const auto *sos0_16 = buffer.data(sos0 + 16);
    const auto *sos0_17 = buffer.data(sos0 + 17);
    const auto *sos0_18 = buffer.data(sos0 + 18);
    const auto *sos0_19 = buffer.data(sos0 + 19);

    const auto *sos1_0 = buffer.data(sos1 + 0);
    const auto *sos1_1 = buffer.data(sos1 + 1);
    const auto *sos1_2 = buffer.data(sos1 + 2);
    const auto *sos1_3 = buffer.data(sos1 + 3);
    const auto *sos1_5 = buffer.data(sos1 + 5);
    const auto *sos1_6 = buffer.data(sos1 + 6);
    const auto *sos1_7 = buffer.data(sos1 + 7);
    const auto *sos1_8 = buffer.data(sos1 + 8);
    const auto *sos1_9 = buffer.data(sos1 + 9);
    const auto *sos1_10 = buffer.data(sos1 + 10);
    const auto *sos1_11 = buffer.data(sos1 + 11);
    const auto *sos1_12 = buffer.data(sos1 + 12);
    const auto *sos1_13 = buffer.data(sos1 + 13);
    const auto *sos1_14 = buffer.data(sos1 + 14);
    const auto *sos1_15 = buffer.data(sos1 + 15);
    const auto *sos1_16 = buffer.data(sos1 + 16);
    const auto *sos1_17 = buffer.data(sos1 + 17);
    const auto *sos1_18 = buffer.data(sos1 + 18);
    const auto *sos1_19 = buffer.data(sos1 + 19);

    const auto *sop_0 = buffer.data(sop + 0);
    const auto *sop_1 = buffer.data(sop + 1);
    const auto *sop_2 = buffer.data(sop + 2);
    const auto *sop_4 = buffer.data(sop + 4);
    const auto *sop_5 = buffer.data(sop + 5);
    const auto *sop_7 = buffer.data(sop + 7);
    const auto *sop_8 = buffer.data(sop + 8);
    const auto *sop_9 = buffer.data(sop + 9);
    const auto *sop_10 = buffer.data(sop + 10);
    const auto *sop_11 = buffer.data(sop + 11);
    const auto *sop_13 = buffer.data(sop + 13);
    const auto *sop_14 = buffer.data(sop + 14);
    const auto *sop_15 = buffer.data(sop + 15);
    const auto *sop_16 = buffer.data(sop + 16);
    const auto *sop_17 = buffer.data(sop + 17);
    const auto *sop_18 = buffer.data(sop + 18);
    const auto *sop_19 = buffer.data(sop + 19);
    const auto *sop_20 = buffer.data(sop + 20);
    const auto *sop_22 = buffer.data(sop + 22);
    const auto *sop_23 = buffer.data(sop + 23);
    const auto *sop_25 = buffer.data(sop + 25);
    const auto *sop_26 = buffer.data(sop + 26);
    const auto *sop_27 = buffer.data(sop + 27);
    const auto *sop_28 = buffer.data(sop + 28);
    const auto *sop_29 = buffer.data(sop + 29);
    const auto *sop_30 = buffer.data(sop + 30);
    const auto *sop_31 = buffer.data(sop + 31);
    const auto *sop_32 = buffer.data(sop + 32);
    const auto *sop_34 = buffer.data(sop + 34);
    const auto *sop_35 = buffer.data(sop + 35);
    const auto *sop_36 = buffer.data(sop + 36);
    const auto *sop_37 = buffer.data(sop + 37);
    const auto *sop_38 = buffer.data(sop + 38);
    const auto *sop_40 = buffer.data(sop + 40);
    const auto *sop_41 = buffer.data(sop + 41);
    const auto *sop_42 = buffer.data(sop + 42);
    const auto *sop_43 = buffer.data(sop + 43);
    const auto *sop_44 = buffer.data(sop + 44);
    const auto *sop_45 = buffer.data(sop + 45);
    const auto *sop_46 = buffer.data(sop + 46);
    const auto *sop_47 = buffer.data(sop + 47);
    const auto *sop_49 = buffer.data(sop + 49);
    const auto *sop_50 = buffer.data(sop + 50);
    const auto *sop_51 = buffer.data(sop + 51);
    const auto *sop_52 = buffer.data(sop + 52);
    const auto *sop_53 = buffer.data(sop + 53);
    const auto *sop_54 = buffer.data(sop + 54);
    const auto *sop_55 = buffer.data(sop + 55);
    const auto *sop_56 = buffer.data(sop + 56);
    const auto *sop_58 = buffer.data(sop + 58);
    const auto *sop_59 = buffer.data(sop + 59);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pc_x, pc_y, snp_0, snp_1, snp_2, sos0_0, \
                         sos1_0, sop_0, sop_1, sop_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * snp_0[k]
                 + f_1 * sos0_0[k]
                 - f_2 * sos1_0[k]
                 + f_3 * pc_x[k] * sop_0[k];

        t_1[k] = f_0 * snp_1[k]
                 + f_3 * pc_x[k] * sop_1[k];

        t_2[k] = f_0 * snp_2[k]
                 + f_3 * pc_x[k] * sop_2[k];

        t_3[k] = f_1 * sos0_0[k]
                 - f_2 * sos1_0[k]
                 + f_3 * pc_y[k] * sop_1[k];

        t_4[k] = f_3 * pc_y[k] * sop_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pb_y, pc_x, pc_y, pc_z, snd0_0, snp_4, snd1_0, sos0_0, \
                         sos1_0, sop_2, sop_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * sos0_0[k]
                 - f_2 * sos1_0[k]
                 + f_3 * pc_z[k] * sop_2[k];

        t_6[k] = pb_y[k] * snd0_0[k]
                 - f_4 * pc_y[k] * snd1_0[k];

        t_7[k] = f_5 * snp_4[k]
                 + f_3 * pc_x[k] * sop_4[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pb_y, pc_x, pc_y, snd0_5, snp_1, snp_2, snp_5, \
                         snd1_5, sos0_1, sos1_1, sop_4, sop_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_5 * snp_5[k]
                 + f_3 * pc_x[k] * sop_5[k];

        t_9[k] = f_6 * snp_1[k]
                 + f_1 * sos0_1[k]
                 - f_2 * sos1_1[k]
                 + f_3 * pc_y[k] * sop_4[k];

        t_10[k] = f_6 * snp_2[k]
                  + f_3 * pc_y[k] * sop_5[k];

        t_11[k] = pb_y[k] * snd0_5[k]
                  - f_4 * pc_y[k] * snd1_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pb_z, pc_x, pc_z, snd0_0, snd0_3, snp_7, \
                         snp_8, snd1_0, snd1_3, sop_7, sop_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pb_z[k] * snd0_0[k]
                  - f_4 * pc_z[k] * snd1_0[k];

        t_13[k] = f_5 * snp_7[k]
                  + f_3 * pc_x[k] * sop_7[k];

        t_14[k] = f_5 * snp_8[k]
                  + f_3 * pc_x[k] * sop_8[k];

        t_15[k] = pb_z[k] * snd0_3[k]
                  - f_4 * pc_z[k] * snd1_3[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pc_x, pc_y, pc_z, snp_2, snp_9, sos0_2, sos0_3, \
                         sos1_2, sos1_3, sop_8, sop_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * pc_y[k] * sop_8[k];

        t_17[k] = f_6 * snp_2[k]
                  + f_1 * sos0_2[k]
                  - f_2 * sos1_2[k]
                  + f_3 * pc_z[k] * sop_8[k];

        t_18[k] = f_7 * snp_9[k]
                  + f_1 * sos0_3[k]
                  - f_2 * sos1_3[k]
                  + f_3 * pc_x[k] * sop_9[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pc_x, pc_y, pc_z, snp_4, snp_5, snp_10, \
                         snp_11, sos0_3, sos1_3, sop_10, sop_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_7 * snp_10[k]
                  + f_3 * pc_x[k] * sop_10[k];

        t_20[k] = f_7 * snp_11[k]
                  + f_3 * pc_x[k] * sop_11[k];

        t_21[k] = f_8 * snp_4[k]
                  + f_1 * sos0_3[k]
                  - f_2 * sos1_3[k]
                  + f_3 * pc_y[k] * sop_10[k];

        t_22[k] = f_8 * snp_5[k]
                  + f_3 * pc_y[k] * sop_11[k];

        t_23[k] = f_1 * sos0_3[k]
                  - f_2 * sos1_3[k]
                  + f_3 * pc_z[k] * sop_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pb_y, pc_x, pc_y, snd0_12, snp_13, snp_14, snd1_12, \
                         sop_13, sop_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pb_y[k] * snd0_12[k]
                  - f_4 * pc_y[k] * snd1_12[k];

        t_25[k] = f_7 * snp_13[k]
                  + f_3 * pc_x[k] * sop_13[k];

        t_26[k] = f_7 * snp_14[k]
                  + f_3 * pc_x[k] * sop_14[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pb_y, pb_z, pc_y, pc_z, snd0_9, snd0_17, snp_8, \
                         snd1_9, snd1_17, sop_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = pb_z[k] * snd0_9[k]
                  - f_4 * pc_z[k] * snd1_9[k];

        t_28[k] = f_6 * snp_8[k]
                  + f_3 * pc_y[k] * sop_14[k];

        t_29[k] = pb_y[k] * snd0_17[k]
                  - f_4 * pc_y[k] * snd1_17[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pc_x, pc_y, snp_15, snp_16, snp_17, \
                         sos0_5, sos1_5, sop_15, sop_16, sop_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_7 * snp_15[k]
                  + f_1 * sos0_5[k]
                  - f_2 * sos1_5[k]
                  + f_3 * pc_x[k] * sop_15[k];

        t_31[k] = f_7 * snp_16[k]
                  + f_3 * pc_x[k] * sop_16[k];

        t_32[k] = f_7 * snp_17[k]
                  + f_3 * pc_x[k] * sop_17[k];

        t_33[k] = f_1 * sos0_5[k]
                  - f_2 * sos1_5[k]
                  + f_3 * pc_y[k] * sop_16[k];

        t_34[k] = f_3 * pc_y[k] * sop_17[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pc_x, pc_z, snp_8, snp_18, snp_19, sos0_5, sos0_6, \
                         sos1_5, sos1_6, sop_17, sop_18, sop_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_8 * snp_8[k]
                  + f_1 * sos0_5[k]
                  - f_2 * sos1_5[k]
                  + f_3 * pc_z[k] * sop_17[k];

        t_36[k] = f_9 * snp_18[k]
                  + f_1 * sos0_6[k]
                  - f_2 * sos1_6[k]
                  + f_3 * pc_x[k] * sop_18[k];

        t_37[k] = f_9 * snp_19[k]
                  + f_3 * pc_x[k] * sop_19[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pc_x, pc_y, pc_z, snp_10, snp_11, snp_20, \
                         sos0_6, sos1_6, sop_19, sop_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_9 * snp_20[k]
                  + f_3 * pc_x[k] * sop_20[k];

        t_39[k] = f_10 * snp_10[k]
                  + f_1 * sos0_6[k]
                  - f_2 * sos1_6[k]
                  + f_3 * pc_y[k] * sop_19[k];

        t_40[k] = f_10 * snp_11[k]
                  + f_3 * pc_y[k] * sop_20[k];

        t_41[k] = f_1 * sos0_6[k]
                  - f_2 * sos1_6[k]
                  + f_3 * pc_z[k] * sop_20[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pb_z, pc_x, pc_z, snd0_18, snd0_21, snp_22, \
                         snp_23, snd1_18, snd1_21, sop_22, sop_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_z[k] * snd0_18[k]
                  - f_4 * pc_z[k] * snd1_18[k];

        t_43[k] = f_9 * snp_22[k]
                  + f_3 * pc_x[k] * sop_22[k];

        t_44[k] = f_9 * snp_23[k]
                  + f_3 * pc_x[k] * sop_23[k];

        t_45[k] = pb_z[k] * snd0_21[k]
                  - f_4 * pc_z[k] * snd1_21[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, pb_y, pc_y, pc_z, snd0_30, snp_11, snp_14, snd1_30, \
                         sos0_7, sos1_7, sop_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_8 * snp_14[k]
                  + f_3 * pc_y[k] * sop_23[k];

        t_47[k] = f_6 * snp_11[k]
                  + f_1 * sos0_7[k]
                  - f_2 * sos1_7[k]
                  + f_3 * pc_z[k] * sop_23[k];

        t_48[k] = pb_y[k] * snd0_30[k]
                  - f_4 * pc_y[k] * snd1_30[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pc_x, pc_y, snp_16, snp_17, snp_25, snp_26, \
                         sos0_8, sos1_8, sop_25, sop_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_9 * snp_25[k]
                  + f_3 * pc_x[k] * sop_25[k];

        t_50[k] = f_9 * snp_26[k]
                  + f_3 * pc_x[k] * sop_26[k];

        t_51[k] = f_6 * snp_16[k]
                  + f_1 * sos0_8[k]
                  - f_2 * sos1_8[k]
                  + f_3 * pc_y[k] * sop_25[k];

        t_52[k] = f_6 * snp_17[k]
                  + f_3 * pc_y[k] * sop_26[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pb_y, pc_x, pc_y, snd0_35, snp_27, snp_28, snd1_35, \
                         sos0_9, sos1_9, sop_27, sop_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = pb_y[k] * snd0_35[k]
                  - f_4 * pc_y[k] * snd1_35[k];

        t_54[k] = f_9 * snp_27[k]
                  + f_1 * sos0_9[k]
                  - f_2 * sos1_9[k]
                  + f_3 * pc_x[k] * sop_27[k];

        t_55[k] = f_9 * snp_28[k]
                  + f_3 * pc_x[k] * sop_28[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pc_x, pc_y, pc_z, snp_17, snp_29, sos0_9, \
                         sos1_9, sop_28, sop_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_9 * snp_29[k]
                  + f_3 * pc_x[k] * sop_29[k];

        t_57[k] = f_1 * sos0_9[k]
                  - f_2 * sos1_9[k]
                  + f_3 * pc_y[k] * sop_28[k];

        t_58[k] = f_3 * pc_y[k] * sop_29[k];

        t_59[k] = f_10 * snp_17[k]
                  + f_1 * sos0_9[k]
                  - f_2 * sos1_9[k]
                  + f_3 * pc_z[k] * sop_29[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pc_x, pc_y, snp_19, snp_30, snp_31, snp_32, \
                         sos0_10, sos1_10, sop_30, sop_31, sop_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_11 * snp_30[k]
                  + f_1 * sos0_10[k]
                  - f_2 * sos1_10[k]
                  + f_3 * pc_x[k] * sop_30[k];

        t_61[k] = f_11 * snp_31[k]
                  + f_3 * pc_x[k] * sop_31[k];

        t_62[k] = f_11 * snp_32[k]
                  + f_3 * pc_x[k] * sop_32[k];

        t_63[k] = f_12 * snp_19[k]
                  + f_1 * sos0_10[k]
                  - f_2 * sos1_10[k]
                  + f_3 * pc_y[k] * sop_31[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pb_z, pc_x, pc_y, pc_z, snd0_36, snp_20, \
                         snp_34, snd1_36, sos0_10, sos1_10, sop_32, \
                         sop_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_12 * snp_20[k]
                  + f_3 * pc_y[k] * sop_32[k];

        t_65[k] = f_1 * sos0_10[k]
                  - f_2 * sos1_10[k]
                  + f_3 * pc_z[k] * sop_32[k];

        t_66[k] = pb_z[k] * snd0_36[k]
                  - f_4 * pc_z[k] * snd1_36[k];

        t_67[k] = f_11 * snp_34[k]
                  + f_3 * pc_x[k] * sop_34[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pb_z, pc_x, pc_y, pc_z, snd0_39, snp_20, \
                         snp_23, snp_35, snd1_39, sos0_11, sos1_11, \
                         sop_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_11 * snp_35[k]
                  + f_3 * pc_x[k] * sop_35[k];

        t_69[k] = pb_z[k] * snd0_39[k]
                  - f_4 * pc_z[k] * snd1_39[k];

        t_70[k] = f_10 * snp_23[k]
                  + f_3 * pc_y[k] * sop_35[k];

        t_71[k] = f_6 * snp_20[k]
                  + f_1 * sos0_11[k]
                  - f_2 * sos1_11[k]
                  + f_3 * pc_z[k] * sop_35[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pc_x, pc_y, snp_25, snp_36, snp_37, snp_38, \
                         sos0_12, sos1_12, sop_36, sop_37, sop_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_11 * snp_36[k]
                  + f_1 * sos0_12[k]
                  - f_2 * sos1_12[k]
                  + f_3 * pc_x[k] * sop_36[k];

        t_73[k] = f_11 * snp_37[k]
                  + f_3 * pc_x[k] * sop_37[k];

        t_74[k] = f_11 * snp_38[k]
                  + f_3 * pc_x[k] * sop_38[k];

        t_75[k] = f_8 * snp_25[k]
                  + f_1 * sos0_12[k]
                  - f_2 * sos1_12[k]
                  + f_3 * pc_y[k] * sop_37[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, pb_y, pc_y, pc_z, snd0_54, snp_23, snp_26, snd1_54, \
                         sos0_12, sos1_12, sop_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_8 * snp_26[k]
                  + f_3 * pc_y[k] * sop_38[k];

        t_77[k] = f_8 * snp_23[k]
                  + f_1 * sos0_12[k]
                  - f_2 * sos1_12[k]
                  + f_3 * pc_z[k] * sop_38[k];

        t_78[k] = pb_y[k] * snd0_54[k]
                  - f_4 * pc_y[k] * snd1_54[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pc_x, pc_y, snp_28, snp_29, snp_40, snp_41, \
                         sos0_13, sos1_13, sop_40, sop_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_11 * snp_40[k]
                  + f_3 * pc_x[k] * sop_40[k];

        t_80[k] = f_11 * snp_41[k]
                  + f_3 * pc_x[k] * sop_41[k];

        t_81[k] = f_6 * snp_28[k]
                  + f_1 * sos0_13[k]
                  - f_2 * sos1_13[k]
                  + f_3 * pc_y[k] * sop_40[k];

        t_82[k] = f_6 * snp_29[k]
                  + f_3 * pc_y[k] * sop_41[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, pb_y, pc_x, pc_y, snd0_59, snp_42, snp_43, snd1_59, \
                         sos0_14, sos1_14, sop_42, sop_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = pb_y[k] * snd0_59[k]
                  - f_4 * pc_y[k] * snd1_59[k];

        t_84[k] = f_11 * snp_42[k]
                  + f_1 * sos0_14[k]
                  - f_2 * sos1_14[k]
                  + f_3 * pc_x[k] * sop_42[k];

        t_85[k] = f_11 * snp_43[k]
                  + f_3 * pc_x[k] * sop_43[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pc_x, pc_y, pc_z, snp_29, snp_44, sos0_14, \
                         sos1_14, sop_43, sop_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_11 * snp_44[k]
                  + f_3 * pc_x[k] * sop_44[k];

        t_87[k] = f_1 * sos0_14[k]
                  - f_2 * sos1_14[k]
                  + f_3 * pc_y[k] * sop_43[k];

        t_88[k] = f_3 * pc_y[k] * sop_44[k];

        t_89[k] = f_12 * snp_29[k]
                  + f_1 * sos0_14[k]
                  - f_2 * sos1_14[k]
                  + f_3 * pc_z[k] * sop_44[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, pc_x, pc_y, snp_31, snp_45, snp_46, snp_47, \
                         sos0_15, sos1_15, sop_45, sop_46, sop_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_13 * snp_45[k]
                  + f_1 * sos0_15[k]
                  - f_2 * sos1_15[k]
                  + f_3 * pc_x[k] * sop_45[k];

        t_91[k] = f_13 * snp_46[k]
                  + f_3 * pc_x[k] * sop_46[k];

        t_92[k] = f_13 * snp_47[k]
                  + f_3 * pc_x[k] * sop_47[k];

        t_93[k] = f_14 * snp_31[k]
                  + f_1 * sos0_15[k]
                  - f_2 * sos1_15[k]
                  + f_3 * pc_y[k] * sop_46[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, pb_z, pc_x, pc_y, pc_z, snd0_60, snp_32, \
                         snp_49, snd1_60, sos0_15, sos1_15, sop_47, \
                         sop_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_14 * snp_32[k]
                  + f_3 * pc_y[k] * sop_47[k];

        t_95[k] = f_1 * sos0_15[k]
                  - f_2 * sos1_15[k]
                  + f_3 * pc_z[k] * sop_47[k];

        t_96[k] = pb_z[k] * snd0_60[k]
                  - f_4 * pc_z[k] * snd1_60[k];

        t_97[k] = f_13 * snp_49[k]
                  + f_3 * pc_x[k] * sop_49[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, pb_z, pc_x, pc_y, pc_z, snd0_63, snp_32, \
                         snp_35, snp_50, snd1_63, sos0_16, sos1_16, \
                         sop_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_13 * snp_50[k]
                  + f_3 * pc_x[k] * sop_50[k];

        t_99[k] = pb_z[k] * snd0_63[k]
                  - f_4 * pc_z[k] * snd1_63[k];

        t_100[k] = f_12 * snp_35[k]
                   + f_3 * pc_y[k] * sop_50[k];

        t_101[k] = f_6 * snp_32[k]
                   + f_1 * sos0_16[k]
                   - f_2 * sos1_16[k]
                   + f_3 * pc_z[k] * sop_50[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, pc_x, pc_y, snp_37, snp_51, snp_52, \
                         snp_53, sos0_17, sos1_17, sop_51, sop_52, \
                         sop_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_13 * snp_51[k]
                   + f_1 * sos0_17[k]
                   - f_2 * sos1_17[k]
                   + f_3 * pc_x[k] * sop_51[k];

        t_103[k] = f_13 * snp_52[k]
                   + f_3 * pc_x[k] * sop_52[k];

        t_104[k] = f_13 * snp_53[k]
                   + f_3 * pc_x[k] * sop_53[k];

        t_105[k] = f_10 * snp_37[k]
                   + f_1 * sos0_17[k]
                   - f_2 * sos1_17[k]
                   + f_3 * pc_y[k] * sop_52[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, pc_x, pc_y, pc_z, snp_35, snp_38, snp_54, \
                         sos0_17, sos0_18, sos1_17, sos1_18, sop_53, \
                         sop_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_10 * snp_38[k]
                   + f_3 * pc_y[k] * sop_53[k];

        t_107[k] = f_8 * snp_35[k]
                   + f_1 * sos0_17[k]
                   - f_2 * sos1_17[k]
                   + f_3 * pc_z[k] * sop_53[k];

        t_108[k] = f_13 * snp_54[k]
                   + f_1 * sos0_18[k]
                   - f_2 * sos1_18[k]
                   + f_3 * pc_x[k] * sop_54[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pc_x, pc_y, snp_40, snp_41, snp_55, \
                         snp_56, sos0_18, sos1_18, sop_55, sop_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_13 * snp_55[k]
                   + f_3 * pc_x[k] * sop_55[k];

        t_110[k] = f_13 * snp_56[k]
                   + f_3 * pc_x[k] * sop_56[k];

        t_111[k] = f_8 * snp_40[k]
                   + f_1 * sos0_18[k]
                   - f_2 * sos1_18[k]
                   + f_3 * pc_y[k] * sop_55[k];

        t_112[k] = f_8 * snp_41[k]
                   + f_3 * pc_y[k] * sop_56[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, pb_y, pc_x, pc_y, pc_z, snd0_84, snp_38, snp_58, \
                         snd1_84, sos0_18, sos1_18, sop_56, sop_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_10 * snp_38[k]
                   + f_1 * sos0_18[k]
                   - f_2 * sos1_18[k]
                   + f_3 * pc_z[k] * sop_56[k];

        t_114[k] = pb_y[k] * snd0_84[k]
                   - f_4 * pc_y[k] * snd1_84[k];

        t_115[k] = f_13 * snp_58[k]
                   + f_3 * pc_x[k] * sop_58[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, pb_y, pc_x, pc_y, snd0_89, snp_43, \
                         snp_44, snp_59, snd1_89, sos0_19, sos1_19, sop_58, \
                         sop_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_13 * snp_59[k]
                   + f_3 * pc_x[k] * sop_59[k];

        t_117[k] = f_6 * snp_43[k]
                   + f_1 * sos0_19[k]
                   - f_2 * sos1_19[k]
                   + f_3 * pc_y[k] * sop_58[k];

        t_118[k] = f_6 * snp_44[k]
                   + f_3 * pc_y[k] * sop_59[k];

        t_119[k] = pb_y[k] * snd0_89[k]
                   - f_4 * pc_y[k] * snd1_89[k];
    }
}

static auto
compute_prim_sod_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t snd0,
                                                          const size_t snp, const size_t snd1,
                                                          const size_t sos0, const size_t sos1,
                                                          const size_t sop, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 0.5 / gamma;
    const auto f_2 = 0.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_6 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 4.0 / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 3.5 / q;
    const auto f_12 = 2.0 / q;
    const auto f_13 = 3.0 / q;
    const auto f_14 = 2.5 / q;

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
    auto *t_233 = buffer.data(target + 233);
    auto *t_234 = buffer.data(target + 234);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *snd0_90 = buffer.data(snd0 + 90);
    const auto *snd0_93 = buffer.data(snd0 + 93);
    const auto *snd0_120 = buffer.data(snd0 + 120);
    const auto *snd0_125 = buffer.data(snd0 + 125);
    const auto *snd0_126 = buffer.data(snd0 + 126);
    const auto *snd0_129 = buffer.data(snd0 + 129);
    const auto *snd0_162 = buffer.data(snd0 + 162);
    const auto *snd0_167 = buffer.data(snd0 + 167);
    const auto *snd0_168 = buffer.data(snd0 + 168);
    const auto *snd0_171 = buffer.data(snd0 + 171);

    const auto *snp_44 = buffer.data(snp + 44);
    const auto *snp_46 = buffer.data(snp + 46);
    const auto *snp_47 = buffer.data(snp + 47);
    const auto *snp_50 = buffer.data(snp + 50);
    const auto *snp_52 = buffer.data(snp + 52);
    const auto *snp_53 = buffer.data(snp + 53);
    const auto *snp_55 = buffer.data(snp + 55);
    const auto *snp_56 = buffer.data(snp + 56);
    const auto *snp_58 = buffer.data(snp + 58);
    const auto *snp_59 = buffer.data(snp + 59);
    const auto *snp_60 = buffer.data(snp + 60);
    const auto *snp_61 = buffer.data(snp + 61);
    const auto *snp_62 = buffer.data(snp + 62);
    const auto *snp_63 = buffer.data(snp + 63);
    const auto *snp_64 = buffer.data(snp + 64);
    const auto *snp_65 = buffer.data(snp + 65);
    const auto *snp_67 = buffer.data(snp + 67);
    const auto *snp_68 = buffer.data(snp + 68);
    const auto *snp_69 = buffer.data(snp + 69);
    const auto *snp_70 = buffer.data(snp + 70);
    const auto *snp_71 = buffer.data(snp + 71);
    const auto *snp_72 = buffer.data(snp + 72);
    const auto *snp_73 = buffer.data(snp + 73);
    const auto *snp_74 = buffer.data(snp + 74);
    const auto *snp_75 = buffer.data(snp + 75);
    const auto *snp_76 = buffer.data(snp + 76);
    const auto *snp_77 = buffer.data(snp + 77);
    const auto *snp_79 = buffer.data(snp + 79);
    const auto *snp_80 = buffer.data(snp + 80);
    const auto *snp_81 = buffer.data(snp + 81);
    const auto *snp_82 = buffer.data(snp + 82);
    const auto *snp_83 = buffer.data(snp + 83);
    const auto *snp_84 = buffer.data(snp + 84);
    const auto *snp_85 = buffer.data(snp + 85);
    const auto *snp_86 = buffer.data(snp + 86);
    const auto *snp_88 = buffer.data(snp + 88);
    const auto *snp_89 = buffer.data(snp + 89);
    const auto *snp_90 = buffer.data(snp + 90);
    const auto *snp_91 = buffer.data(snp + 91);
    const auto *snp_92 = buffer.data(snp + 92);
    const auto *snp_93 = buffer.data(snp + 93);
    const auto *snp_94 = buffer.data(snp + 94);
    const auto *snp_95 = buffer.data(snp + 95);
    const auto *snp_96 = buffer.data(snp + 96);
    const auto *snp_97 = buffer.data(snp + 97);
    const auto *snp_98 = buffer.data(snp + 98);
    const auto *snp_99 = buffer.data(snp + 99);
    const auto *snp_100 = buffer.data(snp + 100);
    const auto *snp_101 = buffer.data(snp + 101);
    const auto *snp_103 = buffer.data(snp + 103);
    const auto *snp_104 = buffer.data(snp + 104);
    const auto *snp_105 = buffer.data(snp + 105);
    const auto *snp_106 = buffer.data(snp + 106);
    const auto *snp_107 = buffer.data(snp + 107);
    const auto *snp_108 = buffer.data(snp + 108);
    const auto *snp_109 = buffer.data(snp + 109);
    const auto *snp_110 = buffer.data(snp + 110);
    const auto *snp_112 = buffer.data(snp + 112);
    const auto *snp_113 = buffer.data(snp + 113);
    const auto *snp_114 = buffer.data(snp + 114);
    const auto *snp_115 = buffer.data(snp + 115);
    const auto *snp_116 = buffer.data(snp + 116);
    const auto *snp_117 = buffer.data(snp + 117);

    const auto *snd1_90 = buffer.data(snd1 + 90);
    const auto *snd1_93 = buffer.data(snd1 + 93);
    const auto *snd1_120 = buffer.data(snd1 + 120);
    const auto *snd1_125 = buffer.data(snd1 + 125);
    const auto *snd1_126 = buffer.data(snd1 + 126);
    const auto *snd1_129 = buffer.data(snd1 + 129);
    const auto *snd1_162 = buffer.data(snd1 + 162);
    const auto *snd1_167 = buffer.data(snd1 + 167);
    const auto *snd1_168 = buffer.data(snd1 + 168);
    const auto *snd1_171 = buffer.data(snd1 + 171);

    const auto *sos0_20 = buffer.data(sos0 + 20);
    const auto *sos0_21 = buffer.data(sos0 + 21);
    const auto *sos0_22 = buffer.data(sos0 + 22);
    const auto *sos0_23 = buffer.data(sos0 + 23);
    const auto *sos0_24 = buffer.data(sos0 + 24);
    const auto *sos0_25 = buffer.data(sos0 + 25);
    const auto *sos0_26 = buffer.data(sos0 + 26);
    const auto *sos0_27 = buffer.data(sos0 + 27);
    const auto *sos0_28 = buffer.data(sos0 + 28);
    const auto *sos0_29 = buffer.data(sos0 + 29);
    const auto *sos0_30 = buffer.data(sos0 + 30);
    const auto *sos0_31 = buffer.data(sos0 + 31);
    const auto *sos0_32 = buffer.data(sos0 + 32);
    const auto *sos0_33 = buffer.data(sos0 + 33);
    const auto *sos0_34 = buffer.data(sos0 + 34);
    const auto *sos0_35 = buffer.data(sos0 + 35);
    const auto *sos0_36 = buffer.data(sos0 + 36);
    const auto *sos0_37 = buffer.data(sos0 + 37);
    const auto *sos0_38 = buffer.data(sos0 + 38);
    const auto *sos0_39 = buffer.data(sos0 + 39);

    const auto *sos1_20 = buffer.data(sos1 + 20);
    const auto *sos1_21 = buffer.data(sos1 + 21);
    const auto *sos1_22 = buffer.data(sos1 + 22);
    const auto *sos1_23 = buffer.data(sos1 + 23);
    const auto *sos1_24 = buffer.data(sos1 + 24);
    const auto *sos1_25 = buffer.data(sos1 + 25);
    const auto *sos1_26 = buffer.data(sos1 + 26);
    const auto *sos1_27 = buffer.data(sos1 + 27);
    const auto *sos1_28 = buffer.data(sos1 + 28);
    const auto *sos1_29 = buffer.data(sos1 + 29);
    const auto *sos1_30 = buffer.data(sos1 + 30);
    const auto *sos1_31 = buffer.data(sos1 + 31);
    const auto *sos1_32 = buffer.data(sos1 + 32);
    const auto *sos1_33 = buffer.data(sos1 + 33);
    const auto *sos1_34 = buffer.data(sos1 + 34);
    const auto *sos1_35 = buffer.data(sos1 + 35);
    const auto *sos1_36 = buffer.data(sos1 + 36);
    const auto *sos1_37 = buffer.data(sos1 + 37);
    const auto *sos1_38 = buffer.data(sos1 + 38);
    const auto *sos1_39 = buffer.data(sos1 + 39);

    const auto *sop_60 = buffer.data(sop + 60);
    const auto *sop_61 = buffer.data(sop + 61);
    const auto *sop_62 = buffer.data(sop + 62);
    const auto *sop_63 = buffer.data(sop + 63);
    const auto *sop_64 = buffer.data(sop + 64);
    const auto *sop_65 = buffer.data(sop + 65);
    const auto *sop_67 = buffer.data(sop + 67);
    const auto *sop_68 = buffer.data(sop + 68);
    const auto *sop_69 = buffer.data(sop + 69);
    const auto *sop_70 = buffer.data(sop + 70);
    const auto *sop_71 = buffer.data(sop + 71);
    const auto *sop_72 = buffer.data(sop + 72);
    const auto *sop_73 = buffer.data(sop + 73);
    const auto *sop_74 = buffer.data(sop + 74);
    const auto *sop_75 = buffer.data(sop + 75);
    const auto *sop_76 = buffer.data(sop + 76);
    const auto *sop_77 = buffer.data(sop + 77);
    const auto *sop_79 = buffer.data(sop + 79);
    const auto *sop_80 = buffer.data(sop + 80);
    const auto *sop_81 = buffer.data(sop + 81);
    const auto *sop_82 = buffer.data(sop + 82);
    const auto *sop_83 = buffer.data(sop + 83);
    const auto *sop_84 = buffer.data(sop + 84);
    const auto *sop_85 = buffer.data(sop + 85);
    const auto *sop_86 = buffer.data(sop + 86);
    const auto *sop_88 = buffer.data(sop + 88);
    const auto *sop_89 = buffer.data(sop + 89);
    const auto *sop_90 = buffer.data(sop + 90);
    const auto *sop_91 = buffer.data(sop + 91);
    const auto *sop_92 = buffer.data(sop + 92);
    const auto *sop_93 = buffer.data(sop + 93);
    const auto *sop_94 = buffer.data(sop + 94);
    const auto *sop_95 = buffer.data(sop + 95);
    const auto *sop_96 = buffer.data(sop + 96);
    const auto *sop_97 = buffer.data(sop + 97);
    const auto *sop_98 = buffer.data(sop + 98);
    const auto *sop_99 = buffer.data(sop + 99);
    const auto *sop_100 = buffer.data(sop + 100);
    const auto *sop_101 = buffer.data(sop + 101);
    const auto *sop_103 = buffer.data(sop + 103);
    const auto *sop_104 = buffer.data(sop + 104);
    const auto *sop_105 = buffer.data(sop + 105);
    const auto *sop_106 = buffer.data(sop + 106);
    const auto *sop_107 = buffer.data(sop + 107);
    const auto *sop_108 = buffer.data(sop + 108);
    const auto *sop_109 = buffer.data(sop + 109);
    const auto *sop_110 = buffer.data(sop + 110);
    const auto *sop_112 = buffer.data(sop + 112);
    const auto *sop_113 = buffer.data(sop + 113);
    const auto *sop_114 = buffer.data(sop + 114);
    const auto *sop_115 = buffer.data(sop + 115);
    const auto *sop_116 = buffer.data(sop + 116);
    const auto *sop_117 = buffer.data(sop + 117);

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, pc_x, pc_y, snp_60, snp_61, \
                         snp_62, sos0_20, sos1_20, sop_60, sop_61, \
                         sop_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_13 * snp_60[k]
                   + f_1 * sos0_20[k]
                   - f_2 * sos1_20[k]
                   + f_3 * pc_x[k] * sop_60[k];

        t_121[k] = f_13 * snp_61[k]
                   + f_3 * pc_x[k] * sop_61[k];

        t_122[k] = f_13 * snp_62[k]
                   + f_3 * pc_x[k] * sop_62[k];

        t_123[k] = f_1 * sos0_20[k]
                   - f_2 * sos1_20[k]
                   + f_3 * pc_y[k] * sop_61[k];

        t_124[k] = f_3 * pc_y[k] * sop_62[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, pc_x, pc_z, snp_44, snp_63, snp_64, sos0_20, \
                         sos0_21, sos1_20, sos1_21, sop_62, sop_63, \
                         sop_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_14 * snp_44[k]
                   + f_1 * sos0_20[k]
                   - f_2 * sos1_20[k]
                   + f_3 * pc_z[k] * sop_62[k];

        t_126[k] = f_14 * snp_63[k]
                   + f_1 * sos0_21[k]
                   - f_2 * sos1_21[k]
                   + f_3 * pc_x[k] * sop_63[k];

        t_127[k] = f_14 * snp_64[k]
                   + f_3 * pc_x[k] * sop_64[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, pc_x, pc_y, pc_z, snp_46, snp_47, snp_65, \
                         sos0_21, sos1_21, sop_64, sop_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_14 * snp_65[k]
                   + f_3 * pc_x[k] * sop_65[k];

        t_129[k] = f_13 * snp_46[k]
                   + f_1 * sos0_21[k]
                   - f_2 * sos1_21[k]
                   + f_3 * pc_y[k] * sop_64[k];

        t_130[k] = f_13 * snp_47[k]
                   + f_3 * pc_y[k] * sop_65[k];

        t_131[k] = f_1 * sos0_21[k]
                   - f_2 * sos1_21[k]
                   + f_3 * pc_z[k] * sop_65[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pb_z, pc_x, pc_z, snd0_90, snd0_93, \
                         snp_67, snp_68, snd1_90, snd1_93, sop_67, \
                         sop_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = pb_z[k] * snd0_90[k]
                   - f_4 * pc_z[k] * snd1_90[k];

        t_133[k] = f_14 * snp_67[k]
                   + f_3 * pc_x[k] * sop_67[k];

        t_134[k] = f_14 * snp_68[k]
                   + f_3 * pc_x[k] * sop_68[k];

        t_135[k] = pb_z[k] * snd0_93[k]
                   - f_4 * pc_z[k] * snd1_93[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, pc_x, pc_y, pc_z, snp_47, snp_50, snp_69, \
                         sos0_22, sos0_23, sos1_22, sos1_23, sop_68, \
                         sop_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_14 * snp_50[k]
                   + f_3 * pc_y[k] * sop_68[k];

        t_137[k] = f_6 * snp_47[k]
                   + f_1 * sos0_22[k]
                   - f_2 * sos1_22[k]
                   + f_3 * pc_z[k] * sop_68[k];

        t_138[k] = f_14 * snp_69[k]
                   + f_1 * sos0_23[k]
                   - f_2 * sos1_23[k]
                   + f_3 * pc_x[k] * sop_69[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, pc_x, pc_y, snp_52, snp_53, snp_70, \
                         snp_71, sos0_23, sos1_23, sop_70, sop_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_14 * snp_70[k]
                   + f_3 * pc_x[k] * sop_70[k];

        t_140[k] = f_14 * snp_71[k]
                   + f_3 * pc_x[k] * sop_71[k];

        t_141[k] = f_12 * snp_52[k]
                   + f_1 * sos0_23[k]
                   - f_2 * sos1_23[k]
                   + f_3 * pc_y[k] * sop_70[k];

        t_142[k] = f_12 * snp_53[k]
                   + f_3 * pc_y[k] * sop_71[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, pc_x, pc_z, snp_50, snp_72, snp_73, sos0_23, \
                         sos0_24, sos1_23, sos1_24, sop_71, sop_72, \
                         sop_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_8 * snp_50[k]
                   + f_1 * sos0_23[k]
                   - f_2 * sos1_23[k]
                   + f_3 * pc_z[k] * sop_71[k];

        t_144[k] = f_14 * snp_72[k]
                   + f_1 * sos0_24[k]
                   - f_2 * sos1_24[k]
                   + f_3 * pc_x[k] * sop_72[k];

        t_145[k] = f_14 * snp_73[k]
                   + f_3 * pc_x[k] * sop_73[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pc_x, pc_y, pc_z, snp_53, snp_55, snp_56, \
                         snp_74, sos0_24, sos1_24, sop_73, sop_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_14 * snp_74[k]
                   + f_3 * pc_x[k] * sop_74[k];

        t_147[k] = f_10 * snp_55[k]
                   + f_1 * sos0_24[k]
                   - f_2 * sos1_24[k]
                   + f_3 * pc_y[k] * sop_73[k];

        t_148[k] = f_10 * snp_56[k]
                   + f_3 * pc_y[k] * sop_74[k];

        t_149[k] = f_10 * snp_53[k]
                   + f_1 * sos0_24[k]
                   - f_2 * sos1_24[k]
                   + f_3 * pc_z[k] * sop_74[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, pc_x, pc_y, snp_58, snp_75, snp_76, \
                         snp_77, sos0_25, sos1_25, sop_75, sop_76, \
                         sop_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_14 * snp_75[k]
                   + f_1 * sos0_25[k]
                   - f_2 * sos1_25[k]
                   + f_3 * pc_x[k] * sop_75[k];

        t_151[k] = f_14 * snp_76[k]
                   + f_3 * pc_x[k] * sop_76[k];

        t_152[k] = f_14 * snp_77[k]
                   + f_3 * pc_x[k] * sop_77[k];

        t_153[k] = f_8 * snp_58[k]
                   + f_1 * sos0_25[k]
                   - f_2 * sos1_25[k]
                   + f_3 * pc_y[k] * sop_76[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, pb_y, pc_y, pc_z, snd0_120, snp_56, snp_59, \
                         snd1_120, sos0_25, sos1_25, sop_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = f_8 * snp_59[k]
                   + f_3 * pc_y[k] * sop_77[k];

        t_155[k] = f_12 * snp_56[k]
                   + f_1 * sos0_25[k]
                   - f_2 * sos1_25[k]
                   + f_3 * pc_z[k] * sop_77[k];

        t_156[k] = pb_y[k] * snd0_120[k]
                   - f_4 * pc_y[k] * snd1_120[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, pc_x, pc_y, snp_61, snp_62, snp_79, \
                         snp_80, sos0_26, sos1_26, sop_79, sop_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = f_14 * snp_79[k]
                   + f_3 * pc_x[k] * sop_79[k];

        t_158[k] = f_14 * snp_80[k]
                   + f_3 * pc_x[k] * sop_80[k];

        t_159[k] = f_6 * snp_61[k]
                   + f_1 * sos0_26[k]
                   - f_2 * sos1_26[k]
                   + f_3 * pc_y[k] * sop_79[k];

        t_160[k] = f_6 * snp_62[k]
                   + f_3 * pc_y[k] * sop_80[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, pb_y, pc_x, pc_y, snd0_125, snp_81, snp_82, \
                         snd1_125, sos0_27, sos1_27, sop_81, sop_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = pb_y[k] * snd0_125[k]
                   - f_4 * pc_y[k] * snd1_125[k];

        t_162[k] = f_14 * snp_81[k]
                   + f_1 * sos0_27[k]
                   - f_2 * sos1_27[k]
                   + f_3 * pc_x[k] * sop_81[k];

        t_163[k] = f_14 * snp_82[k]
                   + f_3 * pc_x[k] * sop_82[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, pc_x, pc_y, pc_z, snp_62, snp_83, \
                         sos0_27, sos1_27, sop_82, sop_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_14 * snp_83[k]
                   + f_3 * pc_x[k] * sop_83[k];

        t_165[k] = f_1 * sos0_27[k]
                   - f_2 * sos1_27[k]
                   + f_3 * pc_y[k] * sop_82[k];

        t_166[k] = f_3 * pc_y[k] * sop_83[k];

        t_167[k] = f_13 * snp_62[k]
                   + f_1 * sos0_27[k]
                   - f_2 * sos1_27[k]
                   + f_3 * pc_z[k] * sop_83[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, pc_x, pc_y, snp_64, snp_84, snp_85, \
                         snp_86, sos0_28, sos1_28, sop_84, sop_85, \
                         sop_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_12 * snp_84[k]
                   + f_1 * sos0_28[k]
                   - f_2 * sos1_28[k]
                   + f_3 * pc_x[k] * sop_84[k];

        t_169[k] = f_12 * snp_85[k]
                   + f_3 * pc_x[k] * sop_85[k];

        t_170[k] = f_12 * snp_86[k]
                   + f_3 * pc_x[k] * sop_86[k];

        t_171[k] = f_11 * snp_64[k]
                   + f_1 * sos0_28[k]
                   - f_2 * sos1_28[k]
                   + f_3 * pc_y[k] * sop_85[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, pb_z, pc_x, pc_y, pc_z, snd0_126, snp_65, \
                         snp_88, snd1_126, sos0_28, sos1_28, sop_86, \
                         sop_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_11 * snp_65[k]
                   + f_3 * pc_y[k] * sop_86[k];

        t_173[k] = f_1 * sos0_28[k]
                   - f_2 * sos1_28[k]
                   + f_3 * pc_z[k] * sop_86[k];

        t_174[k] = pb_z[k] * snd0_126[k]
                   - f_4 * pc_z[k] * snd1_126[k];

        t_175[k] = f_12 * snp_88[k]
                   + f_3 * pc_x[k] * sop_88[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pb_z, pc_x, pc_y, pc_z, snd0_129, snp_65, \
                         snp_68, snp_89, snd1_129, sos0_29, sos1_29, \
                         sop_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_12 * snp_89[k]
                   + f_3 * pc_x[k] * sop_89[k];

        t_177[k] = pb_z[k] * snd0_129[k]
                   - f_4 * pc_z[k] * snd1_129[k];

        t_178[k] = f_13 * snp_68[k]
                   + f_3 * pc_y[k] * sop_89[k];

        t_179[k] = f_6 * snp_65[k]
                   + f_1 * sos0_29[k]
                   - f_2 * sos1_29[k]
                   + f_3 * pc_z[k] * sop_89[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pc_x, pc_y, snp_70, snp_90, snp_91, \
                         snp_92, sos0_30, sos1_30, sop_90, sop_91, \
                         sop_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_12 * snp_90[k]
                   + f_1 * sos0_30[k]
                   - f_2 * sos1_30[k]
                   + f_3 * pc_x[k] * sop_90[k];

        t_181[k] = f_12 * snp_91[k]
                   + f_3 * pc_x[k] * sop_91[k];

        t_182[k] = f_12 * snp_92[k]
                   + f_3 * pc_x[k] * sop_92[k];

        t_183[k] = f_14 * snp_70[k]
                   + f_1 * sos0_30[k]
                   - f_2 * sos1_30[k]
                   + f_3 * pc_y[k] * sop_91[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, pc_x, pc_y, pc_z, snp_68, snp_71, snp_93, \
                         sos0_30, sos0_31, sos1_30, sos1_31, sop_92, \
                         sop_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_14 * snp_71[k]
                   + f_3 * pc_y[k] * sop_92[k];

        t_185[k] = f_8 * snp_68[k]
                   + f_1 * sos0_30[k]
                   - f_2 * sos1_30[k]
                   + f_3 * pc_z[k] * sop_92[k];

        t_186[k] = f_12 * snp_93[k]
                   + f_1 * sos0_31[k]
                   - f_2 * sos1_31[k]
                   + f_3 * pc_x[k] * sop_93[k];
    }

#pragma omp simd aligned(t_187, t_188, t_189, t_190, pc_x, pc_y, snp_73, snp_74, snp_94, \
                         snp_95, sos0_31, sos1_31, sop_94, sop_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_187[k] = f_12 * snp_94[k]
                   + f_3 * pc_x[k] * sop_94[k];

        t_188[k] = f_12 * snp_95[k]
                   + f_3 * pc_x[k] * sop_95[k];

        t_189[k] = f_12 * snp_73[k]
                   + f_1 * sos0_31[k]
                   - f_2 * sos1_31[k]
                   + f_3 * pc_y[k] * sop_94[k];

        t_190[k] = f_12 * snp_74[k]
                   + f_3 * pc_y[k] * sop_95[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, pc_x, pc_z, snp_71, snp_96, snp_97, sos0_31, \
                         sos0_32, sos1_31, sos1_32, sop_95, sop_96, \
                         sop_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = f_10 * snp_71[k]
                   + f_1 * sos0_31[k]
                   - f_2 * sos1_31[k]
                   + f_3 * pc_z[k] * sop_95[k];

        t_192[k] = f_12 * snp_96[k]
                   + f_1 * sos0_32[k]
                   - f_2 * sos1_32[k]
                   + f_3 * pc_x[k] * sop_96[k];

        t_193[k] = f_12 * snp_97[k]
                   + f_3 * pc_x[k] * sop_97[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, t_197, pc_x, pc_y, pc_z, snp_74, snp_76, snp_77, \
                         snp_98, sos0_32, sos1_32, sop_97, sop_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = f_12 * snp_98[k]
                   + f_3 * pc_x[k] * sop_98[k];

        t_195[k] = f_10 * snp_76[k]
                   + f_1 * sos0_32[k]
                   - f_2 * sos1_32[k]
                   + f_3 * pc_y[k] * sop_97[k];

        t_196[k] = f_10 * snp_77[k]
                   + f_3 * pc_y[k] * sop_98[k];

        t_197[k] = f_12 * snp_74[k]
                   + f_1 * sos0_32[k]
                   - f_2 * sos1_32[k]
                   + f_3 * pc_z[k] * sop_98[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, pc_x, pc_y, snp_79, snp_99, snp_100, \
                         snp_101, sos0_33, sos1_33, sop_99, sop_100, \
                         sop_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_12 * snp_99[k]
                   + f_1 * sos0_33[k]
                   - f_2 * sos1_33[k]
                   + f_3 * pc_x[k] * sop_99[k];

        t_199[k] = f_12 * snp_100[k]
                   + f_3 * pc_x[k] * sop_100[k];

        t_200[k] = f_12 * snp_101[k]
                   + f_3 * pc_x[k] * sop_101[k];

        t_201[k] = f_8 * snp_79[k]
                   + f_1 * sos0_33[k]
                   - f_2 * sos1_33[k]
                   + f_3 * pc_y[k] * sop_100[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, pb_y, pc_y, pc_z, snd0_162, snp_77, snp_80, \
                         snd1_162, sos0_33, sos1_33, sop_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = f_8 * snp_80[k]
                   + f_3 * pc_y[k] * sop_101[k];

        t_203[k] = f_14 * snp_77[k]
                   + f_1 * sos0_33[k]
                   - f_2 * sos1_33[k]
                   + f_3 * pc_z[k] * sop_101[k];

        t_204[k] = pb_y[k] * snd0_162[k]
                   - f_4 * pc_y[k] * snd1_162[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, pc_x, pc_y, snp_82, snp_83, snp_103, \
                         snp_104, sos0_34, sos1_34, sop_103, sop_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = f_12 * snp_103[k]
                   + f_3 * pc_x[k] * sop_103[k];

        t_206[k] = f_12 * snp_104[k]
                   + f_3 * pc_x[k] * sop_104[k];

        t_207[k] = f_6 * snp_82[k]
                   + f_1 * sos0_34[k]
                   - f_2 * sos1_34[k]
                   + f_3 * pc_y[k] * sop_103[k];

        t_208[k] = f_6 * snp_83[k]
                   + f_3 * pc_y[k] * sop_104[k];
    }

#pragma omp simd aligned(t_209, t_210, t_211, pb_y, pc_x, pc_y, snd0_167, snp_105, snp_106, \
                         snd1_167, sos0_35, sos1_35, sop_105, sop_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_209[k] = pb_y[k] * snd0_167[k]
                   - f_4 * pc_y[k] * snd1_167[k];

        t_210[k] = f_12 * snp_105[k]
                   + f_1 * sos0_35[k]
                   - f_2 * sos1_35[k]
                   + f_3 * pc_x[k] * sop_105[k];

        t_211[k] = f_12 * snp_106[k]
                   + f_3 * pc_x[k] * sop_106[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, pc_x, pc_y, pc_z, snp_83, snp_107, \
                         sos0_35, sos1_35, sop_106, sop_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_12 * snp_107[k]
                   + f_3 * pc_x[k] * sop_107[k];

        t_213[k] = f_1 * sos0_35[k]
                   - f_2 * sos1_35[k]
                   + f_3 * pc_y[k] * sop_106[k];

        t_214[k] = f_3 * pc_y[k] * sop_107[k];

        t_215[k] = f_11 * snp_83[k]
                   + f_1 * sos0_35[k]
                   - f_2 * sos1_35[k]
                   + f_3 * pc_z[k] * sop_107[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, pc_x, pc_y, snp_85, snp_108, snp_109, \
                         snp_110, sos0_36, sos1_36, sop_108, sop_109, \
                         sop_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = f_10 * snp_108[k]
                   + f_1 * sos0_36[k]
                   - f_2 * sos1_36[k]
                   + f_3 * pc_x[k] * sop_108[k];

        t_217[k] = f_10 * snp_109[k]
                   + f_3 * pc_x[k] * sop_109[k];

        t_218[k] = f_10 * snp_110[k]
                   + f_3 * pc_x[k] * sop_110[k];

        t_219[k] = f_9 * snp_85[k]
                   + f_1 * sos0_36[k]
                   - f_2 * sos1_36[k]
                   + f_3 * pc_y[k] * sop_109[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, pb_z, pc_x, pc_y, pc_z, snd0_168, snp_86, \
                         snp_112, snd1_168, sos0_36, sos1_36, sop_110, \
                         sop_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = f_9 * snp_86[k]
                   + f_3 * pc_y[k] * sop_110[k];

        t_221[k] = f_1 * sos0_36[k]
                   - f_2 * sos1_36[k]
                   + f_3 * pc_z[k] * sop_110[k];

        t_222[k] = pb_z[k] * snd0_168[k]
                   - f_4 * pc_z[k] * snd1_168[k];

        t_223[k] = f_10 * snp_112[k]
                   + f_3 * pc_x[k] * sop_112[k];
    }

#pragma omp simd aligned(t_224, t_225, t_226, t_227, pb_z, pc_x, pc_y, pc_z, snd0_171, snp_86, \
                         snp_89, snp_113, snd1_171, sos0_37, sos1_37, \
                         sop_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_224[k] = f_10 * snp_113[k]
                   + f_3 * pc_x[k] * sop_113[k];

        t_225[k] = pb_z[k] * snd0_171[k]
                   - f_4 * pc_z[k] * snd1_171[k];

        t_226[k] = f_11 * snp_89[k]
                   + f_3 * pc_y[k] * sop_113[k];

        t_227[k] = f_6 * snp_86[k]
                   + f_1 * sos0_37[k]
                   - f_2 * sos1_37[k]
                   + f_3 * pc_z[k] * sop_113[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, t_231, pc_x, pc_y, snp_91, snp_114, snp_115, \
                         snp_116, sos0_38, sos1_38, sop_114, sop_115, \
                         sop_116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_10 * snp_114[k]
                   + f_1 * sos0_38[k]
                   - f_2 * sos1_38[k]
                   + f_3 * pc_x[k] * sop_114[k];

        t_229[k] = f_10 * snp_115[k]
                   + f_3 * pc_x[k] * sop_115[k];

        t_230[k] = f_10 * snp_116[k]
                   + f_3 * pc_x[k] * sop_116[k];

        t_231[k] = f_13 * snp_91[k]
                   + f_1 * sos0_38[k]
                   - f_2 * sos1_38[k]
                   + f_3 * pc_y[k] * sop_115[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, pc_x, pc_y, pc_z, snp_89, snp_92, snp_117, \
                         sos0_38, sos0_39, sos1_38, sos1_39, sop_116, \
                         sop_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = f_13 * snp_92[k]
                   + f_3 * pc_y[k] * sop_116[k];

        t_233[k] = f_8 * snp_89[k]
                   + f_1 * sos0_38[k]
                   - f_2 * sos1_38[k]
                   + f_3 * pc_z[k] * sop_116[k];

        t_234[k] = f_10 * snp_117[k]
                   + f_1 * sos0_39[k]
                   - f_2 * sos1_39[k]
                   + f_3 * pc_x[k] * sop_117[k];
    }
}

static auto
compute_prim_sod_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t snd0,
                                                          const size_t snp, const size_t snd1,
                                                          const size_t sos0, const size_t sos1,
                                                          const size_t sop, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 0.5 / gamma;
    const auto f_2 = 0.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 5.0 / q;
    const auto f_6 = 0.5 / q;
    const auto f_7 = 4.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 4.0 / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 3.5 / q;
    const auto f_12 = 2.0 / q;
    const auto f_13 = 3.0 / q;
    const auto f_14 = 2.5 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *snd0_210 = buffer.data(snd0 + 210);
    const auto *snd0_215 = buffer.data(snd0 + 215);
    const auto *snd0_216 = buffer.data(snd0 + 216);
    const auto *snd0_219 = buffer.data(snd0 + 219);
    const auto *snd0_264 = buffer.data(snd0 + 264);
    const auto *snd0_269 = buffer.data(snd0 + 269);
    const auto *snd0_270 = buffer.data(snd0 + 270);
    const auto *snd0_330 = buffer.data(snd0 + 330);
    const auto *snd0_333 = buffer.data(snd0 + 333);
    const auto *snd0_335 = buffer.data(snd0 + 335);
    const auto *snd0_339 = buffer.data(snd0 + 339);
    const auto *snd0_341 = buffer.data(snd0 + 341);
    const auto *snd0_342 = buffer.data(snd0 + 342);
    const auto *snd0_345 = buffer.data(snd0 + 345);
    const auto *snd0_347 = buffer.data(snd0 + 347);
    const auto *snd0_348 = buffer.data(snd0 + 348);
    const auto *snd0_351 = buffer.data(snd0 + 351);

    const auto *snp_92 = buffer.data(snp + 92);
    const auto *snp_94 = buffer.data(snp + 94);
    const auto *snp_95 = buffer.data(snp + 95);
    const auto *snp_97 = buffer.data(snp + 97);
    const auto *snp_98 = buffer.data(snp + 98);
    const auto *snp_100 = buffer.data(snp + 100);
    const auto *snp_101 = buffer.data(snp + 101);
    const auto *snp_103 = buffer.data(snp + 103);
    const auto *snp_104 = buffer.data(snp + 104);
    const auto *snp_106 = buffer.data(snp + 106);
    const auto *snp_107 = buffer.data(snp + 107);
    const auto *snp_109 = buffer.data(snp + 109);
    const auto *snp_110 = buffer.data(snp + 110);
    const auto *snp_113 = buffer.data(snp + 113);
    const auto *snp_115 = buffer.data(snp + 115);
    const auto *snp_116 = buffer.data(snp + 116);
    const auto *snp_118 = buffer.data(snp + 118);
    const auto *snp_119 = buffer.data(snp + 119);
    const auto *snp_120 = buffer.data(snp + 120);
    const auto *snp_121 = buffer.data(snp + 121);
    const auto *snp_122 = buffer.data(snp + 122);
    const auto *snp_123 = buffer.data(snp + 123);
    const auto *snp_124 = buffer.data(snp + 124);
    const auto *snp_125 = buffer.data(snp + 125);
    const auto *snp_126 = buffer.data(snp + 126);
    const auto *snp_127 = buffer.data(snp + 127);
    const auto *snp_128 = buffer.data(snp + 128);
    const auto *snp_130 = buffer.data(snp + 130);
    const auto *snp_131 = buffer.data(snp + 131);
    const auto *snp_132 = buffer.data(snp + 132);
    const auto *snp_133 = buffer.data(snp + 133);
    const auto *snp_134 = buffer.data(snp + 134);
    const auto *snp_135 = buffer.data(snp + 135);
    const auto *snp_136 = buffer.data(snp + 136);
    const auto *snp_137 = buffer.data(snp + 137);
    const auto *snp_139 = buffer.data(snp + 139);
    const auto *snp_140 = buffer.data(snp + 140);
    const auto *snp_141 = buffer.data(snp + 141);
    const auto *snp_142 = buffer.data(snp + 142);
    const auto *snp_143 = buffer.data(snp + 143);
    const auto *snp_144 = buffer.data(snp + 144);
    const auto *snp_145 = buffer.data(snp + 145);
    const auto *snp_146 = buffer.data(snp + 146);
    const auto *snp_147 = buffer.data(snp + 147);
    const auto *snp_148 = buffer.data(snp + 148);
    const auto *snp_149 = buffer.data(snp + 149);
    const auto *snp_150 = buffer.data(snp + 150);
    const auto *snp_151 = buffer.data(snp + 151);
    const auto *snp_152 = buffer.data(snp + 152);
    const auto *snp_153 = buffer.data(snp + 153);
    const auto *snp_154 = buffer.data(snp + 154);
    const auto *snp_155 = buffer.data(snp + 155);
    const auto *snp_156 = buffer.data(snp + 156);
    const auto *snp_157 = buffer.data(snp + 157);
    const auto *snp_158 = buffer.data(snp + 158);
    const auto *snp_160 = buffer.data(snp + 160);
    const auto *snp_161 = buffer.data(snp + 161);
    const auto *snp_162 = buffer.data(snp + 162);
    const auto *snp_163 = buffer.data(snp + 163);
    const auto *snp_164 = buffer.data(snp + 164);
    const auto *snp_165 = buffer.data(snp + 165);
    const auto *snp_166 = buffer.data(snp + 166);
    const auto *snp_167 = buffer.data(snp + 167);
    const auto *snp_169 = buffer.data(snp + 169);
    const auto *snp_170 = buffer.data(snp + 170);
    const auto *snp_171 = buffer.data(snp + 171);
    const auto *snp_172 = buffer.data(snp + 172);
    const auto *snp_173 = buffer.data(snp + 173);
    const auto *snp_174 = buffer.data(snp + 174);
    const auto *snp_175 = buffer.data(snp + 175);
    const auto *snp_176 = buffer.data(snp + 176);

    const auto *snd1_210 = buffer.data(snd1 + 210);
    const auto *snd1_215 = buffer.data(snd1 + 215);
    const auto *snd1_216 = buffer.data(snd1 + 216);
    const auto *snd1_219 = buffer.data(snd1 + 219);
    const auto *snd1_264 = buffer.data(snd1 + 264);
    const auto *snd1_269 = buffer.data(snd1 + 269);
    const auto *snd1_270 = buffer.data(snd1 + 270);
    const auto *snd1_330 = buffer.data(snd1 + 330);
    const auto *snd1_333 = buffer.data(snd1 + 333);
    const auto *snd1_335 = buffer.data(snd1 + 335);
    const auto *snd1_339 = buffer.data(snd1 + 339);
    const auto *snd1_341 = buffer.data(snd1 + 341);
    const auto *snd1_342 = buffer.data(snd1 + 342);
    const auto *snd1_345 = buffer.data(snd1 + 345);
    const auto *snd1_347 = buffer.data(snd1 + 347);
    const auto *snd1_348 = buffer.data(snd1 + 348);
    const auto *snd1_351 = buffer.data(snd1 + 351);

    const auto *sos0_39 = buffer.data(sos0 + 39);
    const auto *sos0_40 = buffer.data(sos0 + 40);
    const auto *sos0_41 = buffer.data(sos0 + 41);
    const auto *sos0_42 = buffer.data(sos0 + 42);
    const auto *sos0_43 = buffer.data(sos0 + 43);
    const auto *sos0_44 = buffer.data(sos0 + 44);
    const auto *sos0_45 = buffer.data(sos0 + 45);
    const auto *sos0_46 = buffer.data(sos0 + 46);
    const auto *sos0_47 = buffer.data(sos0 + 47);
    const auto *sos0_48 = buffer.data(sos0 + 48);
    const auto *sos0_49 = buffer.data(sos0 + 49);
    const auto *sos0_50 = buffer.data(sos0 + 50);
    const auto *sos0_51 = buffer.data(sos0 + 51);
    const auto *sos0_52 = buffer.data(sos0 + 52);
    const auto *sos0_53 = buffer.data(sos0 + 53);
    const auto *sos0_54 = buffer.data(sos0 + 54);

    const auto *sos1_39 = buffer.data(sos1 + 39);
    const auto *sos1_40 = buffer.data(sos1 + 40);
    const auto *sos1_41 = buffer.data(sos1 + 41);
    const auto *sos1_42 = buffer.data(sos1 + 42);
    const auto *sos1_43 = buffer.data(sos1 + 43);
    const auto *sos1_44 = buffer.data(sos1 + 44);
    const auto *sos1_45 = buffer.data(sos1 + 45);
    const auto *sos1_46 = buffer.data(sos1 + 46);
    const auto *sos1_47 = buffer.data(sos1 + 47);
    const auto *sos1_48 = buffer.data(sos1 + 48);
    const auto *sos1_49 = buffer.data(sos1 + 49);
    const auto *sos1_50 = buffer.data(sos1 + 50);
    const auto *sos1_51 = buffer.data(sos1 + 51);
    const auto *sos1_52 = buffer.data(sos1 + 52);
    const auto *sos1_53 = buffer.data(sos1 + 53);
    const auto *sos1_54 = buffer.data(sos1 + 54);

    const auto *sop_118 = buffer.data(sop + 118);
    const auto *sop_119 = buffer.data(sop + 119);
    const auto *sop_120 = buffer.data(sop + 120);
    const auto *sop_121 = buffer.data(sop + 121);
    const auto *sop_122 = buffer.data(sop + 122);
    const auto *sop_123 = buffer.data(sop + 123);
    const auto *sop_124 = buffer.data(sop + 124);
    const auto *sop_125 = buffer.data(sop + 125);
    const auto *sop_126 = buffer.data(sop + 126);
    const auto *sop_127 = buffer.data(sop + 127);
    const auto *sop_128 = buffer.data(sop + 128);
    const auto *sop_130 = buffer.data(sop + 130);
    const auto *sop_131 = buffer.data(sop + 131);
    const auto *sop_132 = buffer.data(sop + 132);
    const auto *sop_133 = buffer.data(sop + 133);
    const auto *sop_134 = buffer.data(sop + 134);
    const auto *sop_135 = buffer.data(sop + 135);
    const auto *sop_136 = buffer.data(sop + 136);
    const auto *sop_137 = buffer.data(sop + 137);
    const auto *sop_139 = buffer.data(sop + 139);
    const auto *sop_140 = buffer.data(sop + 140);
    const auto *sop_141 = buffer.data(sop + 141);
    const auto *sop_142 = buffer.data(sop + 142);
    const auto *sop_143 = buffer.data(sop + 143);
    const auto *sop_144 = buffer.data(sop + 144);
    const auto *sop_145 = buffer.data(sop + 145);
    const auto *sop_146 = buffer.data(sop + 146);
    const auto *sop_147 = buffer.data(sop + 147);
    const auto *sop_148 = buffer.data(sop + 148);
    const auto *sop_149 = buffer.data(sop + 149);
    const auto *sop_150 = buffer.data(sop + 150);
    const auto *sop_151 = buffer.data(sop + 151);
    const auto *sop_152 = buffer.data(sop + 152);
    const auto *sop_153 = buffer.data(sop + 153);
    const auto *sop_154 = buffer.data(sop + 154);
    const auto *sop_155 = buffer.data(sop + 155);
    const auto *sop_156 = buffer.data(sop + 156);
    const auto *sop_157 = buffer.data(sop + 157);
    const auto *sop_158 = buffer.data(sop + 158);
    const auto *sop_160 = buffer.data(sop + 160);
    const auto *sop_161 = buffer.data(sop + 161);
    const auto *sop_162 = buffer.data(sop + 162);
    const auto *sop_163 = buffer.data(sop + 163);
    const auto *sop_164 = buffer.data(sop + 164);
    const auto *sop_166 = buffer.data(sop + 166);
    const auto *sop_167 = buffer.data(sop + 167);
    const auto *sop_169 = buffer.data(sop + 169);
    const auto *sop_170 = buffer.data(sop + 170);
    const auto *sop_172 = buffer.data(sop + 172);
    const auto *sop_173 = buffer.data(sop + 173);
    const auto *sop_175 = buffer.data(sop + 175);
    const auto *sop_176 = buffer.data(sop + 176);

#pragma omp simd aligned(t_235, t_236, t_237, t_238, pc_x, pc_y, snp_94, snp_95, snp_118, \
                         snp_119, sos0_39, sos1_39, sop_118, sop_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = f_10 * snp_118[k]
                   + f_3 * pc_x[k] * sop_118[k];

        t_236[k] = f_10 * snp_119[k]
                   + f_3 * pc_x[k] * sop_119[k];

        t_237[k] = f_14 * snp_94[k]
                   + f_1 * sos0_39[k]
                   - f_2 * sos1_39[k]
                   + f_3 * pc_y[k] * sop_118[k];

        t_238[k] = f_14 * snp_95[k]
                   + f_3 * pc_y[k] * sop_119[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, pc_x, pc_z, snp_92, snp_120, snp_121, sos0_39, \
                         sos0_40, sos1_39, sos1_40, sop_119, sop_120, \
                         sop_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_10 * snp_92[k]
                   + f_1 * sos0_39[k]
                   - f_2 * sos1_39[k]
                   + f_3 * pc_z[k] * sop_119[k];

        t_240[k] = f_10 * snp_120[k]
                   + f_1 * sos0_40[k]
                   - f_2 * sos1_40[k]
                   + f_3 * pc_x[k] * sop_120[k];

        t_241[k] = f_10 * snp_121[k]
                   + f_3 * pc_x[k] * sop_121[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, t_245, pc_x, pc_y, pc_z, snp_95, snp_97, snp_98, \
                         snp_122, sos0_40, sos1_40, sop_121, sop_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = f_10 * snp_122[k]
                   + f_3 * pc_x[k] * sop_122[k];

        t_243[k] = f_12 * snp_97[k]
                   + f_1 * sos0_40[k]
                   - f_2 * sos1_40[k]
                   + f_3 * pc_y[k] * sop_121[k];

        t_244[k] = f_12 * snp_98[k]
                   + f_3 * pc_y[k] * sop_122[k];

        t_245[k] = f_12 * snp_95[k]
                   + f_1 * sos0_40[k]
                   - f_2 * sos1_40[k]
                   + f_3 * pc_z[k] * sop_122[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, t_249, pc_x, pc_y, snp_100, snp_123, snp_124, \
                         snp_125, sos0_41, sos1_41, sop_123, sop_124, \
                         sop_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = f_10 * snp_123[k]
                   + f_1 * sos0_41[k]
                   - f_2 * sos1_41[k]
                   + f_3 * pc_x[k] * sop_123[k];

        t_247[k] = f_10 * snp_124[k]
                   + f_3 * pc_x[k] * sop_124[k];

        t_248[k] = f_10 * snp_125[k]
                   + f_3 * pc_x[k] * sop_125[k];

        t_249[k] = f_10 * snp_100[k]
                   + f_1 * sos0_41[k]
                   - f_2 * sos1_41[k]
                   + f_3 * pc_y[k] * sop_124[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, pc_x, pc_y, pc_z, snp_98, snp_101, snp_126, \
                         sos0_41, sos0_42, sos1_41, sos1_42, sop_125, \
                         sop_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = f_10 * snp_101[k]
                   + f_3 * pc_y[k] * sop_125[k];

        t_251[k] = f_14 * snp_98[k]
                   + f_1 * sos0_41[k]
                   - f_2 * sos1_41[k]
                   + f_3 * pc_z[k] * sop_125[k];

        t_252[k] = f_10 * snp_126[k]
                   + f_1 * sos0_42[k]
                   - f_2 * sos1_42[k]
                   + f_3 * pc_x[k] * sop_126[k];
    }

#pragma omp simd aligned(t_253, t_254, t_255, t_256, pc_x, pc_y, snp_103, snp_104, snp_127, \
                         snp_128, sos0_42, sos1_42, sop_127, sop_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_253[k] = f_10 * snp_127[k]
                   + f_3 * pc_x[k] * sop_127[k];

        t_254[k] = f_10 * snp_128[k]
                   + f_3 * pc_x[k] * sop_128[k];

        t_255[k] = f_8 * snp_103[k]
                   + f_1 * sos0_42[k]
                   - f_2 * sos1_42[k]
                   + f_3 * pc_y[k] * sop_127[k];

        t_256[k] = f_8 * snp_104[k]
                   + f_3 * pc_y[k] * sop_128[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, pb_y, pc_x, pc_y, pc_z, snd0_210, snp_101, \
                         snp_130, snd1_210, sos0_42, sos1_42, sop_128, \
                         sop_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_13 * snp_101[k]
                   + f_1 * sos0_42[k]
                   - f_2 * sos1_42[k]
                   + f_3 * pc_z[k] * sop_128[k];

        t_258[k] = pb_y[k] * snd0_210[k]
                   - f_4 * pc_y[k] * snd1_210[k];

        t_259[k] = f_10 * snp_130[k]
                   + f_3 * pc_x[k] * sop_130[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pb_y, pc_x, pc_y, snd0_215, snp_106, \
                         snp_107, snp_131, snd1_215, sos0_43, sos1_43, sop_130, \
                         sop_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_10 * snp_131[k]
                   + f_3 * pc_x[k] * sop_131[k];

        t_261[k] = f_6 * snp_106[k]
                   + f_1 * sos0_43[k]
                   - f_2 * sos1_43[k]
                   + f_3 * pc_y[k] * sop_130[k];

        t_262[k] = f_6 * snp_107[k]
                   + f_3 * pc_y[k] * sop_131[k];

        t_263[k] = pb_y[k] * snd0_215[k]
                   - f_4 * pc_y[k] * snd1_215[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, t_268, pc_x, pc_y, snp_132, snp_133, \
                         snp_134, sos0_44, sos1_44, sop_132, sop_133, \
                         sop_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_10 * snp_132[k]
                   + f_1 * sos0_44[k]
                   - f_2 * sos1_44[k]
                   + f_3 * pc_x[k] * sop_132[k];

        t_265[k] = f_10 * snp_133[k]
                   + f_3 * pc_x[k] * sop_133[k];

        t_266[k] = f_10 * snp_134[k]
                   + f_3 * pc_x[k] * sop_134[k];

        t_267[k] = f_1 * sos0_44[k]
                   - f_2 * sos1_44[k]
                   + f_3 * pc_y[k] * sop_133[k];

        t_268[k] = f_3 * pc_y[k] * sop_134[k];
    }

#pragma omp simd aligned(t_269, t_270, t_271, pc_x, pc_z, snp_107, snp_135, snp_136, sos0_44, \
                         sos0_45, sos1_44, sos1_45, sop_134, sop_135, \
                         sop_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_269[k] = f_9 * snp_107[k]
                   + f_1 * sos0_44[k]
                   - f_2 * sos1_44[k]
                   + f_3 * pc_z[k] * sop_134[k];

        t_270[k] = f_8 * snp_135[k]
                   + f_1 * sos0_45[k]
                   - f_2 * sos1_45[k]
                   + f_3 * pc_x[k] * sop_135[k];

        t_271[k] = f_8 * snp_136[k]
                   + f_3 * pc_x[k] * sop_136[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, pc_x, pc_y, pc_z, snp_109, snp_110, \
                         snp_137, sos0_45, sos1_45, sop_136, sop_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_8 * snp_137[k]
                   + f_3 * pc_x[k] * sop_137[k];

        t_273[k] = f_7 * snp_109[k]
                   + f_1 * sos0_45[k]
                   - f_2 * sos1_45[k]
                   + f_3 * pc_y[k] * sop_136[k];

        t_274[k] = f_7 * snp_110[k]
                   + f_3 * pc_y[k] * sop_137[k];

        t_275[k] = f_1 * sos0_45[k]
                   - f_2 * sos1_45[k]
                   + f_3 * pc_z[k] * sop_137[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pb_z, pc_x, pc_z, snd0_216, snd0_219, \
                         snp_139, snp_140, snd1_216, snd1_219, sop_139, \
                         sop_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = pb_z[k] * snd0_216[k]
                   - f_4 * pc_z[k] * snd1_216[k];

        t_277[k] = f_8 * snp_139[k]
                   + f_3 * pc_x[k] * sop_139[k];

        t_278[k] = f_8 * snp_140[k]
                   + f_3 * pc_x[k] * sop_140[k];

        t_279[k] = pb_z[k] * snd0_219[k]
                   - f_4 * pc_z[k] * snd1_219[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, pc_x, pc_y, pc_z, snp_110, snp_113, snp_141, \
                         sos0_46, sos0_47, sos1_46, sos1_47, sop_140, \
                         sop_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_9 * snp_113[k]
                   + f_3 * pc_y[k] * sop_140[k];

        t_281[k] = f_6 * snp_110[k]
                   + f_1 * sos0_46[k]
                   - f_2 * sos1_46[k]
                   + f_3 * pc_z[k] * sop_140[k];

        t_282[k] = f_8 * snp_141[k]
                   + f_1 * sos0_47[k]
                   - f_2 * sos1_47[k]
                   + f_3 * pc_x[k] * sop_141[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, t_286, pc_x, pc_y, snp_115, snp_116, snp_142, \
                         snp_143, sos0_47, sos1_47, sop_142, sop_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_8 * snp_142[k]
                   + f_3 * pc_x[k] * sop_142[k];

        t_284[k] = f_8 * snp_143[k]
                   + f_3 * pc_x[k] * sop_143[k];

        t_285[k] = f_11 * snp_115[k]
                   + f_1 * sos0_47[k]
                   - f_2 * sos1_47[k]
                   + f_3 * pc_y[k] * sop_142[k];

        t_286[k] = f_11 * snp_116[k]
                   + f_3 * pc_y[k] * sop_143[k];
    }

#pragma omp simd aligned(t_287, t_288, t_289, pc_x, pc_z, snp_113, snp_144, snp_145, sos0_47, \
                         sos0_48, sos1_47, sos1_48, sop_143, sop_144, \
                         sop_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_287[k] = f_8 * snp_113[k]
                   + f_1 * sos0_47[k]
                   - f_2 * sos1_47[k]
                   + f_3 * pc_z[k] * sop_143[k];

        t_288[k] = f_8 * snp_144[k]
                   + f_1 * sos0_48[k]
                   - f_2 * sos1_48[k]
                   + f_3 * pc_x[k] * sop_144[k];

        t_289[k] = f_8 * snp_145[k]
                   + f_3 * pc_x[k] * sop_145[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, pc_x, pc_y, pc_z, snp_116, snp_118, \
                         snp_119, snp_146, sos0_48, sos1_48, sop_145, \
                         sop_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = f_8 * snp_146[k]
                   + f_3 * pc_x[k] * sop_146[k];

        t_291[k] = f_13 * snp_118[k]
                   + f_1 * sos0_48[k]
                   - f_2 * sos1_48[k]
                   + f_3 * pc_y[k] * sop_145[k];

        t_292[k] = f_13 * snp_119[k]
                   + f_3 * pc_y[k] * sop_146[k];

        t_293[k] = f_10 * snp_116[k]
                   + f_1 * sos0_48[k]
                   - f_2 * sos1_48[k]
                   + f_3 * pc_z[k] * sop_146[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, pc_x, pc_y, snp_121, snp_147, snp_148, \
                         snp_149, sos0_49, sos1_49, sop_147, sop_148, \
                         sop_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_8 * snp_147[k]
                   + f_1 * sos0_49[k]
                   - f_2 * sos1_49[k]
                   + f_3 * pc_x[k] * sop_147[k];

        t_295[k] = f_8 * snp_148[k]
                   + f_3 * pc_x[k] * sop_148[k];

        t_296[k] = f_8 * snp_149[k]
                   + f_3 * pc_x[k] * sop_149[k];

        t_297[k] = f_14 * snp_121[k]
                   + f_1 * sos0_49[k]
                   - f_2 * sos1_49[k]
                   + f_3 * pc_y[k] * sop_148[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, pc_x, pc_y, pc_z, snp_119, snp_122, snp_150, \
                         sos0_49, sos0_50, sos1_49, sos1_50, sop_149, \
                         sop_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = f_14 * snp_122[k]
                   + f_3 * pc_y[k] * sop_149[k];

        t_299[k] = f_12 * snp_119[k]
                   + f_1 * sos0_49[k]
                   - f_2 * sos1_49[k]
                   + f_3 * pc_z[k] * sop_149[k];

        t_300[k] = f_8 * snp_150[k]
                   + f_1 * sos0_50[k]
                   - f_2 * sos1_50[k]
                   + f_3 * pc_x[k] * sop_150[k];
    }

#pragma omp simd aligned(t_301, t_302, t_303, t_304, pc_x, pc_y, snp_124, snp_125, snp_151, \
                         snp_152, sos0_50, sos1_50, sop_151, sop_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_301[k] = f_8 * snp_151[k]
                   + f_3 * pc_x[k] * sop_151[k];

        t_302[k] = f_8 * snp_152[k]
                   + f_3 * pc_x[k] * sop_152[k];

        t_303[k] = f_12 * snp_124[k]
                   + f_1 * sos0_50[k]
                   - f_2 * sos1_50[k]
                   + f_3 * pc_y[k] * sop_151[k];

        t_304[k] = f_12 * snp_125[k]
                   + f_3 * pc_y[k] * sop_152[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, pc_x, pc_z, snp_122, snp_153, snp_154, sos0_50, \
                         sos0_51, sos1_50, sos1_51, sop_152, sop_153, \
                         sop_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = f_14 * snp_122[k]
                   + f_1 * sos0_50[k]
                   - f_2 * sos1_50[k]
                   + f_3 * pc_z[k] * sop_152[k];

        t_306[k] = f_8 * snp_153[k]
                   + f_1 * sos0_51[k]
                   - f_2 * sos1_51[k]
                   + f_3 * pc_x[k] * sop_153[k];

        t_307[k] = f_8 * snp_154[k]
                   + f_3 * pc_x[k] * sop_154[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, t_311, pc_x, pc_y, pc_z, snp_125, snp_127, \
                         snp_128, snp_155, sos0_51, sos1_51, sop_154, \
                         sop_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = f_8 * snp_155[k]
                   + f_3 * pc_x[k] * sop_155[k];

        t_309[k] = f_10 * snp_127[k]
                   + f_1 * sos0_51[k]
                   - f_2 * sos1_51[k]
                   + f_3 * pc_y[k] * sop_154[k];

        t_310[k] = f_10 * snp_128[k]
                   + f_3 * pc_y[k] * sop_155[k];

        t_311[k] = f_13 * snp_125[k]
                   + f_1 * sos0_51[k]
                   - f_2 * sos1_51[k]
                   + f_3 * pc_z[k] * sop_155[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, t_315, pc_x, pc_y, snp_130, snp_156, snp_157, \
                         snp_158, sos0_52, sos1_52, sop_156, sop_157, \
                         sop_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = f_8 * snp_156[k]
                   + f_1 * sos0_52[k]
                   - f_2 * sos1_52[k]
                   + f_3 * pc_x[k] * sop_156[k];

        t_313[k] = f_8 * snp_157[k]
                   + f_3 * pc_x[k] * sop_157[k];

        t_314[k] = f_8 * snp_158[k]
                   + f_3 * pc_x[k] * sop_158[k];

        t_315[k] = f_8 * snp_130[k]
                   + f_1 * sos0_52[k]
                   - f_2 * sos1_52[k]
                   + f_3 * pc_y[k] * sop_157[k];
    }

#pragma omp simd aligned(t_316, t_317, t_318, pb_y, pc_y, pc_z, snd0_264, snp_128, snp_131, \
                         snd1_264, sos0_52, sos1_52, sop_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_316[k] = f_8 * snp_131[k]
                   + f_3 * pc_y[k] * sop_158[k];

        t_317[k] = f_11 * snp_128[k]
                   + f_1 * sos0_52[k]
                   - f_2 * sos1_52[k]
                   + f_3 * pc_z[k] * sop_158[k];

        t_318[k] = pb_y[k] * snd0_264[k]
                   - f_4 * pc_y[k] * snd1_264[k];
    }

#pragma omp simd aligned(t_319, t_320, t_321, t_322, pc_x, pc_y, snp_133, snp_134, snp_160, \
                         snp_161, sos0_53, sos1_53, sop_160, sop_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_319[k] = f_8 * snp_160[k]
                   + f_3 * pc_x[k] * sop_160[k];

        t_320[k] = f_8 * snp_161[k]
                   + f_3 * pc_x[k] * sop_161[k];

        t_321[k] = f_6 * snp_133[k]
                   + f_1 * sos0_53[k]
                   - f_2 * sos1_53[k]
                   + f_3 * pc_y[k] * sop_160[k];

        t_322[k] = f_6 * snp_134[k]
                   + f_3 * pc_y[k] * sop_161[k];
    }

#pragma omp simd aligned(t_323, t_324, t_325, pb_y, pc_x, pc_y, snd0_269, snp_162, snp_163, \
                         snd1_269, sos0_54, sos1_54, sop_162, sop_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_323[k] = pb_y[k] * snd0_269[k]
                   - f_4 * pc_y[k] * snd1_269[k];

        t_324[k] = f_8 * snp_162[k]
                   + f_1 * sos0_54[k]
                   - f_2 * sos1_54[k]
                   + f_3 * pc_x[k] * sop_162[k];

        t_325[k] = f_8 * snp_163[k]
                   + f_3 * pc_x[k] * sop_163[k];
    }

#pragma omp simd aligned(t_326, t_327, t_328, t_329, pc_x, pc_y, pc_z, snp_134, snp_164, \
                         sos0_54, sos1_54, sop_163, sop_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_326[k] = f_8 * snp_164[k]
                   + f_3 * pc_x[k] * sop_164[k];

        t_327[k] = f_1 * sos0_54[k]
                   - f_2 * sos1_54[k]
                   + f_3 * pc_y[k] * sop_163[k];

        t_328[k] = f_3 * pc_y[k] * sop_164[k];

        t_329[k] = f_7 * snp_134[k]
                   + f_1 * sos0_54[k]
                   - f_2 * sos1_54[k]
                   + f_3 * pc_z[k] * sop_164[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, pb_x, pc_x, snd0_330, snd0_333, snp_165, \
                         snp_166, snp_167, snd1_330, snd1_333, sop_166, \
                         sop_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = pb_x[k] * snd0_330[k]
                   + f_8 * snp_165[k]
                   - f_4 * pc_x[k] * snd1_330[k];

        t_331[k] = f_6 * snp_166[k]
                   + f_3 * pc_x[k] * sop_166[k];

        t_332[k] = f_6 * snp_167[k]
                   + f_3 * pc_x[k] * sop_167[k];

        t_333[k] = pb_x[k] * snd0_333[k]
                   - f_4 * pc_x[k] * snd1_333[k];
    }

#pragma omp simd aligned(t_334, t_335, t_336, pb_x, pb_z, pc_x, pc_y, pc_z, snd0_270, \
                         snd0_335, snp_137, snd1_270, snd1_335, \
                         sop_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_334[k] = f_5 * snp_137[k]
                   + f_3 * pc_y[k] * sop_167[k];

        t_335[k] = pb_x[k] * snd0_335[k]
                   - f_4 * pc_x[k] * snd1_335[k];

        t_336[k] = pb_z[k] * snd0_270[k]
                   - f_4 * pc_z[k] * snd1_270[k];
    }

#pragma omp simd aligned(t_337, t_338, t_339, t_340, pb_x, pc_x, pc_y, snd0_339, snp_140, \
                         snp_169, snp_170, snd1_339, sop_169, sop_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_337[k] = f_6 * snp_169[k]
                   + f_3 * pc_x[k] * sop_169[k];

        t_338[k] = f_6 * snp_170[k]
                   + f_3 * pc_x[k] * sop_170[k];

        t_339[k] = pb_x[k] * snd0_339[k]
                   - f_4 * pc_x[k] * snd1_339[k];

        t_340[k] = f_7 * snp_140[k]
                   + f_3 * pc_y[k] * sop_170[k];
    }

#pragma omp simd aligned(t_341, t_342, t_343, t_344, pb_x, pc_x, snd0_341, snd0_342, snp_171, \
                         snp_172, snp_173, snd1_341, snd1_342, sop_172, \
                         sop_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_341[k] = pb_x[k] * snd0_341[k]
                   - f_4 * pc_x[k] * snd1_341[k];

        t_342[k] = pb_x[k] * snd0_342[k]
                   + f_8 * snp_171[k]
                   - f_4 * pc_x[k] * snd1_342[k];

        t_343[k] = f_6 * snp_172[k]
                   + f_3 * pc_x[k] * sop_172[k];

        t_344[k] = f_6 * snp_173[k]
                   + f_3 * pc_x[k] * sop_173[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, pb_x, pc_x, pc_y, snd0_345, snd0_347, \
                         snd0_348, snp_143, snp_174, snd1_345, snd1_347, snd1_348, \
                         sop_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = pb_x[k] * snd0_345[k]
                   - f_4 * pc_x[k] * snd1_345[k];

        t_346[k] = f_9 * snp_143[k]
                   + f_3 * pc_y[k] * sop_173[k];

        t_347[k] = pb_x[k] * snd0_347[k]
                   - f_4 * pc_x[k] * snd1_347[k];

        t_348[k] = pb_x[k] * snd0_348[k]
                   + f_8 * snp_174[k]
                   - f_4 * pc_x[k] * snd1_348[k];
    }

#pragma omp simd aligned(t_349, t_350, t_351, t_352, pb_x, pc_x, pc_y, snd0_351, snp_146, \
                         snp_175, snp_176, snd1_351, sop_175, sop_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_349[k] = f_6 * snp_175[k]
                   + f_3 * pc_x[k] * sop_175[k];

        t_350[k] = f_6 * snp_176[k]
                   + f_3 * pc_x[k] * sop_176[k];

        t_351[k] = pb_x[k] * snd0_351[k]
                   - f_4 * pc_x[k] * snd1_351[k];

        t_352[k] = f_11 * snp_146[k]
                   + f_3 * pc_y[k] * sop_176[k];
    }
}

static auto
compute_prim_sod_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t snd0,
                                                          const size_t snp, const size_t snd1,
                                                          const size_t sos0, const size_t sos1,
                                                          const size_t sop, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.5 / q;
    const auto f_1 = 0.5 / gamma;
    const auto f_2 = 0.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 5.0 / q;
    const auto f_6 = 0.5 / q;
    const auto f_7 = 4.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 4.0 / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 3.5 / q;
    const auto f_12 = 2.0 / q;
    const auto f_13 = 3.0 / q;
    const auto f_14 = 2.5 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *snd0_324 = buffer.data(snd0 + 324);
    const auto *snd0_330 = buffer.data(snd0 + 330);
    const auto *snd0_333 = buffer.data(snd0 + 333);
    const auto *snd0_353 = buffer.data(snd0 + 353);
    const auto *snd0_354 = buffer.data(snd0 + 354);
    const auto *snd0_357 = buffer.data(snd0 + 357);
    const auto *snd0_359 = buffer.data(snd0 + 359);
    const auto *snd0_360 = buffer.data(snd0 + 360);
    const auto *snd0_363 = buffer.data(snd0 + 363);
    const auto *snd0_365 = buffer.data(snd0 + 365);
    const auto *snd0_366 = buffer.data(snd0 + 366);
    const auto *snd0_369 = buffer.data(snd0 + 369);
    const auto *snd0_371 = buffer.data(snd0 + 371);
    const auto *snd0_372 = buffer.data(snd0 + 372);
    const auto *snd0_375 = buffer.data(snd0 + 375);
    const auto *snd0_377 = buffer.data(snd0 + 377);
    const auto *snd0_378 = buffer.data(snd0 + 378);
    const auto *snd0_381 = buffer.data(snd0 + 381);
    const auto *snd0_383 = buffer.data(snd0 + 383);
    const auto *snd0_387 = buffer.data(snd0 + 387);
    const auto *snd0_389 = buffer.data(snd0 + 389);
    const auto *snd0_390 = buffer.data(snd0 + 390);
    const auto *snd0_393 = buffer.data(snd0 + 393);
    const auto *snd0_395 = buffer.data(snd0 + 395);

    const auto *snp_149 = buffer.data(snp + 149);
    const auto *snp_152 = buffer.data(snp + 152);
    const auto *snp_155 = buffer.data(snp + 155);
    const auto *snp_158 = buffer.data(snp + 158);
    const auto *snp_161 = buffer.data(snp + 161);
    const auto *snp_164 = buffer.data(snp + 164);
    const auto *snp_166 = buffer.data(snp + 166);
    const auto *snp_167 = buffer.data(snp + 167);
    const auto *snp_170 = buffer.data(snp + 170);
    const auto *snp_172 = buffer.data(snp + 172);
    const auto *snp_173 = buffer.data(snp + 173);
    const auto *snp_175 = buffer.data(snp + 175);
    const auto *snp_176 = buffer.data(snp + 176);
    const auto *snp_177 = buffer.data(snp + 177);
    const auto *snp_178 = buffer.data(snp + 178);
    const auto *snp_179 = buffer.data(snp + 179);
    const auto *snp_180 = buffer.data(snp + 180);
    const auto *snp_181 = buffer.data(snp + 181);
    const auto *snp_182 = buffer.data(snp + 182);
    const auto *snp_183 = buffer.data(snp + 183);
    const auto *snp_184 = buffer.data(snp + 184);
    const auto *snp_185 = buffer.data(snp + 185);
    const auto *snp_186 = buffer.data(snp + 186);
    const auto *snp_187 = buffer.data(snp + 187);
    const auto *snp_188 = buffer.data(snp + 188);
    const auto *snp_189 = buffer.data(snp + 189);
    const auto *snp_190 = buffer.data(snp + 190);
    const auto *snp_191 = buffer.data(snp + 191);
    const auto *snp_193 = buffer.data(snp + 193);
    const auto *snp_194 = buffer.data(snp + 194);
    const auto *snp_195 = buffer.data(snp + 195);
    const auto *snp_196 = buffer.data(snp + 196);
    const auto *snp_197 = buffer.data(snp + 197);

    const auto *snd1_324 = buffer.data(snd1 + 324);
    const auto *snd1_330 = buffer.data(snd1 + 330);
    const auto *snd1_333 = buffer.data(snd1 + 333);
    const auto *snd1_353 = buffer.data(snd1 + 353);
    const auto *snd1_354 = buffer.data(snd1 + 354);
    const auto *snd1_357 = buffer.data(snd1 + 357);
    const auto *snd1_359 = buffer.data(snd1 + 359);
    const auto *snd1_360 = buffer.data(snd1 + 360);
    const auto *snd1_363 = buffer.data(snd1 + 363);
    const auto *snd1_365 = buffer.data(snd1 + 365);
    const auto *snd1_366 = buffer.data(snd1 + 366);
    const auto *snd1_369 = buffer.data(snd1 + 369);
    const auto *snd1_371 = buffer.data(snd1 + 371);
    const auto *snd1_372 = buffer.data(snd1 + 372);
    const auto *snd1_375 = buffer.data(snd1 + 375);
    const auto *snd1_377 = buffer.data(snd1 + 377);
    const auto *snd1_378 = buffer.data(snd1 + 378);
    const auto *snd1_381 = buffer.data(snd1 + 381);
    const auto *snd1_383 = buffer.data(snd1 + 383);
    const auto *snd1_387 = buffer.data(snd1 + 387);
    const auto *snd1_389 = buffer.data(snd1 + 389);
    const auto *snd1_390 = buffer.data(snd1 + 390);
    const auto *snd1_393 = buffer.data(snd1 + 393);
    const auto *snd1_395 = buffer.data(snd1 + 395);

    const auto *sos0_66 = buffer.data(sos0 + 66);
    const auto *sos0_67 = buffer.data(sos0 + 67);
    const auto *sos0_68 = buffer.data(sos0 + 68);
    const auto *sos0_69 = buffer.data(sos0 + 69);
    const auto *sos0_70 = buffer.data(sos0 + 70);
    const auto *sos0_71 = buffer.data(sos0 + 71);
    const auto *sos0_72 = buffer.data(sos0 + 72);
    const auto *sos0_73 = buffer.data(sos0 + 73);
    const auto *sos0_74 = buffer.data(sos0 + 74);
    const auto *sos0_75 = buffer.data(sos0 + 75);
    const auto *sos0_77 = buffer.data(sos0 + 77);

    const auto *sos1_66 = buffer.data(sos1 + 66);
    const auto *sos1_67 = buffer.data(sos1 + 67);
    const auto *sos1_68 = buffer.data(sos1 + 68);
    const auto *sos1_69 = buffer.data(sos1 + 69);
    const auto *sos1_70 = buffer.data(sos1 + 70);
    const auto *sos1_71 = buffer.data(sos1 + 71);
    const auto *sos1_72 = buffer.data(sos1 + 72);
    const auto *sos1_73 = buffer.data(sos1 + 73);
    const auto *sos1_74 = buffer.data(sos1 + 74);
    const auto *sos1_75 = buffer.data(sos1 + 75);
    const auto *sos1_77 = buffer.data(sos1 + 77);

    const auto *sop_178 = buffer.data(sop + 178);
    const auto *sop_179 = buffer.data(sop + 179);
    const auto *sop_181 = buffer.data(sop + 181);
    const auto *sop_182 = buffer.data(sop + 182);
    const auto *sop_184 = buffer.data(sop + 184);
    const auto *sop_185 = buffer.data(sop + 185);
    const auto *sop_187 = buffer.data(sop + 187);
    const auto *sop_188 = buffer.data(sop + 188);
    const auto *sop_190 = buffer.data(sop + 190);
    const auto *sop_191 = buffer.data(sop + 191);
    const auto *sop_193 = buffer.data(sop + 193);
    const auto *sop_194 = buffer.data(sop + 194);
    const auto *sop_196 = buffer.data(sop + 196);
    const auto *sop_197 = buffer.data(sop + 197);
    const auto *sop_198 = buffer.data(sop + 198);
    const auto *sop_199 = buffer.data(sop + 199);
    const auto *sop_200 = buffer.data(sop + 200);
    const auto *sop_202 = buffer.data(sop + 202);
    const auto *sop_203 = buffer.data(sop + 203);
    const auto *sop_204 = buffer.data(sop + 204);
    const auto *sop_205 = buffer.data(sop + 205);
    const auto *sop_206 = buffer.data(sop + 206);
    const auto *sop_207 = buffer.data(sop + 207);
    const auto *sop_208 = buffer.data(sop + 208);
    const auto *sop_209 = buffer.data(sop + 209);
    const auto *sop_210 = buffer.data(sop + 210);
    const auto *sop_211 = buffer.data(sop + 211);
    const auto *sop_212 = buffer.data(sop + 212);
    const auto *sop_213 = buffer.data(sop + 213);
    const auto *sop_214 = buffer.data(sop + 214);
    const auto *sop_215 = buffer.data(sop + 215);
    const auto *sop_216 = buffer.data(sop + 216);
    const auto *sop_217 = buffer.data(sop + 217);
    const auto *sop_218 = buffer.data(sop + 218);
    const auto *sop_219 = buffer.data(sop + 219);
    const auto *sop_220 = buffer.data(sop + 220);
    const auto *sop_221 = buffer.data(sop + 221);
    const auto *sop_222 = buffer.data(sop + 222);
    const auto *sop_223 = buffer.data(sop + 223);
    const auto *sop_224 = buffer.data(sop + 224);
    const auto *sop_225 = buffer.data(sop + 225);
    const auto *sop_226 = buffer.data(sop + 226);
    const auto *sop_227 = buffer.data(sop + 227);
    const auto *sop_229 = buffer.data(sop + 229);
    const auto *sop_230 = buffer.data(sop + 230);
    const auto *sop_231 = buffer.data(sop + 231);
    const auto *sop_232 = buffer.data(sop + 232);
    const auto *sop_233 = buffer.data(sop + 233);

#pragma omp simd aligned(t_353, t_354, t_355, t_356, pb_x, pc_x, snd0_353, snd0_354, snp_177, \
                         snp_178, snp_179, snd1_353, snd1_354, sop_178, \
                         sop_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_353[k] = pb_x[k] * snd0_353[k]
                   - f_4 * pc_x[k] * snd1_353[k];

        t_354[k] = pb_x[k] * snd0_354[k]
                   + f_8 * snp_177[k]
                   - f_4 * pc_x[k] * snd1_354[k];

        t_355[k] = f_6 * snp_178[k]
                   + f_3 * pc_x[k] * sop_178[k];

        t_356[k] = f_6 * snp_179[k]
                   + f_3 * pc_x[k] * sop_179[k];
    }

#pragma omp simd aligned(t_357, t_358, t_359, t_360, pb_x, pc_x, pc_y, snd0_357, snd0_359, \
                         snd0_360, snp_149, snp_180, snd1_357, snd1_359, snd1_360, \
                         sop_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_357[k] = pb_x[k] * snd0_357[k]
                   - f_4 * pc_x[k] * snd1_357[k];

        t_358[k] = f_13 * snp_149[k]
                   + f_3 * pc_y[k] * sop_179[k];

        t_359[k] = pb_x[k] * snd0_359[k]
                   - f_4 * pc_x[k] * snd1_359[k];

        t_360[k] = pb_x[k] * snd0_360[k]
                   + f_8 * snp_180[k]
                   - f_4 * pc_x[k] * snd1_360[k];
    }

#pragma omp simd aligned(t_361, t_362, t_363, t_364, pb_x, pc_x, pc_y, snd0_363, snp_152, \
                         snp_181, snp_182, snd1_363, sop_181, sop_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_361[k] = f_6 * snp_181[k]
                   + f_3 * pc_x[k] * sop_181[k];

        t_362[k] = f_6 * snp_182[k]
                   + f_3 * pc_x[k] * sop_182[k];

        t_363[k] = pb_x[k] * snd0_363[k]
                   - f_4 * pc_x[k] * snd1_363[k];

        t_364[k] = f_14 * snp_152[k]
                   + f_3 * pc_y[k] * sop_182[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, pb_x, pc_x, snd0_365, snd0_366, snp_183, \
                         snp_184, snp_185, snd1_365, snd1_366, sop_184, \
                         sop_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = pb_x[k] * snd0_365[k]
                   - f_4 * pc_x[k] * snd1_365[k];

        t_366[k] = pb_x[k] * snd0_366[k]
                   + f_8 * snp_183[k]
                   - f_4 * pc_x[k] * snd1_366[k];

        t_367[k] = f_6 * snp_184[k]
                   + f_3 * pc_x[k] * sop_184[k];

        t_368[k] = f_6 * snp_185[k]
                   + f_3 * pc_x[k] * sop_185[k];
    }

#pragma omp simd aligned(t_369, t_370, t_371, t_372, pb_x, pc_x, pc_y, snd0_369, snd0_371, \
                         snd0_372, snp_155, snp_186, snd1_369, snd1_371, snd1_372, \
                         sop_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_369[k] = pb_x[k] * snd0_369[k]
                   - f_4 * pc_x[k] * snd1_369[k];

        t_370[k] = f_12 * snp_155[k]
                   + f_3 * pc_y[k] * sop_185[k];

        t_371[k] = pb_x[k] * snd0_371[k]
                   - f_4 * pc_x[k] * snd1_371[k];

        t_372[k] = pb_x[k] * snd0_372[k]
                   + f_8 * snp_186[k]
                   - f_4 * pc_x[k] * snd1_372[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, t_376, pb_x, pc_x, pc_y, snd0_375, snp_158, \
                         snp_187, snp_188, snd1_375, sop_187, sop_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = f_6 * snp_187[k]
                   + f_3 * pc_x[k] * sop_187[k];

        t_374[k] = f_6 * snp_188[k]
                   + f_3 * pc_x[k] * sop_188[k];

        t_375[k] = pb_x[k] * snd0_375[k]
                   - f_4 * pc_x[k] * snd1_375[k];

        t_376[k] = f_10 * snp_158[k]
                   + f_3 * pc_y[k] * sop_188[k];
    }

#pragma omp simd aligned(t_377, t_378, t_379, t_380, pb_x, pc_x, snd0_377, snd0_378, snp_189, \
                         snp_190, snp_191, snd1_377, snd1_378, sop_190, \
                         sop_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_377[k] = pb_x[k] * snd0_377[k]
                   - f_4 * pc_x[k] * snd1_377[k];

        t_378[k] = pb_x[k] * snd0_378[k]
                   + f_8 * snp_189[k]
                   - f_4 * pc_x[k] * snd1_378[k];

        t_379[k] = f_6 * snp_190[k]
                   + f_3 * pc_x[k] * sop_190[k];

        t_380[k] = f_6 * snp_191[k]
                   + f_3 * pc_x[k] * sop_191[k];
    }

#pragma omp simd aligned(t_381, t_382, t_383, t_384, pb_x, pb_y, pc_x, pc_y, snd0_324, \
                         snd0_381, snd0_383, snp_161, snd1_324, snd1_381, snd1_383, \
                         sop_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_381[k] = pb_x[k] * snd0_381[k]
                   - f_4 * pc_x[k] * snd1_381[k];

        t_382[k] = f_8 * snp_161[k]
                   + f_3 * pc_y[k] * sop_191[k];

        t_383[k] = pb_x[k] * snd0_383[k]
                   - f_4 * pc_x[k] * snd1_383[k];

        t_384[k] = pb_y[k] * snd0_324[k]
                   - f_4 * pc_y[k] * snd1_324[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, pb_x, pc_x, pc_y, snd0_387, snp_164, \
                         snp_193, snp_194, snd1_387, sop_193, sop_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = f_6 * snp_193[k]
                   + f_3 * pc_x[k] * sop_193[k];

        t_386[k] = f_6 * snp_194[k]
                   + f_3 * pc_x[k] * sop_194[k];

        t_387[k] = pb_x[k] * snd0_387[k]
                   - f_4 * pc_x[k] * snd1_387[k];

        t_388[k] = f_6 * snp_164[k]
                   + f_3 * pc_y[k] * sop_194[k];
    }

#pragma omp simd aligned(t_389, t_390, t_391, t_392, pb_x, pc_x, snd0_389, snd0_390, snp_195, \
                         snp_196, snp_197, snd1_389, snd1_390, sop_196, \
                         sop_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_389[k] = pb_x[k] * snd0_389[k]
                   - f_4 * pc_x[k] * snd1_389[k];

        t_390[k] = pb_x[k] * snd0_390[k]
                   + f_8 * snp_195[k]
                   - f_4 * pc_x[k] * snd1_390[k];

        t_391[k] = f_6 * snp_196[k]
                   + f_3 * pc_x[k] * sop_196[k];

        t_392[k] = f_6 * snp_197[k]
                   + f_3 * pc_x[k] * sop_197[k];
    }

#pragma omp simd aligned(t_393, t_394, t_395, t_396, pb_x, pc_x, pc_y, snd0_393, snd0_395, \
                         snd1_393, snd1_395, sos0_66, sos1_66, sop_197, \
                         sop_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_393[k] = pb_x[k] * snd0_393[k]
                   - f_4 * pc_x[k] * snd1_393[k];

        t_394[k] = f_3 * pc_y[k] * sop_197[k];

        t_395[k] = pb_x[k] * snd0_395[k]
                   - f_4 * pc_x[k] * snd1_395[k];

        t_396[k] = f_1 * sos0_66[k]
                   - f_2 * sos1_66[k]
                   + f_3 * pc_x[k] * sop_198[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, t_400, t_401, pc_x, pc_y, pc_z, snp_166, \
                         snp_167, sos0_66, sos1_66, sop_199, sop_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = f_3 * pc_x[k] * sop_199[k];

        t_398[k] = f_3 * pc_x[k] * sop_200[k];

        t_399[k] = f_0 * snp_166[k]
                   + f_1 * sos0_66[k]
                   - f_2 * sos1_66[k]
                   + f_3 * pc_y[k] * sop_199[k];

        t_400[k] = f_0 * snp_167[k]
                   + f_3 * pc_y[k] * sop_200[k];

        t_401[k] = f_1 * sos0_66[k]
                   - f_2 * sos1_66[k]
                   + f_3 * pc_z[k] * sop_200[k];
    }

#pragma omp simd aligned(t_402, t_403, t_404, t_405, t_406, pb_z, pc_x, pc_y, pc_z, snd0_330, \
                         snd0_333, snp_170, snd1_330, snd1_333, sop_202, \
                         sop_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_402[k] = pb_z[k] * snd0_330[k]
                   - f_4 * pc_z[k] * snd1_330[k];

        t_403[k] = f_3 * pc_x[k] * sop_202[k];

        t_404[k] = f_3 * pc_x[k] * sop_203[k];

        t_405[k] = pb_z[k] * snd0_333[k]
                   - f_4 * pc_z[k] * snd1_333[k];

        t_406[k] = f_5 * snp_170[k]
                   + f_3 * pc_y[k] * sop_203[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, t_410, pc_x, pc_z, snp_167, sos0_67, sos0_68, \
                         sos1_67, sos1_68, sop_203, sop_204, sop_205, \
                         sop_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = f_6 * snp_167[k]
                   + f_1 * sos0_67[k]
                   - f_2 * sos1_67[k]
                   + f_3 * pc_z[k] * sop_203[k];

        t_408[k] = f_1 * sos0_68[k]
                   - f_2 * sos1_68[k]
                   + f_3 * pc_x[k] * sop_204[k];

        t_409[k] = f_3 * pc_x[k] * sop_205[k];

        t_410[k] = f_3 * pc_x[k] * sop_206[k];
    }

#pragma omp simd aligned(t_411, t_412, t_413, pc_y, pc_z, snp_170, snp_172, snp_173, sos0_68, \
                         sos1_68, sop_205, sop_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_411[k] = f_7 * snp_172[k]
                   + f_1 * sos0_68[k]
                   - f_2 * sos1_68[k]
                   + f_3 * pc_y[k] * sop_205[k];

        t_412[k] = f_7 * snp_173[k]
                   + f_3 * pc_y[k] * sop_206[k];

        t_413[k] = f_8 * snp_170[k]
                   + f_1 * sos0_68[k]
                   - f_2 * sos1_68[k]
                   + f_3 * pc_z[k] * sop_206[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, t_417, t_418, pc_x, pc_y, snp_175, snp_176, \
                         sos0_69, sos1_69, sop_207, sop_208, sop_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_1 * sos0_69[k]
                   - f_2 * sos1_69[k]
                   + f_3 * pc_x[k] * sop_207[k];

        t_415[k] = f_3 * pc_x[k] * sop_208[k];

        t_416[k] = f_3 * pc_x[k] * sop_209[k];

        t_417[k] = f_9 * snp_175[k]
                   + f_1 * sos0_69[k]
                   - f_2 * sos1_69[k]
                   + f_3 * pc_y[k] * sop_208[k];

        t_418[k] = f_9 * snp_176[k]
                   + f_3 * pc_y[k] * sop_209[k];
    }

#pragma omp simd aligned(t_419, t_420, t_421, t_422, pc_x, pc_z, snp_173, sos0_69, sos0_70, \
                         sos1_69, sos1_70, sop_209, sop_210, sop_211, \
                         sop_212 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_419[k] = f_10 * snp_173[k]
                   + f_1 * sos0_69[k]
                   - f_2 * sos1_69[k]
                   + f_3 * pc_z[k] * sop_209[k];

        t_420[k] = f_1 * sos0_70[k]
                   - f_2 * sos1_70[k]
                   + f_3 * pc_x[k] * sop_210[k];

        t_421[k] = f_3 * pc_x[k] * sop_211[k];

        t_422[k] = f_3 * pc_x[k] * sop_212[k];
    }

#pragma omp simd aligned(t_423, t_424, t_425, pc_y, pc_z, snp_176, snp_178, snp_179, sos0_70, \
                         sos1_70, sop_211, sop_212 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_423[k] = f_11 * snp_178[k]
                   + f_1 * sos0_70[k]
                   - f_2 * sos1_70[k]
                   + f_3 * pc_y[k] * sop_211[k];

        t_424[k] = f_11 * snp_179[k]
                   + f_3 * pc_y[k] * sop_212[k];

        t_425[k] = f_12 * snp_176[k]
                   + f_1 * sos0_70[k]
                   - f_2 * sos1_70[k]
                   + f_3 * pc_z[k] * sop_212[k];
    }

#pragma omp simd aligned(t_426, t_427, t_428, t_429, t_430, pc_x, pc_y, snp_181, snp_182, \
                         sos0_71, sos1_71, sop_213, sop_214, sop_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_426[k] = f_1 * sos0_71[k]
                   - f_2 * sos1_71[k]
                   + f_3 * pc_x[k] * sop_213[k];

        t_427[k] = f_3 * pc_x[k] * sop_214[k];

        t_428[k] = f_3 * pc_x[k] * sop_215[k];

        t_429[k] = f_13 * snp_181[k]
                   + f_1 * sos0_71[k]
                   - f_2 * sos1_71[k]
                   + f_3 * pc_y[k] * sop_214[k];

        t_430[k] = f_13 * snp_182[k]
                   + f_3 * pc_y[k] * sop_215[k];
    }

#pragma omp simd aligned(t_431, t_432, t_433, t_434, pc_x, pc_z, snp_179, sos0_71, sos0_72, \
                         sos1_71, sos1_72, sop_215, sop_216, sop_217, \
                         sop_218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_431[k] = f_14 * snp_179[k]
                   + f_1 * sos0_71[k]
                   - f_2 * sos1_71[k]
                   + f_3 * pc_z[k] * sop_215[k];

        t_432[k] = f_1 * sos0_72[k]
                   - f_2 * sos1_72[k]
                   + f_3 * pc_x[k] * sop_216[k];

        t_433[k] = f_3 * pc_x[k] * sop_217[k];

        t_434[k] = f_3 * pc_x[k] * sop_218[k];
    }

#pragma omp simd aligned(t_435, t_436, t_437, pc_y, pc_z, snp_182, snp_184, snp_185, sos0_72, \
                         sos1_72, sop_217, sop_218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_435[k] = f_14 * snp_184[k]
                   + f_1 * sos0_72[k]
                   - f_2 * sos1_72[k]
                   + f_3 * pc_y[k] * sop_217[k];

        t_436[k] = f_14 * snp_185[k]
                   + f_3 * pc_y[k] * sop_218[k];

        t_437[k] = f_13 * snp_182[k]
                   + f_1 * sos0_72[k]
                   - f_2 * sos1_72[k]
                   + f_3 * pc_z[k] * sop_218[k];
    }

#pragma omp simd aligned(t_438, t_439, t_440, t_441, t_442, pc_x, pc_y, snp_187, snp_188, \
                         sos0_73, sos1_73, sop_219, sop_220, sop_221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_438[k] = f_1 * sos0_73[k]
                   - f_2 * sos1_73[k]
                   + f_3 * pc_x[k] * sop_219[k];

        t_439[k] = f_3 * pc_x[k] * sop_220[k];

        t_440[k] = f_3 * pc_x[k] * sop_221[k];

        t_441[k] = f_12 * snp_187[k]
                   + f_1 * sos0_73[k]
                   - f_2 * sos1_73[k]
                   + f_3 * pc_y[k] * sop_220[k];

        t_442[k] = f_12 * snp_188[k]
                   + f_3 * pc_y[k] * sop_221[k];
    }

#pragma omp simd aligned(t_443, t_444, t_445, t_446, pc_x, pc_z, snp_185, sos0_73, sos0_74, \
                         sos1_73, sos1_74, sop_221, sop_222, sop_223, \
                         sop_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_443[k] = f_11 * snp_185[k]
                   + f_1 * sos0_73[k]
                   - f_2 * sos1_73[k]
                   + f_3 * pc_z[k] * sop_221[k];

        t_444[k] = f_1 * sos0_74[k]
                   - f_2 * sos1_74[k]
                   + f_3 * pc_x[k] * sop_222[k];

        t_445[k] = f_3 * pc_x[k] * sop_223[k];

        t_446[k] = f_3 * pc_x[k] * sop_224[k];
    }

#pragma omp simd aligned(t_447, t_448, t_449, pc_y, pc_z, snp_188, snp_190, snp_191, sos0_74, \
                         sos1_74, sop_223, sop_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_447[k] = f_10 * snp_190[k]
                   + f_1 * sos0_74[k]
                   - f_2 * sos1_74[k]
                   + f_3 * pc_y[k] * sop_223[k];

        t_448[k] = f_10 * snp_191[k]
                   + f_3 * pc_y[k] * sop_224[k];

        t_449[k] = f_9 * snp_188[k]
                   + f_1 * sos0_74[k]
                   - f_2 * sos1_74[k]
                   + f_3 * pc_z[k] * sop_224[k];
    }

#pragma omp simd aligned(t_450, t_451, t_452, t_453, t_454, pc_x, pc_y, snp_193, snp_194, \
                         sos0_75, sos1_75, sop_225, sop_226, sop_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_450[k] = f_1 * sos0_75[k]
                   - f_2 * sos1_75[k]
                   + f_3 * pc_x[k] * sop_225[k];

        t_451[k] = f_3 * pc_x[k] * sop_226[k];

        t_452[k] = f_3 * pc_x[k] * sop_227[k];

        t_453[k] = f_8 * snp_193[k]
                   + f_1 * sos0_75[k]
                   - f_2 * sos1_75[k]
                   + f_3 * pc_y[k] * sop_226[k];

        t_454[k] = f_8 * snp_194[k]
                   + f_3 * pc_y[k] * sop_227[k];
    }

#pragma omp simd aligned(t_455, t_456, t_457, t_458, pb_y, pc_x, pc_y, pc_z, snd0_390, \
                         snp_191, snd1_390, sos0_75, sos1_75, sop_227, sop_229, \
                         sop_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_455[k] = f_7 * snp_191[k]
                   + f_1 * sos0_75[k]
                   - f_2 * sos1_75[k]
                   + f_3 * pc_z[k] * sop_227[k];

        t_456[k] = pb_y[k] * snd0_390[k]
                   - f_4 * pc_y[k] * snd1_390[k];

        t_457[k] = f_3 * pc_x[k] * sop_229[k];

        t_458[k] = f_3 * pc_x[k] * sop_230[k];
    }

#pragma omp simd aligned(t_459, t_460, t_461, pb_y, pc_y, snd0_393, snd0_395, snp_196, \
                         snp_197, snd1_393, snd1_395, sop_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_459[k] = pb_y[k] * snd0_393[k]
                   + f_8 * snp_196[k]
                   - f_4 * pc_y[k] * snd1_393[k];

        t_460[k] = f_6 * snp_197[k]
                   + f_3 * pc_y[k] * sop_230[k];

        t_461[k] = pb_y[k] * snd0_395[k]
                   - f_4 * pc_y[k] * snd1_395[k];
    }

#pragma omp simd aligned(t_462, t_463, t_464, t_465, t_466, t_467, pc_x, pc_y, pc_z, snp_197, \
                         sos0_77, sos1_77, sop_231, sop_232, sop_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_462[k] = f_1 * sos0_77[k]
                   - f_2 * sos1_77[k]
                   + f_3 * pc_x[k] * sop_231[k];

        t_463[k] = f_3 * pc_x[k] * sop_232[k];

        t_464[k] = f_3 * pc_x[k] * sop_233[k];

        t_465[k] = f_1 * sos0_77[k]
                   - f_2 * sos1_77[k]
                   + f_3 * pc_y[k] * sop_232[k];

        t_466[k] = f_3 * pc_y[k] * sop_233[k];

        t_467[k] = f_0 * snp_197[k]
                   + f_1 * sos0_77[k]
                   - f_2 * sos1_77[k]
                   + f_3 * pc_z[k] * sop_233[k];
    }
}

auto
compute_prim_sod_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t snd0, const size_t snp,
                                                   const size_t snd1, const size_t sos0,
                                                   const size_t sos1, const size_t sop,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_sod_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, snd0, snp,
                                                              snd1, sos0, sos1, sop, ncols,
                                                              gamma, p, q);

    compute_prim_sod_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, snd0, snp,
                                                              snd1, sos0, sos1, sop, ncols,
                                                              gamma, p, q);

    compute_prim_sod_three_center_electron_repulsion_0_piece2(buffer, target, pb, pc, snd0, snp,
                                                              snd1, sos0, sos1, sop, ncols,
                                                              gamma, p, q);

    compute_prim_sod_three_center_electron_repulsion_0_piece3(buffer, target, pb, pc, snd0, snp,
                                                              snd1, sos0, sos1, sop, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
