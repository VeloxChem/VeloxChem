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


#include "SimdThreeCenterElectronRepulsionVrrRecSLD.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_sld_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t skd0,
                                                          const size_t skp, const size_t skd1,
                                                          const size_t sls0, const size_t sls1,
                                                          const size_t slp, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
    const auto f_1 = 0.5 / gamma;
    const auto f_2 = 0.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 3.5 / q;
    const auto f_6 = 0.5 / q;
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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *skd0_0 = buffer.data(skd0 + 0);
    const auto *skd0_3 = buffer.data(skd0 + 3);
    const auto *skd0_5 = buffer.data(skd0 + 5);
    const auto *skd0_9 = buffer.data(skd0 + 9);
    const auto *skd0_12 = buffer.data(skd0 + 12);
    const auto *skd0_17 = buffer.data(skd0 + 17);
    const auto *skd0_18 = buffer.data(skd0 + 18);
    const auto *skd0_21 = buffer.data(skd0 + 21);
    const auto *skd0_30 = buffer.data(skd0 + 30);
    const auto *skd0_35 = buffer.data(skd0 + 35);
    const auto *skd0_36 = buffer.data(skd0 + 36);
    const auto *skd0_39 = buffer.data(skd0 + 39);
    const auto *skd0_54 = buffer.data(skd0 + 54);
    const auto *skd0_59 = buffer.data(skd0 + 59);
    const auto *skd0_60 = buffer.data(skd0 + 60);
    const auto *skd0_63 = buffer.data(skd0 + 63);
    const auto *skd0_84 = buffer.data(skd0 + 84);
    const auto *skd0_89 = buffer.data(skd0 + 89);

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

    const auto *skd1_0 = buffer.data(skd1 + 0);
    const auto *skd1_3 = buffer.data(skd1 + 3);
    const auto *skd1_5 = buffer.data(skd1 + 5);
    const auto *skd1_9 = buffer.data(skd1 + 9);
    const auto *skd1_12 = buffer.data(skd1 + 12);
    const auto *skd1_17 = buffer.data(skd1 + 17);
    const auto *skd1_18 = buffer.data(skd1 + 18);
    const auto *skd1_21 = buffer.data(skd1 + 21);
    const auto *skd1_30 = buffer.data(skd1 + 30);
    const auto *skd1_35 = buffer.data(skd1 + 35);
    const auto *skd1_36 = buffer.data(skd1 + 36);
    const auto *skd1_39 = buffer.data(skd1 + 39);
    const auto *skd1_54 = buffer.data(skd1 + 54);
    const auto *skd1_59 = buffer.data(skd1 + 59);
    const auto *skd1_60 = buffer.data(skd1 + 60);
    const auto *skd1_63 = buffer.data(skd1 + 63);
    const auto *skd1_84 = buffer.data(skd1 + 84);
    const auto *skd1_89 = buffer.data(skd1 + 89);

    const auto *sls0_0 = buffer.data(sls0 + 0);
    const auto *sls0_1 = buffer.data(sls0 + 1);
    const auto *sls0_2 = buffer.data(sls0 + 2);
    const auto *sls0_3 = buffer.data(sls0 + 3);
    const auto *sls0_5 = buffer.data(sls0 + 5);
    const auto *sls0_6 = buffer.data(sls0 + 6);
    const auto *sls0_7 = buffer.data(sls0 + 7);
    const auto *sls0_8 = buffer.data(sls0 + 8);
    const auto *sls0_9 = buffer.data(sls0 + 9);
    const auto *sls0_10 = buffer.data(sls0 + 10);
    const auto *sls0_11 = buffer.data(sls0 + 11);
    const auto *sls0_12 = buffer.data(sls0 + 12);
    const auto *sls0_13 = buffer.data(sls0 + 13);
    const auto *sls0_14 = buffer.data(sls0 + 14);
    const auto *sls0_15 = buffer.data(sls0 + 15);
    const auto *sls0_16 = buffer.data(sls0 + 16);
    const auto *sls0_17 = buffer.data(sls0 + 17);
    const auto *sls0_18 = buffer.data(sls0 + 18);
    const auto *sls0_19 = buffer.data(sls0 + 19);

    const auto *sls1_0 = buffer.data(sls1 + 0);
    const auto *sls1_1 = buffer.data(sls1 + 1);
    const auto *sls1_2 = buffer.data(sls1 + 2);
    const auto *sls1_3 = buffer.data(sls1 + 3);
    const auto *sls1_5 = buffer.data(sls1 + 5);
    const auto *sls1_6 = buffer.data(sls1 + 6);
    const auto *sls1_7 = buffer.data(sls1 + 7);
    const auto *sls1_8 = buffer.data(sls1 + 8);
    const auto *sls1_9 = buffer.data(sls1 + 9);
    const auto *sls1_10 = buffer.data(sls1 + 10);
    const auto *sls1_11 = buffer.data(sls1 + 11);
    const auto *sls1_12 = buffer.data(sls1 + 12);
    const auto *sls1_13 = buffer.data(sls1 + 13);
    const auto *sls1_14 = buffer.data(sls1 + 14);
    const auto *sls1_15 = buffer.data(sls1 + 15);
    const auto *sls1_16 = buffer.data(sls1 + 16);
    const auto *sls1_17 = buffer.data(sls1 + 17);
    const auto *sls1_18 = buffer.data(sls1 + 18);
    const auto *sls1_19 = buffer.data(sls1 + 19);

    const auto *slp_0 = buffer.data(slp + 0);
    const auto *slp_1 = buffer.data(slp + 1);
    const auto *slp_2 = buffer.data(slp + 2);
    const auto *slp_4 = buffer.data(slp + 4);
    const auto *slp_5 = buffer.data(slp + 5);
    const auto *slp_7 = buffer.data(slp + 7);
    const auto *slp_8 = buffer.data(slp + 8);
    const auto *slp_9 = buffer.data(slp + 9);
    const auto *slp_10 = buffer.data(slp + 10);
    const auto *slp_11 = buffer.data(slp + 11);
    const auto *slp_13 = buffer.data(slp + 13);
    const auto *slp_14 = buffer.data(slp + 14);
    const auto *slp_15 = buffer.data(slp + 15);
    const auto *slp_16 = buffer.data(slp + 16);
    const auto *slp_17 = buffer.data(slp + 17);
    const auto *slp_18 = buffer.data(slp + 18);
    const auto *slp_19 = buffer.data(slp + 19);
    const auto *slp_20 = buffer.data(slp + 20);
    const auto *slp_22 = buffer.data(slp + 22);
    const auto *slp_23 = buffer.data(slp + 23);
    const auto *slp_25 = buffer.data(slp + 25);
    const auto *slp_26 = buffer.data(slp + 26);
    const auto *slp_27 = buffer.data(slp + 27);
    const auto *slp_28 = buffer.data(slp + 28);
    const auto *slp_29 = buffer.data(slp + 29);
    const auto *slp_30 = buffer.data(slp + 30);
    const auto *slp_31 = buffer.data(slp + 31);
    const auto *slp_32 = buffer.data(slp + 32);
    const auto *slp_34 = buffer.data(slp + 34);
    const auto *slp_35 = buffer.data(slp + 35);
    const auto *slp_36 = buffer.data(slp + 36);
    const auto *slp_37 = buffer.data(slp + 37);
    const auto *slp_38 = buffer.data(slp + 38);
    const auto *slp_40 = buffer.data(slp + 40);
    const auto *slp_41 = buffer.data(slp + 41);
    const auto *slp_42 = buffer.data(slp + 42);
    const auto *slp_43 = buffer.data(slp + 43);
    const auto *slp_44 = buffer.data(slp + 44);
    const auto *slp_45 = buffer.data(slp + 45);
    const auto *slp_46 = buffer.data(slp + 46);
    const auto *slp_47 = buffer.data(slp + 47);
    const auto *slp_49 = buffer.data(slp + 49);
    const auto *slp_50 = buffer.data(slp + 50);
    const auto *slp_51 = buffer.data(slp + 51);
    const auto *slp_52 = buffer.data(slp + 52);
    const auto *slp_53 = buffer.data(slp + 53);
    const auto *slp_54 = buffer.data(slp + 54);
    const auto *slp_55 = buffer.data(slp + 55);
    const auto *slp_56 = buffer.data(slp + 56);
    const auto *slp_58 = buffer.data(slp + 58);
    const auto *slp_59 = buffer.data(slp + 59);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pc_x, pc_y, skp_0, skp_1, skp_2, sls0_0, \
                         sls1_0, slp_0, slp_1, slp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * skp_0[k]
                 + f_1 * sls0_0[k]
                 - f_2 * sls1_0[k]
                 + f_3 * pc_x[k] * slp_0[k];

        t_1[k] = f_0 * skp_1[k]
                 + f_3 * pc_x[k] * slp_1[k];

        t_2[k] = f_0 * skp_2[k]
                 + f_3 * pc_x[k] * slp_2[k];

        t_3[k] = f_1 * sls0_0[k]
                 - f_2 * sls1_0[k]
                 + f_3 * pc_y[k] * slp_1[k];

        t_4[k] = f_3 * pc_y[k] * slp_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pb_y, pc_x, pc_y, pc_z, skd0_0, skp_4, skd1_0, sls0_0, \
                         sls1_0, slp_2, slp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * sls0_0[k]
                 - f_2 * sls1_0[k]
                 + f_3 * pc_z[k] * slp_2[k];

        t_6[k] = pb_y[k] * skd0_0[k]
                 - f_4 * pc_y[k] * skd1_0[k];

        t_7[k] = f_5 * skp_4[k]
                 + f_3 * pc_x[k] * slp_4[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pb_y, pc_x, pc_y, skd0_5, skp_1, skp_2, skp_5, \
                         skd1_5, sls0_1, sls1_1, slp_4, slp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_5 * skp_5[k]
                 + f_3 * pc_x[k] * slp_5[k];

        t_9[k] = f_6 * skp_1[k]
                 + f_1 * sls0_1[k]
                 - f_2 * sls1_1[k]
                 + f_3 * pc_y[k] * slp_4[k];

        t_10[k] = f_6 * skp_2[k]
                  + f_3 * pc_y[k] * slp_5[k];

        t_11[k] = pb_y[k] * skd0_5[k]
                  - f_4 * pc_y[k] * skd1_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pb_z, pc_x, pc_z, skd0_0, skd0_3, skp_7, \
                         skp_8, skd1_0, skd1_3, slp_7, slp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pb_z[k] * skd0_0[k]
                  - f_4 * pc_z[k] * skd1_0[k];

        t_13[k] = f_5 * skp_7[k]
                  + f_3 * pc_x[k] * slp_7[k];

        t_14[k] = f_5 * skp_8[k]
                  + f_3 * pc_x[k] * slp_8[k];

        t_15[k] = pb_z[k] * skd0_3[k]
                  - f_4 * pc_z[k] * skd1_3[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pc_x, pc_y, pc_z, skp_2, skp_9, sls0_2, sls0_3, \
                         sls1_2, sls1_3, slp_8, slp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * pc_y[k] * slp_8[k];

        t_17[k] = f_6 * skp_2[k]
                  + f_1 * sls0_2[k]
                  - f_2 * sls1_2[k]
                  + f_3 * pc_z[k] * slp_8[k];

        t_18[k] = f_7 * skp_9[k]
                  + f_1 * sls0_3[k]
                  - f_2 * sls1_3[k]
                  + f_3 * pc_x[k] * slp_9[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pc_x, pc_y, pc_z, skp_4, skp_5, skp_10, \
                         skp_11, sls0_3, sls1_3, slp_10, slp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_7 * skp_10[k]
                  + f_3 * pc_x[k] * slp_10[k];

        t_20[k] = f_7 * skp_11[k]
                  + f_3 * pc_x[k] * slp_11[k];

        t_21[k] = f_8 * skp_4[k]
                  + f_1 * sls0_3[k]
                  - f_2 * sls1_3[k]
                  + f_3 * pc_y[k] * slp_10[k];

        t_22[k] = f_8 * skp_5[k]
                  + f_3 * pc_y[k] * slp_11[k];

        t_23[k] = f_1 * sls0_3[k]
                  - f_2 * sls1_3[k]
                  + f_3 * pc_z[k] * slp_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pb_y, pc_x, pc_y, skd0_12, skp_13, skp_14, skd1_12, \
                         slp_13, slp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pb_y[k] * skd0_12[k]
                  - f_4 * pc_y[k] * skd1_12[k];

        t_25[k] = f_7 * skp_13[k]
                  + f_3 * pc_x[k] * slp_13[k];

        t_26[k] = f_7 * skp_14[k]
                  + f_3 * pc_x[k] * slp_14[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pb_y, pb_z, pc_y, pc_z, skd0_9, skd0_17, skp_8, \
                         skd1_9, skd1_17, slp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = pb_z[k] * skd0_9[k]
                  - f_4 * pc_z[k] * skd1_9[k];

        t_28[k] = f_6 * skp_8[k]
                  + f_3 * pc_y[k] * slp_14[k];

        t_29[k] = pb_y[k] * skd0_17[k]
                  - f_4 * pc_y[k] * skd1_17[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pc_x, pc_y, skp_15, skp_16, skp_17, \
                         sls0_5, sls1_5, slp_15, slp_16, slp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_7 * skp_15[k]
                  + f_1 * sls0_5[k]
                  - f_2 * sls1_5[k]
                  + f_3 * pc_x[k] * slp_15[k];

        t_31[k] = f_7 * skp_16[k]
                  + f_3 * pc_x[k] * slp_16[k];

        t_32[k] = f_7 * skp_17[k]
                  + f_3 * pc_x[k] * slp_17[k];

        t_33[k] = f_1 * sls0_5[k]
                  - f_2 * sls1_5[k]
                  + f_3 * pc_y[k] * slp_16[k];

        t_34[k] = f_3 * pc_y[k] * slp_17[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pc_x, pc_z, skp_8, skp_18, skp_19, sls0_5, sls0_6, \
                         sls1_5, sls1_6, slp_17, slp_18, slp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_8 * skp_8[k]
                  + f_1 * sls0_5[k]
                  - f_2 * sls1_5[k]
                  + f_3 * pc_z[k] * slp_17[k];

        t_36[k] = f_9 * skp_18[k]
                  + f_1 * sls0_6[k]
                  - f_2 * sls1_6[k]
                  + f_3 * pc_x[k] * slp_18[k];

        t_37[k] = f_9 * skp_19[k]
                  + f_3 * pc_x[k] * slp_19[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pc_x, pc_y, pc_z, skp_10, skp_11, skp_20, \
                         sls0_6, sls1_6, slp_19, slp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_9 * skp_20[k]
                  + f_3 * pc_x[k] * slp_20[k];

        t_39[k] = f_10 * skp_10[k]
                  + f_1 * sls0_6[k]
                  - f_2 * sls1_6[k]
                  + f_3 * pc_y[k] * slp_19[k];

        t_40[k] = f_10 * skp_11[k]
                  + f_3 * pc_y[k] * slp_20[k];

        t_41[k] = f_1 * sls0_6[k]
                  - f_2 * sls1_6[k]
                  + f_3 * pc_z[k] * slp_20[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pb_z, pc_x, pc_z, skd0_18, skd0_21, skp_22, \
                         skp_23, skd1_18, skd1_21, slp_22, slp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_z[k] * skd0_18[k]
                  - f_4 * pc_z[k] * skd1_18[k];

        t_43[k] = f_9 * skp_22[k]
                  + f_3 * pc_x[k] * slp_22[k];

        t_44[k] = f_9 * skp_23[k]
                  + f_3 * pc_x[k] * slp_23[k];

        t_45[k] = pb_z[k] * skd0_21[k]
                  - f_4 * pc_z[k] * skd1_21[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, pb_y, pc_y, pc_z, skd0_30, skp_11, skp_14, skd1_30, \
                         sls0_7, sls1_7, slp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_8 * skp_14[k]
                  + f_3 * pc_y[k] * slp_23[k];

        t_47[k] = f_6 * skp_11[k]
                  + f_1 * sls0_7[k]
                  - f_2 * sls1_7[k]
                  + f_3 * pc_z[k] * slp_23[k];

        t_48[k] = pb_y[k] * skd0_30[k]
                  - f_4 * pc_y[k] * skd1_30[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pc_x, pc_y, skp_16, skp_17, skp_25, skp_26, \
                         sls0_8, sls1_8, slp_25, slp_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_9 * skp_25[k]
                  + f_3 * pc_x[k] * slp_25[k];

        t_50[k] = f_9 * skp_26[k]
                  + f_3 * pc_x[k] * slp_26[k];

        t_51[k] = f_6 * skp_16[k]
                  + f_1 * sls0_8[k]
                  - f_2 * sls1_8[k]
                  + f_3 * pc_y[k] * slp_25[k];

        t_52[k] = f_6 * skp_17[k]
                  + f_3 * pc_y[k] * slp_26[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pb_y, pc_x, pc_y, skd0_35, skp_27, skp_28, skd1_35, \
                         sls0_9, sls1_9, slp_27, slp_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = pb_y[k] * skd0_35[k]
                  - f_4 * pc_y[k] * skd1_35[k];

        t_54[k] = f_9 * skp_27[k]
                  + f_1 * sls0_9[k]
                  - f_2 * sls1_9[k]
                  + f_3 * pc_x[k] * slp_27[k];

        t_55[k] = f_9 * skp_28[k]
                  + f_3 * pc_x[k] * slp_28[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pc_x, pc_y, pc_z, skp_17, skp_29, sls0_9, \
                         sls1_9, slp_28, slp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_9 * skp_29[k]
                  + f_3 * pc_x[k] * slp_29[k];

        t_57[k] = f_1 * sls0_9[k]
                  - f_2 * sls1_9[k]
                  + f_3 * pc_y[k] * slp_28[k];

        t_58[k] = f_3 * pc_y[k] * slp_29[k];

        t_59[k] = f_10 * skp_17[k]
                  + f_1 * sls0_9[k]
                  - f_2 * sls1_9[k]
                  + f_3 * pc_z[k] * slp_29[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pc_x, pc_y, skp_19, skp_30, skp_31, skp_32, \
                         sls0_10, sls1_10, slp_30, slp_31, slp_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_11 * skp_30[k]
                  + f_1 * sls0_10[k]
                  - f_2 * sls1_10[k]
                  + f_3 * pc_x[k] * slp_30[k];

        t_61[k] = f_11 * skp_31[k]
                  + f_3 * pc_x[k] * slp_31[k];

        t_62[k] = f_11 * skp_32[k]
                  + f_3 * pc_x[k] * slp_32[k];

        t_63[k] = f_11 * skp_19[k]
                  + f_1 * sls0_10[k]
                  - f_2 * sls1_10[k]
                  + f_3 * pc_y[k] * slp_31[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pb_z, pc_x, pc_y, pc_z, skd0_36, skp_20, \
                         skp_34, skd1_36, sls0_10, sls1_10, slp_32, \
                         slp_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_11 * skp_20[k]
                  + f_3 * pc_y[k] * slp_32[k];

        t_65[k] = f_1 * sls0_10[k]
                  - f_2 * sls1_10[k]
                  + f_3 * pc_z[k] * slp_32[k];

        t_66[k] = pb_z[k] * skd0_36[k]
                  - f_4 * pc_z[k] * skd1_36[k];

        t_67[k] = f_11 * skp_34[k]
                  + f_3 * pc_x[k] * slp_34[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pb_z, pc_x, pc_y, pc_z, skd0_39, skp_20, \
                         skp_23, skp_35, skd1_39, sls0_11, sls1_11, \
                         slp_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_11 * skp_35[k]
                  + f_3 * pc_x[k] * slp_35[k];

        t_69[k] = pb_z[k] * skd0_39[k]
                  - f_4 * pc_z[k] * skd1_39[k];

        t_70[k] = f_10 * skp_23[k]
                  + f_3 * pc_y[k] * slp_35[k];

        t_71[k] = f_6 * skp_20[k]
                  + f_1 * sls0_11[k]
                  - f_2 * sls1_11[k]
                  + f_3 * pc_z[k] * slp_35[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pc_x, pc_y, skp_25, skp_36, skp_37, skp_38, \
                         sls0_12, sls1_12, slp_36, slp_37, slp_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_11 * skp_36[k]
                  + f_1 * sls0_12[k]
                  - f_2 * sls1_12[k]
                  + f_3 * pc_x[k] * slp_36[k];

        t_73[k] = f_11 * skp_37[k]
                  + f_3 * pc_x[k] * slp_37[k];

        t_74[k] = f_11 * skp_38[k]
                  + f_3 * pc_x[k] * slp_38[k];

        t_75[k] = f_8 * skp_25[k]
                  + f_1 * sls0_12[k]
                  - f_2 * sls1_12[k]
                  + f_3 * pc_y[k] * slp_37[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, pb_y, pc_y, pc_z, skd0_54, skp_23, skp_26, skd1_54, \
                         sls0_12, sls1_12, slp_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_8 * skp_26[k]
                  + f_3 * pc_y[k] * slp_38[k];

        t_77[k] = f_8 * skp_23[k]
                  + f_1 * sls0_12[k]
                  - f_2 * sls1_12[k]
                  + f_3 * pc_z[k] * slp_38[k];

        t_78[k] = pb_y[k] * skd0_54[k]
                  - f_4 * pc_y[k] * skd1_54[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pc_x, pc_y, skp_28, skp_29, skp_40, skp_41, \
                         sls0_13, sls1_13, slp_40, slp_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_11 * skp_40[k]
                  + f_3 * pc_x[k] * slp_40[k];

        t_80[k] = f_11 * skp_41[k]
                  + f_3 * pc_x[k] * slp_41[k];

        t_81[k] = f_6 * skp_28[k]
                  + f_1 * sls0_13[k]
                  - f_2 * sls1_13[k]
                  + f_3 * pc_y[k] * slp_40[k];

        t_82[k] = f_6 * skp_29[k]
                  + f_3 * pc_y[k] * slp_41[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, pb_y, pc_x, pc_y, skd0_59, skp_42, skp_43, skd1_59, \
                         sls0_14, sls1_14, slp_42, slp_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = pb_y[k] * skd0_59[k]
                  - f_4 * pc_y[k] * skd1_59[k];

        t_84[k] = f_11 * skp_42[k]
                  + f_1 * sls0_14[k]
                  - f_2 * sls1_14[k]
                  + f_3 * pc_x[k] * slp_42[k];

        t_85[k] = f_11 * skp_43[k]
                  + f_3 * pc_x[k] * slp_43[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pc_x, pc_y, pc_z, skp_29, skp_44, sls0_14, \
                         sls1_14, slp_43, slp_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_11 * skp_44[k]
                  + f_3 * pc_x[k] * slp_44[k];

        t_87[k] = f_1 * sls0_14[k]
                  - f_2 * sls1_14[k]
                  + f_3 * pc_y[k] * slp_43[k];

        t_88[k] = f_3 * pc_y[k] * slp_44[k];

        t_89[k] = f_11 * skp_29[k]
                  + f_1 * sls0_14[k]
                  - f_2 * sls1_14[k]
                  + f_3 * pc_z[k] * slp_44[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, pc_x, pc_y, skp_31, skp_45, skp_46, skp_47, \
                         sls0_15, sls1_15, slp_45, slp_46, slp_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_10 * skp_45[k]
                  + f_1 * sls0_15[k]
                  - f_2 * sls1_15[k]
                  + f_3 * pc_x[k] * slp_45[k];

        t_91[k] = f_10 * skp_46[k]
                  + f_3 * pc_x[k] * slp_46[k];

        t_92[k] = f_10 * skp_47[k]
                  + f_3 * pc_x[k] * slp_47[k];

        t_93[k] = f_9 * skp_31[k]
                  + f_1 * sls0_15[k]
                  - f_2 * sls1_15[k]
                  + f_3 * pc_y[k] * slp_46[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, pb_z, pc_x, pc_y, pc_z, skd0_60, skp_32, \
                         skp_49, skd1_60, sls0_15, sls1_15, slp_47, \
                         slp_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_9 * skp_32[k]
                  + f_3 * pc_y[k] * slp_47[k];

        t_95[k] = f_1 * sls0_15[k]
                  - f_2 * sls1_15[k]
                  + f_3 * pc_z[k] * slp_47[k];

        t_96[k] = pb_z[k] * skd0_60[k]
                  - f_4 * pc_z[k] * skd1_60[k];

        t_97[k] = f_10 * skp_49[k]
                  + f_3 * pc_x[k] * slp_49[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, pb_z, pc_x, pc_y, pc_z, skd0_63, skp_32, \
                         skp_35, skp_50, skd1_63, sls0_16, sls1_16, \
                         slp_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_10 * skp_50[k]
                  + f_3 * pc_x[k] * slp_50[k];

        t_99[k] = pb_z[k] * skd0_63[k]
                  - f_4 * pc_z[k] * skd1_63[k];

        t_100[k] = f_11 * skp_35[k]
                   + f_3 * pc_y[k] * slp_50[k];

        t_101[k] = f_6 * skp_32[k]
                   + f_1 * sls0_16[k]
                   - f_2 * sls1_16[k]
                   + f_3 * pc_z[k] * slp_50[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, pc_x, pc_y, skp_37, skp_51, skp_52, \
                         skp_53, sls0_17, sls1_17, slp_51, slp_52, \
                         slp_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_10 * skp_51[k]
                   + f_1 * sls0_17[k]
                   - f_2 * sls1_17[k]
                   + f_3 * pc_x[k] * slp_51[k];

        t_103[k] = f_10 * skp_52[k]
                   + f_3 * pc_x[k] * slp_52[k];

        t_104[k] = f_10 * skp_53[k]
                   + f_3 * pc_x[k] * slp_53[k];

        t_105[k] = f_10 * skp_37[k]
                   + f_1 * sls0_17[k]
                   - f_2 * sls1_17[k]
                   + f_3 * pc_y[k] * slp_52[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, pc_x, pc_y, pc_z, skp_35, skp_38, skp_54, \
                         sls0_17, sls0_18, sls1_17, sls1_18, slp_53, \
                         slp_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_10 * skp_38[k]
                   + f_3 * pc_y[k] * slp_53[k];

        t_107[k] = f_8 * skp_35[k]
                   + f_1 * sls0_17[k]
                   - f_2 * sls1_17[k]
                   + f_3 * pc_z[k] * slp_53[k];

        t_108[k] = f_10 * skp_54[k]
                   + f_1 * sls0_18[k]
                   - f_2 * sls1_18[k]
                   + f_3 * pc_x[k] * slp_54[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pc_x, pc_y, skp_40, skp_41, skp_55, \
                         skp_56, sls0_18, sls1_18, slp_55, slp_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_10 * skp_55[k]
                   + f_3 * pc_x[k] * slp_55[k];

        t_110[k] = f_10 * skp_56[k]
                   + f_3 * pc_x[k] * slp_56[k];

        t_111[k] = f_8 * skp_40[k]
                   + f_1 * sls0_18[k]
                   - f_2 * sls1_18[k]
                   + f_3 * pc_y[k] * slp_55[k];

        t_112[k] = f_8 * skp_41[k]
                   + f_3 * pc_y[k] * slp_56[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, pb_y, pc_x, pc_y, pc_z, skd0_84, skp_38, skp_58, \
                         skd1_84, sls0_18, sls1_18, slp_56, slp_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_10 * skp_38[k]
                   + f_1 * sls0_18[k]
                   - f_2 * sls1_18[k]
                   + f_3 * pc_z[k] * slp_56[k];

        t_114[k] = pb_y[k] * skd0_84[k]
                   - f_4 * pc_y[k] * skd1_84[k];

        t_115[k] = f_10 * skp_58[k]
                   + f_3 * pc_x[k] * slp_58[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, pb_y, pc_x, pc_y, skd0_89, skp_43, \
                         skp_44, skp_59, skd1_89, sls0_19, sls1_19, slp_58, \
                         slp_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_10 * skp_59[k]
                   + f_3 * pc_x[k] * slp_59[k];

        t_117[k] = f_6 * skp_43[k]
                   + f_1 * sls0_19[k]
                   - f_2 * sls1_19[k]
                   + f_3 * pc_y[k] * slp_58[k];

        t_118[k] = f_6 * skp_44[k]
                   + f_3 * pc_y[k] * slp_59[k];

        t_119[k] = pb_y[k] * skd0_89[k]
                   - f_4 * pc_y[k] * skd1_89[k];
    }
}

static auto
compute_prim_sld_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t skd0,
                                                          const size_t skp, const size_t skd1,
                                                          const size_t sls0, const size_t sls1,
                                                          const size_t slp, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
    const auto f_1 = 0.5 / gamma;
    const auto f_2 = 0.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 3.5 / q;
    const auto f_6 = 0.5 / q;
    const auto f_7 = 3.0 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 2.5 / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 2.0 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *skd0_90 = buffer.data(skd0 + 90);
    const auto *skd0_93 = buffer.data(skd0 + 93);
    const auto *skd0_120 = buffer.data(skd0 + 120);
    const auto *skd0_125 = buffer.data(skd0 + 125);
    const auto *skd0_126 = buffer.data(skd0 + 126);
    const auto *skd0_162 = buffer.data(skd0 + 162);
    const auto *skd0_168 = buffer.data(skd0 + 168);
    const auto *skd0_171 = buffer.data(skd0 + 171);
    const auto *skd0_173 = buffer.data(skd0 + 173);
    const auto *skd0_177 = buffer.data(skd0 + 177);
    const auto *skd0_179 = buffer.data(skd0 + 179);
    const auto *skd0_180 = buffer.data(skd0 + 180);
    const auto *skd0_183 = buffer.data(skd0 + 183);
    const auto *skd0_185 = buffer.data(skd0 + 185);
    const auto *skd0_186 = buffer.data(skd0 + 186);
    const auto *skd0_189 = buffer.data(skd0 + 189);
    const auto *skd0_191 = buffer.data(skd0 + 191);
    const auto *skd0_192 = buffer.data(skd0 + 192);
    const auto *skd0_195 = buffer.data(skd0 + 195);
    const auto *skd0_197 = buffer.data(skd0 + 197);
    const auto *skd0_198 = buffer.data(skd0 + 198);
    const auto *skd0_201 = buffer.data(skd0 + 201);
    const auto *skd0_203 = buffer.data(skd0 + 203);
    const auto *skd0_207 = buffer.data(skd0 + 207);
    const auto *skd0_209 = buffer.data(skd0 + 209);
    const auto *skd0_210 = buffer.data(skd0 + 210);
    const auto *skd0_213 = buffer.data(skd0 + 213);
    const auto *skd0_215 = buffer.data(skd0 + 215);

    const auto *skp_44 = buffer.data(skp + 44);
    const auto *skp_46 = buffer.data(skp + 46);
    const auto *skp_47 = buffer.data(skp + 47);
    const auto *skp_50 = buffer.data(skp + 50);
    const auto *skp_52 = buffer.data(skp + 52);
    const auto *skp_53 = buffer.data(skp + 53);
    const auto *skp_55 = buffer.data(skp + 55);
    const auto *skp_56 = buffer.data(skp + 56);
    const auto *skp_58 = buffer.data(skp + 58);
    const auto *skp_59 = buffer.data(skp + 59);
    const auto *skp_60 = buffer.data(skp + 60);
    const auto *skp_61 = buffer.data(skp + 61);
    const auto *skp_62 = buffer.data(skp + 62);
    const auto *skp_63 = buffer.data(skp + 63);
    const auto *skp_64 = buffer.data(skp + 64);
    const auto *skp_65 = buffer.data(skp + 65);
    const auto *skp_67 = buffer.data(skp + 67);
    const auto *skp_68 = buffer.data(skp + 68);
    const auto *skp_69 = buffer.data(skp + 69);
    const auto *skp_70 = buffer.data(skp + 70);
    const auto *skp_71 = buffer.data(skp + 71);
    const auto *skp_72 = buffer.data(skp + 72);
    const auto *skp_73 = buffer.data(skp + 73);
    const auto *skp_74 = buffer.data(skp + 74);
    const auto *skp_75 = buffer.data(skp + 75);
    const auto *skp_76 = buffer.data(skp + 76);
    const auto *skp_77 = buffer.data(skp + 77);
    const auto *skp_79 = buffer.data(skp + 79);
    const auto *skp_80 = buffer.data(skp + 80);
    const auto *skp_81 = buffer.data(skp + 81);
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

    const auto *skd1_90 = buffer.data(skd1 + 90);
    const auto *skd1_93 = buffer.data(skd1 + 93);
    const auto *skd1_120 = buffer.data(skd1 + 120);
    const auto *skd1_125 = buffer.data(skd1 + 125);
    const auto *skd1_126 = buffer.data(skd1 + 126);
    const auto *skd1_162 = buffer.data(skd1 + 162);
    const auto *skd1_168 = buffer.data(skd1 + 168);
    const auto *skd1_171 = buffer.data(skd1 + 171);
    const auto *skd1_173 = buffer.data(skd1 + 173);
    const auto *skd1_177 = buffer.data(skd1 + 177);
    const auto *skd1_179 = buffer.data(skd1 + 179);
    const auto *skd1_180 = buffer.data(skd1 + 180);
    const auto *skd1_183 = buffer.data(skd1 + 183);
    const auto *skd1_185 = buffer.data(skd1 + 185);
    const auto *skd1_186 = buffer.data(skd1 + 186);
    const auto *skd1_189 = buffer.data(skd1 + 189);
    const auto *skd1_191 = buffer.data(skd1 + 191);
    const auto *skd1_192 = buffer.data(skd1 + 192);
    const auto *skd1_195 = buffer.data(skd1 + 195);
    const auto *skd1_197 = buffer.data(skd1 + 197);
    const auto *skd1_198 = buffer.data(skd1 + 198);
    const auto *skd1_201 = buffer.data(skd1 + 201);
    const auto *skd1_203 = buffer.data(skd1 + 203);
    const auto *skd1_207 = buffer.data(skd1 + 207);
    const auto *skd1_209 = buffer.data(skd1 + 209);
    const auto *skd1_210 = buffer.data(skd1 + 210);
    const auto *skd1_213 = buffer.data(skd1 + 213);
    const auto *skd1_215 = buffer.data(skd1 + 215);

    const auto *sls0_20 = buffer.data(sls0 + 20);
    const auto *sls0_21 = buffer.data(sls0 + 21);
    const auto *sls0_22 = buffer.data(sls0 + 22);
    const auto *sls0_23 = buffer.data(sls0 + 23);
    const auto *sls0_24 = buffer.data(sls0 + 24);
    const auto *sls0_25 = buffer.data(sls0 + 25);
    const auto *sls0_26 = buffer.data(sls0 + 26);
    const auto *sls0_27 = buffer.data(sls0 + 27);
    const auto *sls0_36 = buffer.data(sls0 + 36);
    const auto *sls0_37 = buffer.data(sls0 + 37);
    const auto *sls0_38 = buffer.data(sls0 + 38);
    const auto *sls0_39 = buffer.data(sls0 + 39);
    const auto *sls0_40 = buffer.data(sls0 + 40);

    const auto *sls1_20 = buffer.data(sls1 + 20);
    const auto *sls1_21 = buffer.data(sls1 + 21);
    const auto *sls1_22 = buffer.data(sls1 + 22);
    const auto *sls1_23 = buffer.data(sls1 + 23);
    const auto *sls1_24 = buffer.data(sls1 + 24);
    const auto *sls1_25 = buffer.data(sls1 + 25);
    const auto *sls1_26 = buffer.data(sls1 + 26);
    const auto *sls1_27 = buffer.data(sls1 + 27);
    const auto *sls1_36 = buffer.data(sls1 + 36);
    const auto *sls1_37 = buffer.data(sls1 + 37);
    const auto *sls1_38 = buffer.data(sls1 + 38);
    const auto *sls1_39 = buffer.data(sls1 + 39);
    const auto *sls1_40 = buffer.data(sls1 + 40);

    const auto *slp_60 = buffer.data(slp + 60);
    const auto *slp_61 = buffer.data(slp + 61);
    const auto *slp_62 = buffer.data(slp + 62);
    const auto *slp_63 = buffer.data(slp + 63);
    const auto *slp_64 = buffer.data(slp + 64);
    const auto *slp_65 = buffer.data(slp + 65);
    const auto *slp_67 = buffer.data(slp + 67);
    const auto *slp_68 = buffer.data(slp + 68);
    const auto *slp_69 = buffer.data(slp + 69);
    const auto *slp_70 = buffer.data(slp + 70);
    const auto *slp_71 = buffer.data(slp + 71);
    const auto *slp_72 = buffer.data(slp + 72);
    const auto *slp_73 = buffer.data(slp + 73);
    const auto *slp_74 = buffer.data(slp + 74);
    const auto *slp_75 = buffer.data(slp + 75);
    const auto *slp_76 = buffer.data(slp + 76);
    const auto *slp_77 = buffer.data(slp + 77);
    const auto *slp_79 = buffer.data(slp + 79);
    const auto *slp_80 = buffer.data(slp + 80);
    const auto *slp_81 = buffer.data(slp + 81);
    const auto *slp_82 = buffer.data(slp + 82);
    const auto *slp_83 = buffer.data(slp + 83);
    const auto *slp_85 = buffer.data(slp + 85);
    const auto *slp_86 = buffer.data(slp + 86);
    const auto *slp_88 = buffer.data(slp + 88);
    const auto *slp_89 = buffer.data(slp + 89);
    const auto *slp_91 = buffer.data(slp + 91);
    const auto *slp_92 = buffer.data(slp + 92);
    const auto *slp_94 = buffer.data(slp + 94);
    const auto *slp_95 = buffer.data(slp + 95);
    const auto *slp_97 = buffer.data(slp + 97);
    const auto *slp_98 = buffer.data(slp + 98);
    const auto *slp_100 = buffer.data(slp + 100);
    const auto *slp_101 = buffer.data(slp + 101);
    const auto *slp_103 = buffer.data(slp + 103);
    const auto *slp_104 = buffer.data(slp + 104);
    const auto *slp_106 = buffer.data(slp + 106);
    const auto *slp_107 = buffer.data(slp + 107);
    const auto *slp_108 = buffer.data(slp + 108);
    const auto *slp_109 = buffer.data(slp + 109);
    const auto *slp_110 = buffer.data(slp + 110);
    const auto *slp_112 = buffer.data(slp + 112);
    const auto *slp_113 = buffer.data(slp + 113);
    const auto *slp_114 = buffer.data(slp + 114);
    const auto *slp_115 = buffer.data(slp + 115);
    const auto *slp_116 = buffer.data(slp + 116);
    const auto *slp_117 = buffer.data(slp + 117);
    const auto *slp_118 = buffer.data(slp + 118);
    const auto *slp_119 = buffer.data(slp + 119);
    const auto *slp_120 = buffer.data(slp + 120);
    const auto *slp_121 = buffer.data(slp + 121);
    const auto *slp_122 = buffer.data(slp + 122);

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, pc_x, pc_y, skp_60, skp_61, \
                         skp_62, sls0_20, sls1_20, slp_60, slp_61, \
                         slp_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_10 * skp_60[k]
                   + f_1 * sls0_20[k]
                   - f_2 * sls1_20[k]
                   + f_3 * pc_x[k] * slp_60[k];

        t_121[k] = f_10 * skp_61[k]
                   + f_3 * pc_x[k] * slp_61[k];

        t_122[k] = f_10 * skp_62[k]
                   + f_3 * pc_x[k] * slp_62[k];

        t_123[k] = f_1 * sls0_20[k]
                   - f_2 * sls1_20[k]
                   + f_3 * pc_y[k] * slp_61[k];

        t_124[k] = f_3 * pc_y[k] * slp_62[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, pc_x, pc_z, skp_44, skp_63, skp_64, sls0_20, \
                         sls0_21, sls1_20, sls1_21, slp_62, slp_63, \
                         slp_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_9 * skp_44[k]
                   + f_1 * sls0_20[k]
                   - f_2 * sls1_20[k]
                   + f_3 * pc_z[k] * slp_62[k];

        t_126[k] = f_8 * skp_63[k]
                   + f_1 * sls0_21[k]
                   - f_2 * sls1_21[k]
                   + f_3 * pc_x[k] * slp_63[k];

        t_127[k] = f_8 * skp_64[k]
                   + f_3 * pc_x[k] * slp_64[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, pc_x, pc_y, pc_z, skp_46, skp_47, skp_65, \
                         sls0_21, sls1_21, slp_64, slp_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_8 * skp_65[k]
                   + f_3 * pc_x[k] * slp_65[k];

        t_129[k] = f_7 * skp_46[k]
                   + f_1 * sls0_21[k]
                   - f_2 * sls1_21[k]
                   + f_3 * pc_y[k] * slp_64[k];

        t_130[k] = f_7 * skp_47[k]
                   + f_3 * pc_y[k] * slp_65[k];

        t_131[k] = f_1 * sls0_21[k]
                   - f_2 * sls1_21[k]
                   + f_3 * pc_z[k] * slp_65[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pb_z, pc_x, pc_z, skd0_90, skd0_93, \
                         skp_67, skp_68, skd1_90, skd1_93, slp_67, \
                         slp_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = pb_z[k] * skd0_90[k]
                   - f_4 * pc_z[k] * skd1_90[k];

        t_133[k] = f_8 * skp_67[k]
                   + f_3 * pc_x[k] * slp_67[k];

        t_134[k] = f_8 * skp_68[k]
                   + f_3 * pc_x[k] * slp_68[k];

        t_135[k] = pb_z[k] * skd0_93[k]
                   - f_4 * pc_z[k] * skd1_93[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, pc_x, pc_y, pc_z, skp_47, skp_50, skp_69, \
                         sls0_22, sls0_23, sls1_22, sls1_23, slp_68, \
                         slp_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_9 * skp_50[k]
                   + f_3 * pc_y[k] * slp_68[k];

        t_137[k] = f_6 * skp_47[k]
                   + f_1 * sls0_22[k]
                   - f_2 * sls1_22[k]
                   + f_3 * pc_z[k] * slp_68[k];

        t_138[k] = f_8 * skp_69[k]
                   + f_1 * sls0_23[k]
                   - f_2 * sls1_23[k]
                   + f_3 * pc_x[k] * slp_69[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, pc_x, pc_y, skp_52, skp_53, skp_70, \
                         skp_71, sls0_23, sls1_23, slp_70, slp_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_8 * skp_70[k]
                   + f_3 * pc_x[k] * slp_70[k];

        t_140[k] = f_8 * skp_71[k]
                   + f_3 * pc_x[k] * slp_71[k];

        t_141[k] = f_11 * skp_52[k]
                   + f_1 * sls0_23[k]
                   - f_2 * sls1_23[k]
                   + f_3 * pc_y[k] * slp_70[k];

        t_142[k] = f_11 * skp_53[k]
                   + f_3 * pc_y[k] * slp_71[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, pc_x, pc_z, skp_50, skp_72, skp_73, sls0_23, \
                         sls0_24, sls1_23, sls1_24, slp_71, slp_72, \
                         slp_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_8 * skp_50[k]
                   + f_1 * sls0_23[k]
                   - f_2 * sls1_23[k]
                   + f_3 * pc_z[k] * slp_71[k];

        t_144[k] = f_8 * skp_72[k]
                   + f_1 * sls0_24[k]
                   - f_2 * sls1_24[k]
                   + f_3 * pc_x[k] * slp_72[k];

        t_145[k] = f_8 * skp_73[k]
                   + f_3 * pc_x[k] * slp_73[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pc_x, pc_y, pc_z, skp_53, skp_55, skp_56, \
                         skp_74, sls0_24, sls1_24, slp_73, slp_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_8 * skp_74[k]
                   + f_3 * pc_x[k] * slp_74[k];

        t_147[k] = f_10 * skp_55[k]
                   + f_1 * sls0_24[k]
                   - f_2 * sls1_24[k]
                   + f_3 * pc_y[k] * slp_73[k];

        t_148[k] = f_10 * skp_56[k]
                   + f_3 * pc_y[k] * slp_74[k];

        t_149[k] = f_10 * skp_53[k]
                   + f_1 * sls0_24[k]
                   - f_2 * sls1_24[k]
                   + f_3 * pc_z[k] * slp_74[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, pc_x, pc_y, skp_58, skp_75, skp_76, \
                         skp_77, sls0_25, sls1_25, slp_75, slp_76, \
                         slp_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_8 * skp_75[k]
                   + f_1 * sls0_25[k]
                   - f_2 * sls1_25[k]
                   + f_3 * pc_x[k] * slp_75[k];

        t_151[k] = f_8 * skp_76[k]
                   + f_3 * pc_x[k] * slp_76[k];

        t_152[k] = f_8 * skp_77[k]
                   + f_3 * pc_x[k] * slp_77[k];

        t_153[k] = f_8 * skp_58[k]
                   + f_1 * sls0_25[k]
                   - f_2 * sls1_25[k]
                   + f_3 * pc_y[k] * slp_76[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, pb_y, pc_y, pc_z, skd0_120, skp_56, skp_59, \
                         skd1_120, sls0_25, sls1_25, slp_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = f_8 * skp_59[k]
                   + f_3 * pc_y[k] * slp_77[k];

        t_155[k] = f_11 * skp_56[k]
                   + f_1 * sls0_25[k]
                   - f_2 * sls1_25[k]
                   + f_3 * pc_z[k] * slp_77[k];

        t_156[k] = pb_y[k] * skd0_120[k]
                   - f_4 * pc_y[k] * skd1_120[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, pc_x, pc_y, skp_61, skp_62, skp_79, \
                         skp_80, sls0_26, sls1_26, slp_79, slp_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = f_8 * skp_79[k]
                   + f_3 * pc_x[k] * slp_79[k];

        t_158[k] = f_8 * skp_80[k]
                   + f_3 * pc_x[k] * slp_80[k];

        t_159[k] = f_6 * skp_61[k]
                   + f_1 * sls0_26[k]
                   - f_2 * sls1_26[k]
                   + f_3 * pc_y[k] * slp_79[k];

        t_160[k] = f_6 * skp_62[k]
                   + f_3 * pc_y[k] * slp_80[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, pb_y, pc_x, pc_y, skd0_125, skp_81, skp_82, \
                         skd1_125, sls0_27, sls1_27, slp_81, slp_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = pb_y[k] * skd0_125[k]
                   - f_4 * pc_y[k] * skd1_125[k];

        t_162[k] = f_8 * skp_81[k]
                   + f_1 * sls0_27[k]
                   - f_2 * sls1_27[k]
                   + f_3 * pc_x[k] * slp_81[k];

        t_163[k] = f_8 * skp_82[k]
                   + f_3 * pc_x[k] * slp_82[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, pc_x, pc_y, pc_z, skp_62, skp_83, \
                         sls0_27, sls1_27, slp_82, slp_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_8 * skp_83[k]
                   + f_3 * pc_x[k] * slp_83[k];

        t_165[k] = f_1 * sls0_27[k]
                   - f_2 * sls1_27[k]
                   + f_3 * pc_y[k] * slp_82[k];

        t_166[k] = f_3 * pc_y[k] * slp_83[k];

        t_167[k] = f_7 * skp_62[k]
                   + f_1 * sls0_27[k]
                   - f_2 * sls1_27[k]
                   + f_3 * pc_z[k] * slp_83[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, pb_x, pc_x, skd0_168, skd0_171, skp_84, \
                         skp_85, skp_86, skd1_168, skd1_171, slp_85, \
                         slp_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = pb_x[k] * skd0_168[k]
                   + f_8 * skp_84[k]
                   - f_4 * pc_x[k] * skd1_168[k];

        t_169[k] = f_6 * skp_85[k]
                   + f_3 * pc_x[k] * slp_85[k];

        t_170[k] = f_6 * skp_86[k]
                   + f_3 * pc_x[k] * slp_86[k];

        t_171[k] = pb_x[k] * skd0_171[k]
                   - f_4 * pc_x[k] * skd1_171[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, pb_x, pb_z, pc_x, pc_y, pc_z, skd0_126, \
                         skd0_173, skp_65, skd1_126, skd1_173, slp_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_5 * skp_65[k]
                   + f_3 * pc_y[k] * slp_86[k];

        t_173[k] = pb_x[k] * skd0_173[k]
                   - f_4 * pc_x[k] * skd1_173[k];

        t_174[k] = pb_z[k] * skd0_126[k]
                   - f_4 * pc_z[k] * skd1_126[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, pb_x, pc_x, pc_y, skd0_177, skp_68, \
                         skp_88, skp_89, skd1_177, slp_88, slp_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = f_6 * skp_88[k]
                   + f_3 * pc_x[k] * slp_88[k];

        t_176[k] = f_6 * skp_89[k]
                   + f_3 * pc_x[k] * slp_89[k];

        t_177[k] = pb_x[k] * skd0_177[k]
                   - f_4 * pc_x[k] * skd1_177[k];

        t_178[k] = f_7 * skp_68[k]
                   + f_3 * pc_y[k] * slp_89[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, t_182, pb_x, pc_x, skd0_179, skd0_180, skp_90, \
                         skp_91, skp_92, skd1_179, skd1_180, slp_91, \
                         slp_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = pb_x[k] * skd0_179[k]
                   - f_4 * pc_x[k] * skd1_179[k];

        t_180[k] = pb_x[k] * skd0_180[k]
                   + f_8 * skp_90[k]
                   - f_4 * pc_x[k] * skd1_180[k];

        t_181[k] = f_6 * skp_91[k]
                   + f_3 * pc_x[k] * slp_91[k];

        t_182[k] = f_6 * skp_92[k]
                   + f_3 * pc_x[k] * slp_92[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, t_186, pb_x, pc_x, pc_y, skd0_183, skd0_185, \
                         skd0_186, skp_71, skp_93, skd1_183, skd1_185, skd1_186, \
                         slp_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = pb_x[k] * skd0_183[k]
                   - f_4 * pc_x[k] * skd1_183[k];

        t_184[k] = f_9 * skp_71[k]
                   + f_3 * pc_y[k] * slp_92[k];

        t_185[k] = pb_x[k] * skd0_185[k]
                   - f_4 * pc_x[k] * skd1_185[k];

        t_186[k] = pb_x[k] * skd0_186[k]
                   + f_8 * skp_93[k]
                   - f_4 * pc_x[k] * skd1_186[k];
    }

#pragma omp simd aligned(t_187, t_188, t_189, t_190, pb_x, pc_x, pc_y, skd0_189, skp_74, \
                         skp_94, skp_95, skd1_189, slp_94, slp_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_187[k] = f_6 * skp_94[k]
                   + f_3 * pc_x[k] * slp_94[k];

        t_188[k] = f_6 * skp_95[k]
                   + f_3 * pc_x[k] * slp_95[k];

        t_189[k] = pb_x[k] * skd0_189[k]
                   - f_4 * pc_x[k] * skd1_189[k];

        t_190[k] = f_11 * skp_74[k]
                   + f_3 * pc_y[k] * slp_95[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, t_194, pb_x, pc_x, skd0_191, skd0_192, skp_96, \
                         skp_97, skp_98, skd1_191, skd1_192, slp_97, \
                         slp_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = pb_x[k] * skd0_191[k]
                   - f_4 * pc_x[k] * skd1_191[k];

        t_192[k] = pb_x[k] * skd0_192[k]
                   + f_8 * skp_96[k]
                   - f_4 * pc_x[k] * skd1_192[k];

        t_193[k] = f_6 * skp_97[k]
                   + f_3 * pc_x[k] * slp_97[k];

        t_194[k] = f_6 * skp_98[k]
                   + f_3 * pc_x[k] * slp_98[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, pb_x, pc_x, pc_y, skd0_195, skd0_197, \
                         skd0_198, skp_77, skp_99, skd1_195, skd1_197, skd1_198, \
                         slp_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = pb_x[k] * skd0_195[k]
                   - f_4 * pc_x[k] * skd1_195[k];

        t_196[k] = f_10 * skp_77[k]
                   + f_3 * pc_y[k] * slp_98[k];

        t_197[k] = pb_x[k] * skd0_197[k]
                   - f_4 * pc_x[k] * skd1_197[k];

        t_198[k] = pb_x[k] * skd0_198[k]
                   + f_8 * skp_99[k]
                   - f_4 * pc_x[k] * skd1_198[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, pb_x, pc_x, pc_y, skd0_201, skp_80, \
                         skp_100, skp_101, skd1_201, slp_100, slp_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = f_6 * skp_100[k]
                   + f_3 * pc_x[k] * slp_100[k];

        t_200[k] = f_6 * skp_101[k]
                   + f_3 * pc_x[k] * slp_101[k];

        t_201[k] = pb_x[k] * skd0_201[k]
                   - f_4 * pc_x[k] * skd1_201[k];

        t_202[k] = f_8 * skp_80[k]
                   + f_3 * pc_y[k] * slp_101[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, pb_x, pb_y, pc_x, pc_y, skd0_162, \
                         skd0_203, skp_103, skp_104, skd1_162, skd1_203, slp_103, \
                         slp_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = pb_x[k] * skd0_203[k]
                   - f_4 * pc_x[k] * skd1_203[k];

        t_204[k] = pb_y[k] * skd0_162[k]
                   - f_4 * pc_y[k] * skd1_162[k];

        t_205[k] = f_6 * skp_103[k]
                   + f_3 * pc_x[k] * slp_103[k];

        t_206[k] = f_6 * skp_104[k]
                   + f_3 * pc_x[k] * slp_104[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, t_210, pb_x, pc_x, pc_y, skd0_207, skd0_209, \
                         skd0_210, skp_83, skp_105, skd1_207, skd1_209, skd1_210, \
                         slp_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = pb_x[k] * skd0_207[k]
                   - f_4 * pc_x[k] * skd1_207[k];

        t_208[k] = f_6 * skp_83[k]
                   + f_3 * pc_y[k] * slp_104[k];

        t_209[k] = pb_x[k] * skd0_209[k]
                   - f_4 * pc_x[k] * skd1_209[k];

        t_210[k] = pb_x[k] * skd0_210[k]
                   + f_8 * skp_105[k]
                   - f_4 * pc_x[k] * skd1_210[k];
    }

#pragma omp simd aligned(t_211, t_212, t_213, t_214, t_215, pb_x, pc_x, pc_y, skd0_213, \
                         skd0_215, skp_106, skp_107, skd1_213, skd1_215, slp_106, \
                         slp_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_211[k] = f_6 * skp_106[k]
                   + f_3 * pc_x[k] * slp_106[k];

        t_212[k] = f_6 * skp_107[k]
                   + f_3 * pc_x[k] * slp_107[k];

        t_213[k] = pb_x[k] * skd0_213[k]
                   - f_4 * pc_x[k] * skd1_213[k];

        t_214[k] = f_3 * pc_y[k] * slp_107[k];

        t_215[k] = pb_x[k] * skd0_215[k]
                   - f_4 * pc_x[k] * skd1_215[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, t_220, t_221, pc_x, pc_y, pc_z, skp_85, \
                         skp_86, sls0_36, sls1_36, slp_108, slp_109, \
                         slp_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = f_1 * sls0_36[k]
                   - f_2 * sls1_36[k]
                   + f_3 * pc_x[k] * slp_108[k];

        t_217[k] = f_3 * pc_x[k] * slp_109[k];

        t_218[k] = f_3 * pc_x[k] * slp_110[k];

        t_219[k] = f_0 * skp_85[k]
                   + f_1 * sls0_36[k]
                   - f_2 * sls1_36[k]
                   + f_3 * pc_y[k] * slp_109[k];

        t_220[k] = f_0 * skp_86[k]
                   + f_3 * pc_y[k] * slp_110[k];

        t_221[k] = f_1 * sls0_36[k]
                   - f_2 * sls1_36[k]
                   + f_3 * pc_z[k] * slp_110[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, t_226, pb_z, pc_x, pc_y, pc_z, skd0_168, \
                         skd0_171, skp_89, skd1_168, skd1_171, slp_112, \
                         slp_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = pb_z[k] * skd0_168[k]
                   - f_4 * pc_z[k] * skd1_168[k];

        t_223[k] = f_3 * pc_x[k] * slp_112[k];

        t_224[k] = f_3 * pc_x[k] * slp_113[k];

        t_225[k] = pb_z[k] * skd0_171[k]
                   - f_4 * pc_z[k] * skd1_171[k];

        t_226[k] = f_5 * skp_89[k]
                   + f_3 * pc_y[k] * slp_113[k];
    }

#pragma omp simd aligned(t_227, t_228, t_229, t_230, pc_x, pc_z, skp_86, sls0_37, sls0_38, \
                         sls1_37, sls1_38, slp_113, slp_114, slp_115, \
                         slp_116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_227[k] = f_6 * skp_86[k]
                   + f_1 * sls0_37[k]
                   - f_2 * sls1_37[k]
                   + f_3 * pc_z[k] * slp_113[k];

        t_228[k] = f_1 * sls0_38[k]
                   - f_2 * sls1_38[k]
                   + f_3 * pc_x[k] * slp_114[k];

        t_229[k] = f_3 * pc_x[k] * slp_115[k];

        t_230[k] = f_3 * pc_x[k] * slp_116[k];
    }

#pragma omp simd aligned(t_231, t_232, t_233, pc_y, pc_z, skp_89, skp_91, skp_92, sls0_38, \
                         sls1_38, slp_115, slp_116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_231[k] = f_7 * skp_91[k]
                   + f_1 * sls0_38[k]
                   - f_2 * sls1_38[k]
                   + f_3 * pc_y[k] * slp_115[k];

        t_232[k] = f_7 * skp_92[k]
                   + f_3 * pc_y[k] * slp_116[k];

        t_233[k] = f_8 * skp_89[k]
                   + f_1 * sls0_38[k]
                   - f_2 * sls1_38[k]
                   + f_3 * pc_z[k] * slp_116[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, t_238, pc_x, pc_y, skp_94, skp_95, \
                         sls0_39, sls1_39, slp_117, slp_118, slp_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_1 * sls0_39[k]
                   - f_2 * sls1_39[k]
                   + f_3 * pc_x[k] * slp_117[k];

        t_235[k] = f_3 * pc_x[k] * slp_118[k];

        t_236[k] = f_3 * pc_x[k] * slp_119[k];

        t_237[k] = f_9 * skp_94[k]
                   + f_1 * sls0_39[k]
                   - f_2 * sls1_39[k]
                   + f_3 * pc_y[k] * slp_118[k];

        t_238[k] = f_9 * skp_95[k]
                   + f_3 * pc_y[k] * slp_119[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, pc_x, pc_z, skp_92, sls0_39, sls0_40, \
                         sls1_39, sls1_40, slp_119, slp_120, slp_121, \
                         slp_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_10 * skp_92[k]
                   + f_1 * sls0_39[k]
                   - f_2 * sls1_39[k]
                   + f_3 * pc_z[k] * slp_119[k];

        t_240[k] = f_1 * sls0_40[k]
                   - f_2 * sls1_40[k]
                   + f_3 * pc_x[k] * slp_120[k];

        t_241[k] = f_3 * pc_x[k] * slp_121[k];

        t_242[k] = f_3 * pc_x[k] * slp_122[k];
    }

#pragma omp simd aligned(t_243, t_244, t_245, pc_y, pc_z, skp_95, skp_97, skp_98, sls0_40, \
                         sls1_40, slp_121, slp_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = f_11 * skp_97[k]
                   + f_1 * sls0_40[k]
                   - f_2 * sls1_40[k]
                   + f_3 * pc_y[k] * slp_121[k];

        t_244[k] = f_11 * skp_98[k]
                   + f_3 * pc_y[k] * slp_122[k];

        t_245[k] = f_11 * skp_95[k]
                   + f_1 * sls0_40[k]
                   - f_2 * sls1_40[k]
                   + f_3 * pc_z[k] * slp_122[k];
    }
}

static auto
compute_prim_sld_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t skd0,
                                                          const size_t skp, const size_t skd1,
                                                          const size_t sls0, const size_t sls1,
                                                          const size_t slp, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
    const auto f_1 = 0.5 / gamma;
    const auto f_2 = 0.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_6 = 0.5 / q;
    const auto f_7 = 3.0 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 2.5 / q;
    const auto f_10 = 1.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *skd0_210 = buffer.data(skd0 + 210);
    const auto *skd0_213 = buffer.data(skd0 + 213);
    const auto *skd0_215 = buffer.data(skd0 + 215);

    const auto *skp_98 = buffer.data(skp + 98);
    const auto *skp_100 = buffer.data(skp + 100);
    const auto *skp_101 = buffer.data(skp + 101);
    const auto *skp_103 = buffer.data(skp + 103);
    const auto *skp_104 = buffer.data(skp + 104);
    const auto *skp_106 = buffer.data(skp + 106);
    const auto *skp_107 = buffer.data(skp + 107);

    const auto *skd1_210 = buffer.data(skd1 + 210);
    const auto *skd1_213 = buffer.data(skd1 + 213);
    const auto *skd1_215 = buffer.data(skd1 + 215);

    const auto *sls0_41 = buffer.data(sls0 + 41);
    const auto *sls0_42 = buffer.data(sls0 + 42);
    const auto *sls0_44 = buffer.data(sls0 + 44);

    const auto *sls1_41 = buffer.data(sls1 + 41);
    const auto *sls1_42 = buffer.data(sls1 + 42);
    const auto *sls1_44 = buffer.data(sls1 + 44);

    const auto *slp_123 = buffer.data(slp + 123);
    const auto *slp_124 = buffer.data(slp + 124);
    const auto *slp_125 = buffer.data(slp + 125);
    const auto *slp_126 = buffer.data(slp + 126);
    const auto *slp_127 = buffer.data(slp + 127);
    const auto *slp_128 = buffer.data(slp + 128);
    const auto *slp_130 = buffer.data(slp + 130);
    const auto *slp_131 = buffer.data(slp + 131);
    const auto *slp_132 = buffer.data(slp + 132);
    const auto *slp_133 = buffer.data(slp + 133);
    const auto *slp_134 = buffer.data(slp + 134);

#pragma omp simd aligned(t_246, t_247, t_248, t_249, t_250, pc_x, pc_y, skp_100, skp_101, \
                         sls0_41, sls1_41, slp_123, slp_124, slp_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = f_1 * sls0_41[k]
                   - f_2 * sls1_41[k]
                   + f_3 * pc_x[k] * slp_123[k];

        t_247[k] = f_3 * pc_x[k] * slp_124[k];

        t_248[k] = f_3 * pc_x[k] * slp_125[k];

        t_249[k] = f_10 * skp_100[k]
                   + f_1 * sls0_41[k]
                   - f_2 * sls1_41[k]
                   + f_3 * pc_y[k] * slp_124[k];

        t_250[k] = f_10 * skp_101[k]
                   + f_3 * pc_y[k] * slp_125[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, pc_x, pc_z, skp_98, sls0_41, sls0_42, \
                         sls1_41, sls1_42, slp_125, slp_126, slp_127, \
                         slp_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = f_9 * skp_98[k]
                   + f_1 * sls0_41[k]
                   - f_2 * sls1_41[k]
                   + f_3 * pc_z[k] * slp_125[k];

        t_252[k] = f_1 * sls0_42[k]
                   - f_2 * sls1_42[k]
                   + f_3 * pc_x[k] * slp_126[k];

        t_253[k] = f_3 * pc_x[k] * slp_127[k];

        t_254[k] = f_3 * pc_x[k] * slp_128[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, pb_y, pc_y, pc_z, skd0_210, skp_101, \
                         skp_103, skp_104, skd1_210, sls0_42, sls1_42, slp_127, \
                         slp_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_8 * skp_103[k]
                   + f_1 * sls0_42[k]
                   - f_2 * sls1_42[k]
                   + f_3 * pc_y[k] * slp_127[k];

        t_256[k] = f_8 * skp_104[k]
                   + f_3 * pc_y[k] * slp_128[k];

        t_257[k] = f_7 * skp_101[k]
                   + f_1 * sls0_42[k]
                   - f_2 * sls1_42[k]
                   + f_3 * pc_z[k] * slp_128[k];

        t_258[k] = pb_y[k] * skd0_210[k]
                   - f_4 * pc_y[k] * skd1_210[k];
    }

#pragma omp simd aligned(t_259, t_260, t_261, t_262, t_263, pb_y, pc_x, pc_y, skd0_213, \
                         skd0_215, skp_106, skp_107, skd1_213, skd1_215, slp_130, \
                         slp_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_259[k] = f_3 * pc_x[k] * slp_130[k];

        t_260[k] = f_3 * pc_x[k] * slp_131[k];

        t_261[k] = pb_y[k] * skd0_213[k]
                   + f_8 * skp_106[k]
                   - f_4 * pc_y[k] * skd1_213[k];

        t_262[k] = f_6 * skp_107[k]
                   + f_3 * pc_y[k] * slp_131[k];

        t_263[k] = pb_y[k] * skd0_215[k]
                   - f_4 * pc_y[k] * skd1_215[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, t_268, t_269, pc_x, pc_y, pc_z, skp_107, \
                         sls0_44, sls1_44, slp_132, slp_133, slp_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_1 * sls0_44[k]
                   - f_2 * sls1_44[k]
                   + f_3 * pc_x[k] * slp_132[k];

        t_265[k] = f_3 * pc_x[k] * slp_133[k];

        t_266[k] = f_3 * pc_x[k] * slp_134[k];

        t_267[k] = f_1 * sls0_44[k]
                   - f_2 * sls1_44[k]
                   + f_3 * pc_y[k] * slp_133[k];

        t_268[k] = f_3 * pc_y[k] * slp_134[k];

        t_269[k] = f_0 * skp_107[k]
                   + f_1 * sls0_44[k]
                   - f_2 * sls1_44[k]
                   + f_3 * pc_z[k] * slp_134[k];
    }
}

auto
compute_prim_sld_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t skd0, const size_t skp,
                                                   const size_t skd1, const size_t sls0,
                                                   const size_t sls1, const size_t slp,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_sld_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, skd0, skp,
                                                              skd1, sls0, sls1, slp, ncols,
                                                              gamma, p, q);

    compute_prim_sld_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, skd0, skp,
                                                              skd1, sls0, sls1, slp, ncols,
                                                              gamma, p, q);

    compute_prim_sld_three_center_electron_repulsion_0_piece2(buffer, target, pb, pc, skd0, skp,
                                                              skd1, sls0, sls1, slp, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
