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


#include "SimdThreeCenterElectronRepulsionVrrRecSKK.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_skk_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sik0,
                                                          const size_t sii, const size_t sik1,
                                                          const size_t skh0, const size_t skh1,
                                                          const size_t ski, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / q;
    const auto f_1 = 3.0 / gamma;
    const auto f_2 = 3.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.0 / gamma;
    const auto f_5 = 2.0 * p / (gamma * q);
    const auto f_6 = 1.5 / gamma;
    const auto f_7 = 1.5 * p / (gamma * q);
    const auto f_8 = 1.0 / gamma;
    const auto f_9 = p / (gamma * q);
    const auto f_10 = 0.5 / gamma;
    const auto f_11 = 0.5 * p / (gamma * q);
    const auto f_12 = gamma / q;
    const auto f_13 = 0.5 / q;
    const auto f_14 = 1.0 / q;
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;
    const auto f_18 = 3.0 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sik0_0 = buffer.data(sik0 + 0);
    const auto *sik0_3 = buffer.data(sik0 + 3);
    const auto *sik0_5 = buffer.data(sik0 + 5);
    const auto *sik0_6 = buffer.data(sik0 + 6);
    const auto *sik0_9 = buffer.data(sik0 + 9);
    const auto *sik0_10 = buffer.data(sik0 + 10);
    const auto *sik0_12 = buffer.data(sik0 + 12);
    const auto *sik0_14 = buffer.data(sik0 + 14);
    const auto *sik0_15 = buffer.data(sik0 + 15);
    const auto *sik0_17 = buffer.data(sik0 + 17);
    const auto *sik0_18 = buffer.data(sik0 + 18);
    const auto *sik0_20 = buffer.data(sik0 + 20);
    const auto *sik0_28 = buffer.data(sik0 + 28);
    const auto *sik0_35 = buffer.data(sik0 + 35);

    const auto *sii_0 = buffer.data(sii + 0);
    const auto *sii_1 = buffer.data(sii + 1);
    const auto *sii_2 = buffer.data(sii + 2);
    const auto *sii_3 = buffer.data(sii + 3);
    const auto *sii_5 = buffer.data(sii + 5);
    const auto *sii_6 = buffer.data(sii + 6);
    const auto *sii_7 = buffer.data(sii + 7);
    const auto *sii_8 = buffer.data(sii + 8);
    const auto *sii_9 = buffer.data(sii + 9);
    const auto *sii_10 = buffer.data(sii + 10);
    const auto *sii_11 = buffer.data(sii + 11);
    const auto *sii_12 = buffer.data(sii + 12);
    const auto *sii_13 = buffer.data(sii + 13);
    const auto *sii_14 = buffer.data(sii + 14);
    const auto *sii_15 = buffer.data(sii + 15);
    const auto *sii_17 = buffer.data(sii + 17);
    const auto *sii_18 = buffer.data(sii + 18);
    const auto *sii_20 = buffer.data(sii + 20);
    const auto *sii_21 = buffer.data(sii + 21);
    const auto *sii_22 = buffer.data(sii + 22);
    const auto *sii_23 = buffer.data(sii + 23);
    const auto *sii_24 = buffer.data(sii + 24);
    const auto *sii_25 = buffer.data(sii + 25);
    const auto *sii_26 = buffer.data(sii + 26);
    const auto *sii_27 = buffer.data(sii + 27);
    const auto *sii_28 = buffer.data(sii + 28);
    const auto *sii_30 = buffer.data(sii + 30);
    const auto *sii_33 = buffer.data(sii + 33);
    const auto *sii_37 = buffer.data(sii + 37);
    const auto *sii_49 = buffer.data(sii + 49);
    const auto *sii_50 = buffer.data(sii + 50);
    const auto *sii_51 = buffer.data(sii + 51);
    const auto *sii_52 = buffer.data(sii + 52);
    const auto *sii_53 = buffer.data(sii + 53);
    const auto *sii_54 = buffer.data(sii + 54);
    const auto *sii_55 = buffer.data(sii + 55);
    const auto *sii_77 = buffer.data(sii + 77);
    const auto *sii_78 = buffer.data(sii + 78);
    const auto *sii_79 = buffer.data(sii + 79);
    const auto *sii_80 = buffer.data(sii + 80);
    const auto *sii_81 = buffer.data(sii + 81);
    const auto *sii_82 = buffer.data(sii + 82);
    const auto *sii_83 = buffer.data(sii + 83);
    const auto *sii_84 = buffer.data(sii + 84);
    const auto *sii_87 = buffer.data(sii + 87);
    const auto *sii_89 = buffer.data(sii + 89);
    const auto *sii_90 = buffer.data(sii + 90);
    const auto *sii_93 = buffer.data(sii + 93);
    const auto *sii_94 = buffer.data(sii + 94);
    const auto *sii_96 = buffer.data(sii + 96);

    const auto *sik1_0 = buffer.data(sik1 + 0);
    const auto *sik1_3 = buffer.data(sik1 + 3);
    const auto *sik1_5 = buffer.data(sik1 + 5);
    const auto *sik1_6 = buffer.data(sik1 + 6);
    const auto *sik1_9 = buffer.data(sik1 + 9);
    const auto *sik1_10 = buffer.data(sik1 + 10);
    const auto *sik1_12 = buffer.data(sik1 + 12);
    const auto *sik1_14 = buffer.data(sik1 + 14);
    const auto *sik1_15 = buffer.data(sik1 + 15);
    const auto *sik1_17 = buffer.data(sik1 + 17);
    const auto *sik1_18 = buffer.data(sik1 + 18);
    const auto *sik1_20 = buffer.data(sik1 + 20);
    const auto *sik1_28 = buffer.data(sik1 + 28);
    const auto *sik1_35 = buffer.data(sik1 + 35);

    const auto *skh0_0 = buffer.data(skh0 + 0);
    const auto *skh0_3 = buffer.data(skh0 + 3);
    const auto *skh0_5 = buffer.data(skh0 + 5);
    const auto *skh0_6 = buffer.data(skh0 + 6);
    const auto *skh0_9 = buffer.data(skh0 + 9);
    const auto *skh0_10 = buffer.data(skh0 + 10);
    const auto *skh0_12 = buffer.data(skh0 + 12);
    const auto *skh0_14 = buffer.data(skh0 + 14);
    const auto *skh0_15 = buffer.data(skh0 + 15);
    const auto *skh0_17 = buffer.data(skh0 + 17);
    const auto *skh0_18 = buffer.data(skh0 + 18);
    const auto *skh0_19 = buffer.data(skh0 + 19);
    const auto *skh0_20 = buffer.data(skh0 + 20);
    const auto *skh0_36 = buffer.data(skh0 + 36);
    const auto *skh0_38 = buffer.data(skh0 + 38);
    const auto *skh0_39 = buffer.data(skh0 + 39);
    const auto *skh0_40 = buffer.data(skh0 + 40);
    const auto *skh0_41 = buffer.data(skh0 + 41);
    const auto *skh0_59 = buffer.data(skh0 + 59);
    const auto *skh0_60 = buffer.data(skh0 + 60);
    const auto *skh0_61 = buffer.data(skh0 + 61);
    const auto *skh0_62 = buffer.data(skh0 + 62);
    const auto *skh0_63 = buffer.data(skh0 + 63);
    const auto *skh0_66 = buffer.data(skh0 + 66);
    const auto *skh0_68 = buffer.data(skh0 + 68);
    const auto *skh0_69 = buffer.data(skh0 + 69);
    const auto *skh0_72 = buffer.data(skh0 + 72);
    const auto *skh0_73 = buffer.data(skh0 + 73);
    const auto *skh0_75 = buffer.data(skh0 + 75);

    const auto *skh1_0 = buffer.data(skh1 + 0);
    const auto *skh1_3 = buffer.data(skh1 + 3);
    const auto *skh1_5 = buffer.data(skh1 + 5);
    const auto *skh1_6 = buffer.data(skh1 + 6);
    const auto *skh1_9 = buffer.data(skh1 + 9);
    const auto *skh1_10 = buffer.data(skh1 + 10);
    const auto *skh1_12 = buffer.data(skh1 + 12);
    const auto *skh1_14 = buffer.data(skh1 + 14);
    const auto *skh1_15 = buffer.data(skh1 + 15);
    const auto *skh1_17 = buffer.data(skh1 + 17);
    const auto *skh1_18 = buffer.data(skh1 + 18);
    const auto *skh1_19 = buffer.data(skh1 + 19);
    const auto *skh1_20 = buffer.data(skh1 + 20);
    const auto *skh1_36 = buffer.data(skh1 + 36);
    const auto *skh1_38 = buffer.data(skh1 + 38);
    const auto *skh1_39 = buffer.data(skh1 + 39);
    const auto *skh1_40 = buffer.data(skh1 + 40);
    const auto *skh1_41 = buffer.data(skh1 + 41);
    const auto *skh1_59 = buffer.data(skh1 + 59);
    const auto *skh1_60 = buffer.data(skh1 + 60);
    const auto *skh1_61 = buffer.data(skh1 + 61);
    const auto *skh1_62 = buffer.data(skh1 + 62);
    const auto *skh1_63 = buffer.data(skh1 + 63);
    const auto *skh1_66 = buffer.data(skh1 + 66);
    const auto *skh1_68 = buffer.data(skh1 + 68);
    const auto *skh1_69 = buffer.data(skh1 + 69);
    const auto *skh1_72 = buffer.data(skh1 + 72);
    const auto *skh1_73 = buffer.data(skh1 + 73);
    const auto *skh1_75 = buffer.data(skh1 + 75);

    const auto *ski_0 = buffer.data(ski + 0);
    const auto *ski_2 = buffer.data(ski + 2);
    const auto *ski_3 = buffer.data(ski + 3);
    const auto *ski_5 = buffer.data(ski + 5);
    const auto *ski_6 = buffer.data(ski + 6);
    const auto *ski_9 = buffer.data(ski + 9);
    const auto *ski_10 = buffer.data(ski + 10);
    const auto *ski_12 = buffer.data(ski + 12);
    const auto *ski_14 = buffer.data(ski + 14);
    const auto *ski_15 = buffer.data(ski + 15);
    const auto *ski_17 = buffer.data(ski + 17);
    const auto *ski_18 = buffer.data(ski + 18);
    const auto *ski_20 = buffer.data(ski + 20);
    const auto *ski_21 = buffer.data(ski + 21);
    const auto *ski_22 = buffer.data(ski + 22);
    const auto *ski_23 = buffer.data(ski + 23);
    const auto *ski_24 = buffer.data(ski + 24);
    const auto *ski_25 = buffer.data(ski + 25);
    const auto *ski_26 = buffer.data(ski + 26);
    const auto *ski_27 = buffer.data(ski + 27);
    const auto *ski_28 = buffer.data(ski + 28);
    const auto *ski_30 = buffer.data(ski + 30);
    const auto *ski_31 = buffer.data(ski + 31);
    const auto *ski_33 = buffer.data(ski + 33);
    const auto *ski_34 = buffer.data(ski + 34);
    const auto *ski_37 = buffer.data(ski + 37);
    const auto *ski_38 = buffer.data(ski + 38);
    const auto *ski_42 = buffer.data(ski + 42);
    const auto *ski_49 = buffer.data(ski + 49);
    const auto *ski_50 = buffer.data(ski + 50);
    const auto *ski_51 = buffer.data(ski + 51);
    const auto *ski_52 = buffer.data(ski + 52);
    const auto *ski_53 = buffer.data(ski + 53);
    const auto *ski_54 = buffer.data(ski + 54);
    const auto *ski_55 = buffer.data(ski + 55);
    const auto *ski_56 = buffer.data(ski + 56);
    const auto *ski_58 = buffer.data(ski + 58);
    const auto *ski_59 = buffer.data(ski + 59);
    const auto *ski_61 = buffer.data(ski + 61);
    const auto *ski_62 = buffer.data(ski + 62);
    const auto *ski_65 = buffer.data(ski + 65);
    const auto *ski_66 = buffer.data(ski + 66);
    const auto *ski_70 = buffer.data(ski + 70);
    const auto *ski_77 = buffer.data(ski + 77);
    const auto *ski_78 = buffer.data(ski + 78);
    const auto *ski_79 = buffer.data(ski + 79);
    const auto *ski_80 = buffer.data(ski + 80);
    const auto *ski_81 = buffer.data(ski + 81);
    const auto *ski_82 = buffer.data(ski + 82);
    const auto *ski_83 = buffer.data(ski + 83);
    const auto *ski_84 = buffer.data(ski + 84);
    const auto *ski_86 = buffer.data(ski + 86);
    const auto *ski_87 = buffer.data(ski + 87);
    const auto *ski_89 = buffer.data(ski + 89);
    const auto *ski_90 = buffer.data(ski + 90);
    const auto *ski_93 = buffer.data(ski + 93);
    const auto *ski_94 = buffer.data(ski + 94);
    const auto *ski_96 = buffer.data(ski + 96);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, sii_0, sii_3, skh0_0, skh0_3, \
                         skh1_0, skh1_3, ski_0, ski_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sii_0[k]
                 + f_1 * skh0_0[k]
                 - f_2 * skh1_0[k]
                 + f_3 * pc_x[k] * ski_0[k];

        t_1[k] = f_3 * pc_y[k] * ski_0[k];

        t_2[k] = f_3 * pc_z[k] * ski_0[k];

        t_3[k] = f_0 * sii_3[k]
                 + f_4 * skh0_3[k]
                 - f_5 * skh1_3[k]
                 + f_3 * pc_x[k] * ski_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pc_x, pc_y, sii_5, sii_6, skh0_5, skh0_6, skh1_5, \
                         skh1_6, ski_2, ski_5, ski_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * pc_y[k] * ski_2[k];

        t_5[k] = f_0 * sii_5[k]
                 + f_4 * skh0_5[k]
                 - f_5 * skh1_5[k]
                 + f_3 * pc_x[k] * ski_5[k];

        t_6[k] = f_0 * sii_6[k]
                 + f_6 * skh0_6[k]
                 - f_7 * skh1_6[k]
                 + f_3 * pc_x[k] * ski_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pc_x, pc_y, pc_z, sii_9, skh0_9, skh1_9, ski_3, ski_5, \
                         ski_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_3 * pc_z[k] * ski_3[k];

        t_8[k] = f_3 * pc_y[k] * ski_5[k];

        t_9[k] = f_0 * sii_9[k]
                 + f_6 * skh0_9[k]
                 - f_7 * skh1_9[k]
                 + f_3 * pc_x[k] * ski_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pc_x, pc_z, sii_10, sii_12, skh0_10, skh0_12, \
                         skh1_10, skh1_12, ski_6, ski_10, ski_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * sii_10[k]
                  + f_8 * skh0_10[k]
                  - f_9 * skh1_10[k]
                  + f_3 * pc_x[k] * ski_10[k];

        t_11[k] = f_3 * pc_z[k] * ski_6[k];

        t_12[k] = f_0 * sii_12[k]
                  + f_8 * skh0_12[k]
                  - f_9 * skh1_12[k]
                  + f_3 * pc_x[k] * ski_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pc_x, pc_y, sii_14, sii_15, skh0_14, skh0_15, \
                         skh1_14, skh1_15, ski_9, ski_14, ski_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_3 * pc_y[k] * ski_9[k];

        t_14[k] = f_0 * sii_14[k]
                  + f_8 * skh0_14[k]
                  - f_9 * skh1_14[k]
                  + f_3 * pc_x[k] * ski_14[k];

        t_15[k] = f_0 * sii_15[k]
                  + f_10 * skh0_15[k]
                  - f_11 * skh1_15[k]
                  + f_3 * pc_x[k] * ski_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pc_x, pc_z, sii_17, sii_18, skh0_17, skh0_18, \
                         skh1_17, skh1_18, ski_10, ski_17, ski_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * pc_z[k] * ski_10[k];

        t_17[k] = f_0 * sii_17[k]
                  + f_10 * skh0_17[k]
                  - f_11 * skh1_17[k]
                  + f_3 * pc_x[k] * ski_17[k];

        t_18[k] = f_0 * sii_18[k]
                  + f_10 * skh0_18[k]
                  - f_11 * skh1_18[k]
                  + f_3 * pc_x[k] * ski_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pc_x, pc_y, sii_20, sii_21, sii_22, skh0_20, \
                         skh1_20, ski_14, ski_20, ski_21, ski_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * pc_y[k] * ski_14[k];

        t_20[k] = f_0 * sii_20[k]
                  + f_10 * skh0_20[k]
                  - f_11 * skh1_20[k]
                  + f_3 * pc_x[k] * ski_20[k];

        t_21[k] = f_0 * sii_21[k]
                  + f_3 * pc_x[k] * ski_21[k];

        t_22[k] = f_0 * sii_22[k]
                  + f_3 * pc_x[k] * ski_22[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, t_27, pc_x, sii_23, sii_24, sii_25, sii_26, \
                         sii_27, ski_23, ski_24, ski_25, ski_26, \
                         ski_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_0 * sii_23[k]
                  + f_3 * pc_x[k] * ski_23[k];

        t_24[k] = f_0 * sii_24[k]
                  + f_3 * pc_x[k] * ski_24[k];

        t_25[k] = f_0 * sii_25[k]
                  + f_3 * pc_x[k] * ski_25[k];

        t_26[k] = f_0 * sii_26[k]
                  + f_3 * pc_x[k] * ski_26[k];

        t_27[k] = f_0 * sii_27[k]
                  + f_3 * pc_x[k] * ski_27[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pc_y, pc_z, skh0_15, skh0_17, skh0_18, \
                         skh1_15, skh1_17, skh1_18, ski_21, ski_23, \
                         ski_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_1 * skh0_15[k]
                  - f_2 * skh1_15[k]
                  + f_3 * pc_y[k] * ski_21[k];

        t_29[k] = f_3 * pc_z[k] * ski_21[k];

        t_30[k] = f_4 * skh0_17[k]
                  - f_5 * skh1_17[k]
                  + f_3 * pc_y[k] * ski_23[k];

        t_31[k] = f_6 * skh0_18[k]
                  - f_7 * skh1_18[k]
                  + f_3 * pc_y[k] * ski_24[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pc_y, pc_z, skh0_19, skh0_20, skh1_19, \
                         skh1_20, ski_25, ski_26, ski_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_8 * skh0_19[k]
                  - f_9 * skh1_19[k]
                  + f_3 * pc_y[k] * ski_25[k];

        t_33[k] = f_10 * skh0_20[k]
                  - f_11 * skh1_20[k]
                  + f_3 * pc_y[k] * ski_26[k];

        t_34[k] = f_3 * pc_y[k] * ski_27[k];

        t_35[k] = f_1 * skh0_20[k]
                  - f_2 * skh1_20[k]
                  + f_3 * pc_z[k] * ski_27[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pb_y, pc_y, pc_z, sik0_0, sik0_3, sii_0, \
                         sii_1, sik1_0, sik1_3, ski_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = pb_y[k] * sik0_0[k]
                  - f_12 * pc_y[k] * sik1_0[k];

        t_37[k] = f_13 * sii_0[k]
                  + f_3 * pc_y[k] * ski_28[k];

        t_38[k] = f_3 * pc_z[k] * ski_28[k];

        t_39[k] = pb_y[k] * sik0_3[k]
                  + f_14 * sii_1[k]
                  - f_12 * pc_y[k] * sik1_3[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pb_y, pc_y, pc_z, sik0_5, sik0_6, sii_2, \
                         sii_3, sik1_5, sik1_6, ski_30, ski_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_13 * sii_2[k]
                  + f_3 * pc_y[k] * ski_30[k];

        t_41[k] = pb_y[k] * sik0_5[k]
                  - f_12 * pc_y[k] * sik1_5[k];

        t_42[k] = pb_y[k] * sik0_6[k]
                  + f_15 * sii_3[k]
                  - f_12 * pc_y[k] * sik1_6[k];

        t_43[k] = f_3 * pc_z[k] * ski_31[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pb_y, pc_y, pc_z, sik0_9, sik0_10, sii_5, \
                         sii_6, sik1_9, sik1_10, ski_33, ski_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_13 * sii_5[k]
                  + f_3 * pc_y[k] * ski_33[k];

        t_45[k] = pb_y[k] * sik0_9[k]
                  - f_12 * pc_y[k] * sik1_9[k];

        t_46[k] = pb_y[k] * sik0_10[k]
                  + f_16 * sii_6[k]
                  - f_12 * pc_y[k] * sik1_10[k];

        t_47[k] = f_3 * pc_z[k] * ski_34[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pb_y, pc_y, sik0_12, sik0_14, sik0_15, sii_8, \
                         sii_9, sii_10, sik1_12, sik1_14, sik1_15, \
                         ski_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = pb_y[k] * sik0_12[k]
                  + f_14 * sii_8[k]
                  - f_12 * pc_y[k] * sik1_12[k];

        t_49[k] = f_13 * sii_9[k]
                  + f_3 * pc_y[k] * ski_37[k];

        t_50[k] = pb_y[k] * sik0_14[k]
                  - f_12 * pc_y[k] * sik1_14[k];

        t_51[k] = pb_y[k] * sik0_15[k]
                  + f_17 * sii_10[k]
                  - f_12 * pc_y[k] * sik1_15[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pb_y, pc_y, pc_z, sik0_17, sik0_18, sii_12, \
                         sii_13, sii_14, sik1_17, sik1_18, ski_38, \
                         ski_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_3 * pc_z[k] * ski_38[k];

        t_53[k] = pb_y[k] * sik0_17[k]
                  + f_15 * sii_12[k]
                  - f_12 * pc_y[k] * sik1_17[k];

        t_54[k] = pb_y[k] * sik0_18[k]
                  + f_14 * sii_13[k]
                  - f_12 * pc_y[k] * sik1_18[k];

        t_55[k] = f_13 * sii_14[k]
                  + f_3 * pc_y[k] * ski_42[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pb_y, pc_x, pc_y, sik0_20, sii_49, sii_50, \
                         sii_51, sik1_20, ski_49, ski_50, ski_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_y[k] * sik0_20[k]
                  - f_12 * pc_y[k] * sik1_20[k];

        t_57[k] = f_18 * sii_49[k]
                  + f_3 * pc_x[k] * ski_49[k];

        t_58[k] = f_18 * sii_50[k]
                  + f_3 * pc_x[k] * ski_50[k];

        t_59[k] = f_18 * sii_51[k]
                  + f_3 * pc_x[k] * ski_51[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pc_x, sii_52, sii_53, sii_54, sii_55, ski_52, \
                         ski_53, ski_54, ski_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_18 * sii_52[k]
                  + f_3 * pc_x[k] * ski_52[k];

        t_61[k] = f_18 * sii_53[k]
                  + f_3 * pc_x[k] * ski_53[k];

        t_62[k] = f_18 * sii_54[k]
                  + f_3 * pc_x[k] * ski_54[k];

        t_63[k] = f_18 * sii_55[k]
                  + f_3 * pc_x[k] * ski_55[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, pc_y, pc_z, sii_21, sii_23, skh0_36, skh0_38, \
                         skh1_36, skh1_38, ski_49, ski_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_13 * sii_21[k]
                  + f_1 * skh0_36[k]
                  - f_2 * skh1_36[k]
                  + f_3 * pc_y[k] * ski_49[k];

        t_65[k] = f_3 * pc_z[k] * ski_49[k];

        t_66[k] = f_13 * sii_23[k]
                  + f_4 * skh0_38[k]
                  - f_5 * skh1_38[k]
                  + f_3 * pc_y[k] * ski_51[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, pc_y, sii_24, sii_25, sii_26, skh0_39, skh0_40, \
                         skh0_41, skh1_39, skh1_40, skh1_41, ski_52, ski_53, \
                         ski_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_13 * sii_24[k]
                  + f_6 * skh0_39[k]
                  - f_7 * skh1_39[k]
                  + f_3 * pc_y[k] * ski_52[k];

        t_68[k] = f_13 * sii_25[k]
                  + f_8 * skh0_40[k]
                  - f_9 * skh1_40[k]
                  + f_3 * pc_y[k] * ski_53[k];

        t_69[k] = f_13 * sii_26[k]
                  + f_10 * skh0_41[k]
                  - f_11 * skh1_41[k]
                  + f_3 * pc_y[k] * ski_54[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pb_y, pb_z, pc_y, pc_z, sik0_0, sik0_35, \
                         sii_27, sik1_0, sik1_35, ski_55, ski_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_13 * sii_27[k]
                  + f_3 * pc_y[k] * ski_55[k];

        t_71[k] = pb_y[k] * sik0_35[k]
                  - f_12 * pc_y[k] * sik1_35[k];

        t_72[k] = pb_z[k] * sik0_0[k]
                  - f_12 * pc_z[k] * sik1_0[k];

        t_73[k] = f_3 * pc_y[k] * ski_56[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pb_z, pc_y, pc_z, sik0_3, sik0_5, sii_0, \
                         sii_2, sik1_3, sik1_5, ski_56, ski_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_13 * sii_0[k]
                  + f_3 * pc_z[k] * ski_56[k];

        t_75[k] = pb_z[k] * sik0_3[k]
                  - f_12 * pc_z[k] * sik1_3[k];

        t_76[k] = f_3 * pc_y[k] * ski_58[k];

        t_77[k] = pb_z[k] * sik0_5[k]
                  + f_14 * sii_2[k]
                  - f_12 * pc_z[k] * sik1_5[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pb_z, pc_y, pc_z, sik0_6, sik0_9, sii_3, \
                         sii_5, sik1_6, sik1_9, ski_59, ski_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = pb_z[k] * sik0_6[k]
                  - f_12 * pc_z[k] * sik1_6[k];

        t_79[k] = f_13 * sii_3[k]
                  + f_3 * pc_z[k] * ski_59[k];

        t_80[k] = f_3 * pc_y[k] * ski_61[k];

        t_81[k] = pb_z[k] * sik0_9[k]
                  + f_15 * sii_5[k]
                  - f_12 * pc_z[k] * sik1_9[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, pb_z, pc_y, pc_z, sik0_10, sik0_12, sii_6, \
                         sii_7, sik1_10, sik1_12, ski_62, ski_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = pb_z[k] * sik0_10[k]
                  - f_12 * pc_z[k] * sik1_10[k];

        t_83[k] = f_13 * sii_6[k]
                  + f_3 * pc_z[k] * ski_62[k];

        t_84[k] = pb_z[k] * sik0_12[k]
                  + f_14 * sii_7[k]
                  - f_12 * pc_z[k] * sik1_12[k];

        t_85[k] = f_3 * pc_y[k] * ski_65[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pb_z, pc_z, sik0_14, sik0_15, sik0_17, sii_9, \
                         sii_10, sii_11, sik1_14, sik1_15, sik1_17, \
                         ski_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = pb_z[k] * sik0_14[k]
                  + f_16 * sii_9[k]
                  - f_12 * pc_z[k] * sik1_14[k];

        t_87[k] = pb_z[k] * sik0_15[k]
                  - f_12 * pc_z[k] * sik1_15[k];

        t_88[k] = f_13 * sii_10[k]
                  + f_3 * pc_z[k] * ski_66[k];

        t_89[k] = pb_z[k] * sik0_17[k]
                  + f_14 * sii_11[k]
                  - f_12 * pc_z[k] * sik1_17[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, pb_z, pc_y, pc_z, sik0_18, sik0_20, sii_12, sii_14, \
                         sik1_18, sik1_20, ski_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = pb_z[k] * sik0_18[k]
                  + f_15 * sii_12[k]
                  - f_12 * pc_z[k] * sik1_18[k];

        t_91[k] = f_3 * pc_y[k] * ski_70[k];

        t_92[k] = pb_z[k] * sik0_20[k]
                  + f_17 * sii_14[k]
                  - f_12 * pc_z[k] * sik1_20[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, t_97, pc_x, sii_77, sii_78, sii_79, sii_80, \
                         sii_81, ski_77, ski_78, ski_79, ski_80, \
                         ski_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_18 * sii_77[k]
                  + f_3 * pc_x[k] * ski_77[k];

        t_94[k] = f_18 * sii_78[k]
                  + f_3 * pc_x[k] * ski_78[k];

        t_95[k] = f_18 * sii_79[k]
                  + f_3 * pc_x[k] * ski_79[k];

        t_96[k] = f_18 * sii_80[k]
                  + f_3 * pc_x[k] * ski_80[k];

        t_97[k] = f_18 * sii_81[k]
                  + f_3 * pc_x[k] * ski_81[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, pb_z, pc_x, pc_z, sik0_28, sii_21, sii_82, \
                         sii_83, sik1_28, ski_77, ski_82, ski_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_18 * sii_82[k]
                  + f_3 * pc_x[k] * ski_82[k];

        t_99[k] = f_18 * sii_83[k]
                  + f_3 * pc_x[k] * ski_83[k];

        t_100[k] = pb_z[k] * sik0_28[k]
                   - f_12 * pc_z[k] * sik1_28[k];

        t_101[k] = f_13 * sii_21[k]
                   + f_3 * pc_z[k] * ski_77[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, pc_y, skh0_59, skh0_60, skh0_61, skh1_59, \
                         skh1_60, skh1_61, ski_79, ski_80, ski_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_4 * skh0_59[k]
                   - f_5 * skh1_59[k]
                   + f_3 * pc_y[k] * ski_79[k];

        t_103[k] = f_6 * skh0_60[k]
                   - f_7 * skh1_60[k]
                   + f_3 * pc_y[k] * ski_80[k];

        t_104[k] = f_8 * skh0_61[k]
                   - f_9 * skh1_61[k]
                   + f_3 * pc_y[k] * ski_81[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, pc_x, pc_y, pc_z, sii_27, sii_84, \
                         skh0_62, skh0_63, skh1_62, skh1_63, ski_82, ski_83, \
                         ski_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_10 * skh0_62[k]
                   - f_11 * skh1_62[k]
                   + f_3 * pc_y[k] * ski_82[k];

        t_106[k] = f_3 * pc_y[k] * ski_83[k];

        t_107[k] = f_13 * sii_27[k]
                   + f_1 * skh0_62[k]
                   - f_2 * skh1_62[k]
                   + f_3 * pc_z[k] * ski_83[k];

        t_108[k] = f_17 * sii_84[k]
                   + f_1 * skh0_63[k]
                   - f_2 * skh1_63[k]
                   + f_3 * pc_x[k] * ski_84[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pc_x, pc_y, pc_z, sii_28, sii_30, sii_87, \
                         skh0_66, skh1_66, ski_84, ski_86, ski_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_14 * sii_28[k]
                   + f_3 * pc_y[k] * ski_84[k];

        t_110[k] = f_3 * pc_z[k] * ski_84[k];

        t_111[k] = f_17 * sii_87[k]
                   + f_4 * skh0_66[k]
                   - f_5 * skh1_66[k]
                   + f_3 * pc_x[k] * ski_87[k];

        t_112[k] = f_14 * sii_30[k]
                   + f_3 * pc_y[k] * ski_86[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, pc_x, pc_z, sii_89, sii_90, skh0_68, skh0_69, \
                         skh1_68, skh1_69, ski_87, ski_89, ski_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_17 * sii_89[k]
                   + f_4 * skh0_68[k]
                   - f_5 * skh1_68[k]
                   + f_3 * pc_x[k] * ski_89[k];

        t_114[k] = f_17 * sii_90[k]
                   + f_6 * skh0_69[k]
                   - f_7 * skh1_69[k]
                   + f_3 * pc_x[k] * ski_90[k];

        t_115[k] = f_3 * pc_z[k] * ski_87[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, pc_x, pc_y, sii_33, sii_93, sii_94, skh0_72, \
                         skh0_73, skh1_72, skh1_73, ski_89, ski_93, \
                         ski_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_14 * sii_33[k]
                   + f_3 * pc_y[k] * ski_89[k];

        t_117[k] = f_17 * sii_93[k]
                   + f_6 * skh0_72[k]
                   - f_7 * skh1_72[k]
                   + f_3 * pc_x[k] * ski_93[k];

        t_118[k] = f_17 * sii_94[k]
                   + f_8 * skh0_73[k]
                   - f_9 * skh1_73[k]
                   + f_3 * pc_x[k] * ski_94[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, pc_x, pc_y, pc_z, sii_37, sii_96, skh0_75, \
                         skh1_75, ski_90, ski_93, ski_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = f_3 * pc_z[k] * ski_90[k];

        t_120[k] = f_17 * sii_96[k]
                   + f_8 * skh0_75[k]
                   - f_9 * skh1_75[k]
                   + f_3 * pc_x[k] * ski_96[k];

        t_121[k] = f_14 * sii_37[k]
                   + f_3 * pc_y[k] * ski_93[k];
    }
}

static auto
compute_prim_skk_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sik0,
                                                          const size_t sii, const size_t sik1,
                                                          const size_t skh0, const size_t skh1,
                                                          const size_t ski, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.0 / gamma;
    const auto f_2 = 3.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.0 / gamma;
    const auto f_5 = 2.0 * p / (gamma * q);
    const auto f_6 = 1.5 / gamma;
    const auto f_7 = 1.5 * p / (gamma * q);
    const auto f_8 = 1.0 / gamma;
    const auto f_9 = p / (gamma * q);
    const auto f_10 = 0.5 / gamma;
    const auto f_11 = 0.5 * p / (gamma * q);
    const auto f_12 = gamma / q;
    const auto f_13 = 0.5 / q;
    const auto f_14 = 1.0 / q;
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;

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

    const auto *sik0_39 = buffer.data(sik0 + 39);
    const auto *sik0_42 = buffer.data(sik0 + 42);
    const auto *sik0_46 = buffer.data(sik0 + 46);
    const auto *sik0_51 = buffer.data(sik0 + 51);
    const auto *sik0_64 = buffer.data(sik0 + 64);
    const auto *sik0_72 = buffer.data(sik0 + 72);
    const auto *sik0_77 = buffer.data(sik0 + 77);
    const auto *sik0_81 = buffer.data(sik0 + 81);
    const auto *sik0_84 = buffer.data(sik0 + 84);
    const auto *sik0_86 = buffer.data(sik0 + 86);
    const auto *sik0_89 = buffer.data(sik0 + 89);
    const auto *sik0_90 = buffer.data(sik0 + 90);
    const auto *sik0_92 = buffer.data(sik0 + 92);
    const auto *sik0_107 = buffer.data(sik0 + 107);

    const auto *sii_28 = buffer.data(sii + 28);
    const auto *sii_31 = buffer.data(sii + 31);
    const auto *sii_34 = buffer.data(sii + 34);
    const auto *sii_38 = buffer.data(sii + 38);
    const auto *sii_42 = buffer.data(sii + 42);
    const auto *sii_49 = buffer.data(sii + 49);
    const auto *sii_51 = buffer.data(sii + 51);
    const auto *sii_52 = buffer.data(sii + 52);
    const auto *sii_53 = buffer.data(sii + 53);
    const auto *sii_54 = buffer.data(sii + 54);
    const auto *sii_55 = buffer.data(sii + 55);
    const auto *sii_56 = buffer.data(sii + 56);
    const auto *sii_58 = buffer.data(sii + 58);
    const auto *sii_59 = buffer.data(sii + 59);
    const auto *sii_61 = buffer.data(sii + 61);
    const auto *sii_62 = buffer.data(sii + 62);
    const auto *sii_64 = buffer.data(sii + 64);
    const auto *sii_65 = buffer.data(sii + 65);
    const auto *sii_66 = buffer.data(sii + 66);
    const auto *sii_68 = buffer.data(sii + 68);
    const auto *sii_69 = buffer.data(sii + 69);
    const auto *sii_70 = buffer.data(sii + 70);
    const auto *sii_77 = buffer.data(sii + 77);
    const auto *sii_79 = buffer.data(sii + 79);
    const auto *sii_80 = buffer.data(sii + 80);
    const auto *sii_81 = buffer.data(sii + 81);
    const auto *sii_82 = buffer.data(sii + 82);
    const auto *sii_83 = buffer.data(sii + 83);
    const auto *sii_84 = buffer.data(sii + 84);
    const auto *sii_86 = buffer.data(sii + 86);
    const auto *sii_89 = buffer.data(sii + 89);
    const auto *sii_93 = buffer.data(sii + 93);
    const auto *sii_98 = buffer.data(sii + 98);
    const auto *sii_99 = buffer.data(sii + 99);
    const auto *sii_101 = buffer.data(sii + 101);
    const auto *sii_102 = buffer.data(sii + 102);
    const auto *sii_104 = buffer.data(sii + 104);
    const auto *sii_105 = buffer.data(sii + 105);
    const auto *sii_106 = buffer.data(sii + 106);
    const auto *sii_107 = buffer.data(sii + 107);
    const auto *sii_108 = buffer.data(sii + 108);
    const auto *sii_109 = buffer.data(sii + 109);
    const auto *sii_110 = buffer.data(sii + 110);
    const auto *sii_111 = buffer.data(sii + 111);
    const auto *sii_133 = buffer.data(sii + 133);
    const auto *sii_134 = buffer.data(sii + 134);
    const auto *sii_135 = buffer.data(sii + 135);
    const auto *sii_136 = buffer.data(sii + 136);
    const auto *sii_137 = buffer.data(sii + 137);
    const auto *sii_138 = buffer.data(sii + 138);
    const auto *sii_139 = buffer.data(sii + 139);
    const auto *sii_140 = buffer.data(sii + 140);
    const auto *sii_143 = buffer.data(sii + 143);
    const auto *sii_145 = buffer.data(sii + 145);
    const auto *sii_146 = buffer.data(sii + 146);
    const auto *sii_149 = buffer.data(sii + 149);
    const auto *sii_150 = buffer.data(sii + 150);
    const auto *sii_152 = buffer.data(sii + 152);
    const auto *sii_154 = buffer.data(sii + 154);
    const auto *sii_155 = buffer.data(sii + 155);
    const auto *sii_157 = buffer.data(sii + 157);
    const auto *sii_158 = buffer.data(sii + 158);
    const auto *sii_160 = buffer.data(sii + 160);
    const auto *sii_161 = buffer.data(sii + 161);
    const auto *sii_162 = buffer.data(sii + 162);
    const auto *sii_163 = buffer.data(sii + 163);
    const auto *sii_164 = buffer.data(sii + 164);
    const auto *sii_165 = buffer.data(sii + 165);
    const auto *sii_166 = buffer.data(sii + 166);
    const auto *sii_167 = buffer.data(sii + 167);
    const auto *sii_168 = buffer.data(sii + 168);
    const auto *sii_171 = buffer.data(sii + 171);
    const auto *sii_173 = buffer.data(sii + 173);
    const auto *sii_174 = buffer.data(sii + 174);
    const auto *sii_177 = buffer.data(sii + 177);
    const auto *sii_178 = buffer.data(sii + 178);
    const auto *sii_180 = buffer.data(sii + 180);
    const auto *sii_182 = buffer.data(sii + 182);
    const auto *sii_183 = buffer.data(sii + 183);
    const auto *sii_185 = buffer.data(sii + 185);
    const auto *sii_186 = buffer.data(sii + 186);

    const auto *sik1_39 = buffer.data(sik1 + 39);
    const auto *sik1_42 = buffer.data(sik1 + 42);
    const auto *sik1_46 = buffer.data(sik1 + 46);
    const auto *sik1_51 = buffer.data(sik1 + 51);
    const auto *sik1_64 = buffer.data(sik1 + 64);
    const auto *sik1_72 = buffer.data(sik1 + 72);
    const auto *sik1_77 = buffer.data(sik1 + 77);
    const auto *sik1_81 = buffer.data(sik1 + 81);
    const auto *sik1_84 = buffer.data(sik1 + 84);
    const auto *sik1_86 = buffer.data(sik1 + 86);
    const auto *sik1_89 = buffer.data(sik1 + 89);
    const auto *sik1_90 = buffer.data(sik1 + 90);
    const auto *sik1_92 = buffer.data(sik1 + 92);
    const auto *sik1_107 = buffer.data(sik1 + 107);

    const auto *skh0_77 = buffer.data(skh0 + 77);
    const auto *skh0_78 = buffer.data(skh0 + 78);
    const auto *skh0_80 = buffer.data(skh0 + 80);
    const auto *skh0_81 = buffer.data(skh0 + 81);
    const auto *skh0_82 = buffer.data(skh0 + 82);
    const auto *skh0_83 = buffer.data(skh0 + 83);
    const auto *skh0_101 = buffer.data(skh0 + 101);
    const auto *skh0_102 = buffer.data(skh0 + 102);
    const auto *skh0_103 = buffer.data(skh0 + 103);
    const auto *skh0_104 = buffer.data(skh0 + 104);
    const auto *skh0_105 = buffer.data(skh0 + 105);
    const auto *skh0_108 = buffer.data(skh0 + 108);
    const auto *skh0_110 = buffer.data(skh0 + 110);
    const auto *skh0_111 = buffer.data(skh0 + 111);
    const auto *skh0_114 = buffer.data(skh0 + 114);
    const auto *skh0_115 = buffer.data(skh0 + 115);
    const auto *skh0_117 = buffer.data(skh0 + 117);
    const auto *skh0_119 = buffer.data(skh0 + 119);
    const auto *skh0_120 = buffer.data(skh0 + 120);
    const auto *skh0_122 = buffer.data(skh0 + 122);
    const auto *skh0_123 = buffer.data(skh0 + 123);
    const auto *skh0_124 = buffer.data(skh0 + 124);
    const auto *skh0_125 = buffer.data(skh0 + 125);
    const auto *skh0_126 = buffer.data(skh0 + 126);
    const auto *skh0_129 = buffer.data(skh0 + 129);
    const auto *skh0_131 = buffer.data(skh0 + 131);
    const auto *skh0_132 = buffer.data(skh0 + 132);
    const auto *skh0_135 = buffer.data(skh0 + 135);
    const auto *skh0_136 = buffer.data(skh0 + 136);
    const auto *skh0_138 = buffer.data(skh0 + 138);
    const auto *skh0_140 = buffer.data(skh0 + 140);
    const auto *skh0_141 = buffer.data(skh0 + 141);
    const auto *skh0_143 = buffer.data(skh0 + 143);
    const auto *skh0_144 = buffer.data(skh0 + 144);

    const auto *skh1_77 = buffer.data(skh1 + 77);
    const auto *skh1_78 = buffer.data(skh1 + 78);
    const auto *skh1_80 = buffer.data(skh1 + 80);
    const auto *skh1_81 = buffer.data(skh1 + 81);
    const auto *skh1_82 = buffer.data(skh1 + 82);
    const auto *skh1_83 = buffer.data(skh1 + 83);
    const auto *skh1_101 = buffer.data(skh1 + 101);
    const auto *skh1_102 = buffer.data(skh1 + 102);
    const auto *skh1_103 = buffer.data(skh1 + 103);
    const auto *skh1_104 = buffer.data(skh1 + 104);
    const auto *skh1_105 = buffer.data(skh1 + 105);
    const auto *skh1_108 = buffer.data(skh1 + 108);
    const auto *skh1_110 = buffer.data(skh1 + 110);
    const auto *skh1_111 = buffer.data(skh1 + 111);
    const auto *skh1_114 = buffer.data(skh1 + 114);
    const auto *skh1_115 = buffer.data(skh1 + 115);
    const auto *skh1_117 = buffer.data(skh1 + 117);
    const auto *skh1_119 = buffer.data(skh1 + 119);
    const auto *skh1_120 = buffer.data(skh1 + 120);
    const auto *skh1_122 = buffer.data(skh1 + 122);
    const auto *skh1_123 = buffer.data(skh1 + 123);
    const auto *skh1_124 = buffer.data(skh1 + 124);
    const auto *skh1_125 = buffer.data(skh1 + 125);
    const auto *skh1_126 = buffer.data(skh1 + 126);
    const auto *skh1_129 = buffer.data(skh1 + 129);
    const auto *skh1_131 = buffer.data(skh1 + 131);
    const auto *skh1_132 = buffer.data(skh1 + 132);
    const auto *skh1_135 = buffer.data(skh1 + 135);
    const auto *skh1_136 = buffer.data(skh1 + 136);
    const auto *skh1_138 = buffer.data(skh1 + 138);
    const auto *skh1_140 = buffer.data(skh1 + 140);
    const auto *skh1_141 = buffer.data(skh1 + 141);
    const auto *skh1_143 = buffer.data(skh1 + 143);
    const auto *skh1_144 = buffer.data(skh1 + 144);

    const auto *ski_94 = buffer.data(ski + 94);
    const auto *ski_98 = buffer.data(ski + 98);
    const auto *ski_99 = buffer.data(ski + 99);
    const auto *ski_101 = buffer.data(ski + 101);
    const auto *ski_102 = buffer.data(ski + 102);
    const auto *ski_104 = buffer.data(ski + 104);
    const auto *ski_105 = buffer.data(ski + 105);
    const auto *ski_106 = buffer.data(ski + 106);
    const auto *ski_107 = buffer.data(ski + 107);
    const auto *ski_108 = buffer.data(ski + 108);
    const auto *ski_109 = buffer.data(ski + 109);
    const auto *ski_110 = buffer.data(ski + 110);
    const auto *ski_111 = buffer.data(ski + 111);
    const auto *ski_112 = buffer.data(ski + 112);
    const auto *ski_114 = buffer.data(ski + 114);
    const auto *ski_115 = buffer.data(ski + 115);
    const auto *ski_117 = buffer.data(ski + 117);
    const auto *ski_118 = buffer.data(ski + 118);
    const auto *ski_121 = buffer.data(ski + 121);
    const auto *ski_122 = buffer.data(ski + 122);
    const auto *ski_126 = buffer.data(ski + 126);
    const auto *ski_133 = buffer.data(ski + 133);
    const auto *ski_134 = buffer.data(ski + 134);
    const auto *ski_135 = buffer.data(ski + 135);
    const auto *ski_136 = buffer.data(ski + 136);
    const auto *ski_137 = buffer.data(ski + 137);
    const auto *ski_138 = buffer.data(ski + 138);
    const auto *ski_139 = buffer.data(ski + 139);
    const auto *ski_140 = buffer.data(ski + 140);
    const auto *ski_142 = buffer.data(ski + 142);
    const auto *ski_143 = buffer.data(ski + 143);
    const auto *ski_145 = buffer.data(ski + 145);
    const auto *ski_146 = buffer.data(ski + 146);
    const auto *ski_149 = buffer.data(ski + 149);
    const auto *ski_150 = buffer.data(ski + 150);
    const auto *ski_152 = buffer.data(ski + 152);
    const auto *ski_154 = buffer.data(ski + 154);
    const auto *ski_155 = buffer.data(ski + 155);
    const auto *ski_157 = buffer.data(ski + 157);
    const auto *ski_158 = buffer.data(ski + 158);
    const auto *ski_160 = buffer.data(ski + 160);
    const auto *ski_161 = buffer.data(ski + 161);
    const auto *ski_162 = buffer.data(ski + 162);
    const auto *ski_163 = buffer.data(ski + 163);
    const auto *ski_164 = buffer.data(ski + 164);
    const auto *ski_165 = buffer.data(ski + 165);
    const auto *ski_166 = buffer.data(ski + 166);
    const auto *ski_167 = buffer.data(ski + 167);
    const auto *ski_168 = buffer.data(ski + 168);
    const auto *ski_170 = buffer.data(ski + 170);
    const auto *ski_171 = buffer.data(ski + 171);
    const auto *ski_173 = buffer.data(ski + 173);
    const auto *ski_174 = buffer.data(ski + 174);
    const auto *ski_177 = buffer.data(ski + 177);
    const auto *ski_178 = buffer.data(ski + 178);
    const auto *ski_180 = buffer.data(ski + 180);
    const auto *ski_182 = buffer.data(ski + 182);
    const auto *ski_183 = buffer.data(ski + 183);
    const auto *ski_185 = buffer.data(ski + 185);
    const auto *ski_186 = buffer.data(ski + 186);

#pragma omp simd aligned(t_122, t_123, t_124, pc_x, pc_z, sii_98, sii_99, skh0_77, skh0_78, \
                         skh1_77, skh1_78, ski_94, ski_98, ski_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_17 * sii_98[k]
                   + f_8 * skh0_77[k]
                   - f_9 * skh1_77[k]
                   + f_3 * pc_x[k] * ski_98[k];

        t_123[k] = f_17 * sii_99[k]
                   + f_10 * skh0_78[k]
                   - f_11 * skh1_78[k]
                   + f_3 * pc_x[k] * ski_99[k];

        t_124[k] = f_3 * pc_z[k] * ski_94[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, pc_x, pc_y, sii_42, sii_101, sii_102, skh0_80, \
                         skh0_81, skh1_80, skh1_81, ski_98, ski_101, \
                         ski_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_17 * sii_101[k]
                   + f_10 * skh0_80[k]
                   - f_11 * skh1_80[k]
                   + f_3 * pc_x[k] * ski_101[k];

        t_126[k] = f_17 * sii_102[k]
                   + f_10 * skh0_81[k]
                   - f_11 * skh1_81[k]
                   + f_3 * pc_x[k] * ski_102[k];

        t_127[k] = f_14 * sii_42[k]
                   + f_3 * pc_y[k] * ski_98[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, pc_x, sii_104, sii_105, sii_106, sii_107, \
                         skh0_83, skh1_83, ski_104, ski_105, ski_106, \
                         ski_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_17 * sii_104[k]
                   + f_10 * skh0_83[k]
                   - f_11 * skh1_83[k]
                   + f_3 * pc_x[k] * ski_104[k];

        t_129[k] = f_17 * sii_105[k]
                   + f_3 * pc_x[k] * ski_105[k];

        t_130[k] = f_17 * sii_106[k]
                   + f_3 * pc_x[k] * ski_106[k];

        t_131[k] = f_17 * sii_107[k]
                   + f_3 * pc_x[k] * ski_107[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pc_x, sii_108, sii_109, sii_110, sii_111, \
                         ski_108, ski_109, ski_110, ski_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_17 * sii_108[k]
                   + f_3 * pc_x[k] * ski_108[k];

        t_133[k] = f_17 * sii_109[k]
                   + f_3 * pc_x[k] * ski_109[k];

        t_134[k] = f_17 * sii_110[k]
                   + f_3 * pc_x[k] * ski_110[k];

        t_135[k] = f_17 * sii_111[k]
                   + f_3 * pc_x[k] * ski_111[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, pc_y, pc_z, sii_49, sii_51, skh0_78, skh0_80, \
                         skh1_78, skh1_80, ski_105, ski_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_14 * sii_49[k]
                   + f_1 * skh0_78[k]
                   - f_2 * skh1_78[k]
                   + f_3 * pc_y[k] * ski_105[k];

        t_137[k] = f_3 * pc_z[k] * ski_105[k];

        t_138[k] = f_14 * sii_51[k]
                   + f_4 * skh0_80[k]
                   - f_5 * skh1_80[k]
                   + f_3 * pc_y[k] * ski_107[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, pc_y, sii_52, sii_53, sii_54, skh0_81, skh0_82, \
                         skh0_83, skh1_81, skh1_82, skh1_83, ski_108, ski_109, \
                         ski_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_14 * sii_52[k]
                   + f_6 * skh0_81[k]
                   - f_7 * skh1_81[k]
                   + f_3 * pc_y[k] * ski_108[k];

        t_140[k] = f_14 * sii_53[k]
                   + f_8 * skh0_82[k]
                   - f_9 * skh1_82[k]
                   + f_3 * pc_y[k] * ski_109[k];

        t_141[k] = f_14 * sii_54[k]
                   + f_10 * skh0_83[k]
                   - f_11 * skh1_83[k]
                   + f_3 * pc_y[k] * ski_110[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, pb_y, pc_y, pc_z, sik0_72, sii_55, \
                         sii_56, sik1_72, skh0_83, skh1_83, ski_111, \
                         ski_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_14 * sii_55[k]
                   + f_3 * pc_y[k] * ski_111[k];

        t_143[k] = f_1 * skh0_83[k]
                   - f_2 * skh1_83[k]
                   + f_3 * pc_z[k] * ski_111[k];

        t_144[k] = pb_y[k] * sik0_72[k]
                   - f_12 * pc_y[k] * sik1_72[k];

        t_145[k] = f_13 * sii_56[k]
                   + f_3 * pc_y[k] * ski_112[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pb_y, pb_z, pc_y, pc_z, sik0_39, sik0_77, \
                         sii_28, sii_58, sik1_39, sik1_77, ski_112, \
                         ski_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_13 * sii_28[k]
                   + f_3 * pc_z[k] * ski_112[k];

        t_147[k] = pb_z[k] * sik0_39[k]
                   - f_12 * pc_z[k] * sik1_39[k];

        t_148[k] = f_13 * sii_58[k]
                   + f_3 * pc_y[k] * ski_114[k];

        t_149[k] = pb_y[k] * sik0_77[k]
                   - f_12 * pc_y[k] * sik1_77[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, pb_y, pb_z, pc_y, pc_z, sik0_42, sik0_81, \
                         sii_31, sii_61, sik1_42, sik1_81, ski_115, \
                         ski_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = pb_z[k] * sik0_42[k]
                   - f_12 * pc_z[k] * sik1_42[k];

        t_151[k] = f_13 * sii_31[k]
                   + f_3 * pc_z[k] * ski_115[k];

        t_152[k] = f_13 * sii_61[k]
                   + f_3 * pc_y[k] * ski_117[k];

        t_153[k] = pb_y[k] * sik0_81[k]
                   - f_12 * pc_y[k] * sik1_81[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, pb_y, pb_z, pc_y, pc_z, sik0_46, sik0_84, \
                         sii_34, sii_64, sik1_46, sik1_84, ski_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = pb_z[k] * sik0_46[k]
                   - f_12 * pc_z[k] * sik1_46[k];

        t_155[k] = f_13 * sii_34[k]
                   + f_3 * pc_z[k] * ski_118[k];

        t_156[k] = pb_y[k] * sik0_84[k]
                   + f_14 * sii_64[k]
                   - f_12 * pc_y[k] * sik1_84[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, pb_y, pb_z, pc_y, pc_z, sik0_51, sik0_86, \
                         sii_38, sii_65, sik1_51, sik1_86, ski_121, \
                         ski_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = f_13 * sii_65[k]
                   + f_3 * pc_y[k] * ski_121[k];

        t_158[k] = pb_y[k] * sik0_86[k]
                   - f_12 * pc_y[k] * sik1_86[k];

        t_159[k] = pb_z[k] * sik0_51[k]
                   - f_12 * pc_z[k] * sik1_51[k];

        t_160[k] = f_13 * sii_38[k]
                   + f_3 * pc_z[k] * ski_122[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, t_164, pb_y, pc_y, sik0_89, sik0_90, sik0_92, \
                         sii_68, sii_69, sii_70, sik1_89, sik1_90, sik1_92, \
                         ski_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = pb_y[k] * sik0_89[k]
                   + f_15 * sii_68[k]
                   - f_12 * pc_y[k] * sik1_89[k];

        t_162[k] = pb_y[k] * sik0_90[k]
                   + f_14 * sii_69[k]
                   - f_12 * pc_y[k] * sik1_90[k];

        t_163[k] = f_13 * sii_70[k]
                   + f_3 * pc_y[k] * ski_126[k];

        t_164[k] = pb_y[k] * sik0_92[k]
                   - f_12 * pc_y[k] * sik1_92[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, pc_x, sii_133, sii_134, sii_135, \
                         sii_136, sii_137, ski_133, ski_134, ski_135, ski_136, \
                         ski_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = f_17 * sii_133[k]
                   + f_3 * pc_x[k] * ski_133[k];

        t_166[k] = f_17 * sii_134[k]
                   + f_3 * pc_x[k] * ski_134[k];

        t_167[k] = f_17 * sii_135[k]
                   + f_3 * pc_x[k] * ski_135[k];

        t_168[k] = f_17 * sii_136[k]
                   + f_3 * pc_x[k] * ski_136[k];

        t_169[k] = f_17 * sii_137[k]
                   + f_3 * pc_x[k] * ski_137[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, pb_z, pc_x, pc_z, sik0_64, sii_49, \
                         sii_138, sii_139, sik1_64, ski_133, ski_138, \
                         ski_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_17 * sii_138[k]
                   + f_3 * pc_x[k] * ski_138[k];

        t_171[k] = f_17 * sii_139[k]
                   + f_3 * pc_x[k] * ski_139[k];

        t_172[k] = pb_z[k] * sik0_64[k]
                   - f_12 * pc_z[k] * sik1_64[k];

        t_173[k] = f_13 * sii_49[k]
                   + f_3 * pc_z[k] * ski_133[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, pc_y, sii_79, sii_80, sii_81, skh0_101, \
                         skh0_102, skh0_103, skh1_101, skh1_102, skh1_103, ski_135, ski_136, \
                         ski_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_13 * sii_79[k]
                   + f_4 * skh0_101[k]
                   - f_5 * skh1_101[k]
                   + f_3 * pc_y[k] * ski_135[k];

        t_175[k] = f_13 * sii_80[k]
                   + f_6 * skh0_102[k]
                   - f_7 * skh1_102[k]
                   + f_3 * pc_y[k] * ski_136[k];

        t_176[k] = f_13 * sii_81[k]
                   + f_8 * skh0_103[k]
                   - f_9 * skh1_103[k]
                   + f_3 * pc_y[k] * ski_137[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, pb_y, pc_y, sik0_107, sii_82, sii_83, sik1_107, \
                         skh0_104, skh1_104, ski_138, ski_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = f_13 * sii_82[k]
                   + f_10 * skh0_104[k]
                   - f_11 * skh1_104[k]
                   + f_3 * pc_y[k] * ski_138[k];

        t_178[k] = f_13 * sii_83[k]
                   + f_3 * pc_y[k] * ski_139[k];

        t_179[k] = pb_y[k] * sik0_107[k]
                   - f_12 * pc_y[k] * sik1_107[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pc_x, pc_y, pc_z, sii_56, sii_140, \
                         sii_143, skh0_105, skh0_108, skh1_105, skh1_108, ski_140, \
                         ski_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_17 * sii_140[k]
                   + f_1 * skh0_105[k]
                   - f_2 * skh1_105[k]
                   + f_3 * pc_x[k] * ski_140[k];

        t_181[k] = f_3 * pc_y[k] * ski_140[k];

        t_182[k] = f_14 * sii_56[k]
                   + f_3 * pc_z[k] * ski_140[k];

        t_183[k] = f_17 * sii_143[k]
                   + f_4 * skh0_108[k]
                   - f_5 * skh1_108[k]
                   + f_3 * pc_x[k] * ski_143[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, pc_x, pc_y, sii_145, sii_146, skh0_110, \
                         skh0_111, skh1_110, skh1_111, ski_142, ski_145, \
                         ski_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_3 * pc_y[k] * ski_142[k];

        t_185[k] = f_17 * sii_145[k]
                   + f_4 * skh0_110[k]
                   - f_5 * skh1_110[k]
                   + f_3 * pc_x[k] * ski_145[k];

        t_186[k] = f_17 * sii_146[k]
                   + f_6 * skh0_111[k]
                   - f_7 * skh1_111[k]
                   + f_3 * pc_x[k] * ski_146[k];
    }

#pragma omp simd aligned(t_187, t_188, t_189, pc_x, pc_y, pc_z, sii_59, sii_149, skh0_114, \
                         skh1_114, ski_143, ski_145, ski_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_187[k] = f_14 * sii_59[k]
                   + f_3 * pc_z[k] * ski_143[k];

        t_188[k] = f_3 * pc_y[k] * ski_145[k];

        t_189[k] = f_17 * sii_149[k]
                   + f_6 * skh0_114[k]
                   - f_7 * skh1_114[k]
                   + f_3 * pc_x[k] * ski_149[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, pc_x, pc_z, sii_62, sii_150, sii_152, skh0_115, \
                         skh0_117, skh1_115, skh1_117, ski_146, ski_150, \
                         ski_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = f_17 * sii_150[k]
                   + f_8 * skh0_115[k]
                   - f_9 * skh1_115[k]
                   + f_3 * pc_x[k] * ski_150[k];

        t_191[k] = f_14 * sii_62[k]
                   + f_3 * pc_z[k] * ski_146[k];

        t_192[k] = f_17 * sii_152[k]
                   + f_8 * skh0_117[k]
                   - f_9 * skh1_117[k]
                   + f_3 * pc_x[k] * ski_152[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, pc_x, pc_y, sii_154, sii_155, skh0_119, \
                         skh0_120, skh1_119, skh1_120, ski_149, ski_154, \
                         ski_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_3 * pc_y[k] * ski_149[k];

        t_194[k] = f_17 * sii_154[k]
                   + f_8 * skh0_119[k]
                   - f_9 * skh1_119[k]
                   + f_3 * pc_x[k] * ski_154[k];

        t_195[k] = f_17 * sii_155[k]
                   + f_10 * skh0_120[k]
                   - f_11 * skh1_120[k]
                   + f_3 * pc_x[k] * ski_155[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, pc_x, pc_z, sii_66, sii_157, sii_158, skh0_122, \
                         skh0_123, skh1_122, skh1_123, ski_150, ski_157, \
                         ski_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = f_14 * sii_66[k]
                   + f_3 * pc_z[k] * ski_150[k];

        t_197[k] = f_17 * sii_157[k]
                   + f_10 * skh0_122[k]
                   - f_11 * skh1_122[k]
                   + f_3 * pc_x[k] * ski_157[k];

        t_198[k] = f_17 * sii_158[k]
                   + f_10 * skh0_123[k]
                   - f_11 * skh1_123[k]
                   + f_3 * pc_x[k] * ski_158[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, pc_x, pc_y, sii_160, sii_161, sii_162, \
                         skh0_125, skh1_125, ski_154, ski_160, ski_161, \
                         ski_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = f_3 * pc_y[k] * ski_154[k];

        t_200[k] = f_17 * sii_160[k]
                   + f_10 * skh0_125[k]
                   - f_11 * skh1_125[k]
                   + f_3 * pc_x[k] * ski_160[k];

        t_201[k] = f_17 * sii_161[k]
                   + f_3 * pc_x[k] * ski_161[k];

        t_202[k] = f_17 * sii_162[k]
                   + f_3 * pc_x[k] * ski_162[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, t_207, pc_x, sii_163, sii_164, sii_165, \
                         sii_166, sii_167, ski_163, ski_164, ski_165, ski_166, \
                         ski_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_17 * sii_163[k]
                   + f_3 * pc_x[k] * ski_163[k];

        t_204[k] = f_17 * sii_164[k]
                   + f_3 * pc_x[k] * ski_164[k];

        t_205[k] = f_17 * sii_165[k]
                   + f_3 * pc_x[k] * ski_165[k];

        t_206[k] = f_17 * sii_166[k]
                   + f_3 * pc_x[k] * ski_166[k];

        t_207[k] = f_17 * sii_167[k]
                   + f_3 * pc_x[k] * ski_167[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, pc_y, pc_z, sii_77, skh0_120, skh0_122, \
                         skh0_123, skh1_120, skh1_122, skh1_123, ski_161, ski_163, \
                         ski_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_1 * skh0_120[k]
                   - f_2 * skh1_120[k]
                   + f_3 * pc_y[k] * ski_161[k];

        t_209[k] = f_14 * sii_77[k]
                   + f_3 * pc_z[k] * ski_161[k];

        t_210[k] = f_4 * skh0_122[k]
                   - f_5 * skh1_122[k]
                   + f_3 * pc_y[k] * ski_163[k];

        t_211[k] = f_6 * skh0_123[k]
                   - f_7 * skh1_123[k]
                   + f_3 * pc_y[k] * ski_164[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, pc_y, pc_z, sii_83, skh0_124, skh0_125, \
                         skh1_124, skh1_125, ski_165, ski_166, \
                         ski_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_8 * skh0_124[k]
                   - f_9 * skh1_124[k]
                   + f_3 * pc_y[k] * ski_165[k];

        t_213[k] = f_10 * skh0_125[k]
                   - f_11 * skh1_125[k]
                   + f_3 * pc_y[k] * ski_166[k];

        t_214[k] = f_3 * pc_y[k] * ski_167[k];

        t_215[k] = f_14 * sii_83[k]
                   + f_1 * skh0_125[k]
                   - f_2 * skh1_125[k]
                   + f_3 * pc_z[k] * ski_167[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, pc_x, pc_y, pc_z, sii_84, sii_168, \
                         sii_171, skh0_126, skh0_129, skh1_126, skh1_129, ski_168, \
                         ski_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = f_16 * sii_168[k]
                   + f_1 * skh0_126[k]
                   - f_2 * skh1_126[k]
                   + f_3 * pc_x[k] * ski_168[k];

        t_217[k] = f_15 * sii_84[k]
                   + f_3 * pc_y[k] * ski_168[k];

        t_218[k] = f_3 * pc_z[k] * ski_168[k];

        t_219[k] = f_16 * sii_171[k]
                   + f_4 * skh0_129[k]
                   - f_5 * skh1_129[k]
                   + f_3 * pc_x[k] * ski_171[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, pc_x, pc_y, sii_86, sii_173, sii_174, skh0_131, \
                         skh0_132, skh1_131, skh1_132, ski_170, ski_173, \
                         ski_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = f_15 * sii_86[k]
                   + f_3 * pc_y[k] * ski_170[k];

        t_221[k] = f_16 * sii_173[k]
                   + f_4 * skh0_131[k]
                   - f_5 * skh1_131[k]
                   + f_3 * pc_x[k] * ski_173[k];

        t_222[k] = f_16 * sii_174[k]
                   + f_6 * skh0_132[k]
                   - f_7 * skh1_132[k]
                   + f_3 * pc_x[k] * ski_174[k];
    }

#pragma omp simd aligned(t_223, t_224, t_225, pc_x, pc_y, pc_z, sii_89, sii_177, skh0_135, \
                         skh1_135, ski_171, ski_173, ski_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = f_3 * pc_z[k] * ski_171[k];

        t_224[k] = f_15 * sii_89[k]
                   + f_3 * pc_y[k] * ski_173[k];

        t_225[k] = f_16 * sii_177[k]
                   + f_6 * skh0_135[k]
                   - f_7 * skh1_135[k]
                   + f_3 * pc_x[k] * ski_177[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, pc_x, pc_z, sii_178, sii_180, skh0_136, \
                         skh0_138, skh1_136, skh1_138, ski_174, ski_178, \
                         ski_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = f_16 * sii_178[k]
                   + f_8 * skh0_136[k]
                   - f_9 * skh1_136[k]
                   + f_3 * pc_x[k] * ski_178[k];

        t_227[k] = f_3 * pc_z[k] * ski_174[k];

        t_228[k] = f_16 * sii_180[k]
                   + f_8 * skh0_138[k]
                   - f_9 * skh1_138[k]
                   + f_3 * pc_x[k] * ski_180[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, pc_x, pc_y, sii_93, sii_182, sii_183, skh0_140, \
                         skh0_141, skh1_140, skh1_141, ski_177, ski_182, \
                         ski_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = f_15 * sii_93[k]
                   + f_3 * pc_y[k] * ski_177[k];

        t_230[k] = f_16 * sii_182[k]
                   + f_8 * skh0_140[k]
                   - f_9 * skh1_140[k]
                   + f_3 * pc_x[k] * ski_182[k];

        t_231[k] = f_16 * sii_183[k]
                   + f_10 * skh0_141[k]
                   - f_11 * skh1_141[k]
                   + f_3 * pc_x[k] * ski_183[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, pc_x, pc_z, sii_185, sii_186, skh0_143, \
                         skh0_144, skh1_143, skh1_144, ski_178, ski_185, \
                         ski_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = f_3 * pc_z[k] * ski_178[k];

        t_233[k] = f_16 * sii_185[k]
                   + f_10 * skh0_143[k]
                   - f_11 * skh1_143[k]
                   + f_3 * pc_x[k] * ski_185[k];

        t_234[k] = f_16 * sii_186[k]
                   + f_10 * skh0_144[k]
                   - f_11 * skh1_144[k]
                   + f_3 * pc_x[k] * ski_186[k];
    }
}

static auto
compute_prim_skk_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sik0,
                                                          const size_t sii, const size_t sik1,
                                                          const size_t skh0, const size_t skh1,
                                                          const size_t ski, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.0 / gamma;
    const auto f_2 = 3.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.0 / gamma;
    const auto f_5 = 2.0 * p / (gamma * q);
    const auto f_6 = 1.5 / gamma;
    const auto f_7 = 1.5 * p / (gamma * q);
    const auto f_8 = 1.0 / gamma;
    const auto f_9 = p / (gamma * q);
    const auto f_10 = 0.5 / gamma;
    const auto f_11 = 0.5 * p / (gamma * q);
    const auto f_12 = gamma / q;
    const auto f_13 = 0.5 / q;
    const auto f_14 = 1.0 / q;
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sik0_108 = buffer.data(sik0 + 108);
    const auto *sik0_111 = buffer.data(sik0 + 111);
    const auto *sik0_114 = buffer.data(sik0 + 114);
    const auto *sik0_118 = buffer.data(sik0 + 118);
    const auto *sik0_120 = buffer.data(sik0 + 120);
    const auto *sik0_123 = buffer.data(sik0 + 123);
    const auto *sik0_125 = buffer.data(sik0 + 125);
    const auto *sik0_126 = buffer.data(sik0 + 126);
    const auto *sik0_136 = buffer.data(sik0 + 136);
    const auto *sik0_180 = buffer.data(sik0 + 180);
    const auto *sik0_183 = buffer.data(sik0 + 183);
    const auto *sik0_185 = buffer.data(sik0 + 185);
    const auto *sik0_186 = buffer.data(sik0 + 186);
    const auto *sik0_189 = buffer.data(sik0 + 189);
    const auto *sik0_190 = buffer.data(sik0 + 190);
    const auto *sik0_192 = buffer.data(sik0 + 192);
    const auto *sik0_194 = buffer.data(sik0 + 194);
    const auto *sik0_195 = buffer.data(sik0 + 195);
    const auto *sik0_197 = buffer.data(sik0 + 197);
    const auto *sik0_198 = buffer.data(sik0 + 198);
    const auto *sik0_200 = buffer.data(sik0 + 200);
    const auto *sik0_215 = buffer.data(sik0 + 215);

    const auto *sii_84 = buffer.data(sii + 84);
    const auto *sii_87 = buffer.data(sii + 87);
    const auto *sii_90 = buffer.data(sii + 90);
    const auto *sii_91 = buffer.data(sii + 91);
    const auto *sii_94 = buffer.data(sii + 94);
    const auto *sii_95 = buffer.data(sii + 95);
    const auto *sii_96 = buffer.data(sii + 96);
    const auto *sii_98 = buffer.data(sii + 98);
    const auto *sii_105 = buffer.data(sii + 105);
    const auto *sii_107 = buffer.data(sii + 107);
    const auto *sii_108 = buffer.data(sii + 108);
    const auto *sii_109 = buffer.data(sii + 109);
    const auto *sii_110 = buffer.data(sii + 110);
    const auto *sii_111 = buffer.data(sii + 111);
    const auto *sii_112 = buffer.data(sii + 112);
    const auto *sii_114 = buffer.data(sii + 114);
    const auto *sii_115 = buffer.data(sii + 115);
    const auto *sii_117 = buffer.data(sii + 117);
    const auto *sii_118 = buffer.data(sii + 118);
    const auto *sii_121 = buffer.data(sii + 121);
    const auto *sii_122 = buffer.data(sii + 122);
    const auto *sii_126 = buffer.data(sii + 126);
    const auto *sii_133 = buffer.data(sii + 133);
    const auto *sii_135 = buffer.data(sii + 135);
    const auto *sii_136 = buffer.data(sii + 136);
    const auto *sii_137 = buffer.data(sii + 137);
    const auto *sii_138 = buffer.data(sii + 138);
    const auto *sii_139 = buffer.data(sii + 139);
    const auto *sii_140 = buffer.data(sii + 140);
    const auto *sii_141 = buffer.data(sii + 141);
    const auto *sii_142 = buffer.data(sii + 142);
    const auto *sii_143 = buffer.data(sii + 143);
    const auto *sii_145 = buffer.data(sii + 145);
    const auto *sii_146 = buffer.data(sii + 146);
    const auto *sii_148 = buffer.data(sii + 148);
    const auto *sii_149 = buffer.data(sii + 149);
    const auto *sii_150 = buffer.data(sii + 150);
    const auto *sii_152 = buffer.data(sii + 152);
    const auto *sii_153 = buffer.data(sii + 153);
    const auto *sii_154 = buffer.data(sii + 154);
    const auto *sii_161 = buffer.data(sii + 161);
    const auto *sii_163 = buffer.data(sii + 163);
    const auto *sii_164 = buffer.data(sii + 164);
    const auto *sii_165 = buffer.data(sii + 165);
    const auto *sii_166 = buffer.data(sii + 166);
    const auto *sii_167 = buffer.data(sii + 167);
    const auto *sii_188 = buffer.data(sii + 188);
    const auto *sii_189 = buffer.data(sii + 189);
    const auto *sii_190 = buffer.data(sii + 190);
    const auto *sii_191 = buffer.data(sii + 191);
    const auto *sii_192 = buffer.data(sii + 192);
    const auto *sii_193 = buffer.data(sii + 193);
    const auto *sii_194 = buffer.data(sii + 194);
    const auto *sii_195 = buffer.data(sii + 195);
    const auto *sii_201 = buffer.data(sii + 201);
    const auto *sii_205 = buffer.data(sii + 205);
    const auto *sii_210 = buffer.data(sii + 210);
    const auto *sii_216 = buffer.data(sii + 216);
    const auto *sii_217 = buffer.data(sii + 217);
    const auto *sii_218 = buffer.data(sii + 218);
    const auto *sii_219 = buffer.data(sii + 219);
    const auto *sii_220 = buffer.data(sii + 220);
    const auto *sii_221 = buffer.data(sii + 221);
    const auto *sii_222 = buffer.data(sii + 222);
    const auto *sii_223 = buffer.data(sii + 223);
    const auto *sii_245 = buffer.data(sii + 245);
    const auto *sii_246 = buffer.data(sii + 246);
    const auto *sii_247 = buffer.data(sii + 247);
    const auto *sii_248 = buffer.data(sii + 248);
    const auto *sii_249 = buffer.data(sii + 249);
    const auto *sii_250 = buffer.data(sii + 250);
    const auto *sii_251 = buffer.data(sii + 251);
    const auto *sii_252 = buffer.data(sii + 252);
    const auto *sii_255 = buffer.data(sii + 255);
    const auto *sii_257 = buffer.data(sii + 257);
    const auto *sii_258 = buffer.data(sii + 258);
    const auto *sii_261 = buffer.data(sii + 261);
    const auto *sii_262 = buffer.data(sii + 262);
    const auto *sii_264 = buffer.data(sii + 264);
    const auto *sii_266 = buffer.data(sii + 266);
    const auto *sii_267 = buffer.data(sii + 267);
    const auto *sii_269 = buffer.data(sii + 269);
    const auto *sii_270 = buffer.data(sii + 270);
    const auto *sii_272 = buffer.data(sii + 272);
    const auto *sii_273 = buffer.data(sii + 273);
    const auto *sii_274 = buffer.data(sii + 274);

    const auto *sik1_108 = buffer.data(sik1 + 108);
    const auto *sik1_111 = buffer.data(sik1 + 111);
    const auto *sik1_114 = buffer.data(sik1 + 114);
    const auto *sik1_118 = buffer.data(sik1 + 118);
    const auto *sik1_120 = buffer.data(sik1 + 120);
    const auto *sik1_123 = buffer.data(sik1 + 123);
    const auto *sik1_125 = buffer.data(sik1 + 125);
    const auto *sik1_126 = buffer.data(sik1 + 126);
    const auto *sik1_136 = buffer.data(sik1 + 136);
    const auto *sik1_180 = buffer.data(sik1 + 180);
    const auto *sik1_183 = buffer.data(sik1 + 183);
    const auto *sik1_185 = buffer.data(sik1 + 185);
    const auto *sik1_186 = buffer.data(sik1 + 186);
    const auto *sik1_189 = buffer.data(sik1 + 189);
    const auto *sik1_190 = buffer.data(sik1 + 190);
    const auto *sik1_192 = buffer.data(sik1 + 192);
    const auto *sik1_194 = buffer.data(sik1 + 194);
    const auto *sik1_195 = buffer.data(sik1 + 195);
    const auto *sik1_197 = buffer.data(sik1 + 197);
    const auto *sik1_198 = buffer.data(sik1 + 198);
    const auto *sik1_200 = buffer.data(sik1 + 200);
    const auto *sik1_215 = buffer.data(sik1 + 215);

    const auto *skh0_141 = buffer.data(skh0 + 141);
    const auto *skh0_143 = buffer.data(skh0 + 143);
    const auto *skh0_144 = buffer.data(skh0 + 144);
    const auto *skh0_145 = buffer.data(skh0 + 145);
    const auto *skh0_146 = buffer.data(skh0 + 146);
    const auto *skh0_152 = buffer.data(skh0 + 152);
    const auto *skh0_156 = buffer.data(skh0 + 156);
    const auto *skh0_161 = buffer.data(skh0 + 161);
    const auto *skh0_164 = buffer.data(skh0 + 164);
    const auto *skh0_165 = buffer.data(skh0 + 165);
    const auto *skh0_166 = buffer.data(skh0 + 166);
    const auto *skh0_167 = buffer.data(skh0 + 167);
    const auto *skh0_183 = buffer.data(skh0 + 183);
    const auto *skh0_185 = buffer.data(skh0 + 185);
    const auto *skh0_186 = buffer.data(skh0 + 186);
    const auto *skh0_187 = buffer.data(skh0 + 187);
    const auto *skh0_188 = buffer.data(skh0 + 188);
    const auto *skh0_189 = buffer.data(skh0 + 189);
    const auto *skh0_192 = buffer.data(skh0 + 192);
    const auto *skh0_194 = buffer.data(skh0 + 194);
    const auto *skh0_195 = buffer.data(skh0 + 195);
    const auto *skh0_198 = buffer.data(skh0 + 198);
    const auto *skh0_199 = buffer.data(skh0 + 199);
    const auto *skh0_201 = buffer.data(skh0 + 201);
    const auto *skh0_203 = buffer.data(skh0 + 203);
    const auto *skh0_204 = buffer.data(skh0 + 204);
    const auto *skh0_206 = buffer.data(skh0 + 206);
    const auto *skh0_207 = buffer.data(skh0 + 207);
    const auto *skh0_209 = buffer.data(skh0 + 209);

    const auto *skh1_141 = buffer.data(skh1 + 141);
    const auto *skh1_143 = buffer.data(skh1 + 143);
    const auto *skh1_144 = buffer.data(skh1 + 144);
    const auto *skh1_145 = buffer.data(skh1 + 145);
    const auto *skh1_146 = buffer.data(skh1 + 146);
    const auto *skh1_152 = buffer.data(skh1 + 152);
    const auto *skh1_156 = buffer.data(skh1 + 156);
    const auto *skh1_161 = buffer.data(skh1 + 161);
    const auto *skh1_164 = buffer.data(skh1 + 164);
    const auto *skh1_165 = buffer.data(skh1 + 165);
    const auto *skh1_166 = buffer.data(skh1 + 166);
    const auto *skh1_167 = buffer.data(skh1 + 167);
    const auto *skh1_183 = buffer.data(skh1 + 183);
    const auto *skh1_185 = buffer.data(skh1 + 185);
    const auto *skh1_186 = buffer.data(skh1 + 186);
    const auto *skh1_187 = buffer.data(skh1 + 187);
    const auto *skh1_188 = buffer.data(skh1 + 188);
    const auto *skh1_189 = buffer.data(skh1 + 189);
    const auto *skh1_192 = buffer.data(skh1 + 192);
    const auto *skh1_194 = buffer.data(skh1 + 194);
    const auto *skh1_195 = buffer.data(skh1 + 195);
    const auto *skh1_198 = buffer.data(skh1 + 198);
    const auto *skh1_199 = buffer.data(skh1 + 199);
    const auto *skh1_201 = buffer.data(skh1 + 201);
    const auto *skh1_203 = buffer.data(skh1 + 203);
    const auto *skh1_204 = buffer.data(skh1 + 204);
    const auto *skh1_206 = buffer.data(skh1 + 206);
    const auto *skh1_207 = buffer.data(skh1 + 207);
    const auto *skh1_209 = buffer.data(skh1 + 209);

    const auto *ski_182 = buffer.data(ski + 182);
    const auto *ski_188 = buffer.data(ski + 188);
    const auto *ski_189 = buffer.data(ski + 189);
    const auto *ski_190 = buffer.data(ski + 190);
    const auto *ski_191 = buffer.data(ski + 191);
    const auto *ski_192 = buffer.data(ski + 192);
    const auto *ski_193 = buffer.data(ski + 193);
    const auto *ski_194 = buffer.data(ski + 194);
    const auto *ski_195 = buffer.data(ski + 195);
    const auto *ski_196 = buffer.data(ski + 196);
    const auto *ski_198 = buffer.data(ski + 198);
    const auto *ski_199 = buffer.data(ski + 199);
    const auto *ski_201 = buffer.data(ski + 201);
    const auto *ski_202 = buffer.data(ski + 202);
    const auto *ski_205 = buffer.data(ski + 205);
    const auto *ski_206 = buffer.data(ski + 206);
    const auto *ski_210 = buffer.data(ski + 210);
    const auto *ski_216 = buffer.data(ski + 216);
    const auto *ski_217 = buffer.data(ski + 217);
    const auto *ski_218 = buffer.data(ski + 218);
    const auto *ski_219 = buffer.data(ski + 219);
    const auto *ski_220 = buffer.data(ski + 220);
    const auto *ski_221 = buffer.data(ski + 221);
    const auto *ski_222 = buffer.data(ski + 222);
    const auto *ski_223 = buffer.data(ski + 223);
    const auto *ski_224 = buffer.data(ski + 224);
    const auto *ski_226 = buffer.data(ski + 226);
    const auto *ski_227 = buffer.data(ski + 227);
    const auto *ski_229 = buffer.data(ski + 229);
    const auto *ski_230 = buffer.data(ski + 230);
    const auto *ski_233 = buffer.data(ski + 233);
    const auto *ski_234 = buffer.data(ski + 234);
    const auto *ski_238 = buffer.data(ski + 238);
    const auto *ski_245 = buffer.data(ski + 245);
    const auto *ski_246 = buffer.data(ski + 246);
    const auto *ski_247 = buffer.data(ski + 247);
    const auto *ski_248 = buffer.data(ski + 248);
    const auto *ski_249 = buffer.data(ski + 249);
    const auto *ski_250 = buffer.data(ski + 250);
    const auto *ski_251 = buffer.data(ski + 251);
    const auto *ski_252 = buffer.data(ski + 252);
    const auto *ski_254 = buffer.data(ski + 254);
    const auto *ski_255 = buffer.data(ski + 255);
    const auto *ski_257 = buffer.data(ski + 257);
    const auto *ski_258 = buffer.data(ski + 258);
    const auto *ski_261 = buffer.data(ski + 261);
    const auto *ski_262 = buffer.data(ski + 262);
    const auto *ski_264 = buffer.data(ski + 264);
    const auto *ski_266 = buffer.data(ski + 266);
    const auto *ski_267 = buffer.data(ski + 267);
    const auto *ski_269 = buffer.data(ski + 269);
    const auto *ski_270 = buffer.data(ski + 270);
    const auto *ski_272 = buffer.data(ski + 272);
    const auto *ski_273 = buffer.data(ski + 273);
    const auto *ski_274 = buffer.data(ski + 274);

#pragma omp simd aligned(t_235, t_236, t_237, t_238, pc_x, pc_y, sii_98, sii_188, sii_189, \
                         sii_190, skh0_146, skh1_146, ski_182, ski_188, ski_189, \
                         ski_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = f_15 * sii_98[k]
                   + f_3 * pc_y[k] * ski_182[k];

        t_236[k] = f_16 * sii_188[k]
                   + f_10 * skh0_146[k]
                   - f_11 * skh1_146[k]
                   + f_3 * pc_x[k] * ski_188[k];

        t_237[k] = f_16 * sii_189[k]
                   + f_3 * pc_x[k] * ski_189[k];

        t_238[k] = f_16 * sii_190[k]
                   + f_3 * pc_x[k] * ski_190[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, t_243, pc_x, sii_191, sii_192, sii_193, \
                         sii_194, sii_195, ski_191, ski_192, ski_193, ski_194, \
                         ski_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_16 * sii_191[k]
                   + f_3 * pc_x[k] * ski_191[k];

        t_240[k] = f_16 * sii_192[k]
                   + f_3 * pc_x[k] * ski_192[k];

        t_241[k] = f_16 * sii_193[k]
                   + f_3 * pc_x[k] * ski_193[k];

        t_242[k] = f_16 * sii_194[k]
                   + f_3 * pc_x[k] * ski_194[k];

        t_243[k] = f_16 * sii_195[k]
                   + f_3 * pc_x[k] * ski_195[k];
    }

#pragma omp simd aligned(t_244, t_245, t_246, pc_y, pc_z, sii_105, sii_107, skh0_141, \
                         skh0_143, skh1_141, skh1_143, ski_189, \
                         ski_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_244[k] = f_15 * sii_105[k]
                   + f_1 * skh0_141[k]
                   - f_2 * skh1_141[k]
                   + f_3 * pc_y[k] * ski_189[k];

        t_245[k] = f_3 * pc_z[k] * ski_189[k];

        t_246[k] = f_15 * sii_107[k]
                   + f_4 * skh0_143[k]
                   - f_5 * skh1_143[k]
                   + f_3 * pc_y[k] * ski_191[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, pc_y, sii_108, sii_109, sii_110, skh0_144, \
                         skh0_145, skh0_146, skh1_144, skh1_145, skh1_146, ski_192, ski_193, \
                         ski_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = f_15 * sii_108[k]
                   + f_6 * skh0_144[k]
                   - f_7 * skh1_144[k]
                   + f_3 * pc_y[k] * ski_192[k];

        t_248[k] = f_15 * sii_109[k]
                   + f_8 * skh0_145[k]
                   - f_9 * skh1_145[k]
                   + f_3 * pc_y[k] * ski_193[k];

        t_249[k] = f_15 * sii_110[k]
                   + f_10 * skh0_146[k]
                   - f_11 * skh1_146[k]
                   + f_3 * pc_y[k] * ski_194[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, pb_z, pc_y, pc_z, sik0_108, sii_111, \
                         sii_112, sik1_108, skh0_146, skh1_146, ski_195, \
                         ski_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = f_15 * sii_111[k]
                   + f_3 * pc_y[k] * ski_195[k];

        t_251[k] = f_1 * skh0_146[k]
                   - f_2 * skh1_146[k]
                   + f_3 * pc_z[k] * ski_195[k];

        t_252[k] = pb_z[k] * sik0_108[k]
                   - f_12 * pc_z[k] * sik1_108[k];

        t_253[k] = f_14 * sii_112[k]
                   + f_3 * pc_y[k] * ski_196[k];
    }

#pragma omp simd aligned(t_254, t_255, t_256, pb_z, pc_y, pc_z, sik0_111, sii_84, sii_114, \
                         sik1_111, ski_196, ski_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_254[k] = f_13 * sii_84[k]
                   + f_3 * pc_z[k] * ski_196[k];

        t_255[k] = pb_z[k] * sik0_111[k]
                   - f_12 * pc_z[k] * sik1_111[k];

        t_256[k] = f_14 * sii_114[k]
                   + f_3 * pc_y[k] * ski_198[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, pb_z, pc_x, pc_z, sik0_114, sii_87, sii_201, \
                         sik1_114, skh0_152, skh1_152, ski_199, \
                         ski_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_16 * sii_201[k]
                   + f_4 * skh0_152[k]
                   - f_5 * skh1_152[k]
                   + f_3 * pc_x[k] * ski_201[k];

        t_258[k] = pb_z[k] * sik0_114[k]
                   - f_12 * pc_z[k] * sik1_114[k];

        t_259[k] = f_13 * sii_87[k]
                   + f_3 * pc_z[k] * ski_199[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, pb_z, pc_x, pc_y, pc_z, sik0_118, sii_117, \
                         sii_205, sik1_118, skh0_156, skh1_156, ski_201, \
                         ski_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_14 * sii_117[k]
                   + f_3 * pc_y[k] * ski_201[k];

        t_261[k] = f_16 * sii_205[k]
                   + f_6 * skh0_156[k]
                   - f_7 * skh1_156[k]
                   + f_3 * pc_x[k] * ski_205[k];

        t_262[k] = pb_z[k] * sik0_118[k]
                   - f_12 * pc_z[k] * sik1_118[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, pb_z, pc_y, pc_z, sik0_120, sii_90, sii_91, \
                         sii_121, sik1_120, ski_202, ski_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = f_13 * sii_90[k]
                   + f_3 * pc_z[k] * ski_202[k];

        t_264[k] = pb_z[k] * sik0_120[k]
                   + f_14 * sii_91[k]
                   - f_12 * pc_z[k] * sik1_120[k];

        t_265[k] = f_14 * sii_121[k]
                   + f_3 * pc_y[k] * ski_205[k];
    }

#pragma omp simd aligned(t_266, t_267, t_268, pb_z, pc_x, pc_z, sik0_123, sii_94, sii_210, \
                         sik1_123, skh0_161, skh1_161, ski_206, \
                         ski_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_266[k] = f_16 * sii_210[k]
                   + f_8 * skh0_161[k]
                   - f_9 * skh1_161[k]
                   + f_3 * pc_x[k] * ski_210[k];

        t_267[k] = pb_z[k] * sik0_123[k]
                   - f_12 * pc_z[k] * sik1_123[k];

        t_268[k] = f_13 * sii_94[k]
                   + f_3 * pc_z[k] * ski_206[k];
    }

#pragma omp simd aligned(t_269, t_270, t_271, pb_z, pc_y, pc_z, sik0_125, sik0_126, sii_95, \
                         sii_96, sii_126, sik1_125, sik1_126, ski_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_269[k] = pb_z[k] * sik0_125[k]
                   + f_14 * sii_95[k]
                   - f_12 * pc_z[k] * sik1_125[k];

        t_270[k] = pb_z[k] * sik0_126[k]
                   + f_15 * sii_96[k]
                   - f_12 * pc_z[k] * sik1_126[k];

        t_271[k] = f_14 * sii_126[k]
                   + f_3 * pc_y[k] * ski_210[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, pc_x, sii_216, sii_217, sii_218, sii_219, \
                         skh0_167, skh1_167, ski_216, ski_217, ski_218, \
                         ski_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_16 * sii_216[k]
                   + f_10 * skh0_167[k]
                   - f_11 * skh1_167[k]
                   + f_3 * pc_x[k] * ski_216[k];

        t_273[k] = f_16 * sii_217[k]
                   + f_3 * pc_x[k] * ski_217[k];

        t_274[k] = f_16 * sii_218[k]
                   + f_3 * pc_x[k] * ski_218[k];

        t_275[k] = f_16 * sii_219[k]
                   + f_3 * pc_x[k] * ski_219[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pc_x, sii_220, sii_221, sii_222, sii_223, \
                         ski_220, ski_221, ski_222, ski_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_16 * sii_220[k]
                   + f_3 * pc_x[k] * ski_220[k];

        t_277[k] = f_16 * sii_221[k]
                   + f_3 * pc_x[k] * ski_221[k];

        t_278[k] = f_16 * sii_222[k]
                   + f_3 * pc_x[k] * ski_222[k];

        t_279[k] = f_16 * sii_223[k]
                   + f_3 * pc_x[k] * ski_223[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, pb_z, pc_y, pc_z, sik0_136, sii_105, sii_135, \
                         sik1_136, skh0_164, skh1_164, ski_217, \
                         ski_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = pb_z[k] * sik0_136[k]
                   - f_12 * pc_z[k] * sik1_136[k];

        t_281[k] = f_13 * sii_105[k]
                   + f_3 * pc_z[k] * ski_217[k];

        t_282[k] = f_14 * sii_135[k]
                   + f_4 * skh0_164[k]
                   - f_5 * skh1_164[k]
                   + f_3 * pc_y[k] * ski_219[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, pc_y, sii_136, sii_137, sii_138, skh0_165, \
                         skh0_166, skh0_167, skh1_165, skh1_166, skh1_167, ski_220, ski_221, \
                         ski_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_14 * sii_136[k]
                   + f_6 * skh0_165[k]
                   - f_7 * skh1_165[k]
                   + f_3 * pc_y[k] * ski_220[k];

        t_284[k] = f_14 * sii_137[k]
                   + f_8 * skh0_166[k]
                   - f_9 * skh1_166[k]
                   + f_3 * pc_y[k] * ski_221[k];

        t_285[k] = f_14 * sii_138[k]
                   + f_10 * skh0_167[k]
                   - f_11 * skh1_167[k]
                   + f_3 * pc_y[k] * ski_222[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, pb_y, pc_y, pc_z, sik0_180, sii_111, \
                         sii_139, sii_140, sik1_180, skh0_167, skh1_167, ski_223, \
                         ski_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = f_14 * sii_139[k]
                   + f_3 * pc_y[k] * ski_223[k];

        t_287[k] = f_13 * sii_111[k]
                   + f_1 * skh0_167[k]
                   - f_2 * skh1_167[k]
                   + f_3 * pc_z[k] * ski_223[k];

        t_288[k] = pb_y[k] * sik0_180[k]
                   - f_12 * pc_y[k] * sik1_180[k];

        t_289[k] = f_13 * sii_140[k]
                   + f_3 * pc_y[k] * ski_224[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, pb_y, pc_y, pc_z, sik0_183, sik0_185, \
                         sii_112, sii_141, sii_142, sik1_183, sik1_185, ski_224, \
                         ski_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = f_14 * sii_112[k]
                   + f_3 * pc_z[k] * ski_224[k];

        t_291[k] = pb_y[k] * sik0_183[k]
                   + f_14 * sii_141[k]
                   - f_12 * pc_y[k] * sik1_183[k];

        t_292[k] = f_13 * sii_142[k]
                   + f_3 * pc_y[k] * ski_226[k];

        t_293[k] = pb_y[k] * sik0_185[k]
                   - f_12 * pc_y[k] * sik1_185[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, pb_y, pc_y, pc_z, sik0_186, sik0_189, \
                         sii_115, sii_143, sii_145, sik1_186, sik1_189, ski_227, \
                         ski_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = pb_y[k] * sik0_186[k]
                   + f_15 * sii_143[k]
                   - f_12 * pc_y[k] * sik1_186[k];

        t_295[k] = f_14 * sii_115[k]
                   + f_3 * pc_z[k] * ski_227[k];

        t_296[k] = f_13 * sii_145[k]
                   + f_3 * pc_y[k] * ski_229[k];

        t_297[k] = pb_y[k] * sik0_189[k]
                   - f_12 * pc_y[k] * sik1_189[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, pb_y, pc_y, pc_z, sik0_190, sik0_192, sii_118, \
                         sii_146, sii_148, sik1_190, sik1_192, \
                         ski_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = pb_y[k] * sik0_190[k]
                   + f_16 * sii_146[k]
                   - f_12 * pc_y[k] * sik1_190[k];

        t_299[k] = f_14 * sii_118[k]
                   + f_3 * pc_z[k] * ski_230[k];

        t_300[k] = pb_y[k] * sik0_192[k]
                   + f_14 * sii_148[k]
                   - f_12 * pc_y[k] * sik1_192[k];
    }

#pragma omp simd aligned(t_301, t_302, t_303, t_304, pb_y, pc_y, pc_z, sik0_194, sik0_195, \
                         sii_122, sii_149, sii_150, sik1_194, sik1_195, ski_233, \
                         ski_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_301[k] = f_13 * sii_149[k]
                   + f_3 * pc_y[k] * ski_233[k];

        t_302[k] = pb_y[k] * sik0_194[k]
                   - f_12 * pc_y[k] * sik1_194[k];

        t_303[k] = pb_y[k] * sik0_195[k]
                   + f_17 * sii_150[k]
                   - f_12 * pc_y[k] * sik1_195[k];

        t_304[k] = f_14 * sii_122[k]
                   + f_3 * pc_z[k] * ski_234[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, pb_y, pc_y, sik0_197, sik0_198, sik0_200, \
                         sii_152, sii_153, sii_154, sik1_197, sik1_198, sik1_200, \
                         ski_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = pb_y[k] * sik0_197[k]
                   + f_15 * sii_152[k]
                   - f_12 * pc_y[k] * sik1_197[k];

        t_306[k] = pb_y[k] * sik0_198[k]
                   + f_14 * sii_153[k]
                   - f_12 * pc_y[k] * sik1_198[k];

        t_307[k] = f_13 * sii_154[k]
                   + f_3 * pc_y[k] * ski_238[k];

        t_308[k] = pb_y[k] * sik0_200[k]
                   - f_12 * pc_y[k] * sik1_200[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, t_313, pc_x, sii_245, sii_246, sii_247, \
                         sii_248, sii_249, ski_245, ski_246, ski_247, ski_248, \
                         ski_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_16 * sii_245[k]
                   + f_3 * pc_x[k] * ski_245[k];

        t_310[k] = f_16 * sii_246[k]
                   + f_3 * pc_x[k] * ski_246[k];

        t_311[k] = f_16 * sii_247[k]
                   + f_3 * pc_x[k] * ski_247[k];

        t_312[k] = f_16 * sii_248[k]
                   + f_3 * pc_x[k] * ski_248[k];

        t_313[k] = f_16 * sii_249[k]
                   + f_3 * pc_x[k] * ski_249[k];
    }

#pragma omp simd aligned(t_314, t_315, t_316, t_317, pc_x, pc_y, pc_z, sii_133, sii_161, \
                         sii_250, sii_251, skh0_183, skh1_183, ski_245, ski_250, \
                         ski_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_314[k] = f_16 * sii_250[k]
                   + f_3 * pc_x[k] * ski_250[k];

        t_315[k] = f_16 * sii_251[k]
                   + f_3 * pc_x[k] * ski_251[k];

        t_316[k] = f_13 * sii_161[k]
                   + f_1 * skh0_183[k]
                   - f_2 * skh1_183[k]
                   + f_3 * pc_y[k] * ski_245[k];

        t_317[k] = f_14 * sii_133[k]
                   + f_3 * pc_z[k] * ski_245[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, pc_y, sii_163, sii_164, sii_165, skh0_185, \
                         skh0_186, skh0_187, skh1_185, skh1_186, skh1_187, ski_247, ski_248, \
                         ski_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = f_13 * sii_163[k]
                   + f_4 * skh0_185[k]
                   - f_5 * skh1_185[k]
                   + f_3 * pc_y[k] * ski_247[k];

        t_319[k] = f_13 * sii_164[k]
                   + f_6 * skh0_186[k]
                   - f_7 * skh1_186[k]
                   + f_3 * pc_y[k] * ski_248[k];

        t_320[k] = f_13 * sii_165[k]
                   + f_8 * skh0_187[k]
                   - f_9 * skh1_187[k]
                   + f_3 * pc_y[k] * ski_249[k];
    }

#pragma omp simd aligned(t_321, t_322, t_323, pb_y, pc_y, sik0_215, sii_166, sii_167, \
                         sik1_215, skh0_188, skh1_188, ski_250, \
                         ski_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_321[k] = f_13 * sii_166[k]
                   + f_10 * skh0_188[k]
                   - f_11 * skh1_188[k]
                   + f_3 * pc_y[k] * ski_250[k];

        t_322[k] = f_13 * sii_167[k]
                   + f_3 * pc_y[k] * ski_251[k];

        t_323[k] = pb_y[k] * sik0_215[k]
                   - f_12 * pc_y[k] * sik1_215[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, pc_x, pc_y, pc_z, sii_140, sii_252, \
                         sii_255, skh0_189, skh0_192, skh1_189, skh1_192, ski_252, \
                         ski_255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = f_16 * sii_252[k]
                   + f_1 * skh0_189[k]
                   - f_2 * skh1_189[k]
                   + f_3 * pc_x[k] * ski_252[k];

        t_325[k] = f_3 * pc_y[k] * ski_252[k];

        t_326[k] = f_15 * sii_140[k]
                   + f_3 * pc_z[k] * ski_252[k];

        t_327[k] = f_16 * sii_255[k]
                   + f_4 * skh0_192[k]
                   - f_5 * skh1_192[k]
                   + f_3 * pc_x[k] * ski_255[k];
    }

#pragma omp simd aligned(t_328, t_329, t_330, pc_x, pc_y, sii_257, sii_258, skh0_194, \
                         skh0_195, skh1_194, skh1_195, ski_254, ski_257, \
                         ski_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_328[k] = f_3 * pc_y[k] * ski_254[k];

        t_329[k] = f_16 * sii_257[k]
                   + f_4 * skh0_194[k]
                   - f_5 * skh1_194[k]
                   + f_3 * pc_x[k] * ski_257[k];

        t_330[k] = f_16 * sii_258[k]
                   + f_6 * skh0_195[k]
                   - f_7 * skh1_195[k]
                   + f_3 * pc_x[k] * ski_258[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, pc_x, pc_y, pc_z, sii_143, sii_261, skh0_198, \
                         skh1_198, ski_255, ski_257, ski_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = f_15 * sii_143[k]
                   + f_3 * pc_z[k] * ski_255[k];

        t_332[k] = f_3 * pc_y[k] * ski_257[k];

        t_333[k] = f_16 * sii_261[k]
                   + f_6 * skh0_198[k]
                   - f_7 * skh1_198[k]
                   + f_3 * pc_x[k] * ski_261[k];
    }

#pragma omp simd aligned(t_334, t_335, t_336, pc_x, pc_z, sii_146, sii_262, sii_264, skh0_199, \
                         skh0_201, skh1_199, skh1_201, ski_258, ski_262, \
                         ski_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_334[k] = f_16 * sii_262[k]
                   + f_8 * skh0_199[k]
                   - f_9 * skh1_199[k]
                   + f_3 * pc_x[k] * ski_262[k];

        t_335[k] = f_15 * sii_146[k]
                   + f_3 * pc_z[k] * ski_258[k];

        t_336[k] = f_16 * sii_264[k]
                   + f_8 * skh0_201[k]
                   - f_9 * skh1_201[k]
                   + f_3 * pc_x[k] * ski_264[k];
    }

#pragma omp simd aligned(t_337, t_338, t_339, pc_x, pc_y, sii_266, sii_267, skh0_203, \
                         skh0_204, skh1_203, skh1_204, ski_261, ski_266, \
                         ski_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_337[k] = f_3 * pc_y[k] * ski_261[k];

        t_338[k] = f_16 * sii_266[k]
                   + f_8 * skh0_203[k]
                   - f_9 * skh1_203[k]
                   + f_3 * pc_x[k] * ski_266[k];

        t_339[k] = f_16 * sii_267[k]
                   + f_10 * skh0_204[k]
                   - f_11 * skh1_204[k]
                   + f_3 * pc_x[k] * ski_267[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, pc_x, pc_z, sii_150, sii_269, sii_270, skh0_206, \
                         skh0_207, skh1_206, skh1_207, ski_262, ski_269, \
                         ski_270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = f_15 * sii_150[k]
                   + f_3 * pc_z[k] * ski_262[k];

        t_341[k] = f_16 * sii_269[k]
                   + f_10 * skh0_206[k]
                   - f_11 * skh1_206[k]
                   + f_3 * pc_x[k] * ski_269[k];

        t_342[k] = f_16 * sii_270[k]
                   + f_10 * skh0_207[k]
                   - f_11 * skh1_207[k]
                   + f_3 * pc_x[k] * ski_270[k];
    }

#pragma omp simd aligned(t_343, t_344, t_345, t_346, pc_x, pc_y, sii_272, sii_273, sii_274, \
                         skh0_209, skh1_209, ski_266, ski_272, ski_273, \
                         ski_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_343[k] = f_3 * pc_y[k] * ski_266[k];

        t_344[k] = f_16 * sii_272[k]
                   + f_10 * skh0_209[k]
                   - f_11 * skh1_209[k]
                   + f_3 * pc_x[k] * ski_272[k];

        t_345[k] = f_16 * sii_273[k]
                   + f_3 * pc_x[k] * ski_273[k];

        t_346[k] = f_16 * sii_274[k]
                   + f_3 * pc_x[k] * ski_274[k];
    }
}

static auto
compute_prim_skk_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sik0,
                                                          const size_t sii, const size_t sik1,
                                                          const size_t skh0, const size_t skh1,
                                                          const size_t ski, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.0 / gamma;
    const auto f_2 = 3.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.0 / gamma;
    const auto f_5 = 2.0 * p / (gamma * q);
    const auto f_6 = 1.5 / gamma;
    const auto f_7 = 1.5 * p / (gamma * q);
    const auto f_8 = 1.0 / gamma;
    const auto f_9 = p / (gamma * q);
    const auto f_10 = 0.5 / gamma;
    const auto f_11 = 0.5 * p / (gamma * q);
    const auto f_12 = gamma / q;
    const auto f_13 = 0.5 / q;
    const auto f_14 = 1.0 / q;
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.0 / q;

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

    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sik0_216 = buffer.data(sik0 + 216);
    const auto *sik0_219 = buffer.data(sik0 + 219);
    const auto *sik0_222 = buffer.data(sik0 + 222);
    const auto *sik0_226 = buffer.data(sik0 + 226);
    const auto *sik0_228 = buffer.data(sik0 + 228);
    const auto *sik0_231 = buffer.data(sik0 + 231);
    const auto *sik0_233 = buffer.data(sik0 + 233);
    const auto *sik0_234 = buffer.data(sik0 + 234);
    const auto *sik0_244 = buffer.data(sik0 + 244);

    const auto *sii_161 = buffer.data(sii + 161);
    const auto *sii_167 = buffer.data(sii + 167);
    const auto *sii_168 = buffer.data(sii + 168);
    const auto *sii_170 = buffer.data(sii + 170);
    const auto *sii_171 = buffer.data(sii + 171);
    const auto *sii_173 = buffer.data(sii + 173);
    const auto *sii_174 = buffer.data(sii + 174);
    const auto *sii_175 = buffer.data(sii + 175);
    const auto *sii_177 = buffer.data(sii + 177);
    const auto *sii_178 = buffer.data(sii + 178);
    const auto *sii_179 = buffer.data(sii + 179);
    const auto *sii_180 = buffer.data(sii + 180);
    const auto *sii_182 = buffer.data(sii + 182);
    const auto *sii_189 = buffer.data(sii + 189);
    const auto *sii_191 = buffer.data(sii + 191);
    const auto *sii_192 = buffer.data(sii + 192);
    const auto *sii_193 = buffer.data(sii + 193);
    const auto *sii_194 = buffer.data(sii + 194);
    const auto *sii_195 = buffer.data(sii + 195);
    const auto *sii_196 = buffer.data(sii + 196);
    const auto *sii_198 = buffer.data(sii + 198);
    const auto *sii_199 = buffer.data(sii + 199);
    const auto *sii_201 = buffer.data(sii + 201);
    const auto *sii_202 = buffer.data(sii + 202);
    const auto *sii_205 = buffer.data(sii + 205);
    const auto *sii_206 = buffer.data(sii + 206);
    const auto *sii_210 = buffer.data(sii + 210);
    const auto *sii_219 = buffer.data(sii + 219);
    const auto *sii_220 = buffer.data(sii + 220);
    const auto *sii_221 = buffer.data(sii + 221);
    const auto *sii_222 = buffer.data(sii + 222);
    const auto *sii_223 = buffer.data(sii + 223);
    const auto *sii_224 = buffer.data(sii + 224);
    const auto *sii_226 = buffer.data(sii + 226);
    const auto *sii_229 = buffer.data(sii + 229);
    const auto *sii_233 = buffer.data(sii + 233);
    const auto *sii_238 = buffer.data(sii + 238);
    const auto *sii_275 = buffer.data(sii + 275);
    const auto *sii_276 = buffer.data(sii + 276);
    const auto *sii_277 = buffer.data(sii + 277);
    const auto *sii_278 = buffer.data(sii + 278);
    const auto *sii_279 = buffer.data(sii + 279);
    const auto *sii_280 = buffer.data(sii + 280);
    const auto *sii_283 = buffer.data(sii + 283);
    const auto *sii_285 = buffer.data(sii + 285);
    const auto *sii_286 = buffer.data(sii + 286);
    const auto *sii_289 = buffer.data(sii + 289);
    const auto *sii_290 = buffer.data(sii + 290);
    const auto *sii_292 = buffer.data(sii + 292);
    const auto *sii_294 = buffer.data(sii + 294);
    const auto *sii_295 = buffer.data(sii + 295);
    const auto *sii_297 = buffer.data(sii + 297);
    const auto *sii_298 = buffer.data(sii + 298);
    const auto *sii_300 = buffer.data(sii + 300);
    const auto *sii_301 = buffer.data(sii + 301);
    const auto *sii_302 = buffer.data(sii + 302);
    const auto *sii_303 = buffer.data(sii + 303);
    const auto *sii_304 = buffer.data(sii + 304);
    const auto *sii_305 = buffer.data(sii + 305);
    const auto *sii_306 = buffer.data(sii + 306);
    const auto *sii_307 = buffer.data(sii + 307);
    const auto *sii_313 = buffer.data(sii + 313);
    const auto *sii_317 = buffer.data(sii + 317);
    const auto *sii_322 = buffer.data(sii + 322);
    const auto *sii_328 = buffer.data(sii + 328);
    const auto *sii_329 = buffer.data(sii + 329);
    const auto *sii_330 = buffer.data(sii + 330);
    const auto *sii_331 = buffer.data(sii + 331);
    const auto *sii_332 = buffer.data(sii + 332);
    const auto *sii_333 = buffer.data(sii + 333);
    const auto *sii_334 = buffer.data(sii + 334);
    const auto *sii_335 = buffer.data(sii + 335);
    const auto *sii_336 = buffer.data(sii + 336);
    const auto *sii_339 = buffer.data(sii + 339);
    const auto *sii_341 = buffer.data(sii + 341);
    const auto *sii_342 = buffer.data(sii + 342);
    const auto *sii_345 = buffer.data(sii + 345);
    const auto *sii_346 = buffer.data(sii + 346);
    const auto *sii_348 = buffer.data(sii + 348);
    const auto *sii_350 = buffer.data(sii + 350);
    const auto *sii_351 = buffer.data(sii + 351);
    const auto *sii_353 = buffer.data(sii + 353);
    const auto *sii_354 = buffer.data(sii + 354);
    const auto *sii_356 = buffer.data(sii + 356);
    const auto *sii_357 = buffer.data(sii + 357);
    const auto *sii_358 = buffer.data(sii + 358);
    const auto *sii_359 = buffer.data(sii + 359);

    const auto *sik1_216 = buffer.data(sik1 + 216);
    const auto *sik1_219 = buffer.data(sik1 + 219);
    const auto *sik1_222 = buffer.data(sik1 + 222);
    const auto *sik1_226 = buffer.data(sik1 + 226);
    const auto *sik1_228 = buffer.data(sik1 + 228);
    const auto *sik1_231 = buffer.data(sik1 + 231);
    const auto *sik1_233 = buffer.data(sik1 + 233);
    const auto *sik1_234 = buffer.data(sik1 + 234);
    const auto *sik1_244 = buffer.data(sik1 + 244);

    const auto *skh0_204 = buffer.data(skh0 + 204);
    const auto *skh0_206 = buffer.data(skh0 + 206);
    const auto *skh0_207 = buffer.data(skh0 + 207);
    const auto *skh0_208 = buffer.data(skh0 + 208);
    const auto *skh0_209 = buffer.data(skh0 + 209);
    const auto *skh0_210 = buffer.data(skh0 + 210);
    const auto *skh0_213 = buffer.data(skh0 + 213);
    const auto *skh0_215 = buffer.data(skh0 + 215);
    const auto *skh0_216 = buffer.data(skh0 + 216);
    const auto *skh0_219 = buffer.data(skh0 + 219);
    const auto *skh0_220 = buffer.data(skh0 + 220);
    const auto *skh0_222 = buffer.data(skh0 + 222);
    const auto *skh0_224 = buffer.data(skh0 + 224);
    const auto *skh0_225 = buffer.data(skh0 + 225);
    const auto *skh0_227 = buffer.data(skh0 + 227);
    const auto *skh0_228 = buffer.data(skh0 + 228);
    const auto *skh0_229 = buffer.data(skh0 + 229);
    const auto *skh0_230 = buffer.data(skh0 + 230);
    const auto *skh0_236 = buffer.data(skh0 + 236);
    const auto *skh0_240 = buffer.data(skh0 + 240);
    const auto *skh0_245 = buffer.data(skh0 + 245);
    const auto *skh0_248 = buffer.data(skh0 + 248);
    const auto *skh0_249 = buffer.data(skh0 + 249);
    const auto *skh0_250 = buffer.data(skh0 + 250);
    const auto *skh0_251 = buffer.data(skh0 + 251);
    const auto *skh0_252 = buffer.data(skh0 + 252);
    const auto *skh0_255 = buffer.data(skh0 + 255);
    const auto *skh0_257 = buffer.data(skh0 + 257);
    const auto *skh0_258 = buffer.data(skh0 + 258);
    const auto *skh0_261 = buffer.data(skh0 + 261);
    const auto *skh0_262 = buffer.data(skh0 + 262);
    const auto *skh0_264 = buffer.data(skh0 + 264);
    const auto *skh0_266 = buffer.data(skh0 + 266);
    const auto *skh0_267 = buffer.data(skh0 + 267);
    const auto *skh0_269 = buffer.data(skh0 + 269);
    const auto *skh0_270 = buffer.data(skh0 + 270);
    const auto *skh0_272 = buffer.data(skh0 + 272);

    const auto *skh1_204 = buffer.data(skh1 + 204);
    const auto *skh1_206 = buffer.data(skh1 + 206);
    const auto *skh1_207 = buffer.data(skh1 + 207);
    const auto *skh1_208 = buffer.data(skh1 + 208);
    const auto *skh1_209 = buffer.data(skh1 + 209);
    const auto *skh1_210 = buffer.data(skh1 + 210);
    const auto *skh1_213 = buffer.data(skh1 + 213);
    const auto *skh1_215 = buffer.data(skh1 + 215);
    const auto *skh1_216 = buffer.data(skh1 + 216);
    const auto *skh1_219 = buffer.data(skh1 + 219);
    const auto *skh1_220 = buffer.data(skh1 + 220);
    const auto *skh1_222 = buffer.data(skh1 + 222);
    const auto *skh1_224 = buffer.data(skh1 + 224);
    const auto *skh1_225 = buffer.data(skh1 + 225);
    const auto *skh1_227 = buffer.data(skh1 + 227);
    const auto *skh1_228 = buffer.data(skh1 + 228);
    const auto *skh1_229 = buffer.data(skh1 + 229);
    const auto *skh1_230 = buffer.data(skh1 + 230);
    const auto *skh1_236 = buffer.data(skh1 + 236);
    const auto *skh1_240 = buffer.data(skh1 + 240);
    const auto *skh1_245 = buffer.data(skh1 + 245);
    const auto *skh1_248 = buffer.data(skh1 + 248);
    const auto *skh1_249 = buffer.data(skh1 + 249);
    const auto *skh1_250 = buffer.data(skh1 + 250);
    const auto *skh1_251 = buffer.data(skh1 + 251);
    const auto *skh1_252 = buffer.data(skh1 + 252);
    const auto *skh1_255 = buffer.data(skh1 + 255);
    const auto *skh1_257 = buffer.data(skh1 + 257);
    const auto *skh1_258 = buffer.data(skh1 + 258);
    const auto *skh1_261 = buffer.data(skh1 + 261);
    const auto *skh1_262 = buffer.data(skh1 + 262);
    const auto *skh1_264 = buffer.data(skh1 + 264);
    const auto *skh1_266 = buffer.data(skh1 + 266);
    const auto *skh1_267 = buffer.data(skh1 + 267);
    const auto *skh1_269 = buffer.data(skh1 + 269);
    const auto *skh1_270 = buffer.data(skh1 + 270);
    const auto *skh1_272 = buffer.data(skh1 + 272);

    const auto *ski_273 = buffer.data(ski + 273);
    const auto *ski_275 = buffer.data(ski + 275);
    const auto *ski_276 = buffer.data(ski + 276);
    const auto *ski_277 = buffer.data(ski + 277);
    const auto *ski_278 = buffer.data(ski + 278);
    const auto *ski_279 = buffer.data(ski + 279);
    const auto *ski_280 = buffer.data(ski + 280);
    const auto *ski_282 = buffer.data(ski + 282);
    const auto *ski_283 = buffer.data(ski + 283);
    const auto *ski_285 = buffer.data(ski + 285);
    const auto *ski_286 = buffer.data(ski + 286);
    const auto *ski_289 = buffer.data(ski + 289);
    const auto *ski_290 = buffer.data(ski + 290);
    const auto *ski_292 = buffer.data(ski + 292);
    const auto *ski_294 = buffer.data(ski + 294);
    const auto *ski_295 = buffer.data(ski + 295);
    const auto *ski_297 = buffer.data(ski + 297);
    const auto *ski_298 = buffer.data(ski + 298);
    const auto *ski_300 = buffer.data(ski + 300);
    const auto *ski_301 = buffer.data(ski + 301);
    const auto *ski_302 = buffer.data(ski + 302);
    const auto *ski_303 = buffer.data(ski + 303);
    const auto *ski_304 = buffer.data(ski + 304);
    const auto *ski_305 = buffer.data(ski + 305);
    const auto *ski_306 = buffer.data(ski + 306);
    const auto *ski_307 = buffer.data(ski + 307);
    const auto *ski_308 = buffer.data(ski + 308);
    const auto *ski_310 = buffer.data(ski + 310);
    const auto *ski_311 = buffer.data(ski + 311);
    const auto *ski_313 = buffer.data(ski + 313);
    const auto *ski_314 = buffer.data(ski + 314);
    const auto *ski_317 = buffer.data(ski + 317);
    const auto *ski_318 = buffer.data(ski + 318);
    const auto *ski_322 = buffer.data(ski + 322);
    const auto *ski_328 = buffer.data(ski + 328);
    const auto *ski_329 = buffer.data(ski + 329);
    const auto *ski_330 = buffer.data(ski + 330);
    const auto *ski_331 = buffer.data(ski + 331);
    const auto *ski_332 = buffer.data(ski + 332);
    const auto *ski_333 = buffer.data(ski + 333);
    const auto *ski_334 = buffer.data(ski + 334);
    const auto *ski_335 = buffer.data(ski + 335);
    const auto *ski_336 = buffer.data(ski + 336);
    const auto *ski_338 = buffer.data(ski + 338);
    const auto *ski_339 = buffer.data(ski + 339);
    const auto *ski_341 = buffer.data(ski + 341);
    const auto *ski_342 = buffer.data(ski + 342);
    const auto *ski_345 = buffer.data(ski + 345);
    const auto *ski_346 = buffer.data(ski + 346);
    const auto *ski_348 = buffer.data(ski + 348);
    const auto *ski_350 = buffer.data(ski + 350);
    const auto *ski_351 = buffer.data(ski + 351);
    const auto *ski_353 = buffer.data(ski + 353);
    const auto *ski_354 = buffer.data(ski + 354);
    const auto *ski_356 = buffer.data(ski + 356);
    const auto *ski_357 = buffer.data(ski + 357);
    const auto *ski_358 = buffer.data(ski + 358);
    const auto *ski_359 = buffer.data(ski + 359);

#pragma omp simd aligned(t_347, t_348, t_349, t_350, t_351, pc_x, sii_275, sii_276, sii_277, \
                         sii_278, sii_279, ski_275, ski_276, ski_277, ski_278, \
                         ski_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_347[k] = f_16 * sii_275[k]
                   + f_3 * pc_x[k] * ski_275[k];

        t_348[k] = f_16 * sii_276[k]
                   + f_3 * pc_x[k] * ski_276[k];

        t_349[k] = f_16 * sii_277[k]
                   + f_3 * pc_x[k] * ski_277[k];

        t_350[k] = f_16 * sii_278[k]
                   + f_3 * pc_x[k] * ski_278[k];

        t_351[k] = f_16 * sii_279[k]
                   + f_3 * pc_x[k] * ski_279[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, t_355, pc_y, pc_z, sii_161, skh0_204, skh0_206, \
                         skh0_207, skh1_204, skh1_206, skh1_207, ski_273, ski_275, \
                         ski_276 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = f_1 * skh0_204[k]
                   - f_2 * skh1_204[k]
                   + f_3 * pc_y[k] * ski_273[k];

        t_353[k] = f_15 * sii_161[k]
                   + f_3 * pc_z[k] * ski_273[k];

        t_354[k] = f_4 * skh0_206[k]
                   - f_5 * skh1_206[k]
                   + f_3 * pc_y[k] * ski_275[k];

        t_355[k] = f_6 * skh0_207[k]
                   - f_7 * skh1_207[k]
                   + f_3 * pc_y[k] * ski_276[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, t_359, pc_y, pc_z, sii_167, skh0_208, skh0_209, \
                         skh1_208, skh1_209, ski_277, ski_278, \
                         ski_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_8 * skh0_208[k]
                   - f_9 * skh1_208[k]
                   + f_3 * pc_y[k] * ski_277[k];

        t_357[k] = f_10 * skh0_209[k]
                   - f_11 * skh1_209[k]
                   + f_3 * pc_y[k] * ski_278[k];

        t_358[k] = f_3 * pc_y[k] * ski_279[k];

        t_359[k] = f_15 * sii_167[k]
                   + f_1 * skh0_209[k]
                   - f_2 * skh1_209[k]
                   + f_3 * pc_z[k] * ski_279[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, pc_x, pc_y, pc_z, sii_168, sii_280, \
                         sii_283, skh0_210, skh0_213, skh1_210, skh1_213, ski_280, \
                         ski_283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = f_15 * sii_280[k]
                   + f_1 * skh0_210[k]
                   - f_2 * skh1_210[k]
                   + f_3 * pc_x[k] * ski_280[k];

        t_361[k] = f_16 * sii_168[k]
                   + f_3 * pc_y[k] * ski_280[k];

        t_362[k] = f_3 * pc_z[k] * ski_280[k];

        t_363[k] = f_15 * sii_283[k]
                   + f_4 * skh0_213[k]
                   - f_5 * skh1_213[k]
                   + f_3 * pc_x[k] * ski_283[k];
    }

#pragma omp simd aligned(t_364, t_365, t_366, pc_x, pc_y, sii_170, sii_285, sii_286, skh0_215, \
                         skh0_216, skh1_215, skh1_216, ski_282, ski_285, \
                         ski_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_364[k] = f_16 * sii_170[k]
                   + f_3 * pc_y[k] * ski_282[k];

        t_365[k] = f_15 * sii_285[k]
                   + f_4 * skh0_215[k]
                   - f_5 * skh1_215[k]
                   + f_3 * pc_x[k] * ski_285[k];

        t_366[k] = f_15 * sii_286[k]
                   + f_6 * skh0_216[k]
                   - f_7 * skh1_216[k]
                   + f_3 * pc_x[k] * ski_286[k];
    }

#pragma omp simd aligned(t_367, t_368, t_369, pc_x, pc_y, pc_z, sii_173, sii_289, skh0_219, \
                         skh1_219, ski_283, ski_285, ski_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_367[k] = f_3 * pc_z[k] * ski_283[k];

        t_368[k] = f_16 * sii_173[k]
                   + f_3 * pc_y[k] * ski_285[k];

        t_369[k] = f_15 * sii_289[k]
                   + f_6 * skh0_219[k]
                   - f_7 * skh1_219[k]
                   + f_3 * pc_x[k] * ski_289[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, pc_x, pc_z, sii_290, sii_292, skh0_220, \
                         skh0_222, skh1_220, skh1_222, ski_286, ski_290, \
                         ski_292 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = f_15 * sii_290[k]
                   + f_8 * skh0_220[k]
                   - f_9 * skh1_220[k]
                   + f_3 * pc_x[k] * ski_290[k];

        t_371[k] = f_3 * pc_z[k] * ski_286[k];

        t_372[k] = f_15 * sii_292[k]
                   + f_8 * skh0_222[k]
                   - f_9 * skh1_222[k]
                   + f_3 * pc_x[k] * ski_292[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, pc_x, pc_y, sii_177, sii_294, sii_295, skh0_224, \
                         skh0_225, skh1_224, skh1_225, ski_289, ski_294, \
                         ski_295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = f_16 * sii_177[k]
                   + f_3 * pc_y[k] * ski_289[k];

        t_374[k] = f_15 * sii_294[k]
                   + f_8 * skh0_224[k]
                   - f_9 * skh1_224[k]
                   + f_3 * pc_x[k] * ski_294[k];

        t_375[k] = f_15 * sii_295[k]
                   + f_10 * skh0_225[k]
                   - f_11 * skh1_225[k]
                   + f_3 * pc_x[k] * ski_295[k];
    }

#pragma omp simd aligned(t_376, t_377, t_378, pc_x, pc_z, sii_297, sii_298, skh0_227, \
                         skh0_228, skh1_227, skh1_228, ski_290, ski_297, \
                         ski_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_376[k] = f_3 * pc_z[k] * ski_290[k];

        t_377[k] = f_15 * sii_297[k]
                   + f_10 * skh0_227[k]
                   - f_11 * skh1_227[k]
                   + f_3 * pc_x[k] * ski_297[k];

        t_378[k] = f_15 * sii_298[k]
                   + f_10 * skh0_228[k]
                   - f_11 * skh1_228[k]
                   + f_3 * pc_x[k] * ski_298[k];
    }

#pragma omp simd aligned(t_379, t_380, t_381, t_382, pc_x, pc_y, sii_182, sii_300, sii_301, \
                         sii_302, skh0_230, skh1_230, ski_294, ski_300, ski_301, \
                         ski_302 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_379[k] = f_16 * sii_182[k]
                   + f_3 * pc_y[k] * ski_294[k];

        t_380[k] = f_15 * sii_300[k]
                   + f_10 * skh0_230[k]
                   - f_11 * skh1_230[k]
                   + f_3 * pc_x[k] * ski_300[k];

        t_381[k] = f_15 * sii_301[k]
                   + f_3 * pc_x[k] * ski_301[k];

        t_382[k] = f_15 * sii_302[k]
                   + f_3 * pc_x[k] * ski_302[k];
    }

#pragma omp simd aligned(t_383, t_384, t_385, t_386, t_387, pc_x, sii_303, sii_304, sii_305, \
                         sii_306, sii_307, ski_303, ski_304, ski_305, ski_306, \
                         ski_307 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_383[k] = f_15 * sii_303[k]
                   + f_3 * pc_x[k] * ski_303[k];

        t_384[k] = f_15 * sii_304[k]
                   + f_3 * pc_x[k] * ski_304[k];

        t_385[k] = f_15 * sii_305[k]
                   + f_3 * pc_x[k] * ski_305[k];

        t_386[k] = f_15 * sii_306[k]
                   + f_3 * pc_x[k] * ski_306[k];

        t_387[k] = f_15 * sii_307[k]
                   + f_3 * pc_x[k] * ski_307[k];
    }

#pragma omp simd aligned(t_388, t_389, t_390, pc_y, pc_z, sii_189, sii_191, skh0_225, \
                         skh0_227, skh1_225, skh1_227, ski_301, \
                         ski_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_388[k] = f_16 * sii_189[k]
                   + f_1 * skh0_225[k]
                   - f_2 * skh1_225[k]
                   + f_3 * pc_y[k] * ski_301[k];

        t_389[k] = f_3 * pc_z[k] * ski_301[k];

        t_390[k] = f_16 * sii_191[k]
                   + f_4 * skh0_227[k]
                   - f_5 * skh1_227[k]
                   + f_3 * pc_y[k] * ski_303[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, pc_y, sii_192, sii_193, sii_194, skh0_228, \
                         skh0_229, skh0_230, skh1_228, skh1_229, skh1_230, ski_304, ski_305, \
                         ski_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = f_16 * sii_192[k]
                   + f_6 * skh0_228[k]
                   - f_7 * skh1_228[k]
                   + f_3 * pc_y[k] * ski_304[k];

        t_392[k] = f_16 * sii_193[k]
                   + f_8 * skh0_229[k]
                   - f_9 * skh1_229[k]
                   + f_3 * pc_y[k] * ski_305[k];

        t_393[k] = f_16 * sii_194[k]
                   + f_10 * skh0_230[k]
                   - f_11 * skh1_230[k]
                   + f_3 * pc_y[k] * ski_306[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, t_397, pb_z, pc_y, pc_z, sik0_216, sii_195, \
                         sii_196, sik1_216, skh0_230, skh1_230, ski_307, \
                         ski_308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = f_16 * sii_195[k]
                   + f_3 * pc_y[k] * ski_307[k];

        t_395[k] = f_1 * skh0_230[k]
                   - f_2 * skh1_230[k]
                   + f_3 * pc_z[k] * ski_307[k];

        t_396[k] = pb_z[k] * sik0_216[k]
                   - f_12 * pc_z[k] * sik1_216[k];

        t_397[k] = f_15 * sii_196[k]
                   + f_3 * pc_y[k] * ski_308[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, pb_z, pc_y, pc_z, sik0_219, sii_168, sii_198, \
                         sik1_219, ski_308, ski_310 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = f_13 * sii_168[k]
                   + f_3 * pc_z[k] * ski_308[k];

        t_399[k] = pb_z[k] * sik0_219[k]
                   - f_12 * pc_z[k] * sik1_219[k];

        t_400[k] = f_15 * sii_198[k]
                   + f_3 * pc_y[k] * ski_310[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, pb_z, pc_x, pc_z, sik0_222, sii_171, sii_313, \
                         sik1_222, skh0_236, skh1_236, ski_311, \
                         ski_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_15 * sii_313[k]
                   + f_4 * skh0_236[k]
                   - f_5 * skh1_236[k]
                   + f_3 * pc_x[k] * ski_313[k];

        t_402[k] = pb_z[k] * sik0_222[k]
                   - f_12 * pc_z[k] * sik1_222[k];

        t_403[k] = f_13 * sii_171[k]
                   + f_3 * pc_z[k] * ski_311[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, pb_z, pc_x, pc_y, pc_z, sik0_226, sii_201, \
                         sii_317, sik1_226, skh0_240, skh1_240, ski_313, \
                         ski_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = f_15 * sii_201[k]
                   + f_3 * pc_y[k] * ski_313[k];

        t_405[k] = f_15 * sii_317[k]
                   + f_6 * skh0_240[k]
                   - f_7 * skh1_240[k]
                   + f_3 * pc_x[k] * ski_317[k];

        t_406[k] = pb_z[k] * sik0_226[k]
                   - f_12 * pc_z[k] * sik1_226[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, pb_z, pc_y, pc_z, sik0_228, sii_174, sii_175, \
                         sii_205, sik1_228, ski_314, ski_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = f_13 * sii_174[k]
                   + f_3 * pc_z[k] * ski_314[k];

        t_408[k] = pb_z[k] * sik0_228[k]
                   + f_14 * sii_175[k]
                   - f_12 * pc_z[k] * sik1_228[k];

        t_409[k] = f_15 * sii_205[k]
                   + f_3 * pc_y[k] * ski_317[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, pb_z, pc_x, pc_z, sik0_231, sii_178, sii_322, \
                         sik1_231, skh0_245, skh1_245, ski_318, \
                         ski_322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = f_15 * sii_322[k]
                   + f_8 * skh0_245[k]
                   - f_9 * skh1_245[k]
                   + f_3 * pc_x[k] * ski_322[k];

        t_411[k] = pb_z[k] * sik0_231[k]
                   - f_12 * pc_z[k] * sik1_231[k];

        t_412[k] = f_13 * sii_178[k]
                   + f_3 * pc_z[k] * ski_318[k];
    }

#pragma omp simd aligned(t_413, t_414, t_415, pb_z, pc_y, pc_z, sik0_233, sik0_234, sii_179, \
                         sii_180, sii_210, sik1_233, sik1_234, \
                         ski_322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_413[k] = pb_z[k] * sik0_233[k]
                   + f_14 * sii_179[k]
                   - f_12 * pc_z[k] * sik1_233[k];

        t_414[k] = pb_z[k] * sik0_234[k]
                   + f_15 * sii_180[k]
                   - f_12 * pc_z[k] * sik1_234[k];

        t_415[k] = f_15 * sii_210[k]
                   + f_3 * pc_y[k] * ski_322[k];
    }

#pragma omp simd aligned(t_416, t_417, t_418, t_419, pc_x, sii_328, sii_329, sii_330, sii_331, \
                         skh0_251, skh1_251, ski_328, ski_329, ski_330, \
                         ski_331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_416[k] = f_15 * sii_328[k]
                   + f_10 * skh0_251[k]
                   - f_11 * skh1_251[k]
                   + f_3 * pc_x[k] * ski_328[k];

        t_417[k] = f_15 * sii_329[k]
                   + f_3 * pc_x[k] * ski_329[k];

        t_418[k] = f_15 * sii_330[k]
                   + f_3 * pc_x[k] * ski_330[k];

        t_419[k] = f_15 * sii_331[k]
                   + f_3 * pc_x[k] * ski_331[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, pc_x, sii_332, sii_333, sii_334, sii_335, \
                         ski_332, ski_333, ski_334, ski_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_15 * sii_332[k]
                   + f_3 * pc_x[k] * ski_332[k];

        t_421[k] = f_15 * sii_333[k]
                   + f_3 * pc_x[k] * ski_333[k];

        t_422[k] = f_15 * sii_334[k]
                   + f_3 * pc_x[k] * ski_334[k];

        t_423[k] = f_15 * sii_335[k]
                   + f_3 * pc_x[k] * ski_335[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, pb_z, pc_y, pc_z, sik0_244, sii_189, sii_219, \
                         sik1_244, skh0_248, skh1_248, ski_329, \
                         ski_331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = pb_z[k] * sik0_244[k]
                   - f_12 * pc_z[k] * sik1_244[k];

        t_425[k] = f_13 * sii_189[k]
                   + f_3 * pc_z[k] * ski_329[k];

        t_426[k] = f_15 * sii_219[k]
                   + f_4 * skh0_248[k]
                   - f_5 * skh1_248[k]
                   + f_3 * pc_y[k] * ski_331[k];
    }

#pragma omp simd aligned(t_427, t_428, t_429, pc_y, sii_220, sii_221, sii_222, skh0_249, \
                         skh0_250, skh0_251, skh1_249, skh1_250, skh1_251, ski_332, ski_333, \
                         ski_334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_427[k] = f_15 * sii_220[k]
                   + f_6 * skh0_249[k]
                   - f_7 * skh1_249[k]
                   + f_3 * pc_y[k] * ski_332[k];

        t_428[k] = f_15 * sii_221[k]
                   + f_8 * skh0_250[k]
                   - f_9 * skh1_250[k]
                   + f_3 * pc_y[k] * ski_333[k];

        t_429[k] = f_15 * sii_222[k]
                   + f_10 * skh0_251[k]
                   - f_11 * skh1_251[k]
                   + f_3 * pc_y[k] * ski_334[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, pc_x, pc_y, pc_z, sii_195, sii_223, sii_336, \
                         skh0_251, skh0_252, skh1_251, skh1_252, ski_335, \
                         ski_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = f_15 * sii_223[k]
                   + f_3 * pc_y[k] * ski_335[k];

        t_431[k] = f_13 * sii_195[k]
                   + f_1 * skh0_251[k]
                   - f_2 * skh1_251[k]
                   + f_3 * pc_z[k] * ski_335[k];

        t_432[k] = f_15 * sii_336[k]
                   + f_1 * skh0_252[k]
                   - f_2 * skh1_252[k]
                   + f_3 * pc_x[k] * ski_336[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, t_436, pc_x, pc_y, pc_z, sii_196, sii_224, \
                         sii_226, sii_339, skh0_255, skh1_255, ski_336, ski_338, \
                         ski_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = f_14 * sii_224[k]
                   + f_3 * pc_y[k] * ski_336[k];

        t_434[k] = f_14 * sii_196[k]
                   + f_3 * pc_z[k] * ski_336[k];

        t_435[k] = f_15 * sii_339[k]
                   + f_4 * skh0_255[k]
                   - f_5 * skh1_255[k]
                   + f_3 * pc_x[k] * ski_339[k];

        t_436[k] = f_14 * sii_226[k]
                   + f_3 * pc_y[k] * ski_338[k];
    }

#pragma omp simd aligned(t_437, t_438, t_439, pc_x, pc_z, sii_199, sii_341, sii_342, skh0_257, \
                         skh0_258, skh1_257, skh1_258, ski_339, ski_341, \
                         ski_342 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_437[k] = f_15 * sii_341[k]
                   + f_4 * skh0_257[k]
                   - f_5 * skh1_257[k]
                   + f_3 * pc_x[k] * ski_341[k];

        t_438[k] = f_15 * sii_342[k]
                   + f_6 * skh0_258[k]
                   - f_7 * skh1_258[k]
                   + f_3 * pc_x[k] * ski_342[k];

        t_439[k] = f_14 * sii_199[k]
                   + f_3 * pc_z[k] * ski_339[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, pc_x, pc_y, sii_229, sii_345, sii_346, skh0_261, \
                         skh0_262, skh1_261, skh1_262, ski_341, ski_345, \
                         ski_346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = f_14 * sii_229[k]
                   + f_3 * pc_y[k] * ski_341[k];

        t_441[k] = f_15 * sii_345[k]
                   + f_6 * skh0_261[k]
                   - f_7 * skh1_261[k]
                   + f_3 * pc_x[k] * ski_345[k];

        t_442[k] = f_15 * sii_346[k]
                   + f_8 * skh0_262[k]
                   - f_9 * skh1_262[k]
                   + f_3 * pc_x[k] * ski_346[k];
    }

#pragma omp simd aligned(t_443, t_444, t_445, pc_x, pc_y, pc_z, sii_202, sii_233, sii_348, \
                         skh0_264, skh1_264, ski_342, ski_345, \
                         ski_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_443[k] = f_14 * sii_202[k]
                   + f_3 * pc_z[k] * ski_342[k];

        t_444[k] = f_15 * sii_348[k]
                   + f_8 * skh0_264[k]
                   - f_9 * skh1_264[k]
                   + f_3 * pc_x[k] * ski_348[k];

        t_445[k] = f_14 * sii_233[k]
                   + f_3 * pc_y[k] * ski_345[k];
    }

#pragma omp simd aligned(t_446, t_447, t_448, pc_x, pc_z, sii_206, sii_350, sii_351, skh0_266, \
                         skh0_267, skh1_266, skh1_267, ski_346, ski_350, \
                         ski_351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_446[k] = f_15 * sii_350[k]
                   + f_8 * skh0_266[k]
                   - f_9 * skh1_266[k]
                   + f_3 * pc_x[k] * ski_350[k];

        t_447[k] = f_15 * sii_351[k]
                   + f_10 * skh0_267[k]
                   - f_11 * skh1_267[k]
                   + f_3 * pc_x[k] * ski_351[k];

        t_448[k] = f_14 * sii_206[k]
                   + f_3 * pc_z[k] * ski_346[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, pc_x, pc_y, sii_238, sii_353, sii_354, skh0_269, \
                         skh0_270, skh1_269, skh1_270, ski_350, ski_353, \
                         ski_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_15 * sii_353[k]
                   + f_10 * skh0_269[k]
                   - f_11 * skh1_269[k]
                   + f_3 * pc_x[k] * ski_353[k];

        t_450[k] = f_15 * sii_354[k]
                   + f_10 * skh0_270[k]
                   - f_11 * skh1_270[k]
                   + f_3 * pc_x[k] * ski_354[k];

        t_451[k] = f_14 * sii_238[k]
                   + f_3 * pc_y[k] * ski_350[k];
    }

#pragma omp simd aligned(t_452, t_453, t_454, t_455, pc_x, sii_356, sii_357, sii_358, sii_359, \
                         skh0_272, skh1_272, ski_356, ski_357, ski_358, \
                         ski_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_452[k] = f_15 * sii_356[k]
                   + f_10 * skh0_272[k]
                   - f_11 * skh1_272[k]
                   + f_3 * pc_x[k] * ski_356[k];

        t_453[k] = f_15 * sii_357[k]
                   + f_3 * pc_x[k] * ski_357[k];

        t_454[k] = f_15 * sii_358[k]
                   + f_3 * pc_x[k] * ski_358[k];

        t_455[k] = f_15 * sii_359[k]
                   + f_3 * pc_x[k] * ski_359[k];
    }
}

static auto
compute_prim_skk_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sik0,
                                                          const size_t sii, const size_t sik1,
                                                          const size_t skh0, const size_t skh1,
                                                          const size_t ski, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.0 / gamma;
    const auto f_2 = 3.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.0 / gamma;
    const auto f_5 = 2.0 * p / (gamma * q);
    const auto f_6 = 1.5 / gamma;
    const auto f_7 = 1.5 * p / (gamma * q);
    const auto f_8 = 1.0 / gamma;
    const auto f_9 = p / (gamma * q);
    const auto f_10 = 0.5 / gamma;
    const auto f_11 = 0.5 * p / (gamma * q);
    const auto f_12 = gamma / q;
    const auto f_13 = 0.5 / q;
    const auto f_14 = 1.0 / q;
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;

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
    auto *t_563 = buffer.data(target + 563);
    auto *t_564 = buffer.data(target + 564);
    auto *t_565 = buffer.data(target + 565);
    auto *t_566 = buffer.data(target + 566);
    auto *t_567 = buffer.data(target + 567);
    auto *t_568 = buffer.data(target + 568);
    auto *t_569 = buffer.data(target + 569);
    auto *t_570 = buffer.data(target + 570);

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sik0_324 = buffer.data(sik0 + 324);
    const auto *sik0_327 = buffer.data(sik0 + 327);
    const auto *sik0_329 = buffer.data(sik0 + 329);
    const auto *sik0_330 = buffer.data(sik0 + 330);
    const auto *sik0_333 = buffer.data(sik0 + 333);
    const auto *sik0_334 = buffer.data(sik0 + 334);
    const auto *sik0_336 = buffer.data(sik0 + 336);
    const auto *sik0_338 = buffer.data(sik0 + 338);
    const auto *sik0_339 = buffer.data(sik0 + 339);
    const auto *sik0_341 = buffer.data(sik0 + 341);
    const auto *sik0_342 = buffer.data(sik0 + 342);
    const auto *sik0_344 = buffer.data(sik0 + 344);
    const auto *sik0_359 = buffer.data(sik0 + 359);

    const auto *sii_217 = buffer.data(sii + 217);
    const auto *sii_223 = buffer.data(sii + 223);
    const auto *sii_224 = buffer.data(sii + 224);
    const auto *sii_227 = buffer.data(sii + 227);
    const auto *sii_230 = buffer.data(sii + 230);
    const auto *sii_234 = buffer.data(sii + 234);
    const auto *sii_245 = buffer.data(sii + 245);
    const auto *sii_247 = buffer.data(sii + 247);
    const auto *sii_248 = buffer.data(sii + 248);
    const auto *sii_249 = buffer.data(sii + 249);
    const auto *sii_250 = buffer.data(sii + 250);
    const auto *sii_251 = buffer.data(sii + 251);
    const auto *sii_252 = buffer.data(sii + 252);
    const auto *sii_253 = buffer.data(sii + 253);
    const auto *sii_254 = buffer.data(sii + 254);
    const auto *sii_255 = buffer.data(sii + 255);
    const auto *sii_257 = buffer.data(sii + 257);
    const auto *sii_258 = buffer.data(sii + 258);
    const auto *sii_260 = buffer.data(sii + 260);
    const auto *sii_261 = buffer.data(sii + 261);
    const auto *sii_262 = buffer.data(sii + 262);
    const auto *sii_264 = buffer.data(sii + 264);
    const auto *sii_265 = buffer.data(sii + 265);
    const auto *sii_266 = buffer.data(sii + 266);
    const auto *sii_273 = buffer.data(sii + 273);
    const auto *sii_275 = buffer.data(sii + 275);
    const auto *sii_276 = buffer.data(sii + 276);
    const auto *sii_277 = buffer.data(sii + 277);
    const auto *sii_278 = buffer.data(sii + 278);
    const auto *sii_279 = buffer.data(sii + 279);
    const auto *sii_280 = buffer.data(sii + 280);
    const auto *sii_282 = buffer.data(sii + 282);
    const auto *sii_285 = buffer.data(sii + 285);
    const auto *sii_289 = buffer.data(sii + 289);
    const auto *sii_294 = buffer.data(sii + 294);
    const auto *sii_301 = buffer.data(sii + 301);
    const auto *sii_303 = buffer.data(sii + 303);
    const auto *sii_360 = buffer.data(sii + 360);
    const auto *sii_361 = buffer.data(sii + 361);
    const auto *sii_362 = buffer.data(sii + 362);
    const auto *sii_363 = buffer.data(sii + 363);
    const auto *sii_385 = buffer.data(sii + 385);
    const auto *sii_386 = buffer.data(sii + 386);
    const auto *sii_387 = buffer.data(sii + 387);
    const auto *sii_388 = buffer.data(sii + 388);
    const auto *sii_389 = buffer.data(sii + 389);
    const auto *sii_390 = buffer.data(sii + 390);
    const auto *sii_391 = buffer.data(sii + 391);
    const auto *sii_392 = buffer.data(sii + 392);
    const auto *sii_395 = buffer.data(sii + 395);
    const auto *sii_397 = buffer.data(sii + 397);
    const auto *sii_398 = buffer.data(sii + 398);
    const auto *sii_401 = buffer.data(sii + 401);
    const auto *sii_402 = buffer.data(sii + 402);
    const auto *sii_404 = buffer.data(sii + 404);
    const auto *sii_406 = buffer.data(sii + 406);
    const auto *sii_407 = buffer.data(sii + 407);
    const auto *sii_409 = buffer.data(sii + 409);
    const auto *sii_410 = buffer.data(sii + 410);
    const auto *sii_412 = buffer.data(sii + 412);
    const auto *sii_413 = buffer.data(sii + 413);
    const auto *sii_414 = buffer.data(sii + 414);
    const auto *sii_415 = buffer.data(sii + 415);
    const auto *sii_416 = buffer.data(sii + 416);
    const auto *sii_417 = buffer.data(sii + 417);
    const auto *sii_418 = buffer.data(sii + 418);
    const auto *sii_419 = buffer.data(sii + 419);
    const auto *sii_420 = buffer.data(sii + 420);
    const auto *sii_423 = buffer.data(sii + 423);
    const auto *sii_425 = buffer.data(sii + 425);
    const auto *sii_426 = buffer.data(sii + 426);
    const auto *sii_429 = buffer.data(sii + 429);
    const auto *sii_430 = buffer.data(sii + 430);
    const auto *sii_432 = buffer.data(sii + 432);
    const auto *sii_434 = buffer.data(sii + 434);
    const auto *sii_435 = buffer.data(sii + 435);
    const auto *sii_437 = buffer.data(sii + 437);
    const auto *sii_438 = buffer.data(sii + 438);
    const auto *sii_440 = buffer.data(sii + 440);
    const auto *sii_441 = buffer.data(sii + 441);
    const auto *sii_442 = buffer.data(sii + 442);
    const auto *sii_443 = buffer.data(sii + 443);
    const auto *sii_444 = buffer.data(sii + 444);
    const auto *sii_445 = buffer.data(sii + 445);
    const auto *sii_446 = buffer.data(sii + 446);
    const auto *sii_447 = buffer.data(sii + 447);

    const auto *sik1_324 = buffer.data(sik1 + 324);
    const auto *sik1_327 = buffer.data(sik1 + 327);
    const auto *sik1_329 = buffer.data(sik1 + 329);
    const auto *sik1_330 = buffer.data(sik1 + 330);
    const auto *sik1_333 = buffer.data(sik1 + 333);
    const auto *sik1_334 = buffer.data(sik1 + 334);
    const auto *sik1_336 = buffer.data(sik1 + 336);
    const auto *sik1_338 = buffer.data(sik1 + 338);
    const auto *sik1_339 = buffer.data(sik1 + 339);
    const auto *sik1_341 = buffer.data(sik1 + 341);
    const auto *sik1_342 = buffer.data(sik1 + 342);
    const auto *sik1_344 = buffer.data(sik1 + 344);
    const auto *sik1_359 = buffer.data(sik1 + 359);

    const auto *skh0_267 = buffer.data(skh0 + 267);
    const auto *skh0_269 = buffer.data(skh0 + 269);
    const auto *skh0_270 = buffer.data(skh0 + 270);
    const auto *skh0_271 = buffer.data(skh0 + 271);
    const auto *skh0_272 = buffer.data(skh0 + 272);
    const auto *skh0_288 = buffer.data(skh0 + 288);
    const auto *skh0_290 = buffer.data(skh0 + 290);
    const auto *skh0_291 = buffer.data(skh0 + 291);
    const auto *skh0_292 = buffer.data(skh0 + 292);
    const auto *skh0_293 = buffer.data(skh0 + 293);
    const auto *skh0_294 = buffer.data(skh0 + 294);
    const auto *skh0_297 = buffer.data(skh0 + 297);
    const auto *skh0_299 = buffer.data(skh0 + 299);
    const auto *skh0_300 = buffer.data(skh0 + 300);
    const auto *skh0_303 = buffer.data(skh0 + 303);
    const auto *skh0_304 = buffer.data(skh0 + 304);
    const auto *skh0_306 = buffer.data(skh0 + 306);
    const auto *skh0_308 = buffer.data(skh0 + 308);
    const auto *skh0_309 = buffer.data(skh0 + 309);
    const auto *skh0_311 = buffer.data(skh0 + 311);
    const auto *skh0_312 = buffer.data(skh0 + 312);
    const auto *skh0_313 = buffer.data(skh0 + 313);
    const auto *skh0_314 = buffer.data(skh0 + 314);
    const auto *skh0_315 = buffer.data(skh0 + 315);
    const auto *skh0_318 = buffer.data(skh0 + 318);
    const auto *skh0_320 = buffer.data(skh0 + 320);
    const auto *skh0_321 = buffer.data(skh0 + 321);
    const auto *skh0_324 = buffer.data(skh0 + 324);
    const auto *skh0_325 = buffer.data(skh0 + 325);
    const auto *skh0_327 = buffer.data(skh0 + 327);
    const auto *skh0_329 = buffer.data(skh0 + 329);
    const auto *skh0_330 = buffer.data(skh0 + 330);
    const auto *skh0_332 = buffer.data(skh0 + 332);
    const auto *skh0_333 = buffer.data(skh0 + 333);
    const auto *skh0_335 = buffer.data(skh0 + 335);

    const auto *skh1_267 = buffer.data(skh1 + 267);
    const auto *skh1_269 = buffer.data(skh1 + 269);
    const auto *skh1_270 = buffer.data(skh1 + 270);
    const auto *skh1_271 = buffer.data(skh1 + 271);
    const auto *skh1_272 = buffer.data(skh1 + 272);
    const auto *skh1_288 = buffer.data(skh1 + 288);
    const auto *skh1_290 = buffer.data(skh1 + 290);
    const auto *skh1_291 = buffer.data(skh1 + 291);
    const auto *skh1_292 = buffer.data(skh1 + 292);
    const auto *skh1_293 = buffer.data(skh1 + 293);
    const auto *skh1_294 = buffer.data(skh1 + 294);
    const auto *skh1_297 = buffer.data(skh1 + 297);
    const auto *skh1_299 = buffer.data(skh1 + 299);
    const auto *skh1_300 = buffer.data(skh1 + 300);
    const auto *skh1_303 = buffer.data(skh1 + 303);
    const auto *skh1_304 = buffer.data(skh1 + 304);
    const auto *skh1_306 = buffer.data(skh1 + 306);
    const auto *skh1_308 = buffer.data(skh1 + 308);
    const auto *skh1_309 = buffer.data(skh1 + 309);
    const auto *skh1_311 = buffer.data(skh1 + 311);
    const auto *skh1_312 = buffer.data(skh1 + 312);
    const auto *skh1_313 = buffer.data(skh1 + 313);
    const auto *skh1_314 = buffer.data(skh1 + 314);
    const auto *skh1_315 = buffer.data(skh1 + 315);
    const auto *skh1_318 = buffer.data(skh1 + 318);
    const auto *skh1_320 = buffer.data(skh1 + 320);
    const auto *skh1_321 = buffer.data(skh1 + 321);
    const auto *skh1_324 = buffer.data(skh1 + 324);
    const auto *skh1_325 = buffer.data(skh1 + 325);
    const auto *skh1_327 = buffer.data(skh1 + 327);
    const auto *skh1_329 = buffer.data(skh1 + 329);
    const auto *skh1_330 = buffer.data(skh1 + 330);
    const auto *skh1_332 = buffer.data(skh1 + 332);
    const auto *skh1_333 = buffer.data(skh1 + 333);
    const auto *skh1_335 = buffer.data(skh1 + 335);

    const auto *ski_357 = buffer.data(ski + 357);
    const auto *ski_359 = buffer.data(ski + 359);
    const auto *ski_360 = buffer.data(ski + 360);
    const auto *ski_361 = buffer.data(ski + 361);
    const auto *ski_362 = buffer.data(ski + 362);
    const auto *ski_363 = buffer.data(ski + 363);
    const auto *ski_364 = buffer.data(ski + 364);
    const auto *ski_366 = buffer.data(ski + 366);
    const auto *ski_367 = buffer.data(ski + 367);
    const auto *ski_369 = buffer.data(ski + 369);
    const auto *ski_370 = buffer.data(ski + 370);
    const auto *ski_373 = buffer.data(ski + 373);
    const auto *ski_374 = buffer.data(ski + 374);
    const auto *ski_378 = buffer.data(ski + 378);
    const auto *ski_385 = buffer.data(ski + 385);
    const auto *ski_386 = buffer.data(ski + 386);
    const auto *ski_387 = buffer.data(ski + 387);
    const auto *ski_388 = buffer.data(ski + 388);
    const auto *ski_389 = buffer.data(ski + 389);
    const auto *ski_390 = buffer.data(ski + 390);
    const auto *ski_391 = buffer.data(ski + 391);
    const auto *ski_392 = buffer.data(ski + 392);
    const auto *ski_394 = buffer.data(ski + 394);
    const auto *ski_395 = buffer.data(ski + 395);
    const auto *ski_397 = buffer.data(ski + 397);
    const auto *ski_398 = buffer.data(ski + 398);
    const auto *ski_401 = buffer.data(ski + 401);
    const auto *ski_402 = buffer.data(ski + 402);
    const auto *ski_404 = buffer.data(ski + 404);
    const auto *ski_406 = buffer.data(ski + 406);
    const auto *ski_407 = buffer.data(ski + 407);
    const auto *ski_409 = buffer.data(ski + 409);
    const auto *ski_410 = buffer.data(ski + 410);
    const auto *ski_412 = buffer.data(ski + 412);
    const auto *ski_413 = buffer.data(ski + 413);
    const auto *ski_414 = buffer.data(ski + 414);
    const auto *ski_415 = buffer.data(ski + 415);
    const auto *ski_416 = buffer.data(ski + 416);
    const auto *ski_417 = buffer.data(ski + 417);
    const auto *ski_418 = buffer.data(ski + 418);
    const auto *ski_419 = buffer.data(ski + 419);
    const auto *ski_420 = buffer.data(ski + 420);
    const auto *ski_422 = buffer.data(ski + 422);
    const auto *ski_423 = buffer.data(ski + 423);
    const auto *ski_425 = buffer.data(ski + 425);
    const auto *ski_426 = buffer.data(ski + 426);
    const auto *ski_429 = buffer.data(ski + 429);
    const auto *ski_430 = buffer.data(ski + 430);
    const auto *ski_432 = buffer.data(ski + 432);
    const auto *ski_434 = buffer.data(ski + 434);
    const auto *ski_435 = buffer.data(ski + 435);
    const auto *ski_437 = buffer.data(ski + 437);
    const auto *ski_438 = buffer.data(ski + 438);
    const auto *ski_440 = buffer.data(ski + 440);
    const auto *ski_441 = buffer.data(ski + 441);
    const auto *ski_442 = buffer.data(ski + 442);
    const auto *ski_443 = buffer.data(ski + 443);
    const auto *ski_444 = buffer.data(ski + 444);
    const auto *ski_445 = buffer.data(ski + 445);
    const auto *ski_446 = buffer.data(ski + 446);
    const auto *ski_447 = buffer.data(ski + 447);

#pragma omp simd aligned(t_456, t_457, t_458, t_459, pc_x, sii_360, sii_361, sii_362, sii_363, \
                         ski_360, ski_361, ski_362, ski_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_456[k] = f_15 * sii_360[k]
                   + f_3 * pc_x[k] * ski_360[k];

        t_457[k] = f_15 * sii_361[k]
                   + f_3 * pc_x[k] * ski_361[k];

        t_458[k] = f_15 * sii_362[k]
                   + f_3 * pc_x[k] * ski_362[k];

        t_459[k] = f_15 * sii_363[k]
                   + f_3 * pc_x[k] * ski_363[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, pc_y, pc_z, sii_217, sii_245, sii_247, skh0_267, \
                         skh0_269, skh1_267, skh1_269, ski_357, \
                         ski_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = f_14 * sii_245[k]
                   + f_1 * skh0_267[k]
                   - f_2 * skh1_267[k]
                   + f_3 * pc_y[k] * ski_357[k];

        t_461[k] = f_14 * sii_217[k]
                   + f_3 * pc_z[k] * ski_357[k];

        t_462[k] = f_14 * sii_247[k]
                   + f_4 * skh0_269[k]
                   - f_5 * skh1_269[k]
                   + f_3 * pc_y[k] * ski_359[k];
    }

#pragma omp simd aligned(t_463, t_464, t_465, pc_y, sii_248, sii_249, sii_250, skh0_270, \
                         skh0_271, skh0_272, skh1_270, skh1_271, skh1_272, ski_360, ski_361, \
                         ski_362 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_463[k] = f_14 * sii_248[k]
                   + f_6 * skh0_270[k]
                   - f_7 * skh1_270[k]
                   + f_3 * pc_y[k] * ski_360[k];

        t_464[k] = f_14 * sii_249[k]
                   + f_8 * skh0_271[k]
                   - f_9 * skh1_271[k]
                   + f_3 * pc_y[k] * ski_361[k];

        t_465[k] = f_14 * sii_250[k]
                   + f_10 * skh0_272[k]
                   - f_11 * skh1_272[k]
                   + f_3 * pc_y[k] * ski_362[k];
    }

#pragma omp simd aligned(t_466, t_467, t_468, t_469, pb_y, pc_y, pc_z, sik0_324, sii_223, \
                         sii_251, sii_252, sik1_324, skh0_272, skh1_272, ski_363, \
                         ski_364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_466[k] = f_14 * sii_251[k]
                   + f_3 * pc_y[k] * ski_363[k];

        t_467[k] = f_14 * sii_223[k]
                   + f_1 * skh0_272[k]
                   - f_2 * skh1_272[k]
                   + f_3 * pc_z[k] * ski_363[k];

        t_468[k] = pb_y[k] * sik0_324[k]
                   - f_12 * pc_y[k] * sik1_324[k];

        t_469[k] = f_13 * sii_252[k]
                   + f_3 * pc_y[k] * ski_364[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, pb_y, pc_y, pc_z, sik0_327, sik0_329, \
                         sii_224, sii_253, sii_254, sik1_327, sik1_329, ski_364, \
                         ski_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = f_15 * sii_224[k]
                   + f_3 * pc_z[k] * ski_364[k];

        t_471[k] = pb_y[k] * sik0_327[k]
                   + f_14 * sii_253[k]
                   - f_12 * pc_y[k] * sik1_327[k];

        t_472[k] = f_13 * sii_254[k]
                   + f_3 * pc_y[k] * ski_366[k];

        t_473[k] = pb_y[k] * sik0_329[k]
                   - f_12 * pc_y[k] * sik1_329[k];
    }

#pragma omp simd aligned(t_474, t_475, t_476, t_477, pb_y, pc_y, pc_z, sik0_330, sik0_333, \
                         sii_227, sii_255, sii_257, sik1_330, sik1_333, ski_367, \
                         ski_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = pb_y[k] * sik0_330[k]
                   + f_15 * sii_255[k]
                   - f_12 * pc_y[k] * sik1_330[k];

        t_475[k] = f_15 * sii_227[k]
                   + f_3 * pc_z[k] * ski_367[k];

        t_476[k] = f_13 * sii_257[k]
                   + f_3 * pc_y[k] * ski_369[k];

        t_477[k] = pb_y[k] * sik0_333[k]
                   - f_12 * pc_y[k] * sik1_333[k];
    }

#pragma omp simd aligned(t_478, t_479, t_480, pb_y, pc_y, pc_z, sik0_334, sik0_336, sii_230, \
                         sii_258, sii_260, sik1_334, sik1_336, \
                         ski_370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_478[k] = pb_y[k] * sik0_334[k]
                   + f_16 * sii_258[k]
                   - f_12 * pc_y[k] * sik1_334[k];

        t_479[k] = f_15 * sii_230[k]
                   + f_3 * pc_z[k] * ski_370[k];

        t_480[k] = pb_y[k] * sik0_336[k]
                   + f_14 * sii_260[k]
                   - f_12 * pc_y[k] * sik1_336[k];
    }

#pragma omp simd aligned(t_481, t_482, t_483, t_484, pb_y, pc_y, pc_z, sik0_338, sik0_339, \
                         sii_234, sii_261, sii_262, sik1_338, sik1_339, ski_373, \
                         ski_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_481[k] = f_13 * sii_261[k]
                   + f_3 * pc_y[k] * ski_373[k];

        t_482[k] = pb_y[k] * sik0_338[k]
                   - f_12 * pc_y[k] * sik1_338[k];

        t_483[k] = pb_y[k] * sik0_339[k]
                   + f_17 * sii_262[k]
                   - f_12 * pc_y[k] * sik1_339[k];

        t_484[k] = f_15 * sii_234[k]
                   + f_3 * pc_z[k] * ski_374[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, t_488, pb_y, pc_y, sik0_341, sik0_342, sik0_344, \
                         sii_264, sii_265, sii_266, sik1_341, sik1_342, sik1_344, \
                         ski_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = pb_y[k] * sik0_341[k]
                   + f_15 * sii_264[k]
                   - f_12 * pc_y[k] * sik1_341[k];

        t_486[k] = pb_y[k] * sik0_342[k]
                   + f_14 * sii_265[k]
                   - f_12 * pc_y[k] * sik1_342[k];

        t_487[k] = f_13 * sii_266[k]
                   + f_3 * pc_y[k] * ski_378[k];

        t_488[k] = pb_y[k] * sik0_344[k]
                   - f_12 * pc_y[k] * sik1_344[k];
    }

#pragma omp simd aligned(t_489, t_490, t_491, t_492, t_493, pc_x, sii_385, sii_386, sii_387, \
                         sii_388, sii_389, ski_385, ski_386, ski_387, ski_388, \
                         ski_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_489[k] = f_15 * sii_385[k]
                   + f_3 * pc_x[k] * ski_385[k];

        t_490[k] = f_15 * sii_386[k]
                   + f_3 * pc_x[k] * ski_386[k];

        t_491[k] = f_15 * sii_387[k]
                   + f_3 * pc_x[k] * ski_387[k];

        t_492[k] = f_15 * sii_388[k]
                   + f_3 * pc_x[k] * ski_388[k];

        t_493[k] = f_15 * sii_389[k]
                   + f_3 * pc_x[k] * ski_389[k];
    }

#pragma omp simd aligned(t_494, t_495, t_496, t_497, pc_x, pc_y, pc_z, sii_245, sii_273, \
                         sii_390, sii_391, skh0_288, skh1_288, ski_385, ski_390, \
                         ski_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_494[k] = f_15 * sii_390[k]
                   + f_3 * pc_x[k] * ski_390[k];

        t_495[k] = f_15 * sii_391[k]
                   + f_3 * pc_x[k] * ski_391[k];

        t_496[k] = f_13 * sii_273[k]
                   + f_1 * skh0_288[k]
                   - f_2 * skh1_288[k]
                   + f_3 * pc_y[k] * ski_385[k];

        t_497[k] = f_15 * sii_245[k]
                   + f_3 * pc_z[k] * ski_385[k];
    }

#pragma omp simd aligned(t_498, t_499, t_500, pc_y, sii_275, sii_276, sii_277, skh0_290, \
                         skh0_291, skh0_292, skh1_290, skh1_291, skh1_292, ski_387, ski_388, \
                         ski_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_498[k] = f_13 * sii_275[k]
                   + f_4 * skh0_290[k]
                   - f_5 * skh1_290[k]
                   + f_3 * pc_y[k] * ski_387[k];

        t_499[k] = f_13 * sii_276[k]
                   + f_6 * skh0_291[k]
                   - f_7 * skh1_291[k]
                   + f_3 * pc_y[k] * ski_388[k];

        t_500[k] = f_13 * sii_277[k]
                   + f_8 * skh0_292[k]
                   - f_9 * skh1_292[k]
                   + f_3 * pc_y[k] * ski_389[k];
    }

#pragma omp simd aligned(t_501, t_502, t_503, pb_y, pc_y, sik0_359, sii_278, sii_279, \
                         sik1_359, skh0_293, skh1_293, ski_390, \
                         ski_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_501[k] = f_13 * sii_278[k]
                   + f_10 * skh0_293[k]
                   - f_11 * skh1_293[k]
                   + f_3 * pc_y[k] * ski_390[k];

        t_502[k] = f_13 * sii_279[k]
                   + f_3 * pc_y[k] * ski_391[k];

        t_503[k] = pb_y[k] * sik0_359[k]
                   - f_12 * pc_y[k] * sik1_359[k];
    }

#pragma omp simd aligned(t_504, t_505, t_506, t_507, pc_x, pc_y, pc_z, sii_252, sii_392, \
                         sii_395, skh0_294, skh0_297, skh1_294, skh1_297, ski_392, \
                         ski_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_504[k] = f_15 * sii_392[k]
                   + f_1 * skh0_294[k]
                   - f_2 * skh1_294[k]
                   + f_3 * pc_x[k] * ski_392[k];

        t_505[k] = f_3 * pc_y[k] * ski_392[k];

        t_506[k] = f_16 * sii_252[k]
                   + f_3 * pc_z[k] * ski_392[k];

        t_507[k] = f_15 * sii_395[k]
                   + f_4 * skh0_297[k]
                   - f_5 * skh1_297[k]
                   + f_3 * pc_x[k] * ski_395[k];
    }

#pragma omp simd aligned(t_508, t_509, t_510, pc_x, pc_y, sii_397, sii_398, skh0_299, \
                         skh0_300, skh1_299, skh1_300, ski_394, ski_397, \
                         ski_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_508[k] = f_3 * pc_y[k] * ski_394[k];

        t_509[k] = f_15 * sii_397[k]
                   + f_4 * skh0_299[k]
                   - f_5 * skh1_299[k]
                   + f_3 * pc_x[k] * ski_397[k];

        t_510[k] = f_15 * sii_398[k]
                   + f_6 * skh0_300[k]
                   - f_7 * skh1_300[k]
                   + f_3 * pc_x[k] * ski_398[k];
    }

#pragma omp simd aligned(t_511, t_512, t_513, pc_x, pc_y, pc_z, sii_255, sii_401, skh0_303, \
                         skh1_303, ski_395, ski_397, ski_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_511[k] = f_16 * sii_255[k]
                   + f_3 * pc_z[k] * ski_395[k];

        t_512[k] = f_3 * pc_y[k] * ski_397[k];

        t_513[k] = f_15 * sii_401[k]
                   + f_6 * skh0_303[k]
                   - f_7 * skh1_303[k]
                   + f_3 * pc_x[k] * ski_401[k];
    }

#pragma omp simd aligned(t_514, t_515, t_516, pc_x, pc_z, sii_258, sii_402, sii_404, skh0_304, \
                         skh0_306, skh1_304, skh1_306, ski_398, ski_402, \
                         ski_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_514[k] = f_15 * sii_402[k]
                   + f_8 * skh0_304[k]
                   - f_9 * skh1_304[k]
                   + f_3 * pc_x[k] * ski_402[k];

        t_515[k] = f_16 * sii_258[k]
                   + f_3 * pc_z[k] * ski_398[k];

        t_516[k] = f_15 * sii_404[k]
                   + f_8 * skh0_306[k]
                   - f_9 * skh1_306[k]
                   + f_3 * pc_x[k] * ski_404[k];
    }

#pragma omp simd aligned(t_517, t_518, t_519, pc_x, pc_y, sii_406, sii_407, skh0_308, \
                         skh0_309, skh1_308, skh1_309, ski_401, ski_406, \
                         ski_407 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_517[k] = f_3 * pc_y[k] * ski_401[k];

        t_518[k] = f_15 * sii_406[k]
                   + f_8 * skh0_308[k]
                   - f_9 * skh1_308[k]
                   + f_3 * pc_x[k] * ski_406[k];

        t_519[k] = f_15 * sii_407[k]
                   + f_10 * skh0_309[k]
                   - f_11 * skh1_309[k]
                   + f_3 * pc_x[k] * ski_407[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, pc_x, pc_z, sii_262, sii_409, sii_410, skh0_311, \
                         skh0_312, skh1_311, skh1_312, ski_402, ski_409, \
                         ski_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = f_16 * sii_262[k]
                   + f_3 * pc_z[k] * ski_402[k];

        t_521[k] = f_15 * sii_409[k]
                   + f_10 * skh0_311[k]
                   - f_11 * skh1_311[k]
                   + f_3 * pc_x[k] * ski_409[k];

        t_522[k] = f_15 * sii_410[k]
                   + f_10 * skh0_312[k]
                   - f_11 * skh1_312[k]
                   + f_3 * pc_x[k] * ski_410[k];
    }

#pragma omp simd aligned(t_523, t_524, t_525, t_526, pc_x, pc_y, sii_412, sii_413, sii_414, \
                         skh0_314, skh1_314, ski_406, ski_412, ski_413, \
                         ski_414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_523[k] = f_3 * pc_y[k] * ski_406[k];

        t_524[k] = f_15 * sii_412[k]
                   + f_10 * skh0_314[k]
                   - f_11 * skh1_314[k]
                   + f_3 * pc_x[k] * ski_412[k];

        t_525[k] = f_15 * sii_413[k]
                   + f_3 * pc_x[k] * ski_413[k];

        t_526[k] = f_15 * sii_414[k]
                   + f_3 * pc_x[k] * ski_414[k];
    }

#pragma omp simd aligned(t_527, t_528, t_529, t_530, t_531, pc_x, sii_415, sii_416, sii_417, \
                         sii_418, sii_419, ski_415, ski_416, ski_417, ski_418, \
                         ski_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_527[k] = f_15 * sii_415[k]
                   + f_3 * pc_x[k] * ski_415[k];

        t_528[k] = f_15 * sii_416[k]
                   + f_3 * pc_x[k] * ski_416[k];

        t_529[k] = f_15 * sii_417[k]
                   + f_3 * pc_x[k] * ski_417[k];

        t_530[k] = f_15 * sii_418[k]
                   + f_3 * pc_x[k] * ski_418[k];

        t_531[k] = f_15 * sii_419[k]
                   + f_3 * pc_x[k] * ski_419[k];
    }

#pragma omp simd aligned(t_532, t_533, t_534, t_535, pc_y, pc_z, sii_273, skh0_309, skh0_311, \
                         skh0_312, skh1_309, skh1_311, skh1_312, ski_413, ski_415, \
                         ski_416 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_532[k] = f_1 * skh0_309[k]
                   - f_2 * skh1_309[k]
                   + f_3 * pc_y[k] * ski_413[k];

        t_533[k] = f_16 * sii_273[k]
                   + f_3 * pc_z[k] * ski_413[k];

        t_534[k] = f_4 * skh0_311[k]
                   - f_5 * skh1_311[k]
                   + f_3 * pc_y[k] * ski_415[k];

        t_535[k] = f_6 * skh0_312[k]
                   - f_7 * skh1_312[k]
                   + f_3 * pc_y[k] * ski_416[k];
    }

#pragma omp simd aligned(t_536, t_537, t_538, t_539, pc_y, pc_z, sii_279, skh0_313, skh0_314, \
                         skh1_313, skh1_314, ski_417, ski_418, \
                         ski_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_536[k] = f_8 * skh0_313[k]
                   - f_9 * skh1_313[k]
                   + f_3 * pc_y[k] * ski_417[k];

        t_537[k] = f_10 * skh0_314[k]
                   - f_11 * skh1_314[k]
                   + f_3 * pc_y[k] * ski_418[k];

        t_538[k] = f_3 * pc_y[k] * ski_419[k];

        t_539[k] = f_16 * sii_279[k]
                   + f_1 * skh0_314[k]
                   - f_2 * skh1_314[k]
                   + f_3 * pc_z[k] * ski_419[k];
    }

#pragma omp simd aligned(t_540, t_541, t_542, t_543, pc_x, pc_y, pc_z, sii_280, sii_420, \
                         sii_423, skh0_315, skh0_318, skh1_315, skh1_318, ski_420, \
                         ski_423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_540[k] = f_14 * sii_420[k]
                   + f_1 * skh0_315[k]
                   - f_2 * skh1_315[k]
                   + f_3 * pc_x[k] * ski_420[k];

        t_541[k] = f_17 * sii_280[k]
                   + f_3 * pc_y[k] * ski_420[k];

        t_542[k] = f_3 * pc_z[k] * ski_420[k];

        t_543[k] = f_14 * sii_423[k]
                   + f_4 * skh0_318[k]
                   - f_5 * skh1_318[k]
                   + f_3 * pc_x[k] * ski_423[k];
    }

#pragma omp simd aligned(t_544, t_545, t_546, pc_x, pc_y, sii_282, sii_425, sii_426, skh0_320, \
                         skh0_321, skh1_320, skh1_321, ski_422, ski_425, \
                         ski_426 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_544[k] = f_17 * sii_282[k]
                   + f_3 * pc_y[k] * ski_422[k];

        t_545[k] = f_14 * sii_425[k]
                   + f_4 * skh0_320[k]
                   - f_5 * skh1_320[k]
                   + f_3 * pc_x[k] * ski_425[k];

        t_546[k] = f_14 * sii_426[k]
                   + f_6 * skh0_321[k]
                   - f_7 * skh1_321[k]
                   + f_3 * pc_x[k] * ski_426[k];
    }

#pragma omp simd aligned(t_547, t_548, t_549, pc_x, pc_y, pc_z, sii_285, sii_429, skh0_324, \
                         skh1_324, ski_423, ski_425, ski_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_547[k] = f_3 * pc_z[k] * ski_423[k];

        t_548[k] = f_17 * sii_285[k]
                   + f_3 * pc_y[k] * ski_425[k];

        t_549[k] = f_14 * sii_429[k]
                   + f_6 * skh0_324[k]
                   - f_7 * skh1_324[k]
                   + f_3 * pc_x[k] * ski_429[k];
    }

#pragma omp simd aligned(t_550, t_551, t_552, pc_x, pc_z, sii_430, sii_432, skh0_325, \
                         skh0_327, skh1_325, skh1_327, ski_426, ski_430, \
                         ski_432 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_550[k] = f_14 * sii_430[k]
                   + f_8 * skh0_325[k]
                   - f_9 * skh1_325[k]
                   + f_3 * pc_x[k] * ski_430[k];

        t_551[k] = f_3 * pc_z[k] * ski_426[k];

        t_552[k] = f_14 * sii_432[k]
                   + f_8 * skh0_327[k]
                   - f_9 * skh1_327[k]
                   + f_3 * pc_x[k] * ski_432[k];
    }

#pragma omp simd aligned(t_553, t_554, t_555, pc_x, pc_y, sii_289, sii_434, sii_435, skh0_329, \
                         skh0_330, skh1_329, skh1_330, ski_429, ski_434, \
                         ski_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_553[k] = f_17 * sii_289[k]
                   + f_3 * pc_y[k] * ski_429[k];

        t_554[k] = f_14 * sii_434[k]
                   + f_8 * skh0_329[k]
                   - f_9 * skh1_329[k]
                   + f_3 * pc_x[k] * ski_434[k];

        t_555[k] = f_14 * sii_435[k]
                   + f_10 * skh0_330[k]
                   - f_11 * skh1_330[k]
                   + f_3 * pc_x[k] * ski_435[k];
    }

#pragma omp simd aligned(t_556, t_557, t_558, pc_x, pc_z, sii_437, sii_438, skh0_332, \
                         skh0_333, skh1_332, skh1_333, ski_430, ski_437, \
                         ski_438 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_556[k] = f_3 * pc_z[k] * ski_430[k];

        t_557[k] = f_14 * sii_437[k]
                   + f_10 * skh0_332[k]
                   - f_11 * skh1_332[k]
                   + f_3 * pc_x[k] * ski_437[k];

        t_558[k] = f_14 * sii_438[k]
                   + f_10 * skh0_333[k]
                   - f_11 * skh1_333[k]
                   + f_3 * pc_x[k] * ski_438[k];
    }

#pragma omp simd aligned(t_559, t_560, t_561, t_562, pc_x, pc_y, sii_294, sii_440, sii_441, \
                         sii_442, skh0_335, skh1_335, ski_434, ski_440, ski_441, \
                         ski_442 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_559[k] = f_17 * sii_294[k]
                   + f_3 * pc_y[k] * ski_434[k];

        t_560[k] = f_14 * sii_440[k]
                   + f_10 * skh0_335[k]
                   - f_11 * skh1_335[k]
                   + f_3 * pc_x[k] * ski_440[k];

        t_561[k] = f_14 * sii_441[k]
                   + f_3 * pc_x[k] * ski_441[k];

        t_562[k] = f_14 * sii_442[k]
                   + f_3 * pc_x[k] * ski_442[k];
    }

#pragma omp simd aligned(t_563, t_564, t_565, t_566, t_567, pc_x, sii_443, sii_444, sii_445, \
                         sii_446, sii_447, ski_443, ski_444, ski_445, ski_446, \
                         ski_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_563[k] = f_14 * sii_443[k]
                   + f_3 * pc_x[k] * ski_443[k];

        t_564[k] = f_14 * sii_444[k]
                   + f_3 * pc_x[k] * ski_444[k];

        t_565[k] = f_14 * sii_445[k]
                   + f_3 * pc_x[k] * ski_445[k];

        t_566[k] = f_14 * sii_446[k]
                   + f_3 * pc_x[k] * ski_446[k];

        t_567[k] = f_14 * sii_447[k]
                   + f_3 * pc_x[k] * ski_447[k];
    }

#pragma omp simd aligned(t_568, t_569, t_570, pc_y, pc_z, sii_301, sii_303, skh0_330, \
                         skh0_332, skh1_330, skh1_332, ski_441, \
                         ski_443 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_568[k] = f_17 * sii_301[k]
                   + f_1 * skh0_330[k]
                   - f_2 * skh1_330[k]
                   + f_3 * pc_y[k] * ski_441[k];

        t_569[k] = f_3 * pc_z[k] * ski_441[k];

        t_570[k] = f_17 * sii_303[k]
                   + f_4 * skh0_332[k]
                   - f_5 * skh1_332[k]
                   + f_3 * pc_y[k] * ski_443[k];
    }
}

static auto
compute_prim_skk_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sik0,
                                                          const size_t sii, const size_t sik1,
                                                          const size_t skh0, const size_t skh1,
                                                          const size_t ski, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.0 / gamma;
    const auto f_2 = 3.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.0 / gamma;
    const auto f_5 = 2.0 * p / (gamma * q);
    const auto f_6 = 1.5 / gamma;
    const auto f_7 = 1.5 * p / (gamma * q);
    const auto f_8 = 1.0 / gamma;
    const auto f_9 = p / (gamma * q);
    const auto f_10 = 0.5 / gamma;
    const auto f_11 = 0.5 * p / (gamma * q);
    const auto f_12 = gamma / q;
    const auto f_13 = 0.5 / q;
    const auto f_14 = 1.0 / q;
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;

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
    auto *t_674 = buffer.data(target + 674);
    auto *t_675 = buffer.data(target + 675);
    auto *t_676 = buffer.data(target + 676);
    auto *t_677 = buffer.data(target + 677);
    auto *t_678 = buffer.data(target + 678);

    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sik0_360 = buffer.data(sik0 + 360);
    const auto *sik0_363 = buffer.data(sik0 + 363);
    const auto *sik0_366 = buffer.data(sik0 + 366);
    const auto *sik0_370 = buffer.data(sik0 + 370);
    const auto *sik0_372 = buffer.data(sik0 + 372);
    const auto *sik0_375 = buffer.data(sik0 + 375);
    const auto *sik0_377 = buffer.data(sik0 + 377);
    const auto *sik0_378 = buffer.data(sik0 + 378);
    const auto *sik0_388 = buffer.data(sik0 + 388);

    const auto *sii_280 = buffer.data(sii + 280);
    const auto *sii_283 = buffer.data(sii + 283);
    const auto *sii_286 = buffer.data(sii + 286);
    const auto *sii_287 = buffer.data(sii + 287);
    const auto *sii_290 = buffer.data(sii + 290);
    const auto *sii_291 = buffer.data(sii + 291);
    const auto *sii_292 = buffer.data(sii + 292);
    const auto *sii_301 = buffer.data(sii + 301);
    const auto *sii_304 = buffer.data(sii + 304);
    const auto *sii_305 = buffer.data(sii + 305);
    const auto *sii_306 = buffer.data(sii + 306);
    const auto *sii_307 = buffer.data(sii + 307);
    const auto *sii_308 = buffer.data(sii + 308);
    const auto *sii_310 = buffer.data(sii + 310);
    const auto *sii_311 = buffer.data(sii + 311);
    const auto *sii_313 = buffer.data(sii + 313);
    const auto *sii_314 = buffer.data(sii + 314);
    const auto *sii_317 = buffer.data(sii + 317);
    const auto *sii_318 = buffer.data(sii + 318);
    const auto *sii_322 = buffer.data(sii + 322);
    const auto *sii_329 = buffer.data(sii + 329);
    const auto *sii_331 = buffer.data(sii + 331);
    const auto *sii_332 = buffer.data(sii + 332);
    const auto *sii_333 = buffer.data(sii + 333);
    const auto *sii_334 = buffer.data(sii + 334);
    const auto *sii_335 = buffer.data(sii + 335);
    const auto *sii_336 = buffer.data(sii + 336);
    const auto *sii_338 = buffer.data(sii + 338);
    const auto *sii_339 = buffer.data(sii + 339);
    const auto *sii_341 = buffer.data(sii + 341);
    const auto *sii_342 = buffer.data(sii + 342);
    const auto *sii_345 = buffer.data(sii + 345);
    const auto *sii_346 = buffer.data(sii + 346);
    const auto *sii_350 = buffer.data(sii + 350);
    const auto *sii_357 = buffer.data(sii + 357);
    const auto *sii_359 = buffer.data(sii + 359);
    const auto *sii_360 = buffer.data(sii + 360);
    const auto *sii_361 = buffer.data(sii + 361);
    const auto *sii_362 = buffer.data(sii + 362);
    const auto *sii_363 = buffer.data(sii + 363);
    const auto *sii_364 = buffer.data(sii + 364);
    const auto *sii_366 = buffer.data(sii + 366);
    const auto *sii_369 = buffer.data(sii + 369);
    const auto *sii_373 = buffer.data(sii + 373);
    const auto *sii_378 = buffer.data(sii + 378);
    const auto *sii_385 = buffer.data(sii + 385);
    const auto *sii_387 = buffer.data(sii + 387);
    const auto *sii_453 = buffer.data(sii + 453);
    const auto *sii_457 = buffer.data(sii + 457);
    const auto *sii_462 = buffer.data(sii + 462);
    const auto *sii_468 = buffer.data(sii + 468);
    const auto *sii_469 = buffer.data(sii + 469);
    const auto *sii_470 = buffer.data(sii + 470);
    const auto *sii_471 = buffer.data(sii + 471);
    const auto *sii_472 = buffer.data(sii + 472);
    const auto *sii_473 = buffer.data(sii + 473);
    const auto *sii_474 = buffer.data(sii + 474);
    const auto *sii_475 = buffer.data(sii + 475);
    const auto *sii_476 = buffer.data(sii + 476);
    const auto *sii_479 = buffer.data(sii + 479);
    const auto *sii_481 = buffer.data(sii + 481);
    const auto *sii_482 = buffer.data(sii + 482);
    const auto *sii_485 = buffer.data(sii + 485);
    const auto *sii_486 = buffer.data(sii + 486);
    const auto *sii_488 = buffer.data(sii + 488);
    const auto *sii_490 = buffer.data(sii + 490);
    const auto *sii_491 = buffer.data(sii + 491);
    const auto *sii_493 = buffer.data(sii + 493);
    const auto *sii_494 = buffer.data(sii + 494);
    const auto *sii_496 = buffer.data(sii + 496);
    const auto *sii_497 = buffer.data(sii + 497);
    const auto *sii_498 = buffer.data(sii + 498);
    const auto *sii_499 = buffer.data(sii + 499);
    const auto *sii_500 = buffer.data(sii + 500);
    const auto *sii_501 = buffer.data(sii + 501);
    const auto *sii_502 = buffer.data(sii + 502);
    const auto *sii_503 = buffer.data(sii + 503);
    const auto *sii_504 = buffer.data(sii + 504);
    const auto *sii_507 = buffer.data(sii + 507);
    const auto *sii_509 = buffer.data(sii + 509);
    const auto *sii_510 = buffer.data(sii + 510);
    const auto *sii_513 = buffer.data(sii + 513);
    const auto *sii_514 = buffer.data(sii + 514);
    const auto *sii_516 = buffer.data(sii + 516);
    const auto *sii_518 = buffer.data(sii + 518);
    const auto *sii_519 = buffer.data(sii + 519);
    const auto *sii_521 = buffer.data(sii + 521);
    const auto *sii_522 = buffer.data(sii + 522);
    const auto *sii_524 = buffer.data(sii + 524);
    const auto *sii_525 = buffer.data(sii + 525);
    const auto *sii_526 = buffer.data(sii + 526);
    const auto *sii_527 = buffer.data(sii + 527);
    const auto *sii_528 = buffer.data(sii + 528);
    const auto *sii_529 = buffer.data(sii + 529);
    const auto *sii_530 = buffer.data(sii + 530);
    const auto *sii_531 = buffer.data(sii + 531);

    const auto *sik1_360 = buffer.data(sik1 + 360);
    const auto *sik1_363 = buffer.data(sik1 + 363);
    const auto *sik1_366 = buffer.data(sik1 + 366);
    const auto *sik1_370 = buffer.data(sik1 + 370);
    const auto *sik1_372 = buffer.data(sik1 + 372);
    const auto *sik1_375 = buffer.data(sik1 + 375);
    const auto *sik1_377 = buffer.data(sik1 + 377);
    const auto *sik1_378 = buffer.data(sik1 + 378);
    const auto *sik1_388 = buffer.data(sik1 + 388);

    const auto *skh0_333 = buffer.data(skh0 + 333);
    const auto *skh0_334 = buffer.data(skh0 + 334);
    const auto *skh0_335 = buffer.data(skh0 + 335);
    const auto *skh0_341 = buffer.data(skh0 + 341);
    const auto *skh0_345 = buffer.data(skh0 + 345);
    const auto *skh0_350 = buffer.data(skh0 + 350);
    const auto *skh0_353 = buffer.data(skh0 + 353);
    const auto *skh0_354 = buffer.data(skh0 + 354);
    const auto *skh0_355 = buffer.data(skh0 + 355);
    const auto *skh0_356 = buffer.data(skh0 + 356);
    const auto *skh0_357 = buffer.data(skh0 + 357);
    const auto *skh0_360 = buffer.data(skh0 + 360);
    const auto *skh0_362 = buffer.data(skh0 + 362);
    const auto *skh0_363 = buffer.data(skh0 + 363);
    const auto *skh0_366 = buffer.data(skh0 + 366);
    const auto *skh0_367 = buffer.data(skh0 + 367);
    const auto *skh0_369 = buffer.data(skh0 + 369);
    const auto *skh0_371 = buffer.data(skh0 + 371);
    const auto *skh0_372 = buffer.data(skh0 + 372);
    const auto *skh0_374 = buffer.data(skh0 + 374);
    const auto *skh0_375 = buffer.data(skh0 + 375);
    const auto *skh0_376 = buffer.data(skh0 + 376);
    const auto *skh0_377 = buffer.data(skh0 + 377);
    const auto *skh0_378 = buffer.data(skh0 + 378);
    const auto *skh0_381 = buffer.data(skh0 + 381);
    const auto *skh0_383 = buffer.data(skh0 + 383);
    const auto *skh0_384 = buffer.data(skh0 + 384);
    const auto *skh0_387 = buffer.data(skh0 + 387);
    const auto *skh0_388 = buffer.data(skh0 + 388);
    const auto *skh0_390 = buffer.data(skh0 + 390);
    const auto *skh0_392 = buffer.data(skh0 + 392);
    const auto *skh0_393 = buffer.data(skh0 + 393);
    const auto *skh0_395 = buffer.data(skh0 + 395);
    const auto *skh0_396 = buffer.data(skh0 + 396);
    const auto *skh0_398 = buffer.data(skh0 + 398);

    const auto *skh1_333 = buffer.data(skh1 + 333);
    const auto *skh1_334 = buffer.data(skh1 + 334);
    const auto *skh1_335 = buffer.data(skh1 + 335);
    const auto *skh1_341 = buffer.data(skh1 + 341);
    const auto *skh1_345 = buffer.data(skh1 + 345);
    const auto *skh1_350 = buffer.data(skh1 + 350);
    const auto *skh1_353 = buffer.data(skh1 + 353);
    const auto *skh1_354 = buffer.data(skh1 + 354);
    const auto *skh1_355 = buffer.data(skh1 + 355);
    const auto *skh1_356 = buffer.data(skh1 + 356);
    const auto *skh1_357 = buffer.data(skh1 + 357);
    const auto *skh1_360 = buffer.data(skh1 + 360);
    const auto *skh1_362 = buffer.data(skh1 + 362);
    const auto *skh1_363 = buffer.data(skh1 + 363);
    const auto *skh1_366 = buffer.data(skh1 + 366);
    const auto *skh1_367 = buffer.data(skh1 + 367);
    const auto *skh1_369 = buffer.data(skh1 + 369);
    const auto *skh1_371 = buffer.data(skh1 + 371);
    const auto *skh1_372 = buffer.data(skh1 + 372);
    const auto *skh1_374 = buffer.data(skh1 + 374);
    const auto *skh1_375 = buffer.data(skh1 + 375);
    const auto *skh1_376 = buffer.data(skh1 + 376);
    const auto *skh1_377 = buffer.data(skh1 + 377);
    const auto *skh1_378 = buffer.data(skh1 + 378);
    const auto *skh1_381 = buffer.data(skh1 + 381);
    const auto *skh1_383 = buffer.data(skh1 + 383);
    const auto *skh1_384 = buffer.data(skh1 + 384);
    const auto *skh1_387 = buffer.data(skh1 + 387);
    const auto *skh1_388 = buffer.data(skh1 + 388);
    const auto *skh1_390 = buffer.data(skh1 + 390);
    const auto *skh1_392 = buffer.data(skh1 + 392);
    const auto *skh1_393 = buffer.data(skh1 + 393);
    const auto *skh1_395 = buffer.data(skh1 + 395);
    const auto *skh1_396 = buffer.data(skh1 + 396);
    const auto *skh1_398 = buffer.data(skh1 + 398);

    const auto *ski_444 = buffer.data(ski + 444);
    const auto *ski_445 = buffer.data(ski + 445);
    const auto *ski_446 = buffer.data(ski + 446);
    const auto *ski_447 = buffer.data(ski + 447);
    const auto *ski_448 = buffer.data(ski + 448);
    const auto *ski_450 = buffer.data(ski + 450);
    const auto *ski_451 = buffer.data(ski + 451);
    const auto *ski_453 = buffer.data(ski + 453);
    const auto *ski_454 = buffer.data(ski + 454);
    const auto *ski_457 = buffer.data(ski + 457);
    const auto *ski_458 = buffer.data(ski + 458);
    const auto *ski_462 = buffer.data(ski + 462);
    const auto *ski_468 = buffer.data(ski + 468);
    const auto *ski_469 = buffer.data(ski + 469);
    const auto *ski_470 = buffer.data(ski + 470);
    const auto *ski_471 = buffer.data(ski + 471);
    const auto *ski_472 = buffer.data(ski + 472);
    const auto *ski_473 = buffer.data(ski + 473);
    const auto *ski_474 = buffer.data(ski + 474);
    const auto *ski_475 = buffer.data(ski + 475);
    const auto *ski_476 = buffer.data(ski + 476);
    const auto *ski_478 = buffer.data(ski + 478);
    const auto *ski_479 = buffer.data(ski + 479);
    const auto *ski_481 = buffer.data(ski + 481);
    const auto *ski_482 = buffer.data(ski + 482);
    const auto *ski_485 = buffer.data(ski + 485);
    const auto *ski_486 = buffer.data(ski + 486);
    const auto *ski_488 = buffer.data(ski + 488);
    const auto *ski_490 = buffer.data(ski + 490);
    const auto *ski_491 = buffer.data(ski + 491);
    const auto *ski_493 = buffer.data(ski + 493);
    const auto *ski_494 = buffer.data(ski + 494);
    const auto *ski_496 = buffer.data(ski + 496);
    const auto *ski_497 = buffer.data(ski + 497);
    const auto *ski_498 = buffer.data(ski + 498);
    const auto *ski_499 = buffer.data(ski + 499);
    const auto *ski_500 = buffer.data(ski + 500);
    const auto *ski_501 = buffer.data(ski + 501);
    const auto *ski_502 = buffer.data(ski + 502);
    const auto *ski_503 = buffer.data(ski + 503);
    const auto *ski_504 = buffer.data(ski + 504);
    const auto *ski_506 = buffer.data(ski + 506);
    const auto *ski_507 = buffer.data(ski + 507);
    const auto *ski_509 = buffer.data(ski + 509);
    const auto *ski_510 = buffer.data(ski + 510);
    const auto *ski_513 = buffer.data(ski + 513);
    const auto *ski_514 = buffer.data(ski + 514);
    const auto *ski_516 = buffer.data(ski + 516);
    const auto *ski_518 = buffer.data(ski + 518);
    const auto *ski_519 = buffer.data(ski + 519);
    const auto *ski_521 = buffer.data(ski + 521);
    const auto *ski_522 = buffer.data(ski + 522);
    const auto *ski_524 = buffer.data(ski + 524);
    const auto *ski_525 = buffer.data(ski + 525);
    const auto *ski_526 = buffer.data(ski + 526);
    const auto *ski_527 = buffer.data(ski + 527);
    const auto *ski_528 = buffer.data(ski + 528);
    const auto *ski_529 = buffer.data(ski + 529);
    const auto *ski_530 = buffer.data(ski + 530);
    const auto *ski_531 = buffer.data(ski + 531);

#pragma omp simd aligned(t_571, t_572, t_573, pc_y, sii_304, sii_305, sii_306, skh0_333, \
                         skh0_334, skh0_335, skh1_333, skh1_334, skh1_335, ski_444, ski_445, \
                         ski_446 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_571[k] = f_17 * sii_304[k]
                   + f_6 * skh0_333[k]
                   - f_7 * skh1_333[k]
                   + f_3 * pc_y[k] * ski_444[k];

        t_572[k] = f_17 * sii_305[k]
                   + f_8 * skh0_334[k]
                   - f_9 * skh1_334[k]
                   + f_3 * pc_y[k] * ski_445[k];

        t_573[k] = f_17 * sii_306[k]
                   + f_10 * skh0_335[k]
                   - f_11 * skh1_335[k]
                   + f_3 * pc_y[k] * ski_446[k];
    }

#pragma omp simd aligned(t_574, t_575, t_576, t_577, pb_z, pc_y, pc_z, sik0_360, sii_307, \
                         sii_308, sik1_360, skh0_335, skh1_335, ski_447, \
                         ski_448 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_574[k] = f_17 * sii_307[k]
                   + f_3 * pc_y[k] * ski_447[k];

        t_575[k] = f_1 * skh0_335[k]
                   - f_2 * skh1_335[k]
                   + f_3 * pc_z[k] * ski_447[k];

        t_576[k] = pb_z[k] * sik0_360[k]
                   - f_12 * pc_z[k] * sik1_360[k];

        t_577[k] = f_16 * sii_308[k]
                   + f_3 * pc_y[k] * ski_448[k];
    }

#pragma omp simd aligned(t_578, t_579, t_580, pb_z, pc_y, pc_z, sik0_363, sii_280, sii_310, \
                         sik1_363, ski_448, ski_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_578[k] = f_13 * sii_280[k]
                   + f_3 * pc_z[k] * ski_448[k];

        t_579[k] = pb_z[k] * sik0_363[k]
                   - f_12 * pc_z[k] * sik1_363[k];

        t_580[k] = f_16 * sii_310[k]
                   + f_3 * pc_y[k] * ski_450[k];
    }

#pragma omp simd aligned(t_581, t_582, t_583, pb_z, pc_x, pc_z, sik0_366, sii_283, sii_453, \
                         sik1_366, skh0_341, skh1_341, ski_451, \
                         ski_453 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_581[k] = f_14 * sii_453[k]
                   + f_4 * skh0_341[k]
                   - f_5 * skh1_341[k]
                   + f_3 * pc_x[k] * ski_453[k];

        t_582[k] = pb_z[k] * sik0_366[k]
                   - f_12 * pc_z[k] * sik1_366[k];

        t_583[k] = f_13 * sii_283[k]
                   + f_3 * pc_z[k] * ski_451[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, pb_z, pc_x, pc_y, pc_z, sik0_370, sii_313, \
                         sii_457, sik1_370, skh0_345, skh1_345, ski_453, \
                         ski_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = f_16 * sii_313[k]
                   + f_3 * pc_y[k] * ski_453[k];

        t_585[k] = f_14 * sii_457[k]
                   + f_6 * skh0_345[k]
                   - f_7 * skh1_345[k]
                   + f_3 * pc_x[k] * ski_457[k];

        t_586[k] = pb_z[k] * sik0_370[k]
                   - f_12 * pc_z[k] * sik1_370[k];
    }

#pragma omp simd aligned(t_587, t_588, t_589, pb_z, pc_y, pc_z, sik0_372, sii_286, sii_287, \
                         sii_317, sik1_372, ski_454, ski_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_587[k] = f_13 * sii_286[k]
                   + f_3 * pc_z[k] * ski_454[k];

        t_588[k] = pb_z[k] * sik0_372[k]
                   + f_14 * sii_287[k]
                   - f_12 * pc_z[k] * sik1_372[k];

        t_589[k] = f_16 * sii_317[k]
                   + f_3 * pc_y[k] * ski_457[k];
    }

#pragma omp simd aligned(t_590, t_591, t_592, pb_z, pc_x, pc_z, sik0_375, sii_290, sii_462, \
                         sik1_375, skh0_350, skh1_350, ski_458, \
                         ski_462 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_590[k] = f_14 * sii_462[k]
                   + f_8 * skh0_350[k]
                   - f_9 * skh1_350[k]
                   + f_3 * pc_x[k] * ski_462[k];

        t_591[k] = pb_z[k] * sik0_375[k]
                   - f_12 * pc_z[k] * sik1_375[k];

        t_592[k] = f_13 * sii_290[k]
                   + f_3 * pc_z[k] * ski_458[k];
    }

#pragma omp simd aligned(t_593, t_594, t_595, pb_z, pc_y, pc_z, sik0_377, sik0_378, sii_291, \
                         sii_292, sii_322, sik1_377, sik1_378, \
                         ski_462 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_593[k] = pb_z[k] * sik0_377[k]
                   + f_14 * sii_291[k]
                   - f_12 * pc_z[k] * sik1_377[k];

        t_594[k] = pb_z[k] * sik0_378[k]
                   + f_15 * sii_292[k]
                   - f_12 * pc_z[k] * sik1_378[k];

        t_595[k] = f_16 * sii_322[k]
                   + f_3 * pc_y[k] * ski_462[k];
    }

#pragma omp simd aligned(t_596, t_597, t_598, t_599, pc_x, sii_468, sii_469, sii_470, sii_471, \
                         skh0_356, skh1_356, ski_468, ski_469, ski_470, \
                         ski_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_596[k] = f_14 * sii_468[k]
                   + f_10 * skh0_356[k]
                   - f_11 * skh1_356[k]
                   + f_3 * pc_x[k] * ski_468[k];

        t_597[k] = f_14 * sii_469[k]
                   + f_3 * pc_x[k] * ski_469[k];

        t_598[k] = f_14 * sii_470[k]
                   + f_3 * pc_x[k] * ski_470[k];

        t_599[k] = f_14 * sii_471[k]
                   + f_3 * pc_x[k] * ski_471[k];
    }

#pragma omp simd aligned(t_600, t_601, t_602, t_603, pc_x, sii_472, sii_473, sii_474, sii_475, \
                         ski_472, ski_473, ski_474, ski_475 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = f_14 * sii_472[k]
                   + f_3 * pc_x[k] * ski_472[k];

        t_601[k] = f_14 * sii_473[k]
                   + f_3 * pc_x[k] * ski_473[k];

        t_602[k] = f_14 * sii_474[k]
                   + f_3 * pc_x[k] * ski_474[k];

        t_603[k] = f_14 * sii_475[k]
                   + f_3 * pc_x[k] * ski_475[k];
    }

#pragma omp simd aligned(t_604, t_605, t_606, pb_z, pc_y, pc_z, sik0_388, sii_301, sii_331, \
                         sik1_388, skh0_353, skh1_353, ski_469, \
                         ski_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_604[k] = pb_z[k] * sik0_388[k]
                   - f_12 * pc_z[k] * sik1_388[k];

        t_605[k] = f_13 * sii_301[k]
                   + f_3 * pc_z[k] * ski_469[k];

        t_606[k] = f_16 * sii_331[k]
                   + f_4 * skh0_353[k]
                   - f_5 * skh1_353[k]
                   + f_3 * pc_y[k] * ski_471[k];
    }

#pragma omp simd aligned(t_607, t_608, t_609, pc_y, sii_332, sii_333, sii_334, skh0_354, \
                         skh0_355, skh0_356, skh1_354, skh1_355, skh1_356, ski_472, ski_473, \
                         ski_474 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_607[k] = f_16 * sii_332[k]
                   + f_6 * skh0_354[k]
                   - f_7 * skh1_354[k]
                   + f_3 * pc_y[k] * ski_472[k];

        t_608[k] = f_16 * sii_333[k]
                   + f_8 * skh0_355[k]
                   - f_9 * skh1_355[k]
                   + f_3 * pc_y[k] * ski_473[k];

        t_609[k] = f_16 * sii_334[k]
                   + f_10 * skh0_356[k]
                   - f_11 * skh1_356[k]
                   + f_3 * pc_y[k] * ski_474[k];
    }

#pragma omp simd aligned(t_610, t_611, t_612, pc_x, pc_y, pc_z, sii_307, sii_335, sii_476, \
                         skh0_356, skh0_357, skh1_356, skh1_357, ski_475, \
                         ski_476 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_610[k] = f_16 * sii_335[k]
                   + f_3 * pc_y[k] * ski_475[k];

        t_611[k] = f_13 * sii_307[k]
                   + f_1 * skh0_356[k]
                   - f_2 * skh1_356[k]
                   + f_3 * pc_z[k] * ski_475[k];

        t_612[k] = f_14 * sii_476[k]
                   + f_1 * skh0_357[k]
                   - f_2 * skh1_357[k]
                   + f_3 * pc_x[k] * ski_476[k];
    }

#pragma omp simd aligned(t_613, t_614, t_615, t_616, pc_x, pc_y, pc_z, sii_308, sii_336, \
                         sii_338, sii_479, skh0_360, skh1_360, ski_476, ski_478, \
                         ski_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_613[k] = f_15 * sii_336[k]
                   + f_3 * pc_y[k] * ski_476[k];

        t_614[k] = f_14 * sii_308[k]
                   + f_3 * pc_z[k] * ski_476[k];

        t_615[k] = f_14 * sii_479[k]
                   + f_4 * skh0_360[k]
                   - f_5 * skh1_360[k]
                   + f_3 * pc_x[k] * ski_479[k];

        t_616[k] = f_15 * sii_338[k]
                   + f_3 * pc_y[k] * ski_478[k];
    }

#pragma omp simd aligned(t_617, t_618, t_619, pc_x, pc_z, sii_311, sii_481, sii_482, skh0_362, \
                         skh0_363, skh1_362, skh1_363, ski_479, ski_481, \
                         ski_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_617[k] = f_14 * sii_481[k]
                   + f_4 * skh0_362[k]
                   - f_5 * skh1_362[k]
                   + f_3 * pc_x[k] * ski_481[k];

        t_618[k] = f_14 * sii_482[k]
                   + f_6 * skh0_363[k]
                   - f_7 * skh1_363[k]
                   + f_3 * pc_x[k] * ski_482[k];

        t_619[k] = f_14 * sii_311[k]
                   + f_3 * pc_z[k] * ski_479[k];
    }

#pragma omp simd aligned(t_620, t_621, t_622, pc_x, pc_y, sii_341, sii_485, sii_486, skh0_366, \
                         skh0_367, skh1_366, skh1_367, ski_481, ski_485, \
                         ski_486 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_620[k] = f_15 * sii_341[k]
                   + f_3 * pc_y[k] * ski_481[k];

        t_621[k] = f_14 * sii_485[k]
                   + f_6 * skh0_366[k]
                   - f_7 * skh1_366[k]
                   + f_3 * pc_x[k] * ski_485[k];

        t_622[k] = f_14 * sii_486[k]
                   + f_8 * skh0_367[k]
                   - f_9 * skh1_367[k]
                   + f_3 * pc_x[k] * ski_486[k];
    }

#pragma omp simd aligned(t_623, t_624, t_625, pc_x, pc_y, pc_z, sii_314, sii_345, sii_488, \
                         skh0_369, skh1_369, ski_482, ski_485, \
                         ski_488 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_623[k] = f_14 * sii_314[k]
                   + f_3 * pc_z[k] * ski_482[k];

        t_624[k] = f_14 * sii_488[k]
                   + f_8 * skh0_369[k]
                   - f_9 * skh1_369[k]
                   + f_3 * pc_x[k] * ski_488[k];

        t_625[k] = f_15 * sii_345[k]
                   + f_3 * pc_y[k] * ski_485[k];
    }

#pragma omp simd aligned(t_626, t_627, t_628, pc_x, pc_z, sii_318, sii_490, sii_491, skh0_371, \
                         skh0_372, skh1_371, skh1_372, ski_486, ski_490, \
                         ski_491 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_626[k] = f_14 * sii_490[k]
                   + f_8 * skh0_371[k]
                   - f_9 * skh1_371[k]
                   + f_3 * pc_x[k] * ski_490[k];

        t_627[k] = f_14 * sii_491[k]
                   + f_10 * skh0_372[k]
                   - f_11 * skh1_372[k]
                   + f_3 * pc_x[k] * ski_491[k];

        t_628[k] = f_14 * sii_318[k]
                   + f_3 * pc_z[k] * ski_486[k];
    }

#pragma omp simd aligned(t_629, t_630, t_631, pc_x, pc_y, sii_350, sii_493, sii_494, skh0_374, \
                         skh0_375, skh1_374, skh1_375, ski_490, ski_493, \
                         ski_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_629[k] = f_14 * sii_493[k]
                   + f_10 * skh0_374[k]
                   - f_11 * skh1_374[k]
                   + f_3 * pc_x[k] * ski_493[k];

        t_630[k] = f_14 * sii_494[k]
                   + f_10 * skh0_375[k]
                   - f_11 * skh1_375[k]
                   + f_3 * pc_x[k] * ski_494[k];

        t_631[k] = f_15 * sii_350[k]
                   + f_3 * pc_y[k] * ski_490[k];
    }

#pragma omp simd aligned(t_632, t_633, t_634, t_635, pc_x, sii_496, sii_497, sii_498, sii_499, \
                         skh0_377, skh1_377, ski_496, ski_497, ski_498, \
                         ski_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_632[k] = f_14 * sii_496[k]
                   + f_10 * skh0_377[k]
                   - f_11 * skh1_377[k]
                   + f_3 * pc_x[k] * ski_496[k];

        t_633[k] = f_14 * sii_497[k]
                   + f_3 * pc_x[k] * ski_497[k];

        t_634[k] = f_14 * sii_498[k]
                   + f_3 * pc_x[k] * ski_498[k];

        t_635[k] = f_14 * sii_499[k]
                   + f_3 * pc_x[k] * ski_499[k];
    }

#pragma omp simd aligned(t_636, t_637, t_638, t_639, pc_x, sii_500, sii_501, sii_502, sii_503, \
                         ski_500, ski_501, ski_502, ski_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_636[k] = f_14 * sii_500[k]
                   + f_3 * pc_x[k] * ski_500[k];

        t_637[k] = f_14 * sii_501[k]
                   + f_3 * pc_x[k] * ski_501[k];

        t_638[k] = f_14 * sii_502[k]
                   + f_3 * pc_x[k] * ski_502[k];

        t_639[k] = f_14 * sii_503[k]
                   + f_3 * pc_x[k] * ski_503[k];
    }

#pragma omp simd aligned(t_640, t_641, t_642, pc_y, pc_z, sii_329, sii_357, sii_359, skh0_372, \
                         skh0_374, skh1_372, skh1_374, ski_497, \
                         ski_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_640[k] = f_15 * sii_357[k]
                   + f_1 * skh0_372[k]
                   - f_2 * skh1_372[k]
                   + f_3 * pc_y[k] * ski_497[k];

        t_641[k] = f_14 * sii_329[k]
                   + f_3 * pc_z[k] * ski_497[k];

        t_642[k] = f_15 * sii_359[k]
                   + f_4 * skh0_374[k]
                   - f_5 * skh1_374[k]
                   + f_3 * pc_y[k] * ski_499[k];
    }

#pragma omp simd aligned(t_643, t_644, t_645, pc_y, sii_360, sii_361, sii_362, skh0_375, \
                         skh0_376, skh0_377, skh1_375, skh1_376, skh1_377, ski_500, ski_501, \
                         ski_502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_643[k] = f_15 * sii_360[k]
                   + f_6 * skh0_375[k]
                   - f_7 * skh1_375[k]
                   + f_3 * pc_y[k] * ski_500[k];

        t_644[k] = f_15 * sii_361[k]
                   + f_8 * skh0_376[k]
                   - f_9 * skh1_376[k]
                   + f_3 * pc_y[k] * ski_501[k];

        t_645[k] = f_15 * sii_362[k]
                   + f_10 * skh0_377[k]
                   - f_11 * skh1_377[k]
                   + f_3 * pc_y[k] * ski_502[k];
    }

#pragma omp simd aligned(t_646, t_647, t_648, pc_x, pc_y, pc_z, sii_335, sii_363, sii_504, \
                         skh0_377, skh0_378, skh1_377, skh1_378, ski_503, \
                         ski_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_646[k] = f_15 * sii_363[k]
                   + f_3 * pc_y[k] * ski_503[k];

        t_647[k] = f_14 * sii_335[k]
                   + f_1 * skh0_377[k]
                   - f_2 * skh1_377[k]
                   + f_3 * pc_z[k] * ski_503[k];

        t_648[k] = f_14 * sii_504[k]
                   + f_1 * skh0_378[k]
                   - f_2 * skh1_378[k]
                   + f_3 * pc_x[k] * ski_504[k];
    }

#pragma omp simd aligned(t_649, t_650, t_651, t_652, pc_x, pc_y, pc_z, sii_336, sii_364, \
                         sii_366, sii_507, skh0_381, skh1_381, ski_504, ski_506, \
                         ski_507 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_649[k] = f_14 * sii_364[k]
                   + f_3 * pc_y[k] * ski_504[k];

        t_650[k] = f_15 * sii_336[k]
                   + f_3 * pc_z[k] * ski_504[k];

        t_651[k] = f_14 * sii_507[k]
                   + f_4 * skh0_381[k]
                   - f_5 * skh1_381[k]
                   + f_3 * pc_x[k] * ski_507[k];

        t_652[k] = f_14 * sii_366[k]
                   + f_3 * pc_y[k] * ski_506[k];
    }

#pragma omp simd aligned(t_653, t_654, t_655, pc_x, pc_z, sii_339, sii_509, sii_510, skh0_383, \
                         skh0_384, skh1_383, skh1_384, ski_507, ski_509, \
                         ski_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_653[k] = f_14 * sii_509[k]
                   + f_4 * skh0_383[k]
                   - f_5 * skh1_383[k]
                   + f_3 * pc_x[k] * ski_509[k];

        t_654[k] = f_14 * sii_510[k]
                   + f_6 * skh0_384[k]
                   - f_7 * skh1_384[k]
                   + f_3 * pc_x[k] * ski_510[k];

        t_655[k] = f_15 * sii_339[k]
                   + f_3 * pc_z[k] * ski_507[k];
    }

#pragma omp simd aligned(t_656, t_657, t_658, pc_x, pc_y, sii_369, sii_513, sii_514, skh0_387, \
                         skh0_388, skh1_387, skh1_388, ski_509, ski_513, \
                         ski_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_656[k] = f_14 * sii_369[k]
                   + f_3 * pc_y[k] * ski_509[k];

        t_657[k] = f_14 * sii_513[k]
                   + f_6 * skh0_387[k]
                   - f_7 * skh1_387[k]
                   + f_3 * pc_x[k] * ski_513[k];

        t_658[k] = f_14 * sii_514[k]
                   + f_8 * skh0_388[k]
                   - f_9 * skh1_388[k]
                   + f_3 * pc_x[k] * ski_514[k];
    }

#pragma omp simd aligned(t_659, t_660, t_661, pc_x, pc_y, pc_z, sii_342, sii_373, sii_516, \
                         skh0_390, skh1_390, ski_510, ski_513, \
                         ski_516 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_659[k] = f_15 * sii_342[k]
                   + f_3 * pc_z[k] * ski_510[k];

        t_660[k] = f_14 * sii_516[k]
                   + f_8 * skh0_390[k]
                   - f_9 * skh1_390[k]
                   + f_3 * pc_x[k] * ski_516[k];

        t_661[k] = f_14 * sii_373[k]
                   + f_3 * pc_y[k] * ski_513[k];
    }

#pragma omp simd aligned(t_662, t_663, t_664, pc_x, pc_z, sii_346, sii_518, sii_519, skh0_392, \
                         skh0_393, skh1_392, skh1_393, ski_514, ski_518, \
                         ski_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_662[k] = f_14 * sii_518[k]
                   + f_8 * skh0_392[k]
                   - f_9 * skh1_392[k]
                   + f_3 * pc_x[k] * ski_518[k];

        t_663[k] = f_14 * sii_519[k]
                   + f_10 * skh0_393[k]
                   - f_11 * skh1_393[k]
                   + f_3 * pc_x[k] * ski_519[k];

        t_664[k] = f_15 * sii_346[k]
                   + f_3 * pc_z[k] * ski_514[k];
    }

#pragma omp simd aligned(t_665, t_666, t_667, pc_x, pc_y, sii_378, sii_521, sii_522, skh0_395, \
                         skh0_396, skh1_395, skh1_396, ski_518, ski_521, \
                         ski_522 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_665[k] = f_14 * sii_521[k]
                   + f_10 * skh0_395[k]
                   - f_11 * skh1_395[k]
                   + f_3 * pc_x[k] * ski_521[k];

        t_666[k] = f_14 * sii_522[k]
                   + f_10 * skh0_396[k]
                   - f_11 * skh1_396[k]
                   + f_3 * pc_x[k] * ski_522[k];

        t_667[k] = f_14 * sii_378[k]
                   + f_3 * pc_y[k] * ski_518[k];
    }

#pragma omp simd aligned(t_668, t_669, t_670, t_671, pc_x, sii_524, sii_525, sii_526, sii_527, \
                         skh0_398, skh1_398, ski_524, ski_525, ski_526, \
                         ski_527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_668[k] = f_14 * sii_524[k]
                   + f_10 * skh0_398[k]
                   - f_11 * skh1_398[k]
                   + f_3 * pc_x[k] * ski_524[k];

        t_669[k] = f_14 * sii_525[k]
                   + f_3 * pc_x[k] * ski_525[k];

        t_670[k] = f_14 * sii_526[k]
                   + f_3 * pc_x[k] * ski_526[k];

        t_671[k] = f_14 * sii_527[k]
                   + f_3 * pc_x[k] * ski_527[k];
    }

#pragma omp simd aligned(t_672, t_673, t_674, t_675, pc_x, sii_528, sii_529, sii_530, sii_531, \
                         ski_528, ski_529, ski_530, ski_531 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_672[k] = f_14 * sii_528[k]
                   + f_3 * pc_x[k] * ski_528[k];

        t_673[k] = f_14 * sii_529[k]
                   + f_3 * pc_x[k] * ski_529[k];

        t_674[k] = f_14 * sii_530[k]
                   + f_3 * pc_x[k] * ski_530[k];

        t_675[k] = f_14 * sii_531[k]
                   + f_3 * pc_x[k] * ski_531[k];
    }

#pragma omp simd aligned(t_676, t_677, t_678, pc_y, pc_z, sii_357, sii_385, sii_387, skh0_393, \
                         skh0_395, skh1_393, skh1_395, ski_525, \
                         ski_527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_676[k] = f_14 * sii_385[k]
                   + f_1 * skh0_393[k]
                   - f_2 * skh1_393[k]
                   + f_3 * pc_y[k] * ski_525[k];

        t_677[k] = f_15 * sii_357[k]
                   + f_3 * pc_z[k] * ski_525[k];

        t_678[k] = f_14 * sii_387[k]
                   + f_4 * skh0_395[k]
                   - f_5 * skh1_395[k]
                   + f_3 * pc_y[k] * ski_527[k];
    }
}

static auto
compute_prim_skk_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sik0,
                                                          const size_t sii, const size_t sik1,
                                                          const size_t skh0, const size_t skh1,
                                                          const size_t ski, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / q;
    const auto f_1 = 3.0 / gamma;
    const auto f_2 = 3.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.0 / gamma;
    const auto f_5 = 2.0 * p / (gamma * q);
    const auto f_6 = 1.5 / gamma;
    const auto f_7 = 1.5 * p / (gamma * q);
    const auto f_8 = 1.0 / gamma;
    const auto f_9 = p / (gamma * q);
    const auto f_10 = 0.5 / gamma;
    const auto f_11 = 0.5 * p / (gamma * q);
    const auto f_12 = gamma / q;
    const auto f_13 = 0.5 / q;
    const auto f_14 = 1.0 / q;
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;
    const auto f_18 = 3.0 / q;

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
    auto *t_791 = buffer.data(target + 791);
    auto *t_792 = buffer.data(target + 792);
    auto *t_793 = buffer.data(target + 793);
    auto *t_794 = buffer.data(target + 794);
    auto *t_795 = buffer.data(target + 795);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sik0_504 = buffer.data(sik0 + 504);
    const auto *sik0_507 = buffer.data(sik0 + 507);
    const auto *sik0_509 = buffer.data(sik0 + 509);
    const auto *sik0_510 = buffer.data(sik0 + 510);
    const auto *sik0_513 = buffer.data(sik0 + 513);
    const auto *sik0_514 = buffer.data(sik0 + 514);
    const auto *sik0_516 = buffer.data(sik0 + 516);
    const auto *sik0_518 = buffer.data(sik0 + 518);
    const auto *sik0_519 = buffer.data(sik0 + 519);
    const auto *sik0_521 = buffer.data(sik0 + 521);
    const auto *sik0_522 = buffer.data(sik0 + 522);
    const auto *sik0_524 = buffer.data(sik0 + 524);
    const auto *sik0_539 = buffer.data(sik0 + 539);
    const auto *sik0_540 = buffer.data(sik0 + 540);
    const auto *sik0_543 = buffer.data(sik0 + 543);
    const auto *sik0_756 = buffer.data(sik0 + 756);
    const auto *sik0_759 = buffer.data(sik0 + 759);
    const auto *sik0_761 = buffer.data(sik0 + 761);
    const auto *sik0_762 = buffer.data(sik0 + 762);
    const auto *sik0_765 = buffer.data(sik0 + 765);
    const auto *sik0_766 = buffer.data(sik0 + 766);
    const auto *sik0_768 = buffer.data(sik0 + 768);
    const auto *sik0_770 = buffer.data(sik0 + 770);
    const auto *sik0_771 = buffer.data(sik0 + 771);
    const auto *sik0_773 = buffer.data(sik0 + 773);
    const auto *sik0_774 = buffer.data(sik0 + 774);
    const auto *sik0_776 = buffer.data(sik0 + 776);
    const auto *sik0_784 = buffer.data(sik0 + 784);
    const auto *sik0_786 = buffer.data(sik0 + 786);
    const auto *sik0_787 = buffer.data(sik0 + 787);
    const auto *sik0_788 = buffer.data(sik0 + 788);
    const auto *sik0_789 = buffer.data(sik0 + 789);
    const auto *sik0_791 = buffer.data(sik0 + 791);

    const auto *sii_363 = buffer.data(sii + 363);
    const auto *sii_364 = buffer.data(sii + 364);
    const auto *sii_367 = buffer.data(sii + 367);
    const auto *sii_370 = buffer.data(sii + 370);
    const auto *sii_374 = buffer.data(sii + 374);
    const auto *sii_385 = buffer.data(sii + 385);
    const auto *sii_388 = buffer.data(sii + 388);
    const auto *sii_389 = buffer.data(sii + 389);
    const auto *sii_390 = buffer.data(sii + 390);
    const auto *sii_391 = buffer.data(sii + 391);
    const auto *sii_392 = buffer.data(sii + 392);
    const auto *sii_393 = buffer.data(sii + 393);
    const auto *sii_394 = buffer.data(sii + 394);
    const auto *sii_395 = buffer.data(sii + 395);
    const auto *sii_397 = buffer.data(sii + 397);
    const auto *sii_398 = buffer.data(sii + 398);
    const auto *sii_400 = buffer.data(sii + 400);
    const auto *sii_401 = buffer.data(sii + 401);
    const auto *sii_402 = buffer.data(sii + 402);
    const auto *sii_404 = buffer.data(sii + 404);
    const auto *sii_405 = buffer.data(sii + 405);
    const auto *sii_406 = buffer.data(sii + 406);
    const auto *sii_413 = buffer.data(sii + 413);
    const auto *sii_415 = buffer.data(sii + 415);
    const auto *sii_416 = buffer.data(sii + 416);
    const auto *sii_417 = buffer.data(sii + 417);
    const auto *sii_418 = buffer.data(sii + 418);
    const auto *sii_419 = buffer.data(sii + 419);
    const auto *sii_420 = buffer.data(sii + 420);
    const auto *sii_422 = buffer.data(sii + 422);
    const auto *sii_425 = buffer.data(sii + 425);
    const auto *sii_429 = buffer.data(sii + 429);
    const auto *sii_434 = buffer.data(sii + 434);
    const auto *sii_447 = buffer.data(sii + 447);
    const auto *sii_448 = buffer.data(sii + 448);
    const auto *sii_553 = buffer.data(sii + 553);
    const auto *sii_554 = buffer.data(sii + 554);
    const auto *sii_555 = buffer.data(sii + 555);
    const auto *sii_556 = buffer.data(sii + 556);
    const auto *sii_557 = buffer.data(sii + 557);
    const auto *sii_558 = buffer.data(sii + 558);
    const auto *sii_559 = buffer.data(sii + 559);
    const auto *sii_560 = buffer.data(sii + 560);
    const auto *sii_563 = buffer.data(sii + 563);
    const auto *sii_565 = buffer.data(sii + 565);
    const auto *sii_566 = buffer.data(sii + 566);
    const auto *sii_569 = buffer.data(sii + 569);
    const auto *sii_570 = buffer.data(sii + 570);
    const auto *sii_572 = buffer.data(sii + 572);
    const auto *sii_574 = buffer.data(sii + 574);
    const auto *sii_575 = buffer.data(sii + 575);
    const auto *sii_577 = buffer.data(sii + 577);
    const auto *sii_578 = buffer.data(sii + 578);
    const auto *sii_580 = buffer.data(sii + 580);
    const auto *sii_581 = buffer.data(sii + 581);
    const auto *sii_582 = buffer.data(sii + 582);
    const auto *sii_583 = buffer.data(sii + 583);
    const auto *sii_584 = buffer.data(sii + 584);
    const auto *sii_585 = buffer.data(sii + 585);
    const auto *sii_586 = buffer.data(sii + 586);
    const auto *sii_587 = buffer.data(sii + 587);
    const auto *sii_588 = buffer.data(sii + 588);
    const auto *sii_591 = buffer.data(sii + 591);
    const auto *sii_593 = buffer.data(sii + 593);
    const auto *sii_594 = buffer.data(sii + 594);
    const auto *sii_597 = buffer.data(sii + 597);
    const auto *sii_598 = buffer.data(sii + 598);
    const auto *sii_600 = buffer.data(sii + 600);
    const auto *sii_602 = buffer.data(sii + 602);
    const auto *sii_603 = buffer.data(sii + 603);
    const auto *sii_605 = buffer.data(sii + 605);
    const auto *sii_606 = buffer.data(sii + 606);
    const auto *sii_608 = buffer.data(sii + 608);
    const auto *sii_609 = buffer.data(sii + 609);
    const auto *sii_610 = buffer.data(sii + 610);
    const auto *sii_611 = buffer.data(sii + 611);
    const auto *sii_612 = buffer.data(sii + 612);
    const auto *sii_613 = buffer.data(sii + 613);
    const auto *sii_614 = buffer.data(sii + 614);
    const auto *sii_615 = buffer.data(sii + 615);

    const auto *sik1_504 = buffer.data(sik1 + 504);
    const auto *sik1_507 = buffer.data(sik1 + 507);
    const auto *sik1_509 = buffer.data(sik1 + 509);
    const auto *sik1_510 = buffer.data(sik1 + 510);
    const auto *sik1_513 = buffer.data(sik1 + 513);
    const auto *sik1_514 = buffer.data(sik1 + 514);
    const auto *sik1_516 = buffer.data(sik1 + 516);
    const auto *sik1_518 = buffer.data(sik1 + 518);
    const auto *sik1_519 = buffer.data(sik1 + 519);
    const auto *sik1_521 = buffer.data(sik1 + 521);
    const auto *sik1_522 = buffer.data(sik1 + 522);
    const auto *sik1_524 = buffer.data(sik1 + 524);
    const auto *sik1_539 = buffer.data(sik1 + 539);
    const auto *sik1_540 = buffer.data(sik1 + 540);
    const auto *sik1_543 = buffer.data(sik1 + 543);
    const auto *sik1_756 = buffer.data(sik1 + 756);
    const auto *sik1_759 = buffer.data(sik1 + 759);
    const auto *sik1_761 = buffer.data(sik1 + 761);
    const auto *sik1_762 = buffer.data(sik1 + 762);
    const auto *sik1_765 = buffer.data(sik1 + 765);
    const auto *sik1_766 = buffer.data(sik1 + 766);
    const auto *sik1_768 = buffer.data(sik1 + 768);
    const auto *sik1_770 = buffer.data(sik1 + 770);
    const auto *sik1_771 = buffer.data(sik1 + 771);
    const auto *sik1_773 = buffer.data(sik1 + 773);
    const auto *sik1_774 = buffer.data(sik1 + 774);
    const auto *sik1_776 = buffer.data(sik1 + 776);
    const auto *sik1_784 = buffer.data(sik1 + 784);
    const auto *sik1_786 = buffer.data(sik1 + 786);
    const auto *sik1_787 = buffer.data(sik1 + 787);
    const auto *sik1_788 = buffer.data(sik1 + 788);
    const auto *sik1_789 = buffer.data(sik1 + 789);
    const auto *sik1_791 = buffer.data(sik1 + 791);

    const auto *skh0_396 = buffer.data(skh0 + 396);
    const auto *skh0_397 = buffer.data(skh0 + 397);
    const auto *skh0_398 = buffer.data(skh0 + 398);
    const auto *skh0_414 = buffer.data(skh0 + 414);
    const auto *skh0_416 = buffer.data(skh0 + 416);
    const auto *skh0_417 = buffer.data(skh0 + 417);
    const auto *skh0_418 = buffer.data(skh0 + 418);
    const auto *skh0_419 = buffer.data(skh0 + 419);
    const auto *skh0_420 = buffer.data(skh0 + 420);
    const auto *skh0_423 = buffer.data(skh0 + 423);
    const auto *skh0_425 = buffer.data(skh0 + 425);
    const auto *skh0_426 = buffer.data(skh0 + 426);
    const auto *skh0_429 = buffer.data(skh0 + 429);
    const auto *skh0_430 = buffer.data(skh0 + 430);
    const auto *skh0_432 = buffer.data(skh0 + 432);
    const auto *skh0_434 = buffer.data(skh0 + 434);
    const auto *skh0_435 = buffer.data(skh0 + 435);
    const auto *skh0_437 = buffer.data(skh0 + 437);
    const auto *skh0_438 = buffer.data(skh0 + 438);
    const auto *skh0_439 = buffer.data(skh0 + 439);
    const auto *skh0_440 = buffer.data(skh0 + 440);

    const auto *skh1_396 = buffer.data(skh1 + 396);
    const auto *skh1_397 = buffer.data(skh1 + 397);
    const auto *skh1_398 = buffer.data(skh1 + 398);
    const auto *skh1_414 = buffer.data(skh1 + 414);
    const auto *skh1_416 = buffer.data(skh1 + 416);
    const auto *skh1_417 = buffer.data(skh1 + 417);
    const auto *skh1_418 = buffer.data(skh1 + 418);
    const auto *skh1_419 = buffer.data(skh1 + 419);
    const auto *skh1_420 = buffer.data(skh1 + 420);
    const auto *skh1_423 = buffer.data(skh1 + 423);
    const auto *skh1_425 = buffer.data(skh1 + 425);
    const auto *skh1_426 = buffer.data(skh1 + 426);
    const auto *skh1_429 = buffer.data(skh1 + 429);
    const auto *skh1_430 = buffer.data(skh1 + 430);
    const auto *skh1_432 = buffer.data(skh1 + 432);
    const auto *skh1_434 = buffer.data(skh1 + 434);
    const auto *skh1_435 = buffer.data(skh1 + 435);
    const auto *skh1_437 = buffer.data(skh1 + 437);
    const auto *skh1_438 = buffer.data(skh1 + 438);
    const auto *skh1_439 = buffer.data(skh1 + 439);
    const auto *skh1_440 = buffer.data(skh1 + 440);

    const auto *ski_528 = buffer.data(ski + 528);
    const auto *ski_529 = buffer.data(ski + 529);
    const auto *ski_530 = buffer.data(ski + 530);
    const auto *ski_531 = buffer.data(ski + 531);
    const auto *ski_532 = buffer.data(ski + 532);
    const auto *ski_534 = buffer.data(ski + 534);
    const auto *ski_535 = buffer.data(ski + 535);
    const auto *ski_537 = buffer.data(ski + 537);
    const auto *ski_538 = buffer.data(ski + 538);
    const auto *ski_541 = buffer.data(ski + 541);
    const auto *ski_542 = buffer.data(ski + 542);
    const auto *ski_546 = buffer.data(ski + 546);
    const auto *ski_553 = buffer.data(ski + 553);
    const auto *ski_554 = buffer.data(ski + 554);
    const auto *ski_555 = buffer.data(ski + 555);
    const auto *ski_556 = buffer.data(ski + 556);
    const auto *ski_557 = buffer.data(ski + 557);
    const auto *ski_558 = buffer.data(ski + 558);
    const auto *ski_559 = buffer.data(ski + 559);
    const auto *ski_560 = buffer.data(ski + 560);
    const auto *ski_562 = buffer.data(ski + 562);
    const auto *ski_563 = buffer.data(ski + 563);
    const auto *ski_565 = buffer.data(ski + 565);
    const auto *ski_566 = buffer.data(ski + 566);
    const auto *ski_569 = buffer.data(ski + 569);
    const auto *ski_570 = buffer.data(ski + 570);
    const auto *ski_572 = buffer.data(ski + 572);
    const auto *ski_574 = buffer.data(ski + 574);
    const auto *ski_575 = buffer.data(ski + 575);
    const auto *ski_577 = buffer.data(ski + 577);
    const auto *ski_578 = buffer.data(ski + 578);
    const auto *ski_580 = buffer.data(ski + 580);
    const auto *ski_581 = buffer.data(ski + 581);
    const auto *ski_582 = buffer.data(ski + 582);
    const auto *ski_583 = buffer.data(ski + 583);
    const auto *ski_584 = buffer.data(ski + 584);
    const auto *ski_585 = buffer.data(ski + 585);
    const auto *ski_586 = buffer.data(ski + 586);
    const auto *ski_587 = buffer.data(ski + 587);
    const auto *ski_588 = buffer.data(ski + 588);
    const auto *ski_590 = buffer.data(ski + 590);
    const auto *ski_591 = buffer.data(ski + 591);
    const auto *ski_593 = buffer.data(ski + 593);
    const auto *ski_594 = buffer.data(ski + 594);
    const auto *ski_597 = buffer.data(ski + 597);
    const auto *ski_598 = buffer.data(ski + 598);
    const auto *ski_602 = buffer.data(ski + 602);
    const auto *ski_609 = buffer.data(ski + 609);
    const auto *ski_610 = buffer.data(ski + 610);
    const auto *ski_611 = buffer.data(ski + 611);
    const auto *ski_612 = buffer.data(ski + 612);
    const auto *ski_613 = buffer.data(ski + 613);
    const auto *ski_614 = buffer.data(ski + 614);
    const auto *ski_615 = buffer.data(ski + 615);
    const auto *ski_616 = buffer.data(ski + 616);

#pragma omp simd aligned(t_679, t_680, t_681, pc_y, sii_388, sii_389, sii_390, skh0_396, \
                         skh0_397, skh0_398, skh1_396, skh1_397, skh1_398, ski_528, ski_529, \
                         ski_530 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_679[k] = f_14 * sii_388[k]
                   + f_6 * skh0_396[k]
                   - f_7 * skh1_396[k]
                   + f_3 * pc_y[k] * ski_528[k];

        t_680[k] = f_14 * sii_389[k]
                   + f_8 * skh0_397[k]
                   - f_9 * skh1_397[k]
                   + f_3 * pc_y[k] * ski_529[k];

        t_681[k] = f_14 * sii_390[k]
                   + f_10 * skh0_398[k]
                   - f_11 * skh1_398[k]
                   + f_3 * pc_y[k] * ski_530[k];
    }

#pragma omp simd aligned(t_682, t_683, t_684, t_685, pb_y, pc_y, pc_z, sik0_504, sii_363, \
                         sii_391, sii_392, sik1_504, skh0_398, skh1_398, ski_531, \
                         ski_532 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_682[k] = f_14 * sii_391[k]
                   + f_3 * pc_y[k] * ski_531[k];

        t_683[k] = f_15 * sii_363[k]
                   + f_1 * skh0_398[k]
                   - f_2 * skh1_398[k]
                   + f_3 * pc_z[k] * ski_531[k];

        t_684[k] = pb_y[k] * sik0_504[k]
                   - f_12 * pc_y[k] * sik1_504[k];

        t_685[k] = f_13 * sii_392[k]
                   + f_3 * pc_y[k] * ski_532[k];
    }

#pragma omp simd aligned(t_686, t_687, t_688, t_689, pb_y, pc_y, pc_z, sik0_507, sik0_509, \
                         sii_364, sii_393, sii_394, sik1_507, sik1_509, ski_532, \
                         ski_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_686[k] = f_16 * sii_364[k]
                   + f_3 * pc_z[k] * ski_532[k];

        t_687[k] = pb_y[k] * sik0_507[k]
                   + f_14 * sii_393[k]
                   - f_12 * pc_y[k] * sik1_507[k];

        t_688[k] = f_13 * sii_394[k]
                   + f_3 * pc_y[k] * ski_534[k];

        t_689[k] = pb_y[k] * sik0_509[k]
                   - f_12 * pc_y[k] * sik1_509[k];
    }

#pragma omp simd aligned(t_690, t_691, t_692, t_693, pb_y, pc_y, pc_z, sik0_510, sik0_513, \
                         sii_367, sii_395, sii_397, sik1_510, sik1_513, ski_535, \
                         ski_537 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_690[k] = pb_y[k] * sik0_510[k]
                   + f_15 * sii_395[k]
                   - f_12 * pc_y[k] * sik1_510[k];

        t_691[k] = f_16 * sii_367[k]
                   + f_3 * pc_z[k] * ski_535[k];

        t_692[k] = f_13 * sii_397[k]
                   + f_3 * pc_y[k] * ski_537[k];

        t_693[k] = pb_y[k] * sik0_513[k]
                   - f_12 * pc_y[k] * sik1_513[k];
    }

#pragma omp simd aligned(t_694, t_695, t_696, pb_y, pc_y, pc_z, sik0_514, sik0_516, sii_370, \
                         sii_398, sii_400, sik1_514, sik1_516, \
                         ski_538 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_694[k] = pb_y[k] * sik0_514[k]
                   + f_16 * sii_398[k]
                   - f_12 * pc_y[k] * sik1_514[k];

        t_695[k] = f_16 * sii_370[k]
                   + f_3 * pc_z[k] * ski_538[k];

        t_696[k] = pb_y[k] * sik0_516[k]
                   + f_14 * sii_400[k]
                   - f_12 * pc_y[k] * sik1_516[k];
    }

#pragma omp simd aligned(t_697, t_698, t_699, t_700, pb_y, pc_y, pc_z, sik0_518, sik0_519, \
                         sii_374, sii_401, sii_402, sik1_518, sik1_519, ski_541, \
                         ski_542 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_697[k] = f_13 * sii_401[k]
                   + f_3 * pc_y[k] * ski_541[k];

        t_698[k] = pb_y[k] * sik0_518[k]
                   - f_12 * pc_y[k] * sik1_518[k];

        t_699[k] = pb_y[k] * sik0_519[k]
                   + f_17 * sii_402[k]
                   - f_12 * pc_y[k] * sik1_519[k];

        t_700[k] = f_16 * sii_374[k]
                   + f_3 * pc_z[k] * ski_542[k];
    }

#pragma omp simd aligned(t_701, t_702, t_703, t_704, pb_y, pc_y, sik0_521, sik0_522, sik0_524, \
                         sii_404, sii_405, sii_406, sik1_521, sik1_522, sik1_524, \
                         ski_546 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_701[k] = pb_y[k] * sik0_521[k]
                   + f_15 * sii_404[k]
                   - f_12 * pc_y[k] * sik1_521[k];

        t_702[k] = pb_y[k] * sik0_522[k]
                   + f_14 * sii_405[k]
                   - f_12 * pc_y[k] * sik1_522[k];

        t_703[k] = f_13 * sii_406[k]
                   + f_3 * pc_y[k] * ski_546[k];

        t_704[k] = pb_y[k] * sik0_524[k]
                   - f_12 * pc_y[k] * sik1_524[k];
    }

#pragma omp simd aligned(t_705, t_706, t_707, t_708, t_709, pc_x, sii_553, sii_554, sii_555, \
                         sii_556, sii_557, ski_553, ski_554, ski_555, ski_556, \
                         ski_557 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_705[k] = f_14 * sii_553[k]
                   + f_3 * pc_x[k] * ski_553[k];

        t_706[k] = f_14 * sii_554[k]
                   + f_3 * pc_x[k] * ski_554[k];

        t_707[k] = f_14 * sii_555[k]
                   + f_3 * pc_x[k] * ski_555[k];

        t_708[k] = f_14 * sii_556[k]
                   + f_3 * pc_x[k] * ski_556[k];

        t_709[k] = f_14 * sii_557[k]
                   + f_3 * pc_x[k] * ski_557[k];
    }

#pragma omp simd aligned(t_710, t_711, t_712, t_713, pc_x, pc_y, pc_z, sii_385, sii_413, \
                         sii_558, sii_559, skh0_414, skh1_414, ski_553, ski_558, \
                         ski_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_710[k] = f_14 * sii_558[k]
                   + f_3 * pc_x[k] * ski_558[k];

        t_711[k] = f_14 * sii_559[k]
                   + f_3 * pc_x[k] * ski_559[k];

        t_712[k] = f_13 * sii_413[k]
                   + f_1 * skh0_414[k]
                   - f_2 * skh1_414[k]
                   + f_3 * pc_y[k] * ski_553[k];

        t_713[k] = f_16 * sii_385[k]
                   + f_3 * pc_z[k] * ski_553[k];
    }

#pragma omp simd aligned(t_714, t_715, t_716, pc_y, sii_415, sii_416, sii_417, skh0_416, \
                         skh0_417, skh0_418, skh1_416, skh1_417, skh1_418, ski_555, ski_556, \
                         ski_557 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_714[k] = f_13 * sii_415[k]
                   + f_4 * skh0_416[k]
                   - f_5 * skh1_416[k]
                   + f_3 * pc_y[k] * ski_555[k];

        t_715[k] = f_13 * sii_416[k]
                   + f_6 * skh0_417[k]
                   - f_7 * skh1_417[k]
                   + f_3 * pc_y[k] * ski_556[k];

        t_716[k] = f_13 * sii_417[k]
                   + f_8 * skh0_418[k]
                   - f_9 * skh1_418[k]
                   + f_3 * pc_y[k] * ski_557[k];
    }

#pragma omp simd aligned(t_717, t_718, t_719, pb_y, pc_y, sik0_539, sii_418, sii_419, \
                         sik1_539, skh0_419, skh1_419, ski_558, \
                         ski_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_717[k] = f_13 * sii_418[k]
                   + f_10 * skh0_419[k]
                   - f_11 * skh1_419[k]
                   + f_3 * pc_y[k] * ski_558[k];

        t_718[k] = f_13 * sii_419[k]
                   + f_3 * pc_y[k] * ski_559[k];

        t_719[k] = pb_y[k] * sik0_539[k]
                   - f_12 * pc_y[k] * sik1_539[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, t_723, pc_x, pc_y, pc_z, sii_392, sii_560, \
                         sii_563, skh0_420, skh0_423, skh1_420, skh1_423, ski_560, \
                         ski_563 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = f_14 * sii_560[k]
                   + f_1 * skh0_420[k]
                   - f_2 * skh1_420[k]
                   + f_3 * pc_x[k] * ski_560[k];

        t_721[k] = f_3 * pc_y[k] * ski_560[k];

        t_722[k] = f_17 * sii_392[k]
                   + f_3 * pc_z[k] * ski_560[k];

        t_723[k] = f_14 * sii_563[k]
                   + f_4 * skh0_423[k]
                   - f_5 * skh1_423[k]
                   + f_3 * pc_x[k] * ski_563[k];
    }

#pragma omp simd aligned(t_724, t_725, t_726, pc_x, pc_y, sii_565, sii_566, skh0_425, \
                         skh0_426, skh1_425, skh1_426, ski_562, ski_565, \
                         ski_566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_724[k] = f_3 * pc_y[k] * ski_562[k];

        t_725[k] = f_14 * sii_565[k]
                   + f_4 * skh0_425[k]
                   - f_5 * skh1_425[k]
                   + f_3 * pc_x[k] * ski_565[k];

        t_726[k] = f_14 * sii_566[k]
                   + f_6 * skh0_426[k]
                   - f_7 * skh1_426[k]
                   + f_3 * pc_x[k] * ski_566[k];
    }

#pragma omp simd aligned(t_727, t_728, t_729, pc_x, pc_y, pc_z, sii_395, sii_569, skh0_429, \
                         skh1_429, ski_563, ski_565, ski_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_727[k] = f_17 * sii_395[k]
                   + f_3 * pc_z[k] * ski_563[k];

        t_728[k] = f_3 * pc_y[k] * ski_565[k];

        t_729[k] = f_14 * sii_569[k]
                   + f_6 * skh0_429[k]
                   - f_7 * skh1_429[k]
                   + f_3 * pc_x[k] * ski_569[k];
    }

#pragma omp simd aligned(t_730, t_731, t_732, pc_x, pc_z, sii_398, sii_570, sii_572, skh0_430, \
                         skh0_432, skh1_430, skh1_432, ski_566, ski_570, \
                         ski_572 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_730[k] = f_14 * sii_570[k]
                   + f_8 * skh0_430[k]
                   - f_9 * skh1_430[k]
                   + f_3 * pc_x[k] * ski_570[k];

        t_731[k] = f_17 * sii_398[k]
                   + f_3 * pc_z[k] * ski_566[k];

        t_732[k] = f_14 * sii_572[k]
                   + f_8 * skh0_432[k]
                   - f_9 * skh1_432[k]
                   + f_3 * pc_x[k] * ski_572[k];
    }

#pragma omp simd aligned(t_733, t_734, t_735, pc_x, pc_y, sii_574, sii_575, skh0_434, \
                         skh0_435, skh1_434, skh1_435, ski_569, ski_574, \
                         ski_575 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_733[k] = f_3 * pc_y[k] * ski_569[k];

        t_734[k] = f_14 * sii_574[k]
                   + f_8 * skh0_434[k]
                   - f_9 * skh1_434[k]
                   + f_3 * pc_x[k] * ski_574[k];

        t_735[k] = f_14 * sii_575[k]
                   + f_10 * skh0_435[k]
                   - f_11 * skh1_435[k]
                   + f_3 * pc_x[k] * ski_575[k];
    }

#pragma omp simd aligned(t_736, t_737, t_738, pc_x, pc_z, sii_402, sii_577, sii_578, skh0_437, \
                         skh0_438, skh1_437, skh1_438, ski_570, ski_577, \
                         ski_578 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_736[k] = f_17 * sii_402[k]
                   + f_3 * pc_z[k] * ski_570[k];

        t_737[k] = f_14 * sii_577[k]
                   + f_10 * skh0_437[k]
                   - f_11 * skh1_437[k]
                   + f_3 * pc_x[k] * ski_577[k];

        t_738[k] = f_14 * sii_578[k]
                   + f_10 * skh0_438[k]
                   - f_11 * skh1_438[k]
                   + f_3 * pc_x[k] * ski_578[k];
    }

#pragma omp simd aligned(t_739, t_740, t_741, t_742, pc_x, pc_y, sii_580, sii_581, sii_582, \
                         skh0_440, skh1_440, ski_574, ski_580, ski_581, \
                         ski_582 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_739[k] = f_3 * pc_y[k] * ski_574[k];

        t_740[k] = f_14 * sii_580[k]
                   + f_10 * skh0_440[k]
                   - f_11 * skh1_440[k]
                   + f_3 * pc_x[k] * ski_580[k];

        t_741[k] = f_14 * sii_581[k]
                   + f_3 * pc_x[k] * ski_581[k];

        t_742[k] = f_14 * sii_582[k]
                   + f_3 * pc_x[k] * ski_582[k];
    }

#pragma omp simd aligned(t_743, t_744, t_745, t_746, t_747, pc_x, sii_583, sii_584, sii_585, \
                         sii_586, sii_587, ski_583, ski_584, ski_585, ski_586, \
                         ski_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_743[k] = f_14 * sii_583[k]
                   + f_3 * pc_x[k] * ski_583[k];

        t_744[k] = f_14 * sii_584[k]
                   + f_3 * pc_x[k] * ski_584[k];

        t_745[k] = f_14 * sii_585[k]
                   + f_3 * pc_x[k] * ski_585[k];

        t_746[k] = f_14 * sii_586[k]
                   + f_3 * pc_x[k] * ski_586[k];

        t_747[k] = f_14 * sii_587[k]
                   + f_3 * pc_x[k] * ski_587[k];
    }

#pragma omp simd aligned(t_748, t_749, t_750, t_751, pc_y, pc_z, sii_413, skh0_435, skh0_437, \
                         skh0_438, skh1_435, skh1_437, skh1_438, ski_581, ski_583, \
                         ski_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_748[k] = f_1 * skh0_435[k]
                   - f_2 * skh1_435[k]
                   + f_3 * pc_y[k] * ski_581[k];

        t_749[k] = f_17 * sii_413[k]
                   + f_3 * pc_z[k] * ski_581[k];

        t_750[k] = f_4 * skh0_437[k]
                   - f_5 * skh1_437[k]
                   + f_3 * pc_y[k] * ski_583[k];

        t_751[k] = f_6 * skh0_438[k]
                   - f_7 * skh1_438[k]
                   + f_3 * pc_y[k] * ski_584[k];
    }

#pragma omp simd aligned(t_752, t_753, t_754, t_755, pc_y, pc_z, sii_419, skh0_439, skh0_440, \
                         skh1_439, skh1_440, ski_585, ski_586, \
                         ski_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_752[k] = f_8 * skh0_439[k]
                   - f_9 * skh1_439[k]
                   + f_3 * pc_y[k] * ski_585[k];

        t_753[k] = f_10 * skh0_440[k]
                   - f_11 * skh1_440[k]
                   + f_3 * pc_y[k] * ski_586[k];

        t_754[k] = f_3 * pc_y[k] * ski_587[k];

        t_755[k] = f_17 * sii_419[k]
                   + f_1 * skh0_440[k]
                   - f_2 * skh1_440[k]
                   + f_3 * pc_z[k] * ski_587[k];
    }

#pragma omp simd aligned(t_756, t_757, t_758, t_759, pb_x, pc_x, pc_y, pc_z, sik0_756, \
                         sik0_759, sii_420, sii_588, sii_591, sik1_756, sik1_759, \
                         ski_588 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_756[k] = pb_x[k] * sik0_756[k]
                   + f_0 * sii_588[k]
                   - f_12 * pc_x[k] * sik1_756[k];

        t_757[k] = f_18 * sii_420[k]
                   + f_3 * pc_y[k] * ski_588[k];

        t_758[k] = f_3 * pc_z[k] * ski_588[k];

        t_759[k] = pb_x[k] * sik0_759[k]
                   + f_17 * sii_591[k]
                   - f_12 * pc_x[k] * sik1_759[k];
    }

#pragma omp simd aligned(t_760, t_761, t_762, pb_x, pc_x, pc_y, sik0_761, sik0_762, sii_422, \
                         sii_593, sii_594, sik1_761, sik1_762, \
                         ski_590 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_760[k] = f_18 * sii_422[k]
                   + f_3 * pc_y[k] * ski_590[k];

        t_761[k] = pb_x[k] * sik0_761[k]
                   + f_17 * sii_593[k]
                   - f_12 * pc_x[k] * sik1_761[k];

        t_762[k] = pb_x[k] * sik0_762[k]
                   + f_16 * sii_594[k]
                   - f_12 * pc_x[k] * sik1_762[k];
    }

#pragma omp simd aligned(t_763, t_764, t_765, pb_x, pc_x, pc_y, pc_z, sik0_765, sii_425, \
                         sii_597, sik1_765, ski_591, ski_593 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_763[k] = f_3 * pc_z[k] * ski_591[k];

        t_764[k] = f_18 * sii_425[k]
                   + f_3 * pc_y[k] * ski_593[k];

        t_765[k] = pb_x[k] * sik0_765[k]
                   + f_16 * sii_597[k]
                   - f_12 * pc_x[k] * sik1_765[k];
    }

#pragma omp simd aligned(t_766, t_767, t_768, pb_x, pc_x, pc_z, sik0_766, sik0_768, sii_598, \
                         sii_600, sik1_766, sik1_768, ski_594 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_766[k] = pb_x[k] * sik0_766[k]
                   + f_15 * sii_598[k]
                   - f_12 * pc_x[k] * sik1_766[k];

        t_767[k] = f_3 * pc_z[k] * ski_594[k];

        t_768[k] = pb_x[k] * sik0_768[k]
                   + f_15 * sii_600[k]
                   - f_12 * pc_x[k] * sik1_768[k];
    }

#pragma omp simd aligned(t_769, t_770, t_771, pb_x, pc_x, pc_y, sik0_770, sik0_771, sii_429, \
                         sii_602, sii_603, sik1_770, sik1_771, \
                         ski_597 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_769[k] = f_18 * sii_429[k]
                   + f_3 * pc_y[k] * ski_597[k];

        t_770[k] = pb_x[k] * sik0_770[k]
                   + f_15 * sii_602[k]
                   - f_12 * pc_x[k] * sik1_770[k];

        t_771[k] = pb_x[k] * sik0_771[k]
                   + f_14 * sii_603[k]
                   - f_12 * pc_x[k] * sik1_771[k];
    }

#pragma omp simd aligned(t_772, t_773, t_774, pb_x, pc_x, pc_z, sik0_773, sik0_774, sii_605, \
                         sii_606, sik1_773, sik1_774, ski_598 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_772[k] = f_3 * pc_z[k] * ski_598[k];

        t_773[k] = pb_x[k] * sik0_773[k]
                   + f_14 * sii_605[k]
                   - f_12 * pc_x[k] * sik1_773[k];

        t_774[k] = pb_x[k] * sik0_774[k]
                   + f_14 * sii_606[k]
                   - f_12 * pc_x[k] * sik1_774[k];
    }

#pragma omp simd aligned(t_775, t_776, t_777, t_778, pb_x, pc_x, pc_y, sik0_776, sii_434, \
                         sii_608, sii_609, sii_610, sik1_776, ski_602, ski_609, \
                         ski_610 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_775[k] = f_18 * sii_434[k]
                   + f_3 * pc_y[k] * ski_602[k];

        t_776[k] = pb_x[k] * sik0_776[k]
                   + f_14 * sii_608[k]
                   - f_12 * pc_x[k] * sik1_776[k];

        t_777[k] = f_13 * sii_609[k]
                   + f_3 * pc_x[k] * ski_609[k];

        t_778[k] = f_13 * sii_610[k]
                   + f_3 * pc_x[k] * ski_610[k];
    }

#pragma omp simd aligned(t_779, t_780, t_781, t_782, t_783, pc_x, sii_611, sii_612, sii_613, \
                         sii_614, sii_615, ski_611, ski_612, ski_613, ski_614, \
                         ski_615 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_779[k] = f_13 * sii_611[k]
                   + f_3 * pc_x[k] * ski_611[k];

        t_780[k] = f_13 * sii_612[k]
                   + f_3 * pc_x[k] * ski_612[k];

        t_781[k] = f_13 * sii_613[k]
                   + f_3 * pc_x[k] * ski_613[k];

        t_782[k] = f_13 * sii_614[k]
                   + f_3 * pc_x[k] * ski_614[k];

        t_783[k] = f_13 * sii_615[k]
                   + f_3 * pc_x[k] * ski_615[k];
    }

#pragma omp simd aligned(t_784, t_785, t_786, t_787, pb_x, pc_x, pc_z, sik0_784, sik0_786, \
                         sik0_787, sik1_784, sik1_786, sik1_787, \
                         ski_609 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_784[k] = pb_x[k] * sik0_784[k]
                   - f_12 * pc_x[k] * sik1_784[k];

        t_785[k] = f_3 * pc_z[k] * ski_609[k];

        t_786[k] = pb_x[k] * sik0_786[k]
                   - f_12 * pc_x[k] * sik1_786[k];

        t_787[k] = pb_x[k] * sik0_787[k]
                   - f_12 * pc_x[k] * sik1_787[k];
    }

#pragma omp simd aligned(t_788, t_789, t_790, t_791, pb_x, pc_x, pc_y, sik0_788, sik0_789, \
                         sik0_791, sii_447, sik1_788, sik1_789, sik1_791, \
                         ski_615 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_788[k] = pb_x[k] * sik0_788[k]
                   - f_12 * pc_x[k] * sik1_788[k];

        t_789[k] = pb_x[k] * sik0_789[k]
                   - f_12 * pc_x[k] * sik1_789[k];

        t_790[k] = f_18 * sii_447[k]
                   + f_3 * pc_y[k] * ski_615[k];

        t_791[k] = pb_x[k] * sik0_791[k]
                   - f_12 * pc_x[k] * sik1_791[k];
    }

#pragma omp simd aligned(t_792, t_793, t_794, t_795, pb_z, pc_y, pc_z, sik0_540, sik0_543, \
                         sii_420, sii_448, sik1_540, sik1_543, \
                         ski_616 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_792[k] = pb_z[k] * sik0_540[k]
                   - f_12 * pc_z[k] * sik1_540[k];

        t_793[k] = f_17 * sii_448[k]
                   + f_3 * pc_y[k] * ski_616[k];

        t_794[k] = f_13 * sii_420[k]
                   + f_3 * pc_z[k] * ski_616[k];

        t_795[k] = pb_z[k] * sik0_543[k]
                   - f_12 * pc_z[k] * sik1_543[k];
    }
}

static auto
compute_prim_skk_three_center_electron_repulsion_0_piece7(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sik0,
                                                          const size_t sii, const size_t sik1,
                                                          const size_t ski, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / q;
    const auto f_3 = p / q;
    const auto f_12 = gamma / q;
    const auto f_13 = 0.5 / q;
    const auto f_14 = 1.0 / q;
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;

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
    auto *t_911 = buffer.data(target + 911);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sik0_546 = buffer.data(sik0 + 546);
    const auto *sik0_550 = buffer.data(sik0 + 550);
    const auto *sik0_555 = buffer.data(sik0 + 555);
    const auto *sik0_797 = buffer.data(sik0 + 797);
    const auto *sik0_801 = buffer.data(sik0 + 801);
    const auto *sik0_804 = buffer.data(sik0 + 804);
    const auto *sik0_806 = buffer.data(sik0 + 806);
    const auto *sik0_809 = buffer.data(sik0 + 809);
    const auto *sik0_810 = buffer.data(sik0 + 810);
    const auto *sik0_812 = buffer.data(sik0 + 812);
    const auto *sik0_820 = buffer.data(sik0 + 820);
    const auto *sik0_822 = buffer.data(sik0 + 822);
    const auto *sik0_823 = buffer.data(sik0 + 823);
    const auto *sik0_824 = buffer.data(sik0 + 824);
    const auto *sik0_825 = buffer.data(sik0 + 825);
    const auto *sik0_827 = buffer.data(sik0 + 827);
    const auto *sik0_828 = buffer.data(sik0 + 828);
    const auto *sik0_831 = buffer.data(sik0 + 831);
    const auto *sik0_833 = buffer.data(sik0 + 833);
    const auto *sik0_834 = buffer.data(sik0 + 834);
    const auto *sik0_837 = buffer.data(sik0 + 837);
    const auto *sik0_838 = buffer.data(sik0 + 838);
    const auto *sik0_840 = buffer.data(sik0 + 840);
    const auto *sik0_842 = buffer.data(sik0 + 842);
    const auto *sik0_843 = buffer.data(sik0 + 843);
    const auto *sik0_845 = buffer.data(sik0 + 845);
    const auto *sik0_846 = buffer.data(sik0 + 846);
    const auto *sik0_848 = buffer.data(sik0 + 848);
    const auto *sik0_856 = buffer.data(sik0 + 856);
    const auto *sik0_858 = buffer.data(sik0 + 858);
    const auto *sik0_859 = buffer.data(sik0 + 859);
    const auto *sik0_860 = buffer.data(sik0 + 860);
    const auto *sik0_861 = buffer.data(sik0 + 861);
    const auto *sik0_863 = buffer.data(sik0 + 863);
    const auto *sik0_864 = buffer.data(sik0 + 864);
    const auto *sik0_867 = buffer.data(sik0 + 867);
    const auto *sik0_869 = buffer.data(sik0 + 869);
    const auto *sik0_870 = buffer.data(sik0 + 870);
    const auto *sik0_873 = buffer.data(sik0 + 873);
    const auto *sik0_874 = buffer.data(sik0 + 874);
    const auto *sik0_876 = buffer.data(sik0 + 876);
    const auto *sik0_878 = buffer.data(sik0 + 878);
    const auto *sik0_879 = buffer.data(sik0 + 879);
    const auto *sik0_881 = buffer.data(sik0 + 881);
    const auto *sik0_882 = buffer.data(sik0 + 882);
    const auto *sik0_884 = buffer.data(sik0 + 884);
    const auto *sik0_892 = buffer.data(sik0 + 892);
    const auto *sik0_894 = buffer.data(sik0 + 894);
    const auto *sik0_895 = buffer.data(sik0 + 895);
    const auto *sik0_896 = buffer.data(sik0 + 896);
    const auto *sik0_897 = buffer.data(sik0 + 897);
    const auto *sik0_899 = buffer.data(sik0 + 899);
    const auto *sik0_900 = buffer.data(sik0 + 900);
    const auto *sik0_903 = buffer.data(sik0 + 903);
    const auto *sik0_905 = buffer.data(sik0 + 905);
    const auto *sik0_906 = buffer.data(sik0 + 906);
    const auto *sik0_909 = buffer.data(sik0 + 909);
    const auto *sik0_910 = buffer.data(sik0 + 910);

    const auto *sii_423 = buffer.data(sii + 423);
    const auto *sii_426 = buffer.data(sii + 426);
    const auto *sii_430 = buffer.data(sii + 430);
    const auto *sii_441 = buffer.data(sii + 441);
    const auto *sii_448 = buffer.data(sii + 448);
    const auto *sii_450 = buffer.data(sii + 450);
    const auto *sii_451 = buffer.data(sii + 451);
    const auto *sii_453 = buffer.data(sii + 453);
    const auto *sii_454 = buffer.data(sii + 454);
    const auto *sii_457 = buffer.data(sii + 457);
    const auto *sii_458 = buffer.data(sii + 458);
    const auto *sii_462 = buffer.data(sii + 462);
    const auto *sii_469 = buffer.data(sii + 469);
    const auto *sii_475 = buffer.data(sii + 475);
    const auto *sii_476 = buffer.data(sii + 476);
    const auto *sii_478 = buffer.data(sii + 478);
    const auto *sii_479 = buffer.data(sii + 479);
    const auto *sii_481 = buffer.data(sii + 481);
    const auto *sii_482 = buffer.data(sii + 482);
    const auto *sii_485 = buffer.data(sii + 485);
    const auto *sii_486 = buffer.data(sii + 486);
    const auto *sii_490 = buffer.data(sii + 490);
    const auto *sii_497 = buffer.data(sii + 497);
    const auto *sii_503 = buffer.data(sii + 503);
    const auto *sii_504 = buffer.data(sii + 504);
    const auto *sii_506 = buffer.data(sii + 506);
    const auto *sii_507 = buffer.data(sii + 507);
    const auto *sii_509 = buffer.data(sii + 509);
    const auto *sii_510 = buffer.data(sii + 510);
    const auto *sii_513 = buffer.data(sii + 513);
    const auto *sii_518 = buffer.data(sii + 518);
    const auto *sii_531 = buffer.data(sii + 531);
    const auto *sii_532 = buffer.data(sii + 532);
    const auto *sii_534 = buffer.data(sii + 534);
    const auto *sii_537 = buffer.data(sii + 537);
    const auto *sii_621 = buffer.data(sii + 621);
    const auto *sii_625 = buffer.data(sii + 625);
    const auto *sii_628 = buffer.data(sii + 628);
    const auto *sii_630 = buffer.data(sii + 630);
    const auto *sii_633 = buffer.data(sii + 633);
    const auto *sii_634 = buffer.data(sii + 634);
    const auto *sii_636 = buffer.data(sii + 636);
    const auto *sii_637 = buffer.data(sii + 637);
    const auto *sii_638 = buffer.data(sii + 638);
    const auto *sii_639 = buffer.data(sii + 639);
    const auto *sii_640 = buffer.data(sii + 640);
    const auto *sii_641 = buffer.data(sii + 641);
    const auto *sii_642 = buffer.data(sii + 642);
    const auto *sii_643 = buffer.data(sii + 643);
    const auto *sii_644 = buffer.data(sii + 644);
    const auto *sii_647 = buffer.data(sii + 647);
    const auto *sii_649 = buffer.data(sii + 649);
    const auto *sii_650 = buffer.data(sii + 650);
    const auto *sii_653 = buffer.data(sii + 653);
    const auto *sii_654 = buffer.data(sii + 654);
    const auto *sii_656 = buffer.data(sii + 656);
    const auto *sii_658 = buffer.data(sii + 658);
    const auto *sii_659 = buffer.data(sii + 659);
    const auto *sii_661 = buffer.data(sii + 661);
    const auto *sii_662 = buffer.data(sii + 662);
    const auto *sii_664 = buffer.data(sii + 664);
    const auto *sii_665 = buffer.data(sii + 665);
    const auto *sii_666 = buffer.data(sii + 666);
    const auto *sii_667 = buffer.data(sii + 667);
    const auto *sii_668 = buffer.data(sii + 668);
    const auto *sii_669 = buffer.data(sii + 669);
    const auto *sii_670 = buffer.data(sii + 670);
    const auto *sii_671 = buffer.data(sii + 671);
    const auto *sii_672 = buffer.data(sii + 672);
    const auto *sii_675 = buffer.data(sii + 675);
    const auto *sii_677 = buffer.data(sii + 677);
    const auto *sii_678 = buffer.data(sii + 678);
    const auto *sii_681 = buffer.data(sii + 681);
    const auto *sii_682 = buffer.data(sii + 682);
    const auto *sii_684 = buffer.data(sii + 684);
    const auto *sii_686 = buffer.data(sii + 686);
    const auto *sii_687 = buffer.data(sii + 687);
    const auto *sii_689 = buffer.data(sii + 689);
    const auto *sii_690 = buffer.data(sii + 690);
    const auto *sii_692 = buffer.data(sii + 692);
    const auto *sii_693 = buffer.data(sii + 693);
    const auto *sii_694 = buffer.data(sii + 694);
    const auto *sii_695 = buffer.data(sii + 695);
    const auto *sii_696 = buffer.data(sii + 696);
    const auto *sii_697 = buffer.data(sii + 697);
    const auto *sii_698 = buffer.data(sii + 698);
    const auto *sii_699 = buffer.data(sii + 699);
    const auto *sii_700 = buffer.data(sii + 700);
    const auto *sii_703 = buffer.data(sii + 703);
    const auto *sii_705 = buffer.data(sii + 705);
    const auto *sii_706 = buffer.data(sii + 706);
    const auto *sii_709 = buffer.data(sii + 709);
    const auto *sii_710 = buffer.data(sii + 710);

    const auto *sik1_546 = buffer.data(sik1 + 546);
    const auto *sik1_550 = buffer.data(sik1 + 550);
    const auto *sik1_555 = buffer.data(sik1 + 555);
    const auto *sik1_797 = buffer.data(sik1 + 797);
    const auto *sik1_801 = buffer.data(sik1 + 801);
    const auto *sik1_804 = buffer.data(sik1 + 804);
    const auto *sik1_806 = buffer.data(sik1 + 806);
    const auto *sik1_809 = buffer.data(sik1 + 809);
    const auto *sik1_810 = buffer.data(sik1 + 810);
    const auto *sik1_812 = buffer.data(sik1 + 812);
    const auto *sik1_820 = buffer.data(sik1 + 820);
    const auto *sik1_822 = buffer.data(sik1 + 822);
    const auto *sik1_823 = buffer.data(sik1 + 823);
    const auto *sik1_824 = buffer.data(sik1 + 824);
    const auto *sik1_825 = buffer.data(sik1 + 825);
    const auto *sik1_827 = buffer.data(sik1 + 827);
    const auto *sik1_828 = buffer.data(sik1 + 828);
    const auto *sik1_831 = buffer.data(sik1 + 831);
    const auto *sik1_833 = buffer.data(sik1 + 833);
    const auto *sik1_834 = buffer.data(sik1 + 834);
    const auto *sik1_837 = buffer.data(sik1 + 837);
    const auto *sik1_838 = buffer.data(sik1 + 838);
    const auto *sik1_840 = buffer.data(sik1 + 840);
    const auto *sik1_842 = buffer.data(sik1 + 842);
    const auto *sik1_843 = buffer.data(sik1 + 843);
    const auto *sik1_845 = buffer.data(sik1 + 845);
    const auto *sik1_846 = buffer.data(sik1 + 846);
    const auto *sik1_848 = buffer.data(sik1 + 848);
    const auto *sik1_856 = buffer.data(sik1 + 856);
    const auto *sik1_858 = buffer.data(sik1 + 858);
    const auto *sik1_859 = buffer.data(sik1 + 859);
    const auto *sik1_860 = buffer.data(sik1 + 860);
    const auto *sik1_861 = buffer.data(sik1 + 861);
    const auto *sik1_863 = buffer.data(sik1 + 863);
    const auto *sik1_864 = buffer.data(sik1 + 864);
    const auto *sik1_867 = buffer.data(sik1 + 867);
    const auto *sik1_869 = buffer.data(sik1 + 869);
    const auto *sik1_870 = buffer.data(sik1 + 870);
    const auto *sik1_873 = buffer.data(sik1 + 873);
    const auto *sik1_874 = buffer.data(sik1 + 874);
    const auto *sik1_876 = buffer.data(sik1 + 876);
    const auto *sik1_878 = buffer.data(sik1 + 878);
    const auto *sik1_879 = buffer.data(sik1 + 879);
    const auto *sik1_881 = buffer.data(sik1 + 881);
    const auto *sik1_882 = buffer.data(sik1 + 882);
    const auto *sik1_884 = buffer.data(sik1 + 884);
    const auto *sik1_892 = buffer.data(sik1 + 892);
    const auto *sik1_894 = buffer.data(sik1 + 894);
    const auto *sik1_895 = buffer.data(sik1 + 895);
    const auto *sik1_896 = buffer.data(sik1 + 896);
    const auto *sik1_897 = buffer.data(sik1 + 897);
    const auto *sik1_899 = buffer.data(sik1 + 899);
    const auto *sik1_900 = buffer.data(sik1 + 900);
    const auto *sik1_903 = buffer.data(sik1 + 903);
    const auto *sik1_905 = buffer.data(sik1 + 905);
    const auto *sik1_906 = buffer.data(sik1 + 906);
    const auto *sik1_909 = buffer.data(sik1 + 909);
    const auto *sik1_910 = buffer.data(sik1 + 910);

    const auto *ski_618 = buffer.data(ski + 618);
    const auto *ski_619 = buffer.data(ski + 619);
    const auto *ski_621 = buffer.data(ski + 621);
    const auto *ski_622 = buffer.data(ski + 622);
    const auto *ski_625 = buffer.data(ski + 625);
    const auto *ski_626 = buffer.data(ski + 626);
    const auto *ski_630 = buffer.data(ski + 630);
    const auto *ski_637 = buffer.data(ski + 637);
    const auto *ski_638 = buffer.data(ski + 638);
    const auto *ski_639 = buffer.data(ski + 639);
    const auto *ski_640 = buffer.data(ski + 640);
    const auto *ski_641 = buffer.data(ski + 641);
    const auto *ski_642 = buffer.data(ski + 642);
    const auto *ski_643 = buffer.data(ski + 643);
    const auto *ski_644 = buffer.data(ski + 644);
    const auto *ski_646 = buffer.data(ski + 646);
    const auto *ski_647 = buffer.data(ski + 647);
    const auto *ski_649 = buffer.data(ski + 649);
    const auto *ski_650 = buffer.data(ski + 650);
    const auto *ski_653 = buffer.data(ski + 653);
    const auto *ski_654 = buffer.data(ski + 654);
    const auto *ski_658 = buffer.data(ski + 658);
    const auto *ski_665 = buffer.data(ski + 665);
    const auto *ski_666 = buffer.data(ski + 666);
    const auto *ski_667 = buffer.data(ski + 667);
    const auto *ski_668 = buffer.data(ski + 668);
    const auto *ski_669 = buffer.data(ski + 669);
    const auto *ski_670 = buffer.data(ski + 670);
    const auto *ski_671 = buffer.data(ski + 671);
    const auto *ski_672 = buffer.data(ski + 672);
    const auto *ski_674 = buffer.data(ski + 674);
    const auto *ski_675 = buffer.data(ski + 675);
    const auto *ski_677 = buffer.data(ski + 677);
    const auto *ski_678 = buffer.data(ski + 678);
    const auto *ski_681 = buffer.data(ski + 681);
    const auto *ski_682 = buffer.data(ski + 682);
    const auto *ski_686 = buffer.data(ski + 686);
    const auto *ski_693 = buffer.data(ski + 693);
    const auto *ski_694 = buffer.data(ski + 694);
    const auto *ski_695 = buffer.data(ski + 695);
    const auto *ski_696 = buffer.data(ski + 696);
    const auto *ski_697 = buffer.data(ski + 697);
    const auto *ski_698 = buffer.data(ski + 698);
    const auto *ski_699 = buffer.data(ski + 699);
    const auto *ski_700 = buffer.data(ski + 700);
    const auto *ski_702 = buffer.data(ski + 702);
    const auto *ski_703 = buffer.data(ski + 703);
    const auto *ski_705 = buffer.data(ski + 705);
    const auto *ski_706 = buffer.data(ski + 706);

#pragma omp simd aligned(t_796, t_797, t_798, pb_x, pb_z, pc_x, pc_y, pc_z, sik0_546, \
                         sik0_797, sii_450, sii_621, sik1_546, sik1_797, \
                         ski_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_796[k] = f_17 * sii_450[k]
                   + f_3 * pc_y[k] * ski_618[k];

        t_797[k] = pb_x[k] * sik0_797[k]
                   + f_17 * sii_621[k]
                   - f_12 * pc_x[k] * sik1_797[k];

        t_798[k] = pb_z[k] * sik0_546[k]
                   - f_12 * pc_z[k] * sik1_546[k];
    }

#pragma omp simd aligned(t_799, t_800, t_801, pb_x, pc_x, pc_y, pc_z, sik0_801, sii_423, \
                         sii_453, sii_625, sik1_801, ski_619, ski_621 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_799[k] = f_13 * sii_423[k]
                   + f_3 * pc_z[k] * ski_619[k];

        t_800[k] = f_17 * sii_453[k]
                   + f_3 * pc_y[k] * ski_621[k];

        t_801[k] = pb_x[k] * sik0_801[k]
                   + f_16 * sii_625[k]
                   - f_12 * pc_x[k] * sik1_801[k];
    }

#pragma omp simd aligned(t_802, t_803, t_804, pb_x, pb_z, pc_x, pc_z, sik0_550, sik0_804, \
                         sii_426, sii_628, sik1_550, sik1_804, \
                         ski_622 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_802[k] = pb_z[k] * sik0_550[k]
                   - f_12 * pc_z[k] * sik1_550[k];

        t_803[k] = f_13 * sii_426[k]
                   + f_3 * pc_z[k] * ski_622[k];

        t_804[k] = pb_x[k] * sik0_804[k]
                   + f_15 * sii_628[k]
                   - f_12 * pc_x[k] * sik1_804[k];
    }

#pragma omp simd aligned(t_805, t_806, t_807, pb_x, pb_z, pc_x, pc_y, pc_z, sik0_555, \
                         sik0_806, sii_457, sii_630, sik1_555, sik1_806, \
                         ski_625 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_805[k] = f_17 * sii_457[k]
                   + f_3 * pc_y[k] * ski_625[k];

        t_806[k] = pb_x[k] * sik0_806[k]
                   + f_15 * sii_630[k]
                   - f_12 * pc_x[k] * sik1_806[k];

        t_807[k] = pb_z[k] * sik0_555[k]
                   - f_12 * pc_z[k] * sik1_555[k];
    }

#pragma omp simd aligned(t_808, t_809, t_810, pb_x, pc_x, pc_z, sik0_809, sik0_810, sii_430, \
                         sii_633, sii_634, sik1_809, sik1_810, \
                         ski_626 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_808[k] = f_13 * sii_430[k]
                   + f_3 * pc_z[k] * ski_626[k];

        t_809[k] = pb_x[k] * sik0_809[k]
                   + f_14 * sii_633[k]
                   - f_12 * pc_x[k] * sik1_809[k];

        t_810[k] = pb_x[k] * sik0_810[k]
                   + f_14 * sii_634[k]
                   - f_12 * pc_x[k] * sik1_810[k];
    }

#pragma omp simd aligned(t_811, t_812, t_813, t_814, pb_x, pc_x, pc_y, sik0_812, sii_462, \
                         sii_636, sii_637, sii_638, sik1_812, ski_630, ski_637, \
                         ski_638 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_811[k] = f_17 * sii_462[k]
                   + f_3 * pc_y[k] * ski_630[k];

        t_812[k] = pb_x[k] * sik0_812[k]
                   + f_14 * sii_636[k]
                   - f_12 * pc_x[k] * sik1_812[k];

        t_813[k] = f_13 * sii_637[k]
                   + f_3 * pc_x[k] * ski_637[k];

        t_814[k] = f_13 * sii_638[k]
                   + f_3 * pc_x[k] * ski_638[k];
    }

#pragma omp simd aligned(t_815, t_816, t_817, t_818, t_819, pc_x, sii_639, sii_640, sii_641, \
                         sii_642, sii_643, ski_639, ski_640, ski_641, ski_642, \
                         ski_643 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_815[k] = f_13 * sii_639[k]
                   + f_3 * pc_x[k] * ski_639[k];

        t_816[k] = f_13 * sii_640[k]
                   + f_3 * pc_x[k] * ski_640[k];

        t_817[k] = f_13 * sii_641[k]
                   + f_3 * pc_x[k] * ski_641[k];

        t_818[k] = f_13 * sii_642[k]
                   + f_3 * pc_x[k] * ski_642[k];

        t_819[k] = f_13 * sii_643[k]
                   + f_3 * pc_x[k] * ski_643[k];
    }

#pragma omp simd aligned(t_820, t_821, t_822, t_823, pb_x, pc_x, pc_z, sik0_820, sik0_822, \
                         sik0_823, sii_441, sik1_820, sik1_822, sik1_823, \
                         ski_637 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_820[k] = pb_x[k] * sik0_820[k]
                   - f_12 * pc_x[k] * sik1_820[k];

        t_821[k] = f_13 * sii_441[k]
                   + f_3 * pc_z[k] * ski_637[k];

        t_822[k] = pb_x[k] * sik0_822[k]
                   - f_12 * pc_x[k] * sik1_822[k];

        t_823[k] = pb_x[k] * sik0_823[k]
                   - f_12 * pc_x[k] * sik1_823[k];
    }

#pragma omp simd aligned(t_824, t_825, t_826, t_827, pb_x, pc_x, pc_y, sik0_824, sik0_825, \
                         sik0_827, sii_475, sik1_824, sik1_825, sik1_827, \
                         ski_643 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_824[k] = pb_x[k] * sik0_824[k]
                   - f_12 * pc_x[k] * sik1_824[k];

        t_825[k] = pb_x[k] * sik0_825[k]
                   - f_12 * pc_x[k] * sik1_825[k];

        t_826[k] = f_17 * sii_475[k]
                   + f_3 * pc_y[k] * ski_643[k];

        t_827[k] = pb_x[k] * sik0_827[k]
                   - f_12 * pc_x[k] * sik1_827[k];
    }

#pragma omp simd aligned(t_828, t_829, t_830, pb_x, pc_x, pc_y, pc_z, sik0_828, sii_448, \
                         sii_476, sii_644, sik1_828, ski_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_828[k] = pb_x[k] * sik0_828[k]
                   + f_0 * sii_644[k]
                   - f_12 * pc_x[k] * sik1_828[k];

        t_829[k] = f_16 * sii_476[k]
                   + f_3 * pc_y[k] * ski_644[k];

        t_830[k] = f_14 * sii_448[k]
                   + f_3 * pc_z[k] * ski_644[k];
    }

#pragma omp simd aligned(t_831, t_832, t_833, pb_x, pc_x, pc_y, sik0_831, sik0_833, sii_478, \
                         sii_647, sii_649, sik1_831, sik1_833, \
                         ski_646 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_831[k] = pb_x[k] * sik0_831[k]
                   + f_17 * sii_647[k]
                   - f_12 * pc_x[k] * sik1_831[k];

        t_832[k] = f_16 * sii_478[k]
                   + f_3 * pc_y[k] * ski_646[k];

        t_833[k] = pb_x[k] * sik0_833[k]
                   + f_17 * sii_649[k]
                   - f_12 * pc_x[k] * sik1_833[k];
    }

#pragma omp simd aligned(t_834, t_835, t_836, pb_x, pc_x, pc_y, pc_z, sik0_834, sii_451, \
                         sii_481, sii_650, sik1_834, ski_647, ski_649 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_834[k] = pb_x[k] * sik0_834[k]
                   + f_16 * sii_650[k]
                   - f_12 * pc_x[k] * sik1_834[k];

        t_835[k] = f_14 * sii_451[k]
                   + f_3 * pc_z[k] * ski_647[k];

        t_836[k] = f_16 * sii_481[k]
                   + f_3 * pc_y[k] * ski_649[k];
    }

#pragma omp simd aligned(t_837, t_838, t_839, pb_x, pc_x, pc_z, sik0_837, sik0_838, sii_454, \
                         sii_653, sii_654, sik1_837, sik1_838, \
                         ski_650 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_837[k] = pb_x[k] * sik0_837[k]
                   + f_16 * sii_653[k]
                   - f_12 * pc_x[k] * sik1_837[k];

        t_838[k] = pb_x[k] * sik0_838[k]
                   + f_15 * sii_654[k]
                   - f_12 * pc_x[k] * sik1_838[k];

        t_839[k] = f_14 * sii_454[k]
                   + f_3 * pc_z[k] * ski_650[k];
    }

#pragma omp simd aligned(t_840, t_841, t_842, pb_x, pc_x, pc_y, sik0_840, sik0_842, sii_485, \
                         sii_656, sii_658, sik1_840, sik1_842, \
                         ski_653 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_840[k] = pb_x[k] * sik0_840[k]
                   + f_15 * sii_656[k]
                   - f_12 * pc_x[k] * sik1_840[k];

        t_841[k] = f_16 * sii_485[k]
                   + f_3 * pc_y[k] * ski_653[k];

        t_842[k] = pb_x[k] * sik0_842[k]
                   + f_15 * sii_658[k]
                   - f_12 * pc_x[k] * sik1_842[k];
    }

#pragma omp simd aligned(t_843, t_844, t_845, pb_x, pc_x, pc_z, sik0_843, sik0_845, sii_458, \
                         sii_659, sii_661, sik1_843, sik1_845, \
                         ski_654 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_843[k] = pb_x[k] * sik0_843[k]
                   + f_14 * sii_659[k]
                   - f_12 * pc_x[k] * sik1_843[k];

        t_844[k] = f_14 * sii_458[k]
                   + f_3 * pc_z[k] * ski_654[k];

        t_845[k] = pb_x[k] * sik0_845[k]
                   + f_14 * sii_661[k]
                   - f_12 * pc_x[k] * sik1_845[k];
    }

#pragma omp simd aligned(t_846, t_847, t_848, pb_x, pc_x, pc_y, sik0_846, sik0_848, sii_490, \
                         sii_662, sii_664, sik1_846, sik1_848, \
                         ski_658 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_846[k] = pb_x[k] * sik0_846[k]
                   + f_14 * sii_662[k]
                   - f_12 * pc_x[k] * sik1_846[k];

        t_847[k] = f_16 * sii_490[k]
                   + f_3 * pc_y[k] * ski_658[k];

        t_848[k] = pb_x[k] * sik0_848[k]
                   + f_14 * sii_664[k]
                   - f_12 * pc_x[k] * sik1_848[k];
    }

#pragma omp simd aligned(t_849, t_850, t_851, t_852, t_853, pc_x, sii_665, sii_666, sii_667, \
                         sii_668, sii_669, ski_665, ski_666, ski_667, ski_668, \
                         ski_669 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_849[k] = f_13 * sii_665[k]
                   + f_3 * pc_x[k] * ski_665[k];

        t_850[k] = f_13 * sii_666[k]
                   + f_3 * pc_x[k] * ski_666[k];

        t_851[k] = f_13 * sii_667[k]
                   + f_3 * pc_x[k] * ski_667[k];

        t_852[k] = f_13 * sii_668[k]
                   + f_3 * pc_x[k] * ski_668[k];

        t_853[k] = f_13 * sii_669[k]
                   + f_3 * pc_x[k] * ski_669[k];
    }

#pragma omp simd aligned(t_854, t_855, t_856, t_857, pb_x, pc_x, pc_z, sik0_856, sii_469, \
                         sii_670, sii_671, sik1_856, ski_665, ski_670, \
                         ski_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_854[k] = f_13 * sii_670[k]
                   + f_3 * pc_x[k] * ski_670[k];

        t_855[k] = f_13 * sii_671[k]
                   + f_3 * pc_x[k] * ski_671[k];

        t_856[k] = pb_x[k] * sik0_856[k]
                   - f_12 * pc_x[k] * sik1_856[k];

        t_857[k] = f_14 * sii_469[k]
                   + f_3 * pc_z[k] * ski_665[k];
    }

#pragma omp simd aligned(t_858, t_859, t_860, t_861, pb_x, pc_x, sik0_858, sik0_859, sik0_860, \
                         sik0_861, sik1_858, sik1_859, sik1_860, \
                         sik1_861 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_858[k] = pb_x[k] * sik0_858[k]
                   - f_12 * pc_x[k] * sik1_858[k];

        t_859[k] = pb_x[k] * sik0_859[k]
                   - f_12 * pc_x[k] * sik1_859[k];

        t_860[k] = pb_x[k] * sik0_860[k]
                   - f_12 * pc_x[k] * sik1_860[k];

        t_861[k] = pb_x[k] * sik0_861[k]
                   - f_12 * pc_x[k] * sik1_861[k];
    }

#pragma omp simd aligned(t_862, t_863, t_864, t_865, pb_x, pc_x, pc_y, sik0_863, sik0_864, \
                         sii_503, sii_504, sii_672, sik1_863, sik1_864, ski_671, \
                         ski_672 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_862[k] = f_16 * sii_503[k]
                   + f_3 * pc_y[k] * ski_671[k];

        t_863[k] = pb_x[k] * sik0_863[k]
                   - f_12 * pc_x[k] * sik1_863[k];

        t_864[k] = pb_x[k] * sik0_864[k]
                   + f_0 * sii_672[k]
                   - f_12 * pc_x[k] * sik1_864[k];

        t_865[k] = f_15 * sii_504[k]
                   + f_3 * pc_y[k] * ski_672[k];
    }

#pragma omp simd aligned(t_866, t_867, t_868, pb_x, pc_x, pc_y, pc_z, sik0_867, sii_476, \
                         sii_506, sii_675, sik1_867, ski_672, ski_674 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_866[k] = f_15 * sii_476[k]
                   + f_3 * pc_z[k] * ski_672[k];

        t_867[k] = pb_x[k] * sik0_867[k]
                   + f_17 * sii_675[k]
                   - f_12 * pc_x[k] * sik1_867[k];

        t_868[k] = f_15 * sii_506[k]
                   + f_3 * pc_y[k] * ski_674[k];
    }

#pragma omp simd aligned(t_869, t_870, t_871, pb_x, pc_x, pc_z, sik0_869, sik0_870, sii_479, \
                         sii_677, sii_678, sik1_869, sik1_870, \
                         ski_675 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_869[k] = pb_x[k] * sik0_869[k]
                   + f_17 * sii_677[k]
                   - f_12 * pc_x[k] * sik1_869[k];

        t_870[k] = pb_x[k] * sik0_870[k]
                   + f_16 * sii_678[k]
                   - f_12 * pc_x[k] * sik1_870[k];

        t_871[k] = f_15 * sii_479[k]
                   + f_3 * pc_z[k] * ski_675[k];
    }

#pragma omp simd aligned(t_872, t_873, t_874, pb_x, pc_x, pc_y, sik0_873, sik0_874, sii_509, \
                         sii_681, sii_682, sik1_873, sik1_874, \
                         ski_677 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_872[k] = f_15 * sii_509[k]
                   + f_3 * pc_y[k] * ski_677[k];

        t_873[k] = pb_x[k] * sik0_873[k]
                   + f_16 * sii_681[k]
                   - f_12 * pc_x[k] * sik1_873[k];

        t_874[k] = pb_x[k] * sik0_874[k]
                   + f_15 * sii_682[k]
                   - f_12 * pc_x[k] * sik1_874[k];
    }

#pragma omp simd aligned(t_875, t_876, t_877, pb_x, pc_x, pc_y, pc_z, sik0_876, sii_482, \
                         sii_513, sii_684, sik1_876, ski_678, ski_681 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_875[k] = f_15 * sii_482[k]
                   + f_3 * pc_z[k] * ski_678[k];

        t_876[k] = pb_x[k] * sik0_876[k]
                   + f_15 * sii_684[k]
                   - f_12 * pc_x[k] * sik1_876[k];

        t_877[k] = f_15 * sii_513[k]
                   + f_3 * pc_y[k] * ski_681[k];
    }

#pragma omp simd aligned(t_878, t_879, t_880, pb_x, pc_x, pc_z, sik0_878, sik0_879, sii_486, \
                         sii_686, sii_687, sik1_878, sik1_879, \
                         ski_682 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_878[k] = pb_x[k] * sik0_878[k]
                   + f_15 * sii_686[k]
                   - f_12 * pc_x[k] * sik1_878[k];

        t_879[k] = pb_x[k] * sik0_879[k]
                   + f_14 * sii_687[k]
                   - f_12 * pc_x[k] * sik1_879[k];

        t_880[k] = f_15 * sii_486[k]
                   + f_3 * pc_z[k] * ski_682[k];
    }

#pragma omp simd aligned(t_881, t_882, t_883, pb_x, pc_x, pc_y, sik0_881, sik0_882, sii_518, \
                         sii_689, sii_690, sik1_881, sik1_882, \
                         ski_686 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_881[k] = pb_x[k] * sik0_881[k]
                   + f_14 * sii_689[k]
                   - f_12 * pc_x[k] * sik1_881[k];

        t_882[k] = pb_x[k] * sik0_882[k]
                   + f_14 * sii_690[k]
                   - f_12 * pc_x[k] * sik1_882[k];

        t_883[k] = f_15 * sii_518[k]
                   + f_3 * pc_y[k] * ski_686[k];
    }

#pragma omp simd aligned(t_884, t_885, t_886, t_887, pb_x, pc_x, sik0_884, sii_692, sii_693, \
                         sii_694, sii_695, sik1_884, ski_693, ski_694, \
                         ski_695 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_884[k] = pb_x[k] * sik0_884[k]
                   + f_14 * sii_692[k]
                   - f_12 * pc_x[k] * sik1_884[k];

        t_885[k] = f_13 * sii_693[k]
                   + f_3 * pc_x[k] * ski_693[k];

        t_886[k] = f_13 * sii_694[k]
                   + f_3 * pc_x[k] * ski_694[k];

        t_887[k] = f_13 * sii_695[k]
                   + f_3 * pc_x[k] * ski_695[k];
    }

#pragma omp simd aligned(t_888, t_889, t_890, t_891, pc_x, sii_696, sii_697, sii_698, sii_699, \
                         ski_696, ski_697, ski_698, ski_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_888[k] = f_13 * sii_696[k]
                   + f_3 * pc_x[k] * ski_696[k];

        t_889[k] = f_13 * sii_697[k]
                   + f_3 * pc_x[k] * ski_697[k];

        t_890[k] = f_13 * sii_698[k]
                   + f_3 * pc_x[k] * ski_698[k];

        t_891[k] = f_13 * sii_699[k]
                   + f_3 * pc_x[k] * ski_699[k];
    }

#pragma omp simd aligned(t_892, t_893, t_894, t_895, pb_x, pc_x, pc_z, sik0_892, sik0_894, \
                         sik0_895, sii_497, sik1_892, sik1_894, sik1_895, \
                         ski_693 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_892[k] = pb_x[k] * sik0_892[k]
                   - f_12 * pc_x[k] * sik1_892[k];

        t_893[k] = f_15 * sii_497[k]
                   + f_3 * pc_z[k] * ski_693[k];

        t_894[k] = pb_x[k] * sik0_894[k]
                   - f_12 * pc_x[k] * sik1_894[k];

        t_895[k] = pb_x[k] * sik0_895[k]
                   - f_12 * pc_x[k] * sik1_895[k];
    }

#pragma omp simd aligned(t_896, t_897, t_898, t_899, pb_x, pc_x, pc_y, sik0_896, sik0_897, \
                         sik0_899, sii_531, sik1_896, sik1_897, sik1_899, \
                         ski_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_896[k] = pb_x[k] * sik0_896[k]
                   - f_12 * pc_x[k] * sik1_896[k];

        t_897[k] = pb_x[k] * sik0_897[k]
                   - f_12 * pc_x[k] * sik1_897[k];

        t_898[k] = f_15 * sii_531[k]
                   + f_3 * pc_y[k] * ski_699[k];

        t_899[k] = pb_x[k] * sik0_899[k]
                   - f_12 * pc_x[k] * sik1_899[k];
    }

#pragma omp simd aligned(t_900, t_901, t_902, pb_x, pc_x, pc_y, pc_z, sik0_900, sii_504, \
                         sii_532, sii_700, sik1_900, ski_700 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_900[k] = pb_x[k] * sik0_900[k]
                   + f_0 * sii_700[k]
                   - f_12 * pc_x[k] * sik1_900[k];

        t_901[k] = f_14 * sii_532[k]
                   + f_3 * pc_y[k] * ski_700[k];

        t_902[k] = f_16 * sii_504[k]
                   + f_3 * pc_z[k] * ski_700[k];
    }

#pragma omp simd aligned(t_903, t_904, t_905, pb_x, pc_x, pc_y, sik0_903, sik0_905, sii_534, \
                         sii_703, sii_705, sik1_903, sik1_905, \
                         ski_702 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_903[k] = pb_x[k] * sik0_903[k]
                   + f_17 * sii_703[k]
                   - f_12 * pc_x[k] * sik1_903[k];

        t_904[k] = f_14 * sii_534[k]
                   + f_3 * pc_y[k] * ski_702[k];

        t_905[k] = pb_x[k] * sik0_905[k]
                   + f_17 * sii_705[k]
                   - f_12 * pc_x[k] * sik1_905[k];
    }

#pragma omp simd aligned(t_906, t_907, t_908, pb_x, pc_x, pc_y, pc_z, sik0_906, sii_507, \
                         sii_537, sii_706, sik1_906, ski_703, ski_705 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_906[k] = pb_x[k] * sik0_906[k]
                   + f_16 * sii_706[k]
                   - f_12 * pc_x[k] * sik1_906[k];

        t_907[k] = f_16 * sii_507[k]
                   + f_3 * pc_z[k] * ski_703[k];

        t_908[k] = f_14 * sii_537[k]
                   + f_3 * pc_y[k] * ski_705[k];
    }

#pragma omp simd aligned(t_909, t_910, t_911, pb_x, pc_x, pc_z, sik0_909, sik0_910, sii_510, \
                         sii_709, sii_710, sik1_909, sik1_910, \
                         ski_706 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_909[k] = pb_x[k] * sik0_909[k]
                   + f_16 * sii_709[k]
                   - f_12 * pc_x[k] * sik1_909[k];

        t_910[k] = pb_x[k] * sik0_910[k]
                   + f_15 * sii_710[k]
                   - f_12 * pc_x[k] * sik1_910[k];

        t_911[k] = f_16 * sii_510[k]
                   + f_3 * pc_z[k] * ski_706[k];
    }
}

static auto
compute_prim_skk_three_center_electron_repulsion_0_piece8(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sik0,
                                                          const size_t sii, const size_t sik1,
                                                          const size_t skh0, const size_t skh1,
                                                          const size_t ski, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / q;
    const auto f_1 = 3.0 / gamma;
    const auto f_2 = 3.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.0 / gamma;
    const auto f_5 = 2.0 * p / (gamma * q);
    const auto f_6 = 1.5 / gamma;
    const auto f_7 = 1.5 * p / (gamma * q);
    const auto f_8 = 1.0 / gamma;
    const auto f_9 = p / (gamma * q);
    const auto f_10 = 0.5 / gamma;
    const auto f_11 = 0.5 * p / (gamma * q);
    const auto f_12 = gamma / q;
    const auto f_13 = 0.5 / q;
    const auto f_14 = 1.0 / q;
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;
    const auto f_18 = 3.0 / q;

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
    auto *t_1034 = buffer.data(target + 1034);
    auto *t_1035 = buffer.data(target + 1035);
    auto *t_1036 = buffer.data(target + 1036);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sik0_720 = buffer.data(sik0 + 720);
    const auto *sik0_725 = buffer.data(sik0 + 725);
    const auto *sik0_729 = buffer.data(sik0 + 729);
    const auto *sik0_734 = buffer.data(sik0 + 734);
    const auto *sik0_740 = buffer.data(sik0 + 740);
    const auto *sik0_912 = buffer.data(sik0 + 912);
    const auto *sik0_914 = buffer.data(sik0 + 914);
    const auto *sik0_915 = buffer.data(sik0 + 915);
    const auto *sik0_917 = buffer.data(sik0 + 917);
    const auto *sik0_918 = buffer.data(sik0 + 918);
    const auto *sik0_920 = buffer.data(sik0 + 920);
    const auto *sik0_928 = buffer.data(sik0 + 928);
    const auto *sik0_930 = buffer.data(sik0 + 930);
    const auto *sik0_931 = buffer.data(sik0 + 931);
    const auto *sik0_932 = buffer.data(sik0 + 932);
    const auto *sik0_933 = buffer.data(sik0 + 933);
    const auto *sik0_935 = buffer.data(sik0 + 935);
    const auto *sik0_939 = buffer.data(sik0 + 939);
    const auto *sik0_942 = buffer.data(sik0 + 942);
    const auto *sik0_946 = buffer.data(sik0 + 946);
    const auto *sik0_948 = buffer.data(sik0 + 948);
    const auto *sik0_951 = buffer.data(sik0 + 951);
    const auto *sik0_953 = buffer.data(sik0 + 953);
    const auto *sik0_954 = buffer.data(sik0 + 954);
    const auto *sik0_964 = buffer.data(sik0 + 964);
    const auto *sik0_966 = buffer.data(sik0 + 966);
    const auto *sik0_967 = buffer.data(sik0 + 967);
    const auto *sik0_968 = buffer.data(sik0 + 968);
    const auto *sik0_969 = buffer.data(sik0 + 969);
    const auto *sik0_971 = buffer.data(sik0 + 971);
    const auto *sik0_972 = buffer.data(sik0 + 972);
    const auto *sik0_975 = buffer.data(sik0 + 975);
    const auto *sik0_977 = buffer.data(sik0 + 977);
    const auto *sik0_978 = buffer.data(sik0 + 978);
    const auto *sik0_981 = buffer.data(sik0 + 981);
    const auto *sik0_982 = buffer.data(sik0 + 982);
    const auto *sik0_984 = buffer.data(sik0 + 984);
    const auto *sik0_986 = buffer.data(sik0 + 986);
    const auto *sik0_987 = buffer.data(sik0 + 987);
    const auto *sik0_989 = buffer.data(sik0 + 989);
    const auto *sik0_990 = buffer.data(sik0 + 990);
    const auto *sik0_992 = buffer.data(sik0 + 992);
    const auto *sik0_1000 = buffer.data(sik0 + 1000);
    const auto *sik0_1002 = buffer.data(sik0 + 1002);
    const auto *sik0_1003 = buffer.data(sik0 + 1003);
    const auto *sik0_1004 = buffer.data(sik0 + 1004);
    const auto *sik0_1005 = buffer.data(sik0 + 1005);
    const auto *sik0_1007 = buffer.data(sik0 + 1007);

    const auto *sii_514 = buffer.data(sii + 514);
    const auto *sii_525 = buffer.data(sii + 525);
    const auto *sii_532 = buffer.data(sii + 532);
    const auto *sii_535 = buffer.data(sii + 535);
    const auto *sii_538 = buffer.data(sii + 538);
    const auto *sii_541 = buffer.data(sii + 541);
    const auto *sii_542 = buffer.data(sii + 542);
    const auto *sii_546 = buffer.data(sii + 546);
    const auto *sii_553 = buffer.data(sii + 553);
    const auto *sii_559 = buffer.data(sii + 559);
    const auto *sii_560 = buffer.data(sii + 560);
    const auto *sii_562 = buffer.data(sii + 562);
    const auto *sii_563 = buffer.data(sii + 563);
    const auto *sii_565 = buffer.data(sii + 565);
    const auto *sii_566 = buffer.data(sii + 566);
    const auto *sii_569 = buffer.data(sii + 569);
    const auto *sii_570 = buffer.data(sii + 570);
    const auto *sii_574 = buffer.data(sii + 574);
    const auto *sii_581 = buffer.data(sii + 581);
    const auto *sii_587 = buffer.data(sii + 587);
    const auto *sii_588 = buffer.data(sii + 588);
    const auto *sii_590 = buffer.data(sii + 590);
    const auto *sii_593 = buffer.data(sii + 593);
    const auto *sii_597 = buffer.data(sii + 597);
    const auto *sii_602 = buffer.data(sii + 602);
    const auto *sii_609 = buffer.data(sii + 609);
    const auto *sii_712 = buffer.data(sii + 712);
    const auto *sii_714 = buffer.data(sii + 714);
    const auto *sii_715 = buffer.data(sii + 715);
    const auto *sii_717 = buffer.data(sii + 717);
    const auto *sii_718 = buffer.data(sii + 718);
    const auto *sii_720 = buffer.data(sii + 720);
    const auto *sii_721 = buffer.data(sii + 721);
    const auto *sii_722 = buffer.data(sii + 722);
    const auto *sii_723 = buffer.data(sii + 723);
    const auto *sii_724 = buffer.data(sii + 724);
    const auto *sii_725 = buffer.data(sii + 725);
    const auto *sii_726 = buffer.data(sii + 726);
    const auto *sii_727 = buffer.data(sii + 727);
    const auto *sii_731 = buffer.data(sii + 731);
    const auto *sii_734 = buffer.data(sii + 734);
    const auto *sii_738 = buffer.data(sii + 738);
    const auto *sii_740 = buffer.data(sii + 740);
    const auto *sii_743 = buffer.data(sii + 743);
    const auto *sii_745 = buffer.data(sii + 745);
    const auto *sii_746 = buffer.data(sii + 746);
    const auto *sii_749 = buffer.data(sii + 749);
    const auto *sii_750 = buffer.data(sii + 750);
    const auto *sii_751 = buffer.data(sii + 751);
    const auto *sii_752 = buffer.data(sii + 752);
    const auto *sii_753 = buffer.data(sii + 753);
    const auto *sii_754 = buffer.data(sii + 754);
    const auto *sii_755 = buffer.data(sii + 755);
    const auto *sii_756 = buffer.data(sii + 756);
    const auto *sii_759 = buffer.data(sii + 759);
    const auto *sii_761 = buffer.data(sii + 761);
    const auto *sii_762 = buffer.data(sii + 762);
    const auto *sii_765 = buffer.data(sii + 765);
    const auto *sii_766 = buffer.data(sii + 766);
    const auto *sii_768 = buffer.data(sii + 768);
    const auto *sii_770 = buffer.data(sii + 770);
    const auto *sii_771 = buffer.data(sii + 771);
    const auto *sii_773 = buffer.data(sii + 773);
    const auto *sii_774 = buffer.data(sii + 774);
    const auto *sii_776 = buffer.data(sii + 776);
    const auto *sii_777 = buffer.data(sii + 777);
    const auto *sii_778 = buffer.data(sii + 778);
    const auto *sii_779 = buffer.data(sii + 779);
    const auto *sii_780 = buffer.data(sii + 780);
    const auto *sii_781 = buffer.data(sii + 781);
    const auto *sii_782 = buffer.data(sii + 782);
    const auto *sii_783 = buffer.data(sii + 783);

    const auto *sik1_720 = buffer.data(sik1 + 720);
    const auto *sik1_725 = buffer.data(sik1 + 725);
    const auto *sik1_729 = buffer.data(sik1 + 729);
    const auto *sik1_734 = buffer.data(sik1 + 734);
    const auto *sik1_740 = buffer.data(sik1 + 740);
    const auto *sik1_912 = buffer.data(sik1 + 912);
    const auto *sik1_914 = buffer.data(sik1 + 914);
    const auto *sik1_915 = buffer.data(sik1 + 915);
    const auto *sik1_917 = buffer.data(sik1 + 917);
    const auto *sik1_918 = buffer.data(sik1 + 918);
    const auto *sik1_920 = buffer.data(sik1 + 920);
    const auto *sik1_928 = buffer.data(sik1 + 928);
    const auto *sik1_930 = buffer.data(sik1 + 930);
    const auto *sik1_931 = buffer.data(sik1 + 931);
    const auto *sik1_932 = buffer.data(sik1 + 932);
    const auto *sik1_933 = buffer.data(sik1 + 933);
    const auto *sik1_935 = buffer.data(sik1 + 935);
    const auto *sik1_939 = buffer.data(sik1 + 939);
    const auto *sik1_942 = buffer.data(sik1 + 942);
    const auto *sik1_946 = buffer.data(sik1 + 946);
    const auto *sik1_948 = buffer.data(sik1 + 948);
    const auto *sik1_951 = buffer.data(sik1 + 951);
    const auto *sik1_953 = buffer.data(sik1 + 953);
    const auto *sik1_954 = buffer.data(sik1 + 954);
    const auto *sik1_964 = buffer.data(sik1 + 964);
    const auto *sik1_966 = buffer.data(sik1 + 966);
    const auto *sik1_967 = buffer.data(sik1 + 967);
    const auto *sik1_968 = buffer.data(sik1 + 968);
    const auto *sik1_969 = buffer.data(sik1 + 969);
    const auto *sik1_971 = buffer.data(sik1 + 971);
    const auto *sik1_972 = buffer.data(sik1 + 972);
    const auto *sik1_975 = buffer.data(sik1 + 975);
    const auto *sik1_977 = buffer.data(sik1 + 977);
    const auto *sik1_978 = buffer.data(sik1 + 978);
    const auto *sik1_981 = buffer.data(sik1 + 981);
    const auto *sik1_982 = buffer.data(sik1 + 982);
    const auto *sik1_984 = buffer.data(sik1 + 984);
    const auto *sik1_986 = buffer.data(sik1 + 986);
    const auto *sik1_987 = buffer.data(sik1 + 987);
    const auto *sik1_989 = buffer.data(sik1 + 989);
    const auto *sik1_990 = buffer.data(sik1 + 990);
    const auto *sik1_992 = buffer.data(sik1 + 992);
    const auto *sik1_1000 = buffer.data(sik1 + 1000);
    const auto *sik1_1002 = buffer.data(sik1 + 1002);
    const auto *sik1_1003 = buffer.data(sik1 + 1003);
    const auto *sik1_1004 = buffer.data(sik1 + 1004);
    const auto *sik1_1005 = buffer.data(sik1 + 1005);
    const auto *sik1_1007 = buffer.data(sik1 + 1007);

    const auto *skh0_588 = buffer.data(skh0 + 588);
    const auto *skh0_591 = buffer.data(skh0 + 591);
    const auto *skh0_593 = buffer.data(skh0 + 593);
    const auto *skh0_594 = buffer.data(skh0 + 594);
    const auto *skh0_597 = buffer.data(skh0 + 597);
    const auto *skh0_598 = buffer.data(skh0 + 598);
    const auto *skh0_600 = buffer.data(skh0 + 600);
    const auto *skh0_602 = buffer.data(skh0 + 602);
    const auto *skh0_603 = buffer.data(skh0 + 603);
    const auto *skh0_605 = buffer.data(skh0 + 605);
    const auto *skh0_606 = buffer.data(skh0 + 606);
    const auto *skh0_608 = buffer.data(skh0 + 608);

    const auto *skh1_588 = buffer.data(skh1 + 588);
    const auto *skh1_591 = buffer.data(skh1 + 591);
    const auto *skh1_593 = buffer.data(skh1 + 593);
    const auto *skh1_594 = buffer.data(skh1 + 594);
    const auto *skh1_597 = buffer.data(skh1 + 597);
    const auto *skh1_598 = buffer.data(skh1 + 598);
    const auto *skh1_600 = buffer.data(skh1 + 600);
    const auto *skh1_602 = buffer.data(skh1 + 602);
    const auto *skh1_603 = buffer.data(skh1 + 603);
    const auto *skh1_605 = buffer.data(skh1 + 605);
    const auto *skh1_606 = buffer.data(skh1 + 606);
    const auto *skh1_608 = buffer.data(skh1 + 608);

    const auto *ski_709 = buffer.data(ski + 709);
    const auto *ski_710 = buffer.data(ski + 710);
    const auto *ski_714 = buffer.data(ski + 714);
    const auto *ski_721 = buffer.data(ski + 721);
    const auto *ski_722 = buffer.data(ski + 722);
    const auto *ski_723 = buffer.data(ski + 723);
    const auto *ski_724 = buffer.data(ski + 724);
    const auto *ski_725 = buffer.data(ski + 725);
    const auto *ski_726 = buffer.data(ski + 726);
    const auto *ski_727 = buffer.data(ski + 727);
    const auto *ski_728 = buffer.data(ski + 728);
    const auto *ski_730 = buffer.data(ski + 730);
    const auto *ski_731 = buffer.data(ski + 731);
    const auto *ski_733 = buffer.data(ski + 733);
    const auto *ski_734 = buffer.data(ski + 734);
    const auto *ski_737 = buffer.data(ski + 737);
    const auto *ski_738 = buffer.data(ski + 738);
    const auto *ski_742 = buffer.data(ski + 742);
    const auto *ski_749 = buffer.data(ski + 749);
    const auto *ski_750 = buffer.data(ski + 750);
    const auto *ski_751 = buffer.data(ski + 751);
    const auto *ski_752 = buffer.data(ski + 752);
    const auto *ski_753 = buffer.data(ski + 753);
    const auto *ski_754 = buffer.data(ski + 754);
    const auto *ski_755 = buffer.data(ski + 755);
    const auto *ski_756 = buffer.data(ski + 756);
    const auto *ski_758 = buffer.data(ski + 758);
    const auto *ski_759 = buffer.data(ski + 759);
    const auto *ski_761 = buffer.data(ski + 761);
    const auto *ski_762 = buffer.data(ski + 762);
    const auto *ski_765 = buffer.data(ski + 765);
    const auto *ski_766 = buffer.data(ski + 766);
    const auto *ski_770 = buffer.data(ski + 770);
    const auto *ski_777 = buffer.data(ski + 777);
    const auto *ski_778 = buffer.data(ski + 778);
    const auto *ski_779 = buffer.data(ski + 779);
    const auto *ski_780 = buffer.data(ski + 780);
    const auto *ski_781 = buffer.data(ski + 781);
    const auto *ski_782 = buffer.data(ski + 782);
    const auto *ski_783 = buffer.data(ski + 783);
    const auto *ski_784 = buffer.data(ski + 784);
    const auto *ski_786 = buffer.data(ski + 786);
    const auto *ski_787 = buffer.data(ski + 787);
    const auto *ski_789 = buffer.data(ski + 789);
    const auto *ski_790 = buffer.data(ski + 790);
    const auto *ski_793 = buffer.data(ski + 793);
    const auto *ski_794 = buffer.data(ski + 794);
    const auto *ski_796 = buffer.data(ski + 796);
    const auto *ski_798 = buffer.data(ski + 798);
    const auto *ski_799 = buffer.data(ski + 799);
    const auto *ski_801 = buffer.data(ski + 801);
    const auto *ski_802 = buffer.data(ski + 802);
    const auto *ski_804 = buffer.data(ski + 804);
    const auto *ski_805 = buffer.data(ski + 805);
    const auto *ski_806 = buffer.data(ski + 806);
    const auto *ski_807 = buffer.data(ski + 807);
    const auto *ski_808 = buffer.data(ski + 808);
    const auto *ski_809 = buffer.data(ski + 809);
    const auto *ski_810 = buffer.data(ski + 810);
    const auto *ski_811 = buffer.data(ski + 811);

#pragma omp simd aligned(t_912, t_913, t_914, pb_x, pc_x, pc_y, sik0_912, sik0_914, sii_541, \
                         sii_712, sii_714, sik1_912, sik1_914, \
                         ski_709 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_912[k] = pb_x[k] * sik0_912[k]
                   + f_15 * sii_712[k]
                   - f_12 * pc_x[k] * sik1_912[k];

        t_913[k] = f_14 * sii_541[k]
                   + f_3 * pc_y[k] * ski_709[k];

        t_914[k] = pb_x[k] * sik0_914[k]
                   + f_15 * sii_714[k]
                   - f_12 * pc_x[k] * sik1_914[k];
    }

#pragma omp simd aligned(t_915, t_916, t_917, pb_x, pc_x, pc_z, sik0_915, sik0_917, sii_514, \
                         sii_715, sii_717, sik1_915, sik1_917, \
                         ski_710 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_915[k] = pb_x[k] * sik0_915[k]
                   + f_14 * sii_715[k]
                   - f_12 * pc_x[k] * sik1_915[k];

        t_916[k] = f_16 * sii_514[k]
                   + f_3 * pc_z[k] * ski_710[k];

        t_917[k] = pb_x[k] * sik0_917[k]
                   + f_14 * sii_717[k]
                   - f_12 * pc_x[k] * sik1_917[k];
    }

#pragma omp simd aligned(t_918, t_919, t_920, pb_x, pc_x, pc_y, sik0_918, sik0_920, sii_546, \
                         sii_718, sii_720, sik1_918, sik1_920, \
                         ski_714 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_918[k] = pb_x[k] * sik0_918[k]
                   + f_14 * sii_718[k]
                   - f_12 * pc_x[k] * sik1_918[k];

        t_919[k] = f_14 * sii_546[k]
                   + f_3 * pc_y[k] * ski_714[k];

        t_920[k] = pb_x[k] * sik0_920[k]
                   + f_14 * sii_720[k]
                   - f_12 * pc_x[k] * sik1_920[k];
    }

#pragma omp simd aligned(t_921, t_922, t_923, t_924, t_925, pc_x, sii_721, sii_722, sii_723, \
                         sii_724, sii_725, ski_721, ski_722, ski_723, ski_724, \
                         ski_725 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_921[k] = f_13 * sii_721[k]
                   + f_3 * pc_x[k] * ski_721[k];

        t_922[k] = f_13 * sii_722[k]
                   + f_3 * pc_x[k] * ski_722[k];

        t_923[k] = f_13 * sii_723[k]
                   + f_3 * pc_x[k] * ski_723[k];

        t_924[k] = f_13 * sii_724[k]
                   + f_3 * pc_x[k] * ski_724[k];

        t_925[k] = f_13 * sii_725[k]
                   + f_3 * pc_x[k] * ski_725[k];
    }

#pragma omp simd aligned(t_926, t_927, t_928, t_929, pb_x, pc_x, pc_z, sik0_928, sii_525, \
                         sii_726, sii_727, sik1_928, ski_721, ski_726, \
                         ski_727 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_926[k] = f_13 * sii_726[k]
                   + f_3 * pc_x[k] * ski_726[k];

        t_927[k] = f_13 * sii_727[k]
                   + f_3 * pc_x[k] * ski_727[k];

        t_928[k] = pb_x[k] * sik0_928[k]
                   - f_12 * pc_x[k] * sik1_928[k];

        t_929[k] = f_16 * sii_525[k]
                   + f_3 * pc_z[k] * ski_721[k];
    }

#pragma omp simd aligned(t_930, t_931, t_932, t_933, pb_x, pc_x, sik0_930, sik0_931, sik0_932, \
                         sik0_933, sik1_930, sik1_931, sik1_932, \
                         sik1_933 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_930[k] = pb_x[k] * sik0_930[k]
                   - f_12 * pc_x[k] * sik1_930[k];

        t_931[k] = pb_x[k] * sik0_931[k]
                   - f_12 * pc_x[k] * sik1_931[k];

        t_932[k] = pb_x[k] * sik0_932[k]
                   - f_12 * pc_x[k] * sik1_932[k];

        t_933[k] = pb_x[k] * sik0_933[k]
                   - f_12 * pc_x[k] * sik1_933[k];
    }

#pragma omp simd aligned(t_934, t_935, t_936, t_937, pb_x, pb_y, pc_x, pc_y, sik0_720, \
                         sik0_935, sii_559, sii_560, sik1_720, sik1_935, ski_727, \
                         ski_728 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_934[k] = f_14 * sii_559[k]
                   + f_3 * pc_y[k] * ski_727[k];

        t_935[k] = pb_x[k] * sik0_935[k]
                   - f_12 * pc_x[k] * sik1_935[k];

        t_936[k] = pb_y[k] * sik0_720[k]
                   - f_12 * pc_y[k] * sik1_720[k];

        t_937[k] = f_13 * sii_560[k]
                   + f_3 * pc_y[k] * ski_728[k];
    }

#pragma omp simd aligned(t_938, t_939, t_940, pb_x, pc_x, pc_y, pc_z, sik0_939, sii_532, \
                         sii_562, sii_731, sik1_939, ski_728, ski_730 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_938[k] = f_17 * sii_532[k]
                   + f_3 * pc_z[k] * ski_728[k];

        t_939[k] = pb_x[k] * sik0_939[k]
                   + f_17 * sii_731[k]
                   - f_12 * pc_x[k] * sik1_939[k];

        t_940[k] = f_13 * sii_562[k]
                   + f_3 * pc_y[k] * ski_730[k];
    }

#pragma omp simd aligned(t_941, t_942, t_943, pb_x, pb_y, pc_x, pc_y, pc_z, sik0_725, \
                         sik0_942, sii_535, sii_734, sik1_725, sik1_942, \
                         ski_731 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_941[k] = pb_y[k] * sik0_725[k]
                   - f_12 * pc_y[k] * sik1_725[k];

        t_942[k] = pb_x[k] * sik0_942[k]
                   + f_16 * sii_734[k]
                   - f_12 * pc_x[k] * sik1_942[k];

        t_943[k] = f_17 * sii_535[k]
                   + f_3 * pc_z[k] * ski_731[k];
    }

#pragma omp simd aligned(t_944, t_945, t_946, pb_x, pb_y, pc_x, pc_y, sik0_729, sik0_946, \
                         sii_565, sii_738, sik1_729, sik1_946, \
                         ski_733 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_944[k] = f_13 * sii_565[k]
                   + f_3 * pc_y[k] * ski_733[k];

        t_945[k] = pb_y[k] * sik0_729[k]
                   - f_12 * pc_y[k] * sik1_729[k];

        t_946[k] = pb_x[k] * sik0_946[k]
                   + f_15 * sii_738[k]
                   - f_12 * pc_x[k] * sik1_946[k];
    }

#pragma omp simd aligned(t_947, t_948, t_949, pb_x, pc_x, pc_y, pc_z, sik0_948, sii_538, \
                         sii_569, sii_740, sik1_948, ski_734, ski_737 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_947[k] = f_17 * sii_538[k]
                   + f_3 * pc_z[k] * ski_734[k];

        t_948[k] = pb_x[k] * sik0_948[k]
                   + f_15 * sii_740[k]
                   - f_12 * pc_x[k] * sik1_948[k];

        t_949[k] = f_13 * sii_569[k]
                   + f_3 * pc_y[k] * ski_737[k];
    }

#pragma omp simd aligned(t_950, t_951, t_952, pb_x, pb_y, pc_x, pc_y, pc_z, sik0_734, \
                         sik0_951, sii_542, sii_743, sik1_734, sik1_951, \
                         ski_738 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_950[k] = pb_y[k] * sik0_734[k]
                   - f_12 * pc_y[k] * sik1_734[k];

        t_951[k] = pb_x[k] * sik0_951[k]
                   + f_14 * sii_743[k]
                   - f_12 * pc_x[k] * sik1_951[k];

        t_952[k] = f_17 * sii_542[k]
                   + f_3 * pc_z[k] * ski_738[k];
    }

#pragma omp simd aligned(t_953, t_954, t_955, pb_x, pc_x, pc_y, sik0_953, sik0_954, sii_574, \
                         sii_745, sii_746, sik1_953, sik1_954, \
                         ski_742 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_953[k] = pb_x[k] * sik0_953[k]
                   + f_14 * sii_745[k]
                   - f_12 * pc_x[k] * sik1_953[k];

        t_954[k] = pb_x[k] * sik0_954[k]
                   + f_14 * sii_746[k]
                   - f_12 * pc_x[k] * sik1_954[k];

        t_955[k] = f_13 * sii_574[k]
                   + f_3 * pc_y[k] * ski_742[k];
    }

#pragma omp simd aligned(t_956, t_957, t_958, t_959, pb_y, pc_x, pc_y, sik0_740, sii_749, \
                         sii_750, sii_751, sik1_740, ski_749, ski_750, \
                         ski_751 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_956[k] = pb_y[k] * sik0_740[k]
                   - f_12 * pc_y[k] * sik1_740[k];

        t_957[k] = f_13 * sii_749[k]
                   + f_3 * pc_x[k] * ski_749[k];

        t_958[k] = f_13 * sii_750[k]
                   + f_3 * pc_x[k] * ski_750[k];

        t_959[k] = f_13 * sii_751[k]
                   + f_3 * pc_x[k] * ski_751[k];
    }

#pragma omp simd aligned(t_960, t_961, t_962, t_963, pc_x, sii_752, sii_753, sii_754, sii_755, \
                         ski_752, ski_753, ski_754, ski_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_960[k] = f_13 * sii_752[k]
                   + f_3 * pc_x[k] * ski_752[k];

        t_961[k] = f_13 * sii_753[k]
                   + f_3 * pc_x[k] * ski_753[k];

        t_962[k] = f_13 * sii_754[k]
                   + f_3 * pc_x[k] * ski_754[k];

        t_963[k] = f_13 * sii_755[k]
                   + f_3 * pc_x[k] * ski_755[k];
    }

#pragma omp simd aligned(t_964, t_965, t_966, t_967, pb_x, pc_x, pc_z, sik0_964, sik0_966, \
                         sik0_967, sii_553, sik1_964, sik1_966, sik1_967, \
                         ski_749 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_964[k] = pb_x[k] * sik0_964[k]
                   - f_12 * pc_x[k] * sik1_964[k];

        t_965[k] = f_17 * sii_553[k]
                   + f_3 * pc_z[k] * ski_749[k];

        t_966[k] = pb_x[k] * sik0_966[k]
                   - f_12 * pc_x[k] * sik1_966[k];

        t_967[k] = pb_x[k] * sik0_967[k]
                   - f_12 * pc_x[k] * sik1_967[k];
    }

#pragma omp simd aligned(t_968, t_969, t_970, t_971, pb_x, pc_x, pc_y, sik0_968, sik0_969, \
                         sik0_971, sii_587, sik1_968, sik1_969, sik1_971, \
                         ski_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_968[k] = pb_x[k] * sik0_968[k]
                   - f_12 * pc_x[k] * sik1_968[k];

        t_969[k] = pb_x[k] * sik0_969[k]
                   - f_12 * pc_x[k] * sik1_969[k];

        t_970[k] = f_13 * sii_587[k]
                   + f_3 * pc_y[k] * ski_755[k];

        t_971[k] = pb_x[k] * sik0_971[k]
                   - f_12 * pc_x[k] * sik1_971[k];
    }

#pragma omp simd aligned(t_972, t_973, t_974, t_975, pb_x, pc_x, pc_y, pc_z, sik0_972, \
                         sik0_975, sii_560, sii_756, sii_759, sik1_972, sik1_975, \
                         ski_756 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_972[k] = pb_x[k] * sik0_972[k]
                   + f_0 * sii_756[k]
                   - f_12 * pc_x[k] * sik1_972[k];

        t_973[k] = f_3 * pc_y[k] * ski_756[k];

        t_974[k] = f_18 * sii_560[k]
                   + f_3 * pc_z[k] * ski_756[k];

        t_975[k] = pb_x[k] * sik0_975[k]
                   + f_17 * sii_759[k]
                   - f_12 * pc_x[k] * sik1_975[k];
    }

#pragma omp simd aligned(t_976, t_977, t_978, pb_x, pc_x, pc_y, sik0_977, sik0_978, sii_761, \
                         sii_762, sik1_977, sik1_978, ski_758 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_976[k] = f_3 * pc_y[k] * ski_758[k];

        t_977[k] = pb_x[k] * sik0_977[k]
                   + f_17 * sii_761[k]
                   - f_12 * pc_x[k] * sik1_977[k];

        t_978[k] = pb_x[k] * sik0_978[k]
                   + f_16 * sii_762[k]
                   - f_12 * pc_x[k] * sik1_978[k];
    }

#pragma omp simd aligned(t_979, t_980, t_981, pb_x, pc_x, pc_y, pc_z, sik0_981, sii_563, \
                         sii_765, sik1_981, ski_759, ski_761 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_979[k] = f_18 * sii_563[k]
                   + f_3 * pc_z[k] * ski_759[k];

        t_980[k] = f_3 * pc_y[k] * ski_761[k];

        t_981[k] = pb_x[k] * sik0_981[k]
                   + f_16 * sii_765[k]
                   - f_12 * pc_x[k] * sik1_981[k];
    }

#pragma omp simd aligned(t_982, t_983, t_984, pb_x, pc_x, pc_z, sik0_982, sik0_984, sii_566, \
                         sii_766, sii_768, sik1_982, sik1_984, \
                         ski_762 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_982[k] = pb_x[k] * sik0_982[k]
                   + f_15 * sii_766[k]
                   - f_12 * pc_x[k] * sik1_982[k];

        t_983[k] = f_18 * sii_566[k]
                   + f_3 * pc_z[k] * ski_762[k];

        t_984[k] = pb_x[k] * sik0_984[k]
                   + f_15 * sii_768[k]
                   - f_12 * pc_x[k] * sik1_984[k];
    }

#pragma omp simd aligned(t_985, t_986, t_987, pb_x, pc_x, pc_y, sik0_986, sik0_987, sii_770, \
                         sii_771, sik1_986, sik1_987, ski_765 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_985[k] = f_3 * pc_y[k] * ski_765[k];

        t_986[k] = pb_x[k] * sik0_986[k]
                   + f_15 * sii_770[k]
                   - f_12 * pc_x[k] * sik1_986[k];

        t_987[k] = pb_x[k] * sik0_987[k]
                   + f_14 * sii_771[k]
                   - f_12 * pc_x[k] * sik1_987[k];
    }

#pragma omp simd aligned(t_988, t_989, t_990, pb_x, pc_x, pc_z, sik0_989, sik0_990, sii_570, \
                         sii_773, sii_774, sik1_989, sik1_990, \
                         ski_766 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_988[k] = f_18 * sii_570[k]
                   + f_3 * pc_z[k] * ski_766[k];

        t_989[k] = pb_x[k] * sik0_989[k]
                   + f_14 * sii_773[k]
                   - f_12 * pc_x[k] * sik1_989[k];

        t_990[k] = pb_x[k] * sik0_990[k]
                   + f_14 * sii_774[k]
                   - f_12 * pc_x[k] * sik1_990[k];
    }

#pragma omp simd aligned(t_991, t_992, t_993, t_994, pb_x, pc_x, pc_y, sik0_992, sii_776, \
                         sii_777, sii_778, sik1_992, ski_770, ski_777, \
                         ski_778 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_991[k] = f_3 * pc_y[k] * ski_770[k];

        t_992[k] = pb_x[k] * sik0_992[k]
                   + f_14 * sii_776[k]
                   - f_12 * pc_x[k] * sik1_992[k];

        t_993[k] = f_13 * sii_777[k]
                   + f_3 * pc_x[k] * ski_777[k];

        t_994[k] = f_13 * sii_778[k]
                   + f_3 * pc_x[k] * ski_778[k];
    }

#pragma omp simd aligned(t_995, t_996, t_997, t_998, t_999, pc_x, sii_779, sii_780, sii_781, \
                         sii_782, sii_783, ski_779, ski_780, ski_781, ski_782, \
                         ski_783 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_995[k] = f_13 * sii_779[k]
                   + f_3 * pc_x[k] * ski_779[k];

        t_996[k] = f_13 * sii_780[k]
                   + f_3 * pc_x[k] * ski_780[k];

        t_997[k] = f_13 * sii_781[k]
                   + f_3 * pc_x[k] * ski_781[k];

        t_998[k] = f_13 * sii_782[k]
                   + f_3 * pc_x[k] * ski_782[k];

        t_999[k] = f_13 * sii_783[k]
                   + f_3 * pc_x[k] * ski_783[k];
    }

#pragma omp simd aligned(t_1000, t_1001, t_1002, t_1003, pb_x, pc_x, pc_z, sik0_1000, \
                         sik0_1002, sik0_1003, sii_581, sik1_1000, sik1_1002, sik1_1003, \
                         ski_777 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1000[k] = pb_x[k] * sik0_1000[k]
                    - f_12 * pc_x[k] * sik1_1000[k];

        t_1001[k] = f_18 * sii_581[k]
                    + f_3 * pc_z[k] * ski_777[k];

        t_1002[k] = pb_x[k] * sik0_1002[k]
                    - f_12 * pc_x[k] * sik1_1002[k];

        t_1003[k] = pb_x[k] * sik0_1003[k]
                    - f_12 * pc_x[k] * sik1_1003[k];
    }

#pragma omp simd aligned(t_1004, t_1005, t_1006, t_1007, pb_x, pc_x, pc_y, sik0_1004, \
                         sik0_1005, sik0_1007, sik1_1004, sik1_1005, sik1_1007, \
                         ski_783 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1004[k] = pb_x[k] * sik0_1004[k]
                    - f_12 * pc_x[k] * sik1_1004[k];

        t_1005[k] = pb_x[k] * sik0_1005[k]
                    - f_12 * pc_x[k] * sik1_1005[k];

        t_1006[k] = f_3 * pc_y[k] * ski_783[k];

        t_1007[k] = pb_x[k] * sik0_1007[k]
                    - f_12 * pc_x[k] * sik1_1007[k];
    }

#pragma omp simd aligned(t_1008, t_1009, t_1010, t_1011, pc_x, pc_y, pc_z, sii_588, skh0_588, \
                         skh0_591, skh1_588, skh1_591, ski_784, \
                         ski_787 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1008[k] = f_1 * skh0_588[k]
                    - f_2 * skh1_588[k]
                    + f_3 * pc_x[k] * ski_784[k];

        t_1009[k] = f_0 * sii_588[k]
                    + f_3 * pc_y[k] * ski_784[k];

        t_1010[k] = f_3 * pc_z[k] * ski_784[k];

        t_1011[k] = f_4 * skh0_591[k]
                    - f_5 * skh1_591[k]
                    + f_3 * pc_x[k] * ski_787[k];
    }

#pragma omp simd aligned(t_1012, t_1013, t_1014, t_1015, pc_x, pc_y, pc_z, sii_590, skh0_593, \
                         skh0_594, skh1_593, skh1_594, ski_786, ski_787, ski_789, \
                         ski_790 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1012[k] = f_0 * sii_590[k]
                    + f_3 * pc_y[k] * ski_786[k];

        t_1013[k] = f_4 * skh0_593[k]
                    - f_5 * skh1_593[k]
                    + f_3 * pc_x[k] * ski_789[k];

        t_1014[k] = f_6 * skh0_594[k]
                    - f_7 * skh1_594[k]
                    + f_3 * pc_x[k] * ski_790[k];

        t_1015[k] = f_3 * pc_z[k] * ski_787[k];
    }

#pragma omp simd aligned(t_1016, t_1017, t_1018, t_1019, pc_x, pc_y, pc_z, sii_593, skh0_597, \
                         skh0_598, skh1_597, skh1_598, ski_789, ski_790, ski_793, \
                         ski_794 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1016[k] = f_0 * sii_593[k]
                    + f_3 * pc_y[k] * ski_789[k];

        t_1017[k] = f_6 * skh0_597[k]
                    - f_7 * skh1_597[k]
                    + f_3 * pc_x[k] * ski_793[k];

        t_1018[k] = f_8 * skh0_598[k]
                    - f_9 * skh1_598[k]
                    + f_3 * pc_x[k] * ski_794[k];

        t_1019[k] = f_3 * pc_z[k] * ski_790[k];
    }

#pragma omp simd aligned(t_1020, t_1021, t_1022, pc_x, pc_y, sii_597, skh0_600, skh0_602, \
                         skh1_600, skh1_602, ski_793, ski_796, \
                         ski_798 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1020[k] = f_8 * skh0_600[k]
                    - f_9 * skh1_600[k]
                    + f_3 * pc_x[k] * ski_796[k];

        t_1021[k] = f_0 * sii_597[k]
                    + f_3 * pc_y[k] * ski_793[k];

        t_1022[k] = f_8 * skh0_602[k]
                    - f_9 * skh1_602[k]
                    + f_3 * pc_x[k] * ski_798[k];
    }

#pragma omp simd aligned(t_1023, t_1024, t_1025, t_1026, pc_x, pc_z, skh0_603, skh0_605, \
                         skh0_606, skh1_603, skh1_605, skh1_606, ski_794, ski_799, ski_801, \
                         ski_802 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1023[k] = f_10 * skh0_603[k]
                    - f_11 * skh1_603[k]
                    + f_3 * pc_x[k] * ski_799[k];

        t_1024[k] = f_3 * pc_z[k] * ski_794[k];

        t_1025[k] = f_10 * skh0_605[k]
                    - f_11 * skh1_605[k]
                    + f_3 * pc_x[k] * ski_801[k];

        t_1026[k] = f_10 * skh0_606[k]
                    - f_11 * skh1_606[k]
                    + f_3 * pc_x[k] * ski_802[k];
    }

#pragma omp simd aligned(t_1027, t_1028, t_1029, t_1030, t_1031, pc_x, pc_y, sii_602, \
                         skh0_608, skh1_608, ski_798, ski_804, ski_805, ski_806, \
                         ski_807 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1027[k] = f_0 * sii_602[k]
                    + f_3 * pc_y[k] * ski_798[k];

        t_1028[k] = f_10 * skh0_608[k]
                    - f_11 * skh1_608[k]
                    + f_3 * pc_x[k] * ski_804[k];

        t_1029[k] = f_3 * pc_x[k] * ski_805[k];

        t_1030[k] = f_3 * pc_x[k] * ski_806[k];

        t_1031[k] = f_3 * pc_x[k] * ski_807[k];
    }

#pragma omp simd aligned(t_1032, t_1033, t_1034, t_1035, t_1036, pc_x, pc_y, sii_609, \
                         skh0_603, skh1_603, ski_805, ski_808, ski_809, ski_810, \
                         ski_811 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1032[k] = f_3 * pc_x[k] * ski_808[k];

        t_1033[k] = f_3 * pc_x[k] * ski_809[k];

        t_1034[k] = f_3 * pc_x[k] * ski_810[k];

        t_1035[k] = f_3 * pc_x[k] * ski_811[k];

        t_1036[k] = f_0 * sii_609[k]
                    + f_1 * skh0_603[k]
                    - f_2 * skh1_603[k]
                    + f_3 * pc_y[k] * ski_805[k];
    }
}

static auto
compute_prim_skk_three_center_electron_repulsion_0_piece9(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sik0,
                                                          const size_t sii, const size_t sik1,
                                                          const size_t skh0, const size_t skh1,
                                                          const size_t ski, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / q;
    const auto f_1 = 3.0 / gamma;
    const auto f_2 = 3.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.0 / gamma;
    const auto f_5 = 2.0 * p / (gamma * q);
    const auto f_6 = 1.5 / gamma;
    const auto f_7 = 1.5 * p / (gamma * q);
    const auto f_8 = 1.0 / gamma;
    const auto f_9 = p / (gamma * q);
    const auto f_10 = 0.5 / gamma;
    const auto f_11 = 0.5 * p / (gamma * q);
    const auto f_12 = gamma / q;
    const auto f_13 = 0.5 / q;
    const auto f_14 = 1.0 / q;
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;
    const auto f_18 = 3.0 / q;

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
    auto *t_1149 = buffer.data(target + 1149);
    auto *t_1150 = buffer.data(target + 1150);
    auto *t_1151 = buffer.data(target + 1151);
    auto *t_1152 = buffer.data(target + 1152);
    auto *t_1153 = buffer.data(target + 1153);
    auto *t_1154 = buffer.data(target + 1154);
    auto *t_1155 = buffer.data(target + 1155);
    auto *t_1156 = buffer.data(target + 1156);

    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sik0_756 = buffer.data(sik0 + 756);
    const auto *sik0_759 = buffer.data(sik0 + 759);
    const auto *sik0_762 = buffer.data(sik0 + 762);
    const auto *sik0_766 = buffer.data(sik0 + 766);
    const auto *sik0_771 = buffer.data(sik0 + 771);
    const auto *sik0_784 = buffer.data(sik0 + 784);
    const auto *sik0_786 = buffer.data(sik0 + 786);
    const auto *sik0_787 = buffer.data(sik0 + 787);
    const auto *sik0_788 = buffer.data(sik0 + 788);
    const auto *sik0_789 = buffer.data(sik0 + 789);

    const auto *sii_588 = buffer.data(sii + 588);
    const auto *sii_591 = buffer.data(sii + 591);
    const auto *sii_594 = buffer.data(sii + 594);
    const auto *sii_598 = buffer.data(sii + 598);
    const auto *sii_609 = buffer.data(sii + 609);
    const auto *sii_610 = buffer.data(sii + 610);
    const auto *sii_611 = buffer.data(sii + 611);
    const auto *sii_612 = buffer.data(sii + 612);
    const auto *sii_613 = buffer.data(sii + 613);
    const auto *sii_614 = buffer.data(sii + 614);
    const auto *sii_615 = buffer.data(sii + 615);
    const auto *sii_616 = buffer.data(sii + 616);
    const auto *sii_618 = buffer.data(sii + 618);
    const auto *sii_619 = buffer.data(sii + 619);
    const auto *sii_621 = buffer.data(sii + 621);
    const auto *sii_622 = buffer.data(sii + 622);
    const auto *sii_625 = buffer.data(sii + 625);
    const auto *sii_626 = buffer.data(sii + 626);
    const auto *sii_630 = buffer.data(sii + 630);
    const auto *sii_637 = buffer.data(sii + 637);
    const auto *sii_643 = buffer.data(sii + 643);
    const auto *sii_644 = buffer.data(sii + 644);
    const auto *sii_646 = buffer.data(sii + 646);
    const auto *sii_647 = buffer.data(sii + 647);
    const auto *sii_649 = buffer.data(sii + 649);
    const auto *sii_650 = buffer.data(sii + 650);
    const auto *sii_653 = buffer.data(sii + 653);
    const auto *sii_654 = buffer.data(sii + 654);
    const auto *sii_658 = buffer.data(sii + 658);
    const auto *sii_665 = buffer.data(sii + 665);
    const auto *sii_667 = buffer.data(sii + 667);
    const auto *sii_668 = buffer.data(sii + 668);
    const auto *sii_669 = buffer.data(sii + 669);
    const auto *sii_670 = buffer.data(sii + 670);
    const auto *sii_671 = buffer.data(sii + 671);
    const auto *sii_672 = buffer.data(sii + 672);
    const auto *sii_674 = buffer.data(sii + 674);
    const auto *sii_677 = buffer.data(sii + 677);
    const auto *sii_681 = buffer.data(sii + 681);
    const auto *sii_686 = buffer.data(sii + 686);
    const auto *sii_693 = buffer.data(sii + 693);
    const auto *sii_695 = buffer.data(sii + 695);
    const auto *sii_696 = buffer.data(sii + 696);
    const auto *sii_697 = buffer.data(sii + 697);
    const auto *sii_698 = buffer.data(sii + 698);
    const auto *sii_699 = buffer.data(sii + 699);
    const auto *sii_700 = buffer.data(sii + 700);
    const auto *sii_702 = buffer.data(sii + 702);

    const auto *sik1_756 = buffer.data(sik1 + 756);
    const auto *sik1_759 = buffer.data(sik1 + 759);
    const auto *sik1_762 = buffer.data(sik1 + 762);
    const auto *sik1_766 = buffer.data(sik1 + 766);
    const auto *sik1_771 = buffer.data(sik1 + 771);
    const auto *sik1_784 = buffer.data(sik1 + 784);
    const auto *sik1_786 = buffer.data(sik1 + 786);
    const auto *sik1_787 = buffer.data(sik1 + 787);
    const auto *sik1_788 = buffer.data(sik1 + 788);
    const auto *sik1_789 = buffer.data(sik1 + 789);

    const auto *skh0_605 = buffer.data(skh0 + 605);
    const auto *skh0_606 = buffer.data(skh0 + 606);
    const auto *skh0_607 = buffer.data(skh0 + 607);
    const auto *skh0_608 = buffer.data(skh0 + 608);
    const auto *skh0_614 = buffer.data(skh0 + 614);
    const auto *skh0_618 = buffer.data(skh0 + 618);
    const auto *skh0_621 = buffer.data(skh0 + 621);
    const auto *skh0_623 = buffer.data(skh0 + 623);
    const auto *skh0_626 = buffer.data(skh0 + 626);
    const auto *skh0_627 = buffer.data(skh0 + 627);
    const auto *skh0_629 = buffer.data(skh0 + 629);
    const auto *skh0_630 = buffer.data(skh0 + 630);
    const auto *skh0_633 = buffer.data(skh0 + 633);
    const auto *skh0_635 = buffer.data(skh0 + 635);
    const auto *skh0_636 = buffer.data(skh0 + 636);
    const auto *skh0_639 = buffer.data(skh0 + 639);
    const auto *skh0_640 = buffer.data(skh0 + 640);
    const auto *skh0_642 = buffer.data(skh0 + 642);
    const auto *skh0_644 = buffer.data(skh0 + 644);
    const auto *skh0_645 = buffer.data(skh0 + 645);
    const auto *skh0_647 = buffer.data(skh0 + 647);
    const auto *skh0_648 = buffer.data(skh0 + 648);
    const auto *skh0_649 = buffer.data(skh0 + 649);
    const auto *skh0_650 = buffer.data(skh0 + 650);
    const auto *skh0_651 = buffer.data(skh0 + 651);
    const auto *skh0_654 = buffer.data(skh0 + 654);
    const auto *skh0_656 = buffer.data(skh0 + 656);
    const auto *skh0_657 = buffer.data(skh0 + 657);
    const auto *skh0_660 = buffer.data(skh0 + 660);
    const auto *skh0_661 = buffer.data(skh0 + 661);
    const auto *skh0_663 = buffer.data(skh0 + 663);
    const auto *skh0_665 = buffer.data(skh0 + 665);
    const auto *skh0_666 = buffer.data(skh0 + 666);
    const auto *skh0_668 = buffer.data(skh0 + 668);
    const auto *skh0_669 = buffer.data(skh0 + 669);
    const auto *skh0_670 = buffer.data(skh0 + 670);
    const auto *skh0_671 = buffer.data(skh0 + 671);
    const auto *skh0_672 = buffer.data(skh0 + 672);
    const auto *skh0_675 = buffer.data(skh0 + 675);

    const auto *skh1_605 = buffer.data(skh1 + 605);
    const auto *skh1_606 = buffer.data(skh1 + 606);
    const auto *skh1_607 = buffer.data(skh1 + 607);
    const auto *skh1_608 = buffer.data(skh1 + 608);
    const auto *skh1_614 = buffer.data(skh1 + 614);
    const auto *skh1_618 = buffer.data(skh1 + 618);
    const auto *skh1_621 = buffer.data(skh1 + 621);
    const auto *skh1_623 = buffer.data(skh1 + 623);
    const auto *skh1_626 = buffer.data(skh1 + 626);
    const auto *skh1_627 = buffer.data(skh1 + 627);
    const auto *skh1_629 = buffer.data(skh1 + 629);
    const auto *skh1_630 = buffer.data(skh1 + 630);
    const auto *skh1_633 = buffer.data(skh1 + 633);
    const auto *skh1_635 = buffer.data(skh1 + 635);
    const auto *skh1_636 = buffer.data(skh1 + 636);
    const auto *skh1_639 = buffer.data(skh1 + 639);
    const auto *skh1_640 = buffer.data(skh1 + 640);
    const auto *skh1_642 = buffer.data(skh1 + 642);
    const auto *skh1_644 = buffer.data(skh1 + 644);
    const auto *skh1_645 = buffer.data(skh1 + 645);
    const auto *skh1_647 = buffer.data(skh1 + 647);
    const auto *skh1_648 = buffer.data(skh1 + 648);
    const auto *skh1_649 = buffer.data(skh1 + 649);
    const auto *skh1_650 = buffer.data(skh1 + 650);
    const auto *skh1_651 = buffer.data(skh1 + 651);
    const auto *skh1_654 = buffer.data(skh1 + 654);
    const auto *skh1_656 = buffer.data(skh1 + 656);
    const auto *skh1_657 = buffer.data(skh1 + 657);
    const auto *skh1_660 = buffer.data(skh1 + 660);
    const auto *skh1_661 = buffer.data(skh1 + 661);
    const auto *skh1_663 = buffer.data(skh1 + 663);
    const auto *skh1_665 = buffer.data(skh1 + 665);
    const auto *skh1_666 = buffer.data(skh1 + 666);
    const auto *skh1_668 = buffer.data(skh1 + 668);
    const auto *skh1_669 = buffer.data(skh1 + 669);
    const auto *skh1_670 = buffer.data(skh1 + 670);
    const auto *skh1_671 = buffer.data(skh1 + 671);
    const auto *skh1_672 = buffer.data(skh1 + 672);
    const auto *skh1_675 = buffer.data(skh1 + 675);

    const auto *ski_805 = buffer.data(ski + 805);
    const auto *ski_807 = buffer.data(ski + 807);
    const auto *ski_808 = buffer.data(ski + 808);
    const auto *ski_809 = buffer.data(ski + 809);
    const auto *ski_810 = buffer.data(ski + 810);
    const auto *ski_811 = buffer.data(ski + 811);
    const auto *ski_812 = buffer.data(ski + 812);
    const auto *ski_814 = buffer.data(ski + 814);
    const auto *ski_815 = buffer.data(ski + 815);
    const auto *ski_817 = buffer.data(ski + 817);
    const auto *ski_818 = buffer.data(ski + 818);
    const auto *ski_821 = buffer.data(ski + 821);
    const auto *ski_822 = buffer.data(ski + 822);
    const auto *ski_824 = buffer.data(ski + 824);
    const auto *ski_826 = buffer.data(ski + 826);
    const auto *ski_829 = buffer.data(ski + 829);
    const auto *ski_830 = buffer.data(ski + 830);
    const auto *ski_832 = buffer.data(ski + 832);
    const auto *ski_833 = buffer.data(ski + 833);
    const auto *ski_834 = buffer.data(ski + 834);
    const auto *ski_835 = buffer.data(ski + 835);
    const auto *ski_836 = buffer.data(ski + 836);
    const auto *ski_837 = buffer.data(ski + 837);
    const auto *ski_838 = buffer.data(ski + 838);
    const auto *ski_839 = buffer.data(ski + 839);
    const auto *ski_840 = buffer.data(ski + 840);
    const auto *ski_842 = buffer.data(ski + 842);
    const auto *ski_843 = buffer.data(ski + 843);
    const auto *ski_845 = buffer.data(ski + 845);
    const auto *ski_846 = buffer.data(ski + 846);
    const auto *ski_849 = buffer.data(ski + 849);
    const auto *ski_850 = buffer.data(ski + 850);
    const auto *ski_852 = buffer.data(ski + 852);
    const auto *ski_854 = buffer.data(ski + 854);
    const auto *ski_855 = buffer.data(ski + 855);
    const auto *ski_857 = buffer.data(ski + 857);
    const auto *ski_858 = buffer.data(ski + 858);
    const auto *ski_860 = buffer.data(ski + 860);
    const auto *ski_861 = buffer.data(ski + 861);
    const auto *ski_862 = buffer.data(ski + 862);
    const auto *ski_863 = buffer.data(ski + 863);
    const auto *ski_864 = buffer.data(ski + 864);
    const auto *ski_865 = buffer.data(ski + 865);
    const auto *ski_866 = buffer.data(ski + 866);
    const auto *ski_867 = buffer.data(ski + 867);
    const auto *ski_868 = buffer.data(ski + 868);
    const auto *ski_870 = buffer.data(ski + 870);
    const auto *ski_871 = buffer.data(ski + 871);
    const auto *ski_873 = buffer.data(ski + 873);
    const auto *ski_874 = buffer.data(ski + 874);
    const auto *ski_877 = buffer.data(ski + 877);
    const auto *ski_878 = buffer.data(ski + 878);
    const auto *ski_880 = buffer.data(ski + 880);
    const auto *ski_882 = buffer.data(ski + 882);
    const auto *ski_883 = buffer.data(ski + 883);
    const auto *ski_885 = buffer.data(ski + 885);
    const auto *ski_886 = buffer.data(ski + 886);
    const auto *ski_888 = buffer.data(ski + 888);
    const auto *ski_889 = buffer.data(ski + 889);
    const auto *ski_890 = buffer.data(ski + 890);
    const auto *ski_891 = buffer.data(ski + 891);
    const auto *ski_892 = buffer.data(ski + 892);
    const auto *ski_893 = buffer.data(ski + 893);
    const auto *ski_894 = buffer.data(ski + 894);
    const auto *ski_895 = buffer.data(ski + 895);
    const auto *ski_896 = buffer.data(ski + 896);
    const auto *ski_898 = buffer.data(ski + 898);
    const auto *ski_899 = buffer.data(ski + 899);

#pragma omp simd aligned(t_1037, t_1038, t_1039, pc_y, pc_z, sii_611, sii_612, skh0_605, \
                         skh0_606, skh1_605, skh1_606, ski_805, ski_807, \
                         ski_808 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1037[k] = f_3 * pc_z[k] * ski_805[k];

        t_1038[k] = f_0 * sii_611[k]
                    + f_4 * skh0_605[k]
                    - f_5 * skh1_605[k]
                    + f_3 * pc_y[k] * ski_807[k];

        t_1039[k] = f_0 * sii_612[k]
                    + f_6 * skh0_606[k]
                    - f_7 * skh1_606[k]
                    + f_3 * pc_y[k] * ski_808[k];
    }

#pragma omp simd aligned(t_1040, t_1041, t_1042, t_1043, pc_y, pc_z, sii_613, sii_614, \
                         sii_615, skh0_607, skh0_608, skh1_607, skh1_608, ski_809, ski_810, \
                         ski_811 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1040[k] = f_0 * sii_613[k]
                    + f_8 * skh0_607[k]
                    - f_9 * skh1_607[k]
                    + f_3 * pc_y[k] * ski_809[k];

        t_1041[k] = f_0 * sii_614[k]
                    + f_10 * skh0_608[k]
                    - f_11 * skh1_608[k]
                    + f_3 * pc_y[k] * ski_810[k];

        t_1042[k] = f_0 * sii_615[k]
                    + f_3 * pc_y[k] * ski_811[k];

        t_1043[k] = f_1 * skh0_608[k]
                    - f_2 * skh1_608[k]
                    + f_3 * pc_z[k] * ski_811[k];
    }

#pragma omp simd aligned(t_1044, t_1045, t_1046, t_1047, pb_z, pc_y, pc_z, sik0_756, sik0_759, \
                         sii_588, sii_616, sik1_756, sik1_759, \
                         ski_812 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1044[k] = pb_z[k] * sik0_756[k]
                    - f_12 * pc_z[k] * sik1_756[k];

        t_1045[k] = f_18 * sii_616[k]
                    + f_3 * pc_y[k] * ski_812[k];

        t_1046[k] = f_13 * sii_588[k]
                    + f_3 * pc_z[k] * ski_812[k];

        t_1047[k] = pb_z[k] * sik0_759[k]
                    - f_12 * pc_z[k] * sik1_759[k];
    }

#pragma omp simd aligned(t_1048, t_1049, t_1050, pb_z, pc_x, pc_y, pc_z, sik0_762, sii_618, \
                         sik1_762, skh0_614, skh1_614, ski_814, \
                         ski_817 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1048[k] = f_18 * sii_618[k]
                    + f_3 * pc_y[k] * ski_814[k];

        t_1049[k] = f_4 * skh0_614[k]
                    - f_5 * skh1_614[k]
                    + f_3 * pc_x[k] * ski_817[k];

        t_1050[k] = pb_z[k] * sik0_762[k]
                    - f_12 * pc_z[k] * sik1_762[k];
    }

#pragma omp simd aligned(t_1051, t_1052, t_1053, pc_x, pc_y, pc_z, sii_591, sii_621, skh0_618, \
                         skh1_618, ski_815, ski_817, ski_821 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1051[k] = f_13 * sii_591[k]
                    + f_3 * pc_z[k] * ski_815[k];

        t_1052[k] = f_18 * sii_621[k]
                    + f_3 * pc_y[k] * ski_817[k];

        t_1053[k] = f_6 * skh0_618[k]
                    - f_7 * skh1_618[k]
                    + f_3 * pc_x[k] * ski_821[k];
    }

#pragma omp simd aligned(t_1054, t_1055, t_1056, pb_z, pc_x, pc_z, sik0_766, sii_594, \
                         sik1_766, skh0_621, skh1_621, ski_818, \
                         ski_824 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1054[k] = pb_z[k] * sik0_766[k]
                    - f_12 * pc_z[k] * sik1_766[k];

        t_1055[k] = f_13 * sii_594[k]
                    + f_3 * pc_z[k] * ski_818[k];

        t_1056[k] = f_8 * skh0_621[k]
                    - f_9 * skh1_621[k]
                    + f_3 * pc_x[k] * ski_824[k];
    }

#pragma omp simd aligned(t_1057, t_1058, t_1059, pb_z, pc_x, pc_y, pc_z, sik0_771, sii_625, \
                         sik1_771, skh0_623, skh1_623, ski_821, \
                         ski_826 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1057[k] = f_18 * sii_625[k]
                    + f_3 * pc_y[k] * ski_821[k];

        t_1058[k] = f_8 * skh0_623[k]
                    - f_9 * skh1_623[k]
                    + f_3 * pc_x[k] * ski_826[k];

        t_1059[k] = pb_z[k] * sik0_771[k]
                    - f_12 * pc_z[k] * sik1_771[k];
    }

#pragma omp simd aligned(t_1060, t_1061, t_1062, pc_x, pc_z, sii_598, skh0_626, skh0_627, \
                         skh1_626, skh1_627, ski_822, ski_829, \
                         ski_830 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1060[k] = f_13 * sii_598[k]
                    + f_3 * pc_z[k] * ski_822[k];

        t_1061[k] = f_10 * skh0_626[k]
                    - f_11 * skh1_626[k]
                    + f_3 * pc_x[k] * ski_829[k];

        t_1062[k] = f_10 * skh0_627[k]
                    - f_11 * skh1_627[k]
                    + f_3 * pc_x[k] * ski_830[k];
    }

#pragma omp simd aligned(t_1063, t_1064, t_1065, t_1066, t_1067, pc_x, pc_y, sii_630, \
                         skh0_629, skh1_629, ski_826, ski_832, ski_833, ski_834, \
                         ski_835 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1063[k] = f_18 * sii_630[k]
                    + f_3 * pc_y[k] * ski_826[k];

        t_1064[k] = f_10 * skh0_629[k]
                    - f_11 * skh1_629[k]
                    + f_3 * pc_x[k] * ski_832[k];

        t_1065[k] = f_3 * pc_x[k] * ski_833[k];

        t_1066[k] = f_3 * pc_x[k] * ski_834[k];

        t_1067[k] = f_3 * pc_x[k] * ski_835[k];
    }

#pragma omp simd aligned(t_1068, t_1069, t_1070, t_1071, t_1072, pb_z, pc_x, pc_z, sik0_784, \
                         sik1_784, ski_836, ski_837, ski_838, ski_839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1068[k] = f_3 * pc_x[k] * ski_836[k];

        t_1069[k] = f_3 * pc_x[k] * ski_837[k];

        t_1070[k] = f_3 * pc_x[k] * ski_838[k];

        t_1071[k] = f_3 * pc_x[k] * ski_839[k];

        t_1072[k] = pb_z[k] * sik0_784[k]
                    - f_12 * pc_z[k] * sik1_784[k];
    }

#pragma omp simd aligned(t_1073, t_1074, t_1075, pb_z, pc_z, sik0_786, sik0_787, sii_609, \
                         sii_610, sii_611, sik1_786, sik1_787, \
                         ski_833 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1073[k] = f_13 * sii_609[k]
                    + f_3 * pc_z[k] * ski_833[k];

        t_1074[k] = pb_z[k] * sik0_786[k]
                    + f_14 * sii_610[k]
                    - f_12 * pc_z[k] * sik1_786[k];

        t_1075[k] = pb_z[k] * sik0_787[k]
                    + f_15 * sii_611[k]
                    - f_12 * pc_z[k] * sik1_787[k];
    }

#pragma omp simd aligned(t_1076, t_1077, t_1078, pb_z, pc_y, pc_z, sik0_788, sik0_789, \
                         sii_612, sii_613, sii_643, sik1_788, sik1_789, \
                         ski_839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1076[k] = pb_z[k] * sik0_788[k]
                    + f_16 * sii_612[k]
                    - f_12 * pc_z[k] * sik1_788[k];

        t_1077[k] = pb_z[k] * sik0_789[k]
                    + f_17 * sii_613[k]
                    - f_12 * pc_z[k] * sik1_789[k];

        t_1078[k] = f_18 * sii_643[k]
                    + f_3 * pc_y[k] * ski_839[k];
    }

#pragma omp simd aligned(t_1079, t_1080, t_1081, t_1082, pc_x, pc_y, pc_z, sii_615, sii_616, \
                         sii_644, skh0_629, skh0_630, skh1_629, skh1_630, ski_839, \
                         ski_840 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1079[k] = f_13 * sii_615[k]
                    + f_1 * skh0_629[k]
                    - f_2 * skh1_629[k]
                    + f_3 * pc_z[k] * ski_839[k];

        t_1080[k] = f_1 * skh0_630[k]
                    - f_2 * skh1_630[k]
                    + f_3 * pc_x[k] * ski_840[k];

        t_1081[k] = f_17 * sii_644[k]
                    + f_3 * pc_y[k] * ski_840[k];

        t_1082[k] = f_14 * sii_616[k]
                    + f_3 * pc_z[k] * ski_840[k];
    }

#pragma omp simd aligned(t_1083, t_1084, t_1085, pc_x, pc_y, sii_646, skh0_633, skh0_635, \
                         skh1_633, skh1_635, ski_842, ski_843, \
                         ski_845 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1083[k] = f_4 * skh0_633[k]
                    - f_5 * skh1_633[k]
                    + f_3 * pc_x[k] * ski_843[k];

        t_1084[k] = f_17 * sii_646[k]
                    + f_3 * pc_y[k] * ski_842[k];

        t_1085[k] = f_4 * skh0_635[k]
                    - f_5 * skh1_635[k]
                    + f_3 * pc_x[k] * ski_845[k];
    }

#pragma omp simd aligned(t_1086, t_1087, t_1088, pc_x, pc_y, pc_z, sii_619, sii_649, skh0_636, \
                         skh1_636, ski_843, ski_845, ski_846 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1086[k] = f_6 * skh0_636[k]
                    - f_7 * skh1_636[k]
                    + f_3 * pc_x[k] * ski_846[k];

        t_1087[k] = f_14 * sii_619[k]
                    + f_3 * pc_z[k] * ski_843[k];

        t_1088[k] = f_17 * sii_649[k]
                    + f_3 * pc_y[k] * ski_845[k];
    }

#pragma omp simd aligned(t_1089, t_1090, t_1091, pc_x, pc_z, sii_622, skh0_639, skh0_640, \
                         skh1_639, skh1_640, ski_846, ski_849, \
                         ski_850 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1089[k] = f_6 * skh0_639[k]
                    - f_7 * skh1_639[k]
                    + f_3 * pc_x[k] * ski_849[k];

        t_1090[k] = f_8 * skh0_640[k]
                    - f_9 * skh1_640[k]
                    + f_3 * pc_x[k] * ski_850[k];

        t_1091[k] = f_14 * sii_622[k]
                    + f_3 * pc_z[k] * ski_846[k];
    }

#pragma omp simd aligned(t_1092, t_1093, t_1094, pc_x, pc_y, sii_653, skh0_642, skh0_644, \
                         skh1_642, skh1_644, ski_849, ski_852, \
                         ski_854 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1092[k] = f_8 * skh0_642[k]
                    - f_9 * skh1_642[k]
                    + f_3 * pc_x[k] * ski_852[k];

        t_1093[k] = f_17 * sii_653[k]
                    + f_3 * pc_y[k] * ski_849[k];

        t_1094[k] = f_8 * skh0_644[k]
                    - f_9 * skh1_644[k]
                    + f_3 * pc_x[k] * ski_854[k];
    }

#pragma omp simd aligned(t_1095, t_1096, t_1097, pc_x, pc_z, sii_626, skh0_645, skh0_647, \
                         skh1_645, skh1_647, ski_850, ski_855, \
                         ski_857 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1095[k] = f_10 * skh0_645[k]
                    - f_11 * skh1_645[k]
                    + f_3 * pc_x[k] * ski_855[k];

        t_1096[k] = f_14 * sii_626[k]
                    + f_3 * pc_z[k] * ski_850[k];

        t_1097[k] = f_10 * skh0_647[k]
                    - f_11 * skh1_647[k]
                    + f_3 * pc_x[k] * ski_857[k];
    }

#pragma omp simd aligned(t_1098, t_1099, t_1100, t_1101, pc_x, pc_y, sii_658, skh0_648, \
                         skh0_650, skh1_648, skh1_650, ski_854, ski_858, ski_860, \
                         ski_861 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1098[k] = f_10 * skh0_648[k]
                    - f_11 * skh1_648[k]
                    + f_3 * pc_x[k] * ski_858[k];

        t_1099[k] = f_17 * sii_658[k]
                    + f_3 * pc_y[k] * ski_854[k];

        t_1100[k] = f_10 * skh0_650[k]
                    - f_11 * skh1_650[k]
                    + f_3 * pc_x[k] * ski_860[k];

        t_1101[k] = f_3 * pc_x[k] * ski_861[k];
    }

#pragma omp simd aligned(t_1102, t_1103, t_1104, t_1105, t_1106, t_1107, pc_x, ski_862, \
                         ski_863, ski_864, ski_865, ski_866, ski_867 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1102[k] = f_3 * pc_x[k] * ski_862[k];

        t_1103[k] = f_3 * pc_x[k] * ski_863[k];

        t_1104[k] = f_3 * pc_x[k] * ski_864[k];

        t_1105[k] = f_3 * pc_x[k] * ski_865[k];

        t_1106[k] = f_3 * pc_x[k] * ski_866[k];

        t_1107[k] = f_3 * pc_x[k] * ski_867[k];
    }

#pragma omp simd aligned(t_1108, t_1109, t_1110, pc_y, pc_z, sii_637, sii_665, sii_667, \
                         skh0_645, skh0_647, skh1_645, skh1_647, ski_861, \
                         ski_863 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1108[k] = f_17 * sii_665[k]
                    + f_1 * skh0_645[k]
                    - f_2 * skh1_645[k]
                    + f_3 * pc_y[k] * ski_861[k];

        t_1109[k] = f_14 * sii_637[k]
                    + f_3 * pc_z[k] * ski_861[k];

        t_1110[k] = f_17 * sii_667[k]
                    + f_4 * skh0_647[k]
                    - f_5 * skh1_647[k]
                    + f_3 * pc_y[k] * ski_863[k];
    }

#pragma omp simd aligned(t_1111, t_1112, t_1113, pc_y, sii_668, sii_669, sii_670, skh0_648, \
                         skh0_649, skh0_650, skh1_648, skh1_649, skh1_650, ski_864, ski_865, \
                         ski_866 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1111[k] = f_17 * sii_668[k]
                    + f_6 * skh0_648[k]
                    - f_7 * skh1_648[k]
                    + f_3 * pc_y[k] * ski_864[k];

        t_1112[k] = f_17 * sii_669[k]
                    + f_8 * skh0_649[k]
                    - f_9 * skh1_649[k]
                    + f_3 * pc_y[k] * ski_865[k];

        t_1113[k] = f_17 * sii_670[k]
                    + f_10 * skh0_650[k]
                    - f_11 * skh1_650[k]
                    + f_3 * pc_y[k] * ski_866[k];
    }

#pragma omp simd aligned(t_1114, t_1115, t_1116, t_1117, pc_x, pc_y, pc_z, sii_643, sii_671, \
                         sii_672, skh0_650, skh0_651, skh1_650, skh1_651, ski_867, \
                         ski_868 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1114[k] = f_17 * sii_671[k]
                    + f_3 * pc_y[k] * ski_867[k];

        t_1115[k] = f_14 * sii_643[k]
                    + f_1 * skh0_650[k]
                    - f_2 * skh1_650[k]
                    + f_3 * pc_z[k] * ski_867[k];

        t_1116[k] = f_1 * skh0_651[k]
                    - f_2 * skh1_651[k]
                    + f_3 * pc_x[k] * ski_868[k];

        t_1117[k] = f_16 * sii_672[k]
                    + f_3 * pc_y[k] * ski_868[k];
    }

#pragma omp simd aligned(t_1118, t_1119, t_1120, pc_x, pc_y, pc_z, sii_644, sii_674, skh0_654, \
                         skh1_654, ski_868, ski_870, ski_871 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1118[k] = f_15 * sii_644[k]
                    + f_3 * pc_z[k] * ski_868[k];

        t_1119[k] = f_4 * skh0_654[k]
                    - f_5 * skh1_654[k]
                    + f_3 * pc_x[k] * ski_871[k];

        t_1120[k] = f_16 * sii_674[k]
                    + f_3 * pc_y[k] * ski_870[k];
    }

#pragma omp simd aligned(t_1121, t_1122, t_1123, t_1124, pc_x, pc_y, pc_z, sii_647, sii_677, \
                         skh0_656, skh0_657, skh1_656, skh1_657, ski_871, ski_873, \
                         ski_874 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1121[k] = f_4 * skh0_656[k]
                    - f_5 * skh1_656[k]
                    + f_3 * pc_x[k] * ski_873[k];

        t_1122[k] = f_6 * skh0_657[k]
                    - f_7 * skh1_657[k]
                    + f_3 * pc_x[k] * ski_874[k];

        t_1123[k] = f_15 * sii_647[k]
                    + f_3 * pc_z[k] * ski_871[k];

        t_1124[k] = f_16 * sii_677[k]
                    + f_3 * pc_y[k] * ski_873[k];
    }

#pragma omp simd aligned(t_1125, t_1126, t_1127, pc_x, pc_z, sii_650, skh0_660, skh0_661, \
                         skh1_660, skh1_661, ski_874, ski_877, \
                         ski_878 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1125[k] = f_6 * skh0_660[k]
                    - f_7 * skh1_660[k]
                    + f_3 * pc_x[k] * ski_877[k];

        t_1126[k] = f_8 * skh0_661[k]
                    - f_9 * skh1_661[k]
                    + f_3 * pc_x[k] * ski_878[k];

        t_1127[k] = f_15 * sii_650[k]
                    + f_3 * pc_z[k] * ski_874[k];
    }

#pragma omp simd aligned(t_1128, t_1129, t_1130, pc_x, pc_y, sii_681, skh0_663, skh0_665, \
                         skh1_663, skh1_665, ski_877, ski_880, \
                         ski_882 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1128[k] = f_8 * skh0_663[k]
                    - f_9 * skh1_663[k]
                    + f_3 * pc_x[k] * ski_880[k];

        t_1129[k] = f_16 * sii_681[k]
                    + f_3 * pc_y[k] * ski_877[k];

        t_1130[k] = f_8 * skh0_665[k]
                    - f_9 * skh1_665[k]
                    + f_3 * pc_x[k] * ski_882[k];
    }

#pragma omp simd aligned(t_1131, t_1132, t_1133, pc_x, pc_z, sii_654, skh0_666, skh0_668, \
                         skh1_666, skh1_668, ski_878, ski_883, \
                         ski_885 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1131[k] = f_10 * skh0_666[k]
                    - f_11 * skh1_666[k]
                    + f_3 * pc_x[k] * ski_883[k];

        t_1132[k] = f_15 * sii_654[k]
                    + f_3 * pc_z[k] * ski_878[k];

        t_1133[k] = f_10 * skh0_668[k]
                    - f_11 * skh1_668[k]
                    + f_3 * pc_x[k] * ski_885[k];
    }

#pragma omp simd aligned(t_1134, t_1135, t_1136, t_1137, pc_x, pc_y, sii_686, skh0_669, \
                         skh0_671, skh1_669, skh1_671, ski_882, ski_886, ski_888, \
                         ski_889 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1134[k] = f_10 * skh0_669[k]
                    - f_11 * skh1_669[k]
                    + f_3 * pc_x[k] * ski_886[k];

        t_1135[k] = f_16 * sii_686[k]
                    + f_3 * pc_y[k] * ski_882[k];

        t_1136[k] = f_10 * skh0_671[k]
                    - f_11 * skh1_671[k]
                    + f_3 * pc_x[k] * ski_888[k];

        t_1137[k] = f_3 * pc_x[k] * ski_889[k];
    }

#pragma omp simd aligned(t_1138, t_1139, t_1140, t_1141, t_1142, t_1143, pc_x, ski_890, \
                         ski_891, ski_892, ski_893, ski_894, ski_895 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1138[k] = f_3 * pc_x[k] * ski_890[k];

        t_1139[k] = f_3 * pc_x[k] * ski_891[k];

        t_1140[k] = f_3 * pc_x[k] * ski_892[k];

        t_1141[k] = f_3 * pc_x[k] * ski_893[k];

        t_1142[k] = f_3 * pc_x[k] * ski_894[k];

        t_1143[k] = f_3 * pc_x[k] * ski_895[k];
    }

#pragma omp simd aligned(t_1144, t_1145, t_1146, pc_y, pc_z, sii_665, sii_693, sii_695, \
                         skh0_666, skh0_668, skh1_666, skh1_668, ski_889, \
                         ski_891 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1144[k] = f_16 * sii_693[k]
                    + f_1 * skh0_666[k]
                    - f_2 * skh1_666[k]
                    + f_3 * pc_y[k] * ski_889[k];

        t_1145[k] = f_15 * sii_665[k]
                    + f_3 * pc_z[k] * ski_889[k];

        t_1146[k] = f_16 * sii_695[k]
                    + f_4 * skh0_668[k]
                    - f_5 * skh1_668[k]
                    + f_3 * pc_y[k] * ski_891[k];
    }

#pragma omp simd aligned(t_1147, t_1148, t_1149, pc_y, sii_696, sii_697, sii_698, skh0_669, \
                         skh0_670, skh0_671, skh1_669, skh1_670, skh1_671, ski_892, ski_893, \
                         ski_894 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1147[k] = f_16 * sii_696[k]
                    + f_6 * skh0_669[k]
                    - f_7 * skh1_669[k]
                    + f_3 * pc_y[k] * ski_892[k];

        t_1148[k] = f_16 * sii_697[k]
                    + f_8 * skh0_670[k]
                    - f_9 * skh1_670[k]
                    + f_3 * pc_y[k] * ski_893[k];

        t_1149[k] = f_16 * sii_698[k]
                    + f_10 * skh0_671[k]
                    - f_11 * skh1_671[k]
                    + f_3 * pc_y[k] * ski_894[k];
    }

#pragma omp simd aligned(t_1150, t_1151, t_1152, t_1153, pc_x, pc_y, pc_z, sii_671, sii_699, \
                         sii_700, skh0_671, skh0_672, skh1_671, skh1_672, ski_895, \
                         ski_896 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1150[k] = f_16 * sii_699[k]
                    + f_3 * pc_y[k] * ski_895[k];

        t_1151[k] = f_15 * sii_671[k]
                    + f_1 * skh0_671[k]
                    - f_2 * skh1_671[k]
                    + f_3 * pc_z[k] * ski_895[k];

        t_1152[k] = f_1 * skh0_672[k]
                    - f_2 * skh1_672[k]
                    + f_3 * pc_x[k] * ski_896[k];

        t_1153[k] = f_15 * sii_700[k]
                    + f_3 * pc_y[k] * ski_896[k];
    }

#pragma omp simd aligned(t_1154, t_1155, t_1156, pc_x, pc_y, pc_z, sii_672, sii_702, skh0_675, \
                         skh1_675, ski_896, ski_898, ski_899 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1154[k] = f_16 * sii_672[k]
                    + f_3 * pc_z[k] * ski_896[k];

        t_1155[k] = f_4 * skh0_675[k]
                    - f_5 * skh1_675[k]
                    + f_3 * pc_x[k] * ski_899[k];

        t_1156[k] = f_15 * sii_702[k]
                    + f_3 * pc_y[k] * ski_898[k];
    }
}

static auto
compute_prim_skk_three_center_electron_repulsion_0_piece10(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t sik0,
                                                           const size_t sii, const size_t sik1,
                                                           const size_t skh0, const size_t skh1,
                                                           const size_t ski, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / q;
    const auto f_1 = 3.0 / gamma;
    const auto f_2 = 3.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.0 / gamma;
    const auto f_5 = 2.0 * p / (gamma * q);
    const auto f_6 = 1.5 / gamma;
    const auto f_7 = 1.5 * p / (gamma * q);
    const auto f_8 = 1.0 / gamma;
    const auto f_9 = p / (gamma * q);
    const auto f_10 = 0.5 / gamma;
    const auto f_11 = 0.5 * p / (gamma * q);
    const auto f_12 = gamma / q;
    const auto f_13 = 0.5 / q;
    const auto f_14 = 1.0 / q;
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;
    const auto f_18 = 3.0 / q;

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

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sik0_972 = buffer.data(sik0 + 972);
    const auto *sik0_977 = buffer.data(sik0 + 977);
    const auto *sik0_981 = buffer.data(sik0 + 981);
    const auto *sik0_986 = buffer.data(sik0 + 986);
    const auto *sik0_992 = buffer.data(sik0 + 992);
    const auto *sik0_1000 = buffer.data(sik0 + 1000);
    const auto *sik0_1002 = buffer.data(sik0 + 1002);
    const auto *sik0_1003 = buffer.data(sik0 + 1003);
    const auto *sik0_1004 = buffer.data(sik0 + 1004);
    const auto *sik0_1005 = buffer.data(sik0 + 1005);
    const auto *sik0_1007 = buffer.data(sik0 + 1007);

    const auto *sii_675 = buffer.data(sii + 675);
    const auto *sii_678 = buffer.data(sii + 678);
    const auto *sii_682 = buffer.data(sii + 682);
    const auto *sii_693 = buffer.data(sii + 693);
    const auto *sii_699 = buffer.data(sii + 699);
    const auto *sii_700 = buffer.data(sii + 700);
    const auto *sii_703 = buffer.data(sii + 703);
    const auto *sii_705 = buffer.data(sii + 705);
    const auto *sii_706 = buffer.data(sii + 706);
    const auto *sii_709 = buffer.data(sii + 709);
    const auto *sii_710 = buffer.data(sii + 710);
    const auto *sii_714 = buffer.data(sii + 714);
    const auto *sii_721 = buffer.data(sii + 721);
    const auto *sii_723 = buffer.data(sii + 723);
    const auto *sii_724 = buffer.data(sii + 724);
    const auto *sii_725 = buffer.data(sii + 725);
    const auto *sii_726 = buffer.data(sii + 726);
    const auto *sii_727 = buffer.data(sii + 727);
    const auto *sii_728 = buffer.data(sii + 728);
    const auto *sii_730 = buffer.data(sii + 730);
    const auto *sii_731 = buffer.data(sii + 731);
    const auto *sii_733 = buffer.data(sii + 733);
    const auto *sii_734 = buffer.data(sii + 734);
    const auto *sii_737 = buffer.data(sii + 737);
    const auto *sii_738 = buffer.data(sii + 738);
    const auto *sii_742 = buffer.data(sii + 742);
    const auto *sii_749 = buffer.data(sii + 749);
    const auto *sii_751 = buffer.data(sii + 751);
    const auto *sii_752 = buffer.data(sii + 752);
    const auto *sii_753 = buffer.data(sii + 753);
    const auto *sii_754 = buffer.data(sii + 754);
    const auto *sii_755 = buffer.data(sii + 755);
    const auto *sii_756 = buffer.data(sii + 756);
    const auto *sii_758 = buffer.data(sii + 758);
    const auto *sii_759 = buffer.data(sii + 759);
    const auto *sii_761 = buffer.data(sii + 761);
    const auto *sii_762 = buffer.data(sii + 762);
    const auto *sii_765 = buffer.data(sii + 765);
    const auto *sii_766 = buffer.data(sii + 766);
    const auto *sii_770 = buffer.data(sii + 770);
    const auto *sii_777 = buffer.data(sii + 777);
    const auto *sii_779 = buffer.data(sii + 779);
    const auto *sii_780 = buffer.data(sii + 780);
    const auto *sii_781 = buffer.data(sii + 781);
    const auto *sii_782 = buffer.data(sii + 782);
    const auto *sii_783 = buffer.data(sii + 783);

    const auto *sik1_972 = buffer.data(sik1 + 972);
    const auto *sik1_977 = buffer.data(sik1 + 977);
    const auto *sik1_981 = buffer.data(sik1 + 981);
    const auto *sik1_986 = buffer.data(sik1 + 986);
    const auto *sik1_992 = buffer.data(sik1 + 992);
    const auto *sik1_1000 = buffer.data(sik1 + 1000);
    const auto *sik1_1002 = buffer.data(sik1 + 1002);
    const auto *sik1_1003 = buffer.data(sik1 + 1003);
    const auto *sik1_1004 = buffer.data(sik1 + 1004);
    const auto *sik1_1005 = buffer.data(sik1 + 1005);
    const auto *sik1_1007 = buffer.data(sik1 + 1007);

    const auto *skh0_677 = buffer.data(skh0 + 677);
    const auto *skh0_678 = buffer.data(skh0 + 678);
    const auto *skh0_681 = buffer.data(skh0 + 681);
    const auto *skh0_682 = buffer.data(skh0 + 682);
    const auto *skh0_684 = buffer.data(skh0 + 684);
    const auto *skh0_686 = buffer.data(skh0 + 686);
    const auto *skh0_687 = buffer.data(skh0 + 687);
    const auto *skh0_689 = buffer.data(skh0 + 689);
    const auto *skh0_690 = buffer.data(skh0 + 690);
    const auto *skh0_691 = buffer.data(skh0 + 691);
    const auto *skh0_692 = buffer.data(skh0 + 692);
    const auto *skh0_693 = buffer.data(skh0 + 693);
    const auto *skh0_696 = buffer.data(skh0 + 696);
    const auto *skh0_698 = buffer.data(skh0 + 698);
    const auto *skh0_699 = buffer.data(skh0 + 699);
    const auto *skh0_702 = buffer.data(skh0 + 702);
    const auto *skh0_703 = buffer.data(skh0 + 703);
    const auto *skh0_705 = buffer.data(skh0 + 705);
    const auto *skh0_707 = buffer.data(skh0 + 707);
    const auto *skh0_708 = buffer.data(skh0 + 708);
    const auto *skh0_710 = buffer.data(skh0 + 710);
    const auto *skh0_711 = buffer.data(skh0 + 711);
    const auto *skh0_712 = buffer.data(skh0 + 712);
    const auto *skh0_713 = buffer.data(skh0 + 713);
    const auto *skh0_717 = buffer.data(skh0 + 717);
    const auto *skh0_720 = buffer.data(skh0 + 720);
    const auto *skh0_724 = buffer.data(skh0 + 724);
    const auto *skh0_726 = buffer.data(skh0 + 726);
    const auto *skh0_729 = buffer.data(skh0 + 729);
    const auto *skh0_731 = buffer.data(skh0 + 731);
    const auto *skh0_732 = buffer.data(skh0 + 732);
    const auto *skh0_735 = buffer.data(skh0 + 735);
    const auto *skh0_738 = buffer.data(skh0 + 738);
    const auto *skh0_740 = buffer.data(skh0 + 740);
    const auto *skh0_741 = buffer.data(skh0 + 741);
    const auto *skh0_744 = buffer.data(skh0 + 744);
    const auto *skh0_745 = buffer.data(skh0 + 745);
    const auto *skh0_747 = buffer.data(skh0 + 747);
    const auto *skh0_749 = buffer.data(skh0 + 749);
    const auto *skh0_750 = buffer.data(skh0 + 750);
    const auto *skh0_752 = buffer.data(skh0 + 752);
    const auto *skh0_753 = buffer.data(skh0 + 753);

    const auto *skh1_677 = buffer.data(skh1 + 677);
    const auto *skh1_678 = buffer.data(skh1 + 678);
    const auto *skh1_681 = buffer.data(skh1 + 681);
    const auto *skh1_682 = buffer.data(skh1 + 682);
    const auto *skh1_684 = buffer.data(skh1 + 684);
    const auto *skh1_686 = buffer.data(skh1 + 686);
    const auto *skh1_687 = buffer.data(skh1 + 687);
    const auto *skh1_689 = buffer.data(skh1 + 689);
    const auto *skh1_690 = buffer.data(skh1 + 690);
    const auto *skh1_691 = buffer.data(skh1 + 691);
    const auto *skh1_692 = buffer.data(skh1 + 692);
    const auto *skh1_693 = buffer.data(skh1 + 693);
    const auto *skh1_696 = buffer.data(skh1 + 696);
    const auto *skh1_698 = buffer.data(skh1 + 698);
    const auto *skh1_699 = buffer.data(skh1 + 699);
    const auto *skh1_702 = buffer.data(skh1 + 702);
    const auto *skh1_703 = buffer.data(skh1 + 703);
    const auto *skh1_705 = buffer.data(skh1 + 705);
    const auto *skh1_707 = buffer.data(skh1 + 707);
    const auto *skh1_708 = buffer.data(skh1 + 708);
    const auto *skh1_710 = buffer.data(skh1 + 710);
    const auto *skh1_711 = buffer.data(skh1 + 711);
    const auto *skh1_712 = buffer.data(skh1 + 712);
    const auto *skh1_713 = buffer.data(skh1 + 713);
    const auto *skh1_717 = buffer.data(skh1 + 717);
    const auto *skh1_720 = buffer.data(skh1 + 720);
    const auto *skh1_724 = buffer.data(skh1 + 724);
    const auto *skh1_726 = buffer.data(skh1 + 726);
    const auto *skh1_729 = buffer.data(skh1 + 729);
    const auto *skh1_731 = buffer.data(skh1 + 731);
    const auto *skh1_732 = buffer.data(skh1 + 732);
    const auto *skh1_735 = buffer.data(skh1 + 735);
    const auto *skh1_738 = buffer.data(skh1 + 738);
    const auto *skh1_740 = buffer.data(skh1 + 740);
    const auto *skh1_741 = buffer.data(skh1 + 741);
    const auto *skh1_744 = buffer.data(skh1 + 744);
    const auto *skh1_745 = buffer.data(skh1 + 745);
    const auto *skh1_747 = buffer.data(skh1 + 747);
    const auto *skh1_749 = buffer.data(skh1 + 749);
    const auto *skh1_750 = buffer.data(skh1 + 750);
    const auto *skh1_752 = buffer.data(skh1 + 752);
    const auto *skh1_753 = buffer.data(skh1 + 753);

    const auto *ski_899 = buffer.data(ski + 899);
    const auto *ski_901 = buffer.data(ski + 901);
    const auto *ski_902 = buffer.data(ski + 902);
    const auto *ski_905 = buffer.data(ski + 905);
    const auto *ski_906 = buffer.data(ski + 906);
    const auto *ski_908 = buffer.data(ski + 908);
    const auto *ski_910 = buffer.data(ski + 910);
    const auto *ski_911 = buffer.data(ski + 911);
    const auto *ski_913 = buffer.data(ski + 913);
    const auto *ski_914 = buffer.data(ski + 914);
    const auto *ski_916 = buffer.data(ski + 916);
    const auto *ski_917 = buffer.data(ski + 917);
    const auto *ski_918 = buffer.data(ski + 918);
    const auto *ski_919 = buffer.data(ski + 919);
    const auto *ski_920 = buffer.data(ski + 920);
    const auto *ski_921 = buffer.data(ski + 921);
    const auto *ski_922 = buffer.data(ski + 922);
    const auto *ski_923 = buffer.data(ski + 923);
    const auto *ski_924 = buffer.data(ski + 924);
    const auto *ski_926 = buffer.data(ski + 926);
    const auto *ski_927 = buffer.data(ski + 927);
    const auto *ski_929 = buffer.data(ski + 929);
    const auto *ski_930 = buffer.data(ski + 930);
    const auto *ski_933 = buffer.data(ski + 933);
    const auto *ski_934 = buffer.data(ski + 934);
    const auto *ski_936 = buffer.data(ski + 936);
    const auto *ski_938 = buffer.data(ski + 938);
    const auto *ski_939 = buffer.data(ski + 939);
    const auto *ski_941 = buffer.data(ski + 941);
    const auto *ski_942 = buffer.data(ski + 942);
    const auto *ski_944 = buffer.data(ski + 944);
    const auto *ski_945 = buffer.data(ski + 945);
    const auto *ski_946 = buffer.data(ski + 946);
    const auto *ski_947 = buffer.data(ski + 947);
    const auto *ski_948 = buffer.data(ski + 948);
    const auto *ski_949 = buffer.data(ski + 949);
    const auto *ski_950 = buffer.data(ski + 950);
    const auto *ski_951 = buffer.data(ski + 951);
    const auto *ski_952 = buffer.data(ski + 952);
    const auto *ski_954 = buffer.data(ski + 954);
    const auto *ski_955 = buffer.data(ski + 955);
    const auto *ski_957 = buffer.data(ski + 957);
    const auto *ski_958 = buffer.data(ski + 958);
    const auto *ski_961 = buffer.data(ski + 961);
    const auto *ski_962 = buffer.data(ski + 962);
    const auto *ski_964 = buffer.data(ski + 964);
    const auto *ski_966 = buffer.data(ski + 966);
    const auto *ski_967 = buffer.data(ski + 967);
    const auto *ski_969 = buffer.data(ski + 969);
    const auto *ski_970 = buffer.data(ski + 970);
    const auto *ski_973 = buffer.data(ski + 973);
    const auto *ski_974 = buffer.data(ski + 974);
    const auto *ski_975 = buffer.data(ski + 975);
    const auto *ski_976 = buffer.data(ski + 976);
    const auto *ski_977 = buffer.data(ski + 977);
    const auto *ski_978 = buffer.data(ski + 978);
    const auto *ski_979 = buffer.data(ski + 979);
    const auto *ski_980 = buffer.data(ski + 980);
    const auto *ski_982 = buffer.data(ski + 982);
    const auto *ski_983 = buffer.data(ski + 983);
    const auto *ski_985 = buffer.data(ski + 985);
    const auto *ski_986 = buffer.data(ski + 986);
    const auto *ski_989 = buffer.data(ski + 989);
    const auto *ski_990 = buffer.data(ski + 990);
    const auto *ski_992 = buffer.data(ski + 992);
    const auto *ski_994 = buffer.data(ski + 994);
    const auto *ski_995 = buffer.data(ski + 995);
    const auto *ski_997 = buffer.data(ski + 997);
    const auto *ski_998 = buffer.data(ski + 998);

#pragma omp simd aligned(t_1157, t_1158, t_1159, t_1160, pc_x, pc_y, pc_z, sii_675, sii_705, \
                         skh0_677, skh0_678, skh1_677, skh1_678, ski_899, ski_901, \
                         ski_902 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1157[k] = f_4 * skh0_677[k]
                    - f_5 * skh1_677[k]
                    + f_3 * pc_x[k] * ski_901[k];

        t_1158[k] = f_6 * skh0_678[k]
                    - f_7 * skh1_678[k]
                    + f_3 * pc_x[k] * ski_902[k];

        t_1159[k] = f_16 * sii_675[k]
                    + f_3 * pc_z[k] * ski_899[k];

        t_1160[k] = f_15 * sii_705[k]
                    + f_3 * pc_y[k] * ski_901[k];
    }

#pragma omp simd aligned(t_1161, t_1162, t_1163, pc_x, pc_z, sii_678, skh0_681, skh0_682, \
                         skh1_681, skh1_682, ski_902, ski_905, \
                         ski_906 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1161[k] = f_6 * skh0_681[k]
                    - f_7 * skh1_681[k]
                    + f_3 * pc_x[k] * ski_905[k];

        t_1162[k] = f_8 * skh0_682[k]
                    - f_9 * skh1_682[k]
                    + f_3 * pc_x[k] * ski_906[k];

        t_1163[k] = f_16 * sii_678[k]
                    + f_3 * pc_z[k] * ski_902[k];
    }

#pragma omp simd aligned(t_1164, t_1165, t_1166, pc_x, pc_y, sii_709, skh0_684, skh0_686, \
                         skh1_684, skh1_686, ski_905, ski_908, \
                         ski_910 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1164[k] = f_8 * skh0_684[k]
                    - f_9 * skh1_684[k]
                    + f_3 * pc_x[k] * ski_908[k];

        t_1165[k] = f_15 * sii_709[k]
                    + f_3 * pc_y[k] * ski_905[k];

        t_1166[k] = f_8 * skh0_686[k]
                    - f_9 * skh1_686[k]
                    + f_3 * pc_x[k] * ski_910[k];
    }

#pragma omp simd aligned(t_1167, t_1168, t_1169, pc_x, pc_z, sii_682, skh0_687, skh0_689, \
                         skh1_687, skh1_689, ski_906, ski_911, \
                         ski_913 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1167[k] = f_10 * skh0_687[k]
                    - f_11 * skh1_687[k]
                    + f_3 * pc_x[k] * ski_911[k];

        t_1168[k] = f_16 * sii_682[k]
                    + f_3 * pc_z[k] * ski_906[k];

        t_1169[k] = f_10 * skh0_689[k]
                    - f_11 * skh1_689[k]
                    + f_3 * pc_x[k] * ski_913[k];
    }

#pragma omp simd aligned(t_1170, t_1171, t_1172, t_1173, pc_x, pc_y, sii_714, skh0_690, \
                         skh0_692, skh1_690, skh1_692, ski_910, ski_914, ski_916, \
                         ski_917 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1170[k] = f_10 * skh0_690[k]
                    - f_11 * skh1_690[k]
                    + f_3 * pc_x[k] * ski_914[k];

        t_1171[k] = f_15 * sii_714[k]
                    + f_3 * pc_y[k] * ski_910[k];

        t_1172[k] = f_10 * skh0_692[k]
                    - f_11 * skh1_692[k]
                    + f_3 * pc_x[k] * ski_916[k];

        t_1173[k] = f_3 * pc_x[k] * ski_917[k];
    }

#pragma omp simd aligned(t_1174, t_1175, t_1176, t_1177, t_1178, t_1179, pc_x, ski_918, \
                         ski_919, ski_920, ski_921, ski_922, ski_923 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1174[k] = f_3 * pc_x[k] * ski_918[k];

        t_1175[k] = f_3 * pc_x[k] * ski_919[k];

        t_1176[k] = f_3 * pc_x[k] * ski_920[k];

        t_1177[k] = f_3 * pc_x[k] * ski_921[k];

        t_1178[k] = f_3 * pc_x[k] * ski_922[k];

        t_1179[k] = f_3 * pc_x[k] * ski_923[k];
    }

#pragma omp simd aligned(t_1180, t_1181, t_1182, pc_y, pc_z, sii_693, sii_721, sii_723, \
                         skh0_687, skh0_689, skh1_687, skh1_689, ski_917, \
                         ski_919 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1180[k] = f_15 * sii_721[k]
                    + f_1 * skh0_687[k]
                    - f_2 * skh1_687[k]
                    + f_3 * pc_y[k] * ski_917[k];

        t_1181[k] = f_16 * sii_693[k]
                    + f_3 * pc_z[k] * ski_917[k];

        t_1182[k] = f_15 * sii_723[k]
                    + f_4 * skh0_689[k]
                    - f_5 * skh1_689[k]
                    + f_3 * pc_y[k] * ski_919[k];
    }

#pragma omp simd aligned(t_1183, t_1184, t_1185, pc_y, sii_724, sii_725, sii_726, skh0_690, \
                         skh0_691, skh0_692, skh1_690, skh1_691, skh1_692, ski_920, ski_921, \
                         ski_922 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1183[k] = f_15 * sii_724[k]
                    + f_6 * skh0_690[k]
                    - f_7 * skh1_690[k]
                    + f_3 * pc_y[k] * ski_920[k];

        t_1184[k] = f_15 * sii_725[k]
                    + f_8 * skh0_691[k]
                    - f_9 * skh1_691[k]
                    + f_3 * pc_y[k] * ski_921[k];

        t_1185[k] = f_15 * sii_726[k]
                    + f_10 * skh0_692[k]
                    - f_11 * skh1_692[k]
                    + f_3 * pc_y[k] * ski_922[k];
    }

#pragma omp simd aligned(t_1186, t_1187, t_1188, t_1189, pc_x, pc_y, pc_z, sii_699, sii_727, \
                         sii_728, skh0_692, skh0_693, skh1_692, skh1_693, ski_923, \
                         ski_924 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1186[k] = f_15 * sii_727[k]
                    + f_3 * pc_y[k] * ski_923[k];

        t_1187[k] = f_16 * sii_699[k]
                    + f_1 * skh0_692[k]
                    - f_2 * skh1_692[k]
                    + f_3 * pc_z[k] * ski_923[k];

        t_1188[k] = f_1 * skh0_693[k]
                    - f_2 * skh1_693[k]
                    + f_3 * pc_x[k] * ski_924[k];

        t_1189[k] = f_14 * sii_728[k]
                    + f_3 * pc_y[k] * ski_924[k];
    }

#pragma omp simd aligned(t_1190, t_1191, t_1192, pc_x, pc_y, pc_z, sii_700, sii_730, skh0_696, \
                         skh1_696, ski_924, ski_926, ski_927 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1190[k] = f_17 * sii_700[k]
                    + f_3 * pc_z[k] * ski_924[k];

        t_1191[k] = f_4 * skh0_696[k]
                    - f_5 * skh1_696[k]
                    + f_3 * pc_x[k] * ski_927[k];

        t_1192[k] = f_14 * sii_730[k]
                    + f_3 * pc_y[k] * ski_926[k];
    }

#pragma omp simd aligned(t_1193, t_1194, t_1195, t_1196, pc_x, pc_y, pc_z, sii_703, sii_733, \
                         skh0_698, skh0_699, skh1_698, skh1_699, ski_927, ski_929, \
                         ski_930 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1193[k] = f_4 * skh0_698[k]
                    - f_5 * skh1_698[k]
                    + f_3 * pc_x[k] * ski_929[k];

        t_1194[k] = f_6 * skh0_699[k]
                    - f_7 * skh1_699[k]
                    + f_3 * pc_x[k] * ski_930[k];

        t_1195[k] = f_17 * sii_703[k]
                    + f_3 * pc_z[k] * ski_927[k];

        t_1196[k] = f_14 * sii_733[k]
                    + f_3 * pc_y[k] * ski_929[k];
    }

#pragma omp simd aligned(t_1197, t_1198, t_1199, pc_x, pc_z, sii_706, skh0_702, skh0_703, \
                         skh1_702, skh1_703, ski_930, ski_933, \
                         ski_934 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1197[k] = f_6 * skh0_702[k]
                    - f_7 * skh1_702[k]
                    + f_3 * pc_x[k] * ski_933[k];

        t_1198[k] = f_8 * skh0_703[k]
                    - f_9 * skh1_703[k]
                    + f_3 * pc_x[k] * ski_934[k];

        t_1199[k] = f_17 * sii_706[k]
                    + f_3 * pc_z[k] * ski_930[k];
    }

#pragma omp simd aligned(t_1200, t_1201, t_1202, pc_x, pc_y, sii_737, skh0_705, skh0_707, \
                         skh1_705, skh1_707, ski_933, ski_936, \
                         ski_938 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1200[k] = f_8 * skh0_705[k]
                    - f_9 * skh1_705[k]
                    + f_3 * pc_x[k] * ski_936[k];

        t_1201[k] = f_14 * sii_737[k]
                    + f_3 * pc_y[k] * ski_933[k];

        t_1202[k] = f_8 * skh0_707[k]
                    - f_9 * skh1_707[k]
                    + f_3 * pc_x[k] * ski_938[k];
    }

#pragma omp simd aligned(t_1203, t_1204, t_1205, pc_x, pc_z, sii_710, skh0_708, skh0_710, \
                         skh1_708, skh1_710, ski_934, ski_939, \
                         ski_941 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1203[k] = f_10 * skh0_708[k]
                    - f_11 * skh1_708[k]
                    + f_3 * pc_x[k] * ski_939[k];

        t_1204[k] = f_17 * sii_710[k]
                    + f_3 * pc_z[k] * ski_934[k];

        t_1205[k] = f_10 * skh0_710[k]
                    - f_11 * skh1_710[k]
                    + f_3 * pc_x[k] * ski_941[k];
    }

#pragma omp simd aligned(t_1206, t_1207, t_1208, t_1209, pc_x, pc_y, sii_742, skh0_711, \
                         skh0_713, skh1_711, skh1_713, ski_938, ski_942, ski_944, \
                         ski_945 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1206[k] = f_10 * skh0_711[k]
                    - f_11 * skh1_711[k]
                    + f_3 * pc_x[k] * ski_942[k];

        t_1207[k] = f_14 * sii_742[k]
                    + f_3 * pc_y[k] * ski_938[k];

        t_1208[k] = f_10 * skh0_713[k]
                    - f_11 * skh1_713[k]
                    + f_3 * pc_x[k] * ski_944[k];

        t_1209[k] = f_3 * pc_x[k] * ski_945[k];
    }

#pragma omp simd aligned(t_1210, t_1211, t_1212, t_1213, t_1214, t_1215, pc_x, ski_946, \
                         ski_947, ski_948, ski_949, ski_950, ski_951 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1210[k] = f_3 * pc_x[k] * ski_946[k];

        t_1211[k] = f_3 * pc_x[k] * ski_947[k];

        t_1212[k] = f_3 * pc_x[k] * ski_948[k];

        t_1213[k] = f_3 * pc_x[k] * ski_949[k];

        t_1214[k] = f_3 * pc_x[k] * ski_950[k];

        t_1215[k] = f_3 * pc_x[k] * ski_951[k];
    }

#pragma omp simd aligned(t_1216, t_1217, t_1218, pc_y, pc_z, sii_721, sii_749, sii_751, \
                         skh0_708, skh0_710, skh1_708, skh1_710, ski_945, \
                         ski_947 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1216[k] = f_14 * sii_749[k]
                    + f_1 * skh0_708[k]
                    - f_2 * skh1_708[k]
                    + f_3 * pc_y[k] * ski_945[k];

        t_1217[k] = f_17 * sii_721[k]
                    + f_3 * pc_z[k] * ski_945[k];

        t_1218[k] = f_14 * sii_751[k]
                    + f_4 * skh0_710[k]
                    - f_5 * skh1_710[k]
                    + f_3 * pc_y[k] * ski_947[k];
    }

#pragma omp simd aligned(t_1219, t_1220, t_1221, pc_y, sii_752, sii_753, sii_754, skh0_711, \
                         skh0_712, skh0_713, skh1_711, skh1_712, skh1_713, ski_948, ski_949, \
                         ski_950 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1219[k] = f_14 * sii_752[k]
                    + f_6 * skh0_711[k]
                    - f_7 * skh1_711[k]
                    + f_3 * pc_y[k] * ski_948[k];

        t_1220[k] = f_14 * sii_753[k]
                    + f_8 * skh0_712[k]
                    - f_9 * skh1_712[k]
                    + f_3 * pc_y[k] * ski_949[k];

        t_1221[k] = f_14 * sii_754[k]
                    + f_10 * skh0_713[k]
                    - f_11 * skh1_713[k]
                    + f_3 * pc_y[k] * ski_950[k];
    }

#pragma omp simd aligned(t_1222, t_1223, t_1224, t_1225, pb_y, pc_y, pc_z, sik0_972, sii_727, \
                         sii_755, sii_756, sik1_972, skh0_713, skh1_713, ski_951, \
                         ski_952 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1222[k] = f_14 * sii_755[k]
                    + f_3 * pc_y[k] * ski_951[k];

        t_1223[k] = f_17 * sii_727[k]
                    + f_1 * skh0_713[k]
                    - f_2 * skh1_713[k]
                    + f_3 * pc_z[k] * ski_951[k];

        t_1224[k] = pb_y[k] * sik0_972[k]
                    - f_12 * pc_y[k] * sik1_972[k];

        t_1225[k] = f_13 * sii_756[k]
                    + f_3 * pc_y[k] * ski_952[k];
    }

#pragma omp simd aligned(t_1226, t_1227, t_1228, pc_x, pc_y, pc_z, sii_728, sii_758, skh0_717, \
                         skh1_717, ski_952, ski_954, ski_955 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1226[k] = f_18 * sii_728[k]
                    + f_3 * pc_z[k] * ski_952[k];

        t_1227[k] = f_4 * skh0_717[k]
                    - f_5 * skh1_717[k]
                    + f_3 * pc_x[k] * ski_955[k];

        t_1228[k] = f_13 * sii_758[k]
                    + f_3 * pc_y[k] * ski_954[k];
    }

#pragma omp simd aligned(t_1229, t_1230, t_1231, pb_y, pc_x, pc_y, pc_z, sik0_977, sii_731, \
                         sik1_977, skh0_720, skh1_720, ski_955, \
                         ski_958 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1229[k] = pb_y[k] * sik0_977[k]
                    - f_12 * pc_y[k] * sik1_977[k];

        t_1230[k] = f_6 * skh0_720[k]
                    - f_7 * skh1_720[k]
                    + f_3 * pc_x[k] * ski_958[k];

        t_1231[k] = f_18 * sii_731[k]
                    + f_3 * pc_z[k] * ski_955[k];
    }

#pragma omp simd aligned(t_1232, t_1233, t_1234, pb_y, pc_x, pc_y, sik0_981, sii_761, \
                         sik1_981, skh0_724, skh1_724, ski_957, \
                         ski_962 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1232[k] = f_13 * sii_761[k]
                    + f_3 * pc_y[k] * ski_957[k];

        t_1233[k] = pb_y[k] * sik0_981[k]
                    - f_12 * pc_y[k] * sik1_981[k];

        t_1234[k] = f_8 * skh0_724[k]
                    - f_9 * skh1_724[k]
                    + f_3 * pc_x[k] * ski_962[k];
    }

#pragma omp simd aligned(t_1235, t_1236, t_1237, pc_x, pc_y, pc_z, sii_734, sii_765, skh0_726, \
                         skh1_726, ski_958, ski_961, ski_964 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1235[k] = f_18 * sii_734[k]
                    + f_3 * pc_z[k] * ski_958[k];

        t_1236[k] = f_8 * skh0_726[k]
                    - f_9 * skh1_726[k]
                    + f_3 * pc_x[k] * ski_964[k];

        t_1237[k] = f_13 * sii_765[k]
                    + f_3 * pc_y[k] * ski_961[k];
    }

#pragma omp simd aligned(t_1238, t_1239, t_1240, pb_y, pc_x, pc_y, pc_z, sik0_986, sii_738, \
                         sik1_986, skh0_729, skh1_729, ski_962, \
                         ski_967 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1238[k] = pb_y[k] * sik0_986[k]
                    - f_12 * pc_y[k] * sik1_986[k];

        t_1239[k] = f_10 * skh0_729[k]
                    - f_11 * skh1_729[k]
                    + f_3 * pc_x[k] * ski_967[k];

        t_1240[k] = f_18 * sii_738[k]
                    + f_3 * pc_z[k] * ski_962[k];
    }

#pragma omp simd aligned(t_1241, t_1242, t_1243, pc_x, pc_y, sii_770, skh0_731, skh0_732, \
                         skh1_731, skh1_732, ski_966, ski_969, \
                         ski_970 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1241[k] = f_10 * skh0_731[k]
                    - f_11 * skh1_731[k]
                    + f_3 * pc_x[k] * ski_969[k];

        t_1242[k] = f_10 * skh0_732[k]
                    - f_11 * skh1_732[k]
                    + f_3 * pc_x[k] * ski_970[k];

        t_1243[k] = f_13 * sii_770[k]
                    + f_3 * pc_y[k] * ski_966[k];
    }

#pragma omp simd aligned(t_1244, t_1245, t_1246, t_1247, t_1248, t_1249, pb_y, pc_x, pc_y, \
                         sik0_992, sik1_992, ski_973, ski_974, ski_975, ski_976, \
                         ski_977 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1244[k] = pb_y[k] * sik0_992[k]
                    - f_12 * pc_y[k] * sik1_992[k];

        t_1245[k] = f_3 * pc_x[k] * ski_973[k];

        t_1246[k] = f_3 * pc_x[k] * ski_974[k];

        t_1247[k] = f_3 * pc_x[k] * ski_975[k];

        t_1248[k] = f_3 * pc_x[k] * ski_976[k];

        t_1249[k] = f_3 * pc_x[k] * ski_977[k];
    }

#pragma omp simd aligned(t_1250, t_1251, t_1252, t_1253, pb_y, pc_x, pc_y, pc_z, sik0_1000, \
                         sii_749, sii_777, sik1_1000, ski_973, ski_978, \
                         ski_979 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1250[k] = f_3 * pc_x[k] * ski_978[k];

        t_1251[k] = f_3 * pc_x[k] * ski_979[k];

        t_1252[k] = pb_y[k] * sik0_1000[k]
                    + f_0 * sii_777[k]
                    - f_12 * pc_y[k] * sik1_1000[k];

        t_1253[k] = f_18 * sii_749[k]
                    + f_3 * pc_z[k] * ski_973[k];
    }

#pragma omp simd aligned(t_1254, t_1255, t_1256, pb_y, pc_y, sik0_1002, sik0_1003, sik0_1004, \
                         sii_779, sii_780, sii_781, sik1_1002, sik1_1003, \
                         sik1_1004 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1254[k] = pb_y[k] * sik0_1002[k]
                    + f_17 * sii_779[k]
                    - f_12 * pc_y[k] * sik1_1002[k];

        t_1255[k] = pb_y[k] * sik0_1003[k]
                    + f_16 * sii_780[k]
                    - f_12 * pc_y[k] * sik1_1003[k];

        t_1256[k] = pb_y[k] * sik0_1004[k]
                    + f_15 * sii_781[k]
                    - f_12 * pc_y[k] * sik1_1004[k];
    }

#pragma omp simd aligned(t_1257, t_1258, t_1259, pb_y, pc_y, sik0_1005, sik0_1007, sii_782, \
                         sii_783, sik1_1005, sik1_1007, ski_979 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1257[k] = pb_y[k] * sik0_1005[k]
                    + f_14 * sii_782[k]
                    - f_12 * pc_y[k] * sik1_1005[k];

        t_1258[k] = f_13 * sii_783[k]
                    + f_3 * pc_y[k] * ski_979[k];

        t_1259[k] = pb_y[k] * sik0_1007[k]
                    - f_12 * pc_y[k] * sik1_1007[k];
    }

#pragma omp simd aligned(t_1260, t_1261, t_1262, t_1263, t_1264, pc_x, pc_y, pc_z, sii_756, \
                         skh0_735, skh0_738, skh1_735, skh1_738, ski_980, ski_982, \
                         ski_983 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1260[k] = f_1 * skh0_735[k]
                    - f_2 * skh1_735[k]
                    + f_3 * pc_x[k] * ski_980[k];

        t_1261[k] = f_3 * pc_y[k] * ski_980[k];

        t_1262[k] = f_0 * sii_756[k]
                    + f_3 * pc_z[k] * ski_980[k];

        t_1263[k] = f_4 * skh0_738[k]
                    - f_5 * skh1_738[k]
                    + f_3 * pc_x[k] * ski_983[k];

        t_1264[k] = f_3 * pc_y[k] * ski_982[k];
    }

#pragma omp simd aligned(t_1265, t_1266, t_1267, t_1268, pc_x, pc_y, pc_z, sii_759, skh0_740, \
                         skh0_741, skh1_740, skh1_741, ski_983, ski_985, \
                         ski_986 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1265[k] = f_4 * skh0_740[k]
                    - f_5 * skh1_740[k]
                    + f_3 * pc_x[k] * ski_985[k];

        t_1266[k] = f_6 * skh0_741[k]
                    - f_7 * skh1_741[k]
                    + f_3 * pc_x[k] * ski_986[k];

        t_1267[k] = f_0 * sii_759[k]
                    + f_3 * pc_z[k] * ski_983[k];

        t_1268[k] = f_3 * pc_y[k] * ski_985[k];
    }

#pragma omp simd aligned(t_1269, t_1270, t_1271, pc_x, pc_z, sii_762, skh0_744, skh0_745, \
                         skh1_744, skh1_745, ski_986, ski_989, \
                         ski_990 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1269[k] = f_6 * skh0_744[k]
                    - f_7 * skh1_744[k]
                    + f_3 * pc_x[k] * ski_989[k];

        t_1270[k] = f_8 * skh0_745[k]
                    - f_9 * skh1_745[k]
                    + f_3 * pc_x[k] * ski_990[k];

        t_1271[k] = f_0 * sii_762[k]
                    + f_3 * pc_z[k] * ski_986[k];
    }

#pragma omp simd aligned(t_1272, t_1273, t_1274, t_1275, pc_x, pc_y, skh0_747, skh0_749, \
                         skh0_750, skh1_747, skh1_749, skh1_750, ski_989, ski_992, ski_994, \
                         ski_995 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1272[k] = f_8 * skh0_747[k]
                    - f_9 * skh1_747[k]
                    + f_3 * pc_x[k] * ski_992[k];

        t_1273[k] = f_3 * pc_y[k] * ski_989[k];

        t_1274[k] = f_8 * skh0_749[k]
                    - f_9 * skh1_749[k]
                    + f_3 * pc_x[k] * ski_994[k];

        t_1275[k] = f_10 * skh0_750[k]
                    - f_11 * skh1_750[k]
                    + f_3 * pc_x[k] * ski_995[k];
    }

#pragma omp simd aligned(t_1276, t_1277, t_1278, t_1279, pc_x, pc_y, pc_z, sii_766, skh0_752, \
                         skh0_753, skh1_752, skh1_753, ski_990, ski_994, ski_997, \
                         ski_998 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1276[k] = f_0 * sii_766[k]
                    + f_3 * pc_z[k] * ski_990[k];

        t_1277[k] = f_10 * skh0_752[k]
                    - f_11 * skh1_752[k]
                    + f_3 * pc_x[k] * ski_997[k];

        t_1278[k] = f_10 * skh0_753[k]
                    - f_11 * skh1_753[k]
                    + f_3 * pc_x[k] * ski_998[k];

        t_1279[k] = f_3 * pc_y[k] * ski_994[k];
    }
}

static auto
compute_prim_skk_three_center_electron_repulsion_0_piece11(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pc,
                                                           const size_t sii, const size_t skh0,
                                                           const size_t skh1, const size_t ski,
                                                           const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / q;
    const auto f_1 = 3.0 / gamma;
    const auto f_2 = 3.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.0 / gamma;
    const auto f_5 = 2.0 * p / (gamma * q);
    const auto f_6 = 1.5 / gamma;
    const auto f_7 = 1.5 * p / (gamma * q);
    const auto f_8 = 1.0 / gamma;
    const auto f_9 = p / (gamma * q);
    const auto f_10 = 0.5 / gamma;
    const auto f_11 = 0.5 * p / (gamma * q);

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

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sii_777 = buffer.data(sii + 777);
    const auto *sii_783 = buffer.data(sii + 783);

    const auto *skh0_750 = buffer.data(skh0 + 750);
    const auto *skh0_752 = buffer.data(skh0 + 752);
    const auto *skh0_753 = buffer.data(skh0 + 753);
    const auto *skh0_754 = buffer.data(skh0 + 754);
    const auto *skh0_755 = buffer.data(skh0 + 755);

    const auto *skh1_750 = buffer.data(skh1 + 750);
    const auto *skh1_752 = buffer.data(skh1 + 752);
    const auto *skh1_753 = buffer.data(skh1 + 753);
    const auto *skh1_754 = buffer.data(skh1 + 754);
    const auto *skh1_755 = buffer.data(skh1 + 755);

    const auto *ski_1000 = buffer.data(ski + 1000);
    const auto *ski_1001 = buffer.data(ski + 1001);
    const auto *ski_1002 = buffer.data(ski + 1002);
    const auto *ski_1003 = buffer.data(ski + 1003);
    const auto *ski_1004 = buffer.data(ski + 1004);
    const auto *ski_1005 = buffer.data(ski + 1005);
    const auto *ski_1006 = buffer.data(ski + 1006);
    const auto *ski_1007 = buffer.data(ski + 1007);

#pragma omp simd aligned(t_1280, t_1281, t_1282, t_1283, t_1284, t_1285, pc_x, skh0_755, \
                         skh1_755, ski_1000, ski_1001, ski_1002, ski_1003, ski_1004, \
                         ski_1005 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1280[k] = f_10 * skh0_755[k]
                    - f_11 * skh1_755[k]
                    + f_3 * pc_x[k] * ski_1000[k];

        t_1281[k] = f_3 * pc_x[k] * ski_1001[k];

        t_1282[k] = f_3 * pc_x[k] * ski_1002[k];

        t_1283[k] = f_3 * pc_x[k] * ski_1003[k];

        t_1284[k] = f_3 * pc_x[k] * ski_1004[k];

        t_1285[k] = f_3 * pc_x[k] * ski_1005[k];
    }

#pragma omp simd aligned(t_1286, t_1287, t_1288, t_1289, pc_x, pc_y, pc_z, sii_777, skh0_750, \
                         skh1_750, ski_1001, ski_1006, ski_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1286[k] = f_3 * pc_x[k] * ski_1006[k];

        t_1287[k] = f_3 * pc_x[k] * ski_1007[k];

        t_1288[k] = f_1 * skh0_750[k]
                    - f_2 * skh1_750[k]
                    + f_3 * pc_y[k] * ski_1001[k];

        t_1289[k] = f_0 * sii_777[k]
                    + f_3 * pc_z[k] * ski_1001[k];
    }

#pragma omp simd aligned(t_1290, t_1291, t_1292, pc_y, skh0_752, skh0_753, skh0_754, skh1_752, \
                         skh1_753, skh1_754, ski_1003, ski_1004, \
                         ski_1005 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1290[k] = f_4 * skh0_752[k]
                    - f_5 * skh1_752[k]
                    + f_3 * pc_y[k] * ski_1003[k];

        t_1291[k] = f_6 * skh0_753[k]
                    - f_7 * skh1_753[k]
                    + f_3 * pc_y[k] * ski_1004[k];

        t_1292[k] = f_8 * skh0_754[k]
                    - f_9 * skh1_754[k]
                    + f_3 * pc_y[k] * ski_1005[k];
    }

#pragma omp simd aligned(t_1293, t_1294, t_1295, pc_y, pc_z, sii_783, skh0_755, skh1_755, \
                         ski_1006, ski_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1293[k] = f_10 * skh0_755[k]
                    - f_11 * skh1_755[k]
                    + f_3 * pc_y[k] * ski_1006[k];

        t_1294[k] = f_3 * pc_y[k] * ski_1007[k];

        t_1295[k] = f_0 * sii_783[k]
                    + f_1 * skh0_755[k]
                    - f_2 * skh1_755[k]
                    + f_3 * pc_z[k] * ski_1007[k];
    }
}

auto
compute_prim_skk_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t sik0, const size_t sii,
                                                   const size_t sik1, const size_t skh0,
                                                   const size_t skh1, const size_t ski,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_skk_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, sik0, sii,
                                                              sik1, skh0, skh1, ski, ncols,
                                                              gamma, p, q);

    compute_prim_skk_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, sik0, sii,
                                                              sik1, skh0, skh1, ski, ncols,
                                                              gamma, p, q);

    compute_prim_skk_three_center_electron_repulsion_0_piece2(buffer, target, pb, pc, sik0, sii,
                                                              sik1, skh0, skh1, ski, ncols,
                                                              gamma, p, q);

    compute_prim_skk_three_center_electron_repulsion_0_piece3(buffer, target, pb, pc, sik0, sii,
                                                              sik1, skh0, skh1, ski, ncols,
                                                              gamma, p, q);

    compute_prim_skk_three_center_electron_repulsion_0_piece4(buffer, target, pb, pc, sik0, sii,
                                                              sik1, skh0, skh1, ski, ncols,
                                                              gamma, p, q);

    compute_prim_skk_three_center_electron_repulsion_0_piece5(buffer, target, pb, pc, sik0, sii,
                                                              sik1, skh0, skh1, ski, ncols,
                                                              gamma, p, q);

    compute_prim_skk_three_center_electron_repulsion_0_piece6(buffer, target, pb, pc, sik0, sii,
                                                              sik1, skh0, skh1, ski, ncols,
                                                              gamma, p, q);

    compute_prim_skk_three_center_electron_repulsion_0_piece7(buffer, target, pb, pc, sik0, sii,
                                                              sik1, ski, ncols, gamma, p, q);

    compute_prim_skk_three_center_electron_repulsion_0_piece8(buffer, target, pb, pc, sik0, sii,
                                                              sik1, skh0, skh1, ski, ncols,
                                                              gamma, p, q);

    compute_prim_skk_three_center_electron_repulsion_0_piece9(buffer, target, pb, pc, sik0, sii,
                                                              sik1, skh0, skh1, ski, ncols,
                                                              gamma, p, q);

    compute_prim_skk_three_center_electron_repulsion_0_piece10(buffer, target, pb, pc, sik0,
                                                               sii, sik1, skh0, skh1, ski,
                                                               ncols, gamma, p, q);

    compute_prim_skk_three_center_electron_repulsion_0_piece11(buffer, target, pc, sii, skh0,
                                                               skh1, ski, ncols, gamma, p, q);
}

}  // namespace simdt3ceri
