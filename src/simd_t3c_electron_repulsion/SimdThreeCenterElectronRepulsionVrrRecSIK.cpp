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


#include "SimdThreeCenterElectronRepulsionVrrRecSIK.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_sik_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t shk0,
                                                          const size_t shi, const size_t shk1,
                                                          const size_t sih0, const size_t sih1,
                                                          const size_t sii, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
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

    const auto *shk0_0 = buffer.data(shk0 + 0);
    const auto *shk0_3 = buffer.data(shk0 + 3);
    const auto *shk0_5 = buffer.data(shk0 + 5);
    const auto *shk0_6 = buffer.data(shk0 + 6);
    const auto *shk0_9 = buffer.data(shk0 + 9);
    const auto *shk0_10 = buffer.data(shk0 + 10);
    const auto *shk0_12 = buffer.data(shk0 + 12);
    const auto *shk0_14 = buffer.data(shk0 + 14);
    const auto *shk0_15 = buffer.data(shk0 + 15);
    const auto *shk0_17 = buffer.data(shk0 + 17);
    const auto *shk0_18 = buffer.data(shk0 + 18);
    const auto *shk0_20 = buffer.data(shk0 + 20);
    const auto *shk0_28 = buffer.data(shk0 + 28);
    const auto *shk0_35 = buffer.data(shk0 + 35);

    const auto *shi_0 = buffer.data(shi + 0);
    const auto *shi_1 = buffer.data(shi + 1);
    const auto *shi_2 = buffer.data(shi + 2);
    const auto *shi_3 = buffer.data(shi + 3);
    const auto *shi_5 = buffer.data(shi + 5);
    const auto *shi_6 = buffer.data(shi + 6);
    const auto *shi_7 = buffer.data(shi + 7);
    const auto *shi_8 = buffer.data(shi + 8);
    const auto *shi_9 = buffer.data(shi + 9);
    const auto *shi_10 = buffer.data(shi + 10);
    const auto *shi_11 = buffer.data(shi + 11);
    const auto *shi_12 = buffer.data(shi + 12);
    const auto *shi_13 = buffer.data(shi + 13);
    const auto *shi_14 = buffer.data(shi + 14);
    const auto *shi_15 = buffer.data(shi + 15);
    const auto *shi_17 = buffer.data(shi + 17);
    const auto *shi_18 = buffer.data(shi + 18);
    const auto *shi_20 = buffer.data(shi + 20);
    const auto *shi_21 = buffer.data(shi + 21);
    const auto *shi_22 = buffer.data(shi + 22);
    const auto *shi_23 = buffer.data(shi + 23);
    const auto *shi_24 = buffer.data(shi + 24);
    const auto *shi_25 = buffer.data(shi + 25);
    const auto *shi_26 = buffer.data(shi + 26);
    const auto *shi_27 = buffer.data(shi + 27);
    const auto *shi_28 = buffer.data(shi + 28);
    const auto *shi_30 = buffer.data(shi + 30);
    const auto *shi_33 = buffer.data(shi + 33);
    const auto *shi_37 = buffer.data(shi + 37);
    const auto *shi_49 = buffer.data(shi + 49);
    const auto *shi_50 = buffer.data(shi + 50);
    const auto *shi_51 = buffer.data(shi + 51);
    const auto *shi_52 = buffer.data(shi + 52);
    const auto *shi_53 = buffer.data(shi + 53);
    const auto *shi_54 = buffer.data(shi + 54);
    const auto *shi_55 = buffer.data(shi + 55);
    const auto *shi_77 = buffer.data(shi + 77);
    const auto *shi_78 = buffer.data(shi + 78);
    const auto *shi_79 = buffer.data(shi + 79);
    const auto *shi_80 = buffer.data(shi + 80);
    const auto *shi_81 = buffer.data(shi + 81);
    const auto *shi_82 = buffer.data(shi + 82);
    const auto *shi_83 = buffer.data(shi + 83);
    const auto *shi_84 = buffer.data(shi + 84);
    const auto *shi_87 = buffer.data(shi + 87);
    const auto *shi_89 = buffer.data(shi + 89);
    const auto *shi_90 = buffer.data(shi + 90);
    const auto *shi_93 = buffer.data(shi + 93);
    const auto *shi_94 = buffer.data(shi + 94);
    const auto *shi_96 = buffer.data(shi + 96);

    const auto *shk1_0 = buffer.data(shk1 + 0);
    const auto *shk1_3 = buffer.data(shk1 + 3);
    const auto *shk1_5 = buffer.data(shk1 + 5);
    const auto *shk1_6 = buffer.data(shk1 + 6);
    const auto *shk1_9 = buffer.data(shk1 + 9);
    const auto *shk1_10 = buffer.data(shk1 + 10);
    const auto *shk1_12 = buffer.data(shk1 + 12);
    const auto *shk1_14 = buffer.data(shk1 + 14);
    const auto *shk1_15 = buffer.data(shk1 + 15);
    const auto *shk1_17 = buffer.data(shk1 + 17);
    const auto *shk1_18 = buffer.data(shk1 + 18);
    const auto *shk1_20 = buffer.data(shk1 + 20);
    const auto *shk1_28 = buffer.data(shk1 + 28);
    const auto *shk1_35 = buffer.data(shk1 + 35);

    const auto *sih0_0 = buffer.data(sih0 + 0);
    const auto *sih0_3 = buffer.data(sih0 + 3);
    const auto *sih0_5 = buffer.data(sih0 + 5);
    const auto *sih0_6 = buffer.data(sih0 + 6);
    const auto *sih0_9 = buffer.data(sih0 + 9);
    const auto *sih0_10 = buffer.data(sih0 + 10);
    const auto *sih0_12 = buffer.data(sih0 + 12);
    const auto *sih0_14 = buffer.data(sih0 + 14);
    const auto *sih0_15 = buffer.data(sih0 + 15);
    const auto *sih0_17 = buffer.data(sih0 + 17);
    const auto *sih0_18 = buffer.data(sih0 + 18);
    const auto *sih0_19 = buffer.data(sih0 + 19);
    const auto *sih0_20 = buffer.data(sih0 + 20);
    const auto *sih0_36 = buffer.data(sih0 + 36);
    const auto *sih0_38 = buffer.data(sih0 + 38);
    const auto *sih0_39 = buffer.data(sih0 + 39);
    const auto *sih0_40 = buffer.data(sih0 + 40);
    const auto *sih0_41 = buffer.data(sih0 + 41);
    const auto *sih0_59 = buffer.data(sih0 + 59);
    const auto *sih0_60 = buffer.data(sih0 + 60);
    const auto *sih0_61 = buffer.data(sih0 + 61);
    const auto *sih0_62 = buffer.data(sih0 + 62);
    const auto *sih0_63 = buffer.data(sih0 + 63);
    const auto *sih0_66 = buffer.data(sih0 + 66);
    const auto *sih0_68 = buffer.data(sih0 + 68);
    const auto *sih0_69 = buffer.data(sih0 + 69);
    const auto *sih0_72 = buffer.data(sih0 + 72);
    const auto *sih0_73 = buffer.data(sih0 + 73);
    const auto *sih0_75 = buffer.data(sih0 + 75);

    const auto *sih1_0 = buffer.data(sih1 + 0);
    const auto *sih1_3 = buffer.data(sih1 + 3);
    const auto *sih1_5 = buffer.data(sih1 + 5);
    const auto *sih1_6 = buffer.data(sih1 + 6);
    const auto *sih1_9 = buffer.data(sih1 + 9);
    const auto *sih1_10 = buffer.data(sih1 + 10);
    const auto *sih1_12 = buffer.data(sih1 + 12);
    const auto *sih1_14 = buffer.data(sih1 + 14);
    const auto *sih1_15 = buffer.data(sih1 + 15);
    const auto *sih1_17 = buffer.data(sih1 + 17);
    const auto *sih1_18 = buffer.data(sih1 + 18);
    const auto *sih1_19 = buffer.data(sih1 + 19);
    const auto *sih1_20 = buffer.data(sih1 + 20);
    const auto *sih1_36 = buffer.data(sih1 + 36);
    const auto *sih1_38 = buffer.data(sih1 + 38);
    const auto *sih1_39 = buffer.data(sih1 + 39);
    const auto *sih1_40 = buffer.data(sih1 + 40);
    const auto *sih1_41 = buffer.data(sih1 + 41);
    const auto *sih1_59 = buffer.data(sih1 + 59);
    const auto *sih1_60 = buffer.data(sih1 + 60);
    const auto *sih1_61 = buffer.data(sih1 + 61);
    const auto *sih1_62 = buffer.data(sih1 + 62);
    const auto *sih1_63 = buffer.data(sih1 + 63);
    const auto *sih1_66 = buffer.data(sih1 + 66);
    const auto *sih1_68 = buffer.data(sih1 + 68);
    const auto *sih1_69 = buffer.data(sih1 + 69);
    const auto *sih1_72 = buffer.data(sih1 + 72);
    const auto *sih1_73 = buffer.data(sih1 + 73);
    const auto *sih1_75 = buffer.data(sih1 + 75);

    const auto *sii_0 = buffer.data(sii + 0);
    const auto *sii_2 = buffer.data(sii + 2);
    const auto *sii_3 = buffer.data(sii + 3);
    const auto *sii_5 = buffer.data(sii + 5);
    const auto *sii_6 = buffer.data(sii + 6);
    const auto *sii_9 = buffer.data(sii + 9);
    const auto *sii_10 = buffer.data(sii + 10);
    const auto *sii_12 = buffer.data(sii + 12);
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
    const auto *sii_31 = buffer.data(sii + 31);
    const auto *sii_33 = buffer.data(sii + 33);
    const auto *sii_34 = buffer.data(sii + 34);
    const auto *sii_37 = buffer.data(sii + 37);
    const auto *sii_38 = buffer.data(sii + 38);
    const auto *sii_42 = buffer.data(sii + 42);
    const auto *sii_49 = buffer.data(sii + 49);
    const auto *sii_50 = buffer.data(sii + 50);
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
    const auto *sii_65 = buffer.data(sii + 65);
    const auto *sii_66 = buffer.data(sii + 66);
    const auto *sii_70 = buffer.data(sii + 70);
    const auto *sii_77 = buffer.data(sii + 77);
    const auto *sii_78 = buffer.data(sii + 78);
    const auto *sii_79 = buffer.data(sii + 79);
    const auto *sii_80 = buffer.data(sii + 80);
    const auto *sii_81 = buffer.data(sii + 81);
    const auto *sii_82 = buffer.data(sii + 82);
    const auto *sii_83 = buffer.data(sii + 83);
    const auto *sii_84 = buffer.data(sii + 84);
    const auto *sii_86 = buffer.data(sii + 86);
    const auto *sii_87 = buffer.data(sii + 87);
    const auto *sii_89 = buffer.data(sii + 89);
    const auto *sii_90 = buffer.data(sii + 90);
    const auto *sii_93 = buffer.data(sii + 93);
    const auto *sii_94 = buffer.data(sii + 94);
    const auto *sii_96 = buffer.data(sii + 96);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, shi_0, shi_3, sih0_0, sih0_3, \
                         sih1_0, sih1_3, sii_0, sii_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * shi_0[k]
                 + f_1 * sih0_0[k]
                 - f_2 * sih1_0[k]
                 + f_3 * pc_x[k] * sii_0[k];

        t_1[k] = f_3 * pc_y[k] * sii_0[k];

        t_2[k] = f_3 * pc_z[k] * sii_0[k];

        t_3[k] = f_0 * shi_3[k]
                 + f_4 * sih0_3[k]
                 - f_5 * sih1_3[k]
                 + f_3 * pc_x[k] * sii_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pc_x, pc_y, shi_5, shi_6, sih0_5, sih0_6, sih1_5, \
                         sih1_6, sii_2, sii_5, sii_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * pc_y[k] * sii_2[k];

        t_5[k] = f_0 * shi_5[k]
                 + f_4 * sih0_5[k]
                 - f_5 * sih1_5[k]
                 + f_3 * pc_x[k] * sii_5[k];

        t_6[k] = f_0 * shi_6[k]
                 + f_6 * sih0_6[k]
                 - f_7 * sih1_6[k]
                 + f_3 * pc_x[k] * sii_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pc_x, pc_y, pc_z, shi_9, sih0_9, sih1_9, sii_3, sii_5, \
                         sii_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_3 * pc_z[k] * sii_3[k];

        t_8[k] = f_3 * pc_y[k] * sii_5[k];

        t_9[k] = f_0 * shi_9[k]
                 + f_6 * sih0_9[k]
                 - f_7 * sih1_9[k]
                 + f_3 * pc_x[k] * sii_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pc_x, pc_z, shi_10, shi_12, sih0_10, sih0_12, \
                         sih1_10, sih1_12, sii_6, sii_10, sii_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * shi_10[k]
                  + f_8 * sih0_10[k]
                  - f_9 * sih1_10[k]
                  + f_3 * pc_x[k] * sii_10[k];

        t_11[k] = f_3 * pc_z[k] * sii_6[k];

        t_12[k] = f_0 * shi_12[k]
                  + f_8 * sih0_12[k]
                  - f_9 * sih1_12[k]
                  + f_3 * pc_x[k] * sii_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pc_x, pc_y, shi_14, shi_15, sih0_14, sih0_15, \
                         sih1_14, sih1_15, sii_9, sii_14, sii_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_3 * pc_y[k] * sii_9[k];

        t_14[k] = f_0 * shi_14[k]
                  + f_8 * sih0_14[k]
                  - f_9 * sih1_14[k]
                  + f_3 * pc_x[k] * sii_14[k];

        t_15[k] = f_0 * shi_15[k]
                  + f_10 * sih0_15[k]
                  - f_11 * sih1_15[k]
                  + f_3 * pc_x[k] * sii_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pc_x, pc_z, shi_17, shi_18, sih0_17, sih0_18, \
                         sih1_17, sih1_18, sii_10, sii_17, sii_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * pc_z[k] * sii_10[k];

        t_17[k] = f_0 * shi_17[k]
                  + f_10 * sih0_17[k]
                  - f_11 * sih1_17[k]
                  + f_3 * pc_x[k] * sii_17[k];

        t_18[k] = f_0 * shi_18[k]
                  + f_10 * sih0_18[k]
                  - f_11 * sih1_18[k]
                  + f_3 * pc_x[k] * sii_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pc_x, pc_y, shi_20, shi_21, shi_22, sih0_20, \
                         sih1_20, sii_14, sii_20, sii_21, sii_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * pc_y[k] * sii_14[k];

        t_20[k] = f_0 * shi_20[k]
                  + f_10 * sih0_20[k]
                  - f_11 * sih1_20[k]
                  + f_3 * pc_x[k] * sii_20[k];

        t_21[k] = f_0 * shi_21[k]
                  + f_3 * pc_x[k] * sii_21[k];

        t_22[k] = f_0 * shi_22[k]
                  + f_3 * pc_x[k] * sii_22[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, t_27, pc_x, shi_23, shi_24, shi_25, shi_26, \
                         shi_27, sii_23, sii_24, sii_25, sii_26, \
                         sii_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_0 * shi_23[k]
                  + f_3 * pc_x[k] * sii_23[k];

        t_24[k] = f_0 * shi_24[k]
                  + f_3 * pc_x[k] * sii_24[k];

        t_25[k] = f_0 * shi_25[k]
                  + f_3 * pc_x[k] * sii_25[k];

        t_26[k] = f_0 * shi_26[k]
                  + f_3 * pc_x[k] * sii_26[k];

        t_27[k] = f_0 * shi_27[k]
                  + f_3 * pc_x[k] * sii_27[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pc_y, pc_z, sih0_15, sih0_17, sih0_18, \
                         sih1_15, sih1_17, sih1_18, sii_21, sii_23, \
                         sii_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_1 * sih0_15[k]
                  - f_2 * sih1_15[k]
                  + f_3 * pc_y[k] * sii_21[k];

        t_29[k] = f_3 * pc_z[k] * sii_21[k];

        t_30[k] = f_4 * sih0_17[k]
                  - f_5 * sih1_17[k]
                  + f_3 * pc_y[k] * sii_23[k];

        t_31[k] = f_6 * sih0_18[k]
                  - f_7 * sih1_18[k]
                  + f_3 * pc_y[k] * sii_24[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pc_y, pc_z, sih0_19, sih0_20, sih1_19, \
                         sih1_20, sii_25, sii_26, sii_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_8 * sih0_19[k]
                  - f_9 * sih1_19[k]
                  + f_3 * pc_y[k] * sii_25[k];

        t_33[k] = f_10 * sih0_20[k]
                  - f_11 * sih1_20[k]
                  + f_3 * pc_y[k] * sii_26[k];

        t_34[k] = f_3 * pc_y[k] * sii_27[k];

        t_35[k] = f_1 * sih0_20[k]
                  - f_2 * sih1_20[k]
                  + f_3 * pc_z[k] * sii_27[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pb_y, pc_y, pc_z, shk0_0, shk0_3, shi_0, \
                         shi_1, shk1_0, shk1_3, sii_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = pb_y[k] * shk0_0[k]
                  - f_12 * pc_y[k] * shk1_0[k];

        t_37[k] = f_13 * shi_0[k]
                  + f_3 * pc_y[k] * sii_28[k];

        t_38[k] = f_3 * pc_z[k] * sii_28[k];

        t_39[k] = pb_y[k] * shk0_3[k]
                  + f_14 * shi_1[k]
                  - f_12 * pc_y[k] * shk1_3[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pb_y, pc_y, pc_z, shk0_5, shk0_6, shi_2, \
                         shi_3, shk1_5, shk1_6, sii_30, sii_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_13 * shi_2[k]
                  + f_3 * pc_y[k] * sii_30[k];

        t_41[k] = pb_y[k] * shk0_5[k]
                  - f_12 * pc_y[k] * shk1_5[k];

        t_42[k] = pb_y[k] * shk0_6[k]
                  + f_15 * shi_3[k]
                  - f_12 * pc_y[k] * shk1_6[k];

        t_43[k] = f_3 * pc_z[k] * sii_31[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pb_y, pc_y, pc_z, shk0_9, shk0_10, shi_5, \
                         shi_6, shk1_9, shk1_10, sii_33, sii_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_13 * shi_5[k]
                  + f_3 * pc_y[k] * sii_33[k];

        t_45[k] = pb_y[k] * shk0_9[k]
                  - f_12 * pc_y[k] * shk1_9[k];

        t_46[k] = pb_y[k] * shk0_10[k]
                  + f_16 * shi_6[k]
                  - f_12 * pc_y[k] * shk1_10[k];

        t_47[k] = f_3 * pc_z[k] * sii_34[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pb_y, pc_y, shk0_12, shk0_14, shk0_15, shi_8, \
                         shi_9, shi_10, shk1_12, shk1_14, shk1_15, \
                         sii_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = pb_y[k] * shk0_12[k]
                  + f_14 * shi_8[k]
                  - f_12 * pc_y[k] * shk1_12[k];

        t_49[k] = f_13 * shi_9[k]
                  + f_3 * pc_y[k] * sii_37[k];

        t_50[k] = pb_y[k] * shk0_14[k]
                  - f_12 * pc_y[k] * shk1_14[k];

        t_51[k] = pb_y[k] * shk0_15[k]
                  + f_17 * shi_10[k]
                  - f_12 * pc_y[k] * shk1_15[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pb_y, pc_y, pc_z, shk0_17, shk0_18, shi_12, \
                         shi_13, shi_14, shk1_17, shk1_18, sii_38, \
                         sii_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_3 * pc_z[k] * sii_38[k];

        t_53[k] = pb_y[k] * shk0_17[k]
                  + f_15 * shi_12[k]
                  - f_12 * pc_y[k] * shk1_17[k];

        t_54[k] = pb_y[k] * shk0_18[k]
                  + f_14 * shi_13[k]
                  - f_12 * pc_y[k] * shk1_18[k];

        t_55[k] = f_13 * shi_14[k]
                  + f_3 * pc_y[k] * sii_42[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pb_y, pc_x, pc_y, shk0_20, shi_49, shi_50, \
                         shi_51, shk1_20, sii_49, sii_50, sii_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_y[k] * shk0_20[k]
                  - f_12 * pc_y[k] * shk1_20[k];

        t_57[k] = f_17 * shi_49[k]
                  + f_3 * pc_x[k] * sii_49[k];

        t_58[k] = f_17 * shi_50[k]
                  + f_3 * pc_x[k] * sii_50[k];

        t_59[k] = f_17 * shi_51[k]
                  + f_3 * pc_x[k] * sii_51[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pc_x, shi_52, shi_53, shi_54, shi_55, sii_52, \
                         sii_53, sii_54, sii_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_17 * shi_52[k]
                  + f_3 * pc_x[k] * sii_52[k];

        t_61[k] = f_17 * shi_53[k]
                  + f_3 * pc_x[k] * sii_53[k];

        t_62[k] = f_17 * shi_54[k]
                  + f_3 * pc_x[k] * sii_54[k];

        t_63[k] = f_17 * shi_55[k]
                  + f_3 * pc_x[k] * sii_55[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, pc_y, pc_z, shi_21, shi_23, sih0_36, sih0_38, \
                         sih1_36, sih1_38, sii_49, sii_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_13 * shi_21[k]
                  + f_1 * sih0_36[k]
                  - f_2 * sih1_36[k]
                  + f_3 * pc_y[k] * sii_49[k];

        t_65[k] = f_3 * pc_z[k] * sii_49[k];

        t_66[k] = f_13 * shi_23[k]
                  + f_4 * sih0_38[k]
                  - f_5 * sih1_38[k]
                  + f_3 * pc_y[k] * sii_51[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, pc_y, shi_24, shi_25, shi_26, sih0_39, sih0_40, \
                         sih0_41, sih1_39, sih1_40, sih1_41, sii_52, sii_53, \
                         sii_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_13 * shi_24[k]
                  + f_6 * sih0_39[k]
                  - f_7 * sih1_39[k]
                  + f_3 * pc_y[k] * sii_52[k];

        t_68[k] = f_13 * shi_25[k]
                  + f_8 * sih0_40[k]
                  - f_9 * sih1_40[k]
                  + f_3 * pc_y[k] * sii_53[k];

        t_69[k] = f_13 * shi_26[k]
                  + f_10 * sih0_41[k]
                  - f_11 * sih1_41[k]
                  + f_3 * pc_y[k] * sii_54[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pb_y, pb_z, pc_y, pc_z, shk0_0, shk0_35, \
                         shi_27, shk1_0, shk1_35, sii_55, sii_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_13 * shi_27[k]
                  + f_3 * pc_y[k] * sii_55[k];

        t_71[k] = pb_y[k] * shk0_35[k]
                  - f_12 * pc_y[k] * shk1_35[k];

        t_72[k] = pb_z[k] * shk0_0[k]
                  - f_12 * pc_z[k] * shk1_0[k];

        t_73[k] = f_3 * pc_y[k] * sii_56[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pb_z, pc_y, pc_z, shk0_3, shk0_5, shi_0, \
                         shi_2, shk1_3, shk1_5, sii_56, sii_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_13 * shi_0[k]
                  + f_3 * pc_z[k] * sii_56[k];

        t_75[k] = pb_z[k] * shk0_3[k]
                  - f_12 * pc_z[k] * shk1_3[k];

        t_76[k] = f_3 * pc_y[k] * sii_58[k];

        t_77[k] = pb_z[k] * shk0_5[k]
                  + f_14 * shi_2[k]
                  - f_12 * pc_z[k] * shk1_5[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pb_z, pc_y, pc_z, shk0_6, shk0_9, shi_3, \
                         shi_5, shk1_6, shk1_9, sii_59, sii_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = pb_z[k] * shk0_6[k]
                  - f_12 * pc_z[k] * shk1_6[k];

        t_79[k] = f_13 * shi_3[k]
                  + f_3 * pc_z[k] * sii_59[k];

        t_80[k] = f_3 * pc_y[k] * sii_61[k];

        t_81[k] = pb_z[k] * shk0_9[k]
                  + f_15 * shi_5[k]
                  - f_12 * pc_z[k] * shk1_9[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, pb_z, pc_y, pc_z, shk0_10, shk0_12, shi_6, \
                         shi_7, shk1_10, shk1_12, sii_62, sii_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = pb_z[k] * shk0_10[k]
                  - f_12 * pc_z[k] * shk1_10[k];

        t_83[k] = f_13 * shi_6[k]
                  + f_3 * pc_z[k] * sii_62[k];

        t_84[k] = pb_z[k] * shk0_12[k]
                  + f_14 * shi_7[k]
                  - f_12 * pc_z[k] * shk1_12[k];

        t_85[k] = f_3 * pc_y[k] * sii_65[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pb_z, pc_z, shk0_14, shk0_15, shk0_17, shi_9, \
                         shi_10, shi_11, shk1_14, shk1_15, shk1_17, \
                         sii_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = pb_z[k] * shk0_14[k]
                  + f_16 * shi_9[k]
                  - f_12 * pc_z[k] * shk1_14[k];

        t_87[k] = pb_z[k] * shk0_15[k]
                  - f_12 * pc_z[k] * shk1_15[k];

        t_88[k] = f_13 * shi_10[k]
                  + f_3 * pc_z[k] * sii_66[k];

        t_89[k] = pb_z[k] * shk0_17[k]
                  + f_14 * shi_11[k]
                  - f_12 * pc_z[k] * shk1_17[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, pb_z, pc_y, pc_z, shk0_18, shk0_20, shi_12, shi_14, \
                         shk1_18, shk1_20, sii_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = pb_z[k] * shk0_18[k]
                  + f_15 * shi_12[k]
                  - f_12 * pc_z[k] * shk1_18[k];

        t_91[k] = f_3 * pc_y[k] * sii_70[k];

        t_92[k] = pb_z[k] * shk0_20[k]
                  + f_17 * shi_14[k]
                  - f_12 * pc_z[k] * shk1_20[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, t_97, pc_x, shi_77, shi_78, shi_79, shi_80, \
                         shi_81, sii_77, sii_78, sii_79, sii_80, \
                         sii_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_17 * shi_77[k]
                  + f_3 * pc_x[k] * sii_77[k];

        t_94[k] = f_17 * shi_78[k]
                  + f_3 * pc_x[k] * sii_78[k];

        t_95[k] = f_17 * shi_79[k]
                  + f_3 * pc_x[k] * sii_79[k];

        t_96[k] = f_17 * shi_80[k]
                  + f_3 * pc_x[k] * sii_80[k];

        t_97[k] = f_17 * shi_81[k]
                  + f_3 * pc_x[k] * sii_81[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, pb_z, pc_x, pc_z, shk0_28, shi_21, shi_82, \
                         shi_83, shk1_28, sii_77, sii_82, sii_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_17 * shi_82[k]
                  + f_3 * pc_x[k] * sii_82[k];

        t_99[k] = f_17 * shi_83[k]
                  + f_3 * pc_x[k] * sii_83[k];

        t_100[k] = pb_z[k] * shk0_28[k]
                   - f_12 * pc_z[k] * shk1_28[k];

        t_101[k] = f_13 * shi_21[k]
                   + f_3 * pc_z[k] * sii_77[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, pc_y, sih0_59, sih0_60, sih0_61, sih1_59, \
                         sih1_60, sih1_61, sii_79, sii_80, sii_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_4 * sih0_59[k]
                   - f_5 * sih1_59[k]
                   + f_3 * pc_y[k] * sii_79[k];

        t_103[k] = f_6 * sih0_60[k]
                   - f_7 * sih1_60[k]
                   + f_3 * pc_y[k] * sii_80[k];

        t_104[k] = f_8 * sih0_61[k]
                   - f_9 * sih1_61[k]
                   + f_3 * pc_y[k] * sii_81[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, pc_x, pc_y, pc_z, shi_27, shi_84, \
                         sih0_62, sih0_63, sih1_62, sih1_63, sii_82, sii_83, \
                         sii_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_10 * sih0_62[k]
                   - f_11 * sih1_62[k]
                   + f_3 * pc_y[k] * sii_82[k];

        t_106[k] = f_3 * pc_y[k] * sii_83[k];

        t_107[k] = f_13 * shi_27[k]
                   + f_1 * sih0_62[k]
                   - f_2 * sih1_62[k]
                   + f_3 * pc_z[k] * sii_83[k];

        t_108[k] = f_16 * shi_84[k]
                   + f_1 * sih0_63[k]
                   - f_2 * sih1_63[k]
                   + f_3 * pc_x[k] * sii_84[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pc_x, pc_y, pc_z, shi_28, shi_30, shi_87, \
                         sih0_66, sih1_66, sii_84, sii_86, sii_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_14 * shi_28[k]
                   + f_3 * pc_y[k] * sii_84[k];

        t_110[k] = f_3 * pc_z[k] * sii_84[k];

        t_111[k] = f_16 * shi_87[k]
                   + f_4 * sih0_66[k]
                   - f_5 * sih1_66[k]
                   + f_3 * pc_x[k] * sii_87[k];

        t_112[k] = f_14 * shi_30[k]
                   + f_3 * pc_y[k] * sii_86[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, pc_x, pc_z, shi_89, shi_90, sih0_68, sih0_69, \
                         sih1_68, sih1_69, sii_87, sii_89, sii_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_16 * shi_89[k]
                   + f_4 * sih0_68[k]
                   - f_5 * sih1_68[k]
                   + f_3 * pc_x[k] * sii_89[k];

        t_114[k] = f_16 * shi_90[k]
                   + f_6 * sih0_69[k]
                   - f_7 * sih1_69[k]
                   + f_3 * pc_x[k] * sii_90[k];

        t_115[k] = f_3 * pc_z[k] * sii_87[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, pc_x, pc_y, shi_33, shi_93, shi_94, sih0_72, \
                         sih0_73, sih1_72, sih1_73, sii_89, sii_93, \
                         sii_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_14 * shi_33[k]
                   + f_3 * pc_y[k] * sii_89[k];

        t_117[k] = f_16 * shi_93[k]
                   + f_6 * sih0_72[k]
                   - f_7 * sih1_72[k]
                   + f_3 * pc_x[k] * sii_93[k];

        t_118[k] = f_16 * shi_94[k]
                   + f_8 * sih0_73[k]
                   - f_9 * sih1_73[k]
                   + f_3 * pc_x[k] * sii_94[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, pc_x, pc_y, pc_z, shi_37, shi_96, sih0_75, \
                         sih1_75, sii_90, sii_93, sii_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = f_3 * pc_z[k] * sii_90[k];

        t_120[k] = f_16 * shi_96[k]
                   + f_8 * sih0_75[k]
                   - f_9 * sih1_75[k]
                   + f_3 * pc_x[k] * sii_96[k];

        t_121[k] = f_14 * shi_37[k]
                   + f_3 * pc_y[k] * sii_93[k];
    }
}

static auto
compute_prim_sik_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t shk0,
                                                          const size_t shi, const size_t shk1,
                                                          const size_t sih0, const size_t sih1,
                                                          const size_t sii, const size_t ncols,
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

    const auto *shk0_39 = buffer.data(shk0 + 39);
    const auto *shk0_42 = buffer.data(shk0 + 42);
    const auto *shk0_46 = buffer.data(shk0 + 46);
    const auto *shk0_51 = buffer.data(shk0 + 51);
    const auto *shk0_64 = buffer.data(shk0 + 64);
    const auto *shk0_72 = buffer.data(shk0 + 72);
    const auto *shk0_77 = buffer.data(shk0 + 77);
    const auto *shk0_81 = buffer.data(shk0 + 81);
    const auto *shk0_84 = buffer.data(shk0 + 84);
    const auto *shk0_86 = buffer.data(shk0 + 86);
    const auto *shk0_89 = buffer.data(shk0 + 89);
    const auto *shk0_90 = buffer.data(shk0 + 90);
    const auto *shk0_92 = buffer.data(shk0 + 92);
    const auto *shk0_107 = buffer.data(shk0 + 107);

    const auto *shi_28 = buffer.data(shi + 28);
    const auto *shi_31 = buffer.data(shi + 31);
    const auto *shi_34 = buffer.data(shi + 34);
    const auto *shi_38 = buffer.data(shi + 38);
    const auto *shi_42 = buffer.data(shi + 42);
    const auto *shi_49 = buffer.data(shi + 49);
    const auto *shi_51 = buffer.data(shi + 51);
    const auto *shi_52 = buffer.data(shi + 52);
    const auto *shi_53 = buffer.data(shi + 53);
    const auto *shi_54 = buffer.data(shi + 54);
    const auto *shi_55 = buffer.data(shi + 55);
    const auto *shi_56 = buffer.data(shi + 56);
    const auto *shi_58 = buffer.data(shi + 58);
    const auto *shi_59 = buffer.data(shi + 59);
    const auto *shi_61 = buffer.data(shi + 61);
    const auto *shi_62 = buffer.data(shi + 62);
    const auto *shi_64 = buffer.data(shi + 64);
    const auto *shi_65 = buffer.data(shi + 65);
    const auto *shi_66 = buffer.data(shi + 66);
    const auto *shi_68 = buffer.data(shi + 68);
    const auto *shi_69 = buffer.data(shi + 69);
    const auto *shi_70 = buffer.data(shi + 70);
    const auto *shi_77 = buffer.data(shi + 77);
    const auto *shi_79 = buffer.data(shi + 79);
    const auto *shi_80 = buffer.data(shi + 80);
    const auto *shi_81 = buffer.data(shi + 81);
    const auto *shi_82 = buffer.data(shi + 82);
    const auto *shi_83 = buffer.data(shi + 83);
    const auto *shi_84 = buffer.data(shi + 84);
    const auto *shi_86 = buffer.data(shi + 86);
    const auto *shi_89 = buffer.data(shi + 89);
    const auto *shi_93 = buffer.data(shi + 93);
    const auto *shi_98 = buffer.data(shi + 98);
    const auto *shi_99 = buffer.data(shi + 99);
    const auto *shi_101 = buffer.data(shi + 101);
    const auto *shi_102 = buffer.data(shi + 102);
    const auto *shi_104 = buffer.data(shi + 104);
    const auto *shi_105 = buffer.data(shi + 105);
    const auto *shi_106 = buffer.data(shi + 106);
    const auto *shi_107 = buffer.data(shi + 107);
    const auto *shi_108 = buffer.data(shi + 108);
    const auto *shi_109 = buffer.data(shi + 109);
    const auto *shi_110 = buffer.data(shi + 110);
    const auto *shi_111 = buffer.data(shi + 111);
    const auto *shi_133 = buffer.data(shi + 133);
    const auto *shi_134 = buffer.data(shi + 134);
    const auto *shi_135 = buffer.data(shi + 135);
    const auto *shi_136 = buffer.data(shi + 136);
    const auto *shi_137 = buffer.data(shi + 137);
    const auto *shi_138 = buffer.data(shi + 138);
    const auto *shi_139 = buffer.data(shi + 139);
    const auto *shi_140 = buffer.data(shi + 140);
    const auto *shi_143 = buffer.data(shi + 143);
    const auto *shi_145 = buffer.data(shi + 145);
    const auto *shi_146 = buffer.data(shi + 146);
    const auto *shi_149 = buffer.data(shi + 149);
    const auto *shi_150 = buffer.data(shi + 150);
    const auto *shi_152 = buffer.data(shi + 152);
    const auto *shi_154 = buffer.data(shi + 154);
    const auto *shi_155 = buffer.data(shi + 155);
    const auto *shi_157 = buffer.data(shi + 157);
    const auto *shi_158 = buffer.data(shi + 158);
    const auto *shi_160 = buffer.data(shi + 160);
    const auto *shi_161 = buffer.data(shi + 161);
    const auto *shi_162 = buffer.data(shi + 162);
    const auto *shi_163 = buffer.data(shi + 163);
    const auto *shi_164 = buffer.data(shi + 164);
    const auto *shi_165 = buffer.data(shi + 165);
    const auto *shi_166 = buffer.data(shi + 166);
    const auto *shi_167 = buffer.data(shi + 167);
    const auto *shi_168 = buffer.data(shi + 168);
    const auto *shi_171 = buffer.data(shi + 171);
    const auto *shi_173 = buffer.data(shi + 173);
    const auto *shi_174 = buffer.data(shi + 174);
    const auto *shi_177 = buffer.data(shi + 177);
    const auto *shi_178 = buffer.data(shi + 178);
    const auto *shi_180 = buffer.data(shi + 180);
    const auto *shi_182 = buffer.data(shi + 182);
    const auto *shi_183 = buffer.data(shi + 183);
    const auto *shi_185 = buffer.data(shi + 185);
    const auto *shi_186 = buffer.data(shi + 186);

    const auto *shk1_39 = buffer.data(shk1 + 39);
    const auto *shk1_42 = buffer.data(shk1 + 42);
    const auto *shk1_46 = buffer.data(shk1 + 46);
    const auto *shk1_51 = buffer.data(shk1 + 51);
    const auto *shk1_64 = buffer.data(shk1 + 64);
    const auto *shk1_72 = buffer.data(shk1 + 72);
    const auto *shk1_77 = buffer.data(shk1 + 77);
    const auto *shk1_81 = buffer.data(shk1 + 81);
    const auto *shk1_84 = buffer.data(shk1 + 84);
    const auto *shk1_86 = buffer.data(shk1 + 86);
    const auto *shk1_89 = buffer.data(shk1 + 89);
    const auto *shk1_90 = buffer.data(shk1 + 90);
    const auto *shk1_92 = buffer.data(shk1 + 92);
    const auto *shk1_107 = buffer.data(shk1 + 107);

    const auto *sih0_77 = buffer.data(sih0 + 77);
    const auto *sih0_78 = buffer.data(sih0 + 78);
    const auto *sih0_80 = buffer.data(sih0 + 80);
    const auto *sih0_81 = buffer.data(sih0 + 81);
    const auto *sih0_82 = buffer.data(sih0 + 82);
    const auto *sih0_83 = buffer.data(sih0 + 83);
    const auto *sih0_101 = buffer.data(sih0 + 101);
    const auto *sih0_102 = buffer.data(sih0 + 102);
    const auto *sih0_103 = buffer.data(sih0 + 103);
    const auto *sih0_104 = buffer.data(sih0 + 104);
    const auto *sih0_105 = buffer.data(sih0 + 105);
    const auto *sih0_108 = buffer.data(sih0 + 108);
    const auto *sih0_110 = buffer.data(sih0 + 110);
    const auto *sih0_111 = buffer.data(sih0 + 111);
    const auto *sih0_114 = buffer.data(sih0 + 114);
    const auto *sih0_115 = buffer.data(sih0 + 115);
    const auto *sih0_117 = buffer.data(sih0 + 117);
    const auto *sih0_119 = buffer.data(sih0 + 119);
    const auto *sih0_120 = buffer.data(sih0 + 120);
    const auto *sih0_122 = buffer.data(sih0 + 122);
    const auto *sih0_123 = buffer.data(sih0 + 123);
    const auto *sih0_124 = buffer.data(sih0 + 124);
    const auto *sih0_125 = buffer.data(sih0 + 125);
    const auto *sih0_126 = buffer.data(sih0 + 126);
    const auto *sih0_129 = buffer.data(sih0 + 129);
    const auto *sih0_131 = buffer.data(sih0 + 131);
    const auto *sih0_132 = buffer.data(sih0 + 132);
    const auto *sih0_135 = buffer.data(sih0 + 135);
    const auto *sih0_136 = buffer.data(sih0 + 136);
    const auto *sih0_138 = buffer.data(sih0 + 138);
    const auto *sih0_140 = buffer.data(sih0 + 140);
    const auto *sih0_141 = buffer.data(sih0 + 141);
    const auto *sih0_143 = buffer.data(sih0 + 143);
    const auto *sih0_144 = buffer.data(sih0 + 144);

    const auto *sih1_77 = buffer.data(sih1 + 77);
    const auto *sih1_78 = buffer.data(sih1 + 78);
    const auto *sih1_80 = buffer.data(sih1 + 80);
    const auto *sih1_81 = buffer.data(sih1 + 81);
    const auto *sih1_82 = buffer.data(sih1 + 82);
    const auto *sih1_83 = buffer.data(sih1 + 83);
    const auto *sih1_101 = buffer.data(sih1 + 101);
    const auto *sih1_102 = buffer.data(sih1 + 102);
    const auto *sih1_103 = buffer.data(sih1 + 103);
    const auto *sih1_104 = buffer.data(sih1 + 104);
    const auto *sih1_105 = buffer.data(sih1 + 105);
    const auto *sih1_108 = buffer.data(sih1 + 108);
    const auto *sih1_110 = buffer.data(sih1 + 110);
    const auto *sih1_111 = buffer.data(sih1 + 111);
    const auto *sih1_114 = buffer.data(sih1 + 114);
    const auto *sih1_115 = buffer.data(sih1 + 115);
    const auto *sih1_117 = buffer.data(sih1 + 117);
    const auto *sih1_119 = buffer.data(sih1 + 119);
    const auto *sih1_120 = buffer.data(sih1 + 120);
    const auto *sih1_122 = buffer.data(sih1 + 122);
    const auto *sih1_123 = buffer.data(sih1 + 123);
    const auto *sih1_124 = buffer.data(sih1 + 124);
    const auto *sih1_125 = buffer.data(sih1 + 125);
    const auto *sih1_126 = buffer.data(sih1 + 126);
    const auto *sih1_129 = buffer.data(sih1 + 129);
    const auto *sih1_131 = buffer.data(sih1 + 131);
    const auto *sih1_132 = buffer.data(sih1 + 132);
    const auto *sih1_135 = buffer.data(sih1 + 135);
    const auto *sih1_136 = buffer.data(sih1 + 136);
    const auto *sih1_138 = buffer.data(sih1 + 138);
    const auto *sih1_140 = buffer.data(sih1 + 140);
    const auto *sih1_141 = buffer.data(sih1 + 141);
    const auto *sih1_143 = buffer.data(sih1 + 143);
    const auto *sih1_144 = buffer.data(sih1 + 144);

    const auto *sii_94 = buffer.data(sii + 94);
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
    const auto *sii_112 = buffer.data(sii + 112);
    const auto *sii_114 = buffer.data(sii + 114);
    const auto *sii_115 = buffer.data(sii + 115);
    const auto *sii_117 = buffer.data(sii + 117);
    const auto *sii_118 = buffer.data(sii + 118);
    const auto *sii_121 = buffer.data(sii + 121);
    const auto *sii_122 = buffer.data(sii + 122);
    const auto *sii_126 = buffer.data(sii + 126);
    const auto *sii_133 = buffer.data(sii + 133);
    const auto *sii_134 = buffer.data(sii + 134);
    const auto *sii_135 = buffer.data(sii + 135);
    const auto *sii_136 = buffer.data(sii + 136);
    const auto *sii_137 = buffer.data(sii + 137);
    const auto *sii_138 = buffer.data(sii + 138);
    const auto *sii_139 = buffer.data(sii + 139);
    const auto *sii_140 = buffer.data(sii + 140);
    const auto *sii_142 = buffer.data(sii + 142);
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
    const auto *sii_170 = buffer.data(sii + 170);
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

#pragma omp simd aligned(t_122, t_123, t_124, pc_x, pc_z, shi_98, shi_99, sih0_77, sih0_78, \
                         sih1_77, sih1_78, sii_94, sii_98, sii_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_16 * shi_98[k]
                   + f_8 * sih0_77[k]
                   - f_9 * sih1_77[k]
                   + f_3 * pc_x[k] * sii_98[k];

        t_123[k] = f_16 * shi_99[k]
                   + f_10 * sih0_78[k]
                   - f_11 * sih1_78[k]
                   + f_3 * pc_x[k] * sii_99[k];

        t_124[k] = f_3 * pc_z[k] * sii_94[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, pc_x, pc_y, shi_42, shi_101, shi_102, sih0_80, \
                         sih0_81, sih1_80, sih1_81, sii_98, sii_101, \
                         sii_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_16 * shi_101[k]
                   + f_10 * sih0_80[k]
                   - f_11 * sih1_80[k]
                   + f_3 * pc_x[k] * sii_101[k];

        t_126[k] = f_16 * shi_102[k]
                   + f_10 * sih0_81[k]
                   - f_11 * sih1_81[k]
                   + f_3 * pc_x[k] * sii_102[k];

        t_127[k] = f_14 * shi_42[k]
                   + f_3 * pc_y[k] * sii_98[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, pc_x, shi_104, shi_105, shi_106, shi_107, \
                         sih0_83, sih1_83, sii_104, sii_105, sii_106, \
                         sii_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_16 * shi_104[k]
                   + f_10 * sih0_83[k]
                   - f_11 * sih1_83[k]
                   + f_3 * pc_x[k] * sii_104[k];

        t_129[k] = f_16 * shi_105[k]
                   + f_3 * pc_x[k] * sii_105[k];

        t_130[k] = f_16 * shi_106[k]
                   + f_3 * pc_x[k] * sii_106[k];

        t_131[k] = f_16 * shi_107[k]
                   + f_3 * pc_x[k] * sii_107[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pc_x, shi_108, shi_109, shi_110, shi_111, \
                         sii_108, sii_109, sii_110, sii_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_16 * shi_108[k]
                   + f_3 * pc_x[k] * sii_108[k];

        t_133[k] = f_16 * shi_109[k]
                   + f_3 * pc_x[k] * sii_109[k];

        t_134[k] = f_16 * shi_110[k]
                   + f_3 * pc_x[k] * sii_110[k];

        t_135[k] = f_16 * shi_111[k]
                   + f_3 * pc_x[k] * sii_111[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, pc_y, pc_z, shi_49, shi_51, sih0_78, sih0_80, \
                         sih1_78, sih1_80, sii_105, sii_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_14 * shi_49[k]
                   + f_1 * sih0_78[k]
                   - f_2 * sih1_78[k]
                   + f_3 * pc_y[k] * sii_105[k];

        t_137[k] = f_3 * pc_z[k] * sii_105[k];

        t_138[k] = f_14 * shi_51[k]
                   + f_4 * sih0_80[k]
                   - f_5 * sih1_80[k]
                   + f_3 * pc_y[k] * sii_107[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, pc_y, shi_52, shi_53, shi_54, sih0_81, sih0_82, \
                         sih0_83, sih1_81, sih1_82, sih1_83, sii_108, sii_109, \
                         sii_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_14 * shi_52[k]
                   + f_6 * sih0_81[k]
                   - f_7 * sih1_81[k]
                   + f_3 * pc_y[k] * sii_108[k];

        t_140[k] = f_14 * shi_53[k]
                   + f_8 * sih0_82[k]
                   - f_9 * sih1_82[k]
                   + f_3 * pc_y[k] * sii_109[k];

        t_141[k] = f_14 * shi_54[k]
                   + f_10 * sih0_83[k]
                   - f_11 * sih1_83[k]
                   + f_3 * pc_y[k] * sii_110[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, pb_y, pc_y, pc_z, shk0_72, shi_55, \
                         shi_56, shk1_72, sih0_83, sih1_83, sii_111, \
                         sii_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_14 * shi_55[k]
                   + f_3 * pc_y[k] * sii_111[k];

        t_143[k] = f_1 * sih0_83[k]
                   - f_2 * sih1_83[k]
                   + f_3 * pc_z[k] * sii_111[k];

        t_144[k] = pb_y[k] * shk0_72[k]
                   - f_12 * pc_y[k] * shk1_72[k];

        t_145[k] = f_13 * shi_56[k]
                   + f_3 * pc_y[k] * sii_112[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pb_y, pb_z, pc_y, pc_z, shk0_39, shk0_77, \
                         shi_28, shi_58, shk1_39, shk1_77, sii_112, \
                         sii_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_13 * shi_28[k]
                   + f_3 * pc_z[k] * sii_112[k];

        t_147[k] = pb_z[k] * shk0_39[k]
                   - f_12 * pc_z[k] * shk1_39[k];

        t_148[k] = f_13 * shi_58[k]
                   + f_3 * pc_y[k] * sii_114[k];

        t_149[k] = pb_y[k] * shk0_77[k]
                   - f_12 * pc_y[k] * shk1_77[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, pb_y, pb_z, pc_y, pc_z, shk0_42, shk0_81, \
                         shi_31, shi_61, shk1_42, shk1_81, sii_115, \
                         sii_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = pb_z[k] * shk0_42[k]
                   - f_12 * pc_z[k] * shk1_42[k];

        t_151[k] = f_13 * shi_31[k]
                   + f_3 * pc_z[k] * sii_115[k];

        t_152[k] = f_13 * shi_61[k]
                   + f_3 * pc_y[k] * sii_117[k];

        t_153[k] = pb_y[k] * shk0_81[k]
                   - f_12 * pc_y[k] * shk1_81[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, pb_y, pb_z, pc_y, pc_z, shk0_46, shk0_84, \
                         shi_34, shi_64, shk1_46, shk1_84, sii_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = pb_z[k] * shk0_46[k]
                   - f_12 * pc_z[k] * shk1_46[k];

        t_155[k] = f_13 * shi_34[k]
                   + f_3 * pc_z[k] * sii_118[k];

        t_156[k] = pb_y[k] * shk0_84[k]
                   + f_14 * shi_64[k]
                   - f_12 * pc_y[k] * shk1_84[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, pb_y, pb_z, pc_y, pc_z, shk0_51, shk0_86, \
                         shi_38, shi_65, shk1_51, shk1_86, sii_121, \
                         sii_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = f_13 * shi_65[k]
                   + f_3 * pc_y[k] * sii_121[k];

        t_158[k] = pb_y[k] * shk0_86[k]
                   - f_12 * pc_y[k] * shk1_86[k];

        t_159[k] = pb_z[k] * shk0_51[k]
                   - f_12 * pc_z[k] * shk1_51[k];

        t_160[k] = f_13 * shi_38[k]
                   + f_3 * pc_z[k] * sii_122[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, t_164, pb_y, pc_y, shk0_89, shk0_90, shk0_92, \
                         shi_68, shi_69, shi_70, shk1_89, shk1_90, shk1_92, \
                         sii_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = pb_y[k] * shk0_89[k]
                   + f_15 * shi_68[k]
                   - f_12 * pc_y[k] * shk1_89[k];

        t_162[k] = pb_y[k] * shk0_90[k]
                   + f_14 * shi_69[k]
                   - f_12 * pc_y[k] * shk1_90[k];

        t_163[k] = f_13 * shi_70[k]
                   + f_3 * pc_y[k] * sii_126[k];

        t_164[k] = pb_y[k] * shk0_92[k]
                   - f_12 * pc_y[k] * shk1_92[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, pc_x, shi_133, shi_134, shi_135, \
                         shi_136, shi_137, sii_133, sii_134, sii_135, sii_136, \
                         sii_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = f_16 * shi_133[k]
                   + f_3 * pc_x[k] * sii_133[k];

        t_166[k] = f_16 * shi_134[k]
                   + f_3 * pc_x[k] * sii_134[k];

        t_167[k] = f_16 * shi_135[k]
                   + f_3 * pc_x[k] * sii_135[k];

        t_168[k] = f_16 * shi_136[k]
                   + f_3 * pc_x[k] * sii_136[k];

        t_169[k] = f_16 * shi_137[k]
                   + f_3 * pc_x[k] * sii_137[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, pb_z, pc_x, pc_z, shk0_64, shi_49, \
                         shi_138, shi_139, shk1_64, sii_133, sii_138, \
                         sii_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_16 * shi_138[k]
                   + f_3 * pc_x[k] * sii_138[k];

        t_171[k] = f_16 * shi_139[k]
                   + f_3 * pc_x[k] * sii_139[k];

        t_172[k] = pb_z[k] * shk0_64[k]
                   - f_12 * pc_z[k] * shk1_64[k];

        t_173[k] = f_13 * shi_49[k]
                   + f_3 * pc_z[k] * sii_133[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, pc_y, shi_79, shi_80, shi_81, sih0_101, \
                         sih0_102, sih0_103, sih1_101, sih1_102, sih1_103, sii_135, sii_136, \
                         sii_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_13 * shi_79[k]
                   + f_4 * sih0_101[k]
                   - f_5 * sih1_101[k]
                   + f_3 * pc_y[k] * sii_135[k];

        t_175[k] = f_13 * shi_80[k]
                   + f_6 * sih0_102[k]
                   - f_7 * sih1_102[k]
                   + f_3 * pc_y[k] * sii_136[k];

        t_176[k] = f_13 * shi_81[k]
                   + f_8 * sih0_103[k]
                   - f_9 * sih1_103[k]
                   + f_3 * pc_y[k] * sii_137[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, pb_y, pc_y, shk0_107, shi_82, shi_83, shk1_107, \
                         sih0_104, sih1_104, sii_138, sii_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = f_13 * shi_82[k]
                   + f_10 * sih0_104[k]
                   - f_11 * sih1_104[k]
                   + f_3 * pc_y[k] * sii_138[k];

        t_178[k] = f_13 * shi_83[k]
                   + f_3 * pc_y[k] * sii_139[k];

        t_179[k] = pb_y[k] * shk0_107[k]
                   - f_12 * pc_y[k] * shk1_107[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pc_x, pc_y, pc_z, shi_56, shi_140, \
                         shi_143, sih0_105, sih0_108, sih1_105, sih1_108, sii_140, \
                         sii_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_16 * shi_140[k]
                   + f_1 * sih0_105[k]
                   - f_2 * sih1_105[k]
                   + f_3 * pc_x[k] * sii_140[k];

        t_181[k] = f_3 * pc_y[k] * sii_140[k];

        t_182[k] = f_14 * shi_56[k]
                   + f_3 * pc_z[k] * sii_140[k];

        t_183[k] = f_16 * shi_143[k]
                   + f_4 * sih0_108[k]
                   - f_5 * sih1_108[k]
                   + f_3 * pc_x[k] * sii_143[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, pc_x, pc_y, shi_145, shi_146, sih0_110, \
                         sih0_111, sih1_110, sih1_111, sii_142, sii_145, \
                         sii_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_3 * pc_y[k] * sii_142[k];

        t_185[k] = f_16 * shi_145[k]
                   + f_4 * sih0_110[k]
                   - f_5 * sih1_110[k]
                   + f_3 * pc_x[k] * sii_145[k];

        t_186[k] = f_16 * shi_146[k]
                   + f_6 * sih0_111[k]
                   - f_7 * sih1_111[k]
                   + f_3 * pc_x[k] * sii_146[k];
    }

#pragma omp simd aligned(t_187, t_188, t_189, pc_x, pc_y, pc_z, shi_59, shi_149, sih0_114, \
                         sih1_114, sii_143, sii_145, sii_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_187[k] = f_14 * shi_59[k]
                   + f_3 * pc_z[k] * sii_143[k];

        t_188[k] = f_3 * pc_y[k] * sii_145[k];

        t_189[k] = f_16 * shi_149[k]
                   + f_6 * sih0_114[k]
                   - f_7 * sih1_114[k]
                   + f_3 * pc_x[k] * sii_149[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, pc_x, pc_z, shi_62, shi_150, shi_152, sih0_115, \
                         sih0_117, sih1_115, sih1_117, sii_146, sii_150, \
                         sii_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = f_16 * shi_150[k]
                   + f_8 * sih0_115[k]
                   - f_9 * sih1_115[k]
                   + f_3 * pc_x[k] * sii_150[k];

        t_191[k] = f_14 * shi_62[k]
                   + f_3 * pc_z[k] * sii_146[k];

        t_192[k] = f_16 * shi_152[k]
                   + f_8 * sih0_117[k]
                   - f_9 * sih1_117[k]
                   + f_3 * pc_x[k] * sii_152[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, pc_x, pc_y, shi_154, shi_155, sih0_119, \
                         sih0_120, sih1_119, sih1_120, sii_149, sii_154, \
                         sii_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_3 * pc_y[k] * sii_149[k];

        t_194[k] = f_16 * shi_154[k]
                   + f_8 * sih0_119[k]
                   - f_9 * sih1_119[k]
                   + f_3 * pc_x[k] * sii_154[k];

        t_195[k] = f_16 * shi_155[k]
                   + f_10 * sih0_120[k]
                   - f_11 * sih1_120[k]
                   + f_3 * pc_x[k] * sii_155[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, pc_x, pc_z, shi_66, shi_157, shi_158, sih0_122, \
                         sih0_123, sih1_122, sih1_123, sii_150, sii_157, \
                         sii_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = f_14 * shi_66[k]
                   + f_3 * pc_z[k] * sii_150[k];

        t_197[k] = f_16 * shi_157[k]
                   + f_10 * sih0_122[k]
                   - f_11 * sih1_122[k]
                   + f_3 * pc_x[k] * sii_157[k];

        t_198[k] = f_16 * shi_158[k]
                   + f_10 * sih0_123[k]
                   - f_11 * sih1_123[k]
                   + f_3 * pc_x[k] * sii_158[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, pc_x, pc_y, shi_160, shi_161, shi_162, \
                         sih0_125, sih1_125, sii_154, sii_160, sii_161, \
                         sii_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = f_3 * pc_y[k] * sii_154[k];

        t_200[k] = f_16 * shi_160[k]
                   + f_10 * sih0_125[k]
                   - f_11 * sih1_125[k]
                   + f_3 * pc_x[k] * sii_160[k];

        t_201[k] = f_16 * shi_161[k]
                   + f_3 * pc_x[k] * sii_161[k];

        t_202[k] = f_16 * shi_162[k]
                   + f_3 * pc_x[k] * sii_162[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, t_207, pc_x, shi_163, shi_164, shi_165, \
                         shi_166, shi_167, sii_163, sii_164, sii_165, sii_166, \
                         sii_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_16 * shi_163[k]
                   + f_3 * pc_x[k] * sii_163[k];

        t_204[k] = f_16 * shi_164[k]
                   + f_3 * pc_x[k] * sii_164[k];

        t_205[k] = f_16 * shi_165[k]
                   + f_3 * pc_x[k] * sii_165[k];

        t_206[k] = f_16 * shi_166[k]
                   + f_3 * pc_x[k] * sii_166[k];

        t_207[k] = f_16 * shi_167[k]
                   + f_3 * pc_x[k] * sii_167[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, pc_y, pc_z, shi_77, sih0_120, sih0_122, \
                         sih0_123, sih1_120, sih1_122, sih1_123, sii_161, sii_163, \
                         sii_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_1 * sih0_120[k]
                   - f_2 * sih1_120[k]
                   + f_3 * pc_y[k] * sii_161[k];

        t_209[k] = f_14 * shi_77[k]
                   + f_3 * pc_z[k] * sii_161[k];

        t_210[k] = f_4 * sih0_122[k]
                   - f_5 * sih1_122[k]
                   + f_3 * pc_y[k] * sii_163[k];

        t_211[k] = f_6 * sih0_123[k]
                   - f_7 * sih1_123[k]
                   + f_3 * pc_y[k] * sii_164[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, pc_y, pc_z, shi_83, sih0_124, sih0_125, \
                         sih1_124, sih1_125, sii_165, sii_166, \
                         sii_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_8 * sih0_124[k]
                   - f_9 * sih1_124[k]
                   + f_3 * pc_y[k] * sii_165[k];

        t_213[k] = f_10 * sih0_125[k]
                   - f_11 * sih1_125[k]
                   + f_3 * pc_y[k] * sii_166[k];

        t_214[k] = f_3 * pc_y[k] * sii_167[k];

        t_215[k] = f_14 * shi_83[k]
                   + f_1 * sih0_125[k]
                   - f_2 * sih1_125[k]
                   + f_3 * pc_z[k] * sii_167[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, pc_x, pc_y, pc_z, shi_84, shi_168, \
                         shi_171, sih0_126, sih0_129, sih1_126, sih1_129, sii_168, \
                         sii_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = f_15 * shi_168[k]
                   + f_1 * sih0_126[k]
                   - f_2 * sih1_126[k]
                   + f_3 * pc_x[k] * sii_168[k];

        t_217[k] = f_15 * shi_84[k]
                   + f_3 * pc_y[k] * sii_168[k];

        t_218[k] = f_3 * pc_z[k] * sii_168[k];

        t_219[k] = f_15 * shi_171[k]
                   + f_4 * sih0_129[k]
                   - f_5 * sih1_129[k]
                   + f_3 * pc_x[k] * sii_171[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, pc_x, pc_y, shi_86, shi_173, shi_174, sih0_131, \
                         sih0_132, sih1_131, sih1_132, sii_170, sii_173, \
                         sii_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = f_15 * shi_86[k]
                   + f_3 * pc_y[k] * sii_170[k];

        t_221[k] = f_15 * shi_173[k]
                   + f_4 * sih0_131[k]
                   - f_5 * sih1_131[k]
                   + f_3 * pc_x[k] * sii_173[k];

        t_222[k] = f_15 * shi_174[k]
                   + f_6 * sih0_132[k]
                   - f_7 * sih1_132[k]
                   + f_3 * pc_x[k] * sii_174[k];
    }

#pragma omp simd aligned(t_223, t_224, t_225, pc_x, pc_y, pc_z, shi_89, shi_177, sih0_135, \
                         sih1_135, sii_171, sii_173, sii_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = f_3 * pc_z[k] * sii_171[k];

        t_224[k] = f_15 * shi_89[k]
                   + f_3 * pc_y[k] * sii_173[k];

        t_225[k] = f_15 * shi_177[k]
                   + f_6 * sih0_135[k]
                   - f_7 * sih1_135[k]
                   + f_3 * pc_x[k] * sii_177[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, pc_x, pc_z, shi_178, shi_180, sih0_136, \
                         sih0_138, sih1_136, sih1_138, sii_174, sii_178, \
                         sii_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = f_15 * shi_178[k]
                   + f_8 * sih0_136[k]
                   - f_9 * sih1_136[k]
                   + f_3 * pc_x[k] * sii_178[k];

        t_227[k] = f_3 * pc_z[k] * sii_174[k];

        t_228[k] = f_15 * shi_180[k]
                   + f_8 * sih0_138[k]
                   - f_9 * sih1_138[k]
                   + f_3 * pc_x[k] * sii_180[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, pc_x, pc_y, shi_93, shi_182, shi_183, sih0_140, \
                         sih0_141, sih1_140, sih1_141, sii_177, sii_182, \
                         sii_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = f_15 * shi_93[k]
                   + f_3 * pc_y[k] * sii_177[k];

        t_230[k] = f_15 * shi_182[k]
                   + f_8 * sih0_140[k]
                   - f_9 * sih1_140[k]
                   + f_3 * pc_x[k] * sii_182[k];

        t_231[k] = f_15 * shi_183[k]
                   + f_10 * sih0_141[k]
                   - f_11 * sih1_141[k]
                   + f_3 * pc_x[k] * sii_183[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, pc_x, pc_z, shi_185, shi_186, sih0_143, \
                         sih0_144, sih1_143, sih1_144, sii_178, sii_185, \
                         sii_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = f_3 * pc_z[k] * sii_178[k];

        t_233[k] = f_15 * shi_185[k]
                   + f_10 * sih0_143[k]
                   - f_11 * sih1_143[k]
                   + f_3 * pc_x[k] * sii_185[k];

        t_234[k] = f_15 * shi_186[k]
                   + f_10 * sih0_144[k]
                   - f_11 * sih1_144[k]
                   + f_3 * pc_x[k] * sii_186[k];
    }
}

static auto
compute_prim_sik_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t shk0,
                                                          const size_t shi, const size_t shk1,
                                                          const size_t sih0, const size_t sih1,
                                                          const size_t sii, const size_t ncols,
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

    const auto *shk0_108 = buffer.data(shk0 + 108);
    const auto *shk0_111 = buffer.data(shk0 + 111);
    const auto *shk0_114 = buffer.data(shk0 + 114);
    const auto *shk0_118 = buffer.data(shk0 + 118);
    const auto *shk0_120 = buffer.data(shk0 + 120);
    const auto *shk0_123 = buffer.data(shk0 + 123);
    const auto *shk0_125 = buffer.data(shk0 + 125);
    const auto *shk0_126 = buffer.data(shk0 + 126);
    const auto *shk0_136 = buffer.data(shk0 + 136);
    const auto *shk0_180 = buffer.data(shk0 + 180);
    const auto *shk0_183 = buffer.data(shk0 + 183);
    const auto *shk0_185 = buffer.data(shk0 + 185);
    const auto *shk0_186 = buffer.data(shk0 + 186);
    const auto *shk0_189 = buffer.data(shk0 + 189);
    const auto *shk0_190 = buffer.data(shk0 + 190);
    const auto *shk0_192 = buffer.data(shk0 + 192);
    const auto *shk0_194 = buffer.data(shk0 + 194);
    const auto *shk0_195 = buffer.data(shk0 + 195);
    const auto *shk0_197 = buffer.data(shk0 + 197);
    const auto *shk0_198 = buffer.data(shk0 + 198);
    const auto *shk0_200 = buffer.data(shk0 + 200);
    const auto *shk0_215 = buffer.data(shk0 + 215);

    const auto *shi_84 = buffer.data(shi + 84);
    const auto *shi_87 = buffer.data(shi + 87);
    const auto *shi_90 = buffer.data(shi + 90);
    const auto *shi_91 = buffer.data(shi + 91);
    const auto *shi_94 = buffer.data(shi + 94);
    const auto *shi_95 = buffer.data(shi + 95);
    const auto *shi_96 = buffer.data(shi + 96);
    const auto *shi_98 = buffer.data(shi + 98);
    const auto *shi_105 = buffer.data(shi + 105);
    const auto *shi_107 = buffer.data(shi + 107);
    const auto *shi_108 = buffer.data(shi + 108);
    const auto *shi_109 = buffer.data(shi + 109);
    const auto *shi_110 = buffer.data(shi + 110);
    const auto *shi_111 = buffer.data(shi + 111);
    const auto *shi_112 = buffer.data(shi + 112);
    const auto *shi_114 = buffer.data(shi + 114);
    const auto *shi_115 = buffer.data(shi + 115);
    const auto *shi_117 = buffer.data(shi + 117);
    const auto *shi_118 = buffer.data(shi + 118);
    const auto *shi_121 = buffer.data(shi + 121);
    const auto *shi_122 = buffer.data(shi + 122);
    const auto *shi_126 = buffer.data(shi + 126);
    const auto *shi_133 = buffer.data(shi + 133);
    const auto *shi_135 = buffer.data(shi + 135);
    const auto *shi_136 = buffer.data(shi + 136);
    const auto *shi_137 = buffer.data(shi + 137);
    const auto *shi_138 = buffer.data(shi + 138);
    const auto *shi_139 = buffer.data(shi + 139);
    const auto *shi_140 = buffer.data(shi + 140);
    const auto *shi_141 = buffer.data(shi + 141);
    const auto *shi_142 = buffer.data(shi + 142);
    const auto *shi_143 = buffer.data(shi + 143);
    const auto *shi_145 = buffer.data(shi + 145);
    const auto *shi_146 = buffer.data(shi + 146);
    const auto *shi_148 = buffer.data(shi + 148);
    const auto *shi_149 = buffer.data(shi + 149);
    const auto *shi_150 = buffer.data(shi + 150);
    const auto *shi_152 = buffer.data(shi + 152);
    const auto *shi_153 = buffer.data(shi + 153);
    const auto *shi_154 = buffer.data(shi + 154);
    const auto *shi_161 = buffer.data(shi + 161);
    const auto *shi_163 = buffer.data(shi + 163);
    const auto *shi_164 = buffer.data(shi + 164);
    const auto *shi_165 = buffer.data(shi + 165);
    const auto *shi_166 = buffer.data(shi + 166);
    const auto *shi_167 = buffer.data(shi + 167);
    const auto *shi_188 = buffer.data(shi + 188);
    const auto *shi_189 = buffer.data(shi + 189);
    const auto *shi_190 = buffer.data(shi + 190);
    const auto *shi_191 = buffer.data(shi + 191);
    const auto *shi_192 = buffer.data(shi + 192);
    const auto *shi_193 = buffer.data(shi + 193);
    const auto *shi_194 = buffer.data(shi + 194);
    const auto *shi_195 = buffer.data(shi + 195);
    const auto *shi_201 = buffer.data(shi + 201);
    const auto *shi_205 = buffer.data(shi + 205);
    const auto *shi_210 = buffer.data(shi + 210);
    const auto *shi_216 = buffer.data(shi + 216);
    const auto *shi_217 = buffer.data(shi + 217);
    const auto *shi_218 = buffer.data(shi + 218);
    const auto *shi_219 = buffer.data(shi + 219);
    const auto *shi_220 = buffer.data(shi + 220);
    const auto *shi_221 = buffer.data(shi + 221);
    const auto *shi_222 = buffer.data(shi + 222);
    const auto *shi_223 = buffer.data(shi + 223);
    const auto *shi_245 = buffer.data(shi + 245);
    const auto *shi_246 = buffer.data(shi + 246);
    const auto *shi_247 = buffer.data(shi + 247);
    const auto *shi_248 = buffer.data(shi + 248);
    const auto *shi_249 = buffer.data(shi + 249);
    const auto *shi_250 = buffer.data(shi + 250);
    const auto *shi_251 = buffer.data(shi + 251);
    const auto *shi_252 = buffer.data(shi + 252);
    const auto *shi_255 = buffer.data(shi + 255);
    const auto *shi_257 = buffer.data(shi + 257);
    const auto *shi_258 = buffer.data(shi + 258);
    const auto *shi_261 = buffer.data(shi + 261);
    const auto *shi_262 = buffer.data(shi + 262);
    const auto *shi_264 = buffer.data(shi + 264);
    const auto *shi_266 = buffer.data(shi + 266);
    const auto *shi_267 = buffer.data(shi + 267);
    const auto *shi_269 = buffer.data(shi + 269);
    const auto *shi_270 = buffer.data(shi + 270);
    const auto *shi_272 = buffer.data(shi + 272);
    const auto *shi_273 = buffer.data(shi + 273);
    const auto *shi_274 = buffer.data(shi + 274);

    const auto *shk1_108 = buffer.data(shk1 + 108);
    const auto *shk1_111 = buffer.data(shk1 + 111);
    const auto *shk1_114 = buffer.data(shk1 + 114);
    const auto *shk1_118 = buffer.data(shk1 + 118);
    const auto *shk1_120 = buffer.data(shk1 + 120);
    const auto *shk1_123 = buffer.data(shk1 + 123);
    const auto *shk1_125 = buffer.data(shk1 + 125);
    const auto *shk1_126 = buffer.data(shk1 + 126);
    const auto *shk1_136 = buffer.data(shk1 + 136);
    const auto *shk1_180 = buffer.data(shk1 + 180);
    const auto *shk1_183 = buffer.data(shk1 + 183);
    const auto *shk1_185 = buffer.data(shk1 + 185);
    const auto *shk1_186 = buffer.data(shk1 + 186);
    const auto *shk1_189 = buffer.data(shk1 + 189);
    const auto *shk1_190 = buffer.data(shk1 + 190);
    const auto *shk1_192 = buffer.data(shk1 + 192);
    const auto *shk1_194 = buffer.data(shk1 + 194);
    const auto *shk1_195 = buffer.data(shk1 + 195);
    const auto *shk1_197 = buffer.data(shk1 + 197);
    const auto *shk1_198 = buffer.data(shk1 + 198);
    const auto *shk1_200 = buffer.data(shk1 + 200);
    const auto *shk1_215 = buffer.data(shk1 + 215);

    const auto *sih0_141 = buffer.data(sih0 + 141);
    const auto *sih0_143 = buffer.data(sih0 + 143);
    const auto *sih0_144 = buffer.data(sih0 + 144);
    const auto *sih0_145 = buffer.data(sih0 + 145);
    const auto *sih0_146 = buffer.data(sih0 + 146);
    const auto *sih0_152 = buffer.data(sih0 + 152);
    const auto *sih0_156 = buffer.data(sih0 + 156);
    const auto *sih0_161 = buffer.data(sih0 + 161);
    const auto *sih0_164 = buffer.data(sih0 + 164);
    const auto *sih0_165 = buffer.data(sih0 + 165);
    const auto *sih0_166 = buffer.data(sih0 + 166);
    const auto *sih0_167 = buffer.data(sih0 + 167);
    const auto *sih0_183 = buffer.data(sih0 + 183);
    const auto *sih0_185 = buffer.data(sih0 + 185);
    const auto *sih0_186 = buffer.data(sih0 + 186);
    const auto *sih0_187 = buffer.data(sih0 + 187);
    const auto *sih0_188 = buffer.data(sih0 + 188);
    const auto *sih0_189 = buffer.data(sih0 + 189);
    const auto *sih0_192 = buffer.data(sih0 + 192);
    const auto *sih0_194 = buffer.data(sih0 + 194);
    const auto *sih0_195 = buffer.data(sih0 + 195);
    const auto *sih0_198 = buffer.data(sih0 + 198);
    const auto *sih0_199 = buffer.data(sih0 + 199);
    const auto *sih0_201 = buffer.data(sih0 + 201);
    const auto *sih0_203 = buffer.data(sih0 + 203);
    const auto *sih0_204 = buffer.data(sih0 + 204);
    const auto *sih0_206 = buffer.data(sih0 + 206);
    const auto *sih0_207 = buffer.data(sih0 + 207);
    const auto *sih0_209 = buffer.data(sih0 + 209);

    const auto *sih1_141 = buffer.data(sih1 + 141);
    const auto *sih1_143 = buffer.data(sih1 + 143);
    const auto *sih1_144 = buffer.data(sih1 + 144);
    const auto *sih1_145 = buffer.data(sih1 + 145);
    const auto *sih1_146 = buffer.data(sih1 + 146);
    const auto *sih1_152 = buffer.data(sih1 + 152);
    const auto *sih1_156 = buffer.data(sih1 + 156);
    const auto *sih1_161 = buffer.data(sih1 + 161);
    const auto *sih1_164 = buffer.data(sih1 + 164);
    const auto *sih1_165 = buffer.data(sih1 + 165);
    const auto *sih1_166 = buffer.data(sih1 + 166);
    const auto *sih1_167 = buffer.data(sih1 + 167);
    const auto *sih1_183 = buffer.data(sih1 + 183);
    const auto *sih1_185 = buffer.data(sih1 + 185);
    const auto *sih1_186 = buffer.data(sih1 + 186);
    const auto *sih1_187 = buffer.data(sih1 + 187);
    const auto *sih1_188 = buffer.data(sih1 + 188);
    const auto *sih1_189 = buffer.data(sih1 + 189);
    const auto *sih1_192 = buffer.data(sih1 + 192);
    const auto *sih1_194 = buffer.data(sih1 + 194);
    const auto *sih1_195 = buffer.data(sih1 + 195);
    const auto *sih1_198 = buffer.data(sih1 + 198);
    const auto *sih1_199 = buffer.data(sih1 + 199);
    const auto *sih1_201 = buffer.data(sih1 + 201);
    const auto *sih1_203 = buffer.data(sih1 + 203);
    const auto *sih1_204 = buffer.data(sih1 + 204);
    const auto *sih1_206 = buffer.data(sih1 + 206);
    const auto *sih1_207 = buffer.data(sih1 + 207);
    const auto *sih1_209 = buffer.data(sih1 + 209);

    const auto *sii_182 = buffer.data(sii + 182);
    const auto *sii_188 = buffer.data(sii + 188);
    const auto *sii_189 = buffer.data(sii + 189);
    const auto *sii_190 = buffer.data(sii + 190);
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
    const auto *sii_216 = buffer.data(sii + 216);
    const auto *sii_217 = buffer.data(sii + 217);
    const auto *sii_218 = buffer.data(sii + 218);
    const auto *sii_219 = buffer.data(sii + 219);
    const auto *sii_220 = buffer.data(sii + 220);
    const auto *sii_221 = buffer.data(sii + 221);
    const auto *sii_222 = buffer.data(sii + 222);
    const auto *sii_223 = buffer.data(sii + 223);
    const auto *sii_224 = buffer.data(sii + 224);
    const auto *sii_226 = buffer.data(sii + 226);
    const auto *sii_227 = buffer.data(sii + 227);
    const auto *sii_229 = buffer.data(sii + 229);
    const auto *sii_230 = buffer.data(sii + 230);
    const auto *sii_233 = buffer.data(sii + 233);
    const auto *sii_234 = buffer.data(sii + 234);
    const auto *sii_238 = buffer.data(sii + 238);
    const auto *sii_245 = buffer.data(sii + 245);
    const auto *sii_246 = buffer.data(sii + 246);
    const auto *sii_247 = buffer.data(sii + 247);
    const auto *sii_248 = buffer.data(sii + 248);
    const auto *sii_249 = buffer.data(sii + 249);
    const auto *sii_250 = buffer.data(sii + 250);
    const auto *sii_251 = buffer.data(sii + 251);
    const auto *sii_252 = buffer.data(sii + 252);
    const auto *sii_254 = buffer.data(sii + 254);
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

#pragma omp simd aligned(t_235, t_236, t_237, t_238, pc_x, pc_y, shi_98, shi_188, shi_189, \
                         shi_190, sih0_146, sih1_146, sii_182, sii_188, sii_189, \
                         sii_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = f_15 * shi_98[k]
                   + f_3 * pc_y[k] * sii_182[k];

        t_236[k] = f_15 * shi_188[k]
                   + f_10 * sih0_146[k]
                   - f_11 * sih1_146[k]
                   + f_3 * pc_x[k] * sii_188[k];

        t_237[k] = f_15 * shi_189[k]
                   + f_3 * pc_x[k] * sii_189[k];

        t_238[k] = f_15 * shi_190[k]
                   + f_3 * pc_x[k] * sii_190[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, t_243, pc_x, shi_191, shi_192, shi_193, \
                         shi_194, shi_195, sii_191, sii_192, sii_193, sii_194, \
                         sii_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_15 * shi_191[k]
                   + f_3 * pc_x[k] * sii_191[k];

        t_240[k] = f_15 * shi_192[k]
                   + f_3 * pc_x[k] * sii_192[k];

        t_241[k] = f_15 * shi_193[k]
                   + f_3 * pc_x[k] * sii_193[k];

        t_242[k] = f_15 * shi_194[k]
                   + f_3 * pc_x[k] * sii_194[k];

        t_243[k] = f_15 * shi_195[k]
                   + f_3 * pc_x[k] * sii_195[k];
    }

#pragma omp simd aligned(t_244, t_245, t_246, pc_y, pc_z, shi_105, shi_107, sih0_141, \
                         sih0_143, sih1_141, sih1_143, sii_189, \
                         sii_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_244[k] = f_15 * shi_105[k]
                   + f_1 * sih0_141[k]
                   - f_2 * sih1_141[k]
                   + f_3 * pc_y[k] * sii_189[k];

        t_245[k] = f_3 * pc_z[k] * sii_189[k];

        t_246[k] = f_15 * shi_107[k]
                   + f_4 * sih0_143[k]
                   - f_5 * sih1_143[k]
                   + f_3 * pc_y[k] * sii_191[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, pc_y, shi_108, shi_109, shi_110, sih0_144, \
                         sih0_145, sih0_146, sih1_144, sih1_145, sih1_146, sii_192, sii_193, \
                         sii_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = f_15 * shi_108[k]
                   + f_6 * sih0_144[k]
                   - f_7 * sih1_144[k]
                   + f_3 * pc_y[k] * sii_192[k];

        t_248[k] = f_15 * shi_109[k]
                   + f_8 * sih0_145[k]
                   - f_9 * sih1_145[k]
                   + f_3 * pc_y[k] * sii_193[k];

        t_249[k] = f_15 * shi_110[k]
                   + f_10 * sih0_146[k]
                   - f_11 * sih1_146[k]
                   + f_3 * pc_y[k] * sii_194[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, pb_z, pc_y, pc_z, shk0_108, shi_111, \
                         shi_112, shk1_108, sih0_146, sih1_146, sii_195, \
                         sii_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = f_15 * shi_111[k]
                   + f_3 * pc_y[k] * sii_195[k];

        t_251[k] = f_1 * sih0_146[k]
                   - f_2 * sih1_146[k]
                   + f_3 * pc_z[k] * sii_195[k];

        t_252[k] = pb_z[k] * shk0_108[k]
                   - f_12 * pc_z[k] * shk1_108[k];

        t_253[k] = f_14 * shi_112[k]
                   + f_3 * pc_y[k] * sii_196[k];
    }

#pragma omp simd aligned(t_254, t_255, t_256, pb_z, pc_y, pc_z, shk0_111, shi_84, shi_114, \
                         shk1_111, sii_196, sii_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_254[k] = f_13 * shi_84[k]
                   + f_3 * pc_z[k] * sii_196[k];

        t_255[k] = pb_z[k] * shk0_111[k]
                   - f_12 * pc_z[k] * shk1_111[k];

        t_256[k] = f_14 * shi_114[k]
                   + f_3 * pc_y[k] * sii_198[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, pb_z, pc_x, pc_z, shk0_114, shi_87, shi_201, \
                         shk1_114, sih0_152, sih1_152, sii_199, \
                         sii_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_15 * shi_201[k]
                   + f_4 * sih0_152[k]
                   - f_5 * sih1_152[k]
                   + f_3 * pc_x[k] * sii_201[k];

        t_258[k] = pb_z[k] * shk0_114[k]
                   - f_12 * pc_z[k] * shk1_114[k];

        t_259[k] = f_13 * shi_87[k]
                   + f_3 * pc_z[k] * sii_199[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, pb_z, pc_x, pc_y, pc_z, shk0_118, shi_117, \
                         shi_205, shk1_118, sih0_156, sih1_156, sii_201, \
                         sii_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_14 * shi_117[k]
                   + f_3 * pc_y[k] * sii_201[k];

        t_261[k] = f_15 * shi_205[k]
                   + f_6 * sih0_156[k]
                   - f_7 * sih1_156[k]
                   + f_3 * pc_x[k] * sii_205[k];

        t_262[k] = pb_z[k] * shk0_118[k]
                   - f_12 * pc_z[k] * shk1_118[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, pb_z, pc_y, pc_z, shk0_120, shi_90, shi_91, \
                         shi_121, shk1_120, sii_202, sii_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = f_13 * shi_90[k]
                   + f_3 * pc_z[k] * sii_202[k];

        t_264[k] = pb_z[k] * shk0_120[k]
                   + f_14 * shi_91[k]
                   - f_12 * pc_z[k] * shk1_120[k];

        t_265[k] = f_14 * shi_121[k]
                   + f_3 * pc_y[k] * sii_205[k];
    }

#pragma omp simd aligned(t_266, t_267, t_268, pb_z, pc_x, pc_z, shk0_123, shi_94, shi_210, \
                         shk1_123, sih0_161, sih1_161, sii_206, \
                         sii_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_266[k] = f_15 * shi_210[k]
                   + f_8 * sih0_161[k]
                   - f_9 * sih1_161[k]
                   + f_3 * pc_x[k] * sii_210[k];

        t_267[k] = pb_z[k] * shk0_123[k]
                   - f_12 * pc_z[k] * shk1_123[k];

        t_268[k] = f_13 * shi_94[k]
                   + f_3 * pc_z[k] * sii_206[k];
    }

#pragma omp simd aligned(t_269, t_270, t_271, pb_z, pc_y, pc_z, shk0_125, shk0_126, shi_95, \
                         shi_96, shi_126, shk1_125, shk1_126, sii_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_269[k] = pb_z[k] * shk0_125[k]
                   + f_14 * shi_95[k]
                   - f_12 * pc_z[k] * shk1_125[k];

        t_270[k] = pb_z[k] * shk0_126[k]
                   + f_15 * shi_96[k]
                   - f_12 * pc_z[k] * shk1_126[k];

        t_271[k] = f_14 * shi_126[k]
                   + f_3 * pc_y[k] * sii_210[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, pc_x, shi_216, shi_217, shi_218, shi_219, \
                         sih0_167, sih1_167, sii_216, sii_217, sii_218, \
                         sii_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_15 * shi_216[k]
                   + f_10 * sih0_167[k]
                   - f_11 * sih1_167[k]
                   + f_3 * pc_x[k] * sii_216[k];

        t_273[k] = f_15 * shi_217[k]
                   + f_3 * pc_x[k] * sii_217[k];

        t_274[k] = f_15 * shi_218[k]
                   + f_3 * pc_x[k] * sii_218[k];

        t_275[k] = f_15 * shi_219[k]
                   + f_3 * pc_x[k] * sii_219[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pc_x, shi_220, shi_221, shi_222, shi_223, \
                         sii_220, sii_221, sii_222, sii_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_15 * shi_220[k]
                   + f_3 * pc_x[k] * sii_220[k];

        t_277[k] = f_15 * shi_221[k]
                   + f_3 * pc_x[k] * sii_221[k];

        t_278[k] = f_15 * shi_222[k]
                   + f_3 * pc_x[k] * sii_222[k];

        t_279[k] = f_15 * shi_223[k]
                   + f_3 * pc_x[k] * sii_223[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, pb_z, pc_y, pc_z, shk0_136, shi_105, shi_135, \
                         shk1_136, sih0_164, sih1_164, sii_217, \
                         sii_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = pb_z[k] * shk0_136[k]
                   - f_12 * pc_z[k] * shk1_136[k];

        t_281[k] = f_13 * shi_105[k]
                   + f_3 * pc_z[k] * sii_217[k];

        t_282[k] = f_14 * shi_135[k]
                   + f_4 * sih0_164[k]
                   - f_5 * sih1_164[k]
                   + f_3 * pc_y[k] * sii_219[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, pc_y, shi_136, shi_137, shi_138, sih0_165, \
                         sih0_166, sih0_167, sih1_165, sih1_166, sih1_167, sii_220, sii_221, \
                         sii_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_14 * shi_136[k]
                   + f_6 * sih0_165[k]
                   - f_7 * sih1_165[k]
                   + f_3 * pc_y[k] * sii_220[k];

        t_284[k] = f_14 * shi_137[k]
                   + f_8 * sih0_166[k]
                   - f_9 * sih1_166[k]
                   + f_3 * pc_y[k] * sii_221[k];

        t_285[k] = f_14 * shi_138[k]
                   + f_10 * sih0_167[k]
                   - f_11 * sih1_167[k]
                   + f_3 * pc_y[k] * sii_222[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, pb_y, pc_y, pc_z, shk0_180, shi_111, \
                         shi_139, shi_140, shk1_180, sih0_167, sih1_167, sii_223, \
                         sii_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = f_14 * shi_139[k]
                   + f_3 * pc_y[k] * sii_223[k];

        t_287[k] = f_13 * shi_111[k]
                   + f_1 * sih0_167[k]
                   - f_2 * sih1_167[k]
                   + f_3 * pc_z[k] * sii_223[k];

        t_288[k] = pb_y[k] * shk0_180[k]
                   - f_12 * pc_y[k] * shk1_180[k];

        t_289[k] = f_13 * shi_140[k]
                   + f_3 * pc_y[k] * sii_224[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, pb_y, pc_y, pc_z, shk0_183, shk0_185, \
                         shi_112, shi_141, shi_142, shk1_183, shk1_185, sii_224, \
                         sii_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = f_14 * shi_112[k]
                   + f_3 * pc_z[k] * sii_224[k];

        t_291[k] = pb_y[k] * shk0_183[k]
                   + f_14 * shi_141[k]
                   - f_12 * pc_y[k] * shk1_183[k];

        t_292[k] = f_13 * shi_142[k]
                   + f_3 * pc_y[k] * sii_226[k];

        t_293[k] = pb_y[k] * shk0_185[k]
                   - f_12 * pc_y[k] * shk1_185[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, pb_y, pc_y, pc_z, shk0_186, shk0_189, \
                         shi_115, shi_143, shi_145, shk1_186, shk1_189, sii_227, \
                         sii_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = pb_y[k] * shk0_186[k]
                   + f_15 * shi_143[k]
                   - f_12 * pc_y[k] * shk1_186[k];

        t_295[k] = f_14 * shi_115[k]
                   + f_3 * pc_z[k] * sii_227[k];

        t_296[k] = f_13 * shi_145[k]
                   + f_3 * pc_y[k] * sii_229[k];

        t_297[k] = pb_y[k] * shk0_189[k]
                   - f_12 * pc_y[k] * shk1_189[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, pb_y, pc_y, pc_z, shk0_190, shk0_192, shi_118, \
                         shi_146, shi_148, shk1_190, shk1_192, \
                         sii_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = pb_y[k] * shk0_190[k]
                   + f_16 * shi_146[k]
                   - f_12 * pc_y[k] * shk1_190[k];

        t_299[k] = f_14 * shi_118[k]
                   + f_3 * pc_z[k] * sii_230[k];

        t_300[k] = pb_y[k] * shk0_192[k]
                   + f_14 * shi_148[k]
                   - f_12 * pc_y[k] * shk1_192[k];
    }

#pragma omp simd aligned(t_301, t_302, t_303, t_304, pb_y, pc_y, pc_z, shk0_194, shk0_195, \
                         shi_122, shi_149, shi_150, shk1_194, shk1_195, sii_233, \
                         sii_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_301[k] = f_13 * shi_149[k]
                   + f_3 * pc_y[k] * sii_233[k];

        t_302[k] = pb_y[k] * shk0_194[k]
                   - f_12 * pc_y[k] * shk1_194[k];

        t_303[k] = pb_y[k] * shk0_195[k]
                   + f_17 * shi_150[k]
                   - f_12 * pc_y[k] * shk1_195[k];

        t_304[k] = f_14 * shi_122[k]
                   + f_3 * pc_z[k] * sii_234[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, pb_y, pc_y, shk0_197, shk0_198, shk0_200, \
                         shi_152, shi_153, shi_154, shk1_197, shk1_198, shk1_200, \
                         sii_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = pb_y[k] * shk0_197[k]
                   + f_15 * shi_152[k]
                   - f_12 * pc_y[k] * shk1_197[k];

        t_306[k] = pb_y[k] * shk0_198[k]
                   + f_14 * shi_153[k]
                   - f_12 * pc_y[k] * shk1_198[k];

        t_307[k] = f_13 * shi_154[k]
                   + f_3 * pc_y[k] * sii_238[k];

        t_308[k] = pb_y[k] * shk0_200[k]
                   - f_12 * pc_y[k] * shk1_200[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, t_313, pc_x, shi_245, shi_246, shi_247, \
                         shi_248, shi_249, sii_245, sii_246, sii_247, sii_248, \
                         sii_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_15 * shi_245[k]
                   + f_3 * pc_x[k] * sii_245[k];

        t_310[k] = f_15 * shi_246[k]
                   + f_3 * pc_x[k] * sii_246[k];

        t_311[k] = f_15 * shi_247[k]
                   + f_3 * pc_x[k] * sii_247[k];

        t_312[k] = f_15 * shi_248[k]
                   + f_3 * pc_x[k] * sii_248[k];

        t_313[k] = f_15 * shi_249[k]
                   + f_3 * pc_x[k] * sii_249[k];
    }

#pragma omp simd aligned(t_314, t_315, t_316, t_317, pc_x, pc_y, pc_z, shi_133, shi_161, \
                         shi_250, shi_251, sih0_183, sih1_183, sii_245, sii_250, \
                         sii_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_314[k] = f_15 * shi_250[k]
                   + f_3 * pc_x[k] * sii_250[k];

        t_315[k] = f_15 * shi_251[k]
                   + f_3 * pc_x[k] * sii_251[k];

        t_316[k] = f_13 * shi_161[k]
                   + f_1 * sih0_183[k]
                   - f_2 * sih1_183[k]
                   + f_3 * pc_y[k] * sii_245[k];

        t_317[k] = f_14 * shi_133[k]
                   + f_3 * pc_z[k] * sii_245[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, pc_y, shi_163, shi_164, shi_165, sih0_185, \
                         sih0_186, sih0_187, sih1_185, sih1_186, sih1_187, sii_247, sii_248, \
                         sii_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = f_13 * shi_163[k]
                   + f_4 * sih0_185[k]
                   - f_5 * sih1_185[k]
                   + f_3 * pc_y[k] * sii_247[k];

        t_319[k] = f_13 * shi_164[k]
                   + f_6 * sih0_186[k]
                   - f_7 * sih1_186[k]
                   + f_3 * pc_y[k] * sii_248[k];

        t_320[k] = f_13 * shi_165[k]
                   + f_8 * sih0_187[k]
                   - f_9 * sih1_187[k]
                   + f_3 * pc_y[k] * sii_249[k];
    }

#pragma omp simd aligned(t_321, t_322, t_323, pb_y, pc_y, shk0_215, shi_166, shi_167, \
                         shk1_215, sih0_188, sih1_188, sii_250, \
                         sii_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_321[k] = f_13 * shi_166[k]
                   + f_10 * sih0_188[k]
                   - f_11 * sih1_188[k]
                   + f_3 * pc_y[k] * sii_250[k];

        t_322[k] = f_13 * shi_167[k]
                   + f_3 * pc_y[k] * sii_251[k];

        t_323[k] = pb_y[k] * shk0_215[k]
                   - f_12 * pc_y[k] * shk1_215[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, pc_x, pc_y, pc_z, shi_140, shi_252, \
                         shi_255, sih0_189, sih0_192, sih1_189, sih1_192, sii_252, \
                         sii_255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = f_15 * shi_252[k]
                   + f_1 * sih0_189[k]
                   - f_2 * sih1_189[k]
                   + f_3 * pc_x[k] * sii_252[k];

        t_325[k] = f_3 * pc_y[k] * sii_252[k];

        t_326[k] = f_15 * shi_140[k]
                   + f_3 * pc_z[k] * sii_252[k];

        t_327[k] = f_15 * shi_255[k]
                   + f_4 * sih0_192[k]
                   - f_5 * sih1_192[k]
                   + f_3 * pc_x[k] * sii_255[k];
    }

#pragma omp simd aligned(t_328, t_329, t_330, pc_x, pc_y, shi_257, shi_258, sih0_194, \
                         sih0_195, sih1_194, sih1_195, sii_254, sii_257, \
                         sii_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_328[k] = f_3 * pc_y[k] * sii_254[k];

        t_329[k] = f_15 * shi_257[k]
                   + f_4 * sih0_194[k]
                   - f_5 * sih1_194[k]
                   + f_3 * pc_x[k] * sii_257[k];

        t_330[k] = f_15 * shi_258[k]
                   + f_6 * sih0_195[k]
                   - f_7 * sih1_195[k]
                   + f_3 * pc_x[k] * sii_258[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, pc_x, pc_y, pc_z, shi_143, shi_261, sih0_198, \
                         sih1_198, sii_255, sii_257, sii_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = f_15 * shi_143[k]
                   + f_3 * pc_z[k] * sii_255[k];

        t_332[k] = f_3 * pc_y[k] * sii_257[k];

        t_333[k] = f_15 * shi_261[k]
                   + f_6 * sih0_198[k]
                   - f_7 * sih1_198[k]
                   + f_3 * pc_x[k] * sii_261[k];
    }

#pragma omp simd aligned(t_334, t_335, t_336, pc_x, pc_z, shi_146, shi_262, shi_264, sih0_199, \
                         sih0_201, sih1_199, sih1_201, sii_258, sii_262, \
                         sii_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_334[k] = f_15 * shi_262[k]
                   + f_8 * sih0_199[k]
                   - f_9 * sih1_199[k]
                   + f_3 * pc_x[k] * sii_262[k];

        t_335[k] = f_15 * shi_146[k]
                   + f_3 * pc_z[k] * sii_258[k];

        t_336[k] = f_15 * shi_264[k]
                   + f_8 * sih0_201[k]
                   - f_9 * sih1_201[k]
                   + f_3 * pc_x[k] * sii_264[k];
    }

#pragma omp simd aligned(t_337, t_338, t_339, pc_x, pc_y, shi_266, shi_267, sih0_203, \
                         sih0_204, sih1_203, sih1_204, sii_261, sii_266, \
                         sii_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_337[k] = f_3 * pc_y[k] * sii_261[k];

        t_338[k] = f_15 * shi_266[k]
                   + f_8 * sih0_203[k]
                   - f_9 * sih1_203[k]
                   + f_3 * pc_x[k] * sii_266[k];

        t_339[k] = f_15 * shi_267[k]
                   + f_10 * sih0_204[k]
                   - f_11 * sih1_204[k]
                   + f_3 * pc_x[k] * sii_267[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, pc_x, pc_z, shi_150, shi_269, shi_270, sih0_206, \
                         sih0_207, sih1_206, sih1_207, sii_262, sii_269, \
                         sii_270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = f_15 * shi_150[k]
                   + f_3 * pc_z[k] * sii_262[k];

        t_341[k] = f_15 * shi_269[k]
                   + f_10 * sih0_206[k]
                   - f_11 * sih1_206[k]
                   + f_3 * pc_x[k] * sii_269[k];

        t_342[k] = f_15 * shi_270[k]
                   + f_10 * sih0_207[k]
                   - f_11 * sih1_207[k]
                   + f_3 * pc_x[k] * sii_270[k];
    }

#pragma omp simd aligned(t_343, t_344, t_345, t_346, pc_x, pc_y, shi_272, shi_273, shi_274, \
                         sih0_209, sih1_209, sii_266, sii_272, sii_273, \
                         sii_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_343[k] = f_3 * pc_y[k] * sii_266[k];

        t_344[k] = f_15 * shi_272[k]
                   + f_10 * sih0_209[k]
                   - f_11 * sih1_209[k]
                   + f_3 * pc_x[k] * sii_272[k];

        t_345[k] = f_15 * shi_273[k]
                   + f_3 * pc_x[k] * sii_273[k];

        t_346[k] = f_15 * shi_274[k]
                   + f_3 * pc_x[k] * sii_274[k];
    }
}

static auto
compute_prim_sik_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t shk0,
                                                          const size_t shi, const size_t shk1,
                                                          const size_t sih0, const size_t sih1,
                                                          const size_t sii, const size_t ncols,
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

    const auto *shk0_216 = buffer.data(shk0 + 216);
    const auto *shk0_219 = buffer.data(shk0 + 219);
    const auto *shk0_222 = buffer.data(shk0 + 222);
    const auto *shk0_226 = buffer.data(shk0 + 226);
    const auto *shk0_228 = buffer.data(shk0 + 228);
    const auto *shk0_231 = buffer.data(shk0 + 231);
    const auto *shk0_233 = buffer.data(shk0 + 233);
    const auto *shk0_234 = buffer.data(shk0 + 234);
    const auto *shk0_244 = buffer.data(shk0 + 244);

    const auto *shi_161 = buffer.data(shi + 161);
    const auto *shi_167 = buffer.data(shi + 167);
    const auto *shi_168 = buffer.data(shi + 168);
    const auto *shi_170 = buffer.data(shi + 170);
    const auto *shi_171 = buffer.data(shi + 171);
    const auto *shi_173 = buffer.data(shi + 173);
    const auto *shi_174 = buffer.data(shi + 174);
    const auto *shi_175 = buffer.data(shi + 175);
    const auto *shi_177 = buffer.data(shi + 177);
    const auto *shi_178 = buffer.data(shi + 178);
    const auto *shi_179 = buffer.data(shi + 179);
    const auto *shi_180 = buffer.data(shi + 180);
    const auto *shi_182 = buffer.data(shi + 182);
    const auto *shi_189 = buffer.data(shi + 189);
    const auto *shi_191 = buffer.data(shi + 191);
    const auto *shi_192 = buffer.data(shi + 192);
    const auto *shi_193 = buffer.data(shi + 193);
    const auto *shi_194 = buffer.data(shi + 194);
    const auto *shi_195 = buffer.data(shi + 195);
    const auto *shi_196 = buffer.data(shi + 196);
    const auto *shi_198 = buffer.data(shi + 198);
    const auto *shi_199 = buffer.data(shi + 199);
    const auto *shi_201 = buffer.data(shi + 201);
    const auto *shi_202 = buffer.data(shi + 202);
    const auto *shi_205 = buffer.data(shi + 205);
    const auto *shi_206 = buffer.data(shi + 206);
    const auto *shi_210 = buffer.data(shi + 210);
    const auto *shi_219 = buffer.data(shi + 219);
    const auto *shi_220 = buffer.data(shi + 220);
    const auto *shi_221 = buffer.data(shi + 221);
    const auto *shi_222 = buffer.data(shi + 222);
    const auto *shi_223 = buffer.data(shi + 223);
    const auto *shi_224 = buffer.data(shi + 224);
    const auto *shi_226 = buffer.data(shi + 226);
    const auto *shi_229 = buffer.data(shi + 229);
    const auto *shi_233 = buffer.data(shi + 233);
    const auto *shi_238 = buffer.data(shi + 238);
    const auto *shi_275 = buffer.data(shi + 275);
    const auto *shi_276 = buffer.data(shi + 276);
    const auto *shi_277 = buffer.data(shi + 277);
    const auto *shi_278 = buffer.data(shi + 278);
    const auto *shi_279 = buffer.data(shi + 279);
    const auto *shi_280 = buffer.data(shi + 280);
    const auto *shi_283 = buffer.data(shi + 283);
    const auto *shi_285 = buffer.data(shi + 285);
    const auto *shi_286 = buffer.data(shi + 286);
    const auto *shi_289 = buffer.data(shi + 289);
    const auto *shi_290 = buffer.data(shi + 290);
    const auto *shi_292 = buffer.data(shi + 292);
    const auto *shi_294 = buffer.data(shi + 294);
    const auto *shi_295 = buffer.data(shi + 295);
    const auto *shi_297 = buffer.data(shi + 297);
    const auto *shi_298 = buffer.data(shi + 298);
    const auto *shi_300 = buffer.data(shi + 300);
    const auto *shi_301 = buffer.data(shi + 301);
    const auto *shi_302 = buffer.data(shi + 302);
    const auto *shi_303 = buffer.data(shi + 303);
    const auto *shi_304 = buffer.data(shi + 304);
    const auto *shi_305 = buffer.data(shi + 305);
    const auto *shi_306 = buffer.data(shi + 306);
    const auto *shi_307 = buffer.data(shi + 307);
    const auto *shi_313 = buffer.data(shi + 313);
    const auto *shi_317 = buffer.data(shi + 317);
    const auto *shi_322 = buffer.data(shi + 322);
    const auto *shi_328 = buffer.data(shi + 328);
    const auto *shi_329 = buffer.data(shi + 329);
    const auto *shi_330 = buffer.data(shi + 330);
    const auto *shi_331 = buffer.data(shi + 331);
    const auto *shi_332 = buffer.data(shi + 332);
    const auto *shi_333 = buffer.data(shi + 333);
    const auto *shi_334 = buffer.data(shi + 334);
    const auto *shi_335 = buffer.data(shi + 335);
    const auto *shi_336 = buffer.data(shi + 336);
    const auto *shi_339 = buffer.data(shi + 339);
    const auto *shi_341 = buffer.data(shi + 341);
    const auto *shi_342 = buffer.data(shi + 342);
    const auto *shi_345 = buffer.data(shi + 345);
    const auto *shi_346 = buffer.data(shi + 346);
    const auto *shi_348 = buffer.data(shi + 348);
    const auto *shi_350 = buffer.data(shi + 350);
    const auto *shi_351 = buffer.data(shi + 351);
    const auto *shi_353 = buffer.data(shi + 353);
    const auto *shi_354 = buffer.data(shi + 354);
    const auto *shi_356 = buffer.data(shi + 356);
    const auto *shi_357 = buffer.data(shi + 357);
    const auto *shi_358 = buffer.data(shi + 358);
    const auto *shi_359 = buffer.data(shi + 359);

    const auto *shk1_216 = buffer.data(shk1 + 216);
    const auto *shk1_219 = buffer.data(shk1 + 219);
    const auto *shk1_222 = buffer.data(shk1 + 222);
    const auto *shk1_226 = buffer.data(shk1 + 226);
    const auto *shk1_228 = buffer.data(shk1 + 228);
    const auto *shk1_231 = buffer.data(shk1 + 231);
    const auto *shk1_233 = buffer.data(shk1 + 233);
    const auto *shk1_234 = buffer.data(shk1 + 234);
    const auto *shk1_244 = buffer.data(shk1 + 244);

    const auto *sih0_204 = buffer.data(sih0 + 204);
    const auto *sih0_206 = buffer.data(sih0 + 206);
    const auto *sih0_207 = buffer.data(sih0 + 207);
    const auto *sih0_208 = buffer.data(sih0 + 208);
    const auto *sih0_209 = buffer.data(sih0 + 209);
    const auto *sih0_210 = buffer.data(sih0 + 210);
    const auto *sih0_213 = buffer.data(sih0 + 213);
    const auto *sih0_215 = buffer.data(sih0 + 215);
    const auto *sih0_216 = buffer.data(sih0 + 216);
    const auto *sih0_219 = buffer.data(sih0 + 219);
    const auto *sih0_220 = buffer.data(sih0 + 220);
    const auto *sih0_222 = buffer.data(sih0 + 222);
    const auto *sih0_224 = buffer.data(sih0 + 224);
    const auto *sih0_225 = buffer.data(sih0 + 225);
    const auto *sih0_227 = buffer.data(sih0 + 227);
    const auto *sih0_228 = buffer.data(sih0 + 228);
    const auto *sih0_229 = buffer.data(sih0 + 229);
    const auto *sih0_230 = buffer.data(sih0 + 230);
    const auto *sih0_236 = buffer.data(sih0 + 236);
    const auto *sih0_240 = buffer.data(sih0 + 240);
    const auto *sih0_245 = buffer.data(sih0 + 245);
    const auto *sih0_248 = buffer.data(sih0 + 248);
    const auto *sih0_249 = buffer.data(sih0 + 249);
    const auto *sih0_250 = buffer.data(sih0 + 250);
    const auto *sih0_251 = buffer.data(sih0 + 251);
    const auto *sih0_252 = buffer.data(sih0 + 252);
    const auto *sih0_255 = buffer.data(sih0 + 255);
    const auto *sih0_257 = buffer.data(sih0 + 257);
    const auto *sih0_258 = buffer.data(sih0 + 258);
    const auto *sih0_261 = buffer.data(sih0 + 261);
    const auto *sih0_262 = buffer.data(sih0 + 262);
    const auto *sih0_264 = buffer.data(sih0 + 264);
    const auto *sih0_266 = buffer.data(sih0 + 266);
    const auto *sih0_267 = buffer.data(sih0 + 267);
    const auto *sih0_269 = buffer.data(sih0 + 269);
    const auto *sih0_270 = buffer.data(sih0 + 270);
    const auto *sih0_272 = buffer.data(sih0 + 272);

    const auto *sih1_204 = buffer.data(sih1 + 204);
    const auto *sih1_206 = buffer.data(sih1 + 206);
    const auto *sih1_207 = buffer.data(sih1 + 207);
    const auto *sih1_208 = buffer.data(sih1 + 208);
    const auto *sih1_209 = buffer.data(sih1 + 209);
    const auto *sih1_210 = buffer.data(sih1 + 210);
    const auto *sih1_213 = buffer.data(sih1 + 213);
    const auto *sih1_215 = buffer.data(sih1 + 215);
    const auto *sih1_216 = buffer.data(sih1 + 216);
    const auto *sih1_219 = buffer.data(sih1 + 219);
    const auto *sih1_220 = buffer.data(sih1 + 220);
    const auto *sih1_222 = buffer.data(sih1 + 222);
    const auto *sih1_224 = buffer.data(sih1 + 224);
    const auto *sih1_225 = buffer.data(sih1 + 225);
    const auto *sih1_227 = buffer.data(sih1 + 227);
    const auto *sih1_228 = buffer.data(sih1 + 228);
    const auto *sih1_229 = buffer.data(sih1 + 229);
    const auto *sih1_230 = buffer.data(sih1 + 230);
    const auto *sih1_236 = buffer.data(sih1 + 236);
    const auto *sih1_240 = buffer.data(sih1 + 240);
    const auto *sih1_245 = buffer.data(sih1 + 245);
    const auto *sih1_248 = buffer.data(sih1 + 248);
    const auto *sih1_249 = buffer.data(sih1 + 249);
    const auto *sih1_250 = buffer.data(sih1 + 250);
    const auto *sih1_251 = buffer.data(sih1 + 251);
    const auto *sih1_252 = buffer.data(sih1 + 252);
    const auto *sih1_255 = buffer.data(sih1 + 255);
    const auto *sih1_257 = buffer.data(sih1 + 257);
    const auto *sih1_258 = buffer.data(sih1 + 258);
    const auto *sih1_261 = buffer.data(sih1 + 261);
    const auto *sih1_262 = buffer.data(sih1 + 262);
    const auto *sih1_264 = buffer.data(sih1 + 264);
    const auto *sih1_266 = buffer.data(sih1 + 266);
    const auto *sih1_267 = buffer.data(sih1 + 267);
    const auto *sih1_269 = buffer.data(sih1 + 269);
    const auto *sih1_270 = buffer.data(sih1 + 270);
    const auto *sih1_272 = buffer.data(sih1 + 272);

    const auto *sii_273 = buffer.data(sii + 273);
    const auto *sii_275 = buffer.data(sii + 275);
    const auto *sii_276 = buffer.data(sii + 276);
    const auto *sii_277 = buffer.data(sii + 277);
    const auto *sii_278 = buffer.data(sii + 278);
    const auto *sii_279 = buffer.data(sii + 279);
    const auto *sii_280 = buffer.data(sii + 280);
    const auto *sii_282 = buffer.data(sii + 282);
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
    const auto *sii_308 = buffer.data(sii + 308);
    const auto *sii_310 = buffer.data(sii + 310);
    const auto *sii_311 = buffer.data(sii + 311);
    const auto *sii_313 = buffer.data(sii + 313);
    const auto *sii_314 = buffer.data(sii + 314);
    const auto *sii_317 = buffer.data(sii + 317);
    const auto *sii_318 = buffer.data(sii + 318);
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
    const auto *sii_338 = buffer.data(sii + 338);
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

#pragma omp simd aligned(t_347, t_348, t_349, t_350, t_351, pc_x, shi_275, shi_276, shi_277, \
                         shi_278, shi_279, sii_275, sii_276, sii_277, sii_278, \
                         sii_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_347[k] = f_15 * shi_275[k]
                   + f_3 * pc_x[k] * sii_275[k];

        t_348[k] = f_15 * shi_276[k]
                   + f_3 * pc_x[k] * sii_276[k];

        t_349[k] = f_15 * shi_277[k]
                   + f_3 * pc_x[k] * sii_277[k];

        t_350[k] = f_15 * shi_278[k]
                   + f_3 * pc_x[k] * sii_278[k];

        t_351[k] = f_15 * shi_279[k]
                   + f_3 * pc_x[k] * sii_279[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, t_355, pc_y, pc_z, shi_161, sih0_204, sih0_206, \
                         sih0_207, sih1_204, sih1_206, sih1_207, sii_273, sii_275, \
                         sii_276 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = f_1 * sih0_204[k]
                   - f_2 * sih1_204[k]
                   + f_3 * pc_y[k] * sii_273[k];

        t_353[k] = f_15 * shi_161[k]
                   + f_3 * pc_z[k] * sii_273[k];

        t_354[k] = f_4 * sih0_206[k]
                   - f_5 * sih1_206[k]
                   + f_3 * pc_y[k] * sii_275[k];

        t_355[k] = f_6 * sih0_207[k]
                   - f_7 * sih1_207[k]
                   + f_3 * pc_y[k] * sii_276[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, t_359, pc_y, pc_z, shi_167, sih0_208, sih0_209, \
                         sih1_208, sih1_209, sii_277, sii_278, \
                         sii_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_8 * sih0_208[k]
                   - f_9 * sih1_208[k]
                   + f_3 * pc_y[k] * sii_277[k];

        t_357[k] = f_10 * sih0_209[k]
                   - f_11 * sih1_209[k]
                   + f_3 * pc_y[k] * sii_278[k];

        t_358[k] = f_3 * pc_y[k] * sii_279[k];

        t_359[k] = f_15 * shi_167[k]
                   + f_1 * sih0_209[k]
                   - f_2 * sih1_209[k]
                   + f_3 * pc_z[k] * sii_279[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, pc_x, pc_y, pc_z, shi_168, shi_280, \
                         shi_283, sih0_210, sih0_213, sih1_210, sih1_213, sii_280, \
                         sii_283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = f_14 * shi_280[k]
                   + f_1 * sih0_210[k]
                   - f_2 * sih1_210[k]
                   + f_3 * pc_x[k] * sii_280[k];

        t_361[k] = f_16 * shi_168[k]
                   + f_3 * pc_y[k] * sii_280[k];

        t_362[k] = f_3 * pc_z[k] * sii_280[k];

        t_363[k] = f_14 * shi_283[k]
                   + f_4 * sih0_213[k]
                   - f_5 * sih1_213[k]
                   + f_3 * pc_x[k] * sii_283[k];
    }

#pragma omp simd aligned(t_364, t_365, t_366, pc_x, pc_y, shi_170, shi_285, shi_286, sih0_215, \
                         sih0_216, sih1_215, sih1_216, sii_282, sii_285, \
                         sii_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_364[k] = f_16 * shi_170[k]
                   + f_3 * pc_y[k] * sii_282[k];

        t_365[k] = f_14 * shi_285[k]
                   + f_4 * sih0_215[k]
                   - f_5 * sih1_215[k]
                   + f_3 * pc_x[k] * sii_285[k];

        t_366[k] = f_14 * shi_286[k]
                   + f_6 * sih0_216[k]
                   - f_7 * sih1_216[k]
                   + f_3 * pc_x[k] * sii_286[k];
    }

#pragma omp simd aligned(t_367, t_368, t_369, pc_x, pc_y, pc_z, shi_173, shi_289, sih0_219, \
                         sih1_219, sii_283, sii_285, sii_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_367[k] = f_3 * pc_z[k] * sii_283[k];

        t_368[k] = f_16 * shi_173[k]
                   + f_3 * pc_y[k] * sii_285[k];

        t_369[k] = f_14 * shi_289[k]
                   + f_6 * sih0_219[k]
                   - f_7 * sih1_219[k]
                   + f_3 * pc_x[k] * sii_289[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, pc_x, pc_z, shi_290, shi_292, sih0_220, \
                         sih0_222, sih1_220, sih1_222, sii_286, sii_290, \
                         sii_292 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = f_14 * shi_290[k]
                   + f_8 * sih0_220[k]
                   - f_9 * sih1_220[k]
                   + f_3 * pc_x[k] * sii_290[k];

        t_371[k] = f_3 * pc_z[k] * sii_286[k];

        t_372[k] = f_14 * shi_292[k]
                   + f_8 * sih0_222[k]
                   - f_9 * sih1_222[k]
                   + f_3 * pc_x[k] * sii_292[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, pc_x, pc_y, shi_177, shi_294, shi_295, sih0_224, \
                         sih0_225, sih1_224, sih1_225, sii_289, sii_294, \
                         sii_295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = f_16 * shi_177[k]
                   + f_3 * pc_y[k] * sii_289[k];

        t_374[k] = f_14 * shi_294[k]
                   + f_8 * sih0_224[k]
                   - f_9 * sih1_224[k]
                   + f_3 * pc_x[k] * sii_294[k];

        t_375[k] = f_14 * shi_295[k]
                   + f_10 * sih0_225[k]
                   - f_11 * sih1_225[k]
                   + f_3 * pc_x[k] * sii_295[k];
    }

#pragma omp simd aligned(t_376, t_377, t_378, pc_x, pc_z, shi_297, shi_298, sih0_227, \
                         sih0_228, sih1_227, sih1_228, sii_290, sii_297, \
                         sii_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_376[k] = f_3 * pc_z[k] * sii_290[k];

        t_377[k] = f_14 * shi_297[k]
                   + f_10 * sih0_227[k]
                   - f_11 * sih1_227[k]
                   + f_3 * pc_x[k] * sii_297[k];

        t_378[k] = f_14 * shi_298[k]
                   + f_10 * sih0_228[k]
                   - f_11 * sih1_228[k]
                   + f_3 * pc_x[k] * sii_298[k];
    }

#pragma omp simd aligned(t_379, t_380, t_381, t_382, pc_x, pc_y, shi_182, shi_300, shi_301, \
                         shi_302, sih0_230, sih1_230, sii_294, sii_300, sii_301, \
                         sii_302 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_379[k] = f_16 * shi_182[k]
                   + f_3 * pc_y[k] * sii_294[k];

        t_380[k] = f_14 * shi_300[k]
                   + f_10 * sih0_230[k]
                   - f_11 * sih1_230[k]
                   + f_3 * pc_x[k] * sii_300[k];

        t_381[k] = f_14 * shi_301[k]
                   + f_3 * pc_x[k] * sii_301[k];

        t_382[k] = f_14 * shi_302[k]
                   + f_3 * pc_x[k] * sii_302[k];
    }

#pragma omp simd aligned(t_383, t_384, t_385, t_386, t_387, pc_x, shi_303, shi_304, shi_305, \
                         shi_306, shi_307, sii_303, sii_304, sii_305, sii_306, \
                         sii_307 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_383[k] = f_14 * shi_303[k]
                   + f_3 * pc_x[k] * sii_303[k];

        t_384[k] = f_14 * shi_304[k]
                   + f_3 * pc_x[k] * sii_304[k];

        t_385[k] = f_14 * shi_305[k]
                   + f_3 * pc_x[k] * sii_305[k];

        t_386[k] = f_14 * shi_306[k]
                   + f_3 * pc_x[k] * sii_306[k];

        t_387[k] = f_14 * shi_307[k]
                   + f_3 * pc_x[k] * sii_307[k];
    }

#pragma omp simd aligned(t_388, t_389, t_390, pc_y, pc_z, shi_189, shi_191, sih0_225, \
                         sih0_227, sih1_225, sih1_227, sii_301, \
                         sii_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_388[k] = f_16 * shi_189[k]
                   + f_1 * sih0_225[k]
                   - f_2 * sih1_225[k]
                   + f_3 * pc_y[k] * sii_301[k];

        t_389[k] = f_3 * pc_z[k] * sii_301[k];

        t_390[k] = f_16 * shi_191[k]
                   + f_4 * sih0_227[k]
                   - f_5 * sih1_227[k]
                   + f_3 * pc_y[k] * sii_303[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, pc_y, shi_192, shi_193, shi_194, sih0_228, \
                         sih0_229, sih0_230, sih1_228, sih1_229, sih1_230, sii_304, sii_305, \
                         sii_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = f_16 * shi_192[k]
                   + f_6 * sih0_228[k]
                   - f_7 * sih1_228[k]
                   + f_3 * pc_y[k] * sii_304[k];

        t_392[k] = f_16 * shi_193[k]
                   + f_8 * sih0_229[k]
                   - f_9 * sih1_229[k]
                   + f_3 * pc_y[k] * sii_305[k];

        t_393[k] = f_16 * shi_194[k]
                   + f_10 * sih0_230[k]
                   - f_11 * sih1_230[k]
                   + f_3 * pc_y[k] * sii_306[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, t_397, pb_z, pc_y, pc_z, shk0_216, shi_195, \
                         shi_196, shk1_216, sih0_230, sih1_230, sii_307, \
                         sii_308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = f_16 * shi_195[k]
                   + f_3 * pc_y[k] * sii_307[k];

        t_395[k] = f_1 * sih0_230[k]
                   - f_2 * sih1_230[k]
                   + f_3 * pc_z[k] * sii_307[k];

        t_396[k] = pb_z[k] * shk0_216[k]
                   - f_12 * pc_z[k] * shk1_216[k];

        t_397[k] = f_15 * shi_196[k]
                   + f_3 * pc_y[k] * sii_308[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, pb_z, pc_y, pc_z, shk0_219, shi_168, shi_198, \
                         shk1_219, sii_308, sii_310 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = f_13 * shi_168[k]
                   + f_3 * pc_z[k] * sii_308[k];

        t_399[k] = pb_z[k] * shk0_219[k]
                   - f_12 * pc_z[k] * shk1_219[k];

        t_400[k] = f_15 * shi_198[k]
                   + f_3 * pc_y[k] * sii_310[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, pb_z, pc_x, pc_z, shk0_222, shi_171, shi_313, \
                         shk1_222, sih0_236, sih1_236, sii_311, \
                         sii_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_14 * shi_313[k]
                   + f_4 * sih0_236[k]
                   - f_5 * sih1_236[k]
                   + f_3 * pc_x[k] * sii_313[k];

        t_402[k] = pb_z[k] * shk0_222[k]
                   - f_12 * pc_z[k] * shk1_222[k];

        t_403[k] = f_13 * shi_171[k]
                   + f_3 * pc_z[k] * sii_311[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, pb_z, pc_x, pc_y, pc_z, shk0_226, shi_201, \
                         shi_317, shk1_226, sih0_240, sih1_240, sii_313, \
                         sii_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = f_15 * shi_201[k]
                   + f_3 * pc_y[k] * sii_313[k];

        t_405[k] = f_14 * shi_317[k]
                   + f_6 * sih0_240[k]
                   - f_7 * sih1_240[k]
                   + f_3 * pc_x[k] * sii_317[k];

        t_406[k] = pb_z[k] * shk0_226[k]
                   - f_12 * pc_z[k] * shk1_226[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, pb_z, pc_y, pc_z, shk0_228, shi_174, shi_175, \
                         shi_205, shk1_228, sii_314, sii_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = f_13 * shi_174[k]
                   + f_3 * pc_z[k] * sii_314[k];

        t_408[k] = pb_z[k] * shk0_228[k]
                   + f_14 * shi_175[k]
                   - f_12 * pc_z[k] * shk1_228[k];

        t_409[k] = f_15 * shi_205[k]
                   + f_3 * pc_y[k] * sii_317[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, pb_z, pc_x, pc_z, shk0_231, shi_178, shi_322, \
                         shk1_231, sih0_245, sih1_245, sii_318, \
                         sii_322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = f_14 * shi_322[k]
                   + f_8 * sih0_245[k]
                   - f_9 * sih1_245[k]
                   + f_3 * pc_x[k] * sii_322[k];

        t_411[k] = pb_z[k] * shk0_231[k]
                   - f_12 * pc_z[k] * shk1_231[k];

        t_412[k] = f_13 * shi_178[k]
                   + f_3 * pc_z[k] * sii_318[k];
    }

#pragma omp simd aligned(t_413, t_414, t_415, pb_z, pc_y, pc_z, shk0_233, shk0_234, shi_179, \
                         shi_180, shi_210, shk1_233, shk1_234, \
                         sii_322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_413[k] = pb_z[k] * shk0_233[k]
                   + f_14 * shi_179[k]
                   - f_12 * pc_z[k] * shk1_233[k];

        t_414[k] = pb_z[k] * shk0_234[k]
                   + f_15 * shi_180[k]
                   - f_12 * pc_z[k] * shk1_234[k];

        t_415[k] = f_15 * shi_210[k]
                   + f_3 * pc_y[k] * sii_322[k];
    }

#pragma omp simd aligned(t_416, t_417, t_418, t_419, pc_x, shi_328, shi_329, shi_330, shi_331, \
                         sih0_251, sih1_251, sii_328, sii_329, sii_330, \
                         sii_331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_416[k] = f_14 * shi_328[k]
                   + f_10 * sih0_251[k]
                   - f_11 * sih1_251[k]
                   + f_3 * pc_x[k] * sii_328[k];

        t_417[k] = f_14 * shi_329[k]
                   + f_3 * pc_x[k] * sii_329[k];

        t_418[k] = f_14 * shi_330[k]
                   + f_3 * pc_x[k] * sii_330[k];

        t_419[k] = f_14 * shi_331[k]
                   + f_3 * pc_x[k] * sii_331[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, pc_x, shi_332, shi_333, shi_334, shi_335, \
                         sii_332, sii_333, sii_334, sii_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_14 * shi_332[k]
                   + f_3 * pc_x[k] * sii_332[k];

        t_421[k] = f_14 * shi_333[k]
                   + f_3 * pc_x[k] * sii_333[k];

        t_422[k] = f_14 * shi_334[k]
                   + f_3 * pc_x[k] * sii_334[k];

        t_423[k] = f_14 * shi_335[k]
                   + f_3 * pc_x[k] * sii_335[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, pb_z, pc_y, pc_z, shk0_244, shi_189, shi_219, \
                         shk1_244, sih0_248, sih1_248, sii_329, \
                         sii_331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = pb_z[k] * shk0_244[k]
                   - f_12 * pc_z[k] * shk1_244[k];

        t_425[k] = f_13 * shi_189[k]
                   + f_3 * pc_z[k] * sii_329[k];

        t_426[k] = f_15 * shi_219[k]
                   + f_4 * sih0_248[k]
                   - f_5 * sih1_248[k]
                   + f_3 * pc_y[k] * sii_331[k];
    }

#pragma omp simd aligned(t_427, t_428, t_429, pc_y, shi_220, shi_221, shi_222, sih0_249, \
                         sih0_250, sih0_251, sih1_249, sih1_250, sih1_251, sii_332, sii_333, \
                         sii_334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_427[k] = f_15 * shi_220[k]
                   + f_6 * sih0_249[k]
                   - f_7 * sih1_249[k]
                   + f_3 * pc_y[k] * sii_332[k];

        t_428[k] = f_15 * shi_221[k]
                   + f_8 * sih0_250[k]
                   - f_9 * sih1_250[k]
                   + f_3 * pc_y[k] * sii_333[k];

        t_429[k] = f_15 * shi_222[k]
                   + f_10 * sih0_251[k]
                   - f_11 * sih1_251[k]
                   + f_3 * pc_y[k] * sii_334[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, pc_x, pc_y, pc_z, shi_195, shi_223, shi_336, \
                         sih0_251, sih0_252, sih1_251, sih1_252, sii_335, \
                         sii_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = f_15 * shi_223[k]
                   + f_3 * pc_y[k] * sii_335[k];

        t_431[k] = f_13 * shi_195[k]
                   + f_1 * sih0_251[k]
                   - f_2 * sih1_251[k]
                   + f_3 * pc_z[k] * sii_335[k];

        t_432[k] = f_14 * shi_336[k]
                   + f_1 * sih0_252[k]
                   - f_2 * sih1_252[k]
                   + f_3 * pc_x[k] * sii_336[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, t_436, pc_x, pc_y, pc_z, shi_196, shi_224, \
                         shi_226, shi_339, sih0_255, sih1_255, sii_336, sii_338, \
                         sii_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = f_14 * shi_224[k]
                   + f_3 * pc_y[k] * sii_336[k];

        t_434[k] = f_14 * shi_196[k]
                   + f_3 * pc_z[k] * sii_336[k];

        t_435[k] = f_14 * shi_339[k]
                   + f_4 * sih0_255[k]
                   - f_5 * sih1_255[k]
                   + f_3 * pc_x[k] * sii_339[k];

        t_436[k] = f_14 * shi_226[k]
                   + f_3 * pc_y[k] * sii_338[k];
    }

#pragma omp simd aligned(t_437, t_438, t_439, pc_x, pc_z, shi_199, shi_341, shi_342, sih0_257, \
                         sih0_258, sih1_257, sih1_258, sii_339, sii_341, \
                         sii_342 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_437[k] = f_14 * shi_341[k]
                   + f_4 * sih0_257[k]
                   - f_5 * sih1_257[k]
                   + f_3 * pc_x[k] * sii_341[k];

        t_438[k] = f_14 * shi_342[k]
                   + f_6 * sih0_258[k]
                   - f_7 * sih1_258[k]
                   + f_3 * pc_x[k] * sii_342[k];

        t_439[k] = f_14 * shi_199[k]
                   + f_3 * pc_z[k] * sii_339[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, pc_x, pc_y, shi_229, shi_345, shi_346, sih0_261, \
                         sih0_262, sih1_261, sih1_262, sii_341, sii_345, \
                         sii_346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = f_14 * shi_229[k]
                   + f_3 * pc_y[k] * sii_341[k];

        t_441[k] = f_14 * shi_345[k]
                   + f_6 * sih0_261[k]
                   - f_7 * sih1_261[k]
                   + f_3 * pc_x[k] * sii_345[k];

        t_442[k] = f_14 * shi_346[k]
                   + f_8 * sih0_262[k]
                   - f_9 * sih1_262[k]
                   + f_3 * pc_x[k] * sii_346[k];
    }

#pragma omp simd aligned(t_443, t_444, t_445, pc_x, pc_y, pc_z, shi_202, shi_233, shi_348, \
                         sih0_264, sih1_264, sii_342, sii_345, \
                         sii_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_443[k] = f_14 * shi_202[k]
                   + f_3 * pc_z[k] * sii_342[k];

        t_444[k] = f_14 * shi_348[k]
                   + f_8 * sih0_264[k]
                   - f_9 * sih1_264[k]
                   + f_3 * pc_x[k] * sii_348[k];

        t_445[k] = f_14 * shi_233[k]
                   + f_3 * pc_y[k] * sii_345[k];
    }

#pragma omp simd aligned(t_446, t_447, t_448, pc_x, pc_z, shi_206, shi_350, shi_351, sih0_266, \
                         sih0_267, sih1_266, sih1_267, sii_346, sii_350, \
                         sii_351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_446[k] = f_14 * shi_350[k]
                   + f_8 * sih0_266[k]
                   - f_9 * sih1_266[k]
                   + f_3 * pc_x[k] * sii_350[k];

        t_447[k] = f_14 * shi_351[k]
                   + f_10 * sih0_267[k]
                   - f_11 * sih1_267[k]
                   + f_3 * pc_x[k] * sii_351[k];

        t_448[k] = f_14 * shi_206[k]
                   + f_3 * pc_z[k] * sii_346[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, pc_x, pc_y, shi_238, shi_353, shi_354, sih0_269, \
                         sih0_270, sih1_269, sih1_270, sii_350, sii_353, \
                         sii_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_14 * shi_353[k]
                   + f_10 * sih0_269[k]
                   - f_11 * sih1_269[k]
                   + f_3 * pc_x[k] * sii_353[k];

        t_450[k] = f_14 * shi_354[k]
                   + f_10 * sih0_270[k]
                   - f_11 * sih1_270[k]
                   + f_3 * pc_x[k] * sii_354[k];

        t_451[k] = f_14 * shi_238[k]
                   + f_3 * pc_y[k] * sii_350[k];
    }

#pragma omp simd aligned(t_452, t_453, t_454, t_455, pc_x, shi_356, shi_357, shi_358, shi_359, \
                         sih0_272, sih1_272, sii_356, sii_357, sii_358, \
                         sii_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_452[k] = f_14 * shi_356[k]
                   + f_10 * sih0_272[k]
                   - f_11 * sih1_272[k]
                   + f_3 * pc_x[k] * sii_356[k];

        t_453[k] = f_14 * shi_357[k]
                   + f_3 * pc_x[k] * sii_357[k];

        t_454[k] = f_14 * shi_358[k]
                   + f_3 * pc_x[k] * sii_358[k];

        t_455[k] = f_14 * shi_359[k]
                   + f_3 * pc_x[k] * sii_359[k];
    }
}

static auto
compute_prim_sik_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t shk0,
                                                          const size_t shi, const size_t shk1,
                                                          const size_t sih0, const size_t sih1,
                                                          const size_t sii, const size_t ncols,
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
    const auto f_18 = 3.5 / q;

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
    auto *t_571 = buffer.data(target + 571);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *shk0_324 = buffer.data(shk0 + 324);
    const auto *shk0_327 = buffer.data(shk0 + 327);
    const auto *shk0_329 = buffer.data(shk0 + 329);
    const auto *shk0_330 = buffer.data(shk0 + 330);
    const auto *shk0_333 = buffer.data(shk0 + 333);
    const auto *shk0_334 = buffer.data(shk0 + 334);
    const auto *shk0_336 = buffer.data(shk0 + 336);
    const auto *shk0_338 = buffer.data(shk0 + 338);
    const auto *shk0_339 = buffer.data(shk0 + 339);
    const auto *shk0_341 = buffer.data(shk0 + 341);
    const auto *shk0_342 = buffer.data(shk0 + 342);
    const auto *shk0_344 = buffer.data(shk0 + 344);
    const auto *shk0_359 = buffer.data(shk0 + 359);
    const auto *shk0_540 = buffer.data(shk0 + 540);
    const auto *shk0_543 = buffer.data(shk0 + 543);
    const auto *shk0_545 = buffer.data(shk0 + 545);
    const auto *shk0_546 = buffer.data(shk0 + 546);
    const auto *shk0_549 = buffer.data(shk0 + 549);
    const auto *shk0_550 = buffer.data(shk0 + 550);
    const auto *shk0_552 = buffer.data(shk0 + 552);
    const auto *shk0_554 = buffer.data(shk0 + 554);
    const auto *shk0_555 = buffer.data(shk0 + 555);
    const auto *shk0_557 = buffer.data(shk0 + 557);
    const auto *shk0_558 = buffer.data(shk0 + 558);
    const auto *shk0_560 = buffer.data(shk0 + 560);
    const auto *shk0_568 = buffer.data(shk0 + 568);
    const auto *shk0_570 = buffer.data(shk0 + 570);
    const auto *shk0_571 = buffer.data(shk0 + 571);

    const auto *shi_217 = buffer.data(shi + 217);
    const auto *shi_223 = buffer.data(shi + 223);
    const auto *shi_224 = buffer.data(shi + 224);
    const auto *shi_227 = buffer.data(shi + 227);
    const auto *shi_230 = buffer.data(shi + 230);
    const auto *shi_234 = buffer.data(shi + 234);
    const auto *shi_245 = buffer.data(shi + 245);
    const auto *shi_247 = buffer.data(shi + 247);
    const auto *shi_248 = buffer.data(shi + 248);
    const auto *shi_249 = buffer.data(shi + 249);
    const auto *shi_250 = buffer.data(shi + 250);
    const auto *shi_251 = buffer.data(shi + 251);
    const auto *shi_252 = buffer.data(shi + 252);
    const auto *shi_253 = buffer.data(shi + 253);
    const auto *shi_254 = buffer.data(shi + 254);
    const auto *shi_255 = buffer.data(shi + 255);
    const auto *shi_257 = buffer.data(shi + 257);
    const auto *shi_258 = buffer.data(shi + 258);
    const auto *shi_260 = buffer.data(shi + 260);
    const auto *shi_261 = buffer.data(shi + 261);
    const auto *shi_262 = buffer.data(shi + 262);
    const auto *shi_264 = buffer.data(shi + 264);
    const auto *shi_265 = buffer.data(shi + 265);
    const auto *shi_266 = buffer.data(shi + 266);
    const auto *shi_273 = buffer.data(shi + 273);
    const auto *shi_275 = buffer.data(shi + 275);
    const auto *shi_276 = buffer.data(shi + 276);
    const auto *shi_277 = buffer.data(shi + 277);
    const auto *shi_278 = buffer.data(shi + 278);
    const auto *shi_279 = buffer.data(shi + 279);
    const auto *shi_280 = buffer.data(shi + 280);
    const auto *shi_282 = buffer.data(shi + 282);
    const auto *shi_285 = buffer.data(shi + 285);
    const auto *shi_289 = buffer.data(shi + 289);
    const auto *shi_294 = buffer.data(shi + 294);
    const auto *shi_360 = buffer.data(shi + 360);
    const auto *shi_361 = buffer.data(shi + 361);
    const auto *shi_362 = buffer.data(shi + 362);
    const auto *shi_363 = buffer.data(shi + 363);
    const auto *shi_385 = buffer.data(shi + 385);
    const auto *shi_386 = buffer.data(shi + 386);
    const auto *shi_387 = buffer.data(shi + 387);
    const auto *shi_388 = buffer.data(shi + 388);
    const auto *shi_389 = buffer.data(shi + 389);
    const auto *shi_390 = buffer.data(shi + 390);
    const auto *shi_391 = buffer.data(shi + 391);
    const auto *shi_392 = buffer.data(shi + 392);
    const auto *shi_395 = buffer.data(shi + 395);
    const auto *shi_397 = buffer.data(shi + 397);
    const auto *shi_398 = buffer.data(shi + 398);
    const auto *shi_401 = buffer.data(shi + 401);
    const auto *shi_402 = buffer.data(shi + 402);
    const auto *shi_404 = buffer.data(shi + 404);
    const auto *shi_406 = buffer.data(shi + 406);
    const auto *shi_407 = buffer.data(shi + 407);
    const auto *shi_409 = buffer.data(shi + 409);
    const auto *shi_410 = buffer.data(shi + 410);
    const auto *shi_412 = buffer.data(shi + 412);
    const auto *shi_413 = buffer.data(shi + 413);
    const auto *shi_414 = buffer.data(shi + 414);
    const auto *shi_415 = buffer.data(shi + 415);
    const auto *shi_416 = buffer.data(shi + 416);
    const auto *shi_417 = buffer.data(shi + 417);
    const auto *shi_418 = buffer.data(shi + 418);
    const auto *shi_419 = buffer.data(shi + 419);
    const auto *shi_420 = buffer.data(shi + 420);
    const auto *shi_423 = buffer.data(shi + 423);
    const auto *shi_425 = buffer.data(shi + 425);
    const auto *shi_426 = buffer.data(shi + 426);
    const auto *shi_429 = buffer.data(shi + 429);
    const auto *shi_430 = buffer.data(shi + 430);
    const auto *shi_432 = buffer.data(shi + 432);
    const auto *shi_434 = buffer.data(shi + 434);
    const auto *shi_435 = buffer.data(shi + 435);
    const auto *shi_437 = buffer.data(shi + 437);
    const auto *shi_438 = buffer.data(shi + 438);
    const auto *shi_440 = buffer.data(shi + 440);
    const auto *shi_441 = buffer.data(shi + 441);
    const auto *shi_442 = buffer.data(shi + 442);
    const auto *shi_443 = buffer.data(shi + 443);
    const auto *shi_444 = buffer.data(shi + 444);
    const auto *shi_445 = buffer.data(shi + 445);
    const auto *shi_446 = buffer.data(shi + 446);
    const auto *shi_447 = buffer.data(shi + 447);

    const auto *shk1_324 = buffer.data(shk1 + 324);
    const auto *shk1_327 = buffer.data(shk1 + 327);
    const auto *shk1_329 = buffer.data(shk1 + 329);
    const auto *shk1_330 = buffer.data(shk1 + 330);
    const auto *shk1_333 = buffer.data(shk1 + 333);
    const auto *shk1_334 = buffer.data(shk1 + 334);
    const auto *shk1_336 = buffer.data(shk1 + 336);
    const auto *shk1_338 = buffer.data(shk1 + 338);
    const auto *shk1_339 = buffer.data(shk1 + 339);
    const auto *shk1_341 = buffer.data(shk1 + 341);
    const auto *shk1_342 = buffer.data(shk1 + 342);
    const auto *shk1_344 = buffer.data(shk1 + 344);
    const auto *shk1_359 = buffer.data(shk1 + 359);
    const auto *shk1_540 = buffer.data(shk1 + 540);
    const auto *shk1_543 = buffer.data(shk1 + 543);
    const auto *shk1_545 = buffer.data(shk1 + 545);
    const auto *shk1_546 = buffer.data(shk1 + 546);
    const auto *shk1_549 = buffer.data(shk1 + 549);
    const auto *shk1_550 = buffer.data(shk1 + 550);
    const auto *shk1_552 = buffer.data(shk1 + 552);
    const auto *shk1_554 = buffer.data(shk1 + 554);
    const auto *shk1_555 = buffer.data(shk1 + 555);
    const auto *shk1_557 = buffer.data(shk1 + 557);
    const auto *shk1_558 = buffer.data(shk1 + 558);
    const auto *shk1_560 = buffer.data(shk1 + 560);
    const auto *shk1_568 = buffer.data(shk1 + 568);
    const auto *shk1_570 = buffer.data(shk1 + 570);
    const auto *shk1_571 = buffer.data(shk1 + 571);

    const auto *sih0_267 = buffer.data(sih0 + 267);
    const auto *sih0_269 = buffer.data(sih0 + 269);
    const auto *sih0_270 = buffer.data(sih0 + 270);
    const auto *sih0_271 = buffer.data(sih0 + 271);
    const auto *sih0_272 = buffer.data(sih0 + 272);
    const auto *sih0_288 = buffer.data(sih0 + 288);
    const auto *sih0_290 = buffer.data(sih0 + 290);
    const auto *sih0_291 = buffer.data(sih0 + 291);
    const auto *sih0_292 = buffer.data(sih0 + 292);
    const auto *sih0_293 = buffer.data(sih0 + 293);
    const auto *sih0_294 = buffer.data(sih0 + 294);
    const auto *sih0_297 = buffer.data(sih0 + 297);
    const auto *sih0_299 = buffer.data(sih0 + 299);
    const auto *sih0_300 = buffer.data(sih0 + 300);
    const auto *sih0_303 = buffer.data(sih0 + 303);
    const auto *sih0_304 = buffer.data(sih0 + 304);
    const auto *sih0_306 = buffer.data(sih0 + 306);
    const auto *sih0_308 = buffer.data(sih0 + 308);
    const auto *sih0_309 = buffer.data(sih0 + 309);
    const auto *sih0_311 = buffer.data(sih0 + 311);
    const auto *sih0_312 = buffer.data(sih0 + 312);
    const auto *sih0_313 = buffer.data(sih0 + 313);
    const auto *sih0_314 = buffer.data(sih0 + 314);

    const auto *sih1_267 = buffer.data(sih1 + 267);
    const auto *sih1_269 = buffer.data(sih1 + 269);
    const auto *sih1_270 = buffer.data(sih1 + 270);
    const auto *sih1_271 = buffer.data(sih1 + 271);
    const auto *sih1_272 = buffer.data(sih1 + 272);
    const auto *sih1_288 = buffer.data(sih1 + 288);
    const auto *sih1_290 = buffer.data(sih1 + 290);
    const auto *sih1_291 = buffer.data(sih1 + 291);
    const auto *sih1_292 = buffer.data(sih1 + 292);
    const auto *sih1_293 = buffer.data(sih1 + 293);
    const auto *sih1_294 = buffer.data(sih1 + 294);
    const auto *sih1_297 = buffer.data(sih1 + 297);
    const auto *sih1_299 = buffer.data(sih1 + 299);
    const auto *sih1_300 = buffer.data(sih1 + 300);
    const auto *sih1_303 = buffer.data(sih1 + 303);
    const auto *sih1_304 = buffer.data(sih1 + 304);
    const auto *sih1_306 = buffer.data(sih1 + 306);
    const auto *sih1_308 = buffer.data(sih1 + 308);
    const auto *sih1_309 = buffer.data(sih1 + 309);
    const auto *sih1_311 = buffer.data(sih1 + 311);
    const auto *sih1_312 = buffer.data(sih1 + 312);
    const auto *sih1_313 = buffer.data(sih1 + 313);
    const auto *sih1_314 = buffer.data(sih1 + 314);

    const auto *sii_357 = buffer.data(sii + 357);
    const auto *sii_359 = buffer.data(sii + 359);
    const auto *sii_360 = buffer.data(sii + 360);
    const auto *sii_361 = buffer.data(sii + 361);
    const auto *sii_362 = buffer.data(sii + 362);
    const auto *sii_363 = buffer.data(sii + 363);
    const auto *sii_364 = buffer.data(sii + 364);
    const auto *sii_366 = buffer.data(sii + 366);
    const auto *sii_367 = buffer.data(sii + 367);
    const auto *sii_369 = buffer.data(sii + 369);
    const auto *sii_370 = buffer.data(sii + 370);
    const auto *sii_373 = buffer.data(sii + 373);
    const auto *sii_374 = buffer.data(sii + 374);
    const auto *sii_378 = buffer.data(sii + 378);
    const auto *sii_385 = buffer.data(sii + 385);
    const auto *sii_386 = buffer.data(sii + 386);
    const auto *sii_387 = buffer.data(sii + 387);
    const auto *sii_388 = buffer.data(sii + 388);
    const auto *sii_389 = buffer.data(sii + 389);
    const auto *sii_390 = buffer.data(sii + 390);
    const auto *sii_391 = buffer.data(sii + 391);
    const auto *sii_392 = buffer.data(sii + 392);
    const auto *sii_394 = buffer.data(sii + 394);
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
    const auto *sii_422 = buffer.data(sii + 422);
    const auto *sii_423 = buffer.data(sii + 423);
    const auto *sii_425 = buffer.data(sii + 425);
    const auto *sii_426 = buffer.data(sii + 426);
    const auto *sii_429 = buffer.data(sii + 429);
    const auto *sii_430 = buffer.data(sii + 430);
    const auto *sii_434 = buffer.data(sii + 434);
    const auto *sii_441 = buffer.data(sii + 441);
    const auto *sii_442 = buffer.data(sii + 442);
    const auto *sii_443 = buffer.data(sii + 443);
    const auto *sii_444 = buffer.data(sii + 444);
    const auto *sii_445 = buffer.data(sii + 445);
    const auto *sii_446 = buffer.data(sii + 446);
    const auto *sii_447 = buffer.data(sii + 447);

#pragma omp simd aligned(t_456, t_457, t_458, t_459, pc_x, shi_360, shi_361, shi_362, shi_363, \
                         sii_360, sii_361, sii_362, sii_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_456[k] = f_14 * shi_360[k]
                   + f_3 * pc_x[k] * sii_360[k];

        t_457[k] = f_14 * shi_361[k]
                   + f_3 * pc_x[k] * sii_361[k];

        t_458[k] = f_14 * shi_362[k]
                   + f_3 * pc_x[k] * sii_362[k];

        t_459[k] = f_14 * shi_363[k]
                   + f_3 * pc_x[k] * sii_363[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, pc_y, pc_z, shi_217, shi_245, shi_247, sih0_267, \
                         sih0_269, sih1_267, sih1_269, sii_357, \
                         sii_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = f_14 * shi_245[k]
                   + f_1 * sih0_267[k]
                   - f_2 * sih1_267[k]
                   + f_3 * pc_y[k] * sii_357[k];

        t_461[k] = f_14 * shi_217[k]
                   + f_3 * pc_z[k] * sii_357[k];

        t_462[k] = f_14 * shi_247[k]
                   + f_4 * sih0_269[k]
                   - f_5 * sih1_269[k]
                   + f_3 * pc_y[k] * sii_359[k];
    }

#pragma omp simd aligned(t_463, t_464, t_465, pc_y, shi_248, shi_249, shi_250, sih0_270, \
                         sih0_271, sih0_272, sih1_270, sih1_271, sih1_272, sii_360, sii_361, \
                         sii_362 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_463[k] = f_14 * shi_248[k]
                   + f_6 * sih0_270[k]
                   - f_7 * sih1_270[k]
                   + f_3 * pc_y[k] * sii_360[k];

        t_464[k] = f_14 * shi_249[k]
                   + f_8 * sih0_271[k]
                   - f_9 * sih1_271[k]
                   + f_3 * pc_y[k] * sii_361[k];

        t_465[k] = f_14 * shi_250[k]
                   + f_10 * sih0_272[k]
                   - f_11 * sih1_272[k]
                   + f_3 * pc_y[k] * sii_362[k];
    }

#pragma omp simd aligned(t_466, t_467, t_468, t_469, pb_y, pc_y, pc_z, shk0_324, shi_223, \
                         shi_251, shi_252, shk1_324, sih0_272, sih1_272, sii_363, \
                         sii_364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_466[k] = f_14 * shi_251[k]
                   + f_3 * pc_y[k] * sii_363[k];

        t_467[k] = f_14 * shi_223[k]
                   + f_1 * sih0_272[k]
                   - f_2 * sih1_272[k]
                   + f_3 * pc_z[k] * sii_363[k];

        t_468[k] = pb_y[k] * shk0_324[k]
                   - f_12 * pc_y[k] * shk1_324[k];

        t_469[k] = f_13 * shi_252[k]
                   + f_3 * pc_y[k] * sii_364[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, pb_y, pc_y, pc_z, shk0_327, shk0_329, \
                         shi_224, shi_253, shi_254, shk1_327, shk1_329, sii_364, \
                         sii_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = f_15 * shi_224[k]
                   + f_3 * pc_z[k] * sii_364[k];

        t_471[k] = pb_y[k] * shk0_327[k]
                   + f_14 * shi_253[k]
                   - f_12 * pc_y[k] * shk1_327[k];

        t_472[k] = f_13 * shi_254[k]
                   + f_3 * pc_y[k] * sii_366[k];

        t_473[k] = pb_y[k] * shk0_329[k]
                   - f_12 * pc_y[k] * shk1_329[k];
    }

#pragma omp simd aligned(t_474, t_475, t_476, t_477, pb_y, pc_y, pc_z, shk0_330, shk0_333, \
                         shi_227, shi_255, shi_257, shk1_330, shk1_333, sii_367, \
                         sii_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = pb_y[k] * shk0_330[k]
                   + f_15 * shi_255[k]
                   - f_12 * pc_y[k] * shk1_330[k];

        t_475[k] = f_15 * shi_227[k]
                   + f_3 * pc_z[k] * sii_367[k];

        t_476[k] = f_13 * shi_257[k]
                   + f_3 * pc_y[k] * sii_369[k];

        t_477[k] = pb_y[k] * shk0_333[k]
                   - f_12 * pc_y[k] * shk1_333[k];
    }

#pragma omp simd aligned(t_478, t_479, t_480, pb_y, pc_y, pc_z, shk0_334, shk0_336, shi_230, \
                         shi_258, shi_260, shk1_334, shk1_336, \
                         sii_370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_478[k] = pb_y[k] * shk0_334[k]
                   + f_16 * shi_258[k]
                   - f_12 * pc_y[k] * shk1_334[k];

        t_479[k] = f_15 * shi_230[k]
                   + f_3 * pc_z[k] * sii_370[k];

        t_480[k] = pb_y[k] * shk0_336[k]
                   + f_14 * shi_260[k]
                   - f_12 * pc_y[k] * shk1_336[k];
    }

#pragma omp simd aligned(t_481, t_482, t_483, t_484, pb_y, pc_y, pc_z, shk0_338, shk0_339, \
                         shi_234, shi_261, shi_262, shk1_338, shk1_339, sii_373, \
                         sii_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_481[k] = f_13 * shi_261[k]
                   + f_3 * pc_y[k] * sii_373[k];

        t_482[k] = pb_y[k] * shk0_338[k]
                   - f_12 * pc_y[k] * shk1_338[k];

        t_483[k] = pb_y[k] * shk0_339[k]
                   + f_17 * shi_262[k]
                   - f_12 * pc_y[k] * shk1_339[k];

        t_484[k] = f_15 * shi_234[k]
                   + f_3 * pc_z[k] * sii_374[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, t_488, pb_y, pc_y, shk0_341, shk0_342, shk0_344, \
                         shi_264, shi_265, shi_266, shk1_341, shk1_342, shk1_344, \
                         sii_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = pb_y[k] * shk0_341[k]
                   + f_15 * shi_264[k]
                   - f_12 * pc_y[k] * shk1_341[k];

        t_486[k] = pb_y[k] * shk0_342[k]
                   + f_14 * shi_265[k]
                   - f_12 * pc_y[k] * shk1_342[k];

        t_487[k] = f_13 * shi_266[k]
                   + f_3 * pc_y[k] * sii_378[k];

        t_488[k] = pb_y[k] * shk0_344[k]
                   - f_12 * pc_y[k] * shk1_344[k];
    }

#pragma omp simd aligned(t_489, t_490, t_491, t_492, t_493, pc_x, shi_385, shi_386, shi_387, \
                         shi_388, shi_389, sii_385, sii_386, sii_387, sii_388, \
                         sii_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_489[k] = f_14 * shi_385[k]
                   + f_3 * pc_x[k] * sii_385[k];

        t_490[k] = f_14 * shi_386[k]
                   + f_3 * pc_x[k] * sii_386[k];

        t_491[k] = f_14 * shi_387[k]
                   + f_3 * pc_x[k] * sii_387[k];

        t_492[k] = f_14 * shi_388[k]
                   + f_3 * pc_x[k] * sii_388[k];

        t_493[k] = f_14 * shi_389[k]
                   + f_3 * pc_x[k] * sii_389[k];
    }

#pragma omp simd aligned(t_494, t_495, t_496, t_497, pc_x, pc_y, pc_z, shi_245, shi_273, \
                         shi_390, shi_391, sih0_288, sih1_288, sii_385, sii_390, \
                         sii_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_494[k] = f_14 * shi_390[k]
                   + f_3 * pc_x[k] * sii_390[k];

        t_495[k] = f_14 * shi_391[k]
                   + f_3 * pc_x[k] * sii_391[k];

        t_496[k] = f_13 * shi_273[k]
                   + f_1 * sih0_288[k]
                   - f_2 * sih1_288[k]
                   + f_3 * pc_y[k] * sii_385[k];

        t_497[k] = f_15 * shi_245[k]
                   + f_3 * pc_z[k] * sii_385[k];
    }

#pragma omp simd aligned(t_498, t_499, t_500, pc_y, shi_275, shi_276, shi_277, sih0_290, \
                         sih0_291, sih0_292, sih1_290, sih1_291, sih1_292, sii_387, sii_388, \
                         sii_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_498[k] = f_13 * shi_275[k]
                   + f_4 * sih0_290[k]
                   - f_5 * sih1_290[k]
                   + f_3 * pc_y[k] * sii_387[k];

        t_499[k] = f_13 * shi_276[k]
                   + f_6 * sih0_291[k]
                   - f_7 * sih1_291[k]
                   + f_3 * pc_y[k] * sii_388[k];

        t_500[k] = f_13 * shi_277[k]
                   + f_8 * sih0_292[k]
                   - f_9 * sih1_292[k]
                   + f_3 * pc_y[k] * sii_389[k];
    }

#pragma omp simd aligned(t_501, t_502, t_503, pb_y, pc_y, shk0_359, shi_278, shi_279, \
                         shk1_359, sih0_293, sih1_293, sii_390, \
                         sii_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_501[k] = f_13 * shi_278[k]
                   + f_10 * sih0_293[k]
                   - f_11 * sih1_293[k]
                   + f_3 * pc_y[k] * sii_390[k];

        t_502[k] = f_13 * shi_279[k]
                   + f_3 * pc_y[k] * sii_391[k];

        t_503[k] = pb_y[k] * shk0_359[k]
                   - f_12 * pc_y[k] * shk1_359[k];
    }

#pragma omp simd aligned(t_504, t_505, t_506, t_507, pc_x, pc_y, pc_z, shi_252, shi_392, \
                         shi_395, sih0_294, sih0_297, sih1_294, sih1_297, sii_392, \
                         sii_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_504[k] = f_14 * shi_392[k]
                   + f_1 * sih0_294[k]
                   - f_2 * sih1_294[k]
                   + f_3 * pc_x[k] * sii_392[k];

        t_505[k] = f_3 * pc_y[k] * sii_392[k];

        t_506[k] = f_16 * shi_252[k]
                   + f_3 * pc_z[k] * sii_392[k];

        t_507[k] = f_14 * shi_395[k]
                   + f_4 * sih0_297[k]
                   - f_5 * sih1_297[k]
                   + f_3 * pc_x[k] * sii_395[k];
    }

#pragma omp simd aligned(t_508, t_509, t_510, pc_x, pc_y, shi_397, shi_398, sih0_299, \
                         sih0_300, sih1_299, sih1_300, sii_394, sii_397, \
                         sii_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_508[k] = f_3 * pc_y[k] * sii_394[k];

        t_509[k] = f_14 * shi_397[k]
                   + f_4 * sih0_299[k]
                   - f_5 * sih1_299[k]
                   + f_3 * pc_x[k] * sii_397[k];

        t_510[k] = f_14 * shi_398[k]
                   + f_6 * sih0_300[k]
                   - f_7 * sih1_300[k]
                   + f_3 * pc_x[k] * sii_398[k];
    }

#pragma omp simd aligned(t_511, t_512, t_513, pc_x, pc_y, pc_z, shi_255, shi_401, sih0_303, \
                         sih1_303, sii_395, sii_397, sii_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_511[k] = f_16 * shi_255[k]
                   + f_3 * pc_z[k] * sii_395[k];

        t_512[k] = f_3 * pc_y[k] * sii_397[k];

        t_513[k] = f_14 * shi_401[k]
                   + f_6 * sih0_303[k]
                   - f_7 * sih1_303[k]
                   + f_3 * pc_x[k] * sii_401[k];
    }

#pragma omp simd aligned(t_514, t_515, t_516, pc_x, pc_z, shi_258, shi_402, shi_404, sih0_304, \
                         sih0_306, sih1_304, sih1_306, sii_398, sii_402, \
                         sii_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_514[k] = f_14 * shi_402[k]
                   + f_8 * sih0_304[k]
                   - f_9 * sih1_304[k]
                   + f_3 * pc_x[k] * sii_402[k];

        t_515[k] = f_16 * shi_258[k]
                   + f_3 * pc_z[k] * sii_398[k];

        t_516[k] = f_14 * shi_404[k]
                   + f_8 * sih0_306[k]
                   - f_9 * sih1_306[k]
                   + f_3 * pc_x[k] * sii_404[k];
    }

#pragma omp simd aligned(t_517, t_518, t_519, pc_x, pc_y, shi_406, shi_407, sih0_308, \
                         sih0_309, sih1_308, sih1_309, sii_401, sii_406, \
                         sii_407 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_517[k] = f_3 * pc_y[k] * sii_401[k];

        t_518[k] = f_14 * shi_406[k]
                   + f_8 * sih0_308[k]
                   - f_9 * sih1_308[k]
                   + f_3 * pc_x[k] * sii_406[k];

        t_519[k] = f_14 * shi_407[k]
                   + f_10 * sih0_309[k]
                   - f_11 * sih1_309[k]
                   + f_3 * pc_x[k] * sii_407[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, pc_x, pc_z, shi_262, shi_409, shi_410, sih0_311, \
                         sih0_312, sih1_311, sih1_312, sii_402, sii_409, \
                         sii_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = f_16 * shi_262[k]
                   + f_3 * pc_z[k] * sii_402[k];

        t_521[k] = f_14 * shi_409[k]
                   + f_10 * sih0_311[k]
                   - f_11 * sih1_311[k]
                   + f_3 * pc_x[k] * sii_409[k];

        t_522[k] = f_14 * shi_410[k]
                   + f_10 * sih0_312[k]
                   - f_11 * sih1_312[k]
                   + f_3 * pc_x[k] * sii_410[k];
    }

#pragma omp simd aligned(t_523, t_524, t_525, t_526, pc_x, pc_y, shi_412, shi_413, shi_414, \
                         sih0_314, sih1_314, sii_406, sii_412, sii_413, \
                         sii_414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_523[k] = f_3 * pc_y[k] * sii_406[k];

        t_524[k] = f_14 * shi_412[k]
                   + f_10 * sih0_314[k]
                   - f_11 * sih1_314[k]
                   + f_3 * pc_x[k] * sii_412[k];

        t_525[k] = f_14 * shi_413[k]
                   + f_3 * pc_x[k] * sii_413[k];

        t_526[k] = f_14 * shi_414[k]
                   + f_3 * pc_x[k] * sii_414[k];
    }

#pragma omp simd aligned(t_527, t_528, t_529, t_530, t_531, pc_x, shi_415, shi_416, shi_417, \
                         shi_418, shi_419, sii_415, sii_416, sii_417, sii_418, \
                         sii_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_527[k] = f_14 * shi_415[k]
                   + f_3 * pc_x[k] * sii_415[k];

        t_528[k] = f_14 * shi_416[k]
                   + f_3 * pc_x[k] * sii_416[k];

        t_529[k] = f_14 * shi_417[k]
                   + f_3 * pc_x[k] * sii_417[k];

        t_530[k] = f_14 * shi_418[k]
                   + f_3 * pc_x[k] * sii_418[k];

        t_531[k] = f_14 * shi_419[k]
                   + f_3 * pc_x[k] * sii_419[k];
    }

#pragma omp simd aligned(t_532, t_533, t_534, t_535, pc_y, pc_z, shi_273, sih0_309, sih0_311, \
                         sih0_312, sih1_309, sih1_311, sih1_312, sii_413, sii_415, \
                         sii_416 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_532[k] = f_1 * sih0_309[k]
                   - f_2 * sih1_309[k]
                   + f_3 * pc_y[k] * sii_413[k];

        t_533[k] = f_16 * shi_273[k]
                   + f_3 * pc_z[k] * sii_413[k];

        t_534[k] = f_4 * sih0_311[k]
                   - f_5 * sih1_311[k]
                   + f_3 * pc_y[k] * sii_415[k];

        t_535[k] = f_6 * sih0_312[k]
                   - f_7 * sih1_312[k]
                   + f_3 * pc_y[k] * sii_416[k];
    }

#pragma omp simd aligned(t_536, t_537, t_538, t_539, pc_y, pc_z, shi_279, sih0_313, sih0_314, \
                         sih1_313, sih1_314, sii_417, sii_418, \
                         sii_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_536[k] = f_8 * sih0_313[k]
                   - f_9 * sih1_313[k]
                   + f_3 * pc_y[k] * sii_417[k];

        t_537[k] = f_10 * sih0_314[k]
                   - f_11 * sih1_314[k]
                   + f_3 * pc_y[k] * sii_418[k];

        t_538[k] = f_3 * pc_y[k] * sii_419[k];

        t_539[k] = f_16 * shi_279[k]
                   + f_1 * sih0_314[k]
                   - f_2 * sih1_314[k]
                   + f_3 * pc_z[k] * sii_419[k];
    }

#pragma omp simd aligned(t_540, t_541, t_542, t_543, pb_x, pc_x, pc_y, pc_z, shk0_540, \
                         shk0_543, shi_280, shi_420, shi_423, shk1_540, shk1_543, \
                         sii_420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_540[k] = pb_x[k] * shk0_540[k]
                   + f_18 * shi_420[k]
                   - f_12 * pc_x[k] * shk1_540[k];

        t_541[k] = f_17 * shi_280[k]
                   + f_3 * pc_y[k] * sii_420[k];

        t_542[k] = f_3 * pc_z[k] * sii_420[k];

        t_543[k] = pb_x[k] * shk0_543[k]
                   + f_17 * shi_423[k]
                   - f_12 * pc_x[k] * shk1_543[k];
    }

#pragma omp simd aligned(t_544, t_545, t_546, pb_x, pc_x, pc_y, shk0_545, shk0_546, shi_282, \
                         shi_425, shi_426, shk1_545, shk1_546, \
                         sii_422 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_544[k] = f_17 * shi_282[k]
                   + f_3 * pc_y[k] * sii_422[k];

        t_545[k] = pb_x[k] * shk0_545[k]
                   + f_17 * shi_425[k]
                   - f_12 * pc_x[k] * shk1_545[k];

        t_546[k] = pb_x[k] * shk0_546[k]
                   + f_16 * shi_426[k]
                   - f_12 * pc_x[k] * shk1_546[k];
    }

#pragma omp simd aligned(t_547, t_548, t_549, pb_x, pc_x, pc_y, pc_z, shk0_549, shi_285, \
                         shi_429, shk1_549, sii_423, sii_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_547[k] = f_3 * pc_z[k] * sii_423[k];

        t_548[k] = f_17 * shi_285[k]
                   + f_3 * pc_y[k] * sii_425[k];

        t_549[k] = pb_x[k] * shk0_549[k]
                   + f_16 * shi_429[k]
                   - f_12 * pc_x[k] * shk1_549[k];
    }

#pragma omp simd aligned(t_550, t_551, t_552, pb_x, pc_x, pc_z, shk0_550, shk0_552, shi_430, \
                         shi_432, shk1_550, shk1_552, sii_426 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_550[k] = pb_x[k] * shk0_550[k]
                   + f_15 * shi_430[k]
                   - f_12 * pc_x[k] * shk1_550[k];

        t_551[k] = f_3 * pc_z[k] * sii_426[k];

        t_552[k] = pb_x[k] * shk0_552[k]
                   + f_15 * shi_432[k]
                   - f_12 * pc_x[k] * shk1_552[k];
    }

#pragma omp simd aligned(t_553, t_554, t_555, pb_x, pc_x, pc_y, shk0_554, shk0_555, shi_289, \
                         shi_434, shi_435, shk1_554, shk1_555, \
                         sii_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_553[k] = f_17 * shi_289[k]
                   + f_3 * pc_y[k] * sii_429[k];

        t_554[k] = pb_x[k] * shk0_554[k]
                   + f_15 * shi_434[k]
                   - f_12 * pc_x[k] * shk1_554[k];

        t_555[k] = pb_x[k] * shk0_555[k]
                   + f_14 * shi_435[k]
                   - f_12 * pc_x[k] * shk1_555[k];
    }

#pragma omp simd aligned(t_556, t_557, t_558, pb_x, pc_x, pc_z, shk0_557, shk0_558, shi_437, \
                         shi_438, shk1_557, shk1_558, sii_430 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_556[k] = f_3 * pc_z[k] * sii_430[k];

        t_557[k] = pb_x[k] * shk0_557[k]
                   + f_14 * shi_437[k]
                   - f_12 * pc_x[k] * shk1_557[k];

        t_558[k] = pb_x[k] * shk0_558[k]
                   + f_14 * shi_438[k]
                   - f_12 * pc_x[k] * shk1_558[k];
    }

#pragma omp simd aligned(t_559, t_560, t_561, t_562, pb_x, pc_x, pc_y, shk0_560, shi_294, \
                         shi_440, shi_441, shi_442, shk1_560, sii_434, sii_441, \
                         sii_442 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_559[k] = f_17 * shi_294[k]
                   + f_3 * pc_y[k] * sii_434[k];

        t_560[k] = pb_x[k] * shk0_560[k]
                   + f_14 * shi_440[k]
                   - f_12 * pc_x[k] * shk1_560[k];

        t_561[k] = f_13 * shi_441[k]
                   + f_3 * pc_x[k] * sii_441[k];

        t_562[k] = f_13 * shi_442[k]
                   + f_3 * pc_x[k] * sii_442[k];
    }

#pragma omp simd aligned(t_563, t_564, t_565, t_566, t_567, pc_x, shi_443, shi_444, shi_445, \
                         shi_446, shi_447, sii_443, sii_444, sii_445, sii_446, \
                         sii_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_563[k] = f_13 * shi_443[k]
                   + f_3 * pc_x[k] * sii_443[k];

        t_564[k] = f_13 * shi_444[k]
                   + f_3 * pc_x[k] * sii_444[k];

        t_565[k] = f_13 * shi_445[k]
                   + f_3 * pc_x[k] * sii_445[k];

        t_566[k] = f_13 * shi_446[k]
                   + f_3 * pc_x[k] * sii_446[k];

        t_567[k] = f_13 * shi_447[k]
                   + f_3 * pc_x[k] * sii_447[k];
    }

#pragma omp simd aligned(t_568, t_569, t_570, t_571, pb_x, pc_x, pc_z, shk0_568, shk0_570, \
                         shk0_571, shk1_568, shk1_570, shk1_571, \
                         sii_441 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_568[k] = pb_x[k] * shk0_568[k]
                   - f_12 * pc_x[k] * shk1_568[k];

        t_569[k] = f_3 * pc_z[k] * sii_441[k];

        t_570[k] = pb_x[k] * shk0_570[k]
                   - f_12 * pc_x[k] * shk1_570[k];

        t_571[k] = pb_x[k] * shk0_571[k]
                   - f_12 * pc_x[k] * shk1_571[k];
    }
}

static auto
compute_prim_sik_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t shk0,
                                                          const size_t shi, const size_t shk1,
                                                          const size_t sii, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_3 = p / q;
    const auto f_12 = gamma / q;
    const auto f_13 = 0.5 / q;
    const auto f_14 = 1.0 / q;
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;
    const auto f_18 = 3.5 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *shk0_360 = buffer.data(shk0 + 360);
    const auto *shk0_363 = buffer.data(shk0 + 363);
    const auto *shk0_366 = buffer.data(shk0 + 366);
    const auto *shk0_370 = buffer.data(shk0 + 370);
    const auto *shk0_375 = buffer.data(shk0 + 375);
    const auto *shk0_504 = buffer.data(shk0 + 504);
    const auto *shk0_509 = buffer.data(shk0 + 509);
    const auto *shk0_572 = buffer.data(shk0 + 572);
    const auto *shk0_573 = buffer.data(shk0 + 573);
    const auto *shk0_575 = buffer.data(shk0 + 575);
    const auto *shk0_581 = buffer.data(shk0 + 581);
    const auto *shk0_585 = buffer.data(shk0 + 585);
    const auto *shk0_588 = buffer.data(shk0 + 588);
    const auto *shk0_590 = buffer.data(shk0 + 590);
    const auto *shk0_593 = buffer.data(shk0 + 593);
    const auto *shk0_594 = buffer.data(shk0 + 594);
    const auto *shk0_596 = buffer.data(shk0 + 596);
    const auto *shk0_604 = buffer.data(shk0 + 604);
    const auto *shk0_606 = buffer.data(shk0 + 606);
    const auto *shk0_607 = buffer.data(shk0 + 607);
    const auto *shk0_608 = buffer.data(shk0 + 608);
    const auto *shk0_609 = buffer.data(shk0 + 609);
    const auto *shk0_611 = buffer.data(shk0 + 611);
    const auto *shk0_612 = buffer.data(shk0 + 612);
    const auto *shk0_615 = buffer.data(shk0 + 615);
    const auto *shk0_617 = buffer.data(shk0 + 617);
    const auto *shk0_618 = buffer.data(shk0 + 618);
    const auto *shk0_621 = buffer.data(shk0 + 621);
    const auto *shk0_622 = buffer.data(shk0 + 622);
    const auto *shk0_624 = buffer.data(shk0 + 624);
    const auto *shk0_626 = buffer.data(shk0 + 626);
    const auto *shk0_627 = buffer.data(shk0 + 627);
    const auto *shk0_629 = buffer.data(shk0 + 629);
    const auto *shk0_630 = buffer.data(shk0 + 630);
    const auto *shk0_632 = buffer.data(shk0 + 632);
    const auto *shk0_640 = buffer.data(shk0 + 640);
    const auto *shk0_642 = buffer.data(shk0 + 642);
    const auto *shk0_643 = buffer.data(shk0 + 643);
    const auto *shk0_644 = buffer.data(shk0 + 644);
    const auto *shk0_645 = buffer.data(shk0 + 645);
    const auto *shk0_647 = buffer.data(shk0 + 647);
    const auto *shk0_648 = buffer.data(shk0 + 648);
    const auto *shk0_651 = buffer.data(shk0 + 651);
    const auto *shk0_653 = buffer.data(shk0 + 653);
    const auto *shk0_654 = buffer.data(shk0 + 654);
    const auto *shk0_657 = buffer.data(shk0 + 657);
    const auto *shk0_658 = buffer.data(shk0 + 658);
    const auto *shk0_660 = buffer.data(shk0 + 660);
    const auto *shk0_662 = buffer.data(shk0 + 662);
    const auto *shk0_663 = buffer.data(shk0 + 663);
    const auto *shk0_665 = buffer.data(shk0 + 665);
    const auto *shk0_666 = buffer.data(shk0 + 666);
    const auto *shk0_668 = buffer.data(shk0 + 668);
    const auto *shk0_676 = buffer.data(shk0 + 676);
    const auto *shk0_678 = buffer.data(shk0 + 678);
    const auto *shk0_679 = buffer.data(shk0 + 679);
    const auto *shk0_680 = buffer.data(shk0 + 680);
    const auto *shk0_681 = buffer.data(shk0 + 681);
    const auto *shk0_683 = buffer.data(shk0 + 683);
    const auto *shk0_687 = buffer.data(shk0 + 687);

    const auto *shi_280 = buffer.data(shi + 280);
    const auto *shi_283 = buffer.data(shi + 283);
    const auto *shi_286 = buffer.data(shi + 286);
    const auto *shi_290 = buffer.data(shi + 290);
    const auto *shi_301 = buffer.data(shi + 301);
    const auto *shi_307 = buffer.data(shi + 307);
    const auto *shi_308 = buffer.data(shi + 308);
    const auto *shi_310 = buffer.data(shi + 310);
    const auto *shi_311 = buffer.data(shi + 311);
    const auto *shi_313 = buffer.data(shi + 313);
    const auto *shi_314 = buffer.data(shi + 314);
    const auto *shi_317 = buffer.data(shi + 317);
    const auto *shi_318 = buffer.data(shi + 318);
    const auto *shi_322 = buffer.data(shi + 322);
    const auto *shi_329 = buffer.data(shi + 329);
    const auto *shi_335 = buffer.data(shi + 335);
    const auto *shi_336 = buffer.data(shi + 336);
    const auto *shi_338 = buffer.data(shi + 338);
    const auto *shi_339 = buffer.data(shi + 339);
    const auto *shi_341 = buffer.data(shi + 341);
    const auto *shi_342 = buffer.data(shi + 342);
    const auto *shi_345 = buffer.data(shi + 345);
    const auto *shi_346 = buffer.data(shi + 346);
    const auto *shi_350 = buffer.data(shi + 350);
    const auto *shi_357 = buffer.data(shi + 357);
    const auto *shi_363 = buffer.data(shi + 363);
    const auto *shi_364 = buffer.data(shi + 364);
    const auto *shi_366 = buffer.data(shi + 366);
    const auto *shi_369 = buffer.data(shi + 369);
    const auto *shi_373 = buffer.data(shi + 373);
    const auto *shi_378 = buffer.data(shi + 378);
    const auto *shi_391 = buffer.data(shi + 391);
    const auto *shi_392 = buffer.data(shi + 392);
    const auto *shi_394 = buffer.data(shi + 394);
    const auto *shi_453 = buffer.data(shi + 453);
    const auto *shi_457 = buffer.data(shi + 457);
    const auto *shi_460 = buffer.data(shi + 460);
    const auto *shi_462 = buffer.data(shi + 462);
    const auto *shi_465 = buffer.data(shi + 465);
    const auto *shi_466 = buffer.data(shi + 466);
    const auto *shi_468 = buffer.data(shi + 468);
    const auto *shi_469 = buffer.data(shi + 469);
    const auto *shi_470 = buffer.data(shi + 470);
    const auto *shi_471 = buffer.data(shi + 471);
    const auto *shi_472 = buffer.data(shi + 472);
    const auto *shi_473 = buffer.data(shi + 473);
    const auto *shi_474 = buffer.data(shi + 474);
    const auto *shi_475 = buffer.data(shi + 475);
    const auto *shi_476 = buffer.data(shi + 476);
    const auto *shi_479 = buffer.data(shi + 479);
    const auto *shi_481 = buffer.data(shi + 481);
    const auto *shi_482 = buffer.data(shi + 482);
    const auto *shi_485 = buffer.data(shi + 485);
    const auto *shi_486 = buffer.data(shi + 486);
    const auto *shi_488 = buffer.data(shi + 488);
    const auto *shi_490 = buffer.data(shi + 490);
    const auto *shi_491 = buffer.data(shi + 491);
    const auto *shi_493 = buffer.data(shi + 493);
    const auto *shi_494 = buffer.data(shi + 494);
    const auto *shi_496 = buffer.data(shi + 496);
    const auto *shi_497 = buffer.data(shi + 497);
    const auto *shi_498 = buffer.data(shi + 498);
    const auto *shi_499 = buffer.data(shi + 499);
    const auto *shi_500 = buffer.data(shi + 500);
    const auto *shi_501 = buffer.data(shi + 501);
    const auto *shi_502 = buffer.data(shi + 502);
    const auto *shi_503 = buffer.data(shi + 503);
    const auto *shi_504 = buffer.data(shi + 504);
    const auto *shi_507 = buffer.data(shi + 507);
    const auto *shi_509 = buffer.data(shi + 509);
    const auto *shi_510 = buffer.data(shi + 510);
    const auto *shi_513 = buffer.data(shi + 513);
    const auto *shi_514 = buffer.data(shi + 514);
    const auto *shi_516 = buffer.data(shi + 516);
    const auto *shi_518 = buffer.data(shi + 518);
    const auto *shi_519 = buffer.data(shi + 519);
    const auto *shi_521 = buffer.data(shi + 521);
    const auto *shi_522 = buffer.data(shi + 522);
    const auto *shi_524 = buffer.data(shi + 524);
    const auto *shi_525 = buffer.data(shi + 525);
    const auto *shi_526 = buffer.data(shi + 526);
    const auto *shi_527 = buffer.data(shi + 527);
    const auto *shi_528 = buffer.data(shi + 528);
    const auto *shi_529 = buffer.data(shi + 529);
    const auto *shi_530 = buffer.data(shi + 530);
    const auto *shi_531 = buffer.data(shi + 531);
    const auto *shi_535 = buffer.data(shi + 535);

    const auto *shk1_360 = buffer.data(shk1 + 360);
    const auto *shk1_363 = buffer.data(shk1 + 363);
    const auto *shk1_366 = buffer.data(shk1 + 366);
    const auto *shk1_370 = buffer.data(shk1 + 370);
    const auto *shk1_375 = buffer.data(shk1 + 375);
    const auto *shk1_504 = buffer.data(shk1 + 504);
    const auto *shk1_509 = buffer.data(shk1 + 509);
    const auto *shk1_572 = buffer.data(shk1 + 572);
    const auto *shk1_573 = buffer.data(shk1 + 573);
    const auto *shk1_575 = buffer.data(shk1 + 575);
    const auto *shk1_581 = buffer.data(shk1 + 581);
    const auto *shk1_585 = buffer.data(shk1 + 585);
    const auto *shk1_588 = buffer.data(shk1 + 588);
    const auto *shk1_590 = buffer.data(shk1 + 590);
    const auto *shk1_593 = buffer.data(shk1 + 593);
    const auto *shk1_594 = buffer.data(shk1 + 594);
    const auto *shk1_596 = buffer.data(shk1 + 596);
    const auto *shk1_604 = buffer.data(shk1 + 604);
    const auto *shk1_606 = buffer.data(shk1 + 606);
    const auto *shk1_607 = buffer.data(shk1 + 607);
    const auto *shk1_608 = buffer.data(shk1 + 608);
    const auto *shk1_609 = buffer.data(shk1 + 609);
    const auto *shk1_611 = buffer.data(shk1 + 611);
    const auto *shk1_612 = buffer.data(shk1 + 612);
    const auto *shk1_615 = buffer.data(shk1 + 615);
    const auto *shk1_617 = buffer.data(shk1 + 617);
    const auto *shk1_618 = buffer.data(shk1 + 618);
    const auto *shk1_621 = buffer.data(shk1 + 621);
    const auto *shk1_622 = buffer.data(shk1 + 622);
    const auto *shk1_624 = buffer.data(shk1 + 624);
    const auto *shk1_626 = buffer.data(shk1 + 626);
    const auto *shk1_627 = buffer.data(shk1 + 627);
    const auto *shk1_629 = buffer.data(shk1 + 629);
    const auto *shk1_630 = buffer.data(shk1 + 630);
    const auto *shk1_632 = buffer.data(shk1 + 632);
    const auto *shk1_640 = buffer.data(shk1 + 640);
    const auto *shk1_642 = buffer.data(shk1 + 642);
    const auto *shk1_643 = buffer.data(shk1 + 643);
    const auto *shk1_644 = buffer.data(shk1 + 644);
    const auto *shk1_645 = buffer.data(shk1 + 645);
    const auto *shk1_647 = buffer.data(shk1 + 647);
    const auto *shk1_648 = buffer.data(shk1 + 648);
    const auto *shk1_651 = buffer.data(shk1 + 651);
    const auto *shk1_653 = buffer.data(shk1 + 653);
    const auto *shk1_654 = buffer.data(shk1 + 654);
    const auto *shk1_657 = buffer.data(shk1 + 657);
    const auto *shk1_658 = buffer.data(shk1 + 658);
    const auto *shk1_660 = buffer.data(shk1 + 660);
    const auto *shk1_662 = buffer.data(shk1 + 662);
    const auto *shk1_663 = buffer.data(shk1 + 663);
    const auto *shk1_665 = buffer.data(shk1 + 665);
    const auto *shk1_666 = buffer.data(shk1 + 666);
    const auto *shk1_668 = buffer.data(shk1 + 668);
    const auto *shk1_676 = buffer.data(shk1 + 676);
    const auto *shk1_678 = buffer.data(shk1 + 678);
    const auto *shk1_679 = buffer.data(shk1 + 679);
    const auto *shk1_680 = buffer.data(shk1 + 680);
    const auto *shk1_681 = buffer.data(shk1 + 681);
    const auto *shk1_683 = buffer.data(shk1 + 683);
    const auto *shk1_687 = buffer.data(shk1 + 687);

    const auto *sii_447 = buffer.data(sii + 447);
    const auto *sii_448 = buffer.data(sii + 448);
    const auto *sii_450 = buffer.data(sii + 450);
    const auto *sii_451 = buffer.data(sii + 451);
    const auto *sii_453 = buffer.data(sii + 453);
    const auto *sii_454 = buffer.data(sii + 454);
    const auto *sii_457 = buffer.data(sii + 457);
    const auto *sii_458 = buffer.data(sii + 458);
    const auto *sii_462 = buffer.data(sii + 462);
    const auto *sii_469 = buffer.data(sii + 469);
    const auto *sii_470 = buffer.data(sii + 470);
    const auto *sii_471 = buffer.data(sii + 471);
    const auto *sii_472 = buffer.data(sii + 472);
    const auto *sii_473 = buffer.data(sii + 473);
    const auto *sii_474 = buffer.data(sii + 474);
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
    const auto *sii_498 = buffer.data(sii + 498);
    const auto *sii_499 = buffer.data(sii + 499);
    const auto *sii_500 = buffer.data(sii + 500);
    const auto *sii_501 = buffer.data(sii + 501);
    const auto *sii_502 = buffer.data(sii + 502);
    const auto *sii_503 = buffer.data(sii + 503);
    const auto *sii_504 = buffer.data(sii + 504);
    const auto *sii_506 = buffer.data(sii + 506);
    const auto *sii_507 = buffer.data(sii + 507);
    const auto *sii_509 = buffer.data(sii + 509);
    const auto *sii_510 = buffer.data(sii + 510);
    const auto *sii_513 = buffer.data(sii + 513);
    const auto *sii_514 = buffer.data(sii + 514);
    const auto *sii_518 = buffer.data(sii + 518);
    const auto *sii_525 = buffer.data(sii + 525);
    const auto *sii_526 = buffer.data(sii + 526);
    const auto *sii_527 = buffer.data(sii + 527);
    const auto *sii_528 = buffer.data(sii + 528);
    const auto *sii_529 = buffer.data(sii + 529);
    const auto *sii_530 = buffer.data(sii + 530);
    const auto *sii_531 = buffer.data(sii + 531);
    const auto *sii_532 = buffer.data(sii + 532);
    const auto *sii_534 = buffer.data(sii + 534);

#pragma omp simd aligned(t_572, t_573, t_574, t_575, pb_x, pc_x, pc_y, shk0_572, shk0_573, \
                         shk0_575, shi_307, shk1_572, shk1_573, shk1_575, \
                         sii_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_572[k] = pb_x[k] * shk0_572[k]
                   - f_12 * pc_x[k] * shk1_572[k];

        t_573[k] = pb_x[k] * shk0_573[k]
                   - f_12 * pc_x[k] * shk1_573[k];

        t_574[k] = f_17 * shi_307[k]
                   + f_3 * pc_y[k] * sii_447[k];

        t_575[k] = pb_x[k] * shk0_575[k]
                   - f_12 * pc_x[k] * shk1_575[k];
    }

#pragma omp simd aligned(t_576, t_577, t_578, t_579, pb_z, pc_y, pc_z, shk0_360, shk0_363, \
                         shi_280, shi_308, shk1_360, shk1_363, \
                         sii_448 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_576[k] = pb_z[k] * shk0_360[k]
                   - f_12 * pc_z[k] * shk1_360[k];

        t_577[k] = f_16 * shi_308[k]
                   + f_3 * pc_y[k] * sii_448[k];

        t_578[k] = f_13 * shi_280[k]
                   + f_3 * pc_z[k] * sii_448[k];

        t_579[k] = pb_z[k] * shk0_363[k]
                   - f_12 * pc_z[k] * shk1_363[k];
    }

#pragma omp simd aligned(t_580, t_581, t_582, pb_x, pb_z, pc_x, pc_y, pc_z, shk0_366, \
                         shk0_581, shi_310, shi_453, shk1_366, shk1_581, \
                         sii_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_580[k] = f_16 * shi_310[k]
                   + f_3 * pc_y[k] * sii_450[k];

        t_581[k] = pb_x[k] * shk0_581[k]
                   + f_17 * shi_453[k]
                   - f_12 * pc_x[k] * shk1_581[k];

        t_582[k] = pb_z[k] * shk0_366[k]
                   - f_12 * pc_z[k] * shk1_366[k];
    }

#pragma omp simd aligned(t_583, t_584, t_585, pb_x, pc_x, pc_y, pc_z, shk0_585, shi_283, \
                         shi_313, shi_457, shk1_585, sii_451, sii_453 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_583[k] = f_13 * shi_283[k]
                   + f_3 * pc_z[k] * sii_451[k];

        t_584[k] = f_16 * shi_313[k]
                   + f_3 * pc_y[k] * sii_453[k];

        t_585[k] = pb_x[k] * shk0_585[k]
                   + f_16 * shi_457[k]
                   - f_12 * pc_x[k] * shk1_585[k];
    }

#pragma omp simd aligned(t_586, t_587, t_588, pb_x, pb_z, pc_x, pc_z, shk0_370, shk0_588, \
                         shi_286, shi_460, shk1_370, shk1_588, \
                         sii_454 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_586[k] = pb_z[k] * shk0_370[k]
                   - f_12 * pc_z[k] * shk1_370[k];

        t_587[k] = f_13 * shi_286[k]
                   + f_3 * pc_z[k] * sii_454[k];

        t_588[k] = pb_x[k] * shk0_588[k]
                   + f_15 * shi_460[k]
                   - f_12 * pc_x[k] * shk1_588[k];
    }

#pragma omp simd aligned(t_589, t_590, t_591, pb_x, pb_z, pc_x, pc_y, pc_z, shk0_375, \
                         shk0_590, shi_317, shi_462, shk1_375, shk1_590, \
                         sii_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_589[k] = f_16 * shi_317[k]
                   + f_3 * pc_y[k] * sii_457[k];

        t_590[k] = pb_x[k] * shk0_590[k]
                   + f_15 * shi_462[k]
                   - f_12 * pc_x[k] * shk1_590[k];

        t_591[k] = pb_z[k] * shk0_375[k]
                   - f_12 * pc_z[k] * shk1_375[k];
    }

#pragma omp simd aligned(t_592, t_593, t_594, pb_x, pc_x, pc_z, shk0_593, shk0_594, shi_290, \
                         shi_465, shi_466, shk1_593, shk1_594, \
                         sii_458 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_592[k] = f_13 * shi_290[k]
                   + f_3 * pc_z[k] * sii_458[k];

        t_593[k] = pb_x[k] * shk0_593[k]
                   + f_14 * shi_465[k]
                   - f_12 * pc_x[k] * shk1_593[k];

        t_594[k] = pb_x[k] * shk0_594[k]
                   + f_14 * shi_466[k]
                   - f_12 * pc_x[k] * shk1_594[k];
    }

#pragma omp simd aligned(t_595, t_596, t_597, t_598, pb_x, pc_x, pc_y, shk0_596, shi_322, \
                         shi_468, shi_469, shi_470, shk1_596, sii_462, sii_469, \
                         sii_470 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_595[k] = f_16 * shi_322[k]
                   + f_3 * pc_y[k] * sii_462[k];

        t_596[k] = pb_x[k] * shk0_596[k]
                   + f_14 * shi_468[k]
                   - f_12 * pc_x[k] * shk1_596[k];

        t_597[k] = f_13 * shi_469[k]
                   + f_3 * pc_x[k] * sii_469[k];

        t_598[k] = f_13 * shi_470[k]
                   + f_3 * pc_x[k] * sii_470[k];
    }

#pragma omp simd aligned(t_599, t_600, t_601, t_602, t_603, pc_x, shi_471, shi_472, shi_473, \
                         shi_474, shi_475, sii_471, sii_472, sii_473, sii_474, \
                         sii_475 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_599[k] = f_13 * shi_471[k]
                   + f_3 * pc_x[k] * sii_471[k];

        t_600[k] = f_13 * shi_472[k]
                   + f_3 * pc_x[k] * sii_472[k];

        t_601[k] = f_13 * shi_473[k]
                   + f_3 * pc_x[k] * sii_473[k];

        t_602[k] = f_13 * shi_474[k]
                   + f_3 * pc_x[k] * sii_474[k];

        t_603[k] = f_13 * shi_475[k]
                   + f_3 * pc_x[k] * sii_475[k];
    }

#pragma omp simd aligned(t_604, t_605, t_606, t_607, pb_x, pc_x, pc_z, shk0_604, shk0_606, \
                         shk0_607, shi_301, shk1_604, shk1_606, shk1_607, \
                         sii_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_604[k] = pb_x[k] * shk0_604[k]
                   - f_12 * pc_x[k] * shk1_604[k];

        t_605[k] = f_13 * shi_301[k]
                   + f_3 * pc_z[k] * sii_469[k];

        t_606[k] = pb_x[k] * shk0_606[k]
                   - f_12 * pc_x[k] * shk1_606[k];

        t_607[k] = pb_x[k] * shk0_607[k]
                   - f_12 * pc_x[k] * shk1_607[k];
    }

#pragma omp simd aligned(t_608, t_609, t_610, t_611, pb_x, pc_x, pc_y, shk0_608, shk0_609, \
                         shk0_611, shi_335, shk1_608, shk1_609, shk1_611, \
                         sii_475 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_608[k] = pb_x[k] * shk0_608[k]
                   - f_12 * pc_x[k] * shk1_608[k];

        t_609[k] = pb_x[k] * shk0_609[k]
                   - f_12 * pc_x[k] * shk1_609[k];

        t_610[k] = f_16 * shi_335[k]
                   + f_3 * pc_y[k] * sii_475[k];

        t_611[k] = pb_x[k] * shk0_611[k]
                   - f_12 * pc_x[k] * shk1_611[k];
    }

#pragma omp simd aligned(t_612, t_613, t_614, pb_x, pc_x, pc_y, pc_z, shk0_612, shi_308, \
                         shi_336, shi_476, shk1_612, sii_476 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_612[k] = pb_x[k] * shk0_612[k]
                   + f_18 * shi_476[k]
                   - f_12 * pc_x[k] * shk1_612[k];

        t_613[k] = f_15 * shi_336[k]
                   + f_3 * pc_y[k] * sii_476[k];

        t_614[k] = f_14 * shi_308[k]
                   + f_3 * pc_z[k] * sii_476[k];
    }

#pragma omp simd aligned(t_615, t_616, t_617, pb_x, pc_x, pc_y, shk0_615, shk0_617, shi_338, \
                         shi_479, shi_481, shk1_615, shk1_617, \
                         sii_478 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_615[k] = pb_x[k] * shk0_615[k]
                   + f_17 * shi_479[k]
                   - f_12 * pc_x[k] * shk1_615[k];

        t_616[k] = f_15 * shi_338[k]
                   + f_3 * pc_y[k] * sii_478[k];

        t_617[k] = pb_x[k] * shk0_617[k]
                   + f_17 * shi_481[k]
                   - f_12 * pc_x[k] * shk1_617[k];
    }

#pragma omp simd aligned(t_618, t_619, t_620, pb_x, pc_x, pc_y, pc_z, shk0_618, shi_311, \
                         shi_341, shi_482, shk1_618, sii_479, sii_481 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_618[k] = pb_x[k] * shk0_618[k]
                   + f_16 * shi_482[k]
                   - f_12 * pc_x[k] * shk1_618[k];

        t_619[k] = f_14 * shi_311[k]
                   + f_3 * pc_z[k] * sii_479[k];

        t_620[k] = f_15 * shi_341[k]
                   + f_3 * pc_y[k] * sii_481[k];
    }

#pragma omp simd aligned(t_621, t_622, t_623, pb_x, pc_x, pc_z, shk0_621, shk0_622, shi_314, \
                         shi_485, shi_486, shk1_621, shk1_622, \
                         sii_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_621[k] = pb_x[k] * shk0_621[k]
                   + f_16 * shi_485[k]
                   - f_12 * pc_x[k] * shk1_621[k];

        t_622[k] = pb_x[k] * shk0_622[k]
                   + f_15 * shi_486[k]
                   - f_12 * pc_x[k] * shk1_622[k];

        t_623[k] = f_14 * shi_314[k]
                   + f_3 * pc_z[k] * sii_482[k];
    }

#pragma omp simd aligned(t_624, t_625, t_626, pb_x, pc_x, pc_y, shk0_624, shk0_626, shi_345, \
                         shi_488, shi_490, shk1_624, shk1_626, \
                         sii_485 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_624[k] = pb_x[k] * shk0_624[k]
                   + f_15 * shi_488[k]
                   - f_12 * pc_x[k] * shk1_624[k];

        t_625[k] = f_15 * shi_345[k]
                   + f_3 * pc_y[k] * sii_485[k];

        t_626[k] = pb_x[k] * shk0_626[k]
                   + f_15 * shi_490[k]
                   - f_12 * pc_x[k] * shk1_626[k];
    }

#pragma omp simd aligned(t_627, t_628, t_629, pb_x, pc_x, pc_z, shk0_627, shk0_629, shi_318, \
                         shi_491, shi_493, shk1_627, shk1_629, \
                         sii_486 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_627[k] = pb_x[k] * shk0_627[k]
                   + f_14 * shi_491[k]
                   - f_12 * pc_x[k] * shk1_627[k];

        t_628[k] = f_14 * shi_318[k]
                   + f_3 * pc_z[k] * sii_486[k];

        t_629[k] = pb_x[k] * shk0_629[k]
                   + f_14 * shi_493[k]
                   - f_12 * pc_x[k] * shk1_629[k];
    }

#pragma omp simd aligned(t_630, t_631, t_632, pb_x, pc_x, pc_y, shk0_630, shk0_632, shi_350, \
                         shi_494, shi_496, shk1_630, shk1_632, \
                         sii_490 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_630[k] = pb_x[k] * shk0_630[k]
                   + f_14 * shi_494[k]
                   - f_12 * pc_x[k] * shk1_630[k];

        t_631[k] = f_15 * shi_350[k]
                   + f_3 * pc_y[k] * sii_490[k];

        t_632[k] = pb_x[k] * shk0_632[k]
                   + f_14 * shi_496[k]
                   - f_12 * pc_x[k] * shk1_632[k];
    }

#pragma omp simd aligned(t_633, t_634, t_635, t_636, t_637, pc_x, shi_497, shi_498, shi_499, \
                         shi_500, shi_501, sii_497, sii_498, sii_499, sii_500, \
                         sii_501 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_633[k] = f_13 * shi_497[k]
                   + f_3 * pc_x[k] * sii_497[k];

        t_634[k] = f_13 * shi_498[k]
                   + f_3 * pc_x[k] * sii_498[k];

        t_635[k] = f_13 * shi_499[k]
                   + f_3 * pc_x[k] * sii_499[k];

        t_636[k] = f_13 * shi_500[k]
                   + f_3 * pc_x[k] * sii_500[k];

        t_637[k] = f_13 * shi_501[k]
                   + f_3 * pc_x[k] * sii_501[k];
    }

#pragma omp simd aligned(t_638, t_639, t_640, t_641, pb_x, pc_x, pc_z, shk0_640, shi_329, \
                         shi_502, shi_503, shk1_640, sii_497, sii_502, \
                         sii_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_638[k] = f_13 * shi_502[k]
                   + f_3 * pc_x[k] * sii_502[k];

        t_639[k] = f_13 * shi_503[k]
                   + f_3 * pc_x[k] * sii_503[k];

        t_640[k] = pb_x[k] * shk0_640[k]
                   - f_12 * pc_x[k] * shk1_640[k];

        t_641[k] = f_14 * shi_329[k]
                   + f_3 * pc_z[k] * sii_497[k];
    }

#pragma omp simd aligned(t_642, t_643, t_644, t_645, pb_x, pc_x, shk0_642, shk0_643, shk0_644, \
                         shk0_645, shk1_642, shk1_643, shk1_644, \
                         shk1_645 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_642[k] = pb_x[k] * shk0_642[k]
                   - f_12 * pc_x[k] * shk1_642[k];

        t_643[k] = pb_x[k] * shk0_643[k]
                   - f_12 * pc_x[k] * shk1_643[k];

        t_644[k] = pb_x[k] * shk0_644[k]
                   - f_12 * pc_x[k] * shk1_644[k];

        t_645[k] = pb_x[k] * shk0_645[k]
                   - f_12 * pc_x[k] * shk1_645[k];
    }

#pragma omp simd aligned(t_646, t_647, t_648, t_649, pb_x, pc_x, pc_y, shk0_647, shk0_648, \
                         shi_363, shi_364, shi_504, shk1_647, shk1_648, sii_503, \
                         sii_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_646[k] = f_15 * shi_363[k]
                   + f_3 * pc_y[k] * sii_503[k];

        t_647[k] = pb_x[k] * shk0_647[k]
                   - f_12 * pc_x[k] * shk1_647[k];

        t_648[k] = pb_x[k] * shk0_648[k]
                   + f_18 * shi_504[k]
                   - f_12 * pc_x[k] * shk1_648[k];

        t_649[k] = f_14 * shi_364[k]
                   + f_3 * pc_y[k] * sii_504[k];
    }

#pragma omp simd aligned(t_650, t_651, t_652, pb_x, pc_x, pc_y, pc_z, shk0_651, shi_336, \
                         shi_366, shi_507, shk1_651, sii_504, sii_506 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_650[k] = f_15 * shi_336[k]
                   + f_3 * pc_z[k] * sii_504[k];

        t_651[k] = pb_x[k] * shk0_651[k]
                   + f_17 * shi_507[k]
                   - f_12 * pc_x[k] * shk1_651[k];

        t_652[k] = f_14 * shi_366[k]
                   + f_3 * pc_y[k] * sii_506[k];
    }

#pragma omp simd aligned(t_653, t_654, t_655, pb_x, pc_x, pc_z, shk0_653, shk0_654, shi_339, \
                         shi_509, shi_510, shk1_653, shk1_654, \
                         sii_507 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_653[k] = pb_x[k] * shk0_653[k]
                   + f_17 * shi_509[k]
                   - f_12 * pc_x[k] * shk1_653[k];

        t_654[k] = pb_x[k] * shk0_654[k]
                   + f_16 * shi_510[k]
                   - f_12 * pc_x[k] * shk1_654[k];

        t_655[k] = f_15 * shi_339[k]
                   + f_3 * pc_z[k] * sii_507[k];
    }

#pragma omp simd aligned(t_656, t_657, t_658, pb_x, pc_x, pc_y, shk0_657, shk0_658, shi_369, \
                         shi_513, shi_514, shk1_657, shk1_658, \
                         sii_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_656[k] = f_14 * shi_369[k]
                   + f_3 * pc_y[k] * sii_509[k];

        t_657[k] = pb_x[k] * shk0_657[k]
                   + f_16 * shi_513[k]
                   - f_12 * pc_x[k] * shk1_657[k];

        t_658[k] = pb_x[k] * shk0_658[k]
                   + f_15 * shi_514[k]
                   - f_12 * pc_x[k] * shk1_658[k];
    }

#pragma omp simd aligned(t_659, t_660, t_661, pb_x, pc_x, pc_y, pc_z, shk0_660, shi_342, \
                         shi_373, shi_516, shk1_660, sii_510, sii_513 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_659[k] = f_15 * shi_342[k]
                   + f_3 * pc_z[k] * sii_510[k];

        t_660[k] = pb_x[k] * shk0_660[k]
                   + f_15 * shi_516[k]
                   - f_12 * pc_x[k] * shk1_660[k];

        t_661[k] = f_14 * shi_373[k]
                   + f_3 * pc_y[k] * sii_513[k];
    }

#pragma omp simd aligned(t_662, t_663, t_664, pb_x, pc_x, pc_z, shk0_662, shk0_663, shi_346, \
                         shi_518, shi_519, shk1_662, shk1_663, \
                         sii_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_662[k] = pb_x[k] * shk0_662[k]
                   + f_15 * shi_518[k]
                   - f_12 * pc_x[k] * shk1_662[k];

        t_663[k] = pb_x[k] * shk0_663[k]
                   + f_14 * shi_519[k]
                   - f_12 * pc_x[k] * shk1_663[k];

        t_664[k] = f_15 * shi_346[k]
                   + f_3 * pc_z[k] * sii_514[k];
    }

#pragma omp simd aligned(t_665, t_666, t_667, pb_x, pc_x, pc_y, shk0_665, shk0_666, shi_378, \
                         shi_521, shi_522, shk1_665, shk1_666, \
                         sii_518 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_665[k] = pb_x[k] * shk0_665[k]
                   + f_14 * shi_521[k]
                   - f_12 * pc_x[k] * shk1_665[k];

        t_666[k] = pb_x[k] * shk0_666[k]
                   + f_14 * shi_522[k]
                   - f_12 * pc_x[k] * shk1_666[k];

        t_667[k] = f_14 * shi_378[k]
                   + f_3 * pc_y[k] * sii_518[k];
    }

#pragma omp simd aligned(t_668, t_669, t_670, t_671, pb_x, pc_x, shk0_668, shi_524, shi_525, \
                         shi_526, shi_527, shk1_668, sii_525, sii_526, \
                         sii_527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_668[k] = pb_x[k] * shk0_668[k]
                   + f_14 * shi_524[k]
                   - f_12 * pc_x[k] * shk1_668[k];

        t_669[k] = f_13 * shi_525[k]
                   + f_3 * pc_x[k] * sii_525[k];

        t_670[k] = f_13 * shi_526[k]
                   + f_3 * pc_x[k] * sii_526[k];

        t_671[k] = f_13 * shi_527[k]
                   + f_3 * pc_x[k] * sii_527[k];
    }

#pragma omp simd aligned(t_672, t_673, t_674, t_675, pc_x, shi_528, shi_529, shi_530, shi_531, \
                         sii_528, sii_529, sii_530, sii_531 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_672[k] = f_13 * shi_528[k]
                   + f_3 * pc_x[k] * sii_528[k];

        t_673[k] = f_13 * shi_529[k]
                   + f_3 * pc_x[k] * sii_529[k];

        t_674[k] = f_13 * shi_530[k]
                   + f_3 * pc_x[k] * sii_530[k];

        t_675[k] = f_13 * shi_531[k]
                   + f_3 * pc_x[k] * sii_531[k];
    }

#pragma omp simd aligned(t_676, t_677, t_678, t_679, pb_x, pc_x, pc_z, shk0_676, shk0_678, \
                         shk0_679, shi_357, shk1_676, shk1_678, shk1_679, \
                         sii_525 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_676[k] = pb_x[k] * shk0_676[k]
                   - f_12 * pc_x[k] * shk1_676[k];

        t_677[k] = f_15 * shi_357[k]
                   + f_3 * pc_z[k] * sii_525[k];

        t_678[k] = pb_x[k] * shk0_678[k]
                   - f_12 * pc_x[k] * shk1_678[k];

        t_679[k] = pb_x[k] * shk0_679[k]
                   - f_12 * pc_x[k] * shk1_679[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, t_683, pb_x, pc_x, pc_y, shk0_680, shk0_681, \
                         shk0_683, shi_391, shk1_680, shk1_681, shk1_683, \
                         sii_531 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = pb_x[k] * shk0_680[k]
                   - f_12 * pc_x[k] * shk1_680[k];

        t_681[k] = pb_x[k] * shk0_681[k]
                   - f_12 * pc_x[k] * shk1_681[k];

        t_682[k] = f_14 * shi_391[k]
                   + f_3 * pc_y[k] * sii_531[k];

        t_683[k] = pb_x[k] * shk0_683[k]
                   - f_12 * pc_x[k] * shk1_683[k];
    }

#pragma omp simd aligned(t_684, t_685, t_686, pb_y, pc_y, pc_z, shk0_504, shi_364, shi_392, \
                         shk1_504, sii_532 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_684[k] = pb_y[k] * shk0_504[k]
                   - f_12 * pc_y[k] * shk1_504[k];

        t_685[k] = f_13 * shi_392[k]
                   + f_3 * pc_y[k] * sii_532[k];

        t_686[k] = f_16 * shi_364[k]
                   + f_3 * pc_z[k] * sii_532[k];
    }

#pragma omp simd aligned(t_687, t_688, t_689, pb_x, pb_y, pc_x, pc_y, shk0_509, shk0_687, \
                         shi_394, shi_535, shk1_509, shk1_687, \
                         sii_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_687[k] = pb_x[k] * shk0_687[k]
                   + f_17 * shi_535[k]
                   - f_12 * pc_x[k] * shk1_687[k];

        t_688[k] = f_13 * shi_394[k]
                   + f_3 * pc_y[k] * sii_534[k];

        t_689[k] = pb_y[k] * shk0_509[k]
                   - f_12 * pc_y[k] * shk1_509[k];
    }
}

static auto
compute_prim_sik_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t shk0,
                                                          const size_t shi, const size_t shk1,
                                                          const size_t sih0, const size_t sih1,
                                                          const size_t sii, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
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
    const auto f_18 = 3.5 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *shk0_513 = buffer.data(shk0 + 513);
    const auto *shk0_518 = buffer.data(shk0 + 518);
    const auto *shk0_524 = buffer.data(shk0 + 524);
    const auto *shk0_540 = buffer.data(shk0 + 540);
    const auto *shk0_543 = buffer.data(shk0 + 543);
    const auto *shk0_546 = buffer.data(shk0 + 546);
    const auto *shk0_550 = buffer.data(shk0 + 550);
    const auto *shk0_555 = buffer.data(shk0 + 555);
    const auto *shk0_690 = buffer.data(shk0 + 690);
    const auto *shk0_694 = buffer.data(shk0 + 694);
    const auto *shk0_696 = buffer.data(shk0 + 696);
    const auto *shk0_699 = buffer.data(shk0 + 699);
    const auto *shk0_701 = buffer.data(shk0 + 701);
    const auto *shk0_702 = buffer.data(shk0 + 702);
    const auto *shk0_712 = buffer.data(shk0 + 712);
    const auto *shk0_714 = buffer.data(shk0 + 714);
    const auto *shk0_715 = buffer.data(shk0 + 715);
    const auto *shk0_716 = buffer.data(shk0 + 716);
    const auto *shk0_717 = buffer.data(shk0 + 717);
    const auto *shk0_719 = buffer.data(shk0 + 719);
    const auto *shk0_720 = buffer.data(shk0 + 720);
    const auto *shk0_723 = buffer.data(shk0 + 723);
    const auto *shk0_725 = buffer.data(shk0 + 725);
    const auto *shk0_726 = buffer.data(shk0 + 726);
    const auto *shk0_729 = buffer.data(shk0 + 729);
    const auto *shk0_730 = buffer.data(shk0 + 730);
    const auto *shk0_732 = buffer.data(shk0 + 732);
    const auto *shk0_734 = buffer.data(shk0 + 734);
    const auto *shk0_735 = buffer.data(shk0 + 735);
    const auto *shk0_737 = buffer.data(shk0 + 737);
    const auto *shk0_738 = buffer.data(shk0 + 738);
    const auto *shk0_740 = buffer.data(shk0 + 740);
    const auto *shk0_748 = buffer.data(shk0 + 748);
    const auto *shk0_750 = buffer.data(shk0 + 750);
    const auto *shk0_751 = buffer.data(shk0 + 751);
    const auto *shk0_752 = buffer.data(shk0 + 752);
    const auto *shk0_753 = buffer.data(shk0 + 753);
    const auto *shk0_755 = buffer.data(shk0 + 755);

    const auto *shi_367 = buffer.data(shi + 367);
    const auto *shi_370 = buffer.data(shi + 370);
    const auto *shi_374 = buffer.data(shi + 374);
    const auto *shi_385 = buffer.data(shi + 385);
    const auto *shi_392 = buffer.data(shi + 392);
    const auto *shi_395 = buffer.data(shi + 395);
    const auto *shi_397 = buffer.data(shi + 397);
    const auto *shi_398 = buffer.data(shi + 398);
    const auto *shi_401 = buffer.data(shi + 401);
    const auto *shi_402 = buffer.data(shi + 402);
    const auto *shi_406 = buffer.data(shi + 406);
    const auto *shi_413 = buffer.data(shi + 413);
    const auto *shi_419 = buffer.data(shi + 419);
    const auto *shi_420 = buffer.data(shi + 420);
    const auto *shi_422 = buffer.data(shi + 422);
    const auto *shi_423 = buffer.data(shi + 423);
    const auto *shi_425 = buffer.data(shi + 425);
    const auto *shi_426 = buffer.data(shi + 426);
    const auto *shi_429 = buffer.data(shi + 429);
    const auto *shi_430 = buffer.data(shi + 430);
    const auto *shi_434 = buffer.data(shi + 434);
    const auto *shi_441 = buffer.data(shi + 441);
    const auto *shi_443 = buffer.data(shi + 443);
    const auto *shi_444 = buffer.data(shi + 444);
    const auto *shi_445 = buffer.data(shi + 445);
    const auto *shi_446 = buffer.data(shi + 446);
    const auto *shi_447 = buffer.data(shi + 447);
    const auto *shi_448 = buffer.data(shi + 448);
    const auto *shi_450 = buffer.data(shi + 450);
    const auto *shi_453 = buffer.data(shi + 453);
    const auto *shi_457 = buffer.data(shi + 457);
    const auto *shi_538 = buffer.data(shi + 538);
    const auto *shi_542 = buffer.data(shi + 542);
    const auto *shi_544 = buffer.data(shi + 544);
    const auto *shi_547 = buffer.data(shi + 547);
    const auto *shi_549 = buffer.data(shi + 549);
    const auto *shi_550 = buffer.data(shi + 550);
    const auto *shi_553 = buffer.data(shi + 553);
    const auto *shi_554 = buffer.data(shi + 554);
    const auto *shi_555 = buffer.data(shi + 555);
    const auto *shi_556 = buffer.data(shi + 556);
    const auto *shi_557 = buffer.data(shi + 557);
    const auto *shi_558 = buffer.data(shi + 558);
    const auto *shi_559 = buffer.data(shi + 559);
    const auto *shi_560 = buffer.data(shi + 560);
    const auto *shi_563 = buffer.data(shi + 563);
    const auto *shi_565 = buffer.data(shi + 565);
    const auto *shi_566 = buffer.data(shi + 566);
    const auto *shi_569 = buffer.data(shi + 569);
    const auto *shi_570 = buffer.data(shi + 570);
    const auto *shi_572 = buffer.data(shi + 572);
    const auto *shi_574 = buffer.data(shi + 574);
    const auto *shi_575 = buffer.data(shi + 575);
    const auto *shi_577 = buffer.data(shi + 577);
    const auto *shi_578 = buffer.data(shi + 578);
    const auto *shi_580 = buffer.data(shi + 580);
    const auto *shi_581 = buffer.data(shi + 581);
    const auto *shi_582 = buffer.data(shi + 582);
    const auto *shi_583 = buffer.data(shi + 583);
    const auto *shi_584 = buffer.data(shi + 584);
    const auto *shi_585 = buffer.data(shi + 585);
    const auto *shi_586 = buffer.data(shi + 586);
    const auto *shi_587 = buffer.data(shi + 587);

    const auto *shk1_513 = buffer.data(shk1 + 513);
    const auto *shk1_518 = buffer.data(shk1 + 518);
    const auto *shk1_524 = buffer.data(shk1 + 524);
    const auto *shk1_540 = buffer.data(shk1 + 540);
    const auto *shk1_543 = buffer.data(shk1 + 543);
    const auto *shk1_546 = buffer.data(shk1 + 546);
    const auto *shk1_550 = buffer.data(shk1 + 550);
    const auto *shk1_555 = buffer.data(shk1 + 555);
    const auto *shk1_690 = buffer.data(shk1 + 690);
    const auto *shk1_694 = buffer.data(shk1 + 694);
    const auto *shk1_696 = buffer.data(shk1 + 696);
    const auto *shk1_699 = buffer.data(shk1 + 699);
    const auto *shk1_701 = buffer.data(shk1 + 701);
    const auto *shk1_702 = buffer.data(shk1 + 702);
    const auto *shk1_712 = buffer.data(shk1 + 712);
    const auto *shk1_714 = buffer.data(shk1 + 714);
    const auto *shk1_715 = buffer.data(shk1 + 715);
    const auto *shk1_716 = buffer.data(shk1 + 716);
    const auto *shk1_717 = buffer.data(shk1 + 717);
    const auto *shk1_719 = buffer.data(shk1 + 719);
    const auto *shk1_720 = buffer.data(shk1 + 720);
    const auto *shk1_723 = buffer.data(shk1 + 723);
    const auto *shk1_725 = buffer.data(shk1 + 725);
    const auto *shk1_726 = buffer.data(shk1 + 726);
    const auto *shk1_729 = buffer.data(shk1 + 729);
    const auto *shk1_730 = buffer.data(shk1 + 730);
    const auto *shk1_732 = buffer.data(shk1 + 732);
    const auto *shk1_734 = buffer.data(shk1 + 734);
    const auto *shk1_735 = buffer.data(shk1 + 735);
    const auto *shk1_737 = buffer.data(shk1 + 737);
    const auto *shk1_738 = buffer.data(shk1 + 738);
    const auto *shk1_740 = buffer.data(shk1 + 740);
    const auto *shk1_748 = buffer.data(shk1 + 748);
    const auto *shk1_750 = buffer.data(shk1 + 750);
    const auto *shk1_751 = buffer.data(shk1 + 751);
    const auto *shk1_752 = buffer.data(shk1 + 752);
    const auto *shk1_753 = buffer.data(shk1 + 753);
    const auto *shk1_755 = buffer.data(shk1 + 755);

    const auto *sih0_441 = buffer.data(sih0 + 441);
    const auto *sih0_444 = buffer.data(sih0 + 444);
    const auto *sih0_446 = buffer.data(sih0 + 446);
    const auto *sih0_447 = buffer.data(sih0 + 447);
    const auto *sih0_450 = buffer.data(sih0 + 450);
    const auto *sih0_451 = buffer.data(sih0 + 451);
    const auto *sih0_453 = buffer.data(sih0 + 453);
    const auto *sih0_455 = buffer.data(sih0 + 455);
    const auto *sih0_456 = buffer.data(sih0 + 456);
    const auto *sih0_458 = buffer.data(sih0 + 458);
    const auto *sih0_459 = buffer.data(sih0 + 459);
    const auto *sih0_460 = buffer.data(sih0 + 460);
    const auto *sih0_461 = buffer.data(sih0 + 461);
    const auto *sih0_467 = buffer.data(sih0 + 467);
    const auto *sih0_471 = buffer.data(sih0 + 471);
    const auto *sih0_474 = buffer.data(sih0 + 474);
    const auto *sih0_476 = buffer.data(sih0 + 476);
    const auto *sih0_479 = buffer.data(sih0 + 479);
    const auto *sih0_480 = buffer.data(sih0 + 480);

    const auto *sih1_441 = buffer.data(sih1 + 441);
    const auto *sih1_444 = buffer.data(sih1 + 444);
    const auto *sih1_446 = buffer.data(sih1 + 446);
    const auto *sih1_447 = buffer.data(sih1 + 447);
    const auto *sih1_450 = buffer.data(sih1 + 450);
    const auto *sih1_451 = buffer.data(sih1 + 451);
    const auto *sih1_453 = buffer.data(sih1 + 453);
    const auto *sih1_455 = buffer.data(sih1 + 455);
    const auto *sih1_456 = buffer.data(sih1 + 456);
    const auto *sih1_458 = buffer.data(sih1 + 458);
    const auto *sih1_459 = buffer.data(sih1 + 459);
    const auto *sih1_460 = buffer.data(sih1 + 460);
    const auto *sih1_461 = buffer.data(sih1 + 461);
    const auto *sih1_467 = buffer.data(sih1 + 467);
    const auto *sih1_471 = buffer.data(sih1 + 471);
    const auto *sih1_474 = buffer.data(sih1 + 474);
    const auto *sih1_476 = buffer.data(sih1 + 476);
    const auto *sih1_479 = buffer.data(sih1 + 479);
    const auto *sih1_480 = buffer.data(sih1 + 480);

    const auto *sii_535 = buffer.data(sii + 535);
    const auto *sii_537 = buffer.data(sii + 537);
    const auto *sii_538 = buffer.data(sii + 538);
    const auto *sii_541 = buffer.data(sii + 541);
    const auto *sii_542 = buffer.data(sii + 542);
    const auto *sii_546 = buffer.data(sii + 546);
    const auto *sii_553 = buffer.data(sii + 553);
    const auto *sii_554 = buffer.data(sii + 554);
    const auto *sii_555 = buffer.data(sii + 555);
    const auto *sii_556 = buffer.data(sii + 556);
    const auto *sii_557 = buffer.data(sii + 557);
    const auto *sii_558 = buffer.data(sii + 558);
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
    const auto *sii_582 = buffer.data(sii + 582);
    const auto *sii_583 = buffer.data(sii + 583);
    const auto *sii_584 = buffer.data(sii + 584);
    const auto *sii_585 = buffer.data(sii + 585);
    const auto *sii_586 = buffer.data(sii + 586);
    const auto *sii_587 = buffer.data(sii + 587);
    const auto *sii_588 = buffer.data(sii + 588);
    const auto *sii_590 = buffer.data(sii + 590);
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
    const auto *sii_616 = buffer.data(sii + 616);
    const auto *sii_618 = buffer.data(sii + 618);
    const auto *sii_619 = buffer.data(sii + 619);
    const auto *sii_621 = buffer.data(sii + 621);
    const auto *sii_622 = buffer.data(sii + 622);
    const auto *sii_625 = buffer.data(sii + 625);
    const auto *sii_626 = buffer.data(sii + 626);
    const auto *sii_628 = buffer.data(sii + 628);
    const auto *sii_630 = buffer.data(sii + 630);
    const auto *sii_633 = buffer.data(sii + 633);
    const auto *sii_634 = buffer.data(sii + 634);

#pragma omp simd aligned(t_690, t_691, t_692, pb_x, pc_x, pc_y, pc_z, shk0_690, shi_367, \
                         shi_397, shi_538, shk1_690, sii_535, sii_537 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_690[k] = pb_x[k] * shk0_690[k]
                   + f_16 * shi_538[k]
                   - f_12 * pc_x[k] * shk1_690[k];

        t_691[k] = f_16 * shi_367[k]
                   + f_3 * pc_z[k] * sii_535[k];

        t_692[k] = f_13 * shi_397[k]
                   + f_3 * pc_y[k] * sii_537[k];
    }

#pragma omp simd aligned(t_693, t_694, t_695, pb_x, pb_y, pc_x, pc_y, pc_z, shk0_513, \
                         shk0_694, shi_370, shi_542, shk1_513, shk1_694, \
                         sii_538 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_693[k] = pb_y[k] * shk0_513[k]
                   - f_12 * pc_y[k] * shk1_513[k];

        t_694[k] = pb_x[k] * shk0_694[k]
                   + f_15 * shi_542[k]
                   - f_12 * pc_x[k] * shk1_694[k];

        t_695[k] = f_16 * shi_370[k]
                   + f_3 * pc_z[k] * sii_538[k];
    }

#pragma omp simd aligned(t_696, t_697, t_698, pb_x, pb_y, pc_x, pc_y, shk0_518, shk0_696, \
                         shi_401, shi_544, shk1_518, shk1_696, \
                         sii_541 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_696[k] = pb_x[k] * shk0_696[k]
                   + f_15 * shi_544[k]
                   - f_12 * pc_x[k] * shk1_696[k];

        t_697[k] = f_13 * shi_401[k]
                   + f_3 * pc_y[k] * sii_541[k];

        t_698[k] = pb_y[k] * shk0_518[k]
                   - f_12 * pc_y[k] * shk1_518[k];
    }

#pragma omp simd aligned(t_699, t_700, t_701, pb_x, pc_x, pc_z, shk0_699, shk0_701, shi_374, \
                         shi_547, shi_549, shk1_699, shk1_701, \
                         sii_542 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_699[k] = pb_x[k] * shk0_699[k]
                   + f_14 * shi_547[k]
                   - f_12 * pc_x[k] * shk1_699[k];

        t_700[k] = f_16 * shi_374[k]
                   + f_3 * pc_z[k] * sii_542[k];

        t_701[k] = pb_x[k] * shk0_701[k]
                   + f_14 * shi_549[k]
                   - f_12 * pc_x[k] * shk1_701[k];
    }

#pragma omp simd aligned(t_702, t_703, t_704, pb_x, pb_y, pc_x, pc_y, shk0_524, shk0_702, \
                         shi_406, shi_550, shk1_524, shk1_702, \
                         sii_546 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_702[k] = pb_x[k] * shk0_702[k]
                   + f_14 * shi_550[k]
                   - f_12 * pc_x[k] * shk1_702[k];

        t_703[k] = f_13 * shi_406[k]
                   + f_3 * pc_y[k] * sii_546[k];

        t_704[k] = pb_y[k] * shk0_524[k]
                   - f_12 * pc_y[k] * shk1_524[k];
    }

#pragma omp simd aligned(t_705, t_706, t_707, t_708, t_709, pc_x, shi_553, shi_554, shi_555, \
                         shi_556, shi_557, sii_553, sii_554, sii_555, sii_556, \
                         sii_557 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_705[k] = f_13 * shi_553[k]
                   + f_3 * pc_x[k] * sii_553[k];

        t_706[k] = f_13 * shi_554[k]
                   + f_3 * pc_x[k] * sii_554[k];

        t_707[k] = f_13 * shi_555[k]
                   + f_3 * pc_x[k] * sii_555[k];

        t_708[k] = f_13 * shi_556[k]
                   + f_3 * pc_x[k] * sii_556[k];

        t_709[k] = f_13 * shi_557[k]
                   + f_3 * pc_x[k] * sii_557[k];
    }

#pragma omp simd aligned(t_710, t_711, t_712, t_713, pb_x, pc_x, pc_z, shk0_712, shi_385, \
                         shi_558, shi_559, shk1_712, sii_553, sii_558, \
                         sii_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_710[k] = f_13 * shi_558[k]
                   + f_3 * pc_x[k] * sii_558[k];

        t_711[k] = f_13 * shi_559[k]
                   + f_3 * pc_x[k] * sii_559[k];

        t_712[k] = pb_x[k] * shk0_712[k]
                   - f_12 * pc_x[k] * shk1_712[k];

        t_713[k] = f_16 * shi_385[k]
                   + f_3 * pc_z[k] * sii_553[k];
    }

#pragma omp simd aligned(t_714, t_715, t_716, t_717, pb_x, pc_x, shk0_714, shk0_715, shk0_716, \
                         shk0_717, shk1_714, shk1_715, shk1_716, \
                         shk1_717 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_714[k] = pb_x[k] * shk0_714[k]
                   - f_12 * pc_x[k] * shk1_714[k];

        t_715[k] = pb_x[k] * shk0_715[k]
                   - f_12 * pc_x[k] * shk1_715[k];

        t_716[k] = pb_x[k] * shk0_716[k]
                   - f_12 * pc_x[k] * shk1_716[k];

        t_717[k] = pb_x[k] * shk0_717[k]
                   - f_12 * pc_x[k] * shk1_717[k];
    }

#pragma omp simd aligned(t_718, t_719, t_720, t_721, pb_x, pc_x, pc_y, shk0_719, shk0_720, \
                         shi_419, shi_560, shk1_719, shk1_720, sii_559, \
                         sii_560 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_718[k] = f_13 * shi_419[k]
                   + f_3 * pc_y[k] * sii_559[k];

        t_719[k] = pb_x[k] * shk0_719[k]
                   - f_12 * pc_x[k] * shk1_719[k];

        t_720[k] = pb_x[k] * shk0_720[k]
                   + f_18 * shi_560[k]
                   - f_12 * pc_x[k] * shk1_720[k];

        t_721[k] = f_3 * pc_y[k] * sii_560[k];
    }

#pragma omp simd aligned(t_722, t_723, t_724, pb_x, pc_x, pc_y, pc_z, shk0_723, shi_392, \
                         shi_563, shk1_723, sii_560, sii_562 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_722[k] = f_17 * shi_392[k]
                   + f_3 * pc_z[k] * sii_560[k];

        t_723[k] = pb_x[k] * shk0_723[k]
                   + f_17 * shi_563[k]
                   - f_12 * pc_x[k] * shk1_723[k];

        t_724[k] = f_3 * pc_y[k] * sii_562[k];
    }

#pragma omp simd aligned(t_725, t_726, t_727, pb_x, pc_x, pc_z, shk0_725, shk0_726, shi_395, \
                         shi_565, shi_566, shk1_725, shk1_726, \
                         sii_563 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_725[k] = pb_x[k] * shk0_725[k]
                   + f_17 * shi_565[k]
                   - f_12 * pc_x[k] * shk1_725[k];

        t_726[k] = pb_x[k] * shk0_726[k]
                   + f_16 * shi_566[k]
                   - f_12 * pc_x[k] * shk1_726[k];

        t_727[k] = f_17 * shi_395[k]
                   + f_3 * pc_z[k] * sii_563[k];
    }

#pragma omp simd aligned(t_728, t_729, t_730, pb_x, pc_x, pc_y, shk0_729, shk0_730, shi_569, \
                         shi_570, shk1_729, shk1_730, sii_565 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_728[k] = f_3 * pc_y[k] * sii_565[k];

        t_729[k] = pb_x[k] * shk0_729[k]
                   + f_16 * shi_569[k]
                   - f_12 * pc_x[k] * shk1_729[k];

        t_730[k] = pb_x[k] * shk0_730[k]
                   + f_15 * shi_570[k]
                   - f_12 * pc_x[k] * shk1_730[k];
    }

#pragma omp simd aligned(t_731, t_732, t_733, pb_x, pc_x, pc_y, pc_z, shk0_732, shi_398, \
                         shi_572, shk1_732, sii_566, sii_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_731[k] = f_17 * shi_398[k]
                   + f_3 * pc_z[k] * sii_566[k];

        t_732[k] = pb_x[k] * shk0_732[k]
                   + f_15 * shi_572[k]
                   - f_12 * pc_x[k] * shk1_732[k];

        t_733[k] = f_3 * pc_y[k] * sii_569[k];
    }

#pragma omp simd aligned(t_734, t_735, t_736, pb_x, pc_x, pc_z, shk0_734, shk0_735, shi_402, \
                         shi_574, shi_575, shk1_734, shk1_735, \
                         sii_570 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_734[k] = pb_x[k] * shk0_734[k]
                   + f_15 * shi_574[k]
                   - f_12 * pc_x[k] * shk1_734[k];

        t_735[k] = pb_x[k] * shk0_735[k]
                   + f_14 * shi_575[k]
                   - f_12 * pc_x[k] * shk1_735[k];

        t_736[k] = f_17 * shi_402[k]
                   + f_3 * pc_z[k] * sii_570[k];
    }

#pragma omp simd aligned(t_737, t_738, t_739, pb_x, pc_x, pc_y, shk0_737, shk0_738, shi_577, \
                         shi_578, shk1_737, shk1_738, sii_574 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_737[k] = pb_x[k] * shk0_737[k]
                   + f_14 * shi_577[k]
                   - f_12 * pc_x[k] * shk1_737[k];

        t_738[k] = pb_x[k] * shk0_738[k]
                   + f_14 * shi_578[k]
                   - f_12 * pc_x[k] * shk1_738[k];

        t_739[k] = f_3 * pc_y[k] * sii_574[k];
    }

#pragma omp simd aligned(t_740, t_741, t_742, t_743, pb_x, pc_x, shk0_740, shi_580, shi_581, \
                         shi_582, shi_583, shk1_740, sii_581, sii_582, \
                         sii_583 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_740[k] = pb_x[k] * shk0_740[k]
                   + f_14 * shi_580[k]
                   - f_12 * pc_x[k] * shk1_740[k];

        t_741[k] = f_13 * shi_581[k]
                   + f_3 * pc_x[k] * sii_581[k];

        t_742[k] = f_13 * shi_582[k]
                   + f_3 * pc_x[k] * sii_582[k];

        t_743[k] = f_13 * shi_583[k]
                   + f_3 * pc_x[k] * sii_583[k];
    }

#pragma omp simd aligned(t_744, t_745, t_746, t_747, pc_x, shi_584, shi_585, shi_586, shi_587, \
                         sii_584, sii_585, sii_586, sii_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_744[k] = f_13 * shi_584[k]
                   + f_3 * pc_x[k] * sii_584[k];

        t_745[k] = f_13 * shi_585[k]
                   + f_3 * pc_x[k] * sii_585[k];

        t_746[k] = f_13 * shi_586[k]
                   + f_3 * pc_x[k] * sii_586[k];

        t_747[k] = f_13 * shi_587[k]
                   + f_3 * pc_x[k] * sii_587[k];
    }

#pragma omp simd aligned(t_748, t_749, t_750, t_751, pb_x, pc_x, pc_z, shk0_748, shk0_750, \
                         shk0_751, shi_413, shk1_748, shk1_750, shk1_751, \
                         sii_581 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_748[k] = pb_x[k] * shk0_748[k]
                   - f_12 * pc_x[k] * shk1_748[k];

        t_749[k] = f_17 * shi_413[k]
                   + f_3 * pc_z[k] * sii_581[k];

        t_750[k] = pb_x[k] * shk0_750[k]
                   - f_12 * pc_x[k] * shk1_750[k];

        t_751[k] = pb_x[k] * shk0_751[k]
                   - f_12 * pc_x[k] * shk1_751[k];
    }

#pragma omp simd aligned(t_752, t_753, t_754, t_755, pb_x, pc_x, pc_y, shk0_752, shk0_753, \
                         shk0_755, shk1_752, shk1_753, shk1_755, \
                         sii_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_752[k] = pb_x[k] * shk0_752[k]
                   - f_12 * pc_x[k] * shk1_752[k];

        t_753[k] = pb_x[k] * shk0_753[k]
                   - f_12 * pc_x[k] * shk1_753[k];

        t_754[k] = f_3 * pc_y[k] * sii_587[k];

        t_755[k] = pb_x[k] * shk0_755[k]
                   - f_12 * pc_x[k] * shk1_755[k];
    }

#pragma omp simd aligned(t_756, t_757, t_758, t_759, pc_x, pc_y, pc_z, shi_420, sih0_441, \
                         sih0_444, sih1_441, sih1_444, sii_588, \
                         sii_591 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_756[k] = f_1 * sih0_441[k]
                   - f_2 * sih1_441[k]
                   + f_3 * pc_x[k] * sii_588[k];

        t_757[k] = f_0 * shi_420[k]
                   + f_3 * pc_y[k] * sii_588[k];

        t_758[k] = f_3 * pc_z[k] * sii_588[k];

        t_759[k] = f_4 * sih0_444[k]
                   - f_5 * sih1_444[k]
                   + f_3 * pc_x[k] * sii_591[k];
    }

#pragma omp simd aligned(t_760, t_761, t_762, t_763, pc_x, pc_y, pc_z, shi_422, sih0_446, \
                         sih0_447, sih1_446, sih1_447, sii_590, sii_591, sii_593, \
                         sii_594 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_760[k] = f_0 * shi_422[k]
                   + f_3 * pc_y[k] * sii_590[k];

        t_761[k] = f_4 * sih0_446[k]
                   - f_5 * sih1_446[k]
                   + f_3 * pc_x[k] * sii_593[k];

        t_762[k] = f_6 * sih0_447[k]
                   - f_7 * sih1_447[k]
                   + f_3 * pc_x[k] * sii_594[k];

        t_763[k] = f_3 * pc_z[k] * sii_591[k];
    }

#pragma omp simd aligned(t_764, t_765, t_766, t_767, pc_x, pc_y, pc_z, shi_425, sih0_450, \
                         sih0_451, sih1_450, sih1_451, sii_593, sii_594, sii_597, \
                         sii_598 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_764[k] = f_0 * shi_425[k]
                   + f_3 * pc_y[k] * sii_593[k];

        t_765[k] = f_6 * sih0_450[k]
                   - f_7 * sih1_450[k]
                   + f_3 * pc_x[k] * sii_597[k];

        t_766[k] = f_8 * sih0_451[k]
                   - f_9 * sih1_451[k]
                   + f_3 * pc_x[k] * sii_598[k];

        t_767[k] = f_3 * pc_z[k] * sii_594[k];
    }

#pragma omp simd aligned(t_768, t_769, t_770, pc_x, pc_y, shi_429, sih0_453, sih0_455, \
                         sih1_453, sih1_455, sii_597, sii_600, \
                         sii_602 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_768[k] = f_8 * sih0_453[k]
                   - f_9 * sih1_453[k]
                   + f_3 * pc_x[k] * sii_600[k];

        t_769[k] = f_0 * shi_429[k]
                   + f_3 * pc_y[k] * sii_597[k];

        t_770[k] = f_8 * sih0_455[k]
                   - f_9 * sih1_455[k]
                   + f_3 * pc_x[k] * sii_602[k];
    }

#pragma omp simd aligned(t_771, t_772, t_773, t_774, pc_x, pc_z, sih0_456, sih0_458, sih0_459, \
                         sih1_456, sih1_458, sih1_459, sii_598, sii_603, sii_605, \
                         sii_606 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_771[k] = f_10 * sih0_456[k]
                   - f_11 * sih1_456[k]
                   + f_3 * pc_x[k] * sii_603[k];

        t_772[k] = f_3 * pc_z[k] * sii_598[k];

        t_773[k] = f_10 * sih0_458[k]
                   - f_11 * sih1_458[k]
                   + f_3 * pc_x[k] * sii_605[k];

        t_774[k] = f_10 * sih0_459[k]
                   - f_11 * sih1_459[k]
                   + f_3 * pc_x[k] * sii_606[k];
    }

#pragma omp simd aligned(t_775, t_776, t_777, t_778, t_779, pc_x, pc_y, shi_434, sih0_461, \
                         sih1_461, sii_602, sii_608, sii_609, sii_610, \
                         sii_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_775[k] = f_0 * shi_434[k]
                   + f_3 * pc_y[k] * sii_602[k];

        t_776[k] = f_10 * sih0_461[k]
                   - f_11 * sih1_461[k]
                   + f_3 * pc_x[k] * sii_608[k];

        t_777[k] = f_3 * pc_x[k] * sii_609[k];

        t_778[k] = f_3 * pc_x[k] * sii_610[k];

        t_779[k] = f_3 * pc_x[k] * sii_611[k];
    }

#pragma omp simd aligned(t_780, t_781, t_782, t_783, t_784, pc_x, pc_y, shi_441, sih0_456, \
                         sih1_456, sii_609, sii_612, sii_613, sii_614, \
                         sii_615 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_780[k] = f_3 * pc_x[k] * sii_612[k];

        t_781[k] = f_3 * pc_x[k] * sii_613[k];

        t_782[k] = f_3 * pc_x[k] * sii_614[k];

        t_783[k] = f_3 * pc_x[k] * sii_615[k];

        t_784[k] = f_0 * shi_441[k]
                   + f_1 * sih0_456[k]
                   - f_2 * sih1_456[k]
                   + f_3 * pc_y[k] * sii_609[k];
    }

#pragma omp simd aligned(t_785, t_786, t_787, pc_y, pc_z, shi_443, shi_444, sih0_458, \
                         sih0_459, sih1_458, sih1_459, sii_609, sii_611, \
                         sii_612 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_785[k] = f_3 * pc_z[k] * sii_609[k];

        t_786[k] = f_0 * shi_443[k]
                   + f_4 * sih0_458[k]
                   - f_5 * sih1_458[k]
                   + f_3 * pc_y[k] * sii_611[k];

        t_787[k] = f_0 * shi_444[k]
                   + f_6 * sih0_459[k]
                   - f_7 * sih1_459[k]
                   + f_3 * pc_y[k] * sii_612[k];
    }

#pragma omp simd aligned(t_788, t_789, t_790, t_791, pc_y, pc_z, shi_445, shi_446, shi_447, \
                         sih0_460, sih0_461, sih1_460, sih1_461, sii_613, sii_614, \
                         sii_615 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_788[k] = f_0 * shi_445[k]
                   + f_8 * sih0_460[k]
                   - f_9 * sih1_460[k]
                   + f_3 * pc_y[k] * sii_613[k];

        t_789[k] = f_0 * shi_446[k]
                   + f_10 * sih0_461[k]
                   - f_11 * sih1_461[k]
                   + f_3 * pc_y[k] * sii_614[k];

        t_790[k] = f_0 * shi_447[k]
                   + f_3 * pc_y[k] * sii_615[k];

        t_791[k] = f_1 * sih0_461[k]
                   - f_2 * sih1_461[k]
                   + f_3 * pc_z[k] * sii_615[k];
    }

#pragma omp simd aligned(t_792, t_793, t_794, t_795, pb_z, pc_y, pc_z, shk0_540, shk0_543, \
                         shi_420, shi_448, shk1_540, shk1_543, \
                         sii_616 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_792[k] = pb_z[k] * shk0_540[k]
                   - f_12 * pc_z[k] * shk1_540[k];

        t_793[k] = f_17 * shi_448[k]
                   + f_3 * pc_y[k] * sii_616[k];

        t_794[k] = f_13 * shi_420[k]
                   + f_3 * pc_z[k] * sii_616[k];

        t_795[k] = pb_z[k] * shk0_543[k]
                   - f_12 * pc_z[k] * shk1_543[k];
    }

#pragma omp simd aligned(t_796, t_797, t_798, pb_z, pc_x, pc_y, pc_z, shk0_546, shi_450, \
                         shk1_546, sih0_467, sih1_467, sii_618, \
                         sii_621 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_796[k] = f_17 * shi_450[k]
                   + f_3 * pc_y[k] * sii_618[k];

        t_797[k] = f_4 * sih0_467[k]
                   - f_5 * sih1_467[k]
                   + f_3 * pc_x[k] * sii_621[k];

        t_798[k] = pb_z[k] * shk0_546[k]
                   - f_12 * pc_z[k] * shk1_546[k];
    }

#pragma omp simd aligned(t_799, t_800, t_801, pc_x, pc_y, pc_z, shi_423, shi_453, sih0_471, \
                         sih1_471, sii_619, sii_621, sii_625 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_799[k] = f_13 * shi_423[k]
                   + f_3 * pc_z[k] * sii_619[k];

        t_800[k] = f_17 * shi_453[k]
                   + f_3 * pc_y[k] * sii_621[k];

        t_801[k] = f_6 * sih0_471[k]
                   - f_7 * sih1_471[k]
                   + f_3 * pc_x[k] * sii_625[k];
    }

#pragma omp simd aligned(t_802, t_803, t_804, pb_z, pc_x, pc_z, shk0_550, shi_426, shk1_550, \
                         sih0_474, sih1_474, sii_622, sii_628 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_802[k] = pb_z[k] * shk0_550[k]
                   - f_12 * pc_z[k] * shk1_550[k];

        t_803[k] = f_13 * shi_426[k]
                   + f_3 * pc_z[k] * sii_622[k];

        t_804[k] = f_8 * sih0_474[k]
                   - f_9 * sih1_474[k]
                   + f_3 * pc_x[k] * sii_628[k];
    }

#pragma omp simd aligned(t_805, t_806, t_807, pb_z, pc_x, pc_y, pc_z, shk0_555, shi_457, \
                         shk1_555, sih0_476, sih1_476, sii_625, \
                         sii_630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_805[k] = f_17 * shi_457[k]
                   + f_3 * pc_y[k] * sii_625[k];

        t_806[k] = f_8 * sih0_476[k]
                   - f_9 * sih1_476[k]
                   + f_3 * pc_x[k] * sii_630[k];

        t_807[k] = pb_z[k] * shk0_555[k]
                   - f_12 * pc_z[k] * shk1_555[k];
    }

#pragma omp simd aligned(t_808, t_809, t_810, pc_x, pc_z, shi_430, sih0_479, sih0_480, \
                         sih1_479, sih1_480, sii_626, sii_633, \
                         sii_634 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_808[k] = f_13 * shi_430[k]
                   + f_3 * pc_z[k] * sii_626[k];

        t_809[k] = f_10 * sih0_479[k]
                   - f_11 * sih1_479[k]
                   + f_3 * pc_x[k] * sii_633[k];

        t_810[k] = f_10 * sih0_480[k]
                   - f_11 * sih1_480[k]
                   + f_3 * pc_x[k] * sii_634[k];
    }
}

static auto
compute_prim_sik_three_center_electron_repulsion_0_piece7(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t shk0,
                                                          const size_t shi, const size_t shk1,
                                                          const size_t sih0, const size_t sih1,
                                                          const size_t sii, const size_t ncols,
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

    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *shk0_568 = buffer.data(shk0 + 568);
    const auto *shk0_570 = buffer.data(shk0 + 570);
    const auto *shk0_571 = buffer.data(shk0 + 571);
    const auto *shk0_572 = buffer.data(shk0 + 572);
    const auto *shk0_573 = buffer.data(shk0 + 573);

    const auto *shi_441 = buffer.data(shi + 441);
    const auto *shi_442 = buffer.data(shi + 442);
    const auto *shi_443 = buffer.data(shi + 443);
    const auto *shi_444 = buffer.data(shi + 444);
    const auto *shi_445 = buffer.data(shi + 445);
    const auto *shi_447 = buffer.data(shi + 447);
    const auto *shi_448 = buffer.data(shi + 448);
    const auto *shi_451 = buffer.data(shi + 451);
    const auto *shi_454 = buffer.data(shi + 454);
    const auto *shi_458 = buffer.data(shi + 458);
    const auto *shi_462 = buffer.data(shi + 462);
    const auto *shi_469 = buffer.data(shi + 469);
    const auto *shi_475 = buffer.data(shi + 475);
    const auto *shi_476 = buffer.data(shi + 476);
    const auto *shi_478 = buffer.data(shi + 478);
    const auto *shi_479 = buffer.data(shi + 479);
    const auto *shi_481 = buffer.data(shi + 481);
    const auto *shi_482 = buffer.data(shi + 482);
    const auto *shi_485 = buffer.data(shi + 485);
    const auto *shi_486 = buffer.data(shi + 486);
    const auto *shi_490 = buffer.data(shi + 490);
    const auto *shi_497 = buffer.data(shi + 497);
    const auto *shi_499 = buffer.data(shi + 499);
    const auto *shi_500 = buffer.data(shi + 500);
    const auto *shi_501 = buffer.data(shi + 501);
    const auto *shi_502 = buffer.data(shi + 502);
    const auto *shi_503 = buffer.data(shi + 503);
    const auto *shi_504 = buffer.data(shi + 504);
    const auto *shi_506 = buffer.data(shi + 506);
    const auto *shi_507 = buffer.data(shi + 507);
    const auto *shi_509 = buffer.data(shi + 509);
    const auto *shi_510 = buffer.data(shi + 510);
    const auto *shi_513 = buffer.data(shi + 513);
    const auto *shi_514 = buffer.data(shi + 514);
    const auto *shi_518 = buffer.data(shi + 518);
    const auto *shi_525 = buffer.data(shi + 525);
    const auto *shi_527 = buffer.data(shi + 527);
    const auto *shi_528 = buffer.data(shi + 528);
    const auto *shi_529 = buffer.data(shi + 529);
    const auto *shi_530 = buffer.data(shi + 530);
    const auto *shi_531 = buffer.data(shi + 531);
    const auto *shi_532 = buffer.data(shi + 532);
    const auto *shi_534 = buffer.data(shi + 534);
    const auto *shi_537 = buffer.data(shi + 537);
    const auto *shi_541 = buffer.data(shi + 541);
    const auto *shi_546 = buffer.data(shi + 546);
    const auto *shi_553 = buffer.data(shi + 553);
    const auto *shi_555 = buffer.data(shi + 555);
    const auto *shi_556 = buffer.data(shi + 556);
    const auto *shi_557 = buffer.data(shi + 557);
    const auto *shi_558 = buffer.data(shi + 558);

    const auto *shk1_568 = buffer.data(shk1 + 568);
    const auto *shk1_570 = buffer.data(shk1 + 570);
    const auto *shk1_571 = buffer.data(shk1 + 571);
    const auto *shk1_572 = buffer.data(shk1 + 572);
    const auto *shk1_573 = buffer.data(shk1 + 573);

    const auto *sih0_482 = buffer.data(sih0 + 482);
    const auto *sih0_483 = buffer.data(sih0 + 483);
    const auto *sih0_486 = buffer.data(sih0 + 486);
    const auto *sih0_488 = buffer.data(sih0 + 488);
    const auto *sih0_489 = buffer.data(sih0 + 489);
    const auto *sih0_492 = buffer.data(sih0 + 492);
    const auto *sih0_493 = buffer.data(sih0 + 493);
    const auto *sih0_495 = buffer.data(sih0 + 495);
    const auto *sih0_497 = buffer.data(sih0 + 497);
    const auto *sih0_498 = buffer.data(sih0 + 498);
    const auto *sih0_500 = buffer.data(sih0 + 500);
    const auto *sih0_501 = buffer.data(sih0 + 501);
    const auto *sih0_502 = buffer.data(sih0 + 502);
    const auto *sih0_503 = buffer.data(sih0 + 503);
    const auto *sih0_504 = buffer.data(sih0 + 504);
    const auto *sih0_507 = buffer.data(sih0 + 507);
    const auto *sih0_509 = buffer.data(sih0 + 509);
    const auto *sih0_510 = buffer.data(sih0 + 510);
    const auto *sih0_513 = buffer.data(sih0 + 513);
    const auto *sih0_514 = buffer.data(sih0 + 514);
    const auto *sih0_516 = buffer.data(sih0 + 516);
    const auto *sih0_518 = buffer.data(sih0 + 518);
    const auto *sih0_519 = buffer.data(sih0 + 519);
    const auto *sih0_521 = buffer.data(sih0 + 521);
    const auto *sih0_522 = buffer.data(sih0 + 522);
    const auto *sih0_523 = buffer.data(sih0 + 523);
    const auto *sih0_524 = buffer.data(sih0 + 524);
    const auto *sih0_525 = buffer.data(sih0 + 525);
    const auto *sih0_528 = buffer.data(sih0 + 528);
    const auto *sih0_530 = buffer.data(sih0 + 530);
    const auto *sih0_531 = buffer.data(sih0 + 531);
    const auto *sih0_534 = buffer.data(sih0 + 534);
    const auto *sih0_535 = buffer.data(sih0 + 535);
    const auto *sih0_537 = buffer.data(sih0 + 537);
    const auto *sih0_539 = buffer.data(sih0 + 539);
    const auto *sih0_540 = buffer.data(sih0 + 540);
    const auto *sih0_542 = buffer.data(sih0 + 542);
    const auto *sih0_543 = buffer.data(sih0 + 543);
    const auto *sih0_544 = buffer.data(sih0 + 544);
    const auto *sih0_545 = buffer.data(sih0 + 545);

    const auto *sih1_482 = buffer.data(sih1 + 482);
    const auto *sih1_483 = buffer.data(sih1 + 483);
    const auto *sih1_486 = buffer.data(sih1 + 486);
    const auto *sih1_488 = buffer.data(sih1 + 488);
    const auto *sih1_489 = buffer.data(sih1 + 489);
    const auto *sih1_492 = buffer.data(sih1 + 492);
    const auto *sih1_493 = buffer.data(sih1 + 493);
    const auto *sih1_495 = buffer.data(sih1 + 495);
    const auto *sih1_497 = buffer.data(sih1 + 497);
    const auto *sih1_498 = buffer.data(sih1 + 498);
    const auto *sih1_500 = buffer.data(sih1 + 500);
    const auto *sih1_501 = buffer.data(sih1 + 501);
    const auto *sih1_502 = buffer.data(sih1 + 502);
    const auto *sih1_503 = buffer.data(sih1 + 503);
    const auto *sih1_504 = buffer.data(sih1 + 504);
    const auto *sih1_507 = buffer.data(sih1 + 507);
    const auto *sih1_509 = buffer.data(sih1 + 509);
    const auto *sih1_510 = buffer.data(sih1 + 510);
    const auto *sih1_513 = buffer.data(sih1 + 513);
    const auto *sih1_514 = buffer.data(sih1 + 514);
    const auto *sih1_516 = buffer.data(sih1 + 516);
    const auto *sih1_518 = buffer.data(sih1 + 518);
    const auto *sih1_519 = buffer.data(sih1 + 519);
    const auto *sih1_521 = buffer.data(sih1 + 521);
    const auto *sih1_522 = buffer.data(sih1 + 522);
    const auto *sih1_523 = buffer.data(sih1 + 523);
    const auto *sih1_524 = buffer.data(sih1 + 524);
    const auto *sih1_525 = buffer.data(sih1 + 525);
    const auto *sih1_528 = buffer.data(sih1 + 528);
    const auto *sih1_530 = buffer.data(sih1 + 530);
    const auto *sih1_531 = buffer.data(sih1 + 531);
    const auto *sih1_534 = buffer.data(sih1 + 534);
    const auto *sih1_535 = buffer.data(sih1 + 535);
    const auto *sih1_537 = buffer.data(sih1 + 537);
    const auto *sih1_539 = buffer.data(sih1 + 539);
    const auto *sih1_540 = buffer.data(sih1 + 540);
    const auto *sih1_542 = buffer.data(sih1 + 542);
    const auto *sih1_543 = buffer.data(sih1 + 543);
    const auto *sih1_544 = buffer.data(sih1 + 544);
    const auto *sih1_545 = buffer.data(sih1 + 545);

    const auto *sii_630 = buffer.data(sii + 630);
    const auto *sii_636 = buffer.data(sii + 636);
    const auto *sii_637 = buffer.data(sii + 637);
    const auto *sii_638 = buffer.data(sii + 638);
    const auto *sii_639 = buffer.data(sii + 639);
    const auto *sii_640 = buffer.data(sii + 640);
    const auto *sii_641 = buffer.data(sii + 641);
    const auto *sii_642 = buffer.data(sii + 642);
    const auto *sii_643 = buffer.data(sii + 643);
    const auto *sii_644 = buffer.data(sii + 644);
    const auto *sii_646 = buffer.data(sii + 646);
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
    const auto *sii_674 = buffer.data(sii + 674);
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
    const auto *sii_702 = buffer.data(sii + 702);
    const auto *sii_703 = buffer.data(sii + 703);
    const auto *sii_705 = buffer.data(sii + 705);
    const auto *sii_706 = buffer.data(sii + 706);
    const auto *sii_709 = buffer.data(sii + 709);
    const auto *sii_710 = buffer.data(sii + 710);
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

#pragma omp simd aligned(t_811, t_812, t_813, t_814, t_815, pc_x, pc_y, shi_462, sih0_482, \
                         sih1_482, sii_630, sii_636, sii_637, sii_638, \
                         sii_639 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_811[k] = f_17 * shi_462[k]
                   + f_3 * pc_y[k] * sii_630[k];

        t_812[k] = f_10 * sih0_482[k]
                   - f_11 * sih1_482[k]
                   + f_3 * pc_x[k] * sii_636[k];

        t_813[k] = f_3 * pc_x[k] * sii_637[k];

        t_814[k] = f_3 * pc_x[k] * sii_638[k];

        t_815[k] = f_3 * pc_x[k] * sii_639[k];
    }

#pragma omp simd aligned(t_816, t_817, t_818, t_819, t_820, pb_z, pc_x, pc_z, shk0_568, \
                         shk1_568, sii_640, sii_641, sii_642, sii_643 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_816[k] = f_3 * pc_x[k] * sii_640[k];

        t_817[k] = f_3 * pc_x[k] * sii_641[k];

        t_818[k] = f_3 * pc_x[k] * sii_642[k];

        t_819[k] = f_3 * pc_x[k] * sii_643[k];

        t_820[k] = pb_z[k] * shk0_568[k]
                   - f_12 * pc_z[k] * shk1_568[k];
    }

#pragma omp simd aligned(t_821, t_822, t_823, pb_z, pc_z, shk0_570, shk0_571, shi_441, \
                         shi_442, shi_443, shk1_570, shk1_571, \
                         sii_637 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_821[k] = f_13 * shi_441[k]
                   + f_3 * pc_z[k] * sii_637[k];

        t_822[k] = pb_z[k] * shk0_570[k]
                   + f_14 * shi_442[k]
                   - f_12 * pc_z[k] * shk1_570[k];

        t_823[k] = pb_z[k] * shk0_571[k]
                   + f_15 * shi_443[k]
                   - f_12 * pc_z[k] * shk1_571[k];
    }

#pragma omp simd aligned(t_824, t_825, t_826, pb_z, pc_y, pc_z, shk0_572, shk0_573, shi_444, \
                         shi_445, shi_475, shk1_572, shk1_573, \
                         sii_643 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_824[k] = pb_z[k] * shk0_572[k]
                   + f_16 * shi_444[k]
                   - f_12 * pc_z[k] * shk1_572[k];

        t_825[k] = pb_z[k] * shk0_573[k]
                   + f_17 * shi_445[k]
                   - f_12 * pc_z[k] * shk1_573[k];

        t_826[k] = f_17 * shi_475[k]
                   + f_3 * pc_y[k] * sii_643[k];
    }

#pragma omp simd aligned(t_827, t_828, t_829, t_830, pc_x, pc_y, pc_z, shi_447, shi_448, \
                         shi_476, sih0_482, sih0_483, sih1_482, sih1_483, sii_643, \
                         sii_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_827[k] = f_13 * shi_447[k]
                   + f_1 * sih0_482[k]
                   - f_2 * sih1_482[k]
                   + f_3 * pc_z[k] * sii_643[k];

        t_828[k] = f_1 * sih0_483[k]
                   - f_2 * sih1_483[k]
                   + f_3 * pc_x[k] * sii_644[k];

        t_829[k] = f_16 * shi_476[k]
                   + f_3 * pc_y[k] * sii_644[k];

        t_830[k] = f_14 * shi_448[k]
                   + f_3 * pc_z[k] * sii_644[k];
    }

#pragma omp simd aligned(t_831, t_832, t_833, pc_x, pc_y, shi_478, sih0_486, sih0_488, \
                         sih1_486, sih1_488, sii_646, sii_647, \
                         sii_649 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_831[k] = f_4 * sih0_486[k]
                   - f_5 * sih1_486[k]
                   + f_3 * pc_x[k] * sii_647[k];

        t_832[k] = f_16 * shi_478[k]
                   + f_3 * pc_y[k] * sii_646[k];

        t_833[k] = f_4 * sih0_488[k]
                   - f_5 * sih1_488[k]
                   + f_3 * pc_x[k] * sii_649[k];
    }

#pragma omp simd aligned(t_834, t_835, t_836, pc_x, pc_y, pc_z, shi_451, shi_481, sih0_489, \
                         sih1_489, sii_647, sii_649, sii_650 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_834[k] = f_6 * sih0_489[k]
                   - f_7 * sih1_489[k]
                   + f_3 * pc_x[k] * sii_650[k];

        t_835[k] = f_14 * shi_451[k]
                   + f_3 * pc_z[k] * sii_647[k];

        t_836[k] = f_16 * shi_481[k]
                   + f_3 * pc_y[k] * sii_649[k];
    }

#pragma omp simd aligned(t_837, t_838, t_839, pc_x, pc_z, shi_454, sih0_492, sih0_493, \
                         sih1_492, sih1_493, sii_650, sii_653, \
                         sii_654 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_837[k] = f_6 * sih0_492[k]
                   - f_7 * sih1_492[k]
                   + f_3 * pc_x[k] * sii_653[k];

        t_838[k] = f_8 * sih0_493[k]
                   - f_9 * sih1_493[k]
                   + f_3 * pc_x[k] * sii_654[k];

        t_839[k] = f_14 * shi_454[k]
                   + f_3 * pc_z[k] * sii_650[k];
    }

#pragma omp simd aligned(t_840, t_841, t_842, pc_x, pc_y, shi_485, sih0_495, sih0_497, \
                         sih1_495, sih1_497, sii_653, sii_656, \
                         sii_658 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_840[k] = f_8 * sih0_495[k]
                   - f_9 * sih1_495[k]
                   + f_3 * pc_x[k] * sii_656[k];

        t_841[k] = f_16 * shi_485[k]
                   + f_3 * pc_y[k] * sii_653[k];

        t_842[k] = f_8 * sih0_497[k]
                   - f_9 * sih1_497[k]
                   + f_3 * pc_x[k] * sii_658[k];
    }

#pragma omp simd aligned(t_843, t_844, t_845, pc_x, pc_z, shi_458, sih0_498, sih0_500, \
                         sih1_498, sih1_500, sii_654, sii_659, \
                         sii_661 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_843[k] = f_10 * sih0_498[k]
                   - f_11 * sih1_498[k]
                   + f_3 * pc_x[k] * sii_659[k];

        t_844[k] = f_14 * shi_458[k]
                   + f_3 * pc_z[k] * sii_654[k];

        t_845[k] = f_10 * sih0_500[k]
                   - f_11 * sih1_500[k]
                   + f_3 * pc_x[k] * sii_661[k];
    }

#pragma omp simd aligned(t_846, t_847, t_848, t_849, pc_x, pc_y, shi_490, sih0_501, sih0_503, \
                         sih1_501, sih1_503, sii_658, sii_662, sii_664, \
                         sii_665 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_846[k] = f_10 * sih0_501[k]
                   - f_11 * sih1_501[k]
                   + f_3 * pc_x[k] * sii_662[k];

        t_847[k] = f_16 * shi_490[k]
                   + f_3 * pc_y[k] * sii_658[k];

        t_848[k] = f_10 * sih0_503[k]
                   - f_11 * sih1_503[k]
                   + f_3 * pc_x[k] * sii_664[k];

        t_849[k] = f_3 * pc_x[k] * sii_665[k];
    }

#pragma omp simd aligned(t_850, t_851, t_852, t_853, t_854, t_855, pc_x, sii_666, sii_667, \
                         sii_668, sii_669, sii_670, sii_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_850[k] = f_3 * pc_x[k] * sii_666[k];

        t_851[k] = f_3 * pc_x[k] * sii_667[k];

        t_852[k] = f_3 * pc_x[k] * sii_668[k];

        t_853[k] = f_3 * pc_x[k] * sii_669[k];

        t_854[k] = f_3 * pc_x[k] * sii_670[k];

        t_855[k] = f_3 * pc_x[k] * sii_671[k];
    }

#pragma omp simd aligned(t_856, t_857, t_858, pc_y, pc_z, shi_469, shi_497, shi_499, sih0_498, \
                         sih0_500, sih1_498, sih1_500, sii_665, \
                         sii_667 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_856[k] = f_16 * shi_497[k]
                   + f_1 * sih0_498[k]
                   - f_2 * sih1_498[k]
                   + f_3 * pc_y[k] * sii_665[k];

        t_857[k] = f_14 * shi_469[k]
                   + f_3 * pc_z[k] * sii_665[k];

        t_858[k] = f_16 * shi_499[k]
                   + f_4 * sih0_500[k]
                   - f_5 * sih1_500[k]
                   + f_3 * pc_y[k] * sii_667[k];
    }

#pragma omp simd aligned(t_859, t_860, t_861, pc_y, shi_500, shi_501, shi_502, sih0_501, \
                         sih0_502, sih0_503, sih1_501, sih1_502, sih1_503, sii_668, sii_669, \
                         sii_670 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_859[k] = f_16 * shi_500[k]
                   + f_6 * sih0_501[k]
                   - f_7 * sih1_501[k]
                   + f_3 * pc_y[k] * sii_668[k];

        t_860[k] = f_16 * shi_501[k]
                   + f_8 * sih0_502[k]
                   - f_9 * sih1_502[k]
                   + f_3 * pc_y[k] * sii_669[k];

        t_861[k] = f_16 * shi_502[k]
                   + f_10 * sih0_503[k]
                   - f_11 * sih1_503[k]
                   + f_3 * pc_y[k] * sii_670[k];
    }

#pragma omp simd aligned(t_862, t_863, t_864, t_865, pc_x, pc_y, pc_z, shi_475, shi_503, \
                         shi_504, sih0_503, sih0_504, sih1_503, sih1_504, sii_671, \
                         sii_672 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_862[k] = f_16 * shi_503[k]
                   + f_3 * pc_y[k] * sii_671[k];

        t_863[k] = f_14 * shi_475[k]
                   + f_1 * sih0_503[k]
                   - f_2 * sih1_503[k]
                   + f_3 * pc_z[k] * sii_671[k];

        t_864[k] = f_1 * sih0_504[k]
                   - f_2 * sih1_504[k]
                   + f_3 * pc_x[k] * sii_672[k];

        t_865[k] = f_15 * shi_504[k]
                   + f_3 * pc_y[k] * sii_672[k];
    }

#pragma omp simd aligned(t_866, t_867, t_868, pc_x, pc_y, pc_z, shi_476, shi_506, sih0_507, \
                         sih1_507, sii_672, sii_674, sii_675 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_866[k] = f_15 * shi_476[k]
                   + f_3 * pc_z[k] * sii_672[k];

        t_867[k] = f_4 * sih0_507[k]
                   - f_5 * sih1_507[k]
                   + f_3 * pc_x[k] * sii_675[k];

        t_868[k] = f_15 * shi_506[k]
                   + f_3 * pc_y[k] * sii_674[k];
    }

#pragma omp simd aligned(t_869, t_870, t_871, t_872, pc_x, pc_y, pc_z, shi_479, shi_509, \
                         sih0_509, sih0_510, sih1_509, sih1_510, sii_675, sii_677, \
                         sii_678 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_869[k] = f_4 * sih0_509[k]
                   - f_5 * sih1_509[k]
                   + f_3 * pc_x[k] * sii_677[k];

        t_870[k] = f_6 * sih0_510[k]
                   - f_7 * sih1_510[k]
                   + f_3 * pc_x[k] * sii_678[k];

        t_871[k] = f_15 * shi_479[k]
                   + f_3 * pc_z[k] * sii_675[k];

        t_872[k] = f_15 * shi_509[k]
                   + f_3 * pc_y[k] * sii_677[k];
    }

#pragma omp simd aligned(t_873, t_874, t_875, pc_x, pc_z, shi_482, sih0_513, sih0_514, \
                         sih1_513, sih1_514, sii_678, sii_681, \
                         sii_682 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_873[k] = f_6 * sih0_513[k]
                   - f_7 * sih1_513[k]
                   + f_3 * pc_x[k] * sii_681[k];

        t_874[k] = f_8 * sih0_514[k]
                   - f_9 * sih1_514[k]
                   + f_3 * pc_x[k] * sii_682[k];

        t_875[k] = f_15 * shi_482[k]
                   + f_3 * pc_z[k] * sii_678[k];
    }

#pragma omp simd aligned(t_876, t_877, t_878, pc_x, pc_y, shi_513, sih0_516, sih0_518, \
                         sih1_516, sih1_518, sii_681, sii_684, \
                         sii_686 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_876[k] = f_8 * sih0_516[k]
                   - f_9 * sih1_516[k]
                   + f_3 * pc_x[k] * sii_684[k];

        t_877[k] = f_15 * shi_513[k]
                   + f_3 * pc_y[k] * sii_681[k];

        t_878[k] = f_8 * sih0_518[k]
                   - f_9 * sih1_518[k]
                   + f_3 * pc_x[k] * sii_686[k];
    }

#pragma omp simd aligned(t_879, t_880, t_881, pc_x, pc_z, shi_486, sih0_519, sih0_521, \
                         sih1_519, sih1_521, sii_682, sii_687, \
                         sii_689 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_879[k] = f_10 * sih0_519[k]
                   - f_11 * sih1_519[k]
                   + f_3 * pc_x[k] * sii_687[k];

        t_880[k] = f_15 * shi_486[k]
                   + f_3 * pc_z[k] * sii_682[k];

        t_881[k] = f_10 * sih0_521[k]
                   - f_11 * sih1_521[k]
                   + f_3 * pc_x[k] * sii_689[k];
    }

#pragma omp simd aligned(t_882, t_883, t_884, t_885, pc_x, pc_y, shi_518, sih0_522, sih0_524, \
                         sih1_522, sih1_524, sii_686, sii_690, sii_692, \
                         sii_693 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_882[k] = f_10 * sih0_522[k]
                   - f_11 * sih1_522[k]
                   + f_3 * pc_x[k] * sii_690[k];

        t_883[k] = f_15 * shi_518[k]
                   + f_3 * pc_y[k] * sii_686[k];

        t_884[k] = f_10 * sih0_524[k]
                   - f_11 * sih1_524[k]
                   + f_3 * pc_x[k] * sii_692[k];

        t_885[k] = f_3 * pc_x[k] * sii_693[k];
    }

#pragma omp simd aligned(t_886, t_887, t_888, t_889, t_890, t_891, pc_x, sii_694, sii_695, \
                         sii_696, sii_697, sii_698, sii_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_886[k] = f_3 * pc_x[k] * sii_694[k];

        t_887[k] = f_3 * pc_x[k] * sii_695[k];

        t_888[k] = f_3 * pc_x[k] * sii_696[k];

        t_889[k] = f_3 * pc_x[k] * sii_697[k];

        t_890[k] = f_3 * pc_x[k] * sii_698[k];

        t_891[k] = f_3 * pc_x[k] * sii_699[k];
    }

#pragma omp simd aligned(t_892, t_893, t_894, pc_y, pc_z, shi_497, shi_525, shi_527, sih0_519, \
                         sih0_521, sih1_519, sih1_521, sii_693, \
                         sii_695 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_892[k] = f_15 * shi_525[k]
                   + f_1 * sih0_519[k]
                   - f_2 * sih1_519[k]
                   + f_3 * pc_y[k] * sii_693[k];

        t_893[k] = f_15 * shi_497[k]
                   + f_3 * pc_z[k] * sii_693[k];

        t_894[k] = f_15 * shi_527[k]
                   + f_4 * sih0_521[k]
                   - f_5 * sih1_521[k]
                   + f_3 * pc_y[k] * sii_695[k];
    }

#pragma omp simd aligned(t_895, t_896, t_897, pc_y, shi_528, shi_529, shi_530, sih0_522, \
                         sih0_523, sih0_524, sih1_522, sih1_523, sih1_524, sii_696, sii_697, \
                         sii_698 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_895[k] = f_15 * shi_528[k]
                   + f_6 * sih0_522[k]
                   - f_7 * sih1_522[k]
                   + f_3 * pc_y[k] * sii_696[k];

        t_896[k] = f_15 * shi_529[k]
                   + f_8 * sih0_523[k]
                   - f_9 * sih1_523[k]
                   + f_3 * pc_y[k] * sii_697[k];

        t_897[k] = f_15 * shi_530[k]
                   + f_10 * sih0_524[k]
                   - f_11 * sih1_524[k]
                   + f_3 * pc_y[k] * sii_698[k];
    }

#pragma omp simd aligned(t_898, t_899, t_900, t_901, pc_x, pc_y, pc_z, shi_503, shi_531, \
                         shi_532, sih0_524, sih0_525, sih1_524, sih1_525, sii_699, \
                         sii_700 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_898[k] = f_15 * shi_531[k]
                   + f_3 * pc_y[k] * sii_699[k];

        t_899[k] = f_15 * shi_503[k]
                   + f_1 * sih0_524[k]
                   - f_2 * sih1_524[k]
                   + f_3 * pc_z[k] * sii_699[k];

        t_900[k] = f_1 * sih0_525[k]
                   - f_2 * sih1_525[k]
                   + f_3 * pc_x[k] * sii_700[k];

        t_901[k] = f_14 * shi_532[k]
                   + f_3 * pc_y[k] * sii_700[k];
    }

#pragma omp simd aligned(t_902, t_903, t_904, pc_x, pc_y, pc_z, shi_504, shi_534, sih0_528, \
                         sih1_528, sii_700, sii_702, sii_703 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_902[k] = f_16 * shi_504[k]
                   + f_3 * pc_z[k] * sii_700[k];

        t_903[k] = f_4 * sih0_528[k]
                   - f_5 * sih1_528[k]
                   + f_3 * pc_x[k] * sii_703[k];

        t_904[k] = f_14 * shi_534[k]
                   + f_3 * pc_y[k] * sii_702[k];
    }

#pragma omp simd aligned(t_905, t_906, t_907, t_908, pc_x, pc_y, pc_z, shi_507, shi_537, \
                         sih0_530, sih0_531, sih1_530, sih1_531, sii_703, sii_705, \
                         sii_706 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_905[k] = f_4 * sih0_530[k]
                   - f_5 * sih1_530[k]
                   + f_3 * pc_x[k] * sii_705[k];

        t_906[k] = f_6 * sih0_531[k]
                   - f_7 * sih1_531[k]
                   + f_3 * pc_x[k] * sii_706[k];

        t_907[k] = f_16 * shi_507[k]
                   + f_3 * pc_z[k] * sii_703[k];

        t_908[k] = f_14 * shi_537[k]
                   + f_3 * pc_y[k] * sii_705[k];
    }

#pragma omp simd aligned(t_909, t_910, t_911, pc_x, pc_z, shi_510, sih0_534, sih0_535, \
                         sih1_534, sih1_535, sii_706, sii_709, \
                         sii_710 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_909[k] = f_6 * sih0_534[k]
                   - f_7 * sih1_534[k]
                   + f_3 * pc_x[k] * sii_709[k];

        t_910[k] = f_8 * sih0_535[k]
                   - f_9 * sih1_535[k]
                   + f_3 * pc_x[k] * sii_710[k];

        t_911[k] = f_16 * shi_510[k]
                   + f_3 * pc_z[k] * sii_706[k];
    }

#pragma omp simd aligned(t_912, t_913, t_914, pc_x, pc_y, shi_541, sih0_537, sih0_539, \
                         sih1_537, sih1_539, sii_709, sii_712, \
                         sii_714 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_912[k] = f_8 * sih0_537[k]
                   - f_9 * sih1_537[k]
                   + f_3 * pc_x[k] * sii_712[k];

        t_913[k] = f_14 * shi_541[k]
                   + f_3 * pc_y[k] * sii_709[k];

        t_914[k] = f_8 * sih0_539[k]
                   - f_9 * sih1_539[k]
                   + f_3 * pc_x[k] * sii_714[k];
    }

#pragma omp simd aligned(t_915, t_916, t_917, pc_x, pc_z, shi_514, sih0_540, sih0_542, \
                         sih1_540, sih1_542, sii_710, sii_715, \
                         sii_717 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_915[k] = f_10 * sih0_540[k]
                   - f_11 * sih1_540[k]
                   + f_3 * pc_x[k] * sii_715[k];

        t_916[k] = f_16 * shi_514[k]
                   + f_3 * pc_z[k] * sii_710[k];

        t_917[k] = f_10 * sih0_542[k]
                   - f_11 * sih1_542[k]
                   + f_3 * pc_x[k] * sii_717[k];
    }

#pragma omp simd aligned(t_918, t_919, t_920, t_921, pc_x, pc_y, shi_546, sih0_543, sih0_545, \
                         sih1_543, sih1_545, sii_714, sii_718, sii_720, \
                         sii_721 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_918[k] = f_10 * sih0_543[k]
                   - f_11 * sih1_543[k]
                   + f_3 * pc_x[k] * sii_718[k];

        t_919[k] = f_14 * shi_546[k]
                   + f_3 * pc_y[k] * sii_714[k];

        t_920[k] = f_10 * sih0_545[k]
                   - f_11 * sih1_545[k]
                   + f_3 * pc_x[k] * sii_720[k];

        t_921[k] = f_3 * pc_x[k] * sii_721[k];
    }

#pragma omp simd aligned(t_922, t_923, t_924, t_925, t_926, t_927, pc_x, sii_722, sii_723, \
                         sii_724, sii_725, sii_726, sii_727 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_922[k] = f_3 * pc_x[k] * sii_722[k];

        t_923[k] = f_3 * pc_x[k] * sii_723[k];

        t_924[k] = f_3 * pc_x[k] * sii_724[k];

        t_925[k] = f_3 * pc_x[k] * sii_725[k];

        t_926[k] = f_3 * pc_x[k] * sii_726[k];

        t_927[k] = f_3 * pc_x[k] * sii_727[k];
    }

#pragma omp simd aligned(t_928, t_929, t_930, pc_y, pc_z, shi_525, shi_553, shi_555, sih0_540, \
                         sih0_542, sih1_540, sih1_542, sii_721, \
                         sii_723 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_928[k] = f_14 * shi_553[k]
                   + f_1 * sih0_540[k]
                   - f_2 * sih1_540[k]
                   + f_3 * pc_y[k] * sii_721[k];

        t_929[k] = f_16 * shi_525[k]
                   + f_3 * pc_z[k] * sii_721[k];

        t_930[k] = f_14 * shi_555[k]
                   + f_4 * sih0_542[k]
                   - f_5 * sih1_542[k]
                   + f_3 * pc_y[k] * sii_723[k];
    }

#pragma omp simd aligned(t_931, t_932, t_933, pc_y, shi_556, shi_557, shi_558, sih0_543, \
                         sih0_544, sih0_545, sih1_543, sih1_544, sih1_545, sii_724, sii_725, \
                         sii_726 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_931[k] = f_14 * shi_556[k]
                   + f_6 * sih0_543[k]
                   - f_7 * sih1_543[k]
                   + f_3 * pc_y[k] * sii_724[k];

        t_932[k] = f_14 * shi_557[k]
                   + f_8 * sih0_544[k]
                   - f_9 * sih1_544[k]
                   + f_3 * pc_y[k] * sii_725[k];

        t_933[k] = f_14 * shi_558[k]
                   + f_10 * sih0_545[k]
                   - f_11 * sih1_545[k]
                   + f_3 * pc_y[k] * sii_726[k];
    }
}

static auto
compute_prim_sik_three_center_electron_repulsion_0_piece8(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t shk0,
                                                          const size_t shi, const size_t shk1,
                                                          const size_t sih0, const size_t sih1,
                                                          const size_t sii, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
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
    const auto f_18 = 3.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *shk0_720 = buffer.data(shk0 + 720);
    const auto *shk0_725 = buffer.data(shk0 + 725);
    const auto *shk0_729 = buffer.data(shk0 + 729);
    const auto *shk0_734 = buffer.data(shk0 + 734);
    const auto *shk0_740 = buffer.data(shk0 + 740);
    const auto *shk0_748 = buffer.data(shk0 + 748);
    const auto *shk0_750 = buffer.data(shk0 + 750);
    const auto *shk0_751 = buffer.data(shk0 + 751);
    const auto *shk0_752 = buffer.data(shk0 + 752);
    const auto *shk0_753 = buffer.data(shk0 + 753);
    const auto *shk0_755 = buffer.data(shk0 + 755);

    const auto *shi_531 = buffer.data(shi + 531);
    const auto *shi_532 = buffer.data(shi + 532);
    const auto *shi_535 = buffer.data(shi + 535);
    const auto *shi_538 = buffer.data(shi + 538);
    const auto *shi_542 = buffer.data(shi + 542);
    const auto *shi_553 = buffer.data(shi + 553);
    const auto *shi_559 = buffer.data(shi + 559);
    const auto *shi_560 = buffer.data(shi + 560);
    const auto *shi_562 = buffer.data(shi + 562);
    const auto *shi_563 = buffer.data(shi + 563);
    const auto *shi_565 = buffer.data(shi + 565);
    const auto *shi_566 = buffer.data(shi + 566);
    const auto *shi_569 = buffer.data(shi + 569);
    const auto *shi_570 = buffer.data(shi + 570);
    const auto *shi_574 = buffer.data(shi + 574);
    const auto *shi_581 = buffer.data(shi + 581);
    const auto *shi_583 = buffer.data(shi + 583);
    const auto *shi_584 = buffer.data(shi + 584);
    const auto *shi_585 = buffer.data(shi + 585);
    const auto *shi_586 = buffer.data(shi + 586);
    const auto *shi_587 = buffer.data(shi + 587);

    const auto *shk1_720 = buffer.data(shk1 + 720);
    const auto *shk1_725 = buffer.data(shk1 + 725);
    const auto *shk1_729 = buffer.data(shk1 + 729);
    const auto *shk1_734 = buffer.data(shk1 + 734);
    const auto *shk1_740 = buffer.data(shk1 + 740);
    const auto *shk1_748 = buffer.data(shk1 + 748);
    const auto *shk1_750 = buffer.data(shk1 + 750);
    const auto *shk1_751 = buffer.data(shk1 + 751);
    const auto *shk1_752 = buffer.data(shk1 + 752);
    const auto *shk1_753 = buffer.data(shk1 + 753);
    const auto *shk1_755 = buffer.data(shk1 + 755);

    const auto *sih0_545 = buffer.data(sih0 + 545);
    const auto *sih0_549 = buffer.data(sih0 + 549);
    const auto *sih0_552 = buffer.data(sih0 + 552);
    const auto *sih0_556 = buffer.data(sih0 + 556);
    const auto *sih0_558 = buffer.data(sih0 + 558);
    const auto *sih0_561 = buffer.data(sih0 + 561);
    const auto *sih0_563 = buffer.data(sih0 + 563);
    const auto *sih0_564 = buffer.data(sih0 + 564);
    const auto *sih0_567 = buffer.data(sih0 + 567);
    const auto *sih0_570 = buffer.data(sih0 + 570);
    const auto *sih0_572 = buffer.data(sih0 + 572);
    const auto *sih0_573 = buffer.data(sih0 + 573);
    const auto *sih0_576 = buffer.data(sih0 + 576);
    const auto *sih0_577 = buffer.data(sih0 + 577);
    const auto *sih0_579 = buffer.data(sih0 + 579);
    const auto *sih0_581 = buffer.data(sih0 + 581);
    const auto *sih0_582 = buffer.data(sih0 + 582);
    const auto *sih0_584 = buffer.data(sih0 + 584);
    const auto *sih0_585 = buffer.data(sih0 + 585);
    const auto *sih0_586 = buffer.data(sih0 + 586);
    const auto *sih0_587 = buffer.data(sih0 + 587);

    const auto *sih1_545 = buffer.data(sih1 + 545);
    const auto *sih1_549 = buffer.data(sih1 + 549);
    const auto *sih1_552 = buffer.data(sih1 + 552);
    const auto *sih1_556 = buffer.data(sih1 + 556);
    const auto *sih1_558 = buffer.data(sih1 + 558);
    const auto *sih1_561 = buffer.data(sih1 + 561);
    const auto *sih1_563 = buffer.data(sih1 + 563);
    const auto *sih1_564 = buffer.data(sih1 + 564);
    const auto *sih1_567 = buffer.data(sih1 + 567);
    const auto *sih1_570 = buffer.data(sih1 + 570);
    const auto *sih1_572 = buffer.data(sih1 + 572);
    const auto *sih1_573 = buffer.data(sih1 + 573);
    const auto *sih1_576 = buffer.data(sih1 + 576);
    const auto *sih1_577 = buffer.data(sih1 + 577);
    const auto *sih1_579 = buffer.data(sih1 + 579);
    const auto *sih1_581 = buffer.data(sih1 + 581);
    const auto *sih1_582 = buffer.data(sih1 + 582);
    const auto *sih1_584 = buffer.data(sih1 + 584);
    const auto *sih1_585 = buffer.data(sih1 + 585);
    const auto *sih1_586 = buffer.data(sih1 + 586);
    const auto *sih1_587 = buffer.data(sih1 + 587);

    const auto *sii_727 = buffer.data(sii + 727);
    const auto *sii_728 = buffer.data(sii + 728);
    const auto *sii_730 = buffer.data(sii + 730);
    const auto *sii_731 = buffer.data(sii + 731);
    const auto *sii_733 = buffer.data(sii + 733);
    const auto *sii_734 = buffer.data(sii + 734);
    const auto *sii_737 = buffer.data(sii + 737);
    const auto *sii_738 = buffer.data(sii + 738);
    const auto *sii_740 = buffer.data(sii + 740);
    const auto *sii_742 = buffer.data(sii + 742);
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
    const auto *sii_758 = buffer.data(sii + 758);
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

#pragma omp simd aligned(t_934, t_935, t_936, t_937, pb_y, pc_y, pc_z, shk0_720, shi_531, \
                         shi_559, shi_560, shk1_720, sih0_545, sih1_545, sii_727, \
                         sii_728 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_934[k] = f_14 * shi_559[k]
                   + f_3 * pc_y[k] * sii_727[k];

        t_935[k] = f_16 * shi_531[k]
                   + f_1 * sih0_545[k]
                   - f_2 * sih1_545[k]
                   + f_3 * pc_z[k] * sii_727[k];

        t_936[k] = pb_y[k] * shk0_720[k]
                   - f_12 * pc_y[k] * shk1_720[k];

        t_937[k] = f_13 * shi_560[k]
                   + f_3 * pc_y[k] * sii_728[k];
    }

#pragma omp simd aligned(t_938, t_939, t_940, pc_x, pc_y, pc_z, shi_532, shi_562, sih0_549, \
                         sih1_549, sii_728, sii_730, sii_731 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_938[k] = f_17 * shi_532[k]
                   + f_3 * pc_z[k] * sii_728[k];

        t_939[k] = f_4 * sih0_549[k]
                   - f_5 * sih1_549[k]
                   + f_3 * pc_x[k] * sii_731[k];

        t_940[k] = f_13 * shi_562[k]
                   + f_3 * pc_y[k] * sii_730[k];
    }

#pragma omp simd aligned(t_941, t_942, t_943, pb_y, pc_x, pc_y, pc_z, shk0_725, shi_535, \
                         shk1_725, sih0_552, sih1_552, sii_731, \
                         sii_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_941[k] = pb_y[k] * shk0_725[k]
                   - f_12 * pc_y[k] * shk1_725[k];

        t_942[k] = f_6 * sih0_552[k]
                   - f_7 * sih1_552[k]
                   + f_3 * pc_x[k] * sii_734[k];

        t_943[k] = f_17 * shi_535[k]
                   + f_3 * pc_z[k] * sii_731[k];
    }

#pragma omp simd aligned(t_944, t_945, t_946, pb_y, pc_x, pc_y, shk0_729, shi_565, shk1_729, \
                         sih0_556, sih1_556, sii_733, sii_738 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_944[k] = f_13 * shi_565[k]
                   + f_3 * pc_y[k] * sii_733[k];

        t_945[k] = pb_y[k] * shk0_729[k]
                   - f_12 * pc_y[k] * shk1_729[k];

        t_946[k] = f_8 * sih0_556[k]
                   - f_9 * sih1_556[k]
                   + f_3 * pc_x[k] * sii_738[k];
    }

#pragma omp simd aligned(t_947, t_948, t_949, pc_x, pc_y, pc_z, shi_538, shi_569, sih0_558, \
                         sih1_558, sii_734, sii_737, sii_740 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_947[k] = f_17 * shi_538[k]
                   + f_3 * pc_z[k] * sii_734[k];

        t_948[k] = f_8 * sih0_558[k]
                   - f_9 * sih1_558[k]
                   + f_3 * pc_x[k] * sii_740[k];

        t_949[k] = f_13 * shi_569[k]
                   + f_3 * pc_y[k] * sii_737[k];
    }

#pragma omp simd aligned(t_950, t_951, t_952, pb_y, pc_x, pc_y, pc_z, shk0_734, shi_542, \
                         shk1_734, sih0_561, sih1_561, sii_738, \
                         sii_743 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_950[k] = pb_y[k] * shk0_734[k]
                   - f_12 * pc_y[k] * shk1_734[k];

        t_951[k] = f_10 * sih0_561[k]
                   - f_11 * sih1_561[k]
                   + f_3 * pc_x[k] * sii_743[k];

        t_952[k] = f_17 * shi_542[k]
                   + f_3 * pc_z[k] * sii_738[k];
    }

#pragma omp simd aligned(t_953, t_954, t_955, pc_x, pc_y, shi_574, sih0_563, sih0_564, \
                         sih1_563, sih1_564, sii_742, sii_745, \
                         sii_746 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_953[k] = f_10 * sih0_563[k]
                   - f_11 * sih1_563[k]
                   + f_3 * pc_x[k] * sii_745[k];

        t_954[k] = f_10 * sih0_564[k]
                   - f_11 * sih1_564[k]
                   + f_3 * pc_x[k] * sii_746[k];

        t_955[k] = f_13 * shi_574[k]
                   + f_3 * pc_y[k] * sii_742[k];
    }

#pragma omp simd aligned(t_956, t_957, t_958, t_959, t_960, t_961, pb_y, pc_x, pc_y, shk0_740, \
                         shk1_740, sii_749, sii_750, sii_751, sii_752, \
                         sii_753 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_956[k] = pb_y[k] * shk0_740[k]
                   - f_12 * pc_y[k] * shk1_740[k];

        t_957[k] = f_3 * pc_x[k] * sii_749[k];

        t_958[k] = f_3 * pc_x[k] * sii_750[k];

        t_959[k] = f_3 * pc_x[k] * sii_751[k];

        t_960[k] = f_3 * pc_x[k] * sii_752[k];

        t_961[k] = f_3 * pc_x[k] * sii_753[k];
    }

#pragma omp simd aligned(t_962, t_963, t_964, t_965, pb_y, pc_x, pc_y, pc_z, shk0_748, \
                         shi_553, shi_581, shk1_748, sii_749, sii_754, \
                         sii_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_962[k] = f_3 * pc_x[k] * sii_754[k];

        t_963[k] = f_3 * pc_x[k] * sii_755[k];

        t_964[k] = pb_y[k] * shk0_748[k]
                   + f_18 * shi_581[k]
                   - f_12 * pc_y[k] * shk1_748[k];

        t_965[k] = f_17 * shi_553[k]
                   + f_3 * pc_z[k] * sii_749[k];
    }

#pragma omp simd aligned(t_966, t_967, t_968, pb_y, pc_y, shk0_750, shk0_751, shk0_752, \
                         shi_583, shi_584, shi_585, shk1_750, shk1_751, \
                         shk1_752 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_966[k] = pb_y[k] * shk0_750[k]
                   + f_17 * shi_583[k]
                   - f_12 * pc_y[k] * shk1_750[k];

        t_967[k] = pb_y[k] * shk0_751[k]
                   + f_16 * shi_584[k]
                   - f_12 * pc_y[k] * shk1_751[k];

        t_968[k] = pb_y[k] * shk0_752[k]
                   + f_15 * shi_585[k]
                   - f_12 * pc_y[k] * shk1_752[k];
    }

#pragma omp simd aligned(t_969, t_970, t_971, pb_y, pc_y, shk0_753, shk0_755, shi_586, \
                         shi_587, shk1_753, shk1_755, sii_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_969[k] = pb_y[k] * shk0_753[k]
                   + f_14 * shi_586[k]
                   - f_12 * pc_y[k] * shk1_753[k];

        t_970[k] = f_13 * shi_587[k]
                   + f_3 * pc_y[k] * sii_755[k];

        t_971[k] = pb_y[k] * shk0_755[k]
                   - f_12 * pc_y[k] * shk1_755[k];
    }

#pragma omp simd aligned(t_972, t_973, t_974, t_975, t_976, pc_x, pc_y, pc_z, shi_560, \
                         sih0_567, sih0_570, sih1_567, sih1_570, sii_756, sii_758, \
                         sii_759 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_972[k] = f_1 * sih0_567[k]
                   - f_2 * sih1_567[k]
                   + f_3 * pc_x[k] * sii_756[k];

        t_973[k] = f_3 * pc_y[k] * sii_756[k];

        t_974[k] = f_0 * shi_560[k]
                   + f_3 * pc_z[k] * sii_756[k];

        t_975[k] = f_4 * sih0_570[k]
                   - f_5 * sih1_570[k]
                   + f_3 * pc_x[k] * sii_759[k];

        t_976[k] = f_3 * pc_y[k] * sii_758[k];
    }

#pragma omp simd aligned(t_977, t_978, t_979, t_980, pc_x, pc_y, pc_z, shi_563, sih0_572, \
                         sih0_573, sih1_572, sih1_573, sii_759, sii_761, \
                         sii_762 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_977[k] = f_4 * sih0_572[k]
                   - f_5 * sih1_572[k]
                   + f_3 * pc_x[k] * sii_761[k];

        t_978[k] = f_6 * sih0_573[k]
                   - f_7 * sih1_573[k]
                   + f_3 * pc_x[k] * sii_762[k];

        t_979[k] = f_0 * shi_563[k]
                   + f_3 * pc_z[k] * sii_759[k];

        t_980[k] = f_3 * pc_y[k] * sii_761[k];
    }

#pragma omp simd aligned(t_981, t_982, t_983, pc_x, pc_z, shi_566, sih0_576, sih0_577, \
                         sih1_576, sih1_577, sii_762, sii_765, \
                         sii_766 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_981[k] = f_6 * sih0_576[k]
                   - f_7 * sih1_576[k]
                   + f_3 * pc_x[k] * sii_765[k];

        t_982[k] = f_8 * sih0_577[k]
                   - f_9 * sih1_577[k]
                   + f_3 * pc_x[k] * sii_766[k];

        t_983[k] = f_0 * shi_566[k]
                   + f_3 * pc_z[k] * sii_762[k];
    }

#pragma omp simd aligned(t_984, t_985, t_986, t_987, pc_x, pc_y, sih0_579, sih0_581, sih0_582, \
                         sih1_579, sih1_581, sih1_582, sii_765, sii_768, sii_770, \
                         sii_771 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_984[k] = f_8 * sih0_579[k]
                   - f_9 * sih1_579[k]
                   + f_3 * pc_x[k] * sii_768[k];

        t_985[k] = f_3 * pc_y[k] * sii_765[k];

        t_986[k] = f_8 * sih0_581[k]
                   - f_9 * sih1_581[k]
                   + f_3 * pc_x[k] * sii_770[k];

        t_987[k] = f_10 * sih0_582[k]
                   - f_11 * sih1_582[k]
                   + f_3 * pc_x[k] * sii_771[k];
    }

#pragma omp simd aligned(t_988, t_989, t_990, t_991, pc_x, pc_y, pc_z, shi_570, sih0_584, \
                         sih0_585, sih1_584, sih1_585, sii_766, sii_770, sii_773, \
                         sii_774 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_988[k] = f_0 * shi_570[k]
                   + f_3 * pc_z[k] * sii_766[k];

        t_989[k] = f_10 * sih0_584[k]
                   - f_11 * sih1_584[k]
                   + f_3 * pc_x[k] * sii_773[k];

        t_990[k] = f_10 * sih0_585[k]
                   - f_11 * sih1_585[k]
                   + f_3 * pc_x[k] * sii_774[k];

        t_991[k] = f_3 * pc_y[k] * sii_770[k];
    }

#pragma omp simd aligned(t_992, t_993, t_994, t_995, t_996, t_997, pc_x, sih0_587, sih1_587, \
                         sii_776, sii_777, sii_778, sii_779, sii_780, \
                         sii_781 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_992[k] = f_10 * sih0_587[k]
                   - f_11 * sih1_587[k]
                   + f_3 * pc_x[k] * sii_776[k];

        t_993[k] = f_3 * pc_x[k] * sii_777[k];

        t_994[k] = f_3 * pc_x[k] * sii_778[k];

        t_995[k] = f_3 * pc_x[k] * sii_779[k];

        t_996[k] = f_3 * pc_x[k] * sii_780[k];

        t_997[k] = f_3 * pc_x[k] * sii_781[k];
    }

#pragma omp simd aligned(t_998, t_999, t_1000, t_1001, pc_x, pc_y, pc_z, shi_581, sih0_582, \
                         sih1_582, sii_777, sii_782, sii_783 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_998[k] = f_3 * pc_x[k] * sii_782[k];

        t_999[k] = f_3 * pc_x[k] * sii_783[k];

        t_1000[k] = f_1 * sih0_582[k]
                    - f_2 * sih1_582[k]
                    + f_3 * pc_y[k] * sii_777[k];

        t_1001[k] = f_0 * shi_581[k]
                    + f_3 * pc_z[k] * sii_777[k];
    }

#pragma omp simd aligned(t_1002, t_1003, t_1004, pc_y, sih0_584, sih0_585, sih0_586, sih1_584, \
                         sih1_585, sih1_586, sii_779, sii_780, \
                         sii_781 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1002[k] = f_4 * sih0_584[k]
                    - f_5 * sih1_584[k]
                    + f_3 * pc_y[k] * sii_779[k];

        t_1003[k] = f_6 * sih0_585[k]
                    - f_7 * sih1_585[k]
                    + f_3 * pc_y[k] * sii_780[k];

        t_1004[k] = f_8 * sih0_586[k]
                    - f_9 * sih1_586[k]
                    + f_3 * pc_y[k] * sii_781[k];
    }

#pragma omp simd aligned(t_1005, t_1006, t_1007, pc_y, pc_z, shi_587, sih0_587, sih1_587, \
                         sii_782, sii_783 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1005[k] = f_10 * sih0_587[k]
                    - f_11 * sih1_587[k]
                    + f_3 * pc_y[k] * sii_782[k];

        t_1006[k] = f_3 * pc_y[k] * sii_783[k];

        t_1007[k] = f_0 * shi_587[k]
                    + f_1 * sih0_587[k]
                    - f_2 * sih1_587[k]
                    + f_3 * pc_z[k] * sii_783[k];
    }
}

auto
compute_prim_sik_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t shk0, const size_t shi,
                                                   const size_t shk1, const size_t sih0,
                                                   const size_t sih1, const size_t sii,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_sik_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, shk0, shi,
                                                              shk1, sih0, sih1, sii, ncols,
                                                              gamma, p, q);

    compute_prim_sik_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, shk0, shi,
                                                              shk1, sih0, sih1, sii, ncols,
                                                              gamma, p, q);

    compute_prim_sik_three_center_electron_repulsion_0_piece2(buffer, target, pb, pc, shk0, shi,
                                                              shk1, sih0, sih1, sii, ncols,
                                                              gamma, p, q);

    compute_prim_sik_three_center_electron_repulsion_0_piece3(buffer, target, pb, pc, shk0, shi,
                                                              shk1, sih0, sih1, sii, ncols,
                                                              gamma, p, q);

    compute_prim_sik_three_center_electron_repulsion_0_piece4(buffer, target, pb, pc, shk0, shi,
                                                              shk1, sih0, sih1, sii, ncols,
                                                              gamma, p, q);

    compute_prim_sik_three_center_electron_repulsion_0_piece5(buffer, target, pb, pc, shk0, shi,
                                                              shk1, sii, ncols, gamma, p, q);

    compute_prim_sik_three_center_electron_repulsion_0_piece6(buffer, target, pb, pc, shk0, shi,
                                                              shk1, sih0, sih1, sii, ncols,
                                                              gamma, p, q);

    compute_prim_sik_three_center_electron_repulsion_0_piece7(buffer, target, pb, pc, shk0, shi,
                                                              shk1, sih0, sih1, sii, ncols,
                                                              gamma, p, q);

    compute_prim_sik_three_center_electron_repulsion_0_piece8(buffer, target, pb, pc, shk0, shi,
                                                              shk1, sih0, sih1, sii, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
