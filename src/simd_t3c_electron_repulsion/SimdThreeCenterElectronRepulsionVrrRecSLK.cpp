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


#include "SimdThreeCenterElectronRepulsionVrrRecSLK.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_slk_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t skk0,
                                                          const size_t ski, const size_t skk1,
                                                          const size_t slh0, const size_t slh1,
                                                          const size_t sli, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
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
    const auto f_19 = 3.0 / q;

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

    const auto *skk0_0 = buffer.data(skk0 + 0);
    const auto *skk0_3 = buffer.data(skk0 + 3);
    const auto *skk0_5 = buffer.data(skk0 + 5);
    const auto *skk0_6 = buffer.data(skk0 + 6);
    const auto *skk0_9 = buffer.data(skk0 + 9);
    const auto *skk0_10 = buffer.data(skk0 + 10);
    const auto *skk0_12 = buffer.data(skk0 + 12);
    const auto *skk0_14 = buffer.data(skk0 + 14);
    const auto *skk0_15 = buffer.data(skk0 + 15);
    const auto *skk0_17 = buffer.data(skk0 + 17);
    const auto *skk0_18 = buffer.data(skk0 + 18);
    const auto *skk0_20 = buffer.data(skk0 + 20);
    const auto *skk0_28 = buffer.data(skk0 + 28);
    const auto *skk0_35 = buffer.data(skk0 + 35);

    const auto *ski_0 = buffer.data(ski + 0);
    const auto *ski_1 = buffer.data(ski + 1);
    const auto *ski_2 = buffer.data(ski + 2);
    const auto *ski_3 = buffer.data(ski + 3);
    const auto *ski_5 = buffer.data(ski + 5);
    const auto *ski_6 = buffer.data(ski + 6);
    const auto *ski_7 = buffer.data(ski + 7);
    const auto *ski_8 = buffer.data(ski + 8);
    const auto *ski_9 = buffer.data(ski + 9);
    const auto *ski_10 = buffer.data(ski + 10);
    const auto *ski_11 = buffer.data(ski + 11);
    const auto *ski_12 = buffer.data(ski + 12);
    const auto *ski_13 = buffer.data(ski + 13);
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
    const auto *ski_33 = buffer.data(ski + 33);
    const auto *ski_37 = buffer.data(ski + 37);
    const auto *ski_49 = buffer.data(ski + 49);
    const auto *ski_50 = buffer.data(ski + 50);
    const auto *ski_51 = buffer.data(ski + 51);
    const auto *ski_52 = buffer.data(ski + 52);
    const auto *ski_53 = buffer.data(ski + 53);
    const auto *ski_54 = buffer.data(ski + 54);
    const auto *ski_55 = buffer.data(ski + 55);
    const auto *ski_77 = buffer.data(ski + 77);
    const auto *ski_78 = buffer.data(ski + 78);
    const auto *ski_79 = buffer.data(ski + 79);
    const auto *ski_80 = buffer.data(ski + 80);
    const auto *ski_81 = buffer.data(ski + 81);
    const auto *ski_82 = buffer.data(ski + 82);
    const auto *ski_83 = buffer.data(ski + 83);
    const auto *ski_84 = buffer.data(ski + 84);
    const auto *ski_87 = buffer.data(ski + 87);
    const auto *ski_89 = buffer.data(ski + 89);
    const auto *ski_90 = buffer.data(ski + 90);
    const auto *ski_93 = buffer.data(ski + 93);
    const auto *ski_94 = buffer.data(ski + 94);
    const auto *ski_96 = buffer.data(ski + 96);

    const auto *skk1_0 = buffer.data(skk1 + 0);
    const auto *skk1_3 = buffer.data(skk1 + 3);
    const auto *skk1_5 = buffer.data(skk1 + 5);
    const auto *skk1_6 = buffer.data(skk1 + 6);
    const auto *skk1_9 = buffer.data(skk1 + 9);
    const auto *skk1_10 = buffer.data(skk1 + 10);
    const auto *skk1_12 = buffer.data(skk1 + 12);
    const auto *skk1_14 = buffer.data(skk1 + 14);
    const auto *skk1_15 = buffer.data(skk1 + 15);
    const auto *skk1_17 = buffer.data(skk1 + 17);
    const auto *skk1_18 = buffer.data(skk1 + 18);
    const auto *skk1_20 = buffer.data(skk1 + 20);
    const auto *skk1_28 = buffer.data(skk1 + 28);
    const auto *skk1_35 = buffer.data(skk1 + 35);

    const auto *slh0_0 = buffer.data(slh0 + 0);
    const auto *slh0_3 = buffer.data(slh0 + 3);
    const auto *slh0_5 = buffer.data(slh0 + 5);
    const auto *slh0_6 = buffer.data(slh0 + 6);
    const auto *slh0_9 = buffer.data(slh0 + 9);
    const auto *slh0_10 = buffer.data(slh0 + 10);
    const auto *slh0_12 = buffer.data(slh0 + 12);
    const auto *slh0_14 = buffer.data(slh0 + 14);
    const auto *slh0_15 = buffer.data(slh0 + 15);
    const auto *slh0_17 = buffer.data(slh0 + 17);
    const auto *slh0_18 = buffer.data(slh0 + 18);
    const auto *slh0_19 = buffer.data(slh0 + 19);
    const auto *slh0_20 = buffer.data(slh0 + 20);
    const auto *slh0_36 = buffer.data(slh0 + 36);
    const auto *slh0_38 = buffer.data(slh0 + 38);
    const auto *slh0_39 = buffer.data(slh0 + 39);
    const auto *slh0_40 = buffer.data(slh0 + 40);
    const auto *slh0_41 = buffer.data(slh0 + 41);
    const auto *slh0_59 = buffer.data(slh0 + 59);
    const auto *slh0_60 = buffer.data(slh0 + 60);
    const auto *slh0_61 = buffer.data(slh0 + 61);
    const auto *slh0_62 = buffer.data(slh0 + 62);
    const auto *slh0_63 = buffer.data(slh0 + 63);
    const auto *slh0_66 = buffer.data(slh0 + 66);
    const auto *slh0_68 = buffer.data(slh0 + 68);
    const auto *slh0_69 = buffer.data(slh0 + 69);
    const auto *slh0_72 = buffer.data(slh0 + 72);
    const auto *slh0_73 = buffer.data(slh0 + 73);
    const auto *slh0_75 = buffer.data(slh0 + 75);

    const auto *slh1_0 = buffer.data(slh1 + 0);
    const auto *slh1_3 = buffer.data(slh1 + 3);
    const auto *slh1_5 = buffer.data(slh1 + 5);
    const auto *slh1_6 = buffer.data(slh1 + 6);
    const auto *slh1_9 = buffer.data(slh1 + 9);
    const auto *slh1_10 = buffer.data(slh1 + 10);
    const auto *slh1_12 = buffer.data(slh1 + 12);
    const auto *slh1_14 = buffer.data(slh1 + 14);
    const auto *slh1_15 = buffer.data(slh1 + 15);
    const auto *slh1_17 = buffer.data(slh1 + 17);
    const auto *slh1_18 = buffer.data(slh1 + 18);
    const auto *slh1_19 = buffer.data(slh1 + 19);
    const auto *slh1_20 = buffer.data(slh1 + 20);
    const auto *slh1_36 = buffer.data(slh1 + 36);
    const auto *slh1_38 = buffer.data(slh1 + 38);
    const auto *slh1_39 = buffer.data(slh1 + 39);
    const auto *slh1_40 = buffer.data(slh1 + 40);
    const auto *slh1_41 = buffer.data(slh1 + 41);
    const auto *slh1_59 = buffer.data(slh1 + 59);
    const auto *slh1_60 = buffer.data(slh1 + 60);
    const auto *slh1_61 = buffer.data(slh1 + 61);
    const auto *slh1_62 = buffer.data(slh1 + 62);
    const auto *slh1_63 = buffer.data(slh1 + 63);
    const auto *slh1_66 = buffer.data(slh1 + 66);
    const auto *slh1_68 = buffer.data(slh1 + 68);
    const auto *slh1_69 = buffer.data(slh1 + 69);
    const auto *slh1_72 = buffer.data(slh1 + 72);
    const auto *slh1_73 = buffer.data(slh1 + 73);
    const auto *slh1_75 = buffer.data(slh1 + 75);

    const auto *sli_0 = buffer.data(sli + 0);
    const auto *sli_2 = buffer.data(sli + 2);
    const auto *sli_3 = buffer.data(sli + 3);
    const auto *sli_5 = buffer.data(sli + 5);
    const auto *sli_6 = buffer.data(sli + 6);
    const auto *sli_9 = buffer.data(sli + 9);
    const auto *sli_10 = buffer.data(sli + 10);
    const auto *sli_12 = buffer.data(sli + 12);
    const auto *sli_14 = buffer.data(sli + 14);
    const auto *sli_15 = buffer.data(sli + 15);
    const auto *sli_17 = buffer.data(sli + 17);
    const auto *sli_18 = buffer.data(sli + 18);
    const auto *sli_20 = buffer.data(sli + 20);
    const auto *sli_21 = buffer.data(sli + 21);
    const auto *sli_22 = buffer.data(sli + 22);
    const auto *sli_23 = buffer.data(sli + 23);
    const auto *sli_24 = buffer.data(sli + 24);
    const auto *sli_25 = buffer.data(sli + 25);
    const auto *sli_26 = buffer.data(sli + 26);
    const auto *sli_27 = buffer.data(sli + 27);
    const auto *sli_28 = buffer.data(sli + 28);
    const auto *sli_30 = buffer.data(sli + 30);
    const auto *sli_31 = buffer.data(sli + 31);
    const auto *sli_33 = buffer.data(sli + 33);
    const auto *sli_34 = buffer.data(sli + 34);
    const auto *sli_37 = buffer.data(sli + 37);
    const auto *sli_38 = buffer.data(sli + 38);
    const auto *sli_42 = buffer.data(sli + 42);
    const auto *sli_49 = buffer.data(sli + 49);
    const auto *sli_50 = buffer.data(sli + 50);
    const auto *sli_51 = buffer.data(sli + 51);
    const auto *sli_52 = buffer.data(sli + 52);
    const auto *sli_53 = buffer.data(sli + 53);
    const auto *sli_54 = buffer.data(sli + 54);
    const auto *sli_55 = buffer.data(sli + 55);
    const auto *sli_56 = buffer.data(sli + 56);
    const auto *sli_58 = buffer.data(sli + 58);
    const auto *sli_59 = buffer.data(sli + 59);
    const auto *sli_61 = buffer.data(sli + 61);
    const auto *sli_62 = buffer.data(sli + 62);
    const auto *sli_65 = buffer.data(sli + 65);
    const auto *sli_66 = buffer.data(sli + 66);
    const auto *sli_70 = buffer.data(sli + 70);
    const auto *sli_77 = buffer.data(sli + 77);
    const auto *sli_78 = buffer.data(sli + 78);
    const auto *sli_79 = buffer.data(sli + 79);
    const auto *sli_80 = buffer.data(sli + 80);
    const auto *sli_81 = buffer.data(sli + 81);
    const auto *sli_82 = buffer.data(sli + 82);
    const auto *sli_83 = buffer.data(sli + 83);
    const auto *sli_84 = buffer.data(sli + 84);
    const auto *sli_86 = buffer.data(sli + 86);
    const auto *sli_87 = buffer.data(sli + 87);
    const auto *sli_89 = buffer.data(sli + 89);
    const auto *sli_90 = buffer.data(sli + 90);
    const auto *sli_93 = buffer.data(sli + 93);
    const auto *sli_94 = buffer.data(sli + 94);
    const auto *sli_96 = buffer.data(sli + 96);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, ski_0, ski_3, slh0_0, slh0_3, \
                         slh1_0, slh1_3, sli_0, sli_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ski_0[k]
                 + f_1 * slh0_0[k]
                 - f_2 * slh1_0[k]
                 + f_3 * pc_x[k] * sli_0[k];

        t_1[k] = f_3 * pc_y[k] * sli_0[k];

        t_2[k] = f_3 * pc_z[k] * sli_0[k];

        t_3[k] = f_0 * ski_3[k]
                 + f_4 * slh0_3[k]
                 - f_5 * slh1_3[k]
                 + f_3 * pc_x[k] * sli_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pc_x, pc_y, ski_5, ski_6, slh0_5, slh0_6, slh1_5, \
                         slh1_6, sli_2, sli_5, sli_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * pc_y[k] * sli_2[k];

        t_5[k] = f_0 * ski_5[k]
                 + f_4 * slh0_5[k]
                 - f_5 * slh1_5[k]
                 + f_3 * pc_x[k] * sli_5[k];

        t_6[k] = f_0 * ski_6[k]
                 + f_6 * slh0_6[k]
                 - f_7 * slh1_6[k]
                 + f_3 * pc_x[k] * sli_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pc_x, pc_y, pc_z, ski_9, slh0_9, slh1_9, sli_3, sli_5, \
                         sli_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_3 * pc_z[k] * sli_3[k];

        t_8[k] = f_3 * pc_y[k] * sli_5[k];

        t_9[k] = f_0 * ski_9[k]
                 + f_6 * slh0_9[k]
                 - f_7 * slh1_9[k]
                 + f_3 * pc_x[k] * sli_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pc_x, pc_z, ski_10, ski_12, slh0_10, slh0_12, \
                         slh1_10, slh1_12, sli_6, sli_10, sli_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * ski_10[k]
                  + f_8 * slh0_10[k]
                  - f_9 * slh1_10[k]
                  + f_3 * pc_x[k] * sli_10[k];

        t_11[k] = f_3 * pc_z[k] * sli_6[k];

        t_12[k] = f_0 * ski_12[k]
                  + f_8 * slh0_12[k]
                  - f_9 * slh1_12[k]
                  + f_3 * pc_x[k] * sli_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pc_x, pc_y, ski_14, ski_15, slh0_14, slh0_15, \
                         slh1_14, slh1_15, sli_9, sli_14, sli_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_3 * pc_y[k] * sli_9[k];

        t_14[k] = f_0 * ski_14[k]
                  + f_8 * slh0_14[k]
                  - f_9 * slh1_14[k]
                  + f_3 * pc_x[k] * sli_14[k];

        t_15[k] = f_0 * ski_15[k]
                  + f_10 * slh0_15[k]
                  - f_11 * slh1_15[k]
                  + f_3 * pc_x[k] * sli_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pc_x, pc_z, ski_17, ski_18, slh0_17, slh0_18, \
                         slh1_17, slh1_18, sli_10, sli_17, sli_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * pc_z[k] * sli_10[k];

        t_17[k] = f_0 * ski_17[k]
                  + f_10 * slh0_17[k]
                  - f_11 * slh1_17[k]
                  + f_3 * pc_x[k] * sli_17[k];

        t_18[k] = f_0 * ski_18[k]
                  + f_10 * slh0_18[k]
                  - f_11 * slh1_18[k]
                  + f_3 * pc_x[k] * sli_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pc_x, pc_y, ski_20, ski_21, ski_22, slh0_20, \
                         slh1_20, sli_14, sli_20, sli_21, sli_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * pc_y[k] * sli_14[k];

        t_20[k] = f_0 * ski_20[k]
                  + f_10 * slh0_20[k]
                  - f_11 * slh1_20[k]
                  + f_3 * pc_x[k] * sli_20[k];

        t_21[k] = f_0 * ski_21[k]
                  + f_3 * pc_x[k] * sli_21[k];

        t_22[k] = f_0 * ski_22[k]
                  + f_3 * pc_x[k] * sli_22[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, t_27, pc_x, ski_23, ski_24, ski_25, ski_26, \
                         ski_27, sli_23, sli_24, sli_25, sli_26, \
                         sli_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_0 * ski_23[k]
                  + f_3 * pc_x[k] * sli_23[k];

        t_24[k] = f_0 * ski_24[k]
                  + f_3 * pc_x[k] * sli_24[k];

        t_25[k] = f_0 * ski_25[k]
                  + f_3 * pc_x[k] * sli_25[k];

        t_26[k] = f_0 * ski_26[k]
                  + f_3 * pc_x[k] * sli_26[k];

        t_27[k] = f_0 * ski_27[k]
                  + f_3 * pc_x[k] * sli_27[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pc_y, pc_z, slh0_15, slh0_17, slh0_18, \
                         slh1_15, slh1_17, slh1_18, sli_21, sli_23, \
                         sli_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_1 * slh0_15[k]
                  - f_2 * slh1_15[k]
                  + f_3 * pc_y[k] * sli_21[k];

        t_29[k] = f_3 * pc_z[k] * sli_21[k];

        t_30[k] = f_4 * slh0_17[k]
                  - f_5 * slh1_17[k]
                  + f_3 * pc_y[k] * sli_23[k];

        t_31[k] = f_6 * slh0_18[k]
                  - f_7 * slh1_18[k]
                  + f_3 * pc_y[k] * sli_24[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pc_y, pc_z, slh0_19, slh0_20, slh1_19, \
                         slh1_20, sli_25, sli_26, sli_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_8 * slh0_19[k]
                  - f_9 * slh1_19[k]
                  + f_3 * pc_y[k] * sli_25[k];

        t_33[k] = f_10 * slh0_20[k]
                  - f_11 * slh1_20[k]
                  + f_3 * pc_y[k] * sli_26[k];

        t_34[k] = f_3 * pc_y[k] * sli_27[k];

        t_35[k] = f_1 * slh0_20[k]
                  - f_2 * slh1_20[k]
                  + f_3 * pc_z[k] * sli_27[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pb_y, pc_y, pc_z, skk0_0, skk0_3, ski_0, \
                         ski_1, skk1_0, skk1_3, sli_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = pb_y[k] * skk0_0[k]
                  - f_12 * pc_y[k] * skk1_0[k];

        t_37[k] = f_13 * ski_0[k]
                  + f_3 * pc_y[k] * sli_28[k];

        t_38[k] = f_3 * pc_z[k] * sli_28[k];

        t_39[k] = pb_y[k] * skk0_3[k]
                  + f_14 * ski_1[k]
                  - f_12 * pc_y[k] * skk1_3[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pb_y, pc_y, pc_z, skk0_5, skk0_6, ski_2, \
                         ski_3, skk1_5, skk1_6, sli_30, sli_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_13 * ski_2[k]
                  + f_3 * pc_y[k] * sli_30[k];

        t_41[k] = pb_y[k] * skk0_5[k]
                  - f_12 * pc_y[k] * skk1_5[k];

        t_42[k] = pb_y[k] * skk0_6[k]
                  + f_15 * ski_3[k]
                  - f_12 * pc_y[k] * skk1_6[k];

        t_43[k] = f_3 * pc_z[k] * sli_31[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pb_y, pc_y, pc_z, skk0_9, skk0_10, ski_5, \
                         ski_6, skk1_9, skk1_10, sli_33, sli_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_13 * ski_5[k]
                  + f_3 * pc_y[k] * sli_33[k];

        t_45[k] = pb_y[k] * skk0_9[k]
                  - f_12 * pc_y[k] * skk1_9[k];

        t_46[k] = pb_y[k] * skk0_10[k]
                  + f_16 * ski_6[k]
                  - f_12 * pc_y[k] * skk1_10[k];

        t_47[k] = f_3 * pc_z[k] * sli_34[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pb_y, pc_y, skk0_12, skk0_14, skk0_15, ski_8, \
                         ski_9, ski_10, skk1_12, skk1_14, skk1_15, \
                         sli_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = pb_y[k] * skk0_12[k]
                  + f_14 * ski_8[k]
                  - f_12 * pc_y[k] * skk1_12[k];

        t_49[k] = f_13 * ski_9[k]
                  + f_3 * pc_y[k] * sli_37[k];

        t_50[k] = pb_y[k] * skk0_14[k]
                  - f_12 * pc_y[k] * skk1_14[k];

        t_51[k] = pb_y[k] * skk0_15[k]
                  + f_17 * ski_10[k]
                  - f_12 * pc_y[k] * skk1_15[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pb_y, pc_y, pc_z, skk0_17, skk0_18, ski_12, \
                         ski_13, ski_14, skk1_17, skk1_18, sli_38, \
                         sli_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_3 * pc_z[k] * sli_38[k];

        t_53[k] = pb_y[k] * skk0_17[k]
                  + f_15 * ski_12[k]
                  - f_12 * pc_y[k] * skk1_17[k];

        t_54[k] = pb_y[k] * skk0_18[k]
                  + f_14 * ski_13[k]
                  - f_12 * pc_y[k] * skk1_18[k];

        t_55[k] = f_13 * ski_14[k]
                  + f_3 * pc_y[k] * sli_42[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pb_y, pc_x, pc_y, skk0_20, ski_49, ski_50, \
                         ski_51, skk1_20, sli_49, sli_50, sli_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_y[k] * skk0_20[k]
                  - f_12 * pc_y[k] * skk1_20[k];

        t_57[k] = f_18 * ski_49[k]
                  + f_3 * pc_x[k] * sli_49[k];

        t_58[k] = f_18 * ski_50[k]
                  + f_3 * pc_x[k] * sli_50[k];

        t_59[k] = f_18 * ski_51[k]
                  + f_3 * pc_x[k] * sli_51[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pc_x, ski_52, ski_53, ski_54, ski_55, sli_52, \
                         sli_53, sli_54, sli_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_18 * ski_52[k]
                  + f_3 * pc_x[k] * sli_52[k];

        t_61[k] = f_18 * ski_53[k]
                  + f_3 * pc_x[k] * sli_53[k];

        t_62[k] = f_18 * ski_54[k]
                  + f_3 * pc_x[k] * sli_54[k];

        t_63[k] = f_18 * ski_55[k]
                  + f_3 * pc_x[k] * sli_55[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, pc_y, pc_z, ski_21, ski_23, slh0_36, slh0_38, \
                         slh1_36, slh1_38, sli_49, sli_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_13 * ski_21[k]
                  + f_1 * slh0_36[k]
                  - f_2 * slh1_36[k]
                  + f_3 * pc_y[k] * sli_49[k];

        t_65[k] = f_3 * pc_z[k] * sli_49[k];

        t_66[k] = f_13 * ski_23[k]
                  + f_4 * slh0_38[k]
                  - f_5 * slh1_38[k]
                  + f_3 * pc_y[k] * sli_51[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, pc_y, ski_24, ski_25, ski_26, slh0_39, slh0_40, \
                         slh0_41, slh1_39, slh1_40, slh1_41, sli_52, sli_53, \
                         sli_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_13 * ski_24[k]
                  + f_6 * slh0_39[k]
                  - f_7 * slh1_39[k]
                  + f_3 * pc_y[k] * sli_52[k];

        t_68[k] = f_13 * ski_25[k]
                  + f_8 * slh0_40[k]
                  - f_9 * slh1_40[k]
                  + f_3 * pc_y[k] * sli_53[k];

        t_69[k] = f_13 * ski_26[k]
                  + f_10 * slh0_41[k]
                  - f_11 * slh1_41[k]
                  + f_3 * pc_y[k] * sli_54[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pb_y, pb_z, pc_y, pc_z, skk0_0, skk0_35, \
                         ski_27, skk1_0, skk1_35, sli_55, sli_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_13 * ski_27[k]
                  + f_3 * pc_y[k] * sli_55[k];

        t_71[k] = pb_y[k] * skk0_35[k]
                  - f_12 * pc_y[k] * skk1_35[k];

        t_72[k] = pb_z[k] * skk0_0[k]
                  - f_12 * pc_z[k] * skk1_0[k];

        t_73[k] = f_3 * pc_y[k] * sli_56[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pb_z, pc_y, pc_z, skk0_3, skk0_5, ski_0, \
                         ski_2, skk1_3, skk1_5, sli_56, sli_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_13 * ski_0[k]
                  + f_3 * pc_z[k] * sli_56[k];

        t_75[k] = pb_z[k] * skk0_3[k]
                  - f_12 * pc_z[k] * skk1_3[k];

        t_76[k] = f_3 * pc_y[k] * sli_58[k];

        t_77[k] = pb_z[k] * skk0_5[k]
                  + f_14 * ski_2[k]
                  - f_12 * pc_z[k] * skk1_5[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pb_z, pc_y, pc_z, skk0_6, skk0_9, ski_3, \
                         ski_5, skk1_6, skk1_9, sli_59, sli_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = pb_z[k] * skk0_6[k]
                  - f_12 * pc_z[k] * skk1_6[k];

        t_79[k] = f_13 * ski_3[k]
                  + f_3 * pc_z[k] * sli_59[k];

        t_80[k] = f_3 * pc_y[k] * sli_61[k];

        t_81[k] = pb_z[k] * skk0_9[k]
                  + f_15 * ski_5[k]
                  - f_12 * pc_z[k] * skk1_9[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, pb_z, pc_y, pc_z, skk0_10, skk0_12, ski_6, \
                         ski_7, skk1_10, skk1_12, sli_62, sli_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = pb_z[k] * skk0_10[k]
                  - f_12 * pc_z[k] * skk1_10[k];

        t_83[k] = f_13 * ski_6[k]
                  + f_3 * pc_z[k] * sli_62[k];

        t_84[k] = pb_z[k] * skk0_12[k]
                  + f_14 * ski_7[k]
                  - f_12 * pc_z[k] * skk1_12[k];

        t_85[k] = f_3 * pc_y[k] * sli_65[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pb_z, pc_z, skk0_14, skk0_15, skk0_17, ski_9, \
                         ski_10, ski_11, skk1_14, skk1_15, skk1_17, \
                         sli_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = pb_z[k] * skk0_14[k]
                  + f_16 * ski_9[k]
                  - f_12 * pc_z[k] * skk1_14[k];

        t_87[k] = pb_z[k] * skk0_15[k]
                  - f_12 * pc_z[k] * skk1_15[k];

        t_88[k] = f_13 * ski_10[k]
                  + f_3 * pc_z[k] * sli_66[k];

        t_89[k] = pb_z[k] * skk0_17[k]
                  + f_14 * ski_11[k]
                  - f_12 * pc_z[k] * skk1_17[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, pb_z, pc_y, pc_z, skk0_18, skk0_20, ski_12, ski_14, \
                         skk1_18, skk1_20, sli_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = pb_z[k] * skk0_18[k]
                  + f_15 * ski_12[k]
                  - f_12 * pc_z[k] * skk1_18[k];

        t_91[k] = f_3 * pc_y[k] * sli_70[k];

        t_92[k] = pb_z[k] * skk0_20[k]
                  + f_17 * ski_14[k]
                  - f_12 * pc_z[k] * skk1_20[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, t_97, pc_x, ski_77, ski_78, ski_79, ski_80, \
                         ski_81, sli_77, sli_78, sli_79, sli_80, \
                         sli_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_18 * ski_77[k]
                  + f_3 * pc_x[k] * sli_77[k];

        t_94[k] = f_18 * ski_78[k]
                  + f_3 * pc_x[k] * sli_78[k];

        t_95[k] = f_18 * ski_79[k]
                  + f_3 * pc_x[k] * sli_79[k];

        t_96[k] = f_18 * ski_80[k]
                  + f_3 * pc_x[k] * sli_80[k];

        t_97[k] = f_18 * ski_81[k]
                  + f_3 * pc_x[k] * sli_81[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, pb_z, pc_x, pc_z, skk0_28, ski_21, ski_82, \
                         ski_83, skk1_28, sli_77, sli_82, sli_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_18 * ski_82[k]
                  + f_3 * pc_x[k] * sli_82[k];

        t_99[k] = f_18 * ski_83[k]
                  + f_3 * pc_x[k] * sli_83[k];

        t_100[k] = pb_z[k] * skk0_28[k]
                   - f_12 * pc_z[k] * skk1_28[k];

        t_101[k] = f_13 * ski_21[k]
                   + f_3 * pc_z[k] * sli_77[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, pc_y, slh0_59, slh0_60, slh0_61, slh1_59, \
                         slh1_60, slh1_61, sli_79, sli_80, sli_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_4 * slh0_59[k]
                   - f_5 * slh1_59[k]
                   + f_3 * pc_y[k] * sli_79[k];

        t_103[k] = f_6 * slh0_60[k]
                   - f_7 * slh1_60[k]
                   + f_3 * pc_y[k] * sli_80[k];

        t_104[k] = f_8 * slh0_61[k]
                   - f_9 * slh1_61[k]
                   + f_3 * pc_y[k] * sli_81[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, pc_x, pc_y, pc_z, ski_27, ski_84, \
                         slh0_62, slh0_63, slh1_62, slh1_63, sli_82, sli_83, \
                         sli_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_10 * slh0_62[k]
                   - f_11 * slh1_62[k]
                   + f_3 * pc_y[k] * sli_82[k];

        t_106[k] = f_3 * pc_y[k] * sli_83[k];

        t_107[k] = f_13 * ski_27[k]
                   + f_1 * slh0_62[k]
                   - f_2 * slh1_62[k]
                   + f_3 * pc_z[k] * sli_83[k];

        t_108[k] = f_19 * ski_84[k]
                   + f_1 * slh0_63[k]
                   - f_2 * slh1_63[k]
                   + f_3 * pc_x[k] * sli_84[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pc_x, pc_y, pc_z, ski_28, ski_30, ski_87, \
                         slh0_66, slh1_66, sli_84, sli_86, sli_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_14 * ski_28[k]
                   + f_3 * pc_y[k] * sli_84[k];

        t_110[k] = f_3 * pc_z[k] * sli_84[k];

        t_111[k] = f_19 * ski_87[k]
                   + f_4 * slh0_66[k]
                   - f_5 * slh1_66[k]
                   + f_3 * pc_x[k] * sli_87[k];

        t_112[k] = f_14 * ski_30[k]
                   + f_3 * pc_y[k] * sli_86[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, pc_x, pc_z, ski_89, ski_90, slh0_68, slh0_69, \
                         slh1_68, slh1_69, sli_87, sli_89, sli_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_19 * ski_89[k]
                   + f_4 * slh0_68[k]
                   - f_5 * slh1_68[k]
                   + f_3 * pc_x[k] * sli_89[k];

        t_114[k] = f_19 * ski_90[k]
                   + f_6 * slh0_69[k]
                   - f_7 * slh1_69[k]
                   + f_3 * pc_x[k] * sli_90[k];

        t_115[k] = f_3 * pc_z[k] * sli_87[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, pc_x, pc_y, ski_33, ski_93, ski_94, slh0_72, \
                         slh0_73, slh1_72, slh1_73, sli_89, sli_93, \
                         sli_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_14 * ski_33[k]
                   + f_3 * pc_y[k] * sli_89[k];

        t_117[k] = f_19 * ski_93[k]
                   + f_6 * slh0_72[k]
                   - f_7 * slh1_72[k]
                   + f_3 * pc_x[k] * sli_93[k];

        t_118[k] = f_19 * ski_94[k]
                   + f_8 * slh0_73[k]
                   - f_9 * slh1_73[k]
                   + f_3 * pc_x[k] * sli_94[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, pc_x, pc_y, pc_z, ski_37, ski_96, slh0_75, \
                         slh1_75, sli_90, sli_93, sli_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = f_3 * pc_z[k] * sli_90[k];

        t_120[k] = f_19 * ski_96[k]
                   + f_8 * slh0_75[k]
                   - f_9 * slh1_75[k]
                   + f_3 * pc_x[k] * sli_96[k];

        t_121[k] = f_14 * ski_37[k]
                   + f_3 * pc_y[k] * sli_93[k];
    }
}

static auto
compute_prim_slk_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t skk0,
                                                          const size_t ski, const size_t skk1,
                                                          const size_t slh0, const size_t slh1,
                                                          const size_t sli, const size_t ncols,
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
    const auto f_17 = 2.5 / q;
    const auto f_19 = 3.0 / q;

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

    const auto *skk0_39 = buffer.data(skk0 + 39);
    const auto *skk0_42 = buffer.data(skk0 + 42);
    const auto *skk0_46 = buffer.data(skk0 + 46);
    const auto *skk0_51 = buffer.data(skk0 + 51);
    const auto *skk0_64 = buffer.data(skk0 + 64);
    const auto *skk0_72 = buffer.data(skk0 + 72);
    const auto *skk0_77 = buffer.data(skk0 + 77);
    const auto *skk0_81 = buffer.data(skk0 + 81);
    const auto *skk0_84 = buffer.data(skk0 + 84);
    const auto *skk0_86 = buffer.data(skk0 + 86);
    const auto *skk0_89 = buffer.data(skk0 + 89);
    const auto *skk0_90 = buffer.data(skk0 + 90);
    const auto *skk0_92 = buffer.data(skk0 + 92);
    const auto *skk0_107 = buffer.data(skk0 + 107);

    const auto *ski_28 = buffer.data(ski + 28);
    const auto *ski_31 = buffer.data(ski + 31);
    const auto *ski_34 = buffer.data(ski + 34);
    const auto *ski_38 = buffer.data(ski + 38);
    const auto *ski_42 = buffer.data(ski + 42);
    const auto *ski_49 = buffer.data(ski + 49);
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
    const auto *ski_64 = buffer.data(ski + 64);
    const auto *ski_65 = buffer.data(ski + 65);
    const auto *ski_66 = buffer.data(ski + 66);
    const auto *ski_68 = buffer.data(ski + 68);
    const auto *ski_69 = buffer.data(ski + 69);
    const auto *ski_70 = buffer.data(ski + 70);
    const auto *ski_77 = buffer.data(ski + 77);
    const auto *ski_79 = buffer.data(ski + 79);
    const auto *ski_80 = buffer.data(ski + 80);
    const auto *ski_81 = buffer.data(ski + 81);
    const auto *ski_82 = buffer.data(ski + 82);
    const auto *ski_83 = buffer.data(ski + 83);
    const auto *ski_84 = buffer.data(ski + 84);
    const auto *ski_86 = buffer.data(ski + 86);
    const auto *ski_89 = buffer.data(ski + 89);
    const auto *ski_93 = buffer.data(ski + 93);
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
    const auto *ski_133 = buffer.data(ski + 133);
    const auto *ski_134 = buffer.data(ski + 134);
    const auto *ski_135 = buffer.data(ski + 135);
    const auto *ski_136 = buffer.data(ski + 136);
    const auto *ski_137 = buffer.data(ski + 137);
    const auto *ski_138 = buffer.data(ski + 138);
    const auto *ski_139 = buffer.data(ski + 139);
    const auto *ski_140 = buffer.data(ski + 140);
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

    const auto *skk1_39 = buffer.data(skk1 + 39);
    const auto *skk1_42 = buffer.data(skk1 + 42);
    const auto *skk1_46 = buffer.data(skk1 + 46);
    const auto *skk1_51 = buffer.data(skk1 + 51);
    const auto *skk1_64 = buffer.data(skk1 + 64);
    const auto *skk1_72 = buffer.data(skk1 + 72);
    const auto *skk1_77 = buffer.data(skk1 + 77);
    const auto *skk1_81 = buffer.data(skk1 + 81);
    const auto *skk1_84 = buffer.data(skk1 + 84);
    const auto *skk1_86 = buffer.data(skk1 + 86);
    const auto *skk1_89 = buffer.data(skk1 + 89);
    const auto *skk1_90 = buffer.data(skk1 + 90);
    const auto *skk1_92 = buffer.data(skk1 + 92);
    const auto *skk1_107 = buffer.data(skk1 + 107);

    const auto *slh0_77 = buffer.data(slh0 + 77);
    const auto *slh0_78 = buffer.data(slh0 + 78);
    const auto *slh0_80 = buffer.data(slh0 + 80);
    const auto *slh0_81 = buffer.data(slh0 + 81);
    const auto *slh0_82 = buffer.data(slh0 + 82);
    const auto *slh0_83 = buffer.data(slh0 + 83);
    const auto *slh0_101 = buffer.data(slh0 + 101);
    const auto *slh0_102 = buffer.data(slh0 + 102);
    const auto *slh0_103 = buffer.data(slh0 + 103);
    const auto *slh0_104 = buffer.data(slh0 + 104);
    const auto *slh0_105 = buffer.data(slh0 + 105);
    const auto *slh0_108 = buffer.data(slh0 + 108);
    const auto *slh0_110 = buffer.data(slh0 + 110);
    const auto *slh0_111 = buffer.data(slh0 + 111);
    const auto *slh0_114 = buffer.data(slh0 + 114);
    const auto *slh0_115 = buffer.data(slh0 + 115);
    const auto *slh0_117 = buffer.data(slh0 + 117);
    const auto *slh0_119 = buffer.data(slh0 + 119);
    const auto *slh0_120 = buffer.data(slh0 + 120);
    const auto *slh0_122 = buffer.data(slh0 + 122);
    const auto *slh0_123 = buffer.data(slh0 + 123);
    const auto *slh0_124 = buffer.data(slh0 + 124);
    const auto *slh0_125 = buffer.data(slh0 + 125);
    const auto *slh0_126 = buffer.data(slh0 + 126);
    const auto *slh0_129 = buffer.data(slh0 + 129);
    const auto *slh0_131 = buffer.data(slh0 + 131);
    const auto *slh0_132 = buffer.data(slh0 + 132);
    const auto *slh0_135 = buffer.data(slh0 + 135);
    const auto *slh0_136 = buffer.data(slh0 + 136);
    const auto *slh0_138 = buffer.data(slh0 + 138);
    const auto *slh0_140 = buffer.data(slh0 + 140);
    const auto *slh0_141 = buffer.data(slh0 + 141);
    const auto *slh0_143 = buffer.data(slh0 + 143);
    const auto *slh0_144 = buffer.data(slh0 + 144);

    const auto *slh1_77 = buffer.data(slh1 + 77);
    const auto *slh1_78 = buffer.data(slh1 + 78);
    const auto *slh1_80 = buffer.data(slh1 + 80);
    const auto *slh1_81 = buffer.data(slh1 + 81);
    const auto *slh1_82 = buffer.data(slh1 + 82);
    const auto *slh1_83 = buffer.data(slh1 + 83);
    const auto *slh1_101 = buffer.data(slh1 + 101);
    const auto *slh1_102 = buffer.data(slh1 + 102);
    const auto *slh1_103 = buffer.data(slh1 + 103);
    const auto *slh1_104 = buffer.data(slh1 + 104);
    const auto *slh1_105 = buffer.data(slh1 + 105);
    const auto *slh1_108 = buffer.data(slh1 + 108);
    const auto *slh1_110 = buffer.data(slh1 + 110);
    const auto *slh1_111 = buffer.data(slh1 + 111);
    const auto *slh1_114 = buffer.data(slh1 + 114);
    const auto *slh1_115 = buffer.data(slh1 + 115);
    const auto *slh1_117 = buffer.data(slh1 + 117);
    const auto *slh1_119 = buffer.data(slh1 + 119);
    const auto *slh1_120 = buffer.data(slh1 + 120);
    const auto *slh1_122 = buffer.data(slh1 + 122);
    const auto *slh1_123 = buffer.data(slh1 + 123);
    const auto *slh1_124 = buffer.data(slh1 + 124);
    const auto *slh1_125 = buffer.data(slh1 + 125);
    const auto *slh1_126 = buffer.data(slh1 + 126);
    const auto *slh1_129 = buffer.data(slh1 + 129);
    const auto *slh1_131 = buffer.data(slh1 + 131);
    const auto *slh1_132 = buffer.data(slh1 + 132);
    const auto *slh1_135 = buffer.data(slh1 + 135);
    const auto *slh1_136 = buffer.data(slh1 + 136);
    const auto *slh1_138 = buffer.data(slh1 + 138);
    const auto *slh1_140 = buffer.data(slh1 + 140);
    const auto *slh1_141 = buffer.data(slh1 + 141);
    const auto *slh1_143 = buffer.data(slh1 + 143);
    const auto *slh1_144 = buffer.data(slh1 + 144);

    const auto *sli_94 = buffer.data(sli + 94);
    const auto *sli_98 = buffer.data(sli + 98);
    const auto *sli_99 = buffer.data(sli + 99);
    const auto *sli_101 = buffer.data(sli + 101);
    const auto *sli_102 = buffer.data(sli + 102);
    const auto *sli_104 = buffer.data(sli + 104);
    const auto *sli_105 = buffer.data(sli + 105);
    const auto *sli_106 = buffer.data(sli + 106);
    const auto *sli_107 = buffer.data(sli + 107);
    const auto *sli_108 = buffer.data(sli + 108);
    const auto *sli_109 = buffer.data(sli + 109);
    const auto *sli_110 = buffer.data(sli + 110);
    const auto *sli_111 = buffer.data(sli + 111);
    const auto *sli_112 = buffer.data(sli + 112);
    const auto *sli_114 = buffer.data(sli + 114);
    const auto *sli_115 = buffer.data(sli + 115);
    const auto *sli_117 = buffer.data(sli + 117);
    const auto *sli_118 = buffer.data(sli + 118);
    const auto *sli_121 = buffer.data(sli + 121);
    const auto *sli_122 = buffer.data(sli + 122);
    const auto *sli_126 = buffer.data(sli + 126);
    const auto *sli_133 = buffer.data(sli + 133);
    const auto *sli_134 = buffer.data(sli + 134);
    const auto *sli_135 = buffer.data(sli + 135);
    const auto *sli_136 = buffer.data(sli + 136);
    const auto *sli_137 = buffer.data(sli + 137);
    const auto *sli_138 = buffer.data(sli + 138);
    const auto *sli_139 = buffer.data(sli + 139);
    const auto *sli_140 = buffer.data(sli + 140);
    const auto *sli_142 = buffer.data(sli + 142);
    const auto *sli_143 = buffer.data(sli + 143);
    const auto *sli_145 = buffer.data(sli + 145);
    const auto *sli_146 = buffer.data(sli + 146);
    const auto *sli_149 = buffer.data(sli + 149);
    const auto *sli_150 = buffer.data(sli + 150);
    const auto *sli_152 = buffer.data(sli + 152);
    const auto *sli_154 = buffer.data(sli + 154);
    const auto *sli_155 = buffer.data(sli + 155);
    const auto *sli_157 = buffer.data(sli + 157);
    const auto *sli_158 = buffer.data(sli + 158);
    const auto *sli_160 = buffer.data(sli + 160);
    const auto *sli_161 = buffer.data(sli + 161);
    const auto *sli_162 = buffer.data(sli + 162);
    const auto *sli_163 = buffer.data(sli + 163);
    const auto *sli_164 = buffer.data(sli + 164);
    const auto *sli_165 = buffer.data(sli + 165);
    const auto *sli_166 = buffer.data(sli + 166);
    const auto *sli_167 = buffer.data(sli + 167);
    const auto *sli_168 = buffer.data(sli + 168);
    const auto *sli_170 = buffer.data(sli + 170);
    const auto *sli_171 = buffer.data(sli + 171);
    const auto *sli_173 = buffer.data(sli + 173);
    const auto *sli_174 = buffer.data(sli + 174);
    const auto *sli_177 = buffer.data(sli + 177);
    const auto *sli_178 = buffer.data(sli + 178);
    const auto *sli_180 = buffer.data(sli + 180);
    const auto *sli_182 = buffer.data(sli + 182);
    const auto *sli_183 = buffer.data(sli + 183);
    const auto *sli_185 = buffer.data(sli + 185);
    const auto *sli_186 = buffer.data(sli + 186);

#pragma omp simd aligned(t_122, t_123, t_124, pc_x, pc_z, ski_98, ski_99, slh0_77, slh0_78, \
                         slh1_77, slh1_78, sli_94, sli_98, sli_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_19 * ski_98[k]
                   + f_8 * slh0_77[k]
                   - f_9 * slh1_77[k]
                   + f_3 * pc_x[k] * sli_98[k];

        t_123[k] = f_19 * ski_99[k]
                   + f_10 * slh0_78[k]
                   - f_11 * slh1_78[k]
                   + f_3 * pc_x[k] * sli_99[k];

        t_124[k] = f_3 * pc_z[k] * sli_94[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, pc_x, pc_y, ski_42, ski_101, ski_102, slh0_80, \
                         slh0_81, slh1_80, slh1_81, sli_98, sli_101, \
                         sli_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_19 * ski_101[k]
                   + f_10 * slh0_80[k]
                   - f_11 * slh1_80[k]
                   + f_3 * pc_x[k] * sli_101[k];

        t_126[k] = f_19 * ski_102[k]
                   + f_10 * slh0_81[k]
                   - f_11 * slh1_81[k]
                   + f_3 * pc_x[k] * sli_102[k];

        t_127[k] = f_14 * ski_42[k]
                   + f_3 * pc_y[k] * sli_98[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, pc_x, ski_104, ski_105, ski_106, ski_107, \
                         slh0_83, slh1_83, sli_104, sli_105, sli_106, \
                         sli_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_19 * ski_104[k]
                   + f_10 * slh0_83[k]
                   - f_11 * slh1_83[k]
                   + f_3 * pc_x[k] * sli_104[k];

        t_129[k] = f_19 * ski_105[k]
                   + f_3 * pc_x[k] * sli_105[k];

        t_130[k] = f_19 * ski_106[k]
                   + f_3 * pc_x[k] * sli_106[k];

        t_131[k] = f_19 * ski_107[k]
                   + f_3 * pc_x[k] * sli_107[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pc_x, ski_108, ski_109, ski_110, ski_111, \
                         sli_108, sli_109, sli_110, sli_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_19 * ski_108[k]
                   + f_3 * pc_x[k] * sli_108[k];

        t_133[k] = f_19 * ski_109[k]
                   + f_3 * pc_x[k] * sli_109[k];

        t_134[k] = f_19 * ski_110[k]
                   + f_3 * pc_x[k] * sli_110[k];

        t_135[k] = f_19 * ski_111[k]
                   + f_3 * pc_x[k] * sli_111[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, pc_y, pc_z, ski_49, ski_51, slh0_78, slh0_80, \
                         slh1_78, slh1_80, sli_105, sli_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_14 * ski_49[k]
                   + f_1 * slh0_78[k]
                   - f_2 * slh1_78[k]
                   + f_3 * pc_y[k] * sli_105[k];

        t_137[k] = f_3 * pc_z[k] * sli_105[k];

        t_138[k] = f_14 * ski_51[k]
                   + f_4 * slh0_80[k]
                   - f_5 * slh1_80[k]
                   + f_3 * pc_y[k] * sli_107[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, pc_y, ski_52, ski_53, ski_54, slh0_81, slh0_82, \
                         slh0_83, slh1_81, slh1_82, slh1_83, sli_108, sli_109, \
                         sli_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_14 * ski_52[k]
                   + f_6 * slh0_81[k]
                   - f_7 * slh1_81[k]
                   + f_3 * pc_y[k] * sli_108[k];

        t_140[k] = f_14 * ski_53[k]
                   + f_8 * slh0_82[k]
                   - f_9 * slh1_82[k]
                   + f_3 * pc_y[k] * sli_109[k];

        t_141[k] = f_14 * ski_54[k]
                   + f_10 * slh0_83[k]
                   - f_11 * slh1_83[k]
                   + f_3 * pc_y[k] * sli_110[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, pb_y, pc_y, pc_z, skk0_72, ski_55, \
                         ski_56, skk1_72, slh0_83, slh1_83, sli_111, \
                         sli_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_14 * ski_55[k]
                   + f_3 * pc_y[k] * sli_111[k];

        t_143[k] = f_1 * slh0_83[k]
                   - f_2 * slh1_83[k]
                   + f_3 * pc_z[k] * sli_111[k];

        t_144[k] = pb_y[k] * skk0_72[k]
                   - f_12 * pc_y[k] * skk1_72[k];

        t_145[k] = f_13 * ski_56[k]
                   + f_3 * pc_y[k] * sli_112[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pb_y, pb_z, pc_y, pc_z, skk0_39, skk0_77, \
                         ski_28, ski_58, skk1_39, skk1_77, sli_112, \
                         sli_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_13 * ski_28[k]
                   + f_3 * pc_z[k] * sli_112[k];

        t_147[k] = pb_z[k] * skk0_39[k]
                   - f_12 * pc_z[k] * skk1_39[k];

        t_148[k] = f_13 * ski_58[k]
                   + f_3 * pc_y[k] * sli_114[k];

        t_149[k] = pb_y[k] * skk0_77[k]
                   - f_12 * pc_y[k] * skk1_77[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, pb_y, pb_z, pc_y, pc_z, skk0_42, skk0_81, \
                         ski_31, ski_61, skk1_42, skk1_81, sli_115, \
                         sli_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = pb_z[k] * skk0_42[k]
                   - f_12 * pc_z[k] * skk1_42[k];

        t_151[k] = f_13 * ski_31[k]
                   + f_3 * pc_z[k] * sli_115[k];

        t_152[k] = f_13 * ski_61[k]
                   + f_3 * pc_y[k] * sli_117[k];

        t_153[k] = pb_y[k] * skk0_81[k]
                   - f_12 * pc_y[k] * skk1_81[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, pb_y, pb_z, pc_y, pc_z, skk0_46, skk0_84, \
                         ski_34, ski_64, skk1_46, skk1_84, sli_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = pb_z[k] * skk0_46[k]
                   - f_12 * pc_z[k] * skk1_46[k];

        t_155[k] = f_13 * ski_34[k]
                   + f_3 * pc_z[k] * sli_118[k];

        t_156[k] = pb_y[k] * skk0_84[k]
                   + f_14 * ski_64[k]
                   - f_12 * pc_y[k] * skk1_84[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, pb_y, pb_z, pc_y, pc_z, skk0_51, skk0_86, \
                         ski_38, ski_65, skk1_51, skk1_86, sli_121, \
                         sli_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = f_13 * ski_65[k]
                   + f_3 * pc_y[k] * sli_121[k];

        t_158[k] = pb_y[k] * skk0_86[k]
                   - f_12 * pc_y[k] * skk1_86[k];

        t_159[k] = pb_z[k] * skk0_51[k]
                   - f_12 * pc_z[k] * skk1_51[k];

        t_160[k] = f_13 * ski_38[k]
                   + f_3 * pc_z[k] * sli_122[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, t_164, pb_y, pc_y, skk0_89, skk0_90, skk0_92, \
                         ski_68, ski_69, ski_70, skk1_89, skk1_90, skk1_92, \
                         sli_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = pb_y[k] * skk0_89[k]
                   + f_15 * ski_68[k]
                   - f_12 * pc_y[k] * skk1_89[k];

        t_162[k] = pb_y[k] * skk0_90[k]
                   + f_14 * ski_69[k]
                   - f_12 * pc_y[k] * skk1_90[k];

        t_163[k] = f_13 * ski_70[k]
                   + f_3 * pc_y[k] * sli_126[k];

        t_164[k] = pb_y[k] * skk0_92[k]
                   - f_12 * pc_y[k] * skk1_92[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, pc_x, ski_133, ski_134, ski_135, \
                         ski_136, ski_137, sli_133, sli_134, sli_135, sli_136, \
                         sli_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = f_19 * ski_133[k]
                   + f_3 * pc_x[k] * sli_133[k];

        t_166[k] = f_19 * ski_134[k]
                   + f_3 * pc_x[k] * sli_134[k];

        t_167[k] = f_19 * ski_135[k]
                   + f_3 * pc_x[k] * sli_135[k];

        t_168[k] = f_19 * ski_136[k]
                   + f_3 * pc_x[k] * sli_136[k];

        t_169[k] = f_19 * ski_137[k]
                   + f_3 * pc_x[k] * sli_137[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, pb_z, pc_x, pc_z, skk0_64, ski_49, \
                         ski_138, ski_139, skk1_64, sli_133, sli_138, \
                         sli_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_19 * ski_138[k]
                   + f_3 * pc_x[k] * sli_138[k];

        t_171[k] = f_19 * ski_139[k]
                   + f_3 * pc_x[k] * sli_139[k];

        t_172[k] = pb_z[k] * skk0_64[k]
                   - f_12 * pc_z[k] * skk1_64[k];

        t_173[k] = f_13 * ski_49[k]
                   + f_3 * pc_z[k] * sli_133[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, pc_y, ski_79, ski_80, ski_81, slh0_101, \
                         slh0_102, slh0_103, slh1_101, slh1_102, slh1_103, sli_135, sli_136, \
                         sli_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_13 * ski_79[k]
                   + f_4 * slh0_101[k]
                   - f_5 * slh1_101[k]
                   + f_3 * pc_y[k] * sli_135[k];

        t_175[k] = f_13 * ski_80[k]
                   + f_6 * slh0_102[k]
                   - f_7 * slh1_102[k]
                   + f_3 * pc_y[k] * sli_136[k];

        t_176[k] = f_13 * ski_81[k]
                   + f_8 * slh0_103[k]
                   - f_9 * slh1_103[k]
                   + f_3 * pc_y[k] * sli_137[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, pb_y, pc_y, skk0_107, ski_82, ski_83, skk1_107, \
                         slh0_104, slh1_104, sli_138, sli_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = f_13 * ski_82[k]
                   + f_10 * slh0_104[k]
                   - f_11 * slh1_104[k]
                   + f_3 * pc_y[k] * sli_138[k];

        t_178[k] = f_13 * ski_83[k]
                   + f_3 * pc_y[k] * sli_139[k];

        t_179[k] = pb_y[k] * skk0_107[k]
                   - f_12 * pc_y[k] * skk1_107[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pc_x, pc_y, pc_z, ski_56, ski_140, \
                         ski_143, slh0_105, slh0_108, slh1_105, slh1_108, sli_140, \
                         sli_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_19 * ski_140[k]
                   + f_1 * slh0_105[k]
                   - f_2 * slh1_105[k]
                   + f_3 * pc_x[k] * sli_140[k];

        t_181[k] = f_3 * pc_y[k] * sli_140[k];

        t_182[k] = f_14 * ski_56[k]
                   + f_3 * pc_z[k] * sli_140[k];

        t_183[k] = f_19 * ski_143[k]
                   + f_4 * slh0_108[k]
                   - f_5 * slh1_108[k]
                   + f_3 * pc_x[k] * sli_143[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, pc_x, pc_y, ski_145, ski_146, slh0_110, \
                         slh0_111, slh1_110, slh1_111, sli_142, sli_145, \
                         sli_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_3 * pc_y[k] * sli_142[k];

        t_185[k] = f_19 * ski_145[k]
                   + f_4 * slh0_110[k]
                   - f_5 * slh1_110[k]
                   + f_3 * pc_x[k] * sli_145[k];

        t_186[k] = f_19 * ski_146[k]
                   + f_6 * slh0_111[k]
                   - f_7 * slh1_111[k]
                   + f_3 * pc_x[k] * sli_146[k];
    }

#pragma omp simd aligned(t_187, t_188, t_189, pc_x, pc_y, pc_z, ski_59, ski_149, slh0_114, \
                         slh1_114, sli_143, sli_145, sli_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_187[k] = f_14 * ski_59[k]
                   + f_3 * pc_z[k] * sli_143[k];

        t_188[k] = f_3 * pc_y[k] * sli_145[k];

        t_189[k] = f_19 * ski_149[k]
                   + f_6 * slh0_114[k]
                   - f_7 * slh1_114[k]
                   + f_3 * pc_x[k] * sli_149[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, pc_x, pc_z, ski_62, ski_150, ski_152, slh0_115, \
                         slh0_117, slh1_115, slh1_117, sli_146, sli_150, \
                         sli_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = f_19 * ski_150[k]
                   + f_8 * slh0_115[k]
                   - f_9 * slh1_115[k]
                   + f_3 * pc_x[k] * sli_150[k];

        t_191[k] = f_14 * ski_62[k]
                   + f_3 * pc_z[k] * sli_146[k];

        t_192[k] = f_19 * ski_152[k]
                   + f_8 * slh0_117[k]
                   - f_9 * slh1_117[k]
                   + f_3 * pc_x[k] * sli_152[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, pc_x, pc_y, ski_154, ski_155, slh0_119, \
                         slh0_120, slh1_119, slh1_120, sli_149, sli_154, \
                         sli_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_3 * pc_y[k] * sli_149[k];

        t_194[k] = f_19 * ski_154[k]
                   + f_8 * slh0_119[k]
                   - f_9 * slh1_119[k]
                   + f_3 * pc_x[k] * sli_154[k];

        t_195[k] = f_19 * ski_155[k]
                   + f_10 * slh0_120[k]
                   - f_11 * slh1_120[k]
                   + f_3 * pc_x[k] * sli_155[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, pc_x, pc_z, ski_66, ski_157, ski_158, slh0_122, \
                         slh0_123, slh1_122, slh1_123, sli_150, sli_157, \
                         sli_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = f_14 * ski_66[k]
                   + f_3 * pc_z[k] * sli_150[k];

        t_197[k] = f_19 * ski_157[k]
                   + f_10 * slh0_122[k]
                   - f_11 * slh1_122[k]
                   + f_3 * pc_x[k] * sli_157[k];

        t_198[k] = f_19 * ski_158[k]
                   + f_10 * slh0_123[k]
                   - f_11 * slh1_123[k]
                   + f_3 * pc_x[k] * sli_158[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, pc_x, pc_y, ski_160, ski_161, ski_162, \
                         slh0_125, slh1_125, sli_154, sli_160, sli_161, \
                         sli_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = f_3 * pc_y[k] * sli_154[k];

        t_200[k] = f_19 * ski_160[k]
                   + f_10 * slh0_125[k]
                   - f_11 * slh1_125[k]
                   + f_3 * pc_x[k] * sli_160[k];

        t_201[k] = f_19 * ski_161[k]
                   + f_3 * pc_x[k] * sli_161[k];

        t_202[k] = f_19 * ski_162[k]
                   + f_3 * pc_x[k] * sli_162[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, t_207, pc_x, ski_163, ski_164, ski_165, \
                         ski_166, ski_167, sli_163, sli_164, sli_165, sli_166, \
                         sli_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_19 * ski_163[k]
                   + f_3 * pc_x[k] * sli_163[k];

        t_204[k] = f_19 * ski_164[k]
                   + f_3 * pc_x[k] * sli_164[k];

        t_205[k] = f_19 * ski_165[k]
                   + f_3 * pc_x[k] * sli_165[k];

        t_206[k] = f_19 * ski_166[k]
                   + f_3 * pc_x[k] * sli_166[k];

        t_207[k] = f_19 * ski_167[k]
                   + f_3 * pc_x[k] * sli_167[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, pc_y, pc_z, ski_77, slh0_120, slh0_122, \
                         slh0_123, slh1_120, slh1_122, slh1_123, sli_161, sli_163, \
                         sli_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_1 * slh0_120[k]
                   - f_2 * slh1_120[k]
                   + f_3 * pc_y[k] * sli_161[k];

        t_209[k] = f_14 * ski_77[k]
                   + f_3 * pc_z[k] * sli_161[k];

        t_210[k] = f_4 * slh0_122[k]
                   - f_5 * slh1_122[k]
                   + f_3 * pc_y[k] * sli_163[k];

        t_211[k] = f_6 * slh0_123[k]
                   - f_7 * slh1_123[k]
                   + f_3 * pc_y[k] * sli_164[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, pc_y, pc_z, ski_83, slh0_124, slh0_125, \
                         slh1_124, slh1_125, sli_165, sli_166, \
                         sli_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_8 * slh0_124[k]
                   - f_9 * slh1_124[k]
                   + f_3 * pc_y[k] * sli_165[k];

        t_213[k] = f_10 * slh0_125[k]
                   - f_11 * slh1_125[k]
                   + f_3 * pc_y[k] * sli_166[k];

        t_214[k] = f_3 * pc_y[k] * sli_167[k];

        t_215[k] = f_14 * ski_83[k]
                   + f_1 * slh0_125[k]
                   - f_2 * slh1_125[k]
                   + f_3 * pc_z[k] * sli_167[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, pc_x, pc_y, pc_z, ski_84, ski_168, \
                         ski_171, slh0_126, slh0_129, slh1_126, slh1_129, sli_168, \
                         sli_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = f_17 * ski_168[k]
                   + f_1 * slh0_126[k]
                   - f_2 * slh1_126[k]
                   + f_3 * pc_x[k] * sli_168[k];

        t_217[k] = f_15 * ski_84[k]
                   + f_3 * pc_y[k] * sli_168[k];

        t_218[k] = f_3 * pc_z[k] * sli_168[k];

        t_219[k] = f_17 * ski_171[k]
                   + f_4 * slh0_129[k]
                   - f_5 * slh1_129[k]
                   + f_3 * pc_x[k] * sli_171[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, pc_x, pc_y, ski_86, ski_173, ski_174, slh0_131, \
                         slh0_132, slh1_131, slh1_132, sli_170, sli_173, \
                         sli_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = f_15 * ski_86[k]
                   + f_3 * pc_y[k] * sli_170[k];

        t_221[k] = f_17 * ski_173[k]
                   + f_4 * slh0_131[k]
                   - f_5 * slh1_131[k]
                   + f_3 * pc_x[k] * sli_173[k];

        t_222[k] = f_17 * ski_174[k]
                   + f_6 * slh0_132[k]
                   - f_7 * slh1_132[k]
                   + f_3 * pc_x[k] * sli_174[k];
    }

#pragma omp simd aligned(t_223, t_224, t_225, pc_x, pc_y, pc_z, ski_89, ski_177, slh0_135, \
                         slh1_135, sli_171, sli_173, sli_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = f_3 * pc_z[k] * sli_171[k];

        t_224[k] = f_15 * ski_89[k]
                   + f_3 * pc_y[k] * sli_173[k];

        t_225[k] = f_17 * ski_177[k]
                   + f_6 * slh0_135[k]
                   - f_7 * slh1_135[k]
                   + f_3 * pc_x[k] * sli_177[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, pc_x, pc_z, ski_178, ski_180, slh0_136, \
                         slh0_138, slh1_136, slh1_138, sli_174, sli_178, \
                         sli_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = f_17 * ski_178[k]
                   + f_8 * slh0_136[k]
                   - f_9 * slh1_136[k]
                   + f_3 * pc_x[k] * sli_178[k];

        t_227[k] = f_3 * pc_z[k] * sli_174[k];

        t_228[k] = f_17 * ski_180[k]
                   + f_8 * slh0_138[k]
                   - f_9 * slh1_138[k]
                   + f_3 * pc_x[k] * sli_180[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, pc_x, pc_y, ski_93, ski_182, ski_183, slh0_140, \
                         slh0_141, slh1_140, slh1_141, sli_177, sli_182, \
                         sli_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = f_15 * ski_93[k]
                   + f_3 * pc_y[k] * sli_177[k];

        t_230[k] = f_17 * ski_182[k]
                   + f_8 * slh0_140[k]
                   - f_9 * slh1_140[k]
                   + f_3 * pc_x[k] * sli_182[k];

        t_231[k] = f_17 * ski_183[k]
                   + f_10 * slh0_141[k]
                   - f_11 * slh1_141[k]
                   + f_3 * pc_x[k] * sli_183[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, pc_x, pc_z, ski_185, ski_186, slh0_143, \
                         slh0_144, slh1_143, slh1_144, sli_178, sli_185, \
                         sli_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = f_3 * pc_z[k] * sli_178[k];

        t_233[k] = f_17 * ski_185[k]
                   + f_10 * slh0_143[k]
                   - f_11 * slh1_143[k]
                   + f_3 * pc_x[k] * sli_185[k];

        t_234[k] = f_17 * ski_186[k]
                   + f_10 * slh0_144[k]
                   - f_11 * slh1_144[k]
                   + f_3 * pc_x[k] * sli_186[k];
    }
}

static auto
compute_prim_slk_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t skk0,
                                                          const size_t ski, const size_t skk1,
                                                          const size_t slh0, const size_t slh1,
                                                          const size_t sli, const size_t ncols,
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

    const auto *skk0_108 = buffer.data(skk0 + 108);
    const auto *skk0_111 = buffer.data(skk0 + 111);
    const auto *skk0_114 = buffer.data(skk0 + 114);
    const auto *skk0_118 = buffer.data(skk0 + 118);
    const auto *skk0_120 = buffer.data(skk0 + 120);
    const auto *skk0_123 = buffer.data(skk0 + 123);
    const auto *skk0_125 = buffer.data(skk0 + 125);
    const auto *skk0_126 = buffer.data(skk0 + 126);
    const auto *skk0_136 = buffer.data(skk0 + 136);
    const auto *skk0_180 = buffer.data(skk0 + 180);
    const auto *skk0_183 = buffer.data(skk0 + 183);
    const auto *skk0_185 = buffer.data(skk0 + 185);
    const auto *skk0_186 = buffer.data(skk0 + 186);
    const auto *skk0_189 = buffer.data(skk0 + 189);
    const auto *skk0_190 = buffer.data(skk0 + 190);
    const auto *skk0_192 = buffer.data(skk0 + 192);
    const auto *skk0_194 = buffer.data(skk0 + 194);
    const auto *skk0_195 = buffer.data(skk0 + 195);
    const auto *skk0_197 = buffer.data(skk0 + 197);
    const auto *skk0_198 = buffer.data(skk0 + 198);
    const auto *skk0_200 = buffer.data(skk0 + 200);
    const auto *skk0_215 = buffer.data(skk0 + 215);

    const auto *ski_84 = buffer.data(ski + 84);
    const auto *ski_87 = buffer.data(ski + 87);
    const auto *ski_90 = buffer.data(ski + 90);
    const auto *ski_91 = buffer.data(ski + 91);
    const auto *ski_94 = buffer.data(ski + 94);
    const auto *ski_95 = buffer.data(ski + 95);
    const auto *ski_96 = buffer.data(ski + 96);
    const auto *ski_98 = buffer.data(ski + 98);
    const auto *ski_105 = buffer.data(ski + 105);
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
    const auto *ski_135 = buffer.data(ski + 135);
    const auto *ski_136 = buffer.data(ski + 136);
    const auto *ski_137 = buffer.data(ski + 137);
    const auto *ski_138 = buffer.data(ski + 138);
    const auto *ski_139 = buffer.data(ski + 139);
    const auto *ski_140 = buffer.data(ski + 140);
    const auto *ski_141 = buffer.data(ski + 141);
    const auto *ski_142 = buffer.data(ski + 142);
    const auto *ski_143 = buffer.data(ski + 143);
    const auto *ski_145 = buffer.data(ski + 145);
    const auto *ski_146 = buffer.data(ski + 146);
    const auto *ski_148 = buffer.data(ski + 148);
    const auto *ski_149 = buffer.data(ski + 149);
    const auto *ski_150 = buffer.data(ski + 150);
    const auto *ski_152 = buffer.data(ski + 152);
    const auto *ski_153 = buffer.data(ski + 153);
    const auto *ski_154 = buffer.data(ski + 154);
    const auto *ski_161 = buffer.data(ski + 161);
    const auto *ski_163 = buffer.data(ski + 163);
    const auto *ski_164 = buffer.data(ski + 164);
    const auto *ski_165 = buffer.data(ski + 165);
    const auto *ski_166 = buffer.data(ski + 166);
    const auto *ski_167 = buffer.data(ski + 167);
    const auto *ski_188 = buffer.data(ski + 188);
    const auto *ski_189 = buffer.data(ski + 189);
    const auto *ski_190 = buffer.data(ski + 190);
    const auto *ski_191 = buffer.data(ski + 191);
    const auto *ski_192 = buffer.data(ski + 192);
    const auto *ski_193 = buffer.data(ski + 193);
    const auto *ski_194 = buffer.data(ski + 194);
    const auto *ski_195 = buffer.data(ski + 195);
    const auto *ski_201 = buffer.data(ski + 201);
    const auto *ski_205 = buffer.data(ski + 205);
    const auto *ski_210 = buffer.data(ski + 210);
    const auto *ski_216 = buffer.data(ski + 216);
    const auto *ski_217 = buffer.data(ski + 217);
    const auto *ski_218 = buffer.data(ski + 218);
    const auto *ski_219 = buffer.data(ski + 219);
    const auto *ski_220 = buffer.data(ski + 220);
    const auto *ski_221 = buffer.data(ski + 221);
    const auto *ski_222 = buffer.data(ski + 222);
    const auto *ski_223 = buffer.data(ski + 223);
    const auto *ski_245 = buffer.data(ski + 245);
    const auto *ski_246 = buffer.data(ski + 246);
    const auto *ski_247 = buffer.data(ski + 247);
    const auto *ski_248 = buffer.data(ski + 248);
    const auto *ski_249 = buffer.data(ski + 249);
    const auto *ski_250 = buffer.data(ski + 250);
    const auto *ski_251 = buffer.data(ski + 251);
    const auto *ski_252 = buffer.data(ski + 252);
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

    const auto *skk1_108 = buffer.data(skk1 + 108);
    const auto *skk1_111 = buffer.data(skk1 + 111);
    const auto *skk1_114 = buffer.data(skk1 + 114);
    const auto *skk1_118 = buffer.data(skk1 + 118);
    const auto *skk1_120 = buffer.data(skk1 + 120);
    const auto *skk1_123 = buffer.data(skk1 + 123);
    const auto *skk1_125 = buffer.data(skk1 + 125);
    const auto *skk1_126 = buffer.data(skk1 + 126);
    const auto *skk1_136 = buffer.data(skk1 + 136);
    const auto *skk1_180 = buffer.data(skk1 + 180);
    const auto *skk1_183 = buffer.data(skk1 + 183);
    const auto *skk1_185 = buffer.data(skk1 + 185);
    const auto *skk1_186 = buffer.data(skk1 + 186);
    const auto *skk1_189 = buffer.data(skk1 + 189);
    const auto *skk1_190 = buffer.data(skk1 + 190);
    const auto *skk1_192 = buffer.data(skk1 + 192);
    const auto *skk1_194 = buffer.data(skk1 + 194);
    const auto *skk1_195 = buffer.data(skk1 + 195);
    const auto *skk1_197 = buffer.data(skk1 + 197);
    const auto *skk1_198 = buffer.data(skk1 + 198);
    const auto *skk1_200 = buffer.data(skk1 + 200);
    const auto *skk1_215 = buffer.data(skk1 + 215);

    const auto *slh0_141 = buffer.data(slh0 + 141);
    const auto *slh0_143 = buffer.data(slh0 + 143);
    const auto *slh0_144 = buffer.data(slh0 + 144);
    const auto *slh0_145 = buffer.data(slh0 + 145);
    const auto *slh0_146 = buffer.data(slh0 + 146);
    const auto *slh0_152 = buffer.data(slh0 + 152);
    const auto *slh0_156 = buffer.data(slh0 + 156);
    const auto *slh0_161 = buffer.data(slh0 + 161);
    const auto *slh0_164 = buffer.data(slh0 + 164);
    const auto *slh0_165 = buffer.data(slh0 + 165);
    const auto *slh0_166 = buffer.data(slh0 + 166);
    const auto *slh0_167 = buffer.data(slh0 + 167);
    const auto *slh0_183 = buffer.data(slh0 + 183);
    const auto *slh0_185 = buffer.data(slh0 + 185);
    const auto *slh0_186 = buffer.data(slh0 + 186);
    const auto *slh0_187 = buffer.data(slh0 + 187);
    const auto *slh0_188 = buffer.data(slh0 + 188);
    const auto *slh0_189 = buffer.data(slh0 + 189);
    const auto *slh0_192 = buffer.data(slh0 + 192);
    const auto *slh0_194 = buffer.data(slh0 + 194);
    const auto *slh0_195 = buffer.data(slh0 + 195);
    const auto *slh0_198 = buffer.data(slh0 + 198);
    const auto *slh0_199 = buffer.data(slh0 + 199);
    const auto *slh0_201 = buffer.data(slh0 + 201);
    const auto *slh0_203 = buffer.data(slh0 + 203);
    const auto *slh0_204 = buffer.data(slh0 + 204);
    const auto *slh0_206 = buffer.data(slh0 + 206);
    const auto *slh0_207 = buffer.data(slh0 + 207);
    const auto *slh0_209 = buffer.data(slh0 + 209);

    const auto *slh1_141 = buffer.data(slh1 + 141);
    const auto *slh1_143 = buffer.data(slh1 + 143);
    const auto *slh1_144 = buffer.data(slh1 + 144);
    const auto *slh1_145 = buffer.data(slh1 + 145);
    const auto *slh1_146 = buffer.data(slh1 + 146);
    const auto *slh1_152 = buffer.data(slh1 + 152);
    const auto *slh1_156 = buffer.data(slh1 + 156);
    const auto *slh1_161 = buffer.data(slh1 + 161);
    const auto *slh1_164 = buffer.data(slh1 + 164);
    const auto *slh1_165 = buffer.data(slh1 + 165);
    const auto *slh1_166 = buffer.data(slh1 + 166);
    const auto *slh1_167 = buffer.data(slh1 + 167);
    const auto *slh1_183 = buffer.data(slh1 + 183);
    const auto *slh1_185 = buffer.data(slh1 + 185);
    const auto *slh1_186 = buffer.data(slh1 + 186);
    const auto *slh1_187 = buffer.data(slh1 + 187);
    const auto *slh1_188 = buffer.data(slh1 + 188);
    const auto *slh1_189 = buffer.data(slh1 + 189);
    const auto *slh1_192 = buffer.data(slh1 + 192);
    const auto *slh1_194 = buffer.data(slh1 + 194);
    const auto *slh1_195 = buffer.data(slh1 + 195);
    const auto *slh1_198 = buffer.data(slh1 + 198);
    const auto *slh1_199 = buffer.data(slh1 + 199);
    const auto *slh1_201 = buffer.data(slh1 + 201);
    const auto *slh1_203 = buffer.data(slh1 + 203);
    const auto *slh1_204 = buffer.data(slh1 + 204);
    const auto *slh1_206 = buffer.data(slh1 + 206);
    const auto *slh1_207 = buffer.data(slh1 + 207);
    const auto *slh1_209 = buffer.data(slh1 + 209);

    const auto *sli_182 = buffer.data(sli + 182);
    const auto *sli_188 = buffer.data(sli + 188);
    const auto *sli_189 = buffer.data(sli + 189);
    const auto *sli_190 = buffer.data(sli + 190);
    const auto *sli_191 = buffer.data(sli + 191);
    const auto *sli_192 = buffer.data(sli + 192);
    const auto *sli_193 = buffer.data(sli + 193);
    const auto *sli_194 = buffer.data(sli + 194);
    const auto *sli_195 = buffer.data(sli + 195);
    const auto *sli_196 = buffer.data(sli + 196);
    const auto *sli_198 = buffer.data(sli + 198);
    const auto *sli_199 = buffer.data(sli + 199);
    const auto *sli_201 = buffer.data(sli + 201);
    const auto *sli_202 = buffer.data(sli + 202);
    const auto *sli_205 = buffer.data(sli + 205);
    const auto *sli_206 = buffer.data(sli + 206);
    const auto *sli_210 = buffer.data(sli + 210);
    const auto *sli_216 = buffer.data(sli + 216);
    const auto *sli_217 = buffer.data(sli + 217);
    const auto *sli_218 = buffer.data(sli + 218);
    const auto *sli_219 = buffer.data(sli + 219);
    const auto *sli_220 = buffer.data(sli + 220);
    const auto *sli_221 = buffer.data(sli + 221);
    const auto *sli_222 = buffer.data(sli + 222);
    const auto *sli_223 = buffer.data(sli + 223);
    const auto *sli_224 = buffer.data(sli + 224);
    const auto *sli_226 = buffer.data(sli + 226);
    const auto *sli_227 = buffer.data(sli + 227);
    const auto *sli_229 = buffer.data(sli + 229);
    const auto *sli_230 = buffer.data(sli + 230);
    const auto *sli_233 = buffer.data(sli + 233);
    const auto *sli_234 = buffer.data(sli + 234);
    const auto *sli_238 = buffer.data(sli + 238);
    const auto *sli_245 = buffer.data(sli + 245);
    const auto *sli_246 = buffer.data(sli + 246);
    const auto *sli_247 = buffer.data(sli + 247);
    const auto *sli_248 = buffer.data(sli + 248);
    const auto *sli_249 = buffer.data(sli + 249);
    const auto *sli_250 = buffer.data(sli + 250);
    const auto *sli_251 = buffer.data(sli + 251);
    const auto *sli_252 = buffer.data(sli + 252);
    const auto *sli_254 = buffer.data(sli + 254);
    const auto *sli_255 = buffer.data(sli + 255);
    const auto *sli_257 = buffer.data(sli + 257);
    const auto *sli_258 = buffer.data(sli + 258);
    const auto *sli_261 = buffer.data(sli + 261);
    const auto *sli_262 = buffer.data(sli + 262);
    const auto *sli_264 = buffer.data(sli + 264);
    const auto *sli_266 = buffer.data(sli + 266);
    const auto *sli_267 = buffer.data(sli + 267);
    const auto *sli_269 = buffer.data(sli + 269);
    const auto *sli_270 = buffer.data(sli + 270);
    const auto *sli_272 = buffer.data(sli + 272);
    const auto *sli_273 = buffer.data(sli + 273);
    const auto *sli_274 = buffer.data(sli + 274);

#pragma omp simd aligned(t_235, t_236, t_237, t_238, pc_x, pc_y, ski_98, ski_188, ski_189, \
                         ski_190, slh0_146, slh1_146, sli_182, sli_188, sli_189, \
                         sli_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = f_15 * ski_98[k]
                   + f_3 * pc_y[k] * sli_182[k];

        t_236[k] = f_17 * ski_188[k]
                   + f_10 * slh0_146[k]
                   - f_11 * slh1_146[k]
                   + f_3 * pc_x[k] * sli_188[k];

        t_237[k] = f_17 * ski_189[k]
                   + f_3 * pc_x[k] * sli_189[k];

        t_238[k] = f_17 * ski_190[k]
                   + f_3 * pc_x[k] * sli_190[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, t_243, pc_x, ski_191, ski_192, ski_193, \
                         ski_194, ski_195, sli_191, sli_192, sli_193, sli_194, \
                         sli_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_17 * ski_191[k]
                   + f_3 * pc_x[k] * sli_191[k];

        t_240[k] = f_17 * ski_192[k]
                   + f_3 * pc_x[k] * sli_192[k];

        t_241[k] = f_17 * ski_193[k]
                   + f_3 * pc_x[k] * sli_193[k];

        t_242[k] = f_17 * ski_194[k]
                   + f_3 * pc_x[k] * sli_194[k];

        t_243[k] = f_17 * ski_195[k]
                   + f_3 * pc_x[k] * sli_195[k];
    }

#pragma omp simd aligned(t_244, t_245, t_246, pc_y, pc_z, ski_105, ski_107, slh0_141, \
                         slh0_143, slh1_141, slh1_143, sli_189, \
                         sli_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_244[k] = f_15 * ski_105[k]
                   + f_1 * slh0_141[k]
                   - f_2 * slh1_141[k]
                   + f_3 * pc_y[k] * sli_189[k];

        t_245[k] = f_3 * pc_z[k] * sli_189[k];

        t_246[k] = f_15 * ski_107[k]
                   + f_4 * slh0_143[k]
                   - f_5 * slh1_143[k]
                   + f_3 * pc_y[k] * sli_191[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, pc_y, ski_108, ski_109, ski_110, slh0_144, \
                         slh0_145, slh0_146, slh1_144, slh1_145, slh1_146, sli_192, sli_193, \
                         sli_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = f_15 * ski_108[k]
                   + f_6 * slh0_144[k]
                   - f_7 * slh1_144[k]
                   + f_3 * pc_y[k] * sli_192[k];

        t_248[k] = f_15 * ski_109[k]
                   + f_8 * slh0_145[k]
                   - f_9 * slh1_145[k]
                   + f_3 * pc_y[k] * sli_193[k];

        t_249[k] = f_15 * ski_110[k]
                   + f_10 * slh0_146[k]
                   - f_11 * slh1_146[k]
                   + f_3 * pc_y[k] * sli_194[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, pb_z, pc_y, pc_z, skk0_108, ski_111, \
                         ski_112, skk1_108, slh0_146, slh1_146, sli_195, \
                         sli_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = f_15 * ski_111[k]
                   + f_3 * pc_y[k] * sli_195[k];

        t_251[k] = f_1 * slh0_146[k]
                   - f_2 * slh1_146[k]
                   + f_3 * pc_z[k] * sli_195[k];

        t_252[k] = pb_z[k] * skk0_108[k]
                   - f_12 * pc_z[k] * skk1_108[k];

        t_253[k] = f_14 * ski_112[k]
                   + f_3 * pc_y[k] * sli_196[k];
    }

#pragma omp simd aligned(t_254, t_255, t_256, pb_z, pc_y, pc_z, skk0_111, ski_84, ski_114, \
                         skk1_111, sli_196, sli_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_254[k] = f_13 * ski_84[k]
                   + f_3 * pc_z[k] * sli_196[k];

        t_255[k] = pb_z[k] * skk0_111[k]
                   - f_12 * pc_z[k] * skk1_111[k];

        t_256[k] = f_14 * ski_114[k]
                   + f_3 * pc_y[k] * sli_198[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, pb_z, pc_x, pc_z, skk0_114, ski_87, ski_201, \
                         skk1_114, slh0_152, slh1_152, sli_199, \
                         sli_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_17 * ski_201[k]
                   + f_4 * slh0_152[k]
                   - f_5 * slh1_152[k]
                   + f_3 * pc_x[k] * sli_201[k];

        t_258[k] = pb_z[k] * skk0_114[k]
                   - f_12 * pc_z[k] * skk1_114[k];

        t_259[k] = f_13 * ski_87[k]
                   + f_3 * pc_z[k] * sli_199[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, pb_z, pc_x, pc_y, pc_z, skk0_118, ski_117, \
                         ski_205, skk1_118, slh0_156, slh1_156, sli_201, \
                         sli_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_14 * ski_117[k]
                   + f_3 * pc_y[k] * sli_201[k];

        t_261[k] = f_17 * ski_205[k]
                   + f_6 * slh0_156[k]
                   - f_7 * slh1_156[k]
                   + f_3 * pc_x[k] * sli_205[k];

        t_262[k] = pb_z[k] * skk0_118[k]
                   - f_12 * pc_z[k] * skk1_118[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, pb_z, pc_y, pc_z, skk0_120, ski_90, ski_91, \
                         ski_121, skk1_120, sli_202, sli_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = f_13 * ski_90[k]
                   + f_3 * pc_z[k] * sli_202[k];

        t_264[k] = pb_z[k] * skk0_120[k]
                   + f_14 * ski_91[k]
                   - f_12 * pc_z[k] * skk1_120[k];

        t_265[k] = f_14 * ski_121[k]
                   + f_3 * pc_y[k] * sli_205[k];
    }

#pragma omp simd aligned(t_266, t_267, t_268, pb_z, pc_x, pc_z, skk0_123, ski_94, ski_210, \
                         skk1_123, slh0_161, slh1_161, sli_206, \
                         sli_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_266[k] = f_17 * ski_210[k]
                   + f_8 * slh0_161[k]
                   - f_9 * slh1_161[k]
                   + f_3 * pc_x[k] * sli_210[k];

        t_267[k] = pb_z[k] * skk0_123[k]
                   - f_12 * pc_z[k] * skk1_123[k];

        t_268[k] = f_13 * ski_94[k]
                   + f_3 * pc_z[k] * sli_206[k];
    }

#pragma omp simd aligned(t_269, t_270, t_271, pb_z, pc_y, pc_z, skk0_125, skk0_126, ski_95, \
                         ski_96, ski_126, skk1_125, skk1_126, sli_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_269[k] = pb_z[k] * skk0_125[k]
                   + f_14 * ski_95[k]
                   - f_12 * pc_z[k] * skk1_125[k];

        t_270[k] = pb_z[k] * skk0_126[k]
                   + f_15 * ski_96[k]
                   - f_12 * pc_z[k] * skk1_126[k];

        t_271[k] = f_14 * ski_126[k]
                   + f_3 * pc_y[k] * sli_210[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, pc_x, ski_216, ski_217, ski_218, ski_219, \
                         slh0_167, slh1_167, sli_216, sli_217, sli_218, \
                         sli_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_17 * ski_216[k]
                   + f_10 * slh0_167[k]
                   - f_11 * slh1_167[k]
                   + f_3 * pc_x[k] * sli_216[k];

        t_273[k] = f_17 * ski_217[k]
                   + f_3 * pc_x[k] * sli_217[k];

        t_274[k] = f_17 * ski_218[k]
                   + f_3 * pc_x[k] * sli_218[k];

        t_275[k] = f_17 * ski_219[k]
                   + f_3 * pc_x[k] * sli_219[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pc_x, ski_220, ski_221, ski_222, ski_223, \
                         sli_220, sli_221, sli_222, sli_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_17 * ski_220[k]
                   + f_3 * pc_x[k] * sli_220[k];

        t_277[k] = f_17 * ski_221[k]
                   + f_3 * pc_x[k] * sli_221[k];

        t_278[k] = f_17 * ski_222[k]
                   + f_3 * pc_x[k] * sli_222[k];

        t_279[k] = f_17 * ski_223[k]
                   + f_3 * pc_x[k] * sli_223[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, pb_z, pc_y, pc_z, skk0_136, ski_105, ski_135, \
                         skk1_136, slh0_164, slh1_164, sli_217, \
                         sli_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = pb_z[k] * skk0_136[k]
                   - f_12 * pc_z[k] * skk1_136[k];

        t_281[k] = f_13 * ski_105[k]
                   + f_3 * pc_z[k] * sli_217[k];

        t_282[k] = f_14 * ski_135[k]
                   + f_4 * slh0_164[k]
                   - f_5 * slh1_164[k]
                   + f_3 * pc_y[k] * sli_219[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, pc_y, ski_136, ski_137, ski_138, slh0_165, \
                         slh0_166, slh0_167, slh1_165, slh1_166, slh1_167, sli_220, sli_221, \
                         sli_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_14 * ski_136[k]
                   + f_6 * slh0_165[k]
                   - f_7 * slh1_165[k]
                   + f_3 * pc_y[k] * sli_220[k];

        t_284[k] = f_14 * ski_137[k]
                   + f_8 * slh0_166[k]
                   - f_9 * slh1_166[k]
                   + f_3 * pc_y[k] * sli_221[k];

        t_285[k] = f_14 * ski_138[k]
                   + f_10 * slh0_167[k]
                   - f_11 * slh1_167[k]
                   + f_3 * pc_y[k] * sli_222[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, pb_y, pc_y, pc_z, skk0_180, ski_111, \
                         ski_139, ski_140, skk1_180, slh0_167, slh1_167, sli_223, \
                         sli_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = f_14 * ski_139[k]
                   + f_3 * pc_y[k] * sli_223[k];

        t_287[k] = f_13 * ski_111[k]
                   + f_1 * slh0_167[k]
                   - f_2 * slh1_167[k]
                   + f_3 * pc_z[k] * sli_223[k];

        t_288[k] = pb_y[k] * skk0_180[k]
                   - f_12 * pc_y[k] * skk1_180[k];

        t_289[k] = f_13 * ski_140[k]
                   + f_3 * pc_y[k] * sli_224[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, pb_y, pc_y, pc_z, skk0_183, skk0_185, \
                         ski_112, ski_141, ski_142, skk1_183, skk1_185, sli_224, \
                         sli_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = f_14 * ski_112[k]
                   + f_3 * pc_z[k] * sli_224[k];

        t_291[k] = pb_y[k] * skk0_183[k]
                   + f_14 * ski_141[k]
                   - f_12 * pc_y[k] * skk1_183[k];

        t_292[k] = f_13 * ski_142[k]
                   + f_3 * pc_y[k] * sli_226[k];

        t_293[k] = pb_y[k] * skk0_185[k]
                   - f_12 * pc_y[k] * skk1_185[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, pb_y, pc_y, pc_z, skk0_186, skk0_189, \
                         ski_115, ski_143, ski_145, skk1_186, skk1_189, sli_227, \
                         sli_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = pb_y[k] * skk0_186[k]
                   + f_15 * ski_143[k]
                   - f_12 * pc_y[k] * skk1_186[k];

        t_295[k] = f_14 * ski_115[k]
                   + f_3 * pc_z[k] * sli_227[k];

        t_296[k] = f_13 * ski_145[k]
                   + f_3 * pc_y[k] * sli_229[k];

        t_297[k] = pb_y[k] * skk0_189[k]
                   - f_12 * pc_y[k] * skk1_189[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, pb_y, pc_y, pc_z, skk0_190, skk0_192, ski_118, \
                         ski_146, ski_148, skk1_190, skk1_192, \
                         sli_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = pb_y[k] * skk0_190[k]
                   + f_16 * ski_146[k]
                   - f_12 * pc_y[k] * skk1_190[k];

        t_299[k] = f_14 * ski_118[k]
                   + f_3 * pc_z[k] * sli_230[k];

        t_300[k] = pb_y[k] * skk0_192[k]
                   + f_14 * ski_148[k]
                   - f_12 * pc_y[k] * skk1_192[k];
    }

#pragma omp simd aligned(t_301, t_302, t_303, t_304, pb_y, pc_y, pc_z, skk0_194, skk0_195, \
                         ski_122, ski_149, ski_150, skk1_194, skk1_195, sli_233, \
                         sli_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_301[k] = f_13 * ski_149[k]
                   + f_3 * pc_y[k] * sli_233[k];

        t_302[k] = pb_y[k] * skk0_194[k]
                   - f_12 * pc_y[k] * skk1_194[k];

        t_303[k] = pb_y[k] * skk0_195[k]
                   + f_17 * ski_150[k]
                   - f_12 * pc_y[k] * skk1_195[k];

        t_304[k] = f_14 * ski_122[k]
                   + f_3 * pc_z[k] * sli_234[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, pb_y, pc_y, skk0_197, skk0_198, skk0_200, \
                         ski_152, ski_153, ski_154, skk1_197, skk1_198, skk1_200, \
                         sli_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = pb_y[k] * skk0_197[k]
                   + f_15 * ski_152[k]
                   - f_12 * pc_y[k] * skk1_197[k];

        t_306[k] = pb_y[k] * skk0_198[k]
                   + f_14 * ski_153[k]
                   - f_12 * pc_y[k] * skk1_198[k];

        t_307[k] = f_13 * ski_154[k]
                   + f_3 * pc_y[k] * sli_238[k];

        t_308[k] = pb_y[k] * skk0_200[k]
                   - f_12 * pc_y[k] * skk1_200[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, t_313, pc_x, ski_245, ski_246, ski_247, \
                         ski_248, ski_249, sli_245, sli_246, sli_247, sli_248, \
                         sli_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_17 * ski_245[k]
                   + f_3 * pc_x[k] * sli_245[k];

        t_310[k] = f_17 * ski_246[k]
                   + f_3 * pc_x[k] * sli_246[k];

        t_311[k] = f_17 * ski_247[k]
                   + f_3 * pc_x[k] * sli_247[k];

        t_312[k] = f_17 * ski_248[k]
                   + f_3 * pc_x[k] * sli_248[k];

        t_313[k] = f_17 * ski_249[k]
                   + f_3 * pc_x[k] * sli_249[k];
    }

#pragma omp simd aligned(t_314, t_315, t_316, t_317, pc_x, pc_y, pc_z, ski_133, ski_161, \
                         ski_250, ski_251, slh0_183, slh1_183, sli_245, sli_250, \
                         sli_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_314[k] = f_17 * ski_250[k]
                   + f_3 * pc_x[k] * sli_250[k];

        t_315[k] = f_17 * ski_251[k]
                   + f_3 * pc_x[k] * sli_251[k];

        t_316[k] = f_13 * ski_161[k]
                   + f_1 * slh0_183[k]
                   - f_2 * slh1_183[k]
                   + f_3 * pc_y[k] * sli_245[k];

        t_317[k] = f_14 * ski_133[k]
                   + f_3 * pc_z[k] * sli_245[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, pc_y, ski_163, ski_164, ski_165, slh0_185, \
                         slh0_186, slh0_187, slh1_185, slh1_186, slh1_187, sli_247, sli_248, \
                         sli_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = f_13 * ski_163[k]
                   + f_4 * slh0_185[k]
                   - f_5 * slh1_185[k]
                   + f_3 * pc_y[k] * sli_247[k];

        t_319[k] = f_13 * ski_164[k]
                   + f_6 * slh0_186[k]
                   - f_7 * slh1_186[k]
                   + f_3 * pc_y[k] * sli_248[k];

        t_320[k] = f_13 * ski_165[k]
                   + f_8 * slh0_187[k]
                   - f_9 * slh1_187[k]
                   + f_3 * pc_y[k] * sli_249[k];
    }

#pragma omp simd aligned(t_321, t_322, t_323, pb_y, pc_y, skk0_215, ski_166, ski_167, \
                         skk1_215, slh0_188, slh1_188, sli_250, \
                         sli_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_321[k] = f_13 * ski_166[k]
                   + f_10 * slh0_188[k]
                   - f_11 * slh1_188[k]
                   + f_3 * pc_y[k] * sli_250[k];

        t_322[k] = f_13 * ski_167[k]
                   + f_3 * pc_y[k] * sli_251[k];

        t_323[k] = pb_y[k] * skk0_215[k]
                   - f_12 * pc_y[k] * skk1_215[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, pc_x, pc_y, pc_z, ski_140, ski_252, \
                         ski_255, slh0_189, slh0_192, slh1_189, slh1_192, sli_252, \
                         sli_255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = f_17 * ski_252[k]
                   + f_1 * slh0_189[k]
                   - f_2 * slh1_189[k]
                   + f_3 * pc_x[k] * sli_252[k];

        t_325[k] = f_3 * pc_y[k] * sli_252[k];

        t_326[k] = f_15 * ski_140[k]
                   + f_3 * pc_z[k] * sli_252[k];

        t_327[k] = f_17 * ski_255[k]
                   + f_4 * slh0_192[k]
                   - f_5 * slh1_192[k]
                   + f_3 * pc_x[k] * sli_255[k];
    }

#pragma omp simd aligned(t_328, t_329, t_330, pc_x, pc_y, ski_257, ski_258, slh0_194, \
                         slh0_195, slh1_194, slh1_195, sli_254, sli_257, \
                         sli_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_328[k] = f_3 * pc_y[k] * sli_254[k];

        t_329[k] = f_17 * ski_257[k]
                   + f_4 * slh0_194[k]
                   - f_5 * slh1_194[k]
                   + f_3 * pc_x[k] * sli_257[k];

        t_330[k] = f_17 * ski_258[k]
                   + f_6 * slh0_195[k]
                   - f_7 * slh1_195[k]
                   + f_3 * pc_x[k] * sli_258[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, pc_x, pc_y, pc_z, ski_143, ski_261, slh0_198, \
                         slh1_198, sli_255, sli_257, sli_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = f_15 * ski_143[k]
                   + f_3 * pc_z[k] * sli_255[k];

        t_332[k] = f_3 * pc_y[k] * sli_257[k];

        t_333[k] = f_17 * ski_261[k]
                   + f_6 * slh0_198[k]
                   - f_7 * slh1_198[k]
                   + f_3 * pc_x[k] * sli_261[k];
    }

#pragma omp simd aligned(t_334, t_335, t_336, pc_x, pc_z, ski_146, ski_262, ski_264, slh0_199, \
                         slh0_201, slh1_199, slh1_201, sli_258, sli_262, \
                         sli_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_334[k] = f_17 * ski_262[k]
                   + f_8 * slh0_199[k]
                   - f_9 * slh1_199[k]
                   + f_3 * pc_x[k] * sli_262[k];

        t_335[k] = f_15 * ski_146[k]
                   + f_3 * pc_z[k] * sli_258[k];

        t_336[k] = f_17 * ski_264[k]
                   + f_8 * slh0_201[k]
                   - f_9 * slh1_201[k]
                   + f_3 * pc_x[k] * sli_264[k];
    }

#pragma omp simd aligned(t_337, t_338, t_339, pc_x, pc_y, ski_266, ski_267, slh0_203, \
                         slh0_204, slh1_203, slh1_204, sli_261, sli_266, \
                         sli_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_337[k] = f_3 * pc_y[k] * sli_261[k];

        t_338[k] = f_17 * ski_266[k]
                   + f_8 * slh0_203[k]
                   - f_9 * slh1_203[k]
                   + f_3 * pc_x[k] * sli_266[k];

        t_339[k] = f_17 * ski_267[k]
                   + f_10 * slh0_204[k]
                   - f_11 * slh1_204[k]
                   + f_3 * pc_x[k] * sli_267[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, pc_x, pc_z, ski_150, ski_269, ski_270, slh0_206, \
                         slh0_207, slh1_206, slh1_207, sli_262, sli_269, \
                         sli_270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = f_15 * ski_150[k]
                   + f_3 * pc_z[k] * sli_262[k];

        t_341[k] = f_17 * ski_269[k]
                   + f_10 * slh0_206[k]
                   - f_11 * slh1_206[k]
                   + f_3 * pc_x[k] * sli_269[k];

        t_342[k] = f_17 * ski_270[k]
                   + f_10 * slh0_207[k]
                   - f_11 * slh1_207[k]
                   + f_3 * pc_x[k] * sli_270[k];
    }

#pragma omp simd aligned(t_343, t_344, t_345, t_346, pc_x, pc_y, ski_272, ski_273, ski_274, \
                         slh0_209, slh1_209, sli_266, sli_272, sli_273, \
                         sli_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_343[k] = f_3 * pc_y[k] * sli_266[k];

        t_344[k] = f_17 * ski_272[k]
                   + f_10 * slh0_209[k]
                   - f_11 * slh1_209[k]
                   + f_3 * pc_x[k] * sli_272[k];

        t_345[k] = f_17 * ski_273[k]
                   + f_3 * pc_x[k] * sli_273[k];

        t_346[k] = f_17 * ski_274[k]
                   + f_3 * pc_x[k] * sli_274[k];
    }
}

static auto
compute_prim_slk_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t skk0,
                                                          const size_t ski, const size_t skk1,
                                                          const size_t slh0, const size_t slh1,
                                                          const size_t sli, const size_t ncols,
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

    const auto *skk0_216 = buffer.data(skk0 + 216);
    const auto *skk0_219 = buffer.data(skk0 + 219);
    const auto *skk0_222 = buffer.data(skk0 + 222);
    const auto *skk0_226 = buffer.data(skk0 + 226);
    const auto *skk0_228 = buffer.data(skk0 + 228);
    const auto *skk0_231 = buffer.data(skk0 + 231);
    const auto *skk0_233 = buffer.data(skk0 + 233);
    const auto *skk0_234 = buffer.data(skk0 + 234);
    const auto *skk0_244 = buffer.data(skk0 + 244);

    const auto *ski_161 = buffer.data(ski + 161);
    const auto *ski_167 = buffer.data(ski + 167);
    const auto *ski_168 = buffer.data(ski + 168);
    const auto *ski_170 = buffer.data(ski + 170);
    const auto *ski_171 = buffer.data(ski + 171);
    const auto *ski_173 = buffer.data(ski + 173);
    const auto *ski_174 = buffer.data(ski + 174);
    const auto *ski_175 = buffer.data(ski + 175);
    const auto *ski_177 = buffer.data(ski + 177);
    const auto *ski_178 = buffer.data(ski + 178);
    const auto *ski_179 = buffer.data(ski + 179);
    const auto *ski_180 = buffer.data(ski + 180);
    const auto *ski_182 = buffer.data(ski + 182);
    const auto *ski_189 = buffer.data(ski + 189);
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
    const auto *ski_219 = buffer.data(ski + 219);
    const auto *ski_220 = buffer.data(ski + 220);
    const auto *ski_221 = buffer.data(ski + 221);
    const auto *ski_222 = buffer.data(ski + 222);
    const auto *ski_223 = buffer.data(ski + 223);
    const auto *ski_224 = buffer.data(ski + 224);
    const auto *ski_226 = buffer.data(ski + 226);
    const auto *ski_229 = buffer.data(ski + 229);
    const auto *ski_233 = buffer.data(ski + 233);
    const auto *ski_238 = buffer.data(ski + 238);
    const auto *ski_275 = buffer.data(ski + 275);
    const auto *ski_276 = buffer.data(ski + 276);
    const auto *ski_277 = buffer.data(ski + 277);
    const auto *ski_278 = buffer.data(ski + 278);
    const auto *ski_279 = buffer.data(ski + 279);
    const auto *ski_280 = buffer.data(ski + 280);
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
    const auto *ski_313 = buffer.data(ski + 313);
    const auto *ski_317 = buffer.data(ski + 317);
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

    const auto *skk1_216 = buffer.data(skk1 + 216);
    const auto *skk1_219 = buffer.data(skk1 + 219);
    const auto *skk1_222 = buffer.data(skk1 + 222);
    const auto *skk1_226 = buffer.data(skk1 + 226);
    const auto *skk1_228 = buffer.data(skk1 + 228);
    const auto *skk1_231 = buffer.data(skk1 + 231);
    const auto *skk1_233 = buffer.data(skk1 + 233);
    const auto *skk1_234 = buffer.data(skk1 + 234);
    const auto *skk1_244 = buffer.data(skk1 + 244);

    const auto *slh0_204 = buffer.data(slh0 + 204);
    const auto *slh0_206 = buffer.data(slh0 + 206);
    const auto *slh0_207 = buffer.data(slh0 + 207);
    const auto *slh0_208 = buffer.data(slh0 + 208);
    const auto *slh0_209 = buffer.data(slh0 + 209);
    const auto *slh0_210 = buffer.data(slh0 + 210);
    const auto *slh0_213 = buffer.data(slh0 + 213);
    const auto *slh0_215 = buffer.data(slh0 + 215);
    const auto *slh0_216 = buffer.data(slh0 + 216);
    const auto *slh0_219 = buffer.data(slh0 + 219);
    const auto *slh0_220 = buffer.data(slh0 + 220);
    const auto *slh0_222 = buffer.data(slh0 + 222);
    const auto *slh0_224 = buffer.data(slh0 + 224);
    const auto *slh0_225 = buffer.data(slh0 + 225);
    const auto *slh0_227 = buffer.data(slh0 + 227);
    const auto *slh0_228 = buffer.data(slh0 + 228);
    const auto *slh0_229 = buffer.data(slh0 + 229);
    const auto *slh0_230 = buffer.data(slh0 + 230);
    const auto *slh0_236 = buffer.data(slh0 + 236);
    const auto *slh0_240 = buffer.data(slh0 + 240);
    const auto *slh0_245 = buffer.data(slh0 + 245);
    const auto *slh0_248 = buffer.data(slh0 + 248);
    const auto *slh0_249 = buffer.data(slh0 + 249);
    const auto *slh0_250 = buffer.data(slh0 + 250);
    const auto *slh0_251 = buffer.data(slh0 + 251);
    const auto *slh0_252 = buffer.data(slh0 + 252);
    const auto *slh0_255 = buffer.data(slh0 + 255);
    const auto *slh0_257 = buffer.data(slh0 + 257);
    const auto *slh0_258 = buffer.data(slh0 + 258);
    const auto *slh0_261 = buffer.data(slh0 + 261);
    const auto *slh0_262 = buffer.data(slh0 + 262);
    const auto *slh0_264 = buffer.data(slh0 + 264);
    const auto *slh0_266 = buffer.data(slh0 + 266);
    const auto *slh0_267 = buffer.data(slh0 + 267);
    const auto *slh0_269 = buffer.data(slh0 + 269);
    const auto *slh0_270 = buffer.data(slh0 + 270);
    const auto *slh0_272 = buffer.data(slh0 + 272);

    const auto *slh1_204 = buffer.data(slh1 + 204);
    const auto *slh1_206 = buffer.data(slh1 + 206);
    const auto *slh1_207 = buffer.data(slh1 + 207);
    const auto *slh1_208 = buffer.data(slh1 + 208);
    const auto *slh1_209 = buffer.data(slh1 + 209);
    const auto *slh1_210 = buffer.data(slh1 + 210);
    const auto *slh1_213 = buffer.data(slh1 + 213);
    const auto *slh1_215 = buffer.data(slh1 + 215);
    const auto *slh1_216 = buffer.data(slh1 + 216);
    const auto *slh1_219 = buffer.data(slh1 + 219);
    const auto *slh1_220 = buffer.data(slh1 + 220);
    const auto *slh1_222 = buffer.data(slh1 + 222);
    const auto *slh1_224 = buffer.data(slh1 + 224);
    const auto *slh1_225 = buffer.data(slh1 + 225);
    const auto *slh1_227 = buffer.data(slh1 + 227);
    const auto *slh1_228 = buffer.data(slh1 + 228);
    const auto *slh1_229 = buffer.data(slh1 + 229);
    const auto *slh1_230 = buffer.data(slh1 + 230);
    const auto *slh1_236 = buffer.data(slh1 + 236);
    const auto *slh1_240 = buffer.data(slh1 + 240);
    const auto *slh1_245 = buffer.data(slh1 + 245);
    const auto *slh1_248 = buffer.data(slh1 + 248);
    const auto *slh1_249 = buffer.data(slh1 + 249);
    const auto *slh1_250 = buffer.data(slh1 + 250);
    const auto *slh1_251 = buffer.data(slh1 + 251);
    const auto *slh1_252 = buffer.data(slh1 + 252);
    const auto *slh1_255 = buffer.data(slh1 + 255);
    const auto *slh1_257 = buffer.data(slh1 + 257);
    const auto *slh1_258 = buffer.data(slh1 + 258);
    const auto *slh1_261 = buffer.data(slh1 + 261);
    const auto *slh1_262 = buffer.data(slh1 + 262);
    const auto *slh1_264 = buffer.data(slh1 + 264);
    const auto *slh1_266 = buffer.data(slh1 + 266);
    const auto *slh1_267 = buffer.data(slh1 + 267);
    const auto *slh1_269 = buffer.data(slh1 + 269);
    const auto *slh1_270 = buffer.data(slh1 + 270);
    const auto *slh1_272 = buffer.data(slh1 + 272);

    const auto *sli_273 = buffer.data(sli + 273);
    const auto *sli_275 = buffer.data(sli + 275);
    const auto *sli_276 = buffer.data(sli + 276);
    const auto *sli_277 = buffer.data(sli + 277);
    const auto *sli_278 = buffer.data(sli + 278);
    const auto *sli_279 = buffer.data(sli + 279);
    const auto *sli_280 = buffer.data(sli + 280);
    const auto *sli_282 = buffer.data(sli + 282);
    const auto *sli_283 = buffer.data(sli + 283);
    const auto *sli_285 = buffer.data(sli + 285);
    const auto *sli_286 = buffer.data(sli + 286);
    const auto *sli_289 = buffer.data(sli + 289);
    const auto *sli_290 = buffer.data(sli + 290);
    const auto *sli_292 = buffer.data(sli + 292);
    const auto *sli_294 = buffer.data(sli + 294);
    const auto *sli_295 = buffer.data(sli + 295);
    const auto *sli_297 = buffer.data(sli + 297);
    const auto *sli_298 = buffer.data(sli + 298);
    const auto *sli_300 = buffer.data(sli + 300);
    const auto *sli_301 = buffer.data(sli + 301);
    const auto *sli_302 = buffer.data(sli + 302);
    const auto *sli_303 = buffer.data(sli + 303);
    const auto *sli_304 = buffer.data(sli + 304);
    const auto *sli_305 = buffer.data(sli + 305);
    const auto *sli_306 = buffer.data(sli + 306);
    const auto *sli_307 = buffer.data(sli + 307);
    const auto *sli_308 = buffer.data(sli + 308);
    const auto *sli_310 = buffer.data(sli + 310);
    const auto *sli_311 = buffer.data(sli + 311);
    const auto *sli_313 = buffer.data(sli + 313);
    const auto *sli_314 = buffer.data(sli + 314);
    const auto *sli_317 = buffer.data(sli + 317);
    const auto *sli_318 = buffer.data(sli + 318);
    const auto *sli_322 = buffer.data(sli + 322);
    const auto *sli_328 = buffer.data(sli + 328);
    const auto *sli_329 = buffer.data(sli + 329);
    const auto *sli_330 = buffer.data(sli + 330);
    const auto *sli_331 = buffer.data(sli + 331);
    const auto *sli_332 = buffer.data(sli + 332);
    const auto *sli_333 = buffer.data(sli + 333);
    const auto *sli_334 = buffer.data(sli + 334);
    const auto *sli_335 = buffer.data(sli + 335);
    const auto *sli_336 = buffer.data(sli + 336);
    const auto *sli_338 = buffer.data(sli + 338);
    const auto *sli_339 = buffer.data(sli + 339);
    const auto *sli_341 = buffer.data(sli + 341);
    const auto *sli_342 = buffer.data(sli + 342);
    const auto *sli_345 = buffer.data(sli + 345);
    const auto *sli_346 = buffer.data(sli + 346);
    const auto *sli_348 = buffer.data(sli + 348);
    const auto *sli_350 = buffer.data(sli + 350);
    const auto *sli_351 = buffer.data(sli + 351);
    const auto *sli_353 = buffer.data(sli + 353);
    const auto *sli_354 = buffer.data(sli + 354);
    const auto *sli_356 = buffer.data(sli + 356);
    const auto *sli_357 = buffer.data(sli + 357);
    const auto *sli_358 = buffer.data(sli + 358);
    const auto *sli_359 = buffer.data(sli + 359);

#pragma omp simd aligned(t_347, t_348, t_349, t_350, t_351, pc_x, ski_275, ski_276, ski_277, \
                         ski_278, ski_279, sli_275, sli_276, sli_277, sli_278, \
                         sli_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_347[k] = f_17 * ski_275[k]
                   + f_3 * pc_x[k] * sli_275[k];

        t_348[k] = f_17 * ski_276[k]
                   + f_3 * pc_x[k] * sli_276[k];

        t_349[k] = f_17 * ski_277[k]
                   + f_3 * pc_x[k] * sli_277[k];

        t_350[k] = f_17 * ski_278[k]
                   + f_3 * pc_x[k] * sli_278[k];

        t_351[k] = f_17 * ski_279[k]
                   + f_3 * pc_x[k] * sli_279[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, t_355, pc_y, pc_z, ski_161, slh0_204, slh0_206, \
                         slh0_207, slh1_204, slh1_206, slh1_207, sli_273, sli_275, \
                         sli_276 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = f_1 * slh0_204[k]
                   - f_2 * slh1_204[k]
                   + f_3 * pc_y[k] * sli_273[k];

        t_353[k] = f_15 * ski_161[k]
                   + f_3 * pc_z[k] * sli_273[k];

        t_354[k] = f_4 * slh0_206[k]
                   - f_5 * slh1_206[k]
                   + f_3 * pc_y[k] * sli_275[k];

        t_355[k] = f_6 * slh0_207[k]
                   - f_7 * slh1_207[k]
                   + f_3 * pc_y[k] * sli_276[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, t_359, pc_y, pc_z, ski_167, slh0_208, slh0_209, \
                         slh1_208, slh1_209, sli_277, sli_278, \
                         sli_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_8 * slh0_208[k]
                   - f_9 * slh1_208[k]
                   + f_3 * pc_y[k] * sli_277[k];

        t_357[k] = f_10 * slh0_209[k]
                   - f_11 * slh1_209[k]
                   + f_3 * pc_y[k] * sli_278[k];

        t_358[k] = f_3 * pc_y[k] * sli_279[k];

        t_359[k] = f_15 * ski_167[k]
                   + f_1 * slh0_209[k]
                   - f_2 * slh1_209[k]
                   + f_3 * pc_z[k] * sli_279[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, pc_x, pc_y, pc_z, ski_168, ski_280, \
                         ski_283, slh0_210, slh0_213, slh1_210, slh1_213, sli_280, \
                         sli_283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = f_16 * ski_280[k]
                   + f_1 * slh0_210[k]
                   - f_2 * slh1_210[k]
                   + f_3 * pc_x[k] * sli_280[k];

        t_361[k] = f_16 * ski_168[k]
                   + f_3 * pc_y[k] * sli_280[k];

        t_362[k] = f_3 * pc_z[k] * sli_280[k];

        t_363[k] = f_16 * ski_283[k]
                   + f_4 * slh0_213[k]
                   - f_5 * slh1_213[k]
                   + f_3 * pc_x[k] * sli_283[k];
    }

#pragma omp simd aligned(t_364, t_365, t_366, pc_x, pc_y, ski_170, ski_285, ski_286, slh0_215, \
                         slh0_216, slh1_215, slh1_216, sli_282, sli_285, \
                         sli_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_364[k] = f_16 * ski_170[k]
                   + f_3 * pc_y[k] * sli_282[k];

        t_365[k] = f_16 * ski_285[k]
                   + f_4 * slh0_215[k]
                   - f_5 * slh1_215[k]
                   + f_3 * pc_x[k] * sli_285[k];

        t_366[k] = f_16 * ski_286[k]
                   + f_6 * slh0_216[k]
                   - f_7 * slh1_216[k]
                   + f_3 * pc_x[k] * sli_286[k];
    }

#pragma omp simd aligned(t_367, t_368, t_369, pc_x, pc_y, pc_z, ski_173, ski_289, slh0_219, \
                         slh1_219, sli_283, sli_285, sli_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_367[k] = f_3 * pc_z[k] * sli_283[k];

        t_368[k] = f_16 * ski_173[k]
                   + f_3 * pc_y[k] * sli_285[k];

        t_369[k] = f_16 * ski_289[k]
                   + f_6 * slh0_219[k]
                   - f_7 * slh1_219[k]
                   + f_3 * pc_x[k] * sli_289[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, pc_x, pc_z, ski_290, ski_292, slh0_220, \
                         slh0_222, slh1_220, slh1_222, sli_286, sli_290, \
                         sli_292 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = f_16 * ski_290[k]
                   + f_8 * slh0_220[k]
                   - f_9 * slh1_220[k]
                   + f_3 * pc_x[k] * sli_290[k];

        t_371[k] = f_3 * pc_z[k] * sli_286[k];

        t_372[k] = f_16 * ski_292[k]
                   + f_8 * slh0_222[k]
                   - f_9 * slh1_222[k]
                   + f_3 * pc_x[k] * sli_292[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, pc_x, pc_y, ski_177, ski_294, ski_295, slh0_224, \
                         slh0_225, slh1_224, slh1_225, sli_289, sli_294, \
                         sli_295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = f_16 * ski_177[k]
                   + f_3 * pc_y[k] * sli_289[k];

        t_374[k] = f_16 * ski_294[k]
                   + f_8 * slh0_224[k]
                   - f_9 * slh1_224[k]
                   + f_3 * pc_x[k] * sli_294[k];

        t_375[k] = f_16 * ski_295[k]
                   + f_10 * slh0_225[k]
                   - f_11 * slh1_225[k]
                   + f_3 * pc_x[k] * sli_295[k];
    }

#pragma omp simd aligned(t_376, t_377, t_378, pc_x, pc_z, ski_297, ski_298, slh0_227, \
                         slh0_228, slh1_227, slh1_228, sli_290, sli_297, \
                         sli_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_376[k] = f_3 * pc_z[k] * sli_290[k];

        t_377[k] = f_16 * ski_297[k]
                   + f_10 * slh0_227[k]
                   - f_11 * slh1_227[k]
                   + f_3 * pc_x[k] * sli_297[k];

        t_378[k] = f_16 * ski_298[k]
                   + f_10 * slh0_228[k]
                   - f_11 * slh1_228[k]
                   + f_3 * pc_x[k] * sli_298[k];
    }

#pragma omp simd aligned(t_379, t_380, t_381, t_382, pc_x, pc_y, ski_182, ski_300, ski_301, \
                         ski_302, slh0_230, slh1_230, sli_294, sli_300, sli_301, \
                         sli_302 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_379[k] = f_16 * ski_182[k]
                   + f_3 * pc_y[k] * sli_294[k];

        t_380[k] = f_16 * ski_300[k]
                   + f_10 * slh0_230[k]
                   - f_11 * slh1_230[k]
                   + f_3 * pc_x[k] * sli_300[k];

        t_381[k] = f_16 * ski_301[k]
                   + f_3 * pc_x[k] * sli_301[k];

        t_382[k] = f_16 * ski_302[k]
                   + f_3 * pc_x[k] * sli_302[k];
    }

#pragma omp simd aligned(t_383, t_384, t_385, t_386, t_387, pc_x, ski_303, ski_304, ski_305, \
                         ski_306, ski_307, sli_303, sli_304, sli_305, sli_306, \
                         sli_307 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_383[k] = f_16 * ski_303[k]
                   + f_3 * pc_x[k] * sli_303[k];

        t_384[k] = f_16 * ski_304[k]
                   + f_3 * pc_x[k] * sli_304[k];

        t_385[k] = f_16 * ski_305[k]
                   + f_3 * pc_x[k] * sli_305[k];

        t_386[k] = f_16 * ski_306[k]
                   + f_3 * pc_x[k] * sli_306[k];

        t_387[k] = f_16 * ski_307[k]
                   + f_3 * pc_x[k] * sli_307[k];
    }

#pragma omp simd aligned(t_388, t_389, t_390, pc_y, pc_z, ski_189, ski_191, slh0_225, \
                         slh0_227, slh1_225, slh1_227, sli_301, \
                         sli_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_388[k] = f_16 * ski_189[k]
                   + f_1 * slh0_225[k]
                   - f_2 * slh1_225[k]
                   + f_3 * pc_y[k] * sli_301[k];

        t_389[k] = f_3 * pc_z[k] * sli_301[k];

        t_390[k] = f_16 * ski_191[k]
                   + f_4 * slh0_227[k]
                   - f_5 * slh1_227[k]
                   + f_3 * pc_y[k] * sli_303[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, pc_y, ski_192, ski_193, ski_194, slh0_228, \
                         slh0_229, slh0_230, slh1_228, slh1_229, slh1_230, sli_304, sli_305, \
                         sli_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = f_16 * ski_192[k]
                   + f_6 * slh0_228[k]
                   - f_7 * slh1_228[k]
                   + f_3 * pc_y[k] * sli_304[k];

        t_392[k] = f_16 * ski_193[k]
                   + f_8 * slh0_229[k]
                   - f_9 * slh1_229[k]
                   + f_3 * pc_y[k] * sli_305[k];

        t_393[k] = f_16 * ski_194[k]
                   + f_10 * slh0_230[k]
                   - f_11 * slh1_230[k]
                   + f_3 * pc_y[k] * sli_306[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, t_397, pb_z, pc_y, pc_z, skk0_216, ski_195, \
                         ski_196, skk1_216, slh0_230, slh1_230, sli_307, \
                         sli_308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = f_16 * ski_195[k]
                   + f_3 * pc_y[k] * sli_307[k];

        t_395[k] = f_1 * slh0_230[k]
                   - f_2 * slh1_230[k]
                   + f_3 * pc_z[k] * sli_307[k];

        t_396[k] = pb_z[k] * skk0_216[k]
                   - f_12 * pc_z[k] * skk1_216[k];

        t_397[k] = f_15 * ski_196[k]
                   + f_3 * pc_y[k] * sli_308[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, pb_z, pc_y, pc_z, skk0_219, ski_168, ski_198, \
                         skk1_219, sli_308, sli_310 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = f_13 * ski_168[k]
                   + f_3 * pc_z[k] * sli_308[k];

        t_399[k] = pb_z[k] * skk0_219[k]
                   - f_12 * pc_z[k] * skk1_219[k];

        t_400[k] = f_15 * ski_198[k]
                   + f_3 * pc_y[k] * sli_310[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, pb_z, pc_x, pc_z, skk0_222, ski_171, ski_313, \
                         skk1_222, slh0_236, slh1_236, sli_311, \
                         sli_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_16 * ski_313[k]
                   + f_4 * slh0_236[k]
                   - f_5 * slh1_236[k]
                   + f_3 * pc_x[k] * sli_313[k];

        t_402[k] = pb_z[k] * skk0_222[k]
                   - f_12 * pc_z[k] * skk1_222[k];

        t_403[k] = f_13 * ski_171[k]
                   + f_3 * pc_z[k] * sli_311[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, pb_z, pc_x, pc_y, pc_z, skk0_226, ski_201, \
                         ski_317, skk1_226, slh0_240, slh1_240, sli_313, \
                         sli_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = f_15 * ski_201[k]
                   + f_3 * pc_y[k] * sli_313[k];

        t_405[k] = f_16 * ski_317[k]
                   + f_6 * slh0_240[k]
                   - f_7 * slh1_240[k]
                   + f_3 * pc_x[k] * sli_317[k];

        t_406[k] = pb_z[k] * skk0_226[k]
                   - f_12 * pc_z[k] * skk1_226[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, pb_z, pc_y, pc_z, skk0_228, ski_174, ski_175, \
                         ski_205, skk1_228, sli_314, sli_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = f_13 * ski_174[k]
                   + f_3 * pc_z[k] * sli_314[k];

        t_408[k] = pb_z[k] * skk0_228[k]
                   + f_14 * ski_175[k]
                   - f_12 * pc_z[k] * skk1_228[k];

        t_409[k] = f_15 * ski_205[k]
                   + f_3 * pc_y[k] * sli_317[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, pb_z, pc_x, pc_z, skk0_231, ski_178, ski_322, \
                         skk1_231, slh0_245, slh1_245, sli_318, \
                         sli_322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = f_16 * ski_322[k]
                   + f_8 * slh0_245[k]
                   - f_9 * slh1_245[k]
                   + f_3 * pc_x[k] * sli_322[k];

        t_411[k] = pb_z[k] * skk0_231[k]
                   - f_12 * pc_z[k] * skk1_231[k];

        t_412[k] = f_13 * ski_178[k]
                   + f_3 * pc_z[k] * sli_318[k];
    }

#pragma omp simd aligned(t_413, t_414, t_415, pb_z, pc_y, pc_z, skk0_233, skk0_234, ski_179, \
                         ski_180, ski_210, skk1_233, skk1_234, \
                         sli_322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_413[k] = pb_z[k] * skk0_233[k]
                   + f_14 * ski_179[k]
                   - f_12 * pc_z[k] * skk1_233[k];

        t_414[k] = pb_z[k] * skk0_234[k]
                   + f_15 * ski_180[k]
                   - f_12 * pc_z[k] * skk1_234[k];

        t_415[k] = f_15 * ski_210[k]
                   + f_3 * pc_y[k] * sli_322[k];
    }

#pragma omp simd aligned(t_416, t_417, t_418, t_419, pc_x, ski_328, ski_329, ski_330, ski_331, \
                         slh0_251, slh1_251, sli_328, sli_329, sli_330, \
                         sli_331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_416[k] = f_16 * ski_328[k]
                   + f_10 * slh0_251[k]
                   - f_11 * slh1_251[k]
                   + f_3 * pc_x[k] * sli_328[k];

        t_417[k] = f_16 * ski_329[k]
                   + f_3 * pc_x[k] * sli_329[k];

        t_418[k] = f_16 * ski_330[k]
                   + f_3 * pc_x[k] * sli_330[k];

        t_419[k] = f_16 * ski_331[k]
                   + f_3 * pc_x[k] * sli_331[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, pc_x, ski_332, ski_333, ski_334, ski_335, \
                         sli_332, sli_333, sli_334, sli_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_16 * ski_332[k]
                   + f_3 * pc_x[k] * sli_332[k];

        t_421[k] = f_16 * ski_333[k]
                   + f_3 * pc_x[k] * sli_333[k];

        t_422[k] = f_16 * ski_334[k]
                   + f_3 * pc_x[k] * sli_334[k];

        t_423[k] = f_16 * ski_335[k]
                   + f_3 * pc_x[k] * sli_335[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, pb_z, pc_y, pc_z, skk0_244, ski_189, ski_219, \
                         skk1_244, slh0_248, slh1_248, sli_329, \
                         sli_331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = pb_z[k] * skk0_244[k]
                   - f_12 * pc_z[k] * skk1_244[k];

        t_425[k] = f_13 * ski_189[k]
                   + f_3 * pc_z[k] * sli_329[k];

        t_426[k] = f_15 * ski_219[k]
                   + f_4 * slh0_248[k]
                   - f_5 * slh1_248[k]
                   + f_3 * pc_y[k] * sli_331[k];
    }

#pragma omp simd aligned(t_427, t_428, t_429, pc_y, ski_220, ski_221, ski_222, slh0_249, \
                         slh0_250, slh0_251, slh1_249, slh1_250, slh1_251, sli_332, sli_333, \
                         sli_334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_427[k] = f_15 * ski_220[k]
                   + f_6 * slh0_249[k]
                   - f_7 * slh1_249[k]
                   + f_3 * pc_y[k] * sli_332[k];

        t_428[k] = f_15 * ski_221[k]
                   + f_8 * slh0_250[k]
                   - f_9 * slh1_250[k]
                   + f_3 * pc_y[k] * sli_333[k];

        t_429[k] = f_15 * ski_222[k]
                   + f_10 * slh0_251[k]
                   - f_11 * slh1_251[k]
                   + f_3 * pc_y[k] * sli_334[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, pc_x, pc_y, pc_z, ski_195, ski_223, ski_336, \
                         slh0_251, slh0_252, slh1_251, slh1_252, sli_335, \
                         sli_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = f_15 * ski_223[k]
                   + f_3 * pc_y[k] * sli_335[k];

        t_431[k] = f_13 * ski_195[k]
                   + f_1 * slh0_251[k]
                   - f_2 * slh1_251[k]
                   + f_3 * pc_z[k] * sli_335[k];

        t_432[k] = f_16 * ski_336[k]
                   + f_1 * slh0_252[k]
                   - f_2 * slh1_252[k]
                   + f_3 * pc_x[k] * sli_336[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, t_436, pc_x, pc_y, pc_z, ski_196, ski_224, \
                         ski_226, ski_339, slh0_255, slh1_255, sli_336, sli_338, \
                         sli_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = f_14 * ski_224[k]
                   + f_3 * pc_y[k] * sli_336[k];

        t_434[k] = f_14 * ski_196[k]
                   + f_3 * pc_z[k] * sli_336[k];

        t_435[k] = f_16 * ski_339[k]
                   + f_4 * slh0_255[k]
                   - f_5 * slh1_255[k]
                   + f_3 * pc_x[k] * sli_339[k];

        t_436[k] = f_14 * ski_226[k]
                   + f_3 * pc_y[k] * sli_338[k];
    }

#pragma omp simd aligned(t_437, t_438, t_439, pc_x, pc_z, ski_199, ski_341, ski_342, slh0_257, \
                         slh0_258, slh1_257, slh1_258, sli_339, sli_341, \
                         sli_342 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_437[k] = f_16 * ski_341[k]
                   + f_4 * slh0_257[k]
                   - f_5 * slh1_257[k]
                   + f_3 * pc_x[k] * sli_341[k];

        t_438[k] = f_16 * ski_342[k]
                   + f_6 * slh0_258[k]
                   - f_7 * slh1_258[k]
                   + f_3 * pc_x[k] * sli_342[k];

        t_439[k] = f_14 * ski_199[k]
                   + f_3 * pc_z[k] * sli_339[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, pc_x, pc_y, ski_229, ski_345, ski_346, slh0_261, \
                         slh0_262, slh1_261, slh1_262, sli_341, sli_345, \
                         sli_346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = f_14 * ski_229[k]
                   + f_3 * pc_y[k] * sli_341[k];

        t_441[k] = f_16 * ski_345[k]
                   + f_6 * slh0_261[k]
                   - f_7 * slh1_261[k]
                   + f_3 * pc_x[k] * sli_345[k];

        t_442[k] = f_16 * ski_346[k]
                   + f_8 * slh0_262[k]
                   - f_9 * slh1_262[k]
                   + f_3 * pc_x[k] * sli_346[k];
    }

#pragma omp simd aligned(t_443, t_444, t_445, pc_x, pc_y, pc_z, ski_202, ski_233, ski_348, \
                         slh0_264, slh1_264, sli_342, sli_345, \
                         sli_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_443[k] = f_14 * ski_202[k]
                   + f_3 * pc_z[k] * sli_342[k];

        t_444[k] = f_16 * ski_348[k]
                   + f_8 * slh0_264[k]
                   - f_9 * slh1_264[k]
                   + f_3 * pc_x[k] * sli_348[k];

        t_445[k] = f_14 * ski_233[k]
                   + f_3 * pc_y[k] * sli_345[k];
    }

#pragma omp simd aligned(t_446, t_447, t_448, pc_x, pc_z, ski_206, ski_350, ski_351, slh0_266, \
                         slh0_267, slh1_266, slh1_267, sli_346, sli_350, \
                         sli_351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_446[k] = f_16 * ski_350[k]
                   + f_8 * slh0_266[k]
                   - f_9 * slh1_266[k]
                   + f_3 * pc_x[k] * sli_350[k];

        t_447[k] = f_16 * ski_351[k]
                   + f_10 * slh0_267[k]
                   - f_11 * slh1_267[k]
                   + f_3 * pc_x[k] * sli_351[k];

        t_448[k] = f_14 * ski_206[k]
                   + f_3 * pc_z[k] * sli_346[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, pc_x, pc_y, ski_238, ski_353, ski_354, slh0_269, \
                         slh0_270, slh1_269, slh1_270, sli_350, sli_353, \
                         sli_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_16 * ski_353[k]
                   + f_10 * slh0_269[k]
                   - f_11 * slh1_269[k]
                   + f_3 * pc_x[k] * sli_353[k];

        t_450[k] = f_16 * ski_354[k]
                   + f_10 * slh0_270[k]
                   - f_11 * slh1_270[k]
                   + f_3 * pc_x[k] * sli_354[k];

        t_451[k] = f_14 * ski_238[k]
                   + f_3 * pc_y[k] * sli_350[k];
    }

#pragma omp simd aligned(t_452, t_453, t_454, t_455, pc_x, ski_356, ski_357, ski_358, ski_359, \
                         slh0_272, slh1_272, sli_356, sli_357, sli_358, \
                         sli_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_452[k] = f_16 * ski_356[k]
                   + f_10 * slh0_272[k]
                   - f_11 * slh1_272[k]
                   + f_3 * pc_x[k] * sli_356[k];

        t_453[k] = f_16 * ski_357[k]
                   + f_3 * pc_x[k] * sli_357[k];

        t_454[k] = f_16 * ski_358[k]
                   + f_3 * pc_x[k] * sli_358[k];

        t_455[k] = f_16 * ski_359[k]
                   + f_3 * pc_x[k] * sli_359[k];
    }
}

static auto
compute_prim_slk_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t skk0,
                                                          const size_t ski, const size_t skk1,
                                                          const size_t slh0, const size_t slh1,
                                                          const size_t sli, const size_t ncols,
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

    const auto *skk0_324 = buffer.data(skk0 + 324);
    const auto *skk0_327 = buffer.data(skk0 + 327);
    const auto *skk0_329 = buffer.data(skk0 + 329);
    const auto *skk0_330 = buffer.data(skk0 + 330);
    const auto *skk0_333 = buffer.data(skk0 + 333);
    const auto *skk0_334 = buffer.data(skk0 + 334);
    const auto *skk0_336 = buffer.data(skk0 + 336);
    const auto *skk0_338 = buffer.data(skk0 + 338);
    const auto *skk0_339 = buffer.data(skk0 + 339);
    const auto *skk0_341 = buffer.data(skk0 + 341);
    const auto *skk0_342 = buffer.data(skk0 + 342);
    const auto *skk0_344 = buffer.data(skk0 + 344);
    const auto *skk0_359 = buffer.data(skk0 + 359);

    const auto *ski_217 = buffer.data(ski + 217);
    const auto *ski_223 = buffer.data(ski + 223);
    const auto *ski_224 = buffer.data(ski + 224);
    const auto *ski_227 = buffer.data(ski + 227);
    const auto *ski_230 = buffer.data(ski + 230);
    const auto *ski_234 = buffer.data(ski + 234);
    const auto *ski_245 = buffer.data(ski + 245);
    const auto *ski_247 = buffer.data(ski + 247);
    const auto *ski_248 = buffer.data(ski + 248);
    const auto *ski_249 = buffer.data(ski + 249);
    const auto *ski_250 = buffer.data(ski + 250);
    const auto *ski_251 = buffer.data(ski + 251);
    const auto *ski_252 = buffer.data(ski + 252);
    const auto *ski_253 = buffer.data(ski + 253);
    const auto *ski_254 = buffer.data(ski + 254);
    const auto *ski_255 = buffer.data(ski + 255);
    const auto *ski_257 = buffer.data(ski + 257);
    const auto *ski_258 = buffer.data(ski + 258);
    const auto *ski_260 = buffer.data(ski + 260);
    const auto *ski_261 = buffer.data(ski + 261);
    const auto *ski_262 = buffer.data(ski + 262);
    const auto *ski_264 = buffer.data(ski + 264);
    const auto *ski_265 = buffer.data(ski + 265);
    const auto *ski_266 = buffer.data(ski + 266);
    const auto *ski_273 = buffer.data(ski + 273);
    const auto *ski_275 = buffer.data(ski + 275);
    const auto *ski_276 = buffer.data(ski + 276);
    const auto *ski_277 = buffer.data(ski + 277);
    const auto *ski_278 = buffer.data(ski + 278);
    const auto *ski_279 = buffer.data(ski + 279);
    const auto *ski_280 = buffer.data(ski + 280);
    const auto *ski_282 = buffer.data(ski + 282);
    const auto *ski_285 = buffer.data(ski + 285);
    const auto *ski_289 = buffer.data(ski + 289);
    const auto *ski_294 = buffer.data(ski + 294);
    const auto *ski_301 = buffer.data(ski + 301);
    const auto *ski_303 = buffer.data(ski + 303);
    const auto *ski_360 = buffer.data(ski + 360);
    const auto *ski_361 = buffer.data(ski + 361);
    const auto *ski_362 = buffer.data(ski + 362);
    const auto *ski_363 = buffer.data(ski + 363);
    const auto *ski_385 = buffer.data(ski + 385);
    const auto *ski_386 = buffer.data(ski + 386);
    const auto *ski_387 = buffer.data(ski + 387);
    const auto *ski_388 = buffer.data(ski + 388);
    const auto *ski_389 = buffer.data(ski + 389);
    const auto *ski_390 = buffer.data(ski + 390);
    const auto *ski_391 = buffer.data(ski + 391);
    const auto *ski_392 = buffer.data(ski + 392);
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

    const auto *skk1_324 = buffer.data(skk1 + 324);
    const auto *skk1_327 = buffer.data(skk1 + 327);
    const auto *skk1_329 = buffer.data(skk1 + 329);
    const auto *skk1_330 = buffer.data(skk1 + 330);
    const auto *skk1_333 = buffer.data(skk1 + 333);
    const auto *skk1_334 = buffer.data(skk1 + 334);
    const auto *skk1_336 = buffer.data(skk1 + 336);
    const auto *skk1_338 = buffer.data(skk1 + 338);
    const auto *skk1_339 = buffer.data(skk1 + 339);
    const auto *skk1_341 = buffer.data(skk1 + 341);
    const auto *skk1_342 = buffer.data(skk1 + 342);
    const auto *skk1_344 = buffer.data(skk1 + 344);
    const auto *skk1_359 = buffer.data(skk1 + 359);

    const auto *slh0_267 = buffer.data(slh0 + 267);
    const auto *slh0_269 = buffer.data(slh0 + 269);
    const auto *slh0_270 = buffer.data(slh0 + 270);
    const auto *slh0_271 = buffer.data(slh0 + 271);
    const auto *slh0_272 = buffer.data(slh0 + 272);
    const auto *slh0_288 = buffer.data(slh0 + 288);
    const auto *slh0_290 = buffer.data(slh0 + 290);
    const auto *slh0_291 = buffer.data(slh0 + 291);
    const auto *slh0_292 = buffer.data(slh0 + 292);
    const auto *slh0_293 = buffer.data(slh0 + 293);
    const auto *slh0_294 = buffer.data(slh0 + 294);
    const auto *slh0_297 = buffer.data(slh0 + 297);
    const auto *slh0_299 = buffer.data(slh0 + 299);
    const auto *slh0_300 = buffer.data(slh0 + 300);
    const auto *slh0_303 = buffer.data(slh0 + 303);
    const auto *slh0_304 = buffer.data(slh0 + 304);
    const auto *slh0_306 = buffer.data(slh0 + 306);
    const auto *slh0_308 = buffer.data(slh0 + 308);
    const auto *slh0_309 = buffer.data(slh0 + 309);
    const auto *slh0_311 = buffer.data(slh0 + 311);
    const auto *slh0_312 = buffer.data(slh0 + 312);
    const auto *slh0_313 = buffer.data(slh0 + 313);
    const auto *slh0_314 = buffer.data(slh0 + 314);
    const auto *slh0_315 = buffer.data(slh0 + 315);
    const auto *slh0_318 = buffer.data(slh0 + 318);
    const auto *slh0_320 = buffer.data(slh0 + 320);
    const auto *slh0_321 = buffer.data(slh0 + 321);
    const auto *slh0_324 = buffer.data(slh0 + 324);
    const auto *slh0_325 = buffer.data(slh0 + 325);
    const auto *slh0_327 = buffer.data(slh0 + 327);
    const auto *slh0_329 = buffer.data(slh0 + 329);
    const auto *slh0_330 = buffer.data(slh0 + 330);
    const auto *slh0_332 = buffer.data(slh0 + 332);
    const auto *slh0_333 = buffer.data(slh0 + 333);
    const auto *slh0_335 = buffer.data(slh0 + 335);

    const auto *slh1_267 = buffer.data(slh1 + 267);
    const auto *slh1_269 = buffer.data(slh1 + 269);
    const auto *slh1_270 = buffer.data(slh1 + 270);
    const auto *slh1_271 = buffer.data(slh1 + 271);
    const auto *slh1_272 = buffer.data(slh1 + 272);
    const auto *slh1_288 = buffer.data(slh1 + 288);
    const auto *slh1_290 = buffer.data(slh1 + 290);
    const auto *slh1_291 = buffer.data(slh1 + 291);
    const auto *slh1_292 = buffer.data(slh1 + 292);
    const auto *slh1_293 = buffer.data(slh1 + 293);
    const auto *slh1_294 = buffer.data(slh1 + 294);
    const auto *slh1_297 = buffer.data(slh1 + 297);
    const auto *slh1_299 = buffer.data(slh1 + 299);
    const auto *slh1_300 = buffer.data(slh1 + 300);
    const auto *slh1_303 = buffer.data(slh1 + 303);
    const auto *slh1_304 = buffer.data(slh1 + 304);
    const auto *slh1_306 = buffer.data(slh1 + 306);
    const auto *slh1_308 = buffer.data(slh1 + 308);
    const auto *slh1_309 = buffer.data(slh1 + 309);
    const auto *slh1_311 = buffer.data(slh1 + 311);
    const auto *slh1_312 = buffer.data(slh1 + 312);
    const auto *slh1_313 = buffer.data(slh1 + 313);
    const auto *slh1_314 = buffer.data(slh1 + 314);
    const auto *slh1_315 = buffer.data(slh1 + 315);
    const auto *slh1_318 = buffer.data(slh1 + 318);
    const auto *slh1_320 = buffer.data(slh1 + 320);
    const auto *slh1_321 = buffer.data(slh1 + 321);
    const auto *slh1_324 = buffer.data(slh1 + 324);
    const auto *slh1_325 = buffer.data(slh1 + 325);
    const auto *slh1_327 = buffer.data(slh1 + 327);
    const auto *slh1_329 = buffer.data(slh1 + 329);
    const auto *slh1_330 = buffer.data(slh1 + 330);
    const auto *slh1_332 = buffer.data(slh1 + 332);
    const auto *slh1_333 = buffer.data(slh1 + 333);
    const auto *slh1_335 = buffer.data(slh1 + 335);

    const auto *sli_357 = buffer.data(sli + 357);
    const auto *sli_359 = buffer.data(sli + 359);
    const auto *sli_360 = buffer.data(sli + 360);
    const auto *sli_361 = buffer.data(sli + 361);
    const auto *sli_362 = buffer.data(sli + 362);
    const auto *sli_363 = buffer.data(sli + 363);
    const auto *sli_364 = buffer.data(sli + 364);
    const auto *sli_366 = buffer.data(sli + 366);
    const auto *sli_367 = buffer.data(sli + 367);
    const auto *sli_369 = buffer.data(sli + 369);
    const auto *sli_370 = buffer.data(sli + 370);
    const auto *sli_373 = buffer.data(sli + 373);
    const auto *sli_374 = buffer.data(sli + 374);
    const auto *sli_378 = buffer.data(sli + 378);
    const auto *sli_385 = buffer.data(sli + 385);
    const auto *sli_386 = buffer.data(sli + 386);
    const auto *sli_387 = buffer.data(sli + 387);
    const auto *sli_388 = buffer.data(sli + 388);
    const auto *sli_389 = buffer.data(sli + 389);
    const auto *sli_390 = buffer.data(sli + 390);
    const auto *sli_391 = buffer.data(sli + 391);
    const auto *sli_392 = buffer.data(sli + 392);
    const auto *sli_394 = buffer.data(sli + 394);
    const auto *sli_395 = buffer.data(sli + 395);
    const auto *sli_397 = buffer.data(sli + 397);
    const auto *sli_398 = buffer.data(sli + 398);
    const auto *sli_401 = buffer.data(sli + 401);
    const auto *sli_402 = buffer.data(sli + 402);
    const auto *sli_404 = buffer.data(sli + 404);
    const auto *sli_406 = buffer.data(sli + 406);
    const auto *sli_407 = buffer.data(sli + 407);
    const auto *sli_409 = buffer.data(sli + 409);
    const auto *sli_410 = buffer.data(sli + 410);
    const auto *sli_412 = buffer.data(sli + 412);
    const auto *sli_413 = buffer.data(sli + 413);
    const auto *sli_414 = buffer.data(sli + 414);
    const auto *sli_415 = buffer.data(sli + 415);
    const auto *sli_416 = buffer.data(sli + 416);
    const auto *sli_417 = buffer.data(sli + 417);
    const auto *sli_418 = buffer.data(sli + 418);
    const auto *sli_419 = buffer.data(sli + 419);
    const auto *sli_420 = buffer.data(sli + 420);
    const auto *sli_422 = buffer.data(sli + 422);
    const auto *sli_423 = buffer.data(sli + 423);
    const auto *sli_425 = buffer.data(sli + 425);
    const auto *sli_426 = buffer.data(sli + 426);
    const auto *sli_429 = buffer.data(sli + 429);
    const auto *sli_430 = buffer.data(sli + 430);
    const auto *sli_432 = buffer.data(sli + 432);
    const auto *sli_434 = buffer.data(sli + 434);
    const auto *sli_435 = buffer.data(sli + 435);
    const auto *sli_437 = buffer.data(sli + 437);
    const auto *sli_438 = buffer.data(sli + 438);
    const auto *sli_440 = buffer.data(sli + 440);
    const auto *sli_441 = buffer.data(sli + 441);
    const auto *sli_442 = buffer.data(sli + 442);
    const auto *sli_443 = buffer.data(sli + 443);
    const auto *sli_444 = buffer.data(sli + 444);
    const auto *sli_445 = buffer.data(sli + 445);
    const auto *sli_446 = buffer.data(sli + 446);
    const auto *sli_447 = buffer.data(sli + 447);

#pragma omp simd aligned(t_456, t_457, t_458, t_459, pc_x, ski_360, ski_361, ski_362, ski_363, \
                         sli_360, sli_361, sli_362, sli_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_456[k] = f_16 * ski_360[k]
                   + f_3 * pc_x[k] * sli_360[k];

        t_457[k] = f_16 * ski_361[k]
                   + f_3 * pc_x[k] * sli_361[k];

        t_458[k] = f_16 * ski_362[k]
                   + f_3 * pc_x[k] * sli_362[k];

        t_459[k] = f_16 * ski_363[k]
                   + f_3 * pc_x[k] * sli_363[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, pc_y, pc_z, ski_217, ski_245, ski_247, slh0_267, \
                         slh0_269, slh1_267, slh1_269, sli_357, \
                         sli_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = f_14 * ski_245[k]
                   + f_1 * slh0_267[k]
                   - f_2 * slh1_267[k]
                   + f_3 * pc_y[k] * sli_357[k];

        t_461[k] = f_14 * ski_217[k]
                   + f_3 * pc_z[k] * sli_357[k];

        t_462[k] = f_14 * ski_247[k]
                   + f_4 * slh0_269[k]
                   - f_5 * slh1_269[k]
                   + f_3 * pc_y[k] * sli_359[k];
    }

#pragma omp simd aligned(t_463, t_464, t_465, pc_y, ski_248, ski_249, ski_250, slh0_270, \
                         slh0_271, slh0_272, slh1_270, slh1_271, slh1_272, sli_360, sli_361, \
                         sli_362 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_463[k] = f_14 * ski_248[k]
                   + f_6 * slh0_270[k]
                   - f_7 * slh1_270[k]
                   + f_3 * pc_y[k] * sli_360[k];

        t_464[k] = f_14 * ski_249[k]
                   + f_8 * slh0_271[k]
                   - f_9 * slh1_271[k]
                   + f_3 * pc_y[k] * sli_361[k];

        t_465[k] = f_14 * ski_250[k]
                   + f_10 * slh0_272[k]
                   - f_11 * slh1_272[k]
                   + f_3 * pc_y[k] * sli_362[k];
    }

#pragma omp simd aligned(t_466, t_467, t_468, t_469, pb_y, pc_y, pc_z, skk0_324, ski_223, \
                         ski_251, ski_252, skk1_324, slh0_272, slh1_272, sli_363, \
                         sli_364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_466[k] = f_14 * ski_251[k]
                   + f_3 * pc_y[k] * sli_363[k];

        t_467[k] = f_14 * ski_223[k]
                   + f_1 * slh0_272[k]
                   - f_2 * slh1_272[k]
                   + f_3 * pc_z[k] * sli_363[k];

        t_468[k] = pb_y[k] * skk0_324[k]
                   - f_12 * pc_y[k] * skk1_324[k];

        t_469[k] = f_13 * ski_252[k]
                   + f_3 * pc_y[k] * sli_364[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, pb_y, pc_y, pc_z, skk0_327, skk0_329, \
                         ski_224, ski_253, ski_254, skk1_327, skk1_329, sli_364, \
                         sli_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = f_15 * ski_224[k]
                   + f_3 * pc_z[k] * sli_364[k];

        t_471[k] = pb_y[k] * skk0_327[k]
                   + f_14 * ski_253[k]
                   - f_12 * pc_y[k] * skk1_327[k];

        t_472[k] = f_13 * ski_254[k]
                   + f_3 * pc_y[k] * sli_366[k];

        t_473[k] = pb_y[k] * skk0_329[k]
                   - f_12 * pc_y[k] * skk1_329[k];
    }

#pragma omp simd aligned(t_474, t_475, t_476, t_477, pb_y, pc_y, pc_z, skk0_330, skk0_333, \
                         ski_227, ski_255, ski_257, skk1_330, skk1_333, sli_367, \
                         sli_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = pb_y[k] * skk0_330[k]
                   + f_15 * ski_255[k]
                   - f_12 * pc_y[k] * skk1_330[k];

        t_475[k] = f_15 * ski_227[k]
                   + f_3 * pc_z[k] * sli_367[k];

        t_476[k] = f_13 * ski_257[k]
                   + f_3 * pc_y[k] * sli_369[k];

        t_477[k] = pb_y[k] * skk0_333[k]
                   - f_12 * pc_y[k] * skk1_333[k];
    }

#pragma omp simd aligned(t_478, t_479, t_480, pb_y, pc_y, pc_z, skk0_334, skk0_336, ski_230, \
                         ski_258, ski_260, skk1_334, skk1_336, \
                         sli_370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_478[k] = pb_y[k] * skk0_334[k]
                   + f_16 * ski_258[k]
                   - f_12 * pc_y[k] * skk1_334[k];

        t_479[k] = f_15 * ski_230[k]
                   + f_3 * pc_z[k] * sli_370[k];

        t_480[k] = pb_y[k] * skk0_336[k]
                   + f_14 * ski_260[k]
                   - f_12 * pc_y[k] * skk1_336[k];
    }

#pragma omp simd aligned(t_481, t_482, t_483, t_484, pb_y, pc_y, pc_z, skk0_338, skk0_339, \
                         ski_234, ski_261, ski_262, skk1_338, skk1_339, sli_373, \
                         sli_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_481[k] = f_13 * ski_261[k]
                   + f_3 * pc_y[k] * sli_373[k];

        t_482[k] = pb_y[k] * skk0_338[k]
                   - f_12 * pc_y[k] * skk1_338[k];

        t_483[k] = pb_y[k] * skk0_339[k]
                   + f_17 * ski_262[k]
                   - f_12 * pc_y[k] * skk1_339[k];

        t_484[k] = f_15 * ski_234[k]
                   + f_3 * pc_z[k] * sli_374[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, t_488, pb_y, pc_y, skk0_341, skk0_342, skk0_344, \
                         ski_264, ski_265, ski_266, skk1_341, skk1_342, skk1_344, \
                         sli_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = pb_y[k] * skk0_341[k]
                   + f_15 * ski_264[k]
                   - f_12 * pc_y[k] * skk1_341[k];

        t_486[k] = pb_y[k] * skk0_342[k]
                   + f_14 * ski_265[k]
                   - f_12 * pc_y[k] * skk1_342[k];

        t_487[k] = f_13 * ski_266[k]
                   + f_3 * pc_y[k] * sli_378[k];

        t_488[k] = pb_y[k] * skk0_344[k]
                   - f_12 * pc_y[k] * skk1_344[k];
    }

#pragma omp simd aligned(t_489, t_490, t_491, t_492, t_493, pc_x, ski_385, ski_386, ski_387, \
                         ski_388, ski_389, sli_385, sli_386, sli_387, sli_388, \
                         sli_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_489[k] = f_16 * ski_385[k]
                   + f_3 * pc_x[k] * sli_385[k];

        t_490[k] = f_16 * ski_386[k]
                   + f_3 * pc_x[k] * sli_386[k];

        t_491[k] = f_16 * ski_387[k]
                   + f_3 * pc_x[k] * sli_387[k];

        t_492[k] = f_16 * ski_388[k]
                   + f_3 * pc_x[k] * sli_388[k];

        t_493[k] = f_16 * ski_389[k]
                   + f_3 * pc_x[k] * sli_389[k];
    }

#pragma omp simd aligned(t_494, t_495, t_496, t_497, pc_x, pc_y, pc_z, ski_245, ski_273, \
                         ski_390, ski_391, slh0_288, slh1_288, sli_385, sli_390, \
                         sli_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_494[k] = f_16 * ski_390[k]
                   + f_3 * pc_x[k] * sli_390[k];

        t_495[k] = f_16 * ski_391[k]
                   + f_3 * pc_x[k] * sli_391[k];

        t_496[k] = f_13 * ski_273[k]
                   + f_1 * slh0_288[k]
                   - f_2 * slh1_288[k]
                   + f_3 * pc_y[k] * sli_385[k];

        t_497[k] = f_15 * ski_245[k]
                   + f_3 * pc_z[k] * sli_385[k];
    }

#pragma omp simd aligned(t_498, t_499, t_500, pc_y, ski_275, ski_276, ski_277, slh0_290, \
                         slh0_291, slh0_292, slh1_290, slh1_291, slh1_292, sli_387, sli_388, \
                         sli_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_498[k] = f_13 * ski_275[k]
                   + f_4 * slh0_290[k]
                   - f_5 * slh1_290[k]
                   + f_3 * pc_y[k] * sli_387[k];

        t_499[k] = f_13 * ski_276[k]
                   + f_6 * slh0_291[k]
                   - f_7 * slh1_291[k]
                   + f_3 * pc_y[k] * sli_388[k];

        t_500[k] = f_13 * ski_277[k]
                   + f_8 * slh0_292[k]
                   - f_9 * slh1_292[k]
                   + f_3 * pc_y[k] * sli_389[k];
    }

#pragma omp simd aligned(t_501, t_502, t_503, pb_y, pc_y, skk0_359, ski_278, ski_279, \
                         skk1_359, slh0_293, slh1_293, sli_390, \
                         sli_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_501[k] = f_13 * ski_278[k]
                   + f_10 * slh0_293[k]
                   - f_11 * slh1_293[k]
                   + f_3 * pc_y[k] * sli_390[k];

        t_502[k] = f_13 * ski_279[k]
                   + f_3 * pc_y[k] * sli_391[k];

        t_503[k] = pb_y[k] * skk0_359[k]
                   - f_12 * pc_y[k] * skk1_359[k];
    }

#pragma omp simd aligned(t_504, t_505, t_506, t_507, pc_x, pc_y, pc_z, ski_252, ski_392, \
                         ski_395, slh0_294, slh0_297, slh1_294, slh1_297, sli_392, \
                         sli_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_504[k] = f_16 * ski_392[k]
                   + f_1 * slh0_294[k]
                   - f_2 * slh1_294[k]
                   + f_3 * pc_x[k] * sli_392[k];

        t_505[k] = f_3 * pc_y[k] * sli_392[k];

        t_506[k] = f_16 * ski_252[k]
                   + f_3 * pc_z[k] * sli_392[k];

        t_507[k] = f_16 * ski_395[k]
                   + f_4 * slh0_297[k]
                   - f_5 * slh1_297[k]
                   + f_3 * pc_x[k] * sli_395[k];
    }

#pragma omp simd aligned(t_508, t_509, t_510, pc_x, pc_y, ski_397, ski_398, slh0_299, \
                         slh0_300, slh1_299, slh1_300, sli_394, sli_397, \
                         sli_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_508[k] = f_3 * pc_y[k] * sli_394[k];

        t_509[k] = f_16 * ski_397[k]
                   + f_4 * slh0_299[k]
                   - f_5 * slh1_299[k]
                   + f_3 * pc_x[k] * sli_397[k];

        t_510[k] = f_16 * ski_398[k]
                   + f_6 * slh0_300[k]
                   - f_7 * slh1_300[k]
                   + f_3 * pc_x[k] * sli_398[k];
    }

#pragma omp simd aligned(t_511, t_512, t_513, pc_x, pc_y, pc_z, ski_255, ski_401, slh0_303, \
                         slh1_303, sli_395, sli_397, sli_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_511[k] = f_16 * ski_255[k]
                   + f_3 * pc_z[k] * sli_395[k];

        t_512[k] = f_3 * pc_y[k] * sli_397[k];

        t_513[k] = f_16 * ski_401[k]
                   + f_6 * slh0_303[k]
                   - f_7 * slh1_303[k]
                   + f_3 * pc_x[k] * sli_401[k];
    }

#pragma omp simd aligned(t_514, t_515, t_516, pc_x, pc_z, ski_258, ski_402, ski_404, slh0_304, \
                         slh0_306, slh1_304, slh1_306, sli_398, sli_402, \
                         sli_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_514[k] = f_16 * ski_402[k]
                   + f_8 * slh0_304[k]
                   - f_9 * slh1_304[k]
                   + f_3 * pc_x[k] * sli_402[k];

        t_515[k] = f_16 * ski_258[k]
                   + f_3 * pc_z[k] * sli_398[k];

        t_516[k] = f_16 * ski_404[k]
                   + f_8 * slh0_306[k]
                   - f_9 * slh1_306[k]
                   + f_3 * pc_x[k] * sli_404[k];
    }

#pragma omp simd aligned(t_517, t_518, t_519, pc_x, pc_y, ski_406, ski_407, slh0_308, \
                         slh0_309, slh1_308, slh1_309, sli_401, sli_406, \
                         sli_407 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_517[k] = f_3 * pc_y[k] * sli_401[k];

        t_518[k] = f_16 * ski_406[k]
                   + f_8 * slh0_308[k]
                   - f_9 * slh1_308[k]
                   + f_3 * pc_x[k] * sli_406[k];

        t_519[k] = f_16 * ski_407[k]
                   + f_10 * slh0_309[k]
                   - f_11 * slh1_309[k]
                   + f_3 * pc_x[k] * sli_407[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, pc_x, pc_z, ski_262, ski_409, ski_410, slh0_311, \
                         slh0_312, slh1_311, slh1_312, sli_402, sli_409, \
                         sli_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = f_16 * ski_262[k]
                   + f_3 * pc_z[k] * sli_402[k];

        t_521[k] = f_16 * ski_409[k]
                   + f_10 * slh0_311[k]
                   - f_11 * slh1_311[k]
                   + f_3 * pc_x[k] * sli_409[k];

        t_522[k] = f_16 * ski_410[k]
                   + f_10 * slh0_312[k]
                   - f_11 * slh1_312[k]
                   + f_3 * pc_x[k] * sli_410[k];
    }

#pragma omp simd aligned(t_523, t_524, t_525, t_526, pc_x, pc_y, ski_412, ski_413, ski_414, \
                         slh0_314, slh1_314, sli_406, sli_412, sli_413, \
                         sli_414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_523[k] = f_3 * pc_y[k] * sli_406[k];

        t_524[k] = f_16 * ski_412[k]
                   + f_10 * slh0_314[k]
                   - f_11 * slh1_314[k]
                   + f_3 * pc_x[k] * sli_412[k];

        t_525[k] = f_16 * ski_413[k]
                   + f_3 * pc_x[k] * sli_413[k];

        t_526[k] = f_16 * ski_414[k]
                   + f_3 * pc_x[k] * sli_414[k];
    }

#pragma omp simd aligned(t_527, t_528, t_529, t_530, t_531, pc_x, ski_415, ski_416, ski_417, \
                         ski_418, ski_419, sli_415, sli_416, sli_417, sli_418, \
                         sli_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_527[k] = f_16 * ski_415[k]
                   + f_3 * pc_x[k] * sli_415[k];

        t_528[k] = f_16 * ski_416[k]
                   + f_3 * pc_x[k] * sli_416[k];

        t_529[k] = f_16 * ski_417[k]
                   + f_3 * pc_x[k] * sli_417[k];

        t_530[k] = f_16 * ski_418[k]
                   + f_3 * pc_x[k] * sli_418[k];

        t_531[k] = f_16 * ski_419[k]
                   + f_3 * pc_x[k] * sli_419[k];
    }

#pragma omp simd aligned(t_532, t_533, t_534, t_535, pc_y, pc_z, ski_273, slh0_309, slh0_311, \
                         slh0_312, slh1_309, slh1_311, slh1_312, sli_413, sli_415, \
                         sli_416 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_532[k] = f_1 * slh0_309[k]
                   - f_2 * slh1_309[k]
                   + f_3 * pc_y[k] * sli_413[k];

        t_533[k] = f_16 * ski_273[k]
                   + f_3 * pc_z[k] * sli_413[k];

        t_534[k] = f_4 * slh0_311[k]
                   - f_5 * slh1_311[k]
                   + f_3 * pc_y[k] * sli_415[k];

        t_535[k] = f_6 * slh0_312[k]
                   - f_7 * slh1_312[k]
                   + f_3 * pc_y[k] * sli_416[k];
    }

#pragma omp simd aligned(t_536, t_537, t_538, t_539, pc_y, pc_z, ski_279, slh0_313, slh0_314, \
                         slh1_313, slh1_314, sli_417, sli_418, \
                         sli_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_536[k] = f_8 * slh0_313[k]
                   - f_9 * slh1_313[k]
                   + f_3 * pc_y[k] * sli_417[k];

        t_537[k] = f_10 * slh0_314[k]
                   - f_11 * slh1_314[k]
                   + f_3 * pc_y[k] * sli_418[k];

        t_538[k] = f_3 * pc_y[k] * sli_419[k];

        t_539[k] = f_16 * ski_279[k]
                   + f_1 * slh0_314[k]
                   - f_2 * slh1_314[k]
                   + f_3 * pc_z[k] * sli_419[k];
    }

#pragma omp simd aligned(t_540, t_541, t_542, t_543, pc_x, pc_y, pc_z, ski_280, ski_420, \
                         ski_423, slh0_315, slh0_318, slh1_315, slh1_318, sli_420, \
                         sli_423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_540[k] = f_15 * ski_420[k]
                   + f_1 * slh0_315[k]
                   - f_2 * slh1_315[k]
                   + f_3 * pc_x[k] * sli_420[k];

        t_541[k] = f_17 * ski_280[k]
                   + f_3 * pc_y[k] * sli_420[k];

        t_542[k] = f_3 * pc_z[k] * sli_420[k];

        t_543[k] = f_15 * ski_423[k]
                   + f_4 * slh0_318[k]
                   - f_5 * slh1_318[k]
                   + f_3 * pc_x[k] * sli_423[k];
    }

#pragma omp simd aligned(t_544, t_545, t_546, pc_x, pc_y, ski_282, ski_425, ski_426, slh0_320, \
                         slh0_321, slh1_320, slh1_321, sli_422, sli_425, \
                         sli_426 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_544[k] = f_17 * ski_282[k]
                   + f_3 * pc_y[k] * sli_422[k];

        t_545[k] = f_15 * ski_425[k]
                   + f_4 * slh0_320[k]
                   - f_5 * slh1_320[k]
                   + f_3 * pc_x[k] * sli_425[k];

        t_546[k] = f_15 * ski_426[k]
                   + f_6 * slh0_321[k]
                   - f_7 * slh1_321[k]
                   + f_3 * pc_x[k] * sli_426[k];
    }

#pragma omp simd aligned(t_547, t_548, t_549, pc_x, pc_y, pc_z, ski_285, ski_429, slh0_324, \
                         slh1_324, sli_423, sli_425, sli_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_547[k] = f_3 * pc_z[k] * sli_423[k];

        t_548[k] = f_17 * ski_285[k]
                   + f_3 * pc_y[k] * sli_425[k];

        t_549[k] = f_15 * ski_429[k]
                   + f_6 * slh0_324[k]
                   - f_7 * slh1_324[k]
                   + f_3 * pc_x[k] * sli_429[k];
    }

#pragma omp simd aligned(t_550, t_551, t_552, pc_x, pc_z, ski_430, ski_432, slh0_325, \
                         slh0_327, slh1_325, slh1_327, sli_426, sli_430, \
                         sli_432 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_550[k] = f_15 * ski_430[k]
                   + f_8 * slh0_325[k]
                   - f_9 * slh1_325[k]
                   + f_3 * pc_x[k] * sli_430[k];

        t_551[k] = f_3 * pc_z[k] * sli_426[k];

        t_552[k] = f_15 * ski_432[k]
                   + f_8 * slh0_327[k]
                   - f_9 * slh1_327[k]
                   + f_3 * pc_x[k] * sli_432[k];
    }

#pragma omp simd aligned(t_553, t_554, t_555, pc_x, pc_y, ski_289, ski_434, ski_435, slh0_329, \
                         slh0_330, slh1_329, slh1_330, sli_429, sli_434, \
                         sli_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_553[k] = f_17 * ski_289[k]
                   + f_3 * pc_y[k] * sli_429[k];

        t_554[k] = f_15 * ski_434[k]
                   + f_8 * slh0_329[k]
                   - f_9 * slh1_329[k]
                   + f_3 * pc_x[k] * sli_434[k];

        t_555[k] = f_15 * ski_435[k]
                   + f_10 * slh0_330[k]
                   - f_11 * slh1_330[k]
                   + f_3 * pc_x[k] * sli_435[k];
    }

#pragma omp simd aligned(t_556, t_557, t_558, pc_x, pc_z, ski_437, ski_438, slh0_332, \
                         slh0_333, slh1_332, slh1_333, sli_430, sli_437, \
                         sli_438 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_556[k] = f_3 * pc_z[k] * sli_430[k];

        t_557[k] = f_15 * ski_437[k]
                   + f_10 * slh0_332[k]
                   - f_11 * slh1_332[k]
                   + f_3 * pc_x[k] * sli_437[k];

        t_558[k] = f_15 * ski_438[k]
                   + f_10 * slh0_333[k]
                   - f_11 * slh1_333[k]
                   + f_3 * pc_x[k] * sli_438[k];
    }

#pragma omp simd aligned(t_559, t_560, t_561, t_562, pc_x, pc_y, ski_294, ski_440, ski_441, \
                         ski_442, slh0_335, slh1_335, sli_434, sli_440, sli_441, \
                         sli_442 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_559[k] = f_17 * ski_294[k]
                   + f_3 * pc_y[k] * sli_434[k];

        t_560[k] = f_15 * ski_440[k]
                   + f_10 * slh0_335[k]
                   - f_11 * slh1_335[k]
                   + f_3 * pc_x[k] * sli_440[k];

        t_561[k] = f_15 * ski_441[k]
                   + f_3 * pc_x[k] * sli_441[k];

        t_562[k] = f_15 * ski_442[k]
                   + f_3 * pc_x[k] * sli_442[k];
    }

#pragma omp simd aligned(t_563, t_564, t_565, t_566, t_567, pc_x, ski_443, ski_444, ski_445, \
                         ski_446, ski_447, sli_443, sli_444, sli_445, sli_446, \
                         sli_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_563[k] = f_15 * ski_443[k]
                   + f_3 * pc_x[k] * sli_443[k];

        t_564[k] = f_15 * ski_444[k]
                   + f_3 * pc_x[k] * sli_444[k];

        t_565[k] = f_15 * ski_445[k]
                   + f_3 * pc_x[k] * sli_445[k];

        t_566[k] = f_15 * ski_446[k]
                   + f_3 * pc_x[k] * sli_446[k];

        t_567[k] = f_15 * ski_447[k]
                   + f_3 * pc_x[k] * sli_447[k];
    }

#pragma omp simd aligned(t_568, t_569, t_570, pc_y, pc_z, ski_301, ski_303, slh0_330, \
                         slh0_332, slh1_330, slh1_332, sli_441, \
                         sli_443 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_568[k] = f_17 * ski_301[k]
                   + f_1 * slh0_330[k]
                   - f_2 * slh1_330[k]
                   + f_3 * pc_y[k] * sli_441[k];

        t_569[k] = f_3 * pc_z[k] * sli_441[k];

        t_570[k] = f_17 * ski_303[k]
                   + f_4 * slh0_332[k]
                   - f_5 * slh1_332[k]
                   + f_3 * pc_y[k] * sli_443[k];
    }
}

static auto
compute_prim_slk_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t skk0,
                                                          const size_t ski, const size_t skk1,
                                                          const size_t slh0, const size_t slh1,
                                                          const size_t sli, const size_t ncols,
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

    const auto *skk0_360 = buffer.data(skk0 + 360);
    const auto *skk0_363 = buffer.data(skk0 + 363);
    const auto *skk0_366 = buffer.data(skk0 + 366);
    const auto *skk0_370 = buffer.data(skk0 + 370);
    const auto *skk0_372 = buffer.data(skk0 + 372);
    const auto *skk0_375 = buffer.data(skk0 + 375);
    const auto *skk0_377 = buffer.data(skk0 + 377);
    const auto *skk0_378 = buffer.data(skk0 + 378);
    const auto *skk0_388 = buffer.data(skk0 + 388);

    const auto *ski_280 = buffer.data(ski + 280);
    const auto *ski_283 = buffer.data(ski + 283);
    const auto *ski_286 = buffer.data(ski + 286);
    const auto *ski_287 = buffer.data(ski + 287);
    const auto *ski_290 = buffer.data(ski + 290);
    const auto *ski_291 = buffer.data(ski + 291);
    const auto *ski_292 = buffer.data(ski + 292);
    const auto *ski_301 = buffer.data(ski + 301);
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
    const auto *ski_329 = buffer.data(ski + 329);
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
    const auto *ski_350 = buffer.data(ski + 350);
    const auto *ski_357 = buffer.data(ski + 357);
    const auto *ski_359 = buffer.data(ski + 359);
    const auto *ski_360 = buffer.data(ski + 360);
    const auto *ski_361 = buffer.data(ski + 361);
    const auto *ski_362 = buffer.data(ski + 362);
    const auto *ski_363 = buffer.data(ski + 363);
    const auto *ski_364 = buffer.data(ski + 364);
    const auto *ski_366 = buffer.data(ski + 366);
    const auto *ski_369 = buffer.data(ski + 369);
    const auto *ski_373 = buffer.data(ski + 373);
    const auto *ski_378 = buffer.data(ski + 378);
    const auto *ski_385 = buffer.data(ski + 385);
    const auto *ski_387 = buffer.data(ski + 387);
    const auto *ski_453 = buffer.data(ski + 453);
    const auto *ski_457 = buffer.data(ski + 457);
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

    const auto *skk1_360 = buffer.data(skk1 + 360);
    const auto *skk1_363 = buffer.data(skk1 + 363);
    const auto *skk1_366 = buffer.data(skk1 + 366);
    const auto *skk1_370 = buffer.data(skk1 + 370);
    const auto *skk1_372 = buffer.data(skk1 + 372);
    const auto *skk1_375 = buffer.data(skk1 + 375);
    const auto *skk1_377 = buffer.data(skk1 + 377);
    const auto *skk1_378 = buffer.data(skk1 + 378);
    const auto *skk1_388 = buffer.data(skk1 + 388);

    const auto *slh0_333 = buffer.data(slh0 + 333);
    const auto *slh0_334 = buffer.data(slh0 + 334);
    const auto *slh0_335 = buffer.data(slh0 + 335);
    const auto *slh0_341 = buffer.data(slh0 + 341);
    const auto *slh0_345 = buffer.data(slh0 + 345);
    const auto *slh0_350 = buffer.data(slh0 + 350);
    const auto *slh0_353 = buffer.data(slh0 + 353);
    const auto *slh0_354 = buffer.data(slh0 + 354);
    const auto *slh0_355 = buffer.data(slh0 + 355);
    const auto *slh0_356 = buffer.data(slh0 + 356);
    const auto *slh0_357 = buffer.data(slh0 + 357);
    const auto *slh0_360 = buffer.data(slh0 + 360);
    const auto *slh0_362 = buffer.data(slh0 + 362);
    const auto *slh0_363 = buffer.data(slh0 + 363);
    const auto *slh0_366 = buffer.data(slh0 + 366);
    const auto *slh0_367 = buffer.data(slh0 + 367);
    const auto *slh0_369 = buffer.data(slh0 + 369);
    const auto *slh0_371 = buffer.data(slh0 + 371);
    const auto *slh0_372 = buffer.data(slh0 + 372);
    const auto *slh0_374 = buffer.data(slh0 + 374);
    const auto *slh0_375 = buffer.data(slh0 + 375);
    const auto *slh0_376 = buffer.data(slh0 + 376);
    const auto *slh0_377 = buffer.data(slh0 + 377);
    const auto *slh0_378 = buffer.data(slh0 + 378);
    const auto *slh0_381 = buffer.data(slh0 + 381);
    const auto *slh0_383 = buffer.data(slh0 + 383);
    const auto *slh0_384 = buffer.data(slh0 + 384);
    const auto *slh0_387 = buffer.data(slh0 + 387);
    const auto *slh0_388 = buffer.data(slh0 + 388);
    const auto *slh0_390 = buffer.data(slh0 + 390);
    const auto *slh0_392 = buffer.data(slh0 + 392);
    const auto *slh0_393 = buffer.data(slh0 + 393);
    const auto *slh0_395 = buffer.data(slh0 + 395);
    const auto *slh0_396 = buffer.data(slh0 + 396);
    const auto *slh0_398 = buffer.data(slh0 + 398);

    const auto *slh1_333 = buffer.data(slh1 + 333);
    const auto *slh1_334 = buffer.data(slh1 + 334);
    const auto *slh1_335 = buffer.data(slh1 + 335);
    const auto *slh1_341 = buffer.data(slh1 + 341);
    const auto *slh1_345 = buffer.data(slh1 + 345);
    const auto *slh1_350 = buffer.data(slh1 + 350);
    const auto *slh1_353 = buffer.data(slh1 + 353);
    const auto *slh1_354 = buffer.data(slh1 + 354);
    const auto *slh1_355 = buffer.data(slh1 + 355);
    const auto *slh1_356 = buffer.data(slh1 + 356);
    const auto *slh1_357 = buffer.data(slh1 + 357);
    const auto *slh1_360 = buffer.data(slh1 + 360);
    const auto *slh1_362 = buffer.data(slh1 + 362);
    const auto *slh1_363 = buffer.data(slh1 + 363);
    const auto *slh1_366 = buffer.data(slh1 + 366);
    const auto *slh1_367 = buffer.data(slh1 + 367);
    const auto *slh1_369 = buffer.data(slh1 + 369);
    const auto *slh1_371 = buffer.data(slh1 + 371);
    const auto *slh1_372 = buffer.data(slh1 + 372);
    const auto *slh1_374 = buffer.data(slh1 + 374);
    const auto *slh1_375 = buffer.data(slh1 + 375);
    const auto *slh1_376 = buffer.data(slh1 + 376);
    const auto *slh1_377 = buffer.data(slh1 + 377);
    const auto *slh1_378 = buffer.data(slh1 + 378);
    const auto *slh1_381 = buffer.data(slh1 + 381);
    const auto *slh1_383 = buffer.data(slh1 + 383);
    const auto *slh1_384 = buffer.data(slh1 + 384);
    const auto *slh1_387 = buffer.data(slh1 + 387);
    const auto *slh1_388 = buffer.data(slh1 + 388);
    const auto *slh1_390 = buffer.data(slh1 + 390);
    const auto *slh1_392 = buffer.data(slh1 + 392);
    const auto *slh1_393 = buffer.data(slh1 + 393);
    const auto *slh1_395 = buffer.data(slh1 + 395);
    const auto *slh1_396 = buffer.data(slh1 + 396);
    const auto *slh1_398 = buffer.data(slh1 + 398);

    const auto *sli_444 = buffer.data(sli + 444);
    const auto *sli_445 = buffer.data(sli + 445);
    const auto *sli_446 = buffer.data(sli + 446);
    const auto *sli_447 = buffer.data(sli + 447);
    const auto *sli_448 = buffer.data(sli + 448);
    const auto *sli_450 = buffer.data(sli + 450);
    const auto *sli_451 = buffer.data(sli + 451);
    const auto *sli_453 = buffer.data(sli + 453);
    const auto *sli_454 = buffer.data(sli + 454);
    const auto *sli_457 = buffer.data(sli + 457);
    const auto *sli_458 = buffer.data(sli + 458);
    const auto *sli_462 = buffer.data(sli + 462);
    const auto *sli_468 = buffer.data(sli + 468);
    const auto *sli_469 = buffer.data(sli + 469);
    const auto *sli_470 = buffer.data(sli + 470);
    const auto *sli_471 = buffer.data(sli + 471);
    const auto *sli_472 = buffer.data(sli + 472);
    const auto *sli_473 = buffer.data(sli + 473);
    const auto *sli_474 = buffer.data(sli + 474);
    const auto *sli_475 = buffer.data(sli + 475);
    const auto *sli_476 = buffer.data(sli + 476);
    const auto *sli_478 = buffer.data(sli + 478);
    const auto *sli_479 = buffer.data(sli + 479);
    const auto *sli_481 = buffer.data(sli + 481);
    const auto *sli_482 = buffer.data(sli + 482);
    const auto *sli_485 = buffer.data(sli + 485);
    const auto *sli_486 = buffer.data(sli + 486);
    const auto *sli_488 = buffer.data(sli + 488);
    const auto *sli_490 = buffer.data(sli + 490);
    const auto *sli_491 = buffer.data(sli + 491);
    const auto *sli_493 = buffer.data(sli + 493);
    const auto *sli_494 = buffer.data(sli + 494);
    const auto *sli_496 = buffer.data(sli + 496);
    const auto *sli_497 = buffer.data(sli + 497);
    const auto *sli_498 = buffer.data(sli + 498);
    const auto *sli_499 = buffer.data(sli + 499);
    const auto *sli_500 = buffer.data(sli + 500);
    const auto *sli_501 = buffer.data(sli + 501);
    const auto *sli_502 = buffer.data(sli + 502);
    const auto *sli_503 = buffer.data(sli + 503);
    const auto *sli_504 = buffer.data(sli + 504);
    const auto *sli_506 = buffer.data(sli + 506);
    const auto *sli_507 = buffer.data(sli + 507);
    const auto *sli_509 = buffer.data(sli + 509);
    const auto *sli_510 = buffer.data(sli + 510);
    const auto *sli_513 = buffer.data(sli + 513);
    const auto *sli_514 = buffer.data(sli + 514);
    const auto *sli_516 = buffer.data(sli + 516);
    const auto *sli_518 = buffer.data(sli + 518);
    const auto *sli_519 = buffer.data(sli + 519);
    const auto *sli_521 = buffer.data(sli + 521);
    const auto *sli_522 = buffer.data(sli + 522);
    const auto *sli_524 = buffer.data(sli + 524);
    const auto *sli_525 = buffer.data(sli + 525);
    const auto *sli_526 = buffer.data(sli + 526);
    const auto *sli_527 = buffer.data(sli + 527);
    const auto *sli_528 = buffer.data(sli + 528);
    const auto *sli_529 = buffer.data(sli + 529);
    const auto *sli_530 = buffer.data(sli + 530);
    const auto *sli_531 = buffer.data(sli + 531);

#pragma omp simd aligned(t_571, t_572, t_573, pc_y, ski_304, ski_305, ski_306, slh0_333, \
                         slh0_334, slh0_335, slh1_333, slh1_334, slh1_335, sli_444, sli_445, \
                         sli_446 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_571[k] = f_17 * ski_304[k]
                   + f_6 * slh0_333[k]
                   - f_7 * slh1_333[k]
                   + f_3 * pc_y[k] * sli_444[k];

        t_572[k] = f_17 * ski_305[k]
                   + f_8 * slh0_334[k]
                   - f_9 * slh1_334[k]
                   + f_3 * pc_y[k] * sli_445[k];

        t_573[k] = f_17 * ski_306[k]
                   + f_10 * slh0_335[k]
                   - f_11 * slh1_335[k]
                   + f_3 * pc_y[k] * sli_446[k];
    }

#pragma omp simd aligned(t_574, t_575, t_576, t_577, pb_z, pc_y, pc_z, skk0_360, ski_307, \
                         ski_308, skk1_360, slh0_335, slh1_335, sli_447, \
                         sli_448 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_574[k] = f_17 * ski_307[k]
                   + f_3 * pc_y[k] * sli_447[k];

        t_575[k] = f_1 * slh0_335[k]
                   - f_2 * slh1_335[k]
                   + f_3 * pc_z[k] * sli_447[k];

        t_576[k] = pb_z[k] * skk0_360[k]
                   - f_12 * pc_z[k] * skk1_360[k];

        t_577[k] = f_16 * ski_308[k]
                   + f_3 * pc_y[k] * sli_448[k];
    }

#pragma omp simd aligned(t_578, t_579, t_580, pb_z, pc_y, pc_z, skk0_363, ski_280, ski_310, \
                         skk1_363, sli_448, sli_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_578[k] = f_13 * ski_280[k]
                   + f_3 * pc_z[k] * sli_448[k];

        t_579[k] = pb_z[k] * skk0_363[k]
                   - f_12 * pc_z[k] * skk1_363[k];

        t_580[k] = f_16 * ski_310[k]
                   + f_3 * pc_y[k] * sli_450[k];
    }

#pragma omp simd aligned(t_581, t_582, t_583, pb_z, pc_x, pc_z, skk0_366, ski_283, ski_453, \
                         skk1_366, slh0_341, slh1_341, sli_451, \
                         sli_453 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_581[k] = f_15 * ski_453[k]
                   + f_4 * slh0_341[k]
                   - f_5 * slh1_341[k]
                   + f_3 * pc_x[k] * sli_453[k];

        t_582[k] = pb_z[k] * skk0_366[k]
                   - f_12 * pc_z[k] * skk1_366[k];

        t_583[k] = f_13 * ski_283[k]
                   + f_3 * pc_z[k] * sli_451[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, pb_z, pc_x, pc_y, pc_z, skk0_370, ski_313, \
                         ski_457, skk1_370, slh0_345, slh1_345, sli_453, \
                         sli_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = f_16 * ski_313[k]
                   + f_3 * pc_y[k] * sli_453[k];

        t_585[k] = f_15 * ski_457[k]
                   + f_6 * slh0_345[k]
                   - f_7 * slh1_345[k]
                   + f_3 * pc_x[k] * sli_457[k];

        t_586[k] = pb_z[k] * skk0_370[k]
                   - f_12 * pc_z[k] * skk1_370[k];
    }

#pragma omp simd aligned(t_587, t_588, t_589, pb_z, pc_y, pc_z, skk0_372, ski_286, ski_287, \
                         ski_317, skk1_372, sli_454, sli_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_587[k] = f_13 * ski_286[k]
                   + f_3 * pc_z[k] * sli_454[k];

        t_588[k] = pb_z[k] * skk0_372[k]
                   + f_14 * ski_287[k]
                   - f_12 * pc_z[k] * skk1_372[k];

        t_589[k] = f_16 * ski_317[k]
                   + f_3 * pc_y[k] * sli_457[k];
    }

#pragma omp simd aligned(t_590, t_591, t_592, pb_z, pc_x, pc_z, skk0_375, ski_290, ski_462, \
                         skk1_375, slh0_350, slh1_350, sli_458, \
                         sli_462 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_590[k] = f_15 * ski_462[k]
                   + f_8 * slh0_350[k]
                   - f_9 * slh1_350[k]
                   + f_3 * pc_x[k] * sli_462[k];

        t_591[k] = pb_z[k] * skk0_375[k]
                   - f_12 * pc_z[k] * skk1_375[k];

        t_592[k] = f_13 * ski_290[k]
                   + f_3 * pc_z[k] * sli_458[k];
    }

#pragma omp simd aligned(t_593, t_594, t_595, pb_z, pc_y, pc_z, skk0_377, skk0_378, ski_291, \
                         ski_292, ski_322, skk1_377, skk1_378, \
                         sli_462 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_593[k] = pb_z[k] * skk0_377[k]
                   + f_14 * ski_291[k]
                   - f_12 * pc_z[k] * skk1_377[k];

        t_594[k] = pb_z[k] * skk0_378[k]
                   + f_15 * ski_292[k]
                   - f_12 * pc_z[k] * skk1_378[k];

        t_595[k] = f_16 * ski_322[k]
                   + f_3 * pc_y[k] * sli_462[k];
    }

#pragma omp simd aligned(t_596, t_597, t_598, t_599, pc_x, ski_468, ski_469, ski_470, ski_471, \
                         slh0_356, slh1_356, sli_468, sli_469, sli_470, \
                         sli_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_596[k] = f_15 * ski_468[k]
                   + f_10 * slh0_356[k]
                   - f_11 * slh1_356[k]
                   + f_3 * pc_x[k] * sli_468[k];

        t_597[k] = f_15 * ski_469[k]
                   + f_3 * pc_x[k] * sli_469[k];

        t_598[k] = f_15 * ski_470[k]
                   + f_3 * pc_x[k] * sli_470[k];

        t_599[k] = f_15 * ski_471[k]
                   + f_3 * pc_x[k] * sli_471[k];
    }

#pragma omp simd aligned(t_600, t_601, t_602, t_603, pc_x, ski_472, ski_473, ski_474, ski_475, \
                         sli_472, sli_473, sli_474, sli_475 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = f_15 * ski_472[k]
                   + f_3 * pc_x[k] * sli_472[k];

        t_601[k] = f_15 * ski_473[k]
                   + f_3 * pc_x[k] * sli_473[k];

        t_602[k] = f_15 * ski_474[k]
                   + f_3 * pc_x[k] * sli_474[k];

        t_603[k] = f_15 * ski_475[k]
                   + f_3 * pc_x[k] * sli_475[k];
    }

#pragma omp simd aligned(t_604, t_605, t_606, pb_z, pc_y, pc_z, skk0_388, ski_301, ski_331, \
                         skk1_388, slh0_353, slh1_353, sli_469, \
                         sli_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_604[k] = pb_z[k] * skk0_388[k]
                   - f_12 * pc_z[k] * skk1_388[k];

        t_605[k] = f_13 * ski_301[k]
                   + f_3 * pc_z[k] * sli_469[k];

        t_606[k] = f_16 * ski_331[k]
                   + f_4 * slh0_353[k]
                   - f_5 * slh1_353[k]
                   + f_3 * pc_y[k] * sli_471[k];
    }

#pragma omp simd aligned(t_607, t_608, t_609, pc_y, ski_332, ski_333, ski_334, slh0_354, \
                         slh0_355, slh0_356, slh1_354, slh1_355, slh1_356, sli_472, sli_473, \
                         sli_474 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_607[k] = f_16 * ski_332[k]
                   + f_6 * slh0_354[k]
                   - f_7 * slh1_354[k]
                   + f_3 * pc_y[k] * sli_472[k];

        t_608[k] = f_16 * ski_333[k]
                   + f_8 * slh0_355[k]
                   - f_9 * slh1_355[k]
                   + f_3 * pc_y[k] * sli_473[k];

        t_609[k] = f_16 * ski_334[k]
                   + f_10 * slh0_356[k]
                   - f_11 * slh1_356[k]
                   + f_3 * pc_y[k] * sli_474[k];
    }

#pragma omp simd aligned(t_610, t_611, t_612, pc_x, pc_y, pc_z, ski_307, ski_335, ski_476, \
                         slh0_356, slh0_357, slh1_356, slh1_357, sli_475, \
                         sli_476 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_610[k] = f_16 * ski_335[k]
                   + f_3 * pc_y[k] * sli_475[k];

        t_611[k] = f_13 * ski_307[k]
                   + f_1 * slh0_356[k]
                   - f_2 * slh1_356[k]
                   + f_3 * pc_z[k] * sli_475[k];

        t_612[k] = f_15 * ski_476[k]
                   + f_1 * slh0_357[k]
                   - f_2 * slh1_357[k]
                   + f_3 * pc_x[k] * sli_476[k];
    }

#pragma omp simd aligned(t_613, t_614, t_615, t_616, pc_x, pc_y, pc_z, ski_308, ski_336, \
                         ski_338, ski_479, slh0_360, slh1_360, sli_476, sli_478, \
                         sli_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_613[k] = f_15 * ski_336[k]
                   + f_3 * pc_y[k] * sli_476[k];

        t_614[k] = f_14 * ski_308[k]
                   + f_3 * pc_z[k] * sli_476[k];

        t_615[k] = f_15 * ski_479[k]
                   + f_4 * slh0_360[k]
                   - f_5 * slh1_360[k]
                   + f_3 * pc_x[k] * sli_479[k];

        t_616[k] = f_15 * ski_338[k]
                   + f_3 * pc_y[k] * sli_478[k];
    }

#pragma omp simd aligned(t_617, t_618, t_619, pc_x, pc_z, ski_311, ski_481, ski_482, slh0_362, \
                         slh0_363, slh1_362, slh1_363, sli_479, sli_481, \
                         sli_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_617[k] = f_15 * ski_481[k]
                   + f_4 * slh0_362[k]
                   - f_5 * slh1_362[k]
                   + f_3 * pc_x[k] * sli_481[k];

        t_618[k] = f_15 * ski_482[k]
                   + f_6 * slh0_363[k]
                   - f_7 * slh1_363[k]
                   + f_3 * pc_x[k] * sli_482[k];

        t_619[k] = f_14 * ski_311[k]
                   + f_3 * pc_z[k] * sli_479[k];
    }

#pragma omp simd aligned(t_620, t_621, t_622, pc_x, pc_y, ski_341, ski_485, ski_486, slh0_366, \
                         slh0_367, slh1_366, slh1_367, sli_481, sli_485, \
                         sli_486 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_620[k] = f_15 * ski_341[k]
                   + f_3 * pc_y[k] * sli_481[k];

        t_621[k] = f_15 * ski_485[k]
                   + f_6 * slh0_366[k]
                   - f_7 * slh1_366[k]
                   + f_3 * pc_x[k] * sli_485[k];

        t_622[k] = f_15 * ski_486[k]
                   + f_8 * slh0_367[k]
                   - f_9 * slh1_367[k]
                   + f_3 * pc_x[k] * sli_486[k];
    }

#pragma omp simd aligned(t_623, t_624, t_625, pc_x, pc_y, pc_z, ski_314, ski_345, ski_488, \
                         slh0_369, slh1_369, sli_482, sli_485, \
                         sli_488 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_623[k] = f_14 * ski_314[k]
                   + f_3 * pc_z[k] * sli_482[k];

        t_624[k] = f_15 * ski_488[k]
                   + f_8 * slh0_369[k]
                   - f_9 * slh1_369[k]
                   + f_3 * pc_x[k] * sli_488[k];

        t_625[k] = f_15 * ski_345[k]
                   + f_3 * pc_y[k] * sli_485[k];
    }

#pragma omp simd aligned(t_626, t_627, t_628, pc_x, pc_z, ski_318, ski_490, ski_491, slh0_371, \
                         slh0_372, slh1_371, slh1_372, sli_486, sli_490, \
                         sli_491 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_626[k] = f_15 * ski_490[k]
                   + f_8 * slh0_371[k]
                   - f_9 * slh1_371[k]
                   + f_3 * pc_x[k] * sli_490[k];

        t_627[k] = f_15 * ski_491[k]
                   + f_10 * slh0_372[k]
                   - f_11 * slh1_372[k]
                   + f_3 * pc_x[k] * sli_491[k];

        t_628[k] = f_14 * ski_318[k]
                   + f_3 * pc_z[k] * sli_486[k];
    }

#pragma omp simd aligned(t_629, t_630, t_631, pc_x, pc_y, ski_350, ski_493, ski_494, slh0_374, \
                         slh0_375, slh1_374, slh1_375, sli_490, sli_493, \
                         sli_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_629[k] = f_15 * ski_493[k]
                   + f_10 * slh0_374[k]
                   - f_11 * slh1_374[k]
                   + f_3 * pc_x[k] * sli_493[k];

        t_630[k] = f_15 * ski_494[k]
                   + f_10 * slh0_375[k]
                   - f_11 * slh1_375[k]
                   + f_3 * pc_x[k] * sli_494[k];

        t_631[k] = f_15 * ski_350[k]
                   + f_3 * pc_y[k] * sli_490[k];
    }

#pragma omp simd aligned(t_632, t_633, t_634, t_635, pc_x, ski_496, ski_497, ski_498, ski_499, \
                         slh0_377, slh1_377, sli_496, sli_497, sli_498, \
                         sli_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_632[k] = f_15 * ski_496[k]
                   + f_10 * slh0_377[k]
                   - f_11 * slh1_377[k]
                   + f_3 * pc_x[k] * sli_496[k];

        t_633[k] = f_15 * ski_497[k]
                   + f_3 * pc_x[k] * sli_497[k];

        t_634[k] = f_15 * ski_498[k]
                   + f_3 * pc_x[k] * sli_498[k];

        t_635[k] = f_15 * ski_499[k]
                   + f_3 * pc_x[k] * sli_499[k];
    }

#pragma omp simd aligned(t_636, t_637, t_638, t_639, pc_x, ski_500, ski_501, ski_502, ski_503, \
                         sli_500, sli_501, sli_502, sli_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_636[k] = f_15 * ski_500[k]
                   + f_3 * pc_x[k] * sli_500[k];

        t_637[k] = f_15 * ski_501[k]
                   + f_3 * pc_x[k] * sli_501[k];

        t_638[k] = f_15 * ski_502[k]
                   + f_3 * pc_x[k] * sli_502[k];

        t_639[k] = f_15 * ski_503[k]
                   + f_3 * pc_x[k] * sli_503[k];
    }

#pragma omp simd aligned(t_640, t_641, t_642, pc_y, pc_z, ski_329, ski_357, ski_359, slh0_372, \
                         slh0_374, slh1_372, slh1_374, sli_497, \
                         sli_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_640[k] = f_15 * ski_357[k]
                   + f_1 * slh0_372[k]
                   - f_2 * slh1_372[k]
                   + f_3 * pc_y[k] * sli_497[k];

        t_641[k] = f_14 * ski_329[k]
                   + f_3 * pc_z[k] * sli_497[k];

        t_642[k] = f_15 * ski_359[k]
                   + f_4 * slh0_374[k]
                   - f_5 * slh1_374[k]
                   + f_3 * pc_y[k] * sli_499[k];
    }

#pragma omp simd aligned(t_643, t_644, t_645, pc_y, ski_360, ski_361, ski_362, slh0_375, \
                         slh0_376, slh0_377, slh1_375, slh1_376, slh1_377, sli_500, sli_501, \
                         sli_502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_643[k] = f_15 * ski_360[k]
                   + f_6 * slh0_375[k]
                   - f_7 * slh1_375[k]
                   + f_3 * pc_y[k] * sli_500[k];

        t_644[k] = f_15 * ski_361[k]
                   + f_8 * slh0_376[k]
                   - f_9 * slh1_376[k]
                   + f_3 * pc_y[k] * sli_501[k];

        t_645[k] = f_15 * ski_362[k]
                   + f_10 * slh0_377[k]
                   - f_11 * slh1_377[k]
                   + f_3 * pc_y[k] * sli_502[k];
    }

#pragma omp simd aligned(t_646, t_647, t_648, pc_x, pc_y, pc_z, ski_335, ski_363, ski_504, \
                         slh0_377, slh0_378, slh1_377, slh1_378, sli_503, \
                         sli_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_646[k] = f_15 * ski_363[k]
                   + f_3 * pc_y[k] * sli_503[k];

        t_647[k] = f_14 * ski_335[k]
                   + f_1 * slh0_377[k]
                   - f_2 * slh1_377[k]
                   + f_3 * pc_z[k] * sli_503[k];

        t_648[k] = f_15 * ski_504[k]
                   + f_1 * slh0_378[k]
                   - f_2 * slh1_378[k]
                   + f_3 * pc_x[k] * sli_504[k];
    }

#pragma omp simd aligned(t_649, t_650, t_651, t_652, pc_x, pc_y, pc_z, ski_336, ski_364, \
                         ski_366, ski_507, slh0_381, slh1_381, sli_504, sli_506, \
                         sli_507 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_649[k] = f_14 * ski_364[k]
                   + f_3 * pc_y[k] * sli_504[k];

        t_650[k] = f_15 * ski_336[k]
                   + f_3 * pc_z[k] * sli_504[k];

        t_651[k] = f_15 * ski_507[k]
                   + f_4 * slh0_381[k]
                   - f_5 * slh1_381[k]
                   + f_3 * pc_x[k] * sli_507[k];

        t_652[k] = f_14 * ski_366[k]
                   + f_3 * pc_y[k] * sli_506[k];
    }

#pragma omp simd aligned(t_653, t_654, t_655, pc_x, pc_z, ski_339, ski_509, ski_510, slh0_383, \
                         slh0_384, slh1_383, slh1_384, sli_507, sli_509, \
                         sli_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_653[k] = f_15 * ski_509[k]
                   + f_4 * slh0_383[k]
                   - f_5 * slh1_383[k]
                   + f_3 * pc_x[k] * sli_509[k];

        t_654[k] = f_15 * ski_510[k]
                   + f_6 * slh0_384[k]
                   - f_7 * slh1_384[k]
                   + f_3 * pc_x[k] * sli_510[k];

        t_655[k] = f_15 * ski_339[k]
                   + f_3 * pc_z[k] * sli_507[k];
    }

#pragma omp simd aligned(t_656, t_657, t_658, pc_x, pc_y, ski_369, ski_513, ski_514, slh0_387, \
                         slh0_388, slh1_387, slh1_388, sli_509, sli_513, \
                         sli_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_656[k] = f_14 * ski_369[k]
                   + f_3 * pc_y[k] * sli_509[k];

        t_657[k] = f_15 * ski_513[k]
                   + f_6 * slh0_387[k]
                   - f_7 * slh1_387[k]
                   + f_3 * pc_x[k] * sli_513[k];

        t_658[k] = f_15 * ski_514[k]
                   + f_8 * slh0_388[k]
                   - f_9 * slh1_388[k]
                   + f_3 * pc_x[k] * sli_514[k];
    }

#pragma omp simd aligned(t_659, t_660, t_661, pc_x, pc_y, pc_z, ski_342, ski_373, ski_516, \
                         slh0_390, slh1_390, sli_510, sli_513, \
                         sli_516 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_659[k] = f_15 * ski_342[k]
                   + f_3 * pc_z[k] * sli_510[k];

        t_660[k] = f_15 * ski_516[k]
                   + f_8 * slh0_390[k]
                   - f_9 * slh1_390[k]
                   + f_3 * pc_x[k] * sli_516[k];

        t_661[k] = f_14 * ski_373[k]
                   + f_3 * pc_y[k] * sli_513[k];
    }

#pragma omp simd aligned(t_662, t_663, t_664, pc_x, pc_z, ski_346, ski_518, ski_519, slh0_392, \
                         slh0_393, slh1_392, slh1_393, sli_514, sli_518, \
                         sli_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_662[k] = f_15 * ski_518[k]
                   + f_8 * slh0_392[k]
                   - f_9 * slh1_392[k]
                   + f_3 * pc_x[k] * sli_518[k];

        t_663[k] = f_15 * ski_519[k]
                   + f_10 * slh0_393[k]
                   - f_11 * slh1_393[k]
                   + f_3 * pc_x[k] * sli_519[k];

        t_664[k] = f_15 * ski_346[k]
                   + f_3 * pc_z[k] * sli_514[k];
    }

#pragma omp simd aligned(t_665, t_666, t_667, pc_x, pc_y, ski_378, ski_521, ski_522, slh0_395, \
                         slh0_396, slh1_395, slh1_396, sli_518, sli_521, \
                         sli_522 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_665[k] = f_15 * ski_521[k]
                   + f_10 * slh0_395[k]
                   - f_11 * slh1_395[k]
                   + f_3 * pc_x[k] * sli_521[k];

        t_666[k] = f_15 * ski_522[k]
                   + f_10 * slh0_396[k]
                   - f_11 * slh1_396[k]
                   + f_3 * pc_x[k] * sli_522[k];

        t_667[k] = f_14 * ski_378[k]
                   + f_3 * pc_y[k] * sli_518[k];
    }

#pragma omp simd aligned(t_668, t_669, t_670, t_671, pc_x, ski_524, ski_525, ski_526, ski_527, \
                         slh0_398, slh1_398, sli_524, sli_525, sli_526, \
                         sli_527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_668[k] = f_15 * ski_524[k]
                   + f_10 * slh0_398[k]
                   - f_11 * slh1_398[k]
                   + f_3 * pc_x[k] * sli_524[k];

        t_669[k] = f_15 * ski_525[k]
                   + f_3 * pc_x[k] * sli_525[k];

        t_670[k] = f_15 * ski_526[k]
                   + f_3 * pc_x[k] * sli_526[k];

        t_671[k] = f_15 * ski_527[k]
                   + f_3 * pc_x[k] * sli_527[k];
    }

#pragma omp simd aligned(t_672, t_673, t_674, t_675, pc_x, ski_528, ski_529, ski_530, ski_531, \
                         sli_528, sli_529, sli_530, sli_531 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_672[k] = f_15 * ski_528[k]
                   + f_3 * pc_x[k] * sli_528[k];

        t_673[k] = f_15 * ski_529[k]
                   + f_3 * pc_x[k] * sli_529[k];

        t_674[k] = f_15 * ski_530[k]
                   + f_3 * pc_x[k] * sli_530[k];

        t_675[k] = f_15 * ski_531[k]
                   + f_3 * pc_x[k] * sli_531[k];
    }

#pragma omp simd aligned(t_676, t_677, t_678, pc_y, pc_z, ski_357, ski_385, ski_387, slh0_393, \
                         slh0_395, slh1_393, slh1_395, sli_525, \
                         sli_527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_676[k] = f_14 * ski_385[k]
                   + f_1 * slh0_393[k]
                   - f_2 * slh1_393[k]
                   + f_3 * pc_y[k] * sli_525[k];

        t_677[k] = f_15 * ski_357[k]
                   + f_3 * pc_z[k] * sli_525[k];

        t_678[k] = f_14 * ski_387[k]
                   + f_4 * slh0_395[k]
                   - f_5 * slh1_395[k]
                   + f_3 * pc_y[k] * sli_527[k];
    }
}

static auto
compute_prim_slk_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t skk0,
                                                          const size_t ski, const size_t skk1,
                                                          const size_t slh0, const size_t slh1,
                                                          const size_t sli, const size_t ncols,
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
    const auto f_19 = 3.0 / q;

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

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *skk0_504 = buffer.data(skk0 + 504);
    const auto *skk0_507 = buffer.data(skk0 + 507);
    const auto *skk0_509 = buffer.data(skk0 + 509);
    const auto *skk0_510 = buffer.data(skk0 + 510);
    const auto *skk0_513 = buffer.data(skk0 + 513);
    const auto *skk0_514 = buffer.data(skk0 + 514);
    const auto *skk0_516 = buffer.data(skk0 + 516);
    const auto *skk0_518 = buffer.data(skk0 + 518);
    const auto *skk0_519 = buffer.data(skk0 + 519);
    const auto *skk0_521 = buffer.data(skk0 + 521);
    const auto *skk0_522 = buffer.data(skk0 + 522);
    const auto *skk0_524 = buffer.data(skk0 + 524);
    const auto *skk0_539 = buffer.data(skk0 + 539);

    const auto *ski_363 = buffer.data(ski + 363);
    const auto *ski_364 = buffer.data(ski + 364);
    const auto *ski_367 = buffer.data(ski + 367);
    const auto *ski_370 = buffer.data(ski + 370);
    const auto *ski_374 = buffer.data(ski + 374);
    const auto *ski_385 = buffer.data(ski + 385);
    const auto *ski_388 = buffer.data(ski + 388);
    const auto *ski_389 = buffer.data(ski + 389);
    const auto *ski_390 = buffer.data(ski + 390);
    const auto *ski_391 = buffer.data(ski + 391);
    const auto *ski_392 = buffer.data(ski + 392);
    const auto *ski_393 = buffer.data(ski + 393);
    const auto *ski_394 = buffer.data(ski + 394);
    const auto *ski_395 = buffer.data(ski + 395);
    const auto *ski_397 = buffer.data(ski + 397);
    const auto *ski_398 = buffer.data(ski + 398);
    const auto *ski_400 = buffer.data(ski + 400);
    const auto *ski_401 = buffer.data(ski + 401);
    const auto *ski_402 = buffer.data(ski + 402);
    const auto *ski_404 = buffer.data(ski + 404);
    const auto *ski_405 = buffer.data(ski + 405);
    const auto *ski_406 = buffer.data(ski + 406);
    const auto *ski_413 = buffer.data(ski + 413);
    const auto *ski_415 = buffer.data(ski + 415);
    const auto *ski_416 = buffer.data(ski + 416);
    const auto *ski_417 = buffer.data(ski + 417);
    const auto *ski_418 = buffer.data(ski + 418);
    const auto *ski_419 = buffer.data(ski + 419);
    const auto *ski_420 = buffer.data(ski + 420);
    const auto *ski_422 = buffer.data(ski + 422);
    const auto *ski_425 = buffer.data(ski + 425);
    const auto *ski_429 = buffer.data(ski + 429);
    const auto *ski_434 = buffer.data(ski + 434);
    const auto *ski_441 = buffer.data(ski + 441);
    const auto *ski_443 = buffer.data(ski + 443);
    const auto *ski_444 = buffer.data(ski + 444);
    const auto *ski_445 = buffer.data(ski + 445);
    const auto *ski_446 = buffer.data(ski + 446);
    const auto *ski_553 = buffer.data(ski + 553);
    const auto *ski_554 = buffer.data(ski + 554);
    const auto *ski_555 = buffer.data(ski + 555);
    const auto *ski_556 = buffer.data(ski + 556);
    const auto *ski_557 = buffer.data(ski + 557);
    const auto *ski_558 = buffer.data(ski + 558);
    const auto *ski_559 = buffer.data(ski + 559);
    const auto *ski_560 = buffer.data(ski + 560);
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
    const auto *ski_591 = buffer.data(ski + 591);
    const auto *ski_593 = buffer.data(ski + 593);
    const auto *ski_594 = buffer.data(ski + 594);
    const auto *ski_597 = buffer.data(ski + 597);
    const auto *ski_598 = buffer.data(ski + 598);
    const auto *ski_600 = buffer.data(ski + 600);
    const auto *ski_602 = buffer.data(ski + 602);
    const auto *ski_603 = buffer.data(ski + 603);
    const auto *ski_605 = buffer.data(ski + 605);
    const auto *ski_606 = buffer.data(ski + 606);
    const auto *ski_608 = buffer.data(ski + 608);
    const auto *ski_609 = buffer.data(ski + 609);
    const auto *ski_610 = buffer.data(ski + 610);
    const auto *ski_611 = buffer.data(ski + 611);
    const auto *ski_612 = buffer.data(ski + 612);
    const auto *ski_613 = buffer.data(ski + 613);
    const auto *ski_614 = buffer.data(ski + 614);
    const auto *ski_615 = buffer.data(ski + 615);

    const auto *skk1_504 = buffer.data(skk1 + 504);
    const auto *skk1_507 = buffer.data(skk1 + 507);
    const auto *skk1_509 = buffer.data(skk1 + 509);
    const auto *skk1_510 = buffer.data(skk1 + 510);
    const auto *skk1_513 = buffer.data(skk1 + 513);
    const auto *skk1_514 = buffer.data(skk1 + 514);
    const auto *skk1_516 = buffer.data(skk1 + 516);
    const auto *skk1_518 = buffer.data(skk1 + 518);
    const auto *skk1_519 = buffer.data(skk1 + 519);
    const auto *skk1_521 = buffer.data(skk1 + 521);
    const auto *skk1_522 = buffer.data(skk1 + 522);
    const auto *skk1_524 = buffer.data(skk1 + 524);
    const auto *skk1_539 = buffer.data(skk1 + 539);

    const auto *slh0_396 = buffer.data(slh0 + 396);
    const auto *slh0_397 = buffer.data(slh0 + 397);
    const auto *slh0_398 = buffer.data(slh0 + 398);
    const auto *slh0_414 = buffer.data(slh0 + 414);
    const auto *slh0_416 = buffer.data(slh0 + 416);
    const auto *slh0_417 = buffer.data(slh0 + 417);
    const auto *slh0_418 = buffer.data(slh0 + 418);
    const auto *slh0_419 = buffer.data(slh0 + 419);
    const auto *slh0_420 = buffer.data(slh0 + 420);
    const auto *slh0_423 = buffer.data(slh0 + 423);
    const auto *slh0_425 = buffer.data(slh0 + 425);
    const auto *slh0_426 = buffer.data(slh0 + 426);
    const auto *slh0_429 = buffer.data(slh0 + 429);
    const auto *slh0_430 = buffer.data(slh0 + 430);
    const auto *slh0_432 = buffer.data(slh0 + 432);
    const auto *slh0_434 = buffer.data(slh0 + 434);
    const auto *slh0_435 = buffer.data(slh0 + 435);
    const auto *slh0_437 = buffer.data(slh0 + 437);
    const auto *slh0_438 = buffer.data(slh0 + 438);
    const auto *slh0_439 = buffer.data(slh0 + 439);
    const auto *slh0_440 = buffer.data(slh0 + 440);
    const auto *slh0_441 = buffer.data(slh0 + 441);
    const auto *slh0_444 = buffer.data(slh0 + 444);
    const auto *slh0_446 = buffer.data(slh0 + 446);
    const auto *slh0_447 = buffer.data(slh0 + 447);
    const auto *slh0_450 = buffer.data(slh0 + 450);
    const auto *slh0_451 = buffer.data(slh0 + 451);
    const auto *slh0_453 = buffer.data(slh0 + 453);
    const auto *slh0_455 = buffer.data(slh0 + 455);
    const auto *slh0_456 = buffer.data(slh0 + 456);
    const auto *slh0_458 = buffer.data(slh0 + 458);
    const auto *slh0_459 = buffer.data(slh0 + 459);
    const auto *slh0_460 = buffer.data(slh0 + 460);
    const auto *slh0_461 = buffer.data(slh0 + 461);

    const auto *slh1_396 = buffer.data(slh1 + 396);
    const auto *slh1_397 = buffer.data(slh1 + 397);
    const auto *slh1_398 = buffer.data(slh1 + 398);
    const auto *slh1_414 = buffer.data(slh1 + 414);
    const auto *slh1_416 = buffer.data(slh1 + 416);
    const auto *slh1_417 = buffer.data(slh1 + 417);
    const auto *slh1_418 = buffer.data(slh1 + 418);
    const auto *slh1_419 = buffer.data(slh1 + 419);
    const auto *slh1_420 = buffer.data(slh1 + 420);
    const auto *slh1_423 = buffer.data(slh1 + 423);
    const auto *slh1_425 = buffer.data(slh1 + 425);
    const auto *slh1_426 = buffer.data(slh1 + 426);
    const auto *slh1_429 = buffer.data(slh1 + 429);
    const auto *slh1_430 = buffer.data(slh1 + 430);
    const auto *slh1_432 = buffer.data(slh1 + 432);
    const auto *slh1_434 = buffer.data(slh1 + 434);
    const auto *slh1_435 = buffer.data(slh1 + 435);
    const auto *slh1_437 = buffer.data(slh1 + 437);
    const auto *slh1_438 = buffer.data(slh1 + 438);
    const auto *slh1_439 = buffer.data(slh1 + 439);
    const auto *slh1_440 = buffer.data(slh1 + 440);
    const auto *slh1_441 = buffer.data(slh1 + 441);
    const auto *slh1_444 = buffer.data(slh1 + 444);
    const auto *slh1_446 = buffer.data(slh1 + 446);
    const auto *slh1_447 = buffer.data(slh1 + 447);
    const auto *slh1_450 = buffer.data(slh1 + 450);
    const auto *slh1_451 = buffer.data(slh1 + 451);
    const auto *slh1_453 = buffer.data(slh1 + 453);
    const auto *slh1_455 = buffer.data(slh1 + 455);
    const auto *slh1_456 = buffer.data(slh1 + 456);
    const auto *slh1_458 = buffer.data(slh1 + 458);
    const auto *slh1_459 = buffer.data(slh1 + 459);
    const auto *slh1_460 = buffer.data(slh1 + 460);
    const auto *slh1_461 = buffer.data(slh1 + 461);

    const auto *sli_528 = buffer.data(sli + 528);
    const auto *sli_529 = buffer.data(sli + 529);
    const auto *sli_530 = buffer.data(sli + 530);
    const auto *sli_531 = buffer.data(sli + 531);
    const auto *sli_532 = buffer.data(sli + 532);
    const auto *sli_534 = buffer.data(sli + 534);
    const auto *sli_535 = buffer.data(sli + 535);
    const auto *sli_537 = buffer.data(sli + 537);
    const auto *sli_538 = buffer.data(sli + 538);
    const auto *sli_541 = buffer.data(sli + 541);
    const auto *sli_542 = buffer.data(sli + 542);
    const auto *sli_546 = buffer.data(sli + 546);
    const auto *sli_553 = buffer.data(sli + 553);
    const auto *sli_554 = buffer.data(sli + 554);
    const auto *sli_555 = buffer.data(sli + 555);
    const auto *sli_556 = buffer.data(sli + 556);
    const auto *sli_557 = buffer.data(sli + 557);
    const auto *sli_558 = buffer.data(sli + 558);
    const auto *sli_559 = buffer.data(sli + 559);
    const auto *sli_560 = buffer.data(sli + 560);
    const auto *sli_562 = buffer.data(sli + 562);
    const auto *sli_563 = buffer.data(sli + 563);
    const auto *sli_565 = buffer.data(sli + 565);
    const auto *sli_566 = buffer.data(sli + 566);
    const auto *sli_569 = buffer.data(sli + 569);
    const auto *sli_570 = buffer.data(sli + 570);
    const auto *sli_572 = buffer.data(sli + 572);
    const auto *sli_574 = buffer.data(sli + 574);
    const auto *sli_575 = buffer.data(sli + 575);
    const auto *sli_577 = buffer.data(sli + 577);
    const auto *sli_578 = buffer.data(sli + 578);
    const auto *sli_580 = buffer.data(sli + 580);
    const auto *sli_581 = buffer.data(sli + 581);
    const auto *sli_582 = buffer.data(sli + 582);
    const auto *sli_583 = buffer.data(sli + 583);
    const auto *sli_584 = buffer.data(sli + 584);
    const auto *sli_585 = buffer.data(sli + 585);
    const auto *sli_586 = buffer.data(sli + 586);
    const auto *sli_587 = buffer.data(sli + 587);
    const auto *sli_588 = buffer.data(sli + 588);
    const auto *sli_590 = buffer.data(sli + 590);
    const auto *sli_591 = buffer.data(sli + 591);
    const auto *sli_593 = buffer.data(sli + 593);
    const auto *sli_594 = buffer.data(sli + 594);
    const auto *sli_597 = buffer.data(sli + 597);
    const auto *sli_598 = buffer.data(sli + 598);
    const auto *sli_600 = buffer.data(sli + 600);
    const auto *sli_602 = buffer.data(sli + 602);
    const auto *sli_603 = buffer.data(sli + 603);
    const auto *sli_605 = buffer.data(sli + 605);
    const auto *sli_606 = buffer.data(sli + 606);
    const auto *sli_608 = buffer.data(sli + 608);
    const auto *sli_609 = buffer.data(sli + 609);
    const auto *sli_610 = buffer.data(sli + 610);
    const auto *sli_611 = buffer.data(sli + 611);
    const auto *sli_612 = buffer.data(sli + 612);
    const auto *sli_613 = buffer.data(sli + 613);
    const auto *sli_614 = buffer.data(sli + 614);
    const auto *sli_615 = buffer.data(sli + 615);

#pragma omp simd aligned(t_679, t_680, t_681, pc_y, ski_388, ski_389, ski_390, slh0_396, \
                         slh0_397, slh0_398, slh1_396, slh1_397, slh1_398, sli_528, sli_529, \
                         sli_530 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_679[k] = f_14 * ski_388[k]
                   + f_6 * slh0_396[k]
                   - f_7 * slh1_396[k]
                   + f_3 * pc_y[k] * sli_528[k];

        t_680[k] = f_14 * ski_389[k]
                   + f_8 * slh0_397[k]
                   - f_9 * slh1_397[k]
                   + f_3 * pc_y[k] * sli_529[k];

        t_681[k] = f_14 * ski_390[k]
                   + f_10 * slh0_398[k]
                   - f_11 * slh1_398[k]
                   + f_3 * pc_y[k] * sli_530[k];
    }

#pragma omp simd aligned(t_682, t_683, t_684, t_685, pb_y, pc_y, pc_z, skk0_504, ski_363, \
                         ski_391, ski_392, skk1_504, slh0_398, slh1_398, sli_531, \
                         sli_532 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_682[k] = f_14 * ski_391[k]
                   + f_3 * pc_y[k] * sli_531[k];

        t_683[k] = f_15 * ski_363[k]
                   + f_1 * slh0_398[k]
                   - f_2 * slh1_398[k]
                   + f_3 * pc_z[k] * sli_531[k];

        t_684[k] = pb_y[k] * skk0_504[k]
                   - f_12 * pc_y[k] * skk1_504[k];

        t_685[k] = f_13 * ski_392[k]
                   + f_3 * pc_y[k] * sli_532[k];
    }

#pragma omp simd aligned(t_686, t_687, t_688, t_689, pb_y, pc_y, pc_z, skk0_507, skk0_509, \
                         ski_364, ski_393, ski_394, skk1_507, skk1_509, sli_532, \
                         sli_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_686[k] = f_16 * ski_364[k]
                   + f_3 * pc_z[k] * sli_532[k];

        t_687[k] = pb_y[k] * skk0_507[k]
                   + f_14 * ski_393[k]
                   - f_12 * pc_y[k] * skk1_507[k];

        t_688[k] = f_13 * ski_394[k]
                   + f_3 * pc_y[k] * sli_534[k];

        t_689[k] = pb_y[k] * skk0_509[k]
                   - f_12 * pc_y[k] * skk1_509[k];
    }

#pragma omp simd aligned(t_690, t_691, t_692, t_693, pb_y, pc_y, pc_z, skk0_510, skk0_513, \
                         ski_367, ski_395, ski_397, skk1_510, skk1_513, sli_535, \
                         sli_537 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_690[k] = pb_y[k] * skk0_510[k]
                   + f_15 * ski_395[k]
                   - f_12 * pc_y[k] * skk1_510[k];

        t_691[k] = f_16 * ski_367[k]
                   + f_3 * pc_z[k] * sli_535[k];

        t_692[k] = f_13 * ski_397[k]
                   + f_3 * pc_y[k] * sli_537[k];

        t_693[k] = pb_y[k] * skk0_513[k]
                   - f_12 * pc_y[k] * skk1_513[k];
    }

#pragma omp simd aligned(t_694, t_695, t_696, pb_y, pc_y, pc_z, skk0_514, skk0_516, ski_370, \
                         ski_398, ski_400, skk1_514, skk1_516, \
                         sli_538 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_694[k] = pb_y[k] * skk0_514[k]
                   + f_16 * ski_398[k]
                   - f_12 * pc_y[k] * skk1_514[k];

        t_695[k] = f_16 * ski_370[k]
                   + f_3 * pc_z[k] * sli_538[k];

        t_696[k] = pb_y[k] * skk0_516[k]
                   + f_14 * ski_400[k]
                   - f_12 * pc_y[k] * skk1_516[k];
    }

#pragma omp simd aligned(t_697, t_698, t_699, t_700, pb_y, pc_y, pc_z, skk0_518, skk0_519, \
                         ski_374, ski_401, ski_402, skk1_518, skk1_519, sli_541, \
                         sli_542 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_697[k] = f_13 * ski_401[k]
                   + f_3 * pc_y[k] * sli_541[k];

        t_698[k] = pb_y[k] * skk0_518[k]
                   - f_12 * pc_y[k] * skk1_518[k];

        t_699[k] = pb_y[k] * skk0_519[k]
                   + f_17 * ski_402[k]
                   - f_12 * pc_y[k] * skk1_519[k];

        t_700[k] = f_16 * ski_374[k]
                   + f_3 * pc_z[k] * sli_542[k];
    }

#pragma omp simd aligned(t_701, t_702, t_703, t_704, pb_y, pc_y, skk0_521, skk0_522, skk0_524, \
                         ski_404, ski_405, ski_406, skk1_521, skk1_522, skk1_524, \
                         sli_546 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_701[k] = pb_y[k] * skk0_521[k]
                   + f_15 * ski_404[k]
                   - f_12 * pc_y[k] * skk1_521[k];

        t_702[k] = pb_y[k] * skk0_522[k]
                   + f_14 * ski_405[k]
                   - f_12 * pc_y[k] * skk1_522[k];

        t_703[k] = f_13 * ski_406[k]
                   + f_3 * pc_y[k] * sli_546[k];

        t_704[k] = pb_y[k] * skk0_524[k]
                   - f_12 * pc_y[k] * skk1_524[k];
    }

#pragma omp simd aligned(t_705, t_706, t_707, t_708, t_709, pc_x, ski_553, ski_554, ski_555, \
                         ski_556, ski_557, sli_553, sli_554, sli_555, sli_556, \
                         sli_557 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_705[k] = f_15 * ski_553[k]
                   + f_3 * pc_x[k] * sli_553[k];

        t_706[k] = f_15 * ski_554[k]
                   + f_3 * pc_x[k] * sli_554[k];

        t_707[k] = f_15 * ski_555[k]
                   + f_3 * pc_x[k] * sli_555[k];

        t_708[k] = f_15 * ski_556[k]
                   + f_3 * pc_x[k] * sli_556[k];

        t_709[k] = f_15 * ski_557[k]
                   + f_3 * pc_x[k] * sli_557[k];
    }

#pragma omp simd aligned(t_710, t_711, t_712, t_713, pc_x, pc_y, pc_z, ski_385, ski_413, \
                         ski_558, ski_559, slh0_414, slh1_414, sli_553, sli_558, \
                         sli_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_710[k] = f_15 * ski_558[k]
                   + f_3 * pc_x[k] * sli_558[k];

        t_711[k] = f_15 * ski_559[k]
                   + f_3 * pc_x[k] * sli_559[k];

        t_712[k] = f_13 * ski_413[k]
                   + f_1 * slh0_414[k]
                   - f_2 * slh1_414[k]
                   + f_3 * pc_y[k] * sli_553[k];

        t_713[k] = f_16 * ski_385[k]
                   + f_3 * pc_z[k] * sli_553[k];
    }

#pragma omp simd aligned(t_714, t_715, t_716, pc_y, ski_415, ski_416, ski_417, slh0_416, \
                         slh0_417, slh0_418, slh1_416, slh1_417, slh1_418, sli_555, sli_556, \
                         sli_557 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_714[k] = f_13 * ski_415[k]
                   + f_4 * slh0_416[k]
                   - f_5 * slh1_416[k]
                   + f_3 * pc_y[k] * sli_555[k];

        t_715[k] = f_13 * ski_416[k]
                   + f_6 * slh0_417[k]
                   - f_7 * slh1_417[k]
                   + f_3 * pc_y[k] * sli_556[k];

        t_716[k] = f_13 * ski_417[k]
                   + f_8 * slh0_418[k]
                   - f_9 * slh1_418[k]
                   + f_3 * pc_y[k] * sli_557[k];
    }

#pragma omp simd aligned(t_717, t_718, t_719, pb_y, pc_y, skk0_539, ski_418, ski_419, \
                         skk1_539, slh0_419, slh1_419, sli_558, \
                         sli_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_717[k] = f_13 * ski_418[k]
                   + f_10 * slh0_419[k]
                   - f_11 * slh1_419[k]
                   + f_3 * pc_y[k] * sli_558[k];

        t_718[k] = f_13 * ski_419[k]
                   + f_3 * pc_y[k] * sli_559[k];

        t_719[k] = pb_y[k] * skk0_539[k]
                   - f_12 * pc_y[k] * skk1_539[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, t_723, pc_x, pc_y, pc_z, ski_392, ski_560, \
                         ski_563, slh0_420, slh0_423, slh1_420, slh1_423, sli_560, \
                         sli_563 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = f_15 * ski_560[k]
                   + f_1 * slh0_420[k]
                   - f_2 * slh1_420[k]
                   + f_3 * pc_x[k] * sli_560[k];

        t_721[k] = f_3 * pc_y[k] * sli_560[k];

        t_722[k] = f_17 * ski_392[k]
                   + f_3 * pc_z[k] * sli_560[k];

        t_723[k] = f_15 * ski_563[k]
                   + f_4 * slh0_423[k]
                   - f_5 * slh1_423[k]
                   + f_3 * pc_x[k] * sli_563[k];
    }

#pragma omp simd aligned(t_724, t_725, t_726, pc_x, pc_y, ski_565, ski_566, slh0_425, \
                         slh0_426, slh1_425, slh1_426, sli_562, sli_565, \
                         sli_566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_724[k] = f_3 * pc_y[k] * sli_562[k];

        t_725[k] = f_15 * ski_565[k]
                   + f_4 * slh0_425[k]
                   - f_5 * slh1_425[k]
                   + f_3 * pc_x[k] * sli_565[k];

        t_726[k] = f_15 * ski_566[k]
                   + f_6 * slh0_426[k]
                   - f_7 * slh1_426[k]
                   + f_3 * pc_x[k] * sli_566[k];
    }

#pragma omp simd aligned(t_727, t_728, t_729, pc_x, pc_y, pc_z, ski_395, ski_569, slh0_429, \
                         slh1_429, sli_563, sli_565, sli_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_727[k] = f_17 * ski_395[k]
                   + f_3 * pc_z[k] * sli_563[k];

        t_728[k] = f_3 * pc_y[k] * sli_565[k];

        t_729[k] = f_15 * ski_569[k]
                   + f_6 * slh0_429[k]
                   - f_7 * slh1_429[k]
                   + f_3 * pc_x[k] * sli_569[k];
    }

#pragma omp simd aligned(t_730, t_731, t_732, pc_x, pc_z, ski_398, ski_570, ski_572, slh0_430, \
                         slh0_432, slh1_430, slh1_432, sli_566, sli_570, \
                         sli_572 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_730[k] = f_15 * ski_570[k]
                   + f_8 * slh0_430[k]
                   - f_9 * slh1_430[k]
                   + f_3 * pc_x[k] * sli_570[k];

        t_731[k] = f_17 * ski_398[k]
                   + f_3 * pc_z[k] * sli_566[k];

        t_732[k] = f_15 * ski_572[k]
                   + f_8 * slh0_432[k]
                   - f_9 * slh1_432[k]
                   + f_3 * pc_x[k] * sli_572[k];
    }

#pragma omp simd aligned(t_733, t_734, t_735, pc_x, pc_y, ski_574, ski_575, slh0_434, \
                         slh0_435, slh1_434, slh1_435, sli_569, sli_574, \
                         sli_575 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_733[k] = f_3 * pc_y[k] * sli_569[k];

        t_734[k] = f_15 * ski_574[k]
                   + f_8 * slh0_434[k]
                   - f_9 * slh1_434[k]
                   + f_3 * pc_x[k] * sli_574[k];

        t_735[k] = f_15 * ski_575[k]
                   + f_10 * slh0_435[k]
                   - f_11 * slh1_435[k]
                   + f_3 * pc_x[k] * sli_575[k];
    }

#pragma omp simd aligned(t_736, t_737, t_738, pc_x, pc_z, ski_402, ski_577, ski_578, slh0_437, \
                         slh0_438, slh1_437, slh1_438, sli_570, sli_577, \
                         sli_578 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_736[k] = f_17 * ski_402[k]
                   + f_3 * pc_z[k] * sli_570[k];

        t_737[k] = f_15 * ski_577[k]
                   + f_10 * slh0_437[k]
                   - f_11 * slh1_437[k]
                   + f_3 * pc_x[k] * sli_577[k];

        t_738[k] = f_15 * ski_578[k]
                   + f_10 * slh0_438[k]
                   - f_11 * slh1_438[k]
                   + f_3 * pc_x[k] * sli_578[k];
    }

#pragma omp simd aligned(t_739, t_740, t_741, t_742, pc_x, pc_y, ski_580, ski_581, ski_582, \
                         slh0_440, slh1_440, sli_574, sli_580, sli_581, \
                         sli_582 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_739[k] = f_3 * pc_y[k] * sli_574[k];

        t_740[k] = f_15 * ski_580[k]
                   + f_10 * slh0_440[k]
                   - f_11 * slh1_440[k]
                   + f_3 * pc_x[k] * sli_580[k];

        t_741[k] = f_15 * ski_581[k]
                   + f_3 * pc_x[k] * sli_581[k];

        t_742[k] = f_15 * ski_582[k]
                   + f_3 * pc_x[k] * sli_582[k];
    }

#pragma omp simd aligned(t_743, t_744, t_745, t_746, t_747, pc_x, ski_583, ski_584, ski_585, \
                         ski_586, ski_587, sli_583, sli_584, sli_585, sli_586, \
                         sli_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_743[k] = f_15 * ski_583[k]
                   + f_3 * pc_x[k] * sli_583[k];

        t_744[k] = f_15 * ski_584[k]
                   + f_3 * pc_x[k] * sli_584[k];

        t_745[k] = f_15 * ski_585[k]
                   + f_3 * pc_x[k] * sli_585[k];

        t_746[k] = f_15 * ski_586[k]
                   + f_3 * pc_x[k] * sli_586[k];

        t_747[k] = f_15 * ski_587[k]
                   + f_3 * pc_x[k] * sli_587[k];
    }

#pragma omp simd aligned(t_748, t_749, t_750, t_751, pc_y, pc_z, ski_413, slh0_435, slh0_437, \
                         slh0_438, slh1_435, slh1_437, slh1_438, sli_581, sli_583, \
                         sli_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_748[k] = f_1 * slh0_435[k]
                   - f_2 * slh1_435[k]
                   + f_3 * pc_y[k] * sli_581[k];

        t_749[k] = f_17 * ski_413[k]
                   + f_3 * pc_z[k] * sli_581[k];

        t_750[k] = f_4 * slh0_437[k]
                   - f_5 * slh1_437[k]
                   + f_3 * pc_y[k] * sli_583[k];

        t_751[k] = f_6 * slh0_438[k]
                   - f_7 * slh1_438[k]
                   + f_3 * pc_y[k] * sli_584[k];
    }

#pragma omp simd aligned(t_752, t_753, t_754, t_755, pc_y, pc_z, ski_419, slh0_439, slh0_440, \
                         slh1_439, slh1_440, sli_585, sli_586, \
                         sli_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_752[k] = f_8 * slh0_439[k]
                   - f_9 * slh1_439[k]
                   + f_3 * pc_y[k] * sli_585[k];

        t_753[k] = f_10 * slh0_440[k]
                   - f_11 * slh1_440[k]
                   + f_3 * pc_y[k] * sli_586[k];

        t_754[k] = f_3 * pc_y[k] * sli_587[k];

        t_755[k] = f_17 * ski_419[k]
                   + f_1 * slh0_440[k]
                   - f_2 * slh1_440[k]
                   + f_3 * pc_z[k] * sli_587[k];
    }

#pragma omp simd aligned(t_756, t_757, t_758, t_759, pc_x, pc_y, pc_z, ski_420, ski_588, \
                         ski_591, slh0_441, slh0_444, slh1_441, slh1_444, sli_588, \
                         sli_591 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_756[k] = f_14 * ski_588[k]
                   + f_1 * slh0_441[k]
                   - f_2 * slh1_441[k]
                   + f_3 * pc_x[k] * sli_588[k];

        t_757[k] = f_19 * ski_420[k]
                   + f_3 * pc_y[k] * sli_588[k];

        t_758[k] = f_3 * pc_z[k] * sli_588[k];

        t_759[k] = f_14 * ski_591[k]
                   + f_4 * slh0_444[k]
                   - f_5 * slh1_444[k]
                   + f_3 * pc_x[k] * sli_591[k];
    }

#pragma omp simd aligned(t_760, t_761, t_762, pc_x, pc_y, ski_422, ski_593, ski_594, slh0_446, \
                         slh0_447, slh1_446, slh1_447, sli_590, sli_593, \
                         sli_594 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_760[k] = f_19 * ski_422[k]
                   + f_3 * pc_y[k] * sli_590[k];

        t_761[k] = f_14 * ski_593[k]
                   + f_4 * slh0_446[k]
                   - f_5 * slh1_446[k]
                   + f_3 * pc_x[k] * sli_593[k];

        t_762[k] = f_14 * ski_594[k]
                   + f_6 * slh0_447[k]
                   - f_7 * slh1_447[k]
                   + f_3 * pc_x[k] * sli_594[k];
    }

#pragma omp simd aligned(t_763, t_764, t_765, pc_x, pc_y, pc_z, ski_425, ski_597, slh0_450, \
                         slh1_450, sli_591, sli_593, sli_597 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_763[k] = f_3 * pc_z[k] * sli_591[k];

        t_764[k] = f_19 * ski_425[k]
                   + f_3 * pc_y[k] * sli_593[k];

        t_765[k] = f_14 * ski_597[k]
                   + f_6 * slh0_450[k]
                   - f_7 * slh1_450[k]
                   + f_3 * pc_x[k] * sli_597[k];
    }

#pragma omp simd aligned(t_766, t_767, t_768, pc_x, pc_z, ski_598, ski_600, slh0_451, \
                         slh0_453, slh1_451, slh1_453, sli_594, sli_598, \
                         sli_600 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_766[k] = f_14 * ski_598[k]
                   + f_8 * slh0_451[k]
                   - f_9 * slh1_451[k]
                   + f_3 * pc_x[k] * sli_598[k];

        t_767[k] = f_3 * pc_z[k] * sli_594[k];

        t_768[k] = f_14 * ski_600[k]
                   + f_8 * slh0_453[k]
                   - f_9 * slh1_453[k]
                   + f_3 * pc_x[k] * sli_600[k];
    }

#pragma omp simd aligned(t_769, t_770, t_771, pc_x, pc_y, ski_429, ski_602, ski_603, slh0_455, \
                         slh0_456, slh1_455, slh1_456, sli_597, sli_602, \
                         sli_603 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_769[k] = f_19 * ski_429[k]
                   + f_3 * pc_y[k] * sli_597[k];

        t_770[k] = f_14 * ski_602[k]
                   + f_8 * slh0_455[k]
                   - f_9 * slh1_455[k]
                   + f_3 * pc_x[k] * sli_602[k];

        t_771[k] = f_14 * ski_603[k]
                   + f_10 * slh0_456[k]
                   - f_11 * slh1_456[k]
                   + f_3 * pc_x[k] * sli_603[k];
    }

#pragma omp simd aligned(t_772, t_773, t_774, pc_x, pc_z, ski_605, ski_606, slh0_458, \
                         slh0_459, slh1_458, slh1_459, sli_598, sli_605, \
                         sli_606 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_772[k] = f_3 * pc_z[k] * sli_598[k];

        t_773[k] = f_14 * ski_605[k]
                   + f_10 * slh0_458[k]
                   - f_11 * slh1_458[k]
                   + f_3 * pc_x[k] * sli_605[k];

        t_774[k] = f_14 * ski_606[k]
                   + f_10 * slh0_459[k]
                   - f_11 * slh1_459[k]
                   + f_3 * pc_x[k] * sli_606[k];
    }

#pragma omp simd aligned(t_775, t_776, t_777, t_778, pc_x, pc_y, ski_434, ski_608, ski_609, \
                         ski_610, slh0_461, slh1_461, sli_602, sli_608, sli_609, \
                         sli_610 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_775[k] = f_19 * ski_434[k]
                   + f_3 * pc_y[k] * sli_602[k];

        t_776[k] = f_14 * ski_608[k]
                   + f_10 * slh0_461[k]
                   - f_11 * slh1_461[k]
                   + f_3 * pc_x[k] * sli_608[k];

        t_777[k] = f_14 * ski_609[k]
                   + f_3 * pc_x[k] * sli_609[k];

        t_778[k] = f_14 * ski_610[k]
                   + f_3 * pc_x[k] * sli_610[k];
    }

#pragma omp simd aligned(t_779, t_780, t_781, t_782, t_783, pc_x, ski_611, ski_612, ski_613, \
                         ski_614, ski_615, sli_611, sli_612, sli_613, sli_614, \
                         sli_615 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_779[k] = f_14 * ski_611[k]
                   + f_3 * pc_x[k] * sli_611[k];

        t_780[k] = f_14 * ski_612[k]
                   + f_3 * pc_x[k] * sli_612[k];

        t_781[k] = f_14 * ski_613[k]
                   + f_3 * pc_x[k] * sli_613[k];

        t_782[k] = f_14 * ski_614[k]
                   + f_3 * pc_x[k] * sli_614[k];

        t_783[k] = f_14 * ski_615[k]
                   + f_3 * pc_x[k] * sli_615[k];
    }

#pragma omp simd aligned(t_784, t_785, t_786, pc_y, pc_z, ski_441, ski_443, slh0_456, \
                         slh0_458, slh1_456, slh1_458, sli_609, \
                         sli_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_784[k] = f_19 * ski_441[k]
                   + f_1 * slh0_456[k]
                   - f_2 * slh1_456[k]
                   + f_3 * pc_y[k] * sli_609[k];

        t_785[k] = f_3 * pc_z[k] * sli_609[k];

        t_786[k] = f_19 * ski_443[k]
                   + f_4 * slh0_458[k]
                   - f_5 * slh1_458[k]
                   + f_3 * pc_y[k] * sli_611[k];
    }

#pragma omp simd aligned(t_787, t_788, t_789, pc_y, ski_444, ski_445, ski_446, slh0_459, \
                         slh0_460, slh0_461, slh1_459, slh1_460, slh1_461, sli_612, sli_613, \
                         sli_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_787[k] = f_19 * ski_444[k]
                   + f_6 * slh0_459[k]
                   - f_7 * slh1_459[k]
                   + f_3 * pc_y[k] * sli_612[k];

        t_788[k] = f_19 * ski_445[k]
                   + f_8 * slh0_460[k]
                   - f_9 * slh1_460[k]
                   + f_3 * pc_y[k] * sli_613[k];

        t_789[k] = f_19 * ski_446[k]
                   + f_10 * slh0_461[k]
                   - f_11 * slh1_461[k]
                   + f_3 * pc_y[k] * sli_614[k];
    }
}

static auto
compute_prim_slk_three_center_electron_repulsion_0_piece7(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t skk0,
                                                          const size_t ski, const size_t skk1,
                                                          const size_t slh0, const size_t slh1,
                                                          const size_t sli, const size_t ncols,
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
    const auto f_19 = 3.0 / q;

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

    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *skk0_540 = buffer.data(skk0 + 540);
    const auto *skk0_543 = buffer.data(skk0 + 543);
    const auto *skk0_546 = buffer.data(skk0 + 546);
    const auto *skk0_550 = buffer.data(skk0 + 550);
    const auto *skk0_552 = buffer.data(skk0 + 552);
    const auto *skk0_555 = buffer.data(skk0 + 555);
    const auto *skk0_557 = buffer.data(skk0 + 557);
    const auto *skk0_558 = buffer.data(skk0 + 558);
    const auto *skk0_568 = buffer.data(skk0 + 568);

    const auto *ski_420 = buffer.data(ski + 420);
    const auto *ski_423 = buffer.data(ski + 423);
    const auto *ski_426 = buffer.data(ski + 426);
    const auto *ski_427 = buffer.data(ski + 427);
    const auto *ski_430 = buffer.data(ski + 430);
    const auto *ski_431 = buffer.data(ski + 431);
    const auto *ski_432 = buffer.data(ski + 432);
    const auto *ski_441 = buffer.data(ski + 441);
    const auto *ski_447 = buffer.data(ski + 447);
    const auto *ski_448 = buffer.data(ski + 448);
    const auto *ski_450 = buffer.data(ski + 450);
    const auto *ski_451 = buffer.data(ski + 451);
    const auto *ski_453 = buffer.data(ski + 453);
    const auto *ski_454 = buffer.data(ski + 454);
    const auto *ski_457 = buffer.data(ski + 457);
    const auto *ski_458 = buffer.data(ski + 458);
    const auto *ski_462 = buffer.data(ski + 462);
    const auto *ski_469 = buffer.data(ski + 469);
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
    const auto *ski_490 = buffer.data(ski + 490);
    const auto *ski_497 = buffer.data(ski + 497);
    const auto *ski_499 = buffer.data(ski + 499);
    const auto *ski_500 = buffer.data(ski + 500);
    const auto *ski_501 = buffer.data(ski + 501);
    const auto *ski_502 = buffer.data(ski + 502);
    const auto *ski_503 = buffer.data(ski + 503);
    const auto *ski_504 = buffer.data(ski + 504);
    const auto *ski_506 = buffer.data(ski + 506);
    const auto *ski_509 = buffer.data(ski + 509);
    const auto *ski_513 = buffer.data(ski + 513);
    const auto *ski_518 = buffer.data(ski + 518);
    const auto *ski_525 = buffer.data(ski + 525);
    const auto *ski_527 = buffer.data(ski + 527);
    const auto *ski_528 = buffer.data(ski + 528);
    const auto *ski_529 = buffer.data(ski + 529);
    const auto *ski_530 = buffer.data(ski + 530);
    const auto *ski_621 = buffer.data(ski + 621);
    const auto *ski_625 = buffer.data(ski + 625);
    const auto *ski_630 = buffer.data(ski + 630);
    const auto *ski_636 = buffer.data(ski + 636);
    const auto *ski_637 = buffer.data(ski + 637);
    const auto *ski_638 = buffer.data(ski + 638);
    const auto *ski_639 = buffer.data(ski + 639);
    const auto *ski_640 = buffer.data(ski + 640);
    const auto *ski_641 = buffer.data(ski + 641);
    const auto *ski_642 = buffer.data(ski + 642);
    const auto *ski_643 = buffer.data(ski + 643);
    const auto *ski_644 = buffer.data(ski + 644);
    const auto *ski_647 = buffer.data(ski + 647);
    const auto *ski_649 = buffer.data(ski + 649);
    const auto *ski_650 = buffer.data(ski + 650);
    const auto *ski_653 = buffer.data(ski + 653);
    const auto *ski_654 = buffer.data(ski + 654);
    const auto *ski_656 = buffer.data(ski + 656);
    const auto *ski_658 = buffer.data(ski + 658);
    const auto *ski_659 = buffer.data(ski + 659);
    const auto *ski_661 = buffer.data(ski + 661);
    const auto *ski_662 = buffer.data(ski + 662);
    const auto *ski_664 = buffer.data(ski + 664);
    const auto *ski_665 = buffer.data(ski + 665);
    const auto *ski_666 = buffer.data(ski + 666);
    const auto *ski_667 = buffer.data(ski + 667);
    const auto *ski_668 = buffer.data(ski + 668);
    const auto *ski_669 = buffer.data(ski + 669);
    const auto *ski_670 = buffer.data(ski + 670);
    const auto *ski_671 = buffer.data(ski + 671);
    const auto *ski_672 = buffer.data(ski + 672);
    const auto *ski_675 = buffer.data(ski + 675);
    const auto *ski_677 = buffer.data(ski + 677);
    const auto *ski_678 = buffer.data(ski + 678);
    const auto *ski_681 = buffer.data(ski + 681);
    const auto *ski_682 = buffer.data(ski + 682);
    const auto *ski_684 = buffer.data(ski + 684);
    const auto *ski_686 = buffer.data(ski + 686);
    const auto *ski_687 = buffer.data(ski + 687);
    const auto *ski_689 = buffer.data(ski + 689);
    const auto *ski_690 = buffer.data(ski + 690);
    const auto *ski_692 = buffer.data(ski + 692);
    const auto *ski_693 = buffer.data(ski + 693);
    const auto *ski_694 = buffer.data(ski + 694);
    const auto *ski_695 = buffer.data(ski + 695);
    const auto *ski_696 = buffer.data(ski + 696);
    const auto *ski_697 = buffer.data(ski + 697);
    const auto *ski_698 = buffer.data(ski + 698);
    const auto *ski_699 = buffer.data(ski + 699);

    const auto *skk1_540 = buffer.data(skk1 + 540);
    const auto *skk1_543 = buffer.data(skk1 + 543);
    const auto *skk1_546 = buffer.data(skk1 + 546);
    const auto *skk1_550 = buffer.data(skk1 + 550);
    const auto *skk1_552 = buffer.data(skk1 + 552);
    const auto *skk1_555 = buffer.data(skk1 + 555);
    const auto *skk1_557 = buffer.data(skk1 + 557);
    const auto *skk1_558 = buffer.data(skk1 + 558);
    const auto *skk1_568 = buffer.data(skk1 + 568);

    const auto *slh0_461 = buffer.data(slh0 + 461);
    const auto *slh0_467 = buffer.data(slh0 + 467);
    const auto *slh0_471 = buffer.data(slh0 + 471);
    const auto *slh0_476 = buffer.data(slh0 + 476);
    const auto *slh0_479 = buffer.data(slh0 + 479);
    const auto *slh0_480 = buffer.data(slh0 + 480);
    const auto *slh0_481 = buffer.data(slh0 + 481);
    const auto *slh0_482 = buffer.data(slh0 + 482);
    const auto *slh0_483 = buffer.data(slh0 + 483);
    const auto *slh0_486 = buffer.data(slh0 + 486);
    const auto *slh0_488 = buffer.data(slh0 + 488);
    const auto *slh0_489 = buffer.data(slh0 + 489);
    const auto *slh0_492 = buffer.data(slh0 + 492);
    const auto *slh0_493 = buffer.data(slh0 + 493);
    const auto *slh0_495 = buffer.data(slh0 + 495);
    const auto *slh0_497 = buffer.data(slh0 + 497);
    const auto *slh0_498 = buffer.data(slh0 + 498);
    const auto *slh0_500 = buffer.data(slh0 + 500);
    const auto *slh0_501 = buffer.data(slh0 + 501);
    const auto *slh0_502 = buffer.data(slh0 + 502);
    const auto *slh0_503 = buffer.data(slh0 + 503);
    const auto *slh0_504 = buffer.data(slh0 + 504);
    const auto *slh0_507 = buffer.data(slh0 + 507);
    const auto *slh0_509 = buffer.data(slh0 + 509);
    const auto *slh0_510 = buffer.data(slh0 + 510);
    const auto *slh0_513 = buffer.data(slh0 + 513);
    const auto *slh0_514 = buffer.data(slh0 + 514);
    const auto *slh0_516 = buffer.data(slh0 + 516);
    const auto *slh0_518 = buffer.data(slh0 + 518);
    const auto *slh0_519 = buffer.data(slh0 + 519);
    const auto *slh0_521 = buffer.data(slh0 + 521);
    const auto *slh0_522 = buffer.data(slh0 + 522);
    const auto *slh0_523 = buffer.data(slh0 + 523);
    const auto *slh0_524 = buffer.data(slh0 + 524);

    const auto *slh1_461 = buffer.data(slh1 + 461);
    const auto *slh1_467 = buffer.data(slh1 + 467);
    const auto *slh1_471 = buffer.data(slh1 + 471);
    const auto *slh1_476 = buffer.data(slh1 + 476);
    const auto *slh1_479 = buffer.data(slh1 + 479);
    const auto *slh1_480 = buffer.data(slh1 + 480);
    const auto *slh1_481 = buffer.data(slh1 + 481);
    const auto *slh1_482 = buffer.data(slh1 + 482);
    const auto *slh1_483 = buffer.data(slh1 + 483);
    const auto *slh1_486 = buffer.data(slh1 + 486);
    const auto *slh1_488 = buffer.data(slh1 + 488);
    const auto *slh1_489 = buffer.data(slh1 + 489);
    const auto *slh1_492 = buffer.data(slh1 + 492);
    const auto *slh1_493 = buffer.data(slh1 + 493);
    const auto *slh1_495 = buffer.data(slh1 + 495);
    const auto *slh1_497 = buffer.data(slh1 + 497);
    const auto *slh1_498 = buffer.data(slh1 + 498);
    const auto *slh1_500 = buffer.data(slh1 + 500);
    const auto *slh1_501 = buffer.data(slh1 + 501);
    const auto *slh1_502 = buffer.data(slh1 + 502);
    const auto *slh1_503 = buffer.data(slh1 + 503);
    const auto *slh1_504 = buffer.data(slh1 + 504);
    const auto *slh1_507 = buffer.data(slh1 + 507);
    const auto *slh1_509 = buffer.data(slh1 + 509);
    const auto *slh1_510 = buffer.data(slh1 + 510);
    const auto *slh1_513 = buffer.data(slh1 + 513);
    const auto *slh1_514 = buffer.data(slh1 + 514);
    const auto *slh1_516 = buffer.data(slh1 + 516);
    const auto *slh1_518 = buffer.data(slh1 + 518);
    const auto *slh1_519 = buffer.data(slh1 + 519);
    const auto *slh1_521 = buffer.data(slh1 + 521);
    const auto *slh1_522 = buffer.data(slh1 + 522);
    const auto *slh1_523 = buffer.data(slh1 + 523);
    const auto *slh1_524 = buffer.data(slh1 + 524);

    const auto *sli_615 = buffer.data(sli + 615);
    const auto *sli_616 = buffer.data(sli + 616);
    const auto *sli_618 = buffer.data(sli + 618);
    const auto *sli_619 = buffer.data(sli + 619);
    const auto *sli_621 = buffer.data(sli + 621);
    const auto *sli_622 = buffer.data(sli + 622);
    const auto *sli_625 = buffer.data(sli + 625);
    const auto *sli_626 = buffer.data(sli + 626);
    const auto *sli_630 = buffer.data(sli + 630);
    const auto *sli_636 = buffer.data(sli + 636);
    const auto *sli_637 = buffer.data(sli + 637);
    const auto *sli_638 = buffer.data(sli + 638);
    const auto *sli_639 = buffer.data(sli + 639);
    const auto *sli_640 = buffer.data(sli + 640);
    const auto *sli_641 = buffer.data(sli + 641);
    const auto *sli_642 = buffer.data(sli + 642);
    const auto *sli_643 = buffer.data(sli + 643);
    const auto *sli_644 = buffer.data(sli + 644);
    const auto *sli_646 = buffer.data(sli + 646);
    const auto *sli_647 = buffer.data(sli + 647);
    const auto *sli_649 = buffer.data(sli + 649);
    const auto *sli_650 = buffer.data(sli + 650);
    const auto *sli_653 = buffer.data(sli + 653);
    const auto *sli_654 = buffer.data(sli + 654);
    const auto *sli_656 = buffer.data(sli + 656);
    const auto *sli_658 = buffer.data(sli + 658);
    const auto *sli_659 = buffer.data(sli + 659);
    const auto *sli_661 = buffer.data(sli + 661);
    const auto *sli_662 = buffer.data(sli + 662);
    const auto *sli_664 = buffer.data(sli + 664);
    const auto *sli_665 = buffer.data(sli + 665);
    const auto *sli_666 = buffer.data(sli + 666);
    const auto *sli_667 = buffer.data(sli + 667);
    const auto *sli_668 = buffer.data(sli + 668);
    const auto *sli_669 = buffer.data(sli + 669);
    const auto *sli_670 = buffer.data(sli + 670);
    const auto *sli_671 = buffer.data(sli + 671);
    const auto *sli_672 = buffer.data(sli + 672);
    const auto *sli_674 = buffer.data(sli + 674);
    const auto *sli_675 = buffer.data(sli + 675);
    const auto *sli_677 = buffer.data(sli + 677);
    const auto *sli_678 = buffer.data(sli + 678);
    const auto *sli_681 = buffer.data(sli + 681);
    const auto *sli_682 = buffer.data(sli + 682);
    const auto *sli_684 = buffer.data(sli + 684);
    const auto *sli_686 = buffer.data(sli + 686);
    const auto *sli_687 = buffer.data(sli + 687);
    const auto *sli_689 = buffer.data(sli + 689);
    const auto *sli_690 = buffer.data(sli + 690);
    const auto *sli_692 = buffer.data(sli + 692);
    const auto *sli_693 = buffer.data(sli + 693);
    const auto *sli_694 = buffer.data(sli + 694);
    const auto *sli_695 = buffer.data(sli + 695);
    const auto *sli_696 = buffer.data(sli + 696);
    const auto *sli_697 = buffer.data(sli + 697);
    const auto *sli_698 = buffer.data(sli + 698);
    const auto *sli_699 = buffer.data(sli + 699);

#pragma omp simd aligned(t_790, t_791, t_792, t_793, pb_z, pc_y, pc_z, skk0_540, ski_447, \
                         ski_448, skk1_540, slh0_461, slh1_461, sli_615, \
                         sli_616 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_790[k] = f_19 * ski_447[k]
                   + f_3 * pc_y[k] * sli_615[k];

        t_791[k] = f_1 * slh0_461[k]
                   - f_2 * slh1_461[k]
                   + f_3 * pc_z[k] * sli_615[k];

        t_792[k] = pb_z[k] * skk0_540[k]
                   - f_12 * pc_z[k] * skk1_540[k];

        t_793[k] = f_17 * ski_448[k]
                   + f_3 * pc_y[k] * sli_616[k];
    }

#pragma omp simd aligned(t_794, t_795, t_796, pb_z, pc_y, pc_z, skk0_543, ski_420, ski_450, \
                         skk1_543, sli_616, sli_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_794[k] = f_13 * ski_420[k]
                   + f_3 * pc_z[k] * sli_616[k];

        t_795[k] = pb_z[k] * skk0_543[k]
                   - f_12 * pc_z[k] * skk1_543[k];

        t_796[k] = f_17 * ski_450[k]
                   + f_3 * pc_y[k] * sli_618[k];
    }

#pragma omp simd aligned(t_797, t_798, t_799, pb_z, pc_x, pc_z, skk0_546, ski_423, ski_621, \
                         skk1_546, slh0_467, slh1_467, sli_619, \
                         sli_621 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_797[k] = f_14 * ski_621[k]
                   + f_4 * slh0_467[k]
                   - f_5 * slh1_467[k]
                   + f_3 * pc_x[k] * sli_621[k];

        t_798[k] = pb_z[k] * skk0_546[k]
                   - f_12 * pc_z[k] * skk1_546[k];

        t_799[k] = f_13 * ski_423[k]
                   + f_3 * pc_z[k] * sli_619[k];
    }

#pragma omp simd aligned(t_800, t_801, t_802, pb_z, pc_x, pc_y, pc_z, skk0_550, ski_453, \
                         ski_625, skk1_550, slh0_471, slh1_471, sli_621, \
                         sli_625 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_800[k] = f_17 * ski_453[k]
                   + f_3 * pc_y[k] * sli_621[k];

        t_801[k] = f_14 * ski_625[k]
                   + f_6 * slh0_471[k]
                   - f_7 * slh1_471[k]
                   + f_3 * pc_x[k] * sli_625[k];

        t_802[k] = pb_z[k] * skk0_550[k]
                   - f_12 * pc_z[k] * skk1_550[k];
    }

#pragma omp simd aligned(t_803, t_804, t_805, pb_z, pc_y, pc_z, skk0_552, ski_426, ski_427, \
                         ski_457, skk1_552, sli_622, sli_625 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_803[k] = f_13 * ski_426[k]
                   + f_3 * pc_z[k] * sli_622[k];

        t_804[k] = pb_z[k] * skk0_552[k]
                   + f_14 * ski_427[k]
                   - f_12 * pc_z[k] * skk1_552[k];

        t_805[k] = f_17 * ski_457[k]
                   + f_3 * pc_y[k] * sli_625[k];
    }

#pragma omp simd aligned(t_806, t_807, t_808, pb_z, pc_x, pc_z, skk0_555, ski_430, ski_630, \
                         skk1_555, slh0_476, slh1_476, sli_626, \
                         sli_630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_806[k] = f_14 * ski_630[k]
                   + f_8 * slh0_476[k]
                   - f_9 * slh1_476[k]
                   + f_3 * pc_x[k] * sli_630[k];

        t_807[k] = pb_z[k] * skk0_555[k]
                   - f_12 * pc_z[k] * skk1_555[k];

        t_808[k] = f_13 * ski_430[k]
                   + f_3 * pc_z[k] * sli_626[k];
    }

#pragma omp simd aligned(t_809, t_810, t_811, pb_z, pc_y, pc_z, skk0_557, skk0_558, ski_431, \
                         ski_432, ski_462, skk1_557, skk1_558, \
                         sli_630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_809[k] = pb_z[k] * skk0_557[k]
                   + f_14 * ski_431[k]
                   - f_12 * pc_z[k] * skk1_557[k];

        t_810[k] = pb_z[k] * skk0_558[k]
                   + f_15 * ski_432[k]
                   - f_12 * pc_z[k] * skk1_558[k];

        t_811[k] = f_17 * ski_462[k]
                   + f_3 * pc_y[k] * sli_630[k];
    }

#pragma omp simd aligned(t_812, t_813, t_814, t_815, pc_x, ski_636, ski_637, ski_638, ski_639, \
                         slh0_482, slh1_482, sli_636, sli_637, sli_638, \
                         sli_639 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_812[k] = f_14 * ski_636[k]
                   + f_10 * slh0_482[k]
                   - f_11 * slh1_482[k]
                   + f_3 * pc_x[k] * sli_636[k];

        t_813[k] = f_14 * ski_637[k]
                   + f_3 * pc_x[k] * sli_637[k];

        t_814[k] = f_14 * ski_638[k]
                   + f_3 * pc_x[k] * sli_638[k];

        t_815[k] = f_14 * ski_639[k]
                   + f_3 * pc_x[k] * sli_639[k];
    }

#pragma omp simd aligned(t_816, t_817, t_818, t_819, pc_x, ski_640, ski_641, ski_642, ski_643, \
                         sli_640, sli_641, sli_642, sli_643 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_816[k] = f_14 * ski_640[k]
                   + f_3 * pc_x[k] * sli_640[k];

        t_817[k] = f_14 * ski_641[k]
                   + f_3 * pc_x[k] * sli_641[k];

        t_818[k] = f_14 * ski_642[k]
                   + f_3 * pc_x[k] * sli_642[k];

        t_819[k] = f_14 * ski_643[k]
                   + f_3 * pc_x[k] * sli_643[k];
    }

#pragma omp simd aligned(t_820, t_821, t_822, pb_z, pc_y, pc_z, skk0_568, ski_441, ski_471, \
                         skk1_568, slh0_479, slh1_479, sli_637, \
                         sli_639 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_820[k] = pb_z[k] * skk0_568[k]
                   - f_12 * pc_z[k] * skk1_568[k];

        t_821[k] = f_13 * ski_441[k]
                   + f_3 * pc_z[k] * sli_637[k];

        t_822[k] = f_17 * ski_471[k]
                   + f_4 * slh0_479[k]
                   - f_5 * slh1_479[k]
                   + f_3 * pc_y[k] * sli_639[k];
    }

#pragma omp simd aligned(t_823, t_824, t_825, pc_y, ski_472, ski_473, ski_474, slh0_480, \
                         slh0_481, slh0_482, slh1_480, slh1_481, slh1_482, sli_640, sli_641, \
                         sli_642 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_823[k] = f_17 * ski_472[k]
                   + f_6 * slh0_480[k]
                   - f_7 * slh1_480[k]
                   + f_3 * pc_y[k] * sli_640[k];

        t_824[k] = f_17 * ski_473[k]
                   + f_8 * slh0_481[k]
                   - f_9 * slh1_481[k]
                   + f_3 * pc_y[k] * sli_641[k];

        t_825[k] = f_17 * ski_474[k]
                   + f_10 * slh0_482[k]
                   - f_11 * slh1_482[k]
                   + f_3 * pc_y[k] * sli_642[k];
    }

#pragma omp simd aligned(t_826, t_827, t_828, pc_x, pc_y, pc_z, ski_447, ski_475, ski_644, \
                         slh0_482, slh0_483, slh1_482, slh1_483, sli_643, \
                         sli_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_826[k] = f_17 * ski_475[k]
                   + f_3 * pc_y[k] * sli_643[k];

        t_827[k] = f_13 * ski_447[k]
                   + f_1 * slh0_482[k]
                   - f_2 * slh1_482[k]
                   + f_3 * pc_z[k] * sli_643[k];

        t_828[k] = f_14 * ski_644[k]
                   + f_1 * slh0_483[k]
                   - f_2 * slh1_483[k]
                   + f_3 * pc_x[k] * sli_644[k];
    }

#pragma omp simd aligned(t_829, t_830, t_831, t_832, pc_x, pc_y, pc_z, ski_448, ski_476, \
                         ski_478, ski_647, slh0_486, slh1_486, sli_644, sli_646, \
                         sli_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_829[k] = f_16 * ski_476[k]
                   + f_3 * pc_y[k] * sli_644[k];

        t_830[k] = f_14 * ski_448[k]
                   + f_3 * pc_z[k] * sli_644[k];

        t_831[k] = f_14 * ski_647[k]
                   + f_4 * slh0_486[k]
                   - f_5 * slh1_486[k]
                   + f_3 * pc_x[k] * sli_647[k];

        t_832[k] = f_16 * ski_478[k]
                   + f_3 * pc_y[k] * sli_646[k];
    }

#pragma omp simd aligned(t_833, t_834, t_835, pc_x, pc_z, ski_451, ski_649, ski_650, slh0_488, \
                         slh0_489, slh1_488, slh1_489, sli_647, sli_649, \
                         sli_650 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_833[k] = f_14 * ski_649[k]
                   + f_4 * slh0_488[k]
                   - f_5 * slh1_488[k]
                   + f_3 * pc_x[k] * sli_649[k];

        t_834[k] = f_14 * ski_650[k]
                   + f_6 * slh0_489[k]
                   - f_7 * slh1_489[k]
                   + f_3 * pc_x[k] * sli_650[k];

        t_835[k] = f_14 * ski_451[k]
                   + f_3 * pc_z[k] * sli_647[k];
    }

#pragma omp simd aligned(t_836, t_837, t_838, pc_x, pc_y, ski_481, ski_653, ski_654, slh0_492, \
                         slh0_493, slh1_492, slh1_493, sli_649, sli_653, \
                         sli_654 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_836[k] = f_16 * ski_481[k]
                   + f_3 * pc_y[k] * sli_649[k];

        t_837[k] = f_14 * ski_653[k]
                   + f_6 * slh0_492[k]
                   - f_7 * slh1_492[k]
                   + f_3 * pc_x[k] * sli_653[k];

        t_838[k] = f_14 * ski_654[k]
                   + f_8 * slh0_493[k]
                   - f_9 * slh1_493[k]
                   + f_3 * pc_x[k] * sli_654[k];
    }

#pragma omp simd aligned(t_839, t_840, t_841, pc_x, pc_y, pc_z, ski_454, ski_485, ski_656, \
                         slh0_495, slh1_495, sli_650, sli_653, \
                         sli_656 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_839[k] = f_14 * ski_454[k]
                   + f_3 * pc_z[k] * sli_650[k];

        t_840[k] = f_14 * ski_656[k]
                   + f_8 * slh0_495[k]
                   - f_9 * slh1_495[k]
                   + f_3 * pc_x[k] * sli_656[k];

        t_841[k] = f_16 * ski_485[k]
                   + f_3 * pc_y[k] * sli_653[k];
    }

#pragma omp simd aligned(t_842, t_843, t_844, pc_x, pc_z, ski_458, ski_658, ski_659, slh0_497, \
                         slh0_498, slh1_497, slh1_498, sli_654, sli_658, \
                         sli_659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_842[k] = f_14 * ski_658[k]
                   + f_8 * slh0_497[k]
                   - f_9 * slh1_497[k]
                   + f_3 * pc_x[k] * sli_658[k];

        t_843[k] = f_14 * ski_659[k]
                   + f_10 * slh0_498[k]
                   - f_11 * slh1_498[k]
                   + f_3 * pc_x[k] * sli_659[k];

        t_844[k] = f_14 * ski_458[k]
                   + f_3 * pc_z[k] * sli_654[k];
    }

#pragma omp simd aligned(t_845, t_846, t_847, pc_x, pc_y, ski_490, ski_661, ski_662, slh0_500, \
                         slh0_501, slh1_500, slh1_501, sli_658, sli_661, \
                         sli_662 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_845[k] = f_14 * ski_661[k]
                   + f_10 * slh0_500[k]
                   - f_11 * slh1_500[k]
                   + f_3 * pc_x[k] * sli_661[k];

        t_846[k] = f_14 * ski_662[k]
                   + f_10 * slh0_501[k]
                   - f_11 * slh1_501[k]
                   + f_3 * pc_x[k] * sli_662[k];

        t_847[k] = f_16 * ski_490[k]
                   + f_3 * pc_y[k] * sli_658[k];
    }

#pragma omp simd aligned(t_848, t_849, t_850, t_851, pc_x, ski_664, ski_665, ski_666, ski_667, \
                         slh0_503, slh1_503, sli_664, sli_665, sli_666, \
                         sli_667 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_848[k] = f_14 * ski_664[k]
                   + f_10 * slh0_503[k]
                   - f_11 * slh1_503[k]
                   + f_3 * pc_x[k] * sli_664[k];

        t_849[k] = f_14 * ski_665[k]
                   + f_3 * pc_x[k] * sli_665[k];

        t_850[k] = f_14 * ski_666[k]
                   + f_3 * pc_x[k] * sli_666[k];

        t_851[k] = f_14 * ski_667[k]
                   + f_3 * pc_x[k] * sli_667[k];
    }

#pragma omp simd aligned(t_852, t_853, t_854, t_855, pc_x, ski_668, ski_669, ski_670, ski_671, \
                         sli_668, sli_669, sli_670, sli_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_852[k] = f_14 * ski_668[k]
                   + f_3 * pc_x[k] * sli_668[k];

        t_853[k] = f_14 * ski_669[k]
                   + f_3 * pc_x[k] * sli_669[k];

        t_854[k] = f_14 * ski_670[k]
                   + f_3 * pc_x[k] * sli_670[k];

        t_855[k] = f_14 * ski_671[k]
                   + f_3 * pc_x[k] * sli_671[k];
    }

#pragma omp simd aligned(t_856, t_857, t_858, pc_y, pc_z, ski_469, ski_497, ski_499, slh0_498, \
                         slh0_500, slh1_498, slh1_500, sli_665, \
                         sli_667 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_856[k] = f_16 * ski_497[k]
                   + f_1 * slh0_498[k]
                   - f_2 * slh1_498[k]
                   + f_3 * pc_y[k] * sli_665[k];

        t_857[k] = f_14 * ski_469[k]
                   + f_3 * pc_z[k] * sli_665[k];

        t_858[k] = f_16 * ski_499[k]
                   + f_4 * slh0_500[k]
                   - f_5 * slh1_500[k]
                   + f_3 * pc_y[k] * sli_667[k];
    }

#pragma omp simd aligned(t_859, t_860, t_861, pc_y, ski_500, ski_501, ski_502, slh0_501, \
                         slh0_502, slh0_503, slh1_501, slh1_502, slh1_503, sli_668, sli_669, \
                         sli_670 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_859[k] = f_16 * ski_500[k]
                   + f_6 * slh0_501[k]
                   - f_7 * slh1_501[k]
                   + f_3 * pc_y[k] * sli_668[k];

        t_860[k] = f_16 * ski_501[k]
                   + f_8 * slh0_502[k]
                   - f_9 * slh1_502[k]
                   + f_3 * pc_y[k] * sli_669[k];

        t_861[k] = f_16 * ski_502[k]
                   + f_10 * slh0_503[k]
                   - f_11 * slh1_503[k]
                   + f_3 * pc_y[k] * sli_670[k];
    }

#pragma omp simd aligned(t_862, t_863, t_864, pc_x, pc_y, pc_z, ski_475, ski_503, ski_672, \
                         slh0_503, slh0_504, slh1_503, slh1_504, sli_671, \
                         sli_672 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_862[k] = f_16 * ski_503[k]
                   + f_3 * pc_y[k] * sli_671[k];

        t_863[k] = f_14 * ski_475[k]
                   + f_1 * slh0_503[k]
                   - f_2 * slh1_503[k]
                   + f_3 * pc_z[k] * sli_671[k];

        t_864[k] = f_14 * ski_672[k]
                   + f_1 * slh0_504[k]
                   - f_2 * slh1_504[k]
                   + f_3 * pc_x[k] * sli_672[k];
    }

#pragma omp simd aligned(t_865, t_866, t_867, t_868, pc_x, pc_y, pc_z, ski_476, ski_504, \
                         ski_506, ski_675, slh0_507, slh1_507, sli_672, sli_674, \
                         sli_675 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_865[k] = f_15 * ski_504[k]
                   + f_3 * pc_y[k] * sli_672[k];

        t_866[k] = f_15 * ski_476[k]
                   + f_3 * pc_z[k] * sli_672[k];

        t_867[k] = f_14 * ski_675[k]
                   + f_4 * slh0_507[k]
                   - f_5 * slh1_507[k]
                   + f_3 * pc_x[k] * sli_675[k];

        t_868[k] = f_15 * ski_506[k]
                   + f_3 * pc_y[k] * sli_674[k];
    }

#pragma omp simd aligned(t_869, t_870, t_871, pc_x, pc_z, ski_479, ski_677, ski_678, slh0_509, \
                         slh0_510, slh1_509, slh1_510, sli_675, sli_677, \
                         sli_678 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_869[k] = f_14 * ski_677[k]
                   + f_4 * slh0_509[k]
                   - f_5 * slh1_509[k]
                   + f_3 * pc_x[k] * sli_677[k];

        t_870[k] = f_14 * ski_678[k]
                   + f_6 * slh0_510[k]
                   - f_7 * slh1_510[k]
                   + f_3 * pc_x[k] * sli_678[k];

        t_871[k] = f_15 * ski_479[k]
                   + f_3 * pc_z[k] * sli_675[k];
    }

#pragma omp simd aligned(t_872, t_873, t_874, pc_x, pc_y, ski_509, ski_681, ski_682, slh0_513, \
                         slh0_514, slh1_513, slh1_514, sli_677, sli_681, \
                         sli_682 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_872[k] = f_15 * ski_509[k]
                   + f_3 * pc_y[k] * sli_677[k];

        t_873[k] = f_14 * ski_681[k]
                   + f_6 * slh0_513[k]
                   - f_7 * slh1_513[k]
                   + f_3 * pc_x[k] * sli_681[k];

        t_874[k] = f_14 * ski_682[k]
                   + f_8 * slh0_514[k]
                   - f_9 * slh1_514[k]
                   + f_3 * pc_x[k] * sli_682[k];
    }

#pragma omp simd aligned(t_875, t_876, t_877, pc_x, pc_y, pc_z, ski_482, ski_513, ski_684, \
                         slh0_516, slh1_516, sli_678, sli_681, \
                         sli_684 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_875[k] = f_15 * ski_482[k]
                   + f_3 * pc_z[k] * sli_678[k];

        t_876[k] = f_14 * ski_684[k]
                   + f_8 * slh0_516[k]
                   - f_9 * slh1_516[k]
                   + f_3 * pc_x[k] * sli_684[k];

        t_877[k] = f_15 * ski_513[k]
                   + f_3 * pc_y[k] * sli_681[k];
    }

#pragma omp simd aligned(t_878, t_879, t_880, pc_x, pc_z, ski_486, ski_686, ski_687, slh0_518, \
                         slh0_519, slh1_518, slh1_519, sli_682, sli_686, \
                         sli_687 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_878[k] = f_14 * ski_686[k]
                   + f_8 * slh0_518[k]
                   - f_9 * slh1_518[k]
                   + f_3 * pc_x[k] * sli_686[k];

        t_879[k] = f_14 * ski_687[k]
                   + f_10 * slh0_519[k]
                   - f_11 * slh1_519[k]
                   + f_3 * pc_x[k] * sli_687[k];

        t_880[k] = f_15 * ski_486[k]
                   + f_3 * pc_z[k] * sli_682[k];
    }

#pragma omp simd aligned(t_881, t_882, t_883, pc_x, pc_y, ski_518, ski_689, ski_690, slh0_521, \
                         slh0_522, slh1_521, slh1_522, sli_686, sli_689, \
                         sli_690 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_881[k] = f_14 * ski_689[k]
                   + f_10 * slh0_521[k]
                   - f_11 * slh1_521[k]
                   + f_3 * pc_x[k] * sli_689[k];

        t_882[k] = f_14 * ski_690[k]
                   + f_10 * slh0_522[k]
                   - f_11 * slh1_522[k]
                   + f_3 * pc_x[k] * sli_690[k];

        t_883[k] = f_15 * ski_518[k]
                   + f_3 * pc_y[k] * sli_686[k];
    }

#pragma omp simd aligned(t_884, t_885, t_886, t_887, pc_x, ski_692, ski_693, ski_694, ski_695, \
                         slh0_524, slh1_524, sli_692, sli_693, sli_694, \
                         sli_695 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_884[k] = f_14 * ski_692[k]
                   + f_10 * slh0_524[k]
                   - f_11 * slh1_524[k]
                   + f_3 * pc_x[k] * sli_692[k];

        t_885[k] = f_14 * ski_693[k]
                   + f_3 * pc_x[k] * sli_693[k];

        t_886[k] = f_14 * ski_694[k]
                   + f_3 * pc_x[k] * sli_694[k];

        t_887[k] = f_14 * ski_695[k]
                   + f_3 * pc_x[k] * sli_695[k];
    }

#pragma omp simd aligned(t_888, t_889, t_890, t_891, pc_x, ski_696, ski_697, ski_698, ski_699, \
                         sli_696, sli_697, sli_698, sli_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_888[k] = f_14 * ski_696[k]
                   + f_3 * pc_x[k] * sli_696[k];

        t_889[k] = f_14 * ski_697[k]
                   + f_3 * pc_x[k] * sli_697[k];

        t_890[k] = f_14 * ski_698[k]
                   + f_3 * pc_x[k] * sli_698[k];

        t_891[k] = f_14 * ski_699[k]
                   + f_3 * pc_x[k] * sli_699[k];
    }

#pragma omp simd aligned(t_892, t_893, t_894, pc_y, pc_z, ski_497, ski_525, ski_527, slh0_519, \
                         slh0_521, slh1_519, slh1_521, sli_693, \
                         sli_695 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_892[k] = f_15 * ski_525[k]
                   + f_1 * slh0_519[k]
                   - f_2 * slh1_519[k]
                   + f_3 * pc_y[k] * sli_693[k];

        t_893[k] = f_15 * ski_497[k]
                   + f_3 * pc_z[k] * sli_693[k];

        t_894[k] = f_15 * ski_527[k]
                   + f_4 * slh0_521[k]
                   - f_5 * slh1_521[k]
                   + f_3 * pc_y[k] * sli_695[k];
    }

#pragma omp simd aligned(t_895, t_896, t_897, pc_y, ski_528, ski_529, ski_530, slh0_522, \
                         slh0_523, slh0_524, slh1_522, slh1_523, slh1_524, sli_696, sli_697, \
                         sli_698 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_895[k] = f_15 * ski_528[k]
                   + f_6 * slh0_522[k]
                   - f_7 * slh1_522[k]
                   + f_3 * pc_y[k] * sli_696[k];

        t_896[k] = f_15 * ski_529[k]
                   + f_8 * slh0_523[k]
                   - f_9 * slh1_523[k]
                   + f_3 * pc_y[k] * sli_697[k];

        t_897[k] = f_15 * ski_530[k]
                   + f_10 * slh0_524[k]
                   - f_11 * slh1_524[k]
                   + f_3 * pc_y[k] * sli_698[k];
    }
}

static auto
compute_prim_slk_three_center_electron_repulsion_0_piece8(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t skk0,
                                                          const size_t ski, const size_t skk1,
                                                          const size_t slh0, const size_t slh1,
                                                          const size_t sli, const size_t ncols,
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
    const auto f_19 = 3.0 / q;

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

    const auto *skk0_720 = buffer.data(skk0 + 720);
    const auto *skk0_723 = buffer.data(skk0 + 723);
    const auto *skk0_725 = buffer.data(skk0 + 725);
    const auto *skk0_726 = buffer.data(skk0 + 726);
    const auto *skk0_729 = buffer.data(skk0 + 729);
    const auto *skk0_730 = buffer.data(skk0 + 730);
    const auto *skk0_732 = buffer.data(skk0 + 732);
    const auto *skk0_734 = buffer.data(skk0 + 734);
    const auto *skk0_735 = buffer.data(skk0 + 735);
    const auto *skk0_737 = buffer.data(skk0 + 737);
    const auto *skk0_738 = buffer.data(skk0 + 738);
    const auto *skk0_740 = buffer.data(skk0 + 740);
    const auto *skk0_755 = buffer.data(skk0 + 755);

    const auto *ski_503 = buffer.data(ski + 503);
    const auto *ski_504 = buffer.data(ski + 504);
    const auto *ski_507 = buffer.data(ski + 507);
    const auto *ski_510 = buffer.data(ski + 510);
    const auto *ski_514 = buffer.data(ski + 514);
    const auto *ski_525 = buffer.data(ski + 525);
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
    const auto *ski_555 = buffer.data(ski + 555);
    const auto *ski_556 = buffer.data(ski + 556);
    const auto *ski_557 = buffer.data(ski + 557);
    const auto *ski_558 = buffer.data(ski + 558);
    const auto *ski_559 = buffer.data(ski + 559);
    const auto *ski_560 = buffer.data(ski + 560);
    const auto *ski_561 = buffer.data(ski + 561);
    const auto *ski_562 = buffer.data(ski + 562);
    const auto *ski_563 = buffer.data(ski + 563);
    const auto *ski_565 = buffer.data(ski + 565);
    const auto *ski_566 = buffer.data(ski + 566);
    const auto *ski_568 = buffer.data(ski + 568);
    const auto *ski_569 = buffer.data(ski + 569);
    const auto *ski_570 = buffer.data(ski + 570);
    const auto *ski_572 = buffer.data(ski + 572);
    const auto *ski_573 = buffer.data(ski + 573);
    const auto *ski_574 = buffer.data(ski + 574);
    const auto *ski_581 = buffer.data(ski + 581);
    const auto *ski_583 = buffer.data(ski + 583);
    const auto *ski_584 = buffer.data(ski + 584);
    const auto *ski_585 = buffer.data(ski + 585);
    const auto *ski_586 = buffer.data(ski + 586);
    const auto *ski_587 = buffer.data(ski + 587);
    const auto *ski_700 = buffer.data(ski + 700);
    const auto *ski_703 = buffer.data(ski + 703);
    const auto *ski_705 = buffer.data(ski + 705);
    const auto *ski_706 = buffer.data(ski + 706);
    const auto *ski_709 = buffer.data(ski + 709);
    const auto *ski_710 = buffer.data(ski + 710);
    const auto *ski_712 = buffer.data(ski + 712);
    const auto *ski_714 = buffer.data(ski + 714);
    const auto *ski_715 = buffer.data(ski + 715);
    const auto *ski_717 = buffer.data(ski + 717);
    const auto *ski_718 = buffer.data(ski + 718);
    const auto *ski_720 = buffer.data(ski + 720);
    const auto *ski_721 = buffer.data(ski + 721);
    const auto *ski_722 = buffer.data(ski + 722);
    const auto *ski_723 = buffer.data(ski + 723);
    const auto *ski_724 = buffer.data(ski + 724);
    const auto *ski_725 = buffer.data(ski + 725);
    const auto *ski_726 = buffer.data(ski + 726);
    const auto *ski_727 = buffer.data(ski + 727);
    const auto *ski_749 = buffer.data(ski + 749);
    const auto *ski_750 = buffer.data(ski + 750);
    const auto *ski_751 = buffer.data(ski + 751);
    const auto *ski_752 = buffer.data(ski + 752);
    const auto *ski_753 = buffer.data(ski + 753);
    const auto *ski_754 = buffer.data(ski + 754);
    const auto *ski_755 = buffer.data(ski + 755);
    const auto *ski_756 = buffer.data(ski + 756);
    const auto *ski_759 = buffer.data(ski + 759);
    const auto *ski_761 = buffer.data(ski + 761);
    const auto *ski_762 = buffer.data(ski + 762);
    const auto *ski_765 = buffer.data(ski + 765);
    const auto *ski_766 = buffer.data(ski + 766);
    const auto *ski_768 = buffer.data(ski + 768);
    const auto *ski_770 = buffer.data(ski + 770);
    const auto *ski_771 = buffer.data(ski + 771);
    const auto *ski_773 = buffer.data(ski + 773);
    const auto *ski_774 = buffer.data(ski + 774);
    const auto *ski_776 = buffer.data(ski + 776);
    const auto *ski_777 = buffer.data(ski + 777);
    const auto *ski_778 = buffer.data(ski + 778);
    const auto *ski_779 = buffer.data(ski + 779);
    const auto *ski_780 = buffer.data(ski + 780);
    const auto *ski_781 = buffer.data(ski + 781);
    const auto *ski_782 = buffer.data(ski + 782);
    const auto *ski_783 = buffer.data(ski + 783);

    const auto *skk1_720 = buffer.data(skk1 + 720);
    const auto *skk1_723 = buffer.data(skk1 + 723);
    const auto *skk1_725 = buffer.data(skk1 + 725);
    const auto *skk1_726 = buffer.data(skk1 + 726);
    const auto *skk1_729 = buffer.data(skk1 + 729);
    const auto *skk1_730 = buffer.data(skk1 + 730);
    const auto *skk1_732 = buffer.data(skk1 + 732);
    const auto *skk1_734 = buffer.data(skk1 + 734);
    const auto *skk1_735 = buffer.data(skk1 + 735);
    const auto *skk1_737 = buffer.data(skk1 + 737);
    const auto *skk1_738 = buffer.data(skk1 + 738);
    const auto *skk1_740 = buffer.data(skk1 + 740);
    const auto *skk1_755 = buffer.data(skk1 + 755);

    const auto *slh0_524 = buffer.data(slh0 + 524);
    const auto *slh0_525 = buffer.data(slh0 + 525);
    const auto *slh0_528 = buffer.data(slh0 + 528);
    const auto *slh0_530 = buffer.data(slh0 + 530);
    const auto *slh0_531 = buffer.data(slh0 + 531);
    const auto *slh0_534 = buffer.data(slh0 + 534);
    const auto *slh0_535 = buffer.data(slh0 + 535);
    const auto *slh0_537 = buffer.data(slh0 + 537);
    const auto *slh0_539 = buffer.data(slh0 + 539);
    const auto *slh0_540 = buffer.data(slh0 + 540);
    const auto *slh0_542 = buffer.data(slh0 + 542);
    const auto *slh0_543 = buffer.data(slh0 + 543);
    const auto *slh0_544 = buffer.data(slh0 + 544);
    const auto *slh0_545 = buffer.data(slh0 + 545);
    const auto *slh0_561 = buffer.data(slh0 + 561);
    const auto *slh0_563 = buffer.data(slh0 + 563);
    const auto *slh0_564 = buffer.data(slh0 + 564);
    const auto *slh0_565 = buffer.data(slh0 + 565);
    const auto *slh0_566 = buffer.data(slh0 + 566);
    const auto *slh0_567 = buffer.data(slh0 + 567);
    const auto *slh0_570 = buffer.data(slh0 + 570);
    const auto *slh0_572 = buffer.data(slh0 + 572);
    const auto *slh0_573 = buffer.data(slh0 + 573);
    const auto *slh0_576 = buffer.data(slh0 + 576);
    const auto *slh0_577 = buffer.data(slh0 + 577);
    const auto *slh0_579 = buffer.data(slh0 + 579);
    const auto *slh0_581 = buffer.data(slh0 + 581);
    const auto *slh0_582 = buffer.data(slh0 + 582);
    const auto *slh0_584 = buffer.data(slh0 + 584);
    const auto *slh0_585 = buffer.data(slh0 + 585);
    const auto *slh0_586 = buffer.data(slh0 + 586);
    const auto *slh0_587 = buffer.data(slh0 + 587);

    const auto *slh1_524 = buffer.data(slh1 + 524);
    const auto *slh1_525 = buffer.data(slh1 + 525);
    const auto *slh1_528 = buffer.data(slh1 + 528);
    const auto *slh1_530 = buffer.data(slh1 + 530);
    const auto *slh1_531 = buffer.data(slh1 + 531);
    const auto *slh1_534 = buffer.data(slh1 + 534);
    const auto *slh1_535 = buffer.data(slh1 + 535);
    const auto *slh1_537 = buffer.data(slh1 + 537);
    const auto *slh1_539 = buffer.data(slh1 + 539);
    const auto *slh1_540 = buffer.data(slh1 + 540);
    const auto *slh1_542 = buffer.data(slh1 + 542);
    const auto *slh1_543 = buffer.data(slh1 + 543);
    const auto *slh1_544 = buffer.data(slh1 + 544);
    const auto *slh1_545 = buffer.data(slh1 + 545);
    const auto *slh1_561 = buffer.data(slh1 + 561);
    const auto *slh1_563 = buffer.data(slh1 + 563);
    const auto *slh1_564 = buffer.data(slh1 + 564);
    const auto *slh1_565 = buffer.data(slh1 + 565);
    const auto *slh1_566 = buffer.data(slh1 + 566);
    const auto *slh1_567 = buffer.data(slh1 + 567);
    const auto *slh1_570 = buffer.data(slh1 + 570);
    const auto *slh1_572 = buffer.data(slh1 + 572);
    const auto *slh1_573 = buffer.data(slh1 + 573);
    const auto *slh1_576 = buffer.data(slh1 + 576);
    const auto *slh1_577 = buffer.data(slh1 + 577);
    const auto *slh1_579 = buffer.data(slh1 + 579);
    const auto *slh1_581 = buffer.data(slh1 + 581);
    const auto *slh1_582 = buffer.data(slh1 + 582);
    const auto *slh1_584 = buffer.data(slh1 + 584);
    const auto *slh1_585 = buffer.data(slh1 + 585);
    const auto *slh1_586 = buffer.data(slh1 + 586);
    const auto *slh1_587 = buffer.data(slh1 + 587);

    const auto *sli_699 = buffer.data(sli + 699);
    const auto *sli_700 = buffer.data(sli + 700);
    const auto *sli_702 = buffer.data(sli + 702);
    const auto *sli_703 = buffer.data(sli + 703);
    const auto *sli_705 = buffer.data(sli + 705);
    const auto *sli_706 = buffer.data(sli + 706);
    const auto *sli_709 = buffer.data(sli + 709);
    const auto *sli_710 = buffer.data(sli + 710);
    const auto *sli_712 = buffer.data(sli + 712);
    const auto *sli_714 = buffer.data(sli + 714);
    const auto *sli_715 = buffer.data(sli + 715);
    const auto *sli_717 = buffer.data(sli + 717);
    const auto *sli_718 = buffer.data(sli + 718);
    const auto *sli_720 = buffer.data(sli + 720);
    const auto *sli_721 = buffer.data(sli + 721);
    const auto *sli_722 = buffer.data(sli + 722);
    const auto *sli_723 = buffer.data(sli + 723);
    const auto *sli_724 = buffer.data(sli + 724);
    const auto *sli_725 = buffer.data(sli + 725);
    const auto *sli_726 = buffer.data(sli + 726);
    const auto *sli_727 = buffer.data(sli + 727);
    const auto *sli_728 = buffer.data(sli + 728);
    const auto *sli_730 = buffer.data(sli + 730);
    const auto *sli_731 = buffer.data(sli + 731);
    const auto *sli_733 = buffer.data(sli + 733);
    const auto *sli_734 = buffer.data(sli + 734);
    const auto *sli_737 = buffer.data(sli + 737);
    const auto *sli_738 = buffer.data(sli + 738);
    const auto *sli_742 = buffer.data(sli + 742);
    const auto *sli_749 = buffer.data(sli + 749);
    const auto *sli_750 = buffer.data(sli + 750);
    const auto *sli_751 = buffer.data(sli + 751);
    const auto *sli_752 = buffer.data(sli + 752);
    const auto *sli_753 = buffer.data(sli + 753);
    const auto *sli_754 = buffer.data(sli + 754);
    const auto *sli_755 = buffer.data(sli + 755);
    const auto *sli_756 = buffer.data(sli + 756);
    const auto *sli_758 = buffer.data(sli + 758);
    const auto *sli_759 = buffer.data(sli + 759);
    const auto *sli_761 = buffer.data(sli + 761);
    const auto *sli_762 = buffer.data(sli + 762);
    const auto *sli_765 = buffer.data(sli + 765);
    const auto *sli_766 = buffer.data(sli + 766);
    const auto *sli_768 = buffer.data(sli + 768);
    const auto *sli_770 = buffer.data(sli + 770);
    const auto *sli_771 = buffer.data(sli + 771);
    const auto *sli_773 = buffer.data(sli + 773);
    const auto *sli_774 = buffer.data(sli + 774);
    const auto *sli_776 = buffer.data(sli + 776);
    const auto *sli_777 = buffer.data(sli + 777);
    const auto *sli_778 = buffer.data(sli + 778);
    const auto *sli_779 = buffer.data(sli + 779);
    const auto *sli_780 = buffer.data(sli + 780);
    const auto *sli_781 = buffer.data(sli + 781);
    const auto *sli_782 = buffer.data(sli + 782);
    const auto *sli_783 = buffer.data(sli + 783);

#pragma omp simd aligned(t_898, t_899, t_900, pc_x, pc_y, pc_z, ski_503, ski_531, ski_700, \
                         slh0_524, slh0_525, slh1_524, slh1_525, sli_699, \
                         sli_700 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_898[k] = f_15 * ski_531[k]
                   + f_3 * pc_y[k] * sli_699[k];

        t_899[k] = f_15 * ski_503[k]
                   + f_1 * slh0_524[k]
                   - f_2 * slh1_524[k]
                   + f_3 * pc_z[k] * sli_699[k];

        t_900[k] = f_14 * ski_700[k]
                   + f_1 * slh0_525[k]
                   - f_2 * slh1_525[k]
                   + f_3 * pc_x[k] * sli_700[k];
    }

#pragma omp simd aligned(t_901, t_902, t_903, t_904, pc_x, pc_y, pc_z, ski_504, ski_532, \
                         ski_534, ski_703, slh0_528, slh1_528, sli_700, sli_702, \
                         sli_703 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_901[k] = f_14 * ski_532[k]
                   + f_3 * pc_y[k] * sli_700[k];

        t_902[k] = f_16 * ski_504[k]
                   + f_3 * pc_z[k] * sli_700[k];

        t_903[k] = f_14 * ski_703[k]
                   + f_4 * slh0_528[k]
                   - f_5 * slh1_528[k]
                   + f_3 * pc_x[k] * sli_703[k];

        t_904[k] = f_14 * ski_534[k]
                   + f_3 * pc_y[k] * sli_702[k];
    }

#pragma omp simd aligned(t_905, t_906, t_907, pc_x, pc_z, ski_507, ski_705, ski_706, slh0_530, \
                         slh0_531, slh1_530, slh1_531, sli_703, sli_705, \
                         sli_706 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_905[k] = f_14 * ski_705[k]
                   + f_4 * slh0_530[k]
                   - f_5 * slh1_530[k]
                   + f_3 * pc_x[k] * sli_705[k];

        t_906[k] = f_14 * ski_706[k]
                   + f_6 * slh0_531[k]
                   - f_7 * slh1_531[k]
                   + f_3 * pc_x[k] * sli_706[k];

        t_907[k] = f_16 * ski_507[k]
                   + f_3 * pc_z[k] * sli_703[k];
    }

#pragma omp simd aligned(t_908, t_909, t_910, pc_x, pc_y, ski_537, ski_709, ski_710, slh0_534, \
                         slh0_535, slh1_534, slh1_535, sli_705, sli_709, \
                         sli_710 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_908[k] = f_14 * ski_537[k]
                   + f_3 * pc_y[k] * sli_705[k];

        t_909[k] = f_14 * ski_709[k]
                   + f_6 * slh0_534[k]
                   - f_7 * slh1_534[k]
                   + f_3 * pc_x[k] * sli_709[k];

        t_910[k] = f_14 * ski_710[k]
                   + f_8 * slh0_535[k]
                   - f_9 * slh1_535[k]
                   + f_3 * pc_x[k] * sli_710[k];
    }

#pragma omp simd aligned(t_911, t_912, t_913, pc_x, pc_y, pc_z, ski_510, ski_541, ski_712, \
                         slh0_537, slh1_537, sli_706, sli_709, \
                         sli_712 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_911[k] = f_16 * ski_510[k]
                   + f_3 * pc_z[k] * sli_706[k];

        t_912[k] = f_14 * ski_712[k]
                   + f_8 * slh0_537[k]
                   - f_9 * slh1_537[k]
                   + f_3 * pc_x[k] * sli_712[k];

        t_913[k] = f_14 * ski_541[k]
                   + f_3 * pc_y[k] * sli_709[k];
    }

#pragma omp simd aligned(t_914, t_915, t_916, pc_x, pc_z, ski_514, ski_714, ski_715, slh0_539, \
                         slh0_540, slh1_539, slh1_540, sli_710, sli_714, \
                         sli_715 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_914[k] = f_14 * ski_714[k]
                   + f_8 * slh0_539[k]
                   - f_9 * slh1_539[k]
                   + f_3 * pc_x[k] * sli_714[k];

        t_915[k] = f_14 * ski_715[k]
                   + f_10 * slh0_540[k]
                   - f_11 * slh1_540[k]
                   + f_3 * pc_x[k] * sli_715[k];

        t_916[k] = f_16 * ski_514[k]
                   + f_3 * pc_z[k] * sli_710[k];
    }

#pragma omp simd aligned(t_917, t_918, t_919, pc_x, pc_y, ski_546, ski_717, ski_718, slh0_542, \
                         slh0_543, slh1_542, slh1_543, sli_714, sli_717, \
                         sli_718 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_917[k] = f_14 * ski_717[k]
                   + f_10 * slh0_542[k]
                   - f_11 * slh1_542[k]
                   + f_3 * pc_x[k] * sli_717[k];

        t_918[k] = f_14 * ski_718[k]
                   + f_10 * slh0_543[k]
                   - f_11 * slh1_543[k]
                   + f_3 * pc_x[k] * sli_718[k];

        t_919[k] = f_14 * ski_546[k]
                   + f_3 * pc_y[k] * sli_714[k];
    }

#pragma omp simd aligned(t_920, t_921, t_922, t_923, pc_x, ski_720, ski_721, ski_722, ski_723, \
                         slh0_545, slh1_545, sli_720, sli_721, sli_722, \
                         sli_723 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_920[k] = f_14 * ski_720[k]
                   + f_10 * slh0_545[k]
                   - f_11 * slh1_545[k]
                   + f_3 * pc_x[k] * sli_720[k];

        t_921[k] = f_14 * ski_721[k]
                   + f_3 * pc_x[k] * sli_721[k];

        t_922[k] = f_14 * ski_722[k]
                   + f_3 * pc_x[k] * sli_722[k];

        t_923[k] = f_14 * ski_723[k]
                   + f_3 * pc_x[k] * sli_723[k];
    }

#pragma omp simd aligned(t_924, t_925, t_926, t_927, pc_x, ski_724, ski_725, ski_726, ski_727, \
                         sli_724, sli_725, sli_726, sli_727 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_924[k] = f_14 * ski_724[k]
                   + f_3 * pc_x[k] * sli_724[k];

        t_925[k] = f_14 * ski_725[k]
                   + f_3 * pc_x[k] * sli_725[k];

        t_926[k] = f_14 * ski_726[k]
                   + f_3 * pc_x[k] * sli_726[k];

        t_927[k] = f_14 * ski_727[k]
                   + f_3 * pc_x[k] * sli_727[k];
    }

#pragma omp simd aligned(t_928, t_929, t_930, pc_y, pc_z, ski_525, ski_553, ski_555, slh0_540, \
                         slh0_542, slh1_540, slh1_542, sli_721, \
                         sli_723 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_928[k] = f_14 * ski_553[k]
                   + f_1 * slh0_540[k]
                   - f_2 * slh1_540[k]
                   + f_3 * pc_y[k] * sli_721[k];

        t_929[k] = f_16 * ski_525[k]
                   + f_3 * pc_z[k] * sli_721[k];

        t_930[k] = f_14 * ski_555[k]
                   + f_4 * slh0_542[k]
                   - f_5 * slh1_542[k]
                   + f_3 * pc_y[k] * sli_723[k];
    }

#pragma omp simd aligned(t_931, t_932, t_933, pc_y, ski_556, ski_557, ski_558, slh0_543, \
                         slh0_544, slh0_545, slh1_543, slh1_544, slh1_545, sli_724, sli_725, \
                         sli_726 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_931[k] = f_14 * ski_556[k]
                   + f_6 * slh0_543[k]
                   - f_7 * slh1_543[k]
                   + f_3 * pc_y[k] * sli_724[k];

        t_932[k] = f_14 * ski_557[k]
                   + f_8 * slh0_544[k]
                   - f_9 * slh1_544[k]
                   + f_3 * pc_y[k] * sli_725[k];

        t_933[k] = f_14 * ski_558[k]
                   + f_10 * slh0_545[k]
                   - f_11 * slh1_545[k]
                   + f_3 * pc_y[k] * sli_726[k];
    }

#pragma omp simd aligned(t_934, t_935, t_936, t_937, pb_y, pc_y, pc_z, skk0_720, ski_531, \
                         ski_559, ski_560, skk1_720, slh0_545, slh1_545, sli_727, \
                         sli_728 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_934[k] = f_14 * ski_559[k]
                   + f_3 * pc_y[k] * sli_727[k];

        t_935[k] = f_16 * ski_531[k]
                   + f_1 * slh0_545[k]
                   - f_2 * slh1_545[k]
                   + f_3 * pc_z[k] * sli_727[k];

        t_936[k] = pb_y[k] * skk0_720[k]
                   - f_12 * pc_y[k] * skk1_720[k];

        t_937[k] = f_13 * ski_560[k]
                   + f_3 * pc_y[k] * sli_728[k];
    }

#pragma omp simd aligned(t_938, t_939, t_940, t_941, pb_y, pc_y, pc_z, skk0_723, skk0_725, \
                         ski_532, ski_561, ski_562, skk1_723, skk1_725, sli_728, \
                         sli_730 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_938[k] = f_17 * ski_532[k]
                   + f_3 * pc_z[k] * sli_728[k];

        t_939[k] = pb_y[k] * skk0_723[k]
                   + f_14 * ski_561[k]
                   - f_12 * pc_y[k] * skk1_723[k];

        t_940[k] = f_13 * ski_562[k]
                   + f_3 * pc_y[k] * sli_730[k];

        t_941[k] = pb_y[k] * skk0_725[k]
                   - f_12 * pc_y[k] * skk1_725[k];
    }

#pragma omp simd aligned(t_942, t_943, t_944, t_945, pb_y, pc_y, pc_z, skk0_726, skk0_729, \
                         ski_535, ski_563, ski_565, skk1_726, skk1_729, sli_731, \
                         sli_733 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_942[k] = pb_y[k] * skk0_726[k]
                   + f_15 * ski_563[k]
                   - f_12 * pc_y[k] * skk1_726[k];

        t_943[k] = f_17 * ski_535[k]
                   + f_3 * pc_z[k] * sli_731[k];

        t_944[k] = f_13 * ski_565[k]
                   + f_3 * pc_y[k] * sli_733[k];

        t_945[k] = pb_y[k] * skk0_729[k]
                   - f_12 * pc_y[k] * skk1_729[k];
    }

#pragma omp simd aligned(t_946, t_947, t_948, pb_y, pc_y, pc_z, skk0_730, skk0_732, ski_538, \
                         ski_566, ski_568, skk1_730, skk1_732, \
                         sli_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_946[k] = pb_y[k] * skk0_730[k]
                   + f_16 * ski_566[k]
                   - f_12 * pc_y[k] * skk1_730[k];

        t_947[k] = f_17 * ski_538[k]
                   + f_3 * pc_z[k] * sli_734[k];

        t_948[k] = pb_y[k] * skk0_732[k]
                   + f_14 * ski_568[k]
                   - f_12 * pc_y[k] * skk1_732[k];
    }

#pragma omp simd aligned(t_949, t_950, t_951, t_952, pb_y, pc_y, pc_z, skk0_734, skk0_735, \
                         ski_542, ski_569, ski_570, skk1_734, skk1_735, sli_737, \
                         sli_738 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_949[k] = f_13 * ski_569[k]
                   + f_3 * pc_y[k] * sli_737[k];

        t_950[k] = pb_y[k] * skk0_734[k]
                   - f_12 * pc_y[k] * skk1_734[k];

        t_951[k] = pb_y[k] * skk0_735[k]
                   + f_17 * ski_570[k]
                   - f_12 * pc_y[k] * skk1_735[k];

        t_952[k] = f_17 * ski_542[k]
                   + f_3 * pc_z[k] * sli_738[k];
    }

#pragma omp simd aligned(t_953, t_954, t_955, t_956, pb_y, pc_y, skk0_737, skk0_738, skk0_740, \
                         ski_572, ski_573, ski_574, skk1_737, skk1_738, skk1_740, \
                         sli_742 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_953[k] = pb_y[k] * skk0_737[k]
                   + f_15 * ski_572[k]
                   - f_12 * pc_y[k] * skk1_737[k];

        t_954[k] = pb_y[k] * skk0_738[k]
                   + f_14 * ski_573[k]
                   - f_12 * pc_y[k] * skk1_738[k];

        t_955[k] = f_13 * ski_574[k]
                   + f_3 * pc_y[k] * sli_742[k];

        t_956[k] = pb_y[k] * skk0_740[k]
                   - f_12 * pc_y[k] * skk1_740[k];
    }

#pragma omp simd aligned(t_957, t_958, t_959, t_960, t_961, pc_x, ski_749, ski_750, ski_751, \
                         ski_752, ski_753, sli_749, sli_750, sli_751, sli_752, \
                         sli_753 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_957[k] = f_14 * ski_749[k]
                   + f_3 * pc_x[k] * sli_749[k];

        t_958[k] = f_14 * ski_750[k]
                   + f_3 * pc_x[k] * sli_750[k];

        t_959[k] = f_14 * ski_751[k]
                   + f_3 * pc_x[k] * sli_751[k];

        t_960[k] = f_14 * ski_752[k]
                   + f_3 * pc_x[k] * sli_752[k];

        t_961[k] = f_14 * ski_753[k]
                   + f_3 * pc_x[k] * sli_753[k];
    }

#pragma omp simd aligned(t_962, t_963, t_964, t_965, pc_x, pc_y, pc_z, ski_553, ski_581, \
                         ski_754, ski_755, slh0_561, slh1_561, sli_749, sli_754, \
                         sli_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_962[k] = f_14 * ski_754[k]
                   + f_3 * pc_x[k] * sli_754[k];

        t_963[k] = f_14 * ski_755[k]
                   + f_3 * pc_x[k] * sli_755[k];

        t_964[k] = f_13 * ski_581[k]
                   + f_1 * slh0_561[k]
                   - f_2 * slh1_561[k]
                   + f_3 * pc_y[k] * sli_749[k];

        t_965[k] = f_17 * ski_553[k]
                   + f_3 * pc_z[k] * sli_749[k];
    }

#pragma omp simd aligned(t_966, t_967, t_968, pc_y, ski_583, ski_584, ski_585, slh0_563, \
                         slh0_564, slh0_565, slh1_563, slh1_564, slh1_565, sli_751, sli_752, \
                         sli_753 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_966[k] = f_13 * ski_583[k]
                   + f_4 * slh0_563[k]
                   - f_5 * slh1_563[k]
                   + f_3 * pc_y[k] * sli_751[k];

        t_967[k] = f_13 * ski_584[k]
                   + f_6 * slh0_564[k]
                   - f_7 * slh1_564[k]
                   + f_3 * pc_y[k] * sli_752[k];

        t_968[k] = f_13 * ski_585[k]
                   + f_8 * slh0_565[k]
                   - f_9 * slh1_565[k]
                   + f_3 * pc_y[k] * sli_753[k];
    }

#pragma omp simd aligned(t_969, t_970, t_971, pb_y, pc_y, skk0_755, ski_586, ski_587, \
                         skk1_755, slh0_566, slh1_566, sli_754, \
                         sli_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_969[k] = f_13 * ski_586[k]
                   + f_10 * slh0_566[k]
                   - f_11 * slh1_566[k]
                   + f_3 * pc_y[k] * sli_754[k];

        t_970[k] = f_13 * ski_587[k]
                   + f_3 * pc_y[k] * sli_755[k];

        t_971[k] = pb_y[k] * skk0_755[k]
                   - f_12 * pc_y[k] * skk1_755[k];
    }

#pragma omp simd aligned(t_972, t_973, t_974, t_975, pc_x, pc_y, pc_z, ski_560, ski_756, \
                         ski_759, slh0_567, slh0_570, slh1_567, slh1_570, sli_756, \
                         sli_759 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_972[k] = f_14 * ski_756[k]
                   + f_1 * slh0_567[k]
                   - f_2 * slh1_567[k]
                   + f_3 * pc_x[k] * sli_756[k];

        t_973[k] = f_3 * pc_y[k] * sli_756[k];

        t_974[k] = f_19 * ski_560[k]
                   + f_3 * pc_z[k] * sli_756[k];

        t_975[k] = f_14 * ski_759[k]
                   + f_4 * slh0_570[k]
                   - f_5 * slh1_570[k]
                   + f_3 * pc_x[k] * sli_759[k];
    }

#pragma omp simd aligned(t_976, t_977, t_978, pc_x, pc_y, ski_761, ski_762, slh0_572, \
                         slh0_573, slh1_572, slh1_573, sli_758, sli_761, \
                         sli_762 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_976[k] = f_3 * pc_y[k] * sli_758[k];

        t_977[k] = f_14 * ski_761[k]
                   + f_4 * slh0_572[k]
                   - f_5 * slh1_572[k]
                   + f_3 * pc_x[k] * sli_761[k];

        t_978[k] = f_14 * ski_762[k]
                   + f_6 * slh0_573[k]
                   - f_7 * slh1_573[k]
                   + f_3 * pc_x[k] * sli_762[k];
    }

#pragma omp simd aligned(t_979, t_980, t_981, pc_x, pc_y, pc_z, ski_563, ski_765, slh0_576, \
                         slh1_576, sli_759, sli_761, sli_765 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_979[k] = f_19 * ski_563[k]
                   + f_3 * pc_z[k] * sli_759[k];

        t_980[k] = f_3 * pc_y[k] * sli_761[k];

        t_981[k] = f_14 * ski_765[k]
                   + f_6 * slh0_576[k]
                   - f_7 * slh1_576[k]
                   + f_3 * pc_x[k] * sli_765[k];
    }

#pragma omp simd aligned(t_982, t_983, t_984, pc_x, pc_z, ski_566, ski_766, ski_768, slh0_577, \
                         slh0_579, slh1_577, slh1_579, sli_762, sli_766, \
                         sli_768 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_982[k] = f_14 * ski_766[k]
                   + f_8 * slh0_577[k]
                   - f_9 * slh1_577[k]
                   + f_3 * pc_x[k] * sli_766[k];

        t_983[k] = f_19 * ski_566[k]
                   + f_3 * pc_z[k] * sli_762[k];

        t_984[k] = f_14 * ski_768[k]
                   + f_8 * slh0_579[k]
                   - f_9 * slh1_579[k]
                   + f_3 * pc_x[k] * sli_768[k];
    }

#pragma omp simd aligned(t_985, t_986, t_987, pc_x, pc_y, ski_770, ski_771, slh0_581, \
                         slh0_582, slh1_581, slh1_582, sli_765, sli_770, \
                         sli_771 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_985[k] = f_3 * pc_y[k] * sli_765[k];

        t_986[k] = f_14 * ski_770[k]
                   + f_8 * slh0_581[k]
                   - f_9 * slh1_581[k]
                   + f_3 * pc_x[k] * sli_770[k];

        t_987[k] = f_14 * ski_771[k]
                   + f_10 * slh0_582[k]
                   - f_11 * slh1_582[k]
                   + f_3 * pc_x[k] * sli_771[k];
    }

#pragma omp simd aligned(t_988, t_989, t_990, pc_x, pc_z, ski_570, ski_773, ski_774, slh0_584, \
                         slh0_585, slh1_584, slh1_585, sli_766, sli_773, \
                         sli_774 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_988[k] = f_19 * ski_570[k]
                   + f_3 * pc_z[k] * sli_766[k];

        t_989[k] = f_14 * ski_773[k]
                   + f_10 * slh0_584[k]
                   - f_11 * slh1_584[k]
                   + f_3 * pc_x[k] * sli_773[k];

        t_990[k] = f_14 * ski_774[k]
                   + f_10 * slh0_585[k]
                   - f_11 * slh1_585[k]
                   + f_3 * pc_x[k] * sli_774[k];
    }

#pragma omp simd aligned(t_991, t_992, t_993, t_994, pc_x, pc_y, ski_776, ski_777, ski_778, \
                         slh0_587, slh1_587, sli_770, sli_776, sli_777, \
                         sli_778 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_991[k] = f_3 * pc_y[k] * sli_770[k];

        t_992[k] = f_14 * ski_776[k]
                   + f_10 * slh0_587[k]
                   - f_11 * slh1_587[k]
                   + f_3 * pc_x[k] * sli_776[k];

        t_993[k] = f_14 * ski_777[k]
                   + f_3 * pc_x[k] * sli_777[k];

        t_994[k] = f_14 * ski_778[k]
                   + f_3 * pc_x[k] * sli_778[k];
    }

#pragma omp simd aligned(t_995, t_996, t_997, t_998, t_999, pc_x, ski_779, ski_780, ski_781, \
                         ski_782, ski_783, sli_779, sli_780, sli_781, sli_782, \
                         sli_783 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_995[k] = f_14 * ski_779[k]
                   + f_3 * pc_x[k] * sli_779[k];

        t_996[k] = f_14 * ski_780[k]
                   + f_3 * pc_x[k] * sli_780[k];

        t_997[k] = f_14 * ski_781[k]
                   + f_3 * pc_x[k] * sli_781[k];

        t_998[k] = f_14 * ski_782[k]
                   + f_3 * pc_x[k] * sli_782[k];

        t_999[k] = f_14 * ski_783[k]
                   + f_3 * pc_x[k] * sli_783[k];
    }

#pragma omp simd aligned(t_1000, t_1001, t_1002, t_1003, pc_y, pc_z, ski_581, slh0_582, \
                         slh0_584, slh0_585, slh1_582, slh1_584, slh1_585, sli_777, sli_779, \
                         sli_780 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1000[k] = f_1 * slh0_582[k]
                    - f_2 * slh1_582[k]
                    + f_3 * pc_y[k] * sli_777[k];

        t_1001[k] = f_19 * ski_581[k]
                    + f_3 * pc_z[k] * sli_777[k];

        t_1002[k] = f_4 * slh0_584[k]
                    - f_5 * slh1_584[k]
                    + f_3 * pc_y[k] * sli_779[k];

        t_1003[k] = f_6 * slh0_585[k]
                    - f_7 * slh1_585[k]
                    + f_3 * pc_y[k] * sli_780[k];
    }

#pragma omp simd aligned(t_1004, t_1005, t_1006, t_1007, pc_y, pc_z, ski_587, slh0_586, \
                         slh0_587, slh1_586, slh1_587, sli_781, sli_782, \
                         sli_783 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1004[k] = f_8 * slh0_586[k]
                    - f_9 * slh1_586[k]
                    + f_3 * pc_y[k] * sli_781[k];

        t_1005[k] = f_10 * slh0_587[k]
                    - f_11 * slh1_587[k]
                    + f_3 * pc_y[k] * sli_782[k];

        t_1006[k] = f_3 * pc_y[k] * sli_783[k];

        t_1007[k] = f_19 * ski_587[k]
                    + f_1 * slh0_587[k]
                    - f_2 * slh1_587[k]
                    + f_3 * pc_z[k] * sli_783[k];
    }
}

static auto
compute_prim_slk_three_center_electron_repulsion_0_piece9(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t skk0,
                                                          const size_t ski, const size_t skk1,
                                                          const size_t sli, const size_t ncols,
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
    const auto f_19 = 3.0 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *skk0_756 = buffer.data(skk0 + 756);
    const auto *skk0_759 = buffer.data(skk0 + 759);
    const auto *skk0_762 = buffer.data(skk0 + 762);
    const auto *skk0_766 = buffer.data(skk0 + 766);
    const auto *skk0_771 = buffer.data(skk0 + 771);
    const auto *skk0_1008 = buffer.data(skk0 + 1008);
    const auto *skk0_1011 = buffer.data(skk0 + 1011);
    const auto *skk0_1013 = buffer.data(skk0 + 1013);
    const auto *skk0_1014 = buffer.data(skk0 + 1014);
    const auto *skk0_1017 = buffer.data(skk0 + 1017);
    const auto *skk0_1018 = buffer.data(skk0 + 1018);
    const auto *skk0_1020 = buffer.data(skk0 + 1020);
    const auto *skk0_1022 = buffer.data(skk0 + 1022);
    const auto *skk0_1023 = buffer.data(skk0 + 1023);
    const auto *skk0_1025 = buffer.data(skk0 + 1025);
    const auto *skk0_1026 = buffer.data(skk0 + 1026);
    const auto *skk0_1028 = buffer.data(skk0 + 1028);
    const auto *skk0_1036 = buffer.data(skk0 + 1036);
    const auto *skk0_1038 = buffer.data(skk0 + 1038);
    const auto *skk0_1039 = buffer.data(skk0 + 1039);
    const auto *skk0_1040 = buffer.data(skk0 + 1040);
    const auto *skk0_1041 = buffer.data(skk0 + 1041);
    const auto *skk0_1043 = buffer.data(skk0 + 1043);
    const auto *skk0_1049 = buffer.data(skk0 + 1049);
    const auto *skk0_1053 = buffer.data(skk0 + 1053);
    const auto *skk0_1056 = buffer.data(skk0 + 1056);
    const auto *skk0_1058 = buffer.data(skk0 + 1058);
    const auto *skk0_1061 = buffer.data(skk0 + 1061);
    const auto *skk0_1062 = buffer.data(skk0 + 1062);
    const auto *skk0_1064 = buffer.data(skk0 + 1064);
    const auto *skk0_1072 = buffer.data(skk0 + 1072);
    const auto *skk0_1074 = buffer.data(skk0 + 1074);
    const auto *skk0_1075 = buffer.data(skk0 + 1075);
    const auto *skk0_1076 = buffer.data(skk0 + 1076);
    const auto *skk0_1077 = buffer.data(skk0 + 1077);
    const auto *skk0_1079 = buffer.data(skk0 + 1079);
    const auto *skk0_1080 = buffer.data(skk0 + 1080);
    const auto *skk0_1083 = buffer.data(skk0 + 1083);
    const auto *skk0_1085 = buffer.data(skk0 + 1085);
    const auto *skk0_1086 = buffer.data(skk0 + 1086);
    const auto *skk0_1089 = buffer.data(skk0 + 1089);
    const auto *skk0_1090 = buffer.data(skk0 + 1090);
    const auto *skk0_1092 = buffer.data(skk0 + 1092);
    const auto *skk0_1094 = buffer.data(skk0 + 1094);
    const auto *skk0_1095 = buffer.data(skk0 + 1095);
    const auto *skk0_1097 = buffer.data(skk0 + 1097);
    const auto *skk0_1098 = buffer.data(skk0 + 1098);
    const auto *skk0_1100 = buffer.data(skk0 + 1100);
    const auto *skk0_1108 = buffer.data(skk0 + 1108);
    const auto *skk0_1110 = buffer.data(skk0 + 1110);
    const auto *skk0_1111 = buffer.data(skk0 + 1111);
    const auto *skk0_1112 = buffer.data(skk0 + 1112);
    const auto *skk0_1113 = buffer.data(skk0 + 1113);
    const auto *skk0_1115 = buffer.data(skk0 + 1115);
    const auto *skk0_1116 = buffer.data(skk0 + 1116);
    const auto *skk0_1119 = buffer.data(skk0 + 1119);
    const auto *skk0_1121 = buffer.data(skk0 + 1121);
    const auto *skk0_1122 = buffer.data(skk0 + 1122);
    const auto *skk0_1125 = buffer.data(skk0 + 1125);
    const auto *skk0_1126 = buffer.data(skk0 + 1126);

    const auto *ski_588 = buffer.data(ski + 588);
    const auto *ski_590 = buffer.data(ski + 590);
    const auto *ski_591 = buffer.data(ski + 591);
    const auto *ski_593 = buffer.data(ski + 593);
    const auto *ski_594 = buffer.data(ski + 594);
    const auto *ski_597 = buffer.data(ski + 597);
    const auto *ski_598 = buffer.data(ski + 598);
    const auto *ski_602 = buffer.data(ski + 602);
    const auto *ski_609 = buffer.data(ski + 609);
    const auto *ski_615 = buffer.data(ski + 615);
    const auto *ski_616 = buffer.data(ski + 616);
    const auto *ski_618 = buffer.data(ski + 618);
    const auto *ski_619 = buffer.data(ski + 619);
    const auto *ski_621 = buffer.data(ski + 621);
    const auto *ski_622 = buffer.data(ski + 622);
    const auto *ski_625 = buffer.data(ski + 625);
    const auto *ski_626 = buffer.data(ski + 626);
    const auto *ski_630 = buffer.data(ski + 630);
    const auto *ski_637 = buffer.data(ski + 637);
    const auto *ski_643 = buffer.data(ski + 643);
    const auto *ski_644 = buffer.data(ski + 644);
    const auto *ski_646 = buffer.data(ski + 646);
    const auto *ski_647 = buffer.data(ski + 647);
    const auto *ski_649 = buffer.data(ski + 649);
    const auto *ski_653 = buffer.data(ski + 653);
    const auto *ski_658 = buffer.data(ski + 658);
    const auto *ski_671 = buffer.data(ski + 671);
    const auto *ski_672 = buffer.data(ski + 672);
    const auto *ski_674 = buffer.data(ski + 674);
    const auto *ski_677 = buffer.data(ski + 677);
    const auto *ski_784 = buffer.data(ski + 784);
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
    const auto *ski_817 = buffer.data(ski + 817);
    const auto *ski_821 = buffer.data(ski + 821);
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
    const auto *ski_871 = buffer.data(ski + 871);
    const auto *ski_873 = buffer.data(ski + 873);
    const auto *ski_874 = buffer.data(ski + 874);
    const auto *ski_877 = buffer.data(ski + 877);
    const auto *ski_878 = buffer.data(ski + 878);

    const auto *skk1_756 = buffer.data(skk1 + 756);
    const auto *skk1_759 = buffer.data(skk1 + 759);
    const auto *skk1_762 = buffer.data(skk1 + 762);
    const auto *skk1_766 = buffer.data(skk1 + 766);
    const auto *skk1_771 = buffer.data(skk1 + 771);
    const auto *skk1_1008 = buffer.data(skk1 + 1008);
    const auto *skk1_1011 = buffer.data(skk1 + 1011);
    const auto *skk1_1013 = buffer.data(skk1 + 1013);
    const auto *skk1_1014 = buffer.data(skk1 + 1014);
    const auto *skk1_1017 = buffer.data(skk1 + 1017);
    const auto *skk1_1018 = buffer.data(skk1 + 1018);
    const auto *skk1_1020 = buffer.data(skk1 + 1020);
    const auto *skk1_1022 = buffer.data(skk1 + 1022);
    const auto *skk1_1023 = buffer.data(skk1 + 1023);
    const auto *skk1_1025 = buffer.data(skk1 + 1025);
    const auto *skk1_1026 = buffer.data(skk1 + 1026);
    const auto *skk1_1028 = buffer.data(skk1 + 1028);
    const auto *skk1_1036 = buffer.data(skk1 + 1036);
    const auto *skk1_1038 = buffer.data(skk1 + 1038);
    const auto *skk1_1039 = buffer.data(skk1 + 1039);
    const auto *skk1_1040 = buffer.data(skk1 + 1040);
    const auto *skk1_1041 = buffer.data(skk1 + 1041);
    const auto *skk1_1043 = buffer.data(skk1 + 1043);
    const auto *skk1_1049 = buffer.data(skk1 + 1049);
    const auto *skk1_1053 = buffer.data(skk1 + 1053);
    const auto *skk1_1056 = buffer.data(skk1 + 1056);
    const auto *skk1_1058 = buffer.data(skk1 + 1058);
    const auto *skk1_1061 = buffer.data(skk1 + 1061);
    const auto *skk1_1062 = buffer.data(skk1 + 1062);
    const auto *skk1_1064 = buffer.data(skk1 + 1064);
    const auto *skk1_1072 = buffer.data(skk1 + 1072);
    const auto *skk1_1074 = buffer.data(skk1 + 1074);
    const auto *skk1_1075 = buffer.data(skk1 + 1075);
    const auto *skk1_1076 = buffer.data(skk1 + 1076);
    const auto *skk1_1077 = buffer.data(skk1 + 1077);
    const auto *skk1_1079 = buffer.data(skk1 + 1079);
    const auto *skk1_1080 = buffer.data(skk1 + 1080);
    const auto *skk1_1083 = buffer.data(skk1 + 1083);
    const auto *skk1_1085 = buffer.data(skk1 + 1085);
    const auto *skk1_1086 = buffer.data(skk1 + 1086);
    const auto *skk1_1089 = buffer.data(skk1 + 1089);
    const auto *skk1_1090 = buffer.data(skk1 + 1090);
    const auto *skk1_1092 = buffer.data(skk1 + 1092);
    const auto *skk1_1094 = buffer.data(skk1 + 1094);
    const auto *skk1_1095 = buffer.data(skk1 + 1095);
    const auto *skk1_1097 = buffer.data(skk1 + 1097);
    const auto *skk1_1098 = buffer.data(skk1 + 1098);
    const auto *skk1_1100 = buffer.data(skk1 + 1100);
    const auto *skk1_1108 = buffer.data(skk1 + 1108);
    const auto *skk1_1110 = buffer.data(skk1 + 1110);
    const auto *skk1_1111 = buffer.data(skk1 + 1111);
    const auto *skk1_1112 = buffer.data(skk1 + 1112);
    const auto *skk1_1113 = buffer.data(skk1 + 1113);
    const auto *skk1_1115 = buffer.data(skk1 + 1115);
    const auto *skk1_1116 = buffer.data(skk1 + 1116);
    const auto *skk1_1119 = buffer.data(skk1 + 1119);
    const auto *skk1_1121 = buffer.data(skk1 + 1121);
    const auto *skk1_1122 = buffer.data(skk1 + 1122);
    const auto *skk1_1125 = buffer.data(skk1 + 1125);
    const auto *skk1_1126 = buffer.data(skk1 + 1126);

    const auto *sli_784 = buffer.data(sli + 784);
    const auto *sli_786 = buffer.data(sli + 786);
    const auto *sli_787 = buffer.data(sli + 787);
    const auto *sli_789 = buffer.data(sli + 789);
    const auto *sli_790 = buffer.data(sli + 790);
    const auto *sli_793 = buffer.data(sli + 793);
    const auto *sli_794 = buffer.data(sli + 794);
    const auto *sli_798 = buffer.data(sli + 798);
    const auto *sli_805 = buffer.data(sli + 805);
    const auto *sli_806 = buffer.data(sli + 806);
    const auto *sli_807 = buffer.data(sli + 807);
    const auto *sli_808 = buffer.data(sli + 808);
    const auto *sli_809 = buffer.data(sli + 809);
    const auto *sli_810 = buffer.data(sli + 810);
    const auto *sli_811 = buffer.data(sli + 811);
    const auto *sli_812 = buffer.data(sli + 812);
    const auto *sli_814 = buffer.data(sli + 814);
    const auto *sli_815 = buffer.data(sli + 815);
    const auto *sli_817 = buffer.data(sli + 817);
    const auto *sli_818 = buffer.data(sli + 818);
    const auto *sli_821 = buffer.data(sli + 821);
    const auto *sli_822 = buffer.data(sli + 822);
    const auto *sli_826 = buffer.data(sli + 826);
    const auto *sli_833 = buffer.data(sli + 833);
    const auto *sli_834 = buffer.data(sli + 834);
    const auto *sli_835 = buffer.data(sli + 835);
    const auto *sli_836 = buffer.data(sli + 836);
    const auto *sli_837 = buffer.data(sli + 837);
    const auto *sli_838 = buffer.data(sli + 838);
    const auto *sli_839 = buffer.data(sli + 839);
    const auto *sli_840 = buffer.data(sli + 840);
    const auto *sli_842 = buffer.data(sli + 842);
    const auto *sli_843 = buffer.data(sli + 843);
    const auto *sli_845 = buffer.data(sli + 845);
    const auto *sli_846 = buffer.data(sli + 846);
    const auto *sli_849 = buffer.data(sli + 849);
    const auto *sli_850 = buffer.data(sli + 850);
    const auto *sli_854 = buffer.data(sli + 854);
    const auto *sli_861 = buffer.data(sli + 861);
    const auto *sli_862 = buffer.data(sli + 862);
    const auto *sli_863 = buffer.data(sli + 863);
    const auto *sli_864 = buffer.data(sli + 864);
    const auto *sli_865 = buffer.data(sli + 865);
    const auto *sli_866 = buffer.data(sli + 866);
    const auto *sli_867 = buffer.data(sli + 867);
    const auto *sli_868 = buffer.data(sli + 868);
    const auto *sli_870 = buffer.data(sli + 870);
    const auto *sli_871 = buffer.data(sli + 871);
    const auto *sli_873 = buffer.data(sli + 873);

#pragma omp simd aligned(t_1008, t_1009, t_1010, t_1011, pb_x, pc_x, pc_y, pc_z, skk0_1008, \
                         skk0_1011, ski_588, ski_784, ski_787, skk1_1008, skk1_1011, \
                         sli_784 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1008[k] = pb_x[k] * skk0_1008[k]
                    + f_18 * ski_784[k]
                    - f_12 * pc_x[k] * skk1_1008[k];

        t_1009[k] = f_18 * ski_588[k]
                    + f_3 * pc_y[k] * sli_784[k];

        t_1010[k] = f_3 * pc_z[k] * sli_784[k];

        t_1011[k] = pb_x[k] * skk0_1011[k]
                    + f_17 * ski_787[k]
                    - f_12 * pc_x[k] * skk1_1011[k];
    }

#pragma omp simd aligned(t_1012, t_1013, t_1014, pb_x, pc_x, pc_y, skk0_1013, skk0_1014, \
                         ski_590, ski_789, ski_790, skk1_1013, skk1_1014, \
                         sli_786 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1012[k] = f_18 * ski_590[k]
                    + f_3 * pc_y[k] * sli_786[k];

        t_1013[k] = pb_x[k] * skk0_1013[k]
                    + f_17 * ski_789[k]
                    - f_12 * pc_x[k] * skk1_1013[k];

        t_1014[k] = pb_x[k] * skk0_1014[k]
                    + f_16 * ski_790[k]
                    - f_12 * pc_x[k] * skk1_1014[k];
    }

#pragma omp simd aligned(t_1015, t_1016, t_1017, pb_x, pc_x, pc_y, pc_z, skk0_1017, ski_593, \
                         ski_793, skk1_1017, sli_787, sli_789 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1015[k] = f_3 * pc_z[k] * sli_787[k];

        t_1016[k] = f_18 * ski_593[k]
                    + f_3 * pc_y[k] * sli_789[k];

        t_1017[k] = pb_x[k] * skk0_1017[k]
                    + f_16 * ski_793[k]
                    - f_12 * pc_x[k] * skk1_1017[k];
    }

#pragma omp simd aligned(t_1018, t_1019, t_1020, pb_x, pc_x, pc_z, skk0_1018, skk0_1020, \
                         ski_794, ski_796, skk1_1018, skk1_1020, \
                         sli_790 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1018[k] = pb_x[k] * skk0_1018[k]
                    + f_15 * ski_794[k]
                    - f_12 * pc_x[k] * skk1_1018[k];

        t_1019[k] = f_3 * pc_z[k] * sli_790[k];

        t_1020[k] = pb_x[k] * skk0_1020[k]
                    + f_15 * ski_796[k]
                    - f_12 * pc_x[k] * skk1_1020[k];
    }

#pragma omp simd aligned(t_1021, t_1022, t_1023, pb_x, pc_x, pc_y, skk0_1022, skk0_1023, \
                         ski_597, ski_798, ski_799, skk1_1022, skk1_1023, \
                         sli_793 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1021[k] = f_18 * ski_597[k]
                    + f_3 * pc_y[k] * sli_793[k];

        t_1022[k] = pb_x[k] * skk0_1022[k]
                    + f_15 * ski_798[k]
                    - f_12 * pc_x[k] * skk1_1022[k];

        t_1023[k] = pb_x[k] * skk0_1023[k]
                    + f_14 * ski_799[k]
                    - f_12 * pc_x[k] * skk1_1023[k];
    }

#pragma omp simd aligned(t_1024, t_1025, t_1026, pb_x, pc_x, pc_z, skk0_1025, skk0_1026, \
                         ski_801, ski_802, skk1_1025, skk1_1026, \
                         sli_794 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1024[k] = f_3 * pc_z[k] * sli_794[k];

        t_1025[k] = pb_x[k] * skk0_1025[k]
                    + f_14 * ski_801[k]
                    - f_12 * pc_x[k] * skk1_1025[k];

        t_1026[k] = pb_x[k] * skk0_1026[k]
                    + f_14 * ski_802[k]
                    - f_12 * pc_x[k] * skk1_1026[k];
    }

#pragma omp simd aligned(t_1027, t_1028, t_1029, t_1030, pb_x, pc_x, pc_y, skk0_1028, ski_602, \
                         ski_804, ski_805, ski_806, skk1_1028, sli_798, sli_805, \
                         sli_806 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1027[k] = f_18 * ski_602[k]
                    + f_3 * pc_y[k] * sli_798[k];

        t_1028[k] = pb_x[k] * skk0_1028[k]
                    + f_14 * ski_804[k]
                    - f_12 * pc_x[k] * skk1_1028[k];

        t_1029[k] = f_13 * ski_805[k]
                    + f_3 * pc_x[k] * sli_805[k];

        t_1030[k] = f_13 * ski_806[k]
                    + f_3 * pc_x[k] * sli_806[k];
    }

#pragma omp simd aligned(t_1031, t_1032, t_1033, t_1034, t_1035, pc_x, ski_807, ski_808, \
                         ski_809, ski_810, ski_811, sli_807, sli_808, sli_809, sli_810, \
                         sli_811 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1031[k] = f_13 * ski_807[k]
                    + f_3 * pc_x[k] * sli_807[k];

        t_1032[k] = f_13 * ski_808[k]
                    + f_3 * pc_x[k] * sli_808[k];

        t_1033[k] = f_13 * ski_809[k]
                    + f_3 * pc_x[k] * sli_809[k];

        t_1034[k] = f_13 * ski_810[k]
                    + f_3 * pc_x[k] * sli_810[k];

        t_1035[k] = f_13 * ski_811[k]
                    + f_3 * pc_x[k] * sli_811[k];
    }

#pragma omp simd aligned(t_1036, t_1037, t_1038, t_1039, pb_x, pc_x, pc_z, skk0_1036, \
                         skk0_1038, skk0_1039, skk1_1036, skk1_1038, skk1_1039, \
                         sli_805 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1036[k] = pb_x[k] * skk0_1036[k]
                    - f_12 * pc_x[k] * skk1_1036[k];

        t_1037[k] = f_3 * pc_z[k] * sli_805[k];

        t_1038[k] = pb_x[k] * skk0_1038[k]
                    - f_12 * pc_x[k] * skk1_1038[k];

        t_1039[k] = pb_x[k] * skk0_1039[k]
                    - f_12 * pc_x[k] * skk1_1039[k];
    }

#pragma omp simd aligned(t_1040, t_1041, t_1042, t_1043, pb_x, pc_x, pc_y, skk0_1040, \
                         skk0_1041, skk0_1043, ski_615, skk1_1040, skk1_1041, skk1_1043, \
                         sli_811 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1040[k] = pb_x[k] * skk0_1040[k]
                    - f_12 * pc_x[k] * skk1_1040[k];

        t_1041[k] = pb_x[k] * skk0_1041[k]
                    - f_12 * pc_x[k] * skk1_1041[k];

        t_1042[k] = f_18 * ski_615[k]
                    + f_3 * pc_y[k] * sli_811[k];

        t_1043[k] = pb_x[k] * skk0_1043[k]
                    - f_12 * pc_x[k] * skk1_1043[k];
    }

#pragma omp simd aligned(t_1044, t_1045, t_1046, t_1047, pb_z, pc_y, pc_z, skk0_756, skk0_759, \
                         ski_588, ski_616, skk1_756, skk1_759, \
                         sli_812 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1044[k] = pb_z[k] * skk0_756[k]
                    - f_12 * pc_z[k] * skk1_756[k];

        t_1045[k] = f_19 * ski_616[k]
                    + f_3 * pc_y[k] * sli_812[k];

        t_1046[k] = f_13 * ski_588[k]
                    + f_3 * pc_z[k] * sli_812[k];

        t_1047[k] = pb_z[k] * skk0_759[k]
                    - f_12 * pc_z[k] * skk1_759[k];
    }

#pragma omp simd aligned(t_1048, t_1049, t_1050, pb_x, pb_z, pc_x, pc_y, pc_z, skk0_762, \
                         skk0_1049, ski_618, ski_817, skk1_762, skk1_1049, \
                         sli_814 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1048[k] = f_19 * ski_618[k]
                    + f_3 * pc_y[k] * sli_814[k];

        t_1049[k] = pb_x[k] * skk0_1049[k]
                    + f_17 * ski_817[k]
                    - f_12 * pc_x[k] * skk1_1049[k];

        t_1050[k] = pb_z[k] * skk0_762[k]
                    - f_12 * pc_z[k] * skk1_762[k];
    }

#pragma omp simd aligned(t_1051, t_1052, t_1053, pb_x, pc_x, pc_y, pc_z, skk0_1053, ski_591, \
                         ski_621, ski_821, skk1_1053, sli_815, \
                         sli_817 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1051[k] = f_13 * ski_591[k]
                    + f_3 * pc_z[k] * sli_815[k];

        t_1052[k] = f_19 * ski_621[k]
                    + f_3 * pc_y[k] * sli_817[k];

        t_1053[k] = pb_x[k] * skk0_1053[k]
                    + f_16 * ski_821[k]
                    - f_12 * pc_x[k] * skk1_1053[k];
    }

#pragma omp simd aligned(t_1054, t_1055, t_1056, pb_x, pb_z, pc_x, pc_z, skk0_766, skk0_1056, \
                         ski_594, ski_824, skk1_766, skk1_1056, \
                         sli_818 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1054[k] = pb_z[k] * skk0_766[k]
                    - f_12 * pc_z[k] * skk1_766[k];

        t_1055[k] = f_13 * ski_594[k]
                    + f_3 * pc_z[k] * sli_818[k];

        t_1056[k] = pb_x[k] * skk0_1056[k]
                    + f_15 * ski_824[k]
                    - f_12 * pc_x[k] * skk1_1056[k];
    }

#pragma omp simd aligned(t_1057, t_1058, t_1059, pb_x, pb_z, pc_x, pc_y, pc_z, skk0_771, \
                         skk0_1058, ski_625, ski_826, skk1_771, skk1_1058, \
                         sli_821 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1057[k] = f_19 * ski_625[k]
                    + f_3 * pc_y[k] * sli_821[k];

        t_1058[k] = pb_x[k] * skk0_1058[k]
                    + f_15 * ski_826[k]
                    - f_12 * pc_x[k] * skk1_1058[k];

        t_1059[k] = pb_z[k] * skk0_771[k]
                    - f_12 * pc_z[k] * skk1_771[k];
    }

#pragma omp simd aligned(t_1060, t_1061, t_1062, pb_x, pc_x, pc_z, skk0_1061, skk0_1062, \
                         ski_598, ski_829, ski_830, skk1_1061, skk1_1062, \
                         sli_822 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1060[k] = f_13 * ski_598[k]
                    + f_3 * pc_z[k] * sli_822[k];

        t_1061[k] = pb_x[k] * skk0_1061[k]
                    + f_14 * ski_829[k]
                    - f_12 * pc_x[k] * skk1_1061[k];

        t_1062[k] = pb_x[k] * skk0_1062[k]
                    + f_14 * ski_830[k]
                    - f_12 * pc_x[k] * skk1_1062[k];
    }

#pragma omp simd aligned(t_1063, t_1064, t_1065, t_1066, pb_x, pc_x, pc_y, skk0_1064, ski_630, \
                         ski_832, ski_833, ski_834, skk1_1064, sli_826, sli_833, \
                         sli_834 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1063[k] = f_19 * ski_630[k]
                    + f_3 * pc_y[k] * sli_826[k];

        t_1064[k] = pb_x[k] * skk0_1064[k]
                    + f_14 * ski_832[k]
                    - f_12 * pc_x[k] * skk1_1064[k];

        t_1065[k] = f_13 * ski_833[k]
                    + f_3 * pc_x[k] * sli_833[k];

        t_1066[k] = f_13 * ski_834[k]
                    + f_3 * pc_x[k] * sli_834[k];
    }

#pragma omp simd aligned(t_1067, t_1068, t_1069, t_1070, t_1071, pc_x, ski_835, ski_836, \
                         ski_837, ski_838, ski_839, sli_835, sli_836, sli_837, sli_838, \
                         sli_839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1067[k] = f_13 * ski_835[k]
                    + f_3 * pc_x[k] * sli_835[k];

        t_1068[k] = f_13 * ski_836[k]
                    + f_3 * pc_x[k] * sli_836[k];

        t_1069[k] = f_13 * ski_837[k]
                    + f_3 * pc_x[k] * sli_837[k];

        t_1070[k] = f_13 * ski_838[k]
                    + f_3 * pc_x[k] * sli_838[k];

        t_1071[k] = f_13 * ski_839[k]
                    + f_3 * pc_x[k] * sli_839[k];
    }

#pragma omp simd aligned(t_1072, t_1073, t_1074, t_1075, pb_x, pc_x, pc_z, skk0_1072, \
                         skk0_1074, skk0_1075, ski_609, skk1_1072, skk1_1074, skk1_1075, \
                         sli_833 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1072[k] = pb_x[k] * skk0_1072[k]
                    - f_12 * pc_x[k] * skk1_1072[k];

        t_1073[k] = f_13 * ski_609[k]
                    + f_3 * pc_z[k] * sli_833[k];

        t_1074[k] = pb_x[k] * skk0_1074[k]
                    - f_12 * pc_x[k] * skk1_1074[k];

        t_1075[k] = pb_x[k] * skk0_1075[k]
                    - f_12 * pc_x[k] * skk1_1075[k];
    }

#pragma omp simd aligned(t_1076, t_1077, t_1078, t_1079, pb_x, pc_x, pc_y, skk0_1076, \
                         skk0_1077, skk0_1079, ski_643, skk1_1076, skk1_1077, skk1_1079, \
                         sli_839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1076[k] = pb_x[k] * skk0_1076[k]
                    - f_12 * pc_x[k] * skk1_1076[k];

        t_1077[k] = pb_x[k] * skk0_1077[k]
                    - f_12 * pc_x[k] * skk1_1077[k];

        t_1078[k] = f_19 * ski_643[k]
                    + f_3 * pc_y[k] * sli_839[k];

        t_1079[k] = pb_x[k] * skk0_1079[k]
                    - f_12 * pc_x[k] * skk1_1079[k];
    }

#pragma omp simd aligned(t_1080, t_1081, t_1082, pb_x, pc_x, pc_y, pc_z, skk0_1080, ski_616, \
                         ski_644, ski_840, skk1_1080, sli_840 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1080[k] = pb_x[k] * skk0_1080[k]
                    + f_18 * ski_840[k]
                    - f_12 * pc_x[k] * skk1_1080[k];

        t_1081[k] = f_17 * ski_644[k]
                    + f_3 * pc_y[k] * sli_840[k];

        t_1082[k] = f_14 * ski_616[k]
                    + f_3 * pc_z[k] * sli_840[k];
    }

#pragma omp simd aligned(t_1083, t_1084, t_1085, pb_x, pc_x, pc_y, skk0_1083, skk0_1085, \
                         ski_646, ski_843, ski_845, skk1_1083, skk1_1085, \
                         sli_842 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1083[k] = pb_x[k] * skk0_1083[k]
                    + f_17 * ski_843[k]
                    - f_12 * pc_x[k] * skk1_1083[k];

        t_1084[k] = f_17 * ski_646[k]
                    + f_3 * pc_y[k] * sli_842[k];

        t_1085[k] = pb_x[k] * skk0_1085[k]
                    + f_17 * ski_845[k]
                    - f_12 * pc_x[k] * skk1_1085[k];
    }

#pragma omp simd aligned(t_1086, t_1087, t_1088, pb_x, pc_x, pc_y, pc_z, skk0_1086, ski_619, \
                         ski_649, ski_846, skk1_1086, sli_843, \
                         sli_845 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1086[k] = pb_x[k] * skk0_1086[k]
                    + f_16 * ski_846[k]
                    - f_12 * pc_x[k] * skk1_1086[k];

        t_1087[k] = f_14 * ski_619[k]
                    + f_3 * pc_z[k] * sli_843[k];

        t_1088[k] = f_17 * ski_649[k]
                    + f_3 * pc_y[k] * sli_845[k];
    }

#pragma omp simd aligned(t_1089, t_1090, t_1091, pb_x, pc_x, pc_z, skk0_1089, skk0_1090, \
                         ski_622, ski_849, ski_850, skk1_1089, skk1_1090, \
                         sli_846 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1089[k] = pb_x[k] * skk0_1089[k]
                    + f_16 * ski_849[k]
                    - f_12 * pc_x[k] * skk1_1089[k];

        t_1090[k] = pb_x[k] * skk0_1090[k]
                    + f_15 * ski_850[k]
                    - f_12 * pc_x[k] * skk1_1090[k];

        t_1091[k] = f_14 * ski_622[k]
                    + f_3 * pc_z[k] * sli_846[k];
    }

#pragma omp simd aligned(t_1092, t_1093, t_1094, pb_x, pc_x, pc_y, skk0_1092, skk0_1094, \
                         ski_653, ski_852, ski_854, skk1_1092, skk1_1094, \
                         sli_849 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1092[k] = pb_x[k] * skk0_1092[k]
                    + f_15 * ski_852[k]
                    - f_12 * pc_x[k] * skk1_1092[k];

        t_1093[k] = f_17 * ski_653[k]
                    + f_3 * pc_y[k] * sli_849[k];

        t_1094[k] = pb_x[k] * skk0_1094[k]
                    + f_15 * ski_854[k]
                    - f_12 * pc_x[k] * skk1_1094[k];
    }

#pragma omp simd aligned(t_1095, t_1096, t_1097, pb_x, pc_x, pc_z, skk0_1095, skk0_1097, \
                         ski_626, ski_855, ski_857, skk1_1095, skk1_1097, \
                         sli_850 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1095[k] = pb_x[k] * skk0_1095[k]
                    + f_14 * ski_855[k]
                    - f_12 * pc_x[k] * skk1_1095[k];

        t_1096[k] = f_14 * ski_626[k]
                    + f_3 * pc_z[k] * sli_850[k];

        t_1097[k] = pb_x[k] * skk0_1097[k]
                    + f_14 * ski_857[k]
                    - f_12 * pc_x[k] * skk1_1097[k];
    }

#pragma omp simd aligned(t_1098, t_1099, t_1100, pb_x, pc_x, pc_y, skk0_1098, skk0_1100, \
                         ski_658, ski_858, ski_860, skk1_1098, skk1_1100, \
                         sli_854 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1098[k] = pb_x[k] * skk0_1098[k]
                    + f_14 * ski_858[k]
                    - f_12 * pc_x[k] * skk1_1098[k];

        t_1099[k] = f_17 * ski_658[k]
                    + f_3 * pc_y[k] * sli_854[k];

        t_1100[k] = pb_x[k] * skk0_1100[k]
                    + f_14 * ski_860[k]
                    - f_12 * pc_x[k] * skk1_1100[k];
    }

#pragma omp simd aligned(t_1101, t_1102, t_1103, t_1104, t_1105, pc_x, ski_861, ski_862, \
                         ski_863, ski_864, ski_865, sli_861, sli_862, sli_863, sli_864, \
                         sli_865 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1101[k] = f_13 * ski_861[k]
                    + f_3 * pc_x[k] * sli_861[k];

        t_1102[k] = f_13 * ski_862[k]
                    + f_3 * pc_x[k] * sli_862[k];

        t_1103[k] = f_13 * ski_863[k]
                    + f_3 * pc_x[k] * sli_863[k];

        t_1104[k] = f_13 * ski_864[k]
                    + f_3 * pc_x[k] * sli_864[k];

        t_1105[k] = f_13 * ski_865[k]
                    + f_3 * pc_x[k] * sli_865[k];
    }

#pragma omp simd aligned(t_1106, t_1107, t_1108, t_1109, pb_x, pc_x, pc_z, skk0_1108, ski_637, \
                         ski_866, ski_867, skk1_1108, sli_861, sli_866, \
                         sli_867 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1106[k] = f_13 * ski_866[k]
                    + f_3 * pc_x[k] * sli_866[k];

        t_1107[k] = f_13 * ski_867[k]
                    + f_3 * pc_x[k] * sli_867[k];

        t_1108[k] = pb_x[k] * skk0_1108[k]
                    - f_12 * pc_x[k] * skk1_1108[k];

        t_1109[k] = f_14 * ski_637[k]
                    + f_3 * pc_z[k] * sli_861[k];
    }

#pragma omp simd aligned(t_1110, t_1111, t_1112, t_1113, pb_x, pc_x, skk0_1110, skk0_1111, \
                         skk0_1112, skk0_1113, skk1_1110, skk1_1111, skk1_1112, \
                         skk1_1113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1110[k] = pb_x[k] * skk0_1110[k]
                    - f_12 * pc_x[k] * skk1_1110[k];

        t_1111[k] = pb_x[k] * skk0_1111[k]
                    - f_12 * pc_x[k] * skk1_1111[k];

        t_1112[k] = pb_x[k] * skk0_1112[k]
                    - f_12 * pc_x[k] * skk1_1112[k];

        t_1113[k] = pb_x[k] * skk0_1113[k]
                    - f_12 * pc_x[k] * skk1_1113[k];
    }

#pragma omp simd aligned(t_1114, t_1115, t_1116, t_1117, pb_x, pc_x, pc_y, skk0_1115, \
                         skk0_1116, ski_671, ski_672, ski_868, skk1_1115, skk1_1116, sli_867, \
                         sli_868 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1114[k] = f_17 * ski_671[k]
                    + f_3 * pc_y[k] * sli_867[k];

        t_1115[k] = pb_x[k] * skk0_1115[k]
                    - f_12 * pc_x[k] * skk1_1115[k];

        t_1116[k] = pb_x[k] * skk0_1116[k]
                    + f_18 * ski_868[k]
                    - f_12 * pc_x[k] * skk1_1116[k];

        t_1117[k] = f_16 * ski_672[k]
                    + f_3 * pc_y[k] * sli_868[k];
    }

#pragma omp simd aligned(t_1118, t_1119, t_1120, pb_x, pc_x, pc_y, pc_z, skk0_1119, ski_644, \
                         ski_674, ski_871, skk1_1119, sli_868, \
                         sli_870 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1118[k] = f_15 * ski_644[k]
                    + f_3 * pc_z[k] * sli_868[k];

        t_1119[k] = pb_x[k] * skk0_1119[k]
                    + f_17 * ski_871[k]
                    - f_12 * pc_x[k] * skk1_1119[k];

        t_1120[k] = f_16 * ski_674[k]
                    + f_3 * pc_y[k] * sli_870[k];
    }

#pragma omp simd aligned(t_1121, t_1122, t_1123, pb_x, pc_x, pc_z, skk0_1121, skk0_1122, \
                         ski_647, ski_873, ski_874, skk1_1121, skk1_1122, \
                         sli_871 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1121[k] = pb_x[k] * skk0_1121[k]
                    + f_17 * ski_873[k]
                    - f_12 * pc_x[k] * skk1_1121[k];

        t_1122[k] = pb_x[k] * skk0_1122[k]
                    + f_16 * ski_874[k]
                    - f_12 * pc_x[k] * skk1_1122[k];

        t_1123[k] = f_15 * ski_647[k]
                    + f_3 * pc_z[k] * sli_871[k];
    }

#pragma omp simd aligned(t_1124, t_1125, t_1126, pb_x, pc_x, pc_y, skk0_1125, skk0_1126, \
                         ski_677, ski_877, ski_878, skk1_1125, skk1_1126, \
                         sli_873 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1124[k] = f_16 * ski_677[k]
                    + f_3 * pc_y[k] * sli_873[k];

        t_1125[k] = pb_x[k] * skk0_1125[k]
                    + f_16 * ski_877[k]
                    - f_12 * pc_x[k] * skk1_1125[k];

        t_1126[k] = pb_x[k] * skk0_1126[k]
                    + f_15 * ski_878[k]
                    - f_12 * pc_x[k] * skk1_1126[k];
    }
}

static auto
compute_prim_slk_three_center_electron_repulsion_0_piece10(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t skk0,
                                                           const size_t ski, const size_t skk1,
                                                           const size_t sli, const size_t ncols,
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
    const auto f_19 = 3.0 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *skk0_972 = buffer.data(skk0 + 972);
    const auto *skk0_977 = buffer.data(skk0 + 977);
    const auto *skk0_981 = buffer.data(skk0 + 981);
    const auto *skk0_986 = buffer.data(skk0 + 986);
    const auto *skk0_992 = buffer.data(skk0 + 992);
    const auto *skk0_1128 = buffer.data(skk0 + 1128);
    const auto *skk0_1130 = buffer.data(skk0 + 1130);
    const auto *skk0_1131 = buffer.data(skk0 + 1131);
    const auto *skk0_1133 = buffer.data(skk0 + 1133);
    const auto *skk0_1134 = buffer.data(skk0 + 1134);
    const auto *skk0_1136 = buffer.data(skk0 + 1136);
    const auto *skk0_1144 = buffer.data(skk0 + 1144);
    const auto *skk0_1146 = buffer.data(skk0 + 1146);
    const auto *skk0_1147 = buffer.data(skk0 + 1147);
    const auto *skk0_1148 = buffer.data(skk0 + 1148);
    const auto *skk0_1149 = buffer.data(skk0 + 1149);
    const auto *skk0_1151 = buffer.data(skk0 + 1151);
    const auto *skk0_1152 = buffer.data(skk0 + 1152);
    const auto *skk0_1155 = buffer.data(skk0 + 1155);
    const auto *skk0_1157 = buffer.data(skk0 + 1157);
    const auto *skk0_1158 = buffer.data(skk0 + 1158);
    const auto *skk0_1161 = buffer.data(skk0 + 1161);
    const auto *skk0_1162 = buffer.data(skk0 + 1162);
    const auto *skk0_1164 = buffer.data(skk0 + 1164);
    const auto *skk0_1166 = buffer.data(skk0 + 1166);
    const auto *skk0_1167 = buffer.data(skk0 + 1167);
    const auto *skk0_1169 = buffer.data(skk0 + 1169);
    const auto *skk0_1170 = buffer.data(skk0 + 1170);
    const auto *skk0_1172 = buffer.data(skk0 + 1172);
    const auto *skk0_1180 = buffer.data(skk0 + 1180);
    const auto *skk0_1182 = buffer.data(skk0 + 1182);
    const auto *skk0_1183 = buffer.data(skk0 + 1183);
    const auto *skk0_1184 = buffer.data(skk0 + 1184);
    const auto *skk0_1185 = buffer.data(skk0 + 1185);
    const auto *skk0_1187 = buffer.data(skk0 + 1187);
    const auto *skk0_1188 = buffer.data(skk0 + 1188);
    const auto *skk0_1191 = buffer.data(skk0 + 1191);
    const auto *skk0_1193 = buffer.data(skk0 + 1193);
    const auto *skk0_1194 = buffer.data(skk0 + 1194);
    const auto *skk0_1197 = buffer.data(skk0 + 1197);
    const auto *skk0_1198 = buffer.data(skk0 + 1198);
    const auto *skk0_1200 = buffer.data(skk0 + 1200);
    const auto *skk0_1202 = buffer.data(skk0 + 1202);
    const auto *skk0_1203 = buffer.data(skk0 + 1203);
    const auto *skk0_1205 = buffer.data(skk0 + 1205);
    const auto *skk0_1206 = buffer.data(skk0 + 1206);
    const auto *skk0_1208 = buffer.data(skk0 + 1208);
    const auto *skk0_1216 = buffer.data(skk0 + 1216);
    const auto *skk0_1218 = buffer.data(skk0 + 1218);
    const auto *skk0_1219 = buffer.data(skk0 + 1219);
    const auto *skk0_1220 = buffer.data(skk0 + 1220);
    const auto *skk0_1221 = buffer.data(skk0 + 1221);
    const auto *skk0_1223 = buffer.data(skk0 + 1223);
    const auto *skk0_1227 = buffer.data(skk0 + 1227);
    const auto *skk0_1230 = buffer.data(skk0 + 1230);
    const auto *skk0_1234 = buffer.data(skk0 + 1234);
    const auto *skk0_1236 = buffer.data(skk0 + 1236);
    const auto *skk0_1239 = buffer.data(skk0 + 1239);
    const auto *skk0_1241 = buffer.data(skk0 + 1241);
    const auto *skk0_1242 = buffer.data(skk0 + 1242);

    const auto *ski_650 = buffer.data(ski + 650);
    const auto *ski_654 = buffer.data(ski + 654);
    const auto *ski_665 = buffer.data(ski + 665);
    const auto *ski_672 = buffer.data(ski + 672);
    const auto *ski_675 = buffer.data(ski + 675);
    const auto *ski_678 = buffer.data(ski + 678);
    const auto *ski_681 = buffer.data(ski + 681);
    const auto *ski_682 = buffer.data(ski + 682);
    const auto *ski_686 = buffer.data(ski + 686);
    const auto *ski_693 = buffer.data(ski + 693);
    const auto *ski_699 = buffer.data(ski + 699);
    const auto *ski_700 = buffer.data(ski + 700);
    const auto *ski_702 = buffer.data(ski + 702);
    const auto *ski_703 = buffer.data(ski + 703);
    const auto *ski_705 = buffer.data(ski + 705);
    const auto *ski_706 = buffer.data(ski + 706);
    const auto *ski_709 = buffer.data(ski + 709);
    const auto *ski_710 = buffer.data(ski + 710);
    const auto *ski_714 = buffer.data(ski + 714);
    const auto *ski_721 = buffer.data(ski + 721);
    const auto *ski_727 = buffer.data(ski + 727);
    const auto *ski_728 = buffer.data(ski + 728);
    const auto *ski_730 = buffer.data(ski + 730);
    const auto *ski_731 = buffer.data(ski + 731);
    const auto *ski_733 = buffer.data(ski + 733);
    const auto *ski_734 = buffer.data(ski + 734);
    const auto *ski_737 = buffer.data(ski + 737);
    const auto *ski_738 = buffer.data(ski + 738);
    const auto *ski_742 = buffer.data(ski + 742);
    const auto *ski_755 = buffer.data(ski + 755);
    const auto *ski_756 = buffer.data(ski + 756);
    const auto *ski_758 = buffer.data(ski + 758);
    const auto *ski_761 = buffer.data(ski + 761);
    const auto *ski_765 = buffer.data(ski + 765);
    const auto *ski_770 = buffer.data(ski + 770);
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
    const auto *ski_955 = buffer.data(ski + 955);
    const auto *ski_958 = buffer.data(ski + 958);
    const auto *ski_962 = buffer.data(ski + 962);
    const auto *ski_964 = buffer.data(ski + 964);
    const auto *ski_967 = buffer.data(ski + 967);
    const auto *ski_969 = buffer.data(ski + 969);
    const auto *ski_970 = buffer.data(ski + 970);

    const auto *skk1_972 = buffer.data(skk1 + 972);
    const auto *skk1_977 = buffer.data(skk1 + 977);
    const auto *skk1_981 = buffer.data(skk1 + 981);
    const auto *skk1_986 = buffer.data(skk1 + 986);
    const auto *skk1_992 = buffer.data(skk1 + 992);
    const auto *skk1_1128 = buffer.data(skk1 + 1128);
    const auto *skk1_1130 = buffer.data(skk1 + 1130);
    const auto *skk1_1131 = buffer.data(skk1 + 1131);
    const auto *skk1_1133 = buffer.data(skk1 + 1133);
    const auto *skk1_1134 = buffer.data(skk1 + 1134);
    const auto *skk1_1136 = buffer.data(skk1 + 1136);
    const auto *skk1_1144 = buffer.data(skk1 + 1144);
    const auto *skk1_1146 = buffer.data(skk1 + 1146);
    const auto *skk1_1147 = buffer.data(skk1 + 1147);
    const auto *skk1_1148 = buffer.data(skk1 + 1148);
    const auto *skk1_1149 = buffer.data(skk1 + 1149);
    const auto *skk1_1151 = buffer.data(skk1 + 1151);
    const auto *skk1_1152 = buffer.data(skk1 + 1152);
    const auto *skk1_1155 = buffer.data(skk1 + 1155);
    const auto *skk1_1157 = buffer.data(skk1 + 1157);
    const auto *skk1_1158 = buffer.data(skk1 + 1158);
    const auto *skk1_1161 = buffer.data(skk1 + 1161);
    const auto *skk1_1162 = buffer.data(skk1 + 1162);
    const auto *skk1_1164 = buffer.data(skk1 + 1164);
    const auto *skk1_1166 = buffer.data(skk1 + 1166);
    const auto *skk1_1167 = buffer.data(skk1 + 1167);
    const auto *skk1_1169 = buffer.data(skk1 + 1169);
    const auto *skk1_1170 = buffer.data(skk1 + 1170);
    const auto *skk1_1172 = buffer.data(skk1 + 1172);
    const auto *skk1_1180 = buffer.data(skk1 + 1180);
    const auto *skk1_1182 = buffer.data(skk1 + 1182);
    const auto *skk1_1183 = buffer.data(skk1 + 1183);
    const auto *skk1_1184 = buffer.data(skk1 + 1184);
    const auto *skk1_1185 = buffer.data(skk1 + 1185);
    const auto *skk1_1187 = buffer.data(skk1 + 1187);
    const auto *skk1_1188 = buffer.data(skk1 + 1188);
    const auto *skk1_1191 = buffer.data(skk1 + 1191);
    const auto *skk1_1193 = buffer.data(skk1 + 1193);
    const auto *skk1_1194 = buffer.data(skk1 + 1194);
    const auto *skk1_1197 = buffer.data(skk1 + 1197);
    const auto *skk1_1198 = buffer.data(skk1 + 1198);
    const auto *skk1_1200 = buffer.data(skk1 + 1200);
    const auto *skk1_1202 = buffer.data(skk1 + 1202);
    const auto *skk1_1203 = buffer.data(skk1 + 1203);
    const auto *skk1_1205 = buffer.data(skk1 + 1205);
    const auto *skk1_1206 = buffer.data(skk1 + 1206);
    const auto *skk1_1208 = buffer.data(skk1 + 1208);
    const auto *skk1_1216 = buffer.data(skk1 + 1216);
    const auto *skk1_1218 = buffer.data(skk1 + 1218);
    const auto *skk1_1219 = buffer.data(skk1 + 1219);
    const auto *skk1_1220 = buffer.data(skk1 + 1220);
    const auto *skk1_1221 = buffer.data(skk1 + 1221);
    const auto *skk1_1223 = buffer.data(skk1 + 1223);
    const auto *skk1_1227 = buffer.data(skk1 + 1227);
    const auto *skk1_1230 = buffer.data(skk1 + 1230);
    const auto *skk1_1234 = buffer.data(skk1 + 1234);
    const auto *skk1_1236 = buffer.data(skk1 + 1236);
    const auto *skk1_1239 = buffer.data(skk1 + 1239);
    const auto *skk1_1241 = buffer.data(skk1 + 1241);
    const auto *skk1_1242 = buffer.data(skk1 + 1242);

    const auto *sli_874 = buffer.data(sli + 874);
    const auto *sli_877 = buffer.data(sli + 877);
    const auto *sli_878 = buffer.data(sli + 878);
    const auto *sli_882 = buffer.data(sli + 882);
    const auto *sli_889 = buffer.data(sli + 889);
    const auto *sli_890 = buffer.data(sli + 890);
    const auto *sli_891 = buffer.data(sli + 891);
    const auto *sli_892 = buffer.data(sli + 892);
    const auto *sli_893 = buffer.data(sli + 893);
    const auto *sli_894 = buffer.data(sli + 894);
    const auto *sli_895 = buffer.data(sli + 895);
    const auto *sli_896 = buffer.data(sli + 896);
    const auto *sli_898 = buffer.data(sli + 898);
    const auto *sli_899 = buffer.data(sli + 899);
    const auto *sli_901 = buffer.data(sli + 901);
    const auto *sli_902 = buffer.data(sli + 902);
    const auto *sli_905 = buffer.data(sli + 905);
    const auto *sli_906 = buffer.data(sli + 906);
    const auto *sli_910 = buffer.data(sli + 910);
    const auto *sli_917 = buffer.data(sli + 917);
    const auto *sli_918 = buffer.data(sli + 918);
    const auto *sli_919 = buffer.data(sli + 919);
    const auto *sli_920 = buffer.data(sli + 920);
    const auto *sli_921 = buffer.data(sli + 921);
    const auto *sli_922 = buffer.data(sli + 922);
    const auto *sli_923 = buffer.data(sli + 923);
    const auto *sli_924 = buffer.data(sli + 924);
    const auto *sli_926 = buffer.data(sli + 926);
    const auto *sli_927 = buffer.data(sli + 927);
    const auto *sli_929 = buffer.data(sli + 929);
    const auto *sli_930 = buffer.data(sli + 930);
    const auto *sli_933 = buffer.data(sli + 933);
    const auto *sli_934 = buffer.data(sli + 934);
    const auto *sli_938 = buffer.data(sli + 938);
    const auto *sli_945 = buffer.data(sli + 945);
    const auto *sli_946 = buffer.data(sli + 946);
    const auto *sli_947 = buffer.data(sli + 947);
    const auto *sli_948 = buffer.data(sli + 948);
    const auto *sli_949 = buffer.data(sli + 949);
    const auto *sli_950 = buffer.data(sli + 950);
    const auto *sli_951 = buffer.data(sli + 951);
    const auto *sli_952 = buffer.data(sli + 952);
    const auto *sli_954 = buffer.data(sli + 954);
    const auto *sli_955 = buffer.data(sli + 955);
    const auto *sli_957 = buffer.data(sli + 957);
    const auto *sli_958 = buffer.data(sli + 958);
    const auto *sli_961 = buffer.data(sli + 961);
    const auto *sli_962 = buffer.data(sli + 962);
    const auto *sli_966 = buffer.data(sli + 966);

#pragma omp simd aligned(t_1127, t_1128, t_1129, pb_x, pc_x, pc_y, pc_z, skk0_1128, ski_650, \
                         ski_681, ski_880, skk1_1128, sli_874, \
                         sli_877 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1127[k] = f_15 * ski_650[k]
                    + f_3 * pc_z[k] * sli_874[k];

        t_1128[k] = pb_x[k] * skk0_1128[k]
                    + f_15 * ski_880[k]
                    - f_12 * pc_x[k] * skk1_1128[k];

        t_1129[k] = f_16 * ski_681[k]
                    + f_3 * pc_y[k] * sli_877[k];
    }

#pragma omp simd aligned(t_1130, t_1131, t_1132, pb_x, pc_x, pc_z, skk0_1130, skk0_1131, \
                         ski_654, ski_882, ski_883, skk1_1130, skk1_1131, \
                         sli_878 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1130[k] = pb_x[k] * skk0_1130[k]
                    + f_15 * ski_882[k]
                    - f_12 * pc_x[k] * skk1_1130[k];

        t_1131[k] = pb_x[k] * skk0_1131[k]
                    + f_14 * ski_883[k]
                    - f_12 * pc_x[k] * skk1_1131[k];

        t_1132[k] = f_15 * ski_654[k]
                    + f_3 * pc_z[k] * sli_878[k];
    }

#pragma omp simd aligned(t_1133, t_1134, t_1135, pb_x, pc_x, pc_y, skk0_1133, skk0_1134, \
                         ski_686, ski_885, ski_886, skk1_1133, skk1_1134, \
                         sli_882 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1133[k] = pb_x[k] * skk0_1133[k]
                    + f_14 * ski_885[k]
                    - f_12 * pc_x[k] * skk1_1133[k];

        t_1134[k] = pb_x[k] * skk0_1134[k]
                    + f_14 * ski_886[k]
                    - f_12 * pc_x[k] * skk1_1134[k];

        t_1135[k] = f_16 * ski_686[k]
                    + f_3 * pc_y[k] * sli_882[k];
    }

#pragma omp simd aligned(t_1136, t_1137, t_1138, t_1139, pb_x, pc_x, skk0_1136, ski_888, \
                         ski_889, ski_890, ski_891, skk1_1136, sli_889, sli_890, \
                         sli_891 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1136[k] = pb_x[k] * skk0_1136[k]
                    + f_14 * ski_888[k]
                    - f_12 * pc_x[k] * skk1_1136[k];

        t_1137[k] = f_13 * ski_889[k]
                    + f_3 * pc_x[k] * sli_889[k];

        t_1138[k] = f_13 * ski_890[k]
                    + f_3 * pc_x[k] * sli_890[k];

        t_1139[k] = f_13 * ski_891[k]
                    + f_3 * pc_x[k] * sli_891[k];
    }

#pragma omp simd aligned(t_1140, t_1141, t_1142, t_1143, pc_x, ski_892, ski_893, ski_894, \
                         ski_895, sli_892, sli_893, sli_894, sli_895 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1140[k] = f_13 * ski_892[k]
                    + f_3 * pc_x[k] * sli_892[k];

        t_1141[k] = f_13 * ski_893[k]
                    + f_3 * pc_x[k] * sli_893[k];

        t_1142[k] = f_13 * ski_894[k]
                    + f_3 * pc_x[k] * sli_894[k];

        t_1143[k] = f_13 * ski_895[k]
                    + f_3 * pc_x[k] * sli_895[k];
    }

#pragma omp simd aligned(t_1144, t_1145, t_1146, t_1147, pb_x, pc_x, pc_z, skk0_1144, \
                         skk0_1146, skk0_1147, ski_665, skk1_1144, skk1_1146, skk1_1147, \
                         sli_889 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1144[k] = pb_x[k] * skk0_1144[k]
                    - f_12 * pc_x[k] * skk1_1144[k];

        t_1145[k] = f_15 * ski_665[k]
                    + f_3 * pc_z[k] * sli_889[k];

        t_1146[k] = pb_x[k] * skk0_1146[k]
                    - f_12 * pc_x[k] * skk1_1146[k];

        t_1147[k] = pb_x[k] * skk0_1147[k]
                    - f_12 * pc_x[k] * skk1_1147[k];
    }

#pragma omp simd aligned(t_1148, t_1149, t_1150, t_1151, pb_x, pc_x, pc_y, skk0_1148, \
                         skk0_1149, skk0_1151, ski_699, skk1_1148, skk1_1149, skk1_1151, \
                         sli_895 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1148[k] = pb_x[k] * skk0_1148[k]
                    - f_12 * pc_x[k] * skk1_1148[k];

        t_1149[k] = pb_x[k] * skk0_1149[k]
                    - f_12 * pc_x[k] * skk1_1149[k];

        t_1150[k] = f_16 * ski_699[k]
                    + f_3 * pc_y[k] * sli_895[k];

        t_1151[k] = pb_x[k] * skk0_1151[k]
                    - f_12 * pc_x[k] * skk1_1151[k];
    }

#pragma omp simd aligned(t_1152, t_1153, t_1154, pb_x, pc_x, pc_y, pc_z, skk0_1152, ski_672, \
                         ski_700, ski_896, skk1_1152, sli_896 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1152[k] = pb_x[k] * skk0_1152[k]
                    + f_18 * ski_896[k]
                    - f_12 * pc_x[k] * skk1_1152[k];

        t_1153[k] = f_15 * ski_700[k]
                    + f_3 * pc_y[k] * sli_896[k];

        t_1154[k] = f_16 * ski_672[k]
                    + f_3 * pc_z[k] * sli_896[k];
    }

#pragma omp simd aligned(t_1155, t_1156, t_1157, pb_x, pc_x, pc_y, skk0_1155, skk0_1157, \
                         ski_702, ski_899, ski_901, skk1_1155, skk1_1157, \
                         sli_898 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1155[k] = pb_x[k] * skk0_1155[k]
                    + f_17 * ski_899[k]
                    - f_12 * pc_x[k] * skk1_1155[k];

        t_1156[k] = f_15 * ski_702[k]
                    + f_3 * pc_y[k] * sli_898[k];

        t_1157[k] = pb_x[k] * skk0_1157[k]
                    + f_17 * ski_901[k]
                    - f_12 * pc_x[k] * skk1_1157[k];
    }

#pragma omp simd aligned(t_1158, t_1159, t_1160, pb_x, pc_x, pc_y, pc_z, skk0_1158, ski_675, \
                         ski_705, ski_902, skk1_1158, sli_899, \
                         sli_901 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1158[k] = pb_x[k] * skk0_1158[k]
                    + f_16 * ski_902[k]
                    - f_12 * pc_x[k] * skk1_1158[k];

        t_1159[k] = f_16 * ski_675[k]
                    + f_3 * pc_z[k] * sli_899[k];

        t_1160[k] = f_15 * ski_705[k]
                    + f_3 * pc_y[k] * sli_901[k];
    }

#pragma omp simd aligned(t_1161, t_1162, t_1163, pb_x, pc_x, pc_z, skk0_1161, skk0_1162, \
                         ski_678, ski_905, ski_906, skk1_1161, skk1_1162, \
                         sli_902 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1161[k] = pb_x[k] * skk0_1161[k]
                    + f_16 * ski_905[k]
                    - f_12 * pc_x[k] * skk1_1161[k];

        t_1162[k] = pb_x[k] * skk0_1162[k]
                    + f_15 * ski_906[k]
                    - f_12 * pc_x[k] * skk1_1162[k];

        t_1163[k] = f_16 * ski_678[k]
                    + f_3 * pc_z[k] * sli_902[k];
    }

#pragma omp simd aligned(t_1164, t_1165, t_1166, pb_x, pc_x, pc_y, skk0_1164, skk0_1166, \
                         ski_709, ski_908, ski_910, skk1_1164, skk1_1166, \
                         sli_905 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1164[k] = pb_x[k] * skk0_1164[k]
                    + f_15 * ski_908[k]
                    - f_12 * pc_x[k] * skk1_1164[k];

        t_1165[k] = f_15 * ski_709[k]
                    + f_3 * pc_y[k] * sli_905[k];

        t_1166[k] = pb_x[k] * skk0_1166[k]
                    + f_15 * ski_910[k]
                    - f_12 * pc_x[k] * skk1_1166[k];
    }

#pragma omp simd aligned(t_1167, t_1168, t_1169, pb_x, pc_x, pc_z, skk0_1167, skk0_1169, \
                         ski_682, ski_911, ski_913, skk1_1167, skk1_1169, \
                         sli_906 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1167[k] = pb_x[k] * skk0_1167[k]
                    + f_14 * ski_911[k]
                    - f_12 * pc_x[k] * skk1_1167[k];

        t_1168[k] = f_16 * ski_682[k]
                    + f_3 * pc_z[k] * sli_906[k];

        t_1169[k] = pb_x[k] * skk0_1169[k]
                    + f_14 * ski_913[k]
                    - f_12 * pc_x[k] * skk1_1169[k];
    }

#pragma omp simd aligned(t_1170, t_1171, t_1172, pb_x, pc_x, pc_y, skk0_1170, skk0_1172, \
                         ski_714, ski_914, ski_916, skk1_1170, skk1_1172, \
                         sli_910 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1170[k] = pb_x[k] * skk0_1170[k]
                    + f_14 * ski_914[k]
                    - f_12 * pc_x[k] * skk1_1170[k];

        t_1171[k] = f_15 * ski_714[k]
                    + f_3 * pc_y[k] * sli_910[k];

        t_1172[k] = pb_x[k] * skk0_1172[k]
                    + f_14 * ski_916[k]
                    - f_12 * pc_x[k] * skk1_1172[k];
    }

#pragma omp simd aligned(t_1173, t_1174, t_1175, t_1176, t_1177, pc_x, ski_917, ski_918, \
                         ski_919, ski_920, ski_921, sli_917, sli_918, sli_919, sli_920, \
                         sli_921 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1173[k] = f_13 * ski_917[k]
                    + f_3 * pc_x[k] * sli_917[k];

        t_1174[k] = f_13 * ski_918[k]
                    + f_3 * pc_x[k] * sli_918[k];

        t_1175[k] = f_13 * ski_919[k]
                    + f_3 * pc_x[k] * sli_919[k];

        t_1176[k] = f_13 * ski_920[k]
                    + f_3 * pc_x[k] * sli_920[k];

        t_1177[k] = f_13 * ski_921[k]
                    + f_3 * pc_x[k] * sli_921[k];
    }

#pragma omp simd aligned(t_1178, t_1179, t_1180, t_1181, pb_x, pc_x, pc_z, skk0_1180, ski_693, \
                         ski_922, ski_923, skk1_1180, sli_917, sli_922, \
                         sli_923 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1178[k] = f_13 * ski_922[k]
                    + f_3 * pc_x[k] * sli_922[k];

        t_1179[k] = f_13 * ski_923[k]
                    + f_3 * pc_x[k] * sli_923[k];

        t_1180[k] = pb_x[k] * skk0_1180[k]
                    - f_12 * pc_x[k] * skk1_1180[k];

        t_1181[k] = f_16 * ski_693[k]
                    + f_3 * pc_z[k] * sli_917[k];
    }

#pragma omp simd aligned(t_1182, t_1183, t_1184, t_1185, pb_x, pc_x, skk0_1182, skk0_1183, \
                         skk0_1184, skk0_1185, skk1_1182, skk1_1183, skk1_1184, \
                         skk1_1185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1182[k] = pb_x[k] * skk0_1182[k]
                    - f_12 * pc_x[k] * skk1_1182[k];

        t_1183[k] = pb_x[k] * skk0_1183[k]
                    - f_12 * pc_x[k] * skk1_1183[k];

        t_1184[k] = pb_x[k] * skk0_1184[k]
                    - f_12 * pc_x[k] * skk1_1184[k];

        t_1185[k] = pb_x[k] * skk0_1185[k]
                    - f_12 * pc_x[k] * skk1_1185[k];
    }

#pragma omp simd aligned(t_1186, t_1187, t_1188, t_1189, pb_x, pc_x, pc_y, skk0_1187, \
                         skk0_1188, ski_727, ski_728, ski_924, skk1_1187, skk1_1188, sli_923, \
                         sli_924 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1186[k] = f_15 * ski_727[k]
                    + f_3 * pc_y[k] * sli_923[k];

        t_1187[k] = pb_x[k] * skk0_1187[k]
                    - f_12 * pc_x[k] * skk1_1187[k];

        t_1188[k] = pb_x[k] * skk0_1188[k]
                    + f_18 * ski_924[k]
                    - f_12 * pc_x[k] * skk1_1188[k];

        t_1189[k] = f_14 * ski_728[k]
                    + f_3 * pc_y[k] * sli_924[k];
    }

#pragma omp simd aligned(t_1190, t_1191, t_1192, pb_x, pc_x, pc_y, pc_z, skk0_1191, ski_700, \
                         ski_730, ski_927, skk1_1191, sli_924, \
                         sli_926 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1190[k] = f_17 * ski_700[k]
                    + f_3 * pc_z[k] * sli_924[k];

        t_1191[k] = pb_x[k] * skk0_1191[k]
                    + f_17 * ski_927[k]
                    - f_12 * pc_x[k] * skk1_1191[k];

        t_1192[k] = f_14 * ski_730[k]
                    + f_3 * pc_y[k] * sli_926[k];
    }

#pragma omp simd aligned(t_1193, t_1194, t_1195, pb_x, pc_x, pc_z, skk0_1193, skk0_1194, \
                         ski_703, ski_929, ski_930, skk1_1193, skk1_1194, \
                         sli_927 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1193[k] = pb_x[k] * skk0_1193[k]
                    + f_17 * ski_929[k]
                    - f_12 * pc_x[k] * skk1_1193[k];

        t_1194[k] = pb_x[k] * skk0_1194[k]
                    + f_16 * ski_930[k]
                    - f_12 * pc_x[k] * skk1_1194[k];

        t_1195[k] = f_17 * ski_703[k]
                    + f_3 * pc_z[k] * sli_927[k];
    }

#pragma omp simd aligned(t_1196, t_1197, t_1198, pb_x, pc_x, pc_y, skk0_1197, skk0_1198, \
                         ski_733, ski_933, ski_934, skk1_1197, skk1_1198, \
                         sli_929 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1196[k] = f_14 * ski_733[k]
                    + f_3 * pc_y[k] * sli_929[k];

        t_1197[k] = pb_x[k] * skk0_1197[k]
                    + f_16 * ski_933[k]
                    - f_12 * pc_x[k] * skk1_1197[k];

        t_1198[k] = pb_x[k] * skk0_1198[k]
                    + f_15 * ski_934[k]
                    - f_12 * pc_x[k] * skk1_1198[k];
    }

#pragma omp simd aligned(t_1199, t_1200, t_1201, pb_x, pc_x, pc_y, pc_z, skk0_1200, ski_706, \
                         ski_737, ski_936, skk1_1200, sli_930, \
                         sli_933 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1199[k] = f_17 * ski_706[k]
                    + f_3 * pc_z[k] * sli_930[k];

        t_1200[k] = pb_x[k] * skk0_1200[k]
                    + f_15 * ski_936[k]
                    - f_12 * pc_x[k] * skk1_1200[k];

        t_1201[k] = f_14 * ski_737[k]
                    + f_3 * pc_y[k] * sli_933[k];
    }

#pragma omp simd aligned(t_1202, t_1203, t_1204, pb_x, pc_x, pc_z, skk0_1202, skk0_1203, \
                         ski_710, ski_938, ski_939, skk1_1202, skk1_1203, \
                         sli_934 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1202[k] = pb_x[k] * skk0_1202[k]
                    + f_15 * ski_938[k]
                    - f_12 * pc_x[k] * skk1_1202[k];

        t_1203[k] = pb_x[k] * skk0_1203[k]
                    + f_14 * ski_939[k]
                    - f_12 * pc_x[k] * skk1_1203[k];

        t_1204[k] = f_17 * ski_710[k]
                    + f_3 * pc_z[k] * sli_934[k];
    }

#pragma omp simd aligned(t_1205, t_1206, t_1207, pb_x, pc_x, pc_y, skk0_1205, skk0_1206, \
                         ski_742, ski_941, ski_942, skk1_1205, skk1_1206, \
                         sli_938 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1205[k] = pb_x[k] * skk0_1205[k]
                    + f_14 * ski_941[k]
                    - f_12 * pc_x[k] * skk1_1205[k];

        t_1206[k] = pb_x[k] * skk0_1206[k]
                    + f_14 * ski_942[k]
                    - f_12 * pc_x[k] * skk1_1206[k];

        t_1207[k] = f_14 * ski_742[k]
                    + f_3 * pc_y[k] * sli_938[k];
    }

#pragma omp simd aligned(t_1208, t_1209, t_1210, t_1211, pb_x, pc_x, skk0_1208, ski_944, \
                         ski_945, ski_946, ski_947, skk1_1208, sli_945, sli_946, \
                         sli_947 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1208[k] = pb_x[k] * skk0_1208[k]
                    + f_14 * ski_944[k]
                    - f_12 * pc_x[k] * skk1_1208[k];

        t_1209[k] = f_13 * ski_945[k]
                    + f_3 * pc_x[k] * sli_945[k];

        t_1210[k] = f_13 * ski_946[k]
                    + f_3 * pc_x[k] * sli_946[k];

        t_1211[k] = f_13 * ski_947[k]
                    + f_3 * pc_x[k] * sli_947[k];
    }

#pragma omp simd aligned(t_1212, t_1213, t_1214, t_1215, pc_x, ski_948, ski_949, ski_950, \
                         ski_951, sli_948, sli_949, sli_950, sli_951 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1212[k] = f_13 * ski_948[k]
                    + f_3 * pc_x[k] * sli_948[k];

        t_1213[k] = f_13 * ski_949[k]
                    + f_3 * pc_x[k] * sli_949[k];

        t_1214[k] = f_13 * ski_950[k]
                    + f_3 * pc_x[k] * sli_950[k];

        t_1215[k] = f_13 * ski_951[k]
                    + f_3 * pc_x[k] * sli_951[k];
    }

#pragma omp simd aligned(t_1216, t_1217, t_1218, t_1219, pb_x, pc_x, pc_z, skk0_1216, \
                         skk0_1218, skk0_1219, ski_721, skk1_1216, skk1_1218, skk1_1219, \
                         sli_945 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1216[k] = pb_x[k] * skk0_1216[k]
                    - f_12 * pc_x[k] * skk1_1216[k];

        t_1217[k] = f_17 * ski_721[k]
                    + f_3 * pc_z[k] * sli_945[k];

        t_1218[k] = pb_x[k] * skk0_1218[k]
                    - f_12 * pc_x[k] * skk1_1218[k];

        t_1219[k] = pb_x[k] * skk0_1219[k]
                    - f_12 * pc_x[k] * skk1_1219[k];
    }

#pragma omp simd aligned(t_1220, t_1221, t_1222, t_1223, pb_x, pc_x, pc_y, skk0_1220, \
                         skk0_1221, skk0_1223, ski_755, skk1_1220, skk1_1221, skk1_1223, \
                         sli_951 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1220[k] = pb_x[k] * skk0_1220[k]
                    - f_12 * pc_x[k] * skk1_1220[k];

        t_1221[k] = pb_x[k] * skk0_1221[k]
                    - f_12 * pc_x[k] * skk1_1221[k];

        t_1222[k] = f_14 * ski_755[k]
                    + f_3 * pc_y[k] * sli_951[k];

        t_1223[k] = pb_x[k] * skk0_1223[k]
                    - f_12 * pc_x[k] * skk1_1223[k];
    }

#pragma omp simd aligned(t_1224, t_1225, t_1226, pb_y, pc_y, pc_z, skk0_972, ski_728, ski_756, \
                         skk1_972, sli_952 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1224[k] = pb_y[k] * skk0_972[k]
                    - f_12 * pc_y[k] * skk1_972[k];

        t_1225[k] = f_13 * ski_756[k]
                    + f_3 * pc_y[k] * sli_952[k];

        t_1226[k] = f_19 * ski_728[k]
                    + f_3 * pc_z[k] * sli_952[k];
    }

#pragma omp simd aligned(t_1227, t_1228, t_1229, pb_x, pb_y, pc_x, pc_y, skk0_977, skk0_1227, \
                         ski_758, ski_955, skk1_977, skk1_1227, \
                         sli_954 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1227[k] = pb_x[k] * skk0_1227[k]
                    + f_17 * ski_955[k]
                    - f_12 * pc_x[k] * skk1_1227[k];

        t_1228[k] = f_13 * ski_758[k]
                    + f_3 * pc_y[k] * sli_954[k];

        t_1229[k] = pb_y[k] * skk0_977[k]
                    - f_12 * pc_y[k] * skk1_977[k];
    }

#pragma omp simd aligned(t_1230, t_1231, t_1232, pb_x, pc_x, pc_y, pc_z, skk0_1230, ski_731, \
                         ski_761, ski_958, skk1_1230, sli_955, \
                         sli_957 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1230[k] = pb_x[k] * skk0_1230[k]
                    + f_16 * ski_958[k]
                    - f_12 * pc_x[k] * skk1_1230[k];

        t_1231[k] = f_19 * ski_731[k]
                    + f_3 * pc_z[k] * sli_955[k];

        t_1232[k] = f_13 * ski_761[k]
                    + f_3 * pc_y[k] * sli_957[k];
    }

#pragma omp simd aligned(t_1233, t_1234, t_1235, pb_x, pb_y, pc_x, pc_y, pc_z, skk0_981, \
                         skk0_1234, ski_734, ski_962, skk1_981, skk1_1234, \
                         sli_958 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1233[k] = pb_y[k] * skk0_981[k]
                    - f_12 * pc_y[k] * skk1_981[k];

        t_1234[k] = pb_x[k] * skk0_1234[k]
                    + f_15 * ski_962[k]
                    - f_12 * pc_x[k] * skk1_1234[k];

        t_1235[k] = f_19 * ski_734[k]
                    + f_3 * pc_z[k] * sli_958[k];
    }

#pragma omp simd aligned(t_1236, t_1237, t_1238, pb_x, pb_y, pc_x, pc_y, skk0_986, skk0_1236, \
                         ski_765, ski_964, skk1_986, skk1_1236, \
                         sli_961 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1236[k] = pb_x[k] * skk0_1236[k]
                    + f_15 * ski_964[k]
                    - f_12 * pc_x[k] * skk1_1236[k];

        t_1237[k] = f_13 * ski_765[k]
                    + f_3 * pc_y[k] * sli_961[k];

        t_1238[k] = pb_y[k] * skk0_986[k]
                    - f_12 * pc_y[k] * skk1_986[k];
    }

#pragma omp simd aligned(t_1239, t_1240, t_1241, pb_x, pc_x, pc_z, skk0_1239, skk0_1241, \
                         ski_738, ski_967, ski_969, skk1_1239, skk1_1241, \
                         sli_962 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1239[k] = pb_x[k] * skk0_1239[k]
                    + f_14 * ski_967[k]
                    - f_12 * pc_x[k] * skk1_1239[k];

        t_1240[k] = f_19 * ski_738[k]
                    + f_3 * pc_z[k] * sli_962[k];

        t_1241[k] = pb_x[k] * skk0_1241[k]
                    + f_14 * ski_969[k]
                    - f_12 * pc_x[k] * skk1_1241[k];
    }

#pragma omp simd aligned(t_1242, t_1243, t_1244, pb_x, pb_y, pc_x, pc_y, skk0_992, skk0_1242, \
                         ski_770, ski_970, skk1_992, skk1_1242, \
                         sli_966 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1242[k] = pb_x[k] * skk0_1242[k]
                    + f_14 * ski_970[k]
                    - f_12 * pc_x[k] * skk1_1242[k];

        t_1243[k] = f_13 * ski_770[k]
                    + f_3 * pc_y[k] * sli_966[k];

        t_1244[k] = pb_y[k] * skk0_992[k]
                    - f_12 * pc_y[k] * skk1_992[k];
    }
}

static auto
compute_prim_slk_three_center_electron_repulsion_0_piece11(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t skk0,
                                                           const size_t ski, const size_t skk1,
                                                           const size_t slh0, const size_t slh1,
                                                           const size_t sli, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
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
    const auto f_19 = 3.0 / q;

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
    auto *t_1296 = buffer.data(target + 1296);
    auto *t_1297 = buffer.data(target + 1297);
    auto *t_1298 = buffer.data(target + 1298);
    auto *t_1299 = buffer.data(target + 1299);
    auto *t_1300 = buffer.data(target + 1300);
    auto *t_1301 = buffer.data(target + 1301);
    auto *t_1302 = buffer.data(target + 1302);
    auto *t_1303 = buffer.data(target + 1303);
    auto *t_1304 = buffer.data(target + 1304);
    auto *t_1305 = buffer.data(target + 1305);
    auto *t_1306 = buffer.data(target + 1306);
    auto *t_1307 = buffer.data(target + 1307);
    auto *t_1308 = buffer.data(target + 1308);
    auto *t_1309 = buffer.data(target + 1309);
    auto *t_1310 = buffer.data(target + 1310);
    auto *t_1311 = buffer.data(target + 1311);
    auto *t_1312 = buffer.data(target + 1312);
    auto *t_1313 = buffer.data(target + 1313);
    auto *t_1314 = buffer.data(target + 1314);
    auto *t_1315 = buffer.data(target + 1315);
    auto *t_1316 = buffer.data(target + 1316);
    auto *t_1317 = buffer.data(target + 1317);
    auto *t_1318 = buffer.data(target + 1318);
    auto *t_1319 = buffer.data(target + 1319);
    auto *t_1320 = buffer.data(target + 1320);
    auto *t_1321 = buffer.data(target + 1321);
    auto *t_1322 = buffer.data(target + 1322);
    auto *t_1323 = buffer.data(target + 1323);
    auto *t_1324 = buffer.data(target + 1324);
    auto *t_1325 = buffer.data(target + 1325);
    auto *t_1326 = buffer.data(target + 1326);
    auto *t_1327 = buffer.data(target + 1327);
    auto *t_1328 = buffer.data(target + 1328);
    auto *t_1329 = buffer.data(target + 1329);
    auto *t_1330 = buffer.data(target + 1330);
    auto *t_1331 = buffer.data(target + 1331);
    auto *t_1332 = buffer.data(target + 1332);
    auto *t_1333 = buffer.data(target + 1333);
    auto *t_1334 = buffer.data(target + 1334);
    auto *t_1335 = buffer.data(target + 1335);
    auto *t_1336 = buffer.data(target + 1336);
    auto *t_1337 = buffer.data(target + 1337);
    auto *t_1338 = buffer.data(target + 1338);
    auto *t_1339 = buffer.data(target + 1339);
    auto *t_1340 = buffer.data(target + 1340);
    auto *t_1341 = buffer.data(target + 1341);
    auto *t_1342 = buffer.data(target + 1342);
    auto *t_1343 = buffer.data(target + 1343);
    auto *t_1344 = buffer.data(target + 1344);
    auto *t_1345 = buffer.data(target + 1345);
    auto *t_1346 = buffer.data(target + 1346);
    auto *t_1347 = buffer.data(target + 1347);
    auto *t_1348 = buffer.data(target + 1348);
    auto *t_1349 = buffer.data(target + 1349);
    auto *t_1350 = buffer.data(target + 1350);
    auto *t_1351 = buffer.data(target + 1351);
    auto *t_1352 = buffer.data(target + 1352);
    auto *t_1353 = buffer.data(target + 1353);
    auto *t_1354 = buffer.data(target + 1354);
    auto *t_1355 = buffer.data(target + 1355);
    auto *t_1356 = buffer.data(target + 1356);
    auto *t_1357 = buffer.data(target + 1357);
    auto *t_1358 = buffer.data(target + 1358);
    auto *t_1359 = buffer.data(target + 1359);
    auto *t_1360 = buffer.data(target + 1360);
    auto *t_1361 = buffer.data(target + 1361);
    auto *t_1362 = buffer.data(target + 1362);
    auto *t_1363 = buffer.data(target + 1363);
    auto *t_1364 = buffer.data(target + 1364);
    auto *t_1365 = buffer.data(target + 1365);
    auto *t_1366 = buffer.data(target + 1366);
    auto *t_1367 = buffer.data(target + 1367);
    auto *t_1368 = buffer.data(target + 1368);
    auto *t_1369 = buffer.data(target + 1369);
    auto *t_1370 = buffer.data(target + 1370);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *skk0_1008 = buffer.data(skk0 + 1008);
    const auto *skk0_1011 = buffer.data(skk0 + 1011);
    const auto *skk0_1014 = buffer.data(skk0 + 1014);
    const auto *skk0_1018 = buffer.data(skk0 + 1018);
    const auto *skk0_1023 = buffer.data(skk0 + 1023);
    const auto *skk0_1036 = buffer.data(skk0 + 1036);
    const auto *skk0_1038 = buffer.data(skk0 + 1038);
    const auto *skk0_1039 = buffer.data(skk0 + 1039);
    const auto *skk0_1040 = buffer.data(skk0 + 1040);
    const auto *skk0_1041 = buffer.data(skk0 + 1041);
    const auto *skk0_1252 = buffer.data(skk0 + 1252);
    const auto *skk0_1254 = buffer.data(skk0 + 1254);
    const auto *skk0_1255 = buffer.data(skk0 + 1255);
    const auto *skk0_1256 = buffer.data(skk0 + 1256);
    const auto *skk0_1257 = buffer.data(skk0 + 1257);
    const auto *skk0_1259 = buffer.data(skk0 + 1259);
    const auto *skk0_1260 = buffer.data(skk0 + 1260);
    const auto *skk0_1263 = buffer.data(skk0 + 1263);
    const auto *skk0_1265 = buffer.data(skk0 + 1265);
    const auto *skk0_1266 = buffer.data(skk0 + 1266);
    const auto *skk0_1269 = buffer.data(skk0 + 1269);
    const auto *skk0_1270 = buffer.data(skk0 + 1270);
    const auto *skk0_1272 = buffer.data(skk0 + 1272);
    const auto *skk0_1274 = buffer.data(skk0 + 1274);
    const auto *skk0_1275 = buffer.data(skk0 + 1275);
    const auto *skk0_1277 = buffer.data(skk0 + 1277);
    const auto *skk0_1278 = buffer.data(skk0 + 1278);
    const auto *skk0_1280 = buffer.data(skk0 + 1280);
    const auto *skk0_1288 = buffer.data(skk0 + 1288);
    const auto *skk0_1290 = buffer.data(skk0 + 1290);
    const auto *skk0_1291 = buffer.data(skk0 + 1291);
    const auto *skk0_1292 = buffer.data(skk0 + 1292);
    const auto *skk0_1293 = buffer.data(skk0 + 1293);
    const auto *skk0_1295 = buffer.data(skk0 + 1295);

    const auto *ski_749 = buffer.data(ski + 749);
    const auto *ski_756 = buffer.data(ski + 756);
    const auto *ski_759 = buffer.data(ski + 759);
    const auto *ski_762 = buffer.data(ski + 762);
    const auto *ski_766 = buffer.data(ski + 766);
    const auto *ski_777 = buffer.data(ski + 777);
    const auto *ski_783 = buffer.data(ski + 783);
    const auto *ski_784 = buffer.data(ski + 784);
    const auto *ski_786 = buffer.data(ski + 786);
    const auto *ski_787 = buffer.data(ski + 787);
    const auto *ski_789 = buffer.data(ski + 789);
    const auto *ski_790 = buffer.data(ski + 790);
    const auto *ski_793 = buffer.data(ski + 793);
    const auto *ski_794 = buffer.data(ski + 794);
    const auto *ski_798 = buffer.data(ski + 798);
    const auto *ski_805 = buffer.data(ski + 805);
    const auto *ski_806 = buffer.data(ski + 806);
    const auto *ski_807 = buffer.data(ski + 807);
    const auto *ski_808 = buffer.data(ski + 808);
    const auto *ski_809 = buffer.data(ski + 809);
    const auto *ski_810 = buffer.data(ski + 810);
    const auto *ski_811 = buffer.data(ski + 811);
    const auto *ski_812 = buffer.data(ski + 812);
    const auto *ski_814 = buffer.data(ski + 814);
    const auto *ski_817 = buffer.data(ski + 817);
    const auto *ski_821 = buffer.data(ski + 821);
    const auto *ski_826 = buffer.data(ski + 826);
    const auto *ski_839 = buffer.data(ski + 839);
    const auto *ski_840 = buffer.data(ski + 840);
    const auto *ski_973 = buffer.data(ski + 973);
    const auto *ski_974 = buffer.data(ski + 974);
    const auto *ski_975 = buffer.data(ski + 975);
    const auto *ski_976 = buffer.data(ski + 976);
    const auto *ski_977 = buffer.data(ski + 977);
    const auto *ski_978 = buffer.data(ski + 978);
    const auto *ski_979 = buffer.data(ski + 979);
    const auto *ski_980 = buffer.data(ski + 980);
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
    const auto *ski_1000 = buffer.data(ski + 1000);
    const auto *ski_1001 = buffer.data(ski + 1001);
    const auto *ski_1002 = buffer.data(ski + 1002);
    const auto *ski_1003 = buffer.data(ski + 1003);
    const auto *ski_1004 = buffer.data(ski + 1004);
    const auto *ski_1005 = buffer.data(ski + 1005);
    const auto *ski_1006 = buffer.data(ski + 1006);
    const auto *ski_1007 = buffer.data(ski + 1007);

    const auto *skk1_1008 = buffer.data(skk1 + 1008);
    const auto *skk1_1011 = buffer.data(skk1 + 1011);
    const auto *skk1_1014 = buffer.data(skk1 + 1014);
    const auto *skk1_1018 = buffer.data(skk1 + 1018);
    const auto *skk1_1023 = buffer.data(skk1 + 1023);
    const auto *skk1_1036 = buffer.data(skk1 + 1036);
    const auto *skk1_1038 = buffer.data(skk1 + 1038);
    const auto *skk1_1039 = buffer.data(skk1 + 1039);
    const auto *skk1_1040 = buffer.data(skk1 + 1040);
    const auto *skk1_1041 = buffer.data(skk1 + 1041);
    const auto *skk1_1252 = buffer.data(skk1 + 1252);
    const auto *skk1_1254 = buffer.data(skk1 + 1254);
    const auto *skk1_1255 = buffer.data(skk1 + 1255);
    const auto *skk1_1256 = buffer.data(skk1 + 1256);
    const auto *skk1_1257 = buffer.data(skk1 + 1257);
    const auto *skk1_1259 = buffer.data(skk1 + 1259);
    const auto *skk1_1260 = buffer.data(skk1 + 1260);
    const auto *skk1_1263 = buffer.data(skk1 + 1263);
    const auto *skk1_1265 = buffer.data(skk1 + 1265);
    const auto *skk1_1266 = buffer.data(skk1 + 1266);
    const auto *skk1_1269 = buffer.data(skk1 + 1269);
    const auto *skk1_1270 = buffer.data(skk1 + 1270);
    const auto *skk1_1272 = buffer.data(skk1 + 1272);
    const auto *skk1_1274 = buffer.data(skk1 + 1274);
    const auto *skk1_1275 = buffer.data(skk1 + 1275);
    const auto *skk1_1277 = buffer.data(skk1 + 1277);
    const auto *skk1_1278 = buffer.data(skk1 + 1278);
    const auto *skk1_1280 = buffer.data(skk1 + 1280);
    const auto *skk1_1288 = buffer.data(skk1 + 1288);
    const auto *skk1_1290 = buffer.data(skk1 + 1290);
    const auto *skk1_1291 = buffer.data(skk1 + 1291);
    const auto *skk1_1292 = buffer.data(skk1 + 1292);
    const auto *skk1_1293 = buffer.data(skk1 + 1293);
    const auto *skk1_1295 = buffer.data(skk1 + 1295);

    const auto *slh0_756 = buffer.data(slh0 + 756);
    const auto *slh0_759 = buffer.data(slh0 + 759);
    const auto *slh0_761 = buffer.data(slh0 + 761);
    const auto *slh0_762 = buffer.data(slh0 + 762);
    const auto *slh0_765 = buffer.data(slh0 + 765);
    const auto *slh0_766 = buffer.data(slh0 + 766);
    const auto *slh0_768 = buffer.data(slh0 + 768);
    const auto *slh0_770 = buffer.data(slh0 + 770);
    const auto *slh0_771 = buffer.data(slh0 + 771);
    const auto *slh0_773 = buffer.data(slh0 + 773);
    const auto *slh0_774 = buffer.data(slh0 + 774);
    const auto *slh0_775 = buffer.data(slh0 + 775);
    const auto *slh0_776 = buffer.data(slh0 + 776);
    const auto *slh0_782 = buffer.data(slh0 + 782);
    const auto *slh0_786 = buffer.data(slh0 + 786);
    const auto *slh0_789 = buffer.data(slh0 + 789);
    const auto *slh0_791 = buffer.data(slh0 + 791);
    const auto *slh0_794 = buffer.data(slh0 + 794);
    const auto *slh0_795 = buffer.data(slh0 + 795);
    const auto *slh0_797 = buffer.data(slh0 + 797);
    const auto *slh0_798 = buffer.data(slh0 + 798);

    const auto *slh1_756 = buffer.data(slh1 + 756);
    const auto *slh1_759 = buffer.data(slh1 + 759);
    const auto *slh1_761 = buffer.data(slh1 + 761);
    const auto *slh1_762 = buffer.data(slh1 + 762);
    const auto *slh1_765 = buffer.data(slh1 + 765);
    const auto *slh1_766 = buffer.data(slh1 + 766);
    const auto *slh1_768 = buffer.data(slh1 + 768);
    const auto *slh1_770 = buffer.data(slh1 + 770);
    const auto *slh1_771 = buffer.data(slh1 + 771);
    const auto *slh1_773 = buffer.data(slh1 + 773);
    const auto *slh1_774 = buffer.data(slh1 + 774);
    const auto *slh1_775 = buffer.data(slh1 + 775);
    const auto *slh1_776 = buffer.data(slh1 + 776);
    const auto *slh1_782 = buffer.data(slh1 + 782);
    const auto *slh1_786 = buffer.data(slh1 + 786);
    const auto *slh1_789 = buffer.data(slh1 + 789);
    const auto *slh1_791 = buffer.data(slh1 + 791);
    const auto *slh1_794 = buffer.data(slh1 + 794);
    const auto *slh1_795 = buffer.data(slh1 + 795);
    const auto *slh1_797 = buffer.data(slh1 + 797);
    const auto *slh1_798 = buffer.data(slh1 + 798);

    const auto *sli_973 = buffer.data(sli + 973);
    const auto *sli_974 = buffer.data(sli + 974);
    const auto *sli_975 = buffer.data(sli + 975);
    const auto *sli_976 = buffer.data(sli + 976);
    const auto *sli_977 = buffer.data(sli + 977);
    const auto *sli_978 = buffer.data(sli + 978);
    const auto *sli_979 = buffer.data(sli + 979);
    const auto *sli_980 = buffer.data(sli + 980);
    const auto *sli_982 = buffer.data(sli + 982);
    const auto *sli_983 = buffer.data(sli + 983);
    const auto *sli_985 = buffer.data(sli + 985);
    const auto *sli_986 = buffer.data(sli + 986);
    const auto *sli_989 = buffer.data(sli + 989);
    const auto *sli_990 = buffer.data(sli + 990);
    const auto *sli_994 = buffer.data(sli + 994);
    const auto *sli_1001 = buffer.data(sli + 1001);
    const auto *sli_1002 = buffer.data(sli + 1002);
    const auto *sli_1003 = buffer.data(sli + 1003);
    const auto *sli_1004 = buffer.data(sli + 1004);
    const auto *sli_1005 = buffer.data(sli + 1005);
    const auto *sli_1006 = buffer.data(sli + 1006);
    const auto *sli_1007 = buffer.data(sli + 1007);
    const auto *sli_1008 = buffer.data(sli + 1008);
    const auto *sli_1010 = buffer.data(sli + 1010);
    const auto *sli_1011 = buffer.data(sli + 1011);
    const auto *sli_1013 = buffer.data(sli + 1013);
    const auto *sli_1014 = buffer.data(sli + 1014);
    const auto *sli_1017 = buffer.data(sli + 1017);
    const auto *sli_1018 = buffer.data(sli + 1018);
    const auto *sli_1020 = buffer.data(sli + 1020);
    const auto *sli_1022 = buffer.data(sli + 1022);
    const auto *sli_1023 = buffer.data(sli + 1023);
    const auto *sli_1025 = buffer.data(sli + 1025);
    const auto *sli_1026 = buffer.data(sli + 1026);
    const auto *sli_1028 = buffer.data(sli + 1028);
    const auto *sli_1029 = buffer.data(sli + 1029);
    const auto *sli_1030 = buffer.data(sli + 1030);
    const auto *sli_1031 = buffer.data(sli + 1031);
    const auto *sli_1032 = buffer.data(sli + 1032);
    const auto *sli_1033 = buffer.data(sli + 1033);
    const auto *sli_1034 = buffer.data(sli + 1034);
    const auto *sli_1035 = buffer.data(sli + 1035);
    const auto *sli_1036 = buffer.data(sli + 1036);
    const auto *sli_1038 = buffer.data(sli + 1038);
    const auto *sli_1039 = buffer.data(sli + 1039);
    const auto *sli_1041 = buffer.data(sli + 1041);
    const auto *sli_1042 = buffer.data(sli + 1042);
    const auto *sli_1045 = buffer.data(sli + 1045);
    const auto *sli_1046 = buffer.data(sli + 1046);
    const auto *sli_1048 = buffer.data(sli + 1048);
    const auto *sli_1050 = buffer.data(sli + 1050);
    const auto *sli_1053 = buffer.data(sli + 1053);
    const auto *sli_1054 = buffer.data(sli + 1054);
    const auto *sli_1056 = buffer.data(sli + 1056);
    const auto *sli_1057 = buffer.data(sli + 1057);
    const auto *sli_1058 = buffer.data(sli + 1058);
    const auto *sli_1059 = buffer.data(sli + 1059);
    const auto *sli_1060 = buffer.data(sli + 1060);
    const auto *sli_1061 = buffer.data(sli + 1061);
    const auto *sli_1062 = buffer.data(sli + 1062);
    const auto *sli_1063 = buffer.data(sli + 1063);
    const auto *sli_1064 = buffer.data(sli + 1064);

#pragma omp simd aligned(t_1245, t_1246, t_1247, t_1248, t_1249, pc_x, ski_973, ski_974, \
                         ski_975, ski_976, ski_977, sli_973, sli_974, sli_975, sli_976, \
                         sli_977 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1245[k] = f_13 * ski_973[k]
                    + f_3 * pc_x[k] * sli_973[k];

        t_1246[k] = f_13 * ski_974[k]
                    + f_3 * pc_x[k] * sli_974[k];

        t_1247[k] = f_13 * ski_975[k]
                    + f_3 * pc_x[k] * sli_975[k];

        t_1248[k] = f_13 * ski_976[k]
                    + f_3 * pc_x[k] * sli_976[k];

        t_1249[k] = f_13 * ski_977[k]
                    + f_3 * pc_x[k] * sli_977[k];
    }

#pragma omp simd aligned(t_1250, t_1251, t_1252, t_1253, pb_x, pc_x, pc_z, skk0_1252, ski_749, \
                         ski_978, ski_979, skk1_1252, sli_973, sli_978, \
                         sli_979 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1250[k] = f_13 * ski_978[k]
                    + f_3 * pc_x[k] * sli_978[k];

        t_1251[k] = f_13 * ski_979[k]
                    + f_3 * pc_x[k] * sli_979[k];

        t_1252[k] = pb_x[k] * skk0_1252[k]
                    - f_12 * pc_x[k] * skk1_1252[k];

        t_1253[k] = f_19 * ski_749[k]
                    + f_3 * pc_z[k] * sli_973[k];
    }

#pragma omp simd aligned(t_1254, t_1255, t_1256, t_1257, pb_x, pc_x, skk0_1254, skk0_1255, \
                         skk0_1256, skk0_1257, skk1_1254, skk1_1255, skk1_1256, \
                         skk1_1257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1254[k] = pb_x[k] * skk0_1254[k]
                    - f_12 * pc_x[k] * skk1_1254[k];

        t_1255[k] = pb_x[k] * skk0_1255[k]
                    - f_12 * pc_x[k] * skk1_1255[k];

        t_1256[k] = pb_x[k] * skk0_1256[k]
                    - f_12 * pc_x[k] * skk1_1256[k];

        t_1257[k] = pb_x[k] * skk0_1257[k]
                    - f_12 * pc_x[k] * skk1_1257[k];
    }

#pragma omp simd aligned(t_1258, t_1259, t_1260, t_1261, pb_x, pc_x, pc_y, skk0_1259, \
                         skk0_1260, ski_783, ski_980, skk1_1259, skk1_1260, sli_979, \
                         sli_980 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1258[k] = f_13 * ski_783[k]
                    + f_3 * pc_y[k] * sli_979[k];

        t_1259[k] = pb_x[k] * skk0_1259[k]
                    - f_12 * pc_x[k] * skk1_1259[k];

        t_1260[k] = pb_x[k] * skk0_1260[k]
                    + f_18 * ski_980[k]
                    - f_12 * pc_x[k] * skk1_1260[k];

        t_1261[k] = f_3 * pc_y[k] * sli_980[k];
    }

#pragma omp simd aligned(t_1262, t_1263, t_1264, pb_x, pc_x, pc_y, pc_z, skk0_1263, ski_756, \
                         ski_983, skk1_1263, sli_980, sli_982 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1262[k] = f_18 * ski_756[k]
                    + f_3 * pc_z[k] * sli_980[k];

        t_1263[k] = pb_x[k] * skk0_1263[k]
                    + f_17 * ski_983[k]
                    - f_12 * pc_x[k] * skk1_1263[k];

        t_1264[k] = f_3 * pc_y[k] * sli_982[k];
    }

#pragma omp simd aligned(t_1265, t_1266, t_1267, pb_x, pc_x, pc_z, skk0_1265, skk0_1266, \
                         ski_759, ski_985, ski_986, skk1_1265, skk1_1266, \
                         sli_983 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1265[k] = pb_x[k] * skk0_1265[k]
                    + f_17 * ski_985[k]
                    - f_12 * pc_x[k] * skk1_1265[k];

        t_1266[k] = pb_x[k] * skk0_1266[k]
                    + f_16 * ski_986[k]
                    - f_12 * pc_x[k] * skk1_1266[k];

        t_1267[k] = f_18 * ski_759[k]
                    + f_3 * pc_z[k] * sli_983[k];
    }

#pragma omp simd aligned(t_1268, t_1269, t_1270, pb_x, pc_x, pc_y, skk0_1269, skk0_1270, \
                         ski_989, ski_990, skk1_1269, skk1_1270, \
                         sli_985 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1268[k] = f_3 * pc_y[k] * sli_985[k];

        t_1269[k] = pb_x[k] * skk0_1269[k]
                    + f_16 * ski_989[k]
                    - f_12 * pc_x[k] * skk1_1269[k];

        t_1270[k] = pb_x[k] * skk0_1270[k]
                    + f_15 * ski_990[k]
                    - f_12 * pc_x[k] * skk1_1270[k];
    }

#pragma omp simd aligned(t_1271, t_1272, t_1273, pb_x, pc_x, pc_y, pc_z, skk0_1272, ski_762, \
                         ski_992, skk1_1272, sli_986, sli_989 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1271[k] = f_18 * ski_762[k]
                    + f_3 * pc_z[k] * sli_986[k];

        t_1272[k] = pb_x[k] * skk0_1272[k]
                    + f_15 * ski_992[k]
                    - f_12 * pc_x[k] * skk1_1272[k];

        t_1273[k] = f_3 * pc_y[k] * sli_989[k];
    }

#pragma omp simd aligned(t_1274, t_1275, t_1276, pb_x, pc_x, pc_z, skk0_1274, skk0_1275, \
                         ski_766, ski_994, ski_995, skk1_1274, skk1_1275, \
                         sli_990 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1274[k] = pb_x[k] * skk0_1274[k]
                    + f_15 * ski_994[k]
                    - f_12 * pc_x[k] * skk1_1274[k];

        t_1275[k] = pb_x[k] * skk0_1275[k]
                    + f_14 * ski_995[k]
                    - f_12 * pc_x[k] * skk1_1275[k];

        t_1276[k] = f_18 * ski_766[k]
                    + f_3 * pc_z[k] * sli_990[k];
    }

#pragma omp simd aligned(t_1277, t_1278, t_1279, pb_x, pc_x, pc_y, skk0_1277, skk0_1278, \
                         ski_997, ski_998, skk1_1277, skk1_1278, \
                         sli_994 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1277[k] = pb_x[k] * skk0_1277[k]
                    + f_14 * ski_997[k]
                    - f_12 * pc_x[k] * skk1_1277[k];

        t_1278[k] = pb_x[k] * skk0_1278[k]
                    + f_14 * ski_998[k]
                    - f_12 * pc_x[k] * skk1_1278[k];

        t_1279[k] = f_3 * pc_y[k] * sli_994[k];
    }

#pragma omp simd aligned(t_1280, t_1281, t_1282, t_1283, pb_x, pc_x, skk0_1280, ski_1000, \
                         ski_1001, ski_1002, ski_1003, skk1_1280, sli_1001, sli_1002, \
                         sli_1003 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1280[k] = pb_x[k] * skk0_1280[k]
                    + f_14 * ski_1000[k]
                    - f_12 * pc_x[k] * skk1_1280[k];

        t_1281[k] = f_13 * ski_1001[k]
                    + f_3 * pc_x[k] * sli_1001[k];

        t_1282[k] = f_13 * ski_1002[k]
                    + f_3 * pc_x[k] * sli_1002[k];

        t_1283[k] = f_13 * ski_1003[k]
                    + f_3 * pc_x[k] * sli_1003[k];
    }

#pragma omp simd aligned(t_1284, t_1285, t_1286, t_1287, pc_x, ski_1004, ski_1005, ski_1006, \
                         ski_1007, sli_1004, sli_1005, sli_1006, \
                         sli_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1284[k] = f_13 * ski_1004[k]
                    + f_3 * pc_x[k] * sli_1004[k];

        t_1285[k] = f_13 * ski_1005[k]
                    + f_3 * pc_x[k] * sli_1005[k];

        t_1286[k] = f_13 * ski_1006[k]
                    + f_3 * pc_x[k] * sli_1006[k];

        t_1287[k] = f_13 * ski_1007[k]
                    + f_3 * pc_x[k] * sli_1007[k];
    }

#pragma omp simd aligned(t_1288, t_1289, t_1290, t_1291, pb_x, pc_x, pc_z, skk0_1288, \
                         skk0_1290, skk0_1291, ski_777, skk1_1288, skk1_1290, skk1_1291, \
                         sli_1001 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1288[k] = pb_x[k] * skk0_1288[k]
                    - f_12 * pc_x[k] * skk1_1288[k];

        t_1289[k] = f_18 * ski_777[k]
                    + f_3 * pc_z[k] * sli_1001[k];

        t_1290[k] = pb_x[k] * skk0_1290[k]
                    - f_12 * pc_x[k] * skk1_1290[k];

        t_1291[k] = pb_x[k] * skk0_1291[k]
                    - f_12 * pc_x[k] * skk1_1291[k];
    }

#pragma omp simd aligned(t_1292, t_1293, t_1294, t_1295, pb_x, pc_x, pc_y, skk0_1292, \
                         skk0_1293, skk0_1295, skk1_1292, skk1_1293, skk1_1295, \
                         sli_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1292[k] = pb_x[k] * skk0_1292[k]
                    - f_12 * pc_x[k] * skk1_1292[k];

        t_1293[k] = pb_x[k] * skk0_1293[k]
                    - f_12 * pc_x[k] * skk1_1293[k];

        t_1294[k] = f_3 * pc_y[k] * sli_1007[k];

        t_1295[k] = pb_x[k] * skk0_1295[k]
                    - f_12 * pc_x[k] * skk1_1295[k];
    }

#pragma omp simd aligned(t_1296, t_1297, t_1298, t_1299, pc_x, pc_y, pc_z, ski_784, slh0_756, \
                         slh0_759, slh1_756, slh1_759, sli_1008, \
                         sli_1011 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1296[k] = f_1 * slh0_756[k]
                    - f_2 * slh1_756[k]
                    + f_3 * pc_x[k] * sli_1008[k];

        t_1297[k] = f_0 * ski_784[k]
                    + f_3 * pc_y[k] * sli_1008[k];

        t_1298[k] = f_3 * pc_z[k] * sli_1008[k];

        t_1299[k] = f_4 * slh0_759[k]
                    - f_5 * slh1_759[k]
                    + f_3 * pc_x[k] * sli_1011[k];
    }

#pragma omp simd aligned(t_1300, t_1301, t_1302, t_1303, pc_x, pc_y, pc_z, ski_786, slh0_761, \
                         slh0_762, slh1_761, slh1_762, sli_1010, sli_1011, sli_1013, \
                         sli_1014 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1300[k] = f_0 * ski_786[k]
                    + f_3 * pc_y[k] * sli_1010[k];

        t_1301[k] = f_4 * slh0_761[k]
                    - f_5 * slh1_761[k]
                    + f_3 * pc_x[k] * sli_1013[k];

        t_1302[k] = f_6 * slh0_762[k]
                    - f_7 * slh1_762[k]
                    + f_3 * pc_x[k] * sli_1014[k];

        t_1303[k] = f_3 * pc_z[k] * sli_1011[k];
    }

#pragma omp simd aligned(t_1304, t_1305, t_1306, t_1307, pc_x, pc_y, pc_z, ski_789, slh0_765, \
                         slh0_766, slh1_765, slh1_766, sli_1013, sli_1014, sli_1017, \
                         sli_1018 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1304[k] = f_0 * ski_789[k]
                    + f_3 * pc_y[k] * sli_1013[k];

        t_1305[k] = f_6 * slh0_765[k]
                    - f_7 * slh1_765[k]
                    + f_3 * pc_x[k] * sli_1017[k];

        t_1306[k] = f_8 * slh0_766[k]
                    - f_9 * slh1_766[k]
                    + f_3 * pc_x[k] * sli_1018[k];

        t_1307[k] = f_3 * pc_z[k] * sli_1014[k];
    }

#pragma omp simd aligned(t_1308, t_1309, t_1310, pc_x, pc_y, ski_793, slh0_768, slh0_770, \
                         slh1_768, slh1_770, sli_1017, sli_1020, \
                         sli_1022 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1308[k] = f_8 * slh0_768[k]
                    - f_9 * slh1_768[k]
                    + f_3 * pc_x[k] * sli_1020[k];

        t_1309[k] = f_0 * ski_793[k]
                    + f_3 * pc_y[k] * sli_1017[k];

        t_1310[k] = f_8 * slh0_770[k]
                    - f_9 * slh1_770[k]
                    + f_3 * pc_x[k] * sli_1022[k];
    }

#pragma omp simd aligned(t_1311, t_1312, t_1313, t_1314, pc_x, pc_z, slh0_771, slh0_773, \
                         slh0_774, slh1_771, slh1_773, slh1_774, sli_1018, sli_1023, sli_1025, \
                         sli_1026 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1311[k] = f_10 * slh0_771[k]
                    - f_11 * slh1_771[k]
                    + f_3 * pc_x[k] * sli_1023[k];

        t_1312[k] = f_3 * pc_z[k] * sli_1018[k];

        t_1313[k] = f_10 * slh0_773[k]
                    - f_11 * slh1_773[k]
                    + f_3 * pc_x[k] * sli_1025[k];

        t_1314[k] = f_10 * slh0_774[k]
                    - f_11 * slh1_774[k]
                    + f_3 * pc_x[k] * sli_1026[k];
    }

#pragma omp simd aligned(t_1315, t_1316, t_1317, t_1318, t_1319, pc_x, pc_y, ski_798, \
                         slh0_776, slh1_776, sli_1022, sli_1028, sli_1029, sli_1030, \
                         sli_1031 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1315[k] = f_0 * ski_798[k]
                    + f_3 * pc_y[k] * sli_1022[k];

        t_1316[k] = f_10 * slh0_776[k]
                    - f_11 * slh1_776[k]
                    + f_3 * pc_x[k] * sli_1028[k];

        t_1317[k] = f_3 * pc_x[k] * sli_1029[k];

        t_1318[k] = f_3 * pc_x[k] * sli_1030[k];

        t_1319[k] = f_3 * pc_x[k] * sli_1031[k];
    }

#pragma omp simd aligned(t_1320, t_1321, t_1322, t_1323, t_1324, pc_x, pc_y, ski_805, \
                         slh0_771, slh1_771, sli_1029, sli_1032, sli_1033, sli_1034, \
                         sli_1035 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1320[k] = f_3 * pc_x[k] * sli_1032[k];

        t_1321[k] = f_3 * pc_x[k] * sli_1033[k];

        t_1322[k] = f_3 * pc_x[k] * sli_1034[k];

        t_1323[k] = f_3 * pc_x[k] * sli_1035[k];

        t_1324[k] = f_0 * ski_805[k]
                    + f_1 * slh0_771[k]
                    - f_2 * slh1_771[k]
                    + f_3 * pc_y[k] * sli_1029[k];
    }

#pragma omp simd aligned(t_1325, t_1326, t_1327, pc_y, pc_z, ski_807, ski_808, slh0_773, \
                         slh0_774, slh1_773, slh1_774, sli_1029, sli_1031, \
                         sli_1032 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1325[k] = f_3 * pc_z[k] * sli_1029[k];

        t_1326[k] = f_0 * ski_807[k]
                    + f_4 * slh0_773[k]
                    - f_5 * slh1_773[k]
                    + f_3 * pc_y[k] * sli_1031[k];

        t_1327[k] = f_0 * ski_808[k]
                    + f_6 * slh0_774[k]
                    - f_7 * slh1_774[k]
                    + f_3 * pc_y[k] * sli_1032[k];
    }

#pragma omp simd aligned(t_1328, t_1329, t_1330, t_1331, pc_y, pc_z, ski_809, ski_810, \
                         ski_811, slh0_775, slh0_776, slh1_775, slh1_776, sli_1033, sli_1034, \
                         sli_1035 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1328[k] = f_0 * ski_809[k]
                    + f_8 * slh0_775[k]
                    - f_9 * slh1_775[k]
                    + f_3 * pc_y[k] * sli_1033[k];

        t_1329[k] = f_0 * ski_810[k]
                    + f_10 * slh0_776[k]
                    - f_11 * slh1_776[k]
                    + f_3 * pc_y[k] * sli_1034[k];

        t_1330[k] = f_0 * ski_811[k]
                    + f_3 * pc_y[k] * sli_1035[k];

        t_1331[k] = f_1 * slh0_776[k]
                    - f_2 * slh1_776[k]
                    + f_3 * pc_z[k] * sli_1035[k];
    }

#pragma omp simd aligned(t_1332, t_1333, t_1334, t_1335, pb_z, pc_y, pc_z, skk0_1008, \
                         skk0_1011, ski_784, ski_812, skk1_1008, skk1_1011, \
                         sli_1036 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1332[k] = pb_z[k] * skk0_1008[k]
                    - f_12 * pc_z[k] * skk1_1008[k];

        t_1333[k] = f_18 * ski_812[k]
                    + f_3 * pc_y[k] * sli_1036[k];

        t_1334[k] = f_13 * ski_784[k]
                    + f_3 * pc_z[k] * sli_1036[k];

        t_1335[k] = pb_z[k] * skk0_1011[k]
                    - f_12 * pc_z[k] * skk1_1011[k];
    }

#pragma omp simd aligned(t_1336, t_1337, t_1338, pb_z, pc_x, pc_y, pc_z, skk0_1014, ski_814, \
                         skk1_1014, slh0_782, slh1_782, sli_1038, \
                         sli_1041 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1336[k] = f_18 * ski_814[k]
                    + f_3 * pc_y[k] * sli_1038[k];

        t_1337[k] = f_4 * slh0_782[k]
                    - f_5 * slh1_782[k]
                    + f_3 * pc_x[k] * sli_1041[k];

        t_1338[k] = pb_z[k] * skk0_1014[k]
                    - f_12 * pc_z[k] * skk1_1014[k];
    }

#pragma omp simd aligned(t_1339, t_1340, t_1341, pc_x, pc_y, pc_z, ski_787, ski_817, slh0_786, \
                         slh1_786, sli_1039, sli_1041, sli_1045 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1339[k] = f_13 * ski_787[k]
                    + f_3 * pc_z[k] * sli_1039[k];

        t_1340[k] = f_18 * ski_817[k]
                    + f_3 * pc_y[k] * sli_1041[k];

        t_1341[k] = f_6 * slh0_786[k]
                    - f_7 * slh1_786[k]
                    + f_3 * pc_x[k] * sli_1045[k];
    }

#pragma omp simd aligned(t_1342, t_1343, t_1344, pb_z, pc_x, pc_z, skk0_1018, ski_790, \
                         skk1_1018, slh0_789, slh1_789, sli_1042, \
                         sli_1048 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1342[k] = pb_z[k] * skk0_1018[k]
                    - f_12 * pc_z[k] * skk1_1018[k];

        t_1343[k] = f_13 * ski_790[k]
                    + f_3 * pc_z[k] * sli_1042[k];

        t_1344[k] = f_8 * slh0_789[k]
                    - f_9 * slh1_789[k]
                    + f_3 * pc_x[k] * sli_1048[k];
    }

#pragma omp simd aligned(t_1345, t_1346, t_1347, pb_z, pc_x, pc_y, pc_z, skk0_1023, ski_821, \
                         skk1_1023, slh0_791, slh1_791, sli_1045, \
                         sli_1050 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1345[k] = f_18 * ski_821[k]
                    + f_3 * pc_y[k] * sli_1045[k];

        t_1346[k] = f_8 * slh0_791[k]
                    - f_9 * slh1_791[k]
                    + f_3 * pc_x[k] * sli_1050[k];

        t_1347[k] = pb_z[k] * skk0_1023[k]
                    - f_12 * pc_z[k] * skk1_1023[k];
    }

#pragma omp simd aligned(t_1348, t_1349, t_1350, pc_x, pc_z, ski_794, slh0_794, slh0_795, \
                         slh1_794, slh1_795, sli_1046, sli_1053, \
                         sli_1054 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1348[k] = f_13 * ski_794[k]
                    + f_3 * pc_z[k] * sli_1046[k];

        t_1349[k] = f_10 * slh0_794[k]
                    - f_11 * slh1_794[k]
                    + f_3 * pc_x[k] * sli_1053[k];

        t_1350[k] = f_10 * slh0_795[k]
                    - f_11 * slh1_795[k]
                    + f_3 * pc_x[k] * sli_1054[k];
    }

#pragma omp simd aligned(t_1351, t_1352, t_1353, t_1354, t_1355, pc_x, pc_y, ski_826, \
                         slh0_797, slh1_797, sli_1050, sli_1056, sli_1057, sli_1058, \
                         sli_1059 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1351[k] = f_18 * ski_826[k]
                    + f_3 * pc_y[k] * sli_1050[k];

        t_1352[k] = f_10 * slh0_797[k]
                    - f_11 * slh1_797[k]
                    + f_3 * pc_x[k] * sli_1056[k];

        t_1353[k] = f_3 * pc_x[k] * sli_1057[k];

        t_1354[k] = f_3 * pc_x[k] * sli_1058[k];

        t_1355[k] = f_3 * pc_x[k] * sli_1059[k];
    }

#pragma omp simd aligned(t_1356, t_1357, t_1358, t_1359, t_1360, pb_z, pc_x, pc_z, skk0_1036, \
                         skk1_1036, sli_1060, sli_1061, sli_1062, \
                         sli_1063 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1356[k] = f_3 * pc_x[k] * sli_1060[k];

        t_1357[k] = f_3 * pc_x[k] * sli_1061[k];

        t_1358[k] = f_3 * pc_x[k] * sli_1062[k];

        t_1359[k] = f_3 * pc_x[k] * sli_1063[k];

        t_1360[k] = pb_z[k] * skk0_1036[k]
                    - f_12 * pc_z[k] * skk1_1036[k];
    }

#pragma omp simd aligned(t_1361, t_1362, t_1363, pb_z, pc_z, skk0_1038, skk0_1039, ski_805, \
                         ski_806, ski_807, skk1_1038, skk1_1039, \
                         sli_1057 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1361[k] = f_13 * ski_805[k]
                    + f_3 * pc_z[k] * sli_1057[k];

        t_1362[k] = pb_z[k] * skk0_1038[k]
                    + f_14 * ski_806[k]
                    - f_12 * pc_z[k] * skk1_1038[k];

        t_1363[k] = pb_z[k] * skk0_1039[k]
                    + f_15 * ski_807[k]
                    - f_12 * pc_z[k] * skk1_1039[k];
    }

#pragma omp simd aligned(t_1364, t_1365, t_1366, pb_z, pc_y, pc_z, skk0_1040, skk0_1041, \
                         ski_808, ski_809, ski_839, skk1_1040, skk1_1041, \
                         sli_1063 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1364[k] = pb_z[k] * skk0_1040[k]
                    + f_16 * ski_808[k]
                    - f_12 * pc_z[k] * skk1_1040[k];

        t_1365[k] = pb_z[k] * skk0_1041[k]
                    + f_17 * ski_809[k]
                    - f_12 * pc_z[k] * skk1_1041[k];

        t_1366[k] = f_18 * ski_839[k]
                    + f_3 * pc_y[k] * sli_1063[k];
    }

#pragma omp simd aligned(t_1367, t_1368, t_1369, t_1370, pc_x, pc_y, pc_z, ski_811, ski_812, \
                         ski_840, slh0_797, slh0_798, slh1_797, slh1_798, sli_1063, \
                         sli_1064 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1367[k] = f_13 * ski_811[k]
                    + f_1 * slh0_797[k]
                    - f_2 * slh1_797[k]
                    + f_3 * pc_z[k] * sli_1063[k];

        t_1368[k] = f_1 * slh0_798[k]
                    - f_2 * slh1_798[k]
                    + f_3 * pc_x[k] * sli_1064[k];

        t_1369[k] = f_19 * ski_840[k]
                    + f_3 * pc_y[k] * sli_1064[k];

        t_1370[k] = f_14 * ski_812[k]
                    + f_3 * pc_z[k] * sli_1064[k];
    }
}

static auto
compute_prim_slk_three_center_electron_repulsion_0_piece12(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pc,
                                                           const size_t ski, const size_t slh0,
                                                           const size_t slh1, const size_t sli,
                                                           const size_t ncols,
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
    const auto f_14 = 1.0 / q;
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;
    const auto f_19 = 3.0 / q;

    auto *t_1371 = buffer.data(target + 1371);
    auto *t_1372 = buffer.data(target + 1372);
    auto *t_1373 = buffer.data(target + 1373);
    auto *t_1374 = buffer.data(target + 1374);
    auto *t_1375 = buffer.data(target + 1375);
    auto *t_1376 = buffer.data(target + 1376);
    auto *t_1377 = buffer.data(target + 1377);
    auto *t_1378 = buffer.data(target + 1378);
    auto *t_1379 = buffer.data(target + 1379);
    auto *t_1380 = buffer.data(target + 1380);
    auto *t_1381 = buffer.data(target + 1381);
    auto *t_1382 = buffer.data(target + 1382);
    auto *t_1383 = buffer.data(target + 1383);
    auto *t_1384 = buffer.data(target + 1384);
    auto *t_1385 = buffer.data(target + 1385);
    auto *t_1386 = buffer.data(target + 1386);
    auto *t_1387 = buffer.data(target + 1387);
    auto *t_1388 = buffer.data(target + 1388);
    auto *t_1389 = buffer.data(target + 1389);
    auto *t_1390 = buffer.data(target + 1390);
    auto *t_1391 = buffer.data(target + 1391);
    auto *t_1392 = buffer.data(target + 1392);
    auto *t_1393 = buffer.data(target + 1393);
    auto *t_1394 = buffer.data(target + 1394);
    auto *t_1395 = buffer.data(target + 1395);
    auto *t_1396 = buffer.data(target + 1396);
    auto *t_1397 = buffer.data(target + 1397);
    auto *t_1398 = buffer.data(target + 1398);
    auto *t_1399 = buffer.data(target + 1399);
    auto *t_1400 = buffer.data(target + 1400);
    auto *t_1401 = buffer.data(target + 1401);
    auto *t_1402 = buffer.data(target + 1402);
    auto *t_1403 = buffer.data(target + 1403);
    auto *t_1404 = buffer.data(target + 1404);
    auto *t_1405 = buffer.data(target + 1405);
    auto *t_1406 = buffer.data(target + 1406);
    auto *t_1407 = buffer.data(target + 1407);
    auto *t_1408 = buffer.data(target + 1408);
    auto *t_1409 = buffer.data(target + 1409);
    auto *t_1410 = buffer.data(target + 1410);
    auto *t_1411 = buffer.data(target + 1411);
    auto *t_1412 = buffer.data(target + 1412);
    auto *t_1413 = buffer.data(target + 1413);
    auto *t_1414 = buffer.data(target + 1414);
    auto *t_1415 = buffer.data(target + 1415);
    auto *t_1416 = buffer.data(target + 1416);
    auto *t_1417 = buffer.data(target + 1417);
    auto *t_1418 = buffer.data(target + 1418);
    auto *t_1419 = buffer.data(target + 1419);
    auto *t_1420 = buffer.data(target + 1420);
    auto *t_1421 = buffer.data(target + 1421);
    auto *t_1422 = buffer.data(target + 1422);
    auto *t_1423 = buffer.data(target + 1423);
    auto *t_1424 = buffer.data(target + 1424);
    auto *t_1425 = buffer.data(target + 1425);
    auto *t_1426 = buffer.data(target + 1426);
    auto *t_1427 = buffer.data(target + 1427);
    auto *t_1428 = buffer.data(target + 1428);
    auto *t_1429 = buffer.data(target + 1429);
    auto *t_1430 = buffer.data(target + 1430);
    auto *t_1431 = buffer.data(target + 1431);
    auto *t_1432 = buffer.data(target + 1432);
    auto *t_1433 = buffer.data(target + 1433);
    auto *t_1434 = buffer.data(target + 1434);
    auto *t_1435 = buffer.data(target + 1435);
    auto *t_1436 = buffer.data(target + 1436);
    auto *t_1437 = buffer.data(target + 1437);
    auto *t_1438 = buffer.data(target + 1438);
    auto *t_1439 = buffer.data(target + 1439);
    auto *t_1440 = buffer.data(target + 1440);
    auto *t_1441 = buffer.data(target + 1441);
    auto *t_1442 = buffer.data(target + 1442);
    auto *t_1443 = buffer.data(target + 1443);
    auto *t_1444 = buffer.data(target + 1444);
    auto *t_1445 = buffer.data(target + 1445);
    auto *t_1446 = buffer.data(target + 1446);
    auto *t_1447 = buffer.data(target + 1447);
    auto *t_1448 = buffer.data(target + 1448);
    auto *t_1449 = buffer.data(target + 1449);
    auto *t_1450 = buffer.data(target + 1450);
    auto *t_1451 = buffer.data(target + 1451);
    auto *t_1452 = buffer.data(target + 1452);
    auto *t_1453 = buffer.data(target + 1453);
    auto *t_1454 = buffer.data(target + 1454);
    auto *t_1455 = buffer.data(target + 1455);
    auto *t_1456 = buffer.data(target + 1456);
    auto *t_1457 = buffer.data(target + 1457);
    auto *t_1458 = buffer.data(target + 1458);
    auto *t_1459 = buffer.data(target + 1459);
    auto *t_1460 = buffer.data(target + 1460);
    auto *t_1461 = buffer.data(target + 1461);
    auto *t_1462 = buffer.data(target + 1462);
    auto *t_1463 = buffer.data(target + 1463);
    auto *t_1464 = buffer.data(target + 1464);
    auto *t_1465 = buffer.data(target + 1465);
    auto *t_1466 = buffer.data(target + 1466);
    auto *t_1467 = buffer.data(target + 1467);
    auto *t_1468 = buffer.data(target + 1468);
    auto *t_1469 = buffer.data(target + 1469);
    auto *t_1470 = buffer.data(target + 1470);
    auto *t_1471 = buffer.data(target + 1471);
    auto *t_1472 = buffer.data(target + 1472);
    auto *t_1473 = buffer.data(target + 1473);
    auto *t_1474 = buffer.data(target + 1474);
    auto *t_1475 = buffer.data(target + 1475);
    auto *t_1476 = buffer.data(target + 1476);
    auto *t_1477 = buffer.data(target + 1477);
    auto *t_1478 = buffer.data(target + 1478);
    auto *t_1479 = buffer.data(target + 1479);
    auto *t_1480 = buffer.data(target + 1480);
    auto *t_1481 = buffer.data(target + 1481);
    auto *t_1482 = buffer.data(target + 1482);
    auto *t_1483 = buffer.data(target + 1483);
    auto *t_1484 = buffer.data(target + 1484);
    auto *t_1485 = buffer.data(target + 1485);
    auto *t_1486 = buffer.data(target + 1486);
    auto *t_1487 = buffer.data(target + 1487);
    auto *t_1488 = buffer.data(target + 1488);
    auto *t_1489 = buffer.data(target + 1489);
    auto *t_1490 = buffer.data(target + 1490);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ski_815 = buffer.data(ski + 815);
    const auto *ski_818 = buffer.data(ski + 818);
    const auto *ski_822 = buffer.data(ski + 822);
    const auto *ski_833 = buffer.data(ski + 833);
    const auto *ski_839 = buffer.data(ski + 839);
    const auto *ski_840 = buffer.data(ski + 840);
    const auto *ski_842 = buffer.data(ski + 842);
    const auto *ski_843 = buffer.data(ski + 843);
    const auto *ski_845 = buffer.data(ski + 845);
    const auto *ski_846 = buffer.data(ski + 846);
    const auto *ski_849 = buffer.data(ski + 849);
    const auto *ski_850 = buffer.data(ski + 850);
    const auto *ski_854 = buffer.data(ski + 854);
    const auto *ski_861 = buffer.data(ski + 861);
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
    const auto *ski_882 = buffer.data(ski + 882);
    const auto *ski_889 = buffer.data(ski + 889);
    const auto *ski_891 = buffer.data(ski + 891);
    const auto *ski_892 = buffer.data(ski + 892);
    const auto *ski_893 = buffer.data(ski + 893);
    const auto *ski_894 = buffer.data(ski + 894);
    const auto *ski_895 = buffer.data(ski + 895);
    const auto *ski_896 = buffer.data(ski + 896);
    const auto *ski_898 = buffer.data(ski + 898);
    const auto *ski_899 = buffer.data(ski + 899);
    const auto *ski_901 = buffer.data(ski + 901);
    const auto *ski_902 = buffer.data(ski + 902);
    const auto *ski_905 = buffer.data(ski + 905);
    const auto *ski_910 = buffer.data(ski + 910);
    const auto *ski_917 = buffer.data(ski + 917);
    const auto *ski_919 = buffer.data(ski + 919);
    const auto *ski_920 = buffer.data(ski + 920);
    const auto *ski_921 = buffer.data(ski + 921);
    const auto *ski_922 = buffer.data(ski + 922);
    const auto *ski_923 = buffer.data(ski + 923);
    const auto *ski_924 = buffer.data(ski + 924);
    const auto *ski_926 = buffer.data(ski + 926);
    const auto *ski_929 = buffer.data(ski + 929);
    const auto *ski_933 = buffer.data(ski + 933);

    const auto *slh0_801 = buffer.data(slh0 + 801);
    const auto *slh0_803 = buffer.data(slh0 + 803);
    const auto *slh0_804 = buffer.data(slh0 + 804);
    const auto *slh0_807 = buffer.data(slh0 + 807);
    const auto *slh0_808 = buffer.data(slh0 + 808);
    const auto *slh0_810 = buffer.data(slh0 + 810);
    const auto *slh0_812 = buffer.data(slh0 + 812);
    const auto *slh0_813 = buffer.data(slh0 + 813);
    const auto *slh0_815 = buffer.data(slh0 + 815);
    const auto *slh0_816 = buffer.data(slh0 + 816);
    const auto *slh0_817 = buffer.data(slh0 + 817);
    const auto *slh0_818 = buffer.data(slh0 + 818);
    const auto *slh0_819 = buffer.data(slh0 + 819);
    const auto *slh0_822 = buffer.data(slh0 + 822);
    const auto *slh0_824 = buffer.data(slh0 + 824);
    const auto *slh0_825 = buffer.data(slh0 + 825);
    const auto *slh0_828 = buffer.data(slh0 + 828);
    const auto *slh0_829 = buffer.data(slh0 + 829);
    const auto *slh0_831 = buffer.data(slh0 + 831);
    const auto *slh0_833 = buffer.data(slh0 + 833);
    const auto *slh0_834 = buffer.data(slh0 + 834);
    const auto *slh0_836 = buffer.data(slh0 + 836);
    const auto *slh0_837 = buffer.data(slh0 + 837);
    const auto *slh0_838 = buffer.data(slh0 + 838);
    const auto *slh0_839 = buffer.data(slh0 + 839);
    const auto *slh0_840 = buffer.data(slh0 + 840);
    const auto *slh0_843 = buffer.data(slh0 + 843);
    const auto *slh0_845 = buffer.data(slh0 + 845);
    const auto *slh0_846 = buffer.data(slh0 + 846);
    const auto *slh0_849 = buffer.data(slh0 + 849);
    const auto *slh0_850 = buffer.data(slh0 + 850);
    const auto *slh0_852 = buffer.data(slh0 + 852);
    const auto *slh0_854 = buffer.data(slh0 + 854);
    const auto *slh0_855 = buffer.data(slh0 + 855);
    const auto *slh0_857 = buffer.data(slh0 + 857);
    const auto *slh0_858 = buffer.data(slh0 + 858);
    const auto *slh0_859 = buffer.data(slh0 + 859);
    const auto *slh0_860 = buffer.data(slh0 + 860);
    const auto *slh0_861 = buffer.data(slh0 + 861);
    const auto *slh0_864 = buffer.data(slh0 + 864);
    const auto *slh0_866 = buffer.data(slh0 + 866);
    const auto *slh0_867 = buffer.data(slh0 + 867);
    const auto *slh0_870 = buffer.data(slh0 + 870);
    const auto *slh0_871 = buffer.data(slh0 + 871);
    const auto *slh0_873 = buffer.data(slh0 + 873);
    const auto *slh0_875 = buffer.data(slh0 + 875);

    const auto *slh1_801 = buffer.data(slh1 + 801);
    const auto *slh1_803 = buffer.data(slh1 + 803);
    const auto *slh1_804 = buffer.data(slh1 + 804);
    const auto *slh1_807 = buffer.data(slh1 + 807);
    const auto *slh1_808 = buffer.data(slh1 + 808);
    const auto *slh1_810 = buffer.data(slh1 + 810);
    const auto *slh1_812 = buffer.data(slh1 + 812);
    const auto *slh1_813 = buffer.data(slh1 + 813);
    const auto *slh1_815 = buffer.data(slh1 + 815);
    const auto *slh1_816 = buffer.data(slh1 + 816);
    const auto *slh1_817 = buffer.data(slh1 + 817);
    const auto *slh1_818 = buffer.data(slh1 + 818);
    const auto *slh1_819 = buffer.data(slh1 + 819);
    const auto *slh1_822 = buffer.data(slh1 + 822);
    const auto *slh1_824 = buffer.data(slh1 + 824);
    const auto *slh1_825 = buffer.data(slh1 + 825);
    const auto *slh1_828 = buffer.data(slh1 + 828);
    const auto *slh1_829 = buffer.data(slh1 + 829);
    const auto *slh1_831 = buffer.data(slh1 + 831);
    const auto *slh1_833 = buffer.data(slh1 + 833);
    const auto *slh1_834 = buffer.data(slh1 + 834);
    const auto *slh1_836 = buffer.data(slh1 + 836);
    const auto *slh1_837 = buffer.data(slh1 + 837);
    const auto *slh1_838 = buffer.data(slh1 + 838);
    const auto *slh1_839 = buffer.data(slh1 + 839);
    const auto *slh1_840 = buffer.data(slh1 + 840);
    const auto *slh1_843 = buffer.data(slh1 + 843);
    const auto *slh1_845 = buffer.data(slh1 + 845);
    const auto *slh1_846 = buffer.data(slh1 + 846);
    const auto *slh1_849 = buffer.data(slh1 + 849);
    const auto *slh1_850 = buffer.data(slh1 + 850);
    const auto *slh1_852 = buffer.data(slh1 + 852);
    const auto *slh1_854 = buffer.data(slh1 + 854);
    const auto *slh1_855 = buffer.data(slh1 + 855);
    const auto *slh1_857 = buffer.data(slh1 + 857);
    const auto *slh1_858 = buffer.data(slh1 + 858);
    const auto *slh1_859 = buffer.data(slh1 + 859);
    const auto *slh1_860 = buffer.data(slh1 + 860);
    const auto *slh1_861 = buffer.data(slh1 + 861);
    const auto *slh1_864 = buffer.data(slh1 + 864);
    const auto *slh1_866 = buffer.data(slh1 + 866);
    const auto *slh1_867 = buffer.data(slh1 + 867);
    const auto *slh1_870 = buffer.data(slh1 + 870);
    const auto *slh1_871 = buffer.data(slh1 + 871);
    const auto *slh1_873 = buffer.data(slh1 + 873);
    const auto *slh1_875 = buffer.data(slh1 + 875);

    const auto *sli_1066 = buffer.data(sli + 1066);
    const auto *sli_1067 = buffer.data(sli + 1067);
    const auto *sli_1069 = buffer.data(sli + 1069);
    const auto *sli_1070 = buffer.data(sli + 1070);
    const auto *sli_1073 = buffer.data(sli + 1073);
    const auto *sli_1074 = buffer.data(sli + 1074);
    const auto *sli_1076 = buffer.data(sli + 1076);
    const auto *sli_1078 = buffer.data(sli + 1078);
    const auto *sli_1079 = buffer.data(sli + 1079);
    const auto *sli_1081 = buffer.data(sli + 1081);
    const auto *sli_1082 = buffer.data(sli + 1082);
    const auto *sli_1084 = buffer.data(sli + 1084);
    const auto *sli_1085 = buffer.data(sli + 1085);
    const auto *sli_1086 = buffer.data(sli + 1086);
    const auto *sli_1087 = buffer.data(sli + 1087);
    const auto *sli_1088 = buffer.data(sli + 1088);
    const auto *sli_1089 = buffer.data(sli + 1089);
    const auto *sli_1090 = buffer.data(sli + 1090);
    const auto *sli_1091 = buffer.data(sli + 1091);
    const auto *sli_1092 = buffer.data(sli + 1092);
    const auto *sli_1094 = buffer.data(sli + 1094);
    const auto *sli_1095 = buffer.data(sli + 1095);
    const auto *sli_1097 = buffer.data(sli + 1097);
    const auto *sli_1098 = buffer.data(sli + 1098);
    const auto *sli_1101 = buffer.data(sli + 1101);
    const auto *sli_1102 = buffer.data(sli + 1102);
    const auto *sli_1104 = buffer.data(sli + 1104);
    const auto *sli_1106 = buffer.data(sli + 1106);
    const auto *sli_1107 = buffer.data(sli + 1107);
    const auto *sli_1109 = buffer.data(sli + 1109);
    const auto *sli_1110 = buffer.data(sli + 1110);
    const auto *sli_1112 = buffer.data(sli + 1112);
    const auto *sli_1113 = buffer.data(sli + 1113);
    const auto *sli_1114 = buffer.data(sli + 1114);
    const auto *sli_1115 = buffer.data(sli + 1115);
    const auto *sli_1116 = buffer.data(sli + 1116);
    const auto *sli_1117 = buffer.data(sli + 1117);
    const auto *sli_1118 = buffer.data(sli + 1118);
    const auto *sli_1119 = buffer.data(sli + 1119);
    const auto *sli_1120 = buffer.data(sli + 1120);
    const auto *sli_1122 = buffer.data(sli + 1122);
    const auto *sli_1123 = buffer.data(sli + 1123);
    const auto *sli_1125 = buffer.data(sli + 1125);
    const auto *sli_1126 = buffer.data(sli + 1126);
    const auto *sli_1129 = buffer.data(sli + 1129);
    const auto *sli_1130 = buffer.data(sli + 1130);
    const auto *sli_1132 = buffer.data(sli + 1132);
    const auto *sli_1134 = buffer.data(sli + 1134);
    const auto *sli_1135 = buffer.data(sli + 1135);
    const auto *sli_1137 = buffer.data(sli + 1137);
    const auto *sli_1138 = buffer.data(sli + 1138);
    const auto *sli_1140 = buffer.data(sli + 1140);
    const auto *sli_1141 = buffer.data(sli + 1141);
    const auto *sli_1142 = buffer.data(sli + 1142);
    const auto *sli_1143 = buffer.data(sli + 1143);
    const auto *sli_1144 = buffer.data(sli + 1144);
    const auto *sli_1145 = buffer.data(sli + 1145);
    const auto *sli_1146 = buffer.data(sli + 1146);
    const auto *sli_1147 = buffer.data(sli + 1147);
    const auto *sli_1148 = buffer.data(sli + 1148);
    const auto *sli_1150 = buffer.data(sli + 1150);
    const auto *sli_1151 = buffer.data(sli + 1151);
    const auto *sli_1153 = buffer.data(sli + 1153);
    const auto *sli_1154 = buffer.data(sli + 1154);
    const auto *sli_1157 = buffer.data(sli + 1157);
    const auto *sli_1158 = buffer.data(sli + 1158);
    const auto *sli_1160 = buffer.data(sli + 1160);
    const auto *sli_1162 = buffer.data(sli + 1162);

#pragma omp simd aligned(t_1371, t_1372, t_1373, pc_x, pc_y, ski_842, slh0_801, slh0_803, \
                         slh1_801, slh1_803, sli_1066, sli_1067, \
                         sli_1069 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1371[k] = f_4 * slh0_801[k]
                    - f_5 * slh1_801[k]
                    + f_3 * pc_x[k] * sli_1067[k];

        t_1372[k] = f_19 * ski_842[k]
                    + f_3 * pc_y[k] * sli_1066[k];

        t_1373[k] = f_4 * slh0_803[k]
                    - f_5 * slh1_803[k]
                    + f_3 * pc_x[k] * sli_1069[k];
    }

#pragma omp simd aligned(t_1374, t_1375, t_1376, pc_x, pc_y, pc_z, ski_815, ski_845, slh0_804, \
                         slh1_804, sli_1067, sli_1069, sli_1070 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1374[k] = f_6 * slh0_804[k]
                    - f_7 * slh1_804[k]
                    + f_3 * pc_x[k] * sli_1070[k];

        t_1375[k] = f_14 * ski_815[k]
                    + f_3 * pc_z[k] * sli_1067[k];

        t_1376[k] = f_19 * ski_845[k]
                    + f_3 * pc_y[k] * sli_1069[k];
    }

#pragma omp simd aligned(t_1377, t_1378, t_1379, pc_x, pc_z, ski_818, slh0_807, slh0_808, \
                         slh1_807, slh1_808, sli_1070, sli_1073, \
                         sli_1074 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1377[k] = f_6 * slh0_807[k]
                    - f_7 * slh1_807[k]
                    + f_3 * pc_x[k] * sli_1073[k];

        t_1378[k] = f_8 * slh0_808[k]
                    - f_9 * slh1_808[k]
                    + f_3 * pc_x[k] * sli_1074[k];

        t_1379[k] = f_14 * ski_818[k]
                    + f_3 * pc_z[k] * sli_1070[k];
    }

#pragma omp simd aligned(t_1380, t_1381, t_1382, pc_x, pc_y, ski_849, slh0_810, slh0_812, \
                         slh1_810, slh1_812, sli_1073, sli_1076, \
                         sli_1078 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1380[k] = f_8 * slh0_810[k]
                    - f_9 * slh1_810[k]
                    + f_3 * pc_x[k] * sli_1076[k];

        t_1381[k] = f_19 * ski_849[k]
                    + f_3 * pc_y[k] * sli_1073[k];

        t_1382[k] = f_8 * slh0_812[k]
                    - f_9 * slh1_812[k]
                    + f_3 * pc_x[k] * sli_1078[k];
    }

#pragma omp simd aligned(t_1383, t_1384, t_1385, pc_x, pc_z, ski_822, slh0_813, slh0_815, \
                         slh1_813, slh1_815, sli_1074, sli_1079, \
                         sli_1081 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1383[k] = f_10 * slh0_813[k]
                    - f_11 * slh1_813[k]
                    + f_3 * pc_x[k] * sli_1079[k];

        t_1384[k] = f_14 * ski_822[k]
                    + f_3 * pc_z[k] * sli_1074[k];

        t_1385[k] = f_10 * slh0_815[k]
                    - f_11 * slh1_815[k]
                    + f_3 * pc_x[k] * sli_1081[k];
    }

#pragma omp simd aligned(t_1386, t_1387, t_1388, t_1389, pc_x, pc_y, ski_854, slh0_816, \
                         slh0_818, slh1_816, slh1_818, sli_1078, sli_1082, sli_1084, \
                         sli_1085 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1386[k] = f_10 * slh0_816[k]
                    - f_11 * slh1_816[k]
                    + f_3 * pc_x[k] * sli_1082[k];

        t_1387[k] = f_19 * ski_854[k]
                    + f_3 * pc_y[k] * sli_1078[k];

        t_1388[k] = f_10 * slh0_818[k]
                    - f_11 * slh1_818[k]
                    + f_3 * pc_x[k] * sli_1084[k];

        t_1389[k] = f_3 * pc_x[k] * sli_1085[k];
    }

#pragma omp simd aligned(t_1390, t_1391, t_1392, t_1393, t_1394, t_1395, pc_x, sli_1086, \
                         sli_1087, sli_1088, sli_1089, sli_1090, \
                         sli_1091 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1390[k] = f_3 * pc_x[k] * sli_1086[k];

        t_1391[k] = f_3 * pc_x[k] * sli_1087[k];

        t_1392[k] = f_3 * pc_x[k] * sli_1088[k];

        t_1393[k] = f_3 * pc_x[k] * sli_1089[k];

        t_1394[k] = f_3 * pc_x[k] * sli_1090[k];

        t_1395[k] = f_3 * pc_x[k] * sli_1091[k];
    }

#pragma omp simd aligned(t_1396, t_1397, t_1398, pc_y, pc_z, ski_833, ski_861, ski_863, \
                         slh0_813, slh0_815, slh1_813, slh1_815, sli_1085, \
                         sli_1087 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1396[k] = f_19 * ski_861[k]
                    + f_1 * slh0_813[k]
                    - f_2 * slh1_813[k]
                    + f_3 * pc_y[k] * sli_1085[k];

        t_1397[k] = f_14 * ski_833[k]
                    + f_3 * pc_z[k] * sli_1085[k];

        t_1398[k] = f_19 * ski_863[k]
                    + f_4 * slh0_815[k]
                    - f_5 * slh1_815[k]
                    + f_3 * pc_y[k] * sli_1087[k];
    }

#pragma omp simd aligned(t_1399, t_1400, t_1401, pc_y, ski_864, ski_865, ski_866, slh0_816, \
                         slh0_817, slh0_818, slh1_816, slh1_817, slh1_818, sli_1088, sli_1089, \
                         sli_1090 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1399[k] = f_19 * ski_864[k]
                    + f_6 * slh0_816[k]
                    - f_7 * slh1_816[k]
                    + f_3 * pc_y[k] * sli_1088[k];

        t_1400[k] = f_19 * ski_865[k]
                    + f_8 * slh0_817[k]
                    - f_9 * slh1_817[k]
                    + f_3 * pc_y[k] * sli_1089[k];

        t_1401[k] = f_19 * ski_866[k]
                    + f_10 * slh0_818[k]
                    - f_11 * slh1_818[k]
                    + f_3 * pc_y[k] * sli_1090[k];
    }

#pragma omp simd aligned(t_1402, t_1403, t_1404, t_1405, pc_x, pc_y, pc_z, ski_839, ski_867, \
                         ski_868, slh0_818, slh0_819, slh1_818, slh1_819, sli_1091, \
                         sli_1092 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1402[k] = f_19 * ski_867[k]
                    + f_3 * pc_y[k] * sli_1091[k];

        t_1403[k] = f_14 * ski_839[k]
                    + f_1 * slh0_818[k]
                    - f_2 * slh1_818[k]
                    + f_3 * pc_z[k] * sli_1091[k];

        t_1404[k] = f_1 * slh0_819[k]
                    - f_2 * slh1_819[k]
                    + f_3 * pc_x[k] * sli_1092[k];

        t_1405[k] = f_17 * ski_868[k]
                    + f_3 * pc_y[k] * sli_1092[k];
    }

#pragma omp simd aligned(t_1406, t_1407, t_1408, pc_x, pc_y, pc_z, ski_840, ski_870, slh0_822, \
                         slh1_822, sli_1092, sli_1094, sli_1095 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1406[k] = f_15 * ski_840[k]
                    + f_3 * pc_z[k] * sli_1092[k];

        t_1407[k] = f_4 * slh0_822[k]
                    - f_5 * slh1_822[k]
                    + f_3 * pc_x[k] * sli_1095[k];

        t_1408[k] = f_17 * ski_870[k]
                    + f_3 * pc_y[k] * sli_1094[k];
    }

#pragma omp simd aligned(t_1409, t_1410, t_1411, t_1412, pc_x, pc_y, pc_z, ski_843, ski_873, \
                         slh0_824, slh0_825, slh1_824, slh1_825, sli_1095, sli_1097, \
                         sli_1098 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1409[k] = f_4 * slh0_824[k]
                    - f_5 * slh1_824[k]
                    + f_3 * pc_x[k] * sli_1097[k];

        t_1410[k] = f_6 * slh0_825[k]
                    - f_7 * slh1_825[k]
                    + f_3 * pc_x[k] * sli_1098[k];

        t_1411[k] = f_15 * ski_843[k]
                    + f_3 * pc_z[k] * sli_1095[k];

        t_1412[k] = f_17 * ski_873[k]
                    + f_3 * pc_y[k] * sli_1097[k];
    }

#pragma omp simd aligned(t_1413, t_1414, t_1415, pc_x, pc_z, ski_846, slh0_828, slh0_829, \
                         slh1_828, slh1_829, sli_1098, sli_1101, \
                         sli_1102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1413[k] = f_6 * slh0_828[k]
                    - f_7 * slh1_828[k]
                    + f_3 * pc_x[k] * sli_1101[k];

        t_1414[k] = f_8 * slh0_829[k]
                    - f_9 * slh1_829[k]
                    + f_3 * pc_x[k] * sli_1102[k];

        t_1415[k] = f_15 * ski_846[k]
                    + f_3 * pc_z[k] * sli_1098[k];
    }

#pragma omp simd aligned(t_1416, t_1417, t_1418, pc_x, pc_y, ski_877, slh0_831, slh0_833, \
                         slh1_831, slh1_833, sli_1101, sli_1104, \
                         sli_1106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1416[k] = f_8 * slh0_831[k]
                    - f_9 * slh1_831[k]
                    + f_3 * pc_x[k] * sli_1104[k];

        t_1417[k] = f_17 * ski_877[k]
                    + f_3 * pc_y[k] * sli_1101[k];

        t_1418[k] = f_8 * slh0_833[k]
                    - f_9 * slh1_833[k]
                    + f_3 * pc_x[k] * sli_1106[k];
    }

#pragma omp simd aligned(t_1419, t_1420, t_1421, pc_x, pc_z, ski_850, slh0_834, slh0_836, \
                         slh1_834, slh1_836, sli_1102, sli_1107, \
                         sli_1109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1419[k] = f_10 * slh0_834[k]
                    - f_11 * slh1_834[k]
                    + f_3 * pc_x[k] * sli_1107[k];

        t_1420[k] = f_15 * ski_850[k]
                    + f_3 * pc_z[k] * sli_1102[k];

        t_1421[k] = f_10 * slh0_836[k]
                    - f_11 * slh1_836[k]
                    + f_3 * pc_x[k] * sli_1109[k];
    }

#pragma omp simd aligned(t_1422, t_1423, t_1424, t_1425, pc_x, pc_y, ski_882, slh0_837, \
                         slh0_839, slh1_837, slh1_839, sli_1106, sli_1110, sli_1112, \
                         sli_1113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1422[k] = f_10 * slh0_837[k]
                    - f_11 * slh1_837[k]
                    + f_3 * pc_x[k] * sli_1110[k];

        t_1423[k] = f_17 * ski_882[k]
                    + f_3 * pc_y[k] * sli_1106[k];

        t_1424[k] = f_10 * slh0_839[k]
                    - f_11 * slh1_839[k]
                    + f_3 * pc_x[k] * sli_1112[k];

        t_1425[k] = f_3 * pc_x[k] * sli_1113[k];
    }

#pragma omp simd aligned(t_1426, t_1427, t_1428, t_1429, t_1430, t_1431, pc_x, sli_1114, \
                         sli_1115, sli_1116, sli_1117, sli_1118, \
                         sli_1119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1426[k] = f_3 * pc_x[k] * sli_1114[k];

        t_1427[k] = f_3 * pc_x[k] * sli_1115[k];

        t_1428[k] = f_3 * pc_x[k] * sli_1116[k];

        t_1429[k] = f_3 * pc_x[k] * sli_1117[k];

        t_1430[k] = f_3 * pc_x[k] * sli_1118[k];

        t_1431[k] = f_3 * pc_x[k] * sli_1119[k];
    }

#pragma omp simd aligned(t_1432, t_1433, t_1434, pc_y, pc_z, ski_861, ski_889, ski_891, \
                         slh0_834, slh0_836, slh1_834, slh1_836, sli_1113, \
                         sli_1115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1432[k] = f_17 * ski_889[k]
                    + f_1 * slh0_834[k]
                    - f_2 * slh1_834[k]
                    + f_3 * pc_y[k] * sli_1113[k];

        t_1433[k] = f_15 * ski_861[k]
                    + f_3 * pc_z[k] * sli_1113[k];

        t_1434[k] = f_17 * ski_891[k]
                    + f_4 * slh0_836[k]
                    - f_5 * slh1_836[k]
                    + f_3 * pc_y[k] * sli_1115[k];
    }

#pragma omp simd aligned(t_1435, t_1436, t_1437, pc_y, ski_892, ski_893, ski_894, slh0_837, \
                         slh0_838, slh0_839, slh1_837, slh1_838, slh1_839, sli_1116, sli_1117, \
                         sli_1118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1435[k] = f_17 * ski_892[k]
                    + f_6 * slh0_837[k]
                    - f_7 * slh1_837[k]
                    + f_3 * pc_y[k] * sli_1116[k];

        t_1436[k] = f_17 * ski_893[k]
                    + f_8 * slh0_838[k]
                    - f_9 * slh1_838[k]
                    + f_3 * pc_y[k] * sli_1117[k];

        t_1437[k] = f_17 * ski_894[k]
                    + f_10 * slh0_839[k]
                    - f_11 * slh1_839[k]
                    + f_3 * pc_y[k] * sli_1118[k];
    }

#pragma omp simd aligned(t_1438, t_1439, t_1440, t_1441, pc_x, pc_y, pc_z, ski_867, ski_895, \
                         ski_896, slh0_839, slh0_840, slh1_839, slh1_840, sli_1119, \
                         sli_1120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1438[k] = f_17 * ski_895[k]
                    + f_3 * pc_y[k] * sli_1119[k];

        t_1439[k] = f_15 * ski_867[k]
                    + f_1 * slh0_839[k]
                    - f_2 * slh1_839[k]
                    + f_3 * pc_z[k] * sli_1119[k];

        t_1440[k] = f_1 * slh0_840[k]
                    - f_2 * slh1_840[k]
                    + f_3 * pc_x[k] * sli_1120[k];

        t_1441[k] = f_16 * ski_896[k]
                    + f_3 * pc_y[k] * sli_1120[k];
    }

#pragma omp simd aligned(t_1442, t_1443, t_1444, pc_x, pc_y, pc_z, ski_868, ski_898, slh0_843, \
                         slh1_843, sli_1120, sli_1122, sli_1123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1442[k] = f_16 * ski_868[k]
                    + f_3 * pc_z[k] * sli_1120[k];

        t_1443[k] = f_4 * slh0_843[k]
                    - f_5 * slh1_843[k]
                    + f_3 * pc_x[k] * sli_1123[k];

        t_1444[k] = f_16 * ski_898[k]
                    + f_3 * pc_y[k] * sli_1122[k];
    }

#pragma omp simd aligned(t_1445, t_1446, t_1447, t_1448, pc_x, pc_y, pc_z, ski_871, ski_901, \
                         slh0_845, slh0_846, slh1_845, slh1_846, sli_1123, sli_1125, \
                         sli_1126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1445[k] = f_4 * slh0_845[k]
                    - f_5 * slh1_845[k]
                    + f_3 * pc_x[k] * sli_1125[k];

        t_1446[k] = f_6 * slh0_846[k]
                    - f_7 * slh1_846[k]
                    + f_3 * pc_x[k] * sli_1126[k];

        t_1447[k] = f_16 * ski_871[k]
                    + f_3 * pc_z[k] * sli_1123[k];

        t_1448[k] = f_16 * ski_901[k]
                    + f_3 * pc_y[k] * sli_1125[k];
    }

#pragma omp simd aligned(t_1449, t_1450, t_1451, pc_x, pc_z, ski_874, slh0_849, slh0_850, \
                         slh1_849, slh1_850, sli_1126, sli_1129, \
                         sli_1130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1449[k] = f_6 * slh0_849[k]
                    - f_7 * slh1_849[k]
                    + f_3 * pc_x[k] * sli_1129[k];

        t_1450[k] = f_8 * slh0_850[k]
                    - f_9 * slh1_850[k]
                    + f_3 * pc_x[k] * sli_1130[k];

        t_1451[k] = f_16 * ski_874[k]
                    + f_3 * pc_z[k] * sli_1126[k];
    }

#pragma omp simd aligned(t_1452, t_1453, t_1454, pc_x, pc_y, ski_905, slh0_852, slh0_854, \
                         slh1_852, slh1_854, sli_1129, sli_1132, \
                         sli_1134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1452[k] = f_8 * slh0_852[k]
                    - f_9 * slh1_852[k]
                    + f_3 * pc_x[k] * sli_1132[k];

        t_1453[k] = f_16 * ski_905[k]
                    + f_3 * pc_y[k] * sli_1129[k];

        t_1454[k] = f_8 * slh0_854[k]
                    - f_9 * slh1_854[k]
                    + f_3 * pc_x[k] * sli_1134[k];
    }

#pragma omp simd aligned(t_1455, t_1456, t_1457, pc_x, pc_z, ski_878, slh0_855, slh0_857, \
                         slh1_855, slh1_857, sli_1130, sli_1135, \
                         sli_1137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1455[k] = f_10 * slh0_855[k]
                    - f_11 * slh1_855[k]
                    + f_3 * pc_x[k] * sli_1135[k];

        t_1456[k] = f_16 * ski_878[k]
                    + f_3 * pc_z[k] * sli_1130[k];

        t_1457[k] = f_10 * slh0_857[k]
                    - f_11 * slh1_857[k]
                    + f_3 * pc_x[k] * sli_1137[k];
    }

#pragma omp simd aligned(t_1458, t_1459, t_1460, t_1461, pc_x, pc_y, ski_910, slh0_858, \
                         slh0_860, slh1_858, slh1_860, sli_1134, sli_1138, sli_1140, \
                         sli_1141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1458[k] = f_10 * slh0_858[k]
                    - f_11 * slh1_858[k]
                    + f_3 * pc_x[k] * sli_1138[k];

        t_1459[k] = f_16 * ski_910[k]
                    + f_3 * pc_y[k] * sli_1134[k];

        t_1460[k] = f_10 * slh0_860[k]
                    - f_11 * slh1_860[k]
                    + f_3 * pc_x[k] * sli_1140[k];

        t_1461[k] = f_3 * pc_x[k] * sli_1141[k];
    }

#pragma omp simd aligned(t_1462, t_1463, t_1464, t_1465, t_1466, t_1467, pc_x, sli_1142, \
                         sli_1143, sli_1144, sli_1145, sli_1146, \
                         sli_1147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1462[k] = f_3 * pc_x[k] * sli_1142[k];

        t_1463[k] = f_3 * pc_x[k] * sli_1143[k];

        t_1464[k] = f_3 * pc_x[k] * sli_1144[k];

        t_1465[k] = f_3 * pc_x[k] * sli_1145[k];

        t_1466[k] = f_3 * pc_x[k] * sli_1146[k];

        t_1467[k] = f_3 * pc_x[k] * sli_1147[k];
    }

#pragma omp simd aligned(t_1468, t_1469, t_1470, pc_y, pc_z, ski_889, ski_917, ski_919, \
                         slh0_855, slh0_857, slh1_855, slh1_857, sli_1141, \
                         sli_1143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1468[k] = f_16 * ski_917[k]
                    + f_1 * slh0_855[k]
                    - f_2 * slh1_855[k]
                    + f_3 * pc_y[k] * sli_1141[k];

        t_1469[k] = f_16 * ski_889[k]
                    + f_3 * pc_z[k] * sli_1141[k];

        t_1470[k] = f_16 * ski_919[k]
                    + f_4 * slh0_857[k]
                    - f_5 * slh1_857[k]
                    + f_3 * pc_y[k] * sli_1143[k];
    }

#pragma omp simd aligned(t_1471, t_1472, t_1473, pc_y, ski_920, ski_921, ski_922, slh0_858, \
                         slh0_859, slh0_860, slh1_858, slh1_859, slh1_860, sli_1144, sli_1145, \
                         sli_1146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1471[k] = f_16 * ski_920[k]
                    + f_6 * slh0_858[k]
                    - f_7 * slh1_858[k]
                    + f_3 * pc_y[k] * sli_1144[k];

        t_1472[k] = f_16 * ski_921[k]
                    + f_8 * slh0_859[k]
                    - f_9 * slh1_859[k]
                    + f_3 * pc_y[k] * sli_1145[k];

        t_1473[k] = f_16 * ski_922[k]
                    + f_10 * slh0_860[k]
                    - f_11 * slh1_860[k]
                    + f_3 * pc_y[k] * sli_1146[k];
    }

#pragma omp simd aligned(t_1474, t_1475, t_1476, t_1477, pc_x, pc_y, pc_z, ski_895, ski_923, \
                         ski_924, slh0_860, slh0_861, slh1_860, slh1_861, sli_1147, \
                         sli_1148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1474[k] = f_16 * ski_923[k]
                    + f_3 * pc_y[k] * sli_1147[k];

        t_1475[k] = f_16 * ski_895[k]
                    + f_1 * slh0_860[k]
                    - f_2 * slh1_860[k]
                    + f_3 * pc_z[k] * sli_1147[k];

        t_1476[k] = f_1 * slh0_861[k]
                    - f_2 * slh1_861[k]
                    + f_3 * pc_x[k] * sli_1148[k];

        t_1477[k] = f_15 * ski_924[k]
                    + f_3 * pc_y[k] * sli_1148[k];
    }

#pragma omp simd aligned(t_1478, t_1479, t_1480, pc_x, pc_y, pc_z, ski_896, ski_926, slh0_864, \
                         slh1_864, sli_1148, sli_1150, sli_1151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1478[k] = f_17 * ski_896[k]
                    + f_3 * pc_z[k] * sli_1148[k];

        t_1479[k] = f_4 * slh0_864[k]
                    - f_5 * slh1_864[k]
                    + f_3 * pc_x[k] * sli_1151[k];

        t_1480[k] = f_15 * ski_926[k]
                    + f_3 * pc_y[k] * sli_1150[k];
    }

#pragma omp simd aligned(t_1481, t_1482, t_1483, t_1484, pc_x, pc_y, pc_z, ski_899, ski_929, \
                         slh0_866, slh0_867, slh1_866, slh1_867, sli_1151, sli_1153, \
                         sli_1154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1481[k] = f_4 * slh0_866[k]
                    - f_5 * slh1_866[k]
                    + f_3 * pc_x[k] * sli_1153[k];

        t_1482[k] = f_6 * slh0_867[k]
                    - f_7 * slh1_867[k]
                    + f_3 * pc_x[k] * sli_1154[k];

        t_1483[k] = f_17 * ski_899[k]
                    + f_3 * pc_z[k] * sli_1151[k];

        t_1484[k] = f_15 * ski_929[k]
                    + f_3 * pc_y[k] * sli_1153[k];
    }

#pragma omp simd aligned(t_1485, t_1486, t_1487, pc_x, pc_z, ski_902, slh0_870, slh0_871, \
                         slh1_870, slh1_871, sli_1154, sli_1157, \
                         sli_1158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1485[k] = f_6 * slh0_870[k]
                    - f_7 * slh1_870[k]
                    + f_3 * pc_x[k] * sli_1157[k];

        t_1486[k] = f_8 * slh0_871[k]
                    - f_9 * slh1_871[k]
                    + f_3 * pc_x[k] * sli_1158[k];

        t_1487[k] = f_17 * ski_902[k]
                    + f_3 * pc_z[k] * sli_1154[k];
    }

#pragma omp simd aligned(t_1488, t_1489, t_1490, pc_x, pc_y, ski_933, slh0_873, slh0_875, \
                         slh1_873, slh1_875, sli_1157, sli_1160, \
                         sli_1162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1488[k] = f_8 * slh0_873[k]
                    - f_9 * slh1_873[k]
                    + f_3 * pc_x[k] * sli_1160[k];

        t_1489[k] = f_15 * ski_933[k]
                    + f_3 * pc_y[k] * sli_1157[k];

        t_1490[k] = f_8 * slh0_875[k]
                    - f_9 * slh1_875[k]
                    + f_3 * pc_x[k] * sli_1162[k];
    }
}

static auto
compute_prim_slk_three_center_electron_repulsion_0_piece13(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t skk0,
                                                           const size_t ski, const size_t skk1,
                                                           const size_t slh0, const size_t slh1,
                                                           const size_t sli, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
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
    const auto f_19 = 3.0 / q;

    auto *t_1491 = buffer.data(target + 1491);
    auto *t_1492 = buffer.data(target + 1492);
    auto *t_1493 = buffer.data(target + 1493);
    auto *t_1494 = buffer.data(target + 1494);
    auto *t_1495 = buffer.data(target + 1495);
    auto *t_1496 = buffer.data(target + 1496);
    auto *t_1497 = buffer.data(target + 1497);
    auto *t_1498 = buffer.data(target + 1498);
    auto *t_1499 = buffer.data(target + 1499);
    auto *t_1500 = buffer.data(target + 1500);
    auto *t_1501 = buffer.data(target + 1501);
    auto *t_1502 = buffer.data(target + 1502);
    auto *t_1503 = buffer.data(target + 1503);
    auto *t_1504 = buffer.data(target + 1504);
    auto *t_1505 = buffer.data(target + 1505);
    auto *t_1506 = buffer.data(target + 1506);
    auto *t_1507 = buffer.data(target + 1507);
    auto *t_1508 = buffer.data(target + 1508);
    auto *t_1509 = buffer.data(target + 1509);
    auto *t_1510 = buffer.data(target + 1510);
    auto *t_1511 = buffer.data(target + 1511);
    auto *t_1512 = buffer.data(target + 1512);
    auto *t_1513 = buffer.data(target + 1513);
    auto *t_1514 = buffer.data(target + 1514);
    auto *t_1515 = buffer.data(target + 1515);
    auto *t_1516 = buffer.data(target + 1516);
    auto *t_1517 = buffer.data(target + 1517);
    auto *t_1518 = buffer.data(target + 1518);
    auto *t_1519 = buffer.data(target + 1519);
    auto *t_1520 = buffer.data(target + 1520);
    auto *t_1521 = buffer.data(target + 1521);
    auto *t_1522 = buffer.data(target + 1522);
    auto *t_1523 = buffer.data(target + 1523);
    auto *t_1524 = buffer.data(target + 1524);
    auto *t_1525 = buffer.data(target + 1525);
    auto *t_1526 = buffer.data(target + 1526);
    auto *t_1527 = buffer.data(target + 1527);
    auto *t_1528 = buffer.data(target + 1528);
    auto *t_1529 = buffer.data(target + 1529);
    auto *t_1530 = buffer.data(target + 1530);
    auto *t_1531 = buffer.data(target + 1531);
    auto *t_1532 = buffer.data(target + 1532);
    auto *t_1533 = buffer.data(target + 1533);
    auto *t_1534 = buffer.data(target + 1534);
    auto *t_1535 = buffer.data(target + 1535);
    auto *t_1536 = buffer.data(target + 1536);
    auto *t_1537 = buffer.data(target + 1537);
    auto *t_1538 = buffer.data(target + 1538);
    auto *t_1539 = buffer.data(target + 1539);
    auto *t_1540 = buffer.data(target + 1540);
    auto *t_1541 = buffer.data(target + 1541);
    auto *t_1542 = buffer.data(target + 1542);
    auto *t_1543 = buffer.data(target + 1543);
    auto *t_1544 = buffer.data(target + 1544);
    auto *t_1545 = buffer.data(target + 1545);
    auto *t_1546 = buffer.data(target + 1546);
    auto *t_1547 = buffer.data(target + 1547);
    auto *t_1548 = buffer.data(target + 1548);
    auto *t_1549 = buffer.data(target + 1549);
    auto *t_1550 = buffer.data(target + 1550);
    auto *t_1551 = buffer.data(target + 1551);
    auto *t_1552 = buffer.data(target + 1552);
    auto *t_1553 = buffer.data(target + 1553);
    auto *t_1554 = buffer.data(target + 1554);
    auto *t_1555 = buffer.data(target + 1555);
    auto *t_1556 = buffer.data(target + 1556);
    auto *t_1557 = buffer.data(target + 1557);
    auto *t_1558 = buffer.data(target + 1558);
    auto *t_1559 = buffer.data(target + 1559);
    auto *t_1560 = buffer.data(target + 1560);
    auto *t_1561 = buffer.data(target + 1561);
    auto *t_1562 = buffer.data(target + 1562);
    auto *t_1563 = buffer.data(target + 1563);
    auto *t_1564 = buffer.data(target + 1564);
    auto *t_1565 = buffer.data(target + 1565);
    auto *t_1566 = buffer.data(target + 1566);
    auto *t_1567 = buffer.data(target + 1567);
    auto *t_1568 = buffer.data(target + 1568);
    auto *t_1569 = buffer.data(target + 1569);
    auto *t_1570 = buffer.data(target + 1570);
    auto *t_1571 = buffer.data(target + 1571);
    auto *t_1572 = buffer.data(target + 1572);
    auto *t_1573 = buffer.data(target + 1573);
    auto *t_1574 = buffer.data(target + 1574);
    auto *t_1575 = buffer.data(target + 1575);
    auto *t_1576 = buffer.data(target + 1576);
    auto *t_1577 = buffer.data(target + 1577);
    auto *t_1578 = buffer.data(target + 1578);
    auto *t_1579 = buffer.data(target + 1579);
    auto *t_1580 = buffer.data(target + 1580);
    auto *t_1581 = buffer.data(target + 1581);
    auto *t_1582 = buffer.data(target + 1582);
    auto *t_1583 = buffer.data(target + 1583);
    auto *t_1584 = buffer.data(target + 1584);
    auto *t_1585 = buffer.data(target + 1585);
    auto *t_1586 = buffer.data(target + 1586);
    auto *t_1587 = buffer.data(target + 1587);
    auto *t_1588 = buffer.data(target + 1588);
    auto *t_1589 = buffer.data(target + 1589);
    auto *t_1590 = buffer.data(target + 1590);
    auto *t_1591 = buffer.data(target + 1591);
    auto *t_1592 = buffer.data(target + 1592);
    auto *t_1593 = buffer.data(target + 1593);
    auto *t_1594 = buffer.data(target + 1594);
    auto *t_1595 = buffer.data(target + 1595);
    auto *t_1596 = buffer.data(target + 1596);
    auto *t_1597 = buffer.data(target + 1597);
    auto *t_1598 = buffer.data(target + 1598);
    auto *t_1599 = buffer.data(target + 1599);
    auto *t_1600 = buffer.data(target + 1600);
    auto *t_1601 = buffer.data(target + 1601);
    auto *t_1602 = buffer.data(target + 1602);
    auto *t_1603 = buffer.data(target + 1603);
    auto *t_1604 = buffer.data(target + 1604);
    auto *t_1605 = buffer.data(target + 1605);
    auto *t_1606 = buffer.data(target + 1606);
    auto *t_1607 = buffer.data(target + 1607);
    auto *t_1608 = buffer.data(target + 1608);
    auto *t_1609 = buffer.data(target + 1609);
    auto *t_1610 = buffer.data(target + 1610);
    auto *t_1611 = buffer.data(target + 1611);
    auto *t_1612 = buffer.data(target + 1612);
    auto *t_1613 = buffer.data(target + 1613);
    auto *t_1614 = buffer.data(target + 1614);
    auto *t_1615 = buffer.data(target + 1615);
    auto *t_1616 = buffer.data(target + 1616);

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *skk0_1260 = buffer.data(skk0 + 1260);
    const auto *skk0_1265 = buffer.data(skk0 + 1265);
    const auto *skk0_1269 = buffer.data(skk0 + 1269);
    const auto *skk0_1274 = buffer.data(skk0 + 1274);
    const auto *skk0_1280 = buffer.data(skk0 + 1280);
    const auto *skk0_1288 = buffer.data(skk0 + 1288);
    const auto *skk0_1290 = buffer.data(skk0 + 1290);
    const auto *skk0_1291 = buffer.data(skk0 + 1291);
    const auto *skk0_1292 = buffer.data(skk0 + 1292);
    const auto *skk0_1293 = buffer.data(skk0 + 1293);
    const auto *skk0_1295 = buffer.data(skk0 + 1295);

    const auto *ski_906 = buffer.data(ski + 906);
    const auto *ski_917 = buffer.data(ski + 917);
    const auto *ski_923 = buffer.data(ski + 923);
    const auto *ski_924 = buffer.data(ski + 924);
    const auto *ski_927 = buffer.data(ski + 927);
    const auto *ski_930 = buffer.data(ski + 930);
    const auto *ski_934 = buffer.data(ski + 934);
    const auto *ski_938 = buffer.data(ski + 938);
    const auto *ski_945 = buffer.data(ski + 945);
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
    const auto *ski_966 = buffer.data(ski + 966);
    const auto *ski_973 = buffer.data(ski + 973);
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
    const auto *ski_994 = buffer.data(ski + 994);
    const auto *ski_1001 = buffer.data(ski + 1001);
    const auto *ski_1003 = buffer.data(ski + 1003);
    const auto *ski_1004 = buffer.data(ski + 1004);
    const auto *ski_1005 = buffer.data(ski + 1005);
    const auto *ski_1006 = buffer.data(ski + 1006);
    const auto *ski_1007 = buffer.data(ski + 1007);

    const auto *skk1_1260 = buffer.data(skk1 + 1260);
    const auto *skk1_1265 = buffer.data(skk1 + 1265);
    const auto *skk1_1269 = buffer.data(skk1 + 1269);
    const auto *skk1_1274 = buffer.data(skk1 + 1274);
    const auto *skk1_1280 = buffer.data(skk1 + 1280);
    const auto *skk1_1288 = buffer.data(skk1 + 1288);
    const auto *skk1_1290 = buffer.data(skk1 + 1290);
    const auto *skk1_1291 = buffer.data(skk1 + 1291);
    const auto *skk1_1292 = buffer.data(skk1 + 1292);
    const auto *skk1_1293 = buffer.data(skk1 + 1293);
    const auto *skk1_1295 = buffer.data(skk1 + 1295);

    const auto *slh0_876 = buffer.data(slh0 + 876);
    const auto *slh0_878 = buffer.data(slh0 + 878);
    const auto *slh0_879 = buffer.data(slh0 + 879);
    const auto *slh0_880 = buffer.data(slh0 + 880);
    const auto *slh0_881 = buffer.data(slh0 + 881);
    const auto *slh0_882 = buffer.data(slh0 + 882);
    const auto *slh0_885 = buffer.data(slh0 + 885);
    const auto *slh0_887 = buffer.data(slh0 + 887);
    const auto *slh0_888 = buffer.data(slh0 + 888);
    const auto *slh0_891 = buffer.data(slh0 + 891);
    const auto *slh0_892 = buffer.data(slh0 + 892);
    const auto *slh0_894 = buffer.data(slh0 + 894);
    const auto *slh0_896 = buffer.data(slh0 + 896);
    const auto *slh0_897 = buffer.data(slh0 + 897);
    const auto *slh0_899 = buffer.data(slh0 + 899);
    const auto *slh0_900 = buffer.data(slh0 + 900);
    const auto *slh0_901 = buffer.data(slh0 + 901);
    const auto *slh0_902 = buffer.data(slh0 + 902);
    const auto *slh0_906 = buffer.data(slh0 + 906);
    const auto *slh0_909 = buffer.data(slh0 + 909);
    const auto *slh0_913 = buffer.data(slh0 + 913);
    const auto *slh0_915 = buffer.data(slh0 + 915);
    const auto *slh0_918 = buffer.data(slh0 + 918);
    const auto *slh0_920 = buffer.data(slh0 + 920);
    const auto *slh0_921 = buffer.data(slh0 + 921);
    const auto *slh0_924 = buffer.data(slh0 + 924);
    const auto *slh0_927 = buffer.data(slh0 + 927);
    const auto *slh0_929 = buffer.data(slh0 + 929);
    const auto *slh0_930 = buffer.data(slh0 + 930);
    const auto *slh0_933 = buffer.data(slh0 + 933);
    const auto *slh0_934 = buffer.data(slh0 + 934);
    const auto *slh0_936 = buffer.data(slh0 + 936);
    const auto *slh0_938 = buffer.data(slh0 + 938);
    const auto *slh0_939 = buffer.data(slh0 + 939);
    const auto *slh0_941 = buffer.data(slh0 + 941);
    const auto *slh0_942 = buffer.data(slh0 + 942);
    const auto *slh0_943 = buffer.data(slh0 + 943);
    const auto *slh0_944 = buffer.data(slh0 + 944);

    const auto *slh1_876 = buffer.data(slh1 + 876);
    const auto *slh1_878 = buffer.data(slh1 + 878);
    const auto *slh1_879 = buffer.data(slh1 + 879);
    const auto *slh1_880 = buffer.data(slh1 + 880);
    const auto *slh1_881 = buffer.data(slh1 + 881);
    const auto *slh1_882 = buffer.data(slh1 + 882);
    const auto *slh1_885 = buffer.data(slh1 + 885);
    const auto *slh1_887 = buffer.data(slh1 + 887);
    const auto *slh1_888 = buffer.data(slh1 + 888);
    const auto *slh1_891 = buffer.data(slh1 + 891);
    const auto *slh1_892 = buffer.data(slh1 + 892);
    const auto *slh1_894 = buffer.data(slh1 + 894);
    const auto *slh1_896 = buffer.data(slh1 + 896);
    const auto *slh1_897 = buffer.data(slh1 + 897);
    const auto *slh1_899 = buffer.data(slh1 + 899);
    const auto *slh1_900 = buffer.data(slh1 + 900);
    const auto *slh1_901 = buffer.data(slh1 + 901);
    const auto *slh1_902 = buffer.data(slh1 + 902);
    const auto *slh1_906 = buffer.data(slh1 + 906);
    const auto *slh1_909 = buffer.data(slh1 + 909);
    const auto *slh1_913 = buffer.data(slh1 + 913);
    const auto *slh1_915 = buffer.data(slh1 + 915);
    const auto *slh1_918 = buffer.data(slh1 + 918);
    const auto *slh1_920 = buffer.data(slh1 + 920);
    const auto *slh1_921 = buffer.data(slh1 + 921);
    const auto *slh1_924 = buffer.data(slh1 + 924);
    const auto *slh1_927 = buffer.data(slh1 + 927);
    const auto *slh1_929 = buffer.data(slh1 + 929);
    const auto *slh1_930 = buffer.data(slh1 + 930);
    const auto *slh1_933 = buffer.data(slh1 + 933);
    const auto *slh1_934 = buffer.data(slh1 + 934);
    const auto *slh1_936 = buffer.data(slh1 + 936);
    const auto *slh1_938 = buffer.data(slh1 + 938);
    const auto *slh1_939 = buffer.data(slh1 + 939);
    const auto *slh1_941 = buffer.data(slh1 + 941);
    const auto *slh1_942 = buffer.data(slh1 + 942);
    const auto *slh1_943 = buffer.data(slh1 + 943);
    const auto *slh1_944 = buffer.data(slh1 + 944);

    const auto *sli_1158 = buffer.data(sli + 1158);
    const auto *sli_1162 = buffer.data(sli + 1162);
    const auto *sli_1163 = buffer.data(sli + 1163);
    const auto *sli_1165 = buffer.data(sli + 1165);
    const auto *sli_1166 = buffer.data(sli + 1166);
    const auto *sli_1168 = buffer.data(sli + 1168);
    const auto *sli_1169 = buffer.data(sli + 1169);
    const auto *sli_1170 = buffer.data(sli + 1170);
    const auto *sli_1171 = buffer.data(sli + 1171);
    const auto *sli_1172 = buffer.data(sli + 1172);
    const auto *sli_1173 = buffer.data(sli + 1173);
    const auto *sli_1174 = buffer.data(sli + 1174);
    const auto *sli_1175 = buffer.data(sli + 1175);
    const auto *sli_1176 = buffer.data(sli + 1176);
    const auto *sli_1178 = buffer.data(sli + 1178);
    const auto *sli_1179 = buffer.data(sli + 1179);
    const auto *sli_1181 = buffer.data(sli + 1181);
    const auto *sli_1182 = buffer.data(sli + 1182);
    const auto *sli_1185 = buffer.data(sli + 1185);
    const auto *sli_1186 = buffer.data(sli + 1186);
    const auto *sli_1188 = buffer.data(sli + 1188);
    const auto *sli_1190 = buffer.data(sli + 1190);
    const auto *sli_1191 = buffer.data(sli + 1191);
    const auto *sli_1193 = buffer.data(sli + 1193);
    const auto *sli_1194 = buffer.data(sli + 1194);
    const auto *sli_1196 = buffer.data(sli + 1196);
    const auto *sli_1197 = buffer.data(sli + 1197);
    const auto *sli_1198 = buffer.data(sli + 1198);
    const auto *sli_1199 = buffer.data(sli + 1199);
    const auto *sli_1200 = buffer.data(sli + 1200);
    const auto *sli_1201 = buffer.data(sli + 1201);
    const auto *sli_1202 = buffer.data(sli + 1202);
    const auto *sli_1203 = buffer.data(sli + 1203);
    const auto *sli_1204 = buffer.data(sli + 1204);
    const auto *sli_1206 = buffer.data(sli + 1206);
    const auto *sli_1207 = buffer.data(sli + 1207);
    const auto *sli_1209 = buffer.data(sli + 1209);
    const auto *sli_1210 = buffer.data(sli + 1210);
    const auto *sli_1213 = buffer.data(sli + 1213);
    const auto *sli_1214 = buffer.data(sli + 1214);
    const auto *sli_1216 = buffer.data(sli + 1216);
    const auto *sli_1218 = buffer.data(sli + 1218);
    const auto *sli_1219 = buffer.data(sli + 1219);
    const auto *sli_1221 = buffer.data(sli + 1221);
    const auto *sli_1222 = buffer.data(sli + 1222);
    const auto *sli_1225 = buffer.data(sli + 1225);
    const auto *sli_1226 = buffer.data(sli + 1226);
    const auto *sli_1227 = buffer.data(sli + 1227);
    const auto *sli_1228 = buffer.data(sli + 1228);
    const auto *sli_1229 = buffer.data(sli + 1229);
    const auto *sli_1230 = buffer.data(sli + 1230);
    const auto *sli_1231 = buffer.data(sli + 1231);
    const auto *sli_1232 = buffer.data(sli + 1232);
    const auto *sli_1234 = buffer.data(sli + 1234);
    const auto *sli_1235 = buffer.data(sli + 1235);
    const auto *sli_1237 = buffer.data(sli + 1237);
    const auto *sli_1238 = buffer.data(sli + 1238);
    const auto *sli_1241 = buffer.data(sli + 1241);
    const auto *sli_1242 = buffer.data(sli + 1242);
    const auto *sli_1244 = buffer.data(sli + 1244);
    const auto *sli_1246 = buffer.data(sli + 1246);
    const auto *sli_1247 = buffer.data(sli + 1247);
    const auto *sli_1249 = buffer.data(sli + 1249);
    const auto *sli_1250 = buffer.data(sli + 1250);
    const auto *sli_1252 = buffer.data(sli + 1252);
    const auto *sli_1253 = buffer.data(sli + 1253);
    const auto *sli_1254 = buffer.data(sli + 1254);
    const auto *sli_1255 = buffer.data(sli + 1255);
    const auto *sli_1256 = buffer.data(sli + 1256);
    const auto *sli_1257 = buffer.data(sli + 1257);
    const auto *sli_1258 = buffer.data(sli + 1258);
    const auto *sli_1259 = buffer.data(sli + 1259);

#pragma omp simd aligned(t_1491, t_1492, t_1493, pc_x, pc_z, ski_906, slh0_876, slh0_878, \
                         slh1_876, slh1_878, sli_1158, sli_1163, \
                         sli_1165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1491[k] = f_10 * slh0_876[k]
                    - f_11 * slh1_876[k]
                    + f_3 * pc_x[k] * sli_1163[k];

        t_1492[k] = f_17 * ski_906[k]
                    + f_3 * pc_z[k] * sli_1158[k];

        t_1493[k] = f_10 * slh0_878[k]
                    - f_11 * slh1_878[k]
                    + f_3 * pc_x[k] * sli_1165[k];
    }

#pragma omp simd aligned(t_1494, t_1495, t_1496, t_1497, pc_x, pc_y, ski_938, slh0_879, \
                         slh0_881, slh1_879, slh1_881, sli_1162, sli_1166, sli_1168, \
                         sli_1169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1494[k] = f_10 * slh0_879[k]
                    - f_11 * slh1_879[k]
                    + f_3 * pc_x[k] * sli_1166[k];

        t_1495[k] = f_15 * ski_938[k]
                    + f_3 * pc_y[k] * sli_1162[k];

        t_1496[k] = f_10 * slh0_881[k]
                    - f_11 * slh1_881[k]
                    + f_3 * pc_x[k] * sli_1168[k];

        t_1497[k] = f_3 * pc_x[k] * sli_1169[k];
    }

#pragma omp simd aligned(t_1498, t_1499, t_1500, t_1501, t_1502, t_1503, pc_x, sli_1170, \
                         sli_1171, sli_1172, sli_1173, sli_1174, \
                         sli_1175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1498[k] = f_3 * pc_x[k] * sli_1170[k];

        t_1499[k] = f_3 * pc_x[k] * sli_1171[k];

        t_1500[k] = f_3 * pc_x[k] * sli_1172[k];

        t_1501[k] = f_3 * pc_x[k] * sli_1173[k];

        t_1502[k] = f_3 * pc_x[k] * sli_1174[k];

        t_1503[k] = f_3 * pc_x[k] * sli_1175[k];
    }

#pragma omp simd aligned(t_1504, t_1505, t_1506, pc_y, pc_z, ski_917, ski_945, ski_947, \
                         slh0_876, slh0_878, slh1_876, slh1_878, sli_1169, \
                         sli_1171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1504[k] = f_15 * ski_945[k]
                    + f_1 * slh0_876[k]
                    - f_2 * slh1_876[k]
                    + f_3 * pc_y[k] * sli_1169[k];

        t_1505[k] = f_17 * ski_917[k]
                    + f_3 * pc_z[k] * sli_1169[k];

        t_1506[k] = f_15 * ski_947[k]
                    + f_4 * slh0_878[k]
                    - f_5 * slh1_878[k]
                    + f_3 * pc_y[k] * sli_1171[k];
    }

#pragma omp simd aligned(t_1507, t_1508, t_1509, pc_y, ski_948, ski_949, ski_950, slh0_879, \
                         slh0_880, slh0_881, slh1_879, slh1_880, slh1_881, sli_1172, sli_1173, \
                         sli_1174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1507[k] = f_15 * ski_948[k]
                    + f_6 * slh0_879[k]
                    - f_7 * slh1_879[k]
                    + f_3 * pc_y[k] * sli_1172[k];

        t_1508[k] = f_15 * ski_949[k]
                    + f_8 * slh0_880[k]
                    - f_9 * slh1_880[k]
                    + f_3 * pc_y[k] * sli_1173[k];

        t_1509[k] = f_15 * ski_950[k]
                    + f_10 * slh0_881[k]
                    - f_11 * slh1_881[k]
                    + f_3 * pc_y[k] * sli_1174[k];
    }

#pragma omp simd aligned(t_1510, t_1511, t_1512, t_1513, pc_x, pc_y, pc_z, ski_923, ski_951, \
                         ski_952, slh0_881, slh0_882, slh1_881, slh1_882, sli_1175, \
                         sli_1176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1510[k] = f_15 * ski_951[k]
                    + f_3 * pc_y[k] * sli_1175[k];

        t_1511[k] = f_17 * ski_923[k]
                    + f_1 * slh0_881[k]
                    - f_2 * slh1_881[k]
                    + f_3 * pc_z[k] * sli_1175[k];

        t_1512[k] = f_1 * slh0_882[k]
                    - f_2 * slh1_882[k]
                    + f_3 * pc_x[k] * sli_1176[k];

        t_1513[k] = f_14 * ski_952[k]
                    + f_3 * pc_y[k] * sli_1176[k];
    }

#pragma omp simd aligned(t_1514, t_1515, t_1516, pc_x, pc_y, pc_z, ski_924, ski_954, slh0_885, \
                         slh1_885, sli_1176, sli_1178, sli_1179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1514[k] = f_19 * ski_924[k]
                    + f_3 * pc_z[k] * sli_1176[k];

        t_1515[k] = f_4 * slh0_885[k]
                    - f_5 * slh1_885[k]
                    + f_3 * pc_x[k] * sli_1179[k];

        t_1516[k] = f_14 * ski_954[k]
                    + f_3 * pc_y[k] * sli_1178[k];
    }

#pragma omp simd aligned(t_1517, t_1518, t_1519, t_1520, pc_x, pc_y, pc_z, ski_927, ski_957, \
                         slh0_887, slh0_888, slh1_887, slh1_888, sli_1179, sli_1181, \
                         sli_1182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1517[k] = f_4 * slh0_887[k]
                    - f_5 * slh1_887[k]
                    + f_3 * pc_x[k] * sli_1181[k];

        t_1518[k] = f_6 * slh0_888[k]
                    - f_7 * slh1_888[k]
                    + f_3 * pc_x[k] * sli_1182[k];

        t_1519[k] = f_19 * ski_927[k]
                    + f_3 * pc_z[k] * sli_1179[k];

        t_1520[k] = f_14 * ski_957[k]
                    + f_3 * pc_y[k] * sli_1181[k];
    }

#pragma omp simd aligned(t_1521, t_1522, t_1523, pc_x, pc_z, ski_930, slh0_891, slh0_892, \
                         slh1_891, slh1_892, sli_1182, sli_1185, \
                         sli_1186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1521[k] = f_6 * slh0_891[k]
                    - f_7 * slh1_891[k]
                    + f_3 * pc_x[k] * sli_1185[k];

        t_1522[k] = f_8 * slh0_892[k]
                    - f_9 * slh1_892[k]
                    + f_3 * pc_x[k] * sli_1186[k];

        t_1523[k] = f_19 * ski_930[k]
                    + f_3 * pc_z[k] * sli_1182[k];
    }

#pragma omp simd aligned(t_1524, t_1525, t_1526, pc_x, pc_y, ski_961, slh0_894, slh0_896, \
                         slh1_894, slh1_896, sli_1185, sli_1188, \
                         sli_1190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1524[k] = f_8 * slh0_894[k]
                    - f_9 * slh1_894[k]
                    + f_3 * pc_x[k] * sli_1188[k];

        t_1525[k] = f_14 * ski_961[k]
                    + f_3 * pc_y[k] * sli_1185[k];

        t_1526[k] = f_8 * slh0_896[k]
                    - f_9 * slh1_896[k]
                    + f_3 * pc_x[k] * sli_1190[k];
    }

#pragma omp simd aligned(t_1527, t_1528, t_1529, pc_x, pc_z, ski_934, slh0_897, slh0_899, \
                         slh1_897, slh1_899, sli_1186, sli_1191, \
                         sli_1193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1527[k] = f_10 * slh0_897[k]
                    - f_11 * slh1_897[k]
                    + f_3 * pc_x[k] * sli_1191[k];

        t_1528[k] = f_19 * ski_934[k]
                    + f_3 * pc_z[k] * sli_1186[k];

        t_1529[k] = f_10 * slh0_899[k]
                    - f_11 * slh1_899[k]
                    + f_3 * pc_x[k] * sli_1193[k];
    }

#pragma omp simd aligned(t_1530, t_1531, t_1532, t_1533, pc_x, pc_y, ski_966, slh0_900, \
                         slh0_902, slh1_900, slh1_902, sli_1190, sli_1194, sli_1196, \
                         sli_1197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1530[k] = f_10 * slh0_900[k]
                    - f_11 * slh1_900[k]
                    + f_3 * pc_x[k] * sli_1194[k];

        t_1531[k] = f_14 * ski_966[k]
                    + f_3 * pc_y[k] * sli_1190[k];

        t_1532[k] = f_10 * slh0_902[k]
                    - f_11 * slh1_902[k]
                    + f_3 * pc_x[k] * sli_1196[k];

        t_1533[k] = f_3 * pc_x[k] * sli_1197[k];
    }

#pragma omp simd aligned(t_1534, t_1535, t_1536, t_1537, t_1538, t_1539, pc_x, sli_1198, \
                         sli_1199, sli_1200, sli_1201, sli_1202, \
                         sli_1203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1534[k] = f_3 * pc_x[k] * sli_1198[k];

        t_1535[k] = f_3 * pc_x[k] * sli_1199[k];

        t_1536[k] = f_3 * pc_x[k] * sli_1200[k];

        t_1537[k] = f_3 * pc_x[k] * sli_1201[k];

        t_1538[k] = f_3 * pc_x[k] * sli_1202[k];

        t_1539[k] = f_3 * pc_x[k] * sli_1203[k];
    }

#pragma omp simd aligned(t_1540, t_1541, t_1542, pc_y, pc_z, ski_945, ski_973, ski_975, \
                         slh0_897, slh0_899, slh1_897, slh1_899, sli_1197, \
                         sli_1199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1540[k] = f_14 * ski_973[k]
                    + f_1 * slh0_897[k]
                    - f_2 * slh1_897[k]
                    + f_3 * pc_y[k] * sli_1197[k];

        t_1541[k] = f_19 * ski_945[k]
                    + f_3 * pc_z[k] * sli_1197[k];

        t_1542[k] = f_14 * ski_975[k]
                    + f_4 * slh0_899[k]
                    - f_5 * slh1_899[k]
                    + f_3 * pc_y[k] * sli_1199[k];
    }

#pragma omp simd aligned(t_1543, t_1544, t_1545, pc_y, ski_976, ski_977, ski_978, slh0_900, \
                         slh0_901, slh0_902, slh1_900, slh1_901, slh1_902, sli_1200, sli_1201, \
                         sli_1202 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1543[k] = f_14 * ski_976[k]
                    + f_6 * slh0_900[k]
                    - f_7 * slh1_900[k]
                    + f_3 * pc_y[k] * sli_1200[k];

        t_1544[k] = f_14 * ski_977[k]
                    + f_8 * slh0_901[k]
                    - f_9 * slh1_901[k]
                    + f_3 * pc_y[k] * sli_1201[k];

        t_1545[k] = f_14 * ski_978[k]
                    + f_10 * slh0_902[k]
                    - f_11 * slh1_902[k]
                    + f_3 * pc_y[k] * sli_1202[k];
    }

#pragma omp simd aligned(t_1546, t_1547, t_1548, t_1549, pb_y, pc_y, pc_z, skk0_1260, ski_951, \
                         ski_979, ski_980, skk1_1260, slh0_902, slh1_902, sli_1203, \
                         sli_1204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1546[k] = f_14 * ski_979[k]
                    + f_3 * pc_y[k] * sli_1203[k];

        t_1547[k] = f_19 * ski_951[k]
                    + f_1 * slh0_902[k]
                    - f_2 * slh1_902[k]
                    + f_3 * pc_z[k] * sli_1203[k];

        t_1548[k] = pb_y[k] * skk0_1260[k]
                    - f_12 * pc_y[k] * skk1_1260[k];

        t_1549[k] = f_13 * ski_980[k]
                    + f_3 * pc_y[k] * sli_1204[k];
    }

#pragma omp simd aligned(t_1550, t_1551, t_1552, pc_x, pc_y, pc_z, ski_952, ski_982, slh0_906, \
                         slh1_906, sli_1204, sli_1206, sli_1207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1550[k] = f_18 * ski_952[k]
                    + f_3 * pc_z[k] * sli_1204[k];

        t_1551[k] = f_4 * slh0_906[k]
                    - f_5 * slh1_906[k]
                    + f_3 * pc_x[k] * sli_1207[k];

        t_1552[k] = f_13 * ski_982[k]
                    + f_3 * pc_y[k] * sli_1206[k];
    }

#pragma omp simd aligned(t_1553, t_1554, t_1555, pb_y, pc_x, pc_y, pc_z, skk0_1265, ski_955, \
                         skk1_1265, slh0_909, slh1_909, sli_1207, \
                         sli_1210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1553[k] = pb_y[k] * skk0_1265[k]
                    - f_12 * pc_y[k] * skk1_1265[k];

        t_1554[k] = f_6 * slh0_909[k]
                    - f_7 * slh1_909[k]
                    + f_3 * pc_x[k] * sli_1210[k];

        t_1555[k] = f_18 * ski_955[k]
                    + f_3 * pc_z[k] * sli_1207[k];
    }

#pragma omp simd aligned(t_1556, t_1557, t_1558, pb_y, pc_x, pc_y, skk0_1269, ski_985, \
                         skk1_1269, slh0_913, slh1_913, sli_1209, \
                         sli_1214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1556[k] = f_13 * ski_985[k]
                    + f_3 * pc_y[k] * sli_1209[k];

        t_1557[k] = pb_y[k] * skk0_1269[k]
                    - f_12 * pc_y[k] * skk1_1269[k];

        t_1558[k] = f_8 * slh0_913[k]
                    - f_9 * slh1_913[k]
                    + f_3 * pc_x[k] * sli_1214[k];
    }

#pragma omp simd aligned(t_1559, t_1560, t_1561, pc_x, pc_y, pc_z, ski_958, ski_989, slh0_915, \
                         slh1_915, sli_1210, sli_1213, sli_1216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1559[k] = f_18 * ski_958[k]
                    + f_3 * pc_z[k] * sli_1210[k];

        t_1560[k] = f_8 * slh0_915[k]
                    - f_9 * slh1_915[k]
                    + f_3 * pc_x[k] * sli_1216[k];

        t_1561[k] = f_13 * ski_989[k]
                    + f_3 * pc_y[k] * sli_1213[k];
    }

#pragma omp simd aligned(t_1562, t_1563, t_1564, pb_y, pc_x, pc_y, pc_z, skk0_1274, ski_962, \
                         skk1_1274, slh0_918, slh1_918, sli_1214, \
                         sli_1219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1562[k] = pb_y[k] * skk0_1274[k]
                    - f_12 * pc_y[k] * skk1_1274[k];

        t_1563[k] = f_10 * slh0_918[k]
                    - f_11 * slh1_918[k]
                    + f_3 * pc_x[k] * sli_1219[k];

        t_1564[k] = f_18 * ski_962[k]
                    + f_3 * pc_z[k] * sli_1214[k];
    }

#pragma omp simd aligned(t_1565, t_1566, t_1567, pc_x, pc_y, ski_994, slh0_920, slh0_921, \
                         slh1_920, slh1_921, sli_1218, sli_1221, \
                         sli_1222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1565[k] = f_10 * slh0_920[k]
                    - f_11 * slh1_920[k]
                    + f_3 * pc_x[k] * sli_1221[k];

        t_1566[k] = f_10 * slh0_921[k]
                    - f_11 * slh1_921[k]
                    + f_3 * pc_x[k] * sli_1222[k];

        t_1567[k] = f_13 * ski_994[k]
                    + f_3 * pc_y[k] * sli_1218[k];
    }

#pragma omp simd aligned(t_1568, t_1569, t_1570, t_1571, t_1572, t_1573, pb_y, pc_x, pc_y, \
                         skk0_1280, skk1_1280, sli_1225, sli_1226, sli_1227, sli_1228, \
                         sli_1229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1568[k] = pb_y[k] * skk0_1280[k]
                    - f_12 * pc_y[k] * skk1_1280[k];

        t_1569[k] = f_3 * pc_x[k] * sli_1225[k];

        t_1570[k] = f_3 * pc_x[k] * sli_1226[k];

        t_1571[k] = f_3 * pc_x[k] * sli_1227[k];

        t_1572[k] = f_3 * pc_x[k] * sli_1228[k];

        t_1573[k] = f_3 * pc_x[k] * sli_1229[k];
    }

#pragma omp simd aligned(t_1574, t_1575, t_1576, t_1577, pb_y, pc_x, pc_y, pc_z, skk0_1288, \
                         ski_973, ski_1001, skk1_1288, sli_1225, sli_1230, \
                         sli_1231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1574[k] = f_3 * pc_x[k] * sli_1230[k];

        t_1575[k] = f_3 * pc_x[k] * sli_1231[k];

        t_1576[k] = pb_y[k] * skk0_1288[k]
                    + f_18 * ski_1001[k]
                    - f_12 * pc_y[k] * skk1_1288[k];

        t_1577[k] = f_18 * ski_973[k]
                    + f_3 * pc_z[k] * sli_1225[k];
    }

#pragma omp simd aligned(t_1578, t_1579, t_1580, pb_y, pc_y, skk0_1290, skk0_1291, skk0_1292, \
                         ski_1003, ski_1004, ski_1005, skk1_1290, skk1_1291, \
                         skk1_1292 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1578[k] = pb_y[k] * skk0_1290[k]
                    + f_17 * ski_1003[k]
                    - f_12 * pc_y[k] * skk1_1290[k];

        t_1579[k] = pb_y[k] * skk0_1291[k]
                    + f_16 * ski_1004[k]
                    - f_12 * pc_y[k] * skk1_1291[k];

        t_1580[k] = pb_y[k] * skk0_1292[k]
                    + f_15 * ski_1005[k]
                    - f_12 * pc_y[k] * skk1_1292[k];
    }

#pragma omp simd aligned(t_1581, t_1582, t_1583, pb_y, pc_y, skk0_1293, skk0_1295, ski_1006, \
                         ski_1007, skk1_1293, skk1_1295, sli_1231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1581[k] = pb_y[k] * skk0_1293[k]
                    + f_14 * ski_1006[k]
                    - f_12 * pc_y[k] * skk1_1293[k];

        t_1582[k] = f_13 * ski_1007[k]
                    + f_3 * pc_y[k] * sli_1231[k];

        t_1583[k] = pb_y[k] * skk0_1295[k]
                    - f_12 * pc_y[k] * skk1_1295[k];
    }

#pragma omp simd aligned(t_1584, t_1585, t_1586, t_1587, t_1588, pc_x, pc_y, pc_z, ski_980, \
                         slh0_924, slh0_927, slh1_924, slh1_927, sli_1232, sli_1234, \
                         sli_1235 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1584[k] = f_1 * slh0_924[k]
                    - f_2 * slh1_924[k]
                    + f_3 * pc_x[k] * sli_1232[k];

        t_1585[k] = f_3 * pc_y[k] * sli_1232[k];

        t_1586[k] = f_0 * ski_980[k]
                    + f_3 * pc_z[k] * sli_1232[k];

        t_1587[k] = f_4 * slh0_927[k]
                    - f_5 * slh1_927[k]
                    + f_3 * pc_x[k] * sli_1235[k];

        t_1588[k] = f_3 * pc_y[k] * sli_1234[k];
    }

#pragma omp simd aligned(t_1589, t_1590, t_1591, t_1592, pc_x, pc_y, pc_z, ski_983, slh0_929, \
                         slh0_930, slh1_929, slh1_930, sli_1235, sli_1237, \
                         sli_1238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1589[k] = f_4 * slh0_929[k]
                    - f_5 * slh1_929[k]
                    + f_3 * pc_x[k] * sli_1237[k];

        t_1590[k] = f_6 * slh0_930[k]
                    - f_7 * slh1_930[k]
                    + f_3 * pc_x[k] * sli_1238[k];

        t_1591[k] = f_0 * ski_983[k]
                    + f_3 * pc_z[k] * sli_1235[k];

        t_1592[k] = f_3 * pc_y[k] * sli_1237[k];
    }

#pragma omp simd aligned(t_1593, t_1594, t_1595, pc_x, pc_z, ski_986, slh0_933, slh0_934, \
                         slh1_933, slh1_934, sli_1238, sli_1241, \
                         sli_1242 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1593[k] = f_6 * slh0_933[k]
                    - f_7 * slh1_933[k]
                    + f_3 * pc_x[k] * sli_1241[k];

        t_1594[k] = f_8 * slh0_934[k]
                    - f_9 * slh1_934[k]
                    + f_3 * pc_x[k] * sli_1242[k];

        t_1595[k] = f_0 * ski_986[k]
                    + f_3 * pc_z[k] * sli_1238[k];
    }

#pragma omp simd aligned(t_1596, t_1597, t_1598, t_1599, pc_x, pc_y, slh0_936, slh0_938, \
                         slh0_939, slh1_936, slh1_938, slh1_939, sli_1241, sli_1244, sli_1246, \
                         sli_1247 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1596[k] = f_8 * slh0_936[k]
                    - f_9 * slh1_936[k]
                    + f_3 * pc_x[k] * sli_1244[k];

        t_1597[k] = f_3 * pc_y[k] * sli_1241[k];

        t_1598[k] = f_8 * slh0_938[k]
                    - f_9 * slh1_938[k]
                    + f_3 * pc_x[k] * sli_1246[k];

        t_1599[k] = f_10 * slh0_939[k]
                    - f_11 * slh1_939[k]
                    + f_3 * pc_x[k] * sli_1247[k];
    }

#pragma omp simd aligned(t_1600, t_1601, t_1602, t_1603, pc_x, pc_y, pc_z, ski_990, slh0_941, \
                         slh0_942, slh1_941, slh1_942, sli_1242, sli_1246, sli_1249, \
                         sli_1250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1600[k] = f_0 * ski_990[k]
                    + f_3 * pc_z[k] * sli_1242[k];

        t_1601[k] = f_10 * slh0_941[k]
                    - f_11 * slh1_941[k]
                    + f_3 * pc_x[k] * sli_1249[k];

        t_1602[k] = f_10 * slh0_942[k]
                    - f_11 * slh1_942[k]
                    + f_3 * pc_x[k] * sli_1250[k];

        t_1603[k] = f_3 * pc_y[k] * sli_1246[k];
    }

#pragma omp simd aligned(t_1604, t_1605, t_1606, t_1607, t_1608, t_1609, pc_x, slh0_944, \
                         slh1_944, sli_1252, sli_1253, sli_1254, sli_1255, sli_1256, \
                         sli_1257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1604[k] = f_10 * slh0_944[k]
                    - f_11 * slh1_944[k]
                    + f_3 * pc_x[k] * sli_1252[k];

        t_1605[k] = f_3 * pc_x[k] * sli_1253[k];

        t_1606[k] = f_3 * pc_x[k] * sli_1254[k];

        t_1607[k] = f_3 * pc_x[k] * sli_1255[k];

        t_1608[k] = f_3 * pc_x[k] * sli_1256[k];

        t_1609[k] = f_3 * pc_x[k] * sli_1257[k];
    }

#pragma omp simd aligned(t_1610, t_1611, t_1612, t_1613, pc_x, pc_y, pc_z, ski_1001, slh0_939, \
                         slh1_939, sli_1253, sli_1258, sli_1259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1610[k] = f_3 * pc_x[k] * sli_1258[k];

        t_1611[k] = f_3 * pc_x[k] * sli_1259[k];

        t_1612[k] = f_1 * slh0_939[k]
                    - f_2 * slh1_939[k]
                    + f_3 * pc_y[k] * sli_1253[k];

        t_1613[k] = f_0 * ski_1001[k]
                    + f_3 * pc_z[k] * sli_1253[k];
    }

#pragma omp simd aligned(t_1614, t_1615, t_1616, pc_y, slh0_941, slh0_942, slh0_943, slh1_941, \
                         slh1_942, slh1_943, sli_1255, sli_1256, \
                         sli_1257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1614[k] = f_4 * slh0_941[k]
                    - f_5 * slh1_941[k]
                    + f_3 * pc_y[k] * sli_1255[k];

        t_1615[k] = f_6 * slh0_942[k]
                    - f_7 * slh1_942[k]
                    + f_3 * pc_y[k] * sli_1256[k];

        t_1616[k] = f_8 * slh0_943[k]
                    - f_9 * slh1_943[k]
                    + f_3 * pc_y[k] * sli_1257[k];
    }
}

static auto
compute_prim_slk_three_center_electron_repulsion_0_piece14(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pc,
                                                           const size_t ski, const size_t slh0,
                                                           const size_t slh1, const size_t sli,
                                                           const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
    const auto f_1 = 3.0 / gamma;
    const auto f_2 = 3.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_10 = 0.5 / gamma;
    const auto f_11 = 0.5 * p / (gamma * q);

    auto *t_1617 = buffer.data(target + 1617);
    auto *t_1618 = buffer.data(target + 1618);
    auto *t_1619 = buffer.data(target + 1619);

    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ski_1007 = buffer.data(ski + 1007);

    const auto *slh0_944 = buffer.data(slh0 + 944);

    const auto *slh1_944 = buffer.data(slh1 + 944);

    const auto *sli_1258 = buffer.data(sli + 1258);
    const auto *sli_1259 = buffer.data(sli + 1259);

#pragma omp simd aligned(t_1617, t_1618, t_1619, pc_y, pc_z, ski_1007, slh0_944, slh1_944, \
                         sli_1258, sli_1259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1617[k] = f_10 * slh0_944[k]
                    - f_11 * slh1_944[k]
                    + f_3 * pc_y[k] * sli_1258[k];

        t_1618[k] = f_3 * pc_y[k] * sli_1259[k];

        t_1619[k] = f_0 * ski_1007[k]
                    + f_1 * slh0_944[k]
                    - f_2 * slh1_944[k]
                    + f_3 * pc_z[k] * sli_1259[k];
    }
}

auto
compute_prim_slk_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t skk0, const size_t ski,
                                                   const size_t skk1, const size_t slh0,
                                                   const size_t slh1, const size_t sli,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_slk_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, skk0, ski,
                                                              skk1, slh0, slh1, sli, ncols,
                                                              gamma, p, q);

    compute_prim_slk_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, skk0, ski,
                                                              skk1, slh0, slh1, sli, ncols,
                                                              gamma, p, q);

    compute_prim_slk_three_center_electron_repulsion_0_piece2(buffer, target, pb, pc, skk0, ski,
                                                              skk1, slh0, slh1, sli, ncols,
                                                              gamma, p, q);

    compute_prim_slk_three_center_electron_repulsion_0_piece3(buffer, target, pb, pc, skk0, ski,
                                                              skk1, slh0, slh1, sli, ncols,
                                                              gamma, p, q);

    compute_prim_slk_three_center_electron_repulsion_0_piece4(buffer, target, pb, pc, skk0, ski,
                                                              skk1, slh0, slh1, sli, ncols,
                                                              gamma, p, q);

    compute_prim_slk_three_center_electron_repulsion_0_piece5(buffer, target, pb, pc, skk0, ski,
                                                              skk1, slh0, slh1, sli, ncols,
                                                              gamma, p, q);

    compute_prim_slk_three_center_electron_repulsion_0_piece6(buffer, target, pb, pc, skk0, ski,
                                                              skk1, slh0, slh1, sli, ncols,
                                                              gamma, p, q);

    compute_prim_slk_three_center_electron_repulsion_0_piece7(buffer, target, pb, pc, skk0, ski,
                                                              skk1, slh0, slh1, sli, ncols,
                                                              gamma, p, q);

    compute_prim_slk_three_center_electron_repulsion_0_piece8(buffer, target, pb, pc, skk0, ski,
                                                              skk1, slh0, slh1, sli, ncols,
                                                              gamma, p, q);

    compute_prim_slk_three_center_electron_repulsion_0_piece9(buffer, target, pb, pc, skk0, ski,
                                                              skk1, sli, ncols, gamma, p, q);

    compute_prim_slk_three_center_electron_repulsion_0_piece10(buffer, target, pb, pc, skk0,
                                                               ski, skk1, sli, ncols, gamma, p,
                                                               q);

    compute_prim_slk_three_center_electron_repulsion_0_piece11(buffer, target, pb, pc, skk0,
                                                               ski, skk1, slh0, slh1, sli,
                                                               ncols, gamma, p, q);

    compute_prim_slk_three_center_electron_repulsion_0_piece12(buffer, target, pc, ski, slh0,
                                                               slh1, sli, ncols, gamma, p, q);

    compute_prim_slk_three_center_electron_repulsion_0_piece13(buffer, target, pb, pc, skk0,
                                                               ski, skk1, slh0, slh1, sli,
                                                               ncols, gamma, p, q);

    compute_prim_slk_three_center_electron_repulsion_0_piece14(buffer, target, pc, ski, slh0,
                                                               slh1, sli, ncols, gamma, p, q);
}

}  // namespace simdt3ceri
