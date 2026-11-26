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


#include "SimdThreeCenterElectronRepulsionVrrRecSLI.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_sli_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t ski0,
                                                          const size_t skh, const size_t ski1,
                                                          const size_t slg0, const size_t slg1,
                                                          const size_t slh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.5 / gamma;
    const auto f_5 = 1.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 0.5 / gamma;
    const auto f_9 = 0.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 3.5 / q;
    const auto f_16 = 3.0 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ski0_0 = buffer.data(ski0 + 0);
    const auto *ski0_3 = buffer.data(ski0 + 3);
    const auto *ski0_5 = buffer.data(ski0 + 5);
    const auto *ski0_6 = buffer.data(ski0 + 6);
    const auto *ski0_9 = buffer.data(ski0 + 9);
    const auto *ski0_10 = buffer.data(ski0 + 10);
    const auto *ski0_12 = buffer.data(ski0 + 12);
    const auto *ski0_14 = buffer.data(ski0 + 14);
    const auto *ski0_21 = buffer.data(ski0 + 21);
    const auto *ski0_27 = buffer.data(ski0 + 27);
    const auto *ski0_31 = buffer.data(ski0 + 31);
    const auto *ski0_34 = buffer.data(ski0 + 34);
    const auto *ski0_38 = buffer.data(ski0 + 38);
    const auto *ski0_56 = buffer.data(ski0 + 56);
    const auto *ski0_61 = buffer.data(ski0 + 61);
    const auto *ski0_65 = buffer.data(ski0 + 65);

    const auto *skh_0 = buffer.data(skh + 0);
    const auto *skh_1 = buffer.data(skh + 1);
    const auto *skh_2 = buffer.data(skh + 2);
    const auto *skh_3 = buffer.data(skh + 3);
    const auto *skh_5 = buffer.data(skh + 5);
    const auto *skh_6 = buffer.data(skh + 6);
    const auto *skh_7 = buffer.data(skh + 7);
    const auto *skh_8 = buffer.data(skh + 8);
    const auto *skh_9 = buffer.data(skh + 9);
    const auto *skh_10 = buffer.data(skh + 10);
    const auto *skh_12 = buffer.data(skh + 12);
    const auto *skh_14 = buffer.data(skh + 14);
    const auto *skh_15 = buffer.data(skh + 15);
    const auto *skh_16 = buffer.data(skh + 16);
    const auto *skh_17 = buffer.data(skh + 17);
    const auto *skh_18 = buffer.data(skh + 18);
    const auto *skh_19 = buffer.data(skh + 19);
    const auto *skh_20 = buffer.data(skh + 20);
    const auto *skh_21 = buffer.data(skh + 21);
    const auto *skh_23 = buffer.data(skh + 23);
    const auto *skh_24 = buffer.data(skh + 24);
    const auto *skh_26 = buffer.data(skh + 26);
    const auto *skh_27 = buffer.data(skh + 27);
    const auto *skh_30 = buffer.data(skh + 30);
    const auto *skh_36 = buffer.data(skh + 36);
    const auto *skh_37 = buffer.data(skh + 37);
    const auto *skh_38 = buffer.data(skh + 38);
    const auto *skh_39 = buffer.data(skh + 39);
    const auto *skh_40 = buffer.data(skh + 40);
    const auto *skh_41 = buffer.data(skh + 41);
    const auto *skh_42 = buffer.data(skh + 42);
    const auto *skh_44 = buffer.data(skh + 44);
    const auto *skh_47 = buffer.data(skh + 47);
    const auto *skh_57 = buffer.data(skh + 57);
    const auto *skh_58 = buffer.data(skh + 58);
    const auto *skh_59 = buffer.data(skh + 59);
    const auto *skh_60 = buffer.data(skh + 60);
    const auto *skh_61 = buffer.data(skh + 61);
    const auto *skh_62 = buffer.data(skh + 62);
    const auto *skh_63 = buffer.data(skh + 63);
    const auto *skh_66 = buffer.data(skh + 66);
    const auto *skh_68 = buffer.data(skh + 68);
    const auto *skh_69 = buffer.data(skh + 69);
    const auto *skh_72 = buffer.data(skh + 72);
    const auto *skh_73 = buffer.data(skh + 73);
    const auto *skh_75 = buffer.data(skh + 75);
    const auto *skh_77 = buffer.data(skh + 77);
    const auto *skh_78 = buffer.data(skh + 78);
    const auto *skh_79 = buffer.data(skh + 79);
    const auto *skh_80 = buffer.data(skh + 80);
    const auto *skh_81 = buffer.data(skh + 81);
    const auto *skh_82 = buffer.data(skh + 82);
    const auto *skh_83 = buffer.data(skh + 83);

    const auto *ski1_0 = buffer.data(ski1 + 0);
    const auto *ski1_3 = buffer.data(ski1 + 3);
    const auto *ski1_5 = buffer.data(ski1 + 5);
    const auto *ski1_6 = buffer.data(ski1 + 6);
    const auto *ski1_9 = buffer.data(ski1 + 9);
    const auto *ski1_10 = buffer.data(ski1 + 10);
    const auto *ski1_12 = buffer.data(ski1 + 12);
    const auto *ski1_14 = buffer.data(ski1 + 14);
    const auto *ski1_21 = buffer.data(ski1 + 21);
    const auto *ski1_27 = buffer.data(ski1 + 27);
    const auto *ski1_31 = buffer.data(ski1 + 31);
    const auto *ski1_34 = buffer.data(ski1 + 34);
    const auto *ski1_38 = buffer.data(ski1 + 38);
    const auto *ski1_56 = buffer.data(ski1 + 56);
    const auto *ski1_61 = buffer.data(ski1 + 61);
    const auto *ski1_65 = buffer.data(ski1 + 65);

    const auto *slg0_0 = buffer.data(slg0 + 0);
    const auto *slg0_3 = buffer.data(slg0 + 3);
    const auto *slg0_5 = buffer.data(slg0 + 5);
    const auto *slg0_6 = buffer.data(slg0 + 6);
    const auto *slg0_9 = buffer.data(slg0 + 9);
    const auto *slg0_10 = buffer.data(slg0 + 10);
    const auto *slg0_12 = buffer.data(slg0 + 12);
    const auto *slg0_13 = buffer.data(slg0 + 13);
    const auto *slg0_14 = buffer.data(slg0 + 14);
    const auto *slg0_25 = buffer.data(slg0 + 25);
    const auto *slg0_27 = buffer.data(slg0 + 27);
    const auto *slg0_28 = buffer.data(slg0 + 28);
    const auto *slg0_29 = buffer.data(slg0 + 29);
    const auto *slg0_42 = buffer.data(slg0 + 42);
    const auto *slg0_43 = buffer.data(slg0 + 43);
    const auto *slg0_44 = buffer.data(slg0 + 44);
    const auto *slg0_45 = buffer.data(slg0 + 45);
    const auto *slg0_48 = buffer.data(slg0 + 48);
    const auto *slg0_50 = buffer.data(slg0 + 50);
    const auto *slg0_51 = buffer.data(slg0 + 51);
    const auto *slg0_54 = buffer.data(slg0 + 54);
    const auto *slg0_55 = buffer.data(slg0 + 55);
    const auto *slg0_57 = buffer.data(slg0 + 57);
    const auto *slg0_58 = buffer.data(slg0 + 58);
    const auto *slg0_59 = buffer.data(slg0 + 59);

    const auto *slg1_0 = buffer.data(slg1 + 0);
    const auto *slg1_3 = buffer.data(slg1 + 3);
    const auto *slg1_5 = buffer.data(slg1 + 5);
    const auto *slg1_6 = buffer.data(slg1 + 6);
    const auto *slg1_9 = buffer.data(slg1 + 9);
    const auto *slg1_10 = buffer.data(slg1 + 10);
    const auto *slg1_12 = buffer.data(slg1 + 12);
    const auto *slg1_13 = buffer.data(slg1 + 13);
    const auto *slg1_14 = buffer.data(slg1 + 14);
    const auto *slg1_25 = buffer.data(slg1 + 25);
    const auto *slg1_27 = buffer.data(slg1 + 27);
    const auto *slg1_28 = buffer.data(slg1 + 28);
    const auto *slg1_29 = buffer.data(slg1 + 29);
    const auto *slg1_42 = buffer.data(slg1 + 42);
    const auto *slg1_43 = buffer.data(slg1 + 43);
    const auto *slg1_44 = buffer.data(slg1 + 44);
    const auto *slg1_45 = buffer.data(slg1 + 45);
    const auto *slg1_48 = buffer.data(slg1 + 48);
    const auto *slg1_50 = buffer.data(slg1 + 50);
    const auto *slg1_51 = buffer.data(slg1 + 51);
    const auto *slg1_54 = buffer.data(slg1 + 54);
    const auto *slg1_55 = buffer.data(slg1 + 55);
    const auto *slg1_57 = buffer.data(slg1 + 57);
    const auto *slg1_58 = buffer.data(slg1 + 58);
    const auto *slg1_59 = buffer.data(slg1 + 59);

    const auto *slh_0 = buffer.data(slh + 0);
    const auto *slh_2 = buffer.data(slh + 2);
    const auto *slh_3 = buffer.data(slh + 3);
    const auto *slh_5 = buffer.data(slh + 5);
    const auto *slh_6 = buffer.data(slh + 6);
    const auto *slh_9 = buffer.data(slh + 9);
    const auto *slh_10 = buffer.data(slh + 10);
    const auto *slh_12 = buffer.data(slh + 12);
    const auto *slh_14 = buffer.data(slh + 14);
    const auto *slh_15 = buffer.data(slh + 15);
    const auto *slh_16 = buffer.data(slh + 16);
    const auto *slh_17 = buffer.data(slh + 17);
    const auto *slh_18 = buffer.data(slh + 18);
    const auto *slh_19 = buffer.data(slh + 19);
    const auto *slh_20 = buffer.data(slh + 20);
    const auto *slh_21 = buffer.data(slh + 21);
    const auto *slh_23 = buffer.data(slh + 23);
    const auto *slh_24 = buffer.data(slh + 24);
    const auto *slh_26 = buffer.data(slh + 26);
    const auto *slh_27 = buffer.data(slh + 27);
    const auto *slh_30 = buffer.data(slh + 30);
    const auto *slh_36 = buffer.data(slh + 36);
    const auto *slh_37 = buffer.data(slh + 37);
    const auto *slh_38 = buffer.data(slh + 38);
    const auto *slh_39 = buffer.data(slh + 39);
    const auto *slh_40 = buffer.data(slh + 40);
    const auto *slh_41 = buffer.data(slh + 41);
    const auto *slh_42 = buffer.data(slh + 42);
    const auto *slh_44 = buffer.data(slh + 44);
    const auto *slh_45 = buffer.data(slh + 45);
    const auto *slh_47 = buffer.data(slh + 47);
    const auto *slh_48 = buffer.data(slh + 48);
    const auto *slh_51 = buffer.data(slh + 51);
    const auto *slh_57 = buffer.data(slh + 57);
    const auto *slh_58 = buffer.data(slh + 58);
    const auto *slh_59 = buffer.data(slh + 59);
    const auto *slh_60 = buffer.data(slh + 60);
    const auto *slh_61 = buffer.data(slh + 61);
    const auto *slh_62 = buffer.data(slh + 62);
    const auto *slh_63 = buffer.data(slh + 63);
    const auto *slh_65 = buffer.data(slh + 65);
    const auto *slh_66 = buffer.data(slh + 66);
    const auto *slh_68 = buffer.data(slh + 68);
    const auto *slh_69 = buffer.data(slh + 69);
    const auto *slh_72 = buffer.data(slh + 72);
    const auto *slh_73 = buffer.data(slh + 73);
    const auto *slh_75 = buffer.data(slh + 75);
    const auto *slh_77 = buffer.data(slh + 77);
    const auto *slh_78 = buffer.data(slh + 78);
    const auto *slh_79 = buffer.data(slh + 79);
    const auto *slh_80 = buffer.data(slh + 80);
    const auto *slh_81 = buffer.data(slh + 81);
    const auto *slh_82 = buffer.data(slh + 82);
    const auto *slh_83 = buffer.data(slh + 83);
    const auto *slh_84 = buffer.data(slh + 84);
    const auto *slh_86 = buffer.data(slh + 86);
    const auto *slh_87 = buffer.data(slh + 87);
    const auto *slh_89 = buffer.data(slh + 89);
    const auto *slh_90 = buffer.data(slh + 90);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, skh_0, skh_3, slg0_0, slg0_3, \
                         slg1_0, slg1_3, slh_0, slh_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * skh_0[k]
                 + f_1 * slg0_0[k]
                 - f_2 * slg1_0[k]
                 + f_3 * pc_x[k] * slh_0[k];

        t_1[k] = f_3 * pc_y[k] * slh_0[k];

        t_2[k] = f_3 * pc_z[k] * slh_0[k];

        t_3[k] = f_0 * skh_3[k]
                 + f_4 * slg0_3[k]
                 - f_5 * slg1_3[k]
                 + f_3 * pc_x[k] * slh_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pc_x, pc_y, skh_5, skh_6, slg0_5, slg0_6, slg1_5, \
                         slg1_6, slh_2, slh_5, slh_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * pc_y[k] * slh_2[k];

        t_5[k] = f_0 * skh_5[k]
                 + f_4 * slg0_5[k]
                 - f_5 * slg1_5[k]
                 + f_3 * pc_x[k] * slh_5[k];

        t_6[k] = f_0 * skh_6[k]
                 + f_6 * slg0_6[k]
                 - f_7 * slg1_6[k]
                 + f_3 * pc_x[k] * slh_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pc_x, pc_y, pc_z, skh_9, slg0_9, slg1_9, slh_3, slh_5, \
                         slh_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_3 * pc_z[k] * slh_3[k];

        t_8[k] = f_3 * pc_y[k] * slh_5[k];

        t_9[k] = f_0 * skh_9[k]
                 + f_6 * slg0_9[k]
                 - f_7 * slg1_9[k]
                 + f_3 * pc_x[k] * slh_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pc_x, pc_z, skh_10, skh_12, slg0_10, slg0_12, \
                         slg1_10, slg1_12, slh_6, slh_10, slh_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * skh_10[k]
                  + f_8 * slg0_10[k]
                  - f_9 * slg1_10[k]
                  + f_3 * pc_x[k] * slh_10[k];

        t_11[k] = f_3 * pc_z[k] * slh_6[k];

        t_12[k] = f_0 * skh_12[k]
                  + f_8 * slg0_12[k]
                  - f_9 * slg1_12[k]
                  + f_3 * pc_x[k] * slh_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pc_x, pc_y, skh_14, skh_15, skh_16, slg0_14, \
                         slg1_14, slh_9, slh_14, slh_15, slh_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_3 * pc_y[k] * slh_9[k];

        t_14[k] = f_0 * skh_14[k]
                  + f_8 * slg0_14[k]
                  - f_9 * slg1_14[k]
                  + f_3 * pc_x[k] * slh_14[k];

        t_15[k] = f_0 * skh_15[k]
                  + f_3 * pc_x[k] * slh_15[k];

        t_16[k] = f_0 * skh_16[k]
                  + f_3 * pc_x[k] * slh_16[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pc_x, skh_17, skh_18, skh_19, skh_20, slh_17, \
                         slh_18, slh_19, slh_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_0 * skh_17[k]
                  + f_3 * pc_x[k] * slh_17[k];

        t_18[k] = f_0 * skh_18[k]
                  + f_3 * pc_x[k] * slh_18[k];

        t_19[k] = f_0 * skh_19[k]
                  + f_3 * pc_x[k] * slh_19[k];

        t_20[k] = f_0 * skh_20[k]
                  + f_3 * pc_x[k] * slh_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pc_y, pc_z, slg0_10, slg0_12, slg0_13, \
                         slg1_10, slg1_12, slg1_13, slh_15, slh_17, \
                         slh_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * slg0_10[k]
                  - f_2 * slg1_10[k]
                  + f_3 * pc_y[k] * slh_15[k];

        t_22[k] = f_3 * pc_z[k] * slh_15[k];

        t_23[k] = f_4 * slg0_12[k]
                  - f_5 * slg1_12[k]
                  + f_3 * pc_y[k] * slh_17[k];

        t_24[k] = f_6 * slg0_13[k]
                  - f_7 * slg1_13[k]
                  + f_3 * pc_y[k] * slh_18[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pb_y, pc_y, pc_z, ski0_0, skh_0, \
                         ski1_0, slg0_14, slg1_14, slh_19, slh_20, \
                         slh_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_8 * slg0_14[k]
                  - f_9 * slg1_14[k]
                  + f_3 * pc_y[k] * slh_19[k];

        t_26[k] = f_3 * pc_y[k] * slh_20[k];

        t_27[k] = f_1 * slg0_14[k]
                  - f_2 * slg1_14[k]
                  + f_3 * pc_z[k] * slh_20[k];

        t_28[k] = pb_y[k] * ski0_0[k]
                  - f_10 * pc_y[k] * ski1_0[k];

        t_29[k] = f_11 * skh_0[k]
                  + f_3 * pc_y[k] * slh_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pb_y, pc_y, pc_z, ski0_3, ski0_5, skh_1, \
                         skh_2, ski1_3, ski1_5, slh_21, slh_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_3 * pc_z[k] * slh_21[k];

        t_31[k] = pb_y[k] * ski0_3[k]
                  + f_12 * skh_1[k]
                  - f_10 * pc_y[k] * ski1_3[k];

        t_32[k] = f_11 * skh_2[k]
                  + f_3 * pc_y[k] * slh_23[k];

        t_33[k] = pb_y[k] * ski0_5[k]
                  - f_10 * pc_y[k] * ski1_5[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pb_y, pc_y, pc_z, ski0_6, ski0_9, skh_3, \
                         skh_5, ski1_6, ski1_9, slh_24, slh_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pb_y[k] * ski0_6[k]
                  + f_13 * skh_3[k]
                  - f_10 * pc_y[k] * ski1_6[k];

        t_35[k] = f_3 * pc_z[k] * slh_24[k];

        t_36[k] = f_11 * skh_5[k]
                  + f_3 * pc_y[k] * slh_26[k];

        t_37[k] = pb_y[k] * ski0_9[k]
                  - f_10 * pc_y[k] * ski1_9[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pb_y, pc_y, pc_z, ski0_10, ski0_12, skh_6, \
                         skh_8, skh_9, ski1_10, ski1_12, slh_27, \
                         slh_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pb_y[k] * ski0_10[k]
                  + f_14 * skh_6[k]
                  - f_10 * pc_y[k] * ski1_10[k];

        t_39[k] = f_3 * pc_z[k] * slh_27[k];

        t_40[k] = pb_y[k] * ski0_12[k]
                  + f_12 * skh_8[k]
                  - f_10 * pc_y[k] * ski1_12[k];

        t_41[k] = f_11 * skh_9[k]
                  + f_3 * pc_y[k] * slh_30[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pb_y, pc_x, pc_y, ski0_14, skh_36, skh_37, \
                         skh_38, ski1_14, slh_36, slh_37, slh_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_y[k] * ski0_14[k]
                  - f_10 * pc_y[k] * ski1_14[k];

        t_43[k] = f_15 * skh_36[k]
                  + f_3 * pc_x[k] * slh_36[k];

        t_44[k] = f_15 * skh_37[k]
                  + f_3 * pc_x[k] * slh_37[k];

        t_45[k] = f_15 * skh_38[k]
                  + f_3 * pc_x[k] * slh_38[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pc_x, pc_y, skh_15, skh_39, skh_40, skh_41, \
                         slg0_25, slg1_25, slh_36, slh_39, slh_40, \
                         slh_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_15 * skh_39[k]
                  + f_3 * pc_x[k] * slh_39[k];

        t_47[k] = f_15 * skh_40[k]
                  + f_3 * pc_x[k] * slh_40[k];

        t_48[k] = f_15 * skh_41[k]
                  + f_3 * pc_x[k] * slh_41[k];

        t_49[k] = f_11 * skh_15[k]
                  + f_1 * slg0_25[k]
                  - f_2 * slg1_25[k]
                  + f_3 * pc_y[k] * slh_36[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, pc_y, pc_z, skh_17, skh_18, slg0_27, slg0_28, \
                         slg1_27, slg1_28, slh_36, slh_38, slh_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_3 * pc_z[k] * slh_36[k];

        t_51[k] = f_11 * skh_17[k]
                  + f_4 * slg0_27[k]
                  - f_5 * slg1_27[k]
                  + f_3 * pc_y[k] * slh_38[k];

        t_52[k] = f_11 * skh_18[k]
                  + f_6 * slg0_28[k]
                  - f_7 * slg1_28[k]
                  + f_3 * pc_y[k] * slh_39[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pb_y, pc_y, ski0_27, skh_19, skh_20, ski1_27, \
                         slg0_29, slg1_29, slh_40, slh_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_11 * skh_19[k]
                  + f_8 * slg0_29[k]
                  - f_9 * slg1_29[k]
                  + f_3 * pc_y[k] * slh_40[k];

        t_54[k] = f_11 * skh_20[k]
                  + f_3 * pc_y[k] * slh_41[k];

        t_55[k] = pb_y[k] * ski0_27[k]
                  - f_10 * pc_y[k] * ski1_27[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, pb_z, pc_y, pc_z, ski0_0, ski0_3, \
                         skh_0, ski1_0, ski1_3, slh_42, slh_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_z[k] * ski0_0[k]
                  - f_10 * pc_z[k] * ski1_0[k];

        t_57[k] = f_3 * pc_y[k] * slh_42[k];

        t_58[k] = f_11 * skh_0[k]
                  + f_3 * pc_z[k] * slh_42[k];

        t_59[k] = pb_z[k] * ski0_3[k]
                  - f_10 * pc_z[k] * ski1_3[k];

        t_60[k] = f_3 * pc_y[k] * slh_44[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pb_z, pc_y, pc_z, ski0_5, ski0_6, skh_2, \
                         skh_3, ski1_5, ski1_6, slh_45, slh_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = pb_z[k] * ski0_5[k]
                  + f_12 * skh_2[k]
                  - f_10 * pc_z[k] * ski1_5[k];

        t_62[k] = pb_z[k] * ski0_6[k]
                  - f_10 * pc_z[k] * ski1_6[k];

        t_63[k] = f_11 * skh_3[k]
                  + f_3 * pc_z[k] * slh_45[k];

        t_64[k] = f_3 * pc_y[k] * slh_47[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pb_z, pc_z, ski0_9, ski0_10, ski0_12, skh_5, \
                         skh_6, skh_7, ski1_9, ski1_10, ski1_12, \
                         slh_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = pb_z[k] * ski0_9[k]
                  + f_13 * skh_5[k]
                  - f_10 * pc_z[k] * ski1_9[k];

        t_66[k] = pb_z[k] * ski0_10[k]
                  - f_10 * pc_z[k] * ski1_10[k];

        t_67[k] = f_11 * skh_6[k]
                  + f_3 * pc_z[k] * slh_48[k];

        t_68[k] = pb_z[k] * ski0_12[k]
                  + f_12 * skh_7[k]
                  - f_10 * pc_z[k] * ski1_12[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pb_z, pc_x, pc_y, pc_z, ski0_14, skh_9, \
                         skh_57, skh_58, ski1_14, slh_51, slh_57, \
                         slh_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_3 * pc_y[k] * slh_51[k];

        t_70[k] = pb_z[k] * ski0_14[k]
                  + f_14 * skh_9[k]
                  - f_10 * pc_z[k] * ski1_14[k];

        t_71[k] = f_15 * skh_57[k]
                  + f_3 * pc_x[k] * slh_57[k];

        t_72[k] = f_15 * skh_58[k]
                  + f_3 * pc_x[k] * slh_58[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pc_x, skh_59, skh_60, skh_61, skh_62, slh_59, \
                         slh_60, slh_61, slh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_15 * skh_59[k]
                  + f_3 * pc_x[k] * slh_59[k];

        t_74[k] = f_15 * skh_60[k]
                  + f_3 * pc_x[k] * slh_60[k];

        t_75[k] = f_15 * skh_61[k]
                  + f_3 * pc_x[k] * slh_61[k];

        t_76[k] = f_15 * skh_62[k]
                  + f_3 * pc_x[k] * slh_62[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, pb_z, pc_y, pc_z, ski0_21, skh_15, ski1_21, \
                         slg0_42, slg1_42, slh_57, slh_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = pb_z[k] * ski0_21[k]
                  - f_10 * pc_z[k] * ski1_21[k];

        t_78[k] = f_11 * skh_15[k]
                  + f_3 * pc_z[k] * slh_57[k];

        t_79[k] = f_4 * slg0_42[k]
                  - f_5 * slg1_42[k]
                  + f_3 * pc_y[k] * slh_59[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pc_y, pc_z, skh_20, slg0_43, slg0_44, \
                         slg1_43, slg1_44, slh_60, slh_61, slh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_6 * slg0_43[k]
                  - f_7 * slg1_43[k]
                  + f_3 * pc_y[k] * slh_60[k];

        t_81[k] = f_8 * slg0_44[k]
                  - f_9 * slg1_44[k]
                  + f_3 * pc_y[k] * slh_61[k];

        t_82[k] = f_3 * pc_y[k] * slh_62[k];

        t_83[k] = f_11 * skh_20[k]
                  + f_1 * slg0_44[k]
                  - f_2 * slg1_44[k]
                  + f_3 * pc_z[k] * slh_62[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pc_x, pc_y, pc_z, skh_21, skh_63, skh_66, \
                         slg0_45, slg0_48, slg1_45, slg1_48, slh_63, \
                         slh_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_16 * skh_63[k]
                  + f_1 * slg0_45[k]
                  - f_2 * slg1_45[k]
                  + f_3 * pc_x[k] * slh_63[k];

        t_85[k] = f_12 * skh_21[k]
                  + f_3 * pc_y[k] * slh_63[k];

        t_86[k] = f_3 * pc_z[k] * slh_63[k];

        t_87[k] = f_16 * skh_66[k]
                  + f_4 * slg0_48[k]
                  - f_5 * slg1_48[k]
                  + f_3 * pc_x[k] * slh_66[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, pc_x, pc_y, skh_23, skh_68, skh_69, slg0_50, \
                         slg0_51, slg1_50, slg1_51, slh_65, slh_68, \
                         slh_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_12 * skh_23[k]
                  + f_3 * pc_y[k] * slh_65[k];

        t_89[k] = f_16 * skh_68[k]
                  + f_4 * slg0_50[k]
                  - f_5 * slg1_50[k]
                  + f_3 * pc_x[k] * slh_68[k];

        t_90[k] = f_16 * skh_69[k]
                  + f_6 * slg0_51[k]
                  - f_7 * slg1_51[k]
                  + f_3 * pc_x[k] * slh_69[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, pc_x, pc_y, pc_z, skh_26, skh_72, slg0_54, slg1_54, \
                         slh_66, slh_68, slh_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_3 * pc_z[k] * slh_66[k];

        t_92[k] = f_12 * skh_26[k]
                  + f_3 * pc_y[k] * slh_68[k];

        t_93[k] = f_16 * skh_72[k]
                  + f_6 * slg0_54[k]
                  - f_7 * slg1_54[k]
                  + f_3 * pc_x[k] * slh_72[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, pc_x, pc_z, skh_73, skh_75, slg0_55, slg0_57, \
                         slg1_55, slg1_57, slh_69, slh_73, slh_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_16 * skh_73[k]
                  + f_8 * slg0_55[k]
                  - f_9 * slg1_55[k]
                  + f_3 * pc_x[k] * slh_73[k];

        t_95[k] = f_3 * pc_z[k] * slh_69[k];

        t_96[k] = f_16 * skh_75[k]
                  + f_8 * slg0_57[k]
                  - f_9 * slg1_57[k]
                  + f_3 * pc_x[k] * slh_75[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pc_x, pc_y, skh_30, skh_77, skh_78, skh_79, \
                         slg0_59, slg1_59, slh_72, slh_77, slh_78, \
                         slh_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_12 * skh_30[k]
                  + f_3 * pc_y[k] * slh_72[k];

        t_98[k] = f_16 * skh_77[k]
                  + f_8 * slg0_59[k]
                  - f_9 * slg1_59[k]
                  + f_3 * pc_x[k] * slh_77[k];

        t_99[k] = f_16 * skh_78[k]
                  + f_3 * pc_x[k] * slh_78[k];

        t_100[k] = f_16 * skh_79[k]
                   + f_3 * pc_x[k] * slh_79[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, pc_x, skh_80, skh_81, skh_82, skh_83, \
                         slh_80, slh_81, slh_82, slh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_16 * skh_80[k]
                   + f_3 * pc_x[k] * slh_80[k];

        t_102[k] = f_16 * skh_81[k]
                   + f_3 * pc_x[k] * slh_81[k];

        t_103[k] = f_16 * skh_82[k]
                   + f_3 * pc_x[k] * slh_82[k];

        t_104[k] = f_16 * skh_83[k]
                   + f_3 * pc_x[k] * slh_83[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, pc_y, pc_z, skh_36, skh_38, slg0_55, slg0_57, \
                         slg1_55, slg1_57, slh_78, slh_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_12 * skh_36[k]
                   + f_1 * slg0_55[k]
                   - f_2 * slg1_55[k]
                   + f_3 * pc_y[k] * slh_78[k];

        t_106[k] = f_3 * pc_z[k] * slh_78[k];

        t_107[k] = f_12 * skh_38[k]
                   + f_4 * slg0_57[k]
                   - f_5 * slg1_57[k]
                   + f_3 * pc_y[k] * slh_80[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, pc_y, pc_z, skh_39, skh_40, skh_41, \
                         slg0_58, slg0_59, slg1_58, slg1_59, slh_81, slh_82, \
                         slh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_12 * skh_39[k]
                   + f_6 * slg0_58[k]
                   - f_7 * slg1_58[k]
                   + f_3 * pc_y[k] * slh_81[k];

        t_109[k] = f_12 * skh_40[k]
                   + f_8 * slg0_59[k]
                   - f_9 * slg1_59[k]
                   + f_3 * pc_y[k] * slh_82[k];

        t_110[k] = f_12 * skh_41[k]
                   + f_3 * pc_y[k] * slh_83[k];

        t_111[k] = f_1 * slg0_59[k]
                   - f_2 * slg1_59[k]
                   + f_3 * pc_z[k] * slh_83[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, pb_y, pb_z, pc_y, pc_z, ski0_31, ski0_56, \
                         skh_21, skh_42, ski1_31, ski1_56, slh_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = pb_y[k] * ski0_56[k]
                   - f_10 * pc_y[k] * ski1_56[k];

        t_113[k] = f_11 * skh_42[k]
                   + f_3 * pc_y[k] * slh_84[k];

        t_114[k] = f_11 * skh_21[k]
                   + f_3 * pc_z[k] * slh_84[k];

        t_115[k] = pb_z[k] * ski0_31[k]
                   - f_10 * pc_z[k] * ski1_31[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, pb_y, pb_z, pc_y, pc_z, ski0_34, ski0_61, \
                         skh_24, skh_44, ski1_34, ski1_61, slh_86, \
                         slh_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_11 * skh_44[k]
                   + f_3 * pc_y[k] * slh_86[k];

        t_117[k] = pb_y[k] * ski0_61[k]
                   - f_10 * pc_y[k] * ski1_61[k];

        t_118[k] = pb_z[k] * ski0_34[k]
                   - f_10 * pc_z[k] * ski1_34[k];

        t_119[k] = f_11 * skh_24[k]
                   + f_3 * pc_z[k] * slh_87[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, pb_y, pb_z, pc_y, pc_z, ski0_38, ski0_65, \
                         skh_27, skh_47, ski1_38, ski1_65, slh_89, \
                         slh_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_11 * skh_47[k]
                   + f_3 * pc_y[k] * slh_89[k];

        t_121[k] = pb_y[k] * ski0_65[k]
                   - f_10 * pc_y[k] * ski1_65[k];

        t_122[k] = pb_z[k] * ski0_38[k]
                   - f_10 * pc_z[k] * ski1_38[k];

        t_123[k] = f_11 * skh_27[k]
                   + f_3 * pc_z[k] * slh_90[k];
    }
}

static auto
compute_prim_sli_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t ski0,
                                                          const size_t skh, const size_t ski1,
                                                          const size_t slg0, const size_t slg1,
                                                          const size_t slh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.5 / gamma;
    const auto f_5 = 1.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 0.5 / gamma;
    const auto f_9 = 0.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_16 = 3.0 / q;
    const auto f_17 = 2.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ski0_49 = buffer.data(ski0 + 49);
    const auto *ski0_68 = buffer.data(ski0 + 68);
    const auto *ski0_70 = buffer.data(ski0 + 70);
    const auto *ski0_83 = buffer.data(ski0 + 83);
    const auto *ski0_84 = buffer.data(ski0 + 84);
    const auto *ski0_87 = buffer.data(ski0 + 87);
    const auto *ski0_90 = buffer.data(ski0 + 90);
    const auto *ski0_94 = buffer.data(ski0 + 94);
    const auto *ski0_96 = buffer.data(ski0 + 96);
    const auto *ski0_105 = buffer.data(ski0 + 105);
    const auto *ski0_140 = buffer.data(ski0 + 140);
    const auto *ski0_143 = buffer.data(ski0 + 143);
    const auto *ski0_145 = buffer.data(ski0 + 145);
    const auto *ski0_146 = buffer.data(ski0 + 146);
    const auto *ski0_149 = buffer.data(ski0 + 149);
    const auto *ski0_150 = buffer.data(ski0 + 150);
    const auto *ski0_152 = buffer.data(ski0 + 152);
    const auto *ski0_154 = buffer.data(ski0 + 154);

    const auto *skh_36 = buffer.data(skh + 36);
    const auto *skh_42 = buffer.data(skh + 42);
    const auto *skh_45 = buffer.data(skh + 45);
    const auto *skh_48 = buffer.data(skh + 48);
    const auto *skh_50 = buffer.data(skh + 50);
    const auto *skh_51 = buffer.data(skh + 51);
    const auto *skh_57 = buffer.data(skh + 57);
    const auto *skh_59 = buffer.data(skh + 59);
    const auto *skh_60 = buffer.data(skh + 60);
    const auto *skh_61 = buffer.data(skh + 61);
    const auto *skh_62 = buffer.data(skh + 62);
    const auto *skh_63 = buffer.data(skh + 63);
    const auto *skh_65 = buffer.data(skh + 65);
    const auto *skh_66 = buffer.data(skh + 66);
    const auto *skh_68 = buffer.data(skh + 68);
    const auto *skh_69 = buffer.data(skh + 69);
    const auto *skh_70 = buffer.data(skh + 70);
    const auto *skh_72 = buffer.data(skh + 72);
    const auto *skh_78 = buffer.data(skh + 78);
    const auto *skh_80 = buffer.data(skh + 80);
    const auto *skh_81 = buffer.data(skh + 81);
    const auto *skh_82 = buffer.data(skh + 82);
    const auto *skh_83 = buffer.data(skh + 83);
    const auto *skh_84 = buffer.data(skh + 84);
    const auto *skh_86 = buffer.data(skh + 86);
    const auto *skh_87 = buffer.data(skh + 87);
    const auto *skh_89 = buffer.data(skh + 89);
    const auto *skh_90 = buffer.data(skh + 90);
    const auto *skh_93 = buffer.data(skh + 93);
    const auto *skh_99 = buffer.data(skh + 99);
    const auto *skh_100 = buffer.data(skh + 100);
    const auto *skh_101 = buffer.data(skh + 101);
    const auto *skh_102 = buffer.data(skh + 102);
    const auto *skh_103 = buffer.data(skh + 103);
    const auto *skh_104 = buffer.data(skh + 104);
    const auto *skh_105 = buffer.data(skh + 105);
    const auto *skh_106 = buffer.data(skh + 106);
    const auto *skh_107 = buffer.data(skh + 107);
    const auto *skh_108 = buffer.data(skh + 108);
    const auto *skh_110 = buffer.data(skh + 110);
    const auto *skh_111 = buffer.data(skh + 111);
    const auto *skh_113 = buffer.data(skh + 113);
    const auto *skh_114 = buffer.data(skh + 114);
    const auto *skh_115 = buffer.data(skh + 115);
    const auto *skh_117 = buffer.data(skh + 117);
    const auto *skh_119 = buffer.data(skh + 119);
    const auto *skh_120 = buffer.data(skh + 120);
    const auto *skh_121 = buffer.data(skh + 121);
    const auto *skh_122 = buffer.data(skh + 122);
    const auto *skh_123 = buffer.data(skh + 123);
    const auto *skh_124 = buffer.data(skh + 124);
    const auto *skh_125 = buffer.data(skh + 125);
    const auto *skh_126 = buffer.data(skh + 126);
    const auto *skh_129 = buffer.data(skh + 129);
    const auto *skh_131 = buffer.data(skh + 131);
    const auto *skh_132 = buffer.data(skh + 132);
    const auto *skh_135 = buffer.data(skh + 135);
    const auto *skh_136 = buffer.data(skh + 136);
    const auto *skh_138 = buffer.data(skh + 138);
    const auto *skh_140 = buffer.data(skh + 140);
    const auto *skh_141 = buffer.data(skh + 141);
    const auto *skh_142 = buffer.data(skh + 142);
    const auto *skh_143 = buffer.data(skh + 143);
    const auto *skh_144 = buffer.data(skh + 144);
    const auto *skh_145 = buffer.data(skh + 145);
    const auto *skh_146 = buffer.data(skh + 146);
    const auto *skh_152 = buffer.data(skh + 152);
    const auto *skh_156 = buffer.data(skh + 156);
    const auto *skh_161 = buffer.data(skh + 161);
    const auto *skh_162 = buffer.data(skh + 162);
    const auto *skh_163 = buffer.data(skh + 163);
    const auto *skh_164 = buffer.data(skh + 164);
    const auto *skh_165 = buffer.data(skh + 165);
    const auto *skh_166 = buffer.data(skh + 166);
    const auto *skh_167 = buffer.data(skh + 167);
    const auto *skh_183 = buffer.data(skh + 183);

    const auto *ski1_49 = buffer.data(ski1 + 49);
    const auto *ski1_68 = buffer.data(ski1 + 68);
    const auto *ski1_70 = buffer.data(ski1 + 70);
    const auto *ski1_83 = buffer.data(ski1 + 83);
    const auto *ski1_84 = buffer.data(ski1 + 84);
    const auto *ski1_87 = buffer.data(ski1 + 87);
    const auto *ski1_90 = buffer.data(ski1 + 90);
    const auto *ski1_94 = buffer.data(ski1 + 94);
    const auto *ski1_96 = buffer.data(ski1 + 96);
    const auto *ski1_105 = buffer.data(ski1 + 105);
    const auto *ski1_140 = buffer.data(ski1 + 140);
    const auto *ski1_143 = buffer.data(ski1 + 143);
    const auto *ski1_145 = buffer.data(ski1 + 145);
    const auto *ski1_146 = buffer.data(ski1 + 146);
    const auto *ski1_149 = buffer.data(ski1 + 149);
    const auto *ski1_150 = buffer.data(ski1 + 150);
    const auto *ski1_152 = buffer.data(ski1 + 152);
    const auto *ski1_154 = buffer.data(ski1 + 154);

    const auto *slg0_72 = buffer.data(slg0 + 72);
    const auto *slg0_73 = buffer.data(slg0 + 73);
    const auto *slg0_74 = buffer.data(slg0 + 74);
    const auto *slg0_75 = buffer.data(slg0 + 75);
    const auto *slg0_78 = buffer.data(slg0 + 78);
    const auto *slg0_80 = buffer.data(slg0 + 80);
    const auto *slg0_81 = buffer.data(slg0 + 81);
    const auto *slg0_84 = buffer.data(slg0 + 84);
    const auto *slg0_85 = buffer.data(slg0 + 85);
    const auto *slg0_87 = buffer.data(slg0 + 87);
    const auto *slg0_88 = buffer.data(slg0 + 88);
    const auto *slg0_89 = buffer.data(slg0 + 89);
    const auto *slg0_90 = buffer.data(slg0 + 90);
    const auto *slg0_93 = buffer.data(slg0 + 93);
    const auto *slg0_95 = buffer.data(slg0 + 95);
    const auto *slg0_96 = buffer.data(slg0 + 96);
    const auto *slg0_99 = buffer.data(slg0 + 99);
    const auto *slg0_100 = buffer.data(slg0 + 100);
    const auto *slg0_102 = buffer.data(slg0 + 102);
    const auto *slg0_103 = buffer.data(slg0 + 103);
    const auto *slg0_104 = buffer.data(slg0 + 104);
    const auto *slg0_110 = buffer.data(slg0 + 110);
    const auto *slg0_114 = buffer.data(slg0 + 114);
    const auto *slg0_117 = buffer.data(slg0 + 117);
    const auto *slg0_118 = buffer.data(slg0 + 118);
    const auto *slg0_119 = buffer.data(slg0 + 119);

    const auto *slg1_72 = buffer.data(slg1 + 72);
    const auto *slg1_73 = buffer.data(slg1 + 73);
    const auto *slg1_74 = buffer.data(slg1 + 74);
    const auto *slg1_75 = buffer.data(slg1 + 75);
    const auto *slg1_78 = buffer.data(slg1 + 78);
    const auto *slg1_80 = buffer.data(slg1 + 80);
    const auto *slg1_81 = buffer.data(slg1 + 81);
    const auto *slg1_84 = buffer.data(slg1 + 84);
    const auto *slg1_85 = buffer.data(slg1 + 85);
    const auto *slg1_87 = buffer.data(slg1 + 87);
    const auto *slg1_88 = buffer.data(slg1 + 88);
    const auto *slg1_89 = buffer.data(slg1 + 89);
    const auto *slg1_90 = buffer.data(slg1 + 90);
    const auto *slg1_93 = buffer.data(slg1 + 93);
    const auto *slg1_95 = buffer.data(slg1 + 95);
    const auto *slg1_96 = buffer.data(slg1 + 96);
    const auto *slg1_99 = buffer.data(slg1 + 99);
    const auto *slg1_100 = buffer.data(slg1 + 100);
    const auto *slg1_102 = buffer.data(slg1 + 102);
    const auto *slg1_103 = buffer.data(slg1 + 103);
    const auto *slg1_104 = buffer.data(slg1 + 104);
    const auto *slg1_110 = buffer.data(slg1 + 110);
    const auto *slg1_114 = buffer.data(slg1 + 114);
    const auto *slg1_117 = buffer.data(slg1 + 117);
    const auto *slg1_118 = buffer.data(slg1 + 118);
    const auto *slg1_119 = buffer.data(slg1 + 119);

    const auto *slh_93 = buffer.data(slh + 93);
    const auto *slh_99 = buffer.data(slh + 99);
    const auto *slh_100 = buffer.data(slh + 100);
    const auto *slh_101 = buffer.data(slh + 101);
    const auto *slh_102 = buffer.data(slh + 102);
    const auto *slh_103 = buffer.data(slh + 103);
    const auto *slh_104 = buffer.data(slh + 104);
    const auto *slh_105 = buffer.data(slh + 105);
    const auto *slh_107 = buffer.data(slh + 107);
    const auto *slh_108 = buffer.data(slh + 108);
    const auto *slh_110 = buffer.data(slh + 110);
    const auto *slh_111 = buffer.data(slh + 111);
    const auto *slh_114 = buffer.data(slh + 114);
    const auto *slh_115 = buffer.data(slh + 115);
    const auto *slh_117 = buffer.data(slh + 117);
    const auto *slh_119 = buffer.data(slh + 119);
    const auto *slh_120 = buffer.data(slh + 120);
    const auto *slh_121 = buffer.data(slh + 121);
    const auto *slh_122 = buffer.data(slh + 122);
    const auto *slh_123 = buffer.data(slh + 123);
    const auto *slh_124 = buffer.data(slh + 124);
    const auto *slh_125 = buffer.data(slh + 125);
    const auto *slh_126 = buffer.data(slh + 126);
    const auto *slh_128 = buffer.data(slh + 128);
    const auto *slh_129 = buffer.data(slh + 129);
    const auto *slh_131 = buffer.data(slh + 131);
    const auto *slh_132 = buffer.data(slh + 132);
    const auto *slh_135 = buffer.data(slh + 135);
    const auto *slh_136 = buffer.data(slh + 136);
    const auto *slh_138 = buffer.data(slh + 138);
    const auto *slh_140 = buffer.data(slh + 140);
    const auto *slh_141 = buffer.data(slh + 141);
    const auto *slh_142 = buffer.data(slh + 142);
    const auto *slh_143 = buffer.data(slh + 143);
    const auto *slh_144 = buffer.data(slh + 144);
    const auto *slh_145 = buffer.data(slh + 145);
    const auto *slh_146 = buffer.data(slh + 146);
    const auto *slh_147 = buffer.data(slh + 147);
    const auto *slh_149 = buffer.data(slh + 149);
    const auto *slh_150 = buffer.data(slh + 150);
    const auto *slh_152 = buffer.data(slh + 152);
    const auto *slh_153 = buffer.data(slh + 153);
    const auto *slh_156 = buffer.data(slh + 156);
    const auto *slh_161 = buffer.data(slh + 161);
    const auto *slh_162 = buffer.data(slh + 162);
    const auto *slh_163 = buffer.data(slh + 163);
    const auto *slh_164 = buffer.data(slh + 164);
    const auto *slh_165 = buffer.data(slh + 165);
    const auto *slh_166 = buffer.data(slh + 166);
    const auto *slh_167 = buffer.data(slh + 167);
    const auto *slh_168 = buffer.data(slh + 168);
    const auto *slh_170 = buffer.data(slh + 170);
    const auto *slh_171 = buffer.data(slh + 171);
    const auto *slh_173 = buffer.data(slh + 173);
    const auto *slh_174 = buffer.data(slh + 174);
    const auto *slh_177 = buffer.data(slh + 177);
    const auto *slh_183 = buffer.data(slh + 183);

#pragma omp simd aligned(t_124, t_125, t_126, t_127, pb_y, pc_x, pc_y, ski0_68, ski0_70, \
                         skh_50, skh_51, skh_99, ski1_68, ski1_70, slh_93, \
                         slh_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = pb_y[k] * ski0_68[k]
                   + f_12 * skh_50[k]
                   - f_10 * pc_y[k] * ski1_68[k];

        t_125[k] = f_11 * skh_51[k]
                   + f_3 * pc_y[k] * slh_93[k];

        t_126[k] = pb_y[k] * ski0_70[k]
                   - f_10 * pc_y[k] * ski1_70[k];

        t_127[k] = f_16 * skh_99[k]
                   + f_3 * pc_x[k] * slh_99[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, t_132, pc_x, skh_100, skh_101, skh_102, \
                         skh_103, skh_104, slh_100, slh_101, slh_102, slh_103, \
                         slh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_16 * skh_100[k]
                   + f_3 * pc_x[k] * slh_100[k];

        t_129[k] = f_16 * skh_101[k]
                   + f_3 * pc_x[k] * slh_101[k];

        t_130[k] = f_16 * skh_102[k]
                   + f_3 * pc_x[k] * slh_102[k];

        t_131[k] = f_16 * skh_103[k]
                   + f_3 * pc_x[k] * slh_103[k];

        t_132[k] = f_16 * skh_104[k]
                   + f_3 * pc_x[k] * slh_104[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, pb_z, pc_y, pc_z, ski0_49, skh_36, skh_59, \
                         ski1_49, slg0_72, slg1_72, slh_99, slh_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = pb_z[k] * ski0_49[k]
                   - f_10 * pc_z[k] * ski1_49[k];

        t_134[k] = f_11 * skh_36[k]
                   + f_3 * pc_z[k] * slh_99[k];

        t_135[k] = f_11 * skh_59[k]
                   + f_4 * slg0_72[k]
                   - f_5 * slg1_72[k]
                   + f_3 * pc_y[k] * slh_101[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, pc_y, skh_60, skh_61, skh_62, slg0_73, slg0_74, \
                         slg1_73, slg1_74, slh_102, slh_103, slh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_11 * skh_60[k]
                   + f_6 * slg0_73[k]
                   - f_7 * slg1_73[k]
                   + f_3 * pc_y[k] * slh_102[k];

        t_137[k] = f_11 * skh_61[k]
                   + f_8 * slg0_74[k]
                   - f_9 * slg1_74[k]
                   + f_3 * pc_y[k] * slh_103[k];

        t_138[k] = f_11 * skh_62[k]
                   + f_3 * pc_y[k] * slh_104[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, pb_y, pc_x, pc_y, pc_z, ski0_83, skh_42, \
                         skh_105, ski1_83, slg0_75, slg1_75, slh_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = pb_y[k] * ski0_83[k]
                   - f_10 * pc_y[k] * ski1_83[k];

        t_140[k] = f_16 * skh_105[k]
                   + f_1 * slg0_75[k]
                   - f_2 * slg1_75[k]
                   + f_3 * pc_x[k] * slh_105[k];

        t_141[k] = f_3 * pc_y[k] * slh_105[k];

        t_142[k] = f_12 * skh_42[k]
                   + f_3 * pc_z[k] * slh_105[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, pc_x, pc_y, skh_108, skh_110, slg0_78, slg0_80, \
                         slg1_78, slg1_80, slh_107, slh_108, slh_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_16 * skh_108[k]
                   + f_4 * slg0_78[k]
                   - f_5 * slg1_78[k]
                   + f_3 * pc_x[k] * slh_108[k];

        t_144[k] = f_3 * pc_y[k] * slh_107[k];

        t_145[k] = f_16 * skh_110[k]
                   + f_4 * slg0_80[k]
                   - f_5 * slg1_80[k]
                   + f_3 * pc_x[k] * slh_110[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, pc_x, pc_y, pc_z, skh_45, skh_111, slg0_81, \
                         slg1_81, slh_108, slh_110, slh_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_16 * skh_111[k]
                   + f_6 * slg0_81[k]
                   - f_7 * slg1_81[k]
                   + f_3 * pc_x[k] * slh_111[k];

        t_147[k] = f_12 * skh_45[k]
                   + f_3 * pc_z[k] * slh_108[k];

        t_148[k] = f_3 * pc_y[k] * slh_110[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, pc_x, pc_z, skh_48, skh_114, skh_115, slg0_84, \
                         slg0_85, slg1_84, slg1_85, slh_111, slh_114, \
                         slh_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_16 * skh_114[k]
                   + f_6 * slg0_84[k]
                   - f_7 * slg1_84[k]
                   + f_3 * pc_x[k] * slh_114[k];

        t_150[k] = f_16 * skh_115[k]
                   + f_8 * slg0_85[k]
                   - f_9 * slg1_85[k]
                   + f_3 * pc_x[k] * slh_115[k];

        t_151[k] = f_12 * skh_48[k]
                   + f_3 * pc_z[k] * slh_111[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, pc_x, pc_y, skh_117, skh_119, slg0_87, slg0_89, \
                         slg1_87, slg1_89, slh_114, slh_117, slh_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_16 * skh_117[k]
                   + f_8 * slg0_87[k]
                   - f_9 * slg1_87[k]
                   + f_3 * pc_x[k] * slh_117[k];

        t_153[k] = f_3 * pc_y[k] * slh_114[k];

        t_154[k] = f_16 * skh_119[k]
                   + f_8 * slg0_89[k]
                   - f_9 * slg1_89[k]
                   + f_3 * pc_x[k] * slh_119[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, pc_x, skh_120, skh_121, skh_122, \
                         skh_123, skh_124, slh_120, slh_121, slh_122, slh_123, \
                         slh_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_16 * skh_120[k]
                   + f_3 * pc_x[k] * slh_120[k];

        t_156[k] = f_16 * skh_121[k]
                   + f_3 * pc_x[k] * slh_121[k];

        t_157[k] = f_16 * skh_122[k]
                   + f_3 * pc_x[k] * slh_122[k];

        t_158[k] = f_16 * skh_123[k]
                   + f_3 * pc_x[k] * slh_123[k];

        t_159[k] = f_16 * skh_124[k]
                   + f_3 * pc_x[k] * slh_124[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, pc_x, pc_y, pc_z, skh_57, skh_125, \
                         slg0_85, slg0_87, slg1_85, slg1_87, slh_120, slh_122, \
                         slh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = f_16 * skh_125[k]
                   + f_3 * pc_x[k] * slh_125[k];

        t_161[k] = f_1 * slg0_85[k]
                   - f_2 * slg1_85[k]
                   + f_3 * pc_y[k] * slh_120[k];

        t_162[k] = f_12 * skh_57[k]
                   + f_3 * pc_z[k] * slh_120[k];

        t_163[k] = f_4 * slg0_87[k]
                   - f_5 * slg1_87[k]
                   + f_3 * pc_y[k] * slh_122[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, pc_y, pc_z, skh_62, slg0_88, slg0_89, \
                         slg1_88, slg1_89, slh_123, slh_124, slh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_6 * slg0_88[k]
                   - f_7 * slg1_88[k]
                   + f_3 * pc_y[k] * slh_123[k];

        t_165[k] = f_8 * slg0_89[k]
                   - f_9 * slg1_89[k]
                   + f_3 * pc_y[k] * slh_124[k];

        t_166[k] = f_3 * pc_y[k] * slh_125[k];

        t_167[k] = f_12 * skh_62[k]
                   + f_1 * slg0_89[k]
                   - f_2 * slg1_89[k]
                   + f_3 * pc_z[k] * slh_125[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, pc_x, pc_y, pc_z, skh_63, skh_126, \
                         skh_129, slg0_90, slg0_93, slg1_90, slg1_93, slh_126, \
                         slh_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_17 * skh_126[k]
                   + f_1 * slg0_90[k]
                   - f_2 * slg1_90[k]
                   + f_3 * pc_x[k] * slh_126[k];

        t_169[k] = f_13 * skh_63[k]
                   + f_3 * pc_y[k] * slh_126[k];

        t_170[k] = f_3 * pc_z[k] * slh_126[k];

        t_171[k] = f_17 * skh_129[k]
                   + f_4 * slg0_93[k]
                   - f_5 * slg1_93[k]
                   + f_3 * pc_x[k] * slh_129[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, pc_x, pc_y, skh_65, skh_131, skh_132, slg0_95, \
                         slg0_96, slg1_95, slg1_96, slh_128, slh_131, \
                         slh_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_13 * skh_65[k]
                   + f_3 * pc_y[k] * slh_128[k];

        t_173[k] = f_17 * skh_131[k]
                   + f_4 * slg0_95[k]
                   - f_5 * slg1_95[k]
                   + f_3 * pc_x[k] * slh_131[k];

        t_174[k] = f_17 * skh_132[k]
                   + f_6 * slg0_96[k]
                   - f_7 * slg1_96[k]
                   + f_3 * pc_x[k] * slh_132[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, pc_x, pc_y, pc_z, skh_68, skh_135, slg0_99, \
                         slg1_99, slh_129, slh_131, slh_135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = f_3 * pc_z[k] * slh_129[k];

        t_176[k] = f_13 * skh_68[k]
                   + f_3 * pc_y[k] * slh_131[k];

        t_177[k] = f_17 * skh_135[k]
                   + f_6 * slg0_99[k]
                   - f_7 * slg1_99[k]
                   + f_3 * pc_x[k] * slh_135[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, pc_x, pc_z, skh_136, skh_138, slg0_100, \
                         slg0_102, slg1_100, slg1_102, slh_132, slh_136, \
                         slh_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_17 * skh_136[k]
                   + f_8 * slg0_100[k]
                   - f_9 * slg1_100[k]
                   + f_3 * pc_x[k] * slh_136[k];

        t_179[k] = f_3 * pc_z[k] * slh_132[k];

        t_180[k] = f_17 * skh_138[k]
                   + f_8 * slg0_102[k]
                   - f_9 * slg1_102[k]
                   + f_3 * pc_x[k] * slh_138[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pc_x, pc_y, skh_72, skh_140, skh_141, \
                         skh_142, slg0_104, slg1_104, slh_135, slh_140, slh_141, \
                         slh_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_13 * skh_72[k]
                   + f_3 * pc_y[k] * slh_135[k];

        t_182[k] = f_17 * skh_140[k]
                   + f_8 * slg0_104[k]
                   - f_9 * slg1_104[k]
                   + f_3 * pc_x[k] * slh_140[k];

        t_183[k] = f_17 * skh_141[k]
                   + f_3 * pc_x[k] * slh_141[k];

        t_184[k] = f_17 * skh_142[k]
                   + f_3 * pc_x[k] * slh_142[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pc_x, skh_143, skh_144, skh_145, skh_146, \
                         slh_143, slh_144, slh_145, slh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_17 * skh_143[k]
                   + f_3 * pc_x[k] * slh_143[k];

        t_186[k] = f_17 * skh_144[k]
                   + f_3 * pc_x[k] * slh_144[k];

        t_187[k] = f_17 * skh_145[k]
                   + f_3 * pc_x[k] * slh_145[k];

        t_188[k] = f_17 * skh_146[k]
                   + f_3 * pc_x[k] * slh_146[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, pc_y, pc_z, skh_78, skh_80, slg0_100, slg0_102, \
                         slg1_100, slg1_102, slh_141, slh_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_13 * skh_78[k]
                   + f_1 * slg0_100[k]
                   - f_2 * slg1_100[k]
                   + f_3 * pc_y[k] * slh_141[k];

        t_190[k] = f_3 * pc_z[k] * slh_141[k];

        t_191[k] = f_13 * skh_80[k]
                   + f_4 * slg0_102[k]
                   - f_5 * slg1_102[k]
                   + f_3 * pc_y[k] * slh_143[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, pc_y, pc_z, skh_81, skh_82, skh_83, \
                         slg0_103, slg0_104, slg1_103, slg1_104, slh_144, slh_145, \
                         slh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_13 * skh_81[k]
                   + f_6 * slg0_103[k]
                   - f_7 * slg1_103[k]
                   + f_3 * pc_y[k] * slh_144[k];

        t_193[k] = f_13 * skh_82[k]
                   + f_8 * slg0_104[k]
                   - f_9 * slg1_104[k]
                   + f_3 * pc_y[k] * slh_145[k];

        t_194[k] = f_13 * skh_83[k]
                   + f_3 * pc_y[k] * slh_146[k];

        t_195[k] = f_1 * slg0_104[k]
                   - f_2 * slg1_104[k]
                   + f_3 * pc_z[k] * slh_146[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, pb_z, pc_y, pc_z, ski0_84, ski0_87, \
                         skh_63, skh_84, ski1_84, ski1_87, slh_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = pb_z[k] * ski0_84[k]
                   - f_10 * pc_z[k] * ski1_84[k];

        t_197[k] = f_12 * skh_84[k]
                   + f_3 * pc_y[k] * slh_147[k];

        t_198[k] = f_11 * skh_63[k]
                   + f_3 * pc_z[k] * slh_147[k];

        t_199[k] = pb_z[k] * ski0_87[k]
                   - f_10 * pc_z[k] * ski1_87[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, pb_z, pc_x, pc_y, pc_z, ski0_90, skh_86, \
                         skh_152, ski1_90, slg0_110, slg1_110, slh_149, \
                         slh_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_12 * skh_86[k]
                   + f_3 * pc_y[k] * slh_149[k];

        t_201[k] = f_17 * skh_152[k]
                   + f_4 * slg0_110[k]
                   - f_5 * slg1_110[k]
                   + f_3 * pc_x[k] * slh_152[k];

        t_202[k] = pb_z[k] * ski0_90[k]
                   - f_10 * pc_z[k] * ski1_90[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, pc_x, pc_y, pc_z, skh_66, skh_89, skh_156, \
                         slg0_114, slg1_114, slh_150, slh_152, \
                         slh_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_11 * skh_66[k]
                   + f_3 * pc_z[k] * slh_150[k];

        t_204[k] = f_12 * skh_89[k]
                   + f_3 * pc_y[k] * slh_152[k];

        t_205[k] = f_17 * skh_156[k]
                   + f_6 * slg0_114[k]
                   - f_7 * slg1_114[k]
                   + f_3 * pc_x[k] * slh_156[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, pb_z, pc_y, pc_z, ski0_94, ski0_96, \
                         skh_69, skh_70, skh_93, ski1_94, ski1_96, slh_153, \
                         slh_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = pb_z[k] * ski0_94[k]
                   - f_10 * pc_z[k] * ski1_94[k];

        t_207[k] = f_11 * skh_69[k]
                   + f_3 * pc_z[k] * slh_153[k];

        t_208[k] = pb_z[k] * ski0_96[k]
                   + f_12 * skh_70[k]
                   - f_10 * pc_z[k] * ski1_96[k];

        t_209[k] = f_12 * skh_93[k]
                   + f_3 * pc_y[k] * slh_156[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, pc_x, skh_161, skh_162, skh_163, skh_164, \
                         slg0_119, slg1_119, slh_161, slh_162, slh_163, \
                         slh_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_17 * skh_161[k]
                   + f_8 * slg0_119[k]
                   - f_9 * slg1_119[k]
                   + f_3 * pc_x[k] * slh_161[k];

        t_211[k] = f_17 * skh_162[k]
                   + f_3 * pc_x[k] * slh_162[k];

        t_212[k] = f_17 * skh_163[k]
                   + f_3 * pc_x[k] * slh_163[k];

        t_213[k] = f_17 * skh_164[k]
                   + f_3 * pc_x[k] * slh_164[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, pb_z, pc_x, pc_z, ski0_105, skh_165, \
                         skh_166, skh_167, ski1_105, slh_165, slh_166, \
                         slh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = f_17 * skh_165[k]
                   + f_3 * pc_x[k] * slh_165[k];

        t_215[k] = f_17 * skh_166[k]
                   + f_3 * pc_x[k] * slh_166[k];

        t_216[k] = f_17 * skh_167[k]
                   + f_3 * pc_x[k] * slh_167[k];

        t_217[k] = pb_z[k] * ski0_105[k]
                   - f_10 * pc_z[k] * ski1_105[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, pc_y, pc_z, skh_78, skh_101, skh_102, slg0_117, \
                         slg0_118, slg1_117, slg1_118, slh_162, slh_164, \
                         slh_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_11 * skh_78[k]
                   + f_3 * pc_z[k] * slh_162[k];

        t_219[k] = f_12 * skh_101[k]
                   + f_4 * slg0_117[k]
                   - f_5 * slg1_117[k]
                   + f_3 * pc_y[k] * slh_164[k];

        t_220[k] = f_12 * skh_102[k]
                   + f_6 * slg0_118[k]
                   - f_7 * slg1_118[k]
                   + f_3 * pc_y[k] * slh_165[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, pb_y, pc_y, pc_z, ski0_140, skh_83, \
                         skh_103, skh_104, ski1_140, slg0_119, slg1_119, slh_166, \
                         slh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_12 * skh_103[k]
                   + f_8 * slg0_119[k]
                   - f_9 * slg1_119[k]
                   + f_3 * pc_y[k] * slh_166[k];

        t_222[k] = f_12 * skh_104[k]
                   + f_3 * pc_y[k] * slh_167[k];

        t_223[k] = f_11 * skh_83[k]
                   + f_1 * slg0_119[k]
                   - f_2 * slg1_119[k]
                   + f_3 * pc_z[k] * slh_167[k];

        t_224[k] = pb_y[k] * ski0_140[k]
                   - f_10 * pc_y[k] * ski1_140[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, pb_y, pc_y, pc_z, ski0_143, skh_84, \
                         skh_105, skh_106, skh_107, ski1_143, slh_168, \
                         slh_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_11 * skh_105[k]
                   + f_3 * pc_y[k] * slh_168[k];

        t_226[k] = f_12 * skh_84[k]
                   + f_3 * pc_z[k] * slh_168[k];

        t_227[k] = pb_y[k] * ski0_143[k]
                   + f_12 * skh_106[k]
                   - f_10 * pc_y[k] * ski1_143[k];

        t_228[k] = f_11 * skh_107[k]
                   + f_3 * pc_y[k] * slh_170[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, pb_y, pc_y, pc_z, ski0_145, ski0_146, \
                         skh_87, skh_108, skh_110, ski1_145, ski1_146, slh_171, \
                         slh_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = pb_y[k] * ski0_145[k]
                   - f_10 * pc_y[k] * ski1_145[k];

        t_230[k] = pb_y[k] * ski0_146[k]
                   + f_13 * skh_108[k]
                   - f_10 * pc_y[k] * ski1_146[k];

        t_231[k] = f_12 * skh_87[k]
                   + f_3 * pc_z[k] * slh_171[k];

        t_232[k] = f_11 * skh_110[k]
                   + f_3 * pc_y[k] * slh_173[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, pb_y, pc_y, pc_z, ski0_149, ski0_150, skh_90, \
                         skh_111, ski1_149, ski1_150, slh_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = pb_y[k] * ski0_149[k]
                   - f_10 * pc_y[k] * ski1_149[k];

        t_234[k] = pb_y[k] * ski0_150[k]
                   + f_14 * skh_111[k]
                   - f_10 * pc_y[k] * ski1_150[k];

        t_235[k] = f_12 * skh_90[k]
                   + f_3 * pc_z[k] * slh_174[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, t_239, pb_y, pc_x, pc_y, ski0_152, ski0_154, \
                         skh_113, skh_114, skh_183, ski1_152, ski1_154, slh_177, \
                         slh_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = pb_y[k] * ski0_152[k]
                   + f_12 * skh_113[k]
                   - f_10 * pc_y[k] * ski1_152[k];

        t_237[k] = f_11 * skh_114[k]
                   + f_3 * pc_y[k] * slh_177[k];

        t_238[k] = pb_y[k] * ski0_154[k]
                   - f_10 * pc_y[k] * ski1_154[k];

        t_239[k] = f_17 * skh_183[k]
                   + f_3 * pc_x[k] * slh_183[k];
    }
}

static auto
compute_prim_sli_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t ski0,
                                                          const size_t skh, const size_t ski1,
                                                          const size_t slg0, const size_t slg1,
                                                          const size_t slh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.5 / gamma;
    const auto f_5 = 1.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 0.5 / gamma;
    const auto f_9 = 0.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_17 = 2.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ski0_167 = buffer.data(ski0 + 167);
    const auto *ski0_168 = buffer.data(ski0 + 168);
    const auto *ski0_171 = buffer.data(ski0 + 171);
    const auto *ski0_174 = buffer.data(ski0 + 174);
    const auto *ski0_178 = buffer.data(ski0 + 178);
    const auto *ski0_180 = buffer.data(ski0 + 180);
    const auto *ski0_189 = buffer.data(ski0 + 189);

    const auto *skh_99 = buffer.data(skh + 99);
    const auto *skh_105 = buffer.data(skh + 105);
    const auto *skh_108 = buffer.data(skh + 108);
    const auto *skh_111 = buffer.data(skh + 111);
    const auto *skh_120 = buffer.data(skh + 120);
    const auto *skh_122 = buffer.data(skh + 122);
    const auto *skh_123 = buffer.data(skh + 123);
    const auto *skh_124 = buffer.data(skh + 124);
    const auto *skh_125 = buffer.data(skh + 125);
    const auto *skh_126 = buffer.data(skh + 126);
    const auto *skh_128 = buffer.data(skh + 128);
    const auto *skh_129 = buffer.data(skh + 129);
    const auto *skh_131 = buffer.data(skh + 131);
    const auto *skh_132 = buffer.data(skh + 132);
    const auto *skh_133 = buffer.data(skh + 133);
    const auto *skh_135 = buffer.data(skh + 135);
    const auto *skh_141 = buffer.data(skh + 141);
    const auto *skh_143 = buffer.data(skh + 143);
    const auto *skh_144 = buffer.data(skh + 144);
    const auto *skh_145 = buffer.data(skh + 145);
    const auto *skh_146 = buffer.data(skh + 146);
    const auto *skh_147 = buffer.data(skh + 147);
    const auto *skh_149 = buffer.data(skh + 149);
    const auto *skh_150 = buffer.data(skh + 150);
    const auto *skh_152 = buffer.data(skh + 152);
    const auto *skh_153 = buffer.data(skh + 153);
    const auto *skh_156 = buffer.data(skh + 156);
    const auto *skh_164 = buffer.data(skh + 164);
    const auto *skh_165 = buffer.data(skh + 165);
    const auto *skh_166 = buffer.data(skh + 166);
    const auto *skh_167 = buffer.data(skh + 167);
    const auto *skh_168 = buffer.data(skh + 168);
    const auto *skh_170 = buffer.data(skh + 170);
    const auto *skh_173 = buffer.data(skh + 173);
    const auto *skh_177 = buffer.data(skh + 177);
    const auto *skh_184 = buffer.data(skh + 184);
    const auto *skh_185 = buffer.data(skh + 185);
    const auto *skh_186 = buffer.data(skh + 186);
    const auto *skh_187 = buffer.data(skh + 187);
    const auto *skh_188 = buffer.data(skh + 188);
    const auto *skh_189 = buffer.data(skh + 189);
    const auto *skh_192 = buffer.data(skh + 192);
    const auto *skh_194 = buffer.data(skh + 194);
    const auto *skh_195 = buffer.data(skh + 195);
    const auto *skh_198 = buffer.data(skh + 198);
    const auto *skh_199 = buffer.data(skh + 199);
    const auto *skh_201 = buffer.data(skh + 201);
    const auto *skh_203 = buffer.data(skh + 203);
    const auto *skh_204 = buffer.data(skh + 204);
    const auto *skh_205 = buffer.data(skh + 205);
    const auto *skh_206 = buffer.data(skh + 206);
    const auto *skh_207 = buffer.data(skh + 207);
    const auto *skh_208 = buffer.data(skh + 208);
    const auto *skh_209 = buffer.data(skh + 209);
    const auto *skh_210 = buffer.data(skh + 210);
    const auto *skh_213 = buffer.data(skh + 213);
    const auto *skh_215 = buffer.data(skh + 215);
    const auto *skh_216 = buffer.data(skh + 216);
    const auto *skh_219 = buffer.data(skh + 219);
    const auto *skh_220 = buffer.data(skh + 220);
    const auto *skh_222 = buffer.data(skh + 222);
    const auto *skh_224 = buffer.data(skh + 224);
    const auto *skh_225 = buffer.data(skh + 225);
    const auto *skh_226 = buffer.data(skh + 226);
    const auto *skh_227 = buffer.data(skh + 227);
    const auto *skh_228 = buffer.data(skh + 228);
    const auto *skh_229 = buffer.data(skh + 229);
    const auto *skh_230 = buffer.data(skh + 230);
    const auto *skh_236 = buffer.data(skh + 236);
    const auto *skh_240 = buffer.data(skh + 240);
    const auto *skh_245 = buffer.data(skh + 245);
    const auto *skh_246 = buffer.data(skh + 246);
    const auto *skh_247 = buffer.data(skh + 247);
    const auto *skh_248 = buffer.data(skh + 248);
    const auto *skh_249 = buffer.data(skh + 249);
    const auto *skh_250 = buffer.data(skh + 250);
    const auto *skh_251 = buffer.data(skh + 251);
    const auto *skh_252 = buffer.data(skh + 252);
    const auto *skh_255 = buffer.data(skh + 255);
    const auto *skh_257 = buffer.data(skh + 257);
    const auto *skh_258 = buffer.data(skh + 258);
    const auto *skh_261 = buffer.data(skh + 261);
    const auto *skh_262 = buffer.data(skh + 262);
    const auto *skh_264 = buffer.data(skh + 264);
    const auto *skh_266 = buffer.data(skh + 266);

    const auto *ski1_167 = buffer.data(ski1 + 167);
    const auto *ski1_168 = buffer.data(ski1 + 168);
    const auto *ski1_171 = buffer.data(ski1 + 171);
    const auto *ski1_174 = buffer.data(ski1 + 174);
    const auto *ski1_178 = buffer.data(ski1 + 178);
    const auto *ski1_180 = buffer.data(ski1 + 180);
    const auto *ski1_189 = buffer.data(ski1 + 189);

    const auto *slg0_130 = buffer.data(slg0 + 130);
    const auto *slg0_132 = buffer.data(slg0 + 132);
    const auto *slg0_133 = buffer.data(slg0 + 133);
    const auto *slg0_134 = buffer.data(slg0 + 134);
    const auto *slg0_135 = buffer.data(slg0 + 135);
    const auto *slg0_138 = buffer.data(slg0 + 138);
    const auto *slg0_140 = buffer.data(slg0 + 140);
    const auto *slg0_141 = buffer.data(slg0 + 141);
    const auto *slg0_144 = buffer.data(slg0 + 144);
    const auto *slg0_145 = buffer.data(slg0 + 145);
    const auto *slg0_147 = buffer.data(slg0 + 147);
    const auto *slg0_148 = buffer.data(slg0 + 148);
    const auto *slg0_149 = buffer.data(slg0 + 149);
    const auto *slg0_150 = buffer.data(slg0 + 150);
    const auto *slg0_153 = buffer.data(slg0 + 153);
    const auto *slg0_155 = buffer.data(slg0 + 155);
    const auto *slg0_156 = buffer.data(slg0 + 156);
    const auto *slg0_159 = buffer.data(slg0 + 159);
    const auto *slg0_160 = buffer.data(slg0 + 160);
    const auto *slg0_162 = buffer.data(slg0 + 162);
    const auto *slg0_163 = buffer.data(slg0 + 163);
    const auto *slg0_164 = buffer.data(slg0 + 164);
    const auto *slg0_170 = buffer.data(slg0 + 170);
    const auto *slg0_174 = buffer.data(slg0 + 174);
    const auto *slg0_177 = buffer.data(slg0 + 177);
    const auto *slg0_178 = buffer.data(slg0 + 178);
    const auto *slg0_179 = buffer.data(slg0 + 179);
    const auto *slg0_180 = buffer.data(slg0 + 180);
    const auto *slg0_183 = buffer.data(slg0 + 183);
    const auto *slg0_185 = buffer.data(slg0 + 185);
    const auto *slg0_186 = buffer.data(slg0 + 186);
    const auto *slg0_189 = buffer.data(slg0 + 189);
    const auto *slg0_190 = buffer.data(slg0 + 190);
    const auto *slg0_192 = buffer.data(slg0 + 192);
    const auto *slg0_194 = buffer.data(slg0 + 194);

    const auto *slg1_130 = buffer.data(slg1 + 130);
    const auto *slg1_132 = buffer.data(slg1 + 132);
    const auto *slg1_133 = buffer.data(slg1 + 133);
    const auto *slg1_134 = buffer.data(slg1 + 134);
    const auto *slg1_135 = buffer.data(slg1 + 135);
    const auto *slg1_138 = buffer.data(slg1 + 138);
    const auto *slg1_140 = buffer.data(slg1 + 140);
    const auto *slg1_141 = buffer.data(slg1 + 141);
    const auto *slg1_144 = buffer.data(slg1 + 144);
    const auto *slg1_145 = buffer.data(slg1 + 145);
    const auto *slg1_147 = buffer.data(slg1 + 147);
    const auto *slg1_148 = buffer.data(slg1 + 148);
    const auto *slg1_149 = buffer.data(slg1 + 149);
    const auto *slg1_150 = buffer.data(slg1 + 150);
    const auto *slg1_153 = buffer.data(slg1 + 153);
    const auto *slg1_155 = buffer.data(slg1 + 155);
    const auto *slg1_156 = buffer.data(slg1 + 156);
    const auto *slg1_159 = buffer.data(slg1 + 159);
    const auto *slg1_160 = buffer.data(slg1 + 160);
    const auto *slg1_162 = buffer.data(slg1 + 162);
    const auto *slg1_163 = buffer.data(slg1 + 163);
    const auto *slg1_164 = buffer.data(slg1 + 164);
    const auto *slg1_170 = buffer.data(slg1 + 170);
    const auto *slg1_174 = buffer.data(slg1 + 174);
    const auto *slg1_177 = buffer.data(slg1 + 177);
    const auto *slg1_178 = buffer.data(slg1 + 178);
    const auto *slg1_179 = buffer.data(slg1 + 179);
    const auto *slg1_180 = buffer.data(slg1 + 180);
    const auto *slg1_183 = buffer.data(slg1 + 183);
    const auto *slg1_185 = buffer.data(slg1 + 185);
    const auto *slg1_186 = buffer.data(slg1 + 186);
    const auto *slg1_189 = buffer.data(slg1 + 189);
    const auto *slg1_190 = buffer.data(slg1 + 190);
    const auto *slg1_192 = buffer.data(slg1 + 192);
    const auto *slg1_194 = buffer.data(slg1 + 194);

    const auto *slh_183 = buffer.data(slh + 183);
    const auto *slh_184 = buffer.data(slh + 184);
    const auto *slh_185 = buffer.data(slh + 185);
    const auto *slh_186 = buffer.data(slh + 186);
    const auto *slh_187 = buffer.data(slh + 187);
    const auto *slh_188 = buffer.data(slh + 188);
    const auto *slh_189 = buffer.data(slh + 189);
    const auto *slh_191 = buffer.data(slh + 191);
    const auto *slh_192 = buffer.data(slh + 192);
    const auto *slh_194 = buffer.data(slh + 194);
    const auto *slh_195 = buffer.data(slh + 195);
    const auto *slh_198 = buffer.data(slh + 198);
    const auto *slh_199 = buffer.data(slh + 199);
    const auto *slh_201 = buffer.data(slh + 201);
    const auto *slh_203 = buffer.data(slh + 203);
    const auto *slh_204 = buffer.data(slh + 204);
    const auto *slh_205 = buffer.data(slh + 205);
    const auto *slh_206 = buffer.data(slh + 206);
    const auto *slh_207 = buffer.data(slh + 207);
    const auto *slh_208 = buffer.data(slh + 208);
    const auto *slh_209 = buffer.data(slh + 209);
    const auto *slh_210 = buffer.data(slh + 210);
    const auto *slh_212 = buffer.data(slh + 212);
    const auto *slh_213 = buffer.data(slh + 213);
    const auto *slh_215 = buffer.data(slh + 215);
    const auto *slh_216 = buffer.data(slh + 216);
    const auto *slh_219 = buffer.data(slh + 219);
    const auto *slh_220 = buffer.data(slh + 220);
    const auto *slh_222 = buffer.data(slh + 222);
    const auto *slh_224 = buffer.data(slh + 224);
    const auto *slh_225 = buffer.data(slh + 225);
    const auto *slh_226 = buffer.data(slh + 226);
    const auto *slh_227 = buffer.data(slh + 227);
    const auto *slh_228 = buffer.data(slh + 228);
    const auto *slh_229 = buffer.data(slh + 229);
    const auto *slh_230 = buffer.data(slh + 230);
    const auto *slh_231 = buffer.data(slh + 231);
    const auto *slh_233 = buffer.data(slh + 233);
    const auto *slh_234 = buffer.data(slh + 234);
    const auto *slh_236 = buffer.data(slh + 236);
    const auto *slh_237 = buffer.data(slh + 237);
    const auto *slh_240 = buffer.data(slh + 240);
    const auto *slh_245 = buffer.data(slh + 245);
    const auto *slh_246 = buffer.data(slh + 246);
    const auto *slh_247 = buffer.data(slh + 247);
    const auto *slh_248 = buffer.data(slh + 248);
    const auto *slh_249 = buffer.data(slh + 249);
    const auto *slh_250 = buffer.data(slh + 250);
    const auto *slh_251 = buffer.data(slh + 251);
    const auto *slh_252 = buffer.data(slh + 252);
    const auto *slh_254 = buffer.data(slh + 254);
    const auto *slh_255 = buffer.data(slh + 255);
    const auto *slh_257 = buffer.data(slh + 257);
    const auto *slh_258 = buffer.data(slh + 258);
    const auto *slh_261 = buffer.data(slh + 261);
    const auto *slh_262 = buffer.data(slh + 262);
    const auto *slh_264 = buffer.data(slh + 264);
    const auto *slh_266 = buffer.data(slh + 266);

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, pc_x, skh_184, skh_185, skh_186, \
                         skh_187, skh_188, slh_184, slh_185, slh_186, slh_187, \
                         slh_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_17 * skh_184[k]
                   + f_3 * pc_x[k] * slh_184[k];

        t_241[k] = f_17 * skh_185[k]
                   + f_3 * pc_x[k] * slh_185[k];

        t_242[k] = f_17 * skh_186[k]
                   + f_3 * pc_x[k] * slh_186[k];

        t_243[k] = f_17 * skh_187[k]
                   + f_3 * pc_x[k] * slh_187[k];

        t_244[k] = f_17 * skh_188[k]
                   + f_3 * pc_x[k] * slh_188[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, pc_y, pc_z, skh_99, skh_120, skh_122, slg0_130, \
                         slg0_132, slg1_130, slg1_132, slh_183, \
                         slh_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_11 * skh_120[k]
                   + f_1 * slg0_130[k]
                   - f_2 * slg1_130[k]
                   + f_3 * pc_y[k] * slh_183[k];

        t_246[k] = f_12 * skh_99[k]
                   + f_3 * pc_z[k] * slh_183[k];

        t_247[k] = f_11 * skh_122[k]
                   + f_4 * slg0_132[k]
                   - f_5 * slg1_132[k]
                   + f_3 * pc_y[k] * slh_185[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, pc_y, skh_123, skh_124, skh_125, slg0_133, \
                         slg0_134, slg1_133, slg1_134, slh_186, slh_187, \
                         slh_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_11 * skh_123[k]
                   + f_6 * slg0_133[k]
                   - f_7 * slg1_133[k]
                   + f_3 * pc_y[k] * slh_186[k];

        t_249[k] = f_11 * skh_124[k]
                   + f_8 * slg0_134[k]
                   - f_9 * slg1_134[k]
                   + f_3 * pc_y[k] * slh_187[k];

        t_250[k] = f_11 * skh_125[k]
                   + f_3 * pc_y[k] * slh_188[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, pb_y, pc_x, pc_y, pc_z, ski0_167, \
                         skh_105, skh_189, ski1_167, slg0_135, slg1_135, \
                         slh_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = pb_y[k] * ski0_167[k]
                   - f_10 * pc_y[k] * ski1_167[k];

        t_252[k] = f_17 * skh_189[k]
                   + f_1 * slg0_135[k]
                   - f_2 * slg1_135[k]
                   + f_3 * pc_x[k] * slh_189[k];

        t_253[k] = f_3 * pc_y[k] * slh_189[k];

        t_254[k] = f_13 * skh_105[k]
                   + f_3 * pc_z[k] * slh_189[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, pc_x, pc_y, skh_192, skh_194, slg0_138, \
                         slg0_140, slg1_138, slg1_140, slh_191, slh_192, \
                         slh_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_17 * skh_192[k]
                   + f_4 * slg0_138[k]
                   - f_5 * slg1_138[k]
                   + f_3 * pc_x[k] * slh_192[k];

        t_256[k] = f_3 * pc_y[k] * slh_191[k];

        t_257[k] = f_17 * skh_194[k]
                   + f_4 * slg0_140[k]
                   - f_5 * slg1_140[k]
                   + f_3 * pc_x[k] * slh_194[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, pc_x, pc_y, pc_z, skh_108, skh_195, slg0_141, \
                         slg1_141, slh_192, slh_194, slh_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = f_17 * skh_195[k]
                   + f_6 * slg0_141[k]
                   - f_7 * slg1_141[k]
                   + f_3 * pc_x[k] * slh_195[k];

        t_259[k] = f_13 * skh_108[k]
                   + f_3 * pc_z[k] * slh_192[k];

        t_260[k] = f_3 * pc_y[k] * slh_194[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, pc_x, pc_z, skh_111, skh_198, skh_199, slg0_144, \
                         slg0_145, slg1_144, slg1_145, slh_195, slh_198, \
                         slh_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_17 * skh_198[k]
                   + f_6 * slg0_144[k]
                   - f_7 * slg1_144[k]
                   + f_3 * pc_x[k] * slh_198[k];

        t_262[k] = f_17 * skh_199[k]
                   + f_8 * slg0_145[k]
                   - f_9 * slg1_145[k]
                   + f_3 * pc_x[k] * slh_199[k];

        t_263[k] = f_13 * skh_111[k]
                   + f_3 * pc_z[k] * slh_195[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, pc_x, pc_y, skh_201, skh_203, slg0_147, \
                         slg0_149, slg1_147, slg1_149, slh_198, slh_201, \
                         slh_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_17 * skh_201[k]
                   + f_8 * slg0_147[k]
                   - f_9 * slg1_147[k]
                   + f_3 * pc_x[k] * slh_201[k];

        t_265[k] = f_3 * pc_y[k] * slh_198[k];

        t_266[k] = f_17 * skh_203[k]
                   + f_8 * slg0_149[k]
                   - f_9 * slg1_149[k]
                   + f_3 * pc_x[k] * slh_203[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, t_270, t_271, pc_x, skh_204, skh_205, skh_206, \
                         skh_207, skh_208, slh_204, slh_205, slh_206, slh_207, \
                         slh_208 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_17 * skh_204[k]
                   + f_3 * pc_x[k] * slh_204[k];

        t_268[k] = f_17 * skh_205[k]
                   + f_3 * pc_x[k] * slh_205[k];

        t_269[k] = f_17 * skh_206[k]
                   + f_3 * pc_x[k] * slh_206[k];

        t_270[k] = f_17 * skh_207[k]
                   + f_3 * pc_x[k] * slh_207[k];

        t_271[k] = f_17 * skh_208[k]
                   + f_3 * pc_x[k] * slh_208[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, pc_x, pc_y, pc_z, skh_120, skh_209, \
                         slg0_145, slg0_147, slg1_145, slg1_147, slh_204, slh_206, \
                         slh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_17 * skh_209[k]
                   + f_3 * pc_x[k] * slh_209[k];

        t_273[k] = f_1 * slg0_145[k]
                   - f_2 * slg1_145[k]
                   + f_3 * pc_y[k] * slh_204[k];

        t_274[k] = f_13 * skh_120[k]
                   + f_3 * pc_z[k] * slh_204[k];

        t_275[k] = f_4 * slg0_147[k]
                   - f_5 * slg1_147[k]
                   + f_3 * pc_y[k] * slh_206[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pc_y, pc_z, skh_125, slg0_148, slg0_149, \
                         slg1_148, slg1_149, slh_207, slh_208, \
                         slh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_6 * slg0_148[k]
                   - f_7 * slg1_148[k]
                   + f_3 * pc_y[k] * slh_207[k];

        t_277[k] = f_8 * slg0_149[k]
                   - f_9 * slg1_149[k]
                   + f_3 * pc_y[k] * slh_208[k];

        t_278[k] = f_3 * pc_y[k] * slh_209[k];

        t_279[k] = f_13 * skh_125[k]
                   + f_1 * slg0_149[k]
                   - f_2 * slg1_149[k]
                   + f_3 * pc_z[k] * slh_209[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, pc_x, pc_y, pc_z, skh_126, skh_210, \
                         skh_213, slg0_150, slg0_153, slg1_150, slg1_153, slh_210, \
                         slh_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_14 * skh_210[k]
                   + f_1 * slg0_150[k]
                   - f_2 * slg1_150[k]
                   + f_3 * pc_x[k] * slh_210[k];

        t_281[k] = f_14 * skh_126[k]
                   + f_3 * pc_y[k] * slh_210[k];

        t_282[k] = f_3 * pc_z[k] * slh_210[k];

        t_283[k] = f_14 * skh_213[k]
                   + f_4 * slg0_153[k]
                   - f_5 * slg1_153[k]
                   + f_3 * pc_x[k] * slh_213[k];
    }

#pragma omp simd aligned(t_284, t_285, t_286, pc_x, pc_y, skh_128, skh_215, skh_216, slg0_155, \
                         slg0_156, slg1_155, slg1_156, slh_212, slh_215, \
                         slh_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_284[k] = f_14 * skh_128[k]
                   + f_3 * pc_y[k] * slh_212[k];

        t_285[k] = f_14 * skh_215[k]
                   + f_4 * slg0_155[k]
                   - f_5 * slg1_155[k]
                   + f_3 * pc_x[k] * slh_215[k];

        t_286[k] = f_14 * skh_216[k]
                   + f_6 * slg0_156[k]
                   - f_7 * slg1_156[k]
                   + f_3 * pc_x[k] * slh_216[k];
    }

#pragma omp simd aligned(t_287, t_288, t_289, pc_x, pc_y, pc_z, skh_131, skh_219, slg0_159, \
                         slg1_159, slh_213, slh_215, slh_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_287[k] = f_3 * pc_z[k] * slh_213[k];

        t_288[k] = f_14 * skh_131[k]
                   + f_3 * pc_y[k] * slh_215[k];

        t_289[k] = f_14 * skh_219[k]
                   + f_6 * slg0_159[k]
                   - f_7 * slg1_159[k]
                   + f_3 * pc_x[k] * slh_219[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, pc_x, pc_z, skh_220, skh_222, slg0_160, \
                         slg0_162, slg1_160, slg1_162, slh_216, slh_220, \
                         slh_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = f_14 * skh_220[k]
                   + f_8 * slg0_160[k]
                   - f_9 * slg1_160[k]
                   + f_3 * pc_x[k] * slh_220[k];

        t_291[k] = f_3 * pc_z[k] * slh_216[k];

        t_292[k] = f_14 * skh_222[k]
                   + f_8 * slg0_162[k]
                   - f_9 * slg1_162[k]
                   + f_3 * pc_x[k] * slh_222[k];
    }

#pragma omp simd aligned(t_293, t_294, t_295, t_296, pc_x, pc_y, skh_135, skh_224, skh_225, \
                         skh_226, slg0_164, slg1_164, slh_219, slh_224, slh_225, \
                         slh_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_293[k] = f_14 * skh_135[k]
                   + f_3 * pc_y[k] * slh_219[k];

        t_294[k] = f_14 * skh_224[k]
                   + f_8 * slg0_164[k]
                   - f_9 * slg1_164[k]
                   + f_3 * pc_x[k] * slh_224[k];

        t_295[k] = f_14 * skh_225[k]
                   + f_3 * pc_x[k] * slh_225[k];

        t_296[k] = f_14 * skh_226[k]
                   + f_3 * pc_x[k] * slh_226[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, t_300, pc_x, skh_227, skh_228, skh_229, skh_230, \
                         slh_227, slh_228, slh_229, slh_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = f_14 * skh_227[k]
                   + f_3 * pc_x[k] * slh_227[k];

        t_298[k] = f_14 * skh_228[k]
                   + f_3 * pc_x[k] * slh_228[k];

        t_299[k] = f_14 * skh_229[k]
                   + f_3 * pc_x[k] * slh_229[k];

        t_300[k] = f_14 * skh_230[k]
                   + f_3 * pc_x[k] * slh_230[k];
    }

#pragma omp simd aligned(t_301, t_302, t_303, pc_y, pc_z, skh_141, skh_143, slg0_160, \
                         slg0_162, slg1_160, slg1_162, slh_225, \
                         slh_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_301[k] = f_14 * skh_141[k]
                   + f_1 * slg0_160[k]
                   - f_2 * slg1_160[k]
                   + f_3 * pc_y[k] * slh_225[k];

        t_302[k] = f_3 * pc_z[k] * slh_225[k];

        t_303[k] = f_14 * skh_143[k]
                   + f_4 * slg0_162[k]
                   - f_5 * slg1_162[k]
                   + f_3 * pc_y[k] * slh_227[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, pc_y, pc_z, skh_144, skh_145, skh_146, \
                         slg0_163, slg0_164, slg1_163, slg1_164, slh_228, slh_229, \
                         slh_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_304[k] = f_14 * skh_144[k]
                   + f_6 * slg0_163[k]
                   - f_7 * slg1_163[k]
                   + f_3 * pc_y[k] * slh_228[k];

        t_305[k] = f_14 * skh_145[k]
                   + f_8 * slg0_164[k]
                   - f_9 * slg1_164[k]
                   + f_3 * pc_y[k] * slh_229[k];

        t_306[k] = f_14 * skh_146[k]
                   + f_3 * pc_y[k] * slh_230[k];

        t_307[k] = f_1 * slg0_164[k]
                   - f_2 * slg1_164[k]
                   + f_3 * pc_z[k] * slh_230[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, t_311, pb_z, pc_y, pc_z, ski0_168, ski0_171, \
                         skh_126, skh_147, ski1_168, ski1_171, \
                         slh_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = pb_z[k] * ski0_168[k]
                   - f_10 * pc_z[k] * ski1_168[k];

        t_309[k] = f_13 * skh_147[k]
                   + f_3 * pc_y[k] * slh_231[k];

        t_310[k] = f_11 * skh_126[k]
                   + f_3 * pc_z[k] * slh_231[k];

        t_311[k] = pb_z[k] * ski0_171[k]
                   - f_10 * pc_z[k] * ski1_171[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, pb_z, pc_x, pc_y, pc_z, ski0_174, skh_149, \
                         skh_236, ski1_174, slg0_170, slg1_170, slh_233, \
                         slh_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = f_13 * skh_149[k]
                   + f_3 * pc_y[k] * slh_233[k];

        t_313[k] = f_14 * skh_236[k]
                   + f_4 * slg0_170[k]
                   - f_5 * slg1_170[k]
                   + f_3 * pc_x[k] * slh_236[k];

        t_314[k] = pb_z[k] * ski0_174[k]
                   - f_10 * pc_z[k] * ski1_174[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, pc_x, pc_y, pc_z, skh_129, skh_152, skh_240, \
                         slg0_174, slg1_174, slh_234, slh_236, \
                         slh_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = f_11 * skh_129[k]
                   + f_3 * pc_z[k] * slh_234[k];

        t_316[k] = f_13 * skh_152[k]
                   + f_3 * pc_y[k] * slh_236[k];

        t_317[k] = f_14 * skh_240[k]
                   + f_6 * slg0_174[k]
                   - f_7 * slg1_174[k]
                   + f_3 * pc_x[k] * slh_240[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, t_321, pb_z, pc_y, pc_z, ski0_178, ski0_180, \
                         skh_132, skh_133, skh_156, ski1_178, ski1_180, slh_237, \
                         slh_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = pb_z[k] * ski0_178[k]
                   - f_10 * pc_z[k] * ski1_178[k];

        t_319[k] = f_11 * skh_132[k]
                   + f_3 * pc_z[k] * slh_237[k];

        t_320[k] = pb_z[k] * ski0_180[k]
                   + f_12 * skh_133[k]
                   - f_10 * pc_z[k] * ski1_180[k];

        t_321[k] = f_13 * skh_156[k]
                   + f_3 * pc_y[k] * slh_240[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, t_325, pc_x, skh_245, skh_246, skh_247, skh_248, \
                         slg0_179, slg1_179, slh_245, slh_246, slh_247, \
                         slh_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = f_14 * skh_245[k]
                   + f_8 * slg0_179[k]
                   - f_9 * slg1_179[k]
                   + f_3 * pc_x[k] * slh_245[k];

        t_323[k] = f_14 * skh_246[k]
                   + f_3 * pc_x[k] * slh_246[k];

        t_324[k] = f_14 * skh_247[k]
                   + f_3 * pc_x[k] * slh_247[k];

        t_325[k] = f_14 * skh_248[k]
                   + f_3 * pc_x[k] * slh_248[k];
    }

#pragma omp simd aligned(t_326, t_327, t_328, t_329, pb_z, pc_x, pc_z, ski0_189, skh_249, \
                         skh_250, skh_251, ski1_189, slh_249, slh_250, \
                         slh_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_326[k] = f_14 * skh_249[k]
                   + f_3 * pc_x[k] * slh_249[k];

        t_327[k] = f_14 * skh_250[k]
                   + f_3 * pc_x[k] * slh_250[k];

        t_328[k] = f_14 * skh_251[k]
                   + f_3 * pc_x[k] * slh_251[k];

        t_329[k] = pb_z[k] * ski0_189[k]
                   - f_10 * pc_z[k] * ski1_189[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, pc_y, pc_z, skh_141, skh_164, skh_165, slg0_177, \
                         slg0_178, slg1_177, slg1_178, slh_246, slh_248, \
                         slh_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = f_11 * skh_141[k]
                   + f_3 * pc_z[k] * slh_246[k];

        t_331[k] = f_13 * skh_164[k]
                   + f_4 * slg0_177[k]
                   - f_5 * slg1_177[k]
                   + f_3 * pc_y[k] * slh_248[k];

        t_332[k] = f_13 * skh_165[k]
                   + f_6 * slg0_178[k]
                   - f_7 * slg1_178[k]
                   + f_3 * pc_y[k] * slh_249[k];
    }

#pragma omp simd aligned(t_333, t_334, t_335, pc_y, pc_z, skh_146, skh_166, skh_167, slg0_179, \
                         slg1_179, slh_250, slh_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_333[k] = f_13 * skh_166[k]
                   + f_8 * slg0_179[k]
                   - f_9 * slg1_179[k]
                   + f_3 * pc_y[k] * slh_250[k];

        t_334[k] = f_13 * skh_167[k]
                   + f_3 * pc_y[k] * slh_251[k];

        t_335[k] = f_11 * skh_146[k]
                   + f_1 * slg0_179[k]
                   - f_2 * slg1_179[k]
                   + f_3 * pc_z[k] * slh_251[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, pc_x, pc_y, pc_z, skh_147, skh_168, skh_252, \
                         slg0_180, slg1_180, slh_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = f_14 * skh_252[k]
                   + f_1 * slg0_180[k]
                   - f_2 * slg1_180[k]
                   + f_3 * pc_x[k] * slh_252[k];

        t_337[k] = f_12 * skh_168[k]
                   + f_3 * pc_y[k] * slh_252[k];

        t_338[k] = f_12 * skh_147[k]
                   + f_3 * pc_z[k] * slh_252[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, pc_x, pc_y, skh_170, skh_255, skh_257, slg0_183, \
                         slg0_185, slg1_183, slg1_185, slh_254, slh_255, \
                         slh_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = f_14 * skh_255[k]
                   + f_4 * slg0_183[k]
                   - f_5 * slg1_183[k]
                   + f_3 * pc_x[k] * slh_255[k];

        t_340[k] = f_12 * skh_170[k]
                   + f_3 * pc_y[k] * slh_254[k];

        t_341[k] = f_14 * skh_257[k]
                   + f_4 * slg0_185[k]
                   - f_5 * slg1_185[k]
                   + f_3 * pc_x[k] * slh_257[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, pc_x, pc_y, pc_z, skh_150, skh_173, skh_258, \
                         slg0_186, slg1_186, slh_255, slh_257, \
                         slh_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = f_14 * skh_258[k]
                   + f_6 * slg0_186[k]
                   - f_7 * slg1_186[k]
                   + f_3 * pc_x[k] * slh_258[k];

        t_343[k] = f_12 * skh_150[k]
                   + f_3 * pc_z[k] * slh_255[k];

        t_344[k] = f_12 * skh_173[k]
                   + f_3 * pc_y[k] * slh_257[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, pc_x, pc_z, skh_153, skh_261, skh_262, slg0_189, \
                         slg0_190, slg1_189, slg1_190, slh_258, slh_261, \
                         slh_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_14 * skh_261[k]
                   + f_6 * slg0_189[k]
                   - f_7 * slg1_189[k]
                   + f_3 * pc_x[k] * slh_261[k];

        t_346[k] = f_14 * skh_262[k]
                   + f_8 * slg0_190[k]
                   - f_9 * slg1_190[k]
                   + f_3 * pc_x[k] * slh_262[k];

        t_347[k] = f_12 * skh_153[k]
                   + f_3 * pc_z[k] * slh_258[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, pc_x, pc_y, skh_177, skh_264, skh_266, slg0_192, \
                         slg0_194, slg1_192, slg1_194, slh_261, slh_264, \
                         slh_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = f_14 * skh_264[k]
                   + f_8 * slg0_192[k]
                   - f_9 * slg1_192[k]
                   + f_3 * pc_x[k] * slh_264[k];

        t_349[k] = f_12 * skh_177[k]
                   + f_3 * pc_y[k] * slh_261[k];

        t_350[k] = f_14 * skh_266[k]
                   + f_8 * slg0_194[k]
                   - f_9 * slg1_194[k]
                   + f_3 * pc_x[k] * slh_266[k];
    }
}

static auto
compute_prim_sli_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t ski0,
                                                          const size_t skh, const size_t ski1,
                                                          const size_t slg0, const size_t slg1,
                                                          const size_t slh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.5 / gamma;
    const auto f_5 = 1.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 0.5 / gamma;
    const auto f_9 = 0.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_17 = 2.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ski0_252 = buffer.data(ski0 + 252);
    const auto *ski0_255 = buffer.data(ski0 + 255);
    const auto *ski0_257 = buffer.data(ski0 + 257);
    const auto *ski0_258 = buffer.data(ski0 + 258);
    const auto *ski0_261 = buffer.data(ski0 + 261);
    const auto *ski0_262 = buffer.data(ski0 + 262);
    const auto *ski0_264 = buffer.data(ski0 + 264);
    const auto *ski0_266 = buffer.data(ski0 + 266);
    const auto *ski0_279 = buffer.data(ski0 + 279);
    const auto *ski0_280 = buffer.data(ski0 + 280);
    const auto *ski0_283 = buffer.data(ski0 + 283);
    const auto *ski0_286 = buffer.data(ski0 + 286);
    const auto *ski0_290 = buffer.data(ski0 + 290);
    const auto *ski0_292 = buffer.data(ski0 + 292);

    const auto *skh_162 = buffer.data(skh + 162);
    const auto *skh_167 = buffer.data(skh + 167);
    const auto *skh_168 = buffer.data(skh + 168);
    const auto *skh_171 = buffer.data(skh + 171);
    const auto *skh_174 = buffer.data(skh + 174);
    const auto *skh_183 = buffer.data(skh + 183);
    const auto *skh_185 = buffer.data(skh + 185);
    const auto *skh_186 = buffer.data(skh + 186);
    const auto *skh_187 = buffer.data(skh + 187);
    const auto *skh_188 = buffer.data(skh + 188);
    const auto *skh_189 = buffer.data(skh + 189);
    const auto *skh_190 = buffer.data(skh + 190);
    const auto *skh_191 = buffer.data(skh + 191);
    const auto *skh_192 = buffer.data(skh + 192);
    const auto *skh_194 = buffer.data(skh + 194);
    const auto *skh_195 = buffer.data(skh + 195);
    const auto *skh_197 = buffer.data(skh + 197);
    const auto *skh_198 = buffer.data(skh + 198);
    const auto *skh_204 = buffer.data(skh + 204);
    const auto *skh_206 = buffer.data(skh + 206);
    const auto *skh_207 = buffer.data(skh + 207);
    const auto *skh_208 = buffer.data(skh + 208);
    const auto *skh_209 = buffer.data(skh + 209);
    const auto *skh_210 = buffer.data(skh + 210);
    const auto *skh_212 = buffer.data(skh + 212);
    const auto *skh_213 = buffer.data(skh + 213);
    const auto *skh_215 = buffer.data(skh + 215);
    const auto *skh_216 = buffer.data(skh + 216);
    const auto *skh_217 = buffer.data(skh + 217);
    const auto *skh_219 = buffer.data(skh + 219);
    const auto *skh_225 = buffer.data(skh + 225);
    const auto *skh_227 = buffer.data(skh + 227);
    const auto *skh_228 = buffer.data(skh + 228);
    const auto *skh_229 = buffer.data(skh + 229);
    const auto *skh_230 = buffer.data(skh + 230);
    const auto *skh_231 = buffer.data(skh + 231);
    const auto *skh_233 = buffer.data(skh + 233);
    const auto *skh_236 = buffer.data(skh + 236);
    const auto *skh_240 = buffer.data(skh + 240);
    const auto *skh_267 = buffer.data(skh + 267);
    const auto *skh_268 = buffer.data(skh + 268);
    const auto *skh_269 = buffer.data(skh + 269);
    const auto *skh_270 = buffer.data(skh + 270);
    const auto *skh_271 = buffer.data(skh + 271);
    const auto *skh_272 = buffer.data(skh + 272);
    const auto *skh_288 = buffer.data(skh + 288);
    const auto *skh_289 = buffer.data(skh + 289);
    const auto *skh_290 = buffer.data(skh + 290);
    const auto *skh_291 = buffer.data(skh + 291);
    const auto *skh_292 = buffer.data(skh + 292);
    const auto *skh_293 = buffer.data(skh + 293);
    const auto *skh_294 = buffer.data(skh + 294);
    const auto *skh_297 = buffer.data(skh + 297);
    const auto *skh_299 = buffer.data(skh + 299);
    const auto *skh_300 = buffer.data(skh + 300);
    const auto *skh_303 = buffer.data(skh + 303);
    const auto *skh_304 = buffer.data(skh + 304);
    const auto *skh_306 = buffer.data(skh + 306);
    const auto *skh_308 = buffer.data(skh + 308);
    const auto *skh_309 = buffer.data(skh + 309);
    const auto *skh_310 = buffer.data(skh + 310);
    const auto *skh_311 = buffer.data(skh + 311);
    const auto *skh_312 = buffer.data(skh + 312);
    const auto *skh_313 = buffer.data(skh + 313);
    const auto *skh_314 = buffer.data(skh + 314);
    const auto *skh_315 = buffer.data(skh + 315);
    const auto *skh_318 = buffer.data(skh + 318);
    const auto *skh_320 = buffer.data(skh + 320);
    const auto *skh_321 = buffer.data(skh + 321);
    const auto *skh_324 = buffer.data(skh + 324);
    const auto *skh_325 = buffer.data(skh + 325);
    const auto *skh_327 = buffer.data(skh + 327);
    const auto *skh_329 = buffer.data(skh + 329);
    const auto *skh_330 = buffer.data(skh + 330);
    const auto *skh_331 = buffer.data(skh + 331);
    const auto *skh_332 = buffer.data(skh + 332);
    const auto *skh_333 = buffer.data(skh + 333);
    const auto *skh_334 = buffer.data(skh + 334);
    const auto *skh_335 = buffer.data(skh + 335);
    const auto *skh_341 = buffer.data(skh + 341);
    const auto *skh_345 = buffer.data(skh + 345);
    const auto *skh_350 = buffer.data(skh + 350);
    const auto *skh_351 = buffer.data(skh + 351);
    const auto *skh_352 = buffer.data(skh + 352);
    const auto *skh_353 = buffer.data(skh + 353);

    const auto *ski1_252 = buffer.data(ski1 + 252);
    const auto *ski1_255 = buffer.data(ski1 + 255);
    const auto *ski1_257 = buffer.data(ski1 + 257);
    const auto *ski1_258 = buffer.data(ski1 + 258);
    const auto *ski1_261 = buffer.data(ski1 + 261);
    const auto *ski1_262 = buffer.data(ski1 + 262);
    const auto *ski1_264 = buffer.data(ski1 + 264);
    const auto *ski1_266 = buffer.data(ski1 + 266);
    const auto *ski1_279 = buffer.data(ski1 + 279);
    const auto *ski1_280 = buffer.data(ski1 + 280);
    const auto *ski1_283 = buffer.data(ski1 + 283);
    const auto *ski1_286 = buffer.data(ski1 + 286);
    const auto *ski1_290 = buffer.data(ski1 + 290);
    const auto *ski1_292 = buffer.data(ski1 + 292);

    const auto *slg0_190 = buffer.data(slg0 + 190);
    const auto *slg0_192 = buffer.data(slg0 + 192);
    const auto *slg0_193 = buffer.data(slg0 + 193);
    const auto *slg0_194 = buffer.data(slg0 + 194);
    const auto *slg0_205 = buffer.data(slg0 + 205);
    const auto *slg0_207 = buffer.data(slg0 + 207);
    const auto *slg0_208 = buffer.data(slg0 + 208);
    const auto *slg0_209 = buffer.data(slg0 + 209);
    const auto *slg0_210 = buffer.data(slg0 + 210);
    const auto *slg0_213 = buffer.data(slg0 + 213);
    const auto *slg0_215 = buffer.data(slg0 + 215);
    const auto *slg0_216 = buffer.data(slg0 + 216);
    const auto *slg0_219 = buffer.data(slg0 + 219);
    const auto *slg0_220 = buffer.data(slg0 + 220);
    const auto *slg0_222 = buffer.data(slg0 + 222);
    const auto *slg0_223 = buffer.data(slg0 + 223);
    const auto *slg0_224 = buffer.data(slg0 + 224);
    const auto *slg0_225 = buffer.data(slg0 + 225);
    const auto *slg0_228 = buffer.data(slg0 + 228);
    const auto *slg0_230 = buffer.data(slg0 + 230);
    const auto *slg0_231 = buffer.data(slg0 + 231);
    const auto *slg0_234 = buffer.data(slg0 + 234);
    const auto *slg0_235 = buffer.data(slg0 + 235);
    const auto *slg0_237 = buffer.data(slg0 + 237);
    const auto *slg0_238 = buffer.data(slg0 + 238);
    const auto *slg0_239 = buffer.data(slg0 + 239);
    const auto *slg0_245 = buffer.data(slg0 + 245);
    const auto *slg0_249 = buffer.data(slg0 + 249);
    const auto *slg0_254 = buffer.data(slg0 + 254);

    const auto *slg1_190 = buffer.data(slg1 + 190);
    const auto *slg1_192 = buffer.data(slg1 + 192);
    const auto *slg1_193 = buffer.data(slg1 + 193);
    const auto *slg1_194 = buffer.data(slg1 + 194);
    const auto *slg1_205 = buffer.data(slg1 + 205);
    const auto *slg1_207 = buffer.data(slg1 + 207);
    const auto *slg1_208 = buffer.data(slg1 + 208);
    const auto *slg1_209 = buffer.data(slg1 + 209);
    const auto *slg1_210 = buffer.data(slg1 + 210);
    const auto *slg1_213 = buffer.data(slg1 + 213);
    const auto *slg1_215 = buffer.data(slg1 + 215);
    const auto *slg1_216 = buffer.data(slg1 + 216);
    const auto *slg1_219 = buffer.data(slg1 + 219);
    const auto *slg1_220 = buffer.data(slg1 + 220);
    const auto *slg1_222 = buffer.data(slg1 + 222);
    const auto *slg1_223 = buffer.data(slg1 + 223);
    const auto *slg1_224 = buffer.data(slg1 + 224);
    const auto *slg1_225 = buffer.data(slg1 + 225);
    const auto *slg1_228 = buffer.data(slg1 + 228);
    const auto *slg1_230 = buffer.data(slg1 + 230);
    const auto *slg1_231 = buffer.data(slg1 + 231);
    const auto *slg1_234 = buffer.data(slg1 + 234);
    const auto *slg1_235 = buffer.data(slg1 + 235);
    const auto *slg1_237 = buffer.data(slg1 + 237);
    const auto *slg1_238 = buffer.data(slg1 + 238);
    const auto *slg1_239 = buffer.data(slg1 + 239);
    const auto *slg1_245 = buffer.data(slg1 + 245);
    const auto *slg1_249 = buffer.data(slg1 + 249);
    const auto *slg1_254 = buffer.data(slg1 + 254);

    const auto *slh_267 = buffer.data(slh + 267);
    const auto *slh_268 = buffer.data(slh + 268);
    const auto *slh_269 = buffer.data(slh + 269);
    const auto *slh_270 = buffer.data(slh + 270);
    const auto *slh_271 = buffer.data(slh + 271);
    const auto *slh_272 = buffer.data(slh + 272);
    const auto *slh_273 = buffer.data(slh + 273);
    const auto *slh_275 = buffer.data(slh + 275);
    const auto *slh_276 = buffer.data(slh + 276);
    const auto *slh_278 = buffer.data(slh + 278);
    const auto *slh_279 = buffer.data(slh + 279);
    const auto *slh_282 = buffer.data(slh + 282);
    const auto *slh_288 = buffer.data(slh + 288);
    const auto *slh_289 = buffer.data(slh + 289);
    const auto *slh_290 = buffer.data(slh + 290);
    const auto *slh_291 = buffer.data(slh + 291);
    const auto *slh_292 = buffer.data(slh + 292);
    const auto *slh_293 = buffer.data(slh + 293);
    const auto *slh_294 = buffer.data(slh + 294);
    const auto *slh_296 = buffer.data(slh + 296);
    const auto *slh_297 = buffer.data(slh + 297);
    const auto *slh_299 = buffer.data(slh + 299);
    const auto *slh_300 = buffer.data(slh + 300);
    const auto *slh_303 = buffer.data(slh + 303);
    const auto *slh_304 = buffer.data(slh + 304);
    const auto *slh_306 = buffer.data(slh + 306);
    const auto *slh_308 = buffer.data(slh + 308);
    const auto *slh_309 = buffer.data(slh + 309);
    const auto *slh_310 = buffer.data(slh + 310);
    const auto *slh_311 = buffer.data(slh + 311);
    const auto *slh_312 = buffer.data(slh + 312);
    const auto *slh_313 = buffer.data(slh + 313);
    const auto *slh_314 = buffer.data(slh + 314);
    const auto *slh_315 = buffer.data(slh + 315);
    const auto *slh_317 = buffer.data(slh + 317);
    const auto *slh_318 = buffer.data(slh + 318);
    const auto *slh_320 = buffer.data(slh + 320);
    const auto *slh_321 = buffer.data(slh + 321);
    const auto *slh_324 = buffer.data(slh + 324);
    const auto *slh_325 = buffer.data(slh + 325);
    const auto *slh_327 = buffer.data(slh + 327);
    const auto *slh_329 = buffer.data(slh + 329);
    const auto *slh_330 = buffer.data(slh + 330);
    const auto *slh_331 = buffer.data(slh + 331);
    const auto *slh_332 = buffer.data(slh + 332);
    const auto *slh_333 = buffer.data(slh + 333);
    const auto *slh_334 = buffer.data(slh + 334);
    const auto *slh_335 = buffer.data(slh + 335);
    const auto *slh_336 = buffer.data(slh + 336);
    const auto *slh_338 = buffer.data(slh + 338);
    const auto *slh_339 = buffer.data(slh + 339);
    const auto *slh_341 = buffer.data(slh + 341);
    const auto *slh_342 = buffer.data(slh + 342);
    const auto *slh_345 = buffer.data(slh + 345);
    const auto *slh_350 = buffer.data(slh + 350);
    const auto *slh_351 = buffer.data(slh + 351);
    const auto *slh_352 = buffer.data(slh + 352);
    const auto *slh_353 = buffer.data(slh + 353);

#pragma omp simd aligned(t_351, t_352, t_353, t_354, t_355, pc_x, skh_267, skh_268, skh_269, \
                         skh_270, skh_271, slh_267, slh_268, slh_269, slh_270, \
                         slh_271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_351[k] = f_14 * skh_267[k]
                   + f_3 * pc_x[k] * slh_267[k];

        t_352[k] = f_14 * skh_268[k]
                   + f_3 * pc_x[k] * slh_268[k];

        t_353[k] = f_14 * skh_269[k]
                   + f_3 * pc_x[k] * slh_269[k];

        t_354[k] = f_14 * skh_270[k]
                   + f_3 * pc_x[k] * slh_270[k];

        t_355[k] = f_14 * skh_271[k]
                   + f_3 * pc_x[k] * slh_271[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, pc_x, pc_y, pc_z, skh_162, skh_183, skh_272, \
                         slg0_190, slg1_190, slh_267, slh_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_14 * skh_272[k]
                   + f_3 * pc_x[k] * slh_272[k];

        t_357[k] = f_12 * skh_183[k]
                   + f_1 * slg0_190[k]
                   - f_2 * slg1_190[k]
                   + f_3 * pc_y[k] * slh_267[k];

        t_358[k] = f_12 * skh_162[k]
                   + f_3 * pc_z[k] * slh_267[k];
    }

#pragma omp simd aligned(t_359, t_360, t_361, pc_y, skh_185, skh_186, skh_187, slg0_192, \
                         slg0_193, slg0_194, slg1_192, slg1_193, slg1_194, slh_269, slh_270, \
                         slh_271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = f_12 * skh_185[k]
                   + f_4 * slg0_192[k]
                   - f_5 * slg1_192[k]
                   + f_3 * pc_y[k] * slh_269[k];

        t_360[k] = f_12 * skh_186[k]
                   + f_6 * slg0_193[k]
                   - f_7 * slg1_193[k]
                   + f_3 * pc_y[k] * slh_270[k];

        t_361[k] = f_12 * skh_187[k]
                   + f_8 * slg0_194[k]
                   - f_9 * slg1_194[k]
                   + f_3 * pc_y[k] * slh_271[k];
    }

#pragma omp simd aligned(t_362, t_363, t_364, t_365, pb_y, pc_y, pc_z, ski0_252, skh_167, \
                         skh_188, skh_189, ski1_252, slg0_194, slg1_194, slh_272, \
                         slh_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = f_12 * skh_188[k]
                   + f_3 * pc_y[k] * slh_272[k];

        t_363[k] = f_12 * skh_167[k]
                   + f_1 * slg0_194[k]
                   - f_2 * slg1_194[k]
                   + f_3 * pc_z[k] * slh_272[k];

        t_364[k] = pb_y[k] * ski0_252[k]
                   - f_10 * pc_y[k] * ski1_252[k];

        t_365[k] = f_11 * skh_189[k]
                   + f_3 * pc_y[k] * slh_273[k];
    }

#pragma omp simd aligned(t_366, t_367, t_368, t_369, pb_y, pc_y, pc_z, ski0_255, ski0_257, \
                         skh_168, skh_190, skh_191, ski1_255, ski1_257, slh_273, \
                         slh_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = f_13 * skh_168[k]
                   + f_3 * pc_z[k] * slh_273[k];

        t_367[k] = pb_y[k] * ski0_255[k]
                   + f_12 * skh_190[k]
                   - f_10 * pc_y[k] * ski1_255[k];

        t_368[k] = f_11 * skh_191[k]
                   + f_3 * pc_y[k] * slh_275[k];

        t_369[k] = pb_y[k] * ski0_257[k]
                   - f_10 * pc_y[k] * ski1_257[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, pb_y, pc_y, pc_z, ski0_258, ski0_261, \
                         skh_171, skh_192, skh_194, ski1_258, ski1_261, slh_276, \
                         slh_278 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = pb_y[k] * ski0_258[k]
                   + f_13 * skh_192[k]
                   - f_10 * pc_y[k] * ski1_258[k];

        t_371[k] = f_13 * skh_171[k]
                   + f_3 * pc_z[k] * slh_276[k];

        t_372[k] = f_11 * skh_194[k]
                   + f_3 * pc_y[k] * slh_278[k];

        t_373[k] = pb_y[k] * ski0_261[k]
                   - f_10 * pc_y[k] * ski1_261[k];
    }

#pragma omp simd aligned(t_374, t_375, t_376, pb_y, pc_y, pc_z, ski0_262, ski0_264, skh_174, \
                         skh_195, skh_197, ski1_262, ski1_264, \
                         slh_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_374[k] = pb_y[k] * ski0_262[k]
                   + f_14 * skh_195[k]
                   - f_10 * pc_y[k] * ski1_262[k];

        t_375[k] = f_13 * skh_174[k]
                   + f_3 * pc_z[k] * slh_279[k];

        t_376[k] = pb_y[k] * ski0_264[k]
                   + f_12 * skh_197[k]
                   - f_10 * pc_y[k] * ski1_264[k];
    }

#pragma omp simd aligned(t_377, t_378, t_379, t_380, pb_y, pc_x, pc_y, ski0_266, skh_198, \
                         skh_288, skh_289, ski1_266, slh_282, slh_288, \
                         slh_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_377[k] = f_11 * skh_198[k]
                   + f_3 * pc_y[k] * slh_282[k];

        t_378[k] = pb_y[k] * ski0_266[k]
                   - f_10 * pc_y[k] * ski1_266[k];

        t_379[k] = f_14 * skh_288[k]
                   + f_3 * pc_x[k] * slh_288[k];

        t_380[k] = f_14 * skh_289[k]
                   + f_3 * pc_x[k] * slh_289[k];
    }

#pragma omp simd aligned(t_381, t_382, t_383, t_384, pc_x, skh_290, skh_291, skh_292, skh_293, \
                         slh_290, slh_291, slh_292, slh_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_381[k] = f_14 * skh_290[k]
                   + f_3 * pc_x[k] * slh_290[k];

        t_382[k] = f_14 * skh_291[k]
                   + f_3 * pc_x[k] * slh_291[k];

        t_383[k] = f_14 * skh_292[k]
                   + f_3 * pc_x[k] * slh_292[k];

        t_384[k] = f_14 * skh_293[k]
                   + f_3 * pc_x[k] * slh_293[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, pc_y, pc_z, skh_183, skh_204, skh_206, slg0_205, \
                         slg0_207, slg1_205, slg1_207, slh_288, \
                         slh_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = f_11 * skh_204[k]
                   + f_1 * slg0_205[k]
                   - f_2 * slg1_205[k]
                   + f_3 * pc_y[k] * slh_288[k];

        t_386[k] = f_13 * skh_183[k]
                   + f_3 * pc_z[k] * slh_288[k];

        t_387[k] = f_11 * skh_206[k]
                   + f_4 * slg0_207[k]
                   - f_5 * slg1_207[k]
                   + f_3 * pc_y[k] * slh_290[k];
    }

#pragma omp simd aligned(t_388, t_389, t_390, pc_y, skh_207, skh_208, skh_209, slg0_208, \
                         slg0_209, slg1_208, slg1_209, slh_291, slh_292, \
                         slh_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_388[k] = f_11 * skh_207[k]
                   + f_6 * slg0_208[k]
                   - f_7 * slg1_208[k]
                   + f_3 * pc_y[k] * slh_291[k];

        t_389[k] = f_11 * skh_208[k]
                   + f_8 * slg0_209[k]
                   - f_9 * slg1_209[k]
                   + f_3 * pc_y[k] * slh_292[k];

        t_390[k] = f_11 * skh_209[k]
                   + f_3 * pc_y[k] * slh_293[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, t_394, pb_y, pc_x, pc_y, pc_z, ski0_279, \
                         skh_189, skh_294, ski1_279, slg0_210, slg1_210, \
                         slh_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = pb_y[k] * ski0_279[k]
                   - f_10 * pc_y[k] * ski1_279[k];

        t_392[k] = f_14 * skh_294[k]
                   + f_1 * slg0_210[k]
                   - f_2 * slg1_210[k]
                   + f_3 * pc_x[k] * slh_294[k];

        t_393[k] = f_3 * pc_y[k] * slh_294[k];

        t_394[k] = f_14 * skh_189[k]
                   + f_3 * pc_z[k] * slh_294[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, pc_x, pc_y, skh_297, skh_299, slg0_213, \
                         slg0_215, slg1_213, slg1_215, slh_296, slh_297, \
                         slh_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = f_14 * skh_297[k]
                   + f_4 * slg0_213[k]
                   - f_5 * slg1_213[k]
                   + f_3 * pc_x[k] * slh_297[k];

        t_396[k] = f_3 * pc_y[k] * slh_296[k];

        t_397[k] = f_14 * skh_299[k]
                   + f_4 * slg0_215[k]
                   - f_5 * slg1_215[k]
                   + f_3 * pc_x[k] * slh_299[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, pc_x, pc_y, pc_z, skh_192, skh_300, slg0_216, \
                         slg1_216, slh_297, slh_299, slh_300 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = f_14 * skh_300[k]
                   + f_6 * slg0_216[k]
                   - f_7 * slg1_216[k]
                   + f_3 * pc_x[k] * slh_300[k];

        t_399[k] = f_14 * skh_192[k]
                   + f_3 * pc_z[k] * slh_297[k];

        t_400[k] = f_3 * pc_y[k] * slh_299[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, pc_x, pc_z, skh_195, skh_303, skh_304, slg0_219, \
                         slg0_220, slg1_219, slg1_220, slh_300, slh_303, \
                         slh_304 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_14 * skh_303[k]
                   + f_6 * slg0_219[k]
                   - f_7 * slg1_219[k]
                   + f_3 * pc_x[k] * slh_303[k];

        t_402[k] = f_14 * skh_304[k]
                   + f_8 * slg0_220[k]
                   - f_9 * slg1_220[k]
                   + f_3 * pc_x[k] * slh_304[k];

        t_403[k] = f_14 * skh_195[k]
                   + f_3 * pc_z[k] * slh_300[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, pc_x, pc_y, skh_306, skh_308, slg0_222, \
                         slg0_224, slg1_222, slg1_224, slh_303, slh_306, \
                         slh_308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = f_14 * skh_306[k]
                   + f_8 * slg0_222[k]
                   - f_9 * slg1_222[k]
                   + f_3 * pc_x[k] * slh_306[k];

        t_405[k] = f_3 * pc_y[k] * slh_303[k];

        t_406[k] = f_14 * skh_308[k]
                   + f_8 * slg0_224[k]
                   - f_9 * slg1_224[k]
                   + f_3 * pc_x[k] * slh_308[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, t_410, t_411, pc_x, skh_309, skh_310, skh_311, \
                         skh_312, skh_313, slh_309, slh_310, slh_311, slh_312, \
                         slh_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = f_14 * skh_309[k]
                   + f_3 * pc_x[k] * slh_309[k];

        t_408[k] = f_14 * skh_310[k]
                   + f_3 * pc_x[k] * slh_310[k];

        t_409[k] = f_14 * skh_311[k]
                   + f_3 * pc_x[k] * slh_311[k];

        t_410[k] = f_14 * skh_312[k]
                   + f_3 * pc_x[k] * slh_312[k];

        t_411[k] = f_14 * skh_313[k]
                   + f_3 * pc_x[k] * slh_313[k];
    }

#pragma omp simd aligned(t_412, t_413, t_414, t_415, pc_x, pc_y, pc_z, skh_204, skh_314, \
                         slg0_220, slg0_222, slg1_220, slg1_222, slh_309, slh_311, \
                         slh_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_412[k] = f_14 * skh_314[k]
                   + f_3 * pc_x[k] * slh_314[k];

        t_413[k] = f_1 * slg0_220[k]
                   - f_2 * slg1_220[k]
                   + f_3 * pc_y[k] * slh_309[k];

        t_414[k] = f_14 * skh_204[k]
                   + f_3 * pc_z[k] * slh_309[k];

        t_415[k] = f_4 * slg0_222[k]
                   - f_5 * slg1_222[k]
                   + f_3 * pc_y[k] * slh_311[k];
    }

#pragma omp simd aligned(t_416, t_417, t_418, t_419, pc_y, pc_z, skh_209, slg0_223, slg0_224, \
                         slg1_223, slg1_224, slh_312, slh_313, \
                         slh_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_416[k] = f_6 * slg0_223[k]
                   - f_7 * slg1_223[k]
                   + f_3 * pc_y[k] * slh_312[k];

        t_417[k] = f_8 * slg0_224[k]
                   - f_9 * slg1_224[k]
                   + f_3 * pc_y[k] * slh_313[k];

        t_418[k] = f_3 * pc_y[k] * slh_314[k];

        t_419[k] = f_14 * skh_209[k]
                   + f_1 * slg0_224[k]
                   - f_2 * slg1_224[k]
                   + f_3 * pc_z[k] * slh_314[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, pc_x, pc_y, pc_z, skh_210, skh_315, \
                         skh_318, slg0_225, slg0_228, slg1_225, slg1_228, slh_315, \
                         slh_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_13 * skh_315[k]
                   + f_1 * slg0_225[k]
                   - f_2 * slg1_225[k]
                   + f_3 * pc_x[k] * slh_315[k];

        t_421[k] = f_17 * skh_210[k]
                   + f_3 * pc_y[k] * slh_315[k];

        t_422[k] = f_3 * pc_z[k] * slh_315[k];

        t_423[k] = f_13 * skh_318[k]
                   + f_4 * slg0_228[k]
                   - f_5 * slg1_228[k]
                   + f_3 * pc_x[k] * slh_318[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, pc_x, pc_y, skh_212, skh_320, skh_321, slg0_230, \
                         slg0_231, slg1_230, slg1_231, slh_317, slh_320, \
                         slh_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = f_17 * skh_212[k]
                   + f_3 * pc_y[k] * slh_317[k];

        t_425[k] = f_13 * skh_320[k]
                   + f_4 * slg0_230[k]
                   - f_5 * slg1_230[k]
                   + f_3 * pc_x[k] * slh_320[k];

        t_426[k] = f_13 * skh_321[k]
                   + f_6 * slg0_231[k]
                   - f_7 * slg1_231[k]
                   + f_3 * pc_x[k] * slh_321[k];
    }

#pragma omp simd aligned(t_427, t_428, t_429, pc_x, pc_y, pc_z, skh_215, skh_324, slg0_234, \
                         slg1_234, slh_318, slh_320, slh_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_427[k] = f_3 * pc_z[k] * slh_318[k];

        t_428[k] = f_17 * skh_215[k]
                   + f_3 * pc_y[k] * slh_320[k];

        t_429[k] = f_13 * skh_324[k]
                   + f_6 * slg0_234[k]
                   - f_7 * slg1_234[k]
                   + f_3 * pc_x[k] * slh_324[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, pc_x, pc_z, skh_325, skh_327, slg0_235, \
                         slg0_237, slg1_235, slg1_237, slh_321, slh_325, \
                         slh_327 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = f_13 * skh_325[k]
                   + f_8 * slg0_235[k]
                   - f_9 * slg1_235[k]
                   + f_3 * pc_x[k] * slh_325[k];

        t_431[k] = f_3 * pc_z[k] * slh_321[k];

        t_432[k] = f_13 * skh_327[k]
                   + f_8 * slg0_237[k]
                   - f_9 * slg1_237[k]
                   + f_3 * pc_x[k] * slh_327[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, t_436, pc_x, pc_y, skh_219, skh_329, skh_330, \
                         skh_331, slg0_239, slg1_239, slh_324, slh_329, slh_330, \
                         slh_331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = f_17 * skh_219[k]
                   + f_3 * pc_y[k] * slh_324[k];

        t_434[k] = f_13 * skh_329[k]
                   + f_8 * slg0_239[k]
                   - f_9 * slg1_239[k]
                   + f_3 * pc_x[k] * slh_329[k];

        t_435[k] = f_13 * skh_330[k]
                   + f_3 * pc_x[k] * slh_330[k];

        t_436[k] = f_13 * skh_331[k]
                   + f_3 * pc_x[k] * slh_331[k];
    }

#pragma omp simd aligned(t_437, t_438, t_439, t_440, pc_x, skh_332, skh_333, skh_334, skh_335, \
                         slh_332, slh_333, slh_334, slh_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_437[k] = f_13 * skh_332[k]
                   + f_3 * pc_x[k] * slh_332[k];

        t_438[k] = f_13 * skh_333[k]
                   + f_3 * pc_x[k] * slh_333[k];

        t_439[k] = f_13 * skh_334[k]
                   + f_3 * pc_x[k] * slh_334[k];

        t_440[k] = f_13 * skh_335[k]
                   + f_3 * pc_x[k] * slh_335[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, pc_y, pc_z, skh_225, skh_227, slg0_235, \
                         slg0_237, slg1_235, slg1_237, slh_330, \
                         slh_332 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = f_17 * skh_225[k]
                   + f_1 * slg0_235[k]
                   - f_2 * slg1_235[k]
                   + f_3 * pc_y[k] * slh_330[k];

        t_442[k] = f_3 * pc_z[k] * slh_330[k];

        t_443[k] = f_17 * skh_227[k]
                   + f_4 * slg0_237[k]
                   - f_5 * slg1_237[k]
                   + f_3 * pc_y[k] * slh_332[k];
    }

#pragma omp simd aligned(t_444, t_445, t_446, t_447, pc_y, pc_z, skh_228, skh_229, skh_230, \
                         slg0_238, slg0_239, slg1_238, slg1_239, slh_333, slh_334, \
                         slh_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_444[k] = f_17 * skh_228[k]
                   + f_6 * slg0_238[k]
                   - f_7 * slg1_238[k]
                   + f_3 * pc_y[k] * slh_333[k];

        t_445[k] = f_17 * skh_229[k]
                   + f_8 * slg0_239[k]
                   - f_9 * slg1_239[k]
                   + f_3 * pc_y[k] * slh_334[k];

        t_446[k] = f_17 * skh_230[k]
                   + f_3 * pc_y[k] * slh_335[k];

        t_447[k] = f_1 * slg0_239[k]
                   - f_2 * slg1_239[k]
                   + f_3 * pc_z[k] * slh_335[k];
    }

#pragma omp simd aligned(t_448, t_449, t_450, t_451, pb_z, pc_y, pc_z, ski0_280, ski0_283, \
                         skh_210, skh_231, ski1_280, ski1_283, \
                         slh_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_448[k] = pb_z[k] * ski0_280[k]
                   - f_10 * pc_z[k] * ski1_280[k];

        t_449[k] = f_14 * skh_231[k]
                   + f_3 * pc_y[k] * slh_336[k];

        t_450[k] = f_11 * skh_210[k]
                   + f_3 * pc_z[k] * slh_336[k];

        t_451[k] = pb_z[k] * ski0_283[k]
                   - f_10 * pc_z[k] * ski1_283[k];
    }

#pragma omp simd aligned(t_452, t_453, t_454, pb_z, pc_x, pc_y, pc_z, ski0_286, skh_233, \
                         skh_341, ski1_286, slg0_245, slg1_245, slh_338, \
                         slh_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_452[k] = f_14 * skh_233[k]
                   + f_3 * pc_y[k] * slh_338[k];

        t_453[k] = f_13 * skh_341[k]
                   + f_4 * slg0_245[k]
                   - f_5 * slg1_245[k]
                   + f_3 * pc_x[k] * slh_341[k];

        t_454[k] = pb_z[k] * ski0_286[k]
                   - f_10 * pc_z[k] * ski1_286[k];
    }

#pragma omp simd aligned(t_455, t_456, t_457, pc_x, pc_y, pc_z, skh_213, skh_236, skh_345, \
                         slg0_249, slg1_249, slh_339, slh_341, \
                         slh_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_455[k] = f_11 * skh_213[k]
                   + f_3 * pc_z[k] * slh_339[k];

        t_456[k] = f_14 * skh_236[k]
                   + f_3 * pc_y[k] * slh_341[k];

        t_457[k] = f_13 * skh_345[k]
                   + f_6 * slg0_249[k]
                   - f_7 * slg1_249[k]
                   + f_3 * pc_x[k] * slh_345[k];
    }

#pragma omp simd aligned(t_458, t_459, t_460, t_461, pb_z, pc_y, pc_z, ski0_290, ski0_292, \
                         skh_216, skh_217, skh_240, ski1_290, ski1_292, slh_342, \
                         slh_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_458[k] = pb_z[k] * ski0_290[k]
                   - f_10 * pc_z[k] * ski1_290[k];

        t_459[k] = f_11 * skh_216[k]
                   + f_3 * pc_z[k] * slh_342[k];

        t_460[k] = pb_z[k] * ski0_292[k]
                   + f_12 * skh_217[k]
                   - f_10 * pc_z[k] * ski1_292[k];

        t_461[k] = f_14 * skh_240[k]
                   + f_3 * pc_y[k] * slh_345[k];
    }

#pragma omp simd aligned(t_462, t_463, t_464, t_465, pc_x, skh_350, skh_351, skh_352, skh_353, \
                         slg0_254, slg1_254, slh_350, slh_351, slh_352, \
                         slh_353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_462[k] = f_13 * skh_350[k]
                   + f_8 * slg0_254[k]
                   - f_9 * slg1_254[k]
                   + f_3 * pc_x[k] * slh_350[k];

        t_463[k] = f_13 * skh_351[k]
                   + f_3 * pc_x[k] * slh_351[k];

        t_464[k] = f_13 * skh_352[k]
                   + f_3 * pc_x[k] * slh_352[k];

        t_465[k] = f_13 * skh_353[k]
                   + f_3 * pc_x[k] * slh_353[k];
    }
}

static auto
compute_prim_sli_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t ski0,
                                                          const size_t skh, const size_t ski1,
                                                          const size_t slg0, const size_t slg1,
                                                          const size_t slh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.5 / gamma;
    const auto f_5 = 1.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 0.5 / gamma;
    const auto f_9 = 0.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_17 = 2.5 / q;

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
    auto *t_572 = buffer.data(target + 572);
    auto *t_573 = buffer.data(target + 573);
    auto *t_574 = buffer.data(target + 574);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ski0_301 = buffer.data(ski0 + 301);
    const auto *ski0_392 = buffer.data(ski0 + 392);
    const auto *ski0_395 = buffer.data(ski0 + 395);
    const auto *ski0_397 = buffer.data(ski0 + 397);
    const auto *ski0_398 = buffer.data(ski0 + 398);
    const auto *ski0_401 = buffer.data(ski0 + 401);
    const auto *ski0_402 = buffer.data(ski0 + 402);
    const auto *ski0_404 = buffer.data(ski0 + 404);
    const auto *ski0_406 = buffer.data(ski0 + 406);
    const auto *ski0_419 = buffer.data(ski0 + 419);

    const auto *skh_225 = buffer.data(skh + 225);
    const auto *skh_230 = buffer.data(skh + 230);
    const auto *skh_231 = buffer.data(skh + 231);
    const auto *skh_234 = buffer.data(skh + 234);
    const auto *skh_237 = buffer.data(skh + 237);
    const auto *skh_246 = buffer.data(skh + 246);
    const auto *skh_248 = buffer.data(skh + 248);
    const auto *skh_249 = buffer.data(skh + 249);
    const auto *skh_250 = buffer.data(skh + 250);
    const auto *skh_251 = buffer.data(skh + 251);
    const auto *skh_252 = buffer.data(skh + 252);
    const auto *skh_254 = buffer.data(skh + 254);
    const auto *skh_255 = buffer.data(skh + 255);
    const auto *skh_257 = buffer.data(skh + 257);
    const auto *skh_258 = buffer.data(skh + 258);
    const auto *skh_261 = buffer.data(skh + 261);
    const auto *skh_267 = buffer.data(skh + 267);
    const auto *skh_269 = buffer.data(skh + 269);
    const auto *skh_270 = buffer.data(skh + 270);
    const auto *skh_271 = buffer.data(skh + 271);
    const auto *skh_272 = buffer.data(skh + 272);
    const auto *skh_273 = buffer.data(skh + 273);
    const auto *skh_275 = buffer.data(skh + 275);
    const auto *skh_276 = buffer.data(skh + 276);
    const auto *skh_278 = buffer.data(skh + 278);
    const auto *skh_279 = buffer.data(skh + 279);
    const auto *skh_282 = buffer.data(skh + 282);
    const auto *skh_288 = buffer.data(skh + 288);
    const auto *skh_290 = buffer.data(skh + 290);
    const auto *skh_291 = buffer.data(skh + 291);
    const auto *skh_292 = buffer.data(skh + 292);
    const auto *skh_293 = buffer.data(skh + 293);
    const auto *skh_294 = buffer.data(skh + 294);
    const auto *skh_295 = buffer.data(skh + 295);
    const auto *skh_296 = buffer.data(skh + 296);
    const auto *skh_297 = buffer.data(skh + 297);
    const auto *skh_299 = buffer.data(skh + 299);
    const auto *skh_300 = buffer.data(skh + 300);
    const auto *skh_302 = buffer.data(skh + 302);
    const auto *skh_303 = buffer.data(skh + 303);
    const auto *skh_309 = buffer.data(skh + 309);
    const auto *skh_311 = buffer.data(skh + 311);
    const auto *skh_312 = buffer.data(skh + 312);
    const auto *skh_313 = buffer.data(skh + 313);
    const auto *skh_314 = buffer.data(skh + 314);
    const auto *skh_354 = buffer.data(skh + 354);
    const auto *skh_355 = buffer.data(skh + 355);
    const auto *skh_356 = buffer.data(skh + 356);
    const auto *skh_357 = buffer.data(skh + 357);
    const auto *skh_360 = buffer.data(skh + 360);
    const auto *skh_362 = buffer.data(skh + 362);
    const auto *skh_363 = buffer.data(skh + 363);
    const auto *skh_366 = buffer.data(skh + 366);
    const auto *skh_367 = buffer.data(skh + 367);
    const auto *skh_369 = buffer.data(skh + 369);
    const auto *skh_371 = buffer.data(skh + 371);
    const auto *skh_372 = buffer.data(skh + 372);
    const auto *skh_373 = buffer.data(skh + 373);
    const auto *skh_374 = buffer.data(skh + 374);
    const auto *skh_375 = buffer.data(skh + 375);
    const auto *skh_376 = buffer.data(skh + 376);
    const auto *skh_377 = buffer.data(skh + 377);
    const auto *skh_378 = buffer.data(skh + 378);
    const auto *skh_381 = buffer.data(skh + 381);
    const auto *skh_383 = buffer.data(skh + 383);
    const auto *skh_384 = buffer.data(skh + 384);
    const auto *skh_387 = buffer.data(skh + 387);
    const auto *skh_388 = buffer.data(skh + 388);
    const auto *skh_390 = buffer.data(skh + 390);
    const auto *skh_392 = buffer.data(skh + 392);
    const auto *skh_393 = buffer.data(skh + 393);
    const auto *skh_394 = buffer.data(skh + 394);
    const auto *skh_395 = buffer.data(skh + 395);
    const auto *skh_396 = buffer.data(skh + 396);
    const auto *skh_397 = buffer.data(skh + 397);
    const auto *skh_398 = buffer.data(skh + 398);
    const auto *skh_414 = buffer.data(skh + 414);
    const auto *skh_415 = buffer.data(skh + 415);
    const auto *skh_416 = buffer.data(skh + 416);
    const auto *skh_417 = buffer.data(skh + 417);
    const auto *skh_418 = buffer.data(skh + 418);
    const auto *skh_419 = buffer.data(skh + 419);
    const auto *skh_420 = buffer.data(skh + 420);
    const auto *skh_423 = buffer.data(skh + 423);
    const auto *skh_425 = buffer.data(skh + 425);
    const auto *skh_426 = buffer.data(skh + 426);
    const auto *skh_429 = buffer.data(skh + 429);
    const auto *skh_430 = buffer.data(skh + 430);
    const auto *skh_432 = buffer.data(skh + 432);
    const auto *skh_434 = buffer.data(skh + 434);

    const auto *ski1_301 = buffer.data(ski1 + 301);
    const auto *ski1_392 = buffer.data(ski1 + 392);
    const auto *ski1_395 = buffer.data(ski1 + 395);
    const auto *ski1_397 = buffer.data(ski1 + 397);
    const auto *ski1_398 = buffer.data(ski1 + 398);
    const auto *ski1_401 = buffer.data(ski1 + 401);
    const auto *ski1_402 = buffer.data(ski1 + 402);
    const auto *ski1_404 = buffer.data(ski1 + 404);
    const auto *ski1_406 = buffer.data(ski1 + 406);
    const auto *ski1_419 = buffer.data(ski1 + 419);

    const auto *slg0_252 = buffer.data(slg0 + 252);
    const auto *slg0_253 = buffer.data(slg0 + 253);
    const auto *slg0_254 = buffer.data(slg0 + 254);
    const auto *slg0_255 = buffer.data(slg0 + 255);
    const auto *slg0_258 = buffer.data(slg0 + 258);
    const auto *slg0_260 = buffer.data(slg0 + 260);
    const auto *slg0_261 = buffer.data(slg0 + 261);
    const auto *slg0_264 = buffer.data(slg0 + 264);
    const auto *slg0_265 = buffer.data(slg0 + 265);
    const auto *slg0_267 = buffer.data(slg0 + 267);
    const auto *slg0_268 = buffer.data(slg0 + 268);
    const auto *slg0_269 = buffer.data(slg0 + 269);
    const auto *slg0_270 = buffer.data(slg0 + 270);
    const auto *slg0_273 = buffer.data(slg0 + 273);
    const auto *slg0_275 = buffer.data(slg0 + 275);
    const auto *slg0_276 = buffer.data(slg0 + 276);
    const auto *slg0_279 = buffer.data(slg0 + 279);
    const auto *slg0_280 = buffer.data(slg0 + 280);
    const auto *slg0_282 = buffer.data(slg0 + 282);
    const auto *slg0_283 = buffer.data(slg0 + 283);
    const auto *slg0_284 = buffer.data(slg0 + 284);
    const auto *slg0_295 = buffer.data(slg0 + 295);
    const auto *slg0_297 = buffer.data(slg0 + 297);
    const auto *slg0_298 = buffer.data(slg0 + 298);
    const auto *slg0_299 = buffer.data(slg0 + 299);
    const auto *slg0_300 = buffer.data(slg0 + 300);
    const auto *slg0_303 = buffer.data(slg0 + 303);
    const auto *slg0_305 = buffer.data(slg0 + 305);
    const auto *slg0_306 = buffer.data(slg0 + 306);
    const auto *slg0_309 = buffer.data(slg0 + 309);
    const auto *slg0_310 = buffer.data(slg0 + 310);
    const auto *slg0_312 = buffer.data(slg0 + 312);
    const auto *slg0_314 = buffer.data(slg0 + 314);

    const auto *slg1_252 = buffer.data(slg1 + 252);
    const auto *slg1_253 = buffer.data(slg1 + 253);
    const auto *slg1_254 = buffer.data(slg1 + 254);
    const auto *slg1_255 = buffer.data(slg1 + 255);
    const auto *slg1_258 = buffer.data(slg1 + 258);
    const auto *slg1_260 = buffer.data(slg1 + 260);
    const auto *slg1_261 = buffer.data(slg1 + 261);
    const auto *slg1_264 = buffer.data(slg1 + 264);
    const auto *slg1_265 = buffer.data(slg1 + 265);
    const auto *slg1_267 = buffer.data(slg1 + 267);
    const auto *slg1_268 = buffer.data(slg1 + 268);
    const auto *slg1_269 = buffer.data(slg1 + 269);
    const auto *slg1_270 = buffer.data(slg1 + 270);
    const auto *slg1_273 = buffer.data(slg1 + 273);
    const auto *slg1_275 = buffer.data(slg1 + 275);
    const auto *slg1_276 = buffer.data(slg1 + 276);
    const auto *slg1_279 = buffer.data(slg1 + 279);
    const auto *slg1_280 = buffer.data(slg1 + 280);
    const auto *slg1_282 = buffer.data(slg1 + 282);
    const auto *slg1_283 = buffer.data(slg1 + 283);
    const auto *slg1_284 = buffer.data(slg1 + 284);
    const auto *slg1_295 = buffer.data(slg1 + 295);
    const auto *slg1_297 = buffer.data(slg1 + 297);
    const auto *slg1_298 = buffer.data(slg1 + 298);
    const auto *slg1_299 = buffer.data(slg1 + 299);
    const auto *slg1_300 = buffer.data(slg1 + 300);
    const auto *slg1_303 = buffer.data(slg1 + 303);
    const auto *slg1_305 = buffer.data(slg1 + 305);
    const auto *slg1_306 = buffer.data(slg1 + 306);
    const auto *slg1_309 = buffer.data(slg1 + 309);
    const auto *slg1_310 = buffer.data(slg1 + 310);
    const auto *slg1_312 = buffer.data(slg1 + 312);
    const auto *slg1_314 = buffer.data(slg1 + 314);

    const auto *slh_351 = buffer.data(slh + 351);
    const auto *slh_353 = buffer.data(slh + 353);
    const auto *slh_354 = buffer.data(slh + 354);
    const auto *slh_355 = buffer.data(slh + 355);
    const auto *slh_356 = buffer.data(slh + 356);
    const auto *slh_357 = buffer.data(slh + 357);
    const auto *slh_359 = buffer.data(slh + 359);
    const auto *slh_360 = buffer.data(slh + 360);
    const auto *slh_362 = buffer.data(slh + 362);
    const auto *slh_363 = buffer.data(slh + 363);
    const auto *slh_366 = buffer.data(slh + 366);
    const auto *slh_367 = buffer.data(slh + 367);
    const auto *slh_369 = buffer.data(slh + 369);
    const auto *slh_371 = buffer.data(slh + 371);
    const auto *slh_372 = buffer.data(slh + 372);
    const auto *slh_373 = buffer.data(slh + 373);
    const auto *slh_374 = buffer.data(slh + 374);
    const auto *slh_375 = buffer.data(slh + 375);
    const auto *slh_376 = buffer.data(slh + 376);
    const auto *slh_377 = buffer.data(slh + 377);
    const auto *slh_378 = buffer.data(slh + 378);
    const auto *slh_380 = buffer.data(slh + 380);
    const auto *slh_381 = buffer.data(slh + 381);
    const auto *slh_383 = buffer.data(slh + 383);
    const auto *slh_384 = buffer.data(slh + 384);
    const auto *slh_387 = buffer.data(slh + 387);
    const auto *slh_388 = buffer.data(slh + 388);
    const auto *slh_390 = buffer.data(slh + 390);
    const auto *slh_392 = buffer.data(slh + 392);
    const auto *slh_393 = buffer.data(slh + 393);
    const auto *slh_394 = buffer.data(slh + 394);
    const auto *slh_395 = buffer.data(slh + 395);
    const auto *slh_396 = buffer.data(slh + 396);
    const auto *slh_397 = buffer.data(slh + 397);
    const auto *slh_398 = buffer.data(slh + 398);
    const auto *slh_399 = buffer.data(slh + 399);
    const auto *slh_401 = buffer.data(slh + 401);
    const auto *slh_402 = buffer.data(slh + 402);
    const auto *slh_404 = buffer.data(slh + 404);
    const auto *slh_405 = buffer.data(slh + 405);
    const auto *slh_408 = buffer.data(slh + 408);
    const auto *slh_414 = buffer.data(slh + 414);
    const auto *slh_415 = buffer.data(slh + 415);
    const auto *slh_416 = buffer.data(slh + 416);
    const auto *slh_417 = buffer.data(slh + 417);
    const auto *slh_418 = buffer.data(slh + 418);
    const auto *slh_419 = buffer.data(slh + 419);
    const auto *slh_420 = buffer.data(slh + 420);
    const auto *slh_422 = buffer.data(slh + 422);
    const auto *slh_423 = buffer.data(slh + 423);
    const auto *slh_425 = buffer.data(slh + 425);
    const auto *slh_426 = buffer.data(slh + 426);
    const auto *slh_429 = buffer.data(slh + 429);
    const auto *slh_430 = buffer.data(slh + 430);
    const auto *slh_432 = buffer.data(slh + 432);
    const auto *slh_434 = buffer.data(slh + 434);

#pragma omp simd aligned(t_466, t_467, t_468, t_469, pb_z, pc_x, pc_z, ski0_301, skh_354, \
                         skh_355, skh_356, ski1_301, slh_354, slh_355, \
                         slh_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_466[k] = f_13 * skh_354[k]
                   + f_3 * pc_x[k] * slh_354[k];

        t_467[k] = f_13 * skh_355[k]
                   + f_3 * pc_x[k] * slh_355[k];

        t_468[k] = f_13 * skh_356[k]
                   + f_3 * pc_x[k] * slh_356[k];

        t_469[k] = pb_z[k] * ski0_301[k]
                   - f_10 * pc_z[k] * ski1_301[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, pc_y, pc_z, skh_225, skh_248, skh_249, slg0_252, \
                         slg0_253, slg1_252, slg1_253, slh_351, slh_353, \
                         slh_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = f_11 * skh_225[k]
                   + f_3 * pc_z[k] * slh_351[k];

        t_471[k] = f_14 * skh_248[k]
                   + f_4 * slg0_252[k]
                   - f_5 * slg1_252[k]
                   + f_3 * pc_y[k] * slh_353[k];

        t_472[k] = f_14 * skh_249[k]
                   + f_6 * slg0_253[k]
                   - f_7 * slg1_253[k]
                   + f_3 * pc_y[k] * slh_354[k];
    }

#pragma omp simd aligned(t_473, t_474, t_475, pc_y, pc_z, skh_230, skh_250, skh_251, slg0_254, \
                         slg1_254, slh_355, slh_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_473[k] = f_14 * skh_250[k]
                   + f_8 * slg0_254[k]
                   - f_9 * slg1_254[k]
                   + f_3 * pc_y[k] * slh_355[k];

        t_474[k] = f_14 * skh_251[k]
                   + f_3 * pc_y[k] * slh_356[k];

        t_475[k] = f_11 * skh_230[k]
                   + f_1 * slg0_254[k]
                   - f_2 * slg1_254[k]
                   + f_3 * pc_z[k] * slh_356[k];
    }

#pragma omp simd aligned(t_476, t_477, t_478, pc_x, pc_y, pc_z, skh_231, skh_252, skh_357, \
                         slg0_255, slg1_255, slh_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_476[k] = f_13 * skh_357[k]
                   + f_1 * slg0_255[k]
                   - f_2 * slg1_255[k]
                   + f_3 * pc_x[k] * slh_357[k];

        t_477[k] = f_13 * skh_252[k]
                   + f_3 * pc_y[k] * slh_357[k];

        t_478[k] = f_12 * skh_231[k]
                   + f_3 * pc_z[k] * slh_357[k];
    }

#pragma omp simd aligned(t_479, t_480, t_481, pc_x, pc_y, skh_254, skh_360, skh_362, slg0_258, \
                         slg0_260, slg1_258, slg1_260, slh_359, slh_360, \
                         slh_362 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_479[k] = f_13 * skh_360[k]
                   + f_4 * slg0_258[k]
                   - f_5 * slg1_258[k]
                   + f_3 * pc_x[k] * slh_360[k];

        t_480[k] = f_13 * skh_254[k]
                   + f_3 * pc_y[k] * slh_359[k];

        t_481[k] = f_13 * skh_362[k]
                   + f_4 * slg0_260[k]
                   - f_5 * slg1_260[k]
                   + f_3 * pc_x[k] * slh_362[k];
    }

#pragma omp simd aligned(t_482, t_483, t_484, pc_x, pc_y, pc_z, skh_234, skh_257, skh_363, \
                         slg0_261, slg1_261, slh_360, slh_362, \
                         slh_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_482[k] = f_13 * skh_363[k]
                   + f_6 * slg0_261[k]
                   - f_7 * slg1_261[k]
                   + f_3 * pc_x[k] * slh_363[k];

        t_483[k] = f_12 * skh_234[k]
                   + f_3 * pc_z[k] * slh_360[k];

        t_484[k] = f_13 * skh_257[k]
                   + f_3 * pc_y[k] * slh_362[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, pc_x, pc_z, skh_237, skh_366, skh_367, slg0_264, \
                         slg0_265, slg1_264, slg1_265, slh_363, slh_366, \
                         slh_367 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = f_13 * skh_366[k]
                   + f_6 * slg0_264[k]
                   - f_7 * slg1_264[k]
                   + f_3 * pc_x[k] * slh_366[k];

        t_486[k] = f_13 * skh_367[k]
                   + f_8 * slg0_265[k]
                   - f_9 * slg1_265[k]
                   + f_3 * pc_x[k] * slh_367[k];

        t_487[k] = f_12 * skh_237[k]
                   + f_3 * pc_z[k] * slh_363[k];
    }

#pragma omp simd aligned(t_488, t_489, t_490, pc_x, pc_y, skh_261, skh_369, skh_371, slg0_267, \
                         slg0_269, slg1_267, slg1_269, slh_366, slh_369, \
                         slh_371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_488[k] = f_13 * skh_369[k]
                   + f_8 * slg0_267[k]
                   - f_9 * slg1_267[k]
                   + f_3 * pc_x[k] * slh_369[k];

        t_489[k] = f_13 * skh_261[k]
                   + f_3 * pc_y[k] * slh_366[k];

        t_490[k] = f_13 * skh_371[k]
                   + f_8 * slg0_269[k]
                   - f_9 * slg1_269[k]
                   + f_3 * pc_x[k] * slh_371[k];
    }

#pragma omp simd aligned(t_491, t_492, t_493, t_494, t_495, pc_x, skh_372, skh_373, skh_374, \
                         skh_375, skh_376, slh_372, slh_373, slh_374, slh_375, \
                         slh_376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_491[k] = f_13 * skh_372[k]
                   + f_3 * pc_x[k] * slh_372[k];

        t_492[k] = f_13 * skh_373[k]
                   + f_3 * pc_x[k] * slh_373[k];

        t_493[k] = f_13 * skh_374[k]
                   + f_3 * pc_x[k] * slh_374[k];

        t_494[k] = f_13 * skh_375[k]
                   + f_3 * pc_x[k] * slh_375[k];

        t_495[k] = f_13 * skh_376[k]
                   + f_3 * pc_x[k] * slh_376[k];
    }

#pragma omp simd aligned(t_496, t_497, t_498, pc_x, pc_y, pc_z, skh_246, skh_267, skh_377, \
                         slg0_265, slg1_265, slh_372, slh_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_496[k] = f_13 * skh_377[k]
                   + f_3 * pc_x[k] * slh_377[k];

        t_497[k] = f_13 * skh_267[k]
                   + f_1 * slg0_265[k]
                   - f_2 * slg1_265[k]
                   + f_3 * pc_y[k] * slh_372[k];

        t_498[k] = f_12 * skh_246[k]
                   + f_3 * pc_z[k] * slh_372[k];
    }

#pragma omp simd aligned(t_499, t_500, t_501, pc_y, skh_269, skh_270, skh_271, slg0_267, \
                         slg0_268, slg0_269, slg1_267, slg1_268, slg1_269, slh_374, slh_375, \
                         slh_376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_499[k] = f_13 * skh_269[k]
                   + f_4 * slg0_267[k]
                   - f_5 * slg1_267[k]
                   + f_3 * pc_y[k] * slh_374[k];

        t_500[k] = f_13 * skh_270[k]
                   + f_6 * slg0_268[k]
                   - f_7 * slg1_268[k]
                   + f_3 * pc_y[k] * slh_375[k];

        t_501[k] = f_13 * skh_271[k]
                   + f_8 * slg0_269[k]
                   - f_9 * slg1_269[k]
                   + f_3 * pc_y[k] * slh_376[k];
    }

#pragma omp simd aligned(t_502, t_503, t_504, pc_x, pc_y, pc_z, skh_251, skh_272, skh_378, \
                         slg0_269, slg0_270, slg1_269, slg1_270, slh_377, \
                         slh_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_502[k] = f_13 * skh_272[k]
                   + f_3 * pc_y[k] * slh_377[k];

        t_503[k] = f_12 * skh_251[k]
                   + f_1 * slg0_269[k]
                   - f_2 * slg1_269[k]
                   + f_3 * pc_z[k] * slh_377[k];

        t_504[k] = f_13 * skh_378[k]
                   + f_1 * slg0_270[k]
                   - f_2 * slg1_270[k]
                   + f_3 * pc_x[k] * slh_378[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, pc_x, pc_y, pc_z, skh_252, skh_273, \
                         skh_275, skh_381, slg0_273, slg1_273, slh_378, slh_380, \
                         slh_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = f_12 * skh_273[k]
                   + f_3 * pc_y[k] * slh_378[k];

        t_506[k] = f_13 * skh_252[k]
                   + f_3 * pc_z[k] * slh_378[k];

        t_507[k] = f_13 * skh_381[k]
                   + f_4 * slg0_273[k]
                   - f_5 * slg1_273[k]
                   + f_3 * pc_x[k] * slh_381[k];

        t_508[k] = f_12 * skh_275[k]
                   + f_3 * pc_y[k] * slh_380[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, pc_x, pc_z, skh_255, skh_383, skh_384, slg0_275, \
                         slg0_276, slg1_275, slg1_276, slh_381, slh_383, \
                         slh_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = f_13 * skh_383[k]
                   + f_4 * slg0_275[k]
                   - f_5 * slg1_275[k]
                   + f_3 * pc_x[k] * slh_383[k];

        t_510[k] = f_13 * skh_384[k]
                   + f_6 * slg0_276[k]
                   - f_7 * slg1_276[k]
                   + f_3 * pc_x[k] * slh_384[k];

        t_511[k] = f_13 * skh_255[k]
                   + f_3 * pc_z[k] * slh_381[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, pc_x, pc_y, skh_278, skh_387, skh_388, slg0_279, \
                         slg0_280, slg1_279, slg1_280, slh_383, slh_387, \
                         slh_388 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = f_12 * skh_278[k]
                   + f_3 * pc_y[k] * slh_383[k];

        t_513[k] = f_13 * skh_387[k]
                   + f_6 * slg0_279[k]
                   - f_7 * slg1_279[k]
                   + f_3 * pc_x[k] * slh_387[k];

        t_514[k] = f_13 * skh_388[k]
                   + f_8 * slg0_280[k]
                   - f_9 * slg1_280[k]
                   + f_3 * pc_x[k] * slh_388[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, pc_x, pc_y, pc_z, skh_258, skh_282, skh_390, \
                         slg0_282, slg1_282, slh_384, slh_387, \
                         slh_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = f_13 * skh_258[k]
                   + f_3 * pc_z[k] * slh_384[k];

        t_516[k] = f_13 * skh_390[k]
                   + f_8 * slg0_282[k]
                   - f_9 * slg1_282[k]
                   + f_3 * pc_x[k] * slh_390[k];

        t_517[k] = f_12 * skh_282[k]
                   + f_3 * pc_y[k] * slh_387[k];
    }

#pragma omp simd aligned(t_518, t_519, t_520, t_521, pc_x, skh_392, skh_393, skh_394, skh_395, \
                         slg0_284, slg1_284, slh_392, slh_393, slh_394, \
                         slh_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_518[k] = f_13 * skh_392[k]
                   + f_8 * slg0_284[k]
                   - f_9 * slg1_284[k]
                   + f_3 * pc_x[k] * slh_392[k];

        t_519[k] = f_13 * skh_393[k]
                   + f_3 * pc_x[k] * slh_393[k];

        t_520[k] = f_13 * skh_394[k]
                   + f_3 * pc_x[k] * slh_394[k];

        t_521[k] = f_13 * skh_395[k]
                   + f_3 * pc_x[k] * slh_395[k];
    }

#pragma omp simd aligned(t_522, t_523, t_524, t_525, pc_x, pc_y, skh_288, skh_396, skh_397, \
                         skh_398, slg0_280, slg1_280, slh_393, slh_396, slh_397, \
                         slh_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_522[k] = f_13 * skh_396[k]
                   + f_3 * pc_x[k] * slh_396[k];

        t_523[k] = f_13 * skh_397[k]
                   + f_3 * pc_x[k] * slh_397[k];

        t_524[k] = f_13 * skh_398[k]
                   + f_3 * pc_x[k] * slh_398[k];

        t_525[k] = f_12 * skh_288[k]
                   + f_1 * slg0_280[k]
                   - f_2 * slg1_280[k]
                   + f_3 * pc_y[k] * slh_393[k];
    }

#pragma omp simd aligned(t_526, t_527, t_528, pc_y, pc_z, skh_267, skh_290, skh_291, slg0_282, \
                         slg0_283, slg1_282, slg1_283, slh_393, slh_395, \
                         slh_396 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_526[k] = f_13 * skh_267[k]
                   + f_3 * pc_z[k] * slh_393[k];

        t_527[k] = f_12 * skh_290[k]
                   + f_4 * slg0_282[k]
                   - f_5 * slg1_282[k]
                   + f_3 * pc_y[k] * slh_395[k];

        t_528[k] = f_12 * skh_291[k]
                   + f_6 * slg0_283[k]
                   - f_7 * slg1_283[k]
                   + f_3 * pc_y[k] * slh_396[k];
    }

#pragma omp simd aligned(t_529, t_530, t_531, t_532, pb_y, pc_y, pc_z, ski0_392, skh_272, \
                         skh_292, skh_293, ski1_392, slg0_284, slg1_284, slh_397, \
                         slh_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_529[k] = f_12 * skh_292[k]
                   + f_8 * slg0_284[k]
                   - f_9 * slg1_284[k]
                   + f_3 * pc_y[k] * slh_397[k];

        t_530[k] = f_12 * skh_293[k]
                   + f_3 * pc_y[k] * slh_398[k];

        t_531[k] = f_13 * skh_272[k]
                   + f_1 * slg0_284[k]
                   - f_2 * slg1_284[k]
                   + f_3 * pc_z[k] * slh_398[k];

        t_532[k] = pb_y[k] * ski0_392[k]
                   - f_10 * pc_y[k] * ski1_392[k];
    }

#pragma omp simd aligned(t_533, t_534, t_535, t_536, pb_y, pc_y, pc_z, ski0_395, skh_273, \
                         skh_294, skh_295, skh_296, ski1_395, slh_399, \
                         slh_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_533[k] = f_11 * skh_294[k]
                   + f_3 * pc_y[k] * slh_399[k];

        t_534[k] = f_14 * skh_273[k]
                   + f_3 * pc_z[k] * slh_399[k];

        t_535[k] = pb_y[k] * ski0_395[k]
                   + f_12 * skh_295[k]
                   - f_10 * pc_y[k] * ski1_395[k];

        t_536[k] = f_11 * skh_296[k]
                   + f_3 * pc_y[k] * slh_401[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, t_540, pb_y, pc_y, pc_z, ski0_397, ski0_398, \
                         skh_276, skh_297, skh_299, ski1_397, ski1_398, slh_402, \
                         slh_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = pb_y[k] * ski0_397[k]
                   - f_10 * pc_y[k] * ski1_397[k];

        t_538[k] = pb_y[k] * ski0_398[k]
                   + f_13 * skh_297[k]
                   - f_10 * pc_y[k] * ski1_398[k];

        t_539[k] = f_14 * skh_276[k]
                   + f_3 * pc_z[k] * slh_402[k];

        t_540[k] = f_11 * skh_299[k]
                   + f_3 * pc_y[k] * slh_404[k];
    }

#pragma omp simd aligned(t_541, t_542, t_543, pb_y, pc_y, pc_z, ski0_401, ski0_402, skh_279, \
                         skh_300, ski1_401, ski1_402, slh_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_541[k] = pb_y[k] * ski0_401[k]
                   - f_10 * pc_y[k] * ski1_401[k];

        t_542[k] = pb_y[k] * ski0_402[k]
                   + f_14 * skh_300[k]
                   - f_10 * pc_y[k] * ski1_402[k];

        t_543[k] = f_14 * skh_279[k]
                   + f_3 * pc_z[k] * slh_405[k];
    }

#pragma omp simd aligned(t_544, t_545, t_546, t_547, pb_y, pc_x, pc_y, ski0_404, ski0_406, \
                         skh_302, skh_303, skh_414, ski1_404, ski1_406, slh_408, \
                         slh_414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_544[k] = pb_y[k] * ski0_404[k]
                   + f_12 * skh_302[k]
                   - f_10 * pc_y[k] * ski1_404[k];

        t_545[k] = f_11 * skh_303[k]
                   + f_3 * pc_y[k] * slh_408[k];

        t_546[k] = pb_y[k] * ski0_406[k]
                   - f_10 * pc_y[k] * ski1_406[k];

        t_547[k] = f_13 * skh_414[k]
                   + f_3 * pc_x[k] * slh_414[k];
    }

#pragma omp simd aligned(t_548, t_549, t_550, t_551, t_552, pc_x, skh_415, skh_416, skh_417, \
                         skh_418, skh_419, slh_415, slh_416, slh_417, slh_418, \
                         slh_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_548[k] = f_13 * skh_415[k]
                   + f_3 * pc_x[k] * slh_415[k];

        t_549[k] = f_13 * skh_416[k]
                   + f_3 * pc_x[k] * slh_416[k];

        t_550[k] = f_13 * skh_417[k]
                   + f_3 * pc_x[k] * slh_417[k];

        t_551[k] = f_13 * skh_418[k]
                   + f_3 * pc_x[k] * slh_418[k];

        t_552[k] = f_13 * skh_419[k]
                   + f_3 * pc_x[k] * slh_419[k];
    }

#pragma omp simd aligned(t_553, t_554, t_555, pc_y, pc_z, skh_288, skh_309, skh_311, slg0_295, \
                         slg0_297, slg1_295, slg1_297, slh_414, \
                         slh_416 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_553[k] = f_11 * skh_309[k]
                   + f_1 * slg0_295[k]
                   - f_2 * slg1_295[k]
                   + f_3 * pc_y[k] * slh_414[k];

        t_554[k] = f_14 * skh_288[k]
                   + f_3 * pc_z[k] * slh_414[k];

        t_555[k] = f_11 * skh_311[k]
                   + f_4 * slg0_297[k]
                   - f_5 * slg1_297[k]
                   + f_3 * pc_y[k] * slh_416[k];
    }

#pragma omp simd aligned(t_556, t_557, t_558, pc_y, skh_312, skh_313, skh_314, slg0_298, \
                         slg0_299, slg1_298, slg1_299, slh_417, slh_418, \
                         slh_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_556[k] = f_11 * skh_312[k]
                   + f_6 * slg0_298[k]
                   - f_7 * slg1_298[k]
                   + f_3 * pc_y[k] * slh_417[k];

        t_557[k] = f_11 * skh_313[k]
                   + f_8 * slg0_299[k]
                   - f_9 * slg1_299[k]
                   + f_3 * pc_y[k] * slh_418[k];

        t_558[k] = f_11 * skh_314[k]
                   + f_3 * pc_y[k] * slh_419[k];
    }

#pragma omp simd aligned(t_559, t_560, t_561, t_562, pb_y, pc_x, pc_y, pc_z, ski0_419, \
                         skh_294, skh_420, ski1_419, slg0_300, slg1_300, \
                         slh_420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_559[k] = pb_y[k] * ski0_419[k]
                   - f_10 * pc_y[k] * ski1_419[k];

        t_560[k] = f_13 * skh_420[k]
                   + f_1 * slg0_300[k]
                   - f_2 * slg1_300[k]
                   + f_3 * pc_x[k] * slh_420[k];

        t_561[k] = f_3 * pc_y[k] * slh_420[k];

        t_562[k] = f_17 * skh_294[k]
                   + f_3 * pc_z[k] * slh_420[k];
    }

#pragma omp simd aligned(t_563, t_564, t_565, pc_x, pc_y, skh_423, skh_425, slg0_303, \
                         slg0_305, slg1_303, slg1_305, slh_422, slh_423, \
                         slh_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_563[k] = f_13 * skh_423[k]
                   + f_4 * slg0_303[k]
                   - f_5 * slg1_303[k]
                   + f_3 * pc_x[k] * slh_423[k];

        t_564[k] = f_3 * pc_y[k] * slh_422[k];

        t_565[k] = f_13 * skh_425[k]
                   + f_4 * slg0_305[k]
                   - f_5 * slg1_305[k]
                   + f_3 * pc_x[k] * slh_425[k];
    }

#pragma omp simd aligned(t_566, t_567, t_568, pc_x, pc_y, pc_z, skh_297, skh_426, slg0_306, \
                         slg1_306, slh_423, slh_425, slh_426 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_566[k] = f_13 * skh_426[k]
                   + f_6 * slg0_306[k]
                   - f_7 * slg1_306[k]
                   + f_3 * pc_x[k] * slh_426[k];

        t_567[k] = f_17 * skh_297[k]
                   + f_3 * pc_z[k] * slh_423[k];

        t_568[k] = f_3 * pc_y[k] * slh_425[k];
    }

#pragma omp simd aligned(t_569, t_570, t_571, pc_x, pc_z, skh_300, skh_429, skh_430, slg0_309, \
                         slg0_310, slg1_309, slg1_310, slh_426, slh_429, \
                         slh_430 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_569[k] = f_13 * skh_429[k]
                   + f_6 * slg0_309[k]
                   - f_7 * slg1_309[k]
                   + f_3 * pc_x[k] * slh_429[k];

        t_570[k] = f_13 * skh_430[k]
                   + f_8 * slg0_310[k]
                   - f_9 * slg1_310[k]
                   + f_3 * pc_x[k] * slh_430[k];

        t_571[k] = f_17 * skh_300[k]
                   + f_3 * pc_z[k] * slh_426[k];
    }

#pragma omp simd aligned(t_572, t_573, t_574, pc_x, pc_y, skh_432, skh_434, slg0_312, \
                         slg0_314, slg1_312, slg1_314, slh_429, slh_432, \
                         slh_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_572[k] = f_13 * skh_432[k]
                   + f_8 * slg0_312[k]
                   - f_9 * slg1_312[k]
                   + f_3 * pc_x[k] * slh_432[k];

        t_573[k] = f_3 * pc_y[k] * slh_429[k];

        t_574[k] = f_13 * skh_434[k]
                   + f_8 * slg0_314[k]
                   - f_9 * slg1_314[k]
                   + f_3 * pc_x[k] * slh_434[k];
    }
}

static auto
compute_prim_sli_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t ski0,
                                                          const size_t skh, const size_t ski1,
                                                          const size_t slg0, const size_t slg1,
                                                          const size_t slh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.5 / gamma;
    const auto f_5 = 1.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 0.5 / gamma;
    const auto f_9 = 0.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_16 = 3.0 / q;
    const auto f_17 = 2.5 / q;

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

    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ski0_420 = buffer.data(ski0 + 420);
    const auto *ski0_423 = buffer.data(ski0 + 423);
    const auto *ski0_426 = buffer.data(ski0 + 426);
    const auto *ski0_430 = buffer.data(ski0 + 430);
    const auto *ski0_432 = buffer.data(ski0 + 432);
    const auto *ski0_441 = buffer.data(ski0 + 441);

    const auto *skh_309 = buffer.data(skh + 309);
    const auto *skh_314 = buffer.data(skh + 314);
    const auto *skh_315 = buffer.data(skh + 315);
    const auto *skh_317 = buffer.data(skh + 317);
    const auto *skh_318 = buffer.data(skh + 318);
    const auto *skh_320 = buffer.data(skh + 320);
    const auto *skh_321 = buffer.data(skh + 321);
    const auto *skh_322 = buffer.data(skh + 322);
    const auto *skh_324 = buffer.data(skh + 324);
    const auto *skh_330 = buffer.data(skh + 330);
    const auto *skh_332 = buffer.data(skh + 332);
    const auto *skh_333 = buffer.data(skh + 333);
    const auto *skh_334 = buffer.data(skh + 334);
    const auto *skh_335 = buffer.data(skh + 335);
    const auto *skh_336 = buffer.data(skh + 336);
    const auto *skh_338 = buffer.data(skh + 338);
    const auto *skh_339 = buffer.data(skh + 339);
    const auto *skh_341 = buffer.data(skh + 341);
    const auto *skh_342 = buffer.data(skh + 342);
    const auto *skh_345 = buffer.data(skh + 345);
    const auto *skh_351 = buffer.data(skh + 351);
    const auto *skh_353 = buffer.data(skh + 353);
    const auto *skh_354 = buffer.data(skh + 354);
    const auto *skh_355 = buffer.data(skh + 355);
    const auto *skh_356 = buffer.data(skh + 356);
    const auto *skh_357 = buffer.data(skh + 357);
    const auto *skh_359 = buffer.data(skh + 359);
    const auto *skh_360 = buffer.data(skh + 360);
    const auto *skh_362 = buffer.data(skh + 362);
    const auto *skh_363 = buffer.data(skh + 363);
    const auto *skh_366 = buffer.data(skh + 366);
    const auto *skh_372 = buffer.data(skh + 372);
    const auto *skh_374 = buffer.data(skh + 374);
    const auto *skh_375 = buffer.data(skh + 375);
    const auto *skh_376 = buffer.data(skh + 376);
    const auto *skh_377 = buffer.data(skh + 377);
    const auto *skh_378 = buffer.data(skh + 378);
    const auto *skh_380 = buffer.data(skh + 380);
    const auto *skh_383 = buffer.data(skh + 383);
    const auto *skh_387 = buffer.data(skh + 387);
    const auto *skh_435 = buffer.data(skh + 435);
    const auto *skh_436 = buffer.data(skh + 436);
    const auto *skh_437 = buffer.data(skh + 437);
    const auto *skh_438 = buffer.data(skh + 438);
    const auto *skh_439 = buffer.data(skh + 439);
    const auto *skh_440 = buffer.data(skh + 440);
    const auto *skh_441 = buffer.data(skh + 441);
    const auto *skh_444 = buffer.data(skh + 444);
    const auto *skh_446 = buffer.data(skh + 446);
    const auto *skh_447 = buffer.data(skh + 447);
    const auto *skh_450 = buffer.data(skh + 450);
    const auto *skh_451 = buffer.data(skh + 451);
    const auto *skh_453 = buffer.data(skh + 453);
    const auto *skh_455 = buffer.data(skh + 455);
    const auto *skh_456 = buffer.data(skh + 456);
    const auto *skh_457 = buffer.data(skh + 457);
    const auto *skh_458 = buffer.data(skh + 458);
    const auto *skh_459 = buffer.data(skh + 459);
    const auto *skh_460 = buffer.data(skh + 460);
    const auto *skh_461 = buffer.data(skh + 461);
    const auto *skh_467 = buffer.data(skh + 467);
    const auto *skh_471 = buffer.data(skh + 471);
    const auto *skh_476 = buffer.data(skh + 476);
    const auto *skh_477 = buffer.data(skh + 477);
    const auto *skh_478 = buffer.data(skh + 478);
    const auto *skh_479 = buffer.data(skh + 479);
    const auto *skh_480 = buffer.data(skh + 480);
    const auto *skh_481 = buffer.data(skh + 481);
    const auto *skh_482 = buffer.data(skh + 482);
    const auto *skh_483 = buffer.data(skh + 483);
    const auto *skh_486 = buffer.data(skh + 486);
    const auto *skh_488 = buffer.data(skh + 488);
    const auto *skh_489 = buffer.data(skh + 489);
    const auto *skh_492 = buffer.data(skh + 492);
    const auto *skh_493 = buffer.data(skh + 493);
    const auto *skh_495 = buffer.data(skh + 495);
    const auto *skh_497 = buffer.data(skh + 497);
    const auto *skh_498 = buffer.data(skh + 498);
    const auto *skh_499 = buffer.data(skh + 499);
    const auto *skh_500 = buffer.data(skh + 500);
    const auto *skh_501 = buffer.data(skh + 501);
    const auto *skh_502 = buffer.data(skh + 502);
    const auto *skh_503 = buffer.data(skh + 503);
    const auto *skh_504 = buffer.data(skh + 504);
    const auto *skh_507 = buffer.data(skh + 507);
    const auto *skh_509 = buffer.data(skh + 509);
    const auto *skh_510 = buffer.data(skh + 510);
    const auto *skh_513 = buffer.data(skh + 513);
    const auto *skh_514 = buffer.data(skh + 514);
    const auto *skh_516 = buffer.data(skh + 516);

    const auto *ski1_420 = buffer.data(ski1 + 420);
    const auto *ski1_423 = buffer.data(ski1 + 423);
    const auto *ski1_426 = buffer.data(ski1 + 426);
    const auto *ski1_430 = buffer.data(ski1 + 430);
    const auto *ski1_432 = buffer.data(ski1 + 432);
    const auto *ski1_441 = buffer.data(ski1 + 441);

    const auto *slg0_310 = buffer.data(slg0 + 310);
    const auto *slg0_312 = buffer.data(slg0 + 312);
    const auto *slg0_313 = buffer.data(slg0 + 313);
    const auto *slg0_314 = buffer.data(slg0 + 314);
    const auto *slg0_315 = buffer.data(slg0 + 315);
    const auto *slg0_318 = buffer.data(slg0 + 318);
    const auto *slg0_320 = buffer.data(slg0 + 320);
    const auto *slg0_321 = buffer.data(slg0 + 321);
    const auto *slg0_324 = buffer.data(slg0 + 324);
    const auto *slg0_325 = buffer.data(slg0 + 325);
    const auto *slg0_327 = buffer.data(slg0 + 327);
    const auto *slg0_328 = buffer.data(slg0 + 328);
    const auto *slg0_329 = buffer.data(slg0 + 329);
    const auto *slg0_335 = buffer.data(slg0 + 335);
    const auto *slg0_339 = buffer.data(slg0 + 339);
    const auto *slg0_342 = buffer.data(slg0 + 342);
    const auto *slg0_343 = buffer.data(slg0 + 343);
    const auto *slg0_344 = buffer.data(slg0 + 344);
    const auto *slg0_345 = buffer.data(slg0 + 345);
    const auto *slg0_348 = buffer.data(slg0 + 348);
    const auto *slg0_350 = buffer.data(slg0 + 350);
    const auto *slg0_351 = buffer.data(slg0 + 351);
    const auto *slg0_354 = buffer.data(slg0 + 354);
    const auto *slg0_355 = buffer.data(slg0 + 355);
    const auto *slg0_357 = buffer.data(slg0 + 357);
    const auto *slg0_358 = buffer.data(slg0 + 358);
    const auto *slg0_359 = buffer.data(slg0 + 359);
    const auto *slg0_360 = buffer.data(slg0 + 360);
    const auto *slg0_363 = buffer.data(slg0 + 363);
    const auto *slg0_365 = buffer.data(slg0 + 365);
    const auto *slg0_366 = buffer.data(slg0 + 366);
    const auto *slg0_369 = buffer.data(slg0 + 369);
    const auto *slg0_370 = buffer.data(slg0 + 370);
    const auto *slg0_372 = buffer.data(slg0 + 372);

    const auto *slg1_310 = buffer.data(slg1 + 310);
    const auto *slg1_312 = buffer.data(slg1 + 312);
    const auto *slg1_313 = buffer.data(slg1 + 313);
    const auto *slg1_314 = buffer.data(slg1 + 314);
    const auto *slg1_315 = buffer.data(slg1 + 315);
    const auto *slg1_318 = buffer.data(slg1 + 318);
    const auto *slg1_320 = buffer.data(slg1 + 320);
    const auto *slg1_321 = buffer.data(slg1 + 321);
    const auto *slg1_324 = buffer.data(slg1 + 324);
    const auto *slg1_325 = buffer.data(slg1 + 325);
    const auto *slg1_327 = buffer.data(slg1 + 327);
    const auto *slg1_328 = buffer.data(slg1 + 328);
    const auto *slg1_329 = buffer.data(slg1 + 329);
    const auto *slg1_335 = buffer.data(slg1 + 335);
    const auto *slg1_339 = buffer.data(slg1 + 339);
    const auto *slg1_342 = buffer.data(slg1 + 342);
    const auto *slg1_343 = buffer.data(slg1 + 343);
    const auto *slg1_344 = buffer.data(slg1 + 344);
    const auto *slg1_345 = buffer.data(slg1 + 345);
    const auto *slg1_348 = buffer.data(slg1 + 348);
    const auto *slg1_350 = buffer.data(slg1 + 350);
    const auto *slg1_351 = buffer.data(slg1 + 351);
    const auto *slg1_354 = buffer.data(slg1 + 354);
    const auto *slg1_355 = buffer.data(slg1 + 355);
    const auto *slg1_357 = buffer.data(slg1 + 357);
    const auto *slg1_358 = buffer.data(slg1 + 358);
    const auto *slg1_359 = buffer.data(slg1 + 359);
    const auto *slg1_360 = buffer.data(slg1 + 360);
    const auto *slg1_363 = buffer.data(slg1 + 363);
    const auto *slg1_365 = buffer.data(slg1 + 365);
    const auto *slg1_366 = buffer.data(slg1 + 366);
    const auto *slg1_369 = buffer.data(slg1 + 369);
    const auto *slg1_370 = buffer.data(slg1 + 370);
    const auto *slg1_372 = buffer.data(slg1 + 372);

    const auto *slh_435 = buffer.data(slh + 435);
    const auto *slh_436 = buffer.data(slh + 436);
    const auto *slh_437 = buffer.data(slh + 437);
    const auto *slh_438 = buffer.data(slh + 438);
    const auto *slh_439 = buffer.data(slh + 439);
    const auto *slh_440 = buffer.data(slh + 440);
    const auto *slh_441 = buffer.data(slh + 441);
    const auto *slh_443 = buffer.data(slh + 443);
    const auto *slh_444 = buffer.data(slh + 444);
    const auto *slh_446 = buffer.data(slh + 446);
    const auto *slh_447 = buffer.data(slh + 447);
    const auto *slh_450 = buffer.data(slh + 450);
    const auto *slh_451 = buffer.data(slh + 451);
    const auto *slh_453 = buffer.data(slh + 453);
    const auto *slh_455 = buffer.data(slh + 455);
    const auto *slh_456 = buffer.data(slh + 456);
    const auto *slh_457 = buffer.data(slh + 457);
    const auto *slh_458 = buffer.data(slh + 458);
    const auto *slh_459 = buffer.data(slh + 459);
    const auto *slh_460 = buffer.data(slh + 460);
    const auto *slh_461 = buffer.data(slh + 461);
    const auto *slh_462 = buffer.data(slh + 462);
    const auto *slh_464 = buffer.data(slh + 464);
    const auto *slh_465 = buffer.data(slh + 465);
    const auto *slh_467 = buffer.data(slh + 467);
    const auto *slh_468 = buffer.data(slh + 468);
    const auto *slh_471 = buffer.data(slh + 471);
    const auto *slh_476 = buffer.data(slh + 476);
    const auto *slh_477 = buffer.data(slh + 477);
    const auto *slh_478 = buffer.data(slh + 478);
    const auto *slh_479 = buffer.data(slh + 479);
    const auto *slh_480 = buffer.data(slh + 480);
    const auto *slh_481 = buffer.data(slh + 481);
    const auto *slh_482 = buffer.data(slh + 482);
    const auto *slh_483 = buffer.data(slh + 483);
    const auto *slh_485 = buffer.data(slh + 485);
    const auto *slh_486 = buffer.data(slh + 486);
    const auto *slh_488 = buffer.data(slh + 488);
    const auto *slh_489 = buffer.data(slh + 489);
    const auto *slh_492 = buffer.data(slh + 492);
    const auto *slh_493 = buffer.data(slh + 493);
    const auto *slh_495 = buffer.data(slh + 495);
    const auto *slh_497 = buffer.data(slh + 497);
    const auto *slh_498 = buffer.data(slh + 498);
    const auto *slh_499 = buffer.data(slh + 499);
    const auto *slh_500 = buffer.data(slh + 500);
    const auto *slh_501 = buffer.data(slh + 501);
    const auto *slh_502 = buffer.data(slh + 502);
    const auto *slh_503 = buffer.data(slh + 503);
    const auto *slh_504 = buffer.data(slh + 504);
    const auto *slh_506 = buffer.data(slh + 506);
    const auto *slh_507 = buffer.data(slh + 507);
    const auto *slh_509 = buffer.data(slh + 509);
    const auto *slh_510 = buffer.data(slh + 510);
    const auto *slh_513 = buffer.data(slh + 513);
    const auto *slh_514 = buffer.data(slh + 514);
    const auto *slh_516 = buffer.data(slh + 516);

#pragma omp simd aligned(t_575, t_576, t_577, t_578, t_579, pc_x, skh_435, skh_436, skh_437, \
                         skh_438, skh_439, slh_435, slh_436, slh_437, slh_438, \
                         slh_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_575[k] = f_13 * skh_435[k]
                   + f_3 * pc_x[k] * slh_435[k];

        t_576[k] = f_13 * skh_436[k]
                   + f_3 * pc_x[k] * slh_436[k];

        t_577[k] = f_13 * skh_437[k]
                   + f_3 * pc_x[k] * slh_437[k];

        t_578[k] = f_13 * skh_438[k]
                   + f_3 * pc_x[k] * slh_438[k];

        t_579[k] = f_13 * skh_439[k]
                   + f_3 * pc_x[k] * slh_439[k];
    }

#pragma omp simd aligned(t_580, t_581, t_582, t_583, pc_x, pc_y, pc_z, skh_309, skh_440, \
                         slg0_310, slg0_312, slg1_310, slg1_312, slh_435, slh_437, \
                         slh_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_580[k] = f_13 * skh_440[k]
                   + f_3 * pc_x[k] * slh_440[k];

        t_581[k] = f_1 * slg0_310[k]
                   - f_2 * slg1_310[k]
                   + f_3 * pc_y[k] * slh_435[k];

        t_582[k] = f_17 * skh_309[k]
                   + f_3 * pc_z[k] * slh_435[k];

        t_583[k] = f_4 * slg0_312[k]
                   - f_5 * slg1_312[k]
                   + f_3 * pc_y[k] * slh_437[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, t_587, pc_y, pc_z, skh_314, slg0_313, slg0_314, \
                         slg1_313, slg1_314, slh_438, slh_439, \
                         slh_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = f_6 * slg0_313[k]
                   - f_7 * slg1_313[k]
                   + f_3 * pc_y[k] * slh_438[k];

        t_585[k] = f_8 * slg0_314[k]
                   - f_9 * slg1_314[k]
                   + f_3 * pc_y[k] * slh_439[k];

        t_586[k] = f_3 * pc_y[k] * slh_440[k];

        t_587[k] = f_17 * skh_314[k]
                   + f_1 * slg0_314[k]
                   - f_2 * slg1_314[k]
                   + f_3 * pc_z[k] * slh_440[k];
    }

#pragma omp simd aligned(t_588, t_589, t_590, t_591, pc_x, pc_y, pc_z, skh_315, skh_441, \
                         skh_444, slg0_315, slg0_318, slg1_315, slg1_318, slh_441, \
                         slh_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_588[k] = f_12 * skh_441[k]
                   + f_1 * slg0_315[k]
                   - f_2 * slg1_315[k]
                   + f_3 * pc_x[k] * slh_441[k];

        t_589[k] = f_16 * skh_315[k]
                   + f_3 * pc_y[k] * slh_441[k];

        t_590[k] = f_3 * pc_z[k] * slh_441[k];

        t_591[k] = f_12 * skh_444[k]
                   + f_4 * slg0_318[k]
                   - f_5 * slg1_318[k]
                   + f_3 * pc_x[k] * slh_444[k];
    }

#pragma omp simd aligned(t_592, t_593, t_594, pc_x, pc_y, skh_317, skh_446, skh_447, slg0_320, \
                         slg0_321, slg1_320, slg1_321, slh_443, slh_446, \
                         slh_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_592[k] = f_16 * skh_317[k]
                   + f_3 * pc_y[k] * slh_443[k];

        t_593[k] = f_12 * skh_446[k]
                   + f_4 * slg0_320[k]
                   - f_5 * slg1_320[k]
                   + f_3 * pc_x[k] * slh_446[k];

        t_594[k] = f_12 * skh_447[k]
                   + f_6 * slg0_321[k]
                   - f_7 * slg1_321[k]
                   + f_3 * pc_x[k] * slh_447[k];
    }

#pragma omp simd aligned(t_595, t_596, t_597, pc_x, pc_y, pc_z, skh_320, skh_450, slg0_324, \
                         slg1_324, slh_444, slh_446, slh_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_595[k] = f_3 * pc_z[k] * slh_444[k];

        t_596[k] = f_16 * skh_320[k]
                   + f_3 * pc_y[k] * slh_446[k];

        t_597[k] = f_12 * skh_450[k]
                   + f_6 * slg0_324[k]
                   - f_7 * slg1_324[k]
                   + f_3 * pc_x[k] * slh_450[k];
    }

#pragma omp simd aligned(t_598, t_599, t_600, pc_x, pc_z, skh_451, skh_453, slg0_325, \
                         slg0_327, slg1_325, slg1_327, slh_447, slh_451, \
                         slh_453 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_598[k] = f_12 * skh_451[k]
                   + f_8 * slg0_325[k]
                   - f_9 * slg1_325[k]
                   + f_3 * pc_x[k] * slh_451[k];

        t_599[k] = f_3 * pc_z[k] * slh_447[k];

        t_600[k] = f_12 * skh_453[k]
                   + f_8 * slg0_327[k]
                   - f_9 * slg1_327[k]
                   + f_3 * pc_x[k] * slh_453[k];
    }

#pragma omp simd aligned(t_601, t_602, t_603, t_604, pc_x, pc_y, skh_324, skh_455, skh_456, \
                         skh_457, slg0_329, slg1_329, slh_450, slh_455, slh_456, \
                         slh_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_601[k] = f_16 * skh_324[k]
                   + f_3 * pc_y[k] * slh_450[k];

        t_602[k] = f_12 * skh_455[k]
                   + f_8 * slg0_329[k]
                   - f_9 * slg1_329[k]
                   + f_3 * pc_x[k] * slh_455[k];

        t_603[k] = f_12 * skh_456[k]
                   + f_3 * pc_x[k] * slh_456[k];

        t_604[k] = f_12 * skh_457[k]
                   + f_3 * pc_x[k] * slh_457[k];
    }

#pragma omp simd aligned(t_605, t_606, t_607, t_608, pc_x, skh_458, skh_459, skh_460, skh_461, \
                         slh_458, slh_459, slh_460, slh_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_605[k] = f_12 * skh_458[k]
                   + f_3 * pc_x[k] * slh_458[k];

        t_606[k] = f_12 * skh_459[k]
                   + f_3 * pc_x[k] * slh_459[k];

        t_607[k] = f_12 * skh_460[k]
                   + f_3 * pc_x[k] * slh_460[k];

        t_608[k] = f_12 * skh_461[k]
                   + f_3 * pc_x[k] * slh_461[k];
    }

#pragma omp simd aligned(t_609, t_610, t_611, pc_y, pc_z, skh_330, skh_332, slg0_325, \
                         slg0_327, slg1_325, slg1_327, slh_456, \
                         slh_458 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_609[k] = f_16 * skh_330[k]
                   + f_1 * slg0_325[k]
                   - f_2 * slg1_325[k]
                   + f_3 * pc_y[k] * slh_456[k];

        t_610[k] = f_3 * pc_z[k] * slh_456[k];

        t_611[k] = f_16 * skh_332[k]
                   + f_4 * slg0_327[k]
                   - f_5 * slg1_327[k]
                   + f_3 * pc_y[k] * slh_458[k];
    }

#pragma omp simd aligned(t_612, t_613, t_614, t_615, pc_y, pc_z, skh_333, skh_334, skh_335, \
                         slg0_328, slg0_329, slg1_328, slg1_329, slh_459, slh_460, \
                         slh_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_612[k] = f_16 * skh_333[k]
                   + f_6 * slg0_328[k]
                   - f_7 * slg1_328[k]
                   + f_3 * pc_y[k] * slh_459[k];

        t_613[k] = f_16 * skh_334[k]
                   + f_8 * slg0_329[k]
                   - f_9 * slg1_329[k]
                   + f_3 * pc_y[k] * slh_460[k];

        t_614[k] = f_16 * skh_335[k]
                   + f_3 * pc_y[k] * slh_461[k];

        t_615[k] = f_1 * slg0_329[k]
                   - f_2 * slg1_329[k]
                   + f_3 * pc_z[k] * slh_461[k];
    }

#pragma omp simd aligned(t_616, t_617, t_618, t_619, pb_z, pc_y, pc_z, ski0_420, ski0_423, \
                         skh_315, skh_336, ski1_420, ski1_423, \
                         slh_462 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_616[k] = pb_z[k] * ski0_420[k]
                   - f_10 * pc_z[k] * ski1_420[k];

        t_617[k] = f_17 * skh_336[k]
                   + f_3 * pc_y[k] * slh_462[k];

        t_618[k] = f_11 * skh_315[k]
                   + f_3 * pc_z[k] * slh_462[k];

        t_619[k] = pb_z[k] * ski0_423[k]
                   - f_10 * pc_z[k] * ski1_423[k];
    }

#pragma omp simd aligned(t_620, t_621, t_622, pb_z, pc_x, pc_y, pc_z, ski0_426, skh_338, \
                         skh_467, ski1_426, slg0_335, slg1_335, slh_464, \
                         slh_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_620[k] = f_17 * skh_338[k]
                   + f_3 * pc_y[k] * slh_464[k];

        t_621[k] = f_12 * skh_467[k]
                   + f_4 * slg0_335[k]
                   - f_5 * slg1_335[k]
                   + f_3 * pc_x[k] * slh_467[k];

        t_622[k] = pb_z[k] * ski0_426[k]
                   - f_10 * pc_z[k] * ski1_426[k];
    }

#pragma omp simd aligned(t_623, t_624, t_625, pc_x, pc_y, pc_z, skh_318, skh_341, skh_471, \
                         slg0_339, slg1_339, slh_465, slh_467, \
                         slh_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_623[k] = f_11 * skh_318[k]
                   + f_3 * pc_z[k] * slh_465[k];

        t_624[k] = f_17 * skh_341[k]
                   + f_3 * pc_y[k] * slh_467[k];

        t_625[k] = f_12 * skh_471[k]
                   + f_6 * slg0_339[k]
                   - f_7 * slg1_339[k]
                   + f_3 * pc_x[k] * slh_471[k];
    }

#pragma omp simd aligned(t_626, t_627, t_628, t_629, pb_z, pc_y, pc_z, ski0_430, ski0_432, \
                         skh_321, skh_322, skh_345, ski1_430, ski1_432, slh_468, \
                         slh_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_626[k] = pb_z[k] * ski0_430[k]
                   - f_10 * pc_z[k] * ski1_430[k];

        t_627[k] = f_11 * skh_321[k]
                   + f_3 * pc_z[k] * slh_468[k];

        t_628[k] = pb_z[k] * ski0_432[k]
                   + f_12 * skh_322[k]
                   - f_10 * pc_z[k] * ski1_432[k];

        t_629[k] = f_17 * skh_345[k]
                   + f_3 * pc_y[k] * slh_471[k];
    }

#pragma omp simd aligned(t_630, t_631, t_632, t_633, pc_x, skh_476, skh_477, skh_478, skh_479, \
                         slg0_344, slg1_344, slh_476, slh_477, slh_478, \
                         slh_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_630[k] = f_12 * skh_476[k]
                   + f_8 * slg0_344[k]
                   - f_9 * slg1_344[k]
                   + f_3 * pc_x[k] * slh_476[k];

        t_631[k] = f_12 * skh_477[k]
                   + f_3 * pc_x[k] * slh_477[k];

        t_632[k] = f_12 * skh_478[k]
                   + f_3 * pc_x[k] * slh_478[k];

        t_633[k] = f_12 * skh_479[k]
                   + f_3 * pc_x[k] * slh_479[k];
    }

#pragma omp simd aligned(t_634, t_635, t_636, t_637, pb_z, pc_x, pc_z, ski0_441, skh_480, \
                         skh_481, skh_482, ski1_441, slh_480, slh_481, \
                         slh_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_634[k] = f_12 * skh_480[k]
                   + f_3 * pc_x[k] * slh_480[k];

        t_635[k] = f_12 * skh_481[k]
                   + f_3 * pc_x[k] * slh_481[k];

        t_636[k] = f_12 * skh_482[k]
                   + f_3 * pc_x[k] * slh_482[k];

        t_637[k] = pb_z[k] * ski0_441[k]
                   - f_10 * pc_z[k] * ski1_441[k];
    }

#pragma omp simd aligned(t_638, t_639, t_640, pc_y, pc_z, skh_330, skh_353, skh_354, slg0_342, \
                         slg0_343, slg1_342, slg1_343, slh_477, slh_479, \
                         slh_480 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_638[k] = f_11 * skh_330[k]
                   + f_3 * pc_z[k] * slh_477[k];

        t_639[k] = f_17 * skh_353[k]
                   + f_4 * slg0_342[k]
                   - f_5 * slg1_342[k]
                   + f_3 * pc_y[k] * slh_479[k];

        t_640[k] = f_17 * skh_354[k]
                   + f_6 * slg0_343[k]
                   - f_7 * slg1_343[k]
                   + f_3 * pc_y[k] * slh_480[k];
    }

#pragma omp simd aligned(t_641, t_642, t_643, pc_y, pc_z, skh_335, skh_355, skh_356, slg0_344, \
                         slg1_344, slh_481, slh_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_641[k] = f_17 * skh_355[k]
                   + f_8 * slg0_344[k]
                   - f_9 * slg1_344[k]
                   + f_3 * pc_y[k] * slh_481[k];

        t_642[k] = f_17 * skh_356[k]
                   + f_3 * pc_y[k] * slh_482[k];

        t_643[k] = f_11 * skh_335[k]
                   + f_1 * slg0_344[k]
                   - f_2 * slg1_344[k]
                   + f_3 * pc_z[k] * slh_482[k];
    }

#pragma omp simd aligned(t_644, t_645, t_646, pc_x, pc_y, pc_z, skh_336, skh_357, skh_483, \
                         slg0_345, slg1_345, slh_483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_644[k] = f_12 * skh_483[k]
                   + f_1 * slg0_345[k]
                   - f_2 * slg1_345[k]
                   + f_3 * pc_x[k] * slh_483[k];

        t_645[k] = f_14 * skh_357[k]
                   + f_3 * pc_y[k] * slh_483[k];

        t_646[k] = f_12 * skh_336[k]
                   + f_3 * pc_z[k] * slh_483[k];
    }

#pragma omp simd aligned(t_647, t_648, t_649, pc_x, pc_y, skh_359, skh_486, skh_488, slg0_348, \
                         slg0_350, slg1_348, slg1_350, slh_485, slh_486, \
                         slh_488 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_647[k] = f_12 * skh_486[k]
                   + f_4 * slg0_348[k]
                   - f_5 * slg1_348[k]
                   + f_3 * pc_x[k] * slh_486[k];

        t_648[k] = f_14 * skh_359[k]
                   + f_3 * pc_y[k] * slh_485[k];

        t_649[k] = f_12 * skh_488[k]
                   + f_4 * slg0_350[k]
                   - f_5 * slg1_350[k]
                   + f_3 * pc_x[k] * slh_488[k];
    }

#pragma omp simd aligned(t_650, t_651, t_652, pc_x, pc_y, pc_z, skh_339, skh_362, skh_489, \
                         slg0_351, slg1_351, slh_486, slh_488, \
                         slh_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_650[k] = f_12 * skh_489[k]
                   + f_6 * slg0_351[k]
                   - f_7 * slg1_351[k]
                   + f_3 * pc_x[k] * slh_489[k];

        t_651[k] = f_12 * skh_339[k]
                   + f_3 * pc_z[k] * slh_486[k];

        t_652[k] = f_14 * skh_362[k]
                   + f_3 * pc_y[k] * slh_488[k];
    }

#pragma omp simd aligned(t_653, t_654, t_655, pc_x, pc_z, skh_342, skh_492, skh_493, slg0_354, \
                         slg0_355, slg1_354, slg1_355, slh_489, slh_492, \
                         slh_493 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_653[k] = f_12 * skh_492[k]
                   + f_6 * slg0_354[k]
                   - f_7 * slg1_354[k]
                   + f_3 * pc_x[k] * slh_492[k];

        t_654[k] = f_12 * skh_493[k]
                   + f_8 * slg0_355[k]
                   - f_9 * slg1_355[k]
                   + f_3 * pc_x[k] * slh_493[k];

        t_655[k] = f_12 * skh_342[k]
                   + f_3 * pc_z[k] * slh_489[k];
    }

#pragma omp simd aligned(t_656, t_657, t_658, pc_x, pc_y, skh_366, skh_495, skh_497, slg0_357, \
                         slg0_359, slg1_357, slg1_359, slh_492, slh_495, \
                         slh_497 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_656[k] = f_12 * skh_495[k]
                   + f_8 * slg0_357[k]
                   - f_9 * slg1_357[k]
                   + f_3 * pc_x[k] * slh_495[k];

        t_657[k] = f_14 * skh_366[k]
                   + f_3 * pc_y[k] * slh_492[k];

        t_658[k] = f_12 * skh_497[k]
                   + f_8 * slg0_359[k]
                   - f_9 * slg1_359[k]
                   + f_3 * pc_x[k] * slh_497[k];
    }

#pragma omp simd aligned(t_659, t_660, t_661, t_662, t_663, pc_x, skh_498, skh_499, skh_500, \
                         skh_501, skh_502, slh_498, slh_499, slh_500, slh_501, \
                         slh_502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_659[k] = f_12 * skh_498[k]
                   + f_3 * pc_x[k] * slh_498[k];

        t_660[k] = f_12 * skh_499[k]
                   + f_3 * pc_x[k] * slh_499[k];

        t_661[k] = f_12 * skh_500[k]
                   + f_3 * pc_x[k] * slh_500[k];

        t_662[k] = f_12 * skh_501[k]
                   + f_3 * pc_x[k] * slh_501[k];

        t_663[k] = f_12 * skh_502[k]
                   + f_3 * pc_x[k] * slh_502[k];
    }

#pragma omp simd aligned(t_664, t_665, t_666, pc_x, pc_y, pc_z, skh_351, skh_372, skh_503, \
                         slg0_355, slg1_355, slh_498, slh_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_664[k] = f_12 * skh_503[k]
                   + f_3 * pc_x[k] * slh_503[k];

        t_665[k] = f_14 * skh_372[k]
                   + f_1 * slg0_355[k]
                   - f_2 * slg1_355[k]
                   + f_3 * pc_y[k] * slh_498[k];

        t_666[k] = f_12 * skh_351[k]
                   + f_3 * pc_z[k] * slh_498[k];
    }

#pragma omp simd aligned(t_667, t_668, t_669, pc_y, skh_374, skh_375, skh_376, slg0_357, \
                         slg0_358, slg0_359, slg1_357, slg1_358, slg1_359, slh_500, slh_501, \
                         slh_502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_667[k] = f_14 * skh_374[k]
                   + f_4 * slg0_357[k]
                   - f_5 * slg1_357[k]
                   + f_3 * pc_y[k] * slh_500[k];

        t_668[k] = f_14 * skh_375[k]
                   + f_6 * slg0_358[k]
                   - f_7 * slg1_358[k]
                   + f_3 * pc_y[k] * slh_501[k];

        t_669[k] = f_14 * skh_376[k]
                   + f_8 * slg0_359[k]
                   - f_9 * slg1_359[k]
                   + f_3 * pc_y[k] * slh_502[k];
    }

#pragma omp simd aligned(t_670, t_671, t_672, pc_x, pc_y, pc_z, skh_356, skh_377, skh_504, \
                         slg0_359, slg0_360, slg1_359, slg1_360, slh_503, \
                         slh_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = f_14 * skh_377[k]
                   + f_3 * pc_y[k] * slh_503[k];

        t_671[k] = f_12 * skh_356[k]
                   + f_1 * slg0_359[k]
                   - f_2 * slg1_359[k]
                   + f_3 * pc_z[k] * slh_503[k];

        t_672[k] = f_12 * skh_504[k]
                   + f_1 * slg0_360[k]
                   - f_2 * slg1_360[k]
                   + f_3 * pc_x[k] * slh_504[k];
    }

#pragma omp simd aligned(t_673, t_674, t_675, t_676, pc_x, pc_y, pc_z, skh_357, skh_378, \
                         skh_380, skh_507, slg0_363, slg1_363, slh_504, slh_506, \
                         slh_507 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_673[k] = f_13 * skh_378[k]
                   + f_3 * pc_y[k] * slh_504[k];

        t_674[k] = f_13 * skh_357[k]
                   + f_3 * pc_z[k] * slh_504[k];

        t_675[k] = f_12 * skh_507[k]
                   + f_4 * slg0_363[k]
                   - f_5 * slg1_363[k]
                   + f_3 * pc_x[k] * slh_507[k];

        t_676[k] = f_13 * skh_380[k]
                   + f_3 * pc_y[k] * slh_506[k];
    }

#pragma omp simd aligned(t_677, t_678, t_679, pc_x, pc_z, skh_360, skh_509, skh_510, slg0_365, \
                         slg0_366, slg1_365, slg1_366, slh_507, slh_509, \
                         slh_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_677[k] = f_12 * skh_509[k]
                   + f_4 * slg0_365[k]
                   - f_5 * slg1_365[k]
                   + f_3 * pc_x[k] * slh_509[k];

        t_678[k] = f_12 * skh_510[k]
                   + f_6 * slg0_366[k]
                   - f_7 * slg1_366[k]
                   + f_3 * pc_x[k] * slh_510[k];

        t_679[k] = f_13 * skh_360[k]
                   + f_3 * pc_z[k] * slh_507[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, pc_x, pc_y, skh_383, skh_513, skh_514, slg0_369, \
                         slg0_370, slg1_369, slg1_370, slh_509, slh_513, \
                         slh_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = f_13 * skh_383[k]
                   + f_3 * pc_y[k] * slh_509[k];

        t_681[k] = f_12 * skh_513[k]
                   + f_6 * slg0_369[k]
                   - f_7 * slg1_369[k]
                   + f_3 * pc_x[k] * slh_513[k];

        t_682[k] = f_12 * skh_514[k]
                   + f_8 * slg0_370[k]
                   - f_9 * slg1_370[k]
                   + f_3 * pc_x[k] * slh_514[k];
    }

#pragma omp simd aligned(t_683, t_684, t_685, pc_x, pc_y, pc_z, skh_363, skh_387, skh_516, \
                         slg0_372, slg1_372, slh_510, slh_513, \
                         slh_516 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_683[k] = f_13 * skh_363[k]
                   + f_3 * pc_z[k] * slh_510[k];

        t_684[k] = f_12 * skh_516[k]
                   + f_8 * slg0_372[k]
                   - f_9 * slg1_372[k]
                   + f_3 * pc_x[k] * slh_516[k];

        t_685[k] = f_13 * skh_387[k]
                   + f_3 * pc_y[k] * slh_513[k];
    }
}

static auto
compute_prim_sli_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t ski0,
                                                          const size_t skh, const size_t ski1,
                                                          const size_t slg0, const size_t slg1,
                                                          const size_t slh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.5 / gamma;
    const auto f_5 = 1.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 0.5 / gamma;
    const auto f_9 = 0.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 3.5 / q;
    const auto f_16 = 3.0 / q;
    const auto f_17 = 2.5 / q;

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
    auto *t_796 = buffer.data(target + 796);
    auto *t_797 = buffer.data(target + 797);
    auto *t_798 = buffer.data(target + 798);
    auto *t_799 = buffer.data(target + 799);
    auto *t_800 = buffer.data(target + 800);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ski0_560 = buffer.data(ski0 + 560);
    const auto *ski0_563 = buffer.data(ski0 + 563);
    const auto *ski0_565 = buffer.data(ski0 + 565);
    const auto *ski0_566 = buffer.data(ski0 + 566);
    const auto *ski0_569 = buffer.data(ski0 + 569);
    const auto *ski0_570 = buffer.data(ski0 + 570);
    const auto *ski0_572 = buffer.data(ski0 + 572);
    const auto *ski0_574 = buffer.data(ski0 + 574);
    const auto *ski0_587 = buffer.data(ski0 + 587);
    const auto *ski0_784 = buffer.data(ski0 + 784);
    const auto *ski0_787 = buffer.data(ski0 + 787);
    const auto *ski0_789 = buffer.data(ski0 + 789);
    const auto *ski0_790 = buffer.data(ski0 + 790);
    const auto *ski0_793 = buffer.data(ski0 + 793);
    const auto *ski0_794 = buffer.data(ski0 + 794);
    const auto *ski0_796 = buffer.data(ski0 + 796);
    const auto *ski0_798 = buffer.data(ski0 + 798);

    const auto *skh_372 = buffer.data(skh + 372);
    const auto *skh_377 = buffer.data(skh + 377);
    const auto *skh_378 = buffer.data(skh + 378);
    const auto *skh_381 = buffer.data(skh + 381);
    const auto *skh_384 = buffer.data(skh + 384);
    const auto *skh_393 = buffer.data(skh + 393);
    const auto *skh_395 = buffer.data(skh + 395);
    const auto *skh_396 = buffer.data(skh + 396);
    const auto *skh_397 = buffer.data(skh + 397);
    const auto *skh_398 = buffer.data(skh + 398);
    const auto *skh_399 = buffer.data(skh + 399);
    const auto *skh_401 = buffer.data(skh + 401);
    const auto *skh_402 = buffer.data(skh + 402);
    const auto *skh_404 = buffer.data(skh + 404);
    const auto *skh_405 = buffer.data(skh + 405);
    const auto *skh_408 = buffer.data(skh + 408);
    const auto *skh_414 = buffer.data(skh + 414);
    const auto *skh_416 = buffer.data(skh + 416);
    const auto *skh_417 = buffer.data(skh + 417);
    const auto *skh_418 = buffer.data(skh + 418);
    const auto *skh_419 = buffer.data(skh + 419);
    const auto *skh_420 = buffer.data(skh + 420);
    const auto *skh_421 = buffer.data(skh + 421);
    const auto *skh_422 = buffer.data(skh + 422);
    const auto *skh_423 = buffer.data(skh + 423);
    const auto *skh_425 = buffer.data(skh + 425);
    const auto *skh_426 = buffer.data(skh + 426);
    const auto *skh_428 = buffer.data(skh + 428);
    const auto *skh_429 = buffer.data(skh + 429);
    const auto *skh_435 = buffer.data(skh + 435);
    const auto *skh_437 = buffer.data(skh + 437);
    const auto *skh_438 = buffer.data(skh + 438);
    const auto *skh_439 = buffer.data(skh + 439);
    const auto *skh_440 = buffer.data(skh + 440);
    const auto *skh_441 = buffer.data(skh + 441);
    const auto *skh_443 = buffer.data(skh + 443);
    const auto *skh_446 = buffer.data(skh + 446);
    const auto *skh_450 = buffer.data(skh + 450);
    const auto *skh_518 = buffer.data(skh + 518);
    const auto *skh_519 = buffer.data(skh + 519);
    const auto *skh_520 = buffer.data(skh + 520);
    const auto *skh_521 = buffer.data(skh + 521);
    const auto *skh_522 = buffer.data(skh + 522);
    const auto *skh_523 = buffer.data(skh + 523);
    const auto *skh_524 = buffer.data(skh + 524);
    const auto *skh_525 = buffer.data(skh + 525);
    const auto *skh_528 = buffer.data(skh + 528);
    const auto *skh_530 = buffer.data(skh + 530);
    const auto *skh_531 = buffer.data(skh + 531);
    const auto *skh_534 = buffer.data(skh + 534);
    const auto *skh_535 = buffer.data(skh + 535);
    const auto *skh_537 = buffer.data(skh + 537);
    const auto *skh_539 = buffer.data(skh + 539);
    const auto *skh_540 = buffer.data(skh + 540);
    const auto *skh_541 = buffer.data(skh + 541);
    const auto *skh_542 = buffer.data(skh + 542);
    const auto *skh_543 = buffer.data(skh + 543);
    const auto *skh_544 = buffer.data(skh + 544);
    const auto *skh_545 = buffer.data(skh + 545);
    const auto *skh_561 = buffer.data(skh + 561);
    const auto *skh_562 = buffer.data(skh + 562);
    const auto *skh_563 = buffer.data(skh + 563);
    const auto *skh_564 = buffer.data(skh + 564);
    const auto *skh_565 = buffer.data(skh + 565);
    const auto *skh_566 = buffer.data(skh + 566);
    const auto *skh_567 = buffer.data(skh + 567);
    const auto *skh_570 = buffer.data(skh + 570);
    const auto *skh_572 = buffer.data(skh + 572);
    const auto *skh_573 = buffer.data(skh + 573);
    const auto *skh_576 = buffer.data(skh + 576);
    const auto *skh_577 = buffer.data(skh + 577);
    const auto *skh_579 = buffer.data(skh + 579);
    const auto *skh_581 = buffer.data(skh + 581);
    const auto *skh_582 = buffer.data(skh + 582);
    const auto *skh_583 = buffer.data(skh + 583);
    const auto *skh_584 = buffer.data(skh + 584);
    const auto *skh_585 = buffer.data(skh + 585);
    const auto *skh_586 = buffer.data(skh + 586);
    const auto *skh_587 = buffer.data(skh + 587);
    const auto *skh_588 = buffer.data(skh + 588);
    const auto *skh_591 = buffer.data(skh + 591);
    const auto *skh_593 = buffer.data(skh + 593);
    const auto *skh_594 = buffer.data(skh + 594);
    const auto *skh_597 = buffer.data(skh + 597);
    const auto *skh_598 = buffer.data(skh + 598);
    const auto *skh_600 = buffer.data(skh + 600);
    const auto *skh_602 = buffer.data(skh + 602);
    const auto *skh_603 = buffer.data(skh + 603);
    const auto *skh_604 = buffer.data(skh + 604);

    const auto *ski1_560 = buffer.data(ski1 + 560);
    const auto *ski1_563 = buffer.data(ski1 + 563);
    const auto *ski1_565 = buffer.data(ski1 + 565);
    const auto *ski1_566 = buffer.data(ski1 + 566);
    const auto *ski1_569 = buffer.data(ski1 + 569);
    const auto *ski1_570 = buffer.data(ski1 + 570);
    const auto *ski1_572 = buffer.data(ski1 + 572);
    const auto *ski1_574 = buffer.data(ski1 + 574);
    const auto *ski1_587 = buffer.data(ski1 + 587);
    const auto *ski1_784 = buffer.data(ski1 + 784);
    const auto *ski1_787 = buffer.data(ski1 + 787);
    const auto *ski1_789 = buffer.data(ski1 + 789);
    const auto *ski1_790 = buffer.data(ski1 + 790);
    const auto *ski1_793 = buffer.data(ski1 + 793);
    const auto *ski1_794 = buffer.data(ski1 + 794);
    const auto *ski1_796 = buffer.data(ski1 + 796);
    const auto *ski1_798 = buffer.data(ski1 + 798);

    const auto *slg0_370 = buffer.data(slg0 + 370);
    const auto *slg0_372 = buffer.data(slg0 + 372);
    const auto *slg0_373 = buffer.data(slg0 + 373);
    const auto *slg0_374 = buffer.data(slg0 + 374);
    const auto *slg0_375 = buffer.data(slg0 + 375);
    const auto *slg0_378 = buffer.data(slg0 + 378);
    const auto *slg0_380 = buffer.data(slg0 + 380);
    const auto *slg0_381 = buffer.data(slg0 + 381);
    const auto *slg0_384 = buffer.data(slg0 + 384);
    const auto *slg0_385 = buffer.data(slg0 + 385);
    const auto *slg0_387 = buffer.data(slg0 + 387);
    const auto *slg0_388 = buffer.data(slg0 + 388);
    const auto *slg0_389 = buffer.data(slg0 + 389);
    const auto *slg0_400 = buffer.data(slg0 + 400);
    const auto *slg0_402 = buffer.data(slg0 + 402);
    const auto *slg0_403 = buffer.data(slg0 + 403);
    const auto *slg0_404 = buffer.data(slg0 + 404);
    const auto *slg0_405 = buffer.data(slg0 + 405);
    const auto *slg0_408 = buffer.data(slg0 + 408);
    const auto *slg0_410 = buffer.data(slg0 + 410);
    const auto *slg0_411 = buffer.data(slg0 + 411);
    const auto *slg0_414 = buffer.data(slg0 + 414);
    const auto *slg0_415 = buffer.data(slg0 + 415);
    const auto *slg0_417 = buffer.data(slg0 + 417);
    const auto *slg0_418 = buffer.data(slg0 + 418);
    const auto *slg0_419 = buffer.data(slg0 + 419);

    const auto *slg1_370 = buffer.data(slg1 + 370);
    const auto *slg1_372 = buffer.data(slg1 + 372);
    const auto *slg1_373 = buffer.data(slg1 + 373);
    const auto *slg1_374 = buffer.data(slg1 + 374);
    const auto *slg1_375 = buffer.data(slg1 + 375);
    const auto *slg1_378 = buffer.data(slg1 + 378);
    const auto *slg1_380 = buffer.data(slg1 + 380);
    const auto *slg1_381 = buffer.data(slg1 + 381);
    const auto *slg1_384 = buffer.data(slg1 + 384);
    const auto *slg1_385 = buffer.data(slg1 + 385);
    const auto *slg1_387 = buffer.data(slg1 + 387);
    const auto *slg1_388 = buffer.data(slg1 + 388);
    const auto *slg1_389 = buffer.data(slg1 + 389);
    const auto *slg1_400 = buffer.data(slg1 + 400);
    const auto *slg1_402 = buffer.data(slg1 + 402);
    const auto *slg1_403 = buffer.data(slg1 + 403);
    const auto *slg1_404 = buffer.data(slg1 + 404);
    const auto *slg1_405 = buffer.data(slg1 + 405);
    const auto *slg1_408 = buffer.data(slg1 + 408);
    const auto *slg1_410 = buffer.data(slg1 + 410);
    const auto *slg1_411 = buffer.data(slg1 + 411);
    const auto *slg1_414 = buffer.data(slg1 + 414);
    const auto *slg1_415 = buffer.data(slg1 + 415);
    const auto *slg1_417 = buffer.data(slg1 + 417);
    const auto *slg1_418 = buffer.data(slg1 + 418);
    const auto *slg1_419 = buffer.data(slg1 + 419);

    const auto *slh_518 = buffer.data(slh + 518);
    const auto *slh_519 = buffer.data(slh + 519);
    const auto *slh_520 = buffer.data(slh + 520);
    const auto *slh_521 = buffer.data(slh + 521);
    const auto *slh_522 = buffer.data(slh + 522);
    const auto *slh_523 = buffer.data(slh + 523);
    const auto *slh_524 = buffer.data(slh + 524);
    const auto *slh_525 = buffer.data(slh + 525);
    const auto *slh_527 = buffer.data(slh + 527);
    const auto *slh_528 = buffer.data(slh + 528);
    const auto *slh_530 = buffer.data(slh + 530);
    const auto *slh_531 = buffer.data(slh + 531);
    const auto *slh_534 = buffer.data(slh + 534);
    const auto *slh_535 = buffer.data(slh + 535);
    const auto *slh_537 = buffer.data(slh + 537);
    const auto *slh_539 = buffer.data(slh + 539);
    const auto *slh_540 = buffer.data(slh + 540);
    const auto *slh_541 = buffer.data(slh + 541);
    const auto *slh_542 = buffer.data(slh + 542);
    const auto *slh_543 = buffer.data(slh + 543);
    const auto *slh_544 = buffer.data(slh + 544);
    const auto *slh_545 = buffer.data(slh + 545);
    const auto *slh_546 = buffer.data(slh + 546);
    const auto *slh_548 = buffer.data(slh + 548);
    const auto *slh_549 = buffer.data(slh + 549);
    const auto *slh_551 = buffer.data(slh + 551);
    const auto *slh_552 = buffer.data(slh + 552);
    const auto *slh_555 = buffer.data(slh + 555);
    const auto *slh_561 = buffer.data(slh + 561);
    const auto *slh_562 = buffer.data(slh + 562);
    const auto *slh_563 = buffer.data(slh + 563);
    const auto *slh_564 = buffer.data(slh + 564);
    const auto *slh_565 = buffer.data(slh + 565);
    const auto *slh_566 = buffer.data(slh + 566);
    const auto *slh_567 = buffer.data(slh + 567);
    const auto *slh_569 = buffer.data(slh + 569);
    const auto *slh_570 = buffer.data(slh + 570);
    const auto *slh_572 = buffer.data(slh + 572);
    const auto *slh_573 = buffer.data(slh + 573);
    const auto *slh_576 = buffer.data(slh + 576);
    const auto *slh_577 = buffer.data(slh + 577);
    const auto *slh_579 = buffer.data(slh + 579);
    const auto *slh_581 = buffer.data(slh + 581);
    const auto *slh_582 = buffer.data(slh + 582);
    const auto *slh_583 = buffer.data(slh + 583);
    const auto *slh_584 = buffer.data(slh + 584);
    const auto *slh_585 = buffer.data(slh + 585);
    const auto *slh_586 = buffer.data(slh + 586);
    const auto *slh_587 = buffer.data(slh + 587);
    const auto *slh_588 = buffer.data(slh + 588);
    const auto *slh_590 = buffer.data(slh + 590);
    const auto *slh_591 = buffer.data(slh + 591);
    const auto *slh_593 = buffer.data(slh + 593);
    const auto *slh_594 = buffer.data(slh + 594);
    const auto *slh_597 = buffer.data(slh + 597);
    const auto *slh_603 = buffer.data(slh + 603);
    const auto *slh_604 = buffer.data(slh + 604);

#pragma omp simd aligned(t_686, t_687, t_688, t_689, pc_x, skh_518, skh_519, skh_520, skh_521, \
                         slg0_374, slg1_374, slh_518, slh_519, slh_520, \
                         slh_521 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_686[k] = f_12 * skh_518[k]
                   + f_8 * slg0_374[k]
                   - f_9 * slg1_374[k]
                   + f_3 * pc_x[k] * slh_518[k];

        t_687[k] = f_12 * skh_519[k]
                   + f_3 * pc_x[k] * slh_519[k];

        t_688[k] = f_12 * skh_520[k]
                   + f_3 * pc_x[k] * slh_520[k];

        t_689[k] = f_12 * skh_521[k]
                   + f_3 * pc_x[k] * slh_521[k];
    }

#pragma omp simd aligned(t_690, t_691, t_692, t_693, pc_x, pc_y, skh_393, skh_522, skh_523, \
                         skh_524, slg0_370, slg1_370, slh_519, slh_522, slh_523, \
                         slh_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_690[k] = f_12 * skh_522[k]
                   + f_3 * pc_x[k] * slh_522[k];

        t_691[k] = f_12 * skh_523[k]
                   + f_3 * pc_x[k] * slh_523[k];

        t_692[k] = f_12 * skh_524[k]
                   + f_3 * pc_x[k] * slh_524[k];

        t_693[k] = f_13 * skh_393[k]
                   + f_1 * slg0_370[k]
                   - f_2 * slg1_370[k]
                   + f_3 * pc_y[k] * slh_519[k];
    }

#pragma omp simd aligned(t_694, t_695, t_696, pc_y, pc_z, skh_372, skh_395, skh_396, slg0_372, \
                         slg0_373, slg1_372, slg1_373, slh_519, slh_521, \
                         slh_522 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_694[k] = f_13 * skh_372[k]
                   + f_3 * pc_z[k] * slh_519[k];

        t_695[k] = f_13 * skh_395[k]
                   + f_4 * slg0_372[k]
                   - f_5 * slg1_372[k]
                   + f_3 * pc_y[k] * slh_521[k];

        t_696[k] = f_13 * skh_396[k]
                   + f_6 * slg0_373[k]
                   - f_7 * slg1_373[k]
                   + f_3 * pc_y[k] * slh_522[k];
    }

#pragma omp simd aligned(t_697, t_698, t_699, pc_y, pc_z, skh_377, skh_397, skh_398, slg0_374, \
                         slg1_374, slh_523, slh_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_697[k] = f_13 * skh_397[k]
                   + f_8 * slg0_374[k]
                   - f_9 * slg1_374[k]
                   + f_3 * pc_y[k] * slh_523[k];

        t_698[k] = f_13 * skh_398[k]
                   + f_3 * pc_y[k] * slh_524[k];

        t_699[k] = f_13 * skh_377[k]
                   + f_1 * slg0_374[k]
                   - f_2 * slg1_374[k]
                   + f_3 * pc_z[k] * slh_524[k];
    }

#pragma omp simd aligned(t_700, t_701, t_702, pc_x, pc_y, pc_z, skh_378, skh_399, skh_525, \
                         slg0_375, slg1_375, slh_525 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_700[k] = f_12 * skh_525[k]
                   + f_1 * slg0_375[k]
                   - f_2 * slg1_375[k]
                   + f_3 * pc_x[k] * slh_525[k];

        t_701[k] = f_12 * skh_399[k]
                   + f_3 * pc_y[k] * slh_525[k];

        t_702[k] = f_14 * skh_378[k]
                   + f_3 * pc_z[k] * slh_525[k];
    }

#pragma omp simd aligned(t_703, t_704, t_705, pc_x, pc_y, skh_401, skh_528, skh_530, slg0_378, \
                         slg0_380, slg1_378, slg1_380, slh_527, slh_528, \
                         slh_530 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_703[k] = f_12 * skh_528[k]
                   + f_4 * slg0_378[k]
                   - f_5 * slg1_378[k]
                   + f_3 * pc_x[k] * slh_528[k];

        t_704[k] = f_12 * skh_401[k]
                   + f_3 * pc_y[k] * slh_527[k];

        t_705[k] = f_12 * skh_530[k]
                   + f_4 * slg0_380[k]
                   - f_5 * slg1_380[k]
                   + f_3 * pc_x[k] * slh_530[k];
    }

#pragma omp simd aligned(t_706, t_707, t_708, pc_x, pc_y, pc_z, skh_381, skh_404, skh_531, \
                         slg0_381, slg1_381, slh_528, slh_530, \
                         slh_531 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_706[k] = f_12 * skh_531[k]
                   + f_6 * slg0_381[k]
                   - f_7 * slg1_381[k]
                   + f_3 * pc_x[k] * slh_531[k];

        t_707[k] = f_14 * skh_381[k]
                   + f_3 * pc_z[k] * slh_528[k];

        t_708[k] = f_12 * skh_404[k]
                   + f_3 * pc_y[k] * slh_530[k];
    }

#pragma omp simd aligned(t_709, t_710, t_711, pc_x, pc_z, skh_384, skh_534, skh_535, slg0_384, \
                         slg0_385, slg1_384, slg1_385, slh_531, slh_534, \
                         slh_535 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_709[k] = f_12 * skh_534[k]
                   + f_6 * slg0_384[k]
                   - f_7 * slg1_384[k]
                   + f_3 * pc_x[k] * slh_534[k];

        t_710[k] = f_12 * skh_535[k]
                   + f_8 * slg0_385[k]
                   - f_9 * slg1_385[k]
                   + f_3 * pc_x[k] * slh_535[k];

        t_711[k] = f_14 * skh_384[k]
                   + f_3 * pc_z[k] * slh_531[k];
    }

#pragma omp simd aligned(t_712, t_713, t_714, pc_x, pc_y, skh_408, skh_537, skh_539, slg0_387, \
                         slg0_389, slg1_387, slg1_389, slh_534, slh_537, \
                         slh_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_712[k] = f_12 * skh_537[k]
                   + f_8 * slg0_387[k]
                   - f_9 * slg1_387[k]
                   + f_3 * pc_x[k] * slh_537[k];

        t_713[k] = f_12 * skh_408[k]
                   + f_3 * pc_y[k] * slh_534[k];

        t_714[k] = f_12 * skh_539[k]
                   + f_8 * slg0_389[k]
                   - f_9 * slg1_389[k]
                   + f_3 * pc_x[k] * slh_539[k];
    }

#pragma omp simd aligned(t_715, t_716, t_717, t_718, t_719, pc_x, skh_540, skh_541, skh_542, \
                         skh_543, skh_544, slh_540, slh_541, slh_542, slh_543, \
                         slh_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_715[k] = f_12 * skh_540[k]
                   + f_3 * pc_x[k] * slh_540[k];

        t_716[k] = f_12 * skh_541[k]
                   + f_3 * pc_x[k] * slh_541[k];

        t_717[k] = f_12 * skh_542[k]
                   + f_3 * pc_x[k] * slh_542[k];

        t_718[k] = f_12 * skh_543[k]
                   + f_3 * pc_x[k] * slh_543[k];

        t_719[k] = f_12 * skh_544[k]
                   + f_3 * pc_x[k] * slh_544[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, pc_x, pc_y, pc_z, skh_393, skh_414, skh_545, \
                         slg0_385, slg1_385, slh_540, slh_545 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = f_12 * skh_545[k]
                   + f_3 * pc_x[k] * slh_545[k];

        t_721[k] = f_12 * skh_414[k]
                   + f_1 * slg0_385[k]
                   - f_2 * slg1_385[k]
                   + f_3 * pc_y[k] * slh_540[k];

        t_722[k] = f_14 * skh_393[k]
                   + f_3 * pc_z[k] * slh_540[k];
    }

#pragma omp simd aligned(t_723, t_724, t_725, pc_y, skh_416, skh_417, skh_418, slg0_387, \
                         slg0_388, slg0_389, slg1_387, slg1_388, slg1_389, slh_542, slh_543, \
                         slh_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_723[k] = f_12 * skh_416[k]
                   + f_4 * slg0_387[k]
                   - f_5 * slg1_387[k]
                   + f_3 * pc_y[k] * slh_542[k];

        t_724[k] = f_12 * skh_417[k]
                   + f_6 * slg0_388[k]
                   - f_7 * slg1_388[k]
                   + f_3 * pc_y[k] * slh_543[k];

        t_725[k] = f_12 * skh_418[k]
                   + f_8 * slg0_389[k]
                   - f_9 * slg1_389[k]
                   + f_3 * pc_y[k] * slh_544[k];
    }

#pragma omp simd aligned(t_726, t_727, t_728, t_729, pb_y, pc_y, pc_z, ski0_560, skh_398, \
                         skh_419, skh_420, ski1_560, slg0_389, slg1_389, slh_545, \
                         slh_546 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_726[k] = f_12 * skh_419[k]
                   + f_3 * pc_y[k] * slh_545[k];

        t_727[k] = f_14 * skh_398[k]
                   + f_1 * slg0_389[k]
                   - f_2 * slg1_389[k]
                   + f_3 * pc_z[k] * slh_545[k];

        t_728[k] = pb_y[k] * ski0_560[k]
                   - f_10 * pc_y[k] * ski1_560[k];

        t_729[k] = f_11 * skh_420[k]
                   + f_3 * pc_y[k] * slh_546[k];
    }

#pragma omp simd aligned(t_730, t_731, t_732, t_733, pb_y, pc_y, pc_z, ski0_563, ski0_565, \
                         skh_399, skh_421, skh_422, ski1_563, ski1_565, slh_546, \
                         slh_548 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_730[k] = f_17 * skh_399[k]
                   + f_3 * pc_z[k] * slh_546[k];

        t_731[k] = pb_y[k] * ski0_563[k]
                   + f_12 * skh_421[k]
                   - f_10 * pc_y[k] * ski1_563[k];

        t_732[k] = f_11 * skh_422[k]
                   + f_3 * pc_y[k] * slh_548[k];

        t_733[k] = pb_y[k] * ski0_565[k]
                   - f_10 * pc_y[k] * ski1_565[k];
    }

#pragma omp simd aligned(t_734, t_735, t_736, t_737, pb_y, pc_y, pc_z, ski0_566, ski0_569, \
                         skh_402, skh_423, skh_425, ski1_566, ski1_569, slh_549, \
                         slh_551 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_734[k] = pb_y[k] * ski0_566[k]
                   + f_13 * skh_423[k]
                   - f_10 * pc_y[k] * ski1_566[k];

        t_735[k] = f_17 * skh_402[k]
                   + f_3 * pc_z[k] * slh_549[k];

        t_736[k] = f_11 * skh_425[k]
                   + f_3 * pc_y[k] * slh_551[k];

        t_737[k] = pb_y[k] * ski0_569[k]
                   - f_10 * pc_y[k] * ski1_569[k];
    }

#pragma omp simd aligned(t_738, t_739, t_740, pb_y, pc_y, pc_z, ski0_570, ski0_572, skh_405, \
                         skh_426, skh_428, ski1_570, ski1_572, \
                         slh_552 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_738[k] = pb_y[k] * ski0_570[k]
                   + f_14 * skh_426[k]
                   - f_10 * pc_y[k] * ski1_570[k];

        t_739[k] = f_17 * skh_405[k]
                   + f_3 * pc_z[k] * slh_552[k];

        t_740[k] = pb_y[k] * ski0_572[k]
                   + f_12 * skh_428[k]
                   - f_10 * pc_y[k] * ski1_572[k];
    }

#pragma omp simd aligned(t_741, t_742, t_743, t_744, pb_y, pc_x, pc_y, ski0_574, skh_429, \
                         skh_561, skh_562, ski1_574, slh_555, slh_561, \
                         slh_562 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_741[k] = f_11 * skh_429[k]
                   + f_3 * pc_y[k] * slh_555[k];

        t_742[k] = pb_y[k] * ski0_574[k]
                   - f_10 * pc_y[k] * ski1_574[k];

        t_743[k] = f_12 * skh_561[k]
                   + f_3 * pc_x[k] * slh_561[k];

        t_744[k] = f_12 * skh_562[k]
                   + f_3 * pc_x[k] * slh_562[k];
    }

#pragma omp simd aligned(t_745, t_746, t_747, t_748, pc_x, skh_563, skh_564, skh_565, skh_566, \
                         slh_563, slh_564, slh_565, slh_566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_745[k] = f_12 * skh_563[k]
                   + f_3 * pc_x[k] * slh_563[k];

        t_746[k] = f_12 * skh_564[k]
                   + f_3 * pc_x[k] * slh_564[k];

        t_747[k] = f_12 * skh_565[k]
                   + f_3 * pc_x[k] * slh_565[k];

        t_748[k] = f_12 * skh_566[k]
                   + f_3 * pc_x[k] * slh_566[k];
    }

#pragma omp simd aligned(t_749, t_750, t_751, pc_y, pc_z, skh_414, skh_435, skh_437, slg0_400, \
                         slg0_402, slg1_400, slg1_402, slh_561, \
                         slh_563 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_749[k] = f_11 * skh_435[k]
                   + f_1 * slg0_400[k]
                   - f_2 * slg1_400[k]
                   + f_3 * pc_y[k] * slh_561[k];

        t_750[k] = f_17 * skh_414[k]
                   + f_3 * pc_z[k] * slh_561[k];

        t_751[k] = f_11 * skh_437[k]
                   + f_4 * slg0_402[k]
                   - f_5 * slg1_402[k]
                   + f_3 * pc_y[k] * slh_563[k];
    }

#pragma omp simd aligned(t_752, t_753, t_754, pc_y, skh_438, skh_439, skh_440, slg0_403, \
                         slg0_404, slg1_403, slg1_404, slh_564, slh_565, \
                         slh_566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_752[k] = f_11 * skh_438[k]
                   + f_6 * slg0_403[k]
                   - f_7 * slg1_403[k]
                   + f_3 * pc_y[k] * slh_564[k];

        t_753[k] = f_11 * skh_439[k]
                   + f_8 * slg0_404[k]
                   - f_9 * slg1_404[k]
                   + f_3 * pc_y[k] * slh_565[k];

        t_754[k] = f_11 * skh_440[k]
                   + f_3 * pc_y[k] * slh_566[k];
    }

#pragma omp simd aligned(t_755, t_756, t_757, t_758, pb_y, pc_x, pc_y, pc_z, ski0_587, \
                         skh_420, skh_567, ski1_587, slg0_405, slg1_405, \
                         slh_567 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_755[k] = pb_y[k] * ski0_587[k]
                   - f_10 * pc_y[k] * ski1_587[k];

        t_756[k] = f_12 * skh_567[k]
                   + f_1 * slg0_405[k]
                   - f_2 * slg1_405[k]
                   + f_3 * pc_x[k] * slh_567[k];

        t_757[k] = f_3 * pc_y[k] * slh_567[k];

        t_758[k] = f_16 * skh_420[k]
                   + f_3 * pc_z[k] * slh_567[k];
    }

#pragma omp simd aligned(t_759, t_760, t_761, pc_x, pc_y, skh_570, skh_572, slg0_408, \
                         slg0_410, slg1_408, slg1_410, slh_569, slh_570, \
                         slh_572 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_759[k] = f_12 * skh_570[k]
                   + f_4 * slg0_408[k]
                   - f_5 * slg1_408[k]
                   + f_3 * pc_x[k] * slh_570[k];

        t_760[k] = f_3 * pc_y[k] * slh_569[k];

        t_761[k] = f_12 * skh_572[k]
                   + f_4 * slg0_410[k]
                   - f_5 * slg1_410[k]
                   + f_3 * pc_x[k] * slh_572[k];
    }

#pragma omp simd aligned(t_762, t_763, t_764, pc_x, pc_y, pc_z, skh_423, skh_573, slg0_411, \
                         slg1_411, slh_570, slh_572, slh_573 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_762[k] = f_12 * skh_573[k]
                   + f_6 * slg0_411[k]
                   - f_7 * slg1_411[k]
                   + f_3 * pc_x[k] * slh_573[k];

        t_763[k] = f_16 * skh_423[k]
                   + f_3 * pc_z[k] * slh_570[k];

        t_764[k] = f_3 * pc_y[k] * slh_572[k];
    }

#pragma omp simd aligned(t_765, t_766, t_767, pc_x, pc_z, skh_426, skh_576, skh_577, slg0_414, \
                         slg0_415, slg1_414, slg1_415, slh_573, slh_576, \
                         slh_577 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_765[k] = f_12 * skh_576[k]
                   + f_6 * slg0_414[k]
                   - f_7 * slg1_414[k]
                   + f_3 * pc_x[k] * slh_576[k];

        t_766[k] = f_12 * skh_577[k]
                   + f_8 * slg0_415[k]
                   - f_9 * slg1_415[k]
                   + f_3 * pc_x[k] * slh_577[k];

        t_767[k] = f_16 * skh_426[k]
                   + f_3 * pc_z[k] * slh_573[k];
    }

#pragma omp simd aligned(t_768, t_769, t_770, pc_x, pc_y, skh_579, skh_581, slg0_417, \
                         slg0_419, slg1_417, slg1_419, slh_576, slh_579, \
                         slh_581 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_768[k] = f_12 * skh_579[k]
                   + f_8 * slg0_417[k]
                   - f_9 * slg1_417[k]
                   + f_3 * pc_x[k] * slh_579[k];

        t_769[k] = f_3 * pc_y[k] * slh_576[k];

        t_770[k] = f_12 * skh_581[k]
                   + f_8 * slg0_419[k]
                   - f_9 * slg1_419[k]
                   + f_3 * pc_x[k] * slh_581[k];
    }

#pragma omp simd aligned(t_771, t_772, t_773, t_774, t_775, pc_x, skh_582, skh_583, skh_584, \
                         skh_585, skh_586, slh_582, slh_583, slh_584, slh_585, \
                         slh_586 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_771[k] = f_12 * skh_582[k]
                   + f_3 * pc_x[k] * slh_582[k];

        t_772[k] = f_12 * skh_583[k]
                   + f_3 * pc_x[k] * slh_583[k];

        t_773[k] = f_12 * skh_584[k]
                   + f_3 * pc_x[k] * slh_584[k];

        t_774[k] = f_12 * skh_585[k]
                   + f_3 * pc_x[k] * slh_585[k];

        t_775[k] = f_12 * skh_586[k]
                   + f_3 * pc_x[k] * slh_586[k];
    }

#pragma omp simd aligned(t_776, t_777, t_778, t_779, pc_x, pc_y, pc_z, skh_435, skh_587, \
                         slg0_415, slg0_417, slg1_415, slg1_417, slh_582, slh_584, \
                         slh_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_776[k] = f_12 * skh_587[k]
                   + f_3 * pc_x[k] * slh_587[k];

        t_777[k] = f_1 * slg0_415[k]
                   - f_2 * slg1_415[k]
                   + f_3 * pc_y[k] * slh_582[k];

        t_778[k] = f_16 * skh_435[k]
                   + f_3 * pc_z[k] * slh_582[k];

        t_779[k] = f_4 * slg0_417[k]
                   - f_5 * slg1_417[k]
                   + f_3 * pc_y[k] * slh_584[k];
    }

#pragma omp simd aligned(t_780, t_781, t_782, t_783, pc_y, pc_z, skh_440, slg0_418, slg0_419, \
                         slg1_418, slg1_419, slh_585, slh_586, \
                         slh_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_780[k] = f_6 * slg0_418[k]
                   - f_7 * slg1_418[k]
                   + f_3 * pc_y[k] * slh_585[k];

        t_781[k] = f_8 * slg0_419[k]
                   - f_9 * slg1_419[k]
                   + f_3 * pc_y[k] * slh_586[k];

        t_782[k] = f_3 * pc_y[k] * slh_587[k];

        t_783[k] = f_16 * skh_440[k]
                   + f_1 * slg0_419[k]
                   - f_2 * slg1_419[k]
                   + f_3 * pc_z[k] * slh_587[k];
    }

#pragma omp simd aligned(t_784, t_785, t_786, t_787, pb_x, pc_x, pc_y, pc_z, ski0_784, \
                         ski0_787, skh_441, skh_588, skh_591, ski1_784, ski1_787, \
                         slh_588 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_784[k] = pb_x[k] * ski0_784[k]
                   + f_16 * skh_588[k]
                   - f_10 * pc_x[k] * ski1_784[k];

        t_785[k] = f_15 * skh_441[k]
                   + f_3 * pc_y[k] * slh_588[k];

        t_786[k] = f_3 * pc_z[k] * slh_588[k];

        t_787[k] = pb_x[k] * ski0_787[k]
                   + f_14 * skh_591[k]
                   - f_10 * pc_x[k] * ski1_787[k];
    }

#pragma omp simd aligned(t_788, t_789, t_790, pb_x, pc_x, pc_y, ski0_789, ski0_790, skh_443, \
                         skh_593, skh_594, ski1_789, ski1_790, \
                         slh_590 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_788[k] = f_15 * skh_443[k]
                   + f_3 * pc_y[k] * slh_590[k];

        t_789[k] = pb_x[k] * ski0_789[k]
                   + f_14 * skh_593[k]
                   - f_10 * pc_x[k] * ski1_789[k];

        t_790[k] = pb_x[k] * ski0_790[k]
                   + f_13 * skh_594[k]
                   - f_10 * pc_x[k] * ski1_790[k];
    }

#pragma omp simd aligned(t_791, t_792, t_793, pb_x, pc_x, pc_y, pc_z, ski0_793, skh_446, \
                         skh_597, ski1_793, slh_591, slh_593 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_791[k] = f_3 * pc_z[k] * slh_591[k];

        t_792[k] = f_15 * skh_446[k]
                   + f_3 * pc_y[k] * slh_593[k];

        t_793[k] = pb_x[k] * ski0_793[k]
                   + f_13 * skh_597[k]
                   - f_10 * pc_x[k] * ski1_793[k];
    }

#pragma omp simd aligned(t_794, t_795, t_796, pb_x, pc_x, pc_z, ski0_794, ski0_796, skh_598, \
                         skh_600, ski1_794, ski1_796, slh_594 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_794[k] = pb_x[k] * ski0_794[k]
                   + f_12 * skh_598[k]
                   - f_10 * pc_x[k] * ski1_794[k];

        t_795[k] = f_3 * pc_z[k] * slh_594[k];

        t_796[k] = pb_x[k] * ski0_796[k]
                   + f_12 * skh_600[k]
                   - f_10 * pc_x[k] * ski1_796[k];
    }

#pragma omp simd aligned(t_797, t_798, t_799, t_800, pb_x, pc_x, pc_y, ski0_798, skh_450, \
                         skh_602, skh_603, skh_604, ski1_798, slh_597, slh_603, \
                         slh_604 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_797[k] = f_15 * skh_450[k]
                   + f_3 * pc_y[k] * slh_597[k];

        t_798[k] = pb_x[k] * ski0_798[k]
                   + f_12 * skh_602[k]
                   - f_10 * pc_x[k] * ski1_798[k];

        t_799[k] = f_11 * skh_603[k]
                   + f_3 * pc_x[k] * slh_603[k];

        t_800[k] = f_11 * skh_604[k]
                   + f_3 * pc_x[k] * slh_604[k];
    }
}

static auto
compute_prim_sli_three_center_electron_repulsion_0_piece7(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t ski0,
                                                          const size_t skh, const size_t ski1,
                                                          const size_t slh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_3 = p / q;
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 3.5 / q;
    const auto f_16 = 3.0 / q;
    const auto f_17 = 2.5 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ski0_588 = buffer.data(ski0 + 588);
    const auto *ski0_591 = buffer.data(ski0 + 591);
    const auto *ski0_594 = buffer.data(ski0 + 594);
    const auto *ski0_598 = buffer.data(ski0 + 598);
    const auto *ski0_805 = buffer.data(ski0 + 805);
    const auto *ski0_807 = buffer.data(ski0 + 807);
    const auto *ski0_808 = buffer.data(ski0 + 808);
    const auto *ski0_809 = buffer.data(ski0 + 809);
    const auto *ski0_811 = buffer.data(ski0 + 811);
    const auto *ski0_817 = buffer.data(ski0 + 817);
    const auto *ski0_821 = buffer.data(ski0 + 821);
    const auto *ski0_824 = buffer.data(ski0 + 824);
    const auto *ski0_826 = buffer.data(ski0 + 826);
    const auto *ski0_833 = buffer.data(ski0 + 833);
    const auto *ski0_835 = buffer.data(ski0 + 835);
    const auto *ski0_836 = buffer.data(ski0 + 836);
    const auto *ski0_837 = buffer.data(ski0 + 837);
    const auto *ski0_839 = buffer.data(ski0 + 839);
    const auto *ski0_840 = buffer.data(ski0 + 840);
    const auto *ski0_843 = buffer.data(ski0 + 843);
    const auto *ski0_845 = buffer.data(ski0 + 845);
    const auto *ski0_846 = buffer.data(ski0 + 846);
    const auto *ski0_849 = buffer.data(ski0 + 849);
    const auto *ski0_850 = buffer.data(ski0 + 850);
    const auto *ski0_852 = buffer.data(ski0 + 852);
    const auto *ski0_854 = buffer.data(ski0 + 854);
    const auto *ski0_861 = buffer.data(ski0 + 861);
    const auto *ski0_863 = buffer.data(ski0 + 863);
    const auto *ski0_864 = buffer.data(ski0 + 864);
    const auto *ski0_865 = buffer.data(ski0 + 865);
    const auto *ski0_867 = buffer.data(ski0 + 867);
    const auto *ski0_868 = buffer.data(ski0 + 868);
    const auto *ski0_871 = buffer.data(ski0 + 871);
    const auto *ski0_873 = buffer.data(ski0 + 873);
    const auto *ski0_874 = buffer.data(ski0 + 874);
    const auto *ski0_877 = buffer.data(ski0 + 877);
    const auto *ski0_878 = buffer.data(ski0 + 878);
    const auto *ski0_880 = buffer.data(ski0 + 880);
    const auto *ski0_882 = buffer.data(ski0 + 882);
    const auto *ski0_889 = buffer.data(ski0 + 889);
    const auto *ski0_891 = buffer.data(ski0 + 891);
    const auto *ski0_892 = buffer.data(ski0 + 892);
    const auto *ski0_893 = buffer.data(ski0 + 893);
    const auto *ski0_895 = buffer.data(ski0 + 895);
    const auto *ski0_896 = buffer.data(ski0 + 896);
    const auto *ski0_899 = buffer.data(ski0 + 899);
    const auto *ski0_901 = buffer.data(ski0 + 901);
    const auto *ski0_902 = buffer.data(ski0 + 902);
    const auto *ski0_905 = buffer.data(ski0 + 905);
    const auto *ski0_906 = buffer.data(ski0 + 906);
    const auto *ski0_908 = buffer.data(ski0 + 908);
    const auto *ski0_910 = buffer.data(ski0 + 910);
    const auto *ski0_917 = buffer.data(ski0 + 917);
    const auto *ski0_919 = buffer.data(ski0 + 919);
    const auto *ski0_920 = buffer.data(ski0 + 920);
    const auto *ski0_921 = buffer.data(ski0 + 921);

    const auto *skh_441 = buffer.data(skh + 441);
    const auto *skh_444 = buffer.data(skh + 444);
    const auto *skh_447 = buffer.data(skh + 447);
    const auto *skh_456 = buffer.data(skh + 456);
    const auto *skh_461 = buffer.data(skh + 461);
    const auto *skh_462 = buffer.data(skh + 462);
    const auto *skh_464 = buffer.data(skh + 464);
    const auto *skh_465 = buffer.data(skh + 465);
    const auto *skh_467 = buffer.data(skh + 467);
    const auto *skh_468 = buffer.data(skh + 468);
    const auto *skh_471 = buffer.data(skh + 471);
    const auto *skh_477 = buffer.data(skh + 477);
    const auto *skh_482 = buffer.data(skh + 482);
    const auto *skh_483 = buffer.data(skh + 483);
    const auto *skh_485 = buffer.data(skh + 485);
    const auto *skh_486 = buffer.data(skh + 486);
    const auto *skh_488 = buffer.data(skh + 488);
    const auto *skh_489 = buffer.data(skh + 489);
    const auto *skh_492 = buffer.data(skh + 492);
    const auto *skh_498 = buffer.data(skh + 498);
    const auto *skh_503 = buffer.data(skh + 503);
    const auto *skh_504 = buffer.data(skh + 504);
    const auto *skh_506 = buffer.data(skh + 506);
    const auto *skh_507 = buffer.data(skh + 507);
    const auto *skh_509 = buffer.data(skh + 509);
    const auto *skh_510 = buffer.data(skh + 510);
    const auto *skh_513 = buffer.data(skh + 513);
    const auto *skh_519 = buffer.data(skh + 519);
    const auto *skh_524 = buffer.data(skh + 524);
    const auto *skh_525 = buffer.data(skh + 525);
    const auto *skh_527 = buffer.data(skh + 527);
    const auto *skh_530 = buffer.data(skh + 530);
    const auto *skh_534 = buffer.data(skh + 534);
    const auto *skh_605 = buffer.data(skh + 605);
    const auto *skh_606 = buffer.data(skh + 606);
    const auto *skh_607 = buffer.data(skh + 607);
    const auto *skh_608 = buffer.data(skh + 608);
    const auto *skh_614 = buffer.data(skh + 614);
    const auto *skh_618 = buffer.data(skh + 618);
    const auto *skh_621 = buffer.data(skh + 621);
    const auto *skh_623 = buffer.data(skh + 623);
    const auto *skh_624 = buffer.data(skh + 624);
    const auto *skh_625 = buffer.data(skh + 625);
    const auto *skh_626 = buffer.data(skh + 626);
    const auto *skh_627 = buffer.data(skh + 627);
    const auto *skh_628 = buffer.data(skh + 628);
    const auto *skh_629 = buffer.data(skh + 629);
    const auto *skh_630 = buffer.data(skh + 630);
    const auto *skh_633 = buffer.data(skh + 633);
    const auto *skh_635 = buffer.data(skh + 635);
    const auto *skh_636 = buffer.data(skh + 636);
    const auto *skh_639 = buffer.data(skh + 639);
    const auto *skh_640 = buffer.data(skh + 640);
    const auto *skh_642 = buffer.data(skh + 642);
    const auto *skh_644 = buffer.data(skh + 644);
    const auto *skh_645 = buffer.data(skh + 645);
    const auto *skh_646 = buffer.data(skh + 646);
    const auto *skh_647 = buffer.data(skh + 647);
    const auto *skh_648 = buffer.data(skh + 648);
    const auto *skh_649 = buffer.data(skh + 649);
    const auto *skh_650 = buffer.data(skh + 650);
    const auto *skh_651 = buffer.data(skh + 651);
    const auto *skh_654 = buffer.data(skh + 654);
    const auto *skh_656 = buffer.data(skh + 656);
    const auto *skh_657 = buffer.data(skh + 657);
    const auto *skh_660 = buffer.data(skh + 660);
    const auto *skh_661 = buffer.data(skh + 661);
    const auto *skh_663 = buffer.data(skh + 663);
    const auto *skh_665 = buffer.data(skh + 665);
    const auto *skh_666 = buffer.data(skh + 666);
    const auto *skh_667 = buffer.data(skh + 667);
    const auto *skh_668 = buffer.data(skh + 668);
    const auto *skh_669 = buffer.data(skh + 669);
    const auto *skh_670 = buffer.data(skh + 670);
    const auto *skh_671 = buffer.data(skh + 671);
    const auto *skh_672 = buffer.data(skh + 672);
    const auto *skh_675 = buffer.data(skh + 675);
    const auto *skh_677 = buffer.data(skh + 677);
    const auto *skh_678 = buffer.data(skh + 678);
    const auto *skh_681 = buffer.data(skh + 681);
    const auto *skh_682 = buffer.data(skh + 682);
    const auto *skh_684 = buffer.data(skh + 684);
    const auto *skh_686 = buffer.data(skh + 686);
    const auto *skh_687 = buffer.data(skh + 687);
    const auto *skh_688 = buffer.data(skh + 688);
    const auto *skh_689 = buffer.data(skh + 689);
    const auto *skh_690 = buffer.data(skh + 690);
    const auto *skh_691 = buffer.data(skh + 691);
    const auto *skh_692 = buffer.data(skh + 692);

    const auto *ski1_588 = buffer.data(ski1 + 588);
    const auto *ski1_591 = buffer.data(ski1 + 591);
    const auto *ski1_594 = buffer.data(ski1 + 594);
    const auto *ski1_598 = buffer.data(ski1 + 598);
    const auto *ski1_805 = buffer.data(ski1 + 805);
    const auto *ski1_807 = buffer.data(ski1 + 807);
    const auto *ski1_808 = buffer.data(ski1 + 808);
    const auto *ski1_809 = buffer.data(ski1 + 809);
    const auto *ski1_811 = buffer.data(ski1 + 811);
    const auto *ski1_817 = buffer.data(ski1 + 817);
    const auto *ski1_821 = buffer.data(ski1 + 821);
    const auto *ski1_824 = buffer.data(ski1 + 824);
    const auto *ski1_826 = buffer.data(ski1 + 826);
    const auto *ski1_833 = buffer.data(ski1 + 833);
    const auto *ski1_835 = buffer.data(ski1 + 835);
    const auto *ski1_836 = buffer.data(ski1 + 836);
    const auto *ski1_837 = buffer.data(ski1 + 837);
    const auto *ski1_839 = buffer.data(ski1 + 839);
    const auto *ski1_840 = buffer.data(ski1 + 840);
    const auto *ski1_843 = buffer.data(ski1 + 843);
    const auto *ski1_845 = buffer.data(ski1 + 845);
    const auto *ski1_846 = buffer.data(ski1 + 846);
    const auto *ski1_849 = buffer.data(ski1 + 849);
    const auto *ski1_850 = buffer.data(ski1 + 850);
    const auto *ski1_852 = buffer.data(ski1 + 852);
    const auto *ski1_854 = buffer.data(ski1 + 854);
    const auto *ski1_861 = buffer.data(ski1 + 861);
    const auto *ski1_863 = buffer.data(ski1 + 863);
    const auto *ski1_864 = buffer.data(ski1 + 864);
    const auto *ski1_865 = buffer.data(ski1 + 865);
    const auto *ski1_867 = buffer.data(ski1 + 867);
    const auto *ski1_868 = buffer.data(ski1 + 868);
    const auto *ski1_871 = buffer.data(ski1 + 871);
    const auto *ski1_873 = buffer.data(ski1 + 873);
    const auto *ski1_874 = buffer.data(ski1 + 874);
    const auto *ski1_877 = buffer.data(ski1 + 877);
    const auto *ski1_878 = buffer.data(ski1 + 878);
    const auto *ski1_880 = buffer.data(ski1 + 880);
    const auto *ski1_882 = buffer.data(ski1 + 882);
    const auto *ski1_889 = buffer.data(ski1 + 889);
    const auto *ski1_891 = buffer.data(ski1 + 891);
    const auto *ski1_892 = buffer.data(ski1 + 892);
    const auto *ski1_893 = buffer.data(ski1 + 893);
    const auto *ski1_895 = buffer.data(ski1 + 895);
    const auto *ski1_896 = buffer.data(ski1 + 896);
    const auto *ski1_899 = buffer.data(ski1 + 899);
    const auto *ski1_901 = buffer.data(ski1 + 901);
    const auto *ski1_902 = buffer.data(ski1 + 902);
    const auto *ski1_905 = buffer.data(ski1 + 905);
    const auto *ski1_906 = buffer.data(ski1 + 906);
    const auto *ski1_908 = buffer.data(ski1 + 908);
    const auto *ski1_910 = buffer.data(ski1 + 910);
    const auto *ski1_917 = buffer.data(ski1 + 917);
    const auto *ski1_919 = buffer.data(ski1 + 919);
    const auto *ski1_920 = buffer.data(ski1 + 920);
    const auto *ski1_921 = buffer.data(ski1 + 921);

    const auto *slh_603 = buffer.data(slh + 603);
    const auto *slh_605 = buffer.data(slh + 605);
    const auto *slh_606 = buffer.data(slh + 606);
    const auto *slh_607 = buffer.data(slh + 607);
    const auto *slh_608 = buffer.data(slh + 608);
    const auto *slh_609 = buffer.data(slh + 609);
    const auto *slh_611 = buffer.data(slh + 611);
    const auto *slh_612 = buffer.data(slh + 612);
    const auto *slh_614 = buffer.data(slh + 614);
    const auto *slh_615 = buffer.data(slh + 615);
    const auto *slh_618 = buffer.data(slh + 618);
    const auto *slh_624 = buffer.data(slh + 624);
    const auto *slh_625 = buffer.data(slh + 625);
    const auto *slh_626 = buffer.data(slh + 626);
    const auto *slh_627 = buffer.data(slh + 627);
    const auto *slh_628 = buffer.data(slh + 628);
    const auto *slh_629 = buffer.data(slh + 629);
    const auto *slh_630 = buffer.data(slh + 630);
    const auto *slh_632 = buffer.data(slh + 632);
    const auto *slh_633 = buffer.data(slh + 633);
    const auto *slh_635 = buffer.data(slh + 635);
    const auto *slh_636 = buffer.data(slh + 636);
    const auto *slh_639 = buffer.data(slh + 639);
    const auto *slh_645 = buffer.data(slh + 645);
    const auto *slh_646 = buffer.data(slh + 646);
    const auto *slh_647 = buffer.data(slh + 647);
    const auto *slh_648 = buffer.data(slh + 648);
    const auto *slh_649 = buffer.data(slh + 649);
    const auto *slh_650 = buffer.data(slh + 650);
    const auto *slh_651 = buffer.data(slh + 651);
    const auto *slh_653 = buffer.data(slh + 653);
    const auto *slh_654 = buffer.data(slh + 654);
    const auto *slh_656 = buffer.data(slh + 656);
    const auto *slh_657 = buffer.data(slh + 657);
    const auto *slh_660 = buffer.data(slh + 660);
    const auto *slh_666 = buffer.data(slh + 666);
    const auto *slh_667 = buffer.data(slh + 667);
    const auto *slh_668 = buffer.data(slh + 668);
    const auto *slh_669 = buffer.data(slh + 669);
    const auto *slh_670 = buffer.data(slh + 670);
    const auto *slh_671 = buffer.data(slh + 671);
    const auto *slh_672 = buffer.data(slh + 672);
    const auto *slh_674 = buffer.data(slh + 674);
    const auto *slh_675 = buffer.data(slh + 675);
    const auto *slh_677 = buffer.data(slh + 677);
    const auto *slh_678 = buffer.data(slh + 678);
    const auto *slh_681 = buffer.data(slh + 681);
    const auto *slh_687 = buffer.data(slh + 687);
    const auto *slh_688 = buffer.data(slh + 688);
    const auto *slh_689 = buffer.data(slh + 689);
    const auto *slh_690 = buffer.data(slh + 690);
    const auto *slh_691 = buffer.data(slh + 691);
    const auto *slh_692 = buffer.data(slh + 692);

#pragma omp simd aligned(t_801, t_802, t_803, t_804, pc_x, skh_605, skh_606, skh_607, skh_608, \
                         slh_605, slh_606, slh_607, slh_608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_801[k] = f_11 * skh_605[k]
                   + f_3 * pc_x[k] * slh_605[k];

        t_802[k] = f_11 * skh_606[k]
                   + f_3 * pc_x[k] * slh_606[k];

        t_803[k] = f_11 * skh_607[k]
                   + f_3 * pc_x[k] * slh_607[k];

        t_804[k] = f_11 * skh_608[k]
                   + f_3 * pc_x[k] * slh_608[k];
    }

#pragma omp simd aligned(t_805, t_806, t_807, t_808, pb_x, pc_x, pc_z, ski0_805, ski0_807, \
                         ski0_808, ski1_805, ski1_807, ski1_808, \
                         slh_603 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_805[k] = pb_x[k] * ski0_805[k]
                   - f_10 * pc_x[k] * ski1_805[k];

        t_806[k] = f_3 * pc_z[k] * slh_603[k];

        t_807[k] = pb_x[k] * ski0_807[k]
                   - f_10 * pc_x[k] * ski1_807[k];

        t_808[k] = pb_x[k] * ski0_808[k]
                   - f_10 * pc_x[k] * ski1_808[k];
    }

#pragma omp simd aligned(t_809, t_810, t_811, pb_x, pc_x, pc_y, ski0_809, ski0_811, skh_461, \
                         ski1_809, ski1_811, slh_608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_809[k] = pb_x[k] * ski0_809[k]
                   - f_10 * pc_x[k] * ski1_809[k];

        t_810[k] = f_15 * skh_461[k]
                   + f_3 * pc_y[k] * slh_608[k];

        t_811[k] = pb_x[k] * ski0_811[k]
                   - f_10 * pc_x[k] * ski1_811[k];
    }

#pragma omp simd aligned(t_812, t_813, t_814, t_815, pb_z, pc_y, pc_z, ski0_588, ski0_591, \
                         skh_441, skh_462, ski1_588, ski1_591, \
                         slh_609 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_812[k] = pb_z[k] * ski0_588[k]
                   - f_10 * pc_z[k] * ski1_588[k];

        t_813[k] = f_16 * skh_462[k]
                   + f_3 * pc_y[k] * slh_609[k];

        t_814[k] = f_11 * skh_441[k]
                   + f_3 * pc_z[k] * slh_609[k];

        t_815[k] = pb_z[k] * ski0_591[k]
                   - f_10 * pc_z[k] * ski1_591[k];
    }

#pragma omp simd aligned(t_816, t_817, t_818, pb_x, pb_z, pc_x, pc_y, pc_z, ski0_594, \
                         ski0_817, skh_464, skh_614, ski1_594, ski1_817, \
                         slh_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_816[k] = f_16 * skh_464[k]
                   + f_3 * pc_y[k] * slh_611[k];

        t_817[k] = pb_x[k] * ski0_817[k]
                   + f_14 * skh_614[k]
                   - f_10 * pc_x[k] * ski1_817[k];

        t_818[k] = pb_z[k] * ski0_594[k]
                   - f_10 * pc_z[k] * ski1_594[k];
    }

#pragma omp simd aligned(t_819, t_820, t_821, pb_x, pc_x, pc_y, pc_z, ski0_821, skh_444, \
                         skh_467, skh_618, ski1_821, slh_612, slh_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_819[k] = f_11 * skh_444[k]
                   + f_3 * pc_z[k] * slh_612[k];

        t_820[k] = f_16 * skh_467[k]
                   + f_3 * pc_y[k] * slh_614[k];

        t_821[k] = pb_x[k] * ski0_821[k]
                   + f_13 * skh_618[k]
                   - f_10 * pc_x[k] * ski1_821[k];
    }

#pragma omp simd aligned(t_822, t_823, t_824, pb_x, pb_z, pc_x, pc_z, ski0_598, ski0_824, \
                         skh_447, skh_621, ski1_598, ski1_824, \
                         slh_615 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_822[k] = pb_z[k] * ski0_598[k]
                   - f_10 * pc_z[k] * ski1_598[k];

        t_823[k] = f_11 * skh_447[k]
                   + f_3 * pc_z[k] * slh_615[k];

        t_824[k] = pb_x[k] * ski0_824[k]
                   + f_12 * skh_621[k]
                   - f_10 * pc_x[k] * ski1_824[k];
    }

#pragma omp simd aligned(t_825, t_826, t_827, t_828, pb_x, pc_x, pc_y, ski0_826, skh_471, \
                         skh_623, skh_624, skh_625, ski1_826, slh_618, slh_624, \
                         slh_625 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_825[k] = f_16 * skh_471[k]
                   + f_3 * pc_y[k] * slh_618[k];

        t_826[k] = pb_x[k] * ski0_826[k]
                   + f_12 * skh_623[k]
                   - f_10 * pc_x[k] * ski1_826[k];

        t_827[k] = f_11 * skh_624[k]
                   + f_3 * pc_x[k] * slh_624[k];

        t_828[k] = f_11 * skh_625[k]
                   + f_3 * pc_x[k] * slh_625[k];
    }

#pragma omp simd aligned(t_829, t_830, t_831, t_832, pc_x, skh_626, skh_627, skh_628, skh_629, \
                         slh_626, slh_627, slh_628, slh_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_829[k] = f_11 * skh_626[k]
                   + f_3 * pc_x[k] * slh_626[k];

        t_830[k] = f_11 * skh_627[k]
                   + f_3 * pc_x[k] * slh_627[k];

        t_831[k] = f_11 * skh_628[k]
                   + f_3 * pc_x[k] * slh_628[k];

        t_832[k] = f_11 * skh_629[k]
                   + f_3 * pc_x[k] * slh_629[k];
    }

#pragma omp simd aligned(t_833, t_834, t_835, t_836, pb_x, pc_x, pc_z, ski0_833, ski0_835, \
                         ski0_836, skh_456, ski1_833, ski1_835, ski1_836, \
                         slh_624 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_833[k] = pb_x[k] * ski0_833[k]
                   - f_10 * pc_x[k] * ski1_833[k];

        t_834[k] = f_11 * skh_456[k]
                   + f_3 * pc_z[k] * slh_624[k];

        t_835[k] = pb_x[k] * ski0_835[k]
                   - f_10 * pc_x[k] * ski1_835[k];

        t_836[k] = pb_x[k] * ski0_836[k]
                   - f_10 * pc_x[k] * ski1_836[k];
    }

#pragma omp simd aligned(t_837, t_838, t_839, t_840, pb_x, pc_x, pc_y, ski0_837, ski0_839, \
                         ski0_840, skh_482, skh_630, ski1_837, ski1_839, ski1_840, \
                         slh_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_837[k] = pb_x[k] * ski0_837[k]
                   - f_10 * pc_x[k] * ski1_837[k];

        t_838[k] = f_16 * skh_482[k]
                   + f_3 * pc_y[k] * slh_629[k];

        t_839[k] = pb_x[k] * ski0_839[k]
                   - f_10 * pc_x[k] * ski1_839[k];

        t_840[k] = pb_x[k] * ski0_840[k]
                   + f_16 * skh_630[k]
                   - f_10 * pc_x[k] * ski1_840[k];
    }

#pragma omp simd aligned(t_841, t_842, t_843, t_844, pb_x, pc_x, pc_y, pc_z, ski0_843, \
                         skh_462, skh_483, skh_485, skh_633, ski1_843, slh_630, \
                         slh_632 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_841[k] = f_17 * skh_483[k]
                   + f_3 * pc_y[k] * slh_630[k];

        t_842[k] = f_12 * skh_462[k]
                   + f_3 * pc_z[k] * slh_630[k];

        t_843[k] = pb_x[k] * ski0_843[k]
                   + f_14 * skh_633[k]
                   - f_10 * pc_x[k] * ski1_843[k];

        t_844[k] = f_17 * skh_485[k]
                   + f_3 * pc_y[k] * slh_632[k];
    }

#pragma omp simd aligned(t_845, t_846, t_847, pb_x, pc_x, pc_z, ski0_845, ski0_846, skh_465, \
                         skh_635, skh_636, ski1_845, ski1_846, \
                         slh_633 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_845[k] = pb_x[k] * ski0_845[k]
                   + f_14 * skh_635[k]
                   - f_10 * pc_x[k] * ski1_845[k];

        t_846[k] = pb_x[k] * ski0_846[k]
                   + f_13 * skh_636[k]
                   - f_10 * pc_x[k] * ski1_846[k];

        t_847[k] = f_12 * skh_465[k]
                   + f_3 * pc_z[k] * slh_633[k];
    }

#pragma omp simd aligned(t_848, t_849, t_850, pb_x, pc_x, pc_y, ski0_849, ski0_850, skh_488, \
                         skh_639, skh_640, ski1_849, ski1_850, \
                         slh_635 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_848[k] = f_17 * skh_488[k]
                   + f_3 * pc_y[k] * slh_635[k];

        t_849[k] = pb_x[k] * ski0_849[k]
                   + f_13 * skh_639[k]
                   - f_10 * pc_x[k] * ski1_849[k];

        t_850[k] = pb_x[k] * ski0_850[k]
                   + f_12 * skh_640[k]
                   - f_10 * pc_x[k] * ski1_850[k];
    }

#pragma omp simd aligned(t_851, t_852, t_853, pb_x, pc_x, pc_y, pc_z, ski0_852, skh_468, \
                         skh_492, skh_642, ski1_852, slh_636, slh_639 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_851[k] = f_12 * skh_468[k]
                   + f_3 * pc_z[k] * slh_636[k];

        t_852[k] = pb_x[k] * ski0_852[k]
                   + f_12 * skh_642[k]
                   - f_10 * pc_x[k] * ski1_852[k];

        t_853[k] = f_17 * skh_492[k]
                   + f_3 * pc_y[k] * slh_639[k];
    }

#pragma omp simd aligned(t_854, t_855, t_856, t_857, pb_x, pc_x, ski0_854, skh_644, skh_645, \
                         skh_646, skh_647, ski1_854, slh_645, slh_646, \
                         slh_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_854[k] = pb_x[k] * ski0_854[k]
                   + f_12 * skh_644[k]
                   - f_10 * pc_x[k] * ski1_854[k];

        t_855[k] = f_11 * skh_645[k]
                   + f_3 * pc_x[k] * slh_645[k];

        t_856[k] = f_11 * skh_646[k]
                   + f_3 * pc_x[k] * slh_646[k];

        t_857[k] = f_11 * skh_647[k]
                   + f_3 * pc_x[k] * slh_647[k];
    }

#pragma omp simd aligned(t_858, t_859, t_860, t_861, pb_x, pc_x, ski0_861, skh_648, skh_649, \
                         skh_650, ski1_861, slh_648, slh_649, slh_650 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_858[k] = f_11 * skh_648[k]
                   + f_3 * pc_x[k] * slh_648[k];

        t_859[k] = f_11 * skh_649[k]
                   + f_3 * pc_x[k] * slh_649[k];

        t_860[k] = f_11 * skh_650[k]
                   + f_3 * pc_x[k] * slh_650[k];

        t_861[k] = pb_x[k] * ski0_861[k]
                   - f_10 * pc_x[k] * ski1_861[k];
    }

#pragma omp simd aligned(t_862, t_863, t_864, t_865, pb_x, pc_x, pc_z, ski0_863, ski0_864, \
                         ski0_865, skh_477, ski1_863, ski1_864, ski1_865, \
                         slh_645 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_862[k] = f_12 * skh_477[k]
                   + f_3 * pc_z[k] * slh_645[k];

        t_863[k] = pb_x[k] * ski0_863[k]
                   - f_10 * pc_x[k] * ski1_863[k];

        t_864[k] = pb_x[k] * ski0_864[k]
                   - f_10 * pc_x[k] * ski1_864[k];

        t_865[k] = pb_x[k] * ski0_865[k]
                   - f_10 * pc_x[k] * ski1_865[k];
    }

#pragma omp simd aligned(t_866, t_867, t_868, t_869, pb_x, pc_x, pc_y, ski0_867, ski0_868, \
                         skh_503, skh_504, skh_651, ski1_867, ski1_868, slh_650, \
                         slh_651 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_866[k] = f_17 * skh_503[k]
                   + f_3 * pc_y[k] * slh_650[k];

        t_867[k] = pb_x[k] * ski0_867[k]
                   - f_10 * pc_x[k] * ski1_867[k];

        t_868[k] = pb_x[k] * ski0_868[k]
                   + f_16 * skh_651[k]
                   - f_10 * pc_x[k] * ski1_868[k];

        t_869[k] = f_14 * skh_504[k]
                   + f_3 * pc_y[k] * slh_651[k];
    }

#pragma omp simd aligned(t_870, t_871, t_872, pb_x, pc_x, pc_y, pc_z, ski0_871, skh_483, \
                         skh_506, skh_654, ski1_871, slh_651, slh_653 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_870[k] = f_13 * skh_483[k]
                   + f_3 * pc_z[k] * slh_651[k];

        t_871[k] = pb_x[k] * ski0_871[k]
                   + f_14 * skh_654[k]
                   - f_10 * pc_x[k] * ski1_871[k];

        t_872[k] = f_14 * skh_506[k]
                   + f_3 * pc_y[k] * slh_653[k];
    }

#pragma omp simd aligned(t_873, t_874, t_875, pb_x, pc_x, pc_z, ski0_873, ski0_874, skh_486, \
                         skh_656, skh_657, ski1_873, ski1_874, \
                         slh_654 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_873[k] = pb_x[k] * ski0_873[k]
                   + f_14 * skh_656[k]
                   - f_10 * pc_x[k] * ski1_873[k];

        t_874[k] = pb_x[k] * ski0_874[k]
                   + f_13 * skh_657[k]
                   - f_10 * pc_x[k] * ski1_874[k];

        t_875[k] = f_13 * skh_486[k]
                   + f_3 * pc_z[k] * slh_654[k];
    }

#pragma omp simd aligned(t_876, t_877, t_878, pb_x, pc_x, pc_y, ski0_877, ski0_878, skh_509, \
                         skh_660, skh_661, ski1_877, ski1_878, \
                         slh_656 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_876[k] = f_14 * skh_509[k]
                   + f_3 * pc_y[k] * slh_656[k];

        t_877[k] = pb_x[k] * ski0_877[k]
                   + f_13 * skh_660[k]
                   - f_10 * pc_x[k] * ski1_877[k];

        t_878[k] = pb_x[k] * ski0_878[k]
                   + f_12 * skh_661[k]
                   - f_10 * pc_x[k] * ski1_878[k];
    }

#pragma omp simd aligned(t_879, t_880, t_881, pb_x, pc_x, pc_y, pc_z, ski0_880, skh_489, \
                         skh_513, skh_663, ski1_880, slh_657, slh_660 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_879[k] = f_13 * skh_489[k]
                   + f_3 * pc_z[k] * slh_657[k];

        t_880[k] = pb_x[k] * ski0_880[k]
                   + f_12 * skh_663[k]
                   - f_10 * pc_x[k] * ski1_880[k];

        t_881[k] = f_14 * skh_513[k]
                   + f_3 * pc_y[k] * slh_660[k];
    }

#pragma omp simd aligned(t_882, t_883, t_884, t_885, pb_x, pc_x, ski0_882, skh_665, skh_666, \
                         skh_667, skh_668, ski1_882, slh_666, slh_667, \
                         slh_668 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_882[k] = pb_x[k] * ski0_882[k]
                   + f_12 * skh_665[k]
                   - f_10 * pc_x[k] * ski1_882[k];

        t_883[k] = f_11 * skh_666[k]
                   + f_3 * pc_x[k] * slh_666[k];

        t_884[k] = f_11 * skh_667[k]
                   + f_3 * pc_x[k] * slh_667[k];

        t_885[k] = f_11 * skh_668[k]
                   + f_3 * pc_x[k] * slh_668[k];
    }

#pragma omp simd aligned(t_886, t_887, t_888, t_889, pb_x, pc_x, ski0_889, skh_669, skh_670, \
                         skh_671, ski1_889, slh_669, slh_670, slh_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_886[k] = f_11 * skh_669[k]
                   + f_3 * pc_x[k] * slh_669[k];

        t_887[k] = f_11 * skh_670[k]
                   + f_3 * pc_x[k] * slh_670[k];

        t_888[k] = f_11 * skh_671[k]
                   + f_3 * pc_x[k] * slh_671[k];

        t_889[k] = pb_x[k] * ski0_889[k]
                   - f_10 * pc_x[k] * ski1_889[k];
    }

#pragma omp simd aligned(t_890, t_891, t_892, t_893, pb_x, pc_x, pc_z, ski0_891, ski0_892, \
                         ski0_893, skh_498, ski1_891, ski1_892, ski1_893, \
                         slh_666 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_890[k] = f_13 * skh_498[k]
                   + f_3 * pc_z[k] * slh_666[k];

        t_891[k] = pb_x[k] * ski0_891[k]
                   - f_10 * pc_x[k] * ski1_891[k];

        t_892[k] = pb_x[k] * ski0_892[k]
                   - f_10 * pc_x[k] * ski1_892[k];

        t_893[k] = pb_x[k] * ski0_893[k]
                   - f_10 * pc_x[k] * ski1_893[k];
    }

#pragma omp simd aligned(t_894, t_895, t_896, t_897, pb_x, pc_x, pc_y, ski0_895, ski0_896, \
                         skh_524, skh_525, skh_672, ski1_895, ski1_896, slh_671, \
                         slh_672 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_894[k] = f_14 * skh_524[k]
                   + f_3 * pc_y[k] * slh_671[k];

        t_895[k] = pb_x[k] * ski0_895[k]
                   - f_10 * pc_x[k] * ski1_895[k];

        t_896[k] = pb_x[k] * ski0_896[k]
                   + f_16 * skh_672[k]
                   - f_10 * pc_x[k] * ski1_896[k];

        t_897[k] = f_13 * skh_525[k]
                   + f_3 * pc_y[k] * slh_672[k];
    }

#pragma omp simd aligned(t_898, t_899, t_900, pb_x, pc_x, pc_y, pc_z, ski0_899, skh_504, \
                         skh_527, skh_675, ski1_899, slh_672, slh_674 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_898[k] = f_14 * skh_504[k]
                   + f_3 * pc_z[k] * slh_672[k];

        t_899[k] = pb_x[k] * ski0_899[k]
                   + f_14 * skh_675[k]
                   - f_10 * pc_x[k] * ski1_899[k];

        t_900[k] = f_13 * skh_527[k]
                   + f_3 * pc_y[k] * slh_674[k];
    }

#pragma omp simd aligned(t_901, t_902, t_903, pb_x, pc_x, pc_z, ski0_901, ski0_902, skh_507, \
                         skh_677, skh_678, ski1_901, ski1_902, \
                         slh_675 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_901[k] = pb_x[k] * ski0_901[k]
                   + f_14 * skh_677[k]
                   - f_10 * pc_x[k] * ski1_901[k];

        t_902[k] = pb_x[k] * ski0_902[k]
                   + f_13 * skh_678[k]
                   - f_10 * pc_x[k] * ski1_902[k];

        t_903[k] = f_14 * skh_507[k]
                   + f_3 * pc_z[k] * slh_675[k];
    }

#pragma omp simd aligned(t_904, t_905, t_906, pb_x, pc_x, pc_y, ski0_905, ski0_906, skh_530, \
                         skh_681, skh_682, ski1_905, ski1_906, \
                         slh_677 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_904[k] = f_13 * skh_530[k]
                   + f_3 * pc_y[k] * slh_677[k];

        t_905[k] = pb_x[k] * ski0_905[k]
                   + f_13 * skh_681[k]
                   - f_10 * pc_x[k] * ski1_905[k];

        t_906[k] = pb_x[k] * ski0_906[k]
                   + f_12 * skh_682[k]
                   - f_10 * pc_x[k] * ski1_906[k];
    }

#pragma omp simd aligned(t_907, t_908, t_909, pb_x, pc_x, pc_y, pc_z, ski0_908, skh_510, \
                         skh_534, skh_684, ski1_908, slh_678, slh_681 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_907[k] = f_14 * skh_510[k]
                   + f_3 * pc_z[k] * slh_678[k];

        t_908[k] = pb_x[k] * ski0_908[k]
                   + f_12 * skh_684[k]
                   - f_10 * pc_x[k] * ski1_908[k];

        t_909[k] = f_13 * skh_534[k]
                   + f_3 * pc_y[k] * slh_681[k];
    }

#pragma omp simd aligned(t_910, t_911, t_912, t_913, pb_x, pc_x, ski0_910, skh_686, skh_687, \
                         skh_688, skh_689, ski1_910, slh_687, slh_688, \
                         slh_689 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_910[k] = pb_x[k] * ski0_910[k]
                   + f_12 * skh_686[k]
                   - f_10 * pc_x[k] * ski1_910[k];

        t_911[k] = f_11 * skh_687[k]
                   + f_3 * pc_x[k] * slh_687[k];

        t_912[k] = f_11 * skh_688[k]
                   + f_3 * pc_x[k] * slh_688[k];

        t_913[k] = f_11 * skh_689[k]
                   + f_3 * pc_x[k] * slh_689[k];
    }

#pragma omp simd aligned(t_914, t_915, t_916, t_917, pb_x, pc_x, ski0_917, skh_690, skh_691, \
                         skh_692, ski1_917, slh_690, slh_691, slh_692 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_914[k] = f_11 * skh_690[k]
                   + f_3 * pc_x[k] * slh_690[k];

        t_915[k] = f_11 * skh_691[k]
                   + f_3 * pc_x[k] * slh_691[k];

        t_916[k] = f_11 * skh_692[k]
                   + f_3 * pc_x[k] * slh_692[k];

        t_917[k] = pb_x[k] * ski0_917[k]
                   - f_10 * pc_x[k] * ski1_917[k];
    }

#pragma omp simd aligned(t_918, t_919, t_920, t_921, pb_x, pc_x, pc_z, ski0_919, ski0_920, \
                         ski0_921, skh_519, ski1_919, ski1_920, ski1_921, \
                         slh_687 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_918[k] = f_14 * skh_519[k]
                   + f_3 * pc_z[k] * slh_687[k];

        t_919[k] = pb_x[k] * ski0_919[k]
                   - f_10 * pc_x[k] * ski1_919[k];

        t_920[k] = pb_x[k] * ski0_920[k]
                   - f_10 * pc_x[k] * ski1_920[k];

        t_921[k] = pb_x[k] * ski0_921[k]
                   - f_10 * pc_x[k] * ski1_921[k];
    }
}

static auto
compute_prim_sli_three_center_electron_repulsion_0_piece8(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t ski0,
                                                          const size_t skh, const size_t ski1,
                                                          const size_t slg0, const size_t slg1,
                                                          const size_t slh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.5 / gamma;
    const auto f_5 = 1.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 0.5 / gamma;
    const auto f_9 = 0.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 3.5 / q;
    const auto f_16 = 3.0 / q;
    const auto f_17 = 2.5 / q;

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
    auto *t_1037 = buffer.data(target + 1037);
    auto *t_1038 = buffer.data(target + 1038);
    auto *t_1039 = buffer.data(target + 1039);
    auto *t_1040 = buffer.data(target + 1040);
    auto *t_1041 = buffer.data(target + 1041);
    auto *t_1042 = buffer.data(target + 1042);
    auto *t_1043 = buffer.data(target + 1043);
    auto *t_1044 = buffer.data(target + 1044);
    auto *t_1045 = buffer.data(target + 1045);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ski0_756 = buffer.data(ski0 + 756);
    const auto *ski0_761 = buffer.data(ski0 + 761);
    const auto *ski0_765 = buffer.data(ski0 + 765);
    const auto *ski0_770 = buffer.data(ski0 + 770);
    const auto *ski0_784 = buffer.data(ski0 + 784);
    const auto *ski0_787 = buffer.data(ski0 + 787);
    const auto *ski0_790 = buffer.data(ski0 + 790);
    const auto *ski0_923 = buffer.data(ski0 + 923);
    const auto *ski0_924 = buffer.data(ski0 + 924);
    const auto *ski0_927 = buffer.data(ski0 + 927);
    const auto *ski0_929 = buffer.data(ski0 + 929);
    const auto *ski0_930 = buffer.data(ski0 + 930);
    const auto *ski0_933 = buffer.data(ski0 + 933);
    const auto *ski0_934 = buffer.data(ski0 + 934);
    const auto *ski0_936 = buffer.data(ski0 + 936);
    const auto *ski0_938 = buffer.data(ski0 + 938);
    const auto *ski0_945 = buffer.data(ski0 + 945);
    const auto *ski0_947 = buffer.data(ski0 + 947);
    const auto *ski0_948 = buffer.data(ski0 + 948);
    const auto *ski0_949 = buffer.data(ski0 + 949);
    const auto *ski0_951 = buffer.data(ski0 + 951);
    const auto *ski0_955 = buffer.data(ski0 + 955);
    const auto *ski0_958 = buffer.data(ski0 + 958);
    const auto *ski0_962 = buffer.data(ski0 + 962);
    const auto *ski0_964 = buffer.data(ski0 + 964);
    const auto *ski0_973 = buffer.data(ski0 + 973);
    const auto *ski0_975 = buffer.data(ski0 + 975);
    const auto *ski0_976 = buffer.data(ski0 + 976);
    const auto *ski0_977 = buffer.data(ski0 + 977);
    const auto *ski0_979 = buffer.data(ski0 + 979);
    const auto *ski0_980 = buffer.data(ski0 + 980);
    const auto *ski0_983 = buffer.data(ski0 + 983);
    const auto *ski0_985 = buffer.data(ski0 + 985);
    const auto *ski0_986 = buffer.data(ski0 + 986);
    const auto *ski0_989 = buffer.data(ski0 + 989);
    const auto *ski0_990 = buffer.data(ski0 + 990);
    const auto *ski0_992 = buffer.data(ski0 + 992);
    const auto *ski0_994 = buffer.data(ski0 + 994);
    const auto *ski0_1001 = buffer.data(ski0 + 1001);
    const auto *ski0_1003 = buffer.data(ski0 + 1003);
    const auto *ski0_1004 = buffer.data(ski0 + 1004);
    const auto *ski0_1005 = buffer.data(ski0 + 1005);
    const auto *ski0_1007 = buffer.data(ski0 + 1007);

    const auto *skh_525 = buffer.data(skh + 525);
    const auto *skh_528 = buffer.data(skh + 528);
    const auto *skh_531 = buffer.data(skh + 531);
    const auto *skh_540 = buffer.data(skh + 540);
    const auto *skh_545 = buffer.data(skh + 545);
    const auto *skh_546 = buffer.data(skh + 546);
    const auto *skh_548 = buffer.data(skh + 548);
    const auto *skh_549 = buffer.data(skh + 549);
    const auto *skh_551 = buffer.data(skh + 551);
    const auto *skh_552 = buffer.data(skh + 552);
    const auto *skh_555 = buffer.data(skh + 555);
    const auto *skh_561 = buffer.data(skh + 561);
    const auto *skh_566 = buffer.data(skh + 566);
    const auto *skh_567 = buffer.data(skh + 567);
    const auto *skh_569 = buffer.data(skh + 569);
    const auto *skh_570 = buffer.data(skh + 570);
    const auto *skh_572 = buffer.data(skh + 572);
    const auto *skh_573 = buffer.data(skh + 573);
    const auto *skh_576 = buffer.data(skh + 576);
    const auto *skh_582 = buffer.data(skh + 582);
    const auto *skh_587 = buffer.data(skh + 587);
    const auto *skh_588 = buffer.data(skh + 588);
    const auto *skh_590 = buffer.data(skh + 590);
    const auto *skh_591 = buffer.data(skh + 591);
    const auto *skh_593 = buffer.data(skh + 593);
    const auto *skh_597 = buffer.data(skh + 597);
    const auto *skh_603 = buffer.data(skh + 603);
    const auto *skh_605 = buffer.data(skh + 605);
    const auto *skh_606 = buffer.data(skh + 606);
    const auto *skh_607 = buffer.data(skh + 607);
    const auto *skh_608 = buffer.data(skh + 608);
    const auto *skh_609 = buffer.data(skh + 609);
    const auto *skh_611 = buffer.data(skh + 611);
    const auto *skh_614 = buffer.data(skh + 614);
    const auto *skh_693 = buffer.data(skh + 693);
    const auto *skh_696 = buffer.data(skh + 696);
    const auto *skh_698 = buffer.data(skh + 698);
    const auto *skh_699 = buffer.data(skh + 699);
    const auto *skh_702 = buffer.data(skh + 702);
    const auto *skh_703 = buffer.data(skh + 703);
    const auto *skh_705 = buffer.data(skh + 705);
    const auto *skh_707 = buffer.data(skh + 707);
    const auto *skh_708 = buffer.data(skh + 708);
    const auto *skh_709 = buffer.data(skh + 709);
    const auto *skh_710 = buffer.data(skh + 710);
    const auto *skh_711 = buffer.data(skh + 711);
    const auto *skh_712 = buffer.data(skh + 712);
    const auto *skh_713 = buffer.data(skh + 713);
    const auto *skh_717 = buffer.data(skh + 717);
    const auto *skh_720 = buffer.data(skh + 720);
    const auto *skh_724 = buffer.data(skh + 724);
    const auto *skh_726 = buffer.data(skh + 726);
    const auto *skh_729 = buffer.data(skh + 729);
    const auto *skh_730 = buffer.data(skh + 730);
    const auto *skh_731 = buffer.data(skh + 731);
    const auto *skh_732 = buffer.data(skh + 732);
    const auto *skh_733 = buffer.data(skh + 733);
    const auto *skh_734 = buffer.data(skh + 734);
    const auto *skh_735 = buffer.data(skh + 735);
    const auto *skh_738 = buffer.data(skh + 738);
    const auto *skh_740 = buffer.data(skh + 740);
    const auto *skh_741 = buffer.data(skh + 741);
    const auto *skh_744 = buffer.data(skh + 744);
    const auto *skh_745 = buffer.data(skh + 745);
    const auto *skh_747 = buffer.data(skh + 747);
    const auto *skh_749 = buffer.data(skh + 749);
    const auto *skh_750 = buffer.data(skh + 750);
    const auto *skh_751 = buffer.data(skh + 751);
    const auto *skh_752 = buffer.data(skh + 752);
    const auto *skh_753 = buffer.data(skh + 753);
    const auto *skh_754 = buffer.data(skh + 754);
    const auto *skh_755 = buffer.data(skh + 755);

    const auto *ski1_756 = buffer.data(ski1 + 756);
    const auto *ski1_761 = buffer.data(ski1 + 761);
    const auto *ski1_765 = buffer.data(ski1 + 765);
    const auto *ski1_770 = buffer.data(ski1 + 770);
    const auto *ski1_784 = buffer.data(ski1 + 784);
    const auto *ski1_787 = buffer.data(ski1 + 787);
    const auto *ski1_790 = buffer.data(ski1 + 790);
    const auto *ski1_923 = buffer.data(ski1 + 923);
    const auto *ski1_924 = buffer.data(ski1 + 924);
    const auto *ski1_927 = buffer.data(ski1 + 927);
    const auto *ski1_929 = buffer.data(ski1 + 929);
    const auto *ski1_930 = buffer.data(ski1 + 930);
    const auto *ski1_933 = buffer.data(ski1 + 933);
    const auto *ski1_934 = buffer.data(ski1 + 934);
    const auto *ski1_936 = buffer.data(ski1 + 936);
    const auto *ski1_938 = buffer.data(ski1 + 938);
    const auto *ski1_945 = buffer.data(ski1 + 945);
    const auto *ski1_947 = buffer.data(ski1 + 947);
    const auto *ski1_948 = buffer.data(ski1 + 948);
    const auto *ski1_949 = buffer.data(ski1 + 949);
    const auto *ski1_951 = buffer.data(ski1 + 951);
    const auto *ski1_955 = buffer.data(ski1 + 955);
    const auto *ski1_958 = buffer.data(ski1 + 958);
    const auto *ski1_962 = buffer.data(ski1 + 962);
    const auto *ski1_964 = buffer.data(ski1 + 964);
    const auto *ski1_973 = buffer.data(ski1 + 973);
    const auto *ski1_975 = buffer.data(ski1 + 975);
    const auto *ski1_976 = buffer.data(ski1 + 976);
    const auto *ski1_977 = buffer.data(ski1 + 977);
    const auto *ski1_979 = buffer.data(ski1 + 979);
    const auto *ski1_980 = buffer.data(ski1 + 980);
    const auto *ski1_983 = buffer.data(ski1 + 983);
    const auto *ski1_985 = buffer.data(ski1 + 985);
    const auto *ski1_986 = buffer.data(ski1 + 986);
    const auto *ski1_989 = buffer.data(ski1 + 989);
    const auto *ski1_990 = buffer.data(ski1 + 990);
    const auto *ski1_992 = buffer.data(ski1 + 992);
    const auto *ski1_994 = buffer.data(ski1 + 994);
    const auto *ski1_1001 = buffer.data(ski1 + 1001);
    const auto *ski1_1003 = buffer.data(ski1 + 1003);
    const auto *ski1_1004 = buffer.data(ski1 + 1004);
    const auto *ski1_1005 = buffer.data(ski1 + 1005);
    const auto *ski1_1007 = buffer.data(ski1 + 1007);

    const auto *slg0_540 = buffer.data(slg0 + 540);
    const auto *slg0_543 = buffer.data(slg0 + 543);
    const auto *slg0_545 = buffer.data(slg0 + 545);
    const auto *slg0_546 = buffer.data(slg0 + 546);
    const auto *slg0_549 = buffer.data(slg0 + 549);
    const auto *slg0_550 = buffer.data(slg0 + 550);
    const auto *slg0_552 = buffer.data(slg0 + 552);
    const auto *slg0_553 = buffer.data(slg0 + 553);
    const auto *slg0_554 = buffer.data(slg0 + 554);
    const auto *slg0_560 = buffer.data(slg0 + 560);
    const auto *slg0_564 = buffer.data(slg0 + 564);

    const auto *slg1_540 = buffer.data(slg1 + 540);
    const auto *slg1_543 = buffer.data(slg1 + 543);
    const auto *slg1_545 = buffer.data(slg1 + 545);
    const auto *slg1_546 = buffer.data(slg1 + 546);
    const auto *slg1_549 = buffer.data(slg1 + 549);
    const auto *slg1_550 = buffer.data(slg1 + 550);
    const auto *slg1_552 = buffer.data(slg1 + 552);
    const auto *slg1_553 = buffer.data(slg1 + 553);
    const auto *slg1_554 = buffer.data(slg1 + 554);
    const auto *slg1_560 = buffer.data(slg1 + 560);
    const auto *slg1_564 = buffer.data(slg1 + 564);

    const auto *slh_692 = buffer.data(slh + 692);
    const auto *slh_693 = buffer.data(slh + 693);
    const auto *slh_695 = buffer.data(slh + 695);
    const auto *slh_696 = buffer.data(slh + 696);
    const auto *slh_698 = buffer.data(slh + 698);
    const auto *slh_699 = buffer.data(slh + 699);
    const auto *slh_702 = buffer.data(slh + 702);
    const auto *slh_708 = buffer.data(slh + 708);
    const auto *slh_709 = buffer.data(slh + 709);
    const auto *slh_710 = buffer.data(slh + 710);
    const auto *slh_711 = buffer.data(slh + 711);
    const auto *slh_712 = buffer.data(slh + 712);
    const auto *slh_713 = buffer.data(slh + 713);
    const auto *slh_714 = buffer.data(slh + 714);
    const auto *slh_716 = buffer.data(slh + 716);
    const auto *slh_717 = buffer.data(slh + 717);
    const auto *slh_719 = buffer.data(slh + 719);
    const auto *slh_720 = buffer.data(slh + 720);
    const auto *slh_723 = buffer.data(slh + 723);
    const auto *slh_729 = buffer.data(slh + 729);
    const auto *slh_730 = buffer.data(slh + 730);
    const auto *slh_731 = buffer.data(slh + 731);
    const auto *slh_732 = buffer.data(slh + 732);
    const auto *slh_733 = buffer.data(slh + 733);
    const auto *slh_734 = buffer.data(slh + 734);
    const auto *slh_735 = buffer.data(slh + 735);
    const auto *slh_737 = buffer.data(slh + 737);
    const auto *slh_738 = buffer.data(slh + 738);
    const auto *slh_740 = buffer.data(slh + 740);
    const auto *slh_741 = buffer.data(slh + 741);
    const auto *slh_744 = buffer.data(slh + 744);
    const auto *slh_750 = buffer.data(slh + 750);
    const auto *slh_751 = buffer.data(slh + 751);
    const auto *slh_752 = buffer.data(slh + 752);
    const auto *slh_753 = buffer.data(slh + 753);
    const auto *slh_754 = buffer.data(slh + 754);
    const auto *slh_755 = buffer.data(slh + 755);
    const auto *slh_756 = buffer.data(slh + 756);
    const auto *slh_758 = buffer.data(slh + 758);
    const auto *slh_759 = buffer.data(slh + 759);
    const auto *slh_761 = buffer.data(slh + 761);
    const auto *slh_762 = buffer.data(slh + 762);
    const auto *slh_765 = buffer.data(slh + 765);
    const auto *slh_766 = buffer.data(slh + 766);
    const auto *slh_768 = buffer.data(slh + 768);
    const auto *slh_770 = buffer.data(slh + 770);
    const auto *slh_771 = buffer.data(slh + 771);
    const auto *slh_772 = buffer.data(slh + 772);
    const auto *slh_773 = buffer.data(slh + 773);
    const auto *slh_774 = buffer.data(slh + 774);
    const auto *slh_775 = buffer.data(slh + 775);
    const auto *slh_776 = buffer.data(slh + 776);
    const auto *slh_777 = buffer.data(slh + 777);
    const auto *slh_779 = buffer.data(slh + 779);
    const auto *slh_780 = buffer.data(slh + 780);
    const auto *slh_782 = buffer.data(slh + 782);
    const auto *slh_786 = buffer.data(slh + 786);

#pragma omp simd aligned(t_922, t_923, t_924, t_925, pb_x, pc_x, pc_y, ski0_923, ski0_924, \
                         skh_545, skh_546, skh_693, ski1_923, ski1_924, slh_692, \
                         slh_693 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_922[k] = f_13 * skh_545[k]
                   + f_3 * pc_y[k] * slh_692[k];

        t_923[k] = pb_x[k] * ski0_923[k]
                   - f_10 * pc_x[k] * ski1_923[k];

        t_924[k] = pb_x[k] * ski0_924[k]
                   + f_16 * skh_693[k]
                   - f_10 * pc_x[k] * ski1_924[k];

        t_925[k] = f_12 * skh_546[k]
                   + f_3 * pc_y[k] * slh_693[k];
    }

#pragma omp simd aligned(t_926, t_927, t_928, pb_x, pc_x, pc_y, pc_z, ski0_927, skh_525, \
                         skh_548, skh_696, ski1_927, slh_693, slh_695 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_926[k] = f_17 * skh_525[k]
                   + f_3 * pc_z[k] * slh_693[k];

        t_927[k] = pb_x[k] * ski0_927[k]
                   + f_14 * skh_696[k]
                   - f_10 * pc_x[k] * ski1_927[k];

        t_928[k] = f_12 * skh_548[k]
                   + f_3 * pc_y[k] * slh_695[k];
    }

#pragma omp simd aligned(t_929, t_930, t_931, pb_x, pc_x, pc_z, ski0_929, ski0_930, skh_528, \
                         skh_698, skh_699, ski1_929, ski1_930, \
                         slh_696 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_929[k] = pb_x[k] * ski0_929[k]
                   + f_14 * skh_698[k]
                   - f_10 * pc_x[k] * ski1_929[k];

        t_930[k] = pb_x[k] * ski0_930[k]
                   + f_13 * skh_699[k]
                   - f_10 * pc_x[k] * ski1_930[k];

        t_931[k] = f_17 * skh_528[k]
                   + f_3 * pc_z[k] * slh_696[k];
    }

#pragma omp simd aligned(t_932, t_933, t_934, pb_x, pc_x, pc_y, ski0_933, ski0_934, skh_551, \
                         skh_702, skh_703, ski1_933, ski1_934, \
                         slh_698 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_932[k] = f_12 * skh_551[k]
                   + f_3 * pc_y[k] * slh_698[k];

        t_933[k] = pb_x[k] * ski0_933[k]
                   + f_13 * skh_702[k]
                   - f_10 * pc_x[k] * ski1_933[k];

        t_934[k] = pb_x[k] * ski0_934[k]
                   + f_12 * skh_703[k]
                   - f_10 * pc_x[k] * ski1_934[k];
    }

#pragma omp simd aligned(t_935, t_936, t_937, pb_x, pc_x, pc_y, pc_z, ski0_936, skh_531, \
                         skh_555, skh_705, ski1_936, slh_699, slh_702 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_935[k] = f_17 * skh_531[k]
                   + f_3 * pc_z[k] * slh_699[k];

        t_936[k] = pb_x[k] * ski0_936[k]
                   + f_12 * skh_705[k]
                   - f_10 * pc_x[k] * ski1_936[k];

        t_937[k] = f_12 * skh_555[k]
                   + f_3 * pc_y[k] * slh_702[k];
    }

#pragma omp simd aligned(t_938, t_939, t_940, t_941, pb_x, pc_x, ski0_938, skh_707, skh_708, \
                         skh_709, skh_710, ski1_938, slh_708, slh_709, \
                         slh_710 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_938[k] = pb_x[k] * ski0_938[k]
                   + f_12 * skh_707[k]
                   - f_10 * pc_x[k] * ski1_938[k];

        t_939[k] = f_11 * skh_708[k]
                   + f_3 * pc_x[k] * slh_708[k];

        t_940[k] = f_11 * skh_709[k]
                   + f_3 * pc_x[k] * slh_709[k];

        t_941[k] = f_11 * skh_710[k]
                   + f_3 * pc_x[k] * slh_710[k];
    }

#pragma omp simd aligned(t_942, t_943, t_944, t_945, pb_x, pc_x, ski0_945, skh_711, skh_712, \
                         skh_713, ski1_945, slh_711, slh_712, slh_713 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_942[k] = f_11 * skh_711[k]
                   + f_3 * pc_x[k] * slh_711[k];

        t_943[k] = f_11 * skh_712[k]
                   + f_3 * pc_x[k] * slh_712[k];

        t_944[k] = f_11 * skh_713[k]
                   + f_3 * pc_x[k] * slh_713[k];

        t_945[k] = pb_x[k] * ski0_945[k]
                   - f_10 * pc_x[k] * ski1_945[k];
    }

#pragma omp simd aligned(t_946, t_947, t_948, t_949, pb_x, pc_x, pc_z, ski0_947, ski0_948, \
                         ski0_949, skh_540, ski1_947, ski1_948, ski1_949, \
                         slh_708 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_946[k] = f_17 * skh_540[k]
                   + f_3 * pc_z[k] * slh_708[k];

        t_947[k] = pb_x[k] * ski0_947[k]
                   - f_10 * pc_x[k] * ski1_947[k];

        t_948[k] = pb_x[k] * ski0_948[k]
                   - f_10 * pc_x[k] * ski1_948[k];

        t_949[k] = pb_x[k] * ski0_949[k]
                   - f_10 * pc_x[k] * ski1_949[k];
    }

#pragma omp simd aligned(t_950, t_951, t_952, t_953, pb_x, pb_y, pc_x, pc_y, ski0_756, \
                         ski0_951, skh_566, skh_567, ski1_756, ski1_951, slh_713, \
                         slh_714 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_950[k] = f_12 * skh_566[k]
                   + f_3 * pc_y[k] * slh_713[k];

        t_951[k] = pb_x[k] * ski0_951[k]
                   - f_10 * pc_x[k] * ski1_951[k];

        t_952[k] = pb_y[k] * ski0_756[k]
                   - f_10 * pc_y[k] * ski1_756[k];

        t_953[k] = f_11 * skh_567[k]
                   + f_3 * pc_y[k] * slh_714[k];
    }

#pragma omp simd aligned(t_954, t_955, t_956, pb_x, pc_x, pc_y, pc_z, ski0_955, skh_546, \
                         skh_569, skh_717, ski1_955, slh_714, slh_716 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_954[k] = f_16 * skh_546[k]
                   + f_3 * pc_z[k] * slh_714[k];

        t_955[k] = pb_x[k] * ski0_955[k]
                   + f_14 * skh_717[k]
                   - f_10 * pc_x[k] * ski1_955[k];

        t_956[k] = f_11 * skh_569[k]
                   + f_3 * pc_y[k] * slh_716[k];
    }

#pragma omp simd aligned(t_957, t_958, t_959, pb_x, pb_y, pc_x, pc_y, pc_z, ski0_761, \
                         ski0_958, skh_549, skh_720, ski1_761, ski1_958, \
                         slh_717 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_957[k] = pb_y[k] * ski0_761[k]
                   - f_10 * pc_y[k] * ski1_761[k];

        t_958[k] = pb_x[k] * ski0_958[k]
                   + f_13 * skh_720[k]
                   - f_10 * pc_x[k] * ski1_958[k];

        t_959[k] = f_16 * skh_549[k]
                   + f_3 * pc_z[k] * slh_717[k];
    }

#pragma omp simd aligned(t_960, t_961, t_962, pb_x, pb_y, pc_x, pc_y, ski0_765, ski0_962, \
                         skh_572, skh_724, ski1_765, ski1_962, \
                         slh_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_960[k] = f_11 * skh_572[k]
                   + f_3 * pc_y[k] * slh_719[k];

        t_961[k] = pb_y[k] * ski0_765[k]
                   - f_10 * pc_y[k] * ski1_765[k];

        t_962[k] = pb_x[k] * ski0_962[k]
                   + f_12 * skh_724[k]
                   - f_10 * pc_x[k] * ski1_962[k];
    }

#pragma omp simd aligned(t_963, t_964, t_965, pb_x, pc_x, pc_y, pc_z, ski0_964, skh_552, \
                         skh_576, skh_726, ski1_964, slh_720, slh_723 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_963[k] = f_16 * skh_552[k]
                   + f_3 * pc_z[k] * slh_720[k];

        t_964[k] = pb_x[k] * ski0_964[k]
                   + f_12 * skh_726[k]
                   - f_10 * pc_x[k] * ski1_964[k];

        t_965[k] = f_11 * skh_576[k]
                   + f_3 * pc_y[k] * slh_723[k];
    }

#pragma omp simd aligned(t_966, t_967, t_968, t_969, pb_y, pc_x, pc_y, ski0_770, skh_729, \
                         skh_730, skh_731, ski1_770, slh_729, slh_730, \
                         slh_731 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_966[k] = pb_y[k] * ski0_770[k]
                   - f_10 * pc_y[k] * ski1_770[k];

        t_967[k] = f_11 * skh_729[k]
                   + f_3 * pc_x[k] * slh_729[k];

        t_968[k] = f_11 * skh_730[k]
                   + f_3 * pc_x[k] * slh_730[k];

        t_969[k] = f_11 * skh_731[k]
                   + f_3 * pc_x[k] * slh_731[k];
    }

#pragma omp simd aligned(t_970, t_971, t_972, t_973, pb_x, pc_x, ski0_973, skh_732, skh_733, \
                         skh_734, ski1_973, slh_732, slh_733, slh_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_970[k] = f_11 * skh_732[k]
                   + f_3 * pc_x[k] * slh_732[k];

        t_971[k] = f_11 * skh_733[k]
                   + f_3 * pc_x[k] * slh_733[k];

        t_972[k] = f_11 * skh_734[k]
                   + f_3 * pc_x[k] * slh_734[k];

        t_973[k] = pb_x[k] * ski0_973[k]
                   - f_10 * pc_x[k] * ski1_973[k];
    }

#pragma omp simd aligned(t_974, t_975, t_976, t_977, pb_x, pc_x, pc_z, ski0_975, ski0_976, \
                         ski0_977, skh_561, ski1_975, ski1_976, ski1_977, \
                         slh_729 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_974[k] = f_16 * skh_561[k]
                   + f_3 * pc_z[k] * slh_729[k];

        t_975[k] = pb_x[k] * ski0_975[k]
                   - f_10 * pc_x[k] * ski1_975[k];

        t_976[k] = pb_x[k] * ski0_976[k]
                   - f_10 * pc_x[k] * ski1_976[k];

        t_977[k] = pb_x[k] * ski0_977[k]
                   - f_10 * pc_x[k] * ski1_977[k];
    }

#pragma omp simd aligned(t_978, t_979, t_980, t_981, pb_x, pc_x, pc_y, ski0_979, ski0_980, \
                         skh_587, skh_735, ski1_979, ski1_980, slh_734, \
                         slh_735 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_978[k] = f_11 * skh_587[k]
                   + f_3 * pc_y[k] * slh_734[k];

        t_979[k] = pb_x[k] * ski0_979[k]
                   - f_10 * pc_x[k] * ski1_979[k];

        t_980[k] = pb_x[k] * ski0_980[k]
                   + f_16 * skh_735[k]
                   - f_10 * pc_x[k] * ski1_980[k];

        t_981[k] = f_3 * pc_y[k] * slh_735[k];
    }

#pragma omp simd aligned(t_982, t_983, t_984, pb_x, pc_x, pc_y, pc_z, ski0_983, skh_567, \
                         skh_738, ski1_983, slh_735, slh_737 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_982[k] = f_15 * skh_567[k]
                   + f_3 * pc_z[k] * slh_735[k];

        t_983[k] = pb_x[k] * ski0_983[k]
                   + f_14 * skh_738[k]
                   - f_10 * pc_x[k] * ski1_983[k];

        t_984[k] = f_3 * pc_y[k] * slh_737[k];
    }

#pragma omp simd aligned(t_985, t_986, t_987, pb_x, pc_x, pc_z, ski0_985, ski0_986, skh_570, \
                         skh_740, skh_741, ski1_985, ski1_986, \
                         slh_738 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_985[k] = pb_x[k] * ski0_985[k]
                   + f_14 * skh_740[k]
                   - f_10 * pc_x[k] * ski1_985[k];

        t_986[k] = pb_x[k] * ski0_986[k]
                   + f_13 * skh_741[k]
                   - f_10 * pc_x[k] * ski1_986[k];

        t_987[k] = f_15 * skh_570[k]
                   + f_3 * pc_z[k] * slh_738[k];
    }

#pragma omp simd aligned(t_988, t_989, t_990, pb_x, pc_x, pc_y, ski0_989, ski0_990, skh_744, \
                         skh_745, ski1_989, ski1_990, slh_740 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_988[k] = f_3 * pc_y[k] * slh_740[k];

        t_989[k] = pb_x[k] * ski0_989[k]
                   + f_13 * skh_744[k]
                   - f_10 * pc_x[k] * ski1_989[k];

        t_990[k] = pb_x[k] * ski0_990[k]
                   + f_12 * skh_745[k]
                   - f_10 * pc_x[k] * ski1_990[k];
    }

#pragma omp simd aligned(t_991, t_992, t_993, pb_x, pc_x, pc_y, pc_z, ski0_992, skh_573, \
                         skh_747, ski1_992, slh_741, slh_744 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_991[k] = f_15 * skh_573[k]
                   + f_3 * pc_z[k] * slh_741[k];

        t_992[k] = pb_x[k] * ski0_992[k]
                   + f_12 * skh_747[k]
                   - f_10 * pc_x[k] * ski1_992[k];

        t_993[k] = f_3 * pc_y[k] * slh_744[k];
    }

#pragma omp simd aligned(t_994, t_995, t_996, t_997, pb_x, pc_x, ski0_994, skh_749, skh_750, \
                         skh_751, skh_752, ski1_994, slh_750, slh_751, \
                         slh_752 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_994[k] = pb_x[k] * ski0_994[k]
                   + f_12 * skh_749[k]
                   - f_10 * pc_x[k] * ski1_994[k];

        t_995[k] = f_11 * skh_750[k]
                   + f_3 * pc_x[k] * slh_750[k];

        t_996[k] = f_11 * skh_751[k]
                   + f_3 * pc_x[k] * slh_751[k];

        t_997[k] = f_11 * skh_752[k]
                   + f_3 * pc_x[k] * slh_752[k];
    }

#pragma omp simd aligned(t_998, t_999, t_1000, t_1001, pb_x, pc_x, ski0_1001, skh_753, \
                         skh_754, skh_755, ski1_1001, slh_753, slh_754, \
                         slh_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_998[k] = f_11 * skh_753[k]
                   + f_3 * pc_x[k] * slh_753[k];

        t_999[k] = f_11 * skh_754[k]
                   + f_3 * pc_x[k] * slh_754[k];

        t_1000[k] = f_11 * skh_755[k]
                    + f_3 * pc_x[k] * slh_755[k];

        t_1001[k] = pb_x[k] * ski0_1001[k]
                    - f_10 * pc_x[k] * ski1_1001[k];
    }

#pragma omp simd aligned(t_1002, t_1003, t_1004, t_1005, pb_x, pc_x, pc_z, ski0_1003, \
                         ski0_1004, ski0_1005, skh_582, ski1_1003, ski1_1004, ski1_1005, \
                         slh_750 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1002[k] = f_15 * skh_582[k]
                    + f_3 * pc_z[k] * slh_750[k];

        t_1003[k] = pb_x[k] * ski0_1003[k]
                    - f_10 * pc_x[k] * ski1_1003[k];

        t_1004[k] = pb_x[k] * ski0_1004[k]
                    - f_10 * pc_x[k] * ski1_1004[k];

        t_1005[k] = pb_x[k] * ski0_1005[k]
                    - f_10 * pc_x[k] * ski1_1005[k];
    }

#pragma omp simd aligned(t_1006, t_1007, t_1008, t_1009, t_1010, pb_x, pc_x, pc_y, pc_z, \
                         ski0_1007, skh_588, ski1_1007, slg0_540, slg1_540, slh_755, \
                         slh_756 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1006[k] = f_3 * pc_y[k] * slh_755[k];

        t_1007[k] = pb_x[k] * ski0_1007[k]
                    - f_10 * pc_x[k] * ski1_1007[k];

        t_1008[k] = f_1 * slg0_540[k]
                    - f_2 * slg1_540[k]
                    + f_3 * pc_x[k] * slh_756[k];

        t_1009[k] = f_0 * skh_588[k]
                    + f_3 * pc_y[k] * slh_756[k];

        t_1010[k] = f_3 * pc_z[k] * slh_756[k];
    }

#pragma omp simd aligned(t_1011, t_1012, t_1013, pc_x, pc_y, skh_590, slg0_543, slg0_545, \
                         slg1_543, slg1_545, slh_758, slh_759, \
                         slh_761 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1011[k] = f_4 * slg0_543[k]
                    - f_5 * slg1_543[k]
                    + f_3 * pc_x[k] * slh_759[k];

        t_1012[k] = f_0 * skh_590[k]
                    + f_3 * pc_y[k] * slh_758[k];

        t_1013[k] = f_4 * slg0_545[k]
                    - f_5 * slg1_545[k]
                    + f_3 * pc_x[k] * slh_761[k];
    }

#pragma omp simd aligned(t_1014, t_1015, t_1016, t_1017, pc_x, pc_y, pc_z, skh_593, slg0_546, \
                         slg0_549, slg1_546, slg1_549, slh_759, slh_761, slh_762, \
                         slh_765 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1014[k] = f_6 * slg0_546[k]
                    - f_7 * slg1_546[k]
                    + f_3 * pc_x[k] * slh_762[k];

        t_1015[k] = f_3 * pc_z[k] * slh_759[k];

        t_1016[k] = f_0 * skh_593[k]
                    + f_3 * pc_y[k] * slh_761[k];

        t_1017[k] = f_6 * slg0_549[k]
                    - f_7 * slg1_549[k]
                    + f_3 * pc_x[k] * slh_765[k];
    }

#pragma omp simd aligned(t_1018, t_1019, t_1020, t_1021, pc_x, pc_y, pc_z, skh_597, slg0_550, \
                         slg0_552, slg1_550, slg1_552, slh_762, slh_765, slh_766, \
                         slh_768 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1018[k] = f_8 * slg0_550[k]
                    - f_9 * slg1_550[k]
                    + f_3 * pc_x[k] * slh_766[k];

        t_1019[k] = f_3 * pc_z[k] * slh_762[k];

        t_1020[k] = f_8 * slg0_552[k]
                    - f_9 * slg1_552[k]
                    + f_3 * pc_x[k] * slh_768[k];

        t_1021[k] = f_0 * skh_597[k]
                    + f_3 * pc_y[k] * slh_765[k];
    }

#pragma omp simd aligned(t_1022, t_1023, t_1024, t_1025, t_1026, t_1027, pc_x, slg0_554, \
                         slg1_554, slh_770, slh_771, slh_772, slh_773, slh_774, \
                         slh_775 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1022[k] = f_8 * slg0_554[k]
                    - f_9 * slg1_554[k]
                    + f_3 * pc_x[k] * slh_770[k];

        t_1023[k] = f_3 * pc_x[k] * slh_771[k];

        t_1024[k] = f_3 * pc_x[k] * slh_772[k];

        t_1025[k] = f_3 * pc_x[k] * slh_773[k];

        t_1026[k] = f_3 * pc_x[k] * slh_774[k];

        t_1027[k] = f_3 * pc_x[k] * slh_775[k];
    }

#pragma omp simd aligned(t_1028, t_1029, t_1030, t_1031, pc_x, pc_y, pc_z, skh_603, skh_605, \
                         slg0_550, slg0_552, slg1_550, slg1_552, slh_771, slh_773, \
                         slh_776 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1028[k] = f_3 * pc_x[k] * slh_776[k];

        t_1029[k] = f_0 * skh_603[k]
                    + f_1 * slg0_550[k]
                    - f_2 * slg1_550[k]
                    + f_3 * pc_y[k] * slh_771[k];

        t_1030[k] = f_3 * pc_z[k] * slh_771[k];

        t_1031[k] = f_0 * skh_605[k]
                    + f_4 * slg0_552[k]
                    - f_5 * slg1_552[k]
                    + f_3 * pc_y[k] * slh_773[k];
    }

#pragma omp simd aligned(t_1032, t_1033, t_1034, t_1035, pc_y, pc_z, skh_606, skh_607, \
                         skh_608, slg0_553, slg0_554, slg1_553, slg1_554, slh_774, slh_775, \
                         slh_776 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1032[k] = f_0 * skh_606[k]
                    + f_6 * slg0_553[k]
                    - f_7 * slg1_553[k]
                    + f_3 * pc_y[k] * slh_774[k];

        t_1033[k] = f_0 * skh_607[k]
                    + f_8 * slg0_554[k]
                    - f_9 * slg1_554[k]
                    + f_3 * pc_y[k] * slh_775[k];

        t_1034[k] = f_0 * skh_608[k]
                    + f_3 * pc_y[k] * slh_776[k];

        t_1035[k] = f_1 * slg0_554[k]
                    - f_2 * slg1_554[k]
                    + f_3 * pc_z[k] * slh_776[k];
    }

#pragma omp simd aligned(t_1036, t_1037, t_1038, t_1039, pb_z, pc_y, pc_z, ski0_784, ski0_787, \
                         skh_588, skh_609, ski1_784, ski1_787, \
                         slh_777 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1036[k] = pb_z[k] * ski0_784[k]
                    - f_10 * pc_z[k] * ski1_784[k];

        t_1037[k] = f_15 * skh_609[k]
                    + f_3 * pc_y[k] * slh_777[k];

        t_1038[k] = f_11 * skh_588[k]
                    + f_3 * pc_z[k] * slh_777[k];

        t_1039[k] = pb_z[k] * ski0_787[k]
                    - f_10 * pc_z[k] * ski1_787[k];
    }

#pragma omp simd aligned(t_1040, t_1041, t_1042, pb_z, pc_x, pc_y, pc_z, ski0_790, skh_611, \
                         ski1_790, slg0_560, slg1_560, slh_779, \
                         slh_782 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1040[k] = f_15 * skh_611[k]
                    + f_3 * pc_y[k] * slh_779[k];

        t_1041[k] = f_4 * slg0_560[k]
                    - f_5 * slg1_560[k]
                    + f_3 * pc_x[k] * slh_782[k];

        t_1042[k] = pb_z[k] * ski0_790[k]
                    - f_10 * pc_z[k] * ski1_790[k];
    }

#pragma omp simd aligned(t_1043, t_1044, t_1045, pc_x, pc_y, pc_z, skh_591, skh_614, slg0_564, \
                         slg1_564, slh_780, slh_782, slh_786 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1043[k] = f_11 * skh_591[k]
                    + f_3 * pc_z[k] * slh_780[k];

        t_1044[k] = f_15 * skh_614[k]
                    + f_3 * pc_y[k] * slh_782[k];

        t_1045[k] = f_6 * slg0_564[k]
                    - f_7 * slg1_564[k]
                    + f_3 * pc_x[k] * slh_786[k];
    }
}

static auto
compute_prim_sli_three_center_electron_repulsion_0_piece9(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t ski0,
                                                          const size_t skh, const size_t ski1,
                                                          const size_t slg0, const size_t slg1,
                                                          const size_t slh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.5 / gamma;
    const auto f_5 = 1.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 0.5 / gamma;
    const auto f_9 = 0.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 3.5 / q;
    const auto f_16 = 3.0 / q;
    const auto f_17 = 2.5 / q;

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

    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ski0_794 = buffer.data(ski0 + 794);
    const auto *ski0_805 = buffer.data(ski0 + 805);
    const auto *ski0_807 = buffer.data(ski0 + 807);
    const auto *ski0_808 = buffer.data(ski0 + 808);
    const auto *ski0_809 = buffer.data(ski0 + 809);

    const auto *skh_594 = buffer.data(skh + 594);
    const auto *skh_603 = buffer.data(skh + 603);
    const auto *skh_604 = buffer.data(skh + 604);
    const auto *skh_605 = buffer.data(skh + 605);
    const auto *skh_606 = buffer.data(skh + 606);
    const auto *skh_608 = buffer.data(skh + 608);
    const auto *skh_609 = buffer.data(skh + 609);
    const auto *skh_612 = buffer.data(skh + 612);
    const auto *skh_615 = buffer.data(skh + 615);
    const auto *skh_618 = buffer.data(skh + 618);
    const auto *skh_624 = buffer.data(skh + 624);
    const auto *skh_629 = buffer.data(skh + 629);
    const auto *skh_630 = buffer.data(skh + 630);
    const auto *skh_632 = buffer.data(skh + 632);
    const auto *skh_633 = buffer.data(skh + 633);
    const auto *skh_635 = buffer.data(skh + 635);
    const auto *skh_636 = buffer.data(skh + 636);
    const auto *skh_639 = buffer.data(skh + 639);
    const auto *skh_645 = buffer.data(skh + 645);
    const auto *skh_647 = buffer.data(skh + 647);
    const auto *skh_648 = buffer.data(skh + 648);
    const auto *skh_649 = buffer.data(skh + 649);
    const auto *skh_650 = buffer.data(skh + 650);
    const auto *skh_651 = buffer.data(skh + 651);
    const auto *skh_653 = buffer.data(skh + 653);
    const auto *skh_654 = buffer.data(skh + 654);
    const auto *skh_656 = buffer.data(skh + 656);
    const auto *skh_657 = buffer.data(skh + 657);
    const auto *skh_660 = buffer.data(skh + 660);
    const auto *skh_666 = buffer.data(skh + 666);
    const auto *skh_668 = buffer.data(skh + 668);
    const auto *skh_669 = buffer.data(skh + 669);
    const auto *skh_670 = buffer.data(skh + 670);
    const auto *skh_671 = buffer.data(skh + 671);
    const auto *skh_672 = buffer.data(skh + 672);
    const auto *skh_674 = buffer.data(skh + 674);
    const auto *skh_675 = buffer.data(skh + 675);
    const auto *skh_677 = buffer.data(skh + 677);
    const auto *skh_678 = buffer.data(skh + 678);
    const auto *skh_681 = buffer.data(skh + 681);
    const auto *skh_687 = buffer.data(skh + 687);
    const auto *skh_689 = buffer.data(skh + 689);
    const auto *skh_690 = buffer.data(skh + 690);
    const auto *skh_691 = buffer.data(skh + 691);
    const auto *skh_692 = buffer.data(skh + 692);
    const auto *skh_693 = buffer.data(skh + 693);
    const auto *skh_695 = buffer.data(skh + 695);
    const auto *skh_698 = buffer.data(skh + 698);
    const auto *skh_702 = buffer.data(skh + 702);

    const auto *ski1_794 = buffer.data(ski1 + 794);
    const auto *ski1_805 = buffer.data(ski1 + 805);
    const auto *ski1_807 = buffer.data(ski1 + 807);
    const auto *ski1_808 = buffer.data(ski1 + 808);
    const auto *ski1_809 = buffer.data(ski1 + 809);

    const auto *slg0_567 = buffer.data(slg0 + 567);
    const auto *slg0_569 = buffer.data(slg0 + 569);
    const auto *slg0_570 = buffer.data(slg0 + 570);
    const auto *slg0_573 = buffer.data(slg0 + 573);
    const auto *slg0_575 = buffer.data(slg0 + 575);
    const auto *slg0_576 = buffer.data(slg0 + 576);
    const auto *slg0_579 = buffer.data(slg0 + 579);
    const auto *slg0_580 = buffer.data(slg0 + 580);
    const auto *slg0_582 = buffer.data(slg0 + 582);
    const auto *slg0_583 = buffer.data(slg0 + 583);
    const auto *slg0_584 = buffer.data(slg0 + 584);
    const auto *slg0_585 = buffer.data(slg0 + 585);
    const auto *slg0_588 = buffer.data(slg0 + 588);
    const auto *slg0_590 = buffer.data(slg0 + 590);
    const auto *slg0_591 = buffer.data(slg0 + 591);
    const auto *slg0_594 = buffer.data(slg0 + 594);
    const auto *slg0_595 = buffer.data(slg0 + 595);
    const auto *slg0_597 = buffer.data(slg0 + 597);
    const auto *slg0_598 = buffer.data(slg0 + 598);
    const auto *slg0_599 = buffer.data(slg0 + 599);
    const auto *slg0_600 = buffer.data(slg0 + 600);
    const auto *slg0_603 = buffer.data(slg0 + 603);
    const auto *slg0_605 = buffer.data(slg0 + 605);
    const auto *slg0_606 = buffer.data(slg0 + 606);
    const auto *slg0_609 = buffer.data(slg0 + 609);
    const auto *slg0_610 = buffer.data(slg0 + 610);
    const auto *slg0_612 = buffer.data(slg0 + 612);
    const auto *slg0_613 = buffer.data(slg0 + 613);
    const auto *slg0_614 = buffer.data(slg0 + 614);
    const auto *slg0_615 = buffer.data(slg0 + 615);
    const auto *slg0_618 = buffer.data(slg0 + 618);
    const auto *slg0_620 = buffer.data(slg0 + 620);
    const auto *slg0_621 = buffer.data(slg0 + 621);
    const auto *slg0_624 = buffer.data(slg0 + 624);
    const auto *slg0_625 = buffer.data(slg0 + 625);
    const auto *slg0_627 = buffer.data(slg0 + 627);
    const auto *slg0_629 = buffer.data(slg0 + 629);

    const auto *slg1_567 = buffer.data(slg1 + 567);
    const auto *slg1_569 = buffer.data(slg1 + 569);
    const auto *slg1_570 = buffer.data(slg1 + 570);
    const auto *slg1_573 = buffer.data(slg1 + 573);
    const auto *slg1_575 = buffer.data(slg1 + 575);
    const auto *slg1_576 = buffer.data(slg1 + 576);
    const auto *slg1_579 = buffer.data(slg1 + 579);
    const auto *slg1_580 = buffer.data(slg1 + 580);
    const auto *slg1_582 = buffer.data(slg1 + 582);
    const auto *slg1_583 = buffer.data(slg1 + 583);
    const auto *slg1_584 = buffer.data(slg1 + 584);
    const auto *slg1_585 = buffer.data(slg1 + 585);
    const auto *slg1_588 = buffer.data(slg1 + 588);
    const auto *slg1_590 = buffer.data(slg1 + 590);
    const auto *slg1_591 = buffer.data(slg1 + 591);
    const auto *slg1_594 = buffer.data(slg1 + 594);
    const auto *slg1_595 = buffer.data(slg1 + 595);
    const auto *slg1_597 = buffer.data(slg1 + 597);
    const auto *slg1_598 = buffer.data(slg1 + 598);
    const auto *slg1_599 = buffer.data(slg1 + 599);
    const auto *slg1_600 = buffer.data(slg1 + 600);
    const auto *slg1_603 = buffer.data(slg1 + 603);
    const auto *slg1_605 = buffer.data(slg1 + 605);
    const auto *slg1_606 = buffer.data(slg1 + 606);
    const auto *slg1_609 = buffer.data(slg1 + 609);
    const auto *slg1_610 = buffer.data(slg1 + 610);
    const auto *slg1_612 = buffer.data(slg1 + 612);
    const auto *slg1_613 = buffer.data(slg1 + 613);
    const auto *slg1_614 = buffer.data(slg1 + 614);
    const auto *slg1_615 = buffer.data(slg1 + 615);
    const auto *slg1_618 = buffer.data(slg1 + 618);
    const auto *slg1_620 = buffer.data(slg1 + 620);
    const auto *slg1_621 = buffer.data(slg1 + 621);
    const auto *slg1_624 = buffer.data(slg1 + 624);
    const auto *slg1_625 = buffer.data(slg1 + 625);
    const auto *slg1_627 = buffer.data(slg1 + 627);
    const auto *slg1_629 = buffer.data(slg1 + 629);

    const auto *slh_783 = buffer.data(slh + 783);
    const auto *slh_786 = buffer.data(slh + 786);
    const auto *slh_789 = buffer.data(slh + 789);
    const auto *slh_791 = buffer.data(slh + 791);
    const auto *slh_792 = buffer.data(slh + 792);
    const auto *slh_793 = buffer.data(slh + 793);
    const auto *slh_794 = buffer.data(slh + 794);
    const auto *slh_795 = buffer.data(slh + 795);
    const auto *slh_796 = buffer.data(slh + 796);
    const auto *slh_797 = buffer.data(slh + 797);
    const auto *slh_798 = buffer.data(slh + 798);
    const auto *slh_800 = buffer.data(slh + 800);
    const auto *slh_801 = buffer.data(slh + 801);
    const auto *slh_803 = buffer.data(slh + 803);
    const auto *slh_804 = buffer.data(slh + 804);
    const auto *slh_807 = buffer.data(slh + 807);
    const auto *slh_808 = buffer.data(slh + 808);
    const auto *slh_810 = buffer.data(slh + 810);
    const auto *slh_812 = buffer.data(slh + 812);
    const auto *slh_813 = buffer.data(slh + 813);
    const auto *slh_814 = buffer.data(slh + 814);
    const auto *slh_815 = buffer.data(slh + 815);
    const auto *slh_816 = buffer.data(slh + 816);
    const auto *slh_817 = buffer.data(slh + 817);
    const auto *slh_818 = buffer.data(slh + 818);
    const auto *slh_819 = buffer.data(slh + 819);
    const auto *slh_821 = buffer.data(slh + 821);
    const auto *slh_822 = buffer.data(slh + 822);
    const auto *slh_824 = buffer.data(slh + 824);
    const auto *slh_825 = buffer.data(slh + 825);
    const auto *slh_828 = buffer.data(slh + 828);
    const auto *slh_829 = buffer.data(slh + 829);
    const auto *slh_831 = buffer.data(slh + 831);
    const auto *slh_833 = buffer.data(slh + 833);
    const auto *slh_834 = buffer.data(slh + 834);
    const auto *slh_835 = buffer.data(slh + 835);
    const auto *slh_836 = buffer.data(slh + 836);
    const auto *slh_837 = buffer.data(slh + 837);
    const auto *slh_838 = buffer.data(slh + 838);
    const auto *slh_839 = buffer.data(slh + 839);
    const auto *slh_840 = buffer.data(slh + 840);
    const auto *slh_842 = buffer.data(slh + 842);
    const auto *slh_843 = buffer.data(slh + 843);
    const auto *slh_845 = buffer.data(slh + 845);
    const auto *slh_846 = buffer.data(slh + 846);
    const auto *slh_849 = buffer.data(slh + 849);
    const auto *slh_850 = buffer.data(slh + 850);
    const auto *slh_852 = buffer.data(slh + 852);
    const auto *slh_854 = buffer.data(slh + 854);
    const auto *slh_855 = buffer.data(slh + 855);
    const auto *slh_856 = buffer.data(slh + 856);
    const auto *slh_857 = buffer.data(slh + 857);
    const auto *slh_858 = buffer.data(slh + 858);
    const auto *slh_859 = buffer.data(slh + 859);
    const auto *slh_860 = buffer.data(slh + 860);
    const auto *slh_861 = buffer.data(slh + 861);
    const auto *slh_863 = buffer.data(slh + 863);
    const auto *slh_864 = buffer.data(slh + 864);
    const auto *slh_866 = buffer.data(slh + 866);
    const auto *slh_867 = buffer.data(slh + 867);
    const auto *slh_870 = buffer.data(slh + 870);
    const auto *slh_871 = buffer.data(slh + 871);
    const auto *slh_873 = buffer.data(slh + 873);
    const auto *slh_875 = buffer.data(slh + 875);
    const auto *slh_876 = buffer.data(slh + 876);
    const auto *slh_877 = buffer.data(slh + 877);
    const auto *slh_878 = buffer.data(slh + 878);
    const auto *slh_879 = buffer.data(slh + 879);
    const auto *slh_880 = buffer.data(slh + 880);
    const auto *slh_881 = buffer.data(slh + 881);

#pragma omp simd aligned(t_1046, t_1047, t_1048, pb_z, pc_x, pc_z, ski0_794, skh_594, \
                         ski1_794, slg0_567, slg1_567, slh_783, \
                         slh_789 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1046[k] = pb_z[k] * ski0_794[k]
                    - f_10 * pc_z[k] * ski1_794[k];

        t_1047[k] = f_11 * skh_594[k]
                    + f_3 * pc_z[k] * slh_783[k];

        t_1048[k] = f_8 * slg0_567[k]
                    - f_9 * slg1_567[k]
                    + f_3 * pc_x[k] * slh_789[k];
    }

#pragma omp simd aligned(t_1049, t_1050, t_1051, t_1052, t_1053, pc_x, pc_y, skh_618, \
                         slg0_569, slg1_569, slh_786, slh_791, slh_792, slh_793, \
                         slh_794 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1049[k] = f_15 * skh_618[k]
                    + f_3 * pc_y[k] * slh_786[k];

        t_1050[k] = f_8 * slg0_569[k]
                    - f_9 * slg1_569[k]
                    + f_3 * pc_x[k] * slh_791[k];

        t_1051[k] = f_3 * pc_x[k] * slh_792[k];

        t_1052[k] = f_3 * pc_x[k] * slh_793[k];

        t_1053[k] = f_3 * pc_x[k] * slh_794[k];
    }

#pragma omp simd aligned(t_1054, t_1055, t_1056, t_1057, t_1058, pb_z, pc_x, pc_z, ski0_805, \
                         skh_603, ski1_805, slh_792, slh_795, slh_796, \
                         slh_797 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1054[k] = f_3 * pc_x[k] * slh_795[k];

        t_1055[k] = f_3 * pc_x[k] * slh_796[k];

        t_1056[k] = f_3 * pc_x[k] * slh_797[k];

        t_1057[k] = pb_z[k] * ski0_805[k]
                    - f_10 * pc_z[k] * ski1_805[k];

        t_1058[k] = f_11 * skh_603[k]
                    + f_3 * pc_z[k] * slh_792[k];
    }

#pragma omp simd aligned(t_1059, t_1060, t_1061, pb_z, pc_z, ski0_807, ski0_808, ski0_809, \
                         skh_604, skh_605, skh_606, ski1_807, ski1_808, \
                         ski1_809 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1059[k] = pb_z[k] * ski0_807[k]
                    + f_12 * skh_604[k]
                    - f_10 * pc_z[k] * ski1_807[k];

        t_1060[k] = pb_z[k] * ski0_808[k]
                    + f_13 * skh_605[k]
                    - f_10 * pc_z[k] * ski1_808[k];

        t_1061[k] = pb_z[k] * ski0_809[k]
                    + f_14 * skh_606[k]
                    - f_10 * pc_z[k] * ski1_809[k];
    }

#pragma omp simd aligned(t_1062, t_1063, t_1064, t_1065, pc_x, pc_y, pc_z, skh_608, skh_629, \
                         skh_630, slg0_569, slg0_570, slg1_569, slg1_570, slh_797, \
                         slh_798 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1062[k] = f_15 * skh_629[k]
                    + f_3 * pc_y[k] * slh_797[k];

        t_1063[k] = f_11 * skh_608[k]
                    + f_1 * slg0_569[k]
                    - f_2 * slg1_569[k]
                    + f_3 * pc_z[k] * slh_797[k];

        t_1064[k] = f_1 * slg0_570[k]
                    - f_2 * slg1_570[k]
                    + f_3 * pc_x[k] * slh_798[k];

        t_1065[k] = f_16 * skh_630[k]
                    + f_3 * pc_y[k] * slh_798[k];
    }

#pragma omp simd aligned(t_1066, t_1067, t_1068, pc_x, pc_y, pc_z, skh_609, skh_632, slg0_573, \
                         slg1_573, slh_798, slh_800, slh_801 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1066[k] = f_12 * skh_609[k]
                    + f_3 * pc_z[k] * slh_798[k];

        t_1067[k] = f_4 * slg0_573[k]
                    - f_5 * slg1_573[k]
                    + f_3 * pc_x[k] * slh_801[k];

        t_1068[k] = f_16 * skh_632[k]
                    + f_3 * pc_y[k] * slh_800[k];
    }

#pragma omp simd aligned(t_1069, t_1070, t_1071, t_1072, pc_x, pc_y, pc_z, skh_612, skh_635, \
                         slg0_575, slg0_576, slg1_575, slg1_576, slh_801, slh_803, \
                         slh_804 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1069[k] = f_4 * slg0_575[k]
                    - f_5 * slg1_575[k]
                    + f_3 * pc_x[k] * slh_803[k];

        t_1070[k] = f_6 * slg0_576[k]
                    - f_7 * slg1_576[k]
                    + f_3 * pc_x[k] * slh_804[k];

        t_1071[k] = f_12 * skh_612[k]
                    + f_3 * pc_z[k] * slh_801[k];

        t_1072[k] = f_16 * skh_635[k]
                    + f_3 * pc_y[k] * slh_803[k];
    }

#pragma omp simd aligned(t_1073, t_1074, t_1075, pc_x, pc_z, skh_615, slg0_579, slg0_580, \
                         slg1_579, slg1_580, slh_804, slh_807, \
                         slh_808 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1073[k] = f_6 * slg0_579[k]
                    - f_7 * slg1_579[k]
                    + f_3 * pc_x[k] * slh_807[k];

        t_1074[k] = f_8 * slg0_580[k]
                    - f_9 * slg1_580[k]
                    + f_3 * pc_x[k] * slh_808[k];

        t_1075[k] = f_12 * skh_615[k]
                    + f_3 * pc_z[k] * slh_804[k];
    }

#pragma omp simd aligned(t_1076, t_1077, t_1078, t_1079, pc_x, pc_y, skh_639, slg0_582, \
                         slg0_584, slg1_582, slg1_584, slh_807, slh_810, slh_812, \
                         slh_813 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1076[k] = f_8 * slg0_582[k]
                    - f_9 * slg1_582[k]
                    + f_3 * pc_x[k] * slh_810[k];

        t_1077[k] = f_16 * skh_639[k]
                    + f_3 * pc_y[k] * slh_807[k];

        t_1078[k] = f_8 * slg0_584[k]
                    - f_9 * slg1_584[k]
                    + f_3 * pc_x[k] * slh_812[k];

        t_1079[k] = f_3 * pc_x[k] * slh_813[k];
    }

#pragma omp simd aligned(t_1080, t_1081, t_1082, t_1083, t_1084, pc_x, slh_814, slh_815, \
                         slh_816, slh_817, slh_818 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1080[k] = f_3 * pc_x[k] * slh_814[k];

        t_1081[k] = f_3 * pc_x[k] * slh_815[k];

        t_1082[k] = f_3 * pc_x[k] * slh_816[k];

        t_1083[k] = f_3 * pc_x[k] * slh_817[k];

        t_1084[k] = f_3 * pc_x[k] * slh_818[k];
    }

#pragma omp simd aligned(t_1085, t_1086, t_1087, pc_y, pc_z, skh_624, skh_645, skh_647, \
                         slg0_580, slg0_582, slg1_580, slg1_582, slh_813, \
                         slh_815 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1085[k] = f_16 * skh_645[k]
                    + f_1 * slg0_580[k]
                    - f_2 * slg1_580[k]
                    + f_3 * pc_y[k] * slh_813[k];

        t_1086[k] = f_12 * skh_624[k]
                    + f_3 * pc_z[k] * slh_813[k];

        t_1087[k] = f_16 * skh_647[k]
                    + f_4 * slg0_582[k]
                    - f_5 * slg1_582[k]
                    + f_3 * pc_y[k] * slh_815[k];
    }

#pragma omp simd aligned(t_1088, t_1089, t_1090, pc_y, skh_648, skh_649, skh_650, slg0_583, \
                         slg0_584, slg1_583, slg1_584, slh_816, slh_817, \
                         slh_818 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1088[k] = f_16 * skh_648[k]
                    + f_6 * slg0_583[k]
                    - f_7 * slg1_583[k]
                    + f_3 * pc_y[k] * slh_816[k];

        t_1089[k] = f_16 * skh_649[k]
                    + f_8 * slg0_584[k]
                    - f_9 * slg1_584[k]
                    + f_3 * pc_y[k] * slh_817[k];

        t_1090[k] = f_16 * skh_650[k]
                    + f_3 * pc_y[k] * slh_818[k];
    }

#pragma omp simd aligned(t_1091, t_1092, t_1093, t_1094, pc_x, pc_y, pc_z, skh_629, skh_630, \
                         skh_651, slg0_584, slg0_585, slg1_584, slg1_585, slh_818, \
                         slh_819 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1091[k] = f_12 * skh_629[k]
                    + f_1 * slg0_584[k]
                    - f_2 * slg1_584[k]
                    + f_3 * pc_z[k] * slh_818[k];

        t_1092[k] = f_1 * slg0_585[k]
                    - f_2 * slg1_585[k]
                    + f_3 * pc_x[k] * slh_819[k];

        t_1093[k] = f_17 * skh_651[k]
                    + f_3 * pc_y[k] * slh_819[k];

        t_1094[k] = f_13 * skh_630[k]
                    + f_3 * pc_z[k] * slh_819[k];
    }

#pragma omp simd aligned(t_1095, t_1096, t_1097, pc_x, pc_y, skh_653, slg0_588, slg0_590, \
                         slg1_588, slg1_590, slh_821, slh_822, \
                         slh_824 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1095[k] = f_4 * slg0_588[k]
                    - f_5 * slg1_588[k]
                    + f_3 * pc_x[k] * slh_822[k];

        t_1096[k] = f_17 * skh_653[k]
                    + f_3 * pc_y[k] * slh_821[k];

        t_1097[k] = f_4 * slg0_590[k]
                    - f_5 * slg1_590[k]
                    + f_3 * pc_x[k] * slh_824[k];
    }

#pragma omp simd aligned(t_1098, t_1099, t_1100, pc_x, pc_y, pc_z, skh_633, skh_656, slg0_591, \
                         slg1_591, slh_822, slh_824, slh_825 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1098[k] = f_6 * slg0_591[k]
                    - f_7 * slg1_591[k]
                    + f_3 * pc_x[k] * slh_825[k];

        t_1099[k] = f_13 * skh_633[k]
                    + f_3 * pc_z[k] * slh_822[k];

        t_1100[k] = f_17 * skh_656[k]
                    + f_3 * pc_y[k] * slh_824[k];
    }

#pragma omp simd aligned(t_1101, t_1102, t_1103, pc_x, pc_z, skh_636, slg0_594, slg0_595, \
                         slg1_594, slg1_595, slh_825, slh_828, \
                         slh_829 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1101[k] = f_6 * slg0_594[k]
                    - f_7 * slg1_594[k]
                    + f_3 * pc_x[k] * slh_828[k];

        t_1102[k] = f_8 * slg0_595[k]
                    - f_9 * slg1_595[k]
                    + f_3 * pc_x[k] * slh_829[k];

        t_1103[k] = f_13 * skh_636[k]
                    + f_3 * pc_z[k] * slh_825[k];
    }

#pragma omp simd aligned(t_1104, t_1105, t_1106, t_1107, pc_x, pc_y, skh_660, slg0_597, \
                         slg0_599, slg1_597, slg1_599, slh_828, slh_831, slh_833, \
                         slh_834 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1104[k] = f_8 * slg0_597[k]
                    - f_9 * slg1_597[k]
                    + f_3 * pc_x[k] * slh_831[k];

        t_1105[k] = f_17 * skh_660[k]
                    + f_3 * pc_y[k] * slh_828[k];

        t_1106[k] = f_8 * slg0_599[k]
                    - f_9 * slg1_599[k]
                    + f_3 * pc_x[k] * slh_833[k];

        t_1107[k] = f_3 * pc_x[k] * slh_834[k];
    }

#pragma omp simd aligned(t_1108, t_1109, t_1110, t_1111, t_1112, pc_x, slh_835, slh_836, \
                         slh_837, slh_838, slh_839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1108[k] = f_3 * pc_x[k] * slh_835[k];

        t_1109[k] = f_3 * pc_x[k] * slh_836[k];

        t_1110[k] = f_3 * pc_x[k] * slh_837[k];

        t_1111[k] = f_3 * pc_x[k] * slh_838[k];

        t_1112[k] = f_3 * pc_x[k] * slh_839[k];
    }

#pragma omp simd aligned(t_1113, t_1114, t_1115, pc_y, pc_z, skh_645, skh_666, skh_668, \
                         slg0_595, slg0_597, slg1_595, slg1_597, slh_834, \
                         slh_836 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1113[k] = f_17 * skh_666[k]
                    + f_1 * slg0_595[k]
                    - f_2 * slg1_595[k]
                    + f_3 * pc_y[k] * slh_834[k];

        t_1114[k] = f_13 * skh_645[k]
                    + f_3 * pc_z[k] * slh_834[k];

        t_1115[k] = f_17 * skh_668[k]
                    + f_4 * slg0_597[k]
                    - f_5 * slg1_597[k]
                    + f_3 * pc_y[k] * slh_836[k];
    }

#pragma omp simd aligned(t_1116, t_1117, t_1118, pc_y, skh_669, skh_670, skh_671, slg0_598, \
                         slg0_599, slg1_598, slg1_599, slh_837, slh_838, \
                         slh_839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1116[k] = f_17 * skh_669[k]
                    + f_6 * slg0_598[k]
                    - f_7 * slg1_598[k]
                    + f_3 * pc_y[k] * slh_837[k];

        t_1117[k] = f_17 * skh_670[k]
                    + f_8 * slg0_599[k]
                    - f_9 * slg1_599[k]
                    + f_3 * pc_y[k] * slh_838[k];

        t_1118[k] = f_17 * skh_671[k]
                    + f_3 * pc_y[k] * slh_839[k];
    }

#pragma omp simd aligned(t_1119, t_1120, t_1121, t_1122, pc_x, pc_y, pc_z, skh_650, skh_651, \
                         skh_672, slg0_599, slg0_600, slg1_599, slg1_600, slh_839, \
                         slh_840 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1119[k] = f_13 * skh_650[k]
                    + f_1 * slg0_599[k]
                    - f_2 * slg1_599[k]
                    + f_3 * pc_z[k] * slh_839[k];

        t_1120[k] = f_1 * slg0_600[k]
                    - f_2 * slg1_600[k]
                    + f_3 * pc_x[k] * slh_840[k];

        t_1121[k] = f_14 * skh_672[k]
                    + f_3 * pc_y[k] * slh_840[k];

        t_1122[k] = f_14 * skh_651[k]
                    + f_3 * pc_z[k] * slh_840[k];
    }

#pragma omp simd aligned(t_1123, t_1124, t_1125, pc_x, pc_y, skh_674, slg0_603, slg0_605, \
                         slg1_603, slg1_605, slh_842, slh_843, \
                         slh_845 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1123[k] = f_4 * slg0_603[k]
                    - f_5 * slg1_603[k]
                    + f_3 * pc_x[k] * slh_843[k];

        t_1124[k] = f_14 * skh_674[k]
                    + f_3 * pc_y[k] * slh_842[k];

        t_1125[k] = f_4 * slg0_605[k]
                    - f_5 * slg1_605[k]
                    + f_3 * pc_x[k] * slh_845[k];
    }

#pragma omp simd aligned(t_1126, t_1127, t_1128, pc_x, pc_y, pc_z, skh_654, skh_677, slg0_606, \
                         slg1_606, slh_843, slh_845, slh_846 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1126[k] = f_6 * slg0_606[k]
                    - f_7 * slg1_606[k]
                    + f_3 * pc_x[k] * slh_846[k];

        t_1127[k] = f_14 * skh_654[k]
                    + f_3 * pc_z[k] * slh_843[k];

        t_1128[k] = f_14 * skh_677[k]
                    + f_3 * pc_y[k] * slh_845[k];
    }

#pragma omp simd aligned(t_1129, t_1130, t_1131, pc_x, pc_z, skh_657, slg0_609, slg0_610, \
                         slg1_609, slg1_610, slh_846, slh_849, \
                         slh_850 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1129[k] = f_6 * slg0_609[k]
                    - f_7 * slg1_609[k]
                    + f_3 * pc_x[k] * slh_849[k];

        t_1130[k] = f_8 * slg0_610[k]
                    - f_9 * slg1_610[k]
                    + f_3 * pc_x[k] * slh_850[k];

        t_1131[k] = f_14 * skh_657[k]
                    + f_3 * pc_z[k] * slh_846[k];
    }

#pragma omp simd aligned(t_1132, t_1133, t_1134, t_1135, pc_x, pc_y, skh_681, slg0_612, \
                         slg0_614, slg1_612, slg1_614, slh_849, slh_852, slh_854, \
                         slh_855 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1132[k] = f_8 * slg0_612[k]
                    - f_9 * slg1_612[k]
                    + f_3 * pc_x[k] * slh_852[k];

        t_1133[k] = f_14 * skh_681[k]
                    + f_3 * pc_y[k] * slh_849[k];

        t_1134[k] = f_8 * slg0_614[k]
                    - f_9 * slg1_614[k]
                    + f_3 * pc_x[k] * slh_854[k];

        t_1135[k] = f_3 * pc_x[k] * slh_855[k];
    }

#pragma omp simd aligned(t_1136, t_1137, t_1138, t_1139, t_1140, pc_x, slh_856, slh_857, \
                         slh_858, slh_859, slh_860 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1136[k] = f_3 * pc_x[k] * slh_856[k];

        t_1137[k] = f_3 * pc_x[k] * slh_857[k];

        t_1138[k] = f_3 * pc_x[k] * slh_858[k];

        t_1139[k] = f_3 * pc_x[k] * slh_859[k];

        t_1140[k] = f_3 * pc_x[k] * slh_860[k];
    }

#pragma omp simd aligned(t_1141, t_1142, t_1143, pc_y, pc_z, skh_666, skh_687, skh_689, \
                         slg0_610, slg0_612, slg1_610, slg1_612, slh_855, \
                         slh_857 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1141[k] = f_14 * skh_687[k]
                    + f_1 * slg0_610[k]
                    - f_2 * slg1_610[k]
                    + f_3 * pc_y[k] * slh_855[k];

        t_1142[k] = f_14 * skh_666[k]
                    + f_3 * pc_z[k] * slh_855[k];

        t_1143[k] = f_14 * skh_689[k]
                    + f_4 * slg0_612[k]
                    - f_5 * slg1_612[k]
                    + f_3 * pc_y[k] * slh_857[k];
    }

#pragma omp simd aligned(t_1144, t_1145, t_1146, pc_y, skh_690, skh_691, skh_692, slg0_613, \
                         slg0_614, slg1_613, slg1_614, slh_858, slh_859, \
                         slh_860 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1144[k] = f_14 * skh_690[k]
                    + f_6 * slg0_613[k]
                    - f_7 * slg1_613[k]
                    + f_3 * pc_y[k] * slh_858[k];

        t_1145[k] = f_14 * skh_691[k]
                    + f_8 * slg0_614[k]
                    - f_9 * slg1_614[k]
                    + f_3 * pc_y[k] * slh_859[k];

        t_1146[k] = f_14 * skh_692[k]
                    + f_3 * pc_y[k] * slh_860[k];
    }

#pragma omp simd aligned(t_1147, t_1148, t_1149, t_1150, pc_x, pc_y, pc_z, skh_671, skh_672, \
                         skh_693, slg0_614, slg0_615, slg1_614, slg1_615, slh_860, \
                         slh_861 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1147[k] = f_14 * skh_671[k]
                    + f_1 * slg0_614[k]
                    - f_2 * slg1_614[k]
                    + f_3 * pc_z[k] * slh_860[k];

        t_1148[k] = f_1 * slg0_615[k]
                    - f_2 * slg1_615[k]
                    + f_3 * pc_x[k] * slh_861[k];

        t_1149[k] = f_13 * skh_693[k]
                    + f_3 * pc_y[k] * slh_861[k];

        t_1150[k] = f_17 * skh_672[k]
                    + f_3 * pc_z[k] * slh_861[k];
    }

#pragma omp simd aligned(t_1151, t_1152, t_1153, pc_x, pc_y, skh_695, slg0_618, slg0_620, \
                         slg1_618, slg1_620, slh_863, slh_864, \
                         slh_866 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1151[k] = f_4 * slg0_618[k]
                    - f_5 * slg1_618[k]
                    + f_3 * pc_x[k] * slh_864[k];

        t_1152[k] = f_13 * skh_695[k]
                    + f_3 * pc_y[k] * slh_863[k];

        t_1153[k] = f_4 * slg0_620[k]
                    - f_5 * slg1_620[k]
                    + f_3 * pc_x[k] * slh_866[k];
    }

#pragma omp simd aligned(t_1154, t_1155, t_1156, pc_x, pc_y, pc_z, skh_675, skh_698, slg0_621, \
                         slg1_621, slh_864, slh_866, slh_867 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1154[k] = f_6 * slg0_621[k]
                    - f_7 * slg1_621[k]
                    + f_3 * pc_x[k] * slh_867[k];

        t_1155[k] = f_17 * skh_675[k]
                    + f_3 * pc_z[k] * slh_864[k];

        t_1156[k] = f_13 * skh_698[k]
                    + f_3 * pc_y[k] * slh_866[k];
    }

#pragma omp simd aligned(t_1157, t_1158, t_1159, pc_x, pc_z, skh_678, slg0_624, slg0_625, \
                         slg1_624, slg1_625, slh_867, slh_870, \
                         slh_871 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1157[k] = f_6 * slg0_624[k]
                    - f_7 * slg1_624[k]
                    + f_3 * pc_x[k] * slh_870[k];

        t_1158[k] = f_8 * slg0_625[k]
                    - f_9 * slg1_625[k]
                    + f_3 * pc_x[k] * slh_871[k];

        t_1159[k] = f_17 * skh_678[k]
                    + f_3 * pc_z[k] * slh_867[k];
    }

#pragma omp simd aligned(t_1160, t_1161, t_1162, t_1163, pc_x, pc_y, skh_702, slg0_627, \
                         slg0_629, slg1_627, slg1_629, slh_870, slh_873, slh_875, \
                         slh_876 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1160[k] = f_8 * slg0_627[k]
                    - f_9 * slg1_627[k]
                    + f_3 * pc_x[k] * slh_873[k];

        t_1161[k] = f_13 * skh_702[k]
                    + f_3 * pc_y[k] * slh_870[k];

        t_1162[k] = f_8 * slg0_629[k]
                    - f_9 * slg1_629[k]
                    + f_3 * pc_x[k] * slh_875[k];

        t_1163[k] = f_3 * pc_x[k] * slh_876[k];
    }

#pragma omp simd aligned(t_1164, t_1165, t_1166, t_1167, t_1168, pc_x, slh_877, slh_878, \
                         slh_879, slh_880, slh_881 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1164[k] = f_3 * pc_x[k] * slh_877[k];

        t_1165[k] = f_3 * pc_x[k] * slh_878[k];

        t_1166[k] = f_3 * pc_x[k] * slh_879[k];

        t_1167[k] = f_3 * pc_x[k] * slh_880[k];

        t_1168[k] = f_3 * pc_x[k] * slh_881[k];
    }
}

static auto
compute_prim_sli_three_center_electron_repulsion_0_piece10(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t ski0,
                                                           const size_t skh, const size_t ski1,
                                                           const size_t slg0, const size_t slg1,
                                                           const size_t slh, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.5 / gamma;
    const auto f_5 = 1.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 0.5 / gamma;
    const auto f_9 = 0.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 3.5 / q;
    const auto f_16 = 3.0 / q;
    const auto f_17 = 2.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ski0_980 = buffer.data(ski0 + 980);
    const auto *ski0_985 = buffer.data(ski0 + 985);
    const auto *ski0_989 = buffer.data(ski0 + 989);
    const auto *ski0_994 = buffer.data(ski0 + 994);
    const auto *ski0_1001 = buffer.data(ski0 + 1001);
    const auto *ski0_1003 = buffer.data(ski0 + 1003);
    const auto *ski0_1004 = buffer.data(ski0 + 1004);
    const auto *ski0_1005 = buffer.data(ski0 + 1005);
    const auto *ski0_1007 = buffer.data(ski0 + 1007);

    const auto *skh_687 = buffer.data(skh + 687);
    const auto *skh_692 = buffer.data(skh + 692);
    const auto *skh_693 = buffer.data(skh + 693);
    const auto *skh_696 = buffer.data(skh + 696);
    const auto *skh_699 = buffer.data(skh + 699);
    const auto *skh_708 = buffer.data(skh + 708);
    const auto *skh_710 = buffer.data(skh + 710);
    const auto *skh_711 = buffer.data(skh + 711);
    const auto *skh_712 = buffer.data(skh + 712);
    const auto *skh_713 = buffer.data(skh + 713);
    const auto *skh_714 = buffer.data(skh + 714);
    const auto *skh_716 = buffer.data(skh + 716);
    const auto *skh_717 = buffer.data(skh + 717);
    const auto *skh_719 = buffer.data(skh + 719);
    const auto *skh_720 = buffer.data(skh + 720);
    const auto *skh_723 = buffer.data(skh + 723);
    const auto *skh_729 = buffer.data(skh + 729);
    const auto *skh_731 = buffer.data(skh + 731);
    const auto *skh_732 = buffer.data(skh + 732);
    const auto *skh_733 = buffer.data(skh + 733);
    const auto *skh_734 = buffer.data(skh + 734);
    const auto *skh_735 = buffer.data(skh + 735);
    const auto *skh_737 = buffer.data(skh + 737);
    const auto *skh_738 = buffer.data(skh + 738);
    const auto *skh_740 = buffer.data(skh + 740);
    const auto *skh_741 = buffer.data(skh + 741);
    const auto *skh_744 = buffer.data(skh + 744);
    const auto *skh_750 = buffer.data(skh + 750);
    const auto *skh_752 = buffer.data(skh + 752);
    const auto *skh_753 = buffer.data(skh + 753);
    const auto *skh_754 = buffer.data(skh + 754);
    const auto *skh_755 = buffer.data(skh + 755);

    const auto *ski1_980 = buffer.data(ski1 + 980);
    const auto *ski1_985 = buffer.data(ski1 + 985);
    const auto *ski1_989 = buffer.data(ski1 + 989);
    const auto *ski1_994 = buffer.data(ski1 + 994);
    const auto *ski1_1001 = buffer.data(ski1 + 1001);
    const auto *ski1_1003 = buffer.data(ski1 + 1003);
    const auto *ski1_1004 = buffer.data(ski1 + 1004);
    const auto *ski1_1005 = buffer.data(ski1 + 1005);
    const auto *ski1_1007 = buffer.data(ski1 + 1007);

    const auto *slg0_625 = buffer.data(slg0 + 625);
    const auto *slg0_627 = buffer.data(slg0 + 627);
    const auto *slg0_628 = buffer.data(slg0 + 628);
    const auto *slg0_629 = buffer.data(slg0 + 629);
    const auto *slg0_630 = buffer.data(slg0 + 630);
    const auto *slg0_633 = buffer.data(slg0 + 633);
    const auto *slg0_635 = buffer.data(slg0 + 635);
    const auto *slg0_636 = buffer.data(slg0 + 636);
    const auto *slg0_639 = buffer.data(slg0 + 639);
    const auto *slg0_640 = buffer.data(slg0 + 640);
    const auto *slg0_642 = buffer.data(slg0 + 642);
    const auto *slg0_643 = buffer.data(slg0 + 643);
    const auto *slg0_644 = buffer.data(slg0 + 644);
    const auto *slg0_648 = buffer.data(slg0 + 648);
    const auto *slg0_651 = buffer.data(slg0 + 651);
    const auto *slg0_655 = buffer.data(slg0 + 655);
    const auto *slg0_657 = buffer.data(slg0 + 657);
    const auto *slg0_660 = buffer.data(slg0 + 660);
    const auto *slg0_663 = buffer.data(slg0 + 663);
    const auto *slg0_665 = buffer.data(slg0 + 665);
    const auto *slg0_666 = buffer.data(slg0 + 666);
    const auto *slg0_669 = buffer.data(slg0 + 669);
    const auto *slg0_670 = buffer.data(slg0 + 670);
    const auto *slg0_672 = buffer.data(slg0 + 672);
    const auto *slg0_673 = buffer.data(slg0 + 673);
    const auto *slg0_674 = buffer.data(slg0 + 674);

    const auto *slg1_625 = buffer.data(slg1 + 625);
    const auto *slg1_627 = buffer.data(slg1 + 627);
    const auto *slg1_628 = buffer.data(slg1 + 628);
    const auto *slg1_629 = buffer.data(slg1 + 629);
    const auto *slg1_630 = buffer.data(slg1 + 630);
    const auto *slg1_633 = buffer.data(slg1 + 633);
    const auto *slg1_635 = buffer.data(slg1 + 635);
    const auto *slg1_636 = buffer.data(slg1 + 636);
    const auto *slg1_639 = buffer.data(slg1 + 639);
    const auto *slg1_640 = buffer.data(slg1 + 640);
    const auto *slg1_642 = buffer.data(slg1 + 642);
    const auto *slg1_643 = buffer.data(slg1 + 643);
    const auto *slg1_644 = buffer.data(slg1 + 644);
    const auto *slg1_648 = buffer.data(slg1 + 648);
    const auto *slg1_651 = buffer.data(slg1 + 651);
    const auto *slg1_655 = buffer.data(slg1 + 655);
    const auto *slg1_657 = buffer.data(slg1 + 657);
    const auto *slg1_660 = buffer.data(slg1 + 660);
    const auto *slg1_663 = buffer.data(slg1 + 663);
    const auto *slg1_665 = buffer.data(slg1 + 665);
    const auto *slg1_666 = buffer.data(slg1 + 666);
    const auto *slg1_669 = buffer.data(slg1 + 669);
    const auto *slg1_670 = buffer.data(slg1 + 670);
    const auto *slg1_672 = buffer.data(slg1 + 672);
    const auto *slg1_673 = buffer.data(slg1 + 673);
    const auto *slg1_674 = buffer.data(slg1 + 674);

    const auto *slh_876 = buffer.data(slh + 876);
    const auto *slh_878 = buffer.data(slh + 878);
    const auto *slh_879 = buffer.data(slh + 879);
    const auto *slh_880 = buffer.data(slh + 880);
    const auto *slh_881 = buffer.data(slh + 881);
    const auto *slh_882 = buffer.data(slh + 882);
    const auto *slh_884 = buffer.data(slh + 884);
    const auto *slh_885 = buffer.data(slh + 885);
    const auto *slh_887 = buffer.data(slh + 887);
    const auto *slh_888 = buffer.data(slh + 888);
    const auto *slh_891 = buffer.data(slh + 891);
    const auto *slh_892 = buffer.data(slh + 892);
    const auto *slh_894 = buffer.data(slh + 894);
    const auto *slh_896 = buffer.data(slh + 896);
    const auto *slh_897 = buffer.data(slh + 897);
    const auto *slh_898 = buffer.data(slh + 898);
    const auto *slh_899 = buffer.data(slh + 899);
    const auto *slh_900 = buffer.data(slh + 900);
    const auto *slh_901 = buffer.data(slh + 901);
    const auto *slh_902 = buffer.data(slh + 902);
    const auto *slh_903 = buffer.data(slh + 903);
    const auto *slh_905 = buffer.data(slh + 905);
    const auto *slh_906 = buffer.data(slh + 906);
    const auto *slh_908 = buffer.data(slh + 908);
    const auto *slh_909 = buffer.data(slh + 909);
    const auto *slh_912 = buffer.data(slh + 912);
    const auto *slh_913 = buffer.data(slh + 913);
    const auto *slh_915 = buffer.data(slh + 915);
    const auto *slh_918 = buffer.data(slh + 918);
    const auto *slh_919 = buffer.data(slh + 919);
    const auto *slh_920 = buffer.data(slh + 920);
    const auto *slh_921 = buffer.data(slh + 921);
    const auto *slh_922 = buffer.data(slh + 922);
    const auto *slh_923 = buffer.data(slh + 923);
    const auto *slh_924 = buffer.data(slh + 924);
    const auto *slh_926 = buffer.data(slh + 926);
    const auto *slh_927 = buffer.data(slh + 927);
    const auto *slh_929 = buffer.data(slh + 929);
    const auto *slh_930 = buffer.data(slh + 930);
    const auto *slh_933 = buffer.data(slh + 933);
    const auto *slh_934 = buffer.data(slh + 934);
    const auto *slh_936 = buffer.data(slh + 936);
    const auto *slh_938 = buffer.data(slh + 938);
    const auto *slh_939 = buffer.data(slh + 939);
    const auto *slh_940 = buffer.data(slh + 940);
    const auto *slh_941 = buffer.data(slh + 941);
    const auto *slh_942 = buffer.data(slh + 942);
    const auto *slh_943 = buffer.data(slh + 943);
    const auto *slh_944 = buffer.data(slh + 944);

#pragma omp simd aligned(t_1169, t_1170, t_1171, pc_y, pc_z, skh_687, skh_708, skh_710, \
                         slg0_625, slg0_627, slg1_625, slg1_627, slh_876, \
                         slh_878 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1169[k] = f_13 * skh_708[k]
                    + f_1 * slg0_625[k]
                    - f_2 * slg1_625[k]
                    + f_3 * pc_y[k] * slh_876[k];

        t_1170[k] = f_17 * skh_687[k]
                    + f_3 * pc_z[k] * slh_876[k];

        t_1171[k] = f_13 * skh_710[k]
                    + f_4 * slg0_627[k]
                    - f_5 * slg1_627[k]
                    + f_3 * pc_y[k] * slh_878[k];
    }

#pragma omp simd aligned(t_1172, t_1173, t_1174, pc_y, skh_711, skh_712, skh_713, slg0_628, \
                         slg0_629, slg1_628, slg1_629, slh_879, slh_880, \
                         slh_881 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1172[k] = f_13 * skh_711[k]
                    + f_6 * slg0_628[k]
                    - f_7 * slg1_628[k]
                    + f_3 * pc_y[k] * slh_879[k];

        t_1173[k] = f_13 * skh_712[k]
                    + f_8 * slg0_629[k]
                    - f_9 * slg1_629[k]
                    + f_3 * pc_y[k] * slh_880[k];

        t_1174[k] = f_13 * skh_713[k]
                    + f_3 * pc_y[k] * slh_881[k];
    }

#pragma omp simd aligned(t_1175, t_1176, t_1177, t_1178, pc_x, pc_y, pc_z, skh_692, skh_693, \
                         skh_714, slg0_629, slg0_630, slg1_629, slg1_630, slh_881, \
                         slh_882 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1175[k] = f_17 * skh_692[k]
                    + f_1 * slg0_629[k]
                    - f_2 * slg1_629[k]
                    + f_3 * pc_z[k] * slh_881[k];

        t_1176[k] = f_1 * slg0_630[k]
                    - f_2 * slg1_630[k]
                    + f_3 * pc_x[k] * slh_882[k];

        t_1177[k] = f_12 * skh_714[k]
                    + f_3 * pc_y[k] * slh_882[k];

        t_1178[k] = f_16 * skh_693[k]
                    + f_3 * pc_z[k] * slh_882[k];
    }

#pragma omp simd aligned(t_1179, t_1180, t_1181, pc_x, pc_y, skh_716, slg0_633, slg0_635, \
                         slg1_633, slg1_635, slh_884, slh_885, \
                         slh_887 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1179[k] = f_4 * slg0_633[k]
                    - f_5 * slg1_633[k]
                    + f_3 * pc_x[k] * slh_885[k];

        t_1180[k] = f_12 * skh_716[k]
                    + f_3 * pc_y[k] * slh_884[k];

        t_1181[k] = f_4 * slg0_635[k]
                    - f_5 * slg1_635[k]
                    + f_3 * pc_x[k] * slh_887[k];
    }

#pragma omp simd aligned(t_1182, t_1183, t_1184, pc_x, pc_y, pc_z, skh_696, skh_719, slg0_636, \
                         slg1_636, slh_885, slh_887, slh_888 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1182[k] = f_6 * slg0_636[k]
                    - f_7 * slg1_636[k]
                    + f_3 * pc_x[k] * slh_888[k];

        t_1183[k] = f_16 * skh_696[k]
                    + f_3 * pc_z[k] * slh_885[k];

        t_1184[k] = f_12 * skh_719[k]
                    + f_3 * pc_y[k] * slh_887[k];
    }

#pragma omp simd aligned(t_1185, t_1186, t_1187, pc_x, pc_z, skh_699, slg0_639, slg0_640, \
                         slg1_639, slg1_640, slh_888, slh_891, \
                         slh_892 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1185[k] = f_6 * slg0_639[k]
                    - f_7 * slg1_639[k]
                    + f_3 * pc_x[k] * slh_891[k];

        t_1186[k] = f_8 * slg0_640[k]
                    - f_9 * slg1_640[k]
                    + f_3 * pc_x[k] * slh_892[k];

        t_1187[k] = f_16 * skh_699[k]
                    + f_3 * pc_z[k] * slh_888[k];
    }

#pragma omp simd aligned(t_1188, t_1189, t_1190, t_1191, pc_x, pc_y, skh_723, slg0_642, \
                         slg0_644, slg1_642, slg1_644, slh_891, slh_894, slh_896, \
                         slh_897 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1188[k] = f_8 * slg0_642[k]
                    - f_9 * slg1_642[k]
                    + f_3 * pc_x[k] * slh_894[k];

        t_1189[k] = f_12 * skh_723[k]
                    + f_3 * pc_y[k] * slh_891[k];

        t_1190[k] = f_8 * slg0_644[k]
                    - f_9 * slg1_644[k]
                    + f_3 * pc_x[k] * slh_896[k];

        t_1191[k] = f_3 * pc_x[k] * slh_897[k];
    }

#pragma omp simd aligned(t_1192, t_1193, t_1194, t_1195, t_1196, pc_x, slh_898, slh_899, \
                         slh_900, slh_901, slh_902 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1192[k] = f_3 * pc_x[k] * slh_898[k];

        t_1193[k] = f_3 * pc_x[k] * slh_899[k];

        t_1194[k] = f_3 * pc_x[k] * slh_900[k];

        t_1195[k] = f_3 * pc_x[k] * slh_901[k];

        t_1196[k] = f_3 * pc_x[k] * slh_902[k];
    }

#pragma omp simd aligned(t_1197, t_1198, t_1199, pc_y, pc_z, skh_708, skh_729, skh_731, \
                         slg0_640, slg0_642, slg1_640, slg1_642, slh_897, \
                         slh_899 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1197[k] = f_12 * skh_729[k]
                    + f_1 * slg0_640[k]
                    - f_2 * slg1_640[k]
                    + f_3 * pc_y[k] * slh_897[k];

        t_1198[k] = f_16 * skh_708[k]
                    + f_3 * pc_z[k] * slh_897[k];

        t_1199[k] = f_12 * skh_731[k]
                    + f_4 * slg0_642[k]
                    - f_5 * slg1_642[k]
                    + f_3 * pc_y[k] * slh_899[k];
    }

#pragma omp simd aligned(t_1200, t_1201, t_1202, pc_y, skh_732, skh_733, skh_734, slg0_643, \
                         slg0_644, slg1_643, slg1_644, slh_900, slh_901, \
                         slh_902 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1200[k] = f_12 * skh_732[k]
                    + f_6 * slg0_643[k]
                    - f_7 * slg1_643[k]
                    + f_3 * pc_y[k] * slh_900[k];

        t_1201[k] = f_12 * skh_733[k]
                    + f_8 * slg0_644[k]
                    - f_9 * slg1_644[k]
                    + f_3 * pc_y[k] * slh_901[k];

        t_1202[k] = f_12 * skh_734[k]
                    + f_3 * pc_y[k] * slh_902[k];
    }

#pragma omp simd aligned(t_1203, t_1204, t_1205, t_1206, pb_y, pc_y, pc_z, ski0_980, skh_713, \
                         skh_714, skh_735, ski1_980, slg0_644, slg1_644, slh_902, \
                         slh_903 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1203[k] = f_16 * skh_713[k]
                    + f_1 * slg0_644[k]
                    - f_2 * slg1_644[k]
                    + f_3 * pc_z[k] * slh_902[k];

        t_1204[k] = pb_y[k] * ski0_980[k]
                    - f_10 * pc_y[k] * ski1_980[k];

        t_1205[k] = f_11 * skh_735[k]
                    + f_3 * pc_y[k] * slh_903[k];

        t_1206[k] = f_15 * skh_714[k]
                    + f_3 * pc_z[k] * slh_903[k];
    }

#pragma omp simd aligned(t_1207, t_1208, t_1209, pb_y, pc_x, pc_y, ski0_985, skh_737, \
                         ski1_985, slg0_648, slg1_648, slh_905, \
                         slh_906 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1207[k] = f_4 * slg0_648[k]
                    - f_5 * slg1_648[k]
                    + f_3 * pc_x[k] * slh_906[k];

        t_1208[k] = f_11 * skh_737[k]
                    + f_3 * pc_y[k] * slh_905[k];

        t_1209[k] = pb_y[k] * ski0_985[k]
                    - f_10 * pc_y[k] * ski1_985[k];
    }

#pragma omp simd aligned(t_1210, t_1211, t_1212, pc_x, pc_y, pc_z, skh_717, skh_740, slg0_651, \
                         slg1_651, slh_906, slh_908, slh_909 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1210[k] = f_6 * slg0_651[k]
                    - f_7 * slg1_651[k]
                    + f_3 * pc_x[k] * slh_909[k];

        t_1211[k] = f_15 * skh_717[k]
                    + f_3 * pc_z[k] * slh_906[k];

        t_1212[k] = f_11 * skh_740[k]
                    + f_3 * pc_y[k] * slh_908[k];
    }

#pragma omp simd aligned(t_1213, t_1214, t_1215, pb_y, pc_x, pc_y, pc_z, ski0_989, skh_720, \
                         ski1_989, slg0_655, slg1_655, slh_909, \
                         slh_913 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1213[k] = pb_y[k] * ski0_989[k]
                    - f_10 * pc_y[k] * ski1_989[k];

        t_1214[k] = f_8 * slg0_655[k]
                    - f_9 * slg1_655[k]
                    + f_3 * pc_x[k] * slh_913[k];

        t_1215[k] = f_15 * skh_720[k]
                    + f_3 * pc_z[k] * slh_909[k];
    }

#pragma omp simd aligned(t_1216, t_1217, t_1218, t_1219, pb_y, pc_x, pc_y, ski0_994, skh_744, \
                         ski1_994, slg0_657, slg1_657, slh_912, slh_915, \
                         slh_918 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1216[k] = f_8 * slg0_657[k]
                    - f_9 * slg1_657[k]
                    + f_3 * pc_x[k] * slh_915[k];

        t_1217[k] = f_11 * skh_744[k]
                    + f_3 * pc_y[k] * slh_912[k];

        t_1218[k] = pb_y[k] * ski0_994[k]
                    - f_10 * pc_y[k] * ski1_994[k];

        t_1219[k] = f_3 * pc_x[k] * slh_918[k];
    }

#pragma omp simd aligned(t_1220, t_1221, t_1222, t_1223, t_1224, pc_x, slh_919, slh_920, \
                         slh_921, slh_922, slh_923 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1220[k] = f_3 * pc_x[k] * slh_919[k];

        t_1221[k] = f_3 * pc_x[k] * slh_920[k];

        t_1222[k] = f_3 * pc_x[k] * slh_921[k];

        t_1223[k] = f_3 * pc_x[k] * slh_922[k];

        t_1224[k] = f_3 * pc_x[k] * slh_923[k];
    }

#pragma omp simd aligned(t_1225, t_1226, t_1227, pb_y, pc_y, pc_z, ski0_1001, ski0_1003, \
                         skh_729, skh_750, skh_752, ski1_1001, ski1_1003, \
                         slh_918 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1225[k] = pb_y[k] * ski0_1001[k]
                    + f_16 * skh_750[k]
                    - f_10 * pc_y[k] * ski1_1001[k];

        t_1226[k] = f_15 * skh_729[k]
                    + f_3 * pc_z[k] * slh_918[k];

        t_1227[k] = pb_y[k] * ski0_1003[k]
                    + f_14 * skh_752[k]
                    - f_10 * pc_y[k] * ski1_1003[k];
    }

#pragma omp simd aligned(t_1228, t_1229, t_1230, t_1231, pb_y, pc_y, ski0_1004, ski0_1005, \
                         ski0_1007, skh_753, skh_754, skh_755, ski1_1004, ski1_1005, \
                         ski1_1007, slh_923 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1228[k] = pb_y[k] * ski0_1004[k]
                    + f_13 * skh_753[k]
                    - f_10 * pc_y[k] * ski1_1004[k];

        t_1229[k] = pb_y[k] * ski0_1005[k]
                    + f_12 * skh_754[k]
                    - f_10 * pc_y[k] * ski1_1005[k];

        t_1230[k] = f_11 * skh_755[k]
                    + f_3 * pc_y[k] * slh_923[k];

        t_1231[k] = pb_y[k] * ski0_1007[k]
                    - f_10 * pc_y[k] * ski1_1007[k];
    }

#pragma omp simd aligned(t_1232, t_1233, t_1234, t_1235, t_1236, pc_x, pc_y, pc_z, skh_735, \
                         slg0_660, slg0_663, slg1_660, slg1_663, slh_924, slh_926, \
                         slh_927 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1232[k] = f_1 * slg0_660[k]
                    - f_2 * slg1_660[k]
                    + f_3 * pc_x[k] * slh_924[k];

        t_1233[k] = f_3 * pc_y[k] * slh_924[k];

        t_1234[k] = f_0 * skh_735[k]
                    + f_3 * pc_z[k] * slh_924[k];

        t_1235[k] = f_4 * slg0_663[k]
                    - f_5 * slg1_663[k]
                    + f_3 * pc_x[k] * slh_927[k];

        t_1236[k] = f_3 * pc_y[k] * slh_926[k];
    }

#pragma omp simd aligned(t_1237, t_1238, t_1239, t_1240, pc_x, pc_y, pc_z, skh_738, slg0_665, \
                         slg0_666, slg1_665, slg1_666, slh_927, slh_929, \
                         slh_930 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1237[k] = f_4 * slg0_665[k]
                    - f_5 * slg1_665[k]
                    + f_3 * pc_x[k] * slh_929[k];

        t_1238[k] = f_6 * slg0_666[k]
                    - f_7 * slg1_666[k]
                    + f_3 * pc_x[k] * slh_930[k];

        t_1239[k] = f_0 * skh_738[k]
                    + f_3 * pc_z[k] * slh_927[k];

        t_1240[k] = f_3 * pc_y[k] * slh_929[k];
    }

#pragma omp simd aligned(t_1241, t_1242, t_1243, pc_x, pc_z, skh_741, slg0_669, slg0_670, \
                         slg1_669, slg1_670, slh_930, slh_933, \
                         slh_934 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1241[k] = f_6 * slg0_669[k]
                    - f_7 * slg1_669[k]
                    + f_3 * pc_x[k] * slh_933[k];

        t_1242[k] = f_8 * slg0_670[k]
                    - f_9 * slg1_670[k]
                    + f_3 * pc_x[k] * slh_934[k];

        t_1243[k] = f_0 * skh_741[k]
                    + f_3 * pc_z[k] * slh_930[k];
    }

#pragma omp simd aligned(t_1244, t_1245, t_1246, t_1247, t_1248, pc_x, pc_y, slg0_672, \
                         slg0_674, slg1_672, slg1_674, slh_933, slh_936, slh_938, slh_939, \
                         slh_940 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1244[k] = f_8 * slg0_672[k]
                    - f_9 * slg1_672[k]
                    + f_3 * pc_x[k] * slh_936[k];

        t_1245[k] = f_3 * pc_y[k] * slh_933[k];

        t_1246[k] = f_8 * slg0_674[k]
                    - f_9 * slg1_674[k]
                    + f_3 * pc_x[k] * slh_938[k];

        t_1247[k] = f_3 * pc_x[k] * slh_939[k];

        t_1248[k] = f_3 * pc_x[k] * slh_940[k];
    }

#pragma omp simd aligned(t_1249, t_1250, t_1251, t_1252, t_1253, pc_x, pc_y, slg0_670, \
                         slg1_670, slh_939, slh_941, slh_942, slh_943, \
                         slh_944 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1249[k] = f_3 * pc_x[k] * slh_941[k];

        t_1250[k] = f_3 * pc_x[k] * slh_942[k];

        t_1251[k] = f_3 * pc_x[k] * slh_943[k];

        t_1252[k] = f_3 * pc_x[k] * slh_944[k];

        t_1253[k] = f_1 * slg0_670[k]
                    - f_2 * slg1_670[k]
                    + f_3 * pc_y[k] * slh_939[k];
    }

#pragma omp simd aligned(t_1254, t_1255, t_1256, pc_y, pc_z, skh_750, slg0_672, slg0_673, \
                         slg1_672, slg1_673, slh_939, slh_941, \
                         slh_942 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1254[k] = f_0 * skh_750[k]
                    + f_3 * pc_z[k] * slh_939[k];

        t_1255[k] = f_4 * slg0_672[k]
                    - f_5 * slg1_672[k]
                    + f_3 * pc_y[k] * slh_941[k];

        t_1256[k] = f_6 * slg0_673[k]
                    - f_7 * slg1_673[k]
                    + f_3 * pc_y[k] * slh_942[k];
    }

#pragma omp simd aligned(t_1257, t_1258, t_1259, pc_y, pc_z, skh_755, slg0_674, slg1_674, \
                         slh_943, slh_944 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1257[k] = f_8 * slg0_674[k]
                    - f_9 * slg1_674[k]
                    + f_3 * pc_y[k] * slh_943[k];

        t_1258[k] = f_3 * pc_y[k] * slh_944[k];

        t_1259[k] = f_0 * skh_755[k]
                    + f_1 * slg0_674[k]
                    - f_2 * slg1_674[k]
                    + f_3 * pc_z[k] * slh_944[k];
    }
}

auto
compute_prim_sli_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t ski0, const size_t skh,
                                                   const size_t ski1, const size_t slg0,
                                                   const size_t slg1, const size_t slh,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_sli_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, ski0, skh,
                                                              ski1, slg0, slg1, slh, ncols,
                                                              gamma, p, q);

    compute_prim_sli_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, ski0, skh,
                                                              ski1, slg0, slg1, slh, ncols,
                                                              gamma, p, q);

    compute_prim_sli_three_center_electron_repulsion_0_piece2(buffer, target, pb, pc, ski0, skh,
                                                              ski1, slg0, slg1, slh, ncols,
                                                              gamma, p, q);

    compute_prim_sli_three_center_electron_repulsion_0_piece3(buffer, target, pb, pc, ski0, skh,
                                                              ski1, slg0, slg1, slh, ncols,
                                                              gamma, p, q);

    compute_prim_sli_three_center_electron_repulsion_0_piece4(buffer, target, pb, pc, ski0, skh,
                                                              ski1, slg0, slg1, slh, ncols,
                                                              gamma, p, q);

    compute_prim_sli_three_center_electron_repulsion_0_piece5(buffer, target, pb, pc, ski0, skh,
                                                              ski1, slg0, slg1, slh, ncols,
                                                              gamma, p, q);

    compute_prim_sli_three_center_electron_repulsion_0_piece6(buffer, target, pb, pc, ski0, skh,
                                                              ski1, slg0, slg1, slh, ncols,
                                                              gamma, p, q);

    compute_prim_sli_three_center_electron_repulsion_0_piece7(buffer, target, pb, pc, ski0, skh,
                                                              ski1, slh, ncols, gamma, p, q);

    compute_prim_sli_three_center_electron_repulsion_0_piece8(buffer, target, pb, pc, ski0, skh,
                                                              ski1, slg0, slg1, slh, ncols,
                                                              gamma, p, q);

    compute_prim_sli_three_center_electron_repulsion_0_piece9(buffer, target, pb, pc, ski0, skh,
                                                              ski1, slg0, slg1, slh, ncols,
                                                              gamma, p, q);

    compute_prim_sli_three_center_electron_repulsion_0_piece10(buffer, target, pb, pc, ski0,
                                                               skh, ski1, slg0, slg1, slh,
                                                               ncols, gamma, p, q);
}

}  // namespace simdt3ceri
