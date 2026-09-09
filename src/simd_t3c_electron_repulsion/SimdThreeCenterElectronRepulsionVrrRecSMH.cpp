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


#include "SimdThreeCenterElectronRepulsionVrrRecSMH.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_smh_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t slh0,
                                                          const size_t slg, const size_t slh1,
                                                          const size_t smf0, const size_t smf1,
                                                          const size_t smg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / q;
    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.0 / gamma;
    const auto f_5 = p / (gamma * q);
    const auto f_6 = 0.5 / gamma;
    const auto f_7 = 0.5 * p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_12 = 4.0 / q;
    const auto f_13 = 3.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *slh0_0 = buffer.data(slh0 + 0);
    const auto *slh0_3 = buffer.data(slh0 + 3);
    const auto *slh0_5 = buffer.data(slh0 + 5);
    const auto *slh0_6 = buffer.data(slh0 + 6);
    const auto *slh0_9 = buffer.data(slh0 + 9);
    const auto *slh0_15 = buffer.data(slh0 + 15);
    const auto *slh0_20 = buffer.data(slh0 + 20);
    const auto *slh0_24 = buffer.data(slh0 + 24);
    const auto *slh0_27 = buffer.data(slh0 + 27);
    const auto *slh0_36 = buffer.data(slh0 + 36);
    const auto *slh0_42 = buffer.data(slh0 + 42);
    const auto *slh0_47 = buffer.data(slh0 + 47);
    const auto *slh0_51 = buffer.data(slh0 + 51);
    const auto *slh0_62 = buffer.data(slh0 + 62);

    const auto *slg_0 = buffer.data(slg + 0);
    const auto *slg_1 = buffer.data(slg + 1);
    const auto *slg_2 = buffer.data(slg + 2);
    const auto *slg_3 = buffer.data(slg + 3);
    const auto *slg_5 = buffer.data(slg + 5);
    const auto *slg_6 = buffer.data(slg + 6);
    const auto *slg_9 = buffer.data(slg + 9);
    const auto *slg_10 = buffer.data(slg + 10);
    const auto *slg_11 = buffer.data(slg + 11);
    const auto *slg_12 = buffer.data(slg + 12);
    const auto *slg_13 = buffer.data(slg + 13);
    const auto *slg_14 = buffer.data(slg + 14);
    const auto *slg_15 = buffer.data(slg + 15);
    const auto *slg_17 = buffer.data(slg + 17);
    const auto *slg_18 = buffer.data(slg + 18);
    const auto *slg_20 = buffer.data(slg + 20);
    const auto *slg_25 = buffer.data(slg + 25);
    const auto *slg_26 = buffer.data(slg + 26);
    const auto *slg_27 = buffer.data(slg + 27);
    const auto *slg_28 = buffer.data(slg + 28);
    const auto *slg_29 = buffer.data(slg + 29);
    const auto *slg_30 = buffer.data(slg + 30);
    const auto *slg_32 = buffer.data(slg + 32);
    const auto *slg_33 = buffer.data(slg + 33);
    const auto *slg_35 = buffer.data(slg + 35);
    const auto *slg_40 = buffer.data(slg + 40);
    const auto *slg_41 = buffer.data(slg + 41);
    const auto *slg_42 = buffer.data(slg + 42);
    const auto *slg_43 = buffer.data(slg + 43);
    const auto *slg_44 = buffer.data(slg + 44);
    const auto *slg_45 = buffer.data(slg + 45);
    const auto *slg_48 = buffer.data(slg + 48);
    const auto *slg_50 = buffer.data(slg + 50);
    const auto *slg_51 = buffer.data(slg + 51);
    const auto *slg_54 = buffer.data(slg + 54);
    const auto *slg_55 = buffer.data(slg + 55);
    const auto *slg_56 = buffer.data(slg + 56);
    const auto *slg_57 = buffer.data(slg + 57);
    const auto *slg_58 = buffer.data(slg + 58);
    const auto *slg_59 = buffer.data(slg + 59);
    const auto *slg_70 = buffer.data(slg + 70);
    const auto *slg_71 = buffer.data(slg + 71);
    const auto *slg_72 = buffer.data(slg + 72);
    const auto *slg_73 = buffer.data(slg + 73);
    const auto *slg_74 = buffer.data(slg + 74);
    const auto *slg_75 = buffer.data(slg + 75);
    const auto *slg_78 = buffer.data(slg + 78);
    const auto *slg_80 = buffer.data(slg + 80);
    const auto *slg_81 = buffer.data(slg + 81);
    const auto *slg_84 = buffer.data(slg + 84);
    const auto *slg_85 = buffer.data(slg + 85);
    const auto *slg_86 = buffer.data(slg + 86);
    const auto *slg_87 = buffer.data(slg + 87);
    const auto *slg_88 = buffer.data(slg + 88);
    const auto *slg_89 = buffer.data(slg + 89);

    const auto *slh1_0 = buffer.data(slh1 + 0);
    const auto *slh1_3 = buffer.data(slh1 + 3);
    const auto *slh1_5 = buffer.data(slh1 + 5);
    const auto *slh1_6 = buffer.data(slh1 + 6);
    const auto *slh1_9 = buffer.data(slh1 + 9);
    const auto *slh1_15 = buffer.data(slh1 + 15);
    const auto *slh1_20 = buffer.data(slh1 + 20);
    const auto *slh1_24 = buffer.data(slh1 + 24);
    const auto *slh1_27 = buffer.data(slh1 + 27);
    const auto *slh1_36 = buffer.data(slh1 + 36);
    const auto *slh1_42 = buffer.data(slh1 + 42);
    const auto *slh1_47 = buffer.data(slh1 + 47);
    const auto *slh1_51 = buffer.data(slh1 + 51);
    const auto *slh1_62 = buffer.data(slh1 + 62);

    const auto *smf0_0 = buffer.data(smf0 + 0);
    const auto *smf0_3 = buffer.data(smf0 + 3);
    const auto *smf0_5 = buffer.data(smf0 + 5);
    const auto *smf0_6 = buffer.data(smf0 + 6);
    const auto *smf0_8 = buffer.data(smf0 + 8);
    const auto *smf0_9 = buffer.data(smf0 + 9);
    const auto *smf0_16 = buffer.data(smf0 + 16);
    const auto *smf0_18 = buffer.data(smf0 + 18);
    const auto *smf0_19 = buffer.data(smf0 + 19);
    const auto *smf0_28 = buffer.data(smf0 + 28);
    const auto *smf0_29 = buffer.data(smf0 + 29);
    const auto *smf0_30 = buffer.data(smf0 + 30);
    const auto *smf0_33 = buffer.data(smf0 + 33);
    const auto *smf0_35 = buffer.data(smf0 + 35);
    const auto *smf0_36 = buffer.data(smf0 + 36);
    const auto *smf0_38 = buffer.data(smf0 + 38);
    const auto *smf0_39 = buffer.data(smf0 + 39);
    const auto *smf0_48 = buffer.data(smf0 + 48);
    const auto *smf0_49 = buffer.data(smf0 + 49);
    const auto *smf0_50 = buffer.data(smf0 + 50);
    const auto *smf0_53 = buffer.data(smf0 + 53);
    const auto *smf0_55 = buffer.data(smf0 + 55);
    const auto *smf0_56 = buffer.data(smf0 + 56);
    const auto *smf0_58 = buffer.data(smf0 + 58);
    const auto *smf0_59 = buffer.data(smf0 + 59);

    const auto *smf1_0 = buffer.data(smf1 + 0);
    const auto *smf1_3 = buffer.data(smf1 + 3);
    const auto *smf1_5 = buffer.data(smf1 + 5);
    const auto *smf1_6 = buffer.data(smf1 + 6);
    const auto *smf1_8 = buffer.data(smf1 + 8);
    const auto *smf1_9 = buffer.data(smf1 + 9);
    const auto *smf1_16 = buffer.data(smf1 + 16);
    const auto *smf1_18 = buffer.data(smf1 + 18);
    const auto *smf1_19 = buffer.data(smf1 + 19);
    const auto *smf1_28 = buffer.data(smf1 + 28);
    const auto *smf1_29 = buffer.data(smf1 + 29);
    const auto *smf1_30 = buffer.data(smf1 + 30);
    const auto *smf1_33 = buffer.data(smf1 + 33);
    const auto *smf1_35 = buffer.data(smf1 + 35);
    const auto *smf1_36 = buffer.data(smf1 + 36);
    const auto *smf1_38 = buffer.data(smf1 + 38);
    const auto *smf1_39 = buffer.data(smf1 + 39);
    const auto *smf1_48 = buffer.data(smf1 + 48);
    const auto *smf1_49 = buffer.data(smf1 + 49);
    const auto *smf1_50 = buffer.data(smf1 + 50);
    const auto *smf1_53 = buffer.data(smf1 + 53);
    const auto *smf1_55 = buffer.data(smf1 + 55);
    const auto *smf1_56 = buffer.data(smf1 + 56);
    const auto *smf1_58 = buffer.data(smf1 + 58);
    const auto *smf1_59 = buffer.data(smf1 + 59);

    const auto *smg_0 = buffer.data(smg + 0);
    const auto *smg_2 = buffer.data(smg + 2);
    const auto *smg_3 = buffer.data(smg + 3);
    const auto *smg_5 = buffer.data(smg + 5);
    const auto *smg_6 = buffer.data(smg + 6);
    const auto *smg_9 = buffer.data(smg + 9);
    const auto *smg_10 = buffer.data(smg + 10);
    const auto *smg_11 = buffer.data(smg + 11);
    const auto *smg_12 = buffer.data(smg + 12);
    const auto *smg_13 = buffer.data(smg + 13);
    const auto *smg_14 = buffer.data(smg + 14);
    const auto *smg_15 = buffer.data(smg + 15);
    const auto *smg_17 = buffer.data(smg + 17);
    const auto *smg_18 = buffer.data(smg + 18);
    const auto *smg_20 = buffer.data(smg + 20);
    const auto *smg_25 = buffer.data(smg + 25);
    const auto *smg_26 = buffer.data(smg + 26);
    const auto *smg_27 = buffer.data(smg + 27);
    const auto *smg_28 = buffer.data(smg + 28);
    const auto *smg_29 = buffer.data(smg + 29);
    const auto *smg_30 = buffer.data(smg + 30);
    const auto *smg_32 = buffer.data(smg + 32);
    const auto *smg_33 = buffer.data(smg + 33);
    const auto *smg_35 = buffer.data(smg + 35);
    const auto *smg_40 = buffer.data(smg + 40);
    const auto *smg_41 = buffer.data(smg + 41);
    const auto *smg_42 = buffer.data(smg + 42);
    const auto *smg_43 = buffer.data(smg + 43);
    const auto *smg_44 = buffer.data(smg + 44);
    const auto *smg_45 = buffer.data(smg + 45);
    const auto *smg_47 = buffer.data(smg + 47);
    const auto *smg_48 = buffer.data(smg + 48);
    const auto *smg_50 = buffer.data(smg + 50);
    const auto *smg_51 = buffer.data(smg + 51);
    const auto *smg_54 = buffer.data(smg + 54);
    const auto *smg_55 = buffer.data(smg + 55);
    const auto *smg_56 = buffer.data(smg + 56);
    const auto *smg_57 = buffer.data(smg + 57);
    const auto *smg_58 = buffer.data(smg + 58);
    const auto *smg_59 = buffer.data(smg + 59);
    const auto *smg_60 = buffer.data(smg + 60);
    const auto *smg_62 = buffer.data(smg + 62);
    const auto *smg_63 = buffer.data(smg + 63);
    const auto *smg_65 = buffer.data(smg + 65);
    const auto *smg_70 = buffer.data(smg + 70);
    const auto *smg_71 = buffer.data(smg + 71);
    const auto *smg_72 = buffer.data(smg + 72);
    const auto *smg_73 = buffer.data(smg + 73);
    const auto *smg_74 = buffer.data(smg + 74);
    const auto *smg_75 = buffer.data(smg + 75);
    const auto *smg_77 = buffer.data(smg + 77);
    const auto *smg_78 = buffer.data(smg + 78);
    const auto *smg_80 = buffer.data(smg + 80);
    const auto *smg_81 = buffer.data(smg + 81);
    const auto *smg_84 = buffer.data(smg + 84);
    const auto *smg_85 = buffer.data(smg + 85);
    const auto *smg_86 = buffer.data(smg + 86);
    const auto *smg_87 = buffer.data(smg + 87);
    const auto *smg_88 = buffer.data(smg + 88);
    const auto *smg_89 = buffer.data(smg + 89);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, slg_0, slg_3, smf0_0, smf0_3, \
                         smf1_0, smf1_3, smg_0, smg_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * slg_0[k]
                 + f_1 * smf0_0[k]
                 - f_2 * smf1_0[k]
                 + f_3 * pc_x[k] * smg_0[k];

        t_1[k] = f_3 * pc_y[k] * smg_0[k];

        t_2[k] = f_3 * pc_z[k] * smg_0[k];

        t_3[k] = f_0 * slg_3[k]
                 + f_4 * smf0_3[k]
                 - f_5 * smf1_3[k]
                 + f_3 * pc_x[k] * smg_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pc_x, pc_y, slg_5, slg_6, smf0_5, smf0_6, smf1_5, \
                         smf1_6, smg_2, smg_5, smg_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * pc_y[k] * smg_2[k];

        t_5[k] = f_0 * slg_5[k]
                 + f_4 * smf0_5[k]
                 - f_5 * smf1_5[k]
                 + f_3 * pc_x[k] * smg_5[k];

        t_6[k] = f_0 * slg_6[k]
                 + f_6 * smf0_6[k]
                 - f_7 * smf1_6[k]
                 + f_3 * pc_x[k] * smg_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, pc_x, pc_y, pc_z, slg_9, slg_10, smf0_9, smf1_9, \
                         smg_3, smg_5, smg_9, smg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_3 * pc_z[k] * smg_3[k];

        t_8[k] = f_3 * pc_y[k] * smg_5[k];

        t_9[k] = f_0 * slg_9[k]
                 + f_6 * smf0_9[k]
                 - f_7 * smf1_9[k]
                 + f_3 * pc_x[k] * smg_9[k];

        t_10[k] = f_0 * slg_10[k]
                  + f_3 * pc_x[k] * smg_10[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, pc_x, slg_11, slg_12, slg_13, slg_14, smg_11, \
                         smg_12, smg_13, smg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_0 * slg_11[k]
                  + f_3 * pc_x[k] * smg_11[k];

        t_12[k] = f_0 * slg_12[k]
                  + f_3 * pc_x[k] * smg_12[k];

        t_13[k] = f_0 * slg_13[k]
                  + f_3 * pc_x[k] * smg_13[k];

        t_14[k] = f_0 * slg_14[k]
                  + f_3 * pc_x[k] * smg_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pc_y, pc_z, smf0_6, smf0_8, smf0_9, smf1_6, \
                         smf1_8, smf1_9, smg_10, smg_12, smg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_1 * smf0_6[k]
                  - f_2 * smf1_6[k]
                  + f_3 * pc_y[k] * smg_10[k];

        t_16[k] = f_3 * pc_z[k] * smg_10[k];

        t_17[k] = f_4 * smf0_8[k]
                  - f_5 * smf1_8[k]
                  + f_3 * pc_y[k] * smg_12[k];

        t_18[k] = f_6 * smf0_9[k]
                  - f_7 * smf1_9[k]
                  + f_3 * pc_y[k] * smg_13[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pb_y, pc_y, pc_z, slh0_0, slg_0, \
                         slh1_0, smf0_9, smf1_9, smg_14, smg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * pc_y[k] * smg_14[k];

        t_20[k] = f_1 * smf0_9[k]
                  - f_2 * smf1_9[k]
                  + f_3 * pc_z[k] * smg_14[k];

        t_21[k] = pb_y[k] * slh0_0[k]
                  - f_8 * pc_y[k] * slh1_0[k];

        t_22[k] = f_9 * slg_0[k]
                  + f_3 * pc_y[k] * smg_15[k];

        t_23[k] = f_3 * pc_z[k] * smg_15[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pb_y, pc_y, slh0_3, slh0_5, slh0_6, slg_1, \
                         slg_2, slg_3, slh1_3, slh1_5, slh1_6, smg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pb_y[k] * slh0_3[k]
                  + f_10 * slg_1[k]
                  - f_8 * pc_y[k] * slh1_3[k];

        t_25[k] = f_9 * slg_2[k]
                  + f_3 * pc_y[k] * smg_17[k];

        t_26[k] = pb_y[k] * slh0_5[k]
                  - f_8 * pc_y[k] * slh1_5[k];

        t_27[k] = pb_y[k] * slh0_6[k]
                  + f_11 * slg_3[k]
                  - f_8 * pc_y[k] * slh1_6[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pb_y, pc_x, pc_y, pc_z, slh0_9, slg_5, \
                         slg_25, slh1_9, smg_18, smg_20, smg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_3 * pc_z[k] * smg_18[k];

        t_29[k] = f_9 * slg_5[k]
                  + f_3 * pc_y[k] * smg_20[k];

        t_30[k] = pb_y[k] * slh0_9[k]
                  - f_8 * pc_y[k] * slh1_9[k];

        t_31[k] = f_12 * slg_25[k]
                  + f_3 * pc_x[k] * smg_25[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pc_x, slg_26, slg_27, slg_28, slg_29, smg_26, \
                         smg_27, smg_28, smg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_12 * slg_26[k]
                  + f_3 * pc_x[k] * smg_26[k];

        t_33[k] = f_12 * slg_27[k]
                  + f_3 * pc_x[k] * smg_27[k];

        t_34[k] = f_12 * slg_28[k]
                  + f_3 * pc_x[k] * smg_28[k];

        t_35[k] = f_12 * slg_29[k]
                  + f_3 * pc_x[k] * smg_29[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pc_y, pc_z, slg_10, slg_12, smf0_16, smf0_18, \
                         smf1_16, smf1_18, smg_25, smg_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_9 * slg_10[k]
                  + f_1 * smf0_16[k]
                  - f_2 * smf1_16[k]
                  + f_3 * pc_y[k] * smg_25[k];

        t_37[k] = f_3 * pc_z[k] * smg_25[k];

        t_38[k] = f_9 * slg_12[k]
                  + f_4 * smf0_18[k]
                  - f_5 * smf1_18[k]
                  + f_3 * pc_y[k] * smg_27[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pb_y, pc_y, slh0_20, slg_13, slg_14, slh1_20, \
                         smf0_19, smf1_19, smg_28, smg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_9 * slg_13[k]
                  + f_6 * smf0_19[k]
                  - f_7 * smf1_19[k]
                  + f_3 * pc_y[k] * smg_28[k];

        t_40[k] = f_9 * slg_14[k]
                  + f_3 * pc_y[k] * smg_29[k];

        t_41[k] = pb_y[k] * slh0_20[k]
                  - f_8 * pc_y[k] * slh1_20[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, pb_z, pc_y, pc_z, slh0_0, slh0_3, \
                         slg_0, slh1_0, slh1_3, smg_30, smg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_z[k] * slh0_0[k]
                  - f_8 * pc_z[k] * slh1_0[k];

        t_43[k] = f_3 * pc_y[k] * smg_30[k];

        t_44[k] = f_9 * slg_0[k]
                  + f_3 * pc_z[k] * smg_30[k];

        t_45[k] = pb_z[k] * slh0_3[k]
                  - f_8 * pc_z[k] * slh1_3[k];

        t_46[k] = f_3 * pc_y[k] * smg_32[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pb_z, pc_y, pc_z, slh0_5, slh0_6, slg_2, \
                         slg_3, slh1_5, slh1_6, smg_33, smg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = pb_z[k] * slh0_5[k]
                  + f_10 * slg_2[k]
                  - f_8 * pc_z[k] * slh1_5[k];

        t_48[k] = pb_z[k] * slh0_6[k]
                  - f_8 * pc_z[k] * slh1_6[k];

        t_49[k] = f_9 * slg_3[k]
                  + f_3 * pc_z[k] * smg_33[k];

        t_50[k] = f_3 * pc_y[k] * smg_35[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pb_z, pc_x, pc_z, slh0_9, slg_5, slg_40, \
                         slg_41, slg_42, slh1_9, smg_40, smg_41, \
                         smg_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = pb_z[k] * slh0_9[k]
                  + f_11 * slg_5[k]
                  - f_8 * pc_z[k] * slh1_9[k];

        t_52[k] = f_12 * slg_40[k]
                  + f_3 * pc_x[k] * smg_40[k];

        t_53[k] = f_12 * slg_41[k]
                  + f_3 * pc_x[k] * smg_41[k];

        t_54[k] = f_12 * slg_42[k]
                  + f_3 * pc_x[k] * smg_42[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pb_z, pc_x, pc_z, slh0_15, slg_10, slg_43, \
                         slg_44, slh1_15, smg_40, smg_43, smg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_12 * slg_43[k]
                  + f_3 * pc_x[k] * smg_43[k];

        t_56[k] = f_12 * slg_44[k]
                  + f_3 * pc_x[k] * smg_44[k];

        t_57[k] = pb_z[k] * slh0_15[k]
                  - f_8 * pc_z[k] * slh1_15[k];

        t_58[k] = f_9 * slg_10[k]
                  + f_3 * pc_z[k] * smg_40[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pc_y, pc_z, slg_14, smf0_28, smf0_29, \
                         smf1_28, smf1_29, smg_42, smg_43, smg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_4 * smf0_28[k]
                  - f_5 * smf1_28[k]
                  + f_3 * pc_y[k] * smg_42[k];

        t_60[k] = f_6 * smf0_29[k]
                  - f_7 * smf1_29[k]
                  + f_3 * pc_y[k] * smg_43[k];

        t_61[k] = f_3 * pc_y[k] * smg_44[k];

        t_62[k] = f_9 * slg_14[k]
                  + f_1 * smf0_29[k]
                  - f_2 * smf1_29[k]
                  + f_3 * pc_z[k] * smg_44[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pc_x, pc_y, pc_z, slg_15, slg_45, slg_48, \
                         smf0_30, smf0_33, smf1_30, smf1_33, smg_45, \
                         smg_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_13 * slg_45[k]
                  + f_1 * smf0_30[k]
                  - f_2 * smf1_30[k]
                  + f_3 * pc_x[k] * smg_45[k];

        t_64[k] = f_10 * slg_15[k]
                  + f_3 * pc_y[k] * smg_45[k];

        t_65[k] = f_3 * pc_z[k] * smg_45[k];

        t_66[k] = f_13 * slg_48[k]
                  + f_4 * smf0_33[k]
                  - f_5 * smf1_33[k]
                  + f_3 * pc_x[k] * smg_48[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, pc_x, pc_y, slg_17, slg_50, slg_51, smf0_35, \
                         smf0_36, smf1_35, smf1_36, smg_47, smg_50, \
                         smg_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_10 * slg_17[k]
                  + f_3 * pc_y[k] * smg_47[k];

        t_68[k] = f_13 * slg_50[k]
                  + f_4 * smf0_35[k]
                  - f_5 * smf1_35[k]
                  + f_3 * pc_x[k] * smg_50[k];

        t_69[k] = f_13 * slg_51[k]
                  + f_6 * smf0_36[k]
                  - f_7 * smf1_36[k]
                  + f_3 * pc_x[k] * smg_51[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pc_x, pc_y, pc_z, slg_20, slg_54, slg_55, \
                         smf0_39, smf1_39, smg_48, smg_50, smg_54, \
                         smg_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_3 * pc_z[k] * smg_48[k];

        t_71[k] = f_10 * slg_20[k]
                  + f_3 * pc_y[k] * smg_50[k];

        t_72[k] = f_13 * slg_54[k]
                  + f_6 * smf0_39[k]
                  - f_7 * smf1_39[k]
                  + f_3 * pc_x[k] * smg_54[k];

        t_73[k] = f_13 * slg_55[k]
                  + f_3 * pc_x[k] * smg_55[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pc_x, slg_56, slg_57, slg_58, slg_59, smg_56, \
                         smg_57, smg_58, smg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_13 * slg_56[k]
                  + f_3 * pc_x[k] * smg_56[k];

        t_75[k] = f_13 * slg_57[k]
                  + f_3 * pc_x[k] * smg_57[k];

        t_76[k] = f_13 * slg_58[k]
                  + f_3 * pc_x[k] * smg_58[k];

        t_77[k] = f_13 * slg_59[k]
                  + f_3 * pc_x[k] * smg_59[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pc_y, pc_z, slg_25, slg_27, smf0_36, smf0_38, \
                         smf1_36, smf1_38, smg_55, smg_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_10 * slg_25[k]
                  + f_1 * smf0_36[k]
                  - f_2 * smf1_36[k]
                  + f_3 * pc_y[k] * smg_55[k];

        t_79[k] = f_3 * pc_z[k] * smg_55[k];

        t_80[k] = f_10 * slg_27[k]
                  + f_4 * smf0_38[k]
                  - f_5 * smf1_38[k]
                  + f_3 * pc_y[k] * smg_57[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pb_y, pc_y, pc_z, slh0_42, slg_28, slg_29, \
                         slh1_42, smf0_39, smf1_39, smg_58, smg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_10 * slg_28[k]
                  + f_6 * smf0_39[k]
                  - f_7 * smf1_39[k]
                  + f_3 * pc_y[k] * smg_58[k];

        t_82[k] = f_10 * slg_29[k]
                  + f_3 * pc_y[k] * smg_59[k];

        t_83[k] = f_1 * smf0_39[k]
                  - f_2 * smf1_39[k]
                  + f_3 * pc_z[k] * smg_59[k];

        t_84[k] = pb_y[k] * slh0_42[k]
                  - f_8 * pc_y[k] * slh1_42[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pb_z, pc_y, pc_z, slh0_24, slg_15, slg_30, \
                         slg_32, slh1_24, smg_60, smg_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_9 * slg_30[k]
                  + f_3 * pc_y[k] * smg_60[k];

        t_86[k] = f_9 * slg_15[k]
                  + f_3 * pc_z[k] * smg_60[k];

        t_87[k] = pb_z[k] * slh0_24[k]
                  - f_8 * pc_z[k] * slh1_24[k];

        t_88[k] = f_9 * slg_32[k]
                  + f_3 * pc_y[k] * smg_62[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, pb_y, pb_z, pc_y, pc_z, slh0_27, slh0_47, \
                         slg_18, slg_35, slh1_27, slh1_47, smg_63, \
                         smg_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = pb_y[k] * slh0_47[k]
                  - f_8 * pc_y[k] * slh1_47[k];

        t_90[k] = pb_z[k] * slh0_27[k]
                  - f_8 * pc_z[k] * slh1_27[k];

        t_91[k] = f_9 * slg_18[k]
                  + f_3 * pc_z[k] * smg_63[k];

        t_92[k] = f_9 * slg_35[k]
                  + f_3 * pc_y[k] * smg_65[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, pb_y, pc_x, pc_y, slh0_51, slg_70, slg_71, \
                         slg_72, slh1_51, smg_70, smg_71, smg_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = pb_y[k] * slh0_51[k]
                  - f_8 * pc_y[k] * slh1_51[k];

        t_94[k] = f_13 * slg_70[k]
                  + f_3 * pc_x[k] * smg_70[k];

        t_95[k] = f_13 * slg_71[k]
                  + f_3 * pc_x[k] * smg_71[k];

        t_96[k] = f_13 * slg_72[k]
                  + f_3 * pc_x[k] * smg_72[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pb_z, pc_x, pc_z, slh0_36, slg_25, slg_73, \
                         slg_74, slh1_36, smg_70, smg_73, smg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_13 * slg_73[k]
                  + f_3 * pc_x[k] * smg_73[k];

        t_98[k] = f_13 * slg_74[k]
                  + f_3 * pc_x[k] * smg_74[k];

        t_99[k] = pb_z[k] * slh0_36[k]
                  - f_8 * pc_z[k] * slh1_36[k];

        t_100[k] = f_9 * slg_25[k]
                   + f_3 * pc_z[k] * smg_70[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, pc_y, slg_42, slg_43, slg_44, smf0_48, smf0_49, \
                         smf1_48, smf1_49, smg_72, smg_73, smg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_9 * slg_42[k]
                   + f_4 * smf0_48[k]
                   - f_5 * smf1_48[k]
                   + f_3 * pc_y[k] * smg_72[k];

        t_102[k] = f_9 * slg_43[k]
                   + f_6 * smf0_49[k]
                   - f_7 * smf1_49[k]
                   + f_3 * pc_y[k] * smg_73[k];

        t_103[k] = f_9 * slg_44[k]
                   + f_3 * pc_y[k] * smg_74[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pb_y, pc_x, pc_y, pc_z, slh0_62, slg_30, \
                         slg_75, slh1_62, smf0_50, smf1_50, smg_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pb_y[k] * slh0_62[k]
                   - f_8 * pc_y[k] * slh1_62[k];

        t_105[k] = f_13 * slg_75[k]
                   + f_1 * smf0_50[k]
                   - f_2 * smf1_50[k]
                   + f_3 * pc_x[k] * smg_75[k];

        t_106[k] = f_3 * pc_y[k] * smg_75[k];

        t_107[k] = f_10 * slg_30[k]
                   + f_3 * pc_z[k] * smg_75[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, pc_x, pc_y, slg_78, slg_80, smf0_53, smf0_55, \
                         smf1_53, smf1_55, smg_77, smg_78, smg_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_13 * slg_78[k]
                   + f_4 * smf0_53[k]
                   - f_5 * smf1_53[k]
                   + f_3 * pc_x[k] * smg_78[k];

        t_109[k] = f_3 * pc_y[k] * smg_77[k];

        t_110[k] = f_13 * slg_80[k]
                   + f_4 * smf0_55[k]
                   - f_5 * smf1_55[k]
                   + f_3 * pc_x[k] * smg_80[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, pc_x, pc_y, pc_z, slg_33, slg_81, smf0_56, \
                         smf1_56, smg_78, smg_80, smg_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_13 * slg_81[k]
                   + f_6 * smf0_56[k]
                   - f_7 * smf1_56[k]
                   + f_3 * pc_x[k] * smg_81[k];

        t_112[k] = f_10 * slg_33[k]
                   + f_3 * pc_z[k] * smg_78[k];

        t_113[k] = f_3 * pc_y[k] * smg_80[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pc_x, slg_84, slg_85, slg_86, slg_87, \
                         smf0_59, smf1_59, smg_84, smg_85, smg_86, \
                         smg_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_13 * slg_84[k]
                   + f_6 * smf0_59[k]
                   - f_7 * smf1_59[k]
                   + f_3 * pc_x[k] * smg_84[k];

        t_115[k] = f_13 * slg_85[k]
                   + f_3 * pc_x[k] * smg_85[k];

        t_116[k] = f_13 * slg_86[k]
                   + f_3 * pc_x[k] * smg_86[k];

        t_117[k] = f_13 * slg_87[k]
                   + f_3 * pc_x[k] * smg_87[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, pc_x, pc_y, pc_z, slg_40, slg_88, slg_89, \
                         smf0_56, smf1_56, smg_85, smg_88, smg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_13 * slg_88[k]
                   + f_3 * pc_x[k] * smg_88[k];

        t_119[k] = f_13 * slg_89[k]
                   + f_3 * pc_x[k] * smg_89[k];

        t_120[k] = f_1 * smf0_56[k]
                   - f_2 * smf1_56[k]
                   + f_3 * pc_y[k] * smg_85[k];

        t_121[k] = f_10 * slg_40[k]
                   + f_3 * pc_z[k] * smg_85[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, pc_y, pc_z, slg_44, smf0_58, smf0_59, \
                         smf1_58, smf1_59, smg_87, smg_88, smg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_4 * smf0_58[k]
                   - f_5 * smf1_58[k]
                   + f_3 * pc_y[k] * smg_87[k];

        t_123[k] = f_6 * smf0_59[k]
                   - f_7 * smf1_59[k]
                   + f_3 * pc_y[k] * smg_88[k];

        t_124[k] = f_3 * pc_y[k] * smg_89[k];

        t_125[k] = f_10 * slg_44[k]
                   + f_1 * smf0_59[k]
                   - f_2 * smf1_59[k]
                   + f_3 * pc_z[k] * smg_89[k];
    }
}

static auto
compute_prim_smh_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t slh0,
                                                          const size_t slg, const size_t slh1,
                                                          const size_t smf0, const size_t smf1,
                                                          const size_t smg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.0 / gamma;
    const auto f_5 = p / (gamma * q);
    const auto f_6 = 0.5 / gamma;
    const auto f_7 = 0.5 * p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_14 = 3.0 / q;
    const auto f_15 = 2.5 / q;
    const auto f_16 = 2.0 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *slh0_63 = buffer.data(slh0 + 63);
    const auto *slh0_66 = buffer.data(slh0 + 66);
    const auto *slh0_69 = buffer.data(slh0 + 69);
    const auto *slh0_78 = buffer.data(slh0 + 78);
    const auto *slh0_105 = buffer.data(slh0 + 105);
    const auto *slh0_108 = buffer.data(slh0 + 108);
    const auto *slh0_110 = buffer.data(slh0 + 110);
    const auto *slh0_111 = buffer.data(slh0 + 111);
    const auto *slh0_114 = buffer.data(slh0 + 114);
    const auto *slh0_125 = buffer.data(slh0 + 125);
    const auto *slh0_126 = buffer.data(slh0 + 126);
    const auto *slh0_129 = buffer.data(slh0 + 129);
    const auto *slh0_132 = buffer.data(slh0 + 132);

    const auto *slg_45 = buffer.data(slg + 45);
    const auto *slg_47 = buffer.data(slg + 47);
    const auto *slg_48 = buffer.data(slg + 48);
    const auto *slg_50 = buffer.data(slg + 50);
    const auto *slg_55 = buffer.data(slg + 55);
    const auto *slg_57 = buffer.data(slg + 57);
    const auto *slg_58 = buffer.data(slg + 58);
    const auto *slg_59 = buffer.data(slg + 59);
    const auto *slg_60 = buffer.data(slg + 60);
    const auto *slg_62 = buffer.data(slg + 62);
    const auto *slg_63 = buffer.data(slg + 63);
    const auto *slg_65 = buffer.data(slg + 65);
    const auto *slg_70 = buffer.data(slg + 70);
    const auto *slg_72 = buffer.data(slg + 72);
    const auto *slg_73 = buffer.data(slg + 73);
    const auto *slg_74 = buffer.data(slg + 74);
    const auto *slg_75 = buffer.data(slg + 75);
    const auto *slg_76 = buffer.data(slg + 76);
    const auto *slg_77 = buffer.data(slg + 77);
    const auto *slg_78 = buffer.data(slg + 78);
    const auto *slg_80 = buffer.data(slg + 80);
    const auto *slg_85 = buffer.data(slg + 85);
    const auto *slg_87 = buffer.data(slg + 87);
    const auto *slg_88 = buffer.data(slg + 88);
    const auto *slg_89 = buffer.data(slg + 89);
    const auto *slg_90 = buffer.data(slg + 90);
    const auto *slg_92 = buffer.data(slg + 92);
    const auto *slg_93 = buffer.data(slg + 93);
    const auto *slg_95 = buffer.data(slg + 95);
    const auto *slg_96 = buffer.data(slg + 96);
    const auto *slg_99 = buffer.data(slg + 99);
    const auto *slg_100 = buffer.data(slg + 100);
    const auto *slg_101 = buffer.data(slg + 101);
    const auto *slg_102 = buffer.data(slg + 102);
    const auto *slg_103 = buffer.data(slg + 103);
    const auto *slg_104 = buffer.data(slg + 104);
    const auto *slg_105 = buffer.data(slg + 105);
    const auto *slg_107 = buffer.data(slg + 107);
    const auto *slg_110 = buffer.data(slg + 110);
    const auto *slg_114 = buffer.data(slg + 114);
    const auto *slg_115 = buffer.data(slg + 115);
    const auto *slg_116 = buffer.data(slg + 116);
    const auto *slg_117 = buffer.data(slg + 117);
    const auto *slg_118 = buffer.data(slg + 118);
    const auto *slg_119 = buffer.data(slg + 119);
    const auto *slg_130 = buffer.data(slg + 130);
    const auto *slg_131 = buffer.data(slg + 131);
    const auto *slg_132 = buffer.data(slg + 132);
    const auto *slg_133 = buffer.data(slg + 133);
    const auto *slg_134 = buffer.data(slg + 134);
    const auto *slg_135 = buffer.data(slg + 135);
    const auto *slg_138 = buffer.data(slg + 138);
    const auto *slg_140 = buffer.data(slg + 140);
    const auto *slg_141 = buffer.data(slg + 141);
    const auto *slg_144 = buffer.data(slg + 144);
    const auto *slg_145 = buffer.data(slg + 145);
    const auto *slg_146 = buffer.data(slg + 146);
    const auto *slg_147 = buffer.data(slg + 147);
    const auto *slg_148 = buffer.data(slg + 148);
    const auto *slg_149 = buffer.data(slg + 149);
    const auto *slg_150 = buffer.data(slg + 150);
    const auto *slg_153 = buffer.data(slg + 153);
    const auto *slg_155 = buffer.data(slg + 155);
    const auto *slg_156 = buffer.data(slg + 156);
    const auto *slg_159 = buffer.data(slg + 159);
    const auto *slg_160 = buffer.data(slg + 160);
    const auto *slg_161 = buffer.data(slg + 161);
    const auto *slg_162 = buffer.data(slg + 162);
    const auto *slg_163 = buffer.data(slg + 163);
    const auto *slg_164 = buffer.data(slg + 164);
    const auto *slg_170 = buffer.data(slg + 170);
    const auto *slg_174 = buffer.data(slg + 174);
    const auto *slg_175 = buffer.data(slg + 175);
    const auto *slg_176 = buffer.data(slg + 176);

    const auto *slh1_63 = buffer.data(slh1 + 63);
    const auto *slh1_66 = buffer.data(slh1 + 66);
    const auto *slh1_69 = buffer.data(slh1 + 69);
    const auto *slh1_78 = buffer.data(slh1 + 78);
    const auto *slh1_105 = buffer.data(slh1 + 105);
    const auto *slh1_108 = buffer.data(slh1 + 108);
    const auto *slh1_110 = buffer.data(slh1 + 110);
    const auto *slh1_111 = buffer.data(slh1 + 111);
    const auto *slh1_114 = buffer.data(slh1 + 114);
    const auto *slh1_125 = buffer.data(slh1 + 125);
    const auto *slh1_126 = buffer.data(slh1 + 126);
    const auto *slh1_129 = buffer.data(slh1 + 129);
    const auto *slh1_132 = buffer.data(slh1 + 132);

    const auto *smf0_60 = buffer.data(smf0 + 60);
    const auto *smf0_63 = buffer.data(smf0 + 63);
    const auto *smf0_65 = buffer.data(smf0 + 65);
    const auto *smf0_66 = buffer.data(smf0 + 66);
    const auto *smf0_68 = buffer.data(smf0 + 68);
    const auto *smf0_69 = buffer.data(smf0 + 69);
    const auto *smf0_75 = buffer.data(smf0 + 75);
    const auto *smf0_78 = buffer.data(smf0 + 78);
    const auto *smf0_79 = buffer.data(smf0 + 79);
    const auto *smf0_86 = buffer.data(smf0 + 86);
    const auto *smf0_88 = buffer.data(smf0 + 88);
    const auto *smf0_89 = buffer.data(smf0 + 89);
    const auto *smf0_90 = buffer.data(smf0 + 90);
    const auto *smf0_93 = buffer.data(smf0 + 93);
    const auto *smf0_95 = buffer.data(smf0 + 95);
    const auto *smf0_96 = buffer.data(smf0 + 96);
    const auto *smf0_98 = buffer.data(smf0 + 98);
    const auto *smf0_99 = buffer.data(smf0 + 99);
    const auto *smf0_100 = buffer.data(smf0 + 100);
    const auto *smf0_103 = buffer.data(smf0 + 103);
    const auto *smf0_105 = buffer.data(smf0 + 105);
    const auto *smf0_106 = buffer.data(smf0 + 106);
    const auto *smf0_108 = buffer.data(smf0 + 108);
    const auto *smf0_109 = buffer.data(smf0 + 109);
    const auto *smf0_115 = buffer.data(smf0 + 115);
    const auto *smf0_119 = buffer.data(smf0 + 119);

    const auto *smf1_60 = buffer.data(smf1 + 60);
    const auto *smf1_63 = buffer.data(smf1 + 63);
    const auto *smf1_65 = buffer.data(smf1 + 65);
    const auto *smf1_66 = buffer.data(smf1 + 66);
    const auto *smf1_68 = buffer.data(smf1 + 68);
    const auto *smf1_69 = buffer.data(smf1 + 69);
    const auto *smf1_75 = buffer.data(smf1 + 75);
    const auto *smf1_78 = buffer.data(smf1 + 78);
    const auto *smf1_79 = buffer.data(smf1 + 79);
    const auto *smf1_86 = buffer.data(smf1 + 86);
    const auto *smf1_88 = buffer.data(smf1 + 88);
    const auto *smf1_89 = buffer.data(smf1 + 89);
    const auto *smf1_90 = buffer.data(smf1 + 90);
    const auto *smf1_93 = buffer.data(smf1 + 93);
    const auto *smf1_95 = buffer.data(smf1 + 95);
    const auto *smf1_96 = buffer.data(smf1 + 96);
    const auto *smf1_98 = buffer.data(smf1 + 98);
    const auto *smf1_99 = buffer.data(smf1 + 99);
    const auto *smf1_100 = buffer.data(smf1 + 100);
    const auto *smf1_103 = buffer.data(smf1 + 103);
    const auto *smf1_105 = buffer.data(smf1 + 105);
    const auto *smf1_106 = buffer.data(smf1 + 106);
    const auto *smf1_108 = buffer.data(smf1 + 108);
    const auto *smf1_109 = buffer.data(smf1 + 109);
    const auto *smf1_115 = buffer.data(smf1 + 115);
    const auto *smf1_119 = buffer.data(smf1 + 119);

    const auto *smg_90 = buffer.data(smg + 90);
    const auto *smg_92 = buffer.data(smg + 92);
    const auto *smg_93 = buffer.data(smg + 93);
    const auto *smg_95 = buffer.data(smg + 95);
    const auto *smg_96 = buffer.data(smg + 96);
    const auto *smg_99 = buffer.data(smg + 99);
    const auto *smg_100 = buffer.data(smg + 100);
    const auto *smg_101 = buffer.data(smg + 101);
    const auto *smg_102 = buffer.data(smg + 102);
    const auto *smg_103 = buffer.data(smg + 103);
    const auto *smg_104 = buffer.data(smg + 104);
    const auto *smg_105 = buffer.data(smg + 105);
    const auto *smg_107 = buffer.data(smg + 107);
    const auto *smg_108 = buffer.data(smg + 108);
    const auto *smg_110 = buffer.data(smg + 110);
    const auto *smg_114 = buffer.data(smg + 114);
    const auto *smg_115 = buffer.data(smg + 115);
    const auto *smg_116 = buffer.data(smg + 116);
    const auto *smg_117 = buffer.data(smg + 117);
    const auto *smg_118 = buffer.data(smg + 118);
    const auto *smg_119 = buffer.data(smg + 119);
    const auto *smg_120 = buffer.data(smg + 120);
    const auto *smg_122 = buffer.data(smg + 122);
    const auto *smg_123 = buffer.data(smg + 123);
    const auto *smg_125 = buffer.data(smg + 125);
    const auto *smg_130 = buffer.data(smg + 130);
    const auto *smg_131 = buffer.data(smg + 131);
    const auto *smg_132 = buffer.data(smg + 132);
    const auto *smg_133 = buffer.data(smg + 133);
    const auto *smg_134 = buffer.data(smg + 134);
    const auto *smg_135 = buffer.data(smg + 135);
    const auto *smg_137 = buffer.data(smg + 137);
    const auto *smg_138 = buffer.data(smg + 138);
    const auto *smg_140 = buffer.data(smg + 140);
    const auto *smg_141 = buffer.data(smg + 141);
    const auto *smg_144 = buffer.data(smg + 144);
    const auto *smg_145 = buffer.data(smg + 145);
    const auto *smg_146 = buffer.data(smg + 146);
    const auto *smg_147 = buffer.data(smg + 147);
    const auto *smg_148 = buffer.data(smg + 148);
    const auto *smg_149 = buffer.data(smg + 149);
    const auto *smg_150 = buffer.data(smg + 150);
    const auto *smg_152 = buffer.data(smg + 152);
    const auto *smg_153 = buffer.data(smg + 153);
    const auto *smg_155 = buffer.data(smg + 155);
    const auto *smg_156 = buffer.data(smg + 156);
    const auto *smg_159 = buffer.data(smg + 159);
    const auto *smg_160 = buffer.data(smg + 160);
    const auto *smg_161 = buffer.data(smg + 161);
    const auto *smg_162 = buffer.data(smg + 162);
    const auto *smg_163 = buffer.data(smg + 163);
    const auto *smg_164 = buffer.data(smg + 164);
    const auto *smg_165 = buffer.data(smg + 165);
    const auto *smg_167 = buffer.data(smg + 167);
    const auto *smg_168 = buffer.data(smg + 168);
    const auto *smg_170 = buffer.data(smg + 170);
    const auto *smg_174 = buffer.data(smg + 174);
    const auto *smg_175 = buffer.data(smg + 175);
    const auto *smg_176 = buffer.data(smg + 176);

#pragma omp simd aligned(t_126, t_127, t_128, t_129, pc_x, pc_y, pc_z, slg_45, slg_90, slg_93, \
                         smf0_60, smf0_63, smf1_60, smf1_63, smg_90, \
                         smg_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_14 * slg_90[k]
                   + f_1 * smf0_60[k]
                   - f_2 * smf1_60[k]
                   + f_3 * pc_x[k] * smg_90[k];

        t_127[k] = f_11 * slg_45[k]
                   + f_3 * pc_y[k] * smg_90[k];

        t_128[k] = f_3 * pc_z[k] * smg_90[k];

        t_129[k] = f_14 * slg_93[k]
                   + f_4 * smf0_63[k]
                   - f_5 * smf1_63[k]
                   + f_3 * pc_x[k] * smg_93[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, pc_x, pc_y, slg_47, slg_95, slg_96, smf0_65, \
                         smf0_66, smf1_65, smf1_66, smg_92, smg_95, \
                         smg_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_11 * slg_47[k]
                   + f_3 * pc_y[k] * smg_92[k];

        t_131[k] = f_14 * slg_95[k]
                   + f_4 * smf0_65[k]
                   - f_5 * smf1_65[k]
                   + f_3 * pc_x[k] * smg_95[k];

        t_132[k] = f_14 * slg_96[k]
                   + f_6 * smf0_66[k]
                   - f_7 * smf1_66[k]
                   + f_3 * pc_x[k] * smg_96[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, pc_x, pc_y, pc_z, slg_50, slg_99, \
                         slg_100, smf0_69, smf1_69, smg_93, smg_95, smg_99, \
                         smg_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_3 * pc_z[k] * smg_93[k];

        t_134[k] = f_11 * slg_50[k]
                   + f_3 * pc_y[k] * smg_95[k];

        t_135[k] = f_14 * slg_99[k]
                   + f_6 * smf0_69[k]
                   - f_7 * smf1_69[k]
                   + f_3 * pc_x[k] * smg_99[k];

        t_136[k] = f_14 * slg_100[k]
                   + f_3 * pc_x[k] * smg_100[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, t_140, pc_x, slg_101, slg_102, slg_103, slg_104, \
                         smg_101, smg_102, smg_103, smg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_14 * slg_101[k]
                   + f_3 * pc_x[k] * smg_101[k];

        t_138[k] = f_14 * slg_102[k]
                   + f_3 * pc_x[k] * smg_102[k];

        t_139[k] = f_14 * slg_103[k]
                   + f_3 * pc_x[k] * smg_103[k];

        t_140[k] = f_14 * slg_104[k]
                   + f_3 * pc_x[k] * smg_104[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, pc_y, pc_z, slg_55, slg_57, smf0_66, smf0_68, \
                         smf1_66, smf1_68, smg_100, smg_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_11 * slg_55[k]
                   + f_1 * smf0_66[k]
                   - f_2 * smf1_66[k]
                   + f_3 * pc_y[k] * smg_100[k];

        t_142[k] = f_3 * pc_z[k] * smg_100[k];

        t_143[k] = f_11 * slg_57[k]
                   + f_4 * smf0_68[k]
                   - f_5 * smf1_68[k]
                   + f_3 * pc_y[k] * smg_102[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, pb_z, pc_y, pc_z, slh0_63, slg_58, \
                         slg_59, slh1_63, smf0_69, smf1_69, smg_103, \
                         smg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_11 * slg_58[k]
                   + f_6 * smf0_69[k]
                   - f_7 * smf1_69[k]
                   + f_3 * pc_y[k] * smg_103[k];

        t_145[k] = f_11 * slg_59[k]
                   + f_3 * pc_y[k] * smg_104[k];

        t_146[k] = f_1 * smf0_69[k]
                   - f_2 * smf1_69[k]
                   + f_3 * pc_z[k] * smg_104[k];

        t_147[k] = pb_z[k] * slh0_63[k]
                   - f_8 * pc_z[k] * slh1_63[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, t_151, pb_z, pc_y, pc_z, slh0_66, slg_45, \
                         slg_60, slg_62, slh1_66, smg_105, smg_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_10 * slg_60[k]
                   + f_3 * pc_y[k] * smg_105[k];

        t_149[k] = f_9 * slg_45[k]
                   + f_3 * pc_z[k] * smg_105[k];

        t_150[k] = pb_z[k] * slh0_66[k]
                   - f_8 * pc_z[k] * slh1_66[k];

        t_151[k] = f_10 * slg_62[k]
                   + f_3 * pc_y[k] * smg_107[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, pb_z, pc_x, pc_z, slh0_69, slg_48, slg_110, \
                         slh1_69, smf0_75, smf1_75, smg_108, smg_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_14 * slg_110[k]
                   + f_4 * smf0_75[k]
                   - f_5 * smf1_75[k]
                   + f_3 * pc_x[k] * smg_110[k];

        t_153[k] = pb_z[k] * slh0_69[k]
                   - f_8 * pc_z[k] * slh1_69[k];

        t_154[k] = f_9 * slg_48[k]
                   + f_3 * pc_z[k] * smg_108[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, pc_x, pc_y, slg_65, slg_114, slg_115, \
                         slg_116, smf0_79, smf1_79, smg_110, smg_114, smg_115, \
                         smg_116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_10 * slg_65[k]
                   + f_3 * pc_y[k] * smg_110[k];

        t_156[k] = f_14 * slg_114[k]
                   + f_6 * smf0_79[k]
                   - f_7 * smf1_79[k]
                   + f_3 * pc_x[k] * smg_114[k];

        t_157[k] = f_14 * slg_115[k]
                   + f_3 * pc_x[k] * smg_115[k];

        t_158[k] = f_14 * slg_116[k]
                   + f_3 * pc_x[k] * smg_116[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, pb_z, pc_x, pc_z, slh0_78, slg_117, \
                         slg_118, slg_119, slh1_78, smg_117, smg_118, \
                         smg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_14 * slg_117[k]
                   + f_3 * pc_x[k] * smg_117[k];

        t_160[k] = f_14 * slg_118[k]
                   + f_3 * pc_x[k] * smg_118[k];

        t_161[k] = f_14 * slg_119[k]
                   + f_3 * pc_x[k] * smg_119[k];

        t_162[k] = pb_z[k] * slh0_78[k]
                   - f_8 * pc_z[k] * slh1_78[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, pc_y, pc_z, slg_55, slg_72, slg_73, smf0_78, \
                         smf0_79, smf1_78, smf1_79, smg_115, smg_117, \
                         smg_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_9 * slg_55[k]
                   + f_3 * pc_z[k] * smg_115[k];

        t_164[k] = f_10 * slg_72[k]
                   + f_4 * smf0_78[k]
                   - f_5 * smf1_78[k]
                   + f_3 * pc_y[k] * smg_117[k];

        t_165[k] = f_10 * slg_73[k]
                   + f_6 * smf0_79[k]
                   - f_7 * smf1_79[k]
                   + f_3 * pc_y[k] * smg_118[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, pb_y, pc_y, pc_z, slh0_105, slg_59, \
                         slg_74, slg_75, slh1_105, smf0_79, smf1_79, smg_119, \
                         smg_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_10 * slg_74[k]
                   + f_3 * pc_y[k] * smg_119[k];

        t_167[k] = f_9 * slg_59[k]
                   + f_1 * smf0_79[k]
                   - f_2 * smf1_79[k]
                   + f_3 * pc_z[k] * smg_119[k];

        t_168[k] = pb_y[k] * slh0_105[k]
                   - f_8 * pc_y[k] * slh1_105[k];

        t_169[k] = f_9 * slg_75[k]
                   + f_3 * pc_y[k] * smg_120[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, pb_y, pc_y, pc_z, slh0_108, slh0_110, \
                         slg_60, slg_76, slg_77, slh1_108, slh1_110, smg_120, \
                         smg_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_10 * slg_60[k]
                   + f_3 * pc_z[k] * smg_120[k];

        t_171[k] = pb_y[k] * slh0_108[k]
                   + f_10 * slg_76[k]
                   - f_8 * pc_y[k] * slh1_108[k];

        t_172[k] = f_9 * slg_77[k]
                   + f_3 * pc_y[k] * smg_122[k];

        t_173[k] = pb_y[k] * slh0_110[k]
                   - f_8 * pc_y[k] * slh1_110[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, pb_y, pc_y, pc_z, slh0_111, slh0_114, \
                         slg_63, slg_78, slg_80, slh1_111, slh1_114, smg_123, \
                         smg_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = pb_y[k] * slh0_111[k]
                   + f_11 * slg_78[k]
                   - f_8 * pc_y[k] * slh1_111[k];

        t_175[k] = f_10 * slg_63[k]
                   + f_3 * pc_z[k] * smg_123[k];

        t_176[k] = f_9 * slg_80[k]
                   + f_3 * pc_y[k] * smg_125[k];

        t_177[k] = pb_y[k] * slh0_114[k]
                   - f_8 * pc_y[k] * slh1_114[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, t_181, t_182, pc_x, slg_130, slg_131, slg_132, \
                         slg_133, slg_134, smg_130, smg_131, smg_132, smg_133, \
                         smg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_14 * slg_130[k]
                   + f_3 * pc_x[k] * smg_130[k];

        t_179[k] = f_14 * slg_131[k]
                   + f_3 * pc_x[k] * smg_131[k];

        t_180[k] = f_14 * slg_132[k]
                   + f_3 * pc_x[k] * smg_132[k];

        t_181[k] = f_14 * slg_133[k]
                   + f_3 * pc_x[k] * smg_133[k];

        t_182[k] = f_14 * slg_134[k]
                   + f_3 * pc_x[k] * smg_134[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, pc_y, pc_z, slg_70, slg_85, slg_87, smf0_86, \
                         smf0_88, smf1_86, smf1_88, smg_130, smg_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_9 * slg_85[k]
                   + f_1 * smf0_86[k]
                   - f_2 * smf1_86[k]
                   + f_3 * pc_y[k] * smg_130[k];

        t_184[k] = f_10 * slg_70[k]
                   + f_3 * pc_z[k] * smg_130[k];

        t_185[k] = f_9 * slg_87[k]
                   + f_4 * smf0_88[k]
                   - f_5 * smf1_88[k]
                   + f_3 * pc_y[k] * smg_132[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, pb_y, pc_y, slh0_125, slg_88, slg_89, slh1_125, \
                         smf0_89, smf1_89, smg_133, smg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_9 * slg_88[k]
                   + f_6 * smf0_89[k]
                   - f_7 * smf1_89[k]
                   + f_3 * pc_y[k] * smg_133[k];

        t_187[k] = f_9 * slg_89[k]
                   + f_3 * pc_y[k] * smg_134[k];

        t_188[k] = pb_y[k] * slh0_125[k]
                   - f_8 * pc_y[k] * slh1_125[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, pc_x, pc_y, pc_z, slg_75, slg_135, \
                         slg_138, smf0_90, smf0_93, smf1_90, smf1_93, smg_135, \
                         smg_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_14 * slg_135[k]
                   + f_1 * smf0_90[k]
                   - f_2 * smf1_90[k]
                   + f_3 * pc_x[k] * smg_135[k];

        t_190[k] = f_3 * pc_y[k] * smg_135[k];

        t_191[k] = f_11 * slg_75[k]
                   + f_3 * pc_z[k] * smg_135[k];

        t_192[k] = f_14 * slg_138[k]
                   + f_4 * smf0_93[k]
                   - f_5 * smf1_93[k]
                   + f_3 * pc_x[k] * smg_138[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, pc_x, pc_y, slg_140, slg_141, smf0_95, smf0_96, \
                         smf1_95, smf1_96, smg_137, smg_140, smg_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_3 * pc_y[k] * smg_137[k];

        t_194[k] = f_14 * slg_140[k]
                   + f_4 * smf0_95[k]
                   - f_5 * smf1_95[k]
                   + f_3 * pc_x[k] * smg_140[k];

        t_195[k] = f_14 * slg_141[k]
                   + f_6 * smf0_96[k]
                   - f_7 * smf1_96[k]
                   + f_3 * pc_x[k] * smg_141[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, pc_x, pc_y, pc_z, slg_78, slg_144, \
                         slg_145, smf0_99, smf1_99, smg_138, smg_140, smg_144, \
                         smg_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = f_11 * slg_78[k]
                   + f_3 * pc_z[k] * smg_138[k];

        t_197[k] = f_3 * pc_y[k] * smg_140[k];

        t_198[k] = f_14 * slg_144[k]
                   + f_6 * smf0_99[k]
                   - f_7 * smf1_99[k]
                   + f_3 * pc_x[k] * smg_144[k];

        t_199[k] = f_14 * slg_145[k]
                   + f_3 * pc_x[k] * smg_145[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, pc_x, slg_146, slg_147, slg_148, slg_149, \
                         smg_146, smg_147, smg_148, smg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_14 * slg_146[k]
                   + f_3 * pc_x[k] * smg_146[k];

        t_201[k] = f_14 * slg_147[k]
                   + f_3 * pc_x[k] * smg_147[k];

        t_202[k] = f_14 * slg_148[k]
                   + f_3 * pc_x[k] * smg_148[k];

        t_203[k] = f_14 * slg_149[k]
                   + f_3 * pc_x[k] * smg_149[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, pc_y, pc_z, slg_85, smf0_96, smf0_98, \
                         smf0_99, smf1_96, smf1_98, smf1_99, smg_145, smg_147, \
                         smg_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_1 * smf0_96[k]
                   - f_2 * smf1_96[k]
                   + f_3 * pc_y[k] * smg_145[k];

        t_205[k] = f_11 * slg_85[k]
                   + f_3 * pc_z[k] * smg_145[k];

        t_206[k] = f_4 * smf0_98[k]
                   - f_5 * smf1_98[k]
                   + f_3 * pc_y[k] * smg_147[k];

        t_207[k] = f_6 * smf0_99[k]
                   - f_7 * smf1_99[k]
                   + f_3 * pc_y[k] * smg_148[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, pc_x, pc_y, pc_z, slg_89, slg_90, \
                         slg_150, smf0_99, smf0_100, smf1_99, smf1_100, smg_149, \
                         smg_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_3 * pc_y[k] * smg_149[k];

        t_209[k] = f_11 * slg_89[k]
                   + f_1 * smf0_99[k]
                   - f_2 * smf1_99[k]
                   + f_3 * pc_z[k] * smg_149[k];

        t_210[k] = f_15 * slg_150[k]
                   + f_1 * smf0_100[k]
                   - f_2 * smf1_100[k]
                   + f_3 * pc_x[k] * smg_150[k];

        t_211[k] = f_16 * slg_90[k]
                   + f_3 * pc_y[k] * smg_150[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, pc_x, pc_y, pc_z, slg_92, slg_153, smf0_103, \
                         smf1_103, smg_150, smg_152, smg_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_3 * pc_z[k] * smg_150[k];

        t_213[k] = f_15 * slg_153[k]
                   + f_4 * smf0_103[k]
                   - f_5 * smf1_103[k]
                   + f_3 * pc_x[k] * smg_153[k];

        t_214[k] = f_16 * slg_92[k]
                   + f_3 * pc_y[k] * smg_152[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, pc_x, pc_z, slg_155, slg_156, smf0_105, \
                         smf0_106, smf1_105, smf1_106, smg_153, smg_155, \
                         smg_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = f_15 * slg_155[k]
                   + f_4 * smf0_105[k]
                   - f_5 * smf1_105[k]
                   + f_3 * pc_x[k] * smg_155[k];

        t_216[k] = f_15 * slg_156[k]
                   + f_6 * smf0_106[k]
                   - f_7 * smf1_106[k]
                   + f_3 * pc_x[k] * smg_156[k];

        t_217[k] = f_3 * pc_z[k] * smg_153[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, t_221, pc_x, pc_y, slg_95, slg_159, slg_160, \
                         slg_161, smf0_109, smf1_109, smg_155, smg_159, smg_160, \
                         smg_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_16 * slg_95[k]
                   + f_3 * pc_y[k] * smg_155[k];

        t_219[k] = f_15 * slg_159[k]
                   + f_6 * smf0_109[k]
                   - f_7 * smf1_109[k]
                   + f_3 * pc_x[k] * smg_159[k];

        t_220[k] = f_15 * slg_160[k]
                   + f_3 * pc_x[k] * smg_160[k];

        t_221[k] = f_15 * slg_161[k]
                   + f_3 * pc_x[k] * smg_161[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, pc_x, pc_y, slg_100, slg_162, slg_163, \
                         slg_164, smf0_106, smf1_106, smg_160, smg_162, smg_163, \
                         smg_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_15 * slg_162[k]
                   + f_3 * pc_x[k] * smg_162[k];

        t_223[k] = f_15 * slg_163[k]
                   + f_3 * pc_x[k] * smg_163[k];

        t_224[k] = f_15 * slg_164[k]
                   + f_3 * pc_x[k] * smg_164[k];

        t_225[k] = f_16 * slg_100[k]
                   + f_1 * smf0_106[k]
                   - f_2 * smf1_106[k]
                   + f_3 * pc_y[k] * smg_160[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, pc_y, pc_z, slg_102, slg_103, smf0_108, \
                         smf0_109, smf1_108, smf1_109, smg_160, smg_162, \
                         smg_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = f_3 * pc_z[k] * smg_160[k];

        t_227[k] = f_16 * slg_102[k]
                   + f_4 * smf0_108[k]
                   - f_5 * smf1_108[k]
                   + f_3 * pc_y[k] * smg_162[k];

        t_228[k] = f_16 * slg_103[k]
                   + f_6 * smf0_109[k]
                   - f_7 * smf1_109[k]
                   + f_3 * pc_y[k] * smg_163[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, pb_z, pc_y, pc_z, slh0_126, slg_104, \
                         slg_105, slh1_126, smf0_109, smf1_109, smg_164, \
                         smg_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = f_16 * slg_104[k]
                   + f_3 * pc_y[k] * smg_164[k];

        t_230[k] = f_1 * smf0_109[k]
                   - f_2 * smf1_109[k]
                   + f_3 * pc_z[k] * smg_164[k];

        t_231[k] = pb_z[k] * slh0_126[k]
                   - f_8 * pc_z[k] * slh1_126[k];

        t_232[k] = f_11 * slg_105[k]
                   + f_3 * pc_y[k] * smg_165[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, pb_z, pc_y, pc_z, slh0_129, slg_90, slg_107, \
                         slh1_129, smg_165, smg_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = f_9 * slg_90[k]
                   + f_3 * pc_z[k] * smg_165[k];

        t_234[k] = pb_z[k] * slh0_129[k]
                   - f_8 * pc_z[k] * slh1_129[k];

        t_235[k] = f_11 * slg_107[k]
                   + f_3 * pc_y[k] * smg_167[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, pb_z, pc_x, pc_z, slh0_132, slg_93, slg_170, \
                         slh1_132, smf0_115, smf1_115, smg_168, \
                         smg_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = f_15 * slg_170[k]
                   + f_4 * smf0_115[k]
                   - f_5 * smf1_115[k]
                   + f_3 * pc_x[k] * smg_170[k];

        t_237[k] = pb_z[k] * slh0_132[k]
                   - f_8 * pc_z[k] * slh1_132[k];

        t_238[k] = f_9 * slg_93[k]
                   + f_3 * pc_z[k] * smg_168[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, pc_x, pc_y, slg_110, slg_174, slg_175, \
                         slg_176, smf0_119, smf1_119, smg_170, smg_174, smg_175, \
                         smg_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_11 * slg_110[k]
                   + f_3 * pc_y[k] * smg_170[k];

        t_240[k] = f_15 * slg_174[k]
                   + f_6 * smf0_119[k]
                   - f_7 * smf1_119[k]
                   + f_3 * pc_x[k] * smg_174[k];

        t_241[k] = f_15 * slg_175[k]
                   + f_3 * pc_x[k] * smg_175[k];

        t_242[k] = f_15 * slg_176[k]
                   + f_3 * pc_x[k] * smg_176[k];
    }
}

static auto
compute_prim_smh_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t slh0,
                                                          const size_t slg, const size_t slh1,
                                                          const size_t smf0, const size_t smf1,
                                                          const size_t smg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.0 / gamma;
    const auto f_5 = p / (gamma * q);
    const auto f_6 = 0.5 / gamma;
    const auto f_7 = 0.5 * p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_15 = 2.5 / q;
    const auto f_16 = 2.0 / q;

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
    auto *t_353 = buffer.data(target + 353);
    auto *t_354 = buffer.data(target + 354);
    auto *t_355 = buffer.data(target + 355);
    auto *t_356 = buffer.data(target + 356);
    auto *t_357 = buffer.data(target + 357);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *slh0_141 = buffer.data(slh0 + 141);
    const auto *slh0_189 = buffer.data(slh0 + 189);
    const auto *slh0_192 = buffer.data(slh0 + 192);
    const auto *slh0_194 = buffer.data(slh0 + 194);
    const auto *slh0_195 = buffer.data(slh0 + 195);
    const auto *slh0_198 = buffer.data(slh0 + 198);
    const auto *slh0_209 = buffer.data(slh0 + 209);
    const auto *slh0_210 = buffer.data(slh0 + 210);
    const auto *slh0_213 = buffer.data(slh0 + 213);
    const auto *slh0_216 = buffer.data(slh0 + 216);
    const auto *slh0_225 = buffer.data(slh0 + 225);

    const auto *slg_100 = buffer.data(slg + 100);
    const auto *slg_104 = buffer.data(slg + 104);
    const auto *slg_105 = buffer.data(slg + 105);
    const auto *slg_108 = buffer.data(slg + 108);
    const auto *slg_115 = buffer.data(slg + 115);
    const auto *slg_117 = buffer.data(slg + 117);
    const auto *slg_118 = buffer.data(slg + 118);
    const auto *slg_119 = buffer.data(slg + 119);
    const auto *slg_120 = buffer.data(slg + 120);
    const auto *slg_122 = buffer.data(slg + 122);
    const auto *slg_123 = buffer.data(slg + 123);
    const auto *slg_125 = buffer.data(slg + 125);
    const auto *slg_130 = buffer.data(slg + 130);
    const auto *slg_132 = buffer.data(slg + 132);
    const auto *slg_133 = buffer.data(slg + 133);
    const auto *slg_134 = buffer.data(slg + 134);
    const auto *slg_135 = buffer.data(slg + 135);
    const auto *slg_136 = buffer.data(slg + 136);
    const auto *slg_137 = buffer.data(slg + 137);
    const auto *slg_138 = buffer.data(slg + 138);
    const auto *slg_140 = buffer.data(slg + 140);
    const auto *slg_145 = buffer.data(slg + 145);
    const auto *slg_147 = buffer.data(slg + 147);
    const auto *slg_148 = buffer.data(slg + 148);
    const auto *slg_149 = buffer.data(slg + 149);
    const auto *slg_150 = buffer.data(slg + 150);
    const auto *slg_152 = buffer.data(slg + 152);
    const auto *slg_153 = buffer.data(slg + 153);
    const auto *slg_155 = buffer.data(slg + 155);
    const auto *slg_160 = buffer.data(slg + 160);
    const auto *slg_162 = buffer.data(slg + 162);
    const auto *slg_163 = buffer.data(slg + 163);
    const auto *slg_164 = buffer.data(slg + 164);
    const auto *slg_165 = buffer.data(slg + 165);
    const auto *slg_167 = buffer.data(slg + 167);
    const auto *slg_170 = buffer.data(slg + 170);
    const auto *slg_177 = buffer.data(slg + 177);
    const auto *slg_178 = buffer.data(slg + 178);
    const auto *slg_179 = buffer.data(slg + 179);
    const auto *slg_180 = buffer.data(slg + 180);
    const auto *slg_183 = buffer.data(slg + 183);
    const auto *slg_185 = buffer.data(slg + 185);
    const auto *slg_186 = buffer.data(slg + 186);
    const auto *slg_189 = buffer.data(slg + 189);
    const auto *slg_190 = buffer.data(slg + 190);
    const auto *slg_191 = buffer.data(slg + 191);
    const auto *slg_192 = buffer.data(slg + 192);
    const auto *slg_193 = buffer.data(slg + 193);
    const auto *slg_194 = buffer.data(slg + 194);
    const auto *slg_205 = buffer.data(slg + 205);
    const auto *slg_206 = buffer.data(slg + 206);
    const auto *slg_207 = buffer.data(slg + 207);
    const auto *slg_208 = buffer.data(slg + 208);
    const auto *slg_209 = buffer.data(slg + 209);
    const auto *slg_210 = buffer.data(slg + 210);
    const auto *slg_213 = buffer.data(slg + 213);
    const auto *slg_215 = buffer.data(slg + 215);
    const auto *slg_216 = buffer.data(slg + 216);
    const auto *slg_219 = buffer.data(slg + 219);
    const auto *slg_220 = buffer.data(slg + 220);
    const auto *slg_221 = buffer.data(slg + 221);
    const auto *slg_222 = buffer.data(slg + 222);
    const auto *slg_223 = buffer.data(slg + 223);
    const auto *slg_224 = buffer.data(slg + 224);
    const auto *slg_225 = buffer.data(slg + 225);
    const auto *slg_228 = buffer.data(slg + 228);
    const auto *slg_230 = buffer.data(slg + 230);
    const auto *slg_231 = buffer.data(slg + 231);
    const auto *slg_234 = buffer.data(slg + 234);
    const auto *slg_235 = buffer.data(slg + 235);
    const auto *slg_236 = buffer.data(slg + 236);
    const auto *slg_237 = buffer.data(slg + 237);
    const auto *slg_238 = buffer.data(slg + 238);
    const auto *slg_239 = buffer.data(slg + 239);
    const auto *slg_245 = buffer.data(slg + 245);
    const auto *slg_249 = buffer.data(slg + 249);
    const auto *slg_250 = buffer.data(slg + 250);
    const auto *slg_251 = buffer.data(slg + 251);
    const auto *slg_252 = buffer.data(slg + 252);
    const auto *slg_253 = buffer.data(slg + 253);
    const auto *slg_254 = buffer.data(slg + 254);
    const auto *slg_255 = buffer.data(slg + 255);

    const auto *slh1_141 = buffer.data(slh1 + 141);
    const auto *slh1_189 = buffer.data(slh1 + 189);
    const auto *slh1_192 = buffer.data(slh1 + 192);
    const auto *slh1_194 = buffer.data(slh1 + 194);
    const auto *slh1_195 = buffer.data(slh1 + 195);
    const auto *slh1_198 = buffer.data(slh1 + 198);
    const auto *slh1_209 = buffer.data(slh1 + 209);
    const auto *slh1_210 = buffer.data(slh1 + 210);
    const auto *slh1_213 = buffer.data(slh1 + 213);
    const auto *slh1_216 = buffer.data(slh1 + 216);
    const auto *slh1_225 = buffer.data(slh1 + 225);

    const auto *smf0_118 = buffer.data(smf0 + 118);
    const auto *smf0_119 = buffer.data(smf0 + 119);
    const auto *smf0_120 = buffer.data(smf0 + 120);
    const auto *smf0_123 = buffer.data(smf0 + 123);
    const auto *smf0_125 = buffer.data(smf0 + 125);
    const auto *smf0_126 = buffer.data(smf0 + 126);
    const auto *smf0_128 = buffer.data(smf0 + 128);
    const auto *smf0_129 = buffer.data(smf0 + 129);
    const auto *smf0_136 = buffer.data(smf0 + 136);
    const auto *smf0_138 = buffer.data(smf0 + 138);
    const auto *smf0_139 = buffer.data(smf0 + 139);
    const auto *smf0_140 = buffer.data(smf0 + 140);
    const auto *smf0_143 = buffer.data(smf0 + 143);
    const auto *smf0_145 = buffer.data(smf0 + 145);
    const auto *smf0_146 = buffer.data(smf0 + 146);
    const auto *smf0_148 = buffer.data(smf0 + 148);
    const auto *smf0_149 = buffer.data(smf0 + 149);
    const auto *smf0_150 = buffer.data(smf0 + 150);
    const auto *smf0_153 = buffer.data(smf0 + 153);
    const auto *smf0_155 = buffer.data(smf0 + 155);
    const auto *smf0_156 = buffer.data(smf0 + 156);
    const auto *smf0_158 = buffer.data(smf0 + 158);
    const auto *smf0_159 = buffer.data(smf0 + 159);
    const auto *smf0_165 = buffer.data(smf0 + 165);
    const auto *smf0_168 = buffer.data(smf0 + 168);
    const auto *smf0_169 = buffer.data(smf0 + 169);
    const auto *smf0_170 = buffer.data(smf0 + 170);

    const auto *smf1_118 = buffer.data(smf1 + 118);
    const auto *smf1_119 = buffer.data(smf1 + 119);
    const auto *smf1_120 = buffer.data(smf1 + 120);
    const auto *smf1_123 = buffer.data(smf1 + 123);
    const auto *smf1_125 = buffer.data(smf1 + 125);
    const auto *smf1_126 = buffer.data(smf1 + 126);
    const auto *smf1_128 = buffer.data(smf1 + 128);
    const auto *smf1_129 = buffer.data(smf1 + 129);
    const auto *smf1_136 = buffer.data(smf1 + 136);
    const auto *smf1_138 = buffer.data(smf1 + 138);
    const auto *smf1_139 = buffer.data(smf1 + 139);
    const auto *smf1_140 = buffer.data(smf1 + 140);
    const auto *smf1_143 = buffer.data(smf1 + 143);
    const auto *smf1_145 = buffer.data(smf1 + 145);
    const auto *smf1_146 = buffer.data(smf1 + 146);
    const auto *smf1_148 = buffer.data(smf1 + 148);
    const auto *smf1_149 = buffer.data(smf1 + 149);
    const auto *smf1_150 = buffer.data(smf1 + 150);
    const auto *smf1_153 = buffer.data(smf1 + 153);
    const auto *smf1_155 = buffer.data(smf1 + 155);
    const auto *smf1_156 = buffer.data(smf1 + 156);
    const auto *smf1_158 = buffer.data(smf1 + 158);
    const auto *smf1_159 = buffer.data(smf1 + 159);
    const auto *smf1_165 = buffer.data(smf1 + 165);
    const auto *smf1_168 = buffer.data(smf1 + 168);
    const auto *smf1_169 = buffer.data(smf1 + 169);
    const auto *smf1_170 = buffer.data(smf1 + 170);

    const auto *smg_175 = buffer.data(smg + 175);
    const auto *smg_177 = buffer.data(smg + 177);
    const auto *smg_178 = buffer.data(smg + 178);
    const auto *smg_179 = buffer.data(smg + 179);
    const auto *smg_180 = buffer.data(smg + 180);
    const auto *smg_182 = buffer.data(smg + 182);
    const auto *smg_183 = buffer.data(smg + 183);
    const auto *smg_185 = buffer.data(smg + 185);
    const auto *smg_186 = buffer.data(smg + 186);
    const auto *smg_189 = buffer.data(smg + 189);
    const auto *smg_190 = buffer.data(smg + 190);
    const auto *smg_191 = buffer.data(smg + 191);
    const auto *smg_192 = buffer.data(smg + 192);
    const auto *smg_193 = buffer.data(smg + 193);
    const auto *smg_194 = buffer.data(smg + 194);
    const auto *smg_195 = buffer.data(smg + 195);
    const auto *smg_197 = buffer.data(smg + 197);
    const auto *smg_198 = buffer.data(smg + 198);
    const auto *smg_200 = buffer.data(smg + 200);
    const auto *smg_205 = buffer.data(smg + 205);
    const auto *smg_206 = buffer.data(smg + 206);
    const auto *smg_207 = buffer.data(smg + 207);
    const auto *smg_208 = buffer.data(smg + 208);
    const auto *smg_209 = buffer.data(smg + 209);
    const auto *smg_210 = buffer.data(smg + 210);
    const auto *smg_212 = buffer.data(smg + 212);
    const auto *smg_213 = buffer.data(smg + 213);
    const auto *smg_215 = buffer.data(smg + 215);
    const auto *smg_216 = buffer.data(smg + 216);
    const auto *smg_219 = buffer.data(smg + 219);
    const auto *smg_220 = buffer.data(smg + 220);
    const auto *smg_221 = buffer.data(smg + 221);
    const auto *smg_222 = buffer.data(smg + 222);
    const auto *smg_223 = buffer.data(smg + 223);
    const auto *smg_224 = buffer.data(smg + 224);
    const auto *smg_225 = buffer.data(smg + 225);
    const auto *smg_227 = buffer.data(smg + 227);
    const auto *smg_228 = buffer.data(smg + 228);
    const auto *smg_230 = buffer.data(smg + 230);
    const auto *smg_231 = buffer.data(smg + 231);
    const auto *smg_234 = buffer.data(smg + 234);
    const auto *smg_235 = buffer.data(smg + 235);
    const auto *smg_236 = buffer.data(smg + 236);
    const auto *smg_237 = buffer.data(smg + 237);
    const auto *smg_238 = buffer.data(smg + 238);
    const auto *smg_239 = buffer.data(smg + 239);
    const auto *smg_240 = buffer.data(smg + 240);
    const auto *smg_242 = buffer.data(smg + 242);
    const auto *smg_243 = buffer.data(smg + 243);
    const auto *smg_245 = buffer.data(smg + 245);
    const auto *smg_249 = buffer.data(smg + 249);
    const auto *smg_250 = buffer.data(smg + 250);
    const auto *smg_251 = buffer.data(smg + 251);
    const auto *smg_252 = buffer.data(smg + 252);
    const auto *smg_253 = buffer.data(smg + 253);
    const auto *smg_254 = buffer.data(smg + 254);
    const auto *smg_255 = buffer.data(smg + 255);

#pragma omp simd aligned(t_243, t_244, t_245, t_246, pb_z, pc_x, pc_z, slh0_141, slg_177, \
                         slg_178, slg_179, slh1_141, smg_177, smg_178, \
                         smg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = f_15 * slg_177[k]
                   + f_3 * pc_x[k] * smg_177[k];

        t_244[k] = f_15 * slg_178[k]
                   + f_3 * pc_x[k] * smg_178[k];

        t_245[k] = f_15 * slg_179[k]
                   + f_3 * pc_x[k] * smg_179[k];

        t_246[k] = pb_z[k] * slh0_141[k]
                   - f_8 * pc_z[k] * slh1_141[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, pc_y, pc_z, slg_100, slg_117, slg_118, smf0_118, \
                         smf0_119, smf1_118, smf1_119, smg_175, smg_177, \
                         smg_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = f_9 * slg_100[k]
                   + f_3 * pc_z[k] * smg_175[k];

        t_248[k] = f_11 * slg_117[k]
                   + f_4 * smf0_118[k]
                   - f_5 * smf1_118[k]
                   + f_3 * pc_y[k] * smg_177[k];

        t_249[k] = f_11 * slg_118[k]
                   + f_6 * smf0_119[k]
                   - f_7 * smf1_119[k]
                   + f_3 * pc_y[k] * smg_178[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, pc_x, pc_y, pc_z, slg_104, slg_119, slg_180, \
                         smf0_119, smf0_120, smf1_119, smf1_120, smg_179, \
                         smg_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = f_11 * slg_119[k]
                   + f_3 * pc_y[k] * smg_179[k];

        t_251[k] = f_9 * slg_104[k]
                   + f_1 * smf0_119[k]
                   - f_2 * smf1_119[k]
                   + f_3 * pc_z[k] * smg_179[k];

        t_252[k] = f_15 * slg_180[k]
                   + f_1 * smf0_120[k]
                   - f_2 * smf1_120[k]
                   + f_3 * pc_x[k] * smg_180[k];
    }

#pragma omp simd aligned(t_253, t_254, t_255, t_256, pc_x, pc_y, pc_z, slg_105, slg_120, \
                         slg_122, slg_183, smf0_123, smf1_123, smg_180, smg_182, \
                         smg_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_253[k] = f_10 * slg_120[k]
                   + f_3 * pc_y[k] * smg_180[k];

        t_254[k] = f_10 * slg_105[k]
                   + f_3 * pc_z[k] * smg_180[k];

        t_255[k] = f_15 * slg_183[k]
                   + f_4 * smf0_123[k]
                   - f_5 * smf1_123[k]
                   + f_3 * pc_x[k] * smg_183[k];

        t_256[k] = f_10 * slg_122[k]
                   + f_3 * pc_y[k] * smg_182[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, pc_x, pc_z, slg_108, slg_185, slg_186, smf0_125, \
                         smf0_126, smf1_125, smf1_126, smg_183, smg_185, \
                         smg_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_15 * slg_185[k]
                   + f_4 * smf0_125[k]
                   - f_5 * smf1_125[k]
                   + f_3 * pc_x[k] * smg_185[k];

        t_258[k] = f_15 * slg_186[k]
                   + f_6 * smf0_126[k]
                   - f_7 * smf1_126[k]
                   + f_3 * pc_x[k] * smg_186[k];

        t_259[k] = f_10 * slg_108[k]
                   + f_3 * pc_z[k] * smg_183[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pc_x, pc_y, slg_125, slg_189, slg_190, \
                         slg_191, smf0_129, smf1_129, smg_185, smg_189, smg_190, \
                         smg_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_10 * slg_125[k]
                   + f_3 * pc_y[k] * smg_185[k];

        t_261[k] = f_15 * slg_189[k]
                   + f_6 * smf0_129[k]
                   - f_7 * smf1_129[k]
                   + f_3 * pc_x[k] * smg_189[k];

        t_262[k] = f_15 * slg_190[k]
                   + f_3 * pc_x[k] * smg_190[k];

        t_263[k] = f_15 * slg_191[k]
                   + f_3 * pc_x[k] * smg_191[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, pc_x, pc_y, slg_130, slg_192, slg_193, \
                         slg_194, smf0_126, smf1_126, smg_190, smg_192, smg_193, \
                         smg_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_15 * slg_192[k]
                   + f_3 * pc_x[k] * smg_192[k];

        t_265[k] = f_15 * slg_193[k]
                   + f_3 * pc_x[k] * smg_193[k];

        t_266[k] = f_15 * slg_194[k]
                   + f_3 * pc_x[k] * smg_194[k];

        t_267[k] = f_10 * slg_130[k]
                   + f_1 * smf0_126[k]
                   - f_2 * smf1_126[k]
                   + f_3 * pc_y[k] * smg_190[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, pc_y, pc_z, slg_115, slg_132, slg_133, smf0_128, \
                         smf0_129, smf1_128, smf1_129, smg_190, smg_192, \
                         smg_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_10 * slg_115[k]
                   + f_3 * pc_z[k] * smg_190[k];

        t_269[k] = f_10 * slg_132[k]
                   + f_4 * smf0_128[k]
                   - f_5 * smf1_128[k]
                   + f_3 * pc_y[k] * smg_192[k];

        t_270[k] = f_10 * slg_133[k]
                   + f_6 * smf0_129[k]
                   - f_7 * smf1_129[k]
                   + f_3 * pc_y[k] * smg_193[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, t_274, pb_y, pc_y, pc_z, slh0_189, slg_119, \
                         slg_134, slg_135, slh1_189, smf0_129, smf1_129, smg_194, \
                         smg_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = f_10 * slg_134[k]
                   + f_3 * pc_y[k] * smg_194[k];

        t_272[k] = f_10 * slg_119[k]
                   + f_1 * smf0_129[k]
                   - f_2 * smf1_129[k]
                   + f_3 * pc_z[k] * smg_194[k];

        t_273[k] = pb_y[k] * slh0_189[k]
                   - f_8 * pc_y[k] * slh1_189[k];

        t_274[k] = f_9 * slg_135[k]
                   + f_3 * pc_y[k] * smg_195[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, pb_y, pc_y, pc_z, slh0_192, slh0_194, \
                         slg_120, slg_136, slg_137, slh1_192, slh1_194, smg_195, \
                         smg_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = f_11 * slg_120[k]
                   + f_3 * pc_z[k] * smg_195[k];

        t_276[k] = pb_y[k] * slh0_192[k]
                   + f_10 * slg_136[k]
                   - f_8 * pc_y[k] * slh1_192[k];

        t_277[k] = f_9 * slg_137[k]
                   + f_3 * pc_y[k] * smg_197[k];

        t_278[k] = pb_y[k] * slh0_194[k]
                   - f_8 * pc_y[k] * slh1_194[k];
    }

#pragma omp simd aligned(t_279, t_280, t_281, t_282, pb_y, pc_y, pc_z, slh0_195, slh0_198, \
                         slg_123, slg_138, slg_140, slh1_195, slh1_198, smg_198, \
                         smg_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_279[k] = pb_y[k] * slh0_195[k]
                   + f_11 * slg_138[k]
                   - f_8 * pc_y[k] * slh1_195[k];

        t_280[k] = f_11 * slg_123[k]
                   + f_3 * pc_z[k] * smg_198[k];

        t_281[k] = f_9 * slg_140[k]
                   + f_3 * pc_y[k] * smg_200[k];

        t_282[k] = pb_y[k] * slh0_198[k]
                   - f_8 * pc_y[k] * slh1_198[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, t_286, t_287, pc_x, slg_205, slg_206, slg_207, \
                         slg_208, slg_209, smg_205, smg_206, smg_207, smg_208, \
                         smg_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_15 * slg_205[k]
                   + f_3 * pc_x[k] * smg_205[k];

        t_284[k] = f_15 * slg_206[k]
                   + f_3 * pc_x[k] * smg_206[k];

        t_285[k] = f_15 * slg_207[k]
                   + f_3 * pc_x[k] * smg_207[k];

        t_286[k] = f_15 * slg_208[k]
                   + f_3 * pc_x[k] * smg_208[k];

        t_287[k] = f_15 * slg_209[k]
                   + f_3 * pc_x[k] * smg_209[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, pc_y, pc_z, slg_130, slg_145, slg_147, smf0_136, \
                         smf0_138, smf1_136, smf1_138, smg_205, \
                         smg_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = f_9 * slg_145[k]
                   + f_1 * smf0_136[k]
                   - f_2 * smf1_136[k]
                   + f_3 * pc_y[k] * smg_205[k];

        t_289[k] = f_11 * slg_130[k]
                   + f_3 * pc_z[k] * smg_205[k];

        t_290[k] = f_9 * slg_147[k]
                   + f_4 * smf0_138[k]
                   - f_5 * smf1_138[k]
                   + f_3 * pc_y[k] * smg_207[k];
    }

#pragma omp simd aligned(t_291, t_292, t_293, pb_y, pc_y, slh0_209, slg_148, slg_149, \
                         slh1_209, smf0_139, smf1_139, smg_208, \
                         smg_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_291[k] = f_9 * slg_148[k]
                   + f_6 * smf0_139[k]
                   - f_7 * smf1_139[k]
                   + f_3 * pc_y[k] * smg_208[k];

        t_292[k] = f_9 * slg_149[k]
                   + f_3 * pc_y[k] * smg_209[k];

        t_293[k] = pb_y[k] * slh0_209[k]
                   - f_8 * pc_y[k] * slh1_209[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, pc_x, pc_y, pc_z, slg_135, slg_210, \
                         slg_213, smf0_140, smf0_143, smf1_140, smf1_143, smg_210, \
                         smg_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_15 * slg_210[k]
                   + f_1 * smf0_140[k]
                   - f_2 * smf1_140[k]
                   + f_3 * pc_x[k] * smg_210[k];

        t_295[k] = f_3 * pc_y[k] * smg_210[k];

        t_296[k] = f_16 * slg_135[k]
                   + f_3 * pc_z[k] * smg_210[k];

        t_297[k] = f_15 * slg_213[k]
                   + f_4 * smf0_143[k]
                   - f_5 * smf1_143[k]
                   + f_3 * pc_x[k] * smg_213[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, pc_x, pc_y, slg_215, slg_216, smf0_145, \
                         smf0_146, smf1_145, smf1_146, smg_212, smg_215, \
                         smg_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = f_3 * pc_y[k] * smg_212[k];

        t_299[k] = f_15 * slg_215[k]
                   + f_4 * smf0_145[k]
                   - f_5 * smf1_145[k]
                   + f_3 * pc_x[k] * smg_215[k];

        t_300[k] = f_15 * slg_216[k]
                   + f_6 * smf0_146[k]
                   - f_7 * smf1_146[k]
                   + f_3 * pc_x[k] * smg_216[k];
    }

#pragma omp simd aligned(t_301, t_302, t_303, t_304, pc_x, pc_y, pc_z, slg_138, slg_219, \
                         slg_220, smf0_149, smf1_149, smg_213, smg_215, smg_219, \
                         smg_220 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_301[k] = f_16 * slg_138[k]
                   + f_3 * pc_z[k] * smg_213[k];

        t_302[k] = f_3 * pc_y[k] * smg_215[k];

        t_303[k] = f_15 * slg_219[k]
                   + f_6 * smf0_149[k]
                   - f_7 * smf1_149[k]
                   + f_3 * pc_x[k] * smg_219[k];

        t_304[k] = f_15 * slg_220[k]
                   + f_3 * pc_x[k] * smg_220[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, pc_x, slg_221, slg_222, slg_223, slg_224, \
                         smg_221, smg_222, smg_223, smg_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = f_15 * slg_221[k]
                   + f_3 * pc_x[k] * smg_221[k];

        t_306[k] = f_15 * slg_222[k]
                   + f_3 * pc_x[k] * smg_222[k];

        t_307[k] = f_15 * slg_223[k]
                   + f_3 * pc_x[k] * smg_223[k];

        t_308[k] = f_15 * slg_224[k]
                   + f_3 * pc_x[k] * smg_224[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, pc_y, pc_z, slg_145, smf0_146, smf0_148, \
                         smf0_149, smf1_146, smf1_148, smf1_149, smg_220, smg_222, \
                         smg_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_1 * smf0_146[k]
                   - f_2 * smf1_146[k]
                   + f_3 * pc_y[k] * smg_220[k];

        t_310[k] = f_16 * slg_145[k]
                   + f_3 * pc_z[k] * smg_220[k];

        t_311[k] = f_4 * smf0_148[k]
                   - f_5 * smf1_148[k]
                   + f_3 * pc_y[k] * smg_222[k];

        t_312[k] = f_6 * smf0_149[k]
                   - f_7 * smf1_149[k]
                   + f_3 * pc_y[k] * smg_223[k];
    }

#pragma omp simd aligned(t_313, t_314, t_315, t_316, pc_x, pc_y, pc_z, slg_149, slg_150, \
                         slg_225, smf0_149, smf0_150, smf1_149, smf1_150, smg_224, \
                         smg_225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_313[k] = f_3 * pc_y[k] * smg_224[k];

        t_314[k] = f_16 * slg_149[k]
                   + f_1 * smf0_149[k]
                   - f_2 * smf1_149[k]
                   + f_3 * pc_z[k] * smg_224[k];

        t_315[k] = f_16 * slg_225[k]
                   + f_1 * smf0_150[k]
                   - f_2 * smf1_150[k]
                   + f_3 * pc_x[k] * smg_225[k];

        t_316[k] = f_15 * slg_150[k]
                   + f_3 * pc_y[k] * smg_225[k];
    }

#pragma omp simd aligned(t_317, t_318, t_319, pc_x, pc_y, pc_z, slg_152, slg_228, smf0_153, \
                         smf1_153, smg_225, smg_227, smg_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_317[k] = f_3 * pc_z[k] * smg_225[k];

        t_318[k] = f_16 * slg_228[k]
                   + f_4 * smf0_153[k]
                   - f_5 * smf1_153[k]
                   + f_3 * pc_x[k] * smg_228[k];

        t_319[k] = f_15 * slg_152[k]
                   + f_3 * pc_y[k] * smg_227[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, pc_x, pc_z, slg_230, slg_231, smf0_155, \
                         smf0_156, smf1_155, smf1_156, smg_228, smg_230, \
                         smg_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = f_16 * slg_230[k]
                   + f_4 * smf0_155[k]
                   - f_5 * smf1_155[k]
                   + f_3 * pc_x[k] * smg_230[k];

        t_321[k] = f_16 * slg_231[k]
                   + f_6 * smf0_156[k]
                   - f_7 * smf1_156[k]
                   + f_3 * pc_x[k] * smg_231[k];

        t_322[k] = f_3 * pc_z[k] * smg_228[k];
    }

#pragma omp simd aligned(t_323, t_324, t_325, t_326, pc_x, pc_y, slg_155, slg_234, slg_235, \
                         slg_236, smf0_159, smf1_159, smg_230, smg_234, smg_235, \
                         smg_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_323[k] = f_15 * slg_155[k]
                   + f_3 * pc_y[k] * smg_230[k];

        t_324[k] = f_16 * slg_234[k]
                   + f_6 * smf0_159[k]
                   - f_7 * smf1_159[k]
                   + f_3 * pc_x[k] * smg_234[k];

        t_325[k] = f_16 * slg_235[k]
                   + f_3 * pc_x[k] * smg_235[k];

        t_326[k] = f_16 * slg_236[k]
                   + f_3 * pc_x[k] * smg_236[k];
    }

#pragma omp simd aligned(t_327, t_328, t_329, t_330, pc_x, pc_y, slg_160, slg_237, slg_238, \
                         slg_239, smf0_156, smf1_156, smg_235, smg_237, smg_238, \
                         smg_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_327[k] = f_16 * slg_237[k]
                   + f_3 * pc_x[k] * smg_237[k];

        t_328[k] = f_16 * slg_238[k]
                   + f_3 * pc_x[k] * smg_238[k];

        t_329[k] = f_16 * slg_239[k]
                   + f_3 * pc_x[k] * smg_239[k];

        t_330[k] = f_15 * slg_160[k]
                   + f_1 * smf0_156[k]
                   - f_2 * smf1_156[k]
                   + f_3 * pc_y[k] * smg_235[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, pc_y, pc_z, slg_162, slg_163, smf0_158, \
                         smf0_159, smf1_158, smf1_159, smg_235, smg_237, \
                         smg_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = f_3 * pc_z[k] * smg_235[k];

        t_332[k] = f_15 * slg_162[k]
                   + f_4 * smf0_158[k]
                   - f_5 * smf1_158[k]
                   + f_3 * pc_y[k] * smg_237[k];

        t_333[k] = f_15 * slg_163[k]
                   + f_6 * smf0_159[k]
                   - f_7 * smf1_159[k]
                   + f_3 * pc_y[k] * smg_238[k];
    }

#pragma omp simd aligned(t_334, t_335, t_336, t_337, pb_z, pc_y, pc_z, slh0_210, slg_164, \
                         slg_165, slh1_210, smf0_159, smf1_159, smg_239, \
                         smg_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_334[k] = f_15 * slg_164[k]
                   + f_3 * pc_y[k] * smg_239[k];

        t_335[k] = f_1 * smf0_159[k]
                   - f_2 * smf1_159[k]
                   + f_3 * pc_z[k] * smg_239[k];

        t_336[k] = pb_z[k] * slh0_210[k]
                   - f_8 * pc_z[k] * slh1_210[k];

        t_337[k] = f_16 * slg_165[k]
                   + f_3 * pc_y[k] * smg_240[k];
    }

#pragma omp simd aligned(t_338, t_339, t_340, pb_z, pc_y, pc_z, slh0_213, slg_150, slg_167, \
                         slh1_213, smg_240, smg_242 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_338[k] = f_9 * slg_150[k]
                   + f_3 * pc_z[k] * smg_240[k];

        t_339[k] = pb_z[k] * slh0_213[k]
                   - f_8 * pc_z[k] * slh1_213[k];

        t_340[k] = f_16 * slg_167[k]
                   + f_3 * pc_y[k] * smg_242[k];
    }

#pragma omp simd aligned(t_341, t_342, t_343, pb_z, pc_x, pc_z, slh0_216, slg_153, slg_245, \
                         slh1_216, smf0_165, smf1_165, smg_243, \
                         smg_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_341[k] = f_16 * slg_245[k]
                   + f_4 * smf0_165[k]
                   - f_5 * smf1_165[k]
                   + f_3 * pc_x[k] * smg_245[k];

        t_342[k] = pb_z[k] * slh0_216[k]
                   - f_8 * pc_z[k] * slh1_216[k];

        t_343[k] = f_9 * slg_153[k]
                   + f_3 * pc_z[k] * smg_243[k];
    }

#pragma omp simd aligned(t_344, t_345, t_346, t_347, pc_x, pc_y, slg_170, slg_249, slg_250, \
                         slg_251, smf0_169, smf1_169, smg_245, smg_249, smg_250, \
                         smg_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_344[k] = f_16 * slg_170[k]
                   + f_3 * pc_y[k] * smg_245[k];

        t_345[k] = f_16 * slg_249[k]
                   + f_6 * smf0_169[k]
                   - f_7 * smf1_169[k]
                   + f_3 * pc_x[k] * smg_249[k];

        t_346[k] = f_16 * slg_250[k]
                   + f_3 * pc_x[k] * smg_250[k];

        t_347[k] = f_16 * slg_251[k]
                   + f_3 * pc_x[k] * smg_251[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, t_351, pb_z, pc_x, pc_z, slh0_225, slg_252, \
                         slg_253, slg_254, slh1_225, smg_252, smg_253, \
                         smg_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = f_16 * slg_252[k]
                   + f_3 * pc_x[k] * smg_252[k];

        t_349[k] = f_16 * slg_253[k]
                   + f_3 * pc_x[k] * smg_253[k];

        t_350[k] = f_16 * slg_254[k]
                   + f_3 * pc_x[k] * smg_254[k];

        t_351[k] = pb_z[k] * slh0_225[k]
                   - f_8 * pc_z[k] * slh1_225[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, pc_y, pc_z, slg_160, slg_177, slg_178, smf0_168, \
                         smf0_169, smf1_168, smf1_169, smg_250, smg_252, \
                         smg_253 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = f_9 * slg_160[k]
                   + f_3 * pc_z[k] * smg_250[k];

        t_353[k] = f_16 * slg_177[k]
                   + f_4 * smf0_168[k]
                   - f_5 * smf1_168[k]
                   + f_3 * pc_y[k] * smg_252[k];

        t_354[k] = f_16 * slg_178[k]
                   + f_6 * smf0_169[k]
                   - f_7 * smf1_169[k]
                   + f_3 * pc_y[k] * smg_253[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, pc_x, pc_y, pc_z, slg_164, slg_179, slg_255, \
                         smf0_169, smf0_170, smf1_169, smf1_170, smg_254, \
                         smg_255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = f_16 * slg_179[k]
                   + f_3 * pc_y[k] * smg_254[k];

        t_356[k] = f_9 * slg_164[k]
                   + f_1 * smf0_169[k]
                   - f_2 * smf1_169[k]
                   + f_3 * pc_z[k] * smg_254[k];

        t_357[k] = f_16 * slg_255[k]
                   + f_1 * smf0_170[k]
                   - f_2 * smf1_170[k]
                   + f_3 * pc_x[k] * smg_255[k];
    }
}

static auto
compute_prim_smh_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t slh0,
                                                          const size_t slg, const size_t slh1,
                                                          const size_t smf0, const size_t smf1,
                                                          const size_t smg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.0 / gamma;
    const auto f_5 = p / (gamma * q);
    const auto f_6 = 0.5 / gamma;
    const auto f_7 = 0.5 * p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_14 = 3.0 / q;
    const auto f_15 = 2.5 / q;
    const auto f_16 = 2.0 / q;

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
    auto *t_468 = buffer.data(target + 468);
    auto *t_469 = buffer.data(target + 469);
    auto *t_470 = buffer.data(target + 470);
    auto *t_471 = buffer.data(target + 471);
    auto *t_472 = buffer.data(target + 472);
    auto *t_473 = buffer.data(target + 473);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *slh0_294 = buffer.data(slh0 + 294);
    const auto *slh0_297 = buffer.data(slh0 + 297);
    const auto *slh0_299 = buffer.data(slh0 + 299);
    const auto *slh0_300 = buffer.data(slh0 + 300);
    const auto *slh0_303 = buffer.data(slh0 + 303);
    const auto *slh0_314 = buffer.data(slh0 + 314);
    const auto *slh0_315 = buffer.data(slh0 + 315);
    const auto *slh0_318 = buffer.data(slh0 + 318);
    const auto *slh0_321 = buffer.data(slh0 + 321);

    const auto *slg_165 = buffer.data(slg + 165);
    const auto *slg_168 = buffer.data(slg + 168);
    const auto *slg_175 = buffer.data(slg + 175);
    const auto *slg_179 = buffer.data(slg + 179);
    const auto *slg_180 = buffer.data(slg + 180);
    const auto *slg_182 = buffer.data(slg + 182);
    const auto *slg_183 = buffer.data(slg + 183);
    const auto *slg_185 = buffer.data(slg + 185);
    const auto *slg_190 = buffer.data(slg + 190);
    const auto *slg_192 = buffer.data(slg + 192);
    const auto *slg_193 = buffer.data(slg + 193);
    const auto *slg_194 = buffer.data(slg + 194);
    const auto *slg_195 = buffer.data(slg + 195);
    const auto *slg_197 = buffer.data(slg + 197);
    const auto *slg_198 = buffer.data(slg + 198);
    const auto *slg_200 = buffer.data(slg + 200);
    const auto *slg_205 = buffer.data(slg + 205);
    const auto *slg_207 = buffer.data(slg + 207);
    const auto *slg_208 = buffer.data(slg + 208);
    const auto *slg_209 = buffer.data(slg + 209);
    const auto *slg_210 = buffer.data(slg + 210);
    const auto *slg_211 = buffer.data(slg + 211);
    const auto *slg_212 = buffer.data(slg + 212);
    const auto *slg_213 = buffer.data(slg + 213);
    const auto *slg_215 = buffer.data(slg + 215);
    const auto *slg_220 = buffer.data(slg + 220);
    const auto *slg_222 = buffer.data(slg + 222);
    const auto *slg_223 = buffer.data(slg + 223);
    const auto *slg_224 = buffer.data(slg + 224);
    const auto *slg_225 = buffer.data(slg + 225);
    const auto *slg_227 = buffer.data(slg + 227);
    const auto *slg_228 = buffer.data(slg + 228);
    const auto *slg_230 = buffer.data(slg + 230);
    const auto *slg_235 = buffer.data(slg + 235);
    const auto *slg_237 = buffer.data(slg + 237);
    const auto *slg_238 = buffer.data(slg + 238);
    const auto *slg_239 = buffer.data(slg + 239);
    const auto *slg_240 = buffer.data(slg + 240);
    const auto *slg_242 = buffer.data(slg + 242);
    const auto *slg_245 = buffer.data(slg + 245);
    const auto *slg_258 = buffer.data(slg + 258);
    const auto *slg_260 = buffer.data(slg + 260);
    const auto *slg_261 = buffer.data(slg + 261);
    const auto *slg_264 = buffer.data(slg + 264);
    const auto *slg_265 = buffer.data(slg + 265);
    const auto *slg_266 = buffer.data(slg + 266);
    const auto *slg_267 = buffer.data(slg + 267);
    const auto *slg_268 = buffer.data(slg + 268);
    const auto *slg_269 = buffer.data(slg + 269);
    const auto *slg_270 = buffer.data(slg + 270);
    const auto *slg_273 = buffer.data(slg + 273);
    const auto *slg_275 = buffer.data(slg + 275);
    const auto *slg_276 = buffer.data(slg + 276);
    const auto *slg_279 = buffer.data(slg + 279);
    const auto *slg_280 = buffer.data(slg + 280);
    const auto *slg_281 = buffer.data(slg + 281);
    const auto *slg_282 = buffer.data(slg + 282);
    const auto *slg_283 = buffer.data(slg + 283);
    const auto *slg_284 = buffer.data(slg + 284);
    const auto *slg_295 = buffer.data(slg + 295);
    const auto *slg_296 = buffer.data(slg + 296);
    const auto *slg_297 = buffer.data(slg + 297);
    const auto *slg_298 = buffer.data(slg + 298);
    const auto *slg_299 = buffer.data(slg + 299);
    const auto *slg_300 = buffer.data(slg + 300);
    const auto *slg_303 = buffer.data(slg + 303);
    const auto *slg_305 = buffer.data(slg + 305);
    const auto *slg_306 = buffer.data(slg + 306);
    const auto *slg_309 = buffer.data(slg + 309);
    const auto *slg_310 = buffer.data(slg + 310);
    const auto *slg_311 = buffer.data(slg + 311);
    const auto *slg_312 = buffer.data(slg + 312);
    const auto *slg_313 = buffer.data(slg + 313);
    const auto *slg_314 = buffer.data(slg + 314);
    const auto *slg_315 = buffer.data(slg + 315);
    const auto *slg_318 = buffer.data(slg + 318);
    const auto *slg_320 = buffer.data(slg + 320);
    const auto *slg_321 = buffer.data(slg + 321);
    const auto *slg_324 = buffer.data(slg + 324);
    const auto *slg_325 = buffer.data(slg + 325);
    const auto *slg_326 = buffer.data(slg + 326);
    const auto *slg_327 = buffer.data(slg + 327);
    const auto *slg_328 = buffer.data(slg + 328);
    const auto *slg_329 = buffer.data(slg + 329);
    const auto *slg_335 = buffer.data(slg + 335);
    const auto *slg_339 = buffer.data(slg + 339);
    const auto *slg_340 = buffer.data(slg + 340);
    const auto *slg_341 = buffer.data(slg + 341);

    const auto *slh1_294 = buffer.data(slh1 + 294);
    const auto *slh1_297 = buffer.data(slh1 + 297);
    const auto *slh1_299 = buffer.data(slh1 + 299);
    const auto *slh1_300 = buffer.data(slh1 + 300);
    const auto *slh1_303 = buffer.data(slh1 + 303);
    const auto *slh1_314 = buffer.data(slh1 + 314);
    const auto *slh1_315 = buffer.data(slh1 + 315);
    const auto *slh1_318 = buffer.data(slh1 + 318);
    const auto *slh1_321 = buffer.data(slh1 + 321);

    const auto *smf0_173 = buffer.data(smf0 + 173);
    const auto *smf0_175 = buffer.data(smf0 + 175);
    const auto *smf0_176 = buffer.data(smf0 + 176);
    const auto *smf0_178 = buffer.data(smf0 + 178);
    const auto *smf0_179 = buffer.data(smf0 + 179);
    const auto *smf0_180 = buffer.data(smf0 + 180);
    const auto *smf0_183 = buffer.data(smf0 + 183);
    const auto *smf0_185 = buffer.data(smf0 + 185);
    const auto *smf0_186 = buffer.data(smf0 + 186);
    const auto *smf0_188 = buffer.data(smf0 + 188);
    const auto *smf0_189 = buffer.data(smf0 + 189);
    const auto *smf0_196 = buffer.data(smf0 + 196);
    const auto *smf0_198 = buffer.data(smf0 + 198);
    const auto *smf0_199 = buffer.data(smf0 + 199);
    const auto *smf0_200 = buffer.data(smf0 + 200);
    const auto *smf0_203 = buffer.data(smf0 + 203);
    const auto *smf0_205 = buffer.data(smf0 + 205);
    const auto *smf0_206 = buffer.data(smf0 + 206);
    const auto *smf0_208 = buffer.data(smf0 + 208);
    const auto *smf0_209 = buffer.data(smf0 + 209);
    const auto *smf0_210 = buffer.data(smf0 + 210);
    const auto *smf0_213 = buffer.data(smf0 + 213);
    const auto *smf0_215 = buffer.data(smf0 + 215);
    const auto *smf0_216 = buffer.data(smf0 + 216);
    const auto *smf0_218 = buffer.data(smf0 + 218);
    const auto *smf0_219 = buffer.data(smf0 + 219);
    const auto *smf0_225 = buffer.data(smf0 + 225);
    const auto *smf0_229 = buffer.data(smf0 + 229);

    const auto *smf1_173 = buffer.data(smf1 + 173);
    const auto *smf1_175 = buffer.data(smf1 + 175);
    const auto *smf1_176 = buffer.data(smf1 + 176);
    const auto *smf1_178 = buffer.data(smf1 + 178);
    const auto *smf1_179 = buffer.data(smf1 + 179);
    const auto *smf1_180 = buffer.data(smf1 + 180);
    const auto *smf1_183 = buffer.data(smf1 + 183);
    const auto *smf1_185 = buffer.data(smf1 + 185);
    const auto *smf1_186 = buffer.data(smf1 + 186);
    const auto *smf1_188 = buffer.data(smf1 + 188);
    const auto *smf1_189 = buffer.data(smf1 + 189);
    const auto *smf1_196 = buffer.data(smf1 + 196);
    const auto *smf1_198 = buffer.data(smf1 + 198);
    const auto *smf1_199 = buffer.data(smf1 + 199);
    const auto *smf1_200 = buffer.data(smf1 + 200);
    const auto *smf1_203 = buffer.data(smf1 + 203);
    const auto *smf1_205 = buffer.data(smf1 + 205);
    const auto *smf1_206 = buffer.data(smf1 + 206);
    const auto *smf1_208 = buffer.data(smf1 + 208);
    const auto *smf1_209 = buffer.data(smf1 + 209);
    const auto *smf1_210 = buffer.data(smf1 + 210);
    const auto *smf1_213 = buffer.data(smf1 + 213);
    const auto *smf1_215 = buffer.data(smf1 + 215);
    const auto *smf1_216 = buffer.data(smf1 + 216);
    const auto *smf1_218 = buffer.data(smf1 + 218);
    const auto *smf1_219 = buffer.data(smf1 + 219);
    const auto *smf1_225 = buffer.data(smf1 + 225);
    const auto *smf1_229 = buffer.data(smf1 + 229);

    const auto *smg_255 = buffer.data(smg + 255);
    const auto *smg_257 = buffer.data(smg + 257);
    const auto *smg_258 = buffer.data(smg + 258);
    const auto *smg_260 = buffer.data(smg + 260);
    const auto *smg_261 = buffer.data(smg + 261);
    const auto *smg_264 = buffer.data(smg + 264);
    const auto *smg_265 = buffer.data(smg + 265);
    const auto *smg_266 = buffer.data(smg + 266);
    const auto *smg_267 = buffer.data(smg + 267);
    const auto *smg_268 = buffer.data(smg + 268);
    const auto *smg_269 = buffer.data(smg + 269);
    const auto *smg_270 = buffer.data(smg + 270);
    const auto *smg_272 = buffer.data(smg + 272);
    const auto *smg_273 = buffer.data(smg + 273);
    const auto *smg_275 = buffer.data(smg + 275);
    const auto *smg_276 = buffer.data(smg + 276);
    const auto *smg_279 = buffer.data(smg + 279);
    const auto *smg_280 = buffer.data(smg + 280);
    const auto *smg_281 = buffer.data(smg + 281);
    const auto *smg_282 = buffer.data(smg + 282);
    const auto *smg_283 = buffer.data(smg + 283);
    const auto *smg_284 = buffer.data(smg + 284);
    const auto *smg_285 = buffer.data(smg + 285);
    const auto *smg_287 = buffer.data(smg + 287);
    const auto *smg_288 = buffer.data(smg + 288);
    const auto *smg_290 = buffer.data(smg + 290);
    const auto *smg_295 = buffer.data(smg + 295);
    const auto *smg_296 = buffer.data(smg + 296);
    const auto *smg_297 = buffer.data(smg + 297);
    const auto *smg_298 = buffer.data(smg + 298);
    const auto *smg_299 = buffer.data(smg + 299);
    const auto *smg_300 = buffer.data(smg + 300);
    const auto *smg_302 = buffer.data(smg + 302);
    const auto *smg_303 = buffer.data(smg + 303);
    const auto *smg_305 = buffer.data(smg + 305);
    const auto *smg_306 = buffer.data(smg + 306);
    const auto *smg_309 = buffer.data(smg + 309);
    const auto *smg_310 = buffer.data(smg + 310);
    const auto *smg_311 = buffer.data(smg + 311);
    const auto *smg_312 = buffer.data(smg + 312);
    const auto *smg_313 = buffer.data(smg + 313);
    const auto *smg_314 = buffer.data(smg + 314);
    const auto *smg_315 = buffer.data(smg + 315);
    const auto *smg_317 = buffer.data(smg + 317);
    const auto *smg_318 = buffer.data(smg + 318);
    const auto *smg_320 = buffer.data(smg + 320);
    const auto *smg_321 = buffer.data(smg + 321);
    const auto *smg_324 = buffer.data(smg + 324);
    const auto *smg_325 = buffer.data(smg + 325);
    const auto *smg_326 = buffer.data(smg + 326);
    const auto *smg_327 = buffer.data(smg + 327);
    const auto *smg_328 = buffer.data(smg + 328);
    const auto *smg_329 = buffer.data(smg + 329);
    const auto *smg_330 = buffer.data(smg + 330);
    const auto *smg_332 = buffer.data(smg + 332);
    const auto *smg_333 = buffer.data(smg + 333);
    const auto *smg_335 = buffer.data(smg + 335);
    const auto *smg_339 = buffer.data(smg + 339);
    const auto *smg_340 = buffer.data(smg + 340);
    const auto *smg_341 = buffer.data(smg + 341);

#pragma omp simd aligned(t_358, t_359, t_360, t_361, pc_x, pc_y, pc_z, slg_165, slg_180, \
                         slg_182, slg_258, smf0_173, smf1_173, smg_255, smg_257, \
                         smg_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_358[k] = f_11 * slg_180[k]
                   + f_3 * pc_y[k] * smg_255[k];

        t_359[k] = f_10 * slg_165[k]
                   + f_3 * pc_z[k] * smg_255[k];

        t_360[k] = f_16 * slg_258[k]
                   + f_4 * smf0_173[k]
                   - f_5 * smf1_173[k]
                   + f_3 * pc_x[k] * smg_258[k];

        t_361[k] = f_11 * slg_182[k]
                   + f_3 * pc_y[k] * smg_257[k];
    }

#pragma omp simd aligned(t_362, t_363, t_364, pc_x, pc_z, slg_168, slg_260, slg_261, smf0_175, \
                         smf0_176, smf1_175, smf1_176, smg_258, smg_260, \
                         smg_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = f_16 * slg_260[k]
                   + f_4 * smf0_175[k]
                   - f_5 * smf1_175[k]
                   + f_3 * pc_x[k] * smg_260[k];

        t_363[k] = f_16 * slg_261[k]
                   + f_6 * smf0_176[k]
                   - f_7 * smf1_176[k]
                   + f_3 * pc_x[k] * smg_261[k];

        t_364[k] = f_10 * slg_168[k]
                   + f_3 * pc_z[k] * smg_258[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, pc_x, pc_y, slg_185, slg_264, slg_265, \
                         slg_266, smf0_179, smf1_179, smg_260, smg_264, smg_265, \
                         smg_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = f_11 * slg_185[k]
                   + f_3 * pc_y[k] * smg_260[k];

        t_366[k] = f_16 * slg_264[k]
                   + f_6 * smf0_179[k]
                   - f_7 * smf1_179[k]
                   + f_3 * pc_x[k] * smg_264[k];

        t_367[k] = f_16 * slg_265[k]
                   + f_3 * pc_x[k] * smg_265[k];

        t_368[k] = f_16 * slg_266[k]
                   + f_3 * pc_x[k] * smg_266[k];
    }

#pragma omp simd aligned(t_369, t_370, t_371, t_372, pc_x, pc_y, slg_190, slg_267, slg_268, \
                         slg_269, smf0_176, smf1_176, smg_265, smg_267, smg_268, \
                         smg_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_369[k] = f_16 * slg_267[k]
                   + f_3 * pc_x[k] * smg_267[k];

        t_370[k] = f_16 * slg_268[k]
                   + f_3 * pc_x[k] * smg_268[k];

        t_371[k] = f_16 * slg_269[k]
                   + f_3 * pc_x[k] * smg_269[k];

        t_372[k] = f_11 * slg_190[k]
                   + f_1 * smf0_176[k]
                   - f_2 * smf1_176[k]
                   + f_3 * pc_y[k] * smg_265[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, pc_y, pc_z, slg_175, slg_192, slg_193, smf0_178, \
                         smf0_179, smf1_178, smf1_179, smg_265, smg_267, \
                         smg_268 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = f_10 * slg_175[k]
                   + f_3 * pc_z[k] * smg_265[k];

        t_374[k] = f_11 * slg_192[k]
                   + f_4 * smf0_178[k]
                   - f_5 * smf1_178[k]
                   + f_3 * pc_y[k] * smg_267[k];

        t_375[k] = f_11 * slg_193[k]
                   + f_6 * smf0_179[k]
                   - f_7 * smf1_179[k]
                   + f_3 * pc_y[k] * smg_268[k];
    }

#pragma omp simd aligned(t_376, t_377, t_378, pc_x, pc_y, pc_z, slg_179, slg_194, slg_270, \
                         smf0_179, smf0_180, smf1_179, smf1_180, smg_269, \
                         smg_270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_376[k] = f_11 * slg_194[k]
                   + f_3 * pc_y[k] * smg_269[k];

        t_377[k] = f_10 * slg_179[k]
                   + f_1 * smf0_179[k]
                   - f_2 * smf1_179[k]
                   + f_3 * pc_z[k] * smg_269[k];

        t_378[k] = f_16 * slg_270[k]
                   + f_1 * smf0_180[k]
                   - f_2 * smf1_180[k]
                   + f_3 * pc_x[k] * smg_270[k];
    }

#pragma omp simd aligned(t_379, t_380, t_381, t_382, pc_x, pc_y, pc_z, slg_180, slg_195, \
                         slg_197, slg_273, smf0_183, smf1_183, smg_270, smg_272, \
                         smg_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_379[k] = f_10 * slg_195[k]
                   + f_3 * pc_y[k] * smg_270[k];

        t_380[k] = f_11 * slg_180[k]
                   + f_3 * pc_z[k] * smg_270[k];

        t_381[k] = f_16 * slg_273[k]
                   + f_4 * smf0_183[k]
                   - f_5 * smf1_183[k]
                   + f_3 * pc_x[k] * smg_273[k];

        t_382[k] = f_10 * slg_197[k]
                   + f_3 * pc_y[k] * smg_272[k];
    }

#pragma omp simd aligned(t_383, t_384, t_385, pc_x, pc_z, slg_183, slg_275, slg_276, smf0_185, \
                         smf0_186, smf1_185, smf1_186, smg_273, smg_275, \
                         smg_276 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_383[k] = f_16 * slg_275[k]
                   + f_4 * smf0_185[k]
                   - f_5 * smf1_185[k]
                   + f_3 * pc_x[k] * smg_275[k];

        t_384[k] = f_16 * slg_276[k]
                   + f_6 * smf0_186[k]
                   - f_7 * smf1_186[k]
                   + f_3 * pc_x[k] * smg_276[k];

        t_385[k] = f_11 * slg_183[k]
                   + f_3 * pc_z[k] * smg_273[k];
    }

#pragma omp simd aligned(t_386, t_387, t_388, t_389, pc_x, pc_y, slg_200, slg_279, slg_280, \
                         slg_281, smf0_189, smf1_189, smg_275, smg_279, smg_280, \
                         smg_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_386[k] = f_10 * slg_200[k]
                   + f_3 * pc_y[k] * smg_275[k];

        t_387[k] = f_16 * slg_279[k]
                   + f_6 * smf0_189[k]
                   - f_7 * smf1_189[k]
                   + f_3 * pc_x[k] * smg_279[k];

        t_388[k] = f_16 * slg_280[k]
                   + f_3 * pc_x[k] * smg_280[k];

        t_389[k] = f_16 * slg_281[k]
                   + f_3 * pc_x[k] * smg_281[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, pc_x, pc_y, slg_205, slg_282, slg_283, \
                         slg_284, smf0_186, smf1_186, smg_280, smg_282, smg_283, \
                         smg_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = f_16 * slg_282[k]
                   + f_3 * pc_x[k] * smg_282[k];

        t_391[k] = f_16 * slg_283[k]
                   + f_3 * pc_x[k] * smg_283[k];

        t_392[k] = f_16 * slg_284[k]
                   + f_3 * pc_x[k] * smg_284[k];

        t_393[k] = f_10 * slg_205[k]
                   + f_1 * smf0_186[k]
                   - f_2 * smf1_186[k]
                   + f_3 * pc_y[k] * smg_280[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, pc_y, pc_z, slg_190, slg_207, slg_208, smf0_188, \
                         smf0_189, smf1_188, smf1_189, smg_280, smg_282, \
                         smg_283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = f_11 * slg_190[k]
                   + f_3 * pc_z[k] * smg_280[k];

        t_395[k] = f_10 * slg_207[k]
                   + f_4 * smf0_188[k]
                   - f_5 * smf1_188[k]
                   + f_3 * pc_y[k] * smg_282[k];

        t_396[k] = f_10 * slg_208[k]
                   + f_6 * smf0_189[k]
                   - f_7 * smf1_189[k]
                   + f_3 * pc_y[k] * smg_283[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, t_400, pb_y, pc_y, pc_z, slh0_294, slg_194, \
                         slg_209, slg_210, slh1_294, smf0_189, smf1_189, smg_284, \
                         smg_285 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = f_10 * slg_209[k]
                   + f_3 * pc_y[k] * smg_284[k];

        t_398[k] = f_11 * slg_194[k]
                   + f_1 * smf0_189[k]
                   - f_2 * smf1_189[k]
                   + f_3 * pc_z[k] * smg_284[k];

        t_399[k] = pb_y[k] * slh0_294[k]
                   - f_8 * pc_y[k] * slh1_294[k];

        t_400[k] = f_9 * slg_210[k]
                   + f_3 * pc_y[k] * smg_285[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, t_404, pb_y, pc_y, pc_z, slh0_297, slh0_299, \
                         slg_195, slg_211, slg_212, slh1_297, slh1_299, smg_285, \
                         smg_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_16 * slg_195[k]
                   + f_3 * pc_z[k] * smg_285[k];

        t_402[k] = pb_y[k] * slh0_297[k]
                   + f_10 * slg_211[k]
                   - f_8 * pc_y[k] * slh1_297[k];

        t_403[k] = f_9 * slg_212[k]
                   + f_3 * pc_y[k] * smg_287[k];

        t_404[k] = pb_y[k] * slh0_299[k]
                   - f_8 * pc_y[k] * slh1_299[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, pb_y, pc_y, pc_z, slh0_300, slh0_303, \
                         slg_198, slg_213, slg_215, slh1_300, slh1_303, smg_288, \
                         smg_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = pb_y[k] * slh0_300[k]
                   + f_11 * slg_213[k]
                   - f_8 * pc_y[k] * slh1_300[k];

        t_406[k] = f_16 * slg_198[k]
                   + f_3 * pc_z[k] * smg_288[k];

        t_407[k] = f_9 * slg_215[k]
                   + f_3 * pc_y[k] * smg_290[k];

        t_408[k] = pb_y[k] * slh0_303[k]
                   - f_8 * pc_y[k] * slh1_303[k];
    }

#pragma omp simd aligned(t_409, t_410, t_411, t_412, t_413, pc_x, slg_295, slg_296, slg_297, \
                         slg_298, slg_299, smg_295, smg_296, smg_297, smg_298, \
                         smg_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_409[k] = f_16 * slg_295[k]
                   + f_3 * pc_x[k] * smg_295[k];

        t_410[k] = f_16 * slg_296[k]
                   + f_3 * pc_x[k] * smg_296[k];

        t_411[k] = f_16 * slg_297[k]
                   + f_3 * pc_x[k] * smg_297[k];

        t_412[k] = f_16 * slg_298[k]
                   + f_3 * pc_x[k] * smg_298[k];

        t_413[k] = f_16 * slg_299[k]
                   + f_3 * pc_x[k] * smg_299[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, pc_y, pc_z, slg_205, slg_220, slg_222, smf0_196, \
                         smf0_198, smf1_196, smf1_198, smg_295, \
                         smg_297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_9 * slg_220[k]
                   + f_1 * smf0_196[k]
                   - f_2 * smf1_196[k]
                   + f_3 * pc_y[k] * smg_295[k];

        t_415[k] = f_16 * slg_205[k]
                   + f_3 * pc_z[k] * smg_295[k];

        t_416[k] = f_9 * slg_222[k]
                   + f_4 * smf0_198[k]
                   - f_5 * smf1_198[k]
                   + f_3 * pc_y[k] * smg_297[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, pb_y, pc_y, slh0_314, slg_223, slg_224, \
                         slh1_314, smf0_199, smf1_199, smg_298, \
                         smg_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = f_9 * slg_223[k]
                   + f_6 * smf0_199[k]
                   - f_7 * smf1_199[k]
                   + f_3 * pc_y[k] * smg_298[k];

        t_418[k] = f_9 * slg_224[k]
                   + f_3 * pc_y[k] * smg_299[k];

        t_419[k] = pb_y[k] * slh0_314[k]
                   - f_8 * pc_y[k] * slh1_314[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, pc_x, pc_y, pc_z, slg_210, slg_300, \
                         slg_303, smf0_200, smf0_203, smf1_200, smf1_203, smg_300, \
                         smg_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_16 * slg_300[k]
                   + f_1 * smf0_200[k]
                   - f_2 * smf1_200[k]
                   + f_3 * pc_x[k] * smg_300[k];

        t_421[k] = f_3 * pc_y[k] * smg_300[k];

        t_422[k] = f_15 * slg_210[k]
                   + f_3 * pc_z[k] * smg_300[k];

        t_423[k] = f_16 * slg_303[k]
                   + f_4 * smf0_203[k]
                   - f_5 * smf1_203[k]
                   + f_3 * pc_x[k] * smg_303[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, pc_x, pc_y, slg_305, slg_306, smf0_205, \
                         smf0_206, smf1_205, smf1_206, smg_302, smg_305, \
                         smg_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = f_3 * pc_y[k] * smg_302[k];

        t_425[k] = f_16 * slg_305[k]
                   + f_4 * smf0_205[k]
                   - f_5 * smf1_205[k]
                   + f_3 * pc_x[k] * smg_305[k];

        t_426[k] = f_16 * slg_306[k]
                   + f_6 * smf0_206[k]
                   - f_7 * smf1_206[k]
                   + f_3 * pc_x[k] * smg_306[k];
    }

#pragma omp simd aligned(t_427, t_428, t_429, t_430, pc_x, pc_y, pc_z, slg_213, slg_309, \
                         slg_310, smf0_209, smf1_209, smg_303, smg_305, smg_309, \
                         smg_310 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_427[k] = f_15 * slg_213[k]
                   + f_3 * pc_z[k] * smg_303[k];

        t_428[k] = f_3 * pc_y[k] * smg_305[k];

        t_429[k] = f_16 * slg_309[k]
                   + f_6 * smf0_209[k]
                   - f_7 * smf1_209[k]
                   + f_3 * pc_x[k] * smg_309[k];

        t_430[k] = f_16 * slg_310[k]
                   + f_3 * pc_x[k] * smg_310[k];
    }

#pragma omp simd aligned(t_431, t_432, t_433, t_434, pc_x, slg_311, slg_312, slg_313, slg_314, \
                         smg_311, smg_312, smg_313, smg_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_431[k] = f_16 * slg_311[k]
                   + f_3 * pc_x[k] * smg_311[k];

        t_432[k] = f_16 * slg_312[k]
                   + f_3 * pc_x[k] * smg_312[k];

        t_433[k] = f_16 * slg_313[k]
                   + f_3 * pc_x[k] * smg_313[k];

        t_434[k] = f_16 * slg_314[k]
                   + f_3 * pc_x[k] * smg_314[k];
    }

#pragma omp simd aligned(t_435, t_436, t_437, t_438, pc_y, pc_z, slg_220, smf0_206, smf0_208, \
                         smf0_209, smf1_206, smf1_208, smf1_209, smg_310, smg_312, \
                         smg_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_435[k] = f_1 * smf0_206[k]
                   - f_2 * smf1_206[k]
                   + f_3 * pc_y[k] * smg_310[k];

        t_436[k] = f_15 * slg_220[k]
                   + f_3 * pc_z[k] * smg_310[k];

        t_437[k] = f_4 * smf0_208[k]
                   - f_5 * smf1_208[k]
                   + f_3 * pc_y[k] * smg_312[k];

        t_438[k] = f_6 * smf0_209[k]
                   - f_7 * smf1_209[k]
                   + f_3 * pc_y[k] * smg_313[k];
    }

#pragma omp simd aligned(t_439, t_440, t_441, t_442, pc_x, pc_y, pc_z, slg_224, slg_225, \
                         slg_315, smf0_209, smf0_210, smf1_209, smf1_210, smg_314, \
                         smg_315 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_439[k] = f_3 * pc_y[k] * smg_314[k];

        t_440[k] = f_15 * slg_224[k]
                   + f_1 * smf0_209[k]
                   - f_2 * smf1_209[k]
                   + f_3 * pc_z[k] * smg_314[k];

        t_441[k] = f_11 * slg_315[k]
                   + f_1 * smf0_210[k]
                   - f_2 * smf1_210[k]
                   + f_3 * pc_x[k] * smg_315[k];

        t_442[k] = f_14 * slg_225[k]
                   + f_3 * pc_y[k] * smg_315[k];
    }

#pragma omp simd aligned(t_443, t_444, t_445, pc_x, pc_y, pc_z, slg_227, slg_318, smf0_213, \
                         smf1_213, smg_315, smg_317, smg_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_443[k] = f_3 * pc_z[k] * smg_315[k];

        t_444[k] = f_11 * slg_318[k]
                   + f_4 * smf0_213[k]
                   - f_5 * smf1_213[k]
                   + f_3 * pc_x[k] * smg_318[k];

        t_445[k] = f_14 * slg_227[k]
                   + f_3 * pc_y[k] * smg_317[k];
    }

#pragma omp simd aligned(t_446, t_447, t_448, pc_x, pc_z, slg_320, slg_321, smf0_215, \
                         smf0_216, smf1_215, smf1_216, smg_318, smg_320, \
                         smg_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_446[k] = f_11 * slg_320[k]
                   + f_4 * smf0_215[k]
                   - f_5 * smf1_215[k]
                   + f_3 * pc_x[k] * smg_320[k];

        t_447[k] = f_11 * slg_321[k]
                   + f_6 * smf0_216[k]
                   - f_7 * smf1_216[k]
                   + f_3 * pc_x[k] * smg_321[k];

        t_448[k] = f_3 * pc_z[k] * smg_318[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, t_452, pc_x, pc_y, slg_230, slg_324, slg_325, \
                         slg_326, smf0_219, smf1_219, smg_320, smg_324, smg_325, \
                         smg_326 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_14 * slg_230[k]
                   + f_3 * pc_y[k] * smg_320[k];

        t_450[k] = f_11 * slg_324[k]
                   + f_6 * smf0_219[k]
                   - f_7 * smf1_219[k]
                   + f_3 * pc_x[k] * smg_324[k];

        t_451[k] = f_11 * slg_325[k]
                   + f_3 * pc_x[k] * smg_325[k];

        t_452[k] = f_11 * slg_326[k]
                   + f_3 * pc_x[k] * smg_326[k];
    }

#pragma omp simd aligned(t_453, t_454, t_455, t_456, pc_x, pc_y, slg_235, slg_327, slg_328, \
                         slg_329, smf0_216, smf1_216, smg_325, smg_327, smg_328, \
                         smg_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_453[k] = f_11 * slg_327[k]
                   + f_3 * pc_x[k] * smg_327[k];

        t_454[k] = f_11 * slg_328[k]
                   + f_3 * pc_x[k] * smg_328[k];

        t_455[k] = f_11 * slg_329[k]
                   + f_3 * pc_x[k] * smg_329[k];

        t_456[k] = f_14 * slg_235[k]
                   + f_1 * smf0_216[k]
                   - f_2 * smf1_216[k]
                   + f_3 * pc_y[k] * smg_325[k];
    }

#pragma omp simd aligned(t_457, t_458, t_459, pc_y, pc_z, slg_237, slg_238, smf0_218, \
                         smf0_219, smf1_218, smf1_219, smg_325, smg_327, \
                         smg_328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_457[k] = f_3 * pc_z[k] * smg_325[k];

        t_458[k] = f_14 * slg_237[k]
                   + f_4 * smf0_218[k]
                   - f_5 * smf1_218[k]
                   + f_3 * pc_y[k] * smg_327[k];

        t_459[k] = f_14 * slg_238[k]
                   + f_6 * smf0_219[k]
                   - f_7 * smf1_219[k]
                   + f_3 * pc_y[k] * smg_328[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, t_463, pb_z, pc_y, pc_z, slh0_315, slg_239, \
                         slg_240, slh1_315, smf0_219, smf1_219, smg_329, \
                         smg_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = f_14 * slg_239[k]
                   + f_3 * pc_y[k] * smg_329[k];

        t_461[k] = f_1 * smf0_219[k]
                   - f_2 * smf1_219[k]
                   + f_3 * pc_z[k] * smg_329[k];

        t_462[k] = pb_z[k] * slh0_315[k]
                   - f_8 * pc_z[k] * slh1_315[k];

        t_463[k] = f_15 * slg_240[k]
                   + f_3 * pc_y[k] * smg_330[k];
    }

#pragma omp simd aligned(t_464, t_465, t_466, pb_z, pc_y, pc_z, slh0_318, slg_225, slg_242, \
                         slh1_318, smg_330, smg_332 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_464[k] = f_9 * slg_225[k]
                   + f_3 * pc_z[k] * smg_330[k];

        t_465[k] = pb_z[k] * slh0_318[k]
                   - f_8 * pc_z[k] * slh1_318[k];

        t_466[k] = f_15 * slg_242[k]
                   + f_3 * pc_y[k] * smg_332[k];
    }

#pragma omp simd aligned(t_467, t_468, t_469, pb_z, pc_x, pc_z, slh0_321, slg_228, slg_335, \
                         slh1_321, smf0_225, smf1_225, smg_333, \
                         smg_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_467[k] = f_11 * slg_335[k]
                   + f_4 * smf0_225[k]
                   - f_5 * smf1_225[k]
                   + f_3 * pc_x[k] * smg_335[k];

        t_468[k] = pb_z[k] * slh0_321[k]
                   - f_8 * pc_z[k] * slh1_321[k];

        t_469[k] = f_9 * slg_228[k]
                   + f_3 * pc_z[k] * smg_333[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, pc_x, pc_y, slg_245, slg_339, slg_340, \
                         slg_341, smf0_229, smf1_229, smg_335, smg_339, smg_340, \
                         smg_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = f_15 * slg_245[k]
                   + f_3 * pc_y[k] * smg_335[k];

        t_471[k] = f_11 * slg_339[k]
                   + f_6 * smf0_229[k]
                   - f_7 * smf1_229[k]
                   + f_3 * pc_x[k] * smg_339[k];

        t_472[k] = f_11 * slg_340[k]
                   + f_3 * pc_x[k] * smg_340[k];

        t_473[k] = f_11 * slg_341[k]
                   + f_3 * pc_x[k] * smg_341[k];
    }
}

static auto
compute_prim_smh_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t slh0,
                                                          const size_t slg, const size_t slh1,
                                                          const size_t smf0, const size_t smf1,
                                                          const size_t smg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.0 / gamma;
    const auto f_5 = p / (gamma * q);
    const auto f_6 = 0.5 / gamma;
    const auto f_7 = 0.5 * p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_14 = 3.0 / q;
    const auto f_15 = 2.5 / q;
    const auto f_16 = 2.0 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *slh0_330 = buffer.data(slh0 + 330);
    const auto *slh0_420 = buffer.data(slh0 + 420);
    const auto *slh0_423 = buffer.data(slh0 + 423);
    const auto *slh0_425 = buffer.data(slh0 + 425);
    const auto *slh0_426 = buffer.data(slh0 + 426);
    const auto *slh0_429 = buffer.data(slh0 + 429);
    const auto *slh0_440 = buffer.data(slh0 + 440);

    const auto *slg_235 = buffer.data(slg + 235);
    const auto *slg_239 = buffer.data(slg + 239);
    const auto *slg_240 = buffer.data(slg + 240);
    const auto *slg_243 = buffer.data(slg + 243);
    const auto *slg_250 = buffer.data(slg + 250);
    const auto *slg_252 = buffer.data(slg + 252);
    const auto *slg_253 = buffer.data(slg + 253);
    const auto *slg_254 = buffer.data(slg + 254);
    const auto *slg_255 = buffer.data(slg + 255);
    const auto *slg_257 = buffer.data(slg + 257);
    const auto *slg_258 = buffer.data(slg + 258);
    const auto *slg_260 = buffer.data(slg + 260);
    const auto *slg_265 = buffer.data(slg + 265);
    const auto *slg_267 = buffer.data(slg + 267);
    const auto *slg_268 = buffer.data(slg + 268);
    const auto *slg_269 = buffer.data(slg + 269);
    const auto *slg_270 = buffer.data(slg + 270);
    const auto *slg_272 = buffer.data(slg + 272);
    const auto *slg_273 = buffer.data(slg + 273);
    const auto *slg_275 = buffer.data(slg + 275);
    const auto *slg_280 = buffer.data(slg + 280);
    const auto *slg_282 = buffer.data(slg + 282);
    const auto *slg_283 = buffer.data(slg + 283);
    const auto *slg_284 = buffer.data(slg + 284);
    const auto *slg_285 = buffer.data(slg + 285);
    const auto *slg_287 = buffer.data(slg + 287);
    const auto *slg_288 = buffer.data(slg + 288);
    const auto *slg_290 = buffer.data(slg + 290);
    const auto *slg_295 = buffer.data(slg + 295);
    const auto *slg_297 = buffer.data(slg + 297);
    const auto *slg_298 = buffer.data(slg + 298);
    const auto *slg_299 = buffer.data(slg + 299);
    const auto *slg_300 = buffer.data(slg + 300);
    const auto *slg_301 = buffer.data(slg + 301);
    const auto *slg_302 = buffer.data(slg + 302);
    const auto *slg_303 = buffer.data(slg + 303);
    const auto *slg_305 = buffer.data(slg + 305);
    const auto *slg_310 = buffer.data(slg + 310);
    const auto *slg_312 = buffer.data(slg + 312);
    const auto *slg_313 = buffer.data(slg + 313);
    const auto *slg_314 = buffer.data(slg + 314);
    const auto *slg_342 = buffer.data(slg + 342);
    const auto *slg_343 = buffer.data(slg + 343);
    const auto *slg_344 = buffer.data(slg + 344);
    const auto *slg_345 = buffer.data(slg + 345);
    const auto *slg_348 = buffer.data(slg + 348);
    const auto *slg_350 = buffer.data(slg + 350);
    const auto *slg_351 = buffer.data(slg + 351);
    const auto *slg_354 = buffer.data(slg + 354);
    const auto *slg_355 = buffer.data(slg + 355);
    const auto *slg_356 = buffer.data(slg + 356);
    const auto *slg_357 = buffer.data(slg + 357);
    const auto *slg_358 = buffer.data(slg + 358);
    const auto *slg_359 = buffer.data(slg + 359);
    const auto *slg_360 = buffer.data(slg + 360);
    const auto *slg_363 = buffer.data(slg + 363);
    const auto *slg_365 = buffer.data(slg + 365);
    const auto *slg_366 = buffer.data(slg + 366);
    const auto *slg_369 = buffer.data(slg + 369);
    const auto *slg_370 = buffer.data(slg + 370);
    const auto *slg_371 = buffer.data(slg + 371);
    const auto *slg_372 = buffer.data(slg + 372);
    const auto *slg_373 = buffer.data(slg + 373);
    const auto *slg_374 = buffer.data(slg + 374);
    const auto *slg_375 = buffer.data(slg + 375);
    const auto *slg_378 = buffer.data(slg + 378);
    const auto *slg_380 = buffer.data(slg + 380);
    const auto *slg_381 = buffer.data(slg + 381);
    const auto *slg_384 = buffer.data(slg + 384);
    const auto *slg_385 = buffer.data(slg + 385);
    const auto *slg_386 = buffer.data(slg + 386);
    const auto *slg_387 = buffer.data(slg + 387);
    const auto *slg_388 = buffer.data(slg + 388);
    const auto *slg_389 = buffer.data(slg + 389);
    const auto *slg_400 = buffer.data(slg + 400);
    const auto *slg_401 = buffer.data(slg + 401);
    const auto *slg_402 = buffer.data(slg + 402);
    const auto *slg_403 = buffer.data(slg + 403);
    const auto *slg_404 = buffer.data(slg + 404);
    const auto *slg_405 = buffer.data(slg + 405);
    const auto *slg_408 = buffer.data(slg + 408);
    const auto *slg_410 = buffer.data(slg + 410);
    const auto *slg_411 = buffer.data(slg + 411);
    const auto *slg_414 = buffer.data(slg + 414);
    const auto *slg_415 = buffer.data(slg + 415);
    const auto *slg_416 = buffer.data(slg + 416);
    const auto *slg_417 = buffer.data(slg + 417);
    const auto *slg_418 = buffer.data(slg + 418);
    const auto *slg_419 = buffer.data(slg + 419);

    const auto *slh1_330 = buffer.data(slh1 + 330);
    const auto *slh1_420 = buffer.data(slh1 + 420);
    const auto *slh1_423 = buffer.data(slh1 + 423);
    const auto *slh1_425 = buffer.data(slh1 + 425);
    const auto *slh1_426 = buffer.data(slh1 + 426);
    const auto *slh1_429 = buffer.data(slh1 + 429);
    const auto *slh1_440 = buffer.data(slh1 + 440);

    const auto *smf0_228 = buffer.data(smf0 + 228);
    const auto *smf0_229 = buffer.data(smf0 + 229);
    const auto *smf0_230 = buffer.data(smf0 + 230);
    const auto *smf0_233 = buffer.data(smf0 + 233);
    const auto *smf0_235 = buffer.data(smf0 + 235);
    const auto *smf0_236 = buffer.data(smf0 + 236);
    const auto *smf0_238 = buffer.data(smf0 + 238);
    const auto *smf0_239 = buffer.data(smf0 + 239);
    const auto *smf0_240 = buffer.data(smf0 + 240);
    const auto *smf0_243 = buffer.data(smf0 + 243);
    const auto *smf0_245 = buffer.data(smf0 + 245);
    const auto *smf0_246 = buffer.data(smf0 + 246);
    const auto *smf0_248 = buffer.data(smf0 + 248);
    const auto *smf0_249 = buffer.data(smf0 + 249);
    const auto *smf0_250 = buffer.data(smf0 + 250);
    const auto *smf0_253 = buffer.data(smf0 + 253);
    const auto *smf0_255 = buffer.data(smf0 + 255);
    const auto *smf0_256 = buffer.data(smf0 + 256);
    const auto *smf0_258 = buffer.data(smf0 + 258);
    const auto *smf0_259 = buffer.data(smf0 + 259);
    const auto *smf0_266 = buffer.data(smf0 + 266);
    const auto *smf0_268 = buffer.data(smf0 + 268);
    const auto *smf0_269 = buffer.data(smf0 + 269);
    const auto *smf0_270 = buffer.data(smf0 + 270);
    const auto *smf0_273 = buffer.data(smf0 + 273);
    const auto *smf0_275 = buffer.data(smf0 + 275);
    const auto *smf0_276 = buffer.data(smf0 + 276);
    const auto *smf0_278 = buffer.data(smf0 + 278);
    const auto *smf0_279 = buffer.data(smf0 + 279);

    const auto *smf1_228 = buffer.data(smf1 + 228);
    const auto *smf1_229 = buffer.data(smf1 + 229);
    const auto *smf1_230 = buffer.data(smf1 + 230);
    const auto *smf1_233 = buffer.data(smf1 + 233);
    const auto *smf1_235 = buffer.data(smf1 + 235);
    const auto *smf1_236 = buffer.data(smf1 + 236);
    const auto *smf1_238 = buffer.data(smf1 + 238);
    const auto *smf1_239 = buffer.data(smf1 + 239);
    const auto *smf1_240 = buffer.data(smf1 + 240);
    const auto *smf1_243 = buffer.data(smf1 + 243);
    const auto *smf1_245 = buffer.data(smf1 + 245);
    const auto *smf1_246 = buffer.data(smf1 + 246);
    const auto *smf1_248 = buffer.data(smf1 + 248);
    const auto *smf1_249 = buffer.data(smf1 + 249);
    const auto *smf1_250 = buffer.data(smf1 + 250);
    const auto *smf1_253 = buffer.data(smf1 + 253);
    const auto *smf1_255 = buffer.data(smf1 + 255);
    const auto *smf1_256 = buffer.data(smf1 + 256);
    const auto *smf1_258 = buffer.data(smf1 + 258);
    const auto *smf1_259 = buffer.data(smf1 + 259);
    const auto *smf1_266 = buffer.data(smf1 + 266);
    const auto *smf1_268 = buffer.data(smf1 + 268);
    const auto *smf1_269 = buffer.data(smf1 + 269);
    const auto *smf1_270 = buffer.data(smf1 + 270);
    const auto *smf1_273 = buffer.data(smf1 + 273);
    const auto *smf1_275 = buffer.data(smf1 + 275);
    const auto *smf1_276 = buffer.data(smf1 + 276);
    const auto *smf1_278 = buffer.data(smf1 + 278);
    const auto *smf1_279 = buffer.data(smf1 + 279);

    const auto *smg_340 = buffer.data(smg + 340);
    const auto *smg_342 = buffer.data(smg + 342);
    const auto *smg_343 = buffer.data(smg + 343);
    const auto *smg_344 = buffer.data(smg + 344);
    const auto *smg_345 = buffer.data(smg + 345);
    const auto *smg_347 = buffer.data(smg + 347);
    const auto *smg_348 = buffer.data(smg + 348);
    const auto *smg_350 = buffer.data(smg + 350);
    const auto *smg_351 = buffer.data(smg + 351);
    const auto *smg_354 = buffer.data(smg + 354);
    const auto *smg_355 = buffer.data(smg + 355);
    const auto *smg_356 = buffer.data(smg + 356);
    const auto *smg_357 = buffer.data(smg + 357);
    const auto *smg_358 = buffer.data(smg + 358);
    const auto *smg_359 = buffer.data(smg + 359);
    const auto *smg_360 = buffer.data(smg + 360);
    const auto *smg_362 = buffer.data(smg + 362);
    const auto *smg_363 = buffer.data(smg + 363);
    const auto *smg_365 = buffer.data(smg + 365);
    const auto *smg_366 = buffer.data(smg + 366);
    const auto *smg_369 = buffer.data(smg + 369);
    const auto *smg_370 = buffer.data(smg + 370);
    const auto *smg_371 = buffer.data(smg + 371);
    const auto *smg_372 = buffer.data(smg + 372);
    const auto *smg_373 = buffer.data(smg + 373);
    const auto *smg_374 = buffer.data(smg + 374);
    const auto *smg_375 = buffer.data(smg + 375);
    const auto *smg_377 = buffer.data(smg + 377);
    const auto *smg_378 = buffer.data(smg + 378);
    const auto *smg_380 = buffer.data(smg + 380);
    const auto *smg_381 = buffer.data(smg + 381);
    const auto *smg_384 = buffer.data(smg + 384);
    const auto *smg_385 = buffer.data(smg + 385);
    const auto *smg_386 = buffer.data(smg + 386);
    const auto *smg_387 = buffer.data(smg + 387);
    const auto *smg_388 = buffer.data(smg + 388);
    const auto *smg_389 = buffer.data(smg + 389);
    const auto *smg_390 = buffer.data(smg + 390);
    const auto *smg_392 = buffer.data(smg + 392);
    const auto *smg_393 = buffer.data(smg + 393);
    const auto *smg_395 = buffer.data(smg + 395);
    const auto *smg_400 = buffer.data(smg + 400);
    const auto *smg_401 = buffer.data(smg + 401);
    const auto *smg_402 = buffer.data(smg + 402);
    const auto *smg_403 = buffer.data(smg + 403);
    const auto *smg_404 = buffer.data(smg + 404);
    const auto *smg_405 = buffer.data(smg + 405);
    const auto *smg_407 = buffer.data(smg + 407);
    const auto *smg_408 = buffer.data(smg + 408);
    const auto *smg_410 = buffer.data(smg + 410);
    const auto *smg_411 = buffer.data(smg + 411);
    const auto *smg_414 = buffer.data(smg + 414);
    const auto *smg_415 = buffer.data(smg + 415);
    const auto *smg_416 = buffer.data(smg + 416);
    const auto *smg_417 = buffer.data(smg + 417);
    const auto *smg_418 = buffer.data(smg + 418);
    const auto *smg_419 = buffer.data(smg + 419);

#pragma omp simd aligned(t_474, t_475, t_476, t_477, pb_z, pc_x, pc_z, slh0_330, slg_342, \
                         slg_343, slg_344, slh1_330, smg_342, smg_343, \
                         smg_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = f_11 * slg_342[k]
                   + f_3 * pc_x[k] * smg_342[k];

        t_475[k] = f_11 * slg_343[k]
                   + f_3 * pc_x[k] * smg_343[k];

        t_476[k] = f_11 * slg_344[k]
                   + f_3 * pc_x[k] * smg_344[k];

        t_477[k] = pb_z[k] * slh0_330[k]
                   - f_8 * pc_z[k] * slh1_330[k];
    }

#pragma omp simd aligned(t_478, t_479, t_480, pc_y, pc_z, slg_235, slg_252, slg_253, smf0_228, \
                         smf0_229, smf1_228, smf1_229, smg_340, smg_342, \
                         smg_343 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_478[k] = f_9 * slg_235[k]
                   + f_3 * pc_z[k] * smg_340[k];

        t_479[k] = f_15 * slg_252[k]
                   + f_4 * smf0_228[k]
                   - f_5 * smf1_228[k]
                   + f_3 * pc_y[k] * smg_342[k];

        t_480[k] = f_15 * slg_253[k]
                   + f_6 * smf0_229[k]
                   - f_7 * smf1_229[k]
                   + f_3 * pc_y[k] * smg_343[k];
    }

#pragma omp simd aligned(t_481, t_482, t_483, pc_x, pc_y, pc_z, slg_239, slg_254, slg_345, \
                         smf0_229, smf0_230, smf1_229, smf1_230, smg_344, \
                         smg_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_481[k] = f_15 * slg_254[k]
                   + f_3 * pc_y[k] * smg_344[k];

        t_482[k] = f_9 * slg_239[k]
                   + f_1 * smf0_229[k]
                   - f_2 * smf1_229[k]
                   + f_3 * pc_z[k] * smg_344[k];

        t_483[k] = f_11 * slg_345[k]
                   + f_1 * smf0_230[k]
                   - f_2 * smf1_230[k]
                   + f_3 * pc_x[k] * smg_345[k];
    }

#pragma omp simd aligned(t_484, t_485, t_486, t_487, pc_x, pc_y, pc_z, slg_240, slg_255, \
                         slg_257, slg_348, smf0_233, smf1_233, smg_345, smg_347, \
                         smg_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_484[k] = f_16 * slg_255[k]
                   + f_3 * pc_y[k] * smg_345[k];

        t_485[k] = f_10 * slg_240[k]
                   + f_3 * pc_z[k] * smg_345[k];

        t_486[k] = f_11 * slg_348[k]
                   + f_4 * smf0_233[k]
                   - f_5 * smf1_233[k]
                   + f_3 * pc_x[k] * smg_348[k];

        t_487[k] = f_16 * slg_257[k]
                   + f_3 * pc_y[k] * smg_347[k];
    }

#pragma omp simd aligned(t_488, t_489, t_490, pc_x, pc_z, slg_243, slg_350, slg_351, smf0_235, \
                         smf0_236, smf1_235, smf1_236, smg_348, smg_350, \
                         smg_351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_488[k] = f_11 * slg_350[k]
                   + f_4 * smf0_235[k]
                   - f_5 * smf1_235[k]
                   + f_3 * pc_x[k] * smg_350[k];

        t_489[k] = f_11 * slg_351[k]
                   + f_6 * smf0_236[k]
                   - f_7 * smf1_236[k]
                   + f_3 * pc_x[k] * smg_351[k];

        t_490[k] = f_10 * slg_243[k]
                   + f_3 * pc_z[k] * smg_348[k];
    }

#pragma omp simd aligned(t_491, t_492, t_493, t_494, pc_x, pc_y, slg_260, slg_354, slg_355, \
                         slg_356, smf0_239, smf1_239, smg_350, smg_354, smg_355, \
                         smg_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_491[k] = f_16 * slg_260[k]
                   + f_3 * pc_y[k] * smg_350[k];

        t_492[k] = f_11 * slg_354[k]
                   + f_6 * smf0_239[k]
                   - f_7 * smf1_239[k]
                   + f_3 * pc_x[k] * smg_354[k];

        t_493[k] = f_11 * slg_355[k]
                   + f_3 * pc_x[k] * smg_355[k];

        t_494[k] = f_11 * slg_356[k]
                   + f_3 * pc_x[k] * smg_356[k];
    }

#pragma omp simd aligned(t_495, t_496, t_497, t_498, pc_x, pc_y, slg_265, slg_357, slg_358, \
                         slg_359, smf0_236, smf1_236, smg_355, smg_357, smg_358, \
                         smg_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_495[k] = f_11 * slg_357[k]
                   + f_3 * pc_x[k] * smg_357[k];

        t_496[k] = f_11 * slg_358[k]
                   + f_3 * pc_x[k] * smg_358[k];

        t_497[k] = f_11 * slg_359[k]
                   + f_3 * pc_x[k] * smg_359[k];

        t_498[k] = f_16 * slg_265[k]
                   + f_1 * smf0_236[k]
                   - f_2 * smf1_236[k]
                   + f_3 * pc_y[k] * smg_355[k];
    }

#pragma omp simd aligned(t_499, t_500, t_501, pc_y, pc_z, slg_250, slg_267, slg_268, smf0_238, \
                         smf0_239, smf1_238, smf1_239, smg_355, smg_357, \
                         smg_358 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_499[k] = f_10 * slg_250[k]
                   + f_3 * pc_z[k] * smg_355[k];

        t_500[k] = f_16 * slg_267[k]
                   + f_4 * smf0_238[k]
                   - f_5 * smf1_238[k]
                   + f_3 * pc_y[k] * smg_357[k];

        t_501[k] = f_16 * slg_268[k]
                   + f_6 * smf0_239[k]
                   - f_7 * smf1_239[k]
                   + f_3 * pc_y[k] * smg_358[k];
    }

#pragma omp simd aligned(t_502, t_503, t_504, pc_x, pc_y, pc_z, slg_254, slg_269, slg_360, \
                         smf0_239, smf0_240, smf1_239, smf1_240, smg_359, \
                         smg_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_502[k] = f_16 * slg_269[k]
                   + f_3 * pc_y[k] * smg_359[k];

        t_503[k] = f_10 * slg_254[k]
                   + f_1 * smf0_239[k]
                   - f_2 * smf1_239[k]
                   + f_3 * pc_z[k] * smg_359[k];

        t_504[k] = f_11 * slg_360[k]
                   + f_1 * smf0_240[k]
                   - f_2 * smf1_240[k]
                   + f_3 * pc_x[k] * smg_360[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, pc_x, pc_y, pc_z, slg_255, slg_270, \
                         slg_272, slg_363, smf0_243, smf1_243, smg_360, smg_362, \
                         smg_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = f_11 * slg_270[k]
                   + f_3 * pc_y[k] * smg_360[k];

        t_506[k] = f_11 * slg_255[k]
                   + f_3 * pc_z[k] * smg_360[k];

        t_507[k] = f_11 * slg_363[k]
                   + f_4 * smf0_243[k]
                   - f_5 * smf1_243[k]
                   + f_3 * pc_x[k] * smg_363[k];

        t_508[k] = f_11 * slg_272[k]
                   + f_3 * pc_y[k] * smg_362[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, pc_x, pc_z, slg_258, slg_365, slg_366, smf0_245, \
                         smf0_246, smf1_245, smf1_246, smg_363, smg_365, \
                         smg_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = f_11 * slg_365[k]
                   + f_4 * smf0_245[k]
                   - f_5 * smf1_245[k]
                   + f_3 * pc_x[k] * smg_365[k];

        t_510[k] = f_11 * slg_366[k]
                   + f_6 * smf0_246[k]
                   - f_7 * smf1_246[k]
                   + f_3 * pc_x[k] * smg_366[k];

        t_511[k] = f_11 * slg_258[k]
                   + f_3 * pc_z[k] * smg_363[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, t_515, pc_x, pc_y, slg_275, slg_369, slg_370, \
                         slg_371, smf0_249, smf1_249, smg_365, smg_369, smg_370, \
                         smg_371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = f_11 * slg_275[k]
                   + f_3 * pc_y[k] * smg_365[k];

        t_513[k] = f_11 * slg_369[k]
                   + f_6 * smf0_249[k]
                   - f_7 * smf1_249[k]
                   + f_3 * pc_x[k] * smg_369[k];

        t_514[k] = f_11 * slg_370[k]
                   + f_3 * pc_x[k] * smg_370[k];

        t_515[k] = f_11 * slg_371[k]
                   + f_3 * pc_x[k] * smg_371[k];
    }

#pragma omp simd aligned(t_516, t_517, t_518, t_519, pc_x, pc_y, slg_280, slg_372, slg_373, \
                         slg_374, smf0_246, smf1_246, smg_370, smg_372, smg_373, \
                         smg_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_516[k] = f_11 * slg_372[k]
                   + f_3 * pc_x[k] * smg_372[k];

        t_517[k] = f_11 * slg_373[k]
                   + f_3 * pc_x[k] * smg_373[k];

        t_518[k] = f_11 * slg_374[k]
                   + f_3 * pc_x[k] * smg_374[k];

        t_519[k] = f_11 * slg_280[k]
                   + f_1 * smf0_246[k]
                   - f_2 * smf1_246[k]
                   + f_3 * pc_y[k] * smg_370[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, pc_y, pc_z, slg_265, slg_282, slg_283, smf0_248, \
                         smf0_249, smf1_248, smf1_249, smg_370, smg_372, \
                         smg_373 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = f_11 * slg_265[k]
                   + f_3 * pc_z[k] * smg_370[k];

        t_521[k] = f_11 * slg_282[k]
                   + f_4 * smf0_248[k]
                   - f_5 * smf1_248[k]
                   + f_3 * pc_y[k] * smg_372[k];

        t_522[k] = f_11 * slg_283[k]
                   + f_6 * smf0_249[k]
                   - f_7 * smf1_249[k]
                   + f_3 * pc_y[k] * smg_373[k];
    }

#pragma omp simd aligned(t_523, t_524, t_525, pc_x, pc_y, pc_z, slg_269, slg_284, slg_375, \
                         smf0_249, smf0_250, smf1_249, smf1_250, smg_374, \
                         smg_375 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_523[k] = f_11 * slg_284[k]
                   + f_3 * pc_y[k] * smg_374[k];

        t_524[k] = f_11 * slg_269[k]
                   + f_1 * smf0_249[k]
                   - f_2 * smf1_249[k]
                   + f_3 * pc_z[k] * smg_374[k];

        t_525[k] = f_11 * slg_375[k]
                   + f_1 * smf0_250[k]
                   - f_2 * smf1_250[k]
                   + f_3 * pc_x[k] * smg_375[k];
    }

#pragma omp simd aligned(t_526, t_527, t_528, t_529, pc_x, pc_y, pc_z, slg_270, slg_285, \
                         slg_287, slg_378, smf0_253, smf1_253, smg_375, smg_377, \
                         smg_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_526[k] = f_10 * slg_285[k]
                   + f_3 * pc_y[k] * smg_375[k];

        t_527[k] = f_16 * slg_270[k]
                   + f_3 * pc_z[k] * smg_375[k];

        t_528[k] = f_11 * slg_378[k]
                   + f_4 * smf0_253[k]
                   - f_5 * smf1_253[k]
                   + f_3 * pc_x[k] * smg_378[k];

        t_529[k] = f_10 * slg_287[k]
                   + f_3 * pc_y[k] * smg_377[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, pc_x, pc_z, slg_273, slg_380, slg_381, smf0_255, \
                         smf0_256, smf1_255, smf1_256, smg_378, smg_380, \
                         smg_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = f_11 * slg_380[k]
                   + f_4 * smf0_255[k]
                   - f_5 * smf1_255[k]
                   + f_3 * pc_x[k] * smg_380[k];

        t_531[k] = f_11 * slg_381[k]
                   + f_6 * smf0_256[k]
                   - f_7 * smf1_256[k]
                   + f_3 * pc_x[k] * smg_381[k];

        t_532[k] = f_16 * slg_273[k]
                   + f_3 * pc_z[k] * smg_378[k];
    }

#pragma omp simd aligned(t_533, t_534, t_535, t_536, pc_x, pc_y, slg_290, slg_384, slg_385, \
                         slg_386, smf0_259, smf1_259, smg_380, smg_384, smg_385, \
                         smg_386 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_533[k] = f_10 * slg_290[k]
                   + f_3 * pc_y[k] * smg_380[k];

        t_534[k] = f_11 * slg_384[k]
                   + f_6 * smf0_259[k]
                   - f_7 * smf1_259[k]
                   + f_3 * pc_x[k] * smg_384[k];

        t_535[k] = f_11 * slg_385[k]
                   + f_3 * pc_x[k] * smg_385[k];

        t_536[k] = f_11 * slg_386[k]
                   + f_3 * pc_x[k] * smg_386[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, t_540, pc_x, pc_y, slg_295, slg_387, slg_388, \
                         slg_389, smf0_256, smf1_256, smg_385, smg_387, smg_388, \
                         smg_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = f_11 * slg_387[k]
                   + f_3 * pc_x[k] * smg_387[k];

        t_538[k] = f_11 * slg_388[k]
                   + f_3 * pc_x[k] * smg_388[k];

        t_539[k] = f_11 * slg_389[k]
                   + f_3 * pc_x[k] * smg_389[k];

        t_540[k] = f_10 * slg_295[k]
                   + f_1 * smf0_256[k]
                   - f_2 * smf1_256[k]
                   + f_3 * pc_y[k] * smg_385[k];
    }

#pragma omp simd aligned(t_541, t_542, t_543, pc_y, pc_z, slg_280, slg_297, slg_298, smf0_258, \
                         smf0_259, smf1_258, smf1_259, smg_385, smg_387, \
                         smg_388 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_541[k] = f_16 * slg_280[k]
                   + f_3 * pc_z[k] * smg_385[k];

        t_542[k] = f_10 * slg_297[k]
                   + f_4 * smf0_258[k]
                   - f_5 * smf1_258[k]
                   + f_3 * pc_y[k] * smg_387[k];

        t_543[k] = f_10 * slg_298[k]
                   + f_6 * smf0_259[k]
                   - f_7 * smf1_259[k]
                   + f_3 * pc_y[k] * smg_388[k];
    }

#pragma omp simd aligned(t_544, t_545, t_546, t_547, pb_y, pc_y, pc_z, slh0_420, slg_284, \
                         slg_299, slg_300, slh1_420, smf0_259, smf1_259, smg_389, \
                         smg_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_544[k] = f_10 * slg_299[k]
                   + f_3 * pc_y[k] * smg_389[k];

        t_545[k] = f_16 * slg_284[k]
                   + f_1 * smf0_259[k]
                   - f_2 * smf1_259[k]
                   + f_3 * pc_z[k] * smg_389[k];

        t_546[k] = pb_y[k] * slh0_420[k]
                   - f_8 * pc_y[k] * slh1_420[k];

        t_547[k] = f_9 * slg_300[k]
                   + f_3 * pc_y[k] * smg_390[k];
    }

#pragma omp simd aligned(t_548, t_549, t_550, t_551, pb_y, pc_y, pc_z, slh0_423, slh0_425, \
                         slg_285, slg_301, slg_302, slh1_423, slh1_425, smg_390, \
                         smg_392 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_548[k] = f_15 * slg_285[k]
                   + f_3 * pc_z[k] * smg_390[k];

        t_549[k] = pb_y[k] * slh0_423[k]
                   + f_10 * slg_301[k]
                   - f_8 * pc_y[k] * slh1_423[k];

        t_550[k] = f_9 * slg_302[k]
                   + f_3 * pc_y[k] * smg_392[k];

        t_551[k] = pb_y[k] * slh0_425[k]
                   - f_8 * pc_y[k] * slh1_425[k];
    }

#pragma omp simd aligned(t_552, t_553, t_554, t_555, pb_y, pc_y, pc_z, slh0_426, slh0_429, \
                         slg_288, slg_303, slg_305, slh1_426, slh1_429, smg_393, \
                         smg_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_552[k] = pb_y[k] * slh0_426[k]
                   + f_11 * slg_303[k]
                   - f_8 * pc_y[k] * slh1_426[k];

        t_553[k] = f_15 * slg_288[k]
                   + f_3 * pc_z[k] * smg_393[k];

        t_554[k] = f_9 * slg_305[k]
                   + f_3 * pc_y[k] * smg_395[k];

        t_555[k] = pb_y[k] * slh0_429[k]
                   - f_8 * pc_y[k] * slh1_429[k];
    }

#pragma omp simd aligned(t_556, t_557, t_558, t_559, t_560, pc_x, slg_400, slg_401, slg_402, \
                         slg_403, slg_404, smg_400, smg_401, smg_402, smg_403, \
                         smg_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_556[k] = f_11 * slg_400[k]
                   + f_3 * pc_x[k] * smg_400[k];

        t_557[k] = f_11 * slg_401[k]
                   + f_3 * pc_x[k] * smg_401[k];

        t_558[k] = f_11 * slg_402[k]
                   + f_3 * pc_x[k] * smg_402[k];

        t_559[k] = f_11 * slg_403[k]
                   + f_3 * pc_x[k] * smg_403[k];

        t_560[k] = f_11 * slg_404[k]
                   + f_3 * pc_x[k] * smg_404[k];
    }

#pragma omp simd aligned(t_561, t_562, t_563, pc_y, pc_z, slg_295, slg_310, slg_312, smf0_266, \
                         smf0_268, smf1_266, smf1_268, smg_400, \
                         smg_402 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_561[k] = f_9 * slg_310[k]
                   + f_1 * smf0_266[k]
                   - f_2 * smf1_266[k]
                   + f_3 * pc_y[k] * smg_400[k];

        t_562[k] = f_15 * slg_295[k]
                   + f_3 * pc_z[k] * smg_400[k];

        t_563[k] = f_9 * slg_312[k]
                   + f_4 * smf0_268[k]
                   - f_5 * smf1_268[k]
                   + f_3 * pc_y[k] * smg_402[k];
    }

#pragma omp simd aligned(t_564, t_565, t_566, pb_y, pc_y, slh0_440, slg_313, slg_314, \
                         slh1_440, smf0_269, smf1_269, smg_403, \
                         smg_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_564[k] = f_9 * slg_313[k]
                   + f_6 * smf0_269[k]
                   - f_7 * smf1_269[k]
                   + f_3 * pc_y[k] * smg_403[k];

        t_565[k] = f_9 * slg_314[k]
                   + f_3 * pc_y[k] * smg_404[k];

        t_566[k] = pb_y[k] * slh0_440[k]
                   - f_8 * pc_y[k] * slh1_440[k];
    }

#pragma omp simd aligned(t_567, t_568, t_569, t_570, pc_x, pc_y, pc_z, slg_300, slg_405, \
                         slg_408, smf0_270, smf0_273, smf1_270, smf1_273, smg_405, \
                         smg_408 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_567[k] = f_11 * slg_405[k]
                   + f_1 * smf0_270[k]
                   - f_2 * smf1_270[k]
                   + f_3 * pc_x[k] * smg_405[k];

        t_568[k] = f_3 * pc_y[k] * smg_405[k];

        t_569[k] = f_14 * slg_300[k]
                   + f_3 * pc_z[k] * smg_405[k];

        t_570[k] = f_11 * slg_408[k]
                   + f_4 * smf0_273[k]
                   - f_5 * smf1_273[k]
                   + f_3 * pc_x[k] * smg_408[k];
    }

#pragma omp simd aligned(t_571, t_572, t_573, pc_x, pc_y, slg_410, slg_411, smf0_275, \
                         smf0_276, smf1_275, smf1_276, smg_407, smg_410, \
                         smg_411 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_571[k] = f_3 * pc_y[k] * smg_407[k];

        t_572[k] = f_11 * slg_410[k]
                   + f_4 * smf0_275[k]
                   - f_5 * smf1_275[k]
                   + f_3 * pc_x[k] * smg_410[k];

        t_573[k] = f_11 * slg_411[k]
                   + f_6 * smf0_276[k]
                   - f_7 * smf1_276[k]
                   + f_3 * pc_x[k] * smg_411[k];
    }

#pragma omp simd aligned(t_574, t_575, t_576, t_577, pc_x, pc_y, pc_z, slg_303, slg_414, \
                         slg_415, smf0_279, smf1_279, smg_408, smg_410, smg_414, \
                         smg_415 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_574[k] = f_14 * slg_303[k]
                   + f_3 * pc_z[k] * smg_408[k];

        t_575[k] = f_3 * pc_y[k] * smg_410[k];

        t_576[k] = f_11 * slg_414[k]
                   + f_6 * smf0_279[k]
                   - f_7 * smf1_279[k]
                   + f_3 * pc_x[k] * smg_414[k];

        t_577[k] = f_11 * slg_415[k]
                   + f_3 * pc_x[k] * smg_415[k];
    }

#pragma omp simd aligned(t_578, t_579, t_580, t_581, pc_x, slg_416, slg_417, slg_418, slg_419, \
                         smg_416, smg_417, smg_418, smg_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_578[k] = f_11 * slg_416[k]
                   + f_3 * pc_x[k] * smg_416[k];

        t_579[k] = f_11 * slg_417[k]
                   + f_3 * pc_x[k] * smg_417[k];

        t_580[k] = f_11 * slg_418[k]
                   + f_3 * pc_x[k] * smg_418[k];

        t_581[k] = f_11 * slg_419[k]
                   + f_3 * pc_x[k] * smg_419[k];
    }

#pragma omp simd aligned(t_582, t_583, t_584, t_585, pc_y, pc_z, slg_310, smf0_276, smf0_278, \
                         smf0_279, smf1_276, smf1_278, smf1_279, smg_415, smg_417, \
                         smg_418 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_582[k] = f_1 * smf0_276[k]
                   - f_2 * smf1_276[k]
                   + f_3 * pc_y[k] * smg_415[k];

        t_583[k] = f_14 * slg_310[k]
                   + f_3 * pc_z[k] * smg_415[k];

        t_584[k] = f_4 * smf0_278[k]
                   - f_5 * smf1_278[k]
                   + f_3 * pc_y[k] * smg_417[k];

        t_585[k] = f_6 * smf0_279[k]
                   - f_7 * smf1_279[k]
                   + f_3 * pc_y[k] * smg_418[k];
    }
}

static auto
compute_prim_smh_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t slh0,
                                                          const size_t slg, const size_t slh1,
                                                          const size_t smf0, const size_t smf1,
                                                          const size_t smg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.0 / gamma;
    const auto f_5 = p / (gamma * q);
    const auto f_6 = 0.5 / gamma;
    const auto f_7 = 0.5 * p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_13 = 3.5 / q;
    const auto f_14 = 3.0 / q;
    const auto f_15 = 2.5 / q;
    const auto f_16 = 2.0 / q;

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
    auto *t_690 = buffer.data(target + 690);
    auto *t_691 = buffer.data(target + 691);
    auto *t_692 = buffer.data(target + 692);
    auto *t_693 = buffer.data(target + 693);

    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *slh0_441 = buffer.data(slh0 + 441);
    const auto *slh0_444 = buffer.data(slh0 + 444);
    const auto *slh0_447 = buffer.data(slh0 + 447);
    const auto *slh0_456 = buffer.data(slh0 + 456);

    const auto *slg_314 = buffer.data(slg + 314);
    const auto *slg_315 = buffer.data(slg + 315);
    const auto *slg_317 = buffer.data(slg + 317);
    const auto *slg_318 = buffer.data(slg + 318);
    const auto *slg_320 = buffer.data(slg + 320);
    const auto *slg_325 = buffer.data(slg + 325);
    const auto *slg_327 = buffer.data(slg + 327);
    const auto *slg_328 = buffer.data(slg + 328);
    const auto *slg_329 = buffer.data(slg + 329);
    const auto *slg_330 = buffer.data(slg + 330);
    const auto *slg_332 = buffer.data(slg + 332);
    const auto *slg_333 = buffer.data(slg + 333);
    const auto *slg_335 = buffer.data(slg + 335);
    const auto *slg_340 = buffer.data(slg + 340);
    const auto *slg_342 = buffer.data(slg + 342);
    const auto *slg_343 = buffer.data(slg + 343);
    const auto *slg_344 = buffer.data(slg + 344);
    const auto *slg_345 = buffer.data(slg + 345);
    const auto *slg_347 = buffer.data(slg + 347);
    const auto *slg_348 = buffer.data(slg + 348);
    const auto *slg_350 = buffer.data(slg + 350);
    const auto *slg_355 = buffer.data(slg + 355);
    const auto *slg_357 = buffer.data(slg + 357);
    const auto *slg_358 = buffer.data(slg + 358);
    const auto *slg_359 = buffer.data(slg + 359);
    const auto *slg_360 = buffer.data(slg + 360);
    const auto *slg_362 = buffer.data(slg + 362);
    const auto *slg_363 = buffer.data(slg + 363);
    const auto *slg_365 = buffer.data(slg + 365);
    const auto *slg_370 = buffer.data(slg + 370);
    const auto *slg_372 = buffer.data(slg + 372);
    const auto *slg_373 = buffer.data(slg + 373);
    const auto *slg_374 = buffer.data(slg + 374);
    const auto *slg_375 = buffer.data(slg + 375);
    const auto *slg_377 = buffer.data(slg + 377);
    const auto *slg_380 = buffer.data(slg + 380);
    const auto *slg_385 = buffer.data(slg + 385);
    const auto *slg_387 = buffer.data(slg + 387);
    const auto *slg_388 = buffer.data(slg + 388);
    const auto *slg_389 = buffer.data(slg + 389);
    const auto *slg_420 = buffer.data(slg + 420);
    const auto *slg_423 = buffer.data(slg + 423);
    const auto *slg_425 = buffer.data(slg + 425);
    const auto *slg_426 = buffer.data(slg + 426);
    const auto *slg_429 = buffer.data(slg + 429);
    const auto *slg_430 = buffer.data(slg + 430);
    const auto *slg_431 = buffer.data(slg + 431);
    const auto *slg_432 = buffer.data(slg + 432);
    const auto *slg_433 = buffer.data(slg + 433);
    const auto *slg_434 = buffer.data(slg + 434);
    const auto *slg_440 = buffer.data(slg + 440);
    const auto *slg_444 = buffer.data(slg + 444);
    const auto *slg_445 = buffer.data(slg + 445);
    const auto *slg_446 = buffer.data(slg + 446);
    const auto *slg_447 = buffer.data(slg + 447);
    const auto *slg_448 = buffer.data(slg + 448);
    const auto *slg_449 = buffer.data(slg + 449);
    const auto *slg_450 = buffer.data(slg + 450);
    const auto *slg_453 = buffer.data(slg + 453);
    const auto *slg_455 = buffer.data(slg + 455);
    const auto *slg_456 = buffer.data(slg + 456);
    const auto *slg_459 = buffer.data(slg + 459);
    const auto *slg_460 = buffer.data(slg + 460);
    const auto *slg_461 = buffer.data(slg + 461);
    const auto *slg_462 = buffer.data(slg + 462);
    const auto *slg_463 = buffer.data(slg + 463);
    const auto *slg_464 = buffer.data(slg + 464);
    const auto *slg_465 = buffer.data(slg + 465);
    const auto *slg_468 = buffer.data(slg + 468);
    const auto *slg_470 = buffer.data(slg + 470);
    const auto *slg_471 = buffer.data(slg + 471);
    const auto *slg_474 = buffer.data(slg + 474);
    const auto *slg_475 = buffer.data(slg + 475);
    const auto *slg_476 = buffer.data(slg + 476);
    const auto *slg_477 = buffer.data(slg + 477);
    const auto *slg_478 = buffer.data(slg + 478);
    const auto *slg_479 = buffer.data(slg + 479);
    const auto *slg_480 = buffer.data(slg + 480);
    const auto *slg_483 = buffer.data(slg + 483);
    const auto *slg_485 = buffer.data(slg + 485);
    const auto *slg_486 = buffer.data(slg + 486);
    const auto *slg_489 = buffer.data(slg + 489);
    const auto *slg_490 = buffer.data(slg + 490);
    const auto *slg_491 = buffer.data(slg + 491);
    const auto *slg_492 = buffer.data(slg + 492);
    const auto *slg_493 = buffer.data(slg + 493);
    const auto *slg_494 = buffer.data(slg + 494);
    const auto *slg_495 = buffer.data(slg + 495);

    const auto *slh1_441 = buffer.data(slh1 + 441);
    const auto *slh1_444 = buffer.data(slh1 + 444);
    const auto *slh1_447 = buffer.data(slh1 + 447);
    const auto *slh1_456 = buffer.data(slh1 + 456);

    const auto *smf0_279 = buffer.data(smf0 + 279);
    const auto *smf0_280 = buffer.data(smf0 + 280);
    const auto *smf0_283 = buffer.data(smf0 + 283);
    const auto *smf0_285 = buffer.data(smf0 + 285);
    const auto *smf0_286 = buffer.data(smf0 + 286);
    const auto *smf0_288 = buffer.data(smf0 + 288);
    const auto *smf0_289 = buffer.data(smf0 + 289);
    const auto *smf0_295 = buffer.data(smf0 + 295);
    const auto *smf0_298 = buffer.data(smf0 + 298);
    const auto *smf0_299 = buffer.data(smf0 + 299);
    const auto *smf0_300 = buffer.data(smf0 + 300);
    const auto *smf0_303 = buffer.data(smf0 + 303);
    const auto *smf0_305 = buffer.data(smf0 + 305);
    const auto *smf0_306 = buffer.data(smf0 + 306);
    const auto *smf0_308 = buffer.data(smf0 + 308);
    const auto *smf0_309 = buffer.data(smf0 + 309);
    const auto *smf0_310 = buffer.data(smf0 + 310);
    const auto *smf0_313 = buffer.data(smf0 + 313);
    const auto *smf0_315 = buffer.data(smf0 + 315);
    const auto *smf0_316 = buffer.data(smf0 + 316);
    const auto *smf0_318 = buffer.data(smf0 + 318);
    const auto *smf0_319 = buffer.data(smf0 + 319);
    const auto *smf0_320 = buffer.data(smf0 + 320);
    const auto *smf0_323 = buffer.data(smf0 + 323);
    const auto *smf0_325 = buffer.data(smf0 + 325);
    const auto *smf0_326 = buffer.data(smf0 + 326);
    const auto *smf0_328 = buffer.data(smf0 + 328);
    const auto *smf0_329 = buffer.data(smf0 + 329);
    const auto *smf0_330 = buffer.data(smf0 + 330);

    const auto *smf1_279 = buffer.data(smf1 + 279);
    const auto *smf1_280 = buffer.data(smf1 + 280);
    const auto *smf1_283 = buffer.data(smf1 + 283);
    const auto *smf1_285 = buffer.data(smf1 + 285);
    const auto *smf1_286 = buffer.data(smf1 + 286);
    const auto *smf1_288 = buffer.data(smf1 + 288);
    const auto *smf1_289 = buffer.data(smf1 + 289);
    const auto *smf1_295 = buffer.data(smf1 + 295);
    const auto *smf1_298 = buffer.data(smf1 + 298);
    const auto *smf1_299 = buffer.data(smf1 + 299);
    const auto *smf1_300 = buffer.data(smf1 + 300);
    const auto *smf1_303 = buffer.data(smf1 + 303);
    const auto *smf1_305 = buffer.data(smf1 + 305);
    const auto *smf1_306 = buffer.data(smf1 + 306);
    const auto *smf1_308 = buffer.data(smf1 + 308);
    const auto *smf1_309 = buffer.data(smf1 + 309);
    const auto *smf1_310 = buffer.data(smf1 + 310);
    const auto *smf1_313 = buffer.data(smf1 + 313);
    const auto *smf1_315 = buffer.data(smf1 + 315);
    const auto *smf1_316 = buffer.data(smf1 + 316);
    const auto *smf1_318 = buffer.data(smf1 + 318);
    const auto *smf1_319 = buffer.data(smf1 + 319);
    const auto *smf1_320 = buffer.data(smf1 + 320);
    const auto *smf1_323 = buffer.data(smf1 + 323);
    const auto *smf1_325 = buffer.data(smf1 + 325);
    const auto *smf1_326 = buffer.data(smf1 + 326);
    const auto *smf1_328 = buffer.data(smf1 + 328);
    const auto *smf1_329 = buffer.data(smf1 + 329);
    const auto *smf1_330 = buffer.data(smf1 + 330);

    const auto *smg_419 = buffer.data(smg + 419);
    const auto *smg_420 = buffer.data(smg + 420);
    const auto *smg_422 = buffer.data(smg + 422);
    const auto *smg_423 = buffer.data(smg + 423);
    const auto *smg_425 = buffer.data(smg + 425);
    const auto *smg_426 = buffer.data(smg + 426);
    const auto *smg_429 = buffer.data(smg + 429);
    const auto *smg_430 = buffer.data(smg + 430);
    const auto *smg_431 = buffer.data(smg + 431);
    const auto *smg_432 = buffer.data(smg + 432);
    const auto *smg_433 = buffer.data(smg + 433);
    const auto *smg_434 = buffer.data(smg + 434);
    const auto *smg_435 = buffer.data(smg + 435);
    const auto *smg_437 = buffer.data(smg + 437);
    const auto *smg_438 = buffer.data(smg + 438);
    const auto *smg_440 = buffer.data(smg + 440);
    const auto *smg_444 = buffer.data(smg + 444);
    const auto *smg_445 = buffer.data(smg + 445);
    const auto *smg_446 = buffer.data(smg + 446);
    const auto *smg_447 = buffer.data(smg + 447);
    const auto *smg_448 = buffer.data(smg + 448);
    const auto *smg_449 = buffer.data(smg + 449);
    const auto *smg_450 = buffer.data(smg + 450);
    const auto *smg_452 = buffer.data(smg + 452);
    const auto *smg_453 = buffer.data(smg + 453);
    const auto *smg_455 = buffer.data(smg + 455);
    const auto *smg_456 = buffer.data(smg + 456);
    const auto *smg_459 = buffer.data(smg + 459);
    const auto *smg_460 = buffer.data(smg + 460);
    const auto *smg_461 = buffer.data(smg + 461);
    const auto *smg_462 = buffer.data(smg + 462);
    const auto *smg_463 = buffer.data(smg + 463);
    const auto *smg_464 = buffer.data(smg + 464);
    const auto *smg_465 = buffer.data(smg + 465);
    const auto *smg_467 = buffer.data(smg + 467);
    const auto *smg_468 = buffer.data(smg + 468);
    const auto *smg_470 = buffer.data(smg + 470);
    const auto *smg_471 = buffer.data(smg + 471);
    const auto *smg_474 = buffer.data(smg + 474);
    const auto *smg_475 = buffer.data(smg + 475);
    const auto *smg_476 = buffer.data(smg + 476);
    const auto *smg_477 = buffer.data(smg + 477);
    const auto *smg_478 = buffer.data(smg + 478);
    const auto *smg_479 = buffer.data(smg + 479);
    const auto *smg_480 = buffer.data(smg + 480);
    const auto *smg_482 = buffer.data(smg + 482);
    const auto *smg_483 = buffer.data(smg + 483);
    const auto *smg_485 = buffer.data(smg + 485);
    const auto *smg_486 = buffer.data(smg + 486);
    const auto *smg_489 = buffer.data(smg + 489);
    const auto *smg_490 = buffer.data(smg + 490);
    const auto *smg_491 = buffer.data(smg + 491);
    const auto *smg_492 = buffer.data(smg + 492);
    const auto *smg_493 = buffer.data(smg + 493);
    const auto *smg_494 = buffer.data(smg + 494);
    const auto *smg_495 = buffer.data(smg + 495);

#pragma omp simd aligned(t_586, t_587, t_588, t_589, pc_x, pc_y, pc_z, slg_314, slg_315, \
                         slg_420, smf0_279, smf0_280, smf1_279, smf1_280, smg_419, \
                         smg_420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_586[k] = f_3 * pc_y[k] * smg_419[k];

        t_587[k] = f_14 * slg_314[k]
                   + f_1 * smf0_279[k]
                   - f_2 * smf1_279[k]
                   + f_3 * pc_z[k] * smg_419[k];

        t_588[k] = f_10 * slg_420[k]
                   + f_1 * smf0_280[k]
                   - f_2 * smf1_280[k]
                   + f_3 * pc_x[k] * smg_420[k];

        t_589[k] = f_13 * slg_315[k]
                   + f_3 * pc_y[k] * smg_420[k];
    }

#pragma omp simd aligned(t_590, t_591, t_592, pc_x, pc_y, pc_z, slg_317, slg_423, smf0_283, \
                         smf1_283, smg_420, smg_422, smg_423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_590[k] = f_3 * pc_z[k] * smg_420[k];

        t_591[k] = f_10 * slg_423[k]
                   + f_4 * smf0_283[k]
                   - f_5 * smf1_283[k]
                   + f_3 * pc_x[k] * smg_423[k];

        t_592[k] = f_13 * slg_317[k]
                   + f_3 * pc_y[k] * smg_422[k];
    }

#pragma omp simd aligned(t_593, t_594, t_595, pc_x, pc_z, slg_425, slg_426, smf0_285, \
                         smf0_286, smf1_285, smf1_286, smg_423, smg_425, \
                         smg_426 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_593[k] = f_10 * slg_425[k]
                   + f_4 * smf0_285[k]
                   - f_5 * smf1_285[k]
                   + f_3 * pc_x[k] * smg_425[k];

        t_594[k] = f_10 * slg_426[k]
                   + f_6 * smf0_286[k]
                   - f_7 * smf1_286[k]
                   + f_3 * pc_x[k] * smg_426[k];

        t_595[k] = f_3 * pc_z[k] * smg_423[k];
    }

#pragma omp simd aligned(t_596, t_597, t_598, t_599, pc_x, pc_y, slg_320, slg_429, slg_430, \
                         slg_431, smf0_289, smf1_289, smg_425, smg_429, smg_430, \
                         smg_431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_596[k] = f_13 * slg_320[k]
                   + f_3 * pc_y[k] * smg_425[k];

        t_597[k] = f_10 * slg_429[k]
                   + f_6 * smf0_289[k]
                   - f_7 * smf1_289[k]
                   + f_3 * pc_x[k] * smg_429[k];

        t_598[k] = f_10 * slg_430[k]
                   + f_3 * pc_x[k] * smg_430[k];

        t_599[k] = f_10 * slg_431[k]
                   + f_3 * pc_x[k] * smg_431[k];
    }

#pragma omp simd aligned(t_600, t_601, t_602, t_603, pc_x, pc_y, slg_325, slg_432, slg_433, \
                         slg_434, smf0_286, smf1_286, smg_430, smg_432, smg_433, \
                         smg_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = f_10 * slg_432[k]
                   + f_3 * pc_x[k] * smg_432[k];

        t_601[k] = f_10 * slg_433[k]
                   + f_3 * pc_x[k] * smg_433[k];

        t_602[k] = f_10 * slg_434[k]
                   + f_3 * pc_x[k] * smg_434[k];

        t_603[k] = f_13 * slg_325[k]
                   + f_1 * smf0_286[k]
                   - f_2 * smf1_286[k]
                   + f_3 * pc_y[k] * smg_430[k];
    }

#pragma omp simd aligned(t_604, t_605, t_606, pc_y, pc_z, slg_327, slg_328, smf0_288, \
                         smf0_289, smf1_288, smf1_289, smg_430, smg_432, \
                         smg_433 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_604[k] = f_3 * pc_z[k] * smg_430[k];

        t_605[k] = f_13 * slg_327[k]
                   + f_4 * smf0_288[k]
                   - f_5 * smf1_288[k]
                   + f_3 * pc_y[k] * smg_432[k];

        t_606[k] = f_13 * slg_328[k]
                   + f_6 * smf0_289[k]
                   - f_7 * smf1_289[k]
                   + f_3 * pc_y[k] * smg_433[k];
    }

#pragma omp simd aligned(t_607, t_608, t_609, t_610, pb_z, pc_y, pc_z, slh0_441, slg_329, \
                         slg_330, slh1_441, smf0_289, smf1_289, smg_434, \
                         smg_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_607[k] = f_13 * slg_329[k]
                   + f_3 * pc_y[k] * smg_434[k];

        t_608[k] = f_1 * smf0_289[k]
                   - f_2 * smf1_289[k]
                   + f_3 * pc_z[k] * smg_434[k];

        t_609[k] = pb_z[k] * slh0_441[k]
                   - f_8 * pc_z[k] * slh1_441[k];

        t_610[k] = f_14 * slg_330[k]
                   + f_3 * pc_y[k] * smg_435[k];
    }

#pragma omp simd aligned(t_611, t_612, t_613, pb_z, pc_y, pc_z, slh0_444, slg_315, slg_332, \
                         slh1_444, smg_435, smg_437 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_611[k] = f_9 * slg_315[k]
                   + f_3 * pc_z[k] * smg_435[k];

        t_612[k] = pb_z[k] * slh0_444[k]
                   - f_8 * pc_z[k] * slh1_444[k];

        t_613[k] = f_14 * slg_332[k]
                   + f_3 * pc_y[k] * smg_437[k];
    }

#pragma omp simd aligned(t_614, t_615, t_616, pb_z, pc_x, pc_z, slh0_447, slg_318, slg_440, \
                         slh1_447, smf0_295, smf1_295, smg_438, \
                         smg_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_614[k] = f_10 * slg_440[k]
                   + f_4 * smf0_295[k]
                   - f_5 * smf1_295[k]
                   + f_3 * pc_x[k] * smg_440[k];

        t_615[k] = pb_z[k] * slh0_447[k]
                   - f_8 * pc_z[k] * slh1_447[k];

        t_616[k] = f_9 * slg_318[k]
                   + f_3 * pc_z[k] * smg_438[k];
    }

#pragma omp simd aligned(t_617, t_618, t_619, t_620, pc_x, pc_y, slg_335, slg_444, slg_445, \
                         slg_446, smf0_299, smf1_299, smg_440, smg_444, smg_445, \
                         smg_446 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_617[k] = f_14 * slg_335[k]
                   + f_3 * pc_y[k] * smg_440[k];

        t_618[k] = f_10 * slg_444[k]
                   + f_6 * smf0_299[k]
                   - f_7 * smf1_299[k]
                   + f_3 * pc_x[k] * smg_444[k];

        t_619[k] = f_10 * slg_445[k]
                   + f_3 * pc_x[k] * smg_445[k];

        t_620[k] = f_10 * slg_446[k]
                   + f_3 * pc_x[k] * smg_446[k];
    }

#pragma omp simd aligned(t_621, t_622, t_623, t_624, pb_z, pc_x, pc_z, slh0_456, slg_447, \
                         slg_448, slg_449, slh1_456, smg_447, smg_448, \
                         smg_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_621[k] = f_10 * slg_447[k]
                   + f_3 * pc_x[k] * smg_447[k];

        t_622[k] = f_10 * slg_448[k]
                   + f_3 * pc_x[k] * smg_448[k];

        t_623[k] = f_10 * slg_449[k]
                   + f_3 * pc_x[k] * smg_449[k];

        t_624[k] = pb_z[k] * slh0_456[k]
                   - f_8 * pc_z[k] * slh1_456[k];
    }

#pragma omp simd aligned(t_625, t_626, t_627, pc_y, pc_z, slg_325, slg_342, slg_343, smf0_298, \
                         smf0_299, smf1_298, smf1_299, smg_445, smg_447, \
                         smg_448 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_625[k] = f_9 * slg_325[k]
                   + f_3 * pc_z[k] * smg_445[k];

        t_626[k] = f_14 * slg_342[k]
                   + f_4 * smf0_298[k]
                   - f_5 * smf1_298[k]
                   + f_3 * pc_y[k] * smg_447[k];

        t_627[k] = f_14 * slg_343[k]
                   + f_6 * smf0_299[k]
                   - f_7 * smf1_299[k]
                   + f_3 * pc_y[k] * smg_448[k];
    }

#pragma omp simd aligned(t_628, t_629, t_630, pc_x, pc_y, pc_z, slg_329, slg_344, slg_450, \
                         smf0_299, smf0_300, smf1_299, smf1_300, smg_449, \
                         smg_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_628[k] = f_14 * slg_344[k]
                   + f_3 * pc_y[k] * smg_449[k];

        t_629[k] = f_9 * slg_329[k]
                   + f_1 * smf0_299[k]
                   - f_2 * smf1_299[k]
                   + f_3 * pc_z[k] * smg_449[k];

        t_630[k] = f_10 * slg_450[k]
                   + f_1 * smf0_300[k]
                   - f_2 * smf1_300[k]
                   + f_3 * pc_x[k] * smg_450[k];
    }

#pragma omp simd aligned(t_631, t_632, t_633, t_634, pc_x, pc_y, pc_z, slg_330, slg_345, \
                         slg_347, slg_453, smf0_303, smf1_303, smg_450, smg_452, \
                         smg_453 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_631[k] = f_15 * slg_345[k]
                   + f_3 * pc_y[k] * smg_450[k];

        t_632[k] = f_10 * slg_330[k]
                   + f_3 * pc_z[k] * smg_450[k];

        t_633[k] = f_10 * slg_453[k]
                   + f_4 * smf0_303[k]
                   - f_5 * smf1_303[k]
                   + f_3 * pc_x[k] * smg_453[k];

        t_634[k] = f_15 * slg_347[k]
                   + f_3 * pc_y[k] * smg_452[k];
    }

#pragma omp simd aligned(t_635, t_636, t_637, pc_x, pc_z, slg_333, slg_455, slg_456, smf0_305, \
                         smf0_306, smf1_305, smf1_306, smg_453, smg_455, \
                         smg_456 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_635[k] = f_10 * slg_455[k]
                   + f_4 * smf0_305[k]
                   - f_5 * smf1_305[k]
                   + f_3 * pc_x[k] * smg_455[k];

        t_636[k] = f_10 * slg_456[k]
                   + f_6 * smf0_306[k]
                   - f_7 * smf1_306[k]
                   + f_3 * pc_x[k] * smg_456[k];

        t_637[k] = f_10 * slg_333[k]
                   + f_3 * pc_z[k] * smg_453[k];
    }

#pragma omp simd aligned(t_638, t_639, t_640, t_641, pc_x, pc_y, slg_350, slg_459, slg_460, \
                         slg_461, smf0_309, smf1_309, smg_455, smg_459, smg_460, \
                         smg_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_638[k] = f_15 * slg_350[k]
                   + f_3 * pc_y[k] * smg_455[k];

        t_639[k] = f_10 * slg_459[k]
                   + f_6 * smf0_309[k]
                   - f_7 * smf1_309[k]
                   + f_3 * pc_x[k] * smg_459[k];

        t_640[k] = f_10 * slg_460[k]
                   + f_3 * pc_x[k] * smg_460[k];

        t_641[k] = f_10 * slg_461[k]
                   + f_3 * pc_x[k] * smg_461[k];
    }

#pragma omp simd aligned(t_642, t_643, t_644, t_645, pc_x, pc_y, slg_355, slg_462, slg_463, \
                         slg_464, smf0_306, smf1_306, smg_460, smg_462, smg_463, \
                         smg_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_642[k] = f_10 * slg_462[k]
                   + f_3 * pc_x[k] * smg_462[k];

        t_643[k] = f_10 * slg_463[k]
                   + f_3 * pc_x[k] * smg_463[k];

        t_644[k] = f_10 * slg_464[k]
                   + f_3 * pc_x[k] * smg_464[k];

        t_645[k] = f_15 * slg_355[k]
                   + f_1 * smf0_306[k]
                   - f_2 * smf1_306[k]
                   + f_3 * pc_y[k] * smg_460[k];
    }

#pragma omp simd aligned(t_646, t_647, t_648, pc_y, pc_z, slg_340, slg_357, slg_358, smf0_308, \
                         smf0_309, smf1_308, smf1_309, smg_460, smg_462, \
                         smg_463 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_646[k] = f_10 * slg_340[k]
                   + f_3 * pc_z[k] * smg_460[k];

        t_647[k] = f_15 * slg_357[k]
                   + f_4 * smf0_308[k]
                   - f_5 * smf1_308[k]
                   + f_3 * pc_y[k] * smg_462[k];

        t_648[k] = f_15 * slg_358[k]
                   + f_6 * smf0_309[k]
                   - f_7 * smf1_309[k]
                   + f_3 * pc_y[k] * smg_463[k];
    }

#pragma omp simd aligned(t_649, t_650, t_651, pc_x, pc_y, pc_z, slg_344, slg_359, slg_465, \
                         smf0_309, smf0_310, smf1_309, smf1_310, smg_464, \
                         smg_465 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_649[k] = f_15 * slg_359[k]
                   + f_3 * pc_y[k] * smg_464[k];

        t_650[k] = f_10 * slg_344[k]
                   + f_1 * smf0_309[k]
                   - f_2 * smf1_309[k]
                   + f_3 * pc_z[k] * smg_464[k];

        t_651[k] = f_10 * slg_465[k]
                   + f_1 * smf0_310[k]
                   - f_2 * smf1_310[k]
                   + f_3 * pc_x[k] * smg_465[k];
    }

#pragma omp simd aligned(t_652, t_653, t_654, t_655, pc_x, pc_y, pc_z, slg_345, slg_360, \
                         slg_362, slg_468, smf0_313, smf1_313, smg_465, smg_467, \
                         smg_468 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_652[k] = f_16 * slg_360[k]
                   + f_3 * pc_y[k] * smg_465[k];

        t_653[k] = f_11 * slg_345[k]
                   + f_3 * pc_z[k] * smg_465[k];

        t_654[k] = f_10 * slg_468[k]
                   + f_4 * smf0_313[k]
                   - f_5 * smf1_313[k]
                   + f_3 * pc_x[k] * smg_468[k];

        t_655[k] = f_16 * slg_362[k]
                   + f_3 * pc_y[k] * smg_467[k];
    }

#pragma omp simd aligned(t_656, t_657, t_658, pc_x, pc_z, slg_348, slg_470, slg_471, smf0_315, \
                         smf0_316, smf1_315, smf1_316, smg_468, smg_470, \
                         smg_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_656[k] = f_10 * slg_470[k]
                   + f_4 * smf0_315[k]
                   - f_5 * smf1_315[k]
                   + f_3 * pc_x[k] * smg_470[k];

        t_657[k] = f_10 * slg_471[k]
                   + f_6 * smf0_316[k]
                   - f_7 * smf1_316[k]
                   + f_3 * pc_x[k] * smg_471[k];

        t_658[k] = f_11 * slg_348[k]
                   + f_3 * pc_z[k] * smg_468[k];
    }

#pragma omp simd aligned(t_659, t_660, t_661, t_662, pc_x, pc_y, slg_365, slg_474, slg_475, \
                         slg_476, smf0_319, smf1_319, smg_470, smg_474, smg_475, \
                         smg_476 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_659[k] = f_16 * slg_365[k]
                   + f_3 * pc_y[k] * smg_470[k];

        t_660[k] = f_10 * slg_474[k]
                   + f_6 * smf0_319[k]
                   - f_7 * smf1_319[k]
                   + f_3 * pc_x[k] * smg_474[k];

        t_661[k] = f_10 * slg_475[k]
                   + f_3 * pc_x[k] * smg_475[k];

        t_662[k] = f_10 * slg_476[k]
                   + f_3 * pc_x[k] * smg_476[k];
    }

#pragma omp simd aligned(t_663, t_664, t_665, t_666, pc_x, pc_y, slg_370, slg_477, slg_478, \
                         slg_479, smf0_316, smf1_316, smg_475, smg_477, smg_478, \
                         smg_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_663[k] = f_10 * slg_477[k]
                   + f_3 * pc_x[k] * smg_477[k];

        t_664[k] = f_10 * slg_478[k]
                   + f_3 * pc_x[k] * smg_478[k];

        t_665[k] = f_10 * slg_479[k]
                   + f_3 * pc_x[k] * smg_479[k];

        t_666[k] = f_16 * slg_370[k]
                   + f_1 * smf0_316[k]
                   - f_2 * smf1_316[k]
                   + f_3 * pc_y[k] * smg_475[k];
    }

#pragma omp simd aligned(t_667, t_668, t_669, pc_y, pc_z, slg_355, slg_372, slg_373, smf0_318, \
                         smf0_319, smf1_318, smf1_319, smg_475, smg_477, \
                         smg_478 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_667[k] = f_11 * slg_355[k]
                   + f_3 * pc_z[k] * smg_475[k];

        t_668[k] = f_16 * slg_372[k]
                   + f_4 * smf0_318[k]
                   - f_5 * smf1_318[k]
                   + f_3 * pc_y[k] * smg_477[k];

        t_669[k] = f_16 * slg_373[k]
                   + f_6 * smf0_319[k]
                   - f_7 * smf1_319[k]
                   + f_3 * pc_y[k] * smg_478[k];
    }

#pragma omp simd aligned(t_670, t_671, t_672, pc_x, pc_y, pc_z, slg_359, slg_374, slg_480, \
                         smf0_319, smf0_320, smf1_319, smf1_320, smg_479, \
                         smg_480 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = f_16 * slg_374[k]
                   + f_3 * pc_y[k] * smg_479[k];

        t_671[k] = f_11 * slg_359[k]
                   + f_1 * smf0_319[k]
                   - f_2 * smf1_319[k]
                   + f_3 * pc_z[k] * smg_479[k];

        t_672[k] = f_10 * slg_480[k]
                   + f_1 * smf0_320[k]
                   - f_2 * smf1_320[k]
                   + f_3 * pc_x[k] * smg_480[k];
    }

#pragma omp simd aligned(t_673, t_674, t_675, t_676, pc_x, pc_y, pc_z, slg_360, slg_375, \
                         slg_377, slg_483, smf0_323, smf1_323, smg_480, smg_482, \
                         smg_483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_673[k] = f_11 * slg_375[k]
                   + f_3 * pc_y[k] * smg_480[k];

        t_674[k] = f_16 * slg_360[k]
                   + f_3 * pc_z[k] * smg_480[k];

        t_675[k] = f_10 * slg_483[k]
                   + f_4 * smf0_323[k]
                   - f_5 * smf1_323[k]
                   + f_3 * pc_x[k] * smg_483[k];

        t_676[k] = f_11 * slg_377[k]
                   + f_3 * pc_y[k] * smg_482[k];
    }

#pragma omp simd aligned(t_677, t_678, t_679, pc_x, pc_z, slg_363, slg_485, slg_486, smf0_325, \
                         smf0_326, smf1_325, smf1_326, smg_483, smg_485, \
                         smg_486 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_677[k] = f_10 * slg_485[k]
                   + f_4 * smf0_325[k]
                   - f_5 * smf1_325[k]
                   + f_3 * pc_x[k] * smg_485[k];

        t_678[k] = f_10 * slg_486[k]
                   + f_6 * smf0_326[k]
                   - f_7 * smf1_326[k]
                   + f_3 * pc_x[k] * smg_486[k];

        t_679[k] = f_16 * slg_363[k]
                   + f_3 * pc_z[k] * smg_483[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, t_683, pc_x, pc_y, slg_380, slg_489, slg_490, \
                         slg_491, smf0_329, smf1_329, smg_485, smg_489, smg_490, \
                         smg_491 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = f_11 * slg_380[k]
                   + f_3 * pc_y[k] * smg_485[k];

        t_681[k] = f_10 * slg_489[k]
                   + f_6 * smf0_329[k]
                   - f_7 * smf1_329[k]
                   + f_3 * pc_x[k] * smg_489[k];

        t_682[k] = f_10 * slg_490[k]
                   + f_3 * pc_x[k] * smg_490[k];

        t_683[k] = f_10 * slg_491[k]
                   + f_3 * pc_x[k] * smg_491[k];
    }

#pragma omp simd aligned(t_684, t_685, t_686, t_687, pc_x, pc_y, slg_385, slg_492, slg_493, \
                         slg_494, smf0_326, smf1_326, smg_490, smg_492, smg_493, \
                         smg_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_684[k] = f_10 * slg_492[k]
                   + f_3 * pc_x[k] * smg_492[k];

        t_685[k] = f_10 * slg_493[k]
                   + f_3 * pc_x[k] * smg_493[k];

        t_686[k] = f_10 * slg_494[k]
                   + f_3 * pc_x[k] * smg_494[k];

        t_687[k] = f_11 * slg_385[k]
                   + f_1 * smf0_326[k]
                   - f_2 * smf1_326[k]
                   + f_3 * pc_y[k] * smg_490[k];
    }

#pragma omp simd aligned(t_688, t_689, t_690, pc_y, pc_z, slg_370, slg_387, slg_388, smf0_328, \
                         smf0_329, smf1_328, smf1_329, smg_490, smg_492, \
                         smg_493 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_688[k] = f_16 * slg_370[k]
                   + f_3 * pc_z[k] * smg_490[k];

        t_689[k] = f_11 * slg_387[k]
                   + f_4 * smf0_328[k]
                   - f_5 * smf1_328[k]
                   + f_3 * pc_y[k] * smg_492[k];

        t_690[k] = f_11 * slg_388[k]
                   + f_6 * smf0_329[k]
                   - f_7 * smf1_329[k]
                   + f_3 * pc_y[k] * smg_493[k];
    }

#pragma omp simd aligned(t_691, t_692, t_693, pc_x, pc_y, pc_z, slg_374, slg_389, slg_495, \
                         smf0_329, smf0_330, smf1_329, smf1_330, smg_494, \
                         smg_495 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_691[k] = f_11 * slg_389[k]
                   + f_3 * pc_y[k] * smg_494[k];

        t_692[k] = f_16 * slg_374[k]
                   + f_1 * smf0_329[k]
                   - f_2 * smf1_329[k]
                   + f_3 * pc_z[k] * smg_494[k];

        t_693[k] = f_10 * slg_495[k]
                   + f_1 * smf0_330[k]
                   - f_2 * smf1_330[k]
                   + f_3 * pc_x[k] * smg_495[k];
    }
}

static auto
compute_prim_smh_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t slh0,
                                                          const size_t slg, const size_t slh1,
                                                          const size_t smf0, const size_t smf1,
                                                          const size_t smg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.0 / gamma;
    const auto f_5 = p / (gamma * q);
    const auto f_6 = 0.5 / gamma;
    const auto f_7 = 0.5 * p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_12 = 4.0 / q;
    const auto f_13 = 3.5 / q;
    const auto f_14 = 3.0 / q;
    const auto f_15 = 2.5 / q;

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
    auto *t_811 = buffer.data(target + 811);
    auto *t_812 = buffer.data(target + 812);
    auto *t_813 = buffer.data(target + 813);
    auto *t_814 = buffer.data(target + 814);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *slh0_567 = buffer.data(slh0 + 567);
    const auto *slh0_570 = buffer.data(slh0 + 570);
    const auto *slh0_572 = buffer.data(slh0 + 572);
    const auto *slh0_573 = buffer.data(slh0 + 573);
    const auto *slh0_576 = buffer.data(slh0 + 576);
    const auto *slh0_587 = buffer.data(slh0 + 587);
    const auto *slh0_588 = buffer.data(slh0 + 588);
    const auto *slh0_591 = buffer.data(slh0 + 591);
    const auto *slh0_594 = buffer.data(slh0 + 594);
    const auto *slh0_756 = buffer.data(slh0 + 756);
    const auto *slh0_759 = buffer.data(slh0 + 759);
    const auto *slh0_761 = buffer.data(slh0 + 761);
    const auto *slh0_762 = buffer.data(slh0 + 762);
    const auto *slh0_765 = buffer.data(slh0 + 765);
    const auto *slh0_771 = buffer.data(slh0 + 771);
    const auto *slh0_773 = buffer.data(slh0 + 773);
    const auto *slh0_774 = buffer.data(slh0 + 774);
    const auto *slh0_776 = buffer.data(slh0 + 776);
    const auto *slh0_782 = buffer.data(slh0 + 782);
    const auto *slh0_786 = buffer.data(slh0 + 786);
    const auto *slh0_792 = buffer.data(slh0 + 792);
    const auto *slh0_794 = buffer.data(slh0 + 794);
    const auto *slh0_795 = buffer.data(slh0 + 795);
    const auto *slh0_797 = buffer.data(slh0 + 797);
    const auto *slh0_798 = buffer.data(slh0 + 798);
    const auto *slh0_801 = buffer.data(slh0 + 801);
    const auto *slh0_803 = buffer.data(slh0 + 803);
    const auto *slh0_804 = buffer.data(slh0 + 804);
    const auto *slh0_807 = buffer.data(slh0 + 807);
    const auto *slh0_813 = buffer.data(slh0 + 813);

    const auto *slg_375 = buffer.data(slg + 375);
    const auto *slg_378 = buffer.data(slg + 378);
    const auto *slg_385 = buffer.data(slg + 385);
    const auto *slg_389 = buffer.data(slg + 389);
    const auto *slg_390 = buffer.data(slg + 390);
    const auto *slg_392 = buffer.data(slg + 392);
    const auto *slg_393 = buffer.data(slg + 393);
    const auto *slg_395 = buffer.data(slg + 395);
    const auto *slg_400 = buffer.data(slg + 400);
    const auto *slg_402 = buffer.data(slg + 402);
    const auto *slg_403 = buffer.data(slg + 403);
    const auto *slg_404 = buffer.data(slg + 404);
    const auto *slg_405 = buffer.data(slg + 405);
    const auto *slg_406 = buffer.data(slg + 406);
    const auto *slg_407 = buffer.data(slg + 407);
    const auto *slg_408 = buffer.data(slg + 408);
    const auto *slg_410 = buffer.data(slg + 410);
    const auto *slg_415 = buffer.data(slg + 415);
    const auto *slg_417 = buffer.data(slg + 417);
    const auto *slg_418 = buffer.data(slg + 418);
    const auto *slg_419 = buffer.data(slg + 419);
    const auto *slg_420 = buffer.data(slg + 420);
    const auto *slg_422 = buffer.data(slg + 422);
    const auto *slg_423 = buffer.data(slg + 423);
    const auto *slg_425 = buffer.data(slg + 425);
    const auto *slg_430 = buffer.data(slg + 430);
    const auto *slg_434 = buffer.data(slg + 434);
    const auto *slg_435 = buffer.data(slg + 435);
    const auto *slg_437 = buffer.data(slg + 437);
    const auto *slg_438 = buffer.data(slg + 438);
    const auto *slg_440 = buffer.data(slg + 440);
    const auto *slg_445 = buffer.data(slg + 445);
    const auto *slg_449 = buffer.data(slg + 449);
    const auto *slg_450 = buffer.data(slg + 450);
    const auto *slg_452 = buffer.data(slg + 452);
    const auto *slg_455 = buffer.data(slg + 455);
    const auto *slg_498 = buffer.data(slg + 498);
    const auto *slg_500 = buffer.data(slg + 500);
    const auto *slg_501 = buffer.data(slg + 501);
    const auto *slg_504 = buffer.data(slg + 504);
    const auto *slg_505 = buffer.data(slg + 505);
    const auto *slg_506 = buffer.data(slg + 506);
    const auto *slg_507 = buffer.data(slg + 507);
    const auto *slg_508 = buffer.data(slg + 508);
    const auto *slg_509 = buffer.data(slg + 509);
    const auto *slg_520 = buffer.data(slg + 520);
    const auto *slg_521 = buffer.data(slg + 521);
    const auto *slg_522 = buffer.data(slg + 522);
    const auto *slg_523 = buffer.data(slg + 523);
    const auto *slg_524 = buffer.data(slg + 524);
    const auto *slg_525 = buffer.data(slg + 525);
    const auto *slg_528 = buffer.data(slg + 528);
    const auto *slg_530 = buffer.data(slg + 530);
    const auto *slg_531 = buffer.data(slg + 531);
    const auto *slg_534 = buffer.data(slg + 534);
    const auto *slg_535 = buffer.data(slg + 535);
    const auto *slg_536 = buffer.data(slg + 536);
    const auto *slg_537 = buffer.data(slg + 537);
    const auto *slg_538 = buffer.data(slg + 538);
    const auto *slg_539 = buffer.data(slg + 539);
    const auto *slg_540 = buffer.data(slg + 540);
    const auto *slg_543 = buffer.data(slg + 543);
    const auto *slg_545 = buffer.data(slg + 545);
    const auto *slg_546 = buffer.data(slg + 546);
    const auto *slg_549 = buffer.data(slg + 549);
    const auto *slg_550 = buffer.data(slg + 550);
    const auto *slg_551 = buffer.data(slg + 551);
    const auto *slg_552 = buffer.data(slg + 552);
    const auto *slg_553 = buffer.data(slg + 553);
    const auto *slg_554 = buffer.data(slg + 554);
    const auto *slg_560 = buffer.data(slg + 560);
    const auto *slg_564 = buffer.data(slg + 564);
    const auto *slg_565 = buffer.data(slg + 565);
    const auto *slg_566 = buffer.data(slg + 566);
    const auto *slg_567 = buffer.data(slg + 567);
    const auto *slg_568 = buffer.data(slg + 568);
    const auto *slg_569 = buffer.data(slg + 569);
    const auto *slg_570 = buffer.data(slg + 570);
    const auto *slg_573 = buffer.data(slg + 573);
    const auto *slg_575 = buffer.data(slg + 575);
    const auto *slg_576 = buffer.data(slg + 576);
    const auto *slg_579 = buffer.data(slg + 579);
    const auto *slg_580 = buffer.data(slg + 580);
    const auto *slg_581 = buffer.data(slg + 581);
    const auto *slg_582 = buffer.data(slg + 582);
    const auto *slg_583 = buffer.data(slg + 583);
    const auto *slg_584 = buffer.data(slg + 584);

    const auto *slh1_567 = buffer.data(slh1 + 567);
    const auto *slh1_570 = buffer.data(slh1 + 570);
    const auto *slh1_572 = buffer.data(slh1 + 572);
    const auto *slh1_573 = buffer.data(slh1 + 573);
    const auto *slh1_576 = buffer.data(slh1 + 576);
    const auto *slh1_587 = buffer.data(slh1 + 587);
    const auto *slh1_588 = buffer.data(slh1 + 588);
    const auto *slh1_591 = buffer.data(slh1 + 591);
    const auto *slh1_594 = buffer.data(slh1 + 594);
    const auto *slh1_756 = buffer.data(slh1 + 756);
    const auto *slh1_759 = buffer.data(slh1 + 759);
    const auto *slh1_761 = buffer.data(slh1 + 761);
    const auto *slh1_762 = buffer.data(slh1 + 762);
    const auto *slh1_765 = buffer.data(slh1 + 765);
    const auto *slh1_771 = buffer.data(slh1 + 771);
    const auto *slh1_773 = buffer.data(slh1 + 773);
    const auto *slh1_774 = buffer.data(slh1 + 774);
    const auto *slh1_776 = buffer.data(slh1 + 776);
    const auto *slh1_782 = buffer.data(slh1 + 782);
    const auto *slh1_786 = buffer.data(slh1 + 786);
    const auto *slh1_792 = buffer.data(slh1 + 792);
    const auto *slh1_794 = buffer.data(slh1 + 794);
    const auto *slh1_795 = buffer.data(slh1 + 795);
    const auto *slh1_797 = buffer.data(slh1 + 797);
    const auto *slh1_798 = buffer.data(slh1 + 798);
    const auto *slh1_801 = buffer.data(slh1 + 801);
    const auto *slh1_803 = buffer.data(slh1 + 803);
    const auto *slh1_804 = buffer.data(slh1 + 804);
    const auto *slh1_807 = buffer.data(slh1 + 807);
    const auto *slh1_813 = buffer.data(slh1 + 813);

    const auto *smf0_333 = buffer.data(smf0 + 333);
    const auto *smf0_335 = buffer.data(smf0 + 335);
    const auto *smf0_336 = buffer.data(smf0 + 336);
    const auto *smf0_338 = buffer.data(smf0 + 338);
    const auto *smf0_339 = buffer.data(smf0 + 339);
    const auto *smf0_346 = buffer.data(smf0 + 346);
    const auto *smf0_348 = buffer.data(smf0 + 348);
    const auto *smf0_349 = buffer.data(smf0 + 349);
    const auto *smf0_350 = buffer.data(smf0 + 350);
    const auto *smf0_353 = buffer.data(smf0 + 353);
    const auto *smf0_355 = buffer.data(smf0 + 355);
    const auto *smf0_356 = buffer.data(smf0 + 356);
    const auto *smf0_358 = buffer.data(smf0 + 358);
    const auto *smf0_359 = buffer.data(smf0 + 359);

    const auto *smf1_333 = buffer.data(smf1 + 333);
    const auto *smf1_335 = buffer.data(smf1 + 335);
    const auto *smf1_336 = buffer.data(smf1 + 336);
    const auto *smf1_338 = buffer.data(smf1 + 338);
    const auto *smf1_339 = buffer.data(smf1 + 339);
    const auto *smf1_346 = buffer.data(smf1 + 346);
    const auto *smf1_348 = buffer.data(smf1 + 348);
    const auto *smf1_349 = buffer.data(smf1 + 349);
    const auto *smf1_350 = buffer.data(smf1 + 350);
    const auto *smf1_353 = buffer.data(smf1 + 353);
    const auto *smf1_355 = buffer.data(smf1 + 355);
    const auto *smf1_356 = buffer.data(smf1 + 356);
    const auto *smf1_358 = buffer.data(smf1 + 358);
    const auto *smf1_359 = buffer.data(smf1 + 359);

    const auto *smg_495 = buffer.data(smg + 495);
    const auto *smg_497 = buffer.data(smg + 497);
    const auto *smg_498 = buffer.data(smg + 498);
    const auto *smg_500 = buffer.data(smg + 500);
    const auto *smg_501 = buffer.data(smg + 501);
    const auto *smg_504 = buffer.data(smg + 504);
    const auto *smg_505 = buffer.data(smg + 505);
    const auto *smg_506 = buffer.data(smg + 506);
    const auto *smg_507 = buffer.data(smg + 507);
    const auto *smg_508 = buffer.data(smg + 508);
    const auto *smg_509 = buffer.data(smg + 509);
    const auto *smg_510 = buffer.data(smg + 510);
    const auto *smg_512 = buffer.data(smg + 512);
    const auto *smg_513 = buffer.data(smg + 513);
    const auto *smg_515 = buffer.data(smg + 515);
    const auto *smg_520 = buffer.data(smg + 520);
    const auto *smg_521 = buffer.data(smg + 521);
    const auto *smg_522 = buffer.data(smg + 522);
    const auto *smg_523 = buffer.data(smg + 523);
    const auto *smg_524 = buffer.data(smg + 524);
    const auto *smg_525 = buffer.data(smg + 525);
    const auto *smg_527 = buffer.data(smg + 527);
    const auto *smg_528 = buffer.data(smg + 528);
    const auto *smg_530 = buffer.data(smg + 530);
    const auto *smg_531 = buffer.data(smg + 531);
    const auto *smg_534 = buffer.data(smg + 534);
    const auto *smg_535 = buffer.data(smg + 535);
    const auto *smg_536 = buffer.data(smg + 536);
    const auto *smg_537 = buffer.data(smg + 537);
    const auto *smg_538 = buffer.data(smg + 538);
    const auto *smg_539 = buffer.data(smg + 539);
    const auto *smg_540 = buffer.data(smg + 540);
    const auto *smg_542 = buffer.data(smg + 542);
    const auto *smg_543 = buffer.data(smg + 543);
    const auto *smg_545 = buffer.data(smg + 545);
    const auto *smg_550 = buffer.data(smg + 550);
    const auto *smg_551 = buffer.data(smg + 551);
    const auto *smg_552 = buffer.data(smg + 552);
    const auto *smg_553 = buffer.data(smg + 553);
    const auto *smg_554 = buffer.data(smg + 554);
    const auto *smg_555 = buffer.data(smg + 555);
    const auto *smg_557 = buffer.data(smg + 557);
    const auto *smg_558 = buffer.data(smg + 558);
    const auto *smg_560 = buffer.data(smg + 560);
    const auto *smg_565 = buffer.data(smg + 565);
    const auto *smg_566 = buffer.data(smg + 566);
    const auto *smg_567 = buffer.data(smg + 567);
    const auto *smg_568 = buffer.data(smg + 568);
    const auto *smg_569 = buffer.data(smg + 569);
    const auto *smg_570 = buffer.data(smg + 570);
    const auto *smg_572 = buffer.data(smg + 572);
    const auto *smg_573 = buffer.data(smg + 573);
    const auto *smg_575 = buffer.data(smg + 575);
    const auto *smg_580 = buffer.data(smg + 580);
    const auto *smg_581 = buffer.data(smg + 581);
    const auto *smg_582 = buffer.data(smg + 582);
    const auto *smg_583 = buffer.data(smg + 583);
    const auto *smg_584 = buffer.data(smg + 584);

#pragma omp simd aligned(t_694, t_695, t_696, t_697, pc_x, pc_y, pc_z, slg_375, slg_390, \
                         slg_392, slg_498, smf0_333, smf1_333, smg_495, smg_497, \
                         smg_498 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_694[k] = f_10 * slg_390[k]
                   + f_3 * pc_y[k] * smg_495[k];

        t_695[k] = f_15 * slg_375[k]
                   + f_3 * pc_z[k] * smg_495[k];

        t_696[k] = f_10 * slg_498[k]
                   + f_4 * smf0_333[k]
                   - f_5 * smf1_333[k]
                   + f_3 * pc_x[k] * smg_498[k];

        t_697[k] = f_10 * slg_392[k]
                   + f_3 * pc_y[k] * smg_497[k];
    }

#pragma omp simd aligned(t_698, t_699, t_700, pc_x, pc_z, slg_378, slg_500, slg_501, smf0_335, \
                         smf0_336, smf1_335, smf1_336, smg_498, smg_500, \
                         smg_501 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_698[k] = f_10 * slg_500[k]
                   + f_4 * smf0_335[k]
                   - f_5 * smf1_335[k]
                   + f_3 * pc_x[k] * smg_500[k];

        t_699[k] = f_10 * slg_501[k]
                   + f_6 * smf0_336[k]
                   - f_7 * smf1_336[k]
                   + f_3 * pc_x[k] * smg_501[k];

        t_700[k] = f_15 * slg_378[k]
                   + f_3 * pc_z[k] * smg_498[k];
    }

#pragma omp simd aligned(t_701, t_702, t_703, t_704, pc_x, pc_y, slg_395, slg_504, slg_505, \
                         slg_506, smf0_339, smf1_339, smg_500, smg_504, smg_505, \
                         smg_506 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_701[k] = f_10 * slg_395[k]
                   + f_3 * pc_y[k] * smg_500[k];

        t_702[k] = f_10 * slg_504[k]
                   + f_6 * smf0_339[k]
                   - f_7 * smf1_339[k]
                   + f_3 * pc_x[k] * smg_504[k];

        t_703[k] = f_10 * slg_505[k]
                   + f_3 * pc_x[k] * smg_505[k];

        t_704[k] = f_10 * slg_506[k]
                   + f_3 * pc_x[k] * smg_506[k];
    }

#pragma omp simd aligned(t_705, t_706, t_707, t_708, pc_x, pc_y, slg_400, slg_507, slg_508, \
                         slg_509, smf0_336, smf1_336, smg_505, smg_507, smg_508, \
                         smg_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_705[k] = f_10 * slg_507[k]
                   + f_3 * pc_x[k] * smg_507[k];

        t_706[k] = f_10 * slg_508[k]
                   + f_3 * pc_x[k] * smg_508[k];

        t_707[k] = f_10 * slg_509[k]
                   + f_3 * pc_x[k] * smg_509[k];

        t_708[k] = f_10 * slg_400[k]
                   + f_1 * smf0_336[k]
                   - f_2 * smf1_336[k]
                   + f_3 * pc_y[k] * smg_505[k];
    }

#pragma omp simd aligned(t_709, t_710, t_711, pc_y, pc_z, slg_385, slg_402, slg_403, smf0_338, \
                         smf0_339, smf1_338, smf1_339, smg_505, smg_507, \
                         smg_508 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_709[k] = f_15 * slg_385[k]
                   + f_3 * pc_z[k] * smg_505[k];

        t_710[k] = f_10 * slg_402[k]
                   + f_4 * smf0_338[k]
                   - f_5 * smf1_338[k]
                   + f_3 * pc_y[k] * smg_507[k];

        t_711[k] = f_10 * slg_403[k]
                   + f_6 * smf0_339[k]
                   - f_7 * smf1_339[k]
                   + f_3 * pc_y[k] * smg_508[k];
    }

#pragma omp simd aligned(t_712, t_713, t_714, t_715, pb_y, pc_y, pc_z, slh0_567, slg_389, \
                         slg_404, slg_405, slh1_567, smf0_339, smf1_339, smg_509, \
                         smg_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_712[k] = f_10 * slg_404[k]
                   + f_3 * pc_y[k] * smg_509[k];

        t_713[k] = f_15 * slg_389[k]
                   + f_1 * smf0_339[k]
                   - f_2 * smf1_339[k]
                   + f_3 * pc_z[k] * smg_509[k];

        t_714[k] = pb_y[k] * slh0_567[k]
                   - f_8 * pc_y[k] * slh1_567[k];

        t_715[k] = f_9 * slg_405[k]
                   + f_3 * pc_y[k] * smg_510[k];
    }

#pragma omp simd aligned(t_716, t_717, t_718, t_719, pb_y, pc_y, pc_z, slh0_570, slh0_572, \
                         slg_390, slg_406, slg_407, slh1_570, slh1_572, smg_510, \
                         smg_512 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_716[k] = f_14 * slg_390[k]
                   + f_3 * pc_z[k] * smg_510[k];

        t_717[k] = pb_y[k] * slh0_570[k]
                   + f_10 * slg_406[k]
                   - f_8 * pc_y[k] * slh1_570[k];

        t_718[k] = f_9 * slg_407[k]
                   + f_3 * pc_y[k] * smg_512[k];

        t_719[k] = pb_y[k] * slh0_572[k]
                   - f_8 * pc_y[k] * slh1_572[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, t_723, pb_y, pc_y, pc_z, slh0_573, slh0_576, \
                         slg_393, slg_408, slg_410, slh1_573, slh1_576, smg_513, \
                         smg_515 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = pb_y[k] * slh0_573[k]
                   + f_11 * slg_408[k]
                   - f_8 * pc_y[k] * slh1_573[k];

        t_721[k] = f_14 * slg_393[k]
                   + f_3 * pc_z[k] * smg_513[k];

        t_722[k] = f_9 * slg_410[k]
                   + f_3 * pc_y[k] * smg_515[k];

        t_723[k] = pb_y[k] * slh0_576[k]
                   - f_8 * pc_y[k] * slh1_576[k];
    }

#pragma omp simd aligned(t_724, t_725, t_726, t_727, t_728, pc_x, slg_520, slg_521, slg_522, \
                         slg_523, slg_524, smg_520, smg_521, smg_522, smg_523, \
                         smg_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_724[k] = f_10 * slg_520[k]
                   + f_3 * pc_x[k] * smg_520[k];

        t_725[k] = f_10 * slg_521[k]
                   + f_3 * pc_x[k] * smg_521[k];

        t_726[k] = f_10 * slg_522[k]
                   + f_3 * pc_x[k] * smg_522[k];

        t_727[k] = f_10 * slg_523[k]
                   + f_3 * pc_x[k] * smg_523[k];

        t_728[k] = f_10 * slg_524[k]
                   + f_3 * pc_x[k] * smg_524[k];
    }

#pragma omp simd aligned(t_729, t_730, t_731, pc_y, pc_z, slg_400, slg_415, slg_417, smf0_346, \
                         smf0_348, smf1_346, smf1_348, smg_520, \
                         smg_522 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_729[k] = f_9 * slg_415[k]
                   + f_1 * smf0_346[k]
                   - f_2 * smf1_346[k]
                   + f_3 * pc_y[k] * smg_520[k];

        t_730[k] = f_14 * slg_400[k]
                   + f_3 * pc_z[k] * smg_520[k];

        t_731[k] = f_9 * slg_417[k]
                   + f_4 * smf0_348[k]
                   - f_5 * smf1_348[k]
                   + f_3 * pc_y[k] * smg_522[k];
    }

#pragma omp simd aligned(t_732, t_733, t_734, pb_y, pc_y, slh0_587, slg_418, slg_419, \
                         slh1_587, smf0_349, smf1_349, smg_523, \
                         smg_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_732[k] = f_9 * slg_418[k]
                   + f_6 * smf0_349[k]
                   - f_7 * smf1_349[k]
                   + f_3 * pc_y[k] * smg_523[k];

        t_733[k] = f_9 * slg_419[k]
                   + f_3 * pc_y[k] * smg_524[k];

        t_734[k] = pb_y[k] * slh0_587[k]
                   - f_8 * pc_y[k] * slh1_587[k];
    }

#pragma omp simd aligned(t_735, t_736, t_737, t_738, pc_x, pc_y, pc_z, slg_405, slg_525, \
                         slg_528, smf0_350, smf0_353, smf1_350, smf1_353, smg_525, \
                         smg_528 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_735[k] = f_10 * slg_525[k]
                   + f_1 * smf0_350[k]
                   - f_2 * smf1_350[k]
                   + f_3 * pc_x[k] * smg_525[k];

        t_736[k] = f_3 * pc_y[k] * smg_525[k];

        t_737[k] = f_13 * slg_405[k]
                   + f_3 * pc_z[k] * smg_525[k];

        t_738[k] = f_10 * slg_528[k]
                   + f_4 * smf0_353[k]
                   - f_5 * smf1_353[k]
                   + f_3 * pc_x[k] * smg_528[k];
    }

#pragma omp simd aligned(t_739, t_740, t_741, pc_x, pc_y, slg_530, slg_531, smf0_355, \
                         smf0_356, smf1_355, smf1_356, smg_527, smg_530, \
                         smg_531 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_739[k] = f_3 * pc_y[k] * smg_527[k];

        t_740[k] = f_10 * slg_530[k]
                   + f_4 * smf0_355[k]
                   - f_5 * smf1_355[k]
                   + f_3 * pc_x[k] * smg_530[k];

        t_741[k] = f_10 * slg_531[k]
                   + f_6 * smf0_356[k]
                   - f_7 * smf1_356[k]
                   + f_3 * pc_x[k] * smg_531[k];
    }

#pragma omp simd aligned(t_742, t_743, t_744, t_745, pc_x, pc_y, pc_z, slg_408, slg_534, \
                         slg_535, smf0_359, smf1_359, smg_528, smg_530, smg_534, \
                         smg_535 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_742[k] = f_13 * slg_408[k]
                   + f_3 * pc_z[k] * smg_528[k];

        t_743[k] = f_3 * pc_y[k] * smg_530[k];

        t_744[k] = f_10 * slg_534[k]
                   + f_6 * smf0_359[k]
                   - f_7 * smf1_359[k]
                   + f_3 * pc_x[k] * smg_534[k];

        t_745[k] = f_10 * slg_535[k]
                   + f_3 * pc_x[k] * smg_535[k];
    }

#pragma omp simd aligned(t_746, t_747, t_748, t_749, pc_x, slg_536, slg_537, slg_538, slg_539, \
                         smg_536, smg_537, smg_538, smg_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_746[k] = f_10 * slg_536[k]
                   + f_3 * pc_x[k] * smg_536[k];

        t_747[k] = f_10 * slg_537[k]
                   + f_3 * pc_x[k] * smg_537[k];

        t_748[k] = f_10 * slg_538[k]
                   + f_3 * pc_x[k] * smg_538[k];

        t_749[k] = f_10 * slg_539[k]
                   + f_3 * pc_x[k] * smg_539[k];
    }

#pragma omp simd aligned(t_750, t_751, t_752, t_753, pc_y, pc_z, slg_415, smf0_356, smf0_358, \
                         smf0_359, smf1_356, smf1_358, smf1_359, smg_535, smg_537, \
                         smg_538 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_750[k] = f_1 * smf0_356[k]
                   - f_2 * smf1_356[k]
                   + f_3 * pc_y[k] * smg_535[k];

        t_751[k] = f_13 * slg_415[k]
                   + f_3 * pc_z[k] * smg_535[k];

        t_752[k] = f_4 * smf0_358[k]
                   - f_5 * smf1_358[k]
                   + f_3 * pc_y[k] * smg_537[k];

        t_753[k] = f_6 * smf0_359[k]
                   - f_7 * smf1_359[k]
                   + f_3 * pc_y[k] * smg_538[k];
    }

#pragma omp simd aligned(t_754, t_755, t_756, pb_x, pc_x, pc_y, pc_z, slh0_756, slg_419, \
                         slg_540, slh1_756, smf0_359, smf1_359, \
                         smg_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_754[k] = f_3 * pc_y[k] * smg_539[k];

        t_755[k] = f_13 * slg_419[k]
                   + f_1 * smf0_359[k]
                   - f_2 * smf1_359[k]
                   + f_3 * pc_z[k] * smg_539[k];

        t_756[k] = pb_x[k] * slh0_756[k]
                   + f_15 * slg_540[k]
                   - f_8 * pc_x[k] * slh1_756[k];
    }

#pragma omp simd aligned(t_757, t_758, t_759, t_760, pb_x, pc_x, pc_y, pc_z, slh0_759, \
                         slg_420, slg_422, slg_543, slh1_759, smg_540, \
                         smg_542 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_757[k] = f_12 * slg_420[k]
                   + f_3 * pc_y[k] * smg_540[k];

        t_758[k] = f_3 * pc_z[k] * smg_540[k];

        t_759[k] = pb_x[k] * slh0_759[k]
                   + f_11 * slg_543[k]
                   - f_8 * pc_x[k] * slh1_759[k];

        t_760[k] = f_12 * slg_422[k]
                   + f_3 * pc_y[k] * smg_542[k];
    }

#pragma omp simd aligned(t_761, t_762, t_763, pb_x, pc_x, pc_z, slh0_761, slh0_762, slg_545, \
                         slg_546, slh1_761, slh1_762, smg_543 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_761[k] = pb_x[k] * slh0_761[k]
                   + f_11 * slg_545[k]
                   - f_8 * pc_x[k] * slh1_761[k];

        t_762[k] = pb_x[k] * slh0_762[k]
                   + f_10 * slg_546[k]
                   - f_8 * pc_x[k] * slh1_762[k];

        t_763[k] = f_3 * pc_z[k] * smg_543[k];
    }

#pragma omp simd aligned(t_764, t_765, t_766, t_767, pb_x, pc_x, pc_y, slh0_765, slg_425, \
                         slg_549, slg_550, slg_551, slh1_765, smg_545, smg_550, \
                         smg_551 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_764[k] = f_12 * slg_425[k]
                   + f_3 * pc_y[k] * smg_545[k];

        t_765[k] = pb_x[k] * slh0_765[k]
                   + f_10 * slg_549[k]
                   - f_8 * pc_x[k] * slh1_765[k];

        t_766[k] = f_9 * slg_550[k]
                   + f_3 * pc_x[k] * smg_550[k];

        t_767[k] = f_9 * slg_551[k]
                   + f_3 * pc_x[k] * smg_551[k];
    }

#pragma omp simd aligned(t_768, t_769, t_770, t_771, pb_x, pc_x, slh0_771, slg_552, slg_553, \
                         slg_554, slh1_771, smg_552, smg_553, smg_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_768[k] = f_9 * slg_552[k]
                   + f_3 * pc_x[k] * smg_552[k];

        t_769[k] = f_9 * slg_553[k]
                   + f_3 * pc_x[k] * smg_553[k];

        t_770[k] = f_9 * slg_554[k]
                   + f_3 * pc_x[k] * smg_554[k];

        t_771[k] = pb_x[k] * slh0_771[k]
                   - f_8 * pc_x[k] * slh1_771[k];
    }

#pragma omp simd aligned(t_772, t_773, t_774, t_775, pb_x, pc_x, pc_y, pc_z, slh0_773, \
                         slh0_774, slg_434, slh1_773, slh1_774, smg_550, \
                         smg_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_772[k] = f_3 * pc_z[k] * smg_550[k];

        t_773[k] = pb_x[k] * slh0_773[k]
                   - f_8 * pc_x[k] * slh1_773[k];

        t_774[k] = pb_x[k] * slh0_774[k]
                   - f_8 * pc_x[k] * slh1_774[k];

        t_775[k] = f_12 * slg_434[k]
                   + f_3 * pc_y[k] * smg_554[k];
    }

#pragma omp simd aligned(t_776, t_777, t_778, t_779, pb_x, pb_z, pc_x, pc_y, pc_z, slh0_588, \
                         slh0_776, slg_420, slg_435, slh1_588, slh1_776, \
                         smg_555 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_776[k] = pb_x[k] * slh0_776[k]
                   - f_8 * pc_x[k] * slh1_776[k];

        t_777[k] = pb_z[k] * slh0_588[k]
                   - f_8 * pc_z[k] * slh1_588[k];

        t_778[k] = f_13 * slg_435[k]
                   + f_3 * pc_y[k] * smg_555[k];

        t_779[k] = f_9 * slg_420[k]
                   + f_3 * pc_z[k] * smg_555[k];
    }

#pragma omp simd aligned(t_780, t_781, t_782, pb_x, pb_z, pc_x, pc_y, pc_z, slh0_591, \
                         slh0_782, slg_437, slg_560, slh1_591, slh1_782, \
                         smg_557 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_780[k] = pb_z[k] * slh0_591[k]
                   - f_8 * pc_z[k] * slh1_591[k];

        t_781[k] = f_13 * slg_437[k]
                   + f_3 * pc_y[k] * smg_557[k];

        t_782[k] = pb_x[k] * slh0_782[k]
                   + f_11 * slg_560[k]
                   - f_8 * pc_x[k] * slh1_782[k];
    }

#pragma omp simd aligned(t_783, t_784, t_785, pb_z, pc_y, pc_z, slh0_594, slg_423, slg_440, \
                         slh1_594, smg_558, smg_560 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_783[k] = pb_z[k] * slh0_594[k]
                   - f_8 * pc_z[k] * slh1_594[k];

        t_784[k] = f_9 * slg_423[k]
                   + f_3 * pc_z[k] * smg_558[k];

        t_785[k] = f_13 * slg_440[k]
                   + f_3 * pc_y[k] * smg_560[k];
    }

#pragma omp simd aligned(t_786, t_787, t_788, t_789, pb_x, pc_x, slh0_786, slg_564, slg_565, \
                         slg_566, slg_567, slh1_786, smg_565, smg_566, \
                         smg_567 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_786[k] = pb_x[k] * slh0_786[k]
                   + f_10 * slg_564[k]
                   - f_8 * pc_x[k] * slh1_786[k];

        t_787[k] = f_9 * slg_565[k]
                   + f_3 * pc_x[k] * smg_565[k];

        t_788[k] = f_9 * slg_566[k]
                   + f_3 * pc_x[k] * smg_566[k];

        t_789[k] = f_9 * slg_567[k]
                   + f_3 * pc_x[k] * smg_567[k];
    }

#pragma omp simd aligned(t_790, t_791, t_792, t_793, pb_x, pc_x, pc_z, slh0_792, slg_430, \
                         slg_568, slg_569, slh1_792, smg_565, smg_568, \
                         smg_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_790[k] = f_9 * slg_568[k]
                   + f_3 * pc_x[k] * smg_568[k];

        t_791[k] = f_9 * slg_569[k]
                   + f_3 * pc_x[k] * smg_569[k];

        t_792[k] = pb_x[k] * slh0_792[k]
                   - f_8 * pc_x[k] * slh1_792[k];

        t_793[k] = f_9 * slg_430[k]
                   + f_3 * pc_z[k] * smg_565[k];
    }

#pragma omp simd aligned(t_794, t_795, t_796, t_797, pb_x, pc_x, pc_y, slh0_794, slh0_795, \
                         slh0_797, slg_449, slh1_794, slh1_795, slh1_797, \
                         smg_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_794[k] = pb_x[k] * slh0_794[k]
                   - f_8 * pc_x[k] * slh1_794[k];

        t_795[k] = pb_x[k] * slh0_795[k]
                   - f_8 * pc_x[k] * slh1_795[k];

        t_796[k] = f_13 * slg_449[k]
                   + f_3 * pc_y[k] * smg_569[k];

        t_797[k] = pb_x[k] * slh0_797[k]
                   - f_8 * pc_x[k] * slh1_797[k];
    }

#pragma omp simd aligned(t_798, t_799, t_800, pb_x, pc_x, pc_y, pc_z, slh0_798, slg_435, \
                         slg_450, slg_570, slh1_798, smg_570 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_798[k] = pb_x[k] * slh0_798[k]
                   + f_15 * slg_570[k]
                   - f_8 * pc_x[k] * slh1_798[k];

        t_799[k] = f_14 * slg_450[k]
                   + f_3 * pc_y[k] * smg_570[k];

        t_800[k] = f_10 * slg_435[k]
                   + f_3 * pc_z[k] * smg_570[k];
    }

#pragma omp simd aligned(t_801, t_802, t_803, pb_x, pc_x, pc_y, slh0_801, slh0_803, slg_452, \
                         slg_573, slg_575, slh1_801, slh1_803, \
                         smg_572 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_801[k] = pb_x[k] * slh0_801[k]
                   + f_11 * slg_573[k]
                   - f_8 * pc_x[k] * slh1_801[k];

        t_802[k] = f_14 * slg_452[k]
                   + f_3 * pc_y[k] * smg_572[k];

        t_803[k] = pb_x[k] * slh0_803[k]
                   + f_11 * slg_575[k]
                   - f_8 * pc_x[k] * slh1_803[k];
    }

#pragma omp simd aligned(t_804, t_805, t_806, pb_x, pc_x, pc_y, pc_z, slh0_804, slg_438, \
                         slg_455, slg_576, slh1_804, smg_573, smg_575 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_804[k] = pb_x[k] * slh0_804[k]
                   + f_10 * slg_576[k]
                   - f_8 * pc_x[k] * slh1_804[k];

        t_805[k] = f_10 * slg_438[k]
                   + f_3 * pc_z[k] * smg_573[k];

        t_806[k] = f_14 * slg_455[k]
                   + f_3 * pc_y[k] * smg_575[k];
    }

#pragma omp simd aligned(t_807, t_808, t_809, t_810, pb_x, pc_x, slh0_807, slg_579, slg_580, \
                         slg_581, slg_582, slh1_807, smg_580, smg_581, \
                         smg_582 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_807[k] = pb_x[k] * slh0_807[k]
                   + f_10 * slg_579[k]
                   - f_8 * pc_x[k] * slh1_807[k];

        t_808[k] = f_9 * slg_580[k]
                   + f_3 * pc_x[k] * smg_580[k];

        t_809[k] = f_9 * slg_581[k]
                   + f_3 * pc_x[k] * smg_581[k];

        t_810[k] = f_9 * slg_582[k]
                   + f_3 * pc_x[k] * smg_582[k];
    }

#pragma omp simd aligned(t_811, t_812, t_813, t_814, pb_x, pc_x, pc_z, slh0_813, slg_445, \
                         slg_583, slg_584, slh1_813, smg_580, smg_583, \
                         smg_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_811[k] = f_9 * slg_583[k]
                   + f_3 * pc_x[k] * smg_583[k];

        t_812[k] = f_9 * slg_584[k]
                   + f_3 * pc_x[k] * smg_584[k];

        t_813[k] = pb_x[k] * slh0_813[k]
                   - f_8 * pc_x[k] * slh1_813[k];

        t_814[k] = f_10 * slg_445[k]
                   + f_3 * pc_z[k] * smg_580[k];
    }
}

static auto
compute_prim_smh_three_center_electron_repulsion_0_piece7(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t slh0,
                                                          const size_t slg, const size_t slh1,
                                                          const size_t smg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_3 = p / q;
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_12 = 4.0 / q;
    const auto f_13 = 3.5 / q;
    const auto f_14 = 3.0 / q;
    const auto f_15 = 2.5 / q;
    const auto f_16 = 2.0 / q;

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
    auto *t_934 = buffer.data(target + 934);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *slh0_735 = buffer.data(slh0 + 735);
    const auto *slh0_740 = buffer.data(slh0 + 740);
    const auto *slh0_744 = buffer.data(slh0 + 744);
    const auto *slh0_815 = buffer.data(slh0 + 815);
    const auto *slh0_816 = buffer.data(slh0 + 816);
    const auto *slh0_818 = buffer.data(slh0 + 818);
    const auto *slh0_819 = buffer.data(slh0 + 819);
    const auto *slh0_822 = buffer.data(slh0 + 822);
    const auto *slh0_824 = buffer.data(slh0 + 824);
    const auto *slh0_825 = buffer.data(slh0 + 825);
    const auto *slh0_828 = buffer.data(slh0 + 828);
    const auto *slh0_834 = buffer.data(slh0 + 834);
    const auto *slh0_836 = buffer.data(slh0 + 836);
    const auto *slh0_837 = buffer.data(slh0 + 837);
    const auto *slh0_839 = buffer.data(slh0 + 839);
    const auto *slh0_840 = buffer.data(slh0 + 840);
    const auto *slh0_843 = buffer.data(slh0 + 843);
    const auto *slh0_845 = buffer.data(slh0 + 845);
    const auto *slh0_846 = buffer.data(slh0 + 846);
    const auto *slh0_849 = buffer.data(slh0 + 849);
    const auto *slh0_855 = buffer.data(slh0 + 855);
    const auto *slh0_857 = buffer.data(slh0 + 857);
    const auto *slh0_858 = buffer.data(slh0 + 858);
    const auto *slh0_860 = buffer.data(slh0 + 860);
    const auto *slh0_861 = buffer.data(slh0 + 861);
    const auto *slh0_864 = buffer.data(slh0 + 864);
    const auto *slh0_866 = buffer.data(slh0 + 866);
    const auto *slh0_867 = buffer.data(slh0 + 867);
    const auto *slh0_870 = buffer.data(slh0 + 870);
    const auto *slh0_876 = buffer.data(slh0 + 876);
    const auto *slh0_878 = buffer.data(slh0 + 878);
    const auto *slh0_879 = buffer.data(slh0 + 879);
    const auto *slh0_881 = buffer.data(slh0 + 881);
    const auto *slh0_882 = buffer.data(slh0 + 882);
    const auto *slh0_885 = buffer.data(slh0 + 885);
    const auto *slh0_887 = buffer.data(slh0 + 887);
    const auto *slh0_888 = buffer.data(slh0 + 888);
    const auto *slh0_891 = buffer.data(slh0 + 891);
    const auto *slh0_897 = buffer.data(slh0 + 897);
    const auto *slh0_899 = buffer.data(slh0 + 899);
    const auto *slh0_900 = buffer.data(slh0 + 900);
    const auto *slh0_902 = buffer.data(slh0 + 902);
    const auto *slh0_906 = buffer.data(slh0 + 906);
    const auto *slh0_909 = buffer.data(slh0 + 909);
    const auto *slh0_918 = buffer.data(slh0 + 918);
    const auto *slh0_920 = buffer.data(slh0 + 920);
    const auto *slh0_921 = buffer.data(slh0 + 921);
    const auto *slh0_923 = buffer.data(slh0 + 923);
    const auto *slh0_924 = buffer.data(slh0 + 924);
    const auto *slh0_927 = buffer.data(slh0 + 927);
    const auto *slh0_929 = buffer.data(slh0 + 929);
    const auto *slh0_930 = buffer.data(slh0 + 930);
    const auto *slh0_933 = buffer.data(slh0 + 933);

    const auto *slg_450 = buffer.data(slg + 450);
    const auto *slg_453 = buffer.data(slg + 453);
    const auto *slg_460 = buffer.data(slg + 460);
    const auto *slg_464 = buffer.data(slg + 464);
    const auto *slg_465 = buffer.data(slg + 465);
    const auto *slg_467 = buffer.data(slg + 467);
    const auto *slg_468 = buffer.data(slg + 468);
    const auto *slg_470 = buffer.data(slg + 470);
    const auto *slg_475 = buffer.data(slg + 475);
    const auto *slg_479 = buffer.data(slg + 479);
    const auto *slg_480 = buffer.data(slg + 480);
    const auto *slg_482 = buffer.data(slg + 482);
    const auto *slg_483 = buffer.data(slg + 483);
    const auto *slg_485 = buffer.data(slg + 485);
    const auto *slg_490 = buffer.data(slg + 490);
    const auto *slg_494 = buffer.data(slg + 494);
    const auto *slg_495 = buffer.data(slg + 495);
    const auto *slg_497 = buffer.data(slg + 497);
    const auto *slg_498 = buffer.data(slg + 498);
    const auto *slg_500 = buffer.data(slg + 500);
    const auto *slg_505 = buffer.data(slg + 505);
    const auto *slg_509 = buffer.data(slg + 509);
    const auto *slg_510 = buffer.data(slg + 510);
    const auto *slg_512 = buffer.data(slg + 512);
    const auto *slg_513 = buffer.data(slg + 513);
    const auto *slg_515 = buffer.data(slg + 515);
    const auto *slg_520 = buffer.data(slg + 520);
    const auto *slg_524 = buffer.data(slg + 524);
    const auto *slg_525 = buffer.data(slg + 525);
    const auto *slg_527 = buffer.data(slg + 527);
    const auto *slg_528 = buffer.data(slg + 528);
    const auto *slg_530 = buffer.data(slg + 530);
    const auto *slg_539 = buffer.data(slg + 539);
    const auto *slg_585 = buffer.data(slg + 585);
    const auto *slg_588 = buffer.data(slg + 588);
    const auto *slg_590 = buffer.data(slg + 590);
    const auto *slg_591 = buffer.data(slg + 591);
    const auto *slg_594 = buffer.data(slg + 594);
    const auto *slg_595 = buffer.data(slg + 595);
    const auto *slg_596 = buffer.data(slg + 596);
    const auto *slg_597 = buffer.data(slg + 597);
    const auto *slg_598 = buffer.data(slg + 598);
    const auto *slg_599 = buffer.data(slg + 599);
    const auto *slg_600 = buffer.data(slg + 600);
    const auto *slg_603 = buffer.data(slg + 603);
    const auto *slg_605 = buffer.data(slg + 605);
    const auto *slg_606 = buffer.data(slg + 606);
    const auto *slg_609 = buffer.data(slg + 609);
    const auto *slg_610 = buffer.data(slg + 610);
    const auto *slg_611 = buffer.data(slg + 611);
    const auto *slg_612 = buffer.data(slg + 612);
    const auto *slg_613 = buffer.data(slg + 613);
    const auto *slg_614 = buffer.data(slg + 614);
    const auto *slg_615 = buffer.data(slg + 615);
    const auto *slg_618 = buffer.data(slg + 618);
    const auto *slg_620 = buffer.data(slg + 620);
    const auto *slg_621 = buffer.data(slg + 621);
    const auto *slg_624 = buffer.data(slg + 624);
    const auto *slg_625 = buffer.data(slg + 625);
    const auto *slg_626 = buffer.data(slg + 626);
    const auto *slg_627 = buffer.data(slg + 627);
    const auto *slg_628 = buffer.data(slg + 628);
    const auto *slg_629 = buffer.data(slg + 629);
    const auto *slg_630 = buffer.data(slg + 630);
    const auto *slg_633 = buffer.data(slg + 633);
    const auto *slg_635 = buffer.data(slg + 635);
    const auto *slg_636 = buffer.data(slg + 636);
    const auto *slg_639 = buffer.data(slg + 639);
    const auto *slg_640 = buffer.data(slg + 640);
    const auto *slg_641 = buffer.data(slg + 641);
    const auto *slg_642 = buffer.data(slg + 642);
    const auto *slg_643 = buffer.data(slg + 643);
    const auto *slg_644 = buffer.data(slg + 644);
    const auto *slg_648 = buffer.data(slg + 648);
    const auto *slg_651 = buffer.data(slg + 651);
    const auto *slg_655 = buffer.data(slg + 655);
    const auto *slg_656 = buffer.data(slg + 656);
    const auto *slg_657 = buffer.data(slg + 657);
    const auto *slg_658 = buffer.data(slg + 658);
    const auto *slg_659 = buffer.data(slg + 659);
    const auto *slg_660 = buffer.data(slg + 660);
    const auto *slg_663 = buffer.data(slg + 663);
    const auto *slg_665 = buffer.data(slg + 665);
    const auto *slg_666 = buffer.data(slg + 666);
    const auto *slg_669 = buffer.data(slg + 669);
    const auto *slg_670 = buffer.data(slg + 670);

    const auto *slh1_735 = buffer.data(slh1 + 735);
    const auto *slh1_740 = buffer.data(slh1 + 740);
    const auto *slh1_744 = buffer.data(slh1 + 744);
    const auto *slh1_815 = buffer.data(slh1 + 815);
    const auto *slh1_816 = buffer.data(slh1 + 816);
    const auto *slh1_818 = buffer.data(slh1 + 818);
    const auto *slh1_819 = buffer.data(slh1 + 819);
    const auto *slh1_822 = buffer.data(slh1 + 822);
    const auto *slh1_824 = buffer.data(slh1 + 824);
    const auto *slh1_825 = buffer.data(slh1 + 825);
    const auto *slh1_828 = buffer.data(slh1 + 828);
    const auto *slh1_834 = buffer.data(slh1 + 834);
    const auto *slh1_836 = buffer.data(slh1 + 836);
    const auto *slh1_837 = buffer.data(slh1 + 837);
    const auto *slh1_839 = buffer.data(slh1 + 839);
    const auto *slh1_840 = buffer.data(slh1 + 840);
    const auto *slh1_843 = buffer.data(slh1 + 843);
    const auto *slh1_845 = buffer.data(slh1 + 845);
    const auto *slh1_846 = buffer.data(slh1 + 846);
    const auto *slh1_849 = buffer.data(slh1 + 849);
    const auto *slh1_855 = buffer.data(slh1 + 855);
    const auto *slh1_857 = buffer.data(slh1 + 857);
    const auto *slh1_858 = buffer.data(slh1 + 858);
    const auto *slh1_860 = buffer.data(slh1 + 860);
    const auto *slh1_861 = buffer.data(slh1 + 861);
    const auto *slh1_864 = buffer.data(slh1 + 864);
    const auto *slh1_866 = buffer.data(slh1 + 866);
    const auto *slh1_867 = buffer.data(slh1 + 867);
    const auto *slh1_870 = buffer.data(slh1 + 870);
    const auto *slh1_876 = buffer.data(slh1 + 876);
    const auto *slh1_878 = buffer.data(slh1 + 878);
    const auto *slh1_879 = buffer.data(slh1 + 879);
    const auto *slh1_881 = buffer.data(slh1 + 881);
    const auto *slh1_882 = buffer.data(slh1 + 882);
    const auto *slh1_885 = buffer.data(slh1 + 885);
    const auto *slh1_887 = buffer.data(slh1 + 887);
    const auto *slh1_888 = buffer.data(slh1 + 888);
    const auto *slh1_891 = buffer.data(slh1 + 891);
    const auto *slh1_897 = buffer.data(slh1 + 897);
    const auto *slh1_899 = buffer.data(slh1 + 899);
    const auto *slh1_900 = buffer.data(slh1 + 900);
    const auto *slh1_902 = buffer.data(slh1 + 902);
    const auto *slh1_906 = buffer.data(slh1 + 906);
    const auto *slh1_909 = buffer.data(slh1 + 909);
    const auto *slh1_918 = buffer.data(slh1 + 918);
    const auto *slh1_920 = buffer.data(slh1 + 920);
    const auto *slh1_921 = buffer.data(slh1 + 921);
    const auto *slh1_923 = buffer.data(slh1 + 923);
    const auto *slh1_924 = buffer.data(slh1 + 924);
    const auto *slh1_927 = buffer.data(slh1 + 927);
    const auto *slh1_929 = buffer.data(slh1 + 929);
    const auto *slh1_930 = buffer.data(slh1 + 930);
    const auto *slh1_933 = buffer.data(slh1 + 933);

    const auto *smg_584 = buffer.data(smg + 584);
    const auto *smg_585 = buffer.data(smg + 585);
    const auto *smg_587 = buffer.data(smg + 587);
    const auto *smg_588 = buffer.data(smg + 588);
    const auto *smg_590 = buffer.data(smg + 590);
    const auto *smg_595 = buffer.data(smg + 595);
    const auto *smg_596 = buffer.data(smg + 596);
    const auto *smg_597 = buffer.data(smg + 597);
    const auto *smg_598 = buffer.data(smg + 598);
    const auto *smg_599 = buffer.data(smg + 599);
    const auto *smg_600 = buffer.data(smg + 600);
    const auto *smg_602 = buffer.data(smg + 602);
    const auto *smg_603 = buffer.data(smg + 603);
    const auto *smg_605 = buffer.data(smg + 605);
    const auto *smg_610 = buffer.data(smg + 610);
    const auto *smg_611 = buffer.data(smg + 611);
    const auto *smg_612 = buffer.data(smg + 612);
    const auto *smg_613 = buffer.data(smg + 613);
    const auto *smg_614 = buffer.data(smg + 614);
    const auto *smg_615 = buffer.data(smg + 615);
    const auto *smg_617 = buffer.data(smg + 617);
    const auto *smg_618 = buffer.data(smg + 618);
    const auto *smg_620 = buffer.data(smg + 620);
    const auto *smg_625 = buffer.data(smg + 625);
    const auto *smg_626 = buffer.data(smg + 626);
    const auto *smg_627 = buffer.data(smg + 627);
    const auto *smg_628 = buffer.data(smg + 628);
    const auto *smg_629 = buffer.data(smg + 629);
    const auto *smg_630 = buffer.data(smg + 630);
    const auto *smg_632 = buffer.data(smg + 632);
    const auto *smg_633 = buffer.data(smg + 633);
    const auto *smg_635 = buffer.data(smg + 635);
    const auto *smg_640 = buffer.data(smg + 640);
    const auto *smg_641 = buffer.data(smg + 641);
    const auto *smg_642 = buffer.data(smg + 642);
    const auto *smg_643 = buffer.data(smg + 643);
    const auto *smg_644 = buffer.data(smg + 644);
    const auto *smg_645 = buffer.data(smg + 645);
    const auto *smg_647 = buffer.data(smg + 647);
    const auto *smg_648 = buffer.data(smg + 648);
    const auto *smg_650 = buffer.data(smg + 650);
    const auto *smg_655 = buffer.data(smg + 655);
    const auto *smg_656 = buffer.data(smg + 656);
    const auto *smg_657 = buffer.data(smg + 657);
    const auto *smg_658 = buffer.data(smg + 658);
    const auto *smg_659 = buffer.data(smg + 659);
    const auto *smg_660 = buffer.data(smg + 660);
    const auto *smg_662 = buffer.data(smg + 662);
    const auto *smg_663 = buffer.data(smg + 663);
    const auto *smg_665 = buffer.data(smg + 665);
    const auto *smg_670 = buffer.data(smg + 670);

#pragma omp simd aligned(t_815, t_816, t_817, t_818, pb_x, pc_x, pc_y, slh0_815, slh0_816, \
                         slh0_818, slg_464, slh1_815, slh1_816, slh1_818, \
                         smg_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_815[k] = pb_x[k] * slh0_815[k]
                   - f_8 * pc_x[k] * slh1_815[k];

        t_816[k] = pb_x[k] * slh0_816[k]
                   - f_8 * pc_x[k] * slh1_816[k];

        t_817[k] = f_14 * slg_464[k]
                   + f_3 * pc_y[k] * smg_584[k];

        t_818[k] = pb_x[k] * slh0_818[k]
                   - f_8 * pc_x[k] * slh1_818[k];
    }

#pragma omp simd aligned(t_819, t_820, t_821, pb_x, pc_x, pc_y, pc_z, slh0_819, slg_450, \
                         slg_465, slg_585, slh1_819, smg_585 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_819[k] = pb_x[k] * slh0_819[k]
                   + f_15 * slg_585[k]
                   - f_8 * pc_x[k] * slh1_819[k];

        t_820[k] = f_15 * slg_465[k]
                   + f_3 * pc_y[k] * smg_585[k];

        t_821[k] = f_11 * slg_450[k]
                   + f_3 * pc_z[k] * smg_585[k];
    }

#pragma omp simd aligned(t_822, t_823, t_824, pb_x, pc_x, pc_y, slh0_822, slh0_824, slg_467, \
                         slg_588, slg_590, slh1_822, slh1_824, \
                         smg_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_822[k] = pb_x[k] * slh0_822[k]
                   + f_11 * slg_588[k]
                   - f_8 * pc_x[k] * slh1_822[k];

        t_823[k] = f_15 * slg_467[k]
                   + f_3 * pc_y[k] * smg_587[k];

        t_824[k] = pb_x[k] * slh0_824[k]
                   + f_11 * slg_590[k]
                   - f_8 * pc_x[k] * slh1_824[k];
    }

#pragma omp simd aligned(t_825, t_826, t_827, pb_x, pc_x, pc_y, pc_z, slh0_825, slg_453, \
                         slg_470, slg_591, slh1_825, smg_588, smg_590 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_825[k] = pb_x[k] * slh0_825[k]
                   + f_10 * slg_591[k]
                   - f_8 * pc_x[k] * slh1_825[k];

        t_826[k] = f_11 * slg_453[k]
                   + f_3 * pc_z[k] * smg_588[k];

        t_827[k] = f_15 * slg_470[k]
                   + f_3 * pc_y[k] * smg_590[k];
    }

#pragma omp simd aligned(t_828, t_829, t_830, t_831, pb_x, pc_x, slh0_828, slg_594, slg_595, \
                         slg_596, slg_597, slh1_828, smg_595, smg_596, \
                         smg_597 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_828[k] = pb_x[k] * slh0_828[k]
                   + f_10 * slg_594[k]
                   - f_8 * pc_x[k] * slh1_828[k];

        t_829[k] = f_9 * slg_595[k]
                   + f_3 * pc_x[k] * smg_595[k];

        t_830[k] = f_9 * slg_596[k]
                   + f_3 * pc_x[k] * smg_596[k];

        t_831[k] = f_9 * slg_597[k]
                   + f_3 * pc_x[k] * smg_597[k];
    }

#pragma omp simd aligned(t_832, t_833, t_834, t_835, pb_x, pc_x, pc_z, slh0_834, slg_460, \
                         slg_598, slg_599, slh1_834, smg_595, smg_598, \
                         smg_599 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_832[k] = f_9 * slg_598[k]
                   + f_3 * pc_x[k] * smg_598[k];

        t_833[k] = f_9 * slg_599[k]
                   + f_3 * pc_x[k] * smg_599[k];

        t_834[k] = pb_x[k] * slh0_834[k]
                   - f_8 * pc_x[k] * slh1_834[k];

        t_835[k] = f_11 * slg_460[k]
                   + f_3 * pc_z[k] * smg_595[k];
    }

#pragma omp simd aligned(t_836, t_837, t_838, t_839, pb_x, pc_x, pc_y, slh0_836, slh0_837, \
                         slh0_839, slg_479, slh1_836, slh1_837, slh1_839, \
                         smg_599 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_836[k] = pb_x[k] * slh0_836[k]
                   - f_8 * pc_x[k] * slh1_836[k];

        t_837[k] = pb_x[k] * slh0_837[k]
                   - f_8 * pc_x[k] * slh1_837[k];

        t_838[k] = f_15 * slg_479[k]
                   + f_3 * pc_y[k] * smg_599[k];

        t_839[k] = pb_x[k] * slh0_839[k]
                   - f_8 * pc_x[k] * slh1_839[k];
    }

#pragma omp simd aligned(t_840, t_841, t_842, pb_x, pc_x, pc_y, pc_z, slh0_840, slg_465, \
                         slg_480, slg_600, slh1_840, smg_600 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_840[k] = pb_x[k] * slh0_840[k]
                   + f_15 * slg_600[k]
                   - f_8 * pc_x[k] * slh1_840[k];

        t_841[k] = f_16 * slg_480[k]
                   + f_3 * pc_y[k] * smg_600[k];

        t_842[k] = f_16 * slg_465[k]
                   + f_3 * pc_z[k] * smg_600[k];
    }

#pragma omp simd aligned(t_843, t_844, t_845, pb_x, pc_x, pc_y, slh0_843, slh0_845, slg_482, \
                         slg_603, slg_605, slh1_843, slh1_845, \
                         smg_602 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_843[k] = pb_x[k] * slh0_843[k]
                   + f_11 * slg_603[k]
                   - f_8 * pc_x[k] * slh1_843[k];

        t_844[k] = f_16 * slg_482[k]
                   + f_3 * pc_y[k] * smg_602[k];

        t_845[k] = pb_x[k] * slh0_845[k]
                   + f_11 * slg_605[k]
                   - f_8 * pc_x[k] * slh1_845[k];
    }

#pragma omp simd aligned(t_846, t_847, t_848, pb_x, pc_x, pc_y, pc_z, slh0_846, slg_468, \
                         slg_485, slg_606, slh1_846, smg_603, smg_605 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_846[k] = pb_x[k] * slh0_846[k]
                   + f_10 * slg_606[k]
                   - f_8 * pc_x[k] * slh1_846[k];

        t_847[k] = f_16 * slg_468[k]
                   + f_3 * pc_z[k] * smg_603[k];

        t_848[k] = f_16 * slg_485[k]
                   + f_3 * pc_y[k] * smg_605[k];
    }

#pragma omp simd aligned(t_849, t_850, t_851, t_852, pb_x, pc_x, slh0_849, slg_609, slg_610, \
                         slg_611, slg_612, slh1_849, smg_610, smg_611, \
                         smg_612 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_849[k] = pb_x[k] * slh0_849[k]
                   + f_10 * slg_609[k]
                   - f_8 * pc_x[k] * slh1_849[k];

        t_850[k] = f_9 * slg_610[k]
                   + f_3 * pc_x[k] * smg_610[k];

        t_851[k] = f_9 * slg_611[k]
                   + f_3 * pc_x[k] * smg_611[k];

        t_852[k] = f_9 * slg_612[k]
                   + f_3 * pc_x[k] * smg_612[k];
    }

#pragma omp simd aligned(t_853, t_854, t_855, t_856, pb_x, pc_x, pc_z, slh0_855, slg_475, \
                         slg_613, slg_614, slh1_855, smg_610, smg_613, \
                         smg_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_853[k] = f_9 * slg_613[k]
                   + f_3 * pc_x[k] * smg_613[k];

        t_854[k] = f_9 * slg_614[k]
                   + f_3 * pc_x[k] * smg_614[k];

        t_855[k] = pb_x[k] * slh0_855[k]
                   - f_8 * pc_x[k] * slh1_855[k];

        t_856[k] = f_16 * slg_475[k]
                   + f_3 * pc_z[k] * smg_610[k];
    }

#pragma omp simd aligned(t_857, t_858, t_859, t_860, pb_x, pc_x, pc_y, slh0_857, slh0_858, \
                         slh0_860, slg_494, slh1_857, slh1_858, slh1_860, \
                         smg_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_857[k] = pb_x[k] * slh0_857[k]
                   - f_8 * pc_x[k] * slh1_857[k];

        t_858[k] = pb_x[k] * slh0_858[k]
                   - f_8 * pc_x[k] * slh1_858[k];

        t_859[k] = f_16 * slg_494[k]
                   + f_3 * pc_y[k] * smg_614[k];

        t_860[k] = pb_x[k] * slh0_860[k]
                   - f_8 * pc_x[k] * slh1_860[k];
    }

#pragma omp simd aligned(t_861, t_862, t_863, pb_x, pc_x, pc_y, pc_z, slh0_861, slg_480, \
                         slg_495, slg_615, slh1_861, smg_615 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_861[k] = pb_x[k] * slh0_861[k]
                   + f_15 * slg_615[k]
                   - f_8 * pc_x[k] * slh1_861[k];

        t_862[k] = f_11 * slg_495[k]
                   + f_3 * pc_y[k] * smg_615[k];

        t_863[k] = f_15 * slg_480[k]
                   + f_3 * pc_z[k] * smg_615[k];
    }

#pragma omp simd aligned(t_864, t_865, t_866, pb_x, pc_x, pc_y, slh0_864, slh0_866, slg_497, \
                         slg_618, slg_620, slh1_864, slh1_866, \
                         smg_617 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_864[k] = pb_x[k] * slh0_864[k]
                   + f_11 * slg_618[k]
                   - f_8 * pc_x[k] * slh1_864[k];

        t_865[k] = f_11 * slg_497[k]
                   + f_3 * pc_y[k] * smg_617[k];

        t_866[k] = pb_x[k] * slh0_866[k]
                   + f_11 * slg_620[k]
                   - f_8 * pc_x[k] * slh1_866[k];
    }

#pragma omp simd aligned(t_867, t_868, t_869, pb_x, pc_x, pc_y, pc_z, slh0_867, slg_483, \
                         slg_500, slg_621, slh1_867, smg_618, smg_620 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_867[k] = pb_x[k] * slh0_867[k]
                   + f_10 * slg_621[k]
                   - f_8 * pc_x[k] * slh1_867[k];

        t_868[k] = f_15 * slg_483[k]
                   + f_3 * pc_z[k] * smg_618[k];

        t_869[k] = f_11 * slg_500[k]
                   + f_3 * pc_y[k] * smg_620[k];
    }

#pragma omp simd aligned(t_870, t_871, t_872, t_873, pb_x, pc_x, slh0_870, slg_624, slg_625, \
                         slg_626, slg_627, slh1_870, smg_625, smg_626, \
                         smg_627 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_870[k] = pb_x[k] * slh0_870[k]
                   + f_10 * slg_624[k]
                   - f_8 * pc_x[k] * slh1_870[k];

        t_871[k] = f_9 * slg_625[k]
                   + f_3 * pc_x[k] * smg_625[k];

        t_872[k] = f_9 * slg_626[k]
                   + f_3 * pc_x[k] * smg_626[k];

        t_873[k] = f_9 * slg_627[k]
                   + f_3 * pc_x[k] * smg_627[k];
    }

#pragma omp simd aligned(t_874, t_875, t_876, t_877, pb_x, pc_x, pc_z, slh0_876, slg_490, \
                         slg_628, slg_629, slh1_876, smg_625, smg_628, \
                         smg_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_874[k] = f_9 * slg_628[k]
                   + f_3 * pc_x[k] * smg_628[k];

        t_875[k] = f_9 * slg_629[k]
                   + f_3 * pc_x[k] * smg_629[k];

        t_876[k] = pb_x[k] * slh0_876[k]
                   - f_8 * pc_x[k] * slh1_876[k];

        t_877[k] = f_15 * slg_490[k]
                   + f_3 * pc_z[k] * smg_625[k];
    }

#pragma omp simd aligned(t_878, t_879, t_880, t_881, pb_x, pc_x, pc_y, slh0_878, slh0_879, \
                         slh0_881, slg_509, slh1_878, slh1_879, slh1_881, \
                         smg_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_878[k] = pb_x[k] * slh0_878[k]
                   - f_8 * pc_x[k] * slh1_878[k];

        t_879[k] = pb_x[k] * slh0_879[k]
                   - f_8 * pc_x[k] * slh1_879[k];

        t_880[k] = f_11 * slg_509[k]
                   + f_3 * pc_y[k] * smg_629[k];

        t_881[k] = pb_x[k] * slh0_881[k]
                   - f_8 * pc_x[k] * slh1_881[k];
    }

#pragma omp simd aligned(t_882, t_883, t_884, pb_x, pc_x, pc_y, pc_z, slh0_882, slg_495, \
                         slg_510, slg_630, slh1_882, smg_630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_882[k] = pb_x[k] * slh0_882[k]
                   + f_15 * slg_630[k]
                   - f_8 * pc_x[k] * slh1_882[k];

        t_883[k] = f_10 * slg_510[k]
                   + f_3 * pc_y[k] * smg_630[k];

        t_884[k] = f_14 * slg_495[k]
                   + f_3 * pc_z[k] * smg_630[k];
    }

#pragma omp simd aligned(t_885, t_886, t_887, pb_x, pc_x, pc_y, slh0_885, slh0_887, slg_512, \
                         slg_633, slg_635, slh1_885, slh1_887, \
                         smg_632 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_885[k] = pb_x[k] * slh0_885[k]
                   + f_11 * slg_633[k]
                   - f_8 * pc_x[k] * slh1_885[k];

        t_886[k] = f_10 * slg_512[k]
                   + f_3 * pc_y[k] * smg_632[k];

        t_887[k] = pb_x[k] * slh0_887[k]
                   + f_11 * slg_635[k]
                   - f_8 * pc_x[k] * slh1_887[k];
    }

#pragma omp simd aligned(t_888, t_889, t_890, pb_x, pc_x, pc_y, pc_z, slh0_888, slg_498, \
                         slg_515, slg_636, slh1_888, smg_633, smg_635 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_888[k] = pb_x[k] * slh0_888[k]
                   + f_10 * slg_636[k]
                   - f_8 * pc_x[k] * slh1_888[k];

        t_889[k] = f_14 * slg_498[k]
                   + f_3 * pc_z[k] * smg_633[k];

        t_890[k] = f_10 * slg_515[k]
                   + f_3 * pc_y[k] * smg_635[k];
    }

#pragma omp simd aligned(t_891, t_892, t_893, t_894, pb_x, pc_x, slh0_891, slg_639, slg_640, \
                         slg_641, slg_642, slh1_891, smg_640, smg_641, \
                         smg_642 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_891[k] = pb_x[k] * slh0_891[k]
                   + f_10 * slg_639[k]
                   - f_8 * pc_x[k] * slh1_891[k];

        t_892[k] = f_9 * slg_640[k]
                   + f_3 * pc_x[k] * smg_640[k];

        t_893[k] = f_9 * slg_641[k]
                   + f_3 * pc_x[k] * smg_641[k];

        t_894[k] = f_9 * slg_642[k]
                   + f_3 * pc_x[k] * smg_642[k];
    }

#pragma omp simd aligned(t_895, t_896, t_897, t_898, pb_x, pc_x, pc_z, slh0_897, slg_505, \
                         slg_643, slg_644, slh1_897, smg_640, smg_643, \
                         smg_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_895[k] = f_9 * slg_643[k]
                   + f_3 * pc_x[k] * smg_643[k];

        t_896[k] = f_9 * slg_644[k]
                   + f_3 * pc_x[k] * smg_644[k];

        t_897[k] = pb_x[k] * slh0_897[k]
                   - f_8 * pc_x[k] * slh1_897[k];

        t_898[k] = f_14 * slg_505[k]
                   + f_3 * pc_z[k] * smg_640[k];
    }

#pragma omp simd aligned(t_899, t_900, t_901, t_902, pb_x, pc_x, pc_y, slh0_899, slh0_900, \
                         slh0_902, slg_524, slh1_899, slh1_900, slh1_902, \
                         smg_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_899[k] = pb_x[k] * slh0_899[k]
                   - f_8 * pc_x[k] * slh1_899[k];

        t_900[k] = pb_x[k] * slh0_900[k]
                   - f_8 * pc_x[k] * slh1_900[k];

        t_901[k] = f_10 * slg_524[k]
                   + f_3 * pc_y[k] * smg_644[k];

        t_902[k] = pb_x[k] * slh0_902[k]
                   - f_8 * pc_x[k] * slh1_902[k];
    }

#pragma omp simd aligned(t_903, t_904, t_905, pb_y, pc_y, pc_z, slh0_735, slg_510, slg_525, \
                         slh1_735, smg_645 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_903[k] = pb_y[k] * slh0_735[k]
                   - f_8 * pc_y[k] * slh1_735[k];

        t_904[k] = f_9 * slg_525[k]
                   + f_3 * pc_y[k] * smg_645[k];

        t_905[k] = f_13 * slg_510[k]
                   + f_3 * pc_z[k] * smg_645[k];
    }

#pragma omp simd aligned(t_906, t_907, t_908, pb_x, pb_y, pc_x, pc_y, slh0_740, slh0_906, \
                         slg_527, slg_648, slh1_740, slh1_906, \
                         smg_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_906[k] = pb_x[k] * slh0_906[k]
                   + f_11 * slg_648[k]
                   - f_8 * pc_x[k] * slh1_906[k];

        t_907[k] = f_9 * slg_527[k]
                   + f_3 * pc_y[k] * smg_647[k];

        t_908[k] = pb_y[k] * slh0_740[k]
                   - f_8 * pc_y[k] * slh1_740[k];
    }

#pragma omp simd aligned(t_909, t_910, t_911, pb_x, pc_x, pc_y, pc_z, slh0_909, slg_513, \
                         slg_530, slg_651, slh1_909, smg_648, smg_650 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_909[k] = pb_x[k] * slh0_909[k]
                   + f_10 * slg_651[k]
                   - f_8 * pc_x[k] * slh1_909[k];

        t_910[k] = f_13 * slg_513[k]
                   + f_3 * pc_z[k] * smg_648[k];

        t_911[k] = f_9 * slg_530[k]
                   + f_3 * pc_y[k] * smg_650[k];
    }

#pragma omp simd aligned(t_912, t_913, t_914, t_915, pb_y, pc_x, pc_y, slh0_744, slg_655, \
                         slg_656, slg_657, slh1_744, smg_655, smg_656, \
                         smg_657 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_912[k] = pb_y[k] * slh0_744[k]
                   - f_8 * pc_y[k] * slh1_744[k];

        t_913[k] = f_9 * slg_655[k]
                   + f_3 * pc_x[k] * smg_655[k];

        t_914[k] = f_9 * slg_656[k]
                   + f_3 * pc_x[k] * smg_656[k];

        t_915[k] = f_9 * slg_657[k]
                   + f_3 * pc_x[k] * smg_657[k];
    }

#pragma omp simd aligned(t_916, t_917, t_918, t_919, pb_x, pc_x, pc_z, slh0_918, slg_520, \
                         slg_658, slg_659, slh1_918, smg_655, smg_658, \
                         smg_659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_916[k] = f_9 * slg_658[k]
                   + f_3 * pc_x[k] * smg_658[k];

        t_917[k] = f_9 * slg_659[k]
                   + f_3 * pc_x[k] * smg_659[k];

        t_918[k] = pb_x[k] * slh0_918[k]
                   - f_8 * pc_x[k] * slh1_918[k];

        t_919[k] = f_13 * slg_520[k]
                   + f_3 * pc_z[k] * smg_655[k];
    }

#pragma omp simd aligned(t_920, t_921, t_922, t_923, pb_x, pc_x, pc_y, slh0_920, slh0_921, \
                         slh0_923, slg_539, slh1_920, slh1_921, slh1_923, \
                         smg_659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_920[k] = pb_x[k] * slh0_920[k]
                   - f_8 * pc_x[k] * slh1_920[k];

        t_921[k] = pb_x[k] * slh0_921[k]
                   - f_8 * pc_x[k] * slh1_921[k];

        t_922[k] = f_9 * slg_539[k]
                   + f_3 * pc_y[k] * smg_659[k];

        t_923[k] = pb_x[k] * slh0_923[k]
                   - f_8 * pc_x[k] * slh1_923[k];
    }

#pragma omp simd aligned(t_924, t_925, t_926, t_927, pb_x, pc_x, pc_y, pc_z, slh0_924, \
                         slh0_927, slg_525, slg_660, slg_663, slh1_924, slh1_927, \
                         smg_660 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_924[k] = pb_x[k] * slh0_924[k]
                   + f_15 * slg_660[k]
                   - f_8 * pc_x[k] * slh1_924[k];

        t_925[k] = f_3 * pc_y[k] * smg_660[k];

        t_926[k] = f_12 * slg_525[k]
                   + f_3 * pc_z[k] * smg_660[k];

        t_927[k] = pb_x[k] * slh0_927[k]
                   + f_11 * slg_663[k]
                   - f_8 * pc_x[k] * slh1_927[k];
    }

#pragma omp simd aligned(t_928, t_929, t_930, pb_x, pc_x, pc_y, slh0_929, slh0_930, slg_665, \
                         slg_666, slh1_929, slh1_930, smg_662 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_928[k] = f_3 * pc_y[k] * smg_662[k];

        t_929[k] = pb_x[k] * slh0_929[k]
                   + f_11 * slg_665[k]
                   - f_8 * pc_x[k] * slh1_929[k];

        t_930[k] = pb_x[k] * slh0_930[k]
                   + f_10 * slg_666[k]
                   - f_8 * pc_x[k] * slh1_930[k];
    }

#pragma omp simd aligned(t_931, t_932, t_933, t_934, pb_x, pc_x, pc_y, pc_z, slh0_933, \
                         slg_528, slg_669, slg_670, slh1_933, smg_663, smg_665, \
                         smg_670 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_931[k] = f_12 * slg_528[k]
                   + f_3 * pc_z[k] * smg_663[k];

        t_932[k] = f_3 * pc_y[k] * smg_665[k];

        t_933[k] = pb_x[k] * slh0_933[k]
                   + f_10 * slg_669[k]
                   - f_8 * pc_x[k] * slh1_933[k];

        t_934[k] = f_9 * slg_670[k]
                   + f_3 * pc_x[k] * smg_670[k];
    }
}

static auto
compute_prim_smh_three_center_electron_repulsion_0_piece8(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t slh0,
                                                          const size_t slg, const size_t slh1,
                                                          const size_t smf0, const size_t smf1,
                                                          const size_t smg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / q;
    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.0 / gamma;
    const auto f_5 = p / (gamma * q);
    const auto f_6 = 0.5 / gamma;
    const auto f_7 = 0.5 * p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_12 = 4.0 / q;
    const auto f_13 = 3.5 / q;
    const auto f_14 = 3.0 / q;
    const auto f_15 = 2.5 / q;
    const auto f_16 = 2.0 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *slh0_756 = buffer.data(slh0 + 756);
    const auto *slh0_759 = buffer.data(slh0 + 759);
    const auto *slh0_762 = buffer.data(slh0 + 762);
    const auto *slh0_771 = buffer.data(slh0 + 771);
    const auto *slh0_773 = buffer.data(slh0 + 773);
    const auto *slh0_774 = buffer.data(slh0 + 774);
    const auto *slh0_939 = buffer.data(slh0 + 939);
    const auto *slh0_941 = buffer.data(slh0 + 941);
    const auto *slh0_942 = buffer.data(slh0 + 942);
    const auto *slh0_944 = buffer.data(slh0 + 944);

    const auto *slg_535 = buffer.data(slg + 535);
    const auto *slg_540 = buffer.data(slg + 540);
    const auto *slg_542 = buffer.data(slg + 542);
    const auto *slg_543 = buffer.data(slg + 543);
    const auto *slg_545 = buffer.data(slg + 545);
    const auto *slg_550 = buffer.data(slg + 550);
    const auto *slg_551 = buffer.data(slg + 551);
    const auto *slg_552 = buffer.data(slg + 552);
    const auto *slg_553 = buffer.data(slg + 553);
    const auto *slg_554 = buffer.data(slg + 554);
    const auto *slg_555 = buffer.data(slg + 555);
    const auto *slg_557 = buffer.data(slg + 557);
    const auto *slg_558 = buffer.data(slg + 558);
    const auto *slg_560 = buffer.data(slg + 560);
    const auto *slg_565 = buffer.data(slg + 565);
    const auto *slg_569 = buffer.data(slg + 569);
    const auto *slg_570 = buffer.data(slg + 570);
    const auto *slg_572 = buffer.data(slg + 572);
    const auto *slg_573 = buffer.data(slg + 573);
    const auto *slg_575 = buffer.data(slg + 575);
    const auto *slg_580 = buffer.data(slg + 580);
    const auto *slg_582 = buffer.data(slg + 582);
    const auto *slg_583 = buffer.data(slg + 583);
    const auto *slg_584 = buffer.data(slg + 584);
    const auto *slg_585 = buffer.data(slg + 585);
    const auto *slg_587 = buffer.data(slg + 587);
    const auto *slg_588 = buffer.data(slg + 588);
    const auto *slg_590 = buffer.data(slg + 590);
    const auto *slg_595 = buffer.data(slg + 595);
    const auto *slg_597 = buffer.data(slg + 597);
    const auto *slg_598 = buffer.data(slg + 598);
    const auto *slg_599 = buffer.data(slg + 599);
    const auto *slg_600 = buffer.data(slg + 600);
    const auto *slg_602 = buffer.data(slg + 602);
    const auto *slg_603 = buffer.data(slg + 603);
    const auto *slg_605 = buffer.data(slg + 605);
    const auto *slg_610 = buffer.data(slg + 610);
    const auto *slg_612 = buffer.data(slg + 612);
    const auto *slg_613 = buffer.data(slg + 613);
    const auto *slg_614 = buffer.data(slg + 614);
    const auto *slg_615 = buffer.data(slg + 615);
    const auto *slg_617 = buffer.data(slg + 617);
    const auto *slg_620 = buffer.data(slg + 620);
    const auto *slg_671 = buffer.data(slg + 671);
    const auto *slg_672 = buffer.data(slg + 672);
    const auto *slg_673 = buffer.data(slg + 673);
    const auto *slg_674 = buffer.data(slg + 674);

    const auto *slh1_756 = buffer.data(slh1 + 756);
    const auto *slh1_759 = buffer.data(slh1 + 759);
    const auto *slh1_762 = buffer.data(slh1 + 762);
    const auto *slh1_771 = buffer.data(slh1 + 771);
    const auto *slh1_773 = buffer.data(slh1 + 773);
    const auto *slh1_774 = buffer.data(slh1 + 774);
    const auto *slh1_939 = buffer.data(slh1 + 939);
    const auto *slh1_941 = buffer.data(slh1 + 941);
    const auto *slh1_942 = buffer.data(slh1 + 942);
    const auto *slh1_944 = buffer.data(slh1 + 944);

    const auto *smf0_450 = buffer.data(smf0 + 450);
    const auto *smf0_453 = buffer.data(smf0 + 453);
    const auto *smf0_455 = buffer.data(smf0 + 455);
    const auto *smf0_456 = buffer.data(smf0 + 456);
    const auto *smf0_458 = buffer.data(smf0 + 458);
    const auto *smf0_459 = buffer.data(smf0 + 459);
    const auto *smf0_465 = buffer.data(smf0 + 465);
    const auto *smf0_469 = buffer.data(smf0 + 469);
    const auto *smf0_470 = buffer.data(smf0 + 470);
    const auto *smf0_473 = buffer.data(smf0 + 473);
    const auto *smf0_475 = buffer.data(smf0 + 475);
    const auto *smf0_476 = buffer.data(smf0 + 476);
    const auto *smf0_478 = buffer.data(smf0 + 478);
    const auto *smf0_479 = buffer.data(smf0 + 479);
    const auto *smf0_480 = buffer.data(smf0 + 480);
    const auto *smf0_483 = buffer.data(smf0 + 483);
    const auto *smf0_485 = buffer.data(smf0 + 485);
    const auto *smf0_486 = buffer.data(smf0 + 486);
    const auto *smf0_488 = buffer.data(smf0 + 488);
    const auto *smf0_489 = buffer.data(smf0 + 489);
    const auto *smf0_490 = buffer.data(smf0 + 490);
    const auto *smf0_493 = buffer.data(smf0 + 493);
    const auto *smf0_495 = buffer.data(smf0 + 495);
    const auto *smf0_496 = buffer.data(smf0 + 496);
    const auto *smf0_498 = buffer.data(smf0 + 498);
    const auto *smf0_499 = buffer.data(smf0 + 499);
    const auto *smf0_500 = buffer.data(smf0 + 500);
    const auto *smf0_503 = buffer.data(smf0 + 503);
    const auto *smf0_505 = buffer.data(smf0 + 505);
    const auto *smf0_506 = buffer.data(smf0 + 506);
    const auto *smf0_509 = buffer.data(smf0 + 509);

    const auto *smf1_450 = buffer.data(smf1 + 450);
    const auto *smf1_453 = buffer.data(smf1 + 453);
    const auto *smf1_455 = buffer.data(smf1 + 455);
    const auto *smf1_456 = buffer.data(smf1 + 456);
    const auto *smf1_458 = buffer.data(smf1 + 458);
    const auto *smf1_459 = buffer.data(smf1 + 459);
    const auto *smf1_465 = buffer.data(smf1 + 465);
    const auto *smf1_469 = buffer.data(smf1 + 469);
    const auto *smf1_470 = buffer.data(smf1 + 470);
    const auto *smf1_473 = buffer.data(smf1 + 473);
    const auto *smf1_475 = buffer.data(smf1 + 475);
    const auto *smf1_476 = buffer.data(smf1 + 476);
    const auto *smf1_478 = buffer.data(smf1 + 478);
    const auto *smf1_479 = buffer.data(smf1 + 479);
    const auto *smf1_480 = buffer.data(smf1 + 480);
    const auto *smf1_483 = buffer.data(smf1 + 483);
    const auto *smf1_485 = buffer.data(smf1 + 485);
    const auto *smf1_486 = buffer.data(smf1 + 486);
    const auto *smf1_488 = buffer.data(smf1 + 488);
    const auto *smf1_489 = buffer.data(smf1 + 489);
    const auto *smf1_490 = buffer.data(smf1 + 490);
    const auto *smf1_493 = buffer.data(smf1 + 493);
    const auto *smf1_495 = buffer.data(smf1 + 495);
    const auto *smf1_496 = buffer.data(smf1 + 496);
    const auto *smf1_498 = buffer.data(smf1 + 498);
    const auto *smf1_499 = buffer.data(smf1 + 499);
    const auto *smf1_500 = buffer.data(smf1 + 500);
    const auto *smf1_503 = buffer.data(smf1 + 503);
    const auto *smf1_505 = buffer.data(smf1 + 505);
    const auto *smf1_506 = buffer.data(smf1 + 506);
    const auto *smf1_509 = buffer.data(smf1 + 509);

    const auto *smg_670 = buffer.data(smg + 670);
    const auto *smg_671 = buffer.data(smg + 671);
    const auto *smg_672 = buffer.data(smg + 672);
    const auto *smg_673 = buffer.data(smg + 673);
    const auto *smg_674 = buffer.data(smg + 674);
    const auto *smg_675 = buffer.data(smg + 675);
    const auto *smg_677 = buffer.data(smg + 677);
    const auto *smg_678 = buffer.data(smg + 678);
    const auto *smg_680 = buffer.data(smg + 680);
    const auto *smg_681 = buffer.data(smg + 681);
    const auto *smg_684 = buffer.data(smg + 684);
    const auto *smg_685 = buffer.data(smg + 685);
    const auto *smg_686 = buffer.data(smg + 686);
    const auto *smg_687 = buffer.data(smg + 687);
    const auto *smg_688 = buffer.data(smg + 688);
    const auto *smg_689 = buffer.data(smg + 689);
    const auto *smg_690 = buffer.data(smg + 690);
    const auto *smg_692 = buffer.data(smg + 692);
    const auto *smg_693 = buffer.data(smg + 693);
    const auto *smg_695 = buffer.data(smg + 695);
    const auto *smg_699 = buffer.data(smg + 699);
    const auto *smg_700 = buffer.data(smg + 700);
    const auto *smg_701 = buffer.data(smg + 701);
    const auto *smg_702 = buffer.data(smg + 702);
    const auto *smg_703 = buffer.data(smg + 703);
    const auto *smg_704 = buffer.data(smg + 704);
    const auto *smg_705 = buffer.data(smg + 705);
    const auto *smg_707 = buffer.data(smg + 707);
    const auto *smg_708 = buffer.data(smg + 708);
    const auto *smg_710 = buffer.data(smg + 710);
    const auto *smg_711 = buffer.data(smg + 711);
    const auto *smg_714 = buffer.data(smg + 714);
    const auto *smg_715 = buffer.data(smg + 715);
    const auto *smg_716 = buffer.data(smg + 716);
    const auto *smg_717 = buffer.data(smg + 717);
    const auto *smg_718 = buffer.data(smg + 718);
    const auto *smg_719 = buffer.data(smg + 719);
    const auto *smg_720 = buffer.data(smg + 720);
    const auto *smg_722 = buffer.data(smg + 722);
    const auto *smg_723 = buffer.data(smg + 723);
    const auto *smg_725 = buffer.data(smg + 725);
    const auto *smg_726 = buffer.data(smg + 726);
    const auto *smg_729 = buffer.data(smg + 729);
    const auto *smg_730 = buffer.data(smg + 730);
    const auto *smg_731 = buffer.data(smg + 731);
    const auto *smg_732 = buffer.data(smg + 732);
    const auto *smg_733 = buffer.data(smg + 733);
    const auto *smg_734 = buffer.data(smg + 734);
    const auto *smg_735 = buffer.data(smg + 735);
    const auto *smg_737 = buffer.data(smg + 737);
    const auto *smg_738 = buffer.data(smg + 738);
    const auto *smg_740 = buffer.data(smg + 740);
    const auto *smg_741 = buffer.data(smg + 741);
    const auto *smg_744 = buffer.data(smg + 744);
    const auto *smg_745 = buffer.data(smg + 745);
    const auto *smg_746 = buffer.data(smg + 746);
    const auto *smg_747 = buffer.data(smg + 747);
    const auto *smg_748 = buffer.data(smg + 748);
    const auto *smg_749 = buffer.data(smg + 749);
    const auto *smg_750 = buffer.data(smg + 750);
    const auto *smg_752 = buffer.data(smg + 752);
    const auto *smg_753 = buffer.data(smg + 753);
    const auto *smg_755 = buffer.data(smg + 755);
    const auto *smg_756 = buffer.data(smg + 756);
    const auto *smg_759 = buffer.data(smg + 759);
    const auto *smg_760 = buffer.data(smg + 760);

#pragma omp simd aligned(t_935, t_936, t_937, t_938, pc_x, slg_671, slg_672, slg_673, slg_674, \
                         smg_671, smg_672, smg_673, smg_674 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_935[k] = f_9 * slg_671[k]
                   + f_3 * pc_x[k] * smg_671[k];

        t_936[k] = f_9 * slg_672[k]
                   + f_3 * pc_x[k] * smg_672[k];

        t_937[k] = f_9 * slg_673[k]
                   + f_3 * pc_x[k] * smg_673[k];

        t_938[k] = f_9 * slg_674[k]
                   + f_3 * pc_x[k] * smg_674[k];
    }

#pragma omp simd aligned(t_939, t_940, t_941, t_942, pb_x, pc_x, pc_z, slh0_939, slh0_941, \
                         slh0_942, slg_535, slh1_939, slh1_941, slh1_942, \
                         smg_670 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_939[k] = pb_x[k] * slh0_939[k]
                   - f_8 * pc_x[k] * slh1_939[k];

        t_940[k] = f_12 * slg_535[k]
                   + f_3 * pc_z[k] * smg_670[k];

        t_941[k] = pb_x[k] * slh0_941[k]
                   - f_8 * pc_x[k] * slh1_941[k];

        t_942[k] = pb_x[k] * slh0_942[k]
                   - f_8 * pc_x[k] * slh1_942[k];
    }

#pragma omp simd aligned(t_943, t_944, t_945, t_946, t_947, pb_x, pc_x, pc_y, pc_z, slh0_944, \
                         slg_540, slh1_944, smf0_450, smf1_450, smg_674, \
                         smg_675 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_943[k] = f_3 * pc_y[k] * smg_674[k];

        t_944[k] = pb_x[k] * slh0_944[k]
                   - f_8 * pc_x[k] * slh1_944[k];

        t_945[k] = f_1 * smf0_450[k]
                   - f_2 * smf1_450[k]
                   + f_3 * pc_x[k] * smg_675[k];

        t_946[k] = f_0 * slg_540[k]
                   + f_3 * pc_y[k] * smg_675[k];

        t_947[k] = f_3 * pc_z[k] * smg_675[k];
    }

#pragma omp simd aligned(t_948, t_949, t_950, pc_x, pc_y, slg_542, smf0_453, smf0_455, \
                         smf1_453, smf1_455, smg_677, smg_678, \
                         smg_680 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_948[k] = f_4 * smf0_453[k]
                   - f_5 * smf1_453[k]
                   + f_3 * pc_x[k] * smg_678[k];

        t_949[k] = f_0 * slg_542[k]
                   + f_3 * pc_y[k] * smg_677[k];

        t_950[k] = f_4 * smf0_455[k]
                   - f_5 * smf1_455[k]
                   + f_3 * pc_x[k] * smg_680[k];
    }

#pragma omp simd aligned(t_951, t_952, t_953, t_954, pc_x, pc_y, pc_z, slg_545, smf0_456, \
                         smf0_459, smf1_456, smf1_459, smg_678, smg_680, smg_681, \
                         smg_684 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_951[k] = f_6 * smf0_456[k]
                   - f_7 * smf1_456[k]
                   + f_3 * pc_x[k] * smg_681[k];

        t_952[k] = f_3 * pc_z[k] * smg_678[k];

        t_953[k] = f_0 * slg_545[k]
                   + f_3 * pc_y[k] * smg_680[k];

        t_954[k] = f_6 * smf0_459[k]
                   - f_7 * smf1_459[k]
                   + f_3 * pc_x[k] * smg_684[k];
    }

#pragma omp simd aligned(t_955, t_956, t_957, t_958, t_959, t_960, pc_x, pc_y, slg_550, \
                         smf0_456, smf1_456, smg_685, smg_686, smg_687, smg_688, \
                         smg_689 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_955[k] = f_3 * pc_x[k] * smg_685[k];

        t_956[k] = f_3 * pc_x[k] * smg_686[k];

        t_957[k] = f_3 * pc_x[k] * smg_687[k];

        t_958[k] = f_3 * pc_x[k] * smg_688[k];

        t_959[k] = f_3 * pc_x[k] * smg_689[k];

        t_960[k] = f_0 * slg_550[k]
                   + f_1 * smf0_456[k]
                   - f_2 * smf1_456[k]
                   + f_3 * pc_y[k] * smg_685[k];
    }

#pragma omp simd aligned(t_961, t_962, t_963, pc_y, pc_z, slg_552, slg_553, smf0_458, \
                         smf0_459, smf1_458, smf1_459, smg_685, smg_687, \
                         smg_688 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_961[k] = f_3 * pc_z[k] * smg_685[k];

        t_962[k] = f_0 * slg_552[k]
                   + f_4 * smf0_458[k]
                   - f_5 * smf1_458[k]
                   + f_3 * pc_y[k] * smg_687[k];

        t_963[k] = f_0 * slg_553[k]
                   + f_6 * smf0_459[k]
                   - f_7 * smf1_459[k]
                   + f_3 * pc_y[k] * smg_688[k];
    }

#pragma omp simd aligned(t_964, t_965, t_966, t_967, pb_z, pc_y, pc_z, slh0_756, slg_554, \
                         slg_555, slh1_756, smf0_459, smf1_459, smg_689, \
                         smg_690 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_964[k] = f_0 * slg_554[k]
                   + f_3 * pc_y[k] * smg_689[k];

        t_965[k] = f_1 * smf0_459[k]
                   - f_2 * smf1_459[k]
                   + f_3 * pc_z[k] * smg_689[k];

        t_966[k] = pb_z[k] * slh0_756[k]
                   - f_8 * pc_z[k] * slh1_756[k];

        t_967[k] = f_12 * slg_555[k]
                   + f_3 * pc_y[k] * smg_690[k];
    }

#pragma omp simd aligned(t_968, t_969, t_970, pb_z, pc_y, pc_z, slh0_759, slg_540, slg_557, \
                         slh1_759, smg_690, smg_692 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_968[k] = f_9 * slg_540[k]
                   + f_3 * pc_z[k] * smg_690[k];

        t_969[k] = pb_z[k] * slh0_759[k]
                   - f_8 * pc_z[k] * slh1_759[k];

        t_970[k] = f_12 * slg_557[k]
                   + f_3 * pc_y[k] * smg_692[k];
    }

#pragma omp simd aligned(t_971, t_972, t_973, t_974, pb_z, pc_x, pc_y, pc_z, slh0_762, \
                         slg_543, slg_560, slh1_762, smf0_465, smf1_465, smg_693, \
                         smg_695 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_971[k] = f_4 * smf0_465[k]
                   - f_5 * smf1_465[k]
                   + f_3 * pc_x[k] * smg_695[k];

        t_972[k] = pb_z[k] * slh0_762[k]
                   - f_8 * pc_z[k] * slh1_762[k];

        t_973[k] = f_9 * slg_543[k]
                   + f_3 * pc_z[k] * smg_693[k];

        t_974[k] = f_12 * slg_560[k]
                   + f_3 * pc_y[k] * smg_695[k];
    }

#pragma omp simd aligned(t_975, t_976, t_977, t_978, t_979, t_980, pc_x, smf0_469, smf1_469, \
                         smg_699, smg_700, smg_701, smg_702, smg_703, \
                         smg_704 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_975[k] = f_6 * smf0_469[k]
                   - f_7 * smf1_469[k]
                   + f_3 * pc_x[k] * smg_699[k];

        t_976[k] = f_3 * pc_x[k] * smg_700[k];

        t_977[k] = f_3 * pc_x[k] * smg_701[k];

        t_978[k] = f_3 * pc_x[k] * smg_702[k];

        t_979[k] = f_3 * pc_x[k] * smg_703[k];

        t_980[k] = f_3 * pc_x[k] * smg_704[k];
    }

#pragma omp simd aligned(t_981, t_982, t_983, t_984, pb_z, pc_z, slh0_771, slh0_773, slh0_774, \
                         slg_550, slg_551, slg_552, slh1_771, slh1_773, slh1_774, \
                         smg_700 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_981[k] = pb_z[k] * slh0_771[k]
                   - f_8 * pc_z[k] * slh1_771[k];

        t_982[k] = f_9 * slg_550[k]
                   + f_3 * pc_z[k] * smg_700[k];

        t_983[k] = pb_z[k] * slh0_773[k]
                   + f_10 * slg_551[k]
                   - f_8 * pc_z[k] * slh1_773[k];

        t_984[k] = pb_z[k] * slh0_774[k]
                   + f_11 * slg_552[k]
                   - f_8 * pc_z[k] * slh1_774[k];
    }

#pragma omp simd aligned(t_985, t_986, t_987, t_988, pc_x, pc_y, pc_z, slg_554, slg_569, \
                         slg_570, smf0_469, smf0_470, smf1_469, smf1_470, smg_704, \
                         smg_705 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_985[k] = f_12 * slg_569[k]
                   + f_3 * pc_y[k] * smg_704[k];

        t_986[k] = f_9 * slg_554[k]
                   + f_1 * smf0_469[k]
                   - f_2 * smf1_469[k]
                   + f_3 * pc_z[k] * smg_704[k];

        t_987[k] = f_1 * smf0_470[k]
                   - f_2 * smf1_470[k]
                   + f_3 * pc_x[k] * smg_705[k];

        t_988[k] = f_13 * slg_570[k]
                   + f_3 * pc_y[k] * smg_705[k];
    }

#pragma omp simd aligned(t_989, t_990, t_991, pc_x, pc_y, pc_z, slg_555, slg_572, smf0_473, \
                         smf1_473, smg_705, smg_707, smg_708 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_989[k] = f_10 * slg_555[k]
                   + f_3 * pc_z[k] * smg_705[k];

        t_990[k] = f_4 * smf0_473[k]
                   - f_5 * smf1_473[k]
                   + f_3 * pc_x[k] * smg_708[k];

        t_991[k] = f_13 * slg_572[k]
                   + f_3 * pc_y[k] * smg_707[k];
    }

#pragma omp simd aligned(t_992, t_993, t_994, t_995, pc_x, pc_y, pc_z, slg_558, slg_575, \
                         smf0_475, smf0_476, smf1_475, smf1_476, smg_708, smg_710, \
                         smg_711 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_992[k] = f_4 * smf0_475[k]
                   - f_5 * smf1_475[k]
                   + f_3 * pc_x[k] * smg_710[k];

        t_993[k] = f_6 * smf0_476[k]
                   - f_7 * smf1_476[k]
                   + f_3 * pc_x[k] * smg_711[k];

        t_994[k] = f_10 * slg_558[k]
                   + f_3 * pc_z[k] * smg_708[k];

        t_995[k] = f_13 * slg_575[k]
                   + f_3 * pc_y[k] * smg_710[k];
    }

#pragma omp simd aligned(t_996, t_997, t_998, t_999, t_1000, t_1001, pc_x, smf0_479, smf1_479, \
                         smg_714, smg_715, smg_716, smg_717, smg_718, \
                         smg_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_996[k] = f_6 * smf0_479[k]
                   - f_7 * smf1_479[k]
                   + f_3 * pc_x[k] * smg_714[k];

        t_997[k] = f_3 * pc_x[k] * smg_715[k];

        t_998[k] = f_3 * pc_x[k] * smg_716[k];

        t_999[k] = f_3 * pc_x[k] * smg_717[k];

        t_1000[k] = f_3 * pc_x[k] * smg_718[k];

        t_1001[k] = f_3 * pc_x[k] * smg_719[k];
    }

#pragma omp simd aligned(t_1002, t_1003, t_1004, pc_y, pc_z, slg_565, slg_580, slg_582, \
                         smf0_476, smf0_478, smf1_476, smf1_478, smg_715, \
                         smg_717 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1002[k] = f_13 * slg_580[k]
                    + f_1 * smf0_476[k]
                    - f_2 * smf1_476[k]
                    + f_3 * pc_y[k] * smg_715[k];

        t_1003[k] = f_10 * slg_565[k]
                    + f_3 * pc_z[k] * smg_715[k];

        t_1004[k] = f_13 * slg_582[k]
                    + f_4 * smf0_478[k]
                    - f_5 * smf1_478[k]
                    + f_3 * pc_y[k] * smg_717[k];
    }

#pragma omp simd aligned(t_1005, t_1006, t_1007, pc_y, pc_z, slg_569, slg_583, slg_584, \
                         smf0_479, smf1_479, smg_718, smg_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1005[k] = f_13 * slg_583[k]
                    + f_6 * smf0_479[k]
                    - f_7 * smf1_479[k]
                    + f_3 * pc_y[k] * smg_718[k];

        t_1006[k] = f_13 * slg_584[k]
                    + f_3 * pc_y[k] * smg_719[k];

        t_1007[k] = f_10 * slg_569[k]
                    + f_1 * smf0_479[k]
                    - f_2 * smf1_479[k]
                    + f_3 * pc_z[k] * smg_719[k];
    }

#pragma omp simd aligned(t_1008, t_1009, t_1010, t_1011, pc_x, pc_y, pc_z, slg_570, slg_585, \
                         smf0_480, smf0_483, smf1_480, smf1_483, smg_720, \
                         smg_723 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1008[k] = f_1 * smf0_480[k]
                    - f_2 * smf1_480[k]
                    + f_3 * pc_x[k] * smg_720[k];

        t_1009[k] = f_14 * slg_585[k]
                    + f_3 * pc_y[k] * smg_720[k];

        t_1010[k] = f_11 * slg_570[k]
                    + f_3 * pc_z[k] * smg_720[k];

        t_1011[k] = f_4 * smf0_483[k]
                    - f_5 * smf1_483[k]
                    + f_3 * pc_x[k] * smg_723[k];
    }

#pragma omp simd aligned(t_1012, t_1013, t_1014, pc_x, pc_y, slg_587, smf0_485, smf0_486, \
                         smf1_485, smf1_486, smg_722, smg_725, \
                         smg_726 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1012[k] = f_14 * slg_587[k]
                    + f_3 * pc_y[k] * smg_722[k];

        t_1013[k] = f_4 * smf0_485[k]
                    - f_5 * smf1_485[k]
                    + f_3 * pc_x[k] * smg_725[k];

        t_1014[k] = f_6 * smf0_486[k]
                    - f_7 * smf1_486[k]
                    + f_3 * pc_x[k] * smg_726[k];
    }

#pragma omp simd aligned(t_1015, t_1016, t_1017, t_1018, pc_x, pc_y, pc_z, slg_573, slg_590, \
                         smf0_489, smf1_489, smg_723, smg_725, smg_729, \
                         smg_730 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1015[k] = f_11 * slg_573[k]
                    + f_3 * pc_z[k] * smg_723[k];

        t_1016[k] = f_14 * slg_590[k]
                    + f_3 * pc_y[k] * smg_725[k];

        t_1017[k] = f_6 * smf0_489[k]
                    - f_7 * smf1_489[k]
                    + f_3 * pc_x[k] * smg_729[k];

        t_1018[k] = f_3 * pc_x[k] * smg_730[k];
    }

#pragma omp simd aligned(t_1019, t_1020, t_1021, t_1022, t_1023, pc_x, pc_y, slg_595, \
                         smf0_486, smf1_486, smg_730, smg_731, smg_732, smg_733, \
                         smg_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1019[k] = f_3 * pc_x[k] * smg_731[k];

        t_1020[k] = f_3 * pc_x[k] * smg_732[k];

        t_1021[k] = f_3 * pc_x[k] * smg_733[k];

        t_1022[k] = f_3 * pc_x[k] * smg_734[k];

        t_1023[k] = f_14 * slg_595[k]
                    + f_1 * smf0_486[k]
                    - f_2 * smf1_486[k]
                    + f_3 * pc_y[k] * smg_730[k];
    }

#pragma omp simd aligned(t_1024, t_1025, t_1026, pc_y, pc_z, slg_580, slg_597, slg_598, \
                         smf0_488, smf0_489, smf1_488, smf1_489, smg_730, smg_732, \
                         smg_733 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1024[k] = f_11 * slg_580[k]
                    + f_3 * pc_z[k] * smg_730[k];

        t_1025[k] = f_14 * slg_597[k]
                    + f_4 * smf0_488[k]
                    - f_5 * smf1_488[k]
                    + f_3 * pc_y[k] * smg_732[k];

        t_1026[k] = f_14 * slg_598[k]
                    + f_6 * smf0_489[k]
                    - f_7 * smf1_489[k]
                    + f_3 * pc_y[k] * smg_733[k];
    }

#pragma omp simd aligned(t_1027, t_1028, t_1029, t_1030, pc_x, pc_y, pc_z, slg_584, slg_599, \
                         slg_600, smf0_489, smf0_490, smf1_489, smf1_490, smg_734, \
                         smg_735 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1027[k] = f_14 * slg_599[k]
                    + f_3 * pc_y[k] * smg_734[k];

        t_1028[k] = f_11 * slg_584[k]
                    + f_1 * smf0_489[k]
                    - f_2 * smf1_489[k]
                    + f_3 * pc_z[k] * smg_734[k];

        t_1029[k] = f_1 * smf0_490[k]
                    - f_2 * smf1_490[k]
                    + f_3 * pc_x[k] * smg_735[k];

        t_1030[k] = f_15 * slg_600[k]
                    + f_3 * pc_y[k] * smg_735[k];
    }

#pragma omp simd aligned(t_1031, t_1032, t_1033, pc_x, pc_y, pc_z, slg_585, slg_602, smf0_493, \
                         smf1_493, smg_735, smg_737, smg_738 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1031[k] = f_16 * slg_585[k]
                    + f_3 * pc_z[k] * smg_735[k];

        t_1032[k] = f_4 * smf0_493[k]
                    - f_5 * smf1_493[k]
                    + f_3 * pc_x[k] * smg_738[k];

        t_1033[k] = f_15 * slg_602[k]
                    + f_3 * pc_y[k] * smg_737[k];
    }

#pragma omp simd aligned(t_1034, t_1035, t_1036, t_1037, pc_x, pc_y, pc_z, slg_588, slg_605, \
                         smf0_495, smf0_496, smf1_495, smf1_496, smg_738, smg_740, \
                         smg_741 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1034[k] = f_4 * smf0_495[k]
                    - f_5 * smf1_495[k]
                    + f_3 * pc_x[k] * smg_740[k];

        t_1035[k] = f_6 * smf0_496[k]
                    - f_7 * smf1_496[k]
                    + f_3 * pc_x[k] * smg_741[k];

        t_1036[k] = f_16 * slg_588[k]
                    + f_3 * pc_z[k] * smg_738[k];

        t_1037[k] = f_15 * slg_605[k]
                    + f_3 * pc_y[k] * smg_740[k];
    }

#pragma omp simd aligned(t_1038, t_1039, t_1040, t_1041, t_1042, t_1043, pc_x, smf0_499, \
                         smf1_499, smg_744, smg_745, smg_746, smg_747, smg_748, \
                         smg_749 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1038[k] = f_6 * smf0_499[k]
                    - f_7 * smf1_499[k]
                    + f_3 * pc_x[k] * smg_744[k];

        t_1039[k] = f_3 * pc_x[k] * smg_745[k];

        t_1040[k] = f_3 * pc_x[k] * smg_746[k];

        t_1041[k] = f_3 * pc_x[k] * smg_747[k];

        t_1042[k] = f_3 * pc_x[k] * smg_748[k];

        t_1043[k] = f_3 * pc_x[k] * smg_749[k];
    }

#pragma omp simd aligned(t_1044, t_1045, t_1046, pc_y, pc_z, slg_595, slg_610, slg_612, \
                         smf0_496, smf0_498, smf1_496, smf1_498, smg_745, \
                         smg_747 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1044[k] = f_15 * slg_610[k]
                    + f_1 * smf0_496[k]
                    - f_2 * smf1_496[k]
                    + f_3 * pc_y[k] * smg_745[k];

        t_1045[k] = f_16 * slg_595[k]
                    + f_3 * pc_z[k] * smg_745[k];

        t_1046[k] = f_15 * slg_612[k]
                    + f_4 * smf0_498[k]
                    - f_5 * smf1_498[k]
                    + f_3 * pc_y[k] * smg_747[k];
    }

#pragma omp simd aligned(t_1047, t_1048, t_1049, pc_y, pc_z, slg_599, slg_613, slg_614, \
                         smf0_499, smf1_499, smg_748, smg_749 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1047[k] = f_15 * slg_613[k]
                    + f_6 * smf0_499[k]
                    - f_7 * smf1_499[k]
                    + f_3 * pc_y[k] * smg_748[k];

        t_1048[k] = f_15 * slg_614[k]
                    + f_3 * pc_y[k] * smg_749[k];

        t_1049[k] = f_16 * slg_599[k]
                    + f_1 * smf0_499[k]
                    - f_2 * smf1_499[k]
                    + f_3 * pc_z[k] * smg_749[k];
    }

#pragma omp simd aligned(t_1050, t_1051, t_1052, t_1053, pc_x, pc_y, pc_z, slg_600, slg_615, \
                         smf0_500, smf0_503, smf1_500, smf1_503, smg_750, \
                         smg_753 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1050[k] = f_1 * smf0_500[k]
                    - f_2 * smf1_500[k]
                    + f_3 * pc_x[k] * smg_750[k];

        t_1051[k] = f_16 * slg_615[k]
                    + f_3 * pc_y[k] * smg_750[k];

        t_1052[k] = f_15 * slg_600[k]
                    + f_3 * pc_z[k] * smg_750[k];

        t_1053[k] = f_4 * smf0_503[k]
                    - f_5 * smf1_503[k]
                    + f_3 * pc_x[k] * smg_753[k];
    }

#pragma omp simd aligned(t_1054, t_1055, t_1056, pc_x, pc_y, slg_617, smf0_505, smf0_506, \
                         smf1_505, smf1_506, smg_752, smg_755, \
                         smg_756 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1054[k] = f_16 * slg_617[k]
                    + f_3 * pc_y[k] * smg_752[k];

        t_1055[k] = f_4 * smf0_505[k]
                    - f_5 * smf1_505[k]
                    + f_3 * pc_x[k] * smg_755[k];

        t_1056[k] = f_6 * smf0_506[k]
                    - f_7 * smf1_506[k]
                    + f_3 * pc_x[k] * smg_756[k];
    }

#pragma omp simd aligned(t_1057, t_1058, t_1059, t_1060, pc_x, pc_y, pc_z, slg_603, slg_620, \
                         smf0_509, smf1_509, smg_753, smg_755, smg_759, \
                         smg_760 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1057[k] = f_15 * slg_603[k]
                    + f_3 * pc_z[k] * smg_753[k];

        t_1058[k] = f_16 * slg_620[k]
                    + f_3 * pc_y[k] * smg_755[k];

        t_1059[k] = f_6 * smf0_509[k]
                    - f_7 * smf1_509[k]
                    + f_3 * pc_x[k] * smg_759[k];

        t_1060[k] = f_3 * pc_x[k] * smg_760[k];
    }
}

static auto
compute_prim_smh_three_center_electron_repulsion_0_piece9(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t slh0,
                                                          const size_t slg, const size_t slh1,
                                                          const size_t smf0, const size_t smf1,
                                                          const size_t smg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / q;
    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.0 / gamma;
    const auto f_5 = p / (gamma * q);
    const auto f_6 = 0.5 / gamma;
    const auto f_7 = 0.5 * p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_12 = 4.0 / q;
    const auto f_13 = 3.5 / q;
    const auto f_14 = 3.0 / q;
    const auto f_15 = 2.5 / q;
    const auto f_16 = 2.0 / q;

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

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *slh0_924 = buffer.data(slh0 + 924);
    const auto *slh0_929 = buffer.data(slh0 + 929);
    const auto *slh0_933 = buffer.data(slh0 + 933);
    const auto *slh0_939 = buffer.data(slh0 + 939);
    const auto *slh0_941 = buffer.data(slh0 + 941);
    const auto *slh0_942 = buffer.data(slh0 + 942);
    const auto *slh0_944 = buffer.data(slh0 + 944);

    const auto *slg_610 = buffer.data(slg + 610);
    const auto *slg_614 = buffer.data(slg + 614);
    const auto *slg_615 = buffer.data(slg + 615);
    const auto *slg_618 = buffer.data(slg + 618);
    const auto *slg_625 = buffer.data(slg + 625);
    const auto *slg_627 = buffer.data(slg + 627);
    const auto *slg_628 = buffer.data(slg + 628);
    const auto *slg_629 = buffer.data(slg + 629);
    const auto *slg_630 = buffer.data(slg + 630);
    const auto *slg_632 = buffer.data(slg + 632);
    const auto *slg_633 = buffer.data(slg + 633);
    const auto *slg_635 = buffer.data(slg + 635);
    const auto *slg_640 = buffer.data(slg + 640);
    const auto *slg_642 = buffer.data(slg + 642);
    const auto *slg_643 = buffer.data(slg + 643);
    const auto *slg_644 = buffer.data(slg + 644);
    const auto *slg_645 = buffer.data(slg + 645);
    const auto *slg_647 = buffer.data(slg + 647);
    const auto *slg_648 = buffer.data(slg + 648);
    const auto *slg_650 = buffer.data(slg + 650);
    const auto *slg_655 = buffer.data(slg + 655);
    const auto *slg_657 = buffer.data(slg + 657);
    const auto *slg_658 = buffer.data(slg + 658);
    const auto *slg_659 = buffer.data(slg + 659);
    const auto *slg_660 = buffer.data(slg + 660);
    const auto *slg_662 = buffer.data(slg + 662);
    const auto *slg_663 = buffer.data(slg + 663);
    const auto *slg_665 = buffer.data(slg + 665);
    const auto *slg_670 = buffer.data(slg + 670);
    const auto *slg_672 = buffer.data(slg + 672);
    const auto *slg_673 = buffer.data(slg + 673);
    const auto *slg_674 = buffer.data(slg + 674);

    const auto *slh1_924 = buffer.data(slh1 + 924);
    const auto *slh1_929 = buffer.data(slh1 + 929);
    const auto *slh1_933 = buffer.data(slh1 + 933);
    const auto *slh1_939 = buffer.data(slh1 + 939);
    const auto *slh1_941 = buffer.data(slh1 + 941);
    const auto *slh1_942 = buffer.data(slh1 + 942);
    const auto *slh1_944 = buffer.data(slh1 + 944);

    const auto *smf0_506 = buffer.data(smf0 + 506);
    const auto *smf0_508 = buffer.data(smf0 + 508);
    const auto *smf0_509 = buffer.data(smf0 + 509);
    const auto *smf0_510 = buffer.data(smf0 + 510);
    const auto *smf0_513 = buffer.data(smf0 + 513);
    const auto *smf0_515 = buffer.data(smf0 + 515);
    const auto *smf0_516 = buffer.data(smf0 + 516);
    const auto *smf0_518 = buffer.data(smf0 + 518);
    const auto *smf0_519 = buffer.data(smf0 + 519);
    const auto *smf0_520 = buffer.data(smf0 + 520);
    const auto *smf0_523 = buffer.data(smf0 + 523);
    const auto *smf0_525 = buffer.data(smf0 + 525);
    const auto *smf0_526 = buffer.data(smf0 + 526);
    const auto *smf0_528 = buffer.data(smf0 + 528);
    const auto *smf0_529 = buffer.data(smf0 + 529);
    const auto *smf0_533 = buffer.data(smf0 + 533);
    const auto *smf0_536 = buffer.data(smf0 + 536);
    const auto *smf0_540 = buffer.data(smf0 + 540);
    const auto *smf0_543 = buffer.data(smf0 + 543);
    const auto *smf0_545 = buffer.data(smf0 + 545);
    const auto *smf0_546 = buffer.data(smf0 + 546);
    const auto *smf0_548 = buffer.data(smf0 + 548);
    const auto *smf0_549 = buffer.data(smf0 + 549);

    const auto *smf1_506 = buffer.data(smf1 + 506);
    const auto *smf1_508 = buffer.data(smf1 + 508);
    const auto *smf1_509 = buffer.data(smf1 + 509);
    const auto *smf1_510 = buffer.data(smf1 + 510);
    const auto *smf1_513 = buffer.data(smf1 + 513);
    const auto *smf1_515 = buffer.data(smf1 + 515);
    const auto *smf1_516 = buffer.data(smf1 + 516);
    const auto *smf1_518 = buffer.data(smf1 + 518);
    const auto *smf1_519 = buffer.data(smf1 + 519);
    const auto *smf1_520 = buffer.data(smf1 + 520);
    const auto *smf1_523 = buffer.data(smf1 + 523);
    const auto *smf1_525 = buffer.data(smf1 + 525);
    const auto *smf1_526 = buffer.data(smf1 + 526);
    const auto *smf1_528 = buffer.data(smf1 + 528);
    const auto *smf1_529 = buffer.data(smf1 + 529);
    const auto *smf1_533 = buffer.data(smf1 + 533);
    const auto *smf1_536 = buffer.data(smf1 + 536);
    const auto *smf1_540 = buffer.data(smf1 + 540);
    const auto *smf1_543 = buffer.data(smf1 + 543);
    const auto *smf1_545 = buffer.data(smf1 + 545);
    const auto *smf1_546 = buffer.data(smf1 + 546);
    const auto *smf1_548 = buffer.data(smf1 + 548);
    const auto *smf1_549 = buffer.data(smf1 + 549);

    const auto *smg_760 = buffer.data(smg + 760);
    const auto *smg_761 = buffer.data(smg + 761);
    const auto *smg_762 = buffer.data(smg + 762);
    const auto *smg_763 = buffer.data(smg + 763);
    const auto *smg_764 = buffer.data(smg + 764);
    const auto *smg_765 = buffer.data(smg + 765);
    const auto *smg_767 = buffer.data(smg + 767);
    const auto *smg_768 = buffer.data(smg + 768);
    const auto *smg_770 = buffer.data(smg + 770);
    const auto *smg_771 = buffer.data(smg + 771);
    const auto *smg_774 = buffer.data(smg + 774);
    const auto *smg_775 = buffer.data(smg + 775);
    const auto *smg_776 = buffer.data(smg + 776);
    const auto *smg_777 = buffer.data(smg + 777);
    const auto *smg_778 = buffer.data(smg + 778);
    const auto *smg_779 = buffer.data(smg + 779);
    const auto *smg_780 = buffer.data(smg + 780);
    const auto *smg_782 = buffer.data(smg + 782);
    const auto *smg_783 = buffer.data(smg + 783);
    const auto *smg_785 = buffer.data(smg + 785);
    const auto *smg_786 = buffer.data(smg + 786);
    const auto *smg_789 = buffer.data(smg + 789);
    const auto *smg_790 = buffer.data(smg + 790);
    const auto *smg_791 = buffer.data(smg + 791);
    const auto *smg_792 = buffer.data(smg + 792);
    const auto *smg_793 = buffer.data(smg + 793);
    const auto *smg_794 = buffer.data(smg + 794);
    const auto *smg_795 = buffer.data(smg + 795);
    const auto *smg_797 = buffer.data(smg + 797);
    const auto *smg_798 = buffer.data(smg + 798);
    const auto *smg_800 = buffer.data(smg + 800);
    const auto *smg_801 = buffer.data(smg + 801);
    const auto *smg_805 = buffer.data(smg + 805);
    const auto *smg_806 = buffer.data(smg + 806);
    const auto *smg_807 = buffer.data(smg + 807);
    const auto *smg_808 = buffer.data(smg + 808);
    const auto *smg_809 = buffer.data(smg + 809);
    const auto *smg_810 = buffer.data(smg + 810);
    const auto *smg_812 = buffer.data(smg + 812);
    const auto *smg_813 = buffer.data(smg + 813);
    const auto *smg_815 = buffer.data(smg + 815);
    const auto *smg_816 = buffer.data(smg + 816);
    const auto *smg_819 = buffer.data(smg + 819);
    const auto *smg_820 = buffer.data(smg + 820);
    const auto *smg_821 = buffer.data(smg + 821);
    const auto *smg_822 = buffer.data(smg + 822);
    const auto *smg_823 = buffer.data(smg + 823);
    const auto *smg_824 = buffer.data(smg + 824);

#pragma omp simd aligned(t_1061, t_1062, t_1063, t_1064, t_1065, pc_x, pc_y, slg_625, \
                         smf0_506, smf1_506, smg_760, smg_761, smg_762, smg_763, \
                         smg_764 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1061[k] = f_3 * pc_x[k] * smg_761[k];

        t_1062[k] = f_3 * pc_x[k] * smg_762[k];

        t_1063[k] = f_3 * pc_x[k] * smg_763[k];

        t_1064[k] = f_3 * pc_x[k] * smg_764[k];

        t_1065[k] = f_16 * slg_625[k]
                    + f_1 * smf0_506[k]
                    - f_2 * smf1_506[k]
                    + f_3 * pc_y[k] * smg_760[k];
    }

#pragma omp simd aligned(t_1066, t_1067, t_1068, pc_y, pc_z, slg_610, slg_627, slg_628, \
                         smf0_508, smf0_509, smf1_508, smf1_509, smg_760, smg_762, \
                         smg_763 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1066[k] = f_15 * slg_610[k]
                    + f_3 * pc_z[k] * smg_760[k];

        t_1067[k] = f_16 * slg_627[k]
                    + f_4 * smf0_508[k]
                    - f_5 * smf1_508[k]
                    + f_3 * pc_y[k] * smg_762[k];

        t_1068[k] = f_16 * slg_628[k]
                    + f_6 * smf0_509[k]
                    - f_7 * smf1_509[k]
                    + f_3 * pc_y[k] * smg_763[k];
    }

#pragma omp simd aligned(t_1069, t_1070, t_1071, t_1072, pc_x, pc_y, pc_z, slg_614, slg_629, \
                         slg_630, smf0_509, smf0_510, smf1_509, smf1_510, smg_764, \
                         smg_765 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1069[k] = f_16 * slg_629[k]
                    + f_3 * pc_y[k] * smg_764[k];

        t_1070[k] = f_15 * slg_614[k]
                    + f_1 * smf0_509[k]
                    - f_2 * smf1_509[k]
                    + f_3 * pc_z[k] * smg_764[k];

        t_1071[k] = f_1 * smf0_510[k]
                    - f_2 * smf1_510[k]
                    + f_3 * pc_x[k] * smg_765[k];

        t_1072[k] = f_11 * slg_630[k]
                    + f_3 * pc_y[k] * smg_765[k];
    }

#pragma omp simd aligned(t_1073, t_1074, t_1075, pc_x, pc_y, pc_z, slg_615, slg_632, smf0_513, \
                         smf1_513, smg_765, smg_767, smg_768 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1073[k] = f_14 * slg_615[k]
                    + f_3 * pc_z[k] * smg_765[k];

        t_1074[k] = f_4 * smf0_513[k]
                    - f_5 * smf1_513[k]
                    + f_3 * pc_x[k] * smg_768[k];

        t_1075[k] = f_11 * slg_632[k]
                    + f_3 * pc_y[k] * smg_767[k];
    }

#pragma omp simd aligned(t_1076, t_1077, t_1078, t_1079, pc_x, pc_y, pc_z, slg_618, slg_635, \
                         smf0_515, smf0_516, smf1_515, smf1_516, smg_768, smg_770, \
                         smg_771 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1076[k] = f_4 * smf0_515[k]
                    - f_5 * smf1_515[k]
                    + f_3 * pc_x[k] * smg_770[k];

        t_1077[k] = f_6 * smf0_516[k]
                    - f_7 * smf1_516[k]
                    + f_3 * pc_x[k] * smg_771[k];

        t_1078[k] = f_14 * slg_618[k]
                    + f_3 * pc_z[k] * smg_768[k];

        t_1079[k] = f_11 * slg_635[k]
                    + f_3 * pc_y[k] * smg_770[k];
    }

#pragma omp simd aligned(t_1080, t_1081, t_1082, t_1083, t_1084, t_1085, pc_x, smf0_519, \
                         smf1_519, smg_774, smg_775, smg_776, smg_777, smg_778, \
                         smg_779 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1080[k] = f_6 * smf0_519[k]
                    - f_7 * smf1_519[k]
                    + f_3 * pc_x[k] * smg_774[k];

        t_1081[k] = f_3 * pc_x[k] * smg_775[k];

        t_1082[k] = f_3 * pc_x[k] * smg_776[k];

        t_1083[k] = f_3 * pc_x[k] * smg_777[k];

        t_1084[k] = f_3 * pc_x[k] * smg_778[k];

        t_1085[k] = f_3 * pc_x[k] * smg_779[k];
    }

#pragma omp simd aligned(t_1086, t_1087, t_1088, pc_y, pc_z, slg_625, slg_640, slg_642, \
                         smf0_516, smf0_518, smf1_516, smf1_518, smg_775, \
                         smg_777 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1086[k] = f_11 * slg_640[k]
                    + f_1 * smf0_516[k]
                    - f_2 * smf1_516[k]
                    + f_3 * pc_y[k] * smg_775[k];

        t_1087[k] = f_14 * slg_625[k]
                    + f_3 * pc_z[k] * smg_775[k];

        t_1088[k] = f_11 * slg_642[k]
                    + f_4 * smf0_518[k]
                    - f_5 * smf1_518[k]
                    + f_3 * pc_y[k] * smg_777[k];
    }

#pragma omp simd aligned(t_1089, t_1090, t_1091, pc_y, pc_z, slg_629, slg_643, slg_644, \
                         smf0_519, smf1_519, smg_778, smg_779 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1089[k] = f_11 * slg_643[k]
                    + f_6 * smf0_519[k]
                    - f_7 * smf1_519[k]
                    + f_3 * pc_y[k] * smg_778[k];

        t_1090[k] = f_11 * slg_644[k]
                    + f_3 * pc_y[k] * smg_779[k];

        t_1091[k] = f_14 * slg_629[k]
                    + f_1 * smf0_519[k]
                    - f_2 * smf1_519[k]
                    + f_3 * pc_z[k] * smg_779[k];
    }

#pragma omp simd aligned(t_1092, t_1093, t_1094, t_1095, pc_x, pc_y, pc_z, slg_630, slg_645, \
                         smf0_520, smf0_523, smf1_520, smf1_523, smg_780, \
                         smg_783 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1092[k] = f_1 * smf0_520[k]
                    - f_2 * smf1_520[k]
                    + f_3 * pc_x[k] * smg_780[k];

        t_1093[k] = f_10 * slg_645[k]
                    + f_3 * pc_y[k] * smg_780[k];

        t_1094[k] = f_13 * slg_630[k]
                    + f_3 * pc_z[k] * smg_780[k];

        t_1095[k] = f_4 * smf0_523[k]
                    - f_5 * smf1_523[k]
                    + f_3 * pc_x[k] * smg_783[k];
    }

#pragma omp simd aligned(t_1096, t_1097, t_1098, pc_x, pc_y, slg_647, smf0_525, smf0_526, \
                         smf1_525, smf1_526, smg_782, smg_785, \
                         smg_786 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1096[k] = f_10 * slg_647[k]
                    + f_3 * pc_y[k] * smg_782[k];

        t_1097[k] = f_4 * smf0_525[k]
                    - f_5 * smf1_525[k]
                    + f_3 * pc_x[k] * smg_785[k];

        t_1098[k] = f_6 * smf0_526[k]
                    - f_7 * smf1_526[k]
                    + f_3 * pc_x[k] * smg_786[k];
    }

#pragma omp simd aligned(t_1099, t_1100, t_1101, t_1102, pc_x, pc_y, pc_z, slg_633, slg_650, \
                         smf0_529, smf1_529, smg_783, smg_785, smg_789, \
                         smg_790 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1099[k] = f_13 * slg_633[k]
                    + f_3 * pc_z[k] * smg_783[k];

        t_1100[k] = f_10 * slg_650[k]
                    + f_3 * pc_y[k] * smg_785[k];

        t_1101[k] = f_6 * smf0_529[k]
                    - f_7 * smf1_529[k]
                    + f_3 * pc_x[k] * smg_789[k];

        t_1102[k] = f_3 * pc_x[k] * smg_790[k];
    }

#pragma omp simd aligned(t_1103, t_1104, t_1105, t_1106, t_1107, pc_x, pc_y, slg_655, \
                         smf0_526, smf1_526, smg_790, smg_791, smg_792, smg_793, \
                         smg_794 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1103[k] = f_3 * pc_x[k] * smg_791[k];

        t_1104[k] = f_3 * pc_x[k] * smg_792[k];

        t_1105[k] = f_3 * pc_x[k] * smg_793[k];

        t_1106[k] = f_3 * pc_x[k] * smg_794[k];

        t_1107[k] = f_10 * slg_655[k]
                    + f_1 * smf0_526[k]
                    - f_2 * smf1_526[k]
                    + f_3 * pc_y[k] * smg_790[k];
    }

#pragma omp simd aligned(t_1108, t_1109, t_1110, pc_y, pc_z, slg_640, slg_657, slg_658, \
                         smf0_528, smf0_529, smf1_528, smf1_529, smg_790, smg_792, \
                         smg_793 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1108[k] = f_13 * slg_640[k]
                    + f_3 * pc_z[k] * smg_790[k];

        t_1109[k] = f_10 * slg_657[k]
                    + f_4 * smf0_528[k]
                    - f_5 * smf1_528[k]
                    + f_3 * pc_y[k] * smg_792[k];

        t_1110[k] = f_10 * slg_658[k]
                    + f_6 * smf0_529[k]
                    - f_7 * smf1_529[k]
                    + f_3 * pc_y[k] * smg_793[k];
    }

#pragma omp simd aligned(t_1111, t_1112, t_1113, t_1114, pb_y, pc_y, pc_z, slh0_924, slg_644, \
                         slg_659, slg_660, slh1_924, smf0_529, smf1_529, smg_794, \
                         smg_795 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1111[k] = f_10 * slg_659[k]
                    + f_3 * pc_y[k] * smg_794[k];

        t_1112[k] = f_13 * slg_644[k]
                    + f_1 * smf0_529[k]
                    - f_2 * smf1_529[k]
                    + f_3 * pc_z[k] * smg_794[k];

        t_1113[k] = pb_y[k] * slh0_924[k]
                    - f_8 * pc_y[k] * slh1_924[k];

        t_1114[k] = f_9 * slg_660[k]
                    + f_3 * pc_y[k] * smg_795[k];
    }

#pragma omp simd aligned(t_1115, t_1116, t_1117, pc_x, pc_y, pc_z, slg_645, slg_662, smf0_533, \
                         smf1_533, smg_795, smg_797, smg_798 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1115[k] = f_12 * slg_645[k]
                    + f_3 * pc_z[k] * smg_795[k];

        t_1116[k] = f_4 * smf0_533[k]
                    - f_5 * smf1_533[k]
                    + f_3 * pc_x[k] * smg_798[k];

        t_1117[k] = f_9 * slg_662[k]
                    + f_3 * pc_y[k] * smg_797[k];
    }

#pragma omp simd aligned(t_1118, t_1119, t_1120, pb_y, pc_x, pc_y, pc_z, slh0_929, slg_648, \
                         slh1_929, smf0_536, smf1_536, smg_798, \
                         smg_801 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1118[k] = pb_y[k] * slh0_929[k]
                    - f_8 * pc_y[k] * slh1_929[k];

        t_1119[k] = f_6 * smf0_536[k]
                    - f_7 * smf1_536[k]
                    + f_3 * pc_x[k] * smg_801[k];

        t_1120[k] = f_12 * slg_648[k]
                    + f_3 * pc_z[k] * smg_798[k];
    }

#pragma omp simd aligned(t_1121, t_1122, t_1123, t_1124, t_1125, pb_y, pc_x, pc_y, slh0_933, \
                         slg_665, slh1_933, smg_800, smg_805, smg_806, \
                         smg_807 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1121[k] = f_9 * slg_665[k]
                    + f_3 * pc_y[k] * smg_800[k];

        t_1122[k] = pb_y[k] * slh0_933[k]
                    - f_8 * pc_y[k] * slh1_933[k];

        t_1123[k] = f_3 * pc_x[k] * smg_805[k];

        t_1124[k] = f_3 * pc_x[k] * smg_806[k];

        t_1125[k] = f_3 * pc_x[k] * smg_807[k];
    }

#pragma omp simd aligned(t_1126, t_1127, t_1128, t_1129, pb_y, pc_x, pc_y, pc_z, slh0_939, \
                         slg_655, slg_670, slh1_939, smg_805, smg_808, \
                         smg_809 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1126[k] = f_3 * pc_x[k] * smg_808[k];

        t_1127[k] = f_3 * pc_x[k] * smg_809[k];

        t_1128[k] = pb_y[k] * slh0_939[k]
                    + f_15 * slg_670[k]
                    - f_8 * pc_y[k] * slh1_939[k];

        t_1129[k] = f_12 * slg_655[k]
                    + f_3 * pc_z[k] * smg_805[k];
    }

#pragma omp simd aligned(t_1130, t_1131, t_1132, t_1133, pb_y, pc_y, slh0_941, slh0_942, \
                         slh0_944, slg_672, slg_673, slg_674, slh1_941, slh1_942, slh1_944, \
                         smg_809 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1130[k] = pb_y[k] * slh0_941[k]
                    + f_11 * slg_672[k]
                    - f_8 * pc_y[k] * slh1_941[k];

        t_1131[k] = pb_y[k] * slh0_942[k]
                    + f_10 * slg_673[k]
                    - f_8 * pc_y[k] * slh1_942[k];

        t_1132[k] = f_9 * slg_674[k]
                    + f_3 * pc_y[k] * smg_809[k];

        t_1133[k] = pb_y[k] * slh0_944[k]
                    - f_8 * pc_y[k] * slh1_944[k];
    }

#pragma omp simd aligned(t_1134, t_1135, t_1136, t_1137, t_1138, pc_x, pc_y, pc_z, slg_660, \
                         smf0_540, smf0_543, smf1_540, smf1_543, smg_810, smg_812, \
                         smg_813 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1134[k] = f_1 * smf0_540[k]
                    - f_2 * smf1_540[k]
                    + f_3 * pc_x[k] * smg_810[k];

        t_1135[k] = f_3 * pc_y[k] * smg_810[k];

        t_1136[k] = f_0 * slg_660[k]
                    + f_3 * pc_z[k] * smg_810[k];

        t_1137[k] = f_4 * smf0_543[k]
                    - f_5 * smf1_543[k]
                    + f_3 * pc_x[k] * smg_813[k];

        t_1138[k] = f_3 * pc_y[k] * smg_812[k];
    }

#pragma omp simd aligned(t_1139, t_1140, t_1141, t_1142, pc_x, pc_y, pc_z, slg_663, smf0_545, \
                         smf0_546, smf1_545, smf1_546, smg_813, smg_815, \
                         smg_816 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1139[k] = f_4 * smf0_545[k]
                    - f_5 * smf1_545[k]
                    + f_3 * pc_x[k] * smg_815[k];

        t_1140[k] = f_6 * smf0_546[k]
                    - f_7 * smf1_546[k]
                    + f_3 * pc_x[k] * smg_816[k];

        t_1141[k] = f_0 * slg_663[k]
                    + f_3 * pc_z[k] * smg_813[k];

        t_1142[k] = f_3 * pc_y[k] * smg_815[k];
    }

#pragma omp simd aligned(t_1143, t_1144, t_1145, t_1146, t_1147, t_1148, pc_x, smf0_549, \
                         smf1_549, smg_819, smg_820, smg_821, smg_822, smg_823, \
                         smg_824 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1143[k] = f_6 * smf0_549[k]
                    - f_7 * smf1_549[k]
                    + f_3 * pc_x[k] * smg_819[k];

        t_1144[k] = f_3 * pc_x[k] * smg_820[k];

        t_1145[k] = f_3 * pc_x[k] * smg_821[k];

        t_1146[k] = f_3 * pc_x[k] * smg_822[k];

        t_1147[k] = f_3 * pc_x[k] * smg_823[k];

        t_1148[k] = f_3 * pc_x[k] * smg_824[k];
    }

#pragma omp simd aligned(t_1149, t_1150, t_1151, t_1152, pc_y, pc_z, slg_670, smf0_546, \
                         smf0_548, smf0_549, smf1_546, smf1_548, smf1_549, smg_820, smg_822, \
                         smg_823 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1149[k] = f_1 * smf0_546[k]
                    - f_2 * smf1_546[k]
                    + f_3 * pc_y[k] * smg_820[k];

        t_1150[k] = f_0 * slg_670[k]
                    + f_3 * pc_z[k] * smg_820[k];

        t_1151[k] = f_4 * smf0_548[k]
                    - f_5 * smf1_548[k]
                    + f_3 * pc_y[k] * smg_822[k];

        t_1152[k] = f_6 * smf0_549[k]
                    - f_7 * smf1_549[k]
                    + f_3 * pc_y[k] * smg_823[k];
    }

#pragma omp simd aligned(t_1153, t_1154, pc_y, pc_z, slg_674, smf0_549, smf1_549, \
                         smg_824 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1153[k] = f_3 * pc_y[k] * smg_824[k];

        t_1154[k] = f_0 * slg_674[k]
                    + f_1 * smf0_549[k]
                    - f_2 * smf1_549[k]
                    + f_3 * pc_z[k] * smg_824[k];
    }
}

auto
compute_prim_smh_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t slh0, const size_t slg,
                                                   const size_t slh1, const size_t smf0,
                                                   const size_t smf1, const size_t smg,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_smh_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, slh0, slg,
                                                              slh1, smf0, smf1, smg, ncols,
                                                              gamma, p, q);

    compute_prim_smh_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, slh0, slg,
                                                              slh1, smf0, smf1, smg, ncols,
                                                              gamma, p, q);

    compute_prim_smh_three_center_electron_repulsion_0_piece2(buffer, target, pb, pc, slh0, slg,
                                                              slh1, smf0, smf1, smg, ncols,
                                                              gamma, p, q);

    compute_prim_smh_three_center_electron_repulsion_0_piece3(buffer, target, pb, pc, slh0, slg,
                                                              slh1, smf0, smf1, smg, ncols,
                                                              gamma, p, q);

    compute_prim_smh_three_center_electron_repulsion_0_piece4(buffer, target, pb, pc, slh0, slg,
                                                              slh1, smf0, smf1, smg, ncols,
                                                              gamma, p, q);

    compute_prim_smh_three_center_electron_repulsion_0_piece5(buffer, target, pb, pc, slh0, slg,
                                                              slh1, smf0, smf1, smg, ncols,
                                                              gamma, p, q);

    compute_prim_smh_three_center_electron_repulsion_0_piece6(buffer, target, pb, pc, slh0, slg,
                                                              slh1, smf0, smf1, smg, ncols,
                                                              gamma, p, q);

    compute_prim_smh_three_center_electron_repulsion_0_piece7(buffer, target, pb, pc, slh0, slg,
                                                              slh1, smg, ncols, gamma, p, q);

    compute_prim_smh_three_center_electron_repulsion_0_piece8(buffer, target, pb, pc, slh0, slg,
                                                              slh1, smf0, smf1, smg, ncols,
                                                              gamma, p, q);

    compute_prim_smh_three_center_electron_repulsion_0_piece9(buffer, target, pb, pc, slh0, slg,
                                                              slh1, smf0, smf1, smg, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
