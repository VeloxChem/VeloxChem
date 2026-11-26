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


#include "SimdThreeCenterElectronRepulsionVrrRecSNH.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_snh_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t smh0,
                                                          const size_t smg, const size_t smh1,
                                                          const size_t snf0, const size_t snf1,
                                                          const size_t sng, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.0 / q;
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
    const auto f_12 = 4.5 / q;
    const auto f_13 = 4.0 / q;

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

    const auto *smh0_0 = buffer.data(smh0 + 0);
    const auto *smh0_3 = buffer.data(smh0 + 3);
    const auto *smh0_5 = buffer.data(smh0 + 5);
    const auto *smh0_6 = buffer.data(smh0 + 6);
    const auto *smh0_9 = buffer.data(smh0 + 9);
    const auto *smh0_15 = buffer.data(smh0 + 15);
    const auto *smh0_20 = buffer.data(smh0 + 20);
    const auto *smh0_24 = buffer.data(smh0 + 24);
    const auto *smh0_27 = buffer.data(smh0 + 27);
    const auto *smh0_36 = buffer.data(smh0 + 36);
    const auto *smh0_42 = buffer.data(smh0 + 42);
    const auto *smh0_47 = buffer.data(smh0 + 47);
    const auto *smh0_51 = buffer.data(smh0 + 51);
    const auto *smh0_62 = buffer.data(smh0 + 62);

    const auto *smg_0 = buffer.data(smg + 0);
    const auto *smg_1 = buffer.data(smg + 1);
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
    const auto *smg_48 = buffer.data(smg + 48);
    const auto *smg_50 = buffer.data(smg + 50);
    const auto *smg_51 = buffer.data(smg + 51);
    const auto *smg_54 = buffer.data(smg + 54);
    const auto *smg_55 = buffer.data(smg + 55);
    const auto *smg_56 = buffer.data(smg + 56);
    const auto *smg_57 = buffer.data(smg + 57);
    const auto *smg_58 = buffer.data(smg + 58);
    const auto *smg_59 = buffer.data(smg + 59);
    const auto *smg_70 = buffer.data(smg + 70);
    const auto *smg_71 = buffer.data(smg + 71);
    const auto *smg_72 = buffer.data(smg + 72);
    const auto *smg_73 = buffer.data(smg + 73);
    const auto *smg_74 = buffer.data(smg + 74);
    const auto *smg_75 = buffer.data(smg + 75);
    const auto *smg_78 = buffer.data(smg + 78);
    const auto *smg_80 = buffer.data(smg + 80);
    const auto *smg_81 = buffer.data(smg + 81);
    const auto *smg_84 = buffer.data(smg + 84);
    const auto *smg_85 = buffer.data(smg + 85);
    const auto *smg_86 = buffer.data(smg + 86);
    const auto *smg_87 = buffer.data(smg + 87);
    const auto *smg_88 = buffer.data(smg + 88);
    const auto *smg_89 = buffer.data(smg + 89);

    const auto *smh1_0 = buffer.data(smh1 + 0);
    const auto *smh1_3 = buffer.data(smh1 + 3);
    const auto *smh1_5 = buffer.data(smh1 + 5);
    const auto *smh1_6 = buffer.data(smh1 + 6);
    const auto *smh1_9 = buffer.data(smh1 + 9);
    const auto *smh1_15 = buffer.data(smh1 + 15);
    const auto *smh1_20 = buffer.data(smh1 + 20);
    const auto *smh1_24 = buffer.data(smh1 + 24);
    const auto *smh1_27 = buffer.data(smh1 + 27);
    const auto *smh1_36 = buffer.data(smh1 + 36);
    const auto *smh1_42 = buffer.data(smh1 + 42);
    const auto *smh1_47 = buffer.data(smh1 + 47);
    const auto *smh1_51 = buffer.data(smh1 + 51);
    const auto *smh1_62 = buffer.data(smh1 + 62);

    const auto *snf0_0 = buffer.data(snf0 + 0);
    const auto *snf0_3 = buffer.data(snf0 + 3);
    const auto *snf0_5 = buffer.data(snf0 + 5);
    const auto *snf0_6 = buffer.data(snf0 + 6);
    const auto *snf0_8 = buffer.data(snf0 + 8);
    const auto *snf0_9 = buffer.data(snf0 + 9);
    const auto *snf0_16 = buffer.data(snf0 + 16);
    const auto *snf0_18 = buffer.data(snf0 + 18);
    const auto *snf0_19 = buffer.data(snf0 + 19);
    const auto *snf0_28 = buffer.data(snf0 + 28);
    const auto *snf0_29 = buffer.data(snf0 + 29);
    const auto *snf0_30 = buffer.data(snf0 + 30);
    const auto *snf0_33 = buffer.data(snf0 + 33);
    const auto *snf0_35 = buffer.data(snf0 + 35);
    const auto *snf0_36 = buffer.data(snf0 + 36);
    const auto *snf0_38 = buffer.data(snf0 + 38);
    const auto *snf0_39 = buffer.data(snf0 + 39);
    const auto *snf0_48 = buffer.data(snf0 + 48);
    const auto *snf0_49 = buffer.data(snf0 + 49);
    const auto *snf0_50 = buffer.data(snf0 + 50);
    const auto *snf0_53 = buffer.data(snf0 + 53);
    const auto *snf0_55 = buffer.data(snf0 + 55);
    const auto *snf0_56 = buffer.data(snf0 + 56);
    const auto *snf0_58 = buffer.data(snf0 + 58);
    const auto *snf0_59 = buffer.data(snf0 + 59);

    const auto *snf1_0 = buffer.data(snf1 + 0);
    const auto *snf1_3 = buffer.data(snf1 + 3);
    const auto *snf1_5 = buffer.data(snf1 + 5);
    const auto *snf1_6 = buffer.data(snf1 + 6);
    const auto *snf1_8 = buffer.data(snf1 + 8);
    const auto *snf1_9 = buffer.data(snf1 + 9);
    const auto *snf1_16 = buffer.data(snf1 + 16);
    const auto *snf1_18 = buffer.data(snf1 + 18);
    const auto *snf1_19 = buffer.data(snf1 + 19);
    const auto *snf1_28 = buffer.data(snf1 + 28);
    const auto *snf1_29 = buffer.data(snf1 + 29);
    const auto *snf1_30 = buffer.data(snf1 + 30);
    const auto *snf1_33 = buffer.data(snf1 + 33);
    const auto *snf1_35 = buffer.data(snf1 + 35);
    const auto *snf1_36 = buffer.data(snf1 + 36);
    const auto *snf1_38 = buffer.data(snf1 + 38);
    const auto *snf1_39 = buffer.data(snf1 + 39);
    const auto *snf1_48 = buffer.data(snf1 + 48);
    const auto *snf1_49 = buffer.data(snf1 + 49);
    const auto *snf1_50 = buffer.data(snf1 + 50);
    const auto *snf1_53 = buffer.data(snf1 + 53);
    const auto *snf1_55 = buffer.data(snf1 + 55);
    const auto *snf1_56 = buffer.data(snf1 + 56);
    const auto *snf1_58 = buffer.data(snf1 + 58);
    const auto *snf1_59 = buffer.data(snf1 + 59);

    const auto *sng_0 = buffer.data(sng + 0);
    const auto *sng_2 = buffer.data(sng + 2);
    const auto *sng_3 = buffer.data(sng + 3);
    const auto *sng_5 = buffer.data(sng + 5);
    const auto *sng_6 = buffer.data(sng + 6);
    const auto *sng_9 = buffer.data(sng + 9);
    const auto *sng_10 = buffer.data(sng + 10);
    const auto *sng_11 = buffer.data(sng + 11);
    const auto *sng_12 = buffer.data(sng + 12);
    const auto *sng_13 = buffer.data(sng + 13);
    const auto *sng_14 = buffer.data(sng + 14);
    const auto *sng_15 = buffer.data(sng + 15);
    const auto *sng_17 = buffer.data(sng + 17);
    const auto *sng_18 = buffer.data(sng + 18);
    const auto *sng_20 = buffer.data(sng + 20);
    const auto *sng_25 = buffer.data(sng + 25);
    const auto *sng_26 = buffer.data(sng + 26);
    const auto *sng_27 = buffer.data(sng + 27);
    const auto *sng_28 = buffer.data(sng + 28);
    const auto *sng_29 = buffer.data(sng + 29);
    const auto *sng_30 = buffer.data(sng + 30);
    const auto *sng_32 = buffer.data(sng + 32);
    const auto *sng_33 = buffer.data(sng + 33);
    const auto *sng_35 = buffer.data(sng + 35);
    const auto *sng_40 = buffer.data(sng + 40);
    const auto *sng_41 = buffer.data(sng + 41);
    const auto *sng_42 = buffer.data(sng + 42);
    const auto *sng_43 = buffer.data(sng + 43);
    const auto *sng_44 = buffer.data(sng + 44);
    const auto *sng_45 = buffer.data(sng + 45);
    const auto *sng_47 = buffer.data(sng + 47);
    const auto *sng_48 = buffer.data(sng + 48);
    const auto *sng_50 = buffer.data(sng + 50);
    const auto *sng_51 = buffer.data(sng + 51);
    const auto *sng_54 = buffer.data(sng + 54);
    const auto *sng_55 = buffer.data(sng + 55);
    const auto *sng_56 = buffer.data(sng + 56);
    const auto *sng_57 = buffer.data(sng + 57);
    const auto *sng_58 = buffer.data(sng + 58);
    const auto *sng_59 = buffer.data(sng + 59);
    const auto *sng_60 = buffer.data(sng + 60);
    const auto *sng_62 = buffer.data(sng + 62);
    const auto *sng_63 = buffer.data(sng + 63);
    const auto *sng_65 = buffer.data(sng + 65);
    const auto *sng_70 = buffer.data(sng + 70);
    const auto *sng_71 = buffer.data(sng + 71);
    const auto *sng_72 = buffer.data(sng + 72);
    const auto *sng_73 = buffer.data(sng + 73);
    const auto *sng_74 = buffer.data(sng + 74);
    const auto *sng_75 = buffer.data(sng + 75);
    const auto *sng_77 = buffer.data(sng + 77);
    const auto *sng_78 = buffer.data(sng + 78);
    const auto *sng_80 = buffer.data(sng + 80);
    const auto *sng_81 = buffer.data(sng + 81);
    const auto *sng_84 = buffer.data(sng + 84);
    const auto *sng_85 = buffer.data(sng + 85);
    const auto *sng_86 = buffer.data(sng + 86);
    const auto *sng_87 = buffer.data(sng + 87);
    const auto *sng_88 = buffer.data(sng + 88);
    const auto *sng_89 = buffer.data(sng + 89);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, smg_0, smg_3, snf0_0, snf0_3, \
                         snf1_0, snf1_3, sng_0, sng_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * smg_0[k]
                 + f_1 * snf0_0[k]
                 - f_2 * snf1_0[k]
                 + f_3 * pc_x[k] * sng_0[k];

        t_1[k] = f_3 * pc_y[k] * sng_0[k];

        t_2[k] = f_3 * pc_z[k] * sng_0[k];

        t_3[k] = f_0 * smg_3[k]
                 + f_4 * snf0_3[k]
                 - f_5 * snf1_3[k]
                 + f_3 * pc_x[k] * sng_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pc_x, pc_y, smg_5, smg_6, snf0_5, snf0_6, snf1_5, \
                         snf1_6, sng_2, sng_5, sng_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * pc_y[k] * sng_2[k];

        t_5[k] = f_0 * smg_5[k]
                 + f_4 * snf0_5[k]
                 - f_5 * snf1_5[k]
                 + f_3 * pc_x[k] * sng_5[k];

        t_6[k] = f_0 * smg_6[k]
                 + f_6 * snf0_6[k]
                 - f_7 * snf1_6[k]
                 + f_3 * pc_x[k] * sng_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, pc_x, pc_y, pc_z, smg_9, smg_10, snf0_9, snf1_9, \
                         sng_3, sng_5, sng_9, sng_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_3 * pc_z[k] * sng_3[k];

        t_8[k] = f_3 * pc_y[k] * sng_5[k];

        t_9[k] = f_0 * smg_9[k]
                 + f_6 * snf0_9[k]
                 - f_7 * snf1_9[k]
                 + f_3 * pc_x[k] * sng_9[k];

        t_10[k] = f_0 * smg_10[k]
                  + f_3 * pc_x[k] * sng_10[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, pc_x, smg_11, smg_12, smg_13, smg_14, sng_11, \
                         sng_12, sng_13, sng_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_0 * smg_11[k]
                  + f_3 * pc_x[k] * sng_11[k];

        t_12[k] = f_0 * smg_12[k]
                  + f_3 * pc_x[k] * sng_12[k];

        t_13[k] = f_0 * smg_13[k]
                  + f_3 * pc_x[k] * sng_13[k];

        t_14[k] = f_0 * smg_14[k]
                  + f_3 * pc_x[k] * sng_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pc_y, pc_z, snf0_6, snf0_8, snf0_9, snf1_6, \
                         snf1_8, snf1_9, sng_10, sng_12, sng_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_1 * snf0_6[k]
                  - f_2 * snf1_6[k]
                  + f_3 * pc_y[k] * sng_10[k];

        t_16[k] = f_3 * pc_z[k] * sng_10[k];

        t_17[k] = f_4 * snf0_8[k]
                  - f_5 * snf1_8[k]
                  + f_3 * pc_y[k] * sng_12[k];

        t_18[k] = f_6 * snf0_9[k]
                  - f_7 * snf1_9[k]
                  + f_3 * pc_y[k] * sng_13[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pb_y, pc_y, pc_z, smh0_0, smg_0, \
                         smh1_0, snf0_9, snf1_9, sng_14, sng_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * pc_y[k] * sng_14[k];

        t_20[k] = f_1 * snf0_9[k]
                  - f_2 * snf1_9[k]
                  + f_3 * pc_z[k] * sng_14[k];

        t_21[k] = pb_y[k] * smh0_0[k]
                  - f_8 * pc_y[k] * smh1_0[k];

        t_22[k] = f_9 * smg_0[k]
                  + f_3 * pc_y[k] * sng_15[k];

        t_23[k] = f_3 * pc_z[k] * sng_15[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pb_y, pc_y, smh0_3, smh0_5, smh0_6, smg_1, \
                         smg_2, smg_3, smh1_3, smh1_5, smh1_6, sng_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pb_y[k] * smh0_3[k]
                  + f_10 * smg_1[k]
                  - f_8 * pc_y[k] * smh1_3[k];

        t_25[k] = f_9 * smg_2[k]
                  + f_3 * pc_y[k] * sng_17[k];

        t_26[k] = pb_y[k] * smh0_5[k]
                  - f_8 * pc_y[k] * smh1_5[k];

        t_27[k] = pb_y[k] * smh0_6[k]
                  + f_11 * smg_3[k]
                  - f_8 * pc_y[k] * smh1_6[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pb_y, pc_x, pc_y, pc_z, smh0_9, smg_5, \
                         smg_25, smh1_9, sng_18, sng_20, sng_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_3 * pc_z[k] * sng_18[k];

        t_29[k] = f_9 * smg_5[k]
                  + f_3 * pc_y[k] * sng_20[k];

        t_30[k] = pb_y[k] * smh0_9[k]
                  - f_8 * pc_y[k] * smh1_9[k];

        t_31[k] = f_12 * smg_25[k]
                  + f_3 * pc_x[k] * sng_25[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pc_x, smg_26, smg_27, smg_28, smg_29, sng_26, \
                         sng_27, sng_28, sng_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_12 * smg_26[k]
                  + f_3 * pc_x[k] * sng_26[k];

        t_33[k] = f_12 * smg_27[k]
                  + f_3 * pc_x[k] * sng_27[k];

        t_34[k] = f_12 * smg_28[k]
                  + f_3 * pc_x[k] * sng_28[k];

        t_35[k] = f_12 * smg_29[k]
                  + f_3 * pc_x[k] * sng_29[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pc_y, pc_z, smg_10, smg_12, snf0_16, snf0_18, \
                         snf1_16, snf1_18, sng_25, sng_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_9 * smg_10[k]
                  + f_1 * snf0_16[k]
                  - f_2 * snf1_16[k]
                  + f_3 * pc_y[k] * sng_25[k];

        t_37[k] = f_3 * pc_z[k] * sng_25[k];

        t_38[k] = f_9 * smg_12[k]
                  + f_4 * snf0_18[k]
                  - f_5 * snf1_18[k]
                  + f_3 * pc_y[k] * sng_27[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pb_y, pc_y, smh0_20, smg_13, smg_14, smh1_20, \
                         snf0_19, snf1_19, sng_28, sng_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_9 * smg_13[k]
                  + f_6 * snf0_19[k]
                  - f_7 * snf1_19[k]
                  + f_3 * pc_y[k] * sng_28[k];

        t_40[k] = f_9 * smg_14[k]
                  + f_3 * pc_y[k] * sng_29[k];

        t_41[k] = pb_y[k] * smh0_20[k]
                  - f_8 * pc_y[k] * smh1_20[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, pb_z, pc_y, pc_z, smh0_0, smh0_3, \
                         smg_0, smh1_0, smh1_3, sng_30, sng_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_z[k] * smh0_0[k]
                  - f_8 * pc_z[k] * smh1_0[k];

        t_43[k] = f_3 * pc_y[k] * sng_30[k];

        t_44[k] = f_9 * smg_0[k]
                  + f_3 * pc_z[k] * sng_30[k];

        t_45[k] = pb_z[k] * smh0_3[k]
                  - f_8 * pc_z[k] * smh1_3[k];

        t_46[k] = f_3 * pc_y[k] * sng_32[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pb_z, pc_y, pc_z, smh0_5, smh0_6, smg_2, \
                         smg_3, smh1_5, smh1_6, sng_33, sng_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = pb_z[k] * smh0_5[k]
                  + f_10 * smg_2[k]
                  - f_8 * pc_z[k] * smh1_5[k];

        t_48[k] = pb_z[k] * smh0_6[k]
                  - f_8 * pc_z[k] * smh1_6[k];

        t_49[k] = f_9 * smg_3[k]
                  + f_3 * pc_z[k] * sng_33[k];

        t_50[k] = f_3 * pc_y[k] * sng_35[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pb_z, pc_x, pc_z, smh0_9, smg_5, smg_40, \
                         smg_41, smg_42, smh1_9, sng_40, sng_41, \
                         sng_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = pb_z[k] * smh0_9[k]
                  + f_11 * smg_5[k]
                  - f_8 * pc_z[k] * smh1_9[k];

        t_52[k] = f_12 * smg_40[k]
                  + f_3 * pc_x[k] * sng_40[k];

        t_53[k] = f_12 * smg_41[k]
                  + f_3 * pc_x[k] * sng_41[k];

        t_54[k] = f_12 * smg_42[k]
                  + f_3 * pc_x[k] * sng_42[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pb_z, pc_x, pc_z, smh0_15, smg_10, smg_43, \
                         smg_44, smh1_15, sng_40, sng_43, sng_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_12 * smg_43[k]
                  + f_3 * pc_x[k] * sng_43[k];

        t_56[k] = f_12 * smg_44[k]
                  + f_3 * pc_x[k] * sng_44[k];

        t_57[k] = pb_z[k] * smh0_15[k]
                  - f_8 * pc_z[k] * smh1_15[k];

        t_58[k] = f_9 * smg_10[k]
                  + f_3 * pc_z[k] * sng_40[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pc_y, pc_z, smg_14, snf0_28, snf0_29, \
                         snf1_28, snf1_29, sng_42, sng_43, sng_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_4 * snf0_28[k]
                  - f_5 * snf1_28[k]
                  + f_3 * pc_y[k] * sng_42[k];

        t_60[k] = f_6 * snf0_29[k]
                  - f_7 * snf1_29[k]
                  + f_3 * pc_y[k] * sng_43[k];

        t_61[k] = f_3 * pc_y[k] * sng_44[k];

        t_62[k] = f_9 * smg_14[k]
                  + f_1 * snf0_29[k]
                  - f_2 * snf1_29[k]
                  + f_3 * pc_z[k] * sng_44[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pc_x, pc_y, pc_z, smg_15, smg_45, smg_48, \
                         snf0_30, snf0_33, snf1_30, snf1_33, sng_45, \
                         sng_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_13 * smg_45[k]
                  + f_1 * snf0_30[k]
                  - f_2 * snf1_30[k]
                  + f_3 * pc_x[k] * sng_45[k];

        t_64[k] = f_10 * smg_15[k]
                  + f_3 * pc_y[k] * sng_45[k];

        t_65[k] = f_3 * pc_z[k] * sng_45[k];

        t_66[k] = f_13 * smg_48[k]
                  + f_4 * snf0_33[k]
                  - f_5 * snf1_33[k]
                  + f_3 * pc_x[k] * sng_48[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, pc_x, pc_y, smg_17, smg_50, smg_51, snf0_35, \
                         snf0_36, snf1_35, snf1_36, sng_47, sng_50, \
                         sng_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_10 * smg_17[k]
                  + f_3 * pc_y[k] * sng_47[k];

        t_68[k] = f_13 * smg_50[k]
                  + f_4 * snf0_35[k]
                  - f_5 * snf1_35[k]
                  + f_3 * pc_x[k] * sng_50[k];

        t_69[k] = f_13 * smg_51[k]
                  + f_6 * snf0_36[k]
                  - f_7 * snf1_36[k]
                  + f_3 * pc_x[k] * sng_51[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pc_x, pc_y, pc_z, smg_20, smg_54, smg_55, \
                         snf0_39, snf1_39, sng_48, sng_50, sng_54, \
                         sng_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_3 * pc_z[k] * sng_48[k];

        t_71[k] = f_10 * smg_20[k]
                  + f_3 * pc_y[k] * sng_50[k];

        t_72[k] = f_13 * smg_54[k]
                  + f_6 * snf0_39[k]
                  - f_7 * snf1_39[k]
                  + f_3 * pc_x[k] * sng_54[k];

        t_73[k] = f_13 * smg_55[k]
                  + f_3 * pc_x[k] * sng_55[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pc_x, smg_56, smg_57, smg_58, smg_59, sng_56, \
                         sng_57, sng_58, sng_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_13 * smg_56[k]
                  + f_3 * pc_x[k] * sng_56[k];

        t_75[k] = f_13 * smg_57[k]
                  + f_3 * pc_x[k] * sng_57[k];

        t_76[k] = f_13 * smg_58[k]
                  + f_3 * pc_x[k] * sng_58[k];

        t_77[k] = f_13 * smg_59[k]
                  + f_3 * pc_x[k] * sng_59[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pc_y, pc_z, smg_25, smg_27, snf0_36, snf0_38, \
                         snf1_36, snf1_38, sng_55, sng_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_10 * smg_25[k]
                  + f_1 * snf0_36[k]
                  - f_2 * snf1_36[k]
                  + f_3 * pc_y[k] * sng_55[k];

        t_79[k] = f_3 * pc_z[k] * sng_55[k];

        t_80[k] = f_10 * smg_27[k]
                  + f_4 * snf0_38[k]
                  - f_5 * snf1_38[k]
                  + f_3 * pc_y[k] * sng_57[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pb_y, pc_y, pc_z, smh0_42, smg_28, smg_29, \
                         smh1_42, snf0_39, snf1_39, sng_58, sng_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_10 * smg_28[k]
                  + f_6 * snf0_39[k]
                  - f_7 * snf1_39[k]
                  + f_3 * pc_y[k] * sng_58[k];

        t_82[k] = f_10 * smg_29[k]
                  + f_3 * pc_y[k] * sng_59[k];

        t_83[k] = f_1 * snf0_39[k]
                  - f_2 * snf1_39[k]
                  + f_3 * pc_z[k] * sng_59[k];

        t_84[k] = pb_y[k] * smh0_42[k]
                  - f_8 * pc_y[k] * smh1_42[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pb_z, pc_y, pc_z, smh0_24, smg_15, smg_30, \
                         smg_32, smh1_24, sng_60, sng_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_9 * smg_30[k]
                  + f_3 * pc_y[k] * sng_60[k];

        t_86[k] = f_9 * smg_15[k]
                  + f_3 * pc_z[k] * sng_60[k];

        t_87[k] = pb_z[k] * smh0_24[k]
                  - f_8 * pc_z[k] * smh1_24[k];

        t_88[k] = f_9 * smg_32[k]
                  + f_3 * pc_y[k] * sng_62[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, pb_y, pb_z, pc_y, pc_z, smh0_27, smh0_47, \
                         smg_18, smg_35, smh1_27, smh1_47, sng_63, \
                         sng_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = pb_y[k] * smh0_47[k]
                  - f_8 * pc_y[k] * smh1_47[k];

        t_90[k] = pb_z[k] * smh0_27[k]
                  - f_8 * pc_z[k] * smh1_27[k];

        t_91[k] = f_9 * smg_18[k]
                  + f_3 * pc_z[k] * sng_63[k];

        t_92[k] = f_9 * smg_35[k]
                  + f_3 * pc_y[k] * sng_65[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, pb_y, pc_x, pc_y, smh0_51, smg_70, smg_71, \
                         smg_72, smh1_51, sng_70, sng_71, sng_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = pb_y[k] * smh0_51[k]
                  - f_8 * pc_y[k] * smh1_51[k];

        t_94[k] = f_13 * smg_70[k]
                  + f_3 * pc_x[k] * sng_70[k];

        t_95[k] = f_13 * smg_71[k]
                  + f_3 * pc_x[k] * sng_71[k];

        t_96[k] = f_13 * smg_72[k]
                  + f_3 * pc_x[k] * sng_72[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pb_z, pc_x, pc_z, smh0_36, smg_25, smg_73, \
                         smg_74, smh1_36, sng_70, sng_73, sng_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_13 * smg_73[k]
                  + f_3 * pc_x[k] * sng_73[k];

        t_98[k] = f_13 * smg_74[k]
                  + f_3 * pc_x[k] * sng_74[k];

        t_99[k] = pb_z[k] * smh0_36[k]
                  - f_8 * pc_z[k] * smh1_36[k];

        t_100[k] = f_9 * smg_25[k]
                   + f_3 * pc_z[k] * sng_70[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, pc_y, smg_42, smg_43, smg_44, snf0_48, snf0_49, \
                         snf1_48, snf1_49, sng_72, sng_73, sng_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_9 * smg_42[k]
                   + f_4 * snf0_48[k]
                   - f_5 * snf1_48[k]
                   + f_3 * pc_y[k] * sng_72[k];

        t_102[k] = f_9 * smg_43[k]
                   + f_6 * snf0_49[k]
                   - f_7 * snf1_49[k]
                   + f_3 * pc_y[k] * sng_73[k];

        t_103[k] = f_9 * smg_44[k]
                   + f_3 * pc_y[k] * sng_74[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pb_y, pc_x, pc_y, pc_z, smh0_62, smg_30, \
                         smg_75, smh1_62, snf0_50, snf1_50, sng_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pb_y[k] * smh0_62[k]
                   - f_8 * pc_y[k] * smh1_62[k];

        t_105[k] = f_13 * smg_75[k]
                   + f_1 * snf0_50[k]
                   - f_2 * snf1_50[k]
                   + f_3 * pc_x[k] * sng_75[k];

        t_106[k] = f_3 * pc_y[k] * sng_75[k];

        t_107[k] = f_10 * smg_30[k]
                   + f_3 * pc_z[k] * sng_75[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, pc_x, pc_y, smg_78, smg_80, snf0_53, snf0_55, \
                         snf1_53, snf1_55, sng_77, sng_78, sng_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_13 * smg_78[k]
                   + f_4 * snf0_53[k]
                   - f_5 * snf1_53[k]
                   + f_3 * pc_x[k] * sng_78[k];

        t_109[k] = f_3 * pc_y[k] * sng_77[k];

        t_110[k] = f_13 * smg_80[k]
                   + f_4 * snf0_55[k]
                   - f_5 * snf1_55[k]
                   + f_3 * pc_x[k] * sng_80[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, pc_x, pc_y, pc_z, smg_33, smg_81, snf0_56, \
                         snf1_56, sng_78, sng_80, sng_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_13 * smg_81[k]
                   + f_6 * snf0_56[k]
                   - f_7 * snf1_56[k]
                   + f_3 * pc_x[k] * sng_81[k];

        t_112[k] = f_10 * smg_33[k]
                   + f_3 * pc_z[k] * sng_78[k];

        t_113[k] = f_3 * pc_y[k] * sng_80[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pc_x, smg_84, smg_85, smg_86, smg_87, \
                         snf0_59, snf1_59, sng_84, sng_85, sng_86, \
                         sng_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_13 * smg_84[k]
                   + f_6 * snf0_59[k]
                   - f_7 * snf1_59[k]
                   + f_3 * pc_x[k] * sng_84[k];

        t_115[k] = f_13 * smg_85[k]
                   + f_3 * pc_x[k] * sng_85[k];

        t_116[k] = f_13 * smg_86[k]
                   + f_3 * pc_x[k] * sng_86[k];

        t_117[k] = f_13 * smg_87[k]
                   + f_3 * pc_x[k] * sng_87[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, pc_x, pc_y, pc_z, smg_40, smg_88, smg_89, \
                         snf0_56, snf1_56, sng_85, sng_88, sng_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_13 * smg_88[k]
                   + f_3 * pc_x[k] * sng_88[k];

        t_119[k] = f_13 * smg_89[k]
                   + f_3 * pc_x[k] * sng_89[k];

        t_120[k] = f_1 * snf0_56[k]
                   - f_2 * snf1_56[k]
                   + f_3 * pc_y[k] * sng_85[k];

        t_121[k] = f_10 * smg_40[k]
                   + f_3 * pc_z[k] * sng_85[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, pc_y, pc_z, smg_44, snf0_58, snf0_59, \
                         snf1_58, snf1_59, sng_87, sng_88, sng_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_4 * snf0_58[k]
                   - f_5 * snf1_58[k]
                   + f_3 * pc_y[k] * sng_87[k];

        t_123[k] = f_6 * snf0_59[k]
                   - f_7 * snf1_59[k]
                   + f_3 * pc_y[k] * sng_88[k];

        t_124[k] = f_3 * pc_y[k] * sng_89[k];

        t_125[k] = f_10 * smg_44[k]
                   + f_1 * snf0_59[k]
                   - f_2 * snf1_59[k]
                   + f_3 * pc_z[k] * sng_89[k];
    }
}

static auto
compute_prim_snh_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t smh0,
                                                          const size_t smg, const size_t smh1,
                                                          const size_t snf0, const size_t snf1,
                                                          const size_t sng, const size_t ncols,
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
    const auto f_14 = 3.5 / q;
    const auto f_15 = 3.0 / q;
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

    const auto *smh0_63 = buffer.data(smh0 + 63);
    const auto *smh0_66 = buffer.data(smh0 + 66);
    const auto *smh0_69 = buffer.data(smh0 + 69);
    const auto *smh0_78 = buffer.data(smh0 + 78);
    const auto *smh0_105 = buffer.data(smh0 + 105);
    const auto *smh0_108 = buffer.data(smh0 + 108);
    const auto *smh0_110 = buffer.data(smh0 + 110);
    const auto *smh0_111 = buffer.data(smh0 + 111);
    const auto *smh0_114 = buffer.data(smh0 + 114);
    const auto *smh0_125 = buffer.data(smh0 + 125);
    const auto *smh0_126 = buffer.data(smh0 + 126);
    const auto *smh0_129 = buffer.data(smh0 + 129);
    const auto *smh0_132 = buffer.data(smh0 + 132);

    const auto *smg_45 = buffer.data(smg + 45);
    const auto *smg_47 = buffer.data(smg + 47);
    const auto *smg_48 = buffer.data(smg + 48);
    const auto *smg_50 = buffer.data(smg + 50);
    const auto *smg_55 = buffer.data(smg + 55);
    const auto *smg_57 = buffer.data(smg + 57);
    const auto *smg_58 = buffer.data(smg + 58);
    const auto *smg_59 = buffer.data(smg + 59);
    const auto *smg_60 = buffer.data(smg + 60);
    const auto *smg_62 = buffer.data(smg + 62);
    const auto *smg_63 = buffer.data(smg + 63);
    const auto *smg_65 = buffer.data(smg + 65);
    const auto *smg_70 = buffer.data(smg + 70);
    const auto *smg_72 = buffer.data(smg + 72);
    const auto *smg_73 = buffer.data(smg + 73);
    const auto *smg_74 = buffer.data(smg + 74);
    const auto *smg_75 = buffer.data(smg + 75);
    const auto *smg_76 = buffer.data(smg + 76);
    const auto *smg_77 = buffer.data(smg + 77);
    const auto *smg_78 = buffer.data(smg + 78);
    const auto *smg_80 = buffer.data(smg + 80);
    const auto *smg_85 = buffer.data(smg + 85);
    const auto *smg_87 = buffer.data(smg + 87);
    const auto *smg_88 = buffer.data(smg + 88);
    const auto *smg_89 = buffer.data(smg + 89);
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
    const auto *smg_110 = buffer.data(smg + 110);
    const auto *smg_114 = buffer.data(smg + 114);
    const auto *smg_115 = buffer.data(smg + 115);
    const auto *smg_116 = buffer.data(smg + 116);
    const auto *smg_117 = buffer.data(smg + 117);
    const auto *smg_118 = buffer.data(smg + 118);
    const auto *smg_119 = buffer.data(smg + 119);
    const auto *smg_130 = buffer.data(smg + 130);
    const auto *smg_131 = buffer.data(smg + 131);
    const auto *smg_132 = buffer.data(smg + 132);
    const auto *smg_133 = buffer.data(smg + 133);
    const auto *smg_134 = buffer.data(smg + 134);
    const auto *smg_135 = buffer.data(smg + 135);
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
    const auto *smg_153 = buffer.data(smg + 153);
    const auto *smg_155 = buffer.data(smg + 155);
    const auto *smg_156 = buffer.data(smg + 156);
    const auto *smg_159 = buffer.data(smg + 159);
    const auto *smg_160 = buffer.data(smg + 160);
    const auto *smg_161 = buffer.data(smg + 161);
    const auto *smg_162 = buffer.data(smg + 162);
    const auto *smg_163 = buffer.data(smg + 163);
    const auto *smg_164 = buffer.data(smg + 164);
    const auto *smg_170 = buffer.data(smg + 170);
    const auto *smg_174 = buffer.data(smg + 174);
    const auto *smg_175 = buffer.data(smg + 175);
    const auto *smg_176 = buffer.data(smg + 176);

    const auto *smh1_63 = buffer.data(smh1 + 63);
    const auto *smh1_66 = buffer.data(smh1 + 66);
    const auto *smh1_69 = buffer.data(smh1 + 69);
    const auto *smh1_78 = buffer.data(smh1 + 78);
    const auto *smh1_105 = buffer.data(smh1 + 105);
    const auto *smh1_108 = buffer.data(smh1 + 108);
    const auto *smh1_110 = buffer.data(smh1 + 110);
    const auto *smh1_111 = buffer.data(smh1 + 111);
    const auto *smh1_114 = buffer.data(smh1 + 114);
    const auto *smh1_125 = buffer.data(smh1 + 125);
    const auto *smh1_126 = buffer.data(smh1 + 126);
    const auto *smh1_129 = buffer.data(smh1 + 129);
    const auto *smh1_132 = buffer.data(smh1 + 132);

    const auto *snf0_60 = buffer.data(snf0 + 60);
    const auto *snf0_63 = buffer.data(snf0 + 63);
    const auto *snf0_65 = buffer.data(snf0 + 65);
    const auto *snf0_66 = buffer.data(snf0 + 66);
    const auto *snf0_68 = buffer.data(snf0 + 68);
    const auto *snf0_69 = buffer.data(snf0 + 69);
    const auto *snf0_75 = buffer.data(snf0 + 75);
    const auto *snf0_78 = buffer.data(snf0 + 78);
    const auto *snf0_79 = buffer.data(snf0 + 79);
    const auto *snf0_86 = buffer.data(snf0 + 86);
    const auto *snf0_88 = buffer.data(snf0 + 88);
    const auto *snf0_89 = buffer.data(snf0 + 89);
    const auto *snf0_90 = buffer.data(snf0 + 90);
    const auto *snf0_93 = buffer.data(snf0 + 93);
    const auto *snf0_95 = buffer.data(snf0 + 95);
    const auto *snf0_96 = buffer.data(snf0 + 96);
    const auto *snf0_98 = buffer.data(snf0 + 98);
    const auto *snf0_99 = buffer.data(snf0 + 99);
    const auto *snf0_100 = buffer.data(snf0 + 100);
    const auto *snf0_103 = buffer.data(snf0 + 103);
    const auto *snf0_105 = buffer.data(snf0 + 105);
    const auto *snf0_106 = buffer.data(snf0 + 106);
    const auto *snf0_108 = buffer.data(snf0 + 108);
    const auto *snf0_109 = buffer.data(snf0 + 109);
    const auto *snf0_115 = buffer.data(snf0 + 115);
    const auto *snf0_119 = buffer.data(snf0 + 119);

    const auto *snf1_60 = buffer.data(snf1 + 60);
    const auto *snf1_63 = buffer.data(snf1 + 63);
    const auto *snf1_65 = buffer.data(snf1 + 65);
    const auto *snf1_66 = buffer.data(snf1 + 66);
    const auto *snf1_68 = buffer.data(snf1 + 68);
    const auto *snf1_69 = buffer.data(snf1 + 69);
    const auto *snf1_75 = buffer.data(snf1 + 75);
    const auto *snf1_78 = buffer.data(snf1 + 78);
    const auto *snf1_79 = buffer.data(snf1 + 79);
    const auto *snf1_86 = buffer.data(snf1 + 86);
    const auto *snf1_88 = buffer.data(snf1 + 88);
    const auto *snf1_89 = buffer.data(snf1 + 89);
    const auto *snf1_90 = buffer.data(snf1 + 90);
    const auto *snf1_93 = buffer.data(snf1 + 93);
    const auto *snf1_95 = buffer.data(snf1 + 95);
    const auto *snf1_96 = buffer.data(snf1 + 96);
    const auto *snf1_98 = buffer.data(snf1 + 98);
    const auto *snf1_99 = buffer.data(snf1 + 99);
    const auto *snf1_100 = buffer.data(snf1 + 100);
    const auto *snf1_103 = buffer.data(snf1 + 103);
    const auto *snf1_105 = buffer.data(snf1 + 105);
    const auto *snf1_106 = buffer.data(snf1 + 106);
    const auto *snf1_108 = buffer.data(snf1 + 108);
    const auto *snf1_109 = buffer.data(snf1 + 109);
    const auto *snf1_115 = buffer.data(snf1 + 115);
    const auto *snf1_119 = buffer.data(snf1 + 119);

    const auto *sng_90 = buffer.data(sng + 90);
    const auto *sng_92 = buffer.data(sng + 92);
    const auto *sng_93 = buffer.data(sng + 93);
    const auto *sng_95 = buffer.data(sng + 95);
    const auto *sng_96 = buffer.data(sng + 96);
    const auto *sng_99 = buffer.data(sng + 99);
    const auto *sng_100 = buffer.data(sng + 100);
    const auto *sng_101 = buffer.data(sng + 101);
    const auto *sng_102 = buffer.data(sng + 102);
    const auto *sng_103 = buffer.data(sng + 103);
    const auto *sng_104 = buffer.data(sng + 104);
    const auto *sng_105 = buffer.data(sng + 105);
    const auto *sng_107 = buffer.data(sng + 107);
    const auto *sng_108 = buffer.data(sng + 108);
    const auto *sng_110 = buffer.data(sng + 110);
    const auto *sng_114 = buffer.data(sng + 114);
    const auto *sng_115 = buffer.data(sng + 115);
    const auto *sng_116 = buffer.data(sng + 116);
    const auto *sng_117 = buffer.data(sng + 117);
    const auto *sng_118 = buffer.data(sng + 118);
    const auto *sng_119 = buffer.data(sng + 119);
    const auto *sng_120 = buffer.data(sng + 120);
    const auto *sng_122 = buffer.data(sng + 122);
    const auto *sng_123 = buffer.data(sng + 123);
    const auto *sng_125 = buffer.data(sng + 125);
    const auto *sng_130 = buffer.data(sng + 130);
    const auto *sng_131 = buffer.data(sng + 131);
    const auto *sng_132 = buffer.data(sng + 132);
    const auto *sng_133 = buffer.data(sng + 133);
    const auto *sng_134 = buffer.data(sng + 134);
    const auto *sng_135 = buffer.data(sng + 135);
    const auto *sng_137 = buffer.data(sng + 137);
    const auto *sng_138 = buffer.data(sng + 138);
    const auto *sng_140 = buffer.data(sng + 140);
    const auto *sng_141 = buffer.data(sng + 141);
    const auto *sng_144 = buffer.data(sng + 144);
    const auto *sng_145 = buffer.data(sng + 145);
    const auto *sng_146 = buffer.data(sng + 146);
    const auto *sng_147 = buffer.data(sng + 147);
    const auto *sng_148 = buffer.data(sng + 148);
    const auto *sng_149 = buffer.data(sng + 149);
    const auto *sng_150 = buffer.data(sng + 150);
    const auto *sng_152 = buffer.data(sng + 152);
    const auto *sng_153 = buffer.data(sng + 153);
    const auto *sng_155 = buffer.data(sng + 155);
    const auto *sng_156 = buffer.data(sng + 156);
    const auto *sng_159 = buffer.data(sng + 159);
    const auto *sng_160 = buffer.data(sng + 160);
    const auto *sng_161 = buffer.data(sng + 161);
    const auto *sng_162 = buffer.data(sng + 162);
    const auto *sng_163 = buffer.data(sng + 163);
    const auto *sng_164 = buffer.data(sng + 164);
    const auto *sng_165 = buffer.data(sng + 165);
    const auto *sng_167 = buffer.data(sng + 167);
    const auto *sng_168 = buffer.data(sng + 168);
    const auto *sng_170 = buffer.data(sng + 170);
    const auto *sng_174 = buffer.data(sng + 174);
    const auto *sng_175 = buffer.data(sng + 175);
    const auto *sng_176 = buffer.data(sng + 176);

#pragma omp simd aligned(t_126, t_127, t_128, t_129, pc_x, pc_y, pc_z, smg_45, smg_90, smg_93, \
                         snf0_60, snf0_63, snf1_60, snf1_63, sng_90, \
                         sng_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_14 * smg_90[k]
                   + f_1 * snf0_60[k]
                   - f_2 * snf1_60[k]
                   + f_3 * pc_x[k] * sng_90[k];

        t_127[k] = f_11 * smg_45[k]
                   + f_3 * pc_y[k] * sng_90[k];

        t_128[k] = f_3 * pc_z[k] * sng_90[k];

        t_129[k] = f_14 * smg_93[k]
                   + f_4 * snf0_63[k]
                   - f_5 * snf1_63[k]
                   + f_3 * pc_x[k] * sng_93[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, pc_x, pc_y, smg_47, smg_95, smg_96, snf0_65, \
                         snf0_66, snf1_65, snf1_66, sng_92, sng_95, \
                         sng_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_11 * smg_47[k]
                   + f_3 * pc_y[k] * sng_92[k];

        t_131[k] = f_14 * smg_95[k]
                   + f_4 * snf0_65[k]
                   - f_5 * snf1_65[k]
                   + f_3 * pc_x[k] * sng_95[k];

        t_132[k] = f_14 * smg_96[k]
                   + f_6 * snf0_66[k]
                   - f_7 * snf1_66[k]
                   + f_3 * pc_x[k] * sng_96[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, pc_x, pc_y, pc_z, smg_50, smg_99, \
                         smg_100, snf0_69, snf1_69, sng_93, sng_95, sng_99, \
                         sng_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_3 * pc_z[k] * sng_93[k];

        t_134[k] = f_11 * smg_50[k]
                   + f_3 * pc_y[k] * sng_95[k];

        t_135[k] = f_14 * smg_99[k]
                   + f_6 * snf0_69[k]
                   - f_7 * snf1_69[k]
                   + f_3 * pc_x[k] * sng_99[k];

        t_136[k] = f_14 * smg_100[k]
                   + f_3 * pc_x[k] * sng_100[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, t_140, pc_x, smg_101, smg_102, smg_103, smg_104, \
                         sng_101, sng_102, sng_103, sng_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_14 * smg_101[k]
                   + f_3 * pc_x[k] * sng_101[k];

        t_138[k] = f_14 * smg_102[k]
                   + f_3 * pc_x[k] * sng_102[k];

        t_139[k] = f_14 * smg_103[k]
                   + f_3 * pc_x[k] * sng_103[k];

        t_140[k] = f_14 * smg_104[k]
                   + f_3 * pc_x[k] * sng_104[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, pc_y, pc_z, smg_55, smg_57, snf0_66, snf0_68, \
                         snf1_66, snf1_68, sng_100, sng_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_11 * smg_55[k]
                   + f_1 * snf0_66[k]
                   - f_2 * snf1_66[k]
                   + f_3 * pc_y[k] * sng_100[k];

        t_142[k] = f_3 * pc_z[k] * sng_100[k];

        t_143[k] = f_11 * smg_57[k]
                   + f_4 * snf0_68[k]
                   - f_5 * snf1_68[k]
                   + f_3 * pc_y[k] * sng_102[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, pb_z, pc_y, pc_z, smh0_63, smg_58, \
                         smg_59, smh1_63, snf0_69, snf1_69, sng_103, \
                         sng_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_11 * smg_58[k]
                   + f_6 * snf0_69[k]
                   - f_7 * snf1_69[k]
                   + f_3 * pc_y[k] * sng_103[k];

        t_145[k] = f_11 * smg_59[k]
                   + f_3 * pc_y[k] * sng_104[k];

        t_146[k] = f_1 * snf0_69[k]
                   - f_2 * snf1_69[k]
                   + f_3 * pc_z[k] * sng_104[k];

        t_147[k] = pb_z[k] * smh0_63[k]
                   - f_8 * pc_z[k] * smh1_63[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, t_151, pb_z, pc_y, pc_z, smh0_66, smg_45, \
                         smg_60, smg_62, smh1_66, sng_105, sng_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_10 * smg_60[k]
                   + f_3 * pc_y[k] * sng_105[k];

        t_149[k] = f_9 * smg_45[k]
                   + f_3 * pc_z[k] * sng_105[k];

        t_150[k] = pb_z[k] * smh0_66[k]
                   - f_8 * pc_z[k] * smh1_66[k];

        t_151[k] = f_10 * smg_62[k]
                   + f_3 * pc_y[k] * sng_107[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, pb_z, pc_x, pc_z, smh0_69, smg_48, smg_110, \
                         smh1_69, snf0_75, snf1_75, sng_108, sng_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_14 * smg_110[k]
                   + f_4 * snf0_75[k]
                   - f_5 * snf1_75[k]
                   + f_3 * pc_x[k] * sng_110[k];

        t_153[k] = pb_z[k] * smh0_69[k]
                   - f_8 * pc_z[k] * smh1_69[k];

        t_154[k] = f_9 * smg_48[k]
                   + f_3 * pc_z[k] * sng_108[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, pc_x, pc_y, smg_65, smg_114, smg_115, \
                         smg_116, snf0_79, snf1_79, sng_110, sng_114, sng_115, \
                         sng_116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_10 * smg_65[k]
                   + f_3 * pc_y[k] * sng_110[k];

        t_156[k] = f_14 * smg_114[k]
                   + f_6 * snf0_79[k]
                   - f_7 * snf1_79[k]
                   + f_3 * pc_x[k] * sng_114[k];

        t_157[k] = f_14 * smg_115[k]
                   + f_3 * pc_x[k] * sng_115[k];

        t_158[k] = f_14 * smg_116[k]
                   + f_3 * pc_x[k] * sng_116[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, pb_z, pc_x, pc_z, smh0_78, smg_117, \
                         smg_118, smg_119, smh1_78, sng_117, sng_118, \
                         sng_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_14 * smg_117[k]
                   + f_3 * pc_x[k] * sng_117[k];

        t_160[k] = f_14 * smg_118[k]
                   + f_3 * pc_x[k] * sng_118[k];

        t_161[k] = f_14 * smg_119[k]
                   + f_3 * pc_x[k] * sng_119[k];

        t_162[k] = pb_z[k] * smh0_78[k]
                   - f_8 * pc_z[k] * smh1_78[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, pc_y, pc_z, smg_55, smg_72, smg_73, snf0_78, \
                         snf0_79, snf1_78, snf1_79, sng_115, sng_117, \
                         sng_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_9 * smg_55[k]
                   + f_3 * pc_z[k] * sng_115[k];

        t_164[k] = f_10 * smg_72[k]
                   + f_4 * snf0_78[k]
                   - f_5 * snf1_78[k]
                   + f_3 * pc_y[k] * sng_117[k];

        t_165[k] = f_10 * smg_73[k]
                   + f_6 * snf0_79[k]
                   - f_7 * snf1_79[k]
                   + f_3 * pc_y[k] * sng_118[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, pb_y, pc_y, pc_z, smh0_105, smg_59, \
                         smg_74, smg_75, smh1_105, snf0_79, snf1_79, sng_119, \
                         sng_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_10 * smg_74[k]
                   + f_3 * pc_y[k] * sng_119[k];

        t_167[k] = f_9 * smg_59[k]
                   + f_1 * snf0_79[k]
                   - f_2 * snf1_79[k]
                   + f_3 * pc_z[k] * sng_119[k];

        t_168[k] = pb_y[k] * smh0_105[k]
                   - f_8 * pc_y[k] * smh1_105[k];

        t_169[k] = f_9 * smg_75[k]
                   + f_3 * pc_y[k] * sng_120[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, pb_y, pc_y, pc_z, smh0_108, smh0_110, \
                         smg_60, smg_76, smg_77, smh1_108, smh1_110, sng_120, \
                         sng_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_10 * smg_60[k]
                   + f_3 * pc_z[k] * sng_120[k];

        t_171[k] = pb_y[k] * smh0_108[k]
                   + f_10 * smg_76[k]
                   - f_8 * pc_y[k] * smh1_108[k];

        t_172[k] = f_9 * smg_77[k]
                   + f_3 * pc_y[k] * sng_122[k];

        t_173[k] = pb_y[k] * smh0_110[k]
                   - f_8 * pc_y[k] * smh1_110[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, pb_y, pc_y, pc_z, smh0_111, smh0_114, \
                         smg_63, smg_78, smg_80, smh1_111, smh1_114, sng_123, \
                         sng_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = pb_y[k] * smh0_111[k]
                   + f_11 * smg_78[k]
                   - f_8 * pc_y[k] * smh1_111[k];

        t_175[k] = f_10 * smg_63[k]
                   + f_3 * pc_z[k] * sng_123[k];

        t_176[k] = f_9 * smg_80[k]
                   + f_3 * pc_y[k] * sng_125[k];

        t_177[k] = pb_y[k] * smh0_114[k]
                   - f_8 * pc_y[k] * smh1_114[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, t_181, t_182, pc_x, smg_130, smg_131, smg_132, \
                         smg_133, smg_134, sng_130, sng_131, sng_132, sng_133, \
                         sng_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_14 * smg_130[k]
                   + f_3 * pc_x[k] * sng_130[k];

        t_179[k] = f_14 * smg_131[k]
                   + f_3 * pc_x[k] * sng_131[k];

        t_180[k] = f_14 * smg_132[k]
                   + f_3 * pc_x[k] * sng_132[k];

        t_181[k] = f_14 * smg_133[k]
                   + f_3 * pc_x[k] * sng_133[k];

        t_182[k] = f_14 * smg_134[k]
                   + f_3 * pc_x[k] * sng_134[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, pc_y, pc_z, smg_70, smg_85, smg_87, snf0_86, \
                         snf0_88, snf1_86, snf1_88, sng_130, sng_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_9 * smg_85[k]
                   + f_1 * snf0_86[k]
                   - f_2 * snf1_86[k]
                   + f_3 * pc_y[k] * sng_130[k];

        t_184[k] = f_10 * smg_70[k]
                   + f_3 * pc_z[k] * sng_130[k];

        t_185[k] = f_9 * smg_87[k]
                   + f_4 * snf0_88[k]
                   - f_5 * snf1_88[k]
                   + f_3 * pc_y[k] * sng_132[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, pb_y, pc_y, smh0_125, smg_88, smg_89, smh1_125, \
                         snf0_89, snf1_89, sng_133, sng_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_9 * smg_88[k]
                   + f_6 * snf0_89[k]
                   - f_7 * snf1_89[k]
                   + f_3 * pc_y[k] * sng_133[k];

        t_187[k] = f_9 * smg_89[k]
                   + f_3 * pc_y[k] * sng_134[k];

        t_188[k] = pb_y[k] * smh0_125[k]
                   - f_8 * pc_y[k] * smh1_125[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, pc_x, pc_y, pc_z, smg_75, smg_135, \
                         smg_138, snf0_90, snf0_93, snf1_90, snf1_93, sng_135, \
                         sng_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_14 * smg_135[k]
                   + f_1 * snf0_90[k]
                   - f_2 * snf1_90[k]
                   + f_3 * pc_x[k] * sng_135[k];

        t_190[k] = f_3 * pc_y[k] * sng_135[k];

        t_191[k] = f_11 * smg_75[k]
                   + f_3 * pc_z[k] * sng_135[k];

        t_192[k] = f_14 * smg_138[k]
                   + f_4 * snf0_93[k]
                   - f_5 * snf1_93[k]
                   + f_3 * pc_x[k] * sng_138[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, pc_x, pc_y, smg_140, smg_141, snf0_95, snf0_96, \
                         snf1_95, snf1_96, sng_137, sng_140, sng_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_3 * pc_y[k] * sng_137[k];

        t_194[k] = f_14 * smg_140[k]
                   + f_4 * snf0_95[k]
                   - f_5 * snf1_95[k]
                   + f_3 * pc_x[k] * sng_140[k];

        t_195[k] = f_14 * smg_141[k]
                   + f_6 * snf0_96[k]
                   - f_7 * snf1_96[k]
                   + f_3 * pc_x[k] * sng_141[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, pc_x, pc_y, pc_z, smg_78, smg_144, \
                         smg_145, snf0_99, snf1_99, sng_138, sng_140, sng_144, \
                         sng_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = f_11 * smg_78[k]
                   + f_3 * pc_z[k] * sng_138[k];

        t_197[k] = f_3 * pc_y[k] * sng_140[k];

        t_198[k] = f_14 * smg_144[k]
                   + f_6 * snf0_99[k]
                   - f_7 * snf1_99[k]
                   + f_3 * pc_x[k] * sng_144[k];

        t_199[k] = f_14 * smg_145[k]
                   + f_3 * pc_x[k] * sng_145[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, pc_x, smg_146, smg_147, smg_148, smg_149, \
                         sng_146, sng_147, sng_148, sng_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_14 * smg_146[k]
                   + f_3 * pc_x[k] * sng_146[k];

        t_201[k] = f_14 * smg_147[k]
                   + f_3 * pc_x[k] * sng_147[k];

        t_202[k] = f_14 * smg_148[k]
                   + f_3 * pc_x[k] * sng_148[k];

        t_203[k] = f_14 * smg_149[k]
                   + f_3 * pc_x[k] * sng_149[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, pc_y, pc_z, smg_85, snf0_96, snf0_98, \
                         snf0_99, snf1_96, snf1_98, snf1_99, sng_145, sng_147, \
                         sng_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_1 * snf0_96[k]
                   - f_2 * snf1_96[k]
                   + f_3 * pc_y[k] * sng_145[k];

        t_205[k] = f_11 * smg_85[k]
                   + f_3 * pc_z[k] * sng_145[k];

        t_206[k] = f_4 * snf0_98[k]
                   - f_5 * snf1_98[k]
                   + f_3 * pc_y[k] * sng_147[k];

        t_207[k] = f_6 * snf0_99[k]
                   - f_7 * snf1_99[k]
                   + f_3 * pc_y[k] * sng_148[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, pc_x, pc_y, pc_z, smg_89, smg_90, \
                         smg_150, snf0_99, snf0_100, snf1_99, snf1_100, sng_149, \
                         sng_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_3 * pc_y[k] * sng_149[k];

        t_209[k] = f_11 * smg_89[k]
                   + f_1 * snf0_99[k]
                   - f_2 * snf1_99[k]
                   + f_3 * pc_z[k] * sng_149[k];

        t_210[k] = f_15 * smg_150[k]
                   + f_1 * snf0_100[k]
                   - f_2 * snf1_100[k]
                   + f_3 * pc_x[k] * sng_150[k];

        t_211[k] = f_16 * smg_90[k]
                   + f_3 * pc_y[k] * sng_150[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, pc_x, pc_y, pc_z, smg_92, smg_153, snf0_103, \
                         snf1_103, sng_150, sng_152, sng_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_3 * pc_z[k] * sng_150[k];

        t_213[k] = f_15 * smg_153[k]
                   + f_4 * snf0_103[k]
                   - f_5 * snf1_103[k]
                   + f_3 * pc_x[k] * sng_153[k];

        t_214[k] = f_16 * smg_92[k]
                   + f_3 * pc_y[k] * sng_152[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, pc_x, pc_z, smg_155, smg_156, snf0_105, \
                         snf0_106, snf1_105, snf1_106, sng_153, sng_155, \
                         sng_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = f_15 * smg_155[k]
                   + f_4 * snf0_105[k]
                   - f_5 * snf1_105[k]
                   + f_3 * pc_x[k] * sng_155[k];

        t_216[k] = f_15 * smg_156[k]
                   + f_6 * snf0_106[k]
                   - f_7 * snf1_106[k]
                   + f_3 * pc_x[k] * sng_156[k];

        t_217[k] = f_3 * pc_z[k] * sng_153[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, t_221, pc_x, pc_y, smg_95, smg_159, smg_160, \
                         smg_161, snf0_109, snf1_109, sng_155, sng_159, sng_160, \
                         sng_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_16 * smg_95[k]
                   + f_3 * pc_y[k] * sng_155[k];

        t_219[k] = f_15 * smg_159[k]
                   + f_6 * snf0_109[k]
                   - f_7 * snf1_109[k]
                   + f_3 * pc_x[k] * sng_159[k];

        t_220[k] = f_15 * smg_160[k]
                   + f_3 * pc_x[k] * sng_160[k];

        t_221[k] = f_15 * smg_161[k]
                   + f_3 * pc_x[k] * sng_161[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, pc_x, pc_y, smg_100, smg_162, smg_163, \
                         smg_164, snf0_106, snf1_106, sng_160, sng_162, sng_163, \
                         sng_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_15 * smg_162[k]
                   + f_3 * pc_x[k] * sng_162[k];

        t_223[k] = f_15 * smg_163[k]
                   + f_3 * pc_x[k] * sng_163[k];

        t_224[k] = f_15 * smg_164[k]
                   + f_3 * pc_x[k] * sng_164[k];

        t_225[k] = f_16 * smg_100[k]
                   + f_1 * snf0_106[k]
                   - f_2 * snf1_106[k]
                   + f_3 * pc_y[k] * sng_160[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, pc_y, pc_z, smg_102, smg_103, snf0_108, \
                         snf0_109, snf1_108, snf1_109, sng_160, sng_162, \
                         sng_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = f_3 * pc_z[k] * sng_160[k];

        t_227[k] = f_16 * smg_102[k]
                   + f_4 * snf0_108[k]
                   - f_5 * snf1_108[k]
                   + f_3 * pc_y[k] * sng_162[k];

        t_228[k] = f_16 * smg_103[k]
                   + f_6 * snf0_109[k]
                   - f_7 * snf1_109[k]
                   + f_3 * pc_y[k] * sng_163[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, pb_z, pc_y, pc_z, smh0_126, smg_104, \
                         smg_105, smh1_126, snf0_109, snf1_109, sng_164, \
                         sng_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = f_16 * smg_104[k]
                   + f_3 * pc_y[k] * sng_164[k];

        t_230[k] = f_1 * snf0_109[k]
                   - f_2 * snf1_109[k]
                   + f_3 * pc_z[k] * sng_164[k];

        t_231[k] = pb_z[k] * smh0_126[k]
                   - f_8 * pc_z[k] * smh1_126[k];

        t_232[k] = f_11 * smg_105[k]
                   + f_3 * pc_y[k] * sng_165[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, pb_z, pc_y, pc_z, smh0_129, smg_90, smg_107, \
                         smh1_129, sng_165, sng_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = f_9 * smg_90[k]
                   + f_3 * pc_z[k] * sng_165[k];

        t_234[k] = pb_z[k] * smh0_129[k]
                   - f_8 * pc_z[k] * smh1_129[k];

        t_235[k] = f_11 * smg_107[k]
                   + f_3 * pc_y[k] * sng_167[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, pb_z, pc_x, pc_z, smh0_132, smg_93, smg_170, \
                         smh1_132, snf0_115, snf1_115, sng_168, \
                         sng_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = f_15 * smg_170[k]
                   + f_4 * snf0_115[k]
                   - f_5 * snf1_115[k]
                   + f_3 * pc_x[k] * sng_170[k];

        t_237[k] = pb_z[k] * smh0_132[k]
                   - f_8 * pc_z[k] * smh1_132[k];

        t_238[k] = f_9 * smg_93[k]
                   + f_3 * pc_z[k] * sng_168[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, pc_x, pc_y, smg_110, smg_174, smg_175, \
                         smg_176, snf0_119, snf1_119, sng_170, sng_174, sng_175, \
                         sng_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_11 * smg_110[k]
                   + f_3 * pc_y[k] * sng_170[k];

        t_240[k] = f_15 * smg_174[k]
                   + f_6 * snf0_119[k]
                   - f_7 * snf1_119[k]
                   + f_3 * pc_x[k] * sng_174[k];

        t_241[k] = f_15 * smg_175[k]
                   + f_3 * pc_x[k] * sng_175[k];

        t_242[k] = f_15 * smg_176[k]
                   + f_3 * pc_x[k] * sng_176[k];
    }
}

static auto
compute_prim_snh_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t smh0,
                                                          const size_t smg, const size_t smh1,
                                                          const size_t snf0, const size_t snf1,
                                                          const size_t sng, const size_t ncols,
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
    const auto f_15 = 3.0 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;

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

    const auto *smh0_141 = buffer.data(smh0 + 141);
    const auto *smh0_189 = buffer.data(smh0 + 189);
    const auto *smh0_192 = buffer.data(smh0 + 192);
    const auto *smh0_194 = buffer.data(smh0 + 194);
    const auto *smh0_195 = buffer.data(smh0 + 195);
    const auto *smh0_198 = buffer.data(smh0 + 198);
    const auto *smh0_209 = buffer.data(smh0 + 209);
    const auto *smh0_210 = buffer.data(smh0 + 210);
    const auto *smh0_213 = buffer.data(smh0 + 213);
    const auto *smh0_216 = buffer.data(smh0 + 216);
    const auto *smh0_225 = buffer.data(smh0 + 225);

    const auto *smg_100 = buffer.data(smg + 100);
    const auto *smg_104 = buffer.data(smg + 104);
    const auto *smg_105 = buffer.data(smg + 105);
    const auto *smg_108 = buffer.data(smg + 108);
    const auto *smg_115 = buffer.data(smg + 115);
    const auto *smg_117 = buffer.data(smg + 117);
    const auto *smg_118 = buffer.data(smg + 118);
    const auto *smg_119 = buffer.data(smg + 119);
    const auto *smg_120 = buffer.data(smg + 120);
    const auto *smg_122 = buffer.data(smg + 122);
    const auto *smg_123 = buffer.data(smg + 123);
    const auto *smg_125 = buffer.data(smg + 125);
    const auto *smg_130 = buffer.data(smg + 130);
    const auto *smg_132 = buffer.data(smg + 132);
    const auto *smg_133 = buffer.data(smg + 133);
    const auto *smg_134 = buffer.data(smg + 134);
    const auto *smg_135 = buffer.data(smg + 135);
    const auto *smg_136 = buffer.data(smg + 136);
    const auto *smg_137 = buffer.data(smg + 137);
    const auto *smg_138 = buffer.data(smg + 138);
    const auto *smg_140 = buffer.data(smg + 140);
    const auto *smg_145 = buffer.data(smg + 145);
    const auto *smg_147 = buffer.data(smg + 147);
    const auto *smg_148 = buffer.data(smg + 148);
    const auto *smg_149 = buffer.data(smg + 149);
    const auto *smg_150 = buffer.data(smg + 150);
    const auto *smg_152 = buffer.data(smg + 152);
    const auto *smg_153 = buffer.data(smg + 153);
    const auto *smg_155 = buffer.data(smg + 155);
    const auto *smg_160 = buffer.data(smg + 160);
    const auto *smg_162 = buffer.data(smg + 162);
    const auto *smg_163 = buffer.data(smg + 163);
    const auto *smg_164 = buffer.data(smg + 164);
    const auto *smg_165 = buffer.data(smg + 165);
    const auto *smg_167 = buffer.data(smg + 167);
    const auto *smg_170 = buffer.data(smg + 170);
    const auto *smg_177 = buffer.data(smg + 177);
    const auto *smg_178 = buffer.data(smg + 178);
    const auto *smg_179 = buffer.data(smg + 179);
    const auto *smg_180 = buffer.data(smg + 180);
    const auto *smg_183 = buffer.data(smg + 183);
    const auto *smg_185 = buffer.data(smg + 185);
    const auto *smg_186 = buffer.data(smg + 186);
    const auto *smg_189 = buffer.data(smg + 189);
    const auto *smg_190 = buffer.data(smg + 190);
    const auto *smg_191 = buffer.data(smg + 191);
    const auto *smg_192 = buffer.data(smg + 192);
    const auto *smg_193 = buffer.data(smg + 193);
    const auto *smg_194 = buffer.data(smg + 194);
    const auto *smg_205 = buffer.data(smg + 205);
    const auto *smg_206 = buffer.data(smg + 206);
    const auto *smg_207 = buffer.data(smg + 207);
    const auto *smg_208 = buffer.data(smg + 208);
    const auto *smg_209 = buffer.data(smg + 209);
    const auto *smg_210 = buffer.data(smg + 210);
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
    const auto *smg_228 = buffer.data(smg + 228);
    const auto *smg_230 = buffer.data(smg + 230);
    const auto *smg_231 = buffer.data(smg + 231);
    const auto *smg_234 = buffer.data(smg + 234);
    const auto *smg_235 = buffer.data(smg + 235);
    const auto *smg_236 = buffer.data(smg + 236);
    const auto *smg_237 = buffer.data(smg + 237);
    const auto *smg_238 = buffer.data(smg + 238);
    const auto *smg_239 = buffer.data(smg + 239);
    const auto *smg_245 = buffer.data(smg + 245);
    const auto *smg_249 = buffer.data(smg + 249);
    const auto *smg_250 = buffer.data(smg + 250);
    const auto *smg_251 = buffer.data(smg + 251);
    const auto *smg_252 = buffer.data(smg + 252);
    const auto *smg_253 = buffer.data(smg + 253);
    const auto *smg_254 = buffer.data(smg + 254);
    const auto *smg_255 = buffer.data(smg + 255);

    const auto *smh1_141 = buffer.data(smh1 + 141);
    const auto *smh1_189 = buffer.data(smh1 + 189);
    const auto *smh1_192 = buffer.data(smh1 + 192);
    const auto *smh1_194 = buffer.data(smh1 + 194);
    const auto *smh1_195 = buffer.data(smh1 + 195);
    const auto *smh1_198 = buffer.data(smh1 + 198);
    const auto *smh1_209 = buffer.data(smh1 + 209);
    const auto *smh1_210 = buffer.data(smh1 + 210);
    const auto *smh1_213 = buffer.data(smh1 + 213);
    const auto *smh1_216 = buffer.data(smh1 + 216);
    const auto *smh1_225 = buffer.data(smh1 + 225);

    const auto *snf0_118 = buffer.data(snf0 + 118);
    const auto *snf0_119 = buffer.data(snf0 + 119);
    const auto *snf0_120 = buffer.data(snf0 + 120);
    const auto *snf0_123 = buffer.data(snf0 + 123);
    const auto *snf0_125 = buffer.data(snf0 + 125);
    const auto *snf0_126 = buffer.data(snf0 + 126);
    const auto *snf0_128 = buffer.data(snf0 + 128);
    const auto *snf0_129 = buffer.data(snf0 + 129);
    const auto *snf0_136 = buffer.data(snf0 + 136);
    const auto *snf0_138 = buffer.data(snf0 + 138);
    const auto *snf0_139 = buffer.data(snf0 + 139);
    const auto *snf0_140 = buffer.data(snf0 + 140);
    const auto *snf0_143 = buffer.data(snf0 + 143);
    const auto *snf0_145 = buffer.data(snf0 + 145);
    const auto *snf0_146 = buffer.data(snf0 + 146);
    const auto *snf0_148 = buffer.data(snf0 + 148);
    const auto *snf0_149 = buffer.data(snf0 + 149);
    const auto *snf0_150 = buffer.data(snf0 + 150);
    const auto *snf0_153 = buffer.data(snf0 + 153);
    const auto *snf0_155 = buffer.data(snf0 + 155);
    const auto *snf0_156 = buffer.data(snf0 + 156);
    const auto *snf0_158 = buffer.data(snf0 + 158);
    const auto *snf0_159 = buffer.data(snf0 + 159);
    const auto *snf0_165 = buffer.data(snf0 + 165);
    const auto *snf0_168 = buffer.data(snf0 + 168);
    const auto *snf0_169 = buffer.data(snf0 + 169);
    const auto *snf0_170 = buffer.data(snf0 + 170);

    const auto *snf1_118 = buffer.data(snf1 + 118);
    const auto *snf1_119 = buffer.data(snf1 + 119);
    const auto *snf1_120 = buffer.data(snf1 + 120);
    const auto *snf1_123 = buffer.data(snf1 + 123);
    const auto *snf1_125 = buffer.data(snf1 + 125);
    const auto *snf1_126 = buffer.data(snf1 + 126);
    const auto *snf1_128 = buffer.data(snf1 + 128);
    const auto *snf1_129 = buffer.data(snf1 + 129);
    const auto *snf1_136 = buffer.data(snf1 + 136);
    const auto *snf1_138 = buffer.data(snf1 + 138);
    const auto *snf1_139 = buffer.data(snf1 + 139);
    const auto *snf1_140 = buffer.data(snf1 + 140);
    const auto *snf1_143 = buffer.data(snf1 + 143);
    const auto *snf1_145 = buffer.data(snf1 + 145);
    const auto *snf1_146 = buffer.data(snf1 + 146);
    const auto *snf1_148 = buffer.data(snf1 + 148);
    const auto *snf1_149 = buffer.data(snf1 + 149);
    const auto *snf1_150 = buffer.data(snf1 + 150);
    const auto *snf1_153 = buffer.data(snf1 + 153);
    const auto *snf1_155 = buffer.data(snf1 + 155);
    const auto *snf1_156 = buffer.data(snf1 + 156);
    const auto *snf1_158 = buffer.data(snf1 + 158);
    const auto *snf1_159 = buffer.data(snf1 + 159);
    const auto *snf1_165 = buffer.data(snf1 + 165);
    const auto *snf1_168 = buffer.data(snf1 + 168);
    const auto *snf1_169 = buffer.data(snf1 + 169);
    const auto *snf1_170 = buffer.data(snf1 + 170);

    const auto *sng_175 = buffer.data(sng + 175);
    const auto *sng_177 = buffer.data(sng + 177);
    const auto *sng_178 = buffer.data(sng + 178);
    const auto *sng_179 = buffer.data(sng + 179);
    const auto *sng_180 = buffer.data(sng + 180);
    const auto *sng_182 = buffer.data(sng + 182);
    const auto *sng_183 = buffer.data(sng + 183);
    const auto *sng_185 = buffer.data(sng + 185);
    const auto *sng_186 = buffer.data(sng + 186);
    const auto *sng_189 = buffer.data(sng + 189);
    const auto *sng_190 = buffer.data(sng + 190);
    const auto *sng_191 = buffer.data(sng + 191);
    const auto *sng_192 = buffer.data(sng + 192);
    const auto *sng_193 = buffer.data(sng + 193);
    const auto *sng_194 = buffer.data(sng + 194);
    const auto *sng_195 = buffer.data(sng + 195);
    const auto *sng_197 = buffer.data(sng + 197);
    const auto *sng_198 = buffer.data(sng + 198);
    const auto *sng_200 = buffer.data(sng + 200);
    const auto *sng_205 = buffer.data(sng + 205);
    const auto *sng_206 = buffer.data(sng + 206);
    const auto *sng_207 = buffer.data(sng + 207);
    const auto *sng_208 = buffer.data(sng + 208);
    const auto *sng_209 = buffer.data(sng + 209);
    const auto *sng_210 = buffer.data(sng + 210);
    const auto *sng_212 = buffer.data(sng + 212);
    const auto *sng_213 = buffer.data(sng + 213);
    const auto *sng_215 = buffer.data(sng + 215);
    const auto *sng_216 = buffer.data(sng + 216);
    const auto *sng_219 = buffer.data(sng + 219);
    const auto *sng_220 = buffer.data(sng + 220);
    const auto *sng_221 = buffer.data(sng + 221);
    const auto *sng_222 = buffer.data(sng + 222);
    const auto *sng_223 = buffer.data(sng + 223);
    const auto *sng_224 = buffer.data(sng + 224);
    const auto *sng_225 = buffer.data(sng + 225);
    const auto *sng_227 = buffer.data(sng + 227);
    const auto *sng_228 = buffer.data(sng + 228);
    const auto *sng_230 = buffer.data(sng + 230);
    const auto *sng_231 = buffer.data(sng + 231);
    const auto *sng_234 = buffer.data(sng + 234);
    const auto *sng_235 = buffer.data(sng + 235);
    const auto *sng_236 = buffer.data(sng + 236);
    const auto *sng_237 = buffer.data(sng + 237);
    const auto *sng_238 = buffer.data(sng + 238);
    const auto *sng_239 = buffer.data(sng + 239);
    const auto *sng_240 = buffer.data(sng + 240);
    const auto *sng_242 = buffer.data(sng + 242);
    const auto *sng_243 = buffer.data(sng + 243);
    const auto *sng_245 = buffer.data(sng + 245);
    const auto *sng_249 = buffer.data(sng + 249);
    const auto *sng_250 = buffer.data(sng + 250);
    const auto *sng_251 = buffer.data(sng + 251);
    const auto *sng_252 = buffer.data(sng + 252);
    const auto *sng_253 = buffer.data(sng + 253);
    const auto *sng_254 = buffer.data(sng + 254);
    const auto *sng_255 = buffer.data(sng + 255);

#pragma omp simd aligned(t_243, t_244, t_245, t_246, pb_z, pc_x, pc_z, smh0_141, smg_177, \
                         smg_178, smg_179, smh1_141, sng_177, sng_178, \
                         sng_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = f_15 * smg_177[k]
                   + f_3 * pc_x[k] * sng_177[k];

        t_244[k] = f_15 * smg_178[k]
                   + f_3 * pc_x[k] * sng_178[k];

        t_245[k] = f_15 * smg_179[k]
                   + f_3 * pc_x[k] * sng_179[k];

        t_246[k] = pb_z[k] * smh0_141[k]
                   - f_8 * pc_z[k] * smh1_141[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, pc_y, pc_z, smg_100, smg_117, smg_118, snf0_118, \
                         snf0_119, snf1_118, snf1_119, sng_175, sng_177, \
                         sng_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = f_9 * smg_100[k]
                   + f_3 * pc_z[k] * sng_175[k];

        t_248[k] = f_11 * smg_117[k]
                   + f_4 * snf0_118[k]
                   - f_5 * snf1_118[k]
                   + f_3 * pc_y[k] * sng_177[k];

        t_249[k] = f_11 * smg_118[k]
                   + f_6 * snf0_119[k]
                   - f_7 * snf1_119[k]
                   + f_3 * pc_y[k] * sng_178[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, pc_x, pc_y, pc_z, smg_104, smg_119, smg_180, \
                         snf0_119, snf0_120, snf1_119, snf1_120, sng_179, \
                         sng_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = f_11 * smg_119[k]
                   + f_3 * pc_y[k] * sng_179[k];

        t_251[k] = f_9 * smg_104[k]
                   + f_1 * snf0_119[k]
                   - f_2 * snf1_119[k]
                   + f_3 * pc_z[k] * sng_179[k];

        t_252[k] = f_15 * smg_180[k]
                   + f_1 * snf0_120[k]
                   - f_2 * snf1_120[k]
                   + f_3 * pc_x[k] * sng_180[k];
    }

#pragma omp simd aligned(t_253, t_254, t_255, t_256, pc_x, pc_y, pc_z, smg_105, smg_120, \
                         smg_122, smg_183, snf0_123, snf1_123, sng_180, sng_182, \
                         sng_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_253[k] = f_10 * smg_120[k]
                   + f_3 * pc_y[k] * sng_180[k];

        t_254[k] = f_10 * smg_105[k]
                   + f_3 * pc_z[k] * sng_180[k];

        t_255[k] = f_15 * smg_183[k]
                   + f_4 * snf0_123[k]
                   - f_5 * snf1_123[k]
                   + f_3 * pc_x[k] * sng_183[k];

        t_256[k] = f_10 * smg_122[k]
                   + f_3 * pc_y[k] * sng_182[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, pc_x, pc_z, smg_108, smg_185, smg_186, snf0_125, \
                         snf0_126, snf1_125, snf1_126, sng_183, sng_185, \
                         sng_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_15 * smg_185[k]
                   + f_4 * snf0_125[k]
                   - f_5 * snf1_125[k]
                   + f_3 * pc_x[k] * sng_185[k];

        t_258[k] = f_15 * smg_186[k]
                   + f_6 * snf0_126[k]
                   - f_7 * snf1_126[k]
                   + f_3 * pc_x[k] * sng_186[k];

        t_259[k] = f_10 * smg_108[k]
                   + f_3 * pc_z[k] * sng_183[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pc_x, pc_y, smg_125, smg_189, smg_190, \
                         smg_191, snf0_129, snf1_129, sng_185, sng_189, sng_190, \
                         sng_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_10 * smg_125[k]
                   + f_3 * pc_y[k] * sng_185[k];

        t_261[k] = f_15 * smg_189[k]
                   + f_6 * snf0_129[k]
                   - f_7 * snf1_129[k]
                   + f_3 * pc_x[k] * sng_189[k];

        t_262[k] = f_15 * smg_190[k]
                   + f_3 * pc_x[k] * sng_190[k];

        t_263[k] = f_15 * smg_191[k]
                   + f_3 * pc_x[k] * sng_191[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, pc_x, pc_y, smg_130, smg_192, smg_193, \
                         smg_194, snf0_126, snf1_126, sng_190, sng_192, sng_193, \
                         sng_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_15 * smg_192[k]
                   + f_3 * pc_x[k] * sng_192[k];

        t_265[k] = f_15 * smg_193[k]
                   + f_3 * pc_x[k] * sng_193[k];

        t_266[k] = f_15 * smg_194[k]
                   + f_3 * pc_x[k] * sng_194[k];

        t_267[k] = f_10 * smg_130[k]
                   + f_1 * snf0_126[k]
                   - f_2 * snf1_126[k]
                   + f_3 * pc_y[k] * sng_190[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, pc_y, pc_z, smg_115, smg_132, smg_133, snf0_128, \
                         snf0_129, snf1_128, snf1_129, sng_190, sng_192, \
                         sng_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_10 * smg_115[k]
                   + f_3 * pc_z[k] * sng_190[k];

        t_269[k] = f_10 * smg_132[k]
                   + f_4 * snf0_128[k]
                   - f_5 * snf1_128[k]
                   + f_3 * pc_y[k] * sng_192[k];

        t_270[k] = f_10 * smg_133[k]
                   + f_6 * snf0_129[k]
                   - f_7 * snf1_129[k]
                   + f_3 * pc_y[k] * sng_193[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, t_274, pb_y, pc_y, pc_z, smh0_189, smg_119, \
                         smg_134, smg_135, smh1_189, snf0_129, snf1_129, sng_194, \
                         sng_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = f_10 * smg_134[k]
                   + f_3 * pc_y[k] * sng_194[k];

        t_272[k] = f_10 * smg_119[k]
                   + f_1 * snf0_129[k]
                   - f_2 * snf1_129[k]
                   + f_3 * pc_z[k] * sng_194[k];

        t_273[k] = pb_y[k] * smh0_189[k]
                   - f_8 * pc_y[k] * smh1_189[k];

        t_274[k] = f_9 * smg_135[k]
                   + f_3 * pc_y[k] * sng_195[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, pb_y, pc_y, pc_z, smh0_192, smh0_194, \
                         smg_120, smg_136, smg_137, smh1_192, smh1_194, sng_195, \
                         sng_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = f_11 * smg_120[k]
                   + f_3 * pc_z[k] * sng_195[k];

        t_276[k] = pb_y[k] * smh0_192[k]
                   + f_10 * smg_136[k]
                   - f_8 * pc_y[k] * smh1_192[k];

        t_277[k] = f_9 * smg_137[k]
                   + f_3 * pc_y[k] * sng_197[k];

        t_278[k] = pb_y[k] * smh0_194[k]
                   - f_8 * pc_y[k] * smh1_194[k];
    }

#pragma omp simd aligned(t_279, t_280, t_281, t_282, pb_y, pc_y, pc_z, smh0_195, smh0_198, \
                         smg_123, smg_138, smg_140, smh1_195, smh1_198, sng_198, \
                         sng_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_279[k] = pb_y[k] * smh0_195[k]
                   + f_11 * smg_138[k]
                   - f_8 * pc_y[k] * smh1_195[k];

        t_280[k] = f_11 * smg_123[k]
                   + f_3 * pc_z[k] * sng_198[k];

        t_281[k] = f_9 * smg_140[k]
                   + f_3 * pc_y[k] * sng_200[k];

        t_282[k] = pb_y[k] * smh0_198[k]
                   - f_8 * pc_y[k] * smh1_198[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, t_286, t_287, pc_x, smg_205, smg_206, smg_207, \
                         smg_208, smg_209, sng_205, sng_206, sng_207, sng_208, \
                         sng_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_15 * smg_205[k]
                   + f_3 * pc_x[k] * sng_205[k];

        t_284[k] = f_15 * smg_206[k]
                   + f_3 * pc_x[k] * sng_206[k];

        t_285[k] = f_15 * smg_207[k]
                   + f_3 * pc_x[k] * sng_207[k];

        t_286[k] = f_15 * smg_208[k]
                   + f_3 * pc_x[k] * sng_208[k];

        t_287[k] = f_15 * smg_209[k]
                   + f_3 * pc_x[k] * sng_209[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, pc_y, pc_z, smg_130, smg_145, smg_147, snf0_136, \
                         snf0_138, snf1_136, snf1_138, sng_205, \
                         sng_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = f_9 * smg_145[k]
                   + f_1 * snf0_136[k]
                   - f_2 * snf1_136[k]
                   + f_3 * pc_y[k] * sng_205[k];

        t_289[k] = f_11 * smg_130[k]
                   + f_3 * pc_z[k] * sng_205[k];

        t_290[k] = f_9 * smg_147[k]
                   + f_4 * snf0_138[k]
                   - f_5 * snf1_138[k]
                   + f_3 * pc_y[k] * sng_207[k];
    }

#pragma omp simd aligned(t_291, t_292, t_293, pb_y, pc_y, smh0_209, smg_148, smg_149, \
                         smh1_209, snf0_139, snf1_139, sng_208, \
                         sng_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_291[k] = f_9 * smg_148[k]
                   + f_6 * snf0_139[k]
                   - f_7 * snf1_139[k]
                   + f_3 * pc_y[k] * sng_208[k];

        t_292[k] = f_9 * smg_149[k]
                   + f_3 * pc_y[k] * sng_209[k];

        t_293[k] = pb_y[k] * smh0_209[k]
                   - f_8 * pc_y[k] * smh1_209[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, pc_x, pc_y, pc_z, smg_135, smg_210, \
                         smg_213, snf0_140, snf0_143, snf1_140, snf1_143, sng_210, \
                         sng_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_15 * smg_210[k]
                   + f_1 * snf0_140[k]
                   - f_2 * snf1_140[k]
                   + f_3 * pc_x[k] * sng_210[k];

        t_295[k] = f_3 * pc_y[k] * sng_210[k];

        t_296[k] = f_16 * smg_135[k]
                   + f_3 * pc_z[k] * sng_210[k];

        t_297[k] = f_15 * smg_213[k]
                   + f_4 * snf0_143[k]
                   - f_5 * snf1_143[k]
                   + f_3 * pc_x[k] * sng_213[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, pc_x, pc_y, smg_215, smg_216, snf0_145, \
                         snf0_146, snf1_145, snf1_146, sng_212, sng_215, \
                         sng_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = f_3 * pc_y[k] * sng_212[k];

        t_299[k] = f_15 * smg_215[k]
                   + f_4 * snf0_145[k]
                   - f_5 * snf1_145[k]
                   + f_3 * pc_x[k] * sng_215[k];

        t_300[k] = f_15 * smg_216[k]
                   + f_6 * snf0_146[k]
                   - f_7 * snf1_146[k]
                   + f_3 * pc_x[k] * sng_216[k];
    }

#pragma omp simd aligned(t_301, t_302, t_303, t_304, pc_x, pc_y, pc_z, smg_138, smg_219, \
                         smg_220, snf0_149, snf1_149, sng_213, sng_215, sng_219, \
                         sng_220 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_301[k] = f_16 * smg_138[k]
                   + f_3 * pc_z[k] * sng_213[k];

        t_302[k] = f_3 * pc_y[k] * sng_215[k];

        t_303[k] = f_15 * smg_219[k]
                   + f_6 * snf0_149[k]
                   - f_7 * snf1_149[k]
                   + f_3 * pc_x[k] * sng_219[k];

        t_304[k] = f_15 * smg_220[k]
                   + f_3 * pc_x[k] * sng_220[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, pc_x, smg_221, smg_222, smg_223, smg_224, \
                         sng_221, sng_222, sng_223, sng_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = f_15 * smg_221[k]
                   + f_3 * pc_x[k] * sng_221[k];

        t_306[k] = f_15 * smg_222[k]
                   + f_3 * pc_x[k] * sng_222[k];

        t_307[k] = f_15 * smg_223[k]
                   + f_3 * pc_x[k] * sng_223[k];

        t_308[k] = f_15 * smg_224[k]
                   + f_3 * pc_x[k] * sng_224[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, pc_y, pc_z, smg_145, snf0_146, snf0_148, \
                         snf0_149, snf1_146, snf1_148, snf1_149, sng_220, sng_222, \
                         sng_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_1 * snf0_146[k]
                   - f_2 * snf1_146[k]
                   + f_3 * pc_y[k] * sng_220[k];

        t_310[k] = f_16 * smg_145[k]
                   + f_3 * pc_z[k] * sng_220[k];

        t_311[k] = f_4 * snf0_148[k]
                   - f_5 * snf1_148[k]
                   + f_3 * pc_y[k] * sng_222[k];

        t_312[k] = f_6 * snf0_149[k]
                   - f_7 * snf1_149[k]
                   + f_3 * pc_y[k] * sng_223[k];
    }

#pragma omp simd aligned(t_313, t_314, t_315, t_316, pc_x, pc_y, pc_z, smg_149, smg_150, \
                         smg_225, snf0_149, snf0_150, snf1_149, snf1_150, sng_224, \
                         sng_225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_313[k] = f_3 * pc_y[k] * sng_224[k];

        t_314[k] = f_16 * smg_149[k]
                   + f_1 * snf0_149[k]
                   - f_2 * snf1_149[k]
                   + f_3 * pc_z[k] * sng_224[k];

        t_315[k] = f_17 * smg_225[k]
                   + f_1 * snf0_150[k]
                   - f_2 * snf1_150[k]
                   + f_3 * pc_x[k] * sng_225[k];

        t_316[k] = f_17 * smg_150[k]
                   + f_3 * pc_y[k] * sng_225[k];
    }

#pragma omp simd aligned(t_317, t_318, t_319, pc_x, pc_y, pc_z, smg_152, smg_228, snf0_153, \
                         snf1_153, sng_225, sng_227, sng_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_317[k] = f_3 * pc_z[k] * sng_225[k];

        t_318[k] = f_17 * smg_228[k]
                   + f_4 * snf0_153[k]
                   - f_5 * snf1_153[k]
                   + f_3 * pc_x[k] * sng_228[k];

        t_319[k] = f_17 * smg_152[k]
                   + f_3 * pc_y[k] * sng_227[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, pc_x, pc_z, smg_230, smg_231, snf0_155, \
                         snf0_156, snf1_155, snf1_156, sng_228, sng_230, \
                         sng_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = f_17 * smg_230[k]
                   + f_4 * snf0_155[k]
                   - f_5 * snf1_155[k]
                   + f_3 * pc_x[k] * sng_230[k];

        t_321[k] = f_17 * smg_231[k]
                   + f_6 * snf0_156[k]
                   - f_7 * snf1_156[k]
                   + f_3 * pc_x[k] * sng_231[k];

        t_322[k] = f_3 * pc_z[k] * sng_228[k];
    }

#pragma omp simd aligned(t_323, t_324, t_325, t_326, pc_x, pc_y, smg_155, smg_234, smg_235, \
                         smg_236, snf0_159, snf1_159, sng_230, sng_234, sng_235, \
                         sng_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_323[k] = f_17 * smg_155[k]
                   + f_3 * pc_y[k] * sng_230[k];

        t_324[k] = f_17 * smg_234[k]
                   + f_6 * snf0_159[k]
                   - f_7 * snf1_159[k]
                   + f_3 * pc_x[k] * sng_234[k];

        t_325[k] = f_17 * smg_235[k]
                   + f_3 * pc_x[k] * sng_235[k];

        t_326[k] = f_17 * smg_236[k]
                   + f_3 * pc_x[k] * sng_236[k];
    }

#pragma omp simd aligned(t_327, t_328, t_329, t_330, pc_x, pc_y, smg_160, smg_237, smg_238, \
                         smg_239, snf0_156, snf1_156, sng_235, sng_237, sng_238, \
                         sng_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_327[k] = f_17 * smg_237[k]
                   + f_3 * pc_x[k] * sng_237[k];

        t_328[k] = f_17 * smg_238[k]
                   + f_3 * pc_x[k] * sng_238[k];

        t_329[k] = f_17 * smg_239[k]
                   + f_3 * pc_x[k] * sng_239[k];

        t_330[k] = f_17 * smg_160[k]
                   + f_1 * snf0_156[k]
                   - f_2 * snf1_156[k]
                   + f_3 * pc_y[k] * sng_235[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, pc_y, pc_z, smg_162, smg_163, snf0_158, \
                         snf0_159, snf1_158, snf1_159, sng_235, sng_237, \
                         sng_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = f_3 * pc_z[k] * sng_235[k];

        t_332[k] = f_17 * smg_162[k]
                   + f_4 * snf0_158[k]
                   - f_5 * snf1_158[k]
                   + f_3 * pc_y[k] * sng_237[k];

        t_333[k] = f_17 * smg_163[k]
                   + f_6 * snf0_159[k]
                   - f_7 * snf1_159[k]
                   + f_3 * pc_y[k] * sng_238[k];
    }

#pragma omp simd aligned(t_334, t_335, t_336, t_337, pb_z, pc_y, pc_z, smh0_210, smg_164, \
                         smg_165, smh1_210, snf0_159, snf1_159, sng_239, \
                         sng_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_334[k] = f_17 * smg_164[k]
                   + f_3 * pc_y[k] * sng_239[k];

        t_335[k] = f_1 * snf0_159[k]
                   - f_2 * snf1_159[k]
                   + f_3 * pc_z[k] * sng_239[k];

        t_336[k] = pb_z[k] * smh0_210[k]
                   - f_8 * pc_z[k] * smh1_210[k];

        t_337[k] = f_16 * smg_165[k]
                   + f_3 * pc_y[k] * sng_240[k];
    }

#pragma omp simd aligned(t_338, t_339, t_340, pb_z, pc_y, pc_z, smh0_213, smg_150, smg_167, \
                         smh1_213, sng_240, sng_242 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_338[k] = f_9 * smg_150[k]
                   + f_3 * pc_z[k] * sng_240[k];

        t_339[k] = pb_z[k] * smh0_213[k]
                   - f_8 * pc_z[k] * smh1_213[k];

        t_340[k] = f_16 * smg_167[k]
                   + f_3 * pc_y[k] * sng_242[k];
    }

#pragma omp simd aligned(t_341, t_342, t_343, pb_z, pc_x, pc_z, smh0_216, smg_153, smg_245, \
                         smh1_216, snf0_165, snf1_165, sng_243, \
                         sng_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_341[k] = f_17 * smg_245[k]
                   + f_4 * snf0_165[k]
                   - f_5 * snf1_165[k]
                   + f_3 * pc_x[k] * sng_245[k];

        t_342[k] = pb_z[k] * smh0_216[k]
                   - f_8 * pc_z[k] * smh1_216[k];

        t_343[k] = f_9 * smg_153[k]
                   + f_3 * pc_z[k] * sng_243[k];
    }

#pragma omp simd aligned(t_344, t_345, t_346, t_347, pc_x, pc_y, smg_170, smg_249, smg_250, \
                         smg_251, snf0_169, snf1_169, sng_245, sng_249, sng_250, \
                         sng_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_344[k] = f_16 * smg_170[k]
                   + f_3 * pc_y[k] * sng_245[k];

        t_345[k] = f_17 * smg_249[k]
                   + f_6 * snf0_169[k]
                   - f_7 * snf1_169[k]
                   + f_3 * pc_x[k] * sng_249[k];

        t_346[k] = f_17 * smg_250[k]
                   + f_3 * pc_x[k] * sng_250[k];

        t_347[k] = f_17 * smg_251[k]
                   + f_3 * pc_x[k] * sng_251[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, t_351, pb_z, pc_x, pc_z, smh0_225, smg_252, \
                         smg_253, smg_254, smh1_225, sng_252, sng_253, \
                         sng_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = f_17 * smg_252[k]
                   + f_3 * pc_x[k] * sng_252[k];

        t_349[k] = f_17 * smg_253[k]
                   + f_3 * pc_x[k] * sng_253[k];

        t_350[k] = f_17 * smg_254[k]
                   + f_3 * pc_x[k] * sng_254[k];

        t_351[k] = pb_z[k] * smh0_225[k]
                   - f_8 * pc_z[k] * smh1_225[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, pc_y, pc_z, smg_160, smg_177, smg_178, snf0_168, \
                         snf0_169, snf1_168, snf1_169, sng_250, sng_252, \
                         sng_253 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = f_9 * smg_160[k]
                   + f_3 * pc_z[k] * sng_250[k];

        t_353[k] = f_16 * smg_177[k]
                   + f_4 * snf0_168[k]
                   - f_5 * snf1_168[k]
                   + f_3 * pc_y[k] * sng_252[k];

        t_354[k] = f_16 * smg_178[k]
                   + f_6 * snf0_169[k]
                   - f_7 * snf1_169[k]
                   + f_3 * pc_y[k] * sng_253[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, pc_x, pc_y, pc_z, smg_164, smg_179, smg_255, \
                         snf0_169, snf0_170, snf1_169, snf1_170, sng_254, \
                         sng_255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = f_16 * smg_179[k]
                   + f_3 * pc_y[k] * sng_254[k];

        t_356[k] = f_9 * smg_164[k]
                   + f_1 * snf0_169[k]
                   - f_2 * snf1_169[k]
                   + f_3 * pc_z[k] * sng_254[k];

        t_357[k] = f_17 * smg_255[k]
                   + f_1 * snf0_170[k]
                   - f_2 * snf1_170[k]
                   + f_3 * pc_x[k] * sng_255[k];
    }
}

static auto
compute_prim_snh_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t smh0,
                                                          const size_t smg, const size_t smh1,
                                                          const size_t snf0, const size_t snf1,
                                                          const size_t sng, const size_t ncols,
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
    const auto f_15 = 3.0 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;

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

    const auto *smh0_294 = buffer.data(smh0 + 294);
    const auto *smh0_297 = buffer.data(smh0 + 297);
    const auto *smh0_299 = buffer.data(smh0 + 299);
    const auto *smh0_300 = buffer.data(smh0 + 300);
    const auto *smh0_303 = buffer.data(smh0 + 303);
    const auto *smh0_314 = buffer.data(smh0 + 314);
    const auto *smh0_315 = buffer.data(smh0 + 315);
    const auto *smh0_318 = buffer.data(smh0 + 318);
    const auto *smh0_321 = buffer.data(smh0 + 321);

    const auto *smg_165 = buffer.data(smg + 165);
    const auto *smg_168 = buffer.data(smg + 168);
    const auto *smg_175 = buffer.data(smg + 175);
    const auto *smg_179 = buffer.data(smg + 179);
    const auto *smg_180 = buffer.data(smg + 180);
    const auto *smg_182 = buffer.data(smg + 182);
    const auto *smg_183 = buffer.data(smg + 183);
    const auto *smg_185 = buffer.data(smg + 185);
    const auto *smg_190 = buffer.data(smg + 190);
    const auto *smg_192 = buffer.data(smg + 192);
    const auto *smg_193 = buffer.data(smg + 193);
    const auto *smg_194 = buffer.data(smg + 194);
    const auto *smg_195 = buffer.data(smg + 195);
    const auto *smg_197 = buffer.data(smg + 197);
    const auto *smg_198 = buffer.data(smg + 198);
    const auto *smg_200 = buffer.data(smg + 200);
    const auto *smg_205 = buffer.data(smg + 205);
    const auto *smg_207 = buffer.data(smg + 207);
    const auto *smg_208 = buffer.data(smg + 208);
    const auto *smg_209 = buffer.data(smg + 209);
    const auto *smg_210 = buffer.data(smg + 210);
    const auto *smg_211 = buffer.data(smg + 211);
    const auto *smg_212 = buffer.data(smg + 212);
    const auto *smg_213 = buffer.data(smg + 213);
    const auto *smg_215 = buffer.data(smg + 215);
    const auto *smg_220 = buffer.data(smg + 220);
    const auto *smg_222 = buffer.data(smg + 222);
    const auto *smg_223 = buffer.data(smg + 223);
    const auto *smg_224 = buffer.data(smg + 224);
    const auto *smg_225 = buffer.data(smg + 225);
    const auto *smg_227 = buffer.data(smg + 227);
    const auto *smg_228 = buffer.data(smg + 228);
    const auto *smg_230 = buffer.data(smg + 230);
    const auto *smg_235 = buffer.data(smg + 235);
    const auto *smg_237 = buffer.data(smg + 237);
    const auto *smg_238 = buffer.data(smg + 238);
    const auto *smg_239 = buffer.data(smg + 239);
    const auto *smg_240 = buffer.data(smg + 240);
    const auto *smg_242 = buffer.data(smg + 242);
    const auto *smg_245 = buffer.data(smg + 245);
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
    const auto *smg_273 = buffer.data(smg + 273);
    const auto *smg_275 = buffer.data(smg + 275);
    const auto *smg_276 = buffer.data(smg + 276);
    const auto *smg_279 = buffer.data(smg + 279);
    const auto *smg_280 = buffer.data(smg + 280);
    const auto *smg_281 = buffer.data(smg + 281);
    const auto *smg_282 = buffer.data(smg + 282);
    const auto *smg_283 = buffer.data(smg + 283);
    const auto *smg_284 = buffer.data(smg + 284);
    const auto *smg_295 = buffer.data(smg + 295);
    const auto *smg_296 = buffer.data(smg + 296);
    const auto *smg_297 = buffer.data(smg + 297);
    const auto *smg_298 = buffer.data(smg + 298);
    const auto *smg_299 = buffer.data(smg + 299);
    const auto *smg_300 = buffer.data(smg + 300);
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
    const auto *smg_318 = buffer.data(smg + 318);
    const auto *smg_320 = buffer.data(smg + 320);
    const auto *smg_321 = buffer.data(smg + 321);
    const auto *smg_324 = buffer.data(smg + 324);
    const auto *smg_325 = buffer.data(smg + 325);
    const auto *smg_326 = buffer.data(smg + 326);
    const auto *smg_327 = buffer.data(smg + 327);
    const auto *smg_328 = buffer.data(smg + 328);
    const auto *smg_329 = buffer.data(smg + 329);
    const auto *smg_335 = buffer.data(smg + 335);
    const auto *smg_339 = buffer.data(smg + 339);
    const auto *smg_340 = buffer.data(smg + 340);
    const auto *smg_341 = buffer.data(smg + 341);

    const auto *smh1_294 = buffer.data(smh1 + 294);
    const auto *smh1_297 = buffer.data(smh1 + 297);
    const auto *smh1_299 = buffer.data(smh1 + 299);
    const auto *smh1_300 = buffer.data(smh1 + 300);
    const auto *smh1_303 = buffer.data(smh1 + 303);
    const auto *smh1_314 = buffer.data(smh1 + 314);
    const auto *smh1_315 = buffer.data(smh1 + 315);
    const auto *smh1_318 = buffer.data(smh1 + 318);
    const auto *smh1_321 = buffer.data(smh1 + 321);

    const auto *snf0_173 = buffer.data(snf0 + 173);
    const auto *snf0_175 = buffer.data(snf0 + 175);
    const auto *snf0_176 = buffer.data(snf0 + 176);
    const auto *snf0_178 = buffer.data(snf0 + 178);
    const auto *snf0_179 = buffer.data(snf0 + 179);
    const auto *snf0_180 = buffer.data(snf0 + 180);
    const auto *snf0_183 = buffer.data(snf0 + 183);
    const auto *snf0_185 = buffer.data(snf0 + 185);
    const auto *snf0_186 = buffer.data(snf0 + 186);
    const auto *snf0_188 = buffer.data(snf0 + 188);
    const auto *snf0_189 = buffer.data(snf0 + 189);
    const auto *snf0_196 = buffer.data(snf0 + 196);
    const auto *snf0_198 = buffer.data(snf0 + 198);
    const auto *snf0_199 = buffer.data(snf0 + 199);
    const auto *snf0_200 = buffer.data(snf0 + 200);
    const auto *snf0_203 = buffer.data(snf0 + 203);
    const auto *snf0_205 = buffer.data(snf0 + 205);
    const auto *snf0_206 = buffer.data(snf0 + 206);
    const auto *snf0_208 = buffer.data(snf0 + 208);
    const auto *snf0_209 = buffer.data(snf0 + 209);
    const auto *snf0_210 = buffer.data(snf0 + 210);
    const auto *snf0_213 = buffer.data(snf0 + 213);
    const auto *snf0_215 = buffer.data(snf0 + 215);
    const auto *snf0_216 = buffer.data(snf0 + 216);
    const auto *snf0_218 = buffer.data(snf0 + 218);
    const auto *snf0_219 = buffer.data(snf0 + 219);
    const auto *snf0_225 = buffer.data(snf0 + 225);
    const auto *snf0_229 = buffer.data(snf0 + 229);

    const auto *snf1_173 = buffer.data(snf1 + 173);
    const auto *snf1_175 = buffer.data(snf1 + 175);
    const auto *snf1_176 = buffer.data(snf1 + 176);
    const auto *snf1_178 = buffer.data(snf1 + 178);
    const auto *snf1_179 = buffer.data(snf1 + 179);
    const auto *snf1_180 = buffer.data(snf1 + 180);
    const auto *snf1_183 = buffer.data(snf1 + 183);
    const auto *snf1_185 = buffer.data(snf1 + 185);
    const auto *snf1_186 = buffer.data(snf1 + 186);
    const auto *snf1_188 = buffer.data(snf1 + 188);
    const auto *snf1_189 = buffer.data(snf1 + 189);
    const auto *snf1_196 = buffer.data(snf1 + 196);
    const auto *snf1_198 = buffer.data(snf1 + 198);
    const auto *snf1_199 = buffer.data(snf1 + 199);
    const auto *snf1_200 = buffer.data(snf1 + 200);
    const auto *snf1_203 = buffer.data(snf1 + 203);
    const auto *snf1_205 = buffer.data(snf1 + 205);
    const auto *snf1_206 = buffer.data(snf1 + 206);
    const auto *snf1_208 = buffer.data(snf1 + 208);
    const auto *snf1_209 = buffer.data(snf1 + 209);
    const auto *snf1_210 = buffer.data(snf1 + 210);
    const auto *snf1_213 = buffer.data(snf1 + 213);
    const auto *snf1_215 = buffer.data(snf1 + 215);
    const auto *snf1_216 = buffer.data(snf1 + 216);
    const auto *snf1_218 = buffer.data(snf1 + 218);
    const auto *snf1_219 = buffer.data(snf1 + 219);
    const auto *snf1_225 = buffer.data(snf1 + 225);
    const auto *snf1_229 = buffer.data(snf1 + 229);

    const auto *sng_255 = buffer.data(sng + 255);
    const auto *sng_257 = buffer.data(sng + 257);
    const auto *sng_258 = buffer.data(sng + 258);
    const auto *sng_260 = buffer.data(sng + 260);
    const auto *sng_261 = buffer.data(sng + 261);
    const auto *sng_264 = buffer.data(sng + 264);
    const auto *sng_265 = buffer.data(sng + 265);
    const auto *sng_266 = buffer.data(sng + 266);
    const auto *sng_267 = buffer.data(sng + 267);
    const auto *sng_268 = buffer.data(sng + 268);
    const auto *sng_269 = buffer.data(sng + 269);
    const auto *sng_270 = buffer.data(sng + 270);
    const auto *sng_272 = buffer.data(sng + 272);
    const auto *sng_273 = buffer.data(sng + 273);
    const auto *sng_275 = buffer.data(sng + 275);
    const auto *sng_276 = buffer.data(sng + 276);
    const auto *sng_279 = buffer.data(sng + 279);
    const auto *sng_280 = buffer.data(sng + 280);
    const auto *sng_281 = buffer.data(sng + 281);
    const auto *sng_282 = buffer.data(sng + 282);
    const auto *sng_283 = buffer.data(sng + 283);
    const auto *sng_284 = buffer.data(sng + 284);
    const auto *sng_285 = buffer.data(sng + 285);
    const auto *sng_287 = buffer.data(sng + 287);
    const auto *sng_288 = buffer.data(sng + 288);
    const auto *sng_290 = buffer.data(sng + 290);
    const auto *sng_295 = buffer.data(sng + 295);
    const auto *sng_296 = buffer.data(sng + 296);
    const auto *sng_297 = buffer.data(sng + 297);
    const auto *sng_298 = buffer.data(sng + 298);
    const auto *sng_299 = buffer.data(sng + 299);
    const auto *sng_300 = buffer.data(sng + 300);
    const auto *sng_302 = buffer.data(sng + 302);
    const auto *sng_303 = buffer.data(sng + 303);
    const auto *sng_305 = buffer.data(sng + 305);
    const auto *sng_306 = buffer.data(sng + 306);
    const auto *sng_309 = buffer.data(sng + 309);
    const auto *sng_310 = buffer.data(sng + 310);
    const auto *sng_311 = buffer.data(sng + 311);
    const auto *sng_312 = buffer.data(sng + 312);
    const auto *sng_313 = buffer.data(sng + 313);
    const auto *sng_314 = buffer.data(sng + 314);
    const auto *sng_315 = buffer.data(sng + 315);
    const auto *sng_317 = buffer.data(sng + 317);
    const auto *sng_318 = buffer.data(sng + 318);
    const auto *sng_320 = buffer.data(sng + 320);
    const auto *sng_321 = buffer.data(sng + 321);
    const auto *sng_324 = buffer.data(sng + 324);
    const auto *sng_325 = buffer.data(sng + 325);
    const auto *sng_326 = buffer.data(sng + 326);
    const auto *sng_327 = buffer.data(sng + 327);
    const auto *sng_328 = buffer.data(sng + 328);
    const auto *sng_329 = buffer.data(sng + 329);
    const auto *sng_330 = buffer.data(sng + 330);
    const auto *sng_332 = buffer.data(sng + 332);
    const auto *sng_333 = buffer.data(sng + 333);
    const auto *sng_335 = buffer.data(sng + 335);
    const auto *sng_339 = buffer.data(sng + 339);
    const auto *sng_340 = buffer.data(sng + 340);
    const auto *sng_341 = buffer.data(sng + 341);

#pragma omp simd aligned(t_358, t_359, t_360, t_361, pc_x, pc_y, pc_z, smg_165, smg_180, \
                         smg_182, smg_258, snf0_173, snf1_173, sng_255, sng_257, \
                         sng_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_358[k] = f_11 * smg_180[k]
                   + f_3 * pc_y[k] * sng_255[k];

        t_359[k] = f_10 * smg_165[k]
                   + f_3 * pc_z[k] * sng_255[k];

        t_360[k] = f_17 * smg_258[k]
                   + f_4 * snf0_173[k]
                   - f_5 * snf1_173[k]
                   + f_3 * pc_x[k] * sng_258[k];

        t_361[k] = f_11 * smg_182[k]
                   + f_3 * pc_y[k] * sng_257[k];
    }

#pragma omp simd aligned(t_362, t_363, t_364, pc_x, pc_z, smg_168, smg_260, smg_261, snf0_175, \
                         snf0_176, snf1_175, snf1_176, sng_258, sng_260, \
                         sng_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = f_17 * smg_260[k]
                   + f_4 * snf0_175[k]
                   - f_5 * snf1_175[k]
                   + f_3 * pc_x[k] * sng_260[k];

        t_363[k] = f_17 * smg_261[k]
                   + f_6 * snf0_176[k]
                   - f_7 * snf1_176[k]
                   + f_3 * pc_x[k] * sng_261[k];

        t_364[k] = f_10 * smg_168[k]
                   + f_3 * pc_z[k] * sng_258[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, pc_x, pc_y, smg_185, smg_264, smg_265, \
                         smg_266, snf0_179, snf1_179, sng_260, sng_264, sng_265, \
                         sng_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = f_11 * smg_185[k]
                   + f_3 * pc_y[k] * sng_260[k];

        t_366[k] = f_17 * smg_264[k]
                   + f_6 * snf0_179[k]
                   - f_7 * snf1_179[k]
                   + f_3 * pc_x[k] * sng_264[k];

        t_367[k] = f_17 * smg_265[k]
                   + f_3 * pc_x[k] * sng_265[k];

        t_368[k] = f_17 * smg_266[k]
                   + f_3 * pc_x[k] * sng_266[k];
    }

#pragma omp simd aligned(t_369, t_370, t_371, t_372, pc_x, pc_y, smg_190, smg_267, smg_268, \
                         smg_269, snf0_176, snf1_176, sng_265, sng_267, sng_268, \
                         sng_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_369[k] = f_17 * smg_267[k]
                   + f_3 * pc_x[k] * sng_267[k];

        t_370[k] = f_17 * smg_268[k]
                   + f_3 * pc_x[k] * sng_268[k];

        t_371[k] = f_17 * smg_269[k]
                   + f_3 * pc_x[k] * sng_269[k];

        t_372[k] = f_11 * smg_190[k]
                   + f_1 * snf0_176[k]
                   - f_2 * snf1_176[k]
                   + f_3 * pc_y[k] * sng_265[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, pc_y, pc_z, smg_175, smg_192, smg_193, snf0_178, \
                         snf0_179, snf1_178, snf1_179, sng_265, sng_267, \
                         sng_268 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = f_10 * smg_175[k]
                   + f_3 * pc_z[k] * sng_265[k];

        t_374[k] = f_11 * smg_192[k]
                   + f_4 * snf0_178[k]
                   - f_5 * snf1_178[k]
                   + f_3 * pc_y[k] * sng_267[k];

        t_375[k] = f_11 * smg_193[k]
                   + f_6 * snf0_179[k]
                   - f_7 * snf1_179[k]
                   + f_3 * pc_y[k] * sng_268[k];
    }

#pragma omp simd aligned(t_376, t_377, t_378, pc_x, pc_y, pc_z, smg_179, smg_194, smg_270, \
                         snf0_179, snf0_180, snf1_179, snf1_180, sng_269, \
                         sng_270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_376[k] = f_11 * smg_194[k]
                   + f_3 * pc_y[k] * sng_269[k];

        t_377[k] = f_10 * smg_179[k]
                   + f_1 * snf0_179[k]
                   - f_2 * snf1_179[k]
                   + f_3 * pc_z[k] * sng_269[k];

        t_378[k] = f_17 * smg_270[k]
                   + f_1 * snf0_180[k]
                   - f_2 * snf1_180[k]
                   + f_3 * pc_x[k] * sng_270[k];
    }

#pragma omp simd aligned(t_379, t_380, t_381, t_382, pc_x, pc_y, pc_z, smg_180, smg_195, \
                         smg_197, smg_273, snf0_183, snf1_183, sng_270, sng_272, \
                         sng_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_379[k] = f_10 * smg_195[k]
                   + f_3 * pc_y[k] * sng_270[k];

        t_380[k] = f_11 * smg_180[k]
                   + f_3 * pc_z[k] * sng_270[k];

        t_381[k] = f_17 * smg_273[k]
                   + f_4 * snf0_183[k]
                   - f_5 * snf1_183[k]
                   + f_3 * pc_x[k] * sng_273[k];

        t_382[k] = f_10 * smg_197[k]
                   + f_3 * pc_y[k] * sng_272[k];
    }

#pragma omp simd aligned(t_383, t_384, t_385, pc_x, pc_z, smg_183, smg_275, smg_276, snf0_185, \
                         snf0_186, snf1_185, snf1_186, sng_273, sng_275, \
                         sng_276 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_383[k] = f_17 * smg_275[k]
                   + f_4 * snf0_185[k]
                   - f_5 * snf1_185[k]
                   + f_3 * pc_x[k] * sng_275[k];

        t_384[k] = f_17 * smg_276[k]
                   + f_6 * snf0_186[k]
                   - f_7 * snf1_186[k]
                   + f_3 * pc_x[k] * sng_276[k];

        t_385[k] = f_11 * smg_183[k]
                   + f_3 * pc_z[k] * sng_273[k];
    }

#pragma omp simd aligned(t_386, t_387, t_388, t_389, pc_x, pc_y, smg_200, smg_279, smg_280, \
                         smg_281, snf0_189, snf1_189, sng_275, sng_279, sng_280, \
                         sng_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_386[k] = f_10 * smg_200[k]
                   + f_3 * pc_y[k] * sng_275[k];

        t_387[k] = f_17 * smg_279[k]
                   + f_6 * snf0_189[k]
                   - f_7 * snf1_189[k]
                   + f_3 * pc_x[k] * sng_279[k];

        t_388[k] = f_17 * smg_280[k]
                   + f_3 * pc_x[k] * sng_280[k];

        t_389[k] = f_17 * smg_281[k]
                   + f_3 * pc_x[k] * sng_281[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, pc_x, pc_y, smg_205, smg_282, smg_283, \
                         smg_284, snf0_186, snf1_186, sng_280, sng_282, sng_283, \
                         sng_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = f_17 * smg_282[k]
                   + f_3 * pc_x[k] * sng_282[k];

        t_391[k] = f_17 * smg_283[k]
                   + f_3 * pc_x[k] * sng_283[k];

        t_392[k] = f_17 * smg_284[k]
                   + f_3 * pc_x[k] * sng_284[k];

        t_393[k] = f_10 * smg_205[k]
                   + f_1 * snf0_186[k]
                   - f_2 * snf1_186[k]
                   + f_3 * pc_y[k] * sng_280[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, pc_y, pc_z, smg_190, smg_207, smg_208, snf0_188, \
                         snf0_189, snf1_188, snf1_189, sng_280, sng_282, \
                         sng_283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = f_11 * smg_190[k]
                   + f_3 * pc_z[k] * sng_280[k];

        t_395[k] = f_10 * smg_207[k]
                   + f_4 * snf0_188[k]
                   - f_5 * snf1_188[k]
                   + f_3 * pc_y[k] * sng_282[k];

        t_396[k] = f_10 * smg_208[k]
                   + f_6 * snf0_189[k]
                   - f_7 * snf1_189[k]
                   + f_3 * pc_y[k] * sng_283[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, t_400, pb_y, pc_y, pc_z, smh0_294, smg_194, \
                         smg_209, smg_210, smh1_294, snf0_189, snf1_189, sng_284, \
                         sng_285 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = f_10 * smg_209[k]
                   + f_3 * pc_y[k] * sng_284[k];

        t_398[k] = f_11 * smg_194[k]
                   + f_1 * snf0_189[k]
                   - f_2 * snf1_189[k]
                   + f_3 * pc_z[k] * sng_284[k];

        t_399[k] = pb_y[k] * smh0_294[k]
                   - f_8 * pc_y[k] * smh1_294[k];

        t_400[k] = f_9 * smg_210[k]
                   + f_3 * pc_y[k] * sng_285[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, t_404, pb_y, pc_y, pc_z, smh0_297, smh0_299, \
                         smg_195, smg_211, smg_212, smh1_297, smh1_299, sng_285, \
                         sng_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_16 * smg_195[k]
                   + f_3 * pc_z[k] * sng_285[k];

        t_402[k] = pb_y[k] * smh0_297[k]
                   + f_10 * smg_211[k]
                   - f_8 * pc_y[k] * smh1_297[k];

        t_403[k] = f_9 * smg_212[k]
                   + f_3 * pc_y[k] * sng_287[k];

        t_404[k] = pb_y[k] * smh0_299[k]
                   - f_8 * pc_y[k] * smh1_299[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, pb_y, pc_y, pc_z, smh0_300, smh0_303, \
                         smg_198, smg_213, smg_215, smh1_300, smh1_303, sng_288, \
                         sng_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = pb_y[k] * smh0_300[k]
                   + f_11 * smg_213[k]
                   - f_8 * pc_y[k] * smh1_300[k];

        t_406[k] = f_16 * smg_198[k]
                   + f_3 * pc_z[k] * sng_288[k];

        t_407[k] = f_9 * smg_215[k]
                   + f_3 * pc_y[k] * sng_290[k];

        t_408[k] = pb_y[k] * smh0_303[k]
                   - f_8 * pc_y[k] * smh1_303[k];
    }

#pragma omp simd aligned(t_409, t_410, t_411, t_412, t_413, pc_x, smg_295, smg_296, smg_297, \
                         smg_298, smg_299, sng_295, sng_296, sng_297, sng_298, \
                         sng_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_409[k] = f_17 * smg_295[k]
                   + f_3 * pc_x[k] * sng_295[k];

        t_410[k] = f_17 * smg_296[k]
                   + f_3 * pc_x[k] * sng_296[k];

        t_411[k] = f_17 * smg_297[k]
                   + f_3 * pc_x[k] * sng_297[k];

        t_412[k] = f_17 * smg_298[k]
                   + f_3 * pc_x[k] * sng_298[k];

        t_413[k] = f_17 * smg_299[k]
                   + f_3 * pc_x[k] * sng_299[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, pc_y, pc_z, smg_205, smg_220, smg_222, snf0_196, \
                         snf0_198, snf1_196, snf1_198, sng_295, \
                         sng_297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_9 * smg_220[k]
                   + f_1 * snf0_196[k]
                   - f_2 * snf1_196[k]
                   + f_3 * pc_y[k] * sng_295[k];

        t_415[k] = f_16 * smg_205[k]
                   + f_3 * pc_z[k] * sng_295[k];

        t_416[k] = f_9 * smg_222[k]
                   + f_4 * snf0_198[k]
                   - f_5 * snf1_198[k]
                   + f_3 * pc_y[k] * sng_297[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, pb_y, pc_y, smh0_314, smg_223, smg_224, \
                         smh1_314, snf0_199, snf1_199, sng_298, \
                         sng_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = f_9 * smg_223[k]
                   + f_6 * snf0_199[k]
                   - f_7 * snf1_199[k]
                   + f_3 * pc_y[k] * sng_298[k];

        t_418[k] = f_9 * smg_224[k]
                   + f_3 * pc_y[k] * sng_299[k];

        t_419[k] = pb_y[k] * smh0_314[k]
                   - f_8 * pc_y[k] * smh1_314[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, pc_x, pc_y, pc_z, smg_210, smg_300, \
                         smg_303, snf0_200, snf0_203, snf1_200, snf1_203, sng_300, \
                         sng_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_17 * smg_300[k]
                   + f_1 * snf0_200[k]
                   - f_2 * snf1_200[k]
                   + f_3 * pc_x[k] * sng_300[k];

        t_421[k] = f_3 * pc_y[k] * sng_300[k];

        t_422[k] = f_17 * smg_210[k]
                   + f_3 * pc_z[k] * sng_300[k];

        t_423[k] = f_17 * smg_303[k]
                   + f_4 * snf0_203[k]
                   - f_5 * snf1_203[k]
                   + f_3 * pc_x[k] * sng_303[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, pc_x, pc_y, smg_305, smg_306, snf0_205, \
                         snf0_206, snf1_205, snf1_206, sng_302, sng_305, \
                         sng_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = f_3 * pc_y[k] * sng_302[k];

        t_425[k] = f_17 * smg_305[k]
                   + f_4 * snf0_205[k]
                   - f_5 * snf1_205[k]
                   + f_3 * pc_x[k] * sng_305[k];

        t_426[k] = f_17 * smg_306[k]
                   + f_6 * snf0_206[k]
                   - f_7 * snf1_206[k]
                   + f_3 * pc_x[k] * sng_306[k];
    }

#pragma omp simd aligned(t_427, t_428, t_429, t_430, pc_x, pc_y, pc_z, smg_213, smg_309, \
                         smg_310, snf0_209, snf1_209, sng_303, sng_305, sng_309, \
                         sng_310 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_427[k] = f_17 * smg_213[k]
                   + f_3 * pc_z[k] * sng_303[k];

        t_428[k] = f_3 * pc_y[k] * sng_305[k];

        t_429[k] = f_17 * smg_309[k]
                   + f_6 * snf0_209[k]
                   - f_7 * snf1_209[k]
                   + f_3 * pc_x[k] * sng_309[k];

        t_430[k] = f_17 * smg_310[k]
                   + f_3 * pc_x[k] * sng_310[k];
    }

#pragma omp simd aligned(t_431, t_432, t_433, t_434, pc_x, smg_311, smg_312, smg_313, smg_314, \
                         sng_311, sng_312, sng_313, sng_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_431[k] = f_17 * smg_311[k]
                   + f_3 * pc_x[k] * sng_311[k];

        t_432[k] = f_17 * smg_312[k]
                   + f_3 * pc_x[k] * sng_312[k];

        t_433[k] = f_17 * smg_313[k]
                   + f_3 * pc_x[k] * sng_313[k];

        t_434[k] = f_17 * smg_314[k]
                   + f_3 * pc_x[k] * sng_314[k];
    }

#pragma omp simd aligned(t_435, t_436, t_437, t_438, pc_y, pc_z, smg_220, snf0_206, snf0_208, \
                         snf0_209, snf1_206, snf1_208, snf1_209, sng_310, sng_312, \
                         sng_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_435[k] = f_1 * snf0_206[k]
                   - f_2 * snf1_206[k]
                   + f_3 * pc_y[k] * sng_310[k];

        t_436[k] = f_17 * smg_220[k]
                   + f_3 * pc_z[k] * sng_310[k];

        t_437[k] = f_4 * snf0_208[k]
                   - f_5 * snf1_208[k]
                   + f_3 * pc_y[k] * sng_312[k];

        t_438[k] = f_6 * snf0_209[k]
                   - f_7 * snf1_209[k]
                   + f_3 * pc_y[k] * sng_313[k];
    }

#pragma omp simd aligned(t_439, t_440, t_441, t_442, pc_x, pc_y, pc_z, smg_224, smg_225, \
                         smg_315, snf0_209, snf0_210, snf1_209, snf1_210, sng_314, \
                         sng_315 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_439[k] = f_3 * pc_y[k] * sng_314[k];

        t_440[k] = f_17 * smg_224[k]
                   + f_1 * snf0_209[k]
                   - f_2 * snf1_209[k]
                   + f_3 * pc_z[k] * sng_314[k];

        t_441[k] = f_16 * smg_315[k]
                   + f_1 * snf0_210[k]
                   - f_2 * snf1_210[k]
                   + f_3 * pc_x[k] * sng_315[k];

        t_442[k] = f_15 * smg_225[k]
                   + f_3 * pc_y[k] * sng_315[k];
    }

#pragma omp simd aligned(t_443, t_444, t_445, pc_x, pc_y, pc_z, smg_227, smg_318, snf0_213, \
                         snf1_213, sng_315, sng_317, sng_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_443[k] = f_3 * pc_z[k] * sng_315[k];

        t_444[k] = f_16 * smg_318[k]
                   + f_4 * snf0_213[k]
                   - f_5 * snf1_213[k]
                   + f_3 * pc_x[k] * sng_318[k];

        t_445[k] = f_15 * smg_227[k]
                   + f_3 * pc_y[k] * sng_317[k];
    }

#pragma omp simd aligned(t_446, t_447, t_448, pc_x, pc_z, smg_320, smg_321, snf0_215, \
                         snf0_216, snf1_215, snf1_216, sng_318, sng_320, \
                         sng_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_446[k] = f_16 * smg_320[k]
                   + f_4 * snf0_215[k]
                   - f_5 * snf1_215[k]
                   + f_3 * pc_x[k] * sng_320[k];

        t_447[k] = f_16 * smg_321[k]
                   + f_6 * snf0_216[k]
                   - f_7 * snf1_216[k]
                   + f_3 * pc_x[k] * sng_321[k];

        t_448[k] = f_3 * pc_z[k] * sng_318[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, t_452, pc_x, pc_y, smg_230, smg_324, smg_325, \
                         smg_326, snf0_219, snf1_219, sng_320, sng_324, sng_325, \
                         sng_326 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_15 * smg_230[k]
                   + f_3 * pc_y[k] * sng_320[k];

        t_450[k] = f_16 * smg_324[k]
                   + f_6 * snf0_219[k]
                   - f_7 * snf1_219[k]
                   + f_3 * pc_x[k] * sng_324[k];

        t_451[k] = f_16 * smg_325[k]
                   + f_3 * pc_x[k] * sng_325[k];

        t_452[k] = f_16 * smg_326[k]
                   + f_3 * pc_x[k] * sng_326[k];
    }

#pragma omp simd aligned(t_453, t_454, t_455, t_456, pc_x, pc_y, smg_235, smg_327, smg_328, \
                         smg_329, snf0_216, snf1_216, sng_325, sng_327, sng_328, \
                         sng_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_453[k] = f_16 * smg_327[k]
                   + f_3 * pc_x[k] * sng_327[k];

        t_454[k] = f_16 * smg_328[k]
                   + f_3 * pc_x[k] * sng_328[k];

        t_455[k] = f_16 * smg_329[k]
                   + f_3 * pc_x[k] * sng_329[k];

        t_456[k] = f_15 * smg_235[k]
                   + f_1 * snf0_216[k]
                   - f_2 * snf1_216[k]
                   + f_3 * pc_y[k] * sng_325[k];
    }

#pragma omp simd aligned(t_457, t_458, t_459, pc_y, pc_z, smg_237, smg_238, snf0_218, \
                         snf0_219, snf1_218, snf1_219, sng_325, sng_327, \
                         sng_328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_457[k] = f_3 * pc_z[k] * sng_325[k];

        t_458[k] = f_15 * smg_237[k]
                   + f_4 * snf0_218[k]
                   - f_5 * snf1_218[k]
                   + f_3 * pc_y[k] * sng_327[k];

        t_459[k] = f_15 * smg_238[k]
                   + f_6 * snf0_219[k]
                   - f_7 * snf1_219[k]
                   + f_3 * pc_y[k] * sng_328[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, t_463, pb_z, pc_y, pc_z, smh0_315, smg_239, \
                         smg_240, smh1_315, snf0_219, snf1_219, sng_329, \
                         sng_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = f_15 * smg_239[k]
                   + f_3 * pc_y[k] * sng_329[k];

        t_461[k] = f_1 * snf0_219[k]
                   - f_2 * snf1_219[k]
                   + f_3 * pc_z[k] * sng_329[k];

        t_462[k] = pb_z[k] * smh0_315[k]
                   - f_8 * pc_z[k] * smh1_315[k];

        t_463[k] = f_17 * smg_240[k]
                   + f_3 * pc_y[k] * sng_330[k];
    }

#pragma omp simd aligned(t_464, t_465, t_466, pb_z, pc_y, pc_z, smh0_318, smg_225, smg_242, \
                         smh1_318, sng_330, sng_332 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_464[k] = f_9 * smg_225[k]
                   + f_3 * pc_z[k] * sng_330[k];

        t_465[k] = pb_z[k] * smh0_318[k]
                   - f_8 * pc_z[k] * smh1_318[k];

        t_466[k] = f_17 * smg_242[k]
                   + f_3 * pc_y[k] * sng_332[k];
    }

#pragma omp simd aligned(t_467, t_468, t_469, pb_z, pc_x, pc_z, smh0_321, smg_228, smg_335, \
                         smh1_321, snf0_225, snf1_225, sng_333, \
                         sng_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_467[k] = f_16 * smg_335[k]
                   + f_4 * snf0_225[k]
                   - f_5 * snf1_225[k]
                   + f_3 * pc_x[k] * sng_335[k];

        t_468[k] = pb_z[k] * smh0_321[k]
                   - f_8 * pc_z[k] * smh1_321[k];

        t_469[k] = f_9 * smg_228[k]
                   + f_3 * pc_z[k] * sng_333[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, pc_x, pc_y, smg_245, smg_339, smg_340, \
                         smg_341, snf0_229, snf1_229, sng_335, sng_339, sng_340, \
                         sng_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = f_17 * smg_245[k]
                   + f_3 * pc_y[k] * sng_335[k];

        t_471[k] = f_16 * smg_339[k]
                   + f_6 * snf0_229[k]
                   - f_7 * snf1_229[k]
                   + f_3 * pc_x[k] * sng_339[k];

        t_472[k] = f_16 * smg_340[k]
                   + f_3 * pc_x[k] * sng_340[k];

        t_473[k] = f_16 * smg_341[k]
                   + f_3 * pc_x[k] * sng_341[k];
    }
}

static auto
compute_prim_snh_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t smh0,
                                                          const size_t smg, const size_t smh1,
                                                          const size_t snf0, const size_t snf1,
                                                          const size_t sng, const size_t ncols,
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
    const auto f_15 = 3.0 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;

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

    const auto *smh0_330 = buffer.data(smh0 + 330);
    const auto *smh0_420 = buffer.data(smh0 + 420);
    const auto *smh0_423 = buffer.data(smh0 + 423);
    const auto *smh0_425 = buffer.data(smh0 + 425);
    const auto *smh0_426 = buffer.data(smh0 + 426);
    const auto *smh0_429 = buffer.data(smh0 + 429);
    const auto *smh0_440 = buffer.data(smh0 + 440);

    const auto *smg_235 = buffer.data(smg + 235);
    const auto *smg_239 = buffer.data(smg + 239);
    const auto *smg_240 = buffer.data(smg + 240);
    const auto *smg_243 = buffer.data(smg + 243);
    const auto *smg_250 = buffer.data(smg + 250);
    const auto *smg_252 = buffer.data(smg + 252);
    const auto *smg_253 = buffer.data(smg + 253);
    const auto *smg_254 = buffer.data(smg + 254);
    const auto *smg_255 = buffer.data(smg + 255);
    const auto *smg_257 = buffer.data(smg + 257);
    const auto *smg_258 = buffer.data(smg + 258);
    const auto *smg_260 = buffer.data(smg + 260);
    const auto *smg_265 = buffer.data(smg + 265);
    const auto *smg_267 = buffer.data(smg + 267);
    const auto *smg_268 = buffer.data(smg + 268);
    const auto *smg_269 = buffer.data(smg + 269);
    const auto *smg_270 = buffer.data(smg + 270);
    const auto *smg_272 = buffer.data(smg + 272);
    const auto *smg_273 = buffer.data(smg + 273);
    const auto *smg_275 = buffer.data(smg + 275);
    const auto *smg_280 = buffer.data(smg + 280);
    const auto *smg_282 = buffer.data(smg + 282);
    const auto *smg_283 = buffer.data(smg + 283);
    const auto *smg_284 = buffer.data(smg + 284);
    const auto *smg_285 = buffer.data(smg + 285);
    const auto *smg_287 = buffer.data(smg + 287);
    const auto *smg_288 = buffer.data(smg + 288);
    const auto *smg_290 = buffer.data(smg + 290);
    const auto *smg_295 = buffer.data(smg + 295);
    const auto *smg_297 = buffer.data(smg + 297);
    const auto *smg_298 = buffer.data(smg + 298);
    const auto *smg_299 = buffer.data(smg + 299);
    const auto *smg_300 = buffer.data(smg + 300);
    const auto *smg_301 = buffer.data(smg + 301);
    const auto *smg_302 = buffer.data(smg + 302);
    const auto *smg_303 = buffer.data(smg + 303);
    const auto *smg_305 = buffer.data(smg + 305);
    const auto *smg_310 = buffer.data(smg + 310);
    const auto *smg_312 = buffer.data(smg + 312);
    const auto *smg_313 = buffer.data(smg + 313);
    const auto *smg_314 = buffer.data(smg + 314);
    const auto *smg_342 = buffer.data(smg + 342);
    const auto *smg_343 = buffer.data(smg + 343);
    const auto *smg_344 = buffer.data(smg + 344);
    const auto *smg_345 = buffer.data(smg + 345);
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
    const auto *smg_378 = buffer.data(smg + 378);
    const auto *smg_380 = buffer.data(smg + 380);
    const auto *smg_381 = buffer.data(smg + 381);
    const auto *smg_384 = buffer.data(smg + 384);
    const auto *smg_385 = buffer.data(smg + 385);
    const auto *smg_386 = buffer.data(smg + 386);
    const auto *smg_387 = buffer.data(smg + 387);
    const auto *smg_388 = buffer.data(smg + 388);
    const auto *smg_389 = buffer.data(smg + 389);
    const auto *smg_400 = buffer.data(smg + 400);
    const auto *smg_401 = buffer.data(smg + 401);
    const auto *smg_402 = buffer.data(smg + 402);
    const auto *smg_403 = buffer.data(smg + 403);
    const auto *smg_404 = buffer.data(smg + 404);
    const auto *smg_405 = buffer.data(smg + 405);
    const auto *smg_408 = buffer.data(smg + 408);
    const auto *smg_410 = buffer.data(smg + 410);
    const auto *smg_411 = buffer.data(smg + 411);
    const auto *smg_414 = buffer.data(smg + 414);
    const auto *smg_415 = buffer.data(smg + 415);
    const auto *smg_416 = buffer.data(smg + 416);
    const auto *smg_417 = buffer.data(smg + 417);
    const auto *smg_418 = buffer.data(smg + 418);
    const auto *smg_419 = buffer.data(smg + 419);

    const auto *smh1_330 = buffer.data(smh1 + 330);
    const auto *smh1_420 = buffer.data(smh1 + 420);
    const auto *smh1_423 = buffer.data(smh1 + 423);
    const auto *smh1_425 = buffer.data(smh1 + 425);
    const auto *smh1_426 = buffer.data(smh1 + 426);
    const auto *smh1_429 = buffer.data(smh1 + 429);
    const auto *smh1_440 = buffer.data(smh1 + 440);

    const auto *snf0_228 = buffer.data(snf0 + 228);
    const auto *snf0_229 = buffer.data(snf0 + 229);
    const auto *snf0_230 = buffer.data(snf0 + 230);
    const auto *snf0_233 = buffer.data(snf0 + 233);
    const auto *snf0_235 = buffer.data(snf0 + 235);
    const auto *snf0_236 = buffer.data(snf0 + 236);
    const auto *snf0_238 = buffer.data(snf0 + 238);
    const auto *snf0_239 = buffer.data(snf0 + 239);
    const auto *snf0_240 = buffer.data(snf0 + 240);
    const auto *snf0_243 = buffer.data(snf0 + 243);
    const auto *snf0_245 = buffer.data(snf0 + 245);
    const auto *snf0_246 = buffer.data(snf0 + 246);
    const auto *snf0_248 = buffer.data(snf0 + 248);
    const auto *snf0_249 = buffer.data(snf0 + 249);
    const auto *snf0_250 = buffer.data(snf0 + 250);
    const auto *snf0_253 = buffer.data(snf0 + 253);
    const auto *snf0_255 = buffer.data(snf0 + 255);
    const auto *snf0_256 = buffer.data(snf0 + 256);
    const auto *snf0_258 = buffer.data(snf0 + 258);
    const auto *snf0_259 = buffer.data(snf0 + 259);
    const auto *snf0_266 = buffer.data(snf0 + 266);
    const auto *snf0_268 = buffer.data(snf0 + 268);
    const auto *snf0_269 = buffer.data(snf0 + 269);
    const auto *snf0_270 = buffer.data(snf0 + 270);
    const auto *snf0_273 = buffer.data(snf0 + 273);
    const auto *snf0_275 = buffer.data(snf0 + 275);
    const auto *snf0_276 = buffer.data(snf0 + 276);
    const auto *snf0_278 = buffer.data(snf0 + 278);
    const auto *snf0_279 = buffer.data(snf0 + 279);

    const auto *snf1_228 = buffer.data(snf1 + 228);
    const auto *snf1_229 = buffer.data(snf1 + 229);
    const auto *snf1_230 = buffer.data(snf1 + 230);
    const auto *snf1_233 = buffer.data(snf1 + 233);
    const auto *snf1_235 = buffer.data(snf1 + 235);
    const auto *snf1_236 = buffer.data(snf1 + 236);
    const auto *snf1_238 = buffer.data(snf1 + 238);
    const auto *snf1_239 = buffer.data(snf1 + 239);
    const auto *snf1_240 = buffer.data(snf1 + 240);
    const auto *snf1_243 = buffer.data(snf1 + 243);
    const auto *snf1_245 = buffer.data(snf1 + 245);
    const auto *snf1_246 = buffer.data(snf1 + 246);
    const auto *snf1_248 = buffer.data(snf1 + 248);
    const auto *snf1_249 = buffer.data(snf1 + 249);
    const auto *snf1_250 = buffer.data(snf1 + 250);
    const auto *snf1_253 = buffer.data(snf1 + 253);
    const auto *snf1_255 = buffer.data(snf1 + 255);
    const auto *snf1_256 = buffer.data(snf1 + 256);
    const auto *snf1_258 = buffer.data(snf1 + 258);
    const auto *snf1_259 = buffer.data(snf1 + 259);
    const auto *snf1_266 = buffer.data(snf1 + 266);
    const auto *snf1_268 = buffer.data(snf1 + 268);
    const auto *snf1_269 = buffer.data(snf1 + 269);
    const auto *snf1_270 = buffer.data(snf1 + 270);
    const auto *snf1_273 = buffer.data(snf1 + 273);
    const auto *snf1_275 = buffer.data(snf1 + 275);
    const auto *snf1_276 = buffer.data(snf1 + 276);
    const auto *snf1_278 = buffer.data(snf1 + 278);
    const auto *snf1_279 = buffer.data(snf1 + 279);

    const auto *sng_340 = buffer.data(sng + 340);
    const auto *sng_342 = buffer.data(sng + 342);
    const auto *sng_343 = buffer.data(sng + 343);
    const auto *sng_344 = buffer.data(sng + 344);
    const auto *sng_345 = buffer.data(sng + 345);
    const auto *sng_347 = buffer.data(sng + 347);
    const auto *sng_348 = buffer.data(sng + 348);
    const auto *sng_350 = buffer.data(sng + 350);
    const auto *sng_351 = buffer.data(sng + 351);
    const auto *sng_354 = buffer.data(sng + 354);
    const auto *sng_355 = buffer.data(sng + 355);
    const auto *sng_356 = buffer.data(sng + 356);
    const auto *sng_357 = buffer.data(sng + 357);
    const auto *sng_358 = buffer.data(sng + 358);
    const auto *sng_359 = buffer.data(sng + 359);
    const auto *sng_360 = buffer.data(sng + 360);
    const auto *sng_362 = buffer.data(sng + 362);
    const auto *sng_363 = buffer.data(sng + 363);
    const auto *sng_365 = buffer.data(sng + 365);
    const auto *sng_366 = buffer.data(sng + 366);
    const auto *sng_369 = buffer.data(sng + 369);
    const auto *sng_370 = buffer.data(sng + 370);
    const auto *sng_371 = buffer.data(sng + 371);
    const auto *sng_372 = buffer.data(sng + 372);
    const auto *sng_373 = buffer.data(sng + 373);
    const auto *sng_374 = buffer.data(sng + 374);
    const auto *sng_375 = buffer.data(sng + 375);
    const auto *sng_377 = buffer.data(sng + 377);
    const auto *sng_378 = buffer.data(sng + 378);
    const auto *sng_380 = buffer.data(sng + 380);
    const auto *sng_381 = buffer.data(sng + 381);
    const auto *sng_384 = buffer.data(sng + 384);
    const auto *sng_385 = buffer.data(sng + 385);
    const auto *sng_386 = buffer.data(sng + 386);
    const auto *sng_387 = buffer.data(sng + 387);
    const auto *sng_388 = buffer.data(sng + 388);
    const auto *sng_389 = buffer.data(sng + 389);
    const auto *sng_390 = buffer.data(sng + 390);
    const auto *sng_392 = buffer.data(sng + 392);
    const auto *sng_393 = buffer.data(sng + 393);
    const auto *sng_395 = buffer.data(sng + 395);
    const auto *sng_400 = buffer.data(sng + 400);
    const auto *sng_401 = buffer.data(sng + 401);
    const auto *sng_402 = buffer.data(sng + 402);
    const auto *sng_403 = buffer.data(sng + 403);
    const auto *sng_404 = buffer.data(sng + 404);
    const auto *sng_405 = buffer.data(sng + 405);
    const auto *sng_407 = buffer.data(sng + 407);
    const auto *sng_408 = buffer.data(sng + 408);
    const auto *sng_410 = buffer.data(sng + 410);
    const auto *sng_411 = buffer.data(sng + 411);
    const auto *sng_414 = buffer.data(sng + 414);
    const auto *sng_415 = buffer.data(sng + 415);
    const auto *sng_416 = buffer.data(sng + 416);
    const auto *sng_417 = buffer.data(sng + 417);
    const auto *sng_418 = buffer.data(sng + 418);
    const auto *sng_419 = buffer.data(sng + 419);

#pragma omp simd aligned(t_474, t_475, t_476, t_477, pb_z, pc_x, pc_z, smh0_330, smg_342, \
                         smg_343, smg_344, smh1_330, sng_342, sng_343, \
                         sng_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = f_16 * smg_342[k]
                   + f_3 * pc_x[k] * sng_342[k];

        t_475[k] = f_16 * smg_343[k]
                   + f_3 * pc_x[k] * sng_343[k];

        t_476[k] = f_16 * smg_344[k]
                   + f_3 * pc_x[k] * sng_344[k];

        t_477[k] = pb_z[k] * smh0_330[k]
                   - f_8 * pc_z[k] * smh1_330[k];
    }

#pragma omp simd aligned(t_478, t_479, t_480, pc_y, pc_z, smg_235, smg_252, smg_253, snf0_228, \
                         snf0_229, snf1_228, snf1_229, sng_340, sng_342, \
                         sng_343 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_478[k] = f_9 * smg_235[k]
                   + f_3 * pc_z[k] * sng_340[k];

        t_479[k] = f_17 * smg_252[k]
                   + f_4 * snf0_228[k]
                   - f_5 * snf1_228[k]
                   + f_3 * pc_y[k] * sng_342[k];

        t_480[k] = f_17 * smg_253[k]
                   + f_6 * snf0_229[k]
                   - f_7 * snf1_229[k]
                   + f_3 * pc_y[k] * sng_343[k];
    }

#pragma omp simd aligned(t_481, t_482, t_483, pc_x, pc_y, pc_z, smg_239, smg_254, smg_345, \
                         snf0_229, snf0_230, snf1_229, snf1_230, sng_344, \
                         sng_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_481[k] = f_17 * smg_254[k]
                   + f_3 * pc_y[k] * sng_344[k];

        t_482[k] = f_9 * smg_239[k]
                   + f_1 * snf0_229[k]
                   - f_2 * snf1_229[k]
                   + f_3 * pc_z[k] * sng_344[k];

        t_483[k] = f_16 * smg_345[k]
                   + f_1 * snf0_230[k]
                   - f_2 * snf1_230[k]
                   + f_3 * pc_x[k] * sng_345[k];
    }

#pragma omp simd aligned(t_484, t_485, t_486, t_487, pc_x, pc_y, pc_z, smg_240, smg_255, \
                         smg_257, smg_348, snf0_233, snf1_233, sng_345, sng_347, \
                         sng_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_484[k] = f_16 * smg_255[k]
                   + f_3 * pc_y[k] * sng_345[k];

        t_485[k] = f_10 * smg_240[k]
                   + f_3 * pc_z[k] * sng_345[k];

        t_486[k] = f_16 * smg_348[k]
                   + f_4 * snf0_233[k]
                   - f_5 * snf1_233[k]
                   + f_3 * pc_x[k] * sng_348[k];

        t_487[k] = f_16 * smg_257[k]
                   + f_3 * pc_y[k] * sng_347[k];
    }

#pragma omp simd aligned(t_488, t_489, t_490, pc_x, pc_z, smg_243, smg_350, smg_351, snf0_235, \
                         snf0_236, snf1_235, snf1_236, sng_348, sng_350, \
                         sng_351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_488[k] = f_16 * smg_350[k]
                   + f_4 * snf0_235[k]
                   - f_5 * snf1_235[k]
                   + f_3 * pc_x[k] * sng_350[k];

        t_489[k] = f_16 * smg_351[k]
                   + f_6 * snf0_236[k]
                   - f_7 * snf1_236[k]
                   + f_3 * pc_x[k] * sng_351[k];

        t_490[k] = f_10 * smg_243[k]
                   + f_3 * pc_z[k] * sng_348[k];
    }

#pragma omp simd aligned(t_491, t_492, t_493, t_494, pc_x, pc_y, smg_260, smg_354, smg_355, \
                         smg_356, snf0_239, snf1_239, sng_350, sng_354, sng_355, \
                         sng_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_491[k] = f_16 * smg_260[k]
                   + f_3 * pc_y[k] * sng_350[k];

        t_492[k] = f_16 * smg_354[k]
                   + f_6 * snf0_239[k]
                   - f_7 * snf1_239[k]
                   + f_3 * pc_x[k] * sng_354[k];

        t_493[k] = f_16 * smg_355[k]
                   + f_3 * pc_x[k] * sng_355[k];

        t_494[k] = f_16 * smg_356[k]
                   + f_3 * pc_x[k] * sng_356[k];
    }

#pragma omp simd aligned(t_495, t_496, t_497, t_498, pc_x, pc_y, smg_265, smg_357, smg_358, \
                         smg_359, snf0_236, snf1_236, sng_355, sng_357, sng_358, \
                         sng_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_495[k] = f_16 * smg_357[k]
                   + f_3 * pc_x[k] * sng_357[k];

        t_496[k] = f_16 * smg_358[k]
                   + f_3 * pc_x[k] * sng_358[k];

        t_497[k] = f_16 * smg_359[k]
                   + f_3 * pc_x[k] * sng_359[k];

        t_498[k] = f_16 * smg_265[k]
                   + f_1 * snf0_236[k]
                   - f_2 * snf1_236[k]
                   + f_3 * pc_y[k] * sng_355[k];
    }

#pragma omp simd aligned(t_499, t_500, t_501, pc_y, pc_z, smg_250, smg_267, smg_268, snf0_238, \
                         snf0_239, snf1_238, snf1_239, sng_355, sng_357, \
                         sng_358 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_499[k] = f_10 * smg_250[k]
                   + f_3 * pc_z[k] * sng_355[k];

        t_500[k] = f_16 * smg_267[k]
                   + f_4 * snf0_238[k]
                   - f_5 * snf1_238[k]
                   + f_3 * pc_y[k] * sng_357[k];

        t_501[k] = f_16 * smg_268[k]
                   + f_6 * snf0_239[k]
                   - f_7 * snf1_239[k]
                   + f_3 * pc_y[k] * sng_358[k];
    }

#pragma omp simd aligned(t_502, t_503, t_504, pc_x, pc_y, pc_z, smg_254, smg_269, smg_360, \
                         snf0_239, snf0_240, snf1_239, snf1_240, sng_359, \
                         sng_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_502[k] = f_16 * smg_269[k]
                   + f_3 * pc_y[k] * sng_359[k];

        t_503[k] = f_10 * smg_254[k]
                   + f_1 * snf0_239[k]
                   - f_2 * snf1_239[k]
                   + f_3 * pc_z[k] * sng_359[k];

        t_504[k] = f_16 * smg_360[k]
                   + f_1 * snf0_240[k]
                   - f_2 * snf1_240[k]
                   + f_3 * pc_x[k] * sng_360[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, pc_x, pc_y, pc_z, smg_255, smg_270, \
                         smg_272, smg_363, snf0_243, snf1_243, sng_360, sng_362, \
                         sng_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = f_11 * smg_270[k]
                   + f_3 * pc_y[k] * sng_360[k];

        t_506[k] = f_11 * smg_255[k]
                   + f_3 * pc_z[k] * sng_360[k];

        t_507[k] = f_16 * smg_363[k]
                   + f_4 * snf0_243[k]
                   - f_5 * snf1_243[k]
                   + f_3 * pc_x[k] * sng_363[k];

        t_508[k] = f_11 * smg_272[k]
                   + f_3 * pc_y[k] * sng_362[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, pc_x, pc_z, smg_258, smg_365, smg_366, snf0_245, \
                         snf0_246, snf1_245, snf1_246, sng_363, sng_365, \
                         sng_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = f_16 * smg_365[k]
                   + f_4 * snf0_245[k]
                   - f_5 * snf1_245[k]
                   + f_3 * pc_x[k] * sng_365[k];

        t_510[k] = f_16 * smg_366[k]
                   + f_6 * snf0_246[k]
                   - f_7 * snf1_246[k]
                   + f_3 * pc_x[k] * sng_366[k];

        t_511[k] = f_11 * smg_258[k]
                   + f_3 * pc_z[k] * sng_363[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, t_515, pc_x, pc_y, smg_275, smg_369, smg_370, \
                         smg_371, snf0_249, snf1_249, sng_365, sng_369, sng_370, \
                         sng_371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = f_11 * smg_275[k]
                   + f_3 * pc_y[k] * sng_365[k];

        t_513[k] = f_16 * smg_369[k]
                   + f_6 * snf0_249[k]
                   - f_7 * snf1_249[k]
                   + f_3 * pc_x[k] * sng_369[k];

        t_514[k] = f_16 * smg_370[k]
                   + f_3 * pc_x[k] * sng_370[k];

        t_515[k] = f_16 * smg_371[k]
                   + f_3 * pc_x[k] * sng_371[k];
    }

#pragma omp simd aligned(t_516, t_517, t_518, t_519, pc_x, pc_y, smg_280, smg_372, smg_373, \
                         smg_374, snf0_246, snf1_246, sng_370, sng_372, sng_373, \
                         sng_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_516[k] = f_16 * smg_372[k]
                   + f_3 * pc_x[k] * sng_372[k];

        t_517[k] = f_16 * smg_373[k]
                   + f_3 * pc_x[k] * sng_373[k];

        t_518[k] = f_16 * smg_374[k]
                   + f_3 * pc_x[k] * sng_374[k];

        t_519[k] = f_11 * smg_280[k]
                   + f_1 * snf0_246[k]
                   - f_2 * snf1_246[k]
                   + f_3 * pc_y[k] * sng_370[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, pc_y, pc_z, smg_265, smg_282, smg_283, snf0_248, \
                         snf0_249, snf1_248, snf1_249, sng_370, sng_372, \
                         sng_373 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = f_11 * smg_265[k]
                   + f_3 * pc_z[k] * sng_370[k];

        t_521[k] = f_11 * smg_282[k]
                   + f_4 * snf0_248[k]
                   - f_5 * snf1_248[k]
                   + f_3 * pc_y[k] * sng_372[k];

        t_522[k] = f_11 * smg_283[k]
                   + f_6 * snf0_249[k]
                   - f_7 * snf1_249[k]
                   + f_3 * pc_y[k] * sng_373[k];
    }

#pragma omp simd aligned(t_523, t_524, t_525, pc_x, pc_y, pc_z, smg_269, smg_284, smg_375, \
                         snf0_249, snf0_250, snf1_249, snf1_250, sng_374, \
                         sng_375 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_523[k] = f_11 * smg_284[k]
                   + f_3 * pc_y[k] * sng_374[k];

        t_524[k] = f_11 * smg_269[k]
                   + f_1 * snf0_249[k]
                   - f_2 * snf1_249[k]
                   + f_3 * pc_z[k] * sng_374[k];

        t_525[k] = f_16 * smg_375[k]
                   + f_1 * snf0_250[k]
                   - f_2 * snf1_250[k]
                   + f_3 * pc_x[k] * sng_375[k];
    }

#pragma omp simd aligned(t_526, t_527, t_528, t_529, pc_x, pc_y, pc_z, smg_270, smg_285, \
                         smg_287, smg_378, snf0_253, snf1_253, sng_375, sng_377, \
                         sng_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_526[k] = f_10 * smg_285[k]
                   + f_3 * pc_y[k] * sng_375[k];

        t_527[k] = f_16 * smg_270[k]
                   + f_3 * pc_z[k] * sng_375[k];

        t_528[k] = f_16 * smg_378[k]
                   + f_4 * snf0_253[k]
                   - f_5 * snf1_253[k]
                   + f_3 * pc_x[k] * sng_378[k];

        t_529[k] = f_10 * smg_287[k]
                   + f_3 * pc_y[k] * sng_377[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, pc_x, pc_z, smg_273, smg_380, smg_381, snf0_255, \
                         snf0_256, snf1_255, snf1_256, sng_378, sng_380, \
                         sng_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = f_16 * smg_380[k]
                   + f_4 * snf0_255[k]
                   - f_5 * snf1_255[k]
                   + f_3 * pc_x[k] * sng_380[k];

        t_531[k] = f_16 * smg_381[k]
                   + f_6 * snf0_256[k]
                   - f_7 * snf1_256[k]
                   + f_3 * pc_x[k] * sng_381[k];

        t_532[k] = f_16 * smg_273[k]
                   + f_3 * pc_z[k] * sng_378[k];
    }

#pragma omp simd aligned(t_533, t_534, t_535, t_536, pc_x, pc_y, smg_290, smg_384, smg_385, \
                         smg_386, snf0_259, snf1_259, sng_380, sng_384, sng_385, \
                         sng_386 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_533[k] = f_10 * smg_290[k]
                   + f_3 * pc_y[k] * sng_380[k];

        t_534[k] = f_16 * smg_384[k]
                   + f_6 * snf0_259[k]
                   - f_7 * snf1_259[k]
                   + f_3 * pc_x[k] * sng_384[k];

        t_535[k] = f_16 * smg_385[k]
                   + f_3 * pc_x[k] * sng_385[k];

        t_536[k] = f_16 * smg_386[k]
                   + f_3 * pc_x[k] * sng_386[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, t_540, pc_x, pc_y, smg_295, smg_387, smg_388, \
                         smg_389, snf0_256, snf1_256, sng_385, sng_387, sng_388, \
                         sng_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = f_16 * smg_387[k]
                   + f_3 * pc_x[k] * sng_387[k];

        t_538[k] = f_16 * smg_388[k]
                   + f_3 * pc_x[k] * sng_388[k];

        t_539[k] = f_16 * smg_389[k]
                   + f_3 * pc_x[k] * sng_389[k];

        t_540[k] = f_10 * smg_295[k]
                   + f_1 * snf0_256[k]
                   - f_2 * snf1_256[k]
                   + f_3 * pc_y[k] * sng_385[k];
    }

#pragma omp simd aligned(t_541, t_542, t_543, pc_y, pc_z, smg_280, smg_297, smg_298, snf0_258, \
                         snf0_259, snf1_258, snf1_259, sng_385, sng_387, \
                         sng_388 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_541[k] = f_16 * smg_280[k]
                   + f_3 * pc_z[k] * sng_385[k];

        t_542[k] = f_10 * smg_297[k]
                   + f_4 * snf0_258[k]
                   - f_5 * snf1_258[k]
                   + f_3 * pc_y[k] * sng_387[k];

        t_543[k] = f_10 * smg_298[k]
                   + f_6 * snf0_259[k]
                   - f_7 * snf1_259[k]
                   + f_3 * pc_y[k] * sng_388[k];
    }

#pragma omp simd aligned(t_544, t_545, t_546, t_547, pb_y, pc_y, pc_z, smh0_420, smg_284, \
                         smg_299, smg_300, smh1_420, snf0_259, snf1_259, sng_389, \
                         sng_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_544[k] = f_10 * smg_299[k]
                   + f_3 * pc_y[k] * sng_389[k];

        t_545[k] = f_16 * smg_284[k]
                   + f_1 * snf0_259[k]
                   - f_2 * snf1_259[k]
                   + f_3 * pc_z[k] * sng_389[k];

        t_546[k] = pb_y[k] * smh0_420[k]
                   - f_8 * pc_y[k] * smh1_420[k];

        t_547[k] = f_9 * smg_300[k]
                   + f_3 * pc_y[k] * sng_390[k];
    }

#pragma omp simd aligned(t_548, t_549, t_550, t_551, pb_y, pc_y, pc_z, smh0_423, smh0_425, \
                         smg_285, smg_301, smg_302, smh1_423, smh1_425, sng_390, \
                         sng_392 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_548[k] = f_17 * smg_285[k]
                   + f_3 * pc_z[k] * sng_390[k];

        t_549[k] = pb_y[k] * smh0_423[k]
                   + f_10 * smg_301[k]
                   - f_8 * pc_y[k] * smh1_423[k];

        t_550[k] = f_9 * smg_302[k]
                   + f_3 * pc_y[k] * sng_392[k];

        t_551[k] = pb_y[k] * smh0_425[k]
                   - f_8 * pc_y[k] * smh1_425[k];
    }

#pragma omp simd aligned(t_552, t_553, t_554, t_555, pb_y, pc_y, pc_z, smh0_426, smh0_429, \
                         smg_288, smg_303, smg_305, smh1_426, smh1_429, sng_393, \
                         sng_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_552[k] = pb_y[k] * smh0_426[k]
                   + f_11 * smg_303[k]
                   - f_8 * pc_y[k] * smh1_426[k];

        t_553[k] = f_17 * smg_288[k]
                   + f_3 * pc_z[k] * sng_393[k];

        t_554[k] = f_9 * smg_305[k]
                   + f_3 * pc_y[k] * sng_395[k];

        t_555[k] = pb_y[k] * smh0_429[k]
                   - f_8 * pc_y[k] * smh1_429[k];
    }

#pragma omp simd aligned(t_556, t_557, t_558, t_559, t_560, pc_x, smg_400, smg_401, smg_402, \
                         smg_403, smg_404, sng_400, sng_401, sng_402, sng_403, \
                         sng_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_556[k] = f_16 * smg_400[k]
                   + f_3 * pc_x[k] * sng_400[k];

        t_557[k] = f_16 * smg_401[k]
                   + f_3 * pc_x[k] * sng_401[k];

        t_558[k] = f_16 * smg_402[k]
                   + f_3 * pc_x[k] * sng_402[k];

        t_559[k] = f_16 * smg_403[k]
                   + f_3 * pc_x[k] * sng_403[k];

        t_560[k] = f_16 * smg_404[k]
                   + f_3 * pc_x[k] * sng_404[k];
    }

#pragma omp simd aligned(t_561, t_562, t_563, pc_y, pc_z, smg_295, smg_310, smg_312, snf0_266, \
                         snf0_268, snf1_266, snf1_268, sng_400, \
                         sng_402 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_561[k] = f_9 * smg_310[k]
                   + f_1 * snf0_266[k]
                   - f_2 * snf1_266[k]
                   + f_3 * pc_y[k] * sng_400[k];

        t_562[k] = f_17 * smg_295[k]
                   + f_3 * pc_z[k] * sng_400[k];

        t_563[k] = f_9 * smg_312[k]
                   + f_4 * snf0_268[k]
                   - f_5 * snf1_268[k]
                   + f_3 * pc_y[k] * sng_402[k];
    }

#pragma omp simd aligned(t_564, t_565, t_566, pb_y, pc_y, smh0_440, smg_313, smg_314, \
                         smh1_440, snf0_269, snf1_269, sng_403, \
                         sng_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_564[k] = f_9 * smg_313[k]
                   + f_6 * snf0_269[k]
                   - f_7 * snf1_269[k]
                   + f_3 * pc_y[k] * sng_403[k];

        t_565[k] = f_9 * smg_314[k]
                   + f_3 * pc_y[k] * sng_404[k];

        t_566[k] = pb_y[k] * smh0_440[k]
                   - f_8 * pc_y[k] * smh1_440[k];
    }

#pragma omp simd aligned(t_567, t_568, t_569, t_570, pc_x, pc_y, pc_z, smg_300, smg_405, \
                         smg_408, snf0_270, snf0_273, snf1_270, snf1_273, sng_405, \
                         sng_408 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_567[k] = f_16 * smg_405[k]
                   + f_1 * snf0_270[k]
                   - f_2 * snf1_270[k]
                   + f_3 * pc_x[k] * sng_405[k];

        t_568[k] = f_3 * pc_y[k] * sng_405[k];

        t_569[k] = f_15 * smg_300[k]
                   + f_3 * pc_z[k] * sng_405[k];

        t_570[k] = f_16 * smg_408[k]
                   + f_4 * snf0_273[k]
                   - f_5 * snf1_273[k]
                   + f_3 * pc_x[k] * sng_408[k];
    }

#pragma omp simd aligned(t_571, t_572, t_573, pc_x, pc_y, smg_410, smg_411, snf0_275, \
                         snf0_276, snf1_275, snf1_276, sng_407, sng_410, \
                         sng_411 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_571[k] = f_3 * pc_y[k] * sng_407[k];

        t_572[k] = f_16 * smg_410[k]
                   + f_4 * snf0_275[k]
                   - f_5 * snf1_275[k]
                   + f_3 * pc_x[k] * sng_410[k];

        t_573[k] = f_16 * smg_411[k]
                   + f_6 * snf0_276[k]
                   - f_7 * snf1_276[k]
                   + f_3 * pc_x[k] * sng_411[k];
    }

#pragma omp simd aligned(t_574, t_575, t_576, t_577, pc_x, pc_y, pc_z, smg_303, smg_414, \
                         smg_415, snf0_279, snf1_279, sng_408, sng_410, sng_414, \
                         sng_415 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_574[k] = f_15 * smg_303[k]
                   + f_3 * pc_z[k] * sng_408[k];

        t_575[k] = f_3 * pc_y[k] * sng_410[k];

        t_576[k] = f_16 * smg_414[k]
                   + f_6 * snf0_279[k]
                   - f_7 * snf1_279[k]
                   + f_3 * pc_x[k] * sng_414[k];

        t_577[k] = f_16 * smg_415[k]
                   + f_3 * pc_x[k] * sng_415[k];
    }

#pragma omp simd aligned(t_578, t_579, t_580, t_581, pc_x, smg_416, smg_417, smg_418, smg_419, \
                         sng_416, sng_417, sng_418, sng_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_578[k] = f_16 * smg_416[k]
                   + f_3 * pc_x[k] * sng_416[k];

        t_579[k] = f_16 * smg_417[k]
                   + f_3 * pc_x[k] * sng_417[k];

        t_580[k] = f_16 * smg_418[k]
                   + f_3 * pc_x[k] * sng_418[k];

        t_581[k] = f_16 * smg_419[k]
                   + f_3 * pc_x[k] * sng_419[k];
    }

#pragma omp simd aligned(t_582, t_583, t_584, t_585, pc_y, pc_z, smg_310, snf0_276, snf0_278, \
                         snf0_279, snf1_276, snf1_278, snf1_279, sng_415, sng_417, \
                         sng_418 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_582[k] = f_1 * snf0_276[k]
                   - f_2 * snf1_276[k]
                   + f_3 * pc_y[k] * sng_415[k];

        t_583[k] = f_15 * smg_310[k]
                   + f_3 * pc_z[k] * sng_415[k];

        t_584[k] = f_4 * snf0_278[k]
                   - f_5 * snf1_278[k]
                   + f_3 * pc_y[k] * sng_417[k];

        t_585[k] = f_6 * snf0_279[k]
                   - f_7 * snf1_279[k]
                   + f_3 * pc_y[k] * sng_418[k];
    }
}

static auto
compute_prim_snh_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t smh0,
                                                          const size_t smg, const size_t smh1,
                                                          const size_t snf0, const size_t snf1,
                                                          const size_t sng, const size_t ncols,
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
    const auto f_14 = 3.5 / q;
    const auto f_15 = 3.0 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;

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

    const auto *smh0_441 = buffer.data(smh0 + 441);
    const auto *smh0_444 = buffer.data(smh0 + 444);
    const auto *smh0_447 = buffer.data(smh0 + 447);
    const auto *smh0_456 = buffer.data(smh0 + 456);

    const auto *smg_314 = buffer.data(smg + 314);
    const auto *smg_315 = buffer.data(smg + 315);
    const auto *smg_317 = buffer.data(smg + 317);
    const auto *smg_318 = buffer.data(smg + 318);
    const auto *smg_320 = buffer.data(smg + 320);
    const auto *smg_325 = buffer.data(smg + 325);
    const auto *smg_327 = buffer.data(smg + 327);
    const auto *smg_328 = buffer.data(smg + 328);
    const auto *smg_329 = buffer.data(smg + 329);
    const auto *smg_330 = buffer.data(smg + 330);
    const auto *smg_332 = buffer.data(smg + 332);
    const auto *smg_333 = buffer.data(smg + 333);
    const auto *smg_335 = buffer.data(smg + 335);
    const auto *smg_340 = buffer.data(smg + 340);
    const auto *smg_342 = buffer.data(smg + 342);
    const auto *smg_343 = buffer.data(smg + 343);
    const auto *smg_344 = buffer.data(smg + 344);
    const auto *smg_345 = buffer.data(smg + 345);
    const auto *smg_347 = buffer.data(smg + 347);
    const auto *smg_348 = buffer.data(smg + 348);
    const auto *smg_350 = buffer.data(smg + 350);
    const auto *smg_355 = buffer.data(smg + 355);
    const auto *smg_357 = buffer.data(smg + 357);
    const auto *smg_358 = buffer.data(smg + 358);
    const auto *smg_359 = buffer.data(smg + 359);
    const auto *smg_360 = buffer.data(smg + 360);
    const auto *smg_362 = buffer.data(smg + 362);
    const auto *smg_363 = buffer.data(smg + 363);
    const auto *smg_365 = buffer.data(smg + 365);
    const auto *smg_370 = buffer.data(smg + 370);
    const auto *smg_372 = buffer.data(smg + 372);
    const auto *smg_373 = buffer.data(smg + 373);
    const auto *smg_374 = buffer.data(smg + 374);
    const auto *smg_375 = buffer.data(smg + 375);
    const auto *smg_377 = buffer.data(smg + 377);
    const auto *smg_380 = buffer.data(smg + 380);
    const auto *smg_385 = buffer.data(smg + 385);
    const auto *smg_387 = buffer.data(smg + 387);
    const auto *smg_388 = buffer.data(smg + 388);
    const auto *smg_389 = buffer.data(smg + 389);
    const auto *smg_420 = buffer.data(smg + 420);
    const auto *smg_423 = buffer.data(smg + 423);
    const auto *smg_425 = buffer.data(smg + 425);
    const auto *smg_426 = buffer.data(smg + 426);
    const auto *smg_429 = buffer.data(smg + 429);
    const auto *smg_430 = buffer.data(smg + 430);
    const auto *smg_431 = buffer.data(smg + 431);
    const auto *smg_432 = buffer.data(smg + 432);
    const auto *smg_433 = buffer.data(smg + 433);
    const auto *smg_434 = buffer.data(smg + 434);
    const auto *smg_440 = buffer.data(smg + 440);
    const auto *smg_444 = buffer.data(smg + 444);
    const auto *smg_445 = buffer.data(smg + 445);
    const auto *smg_446 = buffer.data(smg + 446);
    const auto *smg_447 = buffer.data(smg + 447);
    const auto *smg_448 = buffer.data(smg + 448);
    const auto *smg_449 = buffer.data(smg + 449);
    const auto *smg_450 = buffer.data(smg + 450);
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

    const auto *smh1_441 = buffer.data(smh1 + 441);
    const auto *smh1_444 = buffer.data(smh1 + 444);
    const auto *smh1_447 = buffer.data(smh1 + 447);
    const auto *smh1_456 = buffer.data(smh1 + 456);

    const auto *snf0_279 = buffer.data(snf0 + 279);
    const auto *snf0_280 = buffer.data(snf0 + 280);
    const auto *snf0_283 = buffer.data(snf0 + 283);
    const auto *snf0_285 = buffer.data(snf0 + 285);
    const auto *snf0_286 = buffer.data(snf0 + 286);
    const auto *snf0_288 = buffer.data(snf0 + 288);
    const auto *snf0_289 = buffer.data(snf0 + 289);
    const auto *snf0_295 = buffer.data(snf0 + 295);
    const auto *snf0_298 = buffer.data(snf0 + 298);
    const auto *snf0_299 = buffer.data(snf0 + 299);
    const auto *snf0_300 = buffer.data(snf0 + 300);
    const auto *snf0_303 = buffer.data(snf0 + 303);
    const auto *snf0_305 = buffer.data(snf0 + 305);
    const auto *snf0_306 = buffer.data(snf0 + 306);
    const auto *snf0_308 = buffer.data(snf0 + 308);
    const auto *snf0_309 = buffer.data(snf0 + 309);
    const auto *snf0_310 = buffer.data(snf0 + 310);
    const auto *snf0_313 = buffer.data(snf0 + 313);
    const auto *snf0_315 = buffer.data(snf0 + 315);
    const auto *snf0_316 = buffer.data(snf0 + 316);
    const auto *snf0_318 = buffer.data(snf0 + 318);
    const auto *snf0_319 = buffer.data(snf0 + 319);
    const auto *snf0_320 = buffer.data(snf0 + 320);
    const auto *snf0_323 = buffer.data(snf0 + 323);
    const auto *snf0_325 = buffer.data(snf0 + 325);
    const auto *snf0_326 = buffer.data(snf0 + 326);
    const auto *snf0_328 = buffer.data(snf0 + 328);
    const auto *snf0_329 = buffer.data(snf0 + 329);
    const auto *snf0_330 = buffer.data(snf0 + 330);

    const auto *snf1_279 = buffer.data(snf1 + 279);
    const auto *snf1_280 = buffer.data(snf1 + 280);
    const auto *snf1_283 = buffer.data(snf1 + 283);
    const auto *snf1_285 = buffer.data(snf1 + 285);
    const auto *snf1_286 = buffer.data(snf1 + 286);
    const auto *snf1_288 = buffer.data(snf1 + 288);
    const auto *snf1_289 = buffer.data(snf1 + 289);
    const auto *snf1_295 = buffer.data(snf1 + 295);
    const auto *snf1_298 = buffer.data(snf1 + 298);
    const auto *snf1_299 = buffer.data(snf1 + 299);
    const auto *snf1_300 = buffer.data(snf1 + 300);
    const auto *snf1_303 = buffer.data(snf1 + 303);
    const auto *snf1_305 = buffer.data(snf1 + 305);
    const auto *snf1_306 = buffer.data(snf1 + 306);
    const auto *snf1_308 = buffer.data(snf1 + 308);
    const auto *snf1_309 = buffer.data(snf1 + 309);
    const auto *snf1_310 = buffer.data(snf1 + 310);
    const auto *snf1_313 = buffer.data(snf1 + 313);
    const auto *snf1_315 = buffer.data(snf1 + 315);
    const auto *snf1_316 = buffer.data(snf1 + 316);
    const auto *snf1_318 = buffer.data(snf1 + 318);
    const auto *snf1_319 = buffer.data(snf1 + 319);
    const auto *snf1_320 = buffer.data(snf1 + 320);
    const auto *snf1_323 = buffer.data(snf1 + 323);
    const auto *snf1_325 = buffer.data(snf1 + 325);
    const auto *snf1_326 = buffer.data(snf1 + 326);
    const auto *snf1_328 = buffer.data(snf1 + 328);
    const auto *snf1_329 = buffer.data(snf1 + 329);
    const auto *snf1_330 = buffer.data(snf1 + 330);

    const auto *sng_419 = buffer.data(sng + 419);
    const auto *sng_420 = buffer.data(sng + 420);
    const auto *sng_422 = buffer.data(sng + 422);
    const auto *sng_423 = buffer.data(sng + 423);
    const auto *sng_425 = buffer.data(sng + 425);
    const auto *sng_426 = buffer.data(sng + 426);
    const auto *sng_429 = buffer.data(sng + 429);
    const auto *sng_430 = buffer.data(sng + 430);
    const auto *sng_431 = buffer.data(sng + 431);
    const auto *sng_432 = buffer.data(sng + 432);
    const auto *sng_433 = buffer.data(sng + 433);
    const auto *sng_434 = buffer.data(sng + 434);
    const auto *sng_435 = buffer.data(sng + 435);
    const auto *sng_437 = buffer.data(sng + 437);
    const auto *sng_438 = buffer.data(sng + 438);
    const auto *sng_440 = buffer.data(sng + 440);
    const auto *sng_444 = buffer.data(sng + 444);
    const auto *sng_445 = buffer.data(sng + 445);
    const auto *sng_446 = buffer.data(sng + 446);
    const auto *sng_447 = buffer.data(sng + 447);
    const auto *sng_448 = buffer.data(sng + 448);
    const auto *sng_449 = buffer.data(sng + 449);
    const auto *sng_450 = buffer.data(sng + 450);
    const auto *sng_452 = buffer.data(sng + 452);
    const auto *sng_453 = buffer.data(sng + 453);
    const auto *sng_455 = buffer.data(sng + 455);
    const auto *sng_456 = buffer.data(sng + 456);
    const auto *sng_459 = buffer.data(sng + 459);
    const auto *sng_460 = buffer.data(sng + 460);
    const auto *sng_461 = buffer.data(sng + 461);
    const auto *sng_462 = buffer.data(sng + 462);
    const auto *sng_463 = buffer.data(sng + 463);
    const auto *sng_464 = buffer.data(sng + 464);
    const auto *sng_465 = buffer.data(sng + 465);
    const auto *sng_467 = buffer.data(sng + 467);
    const auto *sng_468 = buffer.data(sng + 468);
    const auto *sng_470 = buffer.data(sng + 470);
    const auto *sng_471 = buffer.data(sng + 471);
    const auto *sng_474 = buffer.data(sng + 474);
    const auto *sng_475 = buffer.data(sng + 475);
    const auto *sng_476 = buffer.data(sng + 476);
    const auto *sng_477 = buffer.data(sng + 477);
    const auto *sng_478 = buffer.data(sng + 478);
    const auto *sng_479 = buffer.data(sng + 479);
    const auto *sng_480 = buffer.data(sng + 480);
    const auto *sng_482 = buffer.data(sng + 482);
    const auto *sng_483 = buffer.data(sng + 483);
    const auto *sng_485 = buffer.data(sng + 485);
    const auto *sng_486 = buffer.data(sng + 486);
    const auto *sng_489 = buffer.data(sng + 489);
    const auto *sng_490 = buffer.data(sng + 490);
    const auto *sng_491 = buffer.data(sng + 491);
    const auto *sng_492 = buffer.data(sng + 492);
    const auto *sng_493 = buffer.data(sng + 493);
    const auto *sng_494 = buffer.data(sng + 494);
    const auto *sng_495 = buffer.data(sng + 495);

#pragma omp simd aligned(t_586, t_587, t_588, t_589, pc_x, pc_y, pc_z, smg_314, smg_315, \
                         smg_420, snf0_279, snf0_280, snf1_279, snf1_280, sng_419, \
                         sng_420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_586[k] = f_3 * pc_y[k] * sng_419[k];

        t_587[k] = f_15 * smg_314[k]
                   + f_1 * snf0_279[k]
                   - f_2 * snf1_279[k]
                   + f_3 * pc_z[k] * sng_419[k];

        t_588[k] = f_11 * smg_420[k]
                   + f_1 * snf0_280[k]
                   - f_2 * snf1_280[k]
                   + f_3 * pc_x[k] * sng_420[k];

        t_589[k] = f_14 * smg_315[k]
                   + f_3 * pc_y[k] * sng_420[k];
    }

#pragma omp simd aligned(t_590, t_591, t_592, pc_x, pc_y, pc_z, smg_317, smg_423, snf0_283, \
                         snf1_283, sng_420, sng_422, sng_423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_590[k] = f_3 * pc_z[k] * sng_420[k];

        t_591[k] = f_11 * smg_423[k]
                   + f_4 * snf0_283[k]
                   - f_5 * snf1_283[k]
                   + f_3 * pc_x[k] * sng_423[k];

        t_592[k] = f_14 * smg_317[k]
                   + f_3 * pc_y[k] * sng_422[k];
    }

#pragma omp simd aligned(t_593, t_594, t_595, pc_x, pc_z, smg_425, smg_426, snf0_285, \
                         snf0_286, snf1_285, snf1_286, sng_423, sng_425, \
                         sng_426 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_593[k] = f_11 * smg_425[k]
                   + f_4 * snf0_285[k]
                   - f_5 * snf1_285[k]
                   + f_3 * pc_x[k] * sng_425[k];

        t_594[k] = f_11 * smg_426[k]
                   + f_6 * snf0_286[k]
                   - f_7 * snf1_286[k]
                   + f_3 * pc_x[k] * sng_426[k];

        t_595[k] = f_3 * pc_z[k] * sng_423[k];
    }

#pragma omp simd aligned(t_596, t_597, t_598, t_599, pc_x, pc_y, smg_320, smg_429, smg_430, \
                         smg_431, snf0_289, snf1_289, sng_425, sng_429, sng_430, \
                         sng_431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_596[k] = f_14 * smg_320[k]
                   + f_3 * pc_y[k] * sng_425[k];

        t_597[k] = f_11 * smg_429[k]
                   + f_6 * snf0_289[k]
                   - f_7 * snf1_289[k]
                   + f_3 * pc_x[k] * sng_429[k];

        t_598[k] = f_11 * smg_430[k]
                   + f_3 * pc_x[k] * sng_430[k];

        t_599[k] = f_11 * smg_431[k]
                   + f_3 * pc_x[k] * sng_431[k];
    }

#pragma omp simd aligned(t_600, t_601, t_602, t_603, pc_x, pc_y, smg_325, smg_432, smg_433, \
                         smg_434, snf0_286, snf1_286, sng_430, sng_432, sng_433, \
                         sng_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = f_11 * smg_432[k]
                   + f_3 * pc_x[k] * sng_432[k];

        t_601[k] = f_11 * smg_433[k]
                   + f_3 * pc_x[k] * sng_433[k];

        t_602[k] = f_11 * smg_434[k]
                   + f_3 * pc_x[k] * sng_434[k];

        t_603[k] = f_14 * smg_325[k]
                   + f_1 * snf0_286[k]
                   - f_2 * snf1_286[k]
                   + f_3 * pc_y[k] * sng_430[k];
    }

#pragma omp simd aligned(t_604, t_605, t_606, pc_y, pc_z, smg_327, smg_328, snf0_288, \
                         snf0_289, snf1_288, snf1_289, sng_430, sng_432, \
                         sng_433 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_604[k] = f_3 * pc_z[k] * sng_430[k];

        t_605[k] = f_14 * smg_327[k]
                   + f_4 * snf0_288[k]
                   - f_5 * snf1_288[k]
                   + f_3 * pc_y[k] * sng_432[k];

        t_606[k] = f_14 * smg_328[k]
                   + f_6 * snf0_289[k]
                   - f_7 * snf1_289[k]
                   + f_3 * pc_y[k] * sng_433[k];
    }

#pragma omp simd aligned(t_607, t_608, t_609, t_610, pb_z, pc_y, pc_z, smh0_441, smg_329, \
                         smg_330, smh1_441, snf0_289, snf1_289, sng_434, \
                         sng_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_607[k] = f_14 * smg_329[k]
                   + f_3 * pc_y[k] * sng_434[k];

        t_608[k] = f_1 * snf0_289[k]
                   - f_2 * snf1_289[k]
                   + f_3 * pc_z[k] * sng_434[k];

        t_609[k] = pb_z[k] * smh0_441[k]
                   - f_8 * pc_z[k] * smh1_441[k];

        t_610[k] = f_15 * smg_330[k]
                   + f_3 * pc_y[k] * sng_435[k];
    }

#pragma omp simd aligned(t_611, t_612, t_613, pb_z, pc_y, pc_z, smh0_444, smg_315, smg_332, \
                         smh1_444, sng_435, sng_437 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_611[k] = f_9 * smg_315[k]
                   + f_3 * pc_z[k] * sng_435[k];

        t_612[k] = pb_z[k] * smh0_444[k]
                   - f_8 * pc_z[k] * smh1_444[k];

        t_613[k] = f_15 * smg_332[k]
                   + f_3 * pc_y[k] * sng_437[k];
    }

#pragma omp simd aligned(t_614, t_615, t_616, pb_z, pc_x, pc_z, smh0_447, smg_318, smg_440, \
                         smh1_447, snf0_295, snf1_295, sng_438, \
                         sng_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_614[k] = f_11 * smg_440[k]
                   + f_4 * snf0_295[k]
                   - f_5 * snf1_295[k]
                   + f_3 * pc_x[k] * sng_440[k];

        t_615[k] = pb_z[k] * smh0_447[k]
                   - f_8 * pc_z[k] * smh1_447[k];

        t_616[k] = f_9 * smg_318[k]
                   + f_3 * pc_z[k] * sng_438[k];
    }

#pragma omp simd aligned(t_617, t_618, t_619, t_620, pc_x, pc_y, smg_335, smg_444, smg_445, \
                         smg_446, snf0_299, snf1_299, sng_440, sng_444, sng_445, \
                         sng_446 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_617[k] = f_15 * smg_335[k]
                   + f_3 * pc_y[k] * sng_440[k];

        t_618[k] = f_11 * smg_444[k]
                   + f_6 * snf0_299[k]
                   - f_7 * snf1_299[k]
                   + f_3 * pc_x[k] * sng_444[k];

        t_619[k] = f_11 * smg_445[k]
                   + f_3 * pc_x[k] * sng_445[k];

        t_620[k] = f_11 * smg_446[k]
                   + f_3 * pc_x[k] * sng_446[k];
    }

#pragma omp simd aligned(t_621, t_622, t_623, t_624, pb_z, pc_x, pc_z, smh0_456, smg_447, \
                         smg_448, smg_449, smh1_456, sng_447, sng_448, \
                         sng_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_621[k] = f_11 * smg_447[k]
                   + f_3 * pc_x[k] * sng_447[k];

        t_622[k] = f_11 * smg_448[k]
                   + f_3 * pc_x[k] * sng_448[k];

        t_623[k] = f_11 * smg_449[k]
                   + f_3 * pc_x[k] * sng_449[k];

        t_624[k] = pb_z[k] * smh0_456[k]
                   - f_8 * pc_z[k] * smh1_456[k];
    }

#pragma omp simd aligned(t_625, t_626, t_627, pc_y, pc_z, smg_325, smg_342, smg_343, snf0_298, \
                         snf0_299, snf1_298, snf1_299, sng_445, sng_447, \
                         sng_448 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_625[k] = f_9 * smg_325[k]
                   + f_3 * pc_z[k] * sng_445[k];

        t_626[k] = f_15 * smg_342[k]
                   + f_4 * snf0_298[k]
                   - f_5 * snf1_298[k]
                   + f_3 * pc_y[k] * sng_447[k];

        t_627[k] = f_15 * smg_343[k]
                   + f_6 * snf0_299[k]
                   - f_7 * snf1_299[k]
                   + f_3 * pc_y[k] * sng_448[k];
    }

#pragma omp simd aligned(t_628, t_629, t_630, pc_x, pc_y, pc_z, smg_329, smg_344, smg_450, \
                         snf0_299, snf0_300, snf1_299, snf1_300, sng_449, \
                         sng_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_628[k] = f_15 * smg_344[k]
                   + f_3 * pc_y[k] * sng_449[k];

        t_629[k] = f_9 * smg_329[k]
                   + f_1 * snf0_299[k]
                   - f_2 * snf1_299[k]
                   + f_3 * pc_z[k] * sng_449[k];

        t_630[k] = f_11 * smg_450[k]
                   + f_1 * snf0_300[k]
                   - f_2 * snf1_300[k]
                   + f_3 * pc_x[k] * sng_450[k];
    }

#pragma omp simd aligned(t_631, t_632, t_633, t_634, pc_x, pc_y, pc_z, smg_330, smg_345, \
                         smg_347, smg_453, snf0_303, snf1_303, sng_450, sng_452, \
                         sng_453 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_631[k] = f_17 * smg_345[k]
                   + f_3 * pc_y[k] * sng_450[k];

        t_632[k] = f_10 * smg_330[k]
                   + f_3 * pc_z[k] * sng_450[k];

        t_633[k] = f_11 * smg_453[k]
                   + f_4 * snf0_303[k]
                   - f_5 * snf1_303[k]
                   + f_3 * pc_x[k] * sng_453[k];

        t_634[k] = f_17 * smg_347[k]
                   + f_3 * pc_y[k] * sng_452[k];
    }

#pragma omp simd aligned(t_635, t_636, t_637, pc_x, pc_z, smg_333, smg_455, smg_456, snf0_305, \
                         snf0_306, snf1_305, snf1_306, sng_453, sng_455, \
                         sng_456 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_635[k] = f_11 * smg_455[k]
                   + f_4 * snf0_305[k]
                   - f_5 * snf1_305[k]
                   + f_3 * pc_x[k] * sng_455[k];

        t_636[k] = f_11 * smg_456[k]
                   + f_6 * snf0_306[k]
                   - f_7 * snf1_306[k]
                   + f_3 * pc_x[k] * sng_456[k];

        t_637[k] = f_10 * smg_333[k]
                   + f_3 * pc_z[k] * sng_453[k];
    }

#pragma omp simd aligned(t_638, t_639, t_640, t_641, pc_x, pc_y, smg_350, smg_459, smg_460, \
                         smg_461, snf0_309, snf1_309, sng_455, sng_459, sng_460, \
                         sng_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_638[k] = f_17 * smg_350[k]
                   + f_3 * pc_y[k] * sng_455[k];

        t_639[k] = f_11 * smg_459[k]
                   + f_6 * snf0_309[k]
                   - f_7 * snf1_309[k]
                   + f_3 * pc_x[k] * sng_459[k];

        t_640[k] = f_11 * smg_460[k]
                   + f_3 * pc_x[k] * sng_460[k];

        t_641[k] = f_11 * smg_461[k]
                   + f_3 * pc_x[k] * sng_461[k];
    }

#pragma omp simd aligned(t_642, t_643, t_644, t_645, pc_x, pc_y, smg_355, smg_462, smg_463, \
                         smg_464, snf0_306, snf1_306, sng_460, sng_462, sng_463, \
                         sng_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_642[k] = f_11 * smg_462[k]
                   + f_3 * pc_x[k] * sng_462[k];

        t_643[k] = f_11 * smg_463[k]
                   + f_3 * pc_x[k] * sng_463[k];

        t_644[k] = f_11 * smg_464[k]
                   + f_3 * pc_x[k] * sng_464[k];

        t_645[k] = f_17 * smg_355[k]
                   + f_1 * snf0_306[k]
                   - f_2 * snf1_306[k]
                   + f_3 * pc_y[k] * sng_460[k];
    }

#pragma omp simd aligned(t_646, t_647, t_648, pc_y, pc_z, smg_340, smg_357, smg_358, snf0_308, \
                         snf0_309, snf1_308, snf1_309, sng_460, sng_462, \
                         sng_463 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_646[k] = f_10 * smg_340[k]
                   + f_3 * pc_z[k] * sng_460[k];

        t_647[k] = f_17 * smg_357[k]
                   + f_4 * snf0_308[k]
                   - f_5 * snf1_308[k]
                   + f_3 * pc_y[k] * sng_462[k];

        t_648[k] = f_17 * smg_358[k]
                   + f_6 * snf0_309[k]
                   - f_7 * snf1_309[k]
                   + f_3 * pc_y[k] * sng_463[k];
    }

#pragma omp simd aligned(t_649, t_650, t_651, pc_x, pc_y, pc_z, smg_344, smg_359, smg_465, \
                         snf0_309, snf0_310, snf1_309, snf1_310, sng_464, \
                         sng_465 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_649[k] = f_17 * smg_359[k]
                   + f_3 * pc_y[k] * sng_464[k];

        t_650[k] = f_10 * smg_344[k]
                   + f_1 * snf0_309[k]
                   - f_2 * snf1_309[k]
                   + f_3 * pc_z[k] * sng_464[k];

        t_651[k] = f_11 * smg_465[k]
                   + f_1 * snf0_310[k]
                   - f_2 * snf1_310[k]
                   + f_3 * pc_x[k] * sng_465[k];
    }

#pragma omp simd aligned(t_652, t_653, t_654, t_655, pc_x, pc_y, pc_z, smg_345, smg_360, \
                         smg_362, smg_468, snf0_313, snf1_313, sng_465, sng_467, \
                         sng_468 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_652[k] = f_16 * smg_360[k]
                   + f_3 * pc_y[k] * sng_465[k];

        t_653[k] = f_11 * smg_345[k]
                   + f_3 * pc_z[k] * sng_465[k];

        t_654[k] = f_11 * smg_468[k]
                   + f_4 * snf0_313[k]
                   - f_5 * snf1_313[k]
                   + f_3 * pc_x[k] * sng_468[k];

        t_655[k] = f_16 * smg_362[k]
                   + f_3 * pc_y[k] * sng_467[k];
    }

#pragma omp simd aligned(t_656, t_657, t_658, pc_x, pc_z, smg_348, smg_470, smg_471, snf0_315, \
                         snf0_316, snf1_315, snf1_316, sng_468, sng_470, \
                         sng_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_656[k] = f_11 * smg_470[k]
                   + f_4 * snf0_315[k]
                   - f_5 * snf1_315[k]
                   + f_3 * pc_x[k] * sng_470[k];

        t_657[k] = f_11 * smg_471[k]
                   + f_6 * snf0_316[k]
                   - f_7 * snf1_316[k]
                   + f_3 * pc_x[k] * sng_471[k];

        t_658[k] = f_11 * smg_348[k]
                   + f_3 * pc_z[k] * sng_468[k];
    }

#pragma omp simd aligned(t_659, t_660, t_661, t_662, pc_x, pc_y, smg_365, smg_474, smg_475, \
                         smg_476, snf0_319, snf1_319, sng_470, sng_474, sng_475, \
                         sng_476 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_659[k] = f_16 * smg_365[k]
                   + f_3 * pc_y[k] * sng_470[k];

        t_660[k] = f_11 * smg_474[k]
                   + f_6 * snf0_319[k]
                   - f_7 * snf1_319[k]
                   + f_3 * pc_x[k] * sng_474[k];

        t_661[k] = f_11 * smg_475[k]
                   + f_3 * pc_x[k] * sng_475[k];

        t_662[k] = f_11 * smg_476[k]
                   + f_3 * pc_x[k] * sng_476[k];
    }

#pragma omp simd aligned(t_663, t_664, t_665, t_666, pc_x, pc_y, smg_370, smg_477, smg_478, \
                         smg_479, snf0_316, snf1_316, sng_475, sng_477, sng_478, \
                         sng_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_663[k] = f_11 * smg_477[k]
                   + f_3 * pc_x[k] * sng_477[k];

        t_664[k] = f_11 * smg_478[k]
                   + f_3 * pc_x[k] * sng_478[k];

        t_665[k] = f_11 * smg_479[k]
                   + f_3 * pc_x[k] * sng_479[k];

        t_666[k] = f_16 * smg_370[k]
                   + f_1 * snf0_316[k]
                   - f_2 * snf1_316[k]
                   + f_3 * pc_y[k] * sng_475[k];
    }

#pragma omp simd aligned(t_667, t_668, t_669, pc_y, pc_z, smg_355, smg_372, smg_373, snf0_318, \
                         snf0_319, snf1_318, snf1_319, sng_475, sng_477, \
                         sng_478 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_667[k] = f_11 * smg_355[k]
                   + f_3 * pc_z[k] * sng_475[k];

        t_668[k] = f_16 * smg_372[k]
                   + f_4 * snf0_318[k]
                   - f_5 * snf1_318[k]
                   + f_3 * pc_y[k] * sng_477[k];

        t_669[k] = f_16 * smg_373[k]
                   + f_6 * snf0_319[k]
                   - f_7 * snf1_319[k]
                   + f_3 * pc_y[k] * sng_478[k];
    }

#pragma omp simd aligned(t_670, t_671, t_672, pc_x, pc_y, pc_z, smg_359, smg_374, smg_480, \
                         snf0_319, snf0_320, snf1_319, snf1_320, sng_479, \
                         sng_480 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = f_16 * smg_374[k]
                   + f_3 * pc_y[k] * sng_479[k];

        t_671[k] = f_11 * smg_359[k]
                   + f_1 * snf0_319[k]
                   - f_2 * snf1_319[k]
                   + f_3 * pc_z[k] * sng_479[k];

        t_672[k] = f_11 * smg_480[k]
                   + f_1 * snf0_320[k]
                   - f_2 * snf1_320[k]
                   + f_3 * pc_x[k] * sng_480[k];
    }

#pragma omp simd aligned(t_673, t_674, t_675, t_676, pc_x, pc_y, pc_z, smg_360, smg_375, \
                         smg_377, smg_483, snf0_323, snf1_323, sng_480, sng_482, \
                         sng_483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_673[k] = f_11 * smg_375[k]
                   + f_3 * pc_y[k] * sng_480[k];

        t_674[k] = f_16 * smg_360[k]
                   + f_3 * pc_z[k] * sng_480[k];

        t_675[k] = f_11 * smg_483[k]
                   + f_4 * snf0_323[k]
                   - f_5 * snf1_323[k]
                   + f_3 * pc_x[k] * sng_483[k];

        t_676[k] = f_11 * smg_377[k]
                   + f_3 * pc_y[k] * sng_482[k];
    }

#pragma omp simd aligned(t_677, t_678, t_679, pc_x, pc_z, smg_363, smg_485, smg_486, snf0_325, \
                         snf0_326, snf1_325, snf1_326, sng_483, sng_485, \
                         sng_486 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_677[k] = f_11 * smg_485[k]
                   + f_4 * snf0_325[k]
                   - f_5 * snf1_325[k]
                   + f_3 * pc_x[k] * sng_485[k];

        t_678[k] = f_11 * smg_486[k]
                   + f_6 * snf0_326[k]
                   - f_7 * snf1_326[k]
                   + f_3 * pc_x[k] * sng_486[k];

        t_679[k] = f_16 * smg_363[k]
                   + f_3 * pc_z[k] * sng_483[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, t_683, pc_x, pc_y, smg_380, smg_489, smg_490, \
                         smg_491, snf0_329, snf1_329, sng_485, sng_489, sng_490, \
                         sng_491 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = f_11 * smg_380[k]
                   + f_3 * pc_y[k] * sng_485[k];

        t_681[k] = f_11 * smg_489[k]
                   + f_6 * snf0_329[k]
                   - f_7 * snf1_329[k]
                   + f_3 * pc_x[k] * sng_489[k];

        t_682[k] = f_11 * smg_490[k]
                   + f_3 * pc_x[k] * sng_490[k];

        t_683[k] = f_11 * smg_491[k]
                   + f_3 * pc_x[k] * sng_491[k];
    }

#pragma omp simd aligned(t_684, t_685, t_686, t_687, pc_x, pc_y, smg_385, smg_492, smg_493, \
                         smg_494, snf0_326, snf1_326, sng_490, sng_492, sng_493, \
                         sng_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_684[k] = f_11 * smg_492[k]
                   + f_3 * pc_x[k] * sng_492[k];

        t_685[k] = f_11 * smg_493[k]
                   + f_3 * pc_x[k] * sng_493[k];

        t_686[k] = f_11 * smg_494[k]
                   + f_3 * pc_x[k] * sng_494[k];

        t_687[k] = f_11 * smg_385[k]
                   + f_1 * snf0_326[k]
                   - f_2 * snf1_326[k]
                   + f_3 * pc_y[k] * sng_490[k];
    }

#pragma omp simd aligned(t_688, t_689, t_690, pc_y, pc_z, smg_370, smg_387, smg_388, snf0_328, \
                         snf0_329, snf1_328, snf1_329, sng_490, sng_492, \
                         sng_493 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_688[k] = f_16 * smg_370[k]
                   + f_3 * pc_z[k] * sng_490[k];

        t_689[k] = f_11 * smg_387[k]
                   + f_4 * snf0_328[k]
                   - f_5 * snf1_328[k]
                   + f_3 * pc_y[k] * sng_492[k];

        t_690[k] = f_11 * smg_388[k]
                   + f_6 * snf0_329[k]
                   - f_7 * snf1_329[k]
                   + f_3 * pc_y[k] * sng_493[k];
    }

#pragma omp simd aligned(t_691, t_692, t_693, pc_x, pc_y, pc_z, smg_374, smg_389, smg_495, \
                         snf0_329, snf0_330, snf1_329, snf1_330, sng_494, \
                         sng_495 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_691[k] = f_11 * smg_389[k]
                   + f_3 * pc_y[k] * sng_494[k];

        t_692[k] = f_16 * smg_374[k]
                   + f_1 * snf0_329[k]
                   - f_2 * snf1_329[k]
                   + f_3 * pc_z[k] * sng_494[k];

        t_693[k] = f_11 * smg_495[k]
                   + f_1 * snf0_330[k]
                   - f_2 * snf1_330[k]
                   + f_3 * pc_x[k] * sng_495[k];
    }
}

static auto
compute_prim_snh_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t smh0,
                                                          const size_t smg, const size_t smh1,
                                                          const size_t snf0, const size_t snf1,
                                                          const size_t sng, const size_t ncols,
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
    const auto f_13 = 4.0 / q;
    const auto f_14 = 3.5 / q;
    const auto f_15 = 3.0 / q;
    const auto f_17 = 2.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *smh0_567 = buffer.data(smh0 + 567);
    const auto *smh0_570 = buffer.data(smh0 + 570);
    const auto *smh0_572 = buffer.data(smh0 + 572);
    const auto *smh0_573 = buffer.data(smh0 + 573);
    const auto *smh0_576 = buffer.data(smh0 + 576);
    const auto *smh0_587 = buffer.data(smh0 + 587);
    const auto *smh0_588 = buffer.data(smh0 + 588);
    const auto *smh0_591 = buffer.data(smh0 + 591);
    const auto *smh0_594 = buffer.data(smh0 + 594);
    const auto *smh0_603 = buffer.data(smh0 + 603);

    const auto *smg_375 = buffer.data(smg + 375);
    const auto *smg_378 = buffer.data(smg + 378);
    const auto *smg_385 = buffer.data(smg + 385);
    const auto *smg_389 = buffer.data(smg + 389);
    const auto *smg_390 = buffer.data(smg + 390);
    const auto *smg_392 = buffer.data(smg + 392);
    const auto *smg_393 = buffer.data(smg + 393);
    const auto *smg_395 = buffer.data(smg + 395);
    const auto *smg_400 = buffer.data(smg + 400);
    const auto *smg_402 = buffer.data(smg + 402);
    const auto *smg_403 = buffer.data(smg + 403);
    const auto *smg_404 = buffer.data(smg + 404);
    const auto *smg_405 = buffer.data(smg + 405);
    const auto *smg_406 = buffer.data(smg + 406);
    const auto *smg_407 = buffer.data(smg + 407);
    const auto *smg_408 = buffer.data(smg + 408);
    const auto *smg_410 = buffer.data(smg + 410);
    const auto *smg_415 = buffer.data(smg + 415);
    const auto *smg_417 = buffer.data(smg + 417);
    const auto *smg_418 = buffer.data(smg + 418);
    const auto *smg_419 = buffer.data(smg + 419);
    const auto *smg_420 = buffer.data(smg + 420);
    const auto *smg_422 = buffer.data(smg + 422);
    const auto *smg_423 = buffer.data(smg + 423);
    const auto *smg_425 = buffer.data(smg + 425);
    const auto *smg_430 = buffer.data(smg + 430);
    const auto *smg_432 = buffer.data(smg + 432);
    const auto *smg_433 = buffer.data(smg + 433);
    const auto *smg_434 = buffer.data(smg + 434);
    const auto *smg_435 = buffer.data(smg + 435);
    const auto *smg_437 = buffer.data(smg + 437);
    const auto *smg_438 = buffer.data(smg + 438);
    const auto *smg_440 = buffer.data(smg + 440);
    const auto *smg_447 = buffer.data(smg + 447);
    const auto *smg_448 = buffer.data(smg + 448);
    const auto *smg_449 = buffer.data(smg + 449);
    const auto *smg_450 = buffer.data(smg + 450);
    const auto *smg_452 = buffer.data(smg + 452);
    const auto *smg_455 = buffer.data(smg + 455);
    const auto *smg_498 = buffer.data(smg + 498);
    const auto *smg_500 = buffer.data(smg + 500);
    const auto *smg_501 = buffer.data(smg + 501);
    const auto *smg_504 = buffer.data(smg + 504);
    const auto *smg_505 = buffer.data(smg + 505);
    const auto *smg_506 = buffer.data(smg + 506);
    const auto *smg_507 = buffer.data(smg + 507);
    const auto *smg_508 = buffer.data(smg + 508);
    const auto *smg_509 = buffer.data(smg + 509);
    const auto *smg_520 = buffer.data(smg + 520);
    const auto *smg_521 = buffer.data(smg + 521);
    const auto *smg_522 = buffer.data(smg + 522);
    const auto *smg_523 = buffer.data(smg + 523);
    const auto *smg_524 = buffer.data(smg + 524);
    const auto *smg_525 = buffer.data(smg + 525);
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
    const auto *smg_543 = buffer.data(smg + 543);
    const auto *smg_545 = buffer.data(smg + 545);
    const auto *smg_546 = buffer.data(smg + 546);
    const auto *smg_549 = buffer.data(smg + 549);
    const auto *smg_550 = buffer.data(smg + 550);
    const auto *smg_551 = buffer.data(smg + 551);
    const auto *smg_552 = buffer.data(smg + 552);
    const auto *smg_553 = buffer.data(smg + 553);
    const auto *smg_554 = buffer.data(smg + 554);
    const auto *smg_560 = buffer.data(smg + 560);
    const auto *smg_564 = buffer.data(smg + 564);
    const auto *smg_565 = buffer.data(smg + 565);
    const auto *smg_566 = buffer.data(smg + 566);
    const auto *smg_567 = buffer.data(smg + 567);
    const auto *smg_568 = buffer.data(smg + 568);
    const auto *smg_569 = buffer.data(smg + 569);
    const auto *smg_570 = buffer.data(smg + 570);
    const auto *smg_573 = buffer.data(smg + 573);
    const auto *smg_575 = buffer.data(smg + 575);
    const auto *smg_576 = buffer.data(smg + 576);
    const auto *smg_579 = buffer.data(smg + 579);
    const auto *smg_580 = buffer.data(smg + 580);
    const auto *smg_581 = buffer.data(smg + 581);

    const auto *smh1_567 = buffer.data(smh1 + 567);
    const auto *smh1_570 = buffer.data(smh1 + 570);
    const auto *smh1_572 = buffer.data(smh1 + 572);
    const auto *smh1_573 = buffer.data(smh1 + 573);
    const auto *smh1_576 = buffer.data(smh1 + 576);
    const auto *smh1_587 = buffer.data(smh1 + 587);
    const auto *smh1_588 = buffer.data(smh1 + 588);
    const auto *smh1_591 = buffer.data(smh1 + 591);
    const auto *smh1_594 = buffer.data(smh1 + 594);
    const auto *smh1_603 = buffer.data(smh1 + 603);

    const auto *snf0_333 = buffer.data(snf0 + 333);
    const auto *snf0_335 = buffer.data(snf0 + 335);
    const auto *snf0_336 = buffer.data(snf0 + 336);
    const auto *snf0_338 = buffer.data(snf0 + 338);
    const auto *snf0_339 = buffer.data(snf0 + 339);
    const auto *snf0_346 = buffer.data(snf0 + 346);
    const auto *snf0_348 = buffer.data(snf0 + 348);
    const auto *snf0_349 = buffer.data(snf0 + 349);
    const auto *snf0_350 = buffer.data(snf0 + 350);
    const auto *snf0_353 = buffer.data(snf0 + 353);
    const auto *snf0_355 = buffer.data(snf0 + 355);
    const auto *snf0_356 = buffer.data(snf0 + 356);
    const auto *snf0_358 = buffer.data(snf0 + 358);
    const auto *snf0_359 = buffer.data(snf0 + 359);
    const auto *snf0_360 = buffer.data(snf0 + 360);
    const auto *snf0_363 = buffer.data(snf0 + 363);
    const auto *snf0_365 = buffer.data(snf0 + 365);
    const auto *snf0_366 = buffer.data(snf0 + 366);
    const auto *snf0_368 = buffer.data(snf0 + 368);
    const auto *snf0_369 = buffer.data(snf0 + 369);
    const auto *snf0_375 = buffer.data(snf0 + 375);
    const auto *snf0_378 = buffer.data(snf0 + 378);
    const auto *snf0_379 = buffer.data(snf0 + 379);
    const auto *snf0_380 = buffer.data(snf0 + 380);
    const auto *snf0_383 = buffer.data(snf0 + 383);
    const auto *snf0_385 = buffer.data(snf0 + 385);
    const auto *snf0_386 = buffer.data(snf0 + 386);
    const auto *snf0_389 = buffer.data(snf0 + 389);

    const auto *snf1_333 = buffer.data(snf1 + 333);
    const auto *snf1_335 = buffer.data(snf1 + 335);
    const auto *snf1_336 = buffer.data(snf1 + 336);
    const auto *snf1_338 = buffer.data(snf1 + 338);
    const auto *snf1_339 = buffer.data(snf1 + 339);
    const auto *snf1_346 = buffer.data(snf1 + 346);
    const auto *snf1_348 = buffer.data(snf1 + 348);
    const auto *snf1_349 = buffer.data(snf1 + 349);
    const auto *snf1_350 = buffer.data(snf1 + 350);
    const auto *snf1_353 = buffer.data(snf1 + 353);
    const auto *snf1_355 = buffer.data(snf1 + 355);
    const auto *snf1_356 = buffer.data(snf1 + 356);
    const auto *snf1_358 = buffer.data(snf1 + 358);
    const auto *snf1_359 = buffer.data(snf1 + 359);
    const auto *snf1_360 = buffer.data(snf1 + 360);
    const auto *snf1_363 = buffer.data(snf1 + 363);
    const auto *snf1_365 = buffer.data(snf1 + 365);
    const auto *snf1_366 = buffer.data(snf1 + 366);
    const auto *snf1_368 = buffer.data(snf1 + 368);
    const auto *snf1_369 = buffer.data(snf1 + 369);
    const auto *snf1_375 = buffer.data(snf1 + 375);
    const auto *snf1_378 = buffer.data(snf1 + 378);
    const auto *snf1_379 = buffer.data(snf1 + 379);
    const auto *snf1_380 = buffer.data(snf1 + 380);
    const auto *snf1_383 = buffer.data(snf1 + 383);
    const auto *snf1_385 = buffer.data(snf1 + 385);
    const auto *snf1_386 = buffer.data(snf1 + 386);
    const auto *snf1_389 = buffer.data(snf1 + 389);

    const auto *sng_495 = buffer.data(sng + 495);
    const auto *sng_497 = buffer.data(sng + 497);
    const auto *sng_498 = buffer.data(sng + 498);
    const auto *sng_500 = buffer.data(sng + 500);
    const auto *sng_501 = buffer.data(sng + 501);
    const auto *sng_504 = buffer.data(sng + 504);
    const auto *sng_505 = buffer.data(sng + 505);
    const auto *sng_506 = buffer.data(sng + 506);
    const auto *sng_507 = buffer.data(sng + 507);
    const auto *sng_508 = buffer.data(sng + 508);
    const auto *sng_509 = buffer.data(sng + 509);
    const auto *sng_510 = buffer.data(sng + 510);
    const auto *sng_512 = buffer.data(sng + 512);
    const auto *sng_513 = buffer.data(sng + 513);
    const auto *sng_515 = buffer.data(sng + 515);
    const auto *sng_520 = buffer.data(sng + 520);
    const auto *sng_521 = buffer.data(sng + 521);
    const auto *sng_522 = buffer.data(sng + 522);
    const auto *sng_523 = buffer.data(sng + 523);
    const auto *sng_524 = buffer.data(sng + 524);
    const auto *sng_525 = buffer.data(sng + 525);
    const auto *sng_527 = buffer.data(sng + 527);
    const auto *sng_528 = buffer.data(sng + 528);
    const auto *sng_530 = buffer.data(sng + 530);
    const auto *sng_531 = buffer.data(sng + 531);
    const auto *sng_534 = buffer.data(sng + 534);
    const auto *sng_535 = buffer.data(sng + 535);
    const auto *sng_536 = buffer.data(sng + 536);
    const auto *sng_537 = buffer.data(sng + 537);
    const auto *sng_538 = buffer.data(sng + 538);
    const auto *sng_539 = buffer.data(sng + 539);
    const auto *sng_540 = buffer.data(sng + 540);
    const auto *sng_542 = buffer.data(sng + 542);
    const auto *sng_543 = buffer.data(sng + 543);
    const auto *sng_545 = buffer.data(sng + 545);
    const auto *sng_546 = buffer.data(sng + 546);
    const auto *sng_549 = buffer.data(sng + 549);
    const auto *sng_550 = buffer.data(sng + 550);
    const auto *sng_551 = buffer.data(sng + 551);
    const auto *sng_552 = buffer.data(sng + 552);
    const auto *sng_553 = buffer.data(sng + 553);
    const auto *sng_554 = buffer.data(sng + 554);
    const auto *sng_555 = buffer.data(sng + 555);
    const auto *sng_557 = buffer.data(sng + 557);
    const auto *sng_558 = buffer.data(sng + 558);
    const auto *sng_560 = buffer.data(sng + 560);
    const auto *sng_564 = buffer.data(sng + 564);
    const auto *sng_565 = buffer.data(sng + 565);
    const auto *sng_566 = buffer.data(sng + 566);
    const auto *sng_567 = buffer.data(sng + 567);
    const auto *sng_568 = buffer.data(sng + 568);
    const auto *sng_569 = buffer.data(sng + 569);
    const auto *sng_570 = buffer.data(sng + 570);
    const auto *sng_572 = buffer.data(sng + 572);
    const auto *sng_573 = buffer.data(sng + 573);
    const auto *sng_575 = buffer.data(sng + 575);
    const auto *sng_576 = buffer.data(sng + 576);
    const auto *sng_579 = buffer.data(sng + 579);
    const auto *sng_580 = buffer.data(sng + 580);
    const auto *sng_581 = buffer.data(sng + 581);

#pragma omp simd aligned(t_694, t_695, t_696, t_697, pc_x, pc_y, pc_z, smg_375, smg_390, \
                         smg_392, smg_498, snf0_333, snf1_333, sng_495, sng_497, \
                         sng_498 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_694[k] = f_10 * smg_390[k]
                   + f_3 * pc_y[k] * sng_495[k];

        t_695[k] = f_17 * smg_375[k]
                   + f_3 * pc_z[k] * sng_495[k];

        t_696[k] = f_11 * smg_498[k]
                   + f_4 * snf0_333[k]
                   - f_5 * snf1_333[k]
                   + f_3 * pc_x[k] * sng_498[k];

        t_697[k] = f_10 * smg_392[k]
                   + f_3 * pc_y[k] * sng_497[k];
    }

#pragma omp simd aligned(t_698, t_699, t_700, pc_x, pc_z, smg_378, smg_500, smg_501, snf0_335, \
                         snf0_336, snf1_335, snf1_336, sng_498, sng_500, \
                         sng_501 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_698[k] = f_11 * smg_500[k]
                   + f_4 * snf0_335[k]
                   - f_5 * snf1_335[k]
                   + f_3 * pc_x[k] * sng_500[k];

        t_699[k] = f_11 * smg_501[k]
                   + f_6 * snf0_336[k]
                   - f_7 * snf1_336[k]
                   + f_3 * pc_x[k] * sng_501[k];

        t_700[k] = f_17 * smg_378[k]
                   + f_3 * pc_z[k] * sng_498[k];
    }

#pragma omp simd aligned(t_701, t_702, t_703, t_704, pc_x, pc_y, smg_395, smg_504, smg_505, \
                         smg_506, snf0_339, snf1_339, sng_500, sng_504, sng_505, \
                         sng_506 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_701[k] = f_10 * smg_395[k]
                   + f_3 * pc_y[k] * sng_500[k];

        t_702[k] = f_11 * smg_504[k]
                   + f_6 * snf0_339[k]
                   - f_7 * snf1_339[k]
                   + f_3 * pc_x[k] * sng_504[k];

        t_703[k] = f_11 * smg_505[k]
                   + f_3 * pc_x[k] * sng_505[k];

        t_704[k] = f_11 * smg_506[k]
                   + f_3 * pc_x[k] * sng_506[k];
    }

#pragma omp simd aligned(t_705, t_706, t_707, t_708, pc_x, pc_y, smg_400, smg_507, smg_508, \
                         smg_509, snf0_336, snf1_336, sng_505, sng_507, sng_508, \
                         sng_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_705[k] = f_11 * smg_507[k]
                   + f_3 * pc_x[k] * sng_507[k];

        t_706[k] = f_11 * smg_508[k]
                   + f_3 * pc_x[k] * sng_508[k];

        t_707[k] = f_11 * smg_509[k]
                   + f_3 * pc_x[k] * sng_509[k];

        t_708[k] = f_10 * smg_400[k]
                   + f_1 * snf0_336[k]
                   - f_2 * snf1_336[k]
                   + f_3 * pc_y[k] * sng_505[k];
    }

#pragma omp simd aligned(t_709, t_710, t_711, pc_y, pc_z, smg_385, smg_402, smg_403, snf0_338, \
                         snf0_339, snf1_338, snf1_339, sng_505, sng_507, \
                         sng_508 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_709[k] = f_17 * smg_385[k]
                   + f_3 * pc_z[k] * sng_505[k];

        t_710[k] = f_10 * smg_402[k]
                   + f_4 * snf0_338[k]
                   - f_5 * snf1_338[k]
                   + f_3 * pc_y[k] * sng_507[k];

        t_711[k] = f_10 * smg_403[k]
                   + f_6 * snf0_339[k]
                   - f_7 * snf1_339[k]
                   + f_3 * pc_y[k] * sng_508[k];
    }

#pragma omp simd aligned(t_712, t_713, t_714, t_715, pb_y, pc_y, pc_z, smh0_567, smg_389, \
                         smg_404, smg_405, smh1_567, snf0_339, snf1_339, sng_509, \
                         sng_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_712[k] = f_10 * smg_404[k]
                   + f_3 * pc_y[k] * sng_509[k];

        t_713[k] = f_17 * smg_389[k]
                   + f_1 * snf0_339[k]
                   - f_2 * snf1_339[k]
                   + f_3 * pc_z[k] * sng_509[k];

        t_714[k] = pb_y[k] * smh0_567[k]
                   - f_8 * pc_y[k] * smh1_567[k];

        t_715[k] = f_9 * smg_405[k]
                   + f_3 * pc_y[k] * sng_510[k];
    }

#pragma omp simd aligned(t_716, t_717, t_718, t_719, pb_y, pc_y, pc_z, smh0_570, smh0_572, \
                         smg_390, smg_406, smg_407, smh1_570, smh1_572, sng_510, \
                         sng_512 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_716[k] = f_15 * smg_390[k]
                   + f_3 * pc_z[k] * sng_510[k];

        t_717[k] = pb_y[k] * smh0_570[k]
                   + f_10 * smg_406[k]
                   - f_8 * pc_y[k] * smh1_570[k];

        t_718[k] = f_9 * smg_407[k]
                   + f_3 * pc_y[k] * sng_512[k];

        t_719[k] = pb_y[k] * smh0_572[k]
                   - f_8 * pc_y[k] * smh1_572[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, t_723, pb_y, pc_y, pc_z, smh0_573, smh0_576, \
                         smg_393, smg_408, smg_410, smh1_573, smh1_576, sng_513, \
                         sng_515 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = pb_y[k] * smh0_573[k]
                   + f_11 * smg_408[k]
                   - f_8 * pc_y[k] * smh1_573[k];

        t_721[k] = f_15 * smg_393[k]
                   + f_3 * pc_z[k] * sng_513[k];

        t_722[k] = f_9 * smg_410[k]
                   + f_3 * pc_y[k] * sng_515[k];

        t_723[k] = pb_y[k] * smh0_576[k]
                   - f_8 * pc_y[k] * smh1_576[k];
    }

#pragma omp simd aligned(t_724, t_725, t_726, t_727, t_728, pc_x, smg_520, smg_521, smg_522, \
                         smg_523, smg_524, sng_520, sng_521, sng_522, sng_523, \
                         sng_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_724[k] = f_11 * smg_520[k]
                   + f_3 * pc_x[k] * sng_520[k];

        t_725[k] = f_11 * smg_521[k]
                   + f_3 * pc_x[k] * sng_521[k];

        t_726[k] = f_11 * smg_522[k]
                   + f_3 * pc_x[k] * sng_522[k];

        t_727[k] = f_11 * smg_523[k]
                   + f_3 * pc_x[k] * sng_523[k];

        t_728[k] = f_11 * smg_524[k]
                   + f_3 * pc_x[k] * sng_524[k];
    }

#pragma omp simd aligned(t_729, t_730, t_731, pc_y, pc_z, smg_400, smg_415, smg_417, snf0_346, \
                         snf0_348, snf1_346, snf1_348, sng_520, \
                         sng_522 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_729[k] = f_9 * smg_415[k]
                   + f_1 * snf0_346[k]
                   - f_2 * snf1_346[k]
                   + f_3 * pc_y[k] * sng_520[k];

        t_730[k] = f_15 * smg_400[k]
                   + f_3 * pc_z[k] * sng_520[k];

        t_731[k] = f_9 * smg_417[k]
                   + f_4 * snf0_348[k]
                   - f_5 * snf1_348[k]
                   + f_3 * pc_y[k] * sng_522[k];
    }

#pragma omp simd aligned(t_732, t_733, t_734, pb_y, pc_y, smh0_587, smg_418, smg_419, \
                         smh1_587, snf0_349, snf1_349, sng_523, \
                         sng_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_732[k] = f_9 * smg_418[k]
                   + f_6 * snf0_349[k]
                   - f_7 * snf1_349[k]
                   + f_3 * pc_y[k] * sng_523[k];

        t_733[k] = f_9 * smg_419[k]
                   + f_3 * pc_y[k] * sng_524[k];

        t_734[k] = pb_y[k] * smh0_587[k]
                   - f_8 * pc_y[k] * smh1_587[k];
    }

#pragma omp simd aligned(t_735, t_736, t_737, t_738, pc_x, pc_y, pc_z, smg_405, smg_525, \
                         smg_528, snf0_350, snf0_353, snf1_350, snf1_353, sng_525, \
                         sng_528 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_735[k] = f_11 * smg_525[k]
                   + f_1 * snf0_350[k]
                   - f_2 * snf1_350[k]
                   + f_3 * pc_x[k] * sng_525[k];

        t_736[k] = f_3 * pc_y[k] * sng_525[k];

        t_737[k] = f_14 * smg_405[k]
                   + f_3 * pc_z[k] * sng_525[k];

        t_738[k] = f_11 * smg_528[k]
                   + f_4 * snf0_353[k]
                   - f_5 * snf1_353[k]
                   + f_3 * pc_x[k] * sng_528[k];
    }

#pragma omp simd aligned(t_739, t_740, t_741, pc_x, pc_y, smg_530, smg_531, snf0_355, \
                         snf0_356, snf1_355, snf1_356, sng_527, sng_530, \
                         sng_531 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_739[k] = f_3 * pc_y[k] * sng_527[k];

        t_740[k] = f_11 * smg_530[k]
                   + f_4 * snf0_355[k]
                   - f_5 * snf1_355[k]
                   + f_3 * pc_x[k] * sng_530[k];

        t_741[k] = f_11 * smg_531[k]
                   + f_6 * snf0_356[k]
                   - f_7 * snf1_356[k]
                   + f_3 * pc_x[k] * sng_531[k];
    }

#pragma omp simd aligned(t_742, t_743, t_744, t_745, pc_x, pc_y, pc_z, smg_408, smg_534, \
                         smg_535, snf0_359, snf1_359, sng_528, sng_530, sng_534, \
                         sng_535 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_742[k] = f_14 * smg_408[k]
                   + f_3 * pc_z[k] * sng_528[k];

        t_743[k] = f_3 * pc_y[k] * sng_530[k];

        t_744[k] = f_11 * smg_534[k]
                   + f_6 * snf0_359[k]
                   - f_7 * snf1_359[k]
                   + f_3 * pc_x[k] * sng_534[k];

        t_745[k] = f_11 * smg_535[k]
                   + f_3 * pc_x[k] * sng_535[k];
    }

#pragma omp simd aligned(t_746, t_747, t_748, t_749, pc_x, smg_536, smg_537, smg_538, smg_539, \
                         sng_536, sng_537, sng_538, sng_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_746[k] = f_11 * smg_536[k]
                   + f_3 * pc_x[k] * sng_536[k];

        t_747[k] = f_11 * smg_537[k]
                   + f_3 * pc_x[k] * sng_537[k];

        t_748[k] = f_11 * smg_538[k]
                   + f_3 * pc_x[k] * sng_538[k];

        t_749[k] = f_11 * smg_539[k]
                   + f_3 * pc_x[k] * sng_539[k];
    }

#pragma omp simd aligned(t_750, t_751, t_752, t_753, pc_y, pc_z, smg_415, snf0_356, snf0_358, \
                         snf0_359, snf1_356, snf1_358, snf1_359, sng_535, sng_537, \
                         sng_538 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_750[k] = f_1 * snf0_356[k]
                   - f_2 * snf1_356[k]
                   + f_3 * pc_y[k] * sng_535[k];

        t_751[k] = f_14 * smg_415[k]
                   + f_3 * pc_z[k] * sng_535[k];

        t_752[k] = f_4 * snf0_358[k]
                   - f_5 * snf1_358[k]
                   + f_3 * pc_y[k] * sng_537[k];

        t_753[k] = f_6 * snf0_359[k]
                   - f_7 * snf1_359[k]
                   + f_3 * pc_y[k] * sng_538[k];
    }

#pragma omp simd aligned(t_754, t_755, t_756, t_757, pc_x, pc_y, pc_z, smg_419, smg_420, \
                         smg_540, snf0_359, snf0_360, snf1_359, snf1_360, sng_539, \
                         sng_540 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_754[k] = f_3 * pc_y[k] * sng_539[k];

        t_755[k] = f_14 * smg_419[k]
                   + f_1 * snf0_359[k]
                   - f_2 * snf1_359[k]
                   + f_3 * pc_z[k] * sng_539[k];

        t_756[k] = f_10 * smg_540[k]
                   + f_1 * snf0_360[k]
                   - f_2 * snf1_360[k]
                   + f_3 * pc_x[k] * sng_540[k];

        t_757[k] = f_13 * smg_420[k]
                   + f_3 * pc_y[k] * sng_540[k];
    }

#pragma omp simd aligned(t_758, t_759, t_760, pc_x, pc_y, pc_z, smg_422, smg_543, snf0_363, \
                         snf1_363, sng_540, sng_542, sng_543 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_758[k] = f_3 * pc_z[k] * sng_540[k];

        t_759[k] = f_10 * smg_543[k]
                   + f_4 * snf0_363[k]
                   - f_5 * snf1_363[k]
                   + f_3 * pc_x[k] * sng_543[k];

        t_760[k] = f_13 * smg_422[k]
                   + f_3 * pc_y[k] * sng_542[k];
    }

#pragma omp simd aligned(t_761, t_762, t_763, pc_x, pc_z, smg_545, smg_546, snf0_365, \
                         snf0_366, snf1_365, snf1_366, sng_543, sng_545, \
                         sng_546 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_761[k] = f_10 * smg_545[k]
                   + f_4 * snf0_365[k]
                   - f_5 * snf1_365[k]
                   + f_3 * pc_x[k] * sng_545[k];

        t_762[k] = f_10 * smg_546[k]
                   + f_6 * snf0_366[k]
                   - f_7 * snf1_366[k]
                   + f_3 * pc_x[k] * sng_546[k];

        t_763[k] = f_3 * pc_z[k] * sng_543[k];
    }

#pragma omp simd aligned(t_764, t_765, t_766, t_767, pc_x, pc_y, smg_425, smg_549, smg_550, \
                         smg_551, snf0_369, snf1_369, sng_545, sng_549, sng_550, \
                         sng_551 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_764[k] = f_13 * smg_425[k]
                   + f_3 * pc_y[k] * sng_545[k];

        t_765[k] = f_10 * smg_549[k]
                   + f_6 * snf0_369[k]
                   - f_7 * snf1_369[k]
                   + f_3 * pc_x[k] * sng_549[k];

        t_766[k] = f_10 * smg_550[k]
                   + f_3 * pc_x[k] * sng_550[k];

        t_767[k] = f_10 * smg_551[k]
                   + f_3 * pc_x[k] * sng_551[k];
    }

#pragma omp simd aligned(t_768, t_769, t_770, t_771, pc_x, pc_y, smg_430, smg_552, smg_553, \
                         smg_554, snf0_366, snf1_366, sng_550, sng_552, sng_553, \
                         sng_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_768[k] = f_10 * smg_552[k]
                   + f_3 * pc_x[k] * sng_552[k];

        t_769[k] = f_10 * smg_553[k]
                   + f_3 * pc_x[k] * sng_553[k];

        t_770[k] = f_10 * smg_554[k]
                   + f_3 * pc_x[k] * sng_554[k];

        t_771[k] = f_13 * smg_430[k]
                   + f_1 * snf0_366[k]
                   - f_2 * snf1_366[k]
                   + f_3 * pc_y[k] * sng_550[k];
    }

#pragma omp simd aligned(t_772, t_773, t_774, pc_y, pc_z, smg_432, smg_433, snf0_368, \
                         snf0_369, snf1_368, snf1_369, sng_550, sng_552, \
                         sng_553 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_772[k] = f_3 * pc_z[k] * sng_550[k];

        t_773[k] = f_13 * smg_432[k]
                   + f_4 * snf0_368[k]
                   - f_5 * snf1_368[k]
                   + f_3 * pc_y[k] * sng_552[k];

        t_774[k] = f_13 * smg_433[k]
                   + f_6 * snf0_369[k]
                   - f_7 * snf1_369[k]
                   + f_3 * pc_y[k] * sng_553[k];
    }

#pragma omp simd aligned(t_775, t_776, t_777, t_778, pb_z, pc_y, pc_z, smh0_588, smg_434, \
                         smg_435, smh1_588, snf0_369, snf1_369, sng_554, \
                         sng_555 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_775[k] = f_13 * smg_434[k]
                   + f_3 * pc_y[k] * sng_554[k];

        t_776[k] = f_1 * snf0_369[k]
                   - f_2 * snf1_369[k]
                   + f_3 * pc_z[k] * sng_554[k];

        t_777[k] = pb_z[k] * smh0_588[k]
                   - f_8 * pc_z[k] * smh1_588[k];

        t_778[k] = f_14 * smg_435[k]
                   + f_3 * pc_y[k] * sng_555[k];
    }

#pragma omp simd aligned(t_779, t_780, t_781, pb_z, pc_y, pc_z, smh0_591, smg_420, smg_437, \
                         smh1_591, sng_555, sng_557 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_779[k] = f_9 * smg_420[k]
                   + f_3 * pc_z[k] * sng_555[k];

        t_780[k] = pb_z[k] * smh0_591[k]
                   - f_8 * pc_z[k] * smh1_591[k];

        t_781[k] = f_14 * smg_437[k]
                   + f_3 * pc_y[k] * sng_557[k];
    }

#pragma omp simd aligned(t_782, t_783, t_784, pb_z, pc_x, pc_z, smh0_594, smg_423, smg_560, \
                         smh1_594, snf0_375, snf1_375, sng_558, \
                         sng_560 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_782[k] = f_10 * smg_560[k]
                   + f_4 * snf0_375[k]
                   - f_5 * snf1_375[k]
                   + f_3 * pc_x[k] * sng_560[k];

        t_783[k] = pb_z[k] * smh0_594[k]
                   - f_8 * pc_z[k] * smh1_594[k];

        t_784[k] = f_9 * smg_423[k]
                   + f_3 * pc_z[k] * sng_558[k];
    }

#pragma omp simd aligned(t_785, t_786, t_787, t_788, pc_x, pc_y, smg_440, smg_564, smg_565, \
                         smg_566, snf0_379, snf1_379, sng_560, sng_564, sng_565, \
                         sng_566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_785[k] = f_14 * smg_440[k]
                   + f_3 * pc_y[k] * sng_560[k];

        t_786[k] = f_10 * smg_564[k]
                   + f_6 * snf0_379[k]
                   - f_7 * snf1_379[k]
                   + f_3 * pc_x[k] * sng_564[k];

        t_787[k] = f_10 * smg_565[k]
                   + f_3 * pc_x[k] * sng_565[k];

        t_788[k] = f_10 * smg_566[k]
                   + f_3 * pc_x[k] * sng_566[k];
    }

#pragma omp simd aligned(t_789, t_790, t_791, t_792, pb_z, pc_x, pc_z, smh0_603, smg_567, \
                         smg_568, smg_569, smh1_603, sng_567, sng_568, \
                         sng_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_789[k] = f_10 * smg_567[k]
                   + f_3 * pc_x[k] * sng_567[k];

        t_790[k] = f_10 * smg_568[k]
                   + f_3 * pc_x[k] * sng_568[k];

        t_791[k] = f_10 * smg_569[k]
                   + f_3 * pc_x[k] * sng_569[k];

        t_792[k] = pb_z[k] * smh0_603[k]
                   - f_8 * pc_z[k] * smh1_603[k];
    }

#pragma omp simd aligned(t_793, t_794, t_795, pc_y, pc_z, smg_430, smg_447, smg_448, snf0_378, \
                         snf0_379, snf1_378, snf1_379, sng_565, sng_567, \
                         sng_568 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_793[k] = f_9 * smg_430[k]
                   + f_3 * pc_z[k] * sng_565[k];

        t_794[k] = f_14 * smg_447[k]
                   + f_4 * snf0_378[k]
                   - f_5 * snf1_378[k]
                   + f_3 * pc_y[k] * sng_567[k];

        t_795[k] = f_14 * smg_448[k]
                   + f_6 * snf0_379[k]
                   - f_7 * snf1_379[k]
                   + f_3 * pc_y[k] * sng_568[k];
    }

#pragma omp simd aligned(t_796, t_797, t_798, pc_x, pc_y, pc_z, smg_434, smg_449, smg_570, \
                         snf0_379, snf0_380, snf1_379, snf1_380, sng_569, \
                         sng_570 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_796[k] = f_14 * smg_449[k]
                   + f_3 * pc_y[k] * sng_569[k];

        t_797[k] = f_9 * smg_434[k]
                   + f_1 * snf0_379[k]
                   - f_2 * snf1_379[k]
                   + f_3 * pc_z[k] * sng_569[k];

        t_798[k] = f_10 * smg_570[k]
                   + f_1 * snf0_380[k]
                   - f_2 * snf1_380[k]
                   + f_3 * pc_x[k] * sng_570[k];
    }

#pragma omp simd aligned(t_799, t_800, t_801, t_802, pc_x, pc_y, pc_z, smg_435, smg_450, \
                         smg_452, smg_573, snf0_383, snf1_383, sng_570, sng_572, \
                         sng_573 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_799[k] = f_15 * smg_450[k]
                   + f_3 * pc_y[k] * sng_570[k];

        t_800[k] = f_10 * smg_435[k]
                   + f_3 * pc_z[k] * sng_570[k];

        t_801[k] = f_10 * smg_573[k]
                   + f_4 * snf0_383[k]
                   - f_5 * snf1_383[k]
                   + f_3 * pc_x[k] * sng_573[k];

        t_802[k] = f_15 * smg_452[k]
                   + f_3 * pc_y[k] * sng_572[k];
    }

#pragma omp simd aligned(t_803, t_804, t_805, pc_x, pc_z, smg_438, smg_575, smg_576, snf0_385, \
                         snf0_386, snf1_385, snf1_386, sng_573, sng_575, \
                         sng_576 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_803[k] = f_10 * smg_575[k]
                   + f_4 * snf0_385[k]
                   - f_5 * snf1_385[k]
                   + f_3 * pc_x[k] * sng_575[k];

        t_804[k] = f_10 * smg_576[k]
                   + f_6 * snf0_386[k]
                   - f_7 * snf1_386[k]
                   + f_3 * pc_x[k] * sng_576[k];

        t_805[k] = f_10 * smg_438[k]
                   + f_3 * pc_z[k] * sng_573[k];
    }

#pragma omp simd aligned(t_806, t_807, t_808, t_809, pc_x, pc_y, smg_455, smg_579, smg_580, \
                         smg_581, snf0_389, snf1_389, sng_575, sng_579, sng_580, \
                         sng_581 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_806[k] = f_15 * smg_455[k]
                   + f_3 * pc_y[k] * sng_575[k];

        t_807[k] = f_10 * smg_579[k]
                   + f_6 * snf0_389[k]
                   - f_7 * snf1_389[k]
                   + f_3 * pc_x[k] * sng_579[k];

        t_808[k] = f_10 * smg_580[k]
                   + f_3 * pc_x[k] * sng_580[k];

        t_809[k] = f_10 * smg_581[k]
                   + f_3 * pc_x[k] * sng_581[k];
    }
}

static auto
compute_prim_snh_three_center_electron_repulsion_0_piece7(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t smh0,
                                                          const size_t smg, const size_t smh1,
                                                          const size_t snf0, const size_t snf1,
                                                          const size_t sng, const size_t ncols,
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
    const auto f_14 = 3.5 / q;
    const auto f_15 = 3.0 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *smh0_735 = buffer.data(smh0 + 735);
    const auto *smh0_738 = buffer.data(smh0 + 738);
    const auto *smh0_740 = buffer.data(smh0 + 740);
    const auto *smh0_741 = buffer.data(smh0 + 741);
    const auto *smh0_744 = buffer.data(smh0 + 744);

    const auto *smg_445 = buffer.data(smg + 445);
    const auto *smg_449 = buffer.data(smg + 449);
    const auto *smg_450 = buffer.data(smg + 450);
    const auto *smg_453 = buffer.data(smg + 453);
    const auto *smg_460 = buffer.data(smg + 460);
    const auto *smg_462 = buffer.data(smg + 462);
    const auto *smg_463 = buffer.data(smg + 463);
    const auto *smg_464 = buffer.data(smg + 464);
    const auto *smg_465 = buffer.data(smg + 465);
    const auto *smg_467 = buffer.data(smg + 467);
    const auto *smg_468 = buffer.data(smg + 468);
    const auto *smg_470 = buffer.data(smg + 470);
    const auto *smg_475 = buffer.data(smg + 475);
    const auto *smg_477 = buffer.data(smg + 477);
    const auto *smg_478 = buffer.data(smg + 478);
    const auto *smg_479 = buffer.data(smg + 479);
    const auto *smg_480 = buffer.data(smg + 480);
    const auto *smg_482 = buffer.data(smg + 482);
    const auto *smg_483 = buffer.data(smg + 483);
    const auto *smg_485 = buffer.data(smg + 485);
    const auto *smg_490 = buffer.data(smg + 490);
    const auto *smg_492 = buffer.data(smg + 492);
    const auto *smg_493 = buffer.data(smg + 493);
    const auto *smg_494 = buffer.data(smg + 494);
    const auto *smg_495 = buffer.data(smg + 495);
    const auto *smg_497 = buffer.data(smg + 497);
    const auto *smg_498 = buffer.data(smg + 498);
    const auto *smg_500 = buffer.data(smg + 500);
    const auto *smg_505 = buffer.data(smg + 505);
    const auto *smg_507 = buffer.data(smg + 507);
    const auto *smg_508 = buffer.data(smg + 508);
    const auto *smg_509 = buffer.data(smg + 509);
    const auto *smg_510 = buffer.data(smg + 510);
    const auto *smg_512 = buffer.data(smg + 512);
    const auto *smg_513 = buffer.data(smg + 513);
    const auto *smg_515 = buffer.data(smg + 515);
    const auto *smg_520 = buffer.data(smg + 520);
    const auto *smg_522 = buffer.data(smg + 522);
    const auto *smg_523 = buffer.data(smg + 523);
    const auto *smg_524 = buffer.data(smg + 524);
    const auto *smg_525 = buffer.data(smg + 525);
    const auto *smg_526 = buffer.data(smg + 526);
    const auto *smg_527 = buffer.data(smg + 527);
    const auto *smg_528 = buffer.data(smg + 528);
    const auto *smg_530 = buffer.data(smg + 530);
    const auto *smg_535 = buffer.data(smg + 535);
    const auto *smg_537 = buffer.data(smg + 537);
    const auto *smg_582 = buffer.data(smg + 582);
    const auto *smg_583 = buffer.data(smg + 583);
    const auto *smg_584 = buffer.data(smg + 584);
    const auto *smg_585 = buffer.data(smg + 585);
    const auto *smg_588 = buffer.data(smg + 588);
    const auto *smg_590 = buffer.data(smg + 590);
    const auto *smg_591 = buffer.data(smg + 591);
    const auto *smg_594 = buffer.data(smg + 594);
    const auto *smg_595 = buffer.data(smg + 595);
    const auto *smg_596 = buffer.data(smg + 596);
    const auto *smg_597 = buffer.data(smg + 597);
    const auto *smg_598 = buffer.data(smg + 598);
    const auto *smg_599 = buffer.data(smg + 599);
    const auto *smg_600 = buffer.data(smg + 600);
    const auto *smg_603 = buffer.data(smg + 603);
    const auto *smg_605 = buffer.data(smg + 605);
    const auto *smg_606 = buffer.data(smg + 606);
    const auto *smg_609 = buffer.data(smg + 609);
    const auto *smg_610 = buffer.data(smg + 610);
    const auto *smg_611 = buffer.data(smg + 611);
    const auto *smg_612 = buffer.data(smg + 612);
    const auto *smg_613 = buffer.data(smg + 613);
    const auto *smg_614 = buffer.data(smg + 614);
    const auto *smg_615 = buffer.data(smg + 615);
    const auto *smg_618 = buffer.data(smg + 618);
    const auto *smg_620 = buffer.data(smg + 620);
    const auto *smg_621 = buffer.data(smg + 621);
    const auto *smg_624 = buffer.data(smg + 624);
    const auto *smg_625 = buffer.data(smg + 625);
    const auto *smg_626 = buffer.data(smg + 626);
    const auto *smg_627 = buffer.data(smg + 627);
    const auto *smg_628 = buffer.data(smg + 628);
    const auto *smg_629 = buffer.data(smg + 629);
    const auto *smg_630 = buffer.data(smg + 630);
    const auto *smg_633 = buffer.data(smg + 633);
    const auto *smg_635 = buffer.data(smg + 635);
    const auto *smg_636 = buffer.data(smg + 636);
    const auto *smg_639 = buffer.data(smg + 639);
    const auto *smg_640 = buffer.data(smg + 640);
    const auto *smg_641 = buffer.data(smg + 641);
    const auto *smg_642 = buffer.data(smg + 642);
    const auto *smg_643 = buffer.data(smg + 643);
    const auto *smg_644 = buffer.data(smg + 644);
    const auto *smg_655 = buffer.data(smg + 655);
    const auto *smg_656 = buffer.data(smg + 656);
    const auto *smg_657 = buffer.data(smg + 657);
    const auto *smg_658 = buffer.data(smg + 658);
    const auto *smg_659 = buffer.data(smg + 659);

    const auto *smh1_735 = buffer.data(smh1 + 735);
    const auto *smh1_738 = buffer.data(smh1 + 738);
    const auto *smh1_740 = buffer.data(smh1 + 740);
    const auto *smh1_741 = buffer.data(smh1 + 741);
    const auto *smh1_744 = buffer.data(smh1 + 744);

    const auto *snf0_386 = buffer.data(snf0 + 386);
    const auto *snf0_388 = buffer.data(snf0 + 388);
    const auto *snf0_389 = buffer.data(snf0 + 389);
    const auto *snf0_390 = buffer.data(snf0 + 390);
    const auto *snf0_393 = buffer.data(snf0 + 393);
    const auto *snf0_395 = buffer.data(snf0 + 395);
    const auto *snf0_396 = buffer.data(snf0 + 396);
    const auto *snf0_398 = buffer.data(snf0 + 398);
    const auto *snf0_399 = buffer.data(snf0 + 399);
    const auto *snf0_400 = buffer.data(snf0 + 400);
    const auto *snf0_403 = buffer.data(snf0 + 403);
    const auto *snf0_405 = buffer.data(snf0 + 405);
    const auto *snf0_406 = buffer.data(snf0 + 406);
    const auto *snf0_408 = buffer.data(snf0 + 408);
    const auto *snf0_409 = buffer.data(snf0 + 409);
    const auto *snf0_410 = buffer.data(snf0 + 410);
    const auto *snf0_413 = buffer.data(snf0 + 413);
    const auto *snf0_415 = buffer.data(snf0 + 415);
    const auto *snf0_416 = buffer.data(snf0 + 416);
    const auto *snf0_418 = buffer.data(snf0 + 418);
    const auto *snf0_419 = buffer.data(snf0 + 419);
    const auto *snf0_420 = buffer.data(snf0 + 420);
    const auto *snf0_423 = buffer.data(snf0 + 423);
    const auto *snf0_425 = buffer.data(snf0 + 425);
    const auto *snf0_426 = buffer.data(snf0 + 426);
    const auto *snf0_428 = buffer.data(snf0 + 428);
    const auto *snf0_429 = buffer.data(snf0 + 429);
    const auto *snf0_436 = buffer.data(snf0 + 436);
    const auto *snf0_438 = buffer.data(snf0 + 438);

    const auto *snf1_386 = buffer.data(snf1 + 386);
    const auto *snf1_388 = buffer.data(snf1 + 388);
    const auto *snf1_389 = buffer.data(snf1 + 389);
    const auto *snf1_390 = buffer.data(snf1 + 390);
    const auto *snf1_393 = buffer.data(snf1 + 393);
    const auto *snf1_395 = buffer.data(snf1 + 395);
    const auto *snf1_396 = buffer.data(snf1 + 396);
    const auto *snf1_398 = buffer.data(snf1 + 398);
    const auto *snf1_399 = buffer.data(snf1 + 399);
    const auto *snf1_400 = buffer.data(snf1 + 400);
    const auto *snf1_403 = buffer.data(snf1 + 403);
    const auto *snf1_405 = buffer.data(snf1 + 405);
    const auto *snf1_406 = buffer.data(snf1 + 406);
    const auto *snf1_408 = buffer.data(snf1 + 408);
    const auto *snf1_409 = buffer.data(snf1 + 409);
    const auto *snf1_410 = buffer.data(snf1 + 410);
    const auto *snf1_413 = buffer.data(snf1 + 413);
    const auto *snf1_415 = buffer.data(snf1 + 415);
    const auto *snf1_416 = buffer.data(snf1 + 416);
    const auto *snf1_418 = buffer.data(snf1 + 418);
    const auto *snf1_419 = buffer.data(snf1 + 419);
    const auto *snf1_420 = buffer.data(snf1 + 420);
    const auto *snf1_423 = buffer.data(snf1 + 423);
    const auto *snf1_425 = buffer.data(snf1 + 425);
    const auto *snf1_426 = buffer.data(snf1 + 426);
    const auto *snf1_428 = buffer.data(snf1 + 428);
    const auto *snf1_429 = buffer.data(snf1 + 429);
    const auto *snf1_436 = buffer.data(snf1 + 436);
    const auto *snf1_438 = buffer.data(snf1 + 438);

    const auto *sng_580 = buffer.data(sng + 580);
    const auto *sng_582 = buffer.data(sng + 582);
    const auto *sng_583 = buffer.data(sng + 583);
    const auto *sng_584 = buffer.data(sng + 584);
    const auto *sng_585 = buffer.data(sng + 585);
    const auto *sng_587 = buffer.data(sng + 587);
    const auto *sng_588 = buffer.data(sng + 588);
    const auto *sng_590 = buffer.data(sng + 590);
    const auto *sng_591 = buffer.data(sng + 591);
    const auto *sng_594 = buffer.data(sng + 594);
    const auto *sng_595 = buffer.data(sng + 595);
    const auto *sng_596 = buffer.data(sng + 596);
    const auto *sng_597 = buffer.data(sng + 597);
    const auto *sng_598 = buffer.data(sng + 598);
    const auto *sng_599 = buffer.data(sng + 599);
    const auto *sng_600 = buffer.data(sng + 600);
    const auto *sng_602 = buffer.data(sng + 602);
    const auto *sng_603 = buffer.data(sng + 603);
    const auto *sng_605 = buffer.data(sng + 605);
    const auto *sng_606 = buffer.data(sng + 606);
    const auto *sng_609 = buffer.data(sng + 609);
    const auto *sng_610 = buffer.data(sng + 610);
    const auto *sng_611 = buffer.data(sng + 611);
    const auto *sng_612 = buffer.data(sng + 612);
    const auto *sng_613 = buffer.data(sng + 613);
    const auto *sng_614 = buffer.data(sng + 614);
    const auto *sng_615 = buffer.data(sng + 615);
    const auto *sng_617 = buffer.data(sng + 617);
    const auto *sng_618 = buffer.data(sng + 618);
    const auto *sng_620 = buffer.data(sng + 620);
    const auto *sng_621 = buffer.data(sng + 621);
    const auto *sng_624 = buffer.data(sng + 624);
    const auto *sng_625 = buffer.data(sng + 625);
    const auto *sng_626 = buffer.data(sng + 626);
    const auto *sng_627 = buffer.data(sng + 627);
    const auto *sng_628 = buffer.data(sng + 628);
    const auto *sng_629 = buffer.data(sng + 629);
    const auto *sng_630 = buffer.data(sng + 630);
    const auto *sng_632 = buffer.data(sng + 632);
    const auto *sng_633 = buffer.data(sng + 633);
    const auto *sng_635 = buffer.data(sng + 635);
    const auto *sng_636 = buffer.data(sng + 636);
    const auto *sng_639 = buffer.data(sng + 639);
    const auto *sng_640 = buffer.data(sng + 640);
    const auto *sng_641 = buffer.data(sng + 641);
    const auto *sng_642 = buffer.data(sng + 642);
    const auto *sng_643 = buffer.data(sng + 643);
    const auto *sng_644 = buffer.data(sng + 644);
    const auto *sng_645 = buffer.data(sng + 645);
    const auto *sng_647 = buffer.data(sng + 647);
    const auto *sng_648 = buffer.data(sng + 648);
    const auto *sng_650 = buffer.data(sng + 650);
    const auto *sng_655 = buffer.data(sng + 655);
    const auto *sng_656 = buffer.data(sng + 656);
    const auto *sng_657 = buffer.data(sng + 657);
    const auto *sng_658 = buffer.data(sng + 658);
    const auto *sng_659 = buffer.data(sng + 659);

#pragma omp simd aligned(t_810, t_811, t_812, t_813, pc_x, pc_y, smg_460, smg_582, smg_583, \
                         smg_584, snf0_386, snf1_386, sng_580, sng_582, sng_583, \
                         sng_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_810[k] = f_10 * smg_582[k]
                   + f_3 * pc_x[k] * sng_582[k];

        t_811[k] = f_10 * smg_583[k]
                   + f_3 * pc_x[k] * sng_583[k];

        t_812[k] = f_10 * smg_584[k]
                   + f_3 * pc_x[k] * sng_584[k];

        t_813[k] = f_15 * smg_460[k]
                   + f_1 * snf0_386[k]
                   - f_2 * snf1_386[k]
                   + f_3 * pc_y[k] * sng_580[k];
    }

#pragma omp simd aligned(t_814, t_815, t_816, pc_y, pc_z, smg_445, smg_462, smg_463, snf0_388, \
                         snf0_389, snf1_388, snf1_389, sng_580, sng_582, \
                         sng_583 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_814[k] = f_10 * smg_445[k]
                   + f_3 * pc_z[k] * sng_580[k];

        t_815[k] = f_15 * smg_462[k]
                   + f_4 * snf0_388[k]
                   - f_5 * snf1_388[k]
                   + f_3 * pc_y[k] * sng_582[k];

        t_816[k] = f_15 * smg_463[k]
                   + f_6 * snf0_389[k]
                   - f_7 * snf1_389[k]
                   + f_3 * pc_y[k] * sng_583[k];
    }

#pragma omp simd aligned(t_817, t_818, t_819, pc_x, pc_y, pc_z, smg_449, smg_464, smg_585, \
                         snf0_389, snf0_390, snf1_389, snf1_390, sng_584, \
                         sng_585 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_817[k] = f_15 * smg_464[k]
                   + f_3 * pc_y[k] * sng_584[k];

        t_818[k] = f_10 * smg_449[k]
                   + f_1 * snf0_389[k]
                   - f_2 * snf1_389[k]
                   + f_3 * pc_z[k] * sng_584[k];

        t_819[k] = f_10 * smg_585[k]
                   + f_1 * snf0_390[k]
                   - f_2 * snf1_390[k]
                   + f_3 * pc_x[k] * sng_585[k];
    }

#pragma omp simd aligned(t_820, t_821, t_822, t_823, pc_x, pc_y, pc_z, smg_450, smg_465, \
                         smg_467, smg_588, snf0_393, snf1_393, sng_585, sng_587, \
                         sng_588 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_820[k] = f_17 * smg_465[k]
                   + f_3 * pc_y[k] * sng_585[k];

        t_821[k] = f_11 * smg_450[k]
                   + f_3 * pc_z[k] * sng_585[k];

        t_822[k] = f_10 * smg_588[k]
                   + f_4 * snf0_393[k]
                   - f_5 * snf1_393[k]
                   + f_3 * pc_x[k] * sng_588[k];

        t_823[k] = f_17 * smg_467[k]
                   + f_3 * pc_y[k] * sng_587[k];
    }

#pragma omp simd aligned(t_824, t_825, t_826, pc_x, pc_z, smg_453, smg_590, smg_591, snf0_395, \
                         snf0_396, snf1_395, snf1_396, sng_588, sng_590, \
                         sng_591 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_824[k] = f_10 * smg_590[k]
                   + f_4 * snf0_395[k]
                   - f_5 * snf1_395[k]
                   + f_3 * pc_x[k] * sng_590[k];

        t_825[k] = f_10 * smg_591[k]
                   + f_6 * snf0_396[k]
                   - f_7 * snf1_396[k]
                   + f_3 * pc_x[k] * sng_591[k];

        t_826[k] = f_11 * smg_453[k]
                   + f_3 * pc_z[k] * sng_588[k];
    }

#pragma omp simd aligned(t_827, t_828, t_829, t_830, pc_x, pc_y, smg_470, smg_594, smg_595, \
                         smg_596, snf0_399, snf1_399, sng_590, sng_594, sng_595, \
                         sng_596 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_827[k] = f_17 * smg_470[k]
                   + f_3 * pc_y[k] * sng_590[k];

        t_828[k] = f_10 * smg_594[k]
                   + f_6 * snf0_399[k]
                   - f_7 * snf1_399[k]
                   + f_3 * pc_x[k] * sng_594[k];

        t_829[k] = f_10 * smg_595[k]
                   + f_3 * pc_x[k] * sng_595[k];

        t_830[k] = f_10 * smg_596[k]
                   + f_3 * pc_x[k] * sng_596[k];
    }

#pragma omp simd aligned(t_831, t_832, t_833, t_834, pc_x, pc_y, smg_475, smg_597, smg_598, \
                         smg_599, snf0_396, snf1_396, sng_595, sng_597, sng_598, \
                         sng_599 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_831[k] = f_10 * smg_597[k]
                   + f_3 * pc_x[k] * sng_597[k];

        t_832[k] = f_10 * smg_598[k]
                   + f_3 * pc_x[k] * sng_598[k];

        t_833[k] = f_10 * smg_599[k]
                   + f_3 * pc_x[k] * sng_599[k];

        t_834[k] = f_17 * smg_475[k]
                   + f_1 * snf0_396[k]
                   - f_2 * snf1_396[k]
                   + f_3 * pc_y[k] * sng_595[k];
    }

#pragma omp simd aligned(t_835, t_836, t_837, pc_y, pc_z, smg_460, smg_477, smg_478, snf0_398, \
                         snf0_399, snf1_398, snf1_399, sng_595, sng_597, \
                         sng_598 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_835[k] = f_11 * smg_460[k]
                   + f_3 * pc_z[k] * sng_595[k];

        t_836[k] = f_17 * smg_477[k]
                   + f_4 * snf0_398[k]
                   - f_5 * snf1_398[k]
                   + f_3 * pc_y[k] * sng_597[k];

        t_837[k] = f_17 * smg_478[k]
                   + f_6 * snf0_399[k]
                   - f_7 * snf1_399[k]
                   + f_3 * pc_y[k] * sng_598[k];
    }

#pragma omp simd aligned(t_838, t_839, t_840, pc_x, pc_y, pc_z, smg_464, smg_479, smg_600, \
                         snf0_399, snf0_400, snf1_399, snf1_400, sng_599, \
                         sng_600 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_838[k] = f_17 * smg_479[k]
                   + f_3 * pc_y[k] * sng_599[k];

        t_839[k] = f_11 * smg_464[k]
                   + f_1 * snf0_399[k]
                   - f_2 * snf1_399[k]
                   + f_3 * pc_z[k] * sng_599[k];

        t_840[k] = f_10 * smg_600[k]
                   + f_1 * snf0_400[k]
                   - f_2 * snf1_400[k]
                   + f_3 * pc_x[k] * sng_600[k];
    }

#pragma omp simd aligned(t_841, t_842, t_843, t_844, pc_x, pc_y, pc_z, smg_465, smg_480, \
                         smg_482, smg_603, snf0_403, snf1_403, sng_600, sng_602, \
                         sng_603 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_841[k] = f_16 * smg_480[k]
                   + f_3 * pc_y[k] * sng_600[k];

        t_842[k] = f_16 * smg_465[k]
                   + f_3 * pc_z[k] * sng_600[k];

        t_843[k] = f_10 * smg_603[k]
                   + f_4 * snf0_403[k]
                   - f_5 * snf1_403[k]
                   + f_3 * pc_x[k] * sng_603[k];

        t_844[k] = f_16 * smg_482[k]
                   + f_3 * pc_y[k] * sng_602[k];
    }

#pragma omp simd aligned(t_845, t_846, t_847, pc_x, pc_z, smg_468, smg_605, smg_606, snf0_405, \
                         snf0_406, snf1_405, snf1_406, sng_603, sng_605, \
                         sng_606 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_845[k] = f_10 * smg_605[k]
                   + f_4 * snf0_405[k]
                   - f_5 * snf1_405[k]
                   + f_3 * pc_x[k] * sng_605[k];

        t_846[k] = f_10 * smg_606[k]
                   + f_6 * snf0_406[k]
                   - f_7 * snf1_406[k]
                   + f_3 * pc_x[k] * sng_606[k];

        t_847[k] = f_16 * smg_468[k]
                   + f_3 * pc_z[k] * sng_603[k];
    }

#pragma omp simd aligned(t_848, t_849, t_850, t_851, pc_x, pc_y, smg_485, smg_609, smg_610, \
                         smg_611, snf0_409, snf1_409, sng_605, sng_609, sng_610, \
                         sng_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_848[k] = f_16 * smg_485[k]
                   + f_3 * pc_y[k] * sng_605[k];

        t_849[k] = f_10 * smg_609[k]
                   + f_6 * snf0_409[k]
                   - f_7 * snf1_409[k]
                   + f_3 * pc_x[k] * sng_609[k];

        t_850[k] = f_10 * smg_610[k]
                   + f_3 * pc_x[k] * sng_610[k];

        t_851[k] = f_10 * smg_611[k]
                   + f_3 * pc_x[k] * sng_611[k];
    }

#pragma omp simd aligned(t_852, t_853, t_854, t_855, pc_x, pc_y, smg_490, smg_612, smg_613, \
                         smg_614, snf0_406, snf1_406, sng_610, sng_612, sng_613, \
                         sng_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_852[k] = f_10 * smg_612[k]
                   + f_3 * pc_x[k] * sng_612[k];

        t_853[k] = f_10 * smg_613[k]
                   + f_3 * pc_x[k] * sng_613[k];

        t_854[k] = f_10 * smg_614[k]
                   + f_3 * pc_x[k] * sng_614[k];

        t_855[k] = f_16 * smg_490[k]
                   + f_1 * snf0_406[k]
                   - f_2 * snf1_406[k]
                   + f_3 * pc_y[k] * sng_610[k];
    }

#pragma omp simd aligned(t_856, t_857, t_858, pc_y, pc_z, smg_475, smg_492, smg_493, snf0_408, \
                         snf0_409, snf1_408, snf1_409, sng_610, sng_612, \
                         sng_613 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_856[k] = f_16 * smg_475[k]
                   + f_3 * pc_z[k] * sng_610[k];

        t_857[k] = f_16 * smg_492[k]
                   + f_4 * snf0_408[k]
                   - f_5 * snf1_408[k]
                   + f_3 * pc_y[k] * sng_612[k];

        t_858[k] = f_16 * smg_493[k]
                   + f_6 * snf0_409[k]
                   - f_7 * snf1_409[k]
                   + f_3 * pc_y[k] * sng_613[k];
    }

#pragma omp simd aligned(t_859, t_860, t_861, pc_x, pc_y, pc_z, smg_479, smg_494, smg_615, \
                         snf0_409, snf0_410, snf1_409, snf1_410, sng_614, \
                         sng_615 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_859[k] = f_16 * smg_494[k]
                   + f_3 * pc_y[k] * sng_614[k];

        t_860[k] = f_16 * smg_479[k]
                   + f_1 * snf0_409[k]
                   - f_2 * snf1_409[k]
                   + f_3 * pc_z[k] * sng_614[k];

        t_861[k] = f_10 * smg_615[k]
                   + f_1 * snf0_410[k]
                   - f_2 * snf1_410[k]
                   + f_3 * pc_x[k] * sng_615[k];
    }

#pragma omp simd aligned(t_862, t_863, t_864, t_865, pc_x, pc_y, pc_z, smg_480, smg_495, \
                         smg_497, smg_618, snf0_413, snf1_413, sng_615, sng_617, \
                         sng_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_862[k] = f_11 * smg_495[k]
                   + f_3 * pc_y[k] * sng_615[k];

        t_863[k] = f_17 * smg_480[k]
                   + f_3 * pc_z[k] * sng_615[k];

        t_864[k] = f_10 * smg_618[k]
                   + f_4 * snf0_413[k]
                   - f_5 * snf1_413[k]
                   + f_3 * pc_x[k] * sng_618[k];

        t_865[k] = f_11 * smg_497[k]
                   + f_3 * pc_y[k] * sng_617[k];
    }

#pragma omp simd aligned(t_866, t_867, t_868, pc_x, pc_z, smg_483, smg_620, smg_621, snf0_415, \
                         snf0_416, snf1_415, snf1_416, sng_618, sng_620, \
                         sng_621 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_866[k] = f_10 * smg_620[k]
                   + f_4 * snf0_415[k]
                   - f_5 * snf1_415[k]
                   + f_3 * pc_x[k] * sng_620[k];

        t_867[k] = f_10 * smg_621[k]
                   + f_6 * snf0_416[k]
                   - f_7 * snf1_416[k]
                   + f_3 * pc_x[k] * sng_621[k];

        t_868[k] = f_17 * smg_483[k]
                   + f_3 * pc_z[k] * sng_618[k];
    }

#pragma omp simd aligned(t_869, t_870, t_871, t_872, pc_x, pc_y, smg_500, smg_624, smg_625, \
                         smg_626, snf0_419, snf1_419, sng_620, sng_624, sng_625, \
                         sng_626 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_869[k] = f_11 * smg_500[k]
                   + f_3 * pc_y[k] * sng_620[k];

        t_870[k] = f_10 * smg_624[k]
                   + f_6 * snf0_419[k]
                   - f_7 * snf1_419[k]
                   + f_3 * pc_x[k] * sng_624[k];

        t_871[k] = f_10 * smg_625[k]
                   + f_3 * pc_x[k] * sng_625[k];

        t_872[k] = f_10 * smg_626[k]
                   + f_3 * pc_x[k] * sng_626[k];
    }

#pragma omp simd aligned(t_873, t_874, t_875, t_876, pc_x, pc_y, smg_505, smg_627, smg_628, \
                         smg_629, snf0_416, snf1_416, sng_625, sng_627, sng_628, \
                         sng_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_873[k] = f_10 * smg_627[k]
                   + f_3 * pc_x[k] * sng_627[k];

        t_874[k] = f_10 * smg_628[k]
                   + f_3 * pc_x[k] * sng_628[k];

        t_875[k] = f_10 * smg_629[k]
                   + f_3 * pc_x[k] * sng_629[k];

        t_876[k] = f_11 * smg_505[k]
                   + f_1 * snf0_416[k]
                   - f_2 * snf1_416[k]
                   + f_3 * pc_y[k] * sng_625[k];
    }

#pragma omp simd aligned(t_877, t_878, t_879, pc_y, pc_z, smg_490, smg_507, smg_508, snf0_418, \
                         snf0_419, snf1_418, snf1_419, sng_625, sng_627, \
                         sng_628 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_877[k] = f_17 * smg_490[k]
                   + f_3 * pc_z[k] * sng_625[k];

        t_878[k] = f_11 * smg_507[k]
                   + f_4 * snf0_418[k]
                   - f_5 * snf1_418[k]
                   + f_3 * pc_y[k] * sng_627[k];

        t_879[k] = f_11 * smg_508[k]
                   + f_6 * snf0_419[k]
                   - f_7 * snf1_419[k]
                   + f_3 * pc_y[k] * sng_628[k];
    }

#pragma omp simd aligned(t_880, t_881, t_882, pc_x, pc_y, pc_z, smg_494, smg_509, smg_630, \
                         snf0_419, snf0_420, snf1_419, snf1_420, sng_629, \
                         sng_630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_880[k] = f_11 * smg_509[k]
                   + f_3 * pc_y[k] * sng_629[k];

        t_881[k] = f_17 * smg_494[k]
                   + f_1 * snf0_419[k]
                   - f_2 * snf1_419[k]
                   + f_3 * pc_z[k] * sng_629[k];

        t_882[k] = f_10 * smg_630[k]
                   + f_1 * snf0_420[k]
                   - f_2 * snf1_420[k]
                   + f_3 * pc_x[k] * sng_630[k];
    }

#pragma omp simd aligned(t_883, t_884, t_885, t_886, pc_x, pc_y, pc_z, smg_495, smg_510, \
                         smg_512, smg_633, snf0_423, snf1_423, sng_630, sng_632, \
                         sng_633 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_883[k] = f_10 * smg_510[k]
                   + f_3 * pc_y[k] * sng_630[k];

        t_884[k] = f_15 * smg_495[k]
                   + f_3 * pc_z[k] * sng_630[k];

        t_885[k] = f_10 * smg_633[k]
                   + f_4 * snf0_423[k]
                   - f_5 * snf1_423[k]
                   + f_3 * pc_x[k] * sng_633[k];

        t_886[k] = f_10 * smg_512[k]
                   + f_3 * pc_y[k] * sng_632[k];
    }

#pragma omp simd aligned(t_887, t_888, t_889, pc_x, pc_z, smg_498, smg_635, smg_636, snf0_425, \
                         snf0_426, snf1_425, snf1_426, sng_633, sng_635, \
                         sng_636 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_887[k] = f_10 * smg_635[k]
                   + f_4 * snf0_425[k]
                   - f_5 * snf1_425[k]
                   + f_3 * pc_x[k] * sng_635[k];

        t_888[k] = f_10 * smg_636[k]
                   + f_6 * snf0_426[k]
                   - f_7 * snf1_426[k]
                   + f_3 * pc_x[k] * sng_636[k];

        t_889[k] = f_15 * smg_498[k]
                   + f_3 * pc_z[k] * sng_633[k];
    }

#pragma omp simd aligned(t_890, t_891, t_892, t_893, pc_x, pc_y, smg_515, smg_639, smg_640, \
                         smg_641, snf0_429, snf1_429, sng_635, sng_639, sng_640, \
                         sng_641 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_890[k] = f_10 * smg_515[k]
                   + f_3 * pc_y[k] * sng_635[k];

        t_891[k] = f_10 * smg_639[k]
                   + f_6 * snf0_429[k]
                   - f_7 * snf1_429[k]
                   + f_3 * pc_x[k] * sng_639[k];

        t_892[k] = f_10 * smg_640[k]
                   + f_3 * pc_x[k] * sng_640[k];

        t_893[k] = f_10 * smg_641[k]
                   + f_3 * pc_x[k] * sng_641[k];
    }

#pragma omp simd aligned(t_894, t_895, t_896, t_897, pc_x, pc_y, smg_520, smg_642, smg_643, \
                         smg_644, snf0_426, snf1_426, sng_640, sng_642, sng_643, \
                         sng_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_894[k] = f_10 * smg_642[k]
                   + f_3 * pc_x[k] * sng_642[k];

        t_895[k] = f_10 * smg_643[k]
                   + f_3 * pc_x[k] * sng_643[k];

        t_896[k] = f_10 * smg_644[k]
                   + f_3 * pc_x[k] * sng_644[k];

        t_897[k] = f_10 * smg_520[k]
                   + f_1 * snf0_426[k]
                   - f_2 * snf1_426[k]
                   + f_3 * pc_y[k] * sng_640[k];
    }

#pragma omp simd aligned(t_898, t_899, t_900, pc_y, pc_z, smg_505, smg_522, smg_523, snf0_428, \
                         snf0_429, snf1_428, snf1_429, sng_640, sng_642, \
                         sng_643 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_898[k] = f_15 * smg_505[k]
                   + f_3 * pc_z[k] * sng_640[k];

        t_899[k] = f_10 * smg_522[k]
                   + f_4 * snf0_428[k]
                   - f_5 * snf1_428[k]
                   + f_3 * pc_y[k] * sng_642[k];

        t_900[k] = f_10 * smg_523[k]
                   + f_6 * snf0_429[k]
                   - f_7 * snf1_429[k]
                   + f_3 * pc_y[k] * sng_643[k];
    }

#pragma omp simd aligned(t_901, t_902, t_903, t_904, pb_y, pc_y, pc_z, smh0_735, smg_509, \
                         smg_524, smg_525, smh1_735, snf0_429, snf1_429, sng_644, \
                         sng_645 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_901[k] = f_10 * smg_524[k]
                   + f_3 * pc_y[k] * sng_644[k];

        t_902[k] = f_15 * smg_509[k]
                   + f_1 * snf0_429[k]
                   - f_2 * snf1_429[k]
                   + f_3 * pc_z[k] * sng_644[k];

        t_903[k] = pb_y[k] * smh0_735[k]
                   - f_8 * pc_y[k] * smh1_735[k];

        t_904[k] = f_9 * smg_525[k]
                   + f_3 * pc_y[k] * sng_645[k];
    }

#pragma omp simd aligned(t_905, t_906, t_907, t_908, pb_y, pc_y, pc_z, smh0_738, smh0_740, \
                         smg_510, smg_526, smg_527, smh1_738, smh1_740, sng_645, \
                         sng_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_905[k] = f_14 * smg_510[k]
                   + f_3 * pc_z[k] * sng_645[k];

        t_906[k] = pb_y[k] * smh0_738[k]
                   + f_10 * smg_526[k]
                   - f_8 * pc_y[k] * smh1_738[k];

        t_907[k] = f_9 * smg_527[k]
                   + f_3 * pc_y[k] * sng_647[k];

        t_908[k] = pb_y[k] * smh0_740[k]
                   - f_8 * pc_y[k] * smh1_740[k];
    }

#pragma omp simd aligned(t_909, t_910, t_911, t_912, pb_y, pc_y, pc_z, smh0_741, smh0_744, \
                         smg_513, smg_528, smg_530, smh1_741, smh1_744, sng_648, \
                         sng_650 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_909[k] = pb_y[k] * smh0_741[k]
                   + f_11 * smg_528[k]
                   - f_8 * pc_y[k] * smh1_741[k];

        t_910[k] = f_14 * smg_513[k]
                   + f_3 * pc_z[k] * sng_648[k];

        t_911[k] = f_9 * smg_530[k]
                   + f_3 * pc_y[k] * sng_650[k];

        t_912[k] = pb_y[k] * smh0_744[k]
                   - f_8 * pc_y[k] * smh1_744[k];
    }

#pragma omp simd aligned(t_913, t_914, t_915, t_916, t_917, pc_x, smg_655, smg_656, smg_657, \
                         smg_658, smg_659, sng_655, sng_656, sng_657, sng_658, \
                         sng_659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_913[k] = f_10 * smg_655[k]
                   + f_3 * pc_x[k] * sng_655[k];

        t_914[k] = f_10 * smg_656[k]
                   + f_3 * pc_x[k] * sng_656[k];

        t_915[k] = f_10 * smg_657[k]
                   + f_3 * pc_x[k] * sng_657[k];

        t_916[k] = f_10 * smg_658[k]
                   + f_3 * pc_x[k] * sng_658[k];

        t_917[k] = f_10 * smg_659[k]
                   + f_3 * pc_x[k] * sng_659[k];
    }

#pragma omp simd aligned(t_918, t_919, t_920, pc_y, pc_z, smg_520, smg_535, smg_537, snf0_436, \
                         snf0_438, snf1_436, snf1_438, sng_655, \
                         sng_657 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_918[k] = f_9 * smg_535[k]
                   + f_1 * snf0_436[k]
                   - f_2 * snf1_436[k]
                   + f_3 * pc_y[k] * sng_655[k];

        t_919[k] = f_14 * smg_520[k]
                   + f_3 * pc_z[k] * sng_655[k];

        t_920[k] = f_9 * smg_537[k]
                   + f_4 * snf0_438[k]
                   - f_5 * snf1_438[k]
                   + f_3 * pc_y[k] * sng_657[k];
    }
}

static auto
compute_prim_snh_three_center_electron_repulsion_0_piece8(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t smh0,
                                                          const size_t smg, const size_t smh1,
                                                          const size_t snf0, const size_t snf1,
                                                          const size_t sng, const size_t ncols,
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
    const auto f_12 = 4.5 / q;
    const auto f_13 = 4.0 / q;
    const auto f_14 = 3.5 / q;
    const auto f_15 = 3.0 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;

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
    auto *t_1037 = buffer.data(target + 1037);
    auto *t_1038 = buffer.data(target + 1038);
    auto *t_1039 = buffer.data(target + 1039);
    auto *t_1040 = buffer.data(target + 1040);
    auto *t_1041 = buffer.data(target + 1041);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *smh0_755 = buffer.data(smh0 + 755);
    const auto *smh0_756 = buffer.data(smh0 + 756);
    const auto *smh0_759 = buffer.data(smh0 + 759);
    const auto *smh0_762 = buffer.data(smh0 + 762);
    const auto *smh0_945 = buffer.data(smh0 + 945);
    const auto *smh0_948 = buffer.data(smh0 + 948);
    const auto *smh0_950 = buffer.data(smh0 + 950);
    const auto *smh0_951 = buffer.data(smh0 + 951);
    const auto *smh0_954 = buffer.data(smh0 + 954);
    const auto *smh0_960 = buffer.data(smh0 + 960);
    const auto *smh0_962 = buffer.data(smh0 + 962);
    const auto *smh0_963 = buffer.data(smh0 + 963);
    const auto *smh0_965 = buffer.data(smh0 + 965);
    const auto *smh0_971 = buffer.data(smh0 + 971);
    const auto *smh0_975 = buffer.data(smh0 + 975);
    const auto *smh0_981 = buffer.data(smh0 + 981);
    const auto *smh0_983 = buffer.data(smh0 + 983);
    const auto *smh0_984 = buffer.data(smh0 + 984);
    const auto *smh0_986 = buffer.data(smh0 + 986);
    const auto *smh0_987 = buffer.data(smh0 + 987);
    const auto *smh0_990 = buffer.data(smh0 + 990);
    const auto *smh0_992 = buffer.data(smh0 + 992);
    const auto *smh0_993 = buffer.data(smh0 + 993);
    const auto *smh0_996 = buffer.data(smh0 + 996);
    const auto *smh0_1002 = buffer.data(smh0 + 1002);
    const auto *smh0_1004 = buffer.data(smh0 + 1004);
    const auto *smh0_1005 = buffer.data(smh0 + 1005);
    const auto *smh0_1007 = buffer.data(smh0 + 1007);
    const auto *smh0_1008 = buffer.data(smh0 + 1008);
    const auto *smh0_1011 = buffer.data(smh0 + 1011);
    const auto *smh0_1013 = buffer.data(smh0 + 1013);
    const auto *smh0_1014 = buffer.data(smh0 + 1014);
    const auto *smh0_1017 = buffer.data(smh0 + 1017);
    const auto *smh0_1023 = buffer.data(smh0 + 1023);
    const auto *smh0_1025 = buffer.data(smh0 + 1025);
    const auto *smh0_1026 = buffer.data(smh0 + 1026);
    const auto *smh0_1028 = buffer.data(smh0 + 1028);
    const auto *smh0_1029 = buffer.data(smh0 + 1029);
    const auto *smh0_1032 = buffer.data(smh0 + 1032);
    const auto *smh0_1034 = buffer.data(smh0 + 1034);
    const auto *smh0_1035 = buffer.data(smh0 + 1035);
    const auto *smh0_1038 = buffer.data(smh0 + 1038);

    const auto *smg_525 = buffer.data(smg + 525);
    const auto *smg_528 = buffer.data(smg + 528);
    const auto *smg_535 = buffer.data(smg + 535);
    const auto *smg_538 = buffer.data(smg + 538);
    const auto *smg_539 = buffer.data(smg + 539);
    const auto *smg_540 = buffer.data(smg + 540);
    const auto *smg_542 = buffer.data(smg + 542);
    const auto *smg_543 = buffer.data(smg + 543);
    const auto *smg_545 = buffer.data(smg + 545);
    const auto *smg_550 = buffer.data(smg + 550);
    const auto *smg_554 = buffer.data(smg + 554);
    const auto *smg_555 = buffer.data(smg + 555);
    const auto *smg_557 = buffer.data(smg + 557);
    const auto *smg_558 = buffer.data(smg + 558);
    const auto *smg_560 = buffer.data(smg + 560);
    const auto *smg_565 = buffer.data(smg + 565);
    const auto *smg_569 = buffer.data(smg + 569);
    const auto *smg_570 = buffer.data(smg + 570);
    const auto *smg_572 = buffer.data(smg + 572);
    const auto *smg_573 = buffer.data(smg + 573);
    const auto *smg_575 = buffer.data(smg + 575);
    const auto *smg_580 = buffer.data(smg + 580);
    const auto *smg_584 = buffer.data(smg + 584);
    const auto *smg_585 = buffer.data(smg + 585);
    const auto *smg_587 = buffer.data(smg + 587);
    const auto *smg_588 = buffer.data(smg + 588);
    const auto *smg_590 = buffer.data(smg + 590);
    const auto *smg_599 = buffer.data(smg + 599);
    const auto *smg_600 = buffer.data(smg + 600);
    const auto *smg_602 = buffer.data(smg + 602);
    const auto *smg_605 = buffer.data(smg + 605);
    const auto *smg_660 = buffer.data(smg + 660);
    const auto *smg_663 = buffer.data(smg + 663);
    const auto *smg_665 = buffer.data(smg + 665);
    const auto *smg_666 = buffer.data(smg + 666);
    const auto *smg_669 = buffer.data(smg + 669);
    const auto *smg_670 = buffer.data(smg + 670);
    const auto *smg_671 = buffer.data(smg + 671);
    const auto *smg_672 = buffer.data(smg + 672);
    const auto *smg_673 = buffer.data(smg + 673);
    const auto *smg_674 = buffer.data(smg + 674);
    const auto *smg_675 = buffer.data(smg + 675);
    const auto *smg_678 = buffer.data(smg + 678);
    const auto *smg_680 = buffer.data(smg + 680);
    const auto *smg_681 = buffer.data(smg + 681);
    const auto *smg_684 = buffer.data(smg + 684);
    const auto *smg_685 = buffer.data(smg + 685);
    const auto *smg_686 = buffer.data(smg + 686);
    const auto *smg_687 = buffer.data(smg + 687);
    const auto *smg_688 = buffer.data(smg + 688);
    const auto *smg_689 = buffer.data(smg + 689);
    const auto *smg_695 = buffer.data(smg + 695);
    const auto *smg_699 = buffer.data(smg + 699);
    const auto *smg_700 = buffer.data(smg + 700);
    const auto *smg_701 = buffer.data(smg + 701);
    const auto *smg_702 = buffer.data(smg + 702);
    const auto *smg_703 = buffer.data(smg + 703);
    const auto *smg_704 = buffer.data(smg + 704);
    const auto *smg_705 = buffer.data(smg + 705);
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
    const auto *smg_738 = buffer.data(smg + 738);
    const auto *smg_740 = buffer.data(smg + 740);
    const auto *smg_741 = buffer.data(smg + 741);
    const auto *smg_744 = buffer.data(smg + 744);
    const auto *smg_745 = buffer.data(smg + 745);
    const auto *smg_746 = buffer.data(smg + 746);
    const auto *smg_747 = buffer.data(smg + 747);

    const auto *smh1_755 = buffer.data(smh1 + 755);
    const auto *smh1_756 = buffer.data(smh1 + 756);
    const auto *smh1_759 = buffer.data(smh1 + 759);
    const auto *smh1_762 = buffer.data(smh1 + 762);
    const auto *smh1_945 = buffer.data(smh1 + 945);
    const auto *smh1_948 = buffer.data(smh1 + 948);
    const auto *smh1_950 = buffer.data(smh1 + 950);
    const auto *smh1_951 = buffer.data(smh1 + 951);
    const auto *smh1_954 = buffer.data(smh1 + 954);
    const auto *smh1_960 = buffer.data(smh1 + 960);
    const auto *smh1_962 = buffer.data(smh1 + 962);
    const auto *smh1_963 = buffer.data(smh1 + 963);
    const auto *smh1_965 = buffer.data(smh1 + 965);
    const auto *smh1_971 = buffer.data(smh1 + 971);
    const auto *smh1_975 = buffer.data(smh1 + 975);
    const auto *smh1_981 = buffer.data(smh1 + 981);
    const auto *smh1_983 = buffer.data(smh1 + 983);
    const auto *smh1_984 = buffer.data(smh1 + 984);
    const auto *smh1_986 = buffer.data(smh1 + 986);
    const auto *smh1_987 = buffer.data(smh1 + 987);
    const auto *smh1_990 = buffer.data(smh1 + 990);
    const auto *smh1_992 = buffer.data(smh1 + 992);
    const auto *smh1_993 = buffer.data(smh1 + 993);
    const auto *smh1_996 = buffer.data(smh1 + 996);
    const auto *smh1_1002 = buffer.data(smh1 + 1002);
    const auto *smh1_1004 = buffer.data(smh1 + 1004);
    const auto *smh1_1005 = buffer.data(smh1 + 1005);
    const auto *smh1_1007 = buffer.data(smh1 + 1007);
    const auto *smh1_1008 = buffer.data(smh1 + 1008);
    const auto *smh1_1011 = buffer.data(smh1 + 1011);
    const auto *smh1_1013 = buffer.data(smh1 + 1013);
    const auto *smh1_1014 = buffer.data(smh1 + 1014);
    const auto *smh1_1017 = buffer.data(smh1 + 1017);
    const auto *smh1_1023 = buffer.data(smh1 + 1023);
    const auto *smh1_1025 = buffer.data(smh1 + 1025);
    const auto *smh1_1026 = buffer.data(smh1 + 1026);
    const auto *smh1_1028 = buffer.data(smh1 + 1028);
    const auto *smh1_1029 = buffer.data(smh1 + 1029);
    const auto *smh1_1032 = buffer.data(smh1 + 1032);
    const auto *smh1_1034 = buffer.data(smh1 + 1034);
    const auto *smh1_1035 = buffer.data(smh1 + 1035);
    const auto *smh1_1038 = buffer.data(smh1 + 1038);

    const auto *snf0_439 = buffer.data(snf0 + 439);
    const auto *snf0_440 = buffer.data(snf0 + 440);
    const auto *snf0_443 = buffer.data(snf0 + 443);
    const auto *snf0_445 = buffer.data(snf0 + 445);
    const auto *snf0_446 = buffer.data(snf0 + 446);
    const auto *snf0_448 = buffer.data(snf0 + 448);
    const auto *snf0_449 = buffer.data(snf0 + 449);

    const auto *snf1_439 = buffer.data(snf1 + 439);
    const auto *snf1_440 = buffer.data(snf1 + 440);
    const auto *snf1_443 = buffer.data(snf1 + 443);
    const auto *snf1_445 = buffer.data(snf1 + 445);
    const auto *snf1_446 = buffer.data(snf1 + 446);
    const auto *snf1_448 = buffer.data(snf1 + 448);
    const auto *snf1_449 = buffer.data(snf1 + 449);

    const auto *sng_658 = buffer.data(sng + 658);
    const auto *sng_659 = buffer.data(sng + 659);
    const auto *sng_660 = buffer.data(sng + 660);
    const auto *sng_662 = buffer.data(sng + 662);
    const auto *sng_663 = buffer.data(sng + 663);
    const auto *sng_665 = buffer.data(sng + 665);
    const auto *sng_666 = buffer.data(sng + 666);
    const auto *sng_669 = buffer.data(sng + 669);
    const auto *sng_670 = buffer.data(sng + 670);
    const auto *sng_671 = buffer.data(sng + 671);
    const auto *sng_672 = buffer.data(sng + 672);
    const auto *sng_673 = buffer.data(sng + 673);
    const auto *sng_674 = buffer.data(sng + 674);
    const auto *sng_675 = buffer.data(sng + 675);
    const auto *sng_677 = buffer.data(sng + 677);
    const auto *sng_678 = buffer.data(sng + 678);
    const auto *sng_680 = buffer.data(sng + 680);
    const auto *sng_685 = buffer.data(sng + 685);
    const auto *sng_686 = buffer.data(sng + 686);
    const auto *sng_687 = buffer.data(sng + 687);
    const auto *sng_688 = buffer.data(sng + 688);
    const auto *sng_689 = buffer.data(sng + 689);
    const auto *sng_690 = buffer.data(sng + 690);
    const auto *sng_692 = buffer.data(sng + 692);
    const auto *sng_693 = buffer.data(sng + 693);
    const auto *sng_695 = buffer.data(sng + 695);
    const auto *sng_700 = buffer.data(sng + 700);
    const auto *sng_701 = buffer.data(sng + 701);
    const auto *sng_702 = buffer.data(sng + 702);
    const auto *sng_703 = buffer.data(sng + 703);
    const auto *sng_704 = buffer.data(sng + 704);
    const auto *sng_705 = buffer.data(sng + 705);
    const auto *sng_707 = buffer.data(sng + 707);
    const auto *sng_708 = buffer.data(sng + 708);
    const auto *sng_710 = buffer.data(sng + 710);
    const auto *sng_715 = buffer.data(sng + 715);
    const auto *sng_716 = buffer.data(sng + 716);
    const auto *sng_717 = buffer.data(sng + 717);
    const auto *sng_718 = buffer.data(sng + 718);
    const auto *sng_719 = buffer.data(sng + 719);
    const auto *sng_720 = buffer.data(sng + 720);
    const auto *sng_722 = buffer.data(sng + 722);
    const auto *sng_723 = buffer.data(sng + 723);
    const auto *sng_725 = buffer.data(sng + 725);
    const auto *sng_730 = buffer.data(sng + 730);
    const auto *sng_731 = buffer.data(sng + 731);
    const auto *sng_732 = buffer.data(sng + 732);
    const auto *sng_733 = buffer.data(sng + 733);
    const auto *sng_734 = buffer.data(sng + 734);
    const auto *sng_735 = buffer.data(sng + 735);
    const auto *sng_737 = buffer.data(sng + 737);
    const auto *sng_738 = buffer.data(sng + 738);
    const auto *sng_740 = buffer.data(sng + 740);
    const auto *sng_745 = buffer.data(sng + 745);
    const auto *sng_746 = buffer.data(sng + 746);
    const auto *sng_747 = buffer.data(sng + 747);

#pragma omp simd aligned(t_921, t_922, t_923, pb_y, pc_y, smh0_755, smg_538, smg_539, \
                         smh1_755, snf0_439, snf1_439, sng_658, \
                         sng_659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_921[k] = f_9 * smg_538[k]
                   + f_6 * snf0_439[k]
                   - f_7 * snf1_439[k]
                   + f_3 * pc_y[k] * sng_658[k];

        t_922[k] = f_9 * smg_539[k]
                   + f_3 * pc_y[k] * sng_659[k];

        t_923[k] = pb_y[k] * smh0_755[k]
                   - f_8 * pc_y[k] * smh1_755[k];
    }

#pragma omp simd aligned(t_924, t_925, t_926, t_927, pc_x, pc_y, pc_z, smg_525, smg_660, \
                         smg_663, snf0_440, snf0_443, snf1_440, snf1_443, sng_660, \
                         sng_663 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_924[k] = f_10 * smg_660[k]
                   + f_1 * snf0_440[k]
                   - f_2 * snf1_440[k]
                   + f_3 * pc_x[k] * sng_660[k];

        t_925[k] = f_3 * pc_y[k] * sng_660[k];

        t_926[k] = f_13 * smg_525[k]
                   + f_3 * pc_z[k] * sng_660[k];

        t_927[k] = f_10 * smg_663[k]
                   + f_4 * snf0_443[k]
                   - f_5 * snf1_443[k]
                   + f_3 * pc_x[k] * sng_663[k];
    }

#pragma omp simd aligned(t_928, t_929, t_930, pc_x, pc_y, smg_665, smg_666, snf0_445, \
                         snf0_446, snf1_445, snf1_446, sng_662, sng_665, \
                         sng_666 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_928[k] = f_3 * pc_y[k] * sng_662[k];

        t_929[k] = f_10 * smg_665[k]
                   + f_4 * snf0_445[k]
                   - f_5 * snf1_445[k]
                   + f_3 * pc_x[k] * sng_665[k];

        t_930[k] = f_10 * smg_666[k]
                   + f_6 * snf0_446[k]
                   - f_7 * snf1_446[k]
                   + f_3 * pc_x[k] * sng_666[k];
    }

#pragma omp simd aligned(t_931, t_932, t_933, t_934, pc_x, pc_y, pc_z, smg_528, smg_669, \
                         smg_670, snf0_449, snf1_449, sng_663, sng_665, sng_669, \
                         sng_670 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_931[k] = f_13 * smg_528[k]
                   + f_3 * pc_z[k] * sng_663[k];

        t_932[k] = f_3 * pc_y[k] * sng_665[k];

        t_933[k] = f_10 * smg_669[k]
                   + f_6 * snf0_449[k]
                   - f_7 * snf1_449[k]
                   + f_3 * pc_x[k] * sng_669[k];

        t_934[k] = f_10 * smg_670[k]
                   + f_3 * pc_x[k] * sng_670[k];
    }

#pragma omp simd aligned(t_935, t_936, t_937, t_938, pc_x, smg_671, smg_672, smg_673, smg_674, \
                         sng_671, sng_672, sng_673, sng_674 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_935[k] = f_10 * smg_671[k]
                   + f_3 * pc_x[k] * sng_671[k];

        t_936[k] = f_10 * smg_672[k]
                   + f_3 * pc_x[k] * sng_672[k];

        t_937[k] = f_10 * smg_673[k]
                   + f_3 * pc_x[k] * sng_673[k];

        t_938[k] = f_10 * smg_674[k]
                   + f_3 * pc_x[k] * sng_674[k];
    }

#pragma omp simd aligned(t_939, t_940, t_941, t_942, pc_y, pc_z, smg_535, snf0_446, snf0_448, \
                         snf0_449, snf1_446, snf1_448, snf1_449, sng_670, sng_672, \
                         sng_673 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_939[k] = f_1 * snf0_446[k]
                   - f_2 * snf1_446[k]
                   + f_3 * pc_y[k] * sng_670[k];

        t_940[k] = f_13 * smg_535[k]
                   + f_3 * pc_z[k] * sng_670[k];

        t_941[k] = f_4 * snf0_448[k]
                   - f_5 * snf1_448[k]
                   + f_3 * pc_y[k] * sng_672[k];

        t_942[k] = f_6 * snf0_449[k]
                   - f_7 * snf1_449[k]
                   + f_3 * pc_y[k] * sng_673[k];
    }

#pragma omp simd aligned(t_943, t_944, t_945, pb_x, pc_x, pc_y, pc_z, smh0_945, smg_539, \
                         smg_675, smh1_945, snf0_449, snf1_449, \
                         sng_674 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_943[k] = f_3 * pc_y[k] * sng_674[k];

        t_944[k] = f_13 * smg_539[k]
                   + f_1 * snf0_449[k]
                   - f_2 * snf1_449[k]
                   + f_3 * pc_z[k] * sng_674[k];

        t_945[k] = pb_x[k] * smh0_945[k]
                   + f_17 * smg_675[k]
                   - f_8 * pc_x[k] * smh1_945[k];
    }

#pragma omp simd aligned(t_946, t_947, t_948, t_949, pb_x, pc_x, pc_y, pc_z, smh0_948, \
                         smg_540, smg_542, smg_678, smh1_948, sng_675, \
                         sng_677 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_946[k] = f_12 * smg_540[k]
                   + f_3 * pc_y[k] * sng_675[k];

        t_947[k] = f_3 * pc_z[k] * sng_675[k];

        t_948[k] = pb_x[k] * smh0_948[k]
                   + f_11 * smg_678[k]
                   - f_8 * pc_x[k] * smh1_948[k];

        t_949[k] = f_12 * smg_542[k]
                   + f_3 * pc_y[k] * sng_677[k];
    }

#pragma omp simd aligned(t_950, t_951, t_952, pb_x, pc_x, pc_z, smh0_950, smh0_951, smg_680, \
                         smg_681, smh1_950, smh1_951, sng_678 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_950[k] = pb_x[k] * smh0_950[k]
                   + f_11 * smg_680[k]
                   - f_8 * pc_x[k] * smh1_950[k];

        t_951[k] = pb_x[k] * smh0_951[k]
                   + f_10 * smg_681[k]
                   - f_8 * pc_x[k] * smh1_951[k];

        t_952[k] = f_3 * pc_z[k] * sng_678[k];
    }

#pragma omp simd aligned(t_953, t_954, t_955, t_956, pb_x, pc_x, pc_y, smh0_954, smg_545, \
                         smg_684, smg_685, smg_686, smh1_954, sng_680, sng_685, \
                         sng_686 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_953[k] = f_12 * smg_545[k]
                   + f_3 * pc_y[k] * sng_680[k];

        t_954[k] = pb_x[k] * smh0_954[k]
                   + f_10 * smg_684[k]
                   - f_8 * pc_x[k] * smh1_954[k];

        t_955[k] = f_9 * smg_685[k]
                   + f_3 * pc_x[k] * sng_685[k];

        t_956[k] = f_9 * smg_686[k]
                   + f_3 * pc_x[k] * sng_686[k];
    }

#pragma omp simd aligned(t_957, t_958, t_959, t_960, pb_x, pc_x, smh0_960, smg_687, smg_688, \
                         smg_689, smh1_960, sng_687, sng_688, sng_689 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_957[k] = f_9 * smg_687[k]
                   + f_3 * pc_x[k] * sng_687[k];

        t_958[k] = f_9 * smg_688[k]
                   + f_3 * pc_x[k] * sng_688[k];

        t_959[k] = f_9 * smg_689[k]
                   + f_3 * pc_x[k] * sng_689[k];

        t_960[k] = pb_x[k] * smh0_960[k]
                   - f_8 * pc_x[k] * smh1_960[k];
    }

#pragma omp simd aligned(t_961, t_962, t_963, t_964, pb_x, pc_x, pc_y, pc_z, smh0_962, \
                         smh0_963, smg_554, smh1_962, smh1_963, sng_685, \
                         sng_689 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_961[k] = f_3 * pc_z[k] * sng_685[k];

        t_962[k] = pb_x[k] * smh0_962[k]
                   - f_8 * pc_x[k] * smh1_962[k];

        t_963[k] = pb_x[k] * smh0_963[k]
                   - f_8 * pc_x[k] * smh1_963[k];

        t_964[k] = f_12 * smg_554[k]
                   + f_3 * pc_y[k] * sng_689[k];
    }

#pragma omp simd aligned(t_965, t_966, t_967, t_968, pb_x, pb_z, pc_x, pc_y, pc_z, smh0_756, \
                         smh0_965, smg_540, smg_555, smh1_756, smh1_965, \
                         sng_690 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_965[k] = pb_x[k] * smh0_965[k]
                   - f_8 * pc_x[k] * smh1_965[k];

        t_966[k] = pb_z[k] * smh0_756[k]
                   - f_8 * pc_z[k] * smh1_756[k];

        t_967[k] = f_13 * smg_555[k]
                   + f_3 * pc_y[k] * sng_690[k];

        t_968[k] = f_9 * smg_540[k]
                   + f_3 * pc_z[k] * sng_690[k];
    }

#pragma omp simd aligned(t_969, t_970, t_971, pb_x, pb_z, pc_x, pc_y, pc_z, smh0_759, \
                         smh0_971, smg_557, smg_695, smh1_759, smh1_971, \
                         sng_692 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_969[k] = pb_z[k] * smh0_759[k]
                   - f_8 * pc_z[k] * smh1_759[k];

        t_970[k] = f_13 * smg_557[k]
                   + f_3 * pc_y[k] * sng_692[k];

        t_971[k] = pb_x[k] * smh0_971[k]
                   + f_11 * smg_695[k]
                   - f_8 * pc_x[k] * smh1_971[k];
    }

#pragma omp simd aligned(t_972, t_973, t_974, pb_z, pc_y, pc_z, smh0_762, smg_543, smg_560, \
                         smh1_762, sng_693, sng_695 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_972[k] = pb_z[k] * smh0_762[k]
                   - f_8 * pc_z[k] * smh1_762[k];

        t_973[k] = f_9 * smg_543[k]
                   + f_3 * pc_z[k] * sng_693[k];

        t_974[k] = f_13 * smg_560[k]
                   + f_3 * pc_y[k] * sng_695[k];
    }

#pragma omp simd aligned(t_975, t_976, t_977, t_978, pb_x, pc_x, smh0_975, smg_699, smg_700, \
                         smg_701, smg_702, smh1_975, sng_700, sng_701, \
                         sng_702 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_975[k] = pb_x[k] * smh0_975[k]
                   + f_10 * smg_699[k]
                   - f_8 * pc_x[k] * smh1_975[k];

        t_976[k] = f_9 * smg_700[k]
                   + f_3 * pc_x[k] * sng_700[k];

        t_977[k] = f_9 * smg_701[k]
                   + f_3 * pc_x[k] * sng_701[k];

        t_978[k] = f_9 * smg_702[k]
                   + f_3 * pc_x[k] * sng_702[k];
    }

#pragma omp simd aligned(t_979, t_980, t_981, t_982, pb_x, pc_x, pc_z, smh0_981, smg_550, \
                         smg_703, smg_704, smh1_981, sng_700, sng_703, \
                         sng_704 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_979[k] = f_9 * smg_703[k]
                   + f_3 * pc_x[k] * sng_703[k];

        t_980[k] = f_9 * smg_704[k]
                   + f_3 * pc_x[k] * sng_704[k];

        t_981[k] = pb_x[k] * smh0_981[k]
                   - f_8 * pc_x[k] * smh1_981[k];

        t_982[k] = f_9 * smg_550[k]
                   + f_3 * pc_z[k] * sng_700[k];
    }

#pragma omp simd aligned(t_983, t_984, t_985, t_986, pb_x, pc_x, pc_y, smh0_983, smh0_984, \
                         smh0_986, smg_569, smh1_983, smh1_984, smh1_986, \
                         sng_704 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_983[k] = pb_x[k] * smh0_983[k]
                   - f_8 * pc_x[k] * smh1_983[k];

        t_984[k] = pb_x[k] * smh0_984[k]
                   - f_8 * pc_x[k] * smh1_984[k];

        t_985[k] = f_13 * smg_569[k]
                   + f_3 * pc_y[k] * sng_704[k];

        t_986[k] = pb_x[k] * smh0_986[k]
                   - f_8 * pc_x[k] * smh1_986[k];
    }

#pragma omp simd aligned(t_987, t_988, t_989, pb_x, pc_x, pc_y, pc_z, smh0_987, smg_555, \
                         smg_570, smg_705, smh1_987, sng_705 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_987[k] = pb_x[k] * smh0_987[k]
                   + f_17 * smg_705[k]
                   - f_8 * pc_x[k] * smh1_987[k];

        t_988[k] = f_14 * smg_570[k]
                   + f_3 * pc_y[k] * sng_705[k];

        t_989[k] = f_10 * smg_555[k]
                   + f_3 * pc_z[k] * sng_705[k];
    }

#pragma omp simd aligned(t_990, t_991, t_992, pb_x, pc_x, pc_y, smh0_990, smh0_992, smg_572, \
                         smg_708, smg_710, smh1_990, smh1_992, \
                         sng_707 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_990[k] = pb_x[k] * smh0_990[k]
                   + f_11 * smg_708[k]
                   - f_8 * pc_x[k] * smh1_990[k];

        t_991[k] = f_14 * smg_572[k]
                   + f_3 * pc_y[k] * sng_707[k];

        t_992[k] = pb_x[k] * smh0_992[k]
                   + f_11 * smg_710[k]
                   - f_8 * pc_x[k] * smh1_992[k];
    }

#pragma omp simd aligned(t_993, t_994, t_995, pb_x, pc_x, pc_y, pc_z, smh0_993, smg_558, \
                         smg_575, smg_711, smh1_993, sng_708, sng_710 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_993[k] = pb_x[k] * smh0_993[k]
                   + f_10 * smg_711[k]
                   - f_8 * pc_x[k] * smh1_993[k];

        t_994[k] = f_10 * smg_558[k]
                   + f_3 * pc_z[k] * sng_708[k];

        t_995[k] = f_14 * smg_575[k]
                   + f_3 * pc_y[k] * sng_710[k];
    }

#pragma omp simd aligned(t_996, t_997, t_998, t_999, pb_x, pc_x, smh0_996, smg_714, smg_715, \
                         smg_716, smg_717, smh1_996, sng_715, sng_716, \
                         sng_717 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_996[k] = pb_x[k] * smh0_996[k]
                   + f_10 * smg_714[k]
                   - f_8 * pc_x[k] * smh1_996[k];

        t_997[k] = f_9 * smg_715[k]
                   + f_3 * pc_x[k] * sng_715[k];

        t_998[k] = f_9 * smg_716[k]
                   + f_3 * pc_x[k] * sng_716[k];

        t_999[k] = f_9 * smg_717[k]
                   + f_3 * pc_x[k] * sng_717[k];
    }

#pragma omp simd aligned(t_1000, t_1001, t_1002, t_1003, pb_x, pc_x, pc_z, smh0_1002, smg_565, \
                         smg_718, smg_719, smh1_1002, sng_715, sng_718, \
                         sng_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1000[k] = f_9 * smg_718[k]
                    + f_3 * pc_x[k] * sng_718[k];

        t_1001[k] = f_9 * smg_719[k]
                    + f_3 * pc_x[k] * sng_719[k];

        t_1002[k] = pb_x[k] * smh0_1002[k]
                    - f_8 * pc_x[k] * smh1_1002[k];

        t_1003[k] = f_10 * smg_565[k]
                    + f_3 * pc_z[k] * sng_715[k];
    }

#pragma omp simd aligned(t_1004, t_1005, t_1006, t_1007, pb_x, pc_x, pc_y, smh0_1004, \
                         smh0_1005, smh0_1007, smg_584, smh1_1004, smh1_1005, smh1_1007, \
                         sng_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1004[k] = pb_x[k] * smh0_1004[k]
                    - f_8 * pc_x[k] * smh1_1004[k];

        t_1005[k] = pb_x[k] * smh0_1005[k]
                    - f_8 * pc_x[k] * smh1_1005[k];

        t_1006[k] = f_14 * smg_584[k]
                    + f_3 * pc_y[k] * sng_719[k];

        t_1007[k] = pb_x[k] * smh0_1007[k]
                    - f_8 * pc_x[k] * smh1_1007[k];
    }

#pragma omp simd aligned(t_1008, t_1009, t_1010, pb_x, pc_x, pc_y, pc_z, smh0_1008, smg_570, \
                         smg_585, smg_720, smh1_1008, sng_720 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1008[k] = pb_x[k] * smh0_1008[k]
                    + f_17 * smg_720[k]
                    - f_8 * pc_x[k] * smh1_1008[k];

        t_1009[k] = f_15 * smg_585[k]
                    + f_3 * pc_y[k] * sng_720[k];

        t_1010[k] = f_11 * smg_570[k]
                    + f_3 * pc_z[k] * sng_720[k];
    }

#pragma omp simd aligned(t_1011, t_1012, t_1013, pb_x, pc_x, pc_y, smh0_1011, smh0_1013, \
                         smg_587, smg_723, smg_725, smh1_1011, smh1_1013, \
                         sng_722 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1011[k] = pb_x[k] * smh0_1011[k]
                    + f_11 * smg_723[k]
                    - f_8 * pc_x[k] * smh1_1011[k];

        t_1012[k] = f_15 * smg_587[k]
                    + f_3 * pc_y[k] * sng_722[k];

        t_1013[k] = pb_x[k] * smh0_1013[k]
                    + f_11 * smg_725[k]
                    - f_8 * pc_x[k] * smh1_1013[k];
    }

#pragma omp simd aligned(t_1014, t_1015, t_1016, pb_x, pc_x, pc_y, pc_z, smh0_1014, smg_573, \
                         smg_590, smg_726, smh1_1014, sng_723, \
                         sng_725 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1014[k] = pb_x[k] * smh0_1014[k]
                    + f_10 * smg_726[k]
                    - f_8 * pc_x[k] * smh1_1014[k];

        t_1015[k] = f_11 * smg_573[k]
                    + f_3 * pc_z[k] * sng_723[k];

        t_1016[k] = f_15 * smg_590[k]
                    + f_3 * pc_y[k] * sng_725[k];
    }

#pragma omp simd aligned(t_1017, t_1018, t_1019, t_1020, pb_x, pc_x, smh0_1017, smg_729, \
                         smg_730, smg_731, smg_732, smh1_1017, sng_730, sng_731, \
                         sng_732 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1017[k] = pb_x[k] * smh0_1017[k]
                    + f_10 * smg_729[k]
                    - f_8 * pc_x[k] * smh1_1017[k];

        t_1018[k] = f_9 * smg_730[k]
                    + f_3 * pc_x[k] * sng_730[k];

        t_1019[k] = f_9 * smg_731[k]
                    + f_3 * pc_x[k] * sng_731[k];

        t_1020[k] = f_9 * smg_732[k]
                    + f_3 * pc_x[k] * sng_732[k];
    }

#pragma omp simd aligned(t_1021, t_1022, t_1023, t_1024, pb_x, pc_x, pc_z, smh0_1023, smg_580, \
                         smg_733, smg_734, smh1_1023, sng_730, sng_733, \
                         sng_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1021[k] = f_9 * smg_733[k]
                    + f_3 * pc_x[k] * sng_733[k];

        t_1022[k] = f_9 * smg_734[k]
                    + f_3 * pc_x[k] * sng_734[k];

        t_1023[k] = pb_x[k] * smh0_1023[k]
                    - f_8 * pc_x[k] * smh1_1023[k];

        t_1024[k] = f_11 * smg_580[k]
                    + f_3 * pc_z[k] * sng_730[k];
    }

#pragma omp simd aligned(t_1025, t_1026, t_1027, t_1028, pb_x, pc_x, pc_y, smh0_1025, \
                         smh0_1026, smh0_1028, smg_599, smh1_1025, smh1_1026, smh1_1028, \
                         sng_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1025[k] = pb_x[k] * smh0_1025[k]
                    - f_8 * pc_x[k] * smh1_1025[k];

        t_1026[k] = pb_x[k] * smh0_1026[k]
                    - f_8 * pc_x[k] * smh1_1026[k];

        t_1027[k] = f_15 * smg_599[k]
                    + f_3 * pc_y[k] * sng_734[k];

        t_1028[k] = pb_x[k] * smh0_1028[k]
                    - f_8 * pc_x[k] * smh1_1028[k];
    }

#pragma omp simd aligned(t_1029, t_1030, t_1031, pb_x, pc_x, pc_y, pc_z, smh0_1029, smg_585, \
                         smg_600, smg_735, smh1_1029, sng_735 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1029[k] = pb_x[k] * smh0_1029[k]
                    + f_17 * smg_735[k]
                    - f_8 * pc_x[k] * smh1_1029[k];

        t_1030[k] = f_17 * smg_600[k]
                    + f_3 * pc_y[k] * sng_735[k];

        t_1031[k] = f_16 * smg_585[k]
                    + f_3 * pc_z[k] * sng_735[k];
    }

#pragma omp simd aligned(t_1032, t_1033, t_1034, pb_x, pc_x, pc_y, smh0_1032, smh0_1034, \
                         smg_602, smg_738, smg_740, smh1_1032, smh1_1034, \
                         sng_737 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1032[k] = pb_x[k] * smh0_1032[k]
                    + f_11 * smg_738[k]
                    - f_8 * pc_x[k] * smh1_1032[k];

        t_1033[k] = f_17 * smg_602[k]
                    + f_3 * pc_y[k] * sng_737[k];

        t_1034[k] = pb_x[k] * smh0_1034[k]
                    + f_11 * smg_740[k]
                    - f_8 * pc_x[k] * smh1_1034[k];
    }

#pragma omp simd aligned(t_1035, t_1036, t_1037, pb_x, pc_x, pc_y, pc_z, smh0_1035, smg_588, \
                         smg_605, smg_741, smh1_1035, sng_738, \
                         sng_740 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1035[k] = pb_x[k] * smh0_1035[k]
                    + f_10 * smg_741[k]
                    - f_8 * pc_x[k] * smh1_1035[k];

        t_1036[k] = f_16 * smg_588[k]
                    + f_3 * pc_z[k] * sng_738[k];

        t_1037[k] = f_17 * smg_605[k]
                    + f_3 * pc_y[k] * sng_740[k];
    }

#pragma omp simd aligned(t_1038, t_1039, t_1040, t_1041, pb_x, pc_x, smh0_1038, smg_744, \
                         smg_745, smg_746, smg_747, smh1_1038, sng_745, sng_746, \
                         sng_747 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1038[k] = pb_x[k] * smh0_1038[k]
                    + f_10 * smg_744[k]
                    - f_8 * pc_x[k] * smh1_1038[k];

        t_1039[k] = f_9 * smg_745[k]
                    + f_3 * pc_x[k] * sng_745[k];

        t_1040[k] = f_9 * smg_746[k]
                    + f_3 * pc_x[k] * sng_746[k];

        t_1041[k] = f_9 * smg_747[k]
                    + f_3 * pc_x[k] * sng_747[k];
    }
}

static auto
compute_prim_snh_three_center_electron_repulsion_0_piece9(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t smh0,
                                                          const size_t smg, const size_t smh1,
                                                          const size_t snf0, const size_t snf1,
                                                          const size_t sng, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.0 / q;
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
    const auto f_12 = 4.5 / q;
    const auto f_13 = 4.0 / q;
    const auto f_14 = 3.5 / q;
    const auto f_15 = 3.0 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;

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
    auto *t_1157 = buffer.data(target + 1157);
    auto *t_1158 = buffer.data(target + 1158);
    auto *t_1159 = buffer.data(target + 1159);
    auto *t_1160 = buffer.data(target + 1160);
    auto *t_1161 = buffer.data(target + 1161);
    auto *t_1162 = buffer.data(target + 1162);
    auto *t_1163 = buffer.data(target + 1163);
    auto *t_1164 = buffer.data(target + 1164);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *smh0_924 = buffer.data(smh0 + 924);
    const auto *smh0_929 = buffer.data(smh0 + 929);
    const auto *smh0_933 = buffer.data(smh0 + 933);
    const auto *smh0_1044 = buffer.data(smh0 + 1044);
    const auto *smh0_1046 = buffer.data(smh0 + 1046);
    const auto *smh0_1047 = buffer.data(smh0 + 1047);
    const auto *smh0_1049 = buffer.data(smh0 + 1049);
    const auto *smh0_1050 = buffer.data(smh0 + 1050);
    const auto *smh0_1053 = buffer.data(smh0 + 1053);
    const auto *smh0_1055 = buffer.data(smh0 + 1055);
    const auto *smh0_1056 = buffer.data(smh0 + 1056);
    const auto *smh0_1059 = buffer.data(smh0 + 1059);
    const auto *smh0_1065 = buffer.data(smh0 + 1065);
    const auto *smh0_1067 = buffer.data(smh0 + 1067);
    const auto *smh0_1068 = buffer.data(smh0 + 1068);
    const auto *smh0_1070 = buffer.data(smh0 + 1070);
    const auto *smh0_1071 = buffer.data(smh0 + 1071);
    const auto *smh0_1074 = buffer.data(smh0 + 1074);
    const auto *smh0_1076 = buffer.data(smh0 + 1076);
    const auto *smh0_1077 = buffer.data(smh0 + 1077);
    const auto *smh0_1080 = buffer.data(smh0 + 1080);
    const auto *smh0_1086 = buffer.data(smh0 + 1086);
    const auto *smh0_1088 = buffer.data(smh0 + 1088);
    const auto *smh0_1089 = buffer.data(smh0 + 1089);
    const auto *smh0_1091 = buffer.data(smh0 + 1091);
    const auto *smh0_1092 = buffer.data(smh0 + 1092);
    const auto *smh0_1095 = buffer.data(smh0 + 1095);
    const auto *smh0_1097 = buffer.data(smh0 + 1097);
    const auto *smh0_1098 = buffer.data(smh0 + 1098);
    const auto *smh0_1101 = buffer.data(smh0 + 1101);
    const auto *smh0_1107 = buffer.data(smh0 + 1107);
    const auto *smh0_1109 = buffer.data(smh0 + 1109);
    const auto *smh0_1110 = buffer.data(smh0 + 1110);
    const auto *smh0_1112 = buffer.data(smh0 + 1112);
    const auto *smh0_1116 = buffer.data(smh0 + 1116);
    const auto *smh0_1119 = buffer.data(smh0 + 1119);
    const auto *smh0_1128 = buffer.data(smh0 + 1128);
    const auto *smh0_1130 = buffer.data(smh0 + 1130);
    const auto *smh0_1131 = buffer.data(smh0 + 1131);
    const auto *smh0_1133 = buffer.data(smh0 + 1133);
    const auto *smh0_1134 = buffer.data(smh0 + 1134);
    const auto *smh0_1137 = buffer.data(smh0 + 1137);
    const auto *smh0_1139 = buffer.data(smh0 + 1139);
    const auto *smh0_1140 = buffer.data(smh0 + 1140);
    const auto *smh0_1143 = buffer.data(smh0 + 1143);
    const auto *smh0_1149 = buffer.data(smh0 + 1149);
    const auto *smh0_1151 = buffer.data(smh0 + 1151);
    const auto *smh0_1152 = buffer.data(smh0 + 1152);
    const auto *smh0_1154 = buffer.data(smh0 + 1154);

    const auto *smg_595 = buffer.data(smg + 595);
    const auto *smg_600 = buffer.data(smg + 600);
    const auto *smg_603 = buffer.data(smg + 603);
    const auto *smg_610 = buffer.data(smg + 610);
    const auto *smg_614 = buffer.data(smg + 614);
    const auto *smg_615 = buffer.data(smg + 615);
    const auto *smg_617 = buffer.data(smg + 617);
    const auto *smg_618 = buffer.data(smg + 618);
    const auto *smg_620 = buffer.data(smg + 620);
    const auto *smg_625 = buffer.data(smg + 625);
    const auto *smg_629 = buffer.data(smg + 629);
    const auto *smg_630 = buffer.data(smg + 630);
    const auto *smg_632 = buffer.data(smg + 632);
    const auto *smg_633 = buffer.data(smg + 633);
    const auto *smg_635 = buffer.data(smg + 635);
    const auto *smg_640 = buffer.data(smg + 640);
    const auto *smg_644 = buffer.data(smg + 644);
    const auto *smg_645 = buffer.data(smg + 645);
    const auto *smg_647 = buffer.data(smg + 647);
    const auto *smg_648 = buffer.data(smg + 648);
    const auto *smg_650 = buffer.data(smg + 650);
    const auto *smg_655 = buffer.data(smg + 655);
    const auto *smg_659 = buffer.data(smg + 659);
    const auto *smg_660 = buffer.data(smg + 660);
    const auto *smg_662 = buffer.data(smg + 662);
    const auto *smg_663 = buffer.data(smg + 663);
    const auto *smg_665 = buffer.data(smg + 665);
    const auto *smg_670 = buffer.data(smg + 670);
    const auto *smg_674 = buffer.data(smg + 674);
    const auto *smg_675 = buffer.data(smg + 675);
    const auto *smg_677 = buffer.data(smg + 677);
    const auto *smg_680 = buffer.data(smg + 680);
    const auto *smg_748 = buffer.data(smg + 748);
    const auto *smg_749 = buffer.data(smg + 749);
    const auto *smg_750 = buffer.data(smg + 750);
    const auto *smg_753 = buffer.data(smg + 753);
    const auto *smg_755 = buffer.data(smg + 755);
    const auto *smg_756 = buffer.data(smg + 756);
    const auto *smg_759 = buffer.data(smg + 759);
    const auto *smg_760 = buffer.data(smg + 760);
    const auto *smg_761 = buffer.data(smg + 761);
    const auto *smg_762 = buffer.data(smg + 762);
    const auto *smg_763 = buffer.data(smg + 763);
    const auto *smg_764 = buffer.data(smg + 764);
    const auto *smg_765 = buffer.data(smg + 765);
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
    const auto *smg_783 = buffer.data(smg + 783);
    const auto *smg_785 = buffer.data(smg + 785);
    const auto *smg_786 = buffer.data(smg + 786);
    const auto *smg_789 = buffer.data(smg + 789);
    const auto *smg_790 = buffer.data(smg + 790);
    const auto *smg_791 = buffer.data(smg + 791);
    const auto *smg_792 = buffer.data(smg + 792);
    const auto *smg_793 = buffer.data(smg + 793);
    const auto *smg_794 = buffer.data(smg + 794);
    const auto *smg_798 = buffer.data(smg + 798);
    const auto *smg_801 = buffer.data(smg + 801);
    const auto *smg_805 = buffer.data(smg + 805);
    const auto *smg_806 = buffer.data(smg + 806);
    const auto *smg_807 = buffer.data(smg + 807);
    const auto *smg_808 = buffer.data(smg + 808);
    const auto *smg_809 = buffer.data(smg + 809);
    const auto *smg_810 = buffer.data(smg + 810);
    const auto *smg_813 = buffer.data(smg + 813);
    const auto *smg_815 = buffer.data(smg + 815);
    const auto *smg_816 = buffer.data(smg + 816);
    const auto *smg_819 = buffer.data(smg + 819);
    const auto *smg_820 = buffer.data(smg + 820);
    const auto *smg_821 = buffer.data(smg + 821);
    const auto *smg_822 = buffer.data(smg + 822);
    const auto *smg_823 = buffer.data(smg + 823);
    const auto *smg_824 = buffer.data(smg + 824);

    const auto *smh1_924 = buffer.data(smh1 + 924);
    const auto *smh1_929 = buffer.data(smh1 + 929);
    const auto *smh1_933 = buffer.data(smh1 + 933);
    const auto *smh1_1044 = buffer.data(smh1 + 1044);
    const auto *smh1_1046 = buffer.data(smh1 + 1046);
    const auto *smh1_1047 = buffer.data(smh1 + 1047);
    const auto *smh1_1049 = buffer.data(smh1 + 1049);
    const auto *smh1_1050 = buffer.data(smh1 + 1050);
    const auto *smh1_1053 = buffer.data(smh1 + 1053);
    const auto *smh1_1055 = buffer.data(smh1 + 1055);
    const auto *smh1_1056 = buffer.data(smh1 + 1056);
    const auto *smh1_1059 = buffer.data(smh1 + 1059);
    const auto *smh1_1065 = buffer.data(smh1 + 1065);
    const auto *smh1_1067 = buffer.data(smh1 + 1067);
    const auto *smh1_1068 = buffer.data(smh1 + 1068);
    const auto *smh1_1070 = buffer.data(smh1 + 1070);
    const auto *smh1_1071 = buffer.data(smh1 + 1071);
    const auto *smh1_1074 = buffer.data(smh1 + 1074);
    const auto *smh1_1076 = buffer.data(smh1 + 1076);
    const auto *smh1_1077 = buffer.data(smh1 + 1077);
    const auto *smh1_1080 = buffer.data(smh1 + 1080);
    const auto *smh1_1086 = buffer.data(smh1 + 1086);
    const auto *smh1_1088 = buffer.data(smh1 + 1088);
    const auto *smh1_1089 = buffer.data(smh1 + 1089);
    const auto *smh1_1091 = buffer.data(smh1 + 1091);
    const auto *smh1_1092 = buffer.data(smh1 + 1092);
    const auto *smh1_1095 = buffer.data(smh1 + 1095);
    const auto *smh1_1097 = buffer.data(smh1 + 1097);
    const auto *smh1_1098 = buffer.data(smh1 + 1098);
    const auto *smh1_1101 = buffer.data(smh1 + 1101);
    const auto *smh1_1107 = buffer.data(smh1 + 1107);
    const auto *smh1_1109 = buffer.data(smh1 + 1109);
    const auto *smh1_1110 = buffer.data(smh1 + 1110);
    const auto *smh1_1112 = buffer.data(smh1 + 1112);
    const auto *smh1_1116 = buffer.data(smh1 + 1116);
    const auto *smh1_1119 = buffer.data(smh1 + 1119);
    const auto *smh1_1128 = buffer.data(smh1 + 1128);
    const auto *smh1_1130 = buffer.data(smh1 + 1130);
    const auto *smh1_1131 = buffer.data(smh1 + 1131);
    const auto *smh1_1133 = buffer.data(smh1 + 1133);
    const auto *smh1_1134 = buffer.data(smh1 + 1134);
    const auto *smh1_1137 = buffer.data(smh1 + 1137);
    const auto *smh1_1139 = buffer.data(smh1 + 1139);
    const auto *smh1_1140 = buffer.data(smh1 + 1140);
    const auto *smh1_1143 = buffer.data(smh1 + 1143);
    const auto *smh1_1149 = buffer.data(smh1 + 1149);
    const auto *smh1_1151 = buffer.data(smh1 + 1151);
    const auto *smh1_1152 = buffer.data(smh1 + 1152);
    const auto *smh1_1154 = buffer.data(smh1 + 1154);

    const auto *snf0_550 = buffer.data(snf0 + 550);
    const auto *snf0_553 = buffer.data(snf0 + 553);
    const auto *snf0_555 = buffer.data(snf0 + 555);
    const auto *snf0_556 = buffer.data(snf0 + 556);
    const auto *snf0_559 = buffer.data(snf0 + 559);

    const auto *snf1_550 = buffer.data(snf1 + 550);
    const auto *snf1_553 = buffer.data(snf1 + 553);
    const auto *snf1_555 = buffer.data(snf1 + 555);
    const auto *snf1_556 = buffer.data(snf1 + 556);
    const auto *snf1_559 = buffer.data(snf1 + 559);

    const auto *sng_745 = buffer.data(sng + 745);
    const auto *sng_748 = buffer.data(sng + 748);
    const auto *sng_749 = buffer.data(sng + 749);
    const auto *sng_750 = buffer.data(sng + 750);
    const auto *sng_752 = buffer.data(sng + 752);
    const auto *sng_753 = buffer.data(sng + 753);
    const auto *sng_755 = buffer.data(sng + 755);
    const auto *sng_760 = buffer.data(sng + 760);
    const auto *sng_761 = buffer.data(sng + 761);
    const auto *sng_762 = buffer.data(sng + 762);
    const auto *sng_763 = buffer.data(sng + 763);
    const auto *sng_764 = buffer.data(sng + 764);
    const auto *sng_765 = buffer.data(sng + 765);
    const auto *sng_767 = buffer.data(sng + 767);
    const auto *sng_768 = buffer.data(sng + 768);
    const auto *sng_770 = buffer.data(sng + 770);
    const auto *sng_775 = buffer.data(sng + 775);
    const auto *sng_776 = buffer.data(sng + 776);
    const auto *sng_777 = buffer.data(sng + 777);
    const auto *sng_778 = buffer.data(sng + 778);
    const auto *sng_779 = buffer.data(sng + 779);
    const auto *sng_780 = buffer.data(sng + 780);
    const auto *sng_782 = buffer.data(sng + 782);
    const auto *sng_783 = buffer.data(sng + 783);
    const auto *sng_785 = buffer.data(sng + 785);
    const auto *sng_790 = buffer.data(sng + 790);
    const auto *sng_791 = buffer.data(sng + 791);
    const auto *sng_792 = buffer.data(sng + 792);
    const auto *sng_793 = buffer.data(sng + 793);
    const auto *sng_794 = buffer.data(sng + 794);
    const auto *sng_795 = buffer.data(sng + 795);
    const auto *sng_797 = buffer.data(sng + 797);
    const auto *sng_798 = buffer.data(sng + 798);
    const auto *sng_800 = buffer.data(sng + 800);
    const auto *sng_805 = buffer.data(sng + 805);
    const auto *sng_806 = buffer.data(sng + 806);
    const auto *sng_807 = buffer.data(sng + 807);
    const auto *sng_808 = buffer.data(sng + 808);
    const auto *sng_809 = buffer.data(sng + 809);
    const auto *sng_810 = buffer.data(sng + 810);
    const auto *sng_812 = buffer.data(sng + 812);
    const auto *sng_813 = buffer.data(sng + 813);
    const auto *sng_815 = buffer.data(sng + 815);
    const auto *sng_820 = buffer.data(sng + 820);
    const auto *sng_821 = buffer.data(sng + 821);
    const auto *sng_822 = buffer.data(sng + 822);
    const auto *sng_823 = buffer.data(sng + 823);
    const auto *sng_824 = buffer.data(sng + 824);
    const auto *sng_825 = buffer.data(sng + 825);
    const auto *sng_827 = buffer.data(sng + 827);
    const auto *sng_828 = buffer.data(sng + 828);
    const auto *sng_830 = buffer.data(sng + 830);
    const auto *sng_831 = buffer.data(sng + 831);
    const auto *sng_834 = buffer.data(sng + 834);

#pragma omp simd aligned(t_1042, t_1043, t_1044, t_1045, pb_x, pc_x, pc_z, smh0_1044, smg_595, \
                         smg_748, smg_749, smh1_1044, sng_745, sng_748, \
                         sng_749 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1042[k] = f_9 * smg_748[k]
                    + f_3 * pc_x[k] * sng_748[k];

        t_1043[k] = f_9 * smg_749[k]
                    + f_3 * pc_x[k] * sng_749[k];

        t_1044[k] = pb_x[k] * smh0_1044[k]
                    - f_8 * pc_x[k] * smh1_1044[k];

        t_1045[k] = f_16 * smg_595[k]
                    + f_3 * pc_z[k] * sng_745[k];
    }

#pragma omp simd aligned(t_1046, t_1047, t_1048, t_1049, pb_x, pc_x, pc_y, smh0_1046, \
                         smh0_1047, smh0_1049, smg_614, smh1_1046, smh1_1047, smh1_1049, \
                         sng_749 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1046[k] = pb_x[k] * smh0_1046[k]
                    - f_8 * pc_x[k] * smh1_1046[k];

        t_1047[k] = pb_x[k] * smh0_1047[k]
                    - f_8 * pc_x[k] * smh1_1047[k];

        t_1048[k] = f_17 * smg_614[k]
                    + f_3 * pc_y[k] * sng_749[k];

        t_1049[k] = pb_x[k] * smh0_1049[k]
                    - f_8 * pc_x[k] * smh1_1049[k];
    }

#pragma omp simd aligned(t_1050, t_1051, t_1052, pb_x, pc_x, pc_y, pc_z, smh0_1050, smg_600, \
                         smg_615, smg_750, smh1_1050, sng_750 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1050[k] = pb_x[k] * smh0_1050[k]
                    + f_17 * smg_750[k]
                    - f_8 * pc_x[k] * smh1_1050[k];

        t_1051[k] = f_16 * smg_615[k]
                    + f_3 * pc_y[k] * sng_750[k];

        t_1052[k] = f_17 * smg_600[k]
                    + f_3 * pc_z[k] * sng_750[k];
    }

#pragma omp simd aligned(t_1053, t_1054, t_1055, pb_x, pc_x, pc_y, smh0_1053, smh0_1055, \
                         smg_617, smg_753, smg_755, smh1_1053, smh1_1055, \
                         sng_752 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1053[k] = pb_x[k] * smh0_1053[k]
                    + f_11 * smg_753[k]
                    - f_8 * pc_x[k] * smh1_1053[k];

        t_1054[k] = f_16 * smg_617[k]
                    + f_3 * pc_y[k] * sng_752[k];

        t_1055[k] = pb_x[k] * smh0_1055[k]
                    + f_11 * smg_755[k]
                    - f_8 * pc_x[k] * smh1_1055[k];
    }

#pragma omp simd aligned(t_1056, t_1057, t_1058, pb_x, pc_x, pc_y, pc_z, smh0_1056, smg_603, \
                         smg_620, smg_756, smh1_1056, sng_753, \
                         sng_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1056[k] = pb_x[k] * smh0_1056[k]
                    + f_10 * smg_756[k]
                    - f_8 * pc_x[k] * smh1_1056[k];

        t_1057[k] = f_17 * smg_603[k]
                    + f_3 * pc_z[k] * sng_753[k];

        t_1058[k] = f_16 * smg_620[k]
                    + f_3 * pc_y[k] * sng_755[k];
    }

#pragma omp simd aligned(t_1059, t_1060, t_1061, t_1062, pb_x, pc_x, smh0_1059, smg_759, \
                         smg_760, smg_761, smg_762, smh1_1059, sng_760, sng_761, \
                         sng_762 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1059[k] = pb_x[k] * smh0_1059[k]
                    + f_10 * smg_759[k]
                    - f_8 * pc_x[k] * smh1_1059[k];

        t_1060[k] = f_9 * smg_760[k]
                    + f_3 * pc_x[k] * sng_760[k];

        t_1061[k] = f_9 * smg_761[k]
                    + f_3 * pc_x[k] * sng_761[k];

        t_1062[k] = f_9 * smg_762[k]
                    + f_3 * pc_x[k] * sng_762[k];
    }

#pragma omp simd aligned(t_1063, t_1064, t_1065, t_1066, pb_x, pc_x, pc_z, smh0_1065, smg_610, \
                         smg_763, smg_764, smh1_1065, sng_760, sng_763, \
                         sng_764 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1063[k] = f_9 * smg_763[k]
                    + f_3 * pc_x[k] * sng_763[k];

        t_1064[k] = f_9 * smg_764[k]
                    + f_3 * pc_x[k] * sng_764[k];

        t_1065[k] = pb_x[k] * smh0_1065[k]
                    - f_8 * pc_x[k] * smh1_1065[k];

        t_1066[k] = f_17 * smg_610[k]
                    + f_3 * pc_z[k] * sng_760[k];
    }

#pragma omp simd aligned(t_1067, t_1068, t_1069, t_1070, pb_x, pc_x, pc_y, smh0_1067, \
                         smh0_1068, smh0_1070, smg_629, smh1_1067, smh1_1068, smh1_1070, \
                         sng_764 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1067[k] = pb_x[k] * smh0_1067[k]
                    - f_8 * pc_x[k] * smh1_1067[k];

        t_1068[k] = pb_x[k] * smh0_1068[k]
                    - f_8 * pc_x[k] * smh1_1068[k];

        t_1069[k] = f_16 * smg_629[k]
                    + f_3 * pc_y[k] * sng_764[k];

        t_1070[k] = pb_x[k] * smh0_1070[k]
                    - f_8 * pc_x[k] * smh1_1070[k];
    }

#pragma omp simd aligned(t_1071, t_1072, t_1073, pb_x, pc_x, pc_y, pc_z, smh0_1071, smg_615, \
                         smg_630, smg_765, smh1_1071, sng_765 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1071[k] = pb_x[k] * smh0_1071[k]
                    + f_17 * smg_765[k]
                    - f_8 * pc_x[k] * smh1_1071[k];

        t_1072[k] = f_11 * smg_630[k]
                    + f_3 * pc_y[k] * sng_765[k];

        t_1073[k] = f_15 * smg_615[k]
                    + f_3 * pc_z[k] * sng_765[k];
    }

#pragma omp simd aligned(t_1074, t_1075, t_1076, pb_x, pc_x, pc_y, smh0_1074, smh0_1076, \
                         smg_632, smg_768, smg_770, smh1_1074, smh1_1076, \
                         sng_767 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1074[k] = pb_x[k] * smh0_1074[k]
                    + f_11 * smg_768[k]
                    - f_8 * pc_x[k] * smh1_1074[k];

        t_1075[k] = f_11 * smg_632[k]
                    + f_3 * pc_y[k] * sng_767[k];

        t_1076[k] = pb_x[k] * smh0_1076[k]
                    + f_11 * smg_770[k]
                    - f_8 * pc_x[k] * smh1_1076[k];
    }

#pragma omp simd aligned(t_1077, t_1078, t_1079, pb_x, pc_x, pc_y, pc_z, smh0_1077, smg_618, \
                         smg_635, smg_771, smh1_1077, sng_768, \
                         sng_770 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1077[k] = pb_x[k] * smh0_1077[k]
                    + f_10 * smg_771[k]
                    - f_8 * pc_x[k] * smh1_1077[k];

        t_1078[k] = f_15 * smg_618[k]
                    + f_3 * pc_z[k] * sng_768[k];

        t_1079[k] = f_11 * smg_635[k]
                    + f_3 * pc_y[k] * sng_770[k];
    }

#pragma omp simd aligned(t_1080, t_1081, t_1082, t_1083, pb_x, pc_x, smh0_1080, smg_774, \
                         smg_775, smg_776, smg_777, smh1_1080, sng_775, sng_776, \
                         sng_777 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1080[k] = pb_x[k] * smh0_1080[k]
                    + f_10 * smg_774[k]
                    - f_8 * pc_x[k] * smh1_1080[k];

        t_1081[k] = f_9 * smg_775[k]
                    + f_3 * pc_x[k] * sng_775[k];

        t_1082[k] = f_9 * smg_776[k]
                    + f_3 * pc_x[k] * sng_776[k];

        t_1083[k] = f_9 * smg_777[k]
                    + f_3 * pc_x[k] * sng_777[k];
    }

#pragma omp simd aligned(t_1084, t_1085, t_1086, t_1087, pb_x, pc_x, pc_z, smh0_1086, smg_625, \
                         smg_778, smg_779, smh1_1086, sng_775, sng_778, \
                         sng_779 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1084[k] = f_9 * smg_778[k]
                    + f_3 * pc_x[k] * sng_778[k];

        t_1085[k] = f_9 * smg_779[k]
                    + f_3 * pc_x[k] * sng_779[k];

        t_1086[k] = pb_x[k] * smh0_1086[k]
                    - f_8 * pc_x[k] * smh1_1086[k];

        t_1087[k] = f_15 * smg_625[k]
                    + f_3 * pc_z[k] * sng_775[k];
    }

#pragma omp simd aligned(t_1088, t_1089, t_1090, t_1091, pb_x, pc_x, pc_y, smh0_1088, \
                         smh0_1089, smh0_1091, smg_644, smh1_1088, smh1_1089, smh1_1091, \
                         sng_779 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1088[k] = pb_x[k] * smh0_1088[k]
                    - f_8 * pc_x[k] * smh1_1088[k];

        t_1089[k] = pb_x[k] * smh0_1089[k]
                    - f_8 * pc_x[k] * smh1_1089[k];

        t_1090[k] = f_11 * smg_644[k]
                    + f_3 * pc_y[k] * sng_779[k];

        t_1091[k] = pb_x[k] * smh0_1091[k]
                    - f_8 * pc_x[k] * smh1_1091[k];
    }

#pragma omp simd aligned(t_1092, t_1093, t_1094, pb_x, pc_x, pc_y, pc_z, smh0_1092, smg_630, \
                         smg_645, smg_780, smh1_1092, sng_780 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1092[k] = pb_x[k] * smh0_1092[k]
                    + f_17 * smg_780[k]
                    - f_8 * pc_x[k] * smh1_1092[k];

        t_1093[k] = f_10 * smg_645[k]
                    + f_3 * pc_y[k] * sng_780[k];

        t_1094[k] = f_14 * smg_630[k]
                    + f_3 * pc_z[k] * sng_780[k];
    }

#pragma omp simd aligned(t_1095, t_1096, t_1097, pb_x, pc_x, pc_y, smh0_1095, smh0_1097, \
                         smg_647, smg_783, smg_785, smh1_1095, smh1_1097, \
                         sng_782 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1095[k] = pb_x[k] * smh0_1095[k]
                    + f_11 * smg_783[k]
                    - f_8 * pc_x[k] * smh1_1095[k];

        t_1096[k] = f_10 * smg_647[k]
                    + f_3 * pc_y[k] * sng_782[k];

        t_1097[k] = pb_x[k] * smh0_1097[k]
                    + f_11 * smg_785[k]
                    - f_8 * pc_x[k] * smh1_1097[k];
    }

#pragma omp simd aligned(t_1098, t_1099, t_1100, pb_x, pc_x, pc_y, pc_z, smh0_1098, smg_633, \
                         smg_650, smg_786, smh1_1098, sng_783, \
                         sng_785 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1098[k] = pb_x[k] * smh0_1098[k]
                    + f_10 * smg_786[k]
                    - f_8 * pc_x[k] * smh1_1098[k];

        t_1099[k] = f_14 * smg_633[k]
                    + f_3 * pc_z[k] * sng_783[k];

        t_1100[k] = f_10 * smg_650[k]
                    + f_3 * pc_y[k] * sng_785[k];
    }

#pragma omp simd aligned(t_1101, t_1102, t_1103, t_1104, pb_x, pc_x, smh0_1101, smg_789, \
                         smg_790, smg_791, smg_792, smh1_1101, sng_790, sng_791, \
                         sng_792 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1101[k] = pb_x[k] * smh0_1101[k]
                    + f_10 * smg_789[k]
                    - f_8 * pc_x[k] * smh1_1101[k];

        t_1102[k] = f_9 * smg_790[k]
                    + f_3 * pc_x[k] * sng_790[k];

        t_1103[k] = f_9 * smg_791[k]
                    + f_3 * pc_x[k] * sng_791[k];

        t_1104[k] = f_9 * smg_792[k]
                    + f_3 * pc_x[k] * sng_792[k];
    }

#pragma omp simd aligned(t_1105, t_1106, t_1107, t_1108, pb_x, pc_x, pc_z, smh0_1107, smg_640, \
                         smg_793, smg_794, smh1_1107, sng_790, sng_793, \
                         sng_794 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1105[k] = f_9 * smg_793[k]
                    + f_3 * pc_x[k] * sng_793[k];

        t_1106[k] = f_9 * smg_794[k]
                    + f_3 * pc_x[k] * sng_794[k];

        t_1107[k] = pb_x[k] * smh0_1107[k]
                    - f_8 * pc_x[k] * smh1_1107[k];

        t_1108[k] = f_14 * smg_640[k]
                    + f_3 * pc_z[k] * sng_790[k];
    }

#pragma omp simd aligned(t_1109, t_1110, t_1111, t_1112, pb_x, pc_x, pc_y, smh0_1109, \
                         smh0_1110, smh0_1112, smg_659, smh1_1109, smh1_1110, smh1_1112, \
                         sng_794 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1109[k] = pb_x[k] * smh0_1109[k]
                    - f_8 * pc_x[k] * smh1_1109[k];

        t_1110[k] = pb_x[k] * smh0_1110[k]
                    - f_8 * pc_x[k] * smh1_1110[k];

        t_1111[k] = f_10 * smg_659[k]
                    + f_3 * pc_y[k] * sng_794[k];

        t_1112[k] = pb_x[k] * smh0_1112[k]
                    - f_8 * pc_x[k] * smh1_1112[k];
    }

#pragma omp simd aligned(t_1113, t_1114, t_1115, pb_y, pc_y, pc_z, smh0_924, smg_645, smg_660, \
                         smh1_924, sng_795 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1113[k] = pb_y[k] * smh0_924[k]
                    - f_8 * pc_y[k] * smh1_924[k];

        t_1114[k] = f_9 * smg_660[k]
                    + f_3 * pc_y[k] * sng_795[k];

        t_1115[k] = f_13 * smg_645[k]
                    + f_3 * pc_z[k] * sng_795[k];
    }

#pragma omp simd aligned(t_1116, t_1117, t_1118, pb_x, pb_y, pc_x, pc_y, smh0_929, smh0_1116, \
                         smg_662, smg_798, smh1_929, smh1_1116, \
                         sng_797 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1116[k] = pb_x[k] * smh0_1116[k]
                    + f_11 * smg_798[k]
                    - f_8 * pc_x[k] * smh1_1116[k];

        t_1117[k] = f_9 * smg_662[k]
                    + f_3 * pc_y[k] * sng_797[k];

        t_1118[k] = pb_y[k] * smh0_929[k]
                    - f_8 * pc_y[k] * smh1_929[k];
    }

#pragma omp simd aligned(t_1119, t_1120, t_1121, pb_x, pc_x, pc_y, pc_z, smh0_1119, smg_648, \
                         smg_665, smg_801, smh1_1119, sng_798, \
                         sng_800 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1119[k] = pb_x[k] * smh0_1119[k]
                    + f_10 * smg_801[k]
                    - f_8 * pc_x[k] * smh1_1119[k];

        t_1120[k] = f_13 * smg_648[k]
                    + f_3 * pc_z[k] * sng_798[k];

        t_1121[k] = f_9 * smg_665[k]
                    + f_3 * pc_y[k] * sng_800[k];
    }

#pragma omp simd aligned(t_1122, t_1123, t_1124, t_1125, pb_y, pc_x, pc_y, smh0_933, smg_805, \
                         smg_806, smg_807, smh1_933, sng_805, sng_806, \
                         sng_807 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1122[k] = pb_y[k] * smh0_933[k]
                    - f_8 * pc_y[k] * smh1_933[k];

        t_1123[k] = f_9 * smg_805[k]
                    + f_3 * pc_x[k] * sng_805[k];

        t_1124[k] = f_9 * smg_806[k]
                    + f_3 * pc_x[k] * sng_806[k];

        t_1125[k] = f_9 * smg_807[k]
                    + f_3 * pc_x[k] * sng_807[k];
    }

#pragma omp simd aligned(t_1126, t_1127, t_1128, t_1129, pb_x, pc_x, pc_z, smh0_1128, smg_655, \
                         smg_808, smg_809, smh1_1128, sng_805, sng_808, \
                         sng_809 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1126[k] = f_9 * smg_808[k]
                    + f_3 * pc_x[k] * sng_808[k];

        t_1127[k] = f_9 * smg_809[k]
                    + f_3 * pc_x[k] * sng_809[k];

        t_1128[k] = pb_x[k] * smh0_1128[k]
                    - f_8 * pc_x[k] * smh1_1128[k];

        t_1129[k] = f_13 * smg_655[k]
                    + f_3 * pc_z[k] * sng_805[k];
    }

#pragma omp simd aligned(t_1130, t_1131, t_1132, t_1133, pb_x, pc_x, pc_y, smh0_1130, \
                         smh0_1131, smh0_1133, smg_674, smh1_1130, smh1_1131, smh1_1133, \
                         sng_809 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1130[k] = pb_x[k] * smh0_1130[k]
                    - f_8 * pc_x[k] * smh1_1130[k];

        t_1131[k] = pb_x[k] * smh0_1131[k]
                    - f_8 * pc_x[k] * smh1_1131[k];

        t_1132[k] = f_9 * smg_674[k]
                    + f_3 * pc_y[k] * sng_809[k];

        t_1133[k] = pb_x[k] * smh0_1133[k]
                    - f_8 * pc_x[k] * smh1_1133[k];
    }

#pragma omp simd aligned(t_1134, t_1135, t_1136, t_1137, pb_x, pc_x, pc_y, pc_z, smh0_1134, \
                         smh0_1137, smg_660, smg_810, smg_813, smh1_1134, smh1_1137, \
                         sng_810 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1134[k] = pb_x[k] * smh0_1134[k]
                    + f_17 * smg_810[k]
                    - f_8 * pc_x[k] * smh1_1134[k];

        t_1135[k] = f_3 * pc_y[k] * sng_810[k];

        t_1136[k] = f_12 * smg_660[k]
                    + f_3 * pc_z[k] * sng_810[k];

        t_1137[k] = pb_x[k] * smh0_1137[k]
                    + f_11 * smg_813[k]
                    - f_8 * pc_x[k] * smh1_1137[k];
    }

#pragma omp simd aligned(t_1138, t_1139, t_1140, pb_x, pc_x, pc_y, smh0_1139, smh0_1140, \
                         smg_815, smg_816, smh1_1139, smh1_1140, \
                         sng_812 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1138[k] = f_3 * pc_y[k] * sng_812[k];

        t_1139[k] = pb_x[k] * smh0_1139[k]
                    + f_11 * smg_815[k]
                    - f_8 * pc_x[k] * smh1_1139[k];

        t_1140[k] = pb_x[k] * smh0_1140[k]
                    + f_10 * smg_816[k]
                    - f_8 * pc_x[k] * smh1_1140[k];
    }

#pragma omp simd aligned(t_1141, t_1142, t_1143, t_1144, pb_x, pc_x, pc_y, pc_z, smh0_1143, \
                         smg_663, smg_819, smg_820, smh1_1143, sng_813, sng_815, \
                         sng_820 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1141[k] = f_12 * smg_663[k]
                    + f_3 * pc_z[k] * sng_813[k];

        t_1142[k] = f_3 * pc_y[k] * sng_815[k];

        t_1143[k] = pb_x[k] * smh0_1143[k]
                    + f_10 * smg_819[k]
                    - f_8 * pc_x[k] * smh1_1143[k];

        t_1144[k] = f_9 * smg_820[k]
                    + f_3 * pc_x[k] * sng_820[k];
    }

#pragma omp simd aligned(t_1145, t_1146, t_1147, t_1148, pc_x, smg_821, smg_822, smg_823, \
                         smg_824, sng_821, sng_822, sng_823, sng_824 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1145[k] = f_9 * smg_821[k]
                    + f_3 * pc_x[k] * sng_821[k];

        t_1146[k] = f_9 * smg_822[k]
                    + f_3 * pc_x[k] * sng_822[k];

        t_1147[k] = f_9 * smg_823[k]
                    + f_3 * pc_x[k] * sng_823[k];

        t_1148[k] = f_9 * smg_824[k]
                    + f_3 * pc_x[k] * sng_824[k];
    }

#pragma omp simd aligned(t_1149, t_1150, t_1151, t_1152, pb_x, pc_x, pc_z, smh0_1149, \
                         smh0_1151, smh0_1152, smg_670, smh1_1149, smh1_1151, smh1_1152, \
                         sng_820 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1149[k] = pb_x[k] * smh0_1149[k]
                    - f_8 * pc_x[k] * smh1_1149[k];

        t_1150[k] = f_12 * smg_670[k]
                    + f_3 * pc_z[k] * sng_820[k];

        t_1151[k] = pb_x[k] * smh0_1151[k]
                    - f_8 * pc_x[k] * smh1_1151[k];

        t_1152[k] = pb_x[k] * smh0_1152[k]
                    - f_8 * pc_x[k] * smh1_1152[k];
    }

#pragma omp simd aligned(t_1153, t_1154, t_1155, t_1156, t_1157, pb_x, pc_x, pc_y, pc_z, \
                         smh0_1154, smg_675, smh1_1154, snf0_550, snf1_550, sng_824, \
                         sng_825 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1153[k] = f_3 * pc_y[k] * sng_824[k];

        t_1154[k] = pb_x[k] * smh0_1154[k]
                    - f_8 * pc_x[k] * smh1_1154[k];

        t_1155[k] = f_1 * snf0_550[k]
                    - f_2 * snf1_550[k]
                    + f_3 * pc_x[k] * sng_825[k];

        t_1156[k] = f_0 * smg_675[k]
                    + f_3 * pc_y[k] * sng_825[k];

        t_1157[k] = f_3 * pc_z[k] * sng_825[k];
    }

#pragma omp simd aligned(t_1158, t_1159, t_1160, pc_x, pc_y, smg_677, snf0_553, snf0_555, \
                         snf1_553, snf1_555, sng_827, sng_828, \
                         sng_830 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1158[k] = f_4 * snf0_553[k]
                    - f_5 * snf1_553[k]
                    + f_3 * pc_x[k] * sng_828[k];

        t_1159[k] = f_0 * smg_677[k]
                    + f_3 * pc_y[k] * sng_827[k];

        t_1160[k] = f_4 * snf0_555[k]
                    - f_5 * snf1_555[k]
                    + f_3 * pc_x[k] * sng_830[k];
    }

#pragma omp simd aligned(t_1161, t_1162, t_1163, t_1164, pc_x, pc_y, pc_z, smg_680, snf0_556, \
                         snf0_559, snf1_556, snf1_559, sng_828, sng_830, sng_831, \
                         sng_834 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1161[k] = f_6 * snf0_556[k]
                    - f_7 * snf1_556[k]
                    + f_3 * pc_x[k] * sng_831[k];

        t_1162[k] = f_3 * pc_z[k] * sng_828[k];

        t_1163[k] = f_0 * smg_680[k]
                    + f_3 * pc_y[k] * sng_830[k];

        t_1164[k] = f_6 * snf0_559[k]
                    - f_7 * snf1_559[k]
                    + f_3 * pc_x[k] * sng_834[k];
    }
}

static auto
compute_prim_snh_three_center_electron_repulsion_0_piece10(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t smh0,
                                                           const size_t smg, const size_t smh1,
                                                           const size_t snf0, const size_t snf1,
                                                           const size_t sng, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.0 / q;
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
    const auto f_12 = 4.5 / q;
    const auto f_13 = 4.0 / q;
    const auto f_14 = 3.5 / q;
    const auto f_15 = 3.0 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;

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

    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *smh0_945 = buffer.data(smh0 + 945);
    const auto *smh0_948 = buffer.data(smh0 + 948);
    const auto *smh0_951 = buffer.data(smh0 + 951);
    const auto *smh0_960 = buffer.data(smh0 + 960);
    const auto *smh0_962 = buffer.data(smh0 + 962);
    const auto *smh0_963 = buffer.data(smh0 + 963);

    const auto *smg_675 = buffer.data(smg + 675);
    const auto *smg_678 = buffer.data(smg + 678);
    const auto *smg_685 = buffer.data(smg + 685);
    const auto *smg_686 = buffer.data(smg + 686);
    const auto *smg_687 = buffer.data(smg + 687);
    const auto *smg_688 = buffer.data(smg + 688);
    const auto *smg_689 = buffer.data(smg + 689);
    const auto *smg_690 = buffer.data(smg + 690);
    const auto *smg_692 = buffer.data(smg + 692);
    const auto *smg_693 = buffer.data(smg + 693);
    const auto *smg_695 = buffer.data(smg + 695);
    const auto *smg_700 = buffer.data(smg + 700);
    const auto *smg_704 = buffer.data(smg + 704);
    const auto *smg_705 = buffer.data(smg + 705);
    const auto *smg_707 = buffer.data(smg + 707);
    const auto *smg_708 = buffer.data(smg + 708);
    const auto *smg_710 = buffer.data(smg + 710);
    const auto *smg_715 = buffer.data(smg + 715);
    const auto *smg_717 = buffer.data(smg + 717);
    const auto *smg_718 = buffer.data(smg + 718);
    const auto *smg_719 = buffer.data(smg + 719);
    const auto *smg_720 = buffer.data(smg + 720);
    const auto *smg_722 = buffer.data(smg + 722);
    const auto *smg_723 = buffer.data(smg + 723);
    const auto *smg_725 = buffer.data(smg + 725);
    const auto *smg_730 = buffer.data(smg + 730);
    const auto *smg_732 = buffer.data(smg + 732);
    const auto *smg_733 = buffer.data(smg + 733);
    const auto *smg_734 = buffer.data(smg + 734);
    const auto *smg_735 = buffer.data(smg + 735);
    const auto *smg_737 = buffer.data(smg + 737);
    const auto *smg_738 = buffer.data(smg + 738);
    const auto *smg_740 = buffer.data(smg + 740);
    const auto *smg_745 = buffer.data(smg + 745);
    const auto *smg_747 = buffer.data(smg + 747);
    const auto *smg_748 = buffer.data(smg + 748);
    const auto *smg_749 = buffer.data(smg + 749);
    const auto *smg_750 = buffer.data(smg + 750);
    const auto *smg_752 = buffer.data(smg + 752);
    const auto *smg_753 = buffer.data(smg + 753);
    const auto *smg_755 = buffer.data(smg + 755);
    const auto *smg_760 = buffer.data(smg + 760);
    const auto *smg_762 = buffer.data(smg + 762);
    const auto *smg_763 = buffer.data(smg + 763);
    const auto *smg_764 = buffer.data(smg + 764);
    const auto *smg_765 = buffer.data(smg + 765);
    const auto *smg_767 = buffer.data(smg + 767);
    const auto *smg_770 = buffer.data(smg + 770);

    const auto *smh1_945 = buffer.data(smh1 + 945);
    const auto *smh1_948 = buffer.data(smh1 + 948);
    const auto *smh1_951 = buffer.data(smh1 + 951);
    const auto *smh1_960 = buffer.data(smh1 + 960);
    const auto *smh1_962 = buffer.data(smh1 + 962);
    const auto *smh1_963 = buffer.data(smh1 + 963);

    const auto *snf0_556 = buffer.data(snf0 + 556);
    const auto *snf0_558 = buffer.data(snf0 + 558);
    const auto *snf0_559 = buffer.data(snf0 + 559);
    const auto *snf0_565 = buffer.data(snf0 + 565);
    const auto *snf0_569 = buffer.data(snf0 + 569);
    const auto *snf0_570 = buffer.data(snf0 + 570);
    const auto *snf0_573 = buffer.data(snf0 + 573);
    const auto *snf0_575 = buffer.data(snf0 + 575);
    const auto *snf0_576 = buffer.data(snf0 + 576);
    const auto *snf0_578 = buffer.data(snf0 + 578);
    const auto *snf0_579 = buffer.data(snf0 + 579);
    const auto *snf0_580 = buffer.data(snf0 + 580);
    const auto *snf0_583 = buffer.data(snf0 + 583);
    const auto *snf0_585 = buffer.data(snf0 + 585);
    const auto *snf0_586 = buffer.data(snf0 + 586);
    const auto *snf0_588 = buffer.data(snf0 + 588);
    const auto *snf0_589 = buffer.data(snf0 + 589);
    const auto *snf0_590 = buffer.data(snf0 + 590);
    const auto *snf0_593 = buffer.data(snf0 + 593);
    const auto *snf0_595 = buffer.data(snf0 + 595);
    const auto *snf0_596 = buffer.data(snf0 + 596);
    const auto *snf0_598 = buffer.data(snf0 + 598);
    const auto *snf0_599 = buffer.data(snf0 + 599);
    const auto *snf0_600 = buffer.data(snf0 + 600);
    const auto *snf0_603 = buffer.data(snf0 + 603);
    const auto *snf0_605 = buffer.data(snf0 + 605);
    const auto *snf0_606 = buffer.data(snf0 + 606);
    const auto *snf0_608 = buffer.data(snf0 + 608);
    const auto *snf0_609 = buffer.data(snf0 + 609);
    const auto *snf0_610 = buffer.data(snf0 + 610);
    const auto *snf0_613 = buffer.data(snf0 + 613);
    const auto *snf0_615 = buffer.data(snf0 + 615);
    const auto *snf0_616 = buffer.data(snf0 + 616);

    const auto *snf1_556 = buffer.data(snf1 + 556);
    const auto *snf1_558 = buffer.data(snf1 + 558);
    const auto *snf1_559 = buffer.data(snf1 + 559);
    const auto *snf1_565 = buffer.data(snf1 + 565);
    const auto *snf1_569 = buffer.data(snf1 + 569);
    const auto *snf1_570 = buffer.data(snf1 + 570);
    const auto *snf1_573 = buffer.data(snf1 + 573);
    const auto *snf1_575 = buffer.data(snf1 + 575);
    const auto *snf1_576 = buffer.data(snf1 + 576);
    const auto *snf1_578 = buffer.data(snf1 + 578);
    const auto *snf1_579 = buffer.data(snf1 + 579);
    const auto *snf1_580 = buffer.data(snf1 + 580);
    const auto *snf1_583 = buffer.data(snf1 + 583);
    const auto *snf1_585 = buffer.data(snf1 + 585);
    const auto *snf1_586 = buffer.data(snf1 + 586);
    const auto *snf1_588 = buffer.data(snf1 + 588);
    const auto *snf1_589 = buffer.data(snf1 + 589);
    const auto *snf1_590 = buffer.data(snf1 + 590);
    const auto *snf1_593 = buffer.data(snf1 + 593);
    const auto *snf1_595 = buffer.data(snf1 + 595);
    const auto *snf1_596 = buffer.data(snf1 + 596);
    const auto *snf1_598 = buffer.data(snf1 + 598);
    const auto *snf1_599 = buffer.data(snf1 + 599);
    const auto *snf1_600 = buffer.data(snf1 + 600);
    const auto *snf1_603 = buffer.data(snf1 + 603);
    const auto *snf1_605 = buffer.data(snf1 + 605);
    const auto *snf1_606 = buffer.data(snf1 + 606);
    const auto *snf1_608 = buffer.data(snf1 + 608);
    const auto *snf1_609 = buffer.data(snf1 + 609);
    const auto *snf1_610 = buffer.data(snf1 + 610);
    const auto *snf1_613 = buffer.data(snf1 + 613);
    const auto *snf1_615 = buffer.data(snf1 + 615);
    const auto *snf1_616 = buffer.data(snf1 + 616);

    const auto *sng_835 = buffer.data(sng + 835);
    const auto *sng_836 = buffer.data(sng + 836);
    const auto *sng_837 = buffer.data(sng + 837);
    const auto *sng_838 = buffer.data(sng + 838);
    const auto *sng_839 = buffer.data(sng + 839);
    const auto *sng_840 = buffer.data(sng + 840);
    const auto *sng_842 = buffer.data(sng + 842);
    const auto *sng_843 = buffer.data(sng + 843);
    const auto *sng_845 = buffer.data(sng + 845);
    const auto *sng_849 = buffer.data(sng + 849);
    const auto *sng_850 = buffer.data(sng + 850);
    const auto *sng_851 = buffer.data(sng + 851);
    const auto *sng_852 = buffer.data(sng + 852);
    const auto *sng_853 = buffer.data(sng + 853);
    const auto *sng_854 = buffer.data(sng + 854);
    const auto *sng_855 = buffer.data(sng + 855);
    const auto *sng_857 = buffer.data(sng + 857);
    const auto *sng_858 = buffer.data(sng + 858);
    const auto *sng_860 = buffer.data(sng + 860);
    const auto *sng_861 = buffer.data(sng + 861);
    const auto *sng_864 = buffer.data(sng + 864);
    const auto *sng_865 = buffer.data(sng + 865);
    const auto *sng_866 = buffer.data(sng + 866);
    const auto *sng_867 = buffer.data(sng + 867);
    const auto *sng_868 = buffer.data(sng + 868);
    const auto *sng_869 = buffer.data(sng + 869);
    const auto *sng_870 = buffer.data(sng + 870);
    const auto *sng_872 = buffer.data(sng + 872);
    const auto *sng_873 = buffer.data(sng + 873);
    const auto *sng_875 = buffer.data(sng + 875);
    const auto *sng_876 = buffer.data(sng + 876);
    const auto *sng_879 = buffer.data(sng + 879);
    const auto *sng_880 = buffer.data(sng + 880);
    const auto *sng_881 = buffer.data(sng + 881);
    const auto *sng_882 = buffer.data(sng + 882);
    const auto *sng_883 = buffer.data(sng + 883);
    const auto *sng_884 = buffer.data(sng + 884);
    const auto *sng_885 = buffer.data(sng + 885);
    const auto *sng_887 = buffer.data(sng + 887);
    const auto *sng_888 = buffer.data(sng + 888);
    const auto *sng_890 = buffer.data(sng + 890);
    const auto *sng_891 = buffer.data(sng + 891);
    const auto *sng_894 = buffer.data(sng + 894);
    const auto *sng_895 = buffer.data(sng + 895);
    const auto *sng_896 = buffer.data(sng + 896);
    const auto *sng_897 = buffer.data(sng + 897);
    const auto *sng_898 = buffer.data(sng + 898);
    const auto *sng_899 = buffer.data(sng + 899);
    const auto *sng_900 = buffer.data(sng + 900);
    const auto *sng_902 = buffer.data(sng + 902);
    const auto *sng_903 = buffer.data(sng + 903);
    const auto *sng_905 = buffer.data(sng + 905);
    const auto *sng_906 = buffer.data(sng + 906);
    const auto *sng_909 = buffer.data(sng + 909);
    const auto *sng_910 = buffer.data(sng + 910);
    const auto *sng_911 = buffer.data(sng + 911);
    const auto *sng_912 = buffer.data(sng + 912);
    const auto *sng_913 = buffer.data(sng + 913);
    const auto *sng_914 = buffer.data(sng + 914);
    const auto *sng_915 = buffer.data(sng + 915);
    const auto *sng_917 = buffer.data(sng + 917);
    const auto *sng_918 = buffer.data(sng + 918);
    const auto *sng_920 = buffer.data(sng + 920);
    const auto *sng_921 = buffer.data(sng + 921);

#pragma omp simd aligned(t_1165, t_1166, t_1167, t_1168, t_1169, t_1170, pc_x, pc_y, smg_685, \
                         snf0_556, snf1_556, sng_835, sng_836, sng_837, sng_838, \
                         sng_839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1165[k] = f_3 * pc_x[k] * sng_835[k];

        t_1166[k] = f_3 * pc_x[k] * sng_836[k];

        t_1167[k] = f_3 * pc_x[k] * sng_837[k];

        t_1168[k] = f_3 * pc_x[k] * sng_838[k];

        t_1169[k] = f_3 * pc_x[k] * sng_839[k];

        t_1170[k] = f_0 * smg_685[k]
                    + f_1 * snf0_556[k]
                    - f_2 * snf1_556[k]
                    + f_3 * pc_y[k] * sng_835[k];
    }

#pragma omp simd aligned(t_1171, t_1172, t_1173, pc_y, pc_z, smg_687, smg_688, snf0_558, \
                         snf0_559, snf1_558, snf1_559, sng_835, sng_837, \
                         sng_838 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1171[k] = f_3 * pc_z[k] * sng_835[k];

        t_1172[k] = f_0 * smg_687[k]
                    + f_4 * snf0_558[k]
                    - f_5 * snf1_558[k]
                    + f_3 * pc_y[k] * sng_837[k];

        t_1173[k] = f_0 * smg_688[k]
                    + f_6 * snf0_559[k]
                    - f_7 * snf1_559[k]
                    + f_3 * pc_y[k] * sng_838[k];
    }

#pragma omp simd aligned(t_1174, t_1175, t_1176, t_1177, pb_z, pc_y, pc_z, smh0_945, smg_689, \
                         smg_690, smh1_945, snf0_559, snf1_559, sng_839, \
                         sng_840 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1174[k] = f_0 * smg_689[k]
                    + f_3 * pc_y[k] * sng_839[k];

        t_1175[k] = f_1 * snf0_559[k]
                    - f_2 * snf1_559[k]
                    + f_3 * pc_z[k] * sng_839[k];

        t_1176[k] = pb_z[k] * smh0_945[k]
                    - f_8 * pc_z[k] * smh1_945[k];

        t_1177[k] = f_12 * smg_690[k]
                    + f_3 * pc_y[k] * sng_840[k];
    }

#pragma omp simd aligned(t_1178, t_1179, t_1180, pb_z, pc_y, pc_z, smh0_948, smg_675, smg_692, \
                         smh1_948, sng_840, sng_842 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1178[k] = f_9 * smg_675[k]
                    + f_3 * pc_z[k] * sng_840[k];

        t_1179[k] = pb_z[k] * smh0_948[k]
                    - f_8 * pc_z[k] * smh1_948[k];

        t_1180[k] = f_12 * smg_692[k]
                    + f_3 * pc_y[k] * sng_842[k];
    }

#pragma omp simd aligned(t_1181, t_1182, t_1183, t_1184, pb_z, pc_x, pc_y, pc_z, smh0_951, \
                         smg_678, smg_695, smh1_951, snf0_565, snf1_565, sng_843, \
                         sng_845 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1181[k] = f_4 * snf0_565[k]
                    - f_5 * snf1_565[k]
                    + f_3 * pc_x[k] * sng_845[k];

        t_1182[k] = pb_z[k] * smh0_951[k]
                    - f_8 * pc_z[k] * smh1_951[k];

        t_1183[k] = f_9 * smg_678[k]
                    + f_3 * pc_z[k] * sng_843[k];

        t_1184[k] = f_12 * smg_695[k]
                    + f_3 * pc_y[k] * sng_845[k];
    }

#pragma omp simd aligned(t_1185, t_1186, t_1187, t_1188, t_1189, t_1190, pc_x, snf0_569, \
                         snf1_569, sng_849, sng_850, sng_851, sng_852, sng_853, \
                         sng_854 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1185[k] = f_6 * snf0_569[k]
                    - f_7 * snf1_569[k]
                    + f_3 * pc_x[k] * sng_849[k];

        t_1186[k] = f_3 * pc_x[k] * sng_850[k];

        t_1187[k] = f_3 * pc_x[k] * sng_851[k];

        t_1188[k] = f_3 * pc_x[k] * sng_852[k];

        t_1189[k] = f_3 * pc_x[k] * sng_853[k];

        t_1190[k] = f_3 * pc_x[k] * sng_854[k];
    }

#pragma omp simd aligned(t_1191, t_1192, t_1193, t_1194, pb_z, pc_z, smh0_960, smh0_962, \
                         smh0_963, smg_685, smg_686, smg_687, smh1_960, smh1_962, smh1_963, \
                         sng_850 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1191[k] = pb_z[k] * smh0_960[k]
                    - f_8 * pc_z[k] * smh1_960[k];

        t_1192[k] = f_9 * smg_685[k]
                    + f_3 * pc_z[k] * sng_850[k];

        t_1193[k] = pb_z[k] * smh0_962[k]
                    + f_10 * smg_686[k]
                    - f_8 * pc_z[k] * smh1_962[k];

        t_1194[k] = pb_z[k] * smh0_963[k]
                    + f_11 * smg_687[k]
                    - f_8 * pc_z[k] * smh1_963[k];
    }

#pragma omp simd aligned(t_1195, t_1196, t_1197, t_1198, pc_x, pc_y, pc_z, smg_689, smg_704, \
                         smg_705, snf0_569, snf0_570, snf1_569, snf1_570, sng_854, \
                         sng_855 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1195[k] = f_12 * smg_704[k]
                    + f_3 * pc_y[k] * sng_854[k];

        t_1196[k] = f_9 * smg_689[k]
                    + f_1 * snf0_569[k]
                    - f_2 * snf1_569[k]
                    + f_3 * pc_z[k] * sng_854[k];

        t_1197[k] = f_1 * snf0_570[k]
                    - f_2 * snf1_570[k]
                    + f_3 * pc_x[k] * sng_855[k];

        t_1198[k] = f_13 * smg_705[k]
                    + f_3 * pc_y[k] * sng_855[k];
    }

#pragma omp simd aligned(t_1199, t_1200, t_1201, pc_x, pc_y, pc_z, smg_690, smg_707, snf0_573, \
                         snf1_573, sng_855, sng_857, sng_858 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1199[k] = f_10 * smg_690[k]
                    + f_3 * pc_z[k] * sng_855[k];

        t_1200[k] = f_4 * snf0_573[k]
                    - f_5 * snf1_573[k]
                    + f_3 * pc_x[k] * sng_858[k];

        t_1201[k] = f_13 * smg_707[k]
                    + f_3 * pc_y[k] * sng_857[k];
    }

#pragma omp simd aligned(t_1202, t_1203, t_1204, t_1205, pc_x, pc_y, pc_z, smg_693, smg_710, \
                         snf0_575, snf0_576, snf1_575, snf1_576, sng_858, sng_860, \
                         sng_861 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1202[k] = f_4 * snf0_575[k]
                    - f_5 * snf1_575[k]
                    + f_3 * pc_x[k] * sng_860[k];

        t_1203[k] = f_6 * snf0_576[k]
                    - f_7 * snf1_576[k]
                    + f_3 * pc_x[k] * sng_861[k];

        t_1204[k] = f_10 * smg_693[k]
                    + f_3 * pc_z[k] * sng_858[k];

        t_1205[k] = f_13 * smg_710[k]
                    + f_3 * pc_y[k] * sng_860[k];
    }

#pragma omp simd aligned(t_1206, t_1207, t_1208, t_1209, t_1210, t_1211, pc_x, snf0_579, \
                         snf1_579, sng_864, sng_865, sng_866, sng_867, sng_868, \
                         sng_869 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1206[k] = f_6 * snf0_579[k]
                    - f_7 * snf1_579[k]
                    + f_3 * pc_x[k] * sng_864[k];

        t_1207[k] = f_3 * pc_x[k] * sng_865[k];

        t_1208[k] = f_3 * pc_x[k] * sng_866[k];

        t_1209[k] = f_3 * pc_x[k] * sng_867[k];

        t_1210[k] = f_3 * pc_x[k] * sng_868[k];

        t_1211[k] = f_3 * pc_x[k] * sng_869[k];
    }

#pragma omp simd aligned(t_1212, t_1213, t_1214, pc_y, pc_z, smg_700, smg_715, smg_717, \
                         snf0_576, snf0_578, snf1_576, snf1_578, sng_865, \
                         sng_867 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1212[k] = f_13 * smg_715[k]
                    + f_1 * snf0_576[k]
                    - f_2 * snf1_576[k]
                    + f_3 * pc_y[k] * sng_865[k];

        t_1213[k] = f_10 * smg_700[k]
                    + f_3 * pc_z[k] * sng_865[k];

        t_1214[k] = f_13 * smg_717[k]
                    + f_4 * snf0_578[k]
                    - f_5 * snf1_578[k]
                    + f_3 * pc_y[k] * sng_867[k];
    }

#pragma omp simd aligned(t_1215, t_1216, t_1217, pc_y, pc_z, smg_704, smg_718, smg_719, \
                         snf0_579, snf1_579, sng_868, sng_869 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1215[k] = f_13 * smg_718[k]
                    + f_6 * snf0_579[k]
                    - f_7 * snf1_579[k]
                    + f_3 * pc_y[k] * sng_868[k];

        t_1216[k] = f_13 * smg_719[k]
                    + f_3 * pc_y[k] * sng_869[k];

        t_1217[k] = f_10 * smg_704[k]
                    + f_1 * snf0_579[k]
                    - f_2 * snf1_579[k]
                    + f_3 * pc_z[k] * sng_869[k];
    }

#pragma omp simd aligned(t_1218, t_1219, t_1220, t_1221, pc_x, pc_y, pc_z, smg_705, smg_720, \
                         snf0_580, snf0_583, snf1_580, snf1_583, sng_870, \
                         sng_873 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1218[k] = f_1 * snf0_580[k]
                    - f_2 * snf1_580[k]
                    + f_3 * pc_x[k] * sng_870[k];

        t_1219[k] = f_14 * smg_720[k]
                    + f_3 * pc_y[k] * sng_870[k];

        t_1220[k] = f_11 * smg_705[k]
                    + f_3 * pc_z[k] * sng_870[k];

        t_1221[k] = f_4 * snf0_583[k]
                    - f_5 * snf1_583[k]
                    + f_3 * pc_x[k] * sng_873[k];
    }

#pragma omp simd aligned(t_1222, t_1223, t_1224, pc_x, pc_y, smg_722, snf0_585, snf0_586, \
                         snf1_585, snf1_586, sng_872, sng_875, \
                         sng_876 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1222[k] = f_14 * smg_722[k]
                    + f_3 * pc_y[k] * sng_872[k];

        t_1223[k] = f_4 * snf0_585[k]
                    - f_5 * snf1_585[k]
                    + f_3 * pc_x[k] * sng_875[k];

        t_1224[k] = f_6 * snf0_586[k]
                    - f_7 * snf1_586[k]
                    + f_3 * pc_x[k] * sng_876[k];
    }

#pragma omp simd aligned(t_1225, t_1226, t_1227, t_1228, pc_x, pc_y, pc_z, smg_708, smg_725, \
                         snf0_589, snf1_589, sng_873, sng_875, sng_879, \
                         sng_880 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1225[k] = f_11 * smg_708[k]
                    + f_3 * pc_z[k] * sng_873[k];

        t_1226[k] = f_14 * smg_725[k]
                    + f_3 * pc_y[k] * sng_875[k];

        t_1227[k] = f_6 * snf0_589[k]
                    - f_7 * snf1_589[k]
                    + f_3 * pc_x[k] * sng_879[k];

        t_1228[k] = f_3 * pc_x[k] * sng_880[k];
    }

#pragma omp simd aligned(t_1229, t_1230, t_1231, t_1232, t_1233, pc_x, pc_y, smg_730, \
                         snf0_586, snf1_586, sng_880, sng_881, sng_882, sng_883, \
                         sng_884 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1229[k] = f_3 * pc_x[k] * sng_881[k];

        t_1230[k] = f_3 * pc_x[k] * sng_882[k];

        t_1231[k] = f_3 * pc_x[k] * sng_883[k];

        t_1232[k] = f_3 * pc_x[k] * sng_884[k];

        t_1233[k] = f_14 * smg_730[k]
                    + f_1 * snf0_586[k]
                    - f_2 * snf1_586[k]
                    + f_3 * pc_y[k] * sng_880[k];
    }

#pragma omp simd aligned(t_1234, t_1235, t_1236, pc_y, pc_z, smg_715, smg_732, smg_733, \
                         snf0_588, snf0_589, snf1_588, snf1_589, sng_880, sng_882, \
                         sng_883 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1234[k] = f_11 * smg_715[k]
                    + f_3 * pc_z[k] * sng_880[k];

        t_1235[k] = f_14 * smg_732[k]
                    + f_4 * snf0_588[k]
                    - f_5 * snf1_588[k]
                    + f_3 * pc_y[k] * sng_882[k];

        t_1236[k] = f_14 * smg_733[k]
                    + f_6 * snf0_589[k]
                    - f_7 * snf1_589[k]
                    + f_3 * pc_y[k] * sng_883[k];
    }

#pragma omp simd aligned(t_1237, t_1238, t_1239, t_1240, pc_x, pc_y, pc_z, smg_719, smg_734, \
                         smg_735, snf0_589, snf0_590, snf1_589, snf1_590, sng_884, \
                         sng_885 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1237[k] = f_14 * smg_734[k]
                    + f_3 * pc_y[k] * sng_884[k];

        t_1238[k] = f_11 * smg_719[k]
                    + f_1 * snf0_589[k]
                    - f_2 * snf1_589[k]
                    + f_3 * pc_z[k] * sng_884[k];

        t_1239[k] = f_1 * snf0_590[k]
                    - f_2 * snf1_590[k]
                    + f_3 * pc_x[k] * sng_885[k];

        t_1240[k] = f_15 * smg_735[k]
                    + f_3 * pc_y[k] * sng_885[k];
    }

#pragma omp simd aligned(t_1241, t_1242, t_1243, pc_x, pc_y, pc_z, smg_720, smg_737, snf0_593, \
                         snf1_593, sng_885, sng_887, sng_888 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1241[k] = f_16 * smg_720[k]
                    + f_3 * pc_z[k] * sng_885[k];

        t_1242[k] = f_4 * snf0_593[k]
                    - f_5 * snf1_593[k]
                    + f_3 * pc_x[k] * sng_888[k];

        t_1243[k] = f_15 * smg_737[k]
                    + f_3 * pc_y[k] * sng_887[k];
    }

#pragma omp simd aligned(t_1244, t_1245, t_1246, t_1247, pc_x, pc_y, pc_z, smg_723, smg_740, \
                         snf0_595, snf0_596, snf1_595, snf1_596, sng_888, sng_890, \
                         sng_891 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1244[k] = f_4 * snf0_595[k]
                    - f_5 * snf1_595[k]
                    + f_3 * pc_x[k] * sng_890[k];

        t_1245[k] = f_6 * snf0_596[k]
                    - f_7 * snf1_596[k]
                    + f_3 * pc_x[k] * sng_891[k];

        t_1246[k] = f_16 * smg_723[k]
                    + f_3 * pc_z[k] * sng_888[k];

        t_1247[k] = f_15 * smg_740[k]
                    + f_3 * pc_y[k] * sng_890[k];
    }

#pragma omp simd aligned(t_1248, t_1249, t_1250, t_1251, t_1252, t_1253, pc_x, snf0_599, \
                         snf1_599, sng_894, sng_895, sng_896, sng_897, sng_898, \
                         sng_899 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1248[k] = f_6 * snf0_599[k]
                    - f_7 * snf1_599[k]
                    + f_3 * pc_x[k] * sng_894[k];

        t_1249[k] = f_3 * pc_x[k] * sng_895[k];

        t_1250[k] = f_3 * pc_x[k] * sng_896[k];

        t_1251[k] = f_3 * pc_x[k] * sng_897[k];

        t_1252[k] = f_3 * pc_x[k] * sng_898[k];

        t_1253[k] = f_3 * pc_x[k] * sng_899[k];
    }

#pragma omp simd aligned(t_1254, t_1255, t_1256, pc_y, pc_z, smg_730, smg_745, smg_747, \
                         snf0_596, snf0_598, snf1_596, snf1_598, sng_895, \
                         sng_897 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1254[k] = f_15 * smg_745[k]
                    + f_1 * snf0_596[k]
                    - f_2 * snf1_596[k]
                    + f_3 * pc_y[k] * sng_895[k];

        t_1255[k] = f_16 * smg_730[k]
                    + f_3 * pc_z[k] * sng_895[k];

        t_1256[k] = f_15 * smg_747[k]
                    + f_4 * snf0_598[k]
                    - f_5 * snf1_598[k]
                    + f_3 * pc_y[k] * sng_897[k];
    }

#pragma omp simd aligned(t_1257, t_1258, t_1259, pc_y, pc_z, smg_734, smg_748, smg_749, \
                         snf0_599, snf1_599, sng_898, sng_899 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1257[k] = f_15 * smg_748[k]
                    + f_6 * snf0_599[k]
                    - f_7 * snf1_599[k]
                    + f_3 * pc_y[k] * sng_898[k];

        t_1258[k] = f_15 * smg_749[k]
                    + f_3 * pc_y[k] * sng_899[k];

        t_1259[k] = f_16 * smg_734[k]
                    + f_1 * snf0_599[k]
                    - f_2 * snf1_599[k]
                    + f_3 * pc_z[k] * sng_899[k];
    }

#pragma omp simd aligned(t_1260, t_1261, t_1262, t_1263, pc_x, pc_y, pc_z, smg_735, smg_750, \
                         snf0_600, snf0_603, snf1_600, snf1_603, sng_900, \
                         sng_903 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1260[k] = f_1 * snf0_600[k]
                    - f_2 * snf1_600[k]
                    + f_3 * pc_x[k] * sng_900[k];

        t_1261[k] = f_17 * smg_750[k]
                    + f_3 * pc_y[k] * sng_900[k];

        t_1262[k] = f_17 * smg_735[k]
                    + f_3 * pc_z[k] * sng_900[k];

        t_1263[k] = f_4 * snf0_603[k]
                    - f_5 * snf1_603[k]
                    + f_3 * pc_x[k] * sng_903[k];
    }

#pragma omp simd aligned(t_1264, t_1265, t_1266, pc_x, pc_y, smg_752, snf0_605, snf0_606, \
                         snf1_605, snf1_606, sng_902, sng_905, \
                         sng_906 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1264[k] = f_17 * smg_752[k]
                    + f_3 * pc_y[k] * sng_902[k];

        t_1265[k] = f_4 * snf0_605[k]
                    - f_5 * snf1_605[k]
                    + f_3 * pc_x[k] * sng_905[k];

        t_1266[k] = f_6 * snf0_606[k]
                    - f_7 * snf1_606[k]
                    + f_3 * pc_x[k] * sng_906[k];
    }

#pragma omp simd aligned(t_1267, t_1268, t_1269, t_1270, pc_x, pc_y, pc_z, smg_738, smg_755, \
                         snf0_609, snf1_609, sng_903, sng_905, sng_909, \
                         sng_910 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1267[k] = f_17 * smg_738[k]
                    + f_3 * pc_z[k] * sng_903[k];

        t_1268[k] = f_17 * smg_755[k]
                    + f_3 * pc_y[k] * sng_905[k];

        t_1269[k] = f_6 * snf0_609[k]
                    - f_7 * snf1_609[k]
                    + f_3 * pc_x[k] * sng_909[k];

        t_1270[k] = f_3 * pc_x[k] * sng_910[k];
    }

#pragma omp simd aligned(t_1271, t_1272, t_1273, t_1274, t_1275, pc_x, pc_y, smg_760, \
                         snf0_606, snf1_606, sng_910, sng_911, sng_912, sng_913, \
                         sng_914 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1271[k] = f_3 * pc_x[k] * sng_911[k];

        t_1272[k] = f_3 * pc_x[k] * sng_912[k];

        t_1273[k] = f_3 * pc_x[k] * sng_913[k];

        t_1274[k] = f_3 * pc_x[k] * sng_914[k];

        t_1275[k] = f_17 * smg_760[k]
                    + f_1 * snf0_606[k]
                    - f_2 * snf1_606[k]
                    + f_3 * pc_y[k] * sng_910[k];
    }

#pragma omp simd aligned(t_1276, t_1277, t_1278, pc_y, pc_z, smg_745, smg_762, smg_763, \
                         snf0_608, snf0_609, snf1_608, snf1_609, sng_910, sng_912, \
                         sng_913 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1276[k] = f_17 * smg_745[k]
                    + f_3 * pc_z[k] * sng_910[k];

        t_1277[k] = f_17 * smg_762[k]
                    + f_4 * snf0_608[k]
                    - f_5 * snf1_608[k]
                    + f_3 * pc_y[k] * sng_912[k];

        t_1278[k] = f_17 * smg_763[k]
                    + f_6 * snf0_609[k]
                    - f_7 * snf1_609[k]
                    + f_3 * pc_y[k] * sng_913[k];
    }

#pragma omp simd aligned(t_1279, t_1280, t_1281, t_1282, pc_x, pc_y, pc_z, smg_749, smg_764, \
                         smg_765, snf0_609, snf0_610, snf1_609, snf1_610, sng_914, \
                         sng_915 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1279[k] = f_17 * smg_764[k]
                    + f_3 * pc_y[k] * sng_914[k];

        t_1280[k] = f_17 * smg_749[k]
                    + f_1 * snf0_609[k]
                    - f_2 * snf1_609[k]
                    + f_3 * pc_z[k] * sng_914[k];

        t_1281[k] = f_1 * snf0_610[k]
                    - f_2 * snf1_610[k]
                    + f_3 * pc_x[k] * sng_915[k];

        t_1282[k] = f_16 * smg_765[k]
                    + f_3 * pc_y[k] * sng_915[k];
    }

#pragma omp simd aligned(t_1283, t_1284, t_1285, pc_x, pc_y, pc_z, smg_750, smg_767, snf0_613, \
                         snf1_613, sng_915, sng_917, sng_918 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1283[k] = f_15 * smg_750[k]
                    + f_3 * pc_z[k] * sng_915[k];

        t_1284[k] = f_4 * snf0_613[k]
                    - f_5 * snf1_613[k]
                    + f_3 * pc_x[k] * sng_918[k];

        t_1285[k] = f_16 * smg_767[k]
                    + f_3 * pc_y[k] * sng_917[k];
    }

#pragma omp simd aligned(t_1286, t_1287, t_1288, t_1289, pc_x, pc_y, pc_z, smg_753, smg_770, \
                         snf0_615, snf0_616, snf1_615, snf1_616, sng_918, sng_920, \
                         sng_921 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1286[k] = f_4 * snf0_615[k]
                    - f_5 * snf1_615[k]
                    + f_3 * pc_x[k] * sng_920[k];

        t_1287[k] = f_6 * snf0_616[k]
                    - f_7 * snf1_616[k]
                    + f_3 * pc_x[k] * sng_921[k];

        t_1288[k] = f_15 * smg_753[k]
                    + f_3 * pc_z[k] * sng_918[k];

        t_1289[k] = f_16 * smg_770[k]
                    + f_3 * pc_y[k] * sng_920[k];
    }
}

static auto
compute_prim_snh_three_center_electron_repulsion_0_piece11(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t smh0,
                                                           const size_t smg, const size_t smh1,
                                                           const size_t snf0, const size_t snf1,
                                                           const size_t sng, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.0 / q;
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
    const auto f_12 = 4.5 / q;
    const auto f_13 = 4.0 / q;
    const auto f_14 = 3.5 / q;
    const auto f_15 = 3.0 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *smh0_1134 = buffer.data(smh0 + 1134);
    const auto *smh0_1139 = buffer.data(smh0 + 1139);
    const auto *smh0_1143 = buffer.data(smh0 + 1143);
    const auto *smh0_1149 = buffer.data(smh0 + 1149);
    const auto *smh0_1151 = buffer.data(smh0 + 1151);
    const auto *smh0_1152 = buffer.data(smh0 + 1152);
    const auto *smh0_1154 = buffer.data(smh0 + 1154);

    const auto *smg_760 = buffer.data(smg + 760);
    const auto *smg_764 = buffer.data(smg + 764);
    const auto *smg_765 = buffer.data(smg + 765);
    const auto *smg_768 = buffer.data(smg + 768);
    const auto *smg_775 = buffer.data(smg + 775);
    const auto *smg_777 = buffer.data(smg + 777);
    const auto *smg_778 = buffer.data(smg + 778);
    const auto *smg_779 = buffer.data(smg + 779);
    const auto *smg_780 = buffer.data(smg + 780);
    const auto *smg_782 = buffer.data(smg + 782);
    const auto *smg_783 = buffer.data(smg + 783);
    const auto *smg_785 = buffer.data(smg + 785);
    const auto *smg_790 = buffer.data(smg + 790);
    const auto *smg_792 = buffer.data(smg + 792);
    const auto *smg_793 = buffer.data(smg + 793);
    const auto *smg_794 = buffer.data(smg + 794);
    const auto *smg_795 = buffer.data(smg + 795);
    const auto *smg_797 = buffer.data(smg + 797);
    const auto *smg_798 = buffer.data(smg + 798);
    const auto *smg_800 = buffer.data(smg + 800);
    const auto *smg_805 = buffer.data(smg + 805);
    const auto *smg_807 = buffer.data(smg + 807);
    const auto *smg_808 = buffer.data(smg + 808);
    const auto *smg_809 = buffer.data(smg + 809);
    const auto *smg_810 = buffer.data(smg + 810);
    const auto *smg_812 = buffer.data(smg + 812);
    const auto *smg_813 = buffer.data(smg + 813);
    const auto *smg_815 = buffer.data(smg + 815);
    const auto *smg_820 = buffer.data(smg + 820);
    const auto *smg_822 = buffer.data(smg + 822);
    const auto *smg_823 = buffer.data(smg + 823);
    const auto *smg_824 = buffer.data(smg + 824);

    const auto *smh1_1134 = buffer.data(smh1 + 1134);
    const auto *smh1_1139 = buffer.data(smh1 + 1139);
    const auto *smh1_1143 = buffer.data(smh1 + 1143);
    const auto *smh1_1149 = buffer.data(smh1 + 1149);
    const auto *smh1_1151 = buffer.data(smh1 + 1151);
    const auto *smh1_1152 = buffer.data(smh1 + 1152);
    const auto *smh1_1154 = buffer.data(smh1 + 1154);

    const auto *snf0_616 = buffer.data(snf0 + 616);
    const auto *snf0_618 = buffer.data(snf0 + 618);
    const auto *snf0_619 = buffer.data(snf0 + 619);
    const auto *snf0_620 = buffer.data(snf0 + 620);
    const auto *snf0_623 = buffer.data(snf0 + 623);
    const auto *snf0_625 = buffer.data(snf0 + 625);
    const auto *snf0_626 = buffer.data(snf0 + 626);
    const auto *snf0_628 = buffer.data(snf0 + 628);
    const auto *snf0_629 = buffer.data(snf0 + 629);
    const auto *snf0_630 = buffer.data(snf0 + 630);
    const auto *snf0_633 = buffer.data(snf0 + 633);
    const auto *snf0_635 = buffer.data(snf0 + 635);
    const auto *snf0_636 = buffer.data(snf0 + 636);
    const auto *snf0_638 = buffer.data(snf0 + 638);
    const auto *snf0_639 = buffer.data(snf0 + 639);
    const auto *snf0_643 = buffer.data(snf0 + 643);
    const auto *snf0_646 = buffer.data(snf0 + 646);
    const auto *snf0_650 = buffer.data(snf0 + 650);
    const auto *snf0_653 = buffer.data(snf0 + 653);
    const auto *snf0_655 = buffer.data(snf0 + 655);
    const auto *snf0_656 = buffer.data(snf0 + 656);
    const auto *snf0_658 = buffer.data(snf0 + 658);
    const auto *snf0_659 = buffer.data(snf0 + 659);

    const auto *snf1_616 = buffer.data(snf1 + 616);
    const auto *snf1_618 = buffer.data(snf1 + 618);
    const auto *snf1_619 = buffer.data(snf1 + 619);
    const auto *snf1_620 = buffer.data(snf1 + 620);
    const auto *snf1_623 = buffer.data(snf1 + 623);
    const auto *snf1_625 = buffer.data(snf1 + 625);
    const auto *snf1_626 = buffer.data(snf1 + 626);
    const auto *snf1_628 = buffer.data(snf1 + 628);
    const auto *snf1_629 = buffer.data(snf1 + 629);
    const auto *snf1_630 = buffer.data(snf1 + 630);
    const auto *snf1_633 = buffer.data(snf1 + 633);
    const auto *snf1_635 = buffer.data(snf1 + 635);
    const auto *snf1_636 = buffer.data(snf1 + 636);
    const auto *snf1_638 = buffer.data(snf1 + 638);
    const auto *snf1_639 = buffer.data(snf1 + 639);
    const auto *snf1_643 = buffer.data(snf1 + 643);
    const auto *snf1_646 = buffer.data(snf1 + 646);
    const auto *snf1_650 = buffer.data(snf1 + 650);
    const auto *snf1_653 = buffer.data(snf1 + 653);
    const auto *snf1_655 = buffer.data(snf1 + 655);
    const auto *snf1_656 = buffer.data(snf1 + 656);
    const auto *snf1_658 = buffer.data(snf1 + 658);
    const auto *snf1_659 = buffer.data(snf1 + 659);

    const auto *sng_924 = buffer.data(sng + 924);
    const auto *sng_925 = buffer.data(sng + 925);
    const auto *sng_926 = buffer.data(sng + 926);
    const auto *sng_927 = buffer.data(sng + 927);
    const auto *sng_928 = buffer.data(sng + 928);
    const auto *sng_929 = buffer.data(sng + 929);
    const auto *sng_930 = buffer.data(sng + 930);
    const auto *sng_932 = buffer.data(sng + 932);
    const auto *sng_933 = buffer.data(sng + 933);
    const auto *sng_935 = buffer.data(sng + 935);
    const auto *sng_936 = buffer.data(sng + 936);
    const auto *sng_939 = buffer.data(sng + 939);
    const auto *sng_940 = buffer.data(sng + 940);
    const auto *sng_941 = buffer.data(sng + 941);
    const auto *sng_942 = buffer.data(sng + 942);
    const auto *sng_943 = buffer.data(sng + 943);
    const auto *sng_944 = buffer.data(sng + 944);
    const auto *sng_945 = buffer.data(sng + 945);
    const auto *sng_947 = buffer.data(sng + 947);
    const auto *sng_948 = buffer.data(sng + 948);
    const auto *sng_950 = buffer.data(sng + 950);
    const auto *sng_951 = buffer.data(sng + 951);
    const auto *sng_954 = buffer.data(sng + 954);
    const auto *sng_955 = buffer.data(sng + 955);
    const auto *sng_956 = buffer.data(sng + 956);
    const auto *sng_957 = buffer.data(sng + 957);
    const auto *sng_958 = buffer.data(sng + 958);
    const auto *sng_959 = buffer.data(sng + 959);
    const auto *sng_960 = buffer.data(sng + 960);
    const auto *sng_962 = buffer.data(sng + 962);
    const auto *sng_963 = buffer.data(sng + 963);
    const auto *sng_965 = buffer.data(sng + 965);
    const auto *sng_966 = buffer.data(sng + 966);
    const auto *sng_970 = buffer.data(sng + 970);
    const auto *sng_971 = buffer.data(sng + 971);
    const auto *sng_972 = buffer.data(sng + 972);
    const auto *sng_973 = buffer.data(sng + 973);
    const auto *sng_974 = buffer.data(sng + 974);
    const auto *sng_975 = buffer.data(sng + 975);
    const auto *sng_977 = buffer.data(sng + 977);
    const auto *sng_978 = buffer.data(sng + 978);
    const auto *sng_980 = buffer.data(sng + 980);
    const auto *sng_981 = buffer.data(sng + 981);
    const auto *sng_984 = buffer.data(sng + 984);
    const auto *sng_985 = buffer.data(sng + 985);
    const auto *sng_986 = buffer.data(sng + 986);
    const auto *sng_987 = buffer.data(sng + 987);
    const auto *sng_988 = buffer.data(sng + 988);
    const auto *sng_989 = buffer.data(sng + 989);

#pragma omp simd aligned(t_1290, t_1291, t_1292, t_1293, t_1294, t_1295, pc_x, snf0_619, \
                         snf1_619, sng_924, sng_925, sng_926, sng_927, sng_928, \
                         sng_929 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1290[k] = f_6 * snf0_619[k]
                    - f_7 * snf1_619[k]
                    + f_3 * pc_x[k] * sng_924[k];

        t_1291[k] = f_3 * pc_x[k] * sng_925[k];

        t_1292[k] = f_3 * pc_x[k] * sng_926[k];

        t_1293[k] = f_3 * pc_x[k] * sng_927[k];

        t_1294[k] = f_3 * pc_x[k] * sng_928[k];

        t_1295[k] = f_3 * pc_x[k] * sng_929[k];
    }

#pragma omp simd aligned(t_1296, t_1297, t_1298, pc_y, pc_z, smg_760, smg_775, smg_777, \
                         snf0_616, snf0_618, snf1_616, snf1_618, sng_925, \
                         sng_927 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1296[k] = f_16 * smg_775[k]
                    + f_1 * snf0_616[k]
                    - f_2 * snf1_616[k]
                    + f_3 * pc_y[k] * sng_925[k];

        t_1297[k] = f_15 * smg_760[k]
                    + f_3 * pc_z[k] * sng_925[k];

        t_1298[k] = f_16 * smg_777[k]
                    + f_4 * snf0_618[k]
                    - f_5 * snf1_618[k]
                    + f_3 * pc_y[k] * sng_927[k];
    }

#pragma omp simd aligned(t_1299, t_1300, t_1301, pc_y, pc_z, smg_764, smg_778, smg_779, \
                         snf0_619, snf1_619, sng_928, sng_929 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1299[k] = f_16 * smg_778[k]
                    + f_6 * snf0_619[k]
                    - f_7 * snf1_619[k]
                    + f_3 * pc_y[k] * sng_928[k];

        t_1300[k] = f_16 * smg_779[k]
                    + f_3 * pc_y[k] * sng_929[k];

        t_1301[k] = f_15 * smg_764[k]
                    + f_1 * snf0_619[k]
                    - f_2 * snf1_619[k]
                    + f_3 * pc_z[k] * sng_929[k];
    }

#pragma omp simd aligned(t_1302, t_1303, t_1304, t_1305, pc_x, pc_y, pc_z, smg_765, smg_780, \
                         snf0_620, snf0_623, snf1_620, snf1_623, sng_930, \
                         sng_933 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1302[k] = f_1 * snf0_620[k]
                    - f_2 * snf1_620[k]
                    + f_3 * pc_x[k] * sng_930[k];

        t_1303[k] = f_11 * smg_780[k]
                    + f_3 * pc_y[k] * sng_930[k];

        t_1304[k] = f_14 * smg_765[k]
                    + f_3 * pc_z[k] * sng_930[k];

        t_1305[k] = f_4 * snf0_623[k]
                    - f_5 * snf1_623[k]
                    + f_3 * pc_x[k] * sng_933[k];
    }

#pragma omp simd aligned(t_1306, t_1307, t_1308, pc_x, pc_y, smg_782, snf0_625, snf0_626, \
                         snf1_625, snf1_626, sng_932, sng_935, \
                         sng_936 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1306[k] = f_11 * smg_782[k]
                    + f_3 * pc_y[k] * sng_932[k];

        t_1307[k] = f_4 * snf0_625[k]
                    - f_5 * snf1_625[k]
                    + f_3 * pc_x[k] * sng_935[k];

        t_1308[k] = f_6 * snf0_626[k]
                    - f_7 * snf1_626[k]
                    + f_3 * pc_x[k] * sng_936[k];
    }

#pragma omp simd aligned(t_1309, t_1310, t_1311, t_1312, pc_x, pc_y, pc_z, smg_768, smg_785, \
                         snf0_629, snf1_629, sng_933, sng_935, sng_939, \
                         sng_940 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1309[k] = f_14 * smg_768[k]
                    + f_3 * pc_z[k] * sng_933[k];

        t_1310[k] = f_11 * smg_785[k]
                    + f_3 * pc_y[k] * sng_935[k];

        t_1311[k] = f_6 * snf0_629[k]
                    - f_7 * snf1_629[k]
                    + f_3 * pc_x[k] * sng_939[k];

        t_1312[k] = f_3 * pc_x[k] * sng_940[k];
    }

#pragma omp simd aligned(t_1313, t_1314, t_1315, t_1316, t_1317, pc_x, pc_y, smg_790, \
                         snf0_626, snf1_626, sng_940, sng_941, sng_942, sng_943, \
                         sng_944 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1313[k] = f_3 * pc_x[k] * sng_941[k];

        t_1314[k] = f_3 * pc_x[k] * sng_942[k];

        t_1315[k] = f_3 * pc_x[k] * sng_943[k];

        t_1316[k] = f_3 * pc_x[k] * sng_944[k];

        t_1317[k] = f_11 * smg_790[k]
                    + f_1 * snf0_626[k]
                    - f_2 * snf1_626[k]
                    + f_3 * pc_y[k] * sng_940[k];
    }

#pragma omp simd aligned(t_1318, t_1319, t_1320, pc_y, pc_z, smg_775, smg_792, smg_793, \
                         snf0_628, snf0_629, snf1_628, snf1_629, sng_940, sng_942, \
                         sng_943 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1318[k] = f_14 * smg_775[k]
                    + f_3 * pc_z[k] * sng_940[k];

        t_1319[k] = f_11 * smg_792[k]
                    + f_4 * snf0_628[k]
                    - f_5 * snf1_628[k]
                    + f_3 * pc_y[k] * sng_942[k];

        t_1320[k] = f_11 * smg_793[k]
                    + f_6 * snf0_629[k]
                    - f_7 * snf1_629[k]
                    + f_3 * pc_y[k] * sng_943[k];
    }

#pragma omp simd aligned(t_1321, t_1322, t_1323, t_1324, pc_x, pc_y, pc_z, smg_779, smg_794, \
                         smg_795, snf0_629, snf0_630, snf1_629, snf1_630, sng_944, \
                         sng_945 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1321[k] = f_11 * smg_794[k]
                    + f_3 * pc_y[k] * sng_944[k];

        t_1322[k] = f_14 * smg_779[k]
                    + f_1 * snf0_629[k]
                    - f_2 * snf1_629[k]
                    + f_3 * pc_z[k] * sng_944[k];

        t_1323[k] = f_1 * snf0_630[k]
                    - f_2 * snf1_630[k]
                    + f_3 * pc_x[k] * sng_945[k];

        t_1324[k] = f_10 * smg_795[k]
                    + f_3 * pc_y[k] * sng_945[k];
    }

#pragma omp simd aligned(t_1325, t_1326, t_1327, pc_x, pc_y, pc_z, smg_780, smg_797, snf0_633, \
                         snf1_633, sng_945, sng_947, sng_948 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1325[k] = f_13 * smg_780[k]
                    + f_3 * pc_z[k] * sng_945[k];

        t_1326[k] = f_4 * snf0_633[k]
                    - f_5 * snf1_633[k]
                    + f_3 * pc_x[k] * sng_948[k];

        t_1327[k] = f_10 * smg_797[k]
                    + f_3 * pc_y[k] * sng_947[k];
    }

#pragma omp simd aligned(t_1328, t_1329, t_1330, t_1331, pc_x, pc_y, pc_z, smg_783, smg_800, \
                         snf0_635, snf0_636, snf1_635, snf1_636, sng_948, sng_950, \
                         sng_951 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1328[k] = f_4 * snf0_635[k]
                    - f_5 * snf1_635[k]
                    + f_3 * pc_x[k] * sng_950[k];

        t_1329[k] = f_6 * snf0_636[k]
                    - f_7 * snf1_636[k]
                    + f_3 * pc_x[k] * sng_951[k];

        t_1330[k] = f_13 * smg_783[k]
                    + f_3 * pc_z[k] * sng_948[k];

        t_1331[k] = f_10 * smg_800[k]
                    + f_3 * pc_y[k] * sng_950[k];
    }

#pragma omp simd aligned(t_1332, t_1333, t_1334, t_1335, t_1336, t_1337, pc_x, snf0_639, \
                         snf1_639, sng_954, sng_955, sng_956, sng_957, sng_958, \
                         sng_959 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1332[k] = f_6 * snf0_639[k]
                    - f_7 * snf1_639[k]
                    + f_3 * pc_x[k] * sng_954[k];

        t_1333[k] = f_3 * pc_x[k] * sng_955[k];

        t_1334[k] = f_3 * pc_x[k] * sng_956[k];

        t_1335[k] = f_3 * pc_x[k] * sng_957[k];

        t_1336[k] = f_3 * pc_x[k] * sng_958[k];

        t_1337[k] = f_3 * pc_x[k] * sng_959[k];
    }

#pragma omp simd aligned(t_1338, t_1339, t_1340, pc_y, pc_z, smg_790, smg_805, smg_807, \
                         snf0_636, snf0_638, snf1_636, snf1_638, sng_955, \
                         sng_957 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1338[k] = f_10 * smg_805[k]
                    + f_1 * snf0_636[k]
                    - f_2 * snf1_636[k]
                    + f_3 * pc_y[k] * sng_955[k];

        t_1339[k] = f_13 * smg_790[k]
                    + f_3 * pc_z[k] * sng_955[k];

        t_1340[k] = f_10 * smg_807[k]
                    + f_4 * snf0_638[k]
                    - f_5 * snf1_638[k]
                    + f_3 * pc_y[k] * sng_957[k];
    }

#pragma omp simd aligned(t_1341, t_1342, t_1343, t_1344, pb_y, pc_y, pc_z, smh0_1134, smg_794, \
                         smg_808, smg_809, smh1_1134, snf0_639, snf1_639, sng_958, \
                         sng_959 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1341[k] = f_10 * smg_808[k]
                    + f_6 * snf0_639[k]
                    - f_7 * snf1_639[k]
                    + f_3 * pc_y[k] * sng_958[k];

        t_1342[k] = f_10 * smg_809[k]
                    + f_3 * pc_y[k] * sng_959[k];

        t_1343[k] = f_13 * smg_794[k]
                    + f_1 * snf0_639[k]
                    - f_2 * snf1_639[k]
                    + f_3 * pc_z[k] * sng_959[k];

        t_1344[k] = pb_y[k] * smh0_1134[k]
                    - f_8 * pc_y[k] * smh1_1134[k];
    }

#pragma omp simd aligned(t_1345, t_1346, t_1347, t_1348, pc_x, pc_y, pc_z, smg_795, smg_810, \
                         smg_812, snf0_643, snf1_643, sng_960, sng_962, \
                         sng_963 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1345[k] = f_9 * smg_810[k]
                    + f_3 * pc_y[k] * sng_960[k];

        t_1346[k] = f_12 * smg_795[k]
                    + f_3 * pc_z[k] * sng_960[k];

        t_1347[k] = f_4 * snf0_643[k]
                    - f_5 * snf1_643[k]
                    + f_3 * pc_x[k] * sng_963[k];

        t_1348[k] = f_9 * smg_812[k]
                    + f_3 * pc_y[k] * sng_962[k];
    }

#pragma omp simd aligned(t_1349, t_1350, t_1351, pb_y, pc_x, pc_y, pc_z, smh0_1139, smg_798, \
                         smh1_1139, snf0_646, snf1_646, sng_963, \
                         sng_966 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1349[k] = pb_y[k] * smh0_1139[k]
                    - f_8 * pc_y[k] * smh1_1139[k];

        t_1350[k] = f_6 * snf0_646[k]
                    - f_7 * snf1_646[k]
                    + f_3 * pc_x[k] * sng_966[k];

        t_1351[k] = f_12 * smg_798[k]
                    + f_3 * pc_z[k] * sng_963[k];
    }

#pragma omp simd aligned(t_1352, t_1353, t_1354, t_1355, t_1356, pb_y, pc_x, pc_y, smh0_1143, \
                         smg_815, smh1_1143, sng_965, sng_970, sng_971, \
                         sng_972 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1352[k] = f_9 * smg_815[k]
                    + f_3 * pc_y[k] * sng_965[k];

        t_1353[k] = pb_y[k] * smh0_1143[k]
                    - f_8 * pc_y[k] * smh1_1143[k];

        t_1354[k] = f_3 * pc_x[k] * sng_970[k];

        t_1355[k] = f_3 * pc_x[k] * sng_971[k];

        t_1356[k] = f_3 * pc_x[k] * sng_972[k];
    }

#pragma omp simd aligned(t_1357, t_1358, t_1359, t_1360, pb_y, pc_x, pc_y, pc_z, smh0_1149, \
                         smg_805, smg_820, smh1_1149, sng_970, sng_973, \
                         sng_974 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1357[k] = f_3 * pc_x[k] * sng_973[k];

        t_1358[k] = f_3 * pc_x[k] * sng_974[k];

        t_1359[k] = pb_y[k] * smh0_1149[k]
                    + f_17 * smg_820[k]
                    - f_8 * pc_y[k] * smh1_1149[k];

        t_1360[k] = f_12 * smg_805[k]
                    + f_3 * pc_z[k] * sng_970[k];
    }

#pragma omp simd aligned(t_1361, t_1362, t_1363, t_1364, pb_y, pc_y, smh0_1151, smh0_1152, \
                         smh0_1154, smg_822, smg_823, smg_824, smh1_1151, smh1_1152, \
                         smh1_1154, sng_974 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1361[k] = pb_y[k] * smh0_1151[k]
                    + f_11 * smg_822[k]
                    - f_8 * pc_y[k] * smh1_1151[k];

        t_1362[k] = pb_y[k] * smh0_1152[k]
                    + f_10 * smg_823[k]
                    - f_8 * pc_y[k] * smh1_1152[k];

        t_1363[k] = f_9 * smg_824[k]
                    + f_3 * pc_y[k] * sng_974[k];

        t_1364[k] = pb_y[k] * smh0_1154[k]
                    - f_8 * pc_y[k] * smh1_1154[k];
    }

#pragma omp simd aligned(t_1365, t_1366, t_1367, t_1368, t_1369, pc_x, pc_y, pc_z, smg_810, \
                         snf0_650, snf0_653, snf1_650, snf1_653, sng_975, sng_977, \
                         sng_978 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1365[k] = f_1 * snf0_650[k]
                    - f_2 * snf1_650[k]
                    + f_3 * pc_x[k] * sng_975[k];

        t_1366[k] = f_3 * pc_y[k] * sng_975[k];

        t_1367[k] = f_0 * smg_810[k]
                    + f_3 * pc_z[k] * sng_975[k];

        t_1368[k] = f_4 * snf0_653[k]
                    - f_5 * snf1_653[k]
                    + f_3 * pc_x[k] * sng_978[k];

        t_1369[k] = f_3 * pc_y[k] * sng_977[k];
    }

#pragma omp simd aligned(t_1370, t_1371, t_1372, t_1373, pc_x, pc_y, pc_z, smg_813, snf0_655, \
                         snf0_656, snf1_655, snf1_656, sng_978, sng_980, \
                         sng_981 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1370[k] = f_4 * snf0_655[k]
                    - f_5 * snf1_655[k]
                    + f_3 * pc_x[k] * sng_980[k];

        t_1371[k] = f_6 * snf0_656[k]
                    - f_7 * snf1_656[k]
                    + f_3 * pc_x[k] * sng_981[k];

        t_1372[k] = f_0 * smg_813[k]
                    + f_3 * pc_z[k] * sng_978[k];

        t_1373[k] = f_3 * pc_y[k] * sng_980[k];
    }

#pragma omp simd aligned(t_1374, t_1375, t_1376, t_1377, t_1378, t_1379, pc_x, snf0_659, \
                         snf1_659, sng_984, sng_985, sng_986, sng_987, sng_988, \
                         sng_989 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1374[k] = f_6 * snf0_659[k]
                    - f_7 * snf1_659[k]
                    + f_3 * pc_x[k] * sng_984[k];

        t_1375[k] = f_3 * pc_x[k] * sng_985[k];

        t_1376[k] = f_3 * pc_x[k] * sng_986[k];

        t_1377[k] = f_3 * pc_x[k] * sng_987[k];

        t_1378[k] = f_3 * pc_x[k] * sng_988[k];

        t_1379[k] = f_3 * pc_x[k] * sng_989[k];
    }

#pragma omp simd aligned(t_1380, t_1381, t_1382, t_1383, pc_y, pc_z, smg_820, snf0_656, \
                         snf0_658, snf0_659, snf1_656, snf1_658, snf1_659, sng_985, sng_987, \
                         sng_988 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1380[k] = f_1 * snf0_656[k]
                    - f_2 * snf1_656[k]
                    + f_3 * pc_y[k] * sng_985[k];

        t_1381[k] = f_0 * smg_820[k]
                    + f_3 * pc_z[k] * sng_985[k];

        t_1382[k] = f_4 * snf0_658[k]
                    - f_5 * snf1_658[k]
                    + f_3 * pc_y[k] * sng_987[k];

        t_1383[k] = f_6 * snf0_659[k]
                    - f_7 * snf1_659[k]
                    + f_3 * pc_y[k] * sng_988[k];
    }

#pragma omp simd aligned(t_1384, t_1385, pc_y, pc_z, smg_824, snf0_659, snf1_659, \
                         sng_989 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1384[k] = f_3 * pc_y[k] * sng_989[k];

        t_1385[k] = f_0 * smg_824[k]
                    + f_1 * snf0_659[k]
                    - f_2 * snf1_659[k]
                    + f_3 * pc_z[k] * sng_989[k];
    }
}

auto
compute_prim_snh_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t smh0, const size_t smg,
                                                   const size_t smh1, const size_t snf0,
                                                   const size_t snf1, const size_t sng,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_snh_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, smh0, smg,
                                                              smh1, snf0, snf1, sng, ncols,
                                                              gamma, p, q);

    compute_prim_snh_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, smh0, smg,
                                                              smh1, snf0, snf1, sng, ncols,
                                                              gamma, p, q);

    compute_prim_snh_three_center_electron_repulsion_0_piece2(buffer, target, pb, pc, smh0, smg,
                                                              smh1, snf0, snf1, sng, ncols,
                                                              gamma, p, q);

    compute_prim_snh_three_center_electron_repulsion_0_piece3(buffer, target, pb, pc, smh0, smg,
                                                              smh1, snf0, snf1, sng, ncols,
                                                              gamma, p, q);

    compute_prim_snh_three_center_electron_repulsion_0_piece4(buffer, target, pb, pc, smh0, smg,
                                                              smh1, snf0, snf1, sng, ncols,
                                                              gamma, p, q);

    compute_prim_snh_three_center_electron_repulsion_0_piece5(buffer, target, pb, pc, smh0, smg,
                                                              smh1, snf0, snf1, sng, ncols,
                                                              gamma, p, q);

    compute_prim_snh_three_center_electron_repulsion_0_piece6(buffer, target, pb, pc, smh0, smg,
                                                              smh1, snf0, snf1, sng, ncols,
                                                              gamma, p, q);

    compute_prim_snh_three_center_electron_repulsion_0_piece7(buffer, target, pb, pc, smh0, smg,
                                                              smh1, snf0, snf1, sng, ncols,
                                                              gamma, p, q);

    compute_prim_snh_three_center_electron_repulsion_0_piece8(buffer, target, pb, pc, smh0, smg,
                                                              smh1, snf0, snf1, sng, ncols,
                                                              gamma, p, q);

    compute_prim_snh_three_center_electron_repulsion_0_piece9(buffer, target, pb, pc, smh0, smg,
                                                              smh1, snf0, snf1, sng, ncols,
                                                              gamma, p, q);

    compute_prim_snh_three_center_electron_repulsion_0_piece10(buffer, target, pb, pc, smh0,
                                                               smg, smh1, snf0, snf1, sng,
                                                               ncols, gamma, p, q);

    compute_prim_snh_three_center_electron_repulsion_0_piece11(buffer, target, pb, pc, smh0,
                                                               smg, smh1, snf0, snf1, sng,
                                                               ncols, gamma, p, q);
}

}  // namespace simdt3ceri
