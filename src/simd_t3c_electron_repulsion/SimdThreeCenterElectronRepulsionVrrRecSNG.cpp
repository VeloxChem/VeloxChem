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


#include "SimdThreeCenterElectronRepulsionVrrRecSNG.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_sng_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t smg0,
                                                          const size_t smf, const size_t smg1,
                                                          const size_t snd0, const size_t snd1,
                                                          const size_t snf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.0 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 4.5 / q;
    const auto f_10 = 4.0 / q;
    const auto f_11 = 3.5 / q;
    const auto f_12 = 1.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *smg0_0 = buffer.data(smg0 + 0);
    const auto *smg0_3 = buffer.data(smg0 + 3);
    const auto *smg0_5 = buffer.data(smg0 + 5);
    const auto *smg0_10 = buffer.data(smg0 + 10);
    const auto *smg0_14 = buffer.data(smg0 + 14);
    const auto *smg0_18 = buffer.data(smg0 + 18);
    const auto *smg0_25 = buffer.data(smg0 + 25);
    const auto *smg0_30 = buffer.data(smg0 + 30);
    const auto *smg0_35 = buffer.data(smg0 + 35);
    const auto *smg0_44 = buffer.data(smg0 + 44);
    const auto *smg0_45 = buffer.data(smg0 + 45);
    const auto *smg0_48 = buffer.data(smg0 + 48);
    const auto *smg0_55 = buffer.data(smg0 + 55);
    const auto *smg0_75 = buffer.data(smg0 + 75);
    const auto *smg0_78 = buffer.data(smg0 + 78);

    const auto *smf_0 = buffer.data(smf + 0);
    const auto *smf_1 = buffer.data(smf + 1);
    const auto *smf_2 = buffer.data(smf + 2);
    const auto *smf_3 = buffer.data(smf + 3);
    const auto *smf_5 = buffer.data(smf + 5);
    const auto *smf_6 = buffer.data(smf + 6);
    const auto *smf_7 = buffer.data(smf + 7);
    const auto *smf_8 = buffer.data(smf + 8);
    const auto *smf_9 = buffer.data(smf + 9);
    const auto *smf_10 = buffer.data(smf + 10);
    const auto *smf_12 = buffer.data(smf + 12);
    const auto *smf_16 = buffer.data(smf + 16);
    const auto *smf_17 = buffer.data(smf + 17);
    const auto *smf_18 = buffer.data(smf + 18);
    const auto *smf_19 = buffer.data(smf + 19);
    const auto *smf_20 = buffer.data(smf + 20);
    const auto *smf_22 = buffer.data(smf + 22);
    const auto *smf_26 = buffer.data(smf + 26);
    const auto *smf_27 = buffer.data(smf + 27);
    const auto *smf_28 = buffer.data(smf + 28);
    const auto *smf_29 = buffer.data(smf + 29);
    const auto *smf_30 = buffer.data(smf + 30);
    const auto *smf_32 = buffer.data(smf + 32);
    const auto *smf_33 = buffer.data(smf + 33);
    const auto *smf_35 = buffer.data(smf + 35);
    const auto *smf_36 = buffer.data(smf + 36);
    const auto *smf_37 = buffer.data(smf + 37);
    const auto *smf_38 = buffer.data(smf + 38);
    const auto *smf_39 = buffer.data(smf + 39);
    const auto *smf_40 = buffer.data(smf + 40);
    const auto *smf_42 = buffer.data(smf + 42);
    const auto *smf_46 = buffer.data(smf + 46);
    const auto *smf_47 = buffer.data(smf + 47);
    const auto *smf_48 = buffer.data(smf + 48);
    const auto *smf_49 = buffer.data(smf + 49);
    const auto *smf_50 = buffer.data(smf + 50);
    const auto *smf_51 = buffer.data(smf + 51);
    const auto *smf_52 = buffer.data(smf + 52);
    const auto *smf_53 = buffer.data(smf + 53);
    const auto *smf_55 = buffer.data(smf + 55);
    const auto *smf_56 = buffer.data(smf + 56);
    const auto *smf_57 = buffer.data(smf + 57);
    const auto *smf_58 = buffer.data(smf + 58);
    const auto *smf_59 = buffer.data(smf + 59);
    const auto *smf_60 = buffer.data(smf + 60);
    const auto *smf_63 = buffer.data(smf + 63);
    const auto *smf_65 = buffer.data(smf + 65);
    const auto *smf_66 = buffer.data(smf + 66);
    const auto *smf_67 = buffer.data(smf + 67);
    const auto *smf_68 = buffer.data(smf + 68);
    const auto *smf_69 = buffer.data(smf + 69);
    const auto *smf_75 = buffer.data(smf + 75);
    const auto *smf_76 = buffer.data(smf + 76);
    const auto *smf_77 = buffer.data(smf + 77);
    const auto *smf_78 = buffer.data(smf + 78);
    const auto *smf_79 = buffer.data(smf + 79);

    const auto *smg1_0 = buffer.data(smg1 + 0);
    const auto *smg1_3 = buffer.data(smg1 + 3);
    const auto *smg1_5 = buffer.data(smg1 + 5);
    const auto *smg1_10 = buffer.data(smg1 + 10);
    const auto *smg1_14 = buffer.data(smg1 + 14);
    const auto *smg1_18 = buffer.data(smg1 + 18);
    const auto *smg1_25 = buffer.data(smg1 + 25);
    const auto *smg1_30 = buffer.data(smg1 + 30);
    const auto *smg1_35 = buffer.data(smg1 + 35);
    const auto *smg1_44 = buffer.data(smg1 + 44);
    const auto *smg1_45 = buffer.data(smg1 + 45);
    const auto *smg1_48 = buffer.data(smg1 + 48);
    const auto *smg1_55 = buffer.data(smg1 + 55);
    const auto *smg1_75 = buffer.data(smg1 + 75);
    const auto *smg1_78 = buffer.data(smg1 + 78);

    const auto *snd0_0 = buffer.data(snd0 + 0);
    const auto *snd0_3 = buffer.data(snd0 + 3);
    const auto *snd0_5 = buffer.data(snd0 + 5);
    const auto *snd0_9 = buffer.data(snd0 + 9);
    const auto *snd0_11 = buffer.data(snd0 + 11);
    const auto *snd0_17 = buffer.data(snd0 + 17);
    const auto *snd0_18 = buffer.data(snd0 + 18);
    const auto *snd0_21 = buffer.data(snd0 + 21);
    const auto *snd0_23 = buffer.data(snd0 + 23);
    const auto *snd0_29 = buffer.data(snd0 + 29);
    const auto *snd0_30 = buffer.data(snd0 + 30);
    const auto *snd0_33 = buffer.data(snd0 + 33);
    const auto *snd0_35 = buffer.data(snd0 + 35);
    const auto *snd0_36 = buffer.data(snd0 + 36);
    const auto *snd0_39 = buffer.data(snd0 + 39);
    const auto *snd0_41 = buffer.data(snd0 + 41);
    const auto *snd0_47 = buffer.data(snd0 + 47);

    const auto *snd1_0 = buffer.data(snd1 + 0);
    const auto *snd1_3 = buffer.data(snd1 + 3);
    const auto *snd1_5 = buffer.data(snd1 + 5);
    const auto *snd1_9 = buffer.data(snd1 + 9);
    const auto *snd1_11 = buffer.data(snd1 + 11);
    const auto *snd1_17 = buffer.data(snd1 + 17);
    const auto *snd1_18 = buffer.data(snd1 + 18);
    const auto *snd1_21 = buffer.data(snd1 + 21);
    const auto *snd1_23 = buffer.data(snd1 + 23);
    const auto *snd1_29 = buffer.data(snd1 + 29);
    const auto *snd1_30 = buffer.data(snd1 + 30);
    const auto *snd1_33 = buffer.data(snd1 + 33);
    const auto *snd1_35 = buffer.data(snd1 + 35);
    const auto *snd1_36 = buffer.data(snd1 + 36);
    const auto *snd1_39 = buffer.data(snd1 + 39);
    const auto *snd1_41 = buffer.data(snd1 + 41);
    const auto *snd1_47 = buffer.data(snd1 + 47);

    const auto *snf_0 = buffer.data(snf + 0);
    const auto *snf_2 = buffer.data(snf + 2);
    const auto *snf_3 = buffer.data(snf + 3);
    const auto *snf_5 = buffer.data(snf + 5);
    const auto *snf_6 = buffer.data(snf + 6);
    const auto *snf_7 = buffer.data(snf + 7);
    const auto *snf_8 = buffer.data(snf + 8);
    const auto *snf_9 = buffer.data(snf + 9);
    const auto *snf_10 = buffer.data(snf + 10);
    const auto *snf_12 = buffer.data(snf + 12);
    const auto *snf_16 = buffer.data(snf + 16);
    const auto *snf_17 = buffer.data(snf + 17);
    const auto *snf_18 = buffer.data(snf + 18);
    const auto *snf_19 = buffer.data(snf + 19);
    const auto *snf_20 = buffer.data(snf + 20);
    const auto *snf_22 = buffer.data(snf + 22);
    const auto *snf_26 = buffer.data(snf + 26);
    const auto *snf_27 = buffer.data(snf + 27);
    const auto *snf_28 = buffer.data(snf + 28);
    const auto *snf_29 = buffer.data(snf + 29);
    const auto *snf_30 = buffer.data(snf + 30);
    const auto *snf_32 = buffer.data(snf + 32);
    const auto *snf_33 = buffer.data(snf + 33);
    const auto *snf_35 = buffer.data(snf + 35);
    const auto *snf_36 = buffer.data(snf + 36);
    const auto *snf_37 = buffer.data(snf + 37);
    const auto *snf_38 = buffer.data(snf + 38);
    const auto *snf_39 = buffer.data(snf + 39);
    const auto *snf_40 = buffer.data(snf + 40);
    const auto *snf_42 = buffer.data(snf + 42);
    const auto *snf_46 = buffer.data(snf + 46);
    const auto *snf_47 = buffer.data(snf + 47);
    const auto *snf_48 = buffer.data(snf + 48);
    const auto *snf_49 = buffer.data(snf + 49);
    const auto *snf_50 = buffer.data(snf + 50);
    const auto *snf_52 = buffer.data(snf + 52);
    const auto *snf_53 = buffer.data(snf + 53);
    const auto *snf_55 = buffer.data(snf + 55);
    const auto *snf_56 = buffer.data(snf + 56);
    const auto *snf_57 = buffer.data(snf + 57);
    const auto *snf_58 = buffer.data(snf + 58);
    const auto *snf_59 = buffer.data(snf + 59);
    const auto *snf_60 = buffer.data(snf + 60);
    const auto *snf_62 = buffer.data(snf + 62);
    const auto *snf_63 = buffer.data(snf + 63);
    const auto *snf_65 = buffer.data(snf + 65);
    const auto *snf_66 = buffer.data(snf + 66);
    const auto *snf_67 = buffer.data(snf + 67);
    const auto *snf_68 = buffer.data(snf + 68);
    const auto *snf_69 = buffer.data(snf + 69);
    const auto *snf_70 = buffer.data(snf + 70);
    const auto *snf_72 = buffer.data(snf + 72);
    const auto *snf_75 = buffer.data(snf + 75);
    const auto *snf_76 = buffer.data(snf + 76);
    const auto *snf_77 = buffer.data(snf + 77);
    const auto *snf_78 = buffer.data(snf + 78);
    const auto *snf_79 = buffer.data(snf + 79);
    const auto *snf_80 = buffer.data(snf + 80);
    const auto *snf_82 = buffer.data(snf + 82);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, smf_0, smf_3, snd0_0, snd0_3, \
                         snd1_0, snd1_3, snf_0, snf_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * smf_0[k]
                 + f_1 * snd0_0[k]
                 - f_2 * snd1_0[k]
                 + f_3 * pc_x[k] * snf_0[k];

        t_1[k] = f_3 * pc_y[k] * snf_0[k];

        t_2[k] = f_3 * pc_z[k] * snf_0[k];

        t_3[k] = f_0 * smf_3[k]
                 + f_4 * snd0_3[k]
                 - f_5 * snd1_3[k]
                 + f_3 * pc_x[k] * snf_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pc_x, pc_y, smf_5, smf_6, smf_7, snd0_5, snd1_5, \
                         snf_2, snf_5, snf_6, snf_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * pc_y[k] * snf_2[k];

        t_5[k] = f_0 * smf_5[k]
                 + f_4 * snd0_5[k]
                 - f_5 * snd1_5[k]
                 + f_3 * pc_x[k] * snf_5[k];

        t_6[k] = f_0 * smf_6[k]
                 + f_3 * pc_x[k] * snf_6[k];

        t_7[k] = f_0 * smf_7[k]
                 + f_3 * pc_x[k] * snf_7[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pc_x, pc_y, pc_z, smf_8, smf_9, snd0_3, snd1_3, \
                         snf_6, snf_8, snf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * smf_8[k]
                 + f_3 * pc_x[k] * snf_8[k];

        t_9[k] = f_0 * smf_9[k]
                 + f_3 * pc_x[k] * snf_9[k];

        t_10[k] = f_1 * snd0_3[k]
                  - f_2 * snd1_3[k]
                  + f_3 * pc_y[k] * snf_6[k];

        t_11[k] = f_3 * pc_z[k] * snf_6[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, pb_y, pc_y, pc_z, smg0_0, smf_0, \
                         smg1_0, snd0_5, snd1_5, snf_8, snf_9, snf_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_4 * snd0_5[k]
                  - f_5 * snd1_5[k]
                  + f_3 * pc_y[k] * snf_8[k];

        t_13[k] = f_3 * pc_y[k] * snf_9[k];

        t_14[k] = f_1 * snd0_5[k]
                  - f_2 * snd1_5[k]
                  + f_3 * pc_z[k] * snf_9[k];

        t_15[k] = pb_y[k] * smg0_0[k]
                  - f_6 * pc_y[k] * smg1_0[k];

        t_16[k] = f_7 * smf_0[k]
                  + f_3 * pc_y[k] * snf_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pb_y, pc_y, pc_z, smg0_3, smg0_5, smf_1, \
                         smf_2, smg1_3, smg1_5, snf_10, snf_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_3 * pc_z[k] * snf_10[k];

        t_18[k] = pb_y[k] * smg0_3[k]
                  + f_8 * smf_1[k]
                  - f_6 * pc_y[k] * smg1_3[k];

        t_19[k] = f_7 * smf_2[k]
                  + f_3 * pc_y[k] * snf_12[k];

        t_20[k] = pb_y[k] * smg0_5[k]
                  - f_6 * pc_y[k] * smg1_5[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pc_x, smf_16, smf_17, smf_18, smf_19, snf_16, \
                         snf_17, snf_18, snf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_9 * smf_16[k]
                  + f_3 * pc_x[k] * snf_16[k];

        t_22[k] = f_9 * smf_17[k]
                  + f_3 * pc_x[k] * snf_17[k];

        t_23[k] = f_9 * smf_18[k]
                  + f_3 * pc_x[k] * snf_18[k];

        t_24[k] = f_9 * smf_19[k]
                  + f_3 * pc_x[k] * snf_19[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pc_y, pc_z, smf_6, smf_8, smf_9, snd0_9, \
                         snd0_11, snd1_9, snd1_11, snf_16, snf_18, \
                         snf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_7 * smf_6[k]
                  + f_1 * snd0_9[k]
                  - f_2 * snd1_9[k]
                  + f_3 * pc_y[k] * snf_16[k];

        t_26[k] = f_3 * pc_z[k] * snf_16[k];

        t_27[k] = f_7 * smf_8[k]
                  + f_4 * snd0_11[k]
                  - f_5 * snd1_11[k]
                  + f_3 * pc_y[k] * snf_18[k];

        t_28[k] = f_7 * smf_9[k]
                  + f_3 * pc_y[k] * snf_19[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pb_y, pb_z, pc_y, pc_z, smg0_0, smg0_14, \
                         smf_0, smg1_0, smg1_14, snf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pb_y[k] * smg0_14[k]
                  - f_6 * pc_y[k] * smg1_14[k];

        t_30[k] = pb_z[k] * smg0_0[k]
                  - f_6 * pc_z[k] * smg1_0[k];

        t_31[k] = f_3 * pc_y[k] * snf_20[k];

        t_32[k] = f_7 * smf_0[k]
                  + f_3 * pc_z[k] * snf_20[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pb_z, pc_x, pc_y, pc_z, smg0_3, smg0_5, \
                         smf_2, smf_26, smg1_3, smg1_5, snf_22, \
                         snf_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pb_z[k] * smg0_3[k]
                  - f_6 * pc_z[k] * smg1_3[k];

        t_34[k] = f_3 * pc_y[k] * snf_22[k];

        t_35[k] = pb_z[k] * smg0_5[k]
                  + f_8 * smf_2[k]
                  - f_6 * pc_z[k] * smg1_5[k];

        t_36[k] = f_9 * smf_26[k]
                  + f_3 * pc_x[k] * snf_26[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pb_z, pc_x, pc_z, smg0_10, smf_27, smf_28, \
                         smf_29, smg1_10, snf_27, snf_28, snf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_9 * smf_27[k]
                  + f_3 * pc_x[k] * snf_27[k];

        t_38[k] = f_9 * smf_28[k]
                  + f_3 * pc_x[k] * snf_28[k];

        t_39[k] = f_9 * smf_29[k]
                  + f_3 * pc_x[k] * snf_29[k];

        t_40[k] = pb_z[k] * smg0_10[k]
                  - f_6 * pc_z[k] * smg1_10[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pc_y, pc_z, smf_6, smf_9, snd0_17, snd1_17, \
                         snf_26, snf_28, snf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_7 * smf_6[k]
                  + f_3 * pc_z[k] * snf_26[k];

        t_42[k] = f_4 * snd0_17[k]
                  - f_5 * snd1_17[k]
                  + f_3 * pc_y[k] * snf_28[k];

        t_43[k] = f_3 * pc_y[k] * snf_29[k];

        t_44[k] = f_7 * smf_9[k]
                  + f_1 * snd0_17[k]
                  - f_2 * snd1_17[k]
                  + f_3 * pc_z[k] * snf_29[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pc_x, pc_y, pc_z, smf_10, smf_30, smf_33, \
                         snd0_18, snd0_21, snd1_18, snd1_21, snf_30, \
                         snf_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_10 * smf_30[k]
                  + f_1 * snd0_18[k]
                  - f_2 * snd1_18[k]
                  + f_3 * pc_x[k] * snf_30[k];

        t_46[k] = f_8 * smf_10[k]
                  + f_3 * pc_y[k] * snf_30[k];

        t_47[k] = f_3 * pc_z[k] * snf_30[k];

        t_48[k] = f_10 * smf_33[k]
                  + f_4 * snd0_21[k]
                  - f_5 * snd1_21[k]
                  + f_3 * pc_x[k] * snf_33[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pc_x, pc_y, smf_12, smf_35, smf_36, smf_37, \
                         snd0_23, snd1_23, snf_32, snf_35, snf_36, \
                         snf_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_8 * smf_12[k]
                  + f_3 * pc_y[k] * snf_32[k];

        t_50[k] = f_10 * smf_35[k]
                  + f_4 * snd0_23[k]
                  - f_5 * snd1_23[k]
                  + f_3 * pc_x[k] * snf_35[k];

        t_51[k] = f_10 * smf_36[k]
                  + f_3 * pc_x[k] * snf_36[k];

        t_52[k] = f_10 * smf_37[k]
                  + f_3 * pc_x[k] * snf_37[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, pc_x, pc_y, pc_z, smf_16, smf_38, smf_39, \
                         snd0_21, snd1_21, snf_36, snf_38, snf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_10 * smf_38[k]
                  + f_3 * pc_x[k] * snf_38[k];

        t_54[k] = f_10 * smf_39[k]
                  + f_3 * pc_x[k] * snf_39[k];

        t_55[k] = f_8 * smf_16[k]
                  + f_1 * snd0_21[k]
                  - f_2 * snd1_21[k]
                  + f_3 * pc_y[k] * snf_36[k];

        t_56[k] = f_3 * pc_z[k] * snf_36[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, pb_y, pc_y, pc_z, smg0_30, smf_18, smf_19, \
                         smg1_30, snd0_23, snd1_23, snf_38, snf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_8 * smf_18[k]
                  + f_4 * snd0_23[k]
                  - f_5 * snd1_23[k]
                  + f_3 * pc_y[k] * snf_38[k];

        t_58[k] = f_8 * smf_19[k]
                  + f_3 * pc_y[k] * snf_39[k];

        t_59[k] = f_1 * snd0_23[k]
                  - f_2 * snd1_23[k]
                  + f_3 * pc_z[k] * snf_39[k];

        t_60[k] = pb_y[k] * smg0_30[k]
                  - f_6 * pc_y[k] * smg1_30[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pb_z, pc_y, pc_z, smg0_18, smf_10, smf_20, \
                         smf_22, smg1_18, snf_40, snf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_7 * smf_20[k]
                  + f_3 * pc_y[k] * snf_40[k];

        t_62[k] = f_7 * smf_10[k]
                  + f_3 * pc_z[k] * snf_40[k];

        t_63[k] = pb_z[k] * smg0_18[k]
                  - f_6 * pc_z[k] * smg1_18[k];

        t_64[k] = f_7 * smf_22[k]
                  + f_3 * pc_y[k] * snf_42[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pb_y, pc_x, pc_y, smg0_35, smf_46, smf_47, \
                         smf_48, smg1_35, snf_46, snf_47, snf_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = pb_y[k] * smg0_35[k]
                  - f_6 * pc_y[k] * smg1_35[k];

        t_66[k] = f_10 * smf_46[k]
                  + f_3 * pc_x[k] * snf_46[k];

        t_67[k] = f_10 * smf_47[k]
                  + f_3 * pc_x[k] * snf_47[k];

        t_68[k] = f_10 * smf_48[k]
                  + f_3 * pc_x[k] * snf_48[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pb_z, pc_x, pc_z, smg0_25, smf_16, smf_49, smg1_25, \
                         snf_46, snf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_10 * smf_49[k]
                  + f_3 * pc_x[k] * snf_49[k];

        t_70[k] = pb_z[k] * smg0_25[k]
                  - f_6 * pc_z[k] * smg1_25[k];

        t_71[k] = f_7 * smf_16[k]
                  + f_3 * pc_z[k] * snf_46[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, pb_y, pc_y, smg0_44, smf_28, smf_29, smg1_44, \
                         snd0_29, snd1_29, snf_48, snf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_7 * smf_28[k]
                  + f_4 * snd0_29[k]
                  - f_5 * snd1_29[k]
                  + f_3 * pc_y[k] * snf_48[k];

        t_73[k] = f_7 * smf_29[k]
                  + f_3 * pc_y[k] * snf_49[k];

        t_74[k] = pb_y[k] * smg0_44[k]
                  - f_6 * pc_y[k] * smg1_44[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pc_x, pc_y, pc_z, smf_20, smf_50, smf_53, \
                         snd0_30, snd0_33, snd1_30, snd1_33, snf_50, \
                         snf_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_10 * smf_50[k]
                  + f_1 * snd0_30[k]
                  - f_2 * snd1_30[k]
                  + f_3 * pc_x[k] * snf_50[k];

        t_76[k] = f_3 * pc_y[k] * snf_50[k];

        t_77[k] = f_8 * smf_20[k]
                  + f_3 * pc_z[k] * snf_50[k];

        t_78[k] = f_10 * smf_53[k]
                  + f_4 * snd0_33[k]
                  - f_5 * snd1_33[k]
                  + f_3 * pc_x[k] * snf_53[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pc_x, pc_y, smf_55, smf_56, smf_57, snd0_35, \
                         snd1_35, snf_52, snf_55, snf_56, snf_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_3 * pc_y[k] * snf_52[k];

        t_80[k] = f_10 * smf_55[k]
                  + f_4 * snd0_35[k]
                  - f_5 * snd1_35[k]
                  + f_3 * pc_x[k] * snf_55[k];

        t_81[k] = f_10 * smf_56[k]
                  + f_3 * pc_x[k] * snf_56[k];

        t_82[k] = f_10 * smf_57[k]
                  + f_3 * pc_x[k] * snf_57[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pc_x, pc_y, pc_z, smf_26, smf_58, smf_59, \
                         snd0_33, snd1_33, snf_56, snf_58, snf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_10 * smf_58[k]
                  + f_3 * pc_x[k] * snf_58[k];

        t_84[k] = f_10 * smf_59[k]
                  + f_3 * pc_x[k] * snf_59[k];

        t_85[k] = f_1 * snd0_33[k]
                  - f_2 * snd1_33[k]
                  + f_3 * pc_y[k] * snf_56[k];

        t_86[k] = f_8 * smf_26[k]
                  + f_3 * pc_z[k] * snf_56[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pc_x, pc_y, pc_z, smf_29, smf_60, snd0_35, \
                         snd0_36, snd1_35, snd1_36, snf_58, snf_59, \
                         snf_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_4 * snd0_35[k]
                  - f_5 * snd1_35[k]
                  + f_3 * pc_y[k] * snf_58[k];

        t_88[k] = f_3 * pc_y[k] * snf_59[k];

        t_89[k] = f_8 * smf_29[k]
                  + f_1 * snd0_35[k]
                  - f_2 * snd1_35[k]
                  + f_3 * pc_z[k] * snf_59[k];

        t_90[k] = f_11 * smf_60[k]
                  + f_1 * snd0_36[k]
                  - f_2 * snd1_36[k]
                  + f_3 * pc_x[k] * snf_60[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, pc_x, pc_y, pc_z, smf_30, smf_32, smf_63, \
                         snd0_39, snd1_39, snf_60, snf_62, snf_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_12 * smf_30[k]
                  + f_3 * pc_y[k] * snf_60[k];

        t_92[k] = f_3 * pc_z[k] * snf_60[k];

        t_93[k] = f_11 * smf_63[k]
                  + f_4 * snd0_39[k]
                  - f_5 * snd1_39[k]
                  + f_3 * pc_x[k] * snf_63[k];

        t_94[k] = f_12 * smf_32[k]
                  + f_3 * pc_y[k] * snf_62[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pc_x, smf_65, smf_66, smf_67, smf_68, \
                         snd0_41, snd1_41, snf_65, snf_66, snf_67, \
                         snf_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_11 * smf_65[k]
                  + f_4 * snd0_41[k]
                  - f_5 * snd1_41[k]
                  + f_3 * pc_x[k] * snf_65[k];

        t_96[k] = f_11 * smf_66[k]
                  + f_3 * pc_x[k] * snf_66[k];

        t_97[k] = f_11 * smf_67[k]
                  + f_3 * pc_x[k] * snf_67[k];

        t_98[k] = f_11 * smf_68[k]
                  + f_3 * pc_x[k] * snf_68[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, pc_x, pc_y, pc_z, smf_36, smf_69, snd0_39, \
                         snd1_39, snf_66, snf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_11 * smf_69[k]
                  + f_3 * pc_x[k] * snf_69[k];

        t_100[k] = f_12 * smf_36[k]
                   + f_1 * snd0_39[k]
                   - f_2 * snd1_39[k]
                   + f_3 * pc_y[k] * snf_66[k];

        t_101[k] = f_3 * pc_z[k] * snf_66[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, pb_z, pc_y, pc_z, smg0_45, smf_38, \
                         smf_39, smg1_45, snd0_41, snd1_41, snf_68, \
                         snf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_12 * smf_38[k]
                   + f_4 * snd0_41[k]
                   - f_5 * snd1_41[k]
                   + f_3 * pc_y[k] * snf_68[k];

        t_103[k] = f_12 * smf_39[k]
                   + f_3 * pc_y[k] * snf_69[k];

        t_104[k] = f_1 * snd0_41[k]
                   - f_2 * snd1_41[k]
                   + f_3 * pc_z[k] * snf_69[k];

        t_105[k] = pb_z[k] * smg0_45[k]
                   - f_6 * pc_z[k] * smg1_45[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, pb_z, pc_y, pc_z, smg0_48, smf_30, \
                         smf_40, smf_42, smg1_48, snf_70, snf_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_8 * smf_40[k]
                   + f_3 * pc_y[k] * snf_70[k];

        t_107[k] = f_7 * smf_30[k]
                   + f_3 * pc_z[k] * snf_70[k];

        t_108[k] = pb_z[k] * smg0_48[k]
                   - f_6 * pc_z[k] * smg1_48[k];

        t_109[k] = f_8 * smf_42[k]
                   + f_3 * pc_y[k] * snf_72[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, pc_x, smf_75, smf_76, smf_77, smf_78, \
                         snd0_47, snd1_47, snf_75, snf_76, snf_77, \
                         snf_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = f_11 * smf_75[k]
                   + f_4 * snd0_47[k]
                   - f_5 * snd1_47[k]
                   + f_3 * pc_x[k] * snf_75[k];

        t_111[k] = f_11 * smf_76[k]
                   + f_3 * pc_x[k] * snf_76[k];

        t_112[k] = f_11 * smf_77[k]
                   + f_3 * pc_x[k] * snf_77[k];

        t_113[k] = f_11 * smf_78[k]
                   + f_3 * pc_x[k] * snf_78[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, pb_z, pc_x, pc_z, smg0_55, smf_36, smf_79, \
                         smg1_55, snf_76, snf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_11 * smf_79[k]
                   + f_3 * pc_x[k] * snf_79[k];

        t_115[k] = pb_z[k] * smg0_55[k]
                   - f_6 * pc_z[k] * smg1_55[k];

        t_116[k] = f_7 * smf_36[k]
                   + f_3 * pc_z[k] * snf_76[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pb_y, pc_y, pc_z, smg0_75, smf_39, \
                         smf_48, smf_49, smg1_75, snd0_47, snd1_47, snf_78, \
                         snf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_8 * smf_48[k]
                   + f_4 * snd0_47[k]
                   - f_5 * snd1_47[k]
                   + f_3 * pc_y[k] * snf_78[k];

        t_118[k] = f_8 * smf_49[k]
                   + f_3 * pc_y[k] * snf_79[k];

        t_119[k] = f_7 * smf_39[k]
                   + f_1 * snd0_47[k]
                   - f_2 * snd1_47[k]
                   + f_3 * pc_z[k] * snf_79[k];

        t_120[k] = pb_y[k] * smg0_75[k]
                   - f_6 * pc_y[k] * smg1_75[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pb_y, pc_y, pc_z, smg0_78, smf_40, \
                         smf_50, smf_51, smf_52, smg1_78, snf_80, \
                         snf_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_7 * smf_50[k]
                   + f_3 * pc_y[k] * snf_80[k];

        t_122[k] = f_8 * smf_40[k]
                   + f_3 * pc_z[k] * snf_80[k];

        t_123[k] = pb_y[k] * smg0_78[k]
                   + f_8 * smf_51[k]
                   - f_6 * pc_y[k] * smg1_78[k];

        t_124[k] = f_7 * smf_52[k]
                   + f_3 * pc_y[k] * snf_82[k];
    }
}

static auto
compute_prim_sng_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t smg0,
                                                          const size_t smf, const size_t smg1,
                                                          const size_t snd0, const size_t snd1,
                                                          const size_t snf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_11 = 3.5 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 3.0 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 2.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *smg0_80 = buffer.data(smg0 + 80);
    const auto *smg0_89 = buffer.data(smg0 + 89);
    const auto *smg0_90 = buffer.data(smg0 + 90);
    const auto *smg0_93 = buffer.data(smg0 + 93);
    const auto *smg0_100 = buffer.data(smg0 + 100);
    const auto *smg0_135 = buffer.data(smg0 + 135);
    const auto *smg0_138 = buffer.data(smg0 + 138);
    const auto *smg0_140 = buffer.data(smg0 + 140);
    const auto *smg0_149 = buffer.data(smg0 + 149);
    const auto *smg0_150 = buffer.data(smg0 + 150);
    const auto *smg0_153 = buffer.data(smg0 + 153);

    const auto *smf_46 = buffer.data(smf + 46);
    const auto *smf_50 = buffer.data(smf + 50);
    const auto *smf_56 = buffer.data(smf + 56);
    const auto *smf_58 = buffer.data(smf + 58);
    const auto *smf_59 = buffer.data(smf + 59);
    const auto *smf_60 = buffer.data(smf + 60);
    const auto *smf_62 = buffer.data(smf + 62);
    const auto *smf_66 = buffer.data(smf + 66);
    const auto *smf_68 = buffer.data(smf + 68);
    const auto *smf_69 = buffer.data(smf + 69);
    const auto *smf_70 = buffer.data(smf + 70);
    const auto *smf_72 = buffer.data(smf + 72);
    const auto *smf_76 = buffer.data(smf + 76);
    const auto *smf_78 = buffer.data(smf + 78);
    const auto *smf_79 = buffer.data(smf + 79);
    const auto *smf_80 = buffer.data(smf + 80);
    const auto *smf_82 = buffer.data(smf + 82);
    const auto *smf_86 = buffer.data(smf + 86);
    const auto *smf_87 = buffer.data(smf + 87);
    const auto *smf_88 = buffer.data(smf + 88);
    const auto *smf_89 = buffer.data(smf + 89);
    const auto *smf_90 = buffer.data(smf + 90);
    const auto *smf_91 = buffer.data(smf + 91);
    const auto *smf_92 = buffer.data(smf + 92);
    const auto *smf_93 = buffer.data(smf + 93);
    const auto *smf_95 = buffer.data(smf + 95);
    const auto *smf_96 = buffer.data(smf + 96);
    const auto *smf_97 = buffer.data(smf + 97);
    const auto *smf_98 = buffer.data(smf + 98);
    const auto *smf_99 = buffer.data(smf + 99);
    const auto *smf_100 = buffer.data(smf + 100);
    const auto *smf_102 = buffer.data(smf + 102);
    const auto *smf_103 = buffer.data(smf + 103);
    const auto *smf_105 = buffer.data(smf + 105);
    const auto *smf_106 = buffer.data(smf + 106);
    const auto *smf_107 = buffer.data(smf + 107);
    const auto *smf_108 = buffer.data(smf + 108);
    const auto *smf_109 = buffer.data(smf + 109);
    const auto *smf_110 = buffer.data(smf + 110);
    const auto *smf_112 = buffer.data(smf + 112);
    const auto *smf_115 = buffer.data(smf + 115);
    const auto *smf_116 = buffer.data(smf + 116);
    const auto *smf_117 = buffer.data(smf + 117);
    const auto *smf_118 = buffer.data(smf + 118);
    const auto *smf_119 = buffer.data(smf + 119);
    const auto *smf_120 = buffer.data(smf + 120);
    const auto *smf_123 = buffer.data(smf + 123);
    const auto *smf_125 = buffer.data(smf + 125);
    const auto *smf_126 = buffer.data(smf + 126);
    const auto *smf_127 = buffer.data(smf + 127);
    const auto *smf_128 = buffer.data(smf + 128);
    const auto *smf_129 = buffer.data(smf + 129);
    const auto *smf_136 = buffer.data(smf + 136);
    const auto *smf_137 = buffer.data(smf + 137);
    const auto *smf_138 = buffer.data(smf + 138);
    const auto *smf_139 = buffer.data(smf + 139);
    const auto *smf_140 = buffer.data(smf + 140);
    const auto *smf_143 = buffer.data(smf + 143);
    const auto *smf_145 = buffer.data(smf + 145);
    const auto *smf_146 = buffer.data(smf + 146);
    const auto *smf_147 = buffer.data(smf + 147);
    const auto *smf_148 = buffer.data(smf + 148);
    const auto *smf_149 = buffer.data(smf + 149);
    const auto *smf_150 = buffer.data(smf + 150);
    const auto *smf_153 = buffer.data(smf + 153);
    const auto *smf_155 = buffer.data(smf + 155);
    const auto *smf_156 = buffer.data(smf + 156);
    const auto *smf_157 = buffer.data(smf + 157);
    const auto *smf_158 = buffer.data(smf + 158);
    const auto *smf_159 = buffer.data(smf + 159);

    const auto *smg1_80 = buffer.data(smg1 + 80);
    const auto *smg1_89 = buffer.data(smg1 + 89);
    const auto *smg1_90 = buffer.data(smg1 + 90);
    const auto *smg1_93 = buffer.data(smg1 + 93);
    const auto *smg1_100 = buffer.data(smg1 + 100);
    const auto *smg1_135 = buffer.data(smg1 + 135);
    const auto *smg1_138 = buffer.data(smg1 + 138);
    const auto *smg1_140 = buffer.data(smg1 + 140);
    const auto *smg1_149 = buffer.data(smg1 + 149);
    const auto *smg1_150 = buffer.data(smg1 + 150);
    const auto *smg1_153 = buffer.data(smg1 + 153);

    const auto *snd0_51 = buffer.data(snd0 + 51);
    const auto *snd0_53 = buffer.data(snd0 + 53);
    const auto *snd0_54 = buffer.data(snd0 + 54);
    const auto *snd0_57 = buffer.data(snd0 + 57);
    const auto *snd0_59 = buffer.data(snd0 + 59);
    const auto *snd0_60 = buffer.data(snd0 + 60);
    const auto *snd0_63 = buffer.data(snd0 + 63);
    const auto *snd0_65 = buffer.data(snd0 + 65);
    const auto *snd0_71 = buffer.data(snd0 + 71);
    const auto *snd0_72 = buffer.data(snd0 + 72);
    const auto *snd0_75 = buffer.data(snd0 + 75);
    const auto *snd0_77 = buffer.data(snd0 + 77);
    const auto *snd0_81 = buffer.data(snd0 + 81);
    const auto *snd0_83 = buffer.data(snd0 + 83);
    const auto *snd0_84 = buffer.data(snd0 + 84);
    const auto *snd0_87 = buffer.data(snd0 + 87);
    const auto *snd0_89 = buffer.data(snd0 + 89);
    const auto *snd0_90 = buffer.data(snd0 + 90);
    const auto *snd0_93 = buffer.data(snd0 + 93);
    const auto *snd0_95 = buffer.data(snd0 + 95);

    const auto *snd1_51 = buffer.data(snd1 + 51);
    const auto *snd1_53 = buffer.data(snd1 + 53);
    const auto *snd1_54 = buffer.data(snd1 + 54);
    const auto *snd1_57 = buffer.data(snd1 + 57);
    const auto *snd1_59 = buffer.data(snd1 + 59);
    const auto *snd1_60 = buffer.data(snd1 + 60);
    const auto *snd1_63 = buffer.data(snd1 + 63);
    const auto *snd1_65 = buffer.data(snd1 + 65);
    const auto *snd1_71 = buffer.data(snd1 + 71);
    const auto *snd1_72 = buffer.data(snd1 + 72);
    const auto *snd1_75 = buffer.data(snd1 + 75);
    const auto *snd1_77 = buffer.data(snd1 + 77);
    const auto *snd1_81 = buffer.data(snd1 + 81);
    const auto *snd1_83 = buffer.data(snd1 + 83);
    const auto *snd1_84 = buffer.data(snd1 + 84);
    const auto *snd1_87 = buffer.data(snd1 + 87);
    const auto *snd1_89 = buffer.data(snd1 + 89);
    const auto *snd1_90 = buffer.data(snd1 + 90);
    const auto *snd1_93 = buffer.data(snd1 + 93);
    const auto *snd1_95 = buffer.data(snd1 + 95);

    const auto *snf_86 = buffer.data(snf + 86);
    const auto *snf_87 = buffer.data(snf + 87);
    const auto *snf_88 = buffer.data(snf + 88);
    const auto *snf_89 = buffer.data(snf + 89);
    const auto *snf_90 = buffer.data(snf + 90);
    const auto *snf_92 = buffer.data(snf + 92);
    const auto *snf_93 = buffer.data(snf + 93);
    const auto *snf_95 = buffer.data(snf + 95);
    const auto *snf_96 = buffer.data(snf + 96);
    const auto *snf_97 = buffer.data(snf + 97);
    const auto *snf_98 = buffer.data(snf + 98);
    const auto *snf_99 = buffer.data(snf + 99);
    const auto *snf_100 = buffer.data(snf + 100);
    const auto *snf_102 = buffer.data(snf + 102);
    const auto *snf_103 = buffer.data(snf + 103);
    const auto *snf_105 = buffer.data(snf + 105);
    const auto *snf_106 = buffer.data(snf + 106);
    const auto *snf_107 = buffer.data(snf + 107);
    const auto *snf_108 = buffer.data(snf + 108);
    const auto *snf_109 = buffer.data(snf + 109);
    const auto *snf_110 = buffer.data(snf + 110);
    const auto *snf_112 = buffer.data(snf + 112);
    const auto *snf_115 = buffer.data(snf + 115);
    const auto *snf_116 = buffer.data(snf + 116);
    const auto *snf_117 = buffer.data(snf + 117);
    const auto *snf_118 = buffer.data(snf + 118);
    const auto *snf_119 = buffer.data(snf + 119);
    const auto *snf_120 = buffer.data(snf + 120);
    const auto *snf_122 = buffer.data(snf + 122);
    const auto *snf_123 = buffer.data(snf + 123);
    const auto *snf_125 = buffer.data(snf + 125);
    const auto *snf_126 = buffer.data(snf + 126);
    const auto *snf_127 = buffer.data(snf + 127);
    const auto *snf_128 = buffer.data(snf + 128);
    const auto *snf_129 = buffer.data(snf + 129);
    const auto *snf_130 = buffer.data(snf + 130);
    const auto *snf_132 = buffer.data(snf + 132);
    const auto *snf_136 = buffer.data(snf + 136);
    const auto *snf_137 = buffer.data(snf + 137);
    const auto *snf_138 = buffer.data(snf + 138);
    const auto *snf_139 = buffer.data(snf + 139);
    const auto *snf_140 = buffer.data(snf + 140);
    const auto *snf_142 = buffer.data(snf + 142);
    const auto *snf_143 = buffer.data(snf + 143);
    const auto *snf_145 = buffer.data(snf + 145);
    const auto *snf_146 = buffer.data(snf + 146);
    const auto *snf_147 = buffer.data(snf + 147);
    const auto *snf_148 = buffer.data(snf + 148);
    const auto *snf_149 = buffer.data(snf + 149);
    const auto *snf_150 = buffer.data(snf + 150);
    const auto *snf_152 = buffer.data(snf + 152);
    const auto *snf_153 = buffer.data(snf + 153);
    const auto *snf_155 = buffer.data(snf + 155);
    const auto *snf_156 = buffer.data(snf + 156);
    const auto *snf_157 = buffer.data(snf + 157);
    const auto *snf_158 = buffer.data(snf + 158);
    const auto *snf_159 = buffer.data(snf + 159);
    const auto *snf_160 = buffer.data(snf + 160);
    const auto *snf_162 = buffer.data(snf + 162);

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pb_y, pc_x, pc_y, smg0_80, smf_86, \
                         smf_87, smf_88, smg1_80, snf_86, snf_87, \
                         snf_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = pb_y[k] * smg0_80[k]
                   - f_6 * pc_y[k] * smg1_80[k];

        t_126[k] = f_11 * smf_86[k]
                   + f_3 * pc_x[k] * snf_86[k];

        t_127[k] = f_11 * smf_87[k]
                   + f_3 * pc_x[k] * snf_87[k];

        t_128[k] = f_11 * smf_88[k]
                   + f_3 * pc_x[k] * snf_88[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, pc_x, pc_y, pc_z, smf_46, smf_56, smf_89, \
                         snd0_51, snd1_51, snf_86, snf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_11 * smf_89[k]
                   + f_3 * pc_x[k] * snf_89[k];

        t_130[k] = f_7 * smf_56[k]
                   + f_1 * snd0_51[k]
                   - f_2 * snd1_51[k]
                   + f_3 * pc_y[k] * snf_86[k];

        t_131[k] = f_8 * smf_46[k]
                   + f_3 * pc_z[k] * snf_86[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, pb_y, pc_y, smg0_89, smf_58, smf_59, smg1_89, \
                         snd0_53, snd1_53, snf_88, snf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_7 * smf_58[k]
                   + f_4 * snd0_53[k]
                   - f_5 * snd1_53[k]
                   + f_3 * pc_y[k] * snf_88[k];

        t_133[k] = f_7 * smf_59[k]
                   + f_3 * pc_y[k] * snf_89[k];

        t_134[k] = pb_y[k] * smg0_89[k]
                   - f_6 * pc_y[k] * smg1_89[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, pc_x, pc_y, pc_z, smf_50, smf_90, smf_93, \
                         snd0_54, snd0_57, snd1_54, snd1_57, snf_90, \
                         snf_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_11 * smf_90[k]
                   + f_1 * snd0_54[k]
                   - f_2 * snd1_54[k]
                   + f_3 * pc_x[k] * snf_90[k];

        t_136[k] = f_3 * pc_y[k] * snf_90[k];

        t_137[k] = f_12 * smf_50[k]
                   + f_3 * pc_z[k] * snf_90[k];

        t_138[k] = f_11 * smf_93[k]
                   + f_4 * snd0_57[k]
                   - f_5 * snd1_57[k]
                   + f_3 * pc_x[k] * snf_93[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, pc_x, pc_y, smf_95, smf_96, smf_97, \
                         snd0_59, snd1_59, snf_92, snf_95, snf_96, \
                         snf_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_3 * pc_y[k] * snf_92[k];

        t_140[k] = f_11 * smf_95[k]
                   + f_4 * snd0_59[k]
                   - f_5 * snd1_59[k]
                   + f_3 * pc_x[k] * snf_95[k];

        t_141[k] = f_11 * smf_96[k]
                   + f_3 * pc_x[k] * snf_96[k];

        t_142[k] = f_11 * smf_97[k]
                   + f_3 * pc_x[k] * snf_97[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, t_146, pc_x, pc_y, pc_z, smf_56, smf_98, smf_99, \
                         snd0_57, snd1_57, snf_96, snf_98, snf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_11 * smf_98[k]
                   + f_3 * pc_x[k] * snf_98[k];

        t_144[k] = f_11 * smf_99[k]
                   + f_3 * pc_x[k] * snf_99[k];

        t_145[k] = f_1 * snd0_57[k]
                   - f_2 * snd1_57[k]
                   + f_3 * pc_y[k] * snf_96[k];

        t_146[k] = f_12 * smf_56[k]
                   + f_3 * pc_z[k] * snf_96[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, pc_x, pc_y, pc_z, smf_59, smf_100, \
                         snd0_59, snd0_60, snd1_59, snd1_60, snf_98, snf_99, \
                         snf_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = f_4 * snd0_59[k]
                   - f_5 * snd1_59[k]
                   + f_3 * pc_y[k] * snf_98[k];

        t_148[k] = f_3 * pc_y[k] * snf_99[k];

        t_149[k] = f_12 * smf_59[k]
                   + f_1 * snd0_59[k]
                   - f_2 * snd1_59[k]
                   + f_3 * pc_z[k] * snf_99[k];

        t_150[k] = f_13 * smf_100[k]
                   + f_1 * snd0_60[k]
                   - f_2 * snd1_60[k]
                   + f_3 * pc_x[k] * snf_100[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, pc_x, pc_y, pc_z, smf_60, smf_62, \
                         smf_103, snd0_63, snd1_63, snf_100, snf_102, \
                         snf_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_14 * smf_60[k]
                   + f_3 * pc_y[k] * snf_100[k];

        t_152[k] = f_3 * pc_z[k] * snf_100[k];

        t_153[k] = f_13 * smf_103[k]
                   + f_4 * snd0_63[k]
                   - f_5 * snd1_63[k]
                   + f_3 * pc_x[k] * snf_103[k];

        t_154[k] = f_14 * smf_62[k]
                   + f_3 * pc_y[k] * snf_102[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, pc_x, smf_105, smf_106, smf_107, smf_108, \
                         snd0_65, snd1_65, snf_105, snf_106, snf_107, \
                         snf_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_13 * smf_105[k]
                   + f_4 * snd0_65[k]
                   - f_5 * snd1_65[k]
                   + f_3 * pc_x[k] * snf_105[k];

        t_156[k] = f_13 * smf_106[k]
                   + f_3 * pc_x[k] * snf_106[k];

        t_157[k] = f_13 * smf_107[k]
                   + f_3 * pc_x[k] * snf_107[k];

        t_158[k] = f_13 * smf_108[k]
                   + f_3 * pc_x[k] * snf_108[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, pc_x, pc_y, pc_z, smf_66, smf_109, snd0_63, \
                         snd1_63, snf_106, snf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_13 * smf_109[k]
                   + f_3 * pc_x[k] * snf_109[k];

        t_160[k] = f_14 * smf_66[k]
                   + f_1 * snd0_63[k]
                   - f_2 * snd1_63[k]
                   + f_3 * pc_y[k] * snf_106[k];

        t_161[k] = f_3 * pc_z[k] * snf_106[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, pb_z, pc_y, pc_z, smg0_90, smf_68, \
                         smf_69, smg1_90, snd0_65, snd1_65, snf_108, \
                         snf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_14 * smf_68[k]
                   + f_4 * snd0_65[k]
                   - f_5 * snd1_65[k]
                   + f_3 * pc_y[k] * snf_108[k];

        t_163[k] = f_14 * smf_69[k]
                   + f_3 * pc_y[k] * snf_109[k];

        t_164[k] = f_1 * snd0_65[k]
                   - f_2 * snd1_65[k]
                   + f_3 * pc_z[k] * snf_109[k];

        t_165[k] = pb_z[k] * smg0_90[k]
                   - f_6 * pc_z[k] * smg1_90[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, pb_z, pc_y, pc_z, smg0_93, smf_60, \
                         smf_70, smf_72, smg1_93, snf_110, snf_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_12 * smf_70[k]
                   + f_3 * pc_y[k] * snf_110[k];

        t_167[k] = f_7 * smf_60[k]
                   + f_3 * pc_z[k] * snf_110[k];

        t_168[k] = pb_z[k] * smg0_93[k]
                   - f_6 * pc_z[k] * smg1_93[k];

        t_169[k] = f_12 * smf_72[k]
                   + f_3 * pc_y[k] * snf_112[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, pc_x, smf_115, smf_116, smf_117, smf_118, \
                         snd0_71, snd1_71, snf_115, snf_116, snf_117, \
                         snf_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_13 * smf_115[k]
                   + f_4 * snd0_71[k]
                   - f_5 * snd1_71[k]
                   + f_3 * pc_x[k] * snf_115[k];

        t_171[k] = f_13 * smf_116[k]
                   + f_3 * pc_x[k] * snf_116[k];

        t_172[k] = f_13 * smf_117[k]
                   + f_3 * pc_x[k] * snf_117[k];

        t_173[k] = f_13 * smf_118[k]
                   + f_3 * pc_x[k] * snf_118[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, pb_z, pc_x, pc_z, smg0_100, smf_66, smf_119, \
                         smg1_100, snf_116, snf_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_13 * smf_119[k]
                   + f_3 * pc_x[k] * snf_119[k];

        t_175[k] = pb_z[k] * smg0_100[k]
                   - f_6 * pc_z[k] * smg1_100[k];

        t_176[k] = f_7 * smf_66[k]
                   + f_3 * pc_z[k] * snf_116[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, pc_y, pc_z, smf_69, smf_78, smf_79, snd0_71, \
                         snd1_71, snf_118, snf_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = f_12 * smf_78[k]
                   + f_4 * snd0_71[k]
                   - f_5 * snd1_71[k]
                   + f_3 * pc_y[k] * snf_118[k];

        t_178[k] = f_12 * smf_79[k]
                   + f_3 * pc_y[k] * snf_119[k];

        t_179[k] = f_7 * smf_69[k]
                   + f_1 * snd0_71[k]
                   - f_2 * snd1_71[k]
                   + f_3 * pc_z[k] * snf_119[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, pc_x, pc_y, pc_z, smf_70, smf_80, smf_120, \
                         snd0_72, snd1_72, snf_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_13 * smf_120[k]
                   + f_1 * snd0_72[k]
                   - f_2 * snd1_72[k]
                   + f_3 * pc_x[k] * snf_120[k];

        t_181[k] = f_8 * smf_80[k]
                   + f_3 * pc_y[k] * snf_120[k];

        t_182[k] = f_8 * smf_70[k]
                   + f_3 * pc_z[k] * snf_120[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, pc_x, pc_y, smf_82, smf_123, smf_125, snd0_75, \
                         snd0_77, snd1_75, snd1_77, snf_122, snf_123, \
                         snf_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_13 * smf_123[k]
                   + f_4 * snd0_75[k]
                   - f_5 * snd1_75[k]
                   + f_3 * pc_x[k] * snf_123[k];

        t_184[k] = f_8 * smf_82[k]
                   + f_3 * pc_y[k] * snf_122[k];

        t_185[k] = f_13 * smf_125[k]
                   + f_4 * snd0_77[k]
                   - f_5 * snd1_77[k]
                   + f_3 * pc_x[k] * snf_125[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, t_189, pc_x, smf_126, smf_127, smf_128, smf_129, \
                         snf_126, snf_127, snf_128, snf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_13 * smf_126[k]
                   + f_3 * pc_x[k] * snf_126[k];

        t_187[k] = f_13 * smf_127[k]
                   + f_3 * pc_x[k] * snf_127[k];

        t_188[k] = f_13 * smf_128[k]
                   + f_3 * pc_x[k] * snf_128[k];

        t_189[k] = f_13 * smf_129[k]
                   + f_3 * pc_x[k] * snf_129[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, pc_y, pc_z, smf_76, smf_86, smf_88, snd0_75, \
                         snd0_77, snd1_75, snd1_77, snf_126, snf_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = f_8 * smf_86[k]
                   + f_1 * snd0_75[k]
                   - f_2 * snd1_75[k]
                   + f_3 * pc_y[k] * snf_126[k];

        t_191[k] = f_8 * smf_76[k]
                   + f_3 * pc_z[k] * snf_126[k];

        t_192[k] = f_8 * smf_88[k]
                   + f_4 * snd0_77[k]
                   - f_5 * snd1_77[k]
                   + f_3 * pc_y[k] * snf_128[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, pb_y, pc_y, pc_z, smg0_135, smf_79, \
                         smf_89, smf_90, smg1_135, snd0_77, snd1_77, snf_129, \
                         snf_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_8 * smf_89[k]
                   + f_3 * pc_y[k] * snf_129[k];

        t_194[k] = f_8 * smf_79[k]
                   + f_1 * snd0_77[k]
                   - f_2 * snd1_77[k]
                   + f_3 * pc_z[k] * snf_129[k];

        t_195[k] = pb_y[k] * smg0_135[k]
                   - f_6 * pc_y[k] * smg1_135[k];

        t_196[k] = f_7 * smf_90[k]
                   + f_3 * pc_y[k] * snf_130[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, pb_y, pc_y, pc_z, smg0_138, smg0_140, \
                         smf_80, smf_91, smf_92, smg1_138, smg1_140, snf_130, \
                         snf_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_12 * smf_80[k]
                   + f_3 * pc_z[k] * snf_130[k];

        t_198[k] = pb_y[k] * smg0_138[k]
                   + f_8 * smf_91[k]
                   - f_6 * pc_y[k] * smg1_138[k];

        t_199[k] = f_7 * smf_92[k]
                   + f_3 * pc_y[k] * snf_132[k];

        t_200[k] = pb_y[k] * smg0_140[k]
                   - f_6 * pc_y[k] * smg1_140[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, t_204, pc_x, smf_136, smf_137, smf_138, smf_139, \
                         snf_136, snf_137, snf_138, snf_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = f_13 * smf_136[k]
                   + f_3 * pc_x[k] * snf_136[k];

        t_202[k] = f_13 * smf_137[k]
                   + f_3 * pc_x[k] * snf_137[k];

        t_203[k] = f_13 * smf_138[k]
                   + f_3 * pc_x[k] * snf_138[k];

        t_204[k] = f_13 * smf_139[k]
                   + f_3 * pc_x[k] * snf_139[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, pc_y, pc_z, smf_86, smf_96, smf_98, snd0_81, \
                         snd0_83, snd1_81, snd1_83, snf_136, snf_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = f_7 * smf_96[k]
                   + f_1 * snd0_81[k]
                   - f_2 * snd1_81[k]
                   + f_3 * pc_y[k] * snf_136[k];

        t_206[k] = f_12 * smf_86[k]
                   + f_3 * pc_z[k] * snf_136[k];

        t_207[k] = f_7 * smf_98[k]
                   + f_4 * snd0_83[k]
                   - f_5 * snd1_83[k]
                   + f_3 * pc_y[k] * snf_138[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, pb_y, pc_x, pc_y, smg0_149, smf_99, \
                         smf_140, smg1_149, snd0_84, snd1_84, snf_139, \
                         snf_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_7 * smf_99[k]
                   + f_3 * pc_y[k] * snf_139[k];

        t_209[k] = pb_y[k] * smg0_149[k]
                   - f_6 * pc_y[k] * smg1_149[k];

        t_210[k] = f_13 * smf_140[k]
                   + f_1 * snd0_84[k]
                   - f_2 * snd1_84[k]
                   + f_3 * pc_x[k] * snf_140[k];

        t_211[k] = f_3 * pc_y[k] * snf_140[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, pc_x, pc_y, pc_z, smf_90, smf_143, snd0_87, \
                         snd1_87, snf_140, snf_142, snf_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_14 * smf_90[k]
                   + f_3 * pc_z[k] * snf_140[k];

        t_213[k] = f_13 * smf_143[k]
                   + f_4 * snd0_87[k]
                   - f_5 * snd1_87[k]
                   + f_3 * pc_x[k] * snf_143[k];

        t_214[k] = f_3 * pc_y[k] * snf_142[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, pc_x, smf_145, smf_146, smf_147, smf_148, \
                         snd0_89, snd1_89, snf_145, snf_146, snf_147, \
                         snf_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = f_13 * smf_145[k]
                   + f_4 * snd0_89[k]
                   - f_5 * snd1_89[k]
                   + f_3 * pc_x[k] * snf_145[k];

        t_216[k] = f_13 * smf_146[k]
                   + f_3 * pc_x[k] * snf_146[k];

        t_217[k] = f_13 * smf_147[k]
                   + f_3 * pc_x[k] * snf_147[k];

        t_218[k] = f_13 * smf_148[k]
                   + f_3 * pc_x[k] * snf_148[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, t_222, pc_x, pc_y, pc_z, smf_96, smf_149, \
                         snd0_87, snd0_89, snd1_87, snd1_89, snf_146, snf_148, \
                         snf_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = f_13 * smf_149[k]
                   + f_3 * pc_x[k] * snf_149[k];

        t_220[k] = f_1 * snd0_87[k]
                   - f_2 * snd1_87[k]
                   + f_3 * pc_y[k] * snf_146[k];

        t_221[k] = f_14 * smf_96[k]
                   + f_3 * pc_z[k] * snf_146[k];

        t_222[k] = f_4 * snd0_89[k]
                   - f_5 * snd1_89[k]
                   + f_3 * pc_y[k] * snf_148[k];
    }

#pragma omp simd aligned(t_223, t_224, t_225, t_226, pc_x, pc_y, pc_z, smf_99, smf_100, \
                         smf_150, snd0_89, snd0_90, snd1_89, snd1_90, snf_149, \
                         snf_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = f_3 * pc_y[k] * snf_149[k];

        t_224[k] = f_14 * smf_99[k]
                   + f_1 * snd0_89[k]
                   - f_2 * snd1_89[k]
                   + f_3 * pc_z[k] * snf_149[k];

        t_225[k] = f_15 * smf_150[k]
                   + f_1 * snd0_90[k]
                   - f_2 * snd1_90[k]
                   + f_3 * pc_x[k] * snf_150[k];

        t_226[k] = f_15 * smf_100[k]
                   + f_3 * pc_y[k] * snf_150[k];
    }

#pragma omp simd aligned(t_227, t_228, t_229, pc_x, pc_y, pc_z, smf_102, smf_153, snd0_93, \
                         snd1_93, snf_150, snf_152, snf_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_227[k] = f_3 * pc_z[k] * snf_150[k];

        t_228[k] = f_15 * smf_153[k]
                   + f_4 * snd0_93[k]
                   - f_5 * snd1_93[k]
                   + f_3 * pc_x[k] * snf_153[k];

        t_229[k] = f_15 * smf_102[k]
                   + f_3 * pc_y[k] * snf_152[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, pc_x, smf_155, smf_156, smf_157, smf_158, \
                         snd0_95, snd1_95, snf_155, snf_156, snf_157, \
                         snf_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_15 * smf_155[k]
                   + f_4 * snd0_95[k]
                   - f_5 * snd1_95[k]
                   + f_3 * pc_x[k] * snf_155[k];

        t_231[k] = f_15 * smf_156[k]
                   + f_3 * pc_x[k] * snf_156[k];

        t_232[k] = f_15 * smf_157[k]
                   + f_3 * pc_x[k] * snf_157[k];

        t_233[k] = f_15 * smf_158[k]
                   + f_3 * pc_x[k] * snf_158[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, pc_x, pc_y, pc_z, smf_106, smf_159, snd0_93, \
                         snd1_93, snf_156, snf_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_15 * smf_159[k]
                   + f_3 * pc_x[k] * snf_159[k];

        t_235[k] = f_15 * smf_106[k]
                   + f_1 * snd0_93[k]
                   - f_2 * snd1_93[k]
                   + f_3 * pc_y[k] * snf_156[k];

        t_236[k] = f_3 * pc_z[k] * snf_156[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, t_240, pb_z, pc_y, pc_z, smg0_150, smf_108, \
                         smf_109, smg1_150, snd0_95, snd1_95, snf_158, \
                         snf_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = f_15 * smf_108[k]
                   + f_4 * snd0_95[k]
                   - f_5 * snd1_95[k]
                   + f_3 * pc_y[k] * snf_158[k];

        t_238[k] = f_15 * smf_109[k]
                   + f_3 * pc_y[k] * snf_159[k];

        t_239[k] = f_1 * snd0_95[k]
                   - f_2 * snd1_95[k]
                   + f_3 * pc_z[k] * snf_159[k];

        t_240[k] = pb_z[k] * smg0_150[k]
                   - f_6 * pc_z[k] * smg1_150[k];
    }

#pragma omp simd aligned(t_241, t_242, t_243, t_244, pb_z, pc_y, pc_z, smg0_153, smf_100, \
                         smf_110, smf_112, smg1_153, snf_160, snf_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_241[k] = f_14 * smf_110[k]
                   + f_3 * pc_y[k] * snf_160[k];

        t_242[k] = f_7 * smf_100[k]
                   + f_3 * pc_z[k] * snf_160[k];

        t_243[k] = pb_z[k] * smg0_153[k]
                   - f_6 * pc_z[k] * smg1_153[k];

        t_244[k] = f_14 * smf_112[k]
                   + f_3 * pc_y[k] * snf_162[k];
    }
}

static auto
compute_prim_sng_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t smg0,
                                                          const size_t smf, const size_t smg1,
                                                          const size_t snd0, const size_t snd1,
                                                          const size_t snf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 3.0 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 2.5 / q;

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
    auto *t_358 = buffer.data(target + 358);
    auto *t_359 = buffer.data(target + 359);
    auto *t_360 = buffer.data(target + 360);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *smg0_160 = buffer.data(smg0 + 160);
    const auto *smg0_210 = buffer.data(smg0 + 210);
    const auto *smg0_213 = buffer.data(smg0 + 213);
    const auto *smg0_215 = buffer.data(smg0 + 215);
    const auto *smg0_224 = buffer.data(smg0 + 224);
    const auto *smg0_225 = buffer.data(smg0 + 225);
    const auto *smg0_228 = buffer.data(smg0 + 228);
    const auto *smg0_235 = buffer.data(smg0 + 235);

    const auto *smf_106 = buffer.data(smf + 106);
    const auto *smf_109 = buffer.data(smf + 109);
    const auto *smf_110 = buffer.data(smf + 110);
    const auto *smf_116 = buffer.data(smf + 116);
    const auto *smf_118 = buffer.data(smf + 118);
    const auto *smf_119 = buffer.data(smf + 119);
    const auto *smf_120 = buffer.data(smf + 120);
    const auto *smf_122 = buffer.data(smf + 122);
    const auto *smf_126 = buffer.data(smf + 126);
    const auto *smf_128 = buffer.data(smf + 128);
    const auto *smf_129 = buffer.data(smf + 129);
    const auto *smf_130 = buffer.data(smf + 130);
    const auto *smf_132 = buffer.data(smf + 132);
    const auto *smf_136 = buffer.data(smf + 136);
    const auto *smf_138 = buffer.data(smf + 138);
    const auto *smf_139 = buffer.data(smf + 139);
    const auto *smf_140 = buffer.data(smf + 140);
    const auto *smf_141 = buffer.data(smf + 141);
    const auto *smf_142 = buffer.data(smf + 142);
    const auto *smf_146 = buffer.data(smf + 146);
    const auto *smf_148 = buffer.data(smf + 148);
    const auto *smf_149 = buffer.data(smf + 149);
    const auto *smf_150 = buffer.data(smf + 150);
    const auto *smf_152 = buffer.data(smf + 152);
    const auto *smf_156 = buffer.data(smf + 156);
    const auto *smf_158 = buffer.data(smf + 158);
    const auto *smf_159 = buffer.data(smf + 159);
    const auto *smf_160 = buffer.data(smf + 160);
    const auto *smf_162 = buffer.data(smf + 162);
    const auto *smf_165 = buffer.data(smf + 165);
    const auto *smf_166 = buffer.data(smf + 166);
    const auto *smf_167 = buffer.data(smf + 167);
    const auto *smf_168 = buffer.data(smf + 168);
    const auto *smf_169 = buffer.data(smf + 169);
    const auto *smf_170 = buffer.data(smf + 170);
    const auto *smf_172 = buffer.data(smf + 172);
    const auto *smf_173 = buffer.data(smf + 173);
    const auto *smf_175 = buffer.data(smf + 175);
    const auto *smf_176 = buffer.data(smf + 176);
    const auto *smf_177 = buffer.data(smf + 177);
    const auto *smf_178 = buffer.data(smf + 178);
    const auto *smf_179 = buffer.data(smf + 179);
    const auto *smf_180 = buffer.data(smf + 180);
    const auto *smf_183 = buffer.data(smf + 183);
    const auto *smf_185 = buffer.data(smf + 185);
    const auto *smf_186 = buffer.data(smf + 186);
    const auto *smf_187 = buffer.data(smf + 187);
    const auto *smf_188 = buffer.data(smf + 188);
    const auto *smf_189 = buffer.data(smf + 189);
    const auto *smf_196 = buffer.data(smf + 196);
    const auto *smf_197 = buffer.data(smf + 197);
    const auto *smf_198 = buffer.data(smf + 198);
    const auto *smf_199 = buffer.data(smf + 199);
    const auto *smf_200 = buffer.data(smf + 200);
    const auto *smf_203 = buffer.data(smf + 203);
    const auto *smf_205 = buffer.data(smf + 205);
    const auto *smf_206 = buffer.data(smf + 206);
    const auto *smf_207 = buffer.data(smf + 207);
    const auto *smf_208 = buffer.data(smf + 208);
    const auto *smf_209 = buffer.data(smf + 209);
    const auto *smf_210 = buffer.data(smf + 210);
    const auto *smf_213 = buffer.data(smf + 213);
    const auto *smf_215 = buffer.data(smf + 215);
    const auto *smf_216 = buffer.data(smf + 216);
    const auto *smf_217 = buffer.data(smf + 217);
    const auto *smf_218 = buffer.data(smf + 218);
    const auto *smf_219 = buffer.data(smf + 219);
    const auto *smf_225 = buffer.data(smf + 225);
    const auto *smf_226 = buffer.data(smf + 226);
    const auto *smf_227 = buffer.data(smf + 227);
    const auto *smf_228 = buffer.data(smf + 228);
    const auto *smf_229 = buffer.data(smf + 229);
    const auto *smf_230 = buffer.data(smf + 230);
    const auto *smf_233 = buffer.data(smf + 233);
    const auto *smf_235 = buffer.data(smf + 235);
    const auto *smf_236 = buffer.data(smf + 236);
    const auto *smf_237 = buffer.data(smf + 237);
    const auto *smf_238 = buffer.data(smf + 238);
    const auto *smf_239 = buffer.data(smf + 239);
    const auto *smf_240 = buffer.data(smf + 240);

    const auto *smg1_160 = buffer.data(smg1 + 160);
    const auto *smg1_210 = buffer.data(smg1 + 210);
    const auto *smg1_213 = buffer.data(smg1 + 213);
    const auto *smg1_215 = buffer.data(smg1 + 215);
    const auto *smg1_224 = buffer.data(smg1 + 224);
    const auto *smg1_225 = buffer.data(smg1 + 225);
    const auto *smg1_228 = buffer.data(smg1 + 228);
    const auto *smg1_235 = buffer.data(smg1 + 235);

    const auto *snd0_101 = buffer.data(snd0 + 101);
    const auto *snd0_102 = buffer.data(snd0 + 102);
    const auto *snd0_105 = buffer.data(snd0 + 105);
    const auto *snd0_107 = buffer.data(snd0 + 107);
    const auto *snd0_108 = buffer.data(snd0 + 108);
    const auto *snd0_111 = buffer.data(snd0 + 111);
    const auto *snd0_113 = buffer.data(snd0 + 113);
    const auto *snd0_117 = buffer.data(snd0 + 117);
    const auto *snd0_119 = buffer.data(snd0 + 119);
    const auto *snd0_120 = buffer.data(snd0 + 120);
    const auto *snd0_123 = buffer.data(snd0 + 123);
    const auto *snd0_125 = buffer.data(snd0 + 125);
    const auto *snd0_126 = buffer.data(snd0 + 126);
    const auto *snd0_129 = buffer.data(snd0 + 129);
    const auto *snd0_131 = buffer.data(snd0 + 131);
    const auto *snd0_137 = buffer.data(snd0 + 137);
    const auto *snd0_138 = buffer.data(snd0 + 138);
    const auto *snd0_141 = buffer.data(snd0 + 141);
    const auto *snd0_143 = buffer.data(snd0 + 143);
    const auto *snd0_144 = buffer.data(snd0 + 144);

    const auto *snd1_101 = buffer.data(snd1 + 101);
    const auto *snd1_102 = buffer.data(snd1 + 102);
    const auto *snd1_105 = buffer.data(snd1 + 105);
    const auto *snd1_107 = buffer.data(snd1 + 107);
    const auto *snd1_108 = buffer.data(snd1 + 108);
    const auto *snd1_111 = buffer.data(snd1 + 111);
    const auto *snd1_113 = buffer.data(snd1 + 113);
    const auto *snd1_117 = buffer.data(snd1 + 117);
    const auto *snd1_119 = buffer.data(snd1 + 119);
    const auto *snd1_120 = buffer.data(snd1 + 120);
    const auto *snd1_123 = buffer.data(snd1 + 123);
    const auto *snd1_125 = buffer.data(snd1 + 125);
    const auto *snd1_126 = buffer.data(snd1 + 126);
    const auto *snd1_129 = buffer.data(snd1 + 129);
    const auto *snd1_131 = buffer.data(snd1 + 131);
    const auto *snd1_137 = buffer.data(snd1 + 137);
    const auto *snd1_138 = buffer.data(snd1 + 138);
    const auto *snd1_141 = buffer.data(snd1 + 141);
    const auto *snd1_143 = buffer.data(snd1 + 143);
    const auto *snd1_144 = buffer.data(snd1 + 144);

    const auto *snf_165 = buffer.data(snf + 165);
    const auto *snf_166 = buffer.data(snf + 166);
    const auto *snf_167 = buffer.data(snf + 167);
    const auto *snf_168 = buffer.data(snf + 168);
    const auto *snf_169 = buffer.data(snf + 169);
    const auto *snf_170 = buffer.data(snf + 170);
    const auto *snf_172 = buffer.data(snf + 172);
    const auto *snf_173 = buffer.data(snf + 173);
    const auto *snf_175 = buffer.data(snf + 175);
    const auto *snf_176 = buffer.data(snf + 176);
    const auto *snf_177 = buffer.data(snf + 177);
    const auto *snf_178 = buffer.data(snf + 178);
    const auto *snf_179 = buffer.data(snf + 179);
    const auto *snf_180 = buffer.data(snf + 180);
    const auto *snf_182 = buffer.data(snf + 182);
    const auto *snf_183 = buffer.data(snf + 183);
    const auto *snf_185 = buffer.data(snf + 185);
    const auto *snf_186 = buffer.data(snf + 186);
    const auto *snf_187 = buffer.data(snf + 187);
    const auto *snf_188 = buffer.data(snf + 188);
    const auto *snf_189 = buffer.data(snf + 189);
    const auto *snf_190 = buffer.data(snf + 190);
    const auto *snf_192 = buffer.data(snf + 192);
    const auto *snf_196 = buffer.data(snf + 196);
    const auto *snf_197 = buffer.data(snf + 197);
    const auto *snf_198 = buffer.data(snf + 198);
    const auto *snf_199 = buffer.data(snf + 199);
    const auto *snf_200 = buffer.data(snf + 200);
    const auto *snf_202 = buffer.data(snf + 202);
    const auto *snf_203 = buffer.data(snf + 203);
    const auto *snf_205 = buffer.data(snf + 205);
    const auto *snf_206 = buffer.data(snf + 206);
    const auto *snf_207 = buffer.data(snf + 207);
    const auto *snf_208 = buffer.data(snf + 208);
    const auto *snf_209 = buffer.data(snf + 209);
    const auto *snf_210 = buffer.data(snf + 210);
    const auto *snf_212 = buffer.data(snf + 212);
    const auto *snf_213 = buffer.data(snf + 213);
    const auto *snf_215 = buffer.data(snf + 215);
    const auto *snf_216 = buffer.data(snf + 216);
    const auto *snf_217 = buffer.data(snf + 217);
    const auto *snf_218 = buffer.data(snf + 218);
    const auto *snf_219 = buffer.data(snf + 219);
    const auto *snf_220 = buffer.data(snf + 220);
    const auto *snf_222 = buffer.data(snf + 222);
    const auto *snf_225 = buffer.data(snf + 225);
    const auto *snf_226 = buffer.data(snf + 226);
    const auto *snf_227 = buffer.data(snf + 227);
    const auto *snf_228 = buffer.data(snf + 228);
    const auto *snf_229 = buffer.data(snf + 229);
    const auto *snf_230 = buffer.data(snf + 230);
    const auto *snf_232 = buffer.data(snf + 232);
    const auto *snf_233 = buffer.data(snf + 233);
    const auto *snf_235 = buffer.data(snf + 235);
    const auto *snf_236 = buffer.data(snf + 236);
    const auto *snf_237 = buffer.data(snf + 237);
    const auto *snf_238 = buffer.data(snf + 238);
    const auto *snf_239 = buffer.data(snf + 239);
    const auto *snf_240 = buffer.data(snf + 240);

#pragma omp simd aligned(t_245, t_246, t_247, t_248, pc_x, smf_165, smf_166, smf_167, smf_168, \
                         snd0_101, snd1_101, snf_165, snf_166, snf_167, \
                         snf_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_15 * smf_165[k]
                   + f_4 * snd0_101[k]
                   - f_5 * snd1_101[k]
                   + f_3 * pc_x[k] * snf_165[k];

        t_246[k] = f_15 * smf_166[k]
                   + f_3 * pc_x[k] * snf_166[k];

        t_247[k] = f_15 * smf_167[k]
                   + f_3 * pc_x[k] * snf_167[k];

        t_248[k] = f_15 * smf_168[k]
                   + f_3 * pc_x[k] * snf_168[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, pb_z, pc_x, pc_z, smg0_160, smf_106, smf_169, \
                         smg1_160, snf_166, snf_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = f_15 * smf_169[k]
                   + f_3 * pc_x[k] * snf_169[k];

        t_250[k] = pb_z[k] * smg0_160[k]
                   - f_6 * pc_z[k] * smg1_160[k];

        t_251[k] = f_7 * smf_106[k]
                   + f_3 * pc_z[k] * snf_166[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, pc_y, pc_z, smf_109, smf_118, smf_119, snd0_101, \
                         snd1_101, snf_168, snf_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = f_14 * smf_118[k]
                   + f_4 * snd0_101[k]
                   - f_5 * snd1_101[k]
                   + f_3 * pc_y[k] * snf_168[k];

        t_253[k] = f_14 * smf_119[k]
                   + f_3 * pc_y[k] * snf_169[k];

        t_254[k] = f_7 * smf_109[k]
                   + f_1 * snd0_101[k]
                   - f_2 * snd1_101[k]
                   + f_3 * pc_z[k] * snf_169[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, pc_x, pc_y, pc_z, smf_110, smf_120, smf_170, \
                         snd0_102, snd1_102, snf_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_15 * smf_170[k]
                   + f_1 * snd0_102[k]
                   - f_2 * snd1_102[k]
                   + f_3 * pc_x[k] * snf_170[k];

        t_256[k] = f_12 * smf_120[k]
                   + f_3 * pc_y[k] * snf_170[k];

        t_257[k] = f_8 * smf_110[k]
                   + f_3 * pc_z[k] * snf_170[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, pc_x, pc_y, smf_122, smf_173, smf_175, snd0_105, \
                         snd0_107, snd1_105, snd1_107, snf_172, snf_173, \
                         snf_175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = f_15 * smf_173[k]
                   + f_4 * snd0_105[k]
                   - f_5 * snd1_105[k]
                   + f_3 * pc_x[k] * snf_173[k];

        t_259[k] = f_12 * smf_122[k]
                   + f_3 * pc_y[k] * snf_172[k];

        t_260[k] = f_15 * smf_175[k]
                   + f_4 * snd0_107[k]
                   - f_5 * snd1_107[k]
                   + f_3 * pc_x[k] * snf_175[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, t_264, pc_x, smf_176, smf_177, smf_178, smf_179, \
                         snf_176, snf_177, snf_178, snf_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_15 * smf_176[k]
                   + f_3 * pc_x[k] * snf_176[k];

        t_262[k] = f_15 * smf_177[k]
                   + f_3 * pc_x[k] * snf_177[k];

        t_263[k] = f_15 * smf_178[k]
                   + f_3 * pc_x[k] * snf_178[k];

        t_264[k] = f_15 * smf_179[k]
                   + f_3 * pc_x[k] * snf_179[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, pc_y, pc_z, smf_116, smf_126, smf_128, snd0_105, \
                         snd0_107, snd1_105, snd1_107, snf_176, \
                         snf_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = f_12 * smf_126[k]
                   + f_1 * snd0_105[k]
                   - f_2 * snd1_105[k]
                   + f_3 * pc_y[k] * snf_176[k];

        t_266[k] = f_8 * smf_116[k]
                   + f_3 * pc_z[k] * snf_176[k];

        t_267[k] = f_12 * smf_128[k]
                   + f_4 * snd0_107[k]
                   - f_5 * snd1_107[k]
                   + f_3 * pc_y[k] * snf_178[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, pc_x, pc_y, pc_z, smf_119, smf_129, smf_180, \
                         snd0_107, snd0_108, snd1_107, snd1_108, snf_179, \
                         snf_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_12 * smf_129[k]
                   + f_3 * pc_y[k] * snf_179[k];

        t_269[k] = f_8 * smf_119[k]
                   + f_1 * snd0_107[k]
                   - f_2 * snd1_107[k]
                   + f_3 * pc_z[k] * snf_179[k];

        t_270[k] = f_15 * smf_180[k]
                   + f_1 * snd0_108[k]
                   - f_2 * snd1_108[k]
                   + f_3 * pc_x[k] * snf_180[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, t_274, pc_x, pc_y, pc_z, smf_120, smf_130, \
                         smf_132, smf_183, snd0_111, snd1_111, snf_180, snf_182, \
                         snf_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = f_8 * smf_130[k]
                   + f_3 * pc_y[k] * snf_180[k];

        t_272[k] = f_12 * smf_120[k]
                   + f_3 * pc_z[k] * snf_180[k];

        t_273[k] = f_15 * smf_183[k]
                   + f_4 * snd0_111[k]
                   - f_5 * snd1_111[k]
                   + f_3 * pc_x[k] * snf_183[k];

        t_274[k] = f_8 * smf_132[k]
                   + f_3 * pc_y[k] * snf_182[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, pc_x, smf_185, smf_186, smf_187, smf_188, \
                         snd0_113, snd1_113, snf_185, snf_186, snf_187, \
                         snf_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = f_15 * smf_185[k]
                   + f_4 * snd0_113[k]
                   - f_5 * snd1_113[k]
                   + f_3 * pc_x[k] * snf_185[k];

        t_276[k] = f_15 * smf_186[k]
                   + f_3 * pc_x[k] * snf_186[k];

        t_277[k] = f_15 * smf_187[k]
                   + f_3 * pc_x[k] * snf_187[k];

        t_278[k] = f_15 * smf_188[k]
                   + f_3 * pc_x[k] * snf_188[k];
    }

#pragma omp simd aligned(t_279, t_280, t_281, pc_x, pc_y, pc_z, smf_126, smf_136, smf_189, \
                         snd0_111, snd1_111, snf_186, snf_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_279[k] = f_15 * smf_189[k]
                   + f_3 * pc_x[k] * snf_189[k];

        t_280[k] = f_8 * smf_136[k]
                   + f_1 * snd0_111[k]
                   - f_2 * snd1_111[k]
                   + f_3 * pc_y[k] * snf_186[k];

        t_281[k] = f_12 * smf_126[k]
                   + f_3 * pc_z[k] * snf_186[k];
    }

#pragma omp simd aligned(t_282, t_283, t_284, t_285, pb_y, pc_y, pc_z, smg0_210, smf_129, \
                         smf_138, smf_139, smg1_210, snd0_113, snd1_113, snf_188, \
                         snf_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_282[k] = f_8 * smf_138[k]
                   + f_4 * snd0_113[k]
                   - f_5 * snd1_113[k]
                   + f_3 * pc_y[k] * snf_188[k];

        t_283[k] = f_8 * smf_139[k]
                   + f_3 * pc_y[k] * snf_189[k];

        t_284[k] = f_12 * smf_129[k]
                   + f_1 * snd0_113[k]
                   - f_2 * snd1_113[k]
                   + f_3 * pc_z[k] * snf_189[k];

        t_285[k] = pb_y[k] * smg0_210[k]
                   - f_6 * pc_y[k] * smg1_210[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, pb_y, pc_y, pc_z, smg0_213, smf_130, \
                         smf_140, smf_141, smf_142, smg1_213, snf_190, \
                         snf_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = f_7 * smf_140[k]
                   + f_3 * pc_y[k] * snf_190[k];

        t_287[k] = f_14 * smf_130[k]
                   + f_3 * pc_z[k] * snf_190[k];

        t_288[k] = pb_y[k] * smg0_213[k]
                   + f_8 * smf_141[k]
                   - f_6 * pc_y[k] * smg1_213[k];

        t_289[k] = f_7 * smf_142[k]
                   + f_3 * pc_y[k] * snf_192[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, pb_y, pc_x, pc_y, smg0_215, smf_196, \
                         smf_197, smf_198, smg1_215, snf_196, snf_197, \
                         snf_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = pb_y[k] * smg0_215[k]
                   - f_6 * pc_y[k] * smg1_215[k];

        t_291[k] = f_15 * smf_196[k]
                   + f_3 * pc_x[k] * snf_196[k];

        t_292[k] = f_15 * smf_197[k]
                   + f_3 * pc_x[k] * snf_197[k];

        t_293[k] = f_15 * smf_198[k]
                   + f_3 * pc_x[k] * snf_198[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, pc_x, pc_y, pc_z, smf_136, smf_146, smf_199, \
                         snd0_117, snd1_117, snf_196, snf_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_15 * smf_199[k]
                   + f_3 * pc_x[k] * snf_199[k];

        t_295[k] = f_7 * smf_146[k]
                   + f_1 * snd0_117[k]
                   - f_2 * snd1_117[k]
                   + f_3 * pc_y[k] * snf_196[k];

        t_296[k] = f_14 * smf_136[k]
                   + f_3 * pc_z[k] * snf_196[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, pb_y, pc_y, smg0_224, smf_148, smf_149, \
                         smg1_224, snd0_119, snd1_119, snf_198, \
                         snf_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = f_7 * smf_148[k]
                   + f_4 * snd0_119[k]
                   - f_5 * snd1_119[k]
                   + f_3 * pc_y[k] * snf_198[k];

        t_298[k] = f_7 * smf_149[k]
                   + f_3 * pc_y[k] * snf_199[k];

        t_299[k] = pb_y[k] * smg0_224[k]
                   - f_6 * pc_y[k] * smg1_224[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, pc_x, pc_y, pc_z, smf_140, smf_200, \
                         smf_203, snd0_120, snd0_123, snd1_120, snd1_123, snf_200, \
                         snf_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = f_15 * smf_200[k]
                   + f_1 * snd0_120[k]
                   - f_2 * snd1_120[k]
                   + f_3 * pc_x[k] * snf_200[k];

        t_301[k] = f_3 * pc_y[k] * snf_200[k];

        t_302[k] = f_15 * smf_140[k]
                   + f_3 * pc_z[k] * snf_200[k];

        t_303[k] = f_15 * smf_203[k]
                   + f_4 * snd0_123[k]
                   - f_5 * snd1_123[k]
                   + f_3 * pc_x[k] * snf_203[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, pc_x, pc_y, smf_205, smf_206, smf_207, \
                         snd0_125, snd1_125, snf_202, snf_205, snf_206, \
                         snf_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_304[k] = f_3 * pc_y[k] * snf_202[k];

        t_305[k] = f_15 * smf_205[k]
                   + f_4 * snd0_125[k]
                   - f_5 * snd1_125[k]
                   + f_3 * pc_x[k] * snf_205[k];

        t_306[k] = f_15 * smf_206[k]
                   + f_3 * pc_x[k] * snf_206[k];

        t_307[k] = f_15 * smf_207[k]
                   + f_3 * pc_x[k] * snf_207[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, t_311, pc_x, pc_y, pc_z, smf_146, smf_208, \
                         smf_209, snd0_123, snd1_123, snf_206, snf_208, \
                         snf_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = f_15 * smf_208[k]
                   + f_3 * pc_x[k] * snf_208[k];

        t_309[k] = f_15 * smf_209[k]
                   + f_3 * pc_x[k] * snf_209[k];

        t_310[k] = f_1 * snd0_123[k]
                   - f_2 * snd1_123[k]
                   + f_3 * pc_y[k] * snf_206[k];

        t_311[k] = f_15 * smf_146[k]
                   + f_3 * pc_z[k] * snf_206[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, t_315, pc_x, pc_y, pc_z, smf_149, smf_210, \
                         snd0_125, snd0_126, snd1_125, snd1_126, snf_208, snf_209, \
                         snf_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = f_4 * snd0_125[k]
                   - f_5 * snd1_125[k]
                   + f_3 * pc_y[k] * snf_208[k];

        t_313[k] = f_3 * pc_y[k] * snf_209[k];

        t_314[k] = f_15 * smf_149[k]
                   + f_1 * snd0_125[k]
                   - f_2 * snd1_125[k]
                   + f_3 * pc_z[k] * snf_209[k];

        t_315[k] = f_14 * smf_210[k]
                   + f_1 * snd0_126[k]
                   - f_2 * snd1_126[k]
                   + f_3 * pc_x[k] * snf_210[k];
    }

#pragma omp simd aligned(t_316, t_317, t_318, t_319, pc_x, pc_y, pc_z, smf_150, smf_152, \
                         smf_213, snd0_129, snd1_129, snf_210, snf_212, \
                         snf_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_316[k] = f_13 * smf_150[k]
                   + f_3 * pc_y[k] * snf_210[k];

        t_317[k] = f_3 * pc_z[k] * snf_210[k];

        t_318[k] = f_14 * smf_213[k]
                   + f_4 * snd0_129[k]
                   - f_5 * snd1_129[k]
                   + f_3 * pc_x[k] * snf_213[k];

        t_319[k] = f_13 * smf_152[k]
                   + f_3 * pc_y[k] * snf_212[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, pc_x, smf_215, smf_216, smf_217, smf_218, \
                         snd0_131, snd1_131, snf_215, snf_216, snf_217, \
                         snf_218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = f_14 * smf_215[k]
                   + f_4 * snd0_131[k]
                   - f_5 * snd1_131[k]
                   + f_3 * pc_x[k] * snf_215[k];

        t_321[k] = f_14 * smf_216[k]
                   + f_3 * pc_x[k] * snf_216[k];

        t_322[k] = f_14 * smf_217[k]
                   + f_3 * pc_x[k] * snf_217[k];

        t_323[k] = f_14 * smf_218[k]
                   + f_3 * pc_x[k] * snf_218[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, pc_x, pc_y, pc_z, smf_156, smf_219, snd0_129, \
                         snd1_129, snf_216, snf_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = f_14 * smf_219[k]
                   + f_3 * pc_x[k] * snf_219[k];

        t_325[k] = f_13 * smf_156[k]
                   + f_1 * snd0_129[k]
                   - f_2 * snd1_129[k]
                   + f_3 * pc_y[k] * snf_216[k];

        t_326[k] = f_3 * pc_z[k] * snf_216[k];
    }

#pragma omp simd aligned(t_327, t_328, t_329, t_330, pb_z, pc_y, pc_z, smg0_225, smf_158, \
                         smf_159, smg1_225, snd0_131, snd1_131, snf_218, \
                         snf_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_327[k] = f_13 * smf_158[k]
                   + f_4 * snd0_131[k]
                   - f_5 * snd1_131[k]
                   + f_3 * pc_y[k] * snf_218[k];

        t_328[k] = f_13 * smf_159[k]
                   + f_3 * pc_y[k] * snf_219[k];

        t_329[k] = f_1 * snd0_131[k]
                   - f_2 * snd1_131[k]
                   + f_3 * pc_z[k] * snf_219[k];

        t_330[k] = pb_z[k] * smg0_225[k]
                   - f_6 * pc_z[k] * smg1_225[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, t_334, pb_z, pc_y, pc_z, smg0_228, smf_150, \
                         smf_160, smf_162, smg1_228, snf_220, snf_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = f_15 * smf_160[k]
                   + f_3 * pc_y[k] * snf_220[k];

        t_332[k] = f_7 * smf_150[k]
                   + f_3 * pc_z[k] * snf_220[k];

        t_333[k] = pb_z[k] * smg0_228[k]
                   - f_6 * pc_z[k] * smg1_228[k];

        t_334[k] = f_15 * smf_162[k]
                   + f_3 * pc_y[k] * snf_222[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, pc_x, smf_225, smf_226, smf_227, smf_228, \
                         snd0_137, snd1_137, snf_225, snf_226, snf_227, \
                         snf_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = f_14 * smf_225[k]
                   + f_4 * snd0_137[k]
                   - f_5 * snd1_137[k]
                   + f_3 * pc_x[k] * snf_225[k];

        t_336[k] = f_14 * smf_226[k]
                   + f_3 * pc_x[k] * snf_226[k];

        t_337[k] = f_14 * smf_227[k]
                   + f_3 * pc_x[k] * snf_227[k];

        t_338[k] = f_14 * smf_228[k]
                   + f_3 * pc_x[k] * snf_228[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, pb_z, pc_x, pc_z, smg0_235, smf_156, smf_229, \
                         smg1_235, snf_226, snf_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = f_14 * smf_229[k]
                   + f_3 * pc_x[k] * snf_229[k];

        t_340[k] = pb_z[k] * smg0_235[k]
                   - f_6 * pc_z[k] * smg1_235[k];

        t_341[k] = f_7 * smf_156[k]
                   + f_3 * pc_z[k] * snf_226[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, pc_y, pc_z, smf_159, smf_168, smf_169, snd0_137, \
                         snd1_137, snf_228, snf_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = f_15 * smf_168[k]
                   + f_4 * snd0_137[k]
                   - f_5 * snd1_137[k]
                   + f_3 * pc_y[k] * snf_228[k];

        t_343[k] = f_15 * smf_169[k]
                   + f_3 * pc_y[k] * snf_229[k];

        t_344[k] = f_7 * smf_159[k]
                   + f_1 * snd0_137[k]
                   - f_2 * snd1_137[k]
                   + f_3 * pc_z[k] * snf_229[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, pc_x, pc_y, pc_z, smf_160, smf_170, smf_230, \
                         snd0_138, snd1_138, snf_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_14 * smf_230[k]
                   + f_1 * snd0_138[k]
                   - f_2 * snd1_138[k]
                   + f_3 * pc_x[k] * snf_230[k];

        t_346[k] = f_14 * smf_170[k]
                   + f_3 * pc_y[k] * snf_230[k];

        t_347[k] = f_8 * smf_160[k]
                   + f_3 * pc_z[k] * snf_230[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, pc_x, pc_y, smf_172, smf_233, smf_235, snd0_141, \
                         snd0_143, snd1_141, snd1_143, snf_232, snf_233, \
                         snf_235 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = f_14 * smf_233[k]
                   + f_4 * snd0_141[k]
                   - f_5 * snd1_141[k]
                   + f_3 * pc_x[k] * snf_233[k];

        t_349[k] = f_14 * smf_172[k]
                   + f_3 * pc_y[k] * snf_232[k];

        t_350[k] = f_14 * smf_235[k]
                   + f_4 * snd0_143[k]
                   - f_5 * snd1_143[k]
                   + f_3 * pc_x[k] * snf_235[k];
    }

#pragma omp simd aligned(t_351, t_352, t_353, t_354, pc_x, smf_236, smf_237, smf_238, smf_239, \
                         snf_236, snf_237, snf_238, snf_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_351[k] = f_14 * smf_236[k]
                   + f_3 * pc_x[k] * snf_236[k];

        t_352[k] = f_14 * smf_237[k]
                   + f_3 * pc_x[k] * snf_237[k];

        t_353[k] = f_14 * smf_238[k]
                   + f_3 * pc_x[k] * snf_238[k];

        t_354[k] = f_14 * smf_239[k]
                   + f_3 * pc_x[k] * snf_239[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, pc_y, pc_z, smf_166, smf_176, smf_178, snd0_141, \
                         snd0_143, snd1_141, snd1_143, snf_236, \
                         snf_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = f_14 * smf_176[k]
                   + f_1 * snd0_141[k]
                   - f_2 * snd1_141[k]
                   + f_3 * pc_y[k] * snf_236[k];

        t_356[k] = f_8 * smf_166[k]
                   + f_3 * pc_z[k] * snf_236[k];

        t_357[k] = f_14 * smf_178[k]
                   + f_4 * snd0_143[k]
                   - f_5 * snd1_143[k]
                   + f_3 * pc_y[k] * snf_238[k];
    }

#pragma omp simd aligned(t_358, t_359, t_360, pc_x, pc_y, pc_z, smf_169, smf_179, smf_240, \
                         snd0_143, snd0_144, snd1_143, snd1_144, snf_239, \
                         snf_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_358[k] = f_14 * smf_179[k]
                   + f_3 * pc_y[k] * snf_239[k];

        t_359[k] = f_8 * smf_169[k]
                   + f_1 * snd0_143[k]
                   - f_2 * snd1_143[k]
                   + f_3 * pc_z[k] * snf_239[k];

        t_360[k] = f_14 * smf_240[k]
                   + f_1 * snd0_144[k]
                   - f_2 * snd1_144[k]
                   + f_3 * pc_x[k] * snf_240[k];
    }
}

static auto
compute_prim_sng_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t smg0,
                                                          const size_t smf, const size_t smg1,
                                                          const size_t snd0, const size_t snd1,
                                                          const size_t snf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_11 = 3.5 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 3.0 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 2.5 / q;

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
    auto *t_474 = buffer.data(target + 474);
    auto *t_475 = buffer.data(target + 475);
    auto *t_476 = buffer.data(target + 476);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *smg0_300 = buffer.data(smg0 + 300);
    const auto *smg0_303 = buffer.data(smg0 + 303);
    const auto *smg0_305 = buffer.data(smg0 + 305);
    const auto *smg0_314 = buffer.data(smg0 + 314);
    const auto *smg0_315 = buffer.data(smg0 + 315);
    const auto *smg0_318 = buffer.data(smg0 + 318);
    const auto *smg0_325 = buffer.data(smg0 + 325);

    const auto *smf_170 = buffer.data(smf + 170);
    const auto *smf_176 = buffer.data(smf + 176);
    const auto *smf_179 = buffer.data(smf + 179);
    const auto *smf_180 = buffer.data(smf + 180);
    const auto *smf_182 = buffer.data(smf + 182);
    const auto *smf_186 = buffer.data(smf + 186);
    const auto *smf_188 = buffer.data(smf + 188);
    const auto *smf_189 = buffer.data(smf + 189);
    const auto *smf_190 = buffer.data(smf + 190);
    const auto *smf_192 = buffer.data(smf + 192);
    const auto *smf_196 = buffer.data(smf + 196);
    const auto *smf_198 = buffer.data(smf + 198);
    const auto *smf_199 = buffer.data(smf + 199);
    const auto *smf_200 = buffer.data(smf + 200);
    const auto *smf_201 = buffer.data(smf + 201);
    const auto *smf_202 = buffer.data(smf + 202);
    const auto *smf_206 = buffer.data(smf + 206);
    const auto *smf_208 = buffer.data(smf + 208);
    const auto *smf_209 = buffer.data(smf + 209);
    const auto *smf_210 = buffer.data(smf + 210);
    const auto *smf_212 = buffer.data(smf + 212);
    const auto *smf_216 = buffer.data(smf + 216);
    const auto *smf_218 = buffer.data(smf + 218);
    const auto *smf_219 = buffer.data(smf + 219);
    const auto *smf_220 = buffer.data(smf + 220);
    const auto *smf_222 = buffer.data(smf + 222);
    const auto *smf_226 = buffer.data(smf + 226);
    const auto *smf_228 = buffer.data(smf + 228);
    const auto *smf_229 = buffer.data(smf + 229);
    const auto *smf_230 = buffer.data(smf + 230);
    const auto *smf_232 = buffer.data(smf + 232);
    const auto *smf_236 = buffer.data(smf + 236);
    const auto *smf_238 = buffer.data(smf + 238);
    const auto *smf_239 = buffer.data(smf + 239);
    const auto *smf_240 = buffer.data(smf + 240);
    const auto *smf_242 = buffer.data(smf + 242);
    const auto *smf_243 = buffer.data(smf + 243);
    const auto *smf_245 = buffer.data(smf + 245);
    const auto *smf_246 = buffer.data(smf + 246);
    const auto *smf_247 = buffer.data(smf + 247);
    const auto *smf_248 = buffer.data(smf + 248);
    const auto *smf_249 = buffer.data(smf + 249);
    const auto *smf_250 = buffer.data(smf + 250);
    const auto *smf_253 = buffer.data(smf + 253);
    const auto *smf_255 = buffer.data(smf + 255);
    const auto *smf_256 = buffer.data(smf + 256);
    const auto *smf_257 = buffer.data(smf + 257);
    const auto *smf_258 = buffer.data(smf + 258);
    const auto *smf_259 = buffer.data(smf + 259);
    const auto *smf_266 = buffer.data(smf + 266);
    const auto *smf_267 = buffer.data(smf + 267);
    const auto *smf_268 = buffer.data(smf + 268);
    const auto *smf_269 = buffer.data(smf + 269);
    const auto *smf_270 = buffer.data(smf + 270);
    const auto *smf_273 = buffer.data(smf + 273);
    const auto *smf_275 = buffer.data(smf + 275);
    const auto *smf_276 = buffer.data(smf + 276);
    const auto *smf_277 = buffer.data(smf + 277);
    const auto *smf_278 = buffer.data(smf + 278);
    const auto *smf_279 = buffer.data(smf + 279);
    const auto *smf_280 = buffer.data(smf + 280);
    const auto *smf_283 = buffer.data(smf + 283);
    const auto *smf_285 = buffer.data(smf + 285);
    const auto *smf_286 = buffer.data(smf + 286);
    const auto *smf_287 = buffer.data(smf + 287);
    const auto *smf_288 = buffer.data(smf + 288);
    const auto *smf_289 = buffer.data(smf + 289);
    const auto *smf_295 = buffer.data(smf + 295);
    const auto *smf_296 = buffer.data(smf + 296);
    const auto *smf_297 = buffer.data(smf + 297);
    const auto *smf_298 = buffer.data(smf + 298);
    const auto *smf_299 = buffer.data(smf + 299);
    const auto *smf_300 = buffer.data(smf + 300);
    const auto *smf_303 = buffer.data(smf + 303);
    const auto *smf_305 = buffer.data(smf + 305);
    const auto *smf_306 = buffer.data(smf + 306);
    const auto *smf_307 = buffer.data(smf + 307);
    const auto *smf_308 = buffer.data(smf + 308);
    const auto *smf_309 = buffer.data(smf + 309);
    const auto *smf_310 = buffer.data(smf + 310);
    const auto *smf_313 = buffer.data(smf + 313);
    const auto *smf_315 = buffer.data(smf + 315);
    const auto *smf_316 = buffer.data(smf + 316);
    const auto *smf_317 = buffer.data(smf + 317);
    const auto *smf_318 = buffer.data(smf + 318);
    const auto *smf_319 = buffer.data(smf + 319);

    const auto *smg1_300 = buffer.data(smg1 + 300);
    const auto *smg1_303 = buffer.data(smg1 + 303);
    const auto *smg1_305 = buffer.data(smg1 + 305);
    const auto *smg1_314 = buffer.data(smg1 + 314);
    const auto *smg1_315 = buffer.data(smg1 + 315);
    const auto *smg1_318 = buffer.data(smg1 + 318);
    const auto *smg1_325 = buffer.data(smg1 + 325);

    const auto *snd0_147 = buffer.data(snd0 + 147);
    const auto *snd0_149 = buffer.data(snd0 + 149);
    const auto *snd0_150 = buffer.data(snd0 + 150);
    const auto *snd0_153 = buffer.data(snd0 + 153);
    const auto *snd0_155 = buffer.data(snd0 + 155);
    const auto *snd0_159 = buffer.data(snd0 + 159);
    const auto *snd0_161 = buffer.data(snd0 + 161);
    const auto *snd0_162 = buffer.data(snd0 + 162);
    const auto *snd0_165 = buffer.data(snd0 + 165);
    const auto *snd0_167 = buffer.data(snd0 + 167);
    const auto *snd0_168 = buffer.data(snd0 + 168);
    const auto *snd0_171 = buffer.data(snd0 + 171);
    const auto *snd0_173 = buffer.data(snd0 + 173);
    const auto *snd0_179 = buffer.data(snd0 + 179);
    const auto *snd0_180 = buffer.data(snd0 + 180);
    const auto *snd0_183 = buffer.data(snd0 + 183);
    const auto *snd0_185 = buffer.data(snd0 + 185);
    const auto *snd0_186 = buffer.data(snd0 + 186);
    const auto *snd0_189 = buffer.data(snd0 + 189);
    const auto *snd0_191 = buffer.data(snd0 + 191);

    const auto *snd1_147 = buffer.data(snd1 + 147);
    const auto *snd1_149 = buffer.data(snd1 + 149);
    const auto *snd1_150 = buffer.data(snd1 + 150);
    const auto *snd1_153 = buffer.data(snd1 + 153);
    const auto *snd1_155 = buffer.data(snd1 + 155);
    const auto *snd1_159 = buffer.data(snd1 + 159);
    const auto *snd1_161 = buffer.data(snd1 + 161);
    const auto *snd1_162 = buffer.data(snd1 + 162);
    const auto *snd1_165 = buffer.data(snd1 + 165);
    const auto *snd1_167 = buffer.data(snd1 + 167);
    const auto *snd1_168 = buffer.data(snd1 + 168);
    const auto *snd1_171 = buffer.data(snd1 + 171);
    const auto *snd1_173 = buffer.data(snd1 + 173);
    const auto *snd1_179 = buffer.data(snd1 + 179);
    const auto *snd1_180 = buffer.data(snd1 + 180);
    const auto *snd1_183 = buffer.data(snd1 + 183);
    const auto *snd1_185 = buffer.data(snd1 + 185);
    const auto *snd1_186 = buffer.data(snd1 + 186);
    const auto *snd1_189 = buffer.data(snd1 + 189);
    const auto *snd1_191 = buffer.data(snd1 + 191);

    const auto *snf_240 = buffer.data(snf + 240);
    const auto *snf_242 = buffer.data(snf + 242);
    const auto *snf_243 = buffer.data(snf + 243);
    const auto *snf_245 = buffer.data(snf + 245);
    const auto *snf_246 = buffer.data(snf + 246);
    const auto *snf_247 = buffer.data(snf + 247);
    const auto *snf_248 = buffer.data(snf + 248);
    const auto *snf_249 = buffer.data(snf + 249);
    const auto *snf_250 = buffer.data(snf + 250);
    const auto *snf_252 = buffer.data(snf + 252);
    const auto *snf_253 = buffer.data(snf + 253);
    const auto *snf_255 = buffer.data(snf + 255);
    const auto *snf_256 = buffer.data(snf + 256);
    const auto *snf_257 = buffer.data(snf + 257);
    const auto *snf_258 = buffer.data(snf + 258);
    const auto *snf_259 = buffer.data(snf + 259);
    const auto *snf_260 = buffer.data(snf + 260);
    const auto *snf_262 = buffer.data(snf + 262);
    const auto *snf_266 = buffer.data(snf + 266);
    const auto *snf_267 = buffer.data(snf + 267);
    const auto *snf_268 = buffer.data(snf + 268);
    const auto *snf_269 = buffer.data(snf + 269);
    const auto *snf_270 = buffer.data(snf + 270);
    const auto *snf_272 = buffer.data(snf + 272);
    const auto *snf_273 = buffer.data(snf + 273);
    const auto *snf_275 = buffer.data(snf + 275);
    const auto *snf_276 = buffer.data(snf + 276);
    const auto *snf_277 = buffer.data(snf + 277);
    const auto *snf_278 = buffer.data(snf + 278);
    const auto *snf_279 = buffer.data(snf + 279);
    const auto *snf_280 = buffer.data(snf + 280);
    const auto *snf_282 = buffer.data(snf + 282);
    const auto *snf_283 = buffer.data(snf + 283);
    const auto *snf_285 = buffer.data(snf + 285);
    const auto *snf_286 = buffer.data(snf + 286);
    const auto *snf_287 = buffer.data(snf + 287);
    const auto *snf_288 = buffer.data(snf + 288);
    const auto *snf_289 = buffer.data(snf + 289);
    const auto *snf_290 = buffer.data(snf + 290);
    const auto *snf_292 = buffer.data(snf + 292);
    const auto *snf_295 = buffer.data(snf + 295);
    const auto *snf_296 = buffer.data(snf + 296);
    const auto *snf_297 = buffer.data(snf + 297);
    const auto *snf_298 = buffer.data(snf + 298);
    const auto *snf_299 = buffer.data(snf + 299);
    const auto *snf_300 = buffer.data(snf + 300);
    const auto *snf_302 = buffer.data(snf + 302);
    const auto *snf_303 = buffer.data(snf + 303);
    const auto *snf_305 = buffer.data(snf + 305);
    const auto *snf_306 = buffer.data(snf + 306);
    const auto *snf_307 = buffer.data(snf + 307);
    const auto *snf_308 = buffer.data(snf + 308);
    const auto *snf_309 = buffer.data(snf + 309);
    const auto *snf_310 = buffer.data(snf + 310);
    const auto *snf_312 = buffer.data(snf + 312);
    const auto *snf_313 = buffer.data(snf + 313);
    const auto *snf_315 = buffer.data(snf + 315);
    const auto *snf_316 = buffer.data(snf + 316);
    const auto *snf_317 = buffer.data(snf + 317);
    const auto *snf_318 = buffer.data(snf + 318);
    const auto *snf_319 = buffer.data(snf + 319);

#pragma omp simd aligned(t_361, t_362, t_363, t_364, pc_x, pc_y, pc_z, smf_170, smf_180, \
                         smf_182, smf_243, snd0_147, snd1_147, snf_240, snf_242, \
                         snf_243 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_361[k] = f_12 * smf_180[k]
                   + f_3 * pc_y[k] * snf_240[k];

        t_362[k] = f_12 * smf_170[k]
                   + f_3 * pc_z[k] * snf_240[k];

        t_363[k] = f_14 * smf_243[k]
                   + f_4 * snd0_147[k]
                   - f_5 * snd1_147[k]
                   + f_3 * pc_x[k] * snf_243[k];

        t_364[k] = f_12 * smf_182[k]
                   + f_3 * pc_y[k] * snf_242[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, pc_x, smf_245, smf_246, smf_247, smf_248, \
                         snd0_149, snd1_149, snf_245, snf_246, snf_247, \
                         snf_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = f_14 * smf_245[k]
                   + f_4 * snd0_149[k]
                   - f_5 * snd1_149[k]
                   + f_3 * pc_x[k] * snf_245[k];

        t_366[k] = f_14 * smf_246[k]
                   + f_3 * pc_x[k] * snf_246[k];

        t_367[k] = f_14 * smf_247[k]
                   + f_3 * pc_x[k] * snf_247[k];

        t_368[k] = f_14 * smf_248[k]
                   + f_3 * pc_x[k] * snf_248[k];
    }

#pragma omp simd aligned(t_369, t_370, t_371, pc_x, pc_y, pc_z, smf_176, smf_186, smf_249, \
                         snd0_147, snd1_147, snf_246, snf_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_369[k] = f_14 * smf_249[k]
                   + f_3 * pc_x[k] * snf_249[k];

        t_370[k] = f_12 * smf_186[k]
                   + f_1 * snd0_147[k]
                   - f_2 * snd1_147[k]
                   + f_3 * pc_y[k] * snf_246[k];

        t_371[k] = f_12 * smf_176[k]
                   + f_3 * pc_z[k] * snf_246[k];
    }

#pragma omp simd aligned(t_372, t_373, t_374, pc_y, pc_z, smf_179, smf_188, smf_189, snd0_149, \
                         snd1_149, snf_248, snf_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_372[k] = f_12 * smf_188[k]
                   + f_4 * snd0_149[k]
                   - f_5 * snd1_149[k]
                   + f_3 * pc_y[k] * snf_248[k];

        t_373[k] = f_12 * smf_189[k]
                   + f_3 * pc_y[k] * snf_249[k];

        t_374[k] = f_12 * smf_179[k]
                   + f_1 * snd0_149[k]
                   - f_2 * snd1_149[k]
                   + f_3 * pc_z[k] * snf_249[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, pc_x, pc_y, pc_z, smf_180, smf_190, smf_250, \
                         snd0_150, snd1_150, snf_250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = f_14 * smf_250[k]
                   + f_1 * snd0_150[k]
                   - f_2 * snd1_150[k]
                   + f_3 * pc_x[k] * snf_250[k];

        t_376[k] = f_8 * smf_190[k]
                   + f_3 * pc_y[k] * snf_250[k];

        t_377[k] = f_14 * smf_180[k]
                   + f_3 * pc_z[k] * snf_250[k];
    }

#pragma omp simd aligned(t_378, t_379, t_380, pc_x, pc_y, smf_192, smf_253, smf_255, snd0_153, \
                         snd0_155, snd1_153, snd1_155, snf_252, snf_253, \
                         snf_255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_378[k] = f_14 * smf_253[k]
                   + f_4 * snd0_153[k]
                   - f_5 * snd1_153[k]
                   + f_3 * pc_x[k] * snf_253[k];

        t_379[k] = f_8 * smf_192[k]
                   + f_3 * pc_y[k] * snf_252[k];

        t_380[k] = f_14 * smf_255[k]
                   + f_4 * snd0_155[k]
                   - f_5 * snd1_155[k]
                   + f_3 * pc_x[k] * snf_255[k];
    }

#pragma omp simd aligned(t_381, t_382, t_383, t_384, pc_x, smf_256, smf_257, smf_258, smf_259, \
                         snf_256, snf_257, snf_258, snf_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_381[k] = f_14 * smf_256[k]
                   + f_3 * pc_x[k] * snf_256[k];

        t_382[k] = f_14 * smf_257[k]
                   + f_3 * pc_x[k] * snf_257[k];

        t_383[k] = f_14 * smf_258[k]
                   + f_3 * pc_x[k] * snf_258[k];

        t_384[k] = f_14 * smf_259[k]
                   + f_3 * pc_x[k] * snf_259[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, pc_y, pc_z, smf_186, smf_196, smf_198, snd0_153, \
                         snd0_155, snd1_153, snd1_155, snf_256, \
                         snf_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = f_8 * smf_196[k]
                   + f_1 * snd0_153[k]
                   - f_2 * snd1_153[k]
                   + f_3 * pc_y[k] * snf_256[k];

        t_386[k] = f_14 * smf_186[k]
                   + f_3 * pc_z[k] * snf_256[k];

        t_387[k] = f_8 * smf_198[k]
                   + f_4 * snd0_155[k]
                   - f_5 * snd1_155[k]
                   + f_3 * pc_y[k] * snf_258[k];
    }

#pragma omp simd aligned(t_388, t_389, t_390, t_391, pb_y, pc_y, pc_z, smg0_300, smf_189, \
                         smf_199, smf_200, smg1_300, snd0_155, snd1_155, snf_259, \
                         snf_260 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_388[k] = f_8 * smf_199[k]
                   + f_3 * pc_y[k] * snf_259[k];

        t_389[k] = f_14 * smf_189[k]
                   + f_1 * snd0_155[k]
                   - f_2 * snd1_155[k]
                   + f_3 * pc_z[k] * snf_259[k];

        t_390[k] = pb_y[k] * smg0_300[k]
                   - f_6 * pc_y[k] * smg1_300[k];

        t_391[k] = f_7 * smf_200[k]
                   + f_3 * pc_y[k] * snf_260[k];
    }

#pragma omp simd aligned(t_392, t_393, t_394, t_395, pb_y, pc_y, pc_z, smg0_303, smg0_305, \
                         smf_190, smf_201, smf_202, smg1_303, smg1_305, snf_260, \
                         snf_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_392[k] = f_15 * smf_190[k]
                   + f_3 * pc_z[k] * snf_260[k];

        t_393[k] = pb_y[k] * smg0_303[k]
                   + f_8 * smf_201[k]
                   - f_6 * pc_y[k] * smg1_303[k];

        t_394[k] = f_7 * smf_202[k]
                   + f_3 * pc_y[k] * snf_262[k];

        t_395[k] = pb_y[k] * smg0_305[k]
                   - f_6 * pc_y[k] * smg1_305[k];
    }

#pragma omp simd aligned(t_396, t_397, t_398, t_399, pc_x, smf_266, smf_267, smf_268, smf_269, \
                         snf_266, snf_267, snf_268, snf_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_396[k] = f_14 * smf_266[k]
                   + f_3 * pc_x[k] * snf_266[k];

        t_397[k] = f_14 * smf_267[k]
                   + f_3 * pc_x[k] * snf_267[k];

        t_398[k] = f_14 * smf_268[k]
                   + f_3 * pc_x[k] * snf_268[k];

        t_399[k] = f_14 * smf_269[k]
                   + f_3 * pc_x[k] * snf_269[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, pc_y, pc_z, smf_196, smf_206, smf_208, snd0_159, \
                         snd0_161, snd1_159, snd1_161, snf_266, \
                         snf_268 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = f_7 * smf_206[k]
                   + f_1 * snd0_159[k]
                   - f_2 * snd1_159[k]
                   + f_3 * pc_y[k] * snf_266[k];

        t_401[k] = f_15 * smf_196[k]
                   + f_3 * pc_z[k] * snf_266[k];

        t_402[k] = f_7 * smf_208[k]
                   + f_4 * snd0_161[k]
                   - f_5 * snd1_161[k]
                   + f_3 * pc_y[k] * snf_268[k];
    }

#pragma omp simd aligned(t_403, t_404, t_405, t_406, pb_y, pc_x, pc_y, smg0_314, smf_209, \
                         smf_270, smg1_314, snd0_162, snd1_162, snf_269, \
                         snf_270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_403[k] = f_7 * smf_209[k]
                   + f_3 * pc_y[k] * snf_269[k];

        t_404[k] = pb_y[k] * smg0_314[k]
                   - f_6 * pc_y[k] * smg1_314[k];

        t_405[k] = f_14 * smf_270[k]
                   + f_1 * snd0_162[k]
                   - f_2 * snd1_162[k]
                   + f_3 * pc_x[k] * snf_270[k];

        t_406[k] = f_3 * pc_y[k] * snf_270[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, pc_x, pc_y, pc_z, smf_200, smf_273, snd0_165, \
                         snd1_165, snf_270, snf_272, snf_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = f_13 * smf_200[k]
                   + f_3 * pc_z[k] * snf_270[k];

        t_408[k] = f_14 * smf_273[k]
                   + f_4 * snd0_165[k]
                   - f_5 * snd1_165[k]
                   + f_3 * pc_x[k] * snf_273[k];

        t_409[k] = f_3 * pc_y[k] * snf_272[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, pc_x, smf_275, smf_276, smf_277, smf_278, \
                         snd0_167, snd1_167, snf_275, snf_276, snf_277, \
                         snf_278 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = f_14 * smf_275[k]
                   + f_4 * snd0_167[k]
                   - f_5 * snd1_167[k]
                   + f_3 * pc_x[k] * snf_275[k];

        t_411[k] = f_14 * smf_276[k]
                   + f_3 * pc_x[k] * snf_276[k];

        t_412[k] = f_14 * smf_277[k]
                   + f_3 * pc_x[k] * snf_277[k];

        t_413[k] = f_14 * smf_278[k]
                   + f_3 * pc_x[k] * snf_278[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, t_417, pc_x, pc_y, pc_z, smf_206, smf_279, \
                         snd0_165, snd0_167, snd1_165, snd1_167, snf_276, snf_278, \
                         snf_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_14 * smf_279[k]
                   + f_3 * pc_x[k] * snf_279[k];

        t_415[k] = f_1 * snd0_165[k]
                   - f_2 * snd1_165[k]
                   + f_3 * pc_y[k] * snf_276[k];

        t_416[k] = f_13 * smf_206[k]
                   + f_3 * pc_z[k] * snf_276[k];

        t_417[k] = f_4 * snd0_167[k]
                   - f_5 * snd1_167[k]
                   + f_3 * pc_y[k] * snf_278[k];
    }

#pragma omp simd aligned(t_418, t_419, t_420, t_421, pc_x, pc_y, pc_z, smf_209, smf_210, \
                         smf_280, snd0_167, snd0_168, snd1_167, snd1_168, snf_279, \
                         snf_280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_418[k] = f_3 * pc_y[k] * snf_279[k];

        t_419[k] = f_13 * smf_209[k]
                   + f_1 * snd0_167[k]
                   - f_2 * snd1_167[k]
                   + f_3 * pc_z[k] * snf_279[k];

        t_420[k] = f_12 * smf_280[k]
                   + f_1 * snd0_168[k]
                   - f_2 * snd1_168[k]
                   + f_3 * pc_x[k] * snf_280[k];

        t_421[k] = f_11 * smf_210[k]
                   + f_3 * pc_y[k] * snf_280[k];
    }

#pragma omp simd aligned(t_422, t_423, t_424, pc_x, pc_y, pc_z, smf_212, smf_283, snd0_171, \
                         snd1_171, snf_280, snf_282, snf_283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_422[k] = f_3 * pc_z[k] * snf_280[k];

        t_423[k] = f_12 * smf_283[k]
                   + f_4 * snd0_171[k]
                   - f_5 * snd1_171[k]
                   + f_3 * pc_x[k] * snf_283[k];

        t_424[k] = f_11 * smf_212[k]
                   + f_3 * pc_y[k] * snf_282[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, pc_x, smf_285, smf_286, smf_287, smf_288, \
                         snd0_173, snd1_173, snf_285, snf_286, snf_287, \
                         snf_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = f_12 * smf_285[k]
                   + f_4 * snd0_173[k]
                   - f_5 * snd1_173[k]
                   + f_3 * pc_x[k] * snf_285[k];

        t_426[k] = f_12 * smf_286[k]
                   + f_3 * pc_x[k] * snf_286[k];

        t_427[k] = f_12 * smf_287[k]
                   + f_3 * pc_x[k] * snf_287[k];

        t_428[k] = f_12 * smf_288[k]
                   + f_3 * pc_x[k] * snf_288[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, pc_x, pc_y, pc_z, smf_216, smf_289, snd0_171, \
                         snd1_171, snf_286, snf_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = f_12 * smf_289[k]
                   + f_3 * pc_x[k] * snf_289[k];

        t_430[k] = f_11 * smf_216[k]
                   + f_1 * snd0_171[k]
                   - f_2 * snd1_171[k]
                   + f_3 * pc_y[k] * snf_286[k];

        t_431[k] = f_3 * pc_z[k] * snf_286[k];
    }

#pragma omp simd aligned(t_432, t_433, t_434, t_435, pb_z, pc_y, pc_z, smg0_315, smf_218, \
                         smf_219, smg1_315, snd0_173, snd1_173, snf_288, \
                         snf_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_432[k] = f_11 * smf_218[k]
                   + f_4 * snd0_173[k]
                   - f_5 * snd1_173[k]
                   + f_3 * pc_y[k] * snf_288[k];

        t_433[k] = f_11 * smf_219[k]
                   + f_3 * pc_y[k] * snf_289[k];

        t_434[k] = f_1 * snd0_173[k]
                   - f_2 * snd1_173[k]
                   + f_3 * pc_z[k] * snf_289[k];

        t_435[k] = pb_z[k] * smg0_315[k]
                   - f_6 * pc_z[k] * smg1_315[k];
    }

#pragma omp simd aligned(t_436, t_437, t_438, t_439, pb_z, pc_y, pc_z, smg0_318, smf_210, \
                         smf_220, smf_222, smg1_318, snf_290, snf_292 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_436[k] = f_13 * smf_220[k]
                   + f_3 * pc_y[k] * snf_290[k];

        t_437[k] = f_7 * smf_210[k]
                   + f_3 * pc_z[k] * snf_290[k];

        t_438[k] = pb_z[k] * smg0_318[k]
                   - f_6 * pc_z[k] * smg1_318[k];

        t_439[k] = f_13 * smf_222[k]
                   + f_3 * pc_y[k] * snf_292[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, pc_x, smf_295, smf_296, smf_297, smf_298, \
                         snd0_179, snd1_179, snf_295, snf_296, snf_297, \
                         snf_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = f_12 * smf_295[k]
                   + f_4 * snd0_179[k]
                   - f_5 * snd1_179[k]
                   + f_3 * pc_x[k] * snf_295[k];

        t_441[k] = f_12 * smf_296[k]
                   + f_3 * pc_x[k] * snf_296[k];

        t_442[k] = f_12 * smf_297[k]
                   + f_3 * pc_x[k] * snf_297[k];

        t_443[k] = f_12 * smf_298[k]
                   + f_3 * pc_x[k] * snf_298[k];
    }

#pragma omp simd aligned(t_444, t_445, t_446, pb_z, pc_x, pc_z, smg0_325, smf_216, smf_299, \
                         smg1_325, snf_296, snf_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_444[k] = f_12 * smf_299[k]
                   + f_3 * pc_x[k] * snf_299[k];

        t_445[k] = pb_z[k] * smg0_325[k]
                   - f_6 * pc_z[k] * smg1_325[k];

        t_446[k] = f_7 * smf_216[k]
                   + f_3 * pc_z[k] * snf_296[k];
    }

#pragma omp simd aligned(t_447, t_448, t_449, pc_y, pc_z, smf_219, smf_228, smf_229, snd0_179, \
                         snd1_179, snf_298, snf_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_447[k] = f_13 * smf_228[k]
                   + f_4 * snd0_179[k]
                   - f_5 * snd1_179[k]
                   + f_3 * pc_y[k] * snf_298[k];

        t_448[k] = f_13 * smf_229[k]
                   + f_3 * pc_y[k] * snf_299[k];

        t_449[k] = f_7 * smf_219[k]
                   + f_1 * snd0_179[k]
                   - f_2 * snd1_179[k]
                   + f_3 * pc_z[k] * snf_299[k];
    }

#pragma omp simd aligned(t_450, t_451, t_452, pc_x, pc_y, pc_z, smf_220, smf_230, smf_300, \
                         snd0_180, snd1_180, snf_300 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_450[k] = f_12 * smf_300[k]
                   + f_1 * snd0_180[k]
                   - f_2 * snd1_180[k]
                   + f_3 * pc_x[k] * snf_300[k];

        t_451[k] = f_15 * smf_230[k]
                   + f_3 * pc_y[k] * snf_300[k];

        t_452[k] = f_8 * smf_220[k]
                   + f_3 * pc_z[k] * snf_300[k];
    }

#pragma omp simd aligned(t_453, t_454, t_455, pc_x, pc_y, smf_232, smf_303, smf_305, snd0_183, \
                         snd0_185, snd1_183, snd1_185, snf_302, snf_303, \
                         snf_305 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_453[k] = f_12 * smf_303[k]
                   + f_4 * snd0_183[k]
                   - f_5 * snd1_183[k]
                   + f_3 * pc_x[k] * snf_303[k];

        t_454[k] = f_15 * smf_232[k]
                   + f_3 * pc_y[k] * snf_302[k];

        t_455[k] = f_12 * smf_305[k]
                   + f_4 * snd0_185[k]
                   - f_5 * snd1_185[k]
                   + f_3 * pc_x[k] * snf_305[k];
    }

#pragma omp simd aligned(t_456, t_457, t_458, t_459, pc_x, smf_306, smf_307, smf_308, smf_309, \
                         snf_306, snf_307, snf_308, snf_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_456[k] = f_12 * smf_306[k]
                   + f_3 * pc_x[k] * snf_306[k];

        t_457[k] = f_12 * smf_307[k]
                   + f_3 * pc_x[k] * snf_307[k];

        t_458[k] = f_12 * smf_308[k]
                   + f_3 * pc_x[k] * snf_308[k];

        t_459[k] = f_12 * smf_309[k]
                   + f_3 * pc_x[k] * snf_309[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, pc_y, pc_z, smf_226, smf_236, smf_238, snd0_183, \
                         snd0_185, snd1_183, snd1_185, snf_306, \
                         snf_308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = f_15 * smf_236[k]
                   + f_1 * snd0_183[k]
                   - f_2 * snd1_183[k]
                   + f_3 * pc_y[k] * snf_306[k];

        t_461[k] = f_8 * smf_226[k]
                   + f_3 * pc_z[k] * snf_306[k];

        t_462[k] = f_15 * smf_238[k]
                   + f_4 * snd0_185[k]
                   - f_5 * snd1_185[k]
                   + f_3 * pc_y[k] * snf_308[k];
    }

#pragma omp simd aligned(t_463, t_464, t_465, pc_x, pc_y, pc_z, smf_229, smf_239, smf_310, \
                         snd0_185, snd0_186, snd1_185, snd1_186, snf_309, \
                         snf_310 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_463[k] = f_15 * smf_239[k]
                   + f_3 * pc_y[k] * snf_309[k];

        t_464[k] = f_8 * smf_229[k]
                   + f_1 * snd0_185[k]
                   - f_2 * snd1_185[k]
                   + f_3 * pc_z[k] * snf_309[k];

        t_465[k] = f_12 * smf_310[k]
                   + f_1 * snd0_186[k]
                   - f_2 * snd1_186[k]
                   + f_3 * pc_x[k] * snf_310[k];
    }

#pragma omp simd aligned(t_466, t_467, t_468, t_469, pc_x, pc_y, pc_z, smf_230, smf_240, \
                         smf_242, smf_313, snd0_189, snd1_189, snf_310, snf_312, \
                         snf_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_466[k] = f_14 * smf_240[k]
                   + f_3 * pc_y[k] * snf_310[k];

        t_467[k] = f_12 * smf_230[k]
                   + f_3 * pc_z[k] * snf_310[k];

        t_468[k] = f_12 * smf_313[k]
                   + f_4 * snd0_189[k]
                   - f_5 * snd1_189[k]
                   + f_3 * pc_x[k] * snf_313[k];

        t_469[k] = f_14 * smf_242[k]
                   + f_3 * pc_y[k] * snf_312[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, pc_x, smf_315, smf_316, smf_317, smf_318, \
                         snd0_191, snd1_191, snf_315, snf_316, snf_317, \
                         snf_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = f_12 * smf_315[k]
                   + f_4 * snd0_191[k]
                   - f_5 * snd1_191[k]
                   + f_3 * pc_x[k] * snf_315[k];

        t_471[k] = f_12 * smf_316[k]
                   + f_3 * pc_x[k] * snf_316[k];

        t_472[k] = f_12 * smf_317[k]
                   + f_3 * pc_x[k] * snf_317[k];

        t_473[k] = f_12 * smf_318[k]
                   + f_3 * pc_x[k] * snf_318[k];
    }

#pragma omp simd aligned(t_474, t_475, t_476, pc_x, pc_y, pc_z, smf_236, smf_246, smf_319, \
                         snd0_189, snd1_189, snf_316, snf_319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = f_12 * smf_319[k]
                   + f_3 * pc_x[k] * snf_319[k];

        t_475[k] = f_14 * smf_246[k]
                   + f_1 * snd0_189[k]
                   - f_2 * snd1_189[k]
                   + f_3 * pc_y[k] * snf_316[k];

        t_476[k] = f_12 * smf_236[k]
                   + f_3 * pc_z[k] * snf_316[k];
    }
}

static auto
compute_prim_sng_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t smg0,
                                                          const size_t smf, const size_t smg1,
                                                          const size_t snd0, const size_t snd1,
                                                          const size_t snf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_10 = 4.0 / q;
    const auto f_11 = 3.5 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 3.0 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 2.5 / q;

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
    auto *t_586 = buffer.data(target + 586);
    auto *t_587 = buffer.data(target + 587);
    auto *t_588 = buffer.data(target + 588);
    auto *t_589 = buffer.data(target + 589);
    auto *t_590 = buffer.data(target + 590);
    auto *t_591 = buffer.data(target + 591);
    auto *t_592 = buffer.data(target + 592);
    auto *t_593 = buffer.data(target + 593);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *smg0_405 = buffer.data(smg0 + 405);
    const auto *smg0_408 = buffer.data(smg0 + 408);
    const auto *smg0_410 = buffer.data(smg0 + 410);
    const auto *smg0_419 = buffer.data(smg0 + 419);
    const auto *smg0_420 = buffer.data(smg0 + 420);
    const auto *smg0_423 = buffer.data(smg0 + 423);
    const auto *smg0_430 = buffer.data(smg0 + 430);

    const auto *smf_239 = buffer.data(smf + 239);
    const auto *smf_240 = buffer.data(smf + 240);
    const auto *smf_246 = buffer.data(smf + 246);
    const auto *smf_248 = buffer.data(smf + 248);
    const auto *smf_249 = buffer.data(smf + 249);
    const auto *smf_250 = buffer.data(smf + 250);
    const auto *smf_252 = buffer.data(smf + 252);
    const auto *smf_256 = buffer.data(smf + 256);
    const auto *smf_258 = buffer.data(smf + 258);
    const auto *smf_259 = buffer.data(smf + 259);
    const auto *smf_260 = buffer.data(smf + 260);
    const auto *smf_262 = buffer.data(smf + 262);
    const auto *smf_266 = buffer.data(smf + 266);
    const auto *smf_268 = buffer.data(smf + 268);
    const auto *smf_269 = buffer.data(smf + 269);
    const auto *smf_270 = buffer.data(smf + 270);
    const auto *smf_271 = buffer.data(smf + 271);
    const auto *smf_272 = buffer.data(smf + 272);
    const auto *smf_276 = buffer.data(smf + 276);
    const auto *smf_278 = buffer.data(smf + 278);
    const auto *smf_279 = buffer.data(smf + 279);
    const auto *smf_280 = buffer.data(smf + 280);
    const auto *smf_282 = buffer.data(smf + 282);
    const auto *smf_286 = buffer.data(smf + 286);
    const auto *smf_288 = buffer.data(smf + 288);
    const auto *smf_289 = buffer.data(smf + 289);
    const auto *smf_290 = buffer.data(smf + 290);
    const auto *smf_292 = buffer.data(smf + 292);
    const auto *smf_296 = buffer.data(smf + 296);
    const auto *smf_298 = buffer.data(smf + 298);
    const auto *smf_299 = buffer.data(smf + 299);
    const auto *smf_300 = buffer.data(smf + 300);
    const auto *smf_302 = buffer.data(smf + 302);
    const auto *smf_306 = buffer.data(smf + 306);
    const auto *smf_308 = buffer.data(smf + 308);
    const auto *smf_309 = buffer.data(smf + 309);
    const auto *smf_310 = buffer.data(smf + 310);
    const auto *smf_312 = buffer.data(smf + 312);
    const auto *smf_320 = buffer.data(smf + 320);
    const auto *smf_323 = buffer.data(smf + 323);
    const auto *smf_325 = buffer.data(smf + 325);
    const auto *smf_326 = buffer.data(smf + 326);
    const auto *smf_327 = buffer.data(smf + 327);
    const auto *smf_328 = buffer.data(smf + 328);
    const auto *smf_329 = buffer.data(smf + 329);
    const auto *smf_330 = buffer.data(smf + 330);
    const auto *smf_333 = buffer.data(smf + 333);
    const auto *smf_335 = buffer.data(smf + 335);
    const auto *smf_336 = buffer.data(smf + 336);
    const auto *smf_337 = buffer.data(smf + 337);
    const auto *smf_338 = buffer.data(smf + 338);
    const auto *smf_339 = buffer.data(smf + 339);
    const auto *smf_346 = buffer.data(smf + 346);
    const auto *smf_347 = buffer.data(smf + 347);
    const auto *smf_348 = buffer.data(smf + 348);
    const auto *smf_349 = buffer.data(smf + 349);
    const auto *smf_350 = buffer.data(smf + 350);
    const auto *smf_353 = buffer.data(smf + 353);
    const auto *smf_355 = buffer.data(smf + 355);
    const auto *smf_356 = buffer.data(smf + 356);
    const auto *smf_357 = buffer.data(smf + 357);
    const auto *smf_358 = buffer.data(smf + 358);
    const auto *smf_359 = buffer.data(smf + 359);
    const auto *smf_360 = buffer.data(smf + 360);
    const auto *smf_363 = buffer.data(smf + 363);
    const auto *smf_365 = buffer.data(smf + 365);
    const auto *smf_366 = buffer.data(smf + 366);
    const auto *smf_367 = buffer.data(smf + 367);
    const auto *smf_368 = buffer.data(smf + 368);
    const auto *smf_369 = buffer.data(smf + 369);
    const auto *smf_375 = buffer.data(smf + 375);
    const auto *smf_376 = buffer.data(smf + 376);
    const auto *smf_377 = buffer.data(smf + 377);
    const auto *smf_378 = buffer.data(smf + 378);
    const auto *smf_379 = buffer.data(smf + 379);
    const auto *smf_380 = buffer.data(smf + 380);
    const auto *smf_383 = buffer.data(smf + 383);
    const auto *smf_385 = buffer.data(smf + 385);
    const auto *smf_386 = buffer.data(smf + 386);
    const auto *smf_387 = buffer.data(smf + 387);
    const auto *smf_388 = buffer.data(smf + 388);
    const auto *smf_389 = buffer.data(smf + 389);
    const auto *smf_390 = buffer.data(smf + 390);
    const auto *smf_393 = buffer.data(smf + 393);
    const auto *smf_395 = buffer.data(smf + 395);
    const auto *smf_396 = buffer.data(smf + 396);
    const auto *smf_397 = buffer.data(smf + 397);
    const auto *smf_398 = buffer.data(smf + 398);

    const auto *smg1_405 = buffer.data(smg1 + 405);
    const auto *smg1_408 = buffer.data(smg1 + 408);
    const auto *smg1_410 = buffer.data(smg1 + 410);
    const auto *smg1_419 = buffer.data(smg1 + 419);
    const auto *smg1_420 = buffer.data(smg1 + 420);
    const auto *smg1_423 = buffer.data(smg1 + 423);
    const auto *smg1_430 = buffer.data(smg1 + 430);

    const auto *snd0_191 = buffer.data(snd0 + 191);
    const auto *snd0_192 = buffer.data(snd0 + 192);
    const auto *snd0_195 = buffer.data(snd0 + 195);
    const auto *snd0_197 = buffer.data(snd0 + 197);
    const auto *snd0_198 = buffer.data(snd0 + 198);
    const auto *snd0_201 = buffer.data(snd0 + 201);
    const auto *snd0_203 = buffer.data(snd0 + 203);
    const auto *snd0_207 = buffer.data(snd0 + 207);
    const auto *snd0_209 = buffer.data(snd0 + 209);
    const auto *snd0_210 = buffer.data(snd0 + 210);
    const auto *snd0_213 = buffer.data(snd0 + 213);
    const auto *snd0_215 = buffer.data(snd0 + 215);
    const auto *snd0_216 = buffer.data(snd0 + 216);
    const auto *snd0_219 = buffer.data(snd0 + 219);
    const auto *snd0_221 = buffer.data(snd0 + 221);
    const auto *snd0_227 = buffer.data(snd0 + 227);
    const auto *snd0_228 = buffer.data(snd0 + 228);
    const auto *snd0_231 = buffer.data(snd0 + 231);
    const auto *snd0_233 = buffer.data(snd0 + 233);
    const auto *snd0_234 = buffer.data(snd0 + 234);
    const auto *snd0_237 = buffer.data(snd0 + 237);
    const auto *snd0_239 = buffer.data(snd0 + 239);

    const auto *snd1_191 = buffer.data(snd1 + 191);
    const auto *snd1_192 = buffer.data(snd1 + 192);
    const auto *snd1_195 = buffer.data(snd1 + 195);
    const auto *snd1_197 = buffer.data(snd1 + 197);
    const auto *snd1_198 = buffer.data(snd1 + 198);
    const auto *snd1_201 = buffer.data(snd1 + 201);
    const auto *snd1_203 = buffer.data(snd1 + 203);
    const auto *snd1_207 = buffer.data(snd1 + 207);
    const auto *snd1_209 = buffer.data(snd1 + 209);
    const auto *snd1_210 = buffer.data(snd1 + 210);
    const auto *snd1_213 = buffer.data(snd1 + 213);
    const auto *snd1_215 = buffer.data(snd1 + 215);
    const auto *snd1_216 = buffer.data(snd1 + 216);
    const auto *snd1_219 = buffer.data(snd1 + 219);
    const auto *snd1_221 = buffer.data(snd1 + 221);
    const auto *snd1_227 = buffer.data(snd1 + 227);
    const auto *snd1_228 = buffer.data(snd1 + 228);
    const auto *snd1_231 = buffer.data(snd1 + 231);
    const auto *snd1_233 = buffer.data(snd1 + 233);
    const auto *snd1_234 = buffer.data(snd1 + 234);
    const auto *snd1_237 = buffer.data(snd1 + 237);
    const auto *snd1_239 = buffer.data(snd1 + 239);

    const auto *snf_318 = buffer.data(snf + 318);
    const auto *snf_319 = buffer.data(snf + 319);
    const auto *snf_320 = buffer.data(snf + 320);
    const auto *snf_322 = buffer.data(snf + 322);
    const auto *snf_323 = buffer.data(snf + 323);
    const auto *snf_325 = buffer.data(snf + 325);
    const auto *snf_326 = buffer.data(snf + 326);
    const auto *snf_327 = buffer.data(snf + 327);
    const auto *snf_328 = buffer.data(snf + 328);
    const auto *snf_329 = buffer.data(snf + 329);
    const auto *snf_330 = buffer.data(snf + 330);
    const auto *snf_332 = buffer.data(snf + 332);
    const auto *snf_333 = buffer.data(snf + 333);
    const auto *snf_335 = buffer.data(snf + 335);
    const auto *snf_336 = buffer.data(snf + 336);
    const auto *snf_337 = buffer.data(snf + 337);
    const auto *snf_338 = buffer.data(snf + 338);
    const auto *snf_339 = buffer.data(snf + 339);
    const auto *snf_340 = buffer.data(snf + 340);
    const auto *snf_342 = buffer.data(snf + 342);
    const auto *snf_346 = buffer.data(snf + 346);
    const auto *snf_347 = buffer.data(snf + 347);
    const auto *snf_348 = buffer.data(snf + 348);
    const auto *snf_349 = buffer.data(snf + 349);
    const auto *snf_350 = buffer.data(snf + 350);
    const auto *snf_352 = buffer.data(snf + 352);
    const auto *snf_353 = buffer.data(snf + 353);
    const auto *snf_355 = buffer.data(snf + 355);
    const auto *snf_356 = buffer.data(snf + 356);
    const auto *snf_357 = buffer.data(snf + 357);
    const auto *snf_358 = buffer.data(snf + 358);
    const auto *snf_359 = buffer.data(snf + 359);
    const auto *snf_360 = buffer.data(snf + 360);
    const auto *snf_362 = buffer.data(snf + 362);
    const auto *snf_363 = buffer.data(snf + 363);
    const auto *snf_365 = buffer.data(snf + 365);
    const auto *snf_366 = buffer.data(snf + 366);
    const auto *snf_367 = buffer.data(snf + 367);
    const auto *snf_368 = buffer.data(snf + 368);
    const auto *snf_369 = buffer.data(snf + 369);
    const auto *snf_370 = buffer.data(snf + 370);
    const auto *snf_372 = buffer.data(snf + 372);
    const auto *snf_375 = buffer.data(snf + 375);
    const auto *snf_376 = buffer.data(snf + 376);
    const auto *snf_377 = buffer.data(snf + 377);
    const auto *snf_378 = buffer.data(snf + 378);
    const auto *snf_379 = buffer.data(snf + 379);
    const auto *snf_380 = buffer.data(snf + 380);
    const auto *snf_382 = buffer.data(snf + 382);
    const auto *snf_383 = buffer.data(snf + 383);
    const auto *snf_385 = buffer.data(snf + 385);
    const auto *snf_386 = buffer.data(snf + 386);
    const auto *snf_387 = buffer.data(snf + 387);
    const auto *snf_388 = buffer.data(snf + 388);
    const auto *snf_389 = buffer.data(snf + 389);
    const auto *snf_390 = buffer.data(snf + 390);
    const auto *snf_392 = buffer.data(snf + 392);
    const auto *snf_393 = buffer.data(snf + 393);
    const auto *snf_395 = buffer.data(snf + 395);
    const auto *snf_396 = buffer.data(snf + 396);
    const auto *snf_397 = buffer.data(snf + 397);
    const auto *snf_398 = buffer.data(snf + 398);

#pragma omp simd aligned(t_477, t_478, t_479, pc_y, pc_z, smf_239, smf_248, smf_249, snd0_191, \
                         snd1_191, snf_318, snf_319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_477[k] = f_14 * smf_248[k]
                   + f_4 * snd0_191[k]
                   - f_5 * snd1_191[k]
                   + f_3 * pc_y[k] * snf_318[k];

        t_478[k] = f_14 * smf_249[k]
                   + f_3 * pc_y[k] * snf_319[k];

        t_479[k] = f_12 * smf_239[k]
                   + f_1 * snd0_191[k]
                   - f_2 * snd1_191[k]
                   + f_3 * pc_z[k] * snf_319[k];
    }

#pragma omp simd aligned(t_480, t_481, t_482, pc_x, pc_y, pc_z, smf_240, smf_250, smf_320, \
                         snd0_192, snd1_192, snf_320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_480[k] = f_12 * smf_320[k]
                   + f_1 * snd0_192[k]
                   - f_2 * snd1_192[k]
                   + f_3 * pc_x[k] * snf_320[k];

        t_481[k] = f_12 * smf_250[k]
                   + f_3 * pc_y[k] * snf_320[k];

        t_482[k] = f_14 * smf_240[k]
                   + f_3 * pc_z[k] * snf_320[k];
    }

#pragma omp simd aligned(t_483, t_484, t_485, pc_x, pc_y, smf_252, smf_323, smf_325, snd0_195, \
                         snd0_197, snd1_195, snd1_197, snf_322, snf_323, \
                         snf_325 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_483[k] = f_12 * smf_323[k]
                   + f_4 * snd0_195[k]
                   - f_5 * snd1_195[k]
                   + f_3 * pc_x[k] * snf_323[k];

        t_484[k] = f_12 * smf_252[k]
                   + f_3 * pc_y[k] * snf_322[k];

        t_485[k] = f_12 * smf_325[k]
                   + f_4 * snd0_197[k]
                   - f_5 * snd1_197[k]
                   + f_3 * pc_x[k] * snf_325[k];
    }

#pragma omp simd aligned(t_486, t_487, t_488, t_489, pc_x, smf_326, smf_327, smf_328, smf_329, \
                         snf_326, snf_327, snf_328, snf_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_486[k] = f_12 * smf_326[k]
                   + f_3 * pc_x[k] * snf_326[k];

        t_487[k] = f_12 * smf_327[k]
                   + f_3 * pc_x[k] * snf_327[k];

        t_488[k] = f_12 * smf_328[k]
                   + f_3 * pc_x[k] * snf_328[k];

        t_489[k] = f_12 * smf_329[k]
                   + f_3 * pc_x[k] * snf_329[k];
    }

#pragma omp simd aligned(t_490, t_491, t_492, pc_y, pc_z, smf_246, smf_256, smf_258, snd0_195, \
                         snd0_197, snd1_195, snd1_197, snf_326, \
                         snf_328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_490[k] = f_12 * smf_256[k]
                   + f_1 * snd0_195[k]
                   - f_2 * snd1_195[k]
                   + f_3 * pc_y[k] * snf_326[k];

        t_491[k] = f_14 * smf_246[k]
                   + f_3 * pc_z[k] * snf_326[k];

        t_492[k] = f_12 * smf_258[k]
                   + f_4 * snd0_197[k]
                   - f_5 * snd1_197[k]
                   + f_3 * pc_y[k] * snf_328[k];
    }

#pragma omp simd aligned(t_493, t_494, t_495, pc_x, pc_y, pc_z, smf_249, smf_259, smf_330, \
                         snd0_197, snd0_198, snd1_197, snd1_198, snf_329, \
                         snf_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_493[k] = f_12 * smf_259[k]
                   + f_3 * pc_y[k] * snf_329[k];

        t_494[k] = f_14 * smf_249[k]
                   + f_1 * snd0_197[k]
                   - f_2 * snd1_197[k]
                   + f_3 * pc_z[k] * snf_329[k];

        t_495[k] = f_12 * smf_330[k]
                   + f_1 * snd0_198[k]
                   - f_2 * snd1_198[k]
                   + f_3 * pc_x[k] * snf_330[k];
    }

#pragma omp simd aligned(t_496, t_497, t_498, t_499, pc_x, pc_y, pc_z, smf_250, smf_260, \
                         smf_262, smf_333, snd0_201, snd1_201, snf_330, snf_332, \
                         snf_333 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_496[k] = f_8 * smf_260[k]
                   + f_3 * pc_y[k] * snf_330[k];

        t_497[k] = f_15 * smf_250[k]
                   + f_3 * pc_z[k] * snf_330[k];

        t_498[k] = f_12 * smf_333[k]
                   + f_4 * snd0_201[k]
                   - f_5 * snd1_201[k]
                   + f_3 * pc_x[k] * snf_333[k];

        t_499[k] = f_8 * smf_262[k]
                   + f_3 * pc_y[k] * snf_332[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, pc_x, smf_335, smf_336, smf_337, smf_338, \
                         snd0_203, snd1_203, snf_335, snf_336, snf_337, \
                         snf_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = f_12 * smf_335[k]
                   + f_4 * snd0_203[k]
                   - f_5 * snd1_203[k]
                   + f_3 * pc_x[k] * snf_335[k];

        t_501[k] = f_12 * smf_336[k]
                   + f_3 * pc_x[k] * snf_336[k];

        t_502[k] = f_12 * smf_337[k]
                   + f_3 * pc_x[k] * snf_337[k];

        t_503[k] = f_12 * smf_338[k]
                   + f_3 * pc_x[k] * snf_338[k];
    }

#pragma omp simd aligned(t_504, t_505, t_506, pc_x, pc_y, pc_z, smf_256, smf_266, smf_339, \
                         snd0_201, snd1_201, snf_336, snf_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_504[k] = f_12 * smf_339[k]
                   + f_3 * pc_x[k] * snf_339[k];

        t_505[k] = f_8 * smf_266[k]
                   + f_1 * snd0_201[k]
                   - f_2 * snd1_201[k]
                   + f_3 * pc_y[k] * snf_336[k];

        t_506[k] = f_15 * smf_256[k]
                   + f_3 * pc_z[k] * snf_336[k];
    }

#pragma omp simd aligned(t_507, t_508, t_509, t_510, pb_y, pc_y, pc_z, smg0_405, smf_259, \
                         smf_268, smf_269, smg1_405, snd0_203, snd1_203, snf_338, \
                         snf_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_507[k] = f_8 * smf_268[k]
                   + f_4 * snd0_203[k]
                   - f_5 * snd1_203[k]
                   + f_3 * pc_y[k] * snf_338[k];

        t_508[k] = f_8 * smf_269[k]
                   + f_3 * pc_y[k] * snf_339[k];

        t_509[k] = f_15 * smf_259[k]
                   + f_1 * snd0_203[k]
                   - f_2 * snd1_203[k]
                   + f_3 * pc_z[k] * snf_339[k];

        t_510[k] = pb_y[k] * smg0_405[k]
                   - f_6 * pc_y[k] * smg1_405[k];
    }

#pragma omp simd aligned(t_511, t_512, t_513, t_514, pb_y, pc_y, pc_z, smg0_408, smf_260, \
                         smf_270, smf_271, smf_272, smg1_408, snf_340, \
                         snf_342 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_511[k] = f_7 * smf_270[k]
                   + f_3 * pc_y[k] * snf_340[k];

        t_512[k] = f_13 * smf_260[k]
                   + f_3 * pc_z[k] * snf_340[k];

        t_513[k] = pb_y[k] * smg0_408[k]
                   + f_8 * smf_271[k]
                   - f_6 * pc_y[k] * smg1_408[k];

        t_514[k] = f_7 * smf_272[k]
                   + f_3 * pc_y[k] * snf_342[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, t_518, pb_y, pc_x, pc_y, smg0_410, smf_346, \
                         smf_347, smf_348, smg1_410, snf_346, snf_347, \
                         snf_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = pb_y[k] * smg0_410[k]
                   - f_6 * pc_y[k] * smg1_410[k];

        t_516[k] = f_12 * smf_346[k]
                   + f_3 * pc_x[k] * snf_346[k];

        t_517[k] = f_12 * smf_347[k]
                   + f_3 * pc_x[k] * snf_347[k];

        t_518[k] = f_12 * smf_348[k]
                   + f_3 * pc_x[k] * snf_348[k];
    }

#pragma omp simd aligned(t_519, t_520, t_521, pc_x, pc_y, pc_z, smf_266, smf_276, smf_349, \
                         snd0_207, snd1_207, snf_346, snf_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_519[k] = f_12 * smf_349[k]
                   + f_3 * pc_x[k] * snf_349[k];

        t_520[k] = f_7 * smf_276[k]
                   + f_1 * snd0_207[k]
                   - f_2 * snd1_207[k]
                   + f_3 * pc_y[k] * snf_346[k];

        t_521[k] = f_13 * smf_266[k]
                   + f_3 * pc_z[k] * snf_346[k];
    }

#pragma omp simd aligned(t_522, t_523, t_524, pb_y, pc_y, smg0_419, smf_278, smf_279, \
                         smg1_419, snd0_209, snd1_209, snf_348, \
                         snf_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_522[k] = f_7 * smf_278[k]
                   + f_4 * snd0_209[k]
                   - f_5 * snd1_209[k]
                   + f_3 * pc_y[k] * snf_348[k];

        t_523[k] = f_7 * smf_279[k]
                   + f_3 * pc_y[k] * snf_349[k];

        t_524[k] = pb_y[k] * smg0_419[k]
                   - f_6 * pc_y[k] * smg1_419[k];
    }

#pragma omp simd aligned(t_525, t_526, t_527, t_528, pc_x, pc_y, pc_z, smf_270, smf_350, \
                         smf_353, snd0_210, snd0_213, snd1_210, snd1_213, snf_350, \
                         snf_353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_525[k] = f_12 * smf_350[k]
                   + f_1 * snd0_210[k]
                   - f_2 * snd1_210[k]
                   + f_3 * pc_x[k] * snf_350[k];

        t_526[k] = f_3 * pc_y[k] * snf_350[k];

        t_527[k] = f_11 * smf_270[k]
                   + f_3 * pc_z[k] * snf_350[k];

        t_528[k] = f_12 * smf_353[k]
                   + f_4 * snd0_213[k]
                   - f_5 * snd1_213[k]
                   + f_3 * pc_x[k] * snf_353[k];
    }

#pragma omp simd aligned(t_529, t_530, t_531, t_532, pc_x, pc_y, smf_355, smf_356, smf_357, \
                         snd0_215, snd1_215, snf_352, snf_355, snf_356, \
                         snf_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_529[k] = f_3 * pc_y[k] * snf_352[k];

        t_530[k] = f_12 * smf_355[k]
                   + f_4 * snd0_215[k]
                   - f_5 * snd1_215[k]
                   + f_3 * pc_x[k] * snf_355[k];

        t_531[k] = f_12 * smf_356[k]
                   + f_3 * pc_x[k] * snf_356[k];

        t_532[k] = f_12 * smf_357[k]
                   + f_3 * pc_x[k] * snf_357[k];
    }

#pragma omp simd aligned(t_533, t_534, t_535, t_536, pc_x, pc_y, pc_z, smf_276, smf_358, \
                         smf_359, snd0_213, snd1_213, snf_356, snf_358, \
                         snf_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_533[k] = f_12 * smf_358[k]
                   + f_3 * pc_x[k] * snf_358[k];

        t_534[k] = f_12 * smf_359[k]
                   + f_3 * pc_x[k] * snf_359[k];

        t_535[k] = f_1 * snd0_213[k]
                   - f_2 * snd1_213[k]
                   + f_3 * pc_y[k] * snf_356[k];

        t_536[k] = f_11 * smf_276[k]
                   + f_3 * pc_z[k] * snf_356[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, t_540, pc_x, pc_y, pc_z, smf_279, smf_360, \
                         snd0_215, snd0_216, snd1_215, snd1_216, snf_358, snf_359, \
                         snf_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = f_4 * snd0_215[k]
                   - f_5 * snd1_215[k]
                   + f_3 * pc_y[k] * snf_358[k];

        t_538[k] = f_3 * pc_y[k] * snf_359[k];

        t_539[k] = f_11 * smf_279[k]
                   + f_1 * snd0_215[k]
                   - f_2 * snd1_215[k]
                   + f_3 * pc_z[k] * snf_359[k];

        t_540[k] = f_8 * smf_360[k]
                   + f_1 * snd0_216[k]
                   - f_2 * snd1_216[k]
                   + f_3 * pc_x[k] * snf_360[k];
    }

#pragma omp simd aligned(t_541, t_542, t_543, t_544, pc_x, pc_y, pc_z, smf_280, smf_282, \
                         smf_363, snd0_219, snd1_219, snf_360, snf_362, \
                         snf_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_541[k] = f_10 * smf_280[k]
                   + f_3 * pc_y[k] * snf_360[k];

        t_542[k] = f_3 * pc_z[k] * snf_360[k];

        t_543[k] = f_8 * smf_363[k]
                   + f_4 * snd0_219[k]
                   - f_5 * snd1_219[k]
                   + f_3 * pc_x[k] * snf_363[k];

        t_544[k] = f_10 * smf_282[k]
                   + f_3 * pc_y[k] * snf_362[k];
    }

#pragma omp simd aligned(t_545, t_546, t_547, t_548, pc_x, smf_365, smf_366, smf_367, smf_368, \
                         snd0_221, snd1_221, snf_365, snf_366, snf_367, \
                         snf_368 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_545[k] = f_8 * smf_365[k]
                   + f_4 * snd0_221[k]
                   - f_5 * snd1_221[k]
                   + f_3 * pc_x[k] * snf_365[k];

        t_546[k] = f_8 * smf_366[k]
                   + f_3 * pc_x[k] * snf_366[k];

        t_547[k] = f_8 * smf_367[k]
                   + f_3 * pc_x[k] * snf_367[k];

        t_548[k] = f_8 * smf_368[k]
                   + f_3 * pc_x[k] * snf_368[k];
    }

#pragma omp simd aligned(t_549, t_550, t_551, pc_x, pc_y, pc_z, smf_286, smf_369, snd0_219, \
                         snd1_219, snf_366, snf_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_549[k] = f_8 * smf_369[k]
                   + f_3 * pc_x[k] * snf_369[k];

        t_550[k] = f_10 * smf_286[k]
                   + f_1 * snd0_219[k]
                   - f_2 * snd1_219[k]
                   + f_3 * pc_y[k] * snf_366[k];

        t_551[k] = f_3 * pc_z[k] * snf_366[k];
    }

#pragma omp simd aligned(t_552, t_553, t_554, t_555, pb_z, pc_y, pc_z, smg0_420, smf_288, \
                         smf_289, smg1_420, snd0_221, snd1_221, snf_368, \
                         snf_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_552[k] = f_10 * smf_288[k]
                   + f_4 * snd0_221[k]
                   - f_5 * snd1_221[k]
                   + f_3 * pc_y[k] * snf_368[k];

        t_553[k] = f_10 * smf_289[k]
                   + f_3 * pc_y[k] * snf_369[k];

        t_554[k] = f_1 * snd0_221[k]
                   - f_2 * snd1_221[k]
                   + f_3 * pc_z[k] * snf_369[k];

        t_555[k] = pb_z[k] * smg0_420[k]
                   - f_6 * pc_z[k] * smg1_420[k];
    }

#pragma omp simd aligned(t_556, t_557, t_558, t_559, pb_z, pc_y, pc_z, smg0_423, smf_280, \
                         smf_290, smf_292, smg1_423, snf_370, snf_372 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_556[k] = f_11 * smf_290[k]
                   + f_3 * pc_y[k] * snf_370[k];

        t_557[k] = f_7 * smf_280[k]
                   + f_3 * pc_z[k] * snf_370[k];

        t_558[k] = pb_z[k] * smg0_423[k]
                   - f_6 * pc_z[k] * smg1_423[k];

        t_559[k] = f_11 * smf_292[k]
                   + f_3 * pc_y[k] * snf_372[k];
    }

#pragma omp simd aligned(t_560, t_561, t_562, t_563, pc_x, smf_375, smf_376, smf_377, smf_378, \
                         snd0_227, snd1_227, snf_375, snf_376, snf_377, \
                         snf_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_560[k] = f_8 * smf_375[k]
                   + f_4 * snd0_227[k]
                   - f_5 * snd1_227[k]
                   + f_3 * pc_x[k] * snf_375[k];

        t_561[k] = f_8 * smf_376[k]
                   + f_3 * pc_x[k] * snf_376[k];

        t_562[k] = f_8 * smf_377[k]
                   + f_3 * pc_x[k] * snf_377[k];

        t_563[k] = f_8 * smf_378[k]
                   + f_3 * pc_x[k] * snf_378[k];
    }

#pragma omp simd aligned(t_564, t_565, t_566, pb_z, pc_x, pc_z, smg0_430, smf_286, smf_379, \
                         smg1_430, snf_376, snf_379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_564[k] = f_8 * smf_379[k]
                   + f_3 * pc_x[k] * snf_379[k];

        t_565[k] = pb_z[k] * smg0_430[k]
                   - f_6 * pc_z[k] * smg1_430[k];

        t_566[k] = f_7 * smf_286[k]
                   + f_3 * pc_z[k] * snf_376[k];
    }

#pragma omp simd aligned(t_567, t_568, t_569, pc_y, pc_z, smf_289, smf_298, smf_299, snd0_227, \
                         snd1_227, snf_378, snf_379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_567[k] = f_11 * smf_298[k]
                   + f_4 * snd0_227[k]
                   - f_5 * snd1_227[k]
                   + f_3 * pc_y[k] * snf_378[k];

        t_568[k] = f_11 * smf_299[k]
                   + f_3 * pc_y[k] * snf_379[k];

        t_569[k] = f_7 * smf_289[k]
                   + f_1 * snd0_227[k]
                   - f_2 * snd1_227[k]
                   + f_3 * pc_z[k] * snf_379[k];
    }

#pragma omp simd aligned(t_570, t_571, t_572, pc_x, pc_y, pc_z, smf_290, smf_300, smf_380, \
                         snd0_228, snd1_228, snf_380 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_570[k] = f_8 * smf_380[k]
                   + f_1 * snd0_228[k]
                   - f_2 * snd1_228[k]
                   + f_3 * pc_x[k] * snf_380[k];

        t_571[k] = f_13 * smf_300[k]
                   + f_3 * pc_y[k] * snf_380[k];

        t_572[k] = f_8 * smf_290[k]
                   + f_3 * pc_z[k] * snf_380[k];
    }

#pragma omp simd aligned(t_573, t_574, t_575, pc_x, pc_y, smf_302, smf_383, smf_385, snd0_231, \
                         snd0_233, snd1_231, snd1_233, snf_382, snf_383, \
                         snf_385 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_573[k] = f_8 * smf_383[k]
                   + f_4 * snd0_231[k]
                   - f_5 * snd1_231[k]
                   + f_3 * pc_x[k] * snf_383[k];

        t_574[k] = f_13 * smf_302[k]
                   + f_3 * pc_y[k] * snf_382[k];

        t_575[k] = f_8 * smf_385[k]
                   + f_4 * snd0_233[k]
                   - f_5 * snd1_233[k]
                   + f_3 * pc_x[k] * snf_385[k];
    }

#pragma omp simd aligned(t_576, t_577, t_578, t_579, pc_x, smf_386, smf_387, smf_388, smf_389, \
                         snf_386, snf_387, snf_388, snf_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_576[k] = f_8 * smf_386[k]
                   + f_3 * pc_x[k] * snf_386[k];

        t_577[k] = f_8 * smf_387[k]
                   + f_3 * pc_x[k] * snf_387[k];

        t_578[k] = f_8 * smf_388[k]
                   + f_3 * pc_x[k] * snf_388[k];

        t_579[k] = f_8 * smf_389[k]
                   + f_3 * pc_x[k] * snf_389[k];
    }

#pragma omp simd aligned(t_580, t_581, t_582, pc_y, pc_z, smf_296, smf_306, smf_308, snd0_231, \
                         snd0_233, snd1_231, snd1_233, snf_386, \
                         snf_388 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_580[k] = f_13 * smf_306[k]
                   + f_1 * snd0_231[k]
                   - f_2 * snd1_231[k]
                   + f_3 * pc_y[k] * snf_386[k];

        t_581[k] = f_8 * smf_296[k]
                   + f_3 * pc_z[k] * snf_386[k];

        t_582[k] = f_13 * smf_308[k]
                   + f_4 * snd0_233[k]
                   - f_5 * snd1_233[k]
                   + f_3 * pc_y[k] * snf_388[k];
    }

#pragma omp simd aligned(t_583, t_584, t_585, pc_x, pc_y, pc_z, smf_299, smf_309, smf_390, \
                         snd0_233, snd0_234, snd1_233, snd1_234, snf_389, \
                         snf_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_583[k] = f_13 * smf_309[k]
                   + f_3 * pc_y[k] * snf_389[k];

        t_584[k] = f_8 * smf_299[k]
                   + f_1 * snd0_233[k]
                   - f_2 * snd1_233[k]
                   + f_3 * pc_z[k] * snf_389[k];

        t_585[k] = f_8 * smf_390[k]
                   + f_1 * snd0_234[k]
                   - f_2 * snd1_234[k]
                   + f_3 * pc_x[k] * snf_390[k];
    }

#pragma omp simd aligned(t_586, t_587, t_588, t_589, pc_x, pc_y, pc_z, smf_300, smf_310, \
                         smf_312, smf_393, snd0_237, snd1_237, snf_390, snf_392, \
                         snf_393 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_586[k] = f_15 * smf_310[k]
                   + f_3 * pc_y[k] * snf_390[k];

        t_587[k] = f_12 * smf_300[k]
                   + f_3 * pc_z[k] * snf_390[k];

        t_588[k] = f_8 * smf_393[k]
                   + f_4 * snd0_237[k]
                   - f_5 * snd1_237[k]
                   + f_3 * pc_x[k] * snf_393[k];

        t_589[k] = f_15 * smf_312[k]
                   + f_3 * pc_y[k] * snf_392[k];
    }

#pragma omp simd aligned(t_590, t_591, t_592, t_593, pc_x, smf_395, smf_396, smf_397, smf_398, \
                         snd0_239, snd1_239, snf_395, snf_396, snf_397, \
                         snf_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_590[k] = f_8 * smf_395[k]
                   + f_4 * snd0_239[k]
                   - f_5 * snd1_239[k]
                   + f_3 * pc_x[k] * snf_395[k];

        t_591[k] = f_8 * smf_396[k]
                   + f_3 * pc_x[k] * snf_396[k];

        t_592[k] = f_8 * smf_397[k]
                   + f_3 * pc_x[k] * snf_397[k];

        t_593[k] = f_8 * smf_398[k]
                   + f_3 * pc_x[k] * snf_398[k];
    }
}

static auto
compute_prim_sng_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t smg0,
                                                          const size_t smf, const size_t smg1,
                                                          const size_t snd0, const size_t snd1,
                                                          const size_t snf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 4.5 / q;
    const auto f_10 = 4.0 / q;
    const auto f_11 = 3.5 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 3.0 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 2.5 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *smg0_525 = buffer.data(smg0 + 525);
    const auto *smg0_528 = buffer.data(smg0 + 528);
    const auto *smg0_530 = buffer.data(smg0 + 530);
    const auto *smg0_539 = buffer.data(smg0 + 539);
    const auto *smg0_540 = buffer.data(smg0 + 540);
    const auto *smg0_543 = buffer.data(smg0 + 543);
    const auto *smg0_675 = buffer.data(smg0 + 675);
    const auto *smg0_678 = buffer.data(smg0 + 678);
    const auto *smg0_680 = buffer.data(smg0 + 680);
    const auto *smg0_685 = buffer.data(smg0 + 685);
    const auto *smg0_687 = buffer.data(smg0 + 687);
    const auto *smg0_689 = buffer.data(smg0 + 689);
    const auto *smg0_695 = buffer.data(smg0 + 695);
    const auto *smg0_700 = buffer.data(smg0 + 700);
    const auto *smg0_702 = buffer.data(smg0 + 702);
    const auto *smg0_704 = buffer.data(smg0 + 704);
    const auto *smg0_705 = buffer.data(smg0 + 705);
    const auto *smg0_708 = buffer.data(smg0 + 708);

    const auto *smf_306 = buffer.data(smf + 306);
    const auto *smf_309 = buffer.data(smf + 309);
    const auto *smf_310 = buffer.data(smf + 310);
    const auto *smf_316 = buffer.data(smf + 316);
    const auto *smf_318 = buffer.data(smf + 318);
    const auto *smf_319 = buffer.data(smf + 319);
    const auto *smf_320 = buffer.data(smf + 320);
    const auto *smf_322 = buffer.data(smf + 322);
    const auto *smf_326 = buffer.data(smf + 326);
    const auto *smf_328 = buffer.data(smf + 328);
    const auto *smf_329 = buffer.data(smf + 329);
    const auto *smf_330 = buffer.data(smf + 330);
    const auto *smf_332 = buffer.data(smf + 332);
    const auto *smf_336 = buffer.data(smf + 336);
    const auto *smf_338 = buffer.data(smf + 338);
    const auto *smf_339 = buffer.data(smf + 339);
    const auto *smf_340 = buffer.data(smf + 340);
    const auto *smf_342 = buffer.data(smf + 342);
    const auto *smf_346 = buffer.data(smf + 346);
    const auto *smf_348 = buffer.data(smf + 348);
    const auto *smf_349 = buffer.data(smf + 349);
    const auto *smf_350 = buffer.data(smf + 350);
    const auto *smf_351 = buffer.data(smf + 351);
    const auto *smf_352 = buffer.data(smf + 352);
    const auto *smf_356 = buffer.data(smf + 356);
    const auto *smf_358 = buffer.data(smf + 358);
    const auto *smf_359 = buffer.data(smf + 359);
    const auto *smf_360 = buffer.data(smf + 360);
    const auto *smf_362 = buffer.data(smf + 362);
    const auto *smf_366 = buffer.data(smf + 366);
    const auto *smf_369 = buffer.data(smf + 369);
    const auto *smf_370 = buffer.data(smf + 370);
    const auto *smf_372 = buffer.data(smf + 372);
    const auto *smf_379 = buffer.data(smf + 379);
    const auto *smf_380 = buffer.data(smf + 380);
    const auto *smf_382 = buffer.data(smf + 382);
    const auto *smf_399 = buffer.data(smf + 399);
    const auto *smf_400 = buffer.data(smf + 400);
    const auto *smf_403 = buffer.data(smf + 403);
    const auto *smf_405 = buffer.data(smf + 405);
    const auto *smf_406 = buffer.data(smf + 406);
    const auto *smf_407 = buffer.data(smf + 407);
    const auto *smf_408 = buffer.data(smf + 408);
    const auto *smf_409 = buffer.data(smf + 409);
    const auto *smf_410 = buffer.data(smf + 410);
    const auto *smf_413 = buffer.data(smf + 413);
    const auto *smf_415 = buffer.data(smf + 415);
    const auto *smf_416 = buffer.data(smf + 416);
    const auto *smf_417 = buffer.data(smf + 417);
    const auto *smf_418 = buffer.data(smf + 418);
    const auto *smf_419 = buffer.data(smf + 419);
    const auto *smf_420 = buffer.data(smf + 420);
    const auto *smf_423 = buffer.data(smf + 423);
    const auto *smf_425 = buffer.data(smf + 425);
    const auto *smf_426 = buffer.data(smf + 426);
    const auto *smf_427 = buffer.data(smf + 427);
    const auto *smf_428 = buffer.data(smf + 428);
    const auto *smf_429 = buffer.data(smf + 429);
    const auto *smf_436 = buffer.data(smf + 436);
    const auto *smf_437 = buffer.data(smf + 437);
    const auto *smf_438 = buffer.data(smf + 438);
    const auto *smf_439 = buffer.data(smf + 439);
    const auto *smf_440 = buffer.data(smf + 440);
    const auto *smf_443 = buffer.data(smf + 443);
    const auto *smf_445 = buffer.data(smf + 445);
    const auto *smf_446 = buffer.data(smf + 446);
    const auto *smf_447 = buffer.data(smf + 447);
    const auto *smf_448 = buffer.data(smf + 448);
    const auto *smf_449 = buffer.data(smf + 449);
    const auto *smf_450 = buffer.data(smf + 450);
    const auto *smf_453 = buffer.data(smf + 453);
    const auto *smf_455 = buffer.data(smf + 455);
    const auto *smf_456 = buffer.data(smf + 456);
    const auto *smf_457 = buffer.data(smf + 457);
    const auto *smf_458 = buffer.data(smf + 458);
    const auto *smf_459 = buffer.data(smf + 459);
    const auto *smf_465 = buffer.data(smf + 465);
    const auto *smf_466 = buffer.data(smf + 466);
    const auto *smf_467 = buffer.data(smf + 467);
    const auto *smf_468 = buffer.data(smf + 468);
    const auto *smf_469 = buffer.data(smf + 469);
    const auto *smf_470 = buffer.data(smf + 470);
    const auto *smf_473 = buffer.data(smf + 473);

    const auto *smg1_525 = buffer.data(smg1 + 525);
    const auto *smg1_528 = buffer.data(smg1 + 528);
    const auto *smg1_530 = buffer.data(smg1 + 530);
    const auto *smg1_539 = buffer.data(smg1 + 539);
    const auto *smg1_540 = buffer.data(smg1 + 540);
    const auto *smg1_543 = buffer.data(smg1 + 543);
    const auto *smg1_675 = buffer.data(smg1 + 675);
    const auto *smg1_678 = buffer.data(smg1 + 678);
    const auto *smg1_680 = buffer.data(smg1 + 680);
    const auto *smg1_685 = buffer.data(smg1 + 685);
    const auto *smg1_687 = buffer.data(smg1 + 687);
    const auto *smg1_689 = buffer.data(smg1 + 689);
    const auto *smg1_695 = buffer.data(smg1 + 695);
    const auto *smg1_700 = buffer.data(smg1 + 700);
    const auto *smg1_702 = buffer.data(smg1 + 702);
    const auto *smg1_704 = buffer.data(smg1 + 704);
    const auto *smg1_705 = buffer.data(smg1 + 705);
    const auto *smg1_708 = buffer.data(smg1 + 708);

    const auto *snd0_237 = buffer.data(snd0 + 237);
    const auto *snd0_239 = buffer.data(snd0 + 239);
    const auto *snd0_240 = buffer.data(snd0 + 240);
    const auto *snd0_243 = buffer.data(snd0 + 243);
    const auto *snd0_245 = buffer.data(snd0 + 245);
    const auto *snd0_246 = buffer.data(snd0 + 246);
    const auto *snd0_249 = buffer.data(snd0 + 249);
    const auto *snd0_251 = buffer.data(snd0 + 251);
    const auto *snd0_252 = buffer.data(snd0 + 252);
    const auto *snd0_255 = buffer.data(snd0 + 255);
    const auto *snd0_257 = buffer.data(snd0 + 257);
    const auto *snd0_261 = buffer.data(snd0 + 261);
    const auto *snd0_263 = buffer.data(snd0 + 263);
    const auto *snd0_264 = buffer.data(snd0 + 264);
    const auto *snd0_267 = buffer.data(snd0 + 267);
    const auto *snd0_269 = buffer.data(snd0 + 269);

    const auto *snd1_237 = buffer.data(snd1 + 237);
    const auto *snd1_239 = buffer.data(snd1 + 239);
    const auto *snd1_240 = buffer.data(snd1 + 240);
    const auto *snd1_243 = buffer.data(snd1 + 243);
    const auto *snd1_245 = buffer.data(snd1 + 245);
    const auto *snd1_246 = buffer.data(snd1 + 246);
    const auto *snd1_249 = buffer.data(snd1 + 249);
    const auto *snd1_251 = buffer.data(snd1 + 251);
    const auto *snd1_252 = buffer.data(snd1 + 252);
    const auto *snd1_255 = buffer.data(snd1 + 255);
    const auto *snd1_257 = buffer.data(snd1 + 257);
    const auto *snd1_261 = buffer.data(snd1 + 261);
    const auto *snd1_263 = buffer.data(snd1 + 263);
    const auto *snd1_264 = buffer.data(snd1 + 264);
    const auto *snd1_267 = buffer.data(snd1 + 267);
    const auto *snd1_269 = buffer.data(snd1 + 269);

    const auto *snf_396 = buffer.data(snf + 396);
    const auto *snf_398 = buffer.data(snf + 398);
    const auto *snf_399 = buffer.data(snf + 399);
    const auto *snf_400 = buffer.data(snf + 400);
    const auto *snf_402 = buffer.data(snf + 402);
    const auto *snf_403 = buffer.data(snf + 403);
    const auto *snf_405 = buffer.data(snf + 405);
    const auto *snf_406 = buffer.data(snf + 406);
    const auto *snf_407 = buffer.data(snf + 407);
    const auto *snf_408 = buffer.data(snf + 408);
    const auto *snf_409 = buffer.data(snf + 409);
    const auto *snf_410 = buffer.data(snf + 410);
    const auto *snf_412 = buffer.data(snf + 412);
    const auto *snf_413 = buffer.data(snf + 413);
    const auto *snf_415 = buffer.data(snf + 415);
    const auto *snf_416 = buffer.data(snf + 416);
    const auto *snf_417 = buffer.data(snf + 417);
    const auto *snf_418 = buffer.data(snf + 418);
    const auto *snf_419 = buffer.data(snf + 419);
    const auto *snf_420 = buffer.data(snf + 420);
    const auto *snf_422 = buffer.data(snf + 422);
    const auto *snf_423 = buffer.data(snf + 423);
    const auto *snf_425 = buffer.data(snf + 425);
    const auto *snf_426 = buffer.data(snf + 426);
    const auto *snf_427 = buffer.data(snf + 427);
    const auto *snf_428 = buffer.data(snf + 428);
    const auto *snf_429 = buffer.data(snf + 429);
    const auto *snf_430 = buffer.data(snf + 430);
    const auto *snf_432 = buffer.data(snf + 432);
    const auto *snf_436 = buffer.data(snf + 436);
    const auto *snf_437 = buffer.data(snf + 437);
    const auto *snf_438 = buffer.data(snf + 438);
    const auto *snf_439 = buffer.data(snf + 439);
    const auto *snf_440 = buffer.data(snf + 440);
    const auto *snf_442 = buffer.data(snf + 442);
    const auto *snf_443 = buffer.data(snf + 443);
    const auto *snf_445 = buffer.data(snf + 445);
    const auto *snf_446 = buffer.data(snf + 446);
    const auto *snf_447 = buffer.data(snf + 447);
    const auto *snf_448 = buffer.data(snf + 448);
    const auto *snf_449 = buffer.data(snf + 449);
    const auto *snf_450 = buffer.data(snf + 450);
    const auto *snf_452 = buffer.data(snf + 452);
    const auto *snf_456 = buffer.data(snf + 456);
    const auto *snf_457 = buffer.data(snf + 457);
    const auto *snf_458 = buffer.data(snf + 458);
    const auto *snf_459 = buffer.data(snf + 459);
    const auto *snf_460 = buffer.data(snf + 460);
    const auto *snf_462 = buffer.data(snf + 462);
    const auto *snf_466 = buffer.data(snf + 466);
    const auto *snf_467 = buffer.data(snf + 467);
    const auto *snf_468 = buffer.data(snf + 468);
    const auto *snf_469 = buffer.data(snf + 469);
    const auto *snf_470 = buffer.data(snf + 470);
    const auto *snf_472 = buffer.data(snf + 472);

#pragma omp simd aligned(t_594, t_595, t_596, pc_x, pc_y, pc_z, smf_306, smf_316, smf_399, \
                         snd0_237, snd1_237, snf_396, snf_399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_594[k] = f_8 * smf_399[k]
                   + f_3 * pc_x[k] * snf_399[k];

        t_595[k] = f_15 * smf_316[k]
                   + f_1 * snd0_237[k]
                   - f_2 * snd1_237[k]
                   + f_3 * pc_y[k] * snf_396[k];

        t_596[k] = f_12 * smf_306[k]
                   + f_3 * pc_z[k] * snf_396[k];
    }

#pragma omp simd aligned(t_597, t_598, t_599, pc_y, pc_z, smf_309, smf_318, smf_319, snd0_239, \
                         snd1_239, snf_398, snf_399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_597[k] = f_15 * smf_318[k]
                   + f_4 * snd0_239[k]
                   - f_5 * snd1_239[k]
                   + f_3 * pc_y[k] * snf_398[k];

        t_598[k] = f_15 * smf_319[k]
                   + f_3 * pc_y[k] * snf_399[k];

        t_599[k] = f_12 * smf_309[k]
                   + f_1 * snd0_239[k]
                   - f_2 * snd1_239[k]
                   + f_3 * pc_z[k] * snf_399[k];
    }

#pragma omp simd aligned(t_600, t_601, t_602, pc_x, pc_y, pc_z, smf_310, smf_320, smf_400, \
                         snd0_240, snd1_240, snf_400 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = f_8 * smf_400[k]
                   + f_1 * snd0_240[k]
                   - f_2 * snd1_240[k]
                   + f_3 * pc_x[k] * snf_400[k];

        t_601[k] = f_14 * smf_320[k]
                   + f_3 * pc_y[k] * snf_400[k];

        t_602[k] = f_14 * smf_310[k]
                   + f_3 * pc_z[k] * snf_400[k];
    }

#pragma omp simd aligned(t_603, t_604, t_605, pc_x, pc_y, smf_322, smf_403, smf_405, snd0_243, \
                         snd0_245, snd1_243, snd1_245, snf_402, snf_403, \
                         snf_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_603[k] = f_8 * smf_403[k]
                   + f_4 * snd0_243[k]
                   - f_5 * snd1_243[k]
                   + f_3 * pc_x[k] * snf_403[k];

        t_604[k] = f_14 * smf_322[k]
                   + f_3 * pc_y[k] * snf_402[k];

        t_605[k] = f_8 * smf_405[k]
                   + f_4 * snd0_245[k]
                   - f_5 * snd1_245[k]
                   + f_3 * pc_x[k] * snf_405[k];
    }

#pragma omp simd aligned(t_606, t_607, t_608, t_609, pc_x, smf_406, smf_407, smf_408, smf_409, \
                         snf_406, snf_407, snf_408, snf_409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_606[k] = f_8 * smf_406[k]
                   + f_3 * pc_x[k] * snf_406[k];

        t_607[k] = f_8 * smf_407[k]
                   + f_3 * pc_x[k] * snf_407[k];

        t_608[k] = f_8 * smf_408[k]
                   + f_3 * pc_x[k] * snf_408[k];

        t_609[k] = f_8 * smf_409[k]
                   + f_3 * pc_x[k] * snf_409[k];
    }

#pragma omp simd aligned(t_610, t_611, t_612, pc_y, pc_z, smf_316, smf_326, smf_328, snd0_243, \
                         snd0_245, snd1_243, snd1_245, snf_406, \
                         snf_408 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_610[k] = f_14 * smf_326[k]
                   + f_1 * snd0_243[k]
                   - f_2 * snd1_243[k]
                   + f_3 * pc_y[k] * snf_406[k];

        t_611[k] = f_14 * smf_316[k]
                   + f_3 * pc_z[k] * snf_406[k];

        t_612[k] = f_14 * smf_328[k]
                   + f_4 * snd0_245[k]
                   - f_5 * snd1_245[k]
                   + f_3 * pc_y[k] * snf_408[k];
    }

#pragma omp simd aligned(t_613, t_614, t_615, pc_x, pc_y, pc_z, smf_319, smf_329, smf_410, \
                         snd0_245, snd0_246, snd1_245, snd1_246, snf_409, \
                         snf_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_613[k] = f_14 * smf_329[k]
                   + f_3 * pc_y[k] * snf_409[k];

        t_614[k] = f_14 * smf_319[k]
                   + f_1 * snd0_245[k]
                   - f_2 * snd1_245[k]
                   + f_3 * pc_z[k] * snf_409[k];

        t_615[k] = f_8 * smf_410[k]
                   + f_1 * snd0_246[k]
                   - f_2 * snd1_246[k]
                   + f_3 * pc_x[k] * snf_410[k];
    }

#pragma omp simd aligned(t_616, t_617, t_618, t_619, pc_x, pc_y, pc_z, smf_320, smf_330, \
                         smf_332, smf_413, snd0_249, snd1_249, snf_410, snf_412, \
                         snf_413 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_616[k] = f_12 * smf_330[k]
                   + f_3 * pc_y[k] * snf_410[k];

        t_617[k] = f_15 * smf_320[k]
                   + f_3 * pc_z[k] * snf_410[k];

        t_618[k] = f_8 * smf_413[k]
                   + f_4 * snd0_249[k]
                   - f_5 * snd1_249[k]
                   + f_3 * pc_x[k] * snf_413[k];

        t_619[k] = f_12 * smf_332[k]
                   + f_3 * pc_y[k] * snf_412[k];
    }

#pragma omp simd aligned(t_620, t_621, t_622, t_623, pc_x, smf_415, smf_416, smf_417, smf_418, \
                         snd0_251, snd1_251, snf_415, snf_416, snf_417, \
                         snf_418 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_620[k] = f_8 * smf_415[k]
                   + f_4 * snd0_251[k]
                   - f_5 * snd1_251[k]
                   + f_3 * pc_x[k] * snf_415[k];

        t_621[k] = f_8 * smf_416[k]
                   + f_3 * pc_x[k] * snf_416[k];

        t_622[k] = f_8 * smf_417[k]
                   + f_3 * pc_x[k] * snf_417[k];

        t_623[k] = f_8 * smf_418[k]
                   + f_3 * pc_x[k] * snf_418[k];
    }

#pragma omp simd aligned(t_624, t_625, t_626, pc_x, pc_y, pc_z, smf_326, smf_336, smf_419, \
                         snd0_249, snd1_249, snf_416, snf_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_624[k] = f_8 * smf_419[k]
                   + f_3 * pc_x[k] * snf_419[k];

        t_625[k] = f_12 * smf_336[k]
                   + f_1 * snd0_249[k]
                   - f_2 * snd1_249[k]
                   + f_3 * pc_y[k] * snf_416[k];

        t_626[k] = f_15 * smf_326[k]
                   + f_3 * pc_z[k] * snf_416[k];
    }

#pragma omp simd aligned(t_627, t_628, t_629, pc_y, pc_z, smf_329, smf_338, smf_339, snd0_251, \
                         snd1_251, snf_418, snf_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_627[k] = f_12 * smf_338[k]
                   + f_4 * snd0_251[k]
                   - f_5 * snd1_251[k]
                   + f_3 * pc_y[k] * snf_418[k];

        t_628[k] = f_12 * smf_339[k]
                   + f_3 * pc_y[k] * snf_419[k];

        t_629[k] = f_15 * smf_329[k]
                   + f_1 * snd0_251[k]
                   - f_2 * snd1_251[k]
                   + f_3 * pc_z[k] * snf_419[k];
    }

#pragma omp simd aligned(t_630, t_631, t_632, pc_x, pc_y, pc_z, smf_330, smf_340, smf_420, \
                         snd0_252, snd1_252, snf_420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_630[k] = f_8 * smf_420[k]
                   + f_1 * snd0_252[k]
                   - f_2 * snd1_252[k]
                   + f_3 * pc_x[k] * snf_420[k];

        t_631[k] = f_8 * smf_340[k]
                   + f_3 * pc_y[k] * snf_420[k];

        t_632[k] = f_13 * smf_330[k]
                   + f_3 * pc_z[k] * snf_420[k];
    }

#pragma omp simd aligned(t_633, t_634, t_635, pc_x, pc_y, smf_342, smf_423, smf_425, snd0_255, \
                         snd0_257, snd1_255, snd1_257, snf_422, snf_423, \
                         snf_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_633[k] = f_8 * smf_423[k]
                   + f_4 * snd0_255[k]
                   - f_5 * snd1_255[k]
                   + f_3 * pc_x[k] * snf_423[k];

        t_634[k] = f_8 * smf_342[k]
                   + f_3 * pc_y[k] * snf_422[k];

        t_635[k] = f_8 * smf_425[k]
                   + f_4 * snd0_257[k]
                   - f_5 * snd1_257[k]
                   + f_3 * pc_x[k] * snf_425[k];
    }

#pragma omp simd aligned(t_636, t_637, t_638, t_639, pc_x, smf_426, smf_427, smf_428, smf_429, \
                         snf_426, snf_427, snf_428, snf_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_636[k] = f_8 * smf_426[k]
                   + f_3 * pc_x[k] * snf_426[k];

        t_637[k] = f_8 * smf_427[k]
                   + f_3 * pc_x[k] * snf_427[k];

        t_638[k] = f_8 * smf_428[k]
                   + f_3 * pc_x[k] * snf_428[k];

        t_639[k] = f_8 * smf_429[k]
                   + f_3 * pc_x[k] * snf_429[k];
    }

#pragma omp simd aligned(t_640, t_641, t_642, pc_y, pc_z, smf_336, smf_346, smf_348, snd0_255, \
                         snd0_257, snd1_255, snd1_257, snf_426, \
                         snf_428 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_640[k] = f_8 * smf_346[k]
                   + f_1 * snd0_255[k]
                   - f_2 * snd1_255[k]
                   + f_3 * pc_y[k] * snf_426[k];

        t_641[k] = f_13 * smf_336[k]
                   + f_3 * pc_z[k] * snf_426[k];

        t_642[k] = f_8 * smf_348[k]
                   + f_4 * snd0_257[k]
                   - f_5 * snd1_257[k]
                   + f_3 * pc_y[k] * snf_428[k];
    }

#pragma omp simd aligned(t_643, t_644, t_645, t_646, pb_y, pc_y, pc_z, smg0_525, smf_339, \
                         smf_349, smf_350, smg1_525, snd0_257, snd1_257, snf_429, \
                         snf_430 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_643[k] = f_8 * smf_349[k]
                   + f_3 * pc_y[k] * snf_429[k];

        t_644[k] = f_13 * smf_339[k]
                   + f_1 * snd0_257[k]
                   - f_2 * snd1_257[k]
                   + f_3 * pc_z[k] * snf_429[k];

        t_645[k] = pb_y[k] * smg0_525[k]
                   - f_6 * pc_y[k] * smg1_525[k];

        t_646[k] = f_7 * smf_350[k]
                   + f_3 * pc_y[k] * snf_430[k];
    }

#pragma omp simd aligned(t_647, t_648, t_649, t_650, pb_y, pc_y, pc_z, smg0_528, smg0_530, \
                         smf_340, smf_351, smf_352, smg1_528, smg1_530, snf_430, \
                         snf_432 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_647[k] = f_11 * smf_340[k]
                   + f_3 * pc_z[k] * snf_430[k];

        t_648[k] = pb_y[k] * smg0_528[k]
                   + f_8 * smf_351[k]
                   - f_6 * pc_y[k] * smg1_528[k];

        t_649[k] = f_7 * smf_352[k]
                   + f_3 * pc_y[k] * snf_432[k];

        t_650[k] = pb_y[k] * smg0_530[k]
                   - f_6 * pc_y[k] * smg1_530[k];
    }

#pragma omp simd aligned(t_651, t_652, t_653, t_654, pc_x, smf_436, smf_437, smf_438, smf_439, \
                         snf_436, snf_437, snf_438, snf_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_651[k] = f_8 * smf_436[k]
                   + f_3 * pc_x[k] * snf_436[k];

        t_652[k] = f_8 * smf_437[k]
                   + f_3 * pc_x[k] * snf_437[k];

        t_653[k] = f_8 * smf_438[k]
                   + f_3 * pc_x[k] * snf_438[k];

        t_654[k] = f_8 * smf_439[k]
                   + f_3 * pc_x[k] * snf_439[k];
    }

#pragma omp simd aligned(t_655, t_656, t_657, pc_y, pc_z, smf_346, smf_356, smf_358, snd0_261, \
                         snd0_263, snd1_261, snd1_263, snf_436, \
                         snf_438 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_655[k] = f_7 * smf_356[k]
                   + f_1 * snd0_261[k]
                   - f_2 * snd1_261[k]
                   + f_3 * pc_y[k] * snf_436[k];

        t_656[k] = f_11 * smf_346[k]
                   + f_3 * pc_z[k] * snf_436[k];

        t_657[k] = f_7 * smf_358[k]
                   + f_4 * snd0_263[k]
                   - f_5 * snd1_263[k]
                   + f_3 * pc_y[k] * snf_438[k];
    }

#pragma omp simd aligned(t_658, t_659, t_660, t_661, pb_y, pc_x, pc_y, smg0_539, smf_359, \
                         smf_440, smg1_539, snd0_264, snd1_264, snf_439, \
                         snf_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_658[k] = f_7 * smf_359[k]
                   + f_3 * pc_y[k] * snf_439[k];

        t_659[k] = pb_y[k] * smg0_539[k]
                   - f_6 * pc_y[k] * smg1_539[k];

        t_660[k] = f_8 * smf_440[k]
                   + f_1 * snd0_264[k]
                   - f_2 * snd1_264[k]
                   + f_3 * pc_x[k] * snf_440[k];

        t_661[k] = f_3 * pc_y[k] * snf_440[k];
    }

#pragma omp simd aligned(t_662, t_663, t_664, pc_x, pc_y, pc_z, smf_350, smf_443, snd0_267, \
                         snd1_267, snf_440, snf_442, snf_443 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_662[k] = f_10 * smf_350[k]
                   + f_3 * pc_z[k] * snf_440[k];

        t_663[k] = f_8 * smf_443[k]
                   + f_4 * snd0_267[k]
                   - f_5 * snd1_267[k]
                   + f_3 * pc_x[k] * snf_443[k];

        t_664[k] = f_3 * pc_y[k] * snf_442[k];
    }

#pragma omp simd aligned(t_665, t_666, t_667, t_668, pc_x, smf_445, smf_446, smf_447, smf_448, \
                         snd0_269, snd1_269, snf_445, snf_446, snf_447, \
                         snf_448 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_665[k] = f_8 * smf_445[k]
                   + f_4 * snd0_269[k]
                   - f_5 * snd1_269[k]
                   + f_3 * pc_x[k] * snf_445[k];

        t_666[k] = f_8 * smf_446[k]
                   + f_3 * pc_x[k] * snf_446[k];

        t_667[k] = f_8 * smf_447[k]
                   + f_3 * pc_x[k] * snf_447[k];

        t_668[k] = f_8 * smf_448[k]
                   + f_3 * pc_x[k] * snf_448[k];
    }

#pragma omp simd aligned(t_669, t_670, t_671, t_672, pc_x, pc_y, pc_z, smf_356, smf_449, \
                         snd0_267, snd0_269, snd1_267, snd1_269, snf_446, snf_448, \
                         snf_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_669[k] = f_8 * smf_449[k]
                   + f_3 * pc_x[k] * snf_449[k];

        t_670[k] = f_1 * snd0_267[k]
                   - f_2 * snd1_267[k]
                   + f_3 * pc_y[k] * snf_446[k];

        t_671[k] = f_10 * smf_356[k]
                   + f_3 * pc_z[k] * snf_446[k];

        t_672[k] = f_4 * snd0_269[k]
                   - f_5 * snd1_269[k]
                   + f_3 * pc_y[k] * snf_448[k];
    }

#pragma omp simd aligned(t_673, t_674, t_675, pb_x, pc_x, pc_y, pc_z, smg0_675, smf_359, \
                         smf_450, smg1_675, snd0_269, snd1_269, \
                         snf_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_673[k] = f_3 * pc_y[k] * snf_449[k];

        t_674[k] = f_10 * smf_359[k]
                   + f_1 * snd0_269[k]
                   - f_2 * snd1_269[k]
                   + f_3 * pc_z[k] * snf_449[k];

        t_675[k] = pb_x[k] * smg0_675[k]
                   + f_14 * smf_450[k]
                   - f_6 * pc_x[k] * smg1_675[k];
    }

#pragma omp simd aligned(t_676, t_677, t_678, t_679, pb_x, pc_x, pc_y, pc_z, smg0_678, \
                         smf_360, smf_362, smf_453, smg1_678, snf_450, \
                         snf_452 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_676[k] = f_9 * smf_360[k]
                   + f_3 * pc_y[k] * snf_450[k];

        t_677[k] = f_3 * pc_z[k] * snf_450[k];

        t_678[k] = pb_x[k] * smg0_678[k]
                   + f_8 * smf_453[k]
                   - f_6 * pc_x[k] * smg1_678[k];

        t_679[k] = f_9 * smf_362[k]
                   + f_3 * pc_y[k] * snf_452[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, t_683, pb_x, pc_x, smg0_680, smf_455, smf_456, \
                         smf_457, smf_458, smg1_680, snf_456, snf_457, \
                         snf_458 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = pb_x[k] * smg0_680[k]
                   + f_8 * smf_455[k]
                   - f_6 * pc_x[k] * smg1_680[k];

        t_681[k] = f_7 * smf_456[k]
                   + f_3 * pc_x[k] * snf_456[k];

        t_682[k] = f_7 * smf_457[k]
                   + f_3 * pc_x[k] * snf_457[k];

        t_683[k] = f_7 * smf_458[k]
                   + f_3 * pc_x[k] * snf_458[k];
    }

#pragma omp simd aligned(t_684, t_685, t_686, t_687, pb_x, pc_x, pc_z, smg0_685, smg0_687, \
                         smf_459, smg1_685, smg1_687, snf_456, \
                         snf_459 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_684[k] = f_7 * smf_459[k]
                   + f_3 * pc_x[k] * snf_459[k];

        t_685[k] = pb_x[k] * smg0_685[k]
                   - f_6 * pc_x[k] * smg1_685[k];

        t_686[k] = f_3 * pc_z[k] * snf_456[k];

        t_687[k] = pb_x[k] * smg0_687[k]
                   - f_6 * pc_x[k] * smg1_687[k];
    }

#pragma omp simd aligned(t_688, t_689, t_690, pb_x, pb_z, pc_x, pc_y, pc_z, smg0_540, \
                         smg0_689, smf_369, smg1_540, smg1_689, \
                         snf_459 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_688[k] = f_9 * smf_369[k]
                   + f_3 * pc_y[k] * snf_459[k];

        t_689[k] = pb_x[k] * smg0_689[k]
                   - f_6 * pc_x[k] * smg1_689[k];

        t_690[k] = pb_z[k] * smg0_540[k]
                   - f_6 * pc_z[k] * smg1_540[k];
    }

#pragma omp simd aligned(t_691, t_692, t_693, t_694, pb_z, pc_y, pc_z, smg0_543, smf_360, \
                         smf_370, smf_372, smg1_543, snf_460, snf_462 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_691[k] = f_10 * smf_370[k]
                   + f_3 * pc_y[k] * snf_460[k];

        t_692[k] = f_7 * smf_360[k]
                   + f_3 * pc_z[k] * snf_460[k];

        t_693[k] = pb_z[k] * smg0_543[k]
                   - f_6 * pc_z[k] * smg1_543[k];

        t_694[k] = f_10 * smf_372[k]
                   + f_3 * pc_y[k] * snf_462[k];
    }

#pragma omp simd aligned(t_695, t_696, t_697, t_698, pb_x, pc_x, smg0_695, smf_465, smf_466, \
                         smf_467, smf_468, smg1_695, snf_466, snf_467, \
                         snf_468 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_695[k] = pb_x[k] * smg0_695[k]
                   + f_8 * smf_465[k]
                   - f_6 * pc_x[k] * smg1_695[k];

        t_696[k] = f_7 * smf_466[k]
                   + f_3 * pc_x[k] * snf_466[k];

        t_697[k] = f_7 * smf_467[k]
                   + f_3 * pc_x[k] * snf_467[k];

        t_698[k] = f_7 * smf_468[k]
                   + f_3 * pc_x[k] * snf_468[k];
    }

#pragma omp simd aligned(t_699, t_700, t_701, t_702, pb_x, pc_x, pc_z, smg0_700, smg0_702, \
                         smf_366, smf_469, smg1_700, smg1_702, snf_466, \
                         snf_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_699[k] = f_7 * smf_469[k]
                   + f_3 * pc_x[k] * snf_469[k];

        t_700[k] = pb_x[k] * smg0_700[k]
                   - f_6 * pc_x[k] * smg1_700[k];

        t_701[k] = f_7 * smf_366[k]
                   + f_3 * pc_z[k] * snf_466[k];

        t_702[k] = pb_x[k] * smg0_702[k]
                   - f_6 * pc_x[k] * smg1_702[k];
    }

#pragma omp simd aligned(t_703, t_704, t_705, t_706, pb_x, pc_x, pc_y, smg0_704, smg0_705, \
                         smf_379, smf_380, smf_470, smg1_704, smg1_705, snf_469, \
                         snf_470 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_703[k] = f_10 * smf_379[k]
                   + f_3 * pc_y[k] * snf_469[k];

        t_704[k] = pb_x[k] * smg0_704[k]
                   - f_6 * pc_x[k] * smg1_704[k];

        t_705[k] = pb_x[k] * smg0_705[k]
                   + f_14 * smf_470[k]
                   - f_6 * pc_x[k] * smg1_705[k];

        t_706[k] = f_11 * smf_380[k]
                   + f_3 * pc_y[k] * snf_470[k];
    }

#pragma omp simd aligned(t_707, t_708, t_709, pb_x, pc_x, pc_y, pc_z, smg0_708, smf_370, \
                         smf_382, smf_473, smg1_708, snf_470, snf_472 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_707[k] = f_8 * smf_370[k]
                   + f_3 * pc_z[k] * snf_470[k];

        t_708[k] = pb_x[k] * smg0_708[k]
                   + f_8 * smf_473[k]
                   - f_6 * pc_x[k] * smg1_708[k];

        t_709[k] = f_11 * smf_382[k]
                   + f_3 * pc_y[k] * snf_472[k];
    }
}

static auto
compute_prim_sng_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t smg0,
                                                          const size_t smf, const size_t smg1,
                                                          const size_t snd0, const size_t snd1,
                                                          const size_t snf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.0 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 4.5 / q;
    const auto f_10 = 4.0 / q;
    const auto f_11 = 3.5 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 3.0 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 2.5 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *smg0_660 = buffer.data(smg0 + 660);
    const auto *smg0_665 = buffer.data(smg0 + 665);
    const auto *smg0_710 = buffer.data(smg0 + 710);
    const auto *smg0_715 = buffer.data(smg0 + 715);
    const auto *smg0_717 = buffer.data(smg0 + 717);
    const auto *smg0_719 = buffer.data(smg0 + 719);
    const auto *smg0_720 = buffer.data(smg0 + 720);
    const auto *smg0_723 = buffer.data(smg0 + 723);
    const auto *smg0_725 = buffer.data(smg0 + 725);
    const auto *smg0_730 = buffer.data(smg0 + 730);
    const auto *smg0_732 = buffer.data(smg0 + 732);
    const auto *smg0_734 = buffer.data(smg0 + 734);
    const auto *smg0_735 = buffer.data(smg0 + 735);
    const auto *smg0_738 = buffer.data(smg0 + 738);
    const auto *smg0_740 = buffer.data(smg0 + 740);
    const auto *smg0_745 = buffer.data(smg0 + 745);
    const auto *smg0_747 = buffer.data(smg0 + 747);
    const auto *smg0_749 = buffer.data(smg0 + 749);
    const auto *smg0_750 = buffer.data(smg0 + 750);
    const auto *smg0_753 = buffer.data(smg0 + 753);
    const auto *smg0_755 = buffer.data(smg0 + 755);
    const auto *smg0_760 = buffer.data(smg0 + 760);
    const auto *smg0_762 = buffer.data(smg0 + 762);
    const auto *smg0_764 = buffer.data(smg0 + 764);
    const auto *smg0_765 = buffer.data(smg0 + 765);
    const auto *smg0_768 = buffer.data(smg0 + 768);
    const auto *smg0_770 = buffer.data(smg0 + 770);
    const auto *smg0_775 = buffer.data(smg0 + 775);
    const auto *smg0_777 = buffer.data(smg0 + 777);
    const auto *smg0_779 = buffer.data(smg0 + 779);
    const auto *smg0_780 = buffer.data(smg0 + 780);
    const auto *smg0_783 = buffer.data(smg0 + 783);
    const auto *smg0_785 = buffer.data(smg0 + 785);
    const auto *smg0_790 = buffer.data(smg0 + 790);
    const auto *smg0_792 = buffer.data(smg0 + 792);
    const auto *smg0_794 = buffer.data(smg0 + 794);
    const auto *smg0_798 = buffer.data(smg0 + 798);
    const auto *smg0_805 = buffer.data(smg0 + 805);
    const auto *smg0_807 = buffer.data(smg0 + 807);
    const auto *smg0_809 = buffer.data(smg0 + 809);
    const auto *smg0_810 = buffer.data(smg0 + 810);
    const auto *smg0_813 = buffer.data(smg0 + 813);
    const auto *smg0_815 = buffer.data(smg0 + 815);
    const auto *smg0_820 = buffer.data(smg0 + 820);
    const auto *smg0_822 = buffer.data(smg0 + 822);
    const auto *smg0_824 = buffer.data(smg0 + 824);

    const auto *smf_376 = buffer.data(smf + 376);
    const auto *smf_380 = buffer.data(smf + 380);
    const auto *smf_386 = buffer.data(smf + 386);
    const auto *smf_389 = buffer.data(smf + 389);
    const auto *smf_390 = buffer.data(smf + 390);
    const auto *smf_392 = buffer.data(smf + 392);
    const auto *smf_396 = buffer.data(smf + 396);
    const auto *smf_399 = buffer.data(smf + 399);
    const auto *smf_400 = buffer.data(smf + 400);
    const auto *smf_402 = buffer.data(smf + 402);
    const auto *smf_406 = buffer.data(smf + 406);
    const auto *smf_409 = buffer.data(smf + 409);
    const auto *smf_410 = buffer.data(smf + 410);
    const auto *smf_412 = buffer.data(smf + 412);
    const auto *smf_416 = buffer.data(smf + 416);
    const auto *smf_419 = buffer.data(smf + 419);
    const auto *smf_420 = buffer.data(smf + 420);
    const auto *smf_422 = buffer.data(smf + 422);
    const auto *smf_426 = buffer.data(smf + 426);
    const auto *smf_429 = buffer.data(smf + 429);
    const auto *smf_430 = buffer.data(smf + 430);
    const auto *smf_432 = buffer.data(smf + 432);
    const auto *smf_436 = buffer.data(smf + 436);
    const auto *smf_439 = buffer.data(smf + 439);
    const auto *smf_440 = buffer.data(smf + 440);
    const auto *smf_442 = buffer.data(smf + 442);
    const auto *smf_446 = buffer.data(smf + 446);
    const auto *smf_449 = buffer.data(smf + 449);
    const auto *smf_450 = buffer.data(smf + 450);
    const auto *smf_452 = buffer.data(smf + 452);
    const auto *smf_456 = buffer.data(smf + 456);
    const auto *smf_475 = buffer.data(smf + 475);
    const auto *smf_476 = buffer.data(smf + 476);
    const auto *smf_477 = buffer.data(smf + 477);
    const auto *smf_478 = buffer.data(smf + 478);
    const auto *smf_479 = buffer.data(smf + 479);
    const auto *smf_480 = buffer.data(smf + 480);
    const auto *smf_483 = buffer.data(smf + 483);
    const auto *smf_485 = buffer.data(smf + 485);
    const auto *smf_486 = buffer.data(smf + 486);
    const auto *smf_487 = buffer.data(smf + 487);
    const auto *smf_488 = buffer.data(smf + 488);
    const auto *smf_489 = buffer.data(smf + 489);
    const auto *smf_490 = buffer.data(smf + 490);
    const auto *smf_493 = buffer.data(smf + 493);
    const auto *smf_495 = buffer.data(smf + 495);
    const auto *smf_496 = buffer.data(smf + 496);
    const auto *smf_497 = buffer.data(smf + 497);
    const auto *smf_498 = buffer.data(smf + 498);
    const auto *smf_499 = buffer.data(smf + 499);
    const auto *smf_500 = buffer.data(smf + 500);
    const auto *smf_503 = buffer.data(smf + 503);
    const auto *smf_505 = buffer.data(smf + 505);
    const auto *smf_506 = buffer.data(smf + 506);
    const auto *smf_507 = buffer.data(smf + 507);
    const auto *smf_508 = buffer.data(smf + 508);
    const auto *smf_509 = buffer.data(smf + 509);
    const auto *smf_510 = buffer.data(smf + 510);
    const auto *smf_513 = buffer.data(smf + 513);
    const auto *smf_515 = buffer.data(smf + 515);
    const auto *smf_516 = buffer.data(smf + 516);
    const auto *smf_517 = buffer.data(smf + 517);
    const auto *smf_518 = buffer.data(smf + 518);
    const auto *smf_519 = buffer.data(smf + 519);
    const auto *smf_520 = buffer.data(smf + 520);
    const auto *smf_523 = buffer.data(smf + 523);
    const auto *smf_525 = buffer.data(smf + 525);
    const auto *smf_526 = buffer.data(smf + 526);
    const auto *smf_527 = buffer.data(smf + 527);
    const auto *smf_528 = buffer.data(smf + 528);
    const auto *smf_529 = buffer.data(smf + 529);
    const auto *smf_533 = buffer.data(smf + 533);
    const auto *smf_536 = buffer.data(smf + 536);
    const auto *smf_537 = buffer.data(smf + 537);
    const auto *smf_538 = buffer.data(smf + 538);
    const auto *smf_539 = buffer.data(smf + 539);
    const auto *smf_540 = buffer.data(smf + 540);
    const auto *smf_543 = buffer.data(smf + 543);
    const auto *smf_545 = buffer.data(smf + 545);
    const auto *smf_546 = buffer.data(smf + 546);
    const auto *smf_547 = buffer.data(smf + 547);
    const auto *smf_548 = buffer.data(smf + 548);
    const auto *smf_549 = buffer.data(smf + 549);

    const auto *smg1_660 = buffer.data(smg1 + 660);
    const auto *smg1_665 = buffer.data(smg1 + 665);
    const auto *smg1_710 = buffer.data(smg1 + 710);
    const auto *smg1_715 = buffer.data(smg1 + 715);
    const auto *smg1_717 = buffer.data(smg1 + 717);
    const auto *smg1_719 = buffer.data(smg1 + 719);
    const auto *smg1_720 = buffer.data(smg1 + 720);
    const auto *smg1_723 = buffer.data(smg1 + 723);
    const auto *smg1_725 = buffer.data(smg1 + 725);
    const auto *smg1_730 = buffer.data(smg1 + 730);
    const auto *smg1_732 = buffer.data(smg1 + 732);
    const auto *smg1_734 = buffer.data(smg1 + 734);
    const auto *smg1_735 = buffer.data(smg1 + 735);
    const auto *smg1_738 = buffer.data(smg1 + 738);
    const auto *smg1_740 = buffer.data(smg1 + 740);
    const auto *smg1_745 = buffer.data(smg1 + 745);
    const auto *smg1_747 = buffer.data(smg1 + 747);
    const auto *smg1_749 = buffer.data(smg1 + 749);
    const auto *smg1_750 = buffer.data(smg1 + 750);
    const auto *smg1_753 = buffer.data(smg1 + 753);
    const auto *smg1_755 = buffer.data(smg1 + 755);
    const auto *smg1_760 = buffer.data(smg1 + 760);
    const auto *smg1_762 = buffer.data(smg1 + 762);
    const auto *smg1_764 = buffer.data(smg1 + 764);
    const auto *smg1_765 = buffer.data(smg1 + 765);
    const auto *smg1_768 = buffer.data(smg1 + 768);
    const auto *smg1_770 = buffer.data(smg1 + 770);
    const auto *smg1_775 = buffer.data(smg1 + 775);
    const auto *smg1_777 = buffer.data(smg1 + 777);
    const auto *smg1_779 = buffer.data(smg1 + 779);
    const auto *smg1_780 = buffer.data(smg1 + 780);
    const auto *smg1_783 = buffer.data(smg1 + 783);
    const auto *smg1_785 = buffer.data(smg1 + 785);
    const auto *smg1_790 = buffer.data(smg1 + 790);
    const auto *smg1_792 = buffer.data(smg1 + 792);
    const auto *smg1_794 = buffer.data(smg1 + 794);
    const auto *smg1_798 = buffer.data(smg1 + 798);
    const auto *smg1_805 = buffer.data(smg1 + 805);
    const auto *smg1_807 = buffer.data(smg1 + 807);
    const auto *smg1_809 = buffer.data(smg1 + 809);
    const auto *smg1_810 = buffer.data(smg1 + 810);
    const auto *smg1_813 = buffer.data(smg1 + 813);
    const auto *smg1_815 = buffer.data(smg1 + 815);
    const auto *smg1_820 = buffer.data(smg1 + 820);
    const auto *smg1_822 = buffer.data(smg1 + 822);
    const auto *smg1_824 = buffer.data(smg1 + 824);

    const auto *snd0_330 = buffer.data(snd0 + 330);
    const auto *snd0_333 = buffer.data(snd0 + 333);
    const auto *snd0_335 = buffer.data(snd0 + 335);

    const auto *snd1_330 = buffer.data(snd1 + 330);
    const auto *snd1_333 = buffer.data(snd1 + 333);
    const auto *snd1_335 = buffer.data(snd1 + 335);

    const auto *snf_476 = buffer.data(snf + 476);
    const auto *snf_477 = buffer.data(snf + 477);
    const auto *snf_478 = buffer.data(snf + 478);
    const auto *snf_479 = buffer.data(snf + 479);
    const auto *snf_480 = buffer.data(snf + 480);
    const auto *snf_482 = buffer.data(snf + 482);
    const auto *snf_486 = buffer.data(snf + 486);
    const auto *snf_487 = buffer.data(snf + 487);
    const auto *snf_488 = buffer.data(snf + 488);
    const auto *snf_489 = buffer.data(snf + 489);
    const auto *snf_490 = buffer.data(snf + 490);
    const auto *snf_492 = buffer.data(snf + 492);
    const auto *snf_496 = buffer.data(snf + 496);
    const auto *snf_497 = buffer.data(snf + 497);
    const auto *snf_498 = buffer.data(snf + 498);
    const auto *snf_499 = buffer.data(snf + 499);
    const auto *snf_500 = buffer.data(snf + 500);
    const auto *snf_502 = buffer.data(snf + 502);
    const auto *snf_506 = buffer.data(snf + 506);
    const auto *snf_507 = buffer.data(snf + 507);
    const auto *snf_508 = buffer.data(snf + 508);
    const auto *snf_509 = buffer.data(snf + 509);
    const auto *snf_510 = buffer.data(snf + 510);
    const auto *snf_512 = buffer.data(snf + 512);
    const auto *snf_516 = buffer.data(snf + 516);
    const auto *snf_517 = buffer.data(snf + 517);
    const auto *snf_518 = buffer.data(snf + 518);
    const auto *snf_519 = buffer.data(snf + 519);
    const auto *snf_520 = buffer.data(snf + 520);
    const auto *snf_522 = buffer.data(snf + 522);
    const auto *snf_526 = buffer.data(snf + 526);
    const auto *snf_527 = buffer.data(snf + 527);
    const auto *snf_528 = buffer.data(snf + 528);
    const auto *snf_529 = buffer.data(snf + 529);
    const auto *snf_530 = buffer.data(snf + 530);
    const auto *snf_532 = buffer.data(snf + 532);
    const auto *snf_536 = buffer.data(snf + 536);
    const auto *snf_537 = buffer.data(snf + 537);
    const auto *snf_538 = buffer.data(snf + 538);
    const auto *snf_539 = buffer.data(snf + 539);
    const auto *snf_540 = buffer.data(snf + 540);
    const auto *snf_542 = buffer.data(snf + 542);
    const auto *snf_546 = buffer.data(snf + 546);
    const auto *snf_547 = buffer.data(snf + 547);
    const auto *snf_548 = buffer.data(snf + 548);
    const auto *snf_549 = buffer.data(snf + 549);
    const auto *snf_550 = buffer.data(snf + 550);
    const auto *snf_552 = buffer.data(snf + 552);
    const auto *snf_553 = buffer.data(snf + 553);
    const auto *snf_555 = buffer.data(snf + 555);
    const auto *snf_556 = buffer.data(snf + 556);
    const auto *snf_557 = buffer.data(snf + 557);
    const auto *snf_558 = buffer.data(snf + 558);
    const auto *snf_559 = buffer.data(snf + 559);

#pragma omp simd aligned(t_710, t_711, t_712, t_713, pb_x, pc_x, smg0_710, smf_475, smf_476, \
                         smf_477, smf_478, smg1_710, snf_476, snf_477, \
                         snf_478 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_710[k] = pb_x[k] * smg0_710[k]
                   + f_8 * smf_475[k]
                   - f_6 * pc_x[k] * smg1_710[k];

        t_711[k] = f_7 * smf_476[k]
                   + f_3 * pc_x[k] * snf_476[k];

        t_712[k] = f_7 * smf_477[k]
                   + f_3 * pc_x[k] * snf_477[k];

        t_713[k] = f_7 * smf_478[k]
                   + f_3 * pc_x[k] * snf_478[k];
    }

#pragma omp simd aligned(t_714, t_715, t_716, t_717, pb_x, pc_x, pc_z, smg0_715, smg0_717, \
                         smf_376, smf_479, smg1_715, smg1_717, snf_476, \
                         snf_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_714[k] = f_7 * smf_479[k]
                   + f_3 * pc_x[k] * snf_479[k];

        t_715[k] = pb_x[k] * smg0_715[k]
                   - f_6 * pc_x[k] * smg1_715[k];

        t_716[k] = f_8 * smf_376[k]
                   + f_3 * pc_z[k] * snf_476[k];

        t_717[k] = pb_x[k] * smg0_717[k]
                   - f_6 * pc_x[k] * smg1_717[k];
    }

#pragma omp simd aligned(t_718, t_719, t_720, t_721, pb_x, pc_x, pc_y, smg0_719, smg0_720, \
                         smf_389, smf_390, smf_480, smg1_719, smg1_720, snf_479, \
                         snf_480 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_718[k] = f_11 * smf_389[k]
                   + f_3 * pc_y[k] * snf_479[k];

        t_719[k] = pb_x[k] * smg0_719[k]
                   - f_6 * pc_x[k] * smg1_719[k];

        t_720[k] = pb_x[k] * smg0_720[k]
                   + f_14 * smf_480[k]
                   - f_6 * pc_x[k] * smg1_720[k];

        t_721[k] = f_13 * smf_390[k]
                   + f_3 * pc_y[k] * snf_480[k];
    }

#pragma omp simd aligned(t_722, t_723, t_724, pb_x, pc_x, pc_y, pc_z, smg0_723, smf_380, \
                         smf_392, smf_483, smg1_723, snf_480, snf_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_722[k] = f_12 * smf_380[k]
                   + f_3 * pc_z[k] * snf_480[k];

        t_723[k] = pb_x[k] * smg0_723[k]
                   + f_8 * smf_483[k]
                   - f_6 * pc_x[k] * smg1_723[k];

        t_724[k] = f_13 * smf_392[k]
                   + f_3 * pc_y[k] * snf_482[k];
    }

#pragma omp simd aligned(t_725, t_726, t_727, t_728, pb_x, pc_x, smg0_725, smf_485, smf_486, \
                         smf_487, smf_488, smg1_725, snf_486, snf_487, \
                         snf_488 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_725[k] = pb_x[k] * smg0_725[k]
                   + f_8 * smf_485[k]
                   - f_6 * pc_x[k] * smg1_725[k];

        t_726[k] = f_7 * smf_486[k]
                   + f_3 * pc_x[k] * snf_486[k];

        t_727[k] = f_7 * smf_487[k]
                   + f_3 * pc_x[k] * snf_487[k];

        t_728[k] = f_7 * smf_488[k]
                   + f_3 * pc_x[k] * snf_488[k];
    }

#pragma omp simd aligned(t_729, t_730, t_731, t_732, pb_x, pc_x, pc_z, smg0_730, smg0_732, \
                         smf_386, smf_489, smg1_730, smg1_732, snf_486, \
                         snf_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_729[k] = f_7 * smf_489[k]
                   + f_3 * pc_x[k] * snf_489[k];

        t_730[k] = pb_x[k] * smg0_730[k]
                   - f_6 * pc_x[k] * smg1_730[k];

        t_731[k] = f_12 * smf_386[k]
                   + f_3 * pc_z[k] * snf_486[k];

        t_732[k] = pb_x[k] * smg0_732[k]
                   - f_6 * pc_x[k] * smg1_732[k];
    }

#pragma omp simd aligned(t_733, t_734, t_735, t_736, pb_x, pc_x, pc_y, smg0_734, smg0_735, \
                         smf_399, smf_400, smf_490, smg1_734, smg1_735, snf_489, \
                         snf_490 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_733[k] = f_13 * smf_399[k]
                   + f_3 * pc_y[k] * snf_489[k];

        t_734[k] = pb_x[k] * smg0_734[k]
                   - f_6 * pc_x[k] * smg1_734[k];

        t_735[k] = pb_x[k] * smg0_735[k]
                   + f_14 * smf_490[k]
                   - f_6 * pc_x[k] * smg1_735[k];

        t_736[k] = f_15 * smf_400[k]
                   + f_3 * pc_y[k] * snf_490[k];
    }

#pragma omp simd aligned(t_737, t_738, t_739, pb_x, pc_x, pc_y, pc_z, smg0_738, smf_390, \
                         smf_402, smf_493, smg1_738, snf_490, snf_492 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_737[k] = f_14 * smf_390[k]
                   + f_3 * pc_z[k] * snf_490[k];

        t_738[k] = pb_x[k] * smg0_738[k]
                   + f_8 * smf_493[k]
                   - f_6 * pc_x[k] * smg1_738[k];

        t_739[k] = f_15 * smf_402[k]
                   + f_3 * pc_y[k] * snf_492[k];
    }

#pragma omp simd aligned(t_740, t_741, t_742, t_743, pb_x, pc_x, smg0_740, smf_495, smf_496, \
                         smf_497, smf_498, smg1_740, snf_496, snf_497, \
                         snf_498 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_740[k] = pb_x[k] * smg0_740[k]
                   + f_8 * smf_495[k]
                   - f_6 * pc_x[k] * smg1_740[k];

        t_741[k] = f_7 * smf_496[k]
                   + f_3 * pc_x[k] * snf_496[k];

        t_742[k] = f_7 * smf_497[k]
                   + f_3 * pc_x[k] * snf_497[k];

        t_743[k] = f_7 * smf_498[k]
                   + f_3 * pc_x[k] * snf_498[k];
    }

#pragma omp simd aligned(t_744, t_745, t_746, t_747, pb_x, pc_x, pc_z, smg0_745, smg0_747, \
                         smf_396, smf_499, smg1_745, smg1_747, snf_496, \
                         snf_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_744[k] = f_7 * smf_499[k]
                   + f_3 * pc_x[k] * snf_499[k];

        t_745[k] = pb_x[k] * smg0_745[k]
                   - f_6 * pc_x[k] * smg1_745[k];

        t_746[k] = f_14 * smf_396[k]
                   + f_3 * pc_z[k] * snf_496[k];

        t_747[k] = pb_x[k] * smg0_747[k]
                   - f_6 * pc_x[k] * smg1_747[k];
    }

#pragma omp simd aligned(t_748, t_749, t_750, t_751, pb_x, pc_x, pc_y, smg0_749, smg0_750, \
                         smf_409, smf_410, smf_500, smg1_749, smg1_750, snf_499, \
                         snf_500 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_748[k] = f_15 * smf_409[k]
                   + f_3 * pc_y[k] * snf_499[k];

        t_749[k] = pb_x[k] * smg0_749[k]
                   - f_6 * pc_x[k] * smg1_749[k];

        t_750[k] = pb_x[k] * smg0_750[k]
                   + f_14 * smf_500[k]
                   - f_6 * pc_x[k] * smg1_750[k];

        t_751[k] = f_14 * smf_410[k]
                   + f_3 * pc_y[k] * snf_500[k];
    }

#pragma omp simd aligned(t_752, t_753, t_754, pb_x, pc_x, pc_y, pc_z, smg0_753, smf_400, \
                         smf_412, smf_503, smg1_753, snf_500, snf_502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_752[k] = f_15 * smf_400[k]
                   + f_3 * pc_z[k] * snf_500[k];

        t_753[k] = pb_x[k] * smg0_753[k]
                   + f_8 * smf_503[k]
                   - f_6 * pc_x[k] * smg1_753[k];

        t_754[k] = f_14 * smf_412[k]
                   + f_3 * pc_y[k] * snf_502[k];
    }

#pragma omp simd aligned(t_755, t_756, t_757, t_758, pb_x, pc_x, smg0_755, smf_505, smf_506, \
                         smf_507, smf_508, smg1_755, snf_506, snf_507, \
                         snf_508 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_755[k] = pb_x[k] * smg0_755[k]
                   + f_8 * smf_505[k]
                   - f_6 * pc_x[k] * smg1_755[k];

        t_756[k] = f_7 * smf_506[k]
                   + f_3 * pc_x[k] * snf_506[k];

        t_757[k] = f_7 * smf_507[k]
                   + f_3 * pc_x[k] * snf_507[k];

        t_758[k] = f_7 * smf_508[k]
                   + f_3 * pc_x[k] * snf_508[k];
    }

#pragma omp simd aligned(t_759, t_760, t_761, t_762, pb_x, pc_x, pc_z, smg0_760, smg0_762, \
                         smf_406, smf_509, smg1_760, smg1_762, snf_506, \
                         snf_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_759[k] = f_7 * smf_509[k]
                   + f_3 * pc_x[k] * snf_509[k];

        t_760[k] = pb_x[k] * smg0_760[k]
                   - f_6 * pc_x[k] * smg1_760[k];

        t_761[k] = f_15 * smf_406[k]
                   + f_3 * pc_z[k] * snf_506[k];

        t_762[k] = pb_x[k] * smg0_762[k]
                   - f_6 * pc_x[k] * smg1_762[k];
    }

#pragma omp simd aligned(t_763, t_764, t_765, t_766, pb_x, pc_x, pc_y, smg0_764, smg0_765, \
                         smf_419, smf_420, smf_510, smg1_764, smg1_765, snf_509, \
                         snf_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_763[k] = f_14 * smf_419[k]
                   + f_3 * pc_y[k] * snf_509[k];

        t_764[k] = pb_x[k] * smg0_764[k]
                   - f_6 * pc_x[k] * smg1_764[k];

        t_765[k] = pb_x[k] * smg0_765[k]
                   + f_14 * smf_510[k]
                   - f_6 * pc_x[k] * smg1_765[k];

        t_766[k] = f_12 * smf_420[k]
                   + f_3 * pc_y[k] * snf_510[k];
    }

#pragma omp simd aligned(t_767, t_768, t_769, pb_x, pc_x, pc_y, pc_z, smg0_768, smf_410, \
                         smf_422, smf_513, smg1_768, snf_510, snf_512 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_767[k] = f_13 * smf_410[k]
                   + f_3 * pc_z[k] * snf_510[k];

        t_768[k] = pb_x[k] * smg0_768[k]
                   + f_8 * smf_513[k]
                   - f_6 * pc_x[k] * smg1_768[k];

        t_769[k] = f_12 * smf_422[k]
                   + f_3 * pc_y[k] * snf_512[k];
    }

#pragma omp simd aligned(t_770, t_771, t_772, t_773, pb_x, pc_x, smg0_770, smf_515, smf_516, \
                         smf_517, smf_518, smg1_770, snf_516, snf_517, \
                         snf_518 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_770[k] = pb_x[k] * smg0_770[k]
                   + f_8 * smf_515[k]
                   - f_6 * pc_x[k] * smg1_770[k];

        t_771[k] = f_7 * smf_516[k]
                   + f_3 * pc_x[k] * snf_516[k];

        t_772[k] = f_7 * smf_517[k]
                   + f_3 * pc_x[k] * snf_517[k];

        t_773[k] = f_7 * smf_518[k]
                   + f_3 * pc_x[k] * snf_518[k];
    }

#pragma omp simd aligned(t_774, t_775, t_776, t_777, pb_x, pc_x, pc_z, smg0_775, smg0_777, \
                         smf_416, smf_519, smg1_775, smg1_777, snf_516, \
                         snf_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_774[k] = f_7 * smf_519[k]
                   + f_3 * pc_x[k] * snf_519[k];

        t_775[k] = pb_x[k] * smg0_775[k]
                   - f_6 * pc_x[k] * smg1_775[k];

        t_776[k] = f_13 * smf_416[k]
                   + f_3 * pc_z[k] * snf_516[k];

        t_777[k] = pb_x[k] * smg0_777[k]
                   - f_6 * pc_x[k] * smg1_777[k];
    }

#pragma omp simd aligned(t_778, t_779, t_780, t_781, pb_x, pc_x, pc_y, smg0_779, smg0_780, \
                         smf_429, smf_430, smf_520, smg1_779, smg1_780, snf_519, \
                         snf_520 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_778[k] = f_12 * smf_429[k]
                   + f_3 * pc_y[k] * snf_519[k];

        t_779[k] = pb_x[k] * smg0_779[k]
                   - f_6 * pc_x[k] * smg1_779[k];

        t_780[k] = pb_x[k] * smg0_780[k]
                   + f_14 * smf_520[k]
                   - f_6 * pc_x[k] * smg1_780[k];

        t_781[k] = f_8 * smf_430[k]
                   + f_3 * pc_y[k] * snf_520[k];
    }

#pragma omp simd aligned(t_782, t_783, t_784, pb_x, pc_x, pc_y, pc_z, smg0_783, smf_420, \
                         smf_432, smf_523, smg1_783, snf_520, snf_522 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_782[k] = f_11 * smf_420[k]
                   + f_3 * pc_z[k] * snf_520[k];

        t_783[k] = pb_x[k] * smg0_783[k]
                   + f_8 * smf_523[k]
                   - f_6 * pc_x[k] * smg1_783[k];

        t_784[k] = f_8 * smf_432[k]
                   + f_3 * pc_y[k] * snf_522[k];
    }

#pragma omp simd aligned(t_785, t_786, t_787, t_788, pb_x, pc_x, smg0_785, smf_525, smf_526, \
                         smf_527, smf_528, smg1_785, snf_526, snf_527, \
                         snf_528 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_785[k] = pb_x[k] * smg0_785[k]
                   + f_8 * smf_525[k]
                   - f_6 * pc_x[k] * smg1_785[k];

        t_786[k] = f_7 * smf_526[k]
                   + f_3 * pc_x[k] * snf_526[k];

        t_787[k] = f_7 * smf_527[k]
                   + f_3 * pc_x[k] * snf_527[k];

        t_788[k] = f_7 * smf_528[k]
                   + f_3 * pc_x[k] * snf_528[k];
    }

#pragma omp simd aligned(t_789, t_790, t_791, t_792, pb_x, pc_x, pc_z, smg0_790, smg0_792, \
                         smf_426, smf_529, smg1_790, smg1_792, snf_526, \
                         snf_529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_789[k] = f_7 * smf_529[k]
                   + f_3 * pc_x[k] * snf_529[k];

        t_790[k] = pb_x[k] * smg0_790[k]
                   - f_6 * pc_x[k] * smg1_790[k];

        t_791[k] = f_11 * smf_426[k]
                   + f_3 * pc_z[k] * snf_526[k];

        t_792[k] = pb_x[k] * smg0_792[k]
                   - f_6 * pc_x[k] * smg1_792[k];
    }

#pragma omp simd aligned(t_793, t_794, t_795, t_796, pb_x, pb_y, pc_x, pc_y, smg0_660, \
                         smg0_794, smf_439, smf_440, smg1_660, smg1_794, snf_529, \
                         snf_530 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_793[k] = f_8 * smf_439[k]
                   + f_3 * pc_y[k] * snf_529[k];

        t_794[k] = pb_x[k] * smg0_794[k]
                   - f_6 * pc_x[k] * smg1_794[k];

        t_795[k] = pb_y[k] * smg0_660[k]
                   - f_6 * pc_y[k] * smg1_660[k];

        t_796[k] = f_7 * smf_440[k]
                   + f_3 * pc_y[k] * snf_530[k];
    }

#pragma omp simd aligned(t_797, t_798, t_799, pb_x, pc_x, pc_y, pc_z, smg0_798, smf_430, \
                         smf_442, smf_533, smg1_798, snf_530, snf_532 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_797[k] = f_10 * smf_430[k]
                   + f_3 * pc_z[k] * snf_530[k];

        t_798[k] = pb_x[k] * smg0_798[k]
                   + f_8 * smf_533[k]
                   - f_6 * pc_x[k] * smg1_798[k];

        t_799[k] = f_7 * smf_442[k]
                   + f_3 * pc_y[k] * snf_532[k];
    }

#pragma omp simd aligned(t_800, t_801, t_802, t_803, pb_y, pc_x, pc_y, smg0_665, smf_536, \
                         smf_537, smf_538, smg1_665, snf_536, snf_537, \
                         snf_538 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_800[k] = pb_y[k] * smg0_665[k]
                   - f_6 * pc_y[k] * smg1_665[k];

        t_801[k] = f_7 * smf_536[k]
                   + f_3 * pc_x[k] * snf_536[k];

        t_802[k] = f_7 * smf_537[k]
                   + f_3 * pc_x[k] * snf_537[k];

        t_803[k] = f_7 * smf_538[k]
                   + f_3 * pc_x[k] * snf_538[k];
    }

#pragma omp simd aligned(t_804, t_805, t_806, t_807, pb_x, pc_x, pc_z, smg0_805, smg0_807, \
                         smf_436, smf_539, smg1_805, smg1_807, snf_536, \
                         snf_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_804[k] = f_7 * smf_539[k]
                   + f_3 * pc_x[k] * snf_539[k];

        t_805[k] = pb_x[k] * smg0_805[k]
                   - f_6 * pc_x[k] * smg1_805[k];

        t_806[k] = f_10 * smf_436[k]
                   + f_3 * pc_z[k] * snf_536[k];

        t_807[k] = pb_x[k] * smg0_807[k]
                   - f_6 * pc_x[k] * smg1_807[k];
    }

#pragma omp simd aligned(t_808, t_809, t_810, t_811, pb_x, pc_x, pc_y, smg0_809, smg0_810, \
                         smf_449, smf_540, smg1_809, smg1_810, snf_539, \
                         snf_540 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_808[k] = f_7 * smf_449[k]
                   + f_3 * pc_y[k] * snf_539[k];

        t_809[k] = pb_x[k] * smg0_809[k]
                   - f_6 * pc_x[k] * smg1_809[k];

        t_810[k] = pb_x[k] * smg0_810[k]
                   + f_14 * smf_540[k]
                   - f_6 * pc_x[k] * smg1_810[k];

        t_811[k] = f_3 * pc_y[k] * snf_540[k];
    }

#pragma omp simd aligned(t_812, t_813, t_814, pb_x, pc_x, pc_y, pc_z, smg0_813, smf_440, \
                         smf_543, smg1_813, snf_540, snf_542 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_812[k] = f_9 * smf_440[k]
                   + f_3 * pc_z[k] * snf_540[k];

        t_813[k] = pb_x[k] * smg0_813[k]
                   + f_8 * smf_543[k]
                   - f_6 * pc_x[k] * smg1_813[k];

        t_814[k] = f_3 * pc_y[k] * snf_542[k];
    }

#pragma omp simd aligned(t_815, t_816, t_817, t_818, pb_x, pc_x, smg0_815, smf_545, smf_546, \
                         smf_547, smf_548, smg1_815, snf_546, snf_547, \
                         snf_548 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_815[k] = pb_x[k] * smg0_815[k]
                   + f_8 * smf_545[k]
                   - f_6 * pc_x[k] * smg1_815[k];

        t_816[k] = f_7 * smf_546[k]
                   + f_3 * pc_x[k] * snf_546[k];

        t_817[k] = f_7 * smf_547[k]
                   + f_3 * pc_x[k] * snf_547[k];

        t_818[k] = f_7 * smf_548[k]
                   + f_3 * pc_x[k] * snf_548[k];
    }

#pragma omp simd aligned(t_819, t_820, t_821, t_822, pb_x, pc_x, pc_z, smg0_820, smg0_822, \
                         smf_446, smf_549, smg1_820, smg1_822, snf_546, \
                         snf_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_819[k] = f_7 * smf_549[k]
                   + f_3 * pc_x[k] * snf_549[k];

        t_820[k] = pb_x[k] * smg0_820[k]
                   - f_6 * pc_x[k] * smg1_820[k];

        t_821[k] = f_9 * smf_446[k]
                   + f_3 * pc_z[k] * snf_546[k];

        t_822[k] = pb_x[k] * smg0_822[k]
                   - f_6 * pc_x[k] * smg1_822[k];
    }

#pragma omp simd aligned(t_823, t_824, t_825, t_826, t_827, pb_x, pc_x, pc_y, pc_z, smg0_824, \
                         smf_450, smg1_824, snd0_330, snd1_330, snf_549, \
                         snf_550 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_823[k] = f_3 * pc_y[k] * snf_549[k];

        t_824[k] = pb_x[k] * smg0_824[k]
                   - f_6 * pc_x[k] * smg1_824[k];

        t_825[k] = f_1 * snd0_330[k]
                   - f_2 * snd1_330[k]
                   + f_3 * pc_x[k] * snf_550[k];

        t_826[k] = f_0 * smf_450[k]
                   + f_3 * pc_y[k] * snf_550[k];

        t_827[k] = f_3 * pc_z[k] * snf_550[k];
    }

#pragma omp simd aligned(t_828, t_829, t_830, t_831, pc_x, pc_y, smf_452, snd0_333, snd0_335, \
                         snd1_333, snd1_335, snf_552, snf_553, snf_555, \
                         snf_556 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_828[k] = f_4 * snd0_333[k]
                   - f_5 * snd1_333[k]
                   + f_3 * pc_x[k] * snf_553[k];

        t_829[k] = f_0 * smf_452[k]
                   + f_3 * pc_y[k] * snf_552[k];

        t_830[k] = f_4 * snd0_335[k]
                   - f_5 * snd1_335[k]
                   + f_3 * pc_x[k] * snf_555[k];

        t_831[k] = f_3 * pc_x[k] * snf_556[k];
    }

#pragma omp simd aligned(t_832, t_833, t_834, t_835, t_836, pc_x, pc_y, pc_z, smf_456, \
                         snd0_333, snd1_333, snf_556, snf_557, snf_558, \
                         snf_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_832[k] = f_3 * pc_x[k] * snf_557[k];

        t_833[k] = f_3 * pc_x[k] * snf_558[k];

        t_834[k] = f_3 * pc_x[k] * snf_559[k];

        t_835[k] = f_0 * smf_456[k]
                   + f_1 * snd0_333[k]
                   - f_2 * snd1_333[k]
                   + f_3 * pc_y[k] * snf_556[k];

        t_836[k] = f_3 * pc_z[k] * snf_556[k];
    }
}

static auto
compute_prim_sng_three_center_electron_repulsion_0_piece7(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t smg0,
                                                          const size_t smf, const size_t smg1,
                                                          const size_t snd0, const size_t snd1,
                                                          const size_t snf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.0 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 4.5 / q;
    const auto f_10 = 4.0 / q;
    const auto f_11 = 3.5 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 3.0 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 2.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *smg0_675 = buffer.data(smg0 + 675);
    const auto *smg0_678 = buffer.data(smg0 + 678);
    const auto *smg0_685 = buffer.data(smg0 + 685);
    const auto *smg0_687 = buffer.data(smg0 + 687);
    const auto *smg0_810 = buffer.data(smg0 + 810);

    const auto *smf_450 = buffer.data(smf + 450);
    const auto *smf_456 = buffer.data(smf + 456);
    const auto *smf_457 = buffer.data(smf + 457);
    const auto *smf_458 = buffer.data(smf + 458);
    const auto *smf_459 = buffer.data(smf + 459);
    const auto *smf_460 = buffer.data(smf + 460);
    const auto *smf_462 = buffer.data(smf + 462);
    const auto *smf_466 = buffer.data(smf + 466);
    const auto *smf_469 = buffer.data(smf + 469);
    const auto *smf_470 = buffer.data(smf + 470);
    const auto *smf_472 = buffer.data(smf + 472);
    const auto *smf_476 = buffer.data(smf + 476);
    const auto *smf_478 = buffer.data(smf + 478);
    const auto *smf_479 = buffer.data(smf + 479);
    const auto *smf_480 = buffer.data(smf + 480);
    const auto *smf_482 = buffer.data(smf + 482);
    const auto *smf_486 = buffer.data(smf + 486);
    const auto *smf_488 = buffer.data(smf + 488);
    const auto *smf_489 = buffer.data(smf + 489);
    const auto *smf_490 = buffer.data(smf + 490);
    const auto *smf_492 = buffer.data(smf + 492);
    const auto *smf_496 = buffer.data(smf + 496);
    const auto *smf_498 = buffer.data(smf + 498);
    const auto *smf_499 = buffer.data(smf + 499);
    const auto *smf_500 = buffer.data(smf + 500);
    const auto *smf_502 = buffer.data(smf + 502);
    const auto *smf_506 = buffer.data(smf + 506);
    const auto *smf_508 = buffer.data(smf + 508);
    const auto *smf_509 = buffer.data(smf + 509);
    const auto *smf_510 = buffer.data(smf + 510);
    const auto *smf_512 = buffer.data(smf + 512);
    const auto *smf_516 = buffer.data(smf + 516);
    const auto *smf_518 = buffer.data(smf + 518);
    const auto *smf_519 = buffer.data(smf + 519);
    const auto *smf_520 = buffer.data(smf + 520);
    const auto *smf_522 = buffer.data(smf + 522);
    const auto *smf_526 = buffer.data(smf + 526);
    const auto *smf_528 = buffer.data(smf + 528);
    const auto *smf_529 = buffer.data(smf + 529);
    const auto *smf_530 = buffer.data(smf + 530);
    const auto *smf_532 = buffer.data(smf + 532);
    const auto *smf_536 = buffer.data(smf + 536);
    const auto *smf_538 = buffer.data(smf + 538);
    const auto *smf_539 = buffer.data(smf + 539);
    const auto *smf_540 = buffer.data(smf + 540);
    const auto *smf_542 = buffer.data(smf + 542);

    const auto *smg1_675 = buffer.data(smg1 + 675);
    const auto *smg1_678 = buffer.data(smg1 + 678);
    const auto *smg1_685 = buffer.data(smg1 + 685);
    const auto *smg1_687 = buffer.data(smg1 + 687);
    const auto *smg1_810 = buffer.data(smg1 + 810);

    const auto *snd0_335 = buffer.data(snd0 + 335);
    const auto *snd0_341 = buffer.data(snd0 + 341);
    const auto *snd0_342 = buffer.data(snd0 + 342);
    const auto *snd0_345 = buffer.data(snd0 + 345);
    const auto *snd0_347 = buffer.data(snd0 + 347);
    const auto *snd0_348 = buffer.data(snd0 + 348);
    const auto *snd0_351 = buffer.data(snd0 + 351);
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

    const auto *snd1_335 = buffer.data(snd1 + 335);
    const auto *snd1_341 = buffer.data(snd1 + 341);
    const auto *snd1_342 = buffer.data(snd1 + 342);
    const auto *snd1_345 = buffer.data(snd1 + 345);
    const auto *snd1_347 = buffer.data(snd1 + 347);
    const auto *snd1_348 = buffer.data(snd1 + 348);
    const auto *snd1_351 = buffer.data(snd1 + 351);
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

    const auto *snf_558 = buffer.data(snf + 558);
    const auto *snf_559 = buffer.data(snf + 559);
    const auto *snf_560 = buffer.data(snf + 560);
    const auto *snf_562 = buffer.data(snf + 562);
    const auto *snf_565 = buffer.data(snf + 565);
    const auto *snf_566 = buffer.data(snf + 566);
    const auto *snf_567 = buffer.data(snf + 567);
    const auto *snf_568 = buffer.data(snf + 568);
    const auto *snf_569 = buffer.data(snf + 569);
    const auto *snf_570 = buffer.data(snf + 570);
    const auto *snf_572 = buffer.data(snf + 572);
    const auto *snf_573 = buffer.data(snf + 573);
    const auto *snf_575 = buffer.data(snf + 575);
    const auto *snf_576 = buffer.data(snf + 576);
    const auto *snf_577 = buffer.data(snf + 577);
    const auto *snf_578 = buffer.data(snf + 578);
    const auto *snf_579 = buffer.data(snf + 579);
    const auto *snf_580 = buffer.data(snf + 580);
    const auto *snf_582 = buffer.data(snf + 582);
    const auto *snf_583 = buffer.data(snf + 583);
    const auto *snf_585 = buffer.data(snf + 585);
    const auto *snf_586 = buffer.data(snf + 586);
    const auto *snf_587 = buffer.data(snf + 587);
    const auto *snf_588 = buffer.data(snf + 588);
    const auto *snf_589 = buffer.data(snf + 589);
    const auto *snf_590 = buffer.data(snf + 590);
    const auto *snf_592 = buffer.data(snf + 592);
    const auto *snf_593 = buffer.data(snf + 593);
    const auto *snf_595 = buffer.data(snf + 595);
    const auto *snf_596 = buffer.data(snf + 596);
    const auto *snf_597 = buffer.data(snf + 597);
    const auto *snf_598 = buffer.data(snf + 598);
    const auto *snf_599 = buffer.data(snf + 599);
    const auto *snf_600 = buffer.data(snf + 600);
    const auto *snf_602 = buffer.data(snf + 602);
    const auto *snf_603 = buffer.data(snf + 603);
    const auto *snf_605 = buffer.data(snf + 605);
    const auto *snf_606 = buffer.data(snf + 606);
    const auto *snf_607 = buffer.data(snf + 607);
    const auto *snf_608 = buffer.data(snf + 608);
    const auto *snf_609 = buffer.data(snf + 609);
    const auto *snf_610 = buffer.data(snf + 610);
    const auto *snf_612 = buffer.data(snf + 612);
    const auto *snf_613 = buffer.data(snf + 613);
    const auto *snf_615 = buffer.data(snf + 615);
    const auto *snf_616 = buffer.data(snf + 616);
    const auto *snf_617 = buffer.data(snf + 617);
    const auto *snf_618 = buffer.data(snf + 618);
    const auto *snf_619 = buffer.data(snf + 619);
    const auto *snf_620 = buffer.data(snf + 620);
    const auto *snf_622 = buffer.data(snf + 622);
    const auto *snf_623 = buffer.data(snf + 623);
    const auto *snf_625 = buffer.data(snf + 625);
    const auto *snf_626 = buffer.data(snf + 626);
    const auto *snf_627 = buffer.data(snf + 627);
    const auto *snf_628 = buffer.data(snf + 628);
    const auto *snf_629 = buffer.data(snf + 629);
    const auto *snf_630 = buffer.data(snf + 630);
    const auto *snf_632 = buffer.data(snf + 632);
    const auto *snf_633 = buffer.data(snf + 633);
    const auto *snf_635 = buffer.data(snf + 635);
    const auto *snf_636 = buffer.data(snf + 636);
    const auto *snf_637 = buffer.data(snf + 637);
    const auto *snf_638 = buffer.data(snf + 638);
    const auto *snf_639 = buffer.data(snf + 639);
    const auto *snf_640 = buffer.data(snf + 640);
    const auto *snf_642 = buffer.data(snf + 642);
    const auto *snf_643 = buffer.data(snf + 643);

#pragma omp simd aligned(t_837, t_838, t_839, t_840, pb_z, pc_y, pc_z, smg0_675, smf_458, \
                         smf_459, smg1_675, snd0_335, snd1_335, snf_558, \
                         snf_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_837[k] = f_0 * smf_458[k]
                   + f_4 * snd0_335[k]
                   - f_5 * snd1_335[k]
                   + f_3 * pc_y[k] * snf_558[k];

        t_838[k] = f_0 * smf_459[k]
                   + f_3 * pc_y[k] * snf_559[k];

        t_839[k] = f_1 * snd0_335[k]
                   - f_2 * snd1_335[k]
                   + f_3 * pc_z[k] * snf_559[k];

        t_840[k] = pb_z[k] * smg0_675[k]
                   - f_6 * pc_z[k] * smg1_675[k];
    }

#pragma omp simd aligned(t_841, t_842, t_843, t_844, pb_z, pc_y, pc_z, smg0_678, smf_450, \
                         smf_460, smf_462, smg1_678, snf_560, snf_562 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_841[k] = f_9 * smf_460[k]
                   + f_3 * pc_y[k] * snf_560[k];

        t_842[k] = f_7 * smf_450[k]
                   + f_3 * pc_z[k] * snf_560[k];

        t_843[k] = pb_z[k] * smg0_678[k]
                   - f_6 * pc_z[k] * smg1_678[k];

        t_844[k] = f_9 * smf_462[k]
                   + f_3 * pc_y[k] * snf_562[k];
    }

#pragma omp simd aligned(t_845, t_846, t_847, t_848, t_849, pc_x, snd0_341, snd1_341, snf_565, \
                         snf_566, snf_567, snf_568, snf_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_845[k] = f_4 * snd0_341[k]
                   - f_5 * snd1_341[k]
                   + f_3 * pc_x[k] * snf_565[k];

        t_846[k] = f_3 * pc_x[k] * snf_566[k];

        t_847[k] = f_3 * pc_x[k] * snf_567[k];

        t_848[k] = f_3 * pc_x[k] * snf_568[k];

        t_849[k] = f_3 * pc_x[k] * snf_569[k];
    }

#pragma omp simd aligned(t_850, t_851, t_852, t_853, pb_z, pc_y, pc_z, smg0_685, smg0_687, \
                         smf_456, smf_457, smf_469, smg1_685, smg1_687, snf_566, \
                         snf_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_850[k] = pb_z[k] * smg0_685[k]
                   - f_6 * pc_z[k] * smg1_685[k];

        t_851[k] = f_7 * smf_456[k]
                   + f_3 * pc_z[k] * snf_566[k];

        t_852[k] = pb_z[k] * smg0_687[k]
                   + f_8 * smf_457[k]
                   - f_6 * pc_z[k] * smg1_687[k];

        t_853[k] = f_9 * smf_469[k]
                   + f_3 * pc_y[k] * snf_569[k];
    }

#pragma omp simd aligned(t_854, t_855, t_856, t_857, pc_x, pc_y, pc_z, smf_459, smf_460, \
                         smf_470, snd0_341, snd0_342, snd1_341, snd1_342, snf_569, \
                         snf_570 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_854[k] = f_7 * smf_459[k]
                   + f_1 * snd0_341[k]
                   - f_2 * snd1_341[k]
                   + f_3 * pc_z[k] * snf_569[k];

        t_855[k] = f_1 * snd0_342[k]
                   - f_2 * snd1_342[k]
                   + f_3 * pc_x[k] * snf_570[k];

        t_856[k] = f_10 * smf_470[k]
                   + f_3 * pc_y[k] * snf_570[k];

        t_857[k] = f_8 * smf_460[k]
                   + f_3 * pc_z[k] * snf_570[k];
    }

#pragma omp simd aligned(t_858, t_859, t_860, t_861, pc_x, pc_y, smf_472, snd0_345, snd0_347, \
                         snd1_345, snd1_347, snf_572, snf_573, snf_575, \
                         snf_576 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_858[k] = f_4 * snd0_345[k]
                   - f_5 * snd1_345[k]
                   + f_3 * pc_x[k] * snf_573[k];

        t_859[k] = f_10 * smf_472[k]
                   + f_3 * pc_y[k] * snf_572[k];

        t_860[k] = f_4 * snd0_347[k]
                   - f_5 * snd1_347[k]
                   + f_3 * pc_x[k] * snf_575[k];

        t_861[k] = f_3 * pc_x[k] * snf_576[k];
    }

#pragma omp simd aligned(t_862, t_863, t_864, t_865, t_866, pc_x, pc_y, pc_z, smf_466, \
                         smf_476, snd0_345, snd1_345, snf_576, snf_577, snf_578, \
                         snf_579 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_862[k] = f_3 * pc_x[k] * snf_577[k];

        t_863[k] = f_3 * pc_x[k] * snf_578[k];

        t_864[k] = f_3 * pc_x[k] * snf_579[k];

        t_865[k] = f_10 * smf_476[k]
                   + f_1 * snd0_345[k]
                   - f_2 * snd1_345[k]
                   + f_3 * pc_y[k] * snf_576[k];

        t_866[k] = f_8 * smf_466[k]
                   + f_3 * pc_z[k] * snf_576[k];
    }

#pragma omp simd aligned(t_867, t_868, t_869, pc_y, pc_z, smf_469, smf_478, smf_479, snd0_347, \
                         snd1_347, snf_578, snf_579 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_867[k] = f_10 * smf_478[k]
                   + f_4 * snd0_347[k]
                   - f_5 * snd1_347[k]
                   + f_3 * pc_y[k] * snf_578[k];

        t_868[k] = f_10 * smf_479[k]
                   + f_3 * pc_y[k] * snf_579[k];

        t_869[k] = f_8 * smf_469[k]
                   + f_1 * snd0_347[k]
                   - f_2 * snd1_347[k]
                   + f_3 * pc_z[k] * snf_579[k];
    }

#pragma omp simd aligned(t_870, t_871, t_872, t_873, pc_x, pc_y, pc_z, smf_470, smf_480, \
                         snd0_348, snd0_351, snd1_348, snd1_351, snf_580, \
                         snf_583 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_870[k] = f_1 * snd0_348[k]
                   - f_2 * snd1_348[k]
                   + f_3 * pc_x[k] * snf_580[k];

        t_871[k] = f_11 * smf_480[k]
                   + f_3 * pc_y[k] * snf_580[k];

        t_872[k] = f_12 * smf_470[k]
                   + f_3 * pc_z[k] * snf_580[k];

        t_873[k] = f_4 * snd0_351[k]
                   - f_5 * snd1_351[k]
                   + f_3 * pc_x[k] * snf_583[k];
    }

#pragma omp simd aligned(t_874, t_875, t_876, t_877, t_878, pc_x, pc_y, smf_482, snd0_353, \
                         snd1_353, snf_582, snf_585, snf_586, snf_587, \
                         snf_588 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_874[k] = f_11 * smf_482[k]
                   + f_3 * pc_y[k] * snf_582[k];

        t_875[k] = f_4 * snd0_353[k]
                   - f_5 * snd1_353[k]
                   + f_3 * pc_x[k] * snf_585[k];

        t_876[k] = f_3 * pc_x[k] * snf_586[k];

        t_877[k] = f_3 * pc_x[k] * snf_587[k];

        t_878[k] = f_3 * pc_x[k] * snf_588[k];
    }

#pragma omp simd aligned(t_879, t_880, t_881, pc_x, pc_y, pc_z, smf_476, smf_486, snd0_351, \
                         snd1_351, snf_586, snf_589 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_879[k] = f_3 * pc_x[k] * snf_589[k];

        t_880[k] = f_11 * smf_486[k]
                   + f_1 * snd0_351[k]
                   - f_2 * snd1_351[k]
                   + f_3 * pc_y[k] * snf_586[k];

        t_881[k] = f_12 * smf_476[k]
                   + f_3 * pc_z[k] * snf_586[k];
    }

#pragma omp simd aligned(t_882, t_883, t_884, pc_y, pc_z, smf_479, smf_488, smf_489, snd0_353, \
                         snd1_353, snf_588, snf_589 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_882[k] = f_11 * smf_488[k]
                   + f_4 * snd0_353[k]
                   - f_5 * snd1_353[k]
                   + f_3 * pc_y[k] * snf_588[k];

        t_883[k] = f_11 * smf_489[k]
                   + f_3 * pc_y[k] * snf_589[k];

        t_884[k] = f_12 * smf_479[k]
                   + f_1 * snd0_353[k]
                   - f_2 * snd1_353[k]
                   + f_3 * pc_z[k] * snf_589[k];
    }

#pragma omp simd aligned(t_885, t_886, t_887, t_888, pc_x, pc_y, pc_z, smf_480, smf_490, \
                         snd0_354, snd0_357, snd1_354, snd1_357, snf_590, \
                         snf_593 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_885[k] = f_1 * snd0_354[k]
                   - f_2 * snd1_354[k]
                   + f_3 * pc_x[k] * snf_590[k];

        t_886[k] = f_13 * smf_490[k]
                   + f_3 * pc_y[k] * snf_590[k];

        t_887[k] = f_14 * smf_480[k]
                   + f_3 * pc_z[k] * snf_590[k];

        t_888[k] = f_4 * snd0_357[k]
                   - f_5 * snd1_357[k]
                   + f_3 * pc_x[k] * snf_593[k];
    }

#pragma omp simd aligned(t_889, t_890, t_891, t_892, t_893, pc_x, pc_y, smf_492, snd0_359, \
                         snd1_359, snf_592, snf_595, snf_596, snf_597, \
                         snf_598 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_889[k] = f_13 * smf_492[k]
                   + f_3 * pc_y[k] * snf_592[k];

        t_890[k] = f_4 * snd0_359[k]
                   - f_5 * snd1_359[k]
                   + f_3 * pc_x[k] * snf_595[k];

        t_891[k] = f_3 * pc_x[k] * snf_596[k];

        t_892[k] = f_3 * pc_x[k] * snf_597[k];

        t_893[k] = f_3 * pc_x[k] * snf_598[k];
    }

#pragma omp simd aligned(t_894, t_895, t_896, pc_x, pc_y, pc_z, smf_486, smf_496, snd0_357, \
                         snd1_357, snf_596, snf_599 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_894[k] = f_3 * pc_x[k] * snf_599[k];

        t_895[k] = f_13 * smf_496[k]
                   + f_1 * snd0_357[k]
                   - f_2 * snd1_357[k]
                   + f_3 * pc_y[k] * snf_596[k];

        t_896[k] = f_14 * smf_486[k]
                   + f_3 * pc_z[k] * snf_596[k];
    }

#pragma omp simd aligned(t_897, t_898, t_899, pc_y, pc_z, smf_489, smf_498, smf_499, snd0_359, \
                         snd1_359, snf_598, snf_599 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_897[k] = f_13 * smf_498[k]
                   + f_4 * snd0_359[k]
                   - f_5 * snd1_359[k]
                   + f_3 * pc_y[k] * snf_598[k];

        t_898[k] = f_13 * smf_499[k]
                   + f_3 * pc_y[k] * snf_599[k];

        t_899[k] = f_14 * smf_489[k]
                   + f_1 * snd0_359[k]
                   - f_2 * snd1_359[k]
                   + f_3 * pc_z[k] * snf_599[k];
    }

#pragma omp simd aligned(t_900, t_901, t_902, t_903, pc_x, pc_y, pc_z, smf_490, smf_500, \
                         snd0_360, snd0_363, snd1_360, snd1_363, snf_600, \
                         snf_603 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_900[k] = f_1 * snd0_360[k]
                   - f_2 * snd1_360[k]
                   + f_3 * pc_x[k] * snf_600[k];

        t_901[k] = f_15 * smf_500[k]
                   + f_3 * pc_y[k] * snf_600[k];

        t_902[k] = f_15 * smf_490[k]
                   + f_3 * pc_z[k] * snf_600[k];

        t_903[k] = f_4 * snd0_363[k]
                   - f_5 * snd1_363[k]
                   + f_3 * pc_x[k] * snf_603[k];
    }

#pragma omp simd aligned(t_904, t_905, t_906, t_907, t_908, pc_x, pc_y, smf_502, snd0_365, \
                         snd1_365, snf_602, snf_605, snf_606, snf_607, \
                         snf_608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_904[k] = f_15 * smf_502[k]
                   + f_3 * pc_y[k] * snf_602[k];

        t_905[k] = f_4 * snd0_365[k]
                   - f_5 * snd1_365[k]
                   + f_3 * pc_x[k] * snf_605[k];

        t_906[k] = f_3 * pc_x[k] * snf_606[k];

        t_907[k] = f_3 * pc_x[k] * snf_607[k];

        t_908[k] = f_3 * pc_x[k] * snf_608[k];
    }

#pragma omp simd aligned(t_909, t_910, t_911, pc_x, pc_y, pc_z, smf_496, smf_506, snd0_363, \
                         snd1_363, snf_606, snf_609 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_909[k] = f_3 * pc_x[k] * snf_609[k];

        t_910[k] = f_15 * smf_506[k]
                   + f_1 * snd0_363[k]
                   - f_2 * snd1_363[k]
                   + f_3 * pc_y[k] * snf_606[k];

        t_911[k] = f_15 * smf_496[k]
                   + f_3 * pc_z[k] * snf_606[k];
    }

#pragma omp simd aligned(t_912, t_913, t_914, pc_y, pc_z, smf_499, smf_508, smf_509, snd0_365, \
                         snd1_365, snf_608, snf_609 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_912[k] = f_15 * smf_508[k]
                   + f_4 * snd0_365[k]
                   - f_5 * snd1_365[k]
                   + f_3 * pc_y[k] * snf_608[k];

        t_913[k] = f_15 * smf_509[k]
                   + f_3 * pc_y[k] * snf_609[k];

        t_914[k] = f_15 * smf_499[k]
                   + f_1 * snd0_365[k]
                   - f_2 * snd1_365[k]
                   + f_3 * pc_z[k] * snf_609[k];
    }

#pragma omp simd aligned(t_915, t_916, t_917, t_918, pc_x, pc_y, pc_z, smf_500, smf_510, \
                         snd0_366, snd0_369, snd1_366, snd1_369, snf_610, \
                         snf_613 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_915[k] = f_1 * snd0_366[k]
                   - f_2 * snd1_366[k]
                   + f_3 * pc_x[k] * snf_610[k];

        t_916[k] = f_14 * smf_510[k]
                   + f_3 * pc_y[k] * snf_610[k];

        t_917[k] = f_13 * smf_500[k]
                   + f_3 * pc_z[k] * snf_610[k];

        t_918[k] = f_4 * snd0_369[k]
                   - f_5 * snd1_369[k]
                   + f_3 * pc_x[k] * snf_613[k];
    }

#pragma omp simd aligned(t_919, t_920, t_921, t_922, t_923, pc_x, pc_y, smf_512, snd0_371, \
                         snd1_371, snf_612, snf_615, snf_616, snf_617, \
                         snf_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_919[k] = f_14 * smf_512[k]
                   + f_3 * pc_y[k] * snf_612[k];

        t_920[k] = f_4 * snd0_371[k]
                   - f_5 * snd1_371[k]
                   + f_3 * pc_x[k] * snf_615[k];

        t_921[k] = f_3 * pc_x[k] * snf_616[k];

        t_922[k] = f_3 * pc_x[k] * snf_617[k];

        t_923[k] = f_3 * pc_x[k] * snf_618[k];
    }

#pragma omp simd aligned(t_924, t_925, t_926, pc_x, pc_y, pc_z, smf_506, smf_516, snd0_369, \
                         snd1_369, snf_616, snf_619 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_924[k] = f_3 * pc_x[k] * snf_619[k];

        t_925[k] = f_14 * smf_516[k]
                   + f_1 * snd0_369[k]
                   - f_2 * snd1_369[k]
                   + f_3 * pc_y[k] * snf_616[k];

        t_926[k] = f_13 * smf_506[k]
                   + f_3 * pc_z[k] * snf_616[k];
    }

#pragma omp simd aligned(t_927, t_928, t_929, pc_y, pc_z, smf_509, smf_518, smf_519, snd0_371, \
                         snd1_371, snf_618, snf_619 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_927[k] = f_14 * smf_518[k]
                   + f_4 * snd0_371[k]
                   - f_5 * snd1_371[k]
                   + f_3 * pc_y[k] * snf_618[k];

        t_928[k] = f_14 * smf_519[k]
                   + f_3 * pc_y[k] * snf_619[k];

        t_929[k] = f_13 * smf_509[k]
                   + f_1 * snd0_371[k]
                   - f_2 * snd1_371[k]
                   + f_3 * pc_z[k] * snf_619[k];
    }

#pragma omp simd aligned(t_930, t_931, t_932, t_933, pc_x, pc_y, pc_z, smf_510, smf_520, \
                         snd0_372, snd0_375, snd1_372, snd1_375, snf_620, \
                         snf_623 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_930[k] = f_1 * snd0_372[k]
                   - f_2 * snd1_372[k]
                   + f_3 * pc_x[k] * snf_620[k];

        t_931[k] = f_12 * smf_520[k]
                   + f_3 * pc_y[k] * snf_620[k];

        t_932[k] = f_11 * smf_510[k]
                   + f_3 * pc_z[k] * snf_620[k];

        t_933[k] = f_4 * snd0_375[k]
                   - f_5 * snd1_375[k]
                   + f_3 * pc_x[k] * snf_623[k];
    }

#pragma omp simd aligned(t_934, t_935, t_936, t_937, t_938, pc_x, pc_y, smf_522, snd0_377, \
                         snd1_377, snf_622, snf_625, snf_626, snf_627, \
                         snf_628 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_934[k] = f_12 * smf_522[k]
                   + f_3 * pc_y[k] * snf_622[k];

        t_935[k] = f_4 * snd0_377[k]
                   - f_5 * snd1_377[k]
                   + f_3 * pc_x[k] * snf_625[k];

        t_936[k] = f_3 * pc_x[k] * snf_626[k];

        t_937[k] = f_3 * pc_x[k] * snf_627[k];

        t_938[k] = f_3 * pc_x[k] * snf_628[k];
    }

#pragma omp simd aligned(t_939, t_940, t_941, pc_x, pc_y, pc_z, smf_516, smf_526, snd0_375, \
                         snd1_375, snf_626, snf_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_939[k] = f_3 * pc_x[k] * snf_629[k];

        t_940[k] = f_12 * smf_526[k]
                   + f_1 * snd0_375[k]
                   - f_2 * snd1_375[k]
                   + f_3 * pc_y[k] * snf_626[k];

        t_941[k] = f_11 * smf_516[k]
                   + f_3 * pc_z[k] * snf_626[k];
    }

#pragma omp simd aligned(t_942, t_943, t_944, pc_y, pc_z, smf_519, smf_528, smf_529, snd0_377, \
                         snd1_377, snf_628, snf_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_942[k] = f_12 * smf_528[k]
                   + f_4 * snd0_377[k]
                   - f_5 * snd1_377[k]
                   + f_3 * pc_y[k] * snf_628[k];

        t_943[k] = f_12 * smf_529[k]
                   + f_3 * pc_y[k] * snf_629[k];

        t_944[k] = f_11 * smf_519[k]
                   + f_1 * snd0_377[k]
                   - f_2 * snd1_377[k]
                   + f_3 * pc_z[k] * snf_629[k];
    }

#pragma omp simd aligned(t_945, t_946, t_947, t_948, pc_x, pc_y, pc_z, smf_520, smf_530, \
                         snd0_378, snd0_381, snd1_378, snd1_381, snf_630, \
                         snf_633 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_945[k] = f_1 * snd0_378[k]
                   - f_2 * snd1_378[k]
                   + f_3 * pc_x[k] * snf_630[k];

        t_946[k] = f_8 * smf_530[k]
                   + f_3 * pc_y[k] * snf_630[k];

        t_947[k] = f_10 * smf_520[k]
                   + f_3 * pc_z[k] * snf_630[k];

        t_948[k] = f_4 * snd0_381[k]
                   - f_5 * snd1_381[k]
                   + f_3 * pc_x[k] * snf_633[k];
    }

#pragma omp simd aligned(t_949, t_950, t_951, t_952, t_953, pc_x, pc_y, smf_532, snd0_383, \
                         snd1_383, snf_632, snf_635, snf_636, snf_637, \
                         snf_638 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_949[k] = f_8 * smf_532[k]
                   + f_3 * pc_y[k] * snf_632[k];

        t_950[k] = f_4 * snd0_383[k]
                   - f_5 * snd1_383[k]
                   + f_3 * pc_x[k] * snf_635[k];

        t_951[k] = f_3 * pc_x[k] * snf_636[k];

        t_952[k] = f_3 * pc_x[k] * snf_637[k];

        t_953[k] = f_3 * pc_x[k] * snf_638[k];
    }

#pragma omp simd aligned(t_954, t_955, t_956, pc_x, pc_y, pc_z, smf_526, smf_536, snd0_381, \
                         snd1_381, snf_636, snf_639 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_954[k] = f_3 * pc_x[k] * snf_639[k];

        t_955[k] = f_8 * smf_536[k]
                   + f_1 * snd0_381[k]
                   - f_2 * snd1_381[k]
                   + f_3 * pc_y[k] * snf_636[k];

        t_956[k] = f_10 * smf_526[k]
                   + f_3 * pc_z[k] * snf_636[k];
    }

#pragma omp simd aligned(t_957, t_958, t_959, t_960, pb_y, pc_y, pc_z, smg0_810, smf_529, \
                         smf_538, smf_539, smg1_810, snd0_383, snd1_383, snf_638, \
                         snf_639 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_957[k] = f_8 * smf_538[k]
                   + f_4 * snd0_383[k]
                   - f_5 * snd1_383[k]
                   + f_3 * pc_y[k] * snf_638[k];

        t_958[k] = f_8 * smf_539[k]
                   + f_3 * pc_y[k] * snf_639[k];

        t_959[k] = f_10 * smf_529[k]
                   + f_1 * snd0_383[k]
                   - f_2 * snd1_383[k]
                   + f_3 * pc_z[k] * snf_639[k];

        t_960[k] = pb_y[k] * smg0_810[k]
                   - f_6 * pc_y[k] * smg1_810[k];
    }

#pragma omp simd aligned(t_961, t_962, t_963, t_964, pc_x, pc_y, pc_z, smf_530, smf_540, \
                         smf_542, snd0_387, snd1_387, snf_640, snf_642, \
                         snf_643 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_961[k] = f_7 * smf_540[k]
                   + f_3 * pc_y[k] * snf_640[k];

        t_962[k] = f_9 * smf_530[k]
                   + f_3 * pc_z[k] * snf_640[k];

        t_963[k] = f_4 * snd0_387[k]
                   - f_5 * snd1_387[k]
                   + f_3 * pc_x[k] * snf_643[k];

        t_964[k] = f_7 * smf_542[k]
                   + f_3 * pc_y[k] * snf_642[k];
    }
}

static auto
compute_prim_sng_three_center_electron_repulsion_0_piece8(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t smg0,
                                                          const size_t smf, const size_t smg1,
                                                          const size_t snd0, const size_t snd1,
                                                          const size_t snf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.0 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 4.5 / q;
    const auto f_14 = 2.0 / q;

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

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *smg0_815 = buffer.data(smg0 + 815);
    const auto *smg0_820 = buffer.data(smg0 + 820);
    const auto *smg0_822 = buffer.data(smg0 + 822);
    const auto *smg0_824 = buffer.data(smg0 + 824);

    const auto *smf_536 = buffer.data(smf + 536);
    const auto *smf_540 = buffer.data(smf + 540);
    const auto *smf_546 = buffer.data(smf + 546);
    const auto *smf_548 = buffer.data(smf + 548);
    const auto *smf_549 = buffer.data(smf + 549);

    const auto *smg1_815 = buffer.data(smg1 + 815);
    const auto *smg1_820 = buffer.data(smg1 + 820);
    const auto *smg1_822 = buffer.data(smg1 + 822);
    const auto *smg1_824 = buffer.data(smg1 + 824);

    const auto *snd0_390 = buffer.data(snd0 + 390);
    const auto *snd0_393 = buffer.data(snd0 + 393);
    const auto *snd0_395 = buffer.data(snd0 + 395);

    const auto *snd1_390 = buffer.data(snd1 + 390);
    const auto *snd1_393 = buffer.data(snd1 + 393);
    const auto *snd1_395 = buffer.data(snd1 + 395);

    const auto *snf_646 = buffer.data(snf + 646);
    const auto *snf_647 = buffer.data(snf + 647);
    const auto *snf_648 = buffer.data(snf + 648);
    const auto *snf_649 = buffer.data(snf + 649);
    const auto *snf_650 = buffer.data(snf + 650);
    const auto *snf_652 = buffer.data(snf + 652);
    const auto *snf_653 = buffer.data(snf + 653);
    const auto *snf_655 = buffer.data(snf + 655);
    const auto *snf_656 = buffer.data(snf + 656);
    const auto *snf_657 = buffer.data(snf + 657);
    const auto *snf_658 = buffer.data(snf + 658);
    const auto *snf_659 = buffer.data(snf + 659);

#pragma omp simd aligned(t_965, t_966, t_967, t_968, t_969, pb_y, pc_x, pc_y, smg0_815, \
                         smg1_815, snf_646, snf_647, snf_648, snf_649 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_965[k] = pb_y[k] * smg0_815[k]
                   - f_6 * pc_y[k] * smg1_815[k];

        t_966[k] = f_3 * pc_x[k] * snf_646[k];

        t_967[k] = f_3 * pc_x[k] * snf_647[k];

        t_968[k] = f_3 * pc_x[k] * snf_648[k];

        t_969[k] = f_3 * pc_x[k] * snf_649[k];
    }

#pragma omp simd aligned(t_970, t_971, t_972, pb_y, pc_y, pc_z, smg0_820, smg0_822, smf_536, \
                         smf_546, smf_548, smg1_820, smg1_822, \
                         snf_646 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_970[k] = pb_y[k] * smg0_820[k]
                   + f_14 * smf_546[k]
                   - f_6 * pc_y[k] * smg1_820[k];

        t_971[k] = f_9 * smf_536[k]
                   + f_3 * pc_z[k] * snf_646[k];

        t_972[k] = pb_y[k] * smg0_822[k]
                   + f_8 * smf_548[k]
                   - f_6 * pc_y[k] * smg1_822[k];
    }

#pragma omp simd aligned(t_973, t_974, t_975, t_976, pb_y, pc_x, pc_y, smg0_824, smf_549, \
                         smg1_824, snd0_390, snd1_390, snf_649, \
                         snf_650 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_973[k] = f_7 * smf_549[k]
                   + f_3 * pc_y[k] * snf_649[k];

        t_974[k] = pb_y[k] * smg0_824[k]
                   - f_6 * pc_y[k] * smg1_824[k];

        t_975[k] = f_1 * snd0_390[k]
                   - f_2 * snd1_390[k]
                   + f_3 * pc_x[k] * snf_650[k];

        t_976[k] = f_3 * pc_y[k] * snf_650[k];
    }

#pragma omp simd aligned(t_977, t_978, t_979, t_980, pc_x, pc_y, pc_z, smf_540, snd0_393, \
                         snd0_395, snd1_393, snd1_395, snf_650, snf_652, snf_653, \
                         snf_655 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_977[k] = f_0 * smf_540[k]
                   + f_3 * pc_z[k] * snf_650[k];

        t_978[k] = f_4 * snd0_393[k]
                   - f_5 * snd1_393[k]
                   + f_3 * pc_x[k] * snf_653[k];

        t_979[k] = f_3 * pc_y[k] * snf_652[k];

        t_980[k] = f_4 * snd0_395[k]
                   - f_5 * snd1_395[k]
                   + f_3 * pc_x[k] * snf_655[k];
    }

#pragma omp simd aligned(t_981, t_982, t_983, t_984, t_985, t_986, pc_x, pc_y, pc_z, smf_546, \
                         snd0_393, snd1_393, snf_656, snf_657, snf_658, \
                         snf_659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_981[k] = f_3 * pc_x[k] * snf_656[k];

        t_982[k] = f_3 * pc_x[k] * snf_657[k];

        t_983[k] = f_3 * pc_x[k] * snf_658[k];

        t_984[k] = f_3 * pc_x[k] * snf_659[k];

        t_985[k] = f_1 * snd0_393[k]
                   - f_2 * snd1_393[k]
                   + f_3 * pc_y[k] * snf_656[k];

        t_986[k] = f_0 * smf_546[k]
                   + f_3 * pc_z[k] * snf_656[k];
    }

#pragma omp simd aligned(t_987, t_988, t_989, pc_y, pc_z, smf_549, snd0_395, snd1_395, \
                         snf_658, snf_659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_987[k] = f_4 * snd0_395[k]
                   - f_5 * snd1_395[k]
                   + f_3 * pc_y[k] * snf_658[k];

        t_988[k] = f_3 * pc_y[k] * snf_659[k];

        t_989[k] = f_0 * smf_549[k]
                   + f_1 * snd0_395[k]
                   - f_2 * snd1_395[k]
                   + f_3 * pc_z[k] * snf_659[k];
    }
}

auto
compute_prim_sng_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t smg0, const size_t smf,
                                                   const size_t smg1, const size_t snd0,
                                                   const size_t snd1, const size_t snf,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_sng_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, smg0, smf,
                                                              smg1, snd0, snd1, snf, ncols,
                                                              gamma, p, q);

    compute_prim_sng_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, smg0, smf,
                                                              smg1, snd0, snd1, snf, ncols,
                                                              gamma, p, q);

    compute_prim_sng_three_center_electron_repulsion_0_piece2(buffer, target, pb, pc, smg0, smf,
                                                              smg1, snd0, snd1, snf, ncols,
                                                              gamma, p, q);

    compute_prim_sng_three_center_electron_repulsion_0_piece3(buffer, target, pb, pc, smg0, smf,
                                                              smg1, snd0, snd1, snf, ncols,
                                                              gamma, p, q);

    compute_prim_sng_three_center_electron_repulsion_0_piece4(buffer, target, pb, pc, smg0, smf,
                                                              smg1, snd0, snd1, snf, ncols,
                                                              gamma, p, q);

    compute_prim_sng_three_center_electron_repulsion_0_piece5(buffer, target, pb, pc, smg0, smf,
                                                              smg1, snd0, snd1, snf, ncols,
                                                              gamma, p, q);

    compute_prim_sng_three_center_electron_repulsion_0_piece6(buffer, target, pb, pc, smg0, smf,
                                                              smg1, snd0, snd1, snf, ncols,
                                                              gamma, p, q);

    compute_prim_sng_three_center_electron_repulsion_0_piece7(buffer, target, pb, pc, smg0, smf,
                                                              smg1, snd0, snd1, snf, ncols,
                                                              gamma, p, q);

    compute_prim_sng_three_center_electron_repulsion_0_piece8(buffer, target, pb, pc, smg0, smf,
                                                              smg1, snd0, snd1, snf, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
