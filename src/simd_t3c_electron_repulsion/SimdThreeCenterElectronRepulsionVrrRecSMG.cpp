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


#include "SimdThreeCenterElectronRepulsionVrrRecSMG.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_smg_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t slg0,
                                                          const size_t slf, const size_t slg1,
                                                          const size_t smd0, const size_t smd1,
                                                          const size_t smf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 4.0 / q;
    const auto f_10 = 3.5 / q;
    const auto f_11 = 3.0 / q;
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

    const auto *slg0_0 = buffer.data(slg0 + 0);
    const auto *slg0_3 = buffer.data(slg0 + 3);
    const auto *slg0_5 = buffer.data(slg0 + 5);
    const auto *slg0_10 = buffer.data(slg0 + 10);
    const auto *slg0_14 = buffer.data(slg0 + 14);
    const auto *slg0_18 = buffer.data(slg0 + 18);
    const auto *slg0_25 = buffer.data(slg0 + 25);
    const auto *slg0_30 = buffer.data(slg0 + 30);
    const auto *slg0_35 = buffer.data(slg0 + 35);
    const auto *slg0_44 = buffer.data(slg0 + 44);
    const auto *slg0_45 = buffer.data(slg0 + 45);
    const auto *slg0_48 = buffer.data(slg0 + 48);
    const auto *slg0_55 = buffer.data(slg0 + 55);
    const auto *slg0_75 = buffer.data(slg0 + 75);
    const auto *slg0_78 = buffer.data(slg0 + 78);

    const auto *slf_0 = buffer.data(slf + 0);
    const auto *slf_1 = buffer.data(slf + 1);
    const auto *slf_2 = buffer.data(slf + 2);
    const auto *slf_3 = buffer.data(slf + 3);
    const auto *slf_5 = buffer.data(slf + 5);
    const auto *slf_6 = buffer.data(slf + 6);
    const auto *slf_7 = buffer.data(slf + 7);
    const auto *slf_8 = buffer.data(slf + 8);
    const auto *slf_9 = buffer.data(slf + 9);
    const auto *slf_10 = buffer.data(slf + 10);
    const auto *slf_12 = buffer.data(slf + 12);
    const auto *slf_16 = buffer.data(slf + 16);
    const auto *slf_17 = buffer.data(slf + 17);
    const auto *slf_18 = buffer.data(slf + 18);
    const auto *slf_19 = buffer.data(slf + 19);
    const auto *slf_20 = buffer.data(slf + 20);
    const auto *slf_22 = buffer.data(slf + 22);
    const auto *slf_26 = buffer.data(slf + 26);
    const auto *slf_27 = buffer.data(slf + 27);
    const auto *slf_28 = buffer.data(slf + 28);
    const auto *slf_29 = buffer.data(slf + 29);
    const auto *slf_30 = buffer.data(slf + 30);
    const auto *slf_32 = buffer.data(slf + 32);
    const auto *slf_33 = buffer.data(slf + 33);
    const auto *slf_35 = buffer.data(slf + 35);
    const auto *slf_36 = buffer.data(slf + 36);
    const auto *slf_37 = buffer.data(slf + 37);
    const auto *slf_38 = buffer.data(slf + 38);
    const auto *slf_39 = buffer.data(slf + 39);
    const auto *slf_40 = buffer.data(slf + 40);
    const auto *slf_42 = buffer.data(slf + 42);
    const auto *slf_46 = buffer.data(slf + 46);
    const auto *slf_47 = buffer.data(slf + 47);
    const auto *slf_48 = buffer.data(slf + 48);
    const auto *slf_49 = buffer.data(slf + 49);
    const auto *slf_50 = buffer.data(slf + 50);
    const auto *slf_51 = buffer.data(slf + 51);
    const auto *slf_52 = buffer.data(slf + 52);
    const auto *slf_53 = buffer.data(slf + 53);
    const auto *slf_55 = buffer.data(slf + 55);
    const auto *slf_56 = buffer.data(slf + 56);
    const auto *slf_57 = buffer.data(slf + 57);
    const auto *slf_58 = buffer.data(slf + 58);
    const auto *slf_59 = buffer.data(slf + 59);
    const auto *slf_60 = buffer.data(slf + 60);
    const auto *slf_63 = buffer.data(slf + 63);
    const auto *slf_65 = buffer.data(slf + 65);
    const auto *slf_66 = buffer.data(slf + 66);
    const auto *slf_67 = buffer.data(slf + 67);
    const auto *slf_68 = buffer.data(slf + 68);
    const auto *slf_69 = buffer.data(slf + 69);
    const auto *slf_75 = buffer.data(slf + 75);
    const auto *slf_76 = buffer.data(slf + 76);
    const auto *slf_77 = buffer.data(slf + 77);
    const auto *slf_78 = buffer.data(slf + 78);
    const auto *slf_79 = buffer.data(slf + 79);

    const auto *slg1_0 = buffer.data(slg1 + 0);
    const auto *slg1_3 = buffer.data(slg1 + 3);
    const auto *slg1_5 = buffer.data(slg1 + 5);
    const auto *slg1_10 = buffer.data(slg1 + 10);
    const auto *slg1_14 = buffer.data(slg1 + 14);
    const auto *slg1_18 = buffer.data(slg1 + 18);
    const auto *slg1_25 = buffer.data(slg1 + 25);
    const auto *slg1_30 = buffer.data(slg1 + 30);
    const auto *slg1_35 = buffer.data(slg1 + 35);
    const auto *slg1_44 = buffer.data(slg1 + 44);
    const auto *slg1_45 = buffer.data(slg1 + 45);
    const auto *slg1_48 = buffer.data(slg1 + 48);
    const auto *slg1_55 = buffer.data(slg1 + 55);
    const auto *slg1_75 = buffer.data(slg1 + 75);
    const auto *slg1_78 = buffer.data(slg1 + 78);

    const auto *smd0_0 = buffer.data(smd0 + 0);
    const auto *smd0_3 = buffer.data(smd0 + 3);
    const auto *smd0_5 = buffer.data(smd0 + 5);
    const auto *smd0_9 = buffer.data(smd0 + 9);
    const auto *smd0_11 = buffer.data(smd0 + 11);
    const auto *smd0_17 = buffer.data(smd0 + 17);
    const auto *smd0_18 = buffer.data(smd0 + 18);
    const auto *smd0_21 = buffer.data(smd0 + 21);
    const auto *smd0_23 = buffer.data(smd0 + 23);
    const auto *smd0_29 = buffer.data(smd0 + 29);
    const auto *smd0_30 = buffer.data(smd0 + 30);
    const auto *smd0_33 = buffer.data(smd0 + 33);
    const auto *smd0_35 = buffer.data(smd0 + 35);
    const auto *smd0_36 = buffer.data(smd0 + 36);
    const auto *smd0_39 = buffer.data(smd0 + 39);
    const auto *smd0_41 = buffer.data(smd0 + 41);
    const auto *smd0_47 = buffer.data(smd0 + 47);

    const auto *smd1_0 = buffer.data(smd1 + 0);
    const auto *smd1_3 = buffer.data(smd1 + 3);
    const auto *smd1_5 = buffer.data(smd1 + 5);
    const auto *smd1_9 = buffer.data(smd1 + 9);
    const auto *smd1_11 = buffer.data(smd1 + 11);
    const auto *smd1_17 = buffer.data(smd1 + 17);
    const auto *smd1_18 = buffer.data(smd1 + 18);
    const auto *smd1_21 = buffer.data(smd1 + 21);
    const auto *smd1_23 = buffer.data(smd1 + 23);
    const auto *smd1_29 = buffer.data(smd1 + 29);
    const auto *smd1_30 = buffer.data(smd1 + 30);
    const auto *smd1_33 = buffer.data(smd1 + 33);
    const auto *smd1_35 = buffer.data(smd1 + 35);
    const auto *smd1_36 = buffer.data(smd1 + 36);
    const auto *smd1_39 = buffer.data(smd1 + 39);
    const auto *smd1_41 = buffer.data(smd1 + 41);
    const auto *smd1_47 = buffer.data(smd1 + 47);

    const auto *smf_0 = buffer.data(smf + 0);
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
    const auto *smf_52 = buffer.data(smf + 52);
    const auto *smf_53 = buffer.data(smf + 53);
    const auto *smf_55 = buffer.data(smf + 55);
    const auto *smf_56 = buffer.data(smf + 56);
    const auto *smf_57 = buffer.data(smf + 57);
    const auto *smf_58 = buffer.data(smf + 58);
    const auto *smf_59 = buffer.data(smf + 59);
    const auto *smf_60 = buffer.data(smf + 60);
    const auto *smf_62 = buffer.data(smf + 62);
    const auto *smf_63 = buffer.data(smf + 63);
    const auto *smf_65 = buffer.data(smf + 65);
    const auto *smf_66 = buffer.data(smf + 66);
    const auto *smf_67 = buffer.data(smf + 67);
    const auto *smf_68 = buffer.data(smf + 68);
    const auto *smf_69 = buffer.data(smf + 69);
    const auto *smf_70 = buffer.data(smf + 70);
    const auto *smf_72 = buffer.data(smf + 72);
    const auto *smf_75 = buffer.data(smf + 75);
    const auto *smf_76 = buffer.data(smf + 76);
    const auto *smf_77 = buffer.data(smf + 77);
    const auto *smf_78 = buffer.data(smf + 78);
    const auto *smf_79 = buffer.data(smf + 79);
    const auto *smf_80 = buffer.data(smf + 80);
    const auto *smf_82 = buffer.data(smf + 82);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, slf_0, slf_3, smd0_0, smd0_3, \
                         smd1_0, smd1_3, smf_0, smf_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * slf_0[k]
                 + f_1 * smd0_0[k]
                 - f_2 * smd1_0[k]
                 + f_3 * pc_x[k] * smf_0[k];

        t_1[k] = f_3 * pc_y[k] * smf_0[k];

        t_2[k] = f_3 * pc_z[k] * smf_0[k];

        t_3[k] = f_0 * slf_3[k]
                 + f_4 * smd0_3[k]
                 - f_5 * smd1_3[k]
                 + f_3 * pc_x[k] * smf_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pc_x, pc_y, slf_5, slf_6, slf_7, smd0_5, smd1_5, \
                         smf_2, smf_5, smf_6, smf_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * pc_y[k] * smf_2[k];

        t_5[k] = f_0 * slf_5[k]
                 + f_4 * smd0_5[k]
                 - f_5 * smd1_5[k]
                 + f_3 * pc_x[k] * smf_5[k];

        t_6[k] = f_0 * slf_6[k]
                 + f_3 * pc_x[k] * smf_6[k];

        t_7[k] = f_0 * slf_7[k]
                 + f_3 * pc_x[k] * smf_7[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pc_x, pc_y, pc_z, slf_8, slf_9, smd0_3, smd1_3, \
                         smf_6, smf_8, smf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * slf_8[k]
                 + f_3 * pc_x[k] * smf_8[k];

        t_9[k] = f_0 * slf_9[k]
                 + f_3 * pc_x[k] * smf_9[k];

        t_10[k] = f_1 * smd0_3[k]
                  - f_2 * smd1_3[k]
                  + f_3 * pc_y[k] * smf_6[k];

        t_11[k] = f_3 * pc_z[k] * smf_6[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, pb_y, pc_y, pc_z, slg0_0, slf_0, \
                         slg1_0, smd0_5, smd1_5, smf_8, smf_9, smf_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_4 * smd0_5[k]
                  - f_5 * smd1_5[k]
                  + f_3 * pc_y[k] * smf_8[k];

        t_13[k] = f_3 * pc_y[k] * smf_9[k];

        t_14[k] = f_1 * smd0_5[k]
                  - f_2 * smd1_5[k]
                  + f_3 * pc_z[k] * smf_9[k];

        t_15[k] = pb_y[k] * slg0_0[k]
                  - f_6 * pc_y[k] * slg1_0[k];

        t_16[k] = f_7 * slf_0[k]
                  + f_3 * pc_y[k] * smf_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pb_y, pc_y, pc_z, slg0_3, slg0_5, slf_1, \
                         slf_2, slg1_3, slg1_5, smf_10, smf_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_3 * pc_z[k] * smf_10[k];

        t_18[k] = pb_y[k] * slg0_3[k]
                  + f_8 * slf_1[k]
                  - f_6 * pc_y[k] * slg1_3[k];

        t_19[k] = f_7 * slf_2[k]
                  + f_3 * pc_y[k] * smf_12[k];

        t_20[k] = pb_y[k] * slg0_5[k]
                  - f_6 * pc_y[k] * slg1_5[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pc_x, slf_16, slf_17, slf_18, slf_19, smf_16, \
                         smf_17, smf_18, smf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_9 * slf_16[k]
                  + f_3 * pc_x[k] * smf_16[k];

        t_22[k] = f_9 * slf_17[k]
                  + f_3 * pc_x[k] * smf_17[k];

        t_23[k] = f_9 * slf_18[k]
                  + f_3 * pc_x[k] * smf_18[k];

        t_24[k] = f_9 * slf_19[k]
                  + f_3 * pc_x[k] * smf_19[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pc_y, pc_z, slf_6, slf_8, slf_9, smd0_9, \
                         smd0_11, smd1_9, smd1_11, smf_16, smf_18, \
                         smf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_7 * slf_6[k]
                  + f_1 * smd0_9[k]
                  - f_2 * smd1_9[k]
                  + f_3 * pc_y[k] * smf_16[k];

        t_26[k] = f_3 * pc_z[k] * smf_16[k];

        t_27[k] = f_7 * slf_8[k]
                  + f_4 * smd0_11[k]
                  - f_5 * smd1_11[k]
                  + f_3 * pc_y[k] * smf_18[k];

        t_28[k] = f_7 * slf_9[k]
                  + f_3 * pc_y[k] * smf_19[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pb_y, pb_z, pc_y, pc_z, slg0_0, slg0_14, \
                         slf_0, slg1_0, slg1_14, smf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pb_y[k] * slg0_14[k]
                  - f_6 * pc_y[k] * slg1_14[k];

        t_30[k] = pb_z[k] * slg0_0[k]
                  - f_6 * pc_z[k] * slg1_0[k];

        t_31[k] = f_3 * pc_y[k] * smf_20[k];

        t_32[k] = f_7 * slf_0[k]
                  + f_3 * pc_z[k] * smf_20[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pb_z, pc_x, pc_y, pc_z, slg0_3, slg0_5, \
                         slf_2, slf_26, slg1_3, slg1_5, smf_22, \
                         smf_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pb_z[k] * slg0_3[k]
                  - f_6 * pc_z[k] * slg1_3[k];

        t_34[k] = f_3 * pc_y[k] * smf_22[k];

        t_35[k] = pb_z[k] * slg0_5[k]
                  + f_8 * slf_2[k]
                  - f_6 * pc_z[k] * slg1_5[k];

        t_36[k] = f_9 * slf_26[k]
                  + f_3 * pc_x[k] * smf_26[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pb_z, pc_x, pc_z, slg0_10, slf_27, slf_28, \
                         slf_29, slg1_10, smf_27, smf_28, smf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_9 * slf_27[k]
                  + f_3 * pc_x[k] * smf_27[k];

        t_38[k] = f_9 * slf_28[k]
                  + f_3 * pc_x[k] * smf_28[k];

        t_39[k] = f_9 * slf_29[k]
                  + f_3 * pc_x[k] * smf_29[k];

        t_40[k] = pb_z[k] * slg0_10[k]
                  - f_6 * pc_z[k] * slg1_10[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pc_y, pc_z, slf_6, slf_9, smd0_17, smd1_17, \
                         smf_26, smf_28, smf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_7 * slf_6[k]
                  + f_3 * pc_z[k] * smf_26[k];

        t_42[k] = f_4 * smd0_17[k]
                  - f_5 * smd1_17[k]
                  + f_3 * pc_y[k] * smf_28[k];

        t_43[k] = f_3 * pc_y[k] * smf_29[k];

        t_44[k] = f_7 * slf_9[k]
                  + f_1 * smd0_17[k]
                  - f_2 * smd1_17[k]
                  + f_3 * pc_z[k] * smf_29[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pc_x, pc_y, pc_z, slf_10, slf_30, slf_33, \
                         smd0_18, smd0_21, smd1_18, smd1_21, smf_30, \
                         smf_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_10 * slf_30[k]
                  + f_1 * smd0_18[k]
                  - f_2 * smd1_18[k]
                  + f_3 * pc_x[k] * smf_30[k];

        t_46[k] = f_8 * slf_10[k]
                  + f_3 * pc_y[k] * smf_30[k];

        t_47[k] = f_3 * pc_z[k] * smf_30[k];

        t_48[k] = f_10 * slf_33[k]
                  + f_4 * smd0_21[k]
                  - f_5 * smd1_21[k]
                  + f_3 * pc_x[k] * smf_33[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pc_x, pc_y, slf_12, slf_35, slf_36, slf_37, \
                         smd0_23, smd1_23, smf_32, smf_35, smf_36, \
                         smf_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_8 * slf_12[k]
                  + f_3 * pc_y[k] * smf_32[k];

        t_50[k] = f_10 * slf_35[k]
                  + f_4 * smd0_23[k]
                  - f_5 * smd1_23[k]
                  + f_3 * pc_x[k] * smf_35[k];

        t_51[k] = f_10 * slf_36[k]
                  + f_3 * pc_x[k] * smf_36[k];

        t_52[k] = f_10 * slf_37[k]
                  + f_3 * pc_x[k] * smf_37[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, pc_x, pc_y, pc_z, slf_16, slf_38, slf_39, \
                         smd0_21, smd1_21, smf_36, smf_38, smf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_10 * slf_38[k]
                  + f_3 * pc_x[k] * smf_38[k];

        t_54[k] = f_10 * slf_39[k]
                  + f_3 * pc_x[k] * smf_39[k];

        t_55[k] = f_8 * slf_16[k]
                  + f_1 * smd0_21[k]
                  - f_2 * smd1_21[k]
                  + f_3 * pc_y[k] * smf_36[k];

        t_56[k] = f_3 * pc_z[k] * smf_36[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, pb_y, pc_y, pc_z, slg0_30, slf_18, slf_19, \
                         slg1_30, smd0_23, smd1_23, smf_38, smf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_8 * slf_18[k]
                  + f_4 * smd0_23[k]
                  - f_5 * smd1_23[k]
                  + f_3 * pc_y[k] * smf_38[k];

        t_58[k] = f_8 * slf_19[k]
                  + f_3 * pc_y[k] * smf_39[k];

        t_59[k] = f_1 * smd0_23[k]
                  - f_2 * smd1_23[k]
                  + f_3 * pc_z[k] * smf_39[k];

        t_60[k] = pb_y[k] * slg0_30[k]
                  - f_6 * pc_y[k] * slg1_30[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pb_z, pc_y, pc_z, slg0_18, slf_10, slf_20, \
                         slf_22, slg1_18, smf_40, smf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_7 * slf_20[k]
                  + f_3 * pc_y[k] * smf_40[k];

        t_62[k] = f_7 * slf_10[k]
                  + f_3 * pc_z[k] * smf_40[k];

        t_63[k] = pb_z[k] * slg0_18[k]
                  - f_6 * pc_z[k] * slg1_18[k];

        t_64[k] = f_7 * slf_22[k]
                  + f_3 * pc_y[k] * smf_42[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pb_y, pc_x, pc_y, slg0_35, slf_46, slf_47, \
                         slf_48, slg1_35, smf_46, smf_47, smf_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = pb_y[k] * slg0_35[k]
                  - f_6 * pc_y[k] * slg1_35[k];

        t_66[k] = f_10 * slf_46[k]
                  + f_3 * pc_x[k] * smf_46[k];

        t_67[k] = f_10 * slf_47[k]
                  + f_3 * pc_x[k] * smf_47[k];

        t_68[k] = f_10 * slf_48[k]
                  + f_3 * pc_x[k] * smf_48[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pb_z, pc_x, pc_z, slg0_25, slf_16, slf_49, slg1_25, \
                         smf_46, smf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_10 * slf_49[k]
                  + f_3 * pc_x[k] * smf_49[k];

        t_70[k] = pb_z[k] * slg0_25[k]
                  - f_6 * pc_z[k] * slg1_25[k];

        t_71[k] = f_7 * slf_16[k]
                  + f_3 * pc_z[k] * smf_46[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, pb_y, pc_y, slg0_44, slf_28, slf_29, slg1_44, \
                         smd0_29, smd1_29, smf_48, smf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_7 * slf_28[k]
                  + f_4 * smd0_29[k]
                  - f_5 * smd1_29[k]
                  + f_3 * pc_y[k] * smf_48[k];

        t_73[k] = f_7 * slf_29[k]
                  + f_3 * pc_y[k] * smf_49[k];

        t_74[k] = pb_y[k] * slg0_44[k]
                  - f_6 * pc_y[k] * slg1_44[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pc_x, pc_y, pc_z, slf_20, slf_50, slf_53, \
                         smd0_30, smd0_33, smd1_30, smd1_33, smf_50, \
                         smf_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_10 * slf_50[k]
                  + f_1 * smd0_30[k]
                  - f_2 * smd1_30[k]
                  + f_3 * pc_x[k] * smf_50[k];

        t_76[k] = f_3 * pc_y[k] * smf_50[k];

        t_77[k] = f_8 * slf_20[k]
                  + f_3 * pc_z[k] * smf_50[k];

        t_78[k] = f_10 * slf_53[k]
                  + f_4 * smd0_33[k]
                  - f_5 * smd1_33[k]
                  + f_3 * pc_x[k] * smf_53[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pc_x, pc_y, slf_55, slf_56, slf_57, smd0_35, \
                         smd1_35, smf_52, smf_55, smf_56, smf_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_3 * pc_y[k] * smf_52[k];

        t_80[k] = f_10 * slf_55[k]
                  + f_4 * smd0_35[k]
                  - f_5 * smd1_35[k]
                  + f_3 * pc_x[k] * smf_55[k];

        t_81[k] = f_10 * slf_56[k]
                  + f_3 * pc_x[k] * smf_56[k];

        t_82[k] = f_10 * slf_57[k]
                  + f_3 * pc_x[k] * smf_57[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pc_x, pc_y, pc_z, slf_26, slf_58, slf_59, \
                         smd0_33, smd1_33, smf_56, smf_58, smf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_10 * slf_58[k]
                  + f_3 * pc_x[k] * smf_58[k];

        t_84[k] = f_10 * slf_59[k]
                  + f_3 * pc_x[k] * smf_59[k];

        t_85[k] = f_1 * smd0_33[k]
                  - f_2 * smd1_33[k]
                  + f_3 * pc_y[k] * smf_56[k];

        t_86[k] = f_8 * slf_26[k]
                  + f_3 * pc_z[k] * smf_56[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pc_x, pc_y, pc_z, slf_29, slf_60, smd0_35, \
                         smd0_36, smd1_35, smd1_36, smf_58, smf_59, \
                         smf_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_4 * smd0_35[k]
                  - f_5 * smd1_35[k]
                  + f_3 * pc_y[k] * smf_58[k];

        t_88[k] = f_3 * pc_y[k] * smf_59[k];

        t_89[k] = f_8 * slf_29[k]
                  + f_1 * smd0_35[k]
                  - f_2 * smd1_35[k]
                  + f_3 * pc_z[k] * smf_59[k];

        t_90[k] = f_11 * slf_60[k]
                  + f_1 * smd0_36[k]
                  - f_2 * smd1_36[k]
                  + f_3 * pc_x[k] * smf_60[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, pc_x, pc_y, pc_z, slf_30, slf_32, slf_63, \
                         smd0_39, smd1_39, smf_60, smf_62, smf_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_12 * slf_30[k]
                  + f_3 * pc_y[k] * smf_60[k];

        t_92[k] = f_3 * pc_z[k] * smf_60[k];

        t_93[k] = f_11 * slf_63[k]
                  + f_4 * smd0_39[k]
                  - f_5 * smd1_39[k]
                  + f_3 * pc_x[k] * smf_63[k];

        t_94[k] = f_12 * slf_32[k]
                  + f_3 * pc_y[k] * smf_62[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pc_x, slf_65, slf_66, slf_67, slf_68, \
                         smd0_41, smd1_41, smf_65, smf_66, smf_67, \
                         smf_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_11 * slf_65[k]
                  + f_4 * smd0_41[k]
                  - f_5 * smd1_41[k]
                  + f_3 * pc_x[k] * smf_65[k];

        t_96[k] = f_11 * slf_66[k]
                  + f_3 * pc_x[k] * smf_66[k];

        t_97[k] = f_11 * slf_67[k]
                  + f_3 * pc_x[k] * smf_67[k];

        t_98[k] = f_11 * slf_68[k]
                  + f_3 * pc_x[k] * smf_68[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, pc_x, pc_y, pc_z, slf_36, slf_69, smd0_39, \
                         smd1_39, smf_66, smf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_11 * slf_69[k]
                  + f_3 * pc_x[k] * smf_69[k];

        t_100[k] = f_12 * slf_36[k]
                   + f_1 * smd0_39[k]
                   - f_2 * smd1_39[k]
                   + f_3 * pc_y[k] * smf_66[k];

        t_101[k] = f_3 * pc_z[k] * smf_66[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, pb_z, pc_y, pc_z, slg0_45, slf_38, \
                         slf_39, slg1_45, smd0_41, smd1_41, smf_68, \
                         smf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_12 * slf_38[k]
                   + f_4 * smd0_41[k]
                   - f_5 * smd1_41[k]
                   + f_3 * pc_y[k] * smf_68[k];

        t_103[k] = f_12 * slf_39[k]
                   + f_3 * pc_y[k] * smf_69[k];

        t_104[k] = f_1 * smd0_41[k]
                   - f_2 * smd1_41[k]
                   + f_3 * pc_z[k] * smf_69[k];

        t_105[k] = pb_z[k] * slg0_45[k]
                   - f_6 * pc_z[k] * slg1_45[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, pb_z, pc_y, pc_z, slg0_48, slf_30, \
                         slf_40, slf_42, slg1_48, smf_70, smf_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_8 * slf_40[k]
                   + f_3 * pc_y[k] * smf_70[k];

        t_107[k] = f_7 * slf_30[k]
                   + f_3 * pc_z[k] * smf_70[k];

        t_108[k] = pb_z[k] * slg0_48[k]
                   - f_6 * pc_z[k] * slg1_48[k];

        t_109[k] = f_8 * slf_42[k]
                   + f_3 * pc_y[k] * smf_72[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, pc_x, slf_75, slf_76, slf_77, slf_78, \
                         smd0_47, smd1_47, smf_75, smf_76, smf_77, \
                         smf_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = f_11 * slf_75[k]
                   + f_4 * smd0_47[k]
                   - f_5 * smd1_47[k]
                   + f_3 * pc_x[k] * smf_75[k];

        t_111[k] = f_11 * slf_76[k]
                   + f_3 * pc_x[k] * smf_76[k];

        t_112[k] = f_11 * slf_77[k]
                   + f_3 * pc_x[k] * smf_77[k];

        t_113[k] = f_11 * slf_78[k]
                   + f_3 * pc_x[k] * smf_78[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, pb_z, pc_x, pc_z, slg0_55, slf_36, slf_79, \
                         slg1_55, smf_76, smf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_11 * slf_79[k]
                   + f_3 * pc_x[k] * smf_79[k];

        t_115[k] = pb_z[k] * slg0_55[k]
                   - f_6 * pc_z[k] * slg1_55[k];

        t_116[k] = f_7 * slf_36[k]
                   + f_3 * pc_z[k] * smf_76[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pb_y, pc_y, pc_z, slg0_75, slf_39, \
                         slf_48, slf_49, slg1_75, smd0_47, smd1_47, smf_78, \
                         smf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_8 * slf_48[k]
                   + f_4 * smd0_47[k]
                   - f_5 * smd1_47[k]
                   + f_3 * pc_y[k] * smf_78[k];

        t_118[k] = f_8 * slf_49[k]
                   + f_3 * pc_y[k] * smf_79[k];

        t_119[k] = f_7 * slf_39[k]
                   + f_1 * smd0_47[k]
                   - f_2 * smd1_47[k]
                   + f_3 * pc_z[k] * smf_79[k];

        t_120[k] = pb_y[k] * slg0_75[k]
                   - f_6 * pc_y[k] * slg1_75[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pb_y, pc_y, pc_z, slg0_78, slf_40, \
                         slf_50, slf_51, slf_52, slg1_78, smf_80, \
                         smf_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_7 * slf_50[k]
                   + f_3 * pc_y[k] * smf_80[k];

        t_122[k] = f_8 * slf_40[k]
                   + f_3 * pc_z[k] * smf_80[k];

        t_123[k] = pb_y[k] * slg0_78[k]
                   + f_8 * slf_51[k]
                   - f_6 * pc_y[k] * slg1_78[k];

        t_124[k] = f_7 * slf_52[k]
                   + f_3 * pc_y[k] * smf_82[k];
    }
}

static auto
compute_prim_smg_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t slg0,
                                                          const size_t slf, const size_t slg1,
                                                          const size_t smd0, const size_t smd1,
                                                          const size_t smf, const size_t ncols,
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
    const auto f_11 = 3.0 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 2.5 / q;
    const auto f_14 = 2.0 / q;

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

    const auto *slg0_80 = buffer.data(slg0 + 80);
    const auto *slg0_89 = buffer.data(slg0 + 89);
    const auto *slg0_90 = buffer.data(slg0 + 90);
    const auto *slg0_93 = buffer.data(slg0 + 93);
    const auto *slg0_100 = buffer.data(slg0 + 100);
    const auto *slg0_135 = buffer.data(slg0 + 135);
    const auto *slg0_138 = buffer.data(slg0 + 138);
    const auto *slg0_140 = buffer.data(slg0 + 140);
    const auto *slg0_149 = buffer.data(slg0 + 149);
    const auto *slg0_150 = buffer.data(slg0 + 150);
    const auto *slg0_153 = buffer.data(slg0 + 153);

    const auto *slf_46 = buffer.data(slf + 46);
    const auto *slf_50 = buffer.data(slf + 50);
    const auto *slf_56 = buffer.data(slf + 56);
    const auto *slf_58 = buffer.data(slf + 58);
    const auto *slf_59 = buffer.data(slf + 59);
    const auto *slf_60 = buffer.data(slf + 60);
    const auto *slf_62 = buffer.data(slf + 62);
    const auto *slf_66 = buffer.data(slf + 66);
    const auto *slf_68 = buffer.data(slf + 68);
    const auto *slf_69 = buffer.data(slf + 69);
    const auto *slf_70 = buffer.data(slf + 70);
    const auto *slf_72 = buffer.data(slf + 72);
    const auto *slf_76 = buffer.data(slf + 76);
    const auto *slf_78 = buffer.data(slf + 78);
    const auto *slf_79 = buffer.data(slf + 79);
    const auto *slf_80 = buffer.data(slf + 80);
    const auto *slf_82 = buffer.data(slf + 82);
    const auto *slf_86 = buffer.data(slf + 86);
    const auto *slf_87 = buffer.data(slf + 87);
    const auto *slf_88 = buffer.data(slf + 88);
    const auto *slf_89 = buffer.data(slf + 89);
    const auto *slf_90 = buffer.data(slf + 90);
    const auto *slf_91 = buffer.data(slf + 91);
    const auto *slf_92 = buffer.data(slf + 92);
    const auto *slf_93 = buffer.data(slf + 93);
    const auto *slf_95 = buffer.data(slf + 95);
    const auto *slf_96 = buffer.data(slf + 96);
    const auto *slf_97 = buffer.data(slf + 97);
    const auto *slf_98 = buffer.data(slf + 98);
    const auto *slf_99 = buffer.data(slf + 99);
    const auto *slf_100 = buffer.data(slf + 100);
    const auto *slf_102 = buffer.data(slf + 102);
    const auto *slf_103 = buffer.data(slf + 103);
    const auto *slf_105 = buffer.data(slf + 105);
    const auto *slf_106 = buffer.data(slf + 106);
    const auto *slf_107 = buffer.data(slf + 107);
    const auto *slf_108 = buffer.data(slf + 108);
    const auto *slf_109 = buffer.data(slf + 109);
    const auto *slf_110 = buffer.data(slf + 110);
    const auto *slf_112 = buffer.data(slf + 112);
    const auto *slf_115 = buffer.data(slf + 115);
    const auto *slf_116 = buffer.data(slf + 116);
    const auto *slf_117 = buffer.data(slf + 117);
    const auto *slf_118 = buffer.data(slf + 118);
    const auto *slf_119 = buffer.data(slf + 119);
    const auto *slf_120 = buffer.data(slf + 120);
    const auto *slf_123 = buffer.data(slf + 123);
    const auto *slf_125 = buffer.data(slf + 125);
    const auto *slf_126 = buffer.data(slf + 126);
    const auto *slf_127 = buffer.data(slf + 127);
    const auto *slf_128 = buffer.data(slf + 128);
    const auto *slf_129 = buffer.data(slf + 129);
    const auto *slf_136 = buffer.data(slf + 136);
    const auto *slf_137 = buffer.data(slf + 137);
    const auto *slf_138 = buffer.data(slf + 138);
    const auto *slf_139 = buffer.data(slf + 139);
    const auto *slf_140 = buffer.data(slf + 140);
    const auto *slf_143 = buffer.data(slf + 143);
    const auto *slf_145 = buffer.data(slf + 145);
    const auto *slf_146 = buffer.data(slf + 146);
    const auto *slf_147 = buffer.data(slf + 147);
    const auto *slf_148 = buffer.data(slf + 148);
    const auto *slf_149 = buffer.data(slf + 149);
    const auto *slf_150 = buffer.data(slf + 150);
    const auto *slf_153 = buffer.data(slf + 153);
    const auto *slf_155 = buffer.data(slf + 155);
    const auto *slf_156 = buffer.data(slf + 156);
    const auto *slf_157 = buffer.data(slf + 157);
    const auto *slf_158 = buffer.data(slf + 158);
    const auto *slf_159 = buffer.data(slf + 159);

    const auto *slg1_80 = buffer.data(slg1 + 80);
    const auto *slg1_89 = buffer.data(slg1 + 89);
    const auto *slg1_90 = buffer.data(slg1 + 90);
    const auto *slg1_93 = buffer.data(slg1 + 93);
    const auto *slg1_100 = buffer.data(slg1 + 100);
    const auto *slg1_135 = buffer.data(slg1 + 135);
    const auto *slg1_138 = buffer.data(slg1 + 138);
    const auto *slg1_140 = buffer.data(slg1 + 140);
    const auto *slg1_149 = buffer.data(slg1 + 149);
    const auto *slg1_150 = buffer.data(slg1 + 150);
    const auto *slg1_153 = buffer.data(slg1 + 153);

    const auto *smd0_51 = buffer.data(smd0 + 51);
    const auto *smd0_53 = buffer.data(smd0 + 53);
    const auto *smd0_54 = buffer.data(smd0 + 54);
    const auto *smd0_57 = buffer.data(smd0 + 57);
    const auto *smd0_59 = buffer.data(smd0 + 59);
    const auto *smd0_60 = buffer.data(smd0 + 60);
    const auto *smd0_63 = buffer.data(smd0 + 63);
    const auto *smd0_65 = buffer.data(smd0 + 65);
    const auto *smd0_71 = buffer.data(smd0 + 71);
    const auto *smd0_72 = buffer.data(smd0 + 72);
    const auto *smd0_75 = buffer.data(smd0 + 75);
    const auto *smd0_77 = buffer.data(smd0 + 77);
    const auto *smd0_81 = buffer.data(smd0 + 81);
    const auto *smd0_83 = buffer.data(smd0 + 83);
    const auto *smd0_84 = buffer.data(smd0 + 84);
    const auto *smd0_87 = buffer.data(smd0 + 87);
    const auto *smd0_89 = buffer.data(smd0 + 89);
    const auto *smd0_90 = buffer.data(smd0 + 90);
    const auto *smd0_93 = buffer.data(smd0 + 93);
    const auto *smd0_95 = buffer.data(smd0 + 95);

    const auto *smd1_51 = buffer.data(smd1 + 51);
    const auto *smd1_53 = buffer.data(smd1 + 53);
    const auto *smd1_54 = buffer.data(smd1 + 54);
    const auto *smd1_57 = buffer.data(smd1 + 57);
    const auto *smd1_59 = buffer.data(smd1 + 59);
    const auto *smd1_60 = buffer.data(smd1 + 60);
    const auto *smd1_63 = buffer.data(smd1 + 63);
    const auto *smd1_65 = buffer.data(smd1 + 65);
    const auto *smd1_71 = buffer.data(smd1 + 71);
    const auto *smd1_72 = buffer.data(smd1 + 72);
    const auto *smd1_75 = buffer.data(smd1 + 75);
    const auto *smd1_77 = buffer.data(smd1 + 77);
    const auto *smd1_81 = buffer.data(smd1 + 81);
    const auto *smd1_83 = buffer.data(smd1 + 83);
    const auto *smd1_84 = buffer.data(smd1 + 84);
    const auto *smd1_87 = buffer.data(smd1 + 87);
    const auto *smd1_89 = buffer.data(smd1 + 89);
    const auto *smd1_90 = buffer.data(smd1 + 90);
    const auto *smd1_93 = buffer.data(smd1 + 93);
    const auto *smd1_95 = buffer.data(smd1 + 95);

    const auto *smf_86 = buffer.data(smf + 86);
    const auto *smf_87 = buffer.data(smf + 87);
    const auto *smf_88 = buffer.data(smf + 88);
    const auto *smf_89 = buffer.data(smf + 89);
    const auto *smf_90 = buffer.data(smf + 90);
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
    const auto *smf_122 = buffer.data(smf + 122);
    const auto *smf_123 = buffer.data(smf + 123);
    const auto *smf_125 = buffer.data(smf + 125);
    const auto *smf_126 = buffer.data(smf + 126);
    const auto *smf_127 = buffer.data(smf + 127);
    const auto *smf_128 = buffer.data(smf + 128);
    const auto *smf_129 = buffer.data(smf + 129);
    const auto *smf_130 = buffer.data(smf + 130);
    const auto *smf_132 = buffer.data(smf + 132);
    const auto *smf_136 = buffer.data(smf + 136);
    const auto *smf_137 = buffer.data(smf + 137);
    const auto *smf_138 = buffer.data(smf + 138);
    const auto *smf_139 = buffer.data(smf + 139);
    const auto *smf_140 = buffer.data(smf + 140);
    const auto *smf_142 = buffer.data(smf + 142);
    const auto *smf_143 = buffer.data(smf + 143);
    const auto *smf_145 = buffer.data(smf + 145);
    const auto *smf_146 = buffer.data(smf + 146);
    const auto *smf_147 = buffer.data(smf + 147);
    const auto *smf_148 = buffer.data(smf + 148);
    const auto *smf_149 = buffer.data(smf + 149);
    const auto *smf_150 = buffer.data(smf + 150);
    const auto *smf_152 = buffer.data(smf + 152);
    const auto *smf_153 = buffer.data(smf + 153);
    const auto *smf_155 = buffer.data(smf + 155);
    const auto *smf_156 = buffer.data(smf + 156);
    const auto *smf_157 = buffer.data(smf + 157);
    const auto *smf_158 = buffer.data(smf + 158);
    const auto *smf_159 = buffer.data(smf + 159);
    const auto *smf_160 = buffer.data(smf + 160);
    const auto *smf_162 = buffer.data(smf + 162);

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pb_y, pc_x, pc_y, slg0_80, slf_86, \
                         slf_87, slf_88, slg1_80, smf_86, smf_87, \
                         smf_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = pb_y[k] * slg0_80[k]
                   - f_6 * pc_y[k] * slg1_80[k];

        t_126[k] = f_11 * slf_86[k]
                   + f_3 * pc_x[k] * smf_86[k];

        t_127[k] = f_11 * slf_87[k]
                   + f_3 * pc_x[k] * smf_87[k];

        t_128[k] = f_11 * slf_88[k]
                   + f_3 * pc_x[k] * smf_88[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, pc_x, pc_y, pc_z, slf_46, slf_56, slf_89, \
                         smd0_51, smd1_51, smf_86, smf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_11 * slf_89[k]
                   + f_3 * pc_x[k] * smf_89[k];

        t_130[k] = f_7 * slf_56[k]
                   + f_1 * smd0_51[k]
                   - f_2 * smd1_51[k]
                   + f_3 * pc_y[k] * smf_86[k];

        t_131[k] = f_8 * slf_46[k]
                   + f_3 * pc_z[k] * smf_86[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, pb_y, pc_y, slg0_89, slf_58, slf_59, slg1_89, \
                         smd0_53, smd1_53, smf_88, smf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_7 * slf_58[k]
                   + f_4 * smd0_53[k]
                   - f_5 * smd1_53[k]
                   + f_3 * pc_y[k] * smf_88[k];

        t_133[k] = f_7 * slf_59[k]
                   + f_3 * pc_y[k] * smf_89[k];

        t_134[k] = pb_y[k] * slg0_89[k]
                   - f_6 * pc_y[k] * slg1_89[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, pc_x, pc_y, pc_z, slf_50, slf_90, slf_93, \
                         smd0_54, smd0_57, smd1_54, smd1_57, smf_90, \
                         smf_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_11 * slf_90[k]
                   + f_1 * smd0_54[k]
                   - f_2 * smd1_54[k]
                   + f_3 * pc_x[k] * smf_90[k];

        t_136[k] = f_3 * pc_y[k] * smf_90[k];

        t_137[k] = f_12 * slf_50[k]
                   + f_3 * pc_z[k] * smf_90[k];

        t_138[k] = f_11 * slf_93[k]
                   + f_4 * smd0_57[k]
                   - f_5 * smd1_57[k]
                   + f_3 * pc_x[k] * smf_93[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, pc_x, pc_y, slf_95, slf_96, slf_97, \
                         smd0_59, smd1_59, smf_92, smf_95, smf_96, \
                         smf_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_3 * pc_y[k] * smf_92[k];

        t_140[k] = f_11 * slf_95[k]
                   + f_4 * smd0_59[k]
                   - f_5 * smd1_59[k]
                   + f_3 * pc_x[k] * smf_95[k];

        t_141[k] = f_11 * slf_96[k]
                   + f_3 * pc_x[k] * smf_96[k];

        t_142[k] = f_11 * slf_97[k]
                   + f_3 * pc_x[k] * smf_97[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, t_146, pc_x, pc_y, pc_z, slf_56, slf_98, slf_99, \
                         smd0_57, smd1_57, smf_96, smf_98, smf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_11 * slf_98[k]
                   + f_3 * pc_x[k] * smf_98[k];

        t_144[k] = f_11 * slf_99[k]
                   + f_3 * pc_x[k] * smf_99[k];

        t_145[k] = f_1 * smd0_57[k]
                   - f_2 * smd1_57[k]
                   + f_3 * pc_y[k] * smf_96[k];

        t_146[k] = f_12 * slf_56[k]
                   + f_3 * pc_z[k] * smf_96[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, pc_x, pc_y, pc_z, slf_59, slf_100, \
                         smd0_59, smd0_60, smd1_59, smd1_60, smf_98, smf_99, \
                         smf_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = f_4 * smd0_59[k]
                   - f_5 * smd1_59[k]
                   + f_3 * pc_y[k] * smf_98[k];

        t_148[k] = f_3 * pc_y[k] * smf_99[k];

        t_149[k] = f_12 * slf_59[k]
                   + f_1 * smd0_59[k]
                   - f_2 * smd1_59[k]
                   + f_3 * pc_z[k] * smf_99[k];

        t_150[k] = f_13 * slf_100[k]
                   + f_1 * smd0_60[k]
                   - f_2 * smd1_60[k]
                   + f_3 * pc_x[k] * smf_100[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, pc_x, pc_y, pc_z, slf_60, slf_62, \
                         slf_103, smd0_63, smd1_63, smf_100, smf_102, \
                         smf_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_14 * slf_60[k]
                   + f_3 * pc_y[k] * smf_100[k];

        t_152[k] = f_3 * pc_z[k] * smf_100[k];

        t_153[k] = f_13 * slf_103[k]
                   + f_4 * smd0_63[k]
                   - f_5 * smd1_63[k]
                   + f_3 * pc_x[k] * smf_103[k];

        t_154[k] = f_14 * slf_62[k]
                   + f_3 * pc_y[k] * smf_102[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, pc_x, slf_105, slf_106, slf_107, slf_108, \
                         smd0_65, smd1_65, smf_105, smf_106, smf_107, \
                         smf_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_13 * slf_105[k]
                   + f_4 * smd0_65[k]
                   - f_5 * smd1_65[k]
                   + f_3 * pc_x[k] * smf_105[k];

        t_156[k] = f_13 * slf_106[k]
                   + f_3 * pc_x[k] * smf_106[k];

        t_157[k] = f_13 * slf_107[k]
                   + f_3 * pc_x[k] * smf_107[k];

        t_158[k] = f_13 * slf_108[k]
                   + f_3 * pc_x[k] * smf_108[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, pc_x, pc_y, pc_z, slf_66, slf_109, smd0_63, \
                         smd1_63, smf_106, smf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_13 * slf_109[k]
                   + f_3 * pc_x[k] * smf_109[k];

        t_160[k] = f_14 * slf_66[k]
                   + f_1 * smd0_63[k]
                   - f_2 * smd1_63[k]
                   + f_3 * pc_y[k] * smf_106[k];

        t_161[k] = f_3 * pc_z[k] * smf_106[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, pb_z, pc_y, pc_z, slg0_90, slf_68, \
                         slf_69, slg1_90, smd0_65, smd1_65, smf_108, \
                         smf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_14 * slf_68[k]
                   + f_4 * smd0_65[k]
                   - f_5 * smd1_65[k]
                   + f_3 * pc_y[k] * smf_108[k];

        t_163[k] = f_14 * slf_69[k]
                   + f_3 * pc_y[k] * smf_109[k];

        t_164[k] = f_1 * smd0_65[k]
                   - f_2 * smd1_65[k]
                   + f_3 * pc_z[k] * smf_109[k];

        t_165[k] = pb_z[k] * slg0_90[k]
                   - f_6 * pc_z[k] * slg1_90[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, pb_z, pc_y, pc_z, slg0_93, slf_60, \
                         slf_70, slf_72, slg1_93, smf_110, smf_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_12 * slf_70[k]
                   + f_3 * pc_y[k] * smf_110[k];

        t_167[k] = f_7 * slf_60[k]
                   + f_3 * pc_z[k] * smf_110[k];

        t_168[k] = pb_z[k] * slg0_93[k]
                   - f_6 * pc_z[k] * slg1_93[k];

        t_169[k] = f_12 * slf_72[k]
                   + f_3 * pc_y[k] * smf_112[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, pc_x, slf_115, slf_116, slf_117, slf_118, \
                         smd0_71, smd1_71, smf_115, smf_116, smf_117, \
                         smf_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_13 * slf_115[k]
                   + f_4 * smd0_71[k]
                   - f_5 * smd1_71[k]
                   + f_3 * pc_x[k] * smf_115[k];

        t_171[k] = f_13 * slf_116[k]
                   + f_3 * pc_x[k] * smf_116[k];

        t_172[k] = f_13 * slf_117[k]
                   + f_3 * pc_x[k] * smf_117[k];

        t_173[k] = f_13 * slf_118[k]
                   + f_3 * pc_x[k] * smf_118[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, pb_z, pc_x, pc_z, slg0_100, slf_66, slf_119, \
                         slg1_100, smf_116, smf_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_13 * slf_119[k]
                   + f_3 * pc_x[k] * smf_119[k];

        t_175[k] = pb_z[k] * slg0_100[k]
                   - f_6 * pc_z[k] * slg1_100[k];

        t_176[k] = f_7 * slf_66[k]
                   + f_3 * pc_z[k] * smf_116[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, pc_y, pc_z, slf_69, slf_78, slf_79, smd0_71, \
                         smd1_71, smf_118, smf_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = f_12 * slf_78[k]
                   + f_4 * smd0_71[k]
                   - f_5 * smd1_71[k]
                   + f_3 * pc_y[k] * smf_118[k];

        t_178[k] = f_12 * slf_79[k]
                   + f_3 * pc_y[k] * smf_119[k];

        t_179[k] = f_7 * slf_69[k]
                   + f_1 * smd0_71[k]
                   - f_2 * smd1_71[k]
                   + f_3 * pc_z[k] * smf_119[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, pc_x, pc_y, pc_z, slf_70, slf_80, slf_120, \
                         smd0_72, smd1_72, smf_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_13 * slf_120[k]
                   + f_1 * smd0_72[k]
                   - f_2 * smd1_72[k]
                   + f_3 * pc_x[k] * smf_120[k];

        t_181[k] = f_8 * slf_80[k]
                   + f_3 * pc_y[k] * smf_120[k];

        t_182[k] = f_8 * slf_70[k]
                   + f_3 * pc_z[k] * smf_120[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, pc_x, pc_y, slf_82, slf_123, slf_125, smd0_75, \
                         smd0_77, smd1_75, smd1_77, smf_122, smf_123, \
                         smf_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_13 * slf_123[k]
                   + f_4 * smd0_75[k]
                   - f_5 * smd1_75[k]
                   + f_3 * pc_x[k] * smf_123[k];

        t_184[k] = f_8 * slf_82[k]
                   + f_3 * pc_y[k] * smf_122[k];

        t_185[k] = f_13 * slf_125[k]
                   + f_4 * smd0_77[k]
                   - f_5 * smd1_77[k]
                   + f_3 * pc_x[k] * smf_125[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, t_189, pc_x, slf_126, slf_127, slf_128, slf_129, \
                         smf_126, smf_127, smf_128, smf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_13 * slf_126[k]
                   + f_3 * pc_x[k] * smf_126[k];

        t_187[k] = f_13 * slf_127[k]
                   + f_3 * pc_x[k] * smf_127[k];

        t_188[k] = f_13 * slf_128[k]
                   + f_3 * pc_x[k] * smf_128[k];

        t_189[k] = f_13 * slf_129[k]
                   + f_3 * pc_x[k] * smf_129[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, pc_y, pc_z, slf_76, slf_86, slf_88, smd0_75, \
                         smd0_77, smd1_75, smd1_77, smf_126, smf_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = f_8 * slf_86[k]
                   + f_1 * smd0_75[k]
                   - f_2 * smd1_75[k]
                   + f_3 * pc_y[k] * smf_126[k];

        t_191[k] = f_8 * slf_76[k]
                   + f_3 * pc_z[k] * smf_126[k];

        t_192[k] = f_8 * slf_88[k]
                   + f_4 * smd0_77[k]
                   - f_5 * smd1_77[k]
                   + f_3 * pc_y[k] * smf_128[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, pb_y, pc_y, pc_z, slg0_135, slf_79, \
                         slf_89, slf_90, slg1_135, smd0_77, smd1_77, smf_129, \
                         smf_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_8 * slf_89[k]
                   + f_3 * pc_y[k] * smf_129[k];

        t_194[k] = f_8 * slf_79[k]
                   + f_1 * smd0_77[k]
                   - f_2 * smd1_77[k]
                   + f_3 * pc_z[k] * smf_129[k];

        t_195[k] = pb_y[k] * slg0_135[k]
                   - f_6 * pc_y[k] * slg1_135[k];

        t_196[k] = f_7 * slf_90[k]
                   + f_3 * pc_y[k] * smf_130[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, pb_y, pc_y, pc_z, slg0_138, slg0_140, \
                         slf_80, slf_91, slf_92, slg1_138, slg1_140, smf_130, \
                         smf_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_12 * slf_80[k]
                   + f_3 * pc_z[k] * smf_130[k];

        t_198[k] = pb_y[k] * slg0_138[k]
                   + f_8 * slf_91[k]
                   - f_6 * pc_y[k] * slg1_138[k];

        t_199[k] = f_7 * slf_92[k]
                   + f_3 * pc_y[k] * smf_132[k];

        t_200[k] = pb_y[k] * slg0_140[k]
                   - f_6 * pc_y[k] * slg1_140[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, t_204, pc_x, slf_136, slf_137, slf_138, slf_139, \
                         smf_136, smf_137, smf_138, smf_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = f_13 * slf_136[k]
                   + f_3 * pc_x[k] * smf_136[k];

        t_202[k] = f_13 * slf_137[k]
                   + f_3 * pc_x[k] * smf_137[k];

        t_203[k] = f_13 * slf_138[k]
                   + f_3 * pc_x[k] * smf_138[k];

        t_204[k] = f_13 * slf_139[k]
                   + f_3 * pc_x[k] * smf_139[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, pc_y, pc_z, slf_86, slf_96, slf_98, smd0_81, \
                         smd0_83, smd1_81, smd1_83, smf_136, smf_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = f_7 * slf_96[k]
                   + f_1 * smd0_81[k]
                   - f_2 * smd1_81[k]
                   + f_3 * pc_y[k] * smf_136[k];

        t_206[k] = f_12 * slf_86[k]
                   + f_3 * pc_z[k] * smf_136[k];

        t_207[k] = f_7 * slf_98[k]
                   + f_4 * smd0_83[k]
                   - f_5 * smd1_83[k]
                   + f_3 * pc_y[k] * smf_138[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, pb_y, pc_x, pc_y, slg0_149, slf_99, \
                         slf_140, slg1_149, smd0_84, smd1_84, smf_139, \
                         smf_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_7 * slf_99[k]
                   + f_3 * pc_y[k] * smf_139[k];

        t_209[k] = pb_y[k] * slg0_149[k]
                   - f_6 * pc_y[k] * slg1_149[k];

        t_210[k] = f_13 * slf_140[k]
                   + f_1 * smd0_84[k]
                   - f_2 * smd1_84[k]
                   + f_3 * pc_x[k] * smf_140[k];

        t_211[k] = f_3 * pc_y[k] * smf_140[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, pc_x, pc_y, pc_z, slf_90, slf_143, smd0_87, \
                         smd1_87, smf_140, smf_142, smf_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_14 * slf_90[k]
                   + f_3 * pc_z[k] * smf_140[k];

        t_213[k] = f_13 * slf_143[k]
                   + f_4 * smd0_87[k]
                   - f_5 * smd1_87[k]
                   + f_3 * pc_x[k] * smf_143[k];

        t_214[k] = f_3 * pc_y[k] * smf_142[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, pc_x, slf_145, slf_146, slf_147, slf_148, \
                         smd0_89, smd1_89, smf_145, smf_146, smf_147, \
                         smf_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = f_13 * slf_145[k]
                   + f_4 * smd0_89[k]
                   - f_5 * smd1_89[k]
                   + f_3 * pc_x[k] * smf_145[k];

        t_216[k] = f_13 * slf_146[k]
                   + f_3 * pc_x[k] * smf_146[k];

        t_217[k] = f_13 * slf_147[k]
                   + f_3 * pc_x[k] * smf_147[k];

        t_218[k] = f_13 * slf_148[k]
                   + f_3 * pc_x[k] * smf_148[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, t_222, pc_x, pc_y, pc_z, slf_96, slf_149, \
                         smd0_87, smd0_89, smd1_87, smd1_89, smf_146, smf_148, \
                         smf_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = f_13 * slf_149[k]
                   + f_3 * pc_x[k] * smf_149[k];

        t_220[k] = f_1 * smd0_87[k]
                   - f_2 * smd1_87[k]
                   + f_3 * pc_y[k] * smf_146[k];

        t_221[k] = f_14 * slf_96[k]
                   + f_3 * pc_z[k] * smf_146[k];

        t_222[k] = f_4 * smd0_89[k]
                   - f_5 * smd1_89[k]
                   + f_3 * pc_y[k] * smf_148[k];
    }

#pragma omp simd aligned(t_223, t_224, t_225, t_226, pc_x, pc_y, pc_z, slf_99, slf_100, \
                         slf_150, smd0_89, smd0_90, smd1_89, smd1_90, smf_149, \
                         smf_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = f_3 * pc_y[k] * smf_149[k];

        t_224[k] = f_14 * slf_99[k]
                   + f_1 * smd0_89[k]
                   - f_2 * smd1_89[k]
                   + f_3 * pc_z[k] * smf_149[k];

        t_225[k] = f_14 * slf_150[k]
                   + f_1 * smd0_90[k]
                   - f_2 * smd1_90[k]
                   + f_3 * pc_x[k] * smf_150[k];

        t_226[k] = f_13 * slf_100[k]
                   + f_3 * pc_y[k] * smf_150[k];
    }

#pragma omp simd aligned(t_227, t_228, t_229, pc_x, pc_y, pc_z, slf_102, slf_153, smd0_93, \
                         smd1_93, smf_150, smf_152, smf_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_227[k] = f_3 * pc_z[k] * smf_150[k];

        t_228[k] = f_14 * slf_153[k]
                   + f_4 * smd0_93[k]
                   - f_5 * smd1_93[k]
                   + f_3 * pc_x[k] * smf_153[k];

        t_229[k] = f_13 * slf_102[k]
                   + f_3 * pc_y[k] * smf_152[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, pc_x, slf_155, slf_156, slf_157, slf_158, \
                         smd0_95, smd1_95, smf_155, smf_156, smf_157, \
                         smf_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_14 * slf_155[k]
                   + f_4 * smd0_95[k]
                   - f_5 * smd1_95[k]
                   + f_3 * pc_x[k] * smf_155[k];

        t_231[k] = f_14 * slf_156[k]
                   + f_3 * pc_x[k] * smf_156[k];

        t_232[k] = f_14 * slf_157[k]
                   + f_3 * pc_x[k] * smf_157[k];

        t_233[k] = f_14 * slf_158[k]
                   + f_3 * pc_x[k] * smf_158[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, pc_x, pc_y, pc_z, slf_106, slf_159, smd0_93, \
                         smd1_93, smf_156, smf_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_14 * slf_159[k]
                   + f_3 * pc_x[k] * smf_159[k];

        t_235[k] = f_13 * slf_106[k]
                   + f_1 * smd0_93[k]
                   - f_2 * smd1_93[k]
                   + f_3 * pc_y[k] * smf_156[k];

        t_236[k] = f_3 * pc_z[k] * smf_156[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, t_240, pb_z, pc_y, pc_z, slg0_150, slf_108, \
                         slf_109, slg1_150, smd0_95, smd1_95, smf_158, \
                         smf_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = f_13 * slf_108[k]
                   + f_4 * smd0_95[k]
                   - f_5 * smd1_95[k]
                   + f_3 * pc_y[k] * smf_158[k];

        t_238[k] = f_13 * slf_109[k]
                   + f_3 * pc_y[k] * smf_159[k];

        t_239[k] = f_1 * smd0_95[k]
                   - f_2 * smd1_95[k]
                   + f_3 * pc_z[k] * smf_159[k];

        t_240[k] = pb_z[k] * slg0_150[k]
                   - f_6 * pc_z[k] * slg1_150[k];
    }

#pragma omp simd aligned(t_241, t_242, t_243, t_244, pb_z, pc_y, pc_z, slg0_153, slf_100, \
                         slf_110, slf_112, slg1_153, smf_160, smf_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_241[k] = f_14 * slf_110[k]
                   + f_3 * pc_y[k] * smf_160[k];

        t_242[k] = f_7 * slf_100[k]
                   + f_3 * pc_z[k] * smf_160[k];

        t_243[k] = pb_z[k] * slg0_153[k]
                   - f_6 * pc_z[k] * slg1_153[k];

        t_244[k] = f_14 * slf_112[k]
                   + f_3 * pc_y[k] * smf_162[k];
    }
}

static auto
compute_prim_smg_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t slg0,
                                                          const size_t slf, const size_t slg1,
                                                          const size_t smd0, const size_t smd1,
                                                          const size_t smf, const size_t ncols,
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
    const auto f_11 = 3.0 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 2.5 / q;
    const auto f_14 = 2.0 / q;

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

    const auto *slg0_160 = buffer.data(slg0 + 160);
    const auto *slg0_210 = buffer.data(slg0 + 210);
    const auto *slg0_213 = buffer.data(slg0 + 213);
    const auto *slg0_215 = buffer.data(slg0 + 215);
    const auto *slg0_224 = buffer.data(slg0 + 224);
    const auto *slg0_225 = buffer.data(slg0 + 225);
    const auto *slg0_228 = buffer.data(slg0 + 228);
    const auto *slg0_235 = buffer.data(slg0 + 235);

    const auto *slf_106 = buffer.data(slf + 106);
    const auto *slf_109 = buffer.data(slf + 109);
    const auto *slf_110 = buffer.data(slf + 110);
    const auto *slf_116 = buffer.data(slf + 116);
    const auto *slf_118 = buffer.data(slf + 118);
    const auto *slf_119 = buffer.data(slf + 119);
    const auto *slf_120 = buffer.data(slf + 120);
    const auto *slf_122 = buffer.data(slf + 122);
    const auto *slf_126 = buffer.data(slf + 126);
    const auto *slf_128 = buffer.data(slf + 128);
    const auto *slf_129 = buffer.data(slf + 129);
    const auto *slf_130 = buffer.data(slf + 130);
    const auto *slf_132 = buffer.data(slf + 132);
    const auto *slf_136 = buffer.data(slf + 136);
    const auto *slf_138 = buffer.data(slf + 138);
    const auto *slf_139 = buffer.data(slf + 139);
    const auto *slf_140 = buffer.data(slf + 140);
    const auto *slf_141 = buffer.data(slf + 141);
    const auto *slf_142 = buffer.data(slf + 142);
    const auto *slf_146 = buffer.data(slf + 146);
    const auto *slf_148 = buffer.data(slf + 148);
    const auto *slf_149 = buffer.data(slf + 149);
    const auto *slf_150 = buffer.data(slf + 150);
    const auto *slf_152 = buffer.data(slf + 152);
    const auto *slf_156 = buffer.data(slf + 156);
    const auto *slf_158 = buffer.data(slf + 158);
    const auto *slf_159 = buffer.data(slf + 159);
    const auto *slf_160 = buffer.data(slf + 160);
    const auto *slf_162 = buffer.data(slf + 162);
    const auto *slf_165 = buffer.data(slf + 165);
    const auto *slf_166 = buffer.data(slf + 166);
    const auto *slf_167 = buffer.data(slf + 167);
    const auto *slf_168 = buffer.data(slf + 168);
    const auto *slf_169 = buffer.data(slf + 169);
    const auto *slf_170 = buffer.data(slf + 170);
    const auto *slf_172 = buffer.data(slf + 172);
    const auto *slf_173 = buffer.data(slf + 173);
    const auto *slf_175 = buffer.data(slf + 175);
    const auto *slf_176 = buffer.data(slf + 176);
    const auto *slf_177 = buffer.data(slf + 177);
    const auto *slf_178 = buffer.data(slf + 178);
    const auto *slf_179 = buffer.data(slf + 179);
    const auto *slf_180 = buffer.data(slf + 180);
    const auto *slf_183 = buffer.data(slf + 183);
    const auto *slf_185 = buffer.data(slf + 185);
    const auto *slf_186 = buffer.data(slf + 186);
    const auto *slf_187 = buffer.data(slf + 187);
    const auto *slf_188 = buffer.data(slf + 188);
    const auto *slf_189 = buffer.data(slf + 189);
    const auto *slf_196 = buffer.data(slf + 196);
    const auto *slf_197 = buffer.data(slf + 197);
    const auto *slf_198 = buffer.data(slf + 198);
    const auto *slf_199 = buffer.data(slf + 199);
    const auto *slf_200 = buffer.data(slf + 200);
    const auto *slf_203 = buffer.data(slf + 203);
    const auto *slf_205 = buffer.data(slf + 205);
    const auto *slf_206 = buffer.data(slf + 206);
    const auto *slf_207 = buffer.data(slf + 207);
    const auto *slf_208 = buffer.data(slf + 208);
    const auto *slf_209 = buffer.data(slf + 209);
    const auto *slf_210 = buffer.data(slf + 210);
    const auto *slf_213 = buffer.data(slf + 213);
    const auto *slf_215 = buffer.data(slf + 215);
    const auto *slf_216 = buffer.data(slf + 216);
    const auto *slf_217 = buffer.data(slf + 217);
    const auto *slf_218 = buffer.data(slf + 218);
    const auto *slf_219 = buffer.data(slf + 219);
    const auto *slf_225 = buffer.data(slf + 225);
    const auto *slf_226 = buffer.data(slf + 226);
    const auto *slf_227 = buffer.data(slf + 227);
    const auto *slf_228 = buffer.data(slf + 228);
    const auto *slf_229 = buffer.data(slf + 229);
    const auto *slf_230 = buffer.data(slf + 230);
    const auto *slf_233 = buffer.data(slf + 233);
    const auto *slf_235 = buffer.data(slf + 235);
    const auto *slf_236 = buffer.data(slf + 236);
    const auto *slf_237 = buffer.data(slf + 237);
    const auto *slf_238 = buffer.data(slf + 238);
    const auto *slf_239 = buffer.data(slf + 239);
    const auto *slf_240 = buffer.data(slf + 240);

    const auto *slg1_160 = buffer.data(slg1 + 160);
    const auto *slg1_210 = buffer.data(slg1 + 210);
    const auto *slg1_213 = buffer.data(slg1 + 213);
    const auto *slg1_215 = buffer.data(slg1 + 215);
    const auto *slg1_224 = buffer.data(slg1 + 224);
    const auto *slg1_225 = buffer.data(slg1 + 225);
    const auto *slg1_228 = buffer.data(slg1 + 228);
    const auto *slg1_235 = buffer.data(slg1 + 235);

    const auto *smd0_101 = buffer.data(smd0 + 101);
    const auto *smd0_102 = buffer.data(smd0 + 102);
    const auto *smd0_105 = buffer.data(smd0 + 105);
    const auto *smd0_107 = buffer.data(smd0 + 107);
    const auto *smd0_108 = buffer.data(smd0 + 108);
    const auto *smd0_111 = buffer.data(smd0 + 111);
    const auto *smd0_113 = buffer.data(smd0 + 113);
    const auto *smd0_117 = buffer.data(smd0 + 117);
    const auto *smd0_119 = buffer.data(smd0 + 119);
    const auto *smd0_120 = buffer.data(smd0 + 120);
    const auto *smd0_123 = buffer.data(smd0 + 123);
    const auto *smd0_125 = buffer.data(smd0 + 125);
    const auto *smd0_126 = buffer.data(smd0 + 126);
    const auto *smd0_129 = buffer.data(smd0 + 129);
    const auto *smd0_131 = buffer.data(smd0 + 131);
    const auto *smd0_137 = buffer.data(smd0 + 137);
    const auto *smd0_138 = buffer.data(smd0 + 138);
    const auto *smd0_141 = buffer.data(smd0 + 141);
    const auto *smd0_143 = buffer.data(smd0 + 143);
    const auto *smd0_144 = buffer.data(smd0 + 144);

    const auto *smd1_101 = buffer.data(smd1 + 101);
    const auto *smd1_102 = buffer.data(smd1 + 102);
    const auto *smd1_105 = buffer.data(smd1 + 105);
    const auto *smd1_107 = buffer.data(smd1 + 107);
    const auto *smd1_108 = buffer.data(smd1 + 108);
    const auto *smd1_111 = buffer.data(smd1 + 111);
    const auto *smd1_113 = buffer.data(smd1 + 113);
    const auto *smd1_117 = buffer.data(smd1 + 117);
    const auto *smd1_119 = buffer.data(smd1 + 119);
    const auto *smd1_120 = buffer.data(smd1 + 120);
    const auto *smd1_123 = buffer.data(smd1 + 123);
    const auto *smd1_125 = buffer.data(smd1 + 125);
    const auto *smd1_126 = buffer.data(smd1 + 126);
    const auto *smd1_129 = buffer.data(smd1 + 129);
    const auto *smd1_131 = buffer.data(smd1 + 131);
    const auto *smd1_137 = buffer.data(smd1 + 137);
    const auto *smd1_138 = buffer.data(smd1 + 138);
    const auto *smd1_141 = buffer.data(smd1 + 141);
    const auto *smd1_143 = buffer.data(smd1 + 143);
    const auto *smd1_144 = buffer.data(smd1 + 144);

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
    const auto *smf_182 = buffer.data(smf + 182);
    const auto *smf_183 = buffer.data(smf + 183);
    const auto *smf_185 = buffer.data(smf + 185);
    const auto *smf_186 = buffer.data(smf + 186);
    const auto *smf_187 = buffer.data(smf + 187);
    const auto *smf_188 = buffer.data(smf + 188);
    const auto *smf_189 = buffer.data(smf + 189);
    const auto *smf_190 = buffer.data(smf + 190);
    const auto *smf_192 = buffer.data(smf + 192);
    const auto *smf_196 = buffer.data(smf + 196);
    const auto *smf_197 = buffer.data(smf + 197);
    const auto *smf_198 = buffer.data(smf + 198);
    const auto *smf_199 = buffer.data(smf + 199);
    const auto *smf_200 = buffer.data(smf + 200);
    const auto *smf_202 = buffer.data(smf + 202);
    const auto *smf_203 = buffer.data(smf + 203);
    const auto *smf_205 = buffer.data(smf + 205);
    const auto *smf_206 = buffer.data(smf + 206);
    const auto *smf_207 = buffer.data(smf + 207);
    const auto *smf_208 = buffer.data(smf + 208);
    const auto *smf_209 = buffer.data(smf + 209);
    const auto *smf_210 = buffer.data(smf + 210);
    const auto *smf_212 = buffer.data(smf + 212);
    const auto *smf_213 = buffer.data(smf + 213);
    const auto *smf_215 = buffer.data(smf + 215);
    const auto *smf_216 = buffer.data(smf + 216);
    const auto *smf_217 = buffer.data(smf + 217);
    const auto *smf_218 = buffer.data(smf + 218);
    const auto *smf_219 = buffer.data(smf + 219);
    const auto *smf_220 = buffer.data(smf + 220);
    const auto *smf_222 = buffer.data(smf + 222);
    const auto *smf_225 = buffer.data(smf + 225);
    const auto *smf_226 = buffer.data(smf + 226);
    const auto *smf_227 = buffer.data(smf + 227);
    const auto *smf_228 = buffer.data(smf + 228);
    const auto *smf_229 = buffer.data(smf + 229);
    const auto *smf_230 = buffer.data(smf + 230);
    const auto *smf_232 = buffer.data(smf + 232);
    const auto *smf_233 = buffer.data(smf + 233);
    const auto *smf_235 = buffer.data(smf + 235);
    const auto *smf_236 = buffer.data(smf + 236);
    const auto *smf_237 = buffer.data(smf + 237);
    const auto *smf_238 = buffer.data(smf + 238);
    const auto *smf_239 = buffer.data(smf + 239);
    const auto *smf_240 = buffer.data(smf + 240);

#pragma omp simd aligned(t_245, t_246, t_247, t_248, pc_x, slf_165, slf_166, slf_167, slf_168, \
                         smd0_101, smd1_101, smf_165, smf_166, smf_167, \
                         smf_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_14 * slf_165[k]
                   + f_4 * smd0_101[k]
                   - f_5 * smd1_101[k]
                   + f_3 * pc_x[k] * smf_165[k];

        t_246[k] = f_14 * slf_166[k]
                   + f_3 * pc_x[k] * smf_166[k];

        t_247[k] = f_14 * slf_167[k]
                   + f_3 * pc_x[k] * smf_167[k];

        t_248[k] = f_14 * slf_168[k]
                   + f_3 * pc_x[k] * smf_168[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, pb_z, pc_x, pc_z, slg0_160, slf_106, slf_169, \
                         slg1_160, smf_166, smf_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = f_14 * slf_169[k]
                   + f_3 * pc_x[k] * smf_169[k];

        t_250[k] = pb_z[k] * slg0_160[k]
                   - f_6 * pc_z[k] * slg1_160[k];

        t_251[k] = f_7 * slf_106[k]
                   + f_3 * pc_z[k] * smf_166[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, pc_y, pc_z, slf_109, slf_118, slf_119, smd0_101, \
                         smd1_101, smf_168, smf_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = f_14 * slf_118[k]
                   + f_4 * smd0_101[k]
                   - f_5 * smd1_101[k]
                   + f_3 * pc_y[k] * smf_168[k];

        t_253[k] = f_14 * slf_119[k]
                   + f_3 * pc_y[k] * smf_169[k];

        t_254[k] = f_7 * slf_109[k]
                   + f_1 * smd0_101[k]
                   - f_2 * smd1_101[k]
                   + f_3 * pc_z[k] * smf_169[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, pc_x, pc_y, pc_z, slf_110, slf_120, slf_170, \
                         smd0_102, smd1_102, smf_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_14 * slf_170[k]
                   + f_1 * smd0_102[k]
                   - f_2 * smd1_102[k]
                   + f_3 * pc_x[k] * smf_170[k];

        t_256[k] = f_12 * slf_120[k]
                   + f_3 * pc_y[k] * smf_170[k];

        t_257[k] = f_8 * slf_110[k]
                   + f_3 * pc_z[k] * smf_170[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, pc_x, pc_y, slf_122, slf_173, slf_175, smd0_105, \
                         smd0_107, smd1_105, smd1_107, smf_172, smf_173, \
                         smf_175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = f_14 * slf_173[k]
                   + f_4 * smd0_105[k]
                   - f_5 * smd1_105[k]
                   + f_3 * pc_x[k] * smf_173[k];

        t_259[k] = f_12 * slf_122[k]
                   + f_3 * pc_y[k] * smf_172[k];

        t_260[k] = f_14 * slf_175[k]
                   + f_4 * smd0_107[k]
                   - f_5 * smd1_107[k]
                   + f_3 * pc_x[k] * smf_175[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, t_264, pc_x, slf_176, slf_177, slf_178, slf_179, \
                         smf_176, smf_177, smf_178, smf_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_14 * slf_176[k]
                   + f_3 * pc_x[k] * smf_176[k];

        t_262[k] = f_14 * slf_177[k]
                   + f_3 * pc_x[k] * smf_177[k];

        t_263[k] = f_14 * slf_178[k]
                   + f_3 * pc_x[k] * smf_178[k];

        t_264[k] = f_14 * slf_179[k]
                   + f_3 * pc_x[k] * smf_179[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, pc_y, pc_z, slf_116, slf_126, slf_128, smd0_105, \
                         smd0_107, smd1_105, smd1_107, smf_176, \
                         smf_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = f_12 * slf_126[k]
                   + f_1 * smd0_105[k]
                   - f_2 * smd1_105[k]
                   + f_3 * pc_y[k] * smf_176[k];

        t_266[k] = f_8 * slf_116[k]
                   + f_3 * pc_z[k] * smf_176[k];

        t_267[k] = f_12 * slf_128[k]
                   + f_4 * smd0_107[k]
                   - f_5 * smd1_107[k]
                   + f_3 * pc_y[k] * smf_178[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, pc_x, pc_y, pc_z, slf_119, slf_129, slf_180, \
                         smd0_107, smd0_108, smd1_107, smd1_108, smf_179, \
                         smf_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_12 * slf_129[k]
                   + f_3 * pc_y[k] * smf_179[k];

        t_269[k] = f_8 * slf_119[k]
                   + f_1 * smd0_107[k]
                   - f_2 * smd1_107[k]
                   + f_3 * pc_z[k] * smf_179[k];

        t_270[k] = f_14 * slf_180[k]
                   + f_1 * smd0_108[k]
                   - f_2 * smd1_108[k]
                   + f_3 * pc_x[k] * smf_180[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, t_274, pc_x, pc_y, pc_z, slf_120, slf_130, \
                         slf_132, slf_183, smd0_111, smd1_111, smf_180, smf_182, \
                         smf_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = f_8 * slf_130[k]
                   + f_3 * pc_y[k] * smf_180[k];

        t_272[k] = f_12 * slf_120[k]
                   + f_3 * pc_z[k] * smf_180[k];

        t_273[k] = f_14 * slf_183[k]
                   + f_4 * smd0_111[k]
                   - f_5 * smd1_111[k]
                   + f_3 * pc_x[k] * smf_183[k];

        t_274[k] = f_8 * slf_132[k]
                   + f_3 * pc_y[k] * smf_182[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, pc_x, slf_185, slf_186, slf_187, slf_188, \
                         smd0_113, smd1_113, smf_185, smf_186, smf_187, \
                         smf_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = f_14 * slf_185[k]
                   + f_4 * smd0_113[k]
                   - f_5 * smd1_113[k]
                   + f_3 * pc_x[k] * smf_185[k];

        t_276[k] = f_14 * slf_186[k]
                   + f_3 * pc_x[k] * smf_186[k];

        t_277[k] = f_14 * slf_187[k]
                   + f_3 * pc_x[k] * smf_187[k];

        t_278[k] = f_14 * slf_188[k]
                   + f_3 * pc_x[k] * smf_188[k];
    }

#pragma omp simd aligned(t_279, t_280, t_281, pc_x, pc_y, pc_z, slf_126, slf_136, slf_189, \
                         smd0_111, smd1_111, smf_186, smf_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_279[k] = f_14 * slf_189[k]
                   + f_3 * pc_x[k] * smf_189[k];

        t_280[k] = f_8 * slf_136[k]
                   + f_1 * smd0_111[k]
                   - f_2 * smd1_111[k]
                   + f_3 * pc_y[k] * smf_186[k];

        t_281[k] = f_12 * slf_126[k]
                   + f_3 * pc_z[k] * smf_186[k];
    }

#pragma omp simd aligned(t_282, t_283, t_284, t_285, pb_y, pc_y, pc_z, slg0_210, slf_129, \
                         slf_138, slf_139, slg1_210, smd0_113, smd1_113, smf_188, \
                         smf_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_282[k] = f_8 * slf_138[k]
                   + f_4 * smd0_113[k]
                   - f_5 * smd1_113[k]
                   + f_3 * pc_y[k] * smf_188[k];

        t_283[k] = f_8 * slf_139[k]
                   + f_3 * pc_y[k] * smf_189[k];

        t_284[k] = f_12 * slf_129[k]
                   + f_1 * smd0_113[k]
                   - f_2 * smd1_113[k]
                   + f_3 * pc_z[k] * smf_189[k];

        t_285[k] = pb_y[k] * slg0_210[k]
                   - f_6 * pc_y[k] * slg1_210[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, pb_y, pc_y, pc_z, slg0_213, slf_130, \
                         slf_140, slf_141, slf_142, slg1_213, smf_190, \
                         smf_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = f_7 * slf_140[k]
                   + f_3 * pc_y[k] * smf_190[k];

        t_287[k] = f_14 * slf_130[k]
                   + f_3 * pc_z[k] * smf_190[k];

        t_288[k] = pb_y[k] * slg0_213[k]
                   + f_8 * slf_141[k]
                   - f_6 * pc_y[k] * slg1_213[k];

        t_289[k] = f_7 * slf_142[k]
                   + f_3 * pc_y[k] * smf_192[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, pb_y, pc_x, pc_y, slg0_215, slf_196, \
                         slf_197, slf_198, slg1_215, smf_196, smf_197, \
                         smf_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = pb_y[k] * slg0_215[k]
                   - f_6 * pc_y[k] * slg1_215[k];

        t_291[k] = f_14 * slf_196[k]
                   + f_3 * pc_x[k] * smf_196[k];

        t_292[k] = f_14 * slf_197[k]
                   + f_3 * pc_x[k] * smf_197[k];

        t_293[k] = f_14 * slf_198[k]
                   + f_3 * pc_x[k] * smf_198[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, pc_x, pc_y, pc_z, slf_136, slf_146, slf_199, \
                         smd0_117, smd1_117, smf_196, smf_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_14 * slf_199[k]
                   + f_3 * pc_x[k] * smf_199[k];

        t_295[k] = f_7 * slf_146[k]
                   + f_1 * smd0_117[k]
                   - f_2 * smd1_117[k]
                   + f_3 * pc_y[k] * smf_196[k];

        t_296[k] = f_14 * slf_136[k]
                   + f_3 * pc_z[k] * smf_196[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, pb_y, pc_y, slg0_224, slf_148, slf_149, \
                         slg1_224, smd0_119, smd1_119, smf_198, \
                         smf_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = f_7 * slf_148[k]
                   + f_4 * smd0_119[k]
                   - f_5 * smd1_119[k]
                   + f_3 * pc_y[k] * smf_198[k];

        t_298[k] = f_7 * slf_149[k]
                   + f_3 * pc_y[k] * smf_199[k];

        t_299[k] = pb_y[k] * slg0_224[k]
                   - f_6 * pc_y[k] * slg1_224[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, pc_x, pc_y, pc_z, slf_140, slf_200, \
                         slf_203, smd0_120, smd0_123, smd1_120, smd1_123, smf_200, \
                         smf_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = f_14 * slf_200[k]
                   + f_1 * smd0_120[k]
                   - f_2 * smd1_120[k]
                   + f_3 * pc_x[k] * smf_200[k];

        t_301[k] = f_3 * pc_y[k] * smf_200[k];

        t_302[k] = f_13 * slf_140[k]
                   + f_3 * pc_z[k] * smf_200[k];

        t_303[k] = f_14 * slf_203[k]
                   + f_4 * smd0_123[k]
                   - f_5 * smd1_123[k]
                   + f_3 * pc_x[k] * smf_203[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, pc_x, pc_y, slf_205, slf_206, slf_207, \
                         smd0_125, smd1_125, smf_202, smf_205, smf_206, \
                         smf_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_304[k] = f_3 * pc_y[k] * smf_202[k];

        t_305[k] = f_14 * slf_205[k]
                   + f_4 * smd0_125[k]
                   - f_5 * smd1_125[k]
                   + f_3 * pc_x[k] * smf_205[k];

        t_306[k] = f_14 * slf_206[k]
                   + f_3 * pc_x[k] * smf_206[k];

        t_307[k] = f_14 * slf_207[k]
                   + f_3 * pc_x[k] * smf_207[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, t_311, pc_x, pc_y, pc_z, slf_146, slf_208, \
                         slf_209, smd0_123, smd1_123, smf_206, smf_208, \
                         smf_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = f_14 * slf_208[k]
                   + f_3 * pc_x[k] * smf_208[k];

        t_309[k] = f_14 * slf_209[k]
                   + f_3 * pc_x[k] * smf_209[k];

        t_310[k] = f_1 * smd0_123[k]
                   - f_2 * smd1_123[k]
                   + f_3 * pc_y[k] * smf_206[k];

        t_311[k] = f_13 * slf_146[k]
                   + f_3 * pc_z[k] * smf_206[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, t_315, pc_x, pc_y, pc_z, slf_149, slf_210, \
                         smd0_125, smd0_126, smd1_125, smd1_126, smf_208, smf_209, \
                         smf_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = f_4 * smd0_125[k]
                   - f_5 * smd1_125[k]
                   + f_3 * pc_y[k] * smf_208[k];

        t_313[k] = f_3 * pc_y[k] * smf_209[k];

        t_314[k] = f_13 * slf_149[k]
                   + f_1 * smd0_125[k]
                   - f_2 * smd1_125[k]
                   + f_3 * pc_z[k] * smf_209[k];

        t_315[k] = f_12 * slf_210[k]
                   + f_1 * smd0_126[k]
                   - f_2 * smd1_126[k]
                   + f_3 * pc_x[k] * smf_210[k];
    }

#pragma omp simd aligned(t_316, t_317, t_318, t_319, pc_x, pc_y, pc_z, slf_150, slf_152, \
                         slf_213, smd0_129, smd1_129, smf_210, smf_212, \
                         smf_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_316[k] = f_11 * slf_150[k]
                   + f_3 * pc_y[k] * smf_210[k];

        t_317[k] = f_3 * pc_z[k] * smf_210[k];

        t_318[k] = f_12 * slf_213[k]
                   + f_4 * smd0_129[k]
                   - f_5 * smd1_129[k]
                   + f_3 * pc_x[k] * smf_213[k];

        t_319[k] = f_11 * slf_152[k]
                   + f_3 * pc_y[k] * smf_212[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, pc_x, slf_215, slf_216, slf_217, slf_218, \
                         smd0_131, smd1_131, smf_215, smf_216, smf_217, \
                         smf_218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = f_12 * slf_215[k]
                   + f_4 * smd0_131[k]
                   - f_5 * smd1_131[k]
                   + f_3 * pc_x[k] * smf_215[k];

        t_321[k] = f_12 * slf_216[k]
                   + f_3 * pc_x[k] * smf_216[k];

        t_322[k] = f_12 * slf_217[k]
                   + f_3 * pc_x[k] * smf_217[k];

        t_323[k] = f_12 * slf_218[k]
                   + f_3 * pc_x[k] * smf_218[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, pc_x, pc_y, pc_z, slf_156, slf_219, smd0_129, \
                         smd1_129, smf_216, smf_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = f_12 * slf_219[k]
                   + f_3 * pc_x[k] * smf_219[k];

        t_325[k] = f_11 * slf_156[k]
                   + f_1 * smd0_129[k]
                   - f_2 * smd1_129[k]
                   + f_3 * pc_y[k] * smf_216[k];

        t_326[k] = f_3 * pc_z[k] * smf_216[k];
    }

#pragma omp simd aligned(t_327, t_328, t_329, t_330, pb_z, pc_y, pc_z, slg0_225, slf_158, \
                         slf_159, slg1_225, smd0_131, smd1_131, smf_218, \
                         smf_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_327[k] = f_11 * slf_158[k]
                   + f_4 * smd0_131[k]
                   - f_5 * smd1_131[k]
                   + f_3 * pc_y[k] * smf_218[k];

        t_328[k] = f_11 * slf_159[k]
                   + f_3 * pc_y[k] * smf_219[k];

        t_329[k] = f_1 * smd0_131[k]
                   - f_2 * smd1_131[k]
                   + f_3 * pc_z[k] * smf_219[k];

        t_330[k] = pb_z[k] * slg0_225[k]
                   - f_6 * pc_z[k] * slg1_225[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, t_334, pb_z, pc_y, pc_z, slg0_228, slf_150, \
                         slf_160, slf_162, slg1_228, smf_220, smf_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = f_13 * slf_160[k]
                   + f_3 * pc_y[k] * smf_220[k];

        t_332[k] = f_7 * slf_150[k]
                   + f_3 * pc_z[k] * smf_220[k];

        t_333[k] = pb_z[k] * slg0_228[k]
                   - f_6 * pc_z[k] * slg1_228[k];

        t_334[k] = f_13 * slf_162[k]
                   + f_3 * pc_y[k] * smf_222[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, pc_x, slf_225, slf_226, slf_227, slf_228, \
                         smd0_137, smd1_137, smf_225, smf_226, smf_227, \
                         smf_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = f_12 * slf_225[k]
                   + f_4 * smd0_137[k]
                   - f_5 * smd1_137[k]
                   + f_3 * pc_x[k] * smf_225[k];

        t_336[k] = f_12 * slf_226[k]
                   + f_3 * pc_x[k] * smf_226[k];

        t_337[k] = f_12 * slf_227[k]
                   + f_3 * pc_x[k] * smf_227[k];

        t_338[k] = f_12 * slf_228[k]
                   + f_3 * pc_x[k] * smf_228[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, pb_z, pc_x, pc_z, slg0_235, slf_156, slf_229, \
                         slg1_235, smf_226, smf_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = f_12 * slf_229[k]
                   + f_3 * pc_x[k] * smf_229[k];

        t_340[k] = pb_z[k] * slg0_235[k]
                   - f_6 * pc_z[k] * slg1_235[k];

        t_341[k] = f_7 * slf_156[k]
                   + f_3 * pc_z[k] * smf_226[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, pc_y, pc_z, slf_159, slf_168, slf_169, smd0_137, \
                         smd1_137, smf_228, smf_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = f_13 * slf_168[k]
                   + f_4 * smd0_137[k]
                   - f_5 * smd1_137[k]
                   + f_3 * pc_y[k] * smf_228[k];

        t_343[k] = f_13 * slf_169[k]
                   + f_3 * pc_y[k] * smf_229[k];

        t_344[k] = f_7 * slf_159[k]
                   + f_1 * smd0_137[k]
                   - f_2 * smd1_137[k]
                   + f_3 * pc_z[k] * smf_229[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, pc_x, pc_y, pc_z, slf_160, slf_170, slf_230, \
                         smd0_138, smd1_138, smf_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_12 * slf_230[k]
                   + f_1 * smd0_138[k]
                   - f_2 * smd1_138[k]
                   + f_3 * pc_x[k] * smf_230[k];

        t_346[k] = f_14 * slf_170[k]
                   + f_3 * pc_y[k] * smf_230[k];

        t_347[k] = f_8 * slf_160[k]
                   + f_3 * pc_z[k] * smf_230[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, pc_x, pc_y, slf_172, slf_233, slf_235, smd0_141, \
                         smd0_143, smd1_141, smd1_143, smf_232, smf_233, \
                         smf_235 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = f_12 * slf_233[k]
                   + f_4 * smd0_141[k]
                   - f_5 * smd1_141[k]
                   + f_3 * pc_x[k] * smf_233[k];

        t_349[k] = f_14 * slf_172[k]
                   + f_3 * pc_y[k] * smf_232[k];

        t_350[k] = f_12 * slf_235[k]
                   + f_4 * smd0_143[k]
                   - f_5 * smd1_143[k]
                   + f_3 * pc_x[k] * smf_235[k];
    }

#pragma omp simd aligned(t_351, t_352, t_353, t_354, pc_x, slf_236, slf_237, slf_238, slf_239, \
                         smf_236, smf_237, smf_238, smf_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_351[k] = f_12 * slf_236[k]
                   + f_3 * pc_x[k] * smf_236[k];

        t_352[k] = f_12 * slf_237[k]
                   + f_3 * pc_x[k] * smf_237[k];

        t_353[k] = f_12 * slf_238[k]
                   + f_3 * pc_x[k] * smf_238[k];

        t_354[k] = f_12 * slf_239[k]
                   + f_3 * pc_x[k] * smf_239[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, pc_y, pc_z, slf_166, slf_176, slf_178, smd0_141, \
                         smd0_143, smd1_141, smd1_143, smf_236, \
                         smf_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = f_14 * slf_176[k]
                   + f_1 * smd0_141[k]
                   - f_2 * smd1_141[k]
                   + f_3 * pc_y[k] * smf_236[k];

        t_356[k] = f_8 * slf_166[k]
                   + f_3 * pc_z[k] * smf_236[k];

        t_357[k] = f_14 * slf_178[k]
                   + f_4 * smd0_143[k]
                   - f_5 * smd1_143[k]
                   + f_3 * pc_y[k] * smf_238[k];
    }

#pragma omp simd aligned(t_358, t_359, t_360, pc_x, pc_y, pc_z, slf_169, slf_179, slf_240, \
                         smd0_143, smd0_144, smd1_143, smd1_144, smf_239, \
                         smf_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_358[k] = f_14 * slf_179[k]
                   + f_3 * pc_y[k] * smf_239[k];

        t_359[k] = f_8 * slf_169[k]
                   + f_1 * smd0_143[k]
                   - f_2 * smd1_143[k]
                   + f_3 * pc_z[k] * smf_239[k];

        t_360[k] = f_12 * slf_240[k]
                   + f_1 * smd0_144[k]
                   - f_2 * smd1_144[k]
                   + f_3 * pc_x[k] * smf_240[k];
    }
}

static auto
compute_prim_smg_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t slg0,
                                                          const size_t slf, const size_t slg1,
                                                          const size_t smd0, const size_t smd1,
                                                          const size_t smf, const size_t ncols,
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
    const auto f_10 = 3.5 / q;
    const auto f_11 = 3.0 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 2.5 / q;
    const auto f_14 = 2.0 / q;

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

    const auto *slg0_300 = buffer.data(slg0 + 300);
    const auto *slg0_303 = buffer.data(slg0 + 303);
    const auto *slg0_305 = buffer.data(slg0 + 305);
    const auto *slg0_314 = buffer.data(slg0 + 314);
    const auto *slg0_315 = buffer.data(slg0 + 315);
    const auto *slg0_318 = buffer.data(slg0 + 318);
    const auto *slg0_325 = buffer.data(slg0 + 325);

    const auto *slf_170 = buffer.data(slf + 170);
    const auto *slf_176 = buffer.data(slf + 176);
    const auto *slf_179 = buffer.data(slf + 179);
    const auto *slf_180 = buffer.data(slf + 180);
    const auto *slf_182 = buffer.data(slf + 182);
    const auto *slf_186 = buffer.data(slf + 186);
    const auto *slf_188 = buffer.data(slf + 188);
    const auto *slf_189 = buffer.data(slf + 189);
    const auto *slf_190 = buffer.data(slf + 190);
    const auto *slf_192 = buffer.data(slf + 192);
    const auto *slf_196 = buffer.data(slf + 196);
    const auto *slf_198 = buffer.data(slf + 198);
    const auto *slf_199 = buffer.data(slf + 199);
    const auto *slf_200 = buffer.data(slf + 200);
    const auto *slf_201 = buffer.data(slf + 201);
    const auto *slf_202 = buffer.data(slf + 202);
    const auto *slf_206 = buffer.data(slf + 206);
    const auto *slf_208 = buffer.data(slf + 208);
    const auto *slf_209 = buffer.data(slf + 209);
    const auto *slf_210 = buffer.data(slf + 210);
    const auto *slf_212 = buffer.data(slf + 212);
    const auto *slf_216 = buffer.data(slf + 216);
    const auto *slf_218 = buffer.data(slf + 218);
    const auto *slf_219 = buffer.data(slf + 219);
    const auto *slf_220 = buffer.data(slf + 220);
    const auto *slf_222 = buffer.data(slf + 222);
    const auto *slf_226 = buffer.data(slf + 226);
    const auto *slf_228 = buffer.data(slf + 228);
    const auto *slf_229 = buffer.data(slf + 229);
    const auto *slf_230 = buffer.data(slf + 230);
    const auto *slf_232 = buffer.data(slf + 232);
    const auto *slf_236 = buffer.data(slf + 236);
    const auto *slf_238 = buffer.data(slf + 238);
    const auto *slf_239 = buffer.data(slf + 239);
    const auto *slf_240 = buffer.data(slf + 240);
    const auto *slf_242 = buffer.data(slf + 242);
    const auto *slf_243 = buffer.data(slf + 243);
    const auto *slf_245 = buffer.data(slf + 245);
    const auto *slf_246 = buffer.data(slf + 246);
    const auto *slf_247 = buffer.data(slf + 247);
    const auto *slf_248 = buffer.data(slf + 248);
    const auto *slf_249 = buffer.data(slf + 249);
    const auto *slf_250 = buffer.data(slf + 250);
    const auto *slf_253 = buffer.data(slf + 253);
    const auto *slf_255 = buffer.data(slf + 255);
    const auto *slf_256 = buffer.data(slf + 256);
    const auto *slf_257 = buffer.data(slf + 257);
    const auto *slf_258 = buffer.data(slf + 258);
    const auto *slf_259 = buffer.data(slf + 259);
    const auto *slf_266 = buffer.data(slf + 266);
    const auto *slf_267 = buffer.data(slf + 267);
    const auto *slf_268 = buffer.data(slf + 268);
    const auto *slf_269 = buffer.data(slf + 269);
    const auto *slf_270 = buffer.data(slf + 270);
    const auto *slf_273 = buffer.data(slf + 273);
    const auto *slf_275 = buffer.data(slf + 275);
    const auto *slf_276 = buffer.data(slf + 276);
    const auto *slf_277 = buffer.data(slf + 277);
    const auto *slf_278 = buffer.data(slf + 278);
    const auto *slf_279 = buffer.data(slf + 279);
    const auto *slf_280 = buffer.data(slf + 280);
    const auto *slf_283 = buffer.data(slf + 283);
    const auto *slf_285 = buffer.data(slf + 285);
    const auto *slf_286 = buffer.data(slf + 286);
    const auto *slf_287 = buffer.data(slf + 287);
    const auto *slf_288 = buffer.data(slf + 288);
    const auto *slf_289 = buffer.data(slf + 289);
    const auto *slf_295 = buffer.data(slf + 295);
    const auto *slf_296 = buffer.data(slf + 296);
    const auto *slf_297 = buffer.data(slf + 297);
    const auto *slf_298 = buffer.data(slf + 298);
    const auto *slf_299 = buffer.data(slf + 299);
    const auto *slf_300 = buffer.data(slf + 300);
    const auto *slf_303 = buffer.data(slf + 303);
    const auto *slf_305 = buffer.data(slf + 305);
    const auto *slf_306 = buffer.data(slf + 306);
    const auto *slf_307 = buffer.data(slf + 307);
    const auto *slf_308 = buffer.data(slf + 308);
    const auto *slf_309 = buffer.data(slf + 309);
    const auto *slf_310 = buffer.data(slf + 310);
    const auto *slf_313 = buffer.data(slf + 313);
    const auto *slf_315 = buffer.data(slf + 315);
    const auto *slf_316 = buffer.data(slf + 316);
    const auto *slf_317 = buffer.data(slf + 317);
    const auto *slf_318 = buffer.data(slf + 318);
    const auto *slf_319 = buffer.data(slf + 319);

    const auto *slg1_300 = buffer.data(slg1 + 300);
    const auto *slg1_303 = buffer.data(slg1 + 303);
    const auto *slg1_305 = buffer.data(slg1 + 305);
    const auto *slg1_314 = buffer.data(slg1 + 314);
    const auto *slg1_315 = buffer.data(slg1 + 315);
    const auto *slg1_318 = buffer.data(slg1 + 318);
    const auto *slg1_325 = buffer.data(slg1 + 325);

    const auto *smd0_147 = buffer.data(smd0 + 147);
    const auto *smd0_149 = buffer.data(smd0 + 149);
    const auto *smd0_150 = buffer.data(smd0 + 150);
    const auto *smd0_153 = buffer.data(smd0 + 153);
    const auto *smd0_155 = buffer.data(smd0 + 155);
    const auto *smd0_159 = buffer.data(smd0 + 159);
    const auto *smd0_161 = buffer.data(smd0 + 161);
    const auto *smd0_162 = buffer.data(smd0 + 162);
    const auto *smd0_165 = buffer.data(smd0 + 165);
    const auto *smd0_167 = buffer.data(smd0 + 167);
    const auto *smd0_168 = buffer.data(smd0 + 168);
    const auto *smd0_171 = buffer.data(smd0 + 171);
    const auto *smd0_173 = buffer.data(smd0 + 173);
    const auto *smd0_179 = buffer.data(smd0 + 179);
    const auto *smd0_180 = buffer.data(smd0 + 180);
    const auto *smd0_183 = buffer.data(smd0 + 183);
    const auto *smd0_185 = buffer.data(smd0 + 185);
    const auto *smd0_186 = buffer.data(smd0 + 186);
    const auto *smd0_189 = buffer.data(smd0 + 189);
    const auto *smd0_191 = buffer.data(smd0 + 191);

    const auto *smd1_147 = buffer.data(smd1 + 147);
    const auto *smd1_149 = buffer.data(smd1 + 149);
    const auto *smd1_150 = buffer.data(smd1 + 150);
    const auto *smd1_153 = buffer.data(smd1 + 153);
    const auto *smd1_155 = buffer.data(smd1 + 155);
    const auto *smd1_159 = buffer.data(smd1 + 159);
    const auto *smd1_161 = buffer.data(smd1 + 161);
    const auto *smd1_162 = buffer.data(smd1 + 162);
    const auto *smd1_165 = buffer.data(smd1 + 165);
    const auto *smd1_167 = buffer.data(smd1 + 167);
    const auto *smd1_168 = buffer.data(smd1 + 168);
    const auto *smd1_171 = buffer.data(smd1 + 171);
    const auto *smd1_173 = buffer.data(smd1 + 173);
    const auto *smd1_179 = buffer.data(smd1 + 179);
    const auto *smd1_180 = buffer.data(smd1 + 180);
    const auto *smd1_183 = buffer.data(smd1 + 183);
    const auto *smd1_185 = buffer.data(smd1 + 185);
    const auto *smd1_186 = buffer.data(smd1 + 186);
    const auto *smd1_189 = buffer.data(smd1 + 189);
    const auto *smd1_191 = buffer.data(smd1 + 191);

    const auto *smf_240 = buffer.data(smf + 240);
    const auto *smf_242 = buffer.data(smf + 242);
    const auto *smf_243 = buffer.data(smf + 243);
    const auto *smf_245 = buffer.data(smf + 245);
    const auto *smf_246 = buffer.data(smf + 246);
    const auto *smf_247 = buffer.data(smf + 247);
    const auto *smf_248 = buffer.data(smf + 248);
    const auto *smf_249 = buffer.data(smf + 249);
    const auto *smf_250 = buffer.data(smf + 250);
    const auto *smf_252 = buffer.data(smf + 252);
    const auto *smf_253 = buffer.data(smf + 253);
    const auto *smf_255 = buffer.data(smf + 255);
    const auto *smf_256 = buffer.data(smf + 256);
    const auto *smf_257 = buffer.data(smf + 257);
    const auto *smf_258 = buffer.data(smf + 258);
    const auto *smf_259 = buffer.data(smf + 259);
    const auto *smf_260 = buffer.data(smf + 260);
    const auto *smf_262 = buffer.data(smf + 262);
    const auto *smf_266 = buffer.data(smf + 266);
    const auto *smf_267 = buffer.data(smf + 267);
    const auto *smf_268 = buffer.data(smf + 268);
    const auto *smf_269 = buffer.data(smf + 269);
    const auto *smf_270 = buffer.data(smf + 270);
    const auto *smf_272 = buffer.data(smf + 272);
    const auto *smf_273 = buffer.data(smf + 273);
    const auto *smf_275 = buffer.data(smf + 275);
    const auto *smf_276 = buffer.data(smf + 276);
    const auto *smf_277 = buffer.data(smf + 277);
    const auto *smf_278 = buffer.data(smf + 278);
    const auto *smf_279 = buffer.data(smf + 279);
    const auto *smf_280 = buffer.data(smf + 280);
    const auto *smf_282 = buffer.data(smf + 282);
    const auto *smf_283 = buffer.data(smf + 283);
    const auto *smf_285 = buffer.data(smf + 285);
    const auto *smf_286 = buffer.data(smf + 286);
    const auto *smf_287 = buffer.data(smf + 287);
    const auto *smf_288 = buffer.data(smf + 288);
    const auto *smf_289 = buffer.data(smf + 289);
    const auto *smf_290 = buffer.data(smf + 290);
    const auto *smf_292 = buffer.data(smf + 292);
    const auto *smf_295 = buffer.data(smf + 295);
    const auto *smf_296 = buffer.data(smf + 296);
    const auto *smf_297 = buffer.data(smf + 297);
    const auto *smf_298 = buffer.data(smf + 298);
    const auto *smf_299 = buffer.data(smf + 299);
    const auto *smf_300 = buffer.data(smf + 300);
    const auto *smf_302 = buffer.data(smf + 302);
    const auto *smf_303 = buffer.data(smf + 303);
    const auto *smf_305 = buffer.data(smf + 305);
    const auto *smf_306 = buffer.data(smf + 306);
    const auto *smf_307 = buffer.data(smf + 307);
    const auto *smf_308 = buffer.data(smf + 308);
    const auto *smf_309 = buffer.data(smf + 309);
    const auto *smf_310 = buffer.data(smf + 310);
    const auto *smf_312 = buffer.data(smf + 312);
    const auto *smf_313 = buffer.data(smf + 313);
    const auto *smf_315 = buffer.data(smf + 315);
    const auto *smf_316 = buffer.data(smf + 316);
    const auto *smf_317 = buffer.data(smf + 317);
    const auto *smf_318 = buffer.data(smf + 318);
    const auto *smf_319 = buffer.data(smf + 319);

#pragma omp simd aligned(t_361, t_362, t_363, t_364, pc_x, pc_y, pc_z, slf_170, slf_180, \
                         slf_182, slf_243, smd0_147, smd1_147, smf_240, smf_242, \
                         smf_243 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_361[k] = f_12 * slf_180[k]
                   + f_3 * pc_y[k] * smf_240[k];

        t_362[k] = f_12 * slf_170[k]
                   + f_3 * pc_z[k] * smf_240[k];

        t_363[k] = f_12 * slf_243[k]
                   + f_4 * smd0_147[k]
                   - f_5 * smd1_147[k]
                   + f_3 * pc_x[k] * smf_243[k];

        t_364[k] = f_12 * slf_182[k]
                   + f_3 * pc_y[k] * smf_242[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, pc_x, slf_245, slf_246, slf_247, slf_248, \
                         smd0_149, smd1_149, smf_245, smf_246, smf_247, \
                         smf_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = f_12 * slf_245[k]
                   + f_4 * smd0_149[k]
                   - f_5 * smd1_149[k]
                   + f_3 * pc_x[k] * smf_245[k];

        t_366[k] = f_12 * slf_246[k]
                   + f_3 * pc_x[k] * smf_246[k];

        t_367[k] = f_12 * slf_247[k]
                   + f_3 * pc_x[k] * smf_247[k];

        t_368[k] = f_12 * slf_248[k]
                   + f_3 * pc_x[k] * smf_248[k];
    }

#pragma omp simd aligned(t_369, t_370, t_371, pc_x, pc_y, pc_z, slf_176, slf_186, slf_249, \
                         smd0_147, smd1_147, smf_246, smf_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_369[k] = f_12 * slf_249[k]
                   + f_3 * pc_x[k] * smf_249[k];

        t_370[k] = f_12 * slf_186[k]
                   + f_1 * smd0_147[k]
                   - f_2 * smd1_147[k]
                   + f_3 * pc_y[k] * smf_246[k];

        t_371[k] = f_12 * slf_176[k]
                   + f_3 * pc_z[k] * smf_246[k];
    }

#pragma omp simd aligned(t_372, t_373, t_374, pc_y, pc_z, slf_179, slf_188, slf_189, smd0_149, \
                         smd1_149, smf_248, smf_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_372[k] = f_12 * slf_188[k]
                   + f_4 * smd0_149[k]
                   - f_5 * smd1_149[k]
                   + f_3 * pc_y[k] * smf_248[k];

        t_373[k] = f_12 * slf_189[k]
                   + f_3 * pc_y[k] * smf_249[k];

        t_374[k] = f_12 * slf_179[k]
                   + f_1 * smd0_149[k]
                   - f_2 * smd1_149[k]
                   + f_3 * pc_z[k] * smf_249[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, pc_x, pc_y, pc_z, slf_180, slf_190, slf_250, \
                         smd0_150, smd1_150, smf_250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = f_12 * slf_250[k]
                   + f_1 * smd0_150[k]
                   - f_2 * smd1_150[k]
                   + f_3 * pc_x[k] * smf_250[k];

        t_376[k] = f_8 * slf_190[k]
                   + f_3 * pc_y[k] * smf_250[k];

        t_377[k] = f_14 * slf_180[k]
                   + f_3 * pc_z[k] * smf_250[k];
    }

#pragma omp simd aligned(t_378, t_379, t_380, pc_x, pc_y, slf_192, slf_253, slf_255, smd0_153, \
                         smd0_155, smd1_153, smd1_155, smf_252, smf_253, \
                         smf_255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_378[k] = f_12 * slf_253[k]
                   + f_4 * smd0_153[k]
                   - f_5 * smd1_153[k]
                   + f_3 * pc_x[k] * smf_253[k];

        t_379[k] = f_8 * slf_192[k]
                   + f_3 * pc_y[k] * smf_252[k];

        t_380[k] = f_12 * slf_255[k]
                   + f_4 * smd0_155[k]
                   - f_5 * smd1_155[k]
                   + f_3 * pc_x[k] * smf_255[k];
    }

#pragma omp simd aligned(t_381, t_382, t_383, t_384, pc_x, slf_256, slf_257, slf_258, slf_259, \
                         smf_256, smf_257, smf_258, smf_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_381[k] = f_12 * slf_256[k]
                   + f_3 * pc_x[k] * smf_256[k];

        t_382[k] = f_12 * slf_257[k]
                   + f_3 * pc_x[k] * smf_257[k];

        t_383[k] = f_12 * slf_258[k]
                   + f_3 * pc_x[k] * smf_258[k];

        t_384[k] = f_12 * slf_259[k]
                   + f_3 * pc_x[k] * smf_259[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, pc_y, pc_z, slf_186, slf_196, slf_198, smd0_153, \
                         smd0_155, smd1_153, smd1_155, smf_256, \
                         smf_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = f_8 * slf_196[k]
                   + f_1 * smd0_153[k]
                   - f_2 * smd1_153[k]
                   + f_3 * pc_y[k] * smf_256[k];

        t_386[k] = f_14 * slf_186[k]
                   + f_3 * pc_z[k] * smf_256[k];

        t_387[k] = f_8 * slf_198[k]
                   + f_4 * smd0_155[k]
                   - f_5 * smd1_155[k]
                   + f_3 * pc_y[k] * smf_258[k];
    }

#pragma omp simd aligned(t_388, t_389, t_390, t_391, pb_y, pc_y, pc_z, slg0_300, slf_189, \
                         slf_199, slf_200, slg1_300, smd0_155, smd1_155, smf_259, \
                         smf_260 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_388[k] = f_8 * slf_199[k]
                   + f_3 * pc_y[k] * smf_259[k];

        t_389[k] = f_14 * slf_189[k]
                   + f_1 * smd0_155[k]
                   - f_2 * smd1_155[k]
                   + f_3 * pc_z[k] * smf_259[k];

        t_390[k] = pb_y[k] * slg0_300[k]
                   - f_6 * pc_y[k] * slg1_300[k];

        t_391[k] = f_7 * slf_200[k]
                   + f_3 * pc_y[k] * smf_260[k];
    }

#pragma omp simd aligned(t_392, t_393, t_394, t_395, pb_y, pc_y, pc_z, slg0_303, slg0_305, \
                         slf_190, slf_201, slf_202, slg1_303, slg1_305, smf_260, \
                         smf_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_392[k] = f_13 * slf_190[k]
                   + f_3 * pc_z[k] * smf_260[k];

        t_393[k] = pb_y[k] * slg0_303[k]
                   + f_8 * slf_201[k]
                   - f_6 * pc_y[k] * slg1_303[k];

        t_394[k] = f_7 * slf_202[k]
                   + f_3 * pc_y[k] * smf_262[k];

        t_395[k] = pb_y[k] * slg0_305[k]
                   - f_6 * pc_y[k] * slg1_305[k];
    }

#pragma omp simd aligned(t_396, t_397, t_398, t_399, pc_x, slf_266, slf_267, slf_268, slf_269, \
                         smf_266, smf_267, smf_268, smf_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_396[k] = f_12 * slf_266[k]
                   + f_3 * pc_x[k] * smf_266[k];

        t_397[k] = f_12 * slf_267[k]
                   + f_3 * pc_x[k] * smf_267[k];

        t_398[k] = f_12 * slf_268[k]
                   + f_3 * pc_x[k] * smf_268[k];

        t_399[k] = f_12 * slf_269[k]
                   + f_3 * pc_x[k] * smf_269[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, pc_y, pc_z, slf_196, slf_206, slf_208, smd0_159, \
                         smd0_161, smd1_159, smd1_161, smf_266, \
                         smf_268 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = f_7 * slf_206[k]
                   + f_1 * smd0_159[k]
                   - f_2 * smd1_159[k]
                   + f_3 * pc_y[k] * smf_266[k];

        t_401[k] = f_13 * slf_196[k]
                   + f_3 * pc_z[k] * smf_266[k];

        t_402[k] = f_7 * slf_208[k]
                   + f_4 * smd0_161[k]
                   - f_5 * smd1_161[k]
                   + f_3 * pc_y[k] * smf_268[k];
    }

#pragma omp simd aligned(t_403, t_404, t_405, t_406, pb_y, pc_x, pc_y, slg0_314, slf_209, \
                         slf_270, slg1_314, smd0_162, smd1_162, smf_269, \
                         smf_270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_403[k] = f_7 * slf_209[k]
                   + f_3 * pc_y[k] * smf_269[k];

        t_404[k] = pb_y[k] * slg0_314[k]
                   - f_6 * pc_y[k] * slg1_314[k];

        t_405[k] = f_12 * slf_270[k]
                   + f_1 * smd0_162[k]
                   - f_2 * smd1_162[k]
                   + f_3 * pc_x[k] * smf_270[k];

        t_406[k] = f_3 * pc_y[k] * smf_270[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, pc_x, pc_y, pc_z, slf_200, slf_273, smd0_165, \
                         smd1_165, smf_270, smf_272, smf_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = f_11 * slf_200[k]
                   + f_3 * pc_z[k] * smf_270[k];

        t_408[k] = f_12 * slf_273[k]
                   + f_4 * smd0_165[k]
                   - f_5 * smd1_165[k]
                   + f_3 * pc_x[k] * smf_273[k];

        t_409[k] = f_3 * pc_y[k] * smf_272[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, pc_x, slf_275, slf_276, slf_277, slf_278, \
                         smd0_167, smd1_167, smf_275, smf_276, smf_277, \
                         smf_278 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = f_12 * slf_275[k]
                   + f_4 * smd0_167[k]
                   - f_5 * smd1_167[k]
                   + f_3 * pc_x[k] * smf_275[k];

        t_411[k] = f_12 * slf_276[k]
                   + f_3 * pc_x[k] * smf_276[k];

        t_412[k] = f_12 * slf_277[k]
                   + f_3 * pc_x[k] * smf_277[k];

        t_413[k] = f_12 * slf_278[k]
                   + f_3 * pc_x[k] * smf_278[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, t_417, pc_x, pc_y, pc_z, slf_206, slf_279, \
                         smd0_165, smd0_167, smd1_165, smd1_167, smf_276, smf_278, \
                         smf_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_12 * slf_279[k]
                   + f_3 * pc_x[k] * smf_279[k];

        t_415[k] = f_1 * smd0_165[k]
                   - f_2 * smd1_165[k]
                   + f_3 * pc_y[k] * smf_276[k];

        t_416[k] = f_11 * slf_206[k]
                   + f_3 * pc_z[k] * smf_276[k];

        t_417[k] = f_4 * smd0_167[k]
                   - f_5 * smd1_167[k]
                   + f_3 * pc_y[k] * smf_278[k];
    }

#pragma omp simd aligned(t_418, t_419, t_420, t_421, pc_x, pc_y, pc_z, slf_209, slf_210, \
                         slf_280, smd0_167, smd0_168, smd1_167, smd1_168, smf_279, \
                         smf_280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_418[k] = f_3 * pc_y[k] * smf_279[k];

        t_419[k] = f_11 * slf_209[k]
                   + f_1 * smd0_167[k]
                   - f_2 * smd1_167[k]
                   + f_3 * pc_z[k] * smf_279[k];

        t_420[k] = f_8 * slf_280[k]
                   + f_1 * smd0_168[k]
                   - f_2 * smd1_168[k]
                   + f_3 * pc_x[k] * smf_280[k];

        t_421[k] = f_10 * slf_210[k]
                   + f_3 * pc_y[k] * smf_280[k];
    }

#pragma omp simd aligned(t_422, t_423, t_424, pc_x, pc_y, pc_z, slf_212, slf_283, smd0_171, \
                         smd1_171, smf_280, smf_282, smf_283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_422[k] = f_3 * pc_z[k] * smf_280[k];

        t_423[k] = f_8 * slf_283[k]
                   + f_4 * smd0_171[k]
                   - f_5 * smd1_171[k]
                   + f_3 * pc_x[k] * smf_283[k];

        t_424[k] = f_10 * slf_212[k]
                   + f_3 * pc_y[k] * smf_282[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, pc_x, slf_285, slf_286, slf_287, slf_288, \
                         smd0_173, smd1_173, smf_285, smf_286, smf_287, \
                         smf_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = f_8 * slf_285[k]
                   + f_4 * smd0_173[k]
                   - f_5 * smd1_173[k]
                   + f_3 * pc_x[k] * smf_285[k];

        t_426[k] = f_8 * slf_286[k]
                   + f_3 * pc_x[k] * smf_286[k];

        t_427[k] = f_8 * slf_287[k]
                   + f_3 * pc_x[k] * smf_287[k];

        t_428[k] = f_8 * slf_288[k]
                   + f_3 * pc_x[k] * smf_288[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, pc_x, pc_y, pc_z, slf_216, slf_289, smd0_171, \
                         smd1_171, smf_286, smf_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = f_8 * slf_289[k]
                   + f_3 * pc_x[k] * smf_289[k];

        t_430[k] = f_10 * slf_216[k]
                   + f_1 * smd0_171[k]
                   - f_2 * smd1_171[k]
                   + f_3 * pc_y[k] * smf_286[k];

        t_431[k] = f_3 * pc_z[k] * smf_286[k];
    }

#pragma omp simd aligned(t_432, t_433, t_434, t_435, pb_z, pc_y, pc_z, slg0_315, slf_218, \
                         slf_219, slg1_315, smd0_173, smd1_173, smf_288, \
                         smf_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_432[k] = f_10 * slf_218[k]
                   + f_4 * smd0_173[k]
                   - f_5 * smd1_173[k]
                   + f_3 * pc_y[k] * smf_288[k];

        t_433[k] = f_10 * slf_219[k]
                   + f_3 * pc_y[k] * smf_289[k];

        t_434[k] = f_1 * smd0_173[k]
                   - f_2 * smd1_173[k]
                   + f_3 * pc_z[k] * smf_289[k];

        t_435[k] = pb_z[k] * slg0_315[k]
                   - f_6 * pc_z[k] * slg1_315[k];
    }

#pragma omp simd aligned(t_436, t_437, t_438, t_439, pb_z, pc_y, pc_z, slg0_318, slf_210, \
                         slf_220, slf_222, slg1_318, smf_290, smf_292 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_436[k] = f_11 * slf_220[k]
                   + f_3 * pc_y[k] * smf_290[k];

        t_437[k] = f_7 * slf_210[k]
                   + f_3 * pc_z[k] * smf_290[k];

        t_438[k] = pb_z[k] * slg0_318[k]
                   - f_6 * pc_z[k] * slg1_318[k];

        t_439[k] = f_11 * slf_222[k]
                   + f_3 * pc_y[k] * smf_292[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, pc_x, slf_295, slf_296, slf_297, slf_298, \
                         smd0_179, smd1_179, smf_295, smf_296, smf_297, \
                         smf_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = f_8 * slf_295[k]
                   + f_4 * smd0_179[k]
                   - f_5 * smd1_179[k]
                   + f_3 * pc_x[k] * smf_295[k];

        t_441[k] = f_8 * slf_296[k]
                   + f_3 * pc_x[k] * smf_296[k];

        t_442[k] = f_8 * slf_297[k]
                   + f_3 * pc_x[k] * smf_297[k];

        t_443[k] = f_8 * slf_298[k]
                   + f_3 * pc_x[k] * smf_298[k];
    }

#pragma omp simd aligned(t_444, t_445, t_446, pb_z, pc_x, pc_z, slg0_325, slf_216, slf_299, \
                         slg1_325, smf_296, smf_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_444[k] = f_8 * slf_299[k]
                   + f_3 * pc_x[k] * smf_299[k];

        t_445[k] = pb_z[k] * slg0_325[k]
                   - f_6 * pc_z[k] * slg1_325[k];

        t_446[k] = f_7 * slf_216[k]
                   + f_3 * pc_z[k] * smf_296[k];
    }

#pragma omp simd aligned(t_447, t_448, t_449, pc_y, pc_z, slf_219, slf_228, slf_229, smd0_179, \
                         smd1_179, smf_298, smf_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_447[k] = f_11 * slf_228[k]
                   + f_4 * smd0_179[k]
                   - f_5 * smd1_179[k]
                   + f_3 * pc_y[k] * smf_298[k];

        t_448[k] = f_11 * slf_229[k]
                   + f_3 * pc_y[k] * smf_299[k];

        t_449[k] = f_7 * slf_219[k]
                   + f_1 * smd0_179[k]
                   - f_2 * smd1_179[k]
                   + f_3 * pc_z[k] * smf_299[k];
    }

#pragma omp simd aligned(t_450, t_451, t_452, pc_x, pc_y, pc_z, slf_220, slf_230, slf_300, \
                         smd0_180, smd1_180, smf_300 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_450[k] = f_8 * slf_300[k]
                   + f_1 * smd0_180[k]
                   - f_2 * smd1_180[k]
                   + f_3 * pc_x[k] * smf_300[k];

        t_451[k] = f_13 * slf_230[k]
                   + f_3 * pc_y[k] * smf_300[k];

        t_452[k] = f_8 * slf_220[k]
                   + f_3 * pc_z[k] * smf_300[k];
    }

#pragma omp simd aligned(t_453, t_454, t_455, pc_x, pc_y, slf_232, slf_303, slf_305, smd0_183, \
                         smd0_185, smd1_183, smd1_185, smf_302, smf_303, \
                         smf_305 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_453[k] = f_8 * slf_303[k]
                   + f_4 * smd0_183[k]
                   - f_5 * smd1_183[k]
                   + f_3 * pc_x[k] * smf_303[k];

        t_454[k] = f_13 * slf_232[k]
                   + f_3 * pc_y[k] * smf_302[k];

        t_455[k] = f_8 * slf_305[k]
                   + f_4 * smd0_185[k]
                   - f_5 * smd1_185[k]
                   + f_3 * pc_x[k] * smf_305[k];
    }

#pragma omp simd aligned(t_456, t_457, t_458, t_459, pc_x, slf_306, slf_307, slf_308, slf_309, \
                         smf_306, smf_307, smf_308, smf_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_456[k] = f_8 * slf_306[k]
                   + f_3 * pc_x[k] * smf_306[k];

        t_457[k] = f_8 * slf_307[k]
                   + f_3 * pc_x[k] * smf_307[k];

        t_458[k] = f_8 * slf_308[k]
                   + f_3 * pc_x[k] * smf_308[k];

        t_459[k] = f_8 * slf_309[k]
                   + f_3 * pc_x[k] * smf_309[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, pc_y, pc_z, slf_226, slf_236, slf_238, smd0_183, \
                         smd0_185, smd1_183, smd1_185, smf_306, \
                         smf_308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = f_13 * slf_236[k]
                   + f_1 * smd0_183[k]
                   - f_2 * smd1_183[k]
                   + f_3 * pc_y[k] * smf_306[k];

        t_461[k] = f_8 * slf_226[k]
                   + f_3 * pc_z[k] * smf_306[k];

        t_462[k] = f_13 * slf_238[k]
                   + f_4 * smd0_185[k]
                   - f_5 * smd1_185[k]
                   + f_3 * pc_y[k] * smf_308[k];
    }

#pragma omp simd aligned(t_463, t_464, t_465, pc_x, pc_y, pc_z, slf_229, slf_239, slf_310, \
                         smd0_185, smd0_186, smd1_185, smd1_186, smf_309, \
                         smf_310 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_463[k] = f_13 * slf_239[k]
                   + f_3 * pc_y[k] * smf_309[k];

        t_464[k] = f_8 * slf_229[k]
                   + f_1 * smd0_185[k]
                   - f_2 * smd1_185[k]
                   + f_3 * pc_z[k] * smf_309[k];

        t_465[k] = f_8 * slf_310[k]
                   + f_1 * smd0_186[k]
                   - f_2 * smd1_186[k]
                   + f_3 * pc_x[k] * smf_310[k];
    }

#pragma omp simd aligned(t_466, t_467, t_468, t_469, pc_x, pc_y, pc_z, slf_230, slf_240, \
                         slf_242, slf_313, smd0_189, smd1_189, smf_310, smf_312, \
                         smf_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_466[k] = f_14 * slf_240[k]
                   + f_3 * pc_y[k] * smf_310[k];

        t_467[k] = f_12 * slf_230[k]
                   + f_3 * pc_z[k] * smf_310[k];

        t_468[k] = f_8 * slf_313[k]
                   + f_4 * smd0_189[k]
                   - f_5 * smd1_189[k]
                   + f_3 * pc_x[k] * smf_313[k];

        t_469[k] = f_14 * slf_242[k]
                   + f_3 * pc_y[k] * smf_312[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, pc_x, slf_315, slf_316, slf_317, slf_318, \
                         smd0_191, smd1_191, smf_315, smf_316, smf_317, \
                         smf_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = f_8 * slf_315[k]
                   + f_4 * smd0_191[k]
                   - f_5 * smd1_191[k]
                   + f_3 * pc_x[k] * smf_315[k];

        t_471[k] = f_8 * slf_316[k]
                   + f_3 * pc_x[k] * smf_316[k];

        t_472[k] = f_8 * slf_317[k]
                   + f_3 * pc_x[k] * smf_317[k];

        t_473[k] = f_8 * slf_318[k]
                   + f_3 * pc_x[k] * smf_318[k];
    }

#pragma omp simd aligned(t_474, t_475, t_476, pc_x, pc_y, pc_z, slf_236, slf_246, slf_319, \
                         smd0_189, smd1_189, smf_316, smf_319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = f_8 * slf_319[k]
                   + f_3 * pc_x[k] * smf_319[k];

        t_475[k] = f_14 * slf_246[k]
                   + f_1 * smd0_189[k]
                   - f_2 * smd1_189[k]
                   + f_3 * pc_y[k] * smf_316[k];

        t_476[k] = f_12 * slf_236[k]
                   + f_3 * pc_z[k] * smf_316[k];
    }
}

static auto
compute_prim_smg_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t slg0,
                                                          const size_t slf, const size_t slg1,
                                                          const size_t smd0, const size_t smd1,
                                                          const size_t smf, const size_t ncols,
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
    const auto f_9 = 4.0 / q;
    const auto f_10 = 3.5 / q;
    const auto f_11 = 3.0 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 2.5 / q;
    const auto f_14 = 2.0 / q;

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
    auto *t_594 = buffer.data(target + 594);
    auto *t_595 = buffer.data(target + 595);
    auto *t_596 = buffer.data(target + 596);
    auto *t_597 = buffer.data(target + 597);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *slg0_405 = buffer.data(slg0 + 405);
    const auto *slg0_408 = buffer.data(slg0 + 408);
    const auto *slg0_410 = buffer.data(slg0 + 410);
    const auto *slg0_419 = buffer.data(slg0 + 419);
    const auto *slg0_420 = buffer.data(slg0 + 420);
    const auto *slg0_423 = buffer.data(slg0 + 423);
    const auto *slg0_540 = buffer.data(slg0 + 540);
    const auto *slg0_543 = buffer.data(slg0 + 543);
    const auto *slg0_545 = buffer.data(slg0 + 545);
    const auto *slg0_550 = buffer.data(slg0 + 550);
    const auto *slg0_552 = buffer.data(slg0 + 552);
    const auto *slg0_554 = buffer.data(slg0 + 554);
    const auto *slg0_560 = buffer.data(slg0 + 560);
    const auto *slg0_565 = buffer.data(slg0 + 565);
    const auto *slg0_567 = buffer.data(slg0 + 567);
    const auto *slg0_569 = buffer.data(slg0 + 569);
    const auto *slg0_570 = buffer.data(slg0 + 570);
    const auto *slg0_573 = buffer.data(slg0 + 573);
    const auto *slg0_575 = buffer.data(slg0 + 575);
    const auto *slg0_580 = buffer.data(slg0 + 580);
    const auto *slg0_582 = buffer.data(slg0 + 582);
    const auto *slg0_584 = buffer.data(slg0 + 584);
    const auto *slg0_585 = buffer.data(slg0 + 585);
    const auto *slg0_588 = buffer.data(slg0 + 588);
    const auto *slg0_590 = buffer.data(slg0 + 590);
    const auto *slg0_595 = buffer.data(slg0 + 595);
    const auto *slg0_597 = buffer.data(slg0 + 597);

    const auto *slf_239 = buffer.data(slf + 239);
    const auto *slf_240 = buffer.data(slf + 240);
    const auto *slf_246 = buffer.data(slf + 246);
    const auto *slf_248 = buffer.data(slf + 248);
    const auto *slf_249 = buffer.data(slf + 249);
    const auto *slf_250 = buffer.data(slf + 250);
    const auto *slf_252 = buffer.data(slf + 252);
    const auto *slf_256 = buffer.data(slf + 256);
    const auto *slf_258 = buffer.data(slf + 258);
    const auto *slf_259 = buffer.data(slf + 259);
    const auto *slf_260 = buffer.data(slf + 260);
    const auto *slf_262 = buffer.data(slf + 262);
    const auto *slf_266 = buffer.data(slf + 266);
    const auto *slf_268 = buffer.data(slf + 268);
    const auto *slf_269 = buffer.data(slf + 269);
    const auto *slf_270 = buffer.data(slf + 270);
    const auto *slf_271 = buffer.data(slf + 271);
    const auto *slf_272 = buffer.data(slf + 272);
    const auto *slf_276 = buffer.data(slf + 276);
    const auto *slf_278 = buffer.data(slf + 278);
    const auto *slf_279 = buffer.data(slf + 279);
    const auto *slf_280 = buffer.data(slf + 280);
    const auto *slf_282 = buffer.data(slf + 282);
    const auto *slf_286 = buffer.data(slf + 286);
    const auto *slf_289 = buffer.data(slf + 289);
    const auto *slf_290 = buffer.data(slf + 290);
    const auto *slf_292 = buffer.data(slf + 292);
    const auto *slf_296 = buffer.data(slf + 296);
    const auto *slf_299 = buffer.data(slf + 299);
    const auto *slf_300 = buffer.data(slf + 300);
    const auto *slf_302 = buffer.data(slf + 302);
    const auto *slf_306 = buffer.data(slf + 306);
    const auto *slf_309 = buffer.data(slf + 309);
    const auto *slf_310 = buffer.data(slf + 310);
    const auto *slf_312 = buffer.data(slf + 312);
    const auto *slf_320 = buffer.data(slf + 320);
    const auto *slf_323 = buffer.data(slf + 323);
    const auto *slf_325 = buffer.data(slf + 325);
    const auto *slf_326 = buffer.data(slf + 326);
    const auto *slf_327 = buffer.data(slf + 327);
    const auto *slf_328 = buffer.data(slf + 328);
    const auto *slf_329 = buffer.data(slf + 329);
    const auto *slf_330 = buffer.data(slf + 330);
    const auto *slf_333 = buffer.data(slf + 333);
    const auto *slf_335 = buffer.data(slf + 335);
    const auto *slf_336 = buffer.data(slf + 336);
    const auto *slf_337 = buffer.data(slf + 337);
    const auto *slf_338 = buffer.data(slf + 338);
    const auto *slf_339 = buffer.data(slf + 339);
    const auto *slf_346 = buffer.data(slf + 346);
    const auto *slf_347 = buffer.data(slf + 347);
    const auto *slf_348 = buffer.data(slf + 348);
    const auto *slf_349 = buffer.data(slf + 349);
    const auto *slf_350 = buffer.data(slf + 350);
    const auto *slf_353 = buffer.data(slf + 353);
    const auto *slf_355 = buffer.data(slf + 355);
    const auto *slf_356 = buffer.data(slf + 356);
    const auto *slf_357 = buffer.data(slf + 357);
    const auto *slf_358 = buffer.data(slf + 358);
    const auto *slf_359 = buffer.data(slf + 359);
    const auto *slf_360 = buffer.data(slf + 360);
    const auto *slf_363 = buffer.data(slf + 363);
    const auto *slf_365 = buffer.data(slf + 365);
    const auto *slf_366 = buffer.data(slf + 366);
    const auto *slf_367 = buffer.data(slf + 367);
    const auto *slf_368 = buffer.data(slf + 368);
    const auto *slf_369 = buffer.data(slf + 369);
    const auto *slf_375 = buffer.data(slf + 375);
    const auto *slf_376 = buffer.data(slf + 376);
    const auto *slf_377 = buffer.data(slf + 377);
    const auto *slf_378 = buffer.data(slf + 378);
    const auto *slf_379 = buffer.data(slf + 379);
    const auto *slf_380 = buffer.data(slf + 380);
    const auto *slf_383 = buffer.data(slf + 383);
    const auto *slf_385 = buffer.data(slf + 385);
    const auto *slf_386 = buffer.data(slf + 386);
    const auto *slf_387 = buffer.data(slf + 387);
    const auto *slf_388 = buffer.data(slf + 388);
    const auto *slf_389 = buffer.data(slf + 389);
    const auto *slf_390 = buffer.data(slf + 390);
    const auto *slf_393 = buffer.data(slf + 393);
    const auto *slf_395 = buffer.data(slf + 395);
    const auto *slf_396 = buffer.data(slf + 396);
    const auto *slf_397 = buffer.data(slf + 397);
    const auto *slf_398 = buffer.data(slf + 398);
    const auto *slf_399 = buffer.data(slf + 399);

    const auto *slg1_405 = buffer.data(slg1 + 405);
    const auto *slg1_408 = buffer.data(slg1 + 408);
    const auto *slg1_410 = buffer.data(slg1 + 410);
    const auto *slg1_419 = buffer.data(slg1 + 419);
    const auto *slg1_420 = buffer.data(slg1 + 420);
    const auto *slg1_423 = buffer.data(slg1 + 423);
    const auto *slg1_540 = buffer.data(slg1 + 540);
    const auto *slg1_543 = buffer.data(slg1 + 543);
    const auto *slg1_545 = buffer.data(slg1 + 545);
    const auto *slg1_550 = buffer.data(slg1 + 550);
    const auto *slg1_552 = buffer.data(slg1 + 552);
    const auto *slg1_554 = buffer.data(slg1 + 554);
    const auto *slg1_560 = buffer.data(slg1 + 560);
    const auto *slg1_565 = buffer.data(slg1 + 565);
    const auto *slg1_567 = buffer.data(slg1 + 567);
    const auto *slg1_569 = buffer.data(slg1 + 569);
    const auto *slg1_570 = buffer.data(slg1 + 570);
    const auto *slg1_573 = buffer.data(slg1 + 573);
    const auto *slg1_575 = buffer.data(slg1 + 575);
    const auto *slg1_580 = buffer.data(slg1 + 580);
    const auto *slg1_582 = buffer.data(slg1 + 582);
    const auto *slg1_584 = buffer.data(slg1 + 584);
    const auto *slg1_585 = buffer.data(slg1 + 585);
    const auto *slg1_588 = buffer.data(slg1 + 588);
    const auto *slg1_590 = buffer.data(slg1 + 590);
    const auto *slg1_595 = buffer.data(slg1 + 595);
    const auto *slg1_597 = buffer.data(slg1 + 597);

    const auto *smd0_191 = buffer.data(smd0 + 191);
    const auto *smd0_192 = buffer.data(smd0 + 192);
    const auto *smd0_195 = buffer.data(smd0 + 195);
    const auto *smd0_197 = buffer.data(smd0 + 197);
    const auto *smd0_198 = buffer.data(smd0 + 198);
    const auto *smd0_201 = buffer.data(smd0 + 201);
    const auto *smd0_203 = buffer.data(smd0 + 203);
    const auto *smd0_207 = buffer.data(smd0 + 207);
    const auto *smd0_209 = buffer.data(smd0 + 209);
    const auto *smd0_210 = buffer.data(smd0 + 210);
    const auto *smd0_213 = buffer.data(smd0 + 213);
    const auto *smd0_215 = buffer.data(smd0 + 215);

    const auto *smd1_191 = buffer.data(smd1 + 191);
    const auto *smd1_192 = buffer.data(smd1 + 192);
    const auto *smd1_195 = buffer.data(smd1 + 195);
    const auto *smd1_197 = buffer.data(smd1 + 197);
    const auto *smd1_198 = buffer.data(smd1 + 198);
    const auto *smd1_201 = buffer.data(smd1 + 201);
    const auto *smd1_203 = buffer.data(smd1 + 203);
    const auto *smd1_207 = buffer.data(smd1 + 207);
    const auto *smd1_209 = buffer.data(smd1 + 209);
    const auto *smd1_210 = buffer.data(smd1 + 210);
    const auto *smd1_213 = buffer.data(smd1 + 213);
    const auto *smd1_215 = buffer.data(smd1 + 215);

    const auto *smf_318 = buffer.data(smf + 318);
    const auto *smf_319 = buffer.data(smf + 319);
    const auto *smf_320 = buffer.data(smf + 320);
    const auto *smf_322 = buffer.data(smf + 322);
    const auto *smf_323 = buffer.data(smf + 323);
    const auto *smf_325 = buffer.data(smf + 325);
    const auto *smf_326 = buffer.data(smf + 326);
    const auto *smf_327 = buffer.data(smf + 327);
    const auto *smf_328 = buffer.data(smf + 328);
    const auto *smf_329 = buffer.data(smf + 329);
    const auto *smf_330 = buffer.data(smf + 330);
    const auto *smf_332 = buffer.data(smf + 332);
    const auto *smf_333 = buffer.data(smf + 333);
    const auto *smf_335 = buffer.data(smf + 335);
    const auto *smf_336 = buffer.data(smf + 336);
    const auto *smf_337 = buffer.data(smf + 337);
    const auto *smf_338 = buffer.data(smf + 338);
    const auto *smf_339 = buffer.data(smf + 339);
    const auto *smf_340 = buffer.data(smf + 340);
    const auto *smf_342 = buffer.data(smf + 342);
    const auto *smf_346 = buffer.data(smf + 346);
    const auto *smf_347 = buffer.data(smf + 347);
    const auto *smf_348 = buffer.data(smf + 348);
    const auto *smf_349 = buffer.data(smf + 349);
    const auto *smf_350 = buffer.data(smf + 350);
    const auto *smf_352 = buffer.data(smf + 352);
    const auto *smf_353 = buffer.data(smf + 353);
    const auto *smf_355 = buffer.data(smf + 355);
    const auto *smf_356 = buffer.data(smf + 356);
    const auto *smf_357 = buffer.data(smf + 357);
    const auto *smf_358 = buffer.data(smf + 358);
    const auto *smf_359 = buffer.data(smf + 359);
    const auto *smf_360 = buffer.data(smf + 360);
    const auto *smf_362 = buffer.data(smf + 362);
    const auto *smf_366 = buffer.data(smf + 366);
    const auto *smf_367 = buffer.data(smf + 367);
    const auto *smf_368 = buffer.data(smf + 368);
    const auto *smf_369 = buffer.data(smf + 369);
    const auto *smf_370 = buffer.data(smf + 370);
    const auto *smf_372 = buffer.data(smf + 372);
    const auto *smf_376 = buffer.data(smf + 376);
    const auto *smf_377 = buffer.data(smf + 377);
    const auto *smf_378 = buffer.data(smf + 378);
    const auto *smf_379 = buffer.data(smf + 379);
    const auto *smf_380 = buffer.data(smf + 380);
    const auto *smf_382 = buffer.data(smf + 382);
    const auto *smf_386 = buffer.data(smf + 386);
    const auto *smf_387 = buffer.data(smf + 387);
    const auto *smf_388 = buffer.data(smf + 388);
    const auto *smf_389 = buffer.data(smf + 389);
    const auto *smf_390 = buffer.data(smf + 390);
    const auto *smf_392 = buffer.data(smf + 392);
    const auto *smf_396 = buffer.data(smf + 396);
    const auto *smf_397 = buffer.data(smf + 397);
    const auto *smf_398 = buffer.data(smf + 398);
    const auto *smf_399 = buffer.data(smf + 399);

#pragma omp simd aligned(t_477, t_478, t_479, pc_y, pc_z, slf_239, slf_248, slf_249, smd0_191, \
                         smd1_191, smf_318, smf_319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_477[k] = f_14 * slf_248[k]
                   + f_4 * smd0_191[k]
                   - f_5 * smd1_191[k]
                   + f_3 * pc_y[k] * smf_318[k];

        t_478[k] = f_14 * slf_249[k]
                   + f_3 * pc_y[k] * smf_319[k];

        t_479[k] = f_12 * slf_239[k]
                   + f_1 * smd0_191[k]
                   - f_2 * smd1_191[k]
                   + f_3 * pc_z[k] * smf_319[k];
    }

#pragma omp simd aligned(t_480, t_481, t_482, pc_x, pc_y, pc_z, slf_240, slf_250, slf_320, \
                         smd0_192, smd1_192, smf_320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_480[k] = f_8 * slf_320[k]
                   + f_1 * smd0_192[k]
                   - f_2 * smd1_192[k]
                   + f_3 * pc_x[k] * smf_320[k];

        t_481[k] = f_12 * slf_250[k]
                   + f_3 * pc_y[k] * smf_320[k];

        t_482[k] = f_14 * slf_240[k]
                   + f_3 * pc_z[k] * smf_320[k];
    }

#pragma omp simd aligned(t_483, t_484, t_485, pc_x, pc_y, slf_252, slf_323, slf_325, smd0_195, \
                         smd0_197, smd1_195, smd1_197, smf_322, smf_323, \
                         smf_325 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_483[k] = f_8 * slf_323[k]
                   + f_4 * smd0_195[k]
                   - f_5 * smd1_195[k]
                   + f_3 * pc_x[k] * smf_323[k];

        t_484[k] = f_12 * slf_252[k]
                   + f_3 * pc_y[k] * smf_322[k];

        t_485[k] = f_8 * slf_325[k]
                   + f_4 * smd0_197[k]
                   - f_5 * smd1_197[k]
                   + f_3 * pc_x[k] * smf_325[k];
    }

#pragma omp simd aligned(t_486, t_487, t_488, t_489, pc_x, slf_326, slf_327, slf_328, slf_329, \
                         smf_326, smf_327, smf_328, smf_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_486[k] = f_8 * slf_326[k]
                   + f_3 * pc_x[k] * smf_326[k];

        t_487[k] = f_8 * slf_327[k]
                   + f_3 * pc_x[k] * smf_327[k];

        t_488[k] = f_8 * slf_328[k]
                   + f_3 * pc_x[k] * smf_328[k];

        t_489[k] = f_8 * slf_329[k]
                   + f_3 * pc_x[k] * smf_329[k];
    }

#pragma omp simd aligned(t_490, t_491, t_492, pc_y, pc_z, slf_246, slf_256, slf_258, smd0_195, \
                         smd0_197, smd1_195, smd1_197, smf_326, \
                         smf_328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_490[k] = f_12 * slf_256[k]
                   + f_1 * smd0_195[k]
                   - f_2 * smd1_195[k]
                   + f_3 * pc_y[k] * smf_326[k];

        t_491[k] = f_14 * slf_246[k]
                   + f_3 * pc_z[k] * smf_326[k];

        t_492[k] = f_12 * slf_258[k]
                   + f_4 * smd0_197[k]
                   - f_5 * smd1_197[k]
                   + f_3 * pc_y[k] * smf_328[k];
    }

#pragma omp simd aligned(t_493, t_494, t_495, pc_x, pc_y, pc_z, slf_249, slf_259, slf_330, \
                         smd0_197, smd0_198, smd1_197, smd1_198, smf_329, \
                         smf_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_493[k] = f_12 * slf_259[k]
                   + f_3 * pc_y[k] * smf_329[k];

        t_494[k] = f_14 * slf_249[k]
                   + f_1 * smd0_197[k]
                   - f_2 * smd1_197[k]
                   + f_3 * pc_z[k] * smf_329[k];

        t_495[k] = f_8 * slf_330[k]
                   + f_1 * smd0_198[k]
                   - f_2 * smd1_198[k]
                   + f_3 * pc_x[k] * smf_330[k];
    }

#pragma omp simd aligned(t_496, t_497, t_498, t_499, pc_x, pc_y, pc_z, slf_250, slf_260, \
                         slf_262, slf_333, smd0_201, smd1_201, smf_330, smf_332, \
                         smf_333 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_496[k] = f_8 * slf_260[k]
                   + f_3 * pc_y[k] * smf_330[k];

        t_497[k] = f_13 * slf_250[k]
                   + f_3 * pc_z[k] * smf_330[k];

        t_498[k] = f_8 * slf_333[k]
                   + f_4 * smd0_201[k]
                   - f_5 * smd1_201[k]
                   + f_3 * pc_x[k] * smf_333[k];

        t_499[k] = f_8 * slf_262[k]
                   + f_3 * pc_y[k] * smf_332[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, pc_x, slf_335, slf_336, slf_337, slf_338, \
                         smd0_203, smd1_203, smf_335, smf_336, smf_337, \
                         smf_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = f_8 * slf_335[k]
                   + f_4 * smd0_203[k]
                   - f_5 * smd1_203[k]
                   + f_3 * pc_x[k] * smf_335[k];

        t_501[k] = f_8 * slf_336[k]
                   + f_3 * pc_x[k] * smf_336[k];

        t_502[k] = f_8 * slf_337[k]
                   + f_3 * pc_x[k] * smf_337[k];

        t_503[k] = f_8 * slf_338[k]
                   + f_3 * pc_x[k] * smf_338[k];
    }

#pragma omp simd aligned(t_504, t_505, t_506, pc_x, pc_y, pc_z, slf_256, slf_266, slf_339, \
                         smd0_201, smd1_201, smf_336, smf_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_504[k] = f_8 * slf_339[k]
                   + f_3 * pc_x[k] * smf_339[k];

        t_505[k] = f_8 * slf_266[k]
                   + f_1 * smd0_201[k]
                   - f_2 * smd1_201[k]
                   + f_3 * pc_y[k] * smf_336[k];

        t_506[k] = f_13 * slf_256[k]
                   + f_3 * pc_z[k] * smf_336[k];
    }

#pragma omp simd aligned(t_507, t_508, t_509, t_510, pb_y, pc_y, pc_z, slg0_405, slf_259, \
                         slf_268, slf_269, slg1_405, smd0_203, smd1_203, smf_338, \
                         smf_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_507[k] = f_8 * slf_268[k]
                   + f_4 * smd0_203[k]
                   - f_5 * smd1_203[k]
                   + f_3 * pc_y[k] * smf_338[k];

        t_508[k] = f_8 * slf_269[k]
                   + f_3 * pc_y[k] * smf_339[k];

        t_509[k] = f_13 * slf_259[k]
                   + f_1 * smd0_203[k]
                   - f_2 * smd1_203[k]
                   + f_3 * pc_z[k] * smf_339[k];

        t_510[k] = pb_y[k] * slg0_405[k]
                   - f_6 * pc_y[k] * slg1_405[k];
    }

#pragma omp simd aligned(t_511, t_512, t_513, t_514, pb_y, pc_y, pc_z, slg0_408, slf_260, \
                         slf_270, slf_271, slf_272, slg1_408, smf_340, \
                         smf_342 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_511[k] = f_7 * slf_270[k]
                   + f_3 * pc_y[k] * smf_340[k];

        t_512[k] = f_11 * slf_260[k]
                   + f_3 * pc_z[k] * smf_340[k];

        t_513[k] = pb_y[k] * slg0_408[k]
                   + f_8 * slf_271[k]
                   - f_6 * pc_y[k] * slg1_408[k];

        t_514[k] = f_7 * slf_272[k]
                   + f_3 * pc_y[k] * smf_342[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, t_518, pb_y, pc_x, pc_y, slg0_410, slf_346, \
                         slf_347, slf_348, slg1_410, smf_346, smf_347, \
                         smf_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = pb_y[k] * slg0_410[k]
                   - f_6 * pc_y[k] * slg1_410[k];

        t_516[k] = f_8 * slf_346[k]
                   + f_3 * pc_x[k] * smf_346[k];

        t_517[k] = f_8 * slf_347[k]
                   + f_3 * pc_x[k] * smf_347[k];

        t_518[k] = f_8 * slf_348[k]
                   + f_3 * pc_x[k] * smf_348[k];
    }

#pragma omp simd aligned(t_519, t_520, t_521, pc_x, pc_y, pc_z, slf_266, slf_276, slf_349, \
                         smd0_207, smd1_207, smf_346, smf_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_519[k] = f_8 * slf_349[k]
                   + f_3 * pc_x[k] * smf_349[k];

        t_520[k] = f_7 * slf_276[k]
                   + f_1 * smd0_207[k]
                   - f_2 * smd1_207[k]
                   + f_3 * pc_y[k] * smf_346[k];

        t_521[k] = f_11 * slf_266[k]
                   + f_3 * pc_z[k] * smf_346[k];
    }

#pragma omp simd aligned(t_522, t_523, t_524, pb_y, pc_y, slg0_419, slf_278, slf_279, \
                         slg1_419, smd0_209, smd1_209, smf_348, \
                         smf_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_522[k] = f_7 * slf_278[k]
                   + f_4 * smd0_209[k]
                   - f_5 * smd1_209[k]
                   + f_3 * pc_y[k] * smf_348[k];

        t_523[k] = f_7 * slf_279[k]
                   + f_3 * pc_y[k] * smf_349[k];

        t_524[k] = pb_y[k] * slg0_419[k]
                   - f_6 * pc_y[k] * slg1_419[k];
    }

#pragma omp simd aligned(t_525, t_526, t_527, t_528, pc_x, pc_y, pc_z, slf_270, slf_350, \
                         slf_353, smd0_210, smd0_213, smd1_210, smd1_213, smf_350, \
                         smf_353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_525[k] = f_8 * slf_350[k]
                   + f_1 * smd0_210[k]
                   - f_2 * smd1_210[k]
                   + f_3 * pc_x[k] * smf_350[k];

        t_526[k] = f_3 * pc_y[k] * smf_350[k];

        t_527[k] = f_10 * slf_270[k]
                   + f_3 * pc_z[k] * smf_350[k];

        t_528[k] = f_8 * slf_353[k]
                   + f_4 * smd0_213[k]
                   - f_5 * smd1_213[k]
                   + f_3 * pc_x[k] * smf_353[k];
    }

#pragma omp simd aligned(t_529, t_530, t_531, t_532, pc_x, pc_y, slf_355, slf_356, slf_357, \
                         smd0_215, smd1_215, smf_352, smf_355, smf_356, \
                         smf_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_529[k] = f_3 * pc_y[k] * smf_352[k];

        t_530[k] = f_8 * slf_355[k]
                   + f_4 * smd0_215[k]
                   - f_5 * smd1_215[k]
                   + f_3 * pc_x[k] * smf_355[k];

        t_531[k] = f_8 * slf_356[k]
                   + f_3 * pc_x[k] * smf_356[k];

        t_532[k] = f_8 * slf_357[k]
                   + f_3 * pc_x[k] * smf_357[k];
    }

#pragma omp simd aligned(t_533, t_534, t_535, t_536, pc_x, pc_y, pc_z, slf_276, slf_358, \
                         slf_359, smd0_213, smd1_213, smf_356, smf_358, \
                         smf_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_533[k] = f_8 * slf_358[k]
                   + f_3 * pc_x[k] * smf_358[k];

        t_534[k] = f_8 * slf_359[k]
                   + f_3 * pc_x[k] * smf_359[k];

        t_535[k] = f_1 * smd0_213[k]
                   - f_2 * smd1_213[k]
                   + f_3 * pc_y[k] * smf_356[k];

        t_536[k] = f_10 * slf_276[k]
                   + f_3 * pc_z[k] * smf_356[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, t_540, pb_x, pc_x, pc_y, pc_z, slg0_540, \
                         slf_279, slf_360, slg1_540, smd0_215, smd1_215, smf_358, \
                         smf_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = f_4 * smd0_215[k]
                   - f_5 * smd1_215[k]
                   + f_3 * pc_y[k] * smf_358[k];

        t_538[k] = f_3 * pc_y[k] * smf_359[k];

        t_539[k] = f_10 * slf_279[k]
                   + f_1 * smd0_215[k]
                   - f_2 * smd1_215[k]
                   + f_3 * pc_z[k] * smf_359[k];

        t_540[k] = pb_x[k] * slg0_540[k]
                   + f_14 * slf_360[k]
                   - f_6 * pc_x[k] * slg1_540[k];
    }

#pragma omp simd aligned(t_541, t_542, t_543, t_544, pb_x, pc_x, pc_y, pc_z, slg0_543, \
                         slf_280, slf_282, slf_363, slg1_543, smf_360, \
                         smf_362 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_541[k] = f_9 * slf_280[k]
                   + f_3 * pc_y[k] * smf_360[k];

        t_542[k] = f_3 * pc_z[k] * smf_360[k];

        t_543[k] = pb_x[k] * slg0_543[k]
                   + f_8 * slf_363[k]
                   - f_6 * pc_x[k] * slg1_543[k];

        t_544[k] = f_9 * slf_282[k]
                   + f_3 * pc_y[k] * smf_362[k];
    }

#pragma omp simd aligned(t_545, t_546, t_547, t_548, pb_x, pc_x, slg0_545, slf_365, slf_366, \
                         slf_367, slf_368, slg1_545, smf_366, smf_367, \
                         smf_368 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_545[k] = pb_x[k] * slg0_545[k]
                   + f_8 * slf_365[k]
                   - f_6 * pc_x[k] * slg1_545[k];

        t_546[k] = f_7 * slf_366[k]
                   + f_3 * pc_x[k] * smf_366[k];

        t_547[k] = f_7 * slf_367[k]
                   + f_3 * pc_x[k] * smf_367[k];

        t_548[k] = f_7 * slf_368[k]
                   + f_3 * pc_x[k] * smf_368[k];
    }

#pragma omp simd aligned(t_549, t_550, t_551, t_552, pb_x, pc_x, pc_z, slg0_550, slg0_552, \
                         slf_369, slg1_550, slg1_552, smf_366, \
                         smf_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_549[k] = f_7 * slf_369[k]
                   + f_3 * pc_x[k] * smf_369[k];

        t_550[k] = pb_x[k] * slg0_550[k]
                   - f_6 * pc_x[k] * slg1_550[k];

        t_551[k] = f_3 * pc_z[k] * smf_366[k];

        t_552[k] = pb_x[k] * slg0_552[k]
                   - f_6 * pc_x[k] * slg1_552[k];
    }

#pragma omp simd aligned(t_553, t_554, t_555, pb_x, pb_z, pc_x, pc_y, pc_z, slg0_420, \
                         slg0_554, slf_289, slg1_420, slg1_554, \
                         smf_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_553[k] = f_9 * slf_289[k]
                   + f_3 * pc_y[k] * smf_369[k];

        t_554[k] = pb_x[k] * slg0_554[k]
                   - f_6 * pc_x[k] * slg1_554[k];

        t_555[k] = pb_z[k] * slg0_420[k]
                   - f_6 * pc_z[k] * slg1_420[k];
    }

#pragma omp simd aligned(t_556, t_557, t_558, t_559, pb_z, pc_y, pc_z, slg0_423, slf_280, \
                         slf_290, slf_292, slg1_423, smf_370, smf_372 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_556[k] = f_10 * slf_290[k]
                   + f_3 * pc_y[k] * smf_370[k];

        t_557[k] = f_7 * slf_280[k]
                   + f_3 * pc_z[k] * smf_370[k];

        t_558[k] = pb_z[k] * slg0_423[k]
                   - f_6 * pc_z[k] * slg1_423[k];

        t_559[k] = f_10 * slf_292[k]
                   + f_3 * pc_y[k] * smf_372[k];
    }

#pragma omp simd aligned(t_560, t_561, t_562, t_563, pb_x, pc_x, slg0_560, slf_375, slf_376, \
                         slf_377, slf_378, slg1_560, smf_376, smf_377, \
                         smf_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_560[k] = pb_x[k] * slg0_560[k]
                   + f_8 * slf_375[k]
                   - f_6 * pc_x[k] * slg1_560[k];

        t_561[k] = f_7 * slf_376[k]
                   + f_3 * pc_x[k] * smf_376[k];

        t_562[k] = f_7 * slf_377[k]
                   + f_3 * pc_x[k] * smf_377[k];

        t_563[k] = f_7 * slf_378[k]
                   + f_3 * pc_x[k] * smf_378[k];
    }

#pragma omp simd aligned(t_564, t_565, t_566, t_567, pb_x, pc_x, pc_z, slg0_565, slg0_567, \
                         slf_286, slf_379, slg1_565, slg1_567, smf_376, \
                         smf_379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_564[k] = f_7 * slf_379[k]
                   + f_3 * pc_x[k] * smf_379[k];

        t_565[k] = pb_x[k] * slg0_565[k]
                   - f_6 * pc_x[k] * slg1_565[k];

        t_566[k] = f_7 * slf_286[k]
                   + f_3 * pc_z[k] * smf_376[k];

        t_567[k] = pb_x[k] * slg0_567[k]
                   - f_6 * pc_x[k] * slg1_567[k];
    }

#pragma omp simd aligned(t_568, t_569, t_570, t_571, pb_x, pc_x, pc_y, slg0_569, slg0_570, \
                         slf_299, slf_300, slf_380, slg1_569, slg1_570, smf_379, \
                         smf_380 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_568[k] = f_10 * slf_299[k]
                   + f_3 * pc_y[k] * smf_379[k];

        t_569[k] = pb_x[k] * slg0_569[k]
                   - f_6 * pc_x[k] * slg1_569[k];

        t_570[k] = pb_x[k] * slg0_570[k]
                   + f_14 * slf_380[k]
                   - f_6 * pc_x[k] * slg1_570[k];

        t_571[k] = f_11 * slf_300[k]
                   + f_3 * pc_y[k] * smf_380[k];
    }

#pragma omp simd aligned(t_572, t_573, t_574, pb_x, pc_x, pc_y, pc_z, slg0_573, slf_290, \
                         slf_302, slf_383, slg1_573, smf_380, smf_382 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_572[k] = f_8 * slf_290[k]
                   + f_3 * pc_z[k] * smf_380[k];

        t_573[k] = pb_x[k] * slg0_573[k]
                   + f_8 * slf_383[k]
                   - f_6 * pc_x[k] * slg1_573[k];

        t_574[k] = f_11 * slf_302[k]
                   + f_3 * pc_y[k] * smf_382[k];
    }

#pragma omp simd aligned(t_575, t_576, t_577, t_578, pb_x, pc_x, slg0_575, slf_385, slf_386, \
                         slf_387, slf_388, slg1_575, smf_386, smf_387, \
                         smf_388 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_575[k] = pb_x[k] * slg0_575[k]
                   + f_8 * slf_385[k]
                   - f_6 * pc_x[k] * slg1_575[k];

        t_576[k] = f_7 * slf_386[k]
                   + f_3 * pc_x[k] * smf_386[k];

        t_577[k] = f_7 * slf_387[k]
                   + f_3 * pc_x[k] * smf_387[k];

        t_578[k] = f_7 * slf_388[k]
                   + f_3 * pc_x[k] * smf_388[k];
    }

#pragma omp simd aligned(t_579, t_580, t_581, t_582, pb_x, pc_x, pc_z, slg0_580, slg0_582, \
                         slf_296, slf_389, slg1_580, slg1_582, smf_386, \
                         smf_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_579[k] = f_7 * slf_389[k]
                   + f_3 * pc_x[k] * smf_389[k];

        t_580[k] = pb_x[k] * slg0_580[k]
                   - f_6 * pc_x[k] * slg1_580[k];

        t_581[k] = f_8 * slf_296[k]
                   + f_3 * pc_z[k] * smf_386[k];

        t_582[k] = pb_x[k] * slg0_582[k]
                   - f_6 * pc_x[k] * slg1_582[k];
    }

#pragma omp simd aligned(t_583, t_584, t_585, t_586, pb_x, pc_x, pc_y, slg0_584, slg0_585, \
                         slf_309, slf_310, slf_390, slg1_584, slg1_585, smf_389, \
                         smf_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_583[k] = f_11 * slf_309[k]
                   + f_3 * pc_y[k] * smf_389[k];

        t_584[k] = pb_x[k] * slg0_584[k]
                   - f_6 * pc_x[k] * slg1_584[k];

        t_585[k] = pb_x[k] * slg0_585[k]
                   + f_14 * slf_390[k]
                   - f_6 * pc_x[k] * slg1_585[k];

        t_586[k] = f_13 * slf_310[k]
                   + f_3 * pc_y[k] * smf_390[k];
    }

#pragma omp simd aligned(t_587, t_588, t_589, pb_x, pc_x, pc_y, pc_z, slg0_588, slf_300, \
                         slf_312, slf_393, slg1_588, smf_390, smf_392 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_587[k] = f_12 * slf_300[k]
                   + f_3 * pc_z[k] * smf_390[k];

        t_588[k] = pb_x[k] * slg0_588[k]
                   + f_8 * slf_393[k]
                   - f_6 * pc_x[k] * slg1_588[k];

        t_589[k] = f_13 * slf_312[k]
                   + f_3 * pc_y[k] * smf_392[k];
    }

#pragma omp simd aligned(t_590, t_591, t_592, t_593, pb_x, pc_x, slg0_590, slf_395, slf_396, \
                         slf_397, slf_398, slg1_590, smf_396, smf_397, \
                         smf_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_590[k] = pb_x[k] * slg0_590[k]
                   + f_8 * slf_395[k]
                   - f_6 * pc_x[k] * slg1_590[k];

        t_591[k] = f_7 * slf_396[k]
                   + f_3 * pc_x[k] * smf_396[k];

        t_592[k] = f_7 * slf_397[k]
                   + f_3 * pc_x[k] * smf_397[k];

        t_593[k] = f_7 * slf_398[k]
                   + f_3 * pc_x[k] * smf_398[k];
    }

#pragma omp simd aligned(t_594, t_595, t_596, t_597, pb_x, pc_x, pc_z, slg0_595, slg0_597, \
                         slf_306, slf_399, slg1_595, slg1_597, smf_396, \
                         smf_399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_594[k] = f_7 * slf_399[k]
                   + f_3 * pc_x[k] * smf_399[k];

        t_595[k] = pb_x[k] * slg0_595[k]
                   - f_6 * pc_x[k] * slg1_595[k];

        t_596[k] = f_12 * slf_306[k]
                   + f_3 * pc_z[k] * smf_396[k];

        t_597[k] = pb_x[k] * slg0_597[k]
                   - f_6 * pc_x[k] * slg1_597[k];
    }
}

static auto
compute_prim_smg_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t slg0,
                                                          const size_t slf, const size_t slg1,
                                                          const size_t smd0, const size_t smd1,
                                                          const size_t smf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 4.0 / q;
    const auto f_10 = 3.5 / q;
    const auto f_11 = 3.0 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 2.5 / q;
    const auto f_14 = 2.0 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *slg0_525 = buffer.data(slg0 + 525);
    const auto *slg0_530 = buffer.data(slg0 + 530);
    const auto *slg0_540 = buffer.data(slg0 + 540);
    const auto *slg0_543 = buffer.data(slg0 + 543);
    const auto *slg0_550 = buffer.data(slg0 + 550);
    const auto *slg0_552 = buffer.data(slg0 + 552);
    const auto *slg0_599 = buffer.data(slg0 + 599);
    const auto *slg0_600 = buffer.data(slg0 + 600);
    const auto *slg0_603 = buffer.data(slg0 + 603);
    const auto *slg0_605 = buffer.data(slg0 + 605);
    const auto *slg0_610 = buffer.data(slg0 + 610);
    const auto *slg0_612 = buffer.data(slg0 + 612);
    const auto *slg0_614 = buffer.data(slg0 + 614);
    const auto *slg0_615 = buffer.data(slg0 + 615);
    const auto *slg0_618 = buffer.data(slg0 + 618);
    const auto *slg0_620 = buffer.data(slg0 + 620);
    const auto *slg0_625 = buffer.data(slg0 + 625);
    const auto *slg0_627 = buffer.data(slg0 + 627);
    const auto *slg0_629 = buffer.data(slg0 + 629);
    const auto *slg0_630 = buffer.data(slg0 + 630);
    const auto *slg0_633 = buffer.data(slg0 + 633);
    const auto *slg0_635 = buffer.data(slg0 + 635);
    const auto *slg0_640 = buffer.data(slg0 + 640);
    const auto *slg0_642 = buffer.data(slg0 + 642);
    const auto *slg0_644 = buffer.data(slg0 + 644);
    const auto *slg0_648 = buffer.data(slg0 + 648);
    const auto *slg0_655 = buffer.data(slg0 + 655);
    const auto *slg0_657 = buffer.data(slg0 + 657);
    const auto *slg0_659 = buffer.data(slg0 + 659);
    const auto *slg0_660 = buffer.data(slg0 + 660);
    const auto *slg0_663 = buffer.data(slg0 + 663);
    const auto *slg0_665 = buffer.data(slg0 + 665);
    const auto *slg0_670 = buffer.data(slg0 + 670);
    const auto *slg0_672 = buffer.data(slg0 + 672);
    const auto *slg0_674 = buffer.data(slg0 + 674);

    const auto *slf_310 = buffer.data(slf + 310);
    const auto *slf_316 = buffer.data(slf + 316);
    const auto *slf_319 = buffer.data(slf + 319);
    const auto *slf_320 = buffer.data(slf + 320);
    const auto *slf_322 = buffer.data(slf + 322);
    const auto *slf_326 = buffer.data(slf + 326);
    const auto *slf_329 = buffer.data(slf + 329);
    const auto *slf_330 = buffer.data(slf + 330);
    const auto *slf_332 = buffer.data(slf + 332);
    const auto *slf_336 = buffer.data(slf + 336);
    const auto *slf_339 = buffer.data(slf + 339);
    const auto *slf_340 = buffer.data(slf + 340);
    const auto *slf_342 = buffer.data(slf + 342);
    const auto *slf_346 = buffer.data(slf + 346);
    const auto *slf_349 = buffer.data(slf + 349);
    const auto *slf_350 = buffer.data(slf + 350);
    const auto *slf_352 = buffer.data(slf + 352);
    const auto *slf_356 = buffer.data(slf + 356);
    const auto *slf_359 = buffer.data(slf + 359);
    const auto *slf_360 = buffer.data(slf + 360);
    const auto *slf_362 = buffer.data(slf + 362);
    const auto *slf_366 = buffer.data(slf + 366);
    const auto *slf_367 = buffer.data(slf + 367);
    const auto *slf_368 = buffer.data(slf + 368);
    const auto *slf_369 = buffer.data(slf + 369);
    const auto *slf_370 = buffer.data(slf + 370);
    const auto *slf_372 = buffer.data(slf + 372);
    const auto *slf_376 = buffer.data(slf + 376);
    const auto *slf_379 = buffer.data(slf + 379);
    const auto *slf_380 = buffer.data(slf + 380);
    const auto *slf_382 = buffer.data(slf + 382);
    const auto *slf_386 = buffer.data(slf + 386);
    const auto *slf_388 = buffer.data(slf + 388);
    const auto *slf_389 = buffer.data(slf + 389);
    const auto *slf_390 = buffer.data(slf + 390);
    const auto *slf_392 = buffer.data(slf + 392);
    const auto *slf_400 = buffer.data(slf + 400);
    const auto *slf_403 = buffer.data(slf + 403);
    const auto *slf_405 = buffer.data(slf + 405);
    const auto *slf_406 = buffer.data(slf + 406);
    const auto *slf_407 = buffer.data(slf + 407);
    const auto *slf_408 = buffer.data(slf + 408);
    const auto *slf_409 = buffer.data(slf + 409);
    const auto *slf_410 = buffer.data(slf + 410);
    const auto *slf_413 = buffer.data(slf + 413);
    const auto *slf_415 = buffer.data(slf + 415);
    const auto *slf_416 = buffer.data(slf + 416);
    const auto *slf_417 = buffer.data(slf + 417);
    const auto *slf_418 = buffer.data(slf + 418);
    const auto *slf_419 = buffer.data(slf + 419);
    const auto *slf_420 = buffer.data(slf + 420);
    const auto *slf_423 = buffer.data(slf + 423);
    const auto *slf_425 = buffer.data(slf + 425);
    const auto *slf_426 = buffer.data(slf + 426);
    const auto *slf_427 = buffer.data(slf + 427);
    const auto *slf_428 = buffer.data(slf + 428);
    const auto *slf_429 = buffer.data(slf + 429);
    const auto *slf_433 = buffer.data(slf + 433);
    const auto *slf_436 = buffer.data(slf + 436);
    const auto *slf_437 = buffer.data(slf + 437);
    const auto *slf_438 = buffer.data(slf + 438);
    const auto *slf_439 = buffer.data(slf + 439);
    const auto *slf_440 = buffer.data(slf + 440);
    const auto *slf_443 = buffer.data(slf + 443);
    const auto *slf_445 = buffer.data(slf + 445);
    const auto *slf_446 = buffer.data(slf + 446);
    const auto *slf_447 = buffer.data(slf + 447);
    const auto *slf_448 = buffer.data(slf + 448);
    const auto *slf_449 = buffer.data(slf + 449);

    const auto *slg1_525 = buffer.data(slg1 + 525);
    const auto *slg1_530 = buffer.data(slg1 + 530);
    const auto *slg1_540 = buffer.data(slg1 + 540);
    const auto *slg1_543 = buffer.data(slg1 + 543);
    const auto *slg1_550 = buffer.data(slg1 + 550);
    const auto *slg1_552 = buffer.data(slg1 + 552);
    const auto *slg1_599 = buffer.data(slg1 + 599);
    const auto *slg1_600 = buffer.data(slg1 + 600);
    const auto *slg1_603 = buffer.data(slg1 + 603);
    const auto *slg1_605 = buffer.data(slg1 + 605);
    const auto *slg1_610 = buffer.data(slg1 + 610);
    const auto *slg1_612 = buffer.data(slg1 + 612);
    const auto *slg1_614 = buffer.data(slg1 + 614);
    const auto *slg1_615 = buffer.data(slg1 + 615);
    const auto *slg1_618 = buffer.data(slg1 + 618);
    const auto *slg1_620 = buffer.data(slg1 + 620);
    const auto *slg1_625 = buffer.data(slg1 + 625);
    const auto *slg1_627 = buffer.data(slg1 + 627);
    const auto *slg1_629 = buffer.data(slg1 + 629);
    const auto *slg1_630 = buffer.data(slg1 + 630);
    const auto *slg1_633 = buffer.data(slg1 + 633);
    const auto *slg1_635 = buffer.data(slg1 + 635);
    const auto *slg1_640 = buffer.data(slg1 + 640);
    const auto *slg1_642 = buffer.data(slg1 + 642);
    const auto *slg1_644 = buffer.data(slg1 + 644);
    const auto *slg1_648 = buffer.data(slg1 + 648);
    const auto *slg1_655 = buffer.data(slg1 + 655);
    const auto *slg1_657 = buffer.data(slg1 + 657);
    const auto *slg1_659 = buffer.data(slg1 + 659);
    const auto *slg1_660 = buffer.data(slg1 + 660);
    const auto *slg1_663 = buffer.data(slg1 + 663);
    const auto *slg1_665 = buffer.data(slg1 + 665);
    const auto *slg1_670 = buffer.data(slg1 + 670);
    const auto *slg1_672 = buffer.data(slg1 + 672);
    const auto *slg1_674 = buffer.data(slg1 + 674);

    const auto *smd0_270 = buffer.data(smd0 + 270);
    const auto *smd0_273 = buffer.data(smd0 + 273);
    const auto *smd0_275 = buffer.data(smd0 + 275);
    const auto *smd0_281 = buffer.data(smd0 + 281);
    const auto *smd0_282 = buffer.data(smd0 + 282);
    const auto *smd0_285 = buffer.data(smd0 + 285);
    const auto *smd0_287 = buffer.data(smd0 + 287);
    const auto *smd0_288 = buffer.data(smd0 + 288);
    const auto *smd0_291 = buffer.data(smd0 + 291);
    const auto *smd0_293 = buffer.data(smd0 + 293);

    const auto *smd1_270 = buffer.data(smd1 + 270);
    const auto *smd1_273 = buffer.data(smd1 + 273);
    const auto *smd1_275 = buffer.data(smd1 + 275);
    const auto *smd1_281 = buffer.data(smd1 + 281);
    const auto *smd1_282 = buffer.data(smd1 + 282);
    const auto *smd1_285 = buffer.data(smd1 + 285);
    const auto *smd1_287 = buffer.data(smd1 + 287);
    const auto *smd1_288 = buffer.data(smd1 + 288);
    const auto *smd1_291 = buffer.data(smd1 + 291);
    const auto *smd1_293 = buffer.data(smd1 + 293);

    const auto *smf_399 = buffer.data(smf + 399);
    const auto *smf_400 = buffer.data(smf + 400);
    const auto *smf_402 = buffer.data(smf + 402);
    const auto *smf_406 = buffer.data(smf + 406);
    const auto *smf_407 = buffer.data(smf + 407);
    const auto *smf_408 = buffer.data(smf + 408);
    const auto *smf_409 = buffer.data(smf + 409);
    const auto *smf_410 = buffer.data(smf + 410);
    const auto *smf_412 = buffer.data(smf + 412);
    const auto *smf_416 = buffer.data(smf + 416);
    const auto *smf_417 = buffer.data(smf + 417);
    const auto *smf_418 = buffer.data(smf + 418);
    const auto *smf_419 = buffer.data(smf + 419);
    const auto *smf_420 = buffer.data(smf + 420);
    const auto *smf_422 = buffer.data(smf + 422);
    const auto *smf_426 = buffer.data(smf + 426);
    const auto *smf_427 = buffer.data(smf + 427);
    const auto *smf_428 = buffer.data(smf + 428);
    const auto *smf_429 = buffer.data(smf + 429);
    const auto *smf_430 = buffer.data(smf + 430);
    const auto *smf_432 = buffer.data(smf + 432);
    const auto *smf_436 = buffer.data(smf + 436);
    const auto *smf_437 = buffer.data(smf + 437);
    const auto *smf_438 = buffer.data(smf + 438);
    const auto *smf_439 = buffer.data(smf + 439);
    const auto *smf_440 = buffer.data(smf + 440);
    const auto *smf_442 = buffer.data(smf + 442);
    const auto *smf_446 = buffer.data(smf + 446);
    const auto *smf_447 = buffer.data(smf + 447);
    const auto *smf_448 = buffer.data(smf + 448);
    const auto *smf_449 = buffer.data(smf + 449);
    const auto *smf_450 = buffer.data(smf + 450);
    const auto *smf_452 = buffer.data(smf + 452);
    const auto *smf_453 = buffer.data(smf + 453);
    const auto *smf_455 = buffer.data(smf + 455);
    const auto *smf_456 = buffer.data(smf + 456);
    const auto *smf_457 = buffer.data(smf + 457);
    const auto *smf_458 = buffer.data(smf + 458);
    const auto *smf_459 = buffer.data(smf + 459);
    const auto *smf_460 = buffer.data(smf + 460);
    const auto *smf_462 = buffer.data(smf + 462);
    const auto *smf_465 = buffer.data(smf + 465);
    const auto *smf_466 = buffer.data(smf + 466);
    const auto *smf_467 = buffer.data(smf + 467);
    const auto *smf_468 = buffer.data(smf + 468);
    const auto *smf_469 = buffer.data(smf + 469);
    const auto *smf_470 = buffer.data(smf + 470);
    const auto *smf_472 = buffer.data(smf + 472);
    const auto *smf_473 = buffer.data(smf + 473);
    const auto *smf_475 = buffer.data(smf + 475);
    const auto *smf_476 = buffer.data(smf + 476);
    const auto *smf_477 = buffer.data(smf + 477);
    const auto *smf_478 = buffer.data(smf + 478);
    const auto *smf_479 = buffer.data(smf + 479);
    const auto *smf_480 = buffer.data(smf + 480);
    const auto *smf_482 = buffer.data(smf + 482);
    const auto *smf_483 = buffer.data(smf + 483);
    const auto *smf_485 = buffer.data(smf + 485);
    const auto *smf_486 = buffer.data(smf + 486);
    const auto *smf_487 = buffer.data(smf + 487);
    const auto *smf_488 = buffer.data(smf + 488);

#pragma omp simd aligned(t_598, t_599, t_600, t_601, pb_x, pc_x, pc_y, slg0_599, slg0_600, \
                         slf_319, slf_320, slf_400, slg1_599, slg1_600, smf_399, \
                         smf_400 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_598[k] = f_13 * slf_319[k]
                   + f_3 * pc_y[k] * smf_399[k];

        t_599[k] = pb_x[k] * slg0_599[k]
                   - f_6 * pc_x[k] * slg1_599[k];

        t_600[k] = pb_x[k] * slg0_600[k]
                   + f_14 * slf_400[k]
                   - f_6 * pc_x[k] * slg1_600[k];

        t_601[k] = f_14 * slf_320[k]
                   + f_3 * pc_y[k] * smf_400[k];
    }

#pragma omp simd aligned(t_602, t_603, t_604, pb_x, pc_x, pc_y, pc_z, slg0_603, slf_310, \
                         slf_322, slf_403, slg1_603, smf_400, smf_402 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_602[k] = f_14 * slf_310[k]
                   + f_3 * pc_z[k] * smf_400[k];

        t_603[k] = pb_x[k] * slg0_603[k]
                   + f_8 * slf_403[k]
                   - f_6 * pc_x[k] * slg1_603[k];

        t_604[k] = f_14 * slf_322[k]
                   + f_3 * pc_y[k] * smf_402[k];
    }

#pragma omp simd aligned(t_605, t_606, t_607, t_608, pb_x, pc_x, slg0_605, slf_405, slf_406, \
                         slf_407, slf_408, slg1_605, smf_406, smf_407, \
                         smf_408 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_605[k] = pb_x[k] * slg0_605[k]
                   + f_8 * slf_405[k]
                   - f_6 * pc_x[k] * slg1_605[k];

        t_606[k] = f_7 * slf_406[k]
                   + f_3 * pc_x[k] * smf_406[k];

        t_607[k] = f_7 * slf_407[k]
                   + f_3 * pc_x[k] * smf_407[k];

        t_608[k] = f_7 * slf_408[k]
                   + f_3 * pc_x[k] * smf_408[k];
    }

#pragma omp simd aligned(t_609, t_610, t_611, t_612, pb_x, pc_x, pc_z, slg0_610, slg0_612, \
                         slf_316, slf_409, slg1_610, slg1_612, smf_406, \
                         smf_409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_609[k] = f_7 * slf_409[k]
                   + f_3 * pc_x[k] * smf_409[k];

        t_610[k] = pb_x[k] * slg0_610[k]
                   - f_6 * pc_x[k] * slg1_610[k];

        t_611[k] = f_14 * slf_316[k]
                   + f_3 * pc_z[k] * smf_406[k];

        t_612[k] = pb_x[k] * slg0_612[k]
                   - f_6 * pc_x[k] * slg1_612[k];
    }

#pragma omp simd aligned(t_613, t_614, t_615, t_616, pb_x, pc_x, pc_y, slg0_614, slg0_615, \
                         slf_329, slf_330, slf_410, slg1_614, slg1_615, smf_409, \
                         smf_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_613[k] = f_14 * slf_329[k]
                   + f_3 * pc_y[k] * smf_409[k];

        t_614[k] = pb_x[k] * slg0_614[k]
                   - f_6 * pc_x[k] * slg1_614[k];

        t_615[k] = pb_x[k] * slg0_615[k]
                   + f_14 * slf_410[k]
                   - f_6 * pc_x[k] * slg1_615[k];

        t_616[k] = f_12 * slf_330[k]
                   + f_3 * pc_y[k] * smf_410[k];
    }

#pragma omp simd aligned(t_617, t_618, t_619, pb_x, pc_x, pc_y, pc_z, slg0_618, slf_320, \
                         slf_332, slf_413, slg1_618, smf_410, smf_412 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_617[k] = f_13 * slf_320[k]
                   + f_3 * pc_z[k] * smf_410[k];

        t_618[k] = pb_x[k] * slg0_618[k]
                   + f_8 * slf_413[k]
                   - f_6 * pc_x[k] * slg1_618[k];

        t_619[k] = f_12 * slf_332[k]
                   + f_3 * pc_y[k] * smf_412[k];
    }

#pragma omp simd aligned(t_620, t_621, t_622, t_623, pb_x, pc_x, slg0_620, slf_415, slf_416, \
                         slf_417, slf_418, slg1_620, smf_416, smf_417, \
                         smf_418 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_620[k] = pb_x[k] * slg0_620[k]
                   + f_8 * slf_415[k]
                   - f_6 * pc_x[k] * slg1_620[k];

        t_621[k] = f_7 * slf_416[k]
                   + f_3 * pc_x[k] * smf_416[k];

        t_622[k] = f_7 * slf_417[k]
                   + f_3 * pc_x[k] * smf_417[k];

        t_623[k] = f_7 * slf_418[k]
                   + f_3 * pc_x[k] * smf_418[k];
    }

#pragma omp simd aligned(t_624, t_625, t_626, t_627, pb_x, pc_x, pc_z, slg0_625, slg0_627, \
                         slf_326, slf_419, slg1_625, slg1_627, smf_416, \
                         smf_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_624[k] = f_7 * slf_419[k]
                   + f_3 * pc_x[k] * smf_419[k];

        t_625[k] = pb_x[k] * slg0_625[k]
                   - f_6 * pc_x[k] * slg1_625[k];

        t_626[k] = f_13 * slf_326[k]
                   + f_3 * pc_z[k] * smf_416[k];

        t_627[k] = pb_x[k] * slg0_627[k]
                   - f_6 * pc_x[k] * slg1_627[k];
    }

#pragma omp simd aligned(t_628, t_629, t_630, t_631, pb_x, pc_x, pc_y, slg0_629, slg0_630, \
                         slf_339, slf_340, slf_420, slg1_629, slg1_630, smf_419, \
                         smf_420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_628[k] = f_12 * slf_339[k]
                   + f_3 * pc_y[k] * smf_419[k];

        t_629[k] = pb_x[k] * slg0_629[k]
                   - f_6 * pc_x[k] * slg1_629[k];

        t_630[k] = pb_x[k] * slg0_630[k]
                   + f_14 * slf_420[k]
                   - f_6 * pc_x[k] * slg1_630[k];

        t_631[k] = f_8 * slf_340[k]
                   + f_3 * pc_y[k] * smf_420[k];
    }

#pragma omp simd aligned(t_632, t_633, t_634, pb_x, pc_x, pc_y, pc_z, slg0_633, slf_330, \
                         slf_342, slf_423, slg1_633, smf_420, smf_422 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_632[k] = f_11 * slf_330[k]
                   + f_3 * pc_z[k] * smf_420[k];

        t_633[k] = pb_x[k] * slg0_633[k]
                   + f_8 * slf_423[k]
                   - f_6 * pc_x[k] * slg1_633[k];

        t_634[k] = f_8 * slf_342[k]
                   + f_3 * pc_y[k] * smf_422[k];
    }

#pragma omp simd aligned(t_635, t_636, t_637, t_638, pb_x, pc_x, slg0_635, slf_425, slf_426, \
                         slf_427, slf_428, slg1_635, smf_426, smf_427, \
                         smf_428 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_635[k] = pb_x[k] * slg0_635[k]
                   + f_8 * slf_425[k]
                   - f_6 * pc_x[k] * slg1_635[k];

        t_636[k] = f_7 * slf_426[k]
                   + f_3 * pc_x[k] * smf_426[k];

        t_637[k] = f_7 * slf_427[k]
                   + f_3 * pc_x[k] * smf_427[k];

        t_638[k] = f_7 * slf_428[k]
                   + f_3 * pc_x[k] * smf_428[k];
    }

#pragma omp simd aligned(t_639, t_640, t_641, t_642, pb_x, pc_x, pc_z, slg0_640, slg0_642, \
                         slf_336, slf_429, slg1_640, slg1_642, smf_426, \
                         smf_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_639[k] = f_7 * slf_429[k]
                   + f_3 * pc_x[k] * smf_429[k];

        t_640[k] = pb_x[k] * slg0_640[k]
                   - f_6 * pc_x[k] * slg1_640[k];

        t_641[k] = f_11 * slf_336[k]
                   + f_3 * pc_z[k] * smf_426[k];

        t_642[k] = pb_x[k] * slg0_642[k]
                   - f_6 * pc_x[k] * slg1_642[k];
    }

#pragma omp simd aligned(t_643, t_644, t_645, t_646, pb_x, pb_y, pc_x, pc_y, slg0_525, \
                         slg0_644, slf_349, slf_350, slg1_525, slg1_644, smf_429, \
                         smf_430 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_643[k] = f_8 * slf_349[k]
                   + f_3 * pc_y[k] * smf_429[k];

        t_644[k] = pb_x[k] * slg0_644[k]
                   - f_6 * pc_x[k] * slg1_644[k];

        t_645[k] = pb_y[k] * slg0_525[k]
                   - f_6 * pc_y[k] * slg1_525[k];

        t_646[k] = f_7 * slf_350[k]
                   + f_3 * pc_y[k] * smf_430[k];
    }

#pragma omp simd aligned(t_647, t_648, t_649, pb_x, pc_x, pc_y, pc_z, slg0_648, slf_340, \
                         slf_352, slf_433, slg1_648, smf_430, smf_432 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_647[k] = f_10 * slf_340[k]
                   + f_3 * pc_z[k] * smf_430[k];

        t_648[k] = pb_x[k] * slg0_648[k]
                   + f_8 * slf_433[k]
                   - f_6 * pc_x[k] * slg1_648[k];

        t_649[k] = f_7 * slf_352[k]
                   + f_3 * pc_y[k] * smf_432[k];
    }

#pragma omp simd aligned(t_650, t_651, t_652, t_653, pb_y, pc_x, pc_y, slg0_530, slf_436, \
                         slf_437, slf_438, slg1_530, smf_436, smf_437, \
                         smf_438 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_650[k] = pb_y[k] * slg0_530[k]
                   - f_6 * pc_y[k] * slg1_530[k];

        t_651[k] = f_7 * slf_436[k]
                   + f_3 * pc_x[k] * smf_436[k];

        t_652[k] = f_7 * slf_437[k]
                   + f_3 * pc_x[k] * smf_437[k];

        t_653[k] = f_7 * slf_438[k]
                   + f_3 * pc_x[k] * smf_438[k];
    }

#pragma omp simd aligned(t_654, t_655, t_656, t_657, pb_x, pc_x, pc_z, slg0_655, slg0_657, \
                         slf_346, slf_439, slg1_655, slg1_657, smf_436, \
                         smf_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_654[k] = f_7 * slf_439[k]
                   + f_3 * pc_x[k] * smf_439[k];

        t_655[k] = pb_x[k] * slg0_655[k]
                   - f_6 * pc_x[k] * slg1_655[k];

        t_656[k] = f_10 * slf_346[k]
                   + f_3 * pc_z[k] * smf_436[k];

        t_657[k] = pb_x[k] * slg0_657[k]
                   - f_6 * pc_x[k] * slg1_657[k];
    }

#pragma omp simd aligned(t_658, t_659, t_660, t_661, pb_x, pc_x, pc_y, slg0_659, slg0_660, \
                         slf_359, slf_440, slg1_659, slg1_660, smf_439, \
                         smf_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_658[k] = f_7 * slf_359[k]
                   + f_3 * pc_y[k] * smf_439[k];

        t_659[k] = pb_x[k] * slg0_659[k]
                   - f_6 * pc_x[k] * slg1_659[k];

        t_660[k] = pb_x[k] * slg0_660[k]
                   + f_14 * slf_440[k]
                   - f_6 * pc_x[k] * slg1_660[k];

        t_661[k] = f_3 * pc_y[k] * smf_440[k];
    }

#pragma omp simd aligned(t_662, t_663, t_664, pb_x, pc_x, pc_y, pc_z, slg0_663, slf_350, \
                         slf_443, slg1_663, smf_440, smf_442 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_662[k] = f_9 * slf_350[k]
                   + f_3 * pc_z[k] * smf_440[k];

        t_663[k] = pb_x[k] * slg0_663[k]
                   + f_8 * slf_443[k]
                   - f_6 * pc_x[k] * slg1_663[k];

        t_664[k] = f_3 * pc_y[k] * smf_442[k];
    }

#pragma omp simd aligned(t_665, t_666, t_667, t_668, pb_x, pc_x, slg0_665, slf_445, slf_446, \
                         slf_447, slf_448, slg1_665, smf_446, smf_447, \
                         smf_448 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_665[k] = pb_x[k] * slg0_665[k]
                   + f_8 * slf_445[k]
                   - f_6 * pc_x[k] * slg1_665[k];

        t_666[k] = f_7 * slf_446[k]
                   + f_3 * pc_x[k] * smf_446[k];

        t_667[k] = f_7 * slf_447[k]
                   + f_3 * pc_x[k] * smf_447[k];

        t_668[k] = f_7 * slf_448[k]
                   + f_3 * pc_x[k] * smf_448[k];
    }

#pragma omp simd aligned(t_669, t_670, t_671, t_672, pb_x, pc_x, pc_z, slg0_670, slg0_672, \
                         slf_356, slf_449, slg1_670, slg1_672, smf_446, \
                         smf_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_669[k] = f_7 * slf_449[k]
                   + f_3 * pc_x[k] * smf_449[k];

        t_670[k] = pb_x[k] * slg0_670[k]
                   - f_6 * pc_x[k] * slg1_670[k];

        t_671[k] = f_9 * slf_356[k]
                   + f_3 * pc_z[k] * smf_446[k];

        t_672[k] = pb_x[k] * slg0_672[k]
                   - f_6 * pc_x[k] * slg1_672[k];
    }

#pragma omp simd aligned(t_673, t_674, t_675, t_676, t_677, pb_x, pc_x, pc_y, pc_z, slg0_674, \
                         slf_360, slg1_674, smd0_270, smd1_270, smf_449, \
                         smf_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_673[k] = f_3 * pc_y[k] * smf_449[k];

        t_674[k] = pb_x[k] * slg0_674[k]
                   - f_6 * pc_x[k] * slg1_674[k];

        t_675[k] = f_1 * smd0_270[k]
                   - f_2 * smd1_270[k]
                   + f_3 * pc_x[k] * smf_450[k];

        t_676[k] = f_0 * slf_360[k]
                   + f_3 * pc_y[k] * smf_450[k];

        t_677[k] = f_3 * pc_z[k] * smf_450[k];
    }

#pragma omp simd aligned(t_678, t_679, t_680, t_681, pc_x, pc_y, slf_362, smd0_273, smd0_275, \
                         smd1_273, smd1_275, smf_452, smf_453, smf_455, \
                         smf_456 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_678[k] = f_4 * smd0_273[k]
                   - f_5 * smd1_273[k]
                   + f_3 * pc_x[k] * smf_453[k];

        t_679[k] = f_0 * slf_362[k]
                   + f_3 * pc_y[k] * smf_452[k];

        t_680[k] = f_4 * smd0_275[k]
                   - f_5 * smd1_275[k]
                   + f_3 * pc_x[k] * smf_455[k];

        t_681[k] = f_3 * pc_x[k] * smf_456[k];
    }

#pragma omp simd aligned(t_682, t_683, t_684, t_685, t_686, pc_x, pc_y, pc_z, slf_366, \
                         smd0_273, smd1_273, smf_456, smf_457, smf_458, \
                         smf_459 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_682[k] = f_3 * pc_x[k] * smf_457[k];

        t_683[k] = f_3 * pc_x[k] * smf_458[k];

        t_684[k] = f_3 * pc_x[k] * smf_459[k];

        t_685[k] = f_0 * slf_366[k]
                   + f_1 * smd0_273[k]
                   - f_2 * smd1_273[k]
                   + f_3 * pc_y[k] * smf_456[k];

        t_686[k] = f_3 * pc_z[k] * smf_456[k];
    }

#pragma omp simd aligned(t_687, t_688, t_689, t_690, pb_z, pc_y, pc_z, slg0_540, slf_368, \
                         slf_369, slg1_540, smd0_275, smd1_275, smf_458, \
                         smf_459 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_687[k] = f_0 * slf_368[k]
                   + f_4 * smd0_275[k]
                   - f_5 * smd1_275[k]
                   + f_3 * pc_y[k] * smf_458[k];

        t_688[k] = f_0 * slf_369[k]
                   + f_3 * pc_y[k] * smf_459[k];

        t_689[k] = f_1 * smd0_275[k]
                   - f_2 * smd1_275[k]
                   + f_3 * pc_z[k] * smf_459[k];

        t_690[k] = pb_z[k] * slg0_540[k]
                   - f_6 * pc_z[k] * slg1_540[k];
    }

#pragma omp simd aligned(t_691, t_692, t_693, t_694, pb_z, pc_y, pc_z, slg0_543, slf_360, \
                         slf_370, slf_372, slg1_543, smf_460, smf_462 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_691[k] = f_9 * slf_370[k]
                   + f_3 * pc_y[k] * smf_460[k];

        t_692[k] = f_7 * slf_360[k]
                   + f_3 * pc_z[k] * smf_460[k];

        t_693[k] = pb_z[k] * slg0_543[k]
                   - f_6 * pc_z[k] * slg1_543[k];

        t_694[k] = f_9 * slf_372[k]
                   + f_3 * pc_y[k] * smf_462[k];
    }

#pragma omp simd aligned(t_695, t_696, t_697, t_698, t_699, pc_x, smd0_281, smd1_281, smf_465, \
                         smf_466, smf_467, smf_468, smf_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_695[k] = f_4 * smd0_281[k]
                   - f_5 * smd1_281[k]
                   + f_3 * pc_x[k] * smf_465[k];

        t_696[k] = f_3 * pc_x[k] * smf_466[k];

        t_697[k] = f_3 * pc_x[k] * smf_467[k];

        t_698[k] = f_3 * pc_x[k] * smf_468[k];

        t_699[k] = f_3 * pc_x[k] * smf_469[k];
    }

#pragma omp simd aligned(t_700, t_701, t_702, t_703, pb_z, pc_y, pc_z, slg0_550, slg0_552, \
                         slf_366, slf_367, slf_379, slg1_550, slg1_552, smf_466, \
                         smf_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_700[k] = pb_z[k] * slg0_550[k]
                   - f_6 * pc_z[k] * slg1_550[k];

        t_701[k] = f_7 * slf_366[k]
                   + f_3 * pc_z[k] * smf_466[k];

        t_702[k] = pb_z[k] * slg0_552[k]
                   + f_8 * slf_367[k]
                   - f_6 * pc_z[k] * slg1_552[k];

        t_703[k] = f_9 * slf_379[k]
                   + f_3 * pc_y[k] * smf_469[k];
    }

#pragma omp simd aligned(t_704, t_705, t_706, t_707, pc_x, pc_y, pc_z, slf_369, slf_370, \
                         slf_380, smd0_281, smd0_282, smd1_281, smd1_282, smf_469, \
                         smf_470 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_704[k] = f_7 * slf_369[k]
                   + f_1 * smd0_281[k]
                   - f_2 * smd1_281[k]
                   + f_3 * pc_z[k] * smf_469[k];

        t_705[k] = f_1 * smd0_282[k]
                   - f_2 * smd1_282[k]
                   + f_3 * pc_x[k] * smf_470[k];

        t_706[k] = f_10 * slf_380[k]
                   + f_3 * pc_y[k] * smf_470[k];

        t_707[k] = f_8 * slf_370[k]
                   + f_3 * pc_z[k] * smf_470[k];
    }

#pragma omp simd aligned(t_708, t_709, t_710, t_711, pc_x, pc_y, slf_382, smd0_285, smd0_287, \
                         smd1_285, smd1_287, smf_472, smf_473, smf_475, \
                         smf_476 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_708[k] = f_4 * smd0_285[k]
                   - f_5 * smd1_285[k]
                   + f_3 * pc_x[k] * smf_473[k];

        t_709[k] = f_10 * slf_382[k]
                   + f_3 * pc_y[k] * smf_472[k];

        t_710[k] = f_4 * smd0_287[k]
                   - f_5 * smd1_287[k]
                   + f_3 * pc_x[k] * smf_475[k];

        t_711[k] = f_3 * pc_x[k] * smf_476[k];
    }

#pragma omp simd aligned(t_712, t_713, t_714, t_715, t_716, pc_x, pc_y, pc_z, slf_376, \
                         slf_386, smd0_285, smd1_285, smf_476, smf_477, smf_478, \
                         smf_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_712[k] = f_3 * pc_x[k] * smf_477[k];

        t_713[k] = f_3 * pc_x[k] * smf_478[k];

        t_714[k] = f_3 * pc_x[k] * smf_479[k];

        t_715[k] = f_10 * slf_386[k]
                   + f_1 * smd0_285[k]
                   - f_2 * smd1_285[k]
                   + f_3 * pc_y[k] * smf_476[k];

        t_716[k] = f_8 * slf_376[k]
                   + f_3 * pc_z[k] * smf_476[k];
    }

#pragma omp simd aligned(t_717, t_718, t_719, pc_y, pc_z, slf_379, slf_388, slf_389, smd0_287, \
                         smd1_287, smf_478, smf_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_717[k] = f_10 * slf_388[k]
                   + f_4 * smd0_287[k]
                   - f_5 * smd1_287[k]
                   + f_3 * pc_y[k] * smf_478[k];

        t_718[k] = f_10 * slf_389[k]
                   + f_3 * pc_y[k] * smf_479[k];

        t_719[k] = f_8 * slf_379[k]
                   + f_1 * smd0_287[k]
                   - f_2 * smd1_287[k]
                   + f_3 * pc_z[k] * smf_479[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, t_723, pc_x, pc_y, pc_z, slf_380, slf_390, \
                         smd0_288, smd0_291, smd1_288, smd1_291, smf_480, \
                         smf_483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = f_1 * smd0_288[k]
                   - f_2 * smd1_288[k]
                   + f_3 * pc_x[k] * smf_480[k];

        t_721[k] = f_11 * slf_390[k]
                   + f_3 * pc_y[k] * smf_480[k];

        t_722[k] = f_12 * slf_380[k]
                   + f_3 * pc_z[k] * smf_480[k];

        t_723[k] = f_4 * smd0_291[k]
                   - f_5 * smd1_291[k]
                   + f_3 * pc_x[k] * smf_483[k];
    }

#pragma omp simd aligned(t_724, t_725, t_726, t_727, t_728, pc_x, pc_y, slf_392, smd0_293, \
                         smd1_293, smf_482, smf_485, smf_486, smf_487, \
                         smf_488 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_724[k] = f_11 * slf_392[k]
                   + f_3 * pc_y[k] * smf_482[k];

        t_725[k] = f_4 * smd0_293[k]
                   - f_5 * smd1_293[k]
                   + f_3 * pc_x[k] * smf_485[k];

        t_726[k] = f_3 * pc_x[k] * smf_486[k];

        t_727[k] = f_3 * pc_x[k] * smf_487[k];

        t_728[k] = f_3 * pc_x[k] * smf_488[k];
    }
}

static auto
compute_prim_smg_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t slg0,
                                                          const size_t slf, const size_t slg1,
                                                          const size_t smd0, const size_t smd1,
                                                          const size_t smf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 4.0 / q;
    const auto f_10 = 3.5 / q;
    const auto f_11 = 3.0 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 2.5 / q;
    const auto f_14 = 2.0 / q;

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

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *slg0_660 = buffer.data(slg0 + 660);
    const auto *slg0_665 = buffer.data(slg0 + 665);
    const auto *slg0_670 = buffer.data(slg0 + 670);
    const auto *slg0_672 = buffer.data(slg0 + 672);
    const auto *slg0_674 = buffer.data(slg0 + 674);

    const auto *slf_386 = buffer.data(slf + 386);
    const auto *slf_389 = buffer.data(slf + 389);
    const auto *slf_390 = buffer.data(slf + 390);
    const auto *slf_396 = buffer.data(slf + 396);
    const auto *slf_398 = buffer.data(slf + 398);
    const auto *slf_399 = buffer.data(slf + 399);
    const auto *slf_400 = buffer.data(slf + 400);
    const auto *slf_402 = buffer.data(slf + 402);
    const auto *slf_406 = buffer.data(slf + 406);
    const auto *slf_408 = buffer.data(slf + 408);
    const auto *slf_409 = buffer.data(slf + 409);
    const auto *slf_410 = buffer.data(slf + 410);
    const auto *slf_412 = buffer.data(slf + 412);
    const auto *slf_416 = buffer.data(slf + 416);
    const auto *slf_418 = buffer.data(slf + 418);
    const auto *slf_419 = buffer.data(slf + 419);
    const auto *slf_420 = buffer.data(slf + 420);
    const auto *slf_422 = buffer.data(slf + 422);
    const auto *slf_426 = buffer.data(slf + 426);
    const auto *slf_428 = buffer.data(slf + 428);
    const auto *slf_429 = buffer.data(slf + 429);
    const auto *slf_430 = buffer.data(slf + 430);
    const auto *slf_432 = buffer.data(slf + 432);
    const auto *slf_436 = buffer.data(slf + 436);
    const auto *slf_438 = buffer.data(slf + 438);
    const auto *slf_439 = buffer.data(slf + 439);
    const auto *slf_440 = buffer.data(slf + 440);
    const auto *slf_442 = buffer.data(slf + 442);
    const auto *slf_446 = buffer.data(slf + 446);
    const auto *slf_448 = buffer.data(slf + 448);
    const auto *slf_449 = buffer.data(slf + 449);

    const auto *slg1_660 = buffer.data(slg1 + 660);
    const auto *slg1_665 = buffer.data(slg1 + 665);
    const auto *slg1_670 = buffer.data(slg1 + 670);
    const auto *slg1_672 = buffer.data(slg1 + 672);
    const auto *slg1_674 = buffer.data(slg1 + 674);

    const auto *smd0_291 = buffer.data(smd0 + 291);
    const auto *smd0_293 = buffer.data(smd0 + 293);
    const auto *smd0_294 = buffer.data(smd0 + 294);
    const auto *smd0_297 = buffer.data(smd0 + 297);
    const auto *smd0_299 = buffer.data(smd0 + 299);
    const auto *smd0_300 = buffer.data(smd0 + 300);
    const auto *smd0_303 = buffer.data(smd0 + 303);
    const auto *smd0_305 = buffer.data(smd0 + 305);
    const auto *smd0_306 = buffer.data(smd0 + 306);
    const auto *smd0_309 = buffer.data(smd0 + 309);
    const auto *smd0_311 = buffer.data(smd0 + 311);
    const auto *smd0_312 = buffer.data(smd0 + 312);
    const auto *smd0_315 = buffer.data(smd0 + 315);
    const auto *smd0_317 = buffer.data(smd0 + 317);
    const auto *smd0_321 = buffer.data(smd0 + 321);
    const auto *smd0_324 = buffer.data(smd0 + 324);
    const auto *smd0_327 = buffer.data(smd0 + 327);
    const auto *smd0_329 = buffer.data(smd0 + 329);

    const auto *smd1_291 = buffer.data(smd1 + 291);
    const auto *smd1_293 = buffer.data(smd1 + 293);
    const auto *smd1_294 = buffer.data(smd1 + 294);
    const auto *smd1_297 = buffer.data(smd1 + 297);
    const auto *smd1_299 = buffer.data(smd1 + 299);
    const auto *smd1_300 = buffer.data(smd1 + 300);
    const auto *smd1_303 = buffer.data(smd1 + 303);
    const auto *smd1_305 = buffer.data(smd1 + 305);
    const auto *smd1_306 = buffer.data(smd1 + 306);
    const auto *smd1_309 = buffer.data(smd1 + 309);
    const auto *smd1_311 = buffer.data(smd1 + 311);
    const auto *smd1_312 = buffer.data(smd1 + 312);
    const auto *smd1_315 = buffer.data(smd1 + 315);
    const auto *smd1_317 = buffer.data(smd1 + 317);
    const auto *smd1_321 = buffer.data(smd1 + 321);
    const auto *smd1_324 = buffer.data(smd1 + 324);
    const auto *smd1_327 = buffer.data(smd1 + 327);
    const auto *smd1_329 = buffer.data(smd1 + 329);

    const auto *smf_486 = buffer.data(smf + 486);
    const auto *smf_488 = buffer.data(smf + 488);
    const auto *smf_489 = buffer.data(smf + 489);
    const auto *smf_490 = buffer.data(smf + 490);
    const auto *smf_492 = buffer.data(smf + 492);
    const auto *smf_493 = buffer.data(smf + 493);
    const auto *smf_495 = buffer.data(smf + 495);
    const auto *smf_496 = buffer.data(smf + 496);
    const auto *smf_497 = buffer.data(smf + 497);
    const auto *smf_498 = buffer.data(smf + 498);
    const auto *smf_499 = buffer.data(smf + 499);
    const auto *smf_500 = buffer.data(smf + 500);
    const auto *smf_502 = buffer.data(smf + 502);
    const auto *smf_503 = buffer.data(smf + 503);
    const auto *smf_505 = buffer.data(smf + 505);
    const auto *smf_506 = buffer.data(smf + 506);
    const auto *smf_507 = buffer.data(smf + 507);
    const auto *smf_508 = buffer.data(smf + 508);
    const auto *smf_509 = buffer.data(smf + 509);
    const auto *smf_510 = buffer.data(smf + 510);
    const auto *smf_512 = buffer.data(smf + 512);
    const auto *smf_513 = buffer.data(smf + 513);
    const auto *smf_515 = buffer.data(smf + 515);
    const auto *smf_516 = buffer.data(smf + 516);
    const auto *smf_517 = buffer.data(smf + 517);
    const auto *smf_518 = buffer.data(smf + 518);
    const auto *smf_519 = buffer.data(smf + 519);
    const auto *smf_520 = buffer.data(smf + 520);
    const auto *smf_522 = buffer.data(smf + 522);
    const auto *smf_523 = buffer.data(smf + 523);
    const auto *smf_525 = buffer.data(smf + 525);
    const auto *smf_526 = buffer.data(smf + 526);
    const auto *smf_527 = buffer.data(smf + 527);
    const auto *smf_528 = buffer.data(smf + 528);
    const auto *smf_529 = buffer.data(smf + 529);
    const auto *smf_530 = buffer.data(smf + 530);
    const auto *smf_532 = buffer.data(smf + 532);
    const auto *smf_533 = buffer.data(smf + 533);
    const auto *smf_536 = buffer.data(smf + 536);
    const auto *smf_537 = buffer.data(smf + 537);
    const auto *smf_538 = buffer.data(smf + 538);
    const auto *smf_539 = buffer.data(smf + 539);
    const auto *smf_540 = buffer.data(smf + 540);
    const auto *smf_542 = buffer.data(smf + 542);
    const auto *smf_543 = buffer.data(smf + 543);
    const auto *smf_545 = buffer.data(smf + 545);
    const auto *smf_546 = buffer.data(smf + 546);
    const auto *smf_547 = buffer.data(smf + 547);
    const auto *smf_548 = buffer.data(smf + 548);
    const auto *smf_549 = buffer.data(smf + 549);

#pragma omp simd aligned(t_729, t_730, t_731, pc_x, pc_y, pc_z, slf_386, slf_396, smd0_291, \
                         smd1_291, smf_486, smf_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_729[k] = f_3 * pc_x[k] * smf_489[k];

        t_730[k] = f_11 * slf_396[k]
                   + f_1 * smd0_291[k]
                   - f_2 * smd1_291[k]
                   + f_3 * pc_y[k] * smf_486[k];

        t_731[k] = f_12 * slf_386[k]
                   + f_3 * pc_z[k] * smf_486[k];
    }

#pragma omp simd aligned(t_732, t_733, t_734, pc_y, pc_z, slf_389, slf_398, slf_399, smd0_293, \
                         smd1_293, smf_488, smf_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_732[k] = f_11 * slf_398[k]
                   + f_4 * smd0_293[k]
                   - f_5 * smd1_293[k]
                   + f_3 * pc_y[k] * smf_488[k];

        t_733[k] = f_11 * slf_399[k]
                   + f_3 * pc_y[k] * smf_489[k];

        t_734[k] = f_12 * slf_389[k]
                   + f_1 * smd0_293[k]
                   - f_2 * smd1_293[k]
                   + f_3 * pc_z[k] * smf_489[k];
    }

#pragma omp simd aligned(t_735, t_736, t_737, t_738, pc_x, pc_y, pc_z, slf_390, slf_400, \
                         smd0_294, smd0_297, smd1_294, smd1_297, smf_490, \
                         smf_493 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_735[k] = f_1 * smd0_294[k]
                   - f_2 * smd1_294[k]
                   + f_3 * pc_x[k] * smf_490[k];

        t_736[k] = f_13 * slf_400[k]
                   + f_3 * pc_y[k] * smf_490[k];

        t_737[k] = f_14 * slf_390[k]
                   + f_3 * pc_z[k] * smf_490[k];

        t_738[k] = f_4 * smd0_297[k]
                   - f_5 * smd1_297[k]
                   + f_3 * pc_x[k] * smf_493[k];
    }

#pragma omp simd aligned(t_739, t_740, t_741, t_742, t_743, pc_x, pc_y, slf_402, smd0_299, \
                         smd1_299, smf_492, smf_495, smf_496, smf_497, \
                         smf_498 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_739[k] = f_13 * slf_402[k]
                   + f_3 * pc_y[k] * smf_492[k];

        t_740[k] = f_4 * smd0_299[k]
                   - f_5 * smd1_299[k]
                   + f_3 * pc_x[k] * smf_495[k];

        t_741[k] = f_3 * pc_x[k] * smf_496[k];

        t_742[k] = f_3 * pc_x[k] * smf_497[k];

        t_743[k] = f_3 * pc_x[k] * smf_498[k];
    }

#pragma omp simd aligned(t_744, t_745, t_746, pc_x, pc_y, pc_z, slf_396, slf_406, smd0_297, \
                         smd1_297, smf_496, smf_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_744[k] = f_3 * pc_x[k] * smf_499[k];

        t_745[k] = f_13 * slf_406[k]
                   + f_1 * smd0_297[k]
                   - f_2 * smd1_297[k]
                   + f_3 * pc_y[k] * smf_496[k];

        t_746[k] = f_14 * slf_396[k]
                   + f_3 * pc_z[k] * smf_496[k];
    }

#pragma omp simd aligned(t_747, t_748, t_749, pc_y, pc_z, slf_399, slf_408, slf_409, smd0_299, \
                         smd1_299, smf_498, smf_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_747[k] = f_13 * slf_408[k]
                   + f_4 * smd0_299[k]
                   - f_5 * smd1_299[k]
                   + f_3 * pc_y[k] * smf_498[k];

        t_748[k] = f_13 * slf_409[k]
                   + f_3 * pc_y[k] * smf_499[k];

        t_749[k] = f_14 * slf_399[k]
                   + f_1 * smd0_299[k]
                   - f_2 * smd1_299[k]
                   + f_3 * pc_z[k] * smf_499[k];
    }

#pragma omp simd aligned(t_750, t_751, t_752, t_753, pc_x, pc_y, pc_z, slf_400, slf_410, \
                         smd0_300, smd0_303, smd1_300, smd1_303, smf_500, \
                         smf_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_750[k] = f_1 * smd0_300[k]
                   - f_2 * smd1_300[k]
                   + f_3 * pc_x[k] * smf_500[k];

        t_751[k] = f_14 * slf_410[k]
                   + f_3 * pc_y[k] * smf_500[k];

        t_752[k] = f_13 * slf_400[k]
                   + f_3 * pc_z[k] * smf_500[k];

        t_753[k] = f_4 * smd0_303[k]
                   - f_5 * smd1_303[k]
                   + f_3 * pc_x[k] * smf_503[k];
    }

#pragma omp simd aligned(t_754, t_755, t_756, t_757, t_758, pc_x, pc_y, slf_412, smd0_305, \
                         smd1_305, smf_502, smf_505, smf_506, smf_507, \
                         smf_508 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_754[k] = f_14 * slf_412[k]
                   + f_3 * pc_y[k] * smf_502[k];

        t_755[k] = f_4 * smd0_305[k]
                   - f_5 * smd1_305[k]
                   + f_3 * pc_x[k] * smf_505[k];

        t_756[k] = f_3 * pc_x[k] * smf_506[k];

        t_757[k] = f_3 * pc_x[k] * smf_507[k];

        t_758[k] = f_3 * pc_x[k] * smf_508[k];
    }

#pragma omp simd aligned(t_759, t_760, t_761, pc_x, pc_y, pc_z, slf_406, slf_416, smd0_303, \
                         smd1_303, smf_506, smf_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_759[k] = f_3 * pc_x[k] * smf_509[k];

        t_760[k] = f_14 * slf_416[k]
                   + f_1 * smd0_303[k]
                   - f_2 * smd1_303[k]
                   + f_3 * pc_y[k] * smf_506[k];

        t_761[k] = f_13 * slf_406[k]
                   + f_3 * pc_z[k] * smf_506[k];
    }

#pragma omp simd aligned(t_762, t_763, t_764, pc_y, pc_z, slf_409, slf_418, slf_419, smd0_305, \
                         smd1_305, smf_508, smf_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_762[k] = f_14 * slf_418[k]
                   + f_4 * smd0_305[k]
                   - f_5 * smd1_305[k]
                   + f_3 * pc_y[k] * smf_508[k];

        t_763[k] = f_14 * slf_419[k]
                   + f_3 * pc_y[k] * smf_509[k];

        t_764[k] = f_13 * slf_409[k]
                   + f_1 * smd0_305[k]
                   - f_2 * smd1_305[k]
                   + f_3 * pc_z[k] * smf_509[k];
    }

#pragma omp simd aligned(t_765, t_766, t_767, t_768, pc_x, pc_y, pc_z, slf_410, slf_420, \
                         smd0_306, smd0_309, smd1_306, smd1_309, smf_510, \
                         smf_513 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_765[k] = f_1 * smd0_306[k]
                   - f_2 * smd1_306[k]
                   + f_3 * pc_x[k] * smf_510[k];

        t_766[k] = f_12 * slf_420[k]
                   + f_3 * pc_y[k] * smf_510[k];

        t_767[k] = f_11 * slf_410[k]
                   + f_3 * pc_z[k] * smf_510[k];

        t_768[k] = f_4 * smd0_309[k]
                   - f_5 * smd1_309[k]
                   + f_3 * pc_x[k] * smf_513[k];
    }

#pragma omp simd aligned(t_769, t_770, t_771, t_772, t_773, pc_x, pc_y, slf_422, smd0_311, \
                         smd1_311, smf_512, smf_515, smf_516, smf_517, \
                         smf_518 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_769[k] = f_12 * slf_422[k]
                   + f_3 * pc_y[k] * smf_512[k];

        t_770[k] = f_4 * smd0_311[k]
                   - f_5 * smd1_311[k]
                   + f_3 * pc_x[k] * smf_515[k];

        t_771[k] = f_3 * pc_x[k] * smf_516[k];

        t_772[k] = f_3 * pc_x[k] * smf_517[k];

        t_773[k] = f_3 * pc_x[k] * smf_518[k];
    }

#pragma omp simd aligned(t_774, t_775, t_776, pc_x, pc_y, pc_z, slf_416, slf_426, smd0_309, \
                         smd1_309, smf_516, smf_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_774[k] = f_3 * pc_x[k] * smf_519[k];

        t_775[k] = f_12 * slf_426[k]
                   + f_1 * smd0_309[k]
                   - f_2 * smd1_309[k]
                   + f_3 * pc_y[k] * smf_516[k];

        t_776[k] = f_11 * slf_416[k]
                   + f_3 * pc_z[k] * smf_516[k];
    }

#pragma omp simd aligned(t_777, t_778, t_779, pc_y, pc_z, slf_419, slf_428, slf_429, smd0_311, \
                         smd1_311, smf_518, smf_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_777[k] = f_12 * slf_428[k]
                   + f_4 * smd0_311[k]
                   - f_5 * smd1_311[k]
                   + f_3 * pc_y[k] * smf_518[k];

        t_778[k] = f_12 * slf_429[k]
                   + f_3 * pc_y[k] * smf_519[k];

        t_779[k] = f_11 * slf_419[k]
                   + f_1 * smd0_311[k]
                   - f_2 * smd1_311[k]
                   + f_3 * pc_z[k] * smf_519[k];
    }

#pragma omp simd aligned(t_780, t_781, t_782, t_783, pc_x, pc_y, pc_z, slf_420, slf_430, \
                         smd0_312, smd0_315, smd1_312, smd1_315, smf_520, \
                         smf_523 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_780[k] = f_1 * smd0_312[k]
                   - f_2 * smd1_312[k]
                   + f_3 * pc_x[k] * smf_520[k];

        t_781[k] = f_8 * slf_430[k]
                   + f_3 * pc_y[k] * smf_520[k];

        t_782[k] = f_10 * slf_420[k]
                   + f_3 * pc_z[k] * smf_520[k];

        t_783[k] = f_4 * smd0_315[k]
                   - f_5 * smd1_315[k]
                   + f_3 * pc_x[k] * smf_523[k];
    }

#pragma omp simd aligned(t_784, t_785, t_786, t_787, t_788, pc_x, pc_y, slf_432, smd0_317, \
                         smd1_317, smf_522, smf_525, smf_526, smf_527, \
                         smf_528 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_784[k] = f_8 * slf_432[k]
                   + f_3 * pc_y[k] * smf_522[k];

        t_785[k] = f_4 * smd0_317[k]
                   - f_5 * smd1_317[k]
                   + f_3 * pc_x[k] * smf_525[k];

        t_786[k] = f_3 * pc_x[k] * smf_526[k];

        t_787[k] = f_3 * pc_x[k] * smf_527[k];

        t_788[k] = f_3 * pc_x[k] * smf_528[k];
    }

#pragma omp simd aligned(t_789, t_790, t_791, pc_x, pc_y, pc_z, slf_426, slf_436, smd0_315, \
                         smd1_315, smf_526, smf_529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_789[k] = f_3 * pc_x[k] * smf_529[k];

        t_790[k] = f_8 * slf_436[k]
                   + f_1 * smd0_315[k]
                   - f_2 * smd1_315[k]
                   + f_3 * pc_y[k] * smf_526[k];

        t_791[k] = f_10 * slf_426[k]
                   + f_3 * pc_z[k] * smf_526[k];
    }

#pragma omp simd aligned(t_792, t_793, t_794, t_795, pb_y, pc_y, pc_z, slg0_660, slf_429, \
                         slf_438, slf_439, slg1_660, smd0_317, smd1_317, smf_528, \
                         smf_529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_792[k] = f_8 * slf_438[k]
                   + f_4 * smd0_317[k]
                   - f_5 * smd1_317[k]
                   + f_3 * pc_y[k] * smf_528[k];

        t_793[k] = f_8 * slf_439[k]
                   + f_3 * pc_y[k] * smf_529[k];

        t_794[k] = f_10 * slf_429[k]
                   + f_1 * smd0_317[k]
                   - f_2 * smd1_317[k]
                   + f_3 * pc_z[k] * smf_529[k];

        t_795[k] = pb_y[k] * slg0_660[k]
                   - f_6 * pc_y[k] * slg1_660[k];
    }

#pragma omp simd aligned(t_796, t_797, t_798, t_799, pc_x, pc_y, pc_z, slf_430, slf_440, \
                         slf_442, smd0_321, smd1_321, smf_530, smf_532, \
                         smf_533 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_796[k] = f_7 * slf_440[k]
                   + f_3 * pc_y[k] * smf_530[k];

        t_797[k] = f_9 * slf_430[k]
                   + f_3 * pc_z[k] * smf_530[k];

        t_798[k] = f_4 * smd0_321[k]
                   - f_5 * smd1_321[k]
                   + f_3 * pc_x[k] * smf_533[k];

        t_799[k] = f_7 * slf_442[k]
                   + f_3 * pc_y[k] * smf_532[k];
    }

#pragma omp simd aligned(t_800, t_801, t_802, t_803, t_804, pb_y, pc_x, pc_y, slg0_665, \
                         slg1_665, smf_536, smf_537, smf_538, smf_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_800[k] = pb_y[k] * slg0_665[k]
                   - f_6 * pc_y[k] * slg1_665[k];

        t_801[k] = f_3 * pc_x[k] * smf_536[k];

        t_802[k] = f_3 * pc_x[k] * smf_537[k];

        t_803[k] = f_3 * pc_x[k] * smf_538[k];

        t_804[k] = f_3 * pc_x[k] * smf_539[k];
    }

#pragma omp simd aligned(t_805, t_806, t_807, pb_y, pc_y, pc_z, slg0_670, slg0_672, slf_436, \
                         slf_446, slf_448, slg1_670, slg1_672, \
                         smf_536 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_805[k] = pb_y[k] * slg0_670[k]
                   + f_14 * slf_446[k]
                   - f_6 * pc_y[k] * slg1_670[k];

        t_806[k] = f_9 * slf_436[k]
                   + f_3 * pc_z[k] * smf_536[k];

        t_807[k] = pb_y[k] * slg0_672[k]
                   + f_8 * slf_448[k]
                   - f_6 * pc_y[k] * slg1_672[k];
    }

#pragma omp simd aligned(t_808, t_809, t_810, t_811, pb_y, pc_x, pc_y, slg0_674, slf_449, \
                         slg1_674, smd0_324, smd1_324, smf_539, \
                         smf_540 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_808[k] = f_7 * slf_449[k]
                   + f_3 * pc_y[k] * smf_539[k];

        t_809[k] = pb_y[k] * slg0_674[k]
                   - f_6 * pc_y[k] * slg1_674[k];

        t_810[k] = f_1 * smd0_324[k]
                   - f_2 * smd1_324[k]
                   + f_3 * pc_x[k] * smf_540[k];

        t_811[k] = f_3 * pc_y[k] * smf_540[k];
    }

#pragma omp simd aligned(t_812, t_813, t_814, t_815, pc_x, pc_y, pc_z, slf_440, smd0_327, \
                         smd0_329, smd1_327, smd1_329, smf_540, smf_542, smf_543, \
                         smf_545 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_812[k] = f_0 * slf_440[k]
                   + f_3 * pc_z[k] * smf_540[k];

        t_813[k] = f_4 * smd0_327[k]
                   - f_5 * smd1_327[k]
                   + f_3 * pc_x[k] * smf_543[k];

        t_814[k] = f_3 * pc_y[k] * smf_542[k];

        t_815[k] = f_4 * smd0_329[k]
                   - f_5 * smd1_329[k]
                   + f_3 * pc_x[k] * smf_545[k];
    }

#pragma omp simd aligned(t_816, t_817, t_818, t_819, t_820, t_821, pc_x, pc_y, pc_z, slf_446, \
                         smd0_327, smd1_327, smf_546, smf_547, smf_548, \
                         smf_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_816[k] = f_3 * pc_x[k] * smf_546[k];

        t_817[k] = f_3 * pc_x[k] * smf_547[k];

        t_818[k] = f_3 * pc_x[k] * smf_548[k];

        t_819[k] = f_3 * pc_x[k] * smf_549[k];

        t_820[k] = f_1 * smd0_327[k]
                   - f_2 * smd1_327[k]
                   + f_3 * pc_y[k] * smf_546[k];

        t_821[k] = f_0 * slf_446[k]
                   + f_3 * pc_z[k] * smf_546[k];
    }

#pragma omp simd aligned(t_822, t_823, t_824, pc_y, pc_z, slf_449, smd0_329, smd1_329, \
                         smf_548, smf_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_822[k] = f_4 * smd0_329[k]
                   - f_5 * smd1_329[k]
                   + f_3 * pc_y[k] * smf_548[k];

        t_823[k] = f_3 * pc_y[k] * smf_549[k];

        t_824[k] = f_0 * slf_449[k]
                   + f_1 * smd0_329[k]
                   - f_2 * smd1_329[k]
                   + f_3 * pc_z[k] * smf_549[k];
    }
}

auto
compute_prim_smg_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t slg0, const size_t slf,
                                                   const size_t slg1, const size_t smd0,
                                                   const size_t smd1, const size_t smf,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_smg_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, slg0, slf,
                                                              slg1, smd0, smd1, smf, ncols,
                                                              gamma, p, q);

    compute_prim_smg_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, slg0, slf,
                                                              slg1, smd0, smd1, smf, ncols,
                                                              gamma, p, q);

    compute_prim_smg_three_center_electron_repulsion_0_piece2(buffer, target, pb, pc, slg0, slf,
                                                              slg1, smd0, smd1, smf, ncols,
                                                              gamma, p, q);

    compute_prim_smg_three_center_electron_repulsion_0_piece3(buffer, target, pb, pc, slg0, slf,
                                                              slg1, smd0, smd1, smf, ncols,
                                                              gamma, p, q);

    compute_prim_smg_three_center_electron_repulsion_0_piece4(buffer, target, pb, pc, slg0, slf,
                                                              slg1, smd0, smd1, smf, ncols,
                                                              gamma, p, q);

    compute_prim_smg_three_center_electron_repulsion_0_piece5(buffer, target, pb, pc, slg0, slf,
                                                              slg1, smd0, smd1, smf, ncols,
                                                              gamma, p, q);

    compute_prim_smg_three_center_electron_repulsion_0_piece6(buffer, target, pb, pc, slg0, slf,
                                                              slg1, smd0, smd1, smf, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
