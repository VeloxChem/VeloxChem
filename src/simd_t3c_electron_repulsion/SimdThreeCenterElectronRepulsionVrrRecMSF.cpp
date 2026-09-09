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


#include "SimdThreeCenterElectronRepulsionVrrRecMSF.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_msf_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsf0,
                                                          const size_t lsd, const size_t lsf1,
                                                          const size_t msp0, const size_t msp1,
                                                          const size_t msd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / q;
    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 4.0 / q;
    const auto f_7 = 0.5 / gamma;
    const auto f_8 = 0.5 * p / (gamma * q);
    const auto f_9 = 3.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 3.0 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 2.5 / q;
    const auto f_14 = 2.0 / q;

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
    auto *t_126 = buffer.data(target + 126);
    auto *t_127 = buffer.data(target + 127);
    auto *t_128 = buffer.data(target + 128);
    auto *t_129 = buffer.data(target + 129);
    auto *t_130 = buffer.data(target + 130);
    auto *t_131 = buffer.data(target + 131);
    auto *t_132 = buffer.data(target + 132);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *lsf0_0 = buffer.data(lsf0 + 0);
    const auto *lsf0_6 = buffer.data(lsf0 + 6);
    const auto *lsf0_9 = buffer.data(lsf0 + 9);
    const auto *lsf0_16 = buffer.data(lsf0 + 16);
    const auto *lsf0_20 = buffer.data(lsf0 + 20);
    const auto *lsf0_29 = buffer.data(lsf0 + 29);
    const auto *lsf0_30 = buffer.data(lsf0 + 30);
    const auto *lsf0_36 = buffer.data(lsf0 + 36);
    const auto *lsf0_50 = buffer.data(lsf0 + 50);
    const auto *lsf0_59 = buffer.data(lsf0 + 59);
    const auto *lsf0_60 = buffer.data(lsf0 + 60);
    const auto *lsf0_66 = buffer.data(lsf0 + 66);
    const auto *lsf0_90 = buffer.data(lsf0 + 90);

    const auto *lsd_0 = buffer.data(lsd + 0);
    const auto *lsd_3 = buffer.data(lsd + 3);
    const auto *lsd_5 = buffer.data(lsd + 5);
    const auto *lsd_6 = buffer.data(lsd + 6);
    const auto *lsd_9 = buffer.data(lsd + 9);
    const auto *lsd_11 = buffer.data(lsd + 11);
    const auto *lsd_12 = buffer.data(lsd + 12);
    const auto *lsd_15 = buffer.data(lsd + 15);
    const auto *lsd_17 = buffer.data(lsd + 17);
    const auto *lsd_18 = buffer.data(lsd + 18);
    const auto *lsd_21 = buffer.data(lsd + 21);
    const auto *lsd_23 = buffer.data(lsd + 23);
    const auto *lsd_24 = buffer.data(lsd + 24);
    const auto *lsd_27 = buffer.data(lsd + 27);
    const auto *lsd_28 = buffer.data(lsd + 28);
    const auto *lsd_29 = buffer.data(lsd + 29);
    const auto *lsd_30 = buffer.data(lsd + 30);
    const auto *lsd_33 = buffer.data(lsd + 33);
    const auto *lsd_35 = buffer.data(lsd + 35);
    const auto *lsd_36 = buffer.data(lsd + 36);
    const auto *lsd_39 = buffer.data(lsd + 39);
    const auto *lsd_41 = buffer.data(lsd + 41);
    const auto *lsd_42 = buffer.data(lsd + 42);
    const auto *lsd_45 = buffer.data(lsd + 45);
    const auto *lsd_46 = buffer.data(lsd + 46);
    const auto *lsd_47 = buffer.data(lsd + 47);
    const auto *lsd_48 = buffer.data(lsd + 48);
    const auto *lsd_51 = buffer.data(lsd + 51);
    const auto *lsd_52 = buffer.data(lsd + 52);
    const auto *lsd_53 = buffer.data(lsd + 53);
    const auto *lsd_54 = buffer.data(lsd + 54);
    const auto *lsd_57 = buffer.data(lsd + 57);
    const auto *lsd_59 = buffer.data(lsd + 59);
    const auto *lsd_60 = buffer.data(lsd + 60);
    const auto *lsd_63 = buffer.data(lsd + 63);
    const auto *lsd_65 = buffer.data(lsd + 65);
    const auto *lsd_69 = buffer.data(lsd + 69);
    const auto *lsd_70 = buffer.data(lsd + 70);
    const auto *lsd_71 = buffer.data(lsd + 71);
    const auto *lsd_72 = buffer.data(lsd + 72);
    const auto *lsd_75 = buffer.data(lsd + 75);
    const auto *lsd_76 = buffer.data(lsd + 76);
    const auto *lsd_77 = buffer.data(lsd + 77);

    const auto *lsf1_0 = buffer.data(lsf1 + 0);
    const auto *lsf1_6 = buffer.data(lsf1 + 6);
    const auto *lsf1_9 = buffer.data(lsf1 + 9);
    const auto *lsf1_16 = buffer.data(lsf1 + 16);
    const auto *lsf1_20 = buffer.data(lsf1 + 20);
    const auto *lsf1_29 = buffer.data(lsf1 + 29);
    const auto *lsf1_30 = buffer.data(lsf1 + 30);
    const auto *lsf1_36 = buffer.data(lsf1 + 36);
    const auto *lsf1_50 = buffer.data(lsf1 + 50);
    const auto *lsf1_59 = buffer.data(lsf1 + 59);
    const auto *lsf1_60 = buffer.data(lsf1 + 60);
    const auto *lsf1_66 = buffer.data(lsf1 + 66);
    const auto *lsf1_90 = buffer.data(lsf1 + 90);

    const auto *msp0_0 = buffer.data(msp0 + 0);
    const auto *msp0_1 = buffer.data(msp0 + 1);
    const auto *msp0_2 = buffer.data(msp0 + 2);
    const auto *msp0_4 = buffer.data(msp0 + 4);
    const auto *msp0_8 = buffer.data(msp0 + 8);
    const auto *msp0_9 = buffer.data(msp0 + 9);
    const auto *msp0_10 = buffer.data(msp0 + 10);
    const auto *msp0_11 = buffer.data(msp0 + 11);
    const auto *msp0_15 = buffer.data(msp0 + 15);
    const auto *msp0_16 = buffer.data(msp0 + 16);
    const auto *msp0_17 = buffer.data(msp0 + 17);
    const auto *msp0_18 = buffer.data(msp0 + 18);
    const auto *msp0_19 = buffer.data(msp0 + 19);
    const auto *msp0_20 = buffer.data(msp0 + 20);
    const auto *msp0_23 = buffer.data(msp0 + 23);
    const auto *msp0_25 = buffer.data(msp0 + 25);
    const auto *msp0_27 = buffer.data(msp0 + 27);
    const auto *msp0_28 = buffer.data(msp0 + 28);
    const auto *msp0_29 = buffer.data(msp0 + 29);
    const auto *msp0_30 = buffer.data(msp0 + 30);
    const auto *msp0_31 = buffer.data(msp0 + 31);
    const auto *msp0_32 = buffer.data(msp0 + 32);
    const auto *msp0_35 = buffer.data(msp0 + 35);
    const auto *msp0_36 = buffer.data(msp0 + 36);
    const auto *msp0_37 = buffer.data(msp0 + 37);
    const auto *msp0_38 = buffer.data(msp0 + 38);

    const auto *msp1_0 = buffer.data(msp1 + 0);
    const auto *msp1_1 = buffer.data(msp1 + 1);
    const auto *msp1_2 = buffer.data(msp1 + 2);
    const auto *msp1_4 = buffer.data(msp1 + 4);
    const auto *msp1_8 = buffer.data(msp1 + 8);
    const auto *msp1_9 = buffer.data(msp1 + 9);
    const auto *msp1_10 = buffer.data(msp1 + 10);
    const auto *msp1_11 = buffer.data(msp1 + 11);
    const auto *msp1_15 = buffer.data(msp1 + 15);
    const auto *msp1_16 = buffer.data(msp1 + 16);
    const auto *msp1_17 = buffer.data(msp1 + 17);
    const auto *msp1_18 = buffer.data(msp1 + 18);
    const auto *msp1_19 = buffer.data(msp1 + 19);
    const auto *msp1_20 = buffer.data(msp1 + 20);
    const auto *msp1_23 = buffer.data(msp1 + 23);
    const auto *msp1_25 = buffer.data(msp1 + 25);
    const auto *msp1_27 = buffer.data(msp1 + 27);
    const auto *msp1_28 = buffer.data(msp1 + 28);
    const auto *msp1_29 = buffer.data(msp1 + 29);
    const auto *msp1_30 = buffer.data(msp1 + 30);
    const auto *msp1_31 = buffer.data(msp1 + 31);
    const auto *msp1_32 = buffer.data(msp1 + 32);
    const auto *msp1_35 = buffer.data(msp1 + 35);
    const auto *msp1_36 = buffer.data(msp1 + 36);
    const auto *msp1_37 = buffer.data(msp1 + 37);
    const auto *msp1_38 = buffer.data(msp1 + 38);

    const auto *msd_0 = buffer.data(msd + 0);
    const auto *msd_2 = buffer.data(msd + 2);
    const auto *msd_3 = buffer.data(msd + 3);
    const auto *msd_5 = buffer.data(msd + 5);
    const auto *msd_6 = buffer.data(msd + 6);
    const auto *msd_7 = buffer.data(msd + 7);
    const auto *msd_9 = buffer.data(msd + 9);
    const auto *msd_11 = buffer.data(msd + 11);
    const auto *msd_12 = buffer.data(msd + 12);
    const auto *msd_14 = buffer.data(msd + 14);
    const auto *msd_15 = buffer.data(msd + 15);
    const auto *msd_16 = buffer.data(msd + 16);
    const auto *msd_17 = buffer.data(msd + 17);
    const auto *msd_18 = buffer.data(msd + 18);
    const auto *msd_19 = buffer.data(msd + 19);
    const auto *msd_21 = buffer.data(msd + 21);
    const auto *msd_23 = buffer.data(msd + 23);
    const auto *msd_24 = buffer.data(msd + 24);
    const auto *msd_27 = buffer.data(msd + 27);
    const auto *msd_28 = buffer.data(msd + 28);
    const auto *msd_29 = buffer.data(msd + 29);
    const auto *msd_30 = buffer.data(msd + 30);
    const auto *msd_32 = buffer.data(msd + 32);
    const auto *msd_33 = buffer.data(msd + 33);
    const auto *msd_34 = buffer.data(msd + 34);
    const auto *msd_35 = buffer.data(msd + 35);
    const auto *msd_36 = buffer.data(msd + 36);
    const auto *msd_37 = buffer.data(msd + 37);
    const auto *msd_39 = buffer.data(msd + 39);
    const auto *msd_41 = buffer.data(msd + 41);
    const auto *msd_42 = buffer.data(msd + 42);
    const auto *msd_45 = buffer.data(msd + 45);
    const auto *msd_46 = buffer.data(msd + 46);
    const auto *msd_47 = buffer.data(msd + 47);
    const auto *msd_48 = buffer.data(msd + 48);
    const auto *msd_51 = buffer.data(msd + 51);
    const auto *msd_52 = buffer.data(msd + 52);
    const auto *msd_53 = buffer.data(msd + 53);
    const auto *msd_54 = buffer.data(msd + 54);
    const auto *msd_56 = buffer.data(msd + 56);
    const auto *msd_57 = buffer.data(msd + 57);
    const auto *msd_58 = buffer.data(msd + 58);
    const auto *msd_59 = buffer.data(msd + 59);
    const auto *msd_60 = buffer.data(msd + 60);
    const auto *msd_61 = buffer.data(msd + 61);
    const auto *msd_63 = buffer.data(msd + 63);
    const auto *msd_65 = buffer.data(msd + 65);
    const auto *msd_66 = buffer.data(msd + 66);
    const auto *msd_69 = buffer.data(msd + 69);
    const auto *msd_70 = buffer.data(msd + 70);
    const auto *msd_71 = buffer.data(msd + 71);
    const auto *msd_72 = buffer.data(msd + 72);
    const auto *msd_75 = buffer.data(msd + 75);
    const auto *msd_76 = buffer.data(msd + 76);
    const auto *msd_77 = buffer.data(msd + 77);
    const auto *msd_78 = buffer.data(msd + 78);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pc_x, pc_y, pc_z, lsd_0, lsd_3, msp0_0, \
                         msp1_0, msd_0, msd_2, msd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * lsd_0[k]
                 + f_1 * msp0_0[k]
                 - f_2 * msp1_0[k]
                 + f_3 * pc_x[k] * msd_0[k];

        t_1[k] = f_3 * pc_y[k] * msd_0[k];

        t_2[k] = f_3 * pc_z[k] * msd_0[k];

        t_3[k] = f_0 * lsd_3[k]
                 + f_3 * pc_x[k] * msd_3[k];

        t_4[k] = f_3 * pc_y[k] * msd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pc_x, pc_y, pc_z, lsd_5, msp0_1, msp0_2, \
                         msp1_1, msp1_2, msd_3, msd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * lsd_5[k]
                 + f_3 * pc_x[k] * msd_5[k];

        t_6[k] = f_1 * msp0_1[k]
                 - f_2 * msp1_1[k]
                 + f_3 * pc_y[k] * msd_3[k];

        t_7[k] = f_3 * pc_z[k] * msd_3[k];

        t_8[k] = f_3 * pc_y[k] * msd_5[k];

        t_9[k] = f_1 * msp0_2[k]
                 - f_2 * msp1_2[k]
                 + f_3 * pc_z[k] * msd_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_y, pc_x, pc_y, pc_z, lsf0_0, lsd_0, \
                         lsd_9, lsf1_0, msd_6, msd_7, msd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_y[k] * lsf0_0[k]
                  - f_4 * pc_y[k] * lsf1_0[k];

        t_11[k] = f_5 * lsd_0[k]
                  + f_3 * pc_y[k] * msd_6[k];

        t_12[k] = f_3 * pc_z[k] * msd_6[k];

        t_13[k] = f_6 * lsd_9[k]
                  + f_3 * pc_x[k] * msd_9[k];

        t_14[k] = f_3 * pc_z[k] * msd_7[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pc_x, pc_y, pc_z, lsd_3, lsd_5, lsd_11, \
                         msp0_4, msp1_4, msd_9, msd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_6 * lsd_11[k]
                  + f_3 * pc_x[k] * msd_11[k];

        t_16[k] = f_5 * lsd_3[k]
                  + f_1 * msp0_4[k]
                  - f_2 * msp1_4[k]
                  + f_3 * pc_y[k] * msd_9[k];

        t_17[k] = f_3 * pc_z[k] * msd_9[k];

        t_18[k] = f_5 * lsd_5[k]
                  + f_3 * pc_y[k] * msd_11[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_y, pa_z, pc_y, pc_z, lsf0_0, lsf0_9, \
                         lsd_0, lsf1_0, lsf1_9, msd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pa_y[k] * lsf0_9[k]
                  - f_4 * pc_y[k] * lsf1_9[k];

        t_20[k] = pa_z[k] * lsf0_0[k]
                  - f_4 * pc_z[k] * lsf1_0[k];

        t_21[k] = f_3 * pc_y[k] * msd_12[k];

        t_22[k] = f_5 * lsd_0[k]
                  + f_3 * pc_z[k] * msd_12[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_z, pc_x, pc_y, pc_z, lsf0_6, lsd_15, \
                         lsd_17, lsf1_6, msd_14, msd_15, msd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_6 * lsd_15[k]
                  + f_3 * pc_x[k] * msd_15[k];

        t_24[k] = f_3 * pc_y[k] * msd_14[k];

        t_25[k] = f_6 * lsd_17[k]
                  + f_3 * pc_x[k] * msd_17[k];

        t_26[k] = pa_z[k] * lsf0_6[k]
                  - f_4 * pc_z[k] * lsf1_6[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pc_x, pc_y, pc_z, lsd_5, lsd_18, msp0_8, \
                         msp0_9, msp1_8, msp1_9, msd_16, msd_17, \
                         msd_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_7 * msp0_8[k]
                  - f_8 * msp1_8[k]
                  + f_3 * pc_y[k] * msd_16[k];

        t_28[k] = f_3 * pc_y[k] * msd_17[k];

        t_29[k] = f_5 * lsd_5[k]
                  + f_1 * msp0_8[k]
                  - f_2 * msp1_8[k]
                  + f_3 * pc_z[k] * msd_17[k];

        t_30[k] = f_9 * lsd_18[k]
                  + f_1 * msp0_9[k]
                  - f_2 * msp1_9[k]
                  + f_3 * pc_x[k] * msd_18[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, pc_x, pc_y, pc_z, lsd_6, lsd_21, \
                         lsd_23, msd_18, msd_19, msd_21, msd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_10 * lsd_6[k]
                  + f_3 * pc_y[k] * msd_18[k];

        t_32[k] = f_3 * pc_z[k] * msd_18[k];

        t_33[k] = f_9 * lsd_21[k]
                  + f_3 * pc_x[k] * msd_21[k];

        t_34[k] = f_3 * pc_z[k] * msd_19[k];

        t_35[k] = f_9 * lsd_23[k]
                  + f_3 * pc_x[k] * msd_23[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pc_y, pc_z, lsd_9, lsd_11, msp0_10, msp0_11, \
                         msp1_10, msp1_11, msd_21, msd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_10 * lsd_9[k]
                  + f_1 * msp0_10[k]
                  - f_2 * msp1_10[k]
                  + f_3 * pc_y[k] * msd_21[k];

        t_37[k] = f_3 * pc_z[k] * msd_21[k];

        t_38[k] = f_10 * lsd_11[k]
                  + f_3 * pc_y[k] * msd_23[k];

        t_39[k] = f_1 * msp0_11[k]
                  - f_2 * msp1_11[k]
                  + f_3 * pc_z[k] * msd_23[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_y, pc_x, pc_y, pc_z, lsf0_20, lsd_6, \
                         lsd_12, lsd_27, lsf1_20, msd_24, msd_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pa_y[k] * lsf0_20[k]
                  - f_4 * pc_y[k] * lsf1_20[k];

        t_41[k] = f_5 * lsd_12[k]
                  + f_3 * pc_y[k] * msd_24[k];

        t_42[k] = f_5 * lsd_6[k]
                  + f_3 * pc_z[k] * msd_24[k];

        t_43[k] = f_9 * lsd_27[k]
                  + f_3 * pc_x[k] * msd_27[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_z, pc_x, pc_z, lsf0_16, lsd_9, lsd_28, \
                         lsd_29, lsf1_16, msd_27, msd_28, msd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_9 * lsd_28[k]
                  + f_3 * pc_x[k] * msd_28[k];

        t_45[k] = f_9 * lsd_29[k]
                  + f_3 * pc_x[k] * msd_29[k];

        t_46[k] = pa_z[k] * lsf0_16[k]
                  - f_4 * pc_z[k] * lsf1_16[k];

        t_47[k] = f_5 * lsd_9[k]
                  + f_3 * pc_z[k] * msd_27[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_y, pc_x, pc_y, lsf0_29, lsd_17, lsd_30, \
                         lsf1_29, msp0_15, msp1_15, msd_29, msd_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_5 * lsd_17[k]
                  + f_3 * pc_y[k] * msd_29[k];

        t_49[k] = pa_y[k] * lsf0_29[k]
                  - f_4 * pc_y[k] * lsf1_29[k];

        t_50[k] = f_9 * lsd_30[k]
                  + f_1 * msp0_15[k]
                  - f_2 * msp1_15[k]
                  + f_3 * pc_x[k] * msd_30[k];

        t_51[k] = f_3 * pc_y[k] * msd_30[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pc_x, pc_y, pc_z, lsd_12, lsd_33, lsd_35, \
                         msd_30, msd_32, msd_33, msd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_10 * lsd_12[k]
                  + f_3 * pc_z[k] * msd_30[k];

        t_53[k] = f_9 * lsd_33[k]
                  + f_3 * pc_x[k] * msd_33[k];

        t_54[k] = f_3 * pc_y[k] * msd_32[k];

        t_55[k] = f_9 * lsd_35[k]
                  + f_3 * pc_x[k] * msd_35[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pc_y, pc_z, lsd_17, msp0_16, msp0_17, \
                         msp1_16, msp1_17, msd_33, msd_34, msd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_1 * msp0_16[k]
                  - f_2 * msp1_16[k]
                  + f_3 * pc_y[k] * msd_33[k];

        t_57[k] = f_7 * msp0_17[k]
                  - f_8 * msp1_17[k]
                  + f_3 * pc_y[k] * msd_34[k];

        t_58[k] = f_3 * pc_y[k] * msd_35[k];

        t_59[k] = f_10 * lsd_17[k]
                  + f_1 * msp0_17[k]
                  - f_2 * msp1_17[k]
                  + f_3 * pc_z[k] * msd_35[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, pc_x, pc_y, pc_z, lsd_18, lsd_36, \
                         lsd_39, msp0_18, msp1_18, msd_36, msd_37, \
                         msd_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_11 * lsd_36[k]
                  + f_1 * msp0_18[k]
                  - f_2 * msp1_18[k]
                  + f_3 * pc_x[k] * msd_36[k];

        t_61[k] = f_12 * lsd_18[k]
                  + f_3 * pc_y[k] * msd_36[k];

        t_62[k] = f_3 * pc_z[k] * msd_36[k];

        t_63[k] = f_11 * lsd_39[k]
                  + f_3 * pc_x[k] * msd_39[k];

        t_64[k] = f_3 * pc_z[k] * msd_37[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pc_x, pc_y, pc_z, lsd_21, lsd_23, lsd_41, \
                         msp0_19, msp1_19, msd_39, msd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_11 * lsd_41[k]
                  + f_3 * pc_x[k] * msd_41[k];

        t_66[k] = f_12 * lsd_21[k]
                  + f_1 * msp0_19[k]
                  - f_2 * msp1_19[k]
                  + f_3 * pc_y[k] * msd_39[k];

        t_67[k] = f_3 * pc_z[k] * msd_39[k];

        t_68[k] = f_12 * lsd_23[k]
                  + f_3 * pc_y[k] * msd_41[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pa_z, pc_y, pc_z, lsf0_30, lsd_18, lsd_24, \
                         lsf1_30, msp0_20, msp1_20, msd_41, msd_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_1 * msp0_20[k]
                  - f_2 * msp1_20[k]
                  + f_3 * pc_z[k] * msd_41[k];

        t_70[k] = pa_z[k] * lsf0_30[k]
                  - f_4 * pc_z[k] * lsf1_30[k];

        t_71[k] = f_10 * lsd_24[k]
                  + f_3 * pc_y[k] * msd_42[k];

        t_72[k] = f_5 * lsd_18[k]
                  + f_3 * pc_z[k] * msd_42[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pa_z, pc_x, pc_z, lsf0_36, lsd_45, lsd_46, \
                         lsd_47, lsf1_36, msd_45, msd_46, msd_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_11 * lsd_45[k]
                  + f_3 * pc_x[k] * msd_45[k];

        t_74[k] = f_11 * lsd_46[k]
                  + f_3 * pc_x[k] * msd_46[k];

        t_75[k] = f_11 * lsd_47[k]
                  + f_3 * pc_x[k] * msd_47[k];

        t_76[k] = pa_z[k] * lsf0_36[k]
                  - f_4 * pc_z[k] * lsf1_36[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, pa_y, pc_y, pc_z, lsf0_50, lsd_21, lsd_23, \
                         lsd_29, lsf1_50, msp0_23, msp1_23, msd_45, \
                         msd_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_5 * lsd_21[k]
                  + f_3 * pc_z[k] * msd_45[k];

        t_78[k] = f_10 * lsd_29[k]
                  + f_3 * pc_y[k] * msd_47[k];

        t_79[k] = f_5 * lsd_23[k]
                  + f_1 * msp0_23[k]
                  - f_2 * msp1_23[k]
                  + f_3 * pc_z[k] * msd_47[k];

        t_80[k] = pa_y[k] * lsf0_50[k]
                  - f_4 * pc_y[k] * lsf1_50[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pc_x, pc_y, pc_z, lsd_24, lsd_30, lsd_51, \
                         lsd_52, msd_48, msd_51, msd_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_5 * lsd_30[k]
                  + f_3 * pc_y[k] * msd_48[k];

        t_82[k] = f_10 * lsd_24[k]
                  + f_3 * pc_z[k] * msd_48[k];

        t_83[k] = f_11 * lsd_51[k]
                  + f_3 * pc_x[k] * msd_51[k];

        t_84[k] = f_11 * lsd_52[k]
                  + f_3 * pc_x[k] * msd_52[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pc_x, pc_y, pc_z, lsd_27, lsd_33, lsd_35, \
                         lsd_53, msp0_25, msp1_25, msd_51, msd_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_11 * lsd_53[k]
                  + f_3 * pc_x[k] * msd_53[k];

        t_86[k] = f_5 * lsd_33[k]
                  + f_1 * msp0_25[k]
                  - f_2 * msp1_25[k]
                  + f_3 * pc_y[k] * msd_51[k];

        t_87[k] = f_10 * lsd_27[k]
                  + f_3 * pc_z[k] * msd_51[k];

        t_88[k] = f_5 * lsd_35[k]
                  + f_3 * pc_y[k] * msd_53[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, pa_y, pc_x, pc_y, pc_z, lsf0_59, lsd_30, \
                         lsd_54, lsf1_59, msp0_27, msp1_27, msd_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = pa_y[k] * lsf0_59[k]
                  - f_4 * pc_y[k] * lsf1_59[k];

        t_90[k] = f_11 * lsd_54[k]
                  + f_1 * msp0_27[k]
                  - f_2 * msp1_27[k]
                  + f_3 * pc_x[k] * msd_54[k];

        t_91[k] = f_3 * pc_y[k] * msd_54[k];

        t_92[k] = f_12 * lsd_30[k]
                  + f_3 * pc_z[k] * msd_54[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, pc_x, pc_y, lsd_57, lsd_59, msp0_28, msp1_28, \
                         msd_56, msd_57, msd_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_11 * lsd_57[k]
                  + f_3 * pc_x[k] * msd_57[k];

        t_94[k] = f_3 * pc_y[k] * msd_56[k];

        t_95[k] = f_11 * lsd_59[k]
                  + f_3 * pc_x[k] * msd_59[k];

        t_96[k] = f_1 * msp0_28[k]
                  - f_2 * msp1_28[k]
                  + f_3 * pc_y[k] * msd_57[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pc_x, pc_y, pc_z, lsd_35, lsd_60, msp0_29, \
                         msp0_30, msp1_29, msp1_30, msd_58, msd_59, \
                         msd_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_7 * msp0_29[k]
                  - f_8 * msp1_29[k]
                  + f_3 * pc_y[k] * msd_58[k];

        t_98[k] = f_3 * pc_y[k] * msd_59[k];

        t_99[k] = f_12 * lsd_35[k]
                  + f_1 * msp0_29[k]
                  - f_2 * msp1_29[k]
                  + f_3 * pc_z[k] * msd_59[k];

        t_100[k] = f_13 * lsd_60[k]
                   + f_1 * msp0_30[k]
                   - f_2 * msp1_30[k]
                   + f_3 * pc_x[k] * msd_60[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, t_105, pc_x, pc_y, pc_z, lsd_36, lsd_63, \
                         lsd_65, msd_60, msd_61, msd_63, msd_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_14 * lsd_36[k]
                   + f_3 * pc_y[k] * msd_60[k];

        t_102[k] = f_3 * pc_z[k] * msd_60[k];

        t_103[k] = f_13 * lsd_63[k]
                   + f_3 * pc_x[k] * msd_63[k];

        t_104[k] = f_3 * pc_z[k] * msd_61[k];

        t_105[k] = f_13 * lsd_65[k]
                   + f_3 * pc_x[k] * msd_65[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, pc_y, pc_z, lsd_39, lsd_41, msp0_31, \
                         msp0_32, msp1_31, msp1_32, msd_63, msd_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_14 * lsd_39[k]
                   + f_1 * msp0_31[k]
                   - f_2 * msp1_31[k]
                   + f_3 * pc_y[k] * msd_63[k];

        t_107[k] = f_3 * pc_z[k] * msd_63[k];

        t_108[k] = f_14 * lsd_41[k]
                   + f_3 * pc_y[k] * msd_65[k];

        t_109[k] = f_1 * msp0_32[k]
                   - f_2 * msp1_32[k]
                   + f_3 * pc_z[k] * msd_65[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, pa_z, pc_x, pc_y, pc_z, lsf0_60, lsd_36, \
                         lsd_42, lsd_69, lsf1_60, msd_66, msd_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = pa_z[k] * lsf0_60[k]
                   - f_4 * pc_z[k] * lsf1_60[k];

        t_111[k] = f_12 * lsd_42[k]
                   + f_3 * pc_y[k] * msd_66[k];

        t_112[k] = f_5 * lsd_36[k]
                   + f_3 * pc_z[k] * msd_66[k];

        t_113[k] = f_13 * lsd_69[k]
                   + f_3 * pc_x[k] * msd_69[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pa_z, pc_x, pc_z, lsf0_66, lsd_39, \
                         lsd_70, lsd_71, lsf1_66, msd_69, msd_70, \
                         msd_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_13 * lsd_70[k]
                   + f_3 * pc_x[k] * msd_70[k];

        t_115[k] = f_13 * lsd_71[k]
                   + f_3 * pc_x[k] * msd_71[k];

        t_116[k] = pa_z[k] * lsf0_66[k]
                   - f_4 * pc_z[k] * lsf1_66[k];

        t_117[k] = f_5 * lsd_39[k]
                   + f_3 * pc_z[k] * msd_69[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, pc_x, pc_y, pc_z, lsd_41, lsd_47, lsd_72, \
                         msp0_35, msp0_36, msp1_35, msp1_36, msd_71, \
                         msd_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_12 * lsd_47[k]
                   + f_3 * pc_y[k] * msd_71[k];

        t_119[k] = f_5 * lsd_41[k]
                   + f_1 * msp0_35[k]
                   - f_2 * msp1_35[k]
                   + f_3 * pc_z[k] * msd_71[k];

        t_120[k] = f_13 * lsd_72[k]
                   + f_1 * msp0_36[k]
                   - f_2 * msp1_36[k]
                   + f_3 * pc_x[k] * msd_72[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pc_x, pc_y, pc_z, lsd_42, lsd_48, lsd_75, \
                         lsd_76, msd_72, msd_75, msd_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_10 * lsd_48[k]
                   + f_3 * pc_y[k] * msd_72[k];

        t_122[k] = f_10 * lsd_42[k]
                   + f_3 * pc_z[k] * msd_72[k];

        t_123[k] = f_13 * lsd_75[k]
                   + f_3 * pc_x[k] * msd_75[k];

        t_124[k] = f_13 * lsd_76[k]
                   + f_3 * pc_x[k] * msd_76[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pc_x, pc_y, pc_z, lsd_45, lsd_51, lsd_53, \
                         lsd_77, msp0_37, msp1_37, msd_75, msd_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_13 * lsd_77[k]
                   + f_3 * pc_x[k] * msd_77[k];

        t_126[k] = f_10 * lsd_51[k]
                   + f_1 * msp0_37[k]
                   - f_2 * msp1_37[k]
                   + f_3 * pc_y[k] * msd_75[k];

        t_127[k] = f_10 * lsd_45[k]
                   + f_3 * pc_z[k] * msd_75[k];

        t_128[k] = f_10 * lsd_53[k]
                   + f_3 * pc_y[k] * msd_77[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, pa_y, pc_y, pc_z, lsf0_90, lsd_47, \
                         lsd_48, lsd_54, lsf1_90, msp0_38, msp1_38, msd_77, \
                         msd_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_10 * lsd_47[k]
                   + f_1 * msp0_38[k]
                   - f_2 * msp1_38[k]
                   + f_3 * pc_z[k] * msd_77[k];

        t_130[k] = pa_y[k] * lsf0_90[k]
                   - f_4 * pc_y[k] * lsf1_90[k];

        t_131[k] = f_5 * lsd_54[k]
                   + f_3 * pc_y[k] * msd_78[k];

        t_132[k] = f_12 * lsd_48[k]
                   + f_3 * pc_z[k] * msd_78[k];
    }
}

static auto
compute_prim_msf_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsf0,
                                                          const size_t lsd, const size_t lsf1,
                                                          const size_t msp0, const size_t msp1,
                                                          const size_t msd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_7 = 0.5 / gamma;
    const auto f_8 = 0.5 * p / (gamma * q);
    const auto f_10 = 1.0 / q;
    const auto f_11 = 3.0 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 2.5 / q;
    const auto f_14 = 2.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *lsf0_99 = buffer.data(lsf0 + 99);
    const auto *lsf0_100 = buffer.data(lsf0 + 100);
    const auto *lsf0_106 = buffer.data(lsf0 + 106);
    const auto *lsf0_140 = buffer.data(lsf0 + 140);
    const auto *lsf0_149 = buffer.data(lsf0 + 149);
    const auto *lsf0_150 = buffer.data(lsf0 + 150);
    const auto *lsf0_156 = buffer.data(lsf0 + 156);

    const auto *lsd_51 = buffer.data(lsd + 51);
    const auto *lsd_54 = buffer.data(lsd + 54);
    const auto *lsd_57 = buffer.data(lsd + 57);
    const auto *lsd_59 = buffer.data(lsd + 59);
    const auto *lsd_60 = buffer.data(lsd + 60);
    const auto *lsd_63 = buffer.data(lsd + 63);
    const auto *lsd_65 = buffer.data(lsd + 65);
    const auto *lsd_66 = buffer.data(lsd + 66);
    const auto *lsd_69 = buffer.data(lsd + 69);
    const auto *lsd_71 = buffer.data(lsd + 71);
    const auto *lsd_72 = buffer.data(lsd + 72);
    const auto *lsd_75 = buffer.data(lsd + 75);
    const auto *lsd_77 = buffer.data(lsd + 77);
    const auto *lsd_78 = buffer.data(lsd + 78);
    const auto *lsd_81 = buffer.data(lsd + 81);
    const auto *lsd_82 = buffer.data(lsd + 82);
    const auto *lsd_83 = buffer.data(lsd + 83);
    const auto *lsd_84 = buffer.data(lsd + 84);
    const auto *lsd_87 = buffer.data(lsd + 87);
    const auto *lsd_89 = buffer.data(lsd + 89);
    const auto *lsd_90 = buffer.data(lsd + 90);
    const auto *lsd_93 = buffer.data(lsd + 93);
    const auto *lsd_95 = buffer.data(lsd + 95);
    const auto *lsd_96 = buffer.data(lsd + 96);
    const auto *lsd_99 = buffer.data(lsd + 99);
    const auto *lsd_100 = buffer.data(lsd + 100);
    const auto *lsd_101 = buffer.data(lsd + 101);
    const auto *lsd_102 = buffer.data(lsd + 102);
    const auto *lsd_105 = buffer.data(lsd + 105);
    const auto *lsd_106 = buffer.data(lsd + 106);
    const auto *lsd_107 = buffer.data(lsd + 107);
    const auto *lsd_108 = buffer.data(lsd + 108);
    const auto *lsd_111 = buffer.data(lsd + 111);
    const auto *lsd_112 = buffer.data(lsd + 112);
    const auto *lsd_113 = buffer.data(lsd + 113);
    const auto *lsd_114 = buffer.data(lsd + 114);
    const auto *lsd_117 = buffer.data(lsd + 117);
    const auto *lsd_118 = buffer.data(lsd + 118);
    const auto *lsd_119 = buffer.data(lsd + 119);
    const auto *lsd_120 = buffer.data(lsd + 120);
    const auto *lsd_123 = buffer.data(lsd + 123);
    const auto *lsd_125 = buffer.data(lsd + 125);
    const auto *lsd_126 = buffer.data(lsd + 126);
    const auto *lsd_129 = buffer.data(lsd + 129);
    const auto *lsd_131 = buffer.data(lsd + 131);
    const auto *lsd_135 = buffer.data(lsd + 135);
    const auto *lsd_136 = buffer.data(lsd + 136);
    const auto *lsd_137 = buffer.data(lsd + 137);
    const auto *lsd_138 = buffer.data(lsd + 138);
    const auto *lsd_141 = buffer.data(lsd + 141);
    const auto *lsd_142 = buffer.data(lsd + 142);
    const auto *lsd_143 = buffer.data(lsd + 143);
    const auto *lsd_144 = buffer.data(lsd + 144);
    const auto *lsd_147 = buffer.data(lsd + 147);
    const auto *lsd_148 = buffer.data(lsd + 148);
    const auto *lsd_149 = buffer.data(lsd + 149);
    const auto *lsd_150 = buffer.data(lsd + 150);
    const auto *lsd_153 = buffer.data(lsd + 153);
    const auto *lsd_154 = buffer.data(lsd + 154);
    const auto *lsd_155 = buffer.data(lsd + 155);

    const auto *lsf1_99 = buffer.data(lsf1 + 99);
    const auto *lsf1_100 = buffer.data(lsf1 + 100);
    const auto *lsf1_106 = buffer.data(lsf1 + 106);
    const auto *lsf1_140 = buffer.data(lsf1 + 140);
    const auto *lsf1_149 = buffer.data(lsf1 + 149);
    const auto *lsf1_150 = buffer.data(lsf1 + 150);
    const auto *lsf1_156 = buffer.data(lsf1 + 156);

    const auto *msp0_40 = buffer.data(msp0 + 40);
    const auto *msp0_42 = buffer.data(msp0 + 42);
    const auto *msp0_43 = buffer.data(msp0 + 43);
    const auto *msp0_44 = buffer.data(msp0 + 44);
    const auto *msp0_45 = buffer.data(msp0 + 45);
    const auto *msp0_46 = buffer.data(msp0 + 46);
    const auto *msp0_47 = buffer.data(msp0 + 47);
    const auto *msp0_50 = buffer.data(msp0 + 50);
    const auto *msp0_51 = buffer.data(msp0 + 51);
    const auto *msp0_52 = buffer.data(msp0 + 52);
    const auto *msp0_53 = buffer.data(msp0 + 53);
    const auto *msp0_54 = buffer.data(msp0 + 54);
    const auto *msp0_55 = buffer.data(msp0 + 55);
    const auto *msp0_56 = buffer.data(msp0 + 56);
    const auto *msp0_58 = buffer.data(msp0 + 58);
    const auto *msp0_60 = buffer.data(msp0 + 60);
    const auto *msp0_61 = buffer.data(msp0 + 61);
    const auto *msp0_62 = buffer.data(msp0 + 62);
    const auto *msp0_63 = buffer.data(msp0 + 63);
    const auto *msp0_64 = buffer.data(msp0 + 64);
    const auto *msp0_65 = buffer.data(msp0 + 65);
    const auto *msp0_68 = buffer.data(msp0 + 68);
    const auto *msp0_69 = buffer.data(msp0 + 69);
    const auto *msp0_70 = buffer.data(msp0 + 70);
    const auto *msp0_71 = buffer.data(msp0 + 71);
    const auto *msp0_72 = buffer.data(msp0 + 72);
    const auto *msp0_73 = buffer.data(msp0 + 73);
    const auto *msp0_74 = buffer.data(msp0 + 74);
    const auto *msp0_75 = buffer.data(msp0 + 75);

    const auto *msp1_40 = buffer.data(msp1 + 40);
    const auto *msp1_42 = buffer.data(msp1 + 42);
    const auto *msp1_43 = buffer.data(msp1 + 43);
    const auto *msp1_44 = buffer.data(msp1 + 44);
    const auto *msp1_45 = buffer.data(msp1 + 45);
    const auto *msp1_46 = buffer.data(msp1 + 46);
    const auto *msp1_47 = buffer.data(msp1 + 47);
    const auto *msp1_50 = buffer.data(msp1 + 50);
    const auto *msp1_51 = buffer.data(msp1 + 51);
    const auto *msp1_52 = buffer.data(msp1 + 52);
    const auto *msp1_53 = buffer.data(msp1 + 53);
    const auto *msp1_54 = buffer.data(msp1 + 54);
    const auto *msp1_55 = buffer.data(msp1 + 55);
    const auto *msp1_56 = buffer.data(msp1 + 56);
    const auto *msp1_58 = buffer.data(msp1 + 58);
    const auto *msp1_60 = buffer.data(msp1 + 60);
    const auto *msp1_61 = buffer.data(msp1 + 61);
    const auto *msp1_62 = buffer.data(msp1 + 62);
    const auto *msp1_63 = buffer.data(msp1 + 63);
    const auto *msp1_64 = buffer.data(msp1 + 64);
    const auto *msp1_65 = buffer.data(msp1 + 65);
    const auto *msp1_68 = buffer.data(msp1 + 68);
    const auto *msp1_69 = buffer.data(msp1 + 69);
    const auto *msp1_70 = buffer.data(msp1 + 70);
    const auto *msp1_71 = buffer.data(msp1 + 71);
    const auto *msp1_72 = buffer.data(msp1 + 72);
    const auto *msp1_73 = buffer.data(msp1 + 73);
    const auto *msp1_74 = buffer.data(msp1 + 74);
    const auto *msp1_75 = buffer.data(msp1 + 75);

    const auto *msd_81 = buffer.data(msd + 81);
    const auto *msd_82 = buffer.data(msd + 82);
    const auto *msd_83 = buffer.data(msd + 83);
    const auto *msd_84 = buffer.data(msd + 84);
    const auto *msd_86 = buffer.data(msd + 86);
    const auto *msd_87 = buffer.data(msd + 87);
    const auto *msd_88 = buffer.data(msd + 88);
    const auto *msd_89 = buffer.data(msd + 89);
    const auto *msd_90 = buffer.data(msd + 90);
    const auto *msd_91 = buffer.data(msd + 91);
    const auto *msd_93 = buffer.data(msd + 93);
    const auto *msd_95 = buffer.data(msd + 95);
    const auto *msd_96 = buffer.data(msd + 96);
    const auto *msd_99 = buffer.data(msd + 99);
    const auto *msd_100 = buffer.data(msd + 100);
    const auto *msd_101 = buffer.data(msd + 101);
    const auto *msd_102 = buffer.data(msd + 102);
    const auto *msd_105 = buffer.data(msd + 105);
    const auto *msd_106 = buffer.data(msd + 106);
    const auto *msd_107 = buffer.data(msd + 107);
    const auto *msd_108 = buffer.data(msd + 108);
    const auto *msd_111 = buffer.data(msd + 111);
    const auto *msd_112 = buffer.data(msd + 112);
    const auto *msd_113 = buffer.data(msd + 113);
    const auto *msd_114 = buffer.data(msd + 114);
    const auto *msd_117 = buffer.data(msd + 117);
    const auto *msd_118 = buffer.data(msd + 118);
    const auto *msd_119 = buffer.data(msd + 119);
    const auto *msd_120 = buffer.data(msd + 120);
    const auto *msd_122 = buffer.data(msd + 122);
    const auto *msd_123 = buffer.data(msd + 123);
    const auto *msd_124 = buffer.data(msd + 124);
    const auto *msd_125 = buffer.data(msd + 125);
    const auto *msd_126 = buffer.data(msd + 126);
    const auto *msd_127 = buffer.data(msd + 127);
    const auto *msd_129 = buffer.data(msd + 129);
    const auto *msd_131 = buffer.data(msd + 131);
    const auto *msd_132 = buffer.data(msd + 132);
    const auto *msd_135 = buffer.data(msd + 135);
    const auto *msd_136 = buffer.data(msd + 136);
    const auto *msd_137 = buffer.data(msd + 137);
    const auto *msd_138 = buffer.data(msd + 138);
    const auto *msd_141 = buffer.data(msd + 141);
    const auto *msd_142 = buffer.data(msd + 142);
    const auto *msd_143 = buffer.data(msd + 143);
    const auto *msd_144 = buffer.data(msd + 144);
    const auto *msd_147 = buffer.data(msd + 147);
    const auto *msd_148 = buffer.data(msd + 148);
    const auto *msd_149 = buffer.data(msd + 149);
    const auto *msd_150 = buffer.data(msd + 150);
    const auto *msd_153 = buffer.data(msd + 153);
    const auto *msd_154 = buffer.data(msd + 154);
    const auto *msd_155 = buffer.data(msd + 155);

#pragma omp simd aligned(t_133, t_134, t_135, t_136, pc_x, pc_y, lsd_57, lsd_81, lsd_82, \
                         lsd_83, msp0_40, msp1_40, msd_81, msd_82, \
                         msd_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_13 * lsd_81[k]
                   + f_3 * pc_x[k] * msd_81[k];

        t_134[k] = f_13 * lsd_82[k]
                   + f_3 * pc_x[k] * msd_82[k];

        t_135[k] = f_13 * lsd_83[k]
                   + f_3 * pc_x[k] * msd_83[k];

        t_136[k] = f_5 * lsd_57[k]
                   + f_1 * msp0_40[k]
                   - f_2 * msp1_40[k]
                   + f_3 * pc_y[k] * msd_81[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, pa_y, pc_y, pc_z, lsf0_99, lsd_51, lsd_59, \
                         lsf1_99, msd_81, msd_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_12 * lsd_51[k]
                   + f_3 * pc_z[k] * msd_81[k];

        t_138[k] = f_5 * lsd_59[k]
                   + f_3 * pc_y[k] * msd_83[k];

        t_139[k] = pa_y[k] * lsf0_99[k]
                   - f_4 * pc_y[k] * lsf1_99[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, pc_x, pc_y, pc_z, lsd_54, lsd_84, \
                         lsd_87, msp0_42, msp1_42, msd_84, msd_86, \
                         msd_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_13 * lsd_84[k]
                   + f_1 * msp0_42[k]
                   - f_2 * msp1_42[k]
                   + f_3 * pc_x[k] * msd_84[k];

        t_141[k] = f_3 * pc_y[k] * msd_84[k];

        t_142[k] = f_14 * lsd_54[k]
                   + f_3 * pc_z[k] * msd_84[k];

        t_143[k] = f_13 * lsd_87[k]
                   + f_3 * pc_x[k] * msd_87[k];

        t_144[k] = f_3 * pc_y[k] * msd_86[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, pc_x, pc_y, lsd_89, msp0_43, msp0_44, \
                         msp1_43, msp1_44, msd_87, msd_88, msd_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_13 * lsd_89[k]
                   + f_3 * pc_x[k] * msd_89[k];

        t_146[k] = f_1 * msp0_43[k]
                   - f_2 * msp1_43[k]
                   + f_3 * pc_y[k] * msd_87[k];

        t_147[k] = f_7 * msp0_44[k]
                   - f_8 * msp1_44[k]
                   + f_3 * pc_y[k] * msd_88[k];

        t_148[k] = f_3 * pc_y[k] * msd_89[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, pc_x, pc_y, pc_z, lsd_59, lsd_60, lsd_90, \
                         msp0_44, msp0_45, msp1_44, msp1_45, msd_89, \
                         msd_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_14 * lsd_59[k]
                   + f_1 * msp0_44[k]
                   - f_2 * msp1_44[k]
                   + f_3 * pc_z[k] * msd_89[k];

        t_150[k] = f_14 * lsd_90[k]
                   + f_1 * msp0_45[k]
                   - f_2 * msp1_45[k]
                   + f_3 * pc_x[k] * msd_90[k];

        t_151[k] = f_13 * lsd_60[k]
                   + f_3 * pc_y[k] * msd_90[k];

        t_152[k] = f_3 * pc_z[k] * msd_90[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, t_156, t_157, pc_x, pc_y, pc_z, lsd_63, lsd_93, \
                         lsd_95, msp0_46, msp1_46, msd_91, msd_93, \
                         msd_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_14 * lsd_93[k]
                   + f_3 * pc_x[k] * msd_93[k];

        t_154[k] = f_3 * pc_z[k] * msd_91[k];

        t_155[k] = f_14 * lsd_95[k]
                   + f_3 * pc_x[k] * msd_95[k];

        t_156[k] = f_13 * lsd_63[k]
                   + f_1 * msp0_46[k]
                   - f_2 * msp1_46[k]
                   + f_3 * pc_y[k] * msd_93[k];

        t_157[k] = f_3 * pc_z[k] * msd_93[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, pa_z, pc_y, pc_z, lsf0_100, lsd_65, \
                         lsd_66, lsf1_100, msp0_47, msp1_47, msd_95, \
                         msd_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_13 * lsd_65[k]
                   + f_3 * pc_y[k] * msd_95[k];

        t_159[k] = f_1 * msp0_47[k]
                   - f_2 * msp1_47[k]
                   + f_3 * pc_z[k] * msd_95[k];

        t_160[k] = pa_z[k] * lsf0_100[k]
                   - f_4 * pc_z[k] * lsf1_100[k];

        t_161[k] = f_14 * lsd_66[k]
                   + f_3 * pc_y[k] * msd_96[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, pc_x, pc_z, lsd_60, lsd_99, lsd_100, \
                         lsd_101, msd_96, msd_99, msd_100, msd_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_5 * lsd_60[k]
                   + f_3 * pc_z[k] * msd_96[k];

        t_163[k] = f_14 * lsd_99[k]
                   + f_3 * pc_x[k] * msd_99[k];

        t_164[k] = f_14 * lsd_100[k]
                   + f_3 * pc_x[k] * msd_100[k];

        t_165[k] = f_14 * lsd_101[k]
                   + f_3 * pc_x[k] * msd_101[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, pa_z, pc_y, pc_z, lsf0_106, lsd_63, \
                         lsd_65, lsd_71, lsf1_106, msp0_50, msp1_50, msd_99, \
                         msd_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = pa_z[k] * lsf0_106[k]
                   - f_4 * pc_z[k] * lsf1_106[k];

        t_167[k] = f_5 * lsd_63[k]
                   + f_3 * pc_z[k] * msd_99[k];

        t_168[k] = f_14 * lsd_71[k]
                   + f_3 * pc_y[k] * msd_101[k];

        t_169[k] = f_5 * lsd_65[k]
                   + f_1 * msp0_50[k]
                   - f_2 * msp1_50[k]
                   + f_3 * pc_z[k] * msd_101[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, pc_x, pc_y, pc_z, lsd_66, lsd_72, \
                         lsd_102, lsd_105, msp0_51, msp1_51, msd_102, \
                         msd_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_14 * lsd_102[k]
                   + f_1 * msp0_51[k]
                   - f_2 * msp1_51[k]
                   + f_3 * pc_x[k] * msd_102[k];

        t_171[k] = f_12 * lsd_72[k]
                   + f_3 * pc_y[k] * msd_102[k];

        t_172[k] = f_10 * lsd_66[k]
                   + f_3 * pc_z[k] * msd_102[k];

        t_173[k] = f_14 * lsd_105[k]
                   + f_3 * pc_x[k] * msd_105[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, pc_x, pc_y, pc_z, lsd_69, lsd_75, \
                         lsd_106, lsd_107, msp0_52, msp1_52, msd_105, msd_106, \
                         msd_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_14 * lsd_106[k]
                   + f_3 * pc_x[k] * msd_106[k];

        t_175[k] = f_14 * lsd_107[k]
                   + f_3 * pc_x[k] * msd_107[k];

        t_176[k] = f_12 * lsd_75[k]
                   + f_1 * msp0_52[k]
                   - f_2 * msp1_52[k]
                   + f_3 * pc_y[k] * msd_105[k];

        t_177[k] = f_10 * lsd_69[k]
                   + f_3 * pc_z[k] * msd_105[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, pc_x, pc_y, pc_z, lsd_71, lsd_77, lsd_108, \
                         msp0_53, msp0_54, msp1_53, msp1_54, msd_107, \
                         msd_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_12 * lsd_77[k]
                   + f_3 * pc_y[k] * msd_107[k];

        t_179[k] = f_10 * lsd_71[k]
                   + f_1 * msp0_53[k]
                   - f_2 * msp1_53[k]
                   + f_3 * pc_z[k] * msd_107[k];

        t_180[k] = f_14 * lsd_108[k]
                   + f_1 * msp0_54[k]
                   - f_2 * msp1_54[k]
                   + f_3 * pc_x[k] * msd_108[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pc_x, pc_y, pc_z, lsd_72, lsd_78, \
                         lsd_111, lsd_112, msd_108, msd_111, msd_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_10 * lsd_78[k]
                   + f_3 * pc_y[k] * msd_108[k];

        t_182[k] = f_12 * lsd_72[k]
                   + f_3 * pc_z[k] * msd_108[k];

        t_183[k] = f_14 * lsd_111[k]
                   + f_3 * pc_x[k] * msd_111[k];

        t_184[k] = f_14 * lsd_112[k]
                   + f_3 * pc_x[k] * msd_112[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pc_x, pc_y, pc_z, lsd_75, lsd_81, lsd_83, \
                         lsd_113, msp0_55, msp1_55, msd_111, msd_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_14 * lsd_113[k]
                   + f_3 * pc_x[k] * msd_113[k];

        t_186[k] = f_10 * lsd_81[k]
                   + f_1 * msp0_55[k]
                   - f_2 * msp1_55[k]
                   + f_3 * pc_y[k] * msd_111[k];

        t_187[k] = f_12 * lsd_75[k]
                   + f_3 * pc_z[k] * msd_111[k];

        t_188[k] = f_10 * lsd_83[k]
                   + f_3 * pc_y[k] * msd_113[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, pa_y, pc_y, pc_z, lsf0_140, lsd_77, \
                         lsd_78, lsd_84, lsf1_140, msp0_56, msp1_56, msd_113, \
                         msd_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_12 * lsd_77[k]
                   + f_1 * msp0_56[k]
                   - f_2 * msp1_56[k]
                   + f_3 * pc_z[k] * msd_113[k];

        t_190[k] = pa_y[k] * lsf0_140[k]
                   - f_4 * pc_y[k] * lsf1_140[k];

        t_191[k] = f_5 * lsd_84[k]
                   + f_3 * pc_y[k] * msd_114[k];

        t_192[k] = f_14 * lsd_78[k]
                   + f_3 * pc_z[k] * msd_114[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, pc_x, pc_y, lsd_87, lsd_117, lsd_118, \
                         lsd_119, msp0_58, msp1_58, msd_117, msd_118, \
                         msd_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_14 * lsd_117[k]
                   + f_3 * pc_x[k] * msd_117[k];

        t_194[k] = f_14 * lsd_118[k]
                   + f_3 * pc_x[k] * msd_118[k];

        t_195[k] = f_14 * lsd_119[k]
                   + f_3 * pc_x[k] * msd_119[k];

        t_196[k] = f_5 * lsd_87[k]
                   + f_1 * msp0_58[k]
                   - f_2 * msp1_58[k]
                   + f_3 * pc_y[k] * msd_117[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, pa_y, pc_y, pc_z, lsf0_149, lsd_81, lsd_89, \
                         lsf1_149, msd_117, msd_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_14 * lsd_81[k]
                   + f_3 * pc_z[k] * msd_117[k];

        t_198[k] = f_5 * lsd_89[k]
                   + f_3 * pc_y[k] * msd_119[k];

        t_199[k] = pa_y[k] * lsf0_149[k]
                   - f_4 * pc_y[k] * lsf1_149[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, pc_x, pc_y, pc_z, lsd_84, lsd_120, \
                         lsd_123, msp0_60, msp1_60, msd_120, msd_122, \
                         msd_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_14 * lsd_120[k]
                   + f_1 * msp0_60[k]
                   - f_2 * msp1_60[k]
                   + f_3 * pc_x[k] * msd_120[k];

        t_201[k] = f_3 * pc_y[k] * msd_120[k];

        t_202[k] = f_13 * lsd_84[k]
                   + f_3 * pc_z[k] * msd_120[k];

        t_203[k] = f_14 * lsd_123[k]
                   + f_3 * pc_x[k] * msd_123[k];

        t_204[k] = f_3 * pc_y[k] * msd_122[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, pc_x, pc_y, lsd_125, msp0_61, msp0_62, \
                         msp1_61, msp1_62, msd_123, msd_124, msd_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = f_14 * lsd_125[k]
                   + f_3 * pc_x[k] * msd_125[k];

        t_206[k] = f_1 * msp0_61[k]
                   - f_2 * msp1_61[k]
                   + f_3 * pc_y[k] * msd_123[k];

        t_207[k] = f_7 * msp0_62[k]
                   - f_8 * msp1_62[k]
                   + f_3 * pc_y[k] * msd_124[k];

        t_208[k] = f_3 * pc_y[k] * msd_125[k];
    }

#pragma omp simd aligned(t_209, t_210, t_211, t_212, pc_x, pc_y, pc_z, lsd_89, lsd_90, \
                         lsd_126, msp0_62, msp0_63, msp1_62, msp1_63, msd_125, \
                         msd_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_209[k] = f_13 * lsd_89[k]
                   + f_1 * msp0_62[k]
                   - f_2 * msp1_62[k]
                   + f_3 * pc_z[k] * msd_125[k];

        t_210[k] = f_12 * lsd_126[k]
                   + f_1 * msp0_63[k]
                   - f_2 * msp1_63[k]
                   + f_3 * pc_x[k] * msd_126[k];

        t_211[k] = f_11 * lsd_90[k]
                   + f_3 * pc_y[k] * msd_126[k];

        t_212[k] = f_3 * pc_z[k] * msd_126[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, t_217, pc_x, pc_y, pc_z, lsd_93, lsd_129, \
                         lsd_131, msp0_64, msp1_64, msd_127, msd_129, \
                         msd_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_12 * lsd_129[k]
                   + f_3 * pc_x[k] * msd_129[k];

        t_214[k] = f_3 * pc_z[k] * msd_127[k];

        t_215[k] = f_12 * lsd_131[k]
                   + f_3 * pc_x[k] * msd_131[k];

        t_216[k] = f_11 * lsd_93[k]
                   + f_1 * msp0_64[k]
                   - f_2 * msp1_64[k]
                   + f_3 * pc_y[k] * msd_129[k];

        t_217[k] = f_3 * pc_z[k] * msd_129[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, t_221, pa_z, pc_y, pc_z, lsf0_150, lsd_95, \
                         lsd_96, lsf1_150, msp0_65, msp1_65, msd_131, \
                         msd_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_11 * lsd_95[k]
                   + f_3 * pc_y[k] * msd_131[k];

        t_219[k] = f_1 * msp0_65[k]
                   - f_2 * msp1_65[k]
                   + f_3 * pc_z[k] * msd_131[k];

        t_220[k] = pa_z[k] * lsf0_150[k]
                   - f_4 * pc_z[k] * lsf1_150[k];

        t_221[k] = f_13 * lsd_96[k]
                   + f_3 * pc_y[k] * msd_132[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, pc_x, pc_z, lsd_90, lsd_135, lsd_136, \
                         lsd_137, msd_132, msd_135, msd_136, msd_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_5 * lsd_90[k]
                   + f_3 * pc_z[k] * msd_132[k];

        t_223[k] = f_12 * lsd_135[k]
                   + f_3 * pc_x[k] * msd_135[k];

        t_224[k] = f_12 * lsd_136[k]
                   + f_3 * pc_x[k] * msd_136[k];

        t_225[k] = f_12 * lsd_137[k]
                   + f_3 * pc_x[k] * msd_137[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, t_229, pa_z, pc_y, pc_z, lsf0_156, lsd_93, \
                         lsd_95, lsd_101, lsf1_156, msp0_68, msp1_68, msd_135, \
                         msd_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = pa_z[k] * lsf0_156[k]
                   - f_4 * pc_z[k] * lsf1_156[k];

        t_227[k] = f_5 * lsd_93[k]
                   + f_3 * pc_z[k] * msd_135[k];

        t_228[k] = f_13 * lsd_101[k]
                   + f_3 * pc_y[k] * msd_137[k];

        t_229[k] = f_5 * lsd_95[k]
                   + f_1 * msp0_68[k]
                   - f_2 * msp1_68[k]
                   + f_3 * pc_z[k] * msd_137[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, pc_x, pc_y, pc_z, lsd_96, lsd_102, \
                         lsd_138, lsd_141, msp0_69, msp1_69, msd_138, \
                         msd_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_12 * lsd_138[k]
                   + f_1 * msp0_69[k]
                   - f_2 * msp1_69[k]
                   + f_3 * pc_x[k] * msd_138[k];

        t_231[k] = f_14 * lsd_102[k]
                   + f_3 * pc_y[k] * msd_138[k];

        t_232[k] = f_10 * lsd_96[k]
                   + f_3 * pc_z[k] * msd_138[k];

        t_233[k] = f_12 * lsd_141[k]
                   + f_3 * pc_x[k] * msd_141[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, pc_x, pc_y, pc_z, lsd_99, lsd_105, \
                         lsd_142, lsd_143, msp0_70, msp1_70, msd_141, msd_142, \
                         msd_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_12 * lsd_142[k]
                   + f_3 * pc_x[k] * msd_142[k];

        t_235[k] = f_12 * lsd_143[k]
                   + f_3 * pc_x[k] * msd_143[k];

        t_236[k] = f_14 * lsd_105[k]
                   + f_1 * msp0_70[k]
                   - f_2 * msp1_70[k]
                   + f_3 * pc_y[k] * msd_141[k];

        t_237[k] = f_10 * lsd_99[k]
                   + f_3 * pc_z[k] * msd_141[k];
    }

#pragma omp simd aligned(t_238, t_239, t_240, pc_x, pc_y, pc_z, lsd_101, lsd_107, lsd_144, \
                         msp0_71, msp0_72, msp1_71, msp1_72, msd_143, \
                         msd_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_238[k] = f_14 * lsd_107[k]
                   + f_3 * pc_y[k] * msd_143[k];

        t_239[k] = f_10 * lsd_101[k]
                   + f_1 * msp0_71[k]
                   - f_2 * msp1_71[k]
                   + f_3 * pc_z[k] * msd_143[k];

        t_240[k] = f_12 * lsd_144[k]
                   + f_1 * msp0_72[k]
                   - f_2 * msp1_72[k]
                   + f_3 * pc_x[k] * msd_144[k];
    }

#pragma omp simd aligned(t_241, t_242, t_243, t_244, pc_x, pc_y, pc_z, lsd_102, lsd_108, \
                         lsd_147, lsd_148, msd_144, msd_147, msd_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_241[k] = f_12 * lsd_108[k]
                   + f_3 * pc_y[k] * msd_144[k];

        t_242[k] = f_12 * lsd_102[k]
                   + f_3 * pc_z[k] * msd_144[k];

        t_243[k] = f_12 * lsd_147[k]
                   + f_3 * pc_x[k] * msd_147[k];

        t_244[k] = f_12 * lsd_148[k]
                   + f_3 * pc_x[k] * msd_148[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, pc_x, pc_y, pc_z, lsd_105, lsd_111, \
                         lsd_113, lsd_149, msp0_73, msp1_73, msd_147, \
                         msd_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_12 * lsd_149[k]
                   + f_3 * pc_x[k] * msd_149[k];

        t_246[k] = f_12 * lsd_111[k]
                   + f_1 * msp0_73[k]
                   - f_2 * msp1_73[k]
                   + f_3 * pc_y[k] * msd_147[k];

        t_247[k] = f_12 * lsd_105[k]
                   + f_3 * pc_z[k] * msd_147[k];

        t_248[k] = f_12 * lsd_113[k]
                   + f_3 * pc_y[k] * msd_149[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, pc_x, pc_y, pc_z, lsd_107, lsd_114, lsd_150, \
                         msp0_74, msp0_75, msp1_74, msp1_75, msd_149, \
                         msd_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = f_12 * lsd_107[k]
                   + f_1 * msp0_74[k]
                   - f_2 * msp1_74[k]
                   + f_3 * pc_z[k] * msd_149[k];

        t_250[k] = f_12 * lsd_150[k]
                   + f_1 * msp0_75[k]
                   - f_2 * msp1_75[k]
                   + f_3 * pc_x[k] * msd_150[k];

        t_251[k] = f_10 * lsd_114[k]
                   + f_3 * pc_y[k] * msd_150[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, pc_x, pc_z, lsd_108, lsd_153, lsd_154, \
                         lsd_155, msd_150, msd_153, msd_154, msd_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = f_14 * lsd_108[k]
                   + f_3 * pc_z[k] * msd_150[k];

        t_253[k] = f_12 * lsd_153[k]
                   + f_3 * pc_x[k] * msd_153[k];

        t_254[k] = f_12 * lsd_154[k]
                   + f_3 * pc_x[k] * msd_154[k];

        t_255[k] = f_12 * lsd_155[k]
                   + f_3 * pc_x[k] * msd_155[k];
    }
}

static auto
compute_prim_msf_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsf0,
                                                          const size_t lsd, const size_t lsf1,
                                                          const size_t msp0, const size_t msp1,
                                                          const size_t msd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 4.0 / q;
    const auto f_7 = 0.5 / gamma;
    const auto f_8 = 0.5 * p / (gamma * q);
    const auto f_9 = 3.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 3.0 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 2.5 / q;
    const auto f_14 = 2.0 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *lsf0_200 = buffer.data(lsf0 + 200);
    const auto *lsf0_209 = buffer.data(lsf0 + 209);
    const auto *lsf0_210 = buffer.data(lsf0 + 210);
    const auto *lsf0_216 = buffer.data(lsf0 + 216);
    const auto *lsf0_270 = buffer.data(lsf0 + 270);
    const auto *lsf0_279 = buffer.data(lsf0 + 279);
    const auto *lsf0_280 = buffer.data(lsf0 + 280);
    const auto *lsf0_360 = buffer.data(lsf0 + 360);
    const auto *lsf0_366 = buffer.data(lsf0 + 366);
    const auto *lsf0_369 = buffer.data(lsf0 + 369);
    const auto *lsf0_376 = buffer.data(lsf0 + 376);
    const auto *lsf0_379 = buffer.data(lsf0 + 379);
    const auto *lsf0_380 = buffer.data(lsf0 + 380);

    const auto *lsd_111 = buffer.data(lsd + 111);
    const auto *lsd_113 = buffer.data(lsd + 113);
    const auto *lsd_114 = buffer.data(lsd + 114);
    const auto *lsd_117 = buffer.data(lsd + 117);
    const auto *lsd_119 = buffer.data(lsd + 119);
    const auto *lsd_120 = buffer.data(lsd + 120);
    const auto *lsd_123 = buffer.data(lsd + 123);
    const auto *lsd_125 = buffer.data(lsd + 125);
    const auto *lsd_126 = buffer.data(lsd + 126);
    const auto *lsd_129 = buffer.data(lsd + 129);
    const auto *lsd_131 = buffer.data(lsd + 131);
    const auto *lsd_132 = buffer.data(lsd + 132);
    const auto *lsd_135 = buffer.data(lsd + 135);
    const auto *lsd_137 = buffer.data(lsd + 137);
    const auto *lsd_138 = buffer.data(lsd + 138);
    const auto *lsd_141 = buffer.data(lsd + 141);
    const auto *lsd_143 = buffer.data(lsd + 143);
    const auto *lsd_144 = buffer.data(lsd + 144);
    const auto *lsd_147 = buffer.data(lsd + 147);
    const auto *lsd_149 = buffer.data(lsd + 149);
    const auto *lsd_150 = buffer.data(lsd + 150);
    const auto *lsd_153 = buffer.data(lsd + 153);
    const auto *lsd_155 = buffer.data(lsd + 155);
    const auto *lsd_156 = buffer.data(lsd + 156);
    const auto *lsd_159 = buffer.data(lsd + 159);
    const auto *lsd_160 = buffer.data(lsd + 160);
    const auto *lsd_161 = buffer.data(lsd + 161);
    const auto *lsd_162 = buffer.data(lsd + 162);
    const auto *lsd_165 = buffer.data(lsd + 165);
    const auto *lsd_167 = buffer.data(lsd + 167);
    const auto *lsd_168 = buffer.data(lsd + 168);
    const auto *lsd_171 = buffer.data(lsd + 171);
    const auto *lsd_173 = buffer.data(lsd + 173);
    const auto *lsd_174 = buffer.data(lsd + 174);
    const auto *lsd_177 = buffer.data(lsd + 177);
    const auto *lsd_178 = buffer.data(lsd + 178);
    const auto *lsd_179 = buffer.data(lsd + 179);
    const auto *lsd_180 = buffer.data(lsd + 180);
    const auto *lsd_183 = buffer.data(lsd + 183);
    const auto *lsd_184 = buffer.data(lsd + 184);
    const auto *lsd_185 = buffer.data(lsd + 185);
    const auto *lsd_186 = buffer.data(lsd + 186);
    const auto *lsd_189 = buffer.data(lsd + 189);
    const auto *lsd_190 = buffer.data(lsd + 190);
    const auto *lsd_191 = buffer.data(lsd + 191);
    const auto *lsd_192 = buffer.data(lsd + 192);
    const auto *lsd_195 = buffer.data(lsd + 195);
    const auto *lsd_196 = buffer.data(lsd + 196);
    const auto *lsd_197 = buffer.data(lsd + 197);
    const auto *lsd_198 = buffer.data(lsd + 198);
    const auto *lsd_201 = buffer.data(lsd + 201);
    const auto *lsd_202 = buffer.data(lsd + 202);
    const auto *lsd_203 = buffer.data(lsd + 203);
    const auto *lsd_207 = buffer.data(lsd + 207);
    const auto *lsd_208 = buffer.data(lsd + 208);
    const auto *lsd_209 = buffer.data(lsd + 209);
    const auto *lsd_210 = buffer.data(lsd + 210);
    const auto *lsd_213 = buffer.data(lsd + 213);
    const auto *lsd_215 = buffer.data(lsd + 215);
    const auto *lsd_216 = buffer.data(lsd + 216);
    const auto *lsd_219 = buffer.data(lsd + 219);
    const auto *lsd_221 = buffer.data(lsd + 221);
    const auto *lsd_225 = buffer.data(lsd + 225);
    const auto *lsd_226 = buffer.data(lsd + 226);
    const auto *lsd_227 = buffer.data(lsd + 227);
    const auto *lsd_228 = buffer.data(lsd + 228);

    const auto *lsf1_200 = buffer.data(lsf1 + 200);
    const auto *lsf1_209 = buffer.data(lsf1 + 209);
    const auto *lsf1_210 = buffer.data(lsf1 + 210);
    const auto *lsf1_216 = buffer.data(lsf1 + 216);
    const auto *lsf1_270 = buffer.data(lsf1 + 270);
    const auto *lsf1_279 = buffer.data(lsf1 + 279);
    const auto *lsf1_280 = buffer.data(lsf1 + 280);
    const auto *lsf1_360 = buffer.data(lsf1 + 360);
    const auto *lsf1_366 = buffer.data(lsf1 + 366);
    const auto *lsf1_369 = buffer.data(lsf1 + 369);
    const auto *lsf1_376 = buffer.data(lsf1 + 376);
    const auto *lsf1_379 = buffer.data(lsf1 + 379);
    const auto *lsf1_380 = buffer.data(lsf1 + 380);

    const auto *msp0_76 = buffer.data(msp0 + 76);
    const auto *msp0_77 = buffer.data(msp0 + 77);
    const auto *msp0_79 = buffer.data(msp0 + 79);
    const auto *msp0_81 = buffer.data(msp0 + 81);
    const auto *msp0_82 = buffer.data(msp0 + 82);
    const auto *msp0_83 = buffer.data(msp0 + 83);
    const auto *msp0_84 = buffer.data(msp0 + 84);
    const auto *msp0_85 = buffer.data(msp0 + 85);
    const auto *msp0_86 = buffer.data(msp0 + 86);
    const auto *msp0_89 = buffer.data(msp0 + 89);
    const auto *msp0_90 = buffer.data(msp0 + 90);
    const auto *msp0_91 = buffer.data(msp0 + 91);
    const auto *msp0_92 = buffer.data(msp0 + 92);
    const auto *msp0_93 = buffer.data(msp0 + 93);
    const auto *msp0_94 = buffer.data(msp0 + 94);
    const auto *msp0_95 = buffer.data(msp0 + 95);
    const auto *msp0_96 = buffer.data(msp0 + 96);
    const auto *msp0_97 = buffer.data(msp0 + 97);
    const auto *msp0_98 = buffer.data(msp0 + 98);
    const auto *msp0_99 = buffer.data(msp0 + 99);
    const auto *msp0_100 = buffer.data(msp0 + 100);
    const auto *msp0_101 = buffer.data(msp0 + 101);
    const auto *msp0_103 = buffer.data(msp0 + 103);
    const auto *msp0_105 = buffer.data(msp0 + 105);
    const auto *msp0_106 = buffer.data(msp0 + 106);
    const auto *msp0_107 = buffer.data(msp0 + 107);

    const auto *msp1_76 = buffer.data(msp1 + 76);
    const auto *msp1_77 = buffer.data(msp1 + 77);
    const auto *msp1_79 = buffer.data(msp1 + 79);
    const auto *msp1_81 = buffer.data(msp1 + 81);
    const auto *msp1_82 = buffer.data(msp1 + 82);
    const auto *msp1_83 = buffer.data(msp1 + 83);
    const auto *msp1_84 = buffer.data(msp1 + 84);
    const auto *msp1_85 = buffer.data(msp1 + 85);
    const auto *msp1_86 = buffer.data(msp1 + 86);
    const auto *msp1_89 = buffer.data(msp1 + 89);
    const auto *msp1_90 = buffer.data(msp1 + 90);
    const auto *msp1_91 = buffer.data(msp1 + 91);
    const auto *msp1_92 = buffer.data(msp1 + 92);
    const auto *msp1_93 = buffer.data(msp1 + 93);
    const auto *msp1_94 = buffer.data(msp1 + 94);
    const auto *msp1_95 = buffer.data(msp1 + 95);
    const auto *msp1_96 = buffer.data(msp1 + 96);
    const auto *msp1_97 = buffer.data(msp1 + 97);
    const auto *msp1_98 = buffer.data(msp1 + 98);
    const auto *msp1_99 = buffer.data(msp1 + 99);
    const auto *msp1_100 = buffer.data(msp1 + 100);
    const auto *msp1_101 = buffer.data(msp1 + 101);
    const auto *msp1_103 = buffer.data(msp1 + 103);
    const auto *msp1_105 = buffer.data(msp1 + 105);
    const auto *msp1_106 = buffer.data(msp1 + 106);
    const auto *msp1_107 = buffer.data(msp1 + 107);

    const auto *msd_153 = buffer.data(msd + 153);
    const auto *msd_155 = buffer.data(msd + 155);
    const auto *msd_156 = buffer.data(msd + 156);
    const auto *msd_159 = buffer.data(msd + 159);
    const auto *msd_160 = buffer.data(msd + 160);
    const auto *msd_161 = buffer.data(msd + 161);
    const auto *msd_162 = buffer.data(msd + 162);
    const auto *msd_164 = buffer.data(msd + 164);
    const auto *msd_165 = buffer.data(msd + 165);
    const auto *msd_166 = buffer.data(msd + 166);
    const auto *msd_167 = buffer.data(msd + 167);
    const auto *msd_168 = buffer.data(msd + 168);
    const auto *msd_169 = buffer.data(msd + 169);
    const auto *msd_171 = buffer.data(msd + 171);
    const auto *msd_173 = buffer.data(msd + 173);
    const auto *msd_174 = buffer.data(msd + 174);
    const auto *msd_177 = buffer.data(msd + 177);
    const auto *msd_178 = buffer.data(msd + 178);
    const auto *msd_179 = buffer.data(msd + 179);
    const auto *msd_180 = buffer.data(msd + 180);
    const auto *msd_183 = buffer.data(msd + 183);
    const auto *msd_184 = buffer.data(msd + 184);
    const auto *msd_185 = buffer.data(msd + 185);
    const auto *msd_186 = buffer.data(msd + 186);
    const auto *msd_189 = buffer.data(msd + 189);
    const auto *msd_190 = buffer.data(msd + 190);
    const auto *msd_191 = buffer.data(msd + 191);
    const auto *msd_192 = buffer.data(msd + 192);
    const auto *msd_195 = buffer.data(msd + 195);
    const auto *msd_196 = buffer.data(msd + 196);
    const auto *msd_197 = buffer.data(msd + 197);
    const auto *msd_198 = buffer.data(msd + 198);
    const auto *msd_201 = buffer.data(msd + 201);
    const auto *msd_202 = buffer.data(msd + 202);
    const auto *msd_203 = buffer.data(msd + 203);
    const auto *msd_204 = buffer.data(msd + 204);
    const auto *msd_207 = buffer.data(msd + 207);
    const auto *msd_208 = buffer.data(msd + 208);
    const auto *msd_209 = buffer.data(msd + 209);
    const auto *msd_210 = buffer.data(msd + 210);
    const auto *msd_212 = buffer.data(msd + 212);
    const auto *msd_213 = buffer.data(msd + 213);
    const auto *msd_214 = buffer.data(msd + 214);
    const auto *msd_215 = buffer.data(msd + 215);
    const auto *msd_216 = buffer.data(msd + 216);
    const auto *msd_217 = buffer.data(msd + 217);
    const auto *msd_219 = buffer.data(msd + 219);
    const auto *msd_221 = buffer.data(msd + 221);
    const auto *msd_222 = buffer.data(msd + 222);
    const auto *msd_225 = buffer.data(msd + 225);
    const auto *msd_226 = buffer.data(msd + 226);
    const auto *msd_227 = buffer.data(msd + 227);
    const auto *msd_228 = buffer.data(msd + 228);

#pragma omp simd aligned(t_256, t_257, t_258, t_259, pc_y, pc_z, lsd_111, lsd_113, lsd_117, \
                         lsd_119, msp0_76, msp0_77, msp1_76, msp1_77, msd_153, \
                         msd_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = f_10 * lsd_117[k]
                   + f_1 * msp0_76[k]
                   - f_2 * msp1_76[k]
                   + f_3 * pc_y[k] * msd_153[k];

        t_257[k] = f_14 * lsd_111[k]
                   + f_3 * pc_z[k] * msd_153[k];

        t_258[k] = f_10 * lsd_119[k]
                   + f_3 * pc_y[k] * msd_155[k];

        t_259[k] = f_14 * lsd_113[k]
                   + f_1 * msp0_77[k]
                   - f_2 * msp1_77[k]
                   + f_3 * pc_z[k] * msd_155[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pa_y, pc_x, pc_y, pc_z, lsf0_200, \
                         lsd_114, lsd_120, lsd_159, lsf1_200, msd_156, \
                         msd_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = pa_y[k] * lsf0_200[k]
                   - f_4 * pc_y[k] * lsf1_200[k];

        t_261[k] = f_5 * lsd_120[k]
                   + f_3 * pc_y[k] * msd_156[k];

        t_262[k] = f_13 * lsd_114[k]
                   + f_3 * pc_z[k] * msd_156[k];

        t_263[k] = f_12 * lsd_159[k]
                   + f_3 * pc_x[k] * msd_159[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, pc_x, pc_y, pc_z, lsd_117, lsd_123, \
                         lsd_160, lsd_161, msp0_79, msp1_79, msd_159, msd_160, \
                         msd_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_12 * lsd_160[k]
                   + f_3 * pc_x[k] * msd_160[k];

        t_265[k] = f_12 * lsd_161[k]
                   + f_3 * pc_x[k] * msd_161[k];

        t_266[k] = f_5 * lsd_123[k]
                   + f_1 * msp0_79[k]
                   - f_2 * msp1_79[k]
                   + f_3 * pc_y[k] * msd_159[k];

        t_267[k] = f_13 * lsd_117[k]
                   + f_3 * pc_z[k] * msd_159[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, t_271, pa_y, pc_x, pc_y, lsf0_209, lsd_125, \
                         lsd_162, lsf1_209, msp0_81, msp1_81, msd_161, \
                         msd_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_5 * lsd_125[k]
                   + f_3 * pc_y[k] * msd_161[k];

        t_269[k] = pa_y[k] * lsf0_209[k]
                   - f_4 * pc_y[k] * lsf1_209[k];

        t_270[k] = f_12 * lsd_162[k]
                   + f_1 * msp0_81[k]
                   - f_2 * msp1_81[k]
                   + f_3 * pc_x[k] * msd_162[k];

        t_271[k] = f_3 * pc_y[k] * msd_162[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, pc_x, pc_y, pc_z, lsd_120, lsd_165, \
                         lsd_167, msd_162, msd_164, msd_165, msd_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_11 * lsd_120[k]
                   + f_3 * pc_z[k] * msd_162[k];

        t_273[k] = f_12 * lsd_165[k]
                   + f_3 * pc_x[k] * msd_165[k];

        t_274[k] = f_3 * pc_y[k] * msd_164[k];

        t_275[k] = f_12 * lsd_167[k]
                   + f_3 * pc_x[k] * msd_167[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pc_y, pc_z, lsd_125, msp0_82, msp0_83, \
                         msp1_82, msp1_83, msd_165, msd_166, msd_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_1 * msp0_82[k]
                   - f_2 * msp1_82[k]
                   + f_3 * pc_y[k] * msd_165[k];

        t_277[k] = f_7 * msp0_83[k]
                   - f_8 * msp1_83[k]
                   + f_3 * pc_y[k] * msd_166[k];

        t_278[k] = f_3 * pc_y[k] * msd_167[k];

        t_279[k] = f_11 * lsd_125[k]
                   + f_1 * msp0_83[k]
                   - f_2 * msp1_83[k]
                   + f_3 * pc_z[k] * msd_167[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, pc_x, pc_y, pc_z, lsd_126, \
                         lsd_168, lsd_171, msp0_84, msp1_84, msd_168, msd_169, \
                         msd_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_10 * lsd_168[k]
                   + f_1 * msp0_84[k]
                   - f_2 * msp1_84[k]
                   + f_3 * pc_x[k] * msd_168[k];

        t_281[k] = f_9 * lsd_126[k]
                   + f_3 * pc_y[k] * msd_168[k];

        t_282[k] = f_3 * pc_z[k] * msd_168[k];

        t_283[k] = f_10 * lsd_171[k]
                   + f_3 * pc_x[k] * msd_171[k];

        t_284[k] = f_3 * pc_z[k] * msd_169[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, pc_x, pc_y, pc_z, lsd_129, lsd_131, \
                         lsd_173, msp0_85, msp1_85, msd_171, msd_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = f_10 * lsd_173[k]
                   + f_3 * pc_x[k] * msd_173[k];

        t_286[k] = f_9 * lsd_129[k]
                   + f_1 * msp0_85[k]
                   - f_2 * msp1_85[k]
                   + f_3 * pc_y[k] * msd_171[k];

        t_287[k] = f_3 * pc_z[k] * msd_171[k];

        t_288[k] = f_9 * lsd_131[k]
                   + f_3 * pc_y[k] * msd_173[k];
    }

#pragma omp simd aligned(t_289, t_290, t_291, t_292, pa_z, pc_y, pc_z, lsf0_210, lsd_126, \
                         lsd_132, lsf1_210, msp0_86, msp1_86, msd_173, \
                         msd_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_289[k] = f_1 * msp0_86[k]
                   - f_2 * msp1_86[k]
                   + f_3 * pc_z[k] * msd_173[k];

        t_290[k] = pa_z[k] * lsf0_210[k]
                   - f_4 * pc_z[k] * lsf1_210[k];

        t_291[k] = f_11 * lsd_132[k]
                   + f_3 * pc_y[k] * msd_174[k];

        t_292[k] = f_5 * lsd_126[k]
                   + f_3 * pc_z[k] * msd_174[k];
    }

#pragma omp simd aligned(t_293, t_294, t_295, t_296, pa_z, pc_x, pc_z, lsf0_216, lsd_177, \
                         lsd_178, lsd_179, lsf1_216, msd_177, msd_178, \
                         msd_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_293[k] = f_10 * lsd_177[k]
                   + f_3 * pc_x[k] * msd_177[k];

        t_294[k] = f_10 * lsd_178[k]
                   + f_3 * pc_x[k] * msd_178[k];

        t_295[k] = f_10 * lsd_179[k]
                   + f_3 * pc_x[k] * msd_179[k];

        t_296[k] = pa_z[k] * lsf0_216[k]
                   - f_4 * pc_z[k] * lsf1_216[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, pc_y, pc_z, lsd_129, lsd_131, lsd_137, msp0_89, \
                         msp1_89, msd_177, msd_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = f_5 * lsd_129[k]
                   + f_3 * pc_z[k] * msd_177[k];

        t_298[k] = f_11 * lsd_137[k]
                   + f_3 * pc_y[k] * msd_179[k];

        t_299[k] = f_5 * lsd_131[k]
                   + f_1 * msp0_89[k]
                   - f_2 * msp1_89[k]
                   + f_3 * pc_z[k] * msd_179[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, pc_x, pc_y, pc_z, lsd_132, lsd_138, \
                         lsd_180, lsd_183, msp0_90, msp1_90, msd_180, \
                         msd_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = f_10 * lsd_180[k]
                   + f_1 * msp0_90[k]
                   - f_2 * msp1_90[k]
                   + f_3 * pc_x[k] * msd_180[k];

        t_301[k] = f_13 * lsd_138[k]
                   + f_3 * pc_y[k] * msd_180[k];

        t_302[k] = f_10 * lsd_132[k]
                   + f_3 * pc_z[k] * msd_180[k];

        t_303[k] = f_10 * lsd_183[k]
                   + f_3 * pc_x[k] * msd_183[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, pc_x, pc_y, pc_z, lsd_135, lsd_141, \
                         lsd_184, lsd_185, msp0_91, msp1_91, msd_183, msd_184, \
                         msd_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_304[k] = f_10 * lsd_184[k]
                   + f_3 * pc_x[k] * msd_184[k];

        t_305[k] = f_10 * lsd_185[k]
                   + f_3 * pc_x[k] * msd_185[k];

        t_306[k] = f_13 * lsd_141[k]
                   + f_1 * msp0_91[k]
                   - f_2 * msp1_91[k]
                   + f_3 * pc_y[k] * msd_183[k];

        t_307[k] = f_10 * lsd_135[k]
                   + f_3 * pc_z[k] * msd_183[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, pc_x, pc_y, pc_z, lsd_137, lsd_143, lsd_186, \
                         msp0_92, msp0_93, msp1_92, msp1_93, msd_185, \
                         msd_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = f_13 * lsd_143[k]
                   + f_3 * pc_y[k] * msd_185[k];

        t_309[k] = f_10 * lsd_137[k]
                   + f_1 * msp0_92[k]
                   - f_2 * msp1_92[k]
                   + f_3 * pc_z[k] * msd_185[k];

        t_310[k] = f_10 * lsd_186[k]
                   + f_1 * msp0_93[k]
                   - f_2 * msp1_93[k]
                   + f_3 * pc_x[k] * msd_186[k];
    }

#pragma omp simd aligned(t_311, t_312, t_313, t_314, pc_x, pc_y, pc_z, lsd_138, lsd_144, \
                         lsd_189, lsd_190, msd_186, msd_189, msd_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_311[k] = f_14 * lsd_144[k]
                   + f_3 * pc_y[k] * msd_186[k];

        t_312[k] = f_12 * lsd_138[k]
                   + f_3 * pc_z[k] * msd_186[k];

        t_313[k] = f_10 * lsd_189[k]
                   + f_3 * pc_x[k] * msd_189[k];

        t_314[k] = f_10 * lsd_190[k]
                   + f_3 * pc_x[k] * msd_190[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, pc_x, pc_y, pc_z, lsd_141, lsd_147, \
                         lsd_149, lsd_191, msp0_94, msp1_94, msd_189, \
                         msd_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = f_10 * lsd_191[k]
                   + f_3 * pc_x[k] * msd_191[k];

        t_316[k] = f_14 * lsd_147[k]
                   + f_1 * msp0_94[k]
                   - f_2 * msp1_94[k]
                   + f_3 * pc_y[k] * msd_189[k];

        t_317[k] = f_12 * lsd_141[k]
                   + f_3 * pc_z[k] * msd_189[k];

        t_318[k] = f_14 * lsd_149[k]
                   + f_3 * pc_y[k] * msd_191[k];
    }

#pragma omp simd aligned(t_319, t_320, t_321, pc_x, pc_y, pc_z, lsd_143, lsd_150, lsd_192, \
                         msp0_95, msp0_96, msp1_95, msp1_96, msd_191, \
                         msd_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_319[k] = f_12 * lsd_143[k]
                   + f_1 * msp0_95[k]
                   - f_2 * msp1_95[k]
                   + f_3 * pc_z[k] * msd_191[k];

        t_320[k] = f_10 * lsd_192[k]
                   + f_1 * msp0_96[k]
                   - f_2 * msp1_96[k]
                   + f_3 * pc_x[k] * msd_192[k];

        t_321[k] = f_12 * lsd_150[k]
                   + f_3 * pc_y[k] * msd_192[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, t_325, pc_x, pc_z, lsd_144, lsd_195, lsd_196, \
                         lsd_197, msd_192, msd_195, msd_196, msd_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = f_14 * lsd_144[k]
                   + f_3 * pc_z[k] * msd_192[k];

        t_323[k] = f_10 * lsd_195[k]
                   + f_3 * pc_x[k] * msd_195[k];

        t_324[k] = f_10 * lsd_196[k]
                   + f_3 * pc_x[k] * msd_196[k];

        t_325[k] = f_10 * lsd_197[k]
                   + f_3 * pc_x[k] * msd_197[k];
    }

#pragma omp simd aligned(t_326, t_327, t_328, t_329, pc_y, pc_z, lsd_147, lsd_149, lsd_153, \
                         lsd_155, msp0_97, msp0_98, msp1_97, msp1_98, msd_195, \
                         msd_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_326[k] = f_12 * lsd_153[k]
                   + f_1 * msp0_97[k]
                   - f_2 * msp1_97[k]
                   + f_3 * pc_y[k] * msd_195[k];

        t_327[k] = f_14 * lsd_147[k]
                   + f_3 * pc_z[k] * msd_195[k];

        t_328[k] = f_12 * lsd_155[k]
                   + f_3 * pc_y[k] * msd_197[k];

        t_329[k] = f_14 * lsd_149[k]
                   + f_1 * msp0_98[k]
                   - f_2 * msp1_98[k]
                   + f_3 * pc_z[k] * msd_197[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, pc_x, pc_y, pc_z, lsd_150, lsd_156, \
                         lsd_198, lsd_201, msp0_99, msp1_99, msd_198, \
                         msd_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = f_10 * lsd_198[k]
                   + f_1 * msp0_99[k]
                   - f_2 * msp1_99[k]
                   + f_3 * pc_x[k] * msd_198[k];

        t_331[k] = f_10 * lsd_156[k]
                   + f_3 * pc_y[k] * msd_198[k];

        t_332[k] = f_13 * lsd_150[k]
                   + f_3 * pc_z[k] * msd_198[k];

        t_333[k] = f_10 * lsd_201[k]
                   + f_3 * pc_x[k] * msd_201[k];
    }

#pragma omp simd aligned(t_334, t_335, t_336, t_337, pc_x, pc_y, pc_z, lsd_153, lsd_159, \
                         lsd_202, lsd_203, msp0_100, msp1_100, msd_201, msd_202, \
                         msd_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_334[k] = f_10 * lsd_202[k]
                   + f_3 * pc_x[k] * msd_202[k];

        t_335[k] = f_10 * lsd_203[k]
                   + f_3 * pc_x[k] * msd_203[k];

        t_336[k] = f_10 * lsd_159[k]
                   + f_1 * msp0_100[k]
                   - f_2 * msp1_100[k]
                   + f_3 * pc_y[k] * msd_201[k];

        t_337[k] = f_13 * lsd_153[k]
                   + f_3 * pc_z[k] * msd_201[k];
    }

#pragma omp simd aligned(t_338, t_339, t_340, t_341, pa_y, pc_y, pc_z, lsf0_270, lsd_155, \
                         lsd_161, lsd_162, lsf1_270, msp0_101, msp1_101, msd_203, \
                         msd_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_338[k] = f_10 * lsd_161[k]
                   + f_3 * pc_y[k] * msd_203[k];

        t_339[k] = f_13 * lsd_155[k]
                   + f_1 * msp0_101[k]
                   - f_2 * msp1_101[k]
                   + f_3 * pc_z[k] * msd_203[k];

        t_340[k] = pa_y[k] * lsf0_270[k]
                   - f_4 * pc_y[k] * lsf1_270[k];

        t_341[k] = f_5 * lsd_162[k]
                   + f_3 * pc_y[k] * msd_204[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, t_345, pc_x, pc_z, lsd_156, lsd_207, lsd_208, \
                         lsd_209, msd_204, msd_207, msd_208, msd_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = f_11 * lsd_156[k]
                   + f_3 * pc_z[k] * msd_204[k];

        t_343[k] = f_10 * lsd_207[k]
                   + f_3 * pc_x[k] * msd_207[k];

        t_344[k] = f_10 * lsd_208[k]
                   + f_3 * pc_x[k] * msd_208[k];

        t_345[k] = f_10 * lsd_209[k]
                   + f_3 * pc_x[k] * msd_209[k];
    }

#pragma omp simd aligned(t_346, t_347, t_348, t_349, pa_y, pc_y, pc_z, lsf0_279, lsd_159, \
                         lsd_165, lsd_167, lsf1_279, msp0_103, msp1_103, msd_207, \
                         msd_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_346[k] = f_5 * lsd_165[k]
                   + f_1 * msp0_103[k]
                   - f_2 * msp1_103[k]
                   + f_3 * pc_y[k] * msd_207[k];

        t_347[k] = f_11 * lsd_159[k]
                   + f_3 * pc_z[k] * msd_207[k];

        t_348[k] = f_5 * lsd_167[k]
                   + f_3 * pc_y[k] * msd_209[k];

        t_349[k] = pa_y[k] * lsf0_279[k]
                   - f_4 * pc_y[k] * lsf1_279[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, pc_x, pc_y, pc_z, lsd_162, \
                         lsd_210, lsd_213, msp0_105, msp1_105, msd_210, msd_212, \
                         msd_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = f_10 * lsd_210[k]
                   + f_1 * msp0_105[k]
                   - f_2 * msp1_105[k]
                   + f_3 * pc_x[k] * msd_210[k];

        t_351[k] = f_3 * pc_y[k] * msd_210[k];

        t_352[k] = f_9 * lsd_162[k]
                   + f_3 * pc_z[k] * msd_210[k];

        t_353[k] = f_10 * lsd_213[k]
                   + f_3 * pc_x[k] * msd_213[k];

        t_354[k] = f_3 * pc_y[k] * msd_212[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, pc_x, pc_y, lsd_215, msp0_106, msp0_107, \
                         msp1_106, msp1_107, msd_213, msd_214, \
                         msd_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = f_10 * lsd_215[k]
                   + f_3 * pc_x[k] * msd_215[k];

        t_356[k] = f_1 * msp0_106[k]
                   - f_2 * msp1_106[k]
                   + f_3 * pc_y[k] * msd_213[k];

        t_357[k] = f_7 * msp0_107[k]
                   - f_8 * msp1_107[k]
                   + f_3 * pc_y[k] * msd_214[k];

        t_358[k] = f_3 * pc_y[k] * msd_215[k];
    }

#pragma omp simd aligned(t_359, t_360, t_361, pa_x, pc_x, pc_y, pc_z, lsf0_360, lsd_167, \
                         lsd_168, lsd_216, lsf1_360, msp0_107, msp1_107, msd_215, \
                         msd_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = f_9 * lsd_167[k]
                   + f_1 * msp0_107[k]
                   - f_2 * msp1_107[k]
                   + f_3 * pc_z[k] * msd_215[k];

        t_360[k] = pa_x[k] * lsf0_360[k]
                   + f_12 * lsd_216[k]
                   - f_4 * pc_x[k] * lsf1_360[k];

        t_361[k] = f_6 * lsd_168[k]
                   + f_3 * pc_y[k] * msd_216[k];
    }

#pragma omp simd aligned(t_362, t_363, t_364, t_365, t_366, pa_x, pc_x, pc_z, lsf0_366, \
                         lsd_219, lsd_221, lsf1_366, msd_216, msd_217, msd_219, \
                         msd_221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = f_3 * pc_z[k] * msd_216[k];

        t_363[k] = f_5 * lsd_219[k]
                   + f_3 * pc_x[k] * msd_219[k];

        t_364[k] = f_3 * pc_z[k] * msd_217[k];

        t_365[k] = f_5 * lsd_221[k]
                   + f_3 * pc_x[k] * msd_221[k];

        t_366[k] = pa_x[k] * lsf0_366[k]
                   - f_4 * pc_x[k] * lsf1_366[k];
    }

#pragma omp simd aligned(t_367, t_368, t_369, t_370, pa_x, pa_z, pc_x, pc_y, pc_z, lsf0_280, \
                         lsf0_369, lsd_173, lsf1_280, lsf1_369, msd_219, \
                         msd_221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_367[k] = f_3 * pc_z[k] * msd_219[k];

        t_368[k] = f_6 * lsd_173[k]
                   + f_3 * pc_y[k] * msd_221[k];

        t_369[k] = pa_x[k] * lsf0_369[k]
                   - f_4 * pc_x[k] * lsf1_369[k];

        t_370[k] = pa_z[k] * lsf0_280[k]
                   - f_4 * pc_z[k] * lsf1_280[k];
    }

#pragma omp simd aligned(t_371, t_372, t_373, t_374, pc_x, pc_y, pc_z, lsd_168, lsd_174, \
                         lsd_225, lsd_226, msd_222, msd_225, msd_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_371[k] = f_9 * lsd_174[k]
                   + f_3 * pc_y[k] * msd_222[k];

        t_372[k] = f_5 * lsd_168[k]
                   + f_3 * pc_z[k] * msd_222[k];

        t_373[k] = f_5 * lsd_225[k]
                   + f_3 * pc_x[k] * msd_225[k];

        t_374[k] = f_5 * lsd_226[k]
                   + f_3 * pc_x[k] * msd_226[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, pa_x, pc_x, pc_y, pc_z, lsf0_376, \
                         lsd_171, lsd_179, lsd_227, lsf1_376, msd_225, \
                         msd_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = f_5 * lsd_227[k]
                   + f_3 * pc_x[k] * msd_227[k];

        t_376[k] = pa_x[k] * lsf0_376[k]
                   - f_4 * pc_x[k] * lsf1_376[k];

        t_377[k] = f_5 * lsd_171[k]
                   + f_3 * pc_z[k] * msd_225[k];

        t_378[k] = f_9 * lsd_179[k]
                   + f_3 * pc_y[k] * msd_227[k];
    }

#pragma omp simd aligned(t_379, t_380, t_381, t_382, pa_x, pc_x, pc_y, pc_z, lsf0_379, \
                         lsf0_380, lsd_174, lsd_180, lsd_228, lsf1_379, lsf1_380, \
                         msd_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_379[k] = pa_x[k] * lsf0_379[k]
                   - f_4 * pc_x[k] * lsf1_379[k];

        t_380[k] = pa_x[k] * lsf0_380[k]
                   + f_12 * lsd_228[k]
                   - f_4 * pc_x[k] * lsf1_380[k];

        t_381[k] = f_11 * lsd_180[k]
                   + f_3 * pc_y[k] * msd_228[k];

        t_382[k] = f_10 * lsd_174[k]
                   + f_3 * pc_z[k] * msd_228[k];
    }
}

static auto
compute_prim_msf_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsf0,
                                                          const size_t lsd, const size_t lsf1,
                                                          const size_t msp0, const size_t msp1,
                                                          const size_t msd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / q;
    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 4.0 / q;
    const auto f_7 = 0.5 / gamma;
    const auto f_8 = 0.5 * p / (gamma * q);
    const auto f_9 = 3.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 3.0 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 2.5 / q;
    const auto f_14 = 2.0 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *lsf0_350 = buffer.data(lsf0 + 350);
    const auto *lsf0_360 = buffer.data(lsf0 + 360);
    const auto *lsf0_361 = buffer.data(lsf0 + 361);
    const auto *lsf0_366 = buffer.data(lsf0 + 366);
    const auto *lsf0_386 = buffer.data(lsf0 + 386);
    const auto *lsf0_389 = buffer.data(lsf0 + 389);
    const auto *lsf0_390 = buffer.data(lsf0 + 390);
    const auto *lsf0_396 = buffer.data(lsf0 + 396);
    const auto *lsf0_399 = buffer.data(lsf0 + 399);
    const auto *lsf0_400 = buffer.data(lsf0 + 400);
    const auto *lsf0_406 = buffer.data(lsf0 + 406);
    const auto *lsf0_409 = buffer.data(lsf0 + 409);
    const auto *lsf0_410 = buffer.data(lsf0 + 410);
    const auto *lsf0_416 = buffer.data(lsf0 + 416);
    const auto *lsf0_419 = buffer.data(lsf0 + 419);
    const auto *lsf0_420 = buffer.data(lsf0 + 420);
    const auto *lsf0_426 = buffer.data(lsf0 + 426);
    const auto *lsf0_429 = buffer.data(lsf0 + 429);
    const auto *lsf0_436 = buffer.data(lsf0 + 436);
    const auto *lsf0_439 = buffer.data(lsf0 + 439);
    const auto *lsf0_440 = buffer.data(lsf0 + 440);
    const auto *lsf0_446 = buffer.data(lsf0 + 446);
    const auto *lsf0_447 = buffer.data(lsf0 + 447);
    const auto *lsf0_449 = buffer.data(lsf0 + 449);

    const auto *lsd_177 = buffer.data(lsd + 177);
    const auto *lsd_180 = buffer.data(lsd + 180);
    const auto *lsd_183 = buffer.data(lsd + 183);
    const auto *lsd_185 = buffer.data(lsd + 185);
    const auto *lsd_186 = buffer.data(lsd + 186);
    const auto *lsd_189 = buffer.data(lsd + 189);
    const auto *lsd_191 = buffer.data(lsd + 191);
    const auto *lsd_192 = buffer.data(lsd + 192);
    const auto *lsd_195 = buffer.data(lsd + 195);
    const auto *lsd_197 = buffer.data(lsd + 197);
    const auto *lsd_198 = buffer.data(lsd + 198);
    const auto *lsd_201 = buffer.data(lsd + 201);
    const auto *lsd_203 = buffer.data(lsd + 203);
    const auto *lsd_204 = buffer.data(lsd + 204);
    const auto *lsd_207 = buffer.data(lsd + 207);
    const auto *lsd_209 = buffer.data(lsd + 209);
    const auto *lsd_210 = buffer.data(lsd + 210);
    const auto *lsd_215 = buffer.data(lsd + 215);
    const auto *lsd_219 = buffer.data(lsd + 219);
    const auto *lsd_221 = buffer.data(lsd + 221);
    const auto *lsd_225 = buffer.data(lsd + 225);
    const auto *lsd_227 = buffer.data(lsd + 227);
    const auto *lsd_231 = buffer.data(lsd + 231);
    const auto *lsd_232 = buffer.data(lsd + 232);
    const auto *lsd_233 = buffer.data(lsd + 233);
    const auto *lsd_234 = buffer.data(lsd + 234);
    const auto *lsd_237 = buffer.data(lsd + 237);
    const auto *lsd_238 = buffer.data(lsd + 238);
    const auto *lsd_239 = buffer.data(lsd + 239);
    const auto *lsd_240 = buffer.data(lsd + 240);
    const auto *lsd_243 = buffer.data(lsd + 243);
    const auto *lsd_244 = buffer.data(lsd + 244);
    const auto *lsd_245 = buffer.data(lsd + 245);
    const auto *lsd_246 = buffer.data(lsd + 246);
    const auto *lsd_249 = buffer.data(lsd + 249);
    const auto *lsd_250 = buffer.data(lsd + 250);
    const auto *lsd_251 = buffer.data(lsd + 251);
    const auto *lsd_252 = buffer.data(lsd + 252);
    const auto *lsd_255 = buffer.data(lsd + 255);
    const auto *lsd_256 = buffer.data(lsd + 256);
    const auto *lsd_257 = buffer.data(lsd + 257);
    const auto *lsd_261 = buffer.data(lsd + 261);
    const auto *lsd_262 = buffer.data(lsd + 262);
    const auto *lsd_263 = buffer.data(lsd + 263);
    const auto *lsd_264 = buffer.data(lsd + 264);
    const auto *lsd_267 = buffer.data(lsd + 267);
    const auto *lsd_269 = buffer.data(lsd + 269);

    const auto *lsf1_350 = buffer.data(lsf1 + 350);
    const auto *lsf1_360 = buffer.data(lsf1 + 360);
    const auto *lsf1_361 = buffer.data(lsf1 + 361);
    const auto *lsf1_366 = buffer.data(lsf1 + 366);
    const auto *lsf1_386 = buffer.data(lsf1 + 386);
    const auto *lsf1_389 = buffer.data(lsf1 + 389);
    const auto *lsf1_390 = buffer.data(lsf1 + 390);
    const auto *lsf1_396 = buffer.data(lsf1 + 396);
    const auto *lsf1_399 = buffer.data(lsf1 + 399);
    const auto *lsf1_400 = buffer.data(lsf1 + 400);
    const auto *lsf1_406 = buffer.data(lsf1 + 406);
    const auto *lsf1_409 = buffer.data(lsf1 + 409);
    const auto *lsf1_410 = buffer.data(lsf1 + 410);
    const auto *lsf1_416 = buffer.data(lsf1 + 416);
    const auto *lsf1_419 = buffer.data(lsf1 + 419);
    const auto *lsf1_420 = buffer.data(lsf1 + 420);
    const auto *lsf1_426 = buffer.data(lsf1 + 426);
    const auto *lsf1_429 = buffer.data(lsf1 + 429);
    const auto *lsf1_436 = buffer.data(lsf1 + 436);
    const auto *lsf1_439 = buffer.data(lsf1 + 439);
    const auto *lsf1_440 = buffer.data(lsf1 + 440);
    const auto *lsf1_446 = buffer.data(lsf1 + 446);
    const auto *lsf1_447 = buffer.data(lsf1 + 447);
    const auto *lsf1_449 = buffer.data(lsf1 + 449);

    const auto *msp0_135 = buffer.data(msp0 + 135);
    const auto *msp0_136 = buffer.data(msp0 + 136);
    const auto *msp0_137 = buffer.data(msp0 + 137);
    const auto *msp0_140 = buffer.data(msp0 + 140);
    const auto *msp0_141 = buffer.data(msp0 + 141);
    const auto *msp0_142 = buffer.data(msp0 + 142);
    const auto *msp0_143 = buffer.data(msp0 + 143);
    const auto *msp0_144 = buffer.data(msp0 + 144);
    const auto *msp0_145 = buffer.data(msp0 + 145);
    const auto *msp0_146 = buffer.data(msp0 + 146);
    const auto *msp0_147 = buffer.data(msp0 + 147);
    const auto *msp0_148 = buffer.data(msp0 + 148);
    const auto *msp0_149 = buffer.data(msp0 + 149);
    const auto *msp0_150 = buffer.data(msp0 + 150);
    const auto *msp0_151 = buffer.data(msp0 + 151);
    const auto *msp0_152 = buffer.data(msp0 + 152);
    const auto *msp0_153 = buffer.data(msp0 + 153);
    const auto *msp0_154 = buffer.data(msp0 + 154);

    const auto *msp1_135 = buffer.data(msp1 + 135);
    const auto *msp1_136 = buffer.data(msp1 + 136);
    const auto *msp1_137 = buffer.data(msp1 + 137);
    const auto *msp1_140 = buffer.data(msp1 + 140);
    const auto *msp1_141 = buffer.data(msp1 + 141);
    const auto *msp1_142 = buffer.data(msp1 + 142);
    const auto *msp1_143 = buffer.data(msp1 + 143);
    const auto *msp1_144 = buffer.data(msp1 + 144);
    const auto *msp1_145 = buffer.data(msp1 + 145);
    const auto *msp1_146 = buffer.data(msp1 + 146);
    const auto *msp1_147 = buffer.data(msp1 + 147);
    const auto *msp1_148 = buffer.data(msp1 + 148);
    const auto *msp1_149 = buffer.data(msp1 + 149);
    const auto *msp1_150 = buffer.data(msp1 + 150);
    const auto *msp1_151 = buffer.data(msp1 + 151);
    const auto *msp1_152 = buffer.data(msp1 + 152);
    const auto *msp1_153 = buffer.data(msp1 + 153);
    const auto *msp1_154 = buffer.data(msp1 + 154);

    const auto *msd_231 = buffer.data(msd + 231);
    const auto *msd_232 = buffer.data(msd + 232);
    const auto *msd_233 = buffer.data(msd + 233);
    const auto *msd_234 = buffer.data(msd + 234);
    const auto *msd_237 = buffer.data(msd + 237);
    const auto *msd_238 = buffer.data(msd + 238);
    const auto *msd_239 = buffer.data(msd + 239);
    const auto *msd_240 = buffer.data(msd + 240);
    const auto *msd_243 = buffer.data(msd + 243);
    const auto *msd_244 = buffer.data(msd + 244);
    const auto *msd_245 = buffer.data(msd + 245);
    const auto *msd_246 = buffer.data(msd + 246);
    const auto *msd_249 = buffer.data(msd + 249);
    const auto *msd_250 = buffer.data(msd + 250);
    const auto *msd_251 = buffer.data(msd + 251);
    const auto *msd_252 = buffer.data(msd + 252);
    const auto *msd_255 = buffer.data(msd + 255);
    const auto *msd_256 = buffer.data(msd + 256);
    const auto *msd_257 = buffer.data(msd + 257);
    const auto *msd_258 = buffer.data(msd + 258);
    const auto *msd_261 = buffer.data(msd + 261);
    const auto *msd_262 = buffer.data(msd + 262);
    const auto *msd_263 = buffer.data(msd + 263);
    const auto *msd_264 = buffer.data(msd + 264);
    const auto *msd_266 = buffer.data(msd + 266);
    const auto *msd_267 = buffer.data(msd + 267);
    const auto *msd_269 = buffer.data(msd + 269);
    const auto *msd_270 = buffer.data(msd + 270);
    const auto *msd_271 = buffer.data(msd + 271);
    const auto *msd_273 = buffer.data(msd + 273);
    const auto *msd_274 = buffer.data(msd + 274);
    const auto *msd_275 = buffer.data(msd + 275);
    const auto *msd_278 = buffer.data(msd + 278);
    const auto *msd_279 = buffer.data(msd + 279);
    const auto *msd_280 = buffer.data(msd + 280);
    const auto *msd_281 = buffer.data(msd + 281);
    const auto *msd_282 = buffer.data(msd + 282);
    const auto *msd_283 = buffer.data(msd + 283);
    const auto *msd_284 = buffer.data(msd + 284);
    const auto *msd_285 = buffer.data(msd + 285);
    const auto *msd_286 = buffer.data(msd + 286);
    const auto *msd_287 = buffer.data(msd + 287);
    const auto *msd_288 = buffer.data(msd + 288);
    const auto *msd_289 = buffer.data(msd + 289);
    const auto *msd_290 = buffer.data(msd + 290);
    const auto *msd_291 = buffer.data(msd + 291);
    const auto *msd_292 = buffer.data(msd + 292);
    const auto *msd_293 = buffer.data(msd + 293);
    const auto *msd_294 = buffer.data(msd + 294);
    const auto *msd_295 = buffer.data(msd + 295);
    const auto *msd_296 = buffer.data(msd + 296);
    const auto *msd_297 = buffer.data(msd + 297);
    const auto *msd_298 = buffer.data(msd + 298);
    const auto *msd_299 = buffer.data(msd + 299);
    const auto *msd_300 = buffer.data(msd + 300);
    const auto *msd_301 = buffer.data(msd + 301);
    const auto *msd_302 = buffer.data(msd + 302);
    const auto *msd_303 = buffer.data(msd + 303);
    const auto *msd_304 = buffer.data(msd + 304);
    const auto *msd_305 = buffer.data(msd + 305);
    const auto *msd_306 = buffer.data(msd + 306);
    const auto *msd_307 = buffer.data(msd + 307);

#pragma omp simd aligned(t_383, t_384, t_385, t_386, pa_x, pc_x, lsf0_386, lsd_231, lsd_232, \
                         lsd_233, lsf1_386, msd_231, msd_232, msd_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_383[k] = f_5 * lsd_231[k]
                   + f_3 * pc_x[k] * msd_231[k];

        t_384[k] = f_5 * lsd_232[k]
                   + f_3 * pc_x[k] * msd_232[k];

        t_385[k] = f_5 * lsd_233[k]
                   + f_3 * pc_x[k] * msd_233[k];

        t_386[k] = pa_x[k] * lsf0_386[k]
                   - f_4 * pc_x[k] * lsf1_386[k];
    }

#pragma omp simd aligned(t_387, t_388, t_389, pa_x, pc_x, pc_y, pc_z, lsf0_389, lsd_177, \
                         lsd_185, lsf1_389, msd_231, msd_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_387[k] = f_10 * lsd_177[k]
                   + f_3 * pc_z[k] * msd_231[k];

        t_388[k] = f_11 * lsd_185[k]
                   + f_3 * pc_y[k] * msd_233[k];

        t_389[k] = pa_x[k] * lsf0_389[k]
                   - f_4 * pc_x[k] * lsf1_389[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, pa_x, pc_x, pc_y, pc_z, lsf0_390, \
                         lsd_180, lsd_186, lsd_234, lsd_237, lsf1_390, msd_234, \
                         msd_237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = pa_x[k] * lsf0_390[k]
                   + f_12 * lsd_234[k]
                   - f_4 * pc_x[k] * lsf1_390[k];

        t_391[k] = f_13 * lsd_186[k]
                   + f_3 * pc_y[k] * msd_234[k];

        t_392[k] = f_12 * lsd_180[k]
                   + f_3 * pc_z[k] * msd_234[k];

        t_393[k] = f_5 * lsd_237[k]
                   + f_3 * pc_x[k] * msd_237[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, t_397, pa_x, pc_x, pc_z, lsf0_396, lsd_183, \
                         lsd_238, lsd_239, lsf1_396, msd_237, msd_238, \
                         msd_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = f_5 * lsd_238[k]
                   + f_3 * pc_x[k] * msd_238[k];

        t_395[k] = f_5 * lsd_239[k]
                   + f_3 * pc_x[k] * msd_239[k];

        t_396[k] = pa_x[k] * lsf0_396[k]
                   - f_4 * pc_x[k] * lsf1_396[k];

        t_397[k] = f_12 * lsd_183[k]
                   + f_3 * pc_z[k] * msd_237[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, t_401, pa_x, pc_x, pc_y, lsf0_399, lsf0_400, \
                         lsd_191, lsd_192, lsd_240, lsf1_399, lsf1_400, msd_239, \
                         msd_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = f_13 * lsd_191[k]
                   + f_3 * pc_y[k] * msd_239[k];

        t_399[k] = pa_x[k] * lsf0_399[k]
                   - f_4 * pc_x[k] * lsf1_399[k];

        t_400[k] = pa_x[k] * lsf0_400[k]
                   + f_12 * lsd_240[k]
                   - f_4 * pc_x[k] * lsf1_400[k];

        t_401[k] = f_14 * lsd_192[k]
                   + f_3 * pc_y[k] * msd_240[k];
    }

#pragma omp simd aligned(t_402, t_403, t_404, t_405, pc_x, pc_z, lsd_186, lsd_243, lsd_244, \
                         lsd_245, msd_240, msd_243, msd_244, msd_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_402[k] = f_14 * lsd_186[k]
                   + f_3 * pc_z[k] * msd_240[k];

        t_403[k] = f_5 * lsd_243[k]
                   + f_3 * pc_x[k] * msd_243[k];

        t_404[k] = f_5 * lsd_244[k]
                   + f_3 * pc_x[k] * msd_244[k];

        t_405[k] = f_5 * lsd_245[k]
                   + f_3 * pc_x[k] * msd_245[k];
    }

#pragma omp simd aligned(t_406, t_407, t_408, t_409, pa_x, pc_x, pc_y, pc_z, lsf0_406, \
                         lsf0_409, lsd_189, lsd_197, lsf1_406, lsf1_409, msd_243, \
                         msd_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_406[k] = pa_x[k] * lsf0_406[k]
                   - f_4 * pc_x[k] * lsf1_406[k];

        t_407[k] = f_14 * lsd_189[k]
                   + f_3 * pc_z[k] * msd_243[k];

        t_408[k] = f_14 * lsd_197[k]
                   + f_3 * pc_y[k] * msd_245[k];

        t_409[k] = pa_x[k] * lsf0_409[k]
                   - f_4 * pc_x[k] * lsf1_409[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, pa_x, pc_x, pc_y, pc_z, lsf0_410, \
                         lsd_192, lsd_198, lsd_246, lsd_249, lsf1_410, msd_246, \
                         msd_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = pa_x[k] * lsf0_410[k]
                   + f_12 * lsd_246[k]
                   - f_4 * pc_x[k] * lsf1_410[k];

        t_411[k] = f_12 * lsd_198[k]
                   + f_3 * pc_y[k] * msd_246[k];

        t_412[k] = f_13 * lsd_192[k]
                   + f_3 * pc_z[k] * msd_246[k];

        t_413[k] = f_5 * lsd_249[k]
                   + f_3 * pc_x[k] * msd_249[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, t_417, pa_x, pc_x, pc_z, lsf0_416, lsd_195, \
                         lsd_250, lsd_251, lsf1_416, msd_249, msd_250, \
                         msd_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_5 * lsd_250[k]
                   + f_3 * pc_x[k] * msd_250[k];

        t_415[k] = f_5 * lsd_251[k]
                   + f_3 * pc_x[k] * msd_251[k];

        t_416[k] = pa_x[k] * lsf0_416[k]
                   - f_4 * pc_x[k] * lsf1_416[k];

        t_417[k] = f_13 * lsd_195[k]
                   + f_3 * pc_z[k] * msd_249[k];
    }

#pragma omp simd aligned(t_418, t_419, t_420, t_421, pa_x, pc_x, pc_y, lsf0_419, lsf0_420, \
                         lsd_203, lsd_204, lsd_252, lsf1_419, lsf1_420, msd_251, \
                         msd_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_418[k] = f_12 * lsd_203[k]
                   + f_3 * pc_y[k] * msd_251[k];

        t_419[k] = pa_x[k] * lsf0_419[k]
                   - f_4 * pc_x[k] * lsf1_419[k];

        t_420[k] = pa_x[k] * lsf0_420[k]
                   + f_12 * lsd_252[k]
                   - f_4 * pc_x[k] * lsf1_420[k];

        t_421[k] = f_10 * lsd_204[k]
                   + f_3 * pc_y[k] * msd_252[k];
    }

#pragma omp simd aligned(t_422, t_423, t_424, t_425, pc_x, pc_z, lsd_198, lsd_255, lsd_256, \
                         lsd_257, msd_252, msd_255, msd_256, msd_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_422[k] = f_11 * lsd_198[k]
                   + f_3 * pc_z[k] * msd_252[k];

        t_423[k] = f_5 * lsd_255[k]
                   + f_3 * pc_x[k] * msd_255[k];

        t_424[k] = f_5 * lsd_256[k]
                   + f_3 * pc_x[k] * msd_256[k];

        t_425[k] = f_5 * lsd_257[k]
                   + f_3 * pc_x[k] * msd_257[k];
    }

#pragma omp simd aligned(t_426, t_427, t_428, t_429, pa_x, pc_x, pc_y, pc_z, lsf0_426, \
                         lsf0_429, lsd_201, lsd_209, lsf1_426, lsf1_429, msd_255, \
                         msd_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_426[k] = pa_x[k] * lsf0_426[k]
                   - f_4 * pc_x[k] * lsf1_426[k];

        t_427[k] = f_11 * lsd_201[k]
                   + f_3 * pc_z[k] * msd_255[k];

        t_428[k] = f_10 * lsd_209[k]
                   + f_3 * pc_y[k] * msd_257[k];

        t_429[k] = pa_x[k] * lsf0_429[k]
                   - f_4 * pc_x[k] * lsf1_429[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, pa_y, pc_x, pc_y, pc_z, lsf0_350, \
                         lsd_204, lsd_210, lsd_261, lsf1_350, msd_258, \
                         msd_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = pa_y[k] * lsf0_350[k]
                   - f_4 * pc_y[k] * lsf1_350[k];

        t_431[k] = f_5 * lsd_210[k]
                   + f_3 * pc_y[k] * msd_258[k];

        t_432[k] = f_9 * lsd_204[k]
                   + f_3 * pc_z[k] * msd_258[k];

        t_433[k] = f_5 * lsd_261[k]
                   + f_3 * pc_x[k] * msd_261[k];
    }

#pragma omp simd aligned(t_434, t_435, t_436, t_437, pa_x, pc_x, pc_z, lsf0_436, lsd_207, \
                         lsd_262, lsd_263, lsf1_436, msd_261, msd_262, \
                         msd_263 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_434[k] = f_5 * lsd_262[k]
                   + f_3 * pc_x[k] * msd_262[k];

        t_435[k] = f_5 * lsd_263[k]
                   + f_3 * pc_x[k] * msd_263[k];

        t_436[k] = pa_x[k] * lsf0_436[k]
                   - f_4 * pc_x[k] * lsf1_436[k];

        t_437[k] = f_9 * lsd_207[k]
                   + f_3 * pc_z[k] * msd_261[k];
    }

#pragma omp simd aligned(t_438, t_439, t_440, t_441, pa_x, pc_x, pc_y, lsf0_439, lsf0_440, \
                         lsd_215, lsd_264, lsf1_439, lsf1_440, msd_263, \
                         msd_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_438[k] = f_5 * lsd_215[k]
                   + f_3 * pc_y[k] * msd_263[k];

        t_439[k] = pa_x[k] * lsf0_439[k]
                   - f_4 * pc_x[k] * lsf1_439[k];

        t_440[k] = pa_x[k] * lsf0_440[k]
                   + f_12 * lsd_264[k]
                   - f_4 * pc_x[k] * lsf1_440[k];

        t_441[k] = f_3 * pc_y[k] * msd_264[k];
    }

#pragma omp simd aligned(t_442, t_443, t_444, t_445, pc_x, pc_y, pc_z, lsd_210, lsd_267, \
                         lsd_269, msd_264, msd_266, msd_267, msd_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_442[k] = f_6 * lsd_210[k]
                   + f_3 * pc_z[k] * msd_264[k];

        t_443[k] = f_5 * lsd_267[k]
                   + f_3 * pc_x[k] * msd_267[k];

        t_444[k] = f_3 * pc_y[k] * msd_266[k];

        t_445[k] = f_5 * lsd_269[k]
                   + f_3 * pc_x[k] * msd_269[k];
    }

#pragma omp simd aligned(t_446, t_447, t_448, t_449, pa_x, pc_x, pc_y, lsf0_446, lsf0_447, \
                         lsf0_449, lsf1_446, lsf1_447, lsf1_449, \
                         msd_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_446[k] = pa_x[k] * lsf0_446[k]
                   - f_4 * pc_x[k] * lsf1_446[k];

        t_447[k] = pa_x[k] * lsf0_447[k]
                   - f_4 * pc_x[k] * lsf1_447[k];

        t_448[k] = f_3 * pc_y[k] * msd_269[k];

        t_449[k] = pa_x[k] * lsf0_449[k]
                   - f_4 * pc_x[k] * lsf1_449[k];
    }

#pragma omp simd aligned(t_450, t_451, t_452, t_453, t_454, pc_x, pc_z, msp0_135, msp0_136, \
                         msp1_135, msp1_136, msd_270, msd_271, msd_273, \
                         msd_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_450[k] = f_1 * msp0_135[k]
                   - f_2 * msp1_135[k]
                   + f_3 * pc_x[k] * msd_270[k];

        t_451[k] = f_7 * msp0_136[k]
                   - f_8 * msp1_136[k]
                   + f_3 * pc_x[k] * msd_271[k];

        t_452[k] = f_3 * pc_z[k] * msd_270[k];

        t_453[k] = f_3 * pc_x[k] * msd_273[k];

        t_454[k] = f_3 * pc_x[k] * msd_274[k];
    }

#pragma omp simd aligned(t_455, t_456, t_457, t_458, t_459, pc_x, pc_y, pc_z, lsd_219, \
                         lsd_221, msp0_136, msp0_137, msp1_136, msp1_137, msd_273, \
                         msd_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_455[k] = f_3 * pc_x[k] * msd_275[k];

        t_456[k] = f_0 * lsd_219[k]
                   + f_1 * msp0_136[k]
                   - f_2 * msp1_136[k]
                   + f_3 * pc_y[k] * msd_273[k];

        t_457[k] = f_3 * pc_z[k] * msd_273[k];

        t_458[k] = f_0 * lsd_221[k]
                   + f_3 * pc_y[k] * msd_275[k];

        t_459[k] = f_1 * msp0_137[k]
                   - f_2 * msp1_137[k]
                   + f_3 * pc_z[k] * msd_275[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, t_463, pa_z, pc_x, pc_z, lsf0_360, lsf0_361, \
                         lsf1_360, lsf1_361, msp0_140, msp1_140, msd_278, \
                         msd_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = pa_z[k] * lsf0_360[k]
                   - f_4 * pc_z[k] * lsf1_360[k];

        t_461[k] = pa_z[k] * lsf0_361[k]
                   - f_4 * pc_z[k] * lsf1_361[k];

        t_462[k] = f_7 * msp0_140[k]
                   - f_8 * msp1_140[k]
                   + f_3 * pc_x[k] * msd_278[k];

        t_463[k] = f_3 * pc_x[k] * msd_279[k];
    }

#pragma omp simd aligned(t_464, t_465, t_466, t_467, t_468, pa_z, pc_x, pc_y, pc_z, lsf0_366, \
                         lsd_219, lsd_227, lsf1_366, msd_279, msd_280, \
                         msd_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_464[k] = f_3 * pc_x[k] * msd_280[k];

        t_465[k] = f_3 * pc_x[k] * msd_281[k];

        t_466[k] = pa_z[k] * lsf0_366[k]
                   - f_4 * pc_z[k] * lsf1_366[k];

        t_467[k] = f_5 * lsd_219[k]
                   + f_3 * pc_z[k] * msd_279[k];

        t_468[k] = f_6 * lsd_227[k]
                   + f_3 * pc_y[k] * msd_281[k];
    }

#pragma omp simd aligned(t_469, t_470, t_471, pc_x, pc_z, lsd_221, msp0_140, msp0_141, \
                         msp0_142, msp1_140, msp1_141, msp1_142, msd_281, msd_282, \
                         msd_283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_469[k] = f_5 * lsd_221[k]
                   + f_1 * msp0_140[k]
                   - f_2 * msp1_140[k]
                   + f_3 * pc_z[k] * msd_281[k];

        t_470[k] = f_1 * msp0_141[k]
                   - f_2 * msp1_141[k]
                   + f_3 * pc_x[k] * msd_282[k];

        t_471[k] = f_7 * msp0_142[k]
                   - f_8 * msp1_142[k]
                   + f_3 * pc_x[k] * msd_283[k];
    }

#pragma omp simd aligned(t_472, t_473, t_474, t_475, t_476, pc_x, pc_y, lsd_231, msp0_142, \
                         msp0_143, msp1_142, msp1_143, msd_284, msd_285, msd_286, \
                         msd_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_472[k] = f_7 * msp0_143[k]
                   - f_8 * msp1_143[k]
                   + f_3 * pc_x[k] * msd_284[k];

        t_473[k] = f_3 * pc_x[k] * msd_285[k];

        t_474[k] = f_3 * pc_x[k] * msd_286[k];

        t_475[k] = f_3 * pc_x[k] * msd_287[k];

        t_476[k] = f_9 * lsd_231[k]
                   + f_1 * msp0_142[k]
                   - f_2 * msp1_142[k]
                   + f_3 * pc_y[k] * msd_285[k];
    }

#pragma omp simd aligned(t_477, t_478, t_479, pc_y, pc_z, lsd_225, lsd_227, lsd_233, msp0_143, \
                         msp1_143, msd_285, msd_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_477[k] = f_10 * lsd_225[k]
                   + f_3 * pc_z[k] * msd_285[k];

        t_478[k] = f_9 * lsd_233[k]
                   + f_3 * pc_y[k] * msd_287[k];

        t_479[k] = f_10 * lsd_227[k]
                   + f_1 * msp0_143[k]
                   - f_2 * msp1_143[k]
                   + f_3 * pc_z[k] * msd_287[k];
    }

#pragma omp simd aligned(t_480, t_481, t_482, t_483, pc_x, msp0_144, msp0_145, msp0_146, \
                         msp1_144, msp1_145, msp1_146, msd_288, msd_289, msd_290, \
                         msd_291 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_480[k] = f_1 * msp0_144[k]
                   - f_2 * msp1_144[k]
                   + f_3 * pc_x[k] * msd_288[k];

        t_481[k] = f_7 * msp0_145[k]
                   - f_8 * msp1_145[k]
                   + f_3 * pc_x[k] * msd_289[k];

        t_482[k] = f_7 * msp0_146[k]
                   - f_8 * msp1_146[k]
                   + f_3 * pc_x[k] * msd_290[k];

        t_483[k] = f_3 * pc_x[k] * msd_291[k];
    }

#pragma omp simd aligned(t_484, t_485, t_486, t_487, t_488, pc_x, pc_y, pc_z, lsd_231, \
                         lsd_237, lsd_239, msp0_145, msp1_145, msd_291, msd_292, \
                         msd_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_484[k] = f_3 * pc_x[k] * msd_292[k];

        t_485[k] = f_3 * pc_x[k] * msd_293[k];

        t_486[k] = f_11 * lsd_237[k]
                   + f_1 * msp0_145[k]
                   - f_2 * msp1_145[k]
                   + f_3 * pc_y[k] * msd_291[k];

        t_487[k] = f_12 * lsd_231[k]
                   + f_3 * pc_z[k] * msd_291[k];

        t_488[k] = f_11 * lsd_239[k]
                   + f_3 * pc_y[k] * msd_293[k];
    }

#pragma omp simd aligned(t_489, t_490, t_491, pc_x, pc_z, lsd_233, msp0_146, msp0_147, \
                         msp0_148, msp1_146, msp1_147, msp1_148, msd_293, msd_294, \
                         msd_295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_489[k] = f_12 * lsd_233[k]
                   + f_1 * msp0_146[k]
                   - f_2 * msp1_146[k]
                   + f_3 * pc_z[k] * msd_293[k];

        t_490[k] = f_1 * msp0_147[k]
                   - f_2 * msp1_147[k]
                   + f_3 * pc_x[k] * msd_294[k];

        t_491[k] = f_7 * msp0_148[k]
                   - f_8 * msp1_148[k]
                   + f_3 * pc_x[k] * msd_295[k];
    }

#pragma omp simd aligned(t_492, t_493, t_494, t_495, t_496, pc_x, pc_y, lsd_243, msp0_148, \
                         msp0_149, msp1_148, msp1_149, msd_296, msd_297, msd_298, \
                         msd_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_492[k] = f_7 * msp0_149[k]
                   - f_8 * msp1_149[k]
                   + f_3 * pc_x[k] * msd_296[k];

        t_493[k] = f_3 * pc_x[k] * msd_297[k];

        t_494[k] = f_3 * pc_x[k] * msd_298[k];

        t_495[k] = f_3 * pc_x[k] * msd_299[k];

        t_496[k] = f_13 * lsd_243[k]
                   + f_1 * msp0_148[k]
                   - f_2 * msp1_148[k]
                   + f_3 * pc_y[k] * msd_297[k];
    }

#pragma omp simd aligned(t_497, t_498, t_499, pc_y, pc_z, lsd_237, lsd_239, lsd_245, msp0_149, \
                         msp1_149, msd_297, msd_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_497[k] = f_14 * lsd_237[k]
                   + f_3 * pc_z[k] * msd_297[k];

        t_498[k] = f_13 * lsd_245[k]
                   + f_3 * pc_y[k] * msd_299[k];

        t_499[k] = f_14 * lsd_239[k]
                   + f_1 * msp0_149[k]
                   - f_2 * msp1_149[k]
                   + f_3 * pc_z[k] * msd_299[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, pc_x, msp0_150, msp0_151, msp0_152, \
                         msp1_150, msp1_151, msp1_152, msd_300, msd_301, msd_302, \
                         msd_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = f_1 * msp0_150[k]
                   - f_2 * msp1_150[k]
                   + f_3 * pc_x[k] * msd_300[k];

        t_501[k] = f_7 * msp0_151[k]
                   - f_8 * msp1_151[k]
                   + f_3 * pc_x[k] * msd_301[k];

        t_502[k] = f_7 * msp0_152[k]
                   - f_8 * msp1_152[k]
                   + f_3 * pc_x[k] * msd_302[k];

        t_503[k] = f_3 * pc_x[k] * msd_303[k];
    }

#pragma omp simd aligned(t_504, t_505, t_506, t_507, t_508, pc_x, pc_y, pc_z, lsd_243, \
                         lsd_249, lsd_251, msp0_151, msp1_151, msd_303, msd_304, \
                         msd_305 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_504[k] = f_3 * pc_x[k] * msd_304[k];

        t_505[k] = f_3 * pc_x[k] * msd_305[k];

        t_506[k] = f_14 * lsd_249[k]
                   + f_1 * msp0_151[k]
                   - f_2 * msp1_151[k]
                   + f_3 * pc_y[k] * msd_303[k];

        t_507[k] = f_13 * lsd_243[k]
                   + f_3 * pc_z[k] * msd_303[k];

        t_508[k] = f_14 * lsd_251[k]
                   + f_3 * pc_y[k] * msd_305[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, pc_x, pc_z, lsd_245, msp0_152, msp0_153, \
                         msp0_154, msp1_152, msp1_153, msp1_154, msd_305, msd_306, \
                         msd_307 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = f_13 * lsd_245[k]
                   + f_1 * msp0_152[k]
                   - f_2 * msp1_152[k]
                   + f_3 * pc_z[k] * msd_305[k];

        t_510[k] = f_1 * msp0_153[k]
                   - f_2 * msp1_153[k]
                   + f_3 * pc_x[k] * msd_306[k];

        t_511[k] = f_7 * msp0_154[k]
                   - f_8 * msp1_154[k]
                   + f_3 * pc_x[k] * msd_307[k];
    }
}

static auto
compute_prim_msf_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsf0,
                                                          const size_t lsd, const size_t lsf1,
                                                          const size_t msp0, const size_t msp1,
                                                          const size_t msd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / q;
    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 4.0 / q;
    const auto f_7 = 0.5 / gamma;
    const auto f_8 = 0.5 * p / (gamma * q);
    const auto f_9 = 3.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 3.0 / q;
    const auto f_12 = 1.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *lsf0_440 = buffer.data(lsf0 + 440);
    const auto *lsf0_442 = buffer.data(lsf0 + 442);
    const auto *lsf0_446 = buffer.data(lsf0 + 446);
    const auto *lsf0_449 = buffer.data(lsf0 + 449);

    const auto *lsd_249 = buffer.data(lsd + 249);
    const auto *lsd_251 = buffer.data(lsd + 251);
    const auto *lsd_255 = buffer.data(lsd + 255);
    const auto *lsd_257 = buffer.data(lsd + 257);
    const auto *lsd_261 = buffer.data(lsd + 261);
    const auto *lsd_263 = buffer.data(lsd + 263);
    const auto *lsd_267 = buffer.data(lsd + 267);
    const auto *lsd_269 = buffer.data(lsd + 269);

    const auto *lsf1_440 = buffer.data(lsf1 + 440);
    const auto *lsf1_442 = buffer.data(lsf1 + 442);
    const auto *lsf1_446 = buffer.data(lsf1 + 446);
    const auto *lsf1_449 = buffer.data(lsf1 + 449);

    const auto *msp0_154 = buffer.data(msp0 + 154);
    const auto *msp0_155 = buffer.data(msp0 + 155);
    const auto *msp0_156 = buffer.data(msp0 + 156);
    const auto *msp0_157 = buffer.data(msp0 + 157);
    const auto *msp0_158 = buffer.data(msp0 + 158);
    const auto *msp0_160 = buffer.data(msp0 + 160);
    const auto *msp0_162 = buffer.data(msp0 + 162);
    const auto *msp0_163 = buffer.data(msp0 + 163);
    const auto *msp0_164 = buffer.data(msp0 + 164);

    const auto *msp1_154 = buffer.data(msp1 + 154);
    const auto *msp1_155 = buffer.data(msp1 + 155);
    const auto *msp1_156 = buffer.data(msp1 + 156);
    const auto *msp1_157 = buffer.data(msp1 + 157);
    const auto *msp1_158 = buffer.data(msp1 + 158);
    const auto *msp1_160 = buffer.data(msp1 + 160);
    const auto *msp1_162 = buffer.data(msp1 + 162);
    const auto *msp1_163 = buffer.data(msp1 + 163);
    const auto *msp1_164 = buffer.data(msp1 + 164);

    const auto *msd_308 = buffer.data(msd + 308);
    const auto *msd_309 = buffer.data(msd + 309);
    const auto *msd_310 = buffer.data(msd + 310);
    const auto *msd_311 = buffer.data(msd + 311);
    const auto *msd_312 = buffer.data(msd + 312);
    const auto *msd_313 = buffer.data(msd + 313);
    const auto *msd_314 = buffer.data(msd + 314);
    const auto *msd_315 = buffer.data(msd + 315);
    const auto *msd_316 = buffer.data(msd + 316);
    const auto *msd_317 = buffer.data(msd + 317);
    const auto *msd_319 = buffer.data(msd + 319);
    const auto *msd_321 = buffer.data(msd + 321);
    const auto *msd_322 = buffer.data(msd + 322);
    const auto *msd_323 = buffer.data(msd + 323);
    const auto *msd_324 = buffer.data(msd + 324);
    const auto *msd_326 = buffer.data(msd + 326);
    const auto *msd_327 = buffer.data(msd + 327);
    const auto *msd_328 = buffer.data(msd + 328);
    const auto *msd_329 = buffer.data(msd + 329);

#pragma omp simd aligned(t_512, t_513, t_514, t_515, t_516, pc_x, pc_y, lsd_255, msp0_154, \
                         msp0_155, msp1_154, msp1_155, msd_308, msd_309, msd_310, \
                         msd_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = f_7 * msp0_155[k]
                   - f_8 * msp1_155[k]
                   + f_3 * pc_x[k] * msd_308[k];

        t_513[k] = f_3 * pc_x[k] * msd_309[k];

        t_514[k] = f_3 * pc_x[k] * msd_310[k];

        t_515[k] = f_3 * pc_x[k] * msd_311[k];

        t_516[k] = f_12 * lsd_255[k]
                   + f_1 * msp0_154[k]
                   - f_2 * msp1_154[k]
                   + f_3 * pc_y[k] * msd_309[k];
    }

#pragma omp simd aligned(t_517, t_518, t_519, pc_y, pc_z, lsd_249, lsd_251, lsd_257, msp0_155, \
                         msp1_155, msd_309, msd_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_517[k] = f_11 * lsd_249[k]
                   + f_3 * pc_z[k] * msd_309[k];

        t_518[k] = f_12 * lsd_257[k]
                   + f_3 * pc_y[k] * msd_311[k];

        t_519[k] = f_11 * lsd_251[k]
                   + f_1 * msp0_155[k]
                   - f_2 * msp1_155[k]
                   + f_3 * pc_z[k] * msd_311[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, t_523, pc_x, msp0_156, msp0_157, msp0_158, \
                         msp1_156, msp1_157, msp1_158, msd_312, msd_313, msd_314, \
                         msd_315 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = f_1 * msp0_156[k]
                   - f_2 * msp1_156[k]
                   + f_3 * pc_x[k] * msd_312[k];

        t_521[k] = f_7 * msp0_157[k]
                   - f_8 * msp1_157[k]
                   + f_3 * pc_x[k] * msd_313[k];

        t_522[k] = f_7 * msp0_158[k]
                   - f_8 * msp1_158[k]
                   + f_3 * pc_x[k] * msd_314[k];

        t_523[k] = f_3 * pc_x[k] * msd_315[k];
    }

#pragma omp simd aligned(t_524, t_525, t_526, t_527, t_528, pc_x, pc_y, pc_z, lsd_255, \
                         lsd_261, lsd_263, msp0_157, msp1_157, msd_315, msd_316, \
                         msd_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_524[k] = f_3 * pc_x[k] * msd_316[k];

        t_525[k] = f_3 * pc_x[k] * msd_317[k];

        t_526[k] = f_10 * lsd_261[k]
                   + f_1 * msp0_157[k]
                   - f_2 * msp1_157[k]
                   + f_3 * pc_y[k] * msd_315[k];

        t_527[k] = f_9 * lsd_255[k]
                   + f_3 * pc_z[k] * msd_315[k];

        t_528[k] = f_10 * lsd_263[k]
                   + f_3 * pc_y[k] * msd_317[k];
    }

#pragma omp simd aligned(t_529, t_530, t_531, pa_y, pc_x, pc_y, pc_z, lsf0_440, lsd_257, \
                         lsf1_440, msp0_158, msp0_160, msp1_158, msp1_160, msd_317, \
                         msd_319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_529[k] = f_9 * lsd_257[k]
                   + f_1 * msp0_158[k]
                   - f_2 * msp1_158[k]
                   + f_3 * pc_z[k] * msd_317[k];

        t_530[k] = pa_y[k] * lsf0_440[k]
                   - f_4 * pc_y[k] * lsf1_440[k];

        t_531[k] = f_7 * msp0_160[k]
                   - f_8 * msp1_160[k]
                   + f_3 * pc_x[k] * msd_319[k];
    }

#pragma omp simd aligned(t_532, t_533, t_534, t_535, t_536, pa_y, pc_x, pc_y, lsf0_442, \
                         lsf0_446, lsd_267, lsf1_442, lsf1_446, msd_321, msd_322, \
                         msd_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_532[k] = pa_y[k] * lsf0_442[k]
                   - f_4 * pc_y[k] * lsf1_442[k];

        t_533[k] = f_3 * pc_x[k] * msd_321[k];

        t_534[k] = f_3 * pc_x[k] * msd_322[k];

        t_535[k] = f_3 * pc_x[k] * msd_323[k];

        t_536[k] = pa_y[k] * lsf0_446[k]
                   + f_12 * lsd_267[k]
                   - f_4 * pc_y[k] * lsf1_446[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, pa_y, pc_y, pc_z, lsf0_449, lsd_261, lsd_269, \
                         lsf1_449, msd_321, msd_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = f_6 * lsd_261[k]
                   + f_3 * pc_z[k] * msd_321[k];

        t_538[k] = f_5 * lsd_269[k]
                   + f_3 * pc_y[k] * msd_323[k];

        t_539[k] = pa_y[k] * lsf0_449[k]
                   - f_4 * pc_y[k] * lsf1_449[k];
    }

#pragma omp simd aligned(t_540, t_541, t_542, t_543, t_544, pc_x, pc_y, msp0_162, msp0_164, \
                         msp1_162, msp1_164, msd_324, msd_326, msd_327, \
                         msd_328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_540[k] = f_1 * msp0_162[k]
                   - f_2 * msp1_162[k]
                   + f_3 * pc_x[k] * msd_324[k];

        t_541[k] = f_3 * pc_y[k] * msd_324[k];

        t_542[k] = f_7 * msp0_164[k]
                   - f_8 * msp1_164[k]
                   + f_3 * pc_x[k] * msd_326[k];

        t_543[k] = f_3 * pc_x[k] * msd_327[k];

        t_544[k] = f_3 * pc_x[k] * msd_328[k];
    }

#pragma omp simd aligned(t_545, t_546, t_547, t_548, t_549, pc_x, pc_y, pc_z, lsd_269, \
                         msp0_163, msp0_164, msp1_163, msp1_164, msd_327, msd_328, \
                         msd_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_545[k] = f_3 * pc_x[k] * msd_329[k];

        t_546[k] = f_1 * msp0_163[k]
                   - f_2 * msp1_163[k]
                   + f_3 * pc_y[k] * msd_327[k];

        t_547[k] = f_7 * msp0_164[k]
                   - f_8 * msp1_164[k]
                   + f_3 * pc_y[k] * msd_328[k];

        t_548[k] = f_3 * pc_y[k] * msd_329[k];

        t_549[k] = f_0 * lsd_269[k]
                   + f_1 * msp0_164[k]
                   - f_2 * msp1_164[k]
                   + f_3 * pc_z[k] * msd_329[k];
    }
}

auto
compute_prim_msf_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t lsf0, const size_t lsd,
                                                   const size_t lsf1, const size_t msp0,
                                                   const size_t msp1, const size_t msd,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_msf_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, lsf0, lsd,
                                                              lsf1, msp0, msp1, msd, ncols,
                                                              gamma, p, q);

    compute_prim_msf_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, lsf0, lsd,
                                                              lsf1, msp0, msp1, msd, ncols,
                                                              gamma, p, q);

    compute_prim_msf_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, lsf0, lsd,
                                                              lsf1, msp0, msp1, msd, ncols,
                                                              gamma, p, q);

    compute_prim_msf_three_center_electron_repulsion_0_piece3(buffer, target, pa, pc, lsf0, lsd,
                                                              lsf1, msp0, msp1, msd, ncols,
                                                              gamma, p, q);

    compute_prim_msf_three_center_electron_repulsion_0_piece4(buffer, target, pa, pc, lsf0, lsd,
                                                              lsf1, msp0, msp1, msd, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
