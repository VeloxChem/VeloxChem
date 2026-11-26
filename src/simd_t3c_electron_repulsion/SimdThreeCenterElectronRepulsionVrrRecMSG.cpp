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


#include "SimdThreeCenterElectronRepulsionVrrRecMSG.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_msg_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsg0,
                                                          const size_t lsf, const size_t lsg1,
                                                          const size_t msd0, const size_t msd1,
                                                          const size_t msf, const size_t ncols,
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
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 3.5 / q;
    const auto f_13 = 3.0 / q;
    const auto f_14 = 1.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *lsg0_0 = buffer.data(lsg0 + 0);
    const auto *lsg0_3 = buffer.data(lsg0 + 3);
    const auto *lsg0_5 = buffer.data(lsg0 + 5);
    const auto *lsg0_10 = buffer.data(lsg0 + 10);
    const auto *lsg0_14 = buffer.data(lsg0 + 14);
    const auto *lsg0_18 = buffer.data(lsg0 + 18);
    const auto *lsg0_25 = buffer.data(lsg0 + 25);
    const auto *lsg0_30 = buffer.data(lsg0 + 30);
    const auto *lsg0_35 = buffer.data(lsg0 + 35);
    const auto *lsg0_44 = buffer.data(lsg0 + 44);
    const auto *lsg0_45 = buffer.data(lsg0 + 45);
    const auto *lsg0_48 = buffer.data(lsg0 + 48);
    const auto *lsg0_55 = buffer.data(lsg0 + 55);
    const auto *lsg0_75 = buffer.data(lsg0 + 75);
    const auto *lsg0_78 = buffer.data(lsg0 + 78);
    const auto *lsg0_80 = buffer.data(lsg0 + 80);

    const auto *lsf_0 = buffer.data(lsf + 0);
    const auto *lsf_1 = buffer.data(lsf + 1);
    const auto *lsf_2 = buffer.data(lsf + 2);
    const auto *lsf_6 = buffer.data(lsf + 6);
    const auto *lsf_9 = buffer.data(lsf + 9);
    const auto *lsf_10 = buffer.data(lsf + 10);
    const auto *lsf_16 = buffer.data(lsf + 16);
    const auto *lsf_18 = buffer.data(lsf + 18);
    const auto *lsf_19 = buffer.data(lsf + 19);
    const auto *lsf_20 = buffer.data(lsf + 20);
    const auto *lsf_22 = buffer.data(lsf + 22);
    const auto *lsf_26 = buffer.data(lsf + 26);
    const auto *lsf_27 = buffer.data(lsf + 27);
    const auto *lsf_28 = buffer.data(lsf + 28);
    const auto *lsf_29 = buffer.data(lsf + 29);
    const auto *lsf_30 = buffer.data(lsf + 30);
    const auto *lsf_33 = buffer.data(lsf + 33);
    const auto *lsf_36 = buffer.data(lsf + 36);
    const auto *lsf_38 = buffer.data(lsf + 38);
    const auto *lsf_39 = buffer.data(lsf + 39);
    const auto *lsf_40 = buffer.data(lsf + 40);
    const auto *lsf_42 = buffer.data(lsf + 42);
    const auto *lsf_46 = buffer.data(lsf + 46);
    const auto *lsf_47 = buffer.data(lsf + 47);
    const auto *lsf_48 = buffer.data(lsf + 48);
    const auto *lsf_49 = buffer.data(lsf + 49);
    const auto *lsf_50 = buffer.data(lsf + 50);
    const auto *lsf_51 = buffer.data(lsf + 51);
    const auto *lsf_52 = buffer.data(lsf + 52);
    const auto *lsf_55 = buffer.data(lsf + 55);
    const auto *lsf_56 = buffer.data(lsf + 56);
    const auto *lsf_57 = buffer.data(lsf + 57);
    const auto *lsf_59 = buffer.data(lsf + 59);
    const auto *lsf_60 = buffer.data(lsf + 60);
    const auto *lsf_63 = buffer.data(lsf + 63);
    const auto *lsf_66 = buffer.data(lsf + 66);
    const auto *lsf_68 = buffer.data(lsf + 68);
    const auto *lsf_69 = buffer.data(lsf + 69);
    const auto *lsf_75 = buffer.data(lsf + 75);
    const auto *lsf_76 = buffer.data(lsf + 76);
    const auto *lsf_77 = buffer.data(lsf + 77);
    const auto *lsf_78 = buffer.data(lsf + 78);
    const auto *lsf_79 = buffer.data(lsf + 79);
    const auto *lsf_86 = buffer.data(lsf + 86);
    const auto *lsf_87 = buffer.data(lsf + 87);
    const auto *lsf_88 = buffer.data(lsf + 88);
    const auto *lsf_89 = buffer.data(lsf + 89);

    const auto *lsg1_0 = buffer.data(lsg1 + 0);
    const auto *lsg1_3 = buffer.data(lsg1 + 3);
    const auto *lsg1_5 = buffer.data(lsg1 + 5);
    const auto *lsg1_10 = buffer.data(lsg1 + 10);
    const auto *lsg1_14 = buffer.data(lsg1 + 14);
    const auto *lsg1_18 = buffer.data(lsg1 + 18);
    const auto *lsg1_25 = buffer.data(lsg1 + 25);
    const auto *lsg1_30 = buffer.data(lsg1 + 30);
    const auto *lsg1_35 = buffer.data(lsg1 + 35);
    const auto *lsg1_44 = buffer.data(lsg1 + 44);
    const auto *lsg1_45 = buffer.data(lsg1 + 45);
    const auto *lsg1_48 = buffer.data(lsg1 + 48);
    const auto *lsg1_55 = buffer.data(lsg1 + 55);
    const auto *lsg1_75 = buffer.data(lsg1 + 75);
    const auto *lsg1_78 = buffer.data(lsg1 + 78);
    const auto *lsg1_80 = buffer.data(lsg1 + 80);

    const auto *msd0_0 = buffer.data(msd0 + 0);
    const auto *msd0_3 = buffer.data(msd0 + 3);
    const auto *msd0_5 = buffer.data(msd0 + 5);
    const auto *msd0_9 = buffer.data(msd0 + 9);
    const auto *msd0_16 = buffer.data(msd0 + 16);
    const auto *msd0_17 = buffer.data(msd0 + 17);
    const auto *msd0_18 = buffer.data(msd0 + 18);
    const auto *msd0_21 = buffer.data(msd0 + 21);
    const auto *msd0_23 = buffer.data(msd0 + 23);
    const auto *msd0_29 = buffer.data(msd0 + 29);
    const auto *msd0_30 = buffer.data(msd0 + 30);
    const auto *msd0_33 = buffer.data(msd0 + 33);
    const auto *msd0_34 = buffer.data(msd0 + 34);
    const auto *msd0_35 = buffer.data(msd0 + 35);
    const auto *msd0_36 = buffer.data(msd0 + 36);
    const auto *msd0_39 = buffer.data(msd0 + 39);
    const auto *msd0_41 = buffer.data(msd0 + 41);
    const auto *msd0_47 = buffer.data(msd0 + 47);

    const auto *msd1_0 = buffer.data(msd1 + 0);
    const auto *msd1_3 = buffer.data(msd1 + 3);
    const auto *msd1_5 = buffer.data(msd1 + 5);
    const auto *msd1_9 = buffer.data(msd1 + 9);
    const auto *msd1_16 = buffer.data(msd1 + 16);
    const auto *msd1_17 = buffer.data(msd1 + 17);
    const auto *msd1_18 = buffer.data(msd1 + 18);
    const auto *msd1_21 = buffer.data(msd1 + 21);
    const auto *msd1_23 = buffer.data(msd1 + 23);
    const auto *msd1_29 = buffer.data(msd1 + 29);
    const auto *msd1_30 = buffer.data(msd1 + 30);
    const auto *msd1_33 = buffer.data(msd1 + 33);
    const auto *msd1_34 = buffer.data(msd1 + 34);
    const auto *msd1_35 = buffer.data(msd1 + 35);
    const auto *msd1_36 = buffer.data(msd1 + 36);
    const auto *msd1_39 = buffer.data(msd1 + 39);
    const auto *msd1_41 = buffer.data(msd1 + 41);
    const auto *msd1_47 = buffer.data(msd1 + 47);

    const auto *msf_0 = buffer.data(msf + 0);
    const auto *msf_1 = buffer.data(msf + 1);
    const auto *msf_2 = buffer.data(msf + 2);
    const auto *msf_3 = buffer.data(msf + 3);
    const auto *msf_5 = buffer.data(msf + 5);
    const auto *msf_6 = buffer.data(msf + 6);
    const auto *msf_8 = buffer.data(msf + 8);
    const auto *msf_9 = buffer.data(msf + 9);
    const auto *msf_10 = buffer.data(msf + 10);
    const auto *msf_11 = buffer.data(msf + 11);
    const auto *msf_13 = buffer.data(msf + 13);
    const auto *msf_16 = buffer.data(msf + 16);
    const auto *msf_17 = buffer.data(msf + 17);
    const auto *msf_18 = buffer.data(msf + 18);
    const auto *msf_19 = buffer.data(msf + 19);
    const auto *msf_20 = buffer.data(msf + 20);
    const auto *msf_22 = buffer.data(msf + 22);
    const auto *msf_25 = buffer.data(msf + 25);
    const auto *msf_26 = buffer.data(msf + 26);
    const auto *msf_27 = buffer.data(msf + 27);
    const auto *msf_28 = buffer.data(msf + 28);
    const auto *msf_29 = buffer.data(msf + 29);
    const auto *msf_30 = buffer.data(msf + 30);
    const auto *msf_31 = buffer.data(msf + 31);
    const auto *msf_32 = buffer.data(msf + 32);
    const auto *msf_33 = buffer.data(msf + 33);
    const auto *msf_36 = buffer.data(msf + 36);
    const auto *msf_37 = buffer.data(msf + 37);
    const auto *msf_38 = buffer.data(msf + 38);
    const auto *msf_39 = buffer.data(msf + 39);
    const auto *msf_40 = buffer.data(msf + 40);
    const auto *msf_42 = buffer.data(msf + 42);
    const auto *msf_46 = buffer.data(msf + 46);
    const auto *msf_47 = buffer.data(msf + 47);
    const auto *msf_48 = buffer.data(msf + 48);
    const auto *msf_49 = buffer.data(msf + 49);
    const auto *msf_50 = buffer.data(msf + 50);
    const auto *msf_51 = buffer.data(msf + 51);
    const auto *msf_52 = buffer.data(msf + 52);
    const auto *msf_55 = buffer.data(msf + 55);
    const auto *msf_56 = buffer.data(msf + 56);
    const auto *msf_57 = buffer.data(msf + 57);
    const auto *msf_58 = buffer.data(msf + 58);
    const auto *msf_59 = buffer.data(msf + 59);
    const auto *msf_60 = buffer.data(msf + 60);
    const auto *msf_61 = buffer.data(msf + 61);
    const auto *msf_62 = buffer.data(msf + 62);
    const auto *msf_63 = buffer.data(msf + 63);
    const auto *msf_66 = buffer.data(msf + 66);
    const auto *msf_67 = buffer.data(msf + 67);
    const auto *msf_68 = buffer.data(msf + 68);
    const auto *msf_69 = buffer.data(msf + 69);
    const auto *msf_70 = buffer.data(msf + 70);
    const auto *msf_72 = buffer.data(msf + 72);
    const auto *msf_75 = buffer.data(msf + 75);
    const auto *msf_76 = buffer.data(msf + 76);
    const auto *msf_77 = buffer.data(msf + 77);
    const auto *msf_78 = buffer.data(msf + 78);
    const auto *msf_79 = buffer.data(msf + 79);
    const auto *msf_80 = buffer.data(msf + 80);
    const auto *msf_82 = buffer.data(msf + 82);
    const auto *msf_86 = buffer.data(msf + 86);
    const auto *msf_87 = buffer.data(msf + 87);
    const auto *msf_88 = buffer.data(msf + 88);
    const auto *msf_89 = buffer.data(msf + 89);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, lsf_0, msd0_0, \
                         msd1_0, msf_0, msf_1, msf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * lsf_0[k]
                 + f_1 * msd0_0[k]
                 - f_2 * msd1_0[k]
                 + f_3 * pc_x[k] * msf_0[k];

        t_1[k] = f_3 * pc_y[k] * msf_0[k];

        t_2[k] = f_3 * pc_z[k] * msf_0[k];

        t_3[k] = f_4 * msd0_0[k]
                 - f_5 * msd1_0[k]
                 + f_3 * pc_y[k] * msf_1[k];

        t_4[k] = f_3 * pc_y[k] * msf_2[k];

        t_5[k] = f_4 * msd0_0[k]
                 - f_5 * msd1_0[k]
                 + f_3 * pc_z[k] * msf_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_x, pc_y, pc_z, lsf_6, lsf_9, msd0_3, \
                         msd1_3, msf_3, msf_5, msf_6, msf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * lsf_6[k]
                 + f_3 * pc_x[k] * msf_6[k];

        t_7[k] = f_3 * pc_z[k] * msf_3[k];

        t_8[k] = f_3 * pc_y[k] * msf_5[k];

        t_9[k] = f_0 * lsf_9[k]
                 + f_3 * pc_x[k] * msf_9[k];

        t_10[k] = f_1 * msd0_3[k]
                  - f_2 * msd1_3[k]
                  + f_3 * pc_y[k] * msf_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pa_y, pc_y, pc_z, lsg0_0, lsg1_0, \
                         msd0_5, msd1_5, msf_6, msf_8, msf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * pc_z[k] * msf_6[k];

        t_12[k] = f_4 * msd0_5[k]
                  - f_5 * msd1_5[k]
                  + f_3 * pc_y[k] * msf_8[k];

        t_13[k] = f_3 * pc_y[k] * msf_9[k];

        t_14[k] = f_1 * msd0_5[k]
                  - f_2 * msd1_5[k]
                  + f_3 * pc_z[k] * msf_9[k];

        t_15[k] = pa_y[k] * lsg0_0[k]
                  - f_6 * pc_y[k] * lsg1_0[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pa_y, pc_y, pc_z, lsg0_3, lsg0_5, \
                         lsf_0, lsf_1, lsg1_3, lsg1_5, msf_10, msf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_7 * lsf_0[k]
                  + f_3 * pc_y[k] * msf_10[k];

        t_17[k] = f_3 * pc_z[k] * msf_10[k];

        t_18[k] = pa_y[k] * lsg0_3[k]
                  + f_8 * lsf_1[k]
                  - f_6 * pc_y[k] * lsg1_3[k];

        t_19[k] = f_3 * pc_z[k] * msf_11[k];

        t_20[k] = pa_y[k] * lsg0_5[k]
                  - f_6 * pc_y[k] * lsg1_5[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pc_x, pc_z, lsf_16, lsf_18, lsf_19, msf_13, \
                         msf_16, msf_18, msf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_9 * lsf_16[k]
                  + f_3 * pc_x[k] * msf_16[k];

        t_22[k] = f_3 * pc_z[k] * msf_13[k];

        t_23[k] = f_9 * lsf_18[k]
                  + f_3 * pc_x[k] * msf_18[k];

        t_24[k] = f_9 * lsf_19[k]
                  + f_3 * pc_x[k] * msf_19[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pc_y, pc_z, lsf_6, lsf_9, msd0_9, msd1_9, \
                         msf_16, msf_17, msf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_7 * lsf_6[k]
                  + f_1 * msd0_9[k]
                  - f_2 * msd1_9[k]
                  + f_3 * pc_y[k] * msf_16[k];

        t_26[k] = f_3 * pc_z[k] * msf_16[k];

        t_27[k] = f_4 * msd0_9[k]
                  - f_5 * msd1_9[k]
                  + f_3 * pc_z[k] * msf_17[k];

        t_28[k] = f_7 * lsf_9[k]
                  + f_3 * pc_y[k] * msf_19[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pa_y, pa_z, pc_y, pc_z, lsg0_0, lsg0_14, \
                         lsf_0, lsg1_0, lsg1_14, msf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pa_y[k] * lsg0_14[k]
                  - f_6 * pc_y[k] * lsg1_14[k];

        t_30[k] = pa_z[k] * lsg0_0[k]
                  - f_6 * pc_z[k] * lsg1_0[k];

        t_31[k] = f_3 * pc_y[k] * msf_20[k];

        t_32[k] = f_7 * lsf_0[k]
                  + f_3 * pc_z[k] * msf_20[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pa_z, pc_x, pc_y, pc_z, lsg0_3, lsg0_5, \
                         lsf_2, lsf_26, lsg1_3, lsg1_5, msf_22, \
                         msf_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pa_z[k] * lsg0_3[k]
                  - f_6 * pc_z[k] * lsg1_3[k];

        t_34[k] = f_3 * pc_y[k] * msf_22[k];

        t_35[k] = pa_z[k] * lsg0_5[k]
                  + f_8 * lsf_2[k]
                  - f_6 * pc_z[k] * lsg1_5[k];

        t_36[k] = f_9 * lsf_26[k]
                  + f_3 * pc_x[k] * msf_26[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pa_z, pc_x, pc_y, pc_z, lsg0_10, lsf_27, \
                         lsf_29, lsg1_10, msf_25, msf_27, msf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_9 * lsf_27[k]
                  + f_3 * pc_x[k] * msf_27[k];

        t_38[k] = f_3 * pc_y[k] * msf_25[k];

        t_39[k] = f_9 * lsf_29[k]
                  + f_3 * pc_x[k] * msf_29[k];

        t_40[k] = pa_z[k] * lsg0_10[k]
                  - f_6 * pc_z[k] * lsg1_10[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pc_y, pc_z, lsf_9, msd0_16, msd0_17, msd1_16, \
                         msd1_17, msf_27, msf_28, msf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_10 * msd0_16[k]
                  - f_11 * msd1_16[k]
                  + f_3 * pc_y[k] * msf_27[k];

        t_42[k] = f_4 * msd0_17[k]
                  - f_5 * msd1_17[k]
                  + f_3 * pc_y[k] * msf_28[k];

        t_43[k] = f_3 * pc_y[k] * msf_29[k];

        t_44[k] = f_7 * lsf_9[k]
                  + f_1 * msd0_17[k]
                  - f_2 * msd1_17[k]
                  + f_3 * pc_z[k] * msf_29[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pc_x, pc_y, pc_z, lsf_10, lsf_30, lsf_33, \
                         msd0_18, msd0_21, msd1_18, msd1_21, msf_30, \
                         msf_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_12 * lsf_30[k]
                  + f_1 * msd0_18[k]
                  - f_2 * msd1_18[k]
                  + f_3 * pc_x[k] * msf_30[k];

        t_46[k] = f_8 * lsf_10[k]
                  + f_3 * pc_y[k] * msf_30[k];

        t_47[k] = f_3 * pc_z[k] * msf_30[k];

        t_48[k] = f_12 * lsf_33[k]
                  + f_4 * msd0_21[k]
                  - f_5 * msd1_21[k]
                  + f_3 * pc_x[k] * msf_33[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, pc_x, pc_z, lsf_36, lsf_38, msd0_18, \
                         msd1_18, msf_31, msf_32, msf_33, msf_36, \
                         msf_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_3 * pc_z[k] * msf_31[k];

        t_50[k] = f_4 * msd0_18[k]
                  - f_5 * msd1_18[k]
                  + f_3 * pc_z[k] * msf_32[k];

        t_51[k] = f_12 * lsf_36[k]
                  + f_3 * pc_x[k] * msf_36[k];

        t_52[k] = f_3 * pc_z[k] * msf_33[k];

        t_53[k] = f_12 * lsf_38[k]
                  + f_3 * pc_x[k] * msf_38[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, pc_x, pc_y, pc_z, lsf_16, lsf_19, \
                         lsf_39, msd0_21, msd1_21, msf_36, msf_37, \
                         msf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_12 * lsf_39[k]
                  + f_3 * pc_x[k] * msf_39[k];

        t_55[k] = f_8 * lsf_16[k]
                  + f_1 * msd0_21[k]
                  - f_2 * msd1_21[k]
                  + f_3 * pc_y[k] * msf_36[k];

        t_56[k] = f_3 * pc_z[k] * msf_36[k];

        t_57[k] = f_4 * msd0_21[k]
                  - f_5 * msd1_21[k]
                  + f_3 * pc_z[k] * msf_37[k];

        t_58[k] = f_8 * lsf_19[k]
                  + f_3 * pc_y[k] * msf_39[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pa_y, pc_y, pc_z, lsg0_30, lsf_10, lsf_20, \
                         lsg1_30, msd0_23, msd1_23, msf_39, msf_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_1 * msd0_23[k]
                  - f_2 * msd1_23[k]
                  + f_3 * pc_z[k] * msf_39[k];

        t_60[k] = pa_y[k] * lsg0_30[k]
                  - f_6 * pc_y[k] * lsg1_30[k];

        t_61[k] = f_7 * lsf_20[k]
                  + f_3 * pc_y[k] * msf_40[k];

        t_62[k] = f_7 * lsf_10[k]
                  + f_3 * pc_z[k] * msf_40[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, pa_y, pa_z, pc_y, pc_z, lsg0_18, lsg0_35, lsf_22, \
                         lsg1_18, lsg1_35, msf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = pa_z[k] * lsg0_18[k]
                  - f_6 * pc_z[k] * lsg1_18[k];

        t_64[k] = f_7 * lsf_22[k]
                  + f_3 * pc_y[k] * msf_42[k];

        t_65[k] = pa_y[k] * lsg0_35[k]
                  - f_6 * pc_y[k] * lsg1_35[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pc_x, lsf_46, lsf_47, lsf_48, lsf_49, msf_46, \
                         msf_47, msf_48, msf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_12 * lsf_46[k]
                  + f_3 * pc_x[k] * msf_46[k];

        t_67[k] = f_12 * lsf_47[k]
                  + f_3 * pc_x[k] * msf_47[k];

        t_68[k] = f_12 * lsf_48[k]
                  + f_3 * pc_x[k] * msf_48[k];

        t_69[k] = f_12 * lsf_49[k]
                  + f_3 * pc_x[k] * msf_49[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, pa_z, pc_y, pc_z, lsg0_25, lsf_16, lsf_28, lsg1_25, \
                         msd0_29, msd1_29, msf_46, msf_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = pa_z[k] * lsg0_25[k]
                  - f_6 * pc_z[k] * lsg1_25[k];

        t_71[k] = f_7 * lsf_16[k]
                  + f_3 * pc_z[k] * msf_46[k];

        t_72[k] = f_7 * lsf_28[k]
                  + f_4 * msd0_29[k]
                  - f_5 * msd1_29[k]
                  + f_3 * pc_y[k] * msf_48[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pa_y, pc_x, pc_y, lsg0_44, lsf_29, lsf_50, \
                         lsg1_44, msd0_30, msd1_30, msf_49, msf_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_7 * lsf_29[k]
                  + f_3 * pc_y[k] * msf_49[k];

        t_74[k] = pa_y[k] * lsg0_44[k]
                  - f_6 * pc_y[k] * lsg1_44[k];

        t_75[k] = f_12 * lsf_50[k]
                  + f_1 * msd0_30[k]
                  - f_2 * msd1_30[k]
                  + f_3 * pc_x[k] * msf_50[k];

        t_76[k] = f_3 * pc_y[k] * msf_50[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, pc_y, pc_z, lsf_20, msd0_30, msd1_30, msf_50, \
                         msf_51, msf_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_8 * lsf_20[k]
                  + f_3 * pc_z[k] * msf_50[k];

        t_78[k] = f_4 * msd0_30[k]
                  - f_5 * msd1_30[k]
                  + f_3 * pc_y[k] * msf_51[k];

        t_79[k] = f_3 * pc_y[k] * msf_52[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pc_x, pc_y, lsf_55, lsf_56, lsf_57, msd0_35, \
                         msd1_35, msf_55, msf_56, msf_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_12 * lsf_55[k]
                  + f_4 * msd0_35[k]
                  - f_5 * msd1_35[k]
                  + f_3 * pc_x[k] * msf_55[k];

        t_81[k] = f_12 * lsf_56[k]
                  + f_3 * pc_x[k] * msf_56[k];

        t_82[k] = f_12 * lsf_57[k]
                  + f_3 * pc_x[k] * msf_57[k];

        t_83[k] = f_3 * pc_y[k] * msf_55[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, pc_x, pc_y, lsf_59, msd0_33, msd0_34, msd1_33, \
                         msd1_34, msf_56, msf_57, msf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_12 * lsf_59[k]
                  + f_3 * pc_x[k] * msf_59[k];

        t_85[k] = f_1 * msd0_33[k]
                  - f_2 * msd1_33[k]
                  + f_3 * pc_y[k] * msf_56[k];

        t_86[k] = f_10 * msd0_34[k]
                  - f_11 * msd1_34[k]
                  + f_3 * pc_y[k] * msf_57[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pc_x, pc_y, pc_z, lsf_29, lsf_60, msd0_35, \
                         msd0_36, msd1_35, msd1_36, msf_58, msf_59, \
                         msf_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_4 * msd0_35[k]
                  - f_5 * msd1_35[k]
                  + f_3 * pc_y[k] * msf_58[k];

        t_88[k] = f_3 * pc_y[k] * msf_59[k];

        t_89[k] = f_8 * lsf_29[k]
                  + f_1 * msd0_35[k]
                  - f_2 * msd1_35[k]
                  + f_3 * pc_z[k] * msf_59[k];

        t_90[k] = f_13 * lsf_60[k]
                  + f_1 * msd0_36[k]
                  - f_2 * msd1_36[k]
                  + f_3 * pc_x[k] * msf_60[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, pc_x, pc_y, pc_z, lsf_30, lsf_63, msd0_39, \
                         msd1_39, msf_60, msf_61, msf_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_14 * lsf_30[k]
                  + f_3 * pc_y[k] * msf_60[k];

        t_92[k] = f_3 * pc_z[k] * msf_60[k];

        t_93[k] = f_13 * lsf_63[k]
                  + f_4 * msd0_39[k]
                  - f_5 * msd1_39[k]
                  + f_3 * pc_x[k] * msf_63[k];

        t_94[k] = f_3 * pc_z[k] * msf_61[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pc_x, pc_z, lsf_66, lsf_68, msd0_36, msd1_36, \
                         msf_62, msf_63, msf_66, msf_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_4 * msd0_36[k]
                  - f_5 * msd1_36[k]
                  + f_3 * pc_z[k] * msf_62[k];

        t_96[k] = f_13 * lsf_66[k]
                  + f_3 * pc_x[k] * msf_66[k];

        t_97[k] = f_3 * pc_z[k] * msf_63[k];

        t_98[k] = f_13 * lsf_68[k]
                  + f_3 * pc_x[k] * msf_68[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, pc_x, pc_y, pc_z, lsf_36, lsf_39, \
                         lsf_69, msd0_39, msd1_39, msf_66, msf_67, \
                         msf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_13 * lsf_69[k]
                  + f_3 * pc_x[k] * msf_69[k];

        t_100[k] = f_14 * lsf_36[k]
                   + f_1 * msd0_39[k]
                   - f_2 * msd1_39[k]
                   + f_3 * pc_y[k] * msf_66[k];

        t_101[k] = f_3 * pc_z[k] * msf_66[k];

        t_102[k] = f_4 * msd0_39[k]
                   - f_5 * msd1_39[k]
                   + f_3 * pc_z[k] * msf_67[k];

        t_103[k] = f_14 * lsf_39[k]
                   + f_3 * pc_y[k] * msf_69[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pa_z, pc_y, pc_z, lsg0_45, lsf_30, \
                         lsf_40, lsg1_45, msd0_41, msd1_41, msf_69, \
                         msf_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_1 * msd0_41[k]
                   - f_2 * msd1_41[k]
                   + f_3 * pc_z[k] * msf_69[k];

        t_105[k] = pa_z[k] * lsg0_45[k]
                   - f_6 * pc_z[k] * lsg1_45[k];

        t_106[k] = f_8 * lsf_40[k]
                   + f_3 * pc_y[k] * msf_70[k];

        t_107[k] = f_7 * lsf_30[k]
                   + f_3 * pc_z[k] * msf_70[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, pa_z, pc_x, pc_y, pc_z, lsg0_48, lsf_42, lsf_75, \
                         lsg1_48, msd0_47, msd1_47, msf_72, msf_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = pa_z[k] * lsg0_48[k]
                   - f_6 * pc_z[k] * lsg1_48[k];

        t_109[k] = f_8 * lsf_42[k]
                   + f_3 * pc_y[k] * msf_72[k];

        t_110[k] = f_13 * lsf_75[k]
                   + f_4 * msd0_47[k]
                   - f_5 * msd1_47[k]
                   + f_3 * pc_x[k] * msf_75[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, pc_x, lsf_76, lsf_77, lsf_78, lsf_79, \
                         msf_76, msf_77, msf_78, msf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_13 * lsf_76[k]
                   + f_3 * pc_x[k] * msf_76[k];

        t_112[k] = f_13 * lsf_77[k]
                   + f_3 * pc_x[k] * msf_77[k];

        t_113[k] = f_13 * lsf_78[k]
                   + f_3 * pc_x[k] * msf_78[k];

        t_114[k] = f_13 * lsf_79[k]
                   + f_3 * pc_x[k] * msf_79[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, pa_z, pc_y, pc_z, lsg0_55, lsf_36, lsf_48, \
                         lsg1_55, msd0_47, msd1_47, msf_76, msf_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = pa_z[k] * lsg0_55[k]
                   - f_6 * pc_z[k] * lsg1_55[k];

        t_116[k] = f_7 * lsf_36[k]
                   + f_3 * pc_z[k] * msf_76[k];

        t_117[k] = f_8 * lsf_48[k]
                   + f_4 * msd0_47[k]
                   - f_5 * msd1_47[k]
                   + f_3 * pc_y[k] * msf_78[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, pa_y, pc_y, pc_z, lsg0_75, lsf_39, \
                         lsf_49, lsf_50, lsg1_75, msd0_47, msd1_47, msf_79, \
                         msf_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_8 * lsf_49[k]
                   + f_3 * pc_y[k] * msf_79[k];

        t_119[k] = f_7 * lsf_39[k]
                   + f_1 * msd0_47[k]
                   - f_2 * msd1_47[k]
                   + f_3 * pc_z[k] * msf_79[k];

        t_120[k] = pa_y[k] * lsg0_75[k]
                   - f_6 * pc_y[k] * lsg1_75[k];

        t_121[k] = f_7 * lsf_50[k]
                   + f_3 * pc_y[k] * msf_80[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, pa_y, pc_y, pc_z, lsg0_78, lsg0_80, \
                         lsf_40, lsf_51, lsf_52, lsg1_78, lsg1_80, msf_80, \
                         msf_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_8 * lsf_40[k]
                   + f_3 * pc_z[k] * msf_80[k];

        t_123[k] = pa_y[k] * lsg0_78[k]
                   + f_8 * lsf_51[k]
                   - f_6 * pc_y[k] * lsg1_78[k];

        t_124[k] = f_7 * lsf_52[k]
                   + f_3 * pc_y[k] * msf_82[k];

        t_125[k] = pa_y[k] * lsg0_80[k]
                   - f_6 * pc_y[k] * lsg1_80[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, pc_x, lsf_86, lsf_87, lsf_88, lsf_89, \
                         msf_86, msf_87, msf_88, msf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_13 * lsf_86[k]
                   + f_3 * pc_x[k] * msf_86[k];

        t_127[k] = f_13 * lsf_87[k]
                   + f_3 * pc_x[k] * msf_87[k];

        t_128[k] = f_13 * lsf_88[k]
                   + f_3 * pc_x[k] * msf_88[k];

        t_129[k] = f_13 * lsf_89[k]
                   + f_3 * pc_x[k] * msf_89[k];
    }
}

static auto
compute_prim_msg_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsg0,
                                                          const size_t lsf, const size_t lsg1,
                                                          const size_t msd0, const size_t msd1,
                                                          const size_t msf, const size_t ncols,
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
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_13 = 3.0 / q;
    const auto f_14 = 1.5 / q;
    const auto f_15 = 2.5 / q;
    const auto f_16 = 2.0 / q;

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
    auto *t_246 = buffer.data(target + 246);
    auto *t_247 = buffer.data(target + 247);
    auto *t_248 = buffer.data(target + 248);
    auto *t_249 = buffer.data(target + 249);
    auto *t_250 = buffer.data(target + 250);
    auto *t_251 = buffer.data(target + 251);
    auto *t_252 = buffer.data(target + 252);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *lsg0_89 = buffer.data(lsg0 + 89);
    const auto *lsg0_90 = buffer.data(lsg0 + 90);
    const auto *lsg0_93 = buffer.data(lsg0 + 93);
    const auto *lsg0_100 = buffer.data(lsg0 + 100);
    const auto *lsg0_135 = buffer.data(lsg0 + 135);
    const auto *lsg0_138 = buffer.data(lsg0 + 138);
    const auto *lsg0_140 = buffer.data(lsg0 + 140);
    const auto *lsg0_149 = buffer.data(lsg0 + 149);
    const auto *lsg0_150 = buffer.data(lsg0 + 150);
    const auto *lsg0_153 = buffer.data(lsg0 + 153);
    const auto *lsg0_160 = buffer.data(lsg0 + 160);

    const auto *lsf_46 = buffer.data(lsf + 46);
    const auto *lsf_50 = buffer.data(lsf + 50);
    const auto *lsf_56 = buffer.data(lsf + 56);
    const auto *lsf_58 = buffer.data(lsf + 58);
    const auto *lsf_59 = buffer.data(lsf + 59);
    const auto *lsf_60 = buffer.data(lsf + 60);
    const auto *lsf_66 = buffer.data(lsf + 66);
    const auto *lsf_69 = buffer.data(lsf + 69);
    const auto *lsf_70 = buffer.data(lsf + 70);
    const auto *lsf_72 = buffer.data(lsf + 72);
    const auto *lsf_76 = buffer.data(lsf + 76);
    const auto *lsf_78 = buffer.data(lsf + 78);
    const auto *lsf_79 = buffer.data(lsf + 79);
    const auto *lsf_80 = buffer.data(lsf + 80);
    const auto *lsf_82 = buffer.data(lsf + 82);
    const auto *lsf_86 = buffer.data(lsf + 86);
    const auto *lsf_88 = buffer.data(lsf + 88);
    const auto *lsf_89 = buffer.data(lsf + 89);
    const auto *lsf_90 = buffer.data(lsf + 90);
    const auto *lsf_91 = buffer.data(lsf + 91);
    const auto *lsf_92 = buffer.data(lsf + 92);
    const auto *lsf_95 = buffer.data(lsf + 95);
    const auto *lsf_96 = buffer.data(lsf + 96);
    const auto *lsf_97 = buffer.data(lsf + 97);
    const auto *lsf_98 = buffer.data(lsf + 98);
    const auto *lsf_99 = buffer.data(lsf + 99);
    const auto *lsf_100 = buffer.data(lsf + 100);
    const auto *lsf_103 = buffer.data(lsf + 103);
    const auto *lsf_106 = buffer.data(lsf + 106);
    const auto *lsf_108 = buffer.data(lsf + 108);
    const auto *lsf_109 = buffer.data(lsf + 109);
    const auto *lsf_110 = buffer.data(lsf + 110);
    const auto *lsf_112 = buffer.data(lsf + 112);
    const auto *lsf_115 = buffer.data(lsf + 115);
    const auto *lsf_116 = buffer.data(lsf + 116);
    const auto *lsf_117 = buffer.data(lsf + 117);
    const auto *lsf_118 = buffer.data(lsf + 118);
    const auto *lsf_119 = buffer.data(lsf + 119);
    const auto *lsf_120 = buffer.data(lsf + 120);
    const auto *lsf_123 = buffer.data(lsf + 123);
    const auto *lsf_125 = buffer.data(lsf + 125);
    const auto *lsf_126 = buffer.data(lsf + 126);
    const auto *lsf_127 = buffer.data(lsf + 127);
    const auto *lsf_128 = buffer.data(lsf + 128);
    const auto *lsf_129 = buffer.data(lsf + 129);
    const auto *lsf_136 = buffer.data(lsf + 136);
    const auto *lsf_137 = buffer.data(lsf + 137);
    const auto *lsf_138 = buffer.data(lsf + 138);
    const auto *lsf_139 = buffer.data(lsf + 139);
    const auto *lsf_140 = buffer.data(lsf + 140);
    const auto *lsf_145 = buffer.data(lsf + 145);
    const auto *lsf_146 = buffer.data(lsf + 146);
    const auto *lsf_147 = buffer.data(lsf + 147);
    const auto *lsf_149 = buffer.data(lsf + 149);
    const auto *lsf_150 = buffer.data(lsf + 150);
    const auto *lsf_153 = buffer.data(lsf + 153);
    const auto *lsf_156 = buffer.data(lsf + 156);
    const auto *lsf_158 = buffer.data(lsf + 158);
    const auto *lsf_159 = buffer.data(lsf + 159);
    const auto *lsf_165 = buffer.data(lsf + 165);
    const auto *lsf_166 = buffer.data(lsf + 166);
    const auto *lsf_167 = buffer.data(lsf + 167);
    const auto *lsf_168 = buffer.data(lsf + 168);
    const auto *lsf_169 = buffer.data(lsf + 169);

    const auto *lsg1_89 = buffer.data(lsg1 + 89);
    const auto *lsg1_90 = buffer.data(lsg1 + 90);
    const auto *lsg1_93 = buffer.data(lsg1 + 93);
    const auto *lsg1_100 = buffer.data(lsg1 + 100);
    const auto *lsg1_135 = buffer.data(lsg1 + 135);
    const auto *lsg1_138 = buffer.data(lsg1 + 138);
    const auto *lsg1_140 = buffer.data(lsg1 + 140);
    const auto *lsg1_149 = buffer.data(lsg1 + 149);
    const auto *lsg1_150 = buffer.data(lsg1 + 150);
    const auto *lsg1_153 = buffer.data(lsg1 + 153);
    const auto *lsg1_160 = buffer.data(lsg1 + 160);

    const auto *msd0_51 = buffer.data(msd0 + 51);
    const auto *msd0_53 = buffer.data(msd0 + 53);
    const auto *msd0_54 = buffer.data(msd0 + 54);
    const auto *msd0_57 = buffer.data(msd0 + 57);
    const auto *msd0_58 = buffer.data(msd0 + 58);
    const auto *msd0_59 = buffer.data(msd0 + 59);
    const auto *msd0_60 = buffer.data(msd0 + 60);
    const auto *msd0_63 = buffer.data(msd0 + 63);
    const auto *msd0_65 = buffer.data(msd0 + 65);
    const auto *msd0_71 = buffer.data(msd0 + 71);
    const auto *msd0_72 = buffer.data(msd0 + 72);
    const auto *msd0_75 = buffer.data(msd0 + 75);
    const auto *msd0_77 = buffer.data(msd0 + 77);
    const auto *msd0_81 = buffer.data(msd0 + 81);
    const auto *msd0_83 = buffer.data(msd0 + 83);
    const auto *msd0_84 = buffer.data(msd0 + 84);
    const auto *msd0_87 = buffer.data(msd0 + 87);
    const auto *msd0_88 = buffer.data(msd0 + 88);
    const auto *msd0_89 = buffer.data(msd0 + 89);
    const auto *msd0_90 = buffer.data(msd0 + 90);
    const auto *msd0_93 = buffer.data(msd0 + 93);
    const auto *msd0_95 = buffer.data(msd0 + 95);
    const auto *msd0_101 = buffer.data(msd0 + 101);

    const auto *msd1_51 = buffer.data(msd1 + 51);
    const auto *msd1_53 = buffer.data(msd1 + 53);
    const auto *msd1_54 = buffer.data(msd1 + 54);
    const auto *msd1_57 = buffer.data(msd1 + 57);
    const auto *msd1_58 = buffer.data(msd1 + 58);
    const auto *msd1_59 = buffer.data(msd1 + 59);
    const auto *msd1_60 = buffer.data(msd1 + 60);
    const auto *msd1_63 = buffer.data(msd1 + 63);
    const auto *msd1_65 = buffer.data(msd1 + 65);
    const auto *msd1_71 = buffer.data(msd1 + 71);
    const auto *msd1_72 = buffer.data(msd1 + 72);
    const auto *msd1_75 = buffer.data(msd1 + 75);
    const auto *msd1_77 = buffer.data(msd1 + 77);
    const auto *msd1_81 = buffer.data(msd1 + 81);
    const auto *msd1_83 = buffer.data(msd1 + 83);
    const auto *msd1_84 = buffer.data(msd1 + 84);
    const auto *msd1_87 = buffer.data(msd1 + 87);
    const auto *msd1_88 = buffer.data(msd1 + 88);
    const auto *msd1_89 = buffer.data(msd1 + 89);
    const auto *msd1_90 = buffer.data(msd1 + 90);
    const auto *msd1_93 = buffer.data(msd1 + 93);
    const auto *msd1_95 = buffer.data(msd1 + 95);
    const auto *msd1_101 = buffer.data(msd1 + 101);

    const auto *msf_86 = buffer.data(msf + 86);
    const auto *msf_88 = buffer.data(msf + 88);
    const auto *msf_89 = buffer.data(msf + 89);
    const auto *msf_90 = buffer.data(msf + 90);
    const auto *msf_91 = buffer.data(msf + 91);
    const auto *msf_92 = buffer.data(msf + 92);
    const auto *msf_95 = buffer.data(msf + 95);
    const auto *msf_96 = buffer.data(msf + 96);
    const auto *msf_97 = buffer.data(msf + 97);
    const auto *msf_98 = buffer.data(msf + 98);
    const auto *msf_99 = buffer.data(msf + 99);
    const auto *msf_100 = buffer.data(msf + 100);
    const auto *msf_101 = buffer.data(msf + 101);
    const auto *msf_102 = buffer.data(msf + 102);
    const auto *msf_103 = buffer.data(msf + 103);
    const auto *msf_106 = buffer.data(msf + 106);
    const auto *msf_107 = buffer.data(msf + 107);
    const auto *msf_108 = buffer.data(msf + 108);
    const auto *msf_109 = buffer.data(msf + 109);
    const auto *msf_110 = buffer.data(msf + 110);
    const auto *msf_112 = buffer.data(msf + 112);
    const auto *msf_115 = buffer.data(msf + 115);
    const auto *msf_116 = buffer.data(msf + 116);
    const auto *msf_117 = buffer.data(msf + 117);
    const auto *msf_118 = buffer.data(msf + 118);
    const auto *msf_119 = buffer.data(msf + 119);
    const auto *msf_120 = buffer.data(msf + 120);
    const auto *msf_122 = buffer.data(msf + 122);
    const auto *msf_123 = buffer.data(msf + 123);
    const auto *msf_125 = buffer.data(msf + 125);
    const auto *msf_126 = buffer.data(msf + 126);
    const auto *msf_127 = buffer.data(msf + 127);
    const auto *msf_128 = buffer.data(msf + 128);
    const auto *msf_129 = buffer.data(msf + 129);
    const auto *msf_130 = buffer.data(msf + 130);
    const auto *msf_132 = buffer.data(msf + 132);
    const auto *msf_136 = buffer.data(msf + 136);
    const auto *msf_137 = buffer.data(msf + 137);
    const auto *msf_138 = buffer.data(msf + 138);
    const auto *msf_139 = buffer.data(msf + 139);
    const auto *msf_140 = buffer.data(msf + 140);
    const auto *msf_141 = buffer.data(msf + 141);
    const auto *msf_142 = buffer.data(msf + 142);
    const auto *msf_145 = buffer.data(msf + 145);
    const auto *msf_146 = buffer.data(msf + 146);
    const auto *msf_147 = buffer.data(msf + 147);
    const auto *msf_148 = buffer.data(msf + 148);
    const auto *msf_149 = buffer.data(msf + 149);
    const auto *msf_150 = buffer.data(msf + 150);
    const auto *msf_151 = buffer.data(msf + 151);
    const auto *msf_152 = buffer.data(msf + 152);
    const auto *msf_153 = buffer.data(msf + 153);
    const auto *msf_156 = buffer.data(msf + 156);
    const auto *msf_157 = buffer.data(msf + 157);
    const auto *msf_158 = buffer.data(msf + 158);
    const auto *msf_159 = buffer.data(msf + 159);
    const auto *msf_160 = buffer.data(msf + 160);
    const auto *msf_162 = buffer.data(msf + 162);
    const auto *msf_165 = buffer.data(msf + 165);
    const auto *msf_166 = buffer.data(msf + 166);
    const auto *msf_167 = buffer.data(msf + 167);
    const auto *msf_168 = buffer.data(msf + 168);
    const auto *msf_169 = buffer.data(msf + 169);

#pragma omp simd aligned(t_130, t_131, t_132, pc_y, pc_z, lsf_46, lsf_56, lsf_58, msd0_51, \
                         msd0_53, msd1_51, msd1_53, msf_86, msf_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_7 * lsf_56[k]
                   + f_1 * msd0_51[k]
                   - f_2 * msd1_51[k]
                   + f_3 * pc_y[k] * msf_86[k];

        t_131[k] = f_8 * lsf_46[k]
                   + f_3 * pc_z[k] * msf_86[k];

        t_132[k] = f_7 * lsf_58[k]
                   + f_4 * msd0_53[k]
                   - f_5 * msd1_53[k]
                   + f_3 * pc_y[k] * msf_88[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, pa_y, pc_x, pc_y, lsg0_89, lsf_59, \
                         lsf_90, lsg1_89, msd0_54, msd1_54, msf_89, \
                         msf_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_7 * lsf_59[k]
                   + f_3 * pc_y[k] * msf_89[k];

        t_134[k] = pa_y[k] * lsg0_89[k]
                   - f_6 * pc_y[k] * lsg1_89[k];

        t_135[k] = f_13 * lsf_90[k]
                   + f_1 * msd0_54[k]
                   - f_2 * msd1_54[k]
                   + f_3 * pc_x[k] * msf_90[k];

        t_136[k] = f_3 * pc_y[k] * msf_90[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, pc_y, pc_z, lsf_50, msd0_54, msd1_54, msf_90, \
                         msf_91, msf_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_14 * lsf_50[k]
                   + f_3 * pc_z[k] * msf_90[k];

        t_138[k] = f_4 * msd0_54[k]
                   - f_5 * msd1_54[k]
                   + f_3 * pc_y[k] * msf_91[k];

        t_139[k] = f_3 * pc_y[k] * msf_92[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, pc_x, pc_y, lsf_95, lsf_96, lsf_97, \
                         msd0_59, msd1_59, msf_95, msf_96, msf_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_13 * lsf_95[k]
                   + f_4 * msd0_59[k]
                   - f_5 * msd1_59[k]
                   + f_3 * pc_x[k] * msf_95[k];

        t_141[k] = f_13 * lsf_96[k]
                   + f_3 * pc_x[k] * msf_96[k];

        t_142[k] = f_13 * lsf_97[k]
                   + f_3 * pc_x[k] * msf_97[k];

        t_143[k] = f_3 * pc_y[k] * msf_95[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, pc_x, pc_y, lsf_99, msd0_57, msd0_58, msd1_57, \
                         msd1_58, msf_96, msf_97, msf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_13 * lsf_99[k]
                   + f_3 * pc_x[k] * msf_99[k];

        t_145[k] = f_1 * msd0_57[k]
                   - f_2 * msd1_57[k]
                   + f_3 * pc_y[k] * msf_96[k];

        t_146[k] = f_10 * msd0_58[k]
                   - f_11 * msd1_58[k]
                   + f_3 * pc_y[k] * msf_97[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, pc_x, pc_y, pc_z, lsf_59, lsf_100, \
                         msd0_59, msd0_60, msd1_59, msd1_60, msf_98, msf_99, \
                         msf_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = f_4 * msd0_59[k]
                   - f_5 * msd1_59[k]
                   + f_3 * pc_y[k] * msf_98[k];

        t_148[k] = f_3 * pc_y[k] * msf_99[k];

        t_149[k] = f_14 * lsf_59[k]
                   + f_1 * msd0_59[k]
                   - f_2 * msd1_59[k]
                   + f_3 * pc_z[k] * msf_99[k];

        t_150[k] = f_15 * lsf_100[k]
                   + f_1 * msd0_60[k]
                   - f_2 * msd1_60[k]
                   + f_3 * pc_x[k] * msf_100[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, pc_x, pc_y, pc_z, lsf_60, lsf_103, \
                         msd0_63, msd1_63, msf_100, msf_101, msf_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_16 * lsf_60[k]
                   + f_3 * pc_y[k] * msf_100[k];

        t_152[k] = f_3 * pc_z[k] * msf_100[k];

        t_153[k] = f_15 * lsf_103[k]
                   + f_4 * msd0_63[k]
                   - f_5 * msd1_63[k]
                   + f_3 * pc_x[k] * msf_103[k];

        t_154[k] = f_3 * pc_z[k] * msf_101[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, pc_x, pc_z, lsf_106, lsf_108, msd0_60, \
                         msd1_60, msf_102, msf_103, msf_106, msf_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_4 * msd0_60[k]
                   - f_5 * msd1_60[k]
                   + f_3 * pc_z[k] * msf_102[k];

        t_156[k] = f_15 * lsf_106[k]
                   + f_3 * pc_x[k] * msf_106[k];

        t_157[k] = f_3 * pc_z[k] * msf_103[k];

        t_158[k] = f_15 * lsf_108[k]
                   + f_3 * pc_x[k] * msf_108[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, t_163, pc_x, pc_y, pc_z, lsf_66, lsf_69, \
                         lsf_109, msd0_63, msd1_63, msf_106, msf_107, \
                         msf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_15 * lsf_109[k]
                   + f_3 * pc_x[k] * msf_109[k];

        t_160[k] = f_16 * lsf_66[k]
                   + f_1 * msd0_63[k]
                   - f_2 * msd1_63[k]
                   + f_3 * pc_y[k] * msf_106[k];

        t_161[k] = f_3 * pc_z[k] * msf_106[k];

        t_162[k] = f_4 * msd0_63[k]
                   - f_5 * msd1_63[k]
                   + f_3 * pc_z[k] * msf_107[k];

        t_163[k] = f_16 * lsf_69[k]
                   + f_3 * pc_y[k] * msf_109[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, pa_z, pc_y, pc_z, lsg0_90, lsf_60, \
                         lsf_70, lsg1_90, msd0_65, msd1_65, msf_109, \
                         msf_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_1 * msd0_65[k]
                   - f_2 * msd1_65[k]
                   + f_3 * pc_z[k] * msf_109[k];

        t_165[k] = pa_z[k] * lsg0_90[k]
                   - f_6 * pc_z[k] * lsg1_90[k];

        t_166[k] = f_14 * lsf_70[k]
                   + f_3 * pc_y[k] * msf_110[k];

        t_167[k] = f_7 * lsf_60[k]
                   + f_3 * pc_z[k] * msf_110[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, pa_z, pc_x, pc_y, pc_z, lsg0_93, lsf_72, \
                         lsf_115, lsg1_93, msd0_71, msd1_71, msf_112, \
                         msf_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = pa_z[k] * lsg0_93[k]
                   - f_6 * pc_z[k] * lsg1_93[k];

        t_169[k] = f_14 * lsf_72[k]
                   + f_3 * pc_y[k] * msf_112[k];

        t_170[k] = f_15 * lsf_115[k]
                   + f_4 * msd0_71[k]
                   - f_5 * msd1_71[k]
                   + f_3 * pc_x[k] * msf_115[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, t_174, pc_x, lsf_116, lsf_117, lsf_118, lsf_119, \
                         msf_116, msf_117, msf_118, msf_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = f_15 * lsf_116[k]
                   + f_3 * pc_x[k] * msf_116[k];

        t_172[k] = f_15 * lsf_117[k]
                   + f_3 * pc_x[k] * msf_117[k];

        t_173[k] = f_15 * lsf_118[k]
                   + f_3 * pc_x[k] * msf_118[k];

        t_174[k] = f_15 * lsf_119[k]
                   + f_3 * pc_x[k] * msf_119[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, pa_z, pc_y, pc_z, lsg0_100, lsf_66, lsf_78, \
                         lsg1_100, msd0_71, msd1_71, msf_116, msf_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = pa_z[k] * lsg0_100[k]
                   - f_6 * pc_z[k] * lsg1_100[k];

        t_176[k] = f_7 * lsf_66[k]
                   + f_3 * pc_z[k] * msf_116[k];

        t_177[k] = f_14 * lsf_78[k]
                   + f_4 * msd0_71[k]
                   - f_5 * msd1_71[k]
                   + f_3 * pc_y[k] * msf_118[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, pc_x, pc_y, pc_z, lsf_69, lsf_79, lsf_120, \
                         msd0_71, msd0_72, msd1_71, msd1_72, msf_119, \
                         msf_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_14 * lsf_79[k]
                   + f_3 * pc_y[k] * msf_119[k];

        t_179[k] = f_7 * lsf_69[k]
                   + f_1 * msd0_71[k]
                   - f_2 * msd1_71[k]
                   + f_3 * pc_z[k] * msf_119[k];

        t_180[k] = f_15 * lsf_120[k]
                   + f_1 * msd0_72[k]
                   - f_2 * msd1_72[k]
                   + f_3 * pc_x[k] * msf_120[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pc_x, pc_y, pc_z, lsf_70, lsf_80, lsf_82, \
                         lsf_123, msd0_75, msd1_75, msf_120, msf_122, \
                         msf_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_8 * lsf_80[k]
                   + f_3 * pc_y[k] * msf_120[k];

        t_182[k] = f_8 * lsf_70[k]
                   + f_3 * pc_z[k] * msf_120[k];

        t_183[k] = f_15 * lsf_123[k]
                   + f_4 * msd0_75[k]
                   - f_5 * msd1_75[k]
                   + f_3 * pc_x[k] * msf_123[k];

        t_184[k] = f_8 * lsf_82[k]
                   + f_3 * pc_y[k] * msf_122[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pc_x, lsf_125, lsf_126, lsf_127, lsf_128, \
                         msd0_77, msd1_77, msf_125, msf_126, msf_127, \
                         msf_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_15 * lsf_125[k]
                   + f_4 * msd0_77[k]
                   - f_5 * msd1_77[k]
                   + f_3 * pc_x[k] * msf_125[k];

        t_186[k] = f_15 * lsf_126[k]
                   + f_3 * pc_x[k] * msf_126[k];

        t_187[k] = f_15 * lsf_127[k]
                   + f_3 * pc_x[k] * msf_127[k];

        t_188[k] = f_15 * lsf_128[k]
                   + f_3 * pc_x[k] * msf_128[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, pc_x, pc_y, pc_z, lsf_76, lsf_86, lsf_129, \
                         msd0_75, msd1_75, msf_126, msf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_15 * lsf_129[k]
                   + f_3 * pc_x[k] * msf_129[k];

        t_190[k] = f_8 * lsf_86[k]
                   + f_1 * msd0_75[k]
                   - f_2 * msd1_75[k]
                   + f_3 * pc_y[k] * msf_126[k];

        t_191[k] = f_8 * lsf_76[k]
                   + f_3 * pc_z[k] * msf_126[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, pa_y, pc_y, pc_z, lsg0_135, lsf_79, \
                         lsf_88, lsf_89, lsg1_135, msd0_77, msd1_77, msf_128, \
                         msf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_8 * lsf_88[k]
                   + f_4 * msd0_77[k]
                   - f_5 * msd1_77[k]
                   + f_3 * pc_y[k] * msf_128[k];

        t_193[k] = f_8 * lsf_89[k]
                   + f_3 * pc_y[k] * msf_129[k];

        t_194[k] = f_8 * lsf_79[k]
                   + f_1 * msd0_77[k]
                   - f_2 * msd1_77[k]
                   + f_3 * pc_z[k] * msf_129[k];

        t_195[k] = pa_y[k] * lsg0_135[k]
                   - f_6 * pc_y[k] * lsg1_135[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, pa_y, pc_y, pc_z, lsg0_138, lsf_80, \
                         lsf_90, lsf_91, lsf_92, lsg1_138, msf_130, \
                         msf_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = f_7 * lsf_90[k]
                   + f_3 * pc_y[k] * msf_130[k];

        t_197[k] = f_14 * lsf_80[k]
                   + f_3 * pc_z[k] * msf_130[k];

        t_198[k] = pa_y[k] * lsg0_138[k]
                   + f_8 * lsf_91[k]
                   - f_6 * pc_y[k] * lsg1_138[k];

        t_199[k] = f_7 * lsf_92[k]
                   + f_3 * pc_y[k] * msf_132[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, pa_y, pc_x, pc_y, lsg0_140, lsf_136, \
                         lsf_137, lsf_138, lsg1_140, msf_136, msf_137, \
                         msf_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = pa_y[k] * lsg0_140[k]
                   - f_6 * pc_y[k] * lsg1_140[k];

        t_201[k] = f_15 * lsf_136[k]
                   + f_3 * pc_x[k] * msf_136[k];

        t_202[k] = f_15 * lsf_137[k]
                   + f_3 * pc_x[k] * msf_137[k];

        t_203[k] = f_15 * lsf_138[k]
                   + f_3 * pc_x[k] * msf_138[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, pc_x, pc_y, pc_z, lsf_86, lsf_96, lsf_139, \
                         msd0_81, msd1_81, msf_136, msf_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_15 * lsf_139[k]
                   + f_3 * pc_x[k] * msf_139[k];

        t_205[k] = f_7 * lsf_96[k]
                   + f_1 * msd0_81[k]
                   - f_2 * msd1_81[k]
                   + f_3 * pc_y[k] * msf_136[k];

        t_206[k] = f_14 * lsf_86[k]
                   + f_3 * pc_z[k] * msf_136[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, pa_y, pc_y, lsg0_149, lsf_98, lsf_99, lsg1_149, \
                         msd0_83, msd1_83, msf_138, msf_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = f_7 * lsf_98[k]
                   + f_4 * msd0_83[k]
                   - f_5 * msd1_83[k]
                   + f_3 * pc_y[k] * msf_138[k];

        t_208[k] = f_7 * lsf_99[k]
                   + f_3 * pc_y[k] * msf_139[k];

        t_209[k] = pa_y[k] * lsg0_149[k]
                   - f_6 * pc_y[k] * lsg1_149[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, pc_x, pc_y, pc_z, lsf_90, lsf_140, \
                         msd0_84, msd1_84, msf_140, msf_141, msf_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_15 * lsf_140[k]
                   + f_1 * msd0_84[k]
                   - f_2 * msd1_84[k]
                   + f_3 * pc_x[k] * msf_140[k];

        t_211[k] = f_3 * pc_y[k] * msf_140[k];

        t_212[k] = f_16 * lsf_90[k]
                   + f_3 * pc_z[k] * msf_140[k];

        t_213[k] = f_4 * msd0_84[k]
                   - f_5 * msd1_84[k]
                   + f_3 * pc_y[k] * msf_141[k];

        t_214[k] = f_3 * pc_y[k] * msf_142[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, pc_x, pc_y, lsf_145, lsf_146, lsf_147, \
                         msd0_89, msd1_89, msf_145, msf_146, msf_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = f_15 * lsf_145[k]
                   + f_4 * msd0_89[k]
                   - f_5 * msd1_89[k]
                   + f_3 * pc_x[k] * msf_145[k];

        t_216[k] = f_15 * lsf_146[k]
                   + f_3 * pc_x[k] * msf_146[k];

        t_217[k] = f_15 * lsf_147[k]
                   + f_3 * pc_x[k] * msf_147[k];

        t_218[k] = f_3 * pc_y[k] * msf_145[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, pc_x, pc_y, lsf_149, msd0_87, msd0_88, msd1_87, \
                         msd1_88, msf_146, msf_147, msf_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = f_15 * lsf_149[k]
                   + f_3 * pc_x[k] * msf_149[k];

        t_220[k] = f_1 * msd0_87[k]
                   - f_2 * msd1_87[k]
                   + f_3 * pc_y[k] * msf_146[k];

        t_221[k] = f_10 * msd0_88[k]
                   - f_11 * msd1_88[k]
                   + f_3 * pc_y[k] * msf_147[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, pc_x, pc_y, pc_z, lsf_99, lsf_150, \
                         msd0_89, msd0_90, msd1_89, msd1_90, msf_148, msf_149, \
                         msf_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_4 * msd0_89[k]
                   - f_5 * msd1_89[k]
                   + f_3 * pc_y[k] * msf_148[k];

        t_223[k] = f_3 * pc_y[k] * msf_149[k];

        t_224[k] = f_16 * lsf_99[k]
                   + f_1 * msd0_89[k]
                   - f_2 * msd1_89[k]
                   + f_3 * pc_z[k] * msf_149[k];

        t_225[k] = f_16 * lsf_150[k]
                   + f_1 * msd0_90[k]
                   - f_2 * msd1_90[k]
                   + f_3 * pc_x[k] * msf_150[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, t_229, pc_x, pc_y, pc_z, lsf_100, lsf_153, \
                         msd0_93, msd1_93, msf_150, msf_151, msf_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = f_15 * lsf_100[k]
                   + f_3 * pc_y[k] * msf_150[k];

        t_227[k] = f_3 * pc_z[k] * msf_150[k];

        t_228[k] = f_16 * lsf_153[k]
                   + f_4 * msd0_93[k]
                   - f_5 * msd1_93[k]
                   + f_3 * pc_x[k] * msf_153[k];

        t_229[k] = f_3 * pc_z[k] * msf_151[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, pc_x, pc_z, lsf_156, lsf_158, msd0_90, \
                         msd1_90, msf_152, msf_153, msf_156, msf_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_4 * msd0_90[k]
                   - f_5 * msd1_90[k]
                   + f_3 * pc_z[k] * msf_152[k];

        t_231[k] = f_16 * lsf_156[k]
                   + f_3 * pc_x[k] * msf_156[k];

        t_232[k] = f_3 * pc_z[k] * msf_153[k];

        t_233[k] = f_16 * lsf_158[k]
                   + f_3 * pc_x[k] * msf_158[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, t_238, pc_x, pc_y, pc_z, lsf_106, \
                         lsf_109, lsf_159, msd0_93, msd1_93, msf_156, msf_157, \
                         msf_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_16 * lsf_159[k]
                   + f_3 * pc_x[k] * msf_159[k];

        t_235[k] = f_15 * lsf_106[k]
                   + f_1 * msd0_93[k]
                   - f_2 * msd1_93[k]
                   + f_3 * pc_y[k] * msf_156[k];

        t_236[k] = f_3 * pc_z[k] * msf_156[k];

        t_237[k] = f_4 * msd0_93[k]
                   - f_5 * msd1_93[k]
                   + f_3 * pc_z[k] * msf_157[k];

        t_238[k] = f_15 * lsf_109[k]
                   + f_3 * pc_y[k] * msf_159[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, pa_z, pc_y, pc_z, lsg0_150, lsf_100, \
                         lsf_110, lsg1_150, msd0_95, msd1_95, msf_159, \
                         msf_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_1 * msd0_95[k]
                   - f_2 * msd1_95[k]
                   + f_3 * pc_z[k] * msf_159[k];

        t_240[k] = pa_z[k] * lsg0_150[k]
                   - f_6 * pc_z[k] * lsg1_150[k];

        t_241[k] = f_16 * lsf_110[k]
                   + f_3 * pc_y[k] * msf_160[k];

        t_242[k] = f_7 * lsf_100[k]
                   + f_3 * pc_z[k] * msf_160[k];
    }

#pragma omp simd aligned(t_243, t_244, t_245, pa_z, pc_x, pc_y, pc_z, lsg0_153, lsf_112, \
                         lsf_165, lsg1_153, msd0_101, msd1_101, msf_162, \
                         msf_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = pa_z[k] * lsg0_153[k]
                   - f_6 * pc_z[k] * lsg1_153[k];

        t_244[k] = f_16 * lsf_112[k]
                   + f_3 * pc_y[k] * msf_162[k];

        t_245[k] = f_16 * lsf_165[k]
                   + f_4 * msd0_101[k]
                   - f_5 * msd1_101[k]
                   + f_3 * pc_x[k] * msf_165[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, t_249, pc_x, lsf_166, lsf_167, lsf_168, lsf_169, \
                         msf_166, msf_167, msf_168, msf_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = f_16 * lsf_166[k]
                   + f_3 * pc_x[k] * msf_166[k];

        t_247[k] = f_16 * lsf_167[k]
                   + f_3 * pc_x[k] * msf_167[k];

        t_248[k] = f_16 * lsf_168[k]
                   + f_3 * pc_x[k] * msf_168[k];

        t_249[k] = f_16 * lsf_169[k]
                   + f_3 * pc_x[k] * msf_169[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, pa_z, pc_y, pc_z, lsg0_160, lsf_106, lsf_118, \
                         lsg1_160, msd0_101, msd1_101, msf_166, \
                         msf_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = pa_z[k] * lsg0_160[k]
                   - f_6 * pc_z[k] * lsg1_160[k];

        t_251[k] = f_7 * lsf_106[k]
                   + f_3 * pc_z[k] * msf_166[k];

        t_252[k] = f_16 * lsf_118[k]
                   + f_4 * msd0_101[k]
                   - f_5 * msd1_101[k]
                   + f_3 * pc_y[k] * msf_168[k];
    }
}

static auto
compute_prim_msg_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsg0,
                                                          const size_t lsf, const size_t lsg1,
                                                          const size_t msd0, const size_t msd1,
                                                          const size_t msf, const size_t ncols,
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
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_13 = 3.0 / q;
    const auto f_14 = 1.5 / q;
    const auto f_15 = 2.5 / q;
    const auto f_16 = 2.0 / q;

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
    auto *t_361 = buffer.data(target + 361);
    auto *t_362 = buffer.data(target + 362);
    auto *t_363 = buffer.data(target + 363);
    auto *t_364 = buffer.data(target + 364);
    auto *t_365 = buffer.data(target + 365);
    auto *t_366 = buffer.data(target + 366);
    auto *t_367 = buffer.data(target + 367);
    auto *t_368 = buffer.data(target + 368);
    auto *t_369 = buffer.data(target + 369);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *lsg0_210 = buffer.data(lsg0 + 210);
    const auto *lsg0_213 = buffer.data(lsg0 + 213);
    const auto *lsg0_215 = buffer.data(lsg0 + 215);
    const auto *lsg0_224 = buffer.data(lsg0 + 224);
    const auto *lsg0_225 = buffer.data(lsg0 + 225);
    const auto *lsg0_228 = buffer.data(lsg0 + 228);
    const auto *lsg0_235 = buffer.data(lsg0 + 235);

    const auto *lsf_109 = buffer.data(lsf + 109);
    const auto *lsf_110 = buffer.data(lsf + 110);
    const auto *lsf_116 = buffer.data(lsf + 116);
    const auto *lsf_119 = buffer.data(lsf + 119);
    const auto *lsf_120 = buffer.data(lsf + 120);
    const auto *lsf_122 = buffer.data(lsf + 122);
    const auto *lsf_126 = buffer.data(lsf + 126);
    const auto *lsf_128 = buffer.data(lsf + 128);
    const auto *lsf_129 = buffer.data(lsf + 129);
    const auto *lsf_130 = buffer.data(lsf + 130);
    const auto *lsf_132 = buffer.data(lsf + 132);
    const auto *lsf_136 = buffer.data(lsf + 136);
    const auto *lsf_138 = buffer.data(lsf + 138);
    const auto *lsf_139 = buffer.data(lsf + 139);
    const auto *lsf_140 = buffer.data(lsf + 140);
    const auto *lsf_141 = buffer.data(lsf + 141);
    const auto *lsf_142 = buffer.data(lsf + 142);
    const auto *lsf_146 = buffer.data(lsf + 146);
    const auto *lsf_148 = buffer.data(lsf + 148);
    const auto *lsf_149 = buffer.data(lsf + 149);
    const auto *lsf_150 = buffer.data(lsf + 150);
    const auto *lsf_156 = buffer.data(lsf + 156);
    const auto *lsf_159 = buffer.data(lsf + 159);
    const auto *lsf_160 = buffer.data(lsf + 160);
    const auto *lsf_162 = buffer.data(lsf + 162);
    const auto *lsf_166 = buffer.data(lsf + 166);
    const auto *lsf_168 = buffer.data(lsf + 168);
    const auto *lsf_169 = buffer.data(lsf + 169);
    const auto *lsf_170 = buffer.data(lsf + 170);
    const auto *lsf_172 = buffer.data(lsf + 172);
    const auto *lsf_173 = buffer.data(lsf + 173);
    const auto *lsf_175 = buffer.data(lsf + 175);
    const auto *lsf_176 = buffer.data(lsf + 176);
    const auto *lsf_177 = buffer.data(lsf + 177);
    const auto *lsf_178 = buffer.data(lsf + 178);
    const auto *lsf_179 = buffer.data(lsf + 179);
    const auto *lsf_180 = buffer.data(lsf + 180);
    const auto *lsf_182 = buffer.data(lsf + 182);
    const auto *lsf_183 = buffer.data(lsf + 183);
    const auto *lsf_185 = buffer.data(lsf + 185);
    const auto *lsf_186 = buffer.data(lsf + 186);
    const auto *lsf_187 = buffer.data(lsf + 187);
    const auto *lsf_188 = buffer.data(lsf + 188);
    const auto *lsf_189 = buffer.data(lsf + 189);
    const auto *lsf_196 = buffer.data(lsf + 196);
    const auto *lsf_197 = buffer.data(lsf + 197);
    const auto *lsf_198 = buffer.data(lsf + 198);
    const auto *lsf_199 = buffer.data(lsf + 199);
    const auto *lsf_200 = buffer.data(lsf + 200);
    const auto *lsf_205 = buffer.data(lsf + 205);
    const auto *lsf_206 = buffer.data(lsf + 206);
    const auto *lsf_207 = buffer.data(lsf + 207);
    const auto *lsf_209 = buffer.data(lsf + 209);
    const auto *lsf_210 = buffer.data(lsf + 210);
    const auto *lsf_213 = buffer.data(lsf + 213);
    const auto *lsf_216 = buffer.data(lsf + 216);
    const auto *lsf_218 = buffer.data(lsf + 218);
    const auto *lsf_219 = buffer.data(lsf + 219);
    const auto *lsf_225 = buffer.data(lsf + 225);
    const auto *lsf_226 = buffer.data(lsf + 226);
    const auto *lsf_227 = buffer.data(lsf + 227);
    const auto *lsf_228 = buffer.data(lsf + 228);
    const auto *lsf_229 = buffer.data(lsf + 229);
    const auto *lsf_230 = buffer.data(lsf + 230);
    const auto *lsf_233 = buffer.data(lsf + 233);
    const auto *lsf_235 = buffer.data(lsf + 235);
    const auto *lsf_236 = buffer.data(lsf + 236);
    const auto *lsf_237 = buffer.data(lsf + 237);
    const auto *lsf_238 = buffer.data(lsf + 238);
    const auto *lsf_239 = buffer.data(lsf + 239);
    const auto *lsf_240 = buffer.data(lsf + 240);
    const auto *lsf_243 = buffer.data(lsf + 243);
    const auto *lsf_245 = buffer.data(lsf + 245);
    const auto *lsf_246 = buffer.data(lsf + 246);
    const auto *lsf_247 = buffer.data(lsf + 247);
    const auto *lsf_248 = buffer.data(lsf + 248);
    const auto *lsf_249 = buffer.data(lsf + 249);

    const auto *lsg1_210 = buffer.data(lsg1 + 210);
    const auto *lsg1_213 = buffer.data(lsg1 + 213);
    const auto *lsg1_215 = buffer.data(lsg1 + 215);
    const auto *lsg1_224 = buffer.data(lsg1 + 224);
    const auto *lsg1_225 = buffer.data(lsg1 + 225);
    const auto *lsg1_228 = buffer.data(lsg1 + 228);
    const auto *lsg1_235 = buffer.data(lsg1 + 235);

    const auto *msd0_101 = buffer.data(msd0 + 101);
    const auto *msd0_102 = buffer.data(msd0 + 102);
    const auto *msd0_105 = buffer.data(msd0 + 105);
    const auto *msd0_107 = buffer.data(msd0 + 107);
    const auto *msd0_108 = buffer.data(msd0 + 108);
    const auto *msd0_111 = buffer.data(msd0 + 111);
    const auto *msd0_113 = buffer.data(msd0 + 113);
    const auto *msd0_117 = buffer.data(msd0 + 117);
    const auto *msd0_119 = buffer.data(msd0 + 119);
    const auto *msd0_120 = buffer.data(msd0 + 120);
    const auto *msd0_123 = buffer.data(msd0 + 123);
    const auto *msd0_124 = buffer.data(msd0 + 124);
    const auto *msd0_125 = buffer.data(msd0 + 125);
    const auto *msd0_126 = buffer.data(msd0 + 126);
    const auto *msd0_129 = buffer.data(msd0 + 129);
    const auto *msd0_131 = buffer.data(msd0 + 131);
    const auto *msd0_137 = buffer.data(msd0 + 137);
    const auto *msd0_138 = buffer.data(msd0 + 138);
    const auto *msd0_141 = buffer.data(msd0 + 141);
    const auto *msd0_143 = buffer.data(msd0 + 143);
    const auto *msd0_144 = buffer.data(msd0 + 144);
    const auto *msd0_147 = buffer.data(msd0 + 147);
    const auto *msd0_149 = buffer.data(msd0 + 149);

    const auto *msd1_101 = buffer.data(msd1 + 101);
    const auto *msd1_102 = buffer.data(msd1 + 102);
    const auto *msd1_105 = buffer.data(msd1 + 105);
    const auto *msd1_107 = buffer.data(msd1 + 107);
    const auto *msd1_108 = buffer.data(msd1 + 108);
    const auto *msd1_111 = buffer.data(msd1 + 111);
    const auto *msd1_113 = buffer.data(msd1 + 113);
    const auto *msd1_117 = buffer.data(msd1 + 117);
    const auto *msd1_119 = buffer.data(msd1 + 119);
    const auto *msd1_120 = buffer.data(msd1 + 120);
    const auto *msd1_123 = buffer.data(msd1 + 123);
    const auto *msd1_124 = buffer.data(msd1 + 124);
    const auto *msd1_125 = buffer.data(msd1 + 125);
    const auto *msd1_126 = buffer.data(msd1 + 126);
    const auto *msd1_129 = buffer.data(msd1 + 129);
    const auto *msd1_131 = buffer.data(msd1 + 131);
    const auto *msd1_137 = buffer.data(msd1 + 137);
    const auto *msd1_138 = buffer.data(msd1 + 138);
    const auto *msd1_141 = buffer.data(msd1 + 141);
    const auto *msd1_143 = buffer.data(msd1 + 143);
    const auto *msd1_144 = buffer.data(msd1 + 144);
    const auto *msd1_147 = buffer.data(msd1 + 147);
    const auto *msd1_149 = buffer.data(msd1 + 149);

    const auto *msf_169 = buffer.data(msf + 169);
    const auto *msf_170 = buffer.data(msf + 170);
    const auto *msf_172 = buffer.data(msf + 172);
    const auto *msf_173 = buffer.data(msf + 173);
    const auto *msf_175 = buffer.data(msf + 175);
    const auto *msf_176 = buffer.data(msf + 176);
    const auto *msf_177 = buffer.data(msf + 177);
    const auto *msf_178 = buffer.data(msf + 178);
    const auto *msf_179 = buffer.data(msf + 179);
    const auto *msf_180 = buffer.data(msf + 180);
    const auto *msf_182 = buffer.data(msf + 182);
    const auto *msf_183 = buffer.data(msf + 183);
    const auto *msf_185 = buffer.data(msf + 185);
    const auto *msf_186 = buffer.data(msf + 186);
    const auto *msf_187 = buffer.data(msf + 187);
    const auto *msf_188 = buffer.data(msf + 188);
    const auto *msf_189 = buffer.data(msf + 189);
    const auto *msf_190 = buffer.data(msf + 190);
    const auto *msf_192 = buffer.data(msf + 192);
    const auto *msf_196 = buffer.data(msf + 196);
    const auto *msf_197 = buffer.data(msf + 197);
    const auto *msf_198 = buffer.data(msf + 198);
    const auto *msf_199 = buffer.data(msf + 199);
    const auto *msf_200 = buffer.data(msf + 200);
    const auto *msf_201 = buffer.data(msf + 201);
    const auto *msf_202 = buffer.data(msf + 202);
    const auto *msf_205 = buffer.data(msf + 205);
    const auto *msf_206 = buffer.data(msf + 206);
    const auto *msf_207 = buffer.data(msf + 207);
    const auto *msf_208 = buffer.data(msf + 208);
    const auto *msf_209 = buffer.data(msf + 209);
    const auto *msf_210 = buffer.data(msf + 210);
    const auto *msf_211 = buffer.data(msf + 211);
    const auto *msf_212 = buffer.data(msf + 212);
    const auto *msf_213 = buffer.data(msf + 213);
    const auto *msf_216 = buffer.data(msf + 216);
    const auto *msf_217 = buffer.data(msf + 217);
    const auto *msf_218 = buffer.data(msf + 218);
    const auto *msf_219 = buffer.data(msf + 219);
    const auto *msf_220 = buffer.data(msf + 220);
    const auto *msf_222 = buffer.data(msf + 222);
    const auto *msf_225 = buffer.data(msf + 225);
    const auto *msf_226 = buffer.data(msf + 226);
    const auto *msf_227 = buffer.data(msf + 227);
    const auto *msf_228 = buffer.data(msf + 228);
    const auto *msf_229 = buffer.data(msf + 229);
    const auto *msf_230 = buffer.data(msf + 230);
    const auto *msf_232 = buffer.data(msf + 232);
    const auto *msf_233 = buffer.data(msf + 233);
    const auto *msf_235 = buffer.data(msf + 235);
    const auto *msf_236 = buffer.data(msf + 236);
    const auto *msf_237 = buffer.data(msf + 237);
    const auto *msf_238 = buffer.data(msf + 238);
    const auto *msf_239 = buffer.data(msf + 239);
    const auto *msf_240 = buffer.data(msf + 240);
    const auto *msf_242 = buffer.data(msf + 242);
    const auto *msf_243 = buffer.data(msf + 243);
    const auto *msf_245 = buffer.data(msf + 245);
    const auto *msf_246 = buffer.data(msf + 246);
    const auto *msf_247 = buffer.data(msf + 247);
    const auto *msf_248 = buffer.data(msf + 248);
    const auto *msf_249 = buffer.data(msf + 249);

#pragma omp simd aligned(t_253, t_254, t_255, pc_x, pc_y, pc_z, lsf_109, lsf_119, lsf_170, \
                         msd0_101, msd0_102, msd1_101, msd1_102, msf_169, \
                         msf_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_253[k] = f_16 * lsf_119[k]
                   + f_3 * pc_y[k] * msf_169[k];

        t_254[k] = f_7 * lsf_109[k]
                   + f_1 * msd0_101[k]
                   - f_2 * msd1_101[k]
                   + f_3 * pc_z[k] * msf_169[k];

        t_255[k] = f_16 * lsf_170[k]
                   + f_1 * msd0_102[k]
                   - f_2 * msd1_102[k]
                   + f_3 * pc_x[k] * msf_170[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, t_259, pc_x, pc_y, pc_z, lsf_110, lsf_120, \
                         lsf_122, lsf_173, msd0_105, msd1_105, msf_170, msf_172, \
                         msf_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = f_14 * lsf_120[k]
                   + f_3 * pc_y[k] * msf_170[k];

        t_257[k] = f_8 * lsf_110[k]
                   + f_3 * pc_z[k] * msf_170[k];

        t_258[k] = f_16 * lsf_173[k]
                   + f_4 * msd0_105[k]
                   - f_5 * msd1_105[k]
                   + f_3 * pc_x[k] * msf_173[k];

        t_259[k] = f_14 * lsf_122[k]
                   + f_3 * pc_y[k] * msf_172[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pc_x, lsf_175, lsf_176, lsf_177, lsf_178, \
                         msd0_107, msd1_107, msf_175, msf_176, msf_177, \
                         msf_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_16 * lsf_175[k]
                   + f_4 * msd0_107[k]
                   - f_5 * msd1_107[k]
                   + f_3 * pc_x[k] * msf_175[k];

        t_261[k] = f_16 * lsf_176[k]
                   + f_3 * pc_x[k] * msf_176[k];

        t_262[k] = f_16 * lsf_177[k]
                   + f_3 * pc_x[k] * msf_177[k];

        t_263[k] = f_16 * lsf_178[k]
                   + f_3 * pc_x[k] * msf_178[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, pc_x, pc_y, pc_z, lsf_116, lsf_126, lsf_179, \
                         msd0_105, msd1_105, msf_176, msf_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_16 * lsf_179[k]
                   + f_3 * pc_x[k] * msf_179[k];

        t_265[k] = f_14 * lsf_126[k]
                   + f_1 * msd0_105[k]
                   - f_2 * msd1_105[k]
                   + f_3 * pc_y[k] * msf_176[k];

        t_266[k] = f_8 * lsf_116[k]
                   + f_3 * pc_z[k] * msf_176[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, pc_y, pc_z, lsf_119, lsf_128, lsf_129, msd0_107, \
                         msd1_107, msf_178, msf_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_14 * lsf_128[k]
                   + f_4 * msd0_107[k]
                   - f_5 * msd1_107[k]
                   + f_3 * pc_y[k] * msf_178[k];

        t_268[k] = f_14 * lsf_129[k]
                   + f_3 * pc_y[k] * msf_179[k];

        t_269[k] = f_8 * lsf_119[k]
                   + f_1 * msd0_107[k]
                   - f_2 * msd1_107[k]
                   + f_3 * pc_z[k] * msf_179[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, pc_x, pc_y, pc_z, lsf_120, lsf_130, lsf_180, \
                         msd0_108, msd1_108, msf_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_16 * lsf_180[k]
                   + f_1 * msd0_108[k]
                   - f_2 * msd1_108[k]
                   + f_3 * pc_x[k] * msf_180[k];

        t_271[k] = f_8 * lsf_130[k]
                   + f_3 * pc_y[k] * msf_180[k];

        t_272[k] = f_14 * lsf_120[k]
                   + f_3 * pc_z[k] * msf_180[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, pc_x, pc_y, lsf_132, lsf_183, lsf_185, msd0_111, \
                         msd0_113, msd1_111, msd1_113, msf_182, msf_183, \
                         msf_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = f_16 * lsf_183[k]
                   + f_4 * msd0_111[k]
                   - f_5 * msd1_111[k]
                   + f_3 * pc_x[k] * msf_183[k];

        t_274[k] = f_8 * lsf_132[k]
                   + f_3 * pc_y[k] * msf_182[k];

        t_275[k] = f_16 * lsf_185[k]
                   + f_4 * msd0_113[k]
                   - f_5 * msd1_113[k]
                   + f_3 * pc_x[k] * msf_185[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pc_x, lsf_186, lsf_187, lsf_188, lsf_189, \
                         msf_186, msf_187, msf_188, msf_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_16 * lsf_186[k]
                   + f_3 * pc_x[k] * msf_186[k];

        t_277[k] = f_16 * lsf_187[k]
                   + f_3 * pc_x[k] * msf_187[k];

        t_278[k] = f_16 * lsf_188[k]
                   + f_3 * pc_x[k] * msf_188[k];

        t_279[k] = f_16 * lsf_189[k]
                   + f_3 * pc_x[k] * msf_189[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, pc_y, pc_z, lsf_126, lsf_136, lsf_138, msd0_111, \
                         msd0_113, msd1_111, msd1_113, msf_186, \
                         msf_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_8 * lsf_136[k]
                   + f_1 * msd0_111[k]
                   - f_2 * msd1_111[k]
                   + f_3 * pc_y[k] * msf_186[k];

        t_281[k] = f_14 * lsf_126[k]
                   + f_3 * pc_z[k] * msf_186[k];

        t_282[k] = f_8 * lsf_138[k]
                   + f_4 * msd0_113[k]
                   - f_5 * msd1_113[k]
                   + f_3 * pc_y[k] * msf_188[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, t_286, pa_y, pc_y, pc_z, lsg0_210, lsf_129, \
                         lsf_139, lsf_140, lsg1_210, msd0_113, msd1_113, msf_189, \
                         msf_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_8 * lsf_139[k]
                   + f_3 * pc_y[k] * msf_189[k];

        t_284[k] = f_14 * lsf_129[k]
                   + f_1 * msd0_113[k]
                   - f_2 * msd1_113[k]
                   + f_3 * pc_z[k] * msf_189[k];

        t_285[k] = pa_y[k] * lsg0_210[k]
                   - f_6 * pc_y[k] * lsg1_210[k];

        t_286[k] = f_7 * lsf_140[k]
                   + f_3 * pc_y[k] * msf_190[k];
    }

#pragma omp simd aligned(t_287, t_288, t_289, t_290, pa_y, pc_y, pc_z, lsg0_213, lsg0_215, \
                         lsf_130, lsf_141, lsf_142, lsg1_213, lsg1_215, msf_190, \
                         msf_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_287[k] = f_16 * lsf_130[k]
                   + f_3 * pc_z[k] * msf_190[k];

        t_288[k] = pa_y[k] * lsg0_213[k]
                   + f_8 * lsf_141[k]
                   - f_6 * pc_y[k] * lsg1_213[k];

        t_289[k] = f_7 * lsf_142[k]
                   + f_3 * pc_y[k] * msf_192[k];

        t_290[k] = pa_y[k] * lsg0_215[k]
                   - f_6 * pc_y[k] * lsg1_215[k];
    }

#pragma omp simd aligned(t_291, t_292, t_293, t_294, pc_x, lsf_196, lsf_197, lsf_198, lsf_199, \
                         msf_196, msf_197, msf_198, msf_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_291[k] = f_16 * lsf_196[k]
                   + f_3 * pc_x[k] * msf_196[k];

        t_292[k] = f_16 * lsf_197[k]
                   + f_3 * pc_x[k] * msf_197[k];

        t_293[k] = f_16 * lsf_198[k]
                   + f_3 * pc_x[k] * msf_198[k];

        t_294[k] = f_16 * lsf_199[k]
                   + f_3 * pc_x[k] * msf_199[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, pc_y, pc_z, lsf_136, lsf_146, lsf_148, msd0_117, \
                         msd0_119, msd1_117, msd1_119, msf_196, \
                         msf_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = f_7 * lsf_146[k]
                   + f_1 * msd0_117[k]
                   - f_2 * msd1_117[k]
                   + f_3 * pc_y[k] * msf_196[k];

        t_296[k] = f_16 * lsf_136[k]
                   + f_3 * pc_z[k] * msf_196[k];

        t_297[k] = f_7 * lsf_148[k]
                   + f_4 * msd0_119[k]
                   - f_5 * msd1_119[k]
                   + f_3 * pc_y[k] * msf_198[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, t_301, pa_y, pc_x, pc_y, lsg0_224, lsf_149, \
                         lsf_200, lsg1_224, msd0_120, msd1_120, msf_199, \
                         msf_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = f_7 * lsf_149[k]
                   + f_3 * pc_y[k] * msf_199[k];

        t_299[k] = pa_y[k] * lsg0_224[k]
                   - f_6 * pc_y[k] * lsg1_224[k];

        t_300[k] = f_16 * lsf_200[k]
                   + f_1 * msd0_120[k]
                   - f_2 * msd1_120[k]
                   + f_3 * pc_x[k] * msf_200[k];

        t_301[k] = f_3 * pc_y[k] * msf_200[k];
    }

#pragma omp simd aligned(t_302, t_303, t_304, pc_y, pc_z, lsf_140, msd0_120, msd1_120, \
                         msf_200, msf_201, msf_202 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_302[k] = f_15 * lsf_140[k]
                   + f_3 * pc_z[k] * msf_200[k];

        t_303[k] = f_4 * msd0_120[k]
                   - f_5 * msd1_120[k]
                   + f_3 * pc_y[k] * msf_201[k];

        t_304[k] = f_3 * pc_y[k] * msf_202[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, pc_x, pc_y, lsf_205, lsf_206, lsf_207, \
                         msd0_125, msd1_125, msf_205, msf_206, \
                         msf_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = f_16 * lsf_205[k]
                   + f_4 * msd0_125[k]
                   - f_5 * msd1_125[k]
                   + f_3 * pc_x[k] * msf_205[k];

        t_306[k] = f_16 * lsf_206[k]
                   + f_3 * pc_x[k] * msf_206[k];

        t_307[k] = f_16 * lsf_207[k]
                   + f_3 * pc_x[k] * msf_207[k];

        t_308[k] = f_3 * pc_y[k] * msf_205[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, pc_x, pc_y, lsf_209, msd0_123, msd0_124, \
                         msd1_123, msd1_124, msf_206, msf_207, \
                         msf_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_16 * lsf_209[k]
                   + f_3 * pc_x[k] * msf_209[k];

        t_310[k] = f_1 * msd0_123[k]
                   - f_2 * msd1_123[k]
                   + f_3 * pc_y[k] * msf_206[k];

        t_311[k] = f_10 * msd0_124[k]
                   - f_11 * msd1_124[k]
                   + f_3 * pc_y[k] * msf_207[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, t_315, pc_x, pc_y, pc_z, lsf_149, lsf_210, \
                         msd0_125, msd0_126, msd1_125, msd1_126, msf_208, msf_209, \
                         msf_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = f_4 * msd0_125[k]
                   - f_5 * msd1_125[k]
                   + f_3 * pc_y[k] * msf_208[k];

        t_313[k] = f_3 * pc_y[k] * msf_209[k];

        t_314[k] = f_15 * lsf_149[k]
                   + f_1 * msd0_125[k]
                   - f_2 * msd1_125[k]
                   + f_3 * pc_z[k] * msf_209[k];

        t_315[k] = f_14 * lsf_210[k]
                   + f_1 * msd0_126[k]
                   - f_2 * msd1_126[k]
                   + f_3 * pc_x[k] * msf_210[k];
    }

#pragma omp simd aligned(t_316, t_317, t_318, t_319, pc_x, pc_y, pc_z, lsf_150, lsf_213, \
                         msd0_129, msd1_129, msf_210, msf_211, \
                         msf_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_316[k] = f_13 * lsf_150[k]
                   + f_3 * pc_y[k] * msf_210[k];

        t_317[k] = f_3 * pc_z[k] * msf_210[k];

        t_318[k] = f_14 * lsf_213[k]
                   + f_4 * msd0_129[k]
                   - f_5 * msd1_129[k]
                   + f_3 * pc_x[k] * msf_213[k];

        t_319[k] = f_3 * pc_z[k] * msf_211[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, pc_x, pc_z, lsf_216, lsf_218, msd0_126, \
                         msd1_126, msf_212, msf_213, msf_216, msf_218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = f_4 * msd0_126[k]
                   - f_5 * msd1_126[k]
                   + f_3 * pc_z[k] * msf_212[k];

        t_321[k] = f_14 * lsf_216[k]
                   + f_3 * pc_x[k] * msf_216[k];

        t_322[k] = f_3 * pc_z[k] * msf_213[k];

        t_323[k] = f_14 * lsf_218[k]
                   + f_3 * pc_x[k] * msf_218[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, t_328, pc_x, pc_y, pc_z, lsf_156, \
                         lsf_159, lsf_219, msd0_129, msd1_129, msf_216, msf_217, \
                         msf_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = f_14 * lsf_219[k]
                   + f_3 * pc_x[k] * msf_219[k];

        t_325[k] = f_13 * lsf_156[k]
                   + f_1 * msd0_129[k]
                   - f_2 * msd1_129[k]
                   + f_3 * pc_y[k] * msf_216[k];

        t_326[k] = f_3 * pc_z[k] * msf_216[k];

        t_327[k] = f_4 * msd0_129[k]
                   - f_5 * msd1_129[k]
                   + f_3 * pc_z[k] * msf_217[k];

        t_328[k] = f_13 * lsf_159[k]
                   + f_3 * pc_y[k] * msf_219[k];
    }

#pragma omp simd aligned(t_329, t_330, t_331, t_332, pa_z, pc_y, pc_z, lsg0_225, lsf_150, \
                         lsf_160, lsg1_225, msd0_131, msd1_131, msf_219, \
                         msf_220 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_329[k] = f_1 * msd0_131[k]
                   - f_2 * msd1_131[k]
                   + f_3 * pc_z[k] * msf_219[k];

        t_330[k] = pa_z[k] * lsg0_225[k]
                   - f_6 * pc_z[k] * lsg1_225[k];

        t_331[k] = f_15 * lsf_160[k]
                   + f_3 * pc_y[k] * msf_220[k];

        t_332[k] = f_7 * lsf_150[k]
                   + f_3 * pc_z[k] * msf_220[k];
    }

#pragma omp simd aligned(t_333, t_334, t_335, pa_z, pc_x, pc_y, pc_z, lsg0_228, lsf_162, \
                         lsf_225, lsg1_228, msd0_137, msd1_137, msf_222, \
                         msf_225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_333[k] = pa_z[k] * lsg0_228[k]
                   - f_6 * pc_z[k] * lsg1_228[k];

        t_334[k] = f_15 * lsf_162[k]
                   + f_3 * pc_y[k] * msf_222[k];

        t_335[k] = f_14 * lsf_225[k]
                   + f_4 * msd0_137[k]
                   - f_5 * msd1_137[k]
                   + f_3 * pc_x[k] * msf_225[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, t_339, pc_x, lsf_226, lsf_227, lsf_228, lsf_229, \
                         msf_226, msf_227, msf_228, msf_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = f_14 * lsf_226[k]
                   + f_3 * pc_x[k] * msf_226[k];

        t_337[k] = f_14 * lsf_227[k]
                   + f_3 * pc_x[k] * msf_227[k];

        t_338[k] = f_14 * lsf_228[k]
                   + f_3 * pc_x[k] * msf_228[k];

        t_339[k] = f_14 * lsf_229[k]
                   + f_3 * pc_x[k] * msf_229[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, pa_z, pc_y, pc_z, lsg0_235, lsf_156, lsf_168, \
                         lsg1_235, msd0_137, msd1_137, msf_226, \
                         msf_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = pa_z[k] * lsg0_235[k]
                   - f_6 * pc_z[k] * lsg1_235[k];

        t_341[k] = f_7 * lsf_156[k]
                   + f_3 * pc_z[k] * msf_226[k];

        t_342[k] = f_15 * lsf_168[k]
                   + f_4 * msd0_137[k]
                   - f_5 * msd1_137[k]
                   + f_3 * pc_y[k] * msf_228[k];
    }

#pragma omp simd aligned(t_343, t_344, t_345, pc_x, pc_y, pc_z, lsf_159, lsf_169, lsf_230, \
                         msd0_137, msd0_138, msd1_137, msd1_138, msf_229, \
                         msf_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_343[k] = f_15 * lsf_169[k]
                   + f_3 * pc_y[k] * msf_229[k];

        t_344[k] = f_7 * lsf_159[k]
                   + f_1 * msd0_137[k]
                   - f_2 * msd1_137[k]
                   + f_3 * pc_z[k] * msf_229[k];

        t_345[k] = f_14 * lsf_230[k]
                   + f_1 * msd0_138[k]
                   - f_2 * msd1_138[k]
                   + f_3 * pc_x[k] * msf_230[k];
    }

#pragma omp simd aligned(t_346, t_347, t_348, t_349, pc_x, pc_y, pc_z, lsf_160, lsf_170, \
                         lsf_172, lsf_233, msd0_141, msd1_141, msf_230, msf_232, \
                         msf_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_346[k] = f_16 * lsf_170[k]
                   + f_3 * pc_y[k] * msf_230[k];

        t_347[k] = f_8 * lsf_160[k]
                   + f_3 * pc_z[k] * msf_230[k];

        t_348[k] = f_14 * lsf_233[k]
                   + f_4 * msd0_141[k]
                   - f_5 * msd1_141[k]
                   + f_3 * pc_x[k] * msf_233[k];

        t_349[k] = f_16 * lsf_172[k]
                   + f_3 * pc_y[k] * msf_232[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, pc_x, lsf_235, lsf_236, lsf_237, lsf_238, \
                         msd0_143, msd1_143, msf_235, msf_236, msf_237, \
                         msf_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = f_14 * lsf_235[k]
                   + f_4 * msd0_143[k]
                   - f_5 * msd1_143[k]
                   + f_3 * pc_x[k] * msf_235[k];

        t_351[k] = f_14 * lsf_236[k]
                   + f_3 * pc_x[k] * msf_236[k];

        t_352[k] = f_14 * lsf_237[k]
                   + f_3 * pc_x[k] * msf_237[k];

        t_353[k] = f_14 * lsf_238[k]
                   + f_3 * pc_x[k] * msf_238[k];
    }

#pragma omp simd aligned(t_354, t_355, t_356, pc_x, pc_y, pc_z, lsf_166, lsf_176, lsf_239, \
                         msd0_141, msd1_141, msf_236, msf_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_354[k] = f_14 * lsf_239[k]
                   + f_3 * pc_x[k] * msf_239[k];

        t_355[k] = f_16 * lsf_176[k]
                   + f_1 * msd0_141[k]
                   - f_2 * msd1_141[k]
                   + f_3 * pc_y[k] * msf_236[k];

        t_356[k] = f_8 * lsf_166[k]
                   + f_3 * pc_z[k] * msf_236[k];
    }

#pragma omp simd aligned(t_357, t_358, t_359, pc_y, pc_z, lsf_169, lsf_178, lsf_179, msd0_143, \
                         msd1_143, msf_238, msf_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_357[k] = f_16 * lsf_178[k]
                   + f_4 * msd0_143[k]
                   - f_5 * msd1_143[k]
                   + f_3 * pc_y[k] * msf_238[k];

        t_358[k] = f_16 * lsf_179[k]
                   + f_3 * pc_y[k] * msf_239[k];

        t_359[k] = f_8 * lsf_169[k]
                   + f_1 * msd0_143[k]
                   - f_2 * msd1_143[k]
                   + f_3 * pc_z[k] * msf_239[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, pc_x, pc_y, pc_z, lsf_170, lsf_180, lsf_240, \
                         msd0_144, msd1_144, msf_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = f_14 * lsf_240[k]
                   + f_1 * msd0_144[k]
                   - f_2 * msd1_144[k]
                   + f_3 * pc_x[k] * msf_240[k];

        t_361[k] = f_14 * lsf_180[k]
                   + f_3 * pc_y[k] * msf_240[k];

        t_362[k] = f_14 * lsf_170[k]
                   + f_3 * pc_z[k] * msf_240[k];
    }

#pragma omp simd aligned(t_363, t_364, t_365, pc_x, pc_y, lsf_182, lsf_243, lsf_245, msd0_147, \
                         msd0_149, msd1_147, msd1_149, msf_242, msf_243, \
                         msf_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = f_14 * lsf_243[k]
                   + f_4 * msd0_147[k]
                   - f_5 * msd1_147[k]
                   + f_3 * pc_x[k] * msf_243[k];

        t_364[k] = f_14 * lsf_182[k]
                   + f_3 * pc_y[k] * msf_242[k];

        t_365[k] = f_14 * lsf_245[k]
                   + f_4 * msd0_149[k]
                   - f_5 * msd1_149[k]
                   + f_3 * pc_x[k] * msf_245[k];
    }

#pragma omp simd aligned(t_366, t_367, t_368, t_369, pc_x, lsf_246, lsf_247, lsf_248, lsf_249, \
                         msf_246, msf_247, msf_248, msf_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = f_14 * lsf_246[k]
                   + f_3 * pc_x[k] * msf_246[k];

        t_367[k] = f_14 * lsf_247[k]
                   + f_3 * pc_x[k] * msf_247[k];

        t_368[k] = f_14 * lsf_248[k]
                   + f_3 * pc_x[k] * msf_248[k];

        t_369[k] = f_14 * lsf_249[k]
                   + f_3 * pc_x[k] * msf_249[k];
    }
}

static auto
compute_prim_msg_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsg0,
                                                          const size_t lsf, const size_t lsg1,
                                                          const size_t msd0, const size_t msd1,
                                                          const size_t msf, const size_t ncols,
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
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 3.5 / q;
    const auto f_13 = 3.0 / q;
    const auto f_14 = 1.5 / q;
    const auto f_15 = 2.5 / q;
    const auto f_16 = 2.0 / q;

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
    auto *t_477 = buffer.data(target + 477);
    auto *t_478 = buffer.data(target + 478);
    auto *t_479 = buffer.data(target + 479);
    auto *t_480 = buffer.data(target + 480);
    auto *t_481 = buffer.data(target + 481);
    auto *t_482 = buffer.data(target + 482);
    auto *t_483 = buffer.data(target + 483);
    auto *t_484 = buffer.data(target + 484);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *lsg0_300 = buffer.data(lsg0 + 300);
    const auto *lsg0_303 = buffer.data(lsg0 + 303);
    const auto *lsg0_305 = buffer.data(lsg0 + 305);
    const auto *lsg0_314 = buffer.data(lsg0 + 314);
    const auto *lsg0_315 = buffer.data(lsg0 + 315);
    const auto *lsg0_318 = buffer.data(lsg0 + 318);
    const auto *lsg0_325 = buffer.data(lsg0 + 325);

    const auto *lsf_176 = buffer.data(lsf + 176);
    const auto *lsf_179 = buffer.data(lsf + 179);
    const auto *lsf_180 = buffer.data(lsf + 180);
    const auto *lsf_186 = buffer.data(lsf + 186);
    const auto *lsf_188 = buffer.data(lsf + 188);
    const auto *lsf_189 = buffer.data(lsf + 189);
    const auto *lsf_190 = buffer.data(lsf + 190);
    const auto *lsf_192 = buffer.data(lsf + 192);
    const auto *lsf_196 = buffer.data(lsf + 196);
    const auto *lsf_198 = buffer.data(lsf + 198);
    const auto *lsf_199 = buffer.data(lsf + 199);
    const auto *lsf_200 = buffer.data(lsf + 200);
    const auto *lsf_201 = buffer.data(lsf + 201);
    const auto *lsf_202 = buffer.data(lsf + 202);
    const auto *lsf_206 = buffer.data(lsf + 206);
    const auto *lsf_208 = buffer.data(lsf + 208);
    const auto *lsf_209 = buffer.data(lsf + 209);
    const auto *lsf_210 = buffer.data(lsf + 210);
    const auto *lsf_216 = buffer.data(lsf + 216);
    const auto *lsf_219 = buffer.data(lsf + 219);
    const auto *lsf_220 = buffer.data(lsf + 220);
    const auto *lsf_222 = buffer.data(lsf + 222);
    const auto *lsf_226 = buffer.data(lsf + 226);
    const auto *lsf_228 = buffer.data(lsf + 228);
    const auto *lsf_229 = buffer.data(lsf + 229);
    const auto *lsf_230 = buffer.data(lsf + 230);
    const auto *lsf_232 = buffer.data(lsf + 232);
    const auto *lsf_236 = buffer.data(lsf + 236);
    const auto *lsf_238 = buffer.data(lsf + 238);
    const auto *lsf_239 = buffer.data(lsf + 239);
    const auto *lsf_240 = buffer.data(lsf + 240);
    const auto *lsf_242 = buffer.data(lsf + 242);
    const auto *lsf_246 = buffer.data(lsf + 246);
    const auto *lsf_248 = buffer.data(lsf + 248);
    const auto *lsf_249 = buffer.data(lsf + 249);
    const auto *lsf_250 = buffer.data(lsf + 250);
    const auto *lsf_252 = buffer.data(lsf + 252);
    const auto *lsf_253 = buffer.data(lsf + 253);
    const auto *lsf_255 = buffer.data(lsf + 255);
    const auto *lsf_256 = buffer.data(lsf + 256);
    const auto *lsf_257 = buffer.data(lsf + 257);
    const auto *lsf_258 = buffer.data(lsf + 258);
    const auto *lsf_259 = buffer.data(lsf + 259);
    const auto *lsf_266 = buffer.data(lsf + 266);
    const auto *lsf_267 = buffer.data(lsf + 267);
    const auto *lsf_268 = buffer.data(lsf + 268);
    const auto *lsf_269 = buffer.data(lsf + 269);
    const auto *lsf_270 = buffer.data(lsf + 270);
    const auto *lsf_275 = buffer.data(lsf + 275);
    const auto *lsf_276 = buffer.data(lsf + 276);
    const auto *lsf_277 = buffer.data(lsf + 277);
    const auto *lsf_279 = buffer.data(lsf + 279);
    const auto *lsf_280 = buffer.data(lsf + 280);
    const auto *lsf_283 = buffer.data(lsf + 283);
    const auto *lsf_286 = buffer.data(lsf + 286);
    const auto *lsf_288 = buffer.data(lsf + 288);
    const auto *lsf_289 = buffer.data(lsf + 289);
    const auto *lsf_295 = buffer.data(lsf + 295);
    const auto *lsf_296 = buffer.data(lsf + 296);
    const auto *lsf_297 = buffer.data(lsf + 297);
    const auto *lsf_298 = buffer.data(lsf + 298);
    const auto *lsf_299 = buffer.data(lsf + 299);
    const auto *lsf_300 = buffer.data(lsf + 300);
    const auto *lsf_303 = buffer.data(lsf + 303);
    const auto *lsf_305 = buffer.data(lsf + 305);
    const auto *lsf_306 = buffer.data(lsf + 306);
    const auto *lsf_307 = buffer.data(lsf + 307);
    const auto *lsf_308 = buffer.data(lsf + 308);
    const auto *lsf_309 = buffer.data(lsf + 309);
    const auto *lsf_310 = buffer.data(lsf + 310);
    const auto *lsf_313 = buffer.data(lsf + 313);
    const auto *lsf_315 = buffer.data(lsf + 315);
    const auto *lsf_316 = buffer.data(lsf + 316);
    const auto *lsf_317 = buffer.data(lsf + 317);
    const auto *lsf_318 = buffer.data(lsf + 318);
    const auto *lsf_319 = buffer.data(lsf + 319);
    const auto *lsf_320 = buffer.data(lsf + 320);
    const auto *lsf_323 = buffer.data(lsf + 323);

    const auto *lsg1_300 = buffer.data(lsg1 + 300);
    const auto *lsg1_303 = buffer.data(lsg1 + 303);
    const auto *lsg1_305 = buffer.data(lsg1 + 305);
    const auto *lsg1_314 = buffer.data(lsg1 + 314);
    const auto *lsg1_315 = buffer.data(lsg1 + 315);
    const auto *lsg1_318 = buffer.data(lsg1 + 318);
    const auto *lsg1_325 = buffer.data(lsg1 + 325);

    const auto *msd0_147 = buffer.data(msd0 + 147);
    const auto *msd0_149 = buffer.data(msd0 + 149);
    const auto *msd0_150 = buffer.data(msd0 + 150);
    const auto *msd0_153 = buffer.data(msd0 + 153);
    const auto *msd0_155 = buffer.data(msd0 + 155);
    const auto *msd0_159 = buffer.data(msd0 + 159);
    const auto *msd0_161 = buffer.data(msd0 + 161);
    const auto *msd0_162 = buffer.data(msd0 + 162);
    const auto *msd0_165 = buffer.data(msd0 + 165);
    const auto *msd0_166 = buffer.data(msd0 + 166);
    const auto *msd0_167 = buffer.data(msd0 + 167);
    const auto *msd0_168 = buffer.data(msd0 + 168);
    const auto *msd0_171 = buffer.data(msd0 + 171);
    const auto *msd0_173 = buffer.data(msd0 + 173);
    const auto *msd0_179 = buffer.data(msd0 + 179);
    const auto *msd0_180 = buffer.data(msd0 + 180);
    const auto *msd0_183 = buffer.data(msd0 + 183);
    const auto *msd0_185 = buffer.data(msd0 + 185);
    const auto *msd0_186 = buffer.data(msd0 + 186);
    const auto *msd0_189 = buffer.data(msd0 + 189);
    const auto *msd0_191 = buffer.data(msd0 + 191);
    const auto *msd0_192 = buffer.data(msd0 + 192);
    const auto *msd0_195 = buffer.data(msd0 + 195);

    const auto *msd1_147 = buffer.data(msd1 + 147);
    const auto *msd1_149 = buffer.data(msd1 + 149);
    const auto *msd1_150 = buffer.data(msd1 + 150);
    const auto *msd1_153 = buffer.data(msd1 + 153);
    const auto *msd1_155 = buffer.data(msd1 + 155);
    const auto *msd1_159 = buffer.data(msd1 + 159);
    const auto *msd1_161 = buffer.data(msd1 + 161);
    const auto *msd1_162 = buffer.data(msd1 + 162);
    const auto *msd1_165 = buffer.data(msd1 + 165);
    const auto *msd1_166 = buffer.data(msd1 + 166);
    const auto *msd1_167 = buffer.data(msd1 + 167);
    const auto *msd1_168 = buffer.data(msd1 + 168);
    const auto *msd1_171 = buffer.data(msd1 + 171);
    const auto *msd1_173 = buffer.data(msd1 + 173);
    const auto *msd1_179 = buffer.data(msd1 + 179);
    const auto *msd1_180 = buffer.data(msd1 + 180);
    const auto *msd1_183 = buffer.data(msd1 + 183);
    const auto *msd1_185 = buffer.data(msd1 + 185);
    const auto *msd1_186 = buffer.data(msd1 + 186);
    const auto *msd1_189 = buffer.data(msd1 + 189);
    const auto *msd1_191 = buffer.data(msd1 + 191);
    const auto *msd1_192 = buffer.data(msd1 + 192);
    const auto *msd1_195 = buffer.data(msd1 + 195);

    const auto *msf_246 = buffer.data(msf + 246);
    const auto *msf_248 = buffer.data(msf + 248);
    const auto *msf_249 = buffer.data(msf + 249);
    const auto *msf_250 = buffer.data(msf + 250);
    const auto *msf_252 = buffer.data(msf + 252);
    const auto *msf_253 = buffer.data(msf + 253);
    const auto *msf_255 = buffer.data(msf + 255);
    const auto *msf_256 = buffer.data(msf + 256);
    const auto *msf_257 = buffer.data(msf + 257);
    const auto *msf_258 = buffer.data(msf + 258);
    const auto *msf_259 = buffer.data(msf + 259);
    const auto *msf_260 = buffer.data(msf + 260);
    const auto *msf_262 = buffer.data(msf + 262);
    const auto *msf_266 = buffer.data(msf + 266);
    const auto *msf_267 = buffer.data(msf + 267);
    const auto *msf_268 = buffer.data(msf + 268);
    const auto *msf_269 = buffer.data(msf + 269);
    const auto *msf_270 = buffer.data(msf + 270);
    const auto *msf_271 = buffer.data(msf + 271);
    const auto *msf_272 = buffer.data(msf + 272);
    const auto *msf_275 = buffer.data(msf + 275);
    const auto *msf_276 = buffer.data(msf + 276);
    const auto *msf_277 = buffer.data(msf + 277);
    const auto *msf_278 = buffer.data(msf + 278);
    const auto *msf_279 = buffer.data(msf + 279);
    const auto *msf_280 = buffer.data(msf + 280);
    const auto *msf_281 = buffer.data(msf + 281);
    const auto *msf_282 = buffer.data(msf + 282);
    const auto *msf_283 = buffer.data(msf + 283);
    const auto *msf_286 = buffer.data(msf + 286);
    const auto *msf_287 = buffer.data(msf + 287);
    const auto *msf_288 = buffer.data(msf + 288);
    const auto *msf_289 = buffer.data(msf + 289);
    const auto *msf_290 = buffer.data(msf + 290);
    const auto *msf_292 = buffer.data(msf + 292);
    const auto *msf_295 = buffer.data(msf + 295);
    const auto *msf_296 = buffer.data(msf + 296);
    const auto *msf_297 = buffer.data(msf + 297);
    const auto *msf_298 = buffer.data(msf + 298);
    const auto *msf_299 = buffer.data(msf + 299);
    const auto *msf_300 = buffer.data(msf + 300);
    const auto *msf_302 = buffer.data(msf + 302);
    const auto *msf_303 = buffer.data(msf + 303);
    const auto *msf_305 = buffer.data(msf + 305);
    const auto *msf_306 = buffer.data(msf + 306);
    const auto *msf_307 = buffer.data(msf + 307);
    const auto *msf_308 = buffer.data(msf + 308);
    const auto *msf_309 = buffer.data(msf + 309);
    const auto *msf_310 = buffer.data(msf + 310);
    const auto *msf_312 = buffer.data(msf + 312);
    const auto *msf_313 = buffer.data(msf + 313);
    const auto *msf_315 = buffer.data(msf + 315);
    const auto *msf_316 = buffer.data(msf + 316);
    const auto *msf_317 = buffer.data(msf + 317);
    const auto *msf_318 = buffer.data(msf + 318);
    const auto *msf_319 = buffer.data(msf + 319);
    const auto *msf_320 = buffer.data(msf + 320);
    const auto *msf_322 = buffer.data(msf + 322);
    const auto *msf_323 = buffer.data(msf + 323);

#pragma omp simd aligned(t_370, t_371, t_372, pc_y, pc_z, lsf_176, lsf_186, lsf_188, msd0_147, \
                         msd0_149, msd1_147, msd1_149, msf_246, \
                         msf_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = f_14 * lsf_186[k]
                   + f_1 * msd0_147[k]
                   - f_2 * msd1_147[k]
                   + f_3 * pc_y[k] * msf_246[k];

        t_371[k] = f_14 * lsf_176[k]
                   + f_3 * pc_z[k] * msf_246[k];

        t_372[k] = f_14 * lsf_188[k]
                   + f_4 * msd0_149[k]
                   - f_5 * msd1_149[k]
                   + f_3 * pc_y[k] * msf_248[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, pc_x, pc_y, pc_z, lsf_179, lsf_189, lsf_250, \
                         msd0_149, msd0_150, msd1_149, msd1_150, msf_249, \
                         msf_250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = f_14 * lsf_189[k]
                   + f_3 * pc_y[k] * msf_249[k];

        t_374[k] = f_14 * lsf_179[k]
                   + f_1 * msd0_149[k]
                   - f_2 * msd1_149[k]
                   + f_3 * pc_z[k] * msf_249[k];

        t_375[k] = f_14 * lsf_250[k]
                   + f_1 * msd0_150[k]
                   - f_2 * msd1_150[k]
                   + f_3 * pc_x[k] * msf_250[k];
    }

#pragma omp simd aligned(t_376, t_377, t_378, t_379, pc_x, pc_y, pc_z, lsf_180, lsf_190, \
                         lsf_192, lsf_253, msd0_153, msd1_153, msf_250, msf_252, \
                         msf_253 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_376[k] = f_8 * lsf_190[k]
                   + f_3 * pc_y[k] * msf_250[k];

        t_377[k] = f_16 * lsf_180[k]
                   + f_3 * pc_z[k] * msf_250[k];

        t_378[k] = f_14 * lsf_253[k]
                   + f_4 * msd0_153[k]
                   - f_5 * msd1_153[k]
                   + f_3 * pc_x[k] * msf_253[k];

        t_379[k] = f_8 * lsf_192[k]
                   + f_3 * pc_y[k] * msf_252[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, pc_x, lsf_255, lsf_256, lsf_257, lsf_258, \
                         msd0_155, msd1_155, msf_255, msf_256, msf_257, \
                         msf_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = f_14 * lsf_255[k]
                   + f_4 * msd0_155[k]
                   - f_5 * msd1_155[k]
                   + f_3 * pc_x[k] * msf_255[k];

        t_381[k] = f_14 * lsf_256[k]
                   + f_3 * pc_x[k] * msf_256[k];

        t_382[k] = f_14 * lsf_257[k]
                   + f_3 * pc_x[k] * msf_257[k];

        t_383[k] = f_14 * lsf_258[k]
                   + f_3 * pc_x[k] * msf_258[k];
    }

#pragma omp simd aligned(t_384, t_385, t_386, pc_x, pc_y, pc_z, lsf_186, lsf_196, lsf_259, \
                         msd0_153, msd1_153, msf_256, msf_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_384[k] = f_14 * lsf_259[k]
                   + f_3 * pc_x[k] * msf_259[k];

        t_385[k] = f_8 * lsf_196[k]
                   + f_1 * msd0_153[k]
                   - f_2 * msd1_153[k]
                   + f_3 * pc_y[k] * msf_256[k];

        t_386[k] = f_16 * lsf_186[k]
                   + f_3 * pc_z[k] * msf_256[k];
    }

#pragma omp simd aligned(t_387, t_388, t_389, t_390, pa_y, pc_y, pc_z, lsg0_300, lsf_189, \
                         lsf_198, lsf_199, lsg1_300, msd0_155, msd1_155, msf_258, \
                         msf_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_387[k] = f_8 * lsf_198[k]
                   + f_4 * msd0_155[k]
                   - f_5 * msd1_155[k]
                   + f_3 * pc_y[k] * msf_258[k];

        t_388[k] = f_8 * lsf_199[k]
                   + f_3 * pc_y[k] * msf_259[k];

        t_389[k] = f_16 * lsf_189[k]
                   + f_1 * msd0_155[k]
                   - f_2 * msd1_155[k]
                   + f_3 * pc_z[k] * msf_259[k];

        t_390[k] = pa_y[k] * lsg0_300[k]
                   - f_6 * pc_y[k] * lsg1_300[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, t_394, pa_y, pc_y, pc_z, lsg0_303, lsf_190, \
                         lsf_200, lsf_201, lsf_202, lsg1_303, msf_260, \
                         msf_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = f_7 * lsf_200[k]
                   + f_3 * pc_y[k] * msf_260[k];

        t_392[k] = f_15 * lsf_190[k]
                   + f_3 * pc_z[k] * msf_260[k];

        t_393[k] = pa_y[k] * lsg0_303[k]
                   + f_8 * lsf_201[k]
                   - f_6 * pc_y[k] * lsg1_303[k];

        t_394[k] = f_7 * lsf_202[k]
                   + f_3 * pc_y[k] * msf_262[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, pa_y, pc_x, pc_y, lsg0_305, lsf_266, \
                         lsf_267, lsf_268, lsg1_305, msf_266, msf_267, \
                         msf_268 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = pa_y[k] * lsg0_305[k]
                   - f_6 * pc_y[k] * lsg1_305[k];

        t_396[k] = f_14 * lsf_266[k]
                   + f_3 * pc_x[k] * msf_266[k];

        t_397[k] = f_14 * lsf_267[k]
                   + f_3 * pc_x[k] * msf_267[k];

        t_398[k] = f_14 * lsf_268[k]
                   + f_3 * pc_x[k] * msf_268[k];
    }

#pragma omp simd aligned(t_399, t_400, t_401, pc_x, pc_y, pc_z, lsf_196, lsf_206, lsf_269, \
                         msd0_159, msd1_159, msf_266, msf_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_399[k] = f_14 * lsf_269[k]
                   + f_3 * pc_x[k] * msf_269[k];

        t_400[k] = f_7 * lsf_206[k]
                   + f_1 * msd0_159[k]
                   - f_2 * msd1_159[k]
                   + f_3 * pc_y[k] * msf_266[k];

        t_401[k] = f_15 * lsf_196[k]
                   + f_3 * pc_z[k] * msf_266[k];
    }

#pragma omp simd aligned(t_402, t_403, t_404, pa_y, pc_y, lsg0_314, lsf_208, lsf_209, \
                         lsg1_314, msd0_161, msd1_161, msf_268, \
                         msf_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_402[k] = f_7 * lsf_208[k]
                   + f_4 * msd0_161[k]
                   - f_5 * msd1_161[k]
                   + f_3 * pc_y[k] * msf_268[k];

        t_403[k] = f_7 * lsf_209[k]
                   + f_3 * pc_y[k] * msf_269[k];

        t_404[k] = pa_y[k] * lsg0_314[k]
                   - f_6 * pc_y[k] * lsg1_314[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, pc_x, pc_y, pc_z, lsf_200, \
                         lsf_270, msd0_162, msd1_162, msf_270, msf_271, \
                         msf_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = f_14 * lsf_270[k]
                   + f_1 * msd0_162[k]
                   - f_2 * msd1_162[k]
                   + f_3 * pc_x[k] * msf_270[k];

        t_406[k] = f_3 * pc_y[k] * msf_270[k];

        t_407[k] = f_13 * lsf_200[k]
                   + f_3 * pc_z[k] * msf_270[k];

        t_408[k] = f_4 * msd0_162[k]
                   - f_5 * msd1_162[k]
                   + f_3 * pc_y[k] * msf_271[k];

        t_409[k] = f_3 * pc_y[k] * msf_272[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, pc_x, pc_y, lsf_275, lsf_276, lsf_277, \
                         msd0_167, msd1_167, msf_275, msf_276, \
                         msf_277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = f_14 * lsf_275[k]
                   + f_4 * msd0_167[k]
                   - f_5 * msd1_167[k]
                   + f_3 * pc_x[k] * msf_275[k];

        t_411[k] = f_14 * lsf_276[k]
                   + f_3 * pc_x[k] * msf_276[k];

        t_412[k] = f_14 * lsf_277[k]
                   + f_3 * pc_x[k] * msf_277[k];

        t_413[k] = f_3 * pc_y[k] * msf_275[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, pc_x, pc_y, lsf_279, msd0_165, msd0_166, \
                         msd1_165, msd1_166, msf_276, msf_277, \
                         msf_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_14 * lsf_279[k]
                   + f_3 * pc_x[k] * msf_279[k];

        t_415[k] = f_1 * msd0_165[k]
                   - f_2 * msd1_165[k]
                   + f_3 * pc_y[k] * msf_276[k];

        t_416[k] = f_10 * msd0_166[k]
                   - f_11 * msd1_166[k]
                   + f_3 * pc_y[k] * msf_277[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, t_420, pc_x, pc_y, pc_z, lsf_209, lsf_280, \
                         msd0_167, msd0_168, msd1_167, msd1_168, msf_278, msf_279, \
                         msf_280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = f_4 * msd0_167[k]
                   - f_5 * msd1_167[k]
                   + f_3 * pc_y[k] * msf_278[k];

        t_418[k] = f_3 * pc_y[k] * msf_279[k];

        t_419[k] = f_13 * lsf_209[k]
                   + f_1 * msd0_167[k]
                   - f_2 * msd1_167[k]
                   + f_3 * pc_z[k] * msf_279[k];

        t_420[k] = f_8 * lsf_280[k]
                   + f_1 * msd0_168[k]
                   - f_2 * msd1_168[k]
                   + f_3 * pc_x[k] * msf_280[k];
    }

#pragma omp simd aligned(t_421, t_422, t_423, t_424, pc_x, pc_y, pc_z, lsf_210, lsf_283, \
                         msd0_171, msd1_171, msf_280, msf_281, \
                         msf_283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_421[k] = f_12 * lsf_210[k]
                   + f_3 * pc_y[k] * msf_280[k];

        t_422[k] = f_3 * pc_z[k] * msf_280[k];

        t_423[k] = f_8 * lsf_283[k]
                   + f_4 * msd0_171[k]
                   - f_5 * msd1_171[k]
                   + f_3 * pc_x[k] * msf_283[k];

        t_424[k] = f_3 * pc_z[k] * msf_281[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, pc_x, pc_z, lsf_286, lsf_288, msd0_168, \
                         msd1_168, msf_282, msf_283, msf_286, msf_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = f_4 * msd0_168[k]
                   - f_5 * msd1_168[k]
                   + f_3 * pc_z[k] * msf_282[k];

        t_426[k] = f_8 * lsf_286[k]
                   + f_3 * pc_x[k] * msf_286[k];

        t_427[k] = f_3 * pc_z[k] * msf_283[k];

        t_428[k] = f_8 * lsf_288[k]
                   + f_3 * pc_x[k] * msf_288[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, t_432, t_433, pc_x, pc_y, pc_z, lsf_216, \
                         lsf_219, lsf_289, msd0_171, msd1_171, msf_286, msf_287, \
                         msf_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = f_8 * lsf_289[k]
                   + f_3 * pc_x[k] * msf_289[k];

        t_430[k] = f_12 * lsf_216[k]
                   + f_1 * msd0_171[k]
                   - f_2 * msd1_171[k]
                   + f_3 * pc_y[k] * msf_286[k];

        t_431[k] = f_3 * pc_z[k] * msf_286[k];

        t_432[k] = f_4 * msd0_171[k]
                   - f_5 * msd1_171[k]
                   + f_3 * pc_z[k] * msf_287[k];

        t_433[k] = f_12 * lsf_219[k]
                   + f_3 * pc_y[k] * msf_289[k];
    }

#pragma omp simd aligned(t_434, t_435, t_436, t_437, pa_z, pc_y, pc_z, lsg0_315, lsf_210, \
                         lsf_220, lsg1_315, msd0_173, msd1_173, msf_289, \
                         msf_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_434[k] = f_1 * msd0_173[k]
                   - f_2 * msd1_173[k]
                   + f_3 * pc_z[k] * msf_289[k];

        t_435[k] = pa_z[k] * lsg0_315[k]
                   - f_6 * pc_z[k] * lsg1_315[k];

        t_436[k] = f_13 * lsf_220[k]
                   + f_3 * pc_y[k] * msf_290[k];

        t_437[k] = f_7 * lsf_210[k]
                   + f_3 * pc_z[k] * msf_290[k];
    }

#pragma omp simd aligned(t_438, t_439, t_440, pa_z, pc_x, pc_y, pc_z, lsg0_318, lsf_222, \
                         lsf_295, lsg1_318, msd0_179, msd1_179, msf_292, \
                         msf_295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_438[k] = pa_z[k] * lsg0_318[k]
                   - f_6 * pc_z[k] * lsg1_318[k];

        t_439[k] = f_13 * lsf_222[k]
                   + f_3 * pc_y[k] * msf_292[k];

        t_440[k] = f_8 * lsf_295[k]
                   + f_4 * msd0_179[k]
                   - f_5 * msd1_179[k]
                   + f_3 * pc_x[k] * msf_295[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, t_444, pc_x, lsf_296, lsf_297, lsf_298, lsf_299, \
                         msf_296, msf_297, msf_298, msf_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = f_8 * lsf_296[k]
                   + f_3 * pc_x[k] * msf_296[k];

        t_442[k] = f_8 * lsf_297[k]
                   + f_3 * pc_x[k] * msf_297[k];

        t_443[k] = f_8 * lsf_298[k]
                   + f_3 * pc_x[k] * msf_298[k];

        t_444[k] = f_8 * lsf_299[k]
                   + f_3 * pc_x[k] * msf_299[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, pa_z, pc_y, pc_z, lsg0_325, lsf_216, lsf_228, \
                         lsg1_325, msd0_179, msd1_179, msf_296, \
                         msf_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = pa_z[k] * lsg0_325[k]
                   - f_6 * pc_z[k] * lsg1_325[k];

        t_446[k] = f_7 * lsf_216[k]
                   + f_3 * pc_z[k] * msf_296[k];

        t_447[k] = f_13 * lsf_228[k]
                   + f_4 * msd0_179[k]
                   - f_5 * msd1_179[k]
                   + f_3 * pc_y[k] * msf_298[k];
    }

#pragma omp simd aligned(t_448, t_449, t_450, pc_x, pc_y, pc_z, lsf_219, lsf_229, lsf_300, \
                         msd0_179, msd0_180, msd1_179, msd1_180, msf_299, \
                         msf_300 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_448[k] = f_13 * lsf_229[k]
                   + f_3 * pc_y[k] * msf_299[k];

        t_449[k] = f_7 * lsf_219[k]
                   + f_1 * msd0_179[k]
                   - f_2 * msd1_179[k]
                   + f_3 * pc_z[k] * msf_299[k];

        t_450[k] = f_8 * lsf_300[k]
                   + f_1 * msd0_180[k]
                   - f_2 * msd1_180[k]
                   + f_3 * pc_x[k] * msf_300[k];
    }

#pragma omp simd aligned(t_451, t_452, t_453, t_454, pc_x, pc_y, pc_z, lsf_220, lsf_230, \
                         lsf_232, lsf_303, msd0_183, msd1_183, msf_300, msf_302, \
                         msf_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_451[k] = f_15 * lsf_230[k]
                   + f_3 * pc_y[k] * msf_300[k];

        t_452[k] = f_8 * lsf_220[k]
                   + f_3 * pc_z[k] * msf_300[k];

        t_453[k] = f_8 * lsf_303[k]
                   + f_4 * msd0_183[k]
                   - f_5 * msd1_183[k]
                   + f_3 * pc_x[k] * msf_303[k];

        t_454[k] = f_15 * lsf_232[k]
                   + f_3 * pc_y[k] * msf_302[k];
    }

#pragma omp simd aligned(t_455, t_456, t_457, t_458, pc_x, lsf_305, lsf_306, lsf_307, lsf_308, \
                         msd0_185, msd1_185, msf_305, msf_306, msf_307, \
                         msf_308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_455[k] = f_8 * lsf_305[k]
                   + f_4 * msd0_185[k]
                   - f_5 * msd1_185[k]
                   + f_3 * pc_x[k] * msf_305[k];

        t_456[k] = f_8 * lsf_306[k]
                   + f_3 * pc_x[k] * msf_306[k];

        t_457[k] = f_8 * lsf_307[k]
                   + f_3 * pc_x[k] * msf_307[k];

        t_458[k] = f_8 * lsf_308[k]
                   + f_3 * pc_x[k] * msf_308[k];
    }

#pragma omp simd aligned(t_459, t_460, t_461, pc_x, pc_y, pc_z, lsf_226, lsf_236, lsf_309, \
                         msd0_183, msd1_183, msf_306, msf_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_459[k] = f_8 * lsf_309[k]
                   + f_3 * pc_x[k] * msf_309[k];

        t_460[k] = f_15 * lsf_236[k]
                   + f_1 * msd0_183[k]
                   - f_2 * msd1_183[k]
                   + f_3 * pc_y[k] * msf_306[k];

        t_461[k] = f_8 * lsf_226[k]
                   + f_3 * pc_z[k] * msf_306[k];
    }

#pragma omp simd aligned(t_462, t_463, t_464, pc_y, pc_z, lsf_229, lsf_238, lsf_239, msd0_185, \
                         msd1_185, msf_308, msf_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_462[k] = f_15 * lsf_238[k]
                   + f_4 * msd0_185[k]
                   - f_5 * msd1_185[k]
                   + f_3 * pc_y[k] * msf_308[k];

        t_463[k] = f_15 * lsf_239[k]
                   + f_3 * pc_y[k] * msf_309[k];

        t_464[k] = f_8 * lsf_229[k]
                   + f_1 * msd0_185[k]
                   - f_2 * msd1_185[k]
                   + f_3 * pc_z[k] * msf_309[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, pc_x, pc_y, pc_z, lsf_230, lsf_240, lsf_310, \
                         msd0_186, msd1_186, msf_310 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = f_8 * lsf_310[k]
                   + f_1 * msd0_186[k]
                   - f_2 * msd1_186[k]
                   + f_3 * pc_x[k] * msf_310[k];

        t_466[k] = f_16 * lsf_240[k]
                   + f_3 * pc_y[k] * msf_310[k];

        t_467[k] = f_14 * lsf_230[k]
                   + f_3 * pc_z[k] * msf_310[k];
    }

#pragma omp simd aligned(t_468, t_469, t_470, pc_x, pc_y, lsf_242, lsf_313, lsf_315, msd0_189, \
                         msd0_191, msd1_189, msd1_191, msf_312, msf_313, \
                         msf_315 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_468[k] = f_8 * lsf_313[k]
                   + f_4 * msd0_189[k]
                   - f_5 * msd1_189[k]
                   + f_3 * pc_x[k] * msf_313[k];

        t_469[k] = f_16 * lsf_242[k]
                   + f_3 * pc_y[k] * msf_312[k];

        t_470[k] = f_8 * lsf_315[k]
                   + f_4 * msd0_191[k]
                   - f_5 * msd1_191[k]
                   + f_3 * pc_x[k] * msf_315[k];
    }

#pragma omp simd aligned(t_471, t_472, t_473, t_474, pc_x, lsf_316, lsf_317, lsf_318, lsf_319, \
                         msf_316, msf_317, msf_318, msf_319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_471[k] = f_8 * lsf_316[k]
                   + f_3 * pc_x[k] * msf_316[k];

        t_472[k] = f_8 * lsf_317[k]
                   + f_3 * pc_x[k] * msf_317[k];

        t_473[k] = f_8 * lsf_318[k]
                   + f_3 * pc_x[k] * msf_318[k];

        t_474[k] = f_8 * lsf_319[k]
                   + f_3 * pc_x[k] * msf_319[k];
    }

#pragma omp simd aligned(t_475, t_476, t_477, pc_y, pc_z, lsf_236, lsf_246, lsf_248, msd0_189, \
                         msd0_191, msd1_189, msd1_191, msf_316, \
                         msf_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_475[k] = f_16 * lsf_246[k]
                   + f_1 * msd0_189[k]
                   - f_2 * msd1_189[k]
                   + f_3 * pc_y[k] * msf_316[k];

        t_476[k] = f_14 * lsf_236[k]
                   + f_3 * pc_z[k] * msf_316[k];

        t_477[k] = f_16 * lsf_248[k]
                   + f_4 * msd0_191[k]
                   - f_5 * msd1_191[k]
                   + f_3 * pc_y[k] * msf_318[k];
    }

#pragma omp simd aligned(t_478, t_479, t_480, pc_x, pc_y, pc_z, lsf_239, lsf_249, lsf_320, \
                         msd0_191, msd0_192, msd1_191, msd1_192, msf_319, \
                         msf_320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_478[k] = f_16 * lsf_249[k]
                   + f_3 * pc_y[k] * msf_319[k];

        t_479[k] = f_14 * lsf_239[k]
                   + f_1 * msd0_191[k]
                   - f_2 * msd1_191[k]
                   + f_3 * pc_z[k] * msf_319[k];

        t_480[k] = f_8 * lsf_320[k]
                   + f_1 * msd0_192[k]
                   - f_2 * msd1_192[k]
                   + f_3 * pc_x[k] * msf_320[k];
    }

#pragma omp simd aligned(t_481, t_482, t_483, t_484, pc_x, pc_y, pc_z, lsf_240, lsf_250, \
                         lsf_252, lsf_323, msd0_195, msd1_195, msf_320, msf_322, \
                         msf_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_481[k] = f_14 * lsf_250[k]
                   + f_3 * pc_y[k] * msf_320[k];

        t_482[k] = f_16 * lsf_240[k]
                   + f_3 * pc_z[k] * msf_320[k];

        t_483[k] = f_8 * lsf_323[k]
                   + f_4 * msd0_195[k]
                   - f_5 * msd1_195[k]
                   + f_3 * pc_x[k] * msf_323[k];

        t_484[k] = f_14 * lsf_252[k]
                   + f_3 * pc_y[k] * msf_322[k];
    }
}

static auto
compute_prim_msg_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsg0,
                                                          const size_t lsf, const size_t lsg1,
                                                          const size_t msd0, const size_t msd1,
                                                          const size_t msf, const size_t ncols,
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
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 3.5 / q;
    const auto f_13 = 3.0 / q;
    const auto f_14 = 1.5 / q;
    const auto f_15 = 2.5 / q;
    const auto f_16 = 2.0 / q;

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
    auto *t_598 = buffer.data(target + 598);
    auto *t_599 = buffer.data(target + 599);
    auto *t_600 = buffer.data(target + 600);
    auto *t_601 = buffer.data(target + 601);
    auto *t_602 = buffer.data(target + 602);
    auto *t_603 = buffer.data(target + 603);
    auto *t_604 = buffer.data(target + 604);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *lsg0_405 = buffer.data(lsg0 + 405);
    const auto *lsg0_408 = buffer.data(lsg0 + 408);
    const auto *lsg0_410 = buffer.data(lsg0 + 410);
    const auto *lsg0_419 = buffer.data(lsg0 + 419);
    const auto *lsg0_420 = buffer.data(lsg0 + 420);
    const auto *lsg0_423 = buffer.data(lsg0 + 423);
    const auto *lsg0_540 = buffer.data(lsg0 + 540);
    const auto *lsg0_543 = buffer.data(lsg0 + 543);
    const auto *lsg0_550 = buffer.data(lsg0 + 550);
    const auto *lsg0_552 = buffer.data(lsg0 + 552);
    const auto *lsg0_554 = buffer.data(lsg0 + 554);
    const auto *lsg0_560 = buffer.data(lsg0 + 560);
    const auto *lsg0_565 = buffer.data(lsg0 + 565);
    const auto *lsg0_567 = buffer.data(lsg0 + 567);
    const auto *lsg0_569 = buffer.data(lsg0 + 569);
    const auto *lsg0_570 = buffer.data(lsg0 + 570);
    const auto *lsg0_573 = buffer.data(lsg0 + 573);
    const auto *lsg0_575 = buffer.data(lsg0 + 575);
    const auto *lsg0_580 = buffer.data(lsg0 + 580);
    const auto *lsg0_582 = buffer.data(lsg0 + 582);
    const auto *lsg0_584 = buffer.data(lsg0 + 584);
    const auto *lsg0_585 = buffer.data(lsg0 + 585);
    const auto *lsg0_588 = buffer.data(lsg0 + 588);
    const auto *lsg0_590 = buffer.data(lsg0 + 590);
    const auto *lsg0_595 = buffer.data(lsg0 + 595);
    const auto *lsg0_597 = buffer.data(lsg0 + 597);
    const auto *lsg0_599 = buffer.data(lsg0 + 599);
    const auto *lsg0_600 = buffer.data(lsg0 + 600);
    const auto *lsg0_603 = buffer.data(lsg0 + 603);

    const auto *lsf_246 = buffer.data(lsf + 246);
    const auto *lsf_249 = buffer.data(lsf + 249);
    const auto *lsf_250 = buffer.data(lsf + 250);
    const auto *lsf_256 = buffer.data(lsf + 256);
    const auto *lsf_258 = buffer.data(lsf + 258);
    const auto *lsf_259 = buffer.data(lsf + 259);
    const auto *lsf_260 = buffer.data(lsf + 260);
    const auto *lsf_262 = buffer.data(lsf + 262);
    const auto *lsf_266 = buffer.data(lsf + 266);
    const auto *lsf_268 = buffer.data(lsf + 268);
    const auto *lsf_269 = buffer.data(lsf + 269);
    const auto *lsf_270 = buffer.data(lsf + 270);
    const auto *lsf_271 = buffer.data(lsf + 271);
    const auto *lsf_272 = buffer.data(lsf + 272);
    const auto *lsf_276 = buffer.data(lsf + 276);
    const auto *lsf_278 = buffer.data(lsf + 278);
    const auto *lsf_279 = buffer.data(lsf + 279);
    const auto *lsf_280 = buffer.data(lsf + 280);
    const auto *lsf_286 = buffer.data(lsf + 286);
    const auto *lsf_289 = buffer.data(lsf + 289);
    const auto *lsf_290 = buffer.data(lsf + 290);
    const auto *lsf_292 = buffer.data(lsf + 292);
    const auto *lsf_296 = buffer.data(lsf + 296);
    const auto *lsf_299 = buffer.data(lsf + 299);
    const auto *lsf_300 = buffer.data(lsf + 300);
    const auto *lsf_302 = buffer.data(lsf + 302);
    const auto *lsf_306 = buffer.data(lsf + 306);
    const auto *lsf_309 = buffer.data(lsf + 309);
    const auto *lsf_310 = buffer.data(lsf + 310);
    const auto *lsf_312 = buffer.data(lsf + 312);
    const auto *lsf_319 = buffer.data(lsf + 319);
    const auto *lsf_320 = buffer.data(lsf + 320);
    const auto *lsf_322 = buffer.data(lsf + 322);
    const auto *lsf_325 = buffer.data(lsf + 325);
    const auto *lsf_326 = buffer.data(lsf + 326);
    const auto *lsf_327 = buffer.data(lsf + 327);
    const auto *lsf_328 = buffer.data(lsf + 328);
    const auto *lsf_329 = buffer.data(lsf + 329);
    const auto *lsf_330 = buffer.data(lsf + 330);
    const auto *lsf_333 = buffer.data(lsf + 333);
    const auto *lsf_335 = buffer.data(lsf + 335);
    const auto *lsf_336 = buffer.data(lsf + 336);
    const auto *lsf_337 = buffer.data(lsf + 337);
    const auto *lsf_338 = buffer.data(lsf + 338);
    const auto *lsf_339 = buffer.data(lsf + 339);
    const auto *lsf_346 = buffer.data(lsf + 346);
    const auto *lsf_347 = buffer.data(lsf + 347);
    const auto *lsf_348 = buffer.data(lsf + 348);
    const auto *lsf_349 = buffer.data(lsf + 349);
    const auto *lsf_350 = buffer.data(lsf + 350);
    const auto *lsf_355 = buffer.data(lsf + 355);
    const auto *lsf_356 = buffer.data(lsf + 356);
    const auto *lsf_357 = buffer.data(lsf + 357);
    const auto *lsf_359 = buffer.data(lsf + 359);
    const auto *lsf_360 = buffer.data(lsf + 360);
    const auto *lsf_363 = buffer.data(lsf + 363);
    const auto *lsf_366 = buffer.data(lsf + 366);
    const auto *lsf_368 = buffer.data(lsf + 368);
    const auto *lsf_369 = buffer.data(lsf + 369);
    const auto *lsf_375 = buffer.data(lsf + 375);
    const auto *lsf_376 = buffer.data(lsf + 376);
    const auto *lsf_377 = buffer.data(lsf + 377);
    const auto *lsf_378 = buffer.data(lsf + 378);
    const auto *lsf_379 = buffer.data(lsf + 379);
    const auto *lsf_380 = buffer.data(lsf + 380);
    const auto *lsf_383 = buffer.data(lsf + 383);
    const auto *lsf_385 = buffer.data(lsf + 385);
    const auto *lsf_386 = buffer.data(lsf + 386);
    const auto *lsf_387 = buffer.data(lsf + 387);
    const auto *lsf_388 = buffer.data(lsf + 388);
    const auto *lsf_389 = buffer.data(lsf + 389);
    const auto *lsf_390 = buffer.data(lsf + 390);
    const auto *lsf_393 = buffer.data(lsf + 393);
    const auto *lsf_395 = buffer.data(lsf + 395);
    const auto *lsf_396 = buffer.data(lsf + 396);
    const auto *lsf_397 = buffer.data(lsf + 397);
    const auto *lsf_398 = buffer.data(lsf + 398);
    const auto *lsf_399 = buffer.data(lsf + 399);
    const auto *lsf_400 = buffer.data(lsf + 400);
    const auto *lsf_403 = buffer.data(lsf + 403);

    const auto *lsg1_405 = buffer.data(lsg1 + 405);
    const auto *lsg1_408 = buffer.data(lsg1 + 408);
    const auto *lsg1_410 = buffer.data(lsg1 + 410);
    const auto *lsg1_419 = buffer.data(lsg1 + 419);
    const auto *lsg1_420 = buffer.data(lsg1 + 420);
    const auto *lsg1_423 = buffer.data(lsg1 + 423);
    const auto *lsg1_540 = buffer.data(lsg1 + 540);
    const auto *lsg1_543 = buffer.data(lsg1 + 543);
    const auto *lsg1_550 = buffer.data(lsg1 + 550);
    const auto *lsg1_552 = buffer.data(lsg1 + 552);
    const auto *lsg1_554 = buffer.data(lsg1 + 554);
    const auto *lsg1_560 = buffer.data(lsg1 + 560);
    const auto *lsg1_565 = buffer.data(lsg1 + 565);
    const auto *lsg1_567 = buffer.data(lsg1 + 567);
    const auto *lsg1_569 = buffer.data(lsg1 + 569);
    const auto *lsg1_570 = buffer.data(lsg1 + 570);
    const auto *lsg1_573 = buffer.data(lsg1 + 573);
    const auto *lsg1_575 = buffer.data(lsg1 + 575);
    const auto *lsg1_580 = buffer.data(lsg1 + 580);
    const auto *lsg1_582 = buffer.data(lsg1 + 582);
    const auto *lsg1_584 = buffer.data(lsg1 + 584);
    const auto *lsg1_585 = buffer.data(lsg1 + 585);
    const auto *lsg1_588 = buffer.data(lsg1 + 588);
    const auto *lsg1_590 = buffer.data(lsg1 + 590);
    const auto *lsg1_595 = buffer.data(lsg1 + 595);
    const auto *lsg1_597 = buffer.data(lsg1 + 597);
    const auto *lsg1_599 = buffer.data(lsg1 + 599);
    const auto *lsg1_600 = buffer.data(lsg1 + 600);
    const auto *lsg1_603 = buffer.data(lsg1 + 603);

    const auto *msd0_195 = buffer.data(msd0 + 195);
    const auto *msd0_197 = buffer.data(msd0 + 197);
    const auto *msd0_198 = buffer.data(msd0 + 198);
    const auto *msd0_201 = buffer.data(msd0 + 201);
    const auto *msd0_203 = buffer.data(msd0 + 203);
    const auto *msd0_207 = buffer.data(msd0 + 207);
    const auto *msd0_209 = buffer.data(msd0 + 209);
    const auto *msd0_210 = buffer.data(msd0 + 210);
    const auto *msd0_213 = buffer.data(msd0 + 213);
    const auto *msd0_214 = buffer.data(msd0 + 214);
    const auto *msd0_215 = buffer.data(msd0 + 215);
    const auto *msd0_216 = buffer.data(msd0 + 216);

    const auto *msd1_195 = buffer.data(msd1 + 195);
    const auto *msd1_197 = buffer.data(msd1 + 197);
    const auto *msd1_198 = buffer.data(msd1 + 198);
    const auto *msd1_201 = buffer.data(msd1 + 201);
    const auto *msd1_203 = buffer.data(msd1 + 203);
    const auto *msd1_207 = buffer.data(msd1 + 207);
    const auto *msd1_209 = buffer.data(msd1 + 209);
    const auto *msd1_210 = buffer.data(msd1 + 210);
    const auto *msd1_213 = buffer.data(msd1 + 213);
    const auto *msd1_214 = buffer.data(msd1 + 214);
    const auto *msd1_215 = buffer.data(msd1 + 215);
    const auto *msd1_216 = buffer.data(msd1 + 216);

    const auto *msf_325 = buffer.data(msf + 325);
    const auto *msf_326 = buffer.data(msf + 326);
    const auto *msf_327 = buffer.data(msf + 327);
    const auto *msf_328 = buffer.data(msf + 328);
    const auto *msf_329 = buffer.data(msf + 329);
    const auto *msf_330 = buffer.data(msf + 330);
    const auto *msf_332 = buffer.data(msf + 332);
    const auto *msf_333 = buffer.data(msf + 333);
    const auto *msf_335 = buffer.data(msf + 335);
    const auto *msf_336 = buffer.data(msf + 336);
    const auto *msf_337 = buffer.data(msf + 337);
    const auto *msf_338 = buffer.data(msf + 338);
    const auto *msf_339 = buffer.data(msf + 339);
    const auto *msf_340 = buffer.data(msf + 340);
    const auto *msf_342 = buffer.data(msf + 342);
    const auto *msf_346 = buffer.data(msf + 346);
    const auto *msf_347 = buffer.data(msf + 347);
    const auto *msf_348 = buffer.data(msf + 348);
    const auto *msf_349 = buffer.data(msf + 349);
    const auto *msf_350 = buffer.data(msf + 350);
    const auto *msf_351 = buffer.data(msf + 351);
    const auto *msf_352 = buffer.data(msf + 352);
    const auto *msf_355 = buffer.data(msf + 355);
    const auto *msf_356 = buffer.data(msf + 356);
    const auto *msf_357 = buffer.data(msf + 357);
    const auto *msf_358 = buffer.data(msf + 358);
    const auto *msf_359 = buffer.data(msf + 359);
    const auto *msf_360 = buffer.data(msf + 360);
    const auto *msf_361 = buffer.data(msf + 361);
    const auto *msf_362 = buffer.data(msf + 362);
    const auto *msf_363 = buffer.data(msf + 363);
    const auto *msf_366 = buffer.data(msf + 366);
    const auto *msf_368 = buffer.data(msf + 368);
    const auto *msf_369 = buffer.data(msf + 369);
    const auto *msf_370 = buffer.data(msf + 370);
    const auto *msf_372 = buffer.data(msf + 372);
    const auto *msf_376 = buffer.data(msf + 376);
    const auto *msf_377 = buffer.data(msf + 377);
    const auto *msf_378 = buffer.data(msf + 378);
    const auto *msf_379 = buffer.data(msf + 379);
    const auto *msf_380 = buffer.data(msf + 380);
    const auto *msf_382 = buffer.data(msf + 382);
    const auto *msf_386 = buffer.data(msf + 386);
    const auto *msf_387 = buffer.data(msf + 387);
    const auto *msf_388 = buffer.data(msf + 388);
    const auto *msf_389 = buffer.data(msf + 389);
    const auto *msf_390 = buffer.data(msf + 390);
    const auto *msf_392 = buffer.data(msf + 392);
    const auto *msf_396 = buffer.data(msf + 396);
    const auto *msf_397 = buffer.data(msf + 397);
    const auto *msf_398 = buffer.data(msf + 398);
    const auto *msf_399 = buffer.data(msf + 399);
    const auto *msf_400 = buffer.data(msf + 400);
    const auto *msf_402 = buffer.data(msf + 402);

#pragma omp simd aligned(t_485, t_486, t_487, t_488, pc_x, lsf_325, lsf_326, lsf_327, lsf_328, \
                         msd0_197, msd1_197, msf_325, msf_326, msf_327, \
                         msf_328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = f_8 * lsf_325[k]
                   + f_4 * msd0_197[k]
                   - f_5 * msd1_197[k]
                   + f_3 * pc_x[k] * msf_325[k];

        t_486[k] = f_8 * lsf_326[k]
                   + f_3 * pc_x[k] * msf_326[k];

        t_487[k] = f_8 * lsf_327[k]
                   + f_3 * pc_x[k] * msf_327[k];

        t_488[k] = f_8 * lsf_328[k]
                   + f_3 * pc_x[k] * msf_328[k];
    }

#pragma omp simd aligned(t_489, t_490, t_491, pc_x, pc_y, pc_z, lsf_246, lsf_256, lsf_329, \
                         msd0_195, msd1_195, msf_326, msf_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_489[k] = f_8 * lsf_329[k]
                   + f_3 * pc_x[k] * msf_329[k];

        t_490[k] = f_14 * lsf_256[k]
                   + f_1 * msd0_195[k]
                   - f_2 * msd1_195[k]
                   + f_3 * pc_y[k] * msf_326[k];

        t_491[k] = f_16 * lsf_246[k]
                   + f_3 * pc_z[k] * msf_326[k];
    }

#pragma omp simd aligned(t_492, t_493, t_494, pc_y, pc_z, lsf_249, lsf_258, lsf_259, msd0_197, \
                         msd1_197, msf_328, msf_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_492[k] = f_14 * lsf_258[k]
                   + f_4 * msd0_197[k]
                   - f_5 * msd1_197[k]
                   + f_3 * pc_y[k] * msf_328[k];

        t_493[k] = f_14 * lsf_259[k]
                   + f_3 * pc_y[k] * msf_329[k];

        t_494[k] = f_16 * lsf_249[k]
                   + f_1 * msd0_197[k]
                   - f_2 * msd1_197[k]
                   + f_3 * pc_z[k] * msf_329[k];
    }

#pragma omp simd aligned(t_495, t_496, t_497, pc_x, pc_y, pc_z, lsf_250, lsf_260, lsf_330, \
                         msd0_198, msd1_198, msf_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_495[k] = f_8 * lsf_330[k]
                   + f_1 * msd0_198[k]
                   - f_2 * msd1_198[k]
                   + f_3 * pc_x[k] * msf_330[k];

        t_496[k] = f_8 * lsf_260[k]
                   + f_3 * pc_y[k] * msf_330[k];

        t_497[k] = f_15 * lsf_250[k]
                   + f_3 * pc_z[k] * msf_330[k];
    }

#pragma omp simd aligned(t_498, t_499, t_500, pc_x, pc_y, lsf_262, lsf_333, lsf_335, msd0_201, \
                         msd0_203, msd1_201, msd1_203, msf_332, msf_333, \
                         msf_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_498[k] = f_8 * lsf_333[k]
                   + f_4 * msd0_201[k]
                   - f_5 * msd1_201[k]
                   + f_3 * pc_x[k] * msf_333[k];

        t_499[k] = f_8 * lsf_262[k]
                   + f_3 * pc_y[k] * msf_332[k];

        t_500[k] = f_8 * lsf_335[k]
                   + f_4 * msd0_203[k]
                   - f_5 * msd1_203[k]
                   + f_3 * pc_x[k] * msf_335[k];
    }

#pragma omp simd aligned(t_501, t_502, t_503, t_504, pc_x, lsf_336, lsf_337, lsf_338, lsf_339, \
                         msf_336, msf_337, msf_338, msf_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_501[k] = f_8 * lsf_336[k]
                   + f_3 * pc_x[k] * msf_336[k];

        t_502[k] = f_8 * lsf_337[k]
                   + f_3 * pc_x[k] * msf_337[k];

        t_503[k] = f_8 * lsf_338[k]
                   + f_3 * pc_x[k] * msf_338[k];

        t_504[k] = f_8 * lsf_339[k]
                   + f_3 * pc_x[k] * msf_339[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, pc_y, pc_z, lsf_256, lsf_266, lsf_268, msd0_201, \
                         msd0_203, msd1_201, msd1_203, msf_336, \
                         msf_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = f_8 * lsf_266[k]
                   + f_1 * msd0_201[k]
                   - f_2 * msd1_201[k]
                   + f_3 * pc_y[k] * msf_336[k];

        t_506[k] = f_15 * lsf_256[k]
                   + f_3 * pc_z[k] * msf_336[k];

        t_507[k] = f_8 * lsf_268[k]
                   + f_4 * msd0_203[k]
                   - f_5 * msd1_203[k]
                   + f_3 * pc_y[k] * msf_338[k];
    }

#pragma omp simd aligned(t_508, t_509, t_510, t_511, pa_y, pc_y, pc_z, lsg0_405, lsf_259, \
                         lsf_269, lsf_270, lsg1_405, msd0_203, msd1_203, msf_339, \
                         msf_340 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_508[k] = f_8 * lsf_269[k]
                   + f_3 * pc_y[k] * msf_339[k];

        t_509[k] = f_15 * lsf_259[k]
                   + f_1 * msd0_203[k]
                   - f_2 * msd1_203[k]
                   + f_3 * pc_z[k] * msf_339[k];

        t_510[k] = pa_y[k] * lsg0_405[k]
                   - f_6 * pc_y[k] * lsg1_405[k];

        t_511[k] = f_7 * lsf_270[k]
                   + f_3 * pc_y[k] * msf_340[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, t_515, pa_y, pc_y, pc_z, lsg0_408, lsg0_410, \
                         lsf_260, lsf_271, lsf_272, lsg1_408, lsg1_410, msf_340, \
                         msf_342 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = f_13 * lsf_260[k]
                   + f_3 * pc_z[k] * msf_340[k];

        t_513[k] = pa_y[k] * lsg0_408[k]
                   + f_8 * lsf_271[k]
                   - f_6 * pc_y[k] * lsg1_408[k];

        t_514[k] = f_7 * lsf_272[k]
                   + f_3 * pc_y[k] * msf_342[k];

        t_515[k] = pa_y[k] * lsg0_410[k]
                   - f_6 * pc_y[k] * lsg1_410[k];
    }

#pragma omp simd aligned(t_516, t_517, t_518, t_519, pc_x, lsf_346, lsf_347, lsf_348, lsf_349, \
                         msf_346, msf_347, msf_348, msf_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_516[k] = f_8 * lsf_346[k]
                   + f_3 * pc_x[k] * msf_346[k];

        t_517[k] = f_8 * lsf_347[k]
                   + f_3 * pc_x[k] * msf_347[k];

        t_518[k] = f_8 * lsf_348[k]
                   + f_3 * pc_x[k] * msf_348[k];

        t_519[k] = f_8 * lsf_349[k]
                   + f_3 * pc_x[k] * msf_349[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, pc_y, pc_z, lsf_266, lsf_276, lsf_278, msd0_207, \
                         msd0_209, msd1_207, msd1_209, msf_346, \
                         msf_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = f_7 * lsf_276[k]
                   + f_1 * msd0_207[k]
                   - f_2 * msd1_207[k]
                   + f_3 * pc_y[k] * msf_346[k];

        t_521[k] = f_13 * lsf_266[k]
                   + f_3 * pc_z[k] * msf_346[k];

        t_522[k] = f_7 * lsf_278[k]
                   + f_4 * msd0_209[k]
                   - f_5 * msd1_209[k]
                   + f_3 * pc_y[k] * msf_348[k];
    }

#pragma omp simd aligned(t_523, t_524, t_525, t_526, pa_y, pc_x, pc_y, lsg0_419, lsf_279, \
                         lsf_350, lsg1_419, msd0_210, msd1_210, msf_349, \
                         msf_350 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_523[k] = f_7 * lsf_279[k]
                   + f_3 * pc_y[k] * msf_349[k];

        t_524[k] = pa_y[k] * lsg0_419[k]
                   - f_6 * pc_y[k] * lsg1_419[k];

        t_525[k] = f_8 * lsf_350[k]
                   + f_1 * msd0_210[k]
                   - f_2 * msd1_210[k]
                   + f_3 * pc_x[k] * msf_350[k];

        t_526[k] = f_3 * pc_y[k] * msf_350[k];
    }

#pragma omp simd aligned(t_527, t_528, t_529, pc_y, pc_z, lsf_270, msd0_210, msd1_210, \
                         msf_350, msf_351, msf_352 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_527[k] = f_12 * lsf_270[k]
                   + f_3 * pc_z[k] * msf_350[k];

        t_528[k] = f_4 * msd0_210[k]
                   - f_5 * msd1_210[k]
                   + f_3 * pc_y[k] * msf_351[k];

        t_529[k] = f_3 * pc_y[k] * msf_352[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, t_533, pc_x, pc_y, lsf_355, lsf_356, lsf_357, \
                         msd0_215, msd1_215, msf_355, msf_356, \
                         msf_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = f_8 * lsf_355[k]
                   + f_4 * msd0_215[k]
                   - f_5 * msd1_215[k]
                   + f_3 * pc_x[k] * msf_355[k];

        t_531[k] = f_8 * lsf_356[k]
                   + f_3 * pc_x[k] * msf_356[k];

        t_532[k] = f_8 * lsf_357[k]
                   + f_3 * pc_x[k] * msf_357[k];

        t_533[k] = f_3 * pc_y[k] * msf_355[k];
    }

#pragma omp simd aligned(t_534, t_535, t_536, pc_x, pc_y, lsf_359, msd0_213, msd0_214, \
                         msd1_213, msd1_214, msf_356, msf_357, \
                         msf_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_534[k] = f_8 * lsf_359[k]
                   + f_3 * pc_x[k] * msf_359[k];

        t_535[k] = f_1 * msd0_213[k]
                   - f_2 * msd1_213[k]
                   + f_3 * pc_y[k] * msf_356[k];

        t_536[k] = f_10 * msd0_214[k]
                   - f_11 * msd1_214[k]
                   + f_3 * pc_y[k] * msf_357[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, t_540, pa_x, pc_x, pc_y, pc_z, lsg0_540, \
                         lsf_279, lsf_360, lsg1_540, msd0_215, msd1_215, msf_358, \
                         msf_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = f_4 * msd0_215[k]
                   - f_5 * msd1_215[k]
                   + f_3 * pc_y[k] * msf_358[k];

        t_538[k] = f_3 * pc_y[k] * msf_359[k];

        t_539[k] = f_12 * lsf_279[k]
                   + f_1 * msd0_215[k]
                   - f_2 * msd1_215[k]
                   + f_3 * pc_z[k] * msf_359[k];

        t_540[k] = pa_x[k] * lsg0_540[k]
                   + f_16 * lsf_360[k]
                   - f_6 * pc_x[k] * lsg1_540[k];
    }

#pragma omp simd aligned(t_541, t_542, t_543, t_544, pa_x, pc_x, pc_y, pc_z, lsg0_543, \
                         lsf_280, lsf_363, lsg1_543, msf_360, msf_361 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_541[k] = f_9 * lsf_280[k]
                   + f_3 * pc_y[k] * msf_360[k];

        t_542[k] = f_3 * pc_z[k] * msf_360[k];

        t_543[k] = pa_x[k] * lsg0_543[k]
                   + f_8 * lsf_363[k]
                   - f_6 * pc_x[k] * lsg1_543[k];

        t_544[k] = f_3 * pc_z[k] * msf_361[k];
    }

#pragma omp simd aligned(t_545, t_546, t_547, t_548, pc_x, pc_z, lsf_366, lsf_368, msd0_216, \
                         msd1_216, msf_362, msf_363, msf_366, msf_368 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_545[k] = f_4 * msd0_216[k]
                   - f_5 * msd1_216[k]
                   + f_3 * pc_z[k] * msf_362[k];

        t_546[k] = f_7 * lsf_366[k]
                   + f_3 * pc_x[k] * msf_366[k];

        t_547[k] = f_3 * pc_z[k] * msf_363[k];

        t_548[k] = f_7 * lsf_368[k]
                   + f_3 * pc_x[k] * msf_368[k];
    }

#pragma omp simd aligned(t_549, t_550, t_551, t_552, pa_x, pc_x, pc_z, lsg0_550, lsg0_552, \
                         lsf_369, lsg1_550, lsg1_552, msf_366, \
                         msf_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_549[k] = f_7 * lsf_369[k]
                   + f_3 * pc_x[k] * msf_369[k];

        t_550[k] = pa_x[k] * lsg0_550[k]
                   - f_6 * pc_x[k] * lsg1_550[k];

        t_551[k] = f_3 * pc_z[k] * msf_366[k];

        t_552[k] = pa_x[k] * lsg0_552[k]
                   - f_6 * pc_x[k] * lsg1_552[k];
    }

#pragma omp simd aligned(t_553, t_554, t_555, pa_x, pa_z, pc_x, pc_y, pc_z, lsg0_420, \
                         lsg0_554, lsf_289, lsg1_420, lsg1_554, \
                         msf_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_553[k] = f_9 * lsf_289[k]
                   + f_3 * pc_y[k] * msf_369[k];

        t_554[k] = pa_x[k] * lsg0_554[k]
                   - f_6 * pc_x[k] * lsg1_554[k];

        t_555[k] = pa_z[k] * lsg0_420[k]
                   - f_6 * pc_z[k] * lsg1_420[k];
    }

#pragma omp simd aligned(t_556, t_557, t_558, t_559, pa_z, pc_y, pc_z, lsg0_423, lsf_280, \
                         lsf_290, lsf_292, lsg1_423, msf_370, msf_372 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_556[k] = f_12 * lsf_290[k]
                   + f_3 * pc_y[k] * msf_370[k];

        t_557[k] = f_7 * lsf_280[k]
                   + f_3 * pc_z[k] * msf_370[k];

        t_558[k] = pa_z[k] * lsg0_423[k]
                   - f_6 * pc_z[k] * lsg1_423[k];

        t_559[k] = f_12 * lsf_292[k]
                   + f_3 * pc_y[k] * msf_372[k];
    }

#pragma omp simd aligned(t_560, t_561, t_562, t_563, pa_x, pc_x, lsg0_560, lsf_375, lsf_376, \
                         lsf_377, lsf_378, lsg1_560, msf_376, msf_377, \
                         msf_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_560[k] = pa_x[k] * lsg0_560[k]
                   + f_8 * lsf_375[k]
                   - f_6 * pc_x[k] * lsg1_560[k];

        t_561[k] = f_7 * lsf_376[k]
                   + f_3 * pc_x[k] * msf_376[k];

        t_562[k] = f_7 * lsf_377[k]
                   + f_3 * pc_x[k] * msf_377[k];

        t_563[k] = f_7 * lsf_378[k]
                   + f_3 * pc_x[k] * msf_378[k];
    }

#pragma omp simd aligned(t_564, t_565, t_566, t_567, pa_x, pc_x, pc_z, lsg0_565, lsg0_567, \
                         lsf_286, lsf_379, lsg1_565, lsg1_567, msf_376, \
                         msf_379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_564[k] = f_7 * lsf_379[k]
                   + f_3 * pc_x[k] * msf_379[k];

        t_565[k] = pa_x[k] * lsg0_565[k]
                   - f_6 * pc_x[k] * lsg1_565[k];

        t_566[k] = f_7 * lsf_286[k]
                   + f_3 * pc_z[k] * msf_376[k];

        t_567[k] = pa_x[k] * lsg0_567[k]
                   - f_6 * pc_x[k] * lsg1_567[k];
    }

#pragma omp simd aligned(t_568, t_569, t_570, t_571, pa_x, pc_x, pc_y, lsg0_569, lsg0_570, \
                         lsf_299, lsf_300, lsf_380, lsg1_569, lsg1_570, msf_379, \
                         msf_380 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_568[k] = f_12 * lsf_299[k]
                   + f_3 * pc_y[k] * msf_379[k];

        t_569[k] = pa_x[k] * lsg0_569[k]
                   - f_6 * pc_x[k] * lsg1_569[k];

        t_570[k] = pa_x[k] * lsg0_570[k]
                   + f_16 * lsf_380[k]
                   - f_6 * pc_x[k] * lsg1_570[k];

        t_571[k] = f_13 * lsf_300[k]
                   + f_3 * pc_y[k] * msf_380[k];
    }

#pragma omp simd aligned(t_572, t_573, t_574, pa_x, pc_x, pc_y, pc_z, lsg0_573, lsf_290, \
                         lsf_302, lsf_383, lsg1_573, msf_380, msf_382 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_572[k] = f_8 * lsf_290[k]
                   + f_3 * pc_z[k] * msf_380[k];

        t_573[k] = pa_x[k] * lsg0_573[k]
                   + f_8 * lsf_383[k]
                   - f_6 * pc_x[k] * lsg1_573[k];

        t_574[k] = f_13 * lsf_302[k]
                   + f_3 * pc_y[k] * msf_382[k];
    }

#pragma omp simd aligned(t_575, t_576, t_577, t_578, pa_x, pc_x, lsg0_575, lsf_385, lsf_386, \
                         lsf_387, lsf_388, lsg1_575, msf_386, msf_387, \
                         msf_388 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_575[k] = pa_x[k] * lsg0_575[k]
                   + f_8 * lsf_385[k]
                   - f_6 * pc_x[k] * lsg1_575[k];

        t_576[k] = f_7 * lsf_386[k]
                   + f_3 * pc_x[k] * msf_386[k];

        t_577[k] = f_7 * lsf_387[k]
                   + f_3 * pc_x[k] * msf_387[k];

        t_578[k] = f_7 * lsf_388[k]
                   + f_3 * pc_x[k] * msf_388[k];
    }

#pragma omp simd aligned(t_579, t_580, t_581, t_582, pa_x, pc_x, pc_z, lsg0_580, lsg0_582, \
                         lsf_296, lsf_389, lsg1_580, lsg1_582, msf_386, \
                         msf_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_579[k] = f_7 * lsf_389[k]
                   + f_3 * pc_x[k] * msf_389[k];

        t_580[k] = pa_x[k] * lsg0_580[k]
                   - f_6 * pc_x[k] * lsg1_580[k];

        t_581[k] = f_8 * lsf_296[k]
                   + f_3 * pc_z[k] * msf_386[k];

        t_582[k] = pa_x[k] * lsg0_582[k]
                   - f_6 * pc_x[k] * lsg1_582[k];
    }

#pragma omp simd aligned(t_583, t_584, t_585, t_586, pa_x, pc_x, pc_y, lsg0_584, lsg0_585, \
                         lsf_309, lsf_310, lsf_390, lsg1_584, lsg1_585, msf_389, \
                         msf_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_583[k] = f_13 * lsf_309[k]
                   + f_3 * pc_y[k] * msf_389[k];

        t_584[k] = pa_x[k] * lsg0_584[k]
                   - f_6 * pc_x[k] * lsg1_584[k];

        t_585[k] = pa_x[k] * lsg0_585[k]
                   + f_16 * lsf_390[k]
                   - f_6 * pc_x[k] * lsg1_585[k];

        t_586[k] = f_15 * lsf_310[k]
                   + f_3 * pc_y[k] * msf_390[k];
    }

#pragma omp simd aligned(t_587, t_588, t_589, pa_x, pc_x, pc_y, pc_z, lsg0_588, lsf_300, \
                         lsf_312, lsf_393, lsg1_588, msf_390, msf_392 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_587[k] = f_14 * lsf_300[k]
                   + f_3 * pc_z[k] * msf_390[k];

        t_588[k] = pa_x[k] * lsg0_588[k]
                   + f_8 * lsf_393[k]
                   - f_6 * pc_x[k] * lsg1_588[k];

        t_589[k] = f_15 * lsf_312[k]
                   + f_3 * pc_y[k] * msf_392[k];
    }

#pragma omp simd aligned(t_590, t_591, t_592, t_593, pa_x, pc_x, lsg0_590, lsf_395, lsf_396, \
                         lsf_397, lsf_398, lsg1_590, msf_396, msf_397, \
                         msf_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_590[k] = pa_x[k] * lsg0_590[k]
                   + f_8 * lsf_395[k]
                   - f_6 * pc_x[k] * lsg1_590[k];

        t_591[k] = f_7 * lsf_396[k]
                   + f_3 * pc_x[k] * msf_396[k];

        t_592[k] = f_7 * lsf_397[k]
                   + f_3 * pc_x[k] * msf_397[k];

        t_593[k] = f_7 * lsf_398[k]
                   + f_3 * pc_x[k] * msf_398[k];
    }

#pragma omp simd aligned(t_594, t_595, t_596, t_597, pa_x, pc_x, pc_z, lsg0_595, lsg0_597, \
                         lsf_306, lsf_399, lsg1_595, lsg1_597, msf_396, \
                         msf_399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_594[k] = f_7 * lsf_399[k]
                   + f_3 * pc_x[k] * msf_399[k];

        t_595[k] = pa_x[k] * lsg0_595[k]
                   - f_6 * pc_x[k] * lsg1_595[k];

        t_596[k] = f_14 * lsf_306[k]
                   + f_3 * pc_z[k] * msf_396[k];

        t_597[k] = pa_x[k] * lsg0_597[k]
                   - f_6 * pc_x[k] * lsg1_597[k];
    }

#pragma omp simd aligned(t_598, t_599, t_600, t_601, pa_x, pc_x, pc_y, lsg0_599, lsg0_600, \
                         lsf_319, lsf_320, lsf_400, lsg1_599, lsg1_600, msf_399, \
                         msf_400 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_598[k] = f_15 * lsf_319[k]
                   + f_3 * pc_y[k] * msf_399[k];

        t_599[k] = pa_x[k] * lsg0_599[k]
                   - f_6 * pc_x[k] * lsg1_599[k];

        t_600[k] = pa_x[k] * lsg0_600[k]
                   + f_16 * lsf_400[k]
                   - f_6 * pc_x[k] * lsg1_600[k];

        t_601[k] = f_16 * lsf_320[k]
                   + f_3 * pc_y[k] * msf_400[k];
    }

#pragma omp simd aligned(t_602, t_603, t_604, pa_x, pc_x, pc_y, pc_z, lsg0_603, lsf_310, \
                         lsf_322, lsf_403, lsg1_603, msf_400, msf_402 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_602[k] = f_16 * lsf_310[k]
                   + f_3 * pc_z[k] * msf_400[k];

        t_603[k] = pa_x[k] * lsg0_603[k]
                   + f_8 * lsf_403[k]
                   - f_6 * pc_x[k] * lsg1_603[k];

        t_604[k] = f_16 * lsf_322[k]
                   + f_3 * pc_y[k] * msf_402[k];
    }
}

static auto
compute_prim_msg_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsg0,
                                                          const size_t lsf, const size_t lsg1,
                                                          const size_t msd0, const size_t msd1,
                                                          const size_t msf, const size_t ncols,
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
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 3.5 / q;
    const auto f_13 = 3.0 / q;
    const auto f_14 = 1.5 / q;
    const auto f_15 = 2.5 / q;
    const auto f_16 = 2.0 / q;

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
    auto *t_729 = buffer.data(target + 729);
    auto *t_730 = buffer.data(target + 730);
    auto *t_731 = buffer.data(target + 731);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *lsg0_525 = buffer.data(lsg0 + 525);
    const auto *lsg0_530 = buffer.data(lsg0 + 530);
    const auto *lsg0_540 = buffer.data(lsg0 + 540);
    const auto *lsg0_541 = buffer.data(lsg0 + 541);
    const auto *lsg0_543 = buffer.data(lsg0 + 543);
    const auto *lsg0_550 = buffer.data(lsg0 + 550);
    const auto *lsg0_552 = buffer.data(lsg0 + 552);
    const auto *lsg0_605 = buffer.data(lsg0 + 605);
    const auto *lsg0_610 = buffer.data(lsg0 + 610);
    const auto *lsg0_612 = buffer.data(lsg0 + 612);
    const auto *lsg0_614 = buffer.data(lsg0 + 614);
    const auto *lsg0_615 = buffer.data(lsg0 + 615);
    const auto *lsg0_618 = buffer.data(lsg0 + 618);
    const auto *lsg0_620 = buffer.data(lsg0 + 620);
    const auto *lsg0_625 = buffer.data(lsg0 + 625);
    const auto *lsg0_627 = buffer.data(lsg0 + 627);
    const auto *lsg0_629 = buffer.data(lsg0 + 629);
    const auto *lsg0_630 = buffer.data(lsg0 + 630);
    const auto *lsg0_633 = buffer.data(lsg0 + 633);
    const auto *lsg0_635 = buffer.data(lsg0 + 635);
    const auto *lsg0_640 = buffer.data(lsg0 + 640);
    const auto *lsg0_642 = buffer.data(lsg0 + 642);
    const auto *lsg0_644 = buffer.data(lsg0 + 644);
    const auto *lsg0_648 = buffer.data(lsg0 + 648);
    const auto *lsg0_655 = buffer.data(lsg0 + 655);
    const auto *lsg0_657 = buffer.data(lsg0 + 657);
    const auto *lsg0_659 = buffer.data(lsg0 + 659);
    const auto *lsg0_660 = buffer.data(lsg0 + 660);
    const auto *lsg0_665 = buffer.data(lsg0 + 665);
    const auto *lsg0_670 = buffer.data(lsg0 + 670);
    const auto *lsg0_671 = buffer.data(lsg0 + 671);
    const auto *lsg0_672 = buffer.data(lsg0 + 672);
    const auto *lsg0_674 = buffer.data(lsg0 + 674);

    const auto *lsf_316 = buffer.data(lsf + 316);
    const auto *lsf_320 = buffer.data(lsf + 320);
    const auto *lsf_326 = buffer.data(lsf + 326);
    const auto *lsf_329 = buffer.data(lsf + 329);
    const auto *lsf_330 = buffer.data(lsf + 330);
    const auto *lsf_332 = buffer.data(lsf + 332);
    const auto *lsf_336 = buffer.data(lsf + 336);
    const auto *lsf_339 = buffer.data(lsf + 339);
    const auto *lsf_340 = buffer.data(lsf + 340);
    const auto *lsf_342 = buffer.data(lsf + 342);
    const auto *lsf_346 = buffer.data(lsf + 346);
    const auto *lsf_349 = buffer.data(lsf + 349);
    const auto *lsf_350 = buffer.data(lsf + 350);
    const auto *lsf_352 = buffer.data(lsf + 352);
    const auto *lsf_359 = buffer.data(lsf + 359);
    const auto *lsf_366 = buffer.data(lsf + 366);
    const auto *lsf_367 = buffer.data(lsf + 367);
    const auto *lsf_369 = buffer.data(lsf + 369);
    const auto *lsf_376 = buffer.data(lsf + 376);
    const auto *lsf_379 = buffer.data(lsf + 379);
    const auto *lsf_386 = buffer.data(lsf + 386);
    const auto *lsf_388 = buffer.data(lsf + 388);
    const auto *lsf_389 = buffer.data(lsf + 389);
    const auto *lsf_396 = buffer.data(lsf + 396);
    const auto *lsf_405 = buffer.data(lsf + 405);
    const auto *lsf_406 = buffer.data(lsf + 406);
    const auto *lsf_407 = buffer.data(lsf + 407);
    const auto *lsf_408 = buffer.data(lsf + 408);
    const auto *lsf_409 = buffer.data(lsf + 409);
    const auto *lsf_410 = buffer.data(lsf + 410);
    const auto *lsf_413 = buffer.data(lsf + 413);
    const auto *lsf_415 = buffer.data(lsf + 415);
    const auto *lsf_416 = buffer.data(lsf + 416);
    const auto *lsf_417 = buffer.data(lsf + 417);
    const auto *lsf_418 = buffer.data(lsf + 418);
    const auto *lsf_419 = buffer.data(lsf + 419);
    const auto *lsf_420 = buffer.data(lsf + 420);
    const auto *lsf_423 = buffer.data(lsf + 423);
    const auto *lsf_425 = buffer.data(lsf + 425);
    const auto *lsf_426 = buffer.data(lsf + 426);
    const auto *lsf_427 = buffer.data(lsf + 427);
    const auto *lsf_428 = buffer.data(lsf + 428);
    const auto *lsf_429 = buffer.data(lsf + 429);
    const auto *lsf_433 = buffer.data(lsf + 433);
    const auto *lsf_436 = buffer.data(lsf + 436);
    const auto *lsf_437 = buffer.data(lsf + 437);
    const auto *lsf_438 = buffer.data(lsf + 438);
    const auto *lsf_439 = buffer.data(lsf + 439);
    const auto *lsf_440 = buffer.data(lsf + 440);
    const auto *lsf_445 = buffer.data(lsf + 445);
    const auto *lsf_446 = buffer.data(lsf + 446);
    const auto *lsf_447 = buffer.data(lsf + 447);
    const auto *lsf_449 = buffer.data(lsf + 449);

    const auto *lsg1_525 = buffer.data(lsg1 + 525);
    const auto *lsg1_530 = buffer.data(lsg1 + 530);
    const auto *lsg1_540 = buffer.data(lsg1 + 540);
    const auto *lsg1_541 = buffer.data(lsg1 + 541);
    const auto *lsg1_543 = buffer.data(lsg1 + 543);
    const auto *lsg1_550 = buffer.data(lsg1 + 550);
    const auto *lsg1_552 = buffer.data(lsg1 + 552);
    const auto *lsg1_605 = buffer.data(lsg1 + 605);
    const auto *lsg1_610 = buffer.data(lsg1 + 610);
    const auto *lsg1_612 = buffer.data(lsg1 + 612);
    const auto *lsg1_614 = buffer.data(lsg1 + 614);
    const auto *lsg1_615 = buffer.data(lsg1 + 615);
    const auto *lsg1_618 = buffer.data(lsg1 + 618);
    const auto *lsg1_620 = buffer.data(lsg1 + 620);
    const auto *lsg1_625 = buffer.data(lsg1 + 625);
    const auto *lsg1_627 = buffer.data(lsg1 + 627);
    const auto *lsg1_629 = buffer.data(lsg1 + 629);
    const auto *lsg1_630 = buffer.data(lsg1 + 630);
    const auto *lsg1_633 = buffer.data(lsg1 + 633);
    const auto *lsg1_635 = buffer.data(lsg1 + 635);
    const auto *lsg1_640 = buffer.data(lsg1 + 640);
    const auto *lsg1_642 = buffer.data(lsg1 + 642);
    const auto *lsg1_644 = buffer.data(lsg1 + 644);
    const auto *lsg1_648 = buffer.data(lsg1 + 648);
    const auto *lsg1_655 = buffer.data(lsg1 + 655);
    const auto *lsg1_657 = buffer.data(lsg1 + 657);
    const auto *lsg1_659 = buffer.data(lsg1 + 659);
    const auto *lsg1_660 = buffer.data(lsg1 + 660);
    const auto *lsg1_665 = buffer.data(lsg1 + 665);
    const auto *lsg1_670 = buffer.data(lsg1 + 670);
    const auto *lsg1_671 = buffer.data(lsg1 + 671);
    const auto *lsg1_672 = buffer.data(lsg1 + 672);
    const auto *lsg1_674 = buffer.data(lsg1 + 674);

    const auto *msd0_264 = buffer.data(msd0 + 264);
    const auto *msd0_270 = buffer.data(msd0 + 270);
    const auto *msd0_271 = buffer.data(msd0 + 271);
    const auto *msd0_273 = buffer.data(msd0 + 273);
    const auto *msd0_275 = buffer.data(msd0 + 275);
    const auto *msd0_278 = buffer.data(msd0 + 278);
    const auto *msd0_280 = buffer.data(msd0 + 280);
    const auto *msd0_281 = buffer.data(msd0 + 281);
    const auto *msd0_282 = buffer.data(msd0 + 282);
    const auto *msd0_283 = buffer.data(msd0 + 283);
    const auto *msd0_284 = buffer.data(msd0 + 284);
    const auto *msd0_285 = buffer.data(msd0 + 285);
    const auto *msd0_286 = buffer.data(msd0 + 286);
    const auto *msd0_287 = buffer.data(msd0 + 287);
    const auto *msd0_288 = buffer.data(msd0 + 288);
    const auto *msd0_289 = buffer.data(msd0 + 289);
    const auto *msd0_290 = buffer.data(msd0 + 290);
    const auto *msd0_291 = buffer.data(msd0 + 291);
    const auto *msd0_292 = buffer.data(msd0 + 292);
    const auto *msd0_293 = buffer.data(msd0 + 293);

    const auto *msd1_264 = buffer.data(msd1 + 264);
    const auto *msd1_270 = buffer.data(msd1 + 270);
    const auto *msd1_271 = buffer.data(msd1 + 271);
    const auto *msd1_273 = buffer.data(msd1 + 273);
    const auto *msd1_275 = buffer.data(msd1 + 275);
    const auto *msd1_278 = buffer.data(msd1 + 278);
    const auto *msd1_280 = buffer.data(msd1 + 280);
    const auto *msd1_281 = buffer.data(msd1 + 281);
    const auto *msd1_282 = buffer.data(msd1 + 282);
    const auto *msd1_283 = buffer.data(msd1 + 283);
    const auto *msd1_284 = buffer.data(msd1 + 284);
    const auto *msd1_285 = buffer.data(msd1 + 285);
    const auto *msd1_286 = buffer.data(msd1 + 286);
    const auto *msd1_287 = buffer.data(msd1 + 287);
    const auto *msd1_288 = buffer.data(msd1 + 288);
    const auto *msd1_289 = buffer.data(msd1 + 289);
    const auto *msd1_290 = buffer.data(msd1 + 290);
    const auto *msd1_291 = buffer.data(msd1 + 291);
    const auto *msd1_292 = buffer.data(msd1 + 292);
    const auto *msd1_293 = buffer.data(msd1 + 293);

    const auto *msf_406 = buffer.data(msf + 406);
    const auto *msf_407 = buffer.data(msf + 407);
    const auto *msf_408 = buffer.data(msf + 408);
    const auto *msf_409 = buffer.data(msf + 409);
    const auto *msf_410 = buffer.data(msf + 410);
    const auto *msf_412 = buffer.data(msf + 412);
    const auto *msf_416 = buffer.data(msf + 416);
    const auto *msf_417 = buffer.data(msf + 417);
    const auto *msf_418 = buffer.data(msf + 418);
    const auto *msf_419 = buffer.data(msf + 419);
    const auto *msf_420 = buffer.data(msf + 420);
    const auto *msf_422 = buffer.data(msf + 422);
    const auto *msf_426 = buffer.data(msf + 426);
    const auto *msf_427 = buffer.data(msf + 427);
    const auto *msf_428 = buffer.data(msf + 428);
    const auto *msf_429 = buffer.data(msf + 429);
    const auto *msf_430 = buffer.data(msf + 430);
    const auto *msf_432 = buffer.data(msf + 432);
    const auto *msf_436 = buffer.data(msf + 436);
    const auto *msf_437 = buffer.data(msf + 437);
    const auto *msf_438 = buffer.data(msf + 438);
    const auto *msf_439 = buffer.data(msf + 439);
    const auto *msf_440 = buffer.data(msf + 440);
    const auto *msf_441 = buffer.data(msf + 441);
    const auto *msf_442 = buffer.data(msf + 442);
    const auto *msf_445 = buffer.data(msf + 445);
    const auto *msf_446 = buffer.data(msf + 446);
    const auto *msf_447 = buffer.data(msf + 447);
    const auto *msf_449 = buffer.data(msf + 449);
    const auto *msf_450 = buffer.data(msf + 450);
    const auto *msf_451 = buffer.data(msf + 451);
    const auto *msf_453 = buffer.data(msf + 453);
    const auto *msf_455 = buffer.data(msf + 455);
    const auto *msf_456 = buffer.data(msf + 456);
    const auto *msf_457 = buffer.data(msf + 457);
    const auto *msf_458 = buffer.data(msf + 458);
    const auto *msf_459 = buffer.data(msf + 459);
    const auto *msf_462 = buffer.data(msf + 462);
    const auto *msf_464 = buffer.data(msf + 464);
    const auto *msf_465 = buffer.data(msf + 465);
    const auto *msf_466 = buffer.data(msf + 466);
    const auto *msf_467 = buffer.data(msf + 467);
    const auto *msf_468 = buffer.data(msf + 468);
    const auto *msf_469 = buffer.data(msf + 469);
    const auto *msf_470 = buffer.data(msf + 470);
    const auto *msf_471 = buffer.data(msf + 471);
    const auto *msf_472 = buffer.data(msf + 472);
    const auto *msf_473 = buffer.data(msf + 473);
    const auto *msf_474 = buffer.data(msf + 474);
    const auto *msf_475 = buffer.data(msf + 475);
    const auto *msf_476 = buffer.data(msf + 476);
    const auto *msf_477 = buffer.data(msf + 477);
    const auto *msf_478 = buffer.data(msf + 478);
    const auto *msf_479 = buffer.data(msf + 479);
    const auto *msf_480 = buffer.data(msf + 480);
    const auto *msf_481 = buffer.data(msf + 481);
    const auto *msf_482 = buffer.data(msf + 482);
    const auto *msf_483 = buffer.data(msf + 483);
    const auto *msf_484 = buffer.data(msf + 484);
    const auto *msf_485 = buffer.data(msf + 485);
    const auto *msf_486 = buffer.data(msf + 486);
    const auto *msf_487 = buffer.data(msf + 487);
    const auto *msf_488 = buffer.data(msf + 488);
    const auto *msf_489 = buffer.data(msf + 489);

#pragma omp simd aligned(t_605, t_606, t_607, t_608, pa_x, pc_x, lsg0_605, lsf_405, lsf_406, \
                         lsf_407, lsf_408, lsg1_605, msf_406, msf_407, \
                         msf_408 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_605[k] = pa_x[k] * lsg0_605[k]
                   + f_8 * lsf_405[k]
                   - f_6 * pc_x[k] * lsg1_605[k];

        t_606[k] = f_7 * lsf_406[k]
                   + f_3 * pc_x[k] * msf_406[k];

        t_607[k] = f_7 * lsf_407[k]
                   + f_3 * pc_x[k] * msf_407[k];

        t_608[k] = f_7 * lsf_408[k]
                   + f_3 * pc_x[k] * msf_408[k];
    }

#pragma omp simd aligned(t_609, t_610, t_611, t_612, pa_x, pc_x, pc_z, lsg0_610, lsg0_612, \
                         lsf_316, lsf_409, lsg1_610, lsg1_612, msf_406, \
                         msf_409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_609[k] = f_7 * lsf_409[k]
                   + f_3 * pc_x[k] * msf_409[k];

        t_610[k] = pa_x[k] * lsg0_610[k]
                   - f_6 * pc_x[k] * lsg1_610[k];

        t_611[k] = f_16 * lsf_316[k]
                   + f_3 * pc_z[k] * msf_406[k];

        t_612[k] = pa_x[k] * lsg0_612[k]
                   - f_6 * pc_x[k] * lsg1_612[k];
    }

#pragma omp simd aligned(t_613, t_614, t_615, t_616, pa_x, pc_x, pc_y, lsg0_614, lsg0_615, \
                         lsf_329, lsf_330, lsf_410, lsg1_614, lsg1_615, msf_409, \
                         msf_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_613[k] = f_16 * lsf_329[k]
                   + f_3 * pc_y[k] * msf_409[k];

        t_614[k] = pa_x[k] * lsg0_614[k]
                   - f_6 * pc_x[k] * lsg1_614[k];

        t_615[k] = pa_x[k] * lsg0_615[k]
                   + f_16 * lsf_410[k]
                   - f_6 * pc_x[k] * lsg1_615[k];

        t_616[k] = f_14 * lsf_330[k]
                   + f_3 * pc_y[k] * msf_410[k];
    }

#pragma omp simd aligned(t_617, t_618, t_619, pa_x, pc_x, pc_y, pc_z, lsg0_618, lsf_320, \
                         lsf_332, lsf_413, lsg1_618, msf_410, msf_412 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_617[k] = f_15 * lsf_320[k]
                   + f_3 * pc_z[k] * msf_410[k];

        t_618[k] = pa_x[k] * lsg0_618[k]
                   + f_8 * lsf_413[k]
                   - f_6 * pc_x[k] * lsg1_618[k];

        t_619[k] = f_14 * lsf_332[k]
                   + f_3 * pc_y[k] * msf_412[k];
    }

#pragma omp simd aligned(t_620, t_621, t_622, t_623, pa_x, pc_x, lsg0_620, lsf_415, lsf_416, \
                         lsf_417, lsf_418, lsg1_620, msf_416, msf_417, \
                         msf_418 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_620[k] = pa_x[k] * lsg0_620[k]
                   + f_8 * lsf_415[k]
                   - f_6 * pc_x[k] * lsg1_620[k];

        t_621[k] = f_7 * lsf_416[k]
                   + f_3 * pc_x[k] * msf_416[k];

        t_622[k] = f_7 * lsf_417[k]
                   + f_3 * pc_x[k] * msf_417[k];

        t_623[k] = f_7 * lsf_418[k]
                   + f_3 * pc_x[k] * msf_418[k];
    }

#pragma omp simd aligned(t_624, t_625, t_626, t_627, pa_x, pc_x, pc_z, lsg0_625, lsg0_627, \
                         lsf_326, lsf_419, lsg1_625, lsg1_627, msf_416, \
                         msf_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_624[k] = f_7 * lsf_419[k]
                   + f_3 * pc_x[k] * msf_419[k];

        t_625[k] = pa_x[k] * lsg0_625[k]
                   - f_6 * pc_x[k] * lsg1_625[k];

        t_626[k] = f_15 * lsf_326[k]
                   + f_3 * pc_z[k] * msf_416[k];

        t_627[k] = pa_x[k] * lsg0_627[k]
                   - f_6 * pc_x[k] * lsg1_627[k];
    }

#pragma omp simd aligned(t_628, t_629, t_630, t_631, pa_x, pc_x, pc_y, lsg0_629, lsg0_630, \
                         lsf_339, lsf_340, lsf_420, lsg1_629, lsg1_630, msf_419, \
                         msf_420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_628[k] = f_14 * lsf_339[k]
                   + f_3 * pc_y[k] * msf_419[k];

        t_629[k] = pa_x[k] * lsg0_629[k]
                   - f_6 * pc_x[k] * lsg1_629[k];

        t_630[k] = pa_x[k] * lsg0_630[k]
                   + f_16 * lsf_420[k]
                   - f_6 * pc_x[k] * lsg1_630[k];

        t_631[k] = f_8 * lsf_340[k]
                   + f_3 * pc_y[k] * msf_420[k];
    }

#pragma omp simd aligned(t_632, t_633, t_634, pa_x, pc_x, pc_y, pc_z, lsg0_633, lsf_330, \
                         lsf_342, lsf_423, lsg1_633, msf_420, msf_422 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_632[k] = f_13 * lsf_330[k]
                   + f_3 * pc_z[k] * msf_420[k];

        t_633[k] = pa_x[k] * lsg0_633[k]
                   + f_8 * lsf_423[k]
                   - f_6 * pc_x[k] * lsg1_633[k];

        t_634[k] = f_8 * lsf_342[k]
                   + f_3 * pc_y[k] * msf_422[k];
    }

#pragma omp simd aligned(t_635, t_636, t_637, t_638, pa_x, pc_x, lsg0_635, lsf_425, lsf_426, \
                         lsf_427, lsf_428, lsg1_635, msf_426, msf_427, \
                         msf_428 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_635[k] = pa_x[k] * lsg0_635[k]
                   + f_8 * lsf_425[k]
                   - f_6 * pc_x[k] * lsg1_635[k];

        t_636[k] = f_7 * lsf_426[k]
                   + f_3 * pc_x[k] * msf_426[k];

        t_637[k] = f_7 * lsf_427[k]
                   + f_3 * pc_x[k] * msf_427[k];

        t_638[k] = f_7 * lsf_428[k]
                   + f_3 * pc_x[k] * msf_428[k];
    }

#pragma omp simd aligned(t_639, t_640, t_641, t_642, pa_x, pc_x, pc_z, lsg0_640, lsg0_642, \
                         lsf_336, lsf_429, lsg1_640, lsg1_642, msf_426, \
                         msf_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_639[k] = f_7 * lsf_429[k]
                   + f_3 * pc_x[k] * msf_429[k];

        t_640[k] = pa_x[k] * lsg0_640[k]
                   - f_6 * pc_x[k] * lsg1_640[k];

        t_641[k] = f_13 * lsf_336[k]
                   + f_3 * pc_z[k] * msf_426[k];

        t_642[k] = pa_x[k] * lsg0_642[k]
                   - f_6 * pc_x[k] * lsg1_642[k];
    }

#pragma omp simd aligned(t_643, t_644, t_645, t_646, pa_x, pa_y, pc_x, pc_y, lsg0_525, \
                         lsg0_644, lsf_349, lsf_350, lsg1_525, lsg1_644, msf_429, \
                         msf_430 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_643[k] = f_8 * lsf_349[k]
                   + f_3 * pc_y[k] * msf_429[k];

        t_644[k] = pa_x[k] * lsg0_644[k]
                   - f_6 * pc_x[k] * lsg1_644[k];

        t_645[k] = pa_y[k] * lsg0_525[k]
                   - f_6 * pc_y[k] * lsg1_525[k];

        t_646[k] = f_7 * lsf_350[k]
                   + f_3 * pc_y[k] * msf_430[k];
    }

#pragma omp simd aligned(t_647, t_648, t_649, pa_x, pc_x, pc_y, pc_z, lsg0_648, lsf_340, \
                         lsf_352, lsf_433, lsg1_648, msf_430, msf_432 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_647[k] = f_12 * lsf_340[k]
                   + f_3 * pc_z[k] * msf_430[k];

        t_648[k] = pa_x[k] * lsg0_648[k]
                   + f_8 * lsf_433[k]
                   - f_6 * pc_x[k] * lsg1_648[k];

        t_649[k] = f_7 * lsf_352[k]
                   + f_3 * pc_y[k] * msf_432[k];
    }

#pragma omp simd aligned(t_650, t_651, t_652, t_653, pa_y, pc_x, pc_y, lsg0_530, lsf_436, \
                         lsf_437, lsf_438, lsg1_530, msf_436, msf_437, \
                         msf_438 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_650[k] = pa_y[k] * lsg0_530[k]
                   - f_6 * pc_y[k] * lsg1_530[k];

        t_651[k] = f_7 * lsf_436[k]
                   + f_3 * pc_x[k] * msf_436[k];

        t_652[k] = f_7 * lsf_437[k]
                   + f_3 * pc_x[k] * msf_437[k];

        t_653[k] = f_7 * lsf_438[k]
                   + f_3 * pc_x[k] * msf_438[k];
    }

#pragma omp simd aligned(t_654, t_655, t_656, t_657, pa_x, pc_x, pc_z, lsg0_655, lsg0_657, \
                         lsf_346, lsf_439, lsg1_655, lsg1_657, msf_436, \
                         msf_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_654[k] = f_7 * lsf_439[k]
                   + f_3 * pc_x[k] * msf_439[k];

        t_655[k] = pa_x[k] * lsg0_655[k]
                   - f_6 * pc_x[k] * lsg1_655[k];

        t_656[k] = f_12 * lsf_346[k]
                   + f_3 * pc_z[k] * msf_436[k];

        t_657[k] = pa_x[k] * lsg0_657[k]
                   - f_6 * pc_x[k] * lsg1_657[k];
    }

#pragma omp simd aligned(t_658, t_659, t_660, t_661, pa_x, pc_x, pc_y, lsg0_659, lsg0_660, \
                         lsf_359, lsf_440, lsg1_659, lsg1_660, msf_439, \
                         msf_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_658[k] = f_7 * lsf_359[k]
                   + f_3 * pc_y[k] * msf_439[k];

        t_659[k] = pa_x[k] * lsg0_659[k]
                   - f_6 * pc_x[k] * lsg1_659[k];

        t_660[k] = pa_x[k] * lsg0_660[k]
                   + f_16 * lsf_440[k]
                   - f_6 * pc_x[k] * lsg1_660[k];

        t_661[k] = f_3 * pc_y[k] * msf_440[k];
    }

#pragma omp simd aligned(t_662, t_663, t_664, pc_y, pc_z, lsf_350, msd0_264, msd1_264, \
                         msf_440, msf_441, msf_442 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_662[k] = f_9 * lsf_350[k]
                   + f_3 * pc_z[k] * msf_440[k];

        t_663[k] = f_4 * msd0_264[k]
                   - f_5 * msd1_264[k]
                   + f_3 * pc_y[k] * msf_441[k];

        t_664[k] = f_3 * pc_y[k] * msf_442[k];
    }

#pragma omp simd aligned(t_665, t_666, t_667, t_668, pa_x, pc_x, pc_y, lsg0_665, lsf_445, \
                         lsf_446, lsf_447, lsg1_665, msf_445, msf_446, \
                         msf_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_665[k] = pa_x[k] * lsg0_665[k]
                   + f_8 * lsf_445[k]
                   - f_6 * pc_x[k] * lsg1_665[k];

        t_666[k] = f_7 * lsf_446[k]
                   + f_3 * pc_x[k] * msf_446[k];

        t_667[k] = f_7 * lsf_447[k]
                   + f_3 * pc_x[k] * msf_447[k];

        t_668[k] = f_3 * pc_y[k] * msf_445[k];
    }

#pragma omp simd aligned(t_669, t_670, t_671, t_672, t_673, pa_x, pc_x, pc_y, lsg0_670, \
                         lsg0_671, lsg0_672, lsf_449, lsg1_670, lsg1_671, lsg1_672, \
                         msf_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_669[k] = f_7 * lsf_449[k]
                   + f_3 * pc_x[k] * msf_449[k];

        t_670[k] = pa_x[k] * lsg0_670[k]
                   - f_6 * pc_x[k] * lsg1_670[k];

        t_671[k] = pa_x[k] * lsg0_671[k]
                   - f_6 * pc_x[k] * lsg1_671[k];

        t_672[k] = pa_x[k] * lsg0_672[k]
                   - f_6 * pc_x[k] * lsg1_672[k];

        t_673[k] = f_3 * pc_y[k] * msf_449[k];
    }

#pragma omp simd aligned(t_674, t_675, t_676, t_677, pa_x, pc_x, pc_z, lsg0_674, lsg1_674, \
                         msd0_270, msd0_271, msd1_270, msd1_271, msf_450, \
                         msf_451 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_674[k] = pa_x[k] * lsg0_674[k]
                   - f_6 * pc_x[k] * lsg1_674[k];

        t_675[k] = f_1 * msd0_270[k]
                   - f_2 * msd1_270[k]
                   + f_3 * pc_x[k] * msf_450[k];

        t_676[k] = f_10 * msd0_271[k]
                   - f_11 * msd1_271[k]
                   + f_3 * pc_x[k] * msf_451[k];

        t_677[k] = f_3 * pc_z[k] * msf_450[k];
    }

#pragma omp simd aligned(t_678, t_679, t_680, t_681, t_682, pc_x, pc_z, msd0_273, msd0_275, \
                         msd1_273, msd1_275, msf_451, msf_453, msf_455, msf_456, \
                         msf_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_678[k] = f_4 * msd0_273[k]
                   - f_5 * msd1_273[k]
                   + f_3 * pc_x[k] * msf_453[k];

        t_679[k] = f_3 * pc_z[k] * msf_451[k];

        t_680[k] = f_4 * msd0_275[k]
                   - f_5 * msd1_275[k]
                   + f_3 * pc_x[k] * msf_455[k];

        t_681[k] = f_3 * pc_x[k] * msf_456[k];

        t_682[k] = f_3 * pc_x[k] * msf_457[k];
    }

#pragma omp simd aligned(t_683, t_684, t_685, t_686, t_687, pc_x, pc_y, pc_z, lsf_366, \
                         msd0_273, msd1_273, msf_456, msf_457, msf_458, \
                         msf_459 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_683[k] = f_3 * pc_x[k] * msf_458[k];

        t_684[k] = f_3 * pc_x[k] * msf_459[k];

        t_685[k] = f_0 * lsf_366[k]
                   + f_1 * msd0_273[k]
                   - f_2 * msd1_273[k]
                   + f_3 * pc_y[k] * msf_456[k];

        t_686[k] = f_3 * pc_z[k] * msf_456[k];

        t_687[k] = f_4 * msd0_273[k]
                   - f_5 * msd1_273[k]
                   + f_3 * pc_z[k] * msf_457[k];
    }

#pragma omp simd aligned(t_688, t_689, t_690, t_691, pa_z, pc_y, pc_z, lsg0_540, lsg0_541, \
                         lsf_369, lsg1_540, lsg1_541, msd0_275, msd1_275, \
                         msf_459 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_688[k] = f_0 * lsf_369[k]
                   + f_3 * pc_y[k] * msf_459[k];

        t_689[k] = f_1 * msd0_275[k]
                   - f_2 * msd1_275[k]
                   + f_3 * pc_z[k] * msf_459[k];

        t_690[k] = pa_z[k] * lsg0_540[k]
                   - f_6 * pc_z[k] * lsg1_540[k];

        t_691[k] = pa_z[k] * lsg0_541[k]
                   - f_6 * pc_z[k] * lsg1_541[k];
    }

#pragma omp simd aligned(t_692, t_693, t_694, pa_z, pc_x, pc_z, lsg0_543, lsg1_543, msd0_278, \
                         msd0_280, msd1_278, msd1_280, msf_462, \
                         msf_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_692[k] = f_10 * msd0_278[k]
                   - f_11 * msd1_278[k]
                   + f_3 * pc_x[k] * msf_462[k];

        t_693[k] = pa_z[k] * lsg0_543[k]
                   - f_6 * pc_z[k] * lsg1_543[k];

        t_694[k] = f_4 * msd0_280[k]
                   - f_5 * msd1_280[k]
                   + f_3 * pc_x[k] * msf_464[k];
    }

#pragma omp simd aligned(t_695, t_696, t_697, t_698, t_699, pc_x, msd0_281, msd1_281, msf_465, \
                         msf_466, msf_467, msf_468, msf_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_695[k] = f_4 * msd0_281[k]
                   - f_5 * msd1_281[k]
                   + f_3 * pc_x[k] * msf_465[k];

        t_696[k] = f_3 * pc_x[k] * msf_466[k];

        t_697[k] = f_3 * pc_x[k] * msf_467[k];

        t_698[k] = f_3 * pc_x[k] * msf_468[k];

        t_699[k] = f_3 * pc_x[k] * msf_469[k];
    }

#pragma omp simd aligned(t_700, t_701, t_702, t_703, pa_z, pc_y, pc_z, lsg0_550, lsg0_552, \
                         lsf_366, lsf_367, lsf_379, lsg1_550, lsg1_552, msf_466, \
                         msf_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_700[k] = pa_z[k] * lsg0_550[k]
                   - f_6 * pc_z[k] * lsg1_550[k];

        t_701[k] = f_7 * lsf_366[k]
                   + f_3 * pc_z[k] * msf_466[k];

        t_702[k] = pa_z[k] * lsg0_552[k]
                   + f_8 * lsf_367[k]
                   - f_6 * pc_z[k] * lsg1_552[k];

        t_703[k] = f_9 * lsf_379[k]
                   + f_3 * pc_y[k] * msf_469[k];
    }

#pragma omp simd aligned(t_704, t_705, t_706, pc_x, pc_z, lsf_369, msd0_281, msd0_282, \
                         msd0_283, msd1_281, msd1_282, msd1_283, msf_469, msf_470, \
                         msf_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_704[k] = f_7 * lsf_369[k]
                   + f_1 * msd0_281[k]
                   - f_2 * msd1_281[k]
                   + f_3 * pc_z[k] * msf_469[k];

        t_705[k] = f_1 * msd0_282[k]
                   - f_2 * msd1_282[k]
                   + f_3 * pc_x[k] * msf_470[k];

        t_706[k] = f_10 * msd0_283[k]
                   - f_11 * msd1_283[k]
                   + f_3 * pc_x[k] * msf_471[k];
    }

#pragma omp simd aligned(t_707, t_708, t_709, pc_x, msd0_284, msd0_285, msd0_286, msd1_284, \
                         msd1_285, msd1_286, msf_472, msf_473, \
                         msf_474 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_707[k] = f_10 * msd0_284[k]
                   - f_11 * msd1_284[k]
                   + f_3 * pc_x[k] * msf_472[k];

        t_708[k] = f_4 * msd0_285[k]
                   - f_5 * msd1_285[k]
                   + f_3 * pc_x[k] * msf_473[k];

        t_709[k] = f_4 * msd0_286[k]
                   - f_5 * msd1_286[k]
                   + f_3 * pc_x[k] * msf_474[k];
    }

#pragma omp simd aligned(t_710, t_711, t_712, t_713, t_714, pc_x, msd0_287, msd1_287, msf_475, \
                         msf_476, msf_477, msf_478, msf_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_710[k] = f_4 * msd0_287[k]
                   - f_5 * msd1_287[k]
                   + f_3 * pc_x[k] * msf_475[k];

        t_711[k] = f_3 * pc_x[k] * msf_476[k];

        t_712[k] = f_3 * pc_x[k] * msf_477[k];

        t_713[k] = f_3 * pc_x[k] * msf_478[k];

        t_714[k] = f_3 * pc_x[k] * msf_479[k];
    }

#pragma omp simd aligned(t_715, t_716, t_717, pc_y, pc_z, lsf_376, lsf_386, lsf_388, msd0_285, \
                         msd0_287, msd1_285, msd1_287, msf_476, \
                         msf_478 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_715[k] = f_12 * lsf_386[k]
                   + f_1 * msd0_285[k]
                   - f_2 * msd1_285[k]
                   + f_3 * pc_y[k] * msf_476[k];

        t_716[k] = f_8 * lsf_376[k]
                   + f_3 * pc_z[k] * msf_476[k];

        t_717[k] = f_12 * lsf_388[k]
                   + f_4 * msd0_287[k]
                   - f_5 * msd1_287[k]
                   + f_3 * pc_y[k] * msf_478[k];
    }

#pragma omp simd aligned(t_718, t_719, t_720, pc_x, pc_y, pc_z, lsf_379, lsf_389, msd0_287, \
                         msd0_288, msd1_287, msd1_288, msf_479, \
                         msf_480 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_718[k] = f_12 * lsf_389[k]
                   + f_3 * pc_y[k] * msf_479[k];

        t_719[k] = f_8 * lsf_379[k]
                   + f_1 * msd0_287[k]
                   - f_2 * msd1_287[k]
                   + f_3 * pc_z[k] * msf_479[k];

        t_720[k] = f_1 * msd0_288[k]
                   - f_2 * msd1_288[k]
                   + f_3 * pc_x[k] * msf_480[k];
    }

#pragma omp simd aligned(t_721, t_722, t_723, pc_x, msd0_289, msd0_290, msd0_291, msd1_289, \
                         msd1_290, msd1_291, msf_481, msf_482, \
                         msf_483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_721[k] = f_10 * msd0_289[k]
                   - f_11 * msd1_289[k]
                   + f_3 * pc_x[k] * msf_481[k];

        t_722[k] = f_10 * msd0_290[k]
                   - f_11 * msd1_290[k]
                   + f_3 * pc_x[k] * msf_482[k];

        t_723[k] = f_4 * msd0_291[k]
                   - f_5 * msd1_291[k]
                   + f_3 * pc_x[k] * msf_483[k];
    }

#pragma omp simd aligned(t_724, t_725, t_726, t_727, t_728, pc_x, msd0_292, msd0_293, \
                         msd1_292, msd1_293, msf_484, msf_485, msf_486, msf_487, \
                         msf_488 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_724[k] = f_4 * msd0_292[k]
                   - f_5 * msd1_292[k]
                   + f_3 * pc_x[k] * msf_484[k];

        t_725[k] = f_4 * msd0_293[k]
                   - f_5 * msd1_293[k]
                   + f_3 * pc_x[k] * msf_485[k];

        t_726[k] = f_3 * pc_x[k] * msf_486[k];

        t_727[k] = f_3 * pc_x[k] * msf_487[k];

        t_728[k] = f_3 * pc_x[k] * msf_488[k];
    }

#pragma omp simd aligned(t_729, t_730, t_731, pc_x, pc_y, pc_z, lsf_386, lsf_396, msd0_291, \
                         msd1_291, msf_486, msf_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_729[k] = f_3 * pc_x[k] * msf_489[k];

        t_730[k] = f_13 * lsf_396[k]
                   + f_1 * msd0_291[k]
                   - f_2 * msd1_291[k]
                   + f_3 * pc_y[k] * msf_486[k];

        t_731[k] = f_14 * lsf_386[k]
                   + f_3 * pc_z[k] * msf_486[k];
    }
}

static auto
compute_prim_msg_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsg0,
                                                          const size_t lsf, const size_t lsg1,
                                                          const size_t msd0, const size_t msd1,
                                                          const size_t msf, const size_t ncols,
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
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 3.5 / q;
    const auto f_13 = 3.0 / q;
    const auto f_14 = 1.5 / q;
    const auto f_15 = 2.5 / q;
    const auto f_16 = 2.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *lsg0_660 = buffer.data(lsg0 + 660);
    const auto *lsg0_662 = buffer.data(lsg0 + 662);
    const auto *lsg0_665 = buffer.data(lsg0 + 665);
    const auto *lsg0_670 = buffer.data(lsg0 + 670);
    const auto *lsg0_672 = buffer.data(lsg0 + 672);
    const auto *lsg0_674 = buffer.data(lsg0 + 674);

    const auto *lsf_389 = buffer.data(lsf + 389);
    const auto *lsf_396 = buffer.data(lsf + 396);
    const auto *lsf_398 = buffer.data(lsf + 398);
    const auto *lsf_399 = buffer.data(lsf + 399);
    const auto *lsf_406 = buffer.data(lsf + 406);
    const auto *lsf_408 = buffer.data(lsf + 408);
    const auto *lsf_409 = buffer.data(lsf + 409);
    const auto *lsf_416 = buffer.data(lsf + 416);
    const auto *lsf_418 = buffer.data(lsf + 418);
    const auto *lsf_419 = buffer.data(lsf + 419);
    const auto *lsf_426 = buffer.data(lsf + 426);
    const auto *lsf_428 = buffer.data(lsf + 428);
    const auto *lsf_429 = buffer.data(lsf + 429);
    const auto *lsf_436 = buffer.data(lsf + 436);
    const auto *lsf_438 = buffer.data(lsf + 438);
    const auto *lsf_439 = buffer.data(lsf + 439);
    const auto *lsf_446 = buffer.data(lsf + 446);
    const auto *lsf_448 = buffer.data(lsf + 448);
    const auto *lsf_449 = buffer.data(lsf + 449);

    const auto *lsg1_660 = buffer.data(lsg1 + 660);
    const auto *lsg1_662 = buffer.data(lsg1 + 662);
    const auto *lsg1_665 = buffer.data(lsg1 + 665);
    const auto *lsg1_670 = buffer.data(lsg1 + 670);
    const auto *lsg1_672 = buffer.data(lsg1 + 672);
    const auto *lsg1_674 = buffer.data(lsg1 + 674);

    const auto *msd0_293 = buffer.data(msd0 + 293);
    const auto *msd0_294 = buffer.data(msd0 + 294);
    const auto *msd0_295 = buffer.data(msd0 + 295);
    const auto *msd0_296 = buffer.data(msd0 + 296);
    const auto *msd0_297 = buffer.data(msd0 + 297);
    const auto *msd0_298 = buffer.data(msd0 + 298);
    const auto *msd0_299 = buffer.data(msd0 + 299);
    const auto *msd0_300 = buffer.data(msd0 + 300);
    const auto *msd0_301 = buffer.data(msd0 + 301);
    const auto *msd0_302 = buffer.data(msd0 + 302);
    const auto *msd0_303 = buffer.data(msd0 + 303);
    const auto *msd0_304 = buffer.data(msd0 + 304);
    const auto *msd0_305 = buffer.data(msd0 + 305);
    const auto *msd0_306 = buffer.data(msd0 + 306);
    const auto *msd0_307 = buffer.data(msd0 + 307);
    const auto *msd0_308 = buffer.data(msd0 + 308);
    const auto *msd0_309 = buffer.data(msd0 + 309);
    const auto *msd0_310 = buffer.data(msd0 + 310);
    const auto *msd0_311 = buffer.data(msd0 + 311);
    const auto *msd0_312 = buffer.data(msd0 + 312);
    const auto *msd0_313 = buffer.data(msd0 + 313);
    const auto *msd0_314 = buffer.data(msd0 + 314);
    const auto *msd0_315 = buffer.data(msd0 + 315);
    const auto *msd0_316 = buffer.data(msd0 + 316);
    const auto *msd0_317 = buffer.data(msd0 + 317);
    const auto *msd0_319 = buffer.data(msd0 + 319);
    const auto *msd0_321 = buffer.data(msd0 + 321);
    const auto *msd0_322 = buffer.data(msd0 + 322);
    const auto *msd0_324 = buffer.data(msd0 + 324);
    const auto *msd0_326 = buffer.data(msd0 + 326);
    const auto *msd0_327 = buffer.data(msd0 + 327);
    const auto *msd0_328 = buffer.data(msd0 + 328);
    const auto *msd0_329 = buffer.data(msd0 + 329);

    const auto *msd1_293 = buffer.data(msd1 + 293);
    const auto *msd1_294 = buffer.data(msd1 + 294);
    const auto *msd1_295 = buffer.data(msd1 + 295);
    const auto *msd1_296 = buffer.data(msd1 + 296);
    const auto *msd1_297 = buffer.data(msd1 + 297);
    const auto *msd1_298 = buffer.data(msd1 + 298);
    const auto *msd1_299 = buffer.data(msd1 + 299);
    const auto *msd1_300 = buffer.data(msd1 + 300);
    const auto *msd1_301 = buffer.data(msd1 + 301);
    const auto *msd1_302 = buffer.data(msd1 + 302);
    const auto *msd1_303 = buffer.data(msd1 + 303);
    const auto *msd1_304 = buffer.data(msd1 + 304);
    const auto *msd1_305 = buffer.data(msd1 + 305);
    const auto *msd1_306 = buffer.data(msd1 + 306);
    const auto *msd1_307 = buffer.data(msd1 + 307);
    const auto *msd1_308 = buffer.data(msd1 + 308);
    const auto *msd1_309 = buffer.data(msd1 + 309);
    const auto *msd1_310 = buffer.data(msd1 + 310);
    const auto *msd1_311 = buffer.data(msd1 + 311);
    const auto *msd1_312 = buffer.data(msd1 + 312);
    const auto *msd1_313 = buffer.data(msd1 + 313);
    const auto *msd1_314 = buffer.data(msd1 + 314);
    const auto *msd1_315 = buffer.data(msd1 + 315);
    const auto *msd1_316 = buffer.data(msd1 + 316);
    const auto *msd1_317 = buffer.data(msd1 + 317);
    const auto *msd1_319 = buffer.data(msd1 + 319);
    const auto *msd1_321 = buffer.data(msd1 + 321);
    const auto *msd1_322 = buffer.data(msd1 + 322);
    const auto *msd1_324 = buffer.data(msd1 + 324);
    const auto *msd1_326 = buffer.data(msd1 + 326);
    const auto *msd1_327 = buffer.data(msd1 + 327);
    const auto *msd1_328 = buffer.data(msd1 + 328);
    const auto *msd1_329 = buffer.data(msd1 + 329);

    const auto *msf_488 = buffer.data(msf + 488);
    const auto *msf_489 = buffer.data(msf + 489);
    const auto *msf_490 = buffer.data(msf + 490);
    const auto *msf_491 = buffer.data(msf + 491);
    const auto *msf_492 = buffer.data(msf + 492);
    const auto *msf_493 = buffer.data(msf + 493);
    const auto *msf_494 = buffer.data(msf + 494);
    const auto *msf_495 = buffer.data(msf + 495);
    const auto *msf_496 = buffer.data(msf + 496);
    const auto *msf_497 = buffer.data(msf + 497);
    const auto *msf_498 = buffer.data(msf + 498);
    const auto *msf_499 = buffer.data(msf + 499);
    const auto *msf_500 = buffer.data(msf + 500);
    const auto *msf_501 = buffer.data(msf + 501);
    const auto *msf_502 = buffer.data(msf + 502);
    const auto *msf_503 = buffer.data(msf + 503);
    const auto *msf_504 = buffer.data(msf + 504);
    const auto *msf_505 = buffer.data(msf + 505);
    const auto *msf_506 = buffer.data(msf + 506);
    const auto *msf_507 = buffer.data(msf + 507);
    const auto *msf_508 = buffer.data(msf + 508);
    const auto *msf_509 = buffer.data(msf + 509);
    const auto *msf_510 = buffer.data(msf + 510);
    const auto *msf_511 = buffer.data(msf + 511);
    const auto *msf_512 = buffer.data(msf + 512);
    const auto *msf_513 = buffer.data(msf + 513);
    const auto *msf_514 = buffer.data(msf + 514);
    const auto *msf_515 = buffer.data(msf + 515);
    const auto *msf_516 = buffer.data(msf + 516);
    const auto *msf_517 = buffer.data(msf + 517);
    const auto *msf_518 = buffer.data(msf + 518);
    const auto *msf_519 = buffer.data(msf + 519);
    const auto *msf_520 = buffer.data(msf + 520);
    const auto *msf_521 = buffer.data(msf + 521);
    const auto *msf_522 = buffer.data(msf + 522);
    const auto *msf_523 = buffer.data(msf + 523);
    const auto *msf_524 = buffer.data(msf + 524);
    const auto *msf_525 = buffer.data(msf + 525);
    const auto *msf_526 = buffer.data(msf + 526);
    const auto *msf_527 = buffer.data(msf + 527);
    const auto *msf_528 = buffer.data(msf + 528);
    const auto *msf_529 = buffer.data(msf + 529);
    const auto *msf_531 = buffer.data(msf + 531);
    const auto *msf_533 = buffer.data(msf + 533);
    const auto *msf_534 = buffer.data(msf + 534);
    const auto *msf_536 = buffer.data(msf + 536);
    const auto *msf_537 = buffer.data(msf + 537);
    const auto *msf_538 = buffer.data(msf + 538);
    const auto *msf_539 = buffer.data(msf + 539);
    const auto *msf_540 = buffer.data(msf + 540);
    const auto *msf_542 = buffer.data(msf + 542);
    const auto *msf_543 = buffer.data(msf + 543);
    const auto *msf_545 = buffer.data(msf + 545);
    const auto *msf_546 = buffer.data(msf + 546);
    const auto *msf_547 = buffer.data(msf + 547);
    const auto *msf_548 = buffer.data(msf + 548);
    const auto *msf_549 = buffer.data(msf + 549);

#pragma omp simd aligned(t_732, t_733, t_734, pc_y, pc_z, lsf_389, lsf_398, lsf_399, msd0_293, \
                         msd1_293, msf_488, msf_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_732[k] = f_13 * lsf_398[k]
                   + f_4 * msd0_293[k]
                   - f_5 * msd1_293[k]
                   + f_3 * pc_y[k] * msf_488[k];

        t_733[k] = f_13 * lsf_399[k]
                   + f_3 * pc_y[k] * msf_489[k];

        t_734[k] = f_14 * lsf_389[k]
                   + f_1 * msd0_293[k]
                   - f_2 * msd1_293[k]
                   + f_3 * pc_z[k] * msf_489[k];
    }

#pragma omp simd aligned(t_735, t_736, t_737, pc_x, msd0_294, msd0_295, msd0_296, msd1_294, \
                         msd1_295, msd1_296, msf_490, msf_491, \
                         msf_492 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_735[k] = f_1 * msd0_294[k]
                   - f_2 * msd1_294[k]
                   + f_3 * pc_x[k] * msf_490[k];

        t_736[k] = f_10 * msd0_295[k]
                   - f_11 * msd1_295[k]
                   + f_3 * pc_x[k] * msf_491[k];

        t_737[k] = f_10 * msd0_296[k]
                   - f_11 * msd1_296[k]
                   + f_3 * pc_x[k] * msf_492[k];
    }

#pragma omp simd aligned(t_738, t_739, t_740, t_741, pc_x, msd0_297, msd0_298, msd0_299, \
                         msd1_297, msd1_298, msd1_299, msf_493, msf_494, msf_495, \
                         msf_496 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_738[k] = f_4 * msd0_297[k]
                   - f_5 * msd1_297[k]
                   + f_3 * pc_x[k] * msf_493[k];

        t_739[k] = f_4 * msd0_298[k]
                   - f_5 * msd1_298[k]
                   + f_3 * pc_x[k] * msf_494[k];

        t_740[k] = f_4 * msd0_299[k]
                   - f_5 * msd1_299[k]
                   + f_3 * pc_x[k] * msf_495[k];

        t_741[k] = f_3 * pc_x[k] * msf_496[k];
    }

#pragma omp simd aligned(t_742, t_743, t_744, t_745, t_746, pc_x, pc_y, pc_z, lsf_396, \
                         lsf_406, msd0_297, msd1_297, msf_496, msf_497, msf_498, \
                         msf_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_742[k] = f_3 * pc_x[k] * msf_497[k];

        t_743[k] = f_3 * pc_x[k] * msf_498[k];

        t_744[k] = f_3 * pc_x[k] * msf_499[k];

        t_745[k] = f_15 * lsf_406[k]
                   + f_1 * msd0_297[k]
                   - f_2 * msd1_297[k]
                   + f_3 * pc_y[k] * msf_496[k];

        t_746[k] = f_16 * lsf_396[k]
                   + f_3 * pc_z[k] * msf_496[k];
    }

#pragma omp simd aligned(t_747, t_748, t_749, pc_y, pc_z, lsf_399, lsf_408, lsf_409, msd0_299, \
                         msd1_299, msf_498, msf_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_747[k] = f_15 * lsf_408[k]
                   + f_4 * msd0_299[k]
                   - f_5 * msd1_299[k]
                   + f_3 * pc_y[k] * msf_498[k];

        t_748[k] = f_15 * lsf_409[k]
                   + f_3 * pc_y[k] * msf_499[k];

        t_749[k] = f_16 * lsf_399[k]
                   + f_1 * msd0_299[k]
                   - f_2 * msd1_299[k]
                   + f_3 * pc_z[k] * msf_499[k];
    }

#pragma omp simd aligned(t_750, t_751, t_752, pc_x, msd0_300, msd0_301, msd0_302, msd1_300, \
                         msd1_301, msd1_302, msf_500, msf_501, \
                         msf_502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_750[k] = f_1 * msd0_300[k]
                   - f_2 * msd1_300[k]
                   + f_3 * pc_x[k] * msf_500[k];

        t_751[k] = f_10 * msd0_301[k]
                   - f_11 * msd1_301[k]
                   + f_3 * pc_x[k] * msf_501[k];

        t_752[k] = f_10 * msd0_302[k]
                   - f_11 * msd1_302[k]
                   + f_3 * pc_x[k] * msf_502[k];
    }

#pragma omp simd aligned(t_753, t_754, t_755, t_756, pc_x, msd0_303, msd0_304, msd0_305, \
                         msd1_303, msd1_304, msd1_305, msf_503, msf_504, msf_505, \
                         msf_506 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_753[k] = f_4 * msd0_303[k]
                   - f_5 * msd1_303[k]
                   + f_3 * pc_x[k] * msf_503[k];

        t_754[k] = f_4 * msd0_304[k]
                   - f_5 * msd1_304[k]
                   + f_3 * pc_x[k] * msf_504[k];

        t_755[k] = f_4 * msd0_305[k]
                   - f_5 * msd1_305[k]
                   + f_3 * pc_x[k] * msf_505[k];

        t_756[k] = f_3 * pc_x[k] * msf_506[k];
    }

#pragma omp simd aligned(t_757, t_758, t_759, t_760, t_761, pc_x, pc_y, pc_z, lsf_406, \
                         lsf_416, msd0_303, msd1_303, msf_506, msf_507, msf_508, \
                         msf_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_757[k] = f_3 * pc_x[k] * msf_507[k];

        t_758[k] = f_3 * pc_x[k] * msf_508[k];

        t_759[k] = f_3 * pc_x[k] * msf_509[k];

        t_760[k] = f_16 * lsf_416[k]
                   + f_1 * msd0_303[k]
                   - f_2 * msd1_303[k]
                   + f_3 * pc_y[k] * msf_506[k];

        t_761[k] = f_15 * lsf_406[k]
                   + f_3 * pc_z[k] * msf_506[k];
    }

#pragma omp simd aligned(t_762, t_763, t_764, pc_y, pc_z, lsf_409, lsf_418, lsf_419, msd0_305, \
                         msd1_305, msf_508, msf_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_762[k] = f_16 * lsf_418[k]
                   + f_4 * msd0_305[k]
                   - f_5 * msd1_305[k]
                   + f_3 * pc_y[k] * msf_508[k];

        t_763[k] = f_16 * lsf_419[k]
                   + f_3 * pc_y[k] * msf_509[k];

        t_764[k] = f_15 * lsf_409[k]
                   + f_1 * msd0_305[k]
                   - f_2 * msd1_305[k]
                   + f_3 * pc_z[k] * msf_509[k];
    }

#pragma omp simd aligned(t_765, t_766, t_767, pc_x, msd0_306, msd0_307, msd0_308, msd1_306, \
                         msd1_307, msd1_308, msf_510, msf_511, \
                         msf_512 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_765[k] = f_1 * msd0_306[k]
                   - f_2 * msd1_306[k]
                   + f_3 * pc_x[k] * msf_510[k];

        t_766[k] = f_10 * msd0_307[k]
                   - f_11 * msd1_307[k]
                   + f_3 * pc_x[k] * msf_511[k];

        t_767[k] = f_10 * msd0_308[k]
                   - f_11 * msd1_308[k]
                   + f_3 * pc_x[k] * msf_512[k];
    }

#pragma omp simd aligned(t_768, t_769, t_770, t_771, pc_x, msd0_309, msd0_310, msd0_311, \
                         msd1_309, msd1_310, msd1_311, msf_513, msf_514, msf_515, \
                         msf_516 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_768[k] = f_4 * msd0_309[k]
                   - f_5 * msd1_309[k]
                   + f_3 * pc_x[k] * msf_513[k];

        t_769[k] = f_4 * msd0_310[k]
                   - f_5 * msd1_310[k]
                   + f_3 * pc_x[k] * msf_514[k];

        t_770[k] = f_4 * msd0_311[k]
                   - f_5 * msd1_311[k]
                   + f_3 * pc_x[k] * msf_515[k];

        t_771[k] = f_3 * pc_x[k] * msf_516[k];
    }

#pragma omp simd aligned(t_772, t_773, t_774, t_775, t_776, pc_x, pc_y, pc_z, lsf_416, \
                         lsf_426, msd0_309, msd1_309, msf_516, msf_517, msf_518, \
                         msf_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_772[k] = f_3 * pc_x[k] * msf_517[k];

        t_773[k] = f_3 * pc_x[k] * msf_518[k];

        t_774[k] = f_3 * pc_x[k] * msf_519[k];

        t_775[k] = f_14 * lsf_426[k]
                   + f_1 * msd0_309[k]
                   - f_2 * msd1_309[k]
                   + f_3 * pc_y[k] * msf_516[k];

        t_776[k] = f_13 * lsf_416[k]
                   + f_3 * pc_z[k] * msf_516[k];
    }

#pragma omp simd aligned(t_777, t_778, t_779, pc_y, pc_z, lsf_419, lsf_428, lsf_429, msd0_311, \
                         msd1_311, msf_518, msf_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_777[k] = f_14 * lsf_428[k]
                   + f_4 * msd0_311[k]
                   - f_5 * msd1_311[k]
                   + f_3 * pc_y[k] * msf_518[k];

        t_778[k] = f_14 * lsf_429[k]
                   + f_3 * pc_y[k] * msf_519[k];

        t_779[k] = f_13 * lsf_419[k]
                   + f_1 * msd0_311[k]
                   - f_2 * msd1_311[k]
                   + f_3 * pc_z[k] * msf_519[k];
    }

#pragma omp simd aligned(t_780, t_781, t_782, pc_x, msd0_312, msd0_313, msd0_314, msd1_312, \
                         msd1_313, msd1_314, msf_520, msf_521, \
                         msf_522 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_780[k] = f_1 * msd0_312[k]
                   - f_2 * msd1_312[k]
                   + f_3 * pc_x[k] * msf_520[k];

        t_781[k] = f_10 * msd0_313[k]
                   - f_11 * msd1_313[k]
                   + f_3 * pc_x[k] * msf_521[k];

        t_782[k] = f_10 * msd0_314[k]
                   - f_11 * msd1_314[k]
                   + f_3 * pc_x[k] * msf_522[k];
    }

#pragma omp simd aligned(t_783, t_784, t_785, t_786, pc_x, msd0_315, msd0_316, msd0_317, \
                         msd1_315, msd1_316, msd1_317, msf_523, msf_524, msf_525, \
                         msf_526 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_783[k] = f_4 * msd0_315[k]
                   - f_5 * msd1_315[k]
                   + f_3 * pc_x[k] * msf_523[k];

        t_784[k] = f_4 * msd0_316[k]
                   - f_5 * msd1_316[k]
                   + f_3 * pc_x[k] * msf_524[k];

        t_785[k] = f_4 * msd0_317[k]
                   - f_5 * msd1_317[k]
                   + f_3 * pc_x[k] * msf_525[k];

        t_786[k] = f_3 * pc_x[k] * msf_526[k];
    }

#pragma omp simd aligned(t_787, t_788, t_789, t_790, t_791, pc_x, pc_y, pc_z, lsf_426, \
                         lsf_436, msd0_315, msd1_315, msf_526, msf_527, msf_528, \
                         msf_529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_787[k] = f_3 * pc_x[k] * msf_527[k];

        t_788[k] = f_3 * pc_x[k] * msf_528[k];

        t_789[k] = f_3 * pc_x[k] * msf_529[k];

        t_790[k] = f_8 * lsf_436[k]
                   + f_1 * msd0_315[k]
                   - f_2 * msd1_315[k]
                   + f_3 * pc_y[k] * msf_526[k];

        t_791[k] = f_12 * lsf_426[k]
                   + f_3 * pc_z[k] * msf_526[k];
    }

#pragma omp simd aligned(t_792, t_793, t_794, t_795, pa_y, pc_y, pc_z, lsg0_660, lsf_429, \
                         lsf_438, lsf_439, lsg1_660, msd0_317, msd1_317, msf_528, \
                         msf_529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_792[k] = f_8 * lsf_438[k]
                   + f_4 * msd0_317[k]
                   - f_5 * msd1_317[k]
                   + f_3 * pc_y[k] * msf_528[k];

        t_793[k] = f_8 * lsf_439[k]
                   + f_3 * pc_y[k] * msf_529[k];

        t_794[k] = f_12 * lsf_429[k]
                   + f_1 * msd0_317[k]
                   - f_2 * msd1_317[k]
                   + f_3 * pc_z[k] * msf_529[k];

        t_795[k] = pa_y[k] * lsg0_660[k]
                   - f_6 * pc_y[k] * lsg1_660[k];
    }

#pragma omp simd aligned(t_796, t_797, t_798, pa_y, pc_x, pc_y, lsg0_662, lsg1_662, msd0_319, \
                         msd0_321, msd1_319, msd1_321, msf_531, \
                         msf_533 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_796[k] = f_10 * msd0_319[k]
                   - f_11 * msd1_319[k]
                   + f_3 * pc_x[k] * msf_531[k];

        t_797[k] = pa_y[k] * lsg0_662[k]
                   - f_6 * pc_y[k] * lsg1_662[k];

        t_798[k] = f_4 * msd0_321[k]
                   - f_5 * msd1_321[k]
                   + f_3 * pc_x[k] * msf_533[k];
    }

#pragma omp simd aligned(t_799, t_800, t_801, t_802, t_803, pa_y, pc_x, pc_y, lsg0_665, \
                         lsg1_665, msd0_322, msd1_322, msf_534, msf_536, msf_537, \
                         msf_538 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_799[k] = f_4 * msd0_322[k]
                   - f_5 * msd1_322[k]
                   + f_3 * pc_x[k] * msf_534[k];

        t_800[k] = pa_y[k] * lsg0_665[k]
                   - f_6 * pc_y[k] * lsg1_665[k];

        t_801[k] = f_3 * pc_x[k] * msf_536[k];

        t_802[k] = f_3 * pc_x[k] * msf_537[k];

        t_803[k] = f_3 * pc_x[k] * msf_538[k];
    }

#pragma omp simd aligned(t_804, t_805, t_806, pa_y, pc_x, pc_y, pc_z, lsg0_670, lsf_436, \
                         lsf_446, lsg1_670, msf_536, msf_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_804[k] = f_3 * pc_x[k] * msf_539[k];

        t_805[k] = pa_y[k] * lsg0_670[k]
                   + f_16 * lsf_446[k]
                   - f_6 * pc_y[k] * lsg1_670[k];

        t_806[k] = f_9 * lsf_436[k]
                   + f_3 * pc_z[k] * msf_536[k];
    }

#pragma omp simd aligned(t_807, t_808, t_809, pa_y, pc_y, lsg0_672, lsg0_674, lsf_448, \
                         lsf_449, lsg1_672, lsg1_674, msf_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_807[k] = pa_y[k] * lsg0_672[k]
                   + f_8 * lsf_448[k]
                   - f_6 * pc_y[k] * lsg1_672[k];

        t_808[k] = f_7 * lsf_449[k]
                   + f_3 * pc_y[k] * msf_539[k];

        t_809[k] = pa_y[k] * lsg0_674[k]
                   - f_6 * pc_y[k] * lsg1_674[k];
    }

#pragma omp simd aligned(t_810, t_811, t_812, t_813, t_814, pc_x, pc_y, msd0_324, msd0_326, \
                         msd0_327, msd1_324, msd1_326, msd1_327, msf_540, msf_542, \
                         msf_543 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_810[k] = f_1 * msd0_324[k]
                   - f_2 * msd1_324[k]
                   + f_3 * pc_x[k] * msf_540[k];

        t_811[k] = f_3 * pc_y[k] * msf_540[k];

        t_812[k] = f_10 * msd0_326[k]
                   - f_11 * msd1_326[k]
                   + f_3 * pc_x[k] * msf_542[k];

        t_813[k] = f_4 * msd0_327[k]
                   - f_5 * msd1_327[k]
                   + f_3 * pc_x[k] * msf_543[k];

        t_814[k] = f_3 * pc_y[k] * msf_542[k];
    }

#pragma omp simd aligned(t_815, t_816, t_817, t_818, t_819, pc_x, msd0_329, msd1_329, msf_545, \
                         msf_546, msf_547, msf_548, msf_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_815[k] = f_4 * msd0_329[k]
                   - f_5 * msd1_329[k]
                   + f_3 * pc_x[k] * msf_545[k];

        t_816[k] = f_3 * pc_x[k] * msf_546[k];

        t_817[k] = f_3 * pc_x[k] * msf_547[k];

        t_818[k] = f_3 * pc_x[k] * msf_548[k];

        t_819[k] = f_3 * pc_x[k] * msf_549[k];
    }

#pragma omp simd aligned(t_820, t_821, t_822, t_823, pc_y, msd0_327, msd0_328, msd0_329, \
                         msd1_327, msd1_328, msd1_329, msf_546, msf_547, msf_548, \
                         msf_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_820[k] = f_1 * msd0_327[k]
                   - f_2 * msd1_327[k]
                   + f_3 * pc_y[k] * msf_546[k];

        t_821[k] = f_10 * msd0_328[k]
                   - f_11 * msd1_328[k]
                   + f_3 * pc_y[k] * msf_547[k];

        t_822[k] = f_4 * msd0_329[k]
                   - f_5 * msd1_329[k]
                   + f_3 * pc_y[k] * msf_548[k];

        t_823[k] = f_3 * pc_y[k] * msf_549[k];
    }

#pragma omp simd aligned(t_824, pc_z, lsf_449, msd0_329, msd1_329, \
                         msf_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_824[k] = f_0 * lsf_449[k]
                   + f_1 * msd0_329[k]
                   - f_2 * msd1_329[k]
                   + f_3 * pc_z[k] * msf_549[k];
    }
}

auto
compute_prim_msg_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t lsg0, const size_t lsf,
                                                   const size_t lsg1, const size_t msd0,
                                                   const size_t msd1, const size_t msf,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_msg_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, lsg0, lsf,
                                                              lsg1, msd0, msd1, msf, ncols,
                                                              gamma, p, q);

    compute_prim_msg_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, lsg0, lsf,
                                                              lsg1, msd0, msd1, msf, ncols,
                                                              gamma, p, q);

    compute_prim_msg_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, lsg0, lsf,
                                                              lsg1, msd0, msd1, msf, ncols,
                                                              gamma, p, q);

    compute_prim_msg_three_center_electron_repulsion_0_piece3(buffer, target, pa, pc, lsg0, lsf,
                                                              lsg1, msd0, msd1, msf, ncols,
                                                              gamma, p, q);

    compute_prim_msg_three_center_electron_repulsion_0_piece4(buffer, target, pa, pc, lsg0, lsf,
                                                              lsg1, msd0, msd1, msf, ncols,
                                                              gamma, p, q);

    compute_prim_msg_three_center_electron_repulsion_0_piece5(buffer, target, pa, pc, lsg0, lsf,
                                                              lsg1, msd0, msd1, msf, ncols,
                                                              gamma, p, q);

    compute_prim_msg_three_center_electron_repulsion_0_piece6(buffer, target, pa, pc, lsg0, lsf,
                                                              lsg1, msd0, msd1, msf, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
