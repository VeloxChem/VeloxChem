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


#include "SimdThreeCenterElectronRepulsionVrrRecKSD.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_ksd_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t isd0,
                                                          const size_t isp, const size_t isd1,
                                                          const size_t kss0, const size_t kss1,
                                                          const size_t ksp, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / q;
    const auto f_1 = 0.5 / gamma;
    const auto f_2 = 0.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 3.0 / q;
    const auto f_6 = 0.5 / q;
    const auto f_7 = 2.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 2.0 / q;
    const auto f_10 = 1.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *isd0_0 = buffer.data(isd0 + 0);
    const auto *isd0_3 = buffer.data(isd0 + 3);
    const auto *isd0_5 = buffer.data(isd0 + 5);
    const auto *isd0_9 = buffer.data(isd0 + 9);
    const auto *isd0_12 = buffer.data(isd0 + 12);
    const auto *isd0_17 = buffer.data(isd0 + 17);
    const auto *isd0_18 = buffer.data(isd0 + 18);
    const auto *isd0_21 = buffer.data(isd0 + 21);
    const auto *isd0_30 = buffer.data(isd0 + 30);
    const auto *isd0_35 = buffer.data(isd0 + 35);
    const auto *isd0_36 = buffer.data(isd0 + 36);
    const auto *isd0_39 = buffer.data(isd0 + 39);
    const auto *isd0_54 = buffer.data(isd0 + 54);
    const auto *isd0_59 = buffer.data(isd0 + 59);
    const auto *isd0_60 = buffer.data(isd0 + 60);
    const auto *isd0_63 = buffer.data(isd0 + 63);
    const auto *isd0_84 = buffer.data(isd0 + 84);
    const auto *isd0_89 = buffer.data(isd0 + 89);

    const auto *isp_0 = buffer.data(isp + 0);
    const auto *isp_1 = buffer.data(isp + 1);
    const auto *isp_2 = buffer.data(isp + 2);
    const auto *isp_4 = buffer.data(isp + 4);
    const auto *isp_8 = buffer.data(isp + 8);
    const auto *isp_9 = buffer.data(isp + 9);
    const auto *isp_10 = buffer.data(isp + 10);
    const auto *isp_11 = buffer.data(isp + 11);
    const auto *isp_13 = buffer.data(isp + 13);
    const auto *isp_14 = buffer.data(isp + 14);
    const auto *isp_15 = buffer.data(isp + 15);
    const auto *isp_16 = buffer.data(isp + 16);
    const auto *isp_17 = buffer.data(isp + 17);
    const auto *isp_18 = buffer.data(isp + 18);
    const auto *isp_19 = buffer.data(isp + 19);
    const auto *isp_20 = buffer.data(isp + 20);
    const auto *isp_22 = buffer.data(isp + 22);
    const auto *isp_23 = buffer.data(isp + 23);
    const auto *isp_25 = buffer.data(isp + 25);
    const auto *isp_26 = buffer.data(isp + 26);
    const auto *isp_27 = buffer.data(isp + 27);
    const auto *isp_28 = buffer.data(isp + 28);
    const auto *isp_29 = buffer.data(isp + 29);
    const auto *isp_30 = buffer.data(isp + 30);
    const auto *isp_31 = buffer.data(isp + 31);
    const auto *isp_32 = buffer.data(isp + 32);
    const auto *isp_34 = buffer.data(isp + 34);
    const auto *isp_35 = buffer.data(isp + 35);
    const auto *isp_36 = buffer.data(isp + 36);
    const auto *isp_37 = buffer.data(isp + 37);
    const auto *isp_38 = buffer.data(isp + 38);
    const auto *isp_40 = buffer.data(isp + 40);
    const auto *isp_41 = buffer.data(isp + 41);
    const auto *isp_42 = buffer.data(isp + 42);
    const auto *isp_43 = buffer.data(isp + 43);
    const auto *isp_44 = buffer.data(isp + 44);
    const auto *isp_45 = buffer.data(isp + 45);
    const auto *isp_46 = buffer.data(isp + 46);
    const auto *isp_49 = buffer.data(isp + 49);
    const auto *isp_50 = buffer.data(isp + 50);
    const auto *isp_51 = buffer.data(isp + 51);
    const auto *isp_52 = buffer.data(isp + 52);
    const auto *isp_53 = buffer.data(isp + 53);
    const auto *isp_54 = buffer.data(isp + 54);
    const auto *isp_55 = buffer.data(isp + 55);
    const auto *isp_56 = buffer.data(isp + 56);
    const auto *isp_58 = buffer.data(isp + 58);
    const auto *isp_59 = buffer.data(isp + 59);
    const auto *isp_60 = buffer.data(isp + 60);
    const auto *isp_62 = buffer.data(isp + 62);

    const auto *isd1_0 = buffer.data(isd1 + 0);
    const auto *isd1_3 = buffer.data(isd1 + 3);
    const auto *isd1_5 = buffer.data(isd1 + 5);
    const auto *isd1_9 = buffer.data(isd1 + 9);
    const auto *isd1_12 = buffer.data(isd1 + 12);
    const auto *isd1_17 = buffer.data(isd1 + 17);
    const auto *isd1_18 = buffer.data(isd1 + 18);
    const auto *isd1_21 = buffer.data(isd1 + 21);
    const auto *isd1_30 = buffer.data(isd1 + 30);
    const auto *isd1_35 = buffer.data(isd1 + 35);
    const auto *isd1_36 = buffer.data(isd1 + 36);
    const auto *isd1_39 = buffer.data(isd1 + 39);
    const auto *isd1_54 = buffer.data(isd1 + 54);
    const auto *isd1_59 = buffer.data(isd1 + 59);
    const auto *isd1_60 = buffer.data(isd1 + 60);
    const auto *isd1_63 = buffer.data(isd1 + 63);
    const auto *isd1_84 = buffer.data(isd1 + 84);
    const auto *isd1_89 = buffer.data(isd1 + 89);

    const auto *kss0_0 = buffer.data(kss0 + 0);
    const auto *kss0_1 = buffer.data(kss0 + 1);
    const auto *kss0_2 = buffer.data(kss0 + 2);
    const auto *kss0_3 = buffer.data(kss0 + 3);
    const auto *kss0_5 = buffer.data(kss0 + 5);
    const auto *kss0_6 = buffer.data(kss0 + 6);
    const auto *kss0_7 = buffer.data(kss0 + 7);
    const auto *kss0_8 = buffer.data(kss0 + 8);
    const auto *kss0_9 = buffer.data(kss0 + 9);
    const auto *kss0_10 = buffer.data(kss0 + 10);
    const auto *kss0_11 = buffer.data(kss0 + 11);
    const auto *kss0_12 = buffer.data(kss0 + 12);
    const auto *kss0_13 = buffer.data(kss0 + 13);
    const auto *kss0_14 = buffer.data(kss0 + 14);
    const auto *kss0_15 = buffer.data(kss0 + 15);
    const auto *kss0_16 = buffer.data(kss0 + 16);
    const auto *kss0_17 = buffer.data(kss0 + 17);
    const auto *kss0_18 = buffer.data(kss0 + 18);
    const auto *kss0_19 = buffer.data(kss0 + 19);
    const auto *kss0_20 = buffer.data(kss0 + 20);

    const auto *kss1_0 = buffer.data(kss1 + 0);
    const auto *kss1_1 = buffer.data(kss1 + 1);
    const auto *kss1_2 = buffer.data(kss1 + 2);
    const auto *kss1_3 = buffer.data(kss1 + 3);
    const auto *kss1_5 = buffer.data(kss1 + 5);
    const auto *kss1_6 = buffer.data(kss1 + 6);
    const auto *kss1_7 = buffer.data(kss1 + 7);
    const auto *kss1_8 = buffer.data(kss1 + 8);
    const auto *kss1_9 = buffer.data(kss1 + 9);
    const auto *kss1_10 = buffer.data(kss1 + 10);
    const auto *kss1_11 = buffer.data(kss1 + 11);
    const auto *kss1_12 = buffer.data(kss1 + 12);
    const auto *kss1_13 = buffer.data(kss1 + 13);
    const auto *kss1_14 = buffer.data(kss1 + 14);
    const auto *kss1_15 = buffer.data(kss1 + 15);
    const auto *kss1_16 = buffer.data(kss1 + 16);
    const auto *kss1_17 = buffer.data(kss1 + 17);
    const auto *kss1_18 = buffer.data(kss1 + 18);
    const auto *kss1_19 = buffer.data(kss1 + 19);
    const auto *kss1_20 = buffer.data(kss1 + 20);

    const auto *ksp_0 = buffer.data(ksp + 0);
    const auto *ksp_1 = buffer.data(ksp + 1);
    const auto *ksp_2 = buffer.data(ksp + 2);
    const auto *ksp_3 = buffer.data(ksp + 3);
    const auto *ksp_4 = buffer.data(ksp + 4);
    const auto *ksp_6 = buffer.data(ksp + 6);
    const auto *ksp_8 = buffer.data(ksp + 8);
    const auto *ksp_9 = buffer.data(ksp + 9);
    const auto *ksp_10 = buffer.data(ksp + 10);
    const auto *ksp_11 = buffer.data(ksp + 11);
    const auto *ksp_13 = buffer.data(ksp + 13);
    const auto *ksp_14 = buffer.data(ksp + 14);
    const auto *ksp_15 = buffer.data(ksp + 15);
    const auto *ksp_16 = buffer.data(ksp + 16);
    const auto *ksp_17 = buffer.data(ksp + 17);
    const auto *ksp_18 = buffer.data(ksp + 18);
    const auto *ksp_19 = buffer.data(ksp + 19);
    const auto *ksp_20 = buffer.data(ksp + 20);
    const auto *ksp_22 = buffer.data(ksp + 22);
    const auto *ksp_23 = buffer.data(ksp + 23);
    const auto *ksp_25 = buffer.data(ksp + 25);
    const auto *ksp_26 = buffer.data(ksp + 26);
    const auto *ksp_27 = buffer.data(ksp + 27);
    const auto *ksp_28 = buffer.data(ksp + 28);
    const auto *ksp_29 = buffer.data(ksp + 29);
    const auto *ksp_30 = buffer.data(ksp + 30);
    const auto *ksp_31 = buffer.data(ksp + 31);
    const auto *ksp_32 = buffer.data(ksp + 32);
    const auto *ksp_34 = buffer.data(ksp + 34);
    const auto *ksp_35 = buffer.data(ksp + 35);
    const auto *ksp_36 = buffer.data(ksp + 36);
    const auto *ksp_37 = buffer.data(ksp + 37);
    const auto *ksp_38 = buffer.data(ksp + 38);
    const auto *ksp_40 = buffer.data(ksp + 40);
    const auto *ksp_41 = buffer.data(ksp + 41);
    const auto *ksp_42 = buffer.data(ksp + 42);
    const auto *ksp_43 = buffer.data(ksp + 43);
    const auto *ksp_44 = buffer.data(ksp + 44);
    const auto *ksp_45 = buffer.data(ksp + 45);
    const auto *ksp_46 = buffer.data(ksp + 46);
    const auto *ksp_47 = buffer.data(ksp + 47);
    const auto *ksp_49 = buffer.data(ksp + 49);
    const auto *ksp_50 = buffer.data(ksp + 50);
    const auto *ksp_51 = buffer.data(ksp + 51);
    const auto *ksp_52 = buffer.data(ksp + 52);
    const auto *ksp_53 = buffer.data(ksp + 53);
    const auto *ksp_54 = buffer.data(ksp + 54);
    const auto *ksp_55 = buffer.data(ksp + 55);
    const auto *ksp_56 = buffer.data(ksp + 56);
    const auto *ksp_58 = buffer.data(ksp + 58);
    const auto *ksp_59 = buffer.data(ksp + 59);
    const auto *ksp_60 = buffer.data(ksp + 60);
    const auto *ksp_61 = buffer.data(ksp + 61);
    const auto *ksp_62 = buffer.data(ksp + 62);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, isp_0, kss0_0, \
                         kss1_0, ksp_0, ksp_1, ksp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * isp_0[k]
                 + f_1 * kss0_0[k]
                 - f_2 * kss1_0[k]
                 + f_3 * pc_x[k] * ksp_0[k];

        t_1[k] = f_3 * pc_y[k] * ksp_0[k];

        t_2[k] = f_3 * pc_z[k] * ksp_0[k];

        t_3[k] = f_1 * kss0_0[k]
                 - f_2 * kss1_0[k]
                 + f_3 * pc_y[k] * ksp_1[k];

        t_4[k] = f_3 * pc_y[k] * ksp_2[k];

        t_5[k] = f_1 * kss0_0[k]
                 - f_2 * kss1_0[k]
                 + f_3 * pc_z[k] * ksp_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pa_y, pc_x, pc_y, pc_z, isd0_0, isp_1, isp_4, \
                         isd1_0, kss0_1, kss1_1, ksp_3, ksp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pa_y[k] * isd0_0[k]
                 - f_4 * pc_y[k] * isd1_0[k];

        t_7[k] = f_5 * isp_4[k]
                 + f_3 * pc_x[k] * ksp_4[k];

        t_8[k] = f_3 * pc_z[k] * ksp_3[k];

        t_9[k] = f_6 * isp_1[k]
                 + f_1 * kss0_1[k]
                 - f_2 * kss1_1[k]
                 + f_3 * pc_y[k] * ksp_4[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_y, pa_z, pc_y, pc_z, isd0_0, isd0_5, \
                         isd1_0, isd1_5, ksp_4, ksp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * pc_z[k] * ksp_4[k];

        t_11[k] = pa_y[k] * isd0_5[k]
                  - f_4 * pc_y[k] * isd1_5[k];

        t_12[k] = pa_z[k] * isd0_0[k]
                  - f_4 * pc_z[k] * isd1_0[k];

        t_13[k] = f_3 * pc_y[k] * ksp_6[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_z, pc_x, pc_y, pc_z, isd0_3, isp_2, isp_8, \
                         isd1_3, kss0_2, kss1_2, ksp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_5 * isp_8[k]
                  + f_3 * pc_x[k] * ksp_8[k];

        t_15[k] = pa_z[k] * isd0_3[k]
                  - f_4 * pc_z[k] * isd1_3[k];

        t_16[k] = f_3 * pc_y[k] * ksp_8[k];

        t_17[k] = f_6 * isp_2[k]
                  + f_1 * kss0_2[k]
                  - f_2 * kss1_2[k]
                  + f_3 * pc_z[k] * ksp_8[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, pc_x, pc_y, pc_z, isp_4, isp_9, isp_10, \
                         kss0_3, kss1_3, ksp_9, ksp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_7 * isp_9[k]
                  + f_1 * kss0_3[k]
                  - f_2 * kss1_3[k]
                  + f_3 * pc_x[k] * ksp_9[k];

        t_19[k] = f_7 * isp_10[k]
                  + f_3 * pc_x[k] * ksp_10[k];

        t_20[k] = f_3 * pc_z[k] * ksp_9[k];

        t_21[k] = f_8 * isp_4[k]
                  + f_1 * kss0_3[k]
                  - f_2 * kss1_3[k]
                  + f_3 * pc_y[k] * ksp_10[k];

        t_22[k] = f_3 * pc_z[k] * ksp_10[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_y, pc_x, pc_y, pc_z, isd0_12, isp_13, isd1_12, \
                         kss0_3, kss1_3, ksp_11, ksp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_1 * kss0_3[k]
                  - f_2 * kss1_3[k]
                  + f_3 * pc_z[k] * ksp_11[k];

        t_24[k] = pa_y[k] * isd0_12[k]
                  - f_4 * pc_y[k] * isd1_12[k];

        t_25[k] = f_7 * isp_13[k]
                  + f_3 * pc_x[k] * ksp_13[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_y, pa_z, pc_x, pc_y, pc_z, isd0_9, \
                         isd0_17, isp_8, isp_14, isd1_9, isd1_17, \
                         ksp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_7 * isp_14[k]
                  + f_3 * pc_x[k] * ksp_14[k];

        t_27[k] = pa_z[k] * isd0_9[k]
                  - f_4 * pc_z[k] * isd1_9[k];

        t_28[k] = f_6 * isp_8[k]
                  + f_3 * pc_y[k] * ksp_14[k];

        t_29[k] = pa_y[k] * isd0_17[k]
                  - f_4 * pc_y[k] * isd1_17[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pc_x, pc_y, isp_15, isp_17, kss0_5, \
                         kss1_5, ksp_15, ksp_16, ksp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_7 * isp_15[k]
                  + f_1 * kss0_5[k]
                  - f_2 * kss1_5[k]
                  + f_3 * pc_x[k] * ksp_15[k];

        t_31[k] = f_3 * pc_y[k] * ksp_15[k];

        t_32[k] = f_7 * isp_17[k]
                  + f_3 * pc_x[k] * ksp_17[k];

        t_33[k] = f_1 * kss0_5[k]
                  - f_2 * kss1_5[k]
                  + f_3 * pc_y[k] * ksp_16[k];

        t_34[k] = f_3 * pc_y[k] * ksp_17[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pc_x, pc_z, isp_8, isp_18, isp_19, kss0_5, \
                         kss0_6, kss1_5, kss1_6, ksp_17, ksp_18, \
                         ksp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_8 * isp_8[k]
                  + f_1 * kss0_5[k]
                  - f_2 * kss1_5[k]
                  + f_3 * pc_z[k] * ksp_17[k];

        t_36[k] = f_9 * isp_18[k]
                  + f_1 * kss0_6[k]
                  - f_2 * kss1_6[k]
                  + f_3 * pc_x[k] * ksp_18[k];

        t_37[k] = f_9 * isp_19[k]
                  + f_3 * pc_x[k] * ksp_19[k];

        t_38[k] = f_3 * pc_z[k] * ksp_18[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pa_z, pc_y, pc_z, isd0_18, isp_10, isd1_18, \
                         kss0_6, kss1_6, ksp_19, ksp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_10 * isp_10[k]
                  + f_1 * kss0_6[k]
                  - f_2 * kss1_6[k]
                  + f_3 * pc_y[k] * ksp_19[k];

        t_40[k] = f_3 * pc_z[k] * ksp_19[k];

        t_41[k] = f_1 * kss0_6[k]
                  - f_2 * kss1_6[k]
                  + f_3 * pc_z[k] * ksp_20[k];

        t_42[k] = pa_z[k] * isd0_18[k]
                  - f_4 * pc_z[k] * isd1_18[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, pa_z, pc_x, pc_y, pc_z, isd0_21, isp_14, \
                         isp_22, isp_23, isd1_21, ksp_22, ksp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_9 * isp_22[k]
                  + f_3 * pc_x[k] * ksp_22[k];

        t_44[k] = f_9 * isp_23[k]
                  + f_3 * pc_x[k] * ksp_23[k];

        t_45[k] = pa_z[k] * isd0_21[k]
                  - f_4 * pc_z[k] * isd1_21[k];

        t_46[k] = f_8 * isp_14[k]
                  + f_3 * pc_y[k] * ksp_23[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, pa_y, pc_x, pc_y, pc_z, isd0_30, isp_11, isp_25, \
                         isd1_30, kss0_7, kss1_7, ksp_23, ksp_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_6 * isp_11[k]
                  + f_1 * kss0_7[k]
                  - f_2 * kss1_7[k]
                  + f_3 * pc_z[k] * ksp_23[k];

        t_48[k] = pa_y[k] * isd0_30[k]
                  - f_4 * pc_y[k] * isd1_30[k];

        t_49[k] = f_9 * isp_25[k]
                  + f_3 * pc_x[k] * ksp_25[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_y, pc_x, pc_y, isd0_35, isp_16, isp_17, \
                         isp_26, isd1_35, kss0_8, kss1_8, ksp_25, \
                         ksp_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_9 * isp_26[k]
                  + f_3 * pc_x[k] * ksp_26[k];

        t_51[k] = f_6 * isp_16[k]
                  + f_1 * kss0_8[k]
                  - f_2 * kss1_8[k]
                  + f_3 * pc_y[k] * ksp_25[k];

        t_52[k] = f_6 * isp_17[k]
                  + f_3 * pc_y[k] * ksp_26[k];

        t_53[k] = pa_y[k] * isd0_35[k]
                  - f_4 * pc_y[k] * isd1_35[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, pc_x, pc_y, isp_27, isp_29, kss0_9, \
                         kss1_9, ksp_27, ksp_28, ksp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_9 * isp_27[k]
                  + f_1 * kss0_9[k]
                  - f_2 * kss1_9[k]
                  + f_3 * pc_x[k] * ksp_27[k];

        t_55[k] = f_3 * pc_y[k] * ksp_27[k];

        t_56[k] = f_9 * isp_29[k]
                  + f_3 * pc_x[k] * ksp_29[k];

        t_57[k] = f_1 * kss0_9[k]
                  - f_2 * kss1_9[k]
                  + f_3 * pc_y[k] * ksp_28[k];

        t_58[k] = f_3 * pc_y[k] * ksp_29[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pc_x, pc_z, isp_17, isp_30, isp_31, kss0_9, \
                         kss0_10, kss1_9, kss1_10, ksp_29, ksp_30, \
                         ksp_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_10 * isp_17[k]
                  + f_1 * kss0_9[k]
                  - f_2 * kss1_9[k]
                  + f_3 * pc_z[k] * ksp_29[k];

        t_60[k] = f_10 * isp_30[k]
                  + f_1 * kss0_10[k]
                  - f_2 * kss1_10[k]
                  + f_3 * pc_x[k] * ksp_30[k];

        t_61[k] = f_10 * isp_31[k]
                  + f_3 * pc_x[k] * ksp_31[k];

        t_62[k] = f_3 * pc_z[k] * ksp_30[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pa_z, pc_y, pc_z, isd0_36, isp_19, isd1_36, \
                         kss0_10, kss1_10, ksp_31, ksp_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_9 * isp_19[k]
                  + f_1 * kss0_10[k]
                  - f_2 * kss1_10[k]
                  + f_3 * pc_y[k] * ksp_31[k];

        t_64[k] = f_3 * pc_z[k] * ksp_31[k];

        t_65[k] = f_1 * kss0_10[k]
                  - f_2 * kss1_10[k]
                  + f_3 * pc_z[k] * ksp_32[k];

        t_66[k] = pa_z[k] * isd0_36[k]
                  - f_4 * pc_z[k] * isd1_36[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, pa_z, pc_x, pc_y, pc_z, isd0_39, isp_23, \
                         isp_34, isp_35, isd1_39, ksp_34, ksp_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_10 * isp_34[k]
                  + f_3 * pc_x[k] * ksp_34[k];

        t_68[k] = f_10 * isp_35[k]
                  + f_3 * pc_x[k] * ksp_35[k];

        t_69[k] = pa_z[k] * isd0_39[k]
                  - f_4 * pc_z[k] * isd1_39[k];

        t_70[k] = f_10 * isp_23[k]
                  + f_3 * pc_y[k] * ksp_35[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, pc_x, pc_z, isp_20, isp_36, isp_37, kss0_11, \
                         kss0_12, kss1_11, kss1_12, ksp_35, ksp_36, \
                         ksp_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_6 * isp_20[k]
                  + f_1 * kss0_11[k]
                  - f_2 * kss1_11[k]
                  + f_3 * pc_z[k] * ksp_35[k];

        t_72[k] = f_10 * isp_36[k]
                  + f_1 * kss0_12[k]
                  - f_2 * kss1_12[k]
                  + f_3 * pc_x[k] * ksp_36[k];

        t_73[k] = f_10 * isp_37[k]
                  + f_3 * pc_x[k] * ksp_37[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pc_x, pc_y, pc_z, isp_23, isp_25, isp_26, \
                         isp_38, kss0_12, kss1_12, ksp_37, ksp_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_10 * isp_38[k]
                  + f_3 * pc_x[k] * ksp_38[k];

        t_75[k] = f_8 * isp_25[k]
                  + f_1 * kss0_12[k]
                  - f_2 * kss1_12[k]
                  + f_3 * pc_y[k] * ksp_37[k];

        t_76[k] = f_8 * isp_26[k]
                  + f_3 * pc_y[k] * ksp_38[k];

        t_77[k] = f_8 * isp_23[k]
                  + f_1 * kss0_12[k]
                  - f_2 * kss1_12[k]
                  + f_3 * pc_z[k] * ksp_38[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pa_y, pc_x, pc_y, isd0_54, isp_28, isp_40, \
                         isp_41, isd1_54, kss0_13, kss1_13, ksp_40, \
                         ksp_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = pa_y[k] * isd0_54[k]
                  - f_4 * pc_y[k] * isd1_54[k];

        t_79[k] = f_10 * isp_40[k]
                  + f_3 * pc_x[k] * ksp_40[k];

        t_80[k] = f_10 * isp_41[k]
                  + f_3 * pc_x[k] * ksp_41[k];

        t_81[k] = f_6 * isp_28[k]
                  + f_1 * kss0_13[k]
                  - f_2 * kss1_13[k]
                  + f_3 * pc_y[k] * ksp_40[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, pa_y, pc_x, pc_y, isd0_59, isp_29, isp_42, \
                         isd1_59, kss0_14, kss1_14, ksp_41, ksp_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_6 * isp_29[k]
                  + f_3 * pc_y[k] * ksp_41[k];

        t_83[k] = pa_y[k] * isd0_59[k]
                  - f_4 * pc_y[k] * isd1_59[k];

        t_84[k] = f_10 * isp_42[k]
                  + f_1 * kss0_14[k]
                  - f_2 * kss1_14[k]
                  + f_3 * pc_x[k] * ksp_42[k];

        t_85[k] = f_3 * pc_y[k] * ksp_42[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pc_x, pc_y, pc_z, isp_29, isp_44, kss0_14, \
                         kss1_14, ksp_43, ksp_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_10 * isp_44[k]
                  + f_3 * pc_x[k] * ksp_44[k];

        t_87[k] = f_1 * kss0_14[k]
                  - f_2 * kss1_14[k]
                  + f_3 * pc_y[k] * ksp_43[k];

        t_88[k] = f_3 * pc_y[k] * ksp_44[k];

        t_89[k] = f_9 * isp_29[k]
                  + f_1 * kss0_14[k]
                  - f_2 * kss1_14[k]
                  + f_3 * pc_z[k] * ksp_44[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, pc_x, pc_y, pc_z, isp_31, isp_45, \
                         isp_46, kss0_15, kss1_15, ksp_45, ksp_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_8 * isp_45[k]
                  + f_1 * kss0_15[k]
                  - f_2 * kss1_15[k]
                  + f_3 * pc_x[k] * ksp_45[k];

        t_91[k] = f_8 * isp_46[k]
                  + f_3 * pc_x[k] * ksp_46[k];

        t_92[k] = f_3 * pc_z[k] * ksp_45[k];

        t_93[k] = f_7 * isp_31[k]
                  + f_1 * kss0_15[k]
                  - f_2 * kss1_15[k]
                  + f_3 * pc_y[k] * ksp_46[k];

        t_94[k] = f_3 * pc_z[k] * ksp_46[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pa_z, pc_x, pc_z, isd0_60, isp_49, isp_50, \
                         isd1_60, kss0_15, kss1_15, ksp_47, ksp_49, \
                         ksp_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_1 * kss0_15[k]
                  - f_2 * kss1_15[k]
                  + f_3 * pc_z[k] * ksp_47[k];

        t_96[k] = pa_z[k] * isd0_60[k]
                  - f_4 * pc_z[k] * isd1_60[k];

        t_97[k] = f_8 * isp_49[k]
                  + f_3 * pc_x[k] * ksp_49[k];

        t_98[k] = f_8 * isp_50[k]
                  + f_3 * pc_x[k] * ksp_50[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, pa_z, pc_y, pc_z, isd0_63, isp_32, isp_35, \
                         isd1_63, kss0_16, kss1_16, ksp_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = pa_z[k] * isd0_63[k]
                  - f_4 * pc_z[k] * isd1_63[k];

        t_100[k] = f_9 * isp_35[k]
                   + f_3 * pc_y[k] * ksp_50[k];

        t_101[k] = f_6 * isp_32[k]
                   + f_1 * kss0_16[k]
                   - f_2 * kss1_16[k]
                   + f_3 * pc_z[k] * ksp_50[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, pc_x, pc_y, isp_37, isp_51, isp_52, \
                         isp_53, kss0_17, kss1_17, ksp_51, ksp_52, \
                         ksp_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_8 * isp_51[k]
                   + f_1 * kss0_17[k]
                   - f_2 * kss1_17[k]
                   + f_3 * pc_x[k] * ksp_51[k];

        t_103[k] = f_8 * isp_52[k]
                   + f_3 * pc_x[k] * ksp_52[k];

        t_104[k] = f_8 * isp_53[k]
                   + f_3 * pc_x[k] * ksp_53[k];

        t_105[k] = f_10 * isp_37[k]
                   + f_1 * kss0_17[k]
                   - f_2 * kss1_17[k]
                   + f_3 * pc_y[k] * ksp_52[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, pc_x, pc_y, pc_z, isp_35, isp_38, isp_54, \
                         kss0_17, kss0_18, kss1_17, kss1_18, ksp_53, \
                         ksp_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_10 * isp_38[k]
                   + f_3 * pc_y[k] * ksp_53[k];

        t_107[k] = f_8 * isp_35[k]
                   + f_1 * kss0_17[k]
                   - f_2 * kss1_17[k]
                   + f_3 * pc_z[k] * ksp_53[k];

        t_108[k] = f_8 * isp_54[k]
                   + f_1 * kss0_18[k]
                   - f_2 * kss1_18[k]
                   + f_3 * pc_x[k] * ksp_54[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pc_x, pc_y, isp_40, isp_41, isp_55, \
                         isp_56, kss0_18, kss1_18, ksp_55, ksp_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_8 * isp_55[k]
                   + f_3 * pc_x[k] * ksp_55[k];

        t_110[k] = f_8 * isp_56[k]
                   + f_3 * pc_x[k] * ksp_56[k];

        t_111[k] = f_8 * isp_40[k]
                   + f_1 * kss0_18[k]
                   - f_2 * kss1_18[k]
                   + f_3 * pc_y[k] * ksp_55[k];

        t_112[k] = f_8 * isp_41[k]
                   + f_3 * pc_y[k] * ksp_56[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, pa_y, pc_x, pc_y, pc_z, isd0_84, isp_38, isp_58, \
                         isd1_84, kss0_18, kss1_18, ksp_56, ksp_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_10 * isp_38[k]
                   + f_1 * kss0_18[k]
                   - f_2 * kss1_18[k]
                   + f_3 * pc_z[k] * ksp_56[k];

        t_114[k] = pa_y[k] * isd0_84[k]
                   - f_4 * pc_y[k] * isd1_84[k];

        t_115[k] = f_8 * isp_58[k]
                   + f_3 * pc_x[k] * ksp_58[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, pa_y, pc_x, pc_y, isd0_89, isp_43, \
                         isp_44, isp_59, isd1_89, kss0_19, kss1_19, ksp_58, \
                         ksp_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_8 * isp_59[k]
                   + f_3 * pc_x[k] * ksp_59[k];

        t_117[k] = f_6 * isp_43[k]
                   + f_1 * kss0_19[k]
                   - f_2 * kss1_19[k]
                   + f_3 * pc_y[k] * ksp_58[k];

        t_118[k] = f_6 * isp_44[k]
                   + f_3 * pc_y[k] * ksp_59[k];

        t_119[k] = pa_y[k] * isd0_89[k]
                   - f_4 * pc_y[k] * isd1_89[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, pc_x, pc_y, isp_60, isp_62, \
                         kss0_20, kss1_20, ksp_60, ksp_61, ksp_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_8 * isp_60[k]
                   + f_1 * kss0_20[k]
                   - f_2 * kss1_20[k]
                   + f_3 * pc_x[k] * ksp_60[k];

        t_121[k] = f_3 * pc_y[k] * ksp_60[k];

        t_122[k] = f_8 * isp_62[k]
                   + f_3 * pc_x[k] * ksp_62[k];

        t_123[k] = f_1 * kss0_20[k]
                   - f_2 * kss1_20[k]
                   + f_3 * pc_y[k] * ksp_61[k];

        t_124[k] = f_3 * pc_y[k] * ksp_62[k];
    }
}

static auto
compute_prim_ksd_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t isd0,
                                                          const size_t isp, const size_t isd1,
                                                          const size_t kss0, const size_t kss1,
                                                          const size_t ksp, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / q;
    const auto f_1 = 0.5 / gamma;
    const auto f_2 = 0.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 3.0 / q;
    const auto f_6 = 0.5 / q;
    const auto f_7 = 2.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 2.0 / q;
    const auto f_10 = 1.5 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *isd0_90 = buffer.data(isd0 + 90);
    const auto *isd0_120 = buffer.data(isd0 + 120);
    const auto *isd0_126 = buffer.data(isd0 + 126);
    const auto *isd0_129 = buffer.data(isd0 + 129);
    const auto *isd0_131 = buffer.data(isd0 + 131);
    const auto *isd0_135 = buffer.data(isd0 + 135);
    const auto *isd0_137 = buffer.data(isd0 + 137);
    const auto *isd0_138 = buffer.data(isd0 + 138);
    const auto *isd0_141 = buffer.data(isd0 + 141);
    const auto *isd0_143 = buffer.data(isd0 + 143);
    const auto *isd0_144 = buffer.data(isd0 + 144);
    const auto *isd0_147 = buffer.data(isd0 + 147);
    const auto *isd0_149 = buffer.data(isd0 + 149);
    const auto *isd0_150 = buffer.data(isd0 + 150);
    const auto *isd0_153 = buffer.data(isd0 + 153);
    const auto *isd0_155 = buffer.data(isd0 + 155);
    const auto *isd0_159 = buffer.data(isd0 + 159);
    const auto *isd0_161 = buffer.data(isd0 + 161);
    const auto *isd0_162 = buffer.data(isd0 + 162);
    const auto *isd0_165 = buffer.data(isd0 + 165);
    const auto *isd0_167 = buffer.data(isd0 + 167);

    const auto *isp_44 = buffer.data(isp + 44);
    const auto *isp_50 = buffer.data(isp + 50);
    const auto *isp_53 = buffer.data(isp + 53);
    const auto *isp_56 = buffer.data(isp + 56);
    const auto *isp_59 = buffer.data(isp + 59);
    const auto *isp_62 = buffer.data(isp + 62);
    const auto *isp_63 = buffer.data(isp + 63);
    const auto *isp_64 = buffer.data(isp + 64);
    const auto *isp_65 = buffer.data(isp + 65);
    const auto *isp_67 = buffer.data(isp + 67);
    const auto *isp_68 = buffer.data(isp + 68);
    const auto *isp_69 = buffer.data(isp + 69);
    const auto *isp_70 = buffer.data(isp + 70);
    const auto *isp_71 = buffer.data(isp + 71);
    const auto *isp_72 = buffer.data(isp + 72);
    const auto *isp_73 = buffer.data(isp + 73);
    const auto *isp_74 = buffer.data(isp + 74);
    const auto *isp_75 = buffer.data(isp + 75);
    const auto *isp_76 = buffer.data(isp + 76);
    const auto *isp_77 = buffer.data(isp + 77);
    const auto *isp_79 = buffer.data(isp + 79);
    const auto *isp_80 = buffer.data(isp + 80);
    const auto *isp_81 = buffer.data(isp + 81);
    const auto *isp_82 = buffer.data(isp + 82);
    const auto *isp_83 = buffer.data(isp + 83);

    const auto *isd1_90 = buffer.data(isd1 + 90);
    const auto *isd1_120 = buffer.data(isd1 + 120);
    const auto *isd1_126 = buffer.data(isd1 + 126);
    const auto *isd1_129 = buffer.data(isd1 + 129);
    const auto *isd1_131 = buffer.data(isd1 + 131);
    const auto *isd1_135 = buffer.data(isd1 + 135);
    const auto *isd1_137 = buffer.data(isd1 + 137);
    const auto *isd1_138 = buffer.data(isd1 + 138);
    const auto *isd1_141 = buffer.data(isd1 + 141);
    const auto *isd1_143 = buffer.data(isd1 + 143);
    const auto *isd1_144 = buffer.data(isd1 + 144);
    const auto *isd1_147 = buffer.data(isd1 + 147);
    const auto *isd1_149 = buffer.data(isd1 + 149);
    const auto *isd1_150 = buffer.data(isd1 + 150);
    const auto *isd1_153 = buffer.data(isd1 + 153);
    const auto *isd1_155 = buffer.data(isd1 + 155);
    const auto *isd1_159 = buffer.data(isd1 + 159);
    const auto *isd1_161 = buffer.data(isd1 + 161);
    const auto *isd1_162 = buffer.data(isd1 + 162);
    const auto *isd1_165 = buffer.data(isd1 + 165);
    const auto *isd1_167 = buffer.data(isd1 + 167);

    const auto *kss0_20 = buffer.data(kss0 + 20);
    const auto *kss0_28 = buffer.data(kss0 + 28);
    const auto *kss0_29 = buffer.data(kss0 + 29);
    const auto *kss0_30 = buffer.data(kss0 + 30);
    const auto *kss0_31 = buffer.data(kss0 + 31);
    const auto *kss0_32 = buffer.data(kss0 + 32);
    const auto *kss0_33 = buffer.data(kss0 + 33);
    const auto *kss0_35 = buffer.data(kss0 + 35);

    const auto *kss1_20 = buffer.data(kss1 + 20);
    const auto *kss1_28 = buffer.data(kss1 + 28);
    const auto *kss1_29 = buffer.data(kss1 + 29);
    const auto *kss1_30 = buffer.data(kss1 + 30);
    const auto *kss1_31 = buffer.data(kss1 + 31);
    const auto *kss1_32 = buffer.data(kss1 + 32);
    const auto *kss1_33 = buffer.data(kss1 + 33);
    const auto *kss1_35 = buffer.data(kss1 + 35);

    const auto *ksp_62 = buffer.data(ksp + 62);
    const auto *ksp_63 = buffer.data(ksp + 63);
    const auto *ksp_64 = buffer.data(ksp + 64);
    const auto *ksp_67 = buffer.data(ksp + 67);
    const auto *ksp_68 = buffer.data(ksp + 68);
    const auto *ksp_70 = buffer.data(ksp + 70);
    const auto *ksp_71 = buffer.data(ksp + 71);
    const auto *ksp_73 = buffer.data(ksp + 73);
    const auto *ksp_74 = buffer.data(ksp + 74);
    const auto *ksp_76 = buffer.data(ksp + 76);
    const auto *ksp_77 = buffer.data(ksp + 77);
    const auto *ksp_79 = buffer.data(ksp + 79);
    const auto *ksp_80 = buffer.data(ksp + 80);
    const auto *ksp_81 = buffer.data(ksp + 81);
    const auto *ksp_83 = buffer.data(ksp + 83);
    const auto *ksp_84 = buffer.data(ksp + 84);
    const auto *ksp_85 = buffer.data(ksp + 85);
    const auto *ksp_86 = buffer.data(ksp + 86);
    const auto *ksp_88 = buffer.data(ksp + 88);
    const auto *ksp_89 = buffer.data(ksp + 89);
    const auto *ksp_90 = buffer.data(ksp + 90);
    const auto *ksp_91 = buffer.data(ksp + 91);
    const auto *ksp_92 = buffer.data(ksp + 92);
    const auto *ksp_93 = buffer.data(ksp + 93);
    const auto *ksp_94 = buffer.data(ksp + 94);
    const auto *ksp_95 = buffer.data(ksp + 95);
    const auto *ksp_96 = buffer.data(ksp + 96);
    const auto *ksp_97 = buffer.data(ksp + 97);
    const auto *ksp_98 = buffer.data(ksp + 98);
    const auto *ksp_99 = buffer.data(ksp + 99);
    const auto *ksp_100 = buffer.data(ksp + 100);
    const auto *ksp_101 = buffer.data(ksp + 101);
    const auto *ksp_103 = buffer.data(ksp + 103);
    const auto *ksp_104 = buffer.data(ksp + 104);
    const auto *ksp_105 = buffer.data(ksp + 105);
    const auto *ksp_106 = buffer.data(ksp + 106);
    const auto *ksp_107 = buffer.data(ksp + 107);

#pragma omp simd aligned(t_125, t_126, t_127, pa_x, pc_x, pc_z, isd0_126, isp_44, isp_63, \
                         isp_64, isd1_126, kss0_20, kss1_20, ksp_62, \
                         ksp_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_7 * isp_44[k]
                   + f_1 * kss0_20[k]
                   - f_2 * kss1_20[k]
                   + f_3 * pc_z[k] * ksp_62[k];

        t_126[k] = pa_x[k] * isd0_126[k]
                   + f_8 * isp_63[k]
                   - f_4 * pc_x[k] * isd1_126[k];

        t_127[k] = f_6 * isp_64[k]
                   + f_3 * pc_x[k] * ksp_64[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, pa_x, pc_x, pc_z, isd0_129, isd0_131, \
                         isd1_129, isd1_131, ksp_63, ksp_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_3 * pc_z[k] * ksp_63[k];

        t_129[k] = pa_x[k] * isd0_129[k]
                   - f_4 * pc_x[k] * isd1_129[k];

        t_130[k] = f_3 * pc_z[k] * ksp_64[k];

        t_131[k] = pa_x[k] * isd0_131[k]
                   - f_4 * pc_x[k] * isd1_131[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pa_x, pa_z, pc_x, pc_z, isd0_90, \
                         isd0_135, isp_67, isp_68, isd1_90, isd1_135, ksp_67, \
                         ksp_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = pa_z[k] * isd0_90[k]
                   - f_4 * pc_z[k] * isd1_90[k];

        t_133[k] = f_6 * isp_67[k]
                   + f_3 * pc_x[k] * ksp_67[k];

        t_134[k] = f_6 * isp_68[k]
                   + f_3 * pc_x[k] * ksp_68[k];

        t_135[k] = pa_x[k] * isd0_135[k]
                   - f_4 * pc_x[k] * isd1_135[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, pa_x, pc_x, pc_y, isd0_137, isd0_138, \
                         isp_50, isp_69, isp_70, isd1_137, isd1_138, ksp_68, \
                         ksp_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_7 * isp_50[k]
                   + f_3 * pc_y[k] * ksp_68[k];

        t_137[k] = pa_x[k] * isd0_137[k]
                   - f_4 * pc_x[k] * isd1_137[k];

        t_138[k] = pa_x[k] * isd0_138[k]
                   + f_8 * isp_69[k]
                   - f_4 * pc_x[k] * isd1_138[k];

        t_139[k] = f_6 * isp_70[k]
                   + f_3 * pc_x[k] * ksp_70[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, pa_x, pc_x, pc_y, isd0_141, isd0_143, \
                         isp_53, isp_71, isd1_141, isd1_143, ksp_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_6 * isp_71[k]
                   + f_3 * pc_x[k] * ksp_71[k];

        t_141[k] = pa_x[k] * isd0_141[k]
                   - f_4 * pc_x[k] * isd1_141[k];

        t_142[k] = f_9 * isp_53[k]
                   + f_3 * pc_y[k] * ksp_71[k];

        t_143[k] = pa_x[k] * isd0_143[k]
                   - f_4 * pc_x[k] * isd1_143[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, pa_x, pc_x, isd0_144, isd0_147, isp_72, \
                         isp_73, isp_74, isd1_144, isd1_147, ksp_73, \
                         ksp_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = pa_x[k] * isd0_144[k]
                   + f_8 * isp_72[k]
                   - f_4 * pc_x[k] * isd1_144[k];

        t_145[k] = f_6 * isp_73[k]
                   + f_3 * pc_x[k] * ksp_73[k];

        t_146[k] = f_6 * isp_74[k]
                   + f_3 * pc_x[k] * ksp_74[k];

        t_147[k] = pa_x[k] * isd0_147[k]
                   - f_4 * pc_x[k] * isd1_147[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, t_151, pa_x, pc_x, pc_y, isd0_149, isd0_150, \
                         isp_56, isp_75, isp_76, isd1_149, isd1_150, ksp_74, \
                         ksp_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_10 * isp_56[k]
                   + f_3 * pc_y[k] * ksp_74[k];

        t_149[k] = pa_x[k] * isd0_149[k]
                   - f_4 * pc_x[k] * isd1_149[k];

        t_150[k] = pa_x[k] * isd0_150[k]
                   + f_8 * isp_75[k]
                   - f_4 * pc_x[k] * isd1_150[k];

        t_151[k] = f_6 * isp_76[k]
                   + f_3 * pc_x[k] * ksp_76[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, pa_x, pc_x, pc_y, isd0_153, isd0_155, \
                         isp_59, isp_77, isd1_153, isd1_155, ksp_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_6 * isp_77[k]
                   + f_3 * pc_x[k] * ksp_77[k];

        t_153[k] = pa_x[k] * isd0_153[k]
                   - f_4 * pc_x[k] * isd1_153[k];

        t_154[k] = f_8 * isp_59[k]
                   + f_3 * pc_y[k] * ksp_77[k];

        t_155[k] = pa_x[k] * isd0_155[k]
                   - f_4 * pc_x[k] * isd1_155[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, pa_x, pa_y, pc_x, pc_y, isd0_120, \
                         isd0_159, isp_79, isp_80, isd1_120, isd1_159, ksp_79, \
                         ksp_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = pa_y[k] * isd0_120[k]
                   - f_4 * pc_y[k] * isd1_120[k];

        t_157[k] = f_6 * isp_79[k]
                   + f_3 * pc_x[k] * ksp_79[k];

        t_158[k] = f_6 * isp_80[k]
                   + f_3 * pc_x[k] * ksp_80[k];

        t_159[k] = pa_x[k] * isd0_159[k]
                   - f_4 * pc_x[k] * isd1_159[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, pa_x, pc_x, pc_y, isd0_161, isd0_162, \
                         isp_62, isp_81, isd1_161, isd1_162, ksp_80, \
                         ksp_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = f_6 * isp_62[k]
                   + f_3 * pc_y[k] * ksp_80[k];

        t_161[k] = pa_x[k] * isd0_161[k]
                   - f_4 * pc_x[k] * isd1_161[k];

        t_162[k] = pa_x[k] * isd0_162[k]
                   + f_8 * isp_81[k]
                   - f_4 * pc_x[k] * isd1_162[k];

        t_163[k] = f_3 * pc_y[k] * ksp_81[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, pa_x, pc_x, pc_y, isd0_165, isd0_167, \
                         isp_83, isd1_165, isd1_167, ksp_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_6 * isp_83[k]
                   + f_3 * pc_x[k] * ksp_83[k];

        t_165[k] = pa_x[k] * isd0_165[k]
                   - f_4 * pc_x[k] * isd1_165[k];

        t_166[k] = f_3 * pc_y[k] * ksp_83[k];

        t_167[k] = pa_x[k] * isd0_167[k]
                   - f_4 * pc_x[k] * isd1_167[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, t_172, t_173, pc_x, pc_y, pc_z, isp_64, \
                         kss0_28, kss1_28, ksp_84, ksp_85, ksp_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_1 * kss0_28[k]
                   - f_2 * kss1_28[k]
                   + f_3 * pc_x[k] * ksp_84[k];

        t_169[k] = f_3 * pc_x[k] * ksp_85[k];

        t_170[k] = f_3 * pc_x[k] * ksp_86[k];

        t_171[k] = f_0 * isp_64[k]
                   + f_1 * kss0_28[k]
                   - f_2 * kss1_28[k]
                   + f_3 * pc_y[k] * ksp_85[k];

        t_172[k] = f_3 * pc_z[k] * ksp_85[k];

        t_173[k] = f_1 * kss0_28[k]
                   - f_2 * kss1_28[k]
                   + f_3 * pc_z[k] * ksp_86[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, t_178, pa_z, pc_x, pc_y, pc_z, isd0_126, \
                         isd0_129, isp_68, isd1_126, isd1_129, ksp_88, \
                         ksp_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = pa_z[k] * isd0_126[k]
                   - f_4 * pc_z[k] * isd1_126[k];

        t_175[k] = f_3 * pc_x[k] * ksp_88[k];

        t_176[k] = f_3 * pc_x[k] * ksp_89[k];

        t_177[k] = pa_z[k] * isd0_129[k]
                   - f_4 * pc_z[k] * isd1_129[k];

        t_178[k] = f_5 * isp_68[k]
                   + f_3 * pc_y[k] * ksp_89[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, t_182, pc_x, pc_z, isp_65, kss0_29, kss0_30, \
                         kss1_29, kss1_30, ksp_89, ksp_90, ksp_91, \
                         ksp_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = f_6 * isp_65[k]
                   + f_1 * kss0_29[k]
                   - f_2 * kss1_29[k]
                   + f_3 * pc_z[k] * ksp_89[k];

        t_180[k] = f_1 * kss0_30[k]
                   - f_2 * kss1_30[k]
                   + f_3 * pc_x[k] * ksp_90[k];

        t_181[k] = f_3 * pc_x[k] * ksp_91[k];

        t_182[k] = f_3 * pc_x[k] * ksp_92[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, pc_y, pc_z, isp_68, isp_70, isp_71, kss0_30, \
                         kss1_30, ksp_91, ksp_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_7 * isp_70[k]
                   + f_1 * kss0_30[k]
                   - f_2 * kss1_30[k]
                   + f_3 * pc_y[k] * ksp_91[k];

        t_184[k] = f_7 * isp_71[k]
                   + f_3 * pc_y[k] * ksp_92[k];

        t_185[k] = f_8 * isp_68[k]
                   + f_1 * kss0_30[k]
                   - f_2 * kss1_30[k]
                   + f_3 * pc_z[k] * ksp_92[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, t_189, t_190, pc_x, pc_y, isp_73, isp_74, \
                         kss0_31, kss1_31, ksp_93, ksp_94, ksp_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_1 * kss0_31[k]
                   - f_2 * kss1_31[k]
                   + f_3 * pc_x[k] * ksp_93[k];

        t_187[k] = f_3 * pc_x[k] * ksp_94[k];

        t_188[k] = f_3 * pc_x[k] * ksp_95[k];

        t_189[k] = f_9 * isp_73[k]
                   + f_1 * kss0_31[k]
                   - f_2 * kss1_31[k]
                   + f_3 * pc_y[k] * ksp_94[k];

        t_190[k] = f_9 * isp_74[k]
                   + f_3 * pc_y[k] * ksp_95[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, t_194, pc_x, pc_z, isp_71, kss0_31, kss0_32, \
                         kss1_31, kss1_32, ksp_95, ksp_96, ksp_97, \
                         ksp_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = f_10 * isp_71[k]
                   + f_1 * kss0_31[k]
                   - f_2 * kss1_31[k]
                   + f_3 * pc_z[k] * ksp_95[k];

        t_192[k] = f_1 * kss0_32[k]
                   - f_2 * kss1_32[k]
                   + f_3 * pc_x[k] * ksp_96[k];

        t_193[k] = f_3 * pc_x[k] * ksp_97[k];

        t_194[k] = f_3 * pc_x[k] * ksp_98[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, pc_y, pc_z, isp_74, isp_76, isp_77, kss0_32, \
                         kss1_32, ksp_97, ksp_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = f_10 * isp_76[k]
                   + f_1 * kss0_32[k]
                   - f_2 * kss1_32[k]
                   + f_3 * pc_y[k] * ksp_97[k];

        t_196[k] = f_10 * isp_77[k]
                   + f_3 * pc_y[k] * ksp_98[k];

        t_197[k] = f_9 * isp_74[k]
                   + f_1 * kss0_32[k]
                   - f_2 * kss1_32[k]
                   + f_3 * pc_z[k] * ksp_98[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, t_202, pc_x, pc_y, isp_79, isp_80, \
                         kss0_33, kss1_33, ksp_99, ksp_100, ksp_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_1 * kss0_33[k]
                   - f_2 * kss1_33[k]
                   + f_3 * pc_x[k] * ksp_99[k];

        t_199[k] = f_3 * pc_x[k] * ksp_100[k];

        t_200[k] = f_3 * pc_x[k] * ksp_101[k];

        t_201[k] = f_8 * isp_79[k]
                   + f_1 * kss0_33[k]
                   - f_2 * kss1_33[k]
                   + f_3 * pc_y[k] * ksp_100[k];

        t_202[k] = f_8 * isp_80[k]
                   + f_3 * pc_y[k] * ksp_101[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, pa_y, pc_x, pc_y, pc_z, isd0_162, isp_77, \
                         isd1_162, kss0_33, kss1_33, ksp_101, ksp_103, \
                         ksp_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_7 * isp_77[k]
                   + f_1 * kss0_33[k]
                   - f_2 * kss1_33[k]
                   + f_3 * pc_z[k] * ksp_101[k];

        t_204[k] = pa_y[k] * isd0_162[k]
                   - f_4 * pc_y[k] * isd1_162[k];

        t_205[k] = f_3 * pc_x[k] * ksp_103[k];

        t_206[k] = f_3 * pc_x[k] * ksp_104[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, pa_y, pc_y, isd0_165, isd0_167, isp_82, isp_83, \
                         isd1_165, isd1_167, ksp_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = pa_y[k] * isd0_165[k]
                   + f_8 * isp_82[k]
                   - f_4 * pc_y[k] * isd1_165[k];

        t_208[k] = f_6 * isp_83[k]
                   + f_3 * pc_y[k] * ksp_104[k];

        t_209[k] = pa_y[k] * isd0_167[k]
                   - f_4 * pc_y[k] * isd1_167[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, t_215, pc_x, pc_y, pc_z, isp_83, \
                         kss0_35, kss1_35, ksp_105, ksp_106, ksp_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_1 * kss0_35[k]
                   - f_2 * kss1_35[k]
                   + f_3 * pc_x[k] * ksp_105[k];

        t_211[k] = f_3 * pc_x[k] * ksp_106[k];

        t_212[k] = f_3 * pc_x[k] * ksp_107[k];

        t_213[k] = f_1 * kss0_35[k]
                   - f_2 * kss1_35[k]
                   + f_3 * pc_y[k] * ksp_106[k];

        t_214[k] = f_3 * pc_y[k] * ksp_107[k];

        t_215[k] = f_0 * isp_83[k]
                   + f_1 * kss0_35[k]
                   - f_2 * kss1_35[k]
                   + f_3 * pc_z[k] * ksp_107[k];
    }
}

auto
compute_prim_ksd_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t isd0, const size_t isp,
                                                   const size_t isd1, const size_t kss0,
                                                   const size_t kss1, const size_t ksp,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_ksd_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, isd0, isp,
                                                              isd1, kss0, kss1, ksp, ncols,
                                                              gamma, p, q);

    compute_prim_ksd_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, isd0, isp,
                                                              isd1, kss0, kss1, ksp, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
