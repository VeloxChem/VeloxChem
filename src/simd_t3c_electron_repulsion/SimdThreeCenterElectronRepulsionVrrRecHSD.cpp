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


#include "SimdThreeCenterElectronRepulsionVrrRecHSD.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_hsd_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t gsd0, const size_t gsp,
                                                   const size_t gsd1, const size_t hss0,
                                                   const size_t hss1, const size_t hsp,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
    const auto f_1 = 0.5 / gamma;
    const auto f_2 = 0.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 2.0 / q;
    const auto f_6 = 0.5 / q;
    const auto f_7 = 1.5 / q;
    const auto f_8 = 1.0 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *gsd0_0 = buffer.data(gsd0 + 0);
    const auto *gsd0_3 = buffer.data(gsd0 + 3);
    const auto *gsd0_5 = buffer.data(gsd0 + 5);
    const auto *gsd0_9 = buffer.data(gsd0 + 9);
    const auto *gsd0_12 = buffer.data(gsd0 + 12);
    const auto *gsd0_17 = buffer.data(gsd0 + 17);
    const auto *gsd0_18 = buffer.data(gsd0 + 18);
    const auto *gsd0_21 = buffer.data(gsd0 + 21);
    const auto *gsd0_30 = buffer.data(gsd0 + 30);
    const auto *gsd0_35 = buffer.data(gsd0 + 35);
    const auto *gsd0_36 = buffer.data(gsd0 + 36);
    const auto *gsd0_54 = buffer.data(gsd0 + 54);
    const auto *gsd0_60 = buffer.data(gsd0 + 60);
    const auto *gsd0_63 = buffer.data(gsd0 + 63);
    const auto *gsd0_65 = buffer.data(gsd0 + 65);
    const auto *gsd0_69 = buffer.data(gsd0 + 69);
    const auto *gsd0_71 = buffer.data(gsd0 + 71);
    const auto *gsd0_72 = buffer.data(gsd0 + 72);
    const auto *gsd0_75 = buffer.data(gsd0 + 75);
    const auto *gsd0_77 = buffer.data(gsd0 + 77);
    const auto *gsd0_81 = buffer.data(gsd0 + 81);
    const auto *gsd0_83 = buffer.data(gsd0 + 83);
    const auto *gsd0_84 = buffer.data(gsd0 + 84);
    const auto *gsd0_87 = buffer.data(gsd0 + 87);
    const auto *gsd0_89 = buffer.data(gsd0 + 89);

    const auto *gsp_0 = buffer.data(gsp + 0);
    const auto *gsp_1 = buffer.data(gsp + 1);
    const auto *gsp_2 = buffer.data(gsp + 2);
    const auto *gsp_4 = buffer.data(gsp + 4);
    const auto *gsp_8 = buffer.data(gsp + 8);
    const auto *gsp_9 = buffer.data(gsp + 9);
    const auto *gsp_10 = buffer.data(gsp + 10);
    const auto *gsp_11 = buffer.data(gsp + 11);
    const auto *gsp_13 = buffer.data(gsp + 13);
    const auto *gsp_14 = buffer.data(gsp + 14);
    const auto *gsp_15 = buffer.data(gsp + 15);
    const auto *gsp_16 = buffer.data(gsp + 16);
    const auto *gsp_17 = buffer.data(gsp + 17);
    const auto *gsp_18 = buffer.data(gsp + 18);
    const auto *gsp_19 = buffer.data(gsp + 19);
    const auto *gsp_22 = buffer.data(gsp + 22);
    const auto *gsp_23 = buffer.data(gsp + 23);
    const auto *gsp_25 = buffer.data(gsp + 25);
    const auto *gsp_26 = buffer.data(gsp + 26);
    const auto *gsp_27 = buffer.data(gsp + 27);
    const auto *gsp_29 = buffer.data(gsp + 29);
    const auto *gsp_30 = buffer.data(gsp + 30);
    const auto *gsp_31 = buffer.data(gsp + 31);
    const auto *gsp_32 = buffer.data(gsp + 32);
    const auto *gsp_34 = buffer.data(gsp + 34);
    const auto *gsp_35 = buffer.data(gsp + 35);
    const auto *gsp_36 = buffer.data(gsp + 36);
    const auto *gsp_37 = buffer.data(gsp + 37);
    const auto *gsp_38 = buffer.data(gsp + 38);
    const auto *gsp_40 = buffer.data(gsp + 40);
    const auto *gsp_41 = buffer.data(gsp + 41);
    const auto *gsp_42 = buffer.data(gsp + 42);
    const auto *gsp_43 = buffer.data(gsp + 43);
    const auto *gsp_44 = buffer.data(gsp + 44);

    const auto *gsd1_0 = buffer.data(gsd1 + 0);
    const auto *gsd1_3 = buffer.data(gsd1 + 3);
    const auto *gsd1_5 = buffer.data(gsd1 + 5);
    const auto *gsd1_9 = buffer.data(gsd1 + 9);
    const auto *gsd1_12 = buffer.data(gsd1 + 12);
    const auto *gsd1_17 = buffer.data(gsd1 + 17);
    const auto *gsd1_18 = buffer.data(gsd1 + 18);
    const auto *gsd1_21 = buffer.data(gsd1 + 21);
    const auto *gsd1_30 = buffer.data(gsd1 + 30);
    const auto *gsd1_35 = buffer.data(gsd1 + 35);
    const auto *gsd1_36 = buffer.data(gsd1 + 36);
    const auto *gsd1_54 = buffer.data(gsd1 + 54);
    const auto *gsd1_60 = buffer.data(gsd1 + 60);
    const auto *gsd1_63 = buffer.data(gsd1 + 63);
    const auto *gsd1_65 = buffer.data(gsd1 + 65);
    const auto *gsd1_69 = buffer.data(gsd1 + 69);
    const auto *gsd1_71 = buffer.data(gsd1 + 71);
    const auto *gsd1_72 = buffer.data(gsd1 + 72);
    const auto *gsd1_75 = buffer.data(gsd1 + 75);
    const auto *gsd1_77 = buffer.data(gsd1 + 77);
    const auto *gsd1_81 = buffer.data(gsd1 + 81);
    const auto *gsd1_83 = buffer.data(gsd1 + 83);
    const auto *gsd1_84 = buffer.data(gsd1 + 84);
    const auto *gsd1_87 = buffer.data(gsd1 + 87);
    const auto *gsd1_89 = buffer.data(gsd1 + 89);

    const auto *hss0_0 = buffer.data(hss0 + 0);
    const auto *hss0_1 = buffer.data(hss0 + 1);
    const auto *hss0_2 = buffer.data(hss0 + 2);
    const auto *hss0_3 = buffer.data(hss0 + 3);
    const auto *hss0_5 = buffer.data(hss0 + 5);
    const auto *hss0_6 = buffer.data(hss0 + 6);
    const auto *hss0_7 = buffer.data(hss0 + 7);
    const auto *hss0_8 = buffer.data(hss0 + 8);
    const auto *hss0_9 = buffer.data(hss0 + 9);
    const auto *hss0_15 = buffer.data(hss0 + 15);
    const auto *hss0_16 = buffer.data(hss0 + 16);
    const auto *hss0_17 = buffer.data(hss0 + 17);
    const auto *hss0_18 = buffer.data(hss0 + 18);
    const auto *hss0_20 = buffer.data(hss0 + 20);

    const auto *hss1_0 = buffer.data(hss1 + 0);
    const auto *hss1_1 = buffer.data(hss1 + 1);
    const auto *hss1_2 = buffer.data(hss1 + 2);
    const auto *hss1_3 = buffer.data(hss1 + 3);
    const auto *hss1_5 = buffer.data(hss1 + 5);
    const auto *hss1_6 = buffer.data(hss1 + 6);
    const auto *hss1_7 = buffer.data(hss1 + 7);
    const auto *hss1_8 = buffer.data(hss1 + 8);
    const auto *hss1_9 = buffer.data(hss1 + 9);
    const auto *hss1_15 = buffer.data(hss1 + 15);
    const auto *hss1_16 = buffer.data(hss1 + 16);
    const auto *hss1_17 = buffer.data(hss1 + 17);
    const auto *hss1_18 = buffer.data(hss1 + 18);
    const auto *hss1_20 = buffer.data(hss1 + 20);

    const auto *hsp_0 = buffer.data(hsp + 0);
    const auto *hsp_1 = buffer.data(hsp + 1);
    const auto *hsp_2 = buffer.data(hsp + 2);
    const auto *hsp_3 = buffer.data(hsp + 3);
    const auto *hsp_4 = buffer.data(hsp + 4);
    const auto *hsp_6 = buffer.data(hsp + 6);
    const auto *hsp_8 = buffer.data(hsp + 8);
    const auto *hsp_9 = buffer.data(hsp + 9);
    const auto *hsp_10 = buffer.data(hsp + 10);
    const auto *hsp_11 = buffer.data(hsp + 11);
    const auto *hsp_13 = buffer.data(hsp + 13);
    const auto *hsp_14 = buffer.data(hsp + 14);
    const auto *hsp_15 = buffer.data(hsp + 15);
    const auto *hsp_16 = buffer.data(hsp + 16);
    const auto *hsp_17 = buffer.data(hsp + 17);
    const auto *hsp_18 = buffer.data(hsp + 18);
    const auto *hsp_19 = buffer.data(hsp + 19);
    const auto *hsp_20 = buffer.data(hsp + 20);
    const auto *hsp_22 = buffer.data(hsp + 22);
    const auto *hsp_23 = buffer.data(hsp + 23);
    const auto *hsp_25 = buffer.data(hsp + 25);
    const auto *hsp_26 = buffer.data(hsp + 26);
    const auto *hsp_27 = buffer.data(hsp + 27);
    const auto *hsp_28 = buffer.data(hsp + 28);
    const auto *hsp_29 = buffer.data(hsp + 29);
    const auto *hsp_30 = buffer.data(hsp + 30);
    const auto *hsp_31 = buffer.data(hsp + 31);
    const auto *hsp_34 = buffer.data(hsp + 34);
    const auto *hsp_35 = buffer.data(hsp + 35);
    const auto *hsp_37 = buffer.data(hsp + 37);
    const auto *hsp_38 = buffer.data(hsp + 38);
    const auto *hsp_40 = buffer.data(hsp + 40);
    const auto *hsp_41 = buffer.data(hsp + 41);
    const auto *hsp_42 = buffer.data(hsp + 42);
    const auto *hsp_44 = buffer.data(hsp + 44);
    const auto *hsp_45 = buffer.data(hsp + 45);
    const auto *hsp_46 = buffer.data(hsp + 46);
    const auto *hsp_47 = buffer.data(hsp + 47);
    const auto *hsp_49 = buffer.data(hsp + 49);
    const auto *hsp_50 = buffer.data(hsp + 50);
    const auto *hsp_51 = buffer.data(hsp + 51);
    const auto *hsp_52 = buffer.data(hsp + 52);
    const auto *hsp_53 = buffer.data(hsp + 53);
    const auto *hsp_54 = buffer.data(hsp + 54);
    const auto *hsp_55 = buffer.data(hsp + 55);
    const auto *hsp_56 = buffer.data(hsp + 56);
    const auto *hsp_58 = buffer.data(hsp + 58);
    const auto *hsp_59 = buffer.data(hsp + 59);
    const auto *hsp_60 = buffer.data(hsp + 60);
    const auto *hsp_61 = buffer.data(hsp + 61);
    const auto *hsp_62 = buffer.data(hsp + 62);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, gsp_0, hss0_0, \
                         hss1_0, hsp_0, hsp_1, hsp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gsp_0[k]
                 + f_1 * hss0_0[k]
                 - f_2 * hss1_0[k]
                 + f_3 * pc_x[k] * hsp_0[k];

        t_1[k] = f_3 * pc_y[k] * hsp_0[k];

        t_2[k] = f_3 * pc_z[k] * hsp_0[k];

        t_3[k] = f_1 * hss0_0[k]
                 - f_2 * hss1_0[k]
                 + f_3 * pc_y[k] * hsp_1[k];

        t_4[k] = f_3 * pc_y[k] * hsp_2[k];

        t_5[k] = f_1 * hss0_0[k]
                 - f_2 * hss1_0[k]
                 + f_3 * pc_z[k] * hsp_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pa_y, pc_x, pc_y, pc_z, gsd0_0, gsp_1, gsp_4, \
                         gsd1_0, hss0_1, hss1_1, hsp_3, hsp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pa_y[k] * gsd0_0[k]
                 - f_4 * pc_y[k] * gsd1_0[k];

        t_7[k] = f_5 * gsp_4[k]
                 + f_3 * pc_x[k] * hsp_4[k];

        t_8[k] = f_3 * pc_z[k] * hsp_3[k];

        t_9[k] = f_6 * gsp_1[k]
                 + f_1 * hss0_1[k]
                 - f_2 * hss1_1[k]
                 + f_3 * pc_y[k] * hsp_4[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_y, pa_z, pc_y, pc_z, gsd0_0, gsd0_5, \
                         gsd1_0, gsd1_5, hsp_4, hsp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * pc_z[k] * hsp_4[k];

        t_11[k] = pa_y[k] * gsd0_5[k]
                  - f_4 * pc_y[k] * gsd1_5[k];

        t_12[k] = pa_z[k] * gsd0_0[k]
                  - f_4 * pc_z[k] * gsd1_0[k];

        t_13[k] = f_3 * pc_y[k] * hsp_6[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_z, pc_x, pc_y, pc_z, gsd0_3, gsp_2, gsp_8, \
                         gsd1_3, hss0_2, hss1_2, hsp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_5 * gsp_8[k]
                  + f_3 * pc_x[k] * hsp_8[k];

        t_15[k] = pa_z[k] * gsd0_3[k]
                  - f_4 * pc_z[k] * gsd1_3[k];

        t_16[k] = f_3 * pc_y[k] * hsp_8[k];

        t_17[k] = f_6 * gsp_2[k]
                  + f_1 * hss0_2[k]
                  - f_2 * hss1_2[k]
                  + f_3 * pc_z[k] * hsp_8[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, pc_x, pc_y, pc_z, gsp_4, gsp_9, gsp_10, \
                         hss0_3, hss1_3, hsp_9, hsp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_7 * gsp_9[k]
                  + f_1 * hss0_3[k]
                  - f_2 * hss1_3[k]
                  + f_3 * pc_x[k] * hsp_9[k];

        t_19[k] = f_7 * gsp_10[k]
                  + f_3 * pc_x[k] * hsp_10[k];

        t_20[k] = f_3 * pc_z[k] * hsp_9[k];

        t_21[k] = f_8 * gsp_4[k]
                  + f_1 * hss0_3[k]
                  - f_2 * hss1_3[k]
                  + f_3 * pc_y[k] * hsp_10[k];

        t_22[k] = f_3 * pc_z[k] * hsp_10[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_y, pc_x, pc_y, pc_z, gsd0_12, gsp_13, gsd1_12, \
                         hss0_3, hss1_3, hsp_11, hsp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_1 * hss0_3[k]
                  - f_2 * hss1_3[k]
                  + f_3 * pc_z[k] * hsp_11[k];

        t_24[k] = pa_y[k] * gsd0_12[k]
                  - f_4 * pc_y[k] * gsd1_12[k];

        t_25[k] = f_7 * gsp_13[k]
                  + f_3 * pc_x[k] * hsp_13[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_y, pa_z, pc_x, pc_y, pc_z, gsd0_9, \
                         gsd0_17, gsp_8, gsp_14, gsd1_9, gsd1_17, \
                         hsp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_7 * gsp_14[k]
                  + f_3 * pc_x[k] * hsp_14[k];

        t_27[k] = pa_z[k] * gsd0_9[k]
                  - f_4 * pc_z[k] * gsd1_9[k];

        t_28[k] = f_6 * gsp_8[k]
                  + f_3 * pc_y[k] * hsp_14[k];

        t_29[k] = pa_y[k] * gsd0_17[k]
                  - f_4 * pc_y[k] * gsd1_17[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pc_x, pc_y, gsp_15, gsp_17, hss0_5, \
                         hss1_5, hsp_15, hsp_16, hsp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_7 * gsp_15[k]
                  + f_1 * hss0_5[k]
                  - f_2 * hss1_5[k]
                  + f_3 * pc_x[k] * hsp_15[k];

        t_31[k] = f_3 * pc_y[k] * hsp_15[k];

        t_32[k] = f_7 * gsp_17[k]
                  + f_3 * pc_x[k] * hsp_17[k];

        t_33[k] = f_1 * hss0_5[k]
                  - f_2 * hss1_5[k]
                  + f_3 * pc_y[k] * hsp_16[k];

        t_34[k] = f_3 * pc_y[k] * hsp_17[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pc_x, pc_z, gsp_8, gsp_18, gsp_19, hss0_5, \
                         hss0_6, hss1_5, hss1_6, hsp_17, hsp_18, \
                         hsp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_8 * gsp_8[k]
                  + f_1 * hss0_5[k]
                  - f_2 * hss1_5[k]
                  + f_3 * pc_z[k] * hsp_17[k];

        t_36[k] = f_8 * gsp_18[k]
                  + f_1 * hss0_6[k]
                  - f_2 * hss1_6[k]
                  + f_3 * pc_x[k] * hsp_18[k];

        t_37[k] = f_8 * gsp_19[k]
                  + f_3 * pc_x[k] * hsp_19[k];

        t_38[k] = f_3 * pc_z[k] * hsp_18[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pa_z, pc_y, pc_z, gsd0_18, gsp_10, gsd1_18, \
                         hss0_6, hss1_6, hsp_19, hsp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_7 * gsp_10[k]
                  + f_1 * hss0_6[k]
                  - f_2 * hss1_6[k]
                  + f_3 * pc_y[k] * hsp_19[k];

        t_40[k] = f_3 * pc_z[k] * hsp_19[k];

        t_41[k] = f_1 * hss0_6[k]
                  - f_2 * hss1_6[k]
                  + f_3 * pc_z[k] * hsp_20[k];

        t_42[k] = pa_z[k] * gsd0_18[k]
                  - f_4 * pc_z[k] * gsd1_18[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, pa_z, pc_x, pc_y, pc_z, gsd0_21, gsp_14, \
                         gsp_22, gsp_23, gsd1_21, hsp_22, hsp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_8 * gsp_22[k]
                  + f_3 * pc_x[k] * hsp_22[k];

        t_44[k] = f_8 * gsp_23[k]
                  + f_3 * pc_x[k] * hsp_23[k];

        t_45[k] = pa_z[k] * gsd0_21[k]
                  - f_4 * pc_z[k] * gsd1_21[k];

        t_46[k] = f_8 * gsp_14[k]
                  + f_3 * pc_y[k] * hsp_23[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, pa_y, pc_x, pc_y, pc_z, gsd0_30, gsp_11, gsp_25, \
                         gsd1_30, hss0_7, hss1_7, hsp_23, hsp_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_6 * gsp_11[k]
                  + f_1 * hss0_7[k]
                  - f_2 * hss1_7[k]
                  + f_3 * pc_z[k] * hsp_23[k];

        t_48[k] = pa_y[k] * gsd0_30[k]
                  - f_4 * pc_y[k] * gsd1_30[k];

        t_49[k] = f_8 * gsp_25[k]
                  + f_3 * pc_x[k] * hsp_25[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_y, pc_x, pc_y, gsd0_35, gsp_16, gsp_17, \
                         gsp_26, gsd1_35, hss0_8, hss1_8, hsp_25, \
                         hsp_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_8 * gsp_26[k]
                  + f_3 * pc_x[k] * hsp_26[k];

        t_51[k] = f_6 * gsp_16[k]
                  + f_1 * hss0_8[k]
                  - f_2 * hss1_8[k]
                  + f_3 * pc_y[k] * hsp_25[k];

        t_52[k] = f_6 * gsp_17[k]
                  + f_3 * pc_y[k] * hsp_26[k];

        t_53[k] = pa_y[k] * gsd0_35[k]
                  - f_4 * pc_y[k] * gsd1_35[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, pc_x, pc_y, gsp_27, gsp_29, hss0_9, \
                         hss1_9, hsp_27, hsp_28, hsp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_8 * gsp_27[k]
                  + f_1 * hss0_9[k]
                  - f_2 * hss1_9[k]
                  + f_3 * pc_x[k] * hsp_27[k];

        t_55[k] = f_3 * pc_y[k] * hsp_27[k];

        t_56[k] = f_8 * gsp_29[k]
                  + f_3 * pc_x[k] * hsp_29[k];

        t_57[k] = f_1 * hss0_9[k]
                  - f_2 * hss1_9[k]
                  + f_3 * pc_y[k] * hsp_28[k];

        t_58[k] = f_3 * pc_y[k] * hsp_29[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, pa_x, pc_x, pc_z, gsd0_60, gsp_17, gsp_30, gsp_31, \
                         gsd1_60, hss0_9, hss1_9, hsp_29, hsp_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_7 * gsp_17[k]
                  + f_1 * hss0_9[k]
                  - f_2 * hss1_9[k]
                  + f_3 * pc_z[k] * hsp_29[k];

        t_60[k] = pa_x[k] * gsd0_60[k]
                  + f_8 * gsp_30[k]
                  - f_4 * pc_x[k] * gsd1_60[k];

        t_61[k] = f_6 * gsp_31[k]
                  + f_3 * pc_x[k] * hsp_31[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pa_x, pc_x, pc_z, gsd0_63, gsd0_65, gsd1_63, \
                         gsd1_65, hsp_30, hsp_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_3 * pc_z[k] * hsp_30[k];

        t_63[k] = pa_x[k] * gsd0_63[k]
                  - f_4 * pc_x[k] * gsd1_63[k];

        t_64[k] = f_3 * pc_z[k] * hsp_31[k];

        t_65[k] = pa_x[k] * gsd0_65[k]
                  - f_4 * pc_x[k] * gsd1_65[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pa_x, pa_z, pc_x, pc_z, gsd0_36, gsd0_69, \
                         gsp_34, gsp_35, gsd1_36, gsd1_69, hsp_34, \
                         hsp_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pa_z[k] * gsd0_36[k]
                  - f_4 * pc_z[k] * gsd1_36[k];

        t_67[k] = f_6 * gsp_34[k]
                  + f_3 * pc_x[k] * hsp_34[k];

        t_68[k] = f_6 * gsp_35[k]
                  + f_3 * pc_x[k] * hsp_35[k];

        t_69[k] = pa_x[k] * gsd0_69[k]
                  - f_4 * pc_x[k] * gsd1_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pa_x, pc_x, pc_y, gsd0_71, gsd0_72, gsp_23, \
                         gsp_36, gsp_37, gsd1_71, gsd1_72, hsp_35, \
                         hsp_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_7 * gsp_23[k]
                  + f_3 * pc_y[k] * hsp_35[k];

        t_71[k] = pa_x[k] * gsd0_71[k]
                  - f_4 * pc_x[k] * gsd1_71[k];

        t_72[k] = pa_x[k] * gsd0_72[k]
                  + f_8 * gsp_36[k]
                  - f_4 * pc_x[k] * gsd1_72[k];

        t_73[k] = f_6 * gsp_37[k]
                  + f_3 * pc_x[k] * hsp_37[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pa_x, pc_x, pc_y, gsd0_75, gsd0_77, gsp_26, \
                         gsp_38, gsd1_75, gsd1_77, hsp_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_6 * gsp_38[k]
                  + f_3 * pc_x[k] * hsp_38[k];

        t_75[k] = pa_x[k] * gsd0_75[k]
                  - f_4 * pc_x[k] * gsd1_75[k];

        t_76[k] = f_8 * gsp_26[k]
                  + f_3 * pc_y[k] * hsp_38[k];

        t_77[k] = pa_x[k] * gsd0_77[k]
                  - f_4 * pc_x[k] * gsd1_77[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pa_x, pa_y, pc_x, pc_y, gsd0_54, gsd0_81, \
                         gsp_40, gsp_41, gsd1_54, gsd1_81, hsp_40, \
                         hsp_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = pa_y[k] * gsd0_54[k]
                  - f_4 * pc_y[k] * gsd1_54[k];

        t_79[k] = f_6 * gsp_40[k]
                  + f_3 * pc_x[k] * hsp_40[k];

        t_80[k] = f_6 * gsp_41[k]
                  + f_3 * pc_x[k] * hsp_41[k];

        t_81[k] = pa_x[k] * gsd0_81[k]
                  - f_4 * pc_x[k] * gsd1_81[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, pa_x, pc_x, pc_y, gsd0_83, gsd0_84, gsp_29, \
                         gsp_42, gsd1_83, gsd1_84, hsp_41, hsp_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_6 * gsp_29[k]
                  + f_3 * pc_y[k] * hsp_41[k];

        t_83[k] = pa_x[k] * gsd0_83[k]
                  - f_4 * pc_x[k] * gsd1_83[k];

        t_84[k] = pa_x[k] * gsd0_84[k]
                  + f_8 * gsp_42[k]
                  - f_4 * pc_x[k] * gsd1_84[k];

        t_85[k] = f_3 * pc_y[k] * hsp_42[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pa_x, pc_x, pc_y, gsd0_87, gsd0_89, gsp_44, \
                         gsd1_87, gsd1_89, hsp_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_6 * gsp_44[k]
                  + f_3 * pc_x[k] * hsp_44[k];

        t_87[k] = pa_x[k] * gsd0_87[k]
                  - f_4 * pc_x[k] * gsd1_87[k];

        t_88[k] = f_3 * pc_y[k] * hsp_44[k];

        t_89[k] = pa_x[k] * gsd0_89[k]
                  - f_4 * pc_x[k] * gsd1_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, t_95, pc_x, pc_y, pc_z, gsp_31, \
                         hss0_15, hss1_15, hsp_45, hsp_46, hsp_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_1 * hss0_15[k]
                  - f_2 * hss1_15[k]
                  + f_3 * pc_x[k] * hsp_45[k];

        t_91[k] = f_3 * pc_x[k] * hsp_46[k];

        t_92[k] = f_3 * pc_x[k] * hsp_47[k];

        t_93[k] = f_0 * gsp_31[k]
                  + f_1 * hss0_15[k]
                  - f_2 * hss1_15[k]
                  + f_3 * pc_y[k] * hsp_46[k];

        t_94[k] = f_3 * pc_z[k] * hsp_46[k];

        t_95[k] = f_1 * hss0_15[k]
                  - f_2 * hss1_15[k]
                  + f_3 * pc_z[k] * hsp_47[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, t_100, pa_z, pc_x, pc_y, pc_z, gsd0_60, \
                         gsd0_63, gsp_35, gsd1_60, gsd1_63, hsp_49, \
                         hsp_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = pa_z[k] * gsd0_60[k]
                  - f_4 * pc_z[k] * gsd1_60[k];

        t_97[k] = f_3 * pc_x[k] * hsp_49[k];

        t_98[k] = f_3 * pc_x[k] * hsp_50[k];

        t_99[k] = pa_z[k] * gsd0_63[k]
                  - f_4 * pc_z[k] * gsd1_63[k];

        t_100[k] = f_5 * gsp_35[k]
                   + f_3 * pc_y[k] * hsp_50[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, pc_x, pc_z, gsp_32, hss0_16, hss0_17, \
                         hss1_16, hss1_17, hsp_50, hsp_51, hsp_52, \
                         hsp_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_6 * gsp_32[k]
                   + f_1 * hss0_16[k]
                   - f_2 * hss1_16[k]
                   + f_3 * pc_z[k] * hsp_50[k];

        t_102[k] = f_1 * hss0_17[k]
                   - f_2 * hss1_17[k]
                   + f_3 * pc_x[k] * hsp_51[k];

        t_103[k] = f_3 * pc_x[k] * hsp_52[k];

        t_104[k] = f_3 * pc_x[k] * hsp_53[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, pc_y, pc_z, gsp_35, gsp_37, gsp_38, hss0_17, \
                         hss1_17, hsp_52, hsp_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_7 * gsp_37[k]
                   + f_1 * hss0_17[k]
                   - f_2 * hss1_17[k]
                   + f_3 * pc_y[k] * hsp_52[k];

        t_106[k] = f_7 * gsp_38[k]
                   + f_3 * pc_y[k] * hsp_53[k];

        t_107[k] = f_8 * gsp_35[k]
                   + f_1 * hss0_17[k]
                   - f_2 * hss1_17[k]
                   + f_3 * pc_z[k] * hsp_53[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, t_112, pc_x, pc_y, gsp_40, gsp_41, \
                         hss0_18, hss1_18, hsp_54, hsp_55, hsp_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_1 * hss0_18[k]
                   - f_2 * hss1_18[k]
                   + f_3 * pc_x[k] * hsp_54[k];

        t_109[k] = f_3 * pc_x[k] * hsp_55[k];

        t_110[k] = f_3 * pc_x[k] * hsp_56[k];

        t_111[k] = f_8 * gsp_40[k]
                   + f_1 * hss0_18[k]
                   - f_2 * hss1_18[k]
                   + f_3 * pc_y[k] * hsp_55[k];

        t_112[k] = f_8 * gsp_41[k]
                   + f_3 * pc_y[k] * hsp_56[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pa_y, pc_x, pc_y, pc_z, gsd0_84, gsp_38, \
                         gsd1_84, hss0_18, hss1_18, hsp_56, hsp_58, \
                         hsp_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_7 * gsp_38[k]
                   + f_1 * hss0_18[k]
                   - f_2 * hss1_18[k]
                   + f_3 * pc_z[k] * hsp_56[k];

        t_114[k] = pa_y[k] * gsd0_84[k]
                   - f_4 * pc_y[k] * gsd1_84[k];

        t_115[k] = f_3 * pc_x[k] * hsp_58[k];

        t_116[k] = f_3 * pc_x[k] * hsp_59[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, pa_y, pc_y, gsd0_87, gsd0_89, gsp_43, gsp_44, \
                         gsd1_87, gsd1_89, hsp_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = pa_y[k] * gsd0_87[k]
                   + f_8 * gsp_43[k]
                   - f_4 * pc_y[k] * gsd1_87[k];

        t_118[k] = f_6 * gsp_44[k]
                   + f_3 * pc_y[k] * hsp_59[k];

        t_119[k] = pa_y[k] * gsd0_89[k]
                   - f_4 * pc_y[k] * gsd1_89[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, t_125, pc_x, pc_y, pc_z, gsp_44, \
                         hss0_20, hss1_20, hsp_60, hsp_61, hsp_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_1 * hss0_20[k]
                   - f_2 * hss1_20[k]
                   + f_3 * pc_x[k] * hsp_60[k];

        t_121[k] = f_3 * pc_x[k] * hsp_61[k];

        t_122[k] = f_3 * pc_x[k] * hsp_62[k];

        t_123[k] = f_1 * hss0_20[k]
                   - f_2 * hss1_20[k]
                   + f_3 * pc_y[k] * hsp_61[k];

        t_124[k] = f_3 * pc_y[k] * hsp_62[k];

        t_125[k] = f_0 * gsp_44[k]
                   + f_1 * hss0_20[k]
                   - f_2 * hss1_20[k]
                   + f_3 * pc_z[k] * hsp_62[k];
    }
}

}  // namespace simdt3ceri
