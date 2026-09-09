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


#include "SimdThreeCenterElectronRepulsionVrrRecHSG.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_hsg_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t gsg0,
                                                          const size_t gsf, const size_t gsg1,
                                                          const size_t hsd0, const size_t hsd1,
                                                          const size_t hsf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 2.0 / q;
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
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

    const auto *gsg0_0 = buffer.data(gsg0 + 0);
    const auto *gsg0_3 = buffer.data(gsg0 + 3);
    const auto *gsg0_5 = buffer.data(gsg0 + 5);
    const auto *gsg0_10 = buffer.data(gsg0 + 10);
    const auto *gsg0_14 = buffer.data(gsg0 + 14);
    const auto *gsg0_18 = buffer.data(gsg0 + 18);
    const auto *gsg0_25 = buffer.data(gsg0 + 25);
    const auto *gsg0_30 = buffer.data(gsg0 + 30);
    const auto *gsg0_35 = buffer.data(gsg0 + 35);
    const auto *gsg0_44 = buffer.data(gsg0 + 44);
    const auto *gsg0_45 = buffer.data(gsg0 + 45);
    const auto *gsg0_48 = buffer.data(gsg0 + 48);
    const auto *gsg0_55 = buffer.data(gsg0 + 55);
    const auto *gsg0_75 = buffer.data(gsg0 + 75);
    const auto *gsg0_78 = buffer.data(gsg0 + 78);
    const auto *gsg0_80 = buffer.data(gsg0 + 80);

    const auto *gsf_0 = buffer.data(gsf + 0);
    const auto *gsf_1 = buffer.data(gsf + 1);
    const auto *gsf_2 = buffer.data(gsf + 2);
    const auto *gsf_6 = buffer.data(gsf + 6);
    const auto *gsf_9 = buffer.data(gsf + 9);
    const auto *gsf_10 = buffer.data(gsf + 10);
    const auto *gsf_16 = buffer.data(gsf + 16);
    const auto *gsf_18 = buffer.data(gsf + 18);
    const auto *gsf_19 = buffer.data(gsf + 19);
    const auto *gsf_20 = buffer.data(gsf + 20);
    const auto *gsf_22 = buffer.data(gsf + 22);
    const auto *gsf_26 = buffer.data(gsf + 26);
    const auto *gsf_27 = buffer.data(gsf + 27);
    const auto *gsf_28 = buffer.data(gsf + 28);
    const auto *gsf_29 = buffer.data(gsf + 29);
    const auto *gsf_30 = buffer.data(gsf + 30);
    const auto *gsf_33 = buffer.data(gsf + 33);
    const auto *gsf_36 = buffer.data(gsf + 36);
    const auto *gsf_38 = buffer.data(gsf + 38);
    const auto *gsf_39 = buffer.data(gsf + 39);
    const auto *gsf_40 = buffer.data(gsf + 40);
    const auto *gsf_42 = buffer.data(gsf + 42);
    const auto *gsf_46 = buffer.data(gsf + 46);
    const auto *gsf_47 = buffer.data(gsf + 47);
    const auto *gsf_48 = buffer.data(gsf + 48);
    const auto *gsf_49 = buffer.data(gsf + 49);
    const auto *gsf_50 = buffer.data(gsf + 50);
    const auto *gsf_51 = buffer.data(gsf + 51);
    const auto *gsf_52 = buffer.data(gsf + 52);
    const auto *gsf_55 = buffer.data(gsf + 55);
    const auto *gsf_56 = buffer.data(gsf + 56);
    const auto *gsf_57 = buffer.data(gsf + 57);
    const auto *gsf_59 = buffer.data(gsf + 59);
    const auto *gsf_60 = buffer.data(gsf + 60);
    const auto *gsf_63 = buffer.data(gsf + 63);
    const auto *gsf_66 = buffer.data(gsf + 66);
    const auto *gsf_68 = buffer.data(gsf + 68);
    const auto *gsf_69 = buffer.data(gsf + 69);
    const auto *gsf_75 = buffer.data(gsf + 75);
    const auto *gsf_76 = buffer.data(gsf + 76);
    const auto *gsf_77 = buffer.data(gsf + 77);
    const auto *gsf_78 = buffer.data(gsf + 78);
    const auto *gsf_79 = buffer.data(gsf + 79);
    const auto *gsf_86 = buffer.data(gsf + 86);
    const auto *gsf_87 = buffer.data(gsf + 87);
    const auto *gsf_88 = buffer.data(gsf + 88);
    const auto *gsf_89 = buffer.data(gsf + 89);

    const auto *gsg1_0 = buffer.data(gsg1 + 0);
    const auto *gsg1_3 = buffer.data(gsg1 + 3);
    const auto *gsg1_5 = buffer.data(gsg1 + 5);
    const auto *gsg1_10 = buffer.data(gsg1 + 10);
    const auto *gsg1_14 = buffer.data(gsg1 + 14);
    const auto *gsg1_18 = buffer.data(gsg1 + 18);
    const auto *gsg1_25 = buffer.data(gsg1 + 25);
    const auto *gsg1_30 = buffer.data(gsg1 + 30);
    const auto *gsg1_35 = buffer.data(gsg1 + 35);
    const auto *gsg1_44 = buffer.data(gsg1 + 44);
    const auto *gsg1_45 = buffer.data(gsg1 + 45);
    const auto *gsg1_48 = buffer.data(gsg1 + 48);
    const auto *gsg1_55 = buffer.data(gsg1 + 55);
    const auto *gsg1_75 = buffer.data(gsg1 + 75);
    const auto *gsg1_78 = buffer.data(gsg1 + 78);
    const auto *gsg1_80 = buffer.data(gsg1 + 80);

    const auto *hsd0_0 = buffer.data(hsd0 + 0);
    const auto *hsd0_3 = buffer.data(hsd0 + 3);
    const auto *hsd0_5 = buffer.data(hsd0 + 5);
    const auto *hsd0_9 = buffer.data(hsd0 + 9);
    const auto *hsd0_16 = buffer.data(hsd0 + 16);
    const auto *hsd0_17 = buffer.data(hsd0 + 17);
    const auto *hsd0_18 = buffer.data(hsd0 + 18);
    const auto *hsd0_21 = buffer.data(hsd0 + 21);
    const auto *hsd0_23 = buffer.data(hsd0 + 23);
    const auto *hsd0_29 = buffer.data(hsd0 + 29);
    const auto *hsd0_30 = buffer.data(hsd0 + 30);
    const auto *hsd0_33 = buffer.data(hsd0 + 33);
    const auto *hsd0_34 = buffer.data(hsd0 + 34);
    const auto *hsd0_35 = buffer.data(hsd0 + 35);
    const auto *hsd0_36 = buffer.data(hsd0 + 36);
    const auto *hsd0_39 = buffer.data(hsd0 + 39);
    const auto *hsd0_41 = buffer.data(hsd0 + 41);
    const auto *hsd0_47 = buffer.data(hsd0 + 47);

    const auto *hsd1_0 = buffer.data(hsd1 + 0);
    const auto *hsd1_3 = buffer.data(hsd1 + 3);
    const auto *hsd1_5 = buffer.data(hsd1 + 5);
    const auto *hsd1_9 = buffer.data(hsd1 + 9);
    const auto *hsd1_16 = buffer.data(hsd1 + 16);
    const auto *hsd1_17 = buffer.data(hsd1 + 17);
    const auto *hsd1_18 = buffer.data(hsd1 + 18);
    const auto *hsd1_21 = buffer.data(hsd1 + 21);
    const auto *hsd1_23 = buffer.data(hsd1 + 23);
    const auto *hsd1_29 = buffer.data(hsd1 + 29);
    const auto *hsd1_30 = buffer.data(hsd1 + 30);
    const auto *hsd1_33 = buffer.data(hsd1 + 33);
    const auto *hsd1_34 = buffer.data(hsd1 + 34);
    const auto *hsd1_35 = buffer.data(hsd1 + 35);
    const auto *hsd1_36 = buffer.data(hsd1 + 36);
    const auto *hsd1_39 = buffer.data(hsd1 + 39);
    const auto *hsd1_41 = buffer.data(hsd1 + 41);
    const auto *hsd1_47 = buffer.data(hsd1 + 47);

    const auto *hsf_0 = buffer.data(hsf + 0);
    const auto *hsf_1 = buffer.data(hsf + 1);
    const auto *hsf_2 = buffer.data(hsf + 2);
    const auto *hsf_3 = buffer.data(hsf + 3);
    const auto *hsf_5 = buffer.data(hsf + 5);
    const auto *hsf_6 = buffer.data(hsf + 6);
    const auto *hsf_8 = buffer.data(hsf + 8);
    const auto *hsf_9 = buffer.data(hsf + 9);
    const auto *hsf_10 = buffer.data(hsf + 10);
    const auto *hsf_11 = buffer.data(hsf + 11);
    const auto *hsf_13 = buffer.data(hsf + 13);
    const auto *hsf_16 = buffer.data(hsf + 16);
    const auto *hsf_17 = buffer.data(hsf + 17);
    const auto *hsf_18 = buffer.data(hsf + 18);
    const auto *hsf_19 = buffer.data(hsf + 19);
    const auto *hsf_20 = buffer.data(hsf + 20);
    const auto *hsf_22 = buffer.data(hsf + 22);
    const auto *hsf_25 = buffer.data(hsf + 25);
    const auto *hsf_26 = buffer.data(hsf + 26);
    const auto *hsf_27 = buffer.data(hsf + 27);
    const auto *hsf_28 = buffer.data(hsf + 28);
    const auto *hsf_29 = buffer.data(hsf + 29);
    const auto *hsf_30 = buffer.data(hsf + 30);
    const auto *hsf_31 = buffer.data(hsf + 31);
    const auto *hsf_32 = buffer.data(hsf + 32);
    const auto *hsf_33 = buffer.data(hsf + 33);
    const auto *hsf_36 = buffer.data(hsf + 36);
    const auto *hsf_37 = buffer.data(hsf + 37);
    const auto *hsf_38 = buffer.data(hsf + 38);
    const auto *hsf_39 = buffer.data(hsf + 39);
    const auto *hsf_40 = buffer.data(hsf + 40);
    const auto *hsf_42 = buffer.data(hsf + 42);
    const auto *hsf_46 = buffer.data(hsf + 46);
    const auto *hsf_47 = buffer.data(hsf + 47);
    const auto *hsf_48 = buffer.data(hsf + 48);
    const auto *hsf_49 = buffer.data(hsf + 49);
    const auto *hsf_50 = buffer.data(hsf + 50);
    const auto *hsf_51 = buffer.data(hsf + 51);
    const auto *hsf_52 = buffer.data(hsf + 52);
    const auto *hsf_55 = buffer.data(hsf + 55);
    const auto *hsf_56 = buffer.data(hsf + 56);
    const auto *hsf_57 = buffer.data(hsf + 57);
    const auto *hsf_58 = buffer.data(hsf + 58);
    const auto *hsf_59 = buffer.data(hsf + 59);
    const auto *hsf_60 = buffer.data(hsf + 60);
    const auto *hsf_61 = buffer.data(hsf + 61);
    const auto *hsf_62 = buffer.data(hsf + 62);
    const auto *hsf_63 = buffer.data(hsf + 63);
    const auto *hsf_66 = buffer.data(hsf + 66);
    const auto *hsf_67 = buffer.data(hsf + 67);
    const auto *hsf_68 = buffer.data(hsf + 68);
    const auto *hsf_69 = buffer.data(hsf + 69);
    const auto *hsf_70 = buffer.data(hsf + 70);
    const auto *hsf_72 = buffer.data(hsf + 72);
    const auto *hsf_75 = buffer.data(hsf + 75);
    const auto *hsf_76 = buffer.data(hsf + 76);
    const auto *hsf_77 = buffer.data(hsf + 77);
    const auto *hsf_78 = buffer.data(hsf + 78);
    const auto *hsf_79 = buffer.data(hsf + 79);
    const auto *hsf_80 = buffer.data(hsf + 80);
    const auto *hsf_82 = buffer.data(hsf + 82);
    const auto *hsf_86 = buffer.data(hsf + 86);
    const auto *hsf_87 = buffer.data(hsf + 87);
    const auto *hsf_88 = buffer.data(hsf + 88);
    const auto *hsf_89 = buffer.data(hsf + 89);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, gsf_0, hsd0_0, \
                         hsd1_0, hsf_0, hsf_1, hsf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gsf_0[k]
                 + f_1 * hsd0_0[k]
                 - f_2 * hsd1_0[k]
                 + f_3 * pc_x[k] * hsf_0[k];

        t_1[k] = f_3 * pc_y[k] * hsf_0[k];

        t_2[k] = f_3 * pc_z[k] * hsf_0[k];

        t_3[k] = f_4 * hsd0_0[k]
                 - f_5 * hsd1_0[k]
                 + f_3 * pc_y[k] * hsf_1[k];

        t_4[k] = f_3 * pc_y[k] * hsf_2[k];

        t_5[k] = f_4 * hsd0_0[k]
                 - f_5 * hsd1_0[k]
                 + f_3 * pc_z[k] * hsf_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_x, pc_y, pc_z, gsf_6, gsf_9, hsd0_3, \
                         hsd1_3, hsf_3, hsf_5, hsf_6, hsf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * gsf_6[k]
                 + f_3 * pc_x[k] * hsf_6[k];

        t_7[k] = f_3 * pc_z[k] * hsf_3[k];

        t_8[k] = f_3 * pc_y[k] * hsf_5[k];

        t_9[k] = f_0 * gsf_9[k]
                 + f_3 * pc_x[k] * hsf_9[k];

        t_10[k] = f_1 * hsd0_3[k]
                  - f_2 * hsd1_3[k]
                  + f_3 * pc_y[k] * hsf_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pa_y, pc_y, pc_z, gsg0_0, gsg1_0, \
                         hsd0_5, hsd1_5, hsf_6, hsf_8, hsf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * pc_z[k] * hsf_6[k];

        t_12[k] = f_4 * hsd0_5[k]
                  - f_5 * hsd1_5[k]
                  + f_3 * pc_y[k] * hsf_8[k];

        t_13[k] = f_3 * pc_y[k] * hsf_9[k];

        t_14[k] = f_1 * hsd0_5[k]
                  - f_2 * hsd1_5[k]
                  + f_3 * pc_z[k] * hsf_9[k];

        t_15[k] = pa_y[k] * gsg0_0[k]
                  - f_6 * pc_y[k] * gsg1_0[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pa_y, pc_y, pc_z, gsg0_3, gsg0_5, \
                         gsf_0, gsf_1, gsg1_3, gsg1_5, hsf_10, hsf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_7 * gsf_0[k]
                  + f_3 * pc_y[k] * hsf_10[k];

        t_17[k] = f_3 * pc_z[k] * hsf_10[k];

        t_18[k] = pa_y[k] * gsg0_3[k]
                  + f_8 * gsf_1[k]
                  - f_6 * pc_y[k] * gsg1_3[k];

        t_19[k] = f_3 * pc_z[k] * hsf_11[k];

        t_20[k] = pa_y[k] * gsg0_5[k]
                  - f_6 * pc_y[k] * gsg1_5[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pc_x, pc_z, gsf_16, gsf_18, gsf_19, hsf_13, \
                         hsf_16, hsf_18, hsf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_9 * gsf_16[k]
                  + f_3 * pc_x[k] * hsf_16[k];

        t_22[k] = f_3 * pc_z[k] * hsf_13[k];

        t_23[k] = f_9 * gsf_18[k]
                  + f_3 * pc_x[k] * hsf_18[k];

        t_24[k] = f_9 * gsf_19[k]
                  + f_3 * pc_x[k] * hsf_19[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pc_y, pc_z, gsf_6, gsf_9, hsd0_9, hsd1_9, \
                         hsf_16, hsf_17, hsf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_7 * gsf_6[k]
                  + f_1 * hsd0_9[k]
                  - f_2 * hsd1_9[k]
                  + f_3 * pc_y[k] * hsf_16[k];

        t_26[k] = f_3 * pc_z[k] * hsf_16[k];

        t_27[k] = f_4 * hsd0_9[k]
                  - f_5 * hsd1_9[k]
                  + f_3 * pc_z[k] * hsf_17[k];

        t_28[k] = f_7 * gsf_9[k]
                  + f_3 * pc_y[k] * hsf_19[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pa_y, pa_z, pc_y, pc_z, gsg0_0, gsg0_14, \
                         gsf_0, gsg1_0, gsg1_14, hsf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pa_y[k] * gsg0_14[k]
                  - f_6 * pc_y[k] * gsg1_14[k];

        t_30[k] = pa_z[k] * gsg0_0[k]
                  - f_6 * pc_z[k] * gsg1_0[k];

        t_31[k] = f_3 * pc_y[k] * hsf_20[k];

        t_32[k] = f_7 * gsf_0[k]
                  + f_3 * pc_z[k] * hsf_20[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pa_z, pc_x, pc_y, pc_z, gsg0_3, gsg0_5, \
                         gsf_2, gsf_26, gsg1_3, gsg1_5, hsf_22, \
                         hsf_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pa_z[k] * gsg0_3[k]
                  - f_6 * pc_z[k] * gsg1_3[k];

        t_34[k] = f_3 * pc_y[k] * hsf_22[k];

        t_35[k] = pa_z[k] * gsg0_5[k]
                  + f_8 * gsf_2[k]
                  - f_6 * pc_z[k] * gsg1_5[k];

        t_36[k] = f_9 * gsf_26[k]
                  + f_3 * pc_x[k] * hsf_26[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pa_z, pc_x, pc_y, pc_z, gsg0_10, gsf_27, \
                         gsf_29, gsg1_10, hsf_25, hsf_27, hsf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_9 * gsf_27[k]
                  + f_3 * pc_x[k] * hsf_27[k];

        t_38[k] = f_3 * pc_y[k] * hsf_25[k];

        t_39[k] = f_9 * gsf_29[k]
                  + f_3 * pc_x[k] * hsf_29[k];

        t_40[k] = pa_z[k] * gsg0_10[k]
                  - f_6 * pc_z[k] * gsg1_10[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pc_y, pc_z, gsf_9, hsd0_16, hsd0_17, hsd1_16, \
                         hsd1_17, hsf_27, hsf_28, hsf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_10 * hsd0_16[k]
                  - f_11 * hsd1_16[k]
                  + f_3 * pc_y[k] * hsf_27[k];

        t_42[k] = f_4 * hsd0_17[k]
                  - f_5 * hsd1_17[k]
                  + f_3 * pc_y[k] * hsf_28[k];

        t_43[k] = f_3 * pc_y[k] * hsf_29[k];

        t_44[k] = f_7 * gsf_9[k]
                  + f_1 * hsd0_17[k]
                  - f_2 * hsd1_17[k]
                  + f_3 * pc_z[k] * hsf_29[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pc_x, pc_y, pc_z, gsf_10, gsf_30, gsf_33, \
                         hsd0_18, hsd0_21, hsd1_18, hsd1_21, hsf_30, \
                         hsf_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_12 * gsf_30[k]
                  + f_1 * hsd0_18[k]
                  - f_2 * hsd1_18[k]
                  + f_3 * pc_x[k] * hsf_30[k];

        t_46[k] = f_8 * gsf_10[k]
                  + f_3 * pc_y[k] * hsf_30[k];

        t_47[k] = f_3 * pc_z[k] * hsf_30[k];

        t_48[k] = f_12 * gsf_33[k]
                  + f_4 * hsd0_21[k]
                  - f_5 * hsd1_21[k]
                  + f_3 * pc_x[k] * hsf_33[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, pc_x, pc_z, gsf_36, gsf_38, hsd0_18, \
                         hsd1_18, hsf_31, hsf_32, hsf_33, hsf_36, \
                         hsf_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_3 * pc_z[k] * hsf_31[k];

        t_50[k] = f_4 * hsd0_18[k]
                  - f_5 * hsd1_18[k]
                  + f_3 * pc_z[k] * hsf_32[k];

        t_51[k] = f_12 * gsf_36[k]
                  + f_3 * pc_x[k] * hsf_36[k];

        t_52[k] = f_3 * pc_z[k] * hsf_33[k];

        t_53[k] = f_12 * gsf_38[k]
                  + f_3 * pc_x[k] * hsf_38[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, pc_x, pc_y, pc_z, gsf_16, gsf_19, \
                         gsf_39, hsd0_21, hsd1_21, hsf_36, hsf_37, \
                         hsf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_12 * gsf_39[k]
                  + f_3 * pc_x[k] * hsf_39[k];

        t_55[k] = f_8 * gsf_16[k]
                  + f_1 * hsd0_21[k]
                  - f_2 * hsd1_21[k]
                  + f_3 * pc_y[k] * hsf_36[k];

        t_56[k] = f_3 * pc_z[k] * hsf_36[k];

        t_57[k] = f_4 * hsd0_21[k]
                  - f_5 * hsd1_21[k]
                  + f_3 * pc_z[k] * hsf_37[k];

        t_58[k] = f_8 * gsf_19[k]
                  + f_3 * pc_y[k] * hsf_39[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pa_y, pc_y, pc_z, gsg0_30, gsf_10, gsf_20, \
                         gsg1_30, hsd0_23, hsd1_23, hsf_39, hsf_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_1 * hsd0_23[k]
                  - f_2 * hsd1_23[k]
                  + f_3 * pc_z[k] * hsf_39[k];

        t_60[k] = pa_y[k] * gsg0_30[k]
                  - f_6 * pc_y[k] * gsg1_30[k];

        t_61[k] = f_7 * gsf_20[k]
                  + f_3 * pc_y[k] * hsf_40[k];

        t_62[k] = f_7 * gsf_10[k]
                  + f_3 * pc_z[k] * hsf_40[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, pa_y, pa_z, pc_y, pc_z, gsg0_18, gsg0_35, gsf_22, \
                         gsg1_18, gsg1_35, hsf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = pa_z[k] * gsg0_18[k]
                  - f_6 * pc_z[k] * gsg1_18[k];

        t_64[k] = f_7 * gsf_22[k]
                  + f_3 * pc_y[k] * hsf_42[k];

        t_65[k] = pa_y[k] * gsg0_35[k]
                  - f_6 * pc_y[k] * gsg1_35[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pc_x, gsf_46, gsf_47, gsf_48, gsf_49, hsf_46, \
                         hsf_47, hsf_48, hsf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_12 * gsf_46[k]
                  + f_3 * pc_x[k] * hsf_46[k];

        t_67[k] = f_12 * gsf_47[k]
                  + f_3 * pc_x[k] * hsf_47[k];

        t_68[k] = f_12 * gsf_48[k]
                  + f_3 * pc_x[k] * hsf_48[k];

        t_69[k] = f_12 * gsf_49[k]
                  + f_3 * pc_x[k] * hsf_49[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, pa_z, pc_y, pc_z, gsg0_25, gsf_16, gsf_28, gsg1_25, \
                         hsd0_29, hsd1_29, hsf_46, hsf_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = pa_z[k] * gsg0_25[k]
                  - f_6 * pc_z[k] * gsg1_25[k];

        t_71[k] = f_7 * gsf_16[k]
                  + f_3 * pc_z[k] * hsf_46[k];

        t_72[k] = f_7 * gsf_28[k]
                  + f_4 * hsd0_29[k]
                  - f_5 * hsd1_29[k]
                  + f_3 * pc_y[k] * hsf_48[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pa_y, pc_x, pc_y, gsg0_44, gsf_29, gsf_50, \
                         gsg1_44, hsd0_30, hsd1_30, hsf_49, hsf_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_7 * gsf_29[k]
                  + f_3 * pc_y[k] * hsf_49[k];

        t_74[k] = pa_y[k] * gsg0_44[k]
                  - f_6 * pc_y[k] * gsg1_44[k];

        t_75[k] = f_12 * gsf_50[k]
                  + f_1 * hsd0_30[k]
                  - f_2 * hsd1_30[k]
                  + f_3 * pc_x[k] * hsf_50[k];

        t_76[k] = f_3 * pc_y[k] * hsf_50[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, pc_y, pc_z, gsf_20, hsd0_30, hsd1_30, hsf_50, \
                         hsf_51, hsf_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_8 * gsf_20[k]
                  + f_3 * pc_z[k] * hsf_50[k];

        t_78[k] = f_4 * hsd0_30[k]
                  - f_5 * hsd1_30[k]
                  + f_3 * pc_y[k] * hsf_51[k];

        t_79[k] = f_3 * pc_y[k] * hsf_52[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pc_x, pc_y, gsf_55, gsf_56, gsf_57, hsd0_35, \
                         hsd1_35, hsf_55, hsf_56, hsf_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_12 * gsf_55[k]
                  + f_4 * hsd0_35[k]
                  - f_5 * hsd1_35[k]
                  + f_3 * pc_x[k] * hsf_55[k];

        t_81[k] = f_12 * gsf_56[k]
                  + f_3 * pc_x[k] * hsf_56[k];

        t_82[k] = f_12 * gsf_57[k]
                  + f_3 * pc_x[k] * hsf_57[k];

        t_83[k] = f_3 * pc_y[k] * hsf_55[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, pc_x, pc_y, gsf_59, hsd0_33, hsd0_34, hsd1_33, \
                         hsd1_34, hsf_56, hsf_57, hsf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_12 * gsf_59[k]
                  + f_3 * pc_x[k] * hsf_59[k];

        t_85[k] = f_1 * hsd0_33[k]
                  - f_2 * hsd1_33[k]
                  + f_3 * pc_y[k] * hsf_56[k];

        t_86[k] = f_10 * hsd0_34[k]
                  - f_11 * hsd1_34[k]
                  + f_3 * pc_y[k] * hsf_57[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pc_x, pc_y, pc_z, gsf_29, gsf_60, hsd0_35, \
                         hsd0_36, hsd1_35, hsd1_36, hsf_58, hsf_59, \
                         hsf_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_4 * hsd0_35[k]
                  - f_5 * hsd1_35[k]
                  + f_3 * pc_y[k] * hsf_58[k];

        t_88[k] = f_3 * pc_y[k] * hsf_59[k];

        t_89[k] = f_8 * gsf_29[k]
                  + f_1 * hsd0_35[k]
                  - f_2 * hsd1_35[k]
                  + f_3 * pc_z[k] * hsf_59[k];

        t_90[k] = f_8 * gsf_60[k]
                  + f_1 * hsd0_36[k]
                  - f_2 * hsd1_36[k]
                  + f_3 * pc_x[k] * hsf_60[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, pc_x, pc_y, pc_z, gsf_30, gsf_63, hsd0_39, \
                         hsd1_39, hsf_60, hsf_61, hsf_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_12 * gsf_30[k]
                  + f_3 * pc_y[k] * hsf_60[k];

        t_92[k] = f_3 * pc_z[k] * hsf_60[k];

        t_93[k] = f_8 * gsf_63[k]
                  + f_4 * hsd0_39[k]
                  - f_5 * hsd1_39[k]
                  + f_3 * pc_x[k] * hsf_63[k];

        t_94[k] = f_3 * pc_z[k] * hsf_61[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pc_x, pc_z, gsf_66, gsf_68, hsd0_36, hsd1_36, \
                         hsf_62, hsf_63, hsf_66, hsf_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_4 * hsd0_36[k]
                  - f_5 * hsd1_36[k]
                  + f_3 * pc_z[k] * hsf_62[k];

        t_96[k] = f_8 * gsf_66[k]
                  + f_3 * pc_x[k] * hsf_66[k];

        t_97[k] = f_3 * pc_z[k] * hsf_63[k];

        t_98[k] = f_8 * gsf_68[k]
                  + f_3 * pc_x[k] * hsf_68[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, pc_x, pc_y, pc_z, gsf_36, gsf_39, \
                         gsf_69, hsd0_39, hsd1_39, hsf_66, hsf_67, \
                         hsf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_8 * gsf_69[k]
                  + f_3 * pc_x[k] * hsf_69[k];

        t_100[k] = f_12 * gsf_36[k]
                   + f_1 * hsd0_39[k]
                   - f_2 * hsd1_39[k]
                   + f_3 * pc_y[k] * hsf_66[k];

        t_101[k] = f_3 * pc_z[k] * hsf_66[k];

        t_102[k] = f_4 * hsd0_39[k]
                   - f_5 * hsd1_39[k]
                   + f_3 * pc_z[k] * hsf_67[k];

        t_103[k] = f_12 * gsf_39[k]
                   + f_3 * pc_y[k] * hsf_69[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pa_z, pc_y, pc_z, gsg0_45, gsf_30, \
                         gsf_40, gsg1_45, hsd0_41, hsd1_41, hsf_69, \
                         hsf_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_1 * hsd0_41[k]
                   - f_2 * hsd1_41[k]
                   + f_3 * pc_z[k] * hsf_69[k];

        t_105[k] = pa_z[k] * gsg0_45[k]
                   - f_6 * pc_z[k] * gsg1_45[k];

        t_106[k] = f_8 * gsf_40[k]
                   + f_3 * pc_y[k] * hsf_70[k];

        t_107[k] = f_7 * gsf_30[k]
                   + f_3 * pc_z[k] * hsf_70[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, pa_z, pc_x, pc_y, pc_z, gsg0_48, gsf_42, gsf_75, \
                         gsg1_48, hsd0_47, hsd1_47, hsf_72, hsf_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = pa_z[k] * gsg0_48[k]
                   - f_6 * pc_z[k] * gsg1_48[k];

        t_109[k] = f_8 * gsf_42[k]
                   + f_3 * pc_y[k] * hsf_72[k];

        t_110[k] = f_8 * gsf_75[k]
                   + f_4 * hsd0_47[k]
                   - f_5 * hsd1_47[k]
                   + f_3 * pc_x[k] * hsf_75[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, pc_x, gsf_76, gsf_77, gsf_78, gsf_79, \
                         hsf_76, hsf_77, hsf_78, hsf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_8 * gsf_76[k]
                   + f_3 * pc_x[k] * hsf_76[k];

        t_112[k] = f_8 * gsf_77[k]
                   + f_3 * pc_x[k] * hsf_77[k];

        t_113[k] = f_8 * gsf_78[k]
                   + f_3 * pc_x[k] * hsf_78[k];

        t_114[k] = f_8 * gsf_79[k]
                   + f_3 * pc_x[k] * hsf_79[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, pa_z, pc_y, pc_z, gsg0_55, gsf_36, gsf_48, \
                         gsg1_55, hsd0_47, hsd1_47, hsf_76, hsf_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = pa_z[k] * gsg0_55[k]
                   - f_6 * pc_z[k] * gsg1_55[k];

        t_116[k] = f_7 * gsf_36[k]
                   + f_3 * pc_z[k] * hsf_76[k];

        t_117[k] = f_8 * gsf_48[k]
                   + f_4 * hsd0_47[k]
                   - f_5 * hsd1_47[k]
                   + f_3 * pc_y[k] * hsf_78[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, pa_y, pc_y, pc_z, gsg0_75, gsf_39, \
                         gsf_49, gsf_50, gsg1_75, hsd0_47, hsd1_47, hsf_79, \
                         hsf_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_8 * gsf_49[k]
                   + f_3 * pc_y[k] * hsf_79[k];

        t_119[k] = f_7 * gsf_39[k]
                   + f_1 * hsd0_47[k]
                   - f_2 * hsd1_47[k]
                   + f_3 * pc_z[k] * hsf_79[k];

        t_120[k] = pa_y[k] * gsg0_75[k]
                   - f_6 * pc_y[k] * gsg1_75[k];

        t_121[k] = f_7 * gsf_50[k]
                   + f_3 * pc_y[k] * hsf_80[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, pa_y, pc_y, pc_z, gsg0_78, gsg0_80, \
                         gsf_40, gsf_51, gsf_52, gsg1_78, gsg1_80, hsf_80, \
                         hsf_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_8 * gsf_40[k]
                   + f_3 * pc_z[k] * hsf_80[k];

        t_123[k] = pa_y[k] * gsg0_78[k]
                   + f_8 * gsf_51[k]
                   - f_6 * pc_y[k] * gsg1_78[k];

        t_124[k] = f_7 * gsf_52[k]
                   + f_3 * pc_y[k] * hsf_82[k];

        t_125[k] = pa_y[k] * gsg0_80[k]
                   - f_6 * pc_y[k] * gsg1_80[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, pc_x, gsf_86, gsf_87, gsf_88, gsf_89, \
                         hsf_86, hsf_87, hsf_88, hsf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_8 * gsf_86[k]
                   + f_3 * pc_x[k] * hsf_86[k];

        t_127[k] = f_8 * gsf_87[k]
                   + f_3 * pc_x[k] * hsf_87[k];

        t_128[k] = f_8 * gsf_88[k]
                   + f_3 * pc_x[k] * hsf_88[k];

        t_129[k] = f_8 * gsf_89[k]
                   + f_3 * pc_x[k] * hsf_89[k];
    }
}

static auto
compute_prim_hsg_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t gsg0,
                                                          const size_t gsf, const size_t gsg1,
                                                          const size_t hsd0, const size_t hsd1,
                                                          const size_t hsf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 2.0 / q;
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 1.5 / q;

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
    auto *t_253 = buffer.data(target + 253);
    auto *t_254 = buffer.data(target + 254);
    auto *t_255 = buffer.data(target + 255);
    auto *t_256 = buffer.data(target + 256);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *gsg0_89 = buffer.data(gsg0 + 89);
    const auto *gsg0_90 = buffer.data(gsg0 + 90);
    const auto *gsg0_93 = buffer.data(gsg0 + 93);
    const auto *gsg0_135 = buffer.data(gsg0 + 135);
    const auto *gsg0_140 = buffer.data(gsg0 + 140);
    const auto *gsg0_150 = buffer.data(gsg0 + 150);
    const auto *gsg0_151 = buffer.data(gsg0 + 151);
    const auto *gsg0_153 = buffer.data(gsg0 + 153);
    const auto *gsg0_160 = buffer.data(gsg0 + 160);
    const auto *gsg0_162 = buffer.data(gsg0 + 162);
    const auto *gsg0_164 = buffer.data(gsg0 + 164);
    const auto *gsg0_170 = buffer.data(gsg0 + 170);
    const auto *gsg0_175 = buffer.data(gsg0 + 175);
    const auto *gsg0_177 = buffer.data(gsg0 + 177);
    const auto *gsg0_179 = buffer.data(gsg0 + 179);
    const auto *gsg0_180 = buffer.data(gsg0 + 180);
    const auto *gsg0_183 = buffer.data(gsg0 + 183);
    const auto *gsg0_185 = buffer.data(gsg0 + 185);
    const auto *gsg0_190 = buffer.data(gsg0 + 190);
    const auto *gsg0_192 = buffer.data(gsg0 + 192);
    const auto *gsg0_194 = buffer.data(gsg0 + 194);
    const auto *gsg0_198 = buffer.data(gsg0 + 198);
    const auto *gsg0_205 = buffer.data(gsg0 + 205);
    const auto *gsg0_207 = buffer.data(gsg0 + 207);
    const auto *gsg0_209 = buffer.data(gsg0 + 209);
    const auto *gsg0_210 = buffer.data(gsg0 + 210);
    const auto *gsg0_215 = buffer.data(gsg0 + 215);
    const auto *gsg0_220 = buffer.data(gsg0 + 220);
    const auto *gsg0_221 = buffer.data(gsg0 + 221);
    const auto *gsg0_222 = buffer.data(gsg0 + 222);
    const auto *gsg0_224 = buffer.data(gsg0 + 224);

    const auto *gsf_46 = buffer.data(gsf + 46);
    const auto *gsf_50 = buffer.data(gsf + 50);
    const auto *gsf_56 = buffer.data(gsf + 56);
    const auto *gsf_58 = buffer.data(gsf + 58);
    const auto *gsf_59 = buffer.data(gsf + 59);
    const auto *gsf_60 = buffer.data(gsf + 60);
    const auto *gsf_66 = buffer.data(gsf + 66);
    const auto *gsf_69 = buffer.data(gsf + 69);
    const auto *gsf_70 = buffer.data(gsf + 70);
    const auto *gsf_72 = buffer.data(gsf + 72);
    const auto *gsf_76 = buffer.data(gsf + 76);
    const auto *gsf_79 = buffer.data(gsf + 79);
    const auto *gsf_80 = buffer.data(gsf + 80);
    const auto *gsf_82 = buffer.data(gsf + 82);
    const auto *gsf_86 = buffer.data(gsf + 86);
    const auto *gsf_89 = buffer.data(gsf + 89);
    const auto *gsf_90 = buffer.data(gsf + 90);
    const auto *gsf_92 = buffer.data(gsf + 92);
    const auto *gsf_95 = buffer.data(gsf + 95);
    const auto *gsf_96 = buffer.data(gsf + 96);
    const auto *gsf_97 = buffer.data(gsf + 97);
    const auto *gsf_99 = buffer.data(gsf + 99);
    const auto *gsf_100 = buffer.data(gsf + 100);
    const auto *gsf_103 = buffer.data(gsf + 103);
    const auto *gsf_106 = buffer.data(gsf + 106);
    const auto *gsf_107 = buffer.data(gsf + 107);
    const auto *gsf_108 = buffer.data(gsf + 108);
    const auto *gsf_109 = buffer.data(gsf + 109);
    const auto *gsf_115 = buffer.data(gsf + 115);
    const auto *gsf_116 = buffer.data(gsf + 116);
    const auto *gsf_117 = buffer.data(gsf + 117);
    const auto *gsf_118 = buffer.data(gsf + 118);
    const auto *gsf_119 = buffer.data(gsf + 119);
    const auto *gsf_120 = buffer.data(gsf + 120);
    const auto *gsf_123 = buffer.data(gsf + 123);
    const auto *gsf_125 = buffer.data(gsf + 125);
    const auto *gsf_126 = buffer.data(gsf + 126);
    const auto *gsf_127 = buffer.data(gsf + 127);
    const auto *gsf_128 = buffer.data(gsf + 128);
    const auto *gsf_129 = buffer.data(gsf + 129);
    const auto *gsf_133 = buffer.data(gsf + 133);
    const auto *gsf_136 = buffer.data(gsf + 136);
    const auto *gsf_137 = buffer.data(gsf + 137);
    const auto *gsf_138 = buffer.data(gsf + 138);
    const auto *gsf_139 = buffer.data(gsf + 139);
    const auto *gsf_140 = buffer.data(gsf + 140);
    const auto *gsf_145 = buffer.data(gsf + 145);
    const auto *gsf_146 = buffer.data(gsf + 146);
    const auto *gsf_147 = buffer.data(gsf + 147);
    const auto *gsf_149 = buffer.data(gsf + 149);

    const auto *gsg1_89 = buffer.data(gsg1 + 89);
    const auto *gsg1_90 = buffer.data(gsg1 + 90);
    const auto *gsg1_93 = buffer.data(gsg1 + 93);
    const auto *gsg1_135 = buffer.data(gsg1 + 135);
    const auto *gsg1_140 = buffer.data(gsg1 + 140);
    const auto *gsg1_150 = buffer.data(gsg1 + 150);
    const auto *gsg1_151 = buffer.data(gsg1 + 151);
    const auto *gsg1_153 = buffer.data(gsg1 + 153);
    const auto *gsg1_160 = buffer.data(gsg1 + 160);
    const auto *gsg1_162 = buffer.data(gsg1 + 162);
    const auto *gsg1_164 = buffer.data(gsg1 + 164);
    const auto *gsg1_170 = buffer.data(gsg1 + 170);
    const auto *gsg1_175 = buffer.data(gsg1 + 175);
    const auto *gsg1_177 = buffer.data(gsg1 + 177);
    const auto *gsg1_179 = buffer.data(gsg1 + 179);
    const auto *gsg1_180 = buffer.data(gsg1 + 180);
    const auto *gsg1_183 = buffer.data(gsg1 + 183);
    const auto *gsg1_185 = buffer.data(gsg1 + 185);
    const auto *gsg1_190 = buffer.data(gsg1 + 190);
    const auto *gsg1_192 = buffer.data(gsg1 + 192);
    const auto *gsg1_194 = buffer.data(gsg1 + 194);
    const auto *gsg1_198 = buffer.data(gsg1 + 198);
    const auto *gsg1_205 = buffer.data(gsg1 + 205);
    const auto *gsg1_207 = buffer.data(gsg1 + 207);
    const auto *gsg1_209 = buffer.data(gsg1 + 209);
    const auto *gsg1_210 = buffer.data(gsg1 + 210);
    const auto *gsg1_215 = buffer.data(gsg1 + 215);
    const auto *gsg1_220 = buffer.data(gsg1 + 220);
    const auto *gsg1_221 = buffer.data(gsg1 + 221);
    const auto *gsg1_222 = buffer.data(gsg1 + 222);
    const auto *gsg1_224 = buffer.data(gsg1 + 224);

    const auto *hsd0_51 = buffer.data(hsd0 + 51);
    const auto *hsd0_53 = buffer.data(hsd0 + 53);
    const auto *hsd0_54 = buffer.data(hsd0 + 54);
    const auto *hsd0_57 = buffer.data(hsd0 + 57);
    const auto *hsd0_58 = buffer.data(hsd0 + 58);
    const auto *hsd0_59 = buffer.data(hsd0 + 59);
    const auto *hsd0_60 = buffer.data(hsd0 + 60);
    const auto *hsd0_84 = buffer.data(hsd0 + 84);
    const auto *hsd0_90 = buffer.data(hsd0 + 90);
    const auto *hsd0_91 = buffer.data(hsd0 + 91);
    const auto *hsd0_93 = buffer.data(hsd0 + 93);
    const auto *hsd0_95 = buffer.data(hsd0 + 95);
    const auto *hsd0_98 = buffer.data(hsd0 + 98);
    const auto *hsd0_100 = buffer.data(hsd0 + 100);
    const auto *hsd0_101 = buffer.data(hsd0 + 101);
    const auto *hsd0_102 = buffer.data(hsd0 + 102);
    const auto *hsd0_103 = buffer.data(hsd0 + 103);

    const auto *hsd1_51 = buffer.data(hsd1 + 51);
    const auto *hsd1_53 = buffer.data(hsd1 + 53);
    const auto *hsd1_54 = buffer.data(hsd1 + 54);
    const auto *hsd1_57 = buffer.data(hsd1 + 57);
    const auto *hsd1_58 = buffer.data(hsd1 + 58);
    const auto *hsd1_59 = buffer.data(hsd1 + 59);
    const auto *hsd1_60 = buffer.data(hsd1 + 60);
    const auto *hsd1_84 = buffer.data(hsd1 + 84);
    const auto *hsd1_90 = buffer.data(hsd1 + 90);
    const auto *hsd1_91 = buffer.data(hsd1 + 91);
    const auto *hsd1_93 = buffer.data(hsd1 + 93);
    const auto *hsd1_95 = buffer.data(hsd1 + 95);
    const auto *hsd1_98 = buffer.data(hsd1 + 98);
    const auto *hsd1_100 = buffer.data(hsd1 + 100);
    const auto *hsd1_101 = buffer.data(hsd1 + 101);
    const auto *hsd1_102 = buffer.data(hsd1 + 102);
    const auto *hsd1_103 = buffer.data(hsd1 + 103);

    const auto *hsf_86 = buffer.data(hsf + 86);
    const auto *hsf_88 = buffer.data(hsf + 88);
    const auto *hsf_89 = buffer.data(hsf + 89);
    const auto *hsf_90 = buffer.data(hsf + 90);
    const auto *hsf_91 = buffer.data(hsf + 91);
    const auto *hsf_92 = buffer.data(hsf + 92);
    const auto *hsf_95 = buffer.data(hsf + 95);
    const auto *hsf_96 = buffer.data(hsf + 96);
    const auto *hsf_97 = buffer.data(hsf + 97);
    const auto *hsf_98 = buffer.data(hsf + 98);
    const auto *hsf_99 = buffer.data(hsf + 99);
    const auto *hsf_100 = buffer.data(hsf + 100);
    const auto *hsf_101 = buffer.data(hsf + 101);
    const auto *hsf_102 = buffer.data(hsf + 102);
    const auto *hsf_103 = buffer.data(hsf + 103);
    const auto *hsf_106 = buffer.data(hsf + 106);
    const auto *hsf_108 = buffer.data(hsf + 108);
    const auto *hsf_109 = buffer.data(hsf + 109);
    const auto *hsf_110 = buffer.data(hsf + 110);
    const auto *hsf_112 = buffer.data(hsf + 112);
    const auto *hsf_116 = buffer.data(hsf + 116);
    const auto *hsf_117 = buffer.data(hsf + 117);
    const auto *hsf_118 = buffer.data(hsf + 118);
    const auto *hsf_119 = buffer.data(hsf + 119);
    const auto *hsf_120 = buffer.data(hsf + 120);
    const auto *hsf_122 = buffer.data(hsf + 122);
    const auto *hsf_126 = buffer.data(hsf + 126);
    const auto *hsf_127 = buffer.data(hsf + 127);
    const auto *hsf_128 = buffer.data(hsf + 128);
    const auto *hsf_129 = buffer.data(hsf + 129);
    const auto *hsf_130 = buffer.data(hsf + 130);
    const auto *hsf_132 = buffer.data(hsf + 132);
    const auto *hsf_136 = buffer.data(hsf + 136);
    const auto *hsf_137 = buffer.data(hsf + 137);
    const auto *hsf_138 = buffer.data(hsf + 138);
    const auto *hsf_139 = buffer.data(hsf + 139);
    const auto *hsf_140 = buffer.data(hsf + 140);
    const auto *hsf_141 = buffer.data(hsf + 141);
    const auto *hsf_142 = buffer.data(hsf + 142);
    const auto *hsf_145 = buffer.data(hsf + 145);
    const auto *hsf_146 = buffer.data(hsf + 146);
    const auto *hsf_147 = buffer.data(hsf + 147);
    const auto *hsf_149 = buffer.data(hsf + 149);
    const auto *hsf_150 = buffer.data(hsf + 150);
    const auto *hsf_151 = buffer.data(hsf + 151);
    const auto *hsf_153 = buffer.data(hsf + 153);
    const auto *hsf_155 = buffer.data(hsf + 155);
    const auto *hsf_156 = buffer.data(hsf + 156);
    const auto *hsf_157 = buffer.data(hsf + 157);
    const auto *hsf_158 = buffer.data(hsf + 158);
    const auto *hsf_159 = buffer.data(hsf + 159);
    const auto *hsf_162 = buffer.data(hsf + 162);
    const auto *hsf_164 = buffer.data(hsf + 164);
    const auto *hsf_165 = buffer.data(hsf + 165);
    const auto *hsf_166 = buffer.data(hsf + 166);
    const auto *hsf_167 = buffer.data(hsf + 167);
    const auto *hsf_168 = buffer.data(hsf + 168);
    const auto *hsf_169 = buffer.data(hsf + 169);
    const auto *hsf_170 = buffer.data(hsf + 170);
    const auto *hsf_171 = buffer.data(hsf + 171);

#pragma omp simd aligned(t_130, t_131, t_132, pc_y, pc_z, gsf_46, gsf_56, gsf_58, hsd0_51, \
                         hsd0_53, hsd1_51, hsd1_53, hsf_86, hsf_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_7 * gsf_56[k]
                   + f_1 * hsd0_51[k]
                   - f_2 * hsd1_51[k]
                   + f_3 * pc_y[k] * hsf_86[k];

        t_131[k] = f_8 * gsf_46[k]
                   + f_3 * pc_z[k] * hsf_86[k];

        t_132[k] = f_7 * gsf_58[k]
                   + f_4 * hsd0_53[k]
                   - f_5 * hsd1_53[k]
                   + f_3 * pc_y[k] * hsf_88[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, pa_y, pc_x, pc_y, gsg0_89, gsf_59, \
                         gsf_90, gsg1_89, hsd0_54, hsd1_54, hsf_89, \
                         hsf_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_7 * gsf_59[k]
                   + f_3 * pc_y[k] * hsf_89[k];

        t_134[k] = pa_y[k] * gsg0_89[k]
                   - f_6 * pc_y[k] * gsg1_89[k];

        t_135[k] = f_8 * gsf_90[k]
                   + f_1 * hsd0_54[k]
                   - f_2 * hsd1_54[k]
                   + f_3 * pc_x[k] * hsf_90[k];

        t_136[k] = f_3 * pc_y[k] * hsf_90[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, pc_y, pc_z, gsf_50, hsd0_54, hsd1_54, hsf_90, \
                         hsf_91, hsf_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_12 * gsf_50[k]
                   + f_3 * pc_z[k] * hsf_90[k];

        t_138[k] = f_4 * hsd0_54[k]
                   - f_5 * hsd1_54[k]
                   + f_3 * pc_y[k] * hsf_91[k];

        t_139[k] = f_3 * pc_y[k] * hsf_92[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, pc_x, pc_y, gsf_95, gsf_96, gsf_97, \
                         hsd0_59, hsd1_59, hsf_95, hsf_96, hsf_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_8 * gsf_95[k]
                   + f_4 * hsd0_59[k]
                   - f_5 * hsd1_59[k]
                   + f_3 * pc_x[k] * hsf_95[k];

        t_141[k] = f_8 * gsf_96[k]
                   + f_3 * pc_x[k] * hsf_96[k];

        t_142[k] = f_8 * gsf_97[k]
                   + f_3 * pc_x[k] * hsf_97[k];

        t_143[k] = f_3 * pc_y[k] * hsf_95[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, pc_x, pc_y, gsf_99, hsd0_57, hsd0_58, hsd1_57, \
                         hsd1_58, hsf_96, hsf_97, hsf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_8 * gsf_99[k]
                   + f_3 * pc_x[k] * hsf_99[k];

        t_145[k] = f_1 * hsd0_57[k]
                   - f_2 * hsd1_57[k]
                   + f_3 * pc_y[k] * hsf_96[k];

        t_146[k] = f_10 * hsd0_58[k]
                   - f_11 * hsd1_58[k]
                   + f_3 * pc_y[k] * hsf_97[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, pa_x, pc_x, pc_y, pc_z, gsg0_150, gsf_59, \
                         gsf_100, gsg1_150, hsd0_59, hsd1_59, hsf_98, \
                         hsf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = f_4 * hsd0_59[k]
                   - f_5 * hsd1_59[k]
                   + f_3 * pc_y[k] * hsf_98[k];

        t_148[k] = f_3 * pc_y[k] * hsf_99[k];

        t_149[k] = f_12 * gsf_59[k]
                   + f_1 * hsd0_59[k]
                   - f_2 * hsd1_59[k]
                   + f_3 * pc_z[k] * hsf_99[k];

        t_150[k] = pa_x[k] * gsg0_150[k]
                   + f_9 * gsf_100[k]
                   - f_6 * pc_x[k] * gsg1_150[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, pa_x, pc_x, pc_y, pc_z, gsg0_153, gsf_60, \
                         gsf_103, gsg1_153, hsf_100, hsf_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_9 * gsf_60[k]
                   + f_3 * pc_y[k] * hsf_100[k];

        t_152[k] = f_3 * pc_z[k] * hsf_100[k];

        t_153[k] = pa_x[k] * gsg0_153[k]
                   + f_8 * gsf_103[k]
                   - f_6 * pc_x[k] * gsg1_153[k];

        t_154[k] = f_3 * pc_z[k] * hsf_101[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, pc_x, pc_z, gsf_106, gsf_108, hsd0_60, \
                         hsd1_60, hsf_102, hsf_103, hsf_106, hsf_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_4 * hsd0_60[k]
                   - f_5 * hsd1_60[k]
                   + f_3 * pc_z[k] * hsf_102[k];

        t_156[k] = f_7 * gsf_106[k]
                   + f_3 * pc_x[k] * hsf_106[k];

        t_157[k] = f_3 * pc_z[k] * hsf_103[k];

        t_158[k] = f_7 * gsf_108[k]
                   + f_3 * pc_x[k] * hsf_108[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, pa_x, pc_x, pc_z, gsg0_160, gsg0_162, \
                         gsf_109, gsg1_160, gsg1_162, hsf_106, \
                         hsf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_7 * gsf_109[k]
                   + f_3 * pc_x[k] * hsf_109[k];

        t_160[k] = pa_x[k] * gsg0_160[k]
                   - f_6 * pc_x[k] * gsg1_160[k];

        t_161[k] = f_3 * pc_z[k] * hsf_106[k];

        t_162[k] = pa_x[k] * gsg0_162[k]
                   - f_6 * pc_x[k] * gsg1_162[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, pa_x, pa_z, pc_x, pc_y, pc_z, gsg0_90, gsg0_164, \
                         gsf_69, gsg1_90, gsg1_164, hsf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_9 * gsf_69[k]
                   + f_3 * pc_y[k] * hsf_109[k];

        t_164[k] = pa_x[k] * gsg0_164[k]
                   - f_6 * pc_x[k] * gsg1_164[k];

        t_165[k] = pa_z[k] * gsg0_90[k]
                   - f_6 * pc_z[k] * gsg1_90[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, pa_z, pc_y, pc_z, gsg0_93, gsf_60, \
                         gsf_70, gsf_72, gsg1_93, hsf_110, hsf_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_12 * gsf_70[k]
                   + f_3 * pc_y[k] * hsf_110[k];

        t_167[k] = f_7 * gsf_60[k]
                   + f_3 * pc_z[k] * hsf_110[k];

        t_168[k] = pa_z[k] * gsg0_93[k]
                   - f_6 * pc_z[k] * gsg1_93[k];

        t_169[k] = f_12 * gsf_72[k]
                   + f_3 * pc_y[k] * hsf_112[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, pa_x, pc_x, gsg0_170, gsf_115, gsf_116, \
                         gsf_117, gsf_118, gsg1_170, hsf_116, hsf_117, \
                         hsf_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = pa_x[k] * gsg0_170[k]
                   + f_8 * gsf_115[k]
                   - f_6 * pc_x[k] * gsg1_170[k];

        t_171[k] = f_7 * gsf_116[k]
                   + f_3 * pc_x[k] * hsf_116[k];

        t_172[k] = f_7 * gsf_117[k]
                   + f_3 * pc_x[k] * hsf_117[k];

        t_173[k] = f_7 * gsf_118[k]
                   + f_3 * pc_x[k] * hsf_118[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, pa_x, pc_x, pc_z, gsg0_175, gsg0_177, \
                         gsf_66, gsf_119, gsg1_175, gsg1_177, hsf_116, \
                         hsf_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_7 * gsf_119[k]
                   + f_3 * pc_x[k] * hsf_119[k];

        t_175[k] = pa_x[k] * gsg0_175[k]
                   - f_6 * pc_x[k] * gsg1_175[k];

        t_176[k] = f_7 * gsf_66[k]
                   + f_3 * pc_z[k] * hsf_116[k];

        t_177[k] = pa_x[k] * gsg0_177[k]
                   - f_6 * pc_x[k] * gsg1_177[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, t_181, pa_x, pc_x, pc_y, gsg0_179, gsg0_180, \
                         gsf_79, gsf_80, gsf_120, gsg1_179, gsg1_180, hsf_119, \
                         hsf_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_12 * gsf_79[k]
                   + f_3 * pc_y[k] * hsf_119[k];

        t_179[k] = pa_x[k] * gsg0_179[k]
                   - f_6 * pc_x[k] * gsg1_179[k];

        t_180[k] = pa_x[k] * gsg0_180[k]
                   + f_9 * gsf_120[k]
                   - f_6 * pc_x[k] * gsg1_180[k];

        t_181[k] = f_8 * gsf_80[k]
                   + f_3 * pc_y[k] * hsf_120[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, pa_x, pc_x, pc_y, pc_z, gsg0_183, gsf_70, \
                         gsf_82, gsf_123, gsg1_183, hsf_120, hsf_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = f_8 * gsf_70[k]
                   + f_3 * pc_z[k] * hsf_120[k];

        t_183[k] = pa_x[k] * gsg0_183[k]
                   + f_8 * gsf_123[k]
                   - f_6 * pc_x[k] * gsg1_183[k];

        t_184[k] = f_8 * gsf_82[k]
                   + f_3 * pc_y[k] * hsf_122[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pa_x, pc_x, gsg0_185, gsf_125, gsf_126, \
                         gsf_127, gsf_128, gsg1_185, hsf_126, hsf_127, \
                         hsf_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = pa_x[k] * gsg0_185[k]
                   + f_8 * gsf_125[k]
                   - f_6 * pc_x[k] * gsg1_185[k];

        t_186[k] = f_7 * gsf_126[k]
                   + f_3 * pc_x[k] * hsf_126[k];

        t_187[k] = f_7 * gsf_127[k]
                   + f_3 * pc_x[k] * hsf_127[k];

        t_188[k] = f_7 * gsf_128[k]
                   + f_3 * pc_x[k] * hsf_128[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, pa_x, pc_x, pc_z, gsg0_190, gsg0_192, \
                         gsf_76, gsf_129, gsg1_190, gsg1_192, hsf_126, \
                         hsf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_7 * gsf_129[k]
                   + f_3 * pc_x[k] * hsf_129[k];

        t_190[k] = pa_x[k] * gsg0_190[k]
                   - f_6 * pc_x[k] * gsg1_190[k];

        t_191[k] = f_8 * gsf_76[k]
                   + f_3 * pc_z[k] * hsf_126[k];

        t_192[k] = pa_x[k] * gsg0_192[k]
                   - f_6 * pc_x[k] * gsg1_192[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, pa_x, pa_y, pc_x, pc_y, gsg0_135, \
                         gsg0_194, gsf_89, gsf_90, gsg1_135, gsg1_194, hsf_129, \
                         hsf_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_8 * gsf_89[k]
                   + f_3 * pc_y[k] * hsf_129[k];

        t_194[k] = pa_x[k] * gsg0_194[k]
                   - f_6 * pc_x[k] * gsg1_194[k];

        t_195[k] = pa_y[k] * gsg0_135[k]
                   - f_6 * pc_y[k] * gsg1_135[k];

        t_196[k] = f_7 * gsf_90[k]
                   + f_3 * pc_y[k] * hsf_130[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, pa_x, pc_x, pc_y, pc_z, gsg0_198, gsf_80, \
                         gsf_92, gsf_133, gsg1_198, hsf_130, hsf_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_12 * gsf_80[k]
                   + f_3 * pc_z[k] * hsf_130[k];

        t_198[k] = pa_x[k] * gsg0_198[k]
                   + f_8 * gsf_133[k]
                   - f_6 * pc_x[k] * gsg1_198[k];

        t_199[k] = f_7 * gsf_92[k]
                   + f_3 * pc_y[k] * hsf_132[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, pa_y, pc_x, pc_y, gsg0_140, gsf_136, \
                         gsf_137, gsf_138, gsg1_140, hsf_136, hsf_137, \
                         hsf_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = pa_y[k] * gsg0_140[k]
                   - f_6 * pc_y[k] * gsg1_140[k];

        t_201[k] = f_7 * gsf_136[k]
                   + f_3 * pc_x[k] * hsf_136[k];

        t_202[k] = f_7 * gsf_137[k]
                   + f_3 * pc_x[k] * hsf_137[k];

        t_203[k] = f_7 * gsf_138[k]
                   + f_3 * pc_x[k] * hsf_138[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, pa_x, pc_x, pc_z, gsg0_205, gsg0_207, \
                         gsf_86, gsf_139, gsg1_205, gsg1_207, hsf_136, \
                         hsf_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_7 * gsf_139[k]
                   + f_3 * pc_x[k] * hsf_139[k];

        t_205[k] = pa_x[k] * gsg0_205[k]
                   - f_6 * pc_x[k] * gsg1_205[k];

        t_206[k] = f_12 * gsf_86[k]
                   + f_3 * pc_z[k] * hsf_136[k];

        t_207[k] = pa_x[k] * gsg0_207[k]
                   - f_6 * pc_x[k] * gsg1_207[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, pa_x, pc_x, pc_y, gsg0_209, gsg0_210, \
                         gsf_99, gsf_140, gsg1_209, gsg1_210, hsf_139, \
                         hsf_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_7 * gsf_99[k]
                   + f_3 * pc_y[k] * hsf_139[k];

        t_209[k] = pa_x[k] * gsg0_209[k]
                   - f_6 * pc_x[k] * gsg1_209[k];

        t_210[k] = pa_x[k] * gsg0_210[k]
                   + f_9 * gsf_140[k]
                   - f_6 * pc_x[k] * gsg1_210[k];

        t_211[k] = f_3 * pc_y[k] * hsf_140[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, pc_y, pc_z, gsf_90, hsd0_84, hsd1_84, hsf_140, \
                         hsf_141, hsf_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_9 * gsf_90[k]
                   + f_3 * pc_z[k] * hsf_140[k];

        t_213[k] = f_4 * hsd0_84[k]
                   - f_5 * hsd1_84[k]
                   + f_3 * pc_y[k] * hsf_141[k];

        t_214[k] = f_3 * pc_y[k] * hsf_142[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, pa_x, pc_x, pc_y, gsg0_215, gsf_145, \
                         gsf_146, gsf_147, gsg1_215, hsf_145, hsf_146, \
                         hsf_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = pa_x[k] * gsg0_215[k]
                   + f_8 * gsf_145[k]
                   - f_6 * pc_x[k] * gsg1_215[k];

        t_216[k] = f_7 * gsf_146[k]
                   + f_3 * pc_x[k] * hsf_146[k];

        t_217[k] = f_7 * gsf_147[k]
                   + f_3 * pc_x[k] * hsf_147[k];

        t_218[k] = f_3 * pc_y[k] * hsf_145[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, t_222, t_223, pa_x, pc_x, pc_y, gsg0_220, \
                         gsg0_221, gsg0_222, gsf_149, gsg1_220, gsg1_221, gsg1_222, \
                         hsf_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = f_7 * gsf_149[k]
                   + f_3 * pc_x[k] * hsf_149[k];

        t_220[k] = pa_x[k] * gsg0_220[k]
                   - f_6 * pc_x[k] * gsg1_220[k];

        t_221[k] = pa_x[k] * gsg0_221[k]
                   - f_6 * pc_x[k] * gsg1_221[k];

        t_222[k] = pa_x[k] * gsg0_222[k]
                   - f_6 * pc_x[k] * gsg1_222[k];

        t_223[k] = f_3 * pc_y[k] * hsf_149[k];
    }

#pragma omp simd aligned(t_224, t_225, t_226, t_227, pa_x, pc_x, pc_z, gsg0_224, gsg1_224, \
                         hsd0_90, hsd0_91, hsd1_90, hsd1_91, hsf_150, \
                         hsf_151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_224[k] = pa_x[k] * gsg0_224[k]
                   - f_6 * pc_x[k] * gsg1_224[k];

        t_225[k] = f_1 * hsd0_90[k]
                   - f_2 * hsd1_90[k]
                   + f_3 * pc_x[k] * hsf_150[k];

        t_226[k] = f_10 * hsd0_91[k]
                   - f_11 * hsd1_91[k]
                   + f_3 * pc_x[k] * hsf_151[k];

        t_227[k] = f_3 * pc_z[k] * hsf_150[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, t_231, t_232, pc_x, pc_z, hsd0_93, hsd0_95, \
                         hsd1_93, hsd1_95, hsf_151, hsf_153, hsf_155, hsf_156, \
                         hsf_157 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_4 * hsd0_93[k]
                   - f_5 * hsd1_93[k]
                   + f_3 * pc_x[k] * hsf_153[k];

        t_229[k] = f_3 * pc_z[k] * hsf_151[k];

        t_230[k] = f_4 * hsd0_95[k]
                   - f_5 * hsd1_95[k]
                   + f_3 * pc_x[k] * hsf_155[k];

        t_231[k] = f_3 * pc_x[k] * hsf_156[k];

        t_232[k] = f_3 * pc_x[k] * hsf_157[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, t_236, t_237, pc_x, pc_y, pc_z, gsf_106, \
                         hsd0_93, hsd1_93, hsf_156, hsf_157, hsf_158, \
                         hsf_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = f_3 * pc_x[k] * hsf_158[k];

        t_234[k] = f_3 * pc_x[k] * hsf_159[k];

        t_235[k] = f_0 * gsf_106[k]
                   + f_1 * hsd0_93[k]
                   - f_2 * hsd1_93[k]
                   + f_3 * pc_y[k] * hsf_156[k];

        t_236[k] = f_3 * pc_z[k] * hsf_156[k];

        t_237[k] = f_4 * hsd0_93[k]
                   - f_5 * hsd1_93[k]
                   + f_3 * pc_z[k] * hsf_157[k];
    }

#pragma omp simd aligned(t_238, t_239, t_240, t_241, pa_z, pc_y, pc_z, gsg0_150, gsg0_151, \
                         gsf_109, gsg1_150, gsg1_151, hsd0_95, hsd1_95, \
                         hsf_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_238[k] = f_0 * gsf_109[k]
                   + f_3 * pc_y[k] * hsf_159[k];

        t_239[k] = f_1 * hsd0_95[k]
                   - f_2 * hsd1_95[k]
                   + f_3 * pc_z[k] * hsf_159[k];

        t_240[k] = pa_z[k] * gsg0_150[k]
                   - f_6 * pc_z[k] * gsg1_150[k];

        t_241[k] = pa_z[k] * gsg0_151[k]
                   - f_6 * pc_z[k] * gsg1_151[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, pa_z, pc_x, pc_z, gsg0_153, gsg1_153, hsd0_98, \
                         hsd0_100, hsd1_98, hsd1_100, hsf_162, \
                         hsf_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = f_10 * hsd0_98[k]
                   - f_11 * hsd1_98[k]
                   + f_3 * pc_x[k] * hsf_162[k];

        t_243[k] = pa_z[k] * gsg0_153[k]
                   - f_6 * pc_z[k] * gsg1_153[k];

        t_244[k] = f_4 * hsd0_100[k]
                   - f_5 * hsd1_100[k]
                   + f_3 * pc_x[k] * hsf_164[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, pc_x, hsd0_101, hsd1_101, hsf_165, \
                         hsf_166, hsf_167, hsf_168, hsf_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_4 * hsd0_101[k]
                   - f_5 * hsd1_101[k]
                   + f_3 * pc_x[k] * hsf_165[k];

        t_246[k] = f_3 * pc_x[k] * hsf_166[k];

        t_247[k] = f_3 * pc_x[k] * hsf_167[k];

        t_248[k] = f_3 * pc_x[k] * hsf_168[k];

        t_249[k] = f_3 * pc_x[k] * hsf_169[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, pa_z, pc_y, pc_z, gsg0_160, gsg0_162, \
                         gsf_106, gsf_107, gsf_119, gsg1_160, gsg1_162, hsf_166, \
                         hsf_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = pa_z[k] * gsg0_160[k]
                   - f_6 * pc_z[k] * gsg1_160[k];

        t_251[k] = f_7 * gsf_106[k]
                   + f_3 * pc_z[k] * hsf_166[k];

        t_252[k] = pa_z[k] * gsg0_162[k]
                   + f_8 * gsf_107[k]
                   - f_6 * pc_z[k] * gsg1_162[k];

        t_253[k] = f_9 * gsf_119[k]
                   + f_3 * pc_y[k] * hsf_169[k];
    }

#pragma omp simd aligned(t_254, t_255, t_256, pc_x, pc_z, gsf_109, hsd0_101, hsd0_102, \
                         hsd0_103, hsd1_101, hsd1_102, hsd1_103, hsf_169, hsf_170, \
                         hsf_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_254[k] = f_7 * gsf_109[k]
                   + f_1 * hsd0_101[k]
                   - f_2 * hsd1_101[k]
                   + f_3 * pc_z[k] * hsf_169[k];

        t_255[k] = f_1 * hsd0_102[k]
                   - f_2 * hsd1_102[k]
                   + f_3 * pc_x[k] * hsf_170[k];

        t_256[k] = f_10 * hsd0_103[k]
                   - f_11 * hsd1_103[k]
                   + f_3 * pc_x[k] * hsf_171[k];
    }
}

static auto
compute_prim_hsg_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t gsg0,
                                                          const size_t gsf, const size_t gsg1,
                                                          const size_t hsd0, const size_t hsd1,
                                                          const size_t hsf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 2.0 / q;
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 1.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *gsg0_210 = buffer.data(gsg0 + 210);
    const auto *gsg0_212 = buffer.data(gsg0 + 212);
    const auto *gsg0_215 = buffer.data(gsg0 + 215);
    const auto *gsg0_220 = buffer.data(gsg0 + 220);
    const auto *gsg0_222 = buffer.data(gsg0 + 222);
    const auto *gsg0_224 = buffer.data(gsg0 + 224);

    const auto *gsf_116 = buffer.data(gsf + 116);
    const auto *gsf_119 = buffer.data(gsf + 119);
    const auto *gsf_126 = buffer.data(gsf + 126);
    const auto *gsf_128 = buffer.data(gsf + 128);
    const auto *gsf_129 = buffer.data(gsf + 129);
    const auto *gsf_136 = buffer.data(gsf + 136);
    const auto *gsf_138 = buffer.data(gsf + 138);
    const auto *gsf_139 = buffer.data(gsf + 139);
    const auto *gsf_146 = buffer.data(gsf + 146);
    const auto *gsf_148 = buffer.data(gsf + 148);
    const auto *gsf_149 = buffer.data(gsf + 149);

    const auto *gsg1_210 = buffer.data(gsg1 + 210);
    const auto *gsg1_212 = buffer.data(gsg1 + 212);
    const auto *gsg1_215 = buffer.data(gsg1 + 215);
    const auto *gsg1_220 = buffer.data(gsg1 + 220);
    const auto *gsg1_222 = buffer.data(gsg1 + 222);
    const auto *gsg1_224 = buffer.data(gsg1 + 224);

    const auto *hsd0_104 = buffer.data(hsd0 + 104);
    const auto *hsd0_105 = buffer.data(hsd0 + 105);
    const auto *hsd0_106 = buffer.data(hsd0 + 106);
    const auto *hsd0_107 = buffer.data(hsd0 + 107);
    const auto *hsd0_108 = buffer.data(hsd0 + 108);
    const auto *hsd0_109 = buffer.data(hsd0 + 109);
    const auto *hsd0_110 = buffer.data(hsd0 + 110);
    const auto *hsd0_111 = buffer.data(hsd0 + 111);
    const auto *hsd0_112 = buffer.data(hsd0 + 112);
    const auto *hsd0_113 = buffer.data(hsd0 + 113);
    const auto *hsd0_115 = buffer.data(hsd0 + 115);
    const auto *hsd0_117 = buffer.data(hsd0 + 117);
    const auto *hsd0_118 = buffer.data(hsd0 + 118);
    const auto *hsd0_120 = buffer.data(hsd0 + 120);
    const auto *hsd0_122 = buffer.data(hsd0 + 122);
    const auto *hsd0_123 = buffer.data(hsd0 + 123);
    const auto *hsd0_124 = buffer.data(hsd0 + 124);
    const auto *hsd0_125 = buffer.data(hsd0 + 125);

    const auto *hsd1_104 = buffer.data(hsd1 + 104);
    const auto *hsd1_105 = buffer.data(hsd1 + 105);
    const auto *hsd1_106 = buffer.data(hsd1 + 106);
    const auto *hsd1_107 = buffer.data(hsd1 + 107);
    const auto *hsd1_108 = buffer.data(hsd1 + 108);
    const auto *hsd1_109 = buffer.data(hsd1 + 109);
    const auto *hsd1_110 = buffer.data(hsd1 + 110);
    const auto *hsd1_111 = buffer.data(hsd1 + 111);
    const auto *hsd1_112 = buffer.data(hsd1 + 112);
    const auto *hsd1_113 = buffer.data(hsd1 + 113);
    const auto *hsd1_115 = buffer.data(hsd1 + 115);
    const auto *hsd1_117 = buffer.data(hsd1 + 117);
    const auto *hsd1_118 = buffer.data(hsd1 + 118);
    const auto *hsd1_120 = buffer.data(hsd1 + 120);
    const auto *hsd1_122 = buffer.data(hsd1 + 122);
    const auto *hsd1_123 = buffer.data(hsd1 + 123);
    const auto *hsd1_124 = buffer.data(hsd1 + 124);
    const auto *hsd1_125 = buffer.data(hsd1 + 125);

    const auto *hsf_172 = buffer.data(hsf + 172);
    const auto *hsf_173 = buffer.data(hsf + 173);
    const auto *hsf_174 = buffer.data(hsf + 174);
    const auto *hsf_175 = buffer.data(hsf + 175);
    const auto *hsf_176 = buffer.data(hsf + 176);
    const auto *hsf_177 = buffer.data(hsf + 177);
    const auto *hsf_178 = buffer.data(hsf + 178);
    const auto *hsf_179 = buffer.data(hsf + 179);
    const auto *hsf_180 = buffer.data(hsf + 180);
    const auto *hsf_181 = buffer.data(hsf + 181);
    const auto *hsf_182 = buffer.data(hsf + 182);
    const auto *hsf_183 = buffer.data(hsf + 183);
    const auto *hsf_184 = buffer.data(hsf + 184);
    const auto *hsf_185 = buffer.data(hsf + 185);
    const auto *hsf_186 = buffer.data(hsf + 186);
    const auto *hsf_187 = buffer.data(hsf + 187);
    const auto *hsf_188 = buffer.data(hsf + 188);
    const auto *hsf_189 = buffer.data(hsf + 189);
    const auto *hsf_191 = buffer.data(hsf + 191);
    const auto *hsf_193 = buffer.data(hsf + 193);
    const auto *hsf_194 = buffer.data(hsf + 194);
    const auto *hsf_196 = buffer.data(hsf + 196);
    const auto *hsf_197 = buffer.data(hsf + 197);
    const auto *hsf_198 = buffer.data(hsf + 198);
    const auto *hsf_199 = buffer.data(hsf + 199);
    const auto *hsf_200 = buffer.data(hsf + 200);
    const auto *hsf_202 = buffer.data(hsf + 202);
    const auto *hsf_203 = buffer.data(hsf + 203);
    const auto *hsf_205 = buffer.data(hsf + 205);
    const auto *hsf_206 = buffer.data(hsf + 206);
    const auto *hsf_207 = buffer.data(hsf + 207);
    const auto *hsf_208 = buffer.data(hsf + 208);
    const auto *hsf_209 = buffer.data(hsf + 209);

#pragma omp simd aligned(t_257, t_258, t_259, pc_x, hsd0_104, hsd0_105, hsd0_106, hsd1_104, \
                         hsd1_105, hsd1_106, hsf_172, hsf_173, \
                         hsf_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_10 * hsd0_104[k]
                   - f_11 * hsd1_104[k]
                   + f_3 * pc_x[k] * hsf_172[k];

        t_258[k] = f_4 * hsd0_105[k]
                   - f_5 * hsd1_105[k]
                   + f_3 * pc_x[k] * hsf_173[k];

        t_259[k] = f_4 * hsd0_106[k]
                   - f_5 * hsd1_106[k]
                   + f_3 * pc_x[k] * hsf_174[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, pc_x, hsd0_107, hsd1_107, hsf_175, \
                         hsf_176, hsf_177, hsf_178, hsf_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_4 * hsd0_107[k]
                   - f_5 * hsd1_107[k]
                   + f_3 * pc_x[k] * hsf_175[k];

        t_261[k] = f_3 * pc_x[k] * hsf_176[k];

        t_262[k] = f_3 * pc_x[k] * hsf_177[k];

        t_263[k] = f_3 * pc_x[k] * hsf_178[k];

        t_264[k] = f_3 * pc_x[k] * hsf_179[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, pc_y, pc_z, gsf_116, gsf_126, gsf_128, hsd0_105, \
                         hsd0_107, hsd1_105, hsd1_107, hsf_176, \
                         hsf_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = f_12 * gsf_126[k]
                   + f_1 * hsd0_105[k]
                   - f_2 * hsd1_105[k]
                   + f_3 * pc_y[k] * hsf_176[k];

        t_266[k] = f_8 * gsf_116[k]
                   + f_3 * pc_z[k] * hsf_176[k];

        t_267[k] = f_12 * gsf_128[k]
                   + f_4 * hsd0_107[k]
                   - f_5 * hsd1_107[k]
                   + f_3 * pc_y[k] * hsf_178[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, pc_x, pc_y, pc_z, gsf_119, gsf_129, hsd0_107, \
                         hsd0_108, hsd1_107, hsd1_108, hsf_179, \
                         hsf_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_12 * gsf_129[k]
                   + f_3 * pc_y[k] * hsf_179[k];

        t_269[k] = f_8 * gsf_119[k]
                   + f_1 * hsd0_107[k]
                   - f_2 * hsd1_107[k]
                   + f_3 * pc_z[k] * hsf_179[k];

        t_270[k] = f_1 * hsd0_108[k]
                   - f_2 * hsd1_108[k]
                   + f_3 * pc_x[k] * hsf_180[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, pc_x, hsd0_109, hsd0_110, hsd0_111, hsd1_109, \
                         hsd1_110, hsd1_111, hsf_181, hsf_182, \
                         hsf_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = f_10 * hsd0_109[k]
                   - f_11 * hsd1_109[k]
                   + f_3 * pc_x[k] * hsf_181[k];

        t_272[k] = f_10 * hsd0_110[k]
                   - f_11 * hsd1_110[k]
                   + f_3 * pc_x[k] * hsf_182[k];

        t_273[k] = f_4 * hsd0_111[k]
                   - f_5 * hsd1_111[k]
                   + f_3 * pc_x[k] * hsf_183[k];
    }

#pragma omp simd aligned(t_274, t_275, t_276, t_277, t_278, pc_x, hsd0_112, hsd0_113, \
                         hsd1_112, hsd1_113, hsf_184, hsf_185, hsf_186, hsf_187, \
                         hsf_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_274[k] = f_4 * hsd0_112[k]
                   - f_5 * hsd1_112[k]
                   + f_3 * pc_x[k] * hsf_184[k];

        t_275[k] = f_4 * hsd0_113[k]
                   - f_5 * hsd1_113[k]
                   + f_3 * pc_x[k] * hsf_185[k];

        t_276[k] = f_3 * pc_x[k] * hsf_186[k];

        t_277[k] = f_3 * pc_x[k] * hsf_187[k];

        t_278[k] = f_3 * pc_x[k] * hsf_188[k];
    }

#pragma omp simd aligned(t_279, t_280, t_281, pc_x, pc_y, pc_z, gsf_126, gsf_136, hsd0_111, \
                         hsd1_111, hsf_186, hsf_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_279[k] = f_3 * pc_x[k] * hsf_189[k];

        t_280[k] = f_8 * gsf_136[k]
                   + f_1 * hsd0_111[k]
                   - f_2 * hsd1_111[k]
                   + f_3 * pc_y[k] * hsf_186[k];

        t_281[k] = f_12 * gsf_126[k]
                   + f_3 * pc_z[k] * hsf_186[k];
    }

#pragma omp simd aligned(t_282, t_283, t_284, t_285, pa_y, pc_y, pc_z, gsg0_210, gsf_129, \
                         gsf_138, gsf_139, gsg1_210, hsd0_113, hsd1_113, hsf_188, \
                         hsf_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_282[k] = f_8 * gsf_138[k]
                   + f_4 * hsd0_113[k]
                   - f_5 * hsd1_113[k]
                   + f_3 * pc_y[k] * hsf_188[k];

        t_283[k] = f_8 * gsf_139[k]
                   + f_3 * pc_y[k] * hsf_189[k];

        t_284[k] = f_12 * gsf_129[k]
                   + f_1 * hsd0_113[k]
                   - f_2 * hsd1_113[k]
                   + f_3 * pc_z[k] * hsf_189[k];

        t_285[k] = pa_y[k] * gsg0_210[k]
                   - f_6 * pc_y[k] * gsg1_210[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, pa_y, pc_x, pc_y, gsg0_212, gsg1_212, hsd0_115, \
                         hsd0_117, hsd1_115, hsd1_117, hsf_191, \
                         hsf_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = f_10 * hsd0_115[k]
                   - f_11 * hsd1_115[k]
                   + f_3 * pc_x[k] * hsf_191[k];

        t_287[k] = pa_y[k] * gsg0_212[k]
                   - f_6 * pc_y[k] * gsg1_212[k];

        t_288[k] = f_4 * hsd0_117[k]
                   - f_5 * hsd1_117[k]
                   + f_3 * pc_x[k] * hsf_193[k];
    }

#pragma omp simd aligned(t_289, t_290, t_291, t_292, t_293, pa_y, pc_x, pc_y, gsg0_215, \
                         gsg1_215, hsd0_118, hsd1_118, hsf_194, hsf_196, hsf_197, \
                         hsf_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_289[k] = f_4 * hsd0_118[k]
                   - f_5 * hsd1_118[k]
                   + f_3 * pc_x[k] * hsf_194[k];

        t_290[k] = pa_y[k] * gsg0_215[k]
                   - f_6 * pc_y[k] * gsg1_215[k];

        t_291[k] = f_3 * pc_x[k] * hsf_196[k];

        t_292[k] = f_3 * pc_x[k] * hsf_197[k];

        t_293[k] = f_3 * pc_x[k] * hsf_198[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, pa_y, pc_x, pc_y, pc_z, gsg0_220, gsf_136, \
                         gsf_146, gsg1_220, hsf_196, hsf_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_3 * pc_x[k] * hsf_199[k];

        t_295[k] = pa_y[k] * gsg0_220[k]
                   + f_9 * gsf_146[k]
                   - f_6 * pc_y[k] * gsg1_220[k];

        t_296[k] = f_9 * gsf_136[k]
                   + f_3 * pc_z[k] * hsf_196[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, pa_y, pc_y, gsg0_222, gsg0_224, gsf_148, \
                         gsf_149, gsg1_222, gsg1_224, hsf_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = pa_y[k] * gsg0_222[k]
                   + f_8 * gsf_148[k]
                   - f_6 * pc_y[k] * gsg1_222[k];

        t_298[k] = f_7 * gsf_149[k]
                   + f_3 * pc_y[k] * hsf_199[k];

        t_299[k] = pa_y[k] * gsg0_224[k]
                   - f_6 * pc_y[k] * gsg1_224[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, pc_x, pc_y, hsd0_120, hsd0_122, \
                         hsd0_123, hsd1_120, hsd1_122, hsd1_123, hsf_200, hsf_202, \
                         hsf_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = f_1 * hsd0_120[k]
                   - f_2 * hsd1_120[k]
                   + f_3 * pc_x[k] * hsf_200[k];

        t_301[k] = f_3 * pc_y[k] * hsf_200[k];

        t_302[k] = f_10 * hsd0_122[k]
                   - f_11 * hsd1_122[k]
                   + f_3 * pc_x[k] * hsf_202[k];

        t_303[k] = f_4 * hsd0_123[k]
                   - f_5 * hsd1_123[k]
                   + f_3 * pc_x[k] * hsf_203[k];

        t_304[k] = f_3 * pc_y[k] * hsf_202[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, pc_x, hsd0_125, hsd1_125, hsf_205, \
                         hsf_206, hsf_207, hsf_208, hsf_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = f_4 * hsd0_125[k]
                   - f_5 * hsd1_125[k]
                   + f_3 * pc_x[k] * hsf_205[k];

        t_306[k] = f_3 * pc_x[k] * hsf_206[k];

        t_307[k] = f_3 * pc_x[k] * hsf_207[k];

        t_308[k] = f_3 * pc_x[k] * hsf_208[k];

        t_309[k] = f_3 * pc_x[k] * hsf_209[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, pc_y, hsd0_123, hsd0_124, hsd0_125, \
                         hsd1_123, hsd1_124, hsd1_125, hsf_206, hsf_207, hsf_208, \
                         hsf_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = f_1 * hsd0_123[k]
                   - f_2 * hsd1_123[k]
                   + f_3 * pc_y[k] * hsf_206[k];

        t_311[k] = f_10 * hsd0_124[k]
                   - f_11 * hsd1_124[k]
                   + f_3 * pc_y[k] * hsf_207[k];

        t_312[k] = f_4 * hsd0_125[k]
                   - f_5 * hsd1_125[k]
                   + f_3 * pc_y[k] * hsf_208[k];

        t_313[k] = f_3 * pc_y[k] * hsf_209[k];
    }

#pragma omp simd aligned(t_314, pc_z, gsf_149, hsd0_125, hsd1_125, \
                         hsf_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_314[k] = f_0 * gsf_149[k]
                   + f_1 * hsd0_125[k]
                   - f_2 * hsd1_125[k]
                   + f_3 * pc_z[k] * hsf_209[k];
    }
}

auto
compute_prim_hsg_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t gsg0, const size_t gsf,
                                                   const size_t gsg1, const size_t hsd0,
                                                   const size_t hsd1, const size_t hsf,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_hsg_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, gsg0, gsf,
                                                              gsg1, hsd0, hsd1, hsf, ncols,
                                                              gamma, p, q);

    compute_prim_hsg_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, gsg0, gsf,
                                                              gsg1, hsd0, hsd1, hsf, ncols,
                                                              gamma, p, q);

    compute_prim_hsg_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, gsg0, gsf,
                                                              gsg1, hsd0, hsd1, hsf, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
