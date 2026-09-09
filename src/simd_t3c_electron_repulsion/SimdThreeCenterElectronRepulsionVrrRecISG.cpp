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


#include "SimdThreeCenterElectronRepulsionVrrRecISG.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_isg_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t hsg0,
                                                          const size_t hsf, const size_t hsg1,
                                                          const size_t isd0, const size_t isd1,
                                                          const size_t isf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 2.5 / q;
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 2.0 / q;
    const auto f_13 = 1.5 / q;

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

    const auto *hsg0_0 = buffer.data(hsg0 + 0);
    const auto *hsg0_3 = buffer.data(hsg0 + 3);
    const auto *hsg0_5 = buffer.data(hsg0 + 5);
    const auto *hsg0_10 = buffer.data(hsg0 + 10);
    const auto *hsg0_14 = buffer.data(hsg0 + 14);
    const auto *hsg0_18 = buffer.data(hsg0 + 18);
    const auto *hsg0_25 = buffer.data(hsg0 + 25);
    const auto *hsg0_30 = buffer.data(hsg0 + 30);
    const auto *hsg0_35 = buffer.data(hsg0 + 35);
    const auto *hsg0_44 = buffer.data(hsg0 + 44);
    const auto *hsg0_45 = buffer.data(hsg0 + 45);
    const auto *hsg0_48 = buffer.data(hsg0 + 48);
    const auto *hsg0_55 = buffer.data(hsg0 + 55);
    const auto *hsg0_75 = buffer.data(hsg0 + 75);
    const auto *hsg0_78 = buffer.data(hsg0 + 78);
    const auto *hsg0_80 = buffer.data(hsg0 + 80);

    const auto *hsf_0 = buffer.data(hsf + 0);
    const auto *hsf_1 = buffer.data(hsf + 1);
    const auto *hsf_2 = buffer.data(hsf + 2);
    const auto *hsf_6 = buffer.data(hsf + 6);
    const auto *hsf_9 = buffer.data(hsf + 9);
    const auto *hsf_10 = buffer.data(hsf + 10);
    const auto *hsf_16 = buffer.data(hsf + 16);
    const auto *hsf_18 = buffer.data(hsf + 18);
    const auto *hsf_19 = buffer.data(hsf + 19);
    const auto *hsf_20 = buffer.data(hsf + 20);
    const auto *hsf_22 = buffer.data(hsf + 22);
    const auto *hsf_26 = buffer.data(hsf + 26);
    const auto *hsf_27 = buffer.data(hsf + 27);
    const auto *hsf_28 = buffer.data(hsf + 28);
    const auto *hsf_29 = buffer.data(hsf + 29);
    const auto *hsf_30 = buffer.data(hsf + 30);
    const auto *hsf_33 = buffer.data(hsf + 33);
    const auto *hsf_36 = buffer.data(hsf + 36);
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
    const auto *hsf_59 = buffer.data(hsf + 59);
    const auto *hsf_60 = buffer.data(hsf + 60);
    const auto *hsf_63 = buffer.data(hsf + 63);
    const auto *hsf_66 = buffer.data(hsf + 66);
    const auto *hsf_68 = buffer.data(hsf + 68);
    const auto *hsf_69 = buffer.data(hsf + 69);
    const auto *hsf_75 = buffer.data(hsf + 75);
    const auto *hsf_76 = buffer.data(hsf + 76);
    const auto *hsf_77 = buffer.data(hsf + 77);
    const auto *hsf_78 = buffer.data(hsf + 78);
    const auto *hsf_79 = buffer.data(hsf + 79);
    const auto *hsf_86 = buffer.data(hsf + 86);
    const auto *hsf_87 = buffer.data(hsf + 87);
    const auto *hsf_88 = buffer.data(hsf + 88);
    const auto *hsf_89 = buffer.data(hsf + 89);

    const auto *hsg1_0 = buffer.data(hsg1 + 0);
    const auto *hsg1_3 = buffer.data(hsg1 + 3);
    const auto *hsg1_5 = buffer.data(hsg1 + 5);
    const auto *hsg1_10 = buffer.data(hsg1 + 10);
    const auto *hsg1_14 = buffer.data(hsg1 + 14);
    const auto *hsg1_18 = buffer.data(hsg1 + 18);
    const auto *hsg1_25 = buffer.data(hsg1 + 25);
    const auto *hsg1_30 = buffer.data(hsg1 + 30);
    const auto *hsg1_35 = buffer.data(hsg1 + 35);
    const auto *hsg1_44 = buffer.data(hsg1 + 44);
    const auto *hsg1_45 = buffer.data(hsg1 + 45);
    const auto *hsg1_48 = buffer.data(hsg1 + 48);
    const auto *hsg1_55 = buffer.data(hsg1 + 55);
    const auto *hsg1_75 = buffer.data(hsg1 + 75);
    const auto *hsg1_78 = buffer.data(hsg1 + 78);
    const auto *hsg1_80 = buffer.data(hsg1 + 80);

    const auto *isd0_0 = buffer.data(isd0 + 0);
    const auto *isd0_3 = buffer.data(isd0 + 3);
    const auto *isd0_5 = buffer.data(isd0 + 5);
    const auto *isd0_9 = buffer.data(isd0 + 9);
    const auto *isd0_16 = buffer.data(isd0 + 16);
    const auto *isd0_17 = buffer.data(isd0 + 17);
    const auto *isd0_18 = buffer.data(isd0 + 18);
    const auto *isd0_21 = buffer.data(isd0 + 21);
    const auto *isd0_23 = buffer.data(isd0 + 23);
    const auto *isd0_29 = buffer.data(isd0 + 29);
    const auto *isd0_30 = buffer.data(isd0 + 30);
    const auto *isd0_33 = buffer.data(isd0 + 33);
    const auto *isd0_34 = buffer.data(isd0 + 34);
    const auto *isd0_35 = buffer.data(isd0 + 35);
    const auto *isd0_36 = buffer.data(isd0 + 36);
    const auto *isd0_39 = buffer.data(isd0 + 39);
    const auto *isd0_41 = buffer.data(isd0 + 41);
    const auto *isd0_47 = buffer.data(isd0 + 47);

    const auto *isd1_0 = buffer.data(isd1 + 0);
    const auto *isd1_3 = buffer.data(isd1 + 3);
    const auto *isd1_5 = buffer.data(isd1 + 5);
    const auto *isd1_9 = buffer.data(isd1 + 9);
    const auto *isd1_16 = buffer.data(isd1 + 16);
    const auto *isd1_17 = buffer.data(isd1 + 17);
    const auto *isd1_18 = buffer.data(isd1 + 18);
    const auto *isd1_21 = buffer.data(isd1 + 21);
    const auto *isd1_23 = buffer.data(isd1 + 23);
    const auto *isd1_29 = buffer.data(isd1 + 29);
    const auto *isd1_30 = buffer.data(isd1 + 30);
    const auto *isd1_33 = buffer.data(isd1 + 33);
    const auto *isd1_34 = buffer.data(isd1 + 34);
    const auto *isd1_35 = buffer.data(isd1 + 35);
    const auto *isd1_36 = buffer.data(isd1 + 36);
    const auto *isd1_39 = buffer.data(isd1 + 39);
    const auto *isd1_41 = buffer.data(isd1 + 41);
    const auto *isd1_47 = buffer.data(isd1 + 47);

    const auto *isf_0 = buffer.data(isf + 0);
    const auto *isf_1 = buffer.data(isf + 1);
    const auto *isf_2 = buffer.data(isf + 2);
    const auto *isf_3 = buffer.data(isf + 3);
    const auto *isf_5 = buffer.data(isf + 5);
    const auto *isf_6 = buffer.data(isf + 6);
    const auto *isf_8 = buffer.data(isf + 8);
    const auto *isf_9 = buffer.data(isf + 9);
    const auto *isf_10 = buffer.data(isf + 10);
    const auto *isf_11 = buffer.data(isf + 11);
    const auto *isf_13 = buffer.data(isf + 13);
    const auto *isf_16 = buffer.data(isf + 16);
    const auto *isf_17 = buffer.data(isf + 17);
    const auto *isf_18 = buffer.data(isf + 18);
    const auto *isf_19 = buffer.data(isf + 19);
    const auto *isf_20 = buffer.data(isf + 20);
    const auto *isf_22 = buffer.data(isf + 22);
    const auto *isf_25 = buffer.data(isf + 25);
    const auto *isf_26 = buffer.data(isf + 26);
    const auto *isf_27 = buffer.data(isf + 27);
    const auto *isf_28 = buffer.data(isf + 28);
    const auto *isf_29 = buffer.data(isf + 29);
    const auto *isf_30 = buffer.data(isf + 30);
    const auto *isf_31 = buffer.data(isf + 31);
    const auto *isf_32 = buffer.data(isf + 32);
    const auto *isf_33 = buffer.data(isf + 33);
    const auto *isf_36 = buffer.data(isf + 36);
    const auto *isf_37 = buffer.data(isf + 37);
    const auto *isf_38 = buffer.data(isf + 38);
    const auto *isf_39 = buffer.data(isf + 39);
    const auto *isf_40 = buffer.data(isf + 40);
    const auto *isf_42 = buffer.data(isf + 42);
    const auto *isf_46 = buffer.data(isf + 46);
    const auto *isf_47 = buffer.data(isf + 47);
    const auto *isf_48 = buffer.data(isf + 48);
    const auto *isf_49 = buffer.data(isf + 49);
    const auto *isf_50 = buffer.data(isf + 50);
    const auto *isf_51 = buffer.data(isf + 51);
    const auto *isf_52 = buffer.data(isf + 52);
    const auto *isf_55 = buffer.data(isf + 55);
    const auto *isf_56 = buffer.data(isf + 56);
    const auto *isf_57 = buffer.data(isf + 57);
    const auto *isf_58 = buffer.data(isf + 58);
    const auto *isf_59 = buffer.data(isf + 59);
    const auto *isf_60 = buffer.data(isf + 60);
    const auto *isf_61 = buffer.data(isf + 61);
    const auto *isf_62 = buffer.data(isf + 62);
    const auto *isf_63 = buffer.data(isf + 63);
    const auto *isf_66 = buffer.data(isf + 66);
    const auto *isf_67 = buffer.data(isf + 67);
    const auto *isf_68 = buffer.data(isf + 68);
    const auto *isf_69 = buffer.data(isf + 69);
    const auto *isf_70 = buffer.data(isf + 70);
    const auto *isf_72 = buffer.data(isf + 72);
    const auto *isf_75 = buffer.data(isf + 75);
    const auto *isf_76 = buffer.data(isf + 76);
    const auto *isf_77 = buffer.data(isf + 77);
    const auto *isf_78 = buffer.data(isf + 78);
    const auto *isf_79 = buffer.data(isf + 79);
    const auto *isf_80 = buffer.data(isf + 80);
    const auto *isf_82 = buffer.data(isf + 82);
    const auto *isf_86 = buffer.data(isf + 86);
    const auto *isf_87 = buffer.data(isf + 87);
    const auto *isf_88 = buffer.data(isf + 88);
    const auto *isf_89 = buffer.data(isf + 89);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, hsf_0, isd0_0, \
                         isd1_0, isf_0, isf_1, isf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hsf_0[k]
                 + f_1 * isd0_0[k]
                 - f_2 * isd1_0[k]
                 + f_3 * pc_x[k] * isf_0[k];

        t_1[k] = f_3 * pc_y[k] * isf_0[k];

        t_2[k] = f_3 * pc_z[k] * isf_0[k];

        t_3[k] = f_4 * isd0_0[k]
                 - f_5 * isd1_0[k]
                 + f_3 * pc_y[k] * isf_1[k];

        t_4[k] = f_3 * pc_y[k] * isf_2[k];

        t_5[k] = f_4 * isd0_0[k]
                 - f_5 * isd1_0[k]
                 + f_3 * pc_z[k] * isf_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_x, pc_y, pc_z, hsf_6, hsf_9, isd0_3, \
                         isd1_3, isf_3, isf_5, isf_6, isf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * hsf_6[k]
                 + f_3 * pc_x[k] * isf_6[k];

        t_7[k] = f_3 * pc_z[k] * isf_3[k];

        t_8[k] = f_3 * pc_y[k] * isf_5[k];

        t_9[k] = f_0 * hsf_9[k]
                 + f_3 * pc_x[k] * isf_9[k];

        t_10[k] = f_1 * isd0_3[k]
                  - f_2 * isd1_3[k]
                  + f_3 * pc_y[k] * isf_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pa_y, pc_y, pc_z, hsg0_0, hsg1_0, \
                         isd0_5, isd1_5, isf_6, isf_8, isf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * pc_z[k] * isf_6[k];

        t_12[k] = f_4 * isd0_5[k]
                  - f_5 * isd1_5[k]
                  + f_3 * pc_y[k] * isf_8[k];

        t_13[k] = f_3 * pc_y[k] * isf_9[k];

        t_14[k] = f_1 * isd0_5[k]
                  - f_2 * isd1_5[k]
                  + f_3 * pc_z[k] * isf_9[k];

        t_15[k] = pa_y[k] * hsg0_0[k]
                  - f_6 * pc_y[k] * hsg1_0[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pa_y, pc_y, pc_z, hsg0_3, hsg0_5, \
                         hsf_0, hsf_1, hsg1_3, hsg1_5, isf_10, isf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_7 * hsf_0[k]
                  + f_3 * pc_y[k] * isf_10[k];

        t_17[k] = f_3 * pc_z[k] * isf_10[k];

        t_18[k] = pa_y[k] * hsg0_3[k]
                  + f_8 * hsf_1[k]
                  - f_6 * pc_y[k] * hsg1_3[k];

        t_19[k] = f_3 * pc_z[k] * isf_11[k];

        t_20[k] = pa_y[k] * hsg0_5[k]
                  - f_6 * pc_y[k] * hsg1_5[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pc_x, pc_z, hsf_16, hsf_18, hsf_19, isf_13, \
                         isf_16, isf_18, isf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_9 * hsf_16[k]
                  + f_3 * pc_x[k] * isf_16[k];

        t_22[k] = f_3 * pc_z[k] * isf_13[k];

        t_23[k] = f_9 * hsf_18[k]
                  + f_3 * pc_x[k] * isf_18[k];

        t_24[k] = f_9 * hsf_19[k]
                  + f_3 * pc_x[k] * isf_19[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pc_y, pc_z, hsf_6, hsf_9, isd0_9, isd1_9, \
                         isf_16, isf_17, isf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_7 * hsf_6[k]
                  + f_1 * isd0_9[k]
                  - f_2 * isd1_9[k]
                  + f_3 * pc_y[k] * isf_16[k];

        t_26[k] = f_3 * pc_z[k] * isf_16[k];

        t_27[k] = f_4 * isd0_9[k]
                  - f_5 * isd1_9[k]
                  + f_3 * pc_z[k] * isf_17[k];

        t_28[k] = f_7 * hsf_9[k]
                  + f_3 * pc_y[k] * isf_19[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pa_y, pa_z, pc_y, pc_z, hsg0_0, hsg0_14, \
                         hsf_0, hsg1_0, hsg1_14, isf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pa_y[k] * hsg0_14[k]
                  - f_6 * pc_y[k] * hsg1_14[k];

        t_30[k] = pa_z[k] * hsg0_0[k]
                  - f_6 * pc_z[k] * hsg1_0[k];

        t_31[k] = f_3 * pc_y[k] * isf_20[k];

        t_32[k] = f_7 * hsf_0[k]
                  + f_3 * pc_z[k] * isf_20[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pa_z, pc_x, pc_y, pc_z, hsg0_3, hsg0_5, \
                         hsf_2, hsf_26, hsg1_3, hsg1_5, isf_22, \
                         isf_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pa_z[k] * hsg0_3[k]
                  - f_6 * pc_z[k] * hsg1_3[k];

        t_34[k] = f_3 * pc_y[k] * isf_22[k];

        t_35[k] = pa_z[k] * hsg0_5[k]
                  + f_8 * hsf_2[k]
                  - f_6 * pc_z[k] * hsg1_5[k];

        t_36[k] = f_9 * hsf_26[k]
                  + f_3 * pc_x[k] * isf_26[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pa_z, pc_x, pc_y, pc_z, hsg0_10, hsf_27, \
                         hsf_29, hsg1_10, isf_25, isf_27, isf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_9 * hsf_27[k]
                  + f_3 * pc_x[k] * isf_27[k];

        t_38[k] = f_3 * pc_y[k] * isf_25[k];

        t_39[k] = f_9 * hsf_29[k]
                  + f_3 * pc_x[k] * isf_29[k];

        t_40[k] = pa_z[k] * hsg0_10[k]
                  - f_6 * pc_z[k] * hsg1_10[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pc_y, pc_z, hsf_9, isd0_16, isd0_17, isd1_16, \
                         isd1_17, isf_27, isf_28, isf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_10 * isd0_16[k]
                  - f_11 * isd1_16[k]
                  + f_3 * pc_y[k] * isf_27[k];

        t_42[k] = f_4 * isd0_17[k]
                  - f_5 * isd1_17[k]
                  + f_3 * pc_y[k] * isf_28[k];

        t_43[k] = f_3 * pc_y[k] * isf_29[k];

        t_44[k] = f_7 * hsf_9[k]
                  + f_1 * isd0_17[k]
                  - f_2 * isd1_17[k]
                  + f_3 * pc_z[k] * isf_29[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pc_x, pc_y, pc_z, hsf_10, hsf_30, hsf_33, \
                         isd0_18, isd0_21, isd1_18, isd1_21, isf_30, \
                         isf_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_12 * hsf_30[k]
                  + f_1 * isd0_18[k]
                  - f_2 * isd1_18[k]
                  + f_3 * pc_x[k] * isf_30[k];

        t_46[k] = f_8 * hsf_10[k]
                  + f_3 * pc_y[k] * isf_30[k];

        t_47[k] = f_3 * pc_z[k] * isf_30[k];

        t_48[k] = f_12 * hsf_33[k]
                  + f_4 * isd0_21[k]
                  - f_5 * isd1_21[k]
                  + f_3 * pc_x[k] * isf_33[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, pc_x, pc_z, hsf_36, hsf_38, isd0_18, \
                         isd1_18, isf_31, isf_32, isf_33, isf_36, \
                         isf_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_3 * pc_z[k] * isf_31[k];

        t_50[k] = f_4 * isd0_18[k]
                  - f_5 * isd1_18[k]
                  + f_3 * pc_z[k] * isf_32[k];

        t_51[k] = f_12 * hsf_36[k]
                  + f_3 * pc_x[k] * isf_36[k];

        t_52[k] = f_3 * pc_z[k] * isf_33[k];

        t_53[k] = f_12 * hsf_38[k]
                  + f_3 * pc_x[k] * isf_38[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, pc_x, pc_y, pc_z, hsf_16, hsf_19, \
                         hsf_39, isd0_21, isd1_21, isf_36, isf_37, \
                         isf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_12 * hsf_39[k]
                  + f_3 * pc_x[k] * isf_39[k];

        t_55[k] = f_8 * hsf_16[k]
                  + f_1 * isd0_21[k]
                  - f_2 * isd1_21[k]
                  + f_3 * pc_y[k] * isf_36[k];

        t_56[k] = f_3 * pc_z[k] * isf_36[k];

        t_57[k] = f_4 * isd0_21[k]
                  - f_5 * isd1_21[k]
                  + f_3 * pc_z[k] * isf_37[k];

        t_58[k] = f_8 * hsf_19[k]
                  + f_3 * pc_y[k] * isf_39[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pa_y, pc_y, pc_z, hsg0_30, hsf_10, hsf_20, \
                         hsg1_30, isd0_23, isd1_23, isf_39, isf_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_1 * isd0_23[k]
                  - f_2 * isd1_23[k]
                  + f_3 * pc_z[k] * isf_39[k];

        t_60[k] = pa_y[k] * hsg0_30[k]
                  - f_6 * pc_y[k] * hsg1_30[k];

        t_61[k] = f_7 * hsf_20[k]
                  + f_3 * pc_y[k] * isf_40[k];

        t_62[k] = f_7 * hsf_10[k]
                  + f_3 * pc_z[k] * isf_40[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, pa_y, pa_z, pc_y, pc_z, hsg0_18, hsg0_35, hsf_22, \
                         hsg1_18, hsg1_35, isf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = pa_z[k] * hsg0_18[k]
                  - f_6 * pc_z[k] * hsg1_18[k];

        t_64[k] = f_7 * hsf_22[k]
                  + f_3 * pc_y[k] * isf_42[k];

        t_65[k] = pa_y[k] * hsg0_35[k]
                  - f_6 * pc_y[k] * hsg1_35[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pc_x, hsf_46, hsf_47, hsf_48, hsf_49, isf_46, \
                         isf_47, isf_48, isf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_12 * hsf_46[k]
                  + f_3 * pc_x[k] * isf_46[k];

        t_67[k] = f_12 * hsf_47[k]
                  + f_3 * pc_x[k] * isf_47[k];

        t_68[k] = f_12 * hsf_48[k]
                  + f_3 * pc_x[k] * isf_48[k];

        t_69[k] = f_12 * hsf_49[k]
                  + f_3 * pc_x[k] * isf_49[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, pa_z, pc_y, pc_z, hsg0_25, hsf_16, hsf_28, hsg1_25, \
                         isd0_29, isd1_29, isf_46, isf_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = pa_z[k] * hsg0_25[k]
                  - f_6 * pc_z[k] * hsg1_25[k];

        t_71[k] = f_7 * hsf_16[k]
                  + f_3 * pc_z[k] * isf_46[k];

        t_72[k] = f_7 * hsf_28[k]
                  + f_4 * isd0_29[k]
                  - f_5 * isd1_29[k]
                  + f_3 * pc_y[k] * isf_48[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pa_y, pc_x, pc_y, hsg0_44, hsf_29, hsf_50, \
                         hsg1_44, isd0_30, isd1_30, isf_49, isf_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_7 * hsf_29[k]
                  + f_3 * pc_y[k] * isf_49[k];

        t_74[k] = pa_y[k] * hsg0_44[k]
                  - f_6 * pc_y[k] * hsg1_44[k];

        t_75[k] = f_12 * hsf_50[k]
                  + f_1 * isd0_30[k]
                  - f_2 * isd1_30[k]
                  + f_3 * pc_x[k] * isf_50[k];

        t_76[k] = f_3 * pc_y[k] * isf_50[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, pc_y, pc_z, hsf_20, isd0_30, isd1_30, isf_50, \
                         isf_51, isf_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_8 * hsf_20[k]
                  + f_3 * pc_z[k] * isf_50[k];

        t_78[k] = f_4 * isd0_30[k]
                  - f_5 * isd1_30[k]
                  + f_3 * pc_y[k] * isf_51[k];

        t_79[k] = f_3 * pc_y[k] * isf_52[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pc_x, pc_y, hsf_55, hsf_56, hsf_57, isd0_35, \
                         isd1_35, isf_55, isf_56, isf_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_12 * hsf_55[k]
                  + f_4 * isd0_35[k]
                  - f_5 * isd1_35[k]
                  + f_3 * pc_x[k] * isf_55[k];

        t_81[k] = f_12 * hsf_56[k]
                  + f_3 * pc_x[k] * isf_56[k];

        t_82[k] = f_12 * hsf_57[k]
                  + f_3 * pc_x[k] * isf_57[k];

        t_83[k] = f_3 * pc_y[k] * isf_55[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, pc_x, pc_y, hsf_59, isd0_33, isd0_34, isd1_33, \
                         isd1_34, isf_56, isf_57, isf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_12 * hsf_59[k]
                  + f_3 * pc_x[k] * isf_59[k];

        t_85[k] = f_1 * isd0_33[k]
                  - f_2 * isd1_33[k]
                  + f_3 * pc_y[k] * isf_56[k];

        t_86[k] = f_10 * isd0_34[k]
                  - f_11 * isd1_34[k]
                  + f_3 * pc_y[k] * isf_57[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pc_x, pc_y, pc_z, hsf_29, hsf_60, isd0_35, \
                         isd0_36, isd1_35, isd1_36, isf_58, isf_59, \
                         isf_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_4 * isd0_35[k]
                  - f_5 * isd1_35[k]
                  + f_3 * pc_y[k] * isf_58[k];

        t_88[k] = f_3 * pc_y[k] * isf_59[k];

        t_89[k] = f_8 * hsf_29[k]
                  + f_1 * isd0_35[k]
                  - f_2 * isd1_35[k]
                  + f_3 * pc_z[k] * isf_59[k];

        t_90[k] = f_13 * hsf_60[k]
                  + f_1 * isd0_36[k]
                  - f_2 * isd1_36[k]
                  + f_3 * pc_x[k] * isf_60[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, pc_x, pc_y, pc_z, hsf_30, hsf_63, isd0_39, \
                         isd1_39, isf_60, isf_61, isf_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_13 * hsf_30[k]
                  + f_3 * pc_y[k] * isf_60[k];

        t_92[k] = f_3 * pc_z[k] * isf_60[k];

        t_93[k] = f_13 * hsf_63[k]
                  + f_4 * isd0_39[k]
                  - f_5 * isd1_39[k]
                  + f_3 * pc_x[k] * isf_63[k];

        t_94[k] = f_3 * pc_z[k] * isf_61[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pc_x, pc_z, hsf_66, hsf_68, isd0_36, isd1_36, \
                         isf_62, isf_63, isf_66, isf_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_4 * isd0_36[k]
                  - f_5 * isd1_36[k]
                  + f_3 * pc_z[k] * isf_62[k];

        t_96[k] = f_13 * hsf_66[k]
                  + f_3 * pc_x[k] * isf_66[k];

        t_97[k] = f_3 * pc_z[k] * isf_63[k];

        t_98[k] = f_13 * hsf_68[k]
                  + f_3 * pc_x[k] * isf_68[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, pc_x, pc_y, pc_z, hsf_36, hsf_39, \
                         hsf_69, isd0_39, isd1_39, isf_66, isf_67, \
                         isf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_13 * hsf_69[k]
                  + f_3 * pc_x[k] * isf_69[k];

        t_100[k] = f_13 * hsf_36[k]
                   + f_1 * isd0_39[k]
                   - f_2 * isd1_39[k]
                   + f_3 * pc_y[k] * isf_66[k];

        t_101[k] = f_3 * pc_z[k] * isf_66[k];

        t_102[k] = f_4 * isd0_39[k]
                   - f_5 * isd1_39[k]
                   + f_3 * pc_z[k] * isf_67[k];

        t_103[k] = f_13 * hsf_39[k]
                   + f_3 * pc_y[k] * isf_69[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pa_z, pc_y, pc_z, hsg0_45, hsf_30, \
                         hsf_40, hsg1_45, isd0_41, isd1_41, isf_69, \
                         isf_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_1 * isd0_41[k]
                   - f_2 * isd1_41[k]
                   + f_3 * pc_z[k] * isf_69[k];

        t_105[k] = pa_z[k] * hsg0_45[k]
                   - f_6 * pc_z[k] * hsg1_45[k];

        t_106[k] = f_8 * hsf_40[k]
                   + f_3 * pc_y[k] * isf_70[k];

        t_107[k] = f_7 * hsf_30[k]
                   + f_3 * pc_z[k] * isf_70[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, pa_z, pc_x, pc_y, pc_z, hsg0_48, hsf_42, hsf_75, \
                         hsg1_48, isd0_47, isd1_47, isf_72, isf_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = pa_z[k] * hsg0_48[k]
                   - f_6 * pc_z[k] * hsg1_48[k];

        t_109[k] = f_8 * hsf_42[k]
                   + f_3 * pc_y[k] * isf_72[k];

        t_110[k] = f_13 * hsf_75[k]
                   + f_4 * isd0_47[k]
                   - f_5 * isd1_47[k]
                   + f_3 * pc_x[k] * isf_75[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, pc_x, hsf_76, hsf_77, hsf_78, hsf_79, \
                         isf_76, isf_77, isf_78, isf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_13 * hsf_76[k]
                   + f_3 * pc_x[k] * isf_76[k];

        t_112[k] = f_13 * hsf_77[k]
                   + f_3 * pc_x[k] * isf_77[k];

        t_113[k] = f_13 * hsf_78[k]
                   + f_3 * pc_x[k] * isf_78[k];

        t_114[k] = f_13 * hsf_79[k]
                   + f_3 * pc_x[k] * isf_79[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, pa_z, pc_y, pc_z, hsg0_55, hsf_36, hsf_48, \
                         hsg1_55, isd0_47, isd1_47, isf_76, isf_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = pa_z[k] * hsg0_55[k]
                   - f_6 * pc_z[k] * hsg1_55[k];

        t_116[k] = f_7 * hsf_36[k]
                   + f_3 * pc_z[k] * isf_76[k];

        t_117[k] = f_8 * hsf_48[k]
                   + f_4 * isd0_47[k]
                   - f_5 * isd1_47[k]
                   + f_3 * pc_y[k] * isf_78[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, pa_y, pc_y, pc_z, hsg0_75, hsf_39, \
                         hsf_49, hsf_50, hsg1_75, isd0_47, isd1_47, isf_79, \
                         isf_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_8 * hsf_49[k]
                   + f_3 * pc_y[k] * isf_79[k];

        t_119[k] = f_7 * hsf_39[k]
                   + f_1 * isd0_47[k]
                   - f_2 * isd1_47[k]
                   + f_3 * pc_z[k] * isf_79[k];

        t_120[k] = pa_y[k] * hsg0_75[k]
                   - f_6 * pc_y[k] * hsg1_75[k];

        t_121[k] = f_7 * hsf_50[k]
                   + f_3 * pc_y[k] * isf_80[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, pa_y, pc_y, pc_z, hsg0_78, hsg0_80, \
                         hsf_40, hsf_51, hsf_52, hsg1_78, hsg1_80, isf_80, \
                         isf_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_8 * hsf_40[k]
                   + f_3 * pc_z[k] * isf_80[k];

        t_123[k] = pa_y[k] * hsg0_78[k]
                   + f_8 * hsf_51[k]
                   - f_6 * pc_y[k] * hsg1_78[k];

        t_124[k] = f_7 * hsf_52[k]
                   + f_3 * pc_y[k] * isf_82[k];

        t_125[k] = pa_y[k] * hsg0_80[k]
                   - f_6 * pc_y[k] * hsg1_80[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, pc_x, hsf_86, hsf_87, hsf_88, hsf_89, \
                         isf_86, isf_87, isf_88, isf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_13 * hsf_86[k]
                   + f_3 * pc_x[k] * isf_86[k];

        t_127[k] = f_13 * hsf_87[k]
                   + f_3 * pc_x[k] * isf_87[k];

        t_128[k] = f_13 * hsf_88[k]
                   + f_3 * pc_x[k] * isf_88[k];

        t_129[k] = f_13 * hsf_89[k]
                   + f_3 * pc_x[k] * isf_89[k];
    }
}

static auto
compute_prim_isg_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t hsg0,
                                                          const size_t hsf, const size_t hsg1,
                                                          const size_t isd0, const size_t isd1,
                                                          const size_t isf, const size_t ncols,
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
    const auto f_9 = 2.5 / q;
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 2.0 / q;
    const auto f_13 = 1.5 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *hsg0_89 = buffer.data(hsg0 + 89);
    const auto *hsg0_90 = buffer.data(hsg0 + 90);
    const auto *hsg0_93 = buffer.data(hsg0 + 93);
    const auto *hsg0_100 = buffer.data(hsg0 + 100);
    const auto *hsg0_135 = buffer.data(hsg0 + 135);
    const auto *hsg0_138 = buffer.data(hsg0 + 138);
    const auto *hsg0_140 = buffer.data(hsg0 + 140);
    const auto *hsg0_149 = buffer.data(hsg0 + 149);
    const auto *hsg0_150 = buffer.data(hsg0 + 150);
    const auto *hsg0_153 = buffer.data(hsg0 + 153);
    const auto *hsg0_225 = buffer.data(hsg0 + 225);
    const auto *hsg0_228 = buffer.data(hsg0 + 228);
    const auto *hsg0_235 = buffer.data(hsg0 + 235);
    const auto *hsg0_237 = buffer.data(hsg0 + 237);
    const auto *hsg0_239 = buffer.data(hsg0 + 239);
    const auto *hsg0_245 = buffer.data(hsg0 + 245);
    const auto *hsg0_250 = buffer.data(hsg0 + 250);
    const auto *hsg0_252 = buffer.data(hsg0 + 252);

    const auto *hsf_46 = buffer.data(hsf + 46);
    const auto *hsf_50 = buffer.data(hsf + 50);
    const auto *hsf_56 = buffer.data(hsf + 56);
    const auto *hsf_58 = buffer.data(hsf + 58);
    const auto *hsf_59 = buffer.data(hsf + 59);
    const auto *hsf_60 = buffer.data(hsf + 60);
    const auto *hsf_66 = buffer.data(hsf + 66);
    const auto *hsf_69 = buffer.data(hsf + 69);
    const auto *hsf_70 = buffer.data(hsf + 70);
    const auto *hsf_72 = buffer.data(hsf + 72);
    const auto *hsf_76 = buffer.data(hsf + 76);
    const auto *hsf_78 = buffer.data(hsf + 78);
    const auto *hsf_79 = buffer.data(hsf + 79);
    const auto *hsf_80 = buffer.data(hsf + 80);
    const auto *hsf_82 = buffer.data(hsf + 82);
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
    const auto *hsf_103 = buffer.data(hsf + 103);
    const auto *hsf_106 = buffer.data(hsf + 106);
    const auto *hsf_108 = buffer.data(hsf + 108);
    const auto *hsf_109 = buffer.data(hsf + 109);
    const auto *hsf_110 = buffer.data(hsf + 110);
    const auto *hsf_112 = buffer.data(hsf + 112);
    const auto *hsf_115 = buffer.data(hsf + 115);
    const auto *hsf_116 = buffer.data(hsf + 116);
    const auto *hsf_117 = buffer.data(hsf + 117);
    const auto *hsf_118 = buffer.data(hsf + 118);
    const auto *hsf_119 = buffer.data(hsf + 119);
    const auto *hsf_120 = buffer.data(hsf + 120);
    const auto *hsf_123 = buffer.data(hsf + 123);
    const auto *hsf_125 = buffer.data(hsf + 125);
    const auto *hsf_126 = buffer.data(hsf + 126);
    const auto *hsf_127 = buffer.data(hsf + 127);
    const auto *hsf_128 = buffer.data(hsf + 128);
    const auto *hsf_129 = buffer.data(hsf + 129);
    const auto *hsf_136 = buffer.data(hsf + 136);
    const auto *hsf_137 = buffer.data(hsf + 137);
    const auto *hsf_138 = buffer.data(hsf + 138);
    const auto *hsf_139 = buffer.data(hsf + 139);
    const auto *hsf_140 = buffer.data(hsf + 140);
    const auto *hsf_145 = buffer.data(hsf + 145);
    const auto *hsf_146 = buffer.data(hsf + 146);
    const auto *hsf_147 = buffer.data(hsf + 147);
    const auto *hsf_149 = buffer.data(hsf + 149);
    const auto *hsf_150 = buffer.data(hsf + 150);
    const auto *hsf_153 = buffer.data(hsf + 153);
    const auto *hsf_156 = buffer.data(hsf + 156);
    const auto *hsf_158 = buffer.data(hsf + 158);
    const auto *hsf_159 = buffer.data(hsf + 159);
    const auto *hsf_165 = buffer.data(hsf + 165);
    const auto *hsf_166 = buffer.data(hsf + 166);
    const auto *hsf_167 = buffer.data(hsf + 167);
    const auto *hsf_168 = buffer.data(hsf + 168);
    const auto *hsf_169 = buffer.data(hsf + 169);

    const auto *hsg1_89 = buffer.data(hsg1 + 89);
    const auto *hsg1_90 = buffer.data(hsg1 + 90);
    const auto *hsg1_93 = buffer.data(hsg1 + 93);
    const auto *hsg1_100 = buffer.data(hsg1 + 100);
    const auto *hsg1_135 = buffer.data(hsg1 + 135);
    const auto *hsg1_138 = buffer.data(hsg1 + 138);
    const auto *hsg1_140 = buffer.data(hsg1 + 140);
    const auto *hsg1_149 = buffer.data(hsg1 + 149);
    const auto *hsg1_150 = buffer.data(hsg1 + 150);
    const auto *hsg1_153 = buffer.data(hsg1 + 153);
    const auto *hsg1_225 = buffer.data(hsg1 + 225);
    const auto *hsg1_228 = buffer.data(hsg1 + 228);
    const auto *hsg1_235 = buffer.data(hsg1 + 235);
    const auto *hsg1_237 = buffer.data(hsg1 + 237);
    const auto *hsg1_239 = buffer.data(hsg1 + 239);
    const auto *hsg1_245 = buffer.data(hsg1 + 245);
    const auto *hsg1_250 = buffer.data(hsg1 + 250);
    const auto *hsg1_252 = buffer.data(hsg1 + 252);

    const auto *isd0_51 = buffer.data(isd0 + 51);
    const auto *isd0_53 = buffer.data(isd0 + 53);
    const auto *isd0_54 = buffer.data(isd0 + 54);
    const auto *isd0_57 = buffer.data(isd0 + 57);
    const auto *isd0_58 = buffer.data(isd0 + 58);
    const auto *isd0_59 = buffer.data(isd0 + 59);
    const auto *isd0_60 = buffer.data(isd0 + 60);
    const auto *isd0_63 = buffer.data(isd0 + 63);
    const auto *isd0_65 = buffer.data(isd0 + 65);
    const auto *isd0_71 = buffer.data(isd0 + 71);
    const auto *isd0_72 = buffer.data(isd0 + 72);
    const auto *isd0_75 = buffer.data(isd0 + 75);
    const auto *isd0_77 = buffer.data(isd0 + 77);
    const auto *isd0_81 = buffer.data(isd0 + 81);
    const auto *isd0_83 = buffer.data(isd0 + 83);
    const auto *isd0_84 = buffer.data(isd0 + 84);
    const auto *isd0_87 = buffer.data(isd0 + 87);
    const auto *isd0_88 = buffer.data(isd0 + 88);
    const auto *isd0_89 = buffer.data(isd0 + 89);
    const auto *isd0_90 = buffer.data(isd0 + 90);

    const auto *isd1_51 = buffer.data(isd1 + 51);
    const auto *isd1_53 = buffer.data(isd1 + 53);
    const auto *isd1_54 = buffer.data(isd1 + 54);
    const auto *isd1_57 = buffer.data(isd1 + 57);
    const auto *isd1_58 = buffer.data(isd1 + 58);
    const auto *isd1_59 = buffer.data(isd1 + 59);
    const auto *isd1_60 = buffer.data(isd1 + 60);
    const auto *isd1_63 = buffer.data(isd1 + 63);
    const auto *isd1_65 = buffer.data(isd1 + 65);
    const auto *isd1_71 = buffer.data(isd1 + 71);
    const auto *isd1_72 = buffer.data(isd1 + 72);
    const auto *isd1_75 = buffer.data(isd1 + 75);
    const auto *isd1_77 = buffer.data(isd1 + 77);
    const auto *isd1_81 = buffer.data(isd1 + 81);
    const auto *isd1_83 = buffer.data(isd1 + 83);
    const auto *isd1_84 = buffer.data(isd1 + 84);
    const auto *isd1_87 = buffer.data(isd1 + 87);
    const auto *isd1_88 = buffer.data(isd1 + 88);
    const auto *isd1_89 = buffer.data(isd1 + 89);
    const auto *isd1_90 = buffer.data(isd1 + 90);

    const auto *isf_86 = buffer.data(isf + 86);
    const auto *isf_88 = buffer.data(isf + 88);
    const auto *isf_89 = buffer.data(isf + 89);
    const auto *isf_90 = buffer.data(isf + 90);
    const auto *isf_91 = buffer.data(isf + 91);
    const auto *isf_92 = buffer.data(isf + 92);
    const auto *isf_95 = buffer.data(isf + 95);
    const auto *isf_96 = buffer.data(isf + 96);
    const auto *isf_97 = buffer.data(isf + 97);
    const auto *isf_98 = buffer.data(isf + 98);
    const auto *isf_99 = buffer.data(isf + 99);
    const auto *isf_100 = buffer.data(isf + 100);
    const auto *isf_101 = buffer.data(isf + 101);
    const auto *isf_102 = buffer.data(isf + 102);
    const auto *isf_103 = buffer.data(isf + 103);
    const auto *isf_106 = buffer.data(isf + 106);
    const auto *isf_107 = buffer.data(isf + 107);
    const auto *isf_108 = buffer.data(isf + 108);
    const auto *isf_109 = buffer.data(isf + 109);
    const auto *isf_110 = buffer.data(isf + 110);
    const auto *isf_112 = buffer.data(isf + 112);
    const auto *isf_115 = buffer.data(isf + 115);
    const auto *isf_116 = buffer.data(isf + 116);
    const auto *isf_117 = buffer.data(isf + 117);
    const auto *isf_118 = buffer.data(isf + 118);
    const auto *isf_119 = buffer.data(isf + 119);
    const auto *isf_120 = buffer.data(isf + 120);
    const auto *isf_122 = buffer.data(isf + 122);
    const auto *isf_123 = buffer.data(isf + 123);
    const auto *isf_125 = buffer.data(isf + 125);
    const auto *isf_126 = buffer.data(isf + 126);
    const auto *isf_127 = buffer.data(isf + 127);
    const auto *isf_128 = buffer.data(isf + 128);
    const auto *isf_129 = buffer.data(isf + 129);
    const auto *isf_130 = buffer.data(isf + 130);
    const auto *isf_132 = buffer.data(isf + 132);
    const auto *isf_136 = buffer.data(isf + 136);
    const auto *isf_137 = buffer.data(isf + 137);
    const auto *isf_138 = buffer.data(isf + 138);
    const auto *isf_139 = buffer.data(isf + 139);
    const auto *isf_140 = buffer.data(isf + 140);
    const auto *isf_141 = buffer.data(isf + 141);
    const auto *isf_142 = buffer.data(isf + 142);
    const auto *isf_145 = buffer.data(isf + 145);
    const auto *isf_146 = buffer.data(isf + 146);
    const auto *isf_147 = buffer.data(isf + 147);
    const auto *isf_148 = buffer.data(isf + 148);
    const auto *isf_149 = buffer.data(isf + 149);
    const auto *isf_150 = buffer.data(isf + 150);
    const auto *isf_151 = buffer.data(isf + 151);
    const auto *isf_152 = buffer.data(isf + 152);
    const auto *isf_153 = buffer.data(isf + 153);
    const auto *isf_156 = buffer.data(isf + 156);
    const auto *isf_158 = buffer.data(isf + 158);
    const auto *isf_159 = buffer.data(isf + 159);
    const auto *isf_160 = buffer.data(isf + 160);
    const auto *isf_162 = buffer.data(isf + 162);
    const auto *isf_166 = buffer.data(isf + 166);
    const auto *isf_167 = buffer.data(isf + 167);
    const auto *isf_168 = buffer.data(isf + 168);
    const auto *isf_169 = buffer.data(isf + 169);

#pragma omp simd aligned(t_130, t_131, t_132, pc_y, pc_z, hsf_46, hsf_56, hsf_58, isd0_51, \
                         isd0_53, isd1_51, isd1_53, isf_86, isf_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_7 * hsf_56[k]
                   + f_1 * isd0_51[k]
                   - f_2 * isd1_51[k]
                   + f_3 * pc_y[k] * isf_86[k];

        t_131[k] = f_8 * hsf_46[k]
                   + f_3 * pc_z[k] * isf_86[k];

        t_132[k] = f_7 * hsf_58[k]
                   + f_4 * isd0_53[k]
                   - f_5 * isd1_53[k]
                   + f_3 * pc_y[k] * isf_88[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, pa_y, pc_x, pc_y, hsg0_89, hsf_59, \
                         hsf_90, hsg1_89, isd0_54, isd1_54, isf_89, \
                         isf_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_7 * hsf_59[k]
                   + f_3 * pc_y[k] * isf_89[k];

        t_134[k] = pa_y[k] * hsg0_89[k]
                   - f_6 * pc_y[k] * hsg1_89[k];

        t_135[k] = f_13 * hsf_90[k]
                   + f_1 * isd0_54[k]
                   - f_2 * isd1_54[k]
                   + f_3 * pc_x[k] * isf_90[k];

        t_136[k] = f_3 * pc_y[k] * isf_90[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, pc_y, pc_z, hsf_50, isd0_54, isd1_54, isf_90, \
                         isf_91, isf_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_13 * hsf_50[k]
                   + f_3 * pc_z[k] * isf_90[k];

        t_138[k] = f_4 * isd0_54[k]
                   - f_5 * isd1_54[k]
                   + f_3 * pc_y[k] * isf_91[k];

        t_139[k] = f_3 * pc_y[k] * isf_92[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, pc_x, pc_y, hsf_95, hsf_96, hsf_97, \
                         isd0_59, isd1_59, isf_95, isf_96, isf_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_13 * hsf_95[k]
                   + f_4 * isd0_59[k]
                   - f_5 * isd1_59[k]
                   + f_3 * pc_x[k] * isf_95[k];

        t_141[k] = f_13 * hsf_96[k]
                   + f_3 * pc_x[k] * isf_96[k];

        t_142[k] = f_13 * hsf_97[k]
                   + f_3 * pc_x[k] * isf_97[k];

        t_143[k] = f_3 * pc_y[k] * isf_95[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, pc_x, pc_y, hsf_99, isd0_57, isd0_58, isd1_57, \
                         isd1_58, isf_96, isf_97, isf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_13 * hsf_99[k]
                   + f_3 * pc_x[k] * isf_99[k];

        t_145[k] = f_1 * isd0_57[k]
                   - f_2 * isd1_57[k]
                   + f_3 * pc_y[k] * isf_96[k];

        t_146[k] = f_10 * isd0_58[k]
                   - f_11 * isd1_58[k]
                   + f_3 * pc_y[k] * isf_97[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, pc_x, pc_y, pc_z, hsf_59, hsf_100, \
                         isd0_59, isd0_60, isd1_59, isd1_60, isf_98, isf_99, \
                         isf_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = f_4 * isd0_59[k]
                   - f_5 * isd1_59[k]
                   + f_3 * pc_y[k] * isf_98[k];

        t_148[k] = f_3 * pc_y[k] * isf_99[k];

        t_149[k] = f_13 * hsf_59[k]
                   + f_1 * isd0_59[k]
                   - f_2 * isd1_59[k]
                   + f_3 * pc_z[k] * isf_99[k];

        t_150[k] = f_8 * hsf_100[k]
                   + f_1 * isd0_60[k]
                   - f_2 * isd1_60[k]
                   + f_3 * pc_x[k] * isf_100[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, pc_x, pc_y, pc_z, hsf_60, hsf_103, \
                         isd0_63, isd1_63, isf_100, isf_101, isf_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_12 * hsf_60[k]
                   + f_3 * pc_y[k] * isf_100[k];

        t_152[k] = f_3 * pc_z[k] * isf_100[k];

        t_153[k] = f_8 * hsf_103[k]
                   + f_4 * isd0_63[k]
                   - f_5 * isd1_63[k]
                   + f_3 * pc_x[k] * isf_103[k];

        t_154[k] = f_3 * pc_z[k] * isf_101[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, pc_x, pc_z, hsf_106, hsf_108, isd0_60, \
                         isd1_60, isf_102, isf_103, isf_106, isf_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_4 * isd0_60[k]
                   - f_5 * isd1_60[k]
                   + f_3 * pc_z[k] * isf_102[k];

        t_156[k] = f_8 * hsf_106[k]
                   + f_3 * pc_x[k] * isf_106[k];

        t_157[k] = f_3 * pc_z[k] * isf_103[k];

        t_158[k] = f_8 * hsf_108[k]
                   + f_3 * pc_x[k] * isf_108[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, t_163, pc_x, pc_y, pc_z, hsf_66, hsf_69, \
                         hsf_109, isd0_63, isd1_63, isf_106, isf_107, \
                         isf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_8 * hsf_109[k]
                   + f_3 * pc_x[k] * isf_109[k];

        t_160[k] = f_12 * hsf_66[k]
                   + f_1 * isd0_63[k]
                   - f_2 * isd1_63[k]
                   + f_3 * pc_y[k] * isf_106[k];

        t_161[k] = f_3 * pc_z[k] * isf_106[k];

        t_162[k] = f_4 * isd0_63[k]
                   - f_5 * isd1_63[k]
                   + f_3 * pc_z[k] * isf_107[k];

        t_163[k] = f_12 * hsf_69[k]
                   + f_3 * pc_y[k] * isf_109[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, pa_z, pc_y, pc_z, hsg0_90, hsf_60, \
                         hsf_70, hsg1_90, isd0_65, isd1_65, isf_109, \
                         isf_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_1 * isd0_65[k]
                   - f_2 * isd1_65[k]
                   + f_3 * pc_z[k] * isf_109[k];

        t_165[k] = pa_z[k] * hsg0_90[k]
                   - f_6 * pc_z[k] * hsg1_90[k];

        t_166[k] = f_13 * hsf_70[k]
                   + f_3 * pc_y[k] * isf_110[k];

        t_167[k] = f_7 * hsf_60[k]
                   + f_3 * pc_z[k] * isf_110[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, pa_z, pc_x, pc_y, pc_z, hsg0_93, hsf_72, \
                         hsf_115, hsg1_93, isd0_71, isd1_71, isf_112, \
                         isf_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = pa_z[k] * hsg0_93[k]
                   - f_6 * pc_z[k] * hsg1_93[k];

        t_169[k] = f_13 * hsf_72[k]
                   + f_3 * pc_y[k] * isf_112[k];

        t_170[k] = f_8 * hsf_115[k]
                   + f_4 * isd0_71[k]
                   - f_5 * isd1_71[k]
                   + f_3 * pc_x[k] * isf_115[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, t_174, pc_x, hsf_116, hsf_117, hsf_118, hsf_119, \
                         isf_116, isf_117, isf_118, isf_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = f_8 * hsf_116[k]
                   + f_3 * pc_x[k] * isf_116[k];

        t_172[k] = f_8 * hsf_117[k]
                   + f_3 * pc_x[k] * isf_117[k];

        t_173[k] = f_8 * hsf_118[k]
                   + f_3 * pc_x[k] * isf_118[k];

        t_174[k] = f_8 * hsf_119[k]
                   + f_3 * pc_x[k] * isf_119[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, pa_z, pc_y, pc_z, hsg0_100, hsf_66, hsf_78, \
                         hsg1_100, isd0_71, isd1_71, isf_116, isf_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = pa_z[k] * hsg0_100[k]
                   - f_6 * pc_z[k] * hsg1_100[k];

        t_176[k] = f_7 * hsf_66[k]
                   + f_3 * pc_z[k] * isf_116[k];

        t_177[k] = f_13 * hsf_78[k]
                   + f_4 * isd0_71[k]
                   - f_5 * isd1_71[k]
                   + f_3 * pc_y[k] * isf_118[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, pc_x, pc_y, pc_z, hsf_69, hsf_79, hsf_120, \
                         isd0_71, isd0_72, isd1_71, isd1_72, isf_119, \
                         isf_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_13 * hsf_79[k]
                   + f_3 * pc_y[k] * isf_119[k];

        t_179[k] = f_7 * hsf_69[k]
                   + f_1 * isd0_71[k]
                   - f_2 * isd1_71[k]
                   + f_3 * pc_z[k] * isf_119[k];

        t_180[k] = f_8 * hsf_120[k]
                   + f_1 * isd0_72[k]
                   - f_2 * isd1_72[k]
                   + f_3 * pc_x[k] * isf_120[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pc_x, pc_y, pc_z, hsf_70, hsf_80, hsf_82, \
                         hsf_123, isd0_75, isd1_75, isf_120, isf_122, \
                         isf_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_8 * hsf_80[k]
                   + f_3 * pc_y[k] * isf_120[k];

        t_182[k] = f_8 * hsf_70[k]
                   + f_3 * pc_z[k] * isf_120[k];

        t_183[k] = f_8 * hsf_123[k]
                   + f_4 * isd0_75[k]
                   - f_5 * isd1_75[k]
                   + f_3 * pc_x[k] * isf_123[k];

        t_184[k] = f_8 * hsf_82[k]
                   + f_3 * pc_y[k] * isf_122[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pc_x, hsf_125, hsf_126, hsf_127, hsf_128, \
                         isd0_77, isd1_77, isf_125, isf_126, isf_127, \
                         isf_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_8 * hsf_125[k]
                   + f_4 * isd0_77[k]
                   - f_5 * isd1_77[k]
                   + f_3 * pc_x[k] * isf_125[k];

        t_186[k] = f_8 * hsf_126[k]
                   + f_3 * pc_x[k] * isf_126[k];

        t_187[k] = f_8 * hsf_127[k]
                   + f_3 * pc_x[k] * isf_127[k];

        t_188[k] = f_8 * hsf_128[k]
                   + f_3 * pc_x[k] * isf_128[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, pc_x, pc_y, pc_z, hsf_76, hsf_86, hsf_129, \
                         isd0_75, isd1_75, isf_126, isf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_8 * hsf_129[k]
                   + f_3 * pc_x[k] * isf_129[k];

        t_190[k] = f_8 * hsf_86[k]
                   + f_1 * isd0_75[k]
                   - f_2 * isd1_75[k]
                   + f_3 * pc_y[k] * isf_126[k];

        t_191[k] = f_8 * hsf_76[k]
                   + f_3 * pc_z[k] * isf_126[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, pa_y, pc_y, pc_z, hsg0_135, hsf_79, \
                         hsf_88, hsf_89, hsg1_135, isd0_77, isd1_77, isf_128, \
                         isf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_8 * hsf_88[k]
                   + f_4 * isd0_77[k]
                   - f_5 * isd1_77[k]
                   + f_3 * pc_y[k] * isf_128[k];

        t_193[k] = f_8 * hsf_89[k]
                   + f_3 * pc_y[k] * isf_129[k];

        t_194[k] = f_8 * hsf_79[k]
                   + f_1 * isd0_77[k]
                   - f_2 * isd1_77[k]
                   + f_3 * pc_z[k] * isf_129[k];

        t_195[k] = pa_y[k] * hsg0_135[k]
                   - f_6 * pc_y[k] * hsg1_135[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, pa_y, pc_y, pc_z, hsg0_138, hsf_80, \
                         hsf_90, hsf_91, hsf_92, hsg1_138, isf_130, \
                         isf_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = f_7 * hsf_90[k]
                   + f_3 * pc_y[k] * isf_130[k];

        t_197[k] = f_13 * hsf_80[k]
                   + f_3 * pc_z[k] * isf_130[k];

        t_198[k] = pa_y[k] * hsg0_138[k]
                   + f_8 * hsf_91[k]
                   - f_6 * pc_y[k] * hsg1_138[k];

        t_199[k] = f_7 * hsf_92[k]
                   + f_3 * pc_y[k] * isf_132[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, pa_y, pc_x, pc_y, hsg0_140, hsf_136, \
                         hsf_137, hsf_138, hsg1_140, isf_136, isf_137, \
                         isf_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = pa_y[k] * hsg0_140[k]
                   - f_6 * pc_y[k] * hsg1_140[k];

        t_201[k] = f_8 * hsf_136[k]
                   + f_3 * pc_x[k] * isf_136[k];

        t_202[k] = f_8 * hsf_137[k]
                   + f_3 * pc_x[k] * isf_137[k];

        t_203[k] = f_8 * hsf_138[k]
                   + f_3 * pc_x[k] * isf_138[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, pc_x, pc_y, pc_z, hsf_86, hsf_96, hsf_139, \
                         isd0_81, isd1_81, isf_136, isf_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_8 * hsf_139[k]
                   + f_3 * pc_x[k] * isf_139[k];

        t_205[k] = f_7 * hsf_96[k]
                   + f_1 * isd0_81[k]
                   - f_2 * isd1_81[k]
                   + f_3 * pc_y[k] * isf_136[k];

        t_206[k] = f_13 * hsf_86[k]
                   + f_3 * pc_z[k] * isf_136[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, pa_y, pc_y, hsg0_149, hsf_98, hsf_99, hsg1_149, \
                         isd0_83, isd1_83, isf_138, isf_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = f_7 * hsf_98[k]
                   + f_4 * isd0_83[k]
                   - f_5 * isd1_83[k]
                   + f_3 * pc_y[k] * isf_138[k];

        t_208[k] = f_7 * hsf_99[k]
                   + f_3 * pc_y[k] * isf_139[k];

        t_209[k] = pa_y[k] * hsg0_149[k]
                   - f_6 * pc_y[k] * hsg1_149[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, pc_x, pc_y, pc_z, hsf_90, hsf_140, \
                         isd0_84, isd1_84, isf_140, isf_141, isf_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_8 * hsf_140[k]
                   + f_1 * isd0_84[k]
                   - f_2 * isd1_84[k]
                   + f_3 * pc_x[k] * isf_140[k];

        t_211[k] = f_3 * pc_y[k] * isf_140[k];

        t_212[k] = f_12 * hsf_90[k]
                   + f_3 * pc_z[k] * isf_140[k];

        t_213[k] = f_4 * isd0_84[k]
                   - f_5 * isd1_84[k]
                   + f_3 * pc_y[k] * isf_141[k];

        t_214[k] = f_3 * pc_y[k] * isf_142[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, pc_x, pc_y, hsf_145, hsf_146, hsf_147, \
                         isd0_89, isd1_89, isf_145, isf_146, isf_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = f_8 * hsf_145[k]
                   + f_4 * isd0_89[k]
                   - f_5 * isd1_89[k]
                   + f_3 * pc_x[k] * isf_145[k];

        t_216[k] = f_8 * hsf_146[k]
                   + f_3 * pc_x[k] * isf_146[k];

        t_217[k] = f_8 * hsf_147[k]
                   + f_3 * pc_x[k] * isf_147[k];

        t_218[k] = f_3 * pc_y[k] * isf_145[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, pc_x, pc_y, hsf_149, isd0_87, isd0_88, isd1_87, \
                         isd1_88, isf_146, isf_147, isf_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = f_8 * hsf_149[k]
                   + f_3 * pc_x[k] * isf_149[k];

        t_220[k] = f_1 * isd0_87[k]
                   - f_2 * isd1_87[k]
                   + f_3 * pc_y[k] * isf_146[k];

        t_221[k] = f_10 * isd0_88[k]
                   - f_11 * isd1_88[k]
                   + f_3 * pc_y[k] * isf_147[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, pa_x, pc_x, pc_y, pc_z, hsg0_225, hsf_99, \
                         hsf_150, hsg1_225, isd0_89, isd1_89, isf_148, \
                         isf_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_4 * isd0_89[k]
                   - f_5 * isd1_89[k]
                   + f_3 * pc_y[k] * isf_148[k];

        t_223[k] = f_3 * pc_y[k] * isf_149[k];

        t_224[k] = f_12 * hsf_99[k]
                   + f_1 * isd0_89[k]
                   - f_2 * isd1_89[k]
                   + f_3 * pc_z[k] * isf_149[k];

        t_225[k] = pa_x[k] * hsg0_225[k]
                   + f_12 * hsf_150[k]
                   - f_6 * pc_x[k] * hsg1_225[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, t_229, pa_x, pc_x, pc_y, pc_z, hsg0_228, \
                         hsf_100, hsf_153, hsg1_228, isf_150, isf_151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = f_9 * hsf_100[k]
                   + f_3 * pc_y[k] * isf_150[k];

        t_227[k] = f_3 * pc_z[k] * isf_150[k];

        t_228[k] = pa_x[k] * hsg0_228[k]
                   + f_8 * hsf_153[k]
                   - f_6 * pc_x[k] * hsg1_228[k];

        t_229[k] = f_3 * pc_z[k] * isf_151[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, pc_x, pc_z, hsf_156, hsf_158, isd0_90, \
                         isd1_90, isf_152, isf_153, isf_156, isf_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_4 * isd0_90[k]
                   - f_5 * isd1_90[k]
                   + f_3 * pc_z[k] * isf_152[k];

        t_231[k] = f_7 * hsf_156[k]
                   + f_3 * pc_x[k] * isf_156[k];

        t_232[k] = f_3 * pc_z[k] * isf_153[k];

        t_233[k] = f_7 * hsf_158[k]
                   + f_3 * pc_x[k] * isf_158[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, pa_x, pc_x, pc_z, hsg0_235, hsg0_237, \
                         hsf_159, hsg1_235, hsg1_237, isf_156, \
                         isf_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_7 * hsf_159[k]
                   + f_3 * pc_x[k] * isf_159[k];

        t_235[k] = pa_x[k] * hsg0_235[k]
                   - f_6 * pc_x[k] * hsg1_235[k];

        t_236[k] = f_3 * pc_z[k] * isf_156[k];

        t_237[k] = pa_x[k] * hsg0_237[k]
                   - f_6 * pc_x[k] * hsg1_237[k];
    }

#pragma omp simd aligned(t_238, t_239, t_240, pa_x, pa_z, pc_x, pc_y, pc_z, hsg0_150, \
                         hsg0_239, hsf_109, hsg1_150, hsg1_239, \
                         isf_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_238[k] = f_9 * hsf_109[k]
                   + f_3 * pc_y[k] * isf_159[k];

        t_239[k] = pa_x[k] * hsg0_239[k]
                   - f_6 * pc_x[k] * hsg1_239[k];

        t_240[k] = pa_z[k] * hsg0_150[k]
                   - f_6 * pc_z[k] * hsg1_150[k];
    }

#pragma omp simd aligned(t_241, t_242, t_243, t_244, pa_z, pc_y, pc_z, hsg0_153, hsf_100, \
                         hsf_110, hsf_112, hsg1_153, isf_160, isf_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_241[k] = f_12 * hsf_110[k]
                   + f_3 * pc_y[k] * isf_160[k];

        t_242[k] = f_7 * hsf_100[k]
                   + f_3 * pc_z[k] * isf_160[k];

        t_243[k] = pa_z[k] * hsg0_153[k]
                   - f_6 * pc_z[k] * hsg1_153[k];

        t_244[k] = f_12 * hsf_112[k]
                   + f_3 * pc_y[k] * isf_162[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, pa_x, pc_x, hsg0_245, hsf_165, hsf_166, \
                         hsf_167, hsf_168, hsg1_245, isf_166, isf_167, \
                         isf_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = pa_x[k] * hsg0_245[k]
                   + f_8 * hsf_165[k]
                   - f_6 * pc_x[k] * hsg1_245[k];

        t_246[k] = f_7 * hsf_166[k]
                   + f_3 * pc_x[k] * isf_166[k];

        t_247[k] = f_7 * hsf_167[k]
                   + f_3 * pc_x[k] * isf_167[k];

        t_248[k] = f_7 * hsf_168[k]
                   + f_3 * pc_x[k] * isf_168[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, t_252, pa_x, pc_x, pc_z, hsg0_250, hsg0_252, \
                         hsf_106, hsf_169, hsg1_250, hsg1_252, isf_166, \
                         isf_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = f_7 * hsf_169[k]
                   + f_3 * pc_x[k] * isf_169[k];

        t_250[k] = pa_x[k] * hsg0_250[k]
                   - f_6 * pc_x[k] * hsg1_250[k];

        t_251[k] = f_7 * hsf_106[k]
                   + f_3 * pc_z[k] * isf_166[k];

        t_252[k] = pa_x[k] * hsg0_252[k]
                   - f_6 * pc_x[k] * hsg1_252[k];
    }
}

static auto
compute_prim_isg_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t hsg0,
                                                          const size_t hsf, const size_t hsg1,
                                                          const size_t isd0, const size_t isd1,
                                                          const size_t isf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 2.5 / q;
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 2.0 / q;
    const auto f_13 = 1.5 / q;

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
    auto *t_370 = buffer.data(target + 370);
    auto *t_371 = buffer.data(target + 371);
    auto *t_372 = buffer.data(target + 372);
    auto *t_373 = buffer.data(target + 373);
    auto *t_374 = buffer.data(target + 374);
    auto *t_375 = buffer.data(target + 375);
    auto *t_376 = buffer.data(target + 376);
    auto *t_377 = buffer.data(target + 377);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *hsg0_210 = buffer.data(hsg0 + 210);
    const auto *hsg0_215 = buffer.data(hsg0 + 215);
    const auto *hsg0_225 = buffer.data(hsg0 + 225);
    const auto *hsg0_226 = buffer.data(hsg0 + 226);
    const auto *hsg0_228 = buffer.data(hsg0 + 228);
    const auto *hsg0_235 = buffer.data(hsg0 + 235);
    const auto *hsg0_237 = buffer.data(hsg0 + 237);
    const auto *hsg0_254 = buffer.data(hsg0 + 254);
    const auto *hsg0_255 = buffer.data(hsg0 + 255);
    const auto *hsg0_258 = buffer.data(hsg0 + 258);
    const auto *hsg0_260 = buffer.data(hsg0 + 260);
    const auto *hsg0_265 = buffer.data(hsg0 + 265);
    const auto *hsg0_267 = buffer.data(hsg0 + 267);
    const auto *hsg0_269 = buffer.data(hsg0 + 269);
    const auto *hsg0_270 = buffer.data(hsg0 + 270);
    const auto *hsg0_273 = buffer.data(hsg0 + 273);
    const auto *hsg0_275 = buffer.data(hsg0 + 275);
    const auto *hsg0_280 = buffer.data(hsg0 + 280);
    const auto *hsg0_282 = buffer.data(hsg0 + 282);
    const auto *hsg0_284 = buffer.data(hsg0 + 284);
    const auto *hsg0_288 = buffer.data(hsg0 + 288);
    const auto *hsg0_295 = buffer.data(hsg0 + 295);
    const auto *hsg0_297 = buffer.data(hsg0 + 297);
    const auto *hsg0_299 = buffer.data(hsg0 + 299);
    const auto *hsg0_300 = buffer.data(hsg0 + 300);
    const auto *hsg0_305 = buffer.data(hsg0 + 305);
    const auto *hsg0_310 = buffer.data(hsg0 + 310);
    const auto *hsg0_311 = buffer.data(hsg0 + 311);
    const auto *hsg0_312 = buffer.data(hsg0 + 312);
    const auto *hsg0_314 = buffer.data(hsg0 + 314);

    const auto *hsf_110 = buffer.data(hsf + 110);
    const auto *hsf_116 = buffer.data(hsf + 116);
    const auto *hsf_119 = buffer.data(hsf + 119);
    const auto *hsf_120 = buffer.data(hsf + 120);
    const auto *hsf_122 = buffer.data(hsf + 122);
    const auto *hsf_126 = buffer.data(hsf + 126);
    const auto *hsf_129 = buffer.data(hsf + 129);
    const auto *hsf_130 = buffer.data(hsf + 130);
    const auto *hsf_132 = buffer.data(hsf + 132);
    const auto *hsf_136 = buffer.data(hsf + 136);
    const auto *hsf_139 = buffer.data(hsf + 139);
    const auto *hsf_140 = buffer.data(hsf + 140);
    const auto *hsf_142 = buffer.data(hsf + 142);
    const auto *hsf_149 = buffer.data(hsf + 149);
    const auto *hsf_156 = buffer.data(hsf + 156);
    const auto *hsf_157 = buffer.data(hsf + 157);
    const auto *hsf_159 = buffer.data(hsf + 159);
    const auto *hsf_166 = buffer.data(hsf + 166);
    const auto *hsf_169 = buffer.data(hsf + 169);
    const auto *hsf_170 = buffer.data(hsf + 170);
    const auto *hsf_173 = buffer.data(hsf + 173);
    const auto *hsf_175 = buffer.data(hsf + 175);
    const auto *hsf_176 = buffer.data(hsf + 176);
    const auto *hsf_177 = buffer.data(hsf + 177);
    const auto *hsf_178 = buffer.data(hsf + 178);
    const auto *hsf_179 = buffer.data(hsf + 179);
    const auto *hsf_180 = buffer.data(hsf + 180);
    const auto *hsf_183 = buffer.data(hsf + 183);
    const auto *hsf_185 = buffer.data(hsf + 185);
    const auto *hsf_186 = buffer.data(hsf + 186);
    const auto *hsf_187 = buffer.data(hsf + 187);
    const auto *hsf_188 = buffer.data(hsf + 188);
    const auto *hsf_189 = buffer.data(hsf + 189);
    const auto *hsf_193 = buffer.data(hsf + 193);
    const auto *hsf_196 = buffer.data(hsf + 196);
    const auto *hsf_197 = buffer.data(hsf + 197);
    const auto *hsf_198 = buffer.data(hsf + 198);
    const auto *hsf_199 = buffer.data(hsf + 199);
    const auto *hsf_200 = buffer.data(hsf + 200);
    const auto *hsf_205 = buffer.data(hsf + 205);
    const auto *hsf_206 = buffer.data(hsf + 206);
    const auto *hsf_207 = buffer.data(hsf + 207);
    const auto *hsf_209 = buffer.data(hsf + 209);

    const auto *hsg1_210 = buffer.data(hsg1 + 210);
    const auto *hsg1_215 = buffer.data(hsg1 + 215);
    const auto *hsg1_225 = buffer.data(hsg1 + 225);
    const auto *hsg1_226 = buffer.data(hsg1 + 226);
    const auto *hsg1_228 = buffer.data(hsg1 + 228);
    const auto *hsg1_235 = buffer.data(hsg1 + 235);
    const auto *hsg1_237 = buffer.data(hsg1 + 237);
    const auto *hsg1_254 = buffer.data(hsg1 + 254);
    const auto *hsg1_255 = buffer.data(hsg1 + 255);
    const auto *hsg1_258 = buffer.data(hsg1 + 258);
    const auto *hsg1_260 = buffer.data(hsg1 + 260);
    const auto *hsg1_265 = buffer.data(hsg1 + 265);
    const auto *hsg1_267 = buffer.data(hsg1 + 267);
    const auto *hsg1_269 = buffer.data(hsg1 + 269);
    const auto *hsg1_270 = buffer.data(hsg1 + 270);
    const auto *hsg1_273 = buffer.data(hsg1 + 273);
    const auto *hsg1_275 = buffer.data(hsg1 + 275);
    const auto *hsg1_280 = buffer.data(hsg1 + 280);
    const auto *hsg1_282 = buffer.data(hsg1 + 282);
    const auto *hsg1_284 = buffer.data(hsg1 + 284);
    const auto *hsg1_288 = buffer.data(hsg1 + 288);
    const auto *hsg1_295 = buffer.data(hsg1 + 295);
    const auto *hsg1_297 = buffer.data(hsg1 + 297);
    const auto *hsg1_299 = buffer.data(hsg1 + 299);
    const auto *hsg1_300 = buffer.data(hsg1 + 300);
    const auto *hsg1_305 = buffer.data(hsg1 + 305);
    const auto *hsg1_310 = buffer.data(hsg1 + 310);
    const auto *hsg1_311 = buffer.data(hsg1 + 311);
    const auto *hsg1_312 = buffer.data(hsg1 + 312);
    const auto *hsg1_314 = buffer.data(hsg1 + 314);

    const auto *isd0_120 = buffer.data(isd0 + 120);
    const auto *isd0_126 = buffer.data(isd0 + 126);
    const auto *isd0_127 = buffer.data(isd0 + 127);
    const auto *isd0_129 = buffer.data(isd0 + 129);
    const auto *isd0_131 = buffer.data(isd0 + 131);
    const auto *isd0_134 = buffer.data(isd0 + 134);
    const auto *isd0_136 = buffer.data(isd0 + 136);
    const auto *isd0_137 = buffer.data(isd0 + 137);
    const auto *isd0_138 = buffer.data(isd0 + 138);
    const auto *isd0_139 = buffer.data(isd0 + 139);
    const auto *isd0_140 = buffer.data(isd0 + 140);
    const auto *isd0_141 = buffer.data(isd0 + 141);
    const auto *isd0_142 = buffer.data(isd0 + 142);
    const auto *isd0_143 = buffer.data(isd0 + 143);
    const auto *isd0_144 = buffer.data(isd0 + 144);
    const auto *isd0_145 = buffer.data(isd0 + 145);
    const auto *isd0_146 = buffer.data(isd0 + 146);
    const auto *isd0_147 = buffer.data(isd0 + 147);
    const auto *isd0_148 = buffer.data(isd0 + 148);
    const auto *isd0_149 = buffer.data(isd0 + 149);
    const auto *isd0_150 = buffer.data(isd0 + 150);
    const auto *isd0_151 = buffer.data(isd0 + 151);
    const auto *isd0_152 = buffer.data(isd0 + 152);

    const auto *isd1_120 = buffer.data(isd1 + 120);
    const auto *isd1_126 = buffer.data(isd1 + 126);
    const auto *isd1_127 = buffer.data(isd1 + 127);
    const auto *isd1_129 = buffer.data(isd1 + 129);
    const auto *isd1_131 = buffer.data(isd1 + 131);
    const auto *isd1_134 = buffer.data(isd1 + 134);
    const auto *isd1_136 = buffer.data(isd1 + 136);
    const auto *isd1_137 = buffer.data(isd1 + 137);
    const auto *isd1_138 = buffer.data(isd1 + 138);
    const auto *isd1_139 = buffer.data(isd1 + 139);
    const auto *isd1_140 = buffer.data(isd1 + 140);
    const auto *isd1_141 = buffer.data(isd1 + 141);
    const auto *isd1_142 = buffer.data(isd1 + 142);
    const auto *isd1_143 = buffer.data(isd1 + 143);
    const auto *isd1_144 = buffer.data(isd1 + 144);
    const auto *isd1_145 = buffer.data(isd1 + 145);
    const auto *isd1_146 = buffer.data(isd1 + 146);
    const auto *isd1_147 = buffer.data(isd1 + 147);
    const auto *isd1_148 = buffer.data(isd1 + 148);
    const auto *isd1_149 = buffer.data(isd1 + 149);
    const auto *isd1_150 = buffer.data(isd1 + 150);
    const auto *isd1_151 = buffer.data(isd1 + 151);
    const auto *isd1_152 = buffer.data(isd1 + 152);

    const auto *isf_169 = buffer.data(isf + 169);
    const auto *isf_170 = buffer.data(isf + 170);
    const auto *isf_172 = buffer.data(isf + 172);
    const auto *isf_176 = buffer.data(isf + 176);
    const auto *isf_177 = buffer.data(isf + 177);
    const auto *isf_178 = buffer.data(isf + 178);
    const auto *isf_179 = buffer.data(isf + 179);
    const auto *isf_180 = buffer.data(isf + 180);
    const auto *isf_182 = buffer.data(isf + 182);
    const auto *isf_186 = buffer.data(isf + 186);
    const auto *isf_187 = buffer.data(isf + 187);
    const auto *isf_188 = buffer.data(isf + 188);
    const auto *isf_189 = buffer.data(isf + 189);
    const auto *isf_190 = buffer.data(isf + 190);
    const auto *isf_192 = buffer.data(isf + 192);
    const auto *isf_196 = buffer.data(isf + 196);
    const auto *isf_197 = buffer.data(isf + 197);
    const auto *isf_198 = buffer.data(isf + 198);
    const auto *isf_199 = buffer.data(isf + 199);
    const auto *isf_200 = buffer.data(isf + 200);
    const auto *isf_201 = buffer.data(isf + 201);
    const auto *isf_202 = buffer.data(isf + 202);
    const auto *isf_205 = buffer.data(isf + 205);
    const auto *isf_206 = buffer.data(isf + 206);
    const auto *isf_207 = buffer.data(isf + 207);
    const auto *isf_209 = buffer.data(isf + 209);
    const auto *isf_210 = buffer.data(isf + 210);
    const auto *isf_211 = buffer.data(isf + 211);
    const auto *isf_213 = buffer.data(isf + 213);
    const auto *isf_215 = buffer.data(isf + 215);
    const auto *isf_216 = buffer.data(isf + 216);
    const auto *isf_217 = buffer.data(isf + 217);
    const auto *isf_218 = buffer.data(isf + 218);
    const auto *isf_219 = buffer.data(isf + 219);
    const auto *isf_222 = buffer.data(isf + 222);
    const auto *isf_224 = buffer.data(isf + 224);
    const auto *isf_225 = buffer.data(isf + 225);
    const auto *isf_226 = buffer.data(isf + 226);
    const auto *isf_227 = buffer.data(isf + 227);
    const auto *isf_228 = buffer.data(isf + 228);
    const auto *isf_229 = buffer.data(isf + 229);
    const auto *isf_230 = buffer.data(isf + 230);
    const auto *isf_231 = buffer.data(isf + 231);
    const auto *isf_232 = buffer.data(isf + 232);
    const auto *isf_233 = buffer.data(isf + 233);
    const auto *isf_234 = buffer.data(isf + 234);
    const auto *isf_235 = buffer.data(isf + 235);
    const auto *isf_236 = buffer.data(isf + 236);
    const auto *isf_237 = buffer.data(isf + 237);
    const auto *isf_238 = buffer.data(isf + 238);
    const auto *isf_239 = buffer.data(isf + 239);
    const auto *isf_240 = buffer.data(isf + 240);
    const auto *isf_241 = buffer.data(isf + 241);
    const auto *isf_242 = buffer.data(isf + 242);
    const auto *isf_243 = buffer.data(isf + 243);
    const auto *isf_244 = buffer.data(isf + 244);
    const auto *isf_245 = buffer.data(isf + 245);
    const auto *isf_246 = buffer.data(isf + 246);
    const auto *isf_247 = buffer.data(isf + 247);
    const auto *isf_248 = buffer.data(isf + 248);
    const auto *isf_249 = buffer.data(isf + 249);
    const auto *isf_250 = buffer.data(isf + 250);
    const auto *isf_251 = buffer.data(isf + 251);
    const auto *isf_252 = buffer.data(isf + 252);

#pragma omp simd aligned(t_253, t_254, t_255, t_256, pa_x, pc_x, pc_y, hsg0_254, hsg0_255, \
                         hsf_119, hsf_120, hsf_170, hsg1_254, hsg1_255, isf_169, \
                         isf_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_253[k] = f_12 * hsf_119[k]
                   + f_3 * pc_y[k] * isf_169[k];

        t_254[k] = pa_x[k] * hsg0_254[k]
                   - f_6 * pc_x[k] * hsg1_254[k];

        t_255[k] = pa_x[k] * hsg0_255[k]
                   + f_12 * hsf_170[k]
                   - f_6 * pc_x[k] * hsg1_255[k];

        t_256[k] = f_13 * hsf_120[k]
                   + f_3 * pc_y[k] * isf_170[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, pa_x, pc_x, pc_y, pc_z, hsg0_258, hsf_110, \
                         hsf_122, hsf_173, hsg1_258, isf_170, isf_172 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_8 * hsf_110[k]
                   + f_3 * pc_z[k] * isf_170[k];

        t_258[k] = pa_x[k] * hsg0_258[k]
                   + f_8 * hsf_173[k]
                   - f_6 * pc_x[k] * hsg1_258[k];

        t_259[k] = f_13 * hsf_122[k]
                   + f_3 * pc_y[k] * isf_172[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pa_x, pc_x, hsg0_260, hsf_175, hsf_176, \
                         hsf_177, hsf_178, hsg1_260, isf_176, isf_177, \
                         isf_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = pa_x[k] * hsg0_260[k]
                   + f_8 * hsf_175[k]
                   - f_6 * pc_x[k] * hsg1_260[k];

        t_261[k] = f_7 * hsf_176[k]
                   + f_3 * pc_x[k] * isf_176[k];

        t_262[k] = f_7 * hsf_177[k]
                   + f_3 * pc_x[k] * isf_177[k];

        t_263[k] = f_7 * hsf_178[k]
                   + f_3 * pc_x[k] * isf_178[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, pa_x, pc_x, pc_z, hsg0_265, hsg0_267, \
                         hsf_116, hsf_179, hsg1_265, hsg1_267, isf_176, \
                         isf_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_7 * hsf_179[k]
                   + f_3 * pc_x[k] * isf_179[k];

        t_265[k] = pa_x[k] * hsg0_265[k]
                   - f_6 * pc_x[k] * hsg1_265[k];

        t_266[k] = f_8 * hsf_116[k]
                   + f_3 * pc_z[k] * isf_176[k];

        t_267[k] = pa_x[k] * hsg0_267[k]
                   - f_6 * pc_x[k] * hsg1_267[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, t_271, pa_x, pc_x, pc_y, hsg0_269, hsg0_270, \
                         hsf_129, hsf_130, hsf_180, hsg1_269, hsg1_270, isf_179, \
                         isf_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_13 * hsf_129[k]
                   + f_3 * pc_y[k] * isf_179[k];

        t_269[k] = pa_x[k] * hsg0_269[k]
                   - f_6 * pc_x[k] * hsg1_269[k];

        t_270[k] = pa_x[k] * hsg0_270[k]
                   + f_12 * hsf_180[k]
                   - f_6 * pc_x[k] * hsg1_270[k];

        t_271[k] = f_8 * hsf_130[k]
                   + f_3 * pc_y[k] * isf_180[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, pa_x, pc_x, pc_y, pc_z, hsg0_273, hsf_120, \
                         hsf_132, hsf_183, hsg1_273, isf_180, isf_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_13 * hsf_120[k]
                   + f_3 * pc_z[k] * isf_180[k];

        t_273[k] = pa_x[k] * hsg0_273[k]
                   + f_8 * hsf_183[k]
                   - f_6 * pc_x[k] * hsg1_273[k];

        t_274[k] = f_8 * hsf_132[k]
                   + f_3 * pc_y[k] * isf_182[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, pa_x, pc_x, hsg0_275, hsf_185, hsf_186, \
                         hsf_187, hsf_188, hsg1_275, isf_186, isf_187, \
                         isf_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = pa_x[k] * hsg0_275[k]
                   + f_8 * hsf_185[k]
                   - f_6 * pc_x[k] * hsg1_275[k];

        t_276[k] = f_7 * hsf_186[k]
                   + f_3 * pc_x[k] * isf_186[k];

        t_277[k] = f_7 * hsf_187[k]
                   + f_3 * pc_x[k] * isf_187[k];

        t_278[k] = f_7 * hsf_188[k]
                   + f_3 * pc_x[k] * isf_188[k];
    }

#pragma omp simd aligned(t_279, t_280, t_281, t_282, pa_x, pc_x, pc_z, hsg0_280, hsg0_282, \
                         hsf_126, hsf_189, hsg1_280, hsg1_282, isf_186, \
                         isf_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_279[k] = f_7 * hsf_189[k]
                   + f_3 * pc_x[k] * isf_189[k];

        t_280[k] = pa_x[k] * hsg0_280[k]
                   - f_6 * pc_x[k] * hsg1_280[k];

        t_281[k] = f_13 * hsf_126[k]
                   + f_3 * pc_z[k] * isf_186[k];

        t_282[k] = pa_x[k] * hsg0_282[k]
                   - f_6 * pc_x[k] * hsg1_282[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, t_286, pa_x, pa_y, pc_x, pc_y, hsg0_210, \
                         hsg0_284, hsf_139, hsf_140, hsg1_210, hsg1_284, isf_189, \
                         isf_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_8 * hsf_139[k]
                   + f_3 * pc_y[k] * isf_189[k];

        t_284[k] = pa_x[k] * hsg0_284[k]
                   - f_6 * pc_x[k] * hsg1_284[k];

        t_285[k] = pa_y[k] * hsg0_210[k]
                   - f_6 * pc_y[k] * hsg1_210[k];

        t_286[k] = f_7 * hsf_140[k]
                   + f_3 * pc_y[k] * isf_190[k];
    }

#pragma omp simd aligned(t_287, t_288, t_289, pa_x, pc_x, pc_y, pc_z, hsg0_288, hsf_130, \
                         hsf_142, hsf_193, hsg1_288, isf_190, isf_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_287[k] = f_12 * hsf_130[k]
                   + f_3 * pc_z[k] * isf_190[k];

        t_288[k] = pa_x[k] * hsg0_288[k]
                   + f_8 * hsf_193[k]
                   - f_6 * pc_x[k] * hsg1_288[k];

        t_289[k] = f_7 * hsf_142[k]
                   + f_3 * pc_y[k] * isf_192[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, pa_y, pc_x, pc_y, hsg0_215, hsf_196, \
                         hsf_197, hsf_198, hsg1_215, isf_196, isf_197, \
                         isf_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = pa_y[k] * hsg0_215[k]
                   - f_6 * pc_y[k] * hsg1_215[k];

        t_291[k] = f_7 * hsf_196[k]
                   + f_3 * pc_x[k] * isf_196[k];

        t_292[k] = f_7 * hsf_197[k]
                   + f_3 * pc_x[k] * isf_197[k];

        t_293[k] = f_7 * hsf_198[k]
                   + f_3 * pc_x[k] * isf_198[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, pa_x, pc_x, pc_z, hsg0_295, hsg0_297, \
                         hsf_136, hsf_199, hsg1_295, hsg1_297, isf_196, \
                         isf_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_7 * hsf_199[k]
                   + f_3 * pc_x[k] * isf_199[k];

        t_295[k] = pa_x[k] * hsg0_295[k]
                   - f_6 * pc_x[k] * hsg1_295[k];

        t_296[k] = f_12 * hsf_136[k]
                   + f_3 * pc_z[k] * isf_196[k];

        t_297[k] = pa_x[k] * hsg0_297[k]
                   - f_6 * pc_x[k] * hsg1_297[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, t_301, pa_x, pc_x, pc_y, hsg0_299, hsg0_300, \
                         hsf_149, hsf_200, hsg1_299, hsg1_300, isf_199, \
                         isf_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = f_7 * hsf_149[k]
                   + f_3 * pc_y[k] * isf_199[k];

        t_299[k] = pa_x[k] * hsg0_299[k]
                   - f_6 * pc_x[k] * hsg1_299[k];

        t_300[k] = pa_x[k] * hsg0_300[k]
                   + f_12 * hsf_200[k]
                   - f_6 * pc_x[k] * hsg1_300[k];

        t_301[k] = f_3 * pc_y[k] * isf_200[k];
    }

#pragma omp simd aligned(t_302, t_303, t_304, pc_y, pc_z, hsf_140, isd0_120, isd1_120, \
                         isf_200, isf_201, isf_202 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_302[k] = f_9 * hsf_140[k]
                   + f_3 * pc_z[k] * isf_200[k];

        t_303[k] = f_4 * isd0_120[k]
                   - f_5 * isd1_120[k]
                   + f_3 * pc_y[k] * isf_201[k];

        t_304[k] = f_3 * pc_y[k] * isf_202[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, pa_x, pc_x, pc_y, hsg0_305, hsf_205, \
                         hsf_206, hsf_207, hsg1_305, isf_205, isf_206, \
                         isf_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = pa_x[k] * hsg0_305[k]
                   + f_8 * hsf_205[k]
                   - f_6 * pc_x[k] * hsg1_305[k];

        t_306[k] = f_7 * hsf_206[k]
                   + f_3 * pc_x[k] * isf_206[k];

        t_307[k] = f_7 * hsf_207[k]
                   + f_3 * pc_x[k] * isf_207[k];

        t_308[k] = f_3 * pc_y[k] * isf_205[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, t_313, pa_x, pc_x, pc_y, hsg0_310, \
                         hsg0_311, hsg0_312, hsf_209, hsg1_310, hsg1_311, hsg1_312, \
                         isf_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_7 * hsf_209[k]
                   + f_3 * pc_x[k] * isf_209[k];

        t_310[k] = pa_x[k] * hsg0_310[k]
                   - f_6 * pc_x[k] * hsg1_310[k];

        t_311[k] = pa_x[k] * hsg0_311[k]
                   - f_6 * pc_x[k] * hsg1_311[k];

        t_312[k] = pa_x[k] * hsg0_312[k]
                   - f_6 * pc_x[k] * hsg1_312[k];

        t_313[k] = f_3 * pc_y[k] * isf_209[k];
    }

#pragma omp simd aligned(t_314, t_315, t_316, t_317, pa_x, pc_x, pc_z, hsg0_314, hsg1_314, \
                         isd0_126, isd0_127, isd1_126, isd1_127, isf_210, \
                         isf_211 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_314[k] = pa_x[k] * hsg0_314[k]
                   - f_6 * pc_x[k] * hsg1_314[k];

        t_315[k] = f_1 * isd0_126[k]
                   - f_2 * isd1_126[k]
                   + f_3 * pc_x[k] * isf_210[k];

        t_316[k] = f_10 * isd0_127[k]
                   - f_11 * isd1_127[k]
                   + f_3 * pc_x[k] * isf_211[k];

        t_317[k] = f_3 * pc_z[k] * isf_210[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, t_321, t_322, pc_x, pc_z, isd0_129, isd0_131, \
                         isd1_129, isd1_131, isf_211, isf_213, isf_215, isf_216, \
                         isf_217 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = f_4 * isd0_129[k]
                   - f_5 * isd1_129[k]
                   + f_3 * pc_x[k] * isf_213[k];

        t_319[k] = f_3 * pc_z[k] * isf_211[k];

        t_320[k] = f_4 * isd0_131[k]
                   - f_5 * isd1_131[k]
                   + f_3 * pc_x[k] * isf_215[k];

        t_321[k] = f_3 * pc_x[k] * isf_216[k];

        t_322[k] = f_3 * pc_x[k] * isf_217[k];
    }

#pragma omp simd aligned(t_323, t_324, t_325, t_326, t_327, pc_x, pc_y, pc_z, hsf_156, \
                         isd0_129, isd1_129, isf_216, isf_217, isf_218, \
                         isf_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_323[k] = f_3 * pc_x[k] * isf_218[k];

        t_324[k] = f_3 * pc_x[k] * isf_219[k];

        t_325[k] = f_0 * hsf_156[k]
                   + f_1 * isd0_129[k]
                   - f_2 * isd1_129[k]
                   + f_3 * pc_y[k] * isf_216[k];

        t_326[k] = f_3 * pc_z[k] * isf_216[k];

        t_327[k] = f_4 * isd0_129[k]
                   - f_5 * isd1_129[k]
                   + f_3 * pc_z[k] * isf_217[k];
    }

#pragma omp simd aligned(t_328, t_329, t_330, t_331, pa_z, pc_y, pc_z, hsg0_225, hsg0_226, \
                         hsf_159, hsg1_225, hsg1_226, isd0_131, isd1_131, \
                         isf_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_328[k] = f_0 * hsf_159[k]
                   + f_3 * pc_y[k] * isf_219[k];

        t_329[k] = f_1 * isd0_131[k]
                   - f_2 * isd1_131[k]
                   + f_3 * pc_z[k] * isf_219[k];

        t_330[k] = pa_z[k] * hsg0_225[k]
                   - f_6 * pc_z[k] * hsg1_225[k];

        t_331[k] = pa_z[k] * hsg0_226[k]
                   - f_6 * pc_z[k] * hsg1_226[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, pa_z, pc_x, pc_z, hsg0_228, hsg1_228, isd0_134, \
                         isd0_136, isd1_134, isd1_136, isf_222, \
                         isf_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = f_10 * isd0_134[k]
                   - f_11 * isd1_134[k]
                   + f_3 * pc_x[k] * isf_222[k];

        t_333[k] = pa_z[k] * hsg0_228[k]
                   - f_6 * pc_z[k] * hsg1_228[k];

        t_334[k] = f_4 * isd0_136[k]
                   - f_5 * isd1_136[k]
                   + f_3 * pc_x[k] * isf_224[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, pc_x, isd0_137, isd1_137, isf_225, \
                         isf_226, isf_227, isf_228, isf_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = f_4 * isd0_137[k]
                   - f_5 * isd1_137[k]
                   + f_3 * pc_x[k] * isf_225[k];

        t_336[k] = f_3 * pc_x[k] * isf_226[k];

        t_337[k] = f_3 * pc_x[k] * isf_227[k];

        t_338[k] = f_3 * pc_x[k] * isf_228[k];

        t_339[k] = f_3 * pc_x[k] * isf_229[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, pa_z, pc_y, pc_z, hsg0_235, hsg0_237, \
                         hsf_156, hsf_157, hsf_169, hsg1_235, hsg1_237, isf_226, \
                         isf_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = pa_z[k] * hsg0_235[k]
                   - f_6 * pc_z[k] * hsg1_235[k];

        t_341[k] = f_7 * hsf_156[k]
                   + f_3 * pc_z[k] * isf_226[k];

        t_342[k] = pa_z[k] * hsg0_237[k]
                   + f_8 * hsf_157[k]
                   - f_6 * pc_z[k] * hsg1_237[k];

        t_343[k] = f_9 * hsf_169[k]
                   + f_3 * pc_y[k] * isf_229[k];
    }

#pragma omp simd aligned(t_344, t_345, t_346, pc_x, pc_z, hsf_159, isd0_137, isd0_138, \
                         isd0_139, isd1_137, isd1_138, isd1_139, isf_229, isf_230, \
                         isf_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_344[k] = f_7 * hsf_159[k]
                   + f_1 * isd0_137[k]
                   - f_2 * isd1_137[k]
                   + f_3 * pc_z[k] * isf_229[k];

        t_345[k] = f_1 * isd0_138[k]
                   - f_2 * isd1_138[k]
                   + f_3 * pc_x[k] * isf_230[k];

        t_346[k] = f_10 * isd0_139[k]
                   - f_11 * isd1_139[k]
                   + f_3 * pc_x[k] * isf_231[k];
    }

#pragma omp simd aligned(t_347, t_348, t_349, pc_x, isd0_140, isd0_141, isd0_142, isd1_140, \
                         isd1_141, isd1_142, isf_232, isf_233, \
                         isf_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_347[k] = f_10 * isd0_140[k]
                   - f_11 * isd1_140[k]
                   + f_3 * pc_x[k] * isf_232[k];

        t_348[k] = f_4 * isd0_141[k]
                   - f_5 * isd1_141[k]
                   + f_3 * pc_x[k] * isf_233[k];

        t_349[k] = f_4 * isd0_142[k]
                   - f_5 * isd1_142[k]
                   + f_3 * pc_x[k] * isf_234[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, pc_x, isd0_143, isd1_143, isf_235, \
                         isf_236, isf_237, isf_238, isf_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = f_4 * isd0_143[k]
                   - f_5 * isd1_143[k]
                   + f_3 * pc_x[k] * isf_235[k];

        t_351[k] = f_3 * pc_x[k] * isf_236[k];

        t_352[k] = f_3 * pc_x[k] * isf_237[k];

        t_353[k] = f_3 * pc_x[k] * isf_238[k];

        t_354[k] = f_3 * pc_x[k] * isf_239[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, pc_y, pc_z, hsf_166, hsf_176, hsf_178, isd0_141, \
                         isd0_143, isd1_141, isd1_143, isf_236, \
                         isf_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = f_12 * hsf_176[k]
                   + f_1 * isd0_141[k]
                   - f_2 * isd1_141[k]
                   + f_3 * pc_y[k] * isf_236[k];

        t_356[k] = f_8 * hsf_166[k]
                   + f_3 * pc_z[k] * isf_236[k];

        t_357[k] = f_12 * hsf_178[k]
                   + f_4 * isd0_143[k]
                   - f_5 * isd1_143[k]
                   + f_3 * pc_y[k] * isf_238[k];
    }

#pragma omp simd aligned(t_358, t_359, t_360, pc_x, pc_y, pc_z, hsf_169, hsf_179, isd0_143, \
                         isd0_144, isd1_143, isd1_144, isf_239, \
                         isf_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_358[k] = f_12 * hsf_179[k]
                   + f_3 * pc_y[k] * isf_239[k];

        t_359[k] = f_8 * hsf_169[k]
                   + f_1 * isd0_143[k]
                   - f_2 * isd1_143[k]
                   + f_3 * pc_z[k] * isf_239[k];

        t_360[k] = f_1 * isd0_144[k]
                   - f_2 * isd1_144[k]
                   + f_3 * pc_x[k] * isf_240[k];
    }

#pragma omp simd aligned(t_361, t_362, t_363, pc_x, isd0_145, isd0_146, isd0_147, isd1_145, \
                         isd1_146, isd1_147, isf_241, isf_242, \
                         isf_243 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_361[k] = f_10 * isd0_145[k]
                   - f_11 * isd1_145[k]
                   + f_3 * pc_x[k] * isf_241[k];

        t_362[k] = f_10 * isd0_146[k]
                   - f_11 * isd1_146[k]
                   + f_3 * pc_x[k] * isf_242[k];

        t_363[k] = f_4 * isd0_147[k]
                   - f_5 * isd1_147[k]
                   + f_3 * pc_x[k] * isf_243[k];
    }

#pragma omp simd aligned(t_364, t_365, t_366, t_367, t_368, pc_x, isd0_148, isd0_149, \
                         isd1_148, isd1_149, isf_244, isf_245, isf_246, isf_247, \
                         isf_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_364[k] = f_4 * isd0_148[k]
                   - f_5 * isd1_148[k]
                   + f_3 * pc_x[k] * isf_244[k];

        t_365[k] = f_4 * isd0_149[k]
                   - f_5 * isd1_149[k]
                   + f_3 * pc_x[k] * isf_245[k];

        t_366[k] = f_3 * pc_x[k] * isf_246[k];

        t_367[k] = f_3 * pc_x[k] * isf_247[k];

        t_368[k] = f_3 * pc_x[k] * isf_248[k];
    }

#pragma omp simd aligned(t_369, t_370, t_371, pc_x, pc_y, pc_z, hsf_176, hsf_186, isd0_147, \
                         isd1_147, isf_246, isf_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_369[k] = f_3 * pc_x[k] * isf_249[k];

        t_370[k] = f_13 * hsf_186[k]
                   + f_1 * isd0_147[k]
                   - f_2 * isd1_147[k]
                   + f_3 * pc_y[k] * isf_246[k];

        t_371[k] = f_13 * hsf_176[k]
                   + f_3 * pc_z[k] * isf_246[k];
    }

#pragma omp simd aligned(t_372, t_373, t_374, pc_y, pc_z, hsf_179, hsf_188, hsf_189, isd0_149, \
                         isd1_149, isf_248, isf_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_372[k] = f_13 * hsf_188[k]
                   + f_4 * isd0_149[k]
                   - f_5 * isd1_149[k]
                   + f_3 * pc_y[k] * isf_248[k];

        t_373[k] = f_13 * hsf_189[k]
                   + f_3 * pc_y[k] * isf_249[k];

        t_374[k] = f_13 * hsf_179[k]
                   + f_1 * isd0_149[k]
                   - f_2 * isd1_149[k]
                   + f_3 * pc_z[k] * isf_249[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, pc_x, isd0_150, isd0_151, isd0_152, isd1_150, \
                         isd1_151, isd1_152, isf_250, isf_251, \
                         isf_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = f_1 * isd0_150[k]
                   - f_2 * isd1_150[k]
                   + f_3 * pc_x[k] * isf_250[k];

        t_376[k] = f_10 * isd0_151[k]
                   - f_11 * isd1_151[k]
                   + f_3 * pc_x[k] * isf_251[k];

        t_377[k] = f_10 * isd0_152[k]
                   - f_11 * isd1_152[k]
                   + f_3 * pc_x[k] * isf_252[k];
    }
}

static auto
compute_prim_isg_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t hsg0,
                                                          const size_t hsf, const size_t hsg1,
                                                          const size_t isd0, const size_t isd1,
                                                          const size_t isf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 2.5 / q;
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 2.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *hsg0_300 = buffer.data(hsg0 + 300);
    const auto *hsg0_302 = buffer.data(hsg0 + 302);
    const auto *hsg0_305 = buffer.data(hsg0 + 305);
    const auto *hsg0_310 = buffer.data(hsg0 + 310);
    const auto *hsg0_312 = buffer.data(hsg0 + 312);
    const auto *hsg0_314 = buffer.data(hsg0 + 314);

    const auto *hsf_186 = buffer.data(hsf + 186);
    const auto *hsf_189 = buffer.data(hsf + 189);
    const auto *hsf_196 = buffer.data(hsf + 196);
    const auto *hsf_198 = buffer.data(hsf + 198);
    const auto *hsf_199 = buffer.data(hsf + 199);
    const auto *hsf_206 = buffer.data(hsf + 206);
    const auto *hsf_208 = buffer.data(hsf + 208);
    const auto *hsf_209 = buffer.data(hsf + 209);

    const auto *hsg1_300 = buffer.data(hsg1 + 300);
    const auto *hsg1_302 = buffer.data(hsg1 + 302);
    const auto *hsg1_305 = buffer.data(hsg1 + 305);
    const auto *hsg1_310 = buffer.data(hsg1 + 310);
    const auto *hsg1_312 = buffer.data(hsg1 + 312);
    const auto *hsg1_314 = buffer.data(hsg1 + 314);

    const auto *isd0_153 = buffer.data(isd0 + 153);
    const auto *isd0_154 = buffer.data(isd0 + 154);
    const auto *isd0_155 = buffer.data(isd0 + 155);
    const auto *isd0_157 = buffer.data(isd0 + 157);
    const auto *isd0_159 = buffer.data(isd0 + 159);
    const auto *isd0_160 = buffer.data(isd0 + 160);
    const auto *isd0_162 = buffer.data(isd0 + 162);
    const auto *isd0_164 = buffer.data(isd0 + 164);
    const auto *isd0_165 = buffer.data(isd0 + 165);
    const auto *isd0_166 = buffer.data(isd0 + 166);
    const auto *isd0_167 = buffer.data(isd0 + 167);

    const auto *isd1_153 = buffer.data(isd1 + 153);
    const auto *isd1_154 = buffer.data(isd1 + 154);
    const auto *isd1_155 = buffer.data(isd1 + 155);
    const auto *isd1_157 = buffer.data(isd1 + 157);
    const auto *isd1_159 = buffer.data(isd1 + 159);
    const auto *isd1_160 = buffer.data(isd1 + 160);
    const auto *isd1_162 = buffer.data(isd1 + 162);
    const auto *isd1_164 = buffer.data(isd1 + 164);
    const auto *isd1_165 = buffer.data(isd1 + 165);
    const auto *isd1_166 = buffer.data(isd1 + 166);
    const auto *isd1_167 = buffer.data(isd1 + 167);

    const auto *isf_253 = buffer.data(isf + 253);
    const auto *isf_254 = buffer.data(isf + 254);
    const auto *isf_255 = buffer.data(isf + 255);
    const auto *isf_256 = buffer.data(isf + 256);
    const auto *isf_257 = buffer.data(isf + 257);
    const auto *isf_258 = buffer.data(isf + 258);
    const auto *isf_259 = buffer.data(isf + 259);
    const auto *isf_261 = buffer.data(isf + 261);
    const auto *isf_263 = buffer.data(isf + 263);
    const auto *isf_264 = buffer.data(isf + 264);
    const auto *isf_266 = buffer.data(isf + 266);
    const auto *isf_267 = buffer.data(isf + 267);
    const auto *isf_268 = buffer.data(isf + 268);
    const auto *isf_269 = buffer.data(isf + 269);
    const auto *isf_270 = buffer.data(isf + 270);
    const auto *isf_272 = buffer.data(isf + 272);
    const auto *isf_273 = buffer.data(isf + 273);
    const auto *isf_275 = buffer.data(isf + 275);
    const auto *isf_276 = buffer.data(isf + 276);
    const auto *isf_277 = buffer.data(isf + 277);
    const auto *isf_278 = buffer.data(isf + 278);
    const auto *isf_279 = buffer.data(isf + 279);

#pragma omp simd aligned(t_378, t_379, t_380, t_381, pc_x, isd0_153, isd0_154, isd0_155, \
                         isd1_153, isd1_154, isd1_155, isf_253, isf_254, isf_255, \
                         isf_256 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_378[k] = f_4 * isd0_153[k]
                   - f_5 * isd1_153[k]
                   + f_3 * pc_x[k] * isf_253[k];

        t_379[k] = f_4 * isd0_154[k]
                   - f_5 * isd1_154[k]
                   + f_3 * pc_x[k] * isf_254[k];

        t_380[k] = f_4 * isd0_155[k]
                   - f_5 * isd1_155[k]
                   + f_3 * pc_x[k] * isf_255[k];

        t_381[k] = f_3 * pc_x[k] * isf_256[k];
    }

#pragma omp simd aligned(t_382, t_383, t_384, t_385, t_386, pc_x, pc_y, pc_z, hsf_186, \
                         hsf_196, isd0_153, isd1_153, isf_256, isf_257, isf_258, \
                         isf_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_382[k] = f_3 * pc_x[k] * isf_257[k];

        t_383[k] = f_3 * pc_x[k] * isf_258[k];

        t_384[k] = f_3 * pc_x[k] * isf_259[k];

        t_385[k] = f_8 * hsf_196[k]
                   + f_1 * isd0_153[k]
                   - f_2 * isd1_153[k]
                   + f_3 * pc_y[k] * isf_256[k];

        t_386[k] = f_12 * hsf_186[k]
                   + f_3 * pc_z[k] * isf_256[k];
    }

#pragma omp simd aligned(t_387, t_388, t_389, t_390, pa_y, pc_y, pc_z, hsg0_300, hsf_189, \
                         hsf_198, hsf_199, hsg1_300, isd0_155, isd1_155, isf_258, \
                         isf_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_387[k] = f_8 * hsf_198[k]
                   + f_4 * isd0_155[k]
                   - f_5 * isd1_155[k]
                   + f_3 * pc_y[k] * isf_258[k];

        t_388[k] = f_8 * hsf_199[k]
                   + f_3 * pc_y[k] * isf_259[k];

        t_389[k] = f_12 * hsf_189[k]
                   + f_1 * isd0_155[k]
                   - f_2 * isd1_155[k]
                   + f_3 * pc_z[k] * isf_259[k];

        t_390[k] = pa_y[k] * hsg0_300[k]
                   - f_6 * pc_y[k] * hsg1_300[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, pa_y, pc_x, pc_y, hsg0_302, hsg1_302, isd0_157, \
                         isd0_159, isd1_157, isd1_159, isf_261, \
                         isf_263 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = f_10 * isd0_157[k]
                   - f_11 * isd1_157[k]
                   + f_3 * pc_x[k] * isf_261[k];

        t_392[k] = pa_y[k] * hsg0_302[k]
                   - f_6 * pc_y[k] * hsg1_302[k];

        t_393[k] = f_4 * isd0_159[k]
                   - f_5 * isd1_159[k]
                   + f_3 * pc_x[k] * isf_263[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, t_397, t_398, pa_y, pc_x, pc_y, hsg0_305, \
                         hsg1_305, isd0_160, isd1_160, isf_264, isf_266, isf_267, \
                         isf_268 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = f_4 * isd0_160[k]
                   - f_5 * isd1_160[k]
                   + f_3 * pc_x[k] * isf_264[k];

        t_395[k] = pa_y[k] * hsg0_305[k]
                   - f_6 * pc_y[k] * hsg1_305[k];

        t_396[k] = f_3 * pc_x[k] * isf_266[k];

        t_397[k] = f_3 * pc_x[k] * isf_267[k];

        t_398[k] = f_3 * pc_x[k] * isf_268[k];
    }

#pragma omp simd aligned(t_399, t_400, t_401, pa_y, pc_x, pc_y, pc_z, hsg0_310, hsf_196, \
                         hsf_206, hsg1_310, isf_266, isf_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_399[k] = f_3 * pc_x[k] * isf_269[k];

        t_400[k] = pa_y[k] * hsg0_310[k]
                   + f_12 * hsf_206[k]
                   - f_6 * pc_y[k] * hsg1_310[k];

        t_401[k] = f_9 * hsf_196[k]
                   + f_3 * pc_z[k] * isf_266[k];
    }

#pragma omp simd aligned(t_402, t_403, t_404, pa_y, pc_y, hsg0_312, hsg0_314, hsf_208, \
                         hsf_209, hsg1_312, hsg1_314, isf_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_402[k] = pa_y[k] * hsg0_312[k]
                   + f_8 * hsf_208[k]
                   - f_6 * pc_y[k] * hsg1_312[k];

        t_403[k] = f_7 * hsf_209[k]
                   + f_3 * pc_y[k] * isf_269[k];

        t_404[k] = pa_y[k] * hsg0_314[k]
                   - f_6 * pc_y[k] * hsg1_314[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, pc_x, pc_y, isd0_162, isd0_164, \
                         isd0_165, isd1_162, isd1_164, isd1_165, isf_270, isf_272, \
                         isf_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = f_1 * isd0_162[k]
                   - f_2 * isd1_162[k]
                   + f_3 * pc_x[k] * isf_270[k];

        t_406[k] = f_3 * pc_y[k] * isf_270[k];

        t_407[k] = f_10 * isd0_164[k]
                   - f_11 * isd1_164[k]
                   + f_3 * pc_x[k] * isf_272[k];

        t_408[k] = f_4 * isd0_165[k]
                   - f_5 * isd1_165[k]
                   + f_3 * pc_x[k] * isf_273[k];

        t_409[k] = f_3 * pc_y[k] * isf_272[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, t_414, pc_x, isd0_167, isd1_167, isf_275, \
                         isf_276, isf_277, isf_278, isf_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = f_4 * isd0_167[k]
                   - f_5 * isd1_167[k]
                   + f_3 * pc_x[k] * isf_275[k];

        t_411[k] = f_3 * pc_x[k] * isf_276[k];

        t_412[k] = f_3 * pc_x[k] * isf_277[k];

        t_413[k] = f_3 * pc_x[k] * isf_278[k];

        t_414[k] = f_3 * pc_x[k] * isf_279[k];
    }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, pc_y, isd0_165, isd0_166, isd0_167, \
                         isd1_165, isd1_166, isd1_167, isf_276, isf_277, isf_278, \
                         isf_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_415[k] = f_1 * isd0_165[k]
                   - f_2 * isd1_165[k]
                   + f_3 * pc_y[k] * isf_276[k];

        t_416[k] = f_10 * isd0_166[k]
                   - f_11 * isd1_166[k]
                   + f_3 * pc_y[k] * isf_277[k];

        t_417[k] = f_4 * isd0_167[k]
                   - f_5 * isd1_167[k]
                   + f_3 * pc_y[k] * isf_278[k];

        t_418[k] = f_3 * pc_y[k] * isf_279[k];
    }

#pragma omp simd aligned(t_419, pc_z, hsf_209, isd0_167, isd1_167, \
                         isf_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_419[k] = f_0 * hsf_209[k]
                   + f_1 * isd0_167[k]
                   - f_2 * isd1_167[k]
                   + f_3 * pc_z[k] * isf_279[k];
    }
}

auto
compute_prim_isg_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t hsg0, const size_t hsf,
                                                   const size_t hsg1, const size_t isd0,
                                                   const size_t isd1, const size_t isf,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_isg_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, hsg0, hsf,
                                                              hsg1, isd0, isd1, isf, ncols,
                                                              gamma, p, q);

    compute_prim_isg_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, hsg0, hsf,
                                                              hsg1, isd0, isd1, isf, ncols,
                                                              gamma, p, q);

    compute_prim_isg_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, hsg0, hsf,
                                                              hsg1, isd0, isd1, isf, ncols,
                                                              gamma, p, q);

    compute_prim_isg_three_center_electron_repulsion_0_piece3(buffer, target, pa, pc, hsg0, hsf,
                                                              hsg1, isd0, isd1, isf, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
