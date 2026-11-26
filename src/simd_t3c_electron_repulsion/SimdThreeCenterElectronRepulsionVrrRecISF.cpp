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


#include "SimdThreeCenterElectronRepulsionVrrRecISF.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_isf_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t hsf0,
                                                          const size_t hsd, const size_t hsf1,
                                                          const size_t isp0, const size_t isp1,
                                                          const size_t isd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 2.5 / q;
    const auto f_7 = 0.5 / gamma;
    const auto f_8 = 0.5 * p / (gamma * q);
    const auto f_9 = 2.0 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;

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

    const auto *hsf0_0 = buffer.data(hsf0 + 0);
    const auto *hsf0_6 = buffer.data(hsf0 + 6);
    const auto *hsf0_9 = buffer.data(hsf0 + 9);
    const auto *hsf0_16 = buffer.data(hsf0 + 16);
    const auto *hsf0_20 = buffer.data(hsf0 + 20);
    const auto *hsf0_29 = buffer.data(hsf0 + 29);
    const auto *hsf0_30 = buffer.data(hsf0 + 30);
    const auto *hsf0_36 = buffer.data(hsf0 + 36);
    const auto *hsf0_50 = buffer.data(hsf0 + 50);
    const auto *hsf0_59 = buffer.data(hsf0 + 59);
    const auto *hsf0_60 = buffer.data(hsf0 + 60);
    const auto *hsf0_66 = buffer.data(hsf0 + 66);
    const auto *hsf0_90 = buffer.data(hsf0 + 90);

    const auto *hsd_0 = buffer.data(hsd + 0);
    const auto *hsd_3 = buffer.data(hsd + 3);
    const auto *hsd_5 = buffer.data(hsd + 5);
    const auto *hsd_6 = buffer.data(hsd + 6);
    const auto *hsd_9 = buffer.data(hsd + 9);
    const auto *hsd_11 = buffer.data(hsd + 11);
    const auto *hsd_12 = buffer.data(hsd + 12);
    const auto *hsd_15 = buffer.data(hsd + 15);
    const auto *hsd_17 = buffer.data(hsd + 17);
    const auto *hsd_18 = buffer.data(hsd + 18);
    const auto *hsd_21 = buffer.data(hsd + 21);
    const auto *hsd_23 = buffer.data(hsd + 23);
    const auto *hsd_24 = buffer.data(hsd + 24);
    const auto *hsd_27 = buffer.data(hsd + 27);
    const auto *hsd_28 = buffer.data(hsd + 28);
    const auto *hsd_29 = buffer.data(hsd + 29);
    const auto *hsd_30 = buffer.data(hsd + 30);
    const auto *hsd_33 = buffer.data(hsd + 33);
    const auto *hsd_35 = buffer.data(hsd + 35);
    const auto *hsd_36 = buffer.data(hsd + 36);
    const auto *hsd_39 = buffer.data(hsd + 39);
    const auto *hsd_41 = buffer.data(hsd + 41);
    const auto *hsd_42 = buffer.data(hsd + 42);
    const auto *hsd_45 = buffer.data(hsd + 45);
    const auto *hsd_46 = buffer.data(hsd + 46);
    const auto *hsd_47 = buffer.data(hsd + 47);
    const auto *hsd_48 = buffer.data(hsd + 48);
    const auto *hsd_51 = buffer.data(hsd + 51);
    const auto *hsd_52 = buffer.data(hsd + 52);
    const auto *hsd_53 = buffer.data(hsd + 53);
    const auto *hsd_54 = buffer.data(hsd + 54);
    const auto *hsd_57 = buffer.data(hsd + 57);
    const auto *hsd_59 = buffer.data(hsd + 59);
    const auto *hsd_60 = buffer.data(hsd + 60);
    const auto *hsd_63 = buffer.data(hsd + 63);
    const auto *hsd_65 = buffer.data(hsd + 65);
    const auto *hsd_69 = buffer.data(hsd + 69);
    const auto *hsd_70 = buffer.data(hsd + 70);
    const auto *hsd_71 = buffer.data(hsd + 71);
    const auto *hsd_72 = buffer.data(hsd + 72);
    const auto *hsd_75 = buffer.data(hsd + 75);
    const auto *hsd_76 = buffer.data(hsd + 76);
    const auto *hsd_77 = buffer.data(hsd + 77);

    const auto *hsf1_0 = buffer.data(hsf1 + 0);
    const auto *hsf1_6 = buffer.data(hsf1 + 6);
    const auto *hsf1_9 = buffer.data(hsf1 + 9);
    const auto *hsf1_16 = buffer.data(hsf1 + 16);
    const auto *hsf1_20 = buffer.data(hsf1 + 20);
    const auto *hsf1_29 = buffer.data(hsf1 + 29);
    const auto *hsf1_30 = buffer.data(hsf1 + 30);
    const auto *hsf1_36 = buffer.data(hsf1 + 36);
    const auto *hsf1_50 = buffer.data(hsf1 + 50);
    const auto *hsf1_59 = buffer.data(hsf1 + 59);
    const auto *hsf1_60 = buffer.data(hsf1 + 60);
    const auto *hsf1_66 = buffer.data(hsf1 + 66);
    const auto *hsf1_90 = buffer.data(hsf1 + 90);

    const auto *isp0_0 = buffer.data(isp0 + 0);
    const auto *isp0_1 = buffer.data(isp0 + 1);
    const auto *isp0_2 = buffer.data(isp0 + 2);
    const auto *isp0_4 = buffer.data(isp0 + 4);
    const auto *isp0_8 = buffer.data(isp0 + 8);
    const auto *isp0_9 = buffer.data(isp0 + 9);
    const auto *isp0_10 = buffer.data(isp0 + 10);
    const auto *isp0_11 = buffer.data(isp0 + 11);
    const auto *isp0_15 = buffer.data(isp0 + 15);
    const auto *isp0_16 = buffer.data(isp0 + 16);
    const auto *isp0_17 = buffer.data(isp0 + 17);
    const auto *isp0_18 = buffer.data(isp0 + 18);
    const auto *isp0_19 = buffer.data(isp0 + 19);
    const auto *isp0_20 = buffer.data(isp0 + 20);
    const auto *isp0_23 = buffer.data(isp0 + 23);
    const auto *isp0_25 = buffer.data(isp0 + 25);
    const auto *isp0_27 = buffer.data(isp0 + 27);
    const auto *isp0_28 = buffer.data(isp0 + 28);
    const auto *isp0_29 = buffer.data(isp0 + 29);
    const auto *isp0_30 = buffer.data(isp0 + 30);
    const auto *isp0_31 = buffer.data(isp0 + 31);
    const auto *isp0_32 = buffer.data(isp0 + 32);
    const auto *isp0_35 = buffer.data(isp0 + 35);
    const auto *isp0_36 = buffer.data(isp0 + 36);
    const auto *isp0_37 = buffer.data(isp0 + 37);
    const auto *isp0_38 = buffer.data(isp0 + 38);

    const auto *isp1_0 = buffer.data(isp1 + 0);
    const auto *isp1_1 = buffer.data(isp1 + 1);
    const auto *isp1_2 = buffer.data(isp1 + 2);
    const auto *isp1_4 = buffer.data(isp1 + 4);
    const auto *isp1_8 = buffer.data(isp1 + 8);
    const auto *isp1_9 = buffer.data(isp1 + 9);
    const auto *isp1_10 = buffer.data(isp1 + 10);
    const auto *isp1_11 = buffer.data(isp1 + 11);
    const auto *isp1_15 = buffer.data(isp1 + 15);
    const auto *isp1_16 = buffer.data(isp1 + 16);
    const auto *isp1_17 = buffer.data(isp1 + 17);
    const auto *isp1_18 = buffer.data(isp1 + 18);
    const auto *isp1_19 = buffer.data(isp1 + 19);
    const auto *isp1_20 = buffer.data(isp1 + 20);
    const auto *isp1_23 = buffer.data(isp1 + 23);
    const auto *isp1_25 = buffer.data(isp1 + 25);
    const auto *isp1_27 = buffer.data(isp1 + 27);
    const auto *isp1_28 = buffer.data(isp1 + 28);
    const auto *isp1_29 = buffer.data(isp1 + 29);
    const auto *isp1_30 = buffer.data(isp1 + 30);
    const auto *isp1_31 = buffer.data(isp1 + 31);
    const auto *isp1_32 = buffer.data(isp1 + 32);
    const auto *isp1_35 = buffer.data(isp1 + 35);
    const auto *isp1_36 = buffer.data(isp1 + 36);
    const auto *isp1_37 = buffer.data(isp1 + 37);
    const auto *isp1_38 = buffer.data(isp1 + 38);

    const auto *isd_0 = buffer.data(isd + 0);
    const auto *isd_2 = buffer.data(isd + 2);
    const auto *isd_3 = buffer.data(isd + 3);
    const auto *isd_5 = buffer.data(isd + 5);
    const auto *isd_6 = buffer.data(isd + 6);
    const auto *isd_7 = buffer.data(isd + 7);
    const auto *isd_9 = buffer.data(isd + 9);
    const auto *isd_11 = buffer.data(isd + 11);
    const auto *isd_12 = buffer.data(isd + 12);
    const auto *isd_14 = buffer.data(isd + 14);
    const auto *isd_15 = buffer.data(isd + 15);
    const auto *isd_16 = buffer.data(isd + 16);
    const auto *isd_17 = buffer.data(isd + 17);
    const auto *isd_18 = buffer.data(isd + 18);
    const auto *isd_19 = buffer.data(isd + 19);
    const auto *isd_21 = buffer.data(isd + 21);
    const auto *isd_23 = buffer.data(isd + 23);
    const auto *isd_24 = buffer.data(isd + 24);
    const auto *isd_27 = buffer.data(isd + 27);
    const auto *isd_28 = buffer.data(isd + 28);
    const auto *isd_29 = buffer.data(isd + 29);
    const auto *isd_30 = buffer.data(isd + 30);
    const auto *isd_32 = buffer.data(isd + 32);
    const auto *isd_33 = buffer.data(isd + 33);
    const auto *isd_34 = buffer.data(isd + 34);
    const auto *isd_35 = buffer.data(isd + 35);
    const auto *isd_36 = buffer.data(isd + 36);
    const auto *isd_37 = buffer.data(isd + 37);
    const auto *isd_39 = buffer.data(isd + 39);
    const auto *isd_41 = buffer.data(isd + 41);
    const auto *isd_42 = buffer.data(isd + 42);
    const auto *isd_45 = buffer.data(isd + 45);
    const auto *isd_46 = buffer.data(isd + 46);
    const auto *isd_47 = buffer.data(isd + 47);
    const auto *isd_48 = buffer.data(isd + 48);
    const auto *isd_51 = buffer.data(isd + 51);
    const auto *isd_52 = buffer.data(isd + 52);
    const auto *isd_53 = buffer.data(isd + 53);
    const auto *isd_54 = buffer.data(isd + 54);
    const auto *isd_56 = buffer.data(isd + 56);
    const auto *isd_57 = buffer.data(isd + 57);
    const auto *isd_58 = buffer.data(isd + 58);
    const auto *isd_59 = buffer.data(isd + 59);
    const auto *isd_60 = buffer.data(isd + 60);
    const auto *isd_61 = buffer.data(isd + 61);
    const auto *isd_63 = buffer.data(isd + 63);
    const auto *isd_65 = buffer.data(isd + 65);
    const auto *isd_66 = buffer.data(isd + 66);
    const auto *isd_69 = buffer.data(isd + 69);
    const auto *isd_70 = buffer.data(isd + 70);
    const auto *isd_71 = buffer.data(isd + 71);
    const auto *isd_72 = buffer.data(isd + 72);
    const auto *isd_75 = buffer.data(isd + 75);
    const auto *isd_76 = buffer.data(isd + 76);
    const auto *isd_77 = buffer.data(isd + 77);
    const auto *isd_78 = buffer.data(isd + 78);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pc_x, pc_y, pc_z, hsd_0, hsd_3, isp0_0, \
                         isp1_0, isd_0, isd_2, isd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hsd_0[k]
                 + f_1 * isp0_0[k]
                 - f_2 * isp1_0[k]
                 + f_3 * pc_x[k] * isd_0[k];

        t_1[k] = f_3 * pc_y[k] * isd_0[k];

        t_2[k] = f_3 * pc_z[k] * isd_0[k];

        t_3[k] = f_0 * hsd_3[k]
                 + f_3 * pc_x[k] * isd_3[k];

        t_4[k] = f_3 * pc_y[k] * isd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pc_x, pc_y, pc_z, hsd_5, isp0_1, isp0_2, \
                         isp1_1, isp1_2, isd_3, isd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * hsd_5[k]
                 + f_3 * pc_x[k] * isd_5[k];

        t_6[k] = f_1 * isp0_1[k]
                 - f_2 * isp1_1[k]
                 + f_3 * pc_y[k] * isd_3[k];

        t_7[k] = f_3 * pc_z[k] * isd_3[k];

        t_8[k] = f_3 * pc_y[k] * isd_5[k];

        t_9[k] = f_1 * isp0_2[k]
                 - f_2 * isp1_2[k]
                 + f_3 * pc_z[k] * isd_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_y, pc_x, pc_y, pc_z, hsf0_0, hsd_0, \
                         hsd_9, hsf1_0, isd_6, isd_7, isd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_y[k] * hsf0_0[k]
                  - f_4 * pc_y[k] * hsf1_0[k];

        t_11[k] = f_5 * hsd_0[k]
                  + f_3 * pc_y[k] * isd_6[k];

        t_12[k] = f_3 * pc_z[k] * isd_6[k];

        t_13[k] = f_6 * hsd_9[k]
                  + f_3 * pc_x[k] * isd_9[k];

        t_14[k] = f_3 * pc_z[k] * isd_7[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pc_x, pc_y, pc_z, hsd_3, hsd_5, hsd_11, \
                         isp0_4, isp1_4, isd_9, isd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_6 * hsd_11[k]
                  + f_3 * pc_x[k] * isd_11[k];

        t_16[k] = f_5 * hsd_3[k]
                  + f_1 * isp0_4[k]
                  - f_2 * isp1_4[k]
                  + f_3 * pc_y[k] * isd_9[k];

        t_17[k] = f_3 * pc_z[k] * isd_9[k];

        t_18[k] = f_5 * hsd_5[k]
                  + f_3 * pc_y[k] * isd_11[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_y, pa_z, pc_y, pc_z, hsf0_0, hsf0_9, \
                         hsd_0, hsf1_0, hsf1_9, isd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pa_y[k] * hsf0_9[k]
                  - f_4 * pc_y[k] * hsf1_9[k];

        t_20[k] = pa_z[k] * hsf0_0[k]
                  - f_4 * pc_z[k] * hsf1_0[k];

        t_21[k] = f_3 * pc_y[k] * isd_12[k];

        t_22[k] = f_5 * hsd_0[k]
                  + f_3 * pc_z[k] * isd_12[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_z, pc_x, pc_y, pc_z, hsf0_6, hsd_15, \
                         hsd_17, hsf1_6, isd_14, isd_15, isd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_6 * hsd_15[k]
                  + f_3 * pc_x[k] * isd_15[k];

        t_24[k] = f_3 * pc_y[k] * isd_14[k];

        t_25[k] = f_6 * hsd_17[k]
                  + f_3 * pc_x[k] * isd_17[k];

        t_26[k] = pa_z[k] * hsf0_6[k]
                  - f_4 * pc_z[k] * hsf1_6[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pc_x, pc_y, pc_z, hsd_5, hsd_18, isp0_8, \
                         isp0_9, isp1_8, isp1_9, isd_16, isd_17, \
                         isd_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_7 * isp0_8[k]
                  - f_8 * isp1_8[k]
                  + f_3 * pc_y[k] * isd_16[k];

        t_28[k] = f_3 * pc_y[k] * isd_17[k];

        t_29[k] = f_5 * hsd_5[k]
                  + f_1 * isp0_8[k]
                  - f_2 * isp1_8[k]
                  + f_3 * pc_z[k] * isd_17[k];

        t_30[k] = f_9 * hsd_18[k]
                  + f_1 * isp0_9[k]
                  - f_2 * isp1_9[k]
                  + f_3 * pc_x[k] * isd_18[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, pc_x, pc_y, pc_z, hsd_6, hsd_21, \
                         hsd_23, isd_18, isd_19, isd_21, isd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_10 * hsd_6[k]
                  + f_3 * pc_y[k] * isd_18[k];

        t_32[k] = f_3 * pc_z[k] * isd_18[k];

        t_33[k] = f_9 * hsd_21[k]
                  + f_3 * pc_x[k] * isd_21[k];

        t_34[k] = f_3 * pc_z[k] * isd_19[k];

        t_35[k] = f_9 * hsd_23[k]
                  + f_3 * pc_x[k] * isd_23[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pc_y, pc_z, hsd_9, hsd_11, isp0_10, isp0_11, \
                         isp1_10, isp1_11, isd_21, isd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_10 * hsd_9[k]
                  + f_1 * isp0_10[k]
                  - f_2 * isp1_10[k]
                  + f_3 * pc_y[k] * isd_21[k];

        t_37[k] = f_3 * pc_z[k] * isd_21[k];

        t_38[k] = f_10 * hsd_11[k]
                  + f_3 * pc_y[k] * isd_23[k];

        t_39[k] = f_1 * isp0_11[k]
                  - f_2 * isp1_11[k]
                  + f_3 * pc_z[k] * isd_23[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_y, pc_x, pc_y, pc_z, hsf0_20, hsd_6, \
                         hsd_12, hsd_27, hsf1_20, isd_24, isd_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pa_y[k] * hsf0_20[k]
                  - f_4 * pc_y[k] * hsf1_20[k];

        t_41[k] = f_5 * hsd_12[k]
                  + f_3 * pc_y[k] * isd_24[k];

        t_42[k] = f_5 * hsd_6[k]
                  + f_3 * pc_z[k] * isd_24[k];

        t_43[k] = f_9 * hsd_27[k]
                  + f_3 * pc_x[k] * isd_27[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_z, pc_x, pc_z, hsf0_16, hsd_9, hsd_28, \
                         hsd_29, hsf1_16, isd_27, isd_28, isd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_9 * hsd_28[k]
                  + f_3 * pc_x[k] * isd_28[k];

        t_45[k] = f_9 * hsd_29[k]
                  + f_3 * pc_x[k] * isd_29[k];

        t_46[k] = pa_z[k] * hsf0_16[k]
                  - f_4 * pc_z[k] * hsf1_16[k];

        t_47[k] = f_5 * hsd_9[k]
                  + f_3 * pc_z[k] * isd_27[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_y, pc_x, pc_y, hsf0_29, hsd_17, hsd_30, \
                         hsf1_29, isp0_15, isp1_15, isd_29, isd_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_5 * hsd_17[k]
                  + f_3 * pc_y[k] * isd_29[k];

        t_49[k] = pa_y[k] * hsf0_29[k]
                  - f_4 * pc_y[k] * hsf1_29[k];

        t_50[k] = f_9 * hsd_30[k]
                  + f_1 * isp0_15[k]
                  - f_2 * isp1_15[k]
                  + f_3 * pc_x[k] * isd_30[k];

        t_51[k] = f_3 * pc_y[k] * isd_30[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pc_x, pc_y, pc_z, hsd_12, hsd_33, hsd_35, \
                         isd_30, isd_32, isd_33, isd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_10 * hsd_12[k]
                  + f_3 * pc_z[k] * isd_30[k];

        t_53[k] = f_9 * hsd_33[k]
                  + f_3 * pc_x[k] * isd_33[k];

        t_54[k] = f_3 * pc_y[k] * isd_32[k];

        t_55[k] = f_9 * hsd_35[k]
                  + f_3 * pc_x[k] * isd_35[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pc_y, pc_z, hsd_17, isp0_16, isp0_17, \
                         isp1_16, isp1_17, isd_33, isd_34, isd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_1 * isp0_16[k]
                  - f_2 * isp1_16[k]
                  + f_3 * pc_y[k] * isd_33[k];

        t_57[k] = f_7 * isp0_17[k]
                  - f_8 * isp1_17[k]
                  + f_3 * pc_y[k] * isd_34[k];

        t_58[k] = f_3 * pc_y[k] * isd_35[k];

        t_59[k] = f_10 * hsd_17[k]
                  + f_1 * isp0_17[k]
                  - f_2 * isp1_17[k]
                  + f_3 * pc_z[k] * isd_35[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, pc_x, pc_y, pc_z, hsd_18, hsd_36, \
                         hsd_39, isp0_18, isp1_18, isd_36, isd_37, \
                         isd_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_11 * hsd_36[k]
                  + f_1 * isp0_18[k]
                  - f_2 * isp1_18[k]
                  + f_3 * pc_x[k] * isd_36[k];

        t_61[k] = f_11 * hsd_18[k]
                  + f_3 * pc_y[k] * isd_36[k];

        t_62[k] = f_3 * pc_z[k] * isd_36[k];

        t_63[k] = f_11 * hsd_39[k]
                  + f_3 * pc_x[k] * isd_39[k];

        t_64[k] = f_3 * pc_z[k] * isd_37[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pc_x, pc_y, pc_z, hsd_21, hsd_23, hsd_41, \
                         isp0_19, isp1_19, isd_39, isd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_11 * hsd_41[k]
                  + f_3 * pc_x[k] * isd_41[k];

        t_66[k] = f_11 * hsd_21[k]
                  + f_1 * isp0_19[k]
                  - f_2 * isp1_19[k]
                  + f_3 * pc_y[k] * isd_39[k];

        t_67[k] = f_3 * pc_z[k] * isd_39[k];

        t_68[k] = f_11 * hsd_23[k]
                  + f_3 * pc_y[k] * isd_41[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pa_z, pc_y, pc_z, hsf0_30, hsd_18, hsd_24, \
                         hsf1_30, isp0_20, isp1_20, isd_41, isd_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_1 * isp0_20[k]
                  - f_2 * isp1_20[k]
                  + f_3 * pc_z[k] * isd_41[k];

        t_70[k] = pa_z[k] * hsf0_30[k]
                  - f_4 * pc_z[k] * hsf1_30[k];

        t_71[k] = f_10 * hsd_24[k]
                  + f_3 * pc_y[k] * isd_42[k];

        t_72[k] = f_5 * hsd_18[k]
                  + f_3 * pc_z[k] * isd_42[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pa_z, pc_x, pc_z, hsf0_36, hsd_45, hsd_46, \
                         hsd_47, hsf1_36, isd_45, isd_46, isd_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_11 * hsd_45[k]
                  + f_3 * pc_x[k] * isd_45[k];

        t_74[k] = f_11 * hsd_46[k]
                  + f_3 * pc_x[k] * isd_46[k];

        t_75[k] = f_11 * hsd_47[k]
                  + f_3 * pc_x[k] * isd_47[k];

        t_76[k] = pa_z[k] * hsf0_36[k]
                  - f_4 * pc_z[k] * hsf1_36[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, pa_y, pc_y, pc_z, hsf0_50, hsd_21, hsd_23, \
                         hsd_29, hsf1_50, isp0_23, isp1_23, isd_45, \
                         isd_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_5 * hsd_21[k]
                  + f_3 * pc_z[k] * isd_45[k];

        t_78[k] = f_10 * hsd_29[k]
                  + f_3 * pc_y[k] * isd_47[k];

        t_79[k] = f_5 * hsd_23[k]
                  + f_1 * isp0_23[k]
                  - f_2 * isp1_23[k]
                  + f_3 * pc_z[k] * isd_47[k];

        t_80[k] = pa_y[k] * hsf0_50[k]
                  - f_4 * pc_y[k] * hsf1_50[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pc_x, pc_y, pc_z, hsd_24, hsd_30, hsd_51, \
                         hsd_52, isd_48, isd_51, isd_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_5 * hsd_30[k]
                  + f_3 * pc_y[k] * isd_48[k];

        t_82[k] = f_10 * hsd_24[k]
                  + f_3 * pc_z[k] * isd_48[k];

        t_83[k] = f_11 * hsd_51[k]
                  + f_3 * pc_x[k] * isd_51[k];

        t_84[k] = f_11 * hsd_52[k]
                  + f_3 * pc_x[k] * isd_52[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pc_x, pc_y, pc_z, hsd_27, hsd_33, hsd_35, \
                         hsd_53, isp0_25, isp1_25, isd_51, isd_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_11 * hsd_53[k]
                  + f_3 * pc_x[k] * isd_53[k];

        t_86[k] = f_5 * hsd_33[k]
                  + f_1 * isp0_25[k]
                  - f_2 * isp1_25[k]
                  + f_3 * pc_y[k] * isd_51[k];

        t_87[k] = f_10 * hsd_27[k]
                  + f_3 * pc_z[k] * isd_51[k];

        t_88[k] = f_5 * hsd_35[k]
                  + f_3 * pc_y[k] * isd_53[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, pa_y, pc_x, pc_y, pc_z, hsf0_59, hsd_30, \
                         hsd_54, hsf1_59, isp0_27, isp1_27, isd_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = pa_y[k] * hsf0_59[k]
                  - f_4 * pc_y[k] * hsf1_59[k];

        t_90[k] = f_11 * hsd_54[k]
                  + f_1 * isp0_27[k]
                  - f_2 * isp1_27[k]
                  + f_3 * pc_x[k] * isd_54[k];

        t_91[k] = f_3 * pc_y[k] * isd_54[k];

        t_92[k] = f_11 * hsd_30[k]
                  + f_3 * pc_z[k] * isd_54[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, pc_x, pc_y, hsd_57, hsd_59, isp0_28, isp1_28, \
                         isd_56, isd_57, isd_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_11 * hsd_57[k]
                  + f_3 * pc_x[k] * isd_57[k];

        t_94[k] = f_3 * pc_y[k] * isd_56[k];

        t_95[k] = f_11 * hsd_59[k]
                  + f_3 * pc_x[k] * isd_59[k];

        t_96[k] = f_1 * isp0_28[k]
                  - f_2 * isp1_28[k]
                  + f_3 * pc_y[k] * isd_57[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pc_x, pc_y, pc_z, hsd_35, hsd_60, isp0_29, \
                         isp0_30, isp1_29, isp1_30, isd_58, isd_59, \
                         isd_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_7 * isp0_29[k]
                  - f_8 * isp1_29[k]
                  + f_3 * pc_y[k] * isd_58[k];

        t_98[k] = f_3 * pc_y[k] * isd_59[k];

        t_99[k] = f_11 * hsd_35[k]
                  + f_1 * isp0_29[k]
                  - f_2 * isp1_29[k]
                  + f_3 * pc_z[k] * isd_59[k];

        t_100[k] = f_10 * hsd_60[k]
                   + f_1 * isp0_30[k]
                   - f_2 * isp1_30[k]
                   + f_3 * pc_x[k] * isd_60[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, t_105, pc_x, pc_y, pc_z, hsd_36, hsd_63, \
                         hsd_65, isd_60, isd_61, isd_63, isd_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_9 * hsd_36[k]
                   + f_3 * pc_y[k] * isd_60[k];

        t_102[k] = f_3 * pc_z[k] * isd_60[k];

        t_103[k] = f_10 * hsd_63[k]
                   + f_3 * pc_x[k] * isd_63[k];

        t_104[k] = f_3 * pc_z[k] * isd_61[k];

        t_105[k] = f_10 * hsd_65[k]
                   + f_3 * pc_x[k] * isd_65[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, pc_y, pc_z, hsd_39, hsd_41, isp0_31, \
                         isp0_32, isp1_31, isp1_32, isd_63, isd_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_9 * hsd_39[k]
                   + f_1 * isp0_31[k]
                   - f_2 * isp1_31[k]
                   + f_3 * pc_y[k] * isd_63[k];

        t_107[k] = f_3 * pc_z[k] * isd_63[k];

        t_108[k] = f_9 * hsd_41[k]
                   + f_3 * pc_y[k] * isd_65[k];

        t_109[k] = f_1 * isp0_32[k]
                   - f_2 * isp1_32[k]
                   + f_3 * pc_z[k] * isd_65[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, pa_z, pc_x, pc_y, pc_z, hsf0_60, hsd_36, \
                         hsd_42, hsd_69, hsf1_60, isd_66, isd_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = pa_z[k] * hsf0_60[k]
                   - f_4 * pc_z[k] * hsf1_60[k];

        t_111[k] = f_11 * hsd_42[k]
                   + f_3 * pc_y[k] * isd_66[k];

        t_112[k] = f_5 * hsd_36[k]
                   + f_3 * pc_z[k] * isd_66[k];

        t_113[k] = f_10 * hsd_69[k]
                   + f_3 * pc_x[k] * isd_69[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pa_z, pc_x, pc_z, hsf0_66, hsd_39, \
                         hsd_70, hsd_71, hsf1_66, isd_69, isd_70, \
                         isd_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_10 * hsd_70[k]
                   + f_3 * pc_x[k] * isd_70[k];

        t_115[k] = f_10 * hsd_71[k]
                   + f_3 * pc_x[k] * isd_71[k];

        t_116[k] = pa_z[k] * hsf0_66[k]
                   - f_4 * pc_z[k] * hsf1_66[k];

        t_117[k] = f_5 * hsd_39[k]
                   + f_3 * pc_z[k] * isd_69[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, pc_x, pc_y, pc_z, hsd_41, hsd_47, hsd_72, \
                         isp0_35, isp0_36, isp1_35, isp1_36, isd_71, \
                         isd_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_11 * hsd_47[k]
                   + f_3 * pc_y[k] * isd_71[k];

        t_119[k] = f_5 * hsd_41[k]
                   + f_1 * isp0_35[k]
                   - f_2 * isp1_35[k]
                   + f_3 * pc_z[k] * isd_71[k];

        t_120[k] = f_10 * hsd_72[k]
                   + f_1 * isp0_36[k]
                   - f_2 * isp1_36[k]
                   + f_3 * pc_x[k] * isd_72[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pc_x, pc_y, pc_z, hsd_42, hsd_48, hsd_75, \
                         hsd_76, isd_72, isd_75, isd_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_10 * hsd_48[k]
                   + f_3 * pc_y[k] * isd_72[k];

        t_122[k] = f_10 * hsd_42[k]
                   + f_3 * pc_z[k] * isd_72[k];

        t_123[k] = f_10 * hsd_75[k]
                   + f_3 * pc_x[k] * isd_75[k];

        t_124[k] = f_10 * hsd_76[k]
                   + f_3 * pc_x[k] * isd_76[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pc_x, pc_y, pc_z, hsd_45, hsd_51, hsd_53, \
                         hsd_77, isp0_37, isp1_37, isd_75, isd_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_10 * hsd_77[k]
                   + f_3 * pc_x[k] * isd_77[k];

        t_126[k] = f_10 * hsd_51[k]
                   + f_1 * isp0_37[k]
                   - f_2 * isp1_37[k]
                   + f_3 * pc_y[k] * isd_75[k];

        t_127[k] = f_10 * hsd_45[k]
                   + f_3 * pc_z[k] * isd_75[k];

        t_128[k] = f_10 * hsd_53[k]
                   + f_3 * pc_y[k] * isd_77[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, pa_y, pc_y, pc_z, hsf0_90, hsd_47, \
                         hsd_48, hsd_54, hsf1_90, isp0_38, isp1_38, isd_77, \
                         isd_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_10 * hsd_47[k]
                   + f_1 * isp0_38[k]
                   - f_2 * isp1_38[k]
                   + f_3 * pc_z[k] * isd_77[k];

        t_130[k] = pa_y[k] * hsf0_90[k]
                   - f_4 * pc_y[k] * hsf1_90[k];

        t_131[k] = f_5 * hsd_54[k]
                   + f_3 * pc_y[k] * isd_78[k];

        t_132[k] = f_11 * hsd_48[k]
                   + f_3 * pc_z[k] * isd_78[k];
    }
}

static auto
compute_prim_isf_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t hsf0,
                                                          const size_t hsd, const size_t hsf1,
                                                          const size_t isp0, const size_t isp1,
                                                          const size_t isd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 2.5 / q;
    const auto f_7 = 0.5 / gamma;
    const auto f_8 = 0.5 * p / (gamma * q);
    const auto f_9 = 2.0 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;

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
    auto *t_257 = buffer.data(target + 257);
    auto *t_258 = buffer.data(target + 258);
    auto *t_259 = buffer.data(target + 259);
    auto *t_260 = buffer.data(target + 260);
    auto *t_261 = buffer.data(target + 261);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *hsf0_99 = buffer.data(hsf0 + 99);
    const auto *hsf0_100 = buffer.data(hsf0 + 100);
    const auto *hsf0_140 = buffer.data(hsf0 + 140);
    const auto *hsf0_150 = buffer.data(hsf0 + 150);
    const auto *hsf0_151 = buffer.data(hsf0 + 151);
    const auto *hsf0_156 = buffer.data(hsf0 + 156);
    const auto *hsf0_159 = buffer.data(hsf0 + 159);
    const auto *hsf0_166 = buffer.data(hsf0 + 166);
    const auto *hsf0_169 = buffer.data(hsf0 + 169);
    const auto *hsf0_170 = buffer.data(hsf0 + 170);
    const auto *hsf0_176 = buffer.data(hsf0 + 176);
    const auto *hsf0_179 = buffer.data(hsf0 + 179);
    const auto *hsf0_180 = buffer.data(hsf0 + 180);
    const auto *hsf0_186 = buffer.data(hsf0 + 186);
    const auto *hsf0_189 = buffer.data(hsf0 + 189);
    const auto *hsf0_196 = buffer.data(hsf0 + 196);
    const auto *hsf0_199 = buffer.data(hsf0 + 199);
    const auto *hsf0_200 = buffer.data(hsf0 + 200);
    const auto *hsf0_206 = buffer.data(hsf0 + 206);
    const auto *hsf0_207 = buffer.data(hsf0 + 207);
    const auto *hsf0_209 = buffer.data(hsf0 + 209);

    const auto *hsd_51 = buffer.data(hsd + 51);
    const auto *hsd_54 = buffer.data(hsd + 54);
    const auto *hsd_57 = buffer.data(hsd + 57);
    const auto *hsd_59 = buffer.data(hsd + 59);
    const auto *hsd_60 = buffer.data(hsd + 60);
    const auto *hsd_63 = buffer.data(hsd + 63);
    const auto *hsd_65 = buffer.data(hsd + 65);
    const auto *hsd_66 = buffer.data(hsd + 66);
    const auto *hsd_69 = buffer.data(hsd + 69);
    const auto *hsd_71 = buffer.data(hsd + 71);
    const auto *hsd_72 = buffer.data(hsd + 72);
    const auto *hsd_75 = buffer.data(hsd + 75);
    const auto *hsd_77 = buffer.data(hsd + 77);
    const auto *hsd_78 = buffer.data(hsd + 78);
    const auto *hsd_81 = buffer.data(hsd + 81);
    const auto *hsd_82 = buffer.data(hsd + 82);
    const auto *hsd_83 = buffer.data(hsd + 83);
    const auto *hsd_84 = buffer.data(hsd + 84);
    const auto *hsd_87 = buffer.data(hsd + 87);
    const auto *hsd_89 = buffer.data(hsd + 89);
    const auto *hsd_90 = buffer.data(hsd + 90);
    const auto *hsd_93 = buffer.data(hsd + 93);
    const auto *hsd_95 = buffer.data(hsd + 95);
    const auto *hsd_99 = buffer.data(hsd + 99);
    const auto *hsd_100 = buffer.data(hsd + 100);
    const auto *hsd_101 = buffer.data(hsd + 101);
    const auto *hsd_102 = buffer.data(hsd + 102);
    const auto *hsd_105 = buffer.data(hsd + 105);
    const auto *hsd_106 = buffer.data(hsd + 106);
    const auto *hsd_107 = buffer.data(hsd + 107);
    const auto *hsd_108 = buffer.data(hsd + 108);
    const auto *hsd_111 = buffer.data(hsd + 111);
    const auto *hsd_112 = buffer.data(hsd + 112);
    const auto *hsd_113 = buffer.data(hsd + 113);
    const auto *hsd_117 = buffer.data(hsd + 117);
    const auto *hsd_118 = buffer.data(hsd + 118);
    const auto *hsd_119 = buffer.data(hsd + 119);
    const auto *hsd_120 = buffer.data(hsd + 120);
    const auto *hsd_123 = buffer.data(hsd + 123);
    const auto *hsd_125 = buffer.data(hsd + 125);

    const auto *hsf1_99 = buffer.data(hsf1 + 99);
    const auto *hsf1_100 = buffer.data(hsf1 + 100);
    const auto *hsf1_140 = buffer.data(hsf1 + 140);
    const auto *hsf1_150 = buffer.data(hsf1 + 150);
    const auto *hsf1_151 = buffer.data(hsf1 + 151);
    const auto *hsf1_156 = buffer.data(hsf1 + 156);
    const auto *hsf1_159 = buffer.data(hsf1 + 159);
    const auto *hsf1_166 = buffer.data(hsf1 + 166);
    const auto *hsf1_169 = buffer.data(hsf1 + 169);
    const auto *hsf1_170 = buffer.data(hsf1 + 170);
    const auto *hsf1_176 = buffer.data(hsf1 + 176);
    const auto *hsf1_179 = buffer.data(hsf1 + 179);
    const auto *hsf1_180 = buffer.data(hsf1 + 180);
    const auto *hsf1_186 = buffer.data(hsf1 + 186);
    const auto *hsf1_189 = buffer.data(hsf1 + 189);
    const auto *hsf1_196 = buffer.data(hsf1 + 196);
    const auto *hsf1_199 = buffer.data(hsf1 + 199);
    const auto *hsf1_200 = buffer.data(hsf1 + 200);
    const auto *hsf1_206 = buffer.data(hsf1 + 206);
    const auto *hsf1_207 = buffer.data(hsf1 + 207);
    const auto *hsf1_209 = buffer.data(hsf1 + 209);

    const auto *isp0_40 = buffer.data(isp0 + 40);
    const auto *isp0_42 = buffer.data(isp0 + 42);
    const auto *isp0_43 = buffer.data(isp0 + 43);
    const auto *isp0_44 = buffer.data(isp0 + 44);
    const auto *isp0_63 = buffer.data(isp0 + 63);
    const auto *isp0_64 = buffer.data(isp0 + 64);
    const auto *isp0_65 = buffer.data(isp0 + 65);
    const auto *isp0_68 = buffer.data(isp0 + 68);
    const auto *isp0_69 = buffer.data(isp0 + 69);
    const auto *isp0_70 = buffer.data(isp0 + 70);
    const auto *isp0_71 = buffer.data(isp0 + 71);
    const auto *isp0_72 = buffer.data(isp0 + 72);
    const auto *isp0_73 = buffer.data(isp0 + 73);
    const auto *isp0_74 = buffer.data(isp0 + 74);
    const auto *isp0_75 = buffer.data(isp0 + 75);
    const auto *isp0_76 = buffer.data(isp0 + 76);
    const auto *isp0_77 = buffer.data(isp0 + 77);
    const auto *isp0_79 = buffer.data(isp0 + 79);

    const auto *isp1_40 = buffer.data(isp1 + 40);
    const auto *isp1_42 = buffer.data(isp1 + 42);
    const auto *isp1_43 = buffer.data(isp1 + 43);
    const auto *isp1_44 = buffer.data(isp1 + 44);
    const auto *isp1_63 = buffer.data(isp1 + 63);
    const auto *isp1_64 = buffer.data(isp1 + 64);
    const auto *isp1_65 = buffer.data(isp1 + 65);
    const auto *isp1_68 = buffer.data(isp1 + 68);
    const auto *isp1_69 = buffer.data(isp1 + 69);
    const auto *isp1_70 = buffer.data(isp1 + 70);
    const auto *isp1_71 = buffer.data(isp1 + 71);
    const auto *isp1_72 = buffer.data(isp1 + 72);
    const auto *isp1_73 = buffer.data(isp1 + 73);
    const auto *isp1_74 = buffer.data(isp1 + 74);
    const auto *isp1_75 = buffer.data(isp1 + 75);
    const auto *isp1_76 = buffer.data(isp1 + 76);
    const auto *isp1_77 = buffer.data(isp1 + 77);
    const auto *isp1_79 = buffer.data(isp1 + 79);

    const auto *isd_81 = buffer.data(isd + 81);
    const auto *isd_82 = buffer.data(isd + 82);
    const auto *isd_83 = buffer.data(isd + 83);
    const auto *isd_84 = buffer.data(isd + 84);
    const auto *isd_86 = buffer.data(isd + 86);
    const auto *isd_87 = buffer.data(isd + 87);
    const auto *isd_88 = buffer.data(isd + 88);
    const auto *isd_89 = buffer.data(isd + 89);
    const auto *isd_90 = buffer.data(isd + 90);
    const auto *isd_91 = buffer.data(isd + 91);
    const auto *isd_93 = buffer.data(isd + 93);
    const auto *isd_95 = buffer.data(isd + 95);
    const auto *isd_96 = buffer.data(isd + 96);
    const auto *isd_99 = buffer.data(isd + 99);
    const auto *isd_100 = buffer.data(isd + 100);
    const auto *isd_101 = buffer.data(isd + 101);
    const auto *isd_102 = buffer.data(isd + 102);
    const auto *isd_105 = buffer.data(isd + 105);
    const auto *isd_106 = buffer.data(isd + 106);
    const auto *isd_107 = buffer.data(isd + 107);
    const auto *isd_108 = buffer.data(isd + 108);
    const auto *isd_111 = buffer.data(isd + 111);
    const auto *isd_112 = buffer.data(isd + 112);
    const auto *isd_113 = buffer.data(isd + 113);
    const auto *isd_114 = buffer.data(isd + 114);
    const auto *isd_117 = buffer.data(isd + 117);
    const auto *isd_118 = buffer.data(isd + 118);
    const auto *isd_119 = buffer.data(isd + 119);
    const auto *isd_120 = buffer.data(isd + 120);
    const auto *isd_122 = buffer.data(isd + 122);
    const auto *isd_123 = buffer.data(isd + 123);
    const auto *isd_125 = buffer.data(isd + 125);
    const auto *isd_126 = buffer.data(isd + 126);
    const auto *isd_127 = buffer.data(isd + 127);
    const auto *isd_129 = buffer.data(isd + 129);
    const auto *isd_130 = buffer.data(isd + 130);
    const auto *isd_131 = buffer.data(isd + 131);
    const auto *isd_134 = buffer.data(isd + 134);
    const auto *isd_135 = buffer.data(isd + 135);
    const auto *isd_136 = buffer.data(isd + 136);
    const auto *isd_137 = buffer.data(isd + 137);
    const auto *isd_138 = buffer.data(isd + 138);
    const auto *isd_139 = buffer.data(isd + 139);
    const auto *isd_140 = buffer.data(isd + 140);
    const auto *isd_141 = buffer.data(isd + 141);
    const auto *isd_142 = buffer.data(isd + 142);
    const auto *isd_143 = buffer.data(isd + 143);
    const auto *isd_144 = buffer.data(isd + 144);
    const auto *isd_145 = buffer.data(isd + 145);
    const auto *isd_146 = buffer.data(isd + 146);
    const auto *isd_147 = buffer.data(isd + 147);
    const auto *isd_148 = buffer.data(isd + 148);
    const auto *isd_149 = buffer.data(isd + 149);
    const auto *isd_150 = buffer.data(isd + 150);
    const auto *isd_151 = buffer.data(isd + 151);
    const auto *isd_152 = buffer.data(isd + 152);
    const auto *isd_153 = buffer.data(isd + 153);
    const auto *isd_154 = buffer.data(isd + 154);
    const auto *isd_155 = buffer.data(isd + 155);
    const auto *isd_157 = buffer.data(isd + 157);

#pragma omp simd aligned(t_133, t_134, t_135, t_136, pc_x, pc_y, hsd_57, hsd_81, hsd_82, \
                         hsd_83, isp0_40, isp1_40, isd_81, isd_82, \
                         isd_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_10 * hsd_81[k]
                   + f_3 * pc_x[k] * isd_81[k];

        t_134[k] = f_10 * hsd_82[k]
                   + f_3 * pc_x[k] * isd_82[k];

        t_135[k] = f_10 * hsd_83[k]
                   + f_3 * pc_x[k] * isd_83[k];

        t_136[k] = f_5 * hsd_57[k]
                   + f_1 * isp0_40[k]
                   - f_2 * isp1_40[k]
                   + f_3 * pc_y[k] * isd_81[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, pa_y, pc_y, pc_z, hsf0_99, hsd_51, hsd_59, \
                         hsf1_99, isd_81, isd_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_11 * hsd_51[k]
                   + f_3 * pc_z[k] * isd_81[k];

        t_138[k] = f_5 * hsd_59[k]
                   + f_3 * pc_y[k] * isd_83[k];

        t_139[k] = pa_y[k] * hsf0_99[k]
                   - f_4 * pc_y[k] * hsf1_99[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, pc_x, pc_y, pc_z, hsd_54, hsd_84, \
                         hsd_87, isp0_42, isp1_42, isd_84, isd_86, \
                         isd_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_10 * hsd_84[k]
                   + f_1 * isp0_42[k]
                   - f_2 * isp1_42[k]
                   + f_3 * pc_x[k] * isd_84[k];

        t_141[k] = f_3 * pc_y[k] * isd_84[k];

        t_142[k] = f_9 * hsd_54[k]
                   + f_3 * pc_z[k] * isd_84[k];

        t_143[k] = f_10 * hsd_87[k]
                   + f_3 * pc_x[k] * isd_87[k];

        t_144[k] = f_3 * pc_y[k] * isd_86[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, pc_x, pc_y, hsd_89, isp0_43, isp0_44, \
                         isp1_43, isp1_44, isd_87, isd_88, isd_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_10 * hsd_89[k]
                   + f_3 * pc_x[k] * isd_89[k];

        t_146[k] = f_1 * isp0_43[k]
                   - f_2 * isp1_43[k]
                   + f_3 * pc_y[k] * isd_87[k];

        t_147[k] = f_7 * isp0_44[k]
                   - f_8 * isp1_44[k]
                   + f_3 * pc_y[k] * isd_88[k];

        t_148[k] = f_3 * pc_y[k] * isd_89[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, pa_x, pc_x, pc_y, pc_z, hsf0_150, hsd_59, \
                         hsd_60, hsd_90, hsf1_150, isp0_44, isp1_44, isd_89, \
                         isd_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_9 * hsd_59[k]
                   + f_1 * isp0_44[k]
                   - f_2 * isp1_44[k]
                   + f_3 * pc_z[k] * isd_89[k];

        t_150[k] = pa_x[k] * hsf0_150[k]
                   + f_11 * hsd_90[k]
                   - f_4 * pc_x[k] * hsf1_150[k];

        t_151[k] = f_6 * hsd_60[k]
                   + f_3 * pc_y[k] * isd_90[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, t_156, pa_x, pc_x, pc_z, hsf0_156, \
                         hsd_93, hsd_95, hsf1_156, isd_90, isd_91, isd_93, \
                         isd_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_3 * pc_z[k] * isd_90[k];

        t_153[k] = f_5 * hsd_93[k]
                   + f_3 * pc_x[k] * isd_93[k];

        t_154[k] = f_3 * pc_z[k] * isd_91[k];

        t_155[k] = f_5 * hsd_95[k]
                   + f_3 * pc_x[k] * isd_95[k];

        t_156[k] = pa_x[k] * hsf0_156[k]
                   - f_4 * pc_x[k] * hsf1_156[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, pa_x, pa_z, pc_x, pc_y, pc_z, hsf0_100, \
                         hsf0_159, hsd_65, hsf1_100, hsf1_159, isd_93, \
                         isd_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = f_3 * pc_z[k] * isd_93[k];

        t_158[k] = f_6 * hsd_65[k]
                   + f_3 * pc_y[k] * isd_95[k];

        t_159[k] = pa_x[k] * hsf0_159[k]
                   - f_4 * pc_x[k] * hsf1_159[k];

        t_160[k] = pa_z[k] * hsf0_100[k]
                   - f_4 * pc_z[k] * hsf1_100[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, t_164, pc_x, pc_y, pc_z, hsd_60, hsd_66, hsd_99, \
                         hsd_100, isd_96, isd_99, isd_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = f_9 * hsd_66[k]
                   + f_3 * pc_y[k] * isd_96[k];

        t_162[k] = f_5 * hsd_60[k]
                   + f_3 * pc_z[k] * isd_96[k];

        t_163[k] = f_5 * hsd_99[k]
                   + f_3 * pc_x[k] * isd_99[k];

        t_164[k] = f_5 * hsd_100[k]
                   + f_3 * pc_x[k] * isd_100[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, pa_x, pc_x, pc_y, pc_z, hsf0_166, hsd_63, \
                         hsd_71, hsd_101, hsf1_166, isd_99, isd_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = f_5 * hsd_101[k]
                   + f_3 * pc_x[k] * isd_101[k];

        t_166[k] = pa_x[k] * hsf0_166[k]
                   - f_4 * pc_x[k] * hsf1_166[k];

        t_167[k] = f_5 * hsd_63[k]
                   + f_3 * pc_z[k] * isd_99[k];

        t_168[k] = f_9 * hsd_71[k]
                   + f_3 * pc_y[k] * isd_101[k];
    }

#pragma omp simd aligned(t_169, t_170, t_171, t_172, pa_x, pc_x, pc_y, pc_z, hsf0_169, \
                         hsf0_170, hsd_66, hsd_72, hsd_102, hsf1_169, hsf1_170, \
                         isd_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = pa_x[k] * hsf0_169[k]
                   - f_4 * pc_x[k] * hsf1_169[k];

        t_170[k] = pa_x[k] * hsf0_170[k]
                   + f_11 * hsd_102[k]
                   - f_4 * pc_x[k] * hsf1_170[k];

        t_171[k] = f_11 * hsd_72[k]
                   + f_3 * pc_y[k] * isd_102[k];

        t_172[k] = f_10 * hsd_66[k]
                   + f_3 * pc_z[k] * isd_102[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, pa_x, pc_x, hsf0_176, hsd_105, hsd_106, \
                         hsd_107, hsf1_176, isd_105, isd_106, isd_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = f_5 * hsd_105[k]
                   + f_3 * pc_x[k] * isd_105[k];

        t_174[k] = f_5 * hsd_106[k]
                   + f_3 * pc_x[k] * isd_106[k];

        t_175[k] = f_5 * hsd_107[k]
                   + f_3 * pc_x[k] * isd_107[k];

        t_176[k] = pa_x[k] * hsf0_176[k]
                   - f_4 * pc_x[k] * hsf1_176[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, pa_x, pc_x, pc_y, pc_z, hsf0_179, hsd_69, \
                         hsd_77, hsf1_179, isd_105, isd_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = f_10 * hsd_69[k]
                   + f_3 * pc_z[k] * isd_105[k];

        t_178[k] = f_11 * hsd_77[k]
                   + f_3 * pc_y[k] * isd_107[k];

        t_179[k] = pa_x[k] * hsf0_179[k]
                   - f_4 * pc_x[k] * hsf1_179[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pa_x, pc_x, pc_y, pc_z, hsf0_180, hsd_72, \
                         hsd_78, hsd_108, hsd_111, hsf1_180, isd_108, \
                         isd_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = pa_x[k] * hsf0_180[k]
                   + f_11 * hsd_108[k]
                   - f_4 * pc_x[k] * hsf1_180[k];

        t_181[k] = f_10 * hsd_78[k]
                   + f_3 * pc_y[k] * isd_108[k];

        t_182[k] = f_11 * hsd_72[k]
                   + f_3 * pc_z[k] * isd_108[k];

        t_183[k] = f_5 * hsd_111[k]
                   + f_3 * pc_x[k] * isd_111[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, pa_x, pc_x, pc_z, hsf0_186, hsd_75, \
                         hsd_112, hsd_113, hsf1_186, isd_111, isd_112, \
                         isd_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_5 * hsd_112[k]
                   + f_3 * pc_x[k] * isd_112[k];

        t_185[k] = f_5 * hsd_113[k]
                   + f_3 * pc_x[k] * isd_113[k];

        t_186[k] = pa_x[k] * hsf0_186[k]
                   - f_4 * pc_x[k] * hsf1_186[k];

        t_187[k] = f_11 * hsd_75[k]
                   + f_3 * pc_z[k] * isd_111[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, pa_x, pa_y, pc_x, pc_y, hsf0_140, \
                         hsf0_189, hsd_83, hsd_84, hsf1_140, hsf1_189, isd_113, \
                         isd_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = f_10 * hsd_83[k]
                   + f_3 * pc_y[k] * isd_113[k];

        t_189[k] = pa_x[k] * hsf0_189[k]
                   - f_4 * pc_x[k] * hsf1_189[k];

        t_190[k] = pa_y[k] * hsf0_140[k]
                   - f_4 * pc_y[k] * hsf1_140[k];

        t_191[k] = f_5 * hsd_84[k]
                   + f_3 * pc_y[k] * isd_114[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, pc_x, pc_z, hsd_78, hsd_117, hsd_118, \
                         hsd_119, isd_114, isd_117, isd_118, isd_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_9 * hsd_78[k]
                   + f_3 * pc_z[k] * isd_114[k];

        t_193[k] = f_5 * hsd_117[k]
                   + f_3 * pc_x[k] * isd_117[k];

        t_194[k] = f_5 * hsd_118[k]
                   + f_3 * pc_x[k] * isd_118[k];

        t_195[k] = f_5 * hsd_119[k]
                   + f_3 * pc_x[k] * isd_119[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, pa_x, pc_x, pc_y, pc_z, hsf0_196, \
                         hsf0_199, hsd_81, hsd_89, hsf1_196, hsf1_199, isd_117, \
                         isd_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = pa_x[k] * hsf0_196[k]
                   - f_4 * pc_x[k] * hsf1_196[k];

        t_197[k] = f_9 * hsd_81[k]
                   + f_3 * pc_z[k] * isd_117[k];

        t_198[k] = f_5 * hsd_89[k]
                   + f_3 * pc_y[k] * isd_119[k];

        t_199[k] = pa_x[k] * hsf0_199[k]
                   - f_4 * pc_x[k] * hsf1_199[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, pa_x, pc_x, pc_y, pc_z, hsf0_200, hsd_84, \
                         hsd_120, hsd_123, hsf1_200, isd_120, isd_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = pa_x[k] * hsf0_200[k]
                   + f_11 * hsd_120[k]
                   - f_4 * pc_x[k] * hsf1_200[k];

        t_201[k] = f_3 * pc_y[k] * isd_120[k];

        t_202[k] = f_6 * hsd_84[k]
                   + f_3 * pc_z[k] * isd_120[k];

        t_203[k] = f_5 * hsd_123[k]
                   + f_3 * pc_x[k] * isd_123[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, t_208, pa_x, pc_x, pc_y, hsf0_206, \
                         hsf0_207, hsd_125, hsf1_206, hsf1_207, isd_122, \
                         isd_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_3 * pc_y[k] * isd_122[k];

        t_205[k] = f_5 * hsd_125[k]
                   + f_3 * pc_x[k] * isd_125[k];

        t_206[k] = pa_x[k] * hsf0_206[k]
                   - f_4 * pc_x[k] * hsf1_206[k];

        t_207[k] = pa_x[k] * hsf0_207[k]
                   - f_4 * pc_x[k] * hsf1_207[k];

        t_208[k] = f_3 * pc_y[k] * isd_125[k];
    }

#pragma omp simd aligned(t_209, t_210, t_211, t_212, pa_x, pc_x, pc_z, hsf0_209, hsf1_209, \
                         isp0_63, isp0_64, isp1_63, isp1_64, isd_126, \
                         isd_127 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_209[k] = pa_x[k] * hsf0_209[k]
                   - f_4 * pc_x[k] * hsf1_209[k];

        t_210[k] = f_1 * isp0_63[k]
                   - f_2 * isp1_63[k]
                   + f_3 * pc_x[k] * isd_126[k];

        t_211[k] = f_7 * isp0_64[k]
                   - f_8 * isp1_64[k]
                   + f_3 * pc_x[k] * isd_127[k];

        t_212[k] = f_3 * pc_z[k] * isd_126[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, t_217, t_218, pc_x, pc_y, pc_z, hsd_93, \
                         hsd_95, isp0_64, isp1_64, isd_129, isd_130, \
                         isd_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_3 * pc_x[k] * isd_129[k];

        t_214[k] = f_3 * pc_x[k] * isd_130[k];

        t_215[k] = f_3 * pc_x[k] * isd_131[k];

        t_216[k] = f_0 * hsd_93[k]
                   + f_1 * isp0_64[k]
                   - f_2 * isp1_64[k]
                   + f_3 * pc_y[k] * isd_129[k];

        t_217[k] = f_3 * pc_z[k] * isd_129[k];

        t_218[k] = f_0 * hsd_95[k]
                   + f_3 * pc_y[k] * isd_131[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, pa_z, pc_z, hsf0_150, hsf0_151, hsf1_150, \
                         hsf1_151, isp0_65, isp1_65, isd_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = f_1 * isp0_65[k]
                   - f_2 * isp1_65[k]
                   + f_3 * pc_z[k] * isd_131[k];

        t_220[k] = pa_z[k] * hsf0_150[k]
                   - f_4 * pc_z[k] * hsf1_150[k];

        t_221[k] = pa_z[k] * hsf0_151[k]
                   - f_4 * pc_z[k] * hsf1_151[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, t_226, pa_z, pc_x, pc_z, hsf0_156, \
                         hsf1_156, isp0_68, isp1_68, isd_134, isd_135, isd_136, \
                         isd_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_7 * isp0_68[k]
                   - f_8 * isp1_68[k]
                   + f_3 * pc_x[k] * isd_134[k];

        t_223[k] = f_3 * pc_x[k] * isd_135[k];

        t_224[k] = f_3 * pc_x[k] * isd_136[k];

        t_225[k] = f_3 * pc_x[k] * isd_137[k];

        t_226[k] = pa_z[k] * hsf0_156[k]
                   - f_4 * pc_z[k] * hsf1_156[k];
    }

#pragma omp simd aligned(t_227, t_228, t_229, pc_y, pc_z, hsd_93, hsd_95, hsd_101, isp0_68, \
                         isp1_68, isd_135, isd_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_227[k] = f_5 * hsd_93[k]
                   + f_3 * pc_z[k] * isd_135[k];

        t_228[k] = f_6 * hsd_101[k]
                   + f_3 * pc_y[k] * isd_137[k];

        t_229[k] = f_5 * hsd_95[k]
                   + f_1 * isp0_68[k]
                   - f_2 * isp1_68[k]
                   + f_3 * pc_z[k] * isd_137[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, pc_x, isp0_69, isp0_70, isp0_71, isp1_69, \
                         isp1_70, isp1_71, isd_138, isd_139, isd_140, \
                         isd_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_1 * isp0_69[k]
                   - f_2 * isp1_69[k]
                   + f_3 * pc_x[k] * isd_138[k];

        t_231[k] = f_7 * isp0_70[k]
                   - f_8 * isp1_70[k]
                   + f_3 * pc_x[k] * isd_139[k];

        t_232[k] = f_7 * isp0_71[k]
                   - f_8 * isp1_71[k]
                   + f_3 * pc_x[k] * isd_140[k];

        t_233[k] = f_3 * pc_x[k] * isd_141[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, t_238, pc_x, pc_y, pc_z, hsd_99, hsd_105, \
                         hsd_107, isp0_70, isp1_70, isd_141, isd_142, \
                         isd_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_3 * pc_x[k] * isd_142[k];

        t_235[k] = f_3 * pc_x[k] * isd_143[k];

        t_236[k] = f_9 * hsd_105[k]
                   + f_1 * isp0_70[k]
                   - f_2 * isp1_70[k]
                   + f_3 * pc_y[k] * isd_141[k];

        t_237[k] = f_10 * hsd_99[k]
                   + f_3 * pc_z[k] * isd_141[k];

        t_238[k] = f_9 * hsd_107[k]
                   + f_3 * pc_y[k] * isd_143[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, pc_x, pc_z, hsd_101, isp0_71, isp0_72, isp0_73, \
                         isp1_71, isp1_72, isp1_73, isd_143, isd_144, \
                         isd_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_10 * hsd_101[k]
                   + f_1 * isp0_71[k]
                   - f_2 * isp1_71[k]
                   + f_3 * pc_z[k] * isd_143[k];

        t_240[k] = f_1 * isp0_72[k]
                   - f_2 * isp1_72[k]
                   + f_3 * pc_x[k] * isd_144[k];

        t_241[k] = f_7 * isp0_73[k]
                   - f_8 * isp1_73[k]
                   + f_3 * pc_x[k] * isd_145[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, t_245, t_246, pc_x, pc_y, hsd_111, isp0_73, \
                         isp0_74, isp1_73, isp1_74, isd_146, isd_147, isd_148, \
                         isd_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = f_7 * isp0_74[k]
                   - f_8 * isp1_74[k]
                   + f_3 * pc_x[k] * isd_146[k];

        t_243[k] = f_3 * pc_x[k] * isd_147[k];

        t_244[k] = f_3 * pc_x[k] * isd_148[k];

        t_245[k] = f_3 * pc_x[k] * isd_149[k];

        t_246[k] = f_11 * hsd_111[k]
                   + f_1 * isp0_73[k]
                   - f_2 * isp1_73[k]
                   + f_3 * pc_y[k] * isd_147[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, pc_y, pc_z, hsd_105, hsd_107, hsd_113, isp0_74, \
                         isp1_74, isd_147, isd_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = f_11 * hsd_105[k]
                   + f_3 * pc_z[k] * isd_147[k];

        t_248[k] = f_11 * hsd_113[k]
                   + f_3 * pc_y[k] * isd_149[k];

        t_249[k] = f_11 * hsd_107[k]
                   + f_1 * isp0_74[k]
                   - f_2 * isp1_74[k]
                   + f_3 * pc_z[k] * isd_149[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, pc_x, isp0_75, isp0_76, isp0_77, isp1_75, \
                         isp1_76, isp1_77, isd_150, isd_151, isd_152, \
                         isd_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = f_1 * isp0_75[k]
                   - f_2 * isp1_75[k]
                   + f_3 * pc_x[k] * isd_150[k];

        t_251[k] = f_7 * isp0_76[k]
                   - f_8 * isp1_76[k]
                   + f_3 * pc_x[k] * isd_151[k];

        t_252[k] = f_7 * isp0_77[k]
                   - f_8 * isp1_77[k]
                   + f_3 * pc_x[k] * isd_152[k];

        t_253[k] = f_3 * pc_x[k] * isd_153[k];
    }

#pragma omp simd aligned(t_254, t_255, t_256, t_257, t_258, pc_x, pc_y, pc_z, hsd_111, \
                         hsd_117, hsd_119, isp0_76, isp1_76, isd_153, isd_154, \
                         isd_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_254[k] = f_3 * pc_x[k] * isd_154[k];

        t_255[k] = f_3 * pc_x[k] * isd_155[k];

        t_256[k] = f_10 * hsd_117[k]
                   + f_1 * isp0_76[k]
                   - f_2 * isp1_76[k]
                   + f_3 * pc_y[k] * isd_153[k];

        t_257[k] = f_9 * hsd_111[k]
                   + f_3 * pc_z[k] * isd_153[k];

        t_258[k] = f_10 * hsd_119[k]
                   + f_3 * pc_y[k] * isd_155[k];
    }

#pragma omp simd aligned(t_259, t_260, t_261, pa_y, pc_x, pc_y, pc_z, hsf0_200, hsd_113, \
                         hsf1_200, isp0_77, isp0_79, isp1_77, isp1_79, isd_155, \
                         isd_157 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_259[k] = f_9 * hsd_113[k]
                   + f_1 * isp0_77[k]
                   - f_2 * isp1_77[k]
                   + f_3 * pc_z[k] * isd_155[k];

        t_260[k] = pa_y[k] * hsf0_200[k]
                   - f_4 * pc_y[k] * hsf1_200[k];

        t_261[k] = f_7 * isp0_79[k]
                   - f_8 * isp1_79[k]
                   + f_3 * pc_x[k] * isd_157[k];
    }
}

static auto
compute_prim_isf_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t hsf0,
                                                          const size_t hsd, const size_t hsf1,
                                                          const size_t isp0, const size_t isp1,
                                                          const size_t isd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 2.5 / q;
    const auto f_7 = 0.5 / gamma;
    const auto f_8 = 0.5 * p / (gamma * q);
    const auto f_11 = 1.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *hsf0_202 = buffer.data(hsf0 + 202);
    const auto *hsf0_206 = buffer.data(hsf0 + 206);
    const auto *hsf0_209 = buffer.data(hsf0 + 209);

    const auto *hsd_117 = buffer.data(hsd + 117);
    const auto *hsd_123 = buffer.data(hsd + 123);
    const auto *hsd_125 = buffer.data(hsd + 125);

    const auto *hsf1_202 = buffer.data(hsf1 + 202);
    const auto *hsf1_206 = buffer.data(hsf1 + 206);
    const auto *hsf1_209 = buffer.data(hsf1 + 209);

    const auto *isp0_81 = buffer.data(isp0 + 81);
    const auto *isp0_82 = buffer.data(isp0 + 82);
    const auto *isp0_83 = buffer.data(isp0 + 83);

    const auto *isp1_81 = buffer.data(isp1 + 81);
    const auto *isp1_82 = buffer.data(isp1 + 82);
    const auto *isp1_83 = buffer.data(isp1 + 83);

    const auto *isd_159 = buffer.data(isd + 159);
    const auto *isd_160 = buffer.data(isd + 160);
    const auto *isd_161 = buffer.data(isd + 161);
    const auto *isd_162 = buffer.data(isd + 162);
    const auto *isd_164 = buffer.data(isd + 164);
    const auto *isd_165 = buffer.data(isd + 165);
    const auto *isd_166 = buffer.data(isd + 166);
    const auto *isd_167 = buffer.data(isd + 167);

#pragma omp simd aligned(t_262, t_263, t_264, t_265, t_266, pa_y, pc_x, pc_y, hsf0_202, \
                         hsf0_206, hsd_123, hsf1_202, hsf1_206, isd_159, isd_160, \
                         isd_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_262[k] = pa_y[k] * hsf0_202[k]
                   - f_4 * pc_y[k] * hsf1_202[k];

        t_263[k] = f_3 * pc_x[k] * isd_159[k];

        t_264[k] = f_3 * pc_x[k] * isd_160[k];

        t_265[k] = f_3 * pc_x[k] * isd_161[k];

        t_266[k] = pa_y[k] * hsf0_206[k]
                   + f_11 * hsd_123[k]
                   - f_4 * pc_y[k] * hsf1_206[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, pa_y, pc_y, pc_z, hsf0_209, hsd_117, hsd_125, \
                         hsf1_209, isd_159, isd_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_6 * hsd_117[k]
                   + f_3 * pc_z[k] * isd_159[k];

        t_268[k] = f_5 * hsd_125[k]
                   + f_3 * pc_y[k] * isd_161[k];

        t_269[k] = pa_y[k] * hsf0_209[k]
                   - f_4 * pc_y[k] * hsf1_209[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, pc_x, pc_y, isp0_81, isp0_83, \
                         isp1_81, isp1_83, isd_162, isd_164, isd_165, \
                         isd_166 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_1 * isp0_81[k]
                   - f_2 * isp1_81[k]
                   + f_3 * pc_x[k] * isd_162[k];

        t_271[k] = f_3 * pc_y[k] * isd_162[k];

        t_272[k] = f_7 * isp0_83[k]
                   - f_8 * isp1_83[k]
                   + f_3 * pc_x[k] * isd_164[k];

        t_273[k] = f_3 * pc_x[k] * isd_165[k];

        t_274[k] = f_3 * pc_x[k] * isd_166[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, pc_x, pc_y, pc_z, hsd_125, \
                         isp0_82, isp0_83, isp1_82, isp1_83, isd_165, isd_166, \
                         isd_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = f_3 * pc_x[k] * isd_167[k];

        t_276[k] = f_1 * isp0_82[k]
                   - f_2 * isp1_82[k]
                   + f_3 * pc_y[k] * isd_165[k];

        t_277[k] = f_7 * isp0_83[k]
                   - f_8 * isp1_83[k]
                   + f_3 * pc_y[k] * isd_166[k];

        t_278[k] = f_3 * pc_y[k] * isd_167[k];

        t_279[k] = f_0 * hsd_125[k]
                   + f_1 * isp0_83[k]
                   - f_2 * isp1_83[k]
                   + f_3 * pc_z[k] * isd_167[k];
    }
}

auto
compute_prim_isf_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t hsf0, const size_t hsd,
                                                   const size_t hsf1, const size_t isp0,
                                                   const size_t isp1, const size_t isd,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_isf_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, hsf0, hsd,
                                                              hsf1, isp0, isp1, isd, ncols,
                                                              gamma, p, q);

    compute_prim_isf_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, hsf0, hsd,
                                                              hsf1, isp0, isp1, isd, ncols,
                                                              gamma, p, q);

    compute_prim_isf_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, hsf0, hsd,
                                                              hsf1, isp0, isp1, isd, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
