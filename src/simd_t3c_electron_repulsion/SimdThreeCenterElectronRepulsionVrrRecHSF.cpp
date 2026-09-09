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


#include "SimdThreeCenterElectronRepulsionVrrRecHSF.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_hsf_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t gsf0,
                                                          const size_t gsd, const size_t gsf1,
                                                          const size_t hsp0, const size_t hsp1,
                                                          const size_t hsd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 2.0 / q;
    const auto f_7 = 0.5 / gamma;
    const auto f_8 = 0.5 * p / (gamma * q);
    const auto f_9 = 1.5 / q;
    const auto f_10 = 1.0 / q;

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
    auto *t_133 = buffer.data(target + 133);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *gsf0_0 = buffer.data(gsf0 + 0);
    const auto *gsf0_6 = buffer.data(gsf0 + 6);
    const auto *gsf0_9 = buffer.data(gsf0 + 9);
    const auto *gsf0_16 = buffer.data(gsf0 + 16);
    const auto *gsf0_20 = buffer.data(gsf0 + 20);
    const auto *gsf0_29 = buffer.data(gsf0 + 29);
    const auto *gsf0_30 = buffer.data(gsf0 + 30);
    const auto *gsf0_36 = buffer.data(gsf0 + 36);
    const auto *gsf0_50 = buffer.data(gsf0 + 50);
    const auto *gsf0_59 = buffer.data(gsf0 + 59);
    const auto *gsf0_60 = buffer.data(gsf0 + 60);
    const auto *gsf0_90 = buffer.data(gsf0 + 90);
    const auto *gsf0_100 = buffer.data(gsf0 + 100);
    const auto *gsf0_106 = buffer.data(gsf0 + 106);
    const auto *gsf0_109 = buffer.data(gsf0 + 109);
    const auto *gsf0_116 = buffer.data(gsf0 + 116);
    const auto *gsf0_119 = buffer.data(gsf0 + 119);
    const auto *gsf0_120 = buffer.data(gsf0 + 120);
    const auto *gsf0_126 = buffer.data(gsf0 + 126);
    const auto *gsf0_129 = buffer.data(gsf0 + 129);

    const auto *gsd_0 = buffer.data(gsd + 0);
    const auto *gsd_3 = buffer.data(gsd + 3);
    const auto *gsd_5 = buffer.data(gsd + 5);
    const auto *gsd_6 = buffer.data(gsd + 6);
    const auto *gsd_9 = buffer.data(gsd + 9);
    const auto *gsd_11 = buffer.data(gsd + 11);
    const auto *gsd_12 = buffer.data(gsd + 12);
    const auto *gsd_15 = buffer.data(gsd + 15);
    const auto *gsd_17 = buffer.data(gsd + 17);
    const auto *gsd_18 = buffer.data(gsd + 18);
    const auto *gsd_21 = buffer.data(gsd + 21);
    const auto *gsd_23 = buffer.data(gsd + 23);
    const auto *gsd_24 = buffer.data(gsd + 24);
    const auto *gsd_27 = buffer.data(gsd + 27);
    const auto *gsd_28 = buffer.data(gsd + 28);
    const auto *gsd_29 = buffer.data(gsd + 29);
    const auto *gsd_30 = buffer.data(gsd + 30);
    const auto *gsd_33 = buffer.data(gsd + 33);
    const auto *gsd_35 = buffer.data(gsd + 35);
    const auto *gsd_36 = buffer.data(gsd + 36);
    const auto *gsd_39 = buffer.data(gsd + 39);
    const auto *gsd_41 = buffer.data(gsd + 41);
    const auto *gsd_42 = buffer.data(gsd + 42);
    const auto *gsd_45 = buffer.data(gsd + 45);
    const auto *gsd_46 = buffer.data(gsd + 46);
    const auto *gsd_47 = buffer.data(gsd + 47);
    const auto *gsd_48 = buffer.data(gsd + 48);
    const auto *gsd_51 = buffer.data(gsd + 51);
    const auto *gsd_52 = buffer.data(gsd + 52);
    const auto *gsd_53 = buffer.data(gsd + 53);
    const auto *gsd_54 = buffer.data(gsd + 54);
    const auto *gsd_57 = buffer.data(gsd + 57);
    const auto *gsd_59 = buffer.data(gsd + 59);
    const auto *gsd_60 = buffer.data(gsd + 60);
    const auto *gsd_63 = buffer.data(gsd + 63);
    const auto *gsd_65 = buffer.data(gsd + 65);
    const auto *gsd_69 = buffer.data(gsd + 69);
    const auto *gsd_70 = buffer.data(gsd + 70);
    const auto *gsd_71 = buffer.data(gsd + 71);
    const auto *gsd_72 = buffer.data(gsd + 72);
    const auto *gsd_75 = buffer.data(gsd + 75);
    const auto *gsd_76 = buffer.data(gsd + 76);
    const auto *gsd_77 = buffer.data(gsd + 77);
    const auto *gsd_81 = buffer.data(gsd + 81);

    const auto *gsf1_0 = buffer.data(gsf1 + 0);
    const auto *gsf1_6 = buffer.data(gsf1 + 6);
    const auto *gsf1_9 = buffer.data(gsf1 + 9);
    const auto *gsf1_16 = buffer.data(gsf1 + 16);
    const auto *gsf1_20 = buffer.data(gsf1 + 20);
    const auto *gsf1_29 = buffer.data(gsf1 + 29);
    const auto *gsf1_30 = buffer.data(gsf1 + 30);
    const auto *gsf1_36 = buffer.data(gsf1 + 36);
    const auto *gsf1_50 = buffer.data(gsf1 + 50);
    const auto *gsf1_59 = buffer.data(gsf1 + 59);
    const auto *gsf1_60 = buffer.data(gsf1 + 60);
    const auto *gsf1_90 = buffer.data(gsf1 + 90);
    const auto *gsf1_100 = buffer.data(gsf1 + 100);
    const auto *gsf1_106 = buffer.data(gsf1 + 106);
    const auto *gsf1_109 = buffer.data(gsf1 + 109);
    const auto *gsf1_116 = buffer.data(gsf1 + 116);
    const auto *gsf1_119 = buffer.data(gsf1 + 119);
    const auto *gsf1_120 = buffer.data(gsf1 + 120);
    const auto *gsf1_126 = buffer.data(gsf1 + 126);
    const auto *gsf1_129 = buffer.data(gsf1 + 129);

    const auto *hsp0_0 = buffer.data(hsp0 + 0);
    const auto *hsp0_1 = buffer.data(hsp0 + 1);
    const auto *hsp0_2 = buffer.data(hsp0 + 2);
    const auto *hsp0_4 = buffer.data(hsp0 + 4);
    const auto *hsp0_8 = buffer.data(hsp0 + 8);
    const auto *hsp0_9 = buffer.data(hsp0 + 9);
    const auto *hsp0_10 = buffer.data(hsp0 + 10);
    const auto *hsp0_11 = buffer.data(hsp0 + 11);
    const auto *hsp0_15 = buffer.data(hsp0 + 15);
    const auto *hsp0_16 = buffer.data(hsp0 + 16);
    const auto *hsp0_17 = buffer.data(hsp0 + 17);
    const auto *hsp0_18 = buffer.data(hsp0 + 18);
    const auto *hsp0_19 = buffer.data(hsp0 + 19);
    const auto *hsp0_20 = buffer.data(hsp0 + 20);
    const auto *hsp0_23 = buffer.data(hsp0 + 23);
    const auto *hsp0_25 = buffer.data(hsp0 + 25);
    const auto *hsp0_27 = buffer.data(hsp0 + 27);
    const auto *hsp0_28 = buffer.data(hsp0 + 28);
    const auto *hsp0_29 = buffer.data(hsp0 + 29);

    const auto *hsp1_0 = buffer.data(hsp1 + 0);
    const auto *hsp1_1 = buffer.data(hsp1 + 1);
    const auto *hsp1_2 = buffer.data(hsp1 + 2);
    const auto *hsp1_4 = buffer.data(hsp1 + 4);
    const auto *hsp1_8 = buffer.data(hsp1 + 8);
    const auto *hsp1_9 = buffer.data(hsp1 + 9);
    const auto *hsp1_10 = buffer.data(hsp1 + 10);
    const auto *hsp1_11 = buffer.data(hsp1 + 11);
    const auto *hsp1_15 = buffer.data(hsp1 + 15);
    const auto *hsp1_16 = buffer.data(hsp1 + 16);
    const auto *hsp1_17 = buffer.data(hsp1 + 17);
    const auto *hsp1_18 = buffer.data(hsp1 + 18);
    const auto *hsp1_19 = buffer.data(hsp1 + 19);
    const auto *hsp1_20 = buffer.data(hsp1 + 20);
    const auto *hsp1_23 = buffer.data(hsp1 + 23);
    const auto *hsp1_25 = buffer.data(hsp1 + 25);
    const auto *hsp1_27 = buffer.data(hsp1 + 27);
    const auto *hsp1_28 = buffer.data(hsp1 + 28);
    const auto *hsp1_29 = buffer.data(hsp1 + 29);

    const auto *hsd_0 = buffer.data(hsd + 0);
    const auto *hsd_2 = buffer.data(hsd + 2);
    const auto *hsd_3 = buffer.data(hsd + 3);
    const auto *hsd_5 = buffer.data(hsd + 5);
    const auto *hsd_6 = buffer.data(hsd + 6);
    const auto *hsd_7 = buffer.data(hsd + 7);
    const auto *hsd_9 = buffer.data(hsd + 9);
    const auto *hsd_11 = buffer.data(hsd + 11);
    const auto *hsd_12 = buffer.data(hsd + 12);
    const auto *hsd_14 = buffer.data(hsd + 14);
    const auto *hsd_15 = buffer.data(hsd + 15);
    const auto *hsd_16 = buffer.data(hsd + 16);
    const auto *hsd_17 = buffer.data(hsd + 17);
    const auto *hsd_18 = buffer.data(hsd + 18);
    const auto *hsd_19 = buffer.data(hsd + 19);
    const auto *hsd_21 = buffer.data(hsd + 21);
    const auto *hsd_23 = buffer.data(hsd + 23);
    const auto *hsd_24 = buffer.data(hsd + 24);
    const auto *hsd_27 = buffer.data(hsd + 27);
    const auto *hsd_28 = buffer.data(hsd + 28);
    const auto *hsd_29 = buffer.data(hsd + 29);
    const auto *hsd_30 = buffer.data(hsd + 30);
    const auto *hsd_32 = buffer.data(hsd + 32);
    const auto *hsd_33 = buffer.data(hsd + 33);
    const auto *hsd_34 = buffer.data(hsd + 34);
    const auto *hsd_35 = buffer.data(hsd + 35);
    const auto *hsd_36 = buffer.data(hsd + 36);
    const auto *hsd_37 = buffer.data(hsd + 37);
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
    const auto *hsd_56 = buffer.data(hsd + 56);
    const auto *hsd_57 = buffer.data(hsd + 57);
    const auto *hsd_58 = buffer.data(hsd + 58);
    const auto *hsd_59 = buffer.data(hsd + 59);
    const auto *hsd_60 = buffer.data(hsd + 60);
    const auto *hsd_61 = buffer.data(hsd + 61);
    const auto *hsd_63 = buffer.data(hsd + 63);
    const auto *hsd_65 = buffer.data(hsd + 65);
    const auto *hsd_66 = buffer.data(hsd + 66);
    const auto *hsd_69 = buffer.data(hsd + 69);
    const auto *hsd_70 = buffer.data(hsd + 70);
    const auto *hsd_71 = buffer.data(hsd + 71);
    const auto *hsd_72 = buffer.data(hsd + 72);
    const auto *hsd_75 = buffer.data(hsd + 75);
    const auto *hsd_76 = buffer.data(hsd + 76);
    const auto *hsd_77 = buffer.data(hsd + 77);
    const auto *hsd_78 = buffer.data(hsd + 78);
    const auto *hsd_81 = buffer.data(hsd + 81);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pc_x, pc_y, pc_z, gsd_0, gsd_3, hsp0_0, \
                         hsp1_0, hsd_0, hsd_2, hsd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gsd_0[k]
                 + f_1 * hsp0_0[k]
                 - f_2 * hsp1_0[k]
                 + f_3 * pc_x[k] * hsd_0[k];

        t_1[k] = f_3 * pc_y[k] * hsd_0[k];

        t_2[k] = f_3 * pc_z[k] * hsd_0[k];

        t_3[k] = f_0 * gsd_3[k]
                 + f_3 * pc_x[k] * hsd_3[k];

        t_4[k] = f_3 * pc_y[k] * hsd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pc_x, pc_y, pc_z, gsd_5, hsp0_1, hsp0_2, \
                         hsp1_1, hsp1_2, hsd_3, hsd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * gsd_5[k]
                 + f_3 * pc_x[k] * hsd_5[k];

        t_6[k] = f_1 * hsp0_1[k]
                 - f_2 * hsp1_1[k]
                 + f_3 * pc_y[k] * hsd_3[k];

        t_7[k] = f_3 * pc_z[k] * hsd_3[k];

        t_8[k] = f_3 * pc_y[k] * hsd_5[k];

        t_9[k] = f_1 * hsp0_2[k]
                 - f_2 * hsp1_2[k]
                 + f_3 * pc_z[k] * hsd_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_y, pc_x, pc_y, pc_z, gsf0_0, gsd_0, \
                         gsd_9, gsf1_0, hsd_6, hsd_7, hsd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_y[k] * gsf0_0[k]
                  - f_4 * pc_y[k] * gsf1_0[k];

        t_11[k] = f_5 * gsd_0[k]
                  + f_3 * pc_y[k] * hsd_6[k];

        t_12[k] = f_3 * pc_z[k] * hsd_6[k];

        t_13[k] = f_6 * gsd_9[k]
                  + f_3 * pc_x[k] * hsd_9[k];

        t_14[k] = f_3 * pc_z[k] * hsd_7[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pc_x, pc_y, pc_z, gsd_3, gsd_5, gsd_11, \
                         hsp0_4, hsp1_4, hsd_9, hsd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_6 * gsd_11[k]
                  + f_3 * pc_x[k] * hsd_11[k];

        t_16[k] = f_5 * gsd_3[k]
                  + f_1 * hsp0_4[k]
                  - f_2 * hsp1_4[k]
                  + f_3 * pc_y[k] * hsd_9[k];

        t_17[k] = f_3 * pc_z[k] * hsd_9[k];

        t_18[k] = f_5 * gsd_5[k]
                  + f_3 * pc_y[k] * hsd_11[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_y, pa_z, pc_y, pc_z, gsf0_0, gsf0_9, \
                         gsd_0, gsf1_0, gsf1_9, hsd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pa_y[k] * gsf0_9[k]
                  - f_4 * pc_y[k] * gsf1_9[k];

        t_20[k] = pa_z[k] * gsf0_0[k]
                  - f_4 * pc_z[k] * gsf1_0[k];

        t_21[k] = f_3 * pc_y[k] * hsd_12[k];

        t_22[k] = f_5 * gsd_0[k]
                  + f_3 * pc_z[k] * hsd_12[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_z, pc_x, pc_y, pc_z, gsf0_6, gsd_15, \
                         gsd_17, gsf1_6, hsd_14, hsd_15, hsd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_6 * gsd_15[k]
                  + f_3 * pc_x[k] * hsd_15[k];

        t_24[k] = f_3 * pc_y[k] * hsd_14[k];

        t_25[k] = f_6 * gsd_17[k]
                  + f_3 * pc_x[k] * hsd_17[k];

        t_26[k] = pa_z[k] * gsf0_6[k]
                  - f_4 * pc_z[k] * gsf1_6[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pc_x, pc_y, pc_z, gsd_5, gsd_18, hsp0_8, \
                         hsp0_9, hsp1_8, hsp1_9, hsd_16, hsd_17, \
                         hsd_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_7 * hsp0_8[k]
                  - f_8 * hsp1_8[k]
                  + f_3 * pc_y[k] * hsd_16[k];

        t_28[k] = f_3 * pc_y[k] * hsd_17[k];

        t_29[k] = f_5 * gsd_5[k]
                  + f_1 * hsp0_8[k]
                  - f_2 * hsp1_8[k]
                  + f_3 * pc_z[k] * hsd_17[k];

        t_30[k] = f_9 * gsd_18[k]
                  + f_1 * hsp0_9[k]
                  - f_2 * hsp1_9[k]
                  + f_3 * pc_x[k] * hsd_18[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, pc_x, pc_y, pc_z, gsd_6, gsd_21, \
                         gsd_23, hsd_18, hsd_19, hsd_21, hsd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_10 * gsd_6[k]
                  + f_3 * pc_y[k] * hsd_18[k];

        t_32[k] = f_3 * pc_z[k] * hsd_18[k];

        t_33[k] = f_9 * gsd_21[k]
                  + f_3 * pc_x[k] * hsd_21[k];

        t_34[k] = f_3 * pc_z[k] * hsd_19[k];

        t_35[k] = f_9 * gsd_23[k]
                  + f_3 * pc_x[k] * hsd_23[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pc_y, pc_z, gsd_9, gsd_11, hsp0_10, hsp0_11, \
                         hsp1_10, hsp1_11, hsd_21, hsd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_10 * gsd_9[k]
                  + f_1 * hsp0_10[k]
                  - f_2 * hsp1_10[k]
                  + f_3 * pc_y[k] * hsd_21[k];

        t_37[k] = f_3 * pc_z[k] * hsd_21[k];

        t_38[k] = f_10 * gsd_11[k]
                  + f_3 * pc_y[k] * hsd_23[k];

        t_39[k] = f_1 * hsp0_11[k]
                  - f_2 * hsp1_11[k]
                  + f_3 * pc_z[k] * hsd_23[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_y, pc_x, pc_y, pc_z, gsf0_20, gsd_6, \
                         gsd_12, gsd_27, gsf1_20, hsd_24, hsd_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pa_y[k] * gsf0_20[k]
                  - f_4 * pc_y[k] * gsf1_20[k];

        t_41[k] = f_5 * gsd_12[k]
                  + f_3 * pc_y[k] * hsd_24[k];

        t_42[k] = f_5 * gsd_6[k]
                  + f_3 * pc_z[k] * hsd_24[k];

        t_43[k] = f_9 * gsd_27[k]
                  + f_3 * pc_x[k] * hsd_27[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_z, pc_x, pc_z, gsf0_16, gsd_9, gsd_28, \
                         gsd_29, gsf1_16, hsd_27, hsd_28, hsd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_9 * gsd_28[k]
                  + f_3 * pc_x[k] * hsd_28[k];

        t_45[k] = f_9 * gsd_29[k]
                  + f_3 * pc_x[k] * hsd_29[k];

        t_46[k] = pa_z[k] * gsf0_16[k]
                  - f_4 * pc_z[k] * gsf1_16[k];

        t_47[k] = f_5 * gsd_9[k]
                  + f_3 * pc_z[k] * hsd_27[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_y, pc_x, pc_y, gsf0_29, gsd_17, gsd_30, \
                         gsf1_29, hsp0_15, hsp1_15, hsd_29, hsd_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_5 * gsd_17[k]
                  + f_3 * pc_y[k] * hsd_29[k];

        t_49[k] = pa_y[k] * gsf0_29[k]
                  - f_4 * pc_y[k] * gsf1_29[k];

        t_50[k] = f_9 * gsd_30[k]
                  + f_1 * hsp0_15[k]
                  - f_2 * hsp1_15[k]
                  + f_3 * pc_x[k] * hsd_30[k];

        t_51[k] = f_3 * pc_y[k] * hsd_30[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pc_x, pc_y, pc_z, gsd_12, gsd_33, gsd_35, \
                         hsd_30, hsd_32, hsd_33, hsd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_10 * gsd_12[k]
                  + f_3 * pc_z[k] * hsd_30[k];

        t_53[k] = f_9 * gsd_33[k]
                  + f_3 * pc_x[k] * hsd_33[k];

        t_54[k] = f_3 * pc_y[k] * hsd_32[k];

        t_55[k] = f_9 * gsd_35[k]
                  + f_3 * pc_x[k] * hsd_35[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pc_y, pc_z, gsd_17, hsp0_16, hsp0_17, \
                         hsp1_16, hsp1_17, hsd_33, hsd_34, hsd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_1 * hsp0_16[k]
                  - f_2 * hsp1_16[k]
                  + f_3 * pc_y[k] * hsd_33[k];

        t_57[k] = f_7 * hsp0_17[k]
                  - f_8 * hsp1_17[k]
                  + f_3 * pc_y[k] * hsd_34[k];

        t_58[k] = f_3 * pc_y[k] * hsd_35[k];

        t_59[k] = f_10 * gsd_17[k]
                  + f_1 * hsp0_17[k]
                  - f_2 * hsp1_17[k]
                  + f_3 * pc_z[k] * hsd_35[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, pc_x, pc_y, pc_z, gsd_18, gsd_36, \
                         gsd_39, hsp0_18, hsp1_18, hsd_36, hsd_37, \
                         hsd_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_10 * gsd_36[k]
                  + f_1 * hsp0_18[k]
                  - f_2 * hsp1_18[k]
                  + f_3 * pc_x[k] * hsd_36[k];

        t_61[k] = f_9 * gsd_18[k]
                  + f_3 * pc_y[k] * hsd_36[k];

        t_62[k] = f_3 * pc_z[k] * hsd_36[k];

        t_63[k] = f_10 * gsd_39[k]
                  + f_3 * pc_x[k] * hsd_39[k];

        t_64[k] = f_3 * pc_z[k] * hsd_37[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pc_x, pc_y, pc_z, gsd_21, gsd_23, gsd_41, \
                         hsp0_19, hsp1_19, hsd_39, hsd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_10 * gsd_41[k]
                  + f_3 * pc_x[k] * hsd_41[k];

        t_66[k] = f_9 * gsd_21[k]
                  + f_1 * hsp0_19[k]
                  - f_2 * hsp1_19[k]
                  + f_3 * pc_y[k] * hsd_39[k];

        t_67[k] = f_3 * pc_z[k] * hsd_39[k];

        t_68[k] = f_9 * gsd_23[k]
                  + f_3 * pc_y[k] * hsd_41[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pa_z, pc_y, pc_z, gsf0_30, gsd_18, gsd_24, \
                         gsf1_30, hsp0_20, hsp1_20, hsd_41, hsd_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_1 * hsp0_20[k]
                  - f_2 * hsp1_20[k]
                  + f_3 * pc_z[k] * hsd_41[k];

        t_70[k] = pa_z[k] * gsf0_30[k]
                  - f_4 * pc_z[k] * gsf1_30[k];

        t_71[k] = f_10 * gsd_24[k]
                  + f_3 * pc_y[k] * hsd_42[k];

        t_72[k] = f_5 * gsd_18[k]
                  + f_3 * pc_z[k] * hsd_42[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pa_z, pc_x, pc_z, gsf0_36, gsd_45, gsd_46, \
                         gsd_47, gsf1_36, hsd_45, hsd_46, hsd_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_10 * gsd_45[k]
                  + f_3 * pc_x[k] * hsd_45[k];

        t_74[k] = f_10 * gsd_46[k]
                  + f_3 * pc_x[k] * hsd_46[k];

        t_75[k] = f_10 * gsd_47[k]
                  + f_3 * pc_x[k] * hsd_47[k];

        t_76[k] = pa_z[k] * gsf0_36[k]
                  - f_4 * pc_z[k] * gsf1_36[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, pa_y, pc_y, pc_z, gsf0_50, gsd_21, gsd_23, \
                         gsd_29, gsf1_50, hsp0_23, hsp1_23, hsd_45, \
                         hsd_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_5 * gsd_21[k]
                  + f_3 * pc_z[k] * hsd_45[k];

        t_78[k] = f_10 * gsd_29[k]
                  + f_3 * pc_y[k] * hsd_47[k];

        t_79[k] = f_5 * gsd_23[k]
                  + f_1 * hsp0_23[k]
                  - f_2 * hsp1_23[k]
                  + f_3 * pc_z[k] * hsd_47[k];

        t_80[k] = pa_y[k] * gsf0_50[k]
                  - f_4 * pc_y[k] * gsf1_50[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pc_x, pc_y, pc_z, gsd_24, gsd_30, gsd_51, \
                         gsd_52, hsd_48, hsd_51, hsd_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_5 * gsd_30[k]
                  + f_3 * pc_y[k] * hsd_48[k];

        t_82[k] = f_10 * gsd_24[k]
                  + f_3 * pc_z[k] * hsd_48[k];

        t_83[k] = f_10 * gsd_51[k]
                  + f_3 * pc_x[k] * hsd_51[k];

        t_84[k] = f_10 * gsd_52[k]
                  + f_3 * pc_x[k] * hsd_52[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pc_x, pc_y, pc_z, gsd_27, gsd_33, gsd_35, \
                         gsd_53, hsp0_25, hsp1_25, hsd_51, hsd_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_10 * gsd_53[k]
                  + f_3 * pc_x[k] * hsd_53[k];

        t_86[k] = f_5 * gsd_33[k]
                  + f_1 * hsp0_25[k]
                  - f_2 * hsp1_25[k]
                  + f_3 * pc_y[k] * hsd_51[k];

        t_87[k] = f_10 * gsd_27[k]
                  + f_3 * pc_z[k] * hsd_51[k];

        t_88[k] = f_5 * gsd_35[k]
                  + f_3 * pc_y[k] * hsd_53[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, pa_y, pc_x, pc_y, pc_z, gsf0_59, gsd_30, \
                         gsd_54, gsf1_59, hsp0_27, hsp1_27, hsd_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = pa_y[k] * gsf0_59[k]
                  - f_4 * pc_y[k] * gsf1_59[k];

        t_90[k] = f_10 * gsd_54[k]
                  + f_1 * hsp0_27[k]
                  - f_2 * hsp1_27[k]
                  + f_3 * pc_x[k] * hsd_54[k];

        t_91[k] = f_3 * pc_y[k] * hsd_54[k];

        t_92[k] = f_9 * gsd_30[k]
                  + f_3 * pc_z[k] * hsd_54[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, pc_x, pc_y, gsd_57, gsd_59, hsp0_28, hsp1_28, \
                         hsd_56, hsd_57, hsd_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_10 * gsd_57[k]
                  + f_3 * pc_x[k] * hsd_57[k];

        t_94[k] = f_3 * pc_y[k] * hsd_56[k];

        t_95[k] = f_10 * gsd_59[k]
                  + f_3 * pc_x[k] * hsd_59[k];

        t_96[k] = f_1 * hsp0_28[k]
                  - f_2 * hsp1_28[k]
                  + f_3 * pc_y[k] * hsd_57[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pa_x, pc_x, pc_y, pc_z, gsf0_100, gsd_35, \
                         gsd_60, gsf1_100, hsp0_29, hsp1_29, hsd_58, \
                         hsd_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_7 * hsp0_29[k]
                  - f_8 * hsp1_29[k]
                  + f_3 * pc_y[k] * hsd_58[k];

        t_98[k] = f_3 * pc_y[k] * hsd_59[k];

        t_99[k] = f_9 * gsd_35[k]
                  + f_1 * hsp0_29[k]
                  - f_2 * hsp1_29[k]
                  + f_3 * pc_z[k] * hsd_59[k];

        t_100[k] = pa_x[k] * gsf0_100[k]
                   + f_9 * gsd_60[k]
                   - f_4 * pc_x[k] * gsf1_100[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, t_105, pc_x, pc_y, pc_z, gsd_36, gsd_63, \
                         gsd_65, hsd_60, hsd_61, hsd_63, hsd_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_6 * gsd_36[k]
                   + f_3 * pc_y[k] * hsd_60[k];

        t_102[k] = f_3 * pc_z[k] * hsd_60[k];

        t_103[k] = f_5 * gsd_63[k]
                   + f_3 * pc_x[k] * hsd_63[k];

        t_104[k] = f_3 * pc_z[k] * hsd_61[k];

        t_105[k] = f_5 * gsd_65[k]
                   + f_3 * pc_x[k] * hsd_65[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, pa_x, pc_x, pc_y, pc_z, gsf0_106, \
                         gsf0_109, gsd_41, gsf1_106, gsf1_109, hsd_63, \
                         hsd_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = pa_x[k] * gsf0_106[k]
                   - f_4 * pc_x[k] * gsf1_106[k];

        t_107[k] = f_3 * pc_z[k] * hsd_63[k];

        t_108[k] = f_6 * gsd_41[k]
                   + f_3 * pc_y[k] * hsd_65[k];

        t_109[k] = pa_x[k] * gsf0_109[k]
                   - f_4 * pc_x[k] * gsf1_109[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, pa_z, pc_x, pc_y, pc_z, gsf0_60, gsd_36, \
                         gsd_42, gsd_69, gsf1_60, hsd_66, hsd_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = pa_z[k] * gsf0_60[k]
                   - f_4 * pc_z[k] * gsf1_60[k];

        t_111[k] = f_9 * gsd_42[k]
                   + f_3 * pc_y[k] * hsd_66[k];

        t_112[k] = f_5 * gsd_36[k]
                   + f_3 * pc_z[k] * hsd_66[k];

        t_113[k] = f_5 * gsd_69[k]
                   + f_3 * pc_x[k] * hsd_69[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pa_x, pc_x, pc_z, gsf0_116, gsd_39, \
                         gsd_70, gsd_71, gsf1_116, hsd_69, hsd_70, \
                         hsd_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_5 * gsd_70[k]
                   + f_3 * pc_x[k] * hsd_70[k];

        t_115[k] = f_5 * gsd_71[k]
                   + f_3 * pc_x[k] * hsd_71[k];

        t_116[k] = pa_x[k] * gsf0_116[k]
                   - f_4 * pc_x[k] * gsf1_116[k];

        t_117[k] = f_5 * gsd_39[k]
                   + f_3 * pc_z[k] * hsd_69[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, pa_x, pc_x, pc_y, gsf0_119, gsf0_120, \
                         gsd_47, gsd_48, gsd_72, gsf1_119, gsf1_120, hsd_71, \
                         hsd_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_9 * gsd_47[k]
                   + f_3 * pc_y[k] * hsd_71[k];

        t_119[k] = pa_x[k] * gsf0_119[k]
                   - f_4 * pc_x[k] * gsf1_119[k];

        t_120[k] = pa_x[k] * gsf0_120[k]
                   + f_9 * gsd_72[k]
                   - f_4 * pc_x[k] * gsf1_120[k];

        t_121[k] = f_10 * gsd_48[k]
                   + f_3 * pc_y[k] * hsd_72[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, pc_x, pc_z, gsd_42, gsd_75, gsd_76, \
                         gsd_77, hsd_72, hsd_75, hsd_76, hsd_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_10 * gsd_42[k]
                   + f_3 * pc_z[k] * hsd_72[k];

        t_123[k] = f_5 * gsd_75[k]
                   + f_3 * pc_x[k] * hsd_75[k];

        t_124[k] = f_5 * gsd_76[k]
                   + f_3 * pc_x[k] * hsd_76[k];

        t_125[k] = f_5 * gsd_77[k]
                   + f_3 * pc_x[k] * hsd_77[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, pa_x, pc_x, pc_y, pc_z, gsf0_126, \
                         gsf0_129, gsd_45, gsd_53, gsf1_126, gsf1_129, hsd_75, \
                         hsd_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = pa_x[k] * gsf0_126[k]
                   - f_4 * pc_x[k] * gsf1_126[k];

        t_127[k] = f_10 * gsd_45[k]
                   + f_3 * pc_z[k] * hsd_75[k];

        t_128[k] = f_10 * gsd_53[k]
                   + f_3 * pc_y[k] * hsd_77[k];

        t_129[k] = pa_x[k] * gsf0_129[k]
                   - f_4 * pc_x[k] * gsf1_129[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, pa_y, pc_x, pc_y, pc_z, gsf0_90, gsd_48, \
                         gsd_54, gsd_81, gsf1_90, hsd_78, hsd_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = pa_y[k] * gsf0_90[k]
                   - f_4 * pc_y[k] * gsf1_90[k];

        t_131[k] = f_5 * gsd_54[k]
                   + f_3 * pc_y[k] * hsd_78[k];

        t_132[k] = f_9 * gsd_48[k]
                   + f_3 * pc_z[k] * hsd_78[k];

        t_133[k] = f_5 * gsd_81[k]
                   + f_3 * pc_x[k] * hsd_81[k];
    }
}

static auto
compute_prim_hsf_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t gsf0,
                                                          const size_t gsd, const size_t gsf1,
                                                          const size_t hsp0, const size_t hsp1,
                                                          const size_t hsd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 2.0 / q;
    const auto f_7 = 0.5 / gamma;
    const auto f_8 = 0.5 * p / (gamma * q);
    const auto f_9 = 1.5 / q;
    const auto f_10 = 1.0 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *gsf0_100 = buffer.data(gsf0 + 100);
    const auto *gsf0_101 = buffer.data(gsf0 + 101);
    const auto *gsf0_106 = buffer.data(gsf0 + 106);
    const auto *gsf0_136 = buffer.data(gsf0 + 136);
    const auto *gsf0_139 = buffer.data(gsf0 + 139);
    const auto *gsf0_140 = buffer.data(gsf0 + 140);
    const auto *gsf0_142 = buffer.data(gsf0 + 142);
    const auto *gsf0_146 = buffer.data(gsf0 + 146);
    const auto *gsf0_147 = buffer.data(gsf0 + 147);
    const auto *gsf0_149 = buffer.data(gsf0 + 149);

    const auto *gsd_51 = buffer.data(gsd + 51);
    const auto *gsd_54 = buffer.data(gsd + 54);
    const auto *gsd_59 = buffer.data(gsd + 59);
    const auto *gsd_63 = buffer.data(gsd + 63);
    const auto *gsd_65 = buffer.data(gsd + 65);
    const auto *gsd_69 = buffer.data(gsd + 69);
    const auto *gsd_71 = buffer.data(gsd + 71);
    const auto *gsd_75 = buffer.data(gsd + 75);
    const auto *gsd_77 = buffer.data(gsd + 77);
    const auto *gsd_81 = buffer.data(gsd + 81);
    const auto *gsd_82 = buffer.data(gsd + 82);
    const auto *gsd_83 = buffer.data(gsd + 83);
    const auto *gsd_84 = buffer.data(gsd + 84);
    const auto *gsd_87 = buffer.data(gsd + 87);
    const auto *gsd_89 = buffer.data(gsd + 89);

    const auto *gsf1_100 = buffer.data(gsf1 + 100);
    const auto *gsf1_101 = buffer.data(gsf1 + 101);
    const auto *gsf1_106 = buffer.data(gsf1 + 106);
    const auto *gsf1_136 = buffer.data(gsf1 + 136);
    const auto *gsf1_139 = buffer.data(gsf1 + 139);
    const auto *gsf1_140 = buffer.data(gsf1 + 140);
    const auto *gsf1_142 = buffer.data(gsf1 + 142);
    const auto *gsf1_146 = buffer.data(gsf1 + 146);
    const auto *gsf1_147 = buffer.data(gsf1 + 147);
    const auto *gsf1_149 = buffer.data(gsf1 + 149);

    const auto *hsp0_45 = buffer.data(hsp0 + 45);
    const auto *hsp0_46 = buffer.data(hsp0 + 46);
    const auto *hsp0_47 = buffer.data(hsp0 + 47);
    const auto *hsp0_50 = buffer.data(hsp0 + 50);
    const auto *hsp0_51 = buffer.data(hsp0 + 51);
    const auto *hsp0_52 = buffer.data(hsp0 + 52);
    const auto *hsp0_53 = buffer.data(hsp0 + 53);
    const auto *hsp0_54 = buffer.data(hsp0 + 54);
    const auto *hsp0_55 = buffer.data(hsp0 + 55);
    const auto *hsp0_56 = buffer.data(hsp0 + 56);
    const auto *hsp0_58 = buffer.data(hsp0 + 58);
    const auto *hsp0_60 = buffer.data(hsp0 + 60);
    const auto *hsp0_61 = buffer.data(hsp0 + 61);
    const auto *hsp0_62 = buffer.data(hsp0 + 62);

    const auto *hsp1_45 = buffer.data(hsp1 + 45);
    const auto *hsp1_46 = buffer.data(hsp1 + 46);
    const auto *hsp1_47 = buffer.data(hsp1 + 47);
    const auto *hsp1_50 = buffer.data(hsp1 + 50);
    const auto *hsp1_51 = buffer.data(hsp1 + 51);
    const auto *hsp1_52 = buffer.data(hsp1 + 52);
    const auto *hsp1_53 = buffer.data(hsp1 + 53);
    const auto *hsp1_54 = buffer.data(hsp1 + 54);
    const auto *hsp1_55 = buffer.data(hsp1 + 55);
    const auto *hsp1_56 = buffer.data(hsp1 + 56);
    const auto *hsp1_58 = buffer.data(hsp1 + 58);
    const auto *hsp1_60 = buffer.data(hsp1 + 60);
    const auto *hsp1_61 = buffer.data(hsp1 + 61);
    const auto *hsp1_62 = buffer.data(hsp1 + 62);

    const auto *hsd_81 = buffer.data(hsd + 81);
    const auto *hsd_82 = buffer.data(hsd + 82);
    const auto *hsd_83 = buffer.data(hsd + 83);
    const auto *hsd_84 = buffer.data(hsd + 84);
    const auto *hsd_86 = buffer.data(hsd + 86);
    const auto *hsd_87 = buffer.data(hsd + 87);
    const auto *hsd_89 = buffer.data(hsd + 89);
    const auto *hsd_90 = buffer.data(hsd + 90);
    const auto *hsd_91 = buffer.data(hsd + 91);
    const auto *hsd_93 = buffer.data(hsd + 93);
    const auto *hsd_94 = buffer.data(hsd + 94);
    const auto *hsd_95 = buffer.data(hsd + 95);
    const auto *hsd_98 = buffer.data(hsd + 98);
    const auto *hsd_99 = buffer.data(hsd + 99);
    const auto *hsd_100 = buffer.data(hsd + 100);
    const auto *hsd_101 = buffer.data(hsd + 101);
    const auto *hsd_102 = buffer.data(hsd + 102);
    const auto *hsd_103 = buffer.data(hsd + 103);
    const auto *hsd_104 = buffer.data(hsd + 104);
    const auto *hsd_105 = buffer.data(hsd + 105);
    const auto *hsd_106 = buffer.data(hsd + 106);
    const auto *hsd_107 = buffer.data(hsd + 107);
    const auto *hsd_108 = buffer.data(hsd + 108);
    const auto *hsd_109 = buffer.data(hsd + 109);
    const auto *hsd_110 = buffer.data(hsd + 110);
    const auto *hsd_111 = buffer.data(hsd + 111);
    const auto *hsd_112 = buffer.data(hsd + 112);
    const auto *hsd_113 = buffer.data(hsd + 113);
    const auto *hsd_115 = buffer.data(hsd + 115);
    const auto *hsd_117 = buffer.data(hsd + 117);
    const auto *hsd_118 = buffer.data(hsd + 118);
    const auto *hsd_119 = buffer.data(hsd + 119);
    const auto *hsd_120 = buffer.data(hsd + 120);
    const auto *hsd_122 = buffer.data(hsd + 122);
    const auto *hsd_123 = buffer.data(hsd + 123);
    const auto *hsd_124 = buffer.data(hsd + 124);
    const auto *hsd_125 = buffer.data(hsd + 125);

#pragma omp simd aligned(t_134, t_135, t_136, t_137, pa_x, pc_x, pc_z, gsf0_136, gsd_51, \
                         gsd_82, gsd_83, gsf1_136, hsd_81, hsd_82, \
                         hsd_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_5 * gsd_82[k]
                   + f_3 * pc_x[k] * hsd_82[k];

        t_135[k] = f_5 * gsd_83[k]
                   + f_3 * pc_x[k] * hsd_83[k];

        t_136[k] = pa_x[k] * gsf0_136[k]
                   - f_4 * pc_x[k] * gsf1_136[k];

        t_137[k] = f_9 * gsd_51[k]
                   + f_3 * pc_z[k] * hsd_81[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, pa_x, pc_x, pc_y, gsf0_139, gsf0_140, \
                         gsd_59, gsd_84, gsf1_139, gsf1_140, hsd_83, \
                         hsd_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_5 * gsd_59[k]
                   + f_3 * pc_y[k] * hsd_83[k];

        t_139[k] = pa_x[k] * gsf0_139[k]
                   - f_4 * pc_x[k] * gsf1_139[k];

        t_140[k] = pa_x[k] * gsf0_140[k]
                   + f_9 * gsd_84[k]
                   - f_4 * pc_x[k] * gsf1_140[k];

        t_141[k] = f_3 * pc_y[k] * hsd_84[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, pc_x, pc_y, pc_z, gsd_54, gsd_87, gsd_89, \
                         hsd_84, hsd_86, hsd_87, hsd_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_6 * gsd_54[k]
                   + f_3 * pc_z[k] * hsd_84[k];

        t_143[k] = f_5 * gsd_87[k]
                   + f_3 * pc_x[k] * hsd_87[k];

        t_144[k] = f_3 * pc_y[k] * hsd_86[k];

        t_145[k] = f_5 * gsd_89[k]
                   + f_3 * pc_x[k] * hsd_89[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pa_x, pc_x, pc_y, gsf0_146, gsf0_147, \
                         gsf0_149, gsf1_146, gsf1_147, gsf1_149, \
                         hsd_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = pa_x[k] * gsf0_146[k]
                   - f_4 * pc_x[k] * gsf1_146[k];

        t_147[k] = pa_x[k] * gsf0_147[k]
                   - f_4 * pc_x[k] * gsf1_147[k];

        t_148[k] = f_3 * pc_y[k] * hsd_89[k];

        t_149[k] = pa_x[k] * gsf0_149[k]
                   - f_4 * pc_x[k] * gsf1_149[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, pc_x, pc_z, hsp0_45, hsp0_46, \
                         hsp1_45, hsp1_46, hsd_90, hsd_91, hsd_93, \
                         hsd_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_1 * hsp0_45[k]
                   - f_2 * hsp1_45[k]
                   + f_3 * pc_x[k] * hsd_90[k];

        t_151[k] = f_7 * hsp0_46[k]
                   - f_8 * hsp1_46[k]
                   + f_3 * pc_x[k] * hsd_91[k];

        t_152[k] = f_3 * pc_z[k] * hsd_90[k];

        t_153[k] = f_3 * pc_x[k] * hsd_93[k];

        t_154[k] = f_3 * pc_x[k] * hsd_94[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, pc_x, pc_y, pc_z, gsd_63, gsd_65, \
                         hsp0_46, hsp0_47, hsp1_46, hsp1_47, hsd_93, \
                         hsd_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_3 * pc_x[k] * hsd_95[k];

        t_156[k] = f_0 * gsd_63[k]
                   + f_1 * hsp0_46[k]
                   - f_2 * hsp1_46[k]
                   + f_3 * pc_y[k] * hsd_93[k];

        t_157[k] = f_3 * pc_z[k] * hsd_93[k];

        t_158[k] = f_0 * gsd_65[k]
                   + f_3 * pc_y[k] * hsd_95[k];

        t_159[k] = f_1 * hsp0_47[k]
                   - f_2 * hsp1_47[k]
                   + f_3 * pc_z[k] * hsd_95[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, pa_z, pc_x, pc_z, gsf0_100, gsf0_101, \
                         gsf1_100, gsf1_101, hsp0_50, hsp1_50, hsd_98, \
                         hsd_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = pa_z[k] * gsf0_100[k]
                   - f_4 * pc_z[k] * gsf1_100[k];

        t_161[k] = pa_z[k] * gsf0_101[k]
                   - f_4 * pc_z[k] * gsf1_101[k];

        t_162[k] = f_7 * hsp0_50[k]
                   - f_8 * hsp1_50[k]
                   + f_3 * pc_x[k] * hsd_98[k];

        t_163[k] = f_3 * pc_x[k] * hsd_99[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, t_168, pa_z, pc_x, pc_y, pc_z, gsf0_106, \
                         gsd_63, gsd_71, gsf1_106, hsd_99, hsd_100, \
                         hsd_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_3 * pc_x[k] * hsd_100[k];

        t_165[k] = f_3 * pc_x[k] * hsd_101[k];

        t_166[k] = pa_z[k] * gsf0_106[k]
                   - f_4 * pc_z[k] * gsf1_106[k];

        t_167[k] = f_5 * gsd_63[k]
                   + f_3 * pc_z[k] * hsd_99[k];

        t_168[k] = f_6 * gsd_71[k]
                   + f_3 * pc_y[k] * hsd_101[k];
    }

#pragma omp simd aligned(t_169, t_170, t_171, pc_x, pc_z, gsd_65, hsp0_50, hsp0_51, hsp0_52, \
                         hsp1_50, hsp1_51, hsp1_52, hsd_101, hsd_102, \
                         hsd_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = f_5 * gsd_65[k]
                   + f_1 * hsp0_50[k]
                   - f_2 * hsp1_50[k]
                   + f_3 * pc_z[k] * hsd_101[k];

        t_170[k] = f_1 * hsp0_51[k]
                   - f_2 * hsp1_51[k]
                   + f_3 * pc_x[k] * hsd_102[k];

        t_171[k] = f_7 * hsp0_52[k]
                   - f_8 * hsp1_52[k]
                   + f_3 * pc_x[k] * hsd_103[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, t_176, pc_x, pc_y, gsd_75, hsp0_52, \
                         hsp0_53, hsp1_52, hsp1_53, hsd_104, hsd_105, hsd_106, \
                         hsd_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_7 * hsp0_53[k]
                   - f_8 * hsp1_53[k]
                   + f_3 * pc_x[k] * hsd_104[k];

        t_173[k] = f_3 * pc_x[k] * hsd_105[k];

        t_174[k] = f_3 * pc_x[k] * hsd_106[k];

        t_175[k] = f_3 * pc_x[k] * hsd_107[k];

        t_176[k] = f_9 * gsd_75[k]
                   + f_1 * hsp0_52[k]
                   - f_2 * hsp1_52[k]
                   + f_3 * pc_y[k] * hsd_105[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, pc_y, pc_z, gsd_69, gsd_71, gsd_77, hsp0_53, \
                         hsp1_53, hsd_105, hsd_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = f_10 * gsd_69[k]
                   + f_3 * pc_z[k] * hsd_105[k];

        t_178[k] = f_9 * gsd_77[k]
                   + f_3 * pc_y[k] * hsd_107[k];

        t_179[k] = f_10 * gsd_71[k]
                   + f_1 * hsp0_53[k]
                   - f_2 * hsp1_53[k]
                   + f_3 * pc_z[k] * hsd_107[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pc_x, hsp0_54, hsp0_55, hsp0_56, hsp1_54, \
                         hsp1_55, hsp1_56, hsd_108, hsd_109, hsd_110, \
                         hsd_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_1 * hsp0_54[k]
                   - f_2 * hsp1_54[k]
                   + f_3 * pc_x[k] * hsd_108[k];

        t_181[k] = f_7 * hsp0_55[k]
                   - f_8 * hsp1_55[k]
                   + f_3 * pc_x[k] * hsd_109[k];

        t_182[k] = f_7 * hsp0_56[k]
                   - f_8 * hsp1_56[k]
                   + f_3 * pc_x[k] * hsd_110[k];

        t_183[k] = f_3 * pc_x[k] * hsd_111[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, t_188, pc_x, pc_y, pc_z, gsd_75, gsd_81, \
                         gsd_83, hsp0_55, hsp1_55, hsd_111, hsd_112, \
                         hsd_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_3 * pc_x[k] * hsd_112[k];

        t_185[k] = f_3 * pc_x[k] * hsd_113[k];

        t_186[k] = f_10 * gsd_81[k]
                   + f_1 * hsp0_55[k]
                   - f_2 * hsp1_55[k]
                   + f_3 * pc_y[k] * hsd_111[k];

        t_187[k] = f_9 * gsd_75[k]
                   + f_3 * pc_z[k] * hsd_111[k];

        t_188[k] = f_10 * gsd_83[k]
                   + f_3 * pc_y[k] * hsd_113[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, pa_y, pc_x, pc_y, pc_z, gsf0_140, gsd_77, \
                         gsf1_140, hsp0_56, hsp0_58, hsp1_56, hsp1_58, hsd_113, \
                         hsd_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_9 * gsd_77[k]
                   + f_1 * hsp0_56[k]
                   - f_2 * hsp1_56[k]
                   + f_3 * pc_z[k] * hsd_113[k];

        t_190[k] = pa_y[k] * gsf0_140[k]
                   - f_4 * pc_y[k] * gsf1_140[k];

        t_191[k] = f_7 * hsp0_58[k]
                   - f_8 * hsp1_58[k]
                   + f_3 * pc_x[k] * hsd_115[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, t_196, pa_y, pc_x, pc_y, gsf0_142, \
                         gsf0_146, gsd_87, gsf1_142, gsf1_146, hsd_117, hsd_118, \
                         hsd_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = pa_y[k] * gsf0_142[k]
                   - f_4 * pc_y[k] * gsf1_142[k];

        t_193[k] = f_3 * pc_x[k] * hsd_117[k];

        t_194[k] = f_3 * pc_x[k] * hsd_118[k];

        t_195[k] = f_3 * pc_x[k] * hsd_119[k];

        t_196[k] = pa_y[k] * gsf0_146[k]
                   + f_9 * gsd_87[k]
                   - f_4 * pc_y[k] * gsf1_146[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, pa_y, pc_y, pc_z, gsf0_149, gsd_81, gsd_89, \
                         gsf1_149, hsd_117, hsd_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_6 * gsd_81[k]
                   + f_3 * pc_z[k] * hsd_117[k];

        t_198[k] = f_5 * gsd_89[k]
                   + f_3 * pc_y[k] * hsd_119[k];

        t_199[k] = pa_y[k] * gsf0_149[k]
                   - f_4 * pc_y[k] * gsf1_149[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, pc_x, pc_y, hsp0_60, hsp0_62, \
                         hsp1_60, hsp1_62, hsd_120, hsd_122, hsd_123, \
                         hsd_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_1 * hsp0_60[k]
                   - f_2 * hsp1_60[k]
                   + f_3 * pc_x[k] * hsd_120[k];

        t_201[k] = f_3 * pc_y[k] * hsd_120[k];

        t_202[k] = f_7 * hsp0_62[k]
                   - f_8 * hsp1_62[k]
                   + f_3 * pc_x[k] * hsd_122[k];

        t_203[k] = f_3 * pc_x[k] * hsd_123[k];

        t_204[k] = f_3 * pc_x[k] * hsd_124[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, pc_x, pc_y, pc_z, gsd_89, hsp0_61, \
                         hsp0_62, hsp1_61, hsp1_62, hsd_123, hsd_124, \
                         hsd_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = f_3 * pc_x[k] * hsd_125[k];

        t_206[k] = f_1 * hsp0_61[k]
                   - f_2 * hsp1_61[k]
                   + f_3 * pc_y[k] * hsd_123[k];

        t_207[k] = f_7 * hsp0_62[k]
                   - f_8 * hsp1_62[k]
                   + f_3 * pc_y[k] * hsd_124[k];

        t_208[k] = f_3 * pc_y[k] * hsd_125[k];

        t_209[k] = f_0 * gsd_89[k]
                   + f_1 * hsp0_62[k]
                   - f_2 * hsp1_62[k]
                   + f_3 * pc_z[k] * hsd_125[k];
    }
}

auto
compute_prim_hsf_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t gsf0, const size_t gsd,
                                                   const size_t gsf1, const size_t hsp0,
                                                   const size_t hsp1, const size_t hsd,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_hsf_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, gsf0, gsd,
                                                              gsf1, hsp0, hsp1, hsd, ncols,
                                                              gamma, p, q);

    compute_prim_hsf_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, gsf0, gsd,
                                                              gsf1, hsp0, hsp1, hsd, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
