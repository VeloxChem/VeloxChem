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


#include "SimdThreeCenterElectronRepulsionVrrRecGPP.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_gpp_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t fpp0, const size_t fps,
                                                   const size_t fpp1, const size_t gss,
                                                   const size_t gps, const size_t ncols,
                                                   const double gamma, const double p,
                                                   const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 0.5 / q;
    const auto f_2 = p / q;
    const auto f_3 = gamma / q;
    const auto f_4 = 1.5 / q;
    const auto f_5 = 1.0 / q;

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
    auto *t_134 = buffer.data(target + 134);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *fpp0_0 = buffer.data(fpp0 + 0);
    const auto *fpp0_13 = buffer.data(fpp0 + 13);
    const auto *fpp0_18 = buffer.data(fpp0 + 18);
    const auto *fpp0_26 = buffer.data(fpp0 + 26);
    const auto *fpp0_27 = buffer.data(fpp0 + 27);
    const auto *fpp0_45 = buffer.data(fpp0 + 45);
    const auto *fpp0_58 = buffer.data(fpp0 + 58);
    const auto *fpp0_67 = buffer.data(fpp0 + 67);
    const auto *fpp0_71 = buffer.data(fpp0 + 71);
    const auto *fpp0_76 = buffer.data(fpp0 + 76);
    const auto *fpp0_80 = buffer.data(fpp0 + 80);
    const auto *fpp0_89 = buffer.data(fpp0 + 89);

    const auto *fps_0 = buffer.data(fps + 0);
    const auto *fps_1 = buffer.data(fps + 1);
    const auto *fps_2 = buffer.data(fps + 2);
    const auto *fps_3 = buffer.data(fps + 3);
    const auto *fps_4 = buffer.data(fps + 4);
    const auto *fps_5 = buffer.data(fps + 5);
    const auto *fps_6 = buffer.data(fps + 6);
    const auto *fps_7 = buffer.data(fps + 7);
    const auto *fps_8 = buffer.data(fps + 8);
    const auto *fps_9 = buffer.data(fps + 9);
    const auto *fps_10 = buffer.data(fps + 10);
    const auto *fps_11 = buffer.data(fps + 11);
    const auto *fps_12 = buffer.data(fps + 12);
    const auto *fps_13 = buffer.data(fps + 13);
    const auto *fps_14 = buffer.data(fps + 14);
    const auto *fps_15 = buffer.data(fps + 15);
    const auto *fps_16 = buffer.data(fps + 16);
    const auto *fps_17 = buffer.data(fps + 17);
    const auto *fps_18 = buffer.data(fps + 18);
    const auto *fps_19 = buffer.data(fps + 19);
    const auto *fps_20 = buffer.data(fps + 20);
    const auto *fps_21 = buffer.data(fps + 21);
    const auto *fps_22 = buffer.data(fps + 22);
    const auto *fps_23 = buffer.data(fps + 23);
    const auto *fps_24 = buffer.data(fps + 24);
    const auto *fps_25 = buffer.data(fps + 25);
    const auto *fps_26 = buffer.data(fps + 26);
    const auto *fps_27 = buffer.data(fps + 27);
    const auto *fps_28 = buffer.data(fps + 28);
    const auto *fps_29 = buffer.data(fps + 29);

    const auto *fpp1_0 = buffer.data(fpp1 + 0);
    const auto *fpp1_13 = buffer.data(fpp1 + 13);
    const auto *fpp1_18 = buffer.data(fpp1 + 18);
    const auto *fpp1_26 = buffer.data(fpp1 + 26);
    const auto *fpp1_27 = buffer.data(fpp1 + 27);
    const auto *fpp1_45 = buffer.data(fpp1 + 45);
    const auto *fpp1_58 = buffer.data(fpp1 + 58);
    const auto *fpp1_67 = buffer.data(fpp1 + 67);
    const auto *fpp1_71 = buffer.data(fpp1 + 71);
    const auto *fpp1_76 = buffer.data(fpp1 + 76);
    const auto *fpp1_80 = buffer.data(fpp1 + 80);
    const auto *fpp1_89 = buffer.data(fpp1 + 89);

    const auto *gss_0 = buffer.data(gss + 0);
    const auto *gss_1 = buffer.data(gss + 1);
    const auto *gss_2 = buffer.data(gss + 2);
    const auto *gss_3 = buffer.data(gss + 3);
    const auto *gss_5 = buffer.data(gss + 5);
    const auto *gss_6 = buffer.data(gss + 6);
    const auto *gss_9 = buffer.data(gss + 9);
    const auto *gss_10 = buffer.data(gss + 10);
    const auto *gss_11 = buffer.data(gss + 11);
    const auto *gss_12 = buffer.data(gss + 12);
    const auto *gss_13 = buffer.data(gss + 13);
    const auto *gss_14 = buffer.data(gss + 14);

    const auto *gps_0 = buffer.data(gps + 0);
    const auto *gps_1 = buffer.data(gps + 1);
    const auto *gps_2 = buffer.data(gps + 2);
    const auto *gps_3 = buffer.data(gps + 3);
    const auto *gps_4 = buffer.data(gps + 4);
    const auto *gps_5 = buffer.data(gps + 5);
    const auto *gps_6 = buffer.data(gps + 6);
    const auto *gps_7 = buffer.data(gps + 7);
    const auto *gps_8 = buffer.data(gps + 8);
    const auto *gps_9 = buffer.data(gps + 9);
    const auto *gps_10 = buffer.data(gps + 10);
    const auto *gps_11 = buffer.data(gps + 11);
    const auto *gps_12 = buffer.data(gps + 12);
    const auto *gps_13 = buffer.data(gps + 13);
    const auto *gps_14 = buffer.data(gps + 14);
    const auto *gps_15 = buffer.data(gps + 15);
    const auto *gps_16 = buffer.data(gps + 16);
    const auto *gps_17 = buffer.data(gps + 17);
    const auto *gps_18 = buffer.data(gps + 18);
    const auto *gps_19 = buffer.data(gps + 19);
    const auto *gps_20 = buffer.data(gps + 20);
    const auto *gps_21 = buffer.data(gps + 21);
    const auto *gps_22 = buffer.data(gps + 22);
    const auto *gps_23 = buffer.data(gps + 23);
    const auto *gps_24 = buffer.data(gps + 24);
    const auto *gps_25 = buffer.data(gps + 25);
    const auto *gps_26 = buffer.data(gps + 26);
    const auto *gps_27 = buffer.data(gps + 27);
    const auto *gps_28 = buffer.data(gps + 28);
    const auto *gps_29 = buffer.data(gps + 29);
    const auto *gps_30 = buffer.data(gps + 30);
    const auto *gps_31 = buffer.data(gps + 31);
    const auto *gps_32 = buffer.data(gps + 32);
    const auto *gps_33 = buffer.data(gps + 33);
    const auto *gps_34 = buffer.data(gps + 34);
    const auto *gps_35 = buffer.data(gps + 35);
    const auto *gps_36 = buffer.data(gps + 36);
    const auto *gps_37 = buffer.data(gps + 37);
    const auto *gps_38 = buffer.data(gps + 38);
    const auto *gps_39 = buffer.data(gps + 39);
    const auto *gps_40 = buffer.data(gps + 40);
    const auto *gps_41 = buffer.data(gps + 41);
    const auto *gps_42 = buffer.data(gps + 42);
    const auto *gps_43 = buffer.data(gps + 43);
    const auto *gps_44 = buffer.data(gps + 44);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, fps_0, fps_1, gss_0, \
                         gps_0, gps_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fps_0[k]
                 + f_1 * gss_0[k]
                 + f_2 * pc_x[k] * gps_0[k];

        t_1[k] = f_2 * pc_y[k] * gps_0[k];

        t_2[k] = f_2 * pc_z[k] * gps_0[k];

        t_3[k] = f_0 * fps_1[k]
                 + f_2 * pc_x[k] * gps_1[k];

        t_4[k] = f_1 * gss_0[k]
                 + f_2 * pc_y[k] * gps_1[k];

        t_5[k] = f_2 * pc_z[k] * gps_1[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pa_y, pc_x, pc_y, pc_z, fpp0_0, fps_0, \
                         fps_2, fpp1_0, gss_0, gps_2, gps_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * fps_2[k]
                 + f_2 * pc_x[k] * gps_2[k];

        t_7[k] = f_2 * pc_y[k] * gps_2[k];

        t_8[k] = f_1 * gss_0[k]
                 + f_2 * pc_z[k] * gps_2[k];

        t_9[k] = pa_y[k] * fpp0_0[k]
                 - f_3 * pc_y[k] * fpp1_0[k];

        t_10[k] = f_1 * fps_0[k]
                  + f_2 * pc_y[k] * gps_3[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pc_x, pc_y, pc_z, fps_1, fps_4, fps_5, \
                         gss_1, gps_3, gps_4, gps_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_2 * pc_z[k] * gps_3[k];

        t_12[k] = f_4 * fps_4[k]
                  + f_2 * pc_x[k] * gps_4[k];

        t_13[k] = f_1 * fps_1[k]
                  + f_1 * gss_1[k]
                  + f_2 * pc_y[k] * gps_4[k];

        t_14[k] = f_2 * pc_z[k] * gps_4[k];

        t_15[k] = f_4 * fps_5[k]
                  + f_2 * pc_x[k] * gps_5[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pa_z, pc_y, pc_z, fpp0_0, fps_0, fps_2, \
                         fpp1_0, gss_1, gps_5, gps_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_1 * fps_2[k]
                  + f_2 * pc_y[k] * gps_5[k];

        t_17[k] = f_1 * gss_1[k]
                  + f_2 * pc_z[k] * gps_5[k];

        t_18[k] = pa_z[k] * fpp0_0[k]
                  - f_3 * pc_z[k] * fpp1_0[k];

        t_19[k] = f_2 * pc_y[k] * gps_6[k];

        t_20[k] = f_1 * fps_0[k]
                  + f_2 * pc_z[k] * gps_6[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, t_26, pc_x, pc_y, pc_z, fps_1, fps_2, \
                         fps_7, fps_8, gss_2, gps_7, gps_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_4 * fps_7[k]
                  + f_2 * pc_x[k] * gps_7[k];

        t_22[k] = f_1 * gss_2[k]
                  + f_2 * pc_y[k] * gps_7[k];

        t_23[k] = f_1 * fps_1[k]
                  + f_2 * pc_z[k] * gps_7[k];

        t_24[k] = f_4 * fps_8[k]
                  + f_2 * pc_x[k] * gps_8[k];

        t_25[k] = f_2 * pc_y[k] * gps_8[k];

        t_26[k] = f_1 * fps_2[k]
                  + f_1 * gss_2[k]
                  + f_2 * pc_z[k] * gps_8[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, t_32, pc_x, pc_y, pc_z, fps_3, fps_4, \
                         fps_9, fps_10, gss_3, gps_9, gps_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_5 * fps_9[k]
                  + f_1 * gss_3[k]
                  + f_2 * pc_x[k] * gps_9[k];

        t_28[k] = f_5 * fps_3[k]
                  + f_2 * pc_y[k] * gps_9[k];

        t_29[k] = f_2 * pc_z[k] * gps_9[k];

        t_30[k] = f_5 * fps_10[k]
                  + f_2 * pc_x[k] * gps_10[k];

        t_31[k] = f_5 * fps_4[k]
                  + f_1 * gss_3[k]
                  + f_2 * pc_y[k] * gps_10[k];

        t_32[k] = f_2 * pc_z[k] * gps_10[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pa_y, pc_x, pc_y, pc_z, fpp0_18, fps_5, \
                         fps_11, fpp1_18, gss_3, gps_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_5 * fps_11[k]
                  + f_2 * pc_x[k] * gps_11[k];

        t_34[k] = f_5 * fps_5[k]
                  + f_2 * pc_y[k] * gps_11[k];

        t_35[k] = f_1 * gss_3[k]
                  + f_2 * pc_z[k] * gps_11[k];

        t_36[k] = pa_y[k] * fpp0_18[k]
                  - f_3 * pc_y[k] * fpp1_18[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pa_z, pc_x, pc_y, pc_z, fpp0_13, fps_3, \
                         fps_6, fps_13, fpp1_13, gps_12, gps_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_1 * fps_6[k]
                  + f_2 * pc_y[k] * gps_12[k];

        t_38[k] = f_1 * fps_3[k]
                  + f_2 * pc_z[k] * gps_12[k];

        t_39[k] = f_5 * fps_13[k]
                  + f_2 * pc_x[k] * gps_13[k];

        t_40[k] = pa_z[k] * fpp0_13[k]
                  - f_3 * pc_z[k] * fpp1_13[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pa_y, pc_x, pc_y, pc_z, fpp0_26, fps_4, \
                         fps_8, fps_14, fpp1_26, gps_13, gps_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_1 * fps_4[k]
                  + f_2 * pc_z[k] * gps_13[k];

        t_42[k] = f_5 * fps_14[k]
                  + f_2 * pc_x[k] * gps_14[k];

        t_43[k] = f_1 * fps_8[k]
                  + f_2 * pc_y[k] * gps_14[k];

        t_44[k] = pa_y[k] * fpp0_26[k]
                  - f_3 * pc_y[k] * fpp1_26[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, t_50, pc_x, pc_y, pc_z, fps_6, fps_7, \
                         fps_15, fps_16, gss_5, gps_15, gps_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_5 * fps_15[k]
                  + f_1 * gss_5[k]
                  + f_2 * pc_x[k] * gps_15[k];

        t_46[k] = f_2 * pc_y[k] * gps_15[k];

        t_47[k] = f_5 * fps_6[k]
                  + f_2 * pc_z[k] * gps_15[k];

        t_48[k] = f_5 * fps_16[k]
                  + f_2 * pc_x[k] * gps_16[k];

        t_49[k] = f_1 * gss_5[k]
                  + f_2 * pc_y[k] * gps_16[k];

        t_50[k] = f_5 * fps_7[k]
                  + f_2 * pc_z[k] * gps_16[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, t_55, pc_x, pc_y, pc_z, fps_8, fps_9, fps_17, \
                         fps_18, gss_5, gss_6, gps_17, gps_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_5 * fps_17[k]
                  + f_2 * pc_x[k] * gps_17[k];

        t_52[k] = f_2 * pc_y[k] * gps_17[k];

        t_53[k] = f_5 * fps_8[k]
                  + f_1 * gss_5[k]
                  + f_2 * pc_z[k] * gps_17[k];

        t_54[k] = f_1 * fps_18[k]
                  + f_1 * gss_6[k]
                  + f_2 * pc_x[k] * gps_18[k];

        t_55[k] = f_4 * fps_9[k]
                  + f_2 * pc_y[k] * gps_18[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, pa_x, pc_x, pc_z, fpp0_58, fps_19, \
                         fps_20, fpp1_58, gps_18, gps_19, gps_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_2 * pc_z[k] * gps_18[k];

        t_57[k] = f_1 * fps_19[k]
                  + f_2 * pc_x[k] * gps_19[k];

        t_58[k] = pa_x[k] * fpp0_58[k]
                  - f_3 * pc_x[k] * fpp1_58[k];

        t_59[k] = f_2 * pc_z[k] * gps_19[k];

        t_60[k] = f_1 * fps_20[k]
                  + f_2 * pc_x[k] * gps_20[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, t_65, pa_z, pc_y, pc_z, fpp0_27, fps_9, \
                         fps_11, fps_12, fpp1_27, gss_6, gps_20, \
                         gps_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_4 * fps_11[k]
                  + f_2 * pc_y[k] * gps_20[k];

        t_62[k] = f_1 * gss_6[k]
                  + f_2 * pc_z[k] * gps_20[k];

        t_63[k] = pa_z[k] * fpp0_27[k]
                  - f_3 * pc_z[k] * fpp1_27[k];

        t_64[k] = f_5 * fps_12[k]
                  + f_2 * pc_y[k] * gps_21[k];

        t_65[k] = f_1 * fps_9[k]
                  + f_2 * pc_z[k] * gps_21[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pa_x, pc_x, pc_z, fpp0_67, fps_10, fps_22, \
                         fps_23, fpp1_67, gps_22, gps_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_1 * fps_22[k]
                  + f_2 * pc_x[k] * gps_22[k];

        t_67[k] = pa_x[k] * fpp0_67[k]
                  - f_3 * pc_x[k] * fpp1_67[k];

        t_68[k] = f_1 * fps_10[k]
                  + f_2 * pc_z[k] * gps_22[k];

        t_69[k] = f_1 * fps_23[k]
                  + f_2 * pc_x[k] * gps_23[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pa_x, pa_y, pc_x, pc_y, fpp0_45, fpp0_71, \
                         fps_14, fps_15, fpp1_45, fpp1_71, gps_23, \
                         gps_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_5 * fps_14[k]
                  + f_2 * pc_y[k] * gps_23[k];

        t_71[k] = pa_x[k] * fpp0_71[k]
                  - f_3 * pc_x[k] * fpp1_71[k];

        t_72[k] = pa_y[k] * fpp0_45[k]
                  - f_3 * pc_y[k] * fpp1_45[k];

        t_73[k] = f_1 * fps_15[k]
                  + f_2 * pc_y[k] * gps_24[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pa_x, pc_x, pc_z, fpp0_76, fps_12, fps_13, \
                         fps_25, fpp1_76, gps_24, gps_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_5 * fps_12[k]
                  + f_2 * pc_z[k] * gps_24[k];

        t_75[k] = f_1 * fps_25[k]
                  + f_2 * pc_x[k] * gps_25[k];

        t_76[k] = pa_x[k] * fpp0_76[k]
                  - f_3 * pc_x[k] * fpp1_76[k];

        t_77[k] = f_5 * fps_13[k]
                  + f_2 * pc_z[k] * gps_25[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, pa_x, pc_x, pc_y, fpp0_80, fps_17, \
                         fps_26, fps_27, fpp1_80, gss_9, gps_26, \
                         gps_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_1 * fps_26[k]
                  + f_2 * pc_x[k] * gps_26[k];

        t_79[k] = f_1 * fps_17[k]
                  + f_2 * pc_y[k] * gps_26[k];

        t_80[k] = pa_x[k] * fpp0_80[k]
                  - f_3 * pc_x[k] * fpp1_80[k];

        t_81[k] = f_1 * fps_27[k]
                  + f_1 * gss_9[k]
                  + f_2 * pc_x[k] * gps_27[k];

        t_82[k] = f_2 * pc_y[k] * gps_27[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, t_87, pc_x, pc_y, pc_z, fps_15, fps_16, \
                         fps_28, fps_29, gss_9, gps_27, gps_28, \
                         gps_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_4 * fps_15[k]
                  + f_2 * pc_z[k] * gps_27[k];

        t_84[k] = f_1 * fps_28[k]
                  + f_2 * pc_x[k] * gps_28[k];

        t_85[k] = f_1 * gss_9[k]
                  + f_2 * pc_y[k] * gps_28[k];

        t_86[k] = f_4 * fps_16[k]
                  + f_2 * pc_z[k] * gps_28[k];

        t_87[k] = f_1 * fps_29[k]
                  + f_2 * pc_x[k] * gps_29[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, t_92, pa_x, pc_x, pc_y, pc_z, fpp0_89, \
                         fps_18, fpp1_89, gss_10, gps_29, gps_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_2 * pc_y[k] * gps_29[k];

        t_89[k] = pa_x[k] * fpp0_89[k]
                  - f_3 * pc_x[k] * fpp1_89[k];

        t_90[k] = f_1 * gss_10[k]
                  + f_2 * pc_x[k] * gps_30[k];

        t_91[k] = f_0 * fps_18[k]
                  + f_2 * pc_y[k] * gps_30[k];

        t_92[k] = f_2 * pc_z[k] * gps_30[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, t_97, t_98, pc_x, pc_y, pc_z, fps_19, fps_20, \
                         gss_10, gps_31, gps_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_2 * pc_x[k] * gps_31[k];

        t_94[k] = f_0 * fps_19[k]
                  + f_1 * gss_10[k]
                  + f_2 * pc_y[k] * gps_31[k];

        t_95[k] = f_2 * pc_z[k] * gps_31[k];

        t_96[k] = f_2 * pc_x[k] * gps_32[k];

        t_97[k] = f_0 * fps_20[k]
                  + f_2 * pc_y[k] * gps_32[k];

        t_98[k] = f_1 * gss_10[k]
                  + f_2 * pc_z[k] * gps_32[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, pa_z, pc_x, pc_y, pc_z, fpp0_58, \
                         fps_18, fps_21, fpp1_58, gss_11, gps_33, \
                         gps_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_1 * gss_11[k]
                  + f_2 * pc_x[k] * gps_33[k];

        t_100[k] = f_4 * fps_21[k]
                   + f_2 * pc_y[k] * gps_33[k];

        t_101[k] = f_1 * fps_18[k]
                   + f_2 * pc_z[k] * gps_33[k];

        t_102[k] = f_2 * pc_x[k] * gps_34[k];

        t_103[k] = pa_z[k] * fpp0_58[k]
                   - f_3 * pc_z[k] * fpp1_58[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, t_108, pc_x, pc_y, pc_z, fps_19, fps_20, \
                         fps_23, gss_11, gss_12, gps_34, gps_35, \
                         gps_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_1 * fps_19[k]
                   + f_2 * pc_z[k] * gps_34[k];

        t_105[k] = f_2 * pc_x[k] * gps_35[k];

        t_106[k] = f_4 * fps_23[k]
                   + f_2 * pc_y[k] * gps_35[k];

        t_107[k] = f_1 * fps_20[k]
                   + f_1 * gss_11[k]
                   + f_2 * pc_z[k] * gps_35[k];

        t_108[k] = f_1 * gss_12[k]
                   + f_2 * pc_x[k] * gps_36[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, pc_x, pc_y, pc_z, fps_21, fps_22, \
                         fps_24, fps_25, gss_12, gps_36, gps_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_5 * fps_24[k]
                   + f_2 * pc_y[k] * gps_36[k];

        t_110[k] = f_5 * fps_21[k]
                   + f_2 * pc_z[k] * gps_36[k];

        t_111[k] = f_2 * pc_x[k] * gps_37[k];

        t_112[k] = f_5 * fps_25[k]
                   + f_1 * gss_12[k]
                   + f_2 * pc_y[k] * gps_37[k];

        t_113[k] = f_5 * fps_22[k]
                   + f_2 * pc_z[k] * gps_37[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, pc_x, pc_y, pc_z, fps_23, fps_26, \
                         fps_27, gss_12, gss_13, gps_38, gps_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_2 * pc_x[k] * gps_38[k];

        t_115[k] = f_5 * fps_26[k]
                   + f_2 * pc_y[k] * gps_38[k];

        t_116[k] = f_5 * fps_23[k]
                   + f_1 * gss_12[k]
                   + f_2 * pc_z[k] * gps_38[k];

        t_117[k] = f_1 * gss_13[k]
                   + f_2 * pc_x[k] * gps_39[k];

        t_118[k] = f_1 * fps_27[k]
                   + f_2 * pc_y[k] * gps_39[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, t_123, pc_x, pc_y, pc_z, fps_24, fps_25, \
                         fps_28, gss_13, gps_39, gps_40, gps_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = f_4 * fps_24[k]
                   + f_2 * pc_z[k] * gps_39[k];

        t_120[k] = f_2 * pc_x[k] * gps_40[k];

        t_121[k] = f_1 * fps_28[k]
                   + f_1 * gss_13[k]
                   + f_2 * pc_y[k] * gps_40[k];

        t_122[k] = f_4 * fps_25[k]
                   + f_2 * pc_z[k] * gps_40[k];

        t_123[k] = f_2 * pc_x[k] * gps_41[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, t_128, pa_y, pc_x, pc_y, pc_z, fpp0_89, \
                         fps_27, fps_29, fpp1_89, gss_14, gps_41, \
                         gps_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_1 * fps_29[k]
                   + f_2 * pc_y[k] * gps_41[k];

        t_125[k] = pa_y[k] * fpp0_89[k]
                   - f_3 * pc_y[k] * fpp1_89[k];

        t_126[k] = f_1 * gss_14[k]
                   + f_2 * pc_x[k] * gps_42[k];

        t_127[k] = f_2 * pc_y[k] * gps_42[k];

        t_128[k] = f_0 * fps_27[k]
                   + f_2 * pc_z[k] * gps_42[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, t_133, t_134, pc_x, pc_y, pc_z, fps_28, \
                         fps_29, gss_14, gps_43, gps_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_2 * pc_x[k] * gps_43[k];

        t_130[k] = f_1 * gss_14[k]
                   + f_2 * pc_y[k] * gps_43[k];

        t_131[k] = f_0 * fps_28[k]
                   + f_2 * pc_z[k] * gps_43[k];

        t_132[k] = f_2 * pc_x[k] * gps_44[k];

        t_133[k] = f_2 * pc_y[k] * gps_44[k];

        t_134[k] = f_0 * fps_29[k]
                   + f_1 * gss_14[k]
                   + f_2 * pc_z[k] * gps_44[k];
    }
}

}  // namespace simdt3ceri
