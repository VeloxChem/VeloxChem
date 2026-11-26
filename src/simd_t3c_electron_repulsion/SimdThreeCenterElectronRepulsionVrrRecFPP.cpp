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


#include "SimdThreeCenterElectronRepulsionVrrRecFPP.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_fpp_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t dpp0, const size_t dps,
                                                   const size_t dpp1, const size_t fss,
                                                   const size_t fps, const size_t ncols,
                                                   const double gamma, const double p,
                                                   const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
    const auto f_1 = 0.5 / q;
    const auto f_2 = p / q;
    const auto f_3 = gamma / q;
    const auto f_4 = 1.0 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *dpp0_0 = buffer.data(dpp0 + 0);
    const auto *dpp0_18 = buffer.data(dpp0 + 18);
    const auto *dpp0_31 = buffer.data(dpp0 + 31);
    const auto *dpp0_40 = buffer.data(dpp0 + 40);
    const auto *dpp0_44 = buffer.data(dpp0 + 44);
    const auto *dpp0_53 = buffer.data(dpp0 + 53);

    const auto *dps_0 = buffer.data(dps + 0);
    const auto *dps_1 = buffer.data(dps + 1);
    const auto *dps_2 = buffer.data(dps + 2);
    const auto *dps_3 = buffer.data(dps + 3);
    const auto *dps_4 = buffer.data(dps + 4);
    const auto *dps_5 = buffer.data(dps + 5);
    const auto *dps_6 = buffer.data(dps + 6);
    const auto *dps_7 = buffer.data(dps + 7);
    const auto *dps_8 = buffer.data(dps + 8);
    const auto *dps_9 = buffer.data(dps + 9);
    const auto *dps_10 = buffer.data(dps + 10);
    const auto *dps_11 = buffer.data(dps + 11);
    const auto *dps_12 = buffer.data(dps + 12);
    const auto *dps_13 = buffer.data(dps + 13);
    const auto *dps_14 = buffer.data(dps + 14);
    const auto *dps_15 = buffer.data(dps + 15);
    const auto *dps_16 = buffer.data(dps + 16);
    const auto *dps_17 = buffer.data(dps + 17);

    const auto *dpp1_0 = buffer.data(dpp1 + 0);
    const auto *dpp1_18 = buffer.data(dpp1 + 18);
    const auto *dpp1_31 = buffer.data(dpp1 + 31);
    const auto *dpp1_40 = buffer.data(dpp1 + 40);
    const auto *dpp1_44 = buffer.data(dpp1 + 44);
    const auto *dpp1_53 = buffer.data(dpp1 + 53);

    const auto *fss_0 = buffer.data(fss + 0);
    const auto *fss_1 = buffer.data(fss + 1);
    const auto *fss_2 = buffer.data(fss + 2);
    const auto *fss_3 = buffer.data(fss + 3);
    const auto *fss_5 = buffer.data(fss + 5);
    const auto *fss_6 = buffer.data(fss + 6);
    const auto *fss_7 = buffer.data(fss + 7);
    const auto *fss_8 = buffer.data(fss + 8);
    const auto *fss_9 = buffer.data(fss + 9);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, dps_0, dps_1, fss_0, \
                         fps_0, fps_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dps_0[k]
                 + f_1 * fss_0[k]
                 + f_2 * pc_x[k] * fps_0[k];

        t_1[k] = f_2 * pc_y[k] * fps_0[k];

        t_2[k] = f_2 * pc_z[k] * fps_0[k];

        t_3[k] = f_0 * dps_1[k]
                 + f_2 * pc_x[k] * fps_1[k];

        t_4[k] = f_1 * fss_0[k]
                 + f_2 * pc_y[k] * fps_1[k];

        t_5[k] = f_2 * pc_z[k] * fps_1[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pa_y, pc_x, pc_y, pc_z, dpp0_0, dps_0, \
                         dps_2, dpp1_0, fss_0, fps_2, fps_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * dps_2[k]
                 + f_2 * pc_x[k] * fps_2[k];

        t_7[k] = f_2 * pc_y[k] * fps_2[k];

        t_8[k] = f_1 * fss_0[k]
                 + f_2 * pc_z[k] * fps_2[k];

        t_9[k] = pa_y[k] * dpp0_0[k]
                 - f_3 * pc_y[k] * dpp1_0[k];

        t_10[k] = f_1 * dps_0[k]
                  + f_2 * pc_y[k] * fps_3[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pc_x, pc_y, pc_z, dps_1, dps_4, dps_5, \
                         fss_1, fps_3, fps_4, fps_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_2 * pc_z[k] * fps_3[k];

        t_12[k] = f_4 * dps_4[k]
                  + f_2 * pc_x[k] * fps_4[k];

        t_13[k] = f_1 * dps_1[k]
                  + f_1 * fss_1[k]
                  + f_2 * pc_y[k] * fps_4[k];

        t_14[k] = f_2 * pc_z[k] * fps_4[k];

        t_15[k] = f_4 * dps_5[k]
                  + f_2 * pc_x[k] * fps_5[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pa_z, pc_y, pc_z, dpp0_0, dps_0, dps_2, \
                         dpp1_0, fss_1, fps_5, fps_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_1 * dps_2[k]
                  + f_2 * pc_y[k] * fps_5[k];

        t_17[k] = f_1 * fss_1[k]
                  + f_2 * pc_z[k] * fps_5[k];

        t_18[k] = pa_z[k] * dpp0_0[k]
                  - f_3 * pc_z[k] * dpp1_0[k];

        t_19[k] = f_2 * pc_y[k] * fps_6[k];

        t_20[k] = f_1 * dps_0[k]
                  + f_2 * pc_z[k] * fps_6[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, t_26, pc_x, pc_y, pc_z, dps_1, dps_2, \
                         dps_7, dps_8, fss_2, fps_7, fps_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_4 * dps_7[k]
                  + f_2 * pc_x[k] * fps_7[k];

        t_22[k] = f_1 * fss_2[k]
                  + f_2 * pc_y[k] * fps_7[k];

        t_23[k] = f_1 * dps_1[k]
                  + f_2 * pc_z[k] * fps_7[k];

        t_24[k] = f_4 * dps_8[k]
                  + f_2 * pc_x[k] * fps_8[k];

        t_25[k] = f_2 * pc_y[k] * fps_8[k];

        t_26[k] = f_1 * dps_2[k]
                  + f_1 * fss_2[k]
                  + f_2 * pc_z[k] * fps_8[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pc_x, pc_y, pc_z, dps_3, dps_9, dps_10, \
                         fss_3, fps_9, fps_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_1 * dps_9[k]
                  + f_1 * fss_3[k]
                  + f_2 * pc_x[k] * fps_9[k];

        t_28[k] = f_4 * dps_3[k]
                  + f_2 * pc_y[k] * fps_9[k];

        t_29[k] = f_2 * pc_z[k] * fps_9[k];

        t_30[k] = f_1 * dps_10[k]
                  + f_2 * pc_x[k] * fps_10[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, pa_x, pc_x, pc_y, pc_z, dpp0_31, dps_5, \
                         dps_11, dpp1_31, fss_3, fps_10, fps_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = pa_x[k] * dpp0_31[k]
                  - f_3 * pc_x[k] * dpp1_31[k];

        t_32[k] = f_2 * pc_z[k] * fps_10[k];

        t_33[k] = f_1 * dps_11[k]
                  + f_2 * pc_x[k] * fps_11[k];

        t_34[k] = f_4 * dps_5[k]
                  + f_2 * pc_y[k] * fps_11[k];

        t_35[k] = f_1 * fss_3[k]
                  + f_2 * pc_z[k] * fps_11[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_y, pc_x, pc_y, pc_z, dpp0_18, dps_3, \
                         dps_6, dps_13, dpp1_18, fps_12, fps_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = pa_y[k] * dpp0_18[k]
                  - f_3 * pc_y[k] * dpp1_18[k];

        t_37[k] = f_1 * dps_6[k]
                  + f_2 * pc_y[k] * fps_12[k];

        t_38[k] = f_1 * dps_3[k]
                  + f_2 * pc_z[k] * fps_12[k];

        t_39[k] = f_1 * dps_13[k]
                  + f_2 * pc_x[k] * fps_13[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_x, pc_x, pc_y, pc_z, dpp0_40, dps_4, \
                         dps_8, dps_14, dpp1_40, fps_13, fps_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pa_x[k] * dpp0_40[k]
                  - f_3 * pc_x[k] * dpp1_40[k];

        t_41[k] = f_1 * dps_4[k]
                  + f_2 * pc_z[k] * fps_13[k];

        t_42[k] = f_1 * dps_14[k]
                  + f_2 * pc_x[k] * fps_14[k];

        t_43[k] = f_1 * dps_8[k]
                  + f_2 * pc_y[k] * fps_14[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_x, pc_x, pc_y, pc_z, dpp0_44, dps_6, \
                         dps_15, dpp1_44, fss_5, fps_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = pa_x[k] * dpp0_44[k]
                  - f_3 * pc_x[k] * dpp1_44[k];

        t_45[k] = f_1 * dps_15[k]
                  + f_1 * fss_5[k]
                  + f_2 * pc_x[k] * fps_15[k];

        t_46[k] = f_2 * pc_y[k] * fps_15[k];

        t_47[k] = f_4 * dps_6[k]
                  + f_2 * pc_z[k] * fps_15[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, pc_x, pc_y, pc_z, dps_7, dps_16, \
                         dps_17, fss_5, fps_16, fps_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_1 * dps_16[k]
                  + f_2 * pc_x[k] * fps_16[k];

        t_49[k] = f_1 * fss_5[k]
                  + f_2 * pc_y[k] * fps_16[k];

        t_50[k] = f_4 * dps_7[k]
                  + f_2 * pc_z[k] * fps_16[k];

        t_51[k] = f_1 * dps_17[k]
                  + f_2 * pc_x[k] * fps_17[k];

        t_52[k] = f_2 * pc_y[k] * fps_17[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, pa_x, pc_x, pc_y, pc_z, dpp0_53, dps_9, \
                         dpp1_53, fss_6, fps_18, fps_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = pa_x[k] * dpp0_53[k]
                  - f_3 * pc_x[k] * dpp1_53[k];

        t_54[k] = f_1 * fss_6[k]
                  + f_2 * pc_x[k] * fps_18[k];

        t_55[k] = f_0 * dps_9[k]
                  + f_2 * pc_y[k] * fps_18[k];

        t_56[k] = f_2 * pc_z[k] * fps_18[k];

        t_57[k] = f_2 * pc_x[k] * fps_19[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, t_63, pc_x, pc_y, pc_z, dps_10, dps_11, \
                         fss_6, fss_7, fps_19, fps_20, fps_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_0 * dps_10[k]
                  + f_1 * fss_6[k]
                  + f_2 * pc_y[k] * fps_19[k];

        t_59[k] = f_2 * pc_z[k] * fps_19[k];

        t_60[k] = f_2 * pc_x[k] * fps_20[k];

        t_61[k] = f_0 * dps_11[k]
                  + f_2 * pc_y[k] * fps_20[k];

        t_62[k] = f_1 * fss_6[k]
                  + f_2 * pc_z[k] * fps_20[k];

        t_63[k] = f_1 * fss_7[k]
                  + f_2 * pc_x[k] * fps_21[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, t_68, pa_z, pc_x, pc_y, pc_z, dpp0_31, dps_9, \
                         dps_10, dps_12, dpp1_31, fps_21, fps_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_4 * dps_12[k]
                  + f_2 * pc_y[k] * fps_21[k];

        t_65[k] = f_1 * dps_9[k]
                  + f_2 * pc_z[k] * fps_21[k];

        t_66[k] = f_2 * pc_x[k] * fps_22[k];

        t_67[k] = pa_z[k] * dpp0_31[k]
                  - f_3 * pc_z[k] * dpp1_31[k];

        t_68[k] = f_1 * dps_10[k]
                  + f_2 * pc_z[k] * fps_22[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, t_73, pc_x, pc_y, pc_z, dps_11, dps_14, \
                         dps_15, fss_7, fss_8, fps_23, fps_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_2 * pc_x[k] * fps_23[k];

        t_70[k] = f_4 * dps_14[k]
                  + f_2 * pc_y[k] * fps_23[k];

        t_71[k] = f_1 * dps_11[k]
                  + f_1 * fss_7[k]
                  + f_2 * pc_z[k] * fps_23[k];

        t_72[k] = f_1 * fss_8[k]
                  + f_2 * pc_x[k] * fps_24[k];

        t_73[k] = f_1 * dps_15[k]
                  + f_2 * pc_y[k] * fps_24[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, t_78, pc_x, pc_y, pc_z, dps_12, dps_13, \
                         dps_16, fss_8, fps_24, fps_25, fps_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_4 * dps_12[k]
                  + f_2 * pc_z[k] * fps_24[k];

        t_75[k] = f_2 * pc_x[k] * fps_25[k];

        t_76[k] = f_1 * dps_16[k]
                  + f_1 * fss_8[k]
                  + f_2 * pc_y[k] * fps_25[k];

        t_77[k] = f_4 * dps_13[k]
                  + f_2 * pc_z[k] * fps_25[k];

        t_78[k] = f_2 * pc_x[k] * fps_26[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, t_83, pa_y, pc_x, pc_y, pc_z, dpp0_53, \
                         dps_15, dps_17, dpp1_53, fss_9, fps_26, \
                         fps_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_1 * dps_17[k]
                  + f_2 * pc_y[k] * fps_26[k];

        t_80[k] = pa_y[k] * dpp0_53[k]
                  - f_3 * pc_y[k] * dpp1_53[k];

        t_81[k] = f_1 * fss_9[k]
                  + f_2 * pc_x[k] * fps_27[k];

        t_82[k] = f_2 * pc_y[k] * fps_27[k];

        t_83[k] = f_0 * dps_15[k]
                  + f_2 * pc_z[k] * fps_27[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, t_89, pc_x, pc_y, pc_z, dps_16, dps_17, \
                         fss_9, fps_28, fps_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_2 * pc_x[k] * fps_28[k];

        t_85[k] = f_1 * fss_9[k]
                  + f_2 * pc_y[k] * fps_28[k];

        t_86[k] = f_0 * dps_16[k]
                  + f_2 * pc_z[k] * fps_28[k];

        t_87[k] = f_2 * pc_x[k] * fps_29[k];

        t_88[k] = f_2 * pc_y[k] * fps_29[k];

        t_89[k] = f_0 * dps_17[k]
                  + f_1 * fss_9[k]
                  + f_2 * pc_z[k] * fps_29[k];
    }
}

}  // namespace simdt3ceri
