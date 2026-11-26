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


#include "SimdThreeCenterElectronRepulsionVrrRecGSD.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_gsd_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t fsd0, const size_t fsp,
                                                   const size_t fsd1, const size_t gss0,
                                                   const size_t gss1, const size_t gsp,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 0.5 / gamma;
    const auto f_2 = 0.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 1.5 / q;
    const auto f_6 = 0.5 / q;
    const auto f_7 = 1.0 / q;

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

    const auto *fsd0_0 = buffer.data(fsd0 + 0);
    const auto *fsd0_3 = buffer.data(fsd0 + 3);
    const auto *fsd0_5 = buffer.data(fsd0 + 5);
    const auto *fsd0_9 = buffer.data(fsd0 + 9);
    const auto *fsd0_12 = buffer.data(fsd0 + 12);
    const auto *fsd0_17 = buffer.data(fsd0 + 17);
    const auto *fsd0_18 = buffer.data(fsd0 + 18);
    const auto *fsd0_30 = buffer.data(fsd0 + 30);
    const auto *fsd0_36 = buffer.data(fsd0 + 36);
    const auto *fsd0_39 = buffer.data(fsd0 + 39);
    const auto *fsd0_41 = buffer.data(fsd0 + 41);
    const auto *fsd0_45 = buffer.data(fsd0 + 45);
    const auto *fsd0_47 = buffer.data(fsd0 + 47);
    const auto *fsd0_51 = buffer.data(fsd0 + 51);
    const auto *fsd0_53 = buffer.data(fsd0 + 53);
    const auto *fsd0_54 = buffer.data(fsd0 + 54);
    const auto *fsd0_57 = buffer.data(fsd0 + 57);
    const auto *fsd0_59 = buffer.data(fsd0 + 59);

    const auto *fsp_0 = buffer.data(fsp + 0);
    const auto *fsp_1 = buffer.data(fsp + 1);
    const auto *fsp_2 = buffer.data(fsp + 2);
    const auto *fsp_4 = buffer.data(fsp + 4);
    const auto *fsp_8 = buffer.data(fsp + 8);
    const auto *fsp_9 = buffer.data(fsp + 9);
    const auto *fsp_10 = buffer.data(fsp + 10);
    const auto *fsp_13 = buffer.data(fsp + 13);
    const auto *fsp_14 = buffer.data(fsp + 14);
    const auto *fsp_15 = buffer.data(fsp + 15);
    const auto *fsp_17 = buffer.data(fsp + 17);
    const auto *fsp_18 = buffer.data(fsp + 18);
    const auto *fsp_19 = buffer.data(fsp + 19);
    const auto *fsp_20 = buffer.data(fsp + 20);
    const auto *fsp_22 = buffer.data(fsp + 22);
    const auto *fsp_23 = buffer.data(fsp + 23);
    const auto *fsp_25 = buffer.data(fsp + 25);
    const auto *fsp_26 = buffer.data(fsp + 26);
    const auto *fsp_27 = buffer.data(fsp + 27);
    const auto *fsp_28 = buffer.data(fsp + 28);
    const auto *fsp_29 = buffer.data(fsp + 29);

    const auto *fsd1_0 = buffer.data(fsd1 + 0);
    const auto *fsd1_3 = buffer.data(fsd1 + 3);
    const auto *fsd1_5 = buffer.data(fsd1 + 5);
    const auto *fsd1_9 = buffer.data(fsd1 + 9);
    const auto *fsd1_12 = buffer.data(fsd1 + 12);
    const auto *fsd1_17 = buffer.data(fsd1 + 17);
    const auto *fsd1_18 = buffer.data(fsd1 + 18);
    const auto *fsd1_30 = buffer.data(fsd1 + 30);
    const auto *fsd1_36 = buffer.data(fsd1 + 36);
    const auto *fsd1_39 = buffer.data(fsd1 + 39);
    const auto *fsd1_41 = buffer.data(fsd1 + 41);
    const auto *fsd1_45 = buffer.data(fsd1 + 45);
    const auto *fsd1_47 = buffer.data(fsd1 + 47);
    const auto *fsd1_51 = buffer.data(fsd1 + 51);
    const auto *fsd1_53 = buffer.data(fsd1 + 53);
    const auto *fsd1_54 = buffer.data(fsd1 + 54);
    const auto *fsd1_57 = buffer.data(fsd1 + 57);
    const auto *fsd1_59 = buffer.data(fsd1 + 59);

    const auto *gss0_0 = buffer.data(gss0 + 0);
    const auto *gss0_1 = buffer.data(gss0 + 1);
    const auto *gss0_2 = buffer.data(gss0 + 2);
    const auto *gss0_3 = buffer.data(gss0 + 3);
    const auto *gss0_5 = buffer.data(gss0 + 5);
    const auto *gss0_10 = buffer.data(gss0 + 10);
    const auto *gss0_11 = buffer.data(gss0 + 11);
    const auto *gss0_12 = buffer.data(gss0 + 12);
    const auto *gss0_14 = buffer.data(gss0 + 14);

    const auto *gss1_0 = buffer.data(gss1 + 0);
    const auto *gss1_1 = buffer.data(gss1 + 1);
    const auto *gss1_2 = buffer.data(gss1 + 2);
    const auto *gss1_3 = buffer.data(gss1 + 3);
    const auto *gss1_5 = buffer.data(gss1 + 5);
    const auto *gss1_10 = buffer.data(gss1 + 10);
    const auto *gss1_11 = buffer.data(gss1 + 11);
    const auto *gss1_12 = buffer.data(gss1 + 12);
    const auto *gss1_14 = buffer.data(gss1 + 14);

    const auto *gsp_0 = buffer.data(gsp + 0);
    const auto *gsp_1 = buffer.data(gsp + 1);
    const auto *gsp_2 = buffer.data(gsp + 2);
    const auto *gsp_3 = buffer.data(gsp + 3);
    const auto *gsp_4 = buffer.data(gsp + 4);
    const auto *gsp_6 = buffer.data(gsp + 6);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, fsp_0, gss0_0, \
                         gss1_0, gsp_0, gsp_1, gsp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fsp_0[k]
                 + f_1 * gss0_0[k]
                 - f_2 * gss1_0[k]
                 + f_3 * pc_x[k] * gsp_0[k];

        t_1[k] = f_3 * pc_y[k] * gsp_0[k];

        t_2[k] = f_3 * pc_z[k] * gsp_0[k];

        t_3[k] = f_1 * gss0_0[k]
                 - f_2 * gss1_0[k]
                 + f_3 * pc_y[k] * gsp_1[k];

        t_4[k] = f_3 * pc_y[k] * gsp_2[k];

        t_5[k] = f_1 * gss0_0[k]
                 - f_2 * gss1_0[k]
                 + f_3 * pc_z[k] * gsp_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pa_y, pc_x, pc_y, pc_z, fsd0_0, fsp_1, fsp_4, \
                         fsd1_0, gss0_1, gss1_1, gsp_3, gsp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pa_y[k] * fsd0_0[k]
                 - f_4 * pc_y[k] * fsd1_0[k];

        t_7[k] = f_5 * fsp_4[k]
                 + f_3 * pc_x[k] * gsp_4[k];

        t_8[k] = f_3 * pc_z[k] * gsp_3[k];

        t_9[k] = f_6 * fsp_1[k]
                 + f_1 * gss0_1[k]
                 - f_2 * gss1_1[k]
                 + f_3 * pc_y[k] * gsp_4[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_y, pa_z, pc_y, pc_z, fsd0_0, fsd0_5, \
                         fsd1_0, fsd1_5, gsp_4, gsp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * pc_z[k] * gsp_4[k];

        t_11[k] = pa_y[k] * fsd0_5[k]
                  - f_4 * pc_y[k] * fsd1_5[k];

        t_12[k] = pa_z[k] * fsd0_0[k]
                  - f_4 * pc_z[k] * fsd1_0[k];

        t_13[k] = f_3 * pc_y[k] * gsp_6[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_z, pc_x, pc_y, pc_z, fsd0_3, fsp_2, fsp_8, \
                         fsd1_3, gss0_2, gss1_2, gsp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_5 * fsp_8[k]
                  + f_3 * pc_x[k] * gsp_8[k];

        t_15[k] = pa_z[k] * fsd0_3[k]
                  - f_4 * pc_z[k] * fsd1_3[k];

        t_16[k] = f_3 * pc_y[k] * gsp_8[k];

        t_17[k] = f_6 * fsp_2[k]
                  + f_1 * gss0_2[k]
                  - f_2 * gss1_2[k]
                  + f_3 * pc_z[k] * gsp_8[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, pc_x, pc_y, pc_z, fsp_4, fsp_9, fsp_10, \
                         gss0_3, gss1_3, gsp_9, gsp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_7 * fsp_9[k]
                  + f_1 * gss0_3[k]
                  - f_2 * gss1_3[k]
                  + f_3 * pc_x[k] * gsp_9[k];

        t_19[k] = f_7 * fsp_10[k]
                  + f_3 * pc_x[k] * gsp_10[k];

        t_20[k] = f_3 * pc_z[k] * gsp_9[k];

        t_21[k] = f_7 * fsp_4[k]
                  + f_1 * gss0_3[k]
                  - f_2 * gss1_3[k]
                  + f_3 * pc_y[k] * gsp_10[k];

        t_22[k] = f_3 * pc_z[k] * gsp_10[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_y, pc_x, pc_y, pc_z, fsd0_12, fsp_13, fsd1_12, \
                         gss0_3, gss1_3, gsp_11, gsp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_1 * gss0_3[k]
                  - f_2 * gss1_3[k]
                  + f_3 * pc_z[k] * gsp_11[k];

        t_24[k] = pa_y[k] * fsd0_12[k]
                  - f_4 * pc_y[k] * fsd1_12[k];

        t_25[k] = f_7 * fsp_13[k]
                  + f_3 * pc_x[k] * gsp_13[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_y, pa_z, pc_x, pc_y, pc_z, fsd0_9, \
                         fsd0_17, fsp_8, fsp_14, fsd1_9, fsd1_17, \
                         gsp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_7 * fsp_14[k]
                  + f_3 * pc_x[k] * gsp_14[k];

        t_27[k] = pa_z[k] * fsd0_9[k]
                  - f_4 * pc_z[k] * fsd1_9[k];

        t_28[k] = f_6 * fsp_8[k]
                  + f_3 * pc_y[k] * gsp_14[k];

        t_29[k] = pa_y[k] * fsd0_17[k]
                  - f_4 * pc_y[k] * fsd1_17[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pc_x, pc_y, fsp_15, fsp_17, gss0_5, \
                         gss1_5, gsp_15, gsp_16, gsp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_7 * fsp_15[k]
                  + f_1 * gss0_5[k]
                  - f_2 * gss1_5[k]
                  + f_3 * pc_x[k] * gsp_15[k];

        t_31[k] = f_3 * pc_y[k] * gsp_15[k];

        t_32[k] = f_7 * fsp_17[k]
                  + f_3 * pc_x[k] * gsp_17[k];

        t_33[k] = f_1 * gss0_5[k]
                  - f_2 * gss1_5[k]
                  + f_3 * pc_y[k] * gsp_16[k];

        t_34[k] = f_3 * pc_y[k] * gsp_17[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pa_x, pc_x, pc_z, fsd0_36, fsp_8, fsp_18, fsp_19, \
                         fsd1_36, gss0_5, gss1_5, gsp_17, gsp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_7 * fsp_8[k]
                  + f_1 * gss0_5[k]
                  - f_2 * gss1_5[k]
                  + f_3 * pc_z[k] * gsp_17[k];

        t_36[k] = pa_x[k] * fsd0_36[k]
                  + f_7 * fsp_18[k]
                  - f_4 * pc_x[k] * fsd1_36[k];

        t_37[k] = f_6 * fsp_19[k]
                  + f_3 * pc_x[k] * gsp_19[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_x, pc_x, pc_z, fsd0_39, fsd0_41, fsd1_39, \
                         fsd1_41, gsp_18, gsp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_3 * pc_z[k] * gsp_18[k];

        t_39[k] = pa_x[k] * fsd0_39[k]
                  - f_4 * pc_x[k] * fsd1_39[k];

        t_40[k] = f_3 * pc_z[k] * gsp_19[k];

        t_41[k] = pa_x[k] * fsd0_41[k]
                  - f_4 * pc_x[k] * fsd1_41[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_x, pa_z, pc_x, pc_z, fsd0_18, fsd0_45, \
                         fsp_22, fsp_23, fsd1_18, fsd1_45, gsp_22, \
                         gsp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pa_z[k] * fsd0_18[k]
                  - f_4 * pc_z[k] * fsd1_18[k];

        t_43[k] = f_6 * fsp_22[k]
                  + f_3 * pc_x[k] * gsp_22[k];

        t_44[k] = f_6 * fsp_23[k]
                  + f_3 * pc_x[k] * gsp_23[k];

        t_45[k] = pa_x[k] * fsd0_45[k]
                  - f_4 * pc_x[k] * fsd1_45[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_x, pa_y, pc_x, pc_y, fsd0_30, fsd0_47, \
                         fsp_14, fsp_25, fsd1_30, fsd1_47, gsp_23, \
                         gsp_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_7 * fsp_14[k]
                  + f_3 * pc_y[k] * gsp_23[k];

        t_47[k] = pa_x[k] * fsd0_47[k]
                  - f_4 * pc_x[k] * fsd1_47[k];

        t_48[k] = pa_y[k] * fsd0_30[k]
                  - f_4 * pc_y[k] * fsd1_30[k];

        t_49[k] = f_6 * fsp_25[k]
                  + f_3 * pc_x[k] * gsp_25[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_x, pc_x, pc_y, fsd0_51, fsd0_53, fsp_17, \
                         fsp_26, fsd1_51, fsd1_53, gsp_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_6 * fsp_26[k]
                  + f_3 * pc_x[k] * gsp_26[k];

        t_51[k] = pa_x[k] * fsd0_51[k]
                  - f_4 * pc_x[k] * fsd1_51[k];

        t_52[k] = f_6 * fsp_17[k]
                  + f_3 * pc_y[k] * gsp_26[k];

        t_53[k] = pa_x[k] * fsd0_53[k]
                  - f_4 * pc_x[k] * fsd1_53[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, pa_x, pc_x, pc_y, fsd0_54, fsd0_57, \
                         fsp_27, fsp_29, fsd1_54, fsd1_57, gsp_27, \
                         gsp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = pa_x[k] * fsd0_54[k]
                  + f_7 * fsp_27[k]
                  - f_4 * pc_x[k] * fsd1_54[k];

        t_55[k] = f_3 * pc_y[k] * gsp_27[k];

        t_56[k] = f_6 * fsp_29[k]
                  + f_3 * pc_x[k] * gsp_29[k];

        t_57[k] = pa_x[k] * fsd0_57[k]
                  - f_4 * pc_x[k] * fsd1_57[k];

        t_58[k] = f_3 * pc_y[k] * gsp_29[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, t_63, pa_x, pc_x, pc_y, fsd0_59, fsp_19, \
                         fsd1_59, gss0_10, gss1_10, gsp_30, gsp_31, \
                         gsp_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = pa_x[k] * fsd0_59[k]
                  - f_4 * pc_x[k] * fsd1_59[k];

        t_60[k] = f_1 * gss0_10[k]
                  - f_2 * gss1_10[k]
                  + f_3 * pc_x[k] * gsp_30[k];

        t_61[k] = f_3 * pc_x[k] * gsp_31[k];

        t_62[k] = f_3 * pc_x[k] * gsp_32[k];

        t_63[k] = f_0 * fsp_19[k]
                  + f_1 * gss0_10[k]
                  - f_2 * gss1_10[k]
                  + f_3 * pc_y[k] * gsp_31[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, t_68, pa_z, pc_x, pc_z, fsd0_36, fsd1_36, \
                         gss0_10, gss1_10, gsp_31, gsp_32, gsp_34, \
                         gsp_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_3 * pc_z[k] * gsp_31[k];

        t_65[k] = f_1 * gss0_10[k]
                  - f_2 * gss1_10[k]
                  + f_3 * pc_z[k] * gsp_32[k];

        t_66[k] = pa_z[k] * fsd0_36[k]
                  - f_4 * pc_z[k] * fsd1_36[k];

        t_67[k] = f_3 * pc_x[k] * gsp_34[k];

        t_68[k] = f_3 * pc_x[k] * gsp_35[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pa_z, pc_y, pc_z, fsd0_39, fsp_20, fsp_23, fsd1_39, \
                         gss0_11, gss1_11, gsp_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = pa_z[k] * fsd0_39[k]
                  - f_4 * pc_z[k] * fsd1_39[k];

        t_70[k] = f_5 * fsp_23[k]
                  + f_3 * pc_y[k] * gsp_35[k];

        t_71[k] = f_6 * fsp_20[k]
                  + f_1 * gss0_11[k]
                  - f_2 * gss1_11[k]
                  + f_3 * pc_z[k] * gsp_35[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, pc_x, pc_y, fsp_25, fsp_26, gss0_12, \
                         gss1_12, gsp_36, gsp_37, gsp_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_1 * gss0_12[k]
                  - f_2 * gss1_12[k]
                  + f_3 * pc_x[k] * gsp_36[k];

        t_73[k] = f_3 * pc_x[k] * gsp_37[k];

        t_74[k] = f_3 * pc_x[k] * gsp_38[k];

        t_75[k] = f_7 * fsp_25[k]
                  + f_1 * gss0_12[k]
                  - f_2 * gss1_12[k]
                  + f_3 * pc_y[k] * gsp_37[k];

        t_76[k] = f_7 * fsp_26[k]
                  + f_3 * pc_y[k] * gsp_38[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, pa_y, pc_x, pc_y, pc_z, fsd0_54, fsp_23, \
                         fsd1_54, gss0_12, gss1_12, gsp_38, gsp_40, \
                         gsp_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_7 * fsp_23[k]
                  + f_1 * gss0_12[k]
                  - f_2 * gss1_12[k]
                  + f_3 * pc_z[k] * gsp_38[k];

        t_78[k] = pa_y[k] * fsd0_54[k]
                  - f_4 * pc_y[k] * fsd1_54[k];

        t_79[k] = f_3 * pc_x[k] * gsp_40[k];

        t_80[k] = f_3 * pc_x[k] * gsp_41[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, pa_y, pc_y, fsd0_57, fsd0_59, fsp_28, fsp_29, \
                         fsd1_57, fsd1_59, gsp_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = pa_y[k] * fsd0_57[k]
                  + f_7 * fsp_28[k]
                  - f_4 * pc_y[k] * fsd1_57[k];

        t_82[k] = f_6 * fsp_29[k]
                  + f_3 * pc_y[k] * gsp_41[k];

        t_83[k] = pa_y[k] * fsd0_59[k]
                  - f_4 * pc_y[k] * fsd1_59[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, t_89, pc_x, pc_y, pc_z, fsp_29, \
                         gss0_14, gss1_14, gsp_42, gsp_43, gsp_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_1 * gss0_14[k]
                  - f_2 * gss1_14[k]
                  + f_3 * pc_x[k] * gsp_42[k];

        t_85[k] = f_3 * pc_x[k] * gsp_43[k];

        t_86[k] = f_3 * pc_x[k] * gsp_44[k];

        t_87[k] = f_1 * gss0_14[k]
                  - f_2 * gss1_14[k]
                  + f_3 * pc_y[k] * gsp_43[k];

        t_88[k] = f_3 * pc_y[k] * gsp_44[k];

        t_89[k] = f_0 * fsp_29[k]
                  + f_1 * gss0_14[k]
                  - f_2 * gss1_14[k]
                  + f_3 * pc_z[k] * gsp_44[k];
    }
}

}  // namespace simdt3ceri
