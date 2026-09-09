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


#include "SimdThreeCenterElectronRepulsionVrrRecDSG.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_dsg_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t psg0, const size_t psf,
                                                   const size_t psg1, const size_t dsd0,
                                                   const size_t dsd1, const size_t dsf,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / gamma;
    const auto f_9 = p / (gamma * q);

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

    const auto *psg0_0 = buffer.data(psg0 + 0);
    const auto *psg0_3 = buffer.data(psg0 + 3);
    const auto *psg0_5 = buffer.data(psg0 + 5);
    const auto *psg0_16 = buffer.data(psg0 + 16);
    const auto *psg0_18 = buffer.data(psg0 + 18);
    const auto *psg0_25 = buffer.data(psg0 + 25);
    const auto *psg0_27 = buffer.data(psg0 + 27);
    const auto *psg0_29 = buffer.data(psg0 + 29);
    const auto *psg0_30 = buffer.data(psg0 + 30);
    const auto *psg0_32 = buffer.data(psg0 + 32);
    const auto *psg0_35 = buffer.data(psg0 + 35);
    const auto *psg0_40 = buffer.data(psg0 + 40);
    const auto *psg0_41 = buffer.data(psg0 + 41);
    const auto *psg0_42 = buffer.data(psg0 + 42);
    const auto *psg0_44 = buffer.data(psg0 + 44);

    const auto *psf_0 = buffer.data(psf + 0);
    const auto *psf_6 = buffer.data(psf + 6);
    const auto *psf_9 = buffer.data(psf + 9);
    const auto *psf_13 = buffer.data(psf + 13);
    const auto *psf_16 = buffer.data(psf + 16);
    const auto *psf_18 = buffer.data(psf + 18);
    const auto *psf_19 = buffer.data(psf + 19);
    const auto *psf_25 = buffer.data(psf + 25);
    const auto *psf_26 = buffer.data(psf + 26);
    const auto *psf_27 = buffer.data(psf + 27);
    const auto *psf_28 = buffer.data(psf + 28);
    const auto *psf_29 = buffer.data(psf + 29);

    const auto *psg1_0 = buffer.data(psg1 + 0);
    const auto *psg1_3 = buffer.data(psg1 + 3);
    const auto *psg1_5 = buffer.data(psg1 + 5);
    const auto *psg1_16 = buffer.data(psg1 + 16);
    const auto *psg1_18 = buffer.data(psg1 + 18);
    const auto *psg1_25 = buffer.data(psg1 + 25);
    const auto *psg1_27 = buffer.data(psg1 + 27);
    const auto *psg1_29 = buffer.data(psg1 + 29);
    const auto *psg1_30 = buffer.data(psg1 + 30);
    const auto *psg1_32 = buffer.data(psg1 + 32);
    const auto *psg1_35 = buffer.data(psg1 + 35);
    const auto *psg1_40 = buffer.data(psg1 + 40);
    const auto *psg1_41 = buffer.data(psg1 + 41);
    const auto *psg1_42 = buffer.data(psg1 + 42);
    const auto *psg1_44 = buffer.data(psg1 + 44);

    const auto *dsd0_0 = buffer.data(dsd0 + 0);
    const auto *dsd0_3 = buffer.data(dsd0 + 3);
    const auto *dsd0_5 = buffer.data(dsd0 + 5);
    const auto *dsd0_18 = buffer.data(dsd0 + 18);
    const auto *dsd0_19 = buffer.data(dsd0 + 19);
    const auto *dsd0_21 = buffer.data(dsd0 + 21);
    const auto *dsd0_23 = buffer.data(dsd0 + 23);
    const auto *dsd0_28 = buffer.data(dsd0 + 28);
    const auto *dsd0_30 = buffer.data(dsd0 + 30);
    const auto *dsd0_32 = buffer.data(dsd0 + 32);
    const auto *dsd0_33 = buffer.data(dsd0 + 33);
    const auto *dsd0_34 = buffer.data(dsd0 + 34);
    const auto *dsd0_35 = buffer.data(dsd0 + 35);

    const auto *dsd1_0 = buffer.data(dsd1 + 0);
    const auto *dsd1_3 = buffer.data(dsd1 + 3);
    const auto *dsd1_5 = buffer.data(dsd1 + 5);
    const auto *dsd1_18 = buffer.data(dsd1 + 18);
    const auto *dsd1_19 = buffer.data(dsd1 + 19);
    const auto *dsd1_21 = buffer.data(dsd1 + 21);
    const auto *dsd1_23 = buffer.data(dsd1 + 23);
    const auto *dsd1_28 = buffer.data(dsd1 + 28);
    const auto *dsd1_30 = buffer.data(dsd1 + 30);
    const auto *dsd1_32 = buffer.data(dsd1 + 32);
    const auto *dsd1_33 = buffer.data(dsd1 + 33);
    const auto *dsd1_34 = buffer.data(dsd1 + 34);
    const auto *dsd1_35 = buffer.data(dsd1 + 35);

    const auto *dsf_0 = buffer.data(dsf + 0);
    const auto *dsf_1 = buffer.data(dsf + 1);
    const auto *dsf_2 = buffer.data(dsf + 2);
    const auto *dsf_3 = buffer.data(dsf + 3);
    const auto *dsf_5 = buffer.data(dsf + 5);
    const auto *dsf_6 = buffer.data(dsf + 6);
    const auto *dsf_8 = buffer.data(dsf + 8);
    const auto *dsf_9 = buffer.data(dsf + 9);
    const auto *dsf_10 = buffer.data(dsf + 10);
    const auto *dsf_11 = buffer.data(dsf + 11);
    const auto *dsf_13 = buffer.data(dsf + 13);
    const auto *dsf_16 = buffer.data(dsf + 16);
    const auto *dsf_18 = buffer.data(dsf + 18);
    const auto *dsf_19 = buffer.data(dsf + 19);
    const auto *dsf_20 = buffer.data(dsf + 20);
    const auto *dsf_22 = buffer.data(dsf + 22);
    const auto *dsf_25 = buffer.data(dsf + 25);
    const auto *dsf_26 = buffer.data(dsf + 26);
    const auto *dsf_27 = buffer.data(dsf + 27);
    const auto *dsf_29 = buffer.data(dsf + 29);
    const auto *dsf_30 = buffer.data(dsf + 30);
    const auto *dsf_31 = buffer.data(dsf + 31);
    const auto *dsf_33 = buffer.data(dsf + 33);
    const auto *dsf_35 = buffer.data(dsf + 35);
    const auto *dsf_36 = buffer.data(dsf + 36);
    const auto *dsf_37 = buffer.data(dsf + 37);
    const auto *dsf_38 = buffer.data(dsf + 38);
    const auto *dsf_39 = buffer.data(dsf + 39);
    const auto *dsf_44 = buffer.data(dsf + 44);
    const auto *dsf_46 = buffer.data(dsf + 46);
    const auto *dsf_47 = buffer.data(dsf + 47);
    const auto *dsf_48 = buffer.data(dsf + 48);
    const auto *dsf_49 = buffer.data(dsf + 49);
    const auto *dsf_50 = buffer.data(dsf + 50);
    const auto *dsf_52 = buffer.data(dsf + 52);
    const auto *dsf_53 = buffer.data(dsf + 53);
    const auto *dsf_55 = buffer.data(dsf + 55);
    const auto *dsf_56 = buffer.data(dsf + 56);
    const auto *dsf_57 = buffer.data(dsf + 57);
    const auto *dsf_58 = buffer.data(dsf + 58);
    const auto *dsf_59 = buffer.data(dsf + 59);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, psf_0, dsd0_0, \
                         dsd1_0, dsf_0, dsf_1, dsf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * psf_0[k]
                 + f_1 * dsd0_0[k]
                 - f_2 * dsd1_0[k]
                 + f_3 * pc_x[k] * dsf_0[k];

        t_1[k] = f_3 * pc_y[k] * dsf_0[k];

        t_2[k] = f_3 * pc_z[k] * dsf_0[k];

        t_3[k] = f_4 * dsd0_0[k]
                 - f_5 * dsd1_0[k]
                 + f_3 * pc_y[k] * dsf_1[k];

        t_4[k] = f_3 * pc_y[k] * dsf_2[k];

        t_5[k] = f_4 * dsd0_0[k]
                 - f_5 * dsd1_0[k]
                 + f_3 * pc_z[k] * dsf_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_x, pc_y, pc_z, psf_6, psf_9, dsd0_3, \
                         dsd1_3, dsf_3, dsf_5, dsf_6, dsf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * psf_6[k]
                 + f_3 * pc_x[k] * dsf_6[k];

        t_7[k] = f_3 * pc_z[k] * dsf_3[k];

        t_8[k] = f_3 * pc_y[k] * dsf_5[k];

        t_9[k] = f_0 * psf_9[k]
                 + f_3 * pc_x[k] * dsf_9[k];

        t_10[k] = f_1 * dsd0_3[k]
                  - f_2 * dsd1_3[k]
                  + f_3 * pc_y[k] * dsf_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pa_y, pc_y, pc_z, psg0_0, psg1_0, \
                         dsd0_5, dsd1_5, dsf_6, dsf_8, dsf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * pc_z[k] * dsf_6[k];

        t_12[k] = f_4 * dsd0_5[k]
                  - f_5 * dsd1_5[k]
                  + f_3 * pc_y[k] * dsf_8[k];

        t_13[k] = f_3 * pc_y[k] * dsf_9[k];

        t_14[k] = f_1 * dsd0_5[k]
                  - f_2 * dsd1_5[k]
                  + f_3 * pc_z[k] * dsf_9[k];

        t_15[k] = pa_y[k] * psg0_0[k]
                  - f_6 * pc_y[k] * psg1_0[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pc_x, pc_y, pc_z, psg0_18, psf_0, \
                         psf_13, psg1_18, dsf_10, dsf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_7 * psf_0[k]
                  + f_3 * pc_y[k] * dsf_10[k];

        t_17[k] = f_3 * pc_z[k] * dsf_10[k];

        t_18[k] = pa_x[k] * psg0_18[k]
                  + f_0 * psf_13[k]
                  - f_6 * pc_x[k] * psg1_18[k];

        t_19[k] = f_3 * pc_z[k] * dsf_11[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_y, pc_x, pc_y, pc_z, psg0_5, psf_16, \
                         psf_18, psg1_5, dsf_13, dsf_16, dsf_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_y[k] * psg0_5[k]
                  - f_6 * pc_y[k] * psg1_5[k];

        t_21[k] = f_7 * psf_16[k]
                  + f_3 * pc_x[k] * dsf_16[k];

        t_22[k] = f_3 * pc_z[k] * dsf_13[k];

        t_23[k] = f_7 * psf_18[k]
                  + f_3 * pc_x[k] * dsf_18[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_x, pc_x, pc_z, psg0_25, psg0_27, psf_19, \
                         psg1_25, psg1_27, dsf_16, dsf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_7 * psf_19[k]
                  + f_3 * pc_x[k] * dsf_19[k];

        t_25[k] = pa_x[k] * psg0_25[k]
                  - f_6 * pc_x[k] * psg1_25[k];

        t_26[k] = f_3 * pc_z[k] * dsf_16[k];

        t_27[k] = pa_x[k] * psg0_27[k]
                  - f_6 * pc_x[k] * psg1_27[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_x, pa_z, pc_x, pc_y, pc_z, psg0_0, \
                         psg0_29, psf_9, psg1_0, psg1_29, dsf_19, \
                         dsf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_7 * psf_9[k]
                  + f_3 * pc_y[k] * dsf_19[k];

        t_29[k] = pa_x[k] * psg0_29[k]
                  - f_6 * pc_x[k] * psg1_29[k];

        t_30[k] = pa_z[k] * psg0_0[k]
                  - f_6 * pc_z[k] * psg1_0[k];

        t_31[k] = f_3 * pc_y[k] * dsf_20[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pa_z, pc_y, pc_z, psg0_3, psf_0, psg1_3, dsf_20, \
                         dsf_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_7 * psf_0[k]
                  + f_3 * pc_z[k] * dsf_20[k];

        t_33[k] = pa_z[k] * psg0_3[k]
                  - f_6 * pc_z[k] * psg1_3[k];

        t_34[k] = f_3 * pc_y[k] * dsf_22[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pa_x, pc_x, pc_y, psg0_35, psf_25, psf_26, \
                         psf_27, psg1_35, dsf_25, dsf_26, dsf_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = pa_x[k] * psg0_35[k]
                  + f_0 * psf_25[k]
                  - f_6 * pc_x[k] * psg1_35[k];

        t_36[k] = f_7 * psf_26[k]
                  + f_3 * pc_x[k] * dsf_26[k];

        t_37[k] = f_7 * psf_27[k]
                  + f_3 * pc_x[k] * dsf_27[k];

        t_38[k] = f_3 * pc_y[k] * dsf_25[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, t_43, pa_x, pc_x, pc_y, psg0_40, psg0_41, \
                         psg0_42, psf_29, psg1_40, psg1_41, psg1_42, \
                         dsf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_7 * psf_29[k]
                  + f_3 * pc_x[k] * dsf_29[k];

        t_40[k] = pa_x[k] * psg0_40[k]
                  - f_6 * pc_x[k] * psg1_40[k];

        t_41[k] = pa_x[k] * psg0_41[k]
                  - f_6 * pc_x[k] * psg1_41[k];

        t_42[k] = pa_x[k] * psg0_42[k]
                  - f_6 * pc_x[k] * psg1_42[k];

        t_43[k] = f_3 * pc_y[k] * dsf_29[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_x, pc_x, pc_z, psg0_44, psg1_44, dsd0_18, \
                         dsd0_19, dsd1_18, dsd1_19, dsf_30, dsf_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = pa_x[k] * psg0_44[k]
                  - f_6 * pc_x[k] * psg1_44[k];

        t_45[k] = f_1 * dsd0_18[k]
                  - f_2 * dsd1_18[k]
                  + f_3 * pc_x[k] * dsf_30[k];

        t_46[k] = f_8 * dsd0_19[k]
                  - f_9 * dsd1_19[k]
                  + f_3 * pc_x[k] * dsf_31[k];

        t_47[k] = f_3 * pc_z[k] * dsf_30[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, pc_x, pc_z, dsd0_21, dsd0_23, dsd1_21, \
                         dsd1_23, dsf_31, dsf_33, dsf_35, dsf_36, \
                         dsf_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_4 * dsd0_21[k]
                  - f_5 * dsd1_21[k]
                  + f_3 * pc_x[k] * dsf_33[k];

        t_49[k] = f_3 * pc_z[k] * dsf_31[k];

        t_50[k] = f_4 * dsd0_23[k]
                  - f_5 * dsd1_23[k]
                  + f_3 * pc_x[k] * dsf_35[k];

        t_51[k] = f_3 * pc_x[k] * dsf_36[k];

        t_52[k] = f_3 * pc_x[k] * dsf_37[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, pc_x, pc_y, pc_z, psf_16, dsd0_21, \
                         dsd1_21, dsf_36, dsf_37, dsf_38, dsf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_3 * pc_x[k] * dsf_38[k];

        t_54[k] = f_3 * pc_x[k] * dsf_39[k];

        t_55[k] = f_0 * psf_16[k]
                  + f_1 * dsd0_21[k]
                  - f_2 * dsd1_21[k]
                  + f_3 * pc_y[k] * dsf_36[k];

        t_56[k] = f_3 * pc_z[k] * dsf_36[k];

        t_57[k] = f_4 * dsd0_21[k]
                  - f_5 * dsd1_21[k]
                  + f_3 * pc_z[k] * dsf_37[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, pa_y, pa_z, pc_y, pc_z, psg0_16, psg0_30, \
                         psf_19, psg1_16, psg1_30, dsd0_23, dsd1_23, \
                         dsf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_0 * psf_19[k]
                  + f_3 * pc_y[k] * dsf_39[k];

        t_59[k] = f_1 * dsd0_23[k]
                  - f_2 * dsd1_23[k]
                  + f_3 * pc_z[k] * dsf_39[k];

        t_60[k] = pa_y[k] * psg0_30[k]
                  - f_6 * pc_y[k] * psg1_30[k];

        t_61[k] = pa_z[k] * psg0_16[k]
                  - f_6 * pc_z[k] * psg1_16[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, pa_y, pa_z, pc_x, pc_y, pc_z, psg0_18, psg0_32, \
                         psg1_18, psg1_32, dsd0_28, dsd1_28, dsf_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = pa_y[k] * psg0_32[k]
                  - f_6 * pc_y[k] * psg1_32[k];

        t_63[k] = pa_z[k] * psg0_18[k]
                  - f_6 * pc_z[k] * psg1_18[k];

        t_64[k] = f_4 * dsd0_28[k]
                  - f_5 * dsd1_28[k]
                  + f_3 * pc_x[k] * dsf_44[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, pa_y, pc_x, pc_y, psg0_35, psg1_35, \
                         dsf_46, dsf_47, dsf_48, dsf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = pa_y[k] * psg0_35[k]
                  - f_6 * pc_y[k] * psg1_35[k];

        t_66[k] = f_3 * pc_x[k] * dsf_46[k];

        t_67[k] = f_3 * pc_x[k] * dsf_47[k];

        t_68[k] = f_3 * pc_x[k] * dsf_48[k];

        t_69[k] = f_3 * pc_x[k] * dsf_49[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, pa_y, pa_z, pc_y, pc_z, psg0_25, psg0_42, psf_16, \
                         psf_28, psg1_25, psg1_42, dsf_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = pa_z[k] * psg0_25[k]
                  - f_6 * pc_z[k] * psg1_25[k];

        t_71[k] = f_7 * psf_16[k]
                  + f_3 * pc_z[k] * dsf_46[k];

        t_72[k] = pa_y[k] * psg0_42[k]
                  + f_0 * psf_28[k]
                  - f_6 * pc_y[k] * psg1_42[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pa_y, pc_x, pc_y, psg0_44, psf_29, psg1_44, \
                         dsd0_30, dsd1_30, dsf_49, dsf_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_7 * psf_29[k]
                  + f_3 * pc_y[k] * dsf_49[k];

        t_74[k] = pa_y[k] * psg0_44[k]
                  - f_6 * pc_y[k] * psg1_44[k];

        t_75[k] = f_1 * dsd0_30[k]
                  - f_2 * dsd1_30[k]
                  + f_3 * pc_x[k] * dsf_50[k];

        t_76[k] = f_3 * pc_y[k] * dsf_50[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, pc_x, pc_y, dsd0_32, dsd0_33, dsd0_35, \
                         dsd1_32, dsd1_33, dsd1_35, dsf_52, dsf_53, \
                         dsf_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_8 * dsd0_32[k]
                  - f_9 * dsd1_32[k]
                  + f_3 * pc_x[k] * dsf_52[k];

        t_78[k] = f_4 * dsd0_33[k]
                  - f_5 * dsd1_33[k]
                  + f_3 * pc_x[k] * dsf_53[k];

        t_79[k] = f_3 * pc_y[k] * dsf_52[k];

        t_80[k] = f_4 * dsd0_35[k]
                  - f_5 * dsd1_35[k]
                  + f_3 * pc_x[k] * dsf_55[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, t_85, t_86, pc_x, pc_y, dsd0_33, dsd0_34, \
                         dsd1_33, dsd1_34, dsf_56, dsf_57, dsf_58, \
                         dsf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_3 * pc_x[k] * dsf_56[k];

        t_82[k] = f_3 * pc_x[k] * dsf_57[k];

        t_83[k] = f_3 * pc_x[k] * dsf_58[k];

        t_84[k] = f_3 * pc_x[k] * dsf_59[k];

        t_85[k] = f_1 * dsd0_33[k]
                  - f_2 * dsd1_33[k]
                  + f_3 * pc_y[k] * dsf_56[k];

        t_86[k] = f_8 * dsd0_34[k]
                  - f_9 * dsd1_34[k]
                  + f_3 * pc_y[k] * dsf_57[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, pc_y, pc_z, psf_29, dsd0_35, dsd1_35, dsf_58, \
                         dsf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_4 * dsd0_35[k]
                  - f_5 * dsd1_35[k]
                  + f_3 * pc_y[k] * dsf_58[k];

        t_88[k] = f_3 * pc_y[k] * dsf_59[k];

        t_89[k] = f_0 * psf_29[k]
                  + f_1 * dsd0_35[k]
                  - f_2 * dsd1_35[k]
                  + f_3 * pc_z[k] * dsf_59[k];
    }
}

}  // namespace simdt3ceri
