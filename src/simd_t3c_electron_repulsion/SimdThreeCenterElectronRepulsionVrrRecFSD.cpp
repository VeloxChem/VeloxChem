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


#include "SimdThreeCenterElectronRepulsionVrrRecFSD.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_fsd_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t dsd0, const size_t dsp,
                                                   const size_t dsd1, const size_t fss0,
                                                   const size_t fss1, const size_t fsp,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
    const auto f_1 = 0.5 / gamma;
    const auto f_2 = 0.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 1.0 / q;
    const auto f_6 = 0.5 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *dsd0_0 = buffer.data(dsd0 + 0);
    const auto *dsd0_3 = buffer.data(dsd0 + 3);
    const auto *dsd0_5 = buffer.data(dsd0 + 5);
    const auto *dsd0_12 = buffer.data(dsd0 + 12);
    const auto *dsd0_18 = buffer.data(dsd0 + 18);
    const auto *dsd0_21 = buffer.data(dsd0 + 21);
    const auto *dsd0_23 = buffer.data(dsd0 + 23);
    const auto *dsd0_27 = buffer.data(dsd0 + 27);
    const auto *dsd0_29 = buffer.data(dsd0 + 29);
    const auto *dsd0_30 = buffer.data(dsd0 + 30);
    const auto *dsd0_33 = buffer.data(dsd0 + 33);
    const auto *dsd0_35 = buffer.data(dsd0 + 35);

    const auto *dsp_0 = buffer.data(dsp + 0);
    const auto *dsp_1 = buffer.data(dsp + 1);
    const auto *dsp_2 = buffer.data(dsp + 2);
    const auto *dsp_4 = buffer.data(dsp + 4);
    const auto *dsp_8 = buffer.data(dsp + 8);
    const auto *dsp_9 = buffer.data(dsp + 9);
    const auto *dsp_10 = buffer.data(dsp + 10);
    const auto *dsp_11 = buffer.data(dsp + 11);
    const auto *dsp_13 = buffer.data(dsp + 13);
    const auto *dsp_14 = buffer.data(dsp + 14);
    const auto *dsp_15 = buffer.data(dsp + 15);
    const auto *dsp_16 = buffer.data(dsp + 16);
    const auto *dsp_17 = buffer.data(dsp + 17);

    const auto *dsd1_0 = buffer.data(dsd1 + 0);
    const auto *dsd1_3 = buffer.data(dsd1 + 3);
    const auto *dsd1_5 = buffer.data(dsd1 + 5);
    const auto *dsd1_12 = buffer.data(dsd1 + 12);
    const auto *dsd1_18 = buffer.data(dsd1 + 18);
    const auto *dsd1_21 = buffer.data(dsd1 + 21);
    const auto *dsd1_23 = buffer.data(dsd1 + 23);
    const auto *dsd1_27 = buffer.data(dsd1 + 27);
    const auto *dsd1_29 = buffer.data(dsd1 + 29);
    const auto *dsd1_30 = buffer.data(dsd1 + 30);
    const auto *dsd1_33 = buffer.data(dsd1 + 33);
    const auto *dsd1_35 = buffer.data(dsd1 + 35);

    const auto *fss0_0 = buffer.data(fss0 + 0);
    const auto *fss0_1 = buffer.data(fss0 + 1);
    const auto *fss0_2 = buffer.data(fss0 + 2);
    const auto *fss0_6 = buffer.data(fss0 + 6);
    const auto *fss0_7 = buffer.data(fss0 + 7);
    const auto *fss0_9 = buffer.data(fss0 + 9);

    const auto *fss1_0 = buffer.data(fss1 + 0);
    const auto *fss1_1 = buffer.data(fss1 + 1);
    const auto *fss1_2 = buffer.data(fss1 + 2);
    const auto *fss1_6 = buffer.data(fss1 + 6);
    const auto *fss1_7 = buffer.data(fss1 + 7);
    const auto *fss1_9 = buffer.data(fss1 + 9);

    const auto *fsp_0 = buffer.data(fsp + 0);
    const auto *fsp_1 = buffer.data(fsp + 1);
    const auto *fsp_2 = buffer.data(fsp + 2);
    const auto *fsp_3 = buffer.data(fsp + 3);
    const auto *fsp_4 = buffer.data(fsp + 4);
    const auto *fsp_6 = buffer.data(fsp + 6);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, dsp_0, fss0_0, \
                         fss1_0, fsp_0, fsp_1, fsp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dsp_0[k]
                 + f_1 * fss0_0[k]
                 - f_2 * fss1_0[k]
                 + f_3 * pc_x[k] * fsp_0[k];

        t_1[k] = f_3 * pc_y[k] * fsp_0[k];

        t_2[k] = f_3 * pc_z[k] * fsp_0[k];

        t_3[k] = f_1 * fss0_0[k]
                 - f_2 * fss1_0[k]
                 + f_3 * pc_y[k] * fsp_1[k];

        t_4[k] = f_3 * pc_y[k] * fsp_2[k];

        t_5[k] = f_1 * fss0_0[k]
                 - f_2 * fss1_0[k]
                 + f_3 * pc_z[k] * fsp_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pa_y, pc_x, pc_y, pc_z, dsd0_0, dsp_1, dsp_4, \
                         dsd1_0, fss0_1, fss1_1, fsp_3, fsp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pa_y[k] * dsd0_0[k]
                 - f_4 * pc_y[k] * dsd1_0[k];

        t_7[k] = f_5 * dsp_4[k]
                 + f_3 * pc_x[k] * fsp_4[k];

        t_8[k] = f_3 * pc_z[k] * fsp_3[k];

        t_9[k] = f_6 * dsp_1[k]
                 + f_1 * fss0_1[k]
                 - f_2 * fss1_1[k]
                 + f_3 * pc_y[k] * fsp_4[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_y, pa_z, pc_y, pc_z, dsd0_0, dsd0_5, \
                         dsd1_0, dsd1_5, fsp_4, fsp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * pc_z[k] * fsp_4[k];

        t_11[k] = pa_y[k] * dsd0_5[k]
                  - f_4 * pc_y[k] * dsd1_5[k];

        t_12[k] = pa_z[k] * dsd0_0[k]
                  - f_4 * pc_z[k] * dsd1_0[k];

        t_13[k] = f_3 * pc_y[k] * fsp_6[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_z, pc_x, pc_y, pc_z, dsd0_3, dsp_2, dsp_8, \
                         dsd1_3, fss0_2, fss1_2, fsp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_5 * dsp_8[k]
                  + f_3 * pc_x[k] * fsp_8[k];

        t_15[k] = pa_z[k] * dsd0_3[k]
                  - f_4 * pc_z[k] * dsd1_3[k];

        t_16[k] = f_3 * pc_y[k] * fsp_8[k];

        t_17[k] = f_6 * dsp_2[k]
                  + f_1 * fss0_2[k]
                  - f_2 * fss1_2[k]
                  + f_3 * pc_z[k] * fsp_8[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, pa_x, pc_x, pc_z, dsd0_18, dsd0_21, \
                         dsp_9, dsp_10, dsd1_18, dsd1_21, fsp_9, \
                         fsp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = pa_x[k] * dsd0_18[k]
                  + f_5 * dsp_9[k]
                  - f_4 * pc_x[k] * dsd1_18[k];

        t_19[k] = f_6 * dsp_10[k]
                  + f_3 * pc_x[k] * fsp_10[k];

        t_20[k] = f_3 * pc_z[k] * fsp_9[k];

        t_21[k] = pa_x[k] * dsd0_21[k]
                  - f_4 * pc_x[k] * dsd1_21[k];

        t_22[k] = f_3 * pc_z[k] * fsp_10[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_x, pa_y, pc_x, pc_y, dsd0_12, dsd0_23, \
                         dsp_13, dsp_14, dsd1_12, dsd1_23, fsp_13, \
                         fsp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = pa_x[k] * dsd0_23[k]
                  - f_4 * pc_x[k] * dsd1_23[k];

        t_24[k] = pa_y[k] * dsd0_12[k]
                  - f_4 * pc_y[k] * dsd1_12[k];

        t_25[k] = f_6 * dsp_13[k]
                  + f_3 * pc_x[k] * fsp_13[k];

        t_26[k] = f_6 * dsp_14[k]
                  + f_3 * pc_x[k] * fsp_14[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pa_x, pc_x, pc_y, dsd0_27, dsd0_29, dsd0_30, \
                         dsp_8, dsp_15, dsd1_27, dsd1_29, dsd1_30, \
                         fsp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = pa_x[k] * dsd0_27[k]
                  - f_4 * pc_x[k] * dsd1_27[k];

        t_28[k] = f_6 * dsp_8[k]
                  + f_3 * pc_y[k] * fsp_14[k];

        t_29[k] = pa_x[k] * dsd0_29[k]
                  - f_4 * pc_x[k] * dsd1_29[k];

        t_30[k] = pa_x[k] * dsd0_30[k]
                  + f_5 * dsp_15[k]
                  - f_4 * pc_x[k] * dsd1_30[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, pa_x, pc_x, pc_y, dsd0_33, dsd0_35, \
                         dsp_17, dsd1_33, dsd1_35, fsp_15, fsp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_3 * pc_y[k] * fsp_15[k];

        t_32[k] = f_6 * dsp_17[k]
                  + f_3 * pc_x[k] * fsp_17[k];

        t_33[k] = pa_x[k] * dsd0_33[k]
                  - f_4 * pc_x[k] * dsd1_33[k];

        t_34[k] = f_3 * pc_y[k] * fsp_17[k];

        t_35[k] = pa_x[k] * dsd0_35[k]
                  - f_4 * pc_x[k] * dsd1_35[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, t_41, pc_x, pc_y, pc_z, dsp_10, fss0_6, \
                         fss1_6, fsp_18, fsp_19, fsp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_1 * fss0_6[k]
                  - f_2 * fss1_6[k]
                  + f_3 * pc_x[k] * fsp_18[k];

        t_37[k] = f_3 * pc_x[k] * fsp_19[k];

        t_38[k] = f_3 * pc_x[k] * fsp_20[k];

        t_39[k] = f_0 * dsp_10[k]
                  + f_1 * fss0_6[k]
                  - f_2 * fss1_6[k]
                  + f_3 * pc_y[k] * fsp_19[k];

        t_40[k] = f_3 * pc_z[k] * fsp_19[k];

        t_41[k] = f_1 * fss0_6[k]
                  - f_2 * fss1_6[k]
                  + f_3 * pc_z[k] * fsp_20[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, pa_z, pc_x, pc_y, pc_z, dsd0_18, \
                         dsd0_21, dsp_14, dsd1_18, dsd1_21, fsp_22, \
                         fsp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pa_z[k] * dsd0_18[k]
                  - f_4 * pc_z[k] * dsd1_18[k];

        t_43[k] = f_3 * pc_x[k] * fsp_22[k];

        t_44[k] = f_3 * pc_x[k] * fsp_23[k];

        t_45[k] = pa_z[k] * dsd0_21[k]
                  - f_4 * pc_z[k] * dsd1_21[k];

        t_46[k] = f_5 * dsp_14[k]
                  + f_3 * pc_y[k] * fsp_23[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pa_y, pc_x, pc_y, pc_z, dsd0_30, dsp_11, \
                         dsd1_30, fss0_7, fss1_7, fsp_23, fsp_25, \
                         fsp_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_6 * dsp_11[k]
                  + f_1 * fss0_7[k]
                  - f_2 * fss1_7[k]
                  + f_3 * pc_z[k] * fsp_23[k];

        t_48[k] = pa_y[k] * dsd0_30[k]
                  - f_4 * pc_y[k] * dsd1_30[k];

        t_49[k] = f_3 * pc_x[k] * fsp_25[k];

        t_50[k] = f_3 * pc_x[k] * fsp_26[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, pa_y, pc_y, dsd0_33, dsd0_35, dsp_16, dsp_17, \
                         dsd1_33, dsd1_35, fsp_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = pa_y[k] * dsd0_33[k]
                  + f_5 * dsp_16[k]
                  - f_4 * pc_y[k] * dsd1_33[k];

        t_52[k] = f_6 * dsp_17[k]
                  + f_3 * pc_y[k] * fsp_26[k];

        t_53[k] = pa_y[k] * dsd0_35[k]
                  - f_4 * pc_y[k] * dsd1_35[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, t_59, pc_x, pc_y, pc_z, dsp_17, fss0_9, \
                         fss1_9, fsp_27, fsp_28, fsp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_1 * fss0_9[k]
                  - f_2 * fss1_9[k]
                  + f_3 * pc_x[k] * fsp_27[k];

        t_55[k] = f_3 * pc_x[k] * fsp_28[k];

        t_56[k] = f_3 * pc_x[k] * fsp_29[k];

        t_57[k] = f_1 * fss0_9[k]
                  - f_2 * fss1_9[k]
                  + f_3 * pc_y[k] * fsp_28[k];

        t_58[k] = f_3 * pc_y[k] * fsp_29[k];

        t_59[k] = f_0 * dsp_17[k]
                  + f_1 * fss0_9[k]
                  - f_2 * fss1_9[k]
                  + f_3 * pc_z[k] * fsp_29[k];
    }
}

}  // namespace simdt3ceri
