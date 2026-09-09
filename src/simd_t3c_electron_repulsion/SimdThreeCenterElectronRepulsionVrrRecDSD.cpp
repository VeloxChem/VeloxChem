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


#include "SimdThreeCenterElectronRepulsionVrrRecDSD.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_dsd_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t psd0, const size_t psp,
                                                   const size_t psd1, const size_t dss0,
                                                   const size_t dss1, const size_t dsp,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / q;
    const auto f_1 = 0.5 / gamma;
    const auto f_2 = 0.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *psd0_0 = buffer.data(psd0 + 0);
    const auto *psd0_9 = buffer.data(psd0 + 9);
    const auto *psd0_11 = buffer.data(psd0 + 11);
    const auto *psd0_12 = buffer.data(psd0 + 12);
    const auto *psd0_15 = buffer.data(psd0 + 15);
    const auto *psd0_17 = buffer.data(psd0 + 17);

    const auto *psp_0 = buffer.data(psp + 0);
    const auto *psp_4 = buffer.data(psp + 4);
    const auto *psp_8 = buffer.data(psp + 8);

    const auto *psd1_0 = buffer.data(psd1 + 0);
    const auto *psd1_9 = buffer.data(psd1 + 9);
    const auto *psd1_11 = buffer.data(psd1 + 11);
    const auto *psd1_12 = buffer.data(psd1 + 12);
    const auto *psd1_15 = buffer.data(psd1 + 15);
    const auto *psd1_17 = buffer.data(psd1 + 17);

    const auto *dss0_0 = buffer.data(dss0 + 0);
    const auto *dss0_3 = buffer.data(dss0 + 3);
    const auto *dss0_5 = buffer.data(dss0 + 5);

    const auto *dss1_0 = buffer.data(dss1 + 0);
    const auto *dss1_3 = buffer.data(dss1 + 3);
    const auto *dss1_5 = buffer.data(dss1 + 5);

    const auto *dsp_0 = buffer.data(dsp + 0);
    const auto *dsp_1 = buffer.data(dsp + 1);
    const auto *dsp_2 = buffer.data(dsp + 2);
    const auto *dsp_3 = buffer.data(dsp + 3);
    const auto *dsp_4 = buffer.data(dsp + 4);
    const auto *dsp_6 = buffer.data(dsp + 6);
    const auto *dsp_8 = buffer.data(dsp + 8);
    const auto *dsp_9 = buffer.data(dsp + 9);
    const auto *dsp_10 = buffer.data(dsp + 10);
    const auto *dsp_11 = buffer.data(dsp + 11);
    const auto *dsp_13 = buffer.data(dsp + 13);
    const auto *dsp_14 = buffer.data(dsp + 14);
    const auto *dsp_15 = buffer.data(dsp + 15);
    const auto *dsp_16 = buffer.data(dsp + 16);
    const auto *dsp_17 = buffer.data(dsp + 17);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, psp_0, dss0_0, \
                         dss1_0, dsp_0, dsp_1, dsp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * psp_0[k]
                 + f_1 * dss0_0[k]
                 - f_2 * dss1_0[k]
                 + f_3 * pc_x[k] * dsp_0[k];

        t_1[k] = f_3 * pc_y[k] * dsp_0[k];

        t_2[k] = f_3 * pc_z[k] * dsp_0[k];

        t_3[k] = f_1 * dss0_0[k]
                 - f_2 * dss1_0[k]
                 + f_3 * pc_y[k] * dsp_1[k];

        t_4[k] = f_3 * pc_y[k] * dsp_2[k];

        t_5[k] = f_1 * dss0_0[k]
                 - f_2 * dss1_0[k]
                 + f_3 * pc_z[k] * dsp_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pa_x, pa_y, pc_x, pc_y, pc_z, psd0_0, psd0_9, \
                         psp_4, psd1_0, psd1_9, dsp_3, dsp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pa_y[k] * psd0_0[k]
                 - f_4 * pc_y[k] * psd1_0[k];

        t_7[k] = f_5 * psp_4[k]
                 + f_3 * pc_x[k] * dsp_4[k];

        t_8[k] = f_3 * pc_z[k] * dsp_3[k];

        t_9[k] = pa_x[k] * psd0_9[k]
                 - f_4 * pc_x[k] * psd1_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pa_z, pc_x, pc_y, pc_z, psd0_0, \
                         psd0_11, psd1_0, psd1_11, dsp_4, dsp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * pc_z[k] * dsp_4[k];

        t_11[k] = pa_x[k] * psd0_11[k]
                  - f_4 * pc_x[k] * psd1_11[k];

        t_12[k] = pa_z[k] * psd0_0[k]
                  - f_4 * pc_z[k] * psd1_0[k];

        t_13[k] = f_3 * pc_y[k] * dsp_6[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_x, pc_x, pc_y, psd0_15, psd0_17, psp_8, \
                         psd1_15, psd1_17, dsp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_5 * psp_8[k]
                  + f_3 * pc_x[k] * dsp_8[k];

        t_15[k] = pa_x[k] * psd0_15[k]
                  - f_4 * pc_x[k] * psd1_15[k];

        t_16[k] = f_3 * pc_y[k] * dsp_8[k];

        t_17[k] = pa_x[k] * psd0_17[k]
                  - f_4 * pc_x[k] * psd1_17[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, t_23, pc_x, pc_y, pc_z, psp_4, dss0_3, \
                         dss1_3, dsp_9, dsp_10, dsp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_1 * dss0_3[k]
                  - f_2 * dss1_3[k]
                  + f_3 * pc_x[k] * dsp_9[k];

        t_19[k] = f_3 * pc_x[k] * dsp_10[k];

        t_20[k] = f_3 * pc_x[k] * dsp_11[k];

        t_21[k] = f_0 * psp_4[k]
                  + f_1 * dss0_3[k]
                  - f_2 * dss1_3[k]
                  + f_3 * pc_y[k] * dsp_10[k];

        t_22[k] = f_3 * pc_z[k] * dsp_10[k];

        t_23[k] = f_1 * dss0_3[k]
                  - f_2 * dss1_3[k]
                  + f_3 * pc_z[k] * dsp_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_y, pa_z, pc_x, pc_y, pc_z, psd0_9, \
                         psd0_12, psd1_9, psd1_12, dsp_13, dsp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pa_y[k] * psd0_12[k]
                  - f_4 * pc_y[k] * psd1_12[k];

        t_25[k] = f_3 * pc_x[k] * dsp_13[k];

        t_26[k] = f_3 * pc_x[k] * dsp_14[k];

        t_27[k] = pa_z[k] * psd0_9[k]
                  - f_4 * pc_z[k] * psd1_9[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_y, pc_x, pc_y, psd0_17, psp_8, psd1_17, \
                         dss0_5, dss1_5, dsp_14, dsp_15, dsp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_5 * psp_8[k]
                  + f_3 * pc_y[k] * dsp_14[k];

        t_29[k] = pa_y[k] * psd0_17[k]
                  - f_4 * pc_y[k] * psd1_17[k];

        t_30[k] = f_1 * dss0_5[k]
                  - f_2 * dss1_5[k]
                  + f_3 * pc_x[k] * dsp_15[k];

        t_31[k] = f_3 * pc_x[k] * dsp_16[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pc_x, pc_y, pc_z, psp_8, dss0_5, dss1_5, \
                         dsp_16, dsp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_3 * pc_x[k] * dsp_17[k];

        t_33[k] = f_1 * dss0_5[k]
                  - f_2 * dss1_5[k]
                  + f_3 * pc_y[k] * dsp_16[k];

        t_34[k] = f_3 * pc_y[k] * dsp_17[k];

        t_35[k] = f_0 * psp_8[k]
                  + f_1 * dss0_5[k]
                  - f_2 * dss1_5[k]
                  + f_3 * pc_z[k] * dsp_17[k];
    }
}

}  // namespace simdt3ceri
