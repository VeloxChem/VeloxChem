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


#include "SimdThreeCenterElectronRepulsionVrrRecDSF.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_dsf_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t psf0, const size_t psd,
                                                   const size_t psf1, const size_t dsp0,
                                                   const size_t dsp1, const size_t dsd,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / q;
    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 0.5 / gamma;
    const auto f_7 = 0.5 * p / (gamma * q);

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

    const auto *psf0_0 = buffer.data(psf0 + 0);
    const auto *psf0_11 = buffer.data(psf0 + 11);
    const auto *psf0_16 = buffer.data(psf0 + 16);
    const auto *psf0_19 = buffer.data(psf0 + 19);
    const auto *psf0_20 = buffer.data(psf0 + 20);
    const auto *psf0_22 = buffer.data(psf0 + 22);
    const auto *psf0_26 = buffer.data(psf0 + 26);
    const auto *psf0_27 = buffer.data(psf0 + 27);
    const auto *psf0_29 = buffer.data(psf0 + 29);

    const auto *psd_0 = buffer.data(psd + 0);
    const auto *psd_3 = buffer.data(psd + 3);
    const auto *psd_5 = buffer.data(psd + 5);
    const auto *psd_9 = buffer.data(psd + 9);
    const auto *psd_11 = buffer.data(psd + 11);
    const auto *psd_15 = buffer.data(psd + 15);
    const auto *psd_17 = buffer.data(psd + 17);

    const auto *psf1_0 = buffer.data(psf1 + 0);
    const auto *psf1_11 = buffer.data(psf1 + 11);
    const auto *psf1_16 = buffer.data(psf1 + 16);
    const auto *psf1_19 = buffer.data(psf1 + 19);
    const auto *psf1_20 = buffer.data(psf1 + 20);
    const auto *psf1_22 = buffer.data(psf1 + 22);
    const auto *psf1_26 = buffer.data(psf1 + 26);
    const auto *psf1_27 = buffer.data(psf1 + 27);
    const auto *psf1_29 = buffer.data(psf1 + 29);

    const auto *dsp0_0 = buffer.data(dsp0 + 0);
    const auto *dsp0_1 = buffer.data(dsp0 + 1);
    const auto *dsp0_2 = buffer.data(dsp0 + 2);
    const auto *dsp0_9 = buffer.data(dsp0 + 9);
    const auto *dsp0_10 = buffer.data(dsp0 + 10);
    const auto *dsp0_11 = buffer.data(dsp0 + 11);
    const auto *dsp0_15 = buffer.data(dsp0 + 15);
    const auto *dsp0_16 = buffer.data(dsp0 + 16);
    const auto *dsp0_17 = buffer.data(dsp0 + 17);

    const auto *dsp1_0 = buffer.data(dsp1 + 0);
    const auto *dsp1_1 = buffer.data(dsp1 + 1);
    const auto *dsp1_2 = buffer.data(dsp1 + 2);
    const auto *dsp1_9 = buffer.data(dsp1 + 9);
    const auto *dsp1_10 = buffer.data(dsp1 + 10);
    const auto *dsp1_11 = buffer.data(dsp1 + 11);
    const auto *dsp1_15 = buffer.data(dsp1 + 15);
    const auto *dsp1_16 = buffer.data(dsp1 + 16);
    const auto *dsp1_17 = buffer.data(dsp1 + 17);

    const auto *dsd_0 = buffer.data(dsd + 0);
    const auto *dsd_2 = buffer.data(dsd + 2);
    const auto *dsd_3 = buffer.data(dsd + 3);
    const auto *dsd_5 = buffer.data(dsd + 5);
    const auto *dsd_6 = buffer.data(dsd + 6);
    const auto *dsd_7 = buffer.data(dsd + 7);
    const auto *dsd_9 = buffer.data(dsd + 9);
    const auto *dsd_11 = buffer.data(dsd + 11);
    const auto *dsd_12 = buffer.data(dsd + 12);
    const auto *dsd_14 = buffer.data(dsd + 14);
    const auto *dsd_15 = buffer.data(dsd + 15);
    const auto *dsd_17 = buffer.data(dsd + 17);
    const auto *dsd_18 = buffer.data(dsd + 18);
    const auto *dsd_19 = buffer.data(dsd + 19);
    const auto *dsd_21 = buffer.data(dsd + 21);
    const auto *dsd_22 = buffer.data(dsd + 22);
    const auto *dsd_23 = buffer.data(dsd + 23);
    const auto *dsd_27 = buffer.data(dsd + 27);
    const auto *dsd_28 = buffer.data(dsd + 28);
    const auto *dsd_29 = buffer.data(dsd + 29);
    const auto *dsd_30 = buffer.data(dsd + 30);
    const auto *dsd_32 = buffer.data(dsd + 32);
    const auto *dsd_33 = buffer.data(dsd + 33);
    const auto *dsd_34 = buffer.data(dsd + 34);
    const auto *dsd_35 = buffer.data(dsd + 35);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pc_x, pc_y, pc_z, psd_0, psd_3, dsp0_0, \
                         dsp1_0, dsd_0, dsd_2, dsd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * psd_0[k]
                 + f_1 * dsp0_0[k]
                 - f_2 * dsp1_0[k]
                 + f_3 * pc_x[k] * dsd_0[k];

        t_1[k] = f_3 * pc_y[k] * dsd_0[k];

        t_2[k] = f_3 * pc_z[k] * dsd_0[k];

        t_3[k] = f_0 * psd_3[k]
                 + f_3 * pc_x[k] * dsd_3[k];

        t_4[k] = f_3 * pc_y[k] * dsd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pc_x, pc_y, pc_z, psd_5, dsp0_1, dsp0_2, \
                         dsp1_1, dsp1_2, dsd_3, dsd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * psd_5[k]
                 + f_3 * pc_x[k] * dsd_5[k];

        t_6[k] = f_1 * dsp0_1[k]
                 - f_2 * dsp1_1[k]
                 + f_3 * pc_y[k] * dsd_3[k];

        t_7[k] = f_3 * pc_z[k] * dsd_3[k];

        t_8[k] = f_3 * pc_y[k] * dsd_5[k];

        t_9[k] = f_1 * dsp0_2[k]
                 - f_2 * dsp1_2[k]
                 + f_3 * pc_z[k] * dsd_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_y, pc_x, pc_y, pc_z, psf0_0, psd_0, \
                         psd_9, psf1_0, dsd_6, dsd_7, dsd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_y[k] * psf0_0[k]
                  - f_4 * pc_y[k] * psf1_0[k];

        t_11[k] = f_5 * psd_0[k]
                  + f_3 * pc_y[k] * dsd_6[k];

        t_12[k] = f_3 * pc_z[k] * dsd_6[k];

        t_13[k] = f_5 * psd_9[k]
                  + f_3 * pc_x[k] * dsd_9[k];

        t_14[k] = f_3 * pc_z[k] * dsd_7[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pa_x, pc_x, pc_y, pc_z, psf0_16, psd_5, \
                         psd_11, psf1_16, dsd_9, dsd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_5 * psd_11[k]
                  + f_3 * pc_x[k] * dsd_11[k];

        t_16[k] = pa_x[k] * psf0_16[k]
                  - f_4 * pc_x[k] * psf1_16[k];

        t_17[k] = f_3 * pc_z[k] * dsd_9[k];

        t_18[k] = f_5 * psd_5[k]
                  + f_3 * pc_y[k] * dsd_11[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_x, pa_z, pc_x, pc_y, pc_z, psf0_0, \
                         psf0_19, psd_0, psf1_0, psf1_19, dsd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pa_x[k] * psf0_19[k]
                  - f_4 * pc_x[k] * psf1_19[k];

        t_20[k] = pa_z[k] * psf0_0[k]
                  - f_4 * pc_z[k] * psf1_0[k];

        t_21[k] = f_3 * pc_y[k] * dsd_12[k];

        t_22[k] = f_5 * psd_0[k]
                  + f_3 * pc_z[k] * dsd_12[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_x, pc_x, pc_y, psf0_26, psd_15, psd_17, \
                         psf1_26, dsd_14, dsd_15, dsd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_5 * psd_15[k]
                  + f_3 * pc_x[k] * dsd_15[k];

        t_24[k] = f_3 * pc_y[k] * dsd_14[k];

        t_25[k] = f_5 * psd_17[k]
                  + f_3 * pc_x[k] * dsd_17[k];

        t_26[k] = pa_x[k] * psf0_26[k]
                  - f_4 * pc_x[k] * psf1_26[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pa_x, pc_x, pc_y, psf0_27, psf0_29, psf1_27, \
                         psf1_29, dsp0_9, dsp1_9, dsd_17, dsd_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = pa_x[k] * psf0_27[k]
                  - f_4 * pc_x[k] * psf1_27[k];

        t_28[k] = f_3 * pc_y[k] * dsd_17[k];

        t_29[k] = pa_x[k] * psf0_29[k]
                  - f_4 * pc_x[k] * psf1_29[k];

        t_30[k] = f_1 * dsp0_9[k]
                  - f_2 * dsp1_9[k]
                  + f_3 * pc_x[k] * dsd_18[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, pc_x, pc_z, dsp0_10, dsp1_10, dsd_18, \
                         dsd_19, dsd_21, dsd_22, dsd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_6 * dsp0_10[k]
                  - f_7 * dsp1_10[k]
                  + f_3 * pc_x[k] * dsd_19[k];

        t_32[k] = f_3 * pc_z[k] * dsd_18[k];

        t_33[k] = f_3 * pc_x[k] * dsd_21[k];

        t_34[k] = f_3 * pc_x[k] * dsd_22[k];

        t_35[k] = f_3 * pc_x[k] * dsd_23[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pc_y, pc_z, psd_9, psd_11, dsp0_10, dsp0_11, \
                         dsp1_10, dsp1_11, dsd_21, dsd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_0 * psd_9[k]
                  + f_1 * dsp0_10[k]
                  - f_2 * dsp1_10[k]
                  + f_3 * pc_y[k] * dsd_21[k];

        t_37[k] = f_3 * pc_z[k] * dsd_21[k];

        t_38[k] = f_0 * psd_11[k]
                  + f_3 * pc_y[k] * dsd_23[k];

        t_39[k] = f_1 * dsp0_11[k]
                  - f_2 * dsp1_11[k]
                  + f_3 * pc_z[k] * dsd_23[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_y, pa_z, pc_x, pc_y, pc_z, psf0_11, \
                         psf0_20, psf0_22, psf1_11, psf1_20, psf1_22, \
                         dsd_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pa_y[k] * psf0_20[k]
                  - f_4 * pc_y[k] * psf1_20[k];

        t_41[k] = pa_z[k] * psf0_11[k]
                  - f_4 * pc_z[k] * psf1_11[k];

        t_42[k] = pa_y[k] * psf0_22[k]
                  - f_4 * pc_y[k] * psf1_22[k];

        t_43[k] = f_3 * pc_x[k] * dsd_27[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, pa_z, pc_x, pc_y, pc_z, psf0_16, psd_9, \
                         psd_17, psf1_16, dsd_27, dsd_28, dsd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_3 * pc_x[k] * dsd_28[k];

        t_45[k] = f_3 * pc_x[k] * dsd_29[k];

        t_46[k] = pa_z[k] * psf0_16[k]
                  - f_4 * pc_z[k] * psf1_16[k];

        t_47[k] = f_5 * psd_9[k]
                  + f_3 * pc_z[k] * dsd_27[k];

        t_48[k] = f_5 * psd_17[k]
                  + f_3 * pc_y[k] * dsd_29[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pa_y, pc_x, pc_y, psf0_29, psf1_29, dsp0_15, \
                         dsp0_17, dsp1_15, dsp1_17, dsd_30, dsd_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = pa_y[k] * psf0_29[k]
                  - f_4 * pc_y[k] * psf1_29[k];

        t_50[k] = f_1 * dsp0_15[k]
                  - f_2 * dsp1_15[k]
                  + f_3 * pc_x[k] * dsd_30[k];

        t_51[k] = f_3 * pc_y[k] * dsd_30[k];

        t_52[k] = f_6 * dsp0_17[k]
                  - f_7 * dsp1_17[k]
                  + f_3 * pc_x[k] * dsd_32[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, t_58, pc_x, pc_y, dsp0_16, dsp0_17, \
                         dsp1_16, dsp1_17, dsd_33, dsd_34, dsd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_3 * pc_x[k] * dsd_33[k];

        t_54[k] = f_3 * pc_x[k] * dsd_34[k];

        t_55[k] = f_3 * pc_x[k] * dsd_35[k];

        t_56[k] = f_1 * dsp0_16[k]
                  - f_2 * dsp1_16[k]
                  + f_3 * pc_y[k] * dsd_33[k];

        t_57[k] = f_6 * dsp0_17[k]
                  - f_7 * dsp1_17[k]
                  + f_3 * pc_y[k] * dsd_34[k];

        t_58[k] = f_3 * pc_y[k] * dsd_35[k];
    }

#pragma omp simd aligned(t_59, pc_z, psd_17, dsp0_17, dsp1_17, dsd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_0 * psd_17[k]
                  + f_1 * dsp0_17[k]
                  - f_2 * dsp1_17[k]
                  + f_3 * pc_z[k] * dsd_35[k];
    }
}

}  // namespace simdt3ceri
