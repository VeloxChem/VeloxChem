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


#include "SimdThreeCenterElectronRepulsionVrrRecPSG.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_psg_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t ssg0, const size_t ssf,
                                                   const size_t ssg1, const size_t psd0,
                                                   const size_t psd1, const size_t psf,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = gamma / q;
    const auto f_2 = p / q;
    const auto f_3 = 0.5 / gamma;
    const auto f_4 = 0.5 * p / (gamma * q);
    const auto f_5 = 0.5 / q;
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ssg0_0 = buffer.data(ssg0 + 0);
    const auto *ssg0_3 = buffer.data(ssg0 + 3);
    const auto *ssg0_5 = buffer.data(ssg0 + 5);
    const auto *ssg0_10 = buffer.data(ssg0 + 10);
    const auto *ssg0_12 = buffer.data(ssg0 + 12);
    const auto *ssg0_14 = buffer.data(ssg0 + 14);

    const auto *ssf_0 = buffer.data(ssf + 0);
    const auto *ssf_6 = buffer.data(ssf + 6);
    const auto *ssf_9 = buffer.data(ssf + 9);

    const auto *ssg1_0 = buffer.data(ssg1 + 0);
    const auto *ssg1_3 = buffer.data(ssg1 + 3);
    const auto *ssg1_5 = buffer.data(ssg1 + 5);
    const auto *ssg1_10 = buffer.data(ssg1 + 10);
    const auto *ssg1_12 = buffer.data(ssg1 + 12);
    const auto *ssg1_14 = buffer.data(ssg1 + 14);

    const auto *psd0_0 = buffer.data(psd0 + 0);
    const auto *psd0_7 = buffer.data(psd0 + 7);
    const auto *psd0_9 = buffer.data(psd0 + 9);
    const auto *psd0_14 = buffer.data(psd0 + 14);
    const auto *psd0_16 = buffer.data(psd0 + 16);
    const auto *psd0_17 = buffer.data(psd0 + 17);

    const auto *psd1_0 = buffer.data(psd1 + 0);
    const auto *psd1_7 = buffer.data(psd1 + 7);
    const auto *psd1_9 = buffer.data(psd1 + 9);
    const auto *psd1_14 = buffer.data(psd1 + 14);
    const auto *psd1_16 = buffer.data(psd1 + 16);
    const auto *psd1_17 = buffer.data(psd1 + 17);

    const auto *psf_0 = buffer.data(psf + 0);
    const auto *psf_1 = buffer.data(psf + 1);
    const auto *psf_2 = buffer.data(psf + 2);
    const auto *psf_3 = buffer.data(psf + 3);
    const auto *psf_5 = buffer.data(psf + 5);
    const auto *psf_6 = buffer.data(psf + 6);
    const auto *psf_9 = buffer.data(psf + 9);
    const auto *psf_10 = buffer.data(psf + 10);
    const auto *psf_11 = buffer.data(psf + 11);
    const auto *psf_13 = buffer.data(psf + 13);
    const auto *psf_16 = buffer.data(psf + 16);
    const auto *psf_17 = buffer.data(psf + 17);
    const auto *psf_18 = buffer.data(psf + 18);
    const auto *psf_19 = buffer.data(psf + 19);
    const auto *psf_20 = buffer.data(psf + 20);
    const auto *psf_22 = buffer.data(psf + 22);
    const auto *psf_25 = buffer.data(psf + 25);
    const auto *psf_26 = buffer.data(psf + 26);
    const auto *psf_27 = buffer.data(psf + 27);
    const auto *psf_28 = buffer.data(psf + 28);
    const auto *psf_29 = buffer.data(psf + 29);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pc_x, pc_y, pc_z, ssg0_0, ssf_0, ssg1_0, \
                         psd0_0, psd1_0, psf_0, psf_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = pa_x[k] * ssg0_0[k]
                 + f_0 * ssf_0[k]
                 - f_1 * pc_x[k] * ssg1_0[k];

        t_1[k] = f_2 * pc_y[k] * psf_0[k];

        t_2[k] = f_2 * pc_z[k] * psf_0[k];

        t_3[k] = f_3 * psd0_0[k]
                 - f_4 * psd1_0[k]
                 + f_2 * pc_y[k] * psf_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pc_x, pc_y, pc_z, ssf_6, psd0_0, psd1_0, \
                         psf_2, psf_3, psf_5, psf_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * pc_y[k] * psf_2[k];

        t_5[k] = f_3 * psd0_0[k]
                 - f_4 * psd1_0[k]
                 + f_2 * pc_z[k] * psf_2[k];

        t_6[k] = f_5 * ssf_6[k]
                 + f_2 * pc_x[k] * psf_6[k];

        t_7[k] = f_2 * pc_z[k] * psf_3[k];

        t_8[k] = f_2 * pc_y[k] * psf_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, t_13, pa_x, pc_x, pc_y, pc_z, ssg0_10, \
                         ssg0_12, ssf_9, ssg1_10, ssg1_12, psf_6, \
                         psf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * ssf_9[k]
                 + f_2 * pc_x[k] * psf_9[k];

        t_10[k] = pa_x[k] * ssg0_10[k]
                  - f_1 * pc_x[k] * ssg1_10[k];

        t_11[k] = f_2 * pc_z[k] * psf_6[k];

        t_12[k] = pa_x[k] * ssg0_12[k]
                  - f_1 * pc_x[k] * ssg1_12[k];

        t_13[k] = f_2 * pc_y[k] * psf_9[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_x, pa_y, pc_x, pc_y, ssg0_0, ssg0_14, ssg1_0, \
                         ssg1_14, psd0_7, psd1_7, psf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = pa_x[k] * ssg0_14[k]
                  - f_1 * pc_x[k] * ssg1_14[k];

        t_15[k] = pa_y[k] * ssg0_0[k]
                  - f_1 * pc_y[k] * ssg1_0[k];

        t_16[k] = f_6 * psd0_7[k]
                  - f_7 * psd1_7[k]
                  + f_2 * pc_x[k] * psf_11[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pa_y, pc_x, pc_y, pc_z, ssg0_5, ssg1_5, \
                         psd0_9, psd1_9, psf_10, psf_11, psf_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_2 * pc_z[k] * psf_10[k];

        t_18[k] = f_3 * psd0_9[k]
                  - f_4 * psd1_9[k]
                  + f_2 * pc_x[k] * psf_13[k];

        t_19[k] = f_2 * pc_z[k] * psf_11[k];

        t_20[k] = pa_y[k] * ssg0_5[k]
                  - f_1 * pc_y[k] * ssg1_5[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pa_y, pc_x, pc_y, ssg0_10, ssf_6, \
                         ssg1_10, psf_16, psf_17, psf_18, psf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_2 * pc_x[k] * psf_16[k];

        t_22[k] = f_2 * pc_x[k] * psf_17[k];

        t_23[k] = f_2 * pc_x[k] * psf_18[k];

        t_24[k] = f_2 * pc_x[k] * psf_19[k];

        t_25[k] = pa_y[k] * ssg0_10[k]
                  + f_0 * ssf_6[k]
                  - f_1 * pc_y[k] * ssg1_10[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_y, pc_y, pc_z, ssg0_14, ssf_9, ssg1_14, \
                         psd0_9, psd1_9, psf_16, psf_17, psf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_2 * pc_z[k] * psf_16[k];

        t_27[k] = f_3 * psd0_9[k]
                  - f_4 * psd1_9[k]
                  + f_2 * pc_z[k] * psf_17[k];

        t_28[k] = f_5 * ssf_9[k]
                  + f_2 * pc_y[k] * psf_19[k];

        t_29[k] = pa_y[k] * ssg0_14[k]
                  - f_1 * pc_y[k] * ssg1_14[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_z, pc_x, pc_y, pc_z, ssg0_0, ssg0_3, \
                         ssg1_0, ssg1_3, psd0_14, psd1_14, psf_20, \
                         psf_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pa_z[k] * ssg0_0[k]
                  - f_1 * pc_z[k] * ssg1_0[k];

        t_31[k] = f_2 * pc_y[k] * psf_20[k];

        t_32[k] = f_6 * psd0_14[k]
                  - f_7 * psd1_14[k]
                  + f_2 * pc_x[k] * psf_22[k];

        t_33[k] = pa_z[k] * ssg0_3[k]
                  - f_1 * pc_z[k] * ssg1_3[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, t_38, t_39, pc_x, pc_y, psd0_17, psd1_17, \
                         psf_22, psf_25, psf_26, psf_27, psf_28, \
                         psf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_2 * pc_y[k] * psf_22[k];

        t_35[k] = f_3 * psd0_17[k]
                  - f_4 * psd1_17[k]
                  + f_2 * pc_x[k] * psf_25[k];

        t_36[k] = f_2 * pc_x[k] * psf_26[k];

        t_37[k] = f_2 * pc_x[k] * psf_27[k];

        t_38[k] = f_2 * pc_x[k] * psf_28[k];

        t_39[k] = f_2 * pc_x[k] * psf_29[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_z, pc_y, pc_z, ssg0_10, ssg1_10, psd0_16, \
                         psd0_17, psd1_16, psd1_17, psf_27, psf_28, \
                         psf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pa_z[k] * ssg0_10[k]
                  - f_1 * pc_z[k] * ssg1_10[k];

        t_41[k] = f_6 * psd0_16[k]
                  - f_7 * psd1_16[k]
                  + f_2 * pc_y[k] * psf_27[k];

        t_42[k] = f_3 * psd0_17[k]
                  - f_4 * psd1_17[k]
                  + f_2 * pc_y[k] * psf_28[k];

        t_43[k] = f_2 * pc_y[k] * psf_29[k];
    }

#pragma omp simd aligned(t_44, pa_z, pc_z, ssg0_14, ssf_9, ssg1_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = pa_z[k] * ssg0_14[k]
                  + f_0 * ssf_9[k]
                  - f_1 * pc_z[k] * ssg1_14[k];
    }
}

}  // namespace simdt3ceri
