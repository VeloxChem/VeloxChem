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


#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_psf_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t ssf0, const size_t ssd,
                                                   const size_t ssf1, const size_t psp0,
                                                   const size_t psp1, const size_t psd,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
    const auto f_1 = gamma / q;
    const auto f_2 = p / q;
    const auto f_3 = 0.5 / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ssf0_0 = buffer.data(ssf0 + 0);
    const auto *ssf0_6 = buffer.data(ssf0 + 6);
    const auto *ssf0_9 = buffer.data(ssf0 + 9);

    const auto *ssd_0 = buffer.data(ssd + 0);
    const auto *ssd_3 = buffer.data(ssd + 3);
    const auto *ssd_5 = buffer.data(ssd + 5);

    const auto *ssf1_0 = buffer.data(ssf1 + 0);
    const auto *ssf1_6 = buffer.data(ssf1 + 6);
    const auto *ssf1_9 = buffer.data(ssf1 + 9);

    const auto *psp0_4 = buffer.data(psp0 + 4);
    const auto *psp0_8 = buffer.data(psp0 + 8);

    const auto *psp1_4 = buffer.data(psp1 + 4);
    const auto *psp1_8 = buffer.data(psp1 + 8);

    const auto *psd_0 = buffer.data(psd + 0);
    const auto *psd_2 = buffer.data(psd + 2);
    const auto *psd_3 = buffer.data(psd + 3);
    const auto *psd_5 = buffer.data(psd + 5);
    const auto *psd_6 = buffer.data(psd + 6);
    const auto *psd_7 = buffer.data(psd + 7);
    const auto *psd_9 = buffer.data(psd + 9);
    const auto *psd_10 = buffer.data(psd + 10);
    const auto *psd_11 = buffer.data(psd + 11);
    const auto *psd_12 = buffer.data(psd + 12);
    const auto *psd_14 = buffer.data(psd + 14);
    const auto *psd_15 = buffer.data(psd + 15);
    const auto *psd_16 = buffer.data(psd + 16);
    const auto *psd_17 = buffer.data(psd + 17);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pc_x, pc_y, pc_z, ssf0_0, ssd_0, \
                         ssd_3, ssf1_0, psd_0, psd_2, psd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = pa_x[k] * ssf0_0[k]
                 + f_0 * ssd_0[k]
                 - f_1 * pc_x[k] * ssf1_0[k];

        t_1[k] = f_2 * pc_y[k] * psd_0[k];

        t_2[k] = f_2 * pc_z[k] * psd_0[k];

        t_3[k] = f_3 * ssd_3[k]
                 + f_2 * pc_x[k] * psd_3[k];

        t_4[k] = f_2 * pc_y[k] * psd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_x, pc_x, pc_y, pc_z, ssf0_6, ssf0_9, \
                         ssd_5, ssf1_6, ssf1_9, psd_3, psd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_3 * ssd_5[k]
                 + f_2 * pc_x[k] * psd_5[k];

        t_6[k] = pa_x[k] * ssf0_6[k]
                 - f_1 * pc_x[k] * ssf1_6[k];

        t_7[k] = f_2 * pc_z[k] * psd_3[k];

        t_8[k] = f_2 * pc_y[k] * psd_5[k];

        t_9[k] = pa_x[k] * ssf0_9[k]
                 - f_1 * pc_x[k] * ssf1_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_y, pc_x, pc_y, pc_z, ssf0_0, ssf1_0, \
                         psp0_4, psp1_4, psd_6, psd_7, psd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_y[k] * ssf0_0[k]
                  - f_1 * pc_y[k] * ssf1_0[k];

        t_11[k] = f_4 * psp0_4[k]
                  - f_5 * psp1_4[k]
                  + f_2 * pc_x[k] * psd_7[k];

        t_12[k] = f_2 * pc_z[k] * psd_6[k];

        t_13[k] = f_2 * pc_x[k] * psd_9[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, t_18, pa_y, pc_x, pc_y, pc_z, ssf0_6, ssd_3, \
                         ssd_5, ssf1_6, psd_9, psd_10, psd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_2 * pc_x[k] * psd_10[k];

        t_15[k] = f_2 * pc_x[k] * psd_11[k];

        t_16[k] = pa_y[k] * ssf0_6[k]
                  + f_0 * ssd_3[k]
                  - f_1 * pc_y[k] * ssf1_6[k];

        t_17[k] = f_2 * pc_z[k] * psd_9[k];

        t_18[k] = f_3 * ssd_5[k]
                  + f_2 * pc_y[k] * psd_11[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_y, pa_z, pc_y, pc_z, ssf0_0, ssf0_9, ssf1_0, \
                         ssf1_9, psd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pa_y[k] * ssf0_9[k]
                  - f_1 * pc_y[k] * ssf1_9[k];

        t_20[k] = pa_z[k] * ssf0_0[k]
                  - f_1 * pc_z[k] * ssf1_0[k];

        t_21[k] = f_2 * pc_y[k] * psd_12[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, pa_z, pc_x, pc_z, ssf0_6, ssf1_6, \
                         psp0_8, psp1_8, psd_14, psd_15, psd_16, \
                         psd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_4 * psp0_8[k]
                  - f_5 * psp1_8[k]
                  + f_2 * pc_x[k] * psd_14[k];

        t_23[k] = f_2 * pc_x[k] * psd_15[k];

        t_24[k] = f_2 * pc_x[k] * psd_16[k];

        t_25[k] = f_2 * pc_x[k] * psd_17[k];

        t_26[k] = pa_z[k] * ssf0_6[k]
                  - f_1 * pc_z[k] * ssf1_6[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_z, pc_y, pc_z, ssf0_9, ssd_5, ssf1_9, psp0_8, \
                         psp1_8, psd_16, psd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_4 * psp0_8[k]
                  - f_5 * psp1_8[k]
                  + f_2 * pc_y[k] * psd_16[k];

        t_28[k] = f_2 * pc_y[k] * psd_17[k];

        t_29[k] = pa_z[k] * ssf0_9[k]
                  + f_0 * ssd_5[k]
                  - f_1 * pc_z[k] * ssf1_9[k];
    }
}

}  // namespace simdt3ceri
