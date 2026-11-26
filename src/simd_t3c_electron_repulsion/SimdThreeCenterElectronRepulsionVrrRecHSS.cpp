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


#include "SimdThreeCenterElectronRepulsionVrrRecHSS.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_hss_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t fss0, const size_t fss1,
                                                   const size_t gss0, const size_t gss1,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 2.0 * gamma / (p * q);
    const auto f_2 = gamma / q;
    const auto f_3 = 1.0 / p;
    const auto f_4 = gamma / (p * q);
    const auto f_5 = 0.5 / p;
    const auto f_6 = 0.5 * gamma / (p * q);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *fss0_0 = buffer.data(fss0 + 0);
    const auto *fss0_3 = buffer.data(fss0 + 3);
    const auto *fss0_5 = buffer.data(fss0 + 5);
    const auto *fss0_6 = buffer.data(fss0 + 6);
    const auto *fss0_8 = buffer.data(fss0 + 8);
    const auto *fss0_9 = buffer.data(fss0 + 9);

    const auto *fss1_0 = buffer.data(fss1 + 0);
    const auto *fss1_3 = buffer.data(fss1 + 3);
    const auto *fss1_5 = buffer.data(fss1 + 5);
    const auto *fss1_6 = buffer.data(fss1 + 6);
    const auto *fss1_8 = buffer.data(fss1 + 8);
    const auto *fss1_9 = buffer.data(fss1 + 9);

    const auto *gss0_0 = buffer.data(gss0 + 0);
    const auto *gss0_2 = buffer.data(gss0 + 2);
    const auto *gss0_3 = buffer.data(gss0 + 3);
    const auto *gss0_5 = buffer.data(gss0 + 5);
    const auto *gss0_6 = buffer.data(gss0 + 6);
    const auto *gss0_9 = buffer.data(gss0 + 9);
    const auto *gss0_10 = buffer.data(gss0 + 10);
    const auto *gss0_11 = buffer.data(gss0 + 11);
    const auto *gss0_12 = buffer.data(gss0 + 12);
    const auto *gss0_13 = buffer.data(gss0 + 13);
    const auto *gss0_14 = buffer.data(gss0 + 14);

    const auto *gss1_0 = buffer.data(gss1 + 0);
    const auto *gss1_2 = buffer.data(gss1 + 2);
    const auto *gss1_3 = buffer.data(gss1 + 3);
    const auto *gss1_5 = buffer.data(gss1 + 5);
    const auto *gss1_6 = buffer.data(gss1 + 6);
    const auto *gss1_9 = buffer.data(gss1 + 9);
    const auto *gss1_10 = buffer.data(gss1 + 10);
    const auto *gss1_11 = buffer.data(gss1 + 11);
    const auto *gss1_12 = buffer.data(gss1 + 12);
    const auto *gss1_13 = buffer.data(gss1 + 13);
    const auto *gss1_14 = buffer.data(gss1 + 14);

#pragma omp simd aligned(t_0, t_1, t_2, pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, fss0_0, fss1_0, \
                         gss0_0, gss1_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fss0_0[k]
                 - f_1 * fss1_0[k]
                 + pa_x[k] * gss0_0[k]
                 - f_2 * pc_x[k] * gss1_0[k];

        t_1[k] = pa_y[k] * gss0_0[k]
                 - f_2 * pc_y[k] * gss1_0[k];

        t_2[k] = pa_z[k] * gss0_0[k]
                 - f_2 * pc_z[k] * gss1_0[k];
    }

#pragma omp simd aligned(t_3, t_4, pa_x, pa_y, pc_x, pc_y, fss0_3, fss1_3, gss0_2, gss0_3, \
                         gss1_2, gss1_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_3 * fss0_3[k]
                 - f_4 * fss1_3[k]
                 + pa_x[k] * gss0_3[k]
                 - f_2 * pc_x[k] * gss1_3[k];

        t_4[k] = pa_y[k] * gss0_2[k]
                 - f_2 * pc_y[k] * gss1_2[k];
    }

#pragma omp simd aligned(t_5, t_6, pa_x, pc_x, fss0_5, fss0_6, fss1_5, fss1_6, gss0_5, gss0_6, \
                         gss1_5, gss1_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_3 * fss0_5[k]
                 - f_4 * fss1_5[k]
                 + pa_x[k] * gss0_5[k]
                 - f_2 * pc_x[k] * gss1_5[k];

        t_6[k] = f_5 * fss0_6[k]
                 - f_6 * fss1_6[k]
                 + pa_x[k] * gss0_6[k]
                 - f_2 * pc_x[k] * gss1_6[k];
    }

#pragma omp simd aligned(t_7, t_8, pa_y, pa_z, pc_y, pc_z, gss0_3, gss0_5, gss1_3, \
                         gss1_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = pa_z[k] * gss0_3[k]
                 - f_2 * pc_z[k] * gss1_3[k];

        t_8[k] = pa_y[k] * gss0_5[k]
                 - f_2 * pc_y[k] * gss1_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pa_x, pc_x, fss0_9, fss1_9, gss0_9, gss0_10, \
                         gss0_11, gss0_12, gss1_9, gss1_10, gss1_11, \
                         gss1_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * fss0_9[k]
                 - f_6 * fss1_9[k]
                 + pa_x[k] * gss0_9[k]
                 - f_2 * pc_x[k] * gss1_9[k];

        t_10[k] = pa_x[k] * gss0_10[k]
                  - f_2 * pc_x[k] * gss1_10[k];

        t_11[k] = pa_x[k] * gss0_11[k]
                  - f_2 * pc_x[k] * gss1_11[k];

        t_12[k] = pa_x[k] * gss0_12[k]
                  - f_2 * pc_x[k] * gss1_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_x, pa_y, pc_x, pc_y, fss0_6, fss1_6, gss0_10, \
                         gss0_13, gss0_14, gss1_10, gss1_13, gss1_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_x[k] * gss0_13[k]
                  - f_2 * pc_x[k] * gss1_13[k];

        t_14[k] = pa_x[k] * gss0_14[k]
                  - f_2 * pc_x[k] * gss1_14[k];

        t_15[k] = f_0 * fss0_6[k]
                  - f_1 * fss1_6[k]
                  + pa_y[k] * gss0_10[k]
                  - f_2 * pc_y[k] * gss1_10[k];
    }

#pragma omp simd aligned(t_16, t_17, pa_y, pa_z, pc_y, pc_z, fss0_8, fss1_8, gss0_10, gss0_12, \
                         gss1_10, gss1_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pa_z[k] * gss0_10[k]
                  - f_2 * pc_z[k] * gss1_10[k];

        t_17[k] = f_3 * fss0_8[k]
                  - f_4 * fss1_8[k]
                  + pa_y[k] * gss0_12[k]
                  - f_2 * pc_y[k] * gss1_12[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_y, pa_z, pc_y, pc_z, fss0_9, fss1_9, gss0_13, \
                         gss0_14, gss1_13, gss1_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * fss0_9[k]
                  - f_6 * fss1_9[k]
                  + pa_y[k] * gss0_13[k]
                  - f_2 * pc_y[k] * gss1_13[k];

        t_19[k] = pa_y[k] * gss0_14[k]
                  - f_2 * pc_y[k] * gss1_14[k];

        t_20[k] = f_0 * fss0_9[k]
                  - f_1 * fss1_9[k]
                  + pa_z[k] * gss0_14[k]
                  - f_2 * pc_z[k] * gss1_14[k];
    }
}

}  // namespace simdt3ceri
