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


#include "SimdThreeCenterElectronRepulsionVrrRecISS.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_iss_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t gss0, const size_t gss1,
                                                   const size_t hss0, const size_t hss1,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 2.5 * gamma / (p * q);
    const auto f_2 = gamma / q;
    const auto f_3 = 1.5 / p;
    const auto f_4 = 1.5 * gamma / (p * q);
    const auto f_5 = 1.0 / p;
    const auto f_6 = gamma / (p * q);
    const auto f_7 = 0.5 / p;
    const auto f_8 = 0.5 * gamma / (p * q);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *gss0_0 = buffer.data(gss0 + 0);
    const auto *gss0_3 = buffer.data(gss0 + 3);
    const auto *gss0_5 = buffer.data(gss0 + 5);
    const auto *gss0_6 = buffer.data(gss0 + 6);
    const auto *gss0_9 = buffer.data(gss0 + 9);
    const auto *gss0_10 = buffer.data(gss0 + 10);
    const auto *gss0_12 = buffer.data(gss0 + 12);
    const auto *gss0_13 = buffer.data(gss0 + 13);
    const auto *gss0_14 = buffer.data(gss0 + 14);

    const auto *gss1_0 = buffer.data(gss1 + 0);
    const auto *gss1_3 = buffer.data(gss1 + 3);
    const auto *gss1_5 = buffer.data(gss1 + 5);
    const auto *gss1_6 = buffer.data(gss1 + 6);
    const auto *gss1_9 = buffer.data(gss1 + 9);
    const auto *gss1_10 = buffer.data(gss1 + 10);
    const auto *gss1_12 = buffer.data(gss1 + 12);
    const auto *gss1_13 = buffer.data(gss1 + 13);
    const auto *gss1_14 = buffer.data(gss1 + 14);

    const auto *hss0_0 = buffer.data(hss0 + 0);
    const auto *hss0_2 = buffer.data(hss0 + 2);
    const auto *hss0_3 = buffer.data(hss0 + 3);
    const auto *hss0_5 = buffer.data(hss0 + 5);
    const auto *hss0_6 = buffer.data(hss0 + 6);
    const auto *hss0_9 = buffer.data(hss0 + 9);
    const auto *hss0_10 = buffer.data(hss0 + 10);
    const auto *hss0_12 = buffer.data(hss0 + 12);
    const auto *hss0_14 = buffer.data(hss0 + 14);
    const auto *hss0_15 = buffer.data(hss0 + 15);
    const auto *hss0_16 = buffer.data(hss0 + 16);
    const auto *hss0_17 = buffer.data(hss0 + 17);
    const auto *hss0_18 = buffer.data(hss0 + 18);
    const auto *hss0_19 = buffer.data(hss0 + 19);
    const auto *hss0_20 = buffer.data(hss0 + 20);

    const auto *hss1_0 = buffer.data(hss1 + 0);
    const auto *hss1_2 = buffer.data(hss1 + 2);
    const auto *hss1_3 = buffer.data(hss1 + 3);
    const auto *hss1_5 = buffer.data(hss1 + 5);
    const auto *hss1_6 = buffer.data(hss1 + 6);
    const auto *hss1_9 = buffer.data(hss1 + 9);
    const auto *hss1_10 = buffer.data(hss1 + 10);
    const auto *hss1_12 = buffer.data(hss1 + 12);
    const auto *hss1_14 = buffer.data(hss1 + 14);
    const auto *hss1_15 = buffer.data(hss1 + 15);
    const auto *hss1_16 = buffer.data(hss1 + 16);
    const auto *hss1_17 = buffer.data(hss1 + 17);
    const auto *hss1_18 = buffer.data(hss1 + 18);
    const auto *hss1_19 = buffer.data(hss1 + 19);
    const auto *hss1_20 = buffer.data(hss1 + 20);

#pragma omp simd aligned(t_0, t_1, t_2, pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, gss0_0, gss1_0, \
                         hss0_0, hss1_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gss0_0[k]
                 - f_1 * gss1_0[k]
                 + pa_x[k] * hss0_0[k]
                 - f_2 * pc_x[k] * hss1_0[k];

        t_1[k] = pa_y[k] * hss0_0[k]
                 - f_2 * pc_y[k] * hss1_0[k];

        t_2[k] = pa_z[k] * hss0_0[k]
                 - f_2 * pc_z[k] * hss1_0[k];
    }

#pragma omp simd aligned(t_3, t_4, pa_x, pa_y, pc_x, pc_y, gss0_3, gss1_3, hss0_2, hss0_3, \
                         hss1_2, hss1_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_3 * gss0_3[k]
                 - f_4 * gss1_3[k]
                 + pa_x[k] * hss0_3[k]
                 - f_2 * pc_x[k] * hss1_3[k];

        t_4[k] = pa_y[k] * hss0_2[k]
                 - f_2 * pc_y[k] * hss1_2[k];
    }

#pragma omp simd aligned(t_5, t_6, pa_x, pc_x, gss0_5, gss0_6, gss1_5, gss1_6, hss0_5, hss0_6, \
                         hss1_5, hss1_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_3 * gss0_5[k]
                 - f_4 * gss1_5[k]
                 + pa_x[k] * hss0_5[k]
                 - f_2 * pc_x[k] * hss1_5[k];

        t_6[k] = f_5 * gss0_6[k]
                 - f_6 * gss1_6[k]
                 + pa_x[k] * hss0_6[k]
                 - f_2 * pc_x[k] * hss1_6[k];
    }

#pragma omp simd aligned(t_7, t_8, pa_y, pa_z, pc_y, pc_z, hss0_3, hss0_5, hss1_3, \
                         hss1_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = pa_z[k] * hss0_3[k]
                 - f_2 * pc_z[k] * hss1_3[k];

        t_8[k] = pa_y[k] * hss0_5[k]
                 - f_2 * pc_y[k] * hss1_5[k];
    }

#pragma omp simd aligned(t_9, t_10, pa_x, pc_x, gss0_9, gss0_10, gss1_9, gss1_10, hss0_9, \
                         hss0_10, hss1_9, hss1_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * gss0_9[k]
                 - f_6 * gss1_9[k]
                 + pa_x[k] * hss0_9[k]
                 - f_2 * pc_x[k] * hss1_9[k];

        t_10[k] = f_7 * gss0_10[k]
                  - f_8 * gss1_10[k]
                  + pa_x[k] * hss0_10[k]
                  - f_2 * pc_x[k] * hss1_10[k];
    }

#pragma omp simd aligned(t_11, t_12, pa_x, pa_z, pc_x, pc_z, gss0_12, gss1_12, hss0_6, \
                         hss0_12, hss1_6, hss1_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pa_z[k] * hss0_6[k]
                  - f_2 * pc_z[k] * hss1_6[k];

        t_12[k] = f_7 * gss0_12[k]
                  - f_8 * gss1_12[k]
                  + pa_x[k] * hss0_12[k]
                  - f_2 * pc_x[k] * hss1_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_x, pa_y, pc_x, pc_y, gss0_14, gss1_14, hss0_9, \
                         hss0_14, hss0_15, hss1_9, hss1_14, hss1_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_y[k] * hss0_9[k]
                  - f_2 * pc_y[k] * hss1_9[k];

        t_14[k] = f_7 * gss0_14[k]
                  - f_8 * gss1_14[k]
                  + pa_x[k] * hss0_14[k]
                  - f_2 * pc_x[k] * hss1_14[k];

        t_15[k] = pa_x[k] * hss0_15[k]
                  - f_2 * pc_x[k] * hss1_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pc_x, hss0_16, hss0_17, hss0_18, \
                         hss0_19, hss1_16, hss1_17, hss1_18, hss1_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pa_x[k] * hss0_16[k]
                  - f_2 * pc_x[k] * hss1_16[k];

        t_17[k] = pa_x[k] * hss0_17[k]
                  - f_2 * pc_x[k] * hss1_17[k];

        t_18[k] = pa_x[k] * hss0_18[k]
                  - f_2 * pc_x[k] * hss1_18[k];

        t_19[k] = pa_x[k] * hss0_19[k]
                  - f_2 * pc_x[k] * hss1_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, gss0_10, \
                         gss1_10, hss0_15, hss0_20, hss1_15, hss1_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_x[k] * hss0_20[k]
                  - f_2 * pc_x[k] * hss1_20[k];

        t_21[k] = f_0 * gss0_10[k]
                  - f_1 * gss1_10[k]
                  + pa_y[k] * hss0_15[k]
                  - f_2 * pc_y[k] * hss1_15[k];

        t_22[k] = pa_z[k] * hss0_15[k]
                  - f_2 * pc_z[k] * hss1_15[k];
    }

#pragma omp simd aligned(t_23, t_24, pa_y, pc_y, gss0_12, gss0_13, gss1_12, gss1_13, hss0_17, \
                         hss0_18, hss1_17, hss1_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_3 * gss0_12[k]
                  - f_4 * gss1_12[k]
                  + pa_y[k] * hss0_17[k]
                  - f_2 * pc_y[k] * hss1_17[k];

        t_24[k] = f_5 * gss0_13[k]
                  - f_6 * gss1_13[k]
                  + pa_y[k] * hss0_18[k]
                  - f_2 * pc_y[k] * hss1_18[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pa_y, pa_z, pc_y, pc_z, gss0_14, gss1_14, hss0_19, \
                         hss0_20, hss1_19, hss1_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_7 * gss0_14[k]
                  - f_8 * gss1_14[k]
                  + pa_y[k] * hss0_19[k]
                  - f_2 * pc_y[k] * hss1_19[k];

        t_26[k] = pa_y[k] * hss0_20[k]
                  - f_2 * pc_y[k] * hss1_20[k];

        t_27[k] = f_0 * gss0_14[k]
                  - f_1 * gss1_14[k]
                  + pa_z[k] * hss0_20[k]
                  - f_2 * pc_z[k] * hss1_20[k];
    }
}

}  // namespace simdt3ceri
