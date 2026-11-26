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


#include "SimdThreeCenterElectronRepulsionVrrRecFPS.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_fps_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pb,
                                                   const size_t pc, const size_t dss0,
                                                   const size_t dss1, const size_t dps0,
                                                   const size_t dps1, const size_t fss0,
                                                   const size_t fss1, const size_t ncols,
                                                   const double gamma, const double p,
                                                   const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 1.5 * gamma / (p * q);
    const auto f_2 = gamma / q;
    const auto f_3 = 0.5 / p;
    const auto f_4 = 0.5 * gamma / (p * q);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *dss0_0 = buffer.data(dss0 + 0);
    const auto *dss0_3 = buffer.data(dss0 + 3);
    const auto *dss0_4 = buffer.data(dss0 + 4);
    const auto *dss0_5 = buffer.data(dss0 + 5);

    const auto *dss1_0 = buffer.data(dss1 + 0);
    const auto *dss1_3 = buffer.data(dss1 + 3);
    const auto *dss1_4 = buffer.data(dss1 + 4);
    const auto *dss1_5 = buffer.data(dss1 + 5);

    const auto *dps0_0 = buffer.data(dps0 + 0);
    const auto *dps0_1 = buffer.data(dps0 + 1);
    const auto *dps0_2 = buffer.data(dps0 + 2);
    const auto *dps0_9 = buffer.data(dps0 + 9);
    const auto *dps0_10 = buffer.data(dps0 + 10);
    const auto *dps0_11 = buffer.data(dps0 + 11);
    const auto *dps0_12 = buffer.data(dps0 + 12);
    const auto *dps0_13 = buffer.data(dps0 + 13);
    const auto *dps0_14 = buffer.data(dps0 + 14);
    const auto *dps0_15 = buffer.data(dps0 + 15);
    const auto *dps0_16 = buffer.data(dps0 + 16);
    const auto *dps0_17 = buffer.data(dps0 + 17);

    const auto *dps1_0 = buffer.data(dps1 + 0);
    const auto *dps1_1 = buffer.data(dps1 + 1);
    const auto *dps1_2 = buffer.data(dps1 + 2);
    const auto *dps1_9 = buffer.data(dps1 + 9);
    const auto *dps1_10 = buffer.data(dps1 + 10);
    const auto *dps1_11 = buffer.data(dps1 + 11);
    const auto *dps1_12 = buffer.data(dps1 + 12);
    const auto *dps1_13 = buffer.data(dps1 + 13);
    const auto *dps1_14 = buffer.data(dps1 + 14);
    const auto *dps1_15 = buffer.data(dps1 + 15);
    const auto *dps1_16 = buffer.data(dps1 + 16);
    const auto *dps1_17 = buffer.data(dps1 + 17);

    const auto *fss0_0 = buffer.data(fss0 + 0);
    const auto *fss0_6 = buffer.data(fss0 + 6);
    const auto *fss0_9 = buffer.data(fss0 + 9);

    const auto *fss1_0 = buffer.data(fss1 + 0);
    const auto *fss1_6 = buffer.data(fss1 + 6);
    const auto *fss1_9 = buffer.data(fss1 + 9);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, pc_x, pc_y, pc_z, dss0_0, dss1_0, \
                         fss0_0, fss1_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dss0_0[k]
                 - f_1 * dss1_0[k]
                 + pb_x[k] * fss0_0[k]
                 - f_2 * pc_x[k] * fss1_0[k];

        t_1[k] = pb_y[k] * fss0_0[k]
                 - f_2 * pc_y[k] * fss1_0[k];

        t_2[k] = pb_z[k] * fss0_0[k]
                 - f_2 * pc_z[k] * fss1_0[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, t_6, pa_y, pa_z, pc_y, pc_z, dss0_0, dss1_0, dps0_0, \
                         dps0_1, dps0_2, dps1_0, dps1_1, dps1_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = pa_y[k] * dps0_0[k]
                 - f_2 * pc_y[k] * dps1_0[k];

        t_4[k] = f_3 * dss0_0[k]
                 - f_4 * dss1_0[k]
                 + pa_y[k] * dps0_1[k]
                 - f_2 * pc_y[k] * dps1_1[k];

        t_5[k] = pa_y[k] * dps0_2[k]
                 - f_2 * pc_y[k] * dps1_2[k];

        t_6[k] = pa_z[k] * dps0_0[k]
                 - f_2 * pc_z[k] * dps1_0[k];
    }

#pragma omp simd aligned(t_7, t_8, pa_z, pc_z, dss0_0, dss1_0, dps0_1, dps0_2, dps1_1, \
                         dps1_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = pa_z[k] * dps0_1[k]
                 - f_2 * pc_z[k] * dps1_1[k];

        t_8[k] = f_3 * dss0_0[k]
                 - f_4 * dss1_0[k]
                 + pa_z[k] * dps0_2[k]
                 - f_2 * pc_z[k] * dps1_2[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pc_x, dss0_3, dss1_3, dps0_9, dps0_10, \
                         dps0_11, dps1_9, dps1_10, dps1_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_3 * dss0_3[k]
                 - f_4 * dss1_3[k]
                 + pa_x[k] * dps0_9[k]
                 - f_2 * pc_x[k] * dps1_9[k];

        t_10[k] = pa_x[k] * dps0_10[k]
                  - f_2 * pc_x[k] * dps1_10[k];

        t_11[k] = pa_x[k] * dps0_11[k]
                  - f_2 * pc_x[k] * dps1_11[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pc_x, dss0_4, dss1_4, dps0_12, dps0_13, \
                         dps0_14, dps1_12, dps1_13, dps1_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_3 * dss0_4[k]
                  - f_4 * dss1_4[k]
                  + pa_x[k] * dps0_12[k]
                  - f_2 * pc_x[k] * dps1_12[k];

        t_13[k] = pa_x[k] * dps0_13[k]
                  - f_2 * pc_x[k] * dps1_13[k];

        t_14[k] = pa_x[k] * dps0_14[k]
                  - f_2 * pc_x[k] * dps1_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_x, pc_x, dss0_5, dss1_5, dps0_15, dps0_16, \
                         dps0_17, dps1_15, dps1_16, dps1_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_3 * dss0_5[k]
                  - f_4 * dss1_5[k]
                  + pa_x[k] * dps0_15[k]
                  - f_2 * pc_x[k] * dps1_15[k];

        t_16[k] = pa_x[k] * dps0_16[k]
                  - f_2 * pc_x[k] * dps1_16[k];

        t_17[k] = pa_x[k] * dps0_17[k]
                  - f_2 * pc_x[k] * dps1_17[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pb_x, pb_y, pb_z, pc_x, pc_y, pc_z, dss0_3, dss1_3, \
                         fss0_6, fss1_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = pb_x[k] * fss0_6[k]
                  - f_2 * pc_x[k] * fss1_6[k];

        t_19[k] = f_0 * dss0_3[k]
                  - f_1 * dss1_3[k]
                  + pb_y[k] * fss0_6[k]
                  - f_2 * pc_y[k] * fss1_6[k];

        t_20[k] = pb_z[k] * fss0_6[k]
                  - f_2 * pc_z[k] * fss1_6[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pa_z, pc_z, dss0_3, dss1_3, dps0_9, dps0_10, \
                         dps0_11, dps1_9, dps1_10, dps1_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = pa_z[k] * dps0_9[k]
                  - f_2 * pc_z[k] * dps1_9[k];

        t_22[k] = pa_z[k] * dps0_10[k]
                  - f_2 * pc_z[k] * dps1_10[k];

        t_23[k] = f_3 * dss0_3[k]
                  - f_4 * dss1_3[k]
                  + pa_z[k] * dps0_11[k]
                  - f_2 * pc_z[k] * dps1_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pa_y, pc_y, dss0_5, dss1_5, dps0_15, dps0_16, \
                         dps0_17, dps1_15, dps1_16, dps1_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pa_y[k] * dps0_15[k]
                  - f_2 * pc_y[k] * dps1_15[k];

        t_25[k] = f_3 * dss0_5[k]
                  - f_4 * dss1_5[k]
                  + pa_y[k] * dps0_16[k]
                  - f_2 * pc_y[k] * dps1_16[k];

        t_26[k] = pa_y[k] * dps0_17[k]
                  - f_2 * pc_y[k] * dps1_17[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pb_x, pb_y, pb_z, pc_x, pc_y, pc_z, dss0_5, dss1_5, \
                         fss0_9, fss1_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = pb_x[k] * fss0_9[k]
                  - f_2 * pc_x[k] * fss1_9[k];

        t_28[k] = pb_y[k] * fss0_9[k]
                  - f_2 * pc_y[k] * fss1_9[k];

        t_29[k] = f_0 * dss0_5[k]
                  - f_1 * dss1_5[k]
                  + pb_z[k] * fss0_9[k]
                  - f_2 * pc_z[k] * fss1_9[k];
    }
}

}  // namespace simdt3ceri
