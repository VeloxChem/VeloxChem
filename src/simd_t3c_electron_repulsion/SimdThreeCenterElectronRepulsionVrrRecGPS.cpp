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


#include "SimdThreeCenterElectronRepulsionVrrRecGPS.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_gps_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pb,
                                                   const size_t pc, const size_t dps0,
                                                   const size_t dps1, const size_t fss0,
                                                   const size_t fss1, const size_t fps0,
                                                   const size_t fps1, const size_t gss0,
                                                   const size_t gss1, const size_t ncols,
                                                   const double gamma, const double p,
                                                   const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 2.0 * gamma / (p * q);
    const auto f_2 = gamma / q;
    const auto f_3 = 0.5 / p;
    const auto f_4 = 0.5 * gamma / (p * q);
    const auto f_5 = 1.0 / p;
    const auto f_6 = gamma / (p * q);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *dps0_10 = buffer.data(dps0 + 10);
    const auto *dps0_17 = buffer.data(dps0 + 17);

    const auto *dps1_10 = buffer.data(dps1 + 10);
    const auto *dps1_17 = buffer.data(dps1 + 17);

    const auto *fss0_0 = buffer.data(fss0 + 0);
    const auto *fss0_2 = buffer.data(fss0 + 2);
    const auto *fss0_3 = buffer.data(fss0 + 3);
    const auto *fss0_5 = buffer.data(fss0 + 5);
    const auto *fss0_6 = buffer.data(fss0 + 6);
    const auto *fss0_7 = buffer.data(fss0 + 7);
    const auto *fss0_8 = buffer.data(fss0 + 8);
    const auto *fss0_9 = buffer.data(fss0 + 9);

    const auto *fss1_0 = buffer.data(fss1 + 0);
    const auto *fss1_2 = buffer.data(fss1 + 2);
    const auto *fss1_3 = buffer.data(fss1 + 3);
    const auto *fss1_5 = buffer.data(fss1 + 5);
    const auto *fss1_6 = buffer.data(fss1 + 6);
    const auto *fss1_7 = buffer.data(fss1 + 7);
    const auto *fss1_8 = buffer.data(fss1 + 8);
    const auto *fss1_9 = buffer.data(fss1 + 9);

    const auto *fps0_0 = buffer.data(fps0 + 0);
    const auto *fps0_1 = buffer.data(fps0 + 1);
    const auto *fps0_2 = buffer.data(fps0 + 2);
    const auto *fps0_6 = buffer.data(fps0 + 6);
    const auto *fps0_7 = buffer.data(fps0 + 7);
    const auto *fps0_8 = buffer.data(fps0 + 8);
    const auto *fps0_10 = buffer.data(fps0 + 10);
    const auto *fps0_17 = buffer.data(fps0 + 17);
    const auto *fps0_18 = buffer.data(fps0 + 18);
    const auto *fps0_19 = buffer.data(fps0 + 19);
    const auto *fps0_20 = buffer.data(fps0 + 20);
    const auto *fps0_21 = buffer.data(fps0 + 21);
    const auto *fps0_22 = buffer.data(fps0 + 22);
    const auto *fps0_23 = buffer.data(fps0 + 23);
    const auto *fps0_24 = buffer.data(fps0 + 24);
    const auto *fps0_25 = buffer.data(fps0 + 25);
    const auto *fps0_26 = buffer.data(fps0 + 26);
    const auto *fps0_27 = buffer.data(fps0 + 27);
    const auto *fps0_28 = buffer.data(fps0 + 28);
    const auto *fps0_29 = buffer.data(fps0 + 29);

    const auto *fps1_0 = buffer.data(fps1 + 0);
    const auto *fps1_1 = buffer.data(fps1 + 1);
    const auto *fps1_2 = buffer.data(fps1 + 2);
    const auto *fps1_6 = buffer.data(fps1 + 6);
    const auto *fps1_7 = buffer.data(fps1 + 7);
    const auto *fps1_8 = buffer.data(fps1 + 8);
    const auto *fps1_10 = buffer.data(fps1 + 10);
    const auto *fps1_17 = buffer.data(fps1 + 17);
    const auto *fps1_18 = buffer.data(fps1 + 18);
    const auto *fps1_19 = buffer.data(fps1 + 19);
    const auto *fps1_20 = buffer.data(fps1 + 20);
    const auto *fps1_21 = buffer.data(fps1 + 21);
    const auto *fps1_22 = buffer.data(fps1 + 22);
    const auto *fps1_23 = buffer.data(fps1 + 23);
    const auto *fps1_24 = buffer.data(fps1 + 24);
    const auto *fps1_25 = buffer.data(fps1 + 25);
    const auto *fps1_26 = buffer.data(fps1 + 26);
    const auto *fps1_27 = buffer.data(fps1 + 27);
    const auto *fps1_28 = buffer.data(fps1 + 28);
    const auto *fps1_29 = buffer.data(fps1 + 29);

    const auto *gss0_0 = buffer.data(gss0 + 0);
    const auto *gss0_3 = buffer.data(gss0 + 3);
    const auto *gss0_5 = buffer.data(gss0 + 5);
    const auto *gss0_10 = buffer.data(gss0 + 10);
    const auto *gss0_12 = buffer.data(gss0 + 12);
    const auto *gss0_14 = buffer.data(gss0 + 14);

    const auto *gss1_0 = buffer.data(gss1 + 0);
    const auto *gss1_3 = buffer.data(gss1 + 3);
    const auto *gss1_5 = buffer.data(gss1 + 5);
    const auto *gss1_10 = buffer.data(gss1 + 10);
    const auto *gss1_12 = buffer.data(gss1 + 12);
    const auto *gss1_14 = buffer.data(gss1 + 14);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, pc_x, pc_y, pc_z, fss0_0, fss1_0, \
                         gss0_0, gss1_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fss0_0[k]
                 - f_1 * fss1_0[k]
                 + pb_x[k] * gss0_0[k]
                 - f_2 * pc_x[k] * gss1_0[k];

        t_1[k] = pb_y[k] * gss0_0[k]
                 - f_2 * pc_y[k] * gss1_0[k];

        t_2[k] = pb_z[k] * gss0_0[k]
                 - f_2 * pc_z[k] * gss1_0[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, t_6, pa_y, pa_z, pc_y, pc_z, fss0_0, fss1_0, fps0_0, \
                         fps0_1, fps0_2, fps1_0, fps1_1, fps1_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = pa_y[k] * fps0_0[k]
                 - f_2 * pc_y[k] * fps1_0[k];

        t_4[k] = f_3 * fss0_0[k]
                 - f_4 * fss1_0[k]
                 + pa_y[k] * fps0_1[k]
                 - f_2 * pc_y[k] * fps1_1[k];

        t_5[k] = pa_y[k] * fps0_2[k]
                 - f_2 * pc_y[k] * fps1_2[k];

        t_6[k] = pa_z[k] * fps0_0[k]
                 - f_2 * pc_z[k] * fps1_0[k];
    }

#pragma omp simd aligned(t_7, t_8, pa_z, pc_z, fss0_0, fss1_0, fps0_1, fps0_2, fps1_1, \
                         fps1_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = pa_z[k] * fps0_1[k]
                 - f_2 * pc_z[k] * fps1_1[k];

        t_8[k] = f_3 * fss0_0[k]
                 - f_4 * fss1_0[k]
                 + pa_z[k] * fps0_2[k]
                 - f_2 * pc_z[k] * fps1_2[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pb_x, pb_z, pc_x, pc_z, dps0_10, dps1_10, \
                         fss0_3, fss1_3, fps0_10, fps1_10, gss0_3, \
                         gss1_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * fss0_3[k]
                 - f_6 * fss1_3[k]
                 + pb_x[k] * gss0_3[k]
                 - f_2 * pc_x[k] * gss1_3[k];

        t_10[k] = f_3 * dps0_10[k]
                  - f_4 * dps1_10[k]
                  + pa_x[k] * fps0_10[k]
                  - f_2 * pc_x[k] * fps1_10[k];

        t_11[k] = pb_z[k] * gss0_3[k]
                  - f_2 * pc_z[k] * gss1_3[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_y, pc_y, fss0_2, fss1_2, fps0_6, fps0_7, fps0_8, \
                         fps1_6, fps1_7, fps1_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pa_y[k] * fps0_6[k]
                  - f_2 * pc_y[k] * fps1_6[k];

        t_13[k] = f_3 * fss0_2[k]
                  - f_4 * fss1_2[k]
                  + pa_y[k] * fps0_7[k]
                  - f_2 * pc_y[k] * fps1_7[k];

        t_14[k] = pa_y[k] * fps0_8[k]
                  - f_2 * pc_y[k] * fps1_8[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_x, pb_x, pb_y, pc_x, pc_y, dps0_17, dps1_17, \
                         fss0_5, fss1_5, fps0_17, fps1_17, gss0_5, \
                         gss1_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_5 * fss0_5[k]
                  - f_6 * fss1_5[k]
                  + pb_x[k] * gss0_5[k]
                  - f_2 * pc_x[k] * gss1_5[k];

        t_16[k] = pb_y[k] * gss0_5[k]
                  - f_2 * pc_y[k] * gss1_5[k];

        t_17[k] = f_3 * dps0_17[k]
                  - f_4 * dps1_17[k]
                  + pa_x[k] * fps0_17[k]
                  - f_2 * pc_x[k] * fps1_17[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_x, pc_x, fss0_6, fss1_6, fps0_18, fps0_19, \
                         fps0_20, fps1_18, fps1_19, fps1_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_3 * fss0_6[k]
                  - f_4 * fss1_6[k]
                  + pa_x[k] * fps0_18[k]
                  - f_2 * pc_x[k] * fps1_18[k];

        t_19[k] = pa_x[k] * fps0_19[k]
                  - f_2 * pc_x[k] * fps1_19[k];

        t_20[k] = pa_x[k] * fps0_20[k]
                  - f_2 * pc_x[k] * fps1_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pa_x, pc_x, fss0_7, fss1_7, fps0_21, fps0_22, \
                         fps0_23, fps1_21, fps1_22, fps1_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_3 * fss0_7[k]
                  - f_4 * fss1_7[k]
                  + pa_x[k] * fps0_21[k]
                  - f_2 * pc_x[k] * fps1_21[k];

        t_22[k] = pa_x[k] * fps0_22[k]
                  - f_2 * pc_x[k] * fps1_22[k];

        t_23[k] = pa_x[k] * fps0_23[k]
                  - f_2 * pc_x[k] * fps1_23[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pa_x, pc_x, fss0_8, fss1_8, fps0_24, fps0_25, \
                         fps0_26, fps1_24, fps1_25, fps1_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_3 * fss0_8[k]
                  - f_4 * fss1_8[k]
                  + pa_x[k] * fps0_24[k]
                  - f_2 * pc_x[k] * fps1_24[k];

        t_25[k] = pa_x[k] * fps0_25[k]
                  - f_2 * pc_x[k] * fps1_25[k];

        t_26[k] = pa_x[k] * fps0_26[k]
                  - f_2 * pc_x[k] * fps1_26[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_x, pc_x, fss0_9, fss1_9, fps0_27, fps0_28, \
                         fps0_29, fps1_27, fps1_28, fps1_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_3 * fss0_9[k]
                  - f_4 * fss1_9[k]
                  + pa_x[k] * fps0_27[k]
                  - f_2 * pc_x[k] * fps1_27[k];

        t_28[k] = pa_x[k] * fps0_28[k]
                  - f_2 * pc_x[k] * fps1_28[k];

        t_29[k] = pa_x[k] * fps0_29[k]
                  - f_2 * pc_x[k] * fps1_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pb_x, pb_y, pb_z, pc_x, pc_y, pc_z, fss0_6, fss1_6, \
                         gss0_10, gss1_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pb_x[k] * gss0_10[k]
                  - f_2 * pc_x[k] * gss1_10[k];

        t_31[k] = f_0 * fss0_6[k]
                  - f_1 * fss1_6[k]
                  + pb_y[k] * gss0_10[k]
                  - f_2 * pc_y[k] * gss1_10[k];

        t_32[k] = pb_z[k] * gss0_10[k]
                  - f_2 * pc_z[k] * gss1_10[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, pa_z, pc_z, fss0_6, fss1_6, fps0_18, fps0_19, \
                         fps0_20, fps1_18, fps1_19, fps1_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pa_z[k] * fps0_18[k]
                  - f_2 * pc_z[k] * fps1_18[k];

        t_34[k] = pa_z[k] * fps0_19[k]
                  - f_2 * pc_z[k] * fps1_19[k];

        t_35[k] = f_3 * fss0_6[k]
                  - f_4 * fss1_6[k]
                  + pa_z[k] * fps0_20[k]
                  - f_2 * pc_z[k] * fps1_20[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pa_y, pb_x, pb_y, pc_x, pc_y, dps0_17, dps1_17, \
                         fss0_8, fss1_8, fps0_26, fps1_26, gss0_12, \
                         gss1_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = pb_x[k] * gss0_12[k]
                  - f_2 * pc_x[k] * gss1_12[k];

        t_37[k] = f_5 * fss0_8[k]
                  - f_6 * fss1_8[k]
                  + pb_y[k] * gss0_12[k]
                  - f_2 * pc_y[k] * gss1_12[k];

        t_38[k] = f_3 * dps0_17[k]
                  - f_4 * dps1_17[k]
                  + pa_y[k] * fps0_26[k]
                  - f_2 * pc_y[k] * fps1_26[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pa_y, pc_y, fss0_9, fss1_9, fps0_27, fps0_28, \
                         fps0_29, fps1_27, fps1_28, fps1_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = pa_y[k] * fps0_27[k]
                  - f_2 * pc_y[k] * fps1_27[k];

        t_40[k] = f_3 * fss0_9[k]
                  - f_4 * fss1_9[k]
                  + pa_y[k] * fps0_28[k]
                  - f_2 * pc_y[k] * fps1_28[k];

        t_41[k] = pa_y[k] * fps0_29[k]
                  - f_2 * pc_y[k] * fps1_29[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, pb_x, pb_y, pb_z, pc_x, pc_y, pc_z, fss0_9, fss1_9, \
                         gss0_14, gss1_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_x[k] * gss0_14[k]
                  - f_2 * pc_x[k] * gss1_14[k];

        t_43[k] = pb_y[k] * gss0_14[k]
                  - f_2 * pc_y[k] * gss1_14[k];

        t_44[k] = f_0 * fss0_9[k]
                  - f_1 * fss1_9[k]
                  + pb_z[k] * gss0_14[k]
                  - f_2 * pc_z[k] * gss1_14[k];
    }
}

}  // namespace simdt3ceri
