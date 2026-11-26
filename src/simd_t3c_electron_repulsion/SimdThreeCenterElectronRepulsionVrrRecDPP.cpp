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


#include "SimdThreeCenterElectronRepulsionVrrRecDPP.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_dpp_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t ppp0, const size_t pps,
                                                   const size_t ppp1, const size_t dss,
                                                   const size_t dps, const size_t ncols,
                                                   const double gamma, const double p,
                                                   const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / q;
    const auto f_1 = 0.5 / q;
    const auto f_2 = p / q;
    const auto f_3 = gamma / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ppp0_0 = buffer.data(ppp0 + 0);
    const auto *ppp0_13 = buffer.data(ppp0 + 13);
    const auto *ppp0_26 = buffer.data(ppp0 + 26);

    const auto *pps_0 = buffer.data(pps + 0);
    const auto *pps_1 = buffer.data(pps + 1);
    const auto *pps_2 = buffer.data(pps + 2);
    const auto *pps_3 = buffer.data(pps + 3);
    const auto *pps_4 = buffer.data(pps + 4);
    const auto *pps_5 = buffer.data(pps + 5);
    const auto *pps_6 = buffer.data(pps + 6);
    const auto *pps_7 = buffer.data(pps + 7);
    const auto *pps_8 = buffer.data(pps + 8);

    const auto *ppp1_0 = buffer.data(ppp1 + 0);
    const auto *ppp1_13 = buffer.data(ppp1 + 13);
    const auto *ppp1_26 = buffer.data(ppp1 + 26);

    const auto *dss_0 = buffer.data(dss + 0);
    const auto *dss_1 = buffer.data(dss + 1);
    const auto *dss_2 = buffer.data(dss + 2);
    const auto *dss_3 = buffer.data(dss + 3);
    const auto *dss_4 = buffer.data(dss + 4);
    const auto *dss_5 = buffer.data(dss + 5);

    const auto *dps_0 = buffer.data(dps + 0);
    const auto *dps_1 = buffer.data(dps + 1);
    const auto *dps_2 = buffer.data(dps + 2);
    const auto *dps_3 = buffer.data(dps + 3);
    const auto *dps_4 = buffer.data(dps + 4);
    const auto *dps_5 = buffer.data(dps + 5);
    const auto *dps_6 = buffer.data(dps + 6);
    const auto *dps_7 = buffer.data(dps + 7);
    const auto *dps_8 = buffer.data(dps + 8);
    const auto *dps_9 = buffer.data(dps + 9);
    const auto *dps_10 = buffer.data(dps + 10);
    const auto *dps_11 = buffer.data(dps + 11);
    const auto *dps_12 = buffer.data(dps + 12);
    const auto *dps_13 = buffer.data(dps + 13);
    const auto *dps_14 = buffer.data(dps + 14);
    const auto *dps_15 = buffer.data(dps + 15);
    const auto *dps_16 = buffer.data(dps + 16);
    const auto *dps_17 = buffer.data(dps + 17);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, pps_0, pps_1, dss_0, \
                         dps_0, dps_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pps_0[k]
                 + f_1 * dss_0[k]
                 + f_2 * pc_x[k] * dps_0[k];

        t_1[k] = f_2 * pc_y[k] * dps_0[k];

        t_2[k] = f_2 * pc_z[k] * dps_0[k];

        t_3[k] = f_0 * pps_1[k]
                 + f_2 * pc_x[k] * dps_1[k];

        t_4[k] = f_1 * dss_0[k]
                 + f_2 * pc_y[k] * dps_1[k];

        t_5[k] = f_2 * pc_z[k] * dps_1[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pa_y, pc_x, pc_y, pc_z, ppp0_0, pps_0, \
                         pps_2, ppp1_0, dss_0, dps_2, dps_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * pps_2[k]
                 + f_2 * pc_x[k] * dps_2[k];

        t_7[k] = f_2 * pc_y[k] * dps_2[k];

        t_8[k] = f_1 * dss_0[k]
                 + f_2 * pc_z[k] * dps_2[k];

        t_9[k] = pa_y[k] * ppp0_0[k]
                 - f_3 * pc_y[k] * ppp1_0[k];

        t_10[k] = f_1 * pps_0[k]
                  + f_2 * pc_y[k] * dps_3[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pa_x, pc_x, pc_z, ppp0_13, pps_4, \
                         pps_5, ppp1_13, dps_3, dps_4, dps_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_2 * pc_z[k] * dps_3[k];

        t_12[k] = f_1 * pps_4[k]
                  + f_2 * pc_x[k] * dps_4[k];

        t_13[k] = pa_x[k] * ppp0_13[k]
                  - f_3 * pc_x[k] * ppp1_13[k];

        t_14[k] = f_2 * pc_z[k] * dps_4[k];

        t_15[k] = f_1 * pps_5[k]
                  + f_2 * pc_x[k] * dps_5[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pa_z, pc_y, pc_z, ppp0_0, pps_0, pps_2, \
                         ppp1_0, dss_1, dps_5, dps_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_1 * pps_2[k]
                  + f_2 * pc_y[k] * dps_5[k];

        t_17[k] = f_1 * dss_1[k]
                  + f_2 * pc_z[k] * dps_5[k];

        t_18[k] = pa_z[k] * ppp0_0[k]
                  - f_3 * pc_z[k] * ppp1_0[k];

        t_19[k] = f_2 * pc_y[k] * dps_6[k];

        t_20[k] = f_1 * pps_0[k]
                  + f_2 * pc_z[k] * dps_6[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pc_x, pc_y, pc_z, pps_1, pps_7, pps_8, \
                         dss_2, dps_7, dps_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * pps_7[k]
                  + f_2 * pc_x[k] * dps_7[k];

        t_22[k] = f_1 * dss_2[k]
                  + f_2 * pc_y[k] * dps_7[k];

        t_23[k] = f_1 * pps_1[k]
                  + f_2 * pc_z[k] * dps_7[k];

        t_24[k] = f_1 * pps_8[k]
                  + f_2 * pc_x[k] * dps_8[k];

        t_25[k] = f_2 * pc_y[k] * dps_8[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, pa_x, pc_x, pc_y, pc_z, ppp0_26, pps_3, \
                         ppp1_26, dss_3, dps_9, dps_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pa_x[k] * ppp0_26[k]
                  - f_3 * pc_x[k] * ppp1_26[k];

        t_27[k] = f_1 * dss_3[k]
                  + f_2 * pc_x[k] * dps_9[k];

        t_28[k] = f_0 * pps_3[k]
                  + f_2 * pc_y[k] * dps_9[k];

        t_29[k] = f_2 * pc_z[k] * dps_9[k];

        t_30[k] = f_2 * pc_x[k] * dps_10[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, t_36, pc_x, pc_y, pc_z, pps_4, pps_5, \
                         dss_3, dss_4, dps_10, dps_11, dps_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_0 * pps_4[k]
                  + f_1 * dss_3[k]
                  + f_2 * pc_y[k] * dps_10[k];

        t_32[k] = f_2 * pc_z[k] * dps_10[k];

        t_33[k] = f_2 * pc_x[k] * dps_11[k];

        t_34[k] = f_0 * pps_5[k]
                  + f_2 * pc_y[k] * dps_11[k];

        t_35[k] = f_1 * dss_3[k]
                  + f_2 * pc_z[k] * dps_11[k];

        t_36[k] = f_1 * dss_4[k]
                  + f_2 * pc_x[k] * dps_12[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, t_41, pa_z, pc_x, pc_y, pc_z, ppp0_13, pps_3, \
                         pps_4, pps_6, ppp1_13, dps_12, dps_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_1 * pps_6[k]
                  + f_2 * pc_y[k] * dps_12[k];

        t_38[k] = f_1 * pps_3[k]
                  + f_2 * pc_z[k] * dps_12[k];

        t_39[k] = f_2 * pc_x[k] * dps_13[k];

        t_40[k] = pa_z[k] * ppp0_13[k]
                  - f_3 * pc_z[k] * ppp1_13[k];

        t_41[k] = f_1 * pps_4[k]
                  + f_2 * pc_z[k] * dps_13[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, pa_y, pc_x, pc_y, ppp0_26, pps_8, \
                         ppp1_26, dss_5, dps_14, dps_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_2 * pc_x[k] * dps_14[k];

        t_43[k] = f_1 * pps_8[k]
                  + f_2 * pc_y[k] * dps_14[k];

        t_44[k] = pa_y[k] * ppp0_26[k]
                  - f_3 * pc_y[k] * ppp1_26[k];

        t_45[k] = f_1 * dss_5[k]
                  + f_2 * pc_x[k] * dps_15[k];

        t_46[k] = f_2 * pc_y[k] * dps_15[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, t_51, t_52, pc_x, pc_y, pc_z, pps_6, pps_7, \
                         dss_5, dps_15, dps_16, dps_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_0 * pps_6[k]
                  + f_2 * pc_z[k] * dps_15[k];

        t_48[k] = f_2 * pc_x[k] * dps_16[k];

        t_49[k] = f_1 * dss_5[k]
                  + f_2 * pc_y[k] * dps_16[k];

        t_50[k] = f_0 * pps_7[k]
                  + f_2 * pc_z[k] * dps_16[k];

        t_51[k] = f_2 * pc_x[k] * dps_17[k];

        t_52[k] = f_2 * pc_y[k] * dps_17[k];
    }

#pragma omp simd aligned(t_53, pc_z, pps_8, dss_5, dps_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_0 * pps_8[k]
                  + f_1 * dss_5[k]
                  + f_2 * pc_z[k] * dps_17[k];
    }
}

}  // namespace simdt3ceri
