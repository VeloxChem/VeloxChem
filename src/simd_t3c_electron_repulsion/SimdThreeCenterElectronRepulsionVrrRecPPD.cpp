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


#include "SimdThreeCenterElectronRepulsionVrrRecPPD.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_ppd_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pb,
                                                   const size_t pc, const size_t spd0,
                                                   const size_t spp, const size_t spd1,
                                                   const size_t psd0, const size_t psp,
                                                   const size_t psd1, const size_t pps0,
                                                   const size_t pps1, const size_t ppp,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 0.5 / gamma;
    const auto f_2 = 0.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *spd0_0 = buffer.data(spd0 + 0);
    const auto *spd0_3 = buffer.data(spd0 + 3);
    const auto *spd0_5 = buffer.data(spd0 + 5);
    const auto *spd0_6 = buffer.data(spd0 + 6);
    const auto *spd0_9 = buffer.data(spd0 + 9);
    const auto *spd0_11 = buffer.data(spd0 + 11);
    const auto *spd0_12 = buffer.data(spd0 + 12);
    const auto *spd0_15 = buffer.data(spd0 + 15);
    const auto *spd0_17 = buffer.data(spd0 + 17);

    const auto *spp_0 = buffer.data(spp + 0);
    const auto *spp_4 = buffer.data(spp + 4);
    const auto *spp_8 = buffer.data(spp + 8);

    const auto *spd1_0 = buffer.data(spd1 + 0);
    const auto *spd1_3 = buffer.data(spd1 + 3);
    const auto *spd1_5 = buffer.data(spd1 + 5);
    const auto *spd1_6 = buffer.data(spd1 + 6);
    const auto *spd1_9 = buffer.data(spd1 + 9);
    const auto *spd1_11 = buffer.data(spd1 + 11);
    const auto *spd1_12 = buffer.data(spd1 + 12);
    const auto *spd1_15 = buffer.data(spd1 + 15);
    const auto *spd1_17 = buffer.data(spd1 + 17);

    const auto *psd0_0 = buffer.data(psd0 + 0);
    const auto *psd0_9 = buffer.data(psd0 + 9);
    const auto *psd0_17 = buffer.data(psd0 + 17);

    const auto *psp_0 = buffer.data(psp + 0);
    const auto *psp_2 = buffer.data(psp + 2);
    const auto *psp_4 = buffer.data(psp + 4);
    const auto *psp_5 = buffer.data(psp + 5);
    const auto *psp_7 = buffer.data(psp + 7);
    const auto *psp_8 = buffer.data(psp + 8);

    const auto *psd1_0 = buffer.data(psd1 + 0);
    const auto *psd1_9 = buffer.data(psd1 + 9);
    const auto *psd1_17 = buffer.data(psd1 + 17);

    const auto *pps0_0 = buffer.data(pps0 + 0);
    const auto *pps0_4 = buffer.data(pps0 + 4);
    const auto *pps0_8 = buffer.data(pps0 + 8);

    const auto *pps1_0 = buffer.data(pps1 + 0);
    const auto *pps1_4 = buffer.data(pps1 + 4);
    const auto *pps1_8 = buffer.data(pps1 + 8);

    const auto *ppp_0 = buffer.data(ppp + 0);
    const auto *ppp_1 = buffer.data(ppp + 1);
    const auto *ppp_2 = buffer.data(ppp + 2);
    const auto *ppp_3 = buffer.data(ppp + 3);
    const auto *ppp_5 = buffer.data(ppp + 5);
    const auto *ppp_6 = buffer.data(ppp + 6);
    const auto *ppp_8 = buffer.data(ppp + 8);
    const auto *ppp_10 = buffer.data(ppp + 10);
    const auto *ppp_11 = buffer.data(ppp + 11);
    const auto *ppp_12 = buffer.data(ppp + 12);
    const auto *ppp_13 = buffer.data(ppp + 13);
    const auto *ppp_14 = buffer.data(ppp + 14);
    const auto *ppp_16 = buffer.data(ppp + 16);
    const auto *ppp_17 = buffer.data(ppp + 17);
    const auto *ppp_19 = buffer.data(ppp + 19);
    const auto *ppp_20 = buffer.data(ppp + 20);
    const auto *ppp_22 = buffer.data(ppp + 22);
    const auto *ppp_23 = buffer.data(ppp + 23);
    const auto *ppp_24 = buffer.data(ppp + 24);
    const auto *ppp_25 = buffer.data(ppp + 25);
    const auto *ppp_26 = buffer.data(ppp + 26);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, spp_0, psp_0, pps0_0, \
                         pps1_0, ppp_0, ppp_1, ppp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * spp_0[k]
                 + f_0 * psp_0[k]
                 + f_1 * pps0_0[k]
                 - f_2 * pps1_0[k]
                 + f_3 * pc_x[k] * ppp_0[k];

        t_1[k] = f_3 * pc_y[k] * ppp_0[k];

        t_2[k] = f_3 * pc_z[k] * ppp_0[k];

        t_3[k] = f_1 * pps0_0[k]
                 - f_2 * pps1_0[k]
                 + f_3 * pc_y[k] * ppp_1[k];

        t_4[k] = f_3 * pc_y[k] * ppp_2[k];

        t_5[k] = f_1 * pps0_0[k]
                 - f_2 * pps1_0[k]
                 + f_3 * pc_z[k] * ppp_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pa_x, pb_y, pc_x, pc_y, pc_z, spd0_9, spd1_9, \
                         psd0_0, psp_0, psd1_0, ppp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pb_y[k] * psd0_0[k]
                 - f_4 * pc_y[k] * psd1_0[k];

        t_7[k] = f_0 * psp_0[k]
                 + f_3 * pc_y[k] * ppp_3[k];

        t_8[k] = f_3 * pc_z[k] * ppp_3[k];

        t_9[k] = pa_x[k] * spd0_9[k]
                 - f_4 * pc_x[k] * spd1_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pb_z, pc_x, pc_y, pc_z, spd0_11, \
                         spd1_11, psd0_0, psp_2, psd1_0, ppp_5, ppp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * psp_2[k]
                  + f_3 * pc_y[k] * ppp_5[k];

        t_11[k] = pa_x[k] * spd0_11[k]
                  - f_4 * pc_x[k] * spd1_11[k];

        t_12[k] = pb_z[k] * psd0_0[k]
                  - f_4 * pc_z[k] * psd1_0[k];

        t_13[k] = f_3 * pc_y[k] * ppp_6[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_x, pc_x, pc_y, pc_z, spd0_15, spd0_17, \
                         spd1_15, spd1_17, psp_0, ppp_6, ppp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_0 * psp_0[k]
                  + f_3 * pc_z[k] * ppp_6[k];

        t_15[k] = pa_x[k] * spd0_15[k]
                  - f_4 * pc_x[k] * spd1_15[k];

        t_16[k] = f_3 * pc_y[k] * ppp_8[k];

        t_17[k] = pa_x[k] * spd0_17[k]
                  - f_4 * pc_x[k] * spd1_17[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pa_y, pb_x, pc_x, pc_y, spd0_0, spd1_0, \
                         psd0_9, psp_4, psp_5, psd1_9, ppp_10, ppp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = pa_y[k] * spd0_0[k]
                  - f_4 * pc_y[k] * spd1_0[k];

        t_19[k] = f_0 * psp_4[k]
                  + f_3 * pc_x[k] * ppp_10[k];

        t_20[k] = f_0 * psp_5[k]
                  + f_3 * pc_x[k] * ppp_11[k];

        t_21[k] = pb_x[k] * psd0_9[k]
                  - f_4 * pc_x[k] * psd1_9[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_y, pc_x, pc_y, pc_z, spd0_5, spd1_5, \
                         pps0_4, pps1_4, ppp_10, ppp_12, ppp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_3 * pc_z[k] * ppp_10[k];

        t_23[k] = pa_y[k] * spd0_5[k]
                  - f_4 * pc_y[k] * spd1_5[k];

        t_24[k] = f_1 * pps0_4[k]
                  - f_2 * pps1_4[k]
                  + f_3 * pc_x[k] * ppp_12[k];

        t_25[k] = f_3 * pc_x[k] * ppp_13[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pc_x, pc_y, pc_z, spp_4, psp_4, pps0_4, \
                         pps1_4, ppp_13, ppp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_3 * pc_x[k] * ppp_14[k];

        t_27[k] = f_0 * spp_4[k]
                  + f_0 * psp_4[k]
                  + f_1 * pps0_4[k]
                  - f_2 * pps1_4[k]
                  + f_3 * pc_y[k] * ppp_13[k];

        t_28[k] = f_3 * pc_z[k] * ppp_13[k];

        t_29[k] = f_1 * pps0_4[k]
                  - f_2 * pps1_4[k]
                  + f_3 * pc_z[k] * ppp_14[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_y, pb_z, pc_x, pc_y, pc_z, spd0_12, \
                         spd1_12, psd0_9, psd1_9, ppp_16, ppp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pa_y[k] * spd0_12[k]
                  - f_4 * pc_y[k] * spd1_12[k];

        t_31[k] = f_3 * pc_x[k] * ppp_16[k];

        t_32[k] = f_3 * pc_x[k] * ppp_17[k];

        t_33[k] = pb_z[k] * psd0_9[k]
                  - f_4 * pc_z[k] * psd1_9[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pa_y, pa_z, pc_y, pc_z, spd0_0, spd0_17, spd1_0, \
                         spd1_17, psp_4, ppp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_0 * psp_4[k]
                  + f_3 * pc_z[k] * ppp_16[k];

        t_35[k] = pa_y[k] * spd0_17[k]
                  - f_4 * pc_y[k] * spd1_17[k];

        t_36[k] = pa_z[k] * spd0_0[k]
                  - f_4 * pc_z[k] * spd1_0[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pa_z, pc_x, pc_y, pc_z, spd0_3, spd1_3, \
                         psp_7, psp_8, ppp_19, ppp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_0 * psp_7[k]
                  + f_3 * pc_x[k] * ppp_19[k];

        t_38[k] = f_0 * psp_8[k]
                  + f_3 * pc_x[k] * ppp_20[k];

        t_39[k] = pa_z[k] * spd0_3[k]
                  - f_4 * pc_z[k] * spd1_3[k];

        t_40[k] = f_3 * pc_y[k] * ppp_20[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pa_z, pb_x, pc_x, pc_z, spd0_6, spd1_6, \
                         psd0_17, psd1_17, ppp_22, ppp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = pb_x[k] * psd0_17[k]
                  - f_4 * pc_x[k] * psd1_17[k];

        t_42[k] = pa_z[k] * spd0_6[k]
                  - f_4 * pc_z[k] * spd1_6[k];

        t_43[k] = f_3 * pc_x[k] * ppp_22[k];

        t_44[k] = f_3 * pc_x[k] * ppp_23[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, pa_z, pb_y, pc_y, pc_z, spd0_9, spd1_9, psd0_17, \
                         psp_8, psd1_17, ppp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = pa_z[k] * spd0_9[k]
                  - f_4 * pc_z[k] * spd1_9[k];

        t_46[k] = f_0 * psp_8[k]
                  + f_3 * pc_y[k] * ppp_23[k];

        t_47[k] = pb_y[k] * psd0_17[k]
                  - f_4 * pc_y[k] * psd1_17[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, t_53, pc_x, pc_y, pc_z, spp_8, psp_8, \
                         pps0_8, pps1_8, ppp_24, ppp_25, ppp_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_1 * pps0_8[k]
                  - f_2 * pps1_8[k]
                  + f_3 * pc_x[k] * ppp_24[k];

        t_49[k] = f_3 * pc_x[k] * ppp_25[k];

        t_50[k] = f_3 * pc_x[k] * ppp_26[k];

        t_51[k] = f_1 * pps0_8[k]
                  - f_2 * pps1_8[k]
                  + f_3 * pc_y[k] * ppp_25[k];

        t_52[k] = f_3 * pc_y[k] * ppp_26[k];

        t_53[k] = f_0 * spp_8[k]
                  + f_0 * psp_8[k]
                  + f_1 * pps0_8[k]
                  - f_2 * pps1_8[k]
                  + f_3 * pc_z[k] * ppp_26[k];
    }
}

}  // namespace simdt3ceri
