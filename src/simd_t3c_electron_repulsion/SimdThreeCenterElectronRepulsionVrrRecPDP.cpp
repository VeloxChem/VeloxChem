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


#include "SimdThreeCenterElectronRepulsionVrrRecPDP.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_pdp_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t sds, const size_t ppp0,
                                                   const size_t pps, const size_t ppp1,
                                                   const size_t pds, const size_t ncols,
                                                   const double gamma, const double p,
                                                   const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 1.0 / q;
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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sds_0 = buffer.data(sds + 0);
    const auto *sds_1 = buffer.data(sds + 1);
    const auto *sds_2 = buffer.data(sds + 2);
    const auto *sds_3 = buffer.data(sds + 3);
    const auto *sds_4 = buffer.data(sds + 4);
    const auto *sds_5 = buffer.data(sds + 5);

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

    const auto *pds_0 = buffer.data(pds + 0);
    const auto *pds_1 = buffer.data(pds + 1);
    const auto *pds_2 = buffer.data(pds + 2);
    const auto *pds_3 = buffer.data(pds + 3);
    const auto *pds_4 = buffer.data(pds + 4);
    const auto *pds_5 = buffer.data(pds + 5);
    const auto *pds_6 = buffer.data(pds + 6);
    const auto *pds_7 = buffer.data(pds + 7);
    const auto *pds_8 = buffer.data(pds + 8);
    const auto *pds_9 = buffer.data(pds + 9);
    const auto *pds_10 = buffer.data(pds + 10);
    const auto *pds_11 = buffer.data(pds + 11);
    const auto *pds_12 = buffer.data(pds + 12);
    const auto *pds_13 = buffer.data(pds + 13);
    const auto *pds_14 = buffer.data(pds + 14);
    const auto *pds_15 = buffer.data(pds + 15);
    const auto *pds_16 = buffer.data(pds + 16);
    const auto *pds_17 = buffer.data(pds + 17);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_y, pc_x, pc_y, pc_z, sds_0, ppp0_0, \
                         pps_0, ppp1_0, pds_0, pds_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sds_0[k]
                 + f_1 * pps_0[k]
                 + f_2 * pc_x[k] * pds_0[k];

        t_1[k] = f_2 * pc_y[k] * pds_0[k];

        t_2[k] = f_2 * pc_z[k] * pds_0[k];

        t_3[k] = pb_y[k] * ppp0_0[k]
                 - f_3 * pc_y[k] * ppp1_0[k];

        t_4[k] = f_0 * pps_0[k]
                 + f_2 * pc_y[k] * pds_1[k];

        t_5[k] = f_2 * pc_z[k] * pds_1[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pb_z, pc_x, pc_y, pc_z, sds_3, ppp0_0, \
                         pps_0, pps_1, ppp1_0, pds_2, pds_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pb_z[k] * ppp0_0[k]
                 - f_3 * pc_z[k] * ppp1_0[k];

        t_7[k] = f_2 * pc_y[k] * pds_2[k];

        t_8[k] = f_0 * pps_0[k]
                 + f_2 * pc_z[k] * pds_2[k];

        t_9[k] = f_0 * sds_3[k]
                 + f_2 * pc_x[k] * pds_3[k];

        t_10[k] = f_1 * pps_1[k]
                  + f_2 * pc_y[k] * pds_3[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, pc_x, pc_y, pc_z, sds_4, sds_5, \
                         pps_1, pps_2, pds_3, pds_4, pds_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_2 * pc_z[k] * pds_3[k];

        t_12[k] = f_0 * sds_4[k]
                  + f_2 * pc_x[k] * pds_4[k];

        t_13[k] = f_0 * pps_2[k]
                  + f_2 * pc_y[k] * pds_4[k];

        t_14[k] = f_0 * pps_1[k]
                  + f_2 * pc_z[k] * pds_4[k];

        t_15[k] = f_0 * sds_5[k]
                  + f_2 * pc_x[k] * pds_5[k];

        t_16[k] = f_2 * pc_y[k] * pds_5[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, t_21, pc_x, pc_y, pc_z, sds_0, pps_2, pps_3, \
                         pps_4, pds_5, pds_6, pds_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_1 * pps_2[k]
                  + f_2 * pc_z[k] * pds_5[k];

        t_18[k] = f_1 * pps_3[k]
                  + f_2 * pc_x[k] * pds_6[k];

        t_19[k] = f_0 * sds_0[k]
                  + f_2 * pc_y[k] * pds_6[k];

        t_20[k] = f_2 * pc_z[k] * pds_6[k];

        t_21[k] = f_0 * pps_4[k]
                  + f_2 * pc_x[k] * pds_7[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, pb_x, pc_x, pc_y, pc_z, sds_2, ppp0_13, \
                         pps_3, pps_5, ppp1_13, pds_7, pds_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = pb_x[k] * ppp0_13[k]
                  - f_3 * pc_x[k] * ppp1_13[k];

        t_23[k] = f_2 * pc_z[k] * pds_7[k];

        t_24[k] = f_0 * pps_5[k]
                  + f_2 * pc_x[k] * pds_8[k];

        t_25[k] = f_0 * sds_2[k]
                  + f_2 * pc_y[k] * pds_8[k];

        t_26[k] = f_0 * pps_3[k]
                  + f_2 * pc_z[k] * pds_8[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, t_32, pb_z, pc_x, pc_y, pc_z, sds_3, \
                         ppp0_13, pps_4, ppp1_13, pds_9, pds_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_2 * pc_x[k] * pds_9[k];

        t_28[k] = f_0 * sds_3[k]
                  + f_1 * pps_4[k]
                  + f_2 * pc_y[k] * pds_9[k];

        t_29[k] = f_2 * pc_z[k] * pds_9[k];

        t_30[k] = f_2 * pc_x[k] * pds_10[k];

        t_31[k] = pb_z[k] * ppp0_13[k]
                  - f_3 * pc_z[k] * ppp1_13[k];

        t_32[k] = f_0 * pps_4[k]
                  + f_2 * pc_z[k] * pds_10[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, t_38, pc_x, pc_y, pc_z, sds_0, sds_5, \
                         pps_5, pps_6, pds_11, pds_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_2 * pc_x[k] * pds_11[k];

        t_34[k] = f_0 * sds_5[k]
                  + f_2 * pc_y[k] * pds_11[k];

        t_35[k] = f_1 * pps_5[k]
                  + f_2 * pc_z[k] * pds_11[k];

        t_36[k] = f_1 * pps_6[k]
                  + f_2 * pc_x[k] * pds_12[k];

        t_37[k] = f_2 * pc_y[k] * pds_12[k];

        t_38[k] = f_0 * sds_0[k]
                  + f_2 * pc_z[k] * pds_12[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, t_43, pc_x, pc_y, pc_z, sds_1, pps_6, pps_7, \
                         pps_8, pds_13, pds_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_0 * pps_7[k]
                  + f_2 * pc_x[k] * pds_13[k];

        t_40[k] = f_0 * pps_6[k]
                  + f_2 * pc_y[k] * pds_13[k];

        t_41[k] = f_0 * sds_1[k]
                  + f_2 * pc_z[k] * pds_13[k];

        t_42[k] = f_0 * pps_8[k]
                  + f_2 * pc_x[k] * pds_14[k];

        t_43[k] = f_2 * pc_y[k] * pds_14[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, pb_x, pc_x, pc_y, pc_z, sds_3, ppp0_26, \
                         pps_7, ppp1_26, pds_15, pds_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = pb_x[k] * ppp0_26[k]
                  - f_3 * pc_x[k] * ppp1_26[k];

        t_45[k] = f_2 * pc_x[k] * pds_15[k];

        t_46[k] = f_1 * pps_7[k]
                  + f_2 * pc_y[k] * pds_15[k];

        t_47[k] = f_0 * sds_3[k]
                  + f_2 * pc_z[k] * pds_15[k];

        t_48[k] = f_2 * pc_x[k] * pds_16[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, pb_y, pc_x, pc_y, pc_z, sds_5, ppp0_26, \
                         pps_8, ppp1_26, pds_16, pds_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_0 * pps_8[k]
                  + f_2 * pc_y[k] * pds_16[k];

        t_50[k] = pb_y[k] * ppp0_26[k]
                  - f_3 * pc_y[k] * ppp1_26[k];

        t_51[k] = f_2 * pc_x[k] * pds_17[k];

        t_52[k] = f_2 * pc_y[k] * pds_17[k];

        t_53[k] = f_0 * sds_5[k]
                  + f_1 * pps_8[k]
                  + f_2 * pc_z[k] * pds_17[k];
    }
}

}  // namespace simdt3ceri
