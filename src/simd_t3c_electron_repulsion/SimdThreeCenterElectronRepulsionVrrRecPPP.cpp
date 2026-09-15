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


#include "SimdThreeCenterElectronRepulsionVrrRecPPP.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_ppp_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pc, const size_t sps,
                                                   const size_t pss, const size_t pps,
                                                   const size_t ncols, const double p,
                                                   const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = p / q;

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

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sps_0 = buffer.data(sps + 0);
    const auto *sps_1 = buffer.data(sps + 1);
    const auto *sps_2 = buffer.data(sps + 2);

    const auto *pss_0 = buffer.data(pss + 0);
    const auto *pss_1 = buffer.data(pss + 1);
    const auto *pss_2 = buffer.data(pss + 2);

    const auto *pps_0 = buffer.data(pps + 0);
    const auto *pps_1 = buffer.data(pps + 1);
    const auto *pps_2 = buffer.data(pps + 2);
    const auto *pps_3 = buffer.data(pps + 3);
    const auto *pps_4 = buffer.data(pps + 4);
    const auto *pps_5 = buffer.data(pps + 5);
    const auto *pps_6 = buffer.data(pps + 6);
    const auto *pps_7 = buffer.data(pps + 7);
    const auto *pps_8 = buffer.data(pps + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, sps_0, sps_1, pss_0, \
                         pps_0, pps_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sps_0[k]
                 + f_0 * pss_0[k]
                 + f_1 * pc_x[k] * pps_0[k];

        t_1[k] = f_1 * pc_y[k] * pps_0[k];

        t_2[k] = f_1 * pc_z[k] * pps_0[k];

        t_3[k] = f_0 * sps_1[k]
                 + f_1 * pc_x[k] * pps_1[k];

        t_4[k] = f_0 * pss_0[k]
                 + f_1 * pc_y[k] * pps_1[k];

        t_5[k] = f_1 * pc_z[k] * pps_1[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, t_11, pc_x, pc_y, pc_z, sps_0, sps_2, \
                         pss_0, pss_1, pps_2, pps_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * sps_2[k]
                 + f_1 * pc_x[k] * pps_2[k];

        t_7[k] = f_1 * pc_y[k] * pps_2[k];

        t_8[k] = f_0 * pss_0[k]
                 + f_1 * pc_z[k] * pps_2[k];

        t_9[k] = f_0 * pss_1[k]
                 + f_1 * pc_x[k] * pps_3[k];

        t_10[k] = f_0 * sps_0[k]
                  + f_1 * pc_y[k] * pps_3[k];

        t_11[k] = f_1 * pc_z[k] * pps_3[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, t_17, pc_x, pc_y, pc_z, sps_1, sps_2, \
                         pss_1, pps_4, pps_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_1 * pc_x[k] * pps_4[k];

        t_13[k] = f_0 * sps_1[k]
                  + f_0 * pss_1[k]
                  + f_1 * pc_y[k] * pps_4[k];

        t_14[k] = f_1 * pc_z[k] * pps_4[k];

        t_15[k] = f_1 * pc_x[k] * pps_5[k];

        t_16[k] = f_0 * sps_2[k]
                  + f_1 * pc_y[k] * pps_5[k];

        t_17[k] = f_0 * pss_1[k]
                  + f_1 * pc_z[k] * pps_5[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, t_23, t_24, pc_x, pc_y, pc_z, sps_0, \
                         sps_1, pss_2, pps_6, pps_7, pps_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_0 * pss_2[k]
                  + f_1 * pc_x[k] * pps_6[k];

        t_19[k] = f_1 * pc_y[k] * pps_6[k];

        t_20[k] = f_0 * sps_0[k]
                  + f_1 * pc_z[k] * pps_6[k];

        t_21[k] = f_1 * pc_x[k] * pps_7[k];

        t_22[k] = f_0 * pss_2[k]
                  + f_1 * pc_y[k] * pps_7[k];

        t_23[k] = f_0 * sps_1[k]
                  + f_1 * pc_z[k] * pps_7[k];

        t_24[k] = f_1 * pc_x[k] * pps_8[k];
    }

#pragma omp simd aligned(t_25, t_26, pc_y, pc_z, sps_2, pss_2, pps_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_1 * pc_y[k] * pps_8[k];

        t_26[k] = f_0 * sps_2[k]
                  + f_0 * pss_2[k]
                  + f_1 * pc_z[k] * pps_8[k];
    }
}

}  // namespace simdt3ceri
