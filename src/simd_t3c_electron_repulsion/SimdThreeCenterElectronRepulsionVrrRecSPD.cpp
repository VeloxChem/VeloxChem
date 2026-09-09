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


#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_spd_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t ssd0, const size_t ssp,
                                                   const size_t ssd1, const size_t spp,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / q;
    const auto f_1 = gamma / q;
    const auto f_2 = 0.5 / q;
    const auto f_3 = p / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ssd0_0 = buffer.data(ssd0 + 0);
    const auto *ssd0_3 = buffer.data(ssd0 + 3);
    const auto *ssd0_5 = buffer.data(ssd0 + 5);

    const auto *ssp_0 = buffer.data(ssp + 0);
    const auto *ssp_1 = buffer.data(ssp + 1);
    const auto *ssp_2 = buffer.data(ssp + 2);

    const auto *ssd1_0 = buffer.data(ssd1 + 0);
    const auto *ssd1_3 = buffer.data(ssd1 + 3);
    const auto *ssd1_5 = buffer.data(ssd1 + 5);

    const auto *spp_1 = buffer.data(spp + 1);
    const auto *spp_2 = buffer.data(spp + 2);
    const auto *spp_4 = buffer.data(spp + 4);
    const auto *spp_5 = buffer.data(spp + 5);
    const auto *spp_7 = buffer.data(spp + 7);
    const auto *spp_8 = buffer.data(spp + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pc_x, ssd0_0, ssd0_3, ssp_0, ssp_1, ssp_2, \
                         ssd1_0, ssd1_3, spp_1, spp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = pb_x[k] * ssd0_0[k]
                 + f_0 * ssp_0[k]
                 - f_1 * pc_x[k] * ssd1_0[k];

        t_1[k] = f_2 * ssp_1[k]
                 + f_3 * pc_x[k] * spp_1[k];

        t_2[k] = f_2 * ssp_2[k]
                 + f_3 * pc_x[k] * spp_2[k];

        t_3[k] = pb_x[k] * ssd0_3[k]
                 - f_1 * pc_x[k] * ssd1_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pb_x, pb_y, pc_x, pc_y, ssd0_0, ssd0_5, \
                         ssd1_0, ssd1_5, spp_2, spp_4, spp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * pc_y[k] * spp_2[k];

        t_5[k] = pb_x[k] * ssd0_5[k]
                 - f_1 * pc_x[k] * ssd1_5[k];

        t_6[k] = pb_y[k] * ssd0_0[k]
                 - f_1 * pc_y[k] * ssd1_0[k];

        t_7[k] = f_3 * pc_x[k] * spp_4[k];

        t_8[k] = f_3 * pc_x[k] * spp_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pb_y, pc_y, ssd0_3, ssd0_5, ssp_1, ssp_2, ssd1_3, \
                         ssd1_5, spp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = pb_y[k] * ssd0_3[k]
                 + f_0 * ssp_1[k]
                 - f_1 * pc_y[k] * ssd1_3[k];

        t_10[k] = f_2 * ssp_2[k]
                  + f_3 * pc_y[k] * spp_5[k];

        t_11[k] = pb_y[k] * ssd0_5[k]
                  - f_1 * pc_y[k] * ssd1_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, pb_z, pc_x, pc_y, pc_z, ssd0_0, ssd0_3, \
                         ssd1_0, ssd1_3, spp_7, spp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pb_z[k] * ssd0_0[k]
                  - f_1 * pc_z[k] * ssd1_0[k];

        t_13[k] = f_3 * pc_x[k] * spp_7[k];

        t_14[k] = f_3 * pc_x[k] * spp_8[k];

        t_15[k] = pb_z[k] * ssd0_3[k]
                  - f_1 * pc_z[k] * ssd1_3[k];

        t_16[k] = f_3 * pc_y[k] * spp_8[k];
    }

#pragma omp simd aligned(t_17, pb_z, pc_z, ssd0_5, ssp_2, ssd1_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = pb_z[k] * ssd0_5[k]
                  + f_0 * ssp_2[k]
                  - f_1 * pc_z[k] * ssd1_5[k];
    }
}

}  // namespace simdt3ceri
