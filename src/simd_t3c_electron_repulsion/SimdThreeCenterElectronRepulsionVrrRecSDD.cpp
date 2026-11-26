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


#include "SimdThreeCenterElectronRepulsionVrrRecSDD.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_sdd_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t spd0, const size_t spp,
                                                   const size_t spd1, const size_t sds0,
                                                   const size_t sds1, const size_t sdp,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / q;
    const auto f_1 = 0.5 / gamma;
    const auto f_2 = 0.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *spd0_0 = buffer.data(spd0 + 0);
    const auto *spd0_9 = buffer.data(spd0 + 9);
    const auto *spd0_11 = buffer.data(spd0 + 11);
    const auto *spd0_12 = buffer.data(spd0 + 12);
    const auto *spd0_15 = buffer.data(spd0 + 15);
    const auto *spd0_17 = buffer.data(spd0 + 17);

    const auto *spp_0 = buffer.data(spp + 0);
    const auto *spp_1 = buffer.data(spp + 1);
    const auto *spp_2 = buffer.data(spp + 2);
    const auto *spp_4 = buffer.data(spp + 4);
    const auto *spp_5 = buffer.data(spp + 5);
    const auto *spp_7 = buffer.data(spp + 7);
    const auto *spp_8 = buffer.data(spp + 8);

    const auto *spd1_0 = buffer.data(spd1 + 0);
    const auto *spd1_9 = buffer.data(spd1 + 9);
    const auto *spd1_11 = buffer.data(spd1 + 11);
    const auto *spd1_12 = buffer.data(spd1 + 12);
    const auto *spd1_15 = buffer.data(spd1 + 15);
    const auto *spd1_17 = buffer.data(spd1 + 17);

    const auto *sds0_0 = buffer.data(sds0 + 0);
    const auto *sds0_3 = buffer.data(sds0 + 3);
    const auto *sds0_5 = buffer.data(sds0 + 5);

    const auto *sds1_0 = buffer.data(sds1 + 0);
    const auto *sds1_3 = buffer.data(sds1 + 3);
    const auto *sds1_5 = buffer.data(sds1 + 5);

    const auto *sdp_0 = buffer.data(sdp + 0);
    const auto *sdp_1 = buffer.data(sdp + 1);
    const auto *sdp_2 = buffer.data(sdp + 2);
    const auto *sdp_4 = buffer.data(sdp + 4);
    const auto *sdp_5 = buffer.data(sdp + 5);
    const auto *sdp_7 = buffer.data(sdp + 7);
    const auto *sdp_8 = buffer.data(sdp + 8);
    const auto *sdp_9 = buffer.data(sdp + 9);
    const auto *sdp_10 = buffer.data(sdp + 10);
    const auto *sdp_11 = buffer.data(sdp + 11);
    const auto *sdp_13 = buffer.data(sdp + 13);
    const auto *sdp_14 = buffer.data(sdp + 14);
    const auto *sdp_15 = buffer.data(sdp + 15);
    const auto *sdp_16 = buffer.data(sdp + 16);
    const auto *sdp_17 = buffer.data(sdp + 17);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pc_x, pc_y, spp_0, spp_1, spp_2, sds0_0, \
                         sds1_0, sdp_0, sdp_1, sdp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * spp_0[k]
                 + f_1 * sds0_0[k]
                 - f_2 * sds1_0[k]
                 + f_3 * pc_x[k] * sdp_0[k];

        t_1[k] = f_0 * spp_1[k]
                 + f_3 * pc_x[k] * sdp_1[k];

        t_2[k] = f_0 * spp_2[k]
                 + f_3 * pc_x[k] * sdp_2[k];

        t_3[k] = f_1 * sds0_0[k]
                 - f_2 * sds1_0[k]
                 + f_3 * pc_y[k] * sdp_1[k];

        t_4[k] = f_3 * pc_y[k] * sdp_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pb_y, pc_x, pc_y, pc_z, spd0_0, spp_4, spd1_0, sds0_0, \
                         sds1_0, sdp_2, sdp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * sds0_0[k]
                 - f_2 * sds1_0[k]
                 + f_3 * pc_z[k] * sdp_2[k];

        t_6[k] = pb_y[k] * spd0_0[k]
                 - f_4 * pc_y[k] * spd1_0[k];

        t_7[k] = f_5 * spp_4[k]
                 + f_3 * pc_x[k] * sdp_4[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pb_x, pc_x, pc_y, spd0_9, spd0_11, spp_2, \
                         spp_5, spd1_9, spd1_11, sdp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_5 * spp_5[k]
                 + f_3 * pc_x[k] * sdp_5[k];

        t_9[k] = pb_x[k] * spd0_9[k]
                 - f_4 * pc_x[k] * spd1_9[k];

        t_10[k] = f_5 * spp_2[k]
                  + f_3 * pc_y[k] * sdp_5[k];

        t_11[k] = pb_x[k] * spd0_11[k]
                  - f_4 * pc_x[k] * spd1_11[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pb_x, pb_z, pc_x, pc_z, spd0_0, spd0_15, \
                         spp_7, spp_8, spd1_0, spd1_15, sdp_7, sdp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pb_z[k] * spd0_0[k]
                  - f_4 * pc_z[k] * spd1_0[k];

        t_13[k] = f_5 * spp_7[k]
                  + f_3 * pc_x[k] * sdp_7[k];

        t_14[k] = f_5 * spp_8[k]
                  + f_3 * pc_x[k] * sdp_8[k];

        t_15[k] = pb_x[k] * spd0_15[k]
                  - f_4 * pc_x[k] * spd1_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pb_x, pc_x, pc_y, spd0_17, spd1_17, \
                         sds0_3, sds1_3, sdp_8, sdp_9, sdp_10, sdp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * pc_y[k] * sdp_8[k];

        t_17[k] = pb_x[k] * spd0_17[k]
                  - f_4 * pc_x[k] * spd1_17[k];

        t_18[k] = f_1 * sds0_3[k]
                  - f_2 * sds1_3[k]
                  + f_3 * pc_x[k] * sdp_9[k];

        t_19[k] = f_3 * pc_x[k] * sdp_10[k];

        t_20[k] = f_3 * pc_x[k] * sdp_11[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pb_y, pc_y, pc_z, spd0_12, spp_4, spp_5, \
                         spd1_12, sds0_3, sds1_3, sdp_10, sdp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_0 * spp_4[k]
                  + f_1 * sds0_3[k]
                  - f_2 * sds1_3[k]
                  + f_3 * pc_y[k] * sdp_10[k];

        t_22[k] = f_0 * spp_5[k]
                  + f_3 * pc_y[k] * sdp_11[k];

        t_23[k] = f_1 * sds0_3[k]
                  - f_2 * sds1_3[k]
                  + f_3 * pc_z[k] * sdp_11[k];

        t_24[k] = pb_y[k] * spd0_12[k]
                  - f_4 * pc_y[k] * spd1_12[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pb_z, pc_x, pc_y, pc_z, spd0_9, spp_8, \
                         spd1_9, sdp_13, sdp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_3 * pc_x[k] * sdp_13[k];

        t_26[k] = f_3 * pc_x[k] * sdp_14[k];

        t_27[k] = pb_z[k] * spd0_9[k]
                  - f_4 * pc_z[k] * spd1_9[k];

        t_28[k] = f_5 * spp_8[k]
                  + f_3 * pc_y[k] * sdp_14[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, t_34, pb_y, pc_x, pc_y, spd0_17, \
                         spd1_17, sds0_5, sds1_5, sdp_15, sdp_16, \
                         sdp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pb_y[k] * spd0_17[k]
                  - f_4 * pc_y[k] * spd1_17[k];

        t_30[k] = f_1 * sds0_5[k]
                  - f_2 * sds1_5[k]
                  + f_3 * pc_x[k] * sdp_15[k];

        t_31[k] = f_3 * pc_x[k] * sdp_16[k];

        t_32[k] = f_3 * pc_x[k] * sdp_17[k];

        t_33[k] = f_1 * sds0_5[k]
                  - f_2 * sds1_5[k]
                  + f_3 * pc_y[k] * sdp_16[k];

        t_34[k] = f_3 * pc_y[k] * sdp_17[k];
    }

#pragma omp simd aligned(t_35, pc_z, spp_8, sds0_5, sds1_5, sdp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_0 * spp_8[k]
                  + f_1 * sds0_5[k]
                  - f_2 * sds1_5[k]
                  + f_3 * pc_z[k] * sdp_17[k];
    }
}

}  // namespace simdt3ceri
