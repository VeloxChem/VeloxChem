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


#include "SimdThreeCenterElectronRepulsionVrrRecSDP.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_sdp_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pc, const size_t sps,
                                                   const size_t sds, const size_t ncols,
                                                   const double p, const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / q;
    const auto f_1 = p / q;
    const auto f_2 = 0.5 / q;

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

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sps_0 = buffer.data(sps + 0);
    const auto *sps_1 = buffer.data(sps + 1);
    const auto *sps_2 = buffer.data(sps + 2);

    const auto *sds_0 = buffer.data(sds + 0);
    const auto *sds_1 = buffer.data(sds + 1);
    const auto *sds_2 = buffer.data(sds + 2);
    const auto *sds_3 = buffer.data(sds + 3);
    const auto *sds_4 = buffer.data(sds + 4);
    const auto *sds_5 = buffer.data(sds + 5);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, pc_x, pc_y, pc_z, sps_0, sps_1, \
                         sps_2, sds_0, sds_1, sds_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sps_0[k]
                 + f_1 * pc_x[k] * sds_0[k];

        t_1[k] = f_1 * pc_y[k] * sds_0[k];

        t_2[k] = f_1 * pc_z[k] * sds_0[k];

        t_3[k] = f_2 * sps_1[k]
                 + f_1 * pc_x[k] * sds_1[k];

        t_4[k] = f_2 * sps_0[k]
                 + f_1 * pc_y[k] * sds_1[k];

        t_5[k] = f_1 * pc_z[k] * sds_1[k];

        t_6[k] = f_2 * sps_2[k]
                 + f_1 * pc_x[k] * sds_2[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, t_11, t_12, t_13, pc_x, pc_y, pc_z, sps_0, \
                         sps_1, sps_2, sds_2, sds_3, sds_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_1 * pc_y[k] * sds_2[k];

        t_8[k] = f_2 * sps_0[k]
                 + f_1 * pc_z[k] * sds_2[k];

        t_9[k] = f_1 * pc_x[k] * sds_3[k];

        t_10[k] = f_0 * sps_1[k]
                  + f_1 * pc_y[k] * sds_3[k];

        t_11[k] = f_1 * pc_z[k] * sds_3[k];

        t_12[k] = f_1 * pc_x[k] * sds_4[k];

        t_13[k] = f_2 * sps_2[k]
                  + f_1 * pc_y[k] * sds_4[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pc_x, pc_y, pc_z, sps_1, sps_2, sds_4, \
                         sds_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_2 * sps_1[k]
                  + f_1 * pc_z[k] * sds_4[k];

        t_15[k] = f_1 * pc_x[k] * sds_5[k];

        t_16[k] = f_1 * pc_y[k] * sds_5[k];

        t_17[k] = f_0 * sps_2[k]
                  + f_1 * pc_z[k] * sds_5[k];
    }
}

}  // namespace simdt3ceri
