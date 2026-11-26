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


#include "SimdThreeCenterElectronRepulsionVrrRecPDS.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_pds_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pb,
                                                   const size_t pc, const size_t sps0,
                                                   const size_t sps1, const size_t sds0,
                                                   const size_t sds1, const size_t pps0,
                                                   const size_t pps1, const size_t ncols,
                                                   const double gamma, const double p,
                                                   const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = gamma / (p * q);
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sps0_0 = buffer.data(sps0 + 0);
    const auto *sps0_1 = buffer.data(sps0 + 1);
    const auto *sps0_2 = buffer.data(sps0 + 2);

    const auto *sps1_0 = buffer.data(sps1 + 0);
    const auto *sps1_1 = buffer.data(sps1 + 1);
    const auto *sps1_2 = buffer.data(sps1 + 2);

    const auto *sds0_0 = buffer.data(sds0 + 0);
    const auto *sds0_1 = buffer.data(sds0 + 1);
    const auto *sds0_2 = buffer.data(sds0 + 2);
    const auto *sds0_3 = buffer.data(sds0 + 3);
    const auto *sds0_4 = buffer.data(sds0 + 4);
    const auto *sds0_5 = buffer.data(sds0 + 5);

    const auto *sds1_0 = buffer.data(sds1 + 0);
    const auto *sds1_1 = buffer.data(sds1 + 1);
    const auto *sds1_2 = buffer.data(sds1 + 2);
    const auto *sds1_3 = buffer.data(sds1 + 3);
    const auto *sds1_4 = buffer.data(sds1 + 4);
    const auto *sds1_5 = buffer.data(sds1 + 5);

    const auto *pps0_4 = buffer.data(pps0 + 4);
    const auto *pps0_8 = buffer.data(pps0 + 8);

    const auto *pps1_4 = buffer.data(pps1 + 4);
    const auto *pps1_8 = buffer.data(pps1 + 8);

#pragma omp simd aligned(t_0, t_1, pa_x, pc_x, sps0_0, sps0_1, sps1_0, sps1_1, sds0_0, sds0_1, \
                         sds1_0, sds1_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sps0_0[k]
                 - f_1 * sps1_0[k]
                 + pa_x[k] * sds0_0[k]
                 - f_2 * pc_x[k] * sds1_0[k];

        t_1[k] = f_3 * sps0_1[k]
                 - f_4 * sps1_1[k]
                 + pa_x[k] * sds0_1[k]
                 - f_2 * pc_x[k] * sds1_1[k];
    }

#pragma omp simd aligned(t_2, t_3, t_4, t_5, pa_x, pc_x, sps0_2, sps1_2, sds0_2, sds0_3, \
                         sds0_4, sds0_5, sds1_2, sds1_3, sds1_4, \
                         sds1_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_3 * sps0_2[k]
                 - f_4 * sps1_2[k]
                 + pa_x[k] * sds0_2[k]
                 - f_2 * pc_x[k] * sds1_2[k];

        t_3[k] = pa_x[k] * sds0_3[k]
                 - f_2 * pc_x[k] * sds1_3[k];

        t_4[k] = pa_x[k] * sds0_4[k]
                 - f_2 * pc_x[k] * sds1_4[k];

        t_5[k] = pa_x[k] * sds0_5[k]
                 - f_2 * pc_x[k] * sds1_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_y, pb_x, pc_x, pc_y, sds0_0, sds0_2, sds1_0, \
                         sds1_2, pps0_4, pps1_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pa_y[k] * sds0_0[k]
                 - f_2 * pc_y[k] * sds1_0[k];

        t_7[k] = pb_x[k] * pps0_4[k]
                 - f_2 * pc_x[k] * pps1_4[k];

        t_8[k] = pa_y[k] * sds0_2[k]
                 - f_2 * pc_y[k] * sds1_2[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_y, pc_y, sps0_1, sps0_2, sps1_1, sps1_2, sds0_3, \
                         sds0_4, sds0_5, sds1_3, sds1_4, sds1_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_0 * sps0_1[k]
                 - f_1 * sps1_1[k]
                 + pa_y[k] * sds0_3[k]
                 - f_2 * pc_y[k] * sds1_3[k];

        t_10[k] = f_3 * sps0_2[k]
                  - f_4 * sps1_2[k]
                  + pa_y[k] * sds0_4[k]
                  - f_2 * pc_y[k] * sds1_4[k];

        t_11[k] = pa_y[k] * sds0_5[k]
                  - f_2 * pc_y[k] * sds1_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_z, pb_x, pc_x, pc_z, sds0_0, sds0_1, \
                         sds0_3, sds1_0, sds1_1, sds1_3, pps0_8, \
                         pps1_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pa_z[k] * sds0_0[k]
                  - f_2 * pc_z[k] * sds1_0[k];

        t_13[k] = pa_z[k] * sds0_1[k]
                  - f_2 * pc_z[k] * sds1_1[k];

        t_14[k] = pb_x[k] * pps0_8[k]
                  - f_2 * pc_x[k] * pps1_8[k];

        t_15[k] = pa_z[k] * sds0_3[k]
                  - f_2 * pc_z[k] * sds1_3[k];
    }

#pragma omp simd aligned(t_16, t_17, pa_z, pb_y, pc_y, pc_z, sps0_2, sps1_2, sds0_5, sds1_5, \
                         pps0_8, pps1_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pb_y[k] * pps0_8[k]
                  - f_2 * pc_y[k] * pps1_8[k];

        t_17[k] = f_0 * sps0_2[k]
                  - f_1 * sps1_2[k]
                  + pa_z[k] * sds0_5[k]
                  - f_2 * pc_z[k] * sds1_5[k];
    }
}

}  // namespace simdt3ceri
