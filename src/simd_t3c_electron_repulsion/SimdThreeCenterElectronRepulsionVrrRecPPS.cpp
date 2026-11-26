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


#include "SimdThreeCenterElectronRepulsionVrrRecPPS.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_pps_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t sss0, const size_t sss1,
                                                   const size_t sps0, const size_t sps1,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / p;
    const auto f_1 = 0.5 * gamma / (p * q);
    const auto f_2 = gamma / q;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sss0_0 = buffer.data(sss0 + 0);

    const auto *sss1_0 = buffer.data(sss1 + 0);

    const auto *sps0_0 = buffer.data(sps0 + 0);
    const auto *sps0_1 = buffer.data(sps0 + 1);
    const auto *sps0_2 = buffer.data(sps0 + 2);

    const auto *sps1_0 = buffer.data(sps1 + 0);
    const auto *sps1_1 = buffer.data(sps1 + 1);
    const auto *sps1_2 = buffer.data(sps1 + 2);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, pc_x, pc_y, sss0_0, sss1_0, sps0_0, \
                         sps0_1, sps0_2, sps1_0, sps1_1, sps1_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sss0_0[k]
                 - f_1 * sss1_0[k]
                 + pa_x[k] * sps0_0[k]
                 - f_2 * pc_x[k] * sps1_0[k];

        t_1[k] = pa_x[k] * sps0_1[k]
                 - f_2 * pc_x[k] * sps1_1[k];

        t_2[k] = pa_x[k] * sps0_2[k]
                 - f_2 * pc_x[k] * sps1_2[k];

        t_3[k] = pa_y[k] * sps0_0[k]
                 - f_2 * pc_y[k] * sps1_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_y, pa_z, pc_y, pc_z, sss0_0, sss1_0, sps0_0, \
                         sps0_1, sps0_2, sps1_0, sps1_1, sps1_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_0 * sss0_0[k]
                 - f_1 * sss1_0[k]
                 + pa_y[k] * sps0_1[k]
                 - f_2 * pc_y[k] * sps1_1[k];

        t_5[k] = pa_y[k] * sps0_2[k]
                 - f_2 * pc_y[k] * sps1_2[k];

        t_6[k] = pa_z[k] * sps0_0[k]
                 - f_2 * pc_z[k] * sps1_0[k];

        t_7[k] = pa_z[k] * sps0_1[k]
                 - f_2 * pc_z[k] * sps1_1[k];
    }

#pragma omp simd aligned(t_8, pa_z, pc_z, sss0_0, sss1_0, sps0_2, \
                         sps1_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * sss0_0[k]
                 - f_1 * sss1_0[k]
                 + pa_z[k] * sps0_2[k]
                 - f_2 * pc_z[k] * sps1_2[k];
    }
}

}  // namespace simdt3ceri
