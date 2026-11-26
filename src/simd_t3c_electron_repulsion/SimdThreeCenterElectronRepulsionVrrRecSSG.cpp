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


#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_ssg_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pc, const size_t ssd0,
                                                   const size_t ssd1, const size_t ssf,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / gamma;
    const auto f_1 = 1.5 * p / (gamma * q);
    const auto f_2 = p / q;
    const auto f_3 = 0.5 / gamma;
    const auto f_4 = 0.5 * p / (gamma * q);

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

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ssd0_0 = buffer.data(ssd0 + 0);
    const auto *ssd0_3 = buffer.data(ssd0 + 3);
    const auto *ssd0_5 = buffer.data(ssd0 + 5);

    const auto *ssd1_0 = buffer.data(ssd1 + 0);
    const auto *ssd1_3 = buffer.data(ssd1 + 3);
    const auto *ssd1_5 = buffer.data(ssd1 + 5);

    const auto *ssf_0 = buffer.data(ssf + 0);
    const auto *ssf_2 = buffer.data(ssf + 2);
    const auto *ssf_3 = buffer.data(ssf + 3);
    const auto *ssf_5 = buffer.data(ssf + 5);
    const auto *ssf_6 = buffer.data(ssf + 6);
    const auto *ssf_7 = buffer.data(ssf + 7);
    const auto *ssf_8 = buffer.data(ssf + 8);
    const auto *ssf_9 = buffer.data(ssf + 9);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pc_x, pc_y, pc_z, ssd0_0, ssd0_3, ssd1_0, \
                         ssd1_3, ssf_0, ssf_2, ssf_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ssd0_0[k]
                 - f_1 * ssd1_0[k]
                 + f_2 * pc_x[k] * ssf_0[k];

        t_1[k] = f_2 * pc_y[k] * ssf_0[k];

        t_2[k] = f_2 * pc_z[k] * ssf_0[k];

        t_3[k] = f_3 * ssd0_3[k]
                 - f_4 * ssd1_3[k]
                 + f_2 * pc_x[k] * ssf_3[k];

        t_4[k] = f_2 * pc_y[k] * ssf_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pc_x, ssd0_5, ssd1_5, ssf_5, ssf_6, ssf_7, \
                         ssf_8, ssf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_3 * ssd0_5[k]
                 - f_4 * ssd1_5[k]
                 + f_2 * pc_x[k] * ssf_5[k];

        t_6[k] = f_2 * pc_x[k] * ssf_6[k];

        t_7[k] = f_2 * pc_x[k] * ssf_7[k];

        t_8[k] = f_2 * pc_x[k] * ssf_8[k];

        t_9[k] = f_2 * pc_x[k] * ssf_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pc_y, pc_z, ssd0_3, ssd0_5, ssd1_3, \
                         ssd1_5, ssf_6, ssf_8, ssf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * ssd0_3[k]
                  - f_1 * ssd1_3[k]
                  + f_2 * pc_y[k] * ssf_6[k];

        t_11[k] = f_2 * pc_z[k] * ssf_6[k];

        t_12[k] = f_3 * ssd0_5[k]
                  - f_4 * ssd1_5[k]
                  + f_2 * pc_y[k] * ssf_8[k];

        t_13[k] = f_2 * pc_y[k] * ssf_9[k];

        t_14[k] = f_0 * ssd0_5[k]
                  - f_1 * ssd1_5[k]
                  + f_2 * pc_z[k] * ssf_9[k];
    }
}

}  // namespace simdt3ceri
