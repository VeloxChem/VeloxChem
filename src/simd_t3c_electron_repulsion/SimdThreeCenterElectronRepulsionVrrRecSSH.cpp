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


#include "SimdThreeCenterElectronRepulsionVrrRecSSH.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_ssh_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pc, const size_t ssf0,
                                                   const size_t ssf1, const size_t ssg,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / gamma;
    const auto f_1 = 2.0 * p / (gamma * q);
    const auto f_2 = p / q;
    const auto f_3 = 1.0 / gamma;
    const auto f_4 = p / (gamma * q);
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);

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

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ssf0_0 = buffer.data(ssf0 + 0);
    const auto *ssf0_3 = buffer.data(ssf0 + 3);
    const auto *ssf0_5 = buffer.data(ssf0 + 5);
    const auto *ssf0_6 = buffer.data(ssf0 + 6);
    const auto *ssf0_8 = buffer.data(ssf0 + 8);
    const auto *ssf0_9 = buffer.data(ssf0 + 9);

    const auto *ssf1_0 = buffer.data(ssf1 + 0);
    const auto *ssf1_3 = buffer.data(ssf1 + 3);
    const auto *ssf1_5 = buffer.data(ssf1 + 5);
    const auto *ssf1_6 = buffer.data(ssf1 + 6);
    const auto *ssf1_8 = buffer.data(ssf1 + 8);
    const auto *ssf1_9 = buffer.data(ssf1 + 9);

    const auto *ssg_0 = buffer.data(ssg + 0);
    const auto *ssg_2 = buffer.data(ssg + 2);
    const auto *ssg_3 = buffer.data(ssg + 3);
    const auto *ssg_5 = buffer.data(ssg + 5);
    const auto *ssg_6 = buffer.data(ssg + 6);
    const auto *ssg_9 = buffer.data(ssg + 9);
    const auto *ssg_10 = buffer.data(ssg + 10);
    const auto *ssg_11 = buffer.data(ssg + 11);
    const auto *ssg_12 = buffer.data(ssg + 12);
    const auto *ssg_13 = buffer.data(ssg + 13);
    const auto *ssg_14 = buffer.data(ssg + 14);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pc_x, pc_y, pc_z, ssf0_0, ssf0_3, ssf1_0, \
                         ssf1_3, ssg_0, ssg_2, ssg_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ssf0_0[k]
                 - f_1 * ssf1_0[k]
                 + f_2 * pc_x[k] * ssg_0[k];

        t_1[k] = f_2 * pc_y[k] * ssg_0[k];

        t_2[k] = f_2 * pc_z[k] * ssg_0[k];

        t_3[k] = f_3 * ssf0_3[k]
                 - f_4 * ssf1_3[k]
                 + f_2 * pc_x[k] * ssg_3[k];

        t_4[k] = f_2 * pc_y[k] * ssg_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pc_x, pc_y, pc_z, ssf0_5, ssf0_6, ssf1_5, ssf1_6, \
                         ssg_3, ssg_5, ssg_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_3 * ssf0_5[k]
                 - f_4 * ssf1_5[k]
                 + f_2 * pc_x[k] * ssg_5[k];

        t_6[k] = f_5 * ssf0_6[k]
                 - f_6 * ssf1_6[k]
                 + f_2 * pc_x[k] * ssg_6[k];

        t_7[k] = f_2 * pc_z[k] * ssg_3[k];

        t_8[k] = f_2 * pc_y[k] * ssg_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, t_13, t_14, pc_x, ssf0_9, ssf1_9, ssg_9, \
                         ssg_10, ssg_11, ssg_12, ssg_13, ssg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * ssf0_9[k]
                 - f_6 * ssf1_9[k]
                 + f_2 * pc_x[k] * ssg_9[k];

        t_10[k] = f_2 * pc_x[k] * ssg_10[k];

        t_11[k] = f_2 * pc_x[k] * ssg_11[k];

        t_12[k] = f_2 * pc_x[k] * ssg_12[k];

        t_13[k] = f_2 * pc_x[k] * ssg_13[k];

        t_14[k] = f_2 * pc_x[k] * ssg_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pc_y, pc_z, ssf0_6, ssf0_8, ssf0_9, ssf1_6, \
                         ssf1_8, ssf1_9, ssg_10, ssg_12, ssg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_0 * ssf0_6[k]
                  - f_1 * ssf1_6[k]
                  + f_2 * pc_y[k] * ssg_10[k];

        t_16[k] = f_2 * pc_z[k] * ssg_10[k];

        t_17[k] = f_3 * ssf0_8[k]
                  - f_4 * ssf1_8[k]
                  + f_2 * pc_y[k] * ssg_12[k];

        t_18[k] = f_5 * ssf0_9[k]
                  - f_6 * ssf1_9[k]
                  + f_2 * pc_y[k] * ssg_13[k];
    }

#pragma omp simd aligned(t_19, t_20, pc_y, pc_z, ssf0_9, ssf1_9, \
                         ssg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_2 * pc_y[k] * ssg_14[k];

        t_20[k] = f_0 * ssf0_9[k]
                  - f_1 * ssf1_9[k]
                  + f_2 * pc_z[k] * ssg_14[k];
    }
}

}  // namespace simdt3ceri
