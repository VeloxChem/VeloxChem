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


#include "SimdNuclearPotentialVrrRecSG.hpp"

#include "SimdAlign.hpp"

namespace simdnpot {  // simdnpot namespace

auto
compute_prim_sg_nuclear_potential_0(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                    const size_t pc, const size_t sd0, const size_t sf0,
                                    const size_t sd1, const size_t sf1, const size_t ncols,
                                    const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 0.5 / p;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sd0_0 = buffer.data(sd0 + 0);
    const auto *sd0_3 = buffer.data(sd0 + 3);
    const auto *sd0_5 = buffer.data(sd0 + 5);

    const auto *sf0_0 = buffer.data(sf0 + 0);
    const auto *sf0_2 = buffer.data(sf0 + 2);
    const auto *sf0_3 = buffer.data(sf0 + 3);
    const auto *sf0_5 = buffer.data(sf0 + 5);
    const auto *sf0_6 = buffer.data(sf0 + 6);
    const auto *sf0_7 = buffer.data(sf0 + 7);
    const auto *sf0_8 = buffer.data(sf0 + 8);
    const auto *sf0_9 = buffer.data(sf0 + 9);

    const auto *sd1_0 = buffer.data(sd1 + 0);
    const auto *sd1_3 = buffer.data(sd1 + 3);
    const auto *sd1_5 = buffer.data(sd1 + 5);

    const auto *sf1_0 = buffer.data(sf1 + 0);
    const auto *sf1_2 = buffer.data(sf1 + 2);
    const auto *sf1_3 = buffer.data(sf1 + 3);
    const auto *sf1_5 = buffer.data(sf1 + 5);
    const auto *sf1_6 = buffer.data(sf1 + 6);
    const auto *sf1_7 = buffer.data(sf1 + 7);
    const auto *sf1_8 = buffer.data(sf1 + 8);
    const auto *sf1_9 = buffer.data(sf1 + 9);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, pc_x, pc_y, pc_z, sd0_0, sf0_0, \
                         sd1_0, sf1_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sd0_0[k]
                 + pb_x[k] * sf0_0[k]
                 - f_0 * sd1_0[k]
                 - pc_x[k] * sf1_0[k];

        t_1[k] = pb_y[k] * sf0_0[k]
                 - pc_y[k] * sf1_0[k];

        t_2[k] = pb_z[k] * sf0_0[k]
                 - pc_z[k] * sf1_0[k];
    }

#pragma omp simd aligned(t_3, t_4, pb_x, pb_y, pc_x, pc_y, sd0_3, sf0_2, sf0_3, sd1_3, sf1_2, \
                         sf1_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_1 * sd0_3[k]
                 + pb_x[k] * sf0_3[k]
                 - f_1 * sd1_3[k]
                 - pc_x[k] * sf1_3[k];

        t_4[k] = pb_y[k] * sf0_2[k]
                 - pc_y[k] * sf1_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_x, pc_x, sd0_5, sf0_5, sf0_6, sf0_7, sf0_8, \
                         sd1_5, sf1_5, sf1_6, sf1_7, sf1_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * sd0_5[k]
                 + pb_x[k] * sf0_5[k]
                 - f_1 * sd1_5[k]
                 - pc_x[k] * sf1_5[k];

        t_6[k] = pb_x[k] * sf0_6[k]
                 - pc_x[k] * sf1_6[k];

        t_7[k] = pb_x[k] * sf0_7[k]
                 - pc_x[k] * sf1_7[k];

        t_8[k] = pb_x[k] * sf0_8[k]
                 - pc_x[k] * sf1_8[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pb_x, pb_y, pb_z, pc_x, pc_y, pc_z, sd0_3, sf0_6, \
                         sf0_9, sd1_3, sf1_6, sf1_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = pb_x[k] * sf0_9[k]
                 - pc_x[k] * sf1_9[k];

        t_10[k] = f_0 * sd0_3[k]
                  + pb_y[k] * sf0_6[k]
                  - f_0 * sd1_3[k]
                  - pc_y[k] * sf1_6[k];

        t_11[k] = pb_z[k] * sf0_6[k]
                  - pc_z[k] * sf1_6[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pb_y, pb_z, pc_y, pc_z, sd0_5, sf0_8, sf0_9, sd1_5, \
                         sf1_8, sf1_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_1 * sd0_5[k]
                  + pb_y[k] * sf0_8[k]
                  - f_1 * sd1_5[k]
                  - pc_y[k] * sf1_8[k];

        t_13[k] = pb_y[k] * sf0_9[k]
                  - pc_y[k] * sf1_9[k];

        t_14[k] = f_0 * sd0_5[k]
                  + pb_z[k] * sf0_9[k]
                  - f_0 * sd1_5[k]
                  - pc_z[k] * sf1_9[k];
    }
}

}  // namespace simdnpot
