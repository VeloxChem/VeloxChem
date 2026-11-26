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


#include "SimdNuclearPotentialVrrRecIS.hpp"

#include "SimdAlign.hpp"

namespace simdnpot {  // simdnpot namespace

auto
compute_prim_is_nuclear_potential_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                    const size_t pc, const size_t gs0, const size_t gs1,
                                    const size_t hs0, const size_t hs1, const size_t ncols,
                                    const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 1.5 / p;
    const auto f_2 = 1.0 / p;
    const auto f_3 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *gs0_0 = buffer.data(gs0 + 0);
    const auto *gs0_3 = buffer.data(gs0 + 3);
    const auto *gs0_5 = buffer.data(gs0 + 5);
    const auto *gs0_6 = buffer.data(gs0 + 6);
    const auto *gs0_9 = buffer.data(gs0 + 9);
    const auto *gs0_10 = buffer.data(gs0 + 10);
    const auto *gs0_12 = buffer.data(gs0 + 12);
    const auto *gs0_13 = buffer.data(gs0 + 13);
    const auto *gs0_14 = buffer.data(gs0 + 14);

    const auto *gs1_0 = buffer.data(gs1 + 0);
    const auto *gs1_3 = buffer.data(gs1 + 3);
    const auto *gs1_5 = buffer.data(gs1 + 5);
    const auto *gs1_6 = buffer.data(gs1 + 6);
    const auto *gs1_9 = buffer.data(gs1 + 9);
    const auto *gs1_10 = buffer.data(gs1 + 10);
    const auto *gs1_12 = buffer.data(gs1 + 12);
    const auto *gs1_13 = buffer.data(gs1 + 13);
    const auto *gs1_14 = buffer.data(gs1 + 14);

    const auto *hs0_0 = buffer.data(hs0 + 0);
    const auto *hs0_2 = buffer.data(hs0 + 2);
    const auto *hs0_3 = buffer.data(hs0 + 3);
    const auto *hs0_5 = buffer.data(hs0 + 5);
    const auto *hs0_6 = buffer.data(hs0 + 6);
    const auto *hs0_9 = buffer.data(hs0 + 9);
    const auto *hs0_10 = buffer.data(hs0 + 10);
    const auto *hs0_12 = buffer.data(hs0 + 12);
    const auto *hs0_14 = buffer.data(hs0 + 14);
    const auto *hs0_15 = buffer.data(hs0 + 15);
    const auto *hs0_16 = buffer.data(hs0 + 16);
    const auto *hs0_17 = buffer.data(hs0 + 17);
    const auto *hs0_18 = buffer.data(hs0 + 18);
    const auto *hs0_19 = buffer.data(hs0 + 19);
    const auto *hs0_20 = buffer.data(hs0 + 20);

    const auto *hs1_0 = buffer.data(hs1 + 0);
    const auto *hs1_2 = buffer.data(hs1 + 2);
    const auto *hs1_3 = buffer.data(hs1 + 3);
    const auto *hs1_5 = buffer.data(hs1 + 5);
    const auto *hs1_6 = buffer.data(hs1 + 6);
    const auto *hs1_9 = buffer.data(hs1 + 9);
    const auto *hs1_10 = buffer.data(hs1 + 10);
    const auto *hs1_12 = buffer.data(hs1 + 12);
    const auto *hs1_14 = buffer.data(hs1 + 14);
    const auto *hs1_15 = buffer.data(hs1 + 15);
    const auto *hs1_16 = buffer.data(hs1 + 16);
    const auto *hs1_17 = buffer.data(hs1 + 17);
    const auto *hs1_18 = buffer.data(hs1 + 18);
    const auto *hs1_19 = buffer.data(hs1 + 19);
    const auto *hs1_20 = buffer.data(hs1 + 20);

#pragma omp simd aligned(t_0, t_1, t_2, pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, gs0_0, gs1_0, \
                         hs0_0, hs1_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gs0_0[k]
                 - f_0 * gs1_0[k]
                 + pa_x[k] * hs0_0[k]
                 - pc_x[k] * hs1_0[k];

        t_1[k] = pa_y[k] * hs0_0[k]
                 - pc_y[k] * hs1_0[k];

        t_2[k] = pa_z[k] * hs0_0[k]
                 - pc_z[k] * hs1_0[k];
    }

#pragma omp simd aligned(t_3, t_4, pa_x, pa_y, pc_x, pc_y, gs0_3, gs1_3, hs0_2, hs0_3, hs1_2, \
                         hs1_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_1 * gs0_3[k]
                 - f_1 * gs1_3[k]
                 + pa_x[k] * hs0_3[k]
                 - pc_x[k] * hs1_3[k];

        t_4[k] = pa_y[k] * hs0_2[k]
                 - pc_y[k] * hs1_2[k];
    }

#pragma omp simd aligned(t_5, t_6, pa_x, pc_x, gs0_5, gs0_6, gs1_5, gs1_6, hs0_5, hs0_6, \
                         hs1_5, hs1_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * gs0_5[k]
                 - f_1 * gs1_5[k]
                 + pa_x[k] * hs0_5[k]
                 - pc_x[k] * hs1_5[k];

        t_6[k] = f_2 * gs0_6[k]
                 - f_2 * gs1_6[k]
                 + pa_x[k] * hs0_6[k]
                 - pc_x[k] * hs1_6[k];
    }

#pragma omp simd aligned(t_7, t_8, pa_y, pa_z, pc_y, pc_z, hs0_3, hs0_5, hs1_3, \
                         hs1_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = pa_z[k] * hs0_3[k]
                 - pc_z[k] * hs1_3[k];

        t_8[k] = pa_y[k] * hs0_5[k]
                 - pc_y[k] * hs1_5[k];
    }

#pragma omp simd aligned(t_9, t_10, pa_x, pc_x, gs0_9, gs0_10, gs1_9, gs1_10, hs0_9, hs0_10, \
                         hs1_9, hs1_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_2 * gs0_9[k]
                 - f_2 * gs1_9[k]
                 + pa_x[k] * hs0_9[k]
                 - pc_x[k] * hs1_9[k];

        t_10[k] = f_3 * gs0_10[k]
                  - f_3 * gs1_10[k]
                  + pa_x[k] * hs0_10[k]
                  - pc_x[k] * hs1_10[k];
    }

#pragma omp simd aligned(t_11, t_12, pa_x, pa_z, pc_x, pc_z, gs0_12, gs1_12, hs0_6, hs0_12, \
                         hs1_6, hs1_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pa_z[k] * hs0_6[k]
                  - pc_z[k] * hs1_6[k];

        t_12[k] = f_3 * gs0_12[k]
                  - f_3 * gs1_12[k]
                  + pa_x[k] * hs0_12[k]
                  - pc_x[k] * hs1_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_x, pa_y, pc_x, pc_y, gs0_14, gs1_14, hs0_9, \
                         hs0_14, hs0_15, hs1_9, hs1_14, hs1_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_y[k] * hs0_9[k]
                  - pc_y[k] * hs1_9[k];

        t_14[k] = f_3 * gs0_14[k]
                  - f_3 * gs1_14[k]
                  + pa_x[k] * hs0_14[k]
                  - pc_x[k] * hs1_14[k];

        t_15[k] = pa_x[k] * hs0_15[k]
                  - pc_x[k] * hs1_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pc_x, hs0_16, hs0_17, hs0_18, hs0_19, \
                         hs1_16, hs1_17, hs1_18, hs1_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pa_x[k] * hs0_16[k]
                  - pc_x[k] * hs1_16[k];

        t_17[k] = pa_x[k] * hs0_17[k]
                  - pc_x[k] * hs1_17[k];

        t_18[k] = pa_x[k] * hs0_18[k]
                  - pc_x[k] * hs1_18[k];

        t_19[k] = pa_x[k] * hs0_19[k]
                  - pc_x[k] * hs1_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, gs0_10, gs1_10, \
                         hs0_15, hs0_20, hs1_15, hs1_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_x[k] * hs0_20[k]
                  - pc_x[k] * hs1_20[k];

        t_21[k] = f_0 * gs0_10[k]
                  - f_0 * gs1_10[k]
                  + pa_y[k] * hs0_15[k]
                  - pc_y[k] * hs1_15[k];

        t_22[k] = pa_z[k] * hs0_15[k]
                  - pc_z[k] * hs1_15[k];
    }

#pragma omp simd aligned(t_23, t_24, pa_y, pc_y, gs0_12, gs0_13, gs1_12, gs1_13, hs0_17, \
                         hs0_18, hs1_17, hs1_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_1 * gs0_12[k]
                  - f_1 * gs1_12[k]
                  + pa_y[k] * hs0_17[k]
                  - pc_y[k] * hs1_17[k];

        t_24[k] = f_2 * gs0_13[k]
                  - f_2 * gs1_13[k]
                  + pa_y[k] * hs0_18[k]
                  - pc_y[k] * hs1_18[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pa_y, pa_z, pc_y, pc_z, gs0_14, gs1_14, hs0_19, \
                         hs0_20, hs1_19, hs1_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_3 * gs0_14[k]
                  - f_3 * gs1_14[k]
                  + pa_y[k] * hs0_19[k]
                  - pc_y[k] * hs1_19[k];

        t_26[k] = pa_y[k] * hs0_20[k]
                  - pc_y[k] * hs1_20[k];

        t_27[k] = f_0 * gs0_14[k]
                  - f_0 * gs1_14[k]
                  + pa_z[k] * hs0_20[k]
                  - pc_z[k] * hs1_20[k];
    }
}

}  // namespace simdnpot
