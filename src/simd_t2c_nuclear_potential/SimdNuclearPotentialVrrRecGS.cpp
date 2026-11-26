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


#include "SimdNuclearPotentialVrrRecGS.hpp"

#include "SimdAlign.hpp"

namespace simdnpot {  // simdnpot namespace

auto
compute_prim_gs_nuclear_potential_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                    const size_t pc, const size_t ds0, const size_t ds1,
                                    const size_t fs0, const size_t fs1, const size_t ncols,
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ds0_0 = buffer.data(ds0 + 0);
    const auto *ds0_3 = buffer.data(ds0 + 3);
    const auto *ds0_5 = buffer.data(ds0 + 5);

    const auto *ds1_0 = buffer.data(ds1 + 0);
    const auto *ds1_3 = buffer.data(ds1 + 3);
    const auto *ds1_5 = buffer.data(ds1 + 5);

    const auto *fs0_0 = buffer.data(fs0 + 0);
    const auto *fs0_2 = buffer.data(fs0 + 2);
    const auto *fs0_3 = buffer.data(fs0 + 3);
    const auto *fs0_5 = buffer.data(fs0 + 5);
    const auto *fs0_6 = buffer.data(fs0 + 6);
    const auto *fs0_7 = buffer.data(fs0 + 7);
    const auto *fs0_8 = buffer.data(fs0 + 8);
    const auto *fs0_9 = buffer.data(fs0 + 9);

    const auto *fs1_0 = buffer.data(fs1 + 0);
    const auto *fs1_2 = buffer.data(fs1 + 2);
    const auto *fs1_3 = buffer.data(fs1 + 3);
    const auto *fs1_5 = buffer.data(fs1 + 5);
    const auto *fs1_6 = buffer.data(fs1 + 6);
    const auto *fs1_7 = buffer.data(fs1 + 7);
    const auto *fs1_8 = buffer.data(fs1 + 8);
    const auto *fs1_9 = buffer.data(fs1 + 9);

#pragma omp simd aligned(t_0, t_1, t_2, pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ds0_0, ds1_0, \
                         fs0_0, fs1_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ds0_0[k]
                 - f_0 * ds1_0[k]
                 + pa_x[k] * fs0_0[k]
                 - pc_x[k] * fs1_0[k];

        t_1[k] = pa_y[k] * fs0_0[k]
                 - pc_y[k] * fs1_0[k];

        t_2[k] = pa_z[k] * fs0_0[k]
                 - pc_z[k] * fs1_0[k];
    }

#pragma omp simd aligned(t_3, t_4, pa_x, pa_y, pc_x, pc_y, ds0_3, ds1_3, fs0_2, fs0_3, fs1_2, \
                         fs1_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_1 * ds0_3[k]
                 - f_1 * ds1_3[k]
                 + pa_x[k] * fs0_3[k]
                 - pc_x[k] * fs1_3[k];

        t_4[k] = pa_y[k] * fs0_2[k]
                 - pc_y[k] * fs1_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_x, pc_x, ds0_5, ds1_5, fs0_5, fs0_6, fs0_7, \
                         fs0_8, fs1_5, fs1_6, fs1_7, fs1_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * ds0_5[k]
                 - f_1 * ds1_5[k]
                 + pa_x[k] * fs0_5[k]
                 - pc_x[k] * fs1_5[k];

        t_6[k] = pa_x[k] * fs0_6[k]
                 - pc_x[k] * fs1_6[k];

        t_7[k] = pa_x[k] * fs0_7[k]
                 - pc_x[k] * fs1_7[k];

        t_8[k] = pa_x[k] * fs0_8[k]
                 - pc_x[k] * fs1_8[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ds0_3, ds1_3, \
                         fs0_6, fs0_9, fs1_6, fs1_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = pa_x[k] * fs0_9[k]
                 - pc_x[k] * fs1_9[k];

        t_10[k] = f_0 * ds0_3[k]
                  - f_0 * ds1_3[k]
                  + pa_y[k] * fs0_6[k]
                  - pc_y[k] * fs1_6[k];

        t_11[k] = pa_z[k] * fs0_6[k]
                  - pc_z[k] * fs1_6[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_y, pa_z, pc_y, pc_z, ds0_5, ds1_5, fs0_8, fs0_9, \
                         fs1_8, fs1_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_1 * ds0_5[k]
                  - f_1 * ds1_5[k]
                  + pa_y[k] * fs0_8[k]
                  - pc_y[k] * fs1_8[k];

        t_13[k] = pa_y[k] * fs0_9[k]
                  - pc_y[k] * fs1_9[k];

        t_14[k] = f_0 * ds0_5[k]
                  - f_0 * ds1_5[k]
                  + pa_z[k] * fs0_9[k]
                  - pc_z[k] * fs1_9[k];
    }
}

}  // namespace simdnpot
