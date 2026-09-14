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


#include "SimdNuclearPotentialVrrRecFS.hpp"

#include "SimdAlign.hpp"

namespace simdnpot {  // simdnpot namespace

auto
compute_prim_fs_nuclear_potential_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                    const size_t pc, const size_t ps0, const size_t ps1,
                                    const size_t ds0, const size_t ds1, const size_t ncols,
                                    const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ps0_0 = buffer.data(ps0 + 0);
    const auto *ps0_1 = buffer.data(ps0 + 1);
    const auto *ps0_2 = buffer.data(ps0 + 2);

    const auto *ps1_0 = buffer.data(ps1 + 0);
    const auto *ps1_1 = buffer.data(ps1 + 1);
    const auto *ps1_2 = buffer.data(ps1 + 2);

    const auto *ds0_0 = buffer.data(ds0 + 0);
    const auto *ds0_3 = buffer.data(ds0 + 3);
    const auto *ds0_4 = buffer.data(ds0 + 4);
    const auto *ds0_5 = buffer.data(ds0 + 5);

    const auto *ds1_0 = buffer.data(ds1 + 0);
    const auto *ds1_3 = buffer.data(ds1 + 3);
    const auto *ds1_4 = buffer.data(ds1 + 4);
    const auto *ds1_5 = buffer.data(ds1 + 5);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ps0_0, ps1_0, \
                         ds0_0, ds0_3, ds1_0, ds1_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ps0_0[k]
                 - f_0 * ps1_0[k]
                 + pa_x[k] * ds0_0[k]
                 - pc_x[k] * ds1_0[k];

        t_1[k] = pa_y[k] * ds0_0[k]
                 - pc_y[k] * ds1_0[k];

        t_2[k] = pa_z[k] * ds0_0[k]
                 - pc_z[k] * ds1_0[k];

        t_3[k] = pa_x[k] * ds0_3[k]
                 - pc_x[k] * ds1_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_x, pa_y, pc_x, pc_y, ps0_1, ps1_1, ds0_3, ds0_4, \
                         ds0_5, ds1_3, ds1_4, ds1_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = pa_x[k] * ds0_4[k]
                 - pc_x[k] * ds1_4[k];

        t_5[k] = pa_x[k] * ds0_5[k]
                 - pc_x[k] * ds1_5[k];

        t_6[k] = f_0 * ps0_1[k]
                 - f_0 * ps1_1[k]
                 + pa_y[k] * ds0_3[k]
                 - pc_y[k] * ds1_3[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_y, pa_z, pc_y, pc_z, ps0_2, ps1_2, ds0_3, ds0_5, \
                         ds1_3, ds1_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = pa_z[k] * ds0_3[k]
                 - pc_z[k] * ds1_3[k];

        t_8[k] = pa_y[k] * ds0_5[k]
                 - pc_y[k] * ds1_5[k];

        t_9[k] = f_0 * ps0_2[k]
                 - f_0 * ps1_2[k]
                 + pa_z[k] * ds0_5[k]
                 - pc_z[k] * ds1_5[k];
    }
}

}  // namespace simdnpot
