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


#include "SimdNuclearPotentialVrrRecSD.hpp"

#include "SimdAlign.hpp"

namespace simdnpot {  // simdnpot namespace

auto
compute_prim_sd_nuclear_potential_0(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                    const size_t pc, const size_t ss0, const size_t sp0,
                                    const size_t ss1, const size_t sp1, const size_t ncols,
                                    const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / p;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ss0_0 = buffer.data(ss0 + 0);

    const auto *sp0_0 = buffer.data(sp0 + 0);
    const auto *sp0_1 = buffer.data(sp0 + 1);
    const auto *sp0_2 = buffer.data(sp0 + 2);

    const auto *ss1_0 = buffer.data(ss1 + 0);

    const auto *sp1_0 = buffer.data(sp1 + 0);
    const auto *sp1_1 = buffer.data(sp1 + 1);
    const auto *sp1_2 = buffer.data(sp1 + 2);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pc_x, pc_y, ss0_0, sp0_0, sp0_1, \
                         sp0_2, ss1_0, sp1_0, sp1_1, sp1_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ss0_0[k]
                 + pb_x[k] * sp0_0[k]
                 - f_0 * ss1_0[k]
                 - pc_x[k] * sp1_0[k];

        t_1[k] = pb_x[k] * sp0_1[k]
                 - pc_x[k] * sp1_1[k];

        t_2[k] = pb_x[k] * sp0_2[k]
                 - pc_x[k] * sp1_2[k];

        t_3[k] = f_0 * ss0_0[k]
                 + pb_y[k] * sp0_1[k]
                 - f_0 * ss1_0[k]
                 - pc_y[k] * sp1_1[k];
    }

#pragma omp simd aligned(t_4, t_5, pb_y, pb_z, pc_y, pc_z, ss0_0, sp0_2, ss1_0, \
                         sp1_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = pb_y[k] * sp0_2[k]
                 - pc_y[k] * sp1_2[k];

        t_5[k] = f_0 * ss0_0[k]
                 + pb_z[k] * sp0_2[k]
                 - f_0 * ss1_0[k]
                 - pc_z[k] * sp1_2[k];
    }
}

}  // namespace simdnpot
