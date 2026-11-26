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


#include "SimdKineticEnergyVrrRecSD.hpp"

#include "SimdAlign.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_prim_sd_kinetic_energy_0(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                 const size_t ss_s, const size_t sd_s, const size_t ss,
                                 const size_t sp, const size_t ncols, const double alpha,
                                 const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = alpha / p;
    const auto f_1 = 2.0 * alpha * beta / p;
    const auto f_2 = 0.5 / p;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *ss_s_0 = buffer.data(ss_s + 0);

    const auto *sd_s_0 = buffer.data(sd_s + 0);
    const auto *sd_s_1 = buffer.data(sd_s + 1);
    const auto *sd_s_2 = buffer.data(sd_s + 2);
    const auto *sd_s_3 = buffer.data(sd_s + 3);
    const auto *sd_s_4 = buffer.data(sd_s + 4);
    const auto *sd_s_5 = buffer.data(sd_s + 5);

    const auto *ss_0 = buffer.data(ss + 0);

    const auto *sp_0 = buffer.data(sp + 0);
    const auto *sp_1 = buffer.data(sp + 1);
    const auto *sp_2 = buffer.data(sp + 2);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, ss_s_0, sd_s_0, sd_s_1, sd_s_2, \
                         sd_s_3, ss_0, sp_0, sp_1, sp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -f_0 * ss_s_0[k]
                 + f_1 * sd_s_0[k]
                 + f_2 * ss_0[k]
                 + pb_x[k] * sp_0[k];

        t_1[k] = f_1 * sd_s_1[k]
                 + pb_x[k] * sp_1[k];

        t_2[k] = f_1 * sd_s_2[k]
                 + pb_x[k] * sp_2[k];

        t_3[k] = -f_0 * ss_s_0[k]
                 + f_1 * sd_s_3[k]
                 + f_2 * ss_0[k]
                 + pb_y[k] * sp_1[k];
    }

#pragma omp simd aligned(t_4, t_5, pb_y, pb_z, ss_s_0, sd_s_4, sd_s_5, ss_0, \
                         sp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_1 * sd_s_4[k]
                 + pb_y[k] * sp_2[k];

        t_5[k] = -f_0 * ss_s_0[k]
                 + f_1 * sd_s_5[k]
                 + f_2 * ss_0[k]
                 + pb_z[k] * sp_2[k];
    }
}

}  // namespace simdkin
