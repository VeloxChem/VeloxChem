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


#include "SimdKineticEnergyVrrRecPP.hpp"

#include "SimdAlign.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_prim_pp_kinetic_energy_0(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                 const size_t ss, const size_t pp_s, const size_t ps,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / p;
    const auto f_1 = 2.0 * alpha * beta / p;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *ss_0 = buffer.data(ss + 0);

    const auto *pp_s_0 = buffer.data(pp_s + 0);
    const auto *pp_s_1 = buffer.data(pp_s + 1);
    const auto *pp_s_2 = buffer.data(pp_s + 2);
    const auto *pp_s_3 = buffer.data(pp_s + 3);
    const auto *pp_s_4 = buffer.data(pp_s + 4);
    const auto *pp_s_5 = buffer.data(pp_s + 5);
    const auto *pp_s_6 = buffer.data(pp_s + 6);
    const auto *pp_s_7 = buffer.data(pp_s + 7);
    const auto *pp_s_8 = buffer.data(pp_s + 8);

    const auto *ps_0 = buffer.data(ps + 0);
    const auto *ps_1 = buffer.data(ps + 1);
    const auto *ps_2 = buffer.data(ps + 2);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, ss_0, pp_s_0, pp_s_1, \
                         pp_s_2, pp_s_3, pp_s_4, ps_0, ps_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ss_0[k]
                 + f_1 * pp_s_0[k]
                 + pb_x[k] * ps_0[k];

        t_1[k] = f_1 * pp_s_1[k]
                 + pb_y[k] * ps_0[k];

        t_2[k] = f_1 * pp_s_2[k]
                 + pb_z[k] * ps_0[k];

        t_3[k] = f_1 * pp_s_3[k]
                 + pb_x[k] * ps_1[k];

        t_4[k] = f_0 * ss_0[k]
                 + f_1 * pp_s_4[k]
                 + pb_y[k] * ps_1[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_x, pb_y, pb_z, ss_0, pp_s_5, pp_s_6, pp_s_7, \
                         pp_s_8, ps_1, ps_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * pp_s_5[k]
                 + pb_z[k] * ps_1[k];

        t_6[k] = f_1 * pp_s_6[k]
                 + pb_x[k] * ps_2[k];

        t_7[k] = f_1 * pp_s_7[k]
                 + pb_y[k] * ps_2[k];

        t_8[k] = f_0 * ss_0[k]
                 + f_1 * pp_s_8[k]
                 + pb_z[k] * ps_2[k];
    }
}

}  // namespace simdkin
