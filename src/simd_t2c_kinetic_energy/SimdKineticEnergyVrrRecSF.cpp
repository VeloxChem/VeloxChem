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


#include "SimdKineticEnergyVrrRecSF.hpp"

#include "SimdAlign.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_prim_sf_kinetic_energy_0(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                 const size_t sp_s, const size_t sf_s, const size_t sp,
                                 const size_t sd, const size_t ncols, const double alpha,
                                 const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha / p;
    const auto f_1 = 2.0 * alpha * beta / p;
    const auto f_2 = 1.0 / p;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sp_s_0 = buffer.data(sp_s + 0);
    const auto *sp_s_1 = buffer.data(sp_s + 1);
    const auto *sp_s_2 = buffer.data(sp_s + 2);

    const auto *sf_s_0 = buffer.data(sf_s + 0);
    const auto *sf_s_1 = buffer.data(sf_s + 1);
    const auto *sf_s_2 = buffer.data(sf_s + 2);
    const auto *sf_s_3 = buffer.data(sf_s + 3);
    const auto *sf_s_4 = buffer.data(sf_s + 4);
    const auto *sf_s_5 = buffer.data(sf_s + 5);
    const auto *sf_s_6 = buffer.data(sf_s + 6);
    const auto *sf_s_7 = buffer.data(sf_s + 7);
    const auto *sf_s_8 = buffer.data(sf_s + 8);
    const auto *sf_s_9 = buffer.data(sf_s + 9);

    const auto *sp_0 = buffer.data(sp + 0);
    const auto *sp_1 = buffer.data(sp + 1);
    const auto *sp_2 = buffer.data(sp + 2);

    const auto *sd_0 = buffer.data(sd + 0);
    const auto *sd_3 = buffer.data(sd + 3);
    const auto *sd_4 = buffer.data(sd + 4);
    const auto *sd_5 = buffer.data(sd + 5);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, sp_s_0, sf_s_0, sf_s_1, sf_s_2, \
                         sf_s_3, sp_0, sd_0, sd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -f_0 * sp_s_0[k]
                 + f_1 * sf_s_0[k]
                 + f_2 * sp_0[k]
                 + pb_x[k] * sd_0[k];

        t_1[k] = f_1 * sf_s_1[k]
                 + pb_y[k] * sd_0[k];

        t_2[k] = f_1 * sf_s_2[k]
                 + pb_z[k] * sd_0[k];

        t_3[k] = f_1 * sf_s_3[k]
                 + pb_x[k] * sd_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pb_x, pb_y, pb_z, sp_s_1, sf_s_4, sf_s_5, sf_s_6, \
                         sf_s_7, sp_1, sd_3, sd_4, sd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_1 * sf_s_4[k]
                 + pb_x[k] * sd_4[k];

        t_5[k] = f_1 * sf_s_5[k]
                 + pb_x[k] * sd_5[k];

        t_6[k] = -f_0 * sp_s_1[k]
                 + f_1 * sf_s_6[k]
                 + f_2 * sp_1[k]
                 + pb_y[k] * sd_3[k];

        t_7[k] = f_1 * sf_s_7[k]
                 + pb_z[k] * sd_3[k];
    }

#pragma omp simd aligned(t_8, t_9, pb_y, pb_z, sp_s_2, sf_s_8, sf_s_9, sp_2, \
                         sd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_1 * sf_s_8[k]
                 + pb_y[k] * sd_5[k];

        t_9[k] = -f_0 * sp_s_2[k]
                 + f_1 * sf_s_9[k]
                 + f_2 * sp_2[k]
                 + pb_z[k] * sd_5[k];
    }
}

}  // namespace simdkin
