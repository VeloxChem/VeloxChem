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


#include "SimdTransferPP.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_pp(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t sp, const size_t sd, const size_t nmax) -> void
{
    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

    const auto *sp_0 = buffer.data(sp + 0);
    const auto *sp_1 = buffer.data(sp + 1);
    const auto *sp_2 = buffer.data(sp + 2);

    const auto *sd_0 = buffer.data(sd + 0);
    const auto *sd_1 = buffer.data(sd + 1);
    const auto *sd_2 = buffer.data(sd + 2);
    const auto *sd_3 = buffer.data(sd + 3);
    const auto *sd_4 = buffer.data(sd + 4);
    const auto *sd_5 = buffer.data(sd + 5);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, ab_x, ab_y, sp_0, sp_1, sp_2, sd_0, \
                         sd_1, sd_2, sd_3, sd_4 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_0[k] = -ab_x[k] * sp_0[k]
                 + sd_0[k];

        t_1[k] = -ab_x[k] * sp_1[k]
                 + sd_1[k];

        t_2[k] = -ab_x[k] * sp_2[k]
                 + sd_2[k];

        t_3[k] = -ab_y[k] * sp_0[k]
                 + sd_1[k];

        t_4[k] = -ab_y[k] * sp_1[k]
                 + sd_3[k];

        t_5[k] = -ab_y[k] * sp_2[k]
                 + sd_4[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, ab_z, sp_0, sp_1, sp_2, sd_2, sd_4, \
                         sd_5 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_6[k] = -ab_z[k] * sp_0[k]
                 + sd_2[k];

        t_7[k] = -ab_z[k] * sp_1[k]
                 + sd_4[k];

        t_8[k] = -ab_z[k] * sp_2[k]
                 + sd_5[k];
    }
}

}  // namespace simdtrf
