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


#include "SimdElectronRepulsionGeom10VrrRecPP.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_geom_10_pp_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t sp, const size_t dp,
                                             const size_t ncols, const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *sp_0 = buffer.data(sp + 0);
    const auto *sp_1 = buffer.data(sp + 1);
    const auto *sp_2 = buffer.data(sp + 2);

    const auto *dp_0 = buffer.data(dp + 0);
    const auto *dp_1 = buffer.data(dp + 1);
    const auto *dp_2 = buffer.data(dp + 2);
    const auto *dp_3 = buffer.data(dp + 3);
    const auto *dp_4 = buffer.data(dp + 4);
    const auto *dp_5 = buffer.data(dp + 5);
    const auto *dp_6 = buffer.data(dp + 6);
    const auto *dp_7 = buffer.data(dp + 7);
    const auto *dp_8 = buffer.data(dp + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, sp_0, sp_1, sp_2, dp_0, dp_1, dp_2, \
                         dp_3, dp_4, dp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -sp_0[k]
                 + f_0 * dp_0[k];

        t_1[k] = -sp_1[k]
                 + f_0 * dp_1[k];

        t_2[k] = -sp_2[k]
                 + f_0 * dp_2[k];

        t_3[k] = f_0 * dp_3[k];

        t_4[k] = f_0 * dp_4[k];

        t_5[k] = f_0 * dp_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, dp_6, dp_7, dp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * dp_6[k];

        t_7[k] = f_0 * dp_7[k];

        t_8[k] = f_0 * dp_8[k];
    }
}

auto
compute_prim_geom_10_pp_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t sp, const size_t dp,
                                             const size_t ncols, const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *sp_0 = buffer.data(sp + 0);
    const auto *sp_1 = buffer.data(sp + 1);
    const auto *sp_2 = buffer.data(sp + 2);

    const auto *dp_3 = buffer.data(dp + 3);
    const auto *dp_4 = buffer.data(dp + 4);
    const auto *dp_5 = buffer.data(dp + 5);
    const auto *dp_9 = buffer.data(dp + 9);
    const auto *dp_10 = buffer.data(dp + 10);
    const auto *dp_11 = buffer.data(dp + 11);
    const auto *dp_12 = buffer.data(dp + 12);
    const auto *dp_13 = buffer.data(dp + 13);
    const auto *dp_14 = buffer.data(dp + 14);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, sp_0, sp_1, sp_2, dp_3, dp_4, dp_5, \
                         dp_9, dp_10, dp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dp_3[k];

        t_1[k] = f_0 * dp_4[k];

        t_2[k] = f_0 * dp_5[k];

        t_3[k] = -sp_0[k]
                 + f_0 * dp_9[k];

        t_4[k] = -sp_1[k]
                 + f_0 * dp_10[k];

        t_5[k] = -sp_2[k]
                 + f_0 * dp_11[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, dp_12, dp_13, dp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * dp_12[k];

        t_7[k] = f_0 * dp_13[k];

        t_8[k] = f_0 * dp_14[k];
    }
}

auto
compute_prim_geom_10_pp_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t sp, const size_t dp,
                                             const size_t ncols, const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *sp_0 = buffer.data(sp + 0);
    const auto *sp_1 = buffer.data(sp + 1);
    const auto *sp_2 = buffer.data(sp + 2);

    const auto *dp_6 = buffer.data(dp + 6);
    const auto *dp_7 = buffer.data(dp + 7);
    const auto *dp_8 = buffer.data(dp + 8);
    const auto *dp_12 = buffer.data(dp + 12);
    const auto *dp_13 = buffer.data(dp + 13);
    const auto *dp_14 = buffer.data(dp + 14);
    const auto *dp_15 = buffer.data(dp + 15);
    const auto *dp_16 = buffer.data(dp + 16);
    const auto *dp_17 = buffer.data(dp + 17);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, sp_0, dp_6, dp_7, dp_8, dp_12, \
                         dp_13, dp_14, dp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dp_6[k];

        t_1[k] = f_0 * dp_7[k];

        t_2[k] = f_0 * dp_8[k];

        t_3[k] = f_0 * dp_12[k];

        t_4[k] = f_0 * dp_13[k];

        t_5[k] = f_0 * dp_14[k];

        t_6[k] = -sp_0[k]
                 + f_0 * dp_15[k];
    }

#pragma omp simd aligned(t_7, t_8, sp_1, sp_2, dp_16, dp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = -sp_1[k]
                 + f_0 * dp_16[k];

        t_8[k] = -sp_2[k]
                 + f_0 * dp_17[k];
    }
}

}  // namespace simdt2ceri
