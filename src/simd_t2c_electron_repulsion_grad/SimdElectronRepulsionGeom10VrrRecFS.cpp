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


#include "SimdElectronRepulsionGeom10VrrRecFS.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_geom_10_fs_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t ds, const size_t gs,
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
    auto *t_9 = buffer.data(target + 9);

    const auto *ds_0 = buffer.data(ds + 0);
    const auto *ds_1 = buffer.data(ds + 1);
    const auto *ds_2 = buffer.data(ds + 2);
    const auto *ds_3 = buffer.data(ds + 3);
    const auto *ds_4 = buffer.data(ds + 4);
    const auto *ds_5 = buffer.data(ds + 5);

    const auto *gs_0 = buffer.data(gs + 0);
    const auto *gs_1 = buffer.data(gs + 1);
    const auto *gs_2 = buffer.data(gs + 2);
    const auto *gs_3 = buffer.data(gs + 3);
    const auto *gs_4 = buffer.data(gs + 4);
    const auto *gs_5 = buffer.data(gs + 5);
    const auto *gs_6 = buffer.data(gs + 6);
    const auto *gs_7 = buffer.data(gs + 7);
    const auto *gs_8 = buffer.data(gs + 8);
    const auto *gs_9 = buffer.data(gs + 9);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ds_0, ds_1, ds_2, ds_3, ds_4, gs_0, gs_1, \
                         gs_2, gs_3, gs_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -3.0 * ds_0[k]
                 + f_0 * gs_0[k];

        t_1[k] = -2.0 * ds_1[k]
                 + f_0 * gs_1[k];

        t_2[k] = -2.0 * ds_2[k]
                 + f_0 * gs_2[k];

        t_3[k] = -ds_3[k]
                 + f_0 * gs_3[k];

        t_4[k] = -ds_4[k]
                 + f_0 * gs_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ds_5, gs_5, gs_6, gs_7, gs_8, \
                         gs_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -ds_5[k]
                 + f_0 * gs_5[k];

        t_6[k] = f_0 * gs_6[k];

        t_7[k] = f_0 * gs_7[k];

        t_8[k] = f_0 * gs_8[k];

        t_9[k] = f_0 * gs_9[k];
    }
}

auto
compute_prim_geom_10_fs_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t ds, const size_t gs,
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
    auto *t_9 = buffer.data(target + 9);

    const auto *ds_0 = buffer.data(ds + 0);
    const auto *ds_1 = buffer.data(ds + 1);
    const auto *ds_2 = buffer.data(ds + 2);
    const auto *ds_3 = buffer.data(ds + 3);
    const auto *ds_4 = buffer.data(ds + 4);
    const auto *ds_5 = buffer.data(ds + 5);

    const auto *gs_1 = buffer.data(gs + 1);
    const auto *gs_3 = buffer.data(gs + 3);
    const auto *gs_4 = buffer.data(gs + 4);
    const auto *gs_6 = buffer.data(gs + 6);
    const auto *gs_7 = buffer.data(gs + 7);
    const auto *gs_8 = buffer.data(gs + 8);
    const auto *gs_10 = buffer.data(gs + 10);
    const auto *gs_11 = buffer.data(gs + 11);
    const auto *gs_12 = buffer.data(gs + 12);
    const auto *gs_13 = buffer.data(gs + 13);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, ds_0, ds_1, ds_2, gs_1, gs_3, gs_4, \
                         gs_6, gs_7, gs_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gs_1[k];

        t_1[k] = -ds_0[k]
                 + f_0 * gs_3[k];

        t_2[k] = f_0 * gs_4[k];

        t_3[k] = -2.0 * ds_1[k]
                 + f_0 * gs_6[k];

        t_4[k] = -ds_2[k]
                 + f_0 * gs_7[k];

        t_5[k] = f_0 * gs_8[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, ds_3, ds_4, ds_5, gs_10, gs_11, gs_12, \
                         gs_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -3.0 * ds_3[k]
                 + f_0 * gs_10[k];

        t_7[k] = -2.0 * ds_4[k]
                 + f_0 * gs_11[k];

        t_8[k] = -ds_5[k]
                 + f_0 * gs_12[k];

        t_9[k] = f_0 * gs_13[k];
    }
}

auto
compute_prim_geom_10_fs_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t ds, const size_t gs,
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
    auto *t_9 = buffer.data(target + 9);

    const auto *ds_0 = buffer.data(ds + 0);
    const auto *ds_1 = buffer.data(ds + 1);
    const auto *ds_2 = buffer.data(ds + 2);
    const auto *ds_3 = buffer.data(ds + 3);
    const auto *ds_4 = buffer.data(ds + 4);
    const auto *ds_5 = buffer.data(ds + 5);

    const auto *gs_2 = buffer.data(gs + 2);
    const auto *gs_4 = buffer.data(gs + 4);
    const auto *gs_5 = buffer.data(gs + 5);
    const auto *gs_7 = buffer.data(gs + 7);
    const auto *gs_8 = buffer.data(gs + 8);
    const auto *gs_9 = buffer.data(gs + 9);
    const auto *gs_11 = buffer.data(gs + 11);
    const auto *gs_12 = buffer.data(gs + 12);
    const auto *gs_13 = buffer.data(gs + 13);
    const auto *gs_14 = buffer.data(gs + 14);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, ds_0, ds_1, ds_2, gs_2, gs_4, gs_5, \
                         gs_7, gs_8, gs_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gs_2[k];

        t_1[k] = f_0 * gs_4[k];

        t_2[k] = -ds_0[k]
                 + f_0 * gs_5[k];

        t_3[k] = f_0 * gs_7[k];

        t_4[k] = -ds_1[k]
                 + f_0 * gs_8[k];

        t_5[k] = -2.0 * ds_2[k]
                 + f_0 * gs_9[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, ds_3, ds_4, ds_5, gs_11, gs_12, gs_13, \
                         gs_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * gs_11[k];

        t_7[k] = -ds_3[k]
                 + f_0 * gs_12[k];

        t_8[k] = -2.0 * ds_4[k]
                 + f_0 * gs_13[k];

        t_9[k] = -3.0 * ds_5[k]
                 + f_0 * gs_14[k];
    }
}

}  // namespace simdt2ceri
