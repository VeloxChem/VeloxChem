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


#include "SimdElectronRepulsionGeom10VrrRecGS.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_geom_10_gs_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t fs, const size_t hs,
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
    auto *t_10 = buffer.data(target + 10);
    auto *t_11 = buffer.data(target + 11);
    auto *t_12 = buffer.data(target + 12);
    auto *t_13 = buffer.data(target + 13);
    auto *t_14 = buffer.data(target + 14);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_1 = buffer.data(fs + 1);
    const auto *fs_2 = buffer.data(fs + 2);
    const auto *fs_3 = buffer.data(fs + 3);
    const auto *fs_4 = buffer.data(fs + 4);
    const auto *fs_5 = buffer.data(fs + 5);
    const auto *fs_6 = buffer.data(fs + 6);
    const auto *fs_7 = buffer.data(fs + 7);
    const auto *fs_8 = buffer.data(fs + 8);
    const auto *fs_9 = buffer.data(fs + 9);

    const auto *hs_0 = buffer.data(hs + 0);
    const auto *hs_1 = buffer.data(hs + 1);
    const auto *hs_2 = buffer.data(hs + 2);
    const auto *hs_3 = buffer.data(hs + 3);
    const auto *hs_4 = buffer.data(hs + 4);
    const auto *hs_5 = buffer.data(hs + 5);
    const auto *hs_6 = buffer.data(hs + 6);
    const auto *hs_7 = buffer.data(hs + 7);
    const auto *hs_8 = buffer.data(hs + 8);
    const auto *hs_9 = buffer.data(hs + 9);
    const auto *hs_10 = buffer.data(hs + 10);
    const auto *hs_11 = buffer.data(hs + 11);
    const auto *hs_12 = buffer.data(hs + 12);
    const auto *hs_13 = buffer.data(hs + 13);
    const auto *hs_14 = buffer.data(hs + 14);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, fs_0, fs_1, fs_2, fs_3, fs_4, hs_0, hs_1, \
                         hs_2, hs_3, hs_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -4.0 * fs_0[k]
                 + f_0 * hs_0[k];

        t_1[k] = -3.0 * fs_1[k]
                 + f_0 * hs_1[k];

        t_2[k] = -3.0 * fs_2[k]
                 + f_0 * hs_2[k];

        t_3[k] = -2.0 * fs_3[k]
                 + f_0 * hs_3[k];

        t_4[k] = -2.0 * fs_4[k]
                 + f_0 * hs_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, fs_5, fs_6, fs_7, fs_8, fs_9, hs_5, hs_6, \
                         hs_7, hs_8, hs_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -2.0 * fs_5[k]
                 + f_0 * hs_5[k];

        t_6[k] = -fs_6[k]
                 + f_0 * hs_6[k];

        t_7[k] = -fs_7[k]
                 + f_0 * hs_7[k];

        t_8[k] = -fs_8[k]
                 + f_0 * hs_8[k];

        t_9[k] = -fs_9[k]
                 + f_0 * hs_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, hs_10, hs_11, hs_12, hs_13, \
                         hs_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * hs_10[k];

        t_11[k] = f_0 * hs_11[k];

        t_12[k] = f_0 * hs_12[k];

        t_13[k] = f_0 * hs_13[k];

        t_14[k] = f_0 * hs_14[k];
    }
}

auto
compute_prim_geom_10_gs_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t fs, const size_t hs,
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
    auto *t_10 = buffer.data(target + 10);
    auto *t_11 = buffer.data(target + 11);
    auto *t_12 = buffer.data(target + 12);
    auto *t_13 = buffer.data(target + 13);
    auto *t_14 = buffer.data(target + 14);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_1 = buffer.data(fs + 1);
    const auto *fs_2 = buffer.data(fs + 2);
    const auto *fs_3 = buffer.data(fs + 3);
    const auto *fs_4 = buffer.data(fs + 4);
    const auto *fs_5 = buffer.data(fs + 5);
    const auto *fs_6 = buffer.data(fs + 6);
    const auto *fs_7 = buffer.data(fs + 7);
    const auto *fs_8 = buffer.data(fs + 8);
    const auto *fs_9 = buffer.data(fs + 9);

    const auto *hs_1 = buffer.data(hs + 1);
    const auto *hs_3 = buffer.data(hs + 3);
    const auto *hs_4 = buffer.data(hs + 4);
    const auto *hs_6 = buffer.data(hs + 6);
    const auto *hs_7 = buffer.data(hs + 7);
    const auto *hs_8 = buffer.data(hs + 8);
    const auto *hs_10 = buffer.data(hs + 10);
    const auto *hs_11 = buffer.data(hs + 11);
    const auto *hs_12 = buffer.data(hs + 12);
    const auto *hs_13 = buffer.data(hs + 13);
    const auto *hs_15 = buffer.data(hs + 15);
    const auto *hs_16 = buffer.data(hs + 16);
    const auto *hs_17 = buffer.data(hs + 17);
    const auto *hs_18 = buffer.data(hs + 18);
    const auto *hs_19 = buffer.data(hs + 19);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, fs_0, fs_1, fs_2, hs_1, hs_3, hs_4, \
                         hs_6, hs_7, hs_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hs_1[k];

        t_1[k] = -fs_0[k]
                 + f_0 * hs_3[k];

        t_2[k] = f_0 * hs_4[k];

        t_3[k] = -2.0 * fs_1[k]
                 + f_0 * hs_6[k];

        t_4[k] = -fs_2[k]
                 + f_0 * hs_7[k];

        t_5[k] = f_0 * hs_8[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, fs_3, fs_4, fs_5, fs_6, hs_10, hs_11, \
                         hs_12, hs_13, hs_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -3.0 * fs_3[k]
                 + f_0 * hs_10[k];

        t_7[k] = -2.0 * fs_4[k]
                 + f_0 * hs_11[k];

        t_8[k] = -fs_5[k]
                 + f_0 * hs_12[k];

        t_9[k] = f_0 * hs_13[k];

        t_10[k] = -4.0 * fs_6[k]
                  + f_0 * hs_15[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, fs_7, fs_8, fs_9, hs_16, hs_17, hs_18, \
                         hs_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = -3.0 * fs_7[k]
                  + f_0 * hs_16[k];

        t_12[k] = -2.0 * fs_8[k]
                  + f_0 * hs_17[k];

        t_13[k] = -fs_9[k]
                  + f_0 * hs_18[k];

        t_14[k] = f_0 * hs_19[k];
    }
}

auto
compute_prim_geom_10_gs_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t fs, const size_t hs,
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
    auto *t_10 = buffer.data(target + 10);
    auto *t_11 = buffer.data(target + 11);
    auto *t_12 = buffer.data(target + 12);
    auto *t_13 = buffer.data(target + 13);
    auto *t_14 = buffer.data(target + 14);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_1 = buffer.data(fs + 1);
    const auto *fs_2 = buffer.data(fs + 2);
    const auto *fs_3 = buffer.data(fs + 3);
    const auto *fs_4 = buffer.data(fs + 4);
    const auto *fs_5 = buffer.data(fs + 5);
    const auto *fs_6 = buffer.data(fs + 6);
    const auto *fs_7 = buffer.data(fs + 7);
    const auto *fs_8 = buffer.data(fs + 8);
    const auto *fs_9 = buffer.data(fs + 9);

    const auto *hs_2 = buffer.data(hs + 2);
    const auto *hs_4 = buffer.data(hs + 4);
    const auto *hs_5 = buffer.data(hs + 5);
    const auto *hs_7 = buffer.data(hs + 7);
    const auto *hs_8 = buffer.data(hs + 8);
    const auto *hs_9 = buffer.data(hs + 9);
    const auto *hs_11 = buffer.data(hs + 11);
    const auto *hs_12 = buffer.data(hs + 12);
    const auto *hs_13 = buffer.data(hs + 13);
    const auto *hs_14 = buffer.data(hs + 14);
    const auto *hs_16 = buffer.data(hs + 16);
    const auto *hs_17 = buffer.data(hs + 17);
    const auto *hs_18 = buffer.data(hs + 18);
    const auto *hs_19 = buffer.data(hs + 19);
    const auto *hs_20 = buffer.data(hs + 20);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, fs_0, fs_1, fs_2, hs_2, hs_4, hs_5, \
                         hs_7, hs_8, hs_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hs_2[k];

        t_1[k] = f_0 * hs_4[k];

        t_2[k] = -fs_0[k]
                 + f_0 * hs_5[k];

        t_3[k] = f_0 * hs_7[k];

        t_4[k] = -fs_1[k]
                 + f_0 * hs_8[k];

        t_5[k] = -2.0 * fs_2[k]
                 + f_0 * hs_9[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, t_11, fs_3, fs_4, fs_5, fs_6, hs_11, hs_12, \
                         hs_13, hs_14, hs_16, hs_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * hs_11[k];

        t_7[k] = -fs_3[k]
                 + f_0 * hs_12[k];

        t_8[k] = -2.0 * fs_4[k]
                 + f_0 * hs_13[k];

        t_9[k] = -3.0 * fs_5[k]
                 + f_0 * hs_14[k];

        t_10[k] = f_0 * hs_16[k];

        t_11[k] = -fs_6[k]
                  + f_0 * hs_17[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, fs_7, fs_8, fs_9, hs_18, hs_19, \
                         hs_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = -2.0 * fs_7[k]
                  + f_0 * hs_18[k];

        t_13[k] = -3.0 * fs_8[k]
                  + f_0 * hs_19[k];

        t_14[k] = -4.0 * fs_9[k]
                  + f_0 * hs_20[k];
    }
}

}  // namespace simdt2ceri
