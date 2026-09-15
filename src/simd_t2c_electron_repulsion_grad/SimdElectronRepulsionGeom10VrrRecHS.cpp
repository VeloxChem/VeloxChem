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


#include "SimdElectronRepulsionGeom10VrrRecHS.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_geom_10_hs_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t gs, const size_t is,
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
    auto *t_15 = buffer.data(target + 15);
    auto *t_16 = buffer.data(target + 16);
    auto *t_17 = buffer.data(target + 17);
    auto *t_18 = buffer.data(target + 18);
    auto *t_19 = buffer.data(target + 19);
    auto *t_20 = buffer.data(target + 20);

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
    const auto *gs_10 = buffer.data(gs + 10);
    const auto *gs_11 = buffer.data(gs + 11);
    const auto *gs_12 = buffer.data(gs + 12);
    const auto *gs_13 = buffer.data(gs + 13);
    const auto *gs_14 = buffer.data(gs + 14);

    const auto *is_0 = buffer.data(is + 0);
    const auto *is_1 = buffer.data(is + 1);
    const auto *is_2 = buffer.data(is + 2);
    const auto *is_3 = buffer.data(is + 3);
    const auto *is_4 = buffer.data(is + 4);
    const auto *is_5 = buffer.data(is + 5);
    const auto *is_6 = buffer.data(is + 6);
    const auto *is_7 = buffer.data(is + 7);
    const auto *is_8 = buffer.data(is + 8);
    const auto *is_9 = buffer.data(is + 9);
    const auto *is_10 = buffer.data(is + 10);
    const auto *is_11 = buffer.data(is + 11);
    const auto *is_12 = buffer.data(is + 12);
    const auto *is_13 = buffer.data(is + 13);
    const auto *is_14 = buffer.data(is + 14);
    const auto *is_15 = buffer.data(is + 15);
    const auto *is_16 = buffer.data(is + 16);
    const auto *is_17 = buffer.data(is + 17);
    const auto *is_18 = buffer.data(is + 18);
    const auto *is_19 = buffer.data(is + 19);
    const auto *is_20 = buffer.data(is + 20);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, gs_0, gs_1, gs_2, gs_3, gs_4, is_0, is_1, \
                         is_2, is_3, is_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -5.0 * gs_0[k]
                 + f_0 * is_0[k];

        t_1[k] = -4.0 * gs_1[k]
                 + f_0 * is_1[k];

        t_2[k] = -4.0 * gs_2[k]
                 + f_0 * is_2[k];

        t_3[k] = -3.0 * gs_3[k]
                 + f_0 * is_3[k];

        t_4[k] = -3.0 * gs_4[k]
                 + f_0 * is_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, gs_5, gs_6, gs_7, gs_8, gs_9, is_5, is_6, \
                         is_7, is_8, is_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -3.0 * gs_5[k]
                 + f_0 * is_5[k];

        t_6[k] = -2.0 * gs_6[k]
                 + f_0 * is_6[k];

        t_7[k] = -2.0 * gs_7[k]
                 + f_0 * is_7[k];

        t_8[k] = -2.0 * gs_8[k]
                 + f_0 * is_8[k];

        t_9[k] = -2.0 * gs_9[k]
                 + f_0 * is_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, gs_10, gs_11, gs_12, gs_13, gs_14, \
                         is_10, is_11, is_12, is_13, is_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -gs_10[k]
                  + f_0 * is_10[k];

        t_11[k] = -gs_11[k]
                  + f_0 * is_11[k];

        t_12[k] = -gs_12[k]
                  + f_0 * is_12[k];

        t_13[k] = -gs_13[k]
                  + f_0 * is_13[k];

        t_14[k] = -gs_14[k]
                  + f_0 * is_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, t_20, is_15, is_16, is_17, is_18, \
                         is_19, is_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_0 * is_15[k];

        t_16[k] = f_0 * is_16[k];

        t_17[k] = f_0 * is_17[k];

        t_18[k] = f_0 * is_18[k];

        t_19[k] = f_0 * is_19[k];

        t_20[k] = f_0 * is_20[k];
    }
}

auto
compute_prim_geom_10_hs_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t gs, const size_t is,
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
    auto *t_15 = buffer.data(target + 15);
    auto *t_16 = buffer.data(target + 16);
    auto *t_17 = buffer.data(target + 17);
    auto *t_18 = buffer.data(target + 18);
    auto *t_19 = buffer.data(target + 19);
    auto *t_20 = buffer.data(target + 20);

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
    const auto *gs_10 = buffer.data(gs + 10);
    const auto *gs_11 = buffer.data(gs + 11);
    const auto *gs_12 = buffer.data(gs + 12);
    const auto *gs_13 = buffer.data(gs + 13);
    const auto *gs_14 = buffer.data(gs + 14);

    const auto *is_1 = buffer.data(is + 1);
    const auto *is_3 = buffer.data(is + 3);
    const auto *is_4 = buffer.data(is + 4);
    const auto *is_6 = buffer.data(is + 6);
    const auto *is_7 = buffer.data(is + 7);
    const auto *is_8 = buffer.data(is + 8);
    const auto *is_10 = buffer.data(is + 10);
    const auto *is_11 = buffer.data(is + 11);
    const auto *is_12 = buffer.data(is + 12);
    const auto *is_13 = buffer.data(is + 13);
    const auto *is_15 = buffer.data(is + 15);
    const auto *is_16 = buffer.data(is + 16);
    const auto *is_17 = buffer.data(is + 17);
    const auto *is_18 = buffer.data(is + 18);
    const auto *is_19 = buffer.data(is + 19);
    const auto *is_21 = buffer.data(is + 21);
    const auto *is_22 = buffer.data(is + 22);
    const auto *is_23 = buffer.data(is + 23);
    const auto *is_24 = buffer.data(is + 24);
    const auto *is_25 = buffer.data(is + 25);
    const auto *is_26 = buffer.data(is + 26);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, gs_0, gs_1, gs_2, is_1, is_3, is_4, \
                         is_6, is_7, is_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * is_1[k];

        t_1[k] = -gs_0[k]
                 + f_0 * is_3[k];

        t_2[k] = f_0 * is_4[k];

        t_3[k] = -2.0 * gs_1[k]
                 + f_0 * is_6[k];

        t_4[k] = -gs_2[k]
                 + f_0 * is_7[k];

        t_5[k] = f_0 * is_8[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, gs_3, gs_4, gs_5, gs_6, is_10, is_11, \
                         is_12, is_13, is_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -3.0 * gs_3[k]
                 + f_0 * is_10[k];

        t_7[k] = -2.0 * gs_4[k]
                 + f_0 * is_11[k];

        t_8[k] = -gs_5[k]
                 + f_0 * is_12[k];

        t_9[k] = f_0 * is_13[k];

        t_10[k] = -4.0 * gs_6[k]
                  + f_0 * is_15[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, gs_7, gs_8, gs_9, gs_10, is_16, is_17, \
                         is_18, is_19, is_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = -3.0 * gs_7[k]
                  + f_0 * is_16[k];

        t_12[k] = -2.0 * gs_8[k]
                  + f_0 * is_17[k];

        t_13[k] = -gs_9[k]
                  + f_0 * is_18[k];

        t_14[k] = f_0 * is_19[k];

        t_15[k] = -5.0 * gs_10[k]
                  + f_0 * is_21[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, gs_11, gs_12, gs_13, gs_14, is_22, \
                         is_23, is_24, is_25, is_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = -4.0 * gs_11[k]
                  + f_0 * is_22[k];

        t_17[k] = -3.0 * gs_12[k]
                  + f_0 * is_23[k];

        t_18[k] = -2.0 * gs_13[k]
                  + f_0 * is_24[k];

        t_19[k] = -gs_14[k]
                  + f_0 * is_25[k];

        t_20[k] = f_0 * is_26[k];
    }
}

auto
compute_prim_geom_10_hs_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t gs, const size_t is,
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
    auto *t_15 = buffer.data(target + 15);
    auto *t_16 = buffer.data(target + 16);
    auto *t_17 = buffer.data(target + 17);
    auto *t_18 = buffer.data(target + 18);
    auto *t_19 = buffer.data(target + 19);
    auto *t_20 = buffer.data(target + 20);

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
    const auto *gs_10 = buffer.data(gs + 10);
    const auto *gs_11 = buffer.data(gs + 11);
    const auto *gs_12 = buffer.data(gs + 12);
    const auto *gs_13 = buffer.data(gs + 13);
    const auto *gs_14 = buffer.data(gs + 14);

    const auto *is_2 = buffer.data(is + 2);
    const auto *is_4 = buffer.data(is + 4);
    const auto *is_5 = buffer.data(is + 5);
    const auto *is_7 = buffer.data(is + 7);
    const auto *is_8 = buffer.data(is + 8);
    const auto *is_9 = buffer.data(is + 9);
    const auto *is_11 = buffer.data(is + 11);
    const auto *is_12 = buffer.data(is + 12);
    const auto *is_13 = buffer.data(is + 13);
    const auto *is_14 = buffer.data(is + 14);
    const auto *is_16 = buffer.data(is + 16);
    const auto *is_17 = buffer.data(is + 17);
    const auto *is_18 = buffer.data(is + 18);
    const auto *is_19 = buffer.data(is + 19);
    const auto *is_20 = buffer.data(is + 20);
    const auto *is_22 = buffer.data(is + 22);
    const auto *is_23 = buffer.data(is + 23);
    const auto *is_24 = buffer.data(is + 24);
    const auto *is_25 = buffer.data(is + 25);
    const auto *is_26 = buffer.data(is + 26);
    const auto *is_27 = buffer.data(is + 27);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, gs_0, gs_1, gs_2, is_2, is_4, is_5, \
                         is_7, is_8, is_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * is_2[k];

        t_1[k] = f_0 * is_4[k];

        t_2[k] = -gs_0[k]
                 + f_0 * is_5[k];

        t_3[k] = f_0 * is_7[k];

        t_4[k] = -gs_1[k]
                 + f_0 * is_8[k];

        t_5[k] = -2.0 * gs_2[k]
                 + f_0 * is_9[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, t_11, gs_3, gs_4, gs_5, gs_6, is_11, is_12, \
                         is_13, is_14, is_16, is_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * is_11[k];

        t_7[k] = -gs_3[k]
                 + f_0 * is_12[k];

        t_8[k] = -2.0 * gs_4[k]
                 + f_0 * is_13[k];

        t_9[k] = -3.0 * gs_5[k]
                 + f_0 * is_14[k];

        t_10[k] = f_0 * is_16[k];

        t_11[k] = -gs_6[k]
                  + f_0 * is_17[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, gs_7, gs_8, gs_9, gs_10, is_18, is_19, \
                         is_20, is_22, is_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = -2.0 * gs_7[k]
                  + f_0 * is_18[k];

        t_13[k] = -3.0 * gs_8[k]
                  + f_0 * is_19[k];

        t_14[k] = -4.0 * gs_9[k]
                  + f_0 * is_20[k];

        t_15[k] = f_0 * is_22[k];

        t_16[k] = -gs_10[k]
                  + f_0 * is_23[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, gs_11, gs_12, gs_13, gs_14, is_24, is_25, \
                         is_26, is_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = -2.0 * gs_11[k]
                  + f_0 * is_24[k];

        t_18[k] = -3.0 * gs_12[k]
                  + f_0 * is_25[k];

        t_19[k] = -4.0 * gs_13[k]
                  + f_0 * is_26[k];

        t_20[k] = -5.0 * gs_14[k]
                  + f_0 * is_27[k];
    }
}

}  // namespace simdt2ceri
