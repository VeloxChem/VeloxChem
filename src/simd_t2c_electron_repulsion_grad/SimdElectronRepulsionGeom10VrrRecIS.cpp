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


#include "SimdElectronRepulsionGeom10VrrRecIS.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_geom_10_is_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t hs, const size_t ks,
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
    auto *t_21 = buffer.data(target + 21);
    auto *t_22 = buffer.data(target + 22);
    auto *t_23 = buffer.data(target + 23);
    auto *t_24 = buffer.data(target + 24);
    auto *t_25 = buffer.data(target + 25);
    auto *t_26 = buffer.data(target + 26);
    auto *t_27 = buffer.data(target + 27);

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
    const auto *hs_15 = buffer.data(hs + 15);
    const auto *hs_16 = buffer.data(hs + 16);
    const auto *hs_17 = buffer.data(hs + 17);
    const auto *hs_18 = buffer.data(hs + 18);
    const auto *hs_19 = buffer.data(hs + 19);
    const auto *hs_20 = buffer.data(hs + 20);

    const auto *ks_0 = buffer.data(ks + 0);
    const auto *ks_1 = buffer.data(ks + 1);
    const auto *ks_2 = buffer.data(ks + 2);
    const auto *ks_3 = buffer.data(ks + 3);
    const auto *ks_4 = buffer.data(ks + 4);
    const auto *ks_5 = buffer.data(ks + 5);
    const auto *ks_6 = buffer.data(ks + 6);
    const auto *ks_7 = buffer.data(ks + 7);
    const auto *ks_8 = buffer.data(ks + 8);
    const auto *ks_9 = buffer.data(ks + 9);
    const auto *ks_10 = buffer.data(ks + 10);
    const auto *ks_11 = buffer.data(ks + 11);
    const auto *ks_12 = buffer.data(ks + 12);
    const auto *ks_13 = buffer.data(ks + 13);
    const auto *ks_14 = buffer.data(ks + 14);
    const auto *ks_15 = buffer.data(ks + 15);
    const auto *ks_16 = buffer.data(ks + 16);
    const auto *ks_17 = buffer.data(ks + 17);
    const auto *ks_18 = buffer.data(ks + 18);
    const auto *ks_19 = buffer.data(ks + 19);
    const auto *ks_20 = buffer.data(ks + 20);
    const auto *ks_21 = buffer.data(ks + 21);
    const auto *ks_22 = buffer.data(ks + 22);
    const auto *ks_23 = buffer.data(ks + 23);
    const auto *ks_24 = buffer.data(ks + 24);
    const auto *ks_25 = buffer.data(ks + 25);
    const auto *ks_26 = buffer.data(ks + 26);
    const auto *ks_27 = buffer.data(ks + 27);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, hs_0, hs_1, hs_2, hs_3, hs_4, ks_0, ks_1, \
                         ks_2, ks_3, ks_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -6.0 * hs_0[k]
                 + f_0 * ks_0[k];

        t_1[k] = -5.0 * hs_1[k]
                 + f_0 * ks_1[k];

        t_2[k] = -5.0 * hs_2[k]
                 + f_0 * ks_2[k];

        t_3[k] = -4.0 * hs_3[k]
                 + f_0 * ks_3[k];

        t_4[k] = -4.0 * hs_4[k]
                 + f_0 * ks_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, hs_5, hs_6, hs_7, hs_8, hs_9, ks_5, ks_6, \
                         ks_7, ks_8, ks_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -4.0 * hs_5[k]
                 + f_0 * ks_5[k];

        t_6[k] = -3.0 * hs_6[k]
                 + f_0 * ks_6[k];

        t_7[k] = -3.0 * hs_7[k]
                 + f_0 * ks_7[k];

        t_8[k] = -3.0 * hs_8[k]
                 + f_0 * ks_8[k];

        t_9[k] = -3.0 * hs_9[k]
                 + f_0 * ks_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, hs_10, hs_11, hs_12, hs_13, hs_14, \
                         ks_10, ks_11, ks_12, ks_13, ks_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -2.0 * hs_10[k]
                  + f_0 * ks_10[k];

        t_11[k] = -2.0 * hs_11[k]
                  + f_0 * ks_11[k];

        t_12[k] = -2.0 * hs_12[k]
                  + f_0 * ks_12[k];

        t_13[k] = -2.0 * hs_13[k]
                  + f_0 * ks_13[k];

        t_14[k] = -2.0 * hs_14[k]
                  + f_0 * ks_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, hs_15, hs_16, hs_17, hs_18, hs_19, \
                         ks_15, ks_16, ks_17, ks_18, ks_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -hs_15[k]
                  + f_0 * ks_15[k];

        t_16[k] = -hs_16[k]
                  + f_0 * ks_16[k];

        t_17[k] = -hs_17[k]
                  + f_0 * ks_17[k];

        t_18[k] = -hs_18[k]
                  + f_0 * ks_18[k];

        t_19[k] = -hs_19[k]
                  + f_0 * ks_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, t_25, t_26, hs_20, ks_20, ks_21, ks_22, \
                         ks_23, ks_24, ks_25, ks_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -hs_20[k]
                  + f_0 * ks_20[k];

        t_21[k] = f_0 * ks_21[k];

        t_22[k] = f_0 * ks_22[k];

        t_23[k] = f_0 * ks_23[k];

        t_24[k] = f_0 * ks_24[k];

        t_25[k] = f_0 * ks_25[k];

        t_26[k] = f_0 * ks_26[k];
    }

#pragma omp simd aligned(t_27, ks_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_0 * ks_27[k];
    }
}

auto
compute_prim_geom_10_is_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t hs, const size_t ks,
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
    auto *t_21 = buffer.data(target + 21);
    auto *t_22 = buffer.data(target + 22);
    auto *t_23 = buffer.data(target + 23);
    auto *t_24 = buffer.data(target + 24);
    auto *t_25 = buffer.data(target + 25);
    auto *t_26 = buffer.data(target + 26);
    auto *t_27 = buffer.data(target + 27);

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
    const auto *hs_15 = buffer.data(hs + 15);
    const auto *hs_16 = buffer.data(hs + 16);
    const auto *hs_17 = buffer.data(hs + 17);
    const auto *hs_18 = buffer.data(hs + 18);
    const auto *hs_19 = buffer.data(hs + 19);
    const auto *hs_20 = buffer.data(hs + 20);

    const auto *ks_1 = buffer.data(ks + 1);
    const auto *ks_3 = buffer.data(ks + 3);
    const auto *ks_4 = buffer.data(ks + 4);
    const auto *ks_6 = buffer.data(ks + 6);
    const auto *ks_7 = buffer.data(ks + 7);
    const auto *ks_8 = buffer.data(ks + 8);
    const auto *ks_10 = buffer.data(ks + 10);
    const auto *ks_11 = buffer.data(ks + 11);
    const auto *ks_12 = buffer.data(ks + 12);
    const auto *ks_13 = buffer.data(ks + 13);
    const auto *ks_15 = buffer.data(ks + 15);
    const auto *ks_16 = buffer.data(ks + 16);
    const auto *ks_17 = buffer.data(ks + 17);
    const auto *ks_18 = buffer.data(ks + 18);
    const auto *ks_19 = buffer.data(ks + 19);
    const auto *ks_21 = buffer.data(ks + 21);
    const auto *ks_22 = buffer.data(ks + 22);
    const auto *ks_23 = buffer.data(ks + 23);
    const auto *ks_24 = buffer.data(ks + 24);
    const auto *ks_25 = buffer.data(ks + 25);
    const auto *ks_26 = buffer.data(ks + 26);
    const auto *ks_28 = buffer.data(ks + 28);
    const auto *ks_29 = buffer.data(ks + 29);
    const auto *ks_30 = buffer.data(ks + 30);
    const auto *ks_31 = buffer.data(ks + 31);
    const auto *ks_32 = buffer.data(ks + 32);
    const auto *ks_33 = buffer.data(ks + 33);
    const auto *ks_34 = buffer.data(ks + 34);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, hs_0, hs_1, hs_2, ks_1, ks_3, ks_4, \
                         ks_6, ks_7, ks_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ks_1[k];

        t_1[k] = -hs_0[k]
                 + f_0 * ks_3[k];

        t_2[k] = f_0 * ks_4[k];

        t_3[k] = -2.0 * hs_1[k]
                 + f_0 * ks_6[k];

        t_4[k] = -hs_2[k]
                 + f_0 * ks_7[k];

        t_5[k] = f_0 * ks_8[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, hs_3, hs_4, hs_5, hs_6, ks_10, ks_11, \
                         ks_12, ks_13, ks_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -3.0 * hs_3[k]
                 + f_0 * ks_10[k];

        t_7[k] = -2.0 * hs_4[k]
                 + f_0 * ks_11[k];

        t_8[k] = -hs_5[k]
                 + f_0 * ks_12[k];

        t_9[k] = f_0 * ks_13[k];

        t_10[k] = -4.0 * hs_6[k]
                  + f_0 * ks_15[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, hs_7, hs_8, hs_9, hs_10, ks_16, ks_17, \
                         ks_18, ks_19, ks_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = -3.0 * hs_7[k]
                  + f_0 * ks_16[k];

        t_12[k] = -2.0 * hs_8[k]
                  + f_0 * ks_17[k];

        t_13[k] = -hs_9[k]
                  + f_0 * ks_18[k];

        t_14[k] = f_0 * ks_19[k];

        t_15[k] = -5.0 * hs_10[k]
                  + f_0 * ks_21[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, hs_11, hs_12, hs_13, hs_14, ks_22, \
                         ks_23, ks_24, ks_25, ks_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = -4.0 * hs_11[k]
                  + f_0 * ks_22[k];

        t_17[k] = -3.0 * hs_12[k]
                  + f_0 * ks_23[k];

        t_18[k] = -2.0 * hs_13[k]
                  + f_0 * ks_24[k];

        t_19[k] = -hs_14[k]
                  + f_0 * ks_25[k];

        t_20[k] = f_0 * ks_26[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, hs_15, hs_16, hs_17, hs_18, hs_19, \
                         ks_28, ks_29, ks_30, ks_31, ks_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = -6.0 * hs_15[k]
                  + f_0 * ks_28[k];

        t_22[k] = -5.0 * hs_16[k]
                  + f_0 * ks_29[k];

        t_23[k] = -4.0 * hs_17[k]
                  + f_0 * ks_30[k];

        t_24[k] = -3.0 * hs_18[k]
                  + f_0 * ks_31[k];

        t_25[k] = -2.0 * hs_19[k]
                  + f_0 * ks_32[k];
    }

#pragma omp simd aligned(t_26, t_27, hs_20, ks_33, ks_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = -hs_20[k]
                  + f_0 * ks_33[k];

        t_27[k] = f_0 * ks_34[k];
    }
}

auto
compute_prim_geom_10_is_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t hs, const size_t ks,
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
    auto *t_21 = buffer.data(target + 21);
    auto *t_22 = buffer.data(target + 22);
    auto *t_23 = buffer.data(target + 23);
    auto *t_24 = buffer.data(target + 24);
    auto *t_25 = buffer.data(target + 25);
    auto *t_26 = buffer.data(target + 26);
    auto *t_27 = buffer.data(target + 27);

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
    const auto *hs_15 = buffer.data(hs + 15);
    const auto *hs_16 = buffer.data(hs + 16);
    const auto *hs_17 = buffer.data(hs + 17);
    const auto *hs_18 = buffer.data(hs + 18);
    const auto *hs_19 = buffer.data(hs + 19);
    const auto *hs_20 = buffer.data(hs + 20);

    const auto *ks_2 = buffer.data(ks + 2);
    const auto *ks_4 = buffer.data(ks + 4);
    const auto *ks_5 = buffer.data(ks + 5);
    const auto *ks_7 = buffer.data(ks + 7);
    const auto *ks_8 = buffer.data(ks + 8);
    const auto *ks_9 = buffer.data(ks + 9);
    const auto *ks_11 = buffer.data(ks + 11);
    const auto *ks_12 = buffer.data(ks + 12);
    const auto *ks_13 = buffer.data(ks + 13);
    const auto *ks_14 = buffer.data(ks + 14);
    const auto *ks_16 = buffer.data(ks + 16);
    const auto *ks_17 = buffer.data(ks + 17);
    const auto *ks_18 = buffer.data(ks + 18);
    const auto *ks_19 = buffer.data(ks + 19);
    const auto *ks_20 = buffer.data(ks + 20);
    const auto *ks_22 = buffer.data(ks + 22);
    const auto *ks_23 = buffer.data(ks + 23);
    const auto *ks_24 = buffer.data(ks + 24);
    const auto *ks_25 = buffer.data(ks + 25);
    const auto *ks_26 = buffer.data(ks + 26);
    const auto *ks_27 = buffer.data(ks + 27);
    const auto *ks_29 = buffer.data(ks + 29);
    const auto *ks_30 = buffer.data(ks + 30);
    const auto *ks_31 = buffer.data(ks + 31);
    const auto *ks_32 = buffer.data(ks + 32);
    const auto *ks_33 = buffer.data(ks + 33);
    const auto *ks_34 = buffer.data(ks + 34);
    const auto *ks_35 = buffer.data(ks + 35);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, hs_0, hs_1, hs_2, ks_2, ks_4, ks_5, \
                         ks_7, ks_8, ks_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ks_2[k];

        t_1[k] = f_0 * ks_4[k];

        t_2[k] = -hs_0[k]
                 + f_0 * ks_5[k];

        t_3[k] = f_0 * ks_7[k];

        t_4[k] = -hs_1[k]
                 + f_0 * ks_8[k];

        t_5[k] = -2.0 * hs_2[k]
                 + f_0 * ks_9[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, t_11, hs_3, hs_4, hs_5, hs_6, ks_11, ks_12, \
                         ks_13, ks_14, ks_16, ks_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * ks_11[k];

        t_7[k] = -hs_3[k]
                 + f_0 * ks_12[k];

        t_8[k] = -2.0 * hs_4[k]
                 + f_0 * ks_13[k];

        t_9[k] = -3.0 * hs_5[k]
                 + f_0 * ks_14[k];

        t_10[k] = f_0 * ks_16[k];

        t_11[k] = -hs_6[k]
                  + f_0 * ks_17[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, hs_7, hs_8, hs_9, hs_10, ks_18, ks_19, \
                         ks_20, ks_22, ks_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = -2.0 * hs_7[k]
                  + f_0 * ks_18[k];

        t_13[k] = -3.0 * hs_8[k]
                  + f_0 * ks_19[k];

        t_14[k] = -4.0 * hs_9[k]
                  + f_0 * ks_20[k];

        t_15[k] = f_0 * ks_22[k];

        t_16[k] = -hs_10[k]
                  + f_0 * ks_23[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, t_21, hs_11, hs_12, hs_13, hs_14, ks_24, \
                         ks_25, ks_26, ks_27, ks_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = -2.0 * hs_11[k]
                  + f_0 * ks_24[k];

        t_18[k] = -3.0 * hs_12[k]
                  + f_0 * ks_25[k];

        t_19[k] = -4.0 * hs_13[k]
                  + f_0 * ks_26[k];

        t_20[k] = -5.0 * hs_14[k]
                  + f_0 * ks_27[k];

        t_21[k] = f_0 * ks_29[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, hs_15, hs_16, hs_17, hs_18, hs_19, \
                         ks_30, ks_31, ks_32, ks_33, ks_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = -hs_15[k]
                  + f_0 * ks_30[k];

        t_23[k] = -2.0 * hs_16[k]
                  + f_0 * ks_31[k];

        t_24[k] = -3.0 * hs_17[k]
                  + f_0 * ks_32[k];

        t_25[k] = -4.0 * hs_18[k]
                  + f_0 * ks_33[k];

        t_26[k] = -5.0 * hs_19[k]
                  + f_0 * ks_34[k];
    }

#pragma omp simd aligned(t_27, hs_20, ks_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = -6.0 * hs_20[k]
                  + f_0 * ks_35[k];
    }
}

}  // namespace simdt2ceri
