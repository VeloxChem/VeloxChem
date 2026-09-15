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


#include "SimdElectronRepulsionGeom10VrrRecFP.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_geom_10_fp_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t dp, const size_t gp,
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
    auto *t_28 = buffer.data(target + 28);
    auto *t_29 = buffer.data(target + 29);

    const auto *dp_0 = buffer.data(dp + 0);
    const auto *dp_1 = buffer.data(dp + 1);
    const auto *dp_2 = buffer.data(dp + 2);
    const auto *dp_3 = buffer.data(dp + 3);
    const auto *dp_4 = buffer.data(dp + 4);
    const auto *dp_5 = buffer.data(dp + 5);
    const auto *dp_6 = buffer.data(dp + 6);
    const auto *dp_7 = buffer.data(dp + 7);
    const auto *dp_8 = buffer.data(dp + 8);
    const auto *dp_9 = buffer.data(dp + 9);
    const auto *dp_10 = buffer.data(dp + 10);
    const auto *dp_11 = buffer.data(dp + 11);
    const auto *dp_12 = buffer.data(dp + 12);
    const auto *dp_13 = buffer.data(dp + 13);
    const auto *dp_14 = buffer.data(dp + 14);
    const auto *dp_15 = buffer.data(dp + 15);
    const auto *dp_16 = buffer.data(dp + 16);
    const auto *dp_17 = buffer.data(dp + 17);

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_1 = buffer.data(gp + 1);
    const auto *gp_2 = buffer.data(gp + 2);
    const auto *gp_3 = buffer.data(gp + 3);
    const auto *gp_4 = buffer.data(gp + 4);
    const auto *gp_5 = buffer.data(gp + 5);
    const auto *gp_6 = buffer.data(gp + 6);
    const auto *gp_7 = buffer.data(gp + 7);
    const auto *gp_8 = buffer.data(gp + 8);
    const auto *gp_9 = buffer.data(gp + 9);
    const auto *gp_10 = buffer.data(gp + 10);
    const auto *gp_11 = buffer.data(gp + 11);
    const auto *gp_12 = buffer.data(gp + 12);
    const auto *gp_13 = buffer.data(gp + 13);
    const auto *gp_14 = buffer.data(gp + 14);
    const auto *gp_15 = buffer.data(gp + 15);
    const auto *gp_16 = buffer.data(gp + 16);
    const auto *gp_17 = buffer.data(gp + 17);
    const auto *gp_18 = buffer.data(gp + 18);
    const auto *gp_19 = buffer.data(gp + 19);
    const auto *gp_20 = buffer.data(gp + 20);
    const auto *gp_21 = buffer.data(gp + 21);
    const auto *gp_22 = buffer.data(gp + 22);
    const auto *gp_23 = buffer.data(gp + 23);
    const auto *gp_24 = buffer.data(gp + 24);
    const auto *gp_25 = buffer.data(gp + 25);
    const auto *gp_26 = buffer.data(gp + 26);
    const auto *gp_27 = buffer.data(gp + 27);
    const auto *gp_28 = buffer.data(gp + 28);
    const auto *gp_29 = buffer.data(gp + 29);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, dp_0, dp_1, dp_2, dp_3, dp_4, gp_0, gp_1, \
                         gp_2, gp_3, gp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -3.0 * dp_0[k]
                 + f_0 * gp_0[k];

        t_1[k] = -3.0 * dp_1[k]
                 + f_0 * gp_1[k];

        t_2[k] = -3.0 * dp_2[k]
                 + f_0 * gp_2[k];

        t_3[k] = -2.0 * dp_3[k]
                 + f_0 * gp_3[k];

        t_4[k] = -2.0 * dp_4[k]
                 + f_0 * gp_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, dp_5, dp_6, dp_7, dp_8, dp_9, gp_5, gp_6, \
                         gp_7, gp_8, gp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -2.0 * dp_5[k]
                 + f_0 * gp_5[k];

        t_6[k] = -2.0 * dp_6[k]
                 + f_0 * gp_6[k];

        t_7[k] = -2.0 * dp_7[k]
                 + f_0 * gp_7[k];

        t_8[k] = -2.0 * dp_8[k]
                 + f_0 * gp_8[k];

        t_9[k] = -dp_9[k]
                 + f_0 * gp_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, dp_10, dp_11, dp_12, dp_13, dp_14, \
                         gp_10, gp_11, gp_12, gp_13, gp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -dp_10[k]
                  + f_0 * gp_10[k];

        t_11[k] = -dp_11[k]
                  + f_0 * gp_11[k];

        t_12[k] = -dp_12[k]
                  + f_0 * gp_12[k];

        t_13[k] = -dp_13[k]
                  + f_0 * gp_13[k];

        t_14[k] = -dp_14[k]
                  + f_0 * gp_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, t_20, dp_15, dp_16, dp_17, gp_15, \
                         gp_16, gp_17, gp_18, gp_19, gp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -dp_15[k]
                  + f_0 * gp_15[k];

        t_16[k] = -dp_16[k]
                  + f_0 * gp_16[k];

        t_17[k] = -dp_17[k]
                  + f_0 * gp_17[k];

        t_18[k] = f_0 * gp_18[k];

        t_19[k] = f_0 * gp_19[k];

        t_20[k] = f_0 * gp_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, t_26, t_27, t_28, gp_21, gp_22, gp_23, \
                         gp_24, gp_25, gp_26, gp_27, gp_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_0 * gp_21[k];

        t_22[k] = f_0 * gp_22[k];

        t_23[k] = f_0 * gp_23[k];

        t_24[k] = f_0 * gp_24[k];

        t_25[k] = f_0 * gp_25[k];

        t_26[k] = f_0 * gp_26[k];

        t_27[k] = f_0 * gp_27[k];

        t_28[k] = f_0 * gp_28[k];
    }

#pragma omp simd aligned(t_29, gp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_0 * gp_29[k];
    }
}

auto
compute_prim_geom_10_fp_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t dp, const size_t gp,
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
    auto *t_28 = buffer.data(target + 28);
    auto *t_29 = buffer.data(target + 29);

    const auto *dp_0 = buffer.data(dp + 0);
    const auto *dp_1 = buffer.data(dp + 1);
    const auto *dp_2 = buffer.data(dp + 2);
    const auto *dp_3 = buffer.data(dp + 3);
    const auto *dp_4 = buffer.data(dp + 4);
    const auto *dp_5 = buffer.data(dp + 5);
    const auto *dp_6 = buffer.data(dp + 6);
    const auto *dp_7 = buffer.data(dp + 7);
    const auto *dp_8 = buffer.data(dp + 8);
    const auto *dp_9 = buffer.data(dp + 9);
    const auto *dp_10 = buffer.data(dp + 10);
    const auto *dp_11 = buffer.data(dp + 11);
    const auto *dp_12 = buffer.data(dp + 12);
    const auto *dp_13 = buffer.data(dp + 13);
    const auto *dp_14 = buffer.data(dp + 14);
    const auto *dp_15 = buffer.data(dp + 15);
    const auto *dp_16 = buffer.data(dp + 16);
    const auto *dp_17 = buffer.data(dp + 17);

    const auto *gp_3 = buffer.data(gp + 3);
    const auto *gp_4 = buffer.data(gp + 4);
    const auto *gp_5 = buffer.data(gp + 5);
    const auto *gp_9 = buffer.data(gp + 9);
    const auto *gp_10 = buffer.data(gp + 10);
    const auto *gp_11 = buffer.data(gp + 11);
    const auto *gp_12 = buffer.data(gp + 12);
    const auto *gp_13 = buffer.data(gp + 13);
    const auto *gp_14 = buffer.data(gp + 14);
    const auto *gp_18 = buffer.data(gp + 18);
    const auto *gp_19 = buffer.data(gp + 19);
    const auto *gp_20 = buffer.data(gp + 20);
    const auto *gp_21 = buffer.data(gp + 21);
    const auto *gp_22 = buffer.data(gp + 22);
    const auto *gp_23 = buffer.data(gp + 23);
    const auto *gp_24 = buffer.data(gp + 24);
    const auto *gp_25 = buffer.data(gp + 25);
    const auto *gp_26 = buffer.data(gp + 26);
    const auto *gp_30 = buffer.data(gp + 30);
    const auto *gp_31 = buffer.data(gp + 31);
    const auto *gp_32 = buffer.data(gp + 32);
    const auto *gp_33 = buffer.data(gp + 33);
    const auto *gp_34 = buffer.data(gp + 34);
    const auto *gp_35 = buffer.data(gp + 35);
    const auto *gp_36 = buffer.data(gp + 36);
    const auto *gp_37 = buffer.data(gp + 37);
    const auto *gp_38 = buffer.data(gp + 38);
    const auto *gp_39 = buffer.data(gp + 39);
    const auto *gp_40 = buffer.data(gp + 40);
    const auto *gp_41 = buffer.data(gp + 41);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, dp_0, dp_1, dp_2, gp_3, gp_4, gp_5, \
                         gp_9, gp_10, gp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gp_3[k];

        t_1[k] = f_0 * gp_4[k];

        t_2[k] = f_0 * gp_5[k];

        t_3[k] = -dp_0[k]
                 + f_0 * gp_9[k];

        t_4[k] = -dp_1[k]
                 + f_0 * gp_10[k];

        t_5[k] = -dp_2[k]
                 + f_0 * gp_11[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, t_11, dp_3, dp_4, dp_5, gp_12, gp_13, \
                         gp_14, gp_18, gp_19, gp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * gp_12[k];

        t_7[k] = f_0 * gp_13[k];

        t_8[k] = f_0 * gp_14[k];

        t_9[k] = -2.0 * dp_3[k]
                 + f_0 * gp_18[k];

        t_10[k] = -2.0 * dp_4[k]
                  + f_0 * gp_19[k];

        t_11[k] = -2.0 * dp_5[k]
                  + f_0 * gp_20[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, t_17, dp_6, dp_7, dp_8, gp_21, gp_22, \
                         gp_23, gp_24, gp_25, gp_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = -dp_6[k]
                  + f_0 * gp_21[k];

        t_13[k] = -dp_7[k]
                  + f_0 * gp_22[k];

        t_14[k] = -dp_8[k]
                  + f_0 * gp_23[k];

        t_15[k] = f_0 * gp_24[k];

        t_16[k] = f_0 * gp_25[k];

        t_17[k] = f_0 * gp_26[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, dp_9, dp_10, dp_11, dp_12, dp_13, \
                         gp_30, gp_31, gp_32, gp_33, gp_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = -3.0 * dp_9[k]
                  + f_0 * gp_30[k];

        t_19[k] = -3.0 * dp_10[k]
                  + f_0 * gp_31[k];

        t_20[k] = -3.0 * dp_11[k]
                  + f_0 * gp_32[k];

        t_21[k] = -2.0 * dp_12[k]
                  + f_0 * gp_33[k];

        t_22[k] = -2.0 * dp_13[k]
                  + f_0 * gp_34[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, t_27, t_28, dp_14, dp_15, dp_16, dp_17, \
                         gp_35, gp_36, gp_37, gp_38, gp_39, gp_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = -2.0 * dp_14[k]
                  + f_0 * gp_35[k];

        t_24[k] = -dp_15[k]
                  + f_0 * gp_36[k];

        t_25[k] = -dp_16[k]
                  + f_0 * gp_37[k];

        t_26[k] = -dp_17[k]
                  + f_0 * gp_38[k];

        t_27[k] = f_0 * gp_39[k];

        t_28[k] = f_0 * gp_40[k];
    }

#pragma omp simd aligned(t_29, gp_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_0 * gp_41[k];
    }
}

auto
compute_prim_geom_10_fp_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t dp, const size_t gp,
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
    auto *t_28 = buffer.data(target + 28);
    auto *t_29 = buffer.data(target + 29);

    const auto *dp_0 = buffer.data(dp + 0);
    const auto *dp_1 = buffer.data(dp + 1);
    const auto *dp_2 = buffer.data(dp + 2);
    const auto *dp_3 = buffer.data(dp + 3);
    const auto *dp_4 = buffer.data(dp + 4);
    const auto *dp_5 = buffer.data(dp + 5);
    const auto *dp_6 = buffer.data(dp + 6);
    const auto *dp_7 = buffer.data(dp + 7);
    const auto *dp_8 = buffer.data(dp + 8);
    const auto *dp_9 = buffer.data(dp + 9);
    const auto *dp_10 = buffer.data(dp + 10);
    const auto *dp_11 = buffer.data(dp + 11);
    const auto *dp_12 = buffer.data(dp + 12);
    const auto *dp_13 = buffer.data(dp + 13);
    const auto *dp_14 = buffer.data(dp + 14);
    const auto *dp_15 = buffer.data(dp + 15);
    const auto *dp_16 = buffer.data(dp + 16);
    const auto *dp_17 = buffer.data(dp + 17);

    const auto *gp_6 = buffer.data(gp + 6);
    const auto *gp_7 = buffer.data(gp + 7);
    const auto *gp_8 = buffer.data(gp + 8);
    const auto *gp_12 = buffer.data(gp + 12);
    const auto *gp_13 = buffer.data(gp + 13);
    const auto *gp_14 = buffer.data(gp + 14);
    const auto *gp_15 = buffer.data(gp + 15);
    const auto *gp_16 = buffer.data(gp + 16);
    const auto *gp_17 = buffer.data(gp + 17);
    const auto *gp_21 = buffer.data(gp + 21);
    const auto *gp_22 = buffer.data(gp + 22);
    const auto *gp_23 = buffer.data(gp + 23);
    const auto *gp_24 = buffer.data(gp + 24);
    const auto *gp_25 = buffer.data(gp + 25);
    const auto *gp_26 = buffer.data(gp + 26);
    const auto *gp_27 = buffer.data(gp + 27);
    const auto *gp_28 = buffer.data(gp + 28);
    const auto *gp_29 = buffer.data(gp + 29);
    const auto *gp_33 = buffer.data(gp + 33);
    const auto *gp_34 = buffer.data(gp + 34);
    const auto *gp_35 = buffer.data(gp + 35);
    const auto *gp_36 = buffer.data(gp + 36);
    const auto *gp_37 = buffer.data(gp + 37);
    const auto *gp_38 = buffer.data(gp + 38);
    const auto *gp_39 = buffer.data(gp + 39);
    const auto *gp_40 = buffer.data(gp + 40);
    const auto *gp_41 = buffer.data(gp + 41);
    const auto *gp_42 = buffer.data(gp + 42);
    const auto *gp_43 = buffer.data(gp + 43);
    const auto *gp_44 = buffer.data(gp + 44);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, dp_0, gp_6, gp_7, gp_8, gp_12, \
                         gp_13, gp_14, gp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gp_6[k];

        t_1[k] = f_0 * gp_7[k];

        t_2[k] = f_0 * gp_8[k];

        t_3[k] = f_0 * gp_12[k];

        t_4[k] = f_0 * gp_13[k];

        t_5[k] = f_0 * gp_14[k];

        t_6[k] = -dp_0[k]
                 + f_0 * gp_15[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, t_11, t_12, dp_1, dp_2, dp_3, gp_16, gp_17, \
                         gp_21, gp_22, gp_23, gp_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = -dp_1[k]
                 + f_0 * gp_16[k];

        t_8[k] = -dp_2[k]
                 + f_0 * gp_17[k];

        t_9[k] = f_0 * gp_21[k];

        t_10[k] = f_0 * gp_22[k];

        t_11[k] = f_0 * gp_23[k];

        t_12[k] = -dp_3[k]
                  + f_0 * gp_24[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, dp_4, dp_5, dp_6, dp_7, dp_8, gp_25, \
                         gp_26, gp_27, gp_28, gp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = -dp_4[k]
                  + f_0 * gp_25[k];

        t_14[k] = -dp_5[k]
                  + f_0 * gp_26[k];

        t_15[k] = -2.0 * dp_6[k]
                  + f_0 * gp_27[k];

        t_16[k] = -2.0 * dp_7[k]
                  + f_0 * gp_28[k];

        t_17[k] = -2.0 * dp_8[k]
                  + f_0 * gp_29[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, t_23, dp_9, dp_10, dp_11, gp_33, gp_34, \
                         gp_35, gp_36, gp_37, gp_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_0 * gp_33[k];

        t_19[k] = f_0 * gp_34[k];

        t_20[k] = f_0 * gp_35[k];

        t_21[k] = -dp_9[k]
                  + f_0 * gp_36[k];

        t_22[k] = -dp_10[k]
                  + f_0 * gp_37[k];

        t_23[k] = -dp_11[k]
                  + f_0 * gp_38[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, dp_12, dp_13, dp_14, dp_15, dp_16, \
                         gp_39, gp_40, gp_41, gp_42, gp_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = -2.0 * dp_12[k]
                  + f_0 * gp_39[k];

        t_25[k] = -2.0 * dp_13[k]
                  + f_0 * gp_40[k];

        t_26[k] = -2.0 * dp_14[k]
                  + f_0 * gp_41[k];

        t_27[k] = -3.0 * dp_15[k]
                  + f_0 * gp_42[k];

        t_28[k] = -3.0 * dp_16[k]
                  + f_0 * gp_43[k];
    }

#pragma omp simd aligned(t_29, dp_17, gp_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = -3.0 * dp_17[k]
                  + f_0 * gp_44[k];
    }
}

}  // namespace simdt2ceri
