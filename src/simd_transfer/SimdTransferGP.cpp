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


#include "SimdTransferGP.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_gp_out_of_first(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                            const size_t target, const size_t gs, const size_t hs,
                            const size_t ncomps, const size_t nmax) -> void
{
    // NOTE: what the block carries below the pair reaches this step as a count.
    // The coefficients are powers of the separation and name no index of it, so
    // it is stepped over rather than known.

    for (size_t c = 0; c < ncomps; c++)
    {
        auto *t_0 = buffer.data(target + 0 * ncomps + c);
        auto *t_1 = buffer.data(target + 1 * ncomps + c);
        auto *t_2 = buffer.data(target + 2 * ncomps + c);
        auto *t_3 = buffer.data(target + 3 * ncomps + c);
        auto *t_4 = buffer.data(target + 4 * ncomps + c);
        auto *t_5 = buffer.data(target + 5 * ncomps + c);
        auto *t_6 = buffer.data(target + 6 * ncomps + c);
        auto *t_7 = buffer.data(target + 7 * ncomps + c);
        auto *t_8 = buffer.data(target + 8 * ncomps + c);
        auto *t_9 = buffer.data(target + 9 * ncomps + c);
        auto *t_10 = buffer.data(target + 10 * ncomps + c);
        auto *t_11 = buffer.data(target + 11 * ncomps + c);
        auto *t_12 = buffer.data(target + 12 * ncomps + c);
        auto *t_13 = buffer.data(target + 13 * ncomps + c);
        auto *t_14 = buffer.data(target + 14 * ncomps + c);
        auto *t_15 = buffer.data(target + 15 * ncomps + c);
        auto *t_16 = buffer.data(target + 16 * ncomps + c);
        auto *t_17 = buffer.data(target + 17 * ncomps + c);
        auto *t_18 = buffer.data(target + 18 * ncomps + c);
        auto *t_19 = buffer.data(target + 19 * ncomps + c);
        auto *t_20 = buffer.data(target + 20 * ncomps + c);
        auto *t_21 = buffer.data(target + 21 * ncomps + c);
        auto *t_22 = buffer.data(target + 22 * ncomps + c);
        auto *t_23 = buffer.data(target + 23 * ncomps + c);
        auto *t_24 = buffer.data(target + 24 * ncomps + c);
        auto *t_25 = buffer.data(target + 25 * ncomps + c);
        auto *t_26 = buffer.data(target + 26 * ncomps + c);
        auto *t_27 = buffer.data(target + 27 * ncomps + c);
        auto *t_28 = buffer.data(target + 28 * ncomps + c);
        auto *t_29 = buffer.data(target + 29 * ncomps + c);
        auto *t_30 = buffer.data(target + 30 * ncomps + c);
        auto *t_31 = buffer.data(target + 31 * ncomps + c);
        auto *t_32 = buffer.data(target + 32 * ncomps + c);
        auto *t_33 = buffer.data(target + 33 * ncomps + c);
        auto *t_34 = buffer.data(target + 34 * ncomps + c);
        auto *t_35 = buffer.data(target + 35 * ncomps + c);
        auto *t_36 = buffer.data(target + 36 * ncomps + c);
        auto *t_37 = buffer.data(target + 37 * ncomps + c);
        auto *t_38 = buffer.data(target + 38 * ncomps + c);
        auto *t_39 = buffer.data(target + 39 * ncomps + c);
        auto *t_40 = buffer.data(target + 40 * ncomps + c);
        auto *t_41 = buffer.data(target + 41 * ncomps + c);
        auto *t_42 = buffer.data(target + 42 * ncomps + c);
        auto *t_43 = buffer.data(target + 43 * ncomps + c);
        auto *t_44 = buffer.data(target + 44 * ncomps + c);

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *gs_0 = buffer.data(gs + 0 * ncomps + c);
        const auto *gs_1 = buffer.data(gs + 1 * ncomps + c);
        const auto *gs_2 = buffer.data(gs + 2 * ncomps + c);
        const auto *gs_3 = buffer.data(gs + 3 * ncomps + c);
        const auto *gs_4 = buffer.data(gs + 4 * ncomps + c);
        const auto *gs_5 = buffer.data(gs + 5 * ncomps + c);
        const auto *gs_6 = buffer.data(gs + 6 * ncomps + c);
        const auto *gs_7 = buffer.data(gs + 7 * ncomps + c);
        const auto *gs_8 = buffer.data(gs + 8 * ncomps + c);
        const auto *gs_9 = buffer.data(gs + 9 * ncomps + c);
        const auto *gs_10 = buffer.data(gs + 10 * ncomps + c);
        const auto *gs_11 = buffer.data(gs + 11 * ncomps + c);
        const auto *gs_12 = buffer.data(gs + 12 * ncomps + c);
        const auto *gs_13 = buffer.data(gs + 13 * ncomps + c);
        const auto *gs_14 = buffer.data(gs + 14 * ncomps + c);

        const auto *hs_0 = buffer.data(hs + 0 * ncomps + c);
        const auto *hs_1 = buffer.data(hs + 1 * ncomps + c);
        const auto *hs_2 = buffer.data(hs + 2 * ncomps + c);
        const auto *hs_3 = buffer.data(hs + 3 * ncomps + c);
        const auto *hs_4 = buffer.data(hs + 4 * ncomps + c);
        const auto *hs_5 = buffer.data(hs + 5 * ncomps + c);
        const auto *hs_6 = buffer.data(hs + 6 * ncomps + c);
        const auto *hs_7 = buffer.data(hs + 7 * ncomps + c);
        const auto *hs_8 = buffer.data(hs + 8 * ncomps + c);
        const auto *hs_9 = buffer.data(hs + 9 * ncomps + c);
        const auto *hs_10 = buffer.data(hs + 10 * ncomps + c);
        const auto *hs_11 = buffer.data(hs + 11 * ncomps + c);
        const auto *hs_12 = buffer.data(hs + 12 * ncomps + c);
        const auto *hs_13 = buffer.data(hs + 13 * ncomps + c);
        const auto *hs_14 = buffer.data(hs + 14 * ncomps + c);
        const auto *hs_15 = buffer.data(hs + 15 * ncomps + c);
        const auto *hs_16 = buffer.data(hs + 16 * ncomps + c);
        const auto *hs_17 = buffer.data(hs + 17 * ncomps + c);
        const auto *hs_18 = buffer.data(hs + 18 * ncomps + c);
        const auto *hs_19 = buffer.data(hs + 19 * ncomps + c);
        const auto *hs_20 = buffer.data(hs + 20 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, ab_x, ab_y, ab_z, gs_0, gs_1, hs_0, \
                         hs_1, hs_2, hs_3, hs_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * gs_0[k]
                     + hs_0[k];

            t_1[k] = ab_y[k] * gs_0[k]
                     + hs_1[k];

            t_2[k] = ab_z[k] * gs_0[k]
                     + hs_2[k];

            t_3[k] = ab_x[k] * gs_1[k]
                     + hs_1[k];

            t_4[k] = ab_y[k] * gs_1[k]
                     + hs_3[k];

            t_5[k] = ab_z[k] * gs_1[k]
                     + hs_4[k];
        }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, ab_x, ab_y, ab_z, gs_2, gs_3, hs_2, hs_3, \
                         hs_4, hs_5, hs_6 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_6[k] = ab_x[k] * gs_2[k]
                     + hs_2[k];

            t_7[k] = ab_y[k] * gs_2[k]
                     + hs_4[k];

            t_8[k] = ab_z[k] * gs_2[k]
                     + hs_5[k];

            t_9[k] = ab_x[k] * gs_3[k]
                     + hs_3[k];

            t_10[k] = ab_y[k] * gs_3[k]
                      + hs_6[k];
        }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, ab_x, ab_y, ab_z, gs_3, gs_4, \
                         gs_5, hs_4, hs_5, hs_7, hs_8 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_11[k] = ab_z[k] * gs_3[k]
                      + hs_7[k];

            t_12[k] = ab_x[k] * gs_4[k]
                      + hs_4[k];

            t_13[k] = ab_y[k] * gs_4[k]
                      + hs_7[k];

            t_14[k] = ab_z[k] * gs_4[k]
                      + hs_8[k];

            t_15[k] = ab_x[k] * gs_5[k]
                      + hs_5[k];

            t_16[k] = ab_y[k] * gs_5[k]
                      + hs_8[k];
        }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, t_21, ab_x, ab_y, ab_z, gs_5, gs_6, gs_7, \
                         hs_6, hs_7, hs_9, hs_10, hs_11 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_17[k] = ab_z[k] * gs_5[k]
                      + hs_9[k];

            t_18[k] = ab_x[k] * gs_6[k]
                      + hs_6[k];

            t_19[k] = ab_y[k] * gs_6[k]
                      + hs_10[k];

            t_20[k] = ab_z[k] * gs_6[k]
                      + hs_11[k];

            t_21[k] = ab_x[k] * gs_7[k]
                      + hs_7[k];
        }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, ab_x, ab_y, ab_z, gs_7, gs_8, hs_8, \
                         hs_11, hs_12, hs_13 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_22[k] = ab_y[k] * gs_7[k]
                      + hs_11[k];

            t_23[k] = ab_z[k] * gs_7[k]
                      + hs_12[k];

            t_24[k] = ab_x[k] * gs_8[k]
                      + hs_8[k];

            t_25[k] = ab_y[k] * gs_8[k]
                      + hs_12[k];

            t_26[k] = ab_z[k] * gs_8[k]
                      + hs_13[k];
        }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, ab_x, ab_y, ab_z, gs_9, gs_10, hs_9, \
                         hs_10, hs_13, hs_14, hs_15 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_27[k] = ab_x[k] * gs_9[k]
                      + hs_9[k];

            t_28[k] = ab_y[k] * gs_9[k]
                      + hs_13[k];

            t_29[k] = ab_z[k] * gs_9[k]
                      + hs_14[k];

            t_30[k] = ab_x[k] * gs_10[k]
                      + hs_10[k];

            t_31[k] = ab_y[k] * gs_10[k]
                      + hs_15[k];
        }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, ab_x, ab_y, ab_z, gs_10, gs_11, \
                         gs_12, hs_11, hs_12, hs_16, hs_17 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_32[k] = ab_z[k] * gs_10[k]
                      + hs_16[k];

            t_33[k] = ab_x[k] * gs_11[k]
                      + hs_11[k];

            t_34[k] = ab_y[k] * gs_11[k]
                      + hs_16[k];

            t_35[k] = ab_z[k] * gs_11[k]
                      + hs_17[k];

            t_36[k] = ab_x[k] * gs_12[k]
                      + hs_12[k];

            t_37[k] = ab_y[k] * gs_12[k]
                      + hs_17[k];
        }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, t_43, ab_x, ab_y, ab_z, gs_12, gs_13, \
                         gs_14, hs_13, hs_14, hs_18, hs_19 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_38[k] = ab_z[k] * gs_12[k]
                      + hs_18[k];

            t_39[k] = ab_x[k] * gs_13[k]
                      + hs_13[k];

            t_40[k] = ab_y[k] * gs_13[k]
                      + hs_18[k];

            t_41[k] = ab_z[k] * gs_13[k]
                      + hs_19[k];

            t_42[k] = ab_x[k] * gs_14[k]
                      + hs_14[k];

            t_43[k] = ab_y[k] * gs_14[k]
                      + hs_19[k];
        }

#pragma omp simd aligned(t_44, ab_z, gs_14, hs_20 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_44[k] = ab_z[k] * gs_14[k]
                      + hs_20[k];
        }
    }
}

auto
compute_hrr_gp(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t gs, const size_t hs, const size_t ncomps, const size_t nmax) -> void
{
    // NOTE: what the block carries below the pair reaches this step as a count.
    // The coefficients are powers of the separation and name no index of it, so
    // it is stepped over rather than known.

    for (size_t c = 0; c < ncomps; c++)
    {
        auto *t_0 = buffer.data(target + 0 * ncomps + c);
        auto *t_1 = buffer.data(target + 1 * ncomps + c);
        auto *t_2 = buffer.data(target + 2 * ncomps + c);
        auto *t_3 = buffer.data(target + 3 * ncomps + c);
        auto *t_4 = buffer.data(target + 4 * ncomps + c);
        auto *t_5 = buffer.data(target + 5 * ncomps + c);
        auto *t_6 = buffer.data(target + 6 * ncomps + c);
        auto *t_7 = buffer.data(target + 7 * ncomps + c);
        auto *t_8 = buffer.data(target + 8 * ncomps + c);
        auto *t_9 = buffer.data(target + 9 * ncomps + c);
        auto *t_10 = buffer.data(target + 10 * ncomps + c);
        auto *t_11 = buffer.data(target + 11 * ncomps + c);
        auto *t_12 = buffer.data(target + 12 * ncomps + c);
        auto *t_13 = buffer.data(target + 13 * ncomps + c);
        auto *t_14 = buffer.data(target + 14 * ncomps + c);
        auto *t_15 = buffer.data(target + 15 * ncomps + c);
        auto *t_16 = buffer.data(target + 16 * ncomps + c);
        auto *t_17 = buffer.data(target + 17 * ncomps + c);
        auto *t_18 = buffer.data(target + 18 * ncomps + c);
        auto *t_19 = buffer.data(target + 19 * ncomps + c);
        auto *t_20 = buffer.data(target + 20 * ncomps + c);
        auto *t_21 = buffer.data(target + 21 * ncomps + c);
        auto *t_22 = buffer.data(target + 22 * ncomps + c);
        auto *t_23 = buffer.data(target + 23 * ncomps + c);
        auto *t_24 = buffer.data(target + 24 * ncomps + c);
        auto *t_25 = buffer.data(target + 25 * ncomps + c);
        auto *t_26 = buffer.data(target + 26 * ncomps + c);
        auto *t_27 = buffer.data(target + 27 * ncomps + c);
        auto *t_28 = buffer.data(target + 28 * ncomps + c);
        auto *t_29 = buffer.data(target + 29 * ncomps + c);
        auto *t_30 = buffer.data(target + 30 * ncomps + c);
        auto *t_31 = buffer.data(target + 31 * ncomps + c);
        auto *t_32 = buffer.data(target + 32 * ncomps + c);
        auto *t_33 = buffer.data(target + 33 * ncomps + c);
        auto *t_34 = buffer.data(target + 34 * ncomps + c);
        auto *t_35 = buffer.data(target + 35 * ncomps + c);
        auto *t_36 = buffer.data(target + 36 * ncomps + c);
        auto *t_37 = buffer.data(target + 37 * ncomps + c);
        auto *t_38 = buffer.data(target + 38 * ncomps + c);
        auto *t_39 = buffer.data(target + 39 * ncomps + c);
        auto *t_40 = buffer.data(target + 40 * ncomps + c);
        auto *t_41 = buffer.data(target + 41 * ncomps + c);
        auto *t_42 = buffer.data(target + 42 * ncomps + c);
        auto *t_43 = buffer.data(target + 43 * ncomps + c);
        auto *t_44 = buffer.data(target + 44 * ncomps + c);

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *gs_0 = buffer.data(gs + 0 * ncomps + c);
        const auto *gs_1 = buffer.data(gs + 1 * ncomps + c);
        const auto *gs_2 = buffer.data(gs + 2 * ncomps + c);
        const auto *gs_3 = buffer.data(gs + 3 * ncomps + c);
        const auto *gs_4 = buffer.data(gs + 4 * ncomps + c);
        const auto *gs_5 = buffer.data(gs + 5 * ncomps + c);
        const auto *gs_6 = buffer.data(gs + 6 * ncomps + c);
        const auto *gs_7 = buffer.data(gs + 7 * ncomps + c);
        const auto *gs_8 = buffer.data(gs + 8 * ncomps + c);
        const auto *gs_9 = buffer.data(gs + 9 * ncomps + c);
        const auto *gs_10 = buffer.data(gs + 10 * ncomps + c);
        const auto *gs_11 = buffer.data(gs + 11 * ncomps + c);
        const auto *gs_12 = buffer.data(gs + 12 * ncomps + c);
        const auto *gs_13 = buffer.data(gs + 13 * ncomps + c);
        const auto *gs_14 = buffer.data(gs + 14 * ncomps + c);

        const auto *hs_0 = buffer.data(hs + 0 * ncomps + c);
        const auto *hs_1 = buffer.data(hs + 1 * ncomps + c);
        const auto *hs_2 = buffer.data(hs + 2 * ncomps + c);
        const auto *hs_3 = buffer.data(hs + 3 * ncomps + c);
        const auto *hs_4 = buffer.data(hs + 4 * ncomps + c);
        const auto *hs_5 = buffer.data(hs + 5 * ncomps + c);
        const auto *hs_6 = buffer.data(hs + 6 * ncomps + c);
        const auto *hs_7 = buffer.data(hs + 7 * ncomps + c);
        const auto *hs_8 = buffer.data(hs + 8 * ncomps + c);
        const auto *hs_9 = buffer.data(hs + 9 * ncomps + c);
        const auto *hs_10 = buffer.data(hs + 10 * ncomps + c);
        const auto *hs_11 = buffer.data(hs + 11 * ncomps + c);
        const auto *hs_12 = buffer.data(hs + 12 * ncomps + c);
        const auto *hs_13 = buffer.data(hs + 13 * ncomps + c);
        const auto *hs_14 = buffer.data(hs + 14 * ncomps + c);
        const auto *hs_15 = buffer.data(hs + 15 * ncomps + c);
        const auto *hs_16 = buffer.data(hs + 16 * ncomps + c);
        const auto *hs_17 = buffer.data(hs + 17 * ncomps + c);
        const auto *hs_18 = buffer.data(hs + 18 * ncomps + c);
        const auto *hs_19 = buffer.data(hs + 19 * ncomps + c);
        const auto *hs_20 = buffer.data(hs + 20 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, ab_x, ab_y, ab_z, gs_0, gs_1, hs_0, \
                         hs_1, hs_2, hs_3, hs_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * gs_0[k]
                     + hs_0[k];

            t_1[k] = ab_y[k] * gs_0[k]
                     + hs_1[k];

            t_2[k] = ab_z[k] * gs_0[k]
                     + hs_2[k];

            t_3[k] = ab_x[k] * gs_1[k]
                     + hs_1[k];

            t_4[k] = ab_y[k] * gs_1[k]
                     + hs_3[k];

            t_5[k] = ab_z[k] * gs_1[k]
                     + hs_4[k];
        }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, ab_x, ab_y, ab_z, gs_2, gs_3, hs_2, hs_3, \
                         hs_4, hs_5, hs_6 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_6[k] = ab_x[k] * gs_2[k]
                     + hs_2[k];

            t_7[k] = ab_y[k] * gs_2[k]
                     + hs_4[k];

            t_8[k] = ab_z[k] * gs_2[k]
                     + hs_5[k];

            t_9[k] = ab_x[k] * gs_3[k]
                     + hs_3[k];

            t_10[k] = ab_y[k] * gs_3[k]
                      + hs_6[k];
        }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, ab_x, ab_y, ab_z, gs_3, gs_4, \
                         gs_5, hs_4, hs_5, hs_7, hs_8 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_11[k] = ab_z[k] * gs_3[k]
                      + hs_7[k];

            t_12[k] = ab_x[k] * gs_4[k]
                      + hs_4[k];

            t_13[k] = ab_y[k] * gs_4[k]
                      + hs_7[k];

            t_14[k] = ab_z[k] * gs_4[k]
                      + hs_8[k];

            t_15[k] = ab_x[k] * gs_5[k]
                      + hs_5[k];

            t_16[k] = ab_y[k] * gs_5[k]
                      + hs_8[k];
        }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, t_21, ab_x, ab_y, ab_z, gs_5, gs_6, gs_7, \
                         hs_6, hs_7, hs_9, hs_10, hs_11 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_17[k] = ab_z[k] * gs_5[k]
                      + hs_9[k];

            t_18[k] = ab_x[k] * gs_6[k]
                      + hs_6[k];

            t_19[k] = ab_y[k] * gs_6[k]
                      + hs_10[k];

            t_20[k] = ab_z[k] * gs_6[k]
                      + hs_11[k];

            t_21[k] = ab_x[k] * gs_7[k]
                      + hs_7[k];
        }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, ab_x, ab_y, ab_z, gs_7, gs_8, hs_8, \
                         hs_11, hs_12, hs_13 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_22[k] = ab_y[k] * gs_7[k]
                      + hs_11[k];

            t_23[k] = ab_z[k] * gs_7[k]
                      + hs_12[k];

            t_24[k] = ab_x[k] * gs_8[k]
                      + hs_8[k];

            t_25[k] = ab_y[k] * gs_8[k]
                      + hs_12[k];

            t_26[k] = ab_z[k] * gs_8[k]
                      + hs_13[k];
        }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, ab_x, ab_y, ab_z, gs_9, gs_10, hs_9, \
                         hs_10, hs_13, hs_14, hs_15 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_27[k] = ab_x[k] * gs_9[k]
                      + hs_9[k];

            t_28[k] = ab_y[k] * gs_9[k]
                      + hs_13[k];

            t_29[k] = ab_z[k] * gs_9[k]
                      + hs_14[k];

            t_30[k] = ab_x[k] * gs_10[k]
                      + hs_10[k];

            t_31[k] = ab_y[k] * gs_10[k]
                      + hs_15[k];
        }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, ab_x, ab_y, ab_z, gs_10, gs_11, \
                         gs_12, hs_11, hs_12, hs_16, hs_17 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_32[k] = ab_z[k] * gs_10[k]
                      + hs_16[k];

            t_33[k] = ab_x[k] * gs_11[k]
                      + hs_11[k];

            t_34[k] = ab_y[k] * gs_11[k]
                      + hs_16[k];

            t_35[k] = ab_z[k] * gs_11[k]
                      + hs_17[k];

            t_36[k] = ab_x[k] * gs_12[k]
                      + hs_12[k];

            t_37[k] = ab_y[k] * gs_12[k]
                      + hs_17[k];
        }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, t_43, ab_x, ab_y, ab_z, gs_12, gs_13, \
                         gs_14, hs_13, hs_14, hs_18, hs_19 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_38[k] = ab_z[k] * gs_12[k]
                      + hs_18[k];

            t_39[k] = ab_x[k] * gs_13[k]
                      + hs_13[k];

            t_40[k] = ab_y[k] * gs_13[k]
                      + hs_18[k];

            t_41[k] = ab_z[k] * gs_13[k]
                      + hs_19[k];

            t_42[k] = ab_x[k] * gs_14[k]
                      + hs_14[k];

            t_43[k] = ab_y[k] * gs_14[k]
                      + hs_19[k];
        }

#pragma omp simd aligned(t_44, ab_z, gs_14, hs_20 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_44[k] = ab_z[k] * gs_14[k]
                      + hs_20[k];
        }
    }
}

}  // namespace simdtrf
