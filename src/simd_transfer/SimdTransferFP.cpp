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


#include "SimdTransferFP.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_fp_out_of_first(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                            const size_t target, const size_t fs, const size_t gs,
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

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *fs_0 = buffer.data(fs + 0 * ncomps + c);
        const auto *fs_1 = buffer.data(fs + 1 * ncomps + c);
        const auto *fs_2 = buffer.data(fs + 2 * ncomps + c);
        const auto *fs_3 = buffer.data(fs + 3 * ncomps + c);
        const auto *fs_4 = buffer.data(fs + 4 * ncomps + c);
        const auto *fs_5 = buffer.data(fs + 5 * ncomps + c);
        const auto *fs_6 = buffer.data(fs + 6 * ncomps + c);
        const auto *fs_7 = buffer.data(fs + 7 * ncomps + c);
        const auto *fs_8 = buffer.data(fs + 8 * ncomps + c);
        const auto *fs_9 = buffer.data(fs + 9 * ncomps + c);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, ab_x, ab_y, ab_z, fs_0, fs_1, gs_0, \
                         gs_1, gs_2, gs_3, gs_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * fs_0[k]
                     + gs_0[k];

            t_1[k] = ab_y[k] * fs_0[k]
                     + gs_1[k];

            t_2[k] = ab_z[k] * fs_0[k]
                     + gs_2[k];

            t_3[k] = ab_x[k] * fs_1[k]
                     + gs_1[k];

            t_4[k] = ab_y[k] * fs_1[k]
                     + gs_3[k];

            t_5[k] = ab_z[k] * fs_1[k]
                     + gs_4[k];
        }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, ab_x, ab_y, ab_z, fs_2, fs_3, gs_2, gs_3, \
                         gs_4, gs_5, gs_6 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_6[k] = ab_x[k] * fs_2[k]
                     + gs_2[k];

            t_7[k] = ab_y[k] * fs_2[k]
                     + gs_4[k];

            t_8[k] = ab_z[k] * fs_2[k]
                     + gs_5[k];

            t_9[k] = ab_x[k] * fs_3[k]
                     + gs_3[k];

            t_10[k] = ab_y[k] * fs_3[k]
                      + gs_6[k];
        }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, ab_x, ab_y, ab_z, fs_3, fs_4, \
                         fs_5, gs_4, gs_5, gs_7, gs_8 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_11[k] = ab_z[k] * fs_3[k]
                      + gs_7[k];

            t_12[k] = ab_x[k] * fs_4[k]
                      + gs_4[k];

            t_13[k] = ab_y[k] * fs_4[k]
                      + gs_7[k];

            t_14[k] = ab_z[k] * fs_4[k]
                      + gs_8[k];

            t_15[k] = ab_x[k] * fs_5[k]
                      + gs_5[k];

            t_16[k] = ab_y[k] * fs_5[k]
                      + gs_8[k];
        }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, t_21, ab_x, ab_y, ab_z, fs_5, fs_6, fs_7, \
                         gs_6, gs_7, gs_9, gs_10, gs_11 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_17[k] = ab_z[k] * fs_5[k]
                      + gs_9[k];

            t_18[k] = ab_x[k] * fs_6[k]
                      + gs_6[k];

            t_19[k] = ab_y[k] * fs_6[k]
                      + gs_10[k];

            t_20[k] = ab_z[k] * fs_6[k]
                      + gs_11[k];

            t_21[k] = ab_x[k] * fs_7[k]
                      + gs_7[k];
        }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, ab_x, ab_y, ab_z, fs_7, fs_8, gs_8, \
                         gs_11, gs_12, gs_13 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_22[k] = ab_y[k] * fs_7[k]
                      + gs_11[k];

            t_23[k] = ab_z[k] * fs_7[k]
                      + gs_12[k];

            t_24[k] = ab_x[k] * fs_8[k]
                      + gs_8[k];

            t_25[k] = ab_y[k] * fs_8[k]
                      + gs_12[k];

            t_26[k] = ab_z[k] * fs_8[k]
                      + gs_13[k];
        }

#pragma omp simd aligned(t_27, t_28, t_29, ab_x, ab_y, ab_z, fs_9, gs_9, gs_13, \
                         gs_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_27[k] = ab_x[k] * fs_9[k]
                      + gs_9[k];

            t_28[k] = ab_y[k] * fs_9[k]
                      + gs_13[k];

            t_29[k] = ab_z[k] * fs_9[k]
                      + gs_14[k];
        }
    }
}

auto
compute_hrr_fp_out_of_second(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                             const size_t target, const size_t dp, const size_t dd,
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

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *dp_0 = buffer.data(dp + 0 * ncomps + c);
        const auto *dp_1 = buffer.data(dp + 1 * ncomps + c);
        const auto *dp_2 = buffer.data(dp + 2 * ncomps + c);
        const auto *dp_3 = buffer.data(dp + 3 * ncomps + c);
        const auto *dp_4 = buffer.data(dp + 4 * ncomps + c);
        const auto *dp_5 = buffer.data(dp + 5 * ncomps + c);
        const auto *dp_6 = buffer.data(dp + 6 * ncomps + c);
        const auto *dp_7 = buffer.data(dp + 7 * ncomps + c);
        const auto *dp_8 = buffer.data(dp + 8 * ncomps + c);
        const auto *dp_9 = buffer.data(dp + 9 * ncomps + c);
        const auto *dp_10 = buffer.data(dp + 10 * ncomps + c);
        const auto *dp_11 = buffer.data(dp + 11 * ncomps + c);
        const auto *dp_12 = buffer.data(dp + 12 * ncomps + c);
        const auto *dp_13 = buffer.data(dp + 13 * ncomps + c);
        const auto *dp_14 = buffer.data(dp + 14 * ncomps + c);
        const auto *dp_15 = buffer.data(dp + 15 * ncomps + c);
        const auto *dp_16 = buffer.data(dp + 16 * ncomps + c);
        const auto *dp_17 = buffer.data(dp + 17 * ncomps + c);

        const auto *dd_0 = buffer.data(dd + 0 * ncomps + c);
        const auto *dd_1 = buffer.data(dd + 1 * ncomps + c);
        const auto *dd_2 = buffer.data(dd + 2 * ncomps + c);
        const auto *dd_6 = buffer.data(dd + 6 * ncomps + c);
        const auto *dd_7 = buffer.data(dd + 7 * ncomps + c);
        const auto *dd_8 = buffer.data(dd + 8 * ncomps + c);
        const auto *dd_12 = buffer.data(dd + 12 * ncomps + c);
        const auto *dd_13 = buffer.data(dd + 13 * ncomps + c);
        const auto *dd_14 = buffer.data(dd + 14 * ncomps + c);
        const auto *dd_18 = buffer.data(dd + 18 * ncomps + c);
        const auto *dd_19 = buffer.data(dd + 19 * ncomps + c);
        const auto *dd_20 = buffer.data(dd + 20 * ncomps + c);
        const auto *dd_21 = buffer.data(dd + 21 * ncomps + c);
        const auto *dd_22 = buffer.data(dd + 22 * ncomps + c);
        const auto *dd_24 = buffer.data(dd + 24 * ncomps + c);
        const auto *dd_25 = buffer.data(dd + 25 * ncomps + c);
        const auto *dd_26 = buffer.data(dd + 26 * ncomps + c);
        const auto *dd_27 = buffer.data(dd + 27 * ncomps + c);
        const auto *dd_28 = buffer.data(dd + 28 * ncomps + c);
        const auto *dd_30 = buffer.data(dd + 30 * ncomps + c);
        const auto *dd_31 = buffer.data(dd + 31 * ncomps + c);
        const auto *dd_32 = buffer.data(dd + 32 * ncomps + c);
        const auto *dd_33 = buffer.data(dd + 33 * ncomps + c);
        const auto *dd_34 = buffer.data(dd + 34 * ncomps + c);
        const auto *dd_35 = buffer.data(dd + 35 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, dp_0, dp_1, dp_2, dp_3, dp_4, dd_0, \
                         dd_1, dd_2, dd_6, dd_7 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = -ab_x[k] * dp_0[k]
                     + dd_0[k];

            t_1[k] = -ab_x[k] * dp_1[k]
                     + dd_1[k];

            t_2[k] = -ab_x[k] * dp_2[k]
                     + dd_2[k];

            t_3[k] = -ab_x[k] * dp_3[k]
                     + dd_6[k];

            t_4[k] = -ab_x[k] * dp_4[k]
                     + dd_7[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, dp_5, dp_6, dp_7, dp_8, dp_9, dd_8, \
                         dd_12, dd_13, dd_14, dd_18 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = -ab_x[k] * dp_5[k]
                     + dd_8[k];

            t_6[k] = -ab_x[k] * dp_6[k]
                     + dd_12[k];

            t_7[k] = -ab_x[k] * dp_7[k]
                     + dd_13[k];

            t_8[k] = -ab_x[k] * dp_8[k]
                     + dd_14[k];

            t_9[k] = -ab_x[k] * dp_9[k]
                     + dd_18[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, dp_10, dp_11, dp_12, dp_13, \
                         dp_14, dd_19, dd_20, dd_24, dd_25, dd_26 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = -ab_x[k] * dp_10[k]
                      + dd_19[k];

            t_11[k] = -ab_x[k] * dp_11[k]
                      + dd_20[k];

            t_12[k] = -ab_x[k] * dp_12[k]
                      + dd_24[k];

            t_13[k] = -ab_x[k] * dp_13[k]
                      + dd_25[k];

            t_14[k] = -ab_x[k] * dp_14[k]
                      + dd_26[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, ab_x, ab_y, dp_9, dp_15, dp_16, dp_17, dd_19, \
                         dd_30, dd_31, dd_32 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = -ab_x[k] * dp_15[k]
                      + dd_30[k];

            t_16[k] = -ab_x[k] * dp_16[k]
                      + dd_31[k];

            t_17[k] = -ab_x[k] * dp_17[k]
                      + dd_32[k];

            t_18[k] = -ab_y[k] * dp_9[k]
                      + dd_19[k];
        }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, ab_y, dp_10, dp_11, dp_12, dp_13, \
                         dp_14, dd_21, dd_22, dd_25, dd_27, dd_28 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_19[k] = -ab_y[k] * dp_10[k]
                      + dd_21[k];

            t_20[k] = -ab_y[k] * dp_11[k]
                      + dd_22[k];

            t_21[k] = -ab_y[k] * dp_12[k]
                      + dd_25[k];

            t_22[k] = -ab_y[k] * dp_13[k]
                      + dd_27[k];

            t_23[k] = -ab_y[k] * dp_14[k]
                      + dd_28[k];
        }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, ab_y, ab_z, dp_15, dp_16, dp_17, \
                         dd_31, dd_32, dd_33, dd_34, dd_35 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_24[k] = -ab_y[k] * dp_15[k]
                      + dd_31[k];

            t_25[k] = -ab_y[k] * dp_16[k]
                      + dd_33[k];

            t_26[k] = -ab_y[k] * dp_17[k]
                      + dd_34[k];

            t_27[k] = -ab_z[k] * dp_15[k]
                      + dd_32[k];

            t_28[k] = -ab_z[k] * dp_16[k]
                      + dd_34[k];

            t_29[k] = -ab_z[k] * dp_17[k]
                      + dd_35[k];
        }
    }
}

auto
compute_hrr_fp(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t fs, const size_t gs, const size_t ncomps, const size_t nmax) -> void
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

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *fs_0 = buffer.data(fs + 0 * ncomps + c);
        const auto *fs_1 = buffer.data(fs + 1 * ncomps + c);
        const auto *fs_2 = buffer.data(fs + 2 * ncomps + c);
        const auto *fs_3 = buffer.data(fs + 3 * ncomps + c);
        const auto *fs_4 = buffer.data(fs + 4 * ncomps + c);
        const auto *fs_5 = buffer.data(fs + 5 * ncomps + c);
        const auto *fs_6 = buffer.data(fs + 6 * ncomps + c);
        const auto *fs_7 = buffer.data(fs + 7 * ncomps + c);
        const auto *fs_8 = buffer.data(fs + 8 * ncomps + c);
        const auto *fs_9 = buffer.data(fs + 9 * ncomps + c);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, ab_x, ab_y, ab_z, fs_0, fs_1, gs_0, \
                         gs_1, gs_2, gs_3, gs_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * fs_0[k]
                     + gs_0[k];

            t_1[k] = ab_y[k] * fs_0[k]
                     + gs_1[k];

            t_2[k] = ab_z[k] * fs_0[k]
                     + gs_2[k];

            t_3[k] = ab_x[k] * fs_1[k]
                     + gs_1[k];

            t_4[k] = ab_y[k] * fs_1[k]
                     + gs_3[k];

            t_5[k] = ab_z[k] * fs_1[k]
                     + gs_4[k];
        }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, ab_x, ab_y, ab_z, fs_2, fs_3, gs_2, gs_3, \
                         gs_4, gs_5, gs_6 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_6[k] = ab_x[k] * fs_2[k]
                     + gs_2[k];

            t_7[k] = ab_y[k] * fs_2[k]
                     + gs_4[k];

            t_8[k] = ab_z[k] * fs_2[k]
                     + gs_5[k];

            t_9[k] = ab_x[k] * fs_3[k]
                     + gs_3[k];

            t_10[k] = ab_y[k] * fs_3[k]
                      + gs_6[k];
        }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, ab_x, ab_y, ab_z, fs_3, fs_4, \
                         fs_5, gs_4, gs_5, gs_7, gs_8 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_11[k] = ab_z[k] * fs_3[k]
                      + gs_7[k];

            t_12[k] = ab_x[k] * fs_4[k]
                      + gs_4[k];

            t_13[k] = ab_y[k] * fs_4[k]
                      + gs_7[k];

            t_14[k] = ab_z[k] * fs_4[k]
                      + gs_8[k];

            t_15[k] = ab_x[k] * fs_5[k]
                      + gs_5[k];

            t_16[k] = ab_y[k] * fs_5[k]
                      + gs_8[k];
        }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, t_21, ab_x, ab_y, ab_z, fs_5, fs_6, fs_7, \
                         gs_6, gs_7, gs_9, gs_10, gs_11 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_17[k] = ab_z[k] * fs_5[k]
                      + gs_9[k];

            t_18[k] = ab_x[k] * fs_6[k]
                      + gs_6[k];

            t_19[k] = ab_y[k] * fs_6[k]
                      + gs_10[k];

            t_20[k] = ab_z[k] * fs_6[k]
                      + gs_11[k];

            t_21[k] = ab_x[k] * fs_7[k]
                      + gs_7[k];
        }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, ab_x, ab_y, ab_z, fs_7, fs_8, gs_8, \
                         gs_11, gs_12, gs_13 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_22[k] = ab_y[k] * fs_7[k]
                      + gs_11[k];

            t_23[k] = ab_z[k] * fs_7[k]
                      + gs_12[k];

            t_24[k] = ab_x[k] * fs_8[k]
                      + gs_8[k];

            t_25[k] = ab_y[k] * fs_8[k]
                      + gs_12[k];

            t_26[k] = ab_z[k] * fs_8[k]
                      + gs_13[k];
        }

#pragma omp simd aligned(t_27, t_28, t_29, ab_x, ab_y, ab_z, fs_9, gs_9, gs_13, \
                         gs_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_27[k] = ab_x[k] * fs_9[k]
                      + gs_9[k];

            t_28[k] = ab_y[k] * fs_9[k]
                      + gs_13[k];

            t_29[k] = ab_z[k] * fs_9[k]
                      + gs_14[k];
        }
    }
}

}  // namespace simdtrf
