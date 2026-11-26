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


#include "SimdTransferDD.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_dd_out_of_first(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                            const size_t target, const size_t dp, const size_t fp,
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

        const auto *fp_0 = buffer.data(fp + 0 * ncomps + c);
        const auto *fp_1 = buffer.data(fp + 1 * ncomps + c);
        const auto *fp_2 = buffer.data(fp + 2 * ncomps + c);
        const auto *fp_3 = buffer.data(fp + 3 * ncomps + c);
        const auto *fp_4 = buffer.data(fp + 4 * ncomps + c);
        const auto *fp_5 = buffer.data(fp + 5 * ncomps + c);
        const auto *fp_6 = buffer.data(fp + 6 * ncomps + c);
        const auto *fp_7 = buffer.data(fp + 7 * ncomps + c);
        const auto *fp_8 = buffer.data(fp + 8 * ncomps + c);
        const auto *fp_9 = buffer.data(fp + 9 * ncomps + c);
        const auto *fp_10 = buffer.data(fp + 10 * ncomps + c);
        const auto *fp_11 = buffer.data(fp + 11 * ncomps + c);
        const auto *fp_12 = buffer.data(fp + 12 * ncomps + c);
        const auto *fp_13 = buffer.data(fp + 13 * ncomps + c);
        const auto *fp_14 = buffer.data(fp + 14 * ncomps + c);
        const auto *fp_15 = buffer.data(fp + 15 * ncomps + c);
        const auto *fp_16 = buffer.data(fp + 16 * ncomps + c);
        const auto *fp_17 = buffer.data(fp + 17 * ncomps + c);
        const auto *fp_19 = buffer.data(fp + 19 * ncomps + c);
        const auto *fp_20 = buffer.data(fp + 20 * ncomps + c);
        const auto *fp_22 = buffer.data(fp + 22 * ncomps + c);
        const auto *fp_23 = buffer.data(fp + 23 * ncomps + c);
        const auto *fp_25 = buffer.data(fp + 25 * ncomps + c);
        const auto *fp_26 = buffer.data(fp + 26 * ncomps + c);
        const auto *fp_29 = buffer.data(fp + 29 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, ab_y, dp_0, dp_1, dp_2, fp_0, fp_1, \
                         fp_2, fp_4, fp_5 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * dp_0[k]
                     + fp_0[k];

            t_1[k] = ab_x[k] * dp_1[k]
                     + fp_1[k];

            t_2[k] = ab_x[k] * dp_2[k]
                     + fp_2[k];

            t_3[k] = ab_y[k] * dp_1[k]
                     + fp_4[k];

            t_4[k] = ab_y[k] * dp_2[k]
                     + fp_5[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, ab_x, ab_z, dp_2, dp_3, dp_4, dp_5, fp_3, fp_4, \
                         fp_5, fp_8 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = ab_z[k] * dp_2[k]
                     + fp_8[k];

            t_6[k] = ab_x[k] * dp_3[k]
                     + fp_3[k];

            t_7[k] = ab_x[k] * dp_4[k]
                     + fp_4[k];

            t_8[k] = ab_x[k] * dp_5[k]
                     + fp_5[k];
        }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, ab_x, ab_y, ab_z, dp_4, dp_5, dp_6, fp_6, \
                         fp_10, fp_11, fp_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_9[k] = ab_y[k] * dp_4[k]
                     + fp_10[k];

            t_10[k] = ab_y[k] * dp_5[k]
                      + fp_11[k];

            t_11[k] = ab_z[k] * dp_5[k]
                      + fp_14[k];

            t_12[k] = ab_x[k] * dp_6[k]
                      + fp_6[k];
        }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, ab_x, ab_y, ab_z, dp_7, dp_8, fp_7, \
                         fp_8, fp_13, fp_14, fp_17 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_13[k] = ab_x[k] * dp_7[k]
                      + fp_7[k];

            t_14[k] = ab_x[k] * dp_8[k]
                      + fp_8[k];

            t_15[k] = ab_y[k] * dp_7[k]
                      + fp_13[k];

            t_16[k] = ab_y[k] * dp_8[k]
                      + fp_14[k];

            t_17[k] = ab_z[k] * dp_8[k]
                      + fp_17[k];
        }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, ab_x, ab_y, dp_9, dp_10, dp_11, fp_9, \
                         fp_10, fp_11, fp_19, fp_20 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_18[k] = ab_x[k] * dp_9[k]
                      + fp_9[k];

            t_19[k] = ab_x[k] * dp_10[k]
                      + fp_10[k];

            t_20[k] = ab_x[k] * dp_11[k]
                      + fp_11[k];

            t_21[k] = ab_y[k] * dp_10[k]
                      + fp_19[k];

            t_22[k] = ab_y[k] * dp_11[k]
                      + fp_20[k];
        }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, ab_x, ab_z, dp_11, dp_12, dp_13, dp_14, \
                         fp_12, fp_13, fp_14, fp_23 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_23[k] = ab_z[k] * dp_11[k]
                      + fp_23[k];

            t_24[k] = ab_x[k] * dp_12[k]
                      + fp_12[k];

            t_25[k] = ab_x[k] * dp_13[k]
                      + fp_13[k];

            t_26[k] = ab_x[k] * dp_14[k]
                      + fp_14[k];
        }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, ab_x, ab_y, ab_z, dp_13, dp_14, dp_15, fp_15, \
                         fp_22, fp_23, fp_26 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_27[k] = ab_y[k] * dp_13[k]
                      + fp_22[k];

            t_28[k] = ab_y[k] * dp_14[k]
                      + fp_23[k];

            t_29[k] = ab_z[k] * dp_14[k]
                      + fp_26[k];

            t_30[k] = ab_x[k] * dp_15[k]
                      + fp_15[k];
        }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, ab_x, ab_y, ab_z, dp_16, dp_17, fp_16, \
                         fp_17, fp_25, fp_26, fp_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_31[k] = ab_x[k] * dp_16[k]
                      + fp_16[k];

            t_32[k] = ab_x[k] * dp_17[k]
                      + fp_17[k];

            t_33[k] = ab_y[k] * dp_16[k]
                      + fp_25[k];

            t_34[k] = ab_y[k] * dp_17[k]
                      + fp_26[k];

            t_35[k] = ab_z[k] * dp_17[k]
                      + fp_29[k];
        }
    }
}

auto
compute_hrr_dd(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t pd, const size_t pf, const size_t ncomps, const size_t nmax) -> void
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

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *pd_0 = buffer.data(pd + 0 * ncomps + c);
        const auto *pd_1 = buffer.data(pd + 1 * ncomps + c);
        const auto *pd_2 = buffer.data(pd + 2 * ncomps + c);
        const auto *pd_3 = buffer.data(pd + 3 * ncomps + c);
        const auto *pd_4 = buffer.data(pd + 4 * ncomps + c);
        const auto *pd_5 = buffer.data(pd + 5 * ncomps + c);
        const auto *pd_6 = buffer.data(pd + 6 * ncomps + c);
        const auto *pd_7 = buffer.data(pd + 7 * ncomps + c);
        const auto *pd_8 = buffer.data(pd + 8 * ncomps + c);
        const auto *pd_9 = buffer.data(pd + 9 * ncomps + c);
        const auto *pd_10 = buffer.data(pd + 10 * ncomps + c);
        const auto *pd_11 = buffer.data(pd + 11 * ncomps + c);
        const auto *pd_12 = buffer.data(pd + 12 * ncomps + c);
        const auto *pd_13 = buffer.data(pd + 13 * ncomps + c);
        const auto *pd_14 = buffer.data(pd + 14 * ncomps + c);
        const auto *pd_15 = buffer.data(pd + 15 * ncomps + c);
        const auto *pd_16 = buffer.data(pd + 16 * ncomps + c);
        const auto *pd_17 = buffer.data(pd + 17 * ncomps + c);

        const auto *pf_0 = buffer.data(pf + 0 * ncomps + c);
        const auto *pf_1 = buffer.data(pf + 1 * ncomps + c);
        const auto *pf_2 = buffer.data(pf + 2 * ncomps + c);
        const auto *pf_3 = buffer.data(pf + 3 * ncomps + c);
        const auto *pf_4 = buffer.data(pf + 4 * ncomps + c);
        const auto *pf_5 = buffer.data(pf + 5 * ncomps + c);
        const auto *pf_10 = buffer.data(pf + 10 * ncomps + c);
        const auto *pf_11 = buffer.data(pf + 11 * ncomps + c);
        const auto *pf_12 = buffer.data(pf + 12 * ncomps + c);
        const auto *pf_13 = buffer.data(pf + 13 * ncomps + c);
        const auto *pf_14 = buffer.data(pf + 14 * ncomps + c);
        const auto *pf_15 = buffer.data(pf + 15 * ncomps + c);
        const auto *pf_16 = buffer.data(pf + 16 * ncomps + c);
        const auto *pf_17 = buffer.data(pf + 17 * ncomps + c);
        const auto *pf_18 = buffer.data(pf + 18 * ncomps + c);
        const auto *pf_20 = buffer.data(pf + 20 * ncomps + c);
        const auto *pf_21 = buffer.data(pf + 21 * ncomps + c);
        const auto *pf_22 = buffer.data(pf + 22 * ncomps + c);
        const auto *pf_23 = buffer.data(pf + 23 * ncomps + c);
        const auto *pf_24 = buffer.data(pf + 24 * ncomps + c);
        const auto *pf_25 = buffer.data(pf + 25 * ncomps + c);
        const auto *pf_26 = buffer.data(pf + 26 * ncomps + c);
        const auto *pf_27 = buffer.data(pf + 27 * ncomps + c);
        const auto *pf_28 = buffer.data(pf + 28 * ncomps + c);
        const auto *pf_29 = buffer.data(pf + 29 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, pd_0, pd_1, pd_2, pd_3, pd_4, pf_0, \
                         pf_1, pf_2, pf_3, pf_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = -ab_x[k] * pd_0[k]
                     + pf_0[k];

            t_1[k] = -ab_x[k] * pd_1[k]
                     + pf_1[k];

            t_2[k] = -ab_x[k] * pd_2[k]
                     + pf_2[k];

            t_3[k] = -ab_x[k] * pd_3[k]
                     + pf_3[k];

            t_4[k] = -ab_x[k] * pd_4[k]
                     + pf_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, pd_5, pd_6, pd_7, pd_8, pd_9, pf_5, \
                         pf_10, pf_11, pf_12, pf_13 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = -ab_x[k] * pd_5[k]
                     + pf_5[k];

            t_6[k] = -ab_x[k] * pd_6[k]
                     + pf_10[k];

            t_7[k] = -ab_x[k] * pd_7[k]
                     + pf_11[k];

            t_8[k] = -ab_x[k] * pd_8[k]
                     + pf_12[k];

            t_9[k] = -ab_x[k] * pd_9[k]
                     + pf_13[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, pd_10, pd_11, pd_12, pd_13, \
                         pd_14, pf_14, pf_15, pf_20, pf_21, pf_22 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = -ab_x[k] * pd_10[k]
                      + pf_14[k];

            t_11[k] = -ab_x[k] * pd_11[k]
                      + pf_15[k];

            t_12[k] = -ab_x[k] * pd_12[k]
                      + pf_20[k];

            t_13[k] = -ab_x[k] * pd_13[k]
                      + pf_21[k];

            t_14[k] = -ab_x[k] * pd_14[k]
                      + pf_22[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, ab_x, ab_y, pd_6, pd_15, pd_16, pd_17, pf_11, \
                         pf_23, pf_24, pf_25 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = -ab_x[k] * pd_15[k]
                      + pf_23[k];

            t_16[k] = -ab_x[k] * pd_16[k]
                      + pf_24[k];

            t_17[k] = -ab_x[k] * pd_17[k]
                      + pf_25[k];

            t_18[k] = -ab_y[k] * pd_6[k]
                      + pf_11[k];
        }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, ab_y, pd_7, pd_8, pd_9, pd_10, pd_11, \
                         pf_13, pf_14, pf_16, pf_17, pf_18 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_19[k] = -ab_y[k] * pd_7[k]
                      + pf_13[k];

            t_20[k] = -ab_y[k] * pd_8[k]
                      + pf_14[k];

            t_21[k] = -ab_y[k] * pd_9[k]
                      + pf_16[k];

            t_22[k] = -ab_y[k] * pd_10[k]
                      + pf_17[k];

            t_23[k] = -ab_y[k] * pd_11[k]
                      + pf_18[k];
        }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, ab_y, pd_12, pd_13, pd_14, pd_15, \
                         pd_16, pf_21, pf_23, pf_24, pf_26, pf_27 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_24[k] = -ab_y[k] * pd_12[k]
                      + pf_21[k];

            t_25[k] = -ab_y[k] * pd_13[k]
                      + pf_23[k];

            t_26[k] = -ab_y[k] * pd_14[k]
                      + pf_24[k];

            t_27[k] = -ab_y[k] * pd_15[k]
                      + pf_26[k];

            t_28[k] = -ab_y[k] * pd_16[k]
                      + pf_27[k];
        }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, ab_y, ab_z, pd_12, pd_13, pd_14, pd_17, \
                         pf_22, pf_24, pf_25, pf_28 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_29[k] = -ab_y[k] * pd_17[k]
                      + pf_28[k];

            t_30[k] = -ab_z[k] * pd_12[k]
                      + pf_22[k];

            t_31[k] = -ab_z[k] * pd_13[k]
                      + pf_24[k];

            t_32[k] = -ab_z[k] * pd_14[k]
                      + pf_25[k];
        }

#pragma omp simd aligned(t_33, t_34, t_35, ab_z, pd_15, pd_16, pd_17, pf_27, pf_28, \
                         pf_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_33[k] = -ab_z[k] * pd_15[k]
                      + pf_27[k];

            t_34[k] = -ab_z[k] * pd_16[k]
                      + pf_28[k];

            t_35[k] = -ab_z[k] * pd_17[k]
                      + pf_29[k];
        }
    }
}

}  // namespace simdtrf
