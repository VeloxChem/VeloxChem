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


#include "SimdTransferGeom010XFD.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_geom_010x_fd_out_of_second(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                       const size_t target, const size_t dd_1, const size_t dd_0,
                                       const size_t df_1, const size_t ncomps,
                                       const size_t nmax) -> void
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
        auto *t_45 = buffer.data(target + 45 * ncomps + c);
        auto *t_46 = buffer.data(target + 46 * ncomps + c);
        auto *t_47 = buffer.data(target + 47 * ncomps + c);
        auto *t_48 = buffer.data(target + 48 * ncomps + c);
        auto *t_49 = buffer.data(target + 49 * ncomps + c);
        auto *t_50 = buffer.data(target + 50 * ncomps + c);
        auto *t_51 = buffer.data(target + 51 * ncomps + c);
        auto *t_52 = buffer.data(target + 52 * ncomps + c);
        auto *t_53 = buffer.data(target + 53 * ncomps + c);
        auto *t_54 = buffer.data(target + 54 * ncomps + c);
        auto *t_55 = buffer.data(target + 55 * ncomps + c);
        auto *t_56 = buffer.data(target + 56 * ncomps + c);
        auto *t_57 = buffer.data(target + 57 * ncomps + c);
        auto *t_58 = buffer.data(target + 58 * ncomps + c);
        auto *t_59 = buffer.data(target + 59 * ncomps + c);

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *dd_1_0 = buffer.data(dd_1 + 0 * ncomps + c);
        const auto *dd_1_1 = buffer.data(dd_1 + 1 * ncomps + c);
        const auto *dd_1_2 = buffer.data(dd_1 + 2 * ncomps + c);
        const auto *dd_1_3 = buffer.data(dd_1 + 3 * ncomps + c);
        const auto *dd_1_4 = buffer.data(dd_1 + 4 * ncomps + c);
        const auto *dd_1_5 = buffer.data(dd_1 + 5 * ncomps + c);
        const auto *dd_1_6 = buffer.data(dd_1 + 6 * ncomps + c);
        const auto *dd_1_7 = buffer.data(dd_1 + 7 * ncomps + c);
        const auto *dd_1_8 = buffer.data(dd_1 + 8 * ncomps + c);
        const auto *dd_1_9 = buffer.data(dd_1 + 9 * ncomps + c);
        const auto *dd_1_10 = buffer.data(dd_1 + 10 * ncomps + c);
        const auto *dd_1_11 = buffer.data(dd_1 + 11 * ncomps + c);
        const auto *dd_1_12 = buffer.data(dd_1 + 12 * ncomps + c);
        const auto *dd_1_13 = buffer.data(dd_1 + 13 * ncomps + c);
        const auto *dd_1_14 = buffer.data(dd_1 + 14 * ncomps + c);
        const auto *dd_1_15 = buffer.data(dd_1 + 15 * ncomps + c);
        const auto *dd_1_16 = buffer.data(dd_1 + 16 * ncomps + c);
        const auto *dd_1_17 = buffer.data(dd_1 + 17 * ncomps + c);
        const auto *dd_1_18 = buffer.data(dd_1 + 18 * ncomps + c);
        const auto *dd_1_19 = buffer.data(dd_1 + 19 * ncomps + c);
        const auto *dd_1_20 = buffer.data(dd_1 + 20 * ncomps + c);
        const auto *dd_1_21 = buffer.data(dd_1 + 21 * ncomps + c);
        const auto *dd_1_22 = buffer.data(dd_1 + 22 * ncomps + c);
        const auto *dd_1_23 = buffer.data(dd_1 + 23 * ncomps + c);
        const auto *dd_1_24 = buffer.data(dd_1 + 24 * ncomps + c);
        const auto *dd_1_25 = buffer.data(dd_1 + 25 * ncomps + c);
        const auto *dd_1_26 = buffer.data(dd_1 + 26 * ncomps + c);
        const auto *dd_1_27 = buffer.data(dd_1 + 27 * ncomps + c);
        const auto *dd_1_28 = buffer.data(dd_1 + 28 * ncomps + c);
        const auto *dd_1_29 = buffer.data(dd_1 + 29 * ncomps + c);
        const auto *dd_1_30 = buffer.data(dd_1 + 30 * ncomps + c);
        const auto *dd_1_31 = buffer.data(dd_1 + 31 * ncomps + c);
        const auto *dd_1_32 = buffer.data(dd_1 + 32 * ncomps + c);
        const auto *dd_1_33 = buffer.data(dd_1 + 33 * ncomps + c);
        const auto *dd_1_34 = buffer.data(dd_1 + 34 * ncomps + c);
        const auto *dd_1_35 = buffer.data(dd_1 + 35 * ncomps + c);

        const auto *dd_0_0 = buffer.data(dd_0 + 0 * ncomps + c);
        const auto *dd_0_1 = buffer.data(dd_0 + 1 * ncomps + c);
        const auto *dd_0_2 = buffer.data(dd_0 + 2 * ncomps + c);
        const auto *dd_0_3 = buffer.data(dd_0 + 3 * ncomps + c);
        const auto *dd_0_4 = buffer.data(dd_0 + 4 * ncomps + c);
        const auto *dd_0_5 = buffer.data(dd_0 + 5 * ncomps + c);
        const auto *dd_0_6 = buffer.data(dd_0 + 6 * ncomps + c);
        const auto *dd_0_7 = buffer.data(dd_0 + 7 * ncomps + c);
        const auto *dd_0_8 = buffer.data(dd_0 + 8 * ncomps + c);
        const auto *dd_0_9 = buffer.data(dd_0 + 9 * ncomps + c);
        const auto *dd_0_10 = buffer.data(dd_0 + 10 * ncomps + c);
        const auto *dd_0_11 = buffer.data(dd_0 + 11 * ncomps + c);
        const auto *dd_0_12 = buffer.data(dd_0 + 12 * ncomps + c);
        const auto *dd_0_13 = buffer.data(dd_0 + 13 * ncomps + c);
        const auto *dd_0_14 = buffer.data(dd_0 + 14 * ncomps + c);
        const auto *dd_0_15 = buffer.data(dd_0 + 15 * ncomps + c);
        const auto *dd_0_16 = buffer.data(dd_0 + 16 * ncomps + c);
        const auto *dd_0_17 = buffer.data(dd_0 + 17 * ncomps + c);
        const auto *dd_0_18 = buffer.data(dd_0 + 18 * ncomps + c);
        const auto *dd_0_19 = buffer.data(dd_0 + 19 * ncomps + c);
        const auto *dd_0_20 = buffer.data(dd_0 + 20 * ncomps + c);
        const auto *dd_0_21 = buffer.data(dd_0 + 21 * ncomps + c);
        const auto *dd_0_22 = buffer.data(dd_0 + 22 * ncomps + c);
        const auto *dd_0_23 = buffer.data(dd_0 + 23 * ncomps + c);
        const auto *dd_0_24 = buffer.data(dd_0 + 24 * ncomps + c);
        const auto *dd_0_25 = buffer.data(dd_0 + 25 * ncomps + c);
        const auto *dd_0_26 = buffer.data(dd_0 + 26 * ncomps + c);
        const auto *dd_0_27 = buffer.data(dd_0 + 27 * ncomps + c);
        const auto *dd_0_28 = buffer.data(dd_0 + 28 * ncomps + c);
        const auto *dd_0_29 = buffer.data(dd_0 + 29 * ncomps + c);
        const auto *dd_0_30 = buffer.data(dd_0 + 30 * ncomps + c);
        const auto *dd_0_31 = buffer.data(dd_0 + 31 * ncomps + c);
        const auto *dd_0_32 = buffer.data(dd_0 + 32 * ncomps + c);
        const auto *dd_0_33 = buffer.data(dd_0 + 33 * ncomps + c);
        const auto *dd_0_34 = buffer.data(dd_0 + 34 * ncomps + c);
        const auto *dd_0_35 = buffer.data(dd_0 + 35 * ncomps + c);

        const auto *df_1_0 = buffer.data(df_1 + 0 * ncomps + c);
        const auto *df_1_1 = buffer.data(df_1 + 1 * ncomps + c);
        const auto *df_1_2 = buffer.data(df_1 + 2 * ncomps + c);
        const auto *df_1_3 = buffer.data(df_1 + 3 * ncomps + c);
        const auto *df_1_4 = buffer.data(df_1 + 4 * ncomps + c);
        const auto *df_1_5 = buffer.data(df_1 + 5 * ncomps + c);
        const auto *df_1_10 = buffer.data(df_1 + 10 * ncomps + c);
        const auto *df_1_11 = buffer.data(df_1 + 11 * ncomps + c);
        const auto *df_1_12 = buffer.data(df_1 + 12 * ncomps + c);
        const auto *df_1_13 = buffer.data(df_1 + 13 * ncomps + c);
        const auto *df_1_14 = buffer.data(df_1 + 14 * ncomps + c);
        const auto *df_1_15 = buffer.data(df_1 + 15 * ncomps + c);
        const auto *df_1_20 = buffer.data(df_1 + 20 * ncomps + c);
        const auto *df_1_21 = buffer.data(df_1 + 21 * ncomps + c);
        const auto *df_1_22 = buffer.data(df_1 + 22 * ncomps + c);
        const auto *df_1_23 = buffer.data(df_1 + 23 * ncomps + c);
        const auto *df_1_24 = buffer.data(df_1 + 24 * ncomps + c);
        const auto *df_1_25 = buffer.data(df_1 + 25 * ncomps + c);
        const auto *df_1_30 = buffer.data(df_1 + 30 * ncomps + c);
        const auto *df_1_31 = buffer.data(df_1 + 31 * ncomps + c);
        const auto *df_1_32 = buffer.data(df_1 + 32 * ncomps + c);
        const auto *df_1_33 = buffer.data(df_1 + 33 * ncomps + c);
        const auto *df_1_34 = buffer.data(df_1 + 34 * ncomps + c);
        const auto *df_1_35 = buffer.data(df_1 + 35 * ncomps + c);
        const auto *df_1_36 = buffer.data(df_1 + 36 * ncomps + c);
        const auto *df_1_37 = buffer.data(df_1 + 37 * ncomps + c);
        const auto *df_1_38 = buffer.data(df_1 + 38 * ncomps + c);
        const auto *df_1_40 = buffer.data(df_1 + 40 * ncomps + c);
        const auto *df_1_41 = buffer.data(df_1 + 41 * ncomps + c);
        const auto *df_1_42 = buffer.data(df_1 + 42 * ncomps + c);
        const auto *df_1_43 = buffer.data(df_1 + 43 * ncomps + c);
        const auto *df_1_44 = buffer.data(df_1 + 44 * ncomps + c);
        const auto *df_1_45 = buffer.data(df_1 + 45 * ncomps + c);
        const auto *df_1_46 = buffer.data(df_1 + 46 * ncomps + c);
        const auto *df_1_47 = buffer.data(df_1 + 47 * ncomps + c);
        const auto *df_1_48 = buffer.data(df_1 + 48 * ncomps + c);
        const auto *df_1_50 = buffer.data(df_1 + 50 * ncomps + c);
        const auto *df_1_51 = buffer.data(df_1 + 51 * ncomps + c);
        const auto *df_1_52 = buffer.data(df_1 + 52 * ncomps + c);
        const auto *df_1_53 = buffer.data(df_1 + 53 * ncomps + c);
        const auto *df_1_54 = buffer.data(df_1 + 54 * ncomps + c);
        const auto *df_1_55 = buffer.data(df_1 + 55 * ncomps + c);
        const auto *df_1_56 = buffer.data(df_1 + 56 * ncomps + c);
        const auto *df_1_57 = buffer.data(df_1 + 57 * ncomps + c);
        const auto *df_1_58 = buffer.data(df_1 + 58 * ncomps + c);
        const auto *df_1_59 = buffer.data(df_1 + 59 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, ab_x, dd_1_0, dd_1_1, dd_1_2, dd_0_0, dd_0_1, dd_0_2, \
                         df_1_0, df_1_1, df_1_2 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = -ab_x[k] * dd_1_0[k]
                     + dd_0_0[k]
                     + df_1_0[k];

            t_1[k] = -ab_x[k] * dd_1_1[k]
                     + dd_0_1[k]
                     + df_1_1[k];

            t_2[k] = -ab_x[k] * dd_1_2[k]
                     + dd_0_2[k]
                     + df_1_2[k];
        }

#pragma omp simd aligned(t_3, t_4, t_5, ab_x, dd_1_3, dd_1_4, dd_1_5, dd_0_3, dd_0_4, dd_0_5, \
                         df_1_3, df_1_4, df_1_5 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_3[k] = -ab_x[k] * dd_1_3[k]
                     + dd_0_3[k]
                     + df_1_3[k];

            t_4[k] = -ab_x[k] * dd_1_4[k]
                     + dd_0_4[k]
                     + df_1_4[k];

            t_5[k] = -ab_x[k] * dd_1_5[k]
                     + dd_0_5[k]
                     + df_1_5[k];
        }

#pragma omp simd aligned(t_6, t_7, t_8, ab_x, dd_1_6, dd_1_7, dd_1_8, dd_0_6, dd_0_7, dd_0_8, \
                         df_1_10, df_1_11, df_1_12 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_6[k] = -ab_x[k] * dd_1_6[k]
                     + dd_0_6[k]
                     + df_1_10[k];

            t_7[k] = -ab_x[k] * dd_1_7[k]
                     + dd_0_7[k]
                     + df_1_11[k];

            t_8[k] = -ab_x[k] * dd_1_8[k]
                     + dd_0_8[k]
                     + df_1_12[k];
        }

#pragma omp simd aligned(t_9, t_10, t_11, ab_x, dd_1_9, dd_1_10, dd_1_11, dd_0_9, dd_0_10, \
                         dd_0_11, df_1_13, df_1_14, df_1_15 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_9[k] = -ab_x[k] * dd_1_9[k]
                     + dd_0_9[k]
                     + df_1_13[k];

            t_10[k] = -ab_x[k] * dd_1_10[k]
                      + dd_0_10[k]
                      + df_1_14[k];

            t_11[k] = -ab_x[k] * dd_1_11[k]
                      + dd_0_11[k]
                      + df_1_15[k];
        }

#pragma omp simd aligned(t_12, t_13, t_14, ab_x, dd_1_12, dd_1_13, dd_1_14, dd_0_12, dd_0_13, \
                         dd_0_14, df_1_20, df_1_21, df_1_22 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_12[k] = -ab_x[k] * dd_1_12[k]
                      + dd_0_12[k]
                      + df_1_20[k];

            t_13[k] = -ab_x[k] * dd_1_13[k]
                      + dd_0_13[k]
                      + df_1_21[k];

            t_14[k] = -ab_x[k] * dd_1_14[k]
                      + dd_0_14[k]
                      + df_1_22[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, ab_x, dd_1_15, dd_1_16, dd_1_17, dd_0_15, dd_0_16, \
                         dd_0_17, df_1_23, df_1_24, df_1_25 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = -ab_x[k] * dd_1_15[k]
                      + dd_0_15[k]
                      + df_1_23[k];

            t_16[k] = -ab_x[k] * dd_1_16[k]
                      + dd_0_16[k]
                      + df_1_24[k];

            t_17[k] = -ab_x[k] * dd_1_17[k]
                      + dd_0_17[k]
                      + df_1_25[k];
        }

#pragma omp simd aligned(t_18, t_19, t_20, ab_x, dd_1_18, dd_1_19, dd_1_20, dd_0_18, dd_0_19, \
                         dd_0_20, df_1_30, df_1_31, df_1_32 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_18[k] = -ab_x[k] * dd_1_18[k]
                      + dd_0_18[k]
                      + df_1_30[k];

            t_19[k] = -ab_x[k] * dd_1_19[k]
                      + dd_0_19[k]
                      + df_1_31[k];

            t_20[k] = -ab_x[k] * dd_1_20[k]
                      + dd_0_20[k]
                      + df_1_32[k];
        }

#pragma omp simd aligned(t_21, t_22, t_23, ab_x, dd_1_21, dd_1_22, dd_1_23, dd_0_21, dd_0_22, \
                         dd_0_23, df_1_33, df_1_34, df_1_35 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_21[k] = -ab_x[k] * dd_1_21[k]
                      + dd_0_21[k]
                      + df_1_33[k];

            t_22[k] = -ab_x[k] * dd_1_22[k]
                      + dd_0_22[k]
                      + df_1_34[k];

            t_23[k] = -ab_x[k] * dd_1_23[k]
                      + dd_0_23[k]
                      + df_1_35[k];
        }

#pragma omp simd aligned(t_24, t_25, t_26, ab_x, dd_1_24, dd_1_25, dd_1_26, dd_0_24, dd_0_25, \
                         dd_0_26, df_1_40, df_1_41, df_1_42 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_24[k] = -ab_x[k] * dd_1_24[k]
                      + dd_0_24[k]
                      + df_1_40[k];

            t_25[k] = -ab_x[k] * dd_1_25[k]
                      + dd_0_25[k]
                      + df_1_41[k];

            t_26[k] = -ab_x[k] * dd_1_26[k]
                      + dd_0_26[k]
                      + df_1_42[k];
        }

#pragma omp simd aligned(t_27, t_28, t_29, ab_x, dd_1_27, dd_1_28, dd_1_29, dd_0_27, dd_0_28, \
                         dd_0_29, df_1_43, df_1_44, df_1_45 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_27[k] = -ab_x[k] * dd_1_27[k]
                      + dd_0_27[k]
                      + df_1_43[k];

            t_28[k] = -ab_x[k] * dd_1_28[k]
                      + dd_0_28[k]
                      + df_1_44[k];

            t_29[k] = -ab_x[k] * dd_1_29[k]
                      + dd_0_29[k]
                      + df_1_45[k];
        }

#pragma omp simd aligned(t_30, t_31, t_32, ab_x, dd_1_30, dd_1_31, dd_1_32, dd_0_30, dd_0_31, \
                         dd_0_32, df_1_50, df_1_51, df_1_52 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_30[k] = -ab_x[k] * dd_1_30[k]
                      + dd_0_30[k]
                      + df_1_50[k];

            t_31[k] = -ab_x[k] * dd_1_31[k]
                      + dd_0_31[k]
                      + df_1_51[k];

            t_32[k] = -ab_x[k] * dd_1_32[k]
                      + dd_0_32[k]
                      + df_1_52[k];
        }

#pragma omp simd aligned(t_33, t_34, t_35, ab_x, dd_1_33, dd_1_34, dd_1_35, dd_0_33, dd_0_34, \
                         dd_0_35, df_1_53, df_1_54, df_1_55 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_33[k] = -ab_x[k] * dd_1_33[k]
                      + dd_0_33[k]
                      + df_1_53[k];

            t_34[k] = -ab_x[k] * dd_1_34[k]
                      + dd_0_34[k]
                      + df_1_54[k];

            t_35[k] = -ab_x[k] * dd_1_35[k]
                      + dd_0_35[k]
                      + df_1_55[k];
        }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, ab_y, dd_1_18, dd_1_19, dd_1_20, \
                         dd_1_21, dd_1_22, df_1_31, df_1_33, df_1_34, df_1_36, \
                         df_1_37 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_36[k] = -ab_y[k] * dd_1_18[k]
                      + df_1_31[k];

            t_37[k] = -ab_y[k] * dd_1_19[k]
                      + df_1_33[k];

            t_38[k] = -ab_y[k] * dd_1_20[k]
                      + df_1_34[k];

            t_39[k] = -ab_y[k] * dd_1_21[k]
                      + df_1_36[k];

            t_40[k] = -ab_y[k] * dd_1_22[k]
                      + df_1_37[k];
        }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, t_45, ab_y, dd_1_23, dd_1_24, dd_1_25, \
                         dd_1_26, dd_1_27, df_1_38, df_1_41, df_1_43, df_1_44, \
                         df_1_46 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_41[k] = -ab_y[k] * dd_1_23[k]
                      + df_1_38[k];

            t_42[k] = -ab_y[k] * dd_1_24[k]
                      + df_1_41[k];

            t_43[k] = -ab_y[k] * dd_1_25[k]
                      + df_1_43[k];

            t_44[k] = -ab_y[k] * dd_1_26[k]
                      + df_1_44[k];

            t_45[k] = -ab_y[k] * dd_1_27[k]
                      + df_1_46[k];
        }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, t_50, ab_y, dd_1_28, dd_1_29, dd_1_30, \
                         dd_1_31, dd_1_32, df_1_47, df_1_48, df_1_51, df_1_53, \
                         df_1_54 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_46[k] = -ab_y[k] * dd_1_28[k]
                      + df_1_47[k];

            t_47[k] = -ab_y[k] * dd_1_29[k]
                      + df_1_48[k];

            t_48[k] = -ab_y[k] * dd_1_30[k]
                      + df_1_51[k];

            t_49[k] = -ab_y[k] * dd_1_31[k]
                      + df_1_53[k];

            t_50[k] = -ab_y[k] * dd_1_32[k]
                      + df_1_54[k];
        }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, ab_y, ab_z, dd_1_30, dd_1_33, dd_1_34, \
                         dd_1_35, df_1_52, df_1_56, df_1_57, df_1_58 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_51[k] = -ab_y[k] * dd_1_33[k]
                      + df_1_56[k];

            t_52[k] = -ab_y[k] * dd_1_34[k]
                      + df_1_57[k];

            t_53[k] = -ab_y[k] * dd_1_35[k]
                      + df_1_58[k];

            t_54[k] = -ab_z[k] * dd_1_30[k]
                      + df_1_52[k];
        }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_z, dd_1_31, dd_1_32, dd_1_33, \
                         dd_1_34, dd_1_35, df_1_54, df_1_55, df_1_57, df_1_58, \
                         df_1_59 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_55[k] = -ab_z[k] * dd_1_31[k]
                      + df_1_54[k];

            t_56[k] = -ab_z[k] * dd_1_32[k]
                      + df_1_55[k];

            t_57[k] = -ab_z[k] * dd_1_33[k]
                      + df_1_57[k];

            t_58[k] = -ab_z[k] * dd_1_34[k]
                      + df_1_58[k];

            t_59[k] = -ab_z[k] * dd_1_35[k]
                      + df_1_59[k];
        }
    }
}

}  // namespace simdtrf
