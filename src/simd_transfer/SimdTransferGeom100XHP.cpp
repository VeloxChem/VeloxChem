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


#include "SimdTransferGeom100XHP.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_geom_100x_hp_out_of_first(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                      const size_t target, const size_t hs_1, const size_t hs_0,
                                      const size_t is_1, const size_t ncomps,
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
        auto *t_60 = buffer.data(target + 60 * ncomps + c);
        auto *t_61 = buffer.data(target + 61 * ncomps + c);
        auto *t_62 = buffer.data(target + 62 * ncomps + c);

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *hs_1_0 = buffer.data(hs_1 + 0 * ncomps + c);
        const auto *hs_1_1 = buffer.data(hs_1 + 1 * ncomps + c);
        const auto *hs_1_2 = buffer.data(hs_1 + 2 * ncomps + c);
        const auto *hs_1_3 = buffer.data(hs_1 + 3 * ncomps + c);
        const auto *hs_1_4 = buffer.data(hs_1 + 4 * ncomps + c);
        const auto *hs_1_5 = buffer.data(hs_1 + 5 * ncomps + c);
        const auto *hs_1_6 = buffer.data(hs_1 + 6 * ncomps + c);
        const auto *hs_1_7 = buffer.data(hs_1 + 7 * ncomps + c);
        const auto *hs_1_8 = buffer.data(hs_1 + 8 * ncomps + c);
        const auto *hs_1_9 = buffer.data(hs_1 + 9 * ncomps + c);
        const auto *hs_1_10 = buffer.data(hs_1 + 10 * ncomps + c);
        const auto *hs_1_11 = buffer.data(hs_1 + 11 * ncomps + c);
        const auto *hs_1_12 = buffer.data(hs_1 + 12 * ncomps + c);
        const auto *hs_1_13 = buffer.data(hs_1 + 13 * ncomps + c);
        const auto *hs_1_14 = buffer.data(hs_1 + 14 * ncomps + c);
        const auto *hs_1_15 = buffer.data(hs_1 + 15 * ncomps + c);
        const auto *hs_1_16 = buffer.data(hs_1 + 16 * ncomps + c);
        const auto *hs_1_17 = buffer.data(hs_1 + 17 * ncomps + c);
        const auto *hs_1_18 = buffer.data(hs_1 + 18 * ncomps + c);
        const auto *hs_1_19 = buffer.data(hs_1 + 19 * ncomps + c);
        const auto *hs_1_20 = buffer.data(hs_1 + 20 * ncomps + c);

        const auto *hs_0_0 = buffer.data(hs_0 + 0 * ncomps + c);
        const auto *hs_0_1 = buffer.data(hs_0 + 1 * ncomps + c);
        const auto *hs_0_2 = buffer.data(hs_0 + 2 * ncomps + c);
        const auto *hs_0_3 = buffer.data(hs_0 + 3 * ncomps + c);
        const auto *hs_0_4 = buffer.data(hs_0 + 4 * ncomps + c);
        const auto *hs_0_5 = buffer.data(hs_0 + 5 * ncomps + c);
        const auto *hs_0_6 = buffer.data(hs_0 + 6 * ncomps + c);
        const auto *hs_0_7 = buffer.data(hs_0 + 7 * ncomps + c);
        const auto *hs_0_8 = buffer.data(hs_0 + 8 * ncomps + c);
        const auto *hs_0_9 = buffer.data(hs_0 + 9 * ncomps + c);
        const auto *hs_0_10 = buffer.data(hs_0 + 10 * ncomps + c);
        const auto *hs_0_11 = buffer.data(hs_0 + 11 * ncomps + c);
        const auto *hs_0_12 = buffer.data(hs_0 + 12 * ncomps + c);
        const auto *hs_0_13 = buffer.data(hs_0 + 13 * ncomps + c);
        const auto *hs_0_14 = buffer.data(hs_0 + 14 * ncomps + c);
        const auto *hs_0_15 = buffer.data(hs_0 + 15 * ncomps + c);
        const auto *hs_0_16 = buffer.data(hs_0 + 16 * ncomps + c);
        const auto *hs_0_17 = buffer.data(hs_0 + 17 * ncomps + c);
        const auto *hs_0_18 = buffer.data(hs_0 + 18 * ncomps + c);
        const auto *hs_0_19 = buffer.data(hs_0 + 19 * ncomps + c);
        const auto *hs_0_20 = buffer.data(hs_0 + 20 * ncomps + c);

        const auto *is_1_0 = buffer.data(is_1 + 0 * ncomps + c);
        const auto *is_1_1 = buffer.data(is_1 + 1 * ncomps + c);
        const auto *is_1_2 = buffer.data(is_1 + 2 * ncomps + c);
        const auto *is_1_3 = buffer.data(is_1 + 3 * ncomps + c);
        const auto *is_1_4 = buffer.data(is_1 + 4 * ncomps + c);
        const auto *is_1_5 = buffer.data(is_1 + 5 * ncomps + c);
        const auto *is_1_6 = buffer.data(is_1 + 6 * ncomps + c);
        const auto *is_1_7 = buffer.data(is_1 + 7 * ncomps + c);
        const auto *is_1_8 = buffer.data(is_1 + 8 * ncomps + c);
        const auto *is_1_9 = buffer.data(is_1 + 9 * ncomps + c);
        const auto *is_1_10 = buffer.data(is_1 + 10 * ncomps + c);
        const auto *is_1_11 = buffer.data(is_1 + 11 * ncomps + c);
        const auto *is_1_12 = buffer.data(is_1 + 12 * ncomps + c);
        const auto *is_1_13 = buffer.data(is_1 + 13 * ncomps + c);
        const auto *is_1_14 = buffer.data(is_1 + 14 * ncomps + c);
        const auto *is_1_15 = buffer.data(is_1 + 15 * ncomps + c);
        const auto *is_1_16 = buffer.data(is_1 + 16 * ncomps + c);
        const auto *is_1_17 = buffer.data(is_1 + 17 * ncomps + c);
        const auto *is_1_18 = buffer.data(is_1 + 18 * ncomps + c);
        const auto *is_1_19 = buffer.data(is_1 + 19 * ncomps + c);
        const auto *is_1_20 = buffer.data(is_1 + 20 * ncomps + c);
        const auto *is_1_21 = buffer.data(is_1 + 21 * ncomps + c);
        const auto *is_1_22 = buffer.data(is_1 + 22 * ncomps + c);
        const auto *is_1_23 = buffer.data(is_1 + 23 * ncomps + c);
        const auto *is_1_24 = buffer.data(is_1 + 24 * ncomps + c);
        const auto *is_1_25 = buffer.data(is_1 + 25 * ncomps + c);
        const auto *is_1_26 = buffer.data(is_1 + 26 * ncomps + c);
        const auto *is_1_27 = buffer.data(is_1 + 27 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, ab_y, ab_z, hs_1_0, hs_1_1, hs_0_0, \
                         hs_0_1, is_1_0, is_1_1, is_1_2, is_1_3 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * hs_1_0[k]
                     + hs_0_0[k]
                     + is_1_0[k];

            t_1[k] = ab_y[k] * hs_1_0[k]
                     + is_1_1[k];

            t_2[k] = ab_z[k] * hs_1_0[k]
                     + is_1_2[k];

            t_3[k] = ab_x[k] * hs_1_1[k]
                     + hs_0_1[k]
                     + is_1_1[k];

            t_4[k] = ab_y[k] * hs_1_1[k]
                     + is_1_3[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, ab_x, ab_y, ab_z, hs_1_1, hs_1_2, hs_0_2, is_1_2, \
                         is_1_4, is_1_5 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = ab_z[k] * hs_1_1[k]
                     + is_1_4[k];

            t_6[k] = ab_x[k] * hs_1_2[k]
                     + hs_0_2[k]
                     + is_1_2[k];

            t_7[k] = ab_y[k] * hs_1_2[k]
                     + is_1_4[k];

            t_8[k] = ab_z[k] * hs_1_2[k]
                     + is_1_5[k];
        }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, t_13, ab_x, ab_y, ab_z, hs_1_3, hs_1_4, \
                         hs_0_3, hs_0_4, is_1_3, is_1_4, is_1_6, \
                         is_1_7 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_9[k] = ab_x[k] * hs_1_3[k]
                     + hs_0_3[k]
                     + is_1_3[k];

            t_10[k] = ab_y[k] * hs_1_3[k]
                      + is_1_6[k];

            t_11[k] = ab_z[k] * hs_1_3[k]
                      + is_1_7[k];

            t_12[k] = ab_x[k] * hs_1_4[k]
                      + hs_0_4[k]
                      + is_1_4[k];

            t_13[k] = ab_y[k] * hs_1_4[k]
                      + is_1_7[k];
        }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, ab_x, ab_y, ab_z, hs_1_4, hs_1_5, hs_0_5, \
                         is_1_5, is_1_8, is_1_9 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_14[k] = ab_z[k] * hs_1_4[k]
                      + is_1_8[k];

            t_15[k] = ab_x[k] * hs_1_5[k]
                      + hs_0_5[k]
                      + is_1_5[k];

            t_16[k] = ab_y[k] * hs_1_5[k]
                      + is_1_8[k];

            t_17[k] = ab_z[k] * hs_1_5[k]
                      + is_1_9[k];
        }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, ab_x, ab_y, ab_z, hs_1_6, hs_1_7, \
                         hs_0_6, hs_0_7, is_1_6, is_1_7, is_1_10, \
                         is_1_11 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_18[k] = ab_x[k] * hs_1_6[k]
                      + hs_0_6[k]
                      + is_1_6[k];

            t_19[k] = ab_y[k] * hs_1_6[k]
                      + is_1_10[k];

            t_20[k] = ab_z[k] * hs_1_6[k]
                      + is_1_11[k];

            t_21[k] = ab_x[k] * hs_1_7[k]
                      + hs_0_7[k]
                      + is_1_7[k];

            t_22[k] = ab_y[k] * hs_1_7[k]
                      + is_1_11[k];
        }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, ab_x, ab_y, ab_z, hs_1_7, hs_1_8, hs_0_8, \
                         is_1_8, is_1_12, is_1_13 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_23[k] = ab_z[k] * hs_1_7[k]
                      + is_1_12[k];

            t_24[k] = ab_x[k] * hs_1_8[k]
                      + hs_0_8[k]
                      + is_1_8[k];

            t_25[k] = ab_y[k] * hs_1_8[k]
                      + is_1_12[k];

            t_26[k] = ab_z[k] * hs_1_8[k]
                      + is_1_13[k];
        }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, ab_x, ab_y, ab_z, hs_1_9, hs_1_10, hs_0_9, \
                         hs_0_10, is_1_9, is_1_10, is_1_13, is_1_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_27[k] = ab_x[k] * hs_1_9[k]
                      + hs_0_9[k]
                      + is_1_9[k];

            t_28[k] = ab_y[k] * hs_1_9[k]
                      + is_1_13[k];

            t_29[k] = ab_z[k] * hs_1_9[k]
                      + is_1_14[k];

            t_30[k] = ab_x[k] * hs_1_10[k]
                      + hs_0_10[k]
                      + is_1_10[k];
        }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, ab_x, ab_y, ab_z, hs_1_10, hs_1_11, \
                         hs_0_11, is_1_11, is_1_15, is_1_16, is_1_17 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_31[k] = ab_y[k] * hs_1_10[k]
                      + is_1_15[k];

            t_32[k] = ab_z[k] * hs_1_10[k]
                      + is_1_16[k];

            t_33[k] = ab_x[k] * hs_1_11[k]
                      + hs_0_11[k]
                      + is_1_11[k];

            t_34[k] = ab_y[k] * hs_1_11[k]
                      + is_1_16[k];

            t_35[k] = ab_z[k] * hs_1_11[k]
                      + is_1_17[k];
        }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, ab_x, ab_y, ab_z, hs_1_12, hs_1_13, \
                         hs_0_12, hs_0_13, is_1_12, is_1_13, is_1_17, \
                         is_1_18 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_36[k] = ab_x[k] * hs_1_12[k]
                      + hs_0_12[k]
                      + is_1_12[k];

            t_37[k] = ab_y[k] * hs_1_12[k]
                      + is_1_17[k];

            t_38[k] = ab_z[k] * hs_1_12[k]
                      + is_1_18[k];

            t_39[k] = ab_x[k] * hs_1_13[k]
                      + hs_0_13[k]
                      + is_1_13[k];

            t_40[k] = ab_y[k] * hs_1_13[k]
                      + is_1_18[k];
        }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, ab_x, ab_y, ab_z, hs_1_13, hs_1_14, hs_0_14, \
                         is_1_14, is_1_19, is_1_20 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_41[k] = ab_z[k] * hs_1_13[k]
                      + is_1_19[k];

            t_42[k] = ab_x[k] * hs_1_14[k]
                      + hs_0_14[k]
                      + is_1_14[k];

            t_43[k] = ab_y[k] * hs_1_14[k]
                      + is_1_19[k];

            t_44[k] = ab_z[k] * hs_1_14[k]
                      + is_1_20[k];
        }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, ab_y, ab_z, hs_1_15, hs_1_16, \
                         hs_0_15, hs_0_16, is_1_15, is_1_16, is_1_21, \
                         is_1_22 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_45[k] = ab_x[k] * hs_1_15[k]
                      + hs_0_15[k]
                      + is_1_15[k];

            t_46[k] = ab_y[k] * hs_1_15[k]
                      + is_1_21[k];

            t_47[k] = ab_z[k] * hs_1_15[k]
                      + is_1_22[k];

            t_48[k] = ab_x[k] * hs_1_16[k]
                      + hs_0_16[k]
                      + is_1_16[k];

            t_49[k] = ab_y[k] * hs_1_16[k]
                      + is_1_22[k];
        }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, ab_x, ab_y, ab_z, hs_1_16, hs_1_17, hs_0_17, \
                         is_1_17, is_1_23, is_1_24 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_50[k] = ab_z[k] * hs_1_16[k]
                      + is_1_23[k];

            t_51[k] = ab_x[k] * hs_1_17[k]
                      + hs_0_17[k]
                      + is_1_17[k];

            t_52[k] = ab_y[k] * hs_1_17[k]
                      + is_1_23[k];

            t_53[k] = ab_z[k] * hs_1_17[k]
                      + is_1_24[k];
        }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, ab_x, ab_y, ab_z, hs_1_18, hs_1_19, \
                         hs_0_18, hs_0_19, is_1_18, is_1_19, is_1_24, \
                         is_1_25 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_54[k] = ab_x[k] * hs_1_18[k]
                      + hs_0_18[k]
                      + is_1_18[k];

            t_55[k] = ab_y[k] * hs_1_18[k]
                      + is_1_24[k];

            t_56[k] = ab_z[k] * hs_1_18[k]
                      + is_1_25[k];

            t_57[k] = ab_x[k] * hs_1_19[k]
                      + hs_0_19[k]
                      + is_1_19[k];

            t_58[k] = ab_y[k] * hs_1_19[k]
                      + is_1_25[k];
        }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, ab_x, ab_y, ab_z, hs_1_19, hs_1_20, hs_0_20, \
                         is_1_20, is_1_26, is_1_27 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_59[k] = ab_z[k] * hs_1_19[k]
                      + is_1_26[k];

            t_60[k] = ab_x[k] * hs_1_20[k]
                      + hs_0_20[k]
                      + is_1_20[k];

            t_61[k] = ab_y[k] * hs_1_20[k]
                      + is_1_26[k];

            t_62[k] = ab_z[k] * hs_1_20[k]
                      + is_1_27[k];
        }
    }
}

}  // namespace simdtrf
