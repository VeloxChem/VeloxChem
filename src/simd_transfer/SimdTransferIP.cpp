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


#include "SimdTransferIP.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_ip_out_of_first(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                            const size_t target, const size_t is, const size_t ks,
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
        auto *t_63 = buffer.data(target + 63 * ncomps + c);
        auto *t_64 = buffer.data(target + 64 * ncomps + c);
        auto *t_65 = buffer.data(target + 65 * ncomps + c);
        auto *t_66 = buffer.data(target + 66 * ncomps + c);
        auto *t_67 = buffer.data(target + 67 * ncomps + c);
        auto *t_68 = buffer.data(target + 68 * ncomps + c);
        auto *t_69 = buffer.data(target + 69 * ncomps + c);
        auto *t_70 = buffer.data(target + 70 * ncomps + c);
        auto *t_71 = buffer.data(target + 71 * ncomps + c);
        auto *t_72 = buffer.data(target + 72 * ncomps + c);
        auto *t_73 = buffer.data(target + 73 * ncomps + c);
        auto *t_74 = buffer.data(target + 74 * ncomps + c);
        auto *t_75 = buffer.data(target + 75 * ncomps + c);
        auto *t_76 = buffer.data(target + 76 * ncomps + c);
        auto *t_77 = buffer.data(target + 77 * ncomps + c);
        auto *t_78 = buffer.data(target + 78 * ncomps + c);
        auto *t_79 = buffer.data(target + 79 * ncomps + c);
        auto *t_80 = buffer.data(target + 80 * ncomps + c);
        auto *t_81 = buffer.data(target + 81 * ncomps + c);
        auto *t_82 = buffer.data(target + 82 * ncomps + c);
        auto *t_83 = buffer.data(target + 83 * ncomps + c);

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *is_0 = buffer.data(is + 0 * ncomps + c);
        const auto *is_1 = buffer.data(is + 1 * ncomps + c);
        const auto *is_2 = buffer.data(is + 2 * ncomps + c);
        const auto *is_3 = buffer.data(is + 3 * ncomps + c);
        const auto *is_4 = buffer.data(is + 4 * ncomps + c);
        const auto *is_5 = buffer.data(is + 5 * ncomps + c);
        const auto *is_6 = buffer.data(is + 6 * ncomps + c);
        const auto *is_7 = buffer.data(is + 7 * ncomps + c);
        const auto *is_8 = buffer.data(is + 8 * ncomps + c);
        const auto *is_9 = buffer.data(is + 9 * ncomps + c);
        const auto *is_10 = buffer.data(is + 10 * ncomps + c);
        const auto *is_11 = buffer.data(is + 11 * ncomps + c);
        const auto *is_12 = buffer.data(is + 12 * ncomps + c);
        const auto *is_13 = buffer.data(is + 13 * ncomps + c);
        const auto *is_14 = buffer.data(is + 14 * ncomps + c);
        const auto *is_15 = buffer.data(is + 15 * ncomps + c);
        const auto *is_16 = buffer.data(is + 16 * ncomps + c);
        const auto *is_17 = buffer.data(is + 17 * ncomps + c);
        const auto *is_18 = buffer.data(is + 18 * ncomps + c);
        const auto *is_19 = buffer.data(is + 19 * ncomps + c);
        const auto *is_20 = buffer.data(is + 20 * ncomps + c);
        const auto *is_21 = buffer.data(is + 21 * ncomps + c);
        const auto *is_22 = buffer.data(is + 22 * ncomps + c);
        const auto *is_23 = buffer.data(is + 23 * ncomps + c);
        const auto *is_24 = buffer.data(is + 24 * ncomps + c);
        const auto *is_25 = buffer.data(is + 25 * ncomps + c);
        const auto *is_26 = buffer.data(is + 26 * ncomps + c);
        const auto *is_27 = buffer.data(is + 27 * ncomps + c);

        const auto *ks_0 = buffer.data(ks + 0 * ncomps + c);
        const auto *ks_1 = buffer.data(ks + 1 * ncomps + c);
        const auto *ks_2 = buffer.data(ks + 2 * ncomps + c);
        const auto *ks_3 = buffer.data(ks + 3 * ncomps + c);
        const auto *ks_4 = buffer.data(ks + 4 * ncomps + c);
        const auto *ks_5 = buffer.data(ks + 5 * ncomps + c);
        const auto *ks_6 = buffer.data(ks + 6 * ncomps + c);
        const auto *ks_7 = buffer.data(ks + 7 * ncomps + c);
        const auto *ks_8 = buffer.data(ks + 8 * ncomps + c);
        const auto *ks_9 = buffer.data(ks + 9 * ncomps + c);
        const auto *ks_10 = buffer.data(ks + 10 * ncomps + c);
        const auto *ks_11 = buffer.data(ks + 11 * ncomps + c);
        const auto *ks_12 = buffer.data(ks + 12 * ncomps + c);
        const auto *ks_13 = buffer.data(ks + 13 * ncomps + c);
        const auto *ks_14 = buffer.data(ks + 14 * ncomps + c);
        const auto *ks_15 = buffer.data(ks + 15 * ncomps + c);
        const auto *ks_16 = buffer.data(ks + 16 * ncomps + c);
        const auto *ks_17 = buffer.data(ks + 17 * ncomps + c);
        const auto *ks_18 = buffer.data(ks + 18 * ncomps + c);
        const auto *ks_19 = buffer.data(ks + 19 * ncomps + c);
        const auto *ks_20 = buffer.data(ks + 20 * ncomps + c);
        const auto *ks_21 = buffer.data(ks + 21 * ncomps + c);
        const auto *ks_22 = buffer.data(ks + 22 * ncomps + c);
        const auto *ks_23 = buffer.data(ks + 23 * ncomps + c);
        const auto *ks_24 = buffer.data(ks + 24 * ncomps + c);
        const auto *ks_25 = buffer.data(ks + 25 * ncomps + c);
        const auto *ks_26 = buffer.data(ks + 26 * ncomps + c);
        const auto *ks_27 = buffer.data(ks + 27 * ncomps + c);
        const auto *ks_28 = buffer.data(ks + 28 * ncomps + c);
        const auto *ks_29 = buffer.data(ks + 29 * ncomps + c);
        const auto *ks_30 = buffer.data(ks + 30 * ncomps + c);
        const auto *ks_31 = buffer.data(ks + 31 * ncomps + c);
        const auto *ks_32 = buffer.data(ks + 32 * ncomps + c);
        const auto *ks_33 = buffer.data(ks + 33 * ncomps + c);
        const auto *ks_34 = buffer.data(ks + 34 * ncomps + c);
        const auto *ks_35 = buffer.data(ks + 35 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, ab_x, ab_y, ab_z, is_0, is_1, ks_0, \
                         ks_1, ks_2, ks_3, ks_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * is_0[k]
                     + ks_0[k];

            t_1[k] = ab_y[k] * is_0[k]
                     + ks_1[k];

            t_2[k] = ab_z[k] * is_0[k]
                     + ks_2[k];

            t_3[k] = ab_x[k] * is_1[k]
                     + ks_1[k];

            t_4[k] = ab_y[k] * is_1[k]
                     + ks_3[k];

            t_5[k] = ab_z[k] * is_1[k]
                     + ks_4[k];
        }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, ab_x, ab_y, ab_z, is_2, is_3, ks_2, ks_3, \
                         ks_4, ks_5, ks_6 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_6[k] = ab_x[k] * is_2[k]
                     + ks_2[k];

            t_7[k] = ab_y[k] * is_2[k]
                     + ks_4[k];

            t_8[k] = ab_z[k] * is_2[k]
                     + ks_5[k];

            t_9[k] = ab_x[k] * is_3[k]
                     + ks_3[k];

            t_10[k] = ab_y[k] * is_3[k]
                      + ks_6[k];
        }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, ab_x, ab_y, ab_z, is_3, is_4, \
                         is_5, ks_4, ks_5, ks_7, ks_8 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_11[k] = ab_z[k] * is_3[k]
                      + ks_7[k];

            t_12[k] = ab_x[k] * is_4[k]
                      + ks_4[k];

            t_13[k] = ab_y[k] * is_4[k]
                      + ks_7[k];

            t_14[k] = ab_z[k] * is_4[k]
                      + ks_8[k];

            t_15[k] = ab_x[k] * is_5[k]
                      + ks_5[k];

            t_16[k] = ab_y[k] * is_5[k]
                      + ks_8[k];
        }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, t_21, ab_x, ab_y, ab_z, is_5, is_6, is_7, \
                         ks_6, ks_7, ks_9, ks_10, ks_11 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_17[k] = ab_z[k] * is_5[k]
                      + ks_9[k];

            t_18[k] = ab_x[k] * is_6[k]
                      + ks_6[k];

            t_19[k] = ab_y[k] * is_6[k]
                      + ks_10[k];

            t_20[k] = ab_z[k] * is_6[k]
                      + ks_11[k];

            t_21[k] = ab_x[k] * is_7[k]
                      + ks_7[k];
        }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, ab_x, ab_y, ab_z, is_7, is_8, ks_8, \
                         ks_11, ks_12, ks_13 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_22[k] = ab_y[k] * is_7[k]
                      + ks_11[k];

            t_23[k] = ab_z[k] * is_7[k]
                      + ks_12[k];

            t_24[k] = ab_x[k] * is_8[k]
                      + ks_8[k];

            t_25[k] = ab_y[k] * is_8[k]
                      + ks_12[k];

            t_26[k] = ab_z[k] * is_8[k]
                      + ks_13[k];
        }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, ab_x, ab_y, ab_z, is_9, is_10, ks_9, \
                         ks_10, ks_13, ks_14, ks_15 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_27[k] = ab_x[k] * is_9[k]
                      + ks_9[k];

            t_28[k] = ab_y[k] * is_9[k]
                      + ks_13[k];

            t_29[k] = ab_z[k] * is_9[k]
                      + ks_14[k];

            t_30[k] = ab_x[k] * is_10[k]
                      + ks_10[k];

            t_31[k] = ab_y[k] * is_10[k]
                      + ks_15[k];
        }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, ab_x, ab_y, ab_z, is_10, is_11, \
                         is_12, ks_11, ks_12, ks_16, ks_17 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_32[k] = ab_z[k] * is_10[k]
                      + ks_16[k];

            t_33[k] = ab_x[k] * is_11[k]
                      + ks_11[k];

            t_34[k] = ab_y[k] * is_11[k]
                      + ks_16[k];

            t_35[k] = ab_z[k] * is_11[k]
                      + ks_17[k];

            t_36[k] = ab_x[k] * is_12[k]
                      + ks_12[k];

            t_37[k] = ab_y[k] * is_12[k]
                      + ks_17[k];
        }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, t_43, ab_x, ab_y, ab_z, is_12, is_13, \
                         is_14, ks_13, ks_14, ks_18, ks_19 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_38[k] = ab_z[k] * is_12[k]
                      + ks_18[k];

            t_39[k] = ab_x[k] * is_13[k]
                      + ks_13[k];

            t_40[k] = ab_y[k] * is_13[k]
                      + ks_18[k];

            t_41[k] = ab_z[k] * is_13[k]
                      + ks_19[k];

            t_42[k] = ab_x[k] * is_14[k]
                      + ks_14[k];

            t_43[k] = ab_y[k] * is_14[k]
                      + ks_19[k];
        }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, ab_x, ab_y, ab_z, is_14, is_15, is_16, \
                         ks_15, ks_16, ks_20, ks_21, ks_22 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_44[k] = ab_z[k] * is_14[k]
                      + ks_20[k];

            t_45[k] = ab_x[k] * is_15[k]
                      + ks_15[k];

            t_46[k] = ab_y[k] * is_15[k]
                      + ks_21[k];

            t_47[k] = ab_z[k] * is_15[k]
                      + ks_22[k];

            t_48[k] = ab_x[k] * is_16[k]
                      + ks_16[k];
        }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, ab_x, ab_y, ab_z, is_16, is_17, ks_17, \
                         ks_22, ks_23, ks_24 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_49[k] = ab_y[k] * is_16[k]
                      + ks_22[k];

            t_50[k] = ab_z[k] * is_16[k]
                      + ks_23[k];

            t_51[k] = ab_x[k] * is_17[k]
                      + ks_17[k];

            t_52[k] = ab_y[k] * is_17[k]
                      + ks_23[k];

            t_53[k] = ab_z[k] * is_17[k]
                      + ks_24[k];
        }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, t_59, ab_x, ab_y, ab_z, is_18, is_19, \
                         ks_18, ks_19, ks_24, ks_25, ks_26 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_54[k] = ab_x[k] * is_18[k]
                      + ks_18[k];

            t_55[k] = ab_y[k] * is_18[k]
                      + ks_24[k];

            t_56[k] = ab_z[k] * is_18[k]
                      + ks_25[k];

            t_57[k] = ab_x[k] * is_19[k]
                      + ks_19[k];

            t_58[k] = ab_y[k] * is_19[k]
                      + ks_25[k];

            t_59[k] = ab_z[k] * is_19[k]
                      + ks_26[k];
        }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, ab_y, ab_z, is_20, is_21, ks_20, \
                         ks_21, ks_26, ks_27, ks_28 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_60[k] = ab_x[k] * is_20[k]
                      + ks_20[k];

            t_61[k] = ab_y[k] * is_20[k]
                      + ks_26[k];

            t_62[k] = ab_z[k] * is_20[k]
                      + ks_27[k];

            t_63[k] = ab_x[k] * is_21[k]
                      + ks_21[k];

            t_64[k] = ab_y[k] * is_21[k]
                      + ks_28[k];
        }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, t_70, ab_x, ab_y, ab_z, is_21, is_22, \
                         is_23, ks_22, ks_23, ks_29, ks_30 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_65[k] = ab_z[k] * is_21[k]
                      + ks_29[k];

            t_66[k] = ab_x[k] * is_22[k]
                      + ks_22[k];

            t_67[k] = ab_y[k] * is_22[k]
                      + ks_29[k];

            t_68[k] = ab_z[k] * is_22[k]
                      + ks_30[k];

            t_69[k] = ab_x[k] * is_23[k]
                      + ks_23[k];

            t_70[k] = ab_y[k] * is_23[k]
                      + ks_30[k];
        }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, t_75, t_76, ab_x, ab_y, ab_z, is_23, is_24, \
                         is_25, ks_24, ks_25, ks_31, ks_32 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_71[k] = ab_z[k] * is_23[k]
                      + ks_31[k];

            t_72[k] = ab_x[k] * is_24[k]
                      + ks_24[k];

            t_73[k] = ab_y[k] * is_24[k]
                      + ks_31[k];

            t_74[k] = ab_z[k] * is_24[k]
                      + ks_32[k];

            t_75[k] = ab_x[k] * is_25[k]
                      + ks_25[k];

            t_76[k] = ab_y[k] * is_25[k]
                      + ks_32[k];
        }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, t_81, t_82, ab_x, ab_y, ab_z, is_25, is_26, \
                         is_27, ks_26, ks_27, ks_33, ks_34 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_77[k] = ab_z[k] * is_25[k]
                      + ks_33[k];

            t_78[k] = ab_x[k] * is_26[k]
                      + ks_26[k];

            t_79[k] = ab_y[k] * is_26[k]
                      + ks_33[k];

            t_80[k] = ab_z[k] * is_26[k]
                      + ks_34[k];

            t_81[k] = ab_x[k] * is_27[k]
                      + ks_27[k];

            t_82[k] = ab_y[k] * is_27[k]
                      + ks_34[k];
        }

#pragma omp simd aligned(t_83, ab_z, is_27, ks_35 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_83[k] = ab_z[k] * is_27[k]
                      + ks_35[k];
        }
    }
}

auto
compute_hrr_ip(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t is, const size_t ks, const size_t ncomps, const size_t nmax) -> void
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
        auto *t_63 = buffer.data(target + 63 * ncomps + c);
        auto *t_64 = buffer.data(target + 64 * ncomps + c);
        auto *t_65 = buffer.data(target + 65 * ncomps + c);
        auto *t_66 = buffer.data(target + 66 * ncomps + c);
        auto *t_67 = buffer.data(target + 67 * ncomps + c);
        auto *t_68 = buffer.data(target + 68 * ncomps + c);
        auto *t_69 = buffer.data(target + 69 * ncomps + c);
        auto *t_70 = buffer.data(target + 70 * ncomps + c);
        auto *t_71 = buffer.data(target + 71 * ncomps + c);
        auto *t_72 = buffer.data(target + 72 * ncomps + c);
        auto *t_73 = buffer.data(target + 73 * ncomps + c);
        auto *t_74 = buffer.data(target + 74 * ncomps + c);
        auto *t_75 = buffer.data(target + 75 * ncomps + c);
        auto *t_76 = buffer.data(target + 76 * ncomps + c);
        auto *t_77 = buffer.data(target + 77 * ncomps + c);
        auto *t_78 = buffer.data(target + 78 * ncomps + c);
        auto *t_79 = buffer.data(target + 79 * ncomps + c);
        auto *t_80 = buffer.data(target + 80 * ncomps + c);
        auto *t_81 = buffer.data(target + 81 * ncomps + c);
        auto *t_82 = buffer.data(target + 82 * ncomps + c);
        auto *t_83 = buffer.data(target + 83 * ncomps + c);

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *is_0 = buffer.data(is + 0 * ncomps + c);
        const auto *is_1 = buffer.data(is + 1 * ncomps + c);
        const auto *is_2 = buffer.data(is + 2 * ncomps + c);
        const auto *is_3 = buffer.data(is + 3 * ncomps + c);
        const auto *is_4 = buffer.data(is + 4 * ncomps + c);
        const auto *is_5 = buffer.data(is + 5 * ncomps + c);
        const auto *is_6 = buffer.data(is + 6 * ncomps + c);
        const auto *is_7 = buffer.data(is + 7 * ncomps + c);
        const auto *is_8 = buffer.data(is + 8 * ncomps + c);
        const auto *is_9 = buffer.data(is + 9 * ncomps + c);
        const auto *is_10 = buffer.data(is + 10 * ncomps + c);
        const auto *is_11 = buffer.data(is + 11 * ncomps + c);
        const auto *is_12 = buffer.data(is + 12 * ncomps + c);
        const auto *is_13 = buffer.data(is + 13 * ncomps + c);
        const auto *is_14 = buffer.data(is + 14 * ncomps + c);
        const auto *is_15 = buffer.data(is + 15 * ncomps + c);
        const auto *is_16 = buffer.data(is + 16 * ncomps + c);
        const auto *is_17 = buffer.data(is + 17 * ncomps + c);
        const auto *is_18 = buffer.data(is + 18 * ncomps + c);
        const auto *is_19 = buffer.data(is + 19 * ncomps + c);
        const auto *is_20 = buffer.data(is + 20 * ncomps + c);
        const auto *is_21 = buffer.data(is + 21 * ncomps + c);
        const auto *is_22 = buffer.data(is + 22 * ncomps + c);
        const auto *is_23 = buffer.data(is + 23 * ncomps + c);
        const auto *is_24 = buffer.data(is + 24 * ncomps + c);
        const auto *is_25 = buffer.data(is + 25 * ncomps + c);
        const auto *is_26 = buffer.data(is + 26 * ncomps + c);
        const auto *is_27 = buffer.data(is + 27 * ncomps + c);

        const auto *ks_0 = buffer.data(ks + 0 * ncomps + c);
        const auto *ks_1 = buffer.data(ks + 1 * ncomps + c);
        const auto *ks_2 = buffer.data(ks + 2 * ncomps + c);
        const auto *ks_3 = buffer.data(ks + 3 * ncomps + c);
        const auto *ks_4 = buffer.data(ks + 4 * ncomps + c);
        const auto *ks_5 = buffer.data(ks + 5 * ncomps + c);
        const auto *ks_6 = buffer.data(ks + 6 * ncomps + c);
        const auto *ks_7 = buffer.data(ks + 7 * ncomps + c);
        const auto *ks_8 = buffer.data(ks + 8 * ncomps + c);
        const auto *ks_9 = buffer.data(ks + 9 * ncomps + c);
        const auto *ks_10 = buffer.data(ks + 10 * ncomps + c);
        const auto *ks_11 = buffer.data(ks + 11 * ncomps + c);
        const auto *ks_12 = buffer.data(ks + 12 * ncomps + c);
        const auto *ks_13 = buffer.data(ks + 13 * ncomps + c);
        const auto *ks_14 = buffer.data(ks + 14 * ncomps + c);
        const auto *ks_15 = buffer.data(ks + 15 * ncomps + c);
        const auto *ks_16 = buffer.data(ks + 16 * ncomps + c);
        const auto *ks_17 = buffer.data(ks + 17 * ncomps + c);
        const auto *ks_18 = buffer.data(ks + 18 * ncomps + c);
        const auto *ks_19 = buffer.data(ks + 19 * ncomps + c);
        const auto *ks_20 = buffer.data(ks + 20 * ncomps + c);
        const auto *ks_21 = buffer.data(ks + 21 * ncomps + c);
        const auto *ks_22 = buffer.data(ks + 22 * ncomps + c);
        const auto *ks_23 = buffer.data(ks + 23 * ncomps + c);
        const auto *ks_24 = buffer.data(ks + 24 * ncomps + c);
        const auto *ks_25 = buffer.data(ks + 25 * ncomps + c);
        const auto *ks_26 = buffer.data(ks + 26 * ncomps + c);
        const auto *ks_27 = buffer.data(ks + 27 * ncomps + c);
        const auto *ks_28 = buffer.data(ks + 28 * ncomps + c);
        const auto *ks_29 = buffer.data(ks + 29 * ncomps + c);
        const auto *ks_30 = buffer.data(ks + 30 * ncomps + c);
        const auto *ks_31 = buffer.data(ks + 31 * ncomps + c);
        const auto *ks_32 = buffer.data(ks + 32 * ncomps + c);
        const auto *ks_33 = buffer.data(ks + 33 * ncomps + c);
        const auto *ks_34 = buffer.data(ks + 34 * ncomps + c);
        const auto *ks_35 = buffer.data(ks + 35 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, ab_x, ab_y, ab_z, is_0, is_1, ks_0, \
                         ks_1, ks_2, ks_3, ks_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * is_0[k]
                     + ks_0[k];

            t_1[k] = ab_y[k] * is_0[k]
                     + ks_1[k];

            t_2[k] = ab_z[k] * is_0[k]
                     + ks_2[k];

            t_3[k] = ab_x[k] * is_1[k]
                     + ks_1[k];

            t_4[k] = ab_y[k] * is_1[k]
                     + ks_3[k];

            t_5[k] = ab_z[k] * is_1[k]
                     + ks_4[k];
        }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, ab_x, ab_y, ab_z, is_2, is_3, ks_2, ks_3, \
                         ks_4, ks_5, ks_6 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_6[k] = ab_x[k] * is_2[k]
                     + ks_2[k];

            t_7[k] = ab_y[k] * is_2[k]
                     + ks_4[k];

            t_8[k] = ab_z[k] * is_2[k]
                     + ks_5[k];

            t_9[k] = ab_x[k] * is_3[k]
                     + ks_3[k];

            t_10[k] = ab_y[k] * is_3[k]
                      + ks_6[k];
        }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, ab_x, ab_y, ab_z, is_3, is_4, \
                         is_5, ks_4, ks_5, ks_7, ks_8 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_11[k] = ab_z[k] * is_3[k]
                      + ks_7[k];

            t_12[k] = ab_x[k] * is_4[k]
                      + ks_4[k];

            t_13[k] = ab_y[k] * is_4[k]
                      + ks_7[k];

            t_14[k] = ab_z[k] * is_4[k]
                      + ks_8[k];

            t_15[k] = ab_x[k] * is_5[k]
                      + ks_5[k];

            t_16[k] = ab_y[k] * is_5[k]
                      + ks_8[k];
        }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, t_21, ab_x, ab_y, ab_z, is_5, is_6, is_7, \
                         ks_6, ks_7, ks_9, ks_10, ks_11 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_17[k] = ab_z[k] * is_5[k]
                      + ks_9[k];

            t_18[k] = ab_x[k] * is_6[k]
                      + ks_6[k];

            t_19[k] = ab_y[k] * is_6[k]
                      + ks_10[k];

            t_20[k] = ab_z[k] * is_6[k]
                      + ks_11[k];

            t_21[k] = ab_x[k] * is_7[k]
                      + ks_7[k];
        }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, ab_x, ab_y, ab_z, is_7, is_8, ks_8, \
                         ks_11, ks_12, ks_13 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_22[k] = ab_y[k] * is_7[k]
                      + ks_11[k];

            t_23[k] = ab_z[k] * is_7[k]
                      + ks_12[k];

            t_24[k] = ab_x[k] * is_8[k]
                      + ks_8[k];

            t_25[k] = ab_y[k] * is_8[k]
                      + ks_12[k];

            t_26[k] = ab_z[k] * is_8[k]
                      + ks_13[k];
        }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, ab_x, ab_y, ab_z, is_9, is_10, ks_9, \
                         ks_10, ks_13, ks_14, ks_15 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_27[k] = ab_x[k] * is_9[k]
                      + ks_9[k];

            t_28[k] = ab_y[k] * is_9[k]
                      + ks_13[k];

            t_29[k] = ab_z[k] * is_9[k]
                      + ks_14[k];

            t_30[k] = ab_x[k] * is_10[k]
                      + ks_10[k];

            t_31[k] = ab_y[k] * is_10[k]
                      + ks_15[k];
        }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, ab_x, ab_y, ab_z, is_10, is_11, \
                         is_12, ks_11, ks_12, ks_16, ks_17 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_32[k] = ab_z[k] * is_10[k]
                      + ks_16[k];

            t_33[k] = ab_x[k] * is_11[k]
                      + ks_11[k];

            t_34[k] = ab_y[k] * is_11[k]
                      + ks_16[k];

            t_35[k] = ab_z[k] * is_11[k]
                      + ks_17[k];

            t_36[k] = ab_x[k] * is_12[k]
                      + ks_12[k];

            t_37[k] = ab_y[k] * is_12[k]
                      + ks_17[k];
        }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, t_43, ab_x, ab_y, ab_z, is_12, is_13, \
                         is_14, ks_13, ks_14, ks_18, ks_19 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_38[k] = ab_z[k] * is_12[k]
                      + ks_18[k];

            t_39[k] = ab_x[k] * is_13[k]
                      + ks_13[k];

            t_40[k] = ab_y[k] * is_13[k]
                      + ks_18[k];

            t_41[k] = ab_z[k] * is_13[k]
                      + ks_19[k];

            t_42[k] = ab_x[k] * is_14[k]
                      + ks_14[k];

            t_43[k] = ab_y[k] * is_14[k]
                      + ks_19[k];
        }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, ab_x, ab_y, ab_z, is_14, is_15, is_16, \
                         ks_15, ks_16, ks_20, ks_21, ks_22 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_44[k] = ab_z[k] * is_14[k]
                      + ks_20[k];

            t_45[k] = ab_x[k] * is_15[k]
                      + ks_15[k];

            t_46[k] = ab_y[k] * is_15[k]
                      + ks_21[k];

            t_47[k] = ab_z[k] * is_15[k]
                      + ks_22[k];

            t_48[k] = ab_x[k] * is_16[k]
                      + ks_16[k];
        }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, ab_x, ab_y, ab_z, is_16, is_17, ks_17, \
                         ks_22, ks_23, ks_24 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_49[k] = ab_y[k] * is_16[k]
                      + ks_22[k];

            t_50[k] = ab_z[k] * is_16[k]
                      + ks_23[k];

            t_51[k] = ab_x[k] * is_17[k]
                      + ks_17[k];

            t_52[k] = ab_y[k] * is_17[k]
                      + ks_23[k];

            t_53[k] = ab_z[k] * is_17[k]
                      + ks_24[k];
        }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, t_59, ab_x, ab_y, ab_z, is_18, is_19, \
                         ks_18, ks_19, ks_24, ks_25, ks_26 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_54[k] = ab_x[k] * is_18[k]
                      + ks_18[k];

            t_55[k] = ab_y[k] * is_18[k]
                      + ks_24[k];

            t_56[k] = ab_z[k] * is_18[k]
                      + ks_25[k];

            t_57[k] = ab_x[k] * is_19[k]
                      + ks_19[k];

            t_58[k] = ab_y[k] * is_19[k]
                      + ks_25[k];

            t_59[k] = ab_z[k] * is_19[k]
                      + ks_26[k];
        }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, ab_y, ab_z, is_20, is_21, ks_20, \
                         ks_21, ks_26, ks_27, ks_28 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_60[k] = ab_x[k] * is_20[k]
                      + ks_20[k];

            t_61[k] = ab_y[k] * is_20[k]
                      + ks_26[k];

            t_62[k] = ab_z[k] * is_20[k]
                      + ks_27[k];

            t_63[k] = ab_x[k] * is_21[k]
                      + ks_21[k];

            t_64[k] = ab_y[k] * is_21[k]
                      + ks_28[k];
        }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, t_70, ab_x, ab_y, ab_z, is_21, is_22, \
                         is_23, ks_22, ks_23, ks_29, ks_30 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_65[k] = ab_z[k] * is_21[k]
                      + ks_29[k];

            t_66[k] = ab_x[k] * is_22[k]
                      + ks_22[k];

            t_67[k] = ab_y[k] * is_22[k]
                      + ks_29[k];

            t_68[k] = ab_z[k] * is_22[k]
                      + ks_30[k];

            t_69[k] = ab_x[k] * is_23[k]
                      + ks_23[k];

            t_70[k] = ab_y[k] * is_23[k]
                      + ks_30[k];
        }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, t_75, t_76, ab_x, ab_y, ab_z, is_23, is_24, \
                         is_25, ks_24, ks_25, ks_31, ks_32 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_71[k] = ab_z[k] * is_23[k]
                      + ks_31[k];

            t_72[k] = ab_x[k] * is_24[k]
                      + ks_24[k];

            t_73[k] = ab_y[k] * is_24[k]
                      + ks_31[k];

            t_74[k] = ab_z[k] * is_24[k]
                      + ks_32[k];

            t_75[k] = ab_x[k] * is_25[k]
                      + ks_25[k];

            t_76[k] = ab_y[k] * is_25[k]
                      + ks_32[k];
        }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, t_81, t_82, ab_x, ab_y, ab_z, is_25, is_26, \
                         is_27, ks_26, ks_27, ks_33, ks_34 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_77[k] = ab_z[k] * is_25[k]
                      + ks_33[k];

            t_78[k] = ab_x[k] * is_26[k]
                      + ks_26[k];

            t_79[k] = ab_y[k] * is_26[k]
                      + ks_33[k];

            t_80[k] = ab_z[k] * is_26[k]
                      + ks_34[k];

            t_81[k] = ab_x[k] * is_27[k]
                      + ks_27[k];

            t_82[k] = ab_y[k] * is_27[k]
                      + ks_34[k];
        }

#pragma omp simd aligned(t_83, ab_z, is_27, ks_35 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_83[k] = ab_z[k] * is_27[k]
                      + ks_35[k];
        }
    }
}

}  // namespace simdtrf
