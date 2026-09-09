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


#include "SimdTransferGD.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_gd_out_of_first(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                            const size_t target, const size_t gp, const size_t hp,
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
        auto *t_84 = buffer.data(target + 84 * ncomps + c);
        auto *t_85 = buffer.data(target + 85 * ncomps + c);
        auto *t_86 = buffer.data(target + 86 * ncomps + c);
        auto *t_87 = buffer.data(target + 87 * ncomps + c);
        auto *t_88 = buffer.data(target + 88 * ncomps + c);
        auto *t_89 = buffer.data(target + 89 * ncomps + c);

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *gp_0 = buffer.data(gp + 0 * ncomps + c);
        const auto *gp_1 = buffer.data(gp + 1 * ncomps + c);
        const auto *gp_2 = buffer.data(gp + 2 * ncomps + c);
        const auto *gp_3 = buffer.data(gp + 3 * ncomps + c);
        const auto *gp_4 = buffer.data(gp + 4 * ncomps + c);
        const auto *gp_5 = buffer.data(gp + 5 * ncomps + c);
        const auto *gp_6 = buffer.data(gp + 6 * ncomps + c);
        const auto *gp_7 = buffer.data(gp + 7 * ncomps + c);
        const auto *gp_8 = buffer.data(gp + 8 * ncomps + c);
        const auto *gp_9 = buffer.data(gp + 9 * ncomps + c);
        const auto *gp_10 = buffer.data(gp + 10 * ncomps + c);
        const auto *gp_11 = buffer.data(gp + 11 * ncomps + c);
        const auto *gp_12 = buffer.data(gp + 12 * ncomps + c);
        const auto *gp_13 = buffer.data(gp + 13 * ncomps + c);
        const auto *gp_14 = buffer.data(gp + 14 * ncomps + c);
        const auto *gp_15 = buffer.data(gp + 15 * ncomps + c);
        const auto *gp_16 = buffer.data(gp + 16 * ncomps + c);
        const auto *gp_17 = buffer.data(gp + 17 * ncomps + c);
        const auto *gp_18 = buffer.data(gp + 18 * ncomps + c);
        const auto *gp_19 = buffer.data(gp + 19 * ncomps + c);
        const auto *gp_20 = buffer.data(gp + 20 * ncomps + c);
        const auto *gp_21 = buffer.data(gp + 21 * ncomps + c);
        const auto *gp_22 = buffer.data(gp + 22 * ncomps + c);
        const auto *gp_23 = buffer.data(gp + 23 * ncomps + c);
        const auto *gp_24 = buffer.data(gp + 24 * ncomps + c);
        const auto *gp_25 = buffer.data(gp + 25 * ncomps + c);
        const auto *gp_26 = buffer.data(gp + 26 * ncomps + c);
        const auto *gp_27 = buffer.data(gp + 27 * ncomps + c);
        const auto *gp_28 = buffer.data(gp + 28 * ncomps + c);
        const auto *gp_29 = buffer.data(gp + 29 * ncomps + c);
        const auto *gp_30 = buffer.data(gp + 30 * ncomps + c);
        const auto *gp_31 = buffer.data(gp + 31 * ncomps + c);
        const auto *gp_32 = buffer.data(gp + 32 * ncomps + c);
        const auto *gp_33 = buffer.data(gp + 33 * ncomps + c);
        const auto *gp_34 = buffer.data(gp + 34 * ncomps + c);
        const auto *gp_35 = buffer.data(gp + 35 * ncomps + c);
        const auto *gp_36 = buffer.data(gp + 36 * ncomps + c);
        const auto *gp_37 = buffer.data(gp + 37 * ncomps + c);
        const auto *gp_38 = buffer.data(gp + 38 * ncomps + c);
        const auto *gp_39 = buffer.data(gp + 39 * ncomps + c);
        const auto *gp_40 = buffer.data(gp + 40 * ncomps + c);
        const auto *gp_41 = buffer.data(gp + 41 * ncomps + c);
        const auto *gp_42 = buffer.data(gp + 42 * ncomps + c);
        const auto *gp_43 = buffer.data(gp + 43 * ncomps + c);
        const auto *gp_44 = buffer.data(gp + 44 * ncomps + c);

        const auto *hp_0 = buffer.data(hp + 0 * ncomps + c);
        const auto *hp_1 = buffer.data(hp + 1 * ncomps + c);
        const auto *hp_2 = buffer.data(hp + 2 * ncomps + c);
        const auto *hp_3 = buffer.data(hp + 3 * ncomps + c);
        const auto *hp_4 = buffer.data(hp + 4 * ncomps + c);
        const auto *hp_5 = buffer.data(hp + 5 * ncomps + c);
        const auto *hp_6 = buffer.data(hp + 6 * ncomps + c);
        const auto *hp_7 = buffer.data(hp + 7 * ncomps + c);
        const auto *hp_8 = buffer.data(hp + 8 * ncomps + c);
        const auto *hp_9 = buffer.data(hp + 9 * ncomps + c);
        const auto *hp_10 = buffer.data(hp + 10 * ncomps + c);
        const auto *hp_11 = buffer.data(hp + 11 * ncomps + c);
        const auto *hp_12 = buffer.data(hp + 12 * ncomps + c);
        const auto *hp_13 = buffer.data(hp + 13 * ncomps + c);
        const auto *hp_14 = buffer.data(hp + 14 * ncomps + c);
        const auto *hp_15 = buffer.data(hp + 15 * ncomps + c);
        const auto *hp_16 = buffer.data(hp + 16 * ncomps + c);
        const auto *hp_17 = buffer.data(hp + 17 * ncomps + c);
        const auto *hp_18 = buffer.data(hp + 18 * ncomps + c);
        const auto *hp_19 = buffer.data(hp + 19 * ncomps + c);
        const auto *hp_20 = buffer.data(hp + 20 * ncomps + c);
        const auto *hp_21 = buffer.data(hp + 21 * ncomps + c);
        const auto *hp_22 = buffer.data(hp + 22 * ncomps + c);
        const auto *hp_23 = buffer.data(hp + 23 * ncomps + c);
        const auto *hp_24 = buffer.data(hp + 24 * ncomps + c);
        const auto *hp_25 = buffer.data(hp + 25 * ncomps + c);
        const auto *hp_26 = buffer.data(hp + 26 * ncomps + c);
        const auto *hp_27 = buffer.data(hp + 27 * ncomps + c);
        const auto *hp_28 = buffer.data(hp + 28 * ncomps + c);
        const auto *hp_29 = buffer.data(hp + 29 * ncomps + c);
        const auto *hp_30 = buffer.data(hp + 30 * ncomps + c);
        const auto *hp_31 = buffer.data(hp + 31 * ncomps + c);
        const auto *hp_32 = buffer.data(hp + 32 * ncomps + c);
        const auto *hp_33 = buffer.data(hp + 33 * ncomps + c);
        const auto *hp_34 = buffer.data(hp + 34 * ncomps + c);
        const auto *hp_35 = buffer.data(hp + 35 * ncomps + c);
        const auto *hp_36 = buffer.data(hp + 36 * ncomps + c);
        const auto *hp_37 = buffer.data(hp + 37 * ncomps + c);
        const auto *hp_38 = buffer.data(hp + 38 * ncomps + c);
        const auto *hp_39 = buffer.data(hp + 39 * ncomps + c);
        const auto *hp_40 = buffer.data(hp + 40 * ncomps + c);
        const auto *hp_41 = buffer.data(hp + 41 * ncomps + c);
        const auto *hp_42 = buffer.data(hp + 42 * ncomps + c);
        const auto *hp_43 = buffer.data(hp + 43 * ncomps + c);
        const auto *hp_44 = buffer.data(hp + 44 * ncomps + c);
        const auto *hp_46 = buffer.data(hp + 46 * ncomps + c);
        const auto *hp_47 = buffer.data(hp + 47 * ncomps + c);
        const auto *hp_49 = buffer.data(hp + 49 * ncomps + c);
        const auto *hp_50 = buffer.data(hp + 50 * ncomps + c);
        const auto *hp_52 = buffer.data(hp + 52 * ncomps + c);
        const auto *hp_53 = buffer.data(hp + 53 * ncomps + c);
        const auto *hp_55 = buffer.data(hp + 55 * ncomps + c);
        const auto *hp_56 = buffer.data(hp + 56 * ncomps + c);
        const auto *hp_58 = buffer.data(hp + 58 * ncomps + c);
        const auto *hp_59 = buffer.data(hp + 59 * ncomps + c);
        const auto *hp_62 = buffer.data(hp + 62 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, ab_y, gp_0, gp_1, gp_2, hp_0, hp_1, \
                         hp_2, hp_4, hp_5 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * gp_0[k]
                     + hp_0[k];

            t_1[k] = ab_x[k] * gp_1[k]
                     + hp_1[k];

            t_2[k] = ab_x[k] * gp_2[k]
                     + hp_2[k];

            t_3[k] = ab_y[k] * gp_1[k]
                     + hp_4[k];

            t_4[k] = ab_y[k] * gp_2[k]
                     + hp_5[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, ab_x, ab_z, gp_2, gp_3, gp_4, gp_5, hp_3, hp_4, \
                         hp_5, hp_8 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = ab_z[k] * gp_2[k]
                     + hp_8[k];

            t_6[k] = ab_x[k] * gp_3[k]
                     + hp_3[k];

            t_7[k] = ab_x[k] * gp_4[k]
                     + hp_4[k];

            t_8[k] = ab_x[k] * gp_5[k]
                     + hp_5[k];
        }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, ab_x, ab_y, ab_z, gp_4, gp_5, gp_6, hp_6, \
                         hp_10, hp_11, hp_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_9[k] = ab_y[k] * gp_4[k]
                     + hp_10[k];

            t_10[k] = ab_y[k] * gp_5[k]
                      + hp_11[k];

            t_11[k] = ab_z[k] * gp_5[k]
                      + hp_14[k];

            t_12[k] = ab_x[k] * gp_6[k]
                      + hp_6[k];
        }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, ab_x, ab_y, ab_z, gp_7, gp_8, hp_7, \
                         hp_8, hp_13, hp_14, hp_17 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_13[k] = ab_x[k] * gp_7[k]
                      + hp_7[k];

            t_14[k] = ab_x[k] * gp_8[k]
                      + hp_8[k];

            t_15[k] = ab_y[k] * gp_7[k]
                      + hp_13[k];

            t_16[k] = ab_y[k] * gp_8[k]
                      + hp_14[k];

            t_17[k] = ab_z[k] * gp_8[k]
                      + hp_17[k];
        }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, ab_x, ab_y, gp_9, gp_10, gp_11, hp_9, \
                         hp_10, hp_11, hp_19, hp_20 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_18[k] = ab_x[k] * gp_9[k]
                      + hp_9[k];

            t_19[k] = ab_x[k] * gp_10[k]
                      + hp_10[k];

            t_20[k] = ab_x[k] * gp_11[k]
                      + hp_11[k];

            t_21[k] = ab_y[k] * gp_10[k]
                      + hp_19[k];

            t_22[k] = ab_y[k] * gp_11[k]
                      + hp_20[k];
        }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, ab_x, ab_z, gp_11, gp_12, gp_13, gp_14, \
                         hp_12, hp_13, hp_14, hp_23 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_23[k] = ab_z[k] * gp_11[k]
                      + hp_23[k];

            t_24[k] = ab_x[k] * gp_12[k]
                      + hp_12[k];

            t_25[k] = ab_x[k] * gp_13[k]
                      + hp_13[k];

            t_26[k] = ab_x[k] * gp_14[k]
                      + hp_14[k];
        }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, ab_x, ab_y, ab_z, gp_13, gp_14, gp_15, hp_15, \
                         hp_22, hp_23, hp_26 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_27[k] = ab_y[k] * gp_13[k]
                      + hp_22[k];

            t_28[k] = ab_y[k] * gp_14[k]
                      + hp_23[k];

            t_29[k] = ab_z[k] * gp_14[k]
                      + hp_26[k];

            t_30[k] = ab_x[k] * gp_15[k]
                      + hp_15[k];
        }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, ab_x, ab_y, ab_z, gp_16, gp_17, hp_16, \
                         hp_17, hp_25, hp_26, hp_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_31[k] = ab_x[k] * gp_16[k]
                      + hp_16[k];

            t_32[k] = ab_x[k] * gp_17[k]
                      + hp_17[k];

            t_33[k] = ab_y[k] * gp_16[k]
                      + hp_25[k];

            t_34[k] = ab_y[k] * gp_17[k]
                      + hp_26[k];

            t_35[k] = ab_z[k] * gp_17[k]
                      + hp_29[k];
        }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, ab_x, ab_y, gp_18, gp_19, gp_20, hp_18, \
                         hp_19, hp_20, hp_31, hp_32 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_36[k] = ab_x[k] * gp_18[k]
                      + hp_18[k];

            t_37[k] = ab_x[k] * gp_19[k]
                      + hp_19[k];

            t_38[k] = ab_x[k] * gp_20[k]
                      + hp_20[k];

            t_39[k] = ab_y[k] * gp_19[k]
                      + hp_31[k];

            t_40[k] = ab_y[k] * gp_20[k]
                      + hp_32[k];
        }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, ab_x, ab_z, gp_20, gp_21, gp_22, gp_23, \
                         hp_21, hp_22, hp_23, hp_35 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_41[k] = ab_z[k] * gp_20[k]
                      + hp_35[k];

            t_42[k] = ab_x[k] * gp_21[k]
                      + hp_21[k];

            t_43[k] = ab_x[k] * gp_22[k]
                      + hp_22[k];

            t_44[k] = ab_x[k] * gp_23[k]
                      + hp_23[k];
        }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, ab_x, ab_y, ab_z, gp_22, gp_23, gp_24, hp_24, \
                         hp_34, hp_35, hp_38 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_45[k] = ab_y[k] * gp_22[k]
                      + hp_34[k];

            t_46[k] = ab_y[k] * gp_23[k]
                      + hp_35[k];

            t_47[k] = ab_z[k] * gp_23[k]
                      + hp_38[k];

            t_48[k] = ab_x[k] * gp_24[k]
                      + hp_24[k];
        }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, ab_x, ab_y, ab_z, gp_25, gp_26, hp_25, \
                         hp_26, hp_37, hp_38, hp_41 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_49[k] = ab_x[k] * gp_25[k]
                      + hp_25[k];

            t_50[k] = ab_x[k] * gp_26[k]
                      + hp_26[k];

            t_51[k] = ab_y[k] * gp_25[k]
                      + hp_37[k];

            t_52[k] = ab_y[k] * gp_26[k]
                      + hp_38[k];

            t_53[k] = ab_z[k] * gp_26[k]
                      + hp_41[k];
        }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, ab_x, ab_y, gp_27, gp_28, gp_29, hp_27, \
                         hp_28, hp_29, hp_40, hp_41 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_54[k] = ab_x[k] * gp_27[k]
                      + hp_27[k];

            t_55[k] = ab_x[k] * gp_28[k]
                      + hp_28[k];

            t_56[k] = ab_x[k] * gp_29[k]
                      + hp_29[k];

            t_57[k] = ab_y[k] * gp_28[k]
                      + hp_40[k];

            t_58[k] = ab_y[k] * gp_29[k]
                      + hp_41[k];
        }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, ab_x, ab_z, gp_29, gp_30, gp_31, gp_32, \
                         hp_30, hp_31, hp_32, hp_44 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_59[k] = ab_z[k] * gp_29[k]
                      + hp_44[k];

            t_60[k] = ab_x[k] * gp_30[k]
                      + hp_30[k];

            t_61[k] = ab_x[k] * gp_31[k]
                      + hp_31[k];

            t_62[k] = ab_x[k] * gp_32[k]
                      + hp_32[k];
        }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, ab_x, ab_y, ab_z, gp_31, gp_32, gp_33, hp_33, \
                         hp_46, hp_47, hp_50 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_63[k] = ab_y[k] * gp_31[k]
                      + hp_46[k];

            t_64[k] = ab_y[k] * gp_32[k]
                      + hp_47[k];

            t_65[k] = ab_z[k] * gp_32[k]
                      + hp_50[k];

            t_66[k] = ab_x[k] * gp_33[k]
                      + hp_33[k];
        }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, t_71, ab_x, ab_y, ab_z, gp_34, gp_35, hp_34, \
                         hp_35, hp_49, hp_50, hp_53 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_67[k] = ab_x[k] * gp_34[k]
                      + hp_34[k];

            t_68[k] = ab_x[k] * gp_35[k]
                      + hp_35[k];

            t_69[k] = ab_y[k] * gp_34[k]
                      + hp_49[k];

            t_70[k] = ab_y[k] * gp_35[k]
                      + hp_50[k];

            t_71[k] = ab_z[k] * gp_35[k]
                      + hp_53[k];
        }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, ab_x, ab_y, gp_36, gp_37, gp_38, hp_36, \
                         hp_37, hp_38, hp_52, hp_53 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_72[k] = ab_x[k] * gp_36[k]
                      + hp_36[k];

            t_73[k] = ab_x[k] * gp_37[k]
                      + hp_37[k];

            t_74[k] = ab_x[k] * gp_38[k]
                      + hp_38[k];

            t_75[k] = ab_y[k] * gp_37[k]
                      + hp_52[k];

            t_76[k] = ab_y[k] * gp_38[k]
                      + hp_53[k];
        }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, ab_x, ab_z, gp_38, gp_39, gp_40, gp_41, \
                         hp_39, hp_40, hp_41, hp_56 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_77[k] = ab_z[k] * gp_38[k]
                      + hp_56[k];

            t_78[k] = ab_x[k] * gp_39[k]
                      + hp_39[k];

            t_79[k] = ab_x[k] * gp_40[k]
                      + hp_40[k];

            t_80[k] = ab_x[k] * gp_41[k]
                      + hp_41[k];
        }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, ab_x, ab_y, ab_z, gp_40, gp_41, gp_42, hp_42, \
                         hp_55, hp_56, hp_59 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_81[k] = ab_y[k] * gp_40[k]
                      + hp_55[k];

            t_82[k] = ab_y[k] * gp_41[k]
                      + hp_56[k];

            t_83[k] = ab_z[k] * gp_41[k]
                      + hp_59[k];

            t_84[k] = ab_x[k] * gp_42[k]
                      + hp_42[k];
        }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_x, ab_y, ab_z, gp_43, gp_44, hp_43, \
                         hp_44, hp_58, hp_59, hp_62 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_85[k] = ab_x[k] * gp_43[k]
                      + hp_43[k];

            t_86[k] = ab_x[k] * gp_44[k]
                      + hp_44[k];

            t_87[k] = ab_y[k] * gp_43[k]
                      + hp_58[k];

            t_88[k] = ab_y[k] * gp_44[k]
                      + hp_59[k];

            t_89[k] = ab_z[k] * gp_44[k]
                      + hp_62[k];
        }
    }
}

auto
compute_hrr_gd(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t gp, const size_t hp, const size_t ncomps, const size_t nmax) -> void
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
        auto *t_84 = buffer.data(target + 84 * ncomps + c);
        auto *t_85 = buffer.data(target + 85 * ncomps + c);
        auto *t_86 = buffer.data(target + 86 * ncomps + c);
        auto *t_87 = buffer.data(target + 87 * ncomps + c);
        auto *t_88 = buffer.data(target + 88 * ncomps + c);
        auto *t_89 = buffer.data(target + 89 * ncomps + c);

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *gp_0 = buffer.data(gp + 0 * ncomps + c);
        const auto *gp_1 = buffer.data(gp + 1 * ncomps + c);
        const auto *gp_2 = buffer.data(gp + 2 * ncomps + c);
        const auto *gp_3 = buffer.data(gp + 3 * ncomps + c);
        const auto *gp_4 = buffer.data(gp + 4 * ncomps + c);
        const auto *gp_5 = buffer.data(gp + 5 * ncomps + c);
        const auto *gp_6 = buffer.data(gp + 6 * ncomps + c);
        const auto *gp_7 = buffer.data(gp + 7 * ncomps + c);
        const auto *gp_8 = buffer.data(gp + 8 * ncomps + c);
        const auto *gp_9 = buffer.data(gp + 9 * ncomps + c);
        const auto *gp_10 = buffer.data(gp + 10 * ncomps + c);
        const auto *gp_11 = buffer.data(gp + 11 * ncomps + c);
        const auto *gp_12 = buffer.data(gp + 12 * ncomps + c);
        const auto *gp_13 = buffer.data(gp + 13 * ncomps + c);
        const auto *gp_14 = buffer.data(gp + 14 * ncomps + c);
        const auto *gp_15 = buffer.data(gp + 15 * ncomps + c);
        const auto *gp_16 = buffer.data(gp + 16 * ncomps + c);
        const auto *gp_17 = buffer.data(gp + 17 * ncomps + c);
        const auto *gp_18 = buffer.data(gp + 18 * ncomps + c);
        const auto *gp_19 = buffer.data(gp + 19 * ncomps + c);
        const auto *gp_20 = buffer.data(gp + 20 * ncomps + c);
        const auto *gp_21 = buffer.data(gp + 21 * ncomps + c);
        const auto *gp_22 = buffer.data(gp + 22 * ncomps + c);
        const auto *gp_23 = buffer.data(gp + 23 * ncomps + c);
        const auto *gp_24 = buffer.data(gp + 24 * ncomps + c);
        const auto *gp_25 = buffer.data(gp + 25 * ncomps + c);
        const auto *gp_26 = buffer.data(gp + 26 * ncomps + c);
        const auto *gp_27 = buffer.data(gp + 27 * ncomps + c);
        const auto *gp_28 = buffer.data(gp + 28 * ncomps + c);
        const auto *gp_29 = buffer.data(gp + 29 * ncomps + c);
        const auto *gp_30 = buffer.data(gp + 30 * ncomps + c);
        const auto *gp_31 = buffer.data(gp + 31 * ncomps + c);
        const auto *gp_32 = buffer.data(gp + 32 * ncomps + c);
        const auto *gp_33 = buffer.data(gp + 33 * ncomps + c);
        const auto *gp_34 = buffer.data(gp + 34 * ncomps + c);
        const auto *gp_35 = buffer.data(gp + 35 * ncomps + c);
        const auto *gp_36 = buffer.data(gp + 36 * ncomps + c);
        const auto *gp_37 = buffer.data(gp + 37 * ncomps + c);
        const auto *gp_38 = buffer.data(gp + 38 * ncomps + c);
        const auto *gp_39 = buffer.data(gp + 39 * ncomps + c);
        const auto *gp_40 = buffer.data(gp + 40 * ncomps + c);
        const auto *gp_41 = buffer.data(gp + 41 * ncomps + c);
        const auto *gp_42 = buffer.data(gp + 42 * ncomps + c);
        const auto *gp_43 = buffer.data(gp + 43 * ncomps + c);
        const auto *gp_44 = buffer.data(gp + 44 * ncomps + c);

        const auto *hp_0 = buffer.data(hp + 0 * ncomps + c);
        const auto *hp_1 = buffer.data(hp + 1 * ncomps + c);
        const auto *hp_2 = buffer.data(hp + 2 * ncomps + c);
        const auto *hp_3 = buffer.data(hp + 3 * ncomps + c);
        const auto *hp_4 = buffer.data(hp + 4 * ncomps + c);
        const auto *hp_5 = buffer.data(hp + 5 * ncomps + c);
        const auto *hp_6 = buffer.data(hp + 6 * ncomps + c);
        const auto *hp_7 = buffer.data(hp + 7 * ncomps + c);
        const auto *hp_8 = buffer.data(hp + 8 * ncomps + c);
        const auto *hp_9 = buffer.data(hp + 9 * ncomps + c);
        const auto *hp_10 = buffer.data(hp + 10 * ncomps + c);
        const auto *hp_11 = buffer.data(hp + 11 * ncomps + c);
        const auto *hp_12 = buffer.data(hp + 12 * ncomps + c);
        const auto *hp_13 = buffer.data(hp + 13 * ncomps + c);
        const auto *hp_14 = buffer.data(hp + 14 * ncomps + c);
        const auto *hp_15 = buffer.data(hp + 15 * ncomps + c);
        const auto *hp_16 = buffer.data(hp + 16 * ncomps + c);
        const auto *hp_17 = buffer.data(hp + 17 * ncomps + c);
        const auto *hp_18 = buffer.data(hp + 18 * ncomps + c);
        const auto *hp_19 = buffer.data(hp + 19 * ncomps + c);
        const auto *hp_20 = buffer.data(hp + 20 * ncomps + c);
        const auto *hp_21 = buffer.data(hp + 21 * ncomps + c);
        const auto *hp_22 = buffer.data(hp + 22 * ncomps + c);
        const auto *hp_23 = buffer.data(hp + 23 * ncomps + c);
        const auto *hp_24 = buffer.data(hp + 24 * ncomps + c);
        const auto *hp_25 = buffer.data(hp + 25 * ncomps + c);
        const auto *hp_26 = buffer.data(hp + 26 * ncomps + c);
        const auto *hp_27 = buffer.data(hp + 27 * ncomps + c);
        const auto *hp_28 = buffer.data(hp + 28 * ncomps + c);
        const auto *hp_29 = buffer.data(hp + 29 * ncomps + c);
        const auto *hp_30 = buffer.data(hp + 30 * ncomps + c);
        const auto *hp_31 = buffer.data(hp + 31 * ncomps + c);
        const auto *hp_32 = buffer.data(hp + 32 * ncomps + c);
        const auto *hp_33 = buffer.data(hp + 33 * ncomps + c);
        const auto *hp_34 = buffer.data(hp + 34 * ncomps + c);
        const auto *hp_35 = buffer.data(hp + 35 * ncomps + c);
        const auto *hp_36 = buffer.data(hp + 36 * ncomps + c);
        const auto *hp_37 = buffer.data(hp + 37 * ncomps + c);
        const auto *hp_38 = buffer.data(hp + 38 * ncomps + c);
        const auto *hp_39 = buffer.data(hp + 39 * ncomps + c);
        const auto *hp_40 = buffer.data(hp + 40 * ncomps + c);
        const auto *hp_41 = buffer.data(hp + 41 * ncomps + c);
        const auto *hp_42 = buffer.data(hp + 42 * ncomps + c);
        const auto *hp_43 = buffer.data(hp + 43 * ncomps + c);
        const auto *hp_44 = buffer.data(hp + 44 * ncomps + c);
        const auto *hp_46 = buffer.data(hp + 46 * ncomps + c);
        const auto *hp_47 = buffer.data(hp + 47 * ncomps + c);
        const auto *hp_49 = buffer.data(hp + 49 * ncomps + c);
        const auto *hp_50 = buffer.data(hp + 50 * ncomps + c);
        const auto *hp_52 = buffer.data(hp + 52 * ncomps + c);
        const auto *hp_53 = buffer.data(hp + 53 * ncomps + c);
        const auto *hp_55 = buffer.data(hp + 55 * ncomps + c);
        const auto *hp_56 = buffer.data(hp + 56 * ncomps + c);
        const auto *hp_58 = buffer.data(hp + 58 * ncomps + c);
        const auto *hp_59 = buffer.data(hp + 59 * ncomps + c);
        const auto *hp_62 = buffer.data(hp + 62 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, ab_y, gp_0, gp_1, gp_2, hp_0, hp_1, \
                         hp_2, hp_4, hp_5 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * gp_0[k]
                     + hp_0[k];

            t_1[k] = ab_x[k] * gp_1[k]
                     + hp_1[k];

            t_2[k] = ab_x[k] * gp_2[k]
                     + hp_2[k];

            t_3[k] = ab_y[k] * gp_1[k]
                     + hp_4[k];

            t_4[k] = ab_y[k] * gp_2[k]
                     + hp_5[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, ab_x, ab_z, gp_2, gp_3, gp_4, gp_5, hp_3, hp_4, \
                         hp_5, hp_8 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = ab_z[k] * gp_2[k]
                     + hp_8[k];

            t_6[k] = ab_x[k] * gp_3[k]
                     + hp_3[k];

            t_7[k] = ab_x[k] * gp_4[k]
                     + hp_4[k];

            t_8[k] = ab_x[k] * gp_5[k]
                     + hp_5[k];
        }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, ab_x, ab_y, ab_z, gp_4, gp_5, gp_6, hp_6, \
                         hp_10, hp_11, hp_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_9[k] = ab_y[k] * gp_4[k]
                     + hp_10[k];

            t_10[k] = ab_y[k] * gp_5[k]
                      + hp_11[k];

            t_11[k] = ab_z[k] * gp_5[k]
                      + hp_14[k];

            t_12[k] = ab_x[k] * gp_6[k]
                      + hp_6[k];
        }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, ab_x, ab_y, ab_z, gp_7, gp_8, hp_7, \
                         hp_8, hp_13, hp_14, hp_17 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_13[k] = ab_x[k] * gp_7[k]
                      + hp_7[k];

            t_14[k] = ab_x[k] * gp_8[k]
                      + hp_8[k];

            t_15[k] = ab_y[k] * gp_7[k]
                      + hp_13[k];

            t_16[k] = ab_y[k] * gp_8[k]
                      + hp_14[k];

            t_17[k] = ab_z[k] * gp_8[k]
                      + hp_17[k];
        }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, ab_x, ab_y, gp_9, gp_10, gp_11, hp_9, \
                         hp_10, hp_11, hp_19, hp_20 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_18[k] = ab_x[k] * gp_9[k]
                      + hp_9[k];

            t_19[k] = ab_x[k] * gp_10[k]
                      + hp_10[k];

            t_20[k] = ab_x[k] * gp_11[k]
                      + hp_11[k];

            t_21[k] = ab_y[k] * gp_10[k]
                      + hp_19[k];

            t_22[k] = ab_y[k] * gp_11[k]
                      + hp_20[k];
        }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, ab_x, ab_z, gp_11, gp_12, gp_13, gp_14, \
                         hp_12, hp_13, hp_14, hp_23 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_23[k] = ab_z[k] * gp_11[k]
                      + hp_23[k];

            t_24[k] = ab_x[k] * gp_12[k]
                      + hp_12[k];

            t_25[k] = ab_x[k] * gp_13[k]
                      + hp_13[k];

            t_26[k] = ab_x[k] * gp_14[k]
                      + hp_14[k];
        }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, ab_x, ab_y, ab_z, gp_13, gp_14, gp_15, hp_15, \
                         hp_22, hp_23, hp_26 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_27[k] = ab_y[k] * gp_13[k]
                      + hp_22[k];

            t_28[k] = ab_y[k] * gp_14[k]
                      + hp_23[k];

            t_29[k] = ab_z[k] * gp_14[k]
                      + hp_26[k];

            t_30[k] = ab_x[k] * gp_15[k]
                      + hp_15[k];
        }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, ab_x, ab_y, ab_z, gp_16, gp_17, hp_16, \
                         hp_17, hp_25, hp_26, hp_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_31[k] = ab_x[k] * gp_16[k]
                      + hp_16[k];

            t_32[k] = ab_x[k] * gp_17[k]
                      + hp_17[k];

            t_33[k] = ab_y[k] * gp_16[k]
                      + hp_25[k];

            t_34[k] = ab_y[k] * gp_17[k]
                      + hp_26[k];

            t_35[k] = ab_z[k] * gp_17[k]
                      + hp_29[k];
        }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, ab_x, ab_y, gp_18, gp_19, gp_20, hp_18, \
                         hp_19, hp_20, hp_31, hp_32 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_36[k] = ab_x[k] * gp_18[k]
                      + hp_18[k];

            t_37[k] = ab_x[k] * gp_19[k]
                      + hp_19[k];

            t_38[k] = ab_x[k] * gp_20[k]
                      + hp_20[k];

            t_39[k] = ab_y[k] * gp_19[k]
                      + hp_31[k];

            t_40[k] = ab_y[k] * gp_20[k]
                      + hp_32[k];
        }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, ab_x, ab_z, gp_20, gp_21, gp_22, gp_23, \
                         hp_21, hp_22, hp_23, hp_35 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_41[k] = ab_z[k] * gp_20[k]
                      + hp_35[k];

            t_42[k] = ab_x[k] * gp_21[k]
                      + hp_21[k];

            t_43[k] = ab_x[k] * gp_22[k]
                      + hp_22[k];

            t_44[k] = ab_x[k] * gp_23[k]
                      + hp_23[k];
        }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, ab_x, ab_y, ab_z, gp_22, gp_23, gp_24, hp_24, \
                         hp_34, hp_35, hp_38 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_45[k] = ab_y[k] * gp_22[k]
                      + hp_34[k];

            t_46[k] = ab_y[k] * gp_23[k]
                      + hp_35[k];

            t_47[k] = ab_z[k] * gp_23[k]
                      + hp_38[k];

            t_48[k] = ab_x[k] * gp_24[k]
                      + hp_24[k];
        }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, ab_x, ab_y, ab_z, gp_25, gp_26, hp_25, \
                         hp_26, hp_37, hp_38, hp_41 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_49[k] = ab_x[k] * gp_25[k]
                      + hp_25[k];

            t_50[k] = ab_x[k] * gp_26[k]
                      + hp_26[k];

            t_51[k] = ab_y[k] * gp_25[k]
                      + hp_37[k];

            t_52[k] = ab_y[k] * gp_26[k]
                      + hp_38[k];

            t_53[k] = ab_z[k] * gp_26[k]
                      + hp_41[k];
        }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, ab_x, ab_y, gp_27, gp_28, gp_29, hp_27, \
                         hp_28, hp_29, hp_40, hp_41 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_54[k] = ab_x[k] * gp_27[k]
                      + hp_27[k];

            t_55[k] = ab_x[k] * gp_28[k]
                      + hp_28[k];

            t_56[k] = ab_x[k] * gp_29[k]
                      + hp_29[k];

            t_57[k] = ab_y[k] * gp_28[k]
                      + hp_40[k];

            t_58[k] = ab_y[k] * gp_29[k]
                      + hp_41[k];
        }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, ab_x, ab_z, gp_29, gp_30, gp_31, gp_32, \
                         hp_30, hp_31, hp_32, hp_44 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_59[k] = ab_z[k] * gp_29[k]
                      + hp_44[k];

            t_60[k] = ab_x[k] * gp_30[k]
                      + hp_30[k];

            t_61[k] = ab_x[k] * gp_31[k]
                      + hp_31[k];

            t_62[k] = ab_x[k] * gp_32[k]
                      + hp_32[k];
        }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, ab_x, ab_y, ab_z, gp_31, gp_32, gp_33, hp_33, \
                         hp_46, hp_47, hp_50 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_63[k] = ab_y[k] * gp_31[k]
                      + hp_46[k];

            t_64[k] = ab_y[k] * gp_32[k]
                      + hp_47[k];

            t_65[k] = ab_z[k] * gp_32[k]
                      + hp_50[k];

            t_66[k] = ab_x[k] * gp_33[k]
                      + hp_33[k];
        }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, t_71, ab_x, ab_y, ab_z, gp_34, gp_35, hp_34, \
                         hp_35, hp_49, hp_50, hp_53 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_67[k] = ab_x[k] * gp_34[k]
                      + hp_34[k];

            t_68[k] = ab_x[k] * gp_35[k]
                      + hp_35[k];

            t_69[k] = ab_y[k] * gp_34[k]
                      + hp_49[k];

            t_70[k] = ab_y[k] * gp_35[k]
                      + hp_50[k];

            t_71[k] = ab_z[k] * gp_35[k]
                      + hp_53[k];
        }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, ab_x, ab_y, gp_36, gp_37, gp_38, hp_36, \
                         hp_37, hp_38, hp_52, hp_53 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_72[k] = ab_x[k] * gp_36[k]
                      + hp_36[k];

            t_73[k] = ab_x[k] * gp_37[k]
                      + hp_37[k];

            t_74[k] = ab_x[k] * gp_38[k]
                      + hp_38[k];

            t_75[k] = ab_y[k] * gp_37[k]
                      + hp_52[k];

            t_76[k] = ab_y[k] * gp_38[k]
                      + hp_53[k];
        }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, ab_x, ab_z, gp_38, gp_39, gp_40, gp_41, \
                         hp_39, hp_40, hp_41, hp_56 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_77[k] = ab_z[k] * gp_38[k]
                      + hp_56[k];

            t_78[k] = ab_x[k] * gp_39[k]
                      + hp_39[k];

            t_79[k] = ab_x[k] * gp_40[k]
                      + hp_40[k];

            t_80[k] = ab_x[k] * gp_41[k]
                      + hp_41[k];
        }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, ab_x, ab_y, ab_z, gp_40, gp_41, gp_42, hp_42, \
                         hp_55, hp_56, hp_59 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_81[k] = ab_y[k] * gp_40[k]
                      + hp_55[k];

            t_82[k] = ab_y[k] * gp_41[k]
                      + hp_56[k];

            t_83[k] = ab_z[k] * gp_41[k]
                      + hp_59[k];

            t_84[k] = ab_x[k] * gp_42[k]
                      + hp_42[k];
        }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_x, ab_y, ab_z, gp_43, gp_44, hp_43, \
                         hp_44, hp_58, hp_59, hp_62 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_85[k] = ab_x[k] * gp_43[k]
                      + hp_43[k];

            t_86[k] = ab_x[k] * gp_44[k]
                      + hp_44[k];

            t_87[k] = ab_y[k] * gp_43[k]
                      + hp_58[k];

            t_88[k] = ab_y[k] * gp_44[k]
                      + hp_59[k];

            t_89[k] = ab_z[k] * gp_44[k]
                      + hp_62[k];
        }
    }
}

}  // namespace simdtrf
