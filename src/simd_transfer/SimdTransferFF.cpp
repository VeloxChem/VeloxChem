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


#include "SimdTransferFF.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_ff_out_of_first(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                            const size_t target, const size_t fd, const size_t gd,
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
        auto *t_90 = buffer.data(target + 90 * ncomps + c);
        auto *t_91 = buffer.data(target + 91 * ncomps + c);
        auto *t_92 = buffer.data(target + 92 * ncomps + c);
        auto *t_93 = buffer.data(target + 93 * ncomps + c);
        auto *t_94 = buffer.data(target + 94 * ncomps + c);
        auto *t_95 = buffer.data(target + 95 * ncomps + c);
        auto *t_96 = buffer.data(target + 96 * ncomps + c);
        auto *t_97 = buffer.data(target + 97 * ncomps + c);
        auto *t_98 = buffer.data(target + 98 * ncomps + c);
        auto *t_99 = buffer.data(target + 99 * ncomps + c);

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *fd_0 = buffer.data(fd + 0 * ncomps + c);
        const auto *fd_1 = buffer.data(fd + 1 * ncomps + c);
        const auto *fd_2 = buffer.data(fd + 2 * ncomps + c);
        const auto *fd_3 = buffer.data(fd + 3 * ncomps + c);
        const auto *fd_4 = buffer.data(fd + 4 * ncomps + c);
        const auto *fd_5 = buffer.data(fd + 5 * ncomps + c);
        const auto *fd_6 = buffer.data(fd + 6 * ncomps + c);
        const auto *fd_7 = buffer.data(fd + 7 * ncomps + c);
        const auto *fd_8 = buffer.data(fd + 8 * ncomps + c);
        const auto *fd_9 = buffer.data(fd + 9 * ncomps + c);
        const auto *fd_10 = buffer.data(fd + 10 * ncomps + c);
        const auto *fd_11 = buffer.data(fd + 11 * ncomps + c);
        const auto *fd_12 = buffer.data(fd + 12 * ncomps + c);
        const auto *fd_13 = buffer.data(fd + 13 * ncomps + c);
        const auto *fd_14 = buffer.data(fd + 14 * ncomps + c);
        const auto *fd_15 = buffer.data(fd + 15 * ncomps + c);
        const auto *fd_16 = buffer.data(fd + 16 * ncomps + c);
        const auto *fd_17 = buffer.data(fd + 17 * ncomps + c);
        const auto *fd_18 = buffer.data(fd + 18 * ncomps + c);
        const auto *fd_19 = buffer.data(fd + 19 * ncomps + c);
        const auto *fd_20 = buffer.data(fd + 20 * ncomps + c);
        const auto *fd_21 = buffer.data(fd + 21 * ncomps + c);
        const auto *fd_22 = buffer.data(fd + 22 * ncomps + c);
        const auto *fd_23 = buffer.data(fd + 23 * ncomps + c);
        const auto *fd_24 = buffer.data(fd + 24 * ncomps + c);
        const auto *fd_25 = buffer.data(fd + 25 * ncomps + c);
        const auto *fd_26 = buffer.data(fd + 26 * ncomps + c);
        const auto *fd_27 = buffer.data(fd + 27 * ncomps + c);
        const auto *fd_28 = buffer.data(fd + 28 * ncomps + c);
        const auto *fd_29 = buffer.data(fd + 29 * ncomps + c);
        const auto *fd_30 = buffer.data(fd + 30 * ncomps + c);
        const auto *fd_31 = buffer.data(fd + 31 * ncomps + c);
        const auto *fd_32 = buffer.data(fd + 32 * ncomps + c);
        const auto *fd_33 = buffer.data(fd + 33 * ncomps + c);
        const auto *fd_34 = buffer.data(fd + 34 * ncomps + c);
        const auto *fd_35 = buffer.data(fd + 35 * ncomps + c);
        const auto *fd_36 = buffer.data(fd + 36 * ncomps + c);
        const auto *fd_37 = buffer.data(fd + 37 * ncomps + c);
        const auto *fd_38 = buffer.data(fd + 38 * ncomps + c);
        const auto *fd_39 = buffer.data(fd + 39 * ncomps + c);
        const auto *fd_40 = buffer.data(fd + 40 * ncomps + c);
        const auto *fd_41 = buffer.data(fd + 41 * ncomps + c);
        const auto *fd_42 = buffer.data(fd + 42 * ncomps + c);
        const auto *fd_43 = buffer.data(fd + 43 * ncomps + c);
        const auto *fd_44 = buffer.data(fd + 44 * ncomps + c);
        const auto *fd_45 = buffer.data(fd + 45 * ncomps + c);
        const auto *fd_46 = buffer.data(fd + 46 * ncomps + c);
        const auto *fd_47 = buffer.data(fd + 47 * ncomps + c);
        const auto *fd_48 = buffer.data(fd + 48 * ncomps + c);
        const auto *fd_49 = buffer.data(fd + 49 * ncomps + c);
        const auto *fd_50 = buffer.data(fd + 50 * ncomps + c);
        const auto *fd_51 = buffer.data(fd + 51 * ncomps + c);
        const auto *fd_52 = buffer.data(fd + 52 * ncomps + c);
        const auto *fd_53 = buffer.data(fd + 53 * ncomps + c);
        const auto *fd_54 = buffer.data(fd + 54 * ncomps + c);
        const auto *fd_55 = buffer.data(fd + 55 * ncomps + c);
        const auto *fd_56 = buffer.data(fd + 56 * ncomps + c);
        const auto *fd_57 = buffer.data(fd + 57 * ncomps + c);
        const auto *fd_58 = buffer.data(fd + 58 * ncomps + c);
        const auto *fd_59 = buffer.data(fd + 59 * ncomps + c);

        const auto *gd_0 = buffer.data(gd + 0 * ncomps + c);
        const auto *gd_1 = buffer.data(gd + 1 * ncomps + c);
        const auto *gd_2 = buffer.data(gd + 2 * ncomps + c);
        const auto *gd_3 = buffer.data(gd + 3 * ncomps + c);
        const auto *gd_4 = buffer.data(gd + 4 * ncomps + c);
        const auto *gd_5 = buffer.data(gd + 5 * ncomps + c);
        const auto *gd_6 = buffer.data(gd + 6 * ncomps + c);
        const auto *gd_7 = buffer.data(gd + 7 * ncomps + c);
        const auto *gd_8 = buffer.data(gd + 8 * ncomps + c);
        const auto *gd_9 = buffer.data(gd + 9 * ncomps + c);
        const auto *gd_10 = buffer.data(gd + 10 * ncomps + c);
        const auto *gd_11 = buffer.data(gd + 11 * ncomps + c);
        const auto *gd_12 = buffer.data(gd + 12 * ncomps + c);
        const auto *gd_13 = buffer.data(gd + 13 * ncomps + c);
        const auto *gd_14 = buffer.data(gd + 14 * ncomps + c);
        const auto *gd_15 = buffer.data(gd + 15 * ncomps + c);
        const auto *gd_16 = buffer.data(gd + 16 * ncomps + c);
        const auto *gd_17 = buffer.data(gd + 17 * ncomps + c);
        const auto *gd_18 = buffer.data(gd + 18 * ncomps + c);
        const auto *gd_19 = buffer.data(gd + 19 * ncomps + c);
        const auto *gd_20 = buffer.data(gd + 20 * ncomps + c);
        const auto *gd_21 = buffer.data(gd + 21 * ncomps + c);
        const auto *gd_22 = buffer.data(gd + 22 * ncomps + c);
        const auto *gd_23 = buffer.data(gd + 23 * ncomps + c);
        const auto *gd_24 = buffer.data(gd + 24 * ncomps + c);
        const auto *gd_25 = buffer.data(gd + 25 * ncomps + c);
        const auto *gd_26 = buffer.data(gd + 26 * ncomps + c);
        const auto *gd_27 = buffer.data(gd + 27 * ncomps + c);
        const auto *gd_28 = buffer.data(gd + 28 * ncomps + c);
        const auto *gd_29 = buffer.data(gd + 29 * ncomps + c);
        const auto *gd_30 = buffer.data(gd + 30 * ncomps + c);
        const auto *gd_31 = buffer.data(gd + 31 * ncomps + c);
        const auto *gd_32 = buffer.data(gd + 32 * ncomps + c);
        const auto *gd_33 = buffer.data(gd + 33 * ncomps + c);
        const auto *gd_34 = buffer.data(gd + 34 * ncomps + c);
        const auto *gd_35 = buffer.data(gd + 35 * ncomps + c);
        const auto *gd_36 = buffer.data(gd + 36 * ncomps + c);
        const auto *gd_37 = buffer.data(gd + 37 * ncomps + c);
        const auto *gd_38 = buffer.data(gd + 38 * ncomps + c);
        const auto *gd_39 = buffer.data(gd + 39 * ncomps + c);
        const auto *gd_40 = buffer.data(gd + 40 * ncomps + c);
        const auto *gd_41 = buffer.data(gd + 41 * ncomps + c);
        const auto *gd_42 = buffer.data(gd + 42 * ncomps + c);
        const auto *gd_43 = buffer.data(gd + 43 * ncomps + c);
        const auto *gd_44 = buffer.data(gd + 44 * ncomps + c);
        const auto *gd_45 = buffer.data(gd + 45 * ncomps + c);
        const auto *gd_46 = buffer.data(gd + 46 * ncomps + c);
        const auto *gd_47 = buffer.data(gd + 47 * ncomps + c);
        const auto *gd_48 = buffer.data(gd + 48 * ncomps + c);
        const auto *gd_49 = buffer.data(gd + 49 * ncomps + c);
        const auto *gd_50 = buffer.data(gd + 50 * ncomps + c);
        const auto *gd_51 = buffer.data(gd + 51 * ncomps + c);
        const auto *gd_52 = buffer.data(gd + 52 * ncomps + c);
        const auto *gd_53 = buffer.data(gd + 53 * ncomps + c);
        const auto *gd_54 = buffer.data(gd + 54 * ncomps + c);
        const auto *gd_55 = buffer.data(gd + 55 * ncomps + c);
        const auto *gd_56 = buffer.data(gd + 56 * ncomps + c);
        const auto *gd_57 = buffer.data(gd + 57 * ncomps + c);
        const auto *gd_58 = buffer.data(gd + 58 * ncomps + c);
        const auto *gd_59 = buffer.data(gd + 59 * ncomps + c);
        const auto *gd_63 = buffer.data(gd + 63 * ncomps + c);
        const auto *gd_64 = buffer.data(gd + 64 * ncomps + c);
        const auto *gd_65 = buffer.data(gd + 65 * ncomps + c);
        const auto *gd_69 = buffer.data(gd + 69 * ncomps + c);
        const auto *gd_70 = buffer.data(gd + 70 * ncomps + c);
        const auto *gd_71 = buffer.data(gd + 71 * ncomps + c);
        const auto *gd_75 = buffer.data(gd + 75 * ncomps + c);
        const auto *gd_76 = buffer.data(gd + 76 * ncomps + c);
        const auto *gd_77 = buffer.data(gd + 77 * ncomps + c);
        const auto *gd_81 = buffer.data(gd + 81 * ncomps + c);
        const auto *gd_82 = buffer.data(gd + 82 * ncomps + c);
        const auto *gd_83 = buffer.data(gd + 83 * ncomps + c);
        const auto *gd_89 = buffer.data(gd + 89 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, fd_0, fd_1, fd_2, fd_3, fd_4, gd_0, \
                         gd_1, gd_2, gd_3, gd_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * fd_0[k]
                     + gd_0[k];

            t_1[k] = ab_x[k] * fd_1[k]
                     + gd_1[k];

            t_2[k] = ab_x[k] * fd_2[k]
                     + gd_2[k];

            t_3[k] = ab_x[k] * fd_3[k]
                     + gd_3[k];

            t_4[k] = ab_x[k] * fd_4[k]
                     + gd_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, ab_y, ab_z, fd_3, fd_4, fd_5, gd_5, \
                         gd_9, gd_10, gd_11, gd_17 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = ab_x[k] * fd_5[k]
                     + gd_5[k];

            t_6[k] = ab_y[k] * fd_3[k]
                     + gd_9[k];

            t_7[k] = ab_y[k] * fd_4[k]
                     + gd_10[k];

            t_8[k] = ab_y[k] * fd_5[k]
                     + gd_11[k];

            t_9[k] = ab_z[k] * fd_5[k]
                     + gd_17[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, fd_6, fd_7, fd_8, fd_9, fd_10, \
                         gd_6, gd_7, gd_8, gd_9, gd_10 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = ab_x[k] * fd_6[k]
                      + gd_6[k];

            t_11[k] = ab_x[k] * fd_7[k]
                      + gd_7[k];

            t_12[k] = ab_x[k] * fd_8[k]
                      + gd_8[k];

            t_13[k] = ab_x[k] * fd_9[k]
                      + gd_9[k];

            t_14[k] = ab_x[k] * fd_10[k]
                      + gd_10[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, ab_y, ab_z, fd_9, fd_10, fd_11, \
                         gd_11, gd_21, gd_22, gd_23, gd_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = ab_x[k] * fd_11[k]
                      + gd_11[k];

            t_16[k] = ab_y[k] * fd_9[k]
                      + gd_21[k];

            t_17[k] = ab_y[k] * fd_10[k]
                      + gd_22[k];

            t_18[k] = ab_y[k] * fd_11[k]
                      + gd_23[k];

            t_19[k] = ab_z[k] * fd_11[k]
                      + gd_29[k];
        }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, fd_12, fd_13, fd_14, fd_15, \
                         fd_16, gd_12, gd_13, gd_14, gd_15, gd_16 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_20[k] = ab_x[k] * fd_12[k]
                      + gd_12[k];

            t_21[k] = ab_x[k] * fd_13[k]
                      + gd_13[k];

            t_22[k] = ab_x[k] * fd_14[k]
                      + gd_14[k];

            t_23[k] = ab_x[k] * fd_15[k]
                      + gd_15[k];

            t_24[k] = ab_x[k] * fd_16[k]
                      + gd_16[k];
        }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, ab_y, ab_z, fd_15, fd_16, fd_17, \
                         gd_17, gd_27, gd_28, gd_29, gd_35 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_25[k] = ab_x[k] * fd_17[k]
                      + gd_17[k];

            t_26[k] = ab_y[k] * fd_15[k]
                      + gd_27[k];

            t_27[k] = ab_y[k] * fd_16[k]
                      + gd_28[k];

            t_28[k] = ab_y[k] * fd_17[k]
                      + gd_29[k];

            t_29[k] = ab_z[k] * fd_17[k]
                      + gd_35[k];
        }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, fd_18, fd_19, fd_20, fd_21, \
                         fd_22, gd_18, gd_19, gd_20, gd_21, gd_22 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_30[k] = ab_x[k] * fd_18[k]
                      + gd_18[k];

            t_31[k] = ab_x[k] * fd_19[k]
                      + gd_19[k];

            t_32[k] = ab_x[k] * fd_20[k]
                      + gd_20[k];

            t_33[k] = ab_x[k] * fd_21[k]
                      + gd_21[k];

            t_34[k] = ab_x[k] * fd_22[k]
                      + gd_22[k];
        }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, ab_y, ab_z, fd_21, fd_22, fd_23, \
                         gd_23, gd_39, gd_40, gd_41, gd_47 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_35[k] = ab_x[k] * fd_23[k]
                      + gd_23[k];

            t_36[k] = ab_y[k] * fd_21[k]
                      + gd_39[k];

            t_37[k] = ab_y[k] * fd_22[k]
                      + gd_40[k];

            t_38[k] = ab_y[k] * fd_23[k]
                      + gd_41[k];

            t_39[k] = ab_z[k] * fd_23[k]
                      + gd_47[k];
        }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_x, fd_24, fd_25, fd_26, fd_27, \
                         fd_28, gd_24, gd_25, gd_26, gd_27, gd_28 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_40[k] = ab_x[k] * fd_24[k]
                      + gd_24[k];

            t_41[k] = ab_x[k] * fd_25[k]
                      + gd_25[k];

            t_42[k] = ab_x[k] * fd_26[k]
                      + gd_26[k];

            t_43[k] = ab_x[k] * fd_27[k]
                      + gd_27[k];

            t_44[k] = ab_x[k] * fd_28[k]
                      + gd_28[k];
        }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, ab_y, ab_z, fd_27, fd_28, fd_29, \
                         gd_29, gd_45, gd_46, gd_47, gd_53 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_45[k] = ab_x[k] * fd_29[k]
                      + gd_29[k];

            t_46[k] = ab_y[k] * fd_27[k]
                      + gd_45[k];

            t_47[k] = ab_y[k] * fd_28[k]
                      + gd_46[k];

            t_48[k] = ab_y[k] * fd_29[k]
                      + gd_47[k];

            t_49[k] = ab_z[k] * fd_29[k]
                      + gd_53[k];
        }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, fd_30, fd_31, fd_32, fd_33, \
                         fd_34, gd_30, gd_31, gd_32, gd_33, gd_34 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_50[k] = ab_x[k] * fd_30[k]
                      + gd_30[k];

            t_51[k] = ab_x[k] * fd_31[k]
                      + gd_31[k];

            t_52[k] = ab_x[k] * fd_32[k]
                      + gd_32[k];

            t_53[k] = ab_x[k] * fd_33[k]
                      + gd_33[k];

            t_54[k] = ab_x[k] * fd_34[k]
                      + gd_34[k];
        }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_x, ab_y, ab_z, fd_33, fd_34, fd_35, \
                         gd_35, gd_51, gd_52, gd_53, gd_59 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_55[k] = ab_x[k] * fd_35[k]
                      + gd_35[k];

            t_56[k] = ab_y[k] * fd_33[k]
                      + gd_51[k];

            t_57[k] = ab_y[k] * fd_34[k]
                      + gd_52[k];

            t_58[k] = ab_y[k] * fd_35[k]
                      + gd_53[k];

            t_59[k] = ab_z[k] * fd_35[k]
                      + gd_59[k];
        }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, fd_36, fd_37, fd_38, fd_39, \
                         fd_40, gd_36, gd_37, gd_38, gd_39, gd_40 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_60[k] = ab_x[k] * fd_36[k]
                      + gd_36[k];

            t_61[k] = ab_x[k] * fd_37[k]
                      + gd_37[k];

            t_62[k] = ab_x[k] * fd_38[k]
                      + gd_38[k];

            t_63[k] = ab_x[k] * fd_39[k]
                      + gd_39[k];

            t_64[k] = ab_x[k] * fd_40[k]
                      + gd_40[k];
        }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_x, ab_y, ab_z, fd_39, fd_40, fd_41, \
                         gd_41, gd_63, gd_64, gd_65, gd_71 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_65[k] = ab_x[k] * fd_41[k]
                      + gd_41[k];

            t_66[k] = ab_y[k] * fd_39[k]
                      + gd_63[k];

            t_67[k] = ab_y[k] * fd_40[k]
                      + gd_64[k];

            t_68[k] = ab_y[k] * fd_41[k]
                      + gd_65[k];

            t_69[k] = ab_z[k] * fd_41[k]
                      + gd_71[k];
        }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ab_x, fd_42, fd_43, fd_44, fd_45, \
                         fd_46, gd_42, gd_43, gd_44, gd_45, gd_46 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_70[k] = ab_x[k] * fd_42[k]
                      + gd_42[k];

            t_71[k] = ab_x[k] * fd_43[k]
                      + gd_43[k];

            t_72[k] = ab_x[k] * fd_44[k]
                      + gd_44[k];

            t_73[k] = ab_x[k] * fd_45[k]
                      + gd_45[k];

            t_74[k] = ab_x[k] * fd_46[k]
                      + gd_46[k];
        }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ab_x, ab_y, ab_z, fd_45, fd_46, fd_47, \
                         gd_47, gd_69, gd_70, gd_71, gd_77 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_75[k] = ab_x[k] * fd_47[k]
                      + gd_47[k];

            t_76[k] = ab_y[k] * fd_45[k]
                      + gd_69[k];

            t_77[k] = ab_y[k] * fd_46[k]
                      + gd_70[k];

            t_78[k] = ab_y[k] * fd_47[k]
                      + gd_71[k];

            t_79[k] = ab_z[k] * fd_47[k]
                      + gd_77[k];
        }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ab_x, fd_48, fd_49, fd_50, fd_51, \
                         fd_52, gd_48, gd_49, gd_50, gd_51, gd_52 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_80[k] = ab_x[k] * fd_48[k]
                      + gd_48[k];

            t_81[k] = ab_x[k] * fd_49[k]
                      + gd_49[k];

            t_82[k] = ab_x[k] * fd_50[k]
                      + gd_50[k];

            t_83[k] = ab_x[k] * fd_51[k]
                      + gd_51[k];

            t_84[k] = ab_x[k] * fd_52[k]
                      + gd_52[k];
        }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_x, ab_y, ab_z, fd_51, fd_52, fd_53, \
                         gd_53, gd_75, gd_76, gd_77, gd_83 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_85[k] = ab_x[k] * fd_53[k]
                      + gd_53[k];

            t_86[k] = ab_y[k] * fd_51[k]
                      + gd_75[k];

            t_87[k] = ab_y[k] * fd_52[k]
                      + gd_76[k];

            t_88[k] = ab_y[k] * fd_53[k]
                      + gd_77[k];

            t_89[k] = ab_z[k] * fd_53[k]
                      + gd_83[k];
        }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_x, fd_54, fd_55, fd_56, fd_57, \
                         fd_58, gd_54, gd_55, gd_56, gd_57, gd_58 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_90[k] = ab_x[k] * fd_54[k]
                      + gd_54[k];

            t_91[k] = ab_x[k] * fd_55[k]
                      + gd_55[k];

            t_92[k] = ab_x[k] * fd_56[k]
                      + gd_56[k];

            t_93[k] = ab_x[k] * fd_57[k]
                      + gd_57[k];

            t_94[k] = ab_x[k] * fd_58[k]
                      + gd_58[k];
        }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ab_x, ab_y, ab_z, fd_57, fd_58, fd_59, \
                         gd_59, gd_81, gd_82, gd_83, gd_89 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_95[k] = ab_x[k] * fd_59[k]
                      + gd_59[k];

            t_96[k] = ab_y[k] * fd_57[k]
                      + gd_81[k];

            t_97[k] = ab_y[k] * fd_58[k]
                      + gd_82[k];

            t_98[k] = ab_y[k] * fd_59[k]
                      + gd_83[k];

            t_99[k] = ab_z[k] * fd_59[k]
                      + gd_89[k];
        }
    }
}

auto
compute_hrr_ff(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t df, const size_t dg, const size_t ncomps, const size_t nmax) -> void
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
        auto *t_90 = buffer.data(target + 90 * ncomps + c);
        auto *t_91 = buffer.data(target + 91 * ncomps + c);
        auto *t_92 = buffer.data(target + 92 * ncomps + c);
        auto *t_93 = buffer.data(target + 93 * ncomps + c);
        auto *t_94 = buffer.data(target + 94 * ncomps + c);
        auto *t_95 = buffer.data(target + 95 * ncomps + c);
        auto *t_96 = buffer.data(target + 96 * ncomps + c);
        auto *t_97 = buffer.data(target + 97 * ncomps + c);
        auto *t_98 = buffer.data(target + 98 * ncomps + c);
        auto *t_99 = buffer.data(target + 99 * ncomps + c);

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *df_0 = buffer.data(df + 0 * ncomps + c);
        const auto *df_1 = buffer.data(df + 1 * ncomps + c);
        const auto *df_2 = buffer.data(df + 2 * ncomps + c);
        const auto *df_3 = buffer.data(df + 3 * ncomps + c);
        const auto *df_4 = buffer.data(df + 4 * ncomps + c);
        const auto *df_5 = buffer.data(df + 5 * ncomps + c);
        const auto *df_6 = buffer.data(df + 6 * ncomps + c);
        const auto *df_7 = buffer.data(df + 7 * ncomps + c);
        const auto *df_8 = buffer.data(df + 8 * ncomps + c);
        const auto *df_9 = buffer.data(df + 9 * ncomps + c);
        const auto *df_10 = buffer.data(df + 10 * ncomps + c);
        const auto *df_11 = buffer.data(df + 11 * ncomps + c);
        const auto *df_12 = buffer.data(df + 12 * ncomps + c);
        const auto *df_13 = buffer.data(df + 13 * ncomps + c);
        const auto *df_14 = buffer.data(df + 14 * ncomps + c);
        const auto *df_15 = buffer.data(df + 15 * ncomps + c);
        const auto *df_16 = buffer.data(df + 16 * ncomps + c);
        const auto *df_17 = buffer.data(df + 17 * ncomps + c);
        const auto *df_18 = buffer.data(df + 18 * ncomps + c);
        const auto *df_19 = buffer.data(df + 19 * ncomps + c);
        const auto *df_20 = buffer.data(df + 20 * ncomps + c);
        const auto *df_21 = buffer.data(df + 21 * ncomps + c);
        const auto *df_22 = buffer.data(df + 22 * ncomps + c);
        const auto *df_23 = buffer.data(df + 23 * ncomps + c);
        const auto *df_24 = buffer.data(df + 24 * ncomps + c);
        const auto *df_25 = buffer.data(df + 25 * ncomps + c);
        const auto *df_26 = buffer.data(df + 26 * ncomps + c);
        const auto *df_27 = buffer.data(df + 27 * ncomps + c);
        const auto *df_28 = buffer.data(df + 28 * ncomps + c);
        const auto *df_29 = buffer.data(df + 29 * ncomps + c);
        const auto *df_30 = buffer.data(df + 30 * ncomps + c);
        const auto *df_31 = buffer.data(df + 31 * ncomps + c);
        const auto *df_32 = buffer.data(df + 32 * ncomps + c);
        const auto *df_33 = buffer.data(df + 33 * ncomps + c);
        const auto *df_34 = buffer.data(df + 34 * ncomps + c);
        const auto *df_35 = buffer.data(df + 35 * ncomps + c);
        const auto *df_36 = buffer.data(df + 36 * ncomps + c);
        const auto *df_37 = buffer.data(df + 37 * ncomps + c);
        const auto *df_38 = buffer.data(df + 38 * ncomps + c);
        const auto *df_39 = buffer.data(df + 39 * ncomps + c);
        const auto *df_40 = buffer.data(df + 40 * ncomps + c);
        const auto *df_41 = buffer.data(df + 41 * ncomps + c);
        const auto *df_42 = buffer.data(df + 42 * ncomps + c);
        const auto *df_43 = buffer.data(df + 43 * ncomps + c);
        const auto *df_44 = buffer.data(df + 44 * ncomps + c);
        const auto *df_45 = buffer.data(df + 45 * ncomps + c);
        const auto *df_46 = buffer.data(df + 46 * ncomps + c);
        const auto *df_47 = buffer.data(df + 47 * ncomps + c);
        const auto *df_48 = buffer.data(df + 48 * ncomps + c);
        const auto *df_49 = buffer.data(df + 49 * ncomps + c);
        const auto *df_50 = buffer.data(df + 50 * ncomps + c);
        const auto *df_51 = buffer.data(df + 51 * ncomps + c);
        const auto *df_52 = buffer.data(df + 52 * ncomps + c);
        const auto *df_53 = buffer.data(df + 53 * ncomps + c);
        const auto *df_54 = buffer.data(df + 54 * ncomps + c);
        const auto *df_55 = buffer.data(df + 55 * ncomps + c);
        const auto *df_56 = buffer.data(df + 56 * ncomps + c);
        const auto *df_57 = buffer.data(df + 57 * ncomps + c);
        const auto *df_58 = buffer.data(df + 58 * ncomps + c);
        const auto *df_59 = buffer.data(df + 59 * ncomps + c);

        const auto *dg_0 = buffer.data(dg + 0 * ncomps + c);
        const auto *dg_1 = buffer.data(dg + 1 * ncomps + c);
        const auto *dg_2 = buffer.data(dg + 2 * ncomps + c);
        const auto *dg_3 = buffer.data(dg + 3 * ncomps + c);
        const auto *dg_4 = buffer.data(dg + 4 * ncomps + c);
        const auto *dg_5 = buffer.data(dg + 5 * ncomps + c);
        const auto *dg_6 = buffer.data(dg + 6 * ncomps + c);
        const auto *dg_7 = buffer.data(dg + 7 * ncomps + c);
        const auto *dg_8 = buffer.data(dg + 8 * ncomps + c);
        const auto *dg_9 = buffer.data(dg + 9 * ncomps + c);
        const auto *dg_15 = buffer.data(dg + 15 * ncomps + c);
        const auto *dg_16 = buffer.data(dg + 16 * ncomps + c);
        const auto *dg_17 = buffer.data(dg + 17 * ncomps + c);
        const auto *dg_18 = buffer.data(dg + 18 * ncomps + c);
        const auto *dg_19 = buffer.data(dg + 19 * ncomps + c);
        const auto *dg_20 = buffer.data(dg + 20 * ncomps + c);
        const auto *dg_21 = buffer.data(dg + 21 * ncomps + c);
        const auto *dg_22 = buffer.data(dg + 22 * ncomps + c);
        const auto *dg_23 = buffer.data(dg + 23 * ncomps + c);
        const auto *dg_24 = buffer.data(dg + 24 * ncomps + c);
        const auto *dg_30 = buffer.data(dg + 30 * ncomps + c);
        const auto *dg_31 = buffer.data(dg + 31 * ncomps + c);
        const auto *dg_32 = buffer.data(dg + 32 * ncomps + c);
        const auto *dg_33 = buffer.data(dg + 33 * ncomps + c);
        const auto *dg_34 = buffer.data(dg + 34 * ncomps + c);
        const auto *dg_35 = buffer.data(dg + 35 * ncomps + c);
        const auto *dg_36 = buffer.data(dg + 36 * ncomps + c);
        const auto *dg_37 = buffer.data(dg + 37 * ncomps + c);
        const auto *dg_38 = buffer.data(dg + 38 * ncomps + c);
        const auto *dg_39 = buffer.data(dg + 39 * ncomps + c);
        const auto *dg_45 = buffer.data(dg + 45 * ncomps + c);
        const auto *dg_46 = buffer.data(dg + 46 * ncomps + c);
        const auto *dg_47 = buffer.data(dg + 47 * ncomps + c);
        const auto *dg_48 = buffer.data(dg + 48 * ncomps + c);
        const auto *dg_49 = buffer.data(dg + 49 * ncomps + c);
        const auto *dg_50 = buffer.data(dg + 50 * ncomps + c);
        const auto *dg_51 = buffer.data(dg + 51 * ncomps + c);
        const auto *dg_52 = buffer.data(dg + 52 * ncomps + c);
        const auto *dg_53 = buffer.data(dg + 53 * ncomps + c);
        const auto *dg_54 = buffer.data(dg + 54 * ncomps + c);
        const auto *dg_55 = buffer.data(dg + 55 * ncomps + c);
        const auto *dg_56 = buffer.data(dg + 56 * ncomps + c);
        const auto *dg_57 = buffer.data(dg + 57 * ncomps + c);
        const auto *dg_58 = buffer.data(dg + 58 * ncomps + c);
        const auto *dg_60 = buffer.data(dg + 60 * ncomps + c);
        const auto *dg_61 = buffer.data(dg + 61 * ncomps + c);
        const auto *dg_62 = buffer.data(dg + 62 * ncomps + c);
        const auto *dg_63 = buffer.data(dg + 63 * ncomps + c);
        const auto *dg_64 = buffer.data(dg + 64 * ncomps + c);
        const auto *dg_65 = buffer.data(dg + 65 * ncomps + c);
        const auto *dg_66 = buffer.data(dg + 66 * ncomps + c);
        const auto *dg_67 = buffer.data(dg + 67 * ncomps + c);
        const auto *dg_68 = buffer.data(dg + 68 * ncomps + c);
        const auto *dg_69 = buffer.data(dg + 69 * ncomps + c);
        const auto *dg_70 = buffer.data(dg + 70 * ncomps + c);
        const auto *dg_71 = buffer.data(dg + 71 * ncomps + c);
        const auto *dg_72 = buffer.data(dg + 72 * ncomps + c);
        const auto *dg_73 = buffer.data(dg + 73 * ncomps + c);
        const auto *dg_75 = buffer.data(dg + 75 * ncomps + c);
        const auto *dg_76 = buffer.data(dg + 76 * ncomps + c);
        const auto *dg_77 = buffer.data(dg + 77 * ncomps + c);
        const auto *dg_78 = buffer.data(dg + 78 * ncomps + c);
        const auto *dg_79 = buffer.data(dg + 79 * ncomps + c);
        const auto *dg_80 = buffer.data(dg + 80 * ncomps + c);
        const auto *dg_81 = buffer.data(dg + 81 * ncomps + c);
        const auto *dg_82 = buffer.data(dg + 82 * ncomps + c);
        const auto *dg_83 = buffer.data(dg + 83 * ncomps + c);
        const auto *dg_84 = buffer.data(dg + 84 * ncomps + c);
        const auto *dg_85 = buffer.data(dg + 85 * ncomps + c);
        const auto *dg_86 = buffer.data(dg + 86 * ncomps + c);
        const auto *dg_87 = buffer.data(dg + 87 * ncomps + c);
        const auto *dg_88 = buffer.data(dg + 88 * ncomps + c);
        const auto *dg_89 = buffer.data(dg + 89 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, df_0, df_1, df_2, df_3, df_4, dg_0, \
                         dg_1, dg_2, dg_3, dg_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = -ab_x[k] * df_0[k]
                     + dg_0[k];

            t_1[k] = -ab_x[k] * df_1[k]
                     + dg_1[k];

            t_2[k] = -ab_x[k] * df_2[k]
                     + dg_2[k];

            t_3[k] = -ab_x[k] * df_3[k]
                     + dg_3[k];

            t_4[k] = -ab_x[k] * df_4[k]
                     + dg_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, df_5, df_6, df_7, df_8, df_9, dg_5, \
                         dg_6, dg_7, dg_8, dg_9 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = -ab_x[k] * df_5[k]
                     + dg_5[k];

            t_6[k] = -ab_x[k] * df_6[k]
                     + dg_6[k];

            t_7[k] = -ab_x[k] * df_7[k]
                     + dg_7[k];

            t_8[k] = -ab_x[k] * df_8[k]
                     + dg_8[k];

            t_9[k] = -ab_x[k] * df_9[k]
                     + dg_9[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, df_10, df_11, df_12, df_13, \
                         df_14, dg_15, dg_16, dg_17, dg_18, dg_19 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = -ab_x[k] * df_10[k]
                      + dg_15[k];

            t_11[k] = -ab_x[k] * df_11[k]
                      + dg_16[k];

            t_12[k] = -ab_x[k] * df_12[k]
                      + dg_17[k];

            t_13[k] = -ab_x[k] * df_13[k]
                      + dg_18[k];

            t_14[k] = -ab_x[k] * df_14[k]
                      + dg_19[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, df_15, df_16, df_17, df_18, \
                         df_19, dg_20, dg_21, dg_22, dg_23, dg_24 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = -ab_x[k] * df_15[k]
                      + dg_20[k];

            t_16[k] = -ab_x[k] * df_16[k]
                      + dg_21[k];

            t_17[k] = -ab_x[k] * df_17[k]
                      + dg_22[k];

            t_18[k] = -ab_x[k] * df_18[k]
                      + dg_23[k];

            t_19[k] = -ab_x[k] * df_19[k]
                      + dg_24[k];
        }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, df_20, df_21, df_22, df_23, \
                         df_24, dg_30, dg_31, dg_32, dg_33, dg_34 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_20[k] = -ab_x[k] * df_20[k]
                      + dg_30[k];

            t_21[k] = -ab_x[k] * df_21[k]
                      + dg_31[k];

            t_22[k] = -ab_x[k] * df_22[k]
                      + dg_32[k];

            t_23[k] = -ab_x[k] * df_23[k]
                      + dg_33[k];

            t_24[k] = -ab_x[k] * df_24[k]
                      + dg_34[k];
        }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, df_25, df_26, df_27, df_28, \
                         df_29, dg_35, dg_36, dg_37, dg_38, dg_39 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_25[k] = -ab_x[k] * df_25[k]
                      + dg_35[k];

            t_26[k] = -ab_x[k] * df_26[k]
                      + dg_36[k];

            t_27[k] = -ab_x[k] * df_27[k]
                      + dg_37[k];

            t_28[k] = -ab_x[k] * df_28[k]
                      + dg_38[k];

            t_29[k] = -ab_x[k] * df_29[k]
                      + dg_39[k];
        }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, df_30, df_31, df_32, df_33, \
                         df_34, dg_45, dg_46, dg_47, dg_48, dg_49 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_30[k] = -ab_x[k] * df_30[k]
                      + dg_45[k];

            t_31[k] = -ab_x[k] * df_31[k]
                      + dg_46[k];

            t_32[k] = -ab_x[k] * df_32[k]
                      + dg_47[k];

            t_33[k] = -ab_x[k] * df_33[k]
                      + dg_48[k];

            t_34[k] = -ab_x[k] * df_34[k]
                      + dg_49[k];
        }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, df_35, df_36, df_37, df_38, \
                         df_39, dg_50, dg_51, dg_52, dg_53, dg_54 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_35[k] = -ab_x[k] * df_35[k]
                      + dg_50[k];

            t_36[k] = -ab_x[k] * df_36[k]
                      + dg_51[k];

            t_37[k] = -ab_x[k] * df_37[k]
                      + dg_52[k];

            t_38[k] = -ab_x[k] * df_38[k]
                      + dg_53[k];

            t_39[k] = -ab_x[k] * df_39[k]
                      + dg_54[k];
        }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_x, df_40, df_41, df_42, df_43, \
                         df_44, dg_60, dg_61, dg_62, dg_63, dg_64 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_40[k] = -ab_x[k] * df_40[k]
                      + dg_60[k];

            t_41[k] = -ab_x[k] * df_41[k]
                      + dg_61[k];

            t_42[k] = -ab_x[k] * df_42[k]
                      + dg_62[k];

            t_43[k] = -ab_x[k] * df_43[k]
                      + dg_63[k];

            t_44[k] = -ab_x[k] * df_44[k]
                      + dg_64[k];
        }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, df_45, df_46, df_47, df_48, \
                         df_49, dg_65, dg_66, dg_67, dg_68, dg_69 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_45[k] = -ab_x[k] * df_45[k]
                      + dg_65[k];

            t_46[k] = -ab_x[k] * df_46[k]
                      + dg_66[k];

            t_47[k] = -ab_x[k] * df_47[k]
                      + dg_67[k];

            t_48[k] = -ab_x[k] * df_48[k]
                      + dg_68[k];

            t_49[k] = -ab_x[k] * df_49[k]
                      + dg_69[k];
        }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, df_50, df_51, df_52, df_53, \
                         df_54, dg_75, dg_76, dg_77, dg_78, dg_79 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_50[k] = -ab_x[k] * df_50[k]
                      + dg_75[k];

            t_51[k] = -ab_x[k] * df_51[k]
                      + dg_76[k];

            t_52[k] = -ab_x[k] * df_52[k]
                      + dg_77[k];

            t_53[k] = -ab_x[k] * df_53[k]
                      + dg_78[k];

            t_54[k] = -ab_x[k] * df_54[k]
                      + dg_79[k];
        }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_x, df_55, df_56, df_57, df_58, \
                         df_59, dg_80, dg_81, dg_82, dg_83, dg_84 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_55[k] = -ab_x[k] * df_55[k]
                      + dg_80[k];

            t_56[k] = -ab_x[k] * df_56[k]
                      + dg_81[k];

            t_57[k] = -ab_x[k] * df_57[k]
                      + dg_82[k];

            t_58[k] = -ab_x[k] * df_58[k]
                      + dg_83[k];

            t_59[k] = -ab_x[k] * df_59[k]
                      + dg_84[k];
        }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_y, df_30, df_31, df_32, df_33, \
                         df_34, dg_46, dg_48, dg_49, dg_51, dg_52 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_60[k] = -ab_y[k] * df_30[k]
                      + dg_46[k];

            t_61[k] = -ab_y[k] * df_31[k]
                      + dg_48[k];

            t_62[k] = -ab_y[k] * df_32[k]
                      + dg_49[k];

            t_63[k] = -ab_y[k] * df_33[k]
                      + dg_51[k];

            t_64[k] = -ab_y[k] * df_34[k]
                      + dg_52[k];
        }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_y, df_35, df_36, df_37, df_38, \
                         df_39, dg_53, dg_55, dg_56, dg_57, dg_58 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_65[k] = -ab_y[k] * df_35[k]
                      + dg_53[k];

            t_66[k] = -ab_y[k] * df_36[k]
                      + dg_55[k];

            t_67[k] = -ab_y[k] * df_37[k]
                      + dg_56[k];

            t_68[k] = -ab_y[k] * df_38[k]
                      + dg_57[k];

            t_69[k] = -ab_y[k] * df_39[k]
                      + dg_58[k];
        }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ab_y, df_40, df_41, df_42, df_43, \
                         df_44, dg_61, dg_63, dg_64, dg_66, dg_67 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_70[k] = -ab_y[k] * df_40[k]
                      + dg_61[k];

            t_71[k] = -ab_y[k] * df_41[k]
                      + dg_63[k];

            t_72[k] = -ab_y[k] * df_42[k]
                      + dg_64[k];

            t_73[k] = -ab_y[k] * df_43[k]
                      + dg_66[k];

            t_74[k] = -ab_y[k] * df_44[k]
                      + dg_67[k];
        }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ab_y, df_45, df_46, df_47, df_48, \
                         df_49, dg_68, dg_70, dg_71, dg_72, dg_73 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_75[k] = -ab_y[k] * df_45[k]
                      + dg_68[k];

            t_76[k] = -ab_y[k] * df_46[k]
                      + dg_70[k];

            t_77[k] = -ab_y[k] * df_47[k]
                      + dg_71[k];

            t_78[k] = -ab_y[k] * df_48[k]
                      + dg_72[k];

            t_79[k] = -ab_y[k] * df_49[k]
                      + dg_73[k];
        }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ab_y, df_50, df_51, df_52, df_53, \
                         df_54, dg_76, dg_78, dg_79, dg_81, dg_82 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_80[k] = -ab_y[k] * df_50[k]
                      + dg_76[k];

            t_81[k] = -ab_y[k] * df_51[k]
                      + dg_78[k];

            t_82[k] = -ab_y[k] * df_52[k]
                      + dg_79[k];

            t_83[k] = -ab_y[k] * df_53[k]
                      + dg_81[k];

            t_84[k] = -ab_y[k] * df_54[k]
                      + dg_82[k];
        }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_y, df_55, df_56, df_57, df_58, \
                         df_59, dg_83, dg_85, dg_86, dg_87, dg_88 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_85[k] = -ab_y[k] * df_55[k]
                      + dg_83[k];

            t_86[k] = -ab_y[k] * df_56[k]
                      + dg_85[k];

            t_87[k] = -ab_y[k] * df_57[k]
                      + dg_86[k];

            t_88[k] = -ab_y[k] * df_58[k]
                      + dg_87[k];

            t_89[k] = -ab_y[k] * df_59[k]
                      + dg_88[k];
        }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_z, df_50, df_51, df_52, df_53, \
                         df_54, dg_77, dg_79, dg_80, dg_82, dg_83 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_90[k] = -ab_z[k] * df_50[k]
                      + dg_77[k];

            t_91[k] = -ab_z[k] * df_51[k]
                      + dg_79[k];

            t_92[k] = -ab_z[k] * df_52[k]
                      + dg_80[k];

            t_93[k] = -ab_z[k] * df_53[k]
                      + dg_82[k];

            t_94[k] = -ab_z[k] * df_54[k]
                      + dg_83[k];
        }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ab_z, df_55, df_56, df_57, df_58, \
                         df_59, dg_84, dg_86, dg_87, dg_88, dg_89 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_95[k] = -ab_z[k] * df_55[k]
                      + dg_84[k];

            t_96[k] = -ab_z[k] * df_56[k]
                      + dg_86[k];

            t_97[k] = -ab_z[k] * df_57[k]
                      + dg_87[k];

            t_98[k] = -ab_z[k] * df_58[k]
                      + dg_88[k];

            t_99[k] = -ab_z[k] * df_59[k]
                      + dg_89[k];
        }
    }
}

}  // namespace simdtrf
