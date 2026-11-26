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


#include "SimdTransferGeom010ZPI.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_geom_010z_pi(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                         const size_t target, const size_t si_1, const size_t si_0,
                         const size_t sk_1, const size_t ncomps, const size_t nmax) -> void
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

        const auto *si_1_0 = buffer.data(si_1 + 0 * ncomps + c);
        const auto *si_1_1 = buffer.data(si_1 + 1 * ncomps + c);
        const auto *si_1_2 = buffer.data(si_1 + 2 * ncomps + c);
        const auto *si_1_3 = buffer.data(si_1 + 3 * ncomps + c);
        const auto *si_1_4 = buffer.data(si_1 + 4 * ncomps + c);
        const auto *si_1_5 = buffer.data(si_1 + 5 * ncomps + c);
        const auto *si_1_6 = buffer.data(si_1 + 6 * ncomps + c);
        const auto *si_1_7 = buffer.data(si_1 + 7 * ncomps + c);
        const auto *si_1_8 = buffer.data(si_1 + 8 * ncomps + c);
        const auto *si_1_9 = buffer.data(si_1 + 9 * ncomps + c);
        const auto *si_1_10 = buffer.data(si_1 + 10 * ncomps + c);
        const auto *si_1_11 = buffer.data(si_1 + 11 * ncomps + c);
        const auto *si_1_12 = buffer.data(si_1 + 12 * ncomps + c);
        const auto *si_1_13 = buffer.data(si_1 + 13 * ncomps + c);
        const auto *si_1_14 = buffer.data(si_1 + 14 * ncomps + c);
        const auto *si_1_15 = buffer.data(si_1 + 15 * ncomps + c);
        const auto *si_1_16 = buffer.data(si_1 + 16 * ncomps + c);
        const auto *si_1_17 = buffer.data(si_1 + 17 * ncomps + c);
        const auto *si_1_18 = buffer.data(si_1 + 18 * ncomps + c);
        const auto *si_1_19 = buffer.data(si_1 + 19 * ncomps + c);
        const auto *si_1_20 = buffer.data(si_1 + 20 * ncomps + c);
        const auto *si_1_21 = buffer.data(si_1 + 21 * ncomps + c);
        const auto *si_1_22 = buffer.data(si_1 + 22 * ncomps + c);
        const auto *si_1_23 = buffer.data(si_1 + 23 * ncomps + c);
        const auto *si_1_24 = buffer.data(si_1 + 24 * ncomps + c);
        const auto *si_1_25 = buffer.data(si_1 + 25 * ncomps + c);
        const auto *si_1_26 = buffer.data(si_1 + 26 * ncomps + c);
        const auto *si_1_27 = buffer.data(si_1 + 27 * ncomps + c);

        const auto *si_0_0 = buffer.data(si_0 + 0 * ncomps + c);
        const auto *si_0_1 = buffer.data(si_0 + 1 * ncomps + c);
        const auto *si_0_2 = buffer.data(si_0 + 2 * ncomps + c);
        const auto *si_0_3 = buffer.data(si_0 + 3 * ncomps + c);
        const auto *si_0_4 = buffer.data(si_0 + 4 * ncomps + c);
        const auto *si_0_5 = buffer.data(si_0 + 5 * ncomps + c);
        const auto *si_0_6 = buffer.data(si_0 + 6 * ncomps + c);
        const auto *si_0_7 = buffer.data(si_0 + 7 * ncomps + c);
        const auto *si_0_8 = buffer.data(si_0 + 8 * ncomps + c);
        const auto *si_0_9 = buffer.data(si_0 + 9 * ncomps + c);
        const auto *si_0_10 = buffer.data(si_0 + 10 * ncomps + c);
        const auto *si_0_11 = buffer.data(si_0 + 11 * ncomps + c);
        const auto *si_0_12 = buffer.data(si_0 + 12 * ncomps + c);
        const auto *si_0_13 = buffer.data(si_0 + 13 * ncomps + c);
        const auto *si_0_14 = buffer.data(si_0 + 14 * ncomps + c);
        const auto *si_0_15 = buffer.data(si_0 + 15 * ncomps + c);
        const auto *si_0_16 = buffer.data(si_0 + 16 * ncomps + c);
        const auto *si_0_17 = buffer.data(si_0 + 17 * ncomps + c);
        const auto *si_0_18 = buffer.data(si_0 + 18 * ncomps + c);
        const auto *si_0_19 = buffer.data(si_0 + 19 * ncomps + c);
        const auto *si_0_20 = buffer.data(si_0 + 20 * ncomps + c);
        const auto *si_0_21 = buffer.data(si_0 + 21 * ncomps + c);
        const auto *si_0_22 = buffer.data(si_0 + 22 * ncomps + c);
        const auto *si_0_23 = buffer.data(si_0 + 23 * ncomps + c);
        const auto *si_0_24 = buffer.data(si_0 + 24 * ncomps + c);
        const auto *si_0_25 = buffer.data(si_0 + 25 * ncomps + c);
        const auto *si_0_26 = buffer.data(si_0 + 26 * ncomps + c);
        const auto *si_0_27 = buffer.data(si_0 + 27 * ncomps + c);

        const auto *sk_1_0 = buffer.data(sk_1 + 0 * ncomps + c);
        const auto *sk_1_1 = buffer.data(sk_1 + 1 * ncomps + c);
        const auto *sk_1_2 = buffer.data(sk_1 + 2 * ncomps + c);
        const auto *sk_1_3 = buffer.data(sk_1 + 3 * ncomps + c);
        const auto *sk_1_4 = buffer.data(sk_1 + 4 * ncomps + c);
        const auto *sk_1_5 = buffer.data(sk_1 + 5 * ncomps + c);
        const auto *sk_1_6 = buffer.data(sk_1 + 6 * ncomps + c);
        const auto *sk_1_7 = buffer.data(sk_1 + 7 * ncomps + c);
        const auto *sk_1_8 = buffer.data(sk_1 + 8 * ncomps + c);
        const auto *sk_1_9 = buffer.data(sk_1 + 9 * ncomps + c);
        const auto *sk_1_10 = buffer.data(sk_1 + 10 * ncomps + c);
        const auto *sk_1_11 = buffer.data(sk_1 + 11 * ncomps + c);
        const auto *sk_1_12 = buffer.data(sk_1 + 12 * ncomps + c);
        const auto *sk_1_13 = buffer.data(sk_1 + 13 * ncomps + c);
        const auto *sk_1_14 = buffer.data(sk_1 + 14 * ncomps + c);
        const auto *sk_1_15 = buffer.data(sk_1 + 15 * ncomps + c);
        const auto *sk_1_16 = buffer.data(sk_1 + 16 * ncomps + c);
        const auto *sk_1_17 = buffer.data(sk_1 + 17 * ncomps + c);
        const auto *sk_1_18 = buffer.data(sk_1 + 18 * ncomps + c);
        const auto *sk_1_19 = buffer.data(sk_1 + 19 * ncomps + c);
        const auto *sk_1_20 = buffer.data(sk_1 + 20 * ncomps + c);
        const auto *sk_1_21 = buffer.data(sk_1 + 21 * ncomps + c);
        const auto *sk_1_22 = buffer.data(sk_1 + 22 * ncomps + c);
        const auto *sk_1_23 = buffer.data(sk_1 + 23 * ncomps + c);
        const auto *sk_1_24 = buffer.data(sk_1 + 24 * ncomps + c);
        const auto *sk_1_25 = buffer.data(sk_1 + 25 * ncomps + c);
        const auto *sk_1_26 = buffer.data(sk_1 + 26 * ncomps + c);
        const auto *sk_1_27 = buffer.data(sk_1 + 27 * ncomps + c);
        const auto *sk_1_28 = buffer.data(sk_1 + 28 * ncomps + c);
        const auto *sk_1_29 = buffer.data(sk_1 + 29 * ncomps + c);
        const auto *sk_1_30 = buffer.data(sk_1 + 30 * ncomps + c);
        const auto *sk_1_31 = buffer.data(sk_1 + 31 * ncomps + c);
        const auto *sk_1_32 = buffer.data(sk_1 + 32 * ncomps + c);
        const auto *sk_1_33 = buffer.data(sk_1 + 33 * ncomps + c);
        const auto *sk_1_34 = buffer.data(sk_1 + 34 * ncomps + c);
        const auto *sk_1_35 = buffer.data(sk_1 + 35 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, si_1_0, si_1_1, si_1_2, si_1_3, \
                         si_1_4, sk_1_0, sk_1_1, sk_1_2, sk_1_3, \
                         sk_1_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = -ab_x[k] * si_1_0[k]
                     + sk_1_0[k];

            t_1[k] = -ab_x[k] * si_1_1[k]
                     + sk_1_1[k];

            t_2[k] = -ab_x[k] * si_1_2[k]
                     + sk_1_2[k];

            t_3[k] = -ab_x[k] * si_1_3[k]
                     + sk_1_3[k];

            t_4[k] = -ab_x[k] * si_1_4[k]
                     + sk_1_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, si_1_5, si_1_6, si_1_7, si_1_8, \
                         si_1_9, sk_1_5, sk_1_6, sk_1_7, sk_1_8, \
                         sk_1_9 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = -ab_x[k] * si_1_5[k]
                     + sk_1_5[k];

            t_6[k] = -ab_x[k] * si_1_6[k]
                     + sk_1_6[k];

            t_7[k] = -ab_x[k] * si_1_7[k]
                     + sk_1_7[k];

            t_8[k] = -ab_x[k] * si_1_8[k]
                     + sk_1_8[k];

            t_9[k] = -ab_x[k] * si_1_9[k]
                     + sk_1_9[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, si_1_10, si_1_11, si_1_12, \
                         si_1_13, si_1_14, sk_1_10, sk_1_11, sk_1_12, sk_1_13, \
                         sk_1_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = -ab_x[k] * si_1_10[k]
                      + sk_1_10[k];

            t_11[k] = -ab_x[k] * si_1_11[k]
                      + sk_1_11[k];

            t_12[k] = -ab_x[k] * si_1_12[k]
                      + sk_1_12[k];

            t_13[k] = -ab_x[k] * si_1_13[k]
                      + sk_1_13[k];

            t_14[k] = -ab_x[k] * si_1_14[k]
                      + sk_1_14[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, si_1_15, si_1_16, si_1_17, \
                         si_1_18, si_1_19, sk_1_15, sk_1_16, sk_1_17, sk_1_18, \
                         sk_1_19 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = -ab_x[k] * si_1_15[k]
                      + sk_1_15[k];

            t_16[k] = -ab_x[k] * si_1_16[k]
                      + sk_1_16[k];

            t_17[k] = -ab_x[k] * si_1_17[k]
                      + sk_1_17[k];

            t_18[k] = -ab_x[k] * si_1_18[k]
                      + sk_1_18[k];

            t_19[k] = -ab_x[k] * si_1_19[k]
                      + sk_1_19[k];
        }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, si_1_20, si_1_21, si_1_22, \
                         si_1_23, si_1_24, sk_1_20, sk_1_21, sk_1_22, sk_1_23, \
                         sk_1_24 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_20[k] = -ab_x[k] * si_1_20[k]
                      + sk_1_20[k];

            t_21[k] = -ab_x[k] * si_1_21[k]
                      + sk_1_21[k];

            t_22[k] = -ab_x[k] * si_1_22[k]
                      + sk_1_22[k];

            t_23[k] = -ab_x[k] * si_1_23[k]
                      + sk_1_23[k];

            t_24[k] = -ab_x[k] * si_1_24[k]
                      + sk_1_24[k];
        }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, ab_x, ab_y, si_1_0, si_1_25, si_1_26, \
                         si_1_27, sk_1_1, sk_1_25, sk_1_26, sk_1_27 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_25[k] = -ab_x[k] * si_1_25[k]
                      + sk_1_25[k];

            t_26[k] = -ab_x[k] * si_1_26[k]
                      + sk_1_26[k];

            t_27[k] = -ab_x[k] * si_1_27[k]
                      + sk_1_27[k];

            t_28[k] = -ab_y[k] * si_1_0[k]
                      + sk_1_1[k];
        }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, ab_y, si_1_1, si_1_2, si_1_3, si_1_4, \
                         si_1_5, sk_1_3, sk_1_4, sk_1_6, sk_1_7, \
                         sk_1_8 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_29[k] = -ab_y[k] * si_1_1[k]
                      + sk_1_3[k];

            t_30[k] = -ab_y[k] * si_1_2[k]
                      + sk_1_4[k];

            t_31[k] = -ab_y[k] * si_1_3[k]
                      + sk_1_6[k];

            t_32[k] = -ab_y[k] * si_1_4[k]
                      + sk_1_7[k];

            t_33[k] = -ab_y[k] * si_1_5[k]
                      + sk_1_8[k];
        }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, t_38, ab_y, si_1_6, si_1_7, si_1_8, si_1_9, \
                         si_1_10, sk_1_10, sk_1_11, sk_1_12, sk_1_13, \
                         sk_1_15 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_34[k] = -ab_y[k] * si_1_6[k]
                      + sk_1_10[k];

            t_35[k] = -ab_y[k] * si_1_7[k]
                      + sk_1_11[k];

            t_36[k] = -ab_y[k] * si_1_8[k]
                      + sk_1_12[k];

            t_37[k] = -ab_y[k] * si_1_9[k]
                      + sk_1_13[k];

            t_38[k] = -ab_y[k] * si_1_10[k]
                      + sk_1_15[k];
        }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, t_43, ab_y, si_1_11, si_1_12, si_1_13, \
                         si_1_14, si_1_15, sk_1_16, sk_1_17, sk_1_18, sk_1_19, \
                         sk_1_21 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_39[k] = -ab_y[k] * si_1_11[k]
                      + sk_1_16[k];

            t_40[k] = -ab_y[k] * si_1_12[k]
                      + sk_1_17[k];

            t_41[k] = -ab_y[k] * si_1_13[k]
                      + sk_1_18[k];

            t_42[k] = -ab_y[k] * si_1_14[k]
                      + sk_1_19[k];

            t_43[k] = -ab_y[k] * si_1_15[k]
                      + sk_1_21[k];
        }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, ab_y, si_1_16, si_1_17, si_1_18, \
                         si_1_19, si_1_20, sk_1_22, sk_1_23, sk_1_24, sk_1_25, \
                         sk_1_26 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_44[k] = -ab_y[k] * si_1_16[k]
                      + sk_1_22[k];

            t_45[k] = -ab_y[k] * si_1_17[k]
                      + sk_1_23[k];

            t_46[k] = -ab_y[k] * si_1_18[k]
                      + sk_1_24[k];

            t_47[k] = -ab_y[k] * si_1_19[k]
                      + sk_1_25[k];

            t_48[k] = -ab_y[k] * si_1_20[k]
                      + sk_1_26[k];
        }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, ab_y, si_1_21, si_1_22, si_1_23, \
                         si_1_24, si_1_25, sk_1_28, sk_1_29, sk_1_30, sk_1_31, \
                         sk_1_32 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_49[k] = -ab_y[k] * si_1_21[k]
                      + sk_1_28[k];

            t_50[k] = -ab_y[k] * si_1_22[k]
                      + sk_1_29[k];

            t_51[k] = -ab_y[k] * si_1_23[k]
                      + sk_1_30[k];

            t_52[k] = -ab_y[k] * si_1_24[k]
                      + sk_1_31[k];

            t_53[k] = -ab_y[k] * si_1_25[k]
                      + sk_1_32[k];
        }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, ab_y, ab_z, si_1_0, si_1_1, si_1_26, si_1_27, \
                         si_0_0, si_0_1, sk_1_2, sk_1_4, sk_1_33, \
                         sk_1_34 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_54[k] = -ab_y[k] * si_1_26[k]
                      + sk_1_33[k];

            t_55[k] = -ab_y[k] * si_1_27[k]
                      + sk_1_34[k];

            t_56[k] = -ab_z[k] * si_1_0[k]
                      + si_0_0[k]
                      + sk_1_2[k];

            t_57[k] = -ab_z[k] * si_1_1[k]
                      + si_0_1[k]
                      + sk_1_4[k];
        }

#pragma omp simd aligned(t_58, t_59, t_60, ab_z, si_1_2, si_1_3, si_1_4, si_0_2, si_0_3, \
                         si_0_4, sk_1_5, sk_1_7, sk_1_8 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_58[k] = -ab_z[k] * si_1_2[k]
                      + si_0_2[k]
                      + sk_1_5[k];

            t_59[k] = -ab_z[k] * si_1_3[k]
                      + si_0_3[k]
                      + sk_1_7[k];

            t_60[k] = -ab_z[k] * si_1_4[k]
                      + si_0_4[k]
                      + sk_1_8[k];
        }

#pragma omp simd aligned(t_61, t_62, t_63, ab_z, si_1_5, si_1_6, si_1_7, si_0_5, si_0_6, \
                         si_0_7, sk_1_9, sk_1_11, sk_1_12 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_61[k] = -ab_z[k] * si_1_5[k]
                      + si_0_5[k]
                      + sk_1_9[k];

            t_62[k] = -ab_z[k] * si_1_6[k]
                      + si_0_6[k]
                      + sk_1_11[k];

            t_63[k] = -ab_z[k] * si_1_7[k]
                      + si_0_7[k]
                      + sk_1_12[k];
        }

#pragma omp simd aligned(t_64, t_65, t_66, ab_z, si_1_8, si_1_9, si_1_10, si_0_8, si_0_9, \
                         si_0_10, sk_1_13, sk_1_14, sk_1_16 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_64[k] = -ab_z[k] * si_1_8[k]
                      + si_0_8[k]
                      + sk_1_13[k];

            t_65[k] = -ab_z[k] * si_1_9[k]
                      + si_0_9[k]
                      + sk_1_14[k];

            t_66[k] = -ab_z[k] * si_1_10[k]
                      + si_0_10[k]
                      + sk_1_16[k];
        }

#pragma omp simd aligned(t_67, t_68, t_69, ab_z, si_1_11, si_1_12, si_1_13, si_0_11, si_0_12, \
                         si_0_13, sk_1_17, sk_1_18, sk_1_19 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_67[k] = -ab_z[k] * si_1_11[k]
                      + si_0_11[k]
                      + sk_1_17[k];

            t_68[k] = -ab_z[k] * si_1_12[k]
                      + si_0_12[k]
                      + sk_1_18[k];

            t_69[k] = -ab_z[k] * si_1_13[k]
                      + si_0_13[k]
                      + sk_1_19[k];
        }

#pragma omp simd aligned(t_70, t_71, t_72, ab_z, si_1_14, si_1_15, si_1_16, si_0_14, si_0_15, \
                         si_0_16, sk_1_20, sk_1_22, sk_1_23 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_70[k] = -ab_z[k] * si_1_14[k]
                      + si_0_14[k]
                      + sk_1_20[k];

            t_71[k] = -ab_z[k] * si_1_15[k]
                      + si_0_15[k]
                      + sk_1_22[k];

            t_72[k] = -ab_z[k] * si_1_16[k]
                      + si_0_16[k]
                      + sk_1_23[k];
        }

#pragma omp simd aligned(t_73, t_74, t_75, ab_z, si_1_17, si_1_18, si_1_19, si_0_17, si_0_18, \
                         si_0_19, sk_1_24, sk_1_25, sk_1_26 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_73[k] = -ab_z[k] * si_1_17[k]
                      + si_0_17[k]
                      + sk_1_24[k];

            t_74[k] = -ab_z[k] * si_1_18[k]
                      + si_0_18[k]
                      + sk_1_25[k];

            t_75[k] = -ab_z[k] * si_1_19[k]
                      + si_0_19[k]
                      + sk_1_26[k];
        }

#pragma omp simd aligned(t_76, t_77, t_78, ab_z, si_1_20, si_1_21, si_1_22, si_0_20, si_0_21, \
                         si_0_22, sk_1_27, sk_1_29, sk_1_30 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_76[k] = -ab_z[k] * si_1_20[k]
                      + si_0_20[k]
                      + sk_1_27[k];

            t_77[k] = -ab_z[k] * si_1_21[k]
                      + si_0_21[k]
                      + sk_1_29[k];

            t_78[k] = -ab_z[k] * si_1_22[k]
                      + si_0_22[k]
                      + sk_1_30[k];
        }

#pragma omp simd aligned(t_79, t_80, t_81, ab_z, si_1_23, si_1_24, si_1_25, si_0_23, si_0_24, \
                         si_0_25, sk_1_31, sk_1_32, sk_1_33 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_79[k] = -ab_z[k] * si_1_23[k]
                      + si_0_23[k]
                      + sk_1_31[k];

            t_80[k] = -ab_z[k] * si_1_24[k]
                      + si_0_24[k]
                      + sk_1_32[k];

            t_81[k] = -ab_z[k] * si_1_25[k]
                      + si_0_25[k]
                      + sk_1_33[k];
        }

#pragma omp simd aligned(t_82, t_83, ab_z, si_1_26, si_1_27, si_0_26, si_0_27, sk_1_34, \
                         sk_1_35 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_82[k] = -ab_z[k] * si_1_26[k]
                      + si_0_26[k]
                      + sk_1_34[k];

            t_83[k] = -ab_z[k] * si_1_27[k]
                      + si_0_27[k]
                      + sk_1_35[k];
        }
    }
}

}  // namespace simdtrf
