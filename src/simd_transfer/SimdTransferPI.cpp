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


#include "SimdTransferPI.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_pi(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t si, const size_t sk, const size_t ncomps, const size_t nmax) -> void
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

        const auto *si_0 = buffer.data(si + 0 * ncomps + c);
        const auto *si_1 = buffer.data(si + 1 * ncomps + c);
        const auto *si_2 = buffer.data(si + 2 * ncomps + c);
        const auto *si_3 = buffer.data(si + 3 * ncomps + c);
        const auto *si_4 = buffer.data(si + 4 * ncomps + c);
        const auto *si_5 = buffer.data(si + 5 * ncomps + c);
        const auto *si_6 = buffer.data(si + 6 * ncomps + c);
        const auto *si_7 = buffer.data(si + 7 * ncomps + c);
        const auto *si_8 = buffer.data(si + 8 * ncomps + c);
        const auto *si_9 = buffer.data(si + 9 * ncomps + c);
        const auto *si_10 = buffer.data(si + 10 * ncomps + c);
        const auto *si_11 = buffer.data(si + 11 * ncomps + c);
        const auto *si_12 = buffer.data(si + 12 * ncomps + c);
        const auto *si_13 = buffer.data(si + 13 * ncomps + c);
        const auto *si_14 = buffer.data(si + 14 * ncomps + c);
        const auto *si_15 = buffer.data(si + 15 * ncomps + c);
        const auto *si_16 = buffer.data(si + 16 * ncomps + c);
        const auto *si_17 = buffer.data(si + 17 * ncomps + c);
        const auto *si_18 = buffer.data(si + 18 * ncomps + c);
        const auto *si_19 = buffer.data(si + 19 * ncomps + c);
        const auto *si_20 = buffer.data(si + 20 * ncomps + c);
        const auto *si_21 = buffer.data(si + 21 * ncomps + c);
        const auto *si_22 = buffer.data(si + 22 * ncomps + c);
        const auto *si_23 = buffer.data(si + 23 * ncomps + c);
        const auto *si_24 = buffer.data(si + 24 * ncomps + c);
        const auto *si_25 = buffer.data(si + 25 * ncomps + c);
        const auto *si_26 = buffer.data(si + 26 * ncomps + c);
        const auto *si_27 = buffer.data(si + 27 * ncomps + c);

        const auto *sk_0 = buffer.data(sk + 0 * ncomps + c);
        const auto *sk_1 = buffer.data(sk + 1 * ncomps + c);
        const auto *sk_2 = buffer.data(sk + 2 * ncomps + c);
        const auto *sk_3 = buffer.data(sk + 3 * ncomps + c);
        const auto *sk_4 = buffer.data(sk + 4 * ncomps + c);
        const auto *sk_5 = buffer.data(sk + 5 * ncomps + c);
        const auto *sk_6 = buffer.data(sk + 6 * ncomps + c);
        const auto *sk_7 = buffer.data(sk + 7 * ncomps + c);
        const auto *sk_8 = buffer.data(sk + 8 * ncomps + c);
        const auto *sk_9 = buffer.data(sk + 9 * ncomps + c);
        const auto *sk_10 = buffer.data(sk + 10 * ncomps + c);
        const auto *sk_11 = buffer.data(sk + 11 * ncomps + c);
        const auto *sk_12 = buffer.data(sk + 12 * ncomps + c);
        const auto *sk_13 = buffer.data(sk + 13 * ncomps + c);
        const auto *sk_14 = buffer.data(sk + 14 * ncomps + c);
        const auto *sk_15 = buffer.data(sk + 15 * ncomps + c);
        const auto *sk_16 = buffer.data(sk + 16 * ncomps + c);
        const auto *sk_17 = buffer.data(sk + 17 * ncomps + c);
        const auto *sk_18 = buffer.data(sk + 18 * ncomps + c);
        const auto *sk_19 = buffer.data(sk + 19 * ncomps + c);
        const auto *sk_20 = buffer.data(sk + 20 * ncomps + c);
        const auto *sk_21 = buffer.data(sk + 21 * ncomps + c);
        const auto *sk_22 = buffer.data(sk + 22 * ncomps + c);
        const auto *sk_23 = buffer.data(sk + 23 * ncomps + c);
        const auto *sk_24 = buffer.data(sk + 24 * ncomps + c);
        const auto *sk_25 = buffer.data(sk + 25 * ncomps + c);
        const auto *sk_26 = buffer.data(sk + 26 * ncomps + c);
        const auto *sk_27 = buffer.data(sk + 27 * ncomps + c);
        const auto *sk_28 = buffer.data(sk + 28 * ncomps + c);
        const auto *sk_29 = buffer.data(sk + 29 * ncomps + c);
        const auto *sk_30 = buffer.data(sk + 30 * ncomps + c);
        const auto *sk_31 = buffer.data(sk + 31 * ncomps + c);
        const auto *sk_32 = buffer.data(sk + 32 * ncomps + c);
        const auto *sk_33 = buffer.data(sk + 33 * ncomps + c);
        const auto *sk_34 = buffer.data(sk + 34 * ncomps + c);
        const auto *sk_35 = buffer.data(sk + 35 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, si_0, si_1, si_2, si_3, si_4, sk_0, \
                         sk_1, sk_2, sk_3, sk_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = -ab_x[k] * si_0[k]
                     + sk_0[k];

            t_1[k] = -ab_x[k] * si_1[k]
                     + sk_1[k];

            t_2[k] = -ab_x[k] * si_2[k]
                     + sk_2[k];

            t_3[k] = -ab_x[k] * si_3[k]
                     + sk_3[k];

            t_4[k] = -ab_x[k] * si_4[k]
                     + sk_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, si_5, si_6, si_7, si_8, si_9, sk_5, \
                         sk_6, sk_7, sk_8, sk_9 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = -ab_x[k] * si_5[k]
                     + sk_5[k];

            t_6[k] = -ab_x[k] * si_6[k]
                     + sk_6[k];

            t_7[k] = -ab_x[k] * si_7[k]
                     + sk_7[k];

            t_8[k] = -ab_x[k] * si_8[k]
                     + sk_8[k];

            t_9[k] = -ab_x[k] * si_9[k]
                     + sk_9[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, si_10, si_11, si_12, si_13, \
                         si_14, sk_10, sk_11, sk_12, sk_13, sk_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = -ab_x[k] * si_10[k]
                      + sk_10[k];

            t_11[k] = -ab_x[k] * si_11[k]
                      + sk_11[k];

            t_12[k] = -ab_x[k] * si_12[k]
                      + sk_12[k];

            t_13[k] = -ab_x[k] * si_13[k]
                      + sk_13[k];

            t_14[k] = -ab_x[k] * si_14[k]
                      + sk_14[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, si_15, si_16, si_17, si_18, \
                         si_19, sk_15, sk_16, sk_17, sk_18, sk_19 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = -ab_x[k] * si_15[k]
                      + sk_15[k];

            t_16[k] = -ab_x[k] * si_16[k]
                      + sk_16[k];

            t_17[k] = -ab_x[k] * si_17[k]
                      + sk_17[k];

            t_18[k] = -ab_x[k] * si_18[k]
                      + sk_18[k];

            t_19[k] = -ab_x[k] * si_19[k]
                      + sk_19[k];
        }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, si_20, si_21, si_22, si_23, \
                         si_24, sk_20, sk_21, sk_22, sk_23, sk_24 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_20[k] = -ab_x[k] * si_20[k]
                      + sk_20[k];

            t_21[k] = -ab_x[k] * si_21[k]
                      + sk_21[k];

            t_22[k] = -ab_x[k] * si_22[k]
                      + sk_22[k];

            t_23[k] = -ab_x[k] * si_23[k]
                      + sk_23[k];

            t_24[k] = -ab_x[k] * si_24[k]
                      + sk_24[k];
        }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, ab_x, ab_y, si_0, si_25, si_26, si_27, sk_1, \
                         sk_25, sk_26, sk_27 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_25[k] = -ab_x[k] * si_25[k]
                      + sk_25[k];

            t_26[k] = -ab_x[k] * si_26[k]
                      + sk_26[k];

            t_27[k] = -ab_x[k] * si_27[k]
                      + sk_27[k];

            t_28[k] = -ab_y[k] * si_0[k]
                      + sk_1[k];
        }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, ab_y, si_1, si_2, si_3, si_4, si_5, \
                         sk_3, sk_4, sk_6, sk_7, sk_8 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_29[k] = -ab_y[k] * si_1[k]
                      + sk_3[k];

            t_30[k] = -ab_y[k] * si_2[k]
                      + sk_4[k];

            t_31[k] = -ab_y[k] * si_3[k]
                      + sk_6[k];

            t_32[k] = -ab_y[k] * si_4[k]
                      + sk_7[k];

            t_33[k] = -ab_y[k] * si_5[k]
                      + sk_8[k];
        }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, t_38, ab_y, si_6, si_7, si_8, si_9, si_10, \
                         sk_10, sk_11, sk_12, sk_13, sk_15 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_34[k] = -ab_y[k] * si_6[k]
                      + sk_10[k];

            t_35[k] = -ab_y[k] * si_7[k]
                      + sk_11[k];

            t_36[k] = -ab_y[k] * si_8[k]
                      + sk_12[k];

            t_37[k] = -ab_y[k] * si_9[k]
                      + sk_13[k];

            t_38[k] = -ab_y[k] * si_10[k]
                      + sk_15[k];
        }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, t_43, ab_y, si_11, si_12, si_13, si_14, \
                         si_15, sk_16, sk_17, sk_18, sk_19, sk_21 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_39[k] = -ab_y[k] * si_11[k]
                      + sk_16[k];

            t_40[k] = -ab_y[k] * si_12[k]
                      + sk_17[k];

            t_41[k] = -ab_y[k] * si_13[k]
                      + sk_18[k];

            t_42[k] = -ab_y[k] * si_14[k]
                      + sk_19[k];

            t_43[k] = -ab_y[k] * si_15[k]
                      + sk_21[k];
        }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, ab_y, si_16, si_17, si_18, si_19, \
                         si_20, sk_22, sk_23, sk_24, sk_25, sk_26 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_44[k] = -ab_y[k] * si_16[k]
                      + sk_22[k];

            t_45[k] = -ab_y[k] * si_17[k]
                      + sk_23[k];

            t_46[k] = -ab_y[k] * si_18[k]
                      + sk_24[k];

            t_47[k] = -ab_y[k] * si_19[k]
                      + sk_25[k];

            t_48[k] = -ab_y[k] * si_20[k]
                      + sk_26[k];
        }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, ab_y, si_21, si_22, si_23, si_24, \
                         si_25, sk_28, sk_29, sk_30, sk_31, sk_32 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_49[k] = -ab_y[k] * si_21[k]
                      + sk_28[k];

            t_50[k] = -ab_y[k] * si_22[k]
                      + sk_29[k];

            t_51[k] = -ab_y[k] * si_23[k]
                      + sk_30[k];

            t_52[k] = -ab_y[k] * si_24[k]
                      + sk_31[k];

            t_53[k] = -ab_y[k] * si_25[k]
                      + sk_32[k];
        }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, ab_y, ab_z, si_0, si_1, si_26, si_27, sk_2, \
                         sk_4, sk_33, sk_34 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_54[k] = -ab_y[k] * si_26[k]
                      + sk_33[k];

            t_55[k] = -ab_y[k] * si_27[k]
                      + sk_34[k];

            t_56[k] = -ab_z[k] * si_0[k]
                      + sk_2[k];

            t_57[k] = -ab_z[k] * si_1[k]
                      + sk_4[k];
        }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, ab_z, si_2, si_3, si_4, si_5, si_6, \
                         sk_5, sk_7, sk_8, sk_9, sk_11 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_58[k] = -ab_z[k] * si_2[k]
                      + sk_5[k];

            t_59[k] = -ab_z[k] * si_3[k]
                      + sk_7[k];

            t_60[k] = -ab_z[k] * si_4[k]
                      + sk_8[k];

            t_61[k] = -ab_z[k] * si_5[k]
                      + sk_9[k];

            t_62[k] = -ab_z[k] * si_6[k]
                      + sk_11[k];
        }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, ab_z, si_7, si_8, si_9, si_10, si_11, \
                         sk_12, sk_13, sk_14, sk_16, sk_17 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_63[k] = -ab_z[k] * si_7[k]
                      + sk_12[k];

            t_64[k] = -ab_z[k] * si_8[k]
                      + sk_13[k];

            t_65[k] = -ab_z[k] * si_9[k]
                      + sk_14[k];

            t_66[k] = -ab_z[k] * si_10[k]
                      + sk_16[k];

            t_67[k] = -ab_z[k] * si_11[k]
                      + sk_17[k];
        }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, ab_z, si_12, si_13, si_14, si_15, \
                         si_16, sk_18, sk_19, sk_20, sk_22, sk_23 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_68[k] = -ab_z[k] * si_12[k]
                      + sk_18[k];

            t_69[k] = -ab_z[k] * si_13[k]
                      + sk_19[k];

            t_70[k] = -ab_z[k] * si_14[k]
                      + sk_20[k];

            t_71[k] = -ab_z[k] * si_15[k]
                      + sk_22[k];

            t_72[k] = -ab_z[k] * si_16[k]
                      + sk_23[k];
        }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, ab_z, si_17, si_18, si_19, si_20, \
                         si_21, sk_24, sk_25, sk_26, sk_27, sk_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_73[k] = -ab_z[k] * si_17[k]
                      + sk_24[k];

            t_74[k] = -ab_z[k] * si_18[k]
                      + sk_25[k];

            t_75[k] = -ab_z[k] * si_19[k]
                      + sk_26[k];

            t_76[k] = -ab_z[k] * si_20[k]
                      + sk_27[k];

            t_77[k] = -ab_z[k] * si_21[k]
                      + sk_29[k];
        }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, ab_z, si_22, si_23, si_24, si_25, \
                         si_26, sk_30, sk_31, sk_32, sk_33, sk_34 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_78[k] = -ab_z[k] * si_22[k]
                      + sk_30[k];

            t_79[k] = -ab_z[k] * si_23[k]
                      + sk_31[k];

            t_80[k] = -ab_z[k] * si_24[k]
                      + sk_32[k];

            t_81[k] = -ab_z[k] * si_25[k]
                      + sk_33[k];

            t_82[k] = -ab_z[k] * si_26[k]
                      + sk_34[k];
        }

#pragma omp simd aligned(t_83, ab_z, si_27, sk_35 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_83[k] = -ab_z[k] * si_27[k]
                      + sk_35[k];
        }
    }
}

}  // namespace simdtrf
