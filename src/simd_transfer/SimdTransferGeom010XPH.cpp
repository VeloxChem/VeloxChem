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


#include "SimdTransferGeom010XPH.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_geom_010x_ph(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                         const size_t target, const size_t sh_1, const size_t sh_0,
                         const size_t si_1, const size_t ncomps, const size_t nmax) -> void
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

        const auto *sh_1_0 = buffer.data(sh_1 + 0 * ncomps + c);
        const auto *sh_1_1 = buffer.data(sh_1 + 1 * ncomps + c);
        const auto *sh_1_2 = buffer.data(sh_1 + 2 * ncomps + c);
        const auto *sh_1_3 = buffer.data(sh_1 + 3 * ncomps + c);
        const auto *sh_1_4 = buffer.data(sh_1 + 4 * ncomps + c);
        const auto *sh_1_5 = buffer.data(sh_1 + 5 * ncomps + c);
        const auto *sh_1_6 = buffer.data(sh_1 + 6 * ncomps + c);
        const auto *sh_1_7 = buffer.data(sh_1 + 7 * ncomps + c);
        const auto *sh_1_8 = buffer.data(sh_1 + 8 * ncomps + c);
        const auto *sh_1_9 = buffer.data(sh_1 + 9 * ncomps + c);
        const auto *sh_1_10 = buffer.data(sh_1 + 10 * ncomps + c);
        const auto *sh_1_11 = buffer.data(sh_1 + 11 * ncomps + c);
        const auto *sh_1_12 = buffer.data(sh_1 + 12 * ncomps + c);
        const auto *sh_1_13 = buffer.data(sh_1 + 13 * ncomps + c);
        const auto *sh_1_14 = buffer.data(sh_1 + 14 * ncomps + c);
        const auto *sh_1_15 = buffer.data(sh_1 + 15 * ncomps + c);
        const auto *sh_1_16 = buffer.data(sh_1 + 16 * ncomps + c);
        const auto *sh_1_17 = buffer.data(sh_1 + 17 * ncomps + c);
        const auto *sh_1_18 = buffer.data(sh_1 + 18 * ncomps + c);
        const auto *sh_1_19 = buffer.data(sh_1 + 19 * ncomps + c);
        const auto *sh_1_20 = buffer.data(sh_1 + 20 * ncomps + c);

        const auto *sh_0_0 = buffer.data(sh_0 + 0 * ncomps + c);
        const auto *sh_0_1 = buffer.data(sh_0 + 1 * ncomps + c);
        const auto *sh_0_2 = buffer.data(sh_0 + 2 * ncomps + c);
        const auto *sh_0_3 = buffer.data(sh_0 + 3 * ncomps + c);
        const auto *sh_0_4 = buffer.data(sh_0 + 4 * ncomps + c);
        const auto *sh_0_5 = buffer.data(sh_0 + 5 * ncomps + c);
        const auto *sh_0_6 = buffer.data(sh_0 + 6 * ncomps + c);
        const auto *sh_0_7 = buffer.data(sh_0 + 7 * ncomps + c);
        const auto *sh_0_8 = buffer.data(sh_0 + 8 * ncomps + c);
        const auto *sh_0_9 = buffer.data(sh_0 + 9 * ncomps + c);
        const auto *sh_0_10 = buffer.data(sh_0 + 10 * ncomps + c);
        const auto *sh_0_11 = buffer.data(sh_0 + 11 * ncomps + c);
        const auto *sh_0_12 = buffer.data(sh_0 + 12 * ncomps + c);
        const auto *sh_0_13 = buffer.data(sh_0 + 13 * ncomps + c);
        const auto *sh_0_14 = buffer.data(sh_0 + 14 * ncomps + c);
        const auto *sh_0_15 = buffer.data(sh_0 + 15 * ncomps + c);
        const auto *sh_0_16 = buffer.data(sh_0 + 16 * ncomps + c);
        const auto *sh_0_17 = buffer.data(sh_0 + 17 * ncomps + c);
        const auto *sh_0_18 = buffer.data(sh_0 + 18 * ncomps + c);
        const auto *sh_0_19 = buffer.data(sh_0 + 19 * ncomps + c);
        const auto *sh_0_20 = buffer.data(sh_0 + 20 * ncomps + c);

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

#pragma omp simd aligned(t_0, t_1, t_2, ab_x, sh_1_0, sh_1_1, sh_1_2, sh_0_0, sh_0_1, sh_0_2, \
                         si_1_0, si_1_1, si_1_2 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = -ab_x[k] * sh_1_0[k]
                     + sh_0_0[k]
                     + si_1_0[k];

            t_1[k] = -ab_x[k] * sh_1_1[k]
                     + sh_0_1[k]
                     + si_1_1[k];

            t_2[k] = -ab_x[k] * sh_1_2[k]
                     + sh_0_2[k]
                     + si_1_2[k];
        }

#pragma omp simd aligned(t_3, t_4, t_5, ab_x, sh_1_3, sh_1_4, sh_1_5, sh_0_3, sh_0_4, sh_0_5, \
                         si_1_3, si_1_4, si_1_5 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_3[k] = -ab_x[k] * sh_1_3[k]
                     + sh_0_3[k]
                     + si_1_3[k];

            t_4[k] = -ab_x[k] * sh_1_4[k]
                     + sh_0_4[k]
                     + si_1_4[k];

            t_5[k] = -ab_x[k] * sh_1_5[k]
                     + sh_0_5[k]
                     + si_1_5[k];
        }

#pragma omp simd aligned(t_6, t_7, t_8, ab_x, sh_1_6, sh_1_7, sh_1_8, sh_0_6, sh_0_7, sh_0_8, \
                         si_1_6, si_1_7, si_1_8 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_6[k] = -ab_x[k] * sh_1_6[k]
                     + sh_0_6[k]
                     + si_1_6[k];

            t_7[k] = -ab_x[k] * sh_1_7[k]
                     + sh_0_7[k]
                     + si_1_7[k];

            t_8[k] = -ab_x[k] * sh_1_8[k]
                     + sh_0_8[k]
                     + si_1_8[k];
        }

#pragma omp simd aligned(t_9, t_10, t_11, ab_x, sh_1_9, sh_1_10, sh_1_11, sh_0_9, sh_0_10, \
                         sh_0_11, si_1_9, si_1_10, si_1_11 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_9[k] = -ab_x[k] * sh_1_9[k]
                     + sh_0_9[k]
                     + si_1_9[k];

            t_10[k] = -ab_x[k] * sh_1_10[k]
                      + sh_0_10[k]
                      + si_1_10[k];

            t_11[k] = -ab_x[k] * sh_1_11[k]
                      + sh_0_11[k]
                      + si_1_11[k];
        }

#pragma omp simd aligned(t_12, t_13, t_14, ab_x, sh_1_12, sh_1_13, sh_1_14, sh_0_12, sh_0_13, \
                         sh_0_14, si_1_12, si_1_13, si_1_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_12[k] = -ab_x[k] * sh_1_12[k]
                      + sh_0_12[k]
                      + si_1_12[k];

            t_13[k] = -ab_x[k] * sh_1_13[k]
                      + sh_0_13[k]
                      + si_1_13[k];

            t_14[k] = -ab_x[k] * sh_1_14[k]
                      + sh_0_14[k]
                      + si_1_14[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, ab_x, sh_1_15, sh_1_16, sh_1_17, sh_0_15, sh_0_16, \
                         sh_0_17, si_1_15, si_1_16, si_1_17 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = -ab_x[k] * sh_1_15[k]
                      + sh_0_15[k]
                      + si_1_15[k];

            t_16[k] = -ab_x[k] * sh_1_16[k]
                      + sh_0_16[k]
                      + si_1_16[k];

            t_17[k] = -ab_x[k] * sh_1_17[k]
                      + sh_0_17[k]
                      + si_1_17[k];
        }

#pragma omp simd aligned(t_18, t_19, t_20, ab_x, sh_1_18, sh_1_19, sh_1_20, sh_0_18, sh_0_19, \
                         sh_0_20, si_1_18, si_1_19, si_1_20 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_18[k] = -ab_x[k] * sh_1_18[k]
                      + sh_0_18[k]
                      + si_1_18[k];

            t_19[k] = -ab_x[k] * sh_1_19[k]
                      + sh_0_19[k]
                      + si_1_19[k];

            t_20[k] = -ab_x[k] * sh_1_20[k]
                      + sh_0_20[k]
                      + si_1_20[k];
        }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, ab_y, sh_1_0, sh_1_1, sh_1_2, sh_1_3, \
                         sh_1_4, si_1_1, si_1_3, si_1_4, si_1_6, \
                         si_1_7 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_21[k] = -ab_y[k] * sh_1_0[k]
                      + si_1_1[k];

            t_22[k] = -ab_y[k] * sh_1_1[k]
                      + si_1_3[k];

            t_23[k] = -ab_y[k] * sh_1_2[k]
                      + si_1_4[k];

            t_24[k] = -ab_y[k] * sh_1_3[k]
                      + si_1_6[k];

            t_25[k] = -ab_y[k] * sh_1_4[k]
                      + si_1_7[k];
        }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, ab_y, sh_1_5, sh_1_6, sh_1_7, sh_1_8, \
                         sh_1_9, si_1_8, si_1_10, si_1_11, si_1_12, \
                         si_1_13 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_26[k] = -ab_y[k] * sh_1_5[k]
                      + si_1_8[k];

            t_27[k] = -ab_y[k] * sh_1_6[k]
                      + si_1_10[k];

            t_28[k] = -ab_y[k] * sh_1_7[k]
                      + si_1_11[k];

            t_29[k] = -ab_y[k] * sh_1_8[k]
                      + si_1_12[k];

            t_30[k] = -ab_y[k] * sh_1_9[k]
                      + si_1_13[k];
        }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, ab_y, sh_1_10, sh_1_11, sh_1_12, \
                         sh_1_13, sh_1_14, si_1_15, si_1_16, si_1_17, si_1_18, \
                         si_1_19 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_31[k] = -ab_y[k] * sh_1_10[k]
                      + si_1_15[k];

            t_32[k] = -ab_y[k] * sh_1_11[k]
                      + si_1_16[k];

            t_33[k] = -ab_y[k] * sh_1_12[k]
                      + si_1_17[k];

            t_34[k] = -ab_y[k] * sh_1_13[k]
                      + si_1_18[k];

            t_35[k] = -ab_y[k] * sh_1_14[k]
                      + si_1_19[k];
        }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, ab_y, sh_1_15, sh_1_16, sh_1_17, \
                         sh_1_18, sh_1_19, si_1_21, si_1_22, si_1_23, si_1_24, \
                         si_1_25 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_36[k] = -ab_y[k] * sh_1_15[k]
                      + si_1_21[k];

            t_37[k] = -ab_y[k] * sh_1_16[k]
                      + si_1_22[k];

            t_38[k] = -ab_y[k] * sh_1_17[k]
                      + si_1_23[k];

            t_39[k] = -ab_y[k] * sh_1_18[k]
                      + si_1_24[k];

            t_40[k] = -ab_y[k] * sh_1_19[k]
                      + si_1_25[k];
        }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, ab_y, ab_z, sh_1_0, sh_1_1, sh_1_2, sh_1_20, \
                         si_1_2, si_1_4, si_1_5, si_1_26 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_41[k] = -ab_y[k] * sh_1_20[k]
                      + si_1_26[k];

            t_42[k] = -ab_z[k] * sh_1_0[k]
                      + si_1_2[k];

            t_43[k] = -ab_z[k] * sh_1_1[k]
                      + si_1_4[k];

            t_44[k] = -ab_z[k] * sh_1_2[k]
                      + si_1_5[k];
        }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_z, sh_1_3, sh_1_4, sh_1_5, sh_1_6, \
                         sh_1_7, si_1_7, si_1_8, si_1_9, si_1_11, \
                         si_1_12 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_45[k] = -ab_z[k] * sh_1_3[k]
                      + si_1_7[k];

            t_46[k] = -ab_z[k] * sh_1_4[k]
                      + si_1_8[k];

            t_47[k] = -ab_z[k] * sh_1_5[k]
                      + si_1_9[k];

            t_48[k] = -ab_z[k] * sh_1_6[k]
                      + si_1_11[k];

            t_49[k] = -ab_z[k] * sh_1_7[k]
                      + si_1_12[k];
        }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_z, sh_1_8, sh_1_9, sh_1_10, sh_1_11, \
                         sh_1_12, si_1_13, si_1_14, si_1_16, si_1_17, \
                         si_1_18 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_50[k] = -ab_z[k] * sh_1_8[k]
                      + si_1_13[k];

            t_51[k] = -ab_z[k] * sh_1_9[k]
                      + si_1_14[k];

            t_52[k] = -ab_z[k] * sh_1_10[k]
                      + si_1_16[k];

            t_53[k] = -ab_z[k] * sh_1_11[k]
                      + si_1_17[k];

            t_54[k] = -ab_z[k] * sh_1_12[k]
                      + si_1_18[k];
        }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_z, sh_1_13, sh_1_14, sh_1_15, \
                         sh_1_16, sh_1_17, si_1_19, si_1_20, si_1_22, si_1_23, \
                         si_1_24 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_55[k] = -ab_z[k] * sh_1_13[k]
                      + si_1_19[k];

            t_56[k] = -ab_z[k] * sh_1_14[k]
                      + si_1_20[k];

            t_57[k] = -ab_z[k] * sh_1_15[k]
                      + si_1_22[k];

            t_58[k] = -ab_z[k] * sh_1_16[k]
                      + si_1_23[k];

            t_59[k] = -ab_z[k] * sh_1_17[k]
                      + si_1_24[k];
        }

#pragma omp simd aligned(t_60, t_61, t_62, ab_z, sh_1_18, sh_1_19, sh_1_20, si_1_25, si_1_26, \
                         si_1_27 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_60[k] = -ab_z[k] * sh_1_18[k]
                      + si_1_25[k];

            t_61[k] = -ab_z[k] * sh_1_19[k]
                      + si_1_26[k];

            t_62[k] = -ab_z[k] * sh_1_20[k]
                      + si_1_27[k];
        }
    }
}

}  // namespace simdtrf
