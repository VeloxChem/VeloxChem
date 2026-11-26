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


#include "SimdTransferGeom100XFD.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_geom_100x_fd_out_of_first(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                      const size_t target, const size_t fp_1, const size_t fp_0,
                                      const size_t gp_1, const size_t ncomps,
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

        const auto *fp_1_0 = buffer.data(fp_1 + 0 * ncomps + c);
        const auto *fp_1_1 = buffer.data(fp_1 + 1 * ncomps + c);
        const auto *fp_1_2 = buffer.data(fp_1 + 2 * ncomps + c);
        const auto *fp_1_3 = buffer.data(fp_1 + 3 * ncomps + c);
        const auto *fp_1_4 = buffer.data(fp_1 + 4 * ncomps + c);
        const auto *fp_1_5 = buffer.data(fp_1 + 5 * ncomps + c);
        const auto *fp_1_6 = buffer.data(fp_1 + 6 * ncomps + c);
        const auto *fp_1_7 = buffer.data(fp_1 + 7 * ncomps + c);
        const auto *fp_1_8 = buffer.data(fp_1 + 8 * ncomps + c);
        const auto *fp_1_9 = buffer.data(fp_1 + 9 * ncomps + c);
        const auto *fp_1_10 = buffer.data(fp_1 + 10 * ncomps + c);
        const auto *fp_1_11 = buffer.data(fp_1 + 11 * ncomps + c);
        const auto *fp_1_12 = buffer.data(fp_1 + 12 * ncomps + c);
        const auto *fp_1_13 = buffer.data(fp_1 + 13 * ncomps + c);
        const auto *fp_1_14 = buffer.data(fp_1 + 14 * ncomps + c);
        const auto *fp_1_15 = buffer.data(fp_1 + 15 * ncomps + c);
        const auto *fp_1_16 = buffer.data(fp_1 + 16 * ncomps + c);
        const auto *fp_1_17 = buffer.data(fp_1 + 17 * ncomps + c);
        const auto *fp_1_18 = buffer.data(fp_1 + 18 * ncomps + c);
        const auto *fp_1_19 = buffer.data(fp_1 + 19 * ncomps + c);
        const auto *fp_1_20 = buffer.data(fp_1 + 20 * ncomps + c);
        const auto *fp_1_21 = buffer.data(fp_1 + 21 * ncomps + c);
        const auto *fp_1_22 = buffer.data(fp_1 + 22 * ncomps + c);
        const auto *fp_1_23 = buffer.data(fp_1 + 23 * ncomps + c);
        const auto *fp_1_24 = buffer.data(fp_1 + 24 * ncomps + c);
        const auto *fp_1_25 = buffer.data(fp_1 + 25 * ncomps + c);
        const auto *fp_1_26 = buffer.data(fp_1 + 26 * ncomps + c);
        const auto *fp_1_27 = buffer.data(fp_1 + 27 * ncomps + c);
        const auto *fp_1_28 = buffer.data(fp_1 + 28 * ncomps + c);
        const auto *fp_1_29 = buffer.data(fp_1 + 29 * ncomps + c);

        const auto *fp_0_0 = buffer.data(fp_0 + 0 * ncomps + c);
        const auto *fp_0_1 = buffer.data(fp_0 + 1 * ncomps + c);
        const auto *fp_0_2 = buffer.data(fp_0 + 2 * ncomps + c);
        const auto *fp_0_3 = buffer.data(fp_0 + 3 * ncomps + c);
        const auto *fp_0_4 = buffer.data(fp_0 + 4 * ncomps + c);
        const auto *fp_0_5 = buffer.data(fp_0 + 5 * ncomps + c);
        const auto *fp_0_6 = buffer.data(fp_0 + 6 * ncomps + c);
        const auto *fp_0_7 = buffer.data(fp_0 + 7 * ncomps + c);
        const auto *fp_0_8 = buffer.data(fp_0 + 8 * ncomps + c);
        const auto *fp_0_9 = buffer.data(fp_0 + 9 * ncomps + c);
        const auto *fp_0_10 = buffer.data(fp_0 + 10 * ncomps + c);
        const auto *fp_0_11 = buffer.data(fp_0 + 11 * ncomps + c);
        const auto *fp_0_12 = buffer.data(fp_0 + 12 * ncomps + c);
        const auto *fp_0_13 = buffer.data(fp_0 + 13 * ncomps + c);
        const auto *fp_0_14 = buffer.data(fp_0 + 14 * ncomps + c);
        const auto *fp_0_15 = buffer.data(fp_0 + 15 * ncomps + c);
        const auto *fp_0_16 = buffer.data(fp_0 + 16 * ncomps + c);
        const auto *fp_0_17 = buffer.data(fp_0 + 17 * ncomps + c);
        const auto *fp_0_18 = buffer.data(fp_0 + 18 * ncomps + c);
        const auto *fp_0_19 = buffer.data(fp_0 + 19 * ncomps + c);
        const auto *fp_0_20 = buffer.data(fp_0 + 20 * ncomps + c);
        const auto *fp_0_21 = buffer.data(fp_0 + 21 * ncomps + c);
        const auto *fp_0_22 = buffer.data(fp_0 + 22 * ncomps + c);
        const auto *fp_0_23 = buffer.data(fp_0 + 23 * ncomps + c);
        const auto *fp_0_24 = buffer.data(fp_0 + 24 * ncomps + c);
        const auto *fp_0_25 = buffer.data(fp_0 + 25 * ncomps + c);
        const auto *fp_0_26 = buffer.data(fp_0 + 26 * ncomps + c);
        const auto *fp_0_27 = buffer.data(fp_0 + 27 * ncomps + c);
        const auto *fp_0_28 = buffer.data(fp_0 + 28 * ncomps + c);
        const auto *fp_0_29 = buffer.data(fp_0 + 29 * ncomps + c);

        const auto *gp_1_0 = buffer.data(gp_1 + 0 * ncomps + c);
        const auto *gp_1_1 = buffer.data(gp_1 + 1 * ncomps + c);
        const auto *gp_1_2 = buffer.data(gp_1 + 2 * ncomps + c);
        const auto *gp_1_3 = buffer.data(gp_1 + 3 * ncomps + c);
        const auto *gp_1_4 = buffer.data(gp_1 + 4 * ncomps + c);
        const auto *gp_1_5 = buffer.data(gp_1 + 5 * ncomps + c);
        const auto *gp_1_6 = buffer.data(gp_1 + 6 * ncomps + c);
        const auto *gp_1_7 = buffer.data(gp_1 + 7 * ncomps + c);
        const auto *gp_1_8 = buffer.data(gp_1 + 8 * ncomps + c);
        const auto *gp_1_9 = buffer.data(gp_1 + 9 * ncomps + c);
        const auto *gp_1_10 = buffer.data(gp_1 + 10 * ncomps + c);
        const auto *gp_1_11 = buffer.data(gp_1 + 11 * ncomps + c);
        const auto *gp_1_12 = buffer.data(gp_1 + 12 * ncomps + c);
        const auto *gp_1_13 = buffer.data(gp_1 + 13 * ncomps + c);
        const auto *gp_1_14 = buffer.data(gp_1 + 14 * ncomps + c);
        const auto *gp_1_15 = buffer.data(gp_1 + 15 * ncomps + c);
        const auto *gp_1_16 = buffer.data(gp_1 + 16 * ncomps + c);
        const auto *gp_1_17 = buffer.data(gp_1 + 17 * ncomps + c);
        const auto *gp_1_18 = buffer.data(gp_1 + 18 * ncomps + c);
        const auto *gp_1_19 = buffer.data(gp_1 + 19 * ncomps + c);
        const auto *gp_1_20 = buffer.data(gp_1 + 20 * ncomps + c);
        const auto *gp_1_21 = buffer.data(gp_1 + 21 * ncomps + c);
        const auto *gp_1_22 = buffer.data(gp_1 + 22 * ncomps + c);
        const auto *gp_1_23 = buffer.data(gp_1 + 23 * ncomps + c);
        const auto *gp_1_24 = buffer.data(gp_1 + 24 * ncomps + c);
        const auto *gp_1_25 = buffer.data(gp_1 + 25 * ncomps + c);
        const auto *gp_1_26 = buffer.data(gp_1 + 26 * ncomps + c);
        const auto *gp_1_27 = buffer.data(gp_1 + 27 * ncomps + c);
        const auto *gp_1_28 = buffer.data(gp_1 + 28 * ncomps + c);
        const auto *gp_1_29 = buffer.data(gp_1 + 29 * ncomps + c);
        const auto *gp_1_31 = buffer.data(gp_1 + 31 * ncomps + c);
        const auto *gp_1_32 = buffer.data(gp_1 + 32 * ncomps + c);
        const auto *gp_1_34 = buffer.data(gp_1 + 34 * ncomps + c);
        const auto *gp_1_35 = buffer.data(gp_1 + 35 * ncomps + c);
        const auto *gp_1_37 = buffer.data(gp_1 + 37 * ncomps + c);
        const auto *gp_1_38 = buffer.data(gp_1 + 38 * ncomps + c);
        const auto *gp_1_40 = buffer.data(gp_1 + 40 * ncomps + c);
        const auto *gp_1_41 = buffer.data(gp_1 + 41 * ncomps + c);
        const auto *gp_1_44 = buffer.data(gp_1 + 44 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, ab_x, ab_y, fp_1_0, fp_1_1, fp_1_2, fp_0_0, \
                         fp_0_1, fp_0_2, gp_1_0, gp_1_1, gp_1_2, \
                         gp_1_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * fp_1_0[k]
                     + fp_0_0[k]
                     + gp_1_0[k];

            t_1[k] = ab_x[k] * fp_1_1[k]
                     + fp_0_1[k]
                     + gp_1_1[k];

            t_2[k] = ab_x[k] * fp_1_2[k]
                     + fp_0_2[k]
                     + gp_1_2[k];

            t_3[k] = ab_y[k] * fp_1_1[k]
                     + gp_1_4[k];
        }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, ab_x, ab_y, ab_z, fp_1_2, fp_1_3, fp_1_4, fp_0_3, \
                         fp_0_4, gp_1_3, gp_1_4, gp_1_5, gp_1_8 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_4[k] = ab_y[k] * fp_1_2[k]
                     + gp_1_5[k];

            t_5[k] = ab_z[k] * fp_1_2[k]
                     + gp_1_8[k];

            t_6[k] = ab_x[k] * fp_1_3[k]
                     + fp_0_3[k]
                     + gp_1_3[k];

            t_7[k] = ab_x[k] * fp_1_4[k]
                     + fp_0_4[k]
                     + gp_1_4[k];
        }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, ab_x, ab_y, ab_z, fp_1_4, fp_1_5, fp_0_5, \
                         gp_1_5, gp_1_10, gp_1_11, gp_1_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_8[k] = ab_x[k] * fp_1_5[k]
                     + fp_0_5[k]
                     + gp_1_5[k];

            t_9[k] = ab_y[k] * fp_1_4[k]
                     + gp_1_10[k];

            t_10[k] = ab_y[k] * fp_1_5[k]
                      + gp_1_11[k];

            t_11[k] = ab_z[k] * fp_1_5[k]
                      + gp_1_14[k];
        }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, ab_x, ab_y, fp_1_6, fp_1_7, fp_1_8, fp_0_6, \
                         fp_0_7, fp_0_8, gp_1_6, gp_1_7, gp_1_8, \
                         gp_1_13 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_12[k] = ab_x[k] * fp_1_6[k]
                      + fp_0_6[k]
                      + gp_1_6[k];

            t_13[k] = ab_x[k] * fp_1_7[k]
                      + fp_0_7[k]
                      + gp_1_7[k];

            t_14[k] = ab_x[k] * fp_1_8[k]
                      + fp_0_8[k]
                      + gp_1_8[k];

            t_15[k] = ab_y[k] * fp_1_7[k]
                      + gp_1_13[k];
        }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, ab_x, ab_y, ab_z, fp_1_8, fp_1_9, fp_1_10, \
                         fp_0_9, fp_0_10, gp_1_9, gp_1_10, gp_1_14, \
                         gp_1_17 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_16[k] = ab_y[k] * fp_1_8[k]
                      + gp_1_14[k];

            t_17[k] = ab_z[k] * fp_1_8[k]
                      + gp_1_17[k];

            t_18[k] = ab_x[k] * fp_1_9[k]
                      + fp_0_9[k]
                      + gp_1_9[k];

            t_19[k] = ab_x[k] * fp_1_10[k]
                      + fp_0_10[k]
                      + gp_1_10[k];
        }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, ab_x, ab_y, ab_z, fp_1_10, fp_1_11, fp_0_11, \
                         gp_1_11, gp_1_19, gp_1_20, gp_1_23 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_20[k] = ab_x[k] * fp_1_11[k]
                      + fp_0_11[k]
                      + gp_1_11[k];

            t_21[k] = ab_y[k] * fp_1_10[k]
                      + gp_1_19[k];

            t_22[k] = ab_y[k] * fp_1_11[k]
                      + gp_1_20[k];

            t_23[k] = ab_z[k] * fp_1_11[k]
                      + gp_1_23[k];
        }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, ab_x, ab_y, fp_1_12, fp_1_13, fp_1_14, \
                         fp_0_12, fp_0_13, fp_0_14, gp_1_12, gp_1_13, gp_1_14, \
                         gp_1_22 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_24[k] = ab_x[k] * fp_1_12[k]
                      + fp_0_12[k]
                      + gp_1_12[k];

            t_25[k] = ab_x[k] * fp_1_13[k]
                      + fp_0_13[k]
                      + gp_1_13[k];

            t_26[k] = ab_x[k] * fp_1_14[k]
                      + fp_0_14[k]
                      + gp_1_14[k];

            t_27[k] = ab_y[k] * fp_1_13[k]
                      + gp_1_22[k];
        }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, ab_x, ab_y, ab_z, fp_1_14, fp_1_15, fp_1_16, \
                         fp_0_15, fp_0_16, gp_1_15, gp_1_16, gp_1_23, \
                         gp_1_26 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_28[k] = ab_y[k] * fp_1_14[k]
                      + gp_1_23[k];

            t_29[k] = ab_z[k] * fp_1_14[k]
                      + gp_1_26[k];

            t_30[k] = ab_x[k] * fp_1_15[k]
                      + fp_0_15[k]
                      + gp_1_15[k];

            t_31[k] = ab_x[k] * fp_1_16[k]
                      + fp_0_16[k]
                      + gp_1_16[k];
        }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, ab_x, ab_y, ab_z, fp_1_16, fp_1_17, fp_0_17, \
                         gp_1_17, gp_1_25, gp_1_26, gp_1_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_32[k] = ab_x[k] * fp_1_17[k]
                      + fp_0_17[k]
                      + gp_1_17[k];

            t_33[k] = ab_y[k] * fp_1_16[k]
                      + gp_1_25[k];

            t_34[k] = ab_y[k] * fp_1_17[k]
                      + gp_1_26[k];

            t_35[k] = ab_z[k] * fp_1_17[k]
                      + gp_1_29[k];
        }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, ab_x, ab_y, fp_1_18, fp_1_19, fp_1_20, \
                         fp_0_18, fp_0_19, fp_0_20, gp_1_18, gp_1_19, gp_1_20, \
                         gp_1_31 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_36[k] = ab_x[k] * fp_1_18[k]
                      + fp_0_18[k]
                      + gp_1_18[k];

            t_37[k] = ab_x[k] * fp_1_19[k]
                      + fp_0_19[k]
                      + gp_1_19[k];

            t_38[k] = ab_x[k] * fp_1_20[k]
                      + fp_0_20[k]
                      + gp_1_20[k];

            t_39[k] = ab_y[k] * fp_1_19[k]
                      + gp_1_31[k];
        }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, ab_x, ab_y, ab_z, fp_1_20, fp_1_21, fp_1_22, \
                         fp_0_21, fp_0_22, gp_1_21, gp_1_22, gp_1_32, \
                         gp_1_35 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_40[k] = ab_y[k] * fp_1_20[k]
                      + gp_1_32[k];

            t_41[k] = ab_z[k] * fp_1_20[k]
                      + gp_1_35[k];

            t_42[k] = ab_x[k] * fp_1_21[k]
                      + fp_0_21[k]
                      + gp_1_21[k];

            t_43[k] = ab_x[k] * fp_1_22[k]
                      + fp_0_22[k]
                      + gp_1_22[k];
        }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, ab_x, ab_y, ab_z, fp_1_22, fp_1_23, fp_0_23, \
                         gp_1_23, gp_1_34, gp_1_35, gp_1_38 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_44[k] = ab_x[k] * fp_1_23[k]
                      + fp_0_23[k]
                      + gp_1_23[k];

            t_45[k] = ab_y[k] * fp_1_22[k]
                      + gp_1_34[k];

            t_46[k] = ab_y[k] * fp_1_23[k]
                      + gp_1_35[k];

            t_47[k] = ab_z[k] * fp_1_23[k]
                      + gp_1_38[k];
        }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, ab_x, ab_y, fp_1_24, fp_1_25, fp_1_26, \
                         fp_0_24, fp_0_25, fp_0_26, gp_1_24, gp_1_25, gp_1_26, \
                         gp_1_37 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_48[k] = ab_x[k] * fp_1_24[k]
                      + fp_0_24[k]
                      + gp_1_24[k];

            t_49[k] = ab_x[k] * fp_1_25[k]
                      + fp_0_25[k]
                      + gp_1_25[k];

            t_50[k] = ab_x[k] * fp_1_26[k]
                      + fp_0_26[k]
                      + gp_1_26[k];

            t_51[k] = ab_y[k] * fp_1_25[k]
                      + gp_1_37[k];
        }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, ab_x, ab_y, ab_z, fp_1_26, fp_1_27, fp_1_28, \
                         fp_0_27, fp_0_28, gp_1_27, gp_1_28, gp_1_38, \
                         gp_1_41 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_52[k] = ab_y[k] * fp_1_26[k]
                      + gp_1_38[k];

            t_53[k] = ab_z[k] * fp_1_26[k]
                      + gp_1_41[k];

            t_54[k] = ab_x[k] * fp_1_27[k]
                      + fp_0_27[k]
                      + gp_1_27[k];

            t_55[k] = ab_x[k] * fp_1_28[k]
                      + fp_0_28[k]
                      + gp_1_28[k];
        }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, ab_x, ab_y, ab_z, fp_1_28, fp_1_29, fp_0_29, \
                         gp_1_29, gp_1_40, gp_1_41, gp_1_44 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_56[k] = ab_x[k] * fp_1_29[k]
                      + fp_0_29[k]
                      + gp_1_29[k];

            t_57[k] = ab_y[k] * fp_1_28[k]
                      + gp_1_40[k];

            t_58[k] = ab_y[k] * fp_1_29[k]
                      + gp_1_41[k];

            t_59[k] = ab_z[k] * fp_1_29[k]
                      + gp_1_44[k];
        }
    }
}

}  // namespace simdtrf
