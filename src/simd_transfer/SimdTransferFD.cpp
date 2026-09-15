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


#include "SimdTransferFD.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_fd_out_of_first(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                            const size_t target, const size_t fp, const size_t gp,
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

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

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
        const auto *fp_18 = buffer.data(fp + 18 * ncomps + c);
        const auto *fp_19 = buffer.data(fp + 19 * ncomps + c);
        const auto *fp_20 = buffer.data(fp + 20 * ncomps + c);
        const auto *fp_21 = buffer.data(fp + 21 * ncomps + c);
        const auto *fp_22 = buffer.data(fp + 22 * ncomps + c);
        const auto *fp_23 = buffer.data(fp + 23 * ncomps + c);
        const auto *fp_24 = buffer.data(fp + 24 * ncomps + c);
        const auto *fp_25 = buffer.data(fp + 25 * ncomps + c);
        const auto *fp_26 = buffer.data(fp + 26 * ncomps + c);
        const auto *fp_27 = buffer.data(fp + 27 * ncomps + c);
        const auto *fp_28 = buffer.data(fp + 28 * ncomps + c);
        const auto *fp_29 = buffer.data(fp + 29 * ncomps + c);

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
        const auto *gp_31 = buffer.data(gp + 31 * ncomps + c);
        const auto *gp_32 = buffer.data(gp + 32 * ncomps + c);
        const auto *gp_34 = buffer.data(gp + 34 * ncomps + c);
        const auto *gp_35 = buffer.data(gp + 35 * ncomps + c);
        const auto *gp_37 = buffer.data(gp + 37 * ncomps + c);
        const auto *gp_38 = buffer.data(gp + 38 * ncomps + c);
        const auto *gp_40 = buffer.data(gp + 40 * ncomps + c);
        const auto *gp_41 = buffer.data(gp + 41 * ncomps + c);
        const auto *gp_44 = buffer.data(gp + 44 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, ab_y, fp_0, fp_1, fp_2, gp_0, gp_1, \
                         gp_2, gp_4, gp_5 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * fp_0[k]
                     + gp_0[k];

            t_1[k] = ab_x[k] * fp_1[k]
                     + gp_1[k];

            t_2[k] = ab_x[k] * fp_2[k]
                     + gp_2[k];

            t_3[k] = ab_y[k] * fp_1[k]
                     + gp_4[k];

            t_4[k] = ab_y[k] * fp_2[k]
                     + gp_5[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, ab_x, ab_z, fp_2, fp_3, fp_4, fp_5, gp_3, gp_4, \
                         gp_5, gp_8 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = ab_z[k] * fp_2[k]
                     + gp_8[k];

            t_6[k] = ab_x[k] * fp_3[k]
                     + gp_3[k];

            t_7[k] = ab_x[k] * fp_4[k]
                     + gp_4[k];

            t_8[k] = ab_x[k] * fp_5[k]
                     + gp_5[k];
        }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, ab_x, ab_y, ab_z, fp_4, fp_5, fp_6, gp_6, \
                         gp_10, gp_11, gp_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_9[k] = ab_y[k] * fp_4[k]
                     + gp_10[k];

            t_10[k] = ab_y[k] * fp_5[k]
                      + gp_11[k];

            t_11[k] = ab_z[k] * fp_5[k]
                      + gp_14[k];

            t_12[k] = ab_x[k] * fp_6[k]
                      + gp_6[k];
        }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, ab_x, ab_y, ab_z, fp_7, fp_8, gp_7, \
                         gp_8, gp_13, gp_14, gp_17 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_13[k] = ab_x[k] * fp_7[k]
                      + gp_7[k];

            t_14[k] = ab_x[k] * fp_8[k]
                      + gp_8[k];

            t_15[k] = ab_y[k] * fp_7[k]
                      + gp_13[k];

            t_16[k] = ab_y[k] * fp_8[k]
                      + gp_14[k];

            t_17[k] = ab_z[k] * fp_8[k]
                      + gp_17[k];
        }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, ab_x, ab_y, fp_9, fp_10, fp_11, gp_9, \
                         gp_10, gp_11, gp_19, gp_20 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_18[k] = ab_x[k] * fp_9[k]
                      + gp_9[k];

            t_19[k] = ab_x[k] * fp_10[k]
                      + gp_10[k];

            t_20[k] = ab_x[k] * fp_11[k]
                      + gp_11[k];

            t_21[k] = ab_y[k] * fp_10[k]
                      + gp_19[k];

            t_22[k] = ab_y[k] * fp_11[k]
                      + gp_20[k];
        }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, ab_x, ab_z, fp_11, fp_12, fp_13, fp_14, \
                         gp_12, gp_13, gp_14, gp_23 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_23[k] = ab_z[k] * fp_11[k]
                      + gp_23[k];

            t_24[k] = ab_x[k] * fp_12[k]
                      + gp_12[k];

            t_25[k] = ab_x[k] * fp_13[k]
                      + gp_13[k];

            t_26[k] = ab_x[k] * fp_14[k]
                      + gp_14[k];
        }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, ab_x, ab_y, ab_z, fp_13, fp_14, fp_15, gp_15, \
                         gp_22, gp_23, gp_26 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_27[k] = ab_y[k] * fp_13[k]
                      + gp_22[k];

            t_28[k] = ab_y[k] * fp_14[k]
                      + gp_23[k];

            t_29[k] = ab_z[k] * fp_14[k]
                      + gp_26[k];

            t_30[k] = ab_x[k] * fp_15[k]
                      + gp_15[k];
        }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, ab_x, ab_y, ab_z, fp_16, fp_17, gp_16, \
                         gp_17, gp_25, gp_26, gp_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_31[k] = ab_x[k] * fp_16[k]
                      + gp_16[k];

            t_32[k] = ab_x[k] * fp_17[k]
                      + gp_17[k];

            t_33[k] = ab_y[k] * fp_16[k]
                      + gp_25[k];

            t_34[k] = ab_y[k] * fp_17[k]
                      + gp_26[k];

            t_35[k] = ab_z[k] * fp_17[k]
                      + gp_29[k];
        }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, ab_x, ab_y, fp_18, fp_19, fp_20, gp_18, \
                         gp_19, gp_20, gp_31, gp_32 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_36[k] = ab_x[k] * fp_18[k]
                      + gp_18[k];

            t_37[k] = ab_x[k] * fp_19[k]
                      + gp_19[k];

            t_38[k] = ab_x[k] * fp_20[k]
                      + gp_20[k];

            t_39[k] = ab_y[k] * fp_19[k]
                      + gp_31[k];

            t_40[k] = ab_y[k] * fp_20[k]
                      + gp_32[k];
        }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, ab_x, ab_z, fp_20, fp_21, fp_22, fp_23, \
                         gp_21, gp_22, gp_23, gp_35 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_41[k] = ab_z[k] * fp_20[k]
                      + gp_35[k];

            t_42[k] = ab_x[k] * fp_21[k]
                      + gp_21[k];

            t_43[k] = ab_x[k] * fp_22[k]
                      + gp_22[k];

            t_44[k] = ab_x[k] * fp_23[k]
                      + gp_23[k];
        }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, ab_x, ab_y, ab_z, fp_22, fp_23, fp_24, gp_24, \
                         gp_34, gp_35, gp_38 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_45[k] = ab_y[k] * fp_22[k]
                      + gp_34[k];

            t_46[k] = ab_y[k] * fp_23[k]
                      + gp_35[k];

            t_47[k] = ab_z[k] * fp_23[k]
                      + gp_38[k];

            t_48[k] = ab_x[k] * fp_24[k]
                      + gp_24[k];
        }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, ab_x, ab_y, ab_z, fp_25, fp_26, gp_25, \
                         gp_26, gp_37, gp_38, gp_41 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_49[k] = ab_x[k] * fp_25[k]
                      + gp_25[k];

            t_50[k] = ab_x[k] * fp_26[k]
                      + gp_26[k];

            t_51[k] = ab_y[k] * fp_25[k]
                      + gp_37[k];

            t_52[k] = ab_y[k] * fp_26[k]
                      + gp_38[k];

            t_53[k] = ab_z[k] * fp_26[k]
                      + gp_41[k];
        }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, ab_x, ab_y, fp_27, fp_28, fp_29, gp_27, \
                         gp_28, gp_29, gp_40, gp_41 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_54[k] = ab_x[k] * fp_27[k]
                      + gp_27[k];

            t_55[k] = ab_x[k] * fp_28[k]
                      + gp_28[k];

            t_56[k] = ab_x[k] * fp_29[k]
                      + gp_29[k];

            t_57[k] = ab_y[k] * fp_28[k]
                      + gp_40[k];

            t_58[k] = ab_y[k] * fp_29[k]
                      + gp_41[k];
        }

#pragma omp simd aligned(t_59, ab_z, fp_29, gp_44 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_59[k] = ab_z[k] * fp_29[k]
                      + gp_44[k];
        }
    }
}

auto
compute_hrr_fd_out_of_second(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                             const size_t target, const size_t dd, const size_t df,
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

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *dd_0 = buffer.data(dd + 0 * ncomps + c);
        const auto *dd_1 = buffer.data(dd + 1 * ncomps + c);
        const auto *dd_2 = buffer.data(dd + 2 * ncomps + c);
        const auto *dd_3 = buffer.data(dd + 3 * ncomps + c);
        const auto *dd_4 = buffer.data(dd + 4 * ncomps + c);
        const auto *dd_5 = buffer.data(dd + 5 * ncomps + c);
        const auto *dd_6 = buffer.data(dd + 6 * ncomps + c);
        const auto *dd_7 = buffer.data(dd + 7 * ncomps + c);
        const auto *dd_8 = buffer.data(dd + 8 * ncomps + c);
        const auto *dd_9 = buffer.data(dd + 9 * ncomps + c);
        const auto *dd_10 = buffer.data(dd + 10 * ncomps + c);
        const auto *dd_11 = buffer.data(dd + 11 * ncomps + c);
        const auto *dd_12 = buffer.data(dd + 12 * ncomps + c);
        const auto *dd_13 = buffer.data(dd + 13 * ncomps + c);
        const auto *dd_14 = buffer.data(dd + 14 * ncomps + c);
        const auto *dd_15 = buffer.data(dd + 15 * ncomps + c);
        const auto *dd_16 = buffer.data(dd + 16 * ncomps + c);
        const auto *dd_17 = buffer.data(dd + 17 * ncomps + c);
        const auto *dd_18 = buffer.data(dd + 18 * ncomps + c);
        const auto *dd_19 = buffer.data(dd + 19 * ncomps + c);
        const auto *dd_20 = buffer.data(dd + 20 * ncomps + c);
        const auto *dd_21 = buffer.data(dd + 21 * ncomps + c);
        const auto *dd_22 = buffer.data(dd + 22 * ncomps + c);
        const auto *dd_23 = buffer.data(dd + 23 * ncomps + c);
        const auto *dd_24 = buffer.data(dd + 24 * ncomps + c);
        const auto *dd_25 = buffer.data(dd + 25 * ncomps + c);
        const auto *dd_26 = buffer.data(dd + 26 * ncomps + c);
        const auto *dd_27 = buffer.data(dd + 27 * ncomps + c);
        const auto *dd_28 = buffer.data(dd + 28 * ncomps + c);
        const auto *dd_29 = buffer.data(dd + 29 * ncomps + c);
        const auto *dd_30 = buffer.data(dd + 30 * ncomps + c);
        const auto *dd_31 = buffer.data(dd + 31 * ncomps + c);
        const auto *dd_32 = buffer.data(dd + 32 * ncomps + c);
        const auto *dd_33 = buffer.data(dd + 33 * ncomps + c);
        const auto *dd_34 = buffer.data(dd + 34 * ncomps + c);
        const auto *dd_35 = buffer.data(dd + 35 * ncomps + c);

        const auto *df_0 = buffer.data(df + 0 * ncomps + c);
        const auto *df_1 = buffer.data(df + 1 * ncomps + c);
        const auto *df_2 = buffer.data(df + 2 * ncomps + c);
        const auto *df_3 = buffer.data(df + 3 * ncomps + c);
        const auto *df_4 = buffer.data(df + 4 * ncomps + c);
        const auto *df_5 = buffer.data(df + 5 * ncomps + c);
        const auto *df_10 = buffer.data(df + 10 * ncomps + c);
        const auto *df_11 = buffer.data(df + 11 * ncomps + c);
        const auto *df_12 = buffer.data(df + 12 * ncomps + c);
        const auto *df_13 = buffer.data(df + 13 * ncomps + c);
        const auto *df_14 = buffer.data(df + 14 * ncomps + c);
        const auto *df_15 = buffer.data(df + 15 * ncomps + c);
        const auto *df_20 = buffer.data(df + 20 * ncomps + c);
        const auto *df_21 = buffer.data(df + 21 * ncomps + c);
        const auto *df_22 = buffer.data(df + 22 * ncomps + c);
        const auto *df_23 = buffer.data(df + 23 * ncomps + c);
        const auto *df_24 = buffer.data(df + 24 * ncomps + c);
        const auto *df_25 = buffer.data(df + 25 * ncomps + c);
        const auto *df_30 = buffer.data(df + 30 * ncomps + c);
        const auto *df_31 = buffer.data(df + 31 * ncomps + c);
        const auto *df_32 = buffer.data(df + 32 * ncomps + c);
        const auto *df_33 = buffer.data(df + 33 * ncomps + c);
        const auto *df_34 = buffer.data(df + 34 * ncomps + c);
        const auto *df_35 = buffer.data(df + 35 * ncomps + c);
        const auto *df_36 = buffer.data(df + 36 * ncomps + c);
        const auto *df_37 = buffer.data(df + 37 * ncomps + c);
        const auto *df_38 = buffer.data(df + 38 * ncomps + c);
        const auto *df_40 = buffer.data(df + 40 * ncomps + c);
        const auto *df_41 = buffer.data(df + 41 * ncomps + c);
        const auto *df_42 = buffer.data(df + 42 * ncomps + c);
        const auto *df_43 = buffer.data(df + 43 * ncomps + c);
        const auto *df_44 = buffer.data(df + 44 * ncomps + c);
        const auto *df_45 = buffer.data(df + 45 * ncomps + c);
        const auto *df_46 = buffer.data(df + 46 * ncomps + c);
        const auto *df_47 = buffer.data(df + 47 * ncomps + c);
        const auto *df_48 = buffer.data(df + 48 * ncomps + c);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, dd_0, dd_1, dd_2, dd_3, dd_4, df_0, \
                         df_1, df_2, df_3, df_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = -ab_x[k] * dd_0[k]
                     + df_0[k];

            t_1[k] = -ab_x[k] * dd_1[k]
                     + df_1[k];

            t_2[k] = -ab_x[k] * dd_2[k]
                     + df_2[k];

            t_3[k] = -ab_x[k] * dd_3[k]
                     + df_3[k];

            t_4[k] = -ab_x[k] * dd_4[k]
                     + df_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, dd_5, dd_6, dd_7, dd_8, dd_9, df_5, \
                         df_10, df_11, df_12, df_13 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = -ab_x[k] * dd_5[k]
                     + df_5[k];

            t_6[k] = -ab_x[k] * dd_6[k]
                     + df_10[k];

            t_7[k] = -ab_x[k] * dd_7[k]
                     + df_11[k];

            t_8[k] = -ab_x[k] * dd_8[k]
                     + df_12[k];

            t_9[k] = -ab_x[k] * dd_9[k]
                     + df_13[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, dd_10, dd_11, dd_12, dd_13, \
                         dd_14, df_14, df_15, df_20, df_21, df_22 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = -ab_x[k] * dd_10[k]
                      + df_14[k];

            t_11[k] = -ab_x[k] * dd_11[k]
                      + df_15[k];

            t_12[k] = -ab_x[k] * dd_12[k]
                      + df_20[k];

            t_13[k] = -ab_x[k] * dd_13[k]
                      + df_21[k];

            t_14[k] = -ab_x[k] * dd_14[k]
                      + df_22[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, dd_15, dd_16, dd_17, dd_18, \
                         dd_19, df_23, df_24, df_25, df_30, df_31 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = -ab_x[k] * dd_15[k]
                      + df_23[k];

            t_16[k] = -ab_x[k] * dd_16[k]
                      + df_24[k];

            t_17[k] = -ab_x[k] * dd_17[k]
                      + df_25[k];

            t_18[k] = -ab_x[k] * dd_18[k]
                      + df_30[k];

            t_19[k] = -ab_x[k] * dd_19[k]
                      + df_31[k];
        }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, dd_20, dd_21, dd_22, dd_23, \
                         dd_24, df_32, df_33, df_34, df_35, df_40 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_20[k] = -ab_x[k] * dd_20[k]
                      + df_32[k];

            t_21[k] = -ab_x[k] * dd_21[k]
                      + df_33[k];

            t_22[k] = -ab_x[k] * dd_22[k]
                      + df_34[k];

            t_23[k] = -ab_x[k] * dd_23[k]
                      + df_35[k];

            t_24[k] = -ab_x[k] * dd_24[k]
                      + df_40[k];
        }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, dd_25, dd_26, dd_27, dd_28, \
                         dd_29, df_41, df_42, df_43, df_44, df_45 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_25[k] = -ab_x[k] * dd_25[k]
                      + df_41[k];

            t_26[k] = -ab_x[k] * dd_26[k]
                      + df_42[k];

            t_27[k] = -ab_x[k] * dd_27[k]
                      + df_43[k];

            t_28[k] = -ab_x[k] * dd_28[k]
                      + df_44[k];

            t_29[k] = -ab_x[k] * dd_29[k]
                      + df_45[k];
        }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, dd_30, dd_31, dd_32, dd_33, \
                         dd_34, df_50, df_51, df_52, df_53, df_54 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_30[k] = -ab_x[k] * dd_30[k]
                      + df_50[k];

            t_31[k] = -ab_x[k] * dd_31[k]
                      + df_51[k];

            t_32[k] = -ab_x[k] * dd_32[k]
                      + df_52[k];

            t_33[k] = -ab_x[k] * dd_33[k]
                      + df_53[k];

            t_34[k] = -ab_x[k] * dd_34[k]
                      + df_54[k];
        }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, ab_x, ab_y, dd_18, dd_19, dd_20, dd_35, \
                         df_31, df_33, df_34, df_55 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_35[k] = -ab_x[k] * dd_35[k]
                      + df_55[k];

            t_36[k] = -ab_y[k] * dd_18[k]
                      + df_31[k];

            t_37[k] = -ab_y[k] * dd_19[k]
                      + df_33[k];

            t_38[k] = -ab_y[k] * dd_20[k]
                      + df_34[k];
        }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, t_43, ab_y, dd_21, dd_22, dd_23, dd_24, \
                         dd_25, df_36, df_37, df_38, df_41, df_43 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_39[k] = -ab_y[k] * dd_21[k]
                      + df_36[k];

            t_40[k] = -ab_y[k] * dd_22[k]
                      + df_37[k];

            t_41[k] = -ab_y[k] * dd_23[k]
                      + df_38[k];

            t_42[k] = -ab_y[k] * dd_24[k]
                      + df_41[k];

            t_43[k] = -ab_y[k] * dd_25[k]
                      + df_43[k];
        }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, ab_y, dd_26, dd_27, dd_28, dd_29, \
                         dd_30, df_44, df_46, df_47, df_48, df_51 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_44[k] = -ab_y[k] * dd_26[k]
                      + df_44[k];

            t_45[k] = -ab_y[k] * dd_27[k]
                      + df_46[k];

            t_46[k] = -ab_y[k] * dd_28[k]
                      + df_47[k];

            t_47[k] = -ab_y[k] * dd_29[k]
                      + df_48[k];

            t_48[k] = -ab_y[k] * dd_30[k]
                      + df_51[k];
        }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, ab_y, dd_31, dd_32, dd_33, dd_34, \
                         dd_35, df_53, df_54, df_56, df_57, df_58 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_49[k] = -ab_y[k] * dd_31[k]
                      + df_53[k];

            t_50[k] = -ab_y[k] * dd_32[k]
                      + df_54[k];

            t_51[k] = -ab_y[k] * dd_33[k]
                      + df_56[k];

            t_52[k] = -ab_y[k] * dd_34[k]
                      + df_57[k];

            t_53[k] = -ab_y[k] * dd_35[k]
                      + df_58[k];
        }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, ab_z, dd_30, dd_31, dd_32, dd_33, \
                         dd_34, df_52, df_54, df_55, df_57, df_58 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_54[k] = -ab_z[k] * dd_30[k]
                      + df_52[k];

            t_55[k] = -ab_z[k] * dd_31[k]
                      + df_54[k];

            t_56[k] = -ab_z[k] * dd_32[k]
                      + df_55[k];

            t_57[k] = -ab_z[k] * dd_33[k]
                      + df_57[k];

            t_58[k] = -ab_z[k] * dd_34[k]
                      + df_58[k];
        }

#pragma omp simd aligned(t_59, ab_z, dd_35, df_59 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_59[k] = -ab_z[k] * dd_35[k]
                      + df_59[k];
        }
    }
}

auto
compute_hrr_fd(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t fp, const size_t gp, const size_t ncomps, const size_t nmax) -> void
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
        const auto *fp_18 = buffer.data(fp + 18 * ncomps + c);
        const auto *fp_19 = buffer.data(fp + 19 * ncomps + c);
        const auto *fp_20 = buffer.data(fp + 20 * ncomps + c);
        const auto *fp_21 = buffer.data(fp + 21 * ncomps + c);
        const auto *fp_22 = buffer.data(fp + 22 * ncomps + c);
        const auto *fp_23 = buffer.data(fp + 23 * ncomps + c);
        const auto *fp_24 = buffer.data(fp + 24 * ncomps + c);
        const auto *fp_25 = buffer.data(fp + 25 * ncomps + c);
        const auto *fp_26 = buffer.data(fp + 26 * ncomps + c);
        const auto *fp_27 = buffer.data(fp + 27 * ncomps + c);
        const auto *fp_28 = buffer.data(fp + 28 * ncomps + c);
        const auto *fp_29 = buffer.data(fp + 29 * ncomps + c);

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
        const auto *gp_31 = buffer.data(gp + 31 * ncomps + c);
        const auto *gp_32 = buffer.data(gp + 32 * ncomps + c);
        const auto *gp_34 = buffer.data(gp + 34 * ncomps + c);
        const auto *gp_35 = buffer.data(gp + 35 * ncomps + c);
        const auto *gp_37 = buffer.data(gp + 37 * ncomps + c);
        const auto *gp_38 = buffer.data(gp + 38 * ncomps + c);
        const auto *gp_40 = buffer.data(gp + 40 * ncomps + c);
        const auto *gp_41 = buffer.data(gp + 41 * ncomps + c);
        const auto *gp_44 = buffer.data(gp + 44 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, ab_y, fp_0, fp_1, fp_2, gp_0, gp_1, \
                         gp_2, gp_4, gp_5 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * fp_0[k]
                     + gp_0[k];

            t_1[k] = ab_x[k] * fp_1[k]
                     + gp_1[k];

            t_2[k] = ab_x[k] * fp_2[k]
                     + gp_2[k];

            t_3[k] = ab_y[k] * fp_1[k]
                     + gp_4[k];

            t_4[k] = ab_y[k] * fp_2[k]
                     + gp_5[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, ab_x, ab_z, fp_2, fp_3, fp_4, fp_5, gp_3, gp_4, \
                         gp_5, gp_8 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = ab_z[k] * fp_2[k]
                     + gp_8[k];

            t_6[k] = ab_x[k] * fp_3[k]
                     + gp_3[k];

            t_7[k] = ab_x[k] * fp_4[k]
                     + gp_4[k];

            t_8[k] = ab_x[k] * fp_5[k]
                     + gp_5[k];
        }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, ab_x, ab_y, ab_z, fp_4, fp_5, fp_6, gp_6, \
                         gp_10, gp_11, gp_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_9[k] = ab_y[k] * fp_4[k]
                     + gp_10[k];

            t_10[k] = ab_y[k] * fp_5[k]
                      + gp_11[k];

            t_11[k] = ab_z[k] * fp_5[k]
                      + gp_14[k];

            t_12[k] = ab_x[k] * fp_6[k]
                      + gp_6[k];
        }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, ab_x, ab_y, ab_z, fp_7, fp_8, gp_7, \
                         gp_8, gp_13, gp_14, gp_17 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_13[k] = ab_x[k] * fp_7[k]
                      + gp_7[k];

            t_14[k] = ab_x[k] * fp_8[k]
                      + gp_8[k];

            t_15[k] = ab_y[k] * fp_7[k]
                      + gp_13[k];

            t_16[k] = ab_y[k] * fp_8[k]
                      + gp_14[k];

            t_17[k] = ab_z[k] * fp_8[k]
                      + gp_17[k];
        }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, ab_x, ab_y, fp_9, fp_10, fp_11, gp_9, \
                         gp_10, gp_11, gp_19, gp_20 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_18[k] = ab_x[k] * fp_9[k]
                      + gp_9[k];

            t_19[k] = ab_x[k] * fp_10[k]
                      + gp_10[k];

            t_20[k] = ab_x[k] * fp_11[k]
                      + gp_11[k];

            t_21[k] = ab_y[k] * fp_10[k]
                      + gp_19[k];

            t_22[k] = ab_y[k] * fp_11[k]
                      + gp_20[k];
        }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, ab_x, ab_z, fp_11, fp_12, fp_13, fp_14, \
                         gp_12, gp_13, gp_14, gp_23 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_23[k] = ab_z[k] * fp_11[k]
                      + gp_23[k];

            t_24[k] = ab_x[k] * fp_12[k]
                      + gp_12[k];

            t_25[k] = ab_x[k] * fp_13[k]
                      + gp_13[k];

            t_26[k] = ab_x[k] * fp_14[k]
                      + gp_14[k];
        }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, ab_x, ab_y, ab_z, fp_13, fp_14, fp_15, gp_15, \
                         gp_22, gp_23, gp_26 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_27[k] = ab_y[k] * fp_13[k]
                      + gp_22[k];

            t_28[k] = ab_y[k] * fp_14[k]
                      + gp_23[k];

            t_29[k] = ab_z[k] * fp_14[k]
                      + gp_26[k];

            t_30[k] = ab_x[k] * fp_15[k]
                      + gp_15[k];
        }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, ab_x, ab_y, ab_z, fp_16, fp_17, gp_16, \
                         gp_17, gp_25, gp_26, gp_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_31[k] = ab_x[k] * fp_16[k]
                      + gp_16[k];

            t_32[k] = ab_x[k] * fp_17[k]
                      + gp_17[k];

            t_33[k] = ab_y[k] * fp_16[k]
                      + gp_25[k];

            t_34[k] = ab_y[k] * fp_17[k]
                      + gp_26[k];

            t_35[k] = ab_z[k] * fp_17[k]
                      + gp_29[k];
        }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, ab_x, ab_y, fp_18, fp_19, fp_20, gp_18, \
                         gp_19, gp_20, gp_31, gp_32 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_36[k] = ab_x[k] * fp_18[k]
                      + gp_18[k];

            t_37[k] = ab_x[k] * fp_19[k]
                      + gp_19[k];

            t_38[k] = ab_x[k] * fp_20[k]
                      + gp_20[k];

            t_39[k] = ab_y[k] * fp_19[k]
                      + gp_31[k];

            t_40[k] = ab_y[k] * fp_20[k]
                      + gp_32[k];
        }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, ab_x, ab_z, fp_20, fp_21, fp_22, fp_23, \
                         gp_21, gp_22, gp_23, gp_35 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_41[k] = ab_z[k] * fp_20[k]
                      + gp_35[k];

            t_42[k] = ab_x[k] * fp_21[k]
                      + gp_21[k];

            t_43[k] = ab_x[k] * fp_22[k]
                      + gp_22[k];

            t_44[k] = ab_x[k] * fp_23[k]
                      + gp_23[k];
        }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, ab_x, ab_y, ab_z, fp_22, fp_23, fp_24, gp_24, \
                         gp_34, gp_35, gp_38 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_45[k] = ab_y[k] * fp_22[k]
                      + gp_34[k];

            t_46[k] = ab_y[k] * fp_23[k]
                      + gp_35[k];

            t_47[k] = ab_z[k] * fp_23[k]
                      + gp_38[k];

            t_48[k] = ab_x[k] * fp_24[k]
                      + gp_24[k];
        }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, ab_x, ab_y, ab_z, fp_25, fp_26, gp_25, \
                         gp_26, gp_37, gp_38, gp_41 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_49[k] = ab_x[k] * fp_25[k]
                      + gp_25[k];

            t_50[k] = ab_x[k] * fp_26[k]
                      + gp_26[k];

            t_51[k] = ab_y[k] * fp_25[k]
                      + gp_37[k];

            t_52[k] = ab_y[k] * fp_26[k]
                      + gp_38[k];

            t_53[k] = ab_z[k] * fp_26[k]
                      + gp_41[k];
        }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, ab_x, ab_y, fp_27, fp_28, fp_29, gp_27, \
                         gp_28, gp_29, gp_40, gp_41 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_54[k] = ab_x[k] * fp_27[k]
                      + gp_27[k];

            t_55[k] = ab_x[k] * fp_28[k]
                      + gp_28[k];

            t_56[k] = ab_x[k] * fp_29[k]
                      + gp_29[k];

            t_57[k] = ab_y[k] * fp_28[k]
                      + gp_40[k];

            t_58[k] = ab_y[k] * fp_29[k]
                      + gp_41[k];
        }

#pragma omp simd aligned(t_59, ab_z, fp_29, gp_44 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_59[k] = ab_z[k] * fp_29[k]
                      + gp_44[k];
        }
    }
}

}  // namespace simdtrf
