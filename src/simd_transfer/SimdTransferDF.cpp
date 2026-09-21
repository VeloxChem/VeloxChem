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


#include "SimdTransferDF.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_df(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t pf, const size_t pg, const size_t ncomps, const size_t nmax) -> void
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

        const auto *pf_0 = buffer.data(pf + 0 * ncomps + c);
        const auto *pf_1 = buffer.data(pf + 1 * ncomps + c);
        const auto *pf_2 = buffer.data(pf + 2 * ncomps + c);
        const auto *pf_3 = buffer.data(pf + 3 * ncomps + c);
        const auto *pf_4 = buffer.data(pf + 4 * ncomps + c);
        const auto *pf_5 = buffer.data(pf + 5 * ncomps + c);
        const auto *pf_6 = buffer.data(pf + 6 * ncomps + c);
        const auto *pf_7 = buffer.data(pf + 7 * ncomps + c);
        const auto *pf_8 = buffer.data(pf + 8 * ncomps + c);
        const auto *pf_9 = buffer.data(pf + 9 * ncomps + c);
        const auto *pf_10 = buffer.data(pf + 10 * ncomps + c);
        const auto *pf_11 = buffer.data(pf + 11 * ncomps + c);
        const auto *pf_12 = buffer.data(pf + 12 * ncomps + c);
        const auto *pf_13 = buffer.data(pf + 13 * ncomps + c);
        const auto *pf_14 = buffer.data(pf + 14 * ncomps + c);
        const auto *pf_15 = buffer.data(pf + 15 * ncomps + c);
        const auto *pf_16 = buffer.data(pf + 16 * ncomps + c);
        const auto *pf_17 = buffer.data(pf + 17 * ncomps + c);
        const auto *pf_18 = buffer.data(pf + 18 * ncomps + c);
        const auto *pf_19 = buffer.data(pf + 19 * ncomps + c);
        const auto *pf_20 = buffer.data(pf + 20 * ncomps + c);
        const auto *pf_21 = buffer.data(pf + 21 * ncomps + c);
        const auto *pf_22 = buffer.data(pf + 22 * ncomps + c);
        const auto *pf_23 = buffer.data(pf + 23 * ncomps + c);
        const auto *pf_24 = buffer.data(pf + 24 * ncomps + c);
        const auto *pf_25 = buffer.data(pf + 25 * ncomps + c);
        const auto *pf_26 = buffer.data(pf + 26 * ncomps + c);
        const auto *pf_27 = buffer.data(pf + 27 * ncomps + c);
        const auto *pf_28 = buffer.data(pf + 28 * ncomps + c);
        const auto *pf_29 = buffer.data(pf + 29 * ncomps + c);

        const auto *pg_0 = buffer.data(pg + 0 * ncomps + c);
        const auto *pg_1 = buffer.data(pg + 1 * ncomps + c);
        const auto *pg_2 = buffer.data(pg + 2 * ncomps + c);
        const auto *pg_3 = buffer.data(pg + 3 * ncomps + c);
        const auto *pg_4 = buffer.data(pg + 4 * ncomps + c);
        const auto *pg_5 = buffer.data(pg + 5 * ncomps + c);
        const auto *pg_6 = buffer.data(pg + 6 * ncomps + c);
        const auto *pg_7 = buffer.data(pg + 7 * ncomps + c);
        const auto *pg_8 = buffer.data(pg + 8 * ncomps + c);
        const auto *pg_9 = buffer.data(pg + 9 * ncomps + c);
        const auto *pg_15 = buffer.data(pg + 15 * ncomps + c);
        const auto *pg_16 = buffer.data(pg + 16 * ncomps + c);
        const auto *pg_17 = buffer.data(pg + 17 * ncomps + c);
        const auto *pg_18 = buffer.data(pg + 18 * ncomps + c);
        const auto *pg_19 = buffer.data(pg + 19 * ncomps + c);
        const auto *pg_20 = buffer.data(pg + 20 * ncomps + c);
        const auto *pg_21 = buffer.data(pg + 21 * ncomps + c);
        const auto *pg_22 = buffer.data(pg + 22 * ncomps + c);
        const auto *pg_23 = buffer.data(pg + 23 * ncomps + c);
        const auto *pg_24 = buffer.data(pg + 24 * ncomps + c);
        const auto *pg_25 = buffer.data(pg + 25 * ncomps + c);
        const auto *pg_26 = buffer.data(pg + 26 * ncomps + c);
        const auto *pg_27 = buffer.data(pg + 27 * ncomps + c);
        const auto *pg_28 = buffer.data(pg + 28 * ncomps + c);
        const auto *pg_30 = buffer.data(pg + 30 * ncomps + c);
        const auto *pg_31 = buffer.data(pg + 31 * ncomps + c);
        const auto *pg_32 = buffer.data(pg + 32 * ncomps + c);
        const auto *pg_33 = buffer.data(pg + 33 * ncomps + c);
        const auto *pg_34 = buffer.data(pg + 34 * ncomps + c);
        const auto *pg_35 = buffer.data(pg + 35 * ncomps + c);
        const auto *pg_36 = buffer.data(pg + 36 * ncomps + c);
        const auto *pg_37 = buffer.data(pg + 37 * ncomps + c);
        const auto *pg_38 = buffer.data(pg + 38 * ncomps + c);
        const auto *pg_39 = buffer.data(pg + 39 * ncomps + c);
        const auto *pg_40 = buffer.data(pg + 40 * ncomps + c);
        const auto *pg_41 = buffer.data(pg + 41 * ncomps + c);
        const auto *pg_42 = buffer.data(pg + 42 * ncomps + c);
        const auto *pg_43 = buffer.data(pg + 43 * ncomps + c);
        const auto *pg_44 = buffer.data(pg + 44 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, pf_0, pf_1, pf_2, pf_3, pf_4, pg_0, \
                         pg_1, pg_2, pg_3, pg_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = -ab_x[k] * pf_0[k]
                     + pg_0[k];

            t_1[k] = -ab_x[k] * pf_1[k]
                     + pg_1[k];

            t_2[k] = -ab_x[k] * pf_2[k]
                     + pg_2[k];

            t_3[k] = -ab_x[k] * pf_3[k]
                     + pg_3[k];

            t_4[k] = -ab_x[k] * pf_4[k]
                     + pg_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, pf_5, pf_6, pf_7, pf_8, pf_9, pg_5, \
                         pg_6, pg_7, pg_8, pg_9 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = -ab_x[k] * pf_5[k]
                     + pg_5[k];

            t_6[k] = -ab_x[k] * pf_6[k]
                     + pg_6[k];

            t_7[k] = -ab_x[k] * pf_7[k]
                     + pg_7[k];

            t_8[k] = -ab_x[k] * pf_8[k]
                     + pg_8[k];

            t_9[k] = -ab_x[k] * pf_9[k]
                     + pg_9[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, pf_10, pf_11, pf_12, pf_13, \
                         pf_14, pg_15, pg_16, pg_17, pg_18, pg_19 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = -ab_x[k] * pf_10[k]
                      + pg_15[k];

            t_11[k] = -ab_x[k] * pf_11[k]
                      + pg_16[k];

            t_12[k] = -ab_x[k] * pf_12[k]
                      + pg_17[k];

            t_13[k] = -ab_x[k] * pf_13[k]
                      + pg_18[k];

            t_14[k] = -ab_x[k] * pf_14[k]
                      + pg_19[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, pf_15, pf_16, pf_17, pf_18, \
                         pf_19, pg_20, pg_21, pg_22, pg_23, pg_24 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = -ab_x[k] * pf_15[k]
                      + pg_20[k];

            t_16[k] = -ab_x[k] * pf_16[k]
                      + pg_21[k];

            t_17[k] = -ab_x[k] * pf_17[k]
                      + pg_22[k];

            t_18[k] = -ab_x[k] * pf_18[k]
                      + pg_23[k];

            t_19[k] = -ab_x[k] * pf_19[k]
                      + pg_24[k];
        }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, pf_20, pf_21, pf_22, pf_23, \
                         pf_24, pg_30, pg_31, pg_32, pg_33, pg_34 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_20[k] = -ab_x[k] * pf_20[k]
                      + pg_30[k];

            t_21[k] = -ab_x[k] * pf_21[k]
                      + pg_31[k];

            t_22[k] = -ab_x[k] * pf_22[k]
                      + pg_32[k];

            t_23[k] = -ab_x[k] * pf_23[k]
                      + pg_33[k];

            t_24[k] = -ab_x[k] * pf_24[k]
                      + pg_34[k];
        }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, pf_25, pf_26, pf_27, pf_28, \
                         pf_29, pg_35, pg_36, pg_37, pg_38, pg_39 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_25[k] = -ab_x[k] * pf_25[k]
                      + pg_35[k];

            t_26[k] = -ab_x[k] * pf_26[k]
                      + pg_36[k];

            t_27[k] = -ab_x[k] * pf_27[k]
                      + pg_37[k];

            t_28[k] = -ab_x[k] * pf_28[k]
                      + pg_38[k];

            t_29[k] = -ab_x[k] * pf_29[k]
                      + pg_39[k];
        }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_y, pf_10, pf_11, pf_12, pf_13, \
                         pf_14, pg_16, pg_18, pg_19, pg_21, pg_22 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_30[k] = -ab_y[k] * pf_10[k]
                      + pg_16[k];

            t_31[k] = -ab_y[k] * pf_11[k]
                      + pg_18[k];

            t_32[k] = -ab_y[k] * pf_12[k]
                      + pg_19[k];

            t_33[k] = -ab_y[k] * pf_13[k]
                      + pg_21[k];

            t_34[k] = -ab_y[k] * pf_14[k]
                      + pg_22[k];
        }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_y, pf_15, pf_16, pf_17, pf_18, \
                         pf_19, pg_23, pg_25, pg_26, pg_27, pg_28 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_35[k] = -ab_y[k] * pf_15[k]
                      + pg_23[k];

            t_36[k] = -ab_y[k] * pf_16[k]
                      + pg_25[k];

            t_37[k] = -ab_y[k] * pf_17[k]
                      + pg_26[k];

            t_38[k] = -ab_y[k] * pf_18[k]
                      + pg_27[k];

            t_39[k] = -ab_y[k] * pf_19[k]
                      + pg_28[k];
        }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_y, pf_20, pf_21, pf_22, pf_23, \
                         pf_24, pg_31, pg_33, pg_34, pg_36, pg_37 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_40[k] = -ab_y[k] * pf_20[k]
                      + pg_31[k];

            t_41[k] = -ab_y[k] * pf_21[k]
                      + pg_33[k];

            t_42[k] = -ab_y[k] * pf_22[k]
                      + pg_34[k];

            t_43[k] = -ab_y[k] * pf_23[k]
                      + pg_36[k];

            t_44[k] = -ab_y[k] * pf_24[k]
                      + pg_37[k];
        }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_y, pf_25, pf_26, pf_27, pf_28, \
                         pf_29, pg_38, pg_40, pg_41, pg_42, pg_43 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_45[k] = -ab_y[k] * pf_25[k]
                      + pg_38[k];

            t_46[k] = -ab_y[k] * pf_26[k]
                      + pg_40[k];

            t_47[k] = -ab_y[k] * pf_27[k]
                      + pg_41[k];

            t_48[k] = -ab_y[k] * pf_28[k]
                      + pg_42[k];

            t_49[k] = -ab_y[k] * pf_29[k]
                      + pg_43[k];
        }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_z, pf_20, pf_21, pf_22, pf_23, \
                         pf_24, pg_32, pg_34, pg_35, pg_37, pg_38 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_50[k] = -ab_z[k] * pf_20[k]
                      + pg_32[k];

            t_51[k] = -ab_z[k] * pf_21[k]
                      + pg_34[k];

            t_52[k] = -ab_z[k] * pf_22[k]
                      + pg_35[k];

            t_53[k] = -ab_z[k] * pf_23[k]
                      + pg_37[k];

            t_54[k] = -ab_z[k] * pf_24[k]
                      + pg_38[k];
        }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_z, pf_25, pf_26, pf_27, pf_28, \
                         pf_29, pg_39, pg_41, pg_42, pg_43, pg_44 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_55[k] = -ab_z[k] * pf_25[k]
                      + pg_39[k];

            t_56[k] = -ab_z[k] * pf_26[k]
                      + pg_41[k];

            t_57[k] = -ab_z[k] * pf_27[k]
                      + pg_42[k];

            t_58[k] = -ab_z[k] * pf_28[k]
                      + pg_43[k];

            t_59[k] = -ab_z[k] * pf_29[k]
                      + pg_44[k];
        }
    }
}

auto
compute_hrr_df_out_of_first(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                            const size_t target, const size_t dd, const size_t fd,
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
        const auto *fd_39 = buffer.data(fd + 39 * ncomps + c);
        const auto *fd_40 = buffer.data(fd + 40 * ncomps + c);
        const auto *fd_41 = buffer.data(fd + 41 * ncomps + c);
        const auto *fd_45 = buffer.data(fd + 45 * ncomps + c);
        const auto *fd_46 = buffer.data(fd + 46 * ncomps + c);
        const auto *fd_47 = buffer.data(fd + 47 * ncomps + c);
        const auto *fd_51 = buffer.data(fd + 51 * ncomps + c);
        const auto *fd_52 = buffer.data(fd + 52 * ncomps + c);
        const auto *fd_53 = buffer.data(fd + 53 * ncomps + c);
        const auto *fd_59 = buffer.data(fd + 59 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, dd_0, dd_1, dd_2, dd_3, dd_4, fd_0, \
                         fd_1, fd_2, fd_3, fd_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * dd_0[k]
                     + fd_0[k];

            t_1[k] = ab_x[k] * dd_1[k]
                     + fd_1[k];

            t_2[k] = ab_x[k] * dd_2[k]
                     + fd_2[k];

            t_3[k] = ab_x[k] * dd_3[k]
                     + fd_3[k];

            t_4[k] = ab_x[k] * dd_4[k]
                     + fd_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, ab_y, ab_z, dd_3, dd_4, dd_5, fd_5, \
                         fd_9, fd_10, fd_11, fd_17 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = ab_x[k] * dd_5[k]
                     + fd_5[k];

            t_6[k] = ab_y[k] * dd_3[k]
                     + fd_9[k];

            t_7[k] = ab_y[k] * dd_4[k]
                     + fd_10[k];

            t_8[k] = ab_y[k] * dd_5[k]
                     + fd_11[k];

            t_9[k] = ab_z[k] * dd_5[k]
                     + fd_17[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, dd_6, dd_7, dd_8, dd_9, dd_10, \
                         fd_6, fd_7, fd_8, fd_9, fd_10 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = ab_x[k] * dd_6[k]
                      + fd_6[k];

            t_11[k] = ab_x[k] * dd_7[k]
                      + fd_7[k];

            t_12[k] = ab_x[k] * dd_8[k]
                      + fd_8[k];

            t_13[k] = ab_x[k] * dd_9[k]
                      + fd_9[k];

            t_14[k] = ab_x[k] * dd_10[k]
                      + fd_10[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, ab_y, ab_z, dd_9, dd_10, dd_11, \
                         fd_11, fd_21, fd_22, fd_23, fd_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = ab_x[k] * dd_11[k]
                      + fd_11[k];

            t_16[k] = ab_y[k] * dd_9[k]
                      + fd_21[k];

            t_17[k] = ab_y[k] * dd_10[k]
                      + fd_22[k];

            t_18[k] = ab_y[k] * dd_11[k]
                      + fd_23[k];

            t_19[k] = ab_z[k] * dd_11[k]
                      + fd_29[k];
        }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, dd_12, dd_13, dd_14, dd_15, \
                         dd_16, fd_12, fd_13, fd_14, fd_15, fd_16 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_20[k] = ab_x[k] * dd_12[k]
                      + fd_12[k];

            t_21[k] = ab_x[k] * dd_13[k]
                      + fd_13[k];

            t_22[k] = ab_x[k] * dd_14[k]
                      + fd_14[k];

            t_23[k] = ab_x[k] * dd_15[k]
                      + fd_15[k];

            t_24[k] = ab_x[k] * dd_16[k]
                      + fd_16[k];
        }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, ab_y, ab_z, dd_15, dd_16, dd_17, \
                         fd_17, fd_27, fd_28, fd_29, fd_35 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_25[k] = ab_x[k] * dd_17[k]
                      + fd_17[k];

            t_26[k] = ab_y[k] * dd_15[k]
                      + fd_27[k];

            t_27[k] = ab_y[k] * dd_16[k]
                      + fd_28[k];

            t_28[k] = ab_y[k] * dd_17[k]
                      + fd_29[k];

            t_29[k] = ab_z[k] * dd_17[k]
                      + fd_35[k];
        }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, dd_18, dd_19, dd_20, dd_21, \
                         dd_22, fd_18, fd_19, fd_20, fd_21, fd_22 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_30[k] = ab_x[k] * dd_18[k]
                      + fd_18[k];

            t_31[k] = ab_x[k] * dd_19[k]
                      + fd_19[k];

            t_32[k] = ab_x[k] * dd_20[k]
                      + fd_20[k];

            t_33[k] = ab_x[k] * dd_21[k]
                      + fd_21[k];

            t_34[k] = ab_x[k] * dd_22[k]
                      + fd_22[k];
        }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, ab_y, ab_z, dd_21, dd_22, dd_23, \
                         fd_23, fd_39, fd_40, fd_41, fd_47 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_35[k] = ab_x[k] * dd_23[k]
                      + fd_23[k];

            t_36[k] = ab_y[k] * dd_21[k]
                      + fd_39[k];

            t_37[k] = ab_y[k] * dd_22[k]
                      + fd_40[k];

            t_38[k] = ab_y[k] * dd_23[k]
                      + fd_41[k];

            t_39[k] = ab_z[k] * dd_23[k]
                      + fd_47[k];
        }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_x, dd_24, dd_25, dd_26, dd_27, \
                         dd_28, fd_24, fd_25, fd_26, fd_27, fd_28 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_40[k] = ab_x[k] * dd_24[k]
                      + fd_24[k];

            t_41[k] = ab_x[k] * dd_25[k]
                      + fd_25[k];

            t_42[k] = ab_x[k] * dd_26[k]
                      + fd_26[k];

            t_43[k] = ab_x[k] * dd_27[k]
                      + fd_27[k];

            t_44[k] = ab_x[k] * dd_28[k]
                      + fd_28[k];
        }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, ab_y, ab_z, dd_27, dd_28, dd_29, \
                         fd_29, fd_45, fd_46, fd_47, fd_53 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_45[k] = ab_x[k] * dd_29[k]
                      + fd_29[k];

            t_46[k] = ab_y[k] * dd_27[k]
                      + fd_45[k];

            t_47[k] = ab_y[k] * dd_28[k]
                      + fd_46[k];

            t_48[k] = ab_y[k] * dd_29[k]
                      + fd_47[k];

            t_49[k] = ab_z[k] * dd_29[k]
                      + fd_53[k];
        }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, dd_30, dd_31, dd_32, dd_33, \
                         dd_34, fd_30, fd_31, fd_32, fd_33, fd_34 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_50[k] = ab_x[k] * dd_30[k]
                      + fd_30[k];

            t_51[k] = ab_x[k] * dd_31[k]
                      + fd_31[k];

            t_52[k] = ab_x[k] * dd_32[k]
                      + fd_32[k];

            t_53[k] = ab_x[k] * dd_33[k]
                      + fd_33[k];

            t_54[k] = ab_x[k] * dd_34[k]
                      + fd_34[k];
        }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_x, ab_y, ab_z, dd_33, dd_34, dd_35, \
                         fd_35, fd_51, fd_52, fd_53, fd_59 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_55[k] = ab_x[k] * dd_35[k]
                      + fd_35[k];

            t_56[k] = ab_y[k] * dd_33[k]
                      + fd_51[k];

            t_57[k] = ab_y[k] * dd_34[k]
                      + fd_52[k];

            t_58[k] = ab_y[k] * dd_35[k]
                      + fd_53[k];

            t_59[k] = ab_z[k] * dd_35[k]
                      + fd_59[k];
        }
    }
}

}  // namespace simdtrf
