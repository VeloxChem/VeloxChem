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


#include "SimdTransferGeom010ZDF.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_geom_010z_df(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                         const size_t target, const size_t pf_1, const size_t pf_0,
                         const size_t pg_1, const size_t ncomps, const size_t nmax) -> void
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

        const auto *pf_1_0 = buffer.data(pf_1 + 0 * ncomps + c);
        const auto *pf_1_1 = buffer.data(pf_1 + 1 * ncomps + c);
        const auto *pf_1_2 = buffer.data(pf_1 + 2 * ncomps + c);
        const auto *pf_1_3 = buffer.data(pf_1 + 3 * ncomps + c);
        const auto *pf_1_4 = buffer.data(pf_1 + 4 * ncomps + c);
        const auto *pf_1_5 = buffer.data(pf_1 + 5 * ncomps + c);
        const auto *pf_1_6 = buffer.data(pf_1 + 6 * ncomps + c);
        const auto *pf_1_7 = buffer.data(pf_1 + 7 * ncomps + c);
        const auto *pf_1_8 = buffer.data(pf_1 + 8 * ncomps + c);
        const auto *pf_1_9 = buffer.data(pf_1 + 9 * ncomps + c);
        const auto *pf_1_10 = buffer.data(pf_1 + 10 * ncomps + c);
        const auto *pf_1_11 = buffer.data(pf_1 + 11 * ncomps + c);
        const auto *pf_1_12 = buffer.data(pf_1 + 12 * ncomps + c);
        const auto *pf_1_13 = buffer.data(pf_1 + 13 * ncomps + c);
        const auto *pf_1_14 = buffer.data(pf_1 + 14 * ncomps + c);
        const auto *pf_1_15 = buffer.data(pf_1 + 15 * ncomps + c);
        const auto *pf_1_16 = buffer.data(pf_1 + 16 * ncomps + c);
        const auto *pf_1_17 = buffer.data(pf_1 + 17 * ncomps + c);
        const auto *pf_1_18 = buffer.data(pf_1 + 18 * ncomps + c);
        const auto *pf_1_19 = buffer.data(pf_1 + 19 * ncomps + c);
        const auto *pf_1_20 = buffer.data(pf_1 + 20 * ncomps + c);
        const auto *pf_1_21 = buffer.data(pf_1 + 21 * ncomps + c);
        const auto *pf_1_22 = buffer.data(pf_1 + 22 * ncomps + c);
        const auto *pf_1_23 = buffer.data(pf_1 + 23 * ncomps + c);
        const auto *pf_1_24 = buffer.data(pf_1 + 24 * ncomps + c);
        const auto *pf_1_25 = buffer.data(pf_1 + 25 * ncomps + c);
        const auto *pf_1_26 = buffer.data(pf_1 + 26 * ncomps + c);
        const auto *pf_1_27 = buffer.data(pf_1 + 27 * ncomps + c);
        const auto *pf_1_28 = buffer.data(pf_1 + 28 * ncomps + c);
        const auto *pf_1_29 = buffer.data(pf_1 + 29 * ncomps + c);

        const auto *pf_0_20 = buffer.data(pf_0 + 20 * ncomps + c);
        const auto *pf_0_21 = buffer.data(pf_0 + 21 * ncomps + c);
        const auto *pf_0_22 = buffer.data(pf_0 + 22 * ncomps + c);
        const auto *pf_0_23 = buffer.data(pf_0 + 23 * ncomps + c);
        const auto *pf_0_24 = buffer.data(pf_0 + 24 * ncomps + c);
        const auto *pf_0_25 = buffer.data(pf_0 + 25 * ncomps + c);
        const auto *pf_0_26 = buffer.data(pf_0 + 26 * ncomps + c);
        const auto *pf_0_27 = buffer.data(pf_0 + 27 * ncomps + c);
        const auto *pf_0_28 = buffer.data(pf_0 + 28 * ncomps + c);
        const auto *pf_0_29 = buffer.data(pf_0 + 29 * ncomps + c);

        const auto *pg_1_0 = buffer.data(pg_1 + 0 * ncomps + c);
        const auto *pg_1_1 = buffer.data(pg_1 + 1 * ncomps + c);
        const auto *pg_1_2 = buffer.data(pg_1 + 2 * ncomps + c);
        const auto *pg_1_3 = buffer.data(pg_1 + 3 * ncomps + c);
        const auto *pg_1_4 = buffer.data(pg_1 + 4 * ncomps + c);
        const auto *pg_1_5 = buffer.data(pg_1 + 5 * ncomps + c);
        const auto *pg_1_6 = buffer.data(pg_1 + 6 * ncomps + c);
        const auto *pg_1_7 = buffer.data(pg_1 + 7 * ncomps + c);
        const auto *pg_1_8 = buffer.data(pg_1 + 8 * ncomps + c);
        const auto *pg_1_9 = buffer.data(pg_1 + 9 * ncomps + c);
        const auto *pg_1_15 = buffer.data(pg_1 + 15 * ncomps + c);
        const auto *pg_1_16 = buffer.data(pg_1 + 16 * ncomps + c);
        const auto *pg_1_17 = buffer.data(pg_1 + 17 * ncomps + c);
        const auto *pg_1_18 = buffer.data(pg_1 + 18 * ncomps + c);
        const auto *pg_1_19 = buffer.data(pg_1 + 19 * ncomps + c);
        const auto *pg_1_20 = buffer.data(pg_1 + 20 * ncomps + c);
        const auto *pg_1_21 = buffer.data(pg_1 + 21 * ncomps + c);
        const auto *pg_1_22 = buffer.data(pg_1 + 22 * ncomps + c);
        const auto *pg_1_23 = buffer.data(pg_1 + 23 * ncomps + c);
        const auto *pg_1_24 = buffer.data(pg_1 + 24 * ncomps + c);
        const auto *pg_1_25 = buffer.data(pg_1 + 25 * ncomps + c);
        const auto *pg_1_26 = buffer.data(pg_1 + 26 * ncomps + c);
        const auto *pg_1_27 = buffer.data(pg_1 + 27 * ncomps + c);
        const auto *pg_1_28 = buffer.data(pg_1 + 28 * ncomps + c);
        const auto *pg_1_30 = buffer.data(pg_1 + 30 * ncomps + c);
        const auto *pg_1_31 = buffer.data(pg_1 + 31 * ncomps + c);
        const auto *pg_1_32 = buffer.data(pg_1 + 32 * ncomps + c);
        const auto *pg_1_33 = buffer.data(pg_1 + 33 * ncomps + c);
        const auto *pg_1_34 = buffer.data(pg_1 + 34 * ncomps + c);
        const auto *pg_1_35 = buffer.data(pg_1 + 35 * ncomps + c);
        const auto *pg_1_36 = buffer.data(pg_1 + 36 * ncomps + c);
        const auto *pg_1_37 = buffer.data(pg_1 + 37 * ncomps + c);
        const auto *pg_1_38 = buffer.data(pg_1 + 38 * ncomps + c);
        const auto *pg_1_39 = buffer.data(pg_1 + 39 * ncomps + c);
        const auto *pg_1_40 = buffer.data(pg_1 + 40 * ncomps + c);
        const auto *pg_1_41 = buffer.data(pg_1 + 41 * ncomps + c);
        const auto *pg_1_42 = buffer.data(pg_1 + 42 * ncomps + c);
        const auto *pg_1_43 = buffer.data(pg_1 + 43 * ncomps + c);
        const auto *pg_1_44 = buffer.data(pg_1 + 44 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, pf_1_0, pf_1_1, pf_1_2, pf_1_3, \
                         pf_1_4, pg_1_0, pg_1_1, pg_1_2, pg_1_3, \
                         pg_1_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = -ab_x[k] * pf_1_0[k]
                     + pg_1_0[k];

            t_1[k] = -ab_x[k] * pf_1_1[k]
                     + pg_1_1[k];

            t_2[k] = -ab_x[k] * pf_1_2[k]
                     + pg_1_2[k];

            t_3[k] = -ab_x[k] * pf_1_3[k]
                     + pg_1_3[k];

            t_4[k] = -ab_x[k] * pf_1_4[k]
                     + pg_1_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, pf_1_5, pf_1_6, pf_1_7, pf_1_8, \
                         pf_1_9, pg_1_5, pg_1_6, pg_1_7, pg_1_8, \
                         pg_1_9 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = -ab_x[k] * pf_1_5[k]
                     + pg_1_5[k];

            t_6[k] = -ab_x[k] * pf_1_6[k]
                     + pg_1_6[k];

            t_7[k] = -ab_x[k] * pf_1_7[k]
                     + pg_1_7[k];

            t_8[k] = -ab_x[k] * pf_1_8[k]
                     + pg_1_8[k];

            t_9[k] = -ab_x[k] * pf_1_9[k]
                     + pg_1_9[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, pf_1_10, pf_1_11, pf_1_12, \
                         pf_1_13, pf_1_14, pg_1_15, pg_1_16, pg_1_17, pg_1_18, \
                         pg_1_19 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = -ab_x[k] * pf_1_10[k]
                      + pg_1_15[k];

            t_11[k] = -ab_x[k] * pf_1_11[k]
                      + pg_1_16[k];

            t_12[k] = -ab_x[k] * pf_1_12[k]
                      + pg_1_17[k];

            t_13[k] = -ab_x[k] * pf_1_13[k]
                      + pg_1_18[k];

            t_14[k] = -ab_x[k] * pf_1_14[k]
                      + pg_1_19[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, pf_1_15, pf_1_16, pf_1_17, \
                         pf_1_18, pf_1_19, pg_1_20, pg_1_21, pg_1_22, pg_1_23, \
                         pg_1_24 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = -ab_x[k] * pf_1_15[k]
                      + pg_1_20[k];

            t_16[k] = -ab_x[k] * pf_1_16[k]
                      + pg_1_21[k];

            t_17[k] = -ab_x[k] * pf_1_17[k]
                      + pg_1_22[k];

            t_18[k] = -ab_x[k] * pf_1_18[k]
                      + pg_1_23[k];

            t_19[k] = -ab_x[k] * pf_1_19[k]
                      + pg_1_24[k];
        }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, pf_1_20, pf_1_21, pf_1_22, \
                         pf_1_23, pf_1_24, pg_1_30, pg_1_31, pg_1_32, pg_1_33, \
                         pg_1_34 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_20[k] = -ab_x[k] * pf_1_20[k]
                      + pg_1_30[k];

            t_21[k] = -ab_x[k] * pf_1_21[k]
                      + pg_1_31[k];

            t_22[k] = -ab_x[k] * pf_1_22[k]
                      + pg_1_32[k];

            t_23[k] = -ab_x[k] * pf_1_23[k]
                      + pg_1_33[k];

            t_24[k] = -ab_x[k] * pf_1_24[k]
                      + pg_1_34[k];
        }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, pf_1_25, pf_1_26, pf_1_27, \
                         pf_1_28, pf_1_29, pg_1_35, pg_1_36, pg_1_37, pg_1_38, \
                         pg_1_39 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_25[k] = -ab_x[k] * pf_1_25[k]
                      + pg_1_35[k];

            t_26[k] = -ab_x[k] * pf_1_26[k]
                      + pg_1_36[k];

            t_27[k] = -ab_x[k] * pf_1_27[k]
                      + pg_1_37[k];

            t_28[k] = -ab_x[k] * pf_1_28[k]
                      + pg_1_38[k];

            t_29[k] = -ab_x[k] * pf_1_29[k]
                      + pg_1_39[k];
        }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_y, pf_1_10, pf_1_11, pf_1_12, \
                         pf_1_13, pf_1_14, pg_1_16, pg_1_18, pg_1_19, pg_1_21, \
                         pg_1_22 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_30[k] = -ab_y[k] * pf_1_10[k]
                      + pg_1_16[k];

            t_31[k] = -ab_y[k] * pf_1_11[k]
                      + pg_1_18[k];

            t_32[k] = -ab_y[k] * pf_1_12[k]
                      + pg_1_19[k];

            t_33[k] = -ab_y[k] * pf_1_13[k]
                      + pg_1_21[k];

            t_34[k] = -ab_y[k] * pf_1_14[k]
                      + pg_1_22[k];
        }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_y, pf_1_15, pf_1_16, pf_1_17, \
                         pf_1_18, pf_1_19, pg_1_23, pg_1_25, pg_1_26, pg_1_27, \
                         pg_1_28 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_35[k] = -ab_y[k] * pf_1_15[k]
                      + pg_1_23[k];

            t_36[k] = -ab_y[k] * pf_1_16[k]
                      + pg_1_25[k];

            t_37[k] = -ab_y[k] * pf_1_17[k]
                      + pg_1_26[k];

            t_38[k] = -ab_y[k] * pf_1_18[k]
                      + pg_1_27[k];

            t_39[k] = -ab_y[k] * pf_1_19[k]
                      + pg_1_28[k];
        }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_y, pf_1_20, pf_1_21, pf_1_22, \
                         pf_1_23, pf_1_24, pg_1_31, pg_1_33, pg_1_34, pg_1_36, \
                         pg_1_37 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_40[k] = -ab_y[k] * pf_1_20[k]
                      + pg_1_31[k];

            t_41[k] = -ab_y[k] * pf_1_21[k]
                      + pg_1_33[k];

            t_42[k] = -ab_y[k] * pf_1_22[k]
                      + pg_1_34[k];

            t_43[k] = -ab_y[k] * pf_1_23[k]
                      + pg_1_36[k];

            t_44[k] = -ab_y[k] * pf_1_24[k]
                      + pg_1_37[k];
        }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_y, pf_1_25, pf_1_26, pf_1_27, \
                         pf_1_28, pf_1_29, pg_1_38, pg_1_40, pg_1_41, pg_1_42, \
                         pg_1_43 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_45[k] = -ab_y[k] * pf_1_25[k]
                      + pg_1_38[k];

            t_46[k] = -ab_y[k] * pf_1_26[k]
                      + pg_1_40[k];

            t_47[k] = -ab_y[k] * pf_1_27[k]
                      + pg_1_41[k];

            t_48[k] = -ab_y[k] * pf_1_28[k]
                      + pg_1_42[k];

            t_49[k] = -ab_y[k] * pf_1_29[k]
                      + pg_1_43[k];
        }

#pragma omp simd aligned(t_50, t_51, t_52, ab_z, pf_1_20, pf_1_21, pf_1_22, pf_0_20, pf_0_21, \
                         pf_0_22, pg_1_32, pg_1_34, pg_1_35 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_50[k] = -ab_z[k] * pf_1_20[k]
                      + pf_0_20[k]
                      + pg_1_32[k];

            t_51[k] = -ab_z[k] * pf_1_21[k]
                      + pf_0_21[k]
                      + pg_1_34[k];

            t_52[k] = -ab_z[k] * pf_1_22[k]
                      + pf_0_22[k]
                      + pg_1_35[k];
        }

#pragma omp simd aligned(t_53, t_54, t_55, ab_z, pf_1_23, pf_1_24, pf_1_25, pf_0_23, pf_0_24, \
                         pf_0_25, pg_1_37, pg_1_38, pg_1_39 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_53[k] = -ab_z[k] * pf_1_23[k]
                      + pf_0_23[k]
                      + pg_1_37[k];

            t_54[k] = -ab_z[k] * pf_1_24[k]
                      + pf_0_24[k]
                      + pg_1_38[k];

            t_55[k] = -ab_z[k] * pf_1_25[k]
                      + pf_0_25[k]
                      + pg_1_39[k];
        }

#pragma omp simd aligned(t_56, t_57, t_58, ab_z, pf_1_26, pf_1_27, pf_1_28, pf_0_26, pf_0_27, \
                         pf_0_28, pg_1_41, pg_1_42, pg_1_43 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_56[k] = -ab_z[k] * pf_1_26[k]
                      + pf_0_26[k]
                      + pg_1_41[k];

            t_57[k] = -ab_z[k] * pf_1_27[k]
                      + pf_0_27[k]
                      + pg_1_42[k];

            t_58[k] = -ab_z[k] * pf_1_28[k]
                      + pf_0_28[k]
                      + pg_1_43[k];
        }

#pragma omp simd aligned(t_59, ab_z, pf_1_29, pf_0_29, pg_1_44 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_59[k] = -ab_z[k] * pf_1_29[k]
                      + pf_0_29[k]
                      + pg_1_44[k];
        }
    }
}

}  // namespace simdtrf
