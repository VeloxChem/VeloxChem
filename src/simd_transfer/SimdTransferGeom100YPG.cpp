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


#include "SimdTransferGeom100YPG.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_geom_100y_pg_out_of_first(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                      const size_t target, const size_t pf_1, const size_t pf_0,
                                      const size_t df_1, const size_t ncomps,
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

        const auto *pf_0_6 = buffer.data(pf_0 + 6 * ncomps + c);
        const auto *pf_0_7 = buffer.data(pf_0 + 7 * ncomps + c);
        const auto *pf_0_8 = buffer.data(pf_0 + 8 * ncomps + c);
        const auto *pf_0_9 = buffer.data(pf_0 + 9 * ncomps + c);
        const auto *pf_0_16 = buffer.data(pf_0 + 16 * ncomps + c);
        const auto *pf_0_17 = buffer.data(pf_0 + 17 * ncomps + c);
        const auto *pf_0_18 = buffer.data(pf_0 + 18 * ncomps + c);
        const auto *pf_0_19 = buffer.data(pf_0 + 19 * ncomps + c);
        const auto *pf_0_26 = buffer.data(pf_0 + 26 * ncomps + c);
        const auto *pf_0_27 = buffer.data(pf_0 + 27 * ncomps + c);
        const auto *pf_0_28 = buffer.data(pf_0 + 28 * ncomps + c);
        const auto *pf_0_29 = buffer.data(pf_0 + 29 * ncomps + c);

        const auto *df_1_0 = buffer.data(df_1 + 0 * ncomps + c);
        const auto *df_1_1 = buffer.data(df_1 + 1 * ncomps + c);
        const auto *df_1_2 = buffer.data(df_1 + 2 * ncomps + c);
        const auto *df_1_3 = buffer.data(df_1 + 3 * ncomps + c);
        const auto *df_1_4 = buffer.data(df_1 + 4 * ncomps + c);
        const auto *df_1_5 = buffer.data(df_1 + 5 * ncomps + c);
        const auto *df_1_6 = buffer.data(df_1 + 6 * ncomps + c);
        const auto *df_1_7 = buffer.data(df_1 + 7 * ncomps + c);
        const auto *df_1_8 = buffer.data(df_1 + 8 * ncomps + c);
        const auto *df_1_9 = buffer.data(df_1 + 9 * ncomps + c);
        const auto *df_1_10 = buffer.data(df_1 + 10 * ncomps + c);
        const auto *df_1_11 = buffer.data(df_1 + 11 * ncomps + c);
        const auto *df_1_12 = buffer.data(df_1 + 12 * ncomps + c);
        const auto *df_1_13 = buffer.data(df_1 + 13 * ncomps + c);
        const auto *df_1_14 = buffer.data(df_1 + 14 * ncomps + c);
        const auto *df_1_15 = buffer.data(df_1 + 15 * ncomps + c);
        const auto *df_1_16 = buffer.data(df_1 + 16 * ncomps + c);
        const auto *df_1_17 = buffer.data(df_1 + 17 * ncomps + c);
        const auto *df_1_18 = buffer.data(df_1 + 18 * ncomps + c);
        const auto *df_1_19 = buffer.data(df_1 + 19 * ncomps + c);
        const auto *df_1_20 = buffer.data(df_1 + 20 * ncomps + c);
        const auto *df_1_21 = buffer.data(df_1 + 21 * ncomps + c);
        const auto *df_1_22 = buffer.data(df_1 + 22 * ncomps + c);
        const auto *df_1_23 = buffer.data(df_1 + 23 * ncomps + c);
        const auto *df_1_24 = buffer.data(df_1 + 24 * ncomps + c);
        const auto *df_1_25 = buffer.data(df_1 + 25 * ncomps + c);
        const auto *df_1_26 = buffer.data(df_1 + 26 * ncomps + c);
        const auto *df_1_27 = buffer.data(df_1 + 27 * ncomps + c);
        const auto *df_1_28 = buffer.data(df_1 + 28 * ncomps + c);
        const auto *df_1_29 = buffer.data(df_1 + 29 * ncomps + c);
        const auto *df_1_36 = buffer.data(df_1 + 36 * ncomps + c);
        const auto *df_1_37 = buffer.data(df_1 + 37 * ncomps + c);
        const auto *df_1_38 = buffer.data(df_1 + 38 * ncomps + c);
        const auto *df_1_39 = buffer.data(df_1 + 39 * ncomps + c);
        const auto *df_1_46 = buffer.data(df_1 + 46 * ncomps + c);
        const auto *df_1_47 = buffer.data(df_1 + 47 * ncomps + c);
        const auto *df_1_48 = buffer.data(df_1 + 48 * ncomps + c);
        const auto *df_1_49 = buffer.data(df_1 + 49 * ncomps + c);
        const auto *df_1_59 = buffer.data(df_1 + 59 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, pf_1_0, pf_1_1, pf_1_2, pf_1_3, \
                         pf_1_4, df_1_0, df_1_1, df_1_2, df_1_3, \
                         df_1_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * pf_1_0[k]
                     + df_1_0[k];

            t_1[k] = ab_x[k] * pf_1_1[k]
                     + df_1_1[k];

            t_2[k] = ab_x[k] * pf_1_2[k]
                     + df_1_2[k];

            t_3[k] = ab_x[k] * pf_1_3[k]
                     + df_1_3[k];

            t_4[k] = ab_x[k] * pf_1_4[k]
                     + df_1_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, pf_1_5, pf_1_6, pf_1_7, pf_1_8, \
                         pf_1_9, df_1_5, df_1_6, df_1_7, df_1_8, \
                         df_1_9 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = ab_x[k] * pf_1_5[k]
                     + df_1_5[k];

            t_6[k] = ab_x[k] * pf_1_6[k]
                     + df_1_6[k];

            t_7[k] = ab_x[k] * pf_1_7[k]
                     + df_1_7[k];

            t_8[k] = ab_x[k] * pf_1_8[k]
                     + df_1_8[k];

            t_9[k] = ab_x[k] * pf_1_9[k]
                     + df_1_9[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, ab_y, pf_1_6, pf_1_7, pf_1_8, pf_0_6, pf_0_7, \
                         pf_0_8, df_1_16, df_1_17, df_1_18 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = ab_y[k] * pf_1_6[k]
                      + pf_0_6[k]
                      + df_1_16[k];

            t_11[k] = ab_y[k] * pf_1_7[k]
                      + pf_0_7[k]
                      + df_1_17[k];

            t_12[k] = ab_y[k] * pf_1_8[k]
                      + pf_0_8[k]
                      + df_1_18[k];
        }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, ab_x, ab_y, ab_z, pf_1_9, pf_1_10, pf_1_11, \
                         pf_0_9, df_1_10, df_1_11, df_1_19, df_1_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_13[k] = ab_y[k] * pf_1_9[k]
                      + pf_0_9[k]
                      + df_1_19[k];

            t_14[k] = ab_z[k] * pf_1_9[k]
                      + df_1_29[k];

            t_15[k] = ab_x[k] * pf_1_10[k]
                      + df_1_10[k];

            t_16[k] = ab_x[k] * pf_1_11[k]
                      + df_1_11[k];
        }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, t_21, ab_x, pf_1_12, pf_1_13, pf_1_14, \
                         pf_1_15, pf_1_16, df_1_12, df_1_13, df_1_14, df_1_15, \
                         df_1_16 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_17[k] = ab_x[k] * pf_1_12[k]
                      + df_1_12[k];

            t_18[k] = ab_x[k] * pf_1_13[k]
                      + df_1_13[k];

            t_19[k] = ab_x[k] * pf_1_14[k]
                      + df_1_14[k];

            t_20[k] = ab_x[k] * pf_1_15[k]
                      + df_1_15[k];

            t_21[k] = ab_x[k] * pf_1_16[k]
                      + df_1_16[k];
        }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, ab_x, ab_y, pf_1_16, pf_1_17, pf_1_18, \
                         pf_1_19, pf_0_16, df_1_17, df_1_18, df_1_19, \
                         df_1_36 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_22[k] = ab_x[k] * pf_1_17[k]
                      + df_1_17[k];

            t_23[k] = ab_x[k] * pf_1_18[k]
                      + df_1_18[k];

            t_24[k] = ab_x[k] * pf_1_19[k]
                      + df_1_19[k];

            t_25[k] = ab_y[k] * pf_1_16[k]
                      + pf_0_16[k]
                      + df_1_36[k];
        }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, ab_y, ab_z, pf_1_17, pf_1_18, pf_1_19, \
                         pf_0_17, pf_0_18, pf_0_19, df_1_37, df_1_38, df_1_39, \
                         df_1_49 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_26[k] = ab_y[k] * pf_1_17[k]
                      + pf_0_17[k]
                      + df_1_37[k];

            t_27[k] = ab_y[k] * pf_1_18[k]
                      + pf_0_18[k]
                      + df_1_38[k];

            t_28[k] = ab_y[k] * pf_1_19[k]
                      + pf_0_19[k]
                      + df_1_39[k];

            t_29[k] = ab_z[k] * pf_1_19[k]
                      + df_1_49[k];
        }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, pf_1_20, pf_1_21, pf_1_22, \
                         pf_1_23, pf_1_24, df_1_20, df_1_21, df_1_22, df_1_23, \
                         df_1_24 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_30[k] = ab_x[k] * pf_1_20[k]
                      + df_1_20[k];

            t_31[k] = ab_x[k] * pf_1_21[k]
                      + df_1_21[k];

            t_32[k] = ab_x[k] * pf_1_22[k]
                      + df_1_22[k];

            t_33[k] = ab_x[k] * pf_1_23[k]
                      + df_1_23[k];

            t_34[k] = ab_x[k] * pf_1_24[k]
                      + df_1_24[k];
        }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, pf_1_25, pf_1_26, pf_1_27, \
                         pf_1_28, pf_1_29, df_1_25, df_1_26, df_1_27, df_1_28, \
                         df_1_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_35[k] = ab_x[k] * pf_1_25[k]
                      + df_1_25[k];

            t_36[k] = ab_x[k] * pf_1_26[k]
                      + df_1_26[k];

            t_37[k] = ab_x[k] * pf_1_27[k]
                      + df_1_27[k];

            t_38[k] = ab_x[k] * pf_1_28[k]
                      + df_1_28[k];

            t_39[k] = ab_x[k] * pf_1_29[k]
                      + df_1_29[k];
        }

#pragma omp simd aligned(t_40, t_41, t_42, ab_y, pf_1_26, pf_1_27, pf_1_28, pf_0_26, pf_0_27, \
                         pf_0_28, df_1_46, df_1_47, df_1_48 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_40[k] = ab_y[k] * pf_1_26[k]
                      + pf_0_26[k]
                      + df_1_46[k];

            t_41[k] = ab_y[k] * pf_1_27[k]
                      + pf_0_27[k]
                      + df_1_47[k];

            t_42[k] = ab_y[k] * pf_1_28[k]
                      + pf_0_28[k]
                      + df_1_48[k];
        }

#pragma omp simd aligned(t_43, t_44, ab_y, ab_z, pf_1_29, pf_0_29, df_1_49, \
                         df_1_59 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_43[k] = ab_y[k] * pf_1_29[k]
                      + pf_0_29[k]
                      + df_1_49[k];

            t_44[k] = ab_z[k] * pf_1_29[k]
                      + df_1_59[k];
        }
    }
}

}  // namespace simdtrf
