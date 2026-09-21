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


#include "SimdTransferGeom100YPF.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_geom_100y_pf_out_of_first(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                      const size_t target, const size_t pd_1, const size_t pd_0,
                                      const size_t dd_1, const size_t ncomps,
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

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *pd_1_0 = buffer.data(pd_1 + 0 * ncomps + c);
        const auto *pd_1_1 = buffer.data(pd_1 + 1 * ncomps + c);
        const auto *pd_1_2 = buffer.data(pd_1 + 2 * ncomps + c);
        const auto *pd_1_3 = buffer.data(pd_1 + 3 * ncomps + c);
        const auto *pd_1_4 = buffer.data(pd_1 + 4 * ncomps + c);
        const auto *pd_1_5 = buffer.data(pd_1 + 5 * ncomps + c);
        const auto *pd_1_6 = buffer.data(pd_1 + 6 * ncomps + c);
        const auto *pd_1_7 = buffer.data(pd_1 + 7 * ncomps + c);
        const auto *pd_1_8 = buffer.data(pd_1 + 8 * ncomps + c);
        const auto *pd_1_9 = buffer.data(pd_1 + 9 * ncomps + c);
        const auto *pd_1_10 = buffer.data(pd_1 + 10 * ncomps + c);
        const auto *pd_1_11 = buffer.data(pd_1 + 11 * ncomps + c);
        const auto *pd_1_12 = buffer.data(pd_1 + 12 * ncomps + c);
        const auto *pd_1_13 = buffer.data(pd_1 + 13 * ncomps + c);
        const auto *pd_1_14 = buffer.data(pd_1 + 14 * ncomps + c);
        const auto *pd_1_15 = buffer.data(pd_1 + 15 * ncomps + c);
        const auto *pd_1_16 = buffer.data(pd_1 + 16 * ncomps + c);
        const auto *pd_1_17 = buffer.data(pd_1 + 17 * ncomps + c);

        const auto *pd_0_3 = buffer.data(pd_0 + 3 * ncomps + c);
        const auto *pd_0_4 = buffer.data(pd_0 + 4 * ncomps + c);
        const auto *pd_0_5 = buffer.data(pd_0 + 5 * ncomps + c);
        const auto *pd_0_9 = buffer.data(pd_0 + 9 * ncomps + c);
        const auto *pd_0_10 = buffer.data(pd_0 + 10 * ncomps + c);
        const auto *pd_0_11 = buffer.data(pd_0 + 11 * ncomps + c);
        const auto *pd_0_15 = buffer.data(pd_0 + 15 * ncomps + c);
        const auto *pd_0_16 = buffer.data(pd_0 + 16 * ncomps + c);
        const auto *pd_0_17 = buffer.data(pd_0 + 17 * ncomps + c);

        const auto *dd_1_0 = buffer.data(dd_1 + 0 * ncomps + c);
        const auto *dd_1_1 = buffer.data(dd_1 + 1 * ncomps + c);
        const auto *dd_1_2 = buffer.data(dd_1 + 2 * ncomps + c);
        const auto *dd_1_3 = buffer.data(dd_1 + 3 * ncomps + c);
        const auto *dd_1_4 = buffer.data(dd_1 + 4 * ncomps + c);
        const auto *dd_1_5 = buffer.data(dd_1 + 5 * ncomps + c);
        const auto *dd_1_6 = buffer.data(dd_1 + 6 * ncomps + c);
        const auto *dd_1_7 = buffer.data(dd_1 + 7 * ncomps + c);
        const auto *dd_1_8 = buffer.data(dd_1 + 8 * ncomps + c);
        const auto *dd_1_9 = buffer.data(dd_1 + 9 * ncomps + c);
        const auto *dd_1_10 = buffer.data(dd_1 + 10 * ncomps + c);
        const auto *dd_1_11 = buffer.data(dd_1 + 11 * ncomps + c);
        const auto *dd_1_12 = buffer.data(dd_1 + 12 * ncomps + c);
        const auto *dd_1_13 = buffer.data(dd_1 + 13 * ncomps + c);
        const auto *dd_1_14 = buffer.data(dd_1 + 14 * ncomps + c);
        const auto *dd_1_15 = buffer.data(dd_1 + 15 * ncomps + c);
        const auto *dd_1_16 = buffer.data(dd_1 + 16 * ncomps + c);
        const auto *dd_1_17 = buffer.data(dd_1 + 17 * ncomps + c);
        const auto *dd_1_21 = buffer.data(dd_1 + 21 * ncomps + c);
        const auto *dd_1_22 = buffer.data(dd_1 + 22 * ncomps + c);
        const auto *dd_1_23 = buffer.data(dd_1 + 23 * ncomps + c);
        const auto *dd_1_27 = buffer.data(dd_1 + 27 * ncomps + c);
        const auto *dd_1_28 = buffer.data(dd_1 + 28 * ncomps + c);
        const auto *dd_1_29 = buffer.data(dd_1 + 29 * ncomps + c);
        const auto *dd_1_35 = buffer.data(dd_1 + 35 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, pd_1_0, pd_1_1, pd_1_2, pd_1_3, \
                         pd_1_4, dd_1_0, dd_1_1, dd_1_2, dd_1_3, \
                         dd_1_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * pd_1_0[k]
                     + dd_1_0[k];

            t_1[k] = ab_x[k] * pd_1_1[k]
                     + dd_1_1[k];

            t_2[k] = ab_x[k] * pd_1_2[k]
                     + dd_1_2[k];

            t_3[k] = ab_x[k] * pd_1_3[k]
                     + dd_1_3[k];

            t_4[k] = ab_x[k] * pd_1_4[k]
                     + dd_1_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, ab_x, ab_y, pd_1_3, pd_1_4, pd_1_5, pd_0_3, \
                         pd_0_4, pd_0_5, dd_1_5, dd_1_9, dd_1_10, \
                         dd_1_11 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = ab_x[k] * pd_1_5[k]
                     + dd_1_5[k];

            t_6[k] = ab_y[k] * pd_1_3[k]
                     + pd_0_3[k]
                     + dd_1_9[k];

            t_7[k] = ab_y[k] * pd_1_4[k]
                     + pd_0_4[k]
                     + dd_1_10[k];

            t_8[k] = ab_y[k] * pd_1_5[k]
                     + pd_0_5[k]
                     + dd_1_11[k];
        }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, ab_x, ab_z, pd_1_5, pd_1_6, pd_1_7, pd_1_8, \
                         dd_1_6, dd_1_7, dd_1_8, dd_1_17 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_9[k] = ab_z[k] * pd_1_5[k]
                     + dd_1_17[k];

            t_10[k] = ab_x[k] * pd_1_6[k]
                      + dd_1_6[k];

            t_11[k] = ab_x[k] * pd_1_7[k]
                      + dd_1_7[k];

            t_12[k] = ab_x[k] * pd_1_8[k]
                      + dd_1_8[k];
        }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, ab_x, ab_y, pd_1_9, pd_1_10, pd_1_11, pd_0_9, \
                         dd_1_9, dd_1_10, dd_1_11, dd_1_21 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_13[k] = ab_x[k] * pd_1_9[k]
                      + dd_1_9[k];

            t_14[k] = ab_x[k] * pd_1_10[k]
                      + dd_1_10[k];

            t_15[k] = ab_x[k] * pd_1_11[k]
                      + dd_1_11[k];

            t_16[k] = ab_y[k] * pd_1_9[k]
                      + pd_0_9[k]
                      + dd_1_21[k];
        }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, ab_x, ab_y, ab_z, pd_1_10, pd_1_11, pd_1_12, \
                         pd_0_10, pd_0_11, dd_1_12, dd_1_22, dd_1_23, \
                         dd_1_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_17[k] = ab_y[k] * pd_1_10[k]
                      + pd_0_10[k]
                      + dd_1_22[k];

            t_18[k] = ab_y[k] * pd_1_11[k]
                      + pd_0_11[k]
                      + dd_1_23[k];

            t_19[k] = ab_z[k] * pd_1_11[k]
                      + dd_1_29[k];

            t_20[k] = ab_x[k] * pd_1_12[k]
                      + dd_1_12[k];
        }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, ab_x, pd_1_13, pd_1_14, pd_1_15, \
                         pd_1_16, pd_1_17, dd_1_13, dd_1_14, dd_1_15, dd_1_16, \
                         dd_1_17 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_21[k] = ab_x[k] * pd_1_13[k]
                      + dd_1_13[k];

            t_22[k] = ab_x[k] * pd_1_14[k]
                      + dd_1_14[k];

            t_23[k] = ab_x[k] * pd_1_15[k]
                      + dd_1_15[k];

            t_24[k] = ab_x[k] * pd_1_16[k]
                      + dd_1_16[k];

            t_25[k] = ab_x[k] * pd_1_17[k]
                      + dd_1_17[k];
        }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, ab_y, ab_z, pd_1_15, pd_1_16, pd_1_17, \
                         pd_0_15, pd_0_16, pd_0_17, dd_1_27, dd_1_28, dd_1_29, \
                         dd_1_35 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_26[k] = ab_y[k] * pd_1_15[k]
                      + pd_0_15[k]
                      + dd_1_27[k];

            t_27[k] = ab_y[k] * pd_1_16[k]
                      + pd_0_16[k]
                      + dd_1_28[k];

            t_28[k] = ab_y[k] * pd_1_17[k]
                      + pd_0_17[k]
                      + dd_1_29[k];

            t_29[k] = ab_z[k] * pd_1_17[k]
                      + dd_1_35[k];
        }
    }
}

}  // namespace simdtrf
