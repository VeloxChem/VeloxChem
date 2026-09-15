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


#include "SimdTransferPF.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_pf(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t sf, const size_t sg, const size_t ncomps, const size_t nmax) -> void
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

        const auto *sf_0 = buffer.data(sf + 0 * ncomps + c);
        const auto *sf_1 = buffer.data(sf + 1 * ncomps + c);
        const auto *sf_2 = buffer.data(sf + 2 * ncomps + c);
        const auto *sf_3 = buffer.data(sf + 3 * ncomps + c);
        const auto *sf_4 = buffer.data(sf + 4 * ncomps + c);
        const auto *sf_5 = buffer.data(sf + 5 * ncomps + c);
        const auto *sf_6 = buffer.data(sf + 6 * ncomps + c);
        const auto *sf_7 = buffer.data(sf + 7 * ncomps + c);
        const auto *sf_8 = buffer.data(sf + 8 * ncomps + c);
        const auto *sf_9 = buffer.data(sf + 9 * ncomps + c);

        const auto *sg_0 = buffer.data(sg + 0 * ncomps + c);
        const auto *sg_1 = buffer.data(sg + 1 * ncomps + c);
        const auto *sg_2 = buffer.data(sg + 2 * ncomps + c);
        const auto *sg_3 = buffer.data(sg + 3 * ncomps + c);
        const auto *sg_4 = buffer.data(sg + 4 * ncomps + c);
        const auto *sg_5 = buffer.data(sg + 5 * ncomps + c);
        const auto *sg_6 = buffer.data(sg + 6 * ncomps + c);
        const auto *sg_7 = buffer.data(sg + 7 * ncomps + c);
        const auto *sg_8 = buffer.data(sg + 8 * ncomps + c);
        const auto *sg_9 = buffer.data(sg + 9 * ncomps + c);
        const auto *sg_10 = buffer.data(sg + 10 * ncomps + c);
        const auto *sg_11 = buffer.data(sg + 11 * ncomps + c);
        const auto *sg_12 = buffer.data(sg + 12 * ncomps + c);
        const auto *sg_13 = buffer.data(sg + 13 * ncomps + c);
        const auto *sg_14 = buffer.data(sg + 14 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, sf_0, sf_1, sf_2, sf_3, sf_4, sg_0, \
                         sg_1, sg_2, sg_3, sg_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = -ab_x[k] * sf_0[k]
                     + sg_0[k];

            t_1[k] = -ab_x[k] * sf_1[k]
                     + sg_1[k];

            t_2[k] = -ab_x[k] * sf_2[k]
                     + sg_2[k];

            t_3[k] = -ab_x[k] * sf_3[k]
                     + sg_3[k];

            t_4[k] = -ab_x[k] * sf_4[k]
                     + sg_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, sf_5, sf_6, sf_7, sf_8, sf_9, sg_5, \
                         sg_6, sg_7, sg_8, sg_9 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = -ab_x[k] * sf_5[k]
                     + sg_5[k];

            t_6[k] = -ab_x[k] * sf_6[k]
                     + sg_6[k];

            t_7[k] = -ab_x[k] * sf_7[k]
                     + sg_7[k];

            t_8[k] = -ab_x[k] * sf_8[k]
                     + sg_8[k];

            t_9[k] = -ab_x[k] * sf_9[k]
                     + sg_9[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_y, sf_0, sf_1, sf_2, sf_3, sf_4, \
                         sg_1, sg_3, sg_4, sg_6, sg_7 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = -ab_y[k] * sf_0[k]
                      + sg_1[k];

            t_11[k] = -ab_y[k] * sf_1[k]
                      + sg_3[k];

            t_12[k] = -ab_y[k] * sf_2[k]
                      + sg_4[k];

            t_13[k] = -ab_y[k] * sf_3[k]
                      + sg_6[k];

            t_14[k] = -ab_y[k] * sf_4[k]
                      + sg_7[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_y, sf_5, sf_6, sf_7, sf_8, sf_9, \
                         sg_8, sg_10, sg_11, sg_12, sg_13 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = -ab_y[k] * sf_5[k]
                      + sg_8[k];

            t_16[k] = -ab_y[k] * sf_6[k]
                      + sg_10[k];

            t_17[k] = -ab_y[k] * sf_7[k]
                      + sg_11[k];

            t_18[k] = -ab_y[k] * sf_8[k]
                      + sg_12[k];

            t_19[k] = -ab_y[k] * sf_9[k]
                      + sg_13[k];
        }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_z, sf_0, sf_1, sf_2, sf_3, sf_4, \
                         sg_2, sg_4, sg_5, sg_7, sg_8 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_20[k] = -ab_z[k] * sf_0[k]
                      + sg_2[k];

            t_21[k] = -ab_z[k] * sf_1[k]
                      + sg_4[k];

            t_22[k] = -ab_z[k] * sf_2[k]
                      + sg_5[k];

            t_23[k] = -ab_z[k] * sf_3[k]
                      + sg_7[k];

            t_24[k] = -ab_z[k] * sf_4[k]
                      + sg_8[k];
        }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_z, sf_5, sf_6, sf_7, sf_8, sf_9, \
                         sg_9, sg_11, sg_12, sg_13, sg_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_25[k] = -ab_z[k] * sf_5[k]
                      + sg_9[k];

            t_26[k] = -ab_z[k] * sf_6[k]
                      + sg_11[k];

            t_27[k] = -ab_z[k] * sf_7[k]
                      + sg_12[k];

            t_28[k] = -ab_z[k] * sf_8[k]
                      + sg_13[k];

            t_29[k] = -ab_z[k] * sf_9[k]
                      + sg_14[k];
        }
    }
}

auto
compute_hrr_pf_out_of_first(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                            const size_t target, const size_t pd, const size_t dd,
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

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *pd_0 = buffer.data(pd + 0 * ncomps + c);
        const auto *pd_1 = buffer.data(pd + 1 * ncomps + c);
        const auto *pd_2 = buffer.data(pd + 2 * ncomps + c);
        const auto *pd_3 = buffer.data(pd + 3 * ncomps + c);
        const auto *pd_4 = buffer.data(pd + 4 * ncomps + c);
        const auto *pd_5 = buffer.data(pd + 5 * ncomps + c);
        const auto *pd_6 = buffer.data(pd + 6 * ncomps + c);
        const auto *pd_7 = buffer.data(pd + 7 * ncomps + c);
        const auto *pd_8 = buffer.data(pd + 8 * ncomps + c);
        const auto *pd_9 = buffer.data(pd + 9 * ncomps + c);
        const auto *pd_10 = buffer.data(pd + 10 * ncomps + c);
        const auto *pd_11 = buffer.data(pd + 11 * ncomps + c);
        const auto *pd_12 = buffer.data(pd + 12 * ncomps + c);
        const auto *pd_13 = buffer.data(pd + 13 * ncomps + c);
        const auto *pd_14 = buffer.data(pd + 14 * ncomps + c);
        const auto *pd_15 = buffer.data(pd + 15 * ncomps + c);
        const auto *pd_16 = buffer.data(pd + 16 * ncomps + c);
        const auto *pd_17 = buffer.data(pd + 17 * ncomps + c);

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
        const auto *dd_21 = buffer.data(dd + 21 * ncomps + c);
        const auto *dd_22 = buffer.data(dd + 22 * ncomps + c);
        const auto *dd_23 = buffer.data(dd + 23 * ncomps + c);
        const auto *dd_27 = buffer.data(dd + 27 * ncomps + c);
        const auto *dd_28 = buffer.data(dd + 28 * ncomps + c);
        const auto *dd_29 = buffer.data(dd + 29 * ncomps + c);
        const auto *dd_35 = buffer.data(dd + 35 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, pd_0, pd_1, pd_2, pd_3, pd_4, dd_0, \
                         dd_1, dd_2, dd_3, dd_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * pd_0[k]
                     + dd_0[k];

            t_1[k] = ab_x[k] * pd_1[k]
                     + dd_1[k];

            t_2[k] = ab_x[k] * pd_2[k]
                     + dd_2[k];

            t_3[k] = ab_x[k] * pd_3[k]
                     + dd_3[k];

            t_4[k] = ab_x[k] * pd_4[k]
                     + dd_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, ab_y, ab_z, pd_3, pd_4, pd_5, dd_5, \
                         dd_9, dd_10, dd_11, dd_17 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = ab_x[k] * pd_5[k]
                     + dd_5[k];

            t_6[k] = ab_y[k] * pd_3[k]
                     + dd_9[k];

            t_7[k] = ab_y[k] * pd_4[k]
                     + dd_10[k];

            t_8[k] = ab_y[k] * pd_5[k]
                     + dd_11[k];

            t_9[k] = ab_z[k] * pd_5[k]
                     + dd_17[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, pd_6, pd_7, pd_8, pd_9, pd_10, \
                         dd_6, dd_7, dd_8, dd_9, dd_10 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = ab_x[k] * pd_6[k]
                      + dd_6[k];

            t_11[k] = ab_x[k] * pd_7[k]
                      + dd_7[k];

            t_12[k] = ab_x[k] * pd_8[k]
                      + dd_8[k];

            t_13[k] = ab_x[k] * pd_9[k]
                      + dd_9[k];

            t_14[k] = ab_x[k] * pd_10[k]
                      + dd_10[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, ab_y, ab_z, pd_9, pd_10, pd_11, \
                         dd_11, dd_21, dd_22, dd_23, dd_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = ab_x[k] * pd_11[k]
                      + dd_11[k];

            t_16[k] = ab_y[k] * pd_9[k]
                      + dd_21[k];

            t_17[k] = ab_y[k] * pd_10[k]
                      + dd_22[k];

            t_18[k] = ab_y[k] * pd_11[k]
                      + dd_23[k];

            t_19[k] = ab_z[k] * pd_11[k]
                      + dd_29[k];
        }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, pd_12, pd_13, pd_14, pd_15, \
                         pd_16, dd_12, dd_13, dd_14, dd_15, dd_16 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_20[k] = ab_x[k] * pd_12[k]
                      + dd_12[k];

            t_21[k] = ab_x[k] * pd_13[k]
                      + dd_13[k];

            t_22[k] = ab_x[k] * pd_14[k]
                      + dd_14[k];

            t_23[k] = ab_x[k] * pd_15[k]
                      + dd_15[k];

            t_24[k] = ab_x[k] * pd_16[k]
                      + dd_16[k];
        }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, ab_y, ab_z, pd_15, pd_16, pd_17, \
                         dd_17, dd_27, dd_28, dd_29, dd_35 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_25[k] = ab_x[k] * pd_17[k]
                      + dd_17[k];

            t_26[k] = ab_y[k] * pd_15[k]
                      + dd_27[k];

            t_27[k] = ab_y[k] * pd_16[k]
                      + dd_28[k];

            t_28[k] = ab_y[k] * pd_17[k]
                      + dd_29[k];

            t_29[k] = ab_z[k] * pd_17[k]
                      + dd_35[k];
        }
    }
}

}  // namespace simdtrf
