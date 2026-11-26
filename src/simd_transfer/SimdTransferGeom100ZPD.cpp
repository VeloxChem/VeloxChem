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


#include "SimdTransferGeom100ZPD.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_geom_100z_pd_out_of_first(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                      const size_t target, const size_t pp_1, const size_t pp_0,
                                      const size_t dp_1, const size_t ncomps,
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

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *pp_1_0 = buffer.data(pp_1 + 0 * ncomps + c);
        const auto *pp_1_1 = buffer.data(pp_1 + 1 * ncomps + c);
        const auto *pp_1_2 = buffer.data(pp_1 + 2 * ncomps + c);
        const auto *pp_1_3 = buffer.data(pp_1 + 3 * ncomps + c);
        const auto *pp_1_4 = buffer.data(pp_1 + 4 * ncomps + c);
        const auto *pp_1_5 = buffer.data(pp_1 + 5 * ncomps + c);
        const auto *pp_1_6 = buffer.data(pp_1 + 6 * ncomps + c);
        const auto *pp_1_7 = buffer.data(pp_1 + 7 * ncomps + c);
        const auto *pp_1_8 = buffer.data(pp_1 + 8 * ncomps + c);

        const auto *pp_0_2 = buffer.data(pp_0 + 2 * ncomps + c);
        const auto *pp_0_5 = buffer.data(pp_0 + 5 * ncomps + c);
        const auto *pp_0_8 = buffer.data(pp_0 + 8 * ncomps + c);

        const auto *dp_1_0 = buffer.data(dp_1 + 0 * ncomps + c);
        const auto *dp_1_1 = buffer.data(dp_1 + 1 * ncomps + c);
        const auto *dp_1_2 = buffer.data(dp_1 + 2 * ncomps + c);
        const auto *dp_1_3 = buffer.data(dp_1 + 3 * ncomps + c);
        const auto *dp_1_4 = buffer.data(dp_1 + 4 * ncomps + c);
        const auto *dp_1_5 = buffer.data(dp_1 + 5 * ncomps + c);
        const auto *dp_1_6 = buffer.data(dp_1 + 6 * ncomps + c);
        const auto *dp_1_7 = buffer.data(dp_1 + 7 * ncomps + c);
        const auto *dp_1_8 = buffer.data(dp_1 + 8 * ncomps + c);
        const auto *dp_1_10 = buffer.data(dp_1 + 10 * ncomps + c);
        const auto *dp_1_11 = buffer.data(dp_1 + 11 * ncomps + c);
        const auto *dp_1_13 = buffer.data(dp_1 + 13 * ncomps + c);
        const auto *dp_1_14 = buffer.data(dp_1 + 14 * ncomps + c);
        const auto *dp_1_17 = buffer.data(dp_1 + 17 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, ab_y, pp_1_0, pp_1_1, pp_1_2, dp_1_0, \
                         dp_1_1, dp_1_2, dp_1_4, dp_1_5 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * pp_1_0[k]
                     + dp_1_0[k];

            t_1[k] = ab_x[k] * pp_1_1[k]
                     + dp_1_1[k];

            t_2[k] = ab_x[k] * pp_1_2[k]
                     + dp_1_2[k];

            t_3[k] = ab_y[k] * pp_1_1[k]
                     + dp_1_4[k];

            t_4[k] = ab_y[k] * pp_1_2[k]
                     + dp_1_5[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, ab_x, ab_z, pp_1_2, pp_1_3, pp_1_4, pp_1_5, \
                         pp_0_2, dp_1_3, dp_1_4, dp_1_5, dp_1_8 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = ab_z[k] * pp_1_2[k]
                     + pp_0_2[k]
                     + dp_1_8[k];

            t_6[k] = ab_x[k] * pp_1_3[k]
                     + dp_1_3[k];

            t_7[k] = ab_x[k] * pp_1_4[k]
                     + dp_1_4[k];

            t_8[k] = ab_x[k] * pp_1_5[k]
                     + dp_1_5[k];
        }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, ab_x, ab_y, ab_z, pp_1_4, pp_1_5, pp_1_6, \
                         pp_0_5, dp_1_6, dp_1_10, dp_1_11, dp_1_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_9[k] = ab_y[k] * pp_1_4[k]
                     + dp_1_10[k];

            t_10[k] = ab_y[k] * pp_1_5[k]
                      + dp_1_11[k];

            t_11[k] = ab_z[k] * pp_1_5[k]
                      + pp_0_5[k]
                      + dp_1_14[k];

            t_12[k] = ab_x[k] * pp_1_6[k]
                      + dp_1_6[k];
        }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, ab_x, ab_y, ab_z, pp_1_7, pp_1_8, \
                         pp_0_8, dp_1_7, dp_1_8, dp_1_13, dp_1_14, \
                         dp_1_17 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_13[k] = ab_x[k] * pp_1_7[k]
                      + dp_1_7[k];

            t_14[k] = ab_x[k] * pp_1_8[k]
                      + dp_1_8[k];

            t_15[k] = ab_y[k] * pp_1_7[k]
                      + dp_1_13[k];

            t_16[k] = ab_y[k] * pp_1_8[k]
                      + dp_1_14[k];

            t_17[k] = ab_z[k] * pp_1_8[k]
                      + pp_0_8[k]
                      + dp_1_17[k];
        }
    }
}

}  // namespace simdtrf
