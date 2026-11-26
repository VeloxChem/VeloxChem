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


#include "SimdTransferGeom100XDP.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_geom_100x_dp_out_of_first(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                      const size_t target, const size_t ds_1, const size_t ds_0,
                                      const size_t fs_1, const size_t ncomps,
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

        const auto *ds_1_0 = buffer.data(ds_1 + 0 * ncomps + c);
        const auto *ds_1_1 = buffer.data(ds_1 + 1 * ncomps + c);
        const auto *ds_1_2 = buffer.data(ds_1 + 2 * ncomps + c);
        const auto *ds_1_3 = buffer.data(ds_1 + 3 * ncomps + c);
        const auto *ds_1_4 = buffer.data(ds_1 + 4 * ncomps + c);
        const auto *ds_1_5 = buffer.data(ds_1 + 5 * ncomps + c);

        const auto *ds_0_0 = buffer.data(ds_0 + 0 * ncomps + c);
        const auto *ds_0_1 = buffer.data(ds_0 + 1 * ncomps + c);
        const auto *ds_0_2 = buffer.data(ds_0 + 2 * ncomps + c);
        const auto *ds_0_3 = buffer.data(ds_0 + 3 * ncomps + c);
        const auto *ds_0_4 = buffer.data(ds_0 + 4 * ncomps + c);
        const auto *ds_0_5 = buffer.data(ds_0 + 5 * ncomps + c);

        const auto *fs_1_0 = buffer.data(fs_1 + 0 * ncomps + c);
        const auto *fs_1_1 = buffer.data(fs_1 + 1 * ncomps + c);
        const auto *fs_1_2 = buffer.data(fs_1 + 2 * ncomps + c);
        const auto *fs_1_3 = buffer.data(fs_1 + 3 * ncomps + c);
        const auto *fs_1_4 = buffer.data(fs_1 + 4 * ncomps + c);
        const auto *fs_1_5 = buffer.data(fs_1 + 5 * ncomps + c);
        const auto *fs_1_6 = buffer.data(fs_1 + 6 * ncomps + c);
        const auto *fs_1_7 = buffer.data(fs_1 + 7 * ncomps + c);
        const auto *fs_1_8 = buffer.data(fs_1 + 8 * ncomps + c);
        const auto *fs_1_9 = buffer.data(fs_1 + 9 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, ab_y, ab_z, ds_1_0, ds_1_1, ds_0_0, \
                         ds_0_1, fs_1_0, fs_1_1, fs_1_2, fs_1_3 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * ds_1_0[k]
                     + ds_0_0[k]
                     + fs_1_0[k];

            t_1[k] = ab_y[k] * ds_1_0[k]
                     + fs_1_1[k];

            t_2[k] = ab_z[k] * ds_1_0[k]
                     + fs_1_2[k];

            t_3[k] = ab_x[k] * ds_1_1[k]
                     + ds_0_1[k]
                     + fs_1_1[k];

            t_4[k] = ab_y[k] * ds_1_1[k]
                     + fs_1_3[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, ab_x, ab_y, ab_z, ds_1_1, ds_1_2, ds_0_2, fs_1_2, \
                         fs_1_4, fs_1_5 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = ab_z[k] * ds_1_1[k]
                     + fs_1_4[k];

            t_6[k] = ab_x[k] * ds_1_2[k]
                     + ds_0_2[k]
                     + fs_1_2[k];

            t_7[k] = ab_y[k] * ds_1_2[k]
                     + fs_1_4[k];

            t_8[k] = ab_z[k] * ds_1_2[k]
                     + fs_1_5[k];
        }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, t_13, ab_x, ab_y, ab_z, ds_1_3, ds_1_4, \
                         ds_0_3, ds_0_4, fs_1_3, fs_1_4, fs_1_6, \
                         fs_1_7 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_9[k] = ab_x[k] * ds_1_3[k]
                     + ds_0_3[k]
                     + fs_1_3[k];

            t_10[k] = ab_y[k] * ds_1_3[k]
                      + fs_1_6[k];

            t_11[k] = ab_z[k] * ds_1_3[k]
                      + fs_1_7[k];

            t_12[k] = ab_x[k] * ds_1_4[k]
                      + ds_0_4[k]
                      + fs_1_4[k];

            t_13[k] = ab_y[k] * ds_1_4[k]
                      + fs_1_7[k];
        }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, ab_x, ab_y, ab_z, ds_1_4, ds_1_5, ds_0_5, \
                         fs_1_5, fs_1_8, fs_1_9 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_14[k] = ab_z[k] * ds_1_4[k]
                      + fs_1_8[k];

            t_15[k] = ab_x[k] * ds_1_5[k]
                      + ds_0_5[k]
                      + fs_1_5[k];

            t_16[k] = ab_y[k] * ds_1_5[k]
                      + fs_1_8[k];

            t_17[k] = ab_z[k] * ds_1_5[k]
                      + fs_1_9[k];
        }
    }
}

}  // namespace simdtrf
