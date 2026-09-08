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


#include "SimdTransferDP.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_dp(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t ds, const size_t fs, const size_t nmax) -> void
{
    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);
    auto *t_9 = buffer.data(target + 9);
    auto *t_10 = buffer.data(target + 10);
    auto *t_11 = buffer.data(target + 11);
    auto *t_12 = buffer.data(target + 12);
    auto *t_13 = buffer.data(target + 13);
    auto *t_14 = buffer.data(target + 14);
    auto *t_15 = buffer.data(target + 15);
    auto *t_16 = buffer.data(target + 16);
    auto *t_17 = buffer.data(target + 17);

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

    const auto *ds_0 = buffer.data(ds + 0);
    const auto *ds_1 = buffer.data(ds + 1);
    const auto *ds_2 = buffer.data(ds + 2);
    const auto *ds_3 = buffer.data(ds + 3);
    const auto *ds_4 = buffer.data(ds + 4);
    const auto *ds_5 = buffer.data(ds + 5);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_1 = buffer.data(fs + 1);
    const auto *fs_2 = buffer.data(fs + 2);
    const auto *fs_3 = buffer.data(fs + 3);
    const auto *fs_4 = buffer.data(fs + 4);
    const auto *fs_5 = buffer.data(fs + 5);
    const auto *fs_6 = buffer.data(fs + 6);
    const auto *fs_7 = buffer.data(fs + 7);
    const auto *fs_8 = buffer.data(fs + 8);
    const auto *fs_9 = buffer.data(fs + 9);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, ab_x, ab_y, ab_z, ds_0, ds_1, fs_0, \
                         fs_1, fs_2, fs_3, fs_4 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_0[k] = ab_x[k] * ds_0[k]
                 + fs_0[k];

        t_1[k] = ab_y[k] * ds_0[k]
                 + fs_1[k];

        t_2[k] = ab_z[k] * ds_0[k]
                 + fs_2[k];

        t_3[k] = ab_x[k] * ds_1[k]
                 + fs_1[k];

        t_4[k] = ab_y[k] * ds_1[k]
                 + fs_3[k];

        t_5[k] = ab_z[k] * ds_1[k]
                 + fs_4[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, ab_x, ab_y, ab_z, ds_2, ds_3, fs_2, fs_3, \
                         fs_4, fs_5, fs_6 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_6[k] = ab_x[k] * ds_2[k]
                 + fs_2[k];

        t_7[k] = ab_y[k] * ds_2[k]
                 + fs_4[k];

        t_8[k] = ab_z[k] * ds_2[k]
                 + fs_5[k];

        t_9[k] = ab_x[k] * ds_3[k]
                 + fs_3[k];

        t_10[k] = ab_y[k] * ds_3[k]
                  + fs_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, ab_x, ab_y, ab_z, ds_3, ds_4, \
                         ds_5, fs_4, fs_5, fs_7, fs_8 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_11[k] = ab_z[k] * ds_3[k]
                  + fs_7[k];

        t_12[k] = ab_x[k] * ds_4[k]
                  + fs_4[k];

        t_13[k] = ab_y[k] * ds_4[k]
                  + fs_7[k];

        t_14[k] = ab_z[k] * ds_4[k]
                  + fs_8[k];

        t_15[k] = ab_x[k] * ds_5[k]
                  + fs_5[k];

        t_16[k] = ab_y[k] * ds_5[k]
                  + fs_8[k];
    }

#pragma omp simd aligned(t_17, ab_z, ds_5, fs_9 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_17[k] = ab_z[k] * ds_5[k]
                  + fs_9[k];
    }
}

}  // namespace simdtrf
