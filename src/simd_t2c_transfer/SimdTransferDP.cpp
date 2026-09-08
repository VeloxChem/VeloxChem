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

#include <cmath>
#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_hrr_dp(double *values, const size_t nvalues, CSimdMatrix &buffer,
               const CSimdMatrix &coordinates, const size_t ds, const size_t fs,
               const size_t nmax) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = std::sqrt(3.0);
    const auto f_1 = 0.5 * std::sqrt(3.0);

    // NOTE: the rows of the values are not aligned, starting at this combination's
    // offset in the values block, so they are kept out of the clause below.

    auto *g_0 = values + 0 * nvalues;
    auto *g_1 = values + 1 * nvalues;
    auto *g_2 = values + 2 * nvalues;
    auto *g_3 = values + 3 * nvalues;
    auto *g_4 = values + 4 * nvalues;
    auto *g_5 = values + 5 * nvalues;
    auto *g_6 = values + 6 * nvalues;
    auto *g_7 = values + 7 * nvalues;
    auto *g_8 = values + 8 * nvalues;
    auto *g_9 = values + 9 * nvalues;
    auto *g_10 = values + 10 * nvalues;
    auto *g_11 = values + 11 * nvalues;
    auto *g_12 = values + 12 * nvalues;
    auto *g_13 = values + 13 * nvalues;
    auto *g_14 = values + 14 * nvalues;

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

#pragma omp simd aligned(ab_x, ab_y, ab_z, ds_1, ds_4, fs_1, fs_3, fs_4, fs_7, \
                         fs_8 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * ab_y[k] * ds_1[k]
                 + f_0 * fs_3[k];

        g_1[k] = f_0 * ab_z[k] * ds_1[k]
                 + f_0 * fs_4[k];

        g_2[k] = f_0 * ab_x[k] * ds_1[k]
                 + f_0 * fs_1[k];

        g_3[k] = f_0 * ab_y[k] * ds_4[k]
                 + f_0 * fs_7[k];

        g_4[k] = f_0 * ab_z[k] * ds_4[k]
                 + f_0 * fs_8[k];

        g_5[k] = f_0 * ab_x[k] * ds_4[k]
                 + f_0 * fs_4[k];
    }

#pragma omp simd aligned(ab_y, ab_z, ds_0, ds_3, ds_5, fs_1, fs_2, fs_6, fs_7, fs_8, \
                         fs_9 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = -0.5 * ab_y[k] * ds_0[k]
                 - 0.5 * ab_y[k] * ds_3[k]
                 + ab_y[k] * ds_5[k]
                 - 0.5 * fs_1[k]
                 - 0.5 * fs_6[k]
                 + fs_8[k];

        g_7[k] = -0.5 * ab_z[k] * ds_0[k]
                 - 0.5 * ab_z[k] * ds_3[k]
                 + ab_z[k] * ds_5[k]
                 - 0.5 * fs_2[k]
                 - 0.5 * fs_7[k]
                 + fs_9[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, ds_0, ds_2, ds_3, ds_5, fs_0, fs_2, fs_3, fs_4, \
                         fs_5 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_8[k] = -0.5 * ab_x[k] * ds_0[k]
                 - 0.5 * ab_x[k] * ds_3[k]
                 + ab_x[k] * ds_5[k]
                 - 0.5 * fs_0[k]
                 - 0.5 * fs_3[k]
                 + fs_5[k];

        g_9[k] = f_0 * ab_y[k] * ds_2[k]
                 + f_0 * fs_4[k];

        g_10[k] = f_0 * ab_z[k] * ds_2[k]
                  + f_0 * fs_5[k];

        g_11[k] = f_0 * ab_x[k] * ds_2[k]
                  + f_0 * fs_2[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, ds_0, ds_3, fs_0, fs_1, fs_2, fs_3, fs_6, \
                         fs_7 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_12[k] = f_1 * ab_y[k] * ds_0[k]
                  - f_1 * ab_y[k] * ds_3[k]
                  + f_1 * fs_1[k]
                  - f_1 * fs_6[k];

        g_13[k] = f_1 * ab_z[k] * ds_0[k]
                  - f_1 * ab_z[k] * ds_3[k]
                  + f_1 * fs_2[k]
                  - f_1 * fs_7[k];

        g_14[k] = f_1 * ab_x[k] * ds_0[k]
                  - f_1 * ab_x[k] * ds_3[k]
                  + f_1 * fs_0[k]
                  - f_1 * fs_3[k];
    }
}

}  // namespace simdovl
