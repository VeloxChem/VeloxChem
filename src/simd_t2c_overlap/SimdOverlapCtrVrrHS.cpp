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


#include "SimdOverlapCtrVrrHS.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_ctr_hs_overlap_0(double *values, const size_t nvalues, CSimdMatrix &buffer,
                         const size_t pa, const size_t fs, const size_t gs, const size_t ncols,
                         const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5625 * std::sqrt(14.0) / p;
    const auto f_1 = 0.9375 * std::sqrt(14.0);
    const auto f_2 = 1.875 * std::sqrt(14.0);
    const auto f_3 = 0.1875 * std::sqrt(14.0);
    const auto f_4 = 1.5 * std::sqrt(35.0);
    const auto f_5 = 0.0625 * std::sqrt(70.0) / p;
    const auto f_6 = 0.5 * std::sqrt(70.0) / p;
    const auto f_7 = 0.1875 * std::sqrt(70.0);
    const auto f_8 = 1.5 * std::sqrt(70.0);
    const auto f_9 = 0.125 * std::sqrt(70.0);
    const auto f_10 = 0.0625 * std::sqrt(70.0);
    const auto f_11 = 0.5 * std::sqrt(70.0);
    const auto f_12 = 0.5 * std::sqrt(105.0);
    const auto f_13 = std::sqrt(105.0);
    const auto f_14 = 0.375 * std::sqrt(15.0) / p;
    const auto f_15 = 1.5 * std::sqrt(15.0) / p;
    const auto f_16 = 0.125 * std::sqrt(15.0);
    const auto f_17 = 1.5 * std::sqrt(15.0);
    const auto f_18 = 0.25 * std::sqrt(15.0);
    const auto f_19 = std::sqrt(15.0);
    const auto f_20 = 3.0 / p;
    const auto f_21 = 0.25 * std::sqrt(15.0) / p;
    const auto f_22 = 0.25 * std::sqrt(105.0);
    const auto f_23 = 0.125 * std::sqrt(70.0) / p;
    const auto f_24 = 0.375 * std::sqrt(35.0);
    const auto f_25 = 2.25 * std::sqrt(35.0);
    const auto f_26 = 0.375 * std::sqrt(14.0) / p;
    const auto f_27 = 1.875 * std::sqrt(14.0) / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_1 = buffer.data(fs + 1);
    const auto *fs_2 = buffer.data(fs + 2);
    const auto *fs_3 = buffer.data(fs + 3);
    const auto *fs_4 = buffer.data(fs + 4);
    const auto *fs_5 = buffer.data(fs + 5);

    const auto *gs_0 = buffer.data(gs + 0);
    const auto *gs_1 = buffer.data(gs + 1);
    const auto *gs_2 = buffer.data(gs + 2);
    const auto *gs_3 = buffer.data(gs + 3);
    const auto *gs_4 = buffer.data(gs + 4);
    const auto *gs_5 = buffer.data(gs + 5);
    const auto *gs_6 = buffer.data(gs + 6);
    const auto *gs_7 = buffer.data(gs + 7);
    const auto *gs_8 = buffer.data(gs + 8);
    const auto *gs_9 = buffer.data(gs + 9);
    const auto *gs_10 = buffer.data(gs + 10);

#pragma omp simd aligned(pa_x, pa_y, fs_3, fs_4, gs_0, gs_1, gs_3, gs_4, gs_6, gs_7, gs_8, \
                         gs_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_0[k] += -f_0 * fs_3[k]
                  + f_1 * pa_y[k] * gs_0[k]
                  - f_2 * pa_x[k] * gs_4[k]
                  + f_3 * pa_y[k] * gs_6[k];

        g_1[k] += f_4 * pa_y[k] * gs_1[k]
                  - f_4 * pa_x[k] * gs_7[k];

        g_2[k] += f_5 * fs_3[k]
                  - f_6 * fs_4[k]
                  - f_7 * pa_y[k] * gs_0[k]
                  + f_8 * pa_y[k] * gs_3[k]
                  - f_9 * pa_x[k] * gs_4[k]
                  + f_10 * pa_y[k] * gs_6[k]
                  - f_11 * pa_y[k] * gs_8[k];

        g_3[k] += -f_12 * pa_y[k] * gs_1[k]
                  - f_12 * pa_x[k] * gs_7[k]
                  + f_13 * pa_x[k] * gs_9[k];
    }

#pragma omp simd aligned(pa_x, pa_y, fs_3, fs_4, gs_0, gs_3, gs_4, gs_6, gs_8, \
                         gs_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_4[k] += f_14 * fs_3[k]
                  - f_15 * fs_4[k]
                  + f_16 * pa_y[k] * gs_0[k]
                  - f_17 * pa_y[k] * gs_3[k]
                  + f_18 * pa_x[k] * gs_4[k]
                  + f_16 * pa_y[k] * gs_6[k]
                  - f_17 * pa_y[k] * gs_8[k]
                  + f_19 * pa_y[k] * gs_10[k];
    }

#pragma omp simd aligned(pa_x, pa_y, pa_z, fs_5, gs_0, gs_2, gs_5, gs_6, gs_9, \
                         gs_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_5[k] += -f_20 * fs_5[k]
                  + 1.875 * pa_z[k] * gs_0[k]
                  + 3.75 * pa_z[k] * gs_2[k]
                  - 5.0 * pa_x[k] * gs_5[k]
                  + 1.875 * pa_z[k] * gs_6[k]
                  - 5.0 * pa_y[k] * gs_9[k]
                  + pa_z[k] * gs_10[k];
    }

#pragma omp simd aligned(pa_x, pa_y, pa_z, fs_0, fs_1, fs_2, gs_0, gs_2, gs_3, gs_5, gs_6, \
                         gs_8, gs_9, gs_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_6[k] += f_21 * fs_0[k]
                  + f_21 * fs_1[k]
                  - f_15 * fs_2[k]
                  + f_16 * pa_x[k] * gs_0[k]
                  + f_18 * pa_x[k] * gs_2[k]
                  - f_17 * pa_x[k] * gs_3[k]
                  + f_16 * pa_x[k] * gs_6[k]
                  - f_17 * pa_x[k] * gs_8[k]
                  + f_19 * pa_x[k] * gs_10[k];

        g_7[k] += -f_22 * pa_z[k] * gs_0[k]
                  + f_12 * pa_x[k] * gs_5[k]
                  + f_22 * pa_z[k] * gs_6[k]
                  - f_12 * pa_y[k] * gs_9[k];
    }

#pragma omp simd aligned(pa_x, pa_z, fs_0, fs_1, fs_2, gs_0, gs_2, gs_3, gs_6, \
                         gs_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_8[k] += -f_23 * fs_0[k]
                  + f_23 * fs_1[k]
                  + f_6 * fs_2[k]
                  - f_10 * pa_x[k] * gs_0[k]
                  + f_9 * pa_x[k] * gs_2[k]
                  + f_11 * pa_x[k] * gs_3[k]
                  + f_7 * pa_x[k] * gs_6[k]
                  - f_8 * pa_x[k] * gs_8[k];

        g_9[k] += f_24 * pa_z[k] * gs_0[k]
                  - f_25 * pa_z[k] * gs_2[k]
                  + f_24 * pa_z[k] * gs_6[k];

        g_10[k] += f_26 * fs_0[k]
                   - f_27 * fs_1[k]
                   + f_3 * pa_x[k] * gs_0[k]
                   - f_2 * pa_x[k] * gs_2[k]
                   + f_1 * pa_x[k] * gs_6[k];
    }
}

}  // namespace simdovl
