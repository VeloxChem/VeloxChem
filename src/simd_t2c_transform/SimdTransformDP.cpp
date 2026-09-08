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


#include "SimdTransformDP.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_dp(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t dp,
             const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

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

    const auto *dp_0 = buffer.data(dp + 0);
    const auto *dp_1 = buffer.data(dp + 1);
    const auto *dp_2 = buffer.data(dp + 2);
    const auto *dp_3 = buffer.data(dp + 3);
    const auto *dp_4 = buffer.data(dp + 4);
    const auto *dp_5 = buffer.data(dp + 5);
    const auto *dp_6 = buffer.data(dp + 6);
    const auto *dp_7 = buffer.data(dp + 7);
    const auto *dp_8 = buffer.data(dp + 8);
    const auto *dp_9 = buffer.data(dp + 9);
    const auto *dp_10 = buffer.data(dp + 10);
    const auto *dp_11 = buffer.data(dp + 11);
    const auto *dp_12 = buffer.data(dp + 12);
    const auto *dp_13 = buffer.data(dp + 13);
    const auto *dp_14 = buffer.data(dp + 14);
    const auto *dp_15 = buffer.data(dp + 15);
    const auto *dp_16 = buffer.data(dp + 16);
    const auto *dp_17 = buffer.data(dp + 17);

#pragma omp simd aligned(dp_1, dp_3, dp_4, dp_5, dp_10, dp_12, dp_13, dp_14, \
                         dp_16 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * dp_4[k];

        g_1[k] = f_0 * dp_5[k];

        g_2[k] = f_0 * dp_3[k];

        g_3[k] = f_0 * dp_13[k];

        g_4[k] = f_0 * dp_14[k];

        g_5[k] = f_0 * dp_12[k];

        g_6[k] = -0.5 * dp_1[k]
                 - 0.5 * dp_10[k]
                 + dp_16[k];
    }

#pragma omp simd aligned(dp_0, dp_2, dp_6, dp_7, dp_8, dp_9, dp_11, dp_15, \
                         dp_17 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = -0.5 * dp_2[k]
                 - 0.5 * dp_11[k]
                 + dp_17[k];

        g_8[k] = -0.5 * dp_0[k]
                 - 0.5 * dp_9[k]
                 + dp_15[k];

        g_9[k] = f_0 * dp_7[k];

        g_10[k] = f_0 * dp_8[k];

        g_11[k] = f_0 * dp_6[k];
    }

#pragma omp simd aligned(dp_0, dp_1, dp_2, dp_9, dp_10, dp_11 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_12[k] = f_1 * dp_1[k]
                  - f_1 * dp_10[k];

        g_13[k] = f_1 * dp_2[k]
                  - f_1 * dp_11[k];

        g_14[k] = f_1 * dp_0[k]
                  - f_1 * dp_9[k];
    }
}

}  // namespace simdtrf
