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


#include "SimdTransformDD.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_dd(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t dd,
             const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 0.5 * std::sqrt(3.0);
    const auto f_1 = std::sqrt(3.0);
    const auto f_2 = 0.25 * std::sqrt(3.0);

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
    auto *g_15 = values + 15 * nvalues;
    auto *g_16 = values + 16 * nvalues;
    auto *g_17 = values + 17 * nvalues;
    auto *g_18 = values + 18 * nvalues;
    auto *g_19 = values + 19 * nvalues;
    auto *g_20 = values + 20 * nvalues;
    auto *g_21 = values + 21 * nvalues;
    auto *g_22 = values + 22 * nvalues;
    auto *g_23 = values + 23 * nvalues;
    auto *g_24 = values + 24 * nvalues;

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_1 = buffer.data(dd + 1);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_3 = buffer.data(dd + 3);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_8 = buffer.data(dd + 8);
    const auto *dd_9 = buffer.data(dd + 9);
    const auto *dd_10 = buffer.data(dd + 10);
    const auto *dd_11 = buffer.data(dd + 11);
    const auto *dd_12 = buffer.data(dd + 12);
    const auto *dd_13 = buffer.data(dd + 13);
    const auto *dd_14 = buffer.data(dd + 14);
    const auto *dd_15 = buffer.data(dd + 15);
    const auto *dd_16 = buffer.data(dd + 16);
    const auto *dd_17 = buffer.data(dd + 17);
    const auto *dd_18 = buffer.data(dd + 18);
    const auto *dd_19 = buffer.data(dd + 19);
    const auto *dd_20 = buffer.data(dd + 20);
    const auto *dd_21 = buffer.data(dd + 21);
    const auto *dd_22 = buffer.data(dd + 22);
    const auto *dd_23 = buffer.data(dd + 23);
    const auto *dd_24 = buffer.data(dd + 24);
    const auto *dd_25 = buffer.data(dd + 25);
    const auto *dd_26 = buffer.data(dd + 26);
    const auto *dd_27 = buffer.data(dd + 27);
    const auto *dd_28 = buffer.data(dd + 28);
    const auto *dd_29 = buffer.data(dd + 29);
    const auto *dd_30 = buffer.data(dd + 30);
    const auto *dd_31 = buffer.data(dd + 31);
    const auto *dd_32 = buffer.data(dd + 32);
    const auto *dd_33 = buffer.data(dd + 33);
    const auto *dd_34 = buffer.data(dd + 34);
    const auto *dd_35 = buffer.data(dd + 35);

#pragma omp simd aligned(dd_6, dd_7, dd_8, dd_9, dd_10, dd_11, dd_25, \
                         dd_28 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = 3.0 * dd_7[k];

        g_1[k] = 3.0 * dd_10[k];

        g_2[k] = -f_0 * dd_6[k]
                 - f_0 * dd_9[k]
                 + f_1 * dd_11[k];

        g_3[k] = 3.0 * dd_8[k];

        g_4[k] = 1.5 * dd_6[k]
                 - 1.5 * dd_9[k];

        g_5[k] = 3.0 * dd_25[k];

        g_6[k] = 3.0 * dd_28[k];
    }

#pragma omp simd aligned(dd_1, dd_4, dd_19, dd_22, dd_24, dd_26, dd_27, dd_29, dd_31, \
                         dd_34 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = -f_0 * dd_24[k]
                 - f_0 * dd_27[k]
                 + f_1 * dd_29[k];

        g_8[k] = 3.0 * dd_26[k];

        g_9[k] = 1.5 * dd_24[k]
                 - 1.5 * dd_27[k];

        g_10[k] = -f_0 * dd_1[k]
                  - f_0 * dd_19[k]
                  + f_1 * dd_31[k];

        g_11[k] = -f_0 * dd_4[k]
                  - f_0 * dd_22[k]
                  + f_1 * dd_34[k];
    }

#pragma omp simd aligned(dd_0, dd_2, dd_3, dd_5, dd_18, dd_20, dd_21, dd_23, dd_30, dd_32, \
                         dd_33, dd_35 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_12[k] = 0.25 * dd_0[k]
                  + 0.25 * dd_3[k]
                  - 0.5 * dd_5[k]
                  + 0.25 * dd_18[k]
                  + 0.25 * dd_21[k]
                  - 0.5 * dd_23[k]
                  - 0.5 * dd_30[k]
                  - 0.5 * dd_33[k]
                  + dd_35[k];

        g_13[k] = -f_0 * dd_2[k]
                  - f_0 * dd_20[k]
                  + f_1 * dd_32[k];

        g_14[k] = -f_2 * dd_0[k]
                  + f_2 * dd_3[k]
                  - f_2 * dd_18[k]
                  + f_2 * dd_21[k]
                  + f_0 * dd_30[k]
                  - f_0 * dd_33[k];
    }

#pragma omp simd aligned(dd_1, dd_12, dd_13, dd_14, dd_15, dd_16, dd_17, \
                         dd_19 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_15[k] = 3.0 * dd_13[k];

        g_16[k] = 3.0 * dd_16[k];

        g_17[k] = -f_0 * dd_12[k]
                  - f_0 * dd_15[k]
                  + f_1 * dd_17[k];

        g_18[k] = 3.0 * dd_14[k];

        g_19[k] = 1.5 * dd_12[k]
                  - 1.5 * dd_15[k];

        g_20[k] = 1.5 * dd_1[k]
                  - 1.5 * dd_19[k];
    }

#pragma omp simd aligned(dd_0, dd_2, dd_3, dd_4, dd_5, dd_18, dd_20, dd_21, dd_22, \
                         dd_23 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = 1.5 * dd_4[k]
                  - 1.5 * dd_22[k];

        g_22[k] = -f_2 * dd_0[k]
                  - f_2 * dd_3[k]
                  + f_0 * dd_5[k]
                  + f_2 * dd_18[k]
                  + f_2 * dd_21[k]
                  - f_0 * dd_23[k];

        g_23[k] = 1.5 * dd_2[k]
                  - 1.5 * dd_20[k];

        g_24[k] = 0.75 * dd_0[k]
                  - 0.75 * dd_3[k]
                  - 0.75 * dd_18[k]
                  + 0.75 * dd_21[k];
    }
}

auto
transform_dd_tri(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t dd,
                 const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 0.5 * std::sqrt(3.0);
    const auto f_1 = std::sqrt(3.0);
    const auto f_2 = 0.25 * std::sqrt(3.0);

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
    auto *g_15 = values + 15 * nvalues;
    auto *g_16 = values + 16 * nvalues;
    auto *g_17 = values + 17 * nvalues;
    auto *g_18 = values + 18 * nvalues;
    auto *g_19 = values + 19 * nvalues;
    auto *g_20 = values + 20 * nvalues;
    auto *g_21 = values + 21 * nvalues;
    auto *g_22 = values + 22 * nvalues;
    auto *g_23 = values + 23 * nvalues;
    auto *g_24 = values + 24 * nvalues;

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_3 = buffer.data(dd + 3);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_8 = buffer.data(dd + 8);
    const auto *dd_9 = buffer.data(dd + 9);
    const auto *dd_10 = buffer.data(dd + 10);
    const auto *dd_11 = buffer.data(dd + 11);
    const auto *dd_12 = buffer.data(dd + 12);
    const auto *dd_14 = buffer.data(dd + 14);
    const auto *dd_15 = buffer.data(dd + 15);
    const auto *dd_18 = buffer.data(dd + 18);
    const auto *dd_20 = buffer.data(dd + 20);
    const auto *dd_21 = buffer.data(dd + 21);
    const auto *dd_23 = buffer.data(dd + 23);
    const auto *dd_24 = buffer.data(dd + 24);
    const auto *dd_26 = buffer.data(dd + 26);
    const auto *dd_27 = buffer.data(dd + 27);
    const auto *dd_28 = buffer.data(dd + 28);
    const auto *dd_29 = buffer.data(dd + 29);
    const auto *dd_30 = buffer.data(dd + 30);
    const auto *dd_32 = buffer.data(dd + 32);
    const auto *dd_33 = buffer.data(dd + 33);
    const auto *dd_35 = buffer.data(dd + 35);

#pragma omp simd aligned(dd_6, dd_7, dd_8, dd_9, dd_10, dd_11, dd_28 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = 3.0 * dd_7[k];

        g_1[k] = 3.0 * dd_10[k];
        g_5[k] = g_1[k];

        g_2[k] = -f_0 * dd_6[k]
                 - f_0 * dd_9[k]
                 + f_1 * dd_11[k];
        g_10[k] = g_2[k];

        g_3[k] = 3.0 * dd_8[k];
        g_15[k] = g_3[k];

        g_4[k] = 1.5 * dd_6[k]
                 - 1.5 * dd_9[k];
        g_20[k] = g_4[k];

        g_6[k] = 3.0 * dd_28[k];
    }

#pragma omp simd aligned(dd_24, dd_26, dd_27, dd_29 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = -f_0 * dd_24[k]
                 - f_0 * dd_27[k]
                 + f_1 * dd_29[k];
        g_11[k] = g_7[k];

        g_8[k] = 3.0 * dd_26[k];
        g_16[k] = g_8[k];

        g_9[k] = 1.5 * dd_24[k]
                 - 1.5 * dd_27[k];
        g_21[k] = g_9[k];
    }

#pragma omp simd aligned(dd_0, dd_2, dd_3, dd_5, dd_18, dd_20, dd_21, dd_23, dd_30, dd_32, \
                         dd_33, dd_35 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_12[k] = 0.25 * dd_0[k]
                  + 0.25 * dd_3[k]
                  - 0.5 * dd_5[k]
                  + 0.25 * dd_18[k]
                  + 0.25 * dd_21[k]
                  - 0.5 * dd_23[k]
                  - 0.5 * dd_30[k]
                  - 0.5 * dd_33[k]
                  + dd_35[k];

        g_13[k] = -f_0 * dd_2[k]
                  - f_0 * dd_20[k]
                  + f_1 * dd_32[k];
        g_17[k] = g_13[k];

        g_14[k] = -f_2 * dd_0[k]
                  + f_2 * dd_3[k]
                  - f_2 * dd_18[k]
                  + f_2 * dd_21[k]
                  + f_0 * dd_30[k]
                  - f_0 * dd_33[k];
        g_22[k] = g_14[k];
    }

#pragma omp simd aligned(dd_0, dd_3, dd_12, dd_14, dd_15, dd_18, \
                         dd_21 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = 3.0 * dd_14[k];

        g_19[k] = 1.5 * dd_12[k]
                  - 1.5 * dd_15[k];
        g_23[k] = g_19[k];

        g_24[k] = 0.75 * dd_0[k]
                  - 0.75 * dd_3[k]
                  - 0.75 * dd_18[k]
                  + 0.75 * dd_21[k];
    }
}

}  // namespace simdtrf
