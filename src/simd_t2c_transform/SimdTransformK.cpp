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


#include "SimdTransformK.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_k_inner(CSimdMatrix &buffer, const size_t target, const size_t source,
                  const size_t nrows, const size_t ncols) -> void
{
    // NOTE: the factors are the shell's own, so they are formed once rather than
    // for every atom pair the shell pair reaches.

    const auto f_0 = 0.21875 * std::sqrt(429.0);
    const auto f_1 = 1.09375 * std::sqrt(429.0);
    const auto f_2 = 0.65625 * std::sqrt(429.0);
    const auto f_3 = 0.03125 * std::sqrt(429.0);
    const auto f_4 = 0.1875 * std::sqrt(6006.0);
    const auto f_5 = 0.625 * std::sqrt(6006.0);
    const auto f_6 = 0.15625 * std::sqrt(231.0);
    const auto f_7 = 1.875 * std::sqrt(231.0);
    const auto f_8 = 0.28125 * std::sqrt(231.0);
    const auto f_9 = 3.75 * std::sqrt(231.0);
    const auto f_10 = 0.03125 * std::sqrt(231.0);
    const auto f_11 = 0.375 * std::sqrt(231.0);
    const auto f_12 = 0.75 * std::sqrt(231.0);
    const auto f_13 = 2.5 * std::sqrt(231.0);
    const auto f_14 = 0.28125 * std::sqrt(21.0);
    const auto f_15 = 0.46875 * std::sqrt(21.0);
    const auto f_16 = 5.625 * std::sqrt(21.0);
    const auto f_17 = 0.09375 * std::sqrt(21.0);
    const auto f_18 = 3.75 * std::sqrt(21.0);
    const auto f_19 = 7.5 * std::sqrt(21.0);
    const auto f_20 = 1.875 * std::sqrt(21.0);
    const auto f_21 = 2.5 * std::sqrt(21.0);
    const auto f_22 = 0.9375 * std::sqrt(42.0);
    const auto f_23 = 1.875 * std::sqrt(42.0);
    const auto f_24 = 5.0 * std::sqrt(42.0);
    const auto f_25 = 3.0 * std::sqrt(42.0);
    const auto f_26 = 0.15625 * std::sqrt(7.0);
    const auto f_27 = 0.46875 * std::sqrt(7.0);
    const auto f_28 = 3.75 * std::sqrt(7.0);
    const auto f_29 = 7.5 * std::sqrt(7.0);
    const auto f_30 = 2.0 * std::sqrt(7.0);
    const auto f_31 = 0.46875 * std::sqrt(42.0);
    const auto f_32 = 2.5 * std::sqrt(42.0);
    const auto f_33 = 1.5 * std::sqrt(42.0);
    const auto f_34 = 0.1875 * std::sqrt(231.0);
    const auto f_35 = 0.9375 * std::sqrt(231.0);
    const auto f_36 = 0.625 * std::sqrt(231.0);
    const auto f_37 = 0.03125 * std::sqrt(6006.0);
    const auto f_38 = 0.46875 * std::sqrt(6006.0);

    // NOTE: the other side of the pair reaches this pass as a count of rows and
    // nothing else, its own components running fastest within each of them.

    for (size_t r = 0; r < nrows; r++)
    {
        auto *t_0 = buffer.data(target + r * 15 + 0);
        auto *t_1 = buffer.data(target + r * 15 + 1);
        auto *t_2 = buffer.data(target + r * 15 + 2);
        auto *t_3 = buffer.data(target + r * 15 + 3);
        auto *t_4 = buffer.data(target + r * 15 + 4);
        auto *t_5 = buffer.data(target + r * 15 + 5);
        auto *t_6 = buffer.data(target + r * 15 + 6);
        auto *t_7 = buffer.data(target + r * 15 + 7);
        auto *t_8 = buffer.data(target + r * 15 + 8);
        auto *t_9 = buffer.data(target + r * 15 + 9);
        auto *t_10 = buffer.data(target + r * 15 + 10);
        auto *t_11 = buffer.data(target + r * 15 + 11);
        auto *t_12 = buffer.data(target + r * 15 + 12);
        auto *t_13 = buffer.data(target + r * 15 + 13);
        auto *t_14 = buffer.data(target + r * 15 + 14);

        const auto *s_0 = buffer.data(source + r * 36 + 0);
        const auto *s_1 = buffer.data(source + r * 36 + 1);
        const auto *s_2 = buffer.data(source + r * 36 + 2);
        const auto *s_3 = buffer.data(source + r * 36 + 3);
        const auto *s_4 = buffer.data(source + r * 36 + 4);
        const auto *s_5 = buffer.data(source + r * 36 + 5);
        const auto *s_6 = buffer.data(source + r * 36 + 6);
        const auto *s_7 = buffer.data(source + r * 36 + 7);
        const auto *s_8 = buffer.data(source + r * 36 + 8);
        const auto *s_9 = buffer.data(source + r * 36 + 9);
        const auto *s_10 = buffer.data(source + r * 36 + 10);
        const auto *s_11 = buffer.data(source + r * 36 + 11);
        const auto *s_12 = buffer.data(source + r * 36 + 12);
        const auto *s_13 = buffer.data(source + r * 36 + 13);
        const auto *s_14 = buffer.data(source + r * 36 + 14);
        const auto *s_15 = buffer.data(source + r * 36 + 15);
        const auto *s_16 = buffer.data(source + r * 36 + 16);
        const auto *s_17 = buffer.data(source + r * 36 + 17);
        const auto *s_18 = buffer.data(source + r * 36 + 18);
        const auto *s_19 = buffer.data(source + r * 36 + 19);
        const auto *s_20 = buffer.data(source + r * 36 + 20);
        const auto *s_21 = buffer.data(source + r * 36 + 21);
        const auto *s_22 = buffer.data(source + r * 36 + 22);
        const auto *s_23 = buffer.data(source + r * 36 + 23);
        const auto *s_24 = buffer.data(source + r * 36 + 24);
        const auto *s_25 = buffer.data(source + r * 36 + 25);
        const auto *s_26 = buffer.data(source + r * 36 + 26);
        const auto *s_27 = buffer.data(source + r * 36 + 27);
        const auto *s_28 = buffer.data(source + r * 36 + 28);
        const auto *s_29 = buffer.data(source + r * 36 + 29);
        const auto *s_30 = buffer.data(source + r * 36 + 30);
        const auto *s_31 = buffer.data(source + r * 36 + 31);
        const auto *s_32 = buffer.data(source + r * 36 + 32);
        const auto *s_33 = buffer.data(source + r * 36 + 33);
        const auto *s_34 = buffer.data(source + r * 36 + 34);
        const auto *s_35 = buffer.data(source + r * 36 + 35);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, s_1, s_4, s_6, s_8, s_11, s_13, s_15, s_17, s_22, \
                         s_24, s_28, s_30 : simd::cache_line_size())
        for (size_t k = 0; k < ncols; k++)
        {
            t_0[k] = f_0 * s_1[k]
                     - f_1 * s_6[k]
                     + f_2 * s_15[k]
                     - f_3 * s_28[k];

            t_1[k] = f_4 * s_4[k]
                     - f_5 * s_11[k]
                     + f_4 * s_22[k];

            t_2[k] = -f_6 * s_1[k]
                     + f_6 * s_6[k]
                     + f_7 * s_8[k]
                     + f_8 * s_15[k]
                     - f_9 * s_17[k]
                     - f_10 * s_28[k]
                     + f_11 * s_30[k];

            t_3[k] = -f_12 * s_4[k]
                     + f_13 * s_13[k]
                     + f_12 * s_22[k]
                     - f_13 * s_24[k];
        }

#pragma omp simd aligned(t_4, s_1, s_6, s_8, s_15, s_17, s_19, s_28, s_30, \
                         s_32 : simd::cache_line_size())
        for (size_t k = 0; k < ncols; k++)
        {
            t_4[k] = f_14 * s_1[k]
                     + f_15 * s_6[k]
                     - f_16 * s_8[k]
                     + f_17 * s_15[k]
                     - f_18 * s_17[k]
                     + f_19 * s_19[k]
                     - f_17 * s_28[k]
                     + f_20 * s_30[k]
                     - f_21 * s_32[k];
        }

#pragma omp simd aligned(t_5, s_4, s_11, s_13, s_22, s_24, s_26 : simd::cache_line_size())
        for (size_t k = 0; k < ncols; k++)
        {
            t_5[k] = f_22 * s_4[k]
                     + f_23 * s_11[k]
                     - f_24 * s_13[k]
                     + f_22 * s_22[k]
                     - f_24 * s_24[k]
                     + f_25 * s_26[k];
        }

#pragma omp simd aligned(t_6, s_1, s_6, s_8, s_15, s_17, s_19, s_28, s_30, s_32, \
                         s_34 : simd::cache_line_size())
        for (size_t k = 0; k < ncols; k++)
        {
            t_6[k] = -f_26 * s_1[k]
                     - f_27 * s_6[k]
                     + f_28 * s_8[k]
                     - f_27 * s_15[k]
                     + f_29 * s_17[k]
                     - f_29 * s_19[k]
                     - f_26 * s_28[k]
                     + f_28 * s_30[k]
                     - f_29 * s_32[k]
                     + f_30 * s_34[k];
        }

#pragma omp simd aligned(t_7, s_2, s_7, s_9, s_16, s_18, s_20, s_29, s_31, s_33, \
                         s_35 : simd::cache_line_size())
        for (size_t k = 0; k < ncols; k++)
        {
            t_7[k] = -2.1875 * s_2[k]
                     - 6.5625 * s_7[k]
                     + 13.125 * s_9[k]
                     - 6.5625 * s_16[k]
                     + 26.25 * s_18[k]
                     - 10.5 * s_20[k]
                     - 2.1875 * s_29[k]
                     + 13.125 * s_31[k]
                     - 10.5 * s_33[k]
                     + s_35[k];
        }

#pragma omp simd aligned(t_8, s_0, s_3, s_5, s_10, s_12, s_14, s_21, s_23, s_25, \
                         s_27 : simd::cache_line_size())
        for (size_t k = 0; k < ncols; k++)
        {
            t_8[k] = -f_26 * s_0[k]
                     - f_27 * s_3[k]
                     + f_28 * s_5[k]
                     - f_27 * s_10[k]
                     + f_29 * s_12[k]
                     - f_29 * s_14[k]
                     - f_26 * s_21[k]
                     + f_28 * s_23[k]
                     - f_29 * s_25[k]
                     + f_30 * s_27[k];
        }

#pragma omp simd aligned(t_9, s_2, s_7, s_9, s_16, s_20, s_29, s_31, \
                         s_33 : simd::cache_line_size())
        for (size_t k = 0; k < ncols; k++)
        {
            t_9[k] = f_31 * s_2[k]
                     + f_31 * s_7[k]
                     - f_32 * s_9[k]
                     - f_31 * s_16[k]
                     + f_33 * s_20[k]
                     - f_31 * s_29[k]
                     + f_32 * s_31[k]
                     - f_33 * s_33[k];
        }

#pragma omp simd aligned(t_10, s_0, s_3, s_5, s_10, s_12, s_14, s_21, s_23, \
                         s_25 : simd::cache_line_size())
        for (size_t k = 0; k < ncols; k++)
        {
            t_10[k] = f_17 * s_0[k]
                      - f_17 * s_3[k]
                      - f_20 * s_5[k]
                      - f_15 * s_10[k]
                      + f_18 * s_12[k]
                      + f_21 * s_14[k]
                      - f_14 * s_21[k]
                      + f_16 * s_23[k]
                      - f_19 * s_25[k];
        }

#pragma omp simd aligned(t_11, t_12, s_0, s_2, s_3, s_5, s_7, s_9, s_10, s_12, s_16, s_18, \
                         s_21, s_23, s_29, s_31 : simd::cache_line_size())
        for (size_t k = 0; k < ncols; k++)
        {
            t_11[k] = -f_34 * s_2[k]
                      + f_35 * s_7[k]
                      + f_36 * s_9[k]
                      + f_35 * s_16[k]
                      - f_9 * s_18[k]
                      - f_34 * s_29[k]
                      + f_36 * s_31[k];

            t_12[k] = -f_10 * s_0[k]
                      + f_8 * s_3[k]
                      + f_11 * s_5[k]
                      + f_6 * s_10[k]
                      - f_9 * s_12[k]
                      - f_6 * s_21[k]
                      + f_7 * s_23[k];
        }

#pragma omp simd aligned(t_13, t_14, s_0, s_2, s_3, s_7, s_10, s_16, s_21, \
                         s_29 : simd::cache_line_size())
        for (size_t k = 0; k < ncols; k++)
        {
            t_13[k] = f_37 * s_2[k]
                      - f_38 * s_7[k]
                      + f_38 * s_16[k]
                      - f_37 * s_29[k];

            t_14[k] = f_3 * s_0[k]
                      - f_2 * s_3[k]
                      + f_1 * s_10[k]
                      - f_0 * s_21[k];
        }
    }
}

auto
transform_k_outer(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t source,
                  const size_t ncomps, const size_t nmax) -> void
{
    // NOTE: the factors are the shell's own, so they are formed once rather than
    // for every atom pair the shell pair reaches.

    const auto f_0 = 0.21875 * std::sqrt(429.0);
    const auto f_1 = 1.09375 * std::sqrt(429.0);
    const auto f_2 = 0.65625 * std::sqrt(429.0);
    const auto f_3 = 0.03125 * std::sqrt(429.0);
    const auto f_4 = 0.1875 * std::sqrt(6006.0);
    const auto f_5 = 0.625 * std::sqrt(6006.0);
    const auto f_6 = 0.15625 * std::sqrt(231.0);
    const auto f_7 = 1.875 * std::sqrt(231.0);
    const auto f_8 = 0.28125 * std::sqrt(231.0);
    const auto f_9 = 3.75 * std::sqrt(231.0);
    const auto f_10 = 0.03125 * std::sqrt(231.0);
    const auto f_11 = 0.375 * std::sqrt(231.0);
    const auto f_12 = 0.75 * std::sqrt(231.0);
    const auto f_13 = 2.5 * std::sqrt(231.0);
    const auto f_14 = 0.28125 * std::sqrt(21.0);
    const auto f_15 = 0.46875 * std::sqrt(21.0);
    const auto f_16 = 5.625 * std::sqrt(21.0);
    const auto f_17 = 0.09375 * std::sqrt(21.0);
    const auto f_18 = 3.75 * std::sqrt(21.0);
    const auto f_19 = 7.5 * std::sqrt(21.0);
    const auto f_20 = 1.875 * std::sqrt(21.0);
    const auto f_21 = 2.5 * std::sqrt(21.0);
    const auto f_22 = 0.9375 * std::sqrt(42.0);
    const auto f_23 = 1.875 * std::sqrt(42.0);
    const auto f_24 = 5.0 * std::sqrt(42.0);
    const auto f_25 = 3.0 * std::sqrt(42.0);
    const auto f_26 = 0.15625 * std::sqrt(7.0);
    const auto f_27 = 0.46875 * std::sqrt(7.0);
    const auto f_28 = 3.75 * std::sqrt(7.0);
    const auto f_29 = 7.5 * std::sqrt(7.0);
    const auto f_30 = 2.0 * std::sqrt(7.0);
    const auto f_31 = 0.46875 * std::sqrt(42.0);
    const auto f_32 = 2.5 * std::sqrt(42.0);
    const auto f_33 = 1.5 * std::sqrt(42.0);
    const auto f_34 = 0.1875 * std::sqrt(231.0);
    const auto f_35 = 0.9375 * std::sqrt(231.0);
    const auto f_36 = 0.625 * std::sqrt(231.0);
    const auto f_37 = 0.03125 * std::sqrt(6006.0);
    const auto f_38 = 0.46875 * std::sqrt(6006.0);

    // NOTE: the rows of the values are not aligned, starting at this combination's
    // offset in the values block, so they are kept out of the clause below.

    // NOTE: what the other side has left reaches this pass as a count of
    // components, which is one where that side is a single function.

    for (size_t c = 0; c < ncomps; c++)
    {
        auto *g_0 = values + (0 * ncomps + c) * nvalues;
        auto *g_1 = values + (1 * ncomps + c) * nvalues;
        auto *g_2 = values + (2 * ncomps + c) * nvalues;
        auto *g_3 = values + (3 * ncomps + c) * nvalues;
        auto *g_4 = values + (4 * ncomps + c) * nvalues;
        auto *g_5 = values + (5 * ncomps + c) * nvalues;
        auto *g_6 = values + (6 * ncomps + c) * nvalues;
        auto *g_7 = values + (7 * ncomps + c) * nvalues;
        auto *g_8 = values + (8 * ncomps + c) * nvalues;
        auto *g_9 = values + (9 * ncomps + c) * nvalues;
        auto *g_10 = values + (10 * ncomps + c) * nvalues;
        auto *g_11 = values + (11 * ncomps + c) * nvalues;
        auto *g_12 = values + (12 * ncomps + c) * nvalues;
        auto *g_13 = values + (13 * ncomps + c) * nvalues;
        auto *g_14 = values + (14 * ncomps + c) * nvalues;

        const auto *s_0 = buffer.data(source + 0 * ncomps + c);
        const auto *s_1 = buffer.data(source + 1 * ncomps + c);
        const auto *s_2 = buffer.data(source + 2 * ncomps + c);
        const auto *s_3 = buffer.data(source + 3 * ncomps + c);
        const auto *s_4 = buffer.data(source + 4 * ncomps + c);
        const auto *s_5 = buffer.data(source + 5 * ncomps + c);
        const auto *s_6 = buffer.data(source + 6 * ncomps + c);
        const auto *s_7 = buffer.data(source + 7 * ncomps + c);
        const auto *s_8 = buffer.data(source + 8 * ncomps + c);
        const auto *s_9 = buffer.data(source + 9 * ncomps + c);
        const auto *s_10 = buffer.data(source + 10 * ncomps + c);
        const auto *s_11 = buffer.data(source + 11 * ncomps + c);
        const auto *s_12 = buffer.data(source + 12 * ncomps + c);
        const auto *s_13 = buffer.data(source + 13 * ncomps + c);
        const auto *s_14 = buffer.data(source + 14 * ncomps + c);
        const auto *s_15 = buffer.data(source + 15 * ncomps + c);
        const auto *s_16 = buffer.data(source + 16 * ncomps + c);
        const auto *s_17 = buffer.data(source + 17 * ncomps + c);
        const auto *s_18 = buffer.data(source + 18 * ncomps + c);
        const auto *s_19 = buffer.data(source + 19 * ncomps + c);
        const auto *s_20 = buffer.data(source + 20 * ncomps + c);
        const auto *s_21 = buffer.data(source + 21 * ncomps + c);
        const auto *s_22 = buffer.data(source + 22 * ncomps + c);
        const auto *s_23 = buffer.data(source + 23 * ncomps + c);
        const auto *s_24 = buffer.data(source + 24 * ncomps + c);
        const auto *s_25 = buffer.data(source + 25 * ncomps + c);
        const auto *s_26 = buffer.data(source + 26 * ncomps + c);
        const auto *s_27 = buffer.data(source + 27 * ncomps + c);
        const auto *s_28 = buffer.data(source + 28 * ncomps + c);
        const auto *s_29 = buffer.data(source + 29 * ncomps + c);
        const auto *s_30 = buffer.data(source + 30 * ncomps + c);
        const auto *s_31 = buffer.data(source + 31 * ncomps + c);
        const auto *s_32 = buffer.data(source + 32 * ncomps + c);
        const auto *s_33 = buffer.data(source + 33 * ncomps + c);
        const auto *s_34 = buffer.data(source + 34 * ncomps + c);
        const auto *s_35 = buffer.data(source + 35 * ncomps + c);

#pragma omp simd aligned(s_1, s_4, s_6, s_8, s_11, s_13, s_15, s_17, s_22, s_24, s_28, \
                         s_30 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = f_0 * s_1[k]
                     - f_1 * s_6[k]
                     + f_2 * s_15[k]
                     - f_3 * s_28[k];

            g_1[k] = f_4 * s_4[k]
                     - f_5 * s_11[k]
                     + f_4 * s_22[k];

            g_2[k] = -f_6 * s_1[k]
                     + f_6 * s_6[k]
                     + f_7 * s_8[k]
                     + f_8 * s_15[k]
                     - f_9 * s_17[k]
                     - f_10 * s_28[k]
                     + f_11 * s_30[k];

            g_3[k] = -f_12 * s_4[k]
                     + f_13 * s_13[k]
                     + f_12 * s_22[k]
                     - f_13 * s_24[k];
        }

#pragma omp simd aligned(s_1, s_6, s_8, s_15, s_17, s_19, s_28, s_30, \
                         s_32 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_4[k] = f_14 * s_1[k]
                     + f_15 * s_6[k]
                     - f_16 * s_8[k]
                     + f_17 * s_15[k]
                     - f_18 * s_17[k]
                     + f_19 * s_19[k]
                     - f_17 * s_28[k]
                     + f_20 * s_30[k]
                     - f_21 * s_32[k];
        }

#pragma omp simd aligned(s_4, s_11, s_13, s_22, s_24, s_26 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_5[k] = f_22 * s_4[k]
                     + f_23 * s_11[k]
                     - f_24 * s_13[k]
                     + f_22 * s_22[k]
                     - f_24 * s_24[k]
                     + f_25 * s_26[k];
        }

#pragma omp simd aligned(s_1, s_6, s_8, s_15, s_17, s_19, s_28, s_30, s_32, \
                         s_34 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_6[k] = -f_26 * s_1[k]
                     - f_27 * s_6[k]
                     + f_28 * s_8[k]
                     - f_27 * s_15[k]
                     + f_29 * s_17[k]
                     - f_29 * s_19[k]
                     - f_26 * s_28[k]
                     + f_28 * s_30[k]
                     - f_29 * s_32[k]
                     + f_30 * s_34[k];
        }

#pragma omp simd aligned(s_2, s_7, s_9, s_16, s_18, s_20, s_29, s_31, s_33, \
                         s_35 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_7[k] = -2.1875 * s_2[k]
                     - 6.5625 * s_7[k]
                     + 13.125 * s_9[k]
                     - 6.5625 * s_16[k]
                     + 26.25 * s_18[k]
                     - 10.5 * s_20[k]
                     - 2.1875 * s_29[k]
                     + 13.125 * s_31[k]
                     - 10.5 * s_33[k]
                     + s_35[k];
        }

#pragma omp simd aligned(s_0, s_3, s_5, s_10, s_12, s_14, s_21, s_23, s_25, \
                         s_27 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_8[k] = -f_26 * s_0[k]
                     - f_27 * s_3[k]
                     + f_28 * s_5[k]
                     - f_27 * s_10[k]
                     + f_29 * s_12[k]
                     - f_29 * s_14[k]
                     - f_26 * s_21[k]
                     + f_28 * s_23[k]
                     - f_29 * s_25[k]
                     + f_30 * s_27[k];
        }

#pragma omp simd aligned(s_2, s_7, s_9, s_16, s_20, s_29, s_31, s_33 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_9[k] = f_31 * s_2[k]
                     + f_31 * s_7[k]
                     - f_32 * s_9[k]
                     - f_31 * s_16[k]
                     + f_33 * s_20[k]
                     - f_31 * s_29[k]
                     + f_32 * s_31[k]
                     - f_33 * s_33[k];
        }

#pragma omp simd aligned(s_0, s_3, s_5, s_10, s_12, s_14, s_21, s_23, \
                         s_25 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_10[k] = f_17 * s_0[k]
                      - f_17 * s_3[k]
                      - f_20 * s_5[k]
                      - f_15 * s_10[k]
                      + f_18 * s_12[k]
                      + f_21 * s_14[k]
                      - f_14 * s_21[k]
                      + f_16 * s_23[k]
                      - f_19 * s_25[k];
        }

#pragma omp simd aligned(s_0, s_2, s_3, s_5, s_7, s_9, s_10, s_12, s_16, s_18, s_21, s_23, \
                         s_29, s_31 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_11[k] = -f_34 * s_2[k]
                      + f_35 * s_7[k]
                      + f_36 * s_9[k]
                      + f_35 * s_16[k]
                      - f_9 * s_18[k]
                      - f_34 * s_29[k]
                      + f_36 * s_31[k];

            g_12[k] = -f_10 * s_0[k]
                      + f_8 * s_3[k]
                      + f_11 * s_5[k]
                      + f_6 * s_10[k]
                      - f_9 * s_12[k]
                      - f_6 * s_21[k]
                      + f_7 * s_23[k];
        }

#pragma omp simd aligned(s_0, s_2, s_3, s_7, s_10, s_16, s_21, s_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_13[k] = f_37 * s_2[k]
                      - f_38 * s_7[k]
                      + f_38 * s_16[k]
                      - f_37 * s_29[k];

            g_14[k] = f_3 * s_0[k]
                      - f_2 * s_3[k]
                      + f_1 * s_10[k]
                      - f_0 * s_21[k];
        }
    }
}

auto
transform_k_outer_tri(double *values, const size_t nvalues, CSimdMatrix &buffer,
                      const size_t source, const size_t nmax) -> void
{
    // NOTE: the factors are the shell's own, so they are formed once rather than
    // for every atom pair the shell pair reaches.

    const auto f_0 = 0.21875 * std::sqrt(429.0);
    const auto f_1 = 1.09375 * std::sqrt(429.0);
    const auto f_2 = 0.65625 * std::sqrt(429.0);
    const auto f_3 = 0.03125 * std::sqrt(429.0);
    const auto f_4 = 0.1875 * std::sqrt(6006.0);
    const auto f_5 = 0.625 * std::sqrt(6006.0);
    const auto f_6 = 0.15625 * std::sqrt(231.0);
    const auto f_7 = 1.875 * std::sqrt(231.0);
    const auto f_8 = 0.28125 * std::sqrt(231.0);
    const auto f_9 = 3.75 * std::sqrt(231.0);
    const auto f_10 = 0.03125 * std::sqrt(231.0);
    const auto f_11 = 0.375 * std::sqrt(231.0);
    const auto f_12 = 0.75 * std::sqrt(231.0);
    const auto f_13 = 2.5 * std::sqrt(231.0);
    const auto f_14 = 0.28125 * std::sqrt(21.0);
    const auto f_15 = 0.46875 * std::sqrt(21.0);
    const auto f_16 = 5.625 * std::sqrt(21.0);
    const auto f_17 = 0.09375 * std::sqrt(21.0);
    const auto f_18 = 3.75 * std::sqrt(21.0);
    const auto f_19 = 7.5 * std::sqrt(21.0);
    const auto f_20 = 1.875 * std::sqrt(21.0);
    const auto f_21 = 2.5 * std::sqrt(21.0);
    const auto f_22 = 0.9375 * std::sqrt(42.0);
    const auto f_23 = 1.875 * std::sqrt(42.0);
    const auto f_24 = 5.0 * std::sqrt(42.0);
    const auto f_25 = 3.0 * std::sqrt(42.0);
    const auto f_26 = 0.15625 * std::sqrt(7.0);
    const auto f_27 = 0.46875 * std::sqrt(7.0);
    const auto f_28 = 3.75 * std::sqrt(7.0);
    const auto f_29 = 7.5 * std::sqrt(7.0);
    const auto f_30 = 2.0 * std::sqrt(7.0);
    const auto f_31 = 0.46875 * std::sqrt(42.0);
    const auto f_32 = 2.5 * std::sqrt(42.0);
    const auto f_33 = 1.5 * std::sqrt(42.0);
    const auto f_34 = 0.1875 * std::sqrt(231.0);
    const auto f_35 = 0.9375 * std::sqrt(231.0);
    const auto f_36 = 0.625 * std::sqrt(231.0);
    const auto f_37 = 0.03125 * std::sqrt(6006.0);
    const auto f_38 = 0.46875 * std::sqrt(6006.0);

    // NOTE: the rows of the values are not aligned, starting at this combination's
    // offset in the values block, so they are kept out of the clause below.

    // NOTE: the block is its own transpose, so a row above the diagonal is the
    // one below it read the other way round and is copied rather than computed.

    auto *d_0 = values + 0 * nvalues;
    auto *d_1 = values + 16 * nvalues;
    auto *d_2 = values + 32 * nvalues;
    auto *d_3 = values + 48 * nvalues;
    auto *d_4 = values + 64 * nvalues;
    auto *d_5 = values + 80 * nvalues;
    auto *d_6 = values + 96 * nvalues;
    auto *d_7 = values + 112 * nvalues;
    auto *d_8 = values + 128 * nvalues;
    auto *d_9 = values + 144 * nvalues;
    auto *d_10 = values + 160 * nvalues;
    auto *d_11 = values + 176 * nvalues;
    auto *d_12 = values + 192 * nvalues;
    auto *d_13 = values + 208 * nvalues;
    auto *d_14 = values + 224 * nvalues;

    const auto *q_0_1 = buffer.data(source + 15);
    const auto *q_0_6 = buffer.data(source + 90);
    const auto *q_0_15 = buffer.data(source + 225);
    const auto *q_0_28 = buffer.data(source + 420);
    const auto *q_1_4 = buffer.data(source + 61);
    const auto *q_1_11 = buffer.data(source + 166);
    const auto *q_1_22 = buffer.data(source + 331);
    const auto *q_2_1 = buffer.data(source + 17);
    const auto *q_2_6 = buffer.data(source + 92);
    const auto *q_2_8 = buffer.data(source + 122);
    const auto *q_2_15 = buffer.data(source + 227);
    const auto *q_2_17 = buffer.data(source + 257);
    const auto *q_2_28 = buffer.data(source + 422);
    const auto *q_2_30 = buffer.data(source + 452);
    const auto *q_3_4 = buffer.data(source + 63);
    const auto *q_3_13 = buffer.data(source + 198);
    const auto *q_3_22 = buffer.data(source + 333);
    const auto *q_3_24 = buffer.data(source + 363);
    const auto *q_4_1 = buffer.data(source + 19);
    const auto *q_4_6 = buffer.data(source + 94);
    const auto *q_4_8 = buffer.data(source + 124);
    const auto *q_4_15 = buffer.data(source + 229);
    const auto *q_4_17 = buffer.data(source + 259);
    const auto *q_4_19 = buffer.data(source + 289);
    const auto *q_4_28 = buffer.data(source + 424);
    const auto *q_4_30 = buffer.data(source + 454);
    const auto *q_4_32 = buffer.data(source + 484);
    const auto *q_5_4 = buffer.data(source + 65);
    const auto *q_5_11 = buffer.data(source + 170);
    const auto *q_5_13 = buffer.data(source + 200);
    const auto *q_5_22 = buffer.data(source + 335);
    const auto *q_5_24 = buffer.data(source + 365);
    const auto *q_5_26 = buffer.data(source + 395);
    const auto *q_6_1 = buffer.data(source + 21);
    const auto *q_6_6 = buffer.data(source + 96);
    const auto *q_6_8 = buffer.data(source + 126);
    const auto *q_6_15 = buffer.data(source + 231);
    const auto *q_6_17 = buffer.data(source + 261);
    const auto *q_6_19 = buffer.data(source + 291);
    const auto *q_6_28 = buffer.data(source + 426);
    const auto *q_6_30 = buffer.data(source + 456);
    const auto *q_6_32 = buffer.data(source + 486);
    const auto *q_6_34 = buffer.data(source + 516);
    const auto *q_7_2 = buffer.data(source + 37);
    const auto *q_7_7 = buffer.data(source + 112);
    const auto *q_7_9 = buffer.data(source + 142);
    const auto *q_7_16 = buffer.data(source + 247);
    const auto *q_7_18 = buffer.data(source + 277);
    const auto *q_7_20 = buffer.data(source + 307);
    const auto *q_7_29 = buffer.data(source + 442);
    const auto *q_7_31 = buffer.data(source + 472);
    const auto *q_7_33 = buffer.data(source + 502);
    const auto *q_7_35 = buffer.data(source + 532);
    const auto *q_8_0 = buffer.data(source + 8);
    const auto *q_8_3 = buffer.data(source + 53);
    const auto *q_8_5 = buffer.data(source + 83);
    const auto *q_8_10 = buffer.data(source + 158);
    const auto *q_8_12 = buffer.data(source + 188);
    const auto *q_8_14 = buffer.data(source + 218);
    const auto *q_8_21 = buffer.data(source + 323);
    const auto *q_8_23 = buffer.data(source + 353);
    const auto *q_8_25 = buffer.data(source + 383);
    const auto *q_8_27 = buffer.data(source + 413);
    const auto *q_9_2 = buffer.data(source + 39);
    const auto *q_9_7 = buffer.data(source + 114);
    const auto *q_9_9 = buffer.data(source + 144);
    const auto *q_9_16 = buffer.data(source + 249);
    const auto *q_9_20 = buffer.data(source + 309);
    const auto *q_9_29 = buffer.data(source + 444);
    const auto *q_9_31 = buffer.data(source + 474);
    const auto *q_9_33 = buffer.data(source + 504);
    const auto *q_10_0 = buffer.data(source + 10);
    const auto *q_10_3 = buffer.data(source + 55);
    const auto *q_10_5 = buffer.data(source + 85);
    const auto *q_10_10 = buffer.data(source + 160);
    const auto *q_10_12 = buffer.data(source + 190);
    const auto *q_10_14 = buffer.data(source + 220);
    const auto *q_10_21 = buffer.data(source + 325);
    const auto *q_10_23 = buffer.data(source + 355);
    const auto *q_10_25 = buffer.data(source + 385);
    const auto *q_11_2 = buffer.data(source + 41);
    const auto *q_11_7 = buffer.data(source + 116);
    const auto *q_11_9 = buffer.data(source + 146);
    const auto *q_11_16 = buffer.data(source + 251);
    const auto *q_11_18 = buffer.data(source + 281);
    const auto *q_11_29 = buffer.data(source + 446);
    const auto *q_11_31 = buffer.data(source + 476);
    const auto *q_12_0 = buffer.data(source + 12);
    const auto *q_12_3 = buffer.data(source + 57);
    const auto *q_12_5 = buffer.data(source + 87);
    const auto *q_12_10 = buffer.data(source + 162);
    const auto *q_12_12 = buffer.data(source + 192);
    const auto *q_12_21 = buffer.data(source + 327);
    const auto *q_12_23 = buffer.data(source + 357);
    const auto *q_13_2 = buffer.data(source + 43);
    const auto *q_13_7 = buffer.data(source + 118);
    const auto *q_13_16 = buffer.data(source + 253);
    const auto *q_13_29 = buffer.data(source + 448);
    const auto *q_14_0 = buffer.data(source + 14);
    const auto *q_14_3 = buffer.data(source + 59);
    const auto *q_14_10 = buffer.data(source + 164);
    const auto *q_14_21 = buffer.data(source + 329);

#pragma omp simd aligned(q_0_1, q_0_6, q_0_15, q_0_28, q_1_4, q_1_11, \
                         q_1_22 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        d_0[k] = f_0 * q_0_1[k]
                 - f_1 * q_0_6[k]
                 + f_2 * q_0_15[k]
                 - f_3 * q_0_28[k];

        d_1[k] = f_4 * q_1_4[k]
                 - f_5 * q_1_11[k]
                 + f_4 * q_1_22[k];
    }

#pragma omp simd aligned(q_2_1, q_2_6, q_2_8, q_2_15, q_2_17, q_2_28, q_2_30, q_3_4, q_3_13, \
                         q_3_22, q_3_24 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        d_2[k] = -f_6 * q_2_1[k]
                 + f_6 * q_2_6[k]
                 + f_7 * q_2_8[k]
                 + f_8 * q_2_15[k]
                 - f_9 * q_2_17[k]
                 - f_10 * q_2_28[k]
                 + f_11 * q_2_30[k];

        d_3[k] = -f_12 * q_3_4[k]
                 + f_13 * q_3_13[k]
                 + f_12 * q_3_22[k]
                 - f_13 * q_3_24[k];
    }

#pragma omp simd aligned(q_4_1, q_4_6, q_4_8, q_4_15, q_4_17, q_4_19, q_4_28, q_4_30, \
                         q_4_32 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        d_4[k] = f_14 * q_4_1[k]
                 + f_15 * q_4_6[k]
                 - f_16 * q_4_8[k]
                 + f_17 * q_4_15[k]
                 - f_18 * q_4_17[k]
                 + f_19 * q_4_19[k]
                 - f_17 * q_4_28[k]
                 + f_20 * q_4_30[k]
                 - f_21 * q_4_32[k];
    }

#pragma omp simd aligned(q_5_4, q_5_11, q_5_13, q_5_22, q_5_24, \
                         q_5_26 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        d_5[k] = f_22 * q_5_4[k]
                 + f_23 * q_5_11[k]
                 - f_24 * q_5_13[k]
                 + f_22 * q_5_22[k]
                 - f_24 * q_5_24[k]
                 + f_25 * q_5_26[k];
    }

#pragma omp simd aligned(q_6_1, q_6_6, q_6_8, q_6_15, q_6_17, q_6_19, q_6_28, q_6_30, q_6_32, \
                         q_6_34 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        d_6[k] = -f_26 * q_6_1[k]
                 - f_27 * q_6_6[k]
                 + f_28 * q_6_8[k]
                 - f_27 * q_6_15[k]
                 + f_29 * q_6_17[k]
                 - f_29 * q_6_19[k]
                 - f_26 * q_6_28[k]
                 + f_28 * q_6_30[k]
                 - f_29 * q_6_32[k]
                 + f_30 * q_6_34[k];
    }

#pragma omp simd aligned(q_7_2, q_7_7, q_7_9, q_7_16, q_7_18, q_7_20, q_7_29, q_7_31, q_7_33, \
                         q_7_35 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        d_7[k] = -2.1875 * q_7_2[k]
                 - 6.5625 * q_7_7[k]
                 + 13.125 * q_7_9[k]
                 - 6.5625 * q_7_16[k]
                 + 26.25 * q_7_18[k]
                 - 10.5 * q_7_20[k]
                 - 2.1875 * q_7_29[k]
                 + 13.125 * q_7_31[k]
                 - 10.5 * q_7_33[k]
                 + q_7_35[k];
    }

#pragma omp simd aligned(q_8_0, q_8_3, q_8_5, q_8_10, q_8_12, q_8_14, q_8_21, q_8_23, q_8_25, \
                         q_8_27 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        d_8[k] = -f_26 * q_8_0[k]
                 - f_27 * q_8_3[k]
                 + f_28 * q_8_5[k]
                 - f_27 * q_8_10[k]
                 + f_29 * q_8_12[k]
                 - f_29 * q_8_14[k]
                 - f_26 * q_8_21[k]
                 + f_28 * q_8_23[k]
                 - f_29 * q_8_25[k]
                 + f_30 * q_8_27[k];
    }

#pragma omp simd aligned(q_9_2, q_9_7, q_9_9, q_9_16, q_9_20, q_9_29, q_9_31, \
                         q_9_33 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        d_9[k] = f_31 * q_9_2[k]
                 + f_31 * q_9_7[k]
                 - f_32 * q_9_9[k]
                 - f_31 * q_9_16[k]
                 + f_33 * q_9_20[k]
                 - f_31 * q_9_29[k]
                 + f_32 * q_9_31[k]
                 - f_33 * q_9_33[k];
    }

#pragma omp simd aligned(q_10_0, q_10_3, q_10_5, q_10_10, q_10_12, q_10_14, q_10_21, q_10_23, \
                         q_10_25 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        d_10[k] = f_17 * q_10_0[k]
                  - f_17 * q_10_3[k]
                  - f_20 * q_10_5[k]
                  - f_15 * q_10_10[k]
                  + f_18 * q_10_12[k]
                  + f_21 * q_10_14[k]
                  - f_14 * q_10_21[k]
                  + f_16 * q_10_23[k]
                  - f_19 * q_10_25[k];
    }

#pragma omp simd aligned(q_11_2, q_11_7, q_11_9, q_11_16, q_11_18, q_11_29, q_11_31, q_12_0, \
                         q_12_3, q_12_5, q_12_10, q_12_12, q_12_21, \
                         q_12_23 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        d_11[k] = -f_34 * q_11_2[k]
                  + f_35 * q_11_7[k]
                  + f_36 * q_11_9[k]
                  + f_35 * q_11_16[k]
                  - f_9 * q_11_18[k]
                  - f_34 * q_11_29[k]
                  + f_36 * q_11_31[k];

        d_12[k] = -f_10 * q_12_0[k]
                  + f_8 * q_12_3[k]
                  + f_11 * q_12_5[k]
                  + f_6 * q_12_10[k]
                  - f_9 * q_12_12[k]
                  - f_6 * q_12_21[k]
                  + f_7 * q_12_23[k];
    }

#pragma omp simd aligned(q_13_2, q_13_7, q_13_16, q_13_29, q_14_0, q_14_3, q_14_10, \
                         q_14_21 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        d_13[k] = f_37 * q_13_2[k]
                  - f_38 * q_13_7[k]
                  + f_38 * q_13_16[k]
                  - f_37 * q_13_29[k];

        d_14[k] = f_3 * q_14_0[k]
                  - f_2 * q_14_3[k]
                  + f_1 * q_14_10[k]
                  - f_0 * q_14_21[k];
    }

    for (size_t c = 1; c < 15; c++)
    {
        auto *g_0 = values + (0 + c) * nvalues;
        auto *g_1 = values + (c * 15 + 0) * nvalues;

        const auto *s_1 = buffer.data(source + 15 + c);
        const auto *s_6 = buffer.data(source + 90 + c);
        const auto *s_15 = buffer.data(source + 225 + c);
        const auto *s_28 = buffer.data(source + 420 + c);

#pragma omp simd aligned(s_1, s_6, s_15, s_28 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = f_0 * s_1[k]
                     - f_1 * s_6[k]
                     + f_2 * s_15[k]
                     - f_3 * s_28[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 2; c < 15; c++)
    {
        auto *g_0 = values + (15 + c) * nvalues;
        auto *g_1 = values + (c * 15 + 1) * nvalues;

        const auto *s_4 = buffer.data(source + 60 + c);
        const auto *s_11 = buffer.data(source + 165 + c);
        const auto *s_22 = buffer.data(source + 330 + c);

#pragma omp simd aligned(s_4, s_11, s_22 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = f_4 * s_4[k]
                     - f_5 * s_11[k]
                     + f_4 * s_22[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 3; c < 15; c++)
    {
        auto *g_0 = values + (30 + c) * nvalues;
        auto *g_1 = values + (c * 15 + 2) * nvalues;

        const auto *s_1 = buffer.data(source + 15 + c);
        const auto *s_6 = buffer.data(source + 90 + c);
        const auto *s_8 = buffer.data(source + 120 + c);
        const auto *s_15 = buffer.data(source + 225 + c);
        const auto *s_17 = buffer.data(source + 255 + c);
        const auto *s_28 = buffer.data(source + 420 + c);
        const auto *s_30 = buffer.data(source + 450 + c);

#pragma omp simd aligned(s_1, s_6, s_8, s_15, s_17, s_28, s_30 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = -f_6 * s_1[k]
                     + f_6 * s_6[k]
                     + f_7 * s_8[k]
                     + f_8 * s_15[k]
                     - f_9 * s_17[k]
                     - f_10 * s_28[k]
                     + f_11 * s_30[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 4; c < 15; c++)
    {
        auto *g_0 = values + (45 + c) * nvalues;
        auto *g_1 = values + (c * 15 + 3) * nvalues;

        const auto *s_4 = buffer.data(source + 60 + c);
        const auto *s_13 = buffer.data(source + 195 + c);
        const auto *s_22 = buffer.data(source + 330 + c);
        const auto *s_24 = buffer.data(source + 360 + c);

#pragma omp simd aligned(s_4, s_13, s_22, s_24 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = -f_12 * s_4[k]
                     + f_13 * s_13[k]
                     + f_12 * s_22[k]
                     - f_13 * s_24[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 5; c < 15; c++)
    {
        auto *g_0 = values + (60 + c) * nvalues;
        auto *g_1 = values + (c * 15 + 4) * nvalues;

        const auto *s_1 = buffer.data(source + 15 + c);
        const auto *s_6 = buffer.data(source + 90 + c);
        const auto *s_8 = buffer.data(source + 120 + c);
        const auto *s_15 = buffer.data(source + 225 + c);
        const auto *s_17 = buffer.data(source + 255 + c);
        const auto *s_19 = buffer.data(source + 285 + c);
        const auto *s_28 = buffer.data(source + 420 + c);
        const auto *s_30 = buffer.data(source + 450 + c);
        const auto *s_32 = buffer.data(source + 480 + c);

#pragma omp simd aligned(s_1, s_6, s_8, s_15, s_17, s_19, s_28, s_30, \
                         s_32 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = f_14 * s_1[k]
                     + f_15 * s_6[k]
                     - f_16 * s_8[k]
                     + f_17 * s_15[k]
                     - f_18 * s_17[k]
                     + f_19 * s_19[k]
                     - f_17 * s_28[k]
                     + f_20 * s_30[k]
                     - f_21 * s_32[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 6; c < 15; c++)
    {
        auto *g_0 = values + (75 + c) * nvalues;
        auto *g_1 = values + (c * 15 + 5) * nvalues;

        const auto *s_4 = buffer.data(source + 60 + c);
        const auto *s_11 = buffer.data(source + 165 + c);
        const auto *s_13 = buffer.data(source + 195 + c);
        const auto *s_22 = buffer.data(source + 330 + c);
        const auto *s_24 = buffer.data(source + 360 + c);
        const auto *s_26 = buffer.data(source + 390 + c);

#pragma omp simd aligned(s_4, s_11, s_13, s_22, s_24, s_26 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = f_22 * s_4[k]
                     + f_23 * s_11[k]
                     - f_24 * s_13[k]
                     + f_22 * s_22[k]
                     - f_24 * s_24[k]
                     + f_25 * s_26[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 7; c < 15; c++)
    {
        auto *g_0 = values + (90 + c) * nvalues;
        auto *g_1 = values + (c * 15 + 6) * nvalues;

        const auto *s_1 = buffer.data(source + 15 + c);
        const auto *s_6 = buffer.data(source + 90 + c);
        const auto *s_8 = buffer.data(source + 120 + c);
        const auto *s_15 = buffer.data(source + 225 + c);
        const auto *s_17 = buffer.data(source + 255 + c);
        const auto *s_19 = buffer.data(source + 285 + c);
        const auto *s_28 = buffer.data(source + 420 + c);
        const auto *s_30 = buffer.data(source + 450 + c);
        const auto *s_32 = buffer.data(source + 480 + c);
        const auto *s_34 = buffer.data(source + 510 + c);

#pragma omp simd aligned(s_1, s_6, s_8, s_15, s_17, s_19, s_28, s_30, s_32, \
                         s_34 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = -f_26 * s_1[k]
                     - f_27 * s_6[k]
                     + f_28 * s_8[k]
                     - f_27 * s_15[k]
                     + f_29 * s_17[k]
                     - f_29 * s_19[k]
                     - f_26 * s_28[k]
                     + f_28 * s_30[k]
                     - f_29 * s_32[k]
                     + f_30 * s_34[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 8; c < 15; c++)
    {
        auto *g_0 = values + (105 + c) * nvalues;
        auto *g_1 = values + (c * 15 + 7) * nvalues;

        const auto *s_2 = buffer.data(source + 30 + c);
        const auto *s_7 = buffer.data(source + 105 + c);
        const auto *s_9 = buffer.data(source + 135 + c);
        const auto *s_16 = buffer.data(source + 240 + c);
        const auto *s_18 = buffer.data(source + 270 + c);
        const auto *s_20 = buffer.data(source + 300 + c);
        const auto *s_29 = buffer.data(source + 435 + c);
        const auto *s_31 = buffer.data(source + 465 + c);
        const auto *s_33 = buffer.data(source + 495 + c);
        const auto *s_35 = buffer.data(source + 525 + c);

#pragma omp simd aligned(s_2, s_7, s_9, s_16, s_18, s_20, s_29, s_31, s_33, \
                         s_35 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = -2.1875 * s_2[k]
                     - 6.5625 * s_7[k]
                     + 13.125 * s_9[k]
                     - 6.5625 * s_16[k]
                     + 26.25 * s_18[k]
                     - 10.5 * s_20[k]
                     - 2.1875 * s_29[k]
                     + 13.125 * s_31[k]
                     - 10.5 * s_33[k]
                     + s_35[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 9; c < 15; c++)
    {
        auto *g_0 = values + (120 + c) * nvalues;
        auto *g_1 = values + (c * 15 + 8) * nvalues;

        const auto *s_0 = buffer.data(source + 0 + c);
        const auto *s_3 = buffer.data(source + 45 + c);
        const auto *s_5 = buffer.data(source + 75 + c);
        const auto *s_10 = buffer.data(source + 150 + c);
        const auto *s_12 = buffer.data(source + 180 + c);
        const auto *s_14 = buffer.data(source + 210 + c);
        const auto *s_21 = buffer.data(source + 315 + c);
        const auto *s_23 = buffer.data(source + 345 + c);
        const auto *s_25 = buffer.data(source + 375 + c);
        const auto *s_27 = buffer.data(source + 405 + c);

#pragma omp simd aligned(s_0, s_3, s_5, s_10, s_12, s_14, s_21, s_23, s_25, \
                         s_27 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = -f_26 * s_0[k]
                     - f_27 * s_3[k]
                     + f_28 * s_5[k]
                     - f_27 * s_10[k]
                     + f_29 * s_12[k]
                     - f_29 * s_14[k]
                     - f_26 * s_21[k]
                     + f_28 * s_23[k]
                     - f_29 * s_25[k]
                     + f_30 * s_27[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 10; c < 15; c++)
    {
        auto *g_0 = values + (135 + c) * nvalues;
        auto *g_1 = values + (c * 15 + 9) * nvalues;

        const auto *s_2 = buffer.data(source + 30 + c);
        const auto *s_7 = buffer.data(source + 105 + c);
        const auto *s_9 = buffer.data(source + 135 + c);
        const auto *s_16 = buffer.data(source + 240 + c);
        const auto *s_20 = buffer.data(source + 300 + c);
        const auto *s_29 = buffer.data(source + 435 + c);
        const auto *s_31 = buffer.data(source + 465 + c);
        const auto *s_33 = buffer.data(source + 495 + c);

#pragma omp simd aligned(s_2, s_7, s_9, s_16, s_20, s_29, s_31, s_33 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = f_31 * s_2[k]
                     + f_31 * s_7[k]
                     - f_32 * s_9[k]
                     - f_31 * s_16[k]
                     + f_33 * s_20[k]
                     - f_31 * s_29[k]
                     + f_32 * s_31[k]
                     - f_33 * s_33[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 11; c < 15; c++)
    {
        auto *g_0 = values + (150 + c) * nvalues;
        auto *g_1 = values + (c * 15 + 10) * nvalues;

        const auto *s_0 = buffer.data(source + 0 + c);
        const auto *s_3 = buffer.data(source + 45 + c);
        const auto *s_5 = buffer.data(source + 75 + c);
        const auto *s_10 = buffer.data(source + 150 + c);
        const auto *s_12 = buffer.data(source + 180 + c);
        const auto *s_14 = buffer.data(source + 210 + c);
        const auto *s_21 = buffer.data(source + 315 + c);
        const auto *s_23 = buffer.data(source + 345 + c);
        const auto *s_25 = buffer.data(source + 375 + c);

#pragma omp simd aligned(s_0, s_3, s_5, s_10, s_12, s_14, s_21, s_23, \
                         s_25 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = f_17 * s_0[k]
                     - f_17 * s_3[k]
                     - f_20 * s_5[k]
                     - f_15 * s_10[k]
                     + f_18 * s_12[k]
                     + f_21 * s_14[k]
                     - f_14 * s_21[k]
                     + f_16 * s_23[k]
                     - f_19 * s_25[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 12; c < 15; c++)
    {
        auto *g_0 = values + (165 + c) * nvalues;
        auto *g_1 = values + (c * 15 + 11) * nvalues;

        const auto *s_2 = buffer.data(source + 30 + c);
        const auto *s_7 = buffer.data(source + 105 + c);
        const auto *s_9 = buffer.data(source + 135 + c);
        const auto *s_16 = buffer.data(source + 240 + c);
        const auto *s_18 = buffer.data(source + 270 + c);
        const auto *s_29 = buffer.data(source + 435 + c);
        const auto *s_31 = buffer.data(source + 465 + c);

#pragma omp simd aligned(s_2, s_7, s_9, s_16, s_18, s_29, s_31 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = -f_34 * s_2[k]
                     + f_35 * s_7[k]
                     + f_36 * s_9[k]
                     + f_35 * s_16[k]
                     - f_9 * s_18[k]
                     - f_34 * s_29[k]
                     + f_36 * s_31[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 13; c < 15; c++)
    {
        auto *g_0 = values + (180 + c) * nvalues;
        auto *g_1 = values + (c * 15 + 12) * nvalues;

        const auto *s_0 = buffer.data(source + 0 + c);
        const auto *s_3 = buffer.data(source + 45 + c);
        const auto *s_5 = buffer.data(source + 75 + c);
        const auto *s_10 = buffer.data(source + 150 + c);
        const auto *s_12 = buffer.data(source + 180 + c);
        const auto *s_21 = buffer.data(source + 315 + c);
        const auto *s_23 = buffer.data(source + 345 + c);

#pragma omp simd aligned(s_0, s_3, s_5, s_10, s_12, s_21, s_23 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = -f_10 * s_0[k]
                     + f_8 * s_3[k]
                     + f_11 * s_5[k]
                     + f_6 * s_10[k]
                     - f_9 * s_12[k]
                     - f_6 * s_21[k]
                     + f_7 * s_23[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 14; c < 15; c++)
    {
        auto *g_0 = values + (195 + c) * nvalues;
        auto *g_1 = values + (c * 15 + 13) * nvalues;

        const auto *s_2 = buffer.data(source + 30 + c);
        const auto *s_7 = buffer.data(source + 105 + c);
        const auto *s_16 = buffer.data(source + 240 + c);
        const auto *s_29 = buffer.data(source + 435 + c);

#pragma omp simd aligned(s_2, s_7, s_16, s_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = f_37 * s_2[k]
                     - f_38 * s_7[k]
                     + f_38 * s_16[k]
                     - f_37 * s_29[k];
            g_1[k] = g_0[k];
        }
    }
}

}  // namespace simdtrf
