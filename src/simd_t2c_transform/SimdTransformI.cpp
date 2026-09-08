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


#include "SimdTransformI.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_i_inner(CSimdMatrix &buffer, const size_t target, const size_t source,
                  const size_t nrows, const size_t ncols) -> void
{
    // NOTE: the factors are the shell's own, so they are formed once rather than
    // for every atom pair the shell pair reaches.

    const auto f_0 = 0.1875 * std::sqrt(462.0);
    const auto f_1 = 0.625 * std::sqrt(462.0);
    const auto f_2 = 0.9375 * std::sqrt(154.0);
    const auto f_3 = 1.875 * std::sqrt(154.0);
    const auto f_4 = 0.1875 * std::sqrt(154.0);
    const auto f_5 = 0.75 * std::sqrt(7.0);
    const auto f_6 = 7.5 * std::sqrt(7.0);
    const auto f_7 = 0.5625 * std::sqrt(210.0);
    const auto f_8 = 0.375 * std::sqrt(210.0);
    const auto f_9 = 1.5 * std::sqrt(210.0);
    const auto f_10 = 0.1875 * std::sqrt(210.0);
    const auto f_11 = 0.5 * std::sqrt(210.0);
    const auto f_12 = 0.0625 * std::sqrt(210.0);
    const auto f_13 = 0.125 * std::sqrt(210.0);
    const auto f_14 = std::sqrt(210.0);
    const auto f_15 = 0.625 * std::sqrt(21.0);
    const auto f_16 = 1.25 * std::sqrt(21.0);
    const auto f_17 = 2.5 * std::sqrt(21.0);
    const auto f_18 = std::sqrt(21.0);
    const auto f_19 = 0.03125 * std::sqrt(210.0);
    const auto f_20 = 0.1875 * std::sqrt(7.0);
    const auto f_21 = 0.9375 * std::sqrt(7.0);
    const auto f_22 = 1.875 * std::sqrt(7.0);
    const auto f_23 = 11.25 * std::sqrt(7.0);
    const auto f_24 = 0.03125 * std::sqrt(462.0);
    const auto f_25 = 0.46875 * std::sqrt(462.0);

    // NOTE: the other side of the pair reaches this pass as a count of rows and
    // nothing else, its own components running fastest within each of them.

    for (size_t r = 0; r < nrows; r++)
    {
        auto *t_0 = buffer.data(target + r * 13 + 0);
        auto *t_1 = buffer.data(target + r * 13 + 1);
        auto *t_2 = buffer.data(target + r * 13 + 2);
        auto *t_3 = buffer.data(target + r * 13 + 3);
        auto *t_4 = buffer.data(target + r * 13 + 4);
        auto *t_5 = buffer.data(target + r * 13 + 5);
        auto *t_6 = buffer.data(target + r * 13 + 6);
        auto *t_7 = buffer.data(target + r * 13 + 7);
        auto *t_8 = buffer.data(target + r * 13 + 8);
        auto *t_9 = buffer.data(target + r * 13 + 9);
        auto *t_10 = buffer.data(target + r * 13 + 10);
        auto *t_11 = buffer.data(target + r * 13 + 11);
        auto *t_12 = buffer.data(target + r * 13 + 12);

        const auto *s_0 = buffer.data(source + r * 28 + 0);
        const auto *s_1 = buffer.data(source + r * 28 + 1);
        const auto *s_2 = buffer.data(source + r * 28 + 2);
        const auto *s_3 = buffer.data(source + r * 28 + 3);
        const auto *s_4 = buffer.data(source + r * 28 + 4);
        const auto *s_5 = buffer.data(source + r * 28 + 5);
        const auto *s_6 = buffer.data(source + r * 28 + 6);
        const auto *s_7 = buffer.data(source + r * 28 + 7);
        const auto *s_8 = buffer.data(source + r * 28 + 8);
        const auto *s_9 = buffer.data(source + r * 28 + 9);
        const auto *s_10 = buffer.data(source + r * 28 + 10);
        const auto *s_11 = buffer.data(source + r * 28 + 11);
        const auto *s_12 = buffer.data(source + r * 28 + 12);
        const auto *s_13 = buffer.data(source + r * 28 + 13);
        const auto *s_14 = buffer.data(source + r * 28 + 14);
        const auto *s_15 = buffer.data(source + r * 28 + 15);
        const auto *s_16 = buffer.data(source + r * 28 + 16);
        const auto *s_17 = buffer.data(source + r * 28 + 17);
        const auto *s_18 = buffer.data(source + r * 28 + 18);
        const auto *s_19 = buffer.data(source + r * 28 + 19);
        const auto *s_20 = buffer.data(source + r * 28 + 20);
        const auto *s_21 = buffer.data(source + r * 28 + 21);
        const auto *s_22 = buffer.data(source + r * 28 + 22);
        const auto *s_23 = buffer.data(source + r * 28 + 23);
        const auto *s_24 = buffer.data(source + r * 28 + 24);
        const auto *s_25 = buffer.data(source + r * 28 + 25);
        const auto *s_26 = buffer.data(source + r * 28 + 26);
        const auto *s_27 = buffer.data(source + r * 28 + 27);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, s_1, s_4, s_6, s_8, s_11, s_13, s_15, s_17, \
                         s_19, s_22, s_24 : simd::cache_line_size())
        for (size_t k = 0; k < ncols; k++)
        {
            t_0[k] = f_0 * s_1[k]
                     - f_1 * s_6[k]
                     + f_0 * s_15[k];

            t_1[k] = f_2 * s_4[k]
                     - f_3 * s_11[k]
                     + f_4 * s_22[k];

            t_2[k] = -f_5 * s_1[k]
                     + f_6 * s_8[k]
                     + f_5 * s_15[k]
                     - f_6 * s_17[k];

            t_3[k] = -f_7 * s_4[k]
                     - f_8 * s_11[k]
                     + f_9 * s_13[k]
                     + f_10 * s_22[k]
                     - f_11 * s_24[k];

            t_4[k] = f_12 * s_1[k]
                     + f_13 * s_6[k]
                     - f_14 * s_8[k]
                     + f_12 * s_15[k]
                     - f_14 * s_17[k]
                     + f_14 * s_19[k];
        }

#pragma omp simd aligned(t_5, s_4, s_11, s_13, s_22, s_24, s_26 : simd::cache_line_size())
        for (size_t k = 0; k < ncols; k++)
        {
            t_5[k] = f_15 * s_4[k]
                     + f_16 * s_11[k]
                     - f_17 * s_13[k]
                     + f_15 * s_22[k]
                     - f_17 * s_24[k]
                     + f_18 * s_26[k];
        }

#pragma omp simd aligned(t_6, s_0, s_3, s_5, s_10, s_12, s_14, s_21, s_23, s_25, \
                         s_27 : simd::cache_line_size())
        for (size_t k = 0; k < ncols; k++)
        {
            t_6[k] = -0.3125 * s_0[k]
                     - 0.9375 * s_3[k]
                     + 5.625 * s_5[k]
                     - 0.9375 * s_10[k]
                     + 11.25 * s_12[k]
                     - 7.5 * s_14[k]
                     - 0.3125 * s_21[k]
                     + 5.625 * s_23[k]
                     - 7.5 * s_25[k]
                     + s_27[k];
        }

#pragma omp simd aligned(t_7, t_8, s_0, s_2, s_3, s_5, s_7, s_9, s_10, s_14, s_16, s_18, s_20, \
                         s_21, s_23, s_25 : simd::cache_line_size())
        for (size_t k = 0; k < ncols; k++)
        {
            t_7[k] = f_15 * s_2[k]
                     + f_16 * s_7[k]
                     - f_17 * s_9[k]
                     + f_15 * s_16[k]
                     - f_17 * s_18[k]
                     + f_18 * s_20[k];

            t_8[k] = f_19 * s_0[k]
                     + f_19 * s_3[k]
                     - f_11 * s_5[k]
                     - f_19 * s_10[k]
                     + f_11 * s_14[k]
                     - f_19 * s_21[k]
                     + f_11 * s_23[k]
                     - f_11 * s_25[k];
        }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, s_0, s_2, s_3, s_5, s_7, s_9, s_10, s_12, \
                         s_16, s_18, s_21, s_23 : simd::cache_line_size())
        for (size_t k = 0; k < ncols; k++)
        {
            t_9[k] = -f_10 * s_2[k]
                     + f_8 * s_7[k]
                     + f_11 * s_9[k]
                     + f_7 * s_16[k]
                     - f_9 * s_18[k];

            t_10[k] = -f_20 * s_0[k]
                      + f_21 * s_3[k]
                      + f_22 * s_5[k]
                      + f_21 * s_10[k]
                      - f_23 * s_12[k]
                      - f_20 * s_21[k]
                      + f_22 * s_23[k];

            t_11[k] = f_4 * s_2[k]
                      - f_3 * s_7[k]
                      + f_2 * s_16[k];

            t_12[k] = f_24 * s_0[k]
                      - f_25 * s_3[k]
                      + f_25 * s_10[k]
                      - f_24 * s_21[k];
        }
    }
}

auto
transform_i_outer(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t source,
                  const size_t ncomps, const size_t nmax) -> void
{
    // NOTE: the factors are the shell's own, so they are formed once rather than
    // for every atom pair the shell pair reaches.

    const auto f_0 = 0.1875 * std::sqrt(462.0);
    const auto f_1 = 0.625 * std::sqrt(462.0);
    const auto f_2 = 0.9375 * std::sqrt(154.0);
    const auto f_3 = 1.875 * std::sqrt(154.0);
    const auto f_4 = 0.1875 * std::sqrt(154.0);
    const auto f_5 = 0.75 * std::sqrt(7.0);
    const auto f_6 = 7.5 * std::sqrt(7.0);
    const auto f_7 = 0.5625 * std::sqrt(210.0);
    const auto f_8 = 0.375 * std::sqrt(210.0);
    const auto f_9 = 1.5 * std::sqrt(210.0);
    const auto f_10 = 0.1875 * std::sqrt(210.0);
    const auto f_11 = 0.5 * std::sqrt(210.0);
    const auto f_12 = 0.0625 * std::sqrt(210.0);
    const auto f_13 = 0.125 * std::sqrt(210.0);
    const auto f_14 = std::sqrt(210.0);
    const auto f_15 = 0.625 * std::sqrt(21.0);
    const auto f_16 = 1.25 * std::sqrt(21.0);
    const auto f_17 = 2.5 * std::sqrt(21.0);
    const auto f_18 = std::sqrt(21.0);
    const auto f_19 = 0.03125 * std::sqrt(210.0);
    const auto f_20 = 0.1875 * std::sqrt(7.0);
    const auto f_21 = 0.9375 * std::sqrt(7.0);
    const auto f_22 = 1.875 * std::sqrt(7.0);
    const auto f_23 = 11.25 * std::sqrt(7.0);
    const auto f_24 = 0.03125 * std::sqrt(462.0);
    const auto f_25 = 0.46875 * std::sqrt(462.0);

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

#pragma omp simd aligned(s_1, s_4, s_6, s_8, s_11, s_13, s_15, s_17, s_19, s_22, \
                         s_24 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = f_0 * s_1[k]
                     - f_1 * s_6[k]
                     + f_0 * s_15[k];

            g_1[k] = f_2 * s_4[k]
                     - f_3 * s_11[k]
                     + f_4 * s_22[k];

            g_2[k] = -f_5 * s_1[k]
                     + f_6 * s_8[k]
                     + f_5 * s_15[k]
                     - f_6 * s_17[k];

            g_3[k] = -f_7 * s_4[k]
                     - f_8 * s_11[k]
                     + f_9 * s_13[k]
                     + f_10 * s_22[k]
                     - f_11 * s_24[k];

            g_4[k] = f_12 * s_1[k]
                     + f_13 * s_6[k]
                     - f_14 * s_8[k]
                     + f_12 * s_15[k]
                     - f_14 * s_17[k]
                     + f_14 * s_19[k];
        }

#pragma omp simd aligned(s_4, s_11, s_13, s_22, s_24, s_26 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_5[k] = f_15 * s_4[k]
                     + f_16 * s_11[k]
                     - f_17 * s_13[k]
                     + f_15 * s_22[k]
                     - f_17 * s_24[k]
                     + f_18 * s_26[k];
        }

#pragma omp simd aligned(s_0, s_3, s_5, s_10, s_12, s_14, s_21, s_23, s_25, \
                         s_27 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_6[k] = -0.3125 * s_0[k]
                     - 0.9375 * s_3[k]
                     + 5.625 * s_5[k]
                     - 0.9375 * s_10[k]
                     + 11.25 * s_12[k]
                     - 7.5 * s_14[k]
                     - 0.3125 * s_21[k]
                     + 5.625 * s_23[k]
                     - 7.5 * s_25[k]
                     + s_27[k];
        }

#pragma omp simd aligned(s_0, s_2, s_3, s_5, s_7, s_9, s_10, s_14, s_16, s_18, s_20, s_21, \
                         s_23, s_25 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_7[k] = f_15 * s_2[k]
                     + f_16 * s_7[k]
                     - f_17 * s_9[k]
                     + f_15 * s_16[k]
                     - f_17 * s_18[k]
                     + f_18 * s_20[k];

            g_8[k] = f_19 * s_0[k]
                     + f_19 * s_3[k]
                     - f_11 * s_5[k]
                     - f_19 * s_10[k]
                     + f_11 * s_14[k]
                     - f_19 * s_21[k]
                     + f_11 * s_23[k]
                     - f_11 * s_25[k];
        }

#pragma omp simd aligned(s_0, s_2, s_3, s_5, s_7, s_9, s_10, s_12, s_16, s_18, s_21, \
                         s_23 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_9[k] = -f_10 * s_2[k]
                     + f_8 * s_7[k]
                     + f_11 * s_9[k]
                     + f_7 * s_16[k]
                     - f_9 * s_18[k];

            g_10[k] = -f_20 * s_0[k]
                      + f_21 * s_3[k]
                      + f_22 * s_5[k]
                      + f_21 * s_10[k]
                      - f_23 * s_12[k]
                      - f_20 * s_21[k]
                      + f_22 * s_23[k];

            g_11[k] = f_4 * s_2[k]
                      - f_3 * s_7[k]
                      + f_2 * s_16[k];

            g_12[k] = f_24 * s_0[k]
                      - f_25 * s_3[k]
                      + f_25 * s_10[k]
                      - f_24 * s_21[k];
        }
    }
}

auto
transform_i_outer_tri(double *values, const size_t nvalues, CSimdMatrix &buffer,
                      const size_t source, const size_t nmax) -> void
{
    // NOTE: the factors are the shell's own, so they are formed once rather than
    // for every atom pair the shell pair reaches.

    const auto f_0 = 0.1875 * std::sqrt(462.0);
    const auto f_1 = 0.625 * std::sqrt(462.0);
    const auto f_2 = 0.9375 * std::sqrt(154.0);
    const auto f_3 = 1.875 * std::sqrt(154.0);
    const auto f_4 = 0.1875 * std::sqrt(154.0);
    const auto f_5 = 0.75 * std::sqrt(7.0);
    const auto f_6 = 7.5 * std::sqrt(7.0);
    const auto f_7 = 0.5625 * std::sqrt(210.0);
    const auto f_8 = 0.375 * std::sqrt(210.0);
    const auto f_9 = 1.5 * std::sqrt(210.0);
    const auto f_10 = 0.1875 * std::sqrt(210.0);
    const auto f_11 = 0.5 * std::sqrt(210.0);
    const auto f_12 = 0.0625 * std::sqrt(210.0);
    const auto f_13 = 0.125 * std::sqrt(210.0);
    const auto f_14 = std::sqrt(210.0);
    const auto f_15 = 0.625 * std::sqrt(21.0);
    const auto f_16 = 1.25 * std::sqrt(21.0);
    const auto f_17 = 2.5 * std::sqrt(21.0);
    const auto f_18 = std::sqrt(21.0);
    const auto f_19 = 0.03125 * std::sqrt(210.0);
    const auto f_20 = 0.1875 * std::sqrt(7.0);
    const auto f_21 = 0.9375 * std::sqrt(7.0);
    const auto f_22 = 1.875 * std::sqrt(7.0);
    const auto f_23 = 11.25 * std::sqrt(7.0);
    const auto f_24 = 0.03125 * std::sqrt(462.0);
    const auto f_25 = 0.46875 * std::sqrt(462.0);

    // NOTE: the rows of the values are not aligned, starting at this combination's
    // offset in the values block, so they are kept out of the clause below.

    // NOTE: the block is its own transpose, so a row above the diagonal is the
    // one below it read the other way round and is copied rather than computed.

    auto *d_0 = values + 0 * nvalues;
    auto *d_1 = values + 14 * nvalues;
    auto *d_2 = values + 28 * nvalues;
    auto *d_3 = values + 42 * nvalues;
    auto *d_4 = values + 56 * nvalues;
    auto *d_5 = values + 70 * nvalues;
    auto *d_6 = values + 84 * nvalues;
    auto *d_7 = values + 98 * nvalues;
    auto *d_8 = values + 112 * nvalues;
    auto *d_9 = values + 126 * nvalues;
    auto *d_10 = values + 140 * nvalues;
    auto *d_11 = values + 154 * nvalues;
    auto *d_12 = values + 168 * nvalues;

    const auto *q_0_1 = buffer.data(source + 13);
    const auto *q_0_6 = buffer.data(source + 78);
    const auto *q_0_15 = buffer.data(source + 195);
    const auto *q_1_4 = buffer.data(source + 53);
    const auto *q_1_11 = buffer.data(source + 144);
    const auto *q_1_22 = buffer.data(source + 287);
    const auto *q_2_1 = buffer.data(source + 15);
    const auto *q_2_8 = buffer.data(source + 106);
    const auto *q_2_15 = buffer.data(source + 197);
    const auto *q_2_17 = buffer.data(source + 223);
    const auto *q_3_4 = buffer.data(source + 55);
    const auto *q_3_11 = buffer.data(source + 146);
    const auto *q_3_13 = buffer.data(source + 172);
    const auto *q_3_22 = buffer.data(source + 289);
    const auto *q_3_24 = buffer.data(source + 315);
    const auto *q_4_1 = buffer.data(source + 17);
    const auto *q_4_6 = buffer.data(source + 82);
    const auto *q_4_8 = buffer.data(source + 108);
    const auto *q_4_15 = buffer.data(source + 199);
    const auto *q_4_17 = buffer.data(source + 225);
    const auto *q_4_19 = buffer.data(source + 251);
    const auto *q_5_4 = buffer.data(source + 57);
    const auto *q_5_11 = buffer.data(source + 148);
    const auto *q_5_13 = buffer.data(source + 174);
    const auto *q_5_22 = buffer.data(source + 291);
    const auto *q_5_24 = buffer.data(source + 317);
    const auto *q_5_26 = buffer.data(source + 343);
    const auto *q_6_0 = buffer.data(source + 6);
    const auto *q_6_3 = buffer.data(source + 45);
    const auto *q_6_5 = buffer.data(source + 71);
    const auto *q_6_10 = buffer.data(source + 136);
    const auto *q_6_12 = buffer.data(source + 162);
    const auto *q_6_14 = buffer.data(source + 188);
    const auto *q_6_21 = buffer.data(source + 279);
    const auto *q_6_23 = buffer.data(source + 305);
    const auto *q_6_25 = buffer.data(source + 331);
    const auto *q_6_27 = buffer.data(source + 357);
    const auto *q_7_2 = buffer.data(source + 33);
    const auto *q_7_7 = buffer.data(source + 98);
    const auto *q_7_9 = buffer.data(source + 124);
    const auto *q_7_16 = buffer.data(source + 215);
    const auto *q_7_18 = buffer.data(source + 241);
    const auto *q_7_20 = buffer.data(source + 267);
    const auto *q_8_0 = buffer.data(source + 8);
    const auto *q_8_3 = buffer.data(source + 47);
    const auto *q_8_5 = buffer.data(source + 73);
    const auto *q_8_10 = buffer.data(source + 138);
    const auto *q_8_14 = buffer.data(source + 190);
    const auto *q_8_21 = buffer.data(source + 281);
    const auto *q_8_23 = buffer.data(source + 307);
    const auto *q_8_25 = buffer.data(source + 333);
    const auto *q_9_2 = buffer.data(source + 35);
    const auto *q_9_7 = buffer.data(source + 100);
    const auto *q_9_9 = buffer.data(source + 126);
    const auto *q_9_16 = buffer.data(source + 217);
    const auto *q_9_18 = buffer.data(source + 243);
    const auto *q_10_0 = buffer.data(source + 10);
    const auto *q_10_3 = buffer.data(source + 49);
    const auto *q_10_5 = buffer.data(source + 75);
    const auto *q_10_10 = buffer.data(source + 140);
    const auto *q_10_12 = buffer.data(source + 166);
    const auto *q_10_21 = buffer.data(source + 283);
    const auto *q_10_23 = buffer.data(source + 309);
    const auto *q_11_2 = buffer.data(source + 37);
    const auto *q_11_7 = buffer.data(source + 102);
    const auto *q_11_16 = buffer.data(source + 219);
    const auto *q_12_0 = buffer.data(source + 12);
    const auto *q_12_3 = buffer.data(source + 51);
    const auto *q_12_10 = buffer.data(source + 142);
    const auto *q_12_21 = buffer.data(source + 285);

#pragma omp simd aligned(q_0_1, q_0_6, q_0_15, q_1_4, q_1_11, q_1_22, q_2_1, q_2_8, q_2_15, \
                         q_2_17 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        d_0[k] = f_0 * q_0_1[k]
                 - f_1 * q_0_6[k]
                 + f_0 * q_0_15[k];

        d_1[k] = f_2 * q_1_4[k]
                 - f_3 * q_1_11[k]
                 + f_4 * q_1_22[k];

        d_2[k] = -f_5 * q_2_1[k]
                 + f_6 * q_2_8[k]
                 + f_5 * q_2_15[k]
                 - f_6 * q_2_17[k];
    }

#pragma omp simd aligned(q_3_4, q_3_11, q_3_13, q_3_22, q_3_24, q_4_1, q_4_6, q_4_8, q_4_15, \
                         q_4_17, q_4_19 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        d_3[k] = -f_7 * q_3_4[k]
                 - f_8 * q_3_11[k]
                 + f_9 * q_3_13[k]
                 + f_10 * q_3_22[k]
                 - f_11 * q_3_24[k];

        d_4[k] = f_12 * q_4_1[k]
                 + f_13 * q_4_6[k]
                 - f_14 * q_4_8[k]
                 + f_12 * q_4_15[k]
                 - f_14 * q_4_17[k]
                 + f_14 * q_4_19[k];
    }

#pragma omp simd aligned(q_5_4, q_5_11, q_5_13, q_5_22, q_5_24, \
                         q_5_26 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        d_5[k] = f_15 * q_5_4[k]
                 + f_16 * q_5_11[k]
                 - f_17 * q_5_13[k]
                 + f_15 * q_5_22[k]
                 - f_17 * q_5_24[k]
                 + f_18 * q_5_26[k];
    }

#pragma omp simd aligned(q_6_0, q_6_3, q_6_5, q_6_10, q_6_12, q_6_14, q_6_21, q_6_23, q_6_25, \
                         q_6_27 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        d_6[k] = -0.3125 * q_6_0[k]
                 - 0.9375 * q_6_3[k]
                 + 5.625 * q_6_5[k]
                 - 0.9375 * q_6_10[k]
                 + 11.25 * q_6_12[k]
                 - 7.5 * q_6_14[k]
                 - 0.3125 * q_6_21[k]
                 + 5.625 * q_6_23[k]
                 - 7.5 * q_6_25[k]
                 + q_6_27[k];
    }

#pragma omp simd aligned(q_7_2, q_7_7, q_7_9, q_7_16, q_7_18, q_7_20, q_8_0, q_8_3, q_8_5, \
                         q_8_10, q_8_14, q_8_21, q_8_23, q_8_25 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        d_7[k] = f_15 * q_7_2[k]
                 + f_16 * q_7_7[k]
                 - f_17 * q_7_9[k]
                 + f_15 * q_7_16[k]
                 - f_17 * q_7_18[k]
                 + f_18 * q_7_20[k];

        d_8[k] = f_19 * q_8_0[k]
                 + f_19 * q_8_3[k]
                 - f_11 * q_8_5[k]
                 - f_19 * q_8_10[k]
                 + f_11 * q_8_14[k]
                 - f_19 * q_8_21[k]
                 + f_11 * q_8_23[k]
                 - f_11 * q_8_25[k];
    }

#pragma omp simd aligned(q_9_2, q_9_7, q_9_9, q_9_16, q_9_18, q_10_0, q_10_3, q_10_5, q_10_10, \
                         q_10_12, q_10_21, q_10_23 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        d_9[k] = -f_10 * q_9_2[k]
                 + f_8 * q_9_7[k]
                 + f_11 * q_9_9[k]
                 + f_7 * q_9_16[k]
                 - f_9 * q_9_18[k];

        d_10[k] = -f_20 * q_10_0[k]
                  + f_21 * q_10_3[k]
                  + f_22 * q_10_5[k]
                  + f_21 * q_10_10[k]
                  - f_23 * q_10_12[k]
                  - f_20 * q_10_21[k]
                  + f_22 * q_10_23[k];
    }

#pragma omp simd aligned(q_11_2, q_11_7, q_11_16, q_12_0, q_12_3, q_12_10, \
                         q_12_21 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        d_11[k] = f_4 * q_11_2[k]
                  - f_3 * q_11_7[k]
                  + f_2 * q_11_16[k];

        d_12[k] = f_24 * q_12_0[k]
                  - f_25 * q_12_3[k]
                  + f_25 * q_12_10[k]
                  - f_24 * q_12_21[k];
    }

    for (size_t c = 1; c < 13; c++)
    {
        auto *g_0 = values + (0 + c) * nvalues;
        auto *g_1 = values + (c * 13 + 0) * nvalues;

        const auto *s_1 = buffer.data(source + 13 + c);
        const auto *s_6 = buffer.data(source + 78 + c);
        const auto *s_15 = buffer.data(source + 195 + c);

#pragma omp simd aligned(s_1, s_6, s_15 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = f_0 * s_1[k]
                     - f_1 * s_6[k]
                     + f_0 * s_15[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 2; c < 13; c++)
    {
        auto *g_0 = values + (13 + c) * nvalues;
        auto *g_1 = values + (c * 13 + 1) * nvalues;

        const auto *s_4 = buffer.data(source + 52 + c);
        const auto *s_11 = buffer.data(source + 143 + c);
        const auto *s_22 = buffer.data(source + 286 + c);

#pragma omp simd aligned(s_4, s_11, s_22 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = f_2 * s_4[k]
                     - f_3 * s_11[k]
                     + f_4 * s_22[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 3; c < 13; c++)
    {
        auto *g_0 = values + (26 + c) * nvalues;
        auto *g_1 = values + (c * 13 + 2) * nvalues;

        const auto *s_1 = buffer.data(source + 13 + c);
        const auto *s_8 = buffer.data(source + 104 + c);
        const auto *s_15 = buffer.data(source + 195 + c);
        const auto *s_17 = buffer.data(source + 221 + c);

#pragma omp simd aligned(s_1, s_8, s_15, s_17 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = -f_5 * s_1[k]
                     + f_6 * s_8[k]
                     + f_5 * s_15[k]
                     - f_6 * s_17[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 4; c < 13; c++)
    {
        auto *g_0 = values + (39 + c) * nvalues;
        auto *g_1 = values + (c * 13 + 3) * nvalues;

        const auto *s_4 = buffer.data(source + 52 + c);
        const auto *s_11 = buffer.data(source + 143 + c);
        const auto *s_13 = buffer.data(source + 169 + c);
        const auto *s_22 = buffer.data(source + 286 + c);
        const auto *s_24 = buffer.data(source + 312 + c);

#pragma omp simd aligned(s_4, s_11, s_13, s_22, s_24 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = -f_7 * s_4[k]
                     - f_8 * s_11[k]
                     + f_9 * s_13[k]
                     + f_10 * s_22[k]
                     - f_11 * s_24[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 5; c < 13; c++)
    {
        auto *g_0 = values + (52 + c) * nvalues;
        auto *g_1 = values + (c * 13 + 4) * nvalues;

        const auto *s_1 = buffer.data(source + 13 + c);
        const auto *s_6 = buffer.data(source + 78 + c);
        const auto *s_8 = buffer.data(source + 104 + c);
        const auto *s_15 = buffer.data(source + 195 + c);
        const auto *s_17 = buffer.data(source + 221 + c);
        const auto *s_19 = buffer.data(source + 247 + c);

#pragma omp simd aligned(s_1, s_6, s_8, s_15, s_17, s_19 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = f_12 * s_1[k]
                     + f_13 * s_6[k]
                     - f_14 * s_8[k]
                     + f_12 * s_15[k]
                     - f_14 * s_17[k]
                     + f_14 * s_19[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 6; c < 13; c++)
    {
        auto *g_0 = values + (65 + c) * nvalues;
        auto *g_1 = values + (c * 13 + 5) * nvalues;

        const auto *s_4 = buffer.data(source + 52 + c);
        const auto *s_11 = buffer.data(source + 143 + c);
        const auto *s_13 = buffer.data(source + 169 + c);
        const auto *s_22 = buffer.data(source + 286 + c);
        const auto *s_24 = buffer.data(source + 312 + c);
        const auto *s_26 = buffer.data(source + 338 + c);

#pragma omp simd aligned(s_4, s_11, s_13, s_22, s_24, s_26 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = f_15 * s_4[k]
                     + f_16 * s_11[k]
                     - f_17 * s_13[k]
                     + f_15 * s_22[k]
                     - f_17 * s_24[k]
                     + f_18 * s_26[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 7; c < 13; c++)
    {
        auto *g_0 = values + (78 + c) * nvalues;
        auto *g_1 = values + (c * 13 + 6) * nvalues;

        const auto *s_0 = buffer.data(source + 0 + c);
        const auto *s_3 = buffer.data(source + 39 + c);
        const auto *s_5 = buffer.data(source + 65 + c);
        const auto *s_10 = buffer.data(source + 130 + c);
        const auto *s_12 = buffer.data(source + 156 + c);
        const auto *s_14 = buffer.data(source + 182 + c);
        const auto *s_21 = buffer.data(source + 273 + c);
        const auto *s_23 = buffer.data(source + 299 + c);
        const auto *s_25 = buffer.data(source + 325 + c);
        const auto *s_27 = buffer.data(source + 351 + c);

#pragma omp simd aligned(s_0, s_3, s_5, s_10, s_12, s_14, s_21, s_23, s_25, \
                         s_27 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = -0.3125 * s_0[k]
                     - 0.9375 * s_3[k]
                     + 5.625 * s_5[k]
                     - 0.9375 * s_10[k]
                     + 11.25 * s_12[k]
                     - 7.5 * s_14[k]
                     - 0.3125 * s_21[k]
                     + 5.625 * s_23[k]
                     - 7.5 * s_25[k]
                     + s_27[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 8; c < 13; c++)
    {
        auto *g_0 = values + (91 + c) * nvalues;
        auto *g_1 = values + (c * 13 + 7) * nvalues;

        const auto *s_2 = buffer.data(source + 26 + c);
        const auto *s_7 = buffer.data(source + 91 + c);
        const auto *s_9 = buffer.data(source + 117 + c);
        const auto *s_16 = buffer.data(source + 208 + c);
        const auto *s_18 = buffer.data(source + 234 + c);
        const auto *s_20 = buffer.data(source + 260 + c);

#pragma omp simd aligned(s_2, s_7, s_9, s_16, s_18, s_20 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = f_15 * s_2[k]
                     + f_16 * s_7[k]
                     - f_17 * s_9[k]
                     + f_15 * s_16[k]
                     - f_17 * s_18[k]
                     + f_18 * s_20[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 9; c < 13; c++)
    {
        auto *g_0 = values + (104 + c) * nvalues;
        auto *g_1 = values + (c * 13 + 8) * nvalues;

        const auto *s_0 = buffer.data(source + 0 + c);
        const auto *s_3 = buffer.data(source + 39 + c);
        const auto *s_5 = buffer.data(source + 65 + c);
        const auto *s_10 = buffer.data(source + 130 + c);
        const auto *s_14 = buffer.data(source + 182 + c);
        const auto *s_21 = buffer.data(source + 273 + c);
        const auto *s_23 = buffer.data(source + 299 + c);
        const auto *s_25 = buffer.data(source + 325 + c);

#pragma omp simd aligned(s_0, s_3, s_5, s_10, s_14, s_21, s_23, s_25 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = f_19 * s_0[k]
                     + f_19 * s_3[k]
                     - f_11 * s_5[k]
                     - f_19 * s_10[k]
                     + f_11 * s_14[k]
                     - f_19 * s_21[k]
                     + f_11 * s_23[k]
                     - f_11 * s_25[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 10; c < 13; c++)
    {
        auto *g_0 = values + (117 + c) * nvalues;
        auto *g_1 = values + (c * 13 + 9) * nvalues;

        const auto *s_2 = buffer.data(source + 26 + c);
        const auto *s_7 = buffer.data(source + 91 + c);
        const auto *s_9 = buffer.data(source + 117 + c);
        const auto *s_16 = buffer.data(source + 208 + c);
        const auto *s_18 = buffer.data(source + 234 + c);

#pragma omp simd aligned(s_2, s_7, s_9, s_16, s_18 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = -f_10 * s_2[k]
                     + f_8 * s_7[k]
                     + f_11 * s_9[k]
                     + f_7 * s_16[k]
                     - f_9 * s_18[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 11; c < 13; c++)
    {
        auto *g_0 = values + (130 + c) * nvalues;
        auto *g_1 = values + (c * 13 + 10) * nvalues;

        const auto *s_0 = buffer.data(source + 0 + c);
        const auto *s_3 = buffer.data(source + 39 + c);
        const auto *s_5 = buffer.data(source + 65 + c);
        const auto *s_10 = buffer.data(source + 130 + c);
        const auto *s_12 = buffer.data(source + 156 + c);
        const auto *s_21 = buffer.data(source + 273 + c);
        const auto *s_23 = buffer.data(source + 299 + c);

#pragma omp simd aligned(s_0, s_3, s_5, s_10, s_12, s_21, s_23 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = -f_20 * s_0[k]
                     + f_21 * s_3[k]
                     + f_22 * s_5[k]
                     + f_21 * s_10[k]
                     - f_23 * s_12[k]
                     - f_20 * s_21[k]
                     + f_22 * s_23[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 12; c < 13; c++)
    {
        auto *g_0 = values + (143 + c) * nvalues;
        auto *g_1 = values + (c * 13 + 11) * nvalues;

        const auto *s_2 = buffer.data(source + 26 + c);
        const auto *s_7 = buffer.data(source + 91 + c);
        const auto *s_16 = buffer.data(source + 208 + c);

#pragma omp simd aligned(s_2, s_7, s_16 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = f_4 * s_2[k]
                     - f_3 * s_7[k]
                     + f_2 * s_16[k];
            g_1[k] = g_0[k];
        }
    }
}

}  // namespace simdtrf
