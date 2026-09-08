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


#include "SimdTransformG.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_g_inner(CSimdMatrix &buffer, const size_t target, const size_t source,
                  const size_t nrows, const size_t ncols) -> void
{
    // NOTE: the factors are the shell's own, so they are formed once rather than
    // for every atom pair the shell pair reaches.

    const auto f_0 = 0.5 * std::sqrt(35.0);
    const auto f_1 = 0.75 * std::sqrt(70.0);
    const auto f_2 = 0.25 * std::sqrt(70.0);
    const auto f_3 = 0.5 * std::sqrt(5.0);
    const auto f_4 = 3.0 * std::sqrt(5.0);
    const auto f_5 = 0.75 * std::sqrt(10.0);
    const auto f_6 = std::sqrt(10.0);
    const auto f_7 = 0.25 * std::sqrt(5.0);
    const auto f_8 = 1.5 * std::sqrt(5.0);
    const auto f_9 = 0.125 * std::sqrt(35.0);
    const auto f_10 = 0.75 * std::sqrt(35.0);

    // NOTE: the other side of the pair reaches this pass as a count of rows and
    // nothing else, its own components running fastest within each of them.

    for (size_t r = 0; r < nrows; r++)
    {
        auto *t_0 = buffer.data(target + r * 9 + 0);
        auto *t_1 = buffer.data(target + r * 9 + 1);
        auto *t_2 = buffer.data(target + r * 9 + 2);
        auto *t_3 = buffer.data(target + r * 9 + 3);
        auto *t_4 = buffer.data(target + r * 9 + 4);
        auto *t_5 = buffer.data(target + r * 9 + 5);
        auto *t_6 = buffer.data(target + r * 9 + 6);
        auto *t_7 = buffer.data(target + r * 9 + 7);
        auto *t_8 = buffer.data(target + r * 9 + 8);

        const auto *s_0 = buffer.data(source + r * 15 + 0);
        const auto *s_1 = buffer.data(source + r * 15 + 1);
        const auto *s_2 = buffer.data(source + r * 15 + 2);
        const auto *s_3 = buffer.data(source + r * 15 + 3);
        const auto *s_4 = buffer.data(source + r * 15 + 4);
        const auto *s_5 = buffer.data(source + r * 15 + 5);
        const auto *s_6 = buffer.data(source + r * 15 + 6);
        const auto *s_7 = buffer.data(source + r * 15 + 7);
        const auto *s_8 = buffer.data(source + r * 15 + 8);
        const auto *s_9 = buffer.data(source + r * 15 + 9);
        const auto *s_10 = buffer.data(source + r * 15 + 10);
        const auto *s_11 = buffer.data(source + r * 15 + 11);
        const auto *s_12 = buffer.data(source + r * 15 + 12);
        const auto *s_13 = buffer.data(source + r * 15 + 13);
        const auto *s_14 = buffer.data(source + r * 15 + 14);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, s_1, s_4, s_6, s_8, s_11, \
                         s_13 : simd::cache_line_size())
        for (size_t k = 0; k < ncols; k++)
        {
            t_0[k] = f_0 * s_1[k]
                     - f_0 * s_6[k];

            t_1[k] = f_1 * s_4[k]
                     - f_2 * s_11[k];

            t_2[k] = -f_3 * s_1[k]
                     - f_3 * s_6[k]
                     + f_4 * s_8[k];

            t_3[k] = -f_5 * s_4[k]
                     - f_5 * s_11[k]
                     + f_6 * s_13[k];
        }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, s_0, s_2, s_3, s_5, s_7, s_9, s_10, s_12, \
                         s_14 : simd::cache_line_size())
        for (size_t k = 0; k < ncols; k++)
        {
            t_4[k] = 0.375 * s_0[k]
                     + 0.75 * s_3[k]
                     - 3.0 * s_5[k]
                     + 0.375 * s_10[k]
                     - 3.0 * s_12[k]
                     + s_14[k];

            t_5[k] = -f_5 * s_2[k]
                     - f_5 * s_7[k]
                     + f_6 * s_9[k];

            t_6[k] = -f_7 * s_0[k]
                     + f_8 * s_5[k]
                     + f_7 * s_10[k]
                     - f_8 * s_12[k];

            t_7[k] = f_2 * s_2[k]
                     - f_1 * s_7[k];

            t_8[k] = f_9 * s_0[k]
                     - f_10 * s_3[k]
                     + f_9 * s_10[k];
        }
    }
}

auto
transform_g_outer(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t source,
                  const size_t ncomps, const size_t nmax) -> void
{
    // NOTE: the factors are the shell's own, so they are formed once rather than
    // for every atom pair the shell pair reaches.

    const auto f_0 = 0.5 * std::sqrt(35.0);
    const auto f_1 = 0.75 * std::sqrt(70.0);
    const auto f_2 = 0.25 * std::sqrt(70.0);
    const auto f_3 = 0.5 * std::sqrt(5.0);
    const auto f_4 = 3.0 * std::sqrt(5.0);
    const auto f_5 = 0.75 * std::sqrt(10.0);
    const auto f_6 = std::sqrt(10.0);
    const auto f_7 = 0.25 * std::sqrt(5.0);
    const auto f_8 = 1.5 * std::sqrt(5.0);
    const auto f_9 = 0.125 * std::sqrt(35.0);
    const auto f_10 = 0.75 * std::sqrt(35.0);

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

#pragma omp simd aligned(s_1, s_4, s_6, s_8, s_11, s_13 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = f_0 * s_1[k]
                     - f_0 * s_6[k];

            g_1[k] = f_1 * s_4[k]
                     - f_2 * s_11[k];

            g_2[k] = -f_3 * s_1[k]
                     - f_3 * s_6[k]
                     + f_4 * s_8[k];

            g_3[k] = -f_5 * s_4[k]
                     - f_5 * s_11[k]
                     + f_6 * s_13[k];
        }

#pragma omp simd aligned(s_0, s_2, s_3, s_5, s_7, s_9, s_10, s_12, \
                         s_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_4[k] = 0.375 * s_0[k]
                     + 0.75 * s_3[k]
                     - 3.0 * s_5[k]
                     + 0.375 * s_10[k]
                     - 3.0 * s_12[k]
                     + s_14[k];

            g_5[k] = -f_5 * s_2[k]
                     - f_5 * s_7[k]
                     + f_6 * s_9[k];

            g_6[k] = -f_7 * s_0[k]
                     + f_8 * s_5[k]
                     + f_7 * s_10[k]
                     - f_8 * s_12[k];

            g_7[k] = f_2 * s_2[k]
                     - f_1 * s_7[k];

            g_8[k] = f_9 * s_0[k]
                     - f_10 * s_3[k]
                     + f_9 * s_10[k];
        }
    }
}

auto
transform_g_outer_tri(double *values, const size_t nvalues, CSimdMatrix &buffer,
                      const size_t source, const size_t nmax) -> void
{
    // NOTE: the factors are the shell's own, so they are formed once rather than
    // for every atom pair the shell pair reaches.

    const auto f_0 = 0.5 * std::sqrt(35.0);
    const auto f_1 = 0.75 * std::sqrt(70.0);
    const auto f_2 = 0.25 * std::sqrt(70.0);
    const auto f_3 = 0.5 * std::sqrt(5.0);
    const auto f_4 = 3.0 * std::sqrt(5.0);
    const auto f_5 = 0.75 * std::sqrt(10.0);
    const auto f_6 = std::sqrt(10.0);
    const auto f_7 = 0.25 * std::sqrt(5.0);
    const auto f_8 = 1.5 * std::sqrt(5.0);
    const auto f_9 = 0.125 * std::sqrt(35.0);
    const auto f_10 = 0.75 * std::sqrt(35.0);

    // NOTE: the rows of the values are not aligned, starting at this combination's
    // offset in the values block, so they are kept out of the clause below.

    // NOTE: the block is its own transpose, so a row above the diagonal is the
    // one below it read the other way round and is copied rather than computed.

    auto *d_0 = values + 0 * nvalues;
    auto *d_1 = values + 10 * nvalues;
    auto *d_2 = values + 20 * nvalues;
    auto *d_3 = values + 30 * nvalues;
    auto *d_4 = values + 40 * nvalues;
    auto *d_5 = values + 50 * nvalues;
    auto *d_6 = values + 60 * nvalues;
    auto *d_7 = values + 70 * nvalues;
    auto *d_8 = values + 80 * nvalues;

    const auto *q_0_1 = buffer.data(source + 9);
    const auto *q_0_6 = buffer.data(source + 54);
    const auto *q_1_4 = buffer.data(source + 37);
    const auto *q_1_11 = buffer.data(source + 100);
    const auto *q_2_1 = buffer.data(source + 11);
    const auto *q_2_6 = buffer.data(source + 56);
    const auto *q_2_8 = buffer.data(source + 74);
    const auto *q_3_4 = buffer.data(source + 39);
    const auto *q_3_11 = buffer.data(source + 102);
    const auto *q_3_13 = buffer.data(source + 120);
    const auto *q_4_0 = buffer.data(source + 4);
    const auto *q_4_3 = buffer.data(source + 31);
    const auto *q_4_5 = buffer.data(source + 49);
    const auto *q_4_10 = buffer.data(source + 94);
    const auto *q_4_12 = buffer.data(source + 112);
    const auto *q_4_14 = buffer.data(source + 130);
    const auto *q_5_2 = buffer.data(source + 23);
    const auto *q_5_7 = buffer.data(source + 68);
    const auto *q_5_9 = buffer.data(source + 86);
    const auto *q_6_0 = buffer.data(source + 6);
    const auto *q_6_5 = buffer.data(source + 51);
    const auto *q_6_10 = buffer.data(source + 96);
    const auto *q_6_12 = buffer.data(source + 114);
    const auto *q_7_2 = buffer.data(source + 25);
    const auto *q_7_7 = buffer.data(source + 70);
    const auto *q_8_0 = buffer.data(source + 8);
    const auto *q_8_3 = buffer.data(source + 35);
    const auto *q_8_10 = buffer.data(source + 98);

#pragma omp simd aligned(q_0_1, q_0_6, q_1_4, q_1_11, q_2_1, q_2_6, q_2_8, q_3_4, q_3_11, \
                         q_3_13 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        d_0[k] = f_0 * q_0_1[k]
                 - f_0 * q_0_6[k];

        d_1[k] = f_1 * q_1_4[k]
                 - f_2 * q_1_11[k];

        d_2[k] = -f_3 * q_2_1[k]
                 - f_3 * q_2_6[k]
                 + f_4 * q_2_8[k];

        d_3[k] = -f_5 * q_3_4[k]
                 - f_5 * q_3_11[k]
                 + f_6 * q_3_13[k];
    }

#pragma omp simd aligned(q_4_0, q_4_3, q_4_5, q_4_10, q_4_12, q_4_14, q_5_2, q_5_7, q_5_9, \
                         q_6_0, q_6_5, q_6_10, q_6_12 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        d_4[k] = 0.375 * q_4_0[k]
                 + 0.75 * q_4_3[k]
                 - 3.0 * q_4_5[k]
                 + 0.375 * q_4_10[k]
                 - 3.0 * q_4_12[k]
                 + q_4_14[k];

        d_5[k] = -f_5 * q_5_2[k]
                 - f_5 * q_5_7[k]
                 + f_6 * q_5_9[k];

        d_6[k] = -f_7 * q_6_0[k]
                 + f_8 * q_6_5[k]
                 + f_7 * q_6_10[k]
                 - f_8 * q_6_12[k];
    }

#pragma omp simd aligned(q_7_2, q_7_7, q_8_0, q_8_3, q_8_10 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        d_7[k] = f_2 * q_7_2[k]
                 - f_1 * q_7_7[k];

        d_8[k] = f_9 * q_8_0[k]
                 - f_10 * q_8_3[k]
                 + f_9 * q_8_10[k];
    }

    for (size_t c = 1; c < 9; c++)
    {
        auto *g_0 = values + (0 + c) * nvalues;
        auto *g_1 = values + (c * 9 + 0) * nvalues;

        const auto *s_1 = buffer.data(source + 9 + c);
        const auto *s_6 = buffer.data(source + 54 + c);

#pragma omp simd aligned(s_1, s_6 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = f_0 * s_1[k]
                     - f_0 * s_6[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 2; c < 9; c++)
    {
        auto *g_0 = values + (9 + c) * nvalues;
        auto *g_1 = values + (c * 9 + 1) * nvalues;

        const auto *s_4 = buffer.data(source + 36 + c);
        const auto *s_11 = buffer.data(source + 99 + c);

#pragma omp simd aligned(s_4, s_11 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = f_1 * s_4[k]
                     - f_2 * s_11[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 3; c < 9; c++)
    {
        auto *g_0 = values + (18 + c) * nvalues;
        auto *g_1 = values + (c * 9 + 2) * nvalues;

        const auto *s_1 = buffer.data(source + 9 + c);
        const auto *s_6 = buffer.data(source + 54 + c);
        const auto *s_8 = buffer.data(source + 72 + c);

#pragma omp simd aligned(s_1, s_6, s_8 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = -f_3 * s_1[k]
                     - f_3 * s_6[k]
                     + f_4 * s_8[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 4; c < 9; c++)
    {
        auto *g_0 = values + (27 + c) * nvalues;
        auto *g_1 = values + (c * 9 + 3) * nvalues;

        const auto *s_4 = buffer.data(source + 36 + c);
        const auto *s_11 = buffer.data(source + 99 + c);
        const auto *s_13 = buffer.data(source + 117 + c);

#pragma omp simd aligned(s_4, s_11, s_13 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = -f_5 * s_4[k]
                     - f_5 * s_11[k]
                     + f_6 * s_13[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 5; c < 9; c++)
    {
        auto *g_0 = values + (36 + c) * nvalues;
        auto *g_1 = values + (c * 9 + 4) * nvalues;

        const auto *s_0 = buffer.data(source + 0 + c);
        const auto *s_3 = buffer.data(source + 27 + c);
        const auto *s_5 = buffer.data(source + 45 + c);
        const auto *s_10 = buffer.data(source + 90 + c);
        const auto *s_12 = buffer.data(source + 108 + c);
        const auto *s_14 = buffer.data(source + 126 + c);

#pragma omp simd aligned(s_0, s_3, s_5, s_10, s_12, s_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = 0.375 * s_0[k]
                     + 0.75 * s_3[k]
                     - 3.0 * s_5[k]
                     + 0.375 * s_10[k]
                     - 3.0 * s_12[k]
                     + s_14[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 6; c < 9; c++)
    {
        auto *g_0 = values + (45 + c) * nvalues;
        auto *g_1 = values + (c * 9 + 5) * nvalues;

        const auto *s_2 = buffer.data(source + 18 + c);
        const auto *s_7 = buffer.data(source + 63 + c);
        const auto *s_9 = buffer.data(source + 81 + c);

#pragma omp simd aligned(s_2, s_7, s_9 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = -f_5 * s_2[k]
                     - f_5 * s_7[k]
                     + f_6 * s_9[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 7; c < 9; c++)
    {
        auto *g_0 = values + (54 + c) * nvalues;
        auto *g_1 = values + (c * 9 + 6) * nvalues;

        const auto *s_0 = buffer.data(source + 0 + c);
        const auto *s_5 = buffer.data(source + 45 + c);
        const auto *s_10 = buffer.data(source + 90 + c);
        const auto *s_12 = buffer.data(source + 108 + c);

#pragma omp simd aligned(s_0, s_5, s_10, s_12 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = -f_7 * s_0[k]
                     + f_8 * s_5[k]
                     + f_7 * s_10[k]
                     - f_8 * s_12[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 8; c < 9; c++)
    {
        auto *g_0 = values + (63 + c) * nvalues;
        auto *g_1 = values + (c * 9 + 7) * nvalues;

        const auto *s_2 = buffer.data(source + 18 + c);
        const auto *s_7 = buffer.data(source + 63 + c);

#pragma omp simd aligned(s_2, s_7 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = f_2 * s_2[k]
                     - f_1 * s_7[k];
            g_1[k] = g_0[k];
        }
    }
}

}  // namespace simdtrf
