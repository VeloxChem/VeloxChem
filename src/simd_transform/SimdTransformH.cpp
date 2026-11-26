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


#include "SimdTransformH.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_h_inner(CSimdMatrix &buffer, const size_t target, const size_t source,
                  const size_t nrows, const size_t ncomps, const size_t ncols) -> void
{
    // NOTE: the factors are the shell's own, so they are formed once rather than
    // for every atom pair the shell pair reaches.

    const auto f_0 = 0.9375 * std::sqrt(14.0);
    const auto f_1 = 1.875 * std::sqrt(14.0);
    const auto f_2 = 0.1875 * std::sqrt(14.0);
    const auto f_3 = 1.5 * std::sqrt(35.0);
    const auto f_4 = 0.1875 * std::sqrt(70.0);
    const auto f_5 = 0.125 * std::sqrt(70.0);
    const auto f_6 = 1.5 * std::sqrt(70.0);
    const auto f_7 = 0.0625 * std::sqrt(70.0);
    const auto f_8 = 0.5 * std::sqrt(70.0);
    const auto f_9 = 0.5 * std::sqrt(105.0);
    const auto f_10 = std::sqrt(105.0);
    const auto f_11 = 0.125 * std::sqrt(15.0);
    const auto f_12 = 0.25 * std::sqrt(15.0);
    const auto f_13 = 1.5 * std::sqrt(15.0);
    const auto f_14 = std::sqrt(15.0);
    const auto f_15 = 0.25 * std::sqrt(105.0);
    const auto f_16 = 0.375 * std::sqrt(35.0);
    const auto f_17 = 2.25 * std::sqrt(35.0);

    // NOTE: what sits either side of this index reaches the pass as a count and
    // nothing else -- the rows above it, the components below -- so one routine
    // per shell serves every block the shell appears in.

    for (size_t r = 0; r < nrows; r++)
    {
        for (size_t c = 0; c < ncomps; c++)
        {
            auto *t_0 = buffer.data(target + (r * 11 + 0) * ncomps + c);
            auto *t_1 = buffer.data(target + (r * 11 + 1) * ncomps + c);
            auto *t_2 = buffer.data(target + (r * 11 + 2) * ncomps + c);
            auto *t_3 = buffer.data(target + (r * 11 + 3) * ncomps + c);
            auto *t_4 = buffer.data(target + (r * 11 + 4) * ncomps + c);
            auto *t_5 = buffer.data(target + (r * 11 + 5) * ncomps + c);
            auto *t_6 = buffer.data(target + (r * 11 + 6) * ncomps + c);
            auto *t_7 = buffer.data(target + (r * 11 + 7) * ncomps + c);
            auto *t_8 = buffer.data(target + (r * 11 + 8) * ncomps + c);
            auto *t_9 = buffer.data(target + (r * 11 + 9) * ncomps + c);
            auto *t_10 = buffer.data(target + (r * 11 + 10) * ncomps + c);

            const auto *s_0 = buffer.data(source + (r * 21 + 0) * ncomps + c);
            const auto *s_1 = buffer.data(source + (r * 21 + 1) * ncomps + c);
            const auto *s_2 = buffer.data(source + (r * 21 + 2) * ncomps + c);
            const auto *s_3 = buffer.data(source + (r * 21 + 3) * ncomps + c);
            const auto *s_4 = buffer.data(source + (r * 21 + 4) * ncomps + c);
            const auto *s_5 = buffer.data(source + (r * 21 + 5) * ncomps + c);
            const auto *s_6 = buffer.data(source + (r * 21 + 6) * ncomps + c);
            const auto *s_7 = buffer.data(source + (r * 21 + 7) * ncomps + c);
            const auto *s_8 = buffer.data(source + (r * 21 + 8) * ncomps + c);
            const auto *s_9 = buffer.data(source + (r * 21 + 9) * ncomps + c);
            const auto *s_10 = buffer.data(source + (r * 21 + 10) * ncomps + c);
            const auto *s_11 = buffer.data(source + (r * 21 + 11) * ncomps + c);
            const auto *s_12 = buffer.data(source + (r * 21 + 12) * ncomps + c);
            const auto *s_13 = buffer.data(source + (r * 21 + 13) * ncomps + c);
            const auto *s_14 = buffer.data(source + (r * 21 + 14) * ncomps + c);
            const auto *s_15 = buffer.data(source + (r * 21 + 15) * ncomps + c);
            const auto *s_16 = buffer.data(source + (r * 21 + 16) * ncomps + c);
            const auto *s_17 = buffer.data(source + (r * 21 + 17) * ncomps + c);
            const auto *s_18 = buffer.data(source + (r * 21 + 18) * ncomps + c);
            const auto *s_19 = buffer.data(source + (r * 21 + 19) * ncomps + c);
            const auto *s_20 = buffer.data(source + (r * 21 + 20) * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, s_1, s_4, s_6, s_8, s_11, s_13, s_15, s_17, \
                         s_19 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_0[k] = f_0 * s_1[k]
                         - f_1 * s_6[k]
                         + f_2 * s_15[k];

                t_1[k] = f_3 * s_4[k]
                         - f_3 * s_11[k];

                t_2[k] = -f_4 * s_1[k]
                         - f_5 * s_6[k]
                         + f_6 * s_8[k]
                         + f_7 * s_15[k]
                         - f_8 * s_17[k];

                t_3[k] = -f_9 * s_4[k]
                         - f_9 * s_11[k]
                         + f_10 * s_13[k];

                t_4[k] = f_11 * s_1[k]
                         + f_12 * s_6[k]
                         - f_13 * s_8[k]
                         + f_11 * s_15[k]
                         - f_13 * s_17[k]
                         + f_14 * s_19[k];
            }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, s_0, s_2, s_3, s_5, s_7, s_9, s_10, s_12, s_14, \
                         s_16, s_18, s_20 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_5[k] = 1.875 * s_2[k]
                         + 3.75 * s_7[k]
                         - 5.0 * s_9[k]
                         + 1.875 * s_16[k]
                         - 5.0 * s_18[k]
                         + s_20[k];

                t_6[k] = f_11 * s_0[k]
                         + f_12 * s_3[k]
                         - f_13 * s_5[k]
                         + f_11 * s_10[k]
                         - f_13 * s_12[k]
                         + f_14 * s_14[k];

                t_7[k] = -f_15 * s_2[k]
                         + f_9 * s_9[k]
                         + f_15 * s_16[k]
                         - f_9 * s_18[k];

                t_8[k] = -f_7 * s_0[k]
                         + f_5 * s_3[k]
                         + f_8 * s_5[k]
                         + f_4 * s_10[k]
                         - f_6 * s_12[k];
            }

#pragma omp simd aligned(t_9, t_10, s_0, s_2, s_3, s_7, s_10, s_16 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_9[k] = f_16 * s_2[k]
                         - f_17 * s_7[k]
                         + f_16 * s_16[k];

                t_10[k] = f_2 * s_0[k]
                          - f_1 * s_3[k]
                          + f_0 * s_10[k];
            }
        }
    }
}

auto
transform_h_outer(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t source,
                  const size_t ncomps, const size_t nmax) -> void
{
    // NOTE: the factors are the shell's own, so they are formed once rather than
    // for every atom pair the shell pair reaches.

    const auto f_0 = 0.9375 * std::sqrt(14.0);
    const auto f_1 = 1.875 * std::sqrt(14.0);
    const auto f_2 = 0.1875 * std::sqrt(14.0);
    const auto f_3 = 1.5 * std::sqrt(35.0);
    const auto f_4 = 0.1875 * std::sqrt(70.0);
    const auto f_5 = 0.125 * std::sqrt(70.0);
    const auto f_6 = 1.5 * std::sqrt(70.0);
    const auto f_7 = 0.0625 * std::sqrt(70.0);
    const auto f_8 = 0.5 * std::sqrt(70.0);
    const auto f_9 = 0.5 * std::sqrt(105.0);
    const auto f_10 = std::sqrt(105.0);
    const auto f_11 = 0.125 * std::sqrt(15.0);
    const auto f_12 = 0.25 * std::sqrt(15.0);
    const auto f_13 = 1.5 * std::sqrt(15.0);
    const auto f_14 = std::sqrt(15.0);
    const auto f_15 = 0.25 * std::sqrt(105.0);
    const auto f_16 = 0.375 * std::sqrt(35.0);
    const auto f_17 = 2.25 * std::sqrt(35.0);

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

#pragma omp simd aligned(s_1, s_4, s_6, s_8, s_11, s_13, s_15, s_17, \
                         s_19 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = f_0 * s_1[k]
                     - f_1 * s_6[k]
                     + f_2 * s_15[k];

            g_1[k] = f_3 * s_4[k]
                     - f_3 * s_11[k];

            g_2[k] = -f_4 * s_1[k]
                     - f_5 * s_6[k]
                     + f_6 * s_8[k]
                     + f_7 * s_15[k]
                     - f_8 * s_17[k];

            g_3[k] = -f_9 * s_4[k]
                     - f_9 * s_11[k]
                     + f_10 * s_13[k];

            g_4[k] = f_11 * s_1[k]
                     + f_12 * s_6[k]
                     - f_13 * s_8[k]
                     + f_11 * s_15[k]
                     - f_13 * s_17[k]
                     + f_14 * s_19[k];
        }

#pragma omp simd aligned(s_0, s_2, s_3, s_5, s_7, s_9, s_10, s_12, s_14, s_16, s_18, \
                         s_20 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_5[k] = 1.875 * s_2[k]
                     + 3.75 * s_7[k]
                     - 5.0 * s_9[k]
                     + 1.875 * s_16[k]
                     - 5.0 * s_18[k]
                     + s_20[k];

            g_6[k] = f_11 * s_0[k]
                     + f_12 * s_3[k]
                     - f_13 * s_5[k]
                     + f_11 * s_10[k]
                     - f_13 * s_12[k]
                     + f_14 * s_14[k];

            g_7[k] = -f_15 * s_2[k]
                     + f_9 * s_9[k]
                     + f_15 * s_16[k]
                     - f_9 * s_18[k];

            g_8[k] = -f_7 * s_0[k]
                     + f_5 * s_3[k]
                     + f_8 * s_5[k]
                     + f_4 * s_10[k]
                     - f_6 * s_12[k];
        }

#pragma omp simd aligned(s_0, s_2, s_3, s_7, s_10, s_16 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_9[k] = f_16 * s_2[k]
                     - f_17 * s_7[k]
                     + f_16 * s_16[k];

            g_10[k] = f_2 * s_0[k]
                      - f_1 * s_3[k]
                      + f_0 * s_10[k];
        }
    }
}

auto
transform_h_outer_tri(double *values, const size_t nvalues, CSimdMatrix &buffer,
                      const size_t source, const size_t nmax) -> void
{
    // NOTE: the factors are the shell's own, so they are formed once rather than
    // for every atom pair the shell pair reaches.

    const auto f_0 = 0.9375 * std::sqrt(14.0);
    const auto f_1 = 1.875 * std::sqrt(14.0);
    const auto f_2 = 0.1875 * std::sqrt(14.0);
    const auto f_3 = 1.5 * std::sqrt(35.0);
    const auto f_4 = 0.1875 * std::sqrt(70.0);
    const auto f_5 = 0.125 * std::sqrt(70.0);
    const auto f_6 = 1.5 * std::sqrt(70.0);
    const auto f_7 = 0.0625 * std::sqrt(70.0);
    const auto f_8 = 0.5 * std::sqrt(70.0);
    const auto f_9 = 0.5 * std::sqrt(105.0);
    const auto f_10 = std::sqrt(105.0);
    const auto f_11 = 0.125 * std::sqrt(15.0);
    const auto f_12 = 0.25 * std::sqrt(15.0);
    const auto f_13 = 1.5 * std::sqrt(15.0);
    const auto f_14 = std::sqrt(15.0);
    const auto f_15 = 0.25 * std::sqrt(105.0);
    const auto f_16 = 0.375 * std::sqrt(35.0);
    const auto f_17 = 2.25 * std::sqrt(35.0);

    // NOTE: the rows of the values are not aligned, starting at this combination's
    // offset in the values block, so they are kept out of the clause below.

    // NOTE: the block is its own transpose, so a row above the diagonal is the
    // one below it read the other way round and is copied rather than computed.

    auto *d_0 = values + 0 * nvalues;
    auto *d_1 = values + 12 * nvalues;
    auto *d_2 = values + 24 * nvalues;
    auto *d_3 = values + 36 * nvalues;
    auto *d_4 = values + 48 * nvalues;
    auto *d_5 = values + 60 * nvalues;
    auto *d_6 = values + 72 * nvalues;
    auto *d_7 = values + 84 * nvalues;
    auto *d_8 = values + 96 * nvalues;
    auto *d_9 = values + 108 * nvalues;
    auto *d_10 = values + 120 * nvalues;

    const auto *q_0_1 = buffer.data(source + 11);
    const auto *q_0_6 = buffer.data(source + 66);
    const auto *q_0_15 = buffer.data(source + 165);
    const auto *q_1_4 = buffer.data(source + 45);
    const auto *q_1_11 = buffer.data(source + 122);
    const auto *q_2_1 = buffer.data(source + 13);
    const auto *q_2_6 = buffer.data(source + 68);
    const auto *q_2_8 = buffer.data(source + 90);
    const auto *q_2_15 = buffer.data(source + 167);
    const auto *q_2_17 = buffer.data(source + 189);
    const auto *q_3_4 = buffer.data(source + 47);
    const auto *q_3_11 = buffer.data(source + 124);
    const auto *q_3_13 = buffer.data(source + 146);
    const auto *q_4_1 = buffer.data(source + 15);
    const auto *q_4_6 = buffer.data(source + 70);
    const auto *q_4_8 = buffer.data(source + 92);
    const auto *q_4_15 = buffer.data(source + 169);
    const auto *q_4_17 = buffer.data(source + 191);
    const auto *q_4_19 = buffer.data(source + 213);
    const auto *q_5_2 = buffer.data(source + 27);
    const auto *q_5_7 = buffer.data(source + 82);
    const auto *q_5_9 = buffer.data(source + 104);
    const auto *q_5_16 = buffer.data(source + 181);
    const auto *q_5_18 = buffer.data(source + 203);
    const auto *q_5_20 = buffer.data(source + 225);
    const auto *q_6_0 = buffer.data(source + 6);
    const auto *q_6_3 = buffer.data(source + 39);
    const auto *q_6_5 = buffer.data(source + 61);
    const auto *q_6_10 = buffer.data(source + 116);
    const auto *q_6_12 = buffer.data(source + 138);
    const auto *q_6_14 = buffer.data(source + 160);
    const auto *q_7_2 = buffer.data(source + 29);
    const auto *q_7_9 = buffer.data(source + 106);
    const auto *q_7_16 = buffer.data(source + 183);
    const auto *q_7_18 = buffer.data(source + 205);
    const auto *q_8_0 = buffer.data(source + 8);
    const auto *q_8_3 = buffer.data(source + 41);
    const auto *q_8_5 = buffer.data(source + 63);
    const auto *q_8_10 = buffer.data(source + 118);
    const auto *q_8_12 = buffer.data(source + 140);
    const auto *q_9_2 = buffer.data(source + 31);
    const auto *q_9_7 = buffer.data(source + 86);
    const auto *q_9_16 = buffer.data(source + 185);
    const auto *q_10_0 = buffer.data(source + 10);
    const auto *q_10_3 = buffer.data(source + 43);
    const auto *q_10_10 = buffer.data(source + 120);

#pragma omp simd aligned(q_0_1, q_0_6, q_0_15, q_1_4, q_1_11, q_2_1, q_2_6, q_2_8, q_2_15, \
                         q_2_17 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        d_0[k] = f_0 * q_0_1[k]
                 - f_1 * q_0_6[k]
                 + f_2 * q_0_15[k];

        d_1[k] = f_3 * q_1_4[k]
                 - f_3 * q_1_11[k];

        d_2[k] = -f_4 * q_2_1[k]
                 - f_5 * q_2_6[k]
                 + f_6 * q_2_8[k]
                 + f_7 * q_2_15[k]
                 - f_8 * q_2_17[k];
    }

#pragma omp simd aligned(q_3_4, q_3_11, q_3_13, q_4_1, q_4_6, q_4_8, q_4_15, q_4_17, \
                         q_4_19 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        d_3[k] = -f_9 * q_3_4[k]
                 - f_9 * q_3_11[k]
                 + f_10 * q_3_13[k];

        d_4[k] = f_11 * q_4_1[k]
                 + f_12 * q_4_6[k]
                 - f_13 * q_4_8[k]
                 + f_11 * q_4_15[k]
                 - f_13 * q_4_17[k]
                 + f_14 * q_4_19[k];
    }

#pragma omp simd aligned(q_5_2, q_5_7, q_5_9, q_5_16, q_5_18, q_5_20, q_6_0, q_6_3, q_6_5, \
                         q_6_10, q_6_12, q_6_14 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        d_5[k] = 1.875 * q_5_2[k]
                 + 3.75 * q_5_7[k]
                 - 5.0 * q_5_9[k]
                 + 1.875 * q_5_16[k]
                 - 5.0 * q_5_18[k]
                 + q_5_20[k];

        d_6[k] = f_11 * q_6_0[k]
                 + f_12 * q_6_3[k]
                 - f_13 * q_6_5[k]
                 + f_11 * q_6_10[k]
                 - f_13 * q_6_12[k]
                 + f_14 * q_6_14[k];
    }

#pragma omp simd aligned(q_7_2, q_7_9, q_7_16, q_7_18, q_8_0, q_8_3, q_8_5, q_8_10, q_8_12, \
                         q_9_2, q_9_7, q_9_16 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        d_7[k] = -f_15 * q_7_2[k]
                 + f_9 * q_7_9[k]
                 + f_15 * q_7_16[k]
                 - f_9 * q_7_18[k];

        d_8[k] = -f_7 * q_8_0[k]
                 + f_5 * q_8_3[k]
                 + f_8 * q_8_5[k]
                 + f_4 * q_8_10[k]
                 - f_6 * q_8_12[k];

        d_9[k] = f_16 * q_9_2[k]
                 - f_17 * q_9_7[k]
                 + f_16 * q_9_16[k];
    }

#pragma omp simd aligned(q_10_0, q_10_3, q_10_10 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        d_10[k] = f_2 * q_10_0[k]
                  - f_1 * q_10_3[k]
                  + f_0 * q_10_10[k];
    }

    for (size_t c = 1; c < 11; c++)
    {
        auto *g_0 = values + (0 + c) * nvalues;
        auto *g_1 = values + (c * 11 + 0) * nvalues;

        const auto *s_1 = buffer.data(source + 11 + c);
        const auto *s_6 = buffer.data(source + 66 + c);
        const auto *s_15 = buffer.data(source + 165 + c);

#pragma omp simd aligned(s_1, s_6, s_15 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = f_0 * s_1[k]
                     - f_1 * s_6[k]
                     + f_2 * s_15[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 2; c < 11; c++)
    {
        auto *g_0 = values + (11 + c) * nvalues;
        auto *g_1 = values + (c * 11 + 1) * nvalues;

        const auto *s_4 = buffer.data(source + 44 + c);
        const auto *s_11 = buffer.data(source + 121 + c);

#pragma omp simd aligned(s_4, s_11 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = f_3 * s_4[k]
                     - f_3 * s_11[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 3; c < 11; c++)
    {
        auto *g_0 = values + (22 + c) * nvalues;
        auto *g_1 = values + (c * 11 + 2) * nvalues;

        const auto *s_1 = buffer.data(source + 11 + c);
        const auto *s_6 = buffer.data(source + 66 + c);
        const auto *s_8 = buffer.data(source + 88 + c);
        const auto *s_15 = buffer.data(source + 165 + c);
        const auto *s_17 = buffer.data(source + 187 + c);

#pragma omp simd aligned(s_1, s_6, s_8, s_15, s_17 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = -f_4 * s_1[k]
                     - f_5 * s_6[k]
                     + f_6 * s_8[k]
                     + f_7 * s_15[k]
                     - f_8 * s_17[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 4; c < 11; c++)
    {
        auto *g_0 = values + (33 + c) * nvalues;
        auto *g_1 = values + (c * 11 + 3) * nvalues;

        const auto *s_4 = buffer.data(source + 44 + c);
        const auto *s_11 = buffer.data(source + 121 + c);
        const auto *s_13 = buffer.data(source + 143 + c);

#pragma omp simd aligned(s_4, s_11, s_13 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = -f_9 * s_4[k]
                     - f_9 * s_11[k]
                     + f_10 * s_13[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 5; c < 11; c++)
    {
        auto *g_0 = values + (44 + c) * nvalues;
        auto *g_1 = values + (c * 11 + 4) * nvalues;

        const auto *s_1 = buffer.data(source + 11 + c);
        const auto *s_6 = buffer.data(source + 66 + c);
        const auto *s_8 = buffer.data(source + 88 + c);
        const auto *s_15 = buffer.data(source + 165 + c);
        const auto *s_17 = buffer.data(source + 187 + c);
        const auto *s_19 = buffer.data(source + 209 + c);

#pragma omp simd aligned(s_1, s_6, s_8, s_15, s_17, s_19 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = f_11 * s_1[k]
                     + f_12 * s_6[k]
                     - f_13 * s_8[k]
                     + f_11 * s_15[k]
                     - f_13 * s_17[k]
                     + f_14 * s_19[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 6; c < 11; c++)
    {
        auto *g_0 = values + (55 + c) * nvalues;
        auto *g_1 = values + (c * 11 + 5) * nvalues;

        const auto *s_2 = buffer.data(source + 22 + c);
        const auto *s_7 = buffer.data(source + 77 + c);
        const auto *s_9 = buffer.data(source + 99 + c);
        const auto *s_16 = buffer.data(source + 176 + c);
        const auto *s_18 = buffer.data(source + 198 + c);
        const auto *s_20 = buffer.data(source + 220 + c);

#pragma omp simd aligned(s_2, s_7, s_9, s_16, s_18, s_20 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = 1.875 * s_2[k]
                     + 3.75 * s_7[k]
                     - 5.0 * s_9[k]
                     + 1.875 * s_16[k]
                     - 5.0 * s_18[k]
                     + s_20[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 7; c < 11; c++)
    {
        auto *g_0 = values + (66 + c) * nvalues;
        auto *g_1 = values + (c * 11 + 6) * nvalues;

        const auto *s_0 = buffer.data(source + 0 + c);
        const auto *s_3 = buffer.data(source + 33 + c);
        const auto *s_5 = buffer.data(source + 55 + c);
        const auto *s_10 = buffer.data(source + 110 + c);
        const auto *s_12 = buffer.data(source + 132 + c);
        const auto *s_14 = buffer.data(source + 154 + c);

#pragma omp simd aligned(s_0, s_3, s_5, s_10, s_12, s_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = f_11 * s_0[k]
                     + f_12 * s_3[k]
                     - f_13 * s_5[k]
                     + f_11 * s_10[k]
                     - f_13 * s_12[k]
                     + f_14 * s_14[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 8; c < 11; c++)
    {
        auto *g_0 = values + (77 + c) * nvalues;
        auto *g_1 = values + (c * 11 + 7) * nvalues;

        const auto *s_2 = buffer.data(source + 22 + c);
        const auto *s_9 = buffer.data(source + 99 + c);
        const auto *s_16 = buffer.data(source + 176 + c);
        const auto *s_18 = buffer.data(source + 198 + c);

#pragma omp simd aligned(s_2, s_9, s_16, s_18 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = -f_15 * s_2[k]
                     + f_9 * s_9[k]
                     + f_15 * s_16[k]
                     - f_9 * s_18[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 9; c < 11; c++)
    {
        auto *g_0 = values + (88 + c) * nvalues;
        auto *g_1 = values + (c * 11 + 8) * nvalues;

        const auto *s_0 = buffer.data(source + 0 + c);
        const auto *s_3 = buffer.data(source + 33 + c);
        const auto *s_5 = buffer.data(source + 55 + c);
        const auto *s_10 = buffer.data(source + 110 + c);
        const auto *s_12 = buffer.data(source + 132 + c);

#pragma omp simd aligned(s_0, s_3, s_5, s_10, s_12 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = -f_7 * s_0[k]
                     + f_5 * s_3[k]
                     + f_8 * s_5[k]
                     + f_4 * s_10[k]
                     - f_6 * s_12[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 10; c < 11; c++)
    {
        auto *g_0 = values + (99 + c) * nvalues;
        auto *g_1 = values + (c * 11 + 9) * nvalues;

        const auto *s_2 = buffer.data(source + 22 + c);
        const auto *s_7 = buffer.data(source + 77 + c);
        const auto *s_16 = buffer.data(source + 176 + c);

#pragma omp simd aligned(s_2, s_7, s_16 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = f_16 * s_2[k]
                     - f_17 * s_7[k]
                     + f_16 * s_16[k];
            g_1[k] = g_0[k];
        }
    }
}

}  // namespace simdtrf
