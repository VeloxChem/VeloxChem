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


#include "SimdTransformF.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_f_inner(CSimdMatrix &buffer, const size_t target, const size_t source,
                  const size_t nrows, const size_t ncomps, const size_t ncols) -> void
{
    // NOTE: the factors are the shell's own, so they are formed once rather than
    // for every atom pair the shell pair reaches.

    const auto f_0 = 0.75 * std::sqrt(10.0);
    const auto f_1 = 0.25 * std::sqrt(10.0);
    const auto f_2 = std::sqrt(15.0);
    const auto f_3 = 0.25 * std::sqrt(6.0);
    const auto f_4 = std::sqrt(6.0);
    const auto f_5 = 0.5 * std::sqrt(15.0);

    // NOTE: what sits either side of this index reaches the pass as a count and
    // nothing else -- the rows above it, the components below -- so one routine
    // per shell serves every block the shell appears in.

    for (size_t r = 0; r < nrows; r++)
    {
        for (size_t c = 0; c < ncomps; c++)
        {
            auto *t_0 = buffer.data(target + (r * 7 + 0) * ncomps + c);
            auto *t_1 = buffer.data(target + (r * 7 + 1) * ncomps + c);
            auto *t_2 = buffer.data(target + (r * 7 + 2) * ncomps + c);
            auto *t_3 = buffer.data(target + (r * 7 + 3) * ncomps + c);
            auto *t_4 = buffer.data(target + (r * 7 + 4) * ncomps + c);
            auto *t_5 = buffer.data(target + (r * 7 + 5) * ncomps + c);
            auto *t_6 = buffer.data(target + (r * 7 + 6) * ncomps + c);

            const auto *s_0 = buffer.data(source + (r * 10 + 0) * ncomps + c);
            const auto *s_1 = buffer.data(source + (r * 10 + 1) * ncomps + c);
            const auto *s_2 = buffer.data(source + (r * 10 + 2) * ncomps + c);
            const auto *s_3 = buffer.data(source + (r * 10 + 3) * ncomps + c);
            const auto *s_4 = buffer.data(source + (r * 10 + 4) * ncomps + c);
            const auto *s_5 = buffer.data(source + (r * 10 + 5) * ncomps + c);
            const auto *s_6 = buffer.data(source + (r * 10 + 6) * ncomps + c);
            const auto *s_7 = buffer.data(source + (r * 10 + 7) * ncomps + c);
            const auto *s_8 = buffer.data(source + (r * 10 + 8) * ncomps + c);
            const auto *s_9 = buffer.data(source + (r * 10 + 9) * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, s_0, s_1, s_2, s_3, s_4, s_5, s_6, s_7, \
                         s_8, s_9 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_0[k] = f_0 * s_1[k]
                         - f_1 * s_6[k];

                t_1[k] = f_2 * s_4[k];

                t_2[k] = -f_3 * s_1[k]
                         - f_3 * s_6[k]
                         + f_4 * s_8[k];

                t_3[k] = -1.5 * s_2[k]
                         - 1.5 * s_7[k]
                         + s_9[k];

                t_4[k] = -f_3 * s_0[k]
                         - f_3 * s_3[k]
                         + f_4 * s_5[k];

                t_5[k] = f_5 * s_2[k]
                         - f_5 * s_7[k];
            }

#pragma omp simd aligned(t_6, s_0, s_3 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_6[k] = f_1 * s_0[k]
                         - f_0 * s_3[k];
            }
        }
    }
}

auto
transform_f_outer(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t source,
                  const size_t ncomps, const size_t nmax) -> void
{
    // NOTE: the factors are the shell's own, so they are formed once rather than
    // for every atom pair the shell pair reaches.

    const auto f_0 = 0.75 * std::sqrt(10.0);
    const auto f_1 = 0.25 * std::sqrt(10.0);
    const auto f_2 = std::sqrt(15.0);
    const auto f_3 = 0.25 * std::sqrt(6.0);
    const auto f_4 = std::sqrt(6.0);
    const auto f_5 = 0.5 * std::sqrt(15.0);

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

#pragma omp simd aligned(s_0, s_1, s_2, s_3, s_4, s_5, s_6, s_7, s_8, \
                         s_9 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = f_0 * s_1[k]
                     - f_1 * s_6[k];

            g_1[k] = f_2 * s_4[k];

            g_2[k] = -f_3 * s_1[k]
                     - f_3 * s_6[k]
                     + f_4 * s_8[k];

            g_3[k] = -1.5 * s_2[k]
                     - 1.5 * s_7[k]
                     + s_9[k];

            g_4[k] = -f_3 * s_0[k]
                     - f_3 * s_3[k]
                     + f_4 * s_5[k];

            g_5[k] = f_5 * s_2[k]
                     - f_5 * s_7[k];
        }

#pragma omp simd aligned(s_0, s_3 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_6[k] = f_1 * s_0[k]
                     - f_0 * s_3[k];
        }
    }
}

auto
transform_f_outer_tri(double *values, const size_t nvalues, CSimdMatrix &buffer,
                      const size_t source, const size_t nmax) -> void
{
    // NOTE: the factors are the shell's own, so they are formed once rather than
    // for every atom pair the shell pair reaches.

    const auto f_0 = 0.75 * std::sqrt(10.0);
    const auto f_1 = 0.25 * std::sqrt(10.0);
    const auto f_2 = std::sqrt(15.0);
    const auto f_3 = 0.25 * std::sqrt(6.0);
    const auto f_4 = std::sqrt(6.0);
    const auto f_5 = 0.5 * std::sqrt(15.0);

    // NOTE: the rows of the values are not aligned, starting at this combination's
    // offset in the values block, so they are kept out of the clause below.

    // NOTE: the block is its own transpose, so a row above the diagonal is the
    // one below it read the other way round and is copied rather than computed.

    auto *d_0 = values + 0 * nvalues;
    auto *d_1 = values + 8 * nvalues;
    auto *d_2 = values + 16 * nvalues;
    auto *d_3 = values + 24 * nvalues;
    auto *d_4 = values + 32 * nvalues;
    auto *d_5 = values + 40 * nvalues;
    auto *d_6 = values + 48 * nvalues;

    const auto *q_0_1 = buffer.data(source + 7);
    const auto *q_0_6 = buffer.data(source + 42);
    const auto *q_1_4 = buffer.data(source + 29);
    const auto *q_2_1 = buffer.data(source + 9);
    const auto *q_2_6 = buffer.data(source + 44);
    const auto *q_2_8 = buffer.data(source + 58);
    const auto *q_3_2 = buffer.data(source + 17);
    const auto *q_3_7 = buffer.data(source + 52);
    const auto *q_3_9 = buffer.data(source + 66);
    const auto *q_4_0 = buffer.data(source + 4);
    const auto *q_4_3 = buffer.data(source + 25);
    const auto *q_4_5 = buffer.data(source + 39);
    const auto *q_5_2 = buffer.data(source + 19);
    const auto *q_5_7 = buffer.data(source + 54);
    const auto *q_6_0 = buffer.data(source + 6);
    const auto *q_6_3 = buffer.data(source + 27);

#pragma omp simd aligned(q_0_1, q_0_6, q_1_4, q_2_1, q_2_6, q_2_8, q_3_2, q_3_7, \
                         q_3_9 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        d_0[k] = f_0 * q_0_1[k]
                 - f_1 * q_0_6[k];

        d_1[k] = f_2 * q_1_4[k];

        d_2[k] = -f_3 * q_2_1[k]
                 - f_3 * q_2_6[k]
                 + f_4 * q_2_8[k];

        d_3[k] = -1.5 * q_3_2[k]
                 - 1.5 * q_3_7[k]
                 + q_3_9[k];
    }

#pragma omp simd aligned(q_4_0, q_4_3, q_4_5, q_5_2, q_5_7, q_6_0, \
                         q_6_3 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        d_4[k] = -f_3 * q_4_0[k]
                 - f_3 * q_4_3[k]
                 + f_4 * q_4_5[k];

        d_5[k] = f_5 * q_5_2[k]
                 - f_5 * q_5_7[k];

        d_6[k] = f_1 * q_6_0[k]
                 - f_0 * q_6_3[k];
    }

    for (size_t c = 1; c < 7; c++)
    {
        auto *g_0 = values + (0 + c) * nvalues;
        auto *g_1 = values + (c * 7 + 0) * nvalues;

        const auto *s_1 = buffer.data(source + 7 + c);
        const auto *s_6 = buffer.data(source + 42 + c);

#pragma omp simd aligned(s_1, s_6 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = f_0 * s_1[k]
                     - f_1 * s_6[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 2; c < 7; c++)
    {
        auto *g_0 = values + (7 + c) * nvalues;
        auto *g_1 = values + (c * 7 + 1) * nvalues;

        const auto *s_4 = buffer.data(source + 28 + c);

#pragma omp simd aligned(s_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = f_2 * s_4[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 3; c < 7; c++)
    {
        auto *g_0 = values + (14 + c) * nvalues;
        auto *g_1 = values + (c * 7 + 2) * nvalues;

        const auto *s_1 = buffer.data(source + 7 + c);
        const auto *s_6 = buffer.data(source + 42 + c);
        const auto *s_8 = buffer.data(source + 56 + c);

#pragma omp simd aligned(s_1, s_6, s_8 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = -f_3 * s_1[k]
                     - f_3 * s_6[k]
                     + f_4 * s_8[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 4; c < 7; c++)
    {
        auto *g_0 = values + (21 + c) * nvalues;
        auto *g_1 = values + (c * 7 + 3) * nvalues;

        const auto *s_2 = buffer.data(source + 14 + c);
        const auto *s_7 = buffer.data(source + 49 + c);
        const auto *s_9 = buffer.data(source + 63 + c);

#pragma omp simd aligned(s_2, s_7, s_9 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = -1.5 * s_2[k]
                     - 1.5 * s_7[k]
                     + s_9[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 5; c < 7; c++)
    {
        auto *g_0 = values + (28 + c) * nvalues;
        auto *g_1 = values + (c * 7 + 4) * nvalues;

        const auto *s_0 = buffer.data(source + 0 + c);
        const auto *s_3 = buffer.data(source + 21 + c);
        const auto *s_5 = buffer.data(source + 35 + c);

#pragma omp simd aligned(s_0, s_3, s_5 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = -f_3 * s_0[k]
                     - f_3 * s_3[k]
                     + f_4 * s_5[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 6; c < 7; c++)
    {
        auto *g_0 = values + (35 + c) * nvalues;
        auto *g_1 = values + (c * 7 + 5) * nvalues;

        const auto *s_2 = buffer.data(source + 14 + c);
        const auto *s_7 = buffer.data(source + 49 + c);

#pragma omp simd aligned(s_2, s_7 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = f_5 * s_2[k]
                     - f_5 * s_7[k];
            g_1[k] = g_0[k];
        }
    }
}

}  // namespace simdtrf
