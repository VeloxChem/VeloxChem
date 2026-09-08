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


#include "SimdTransformD.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_d_inner(CSimdMatrix &buffer, const size_t target, const size_t source,
                  const size_t nrows, const size_t ncols) -> void
{
    // NOTE: the factors are the shell's own, so they are formed once rather than
    // for every atom pair the shell pair reaches.

    const auto f_0 = std::sqrt(3.0);
    const auto f_1 = 0.5 * std::sqrt(3.0);

    // NOTE: the other side of the pair reaches this pass as a count of rows and
    // nothing else, its own components running fastest within each of them.

    for (size_t r = 0; r < nrows; r++)
    {
        auto *t_0 = buffer.data(target + r * 5 + 0);
        auto *t_1 = buffer.data(target + r * 5 + 1);
        auto *t_2 = buffer.data(target + r * 5 + 2);
        auto *t_3 = buffer.data(target + r * 5 + 3);
        auto *t_4 = buffer.data(target + r * 5 + 4);

        const auto *s_0 = buffer.data(source + r * 6 + 0);
        const auto *s_1 = buffer.data(source + r * 6 + 1);
        const auto *s_2 = buffer.data(source + r * 6 + 2);
        const auto *s_3 = buffer.data(source + r * 6 + 3);
        const auto *s_4 = buffer.data(source + r * 6 + 4);
        const auto *s_5 = buffer.data(source + r * 6 + 5);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, s_0, s_1, s_2, s_3, s_4, \
                         s_5 : simd::cache_line_size())
        for (size_t k = 0; k < ncols; k++)
        {
            t_0[k] = f_0 * s_1[k];

            t_1[k] = f_0 * s_4[k];

            t_2[k] = -0.5 * s_0[k]
                     - 0.5 * s_3[k]
                     + s_5[k];

            t_3[k] = f_0 * s_2[k];

            t_4[k] = f_1 * s_0[k]
                     - f_1 * s_3[k];
        }
    }
}

auto
transform_d_outer(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t source,
                  const size_t ncomps, const size_t nmax) -> void
{
    // NOTE: the factors are the shell's own, so they are formed once rather than
    // for every atom pair the shell pair reaches.

    const auto f_0 = std::sqrt(3.0);
    const auto f_1 = 0.5 * std::sqrt(3.0);

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

        const auto *s_0 = buffer.data(source + 0 * ncomps + c);
        const auto *s_1 = buffer.data(source + 1 * ncomps + c);
        const auto *s_2 = buffer.data(source + 2 * ncomps + c);
        const auto *s_3 = buffer.data(source + 3 * ncomps + c);
        const auto *s_4 = buffer.data(source + 4 * ncomps + c);
        const auto *s_5 = buffer.data(source + 5 * ncomps + c);

#pragma omp simd aligned(s_0, s_1, s_2, s_3, s_4, s_5 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = f_0 * s_1[k];

            g_1[k] = f_0 * s_4[k];

            g_2[k] = -0.5 * s_0[k]
                     - 0.5 * s_3[k]
                     + s_5[k];

            g_3[k] = f_0 * s_2[k];

            g_4[k] = f_1 * s_0[k]
                     - f_1 * s_3[k];
        }
    }
}

auto
transform_d_outer_tri(double *values, const size_t nvalues, CSimdMatrix &buffer,
                      const size_t source, const size_t nmax) -> void
{
    // NOTE: the factors are the shell's own, so they are formed once rather than
    // for every atom pair the shell pair reaches.

    const auto f_0 = std::sqrt(3.0);
    const auto f_1 = 0.5 * std::sqrt(3.0);

    // NOTE: the rows of the values are not aligned, starting at this combination's
    // offset in the values block, so they are kept out of the clause below.

    // NOTE: the block is its own transpose, so a row above the diagonal is the
    // one below it read the other way round and is copied rather than computed.

    auto *d_0 = values + 0 * nvalues;
    auto *d_1 = values + 6 * nvalues;
    auto *d_2 = values + 12 * nvalues;
    auto *d_3 = values + 18 * nvalues;
    auto *d_4 = values + 24 * nvalues;

    const auto *q_0_1 = buffer.data(source + 5);
    const auto *q_1_4 = buffer.data(source + 21);
    const auto *q_2_0 = buffer.data(source + 2);
    const auto *q_2_3 = buffer.data(source + 17);
    const auto *q_2_5 = buffer.data(source + 27);
    const auto *q_3_2 = buffer.data(source + 13);
    const auto *q_4_0 = buffer.data(source + 4);
    const auto *q_4_3 = buffer.data(source + 19);

#pragma omp simd aligned(q_0_1, q_1_4, q_2_0, q_2_3, q_2_5, q_3_2, q_4_0, \
                         q_4_3 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        d_0[k] = f_0 * q_0_1[k];

        d_1[k] = f_0 * q_1_4[k];

        d_2[k] = -0.5 * q_2_0[k]
                 - 0.5 * q_2_3[k]
                 + q_2_5[k];

        d_3[k] = f_0 * q_3_2[k];

        d_4[k] = f_1 * q_4_0[k]
                 - f_1 * q_4_3[k];
    }

    for (size_t c = 1; c < 5; c++)
    {
        auto *g_0 = values + (0 + c) * nvalues;
        auto *g_1 = values + (c * 5 + 0) * nvalues;

        const auto *s_1 = buffer.data(source + 5 + c);

#pragma omp simd aligned(s_1 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = f_0 * s_1[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 2; c < 5; c++)
    {
        auto *g_0 = values + (5 + c) * nvalues;
        auto *g_1 = values + (c * 5 + 1) * nvalues;

        const auto *s_4 = buffer.data(source + 20 + c);

#pragma omp simd aligned(s_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = f_0 * s_4[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 3; c < 5; c++)
    {
        auto *g_0 = values + (10 + c) * nvalues;
        auto *g_1 = values + (c * 5 + 2) * nvalues;

        const auto *s_0 = buffer.data(source + 0 + c);
        const auto *s_3 = buffer.data(source + 15 + c);
        const auto *s_5 = buffer.data(source + 25 + c);

#pragma omp simd aligned(s_0, s_3, s_5 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = -0.5 * s_0[k]
                     - 0.5 * s_3[k]
                     + s_5[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 4; c < 5; c++)
    {
        auto *g_0 = values + (15 + c) * nvalues;
        auto *g_1 = values + (c * 5 + 3) * nvalues;

        const auto *s_2 = buffer.data(source + 10 + c);

#pragma omp simd aligned(s_2 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = f_0 * s_2[k];
            g_1[k] = g_0[k];
        }
    }
}

}  // namespace simdtrf
