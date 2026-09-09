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


#include "SimdTransformP.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_p_inner(CSimdMatrix &buffer, const size_t target, const size_t source,
                  const size_t nrows, const size_t ncomps, const size_t ncols) -> void
{
    // NOTE: what sits either side of this index reaches the pass as a count and
    // nothing else -- the rows above it, the components below -- so one routine
    // per shell serves every block the shell appears in.

    for (size_t r = 0; r < nrows; r++)
    {
        for (size_t c = 0; c < ncomps; c++)
        {
            auto *t_0 = buffer.data(target + (r * 3 + 0) * ncomps + c);
            auto *t_1 = buffer.data(target + (r * 3 + 1) * ncomps + c);
            auto *t_2 = buffer.data(target + (r * 3 + 2) * ncomps + c);

            const auto *s_0 = buffer.data(source + (r * 3 + 0) * ncomps + c);
            const auto *s_1 = buffer.data(source + (r * 3 + 1) * ncomps + c);
            const auto *s_2 = buffer.data(source + (r * 3 + 2) * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, s_0, s_1, s_2 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_0[k] = s_1[k];

                t_1[k] = s_2[k];

                t_2[k] = s_0[k];
            }
        }
    }
}

auto
transform_p_outer(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t source,
                  const size_t ncomps, const size_t nmax) -> void
{
    // NOTE: the rows of the values are not aligned, starting at this combination's
    // offset in the values block, so they are kept out of the clause below.

    // NOTE: what the other side has left reaches this pass as a count of
    // components, which is one where that side is a single function.

    for (size_t c = 0; c < ncomps; c++)
    {
        auto *g_0 = values + (0 * ncomps + c) * nvalues;
        auto *g_1 = values + (1 * ncomps + c) * nvalues;
        auto *g_2 = values + (2 * ncomps + c) * nvalues;

        const auto *s_0 = buffer.data(source + 0 * ncomps + c);
        const auto *s_1 = buffer.data(source + 1 * ncomps + c);
        const auto *s_2 = buffer.data(source + 2 * ncomps + c);

#pragma omp simd aligned(s_0, s_1, s_2 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = s_1[k];

            g_1[k] = s_2[k];

            g_2[k] = s_0[k];
        }
    }
}

auto
transform_p_outer_tri(double *values, const size_t nvalues, CSimdMatrix &buffer,
                      const size_t source, const size_t nmax) -> void
{
    // NOTE: the rows of the values are not aligned, starting at this combination's
    // offset in the values block, so they are kept out of the clause below.

    // NOTE: the block is its own transpose, so a row above the diagonal is the
    // one below it read the other way round and is copied rather than computed.

    auto *d_0 = values + 0 * nvalues;
    auto *d_1 = values + 4 * nvalues;
    auto *d_2 = values + 8 * nvalues;

    const auto *q_0_1 = buffer.data(source + 3);
    const auto *q_1_2 = buffer.data(source + 7);
    const auto *q_2_0 = buffer.data(source + 2);

#pragma omp simd aligned(q_0_1, q_1_2, q_2_0 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        d_0[k] = q_0_1[k];

        d_1[k] = q_1_2[k];

        d_2[k] = q_2_0[k];
    }

    for (size_t c = 1; c < 3; c++)
    {
        auto *g_0 = values + (0 + c) * nvalues;
        auto *g_1 = values + (c * 3 + 0) * nvalues;

        const auto *s_1 = buffer.data(source + 3 + c);

#pragma omp simd aligned(s_1 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = s_1[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 2; c < 3; c++)
    {
        auto *g_0 = values + (3 + c) * nvalues;
        auto *g_1 = values + (c * 3 + 1) * nvalues;

        const auto *s_2 = buffer.data(source + 6 + c);

#pragma omp simd aligned(s_2 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = s_2[k];
            g_1[k] = g_0[k];
        }
    }
}

}  // namespace simdtrf
