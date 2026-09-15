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


#include "SimdGeometryP1.hpp"

#include "SimdAlign.hpp"

namespace simdgeo {  // simdgeo namespace

auto
geom_p_x(CSimdMatrix &buffer, const size_t target, const size_t s0, const size_t s1,
         const size_t nrows, const size_t ncomps, const size_t ncols,
         const double exponent) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * exponent;

    // NOTE: what sits either side of the differentiated function reaches the pass
    // as a count and nothing else -- the rows above it, the components below -- so
    // one routine per shell serves every block the shell appears in, whatever the
    // other functions carry and whatever the operator is.

    for (size_t r = 0; r < nrows; r++)
    {
        for (size_t c = 0; c < ncomps; c++)
        {
            auto *t_0 = buffer.data(target + (r * 3 + 0) * ncomps + c);
            auto *t_1 = buffer.data(target + (r * 3 + 1) * ncomps + c);
            auto *t_2 = buffer.data(target + (r * 3 + 2) * ncomps + c);
            const auto *s0_0 = buffer.data(s0 + (r * 1 + 0) * ncomps + c);
            const auto *s1_0 = buffer.data(s1 + (r * 6 + 0) * ncomps + c);
            const auto *s1_1 = buffer.data(s1 + (r * 6 + 1) * ncomps + c);
            const auto *s1_2 = buffer.data(s1 + (r * 6 + 2) * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, s0_0, s1_0, s1_1, s1_2 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_0[k] = -s0_0[k]
                         + f_0 * s1_0[k];

                t_1[k] = f_0 * s1_1[k];

                t_2[k] = f_0 * s1_2[k];
            }
        }
    }
}

auto
geom_p_y(CSimdMatrix &buffer, const size_t target, const size_t s0, const size_t s1,
         const size_t nrows, const size_t ncomps, const size_t ncols,
         const double exponent) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * exponent;

    // NOTE: what sits either side of the differentiated function reaches the pass
    // as a count and nothing else -- the rows above it, the components below -- so
    // one routine per shell serves every block the shell appears in, whatever the
    // other functions carry and whatever the operator is.

    for (size_t r = 0; r < nrows; r++)
    {
        for (size_t c = 0; c < ncomps; c++)
        {
            auto *t_0 = buffer.data(target + (r * 3 + 0) * ncomps + c);
            auto *t_1 = buffer.data(target + (r * 3 + 1) * ncomps + c);
            auto *t_2 = buffer.data(target + (r * 3 + 2) * ncomps + c);
            const auto *s0_0 = buffer.data(s0 + (r * 1 + 0) * ncomps + c);
            const auto *s1_1 = buffer.data(s1 + (r * 6 + 1) * ncomps + c);
            const auto *s1_3 = buffer.data(s1 + (r * 6 + 3) * ncomps + c);
            const auto *s1_4 = buffer.data(s1 + (r * 6 + 4) * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, s0_0, s1_1, s1_3, s1_4 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_0[k] = f_0 * s1_1[k];

                t_1[k] = -s0_0[k]
                         + f_0 * s1_3[k];

                t_2[k] = f_0 * s1_4[k];
            }
        }
    }
}

auto
geom_p_z(CSimdMatrix &buffer, const size_t target, const size_t s0, const size_t s1,
         const size_t nrows, const size_t ncomps, const size_t ncols,
         const double exponent) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * exponent;

    // NOTE: what sits either side of the differentiated function reaches the pass
    // as a count and nothing else -- the rows above it, the components below -- so
    // one routine per shell serves every block the shell appears in, whatever the
    // other functions carry and whatever the operator is.

    for (size_t r = 0; r < nrows; r++)
    {
        for (size_t c = 0; c < ncomps; c++)
        {
            auto *t_0 = buffer.data(target + (r * 3 + 0) * ncomps + c);
            auto *t_1 = buffer.data(target + (r * 3 + 1) * ncomps + c);
            auto *t_2 = buffer.data(target + (r * 3 + 2) * ncomps + c);
            const auto *s0_0 = buffer.data(s0 + (r * 1 + 0) * ncomps + c);
            const auto *s1_2 = buffer.data(s1 + (r * 6 + 2) * ncomps + c);
            const auto *s1_4 = buffer.data(s1 + (r * 6 + 4) * ncomps + c);
            const auto *s1_5 = buffer.data(s1 + (r * 6 + 5) * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, s0_0, s1_2, s1_4, s1_5 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_0[k] = f_0 * s1_2[k];

                t_1[k] = f_0 * s1_4[k];

                t_2[k] = -s0_0[k]
                         + f_0 * s1_5[k];
            }
        }
    }
}

}  // namespace simdgeo
