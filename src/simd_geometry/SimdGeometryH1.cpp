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


#include "SimdGeometryH1.hpp"

#include "SimdAlign.hpp"

namespace simdgeo {  // simdgeo namespace

auto
geom_h_x(CSimdMatrix &buffer, const size_t target, const size_t s0, const size_t s1,
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
            auto *t_0 = buffer.data(target + (r * 21 + 0) * ncomps + c);
            auto *t_1 = buffer.data(target + (r * 21 + 1) * ncomps + c);
            auto *t_2 = buffer.data(target + (r * 21 + 2) * ncomps + c);
            auto *t_3 = buffer.data(target + (r * 21 + 3) * ncomps + c);
            auto *t_4 = buffer.data(target + (r * 21 + 4) * ncomps + c);
            auto *t_5 = buffer.data(target + (r * 21 + 5) * ncomps + c);
            auto *t_6 = buffer.data(target + (r * 21 + 6) * ncomps + c);
            auto *t_7 = buffer.data(target + (r * 21 + 7) * ncomps + c);
            auto *t_8 = buffer.data(target + (r * 21 + 8) * ncomps + c);
            auto *t_9 = buffer.data(target + (r * 21 + 9) * ncomps + c);
            auto *t_10 = buffer.data(target + (r * 21 + 10) * ncomps + c);
            auto *t_11 = buffer.data(target + (r * 21 + 11) * ncomps + c);
            auto *t_12 = buffer.data(target + (r * 21 + 12) * ncomps + c);
            auto *t_13 = buffer.data(target + (r * 21 + 13) * ncomps + c);
            auto *t_14 = buffer.data(target + (r * 21 + 14) * ncomps + c);
            auto *t_15 = buffer.data(target + (r * 21 + 15) * ncomps + c);
            auto *t_16 = buffer.data(target + (r * 21 + 16) * ncomps + c);
            auto *t_17 = buffer.data(target + (r * 21 + 17) * ncomps + c);
            auto *t_18 = buffer.data(target + (r * 21 + 18) * ncomps + c);
            auto *t_19 = buffer.data(target + (r * 21 + 19) * ncomps + c);
            auto *t_20 = buffer.data(target + (r * 21 + 20) * ncomps + c);
            const auto *s0_0 = buffer.data(s0 + (r * 15 + 0) * ncomps + c);
            const auto *s0_1 = buffer.data(s0 + (r * 15 + 1) * ncomps + c);
            const auto *s0_2 = buffer.data(s0 + (r * 15 + 2) * ncomps + c);
            const auto *s0_3 = buffer.data(s0 + (r * 15 + 3) * ncomps + c);
            const auto *s0_4 = buffer.data(s0 + (r * 15 + 4) * ncomps + c);
            const auto *s0_5 = buffer.data(s0 + (r * 15 + 5) * ncomps + c);
            const auto *s0_6 = buffer.data(s0 + (r * 15 + 6) * ncomps + c);
            const auto *s0_7 = buffer.data(s0 + (r * 15 + 7) * ncomps + c);
            const auto *s0_8 = buffer.data(s0 + (r * 15 + 8) * ncomps + c);
            const auto *s0_9 = buffer.data(s0 + (r * 15 + 9) * ncomps + c);
            const auto *s0_10 = buffer.data(s0 + (r * 15 + 10) * ncomps + c);
            const auto *s0_11 = buffer.data(s0 + (r * 15 + 11) * ncomps + c);
            const auto *s0_12 = buffer.data(s0 + (r * 15 + 12) * ncomps + c);
            const auto *s0_13 = buffer.data(s0 + (r * 15 + 13) * ncomps + c);
            const auto *s0_14 = buffer.data(s0 + (r * 15 + 14) * ncomps + c);
            const auto *s1_0 = buffer.data(s1 + (r * 28 + 0) * ncomps + c);
            const auto *s1_1 = buffer.data(s1 + (r * 28 + 1) * ncomps + c);
            const auto *s1_2 = buffer.data(s1 + (r * 28 + 2) * ncomps + c);
            const auto *s1_3 = buffer.data(s1 + (r * 28 + 3) * ncomps + c);
            const auto *s1_4 = buffer.data(s1 + (r * 28 + 4) * ncomps + c);
            const auto *s1_5 = buffer.data(s1 + (r * 28 + 5) * ncomps + c);
            const auto *s1_6 = buffer.data(s1 + (r * 28 + 6) * ncomps + c);
            const auto *s1_7 = buffer.data(s1 + (r * 28 + 7) * ncomps + c);
            const auto *s1_8 = buffer.data(s1 + (r * 28 + 8) * ncomps + c);
            const auto *s1_9 = buffer.data(s1 + (r * 28 + 9) * ncomps + c);
            const auto *s1_10 = buffer.data(s1 + (r * 28 + 10) * ncomps + c);
            const auto *s1_11 = buffer.data(s1 + (r * 28 + 11) * ncomps + c);
            const auto *s1_12 = buffer.data(s1 + (r * 28 + 12) * ncomps + c);
            const auto *s1_13 = buffer.data(s1 + (r * 28 + 13) * ncomps + c);
            const auto *s1_14 = buffer.data(s1 + (r * 28 + 14) * ncomps + c);
            const auto *s1_15 = buffer.data(s1 + (r * 28 + 15) * ncomps + c);
            const auto *s1_16 = buffer.data(s1 + (r * 28 + 16) * ncomps + c);
            const auto *s1_17 = buffer.data(s1 + (r * 28 + 17) * ncomps + c);
            const auto *s1_18 = buffer.data(s1 + (r * 28 + 18) * ncomps + c);
            const auto *s1_19 = buffer.data(s1 + (r * 28 + 19) * ncomps + c);
            const auto *s1_20 = buffer.data(s1 + (r * 28 + 20) * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, s0_0, s0_1, s0_2, s0_3, s0_4, s1_0, s1_1, \
                         s1_2, s1_3, s1_4 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_0[k] = -5.0 * s0_0[k]
                         + f_0 * s1_0[k];

                t_1[k] = -4.0 * s0_1[k]
                         + f_0 * s1_1[k];

                t_2[k] = -4.0 * s0_2[k]
                         + f_0 * s1_2[k];

                t_3[k] = -3.0 * s0_3[k]
                         + f_0 * s1_3[k];

                t_4[k] = -3.0 * s0_4[k]
                         + f_0 * s1_4[k];
            }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, s0_5, s0_6, s0_7, s0_8, s0_9, s1_5, s1_6, \
                         s1_7, s1_8, s1_9 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_5[k] = -3.0 * s0_5[k]
                         + f_0 * s1_5[k];

                t_6[k] = -2.0 * s0_6[k]
                         + f_0 * s1_6[k];

                t_7[k] = -2.0 * s0_7[k]
                         + f_0 * s1_7[k];

                t_8[k] = -2.0 * s0_8[k]
                         + f_0 * s1_8[k];

                t_9[k] = -2.0 * s0_9[k]
                         + f_0 * s1_9[k];
            }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, s0_10, s0_11, s0_12, s0_13, s0_14, \
                         s1_10, s1_11, s1_12, s1_13, s1_14 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_10[k] = -s0_10[k]
                          + f_0 * s1_10[k];

                t_11[k] = -s0_11[k]
                          + f_0 * s1_11[k];

                t_12[k] = -s0_12[k]
                          + f_0 * s1_12[k];

                t_13[k] = -s0_13[k]
                          + f_0 * s1_13[k];

                t_14[k] = -s0_14[k]
                          + f_0 * s1_14[k];
            }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, t_20, s1_15, s1_16, s1_17, s1_18, \
                         s1_19, s1_20 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_15[k] = f_0 * s1_15[k];

                t_16[k] = f_0 * s1_16[k];

                t_17[k] = f_0 * s1_17[k];

                t_18[k] = f_0 * s1_18[k];

                t_19[k] = f_0 * s1_19[k];

                t_20[k] = f_0 * s1_20[k];
            }
        }
    }
}

auto
geom_h_y(CSimdMatrix &buffer, const size_t target, const size_t s0, const size_t s1,
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
            auto *t_0 = buffer.data(target + (r * 21 + 0) * ncomps + c);
            auto *t_1 = buffer.data(target + (r * 21 + 1) * ncomps + c);
            auto *t_2 = buffer.data(target + (r * 21 + 2) * ncomps + c);
            auto *t_3 = buffer.data(target + (r * 21 + 3) * ncomps + c);
            auto *t_4 = buffer.data(target + (r * 21 + 4) * ncomps + c);
            auto *t_5 = buffer.data(target + (r * 21 + 5) * ncomps + c);
            auto *t_6 = buffer.data(target + (r * 21 + 6) * ncomps + c);
            auto *t_7 = buffer.data(target + (r * 21 + 7) * ncomps + c);
            auto *t_8 = buffer.data(target + (r * 21 + 8) * ncomps + c);
            auto *t_9 = buffer.data(target + (r * 21 + 9) * ncomps + c);
            auto *t_10 = buffer.data(target + (r * 21 + 10) * ncomps + c);
            auto *t_11 = buffer.data(target + (r * 21 + 11) * ncomps + c);
            auto *t_12 = buffer.data(target + (r * 21 + 12) * ncomps + c);
            auto *t_13 = buffer.data(target + (r * 21 + 13) * ncomps + c);
            auto *t_14 = buffer.data(target + (r * 21 + 14) * ncomps + c);
            auto *t_15 = buffer.data(target + (r * 21 + 15) * ncomps + c);
            auto *t_16 = buffer.data(target + (r * 21 + 16) * ncomps + c);
            auto *t_17 = buffer.data(target + (r * 21 + 17) * ncomps + c);
            auto *t_18 = buffer.data(target + (r * 21 + 18) * ncomps + c);
            auto *t_19 = buffer.data(target + (r * 21 + 19) * ncomps + c);
            auto *t_20 = buffer.data(target + (r * 21 + 20) * ncomps + c);
            const auto *s0_0 = buffer.data(s0 + (r * 15 + 0) * ncomps + c);
            const auto *s0_1 = buffer.data(s0 + (r * 15 + 1) * ncomps + c);
            const auto *s0_2 = buffer.data(s0 + (r * 15 + 2) * ncomps + c);
            const auto *s0_3 = buffer.data(s0 + (r * 15 + 3) * ncomps + c);
            const auto *s0_4 = buffer.data(s0 + (r * 15 + 4) * ncomps + c);
            const auto *s0_5 = buffer.data(s0 + (r * 15 + 5) * ncomps + c);
            const auto *s0_6 = buffer.data(s0 + (r * 15 + 6) * ncomps + c);
            const auto *s0_7 = buffer.data(s0 + (r * 15 + 7) * ncomps + c);
            const auto *s0_8 = buffer.data(s0 + (r * 15 + 8) * ncomps + c);
            const auto *s0_9 = buffer.data(s0 + (r * 15 + 9) * ncomps + c);
            const auto *s0_10 = buffer.data(s0 + (r * 15 + 10) * ncomps + c);
            const auto *s0_11 = buffer.data(s0 + (r * 15 + 11) * ncomps + c);
            const auto *s0_12 = buffer.data(s0 + (r * 15 + 12) * ncomps + c);
            const auto *s0_13 = buffer.data(s0 + (r * 15 + 13) * ncomps + c);
            const auto *s0_14 = buffer.data(s0 + (r * 15 + 14) * ncomps + c);
            const auto *s1_1 = buffer.data(s1 + (r * 28 + 1) * ncomps + c);
            const auto *s1_3 = buffer.data(s1 + (r * 28 + 3) * ncomps + c);
            const auto *s1_4 = buffer.data(s1 + (r * 28 + 4) * ncomps + c);
            const auto *s1_6 = buffer.data(s1 + (r * 28 + 6) * ncomps + c);
            const auto *s1_7 = buffer.data(s1 + (r * 28 + 7) * ncomps + c);
            const auto *s1_8 = buffer.data(s1 + (r * 28 + 8) * ncomps + c);
            const auto *s1_10 = buffer.data(s1 + (r * 28 + 10) * ncomps + c);
            const auto *s1_11 = buffer.data(s1 + (r * 28 + 11) * ncomps + c);
            const auto *s1_12 = buffer.data(s1 + (r * 28 + 12) * ncomps + c);
            const auto *s1_13 = buffer.data(s1 + (r * 28 + 13) * ncomps + c);
            const auto *s1_15 = buffer.data(s1 + (r * 28 + 15) * ncomps + c);
            const auto *s1_16 = buffer.data(s1 + (r * 28 + 16) * ncomps + c);
            const auto *s1_17 = buffer.data(s1 + (r * 28 + 17) * ncomps + c);
            const auto *s1_18 = buffer.data(s1 + (r * 28 + 18) * ncomps + c);
            const auto *s1_19 = buffer.data(s1 + (r * 28 + 19) * ncomps + c);
            const auto *s1_21 = buffer.data(s1 + (r * 28 + 21) * ncomps + c);
            const auto *s1_22 = buffer.data(s1 + (r * 28 + 22) * ncomps + c);
            const auto *s1_23 = buffer.data(s1 + (r * 28 + 23) * ncomps + c);
            const auto *s1_24 = buffer.data(s1 + (r * 28 + 24) * ncomps + c);
            const auto *s1_25 = buffer.data(s1 + (r * 28 + 25) * ncomps + c);
            const auto *s1_26 = buffer.data(s1 + (r * 28 + 26) * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, s0_0, s0_1, s0_2, s1_1, s1_3, s1_4, \
                         s1_6, s1_7, s1_8 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_0[k] = f_0 * s1_1[k];

                t_1[k] = -s0_0[k]
                         + f_0 * s1_3[k];

                t_2[k] = f_0 * s1_4[k];

                t_3[k] = -2.0 * s0_1[k]
                         + f_0 * s1_6[k];

                t_4[k] = -s0_2[k]
                         + f_0 * s1_7[k];

                t_5[k] = f_0 * s1_8[k];
            }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, s0_3, s0_4, s0_5, s0_6, s1_10, s1_11, \
                         s1_12, s1_13, s1_15 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_6[k] = -3.0 * s0_3[k]
                         + f_0 * s1_10[k];

                t_7[k] = -2.0 * s0_4[k]
                         + f_0 * s1_11[k];

                t_8[k] = -s0_5[k]
                         + f_0 * s1_12[k];

                t_9[k] = f_0 * s1_13[k];

                t_10[k] = -4.0 * s0_6[k]
                          + f_0 * s1_15[k];
            }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, s0_7, s0_8, s0_9, s0_10, s1_16, s1_17, \
                         s1_18, s1_19, s1_21 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_11[k] = -3.0 * s0_7[k]
                          + f_0 * s1_16[k];

                t_12[k] = -2.0 * s0_8[k]
                          + f_0 * s1_17[k];

                t_13[k] = -s0_9[k]
                          + f_0 * s1_18[k];

                t_14[k] = f_0 * s1_19[k];

                t_15[k] = -5.0 * s0_10[k]
                          + f_0 * s1_21[k];
            }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, s0_11, s0_12, s0_13, s0_14, s1_22, \
                         s1_23, s1_24, s1_25, s1_26 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_16[k] = -4.0 * s0_11[k]
                          + f_0 * s1_22[k];

                t_17[k] = -3.0 * s0_12[k]
                          + f_0 * s1_23[k];

                t_18[k] = -2.0 * s0_13[k]
                          + f_0 * s1_24[k];

                t_19[k] = -s0_14[k]
                          + f_0 * s1_25[k];

                t_20[k] = f_0 * s1_26[k];
            }
        }
    }
}

auto
geom_h_z(CSimdMatrix &buffer, const size_t target, const size_t s0, const size_t s1,
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
            auto *t_0 = buffer.data(target + (r * 21 + 0) * ncomps + c);
            auto *t_1 = buffer.data(target + (r * 21 + 1) * ncomps + c);
            auto *t_2 = buffer.data(target + (r * 21 + 2) * ncomps + c);
            auto *t_3 = buffer.data(target + (r * 21 + 3) * ncomps + c);
            auto *t_4 = buffer.data(target + (r * 21 + 4) * ncomps + c);
            auto *t_5 = buffer.data(target + (r * 21 + 5) * ncomps + c);
            auto *t_6 = buffer.data(target + (r * 21 + 6) * ncomps + c);
            auto *t_7 = buffer.data(target + (r * 21 + 7) * ncomps + c);
            auto *t_8 = buffer.data(target + (r * 21 + 8) * ncomps + c);
            auto *t_9 = buffer.data(target + (r * 21 + 9) * ncomps + c);
            auto *t_10 = buffer.data(target + (r * 21 + 10) * ncomps + c);
            auto *t_11 = buffer.data(target + (r * 21 + 11) * ncomps + c);
            auto *t_12 = buffer.data(target + (r * 21 + 12) * ncomps + c);
            auto *t_13 = buffer.data(target + (r * 21 + 13) * ncomps + c);
            auto *t_14 = buffer.data(target + (r * 21 + 14) * ncomps + c);
            auto *t_15 = buffer.data(target + (r * 21 + 15) * ncomps + c);
            auto *t_16 = buffer.data(target + (r * 21 + 16) * ncomps + c);
            auto *t_17 = buffer.data(target + (r * 21 + 17) * ncomps + c);
            auto *t_18 = buffer.data(target + (r * 21 + 18) * ncomps + c);
            auto *t_19 = buffer.data(target + (r * 21 + 19) * ncomps + c);
            auto *t_20 = buffer.data(target + (r * 21 + 20) * ncomps + c);
            const auto *s0_0 = buffer.data(s0 + (r * 15 + 0) * ncomps + c);
            const auto *s0_1 = buffer.data(s0 + (r * 15 + 1) * ncomps + c);
            const auto *s0_2 = buffer.data(s0 + (r * 15 + 2) * ncomps + c);
            const auto *s0_3 = buffer.data(s0 + (r * 15 + 3) * ncomps + c);
            const auto *s0_4 = buffer.data(s0 + (r * 15 + 4) * ncomps + c);
            const auto *s0_5 = buffer.data(s0 + (r * 15 + 5) * ncomps + c);
            const auto *s0_6 = buffer.data(s0 + (r * 15 + 6) * ncomps + c);
            const auto *s0_7 = buffer.data(s0 + (r * 15 + 7) * ncomps + c);
            const auto *s0_8 = buffer.data(s0 + (r * 15 + 8) * ncomps + c);
            const auto *s0_9 = buffer.data(s0 + (r * 15 + 9) * ncomps + c);
            const auto *s0_10 = buffer.data(s0 + (r * 15 + 10) * ncomps + c);
            const auto *s0_11 = buffer.data(s0 + (r * 15 + 11) * ncomps + c);
            const auto *s0_12 = buffer.data(s0 + (r * 15 + 12) * ncomps + c);
            const auto *s0_13 = buffer.data(s0 + (r * 15 + 13) * ncomps + c);
            const auto *s0_14 = buffer.data(s0 + (r * 15 + 14) * ncomps + c);
            const auto *s1_2 = buffer.data(s1 + (r * 28 + 2) * ncomps + c);
            const auto *s1_4 = buffer.data(s1 + (r * 28 + 4) * ncomps + c);
            const auto *s1_5 = buffer.data(s1 + (r * 28 + 5) * ncomps + c);
            const auto *s1_7 = buffer.data(s1 + (r * 28 + 7) * ncomps + c);
            const auto *s1_8 = buffer.data(s1 + (r * 28 + 8) * ncomps + c);
            const auto *s1_9 = buffer.data(s1 + (r * 28 + 9) * ncomps + c);
            const auto *s1_11 = buffer.data(s1 + (r * 28 + 11) * ncomps + c);
            const auto *s1_12 = buffer.data(s1 + (r * 28 + 12) * ncomps + c);
            const auto *s1_13 = buffer.data(s1 + (r * 28 + 13) * ncomps + c);
            const auto *s1_14 = buffer.data(s1 + (r * 28 + 14) * ncomps + c);
            const auto *s1_16 = buffer.data(s1 + (r * 28 + 16) * ncomps + c);
            const auto *s1_17 = buffer.data(s1 + (r * 28 + 17) * ncomps + c);
            const auto *s1_18 = buffer.data(s1 + (r * 28 + 18) * ncomps + c);
            const auto *s1_19 = buffer.data(s1 + (r * 28 + 19) * ncomps + c);
            const auto *s1_20 = buffer.data(s1 + (r * 28 + 20) * ncomps + c);
            const auto *s1_22 = buffer.data(s1 + (r * 28 + 22) * ncomps + c);
            const auto *s1_23 = buffer.data(s1 + (r * 28 + 23) * ncomps + c);
            const auto *s1_24 = buffer.data(s1 + (r * 28 + 24) * ncomps + c);
            const auto *s1_25 = buffer.data(s1 + (r * 28 + 25) * ncomps + c);
            const auto *s1_26 = buffer.data(s1 + (r * 28 + 26) * ncomps + c);
            const auto *s1_27 = buffer.data(s1 + (r * 28 + 27) * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, s0_0, s0_1, s0_2, s1_2, s1_4, s1_5, \
                         s1_7, s1_8, s1_9 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_0[k] = f_0 * s1_2[k];

                t_1[k] = f_0 * s1_4[k];

                t_2[k] = -s0_0[k]
                         + f_0 * s1_5[k];

                t_3[k] = f_0 * s1_7[k];

                t_4[k] = -s0_1[k]
                         + f_0 * s1_8[k];

                t_5[k] = -2.0 * s0_2[k]
                         + f_0 * s1_9[k];
            }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, t_11, s0_3, s0_4, s0_5, s0_6, s1_11, s1_12, \
                         s1_13, s1_14, s1_16, s1_17 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_6[k] = f_0 * s1_11[k];

                t_7[k] = -s0_3[k]
                         + f_0 * s1_12[k];

                t_8[k] = -2.0 * s0_4[k]
                         + f_0 * s1_13[k];

                t_9[k] = -3.0 * s0_5[k]
                         + f_0 * s1_14[k];

                t_10[k] = f_0 * s1_16[k];

                t_11[k] = -s0_6[k]
                          + f_0 * s1_17[k];
            }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, s0_7, s0_8, s0_9, s0_10, s1_18, s1_19, \
                         s1_20, s1_22, s1_23 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_12[k] = -2.0 * s0_7[k]
                          + f_0 * s1_18[k];

                t_13[k] = -3.0 * s0_8[k]
                          + f_0 * s1_19[k];

                t_14[k] = -4.0 * s0_9[k]
                          + f_0 * s1_20[k];

                t_15[k] = f_0 * s1_22[k];

                t_16[k] = -s0_10[k]
                          + f_0 * s1_23[k];
            }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, s0_11, s0_12, s0_13, s0_14, s1_24, s1_25, \
                         s1_26, s1_27 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_17[k] = -2.0 * s0_11[k]
                          + f_0 * s1_24[k];

                t_18[k] = -3.0 * s0_12[k]
                          + f_0 * s1_25[k];

                t_19[k] = -4.0 * s0_13[k]
                          + f_0 * s1_26[k];

                t_20[k] = -5.0 * s0_14[k]
                          + f_0 * s1_27[k];
            }
        }
    }
}

}  // namespace simdgeo
