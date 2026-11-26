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


#include "SimdGeometryL1.hpp"

#include "SimdAlign.hpp"

namespace simdgeo {  // simdgeo namespace

auto
geom_l_x(CSimdMatrix &buffer, const size_t target, const size_t s0, const size_t s1,
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
            auto *t_0 = buffer.data(target + (r * 45 + 0) * ncomps + c);
            auto *t_1 = buffer.data(target + (r * 45 + 1) * ncomps + c);
            auto *t_2 = buffer.data(target + (r * 45 + 2) * ncomps + c);
            auto *t_3 = buffer.data(target + (r * 45 + 3) * ncomps + c);
            auto *t_4 = buffer.data(target + (r * 45 + 4) * ncomps + c);
            auto *t_5 = buffer.data(target + (r * 45 + 5) * ncomps + c);
            auto *t_6 = buffer.data(target + (r * 45 + 6) * ncomps + c);
            auto *t_7 = buffer.data(target + (r * 45 + 7) * ncomps + c);
            auto *t_8 = buffer.data(target + (r * 45 + 8) * ncomps + c);
            auto *t_9 = buffer.data(target + (r * 45 + 9) * ncomps + c);
            auto *t_10 = buffer.data(target + (r * 45 + 10) * ncomps + c);
            auto *t_11 = buffer.data(target + (r * 45 + 11) * ncomps + c);
            auto *t_12 = buffer.data(target + (r * 45 + 12) * ncomps + c);
            auto *t_13 = buffer.data(target + (r * 45 + 13) * ncomps + c);
            auto *t_14 = buffer.data(target + (r * 45 + 14) * ncomps + c);
            auto *t_15 = buffer.data(target + (r * 45 + 15) * ncomps + c);
            auto *t_16 = buffer.data(target + (r * 45 + 16) * ncomps + c);
            auto *t_17 = buffer.data(target + (r * 45 + 17) * ncomps + c);
            auto *t_18 = buffer.data(target + (r * 45 + 18) * ncomps + c);
            auto *t_19 = buffer.data(target + (r * 45 + 19) * ncomps + c);
            auto *t_20 = buffer.data(target + (r * 45 + 20) * ncomps + c);
            auto *t_21 = buffer.data(target + (r * 45 + 21) * ncomps + c);
            auto *t_22 = buffer.data(target + (r * 45 + 22) * ncomps + c);
            auto *t_23 = buffer.data(target + (r * 45 + 23) * ncomps + c);
            auto *t_24 = buffer.data(target + (r * 45 + 24) * ncomps + c);
            auto *t_25 = buffer.data(target + (r * 45 + 25) * ncomps + c);
            auto *t_26 = buffer.data(target + (r * 45 + 26) * ncomps + c);
            auto *t_27 = buffer.data(target + (r * 45 + 27) * ncomps + c);
            auto *t_28 = buffer.data(target + (r * 45 + 28) * ncomps + c);
            auto *t_29 = buffer.data(target + (r * 45 + 29) * ncomps + c);
            auto *t_30 = buffer.data(target + (r * 45 + 30) * ncomps + c);
            auto *t_31 = buffer.data(target + (r * 45 + 31) * ncomps + c);
            auto *t_32 = buffer.data(target + (r * 45 + 32) * ncomps + c);
            auto *t_33 = buffer.data(target + (r * 45 + 33) * ncomps + c);
            auto *t_34 = buffer.data(target + (r * 45 + 34) * ncomps + c);
            auto *t_35 = buffer.data(target + (r * 45 + 35) * ncomps + c);
            auto *t_36 = buffer.data(target + (r * 45 + 36) * ncomps + c);
            auto *t_37 = buffer.data(target + (r * 45 + 37) * ncomps + c);
            auto *t_38 = buffer.data(target + (r * 45 + 38) * ncomps + c);
            auto *t_39 = buffer.data(target + (r * 45 + 39) * ncomps + c);
            auto *t_40 = buffer.data(target + (r * 45 + 40) * ncomps + c);
            auto *t_41 = buffer.data(target + (r * 45 + 41) * ncomps + c);
            auto *t_42 = buffer.data(target + (r * 45 + 42) * ncomps + c);
            auto *t_43 = buffer.data(target + (r * 45 + 43) * ncomps + c);
            auto *t_44 = buffer.data(target + (r * 45 + 44) * ncomps + c);
            const auto *s0_0 = buffer.data(s0 + (r * 36 + 0) * ncomps + c);
            const auto *s0_1 = buffer.data(s0 + (r * 36 + 1) * ncomps + c);
            const auto *s0_2 = buffer.data(s0 + (r * 36 + 2) * ncomps + c);
            const auto *s0_3 = buffer.data(s0 + (r * 36 + 3) * ncomps + c);
            const auto *s0_4 = buffer.data(s0 + (r * 36 + 4) * ncomps + c);
            const auto *s0_5 = buffer.data(s0 + (r * 36 + 5) * ncomps + c);
            const auto *s0_6 = buffer.data(s0 + (r * 36 + 6) * ncomps + c);
            const auto *s0_7 = buffer.data(s0 + (r * 36 + 7) * ncomps + c);
            const auto *s0_8 = buffer.data(s0 + (r * 36 + 8) * ncomps + c);
            const auto *s0_9 = buffer.data(s0 + (r * 36 + 9) * ncomps + c);
            const auto *s0_10 = buffer.data(s0 + (r * 36 + 10) * ncomps + c);
            const auto *s0_11 = buffer.data(s0 + (r * 36 + 11) * ncomps + c);
            const auto *s0_12 = buffer.data(s0 + (r * 36 + 12) * ncomps + c);
            const auto *s0_13 = buffer.data(s0 + (r * 36 + 13) * ncomps + c);
            const auto *s0_14 = buffer.data(s0 + (r * 36 + 14) * ncomps + c);
            const auto *s0_15 = buffer.data(s0 + (r * 36 + 15) * ncomps + c);
            const auto *s0_16 = buffer.data(s0 + (r * 36 + 16) * ncomps + c);
            const auto *s0_17 = buffer.data(s0 + (r * 36 + 17) * ncomps + c);
            const auto *s0_18 = buffer.data(s0 + (r * 36 + 18) * ncomps + c);
            const auto *s0_19 = buffer.data(s0 + (r * 36 + 19) * ncomps + c);
            const auto *s0_20 = buffer.data(s0 + (r * 36 + 20) * ncomps + c);
            const auto *s0_21 = buffer.data(s0 + (r * 36 + 21) * ncomps + c);
            const auto *s0_22 = buffer.data(s0 + (r * 36 + 22) * ncomps + c);
            const auto *s0_23 = buffer.data(s0 + (r * 36 + 23) * ncomps + c);
            const auto *s0_24 = buffer.data(s0 + (r * 36 + 24) * ncomps + c);
            const auto *s0_25 = buffer.data(s0 + (r * 36 + 25) * ncomps + c);
            const auto *s0_26 = buffer.data(s0 + (r * 36 + 26) * ncomps + c);
            const auto *s0_27 = buffer.data(s0 + (r * 36 + 27) * ncomps + c);
            const auto *s0_28 = buffer.data(s0 + (r * 36 + 28) * ncomps + c);
            const auto *s0_29 = buffer.data(s0 + (r * 36 + 29) * ncomps + c);
            const auto *s0_30 = buffer.data(s0 + (r * 36 + 30) * ncomps + c);
            const auto *s0_31 = buffer.data(s0 + (r * 36 + 31) * ncomps + c);
            const auto *s0_32 = buffer.data(s0 + (r * 36 + 32) * ncomps + c);
            const auto *s0_33 = buffer.data(s0 + (r * 36 + 33) * ncomps + c);
            const auto *s0_34 = buffer.data(s0 + (r * 36 + 34) * ncomps + c);
            const auto *s0_35 = buffer.data(s0 + (r * 36 + 35) * ncomps + c);
            const auto *s1_0 = buffer.data(s1 + (r * 55 + 0) * ncomps + c);
            const auto *s1_1 = buffer.data(s1 + (r * 55 + 1) * ncomps + c);
            const auto *s1_2 = buffer.data(s1 + (r * 55 + 2) * ncomps + c);
            const auto *s1_3 = buffer.data(s1 + (r * 55 + 3) * ncomps + c);
            const auto *s1_4 = buffer.data(s1 + (r * 55 + 4) * ncomps + c);
            const auto *s1_5 = buffer.data(s1 + (r * 55 + 5) * ncomps + c);
            const auto *s1_6 = buffer.data(s1 + (r * 55 + 6) * ncomps + c);
            const auto *s1_7 = buffer.data(s1 + (r * 55 + 7) * ncomps + c);
            const auto *s1_8 = buffer.data(s1 + (r * 55 + 8) * ncomps + c);
            const auto *s1_9 = buffer.data(s1 + (r * 55 + 9) * ncomps + c);
            const auto *s1_10 = buffer.data(s1 + (r * 55 + 10) * ncomps + c);
            const auto *s1_11 = buffer.data(s1 + (r * 55 + 11) * ncomps + c);
            const auto *s1_12 = buffer.data(s1 + (r * 55 + 12) * ncomps + c);
            const auto *s1_13 = buffer.data(s1 + (r * 55 + 13) * ncomps + c);
            const auto *s1_14 = buffer.data(s1 + (r * 55 + 14) * ncomps + c);
            const auto *s1_15 = buffer.data(s1 + (r * 55 + 15) * ncomps + c);
            const auto *s1_16 = buffer.data(s1 + (r * 55 + 16) * ncomps + c);
            const auto *s1_17 = buffer.data(s1 + (r * 55 + 17) * ncomps + c);
            const auto *s1_18 = buffer.data(s1 + (r * 55 + 18) * ncomps + c);
            const auto *s1_19 = buffer.data(s1 + (r * 55 + 19) * ncomps + c);
            const auto *s1_20 = buffer.data(s1 + (r * 55 + 20) * ncomps + c);
            const auto *s1_21 = buffer.data(s1 + (r * 55 + 21) * ncomps + c);
            const auto *s1_22 = buffer.data(s1 + (r * 55 + 22) * ncomps + c);
            const auto *s1_23 = buffer.data(s1 + (r * 55 + 23) * ncomps + c);
            const auto *s1_24 = buffer.data(s1 + (r * 55 + 24) * ncomps + c);
            const auto *s1_25 = buffer.data(s1 + (r * 55 + 25) * ncomps + c);
            const auto *s1_26 = buffer.data(s1 + (r * 55 + 26) * ncomps + c);
            const auto *s1_27 = buffer.data(s1 + (r * 55 + 27) * ncomps + c);
            const auto *s1_28 = buffer.data(s1 + (r * 55 + 28) * ncomps + c);
            const auto *s1_29 = buffer.data(s1 + (r * 55 + 29) * ncomps + c);
            const auto *s1_30 = buffer.data(s1 + (r * 55 + 30) * ncomps + c);
            const auto *s1_31 = buffer.data(s1 + (r * 55 + 31) * ncomps + c);
            const auto *s1_32 = buffer.data(s1 + (r * 55 + 32) * ncomps + c);
            const auto *s1_33 = buffer.data(s1 + (r * 55 + 33) * ncomps + c);
            const auto *s1_34 = buffer.data(s1 + (r * 55 + 34) * ncomps + c);
            const auto *s1_35 = buffer.data(s1 + (r * 55 + 35) * ncomps + c);
            const auto *s1_36 = buffer.data(s1 + (r * 55 + 36) * ncomps + c);
            const auto *s1_37 = buffer.data(s1 + (r * 55 + 37) * ncomps + c);
            const auto *s1_38 = buffer.data(s1 + (r * 55 + 38) * ncomps + c);
            const auto *s1_39 = buffer.data(s1 + (r * 55 + 39) * ncomps + c);
            const auto *s1_40 = buffer.data(s1 + (r * 55 + 40) * ncomps + c);
            const auto *s1_41 = buffer.data(s1 + (r * 55 + 41) * ncomps + c);
            const auto *s1_42 = buffer.data(s1 + (r * 55 + 42) * ncomps + c);
            const auto *s1_43 = buffer.data(s1 + (r * 55 + 43) * ncomps + c);
            const auto *s1_44 = buffer.data(s1 + (r * 55 + 44) * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, s0_0, s0_1, s0_2, s0_3, s0_4, s1_0, s1_1, \
                         s1_2, s1_3, s1_4 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_0[k] = -8.0 * s0_0[k]
                         + f_0 * s1_0[k];

                t_1[k] = -7.0 * s0_1[k]
                         + f_0 * s1_1[k];

                t_2[k] = -7.0 * s0_2[k]
                         + f_0 * s1_2[k];

                t_3[k] = -6.0 * s0_3[k]
                         + f_0 * s1_3[k];

                t_4[k] = -6.0 * s0_4[k]
                         + f_0 * s1_4[k];
            }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, s0_5, s0_6, s0_7, s0_8, s0_9, s1_5, s1_6, \
                         s1_7, s1_8, s1_9 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_5[k] = -6.0 * s0_5[k]
                         + f_0 * s1_5[k];

                t_6[k] = -5.0 * s0_6[k]
                         + f_0 * s1_6[k];

                t_7[k] = -5.0 * s0_7[k]
                         + f_0 * s1_7[k];

                t_8[k] = -5.0 * s0_8[k]
                         + f_0 * s1_8[k];

                t_9[k] = -5.0 * s0_9[k]
                         + f_0 * s1_9[k];
            }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, s0_10, s0_11, s0_12, s0_13, s0_14, \
                         s1_10, s1_11, s1_12, s1_13, s1_14 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_10[k] = -4.0 * s0_10[k]
                          + f_0 * s1_10[k];

                t_11[k] = -4.0 * s0_11[k]
                          + f_0 * s1_11[k];

                t_12[k] = -4.0 * s0_12[k]
                          + f_0 * s1_12[k];

                t_13[k] = -4.0 * s0_13[k]
                          + f_0 * s1_13[k];

                t_14[k] = -4.0 * s0_14[k]
                          + f_0 * s1_14[k];
            }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, s0_15, s0_16, s0_17, s0_18, s0_19, \
                         s1_15, s1_16, s1_17, s1_18, s1_19 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_15[k] = -3.0 * s0_15[k]
                          + f_0 * s1_15[k];

                t_16[k] = -3.0 * s0_16[k]
                          + f_0 * s1_16[k];

                t_17[k] = -3.0 * s0_17[k]
                          + f_0 * s1_17[k];

                t_18[k] = -3.0 * s0_18[k]
                          + f_0 * s1_18[k];

                t_19[k] = -3.0 * s0_19[k]
                          + f_0 * s1_19[k];
            }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, s0_20, s0_21, s0_22, s0_23, s0_24, \
                         s1_20, s1_21, s1_22, s1_23, s1_24 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_20[k] = -3.0 * s0_20[k]
                          + f_0 * s1_20[k];

                t_21[k] = -2.0 * s0_21[k]
                          + f_0 * s1_21[k];

                t_22[k] = -2.0 * s0_22[k]
                          + f_0 * s1_22[k];

                t_23[k] = -2.0 * s0_23[k]
                          + f_0 * s1_23[k];

                t_24[k] = -2.0 * s0_24[k]
                          + f_0 * s1_24[k];
            }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, s0_25, s0_26, s0_27, s0_28, s0_29, \
                         s1_25, s1_26, s1_27, s1_28, s1_29 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_25[k] = -2.0 * s0_25[k]
                          + f_0 * s1_25[k];

                t_26[k] = -2.0 * s0_26[k]
                          + f_0 * s1_26[k];

                t_27[k] = -2.0 * s0_27[k]
                          + f_0 * s1_27[k];

                t_28[k] = -s0_28[k]
                          + f_0 * s1_28[k];

                t_29[k] = -s0_29[k]
                          + f_0 * s1_29[k];
            }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, s0_30, s0_31, s0_32, s0_33, s0_34, \
                         s1_30, s1_31, s1_32, s1_33, s1_34 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_30[k] = -s0_30[k]
                          + f_0 * s1_30[k];

                t_31[k] = -s0_31[k]
                          + f_0 * s1_31[k];

                t_32[k] = -s0_32[k]
                          + f_0 * s1_32[k];

                t_33[k] = -s0_33[k]
                          + f_0 * s1_33[k];

                t_34[k] = -s0_34[k]
                          + f_0 * s1_34[k];
            }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, t_40, t_41, s0_35, s1_35, s1_36, s1_37, \
                         s1_38, s1_39, s1_40, s1_41 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_35[k] = -s0_35[k]
                          + f_0 * s1_35[k];

                t_36[k] = f_0 * s1_36[k];

                t_37[k] = f_0 * s1_37[k];

                t_38[k] = f_0 * s1_38[k];

                t_39[k] = f_0 * s1_39[k];

                t_40[k] = f_0 * s1_40[k];

                t_41[k] = f_0 * s1_41[k];
            }

#pragma omp simd aligned(t_42, t_43, t_44, s1_42, s1_43, s1_44 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_42[k] = f_0 * s1_42[k];

                t_43[k] = f_0 * s1_43[k];

                t_44[k] = f_0 * s1_44[k];
            }
        }
    }
}

auto
geom_l_y(CSimdMatrix &buffer, const size_t target, const size_t s0, const size_t s1,
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
            auto *t_0 = buffer.data(target + (r * 45 + 0) * ncomps + c);
            auto *t_1 = buffer.data(target + (r * 45 + 1) * ncomps + c);
            auto *t_2 = buffer.data(target + (r * 45 + 2) * ncomps + c);
            auto *t_3 = buffer.data(target + (r * 45 + 3) * ncomps + c);
            auto *t_4 = buffer.data(target + (r * 45 + 4) * ncomps + c);
            auto *t_5 = buffer.data(target + (r * 45 + 5) * ncomps + c);
            auto *t_6 = buffer.data(target + (r * 45 + 6) * ncomps + c);
            auto *t_7 = buffer.data(target + (r * 45 + 7) * ncomps + c);
            auto *t_8 = buffer.data(target + (r * 45 + 8) * ncomps + c);
            auto *t_9 = buffer.data(target + (r * 45 + 9) * ncomps + c);
            auto *t_10 = buffer.data(target + (r * 45 + 10) * ncomps + c);
            auto *t_11 = buffer.data(target + (r * 45 + 11) * ncomps + c);
            auto *t_12 = buffer.data(target + (r * 45 + 12) * ncomps + c);
            auto *t_13 = buffer.data(target + (r * 45 + 13) * ncomps + c);
            auto *t_14 = buffer.data(target + (r * 45 + 14) * ncomps + c);
            auto *t_15 = buffer.data(target + (r * 45 + 15) * ncomps + c);
            auto *t_16 = buffer.data(target + (r * 45 + 16) * ncomps + c);
            auto *t_17 = buffer.data(target + (r * 45 + 17) * ncomps + c);
            auto *t_18 = buffer.data(target + (r * 45 + 18) * ncomps + c);
            auto *t_19 = buffer.data(target + (r * 45 + 19) * ncomps + c);
            auto *t_20 = buffer.data(target + (r * 45 + 20) * ncomps + c);
            auto *t_21 = buffer.data(target + (r * 45 + 21) * ncomps + c);
            auto *t_22 = buffer.data(target + (r * 45 + 22) * ncomps + c);
            auto *t_23 = buffer.data(target + (r * 45 + 23) * ncomps + c);
            auto *t_24 = buffer.data(target + (r * 45 + 24) * ncomps + c);
            auto *t_25 = buffer.data(target + (r * 45 + 25) * ncomps + c);
            auto *t_26 = buffer.data(target + (r * 45 + 26) * ncomps + c);
            auto *t_27 = buffer.data(target + (r * 45 + 27) * ncomps + c);
            auto *t_28 = buffer.data(target + (r * 45 + 28) * ncomps + c);
            auto *t_29 = buffer.data(target + (r * 45 + 29) * ncomps + c);
            auto *t_30 = buffer.data(target + (r * 45 + 30) * ncomps + c);
            auto *t_31 = buffer.data(target + (r * 45 + 31) * ncomps + c);
            auto *t_32 = buffer.data(target + (r * 45 + 32) * ncomps + c);
            auto *t_33 = buffer.data(target + (r * 45 + 33) * ncomps + c);
            auto *t_34 = buffer.data(target + (r * 45 + 34) * ncomps + c);
            auto *t_35 = buffer.data(target + (r * 45 + 35) * ncomps + c);
            auto *t_36 = buffer.data(target + (r * 45 + 36) * ncomps + c);
            auto *t_37 = buffer.data(target + (r * 45 + 37) * ncomps + c);
            auto *t_38 = buffer.data(target + (r * 45 + 38) * ncomps + c);
            auto *t_39 = buffer.data(target + (r * 45 + 39) * ncomps + c);
            auto *t_40 = buffer.data(target + (r * 45 + 40) * ncomps + c);
            auto *t_41 = buffer.data(target + (r * 45 + 41) * ncomps + c);
            auto *t_42 = buffer.data(target + (r * 45 + 42) * ncomps + c);
            auto *t_43 = buffer.data(target + (r * 45 + 43) * ncomps + c);
            auto *t_44 = buffer.data(target + (r * 45 + 44) * ncomps + c);
            const auto *s0_0 = buffer.data(s0 + (r * 36 + 0) * ncomps + c);
            const auto *s0_1 = buffer.data(s0 + (r * 36 + 1) * ncomps + c);
            const auto *s0_2 = buffer.data(s0 + (r * 36 + 2) * ncomps + c);
            const auto *s0_3 = buffer.data(s0 + (r * 36 + 3) * ncomps + c);
            const auto *s0_4 = buffer.data(s0 + (r * 36 + 4) * ncomps + c);
            const auto *s0_5 = buffer.data(s0 + (r * 36 + 5) * ncomps + c);
            const auto *s0_6 = buffer.data(s0 + (r * 36 + 6) * ncomps + c);
            const auto *s0_7 = buffer.data(s0 + (r * 36 + 7) * ncomps + c);
            const auto *s0_8 = buffer.data(s0 + (r * 36 + 8) * ncomps + c);
            const auto *s0_9 = buffer.data(s0 + (r * 36 + 9) * ncomps + c);
            const auto *s0_10 = buffer.data(s0 + (r * 36 + 10) * ncomps + c);
            const auto *s0_11 = buffer.data(s0 + (r * 36 + 11) * ncomps + c);
            const auto *s0_12 = buffer.data(s0 + (r * 36 + 12) * ncomps + c);
            const auto *s0_13 = buffer.data(s0 + (r * 36 + 13) * ncomps + c);
            const auto *s0_14 = buffer.data(s0 + (r * 36 + 14) * ncomps + c);
            const auto *s0_15 = buffer.data(s0 + (r * 36 + 15) * ncomps + c);
            const auto *s0_16 = buffer.data(s0 + (r * 36 + 16) * ncomps + c);
            const auto *s0_17 = buffer.data(s0 + (r * 36 + 17) * ncomps + c);
            const auto *s0_18 = buffer.data(s0 + (r * 36 + 18) * ncomps + c);
            const auto *s0_19 = buffer.data(s0 + (r * 36 + 19) * ncomps + c);
            const auto *s0_20 = buffer.data(s0 + (r * 36 + 20) * ncomps + c);
            const auto *s0_21 = buffer.data(s0 + (r * 36 + 21) * ncomps + c);
            const auto *s0_22 = buffer.data(s0 + (r * 36 + 22) * ncomps + c);
            const auto *s0_23 = buffer.data(s0 + (r * 36 + 23) * ncomps + c);
            const auto *s0_24 = buffer.data(s0 + (r * 36 + 24) * ncomps + c);
            const auto *s0_25 = buffer.data(s0 + (r * 36 + 25) * ncomps + c);
            const auto *s0_26 = buffer.data(s0 + (r * 36 + 26) * ncomps + c);
            const auto *s0_27 = buffer.data(s0 + (r * 36 + 27) * ncomps + c);
            const auto *s0_28 = buffer.data(s0 + (r * 36 + 28) * ncomps + c);
            const auto *s0_29 = buffer.data(s0 + (r * 36 + 29) * ncomps + c);
            const auto *s0_30 = buffer.data(s0 + (r * 36 + 30) * ncomps + c);
            const auto *s0_31 = buffer.data(s0 + (r * 36 + 31) * ncomps + c);
            const auto *s0_32 = buffer.data(s0 + (r * 36 + 32) * ncomps + c);
            const auto *s0_33 = buffer.data(s0 + (r * 36 + 33) * ncomps + c);
            const auto *s0_34 = buffer.data(s0 + (r * 36 + 34) * ncomps + c);
            const auto *s0_35 = buffer.data(s0 + (r * 36 + 35) * ncomps + c);
            const auto *s1_1 = buffer.data(s1 + (r * 55 + 1) * ncomps + c);
            const auto *s1_3 = buffer.data(s1 + (r * 55 + 3) * ncomps + c);
            const auto *s1_4 = buffer.data(s1 + (r * 55 + 4) * ncomps + c);
            const auto *s1_6 = buffer.data(s1 + (r * 55 + 6) * ncomps + c);
            const auto *s1_7 = buffer.data(s1 + (r * 55 + 7) * ncomps + c);
            const auto *s1_8 = buffer.data(s1 + (r * 55 + 8) * ncomps + c);
            const auto *s1_10 = buffer.data(s1 + (r * 55 + 10) * ncomps + c);
            const auto *s1_11 = buffer.data(s1 + (r * 55 + 11) * ncomps + c);
            const auto *s1_12 = buffer.data(s1 + (r * 55 + 12) * ncomps + c);
            const auto *s1_13 = buffer.data(s1 + (r * 55 + 13) * ncomps + c);
            const auto *s1_15 = buffer.data(s1 + (r * 55 + 15) * ncomps + c);
            const auto *s1_16 = buffer.data(s1 + (r * 55 + 16) * ncomps + c);
            const auto *s1_17 = buffer.data(s1 + (r * 55 + 17) * ncomps + c);
            const auto *s1_18 = buffer.data(s1 + (r * 55 + 18) * ncomps + c);
            const auto *s1_19 = buffer.data(s1 + (r * 55 + 19) * ncomps + c);
            const auto *s1_21 = buffer.data(s1 + (r * 55 + 21) * ncomps + c);
            const auto *s1_22 = buffer.data(s1 + (r * 55 + 22) * ncomps + c);
            const auto *s1_23 = buffer.data(s1 + (r * 55 + 23) * ncomps + c);
            const auto *s1_24 = buffer.data(s1 + (r * 55 + 24) * ncomps + c);
            const auto *s1_25 = buffer.data(s1 + (r * 55 + 25) * ncomps + c);
            const auto *s1_26 = buffer.data(s1 + (r * 55 + 26) * ncomps + c);
            const auto *s1_28 = buffer.data(s1 + (r * 55 + 28) * ncomps + c);
            const auto *s1_29 = buffer.data(s1 + (r * 55 + 29) * ncomps + c);
            const auto *s1_30 = buffer.data(s1 + (r * 55 + 30) * ncomps + c);
            const auto *s1_31 = buffer.data(s1 + (r * 55 + 31) * ncomps + c);
            const auto *s1_32 = buffer.data(s1 + (r * 55 + 32) * ncomps + c);
            const auto *s1_33 = buffer.data(s1 + (r * 55 + 33) * ncomps + c);
            const auto *s1_34 = buffer.data(s1 + (r * 55 + 34) * ncomps + c);
            const auto *s1_36 = buffer.data(s1 + (r * 55 + 36) * ncomps + c);
            const auto *s1_37 = buffer.data(s1 + (r * 55 + 37) * ncomps + c);
            const auto *s1_38 = buffer.data(s1 + (r * 55 + 38) * ncomps + c);
            const auto *s1_39 = buffer.data(s1 + (r * 55 + 39) * ncomps + c);
            const auto *s1_40 = buffer.data(s1 + (r * 55 + 40) * ncomps + c);
            const auto *s1_41 = buffer.data(s1 + (r * 55 + 41) * ncomps + c);
            const auto *s1_42 = buffer.data(s1 + (r * 55 + 42) * ncomps + c);
            const auto *s1_43 = buffer.data(s1 + (r * 55 + 43) * ncomps + c);
            const auto *s1_45 = buffer.data(s1 + (r * 55 + 45) * ncomps + c);
            const auto *s1_46 = buffer.data(s1 + (r * 55 + 46) * ncomps + c);
            const auto *s1_47 = buffer.data(s1 + (r * 55 + 47) * ncomps + c);
            const auto *s1_48 = buffer.data(s1 + (r * 55 + 48) * ncomps + c);
            const auto *s1_49 = buffer.data(s1 + (r * 55 + 49) * ncomps + c);
            const auto *s1_50 = buffer.data(s1 + (r * 55 + 50) * ncomps + c);
            const auto *s1_51 = buffer.data(s1 + (r * 55 + 51) * ncomps + c);
            const auto *s1_52 = buffer.data(s1 + (r * 55 + 52) * ncomps + c);
            const auto *s1_53 = buffer.data(s1 + (r * 55 + 53) * ncomps + c);

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

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, s0_15, s0_16, s0_17, s0_18, s0_19, \
                         s1_28, s1_29, s1_30, s1_31, s1_32 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_21[k] = -6.0 * s0_15[k]
                          + f_0 * s1_28[k];

                t_22[k] = -5.0 * s0_16[k]
                          + f_0 * s1_29[k];

                t_23[k] = -4.0 * s0_17[k]
                          + f_0 * s1_30[k];

                t_24[k] = -3.0 * s0_18[k]
                          + f_0 * s1_31[k];

                t_25[k] = -2.0 * s0_19[k]
                          + f_0 * s1_32[k];
            }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, s0_20, s0_21, s0_22, s0_23, s1_33, \
                         s1_34, s1_36, s1_37, s1_38 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_26[k] = -s0_20[k]
                          + f_0 * s1_33[k];

                t_27[k] = f_0 * s1_34[k];

                t_28[k] = -7.0 * s0_21[k]
                          + f_0 * s1_36[k];

                t_29[k] = -6.0 * s0_22[k]
                          + f_0 * s1_37[k];

                t_30[k] = -5.0 * s0_23[k]
                          + f_0 * s1_38[k];
            }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, s0_24, s0_25, s0_26, s0_27, s1_39, \
                         s1_40, s1_41, s1_42, s1_43 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_31[k] = -4.0 * s0_24[k]
                          + f_0 * s1_39[k];

                t_32[k] = -3.0 * s0_25[k]
                          + f_0 * s1_40[k];

                t_33[k] = -2.0 * s0_26[k]
                          + f_0 * s1_41[k];

                t_34[k] = -s0_27[k]
                          + f_0 * s1_42[k];

                t_35[k] = f_0 * s1_43[k];
            }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, s0_28, s0_29, s0_30, s0_31, s0_32, \
                         s1_45, s1_46, s1_47, s1_48, s1_49 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_36[k] = -8.0 * s0_28[k]
                          + f_0 * s1_45[k];

                t_37[k] = -7.0 * s0_29[k]
                          + f_0 * s1_46[k];

                t_38[k] = -6.0 * s0_30[k]
                          + f_0 * s1_47[k];

                t_39[k] = -5.0 * s0_31[k]
                          + f_0 * s1_48[k];

                t_40[k] = -4.0 * s0_32[k]
                          + f_0 * s1_49[k];
            }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, s0_33, s0_34, s0_35, s1_50, s1_51, s1_52, \
                         s1_53 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_41[k] = -3.0 * s0_33[k]
                          + f_0 * s1_50[k];

                t_42[k] = -2.0 * s0_34[k]
                          + f_0 * s1_51[k];

                t_43[k] = -s0_35[k]
                          + f_0 * s1_52[k];

                t_44[k] = f_0 * s1_53[k];
            }
        }
    }
}

auto
geom_l_z(CSimdMatrix &buffer, const size_t target, const size_t s0, const size_t s1,
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
            auto *t_0 = buffer.data(target + (r * 45 + 0) * ncomps + c);
            auto *t_1 = buffer.data(target + (r * 45 + 1) * ncomps + c);
            auto *t_2 = buffer.data(target + (r * 45 + 2) * ncomps + c);
            auto *t_3 = buffer.data(target + (r * 45 + 3) * ncomps + c);
            auto *t_4 = buffer.data(target + (r * 45 + 4) * ncomps + c);
            auto *t_5 = buffer.data(target + (r * 45 + 5) * ncomps + c);
            auto *t_6 = buffer.data(target + (r * 45 + 6) * ncomps + c);
            auto *t_7 = buffer.data(target + (r * 45 + 7) * ncomps + c);
            auto *t_8 = buffer.data(target + (r * 45 + 8) * ncomps + c);
            auto *t_9 = buffer.data(target + (r * 45 + 9) * ncomps + c);
            auto *t_10 = buffer.data(target + (r * 45 + 10) * ncomps + c);
            auto *t_11 = buffer.data(target + (r * 45 + 11) * ncomps + c);
            auto *t_12 = buffer.data(target + (r * 45 + 12) * ncomps + c);
            auto *t_13 = buffer.data(target + (r * 45 + 13) * ncomps + c);
            auto *t_14 = buffer.data(target + (r * 45 + 14) * ncomps + c);
            auto *t_15 = buffer.data(target + (r * 45 + 15) * ncomps + c);
            auto *t_16 = buffer.data(target + (r * 45 + 16) * ncomps + c);
            auto *t_17 = buffer.data(target + (r * 45 + 17) * ncomps + c);
            auto *t_18 = buffer.data(target + (r * 45 + 18) * ncomps + c);
            auto *t_19 = buffer.data(target + (r * 45 + 19) * ncomps + c);
            auto *t_20 = buffer.data(target + (r * 45 + 20) * ncomps + c);
            auto *t_21 = buffer.data(target + (r * 45 + 21) * ncomps + c);
            auto *t_22 = buffer.data(target + (r * 45 + 22) * ncomps + c);
            auto *t_23 = buffer.data(target + (r * 45 + 23) * ncomps + c);
            auto *t_24 = buffer.data(target + (r * 45 + 24) * ncomps + c);
            auto *t_25 = buffer.data(target + (r * 45 + 25) * ncomps + c);
            auto *t_26 = buffer.data(target + (r * 45 + 26) * ncomps + c);
            auto *t_27 = buffer.data(target + (r * 45 + 27) * ncomps + c);
            auto *t_28 = buffer.data(target + (r * 45 + 28) * ncomps + c);
            auto *t_29 = buffer.data(target + (r * 45 + 29) * ncomps + c);
            auto *t_30 = buffer.data(target + (r * 45 + 30) * ncomps + c);
            auto *t_31 = buffer.data(target + (r * 45 + 31) * ncomps + c);
            auto *t_32 = buffer.data(target + (r * 45 + 32) * ncomps + c);
            auto *t_33 = buffer.data(target + (r * 45 + 33) * ncomps + c);
            auto *t_34 = buffer.data(target + (r * 45 + 34) * ncomps + c);
            auto *t_35 = buffer.data(target + (r * 45 + 35) * ncomps + c);
            auto *t_36 = buffer.data(target + (r * 45 + 36) * ncomps + c);
            auto *t_37 = buffer.data(target + (r * 45 + 37) * ncomps + c);
            auto *t_38 = buffer.data(target + (r * 45 + 38) * ncomps + c);
            auto *t_39 = buffer.data(target + (r * 45 + 39) * ncomps + c);
            auto *t_40 = buffer.data(target + (r * 45 + 40) * ncomps + c);
            auto *t_41 = buffer.data(target + (r * 45 + 41) * ncomps + c);
            auto *t_42 = buffer.data(target + (r * 45 + 42) * ncomps + c);
            auto *t_43 = buffer.data(target + (r * 45 + 43) * ncomps + c);
            auto *t_44 = buffer.data(target + (r * 45 + 44) * ncomps + c);
            const auto *s0_0 = buffer.data(s0 + (r * 36 + 0) * ncomps + c);
            const auto *s0_1 = buffer.data(s0 + (r * 36 + 1) * ncomps + c);
            const auto *s0_2 = buffer.data(s0 + (r * 36 + 2) * ncomps + c);
            const auto *s0_3 = buffer.data(s0 + (r * 36 + 3) * ncomps + c);
            const auto *s0_4 = buffer.data(s0 + (r * 36 + 4) * ncomps + c);
            const auto *s0_5 = buffer.data(s0 + (r * 36 + 5) * ncomps + c);
            const auto *s0_6 = buffer.data(s0 + (r * 36 + 6) * ncomps + c);
            const auto *s0_7 = buffer.data(s0 + (r * 36 + 7) * ncomps + c);
            const auto *s0_8 = buffer.data(s0 + (r * 36 + 8) * ncomps + c);
            const auto *s0_9 = buffer.data(s0 + (r * 36 + 9) * ncomps + c);
            const auto *s0_10 = buffer.data(s0 + (r * 36 + 10) * ncomps + c);
            const auto *s0_11 = buffer.data(s0 + (r * 36 + 11) * ncomps + c);
            const auto *s0_12 = buffer.data(s0 + (r * 36 + 12) * ncomps + c);
            const auto *s0_13 = buffer.data(s0 + (r * 36 + 13) * ncomps + c);
            const auto *s0_14 = buffer.data(s0 + (r * 36 + 14) * ncomps + c);
            const auto *s0_15 = buffer.data(s0 + (r * 36 + 15) * ncomps + c);
            const auto *s0_16 = buffer.data(s0 + (r * 36 + 16) * ncomps + c);
            const auto *s0_17 = buffer.data(s0 + (r * 36 + 17) * ncomps + c);
            const auto *s0_18 = buffer.data(s0 + (r * 36 + 18) * ncomps + c);
            const auto *s0_19 = buffer.data(s0 + (r * 36 + 19) * ncomps + c);
            const auto *s0_20 = buffer.data(s0 + (r * 36 + 20) * ncomps + c);
            const auto *s0_21 = buffer.data(s0 + (r * 36 + 21) * ncomps + c);
            const auto *s0_22 = buffer.data(s0 + (r * 36 + 22) * ncomps + c);
            const auto *s0_23 = buffer.data(s0 + (r * 36 + 23) * ncomps + c);
            const auto *s0_24 = buffer.data(s0 + (r * 36 + 24) * ncomps + c);
            const auto *s0_25 = buffer.data(s0 + (r * 36 + 25) * ncomps + c);
            const auto *s0_26 = buffer.data(s0 + (r * 36 + 26) * ncomps + c);
            const auto *s0_27 = buffer.data(s0 + (r * 36 + 27) * ncomps + c);
            const auto *s0_28 = buffer.data(s0 + (r * 36 + 28) * ncomps + c);
            const auto *s0_29 = buffer.data(s0 + (r * 36 + 29) * ncomps + c);
            const auto *s0_30 = buffer.data(s0 + (r * 36 + 30) * ncomps + c);
            const auto *s0_31 = buffer.data(s0 + (r * 36 + 31) * ncomps + c);
            const auto *s0_32 = buffer.data(s0 + (r * 36 + 32) * ncomps + c);
            const auto *s0_33 = buffer.data(s0 + (r * 36 + 33) * ncomps + c);
            const auto *s0_34 = buffer.data(s0 + (r * 36 + 34) * ncomps + c);
            const auto *s0_35 = buffer.data(s0 + (r * 36 + 35) * ncomps + c);
            const auto *s1_2 = buffer.data(s1 + (r * 55 + 2) * ncomps + c);
            const auto *s1_4 = buffer.data(s1 + (r * 55 + 4) * ncomps + c);
            const auto *s1_5 = buffer.data(s1 + (r * 55 + 5) * ncomps + c);
            const auto *s1_7 = buffer.data(s1 + (r * 55 + 7) * ncomps + c);
            const auto *s1_8 = buffer.data(s1 + (r * 55 + 8) * ncomps + c);
            const auto *s1_9 = buffer.data(s1 + (r * 55 + 9) * ncomps + c);
            const auto *s1_11 = buffer.data(s1 + (r * 55 + 11) * ncomps + c);
            const auto *s1_12 = buffer.data(s1 + (r * 55 + 12) * ncomps + c);
            const auto *s1_13 = buffer.data(s1 + (r * 55 + 13) * ncomps + c);
            const auto *s1_14 = buffer.data(s1 + (r * 55 + 14) * ncomps + c);
            const auto *s1_16 = buffer.data(s1 + (r * 55 + 16) * ncomps + c);
            const auto *s1_17 = buffer.data(s1 + (r * 55 + 17) * ncomps + c);
            const auto *s1_18 = buffer.data(s1 + (r * 55 + 18) * ncomps + c);
            const auto *s1_19 = buffer.data(s1 + (r * 55 + 19) * ncomps + c);
            const auto *s1_20 = buffer.data(s1 + (r * 55 + 20) * ncomps + c);
            const auto *s1_22 = buffer.data(s1 + (r * 55 + 22) * ncomps + c);
            const auto *s1_23 = buffer.data(s1 + (r * 55 + 23) * ncomps + c);
            const auto *s1_24 = buffer.data(s1 + (r * 55 + 24) * ncomps + c);
            const auto *s1_25 = buffer.data(s1 + (r * 55 + 25) * ncomps + c);
            const auto *s1_26 = buffer.data(s1 + (r * 55 + 26) * ncomps + c);
            const auto *s1_27 = buffer.data(s1 + (r * 55 + 27) * ncomps + c);
            const auto *s1_29 = buffer.data(s1 + (r * 55 + 29) * ncomps + c);
            const auto *s1_30 = buffer.data(s1 + (r * 55 + 30) * ncomps + c);
            const auto *s1_31 = buffer.data(s1 + (r * 55 + 31) * ncomps + c);
            const auto *s1_32 = buffer.data(s1 + (r * 55 + 32) * ncomps + c);
            const auto *s1_33 = buffer.data(s1 + (r * 55 + 33) * ncomps + c);
            const auto *s1_34 = buffer.data(s1 + (r * 55 + 34) * ncomps + c);
            const auto *s1_35 = buffer.data(s1 + (r * 55 + 35) * ncomps + c);
            const auto *s1_37 = buffer.data(s1 + (r * 55 + 37) * ncomps + c);
            const auto *s1_38 = buffer.data(s1 + (r * 55 + 38) * ncomps + c);
            const auto *s1_39 = buffer.data(s1 + (r * 55 + 39) * ncomps + c);
            const auto *s1_40 = buffer.data(s1 + (r * 55 + 40) * ncomps + c);
            const auto *s1_41 = buffer.data(s1 + (r * 55 + 41) * ncomps + c);
            const auto *s1_42 = buffer.data(s1 + (r * 55 + 42) * ncomps + c);
            const auto *s1_43 = buffer.data(s1 + (r * 55 + 43) * ncomps + c);
            const auto *s1_44 = buffer.data(s1 + (r * 55 + 44) * ncomps + c);
            const auto *s1_46 = buffer.data(s1 + (r * 55 + 46) * ncomps + c);
            const auto *s1_47 = buffer.data(s1 + (r * 55 + 47) * ncomps + c);
            const auto *s1_48 = buffer.data(s1 + (r * 55 + 48) * ncomps + c);
            const auto *s1_49 = buffer.data(s1 + (r * 55 + 49) * ncomps + c);
            const auto *s1_50 = buffer.data(s1 + (r * 55 + 50) * ncomps + c);
            const auto *s1_51 = buffer.data(s1 + (r * 55 + 51) * ncomps + c);
            const auto *s1_52 = buffer.data(s1 + (r * 55 + 52) * ncomps + c);
            const auto *s1_53 = buffer.data(s1 + (r * 55 + 53) * ncomps + c);
            const auto *s1_54 = buffer.data(s1 + (r * 55 + 54) * ncomps + c);

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

#pragma omp simd aligned(t_17, t_18, t_19, t_20, t_21, s0_11, s0_12, s0_13, s0_14, s1_24, \
                         s1_25, s1_26, s1_27, s1_29 : simd::cache_line_size())
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

                t_21[k] = f_0 * s1_29[k];
            }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, s0_15, s0_16, s0_17, s0_18, s0_19, \
                         s1_30, s1_31, s1_32, s1_33, s1_34 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_22[k] = -s0_15[k]
                          + f_0 * s1_30[k];

                t_23[k] = -2.0 * s0_16[k]
                          + f_0 * s1_31[k];

                t_24[k] = -3.0 * s0_17[k]
                          + f_0 * s1_32[k];

                t_25[k] = -4.0 * s0_18[k]
                          + f_0 * s1_33[k];

                t_26[k] = -5.0 * s0_19[k]
                          + f_0 * s1_34[k];
            }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, s0_20, s0_21, s0_22, s0_23, s1_35, \
                         s1_37, s1_38, s1_39, s1_40 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_27[k] = -6.0 * s0_20[k]
                          + f_0 * s1_35[k];

                t_28[k] = f_0 * s1_37[k];

                t_29[k] = -s0_21[k]
                          + f_0 * s1_38[k];

                t_30[k] = -2.0 * s0_22[k]
                          + f_0 * s1_39[k];

                t_31[k] = -3.0 * s0_23[k]
                          + f_0 * s1_40[k];
            }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, s0_24, s0_25, s0_26, s0_27, s1_41, \
                         s1_42, s1_43, s1_44, s1_46 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_32[k] = -4.0 * s0_24[k]
                          + f_0 * s1_41[k];

                t_33[k] = -5.0 * s0_25[k]
                          + f_0 * s1_42[k];

                t_34[k] = -6.0 * s0_26[k]
                          + f_0 * s1_43[k];

                t_35[k] = -7.0 * s0_27[k]
                          + f_0 * s1_44[k];

                t_36[k] = f_0 * s1_46[k];
            }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, t_41, s0_28, s0_29, s0_30, s0_31, s0_32, \
                         s1_47, s1_48, s1_49, s1_50, s1_51 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_37[k] = -s0_28[k]
                          + f_0 * s1_47[k];

                t_38[k] = -2.0 * s0_29[k]
                          + f_0 * s1_48[k];

                t_39[k] = -3.0 * s0_30[k]
                          + f_0 * s1_49[k];

                t_40[k] = -4.0 * s0_31[k]
                          + f_0 * s1_50[k];

                t_41[k] = -5.0 * s0_32[k]
                          + f_0 * s1_51[k];
            }

#pragma omp simd aligned(t_42, t_43, t_44, s0_33, s0_34, s0_35, s1_52, s1_53, \
                         s1_54 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_42[k] = -6.0 * s0_33[k]
                          + f_0 * s1_52[k];

                t_43[k] = -7.0 * s0_34[k]
                          + f_0 * s1_53[k];

                t_44[k] = -8.0 * s0_35[k]
                          + f_0 * s1_54[k];
            }
        }
    }
}

}  // namespace simdgeo
