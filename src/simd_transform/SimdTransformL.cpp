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


#include "SimdTransformL.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_l_inner(CSimdMatrix &buffer, const size_t target, const size_t source,
                  const size_t nrows, const size_t ncomps, const size_t ncols) -> void
{
    // NOTE: the factors are the shell's own, so they are formed once rather than
    // for every atom pair the shell pair reaches.

    const auto f_0 = 0.1875 * std::sqrt(715.0);
    const auto f_1 = 1.3125 * std::sqrt(715.0);
    const auto f_2 = 0.65625 * std::sqrt(715.0);
    const auto f_3 = 3.28125 * std::sqrt(715.0);
    const auto f_4 = 1.96875 * std::sqrt(715.0);
    const auto f_5 = 0.09375 * std::sqrt(715.0);
    const auto f_6 = 0.09375 * std::sqrt(858.0);
    const auto f_7 = 0.21875 * std::sqrt(858.0);
    const auto f_8 = 1.3125 * std::sqrt(858.0);
    const auto f_9 = 4.375 * std::sqrt(858.0);
    const auto f_10 = 0.46875 * std::sqrt(1001.0);
    const auto f_11 = 1.875 * std::sqrt(1001.0);
    const auto f_12 = 0.84375 * std::sqrt(1001.0);
    const auto f_13 = 3.75 * std::sqrt(1001.0);
    const auto f_14 = 0.09375 * std::sqrt(1001.0);
    const auto f_15 = 0.375 * std::sqrt(1001.0);
    const auto f_16 = 0.1875 * std::sqrt(77.0);
    const auto f_17 = 4.5 * std::sqrt(77.0);
    const auto f_18 = 7.5 * std::sqrt(77.0);
    const auto f_19 = 0.28125 * std::sqrt(1155.0);
    const auto f_20 = 0.46875 * std::sqrt(1155.0);
    const auto f_21 = 1.875 * std::sqrt(1155.0);
    const auto f_22 = 0.09375 * std::sqrt(1155.0);
    const auto f_23 = 1.25 * std::sqrt(1155.0);
    const auto f_24 = 1.5 * std::sqrt(1155.0);
    const auto f_25 = 0.625 * std::sqrt(1155.0);
    const auto f_26 = 0.5 * std::sqrt(1155.0);
    const auto f_27 = 0.09375 * std::sqrt(70.0);
    const auto f_28 = 0.28125 * std::sqrt(70.0);
    const auto f_29 = 2.8125 * std::sqrt(70.0);
    const auto f_30 = 5.625 * std::sqrt(70.0);
    const auto f_31 = 7.5 * std::sqrt(70.0);
    const auto f_32 = 3.0 * std::sqrt(70.0);
    const auto f_33 = 0.046875 * std::sqrt(70.0);
    const auto f_34 = 1.40625 * std::sqrt(70.0);
    const auto f_35 = 3.75 * std::sqrt(70.0);
    const auto f_36 = 1.5 * std::sqrt(70.0);
    const auto f_37 = 0.046875 * std::sqrt(77.0);
    const auto f_38 = 1.125 * std::sqrt(77.0);
    const auto f_39 = 0.46875 * std::sqrt(77.0);
    const auto f_40 = 5.625 * std::sqrt(77.0);
    const auto f_41 = 1.875 * std::sqrt(77.0);
    const auto f_42 = 11.25 * std::sqrt(77.0);
    const auto f_43 = 0.015625 * std::sqrt(858.0);
    const auto f_44 = 3.28125 * std::sqrt(858.0);
    const auto f_45 = 0.0234375 * std::sqrt(715.0);
    const auto f_46 = 1.640625 * std::sqrt(715.0);

    // NOTE: what sits either side of this index reaches the pass as a count and
    // nothing else -- the rows above it, the components below -- so one routine
    // per shell serves every block the shell appears in.

    for (size_t r = 0; r < nrows; r++)
    {
        for (size_t c = 0; c < ncomps; c++)
        {
            auto *t_0 = buffer.data(target + (r * 17 + 0) * ncomps + c);
            auto *t_1 = buffer.data(target + (r * 17 + 1) * ncomps + c);
            auto *t_2 = buffer.data(target + (r * 17 + 2) * ncomps + c);
            auto *t_3 = buffer.data(target + (r * 17 + 3) * ncomps + c);
            auto *t_4 = buffer.data(target + (r * 17 + 4) * ncomps + c);
            auto *t_5 = buffer.data(target + (r * 17 + 5) * ncomps + c);
            auto *t_6 = buffer.data(target + (r * 17 + 6) * ncomps + c);
            auto *t_7 = buffer.data(target + (r * 17 + 7) * ncomps + c);
            auto *t_8 = buffer.data(target + (r * 17 + 8) * ncomps + c);
            auto *t_9 = buffer.data(target + (r * 17 + 9) * ncomps + c);
            auto *t_10 = buffer.data(target + (r * 17 + 10) * ncomps + c);
            auto *t_11 = buffer.data(target + (r * 17 + 11) * ncomps + c);
            auto *t_12 = buffer.data(target + (r * 17 + 12) * ncomps + c);
            auto *t_13 = buffer.data(target + (r * 17 + 13) * ncomps + c);
            auto *t_14 = buffer.data(target + (r * 17 + 14) * ncomps + c);
            auto *t_15 = buffer.data(target + (r * 17 + 15) * ncomps + c);
            auto *t_16 = buffer.data(target + (r * 17 + 16) * ncomps + c);

            const auto *s_0 = buffer.data(source + (r * 45 + 0) * ncomps + c);
            const auto *s_1 = buffer.data(source + (r * 45 + 1) * ncomps + c);
            const auto *s_2 = buffer.data(source + (r * 45 + 2) * ncomps + c);
            const auto *s_3 = buffer.data(source + (r * 45 + 3) * ncomps + c);
            const auto *s_4 = buffer.data(source + (r * 45 + 4) * ncomps + c);
            const auto *s_5 = buffer.data(source + (r * 45 + 5) * ncomps + c);
            const auto *s_6 = buffer.data(source + (r * 45 + 6) * ncomps + c);
            const auto *s_7 = buffer.data(source + (r * 45 + 7) * ncomps + c);
            const auto *s_8 = buffer.data(source + (r * 45 + 8) * ncomps + c);
            const auto *s_9 = buffer.data(source + (r * 45 + 9) * ncomps + c);
            const auto *s_10 = buffer.data(source + (r * 45 + 10) * ncomps + c);
            const auto *s_11 = buffer.data(source + (r * 45 + 11) * ncomps + c);
            const auto *s_12 = buffer.data(source + (r * 45 + 12) * ncomps + c);
            const auto *s_13 = buffer.data(source + (r * 45 + 13) * ncomps + c);
            const auto *s_14 = buffer.data(source + (r * 45 + 14) * ncomps + c);
            const auto *s_15 = buffer.data(source + (r * 45 + 15) * ncomps + c);
            const auto *s_16 = buffer.data(source + (r * 45 + 16) * ncomps + c);
            const auto *s_17 = buffer.data(source + (r * 45 + 17) * ncomps + c);
            const auto *s_18 = buffer.data(source + (r * 45 + 18) * ncomps + c);
            const auto *s_19 = buffer.data(source + (r * 45 + 19) * ncomps + c);
            const auto *s_20 = buffer.data(source + (r * 45 + 20) * ncomps + c);
            const auto *s_21 = buffer.data(source + (r * 45 + 21) * ncomps + c);
            const auto *s_22 = buffer.data(source + (r * 45 + 22) * ncomps + c);
            const auto *s_23 = buffer.data(source + (r * 45 + 23) * ncomps + c);
            const auto *s_24 = buffer.data(source + (r * 45 + 24) * ncomps + c);
            const auto *s_25 = buffer.data(source + (r * 45 + 25) * ncomps + c);
            const auto *s_26 = buffer.data(source + (r * 45 + 26) * ncomps + c);
            const auto *s_27 = buffer.data(source + (r * 45 + 27) * ncomps + c);
            const auto *s_28 = buffer.data(source + (r * 45 + 28) * ncomps + c);
            const auto *s_29 = buffer.data(source + (r * 45 + 29) * ncomps + c);
            const auto *s_30 = buffer.data(source + (r * 45 + 30) * ncomps + c);
            const auto *s_31 = buffer.data(source + (r * 45 + 31) * ncomps + c);
            const auto *s_32 = buffer.data(source + (r * 45 + 32) * ncomps + c);
            const auto *s_33 = buffer.data(source + (r * 45 + 33) * ncomps + c);
            const auto *s_34 = buffer.data(source + (r * 45 + 34) * ncomps + c);
            const auto *s_35 = buffer.data(source + (r * 45 + 35) * ncomps + c);
            const auto *s_36 = buffer.data(source + (r * 45 + 36) * ncomps + c);
            const auto *s_37 = buffer.data(source + (r * 45 + 37) * ncomps + c);
            const auto *s_38 = buffer.data(source + (r * 45 + 38) * ncomps + c);
            const auto *s_39 = buffer.data(source + (r * 45 + 39) * ncomps + c);
            const auto *s_40 = buffer.data(source + (r * 45 + 40) * ncomps + c);
            const auto *s_41 = buffer.data(source + (r * 45 + 41) * ncomps + c);
            const auto *s_42 = buffer.data(source + (r * 45 + 42) * ncomps + c);
            const auto *s_43 = buffer.data(source + (r * 45 + 43) * ncomps + c);
            const auto *s_44 = buffer.data(source + (r * 45 + 44) * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, s_1, s_4, s_6, s_8, s_11, s_15, s_17, s_22, s_28, \
                         s_30, s_37 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_0[k] = f_0 * s_1[k]
                         - f_1 * s_6[k]
                         + f_1 * s_15[k]
                         - f_0 * s_28[k];

                t_1[k] = f_2 * s_4[k]
                         - f_3 * s_11[k]
                         + f_4 * s_22[k]
                         - f_5 * s_37[k];

                t_2[k] = -f_6 * s_1[k]
                         + f_7 * s_6[k]
                         + f_8 * s_8[k]
                         + f_7 * s_15[k]
                         - f_9 * s_17[k]
                         - f_6 * s_28[k]
                         + f_8 * s_30[k];
            }

#pragma omp simd aligned(t_3, s_4, s_11, s_13, s_22, s_24, s_37, s_39 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_3[k] = -f_10 * s_4[k]
                         + f_10 * s_11[k]
                         + f_11 * s_13[k]
                         + f_12 * s_22[k]
                         - f_13 * s_24[k]
                         - f_14 * s_37[k]
                         + f_15 * s_39[k];
            }

#pragma omp simd aligned(t_4, s_1, s_6, s_8, s_15, s_19, s_28, s_30, \
                         s_32 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_4[k] = f_16 * s_1[k]
                         + f_16 * s_6[k]
                         - f_17 * s_8[k]
                         - f_16 * s_15[k]
                         + f_18 * s_19[k]
                         - f_16 * s_28[k]
                         + f_17 * s_30[k]
                         - f_18 * s_32[k];
            }

#pragma omp simd aligned(t_5, s_4, s_11, s_13, s_22, s_24, s_26, s_37, s_39, \
                         s_41 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_5[k] = f_19 * s_4[k]
                         + f_20 * s_11[k]
                         - f_21 * s_13[k]
                         + f_22 * s_22[k]
                         - f_23 * s_24[k]
                         + f_24 * s_26[k]
                         - f_22 * s_37[k]
                         + f_25 * s_39[k]
                         - f_26 * s_41[k];
            }

#pragma omp simd aligned(t_6, s_1, s_6, s_8, s_15, s_17, s_19, s_28, s_30, s_32, \
                         s_34 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_6[k] = -f_27 * s_1[k]
                         - f_28 * s_6[k]
                         + f_29 * s_8[k]
                         - f_28 * s_15[k]
                         + f_30 * s_17[k]
                         - f_31 * s_19[k]
                         - f_27 * s_28[k]
                         + f_29 * s_30[k]
                         - f_31 * s_32[k]
                         + f_32 * s_34[k];
            }

#pragma omp simd aligned(t_7, s_4, s_11, s_13, s_22, s_24, s_26, s_37, s_39, s_41, \
                         s_43 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_7[k] = -3.28125 * s_4[k]
                         - 9.84375 * s_11[k]
                         + 26.25 * s_13[k]
                         - 9.84375 * s_22[k]
                         + 52.5 * s_24[k]
                         - 31.5 * s_26[k]
                         - 3.28125 * s_37[k]
                         + 26.25 * s_39[k]
                         - 31.5 * s_41[k]
                         + 6.0 * s_43[k];
            }

#pragma omp simd aligned(t_8, s_0, s_3, s_5, s_10, s_12, s_14, s_21, s_23, s_25, s_27, s_36, \
                         s_38, s_40, s_42, s_44 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_8[k] = 0.2734375 * s_0[k]
                         + 1.09375 * s_3[k]
                         - 8.75 * s_5[k]
                         + 1.640625 * s_10[k]
                         - 26.25 * s_12[k]
                         + 26.25 * s_14[k]
                         + 1.09375 * s_21[k]
                         - 26.25 * s_23[k]
                         + 52.5 * s_25[k]
                         - 14.0 * s_27[k]
                         + 0.2734375 * s_36[k]
                         - 8.75 * s_38[k]
                         + 26.25 * s_40[k]
                         - 14.0 * s_42[k]
                         + s_44[k];
            }

#pragma omp simd aligned(t_9, s_2, s_7, s_9, s_16, s_18, s_20, s_29, s_31, s_33, \
                         s_35 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_9[k] = -3.28125 * s_2[k]
                         - 9.84375 * s_7[k]
                         + 26.25 * s_9[k]
                         - 9.84375 * s_16[k]
                         + 52.5 * s_18[k]
                         - 31.5 * s_20[k]
                         - 3.28125 * s_29[k]
                         + 26.25 * s_31[k]
                         - 31.5 * s_33[k]
                         + 6.0 * s_35[k];
            }

#pragma omp simd aligned(t_10, s_0, s_3, s_5, s_12, s_14, s_21, s_23, s_27, s_36, s_38, s_40, \
                         s_42 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_10[k] = -f_33 * s_0[k]
                          - f_27 * s_3[k]
                          + f_34 * s_5[k]
                          + f_34 * s_12[k]
                          - f_35 * s_14[k]
                          + f_27 * s_21[k]
                          - f_34 * s_23[k]
                          + f_36 * s_27[k]
                          + f_33 * s_36[k]
                          - f_34 * s_38[k]
                          + f_35 * s_40[k]
                          - f_36 * s_42[k];
            }

#pragma omp simd aligned(t_11, s_2, s_7, s_9, s_16, s_18, s_20, s_29, s_31, \
                         s_33 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_11[k] = f_22 * s_2[k]
                          - f_22 * s_7[k]
                          - f_25 * s_9[k]
                          - f_20 * s_16[k]
                          + f_23 * s_18[k]
                          + f_26 * s_20[k]
                          - f_19 * s_29[k]
                          + f_21 * s_31[k]
                          - f_24 * s_33[k];
            }

#pragma omp simd aligned(t_12, s_0, s_3, s_5, s_10, s_12, s_14, s_21, s_23, s_25, s_36, s_38, \
                         s_40 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_12[k] = f_37 * s_0[k]
                          - f_16 * s_3[k]
                          - f_38 * s_5[k]
                          - f_39 * s_10[k]
                          + f_40 * s_12[k]
                          + f_41 * s_14[k]
                          - f_16 * s_21[k]
                          + f_40 * s_23[k]
                          - f_42 * s_25[k]
                          + f_37 * s_36[k]
                          - f_38 * s_38[k]
                          + f_41 * s_40[k];
            }

#pragma omp simd aligned(t_13, s_2, s_7, s_9, s_16, s_18, s_29, s_31 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_13[k] = -f_14 * s_2[k]
                          + f_12 * s_7[k]
                          + f_15 * s_9[k]
                          + f_10 * s_16[k]
                          - f_13 * s_18[k]
                          - f_10 * s_29[k]
                          + f_11 * s_31[k];
            }

#pragma omp simd aligned(t_14, t_15, t_16, s_0, s_2, s_3, s_5, s_7, s_10, s_12, s_16, s_21, \
                         s_23, s_29, s_36, s_38 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                t_14[k] = -f_43 * s_0[k]
                          + f_7 * s_3[k]
                          + f_7 * s_5[k]
                          - f_44 * s_12[k]
                          - f_7 * s_21[k]
                          + f_44 * s_23[k]
                          + f_43 * s_36[k]
                          - f_7 * s_38[k];

                t_15[k] = f_5 * s_2[k]
                          - f_4 * s_7[k]
                          + f_3 * s_16[k]
                          - f_2 * s_29[k];

                t_16[k] = f_45 * s_0[k]
                          - f_2 * s_3[k]
                          + f_46 * s_10[k]
                          - f_2 * s_21[k]
                          + f_45 * s_36[k];
            }
        }
    }
}

auto
transform_l_outer(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t source,
                  const size_t ncomps, const size_t nmax) -> void
{
    // NOTE: the factors are the shell's own, so they are formed once rather than
    // for every atom pair the shell pair reaches.

    const auto f_0 = 0.1875 * std::sqrt(715.0);
    const auto f_1 = 1.3125 * std::sqrt(715.0);
    const auto f_2 = 0.65625 * std::sqrt(715.0);
    const auto f_3 = 3.28125 * std::sqrt(715.0);
    const auto f_4 = 1.96875 * std::sqrt(715.0);
    const auto f_5 = 0.09375 * std::sqrt(715.0);
    const auto f_6 = 0.09375 * std::sqrt(858.0);
    const auto f_7 = 0.21875 * std::sqrt(858.0);
    const auto f_8 = 1.3125 * std::sqrt(858.0);
    const auto f_9 = 4.375 * std::sqrt(858.0);
    const auto f_10 = 0.46875 * std::sqrt(1001.0);
    const auto f_11 = 1.875 * std::sqrt(1001.0);
    const auto f_12 = 0.84375 * std::sqrt(1001.0);
    const auto f_13 = 3.75 * std::sqrt(1001.0);
    const auto f_14 = 0.09375 * std::sqrt(1001.0);
    const auto f_15 = 0.375 * std::sqrt(1001.0);
    const auto f_16 = 0.1875 * std::sqrt(77.0);
    const auto f_17 = 4.5 * std::sqrt(77.0);
    const auto f_18 = 7.5 * std::sqrt(77.0);
    const auto f_19 = 0.28125 * std::sqrt(1155.0);
    const auto f_20 = 0.46875 * std::sqrt(1155.0);
    const auto f_21 = 1.875 * std::sqrt(1155.0);
    const auto f_22 = 0.09375 * std::sqrt(1155.0);
    const auto f_23 = 1.25 * std::sqrt(1155.0);
    const auto f_24 = 1.5 * std::sqrt(1155.0);
    const auto f_25 = 0.625 * std::sqrt(1155.0);
    const auto f_26 = 0.5 * std::sqrt(1155.0);
    const auto f_27 = 0.09375 * std::sqrt(70.0);
    const auto f_28 = 0.28125 * std::sqrt(70.0);
    const auto f_29 = 2.8125 * std::sqrt(70.0);
    const auto f_30 = 5.625 * std::sqrt(70.0);
    const auto f_31 = 7.5 * std::sqrt(70.0);
    const auto f_32 = 3.0 * std::sqrt(70.0);
    const auto f_33 = 0.046875 * std::sqrt(70.0);
    const auto f_34 = 1.40625 * std::sqrt(70.0);
    const auto f_35 = 3.75 * std::sqrt(70.0);
    const auto f_36 = 1.5 * std::sqrt(70.0);
    const auto f_37 = 0.046875 * std::sqrt(77.0);
    const auto f_38 = 1.125 * std::sqrt(77.0);
    const auto f_39 = 0.46875 * std::sqrt(77.0);
    const auto f_40 = 5.625 * std::sqrt(77.0);
    const auto f_41 = 1.875 * std::sqrt(77.0);
    const auto f_42 = 11.25 * std::sqrt(77.0);
    const auto f_43 = 0.015625 * std::sqrt(858.0);
    const auto f_44 = 3.28125 * std::sqrt(858.0);
    const auto f_45 = 0.0234375 * std::sqrt(715.0);
    const auto f_46 = 1.640625 * std::sqrt(715.0);

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
        auto *g_15 = values + (15 * ncomps + c) * nvalues;
        auto *g_16 = values + (16 * ncomps + c) * nvalues;

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
        const auto *s_36 = buffer.data(source + 36 * ncomps + c);
        const auto *s_37 = buffer.data(source + 37 * ncomps + c);
        const auto *s_38 = buffer.data(source + 38 * ncomps + c);
        const auto *s_39 = buffer.data(source + 39 * ncomps + c);
        const auto *s_40 = buffer.data(source + 40 * ncomps + c);
        const auto *s_41 = buffer.data(source + 41 * ncomps + c);
        const auto *s_42 = buffer.data(source + 42 * ncomps + c);
        const auto *s_43 = buffer.data(source + 43 * ncomps + c);
        const auto *s_44 = buffer.data(source + 44 * ncomps + c);

#pragma omp simd aligned(s_1, s_4, s_6, s_8, s_11, s_15, s_17, s_22, s_28, s_30, \
                         s_37 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = f_0 * s_1[k]
                     - f_1 * s_6[k]
                     + f_1 * s_15[k]
                     - f_0 * s_28[k];

            g_1[k] = f_2 * s_4[k]
                     - f_3 * s_11[k]
                     + f_4 * s_22[k]
                     - f_5 * s_37[k];

            g_2[k] = -f_6 * s_1[k]
                     + f_7 * s_6[k]
                     + f_8 * s_8[k]
                     + f_7 * s_15[k]
                     - f_9 * s_17[k]
                     - f_6 * s_28[k]
                     + f_8 * s_30[k];
        }

#pragma omp simd aligned(s_4, s_11, s_13, s_22, s_24, s_37, s_39 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_3[k] = -f_10 * s_4[k]
                     + f_10 * s_11[k]
                     + f_11 * s_13[k]
                     + f_12 * s_22[k]
                     - f_13 * s_24[k]
                     - f_14 * s_37[k]
                     + f_15 * s_39[k];
        }

#pragma omp simd aligned(s_1, s_6, s_8, s_15, s_19, s_28, s_30, s_32 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_4[k] = f_16 * s_1[k]
                     + f_16 * s_6[k]
                     - f_17 * s_8[k]
                     - f_16 * s_15[k]
                     + f_18 * s_19[k]
                     - f_16 * s_28[k]
                     + f_17 * s_30[k]
                     - f_18 * s_32[k];
        }

#pragma omp simd aligned(s_4, s_11, s_13, s_22, s_24, s_26, s_37, s_39, \
                         s_41 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_5[k] = f_19 * s_4[k]
                     + f_20 * s_11[k]
                     - f_21 * s_13[k]
                     + f_22 * s_22[k]
                     - f_23 * s_24[k]
                     + f_24 * s_26[k]
                     - f_22 * s_37[k]
                     + f_25 * s_39[k]
                     - f_26 * s_41[k];
        }

#pragma omp simd aligned(s_1, s_6, s_8, s_15, s_17, s_19, s_28, s_30, s_32, \
                         s_34 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_6[k] = -f_27 * s_1[k]
                     - f_28 * s_6[k]
                     + f_29 * s_8[k]
                     - f_28 * s_15[k]
                     + f_30 * s_17[k]
                     - f_31 * s_19[k]
                     - f_27 * s_28[k]
                     + f_29 * s_30[k]
                     - f_31 * s_32[k]
                     + f_32 * s_34[k];
        }

#pragma omp simd aligned(s_4, s_11, s_13, s_22, s_24, s_26, s_37, s_39, s_41, \
                         s_43 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_7[k] = -3.28125 * s_4[k]
                     - 9.84375 * s_11[k]
                     + 26.25 * s_13[k]
                     - 9.84375 * s_22[k]
                     + 52.5 * s_24[k]
                     - 31.5 * s_26[k]
                     - 3.28125 * s_37[k]
                     + 26.25 * s_39[k]
                     - 31.5 * s_41[k]
                     + 6.0 * s_43[k];
        }

#pragma omp simd aligned(s_0, s_3, s_5, s_10, s_12, s_14, s_21, s_23, s_25, s_27, s_36, s_38, \
                         s_40, s_42, s_44 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_8[k] = 0.2734375 * s_0[k]
                     + 1.09375 * s_3[k]
                     - 8.75 * s_5[k]
                     + 1.640625 * s_10[k]
                     - 26.25 * s_12[k]
                     + 26.25 * s_14[k]
                     + 1.09375 * s_21[k]
                     - 26.25 * s_23[k]
                     + 52.5 * s_25[k]
                     - 14.0 * s_27[k]
                     + 0.2734375 * s_36[k]
                     - 8.75 * s_38[k]
                     + 26.25 * s_40[k]
                     - 14.0 * s_42[k]
                     + s_44[k];
        }

#pragma omp simd aligned(s_2, s_7, s_9, s_16, s_18, s_20, s_29, s_31, s_33, \
                         s_35 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_9[k] = -3.28125 * s_2[k]
                     - 9.84375 * s_7[k]
                     + 26.25 * s_9[k]
                     - 9.84375 * s_16[k]
                     + 52.5 * s_18[k]
                     - 31.5 * s_20[k]
                     - 3.28125 * s_29[k]
                     + 26.25 * s_31[k]
                     - 31.5 * s_33[k]
                     + 6.0 * s_35[k];
        }

#pragma omp simd aligned(s_0, s_3, s_5, s_12, s_14, s_21, s_23, s_27, s_36, s_38, s_40, \
                         s_42 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_10[k] = -f_33 * s_0[k]
                      - f_27 * s_3[k]
                      + f_34 * s_5[k]
                      + f_34 * s_12[k]
                      - f_35 * s_14[k]
                      + f_27 * s_21[k]
                      - f_34 * s_23[k]
                      + f_36 * s_27[k]
                      + f_33 * s_36[k]
                      - f_34 * s_38[k]
                      + f_35 * s_40[k]
                      - f_36 * s_42[k];
        }

#pragma omp simd aligned(s_2, s_7, s_9, s_16, s_18, s_20, s_29, s_31, \
                         s_33 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_11[k] = f_22 * s_2[k]
                      - f_22 * s_7[k]
                      - f_25 * s_9[k]
                      - f_20 * s_16[k]
                      + f_23 * s_18[k]
                      + f_26 * s_20[k]
                      - f_19 * s_29[k]
                      + f_21 * s_31[k]
                      - f_24 * s_33[k];
        }

#pragma omp simd aligned(s_0, s_3, s_5, s_10, s_12, s_14, s_21, s_23, s_25, s_36, s_38, \
                         s_40 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_12[k] = f_37 * s_0[k]
                      - f_16 * s_3[k]
                      - f_38 * s_5[k]
                      - f_39 * s_10[k]
                      + f_40 * s_12[k]
                      + f_41 * s_14[k]
                      - f_16 * s_21[k]
                      + f_40 * s_23[k]
                      - f_42 * s_25[k]
                      + f_37 * s_36[k]
                      - f_38 * s_38[k]
                      + f_41 * s_40[k];
        }

#pragma omp simd aligned(s_2, s_7, s_9, s_16, s_18, s_29, s_31 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_13[k] = -f_14 * s_2[k]
                      + f_12 * s_7[k]
                      + f_15 * s_9[k]
                      + f_10 * s_16[k]
                      - f_13 * s_18[k]
                      - f_10 * s_29[k]
                      + f_11 * s_31[k];
        }

#pragma omp simd aligned(s_0, s_2, s_3, s_5, s_7, s_10, s_12, s_16, s_21, s_23, s_29, s_36, \
                         s_38 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_14[k] = -f_43 * s_0[k]
                      + f_7 * s_3[k]
                      + f_7 * s_5[k]
                      - f_44 * s_12[k]
                      - f_7 * s_21[k]
                      + f_44 * s_23[k]
                      + f_43 * s_36[k]
                      - f_7 * s_38[k];

            g_15[k] = f_5 * s_2[k]
                      - f_4 * s_7[k]
                      + f_3 * s_16[k]
                      - f_2 * s_29[k];

            g_16[k] = f_45 * s_0[k]
                      - f_2 * s_3[k]
                      + f_46 * s_10[k]
                      - f_2 * s_21[k]
                      + f_45 * s_36[k];
        }
    }
}

auto
transform_l_outer_tri(double *values, const size_t nvalues, CSimdMatrix &buffer,
                      const size_t source, const size_t nmax) -> void
{
    // NOTE: the factors are the shell's own, so they are formed once rather than
    // for every atom pair the shell pair reaches.

    const auto f_0 = 0.1875 * std::sqrt(715.0);
    const auto f_1 = 1.3125 * std::sqrt(715.0);
    const auto f_2 = 0.65625 * std::sqrt(715.0);
    const auto f_3 = 3.28125 * std::sqrt(715.0);
    const auto f_4 = 1.96875 * std::sqrt(715.0);
    const auto f_5 = 0.09375 * std::sqrt(715.0);
    const auto f_6 = 0.09375 * std::sqrt(858.0);
    const auto f_7 = 0.21875 * std::sqrt(858.0);
    const auto f_8 = 1.3125 * std::sqrt(858.0);
    const auto f_9 = 4.375 * std::sqrt(858.0);
    const auto f_10 = 0.46875 * std::sqrt(1001.0);
    const auto f_11 = 1.875 * std::sqrt(1001.0);
    const auto f_12 = 0.84375 * std::sqrt(1001.0);
    const auto f_13 = 3.75 * std::sqrt(1001.0);
    const auto f_14 = 0.09375 * std::sqrt(1001.0);
    const auto f_15 = 0.375 * std::sqrt(1001.0);
    const auto f_16 = 0.1875 * std::sqrt(77.0);
    const auto f_17 = 4.5 * std::sqrt(77.0);
    const auto f_18 = 7.5 * std::sqrt(77.0);
    const auto f_19 = 0.28125 * std::sqrt(1155.0);
    const auto f_20 = 0.46875 * std::sqrt(1155.0);
    const auto f_21 = 1.875 * std::sqrt(1155.0);
    const auto f_22 = 0.09375 * std::sqrt(1155.0);
    const auto f_23 = 1.25 * std::sqrt(1155.0);
    const auto f_24 = 1.5 * std::sqrt(1155.0);
    const auto f_25 = 0.625 * std::sqrt(1155.0);
    const auto f_26 = 0.5 * std::sqrt(1155.0);
    const auto f_27 = 0.09375 * std::sqrt(70.0);
    const auto f_28 = 0.28125 * std::sqrt(70.0);
    const auto f_29 = 2.8125 * std::sqrt(70.0);
    const auto f_30 = 5.625 * std::sqrt(70.0);
    const auto f_31 = 7.5 * std::sqrt(70.0);
    const auto f_32 = 3.0 * std::sqrt(70.0);
    const auto f_33 = 0.046875 * std::sqrt(70.0);
    const auto f_34 = 1.40625 * std::sqrt(70.0);
    const auto f_35 = 3.75 * std::sqrt(70.0);
    const auto f_36 = 1.5 * std::sqrt(70.0);
    const auto f_37 = 0.046875 * std::sqrt(77.0);
    const auto f_38 = 1.125 * std::sqrt(77.0);
    const auto f_39 = 0.46875 * std::sqrt(77.0);
    const auto f_40 = 5.625 * std::sqrt(77.0);
    const auto f_41 = 1.875 * std::sqrt(77.0);
    const auto f_42 = 11.25 * std::sqrt(77.0);
    const auto f_43 = 0.015625 * std::sqrt(858.0);
    const auto f_44 = 3.28125 * std::sqrt(858.0);
    const auto f_45 = 0.0234375 * std::sqrt(715.0);
    const auto f_46 = 1.640625 * std::sqrt(715.0);

    // NOTE: the rows of the values are not aligned, starting at this combination's
    // offset in the values block, so they are kept out of the clause below.

    // NOTE: the block is its own transpose, so a row above the diagonal is the
    // one below it read the other way round and is copied rather than computed.

    auto *d_0 = values + 0 * nvalues;
    auto *d_1 = values + 18 * nvalues;
    auto *d_2 = values + 36 * nvalues;
    auto *d_3 = values + 54 * nvalues;
    auto *d_4 = values + 72 * nvalues;
    auto *d_5 = values + 90 * nvalues;
    auto *d_6 = values + 108 * nvalues;
    auto *d_7 = values + 126 * nvalues;
    auto *d_8 = values + 144 * nvalues;
    auto *d_9 = values + 162 * nvalues;
    auto *d_10 = values + 180 * nvalues;
    auto *d_11 = values + 198 * nvalues;
    auto *d_12 = values + 216 * nvalues;
    auto *d_13 = values + 234 * nvalues;
    auto *d_14 = values + 252 * nvalues;
    auto *d_15 = values + 270 * nvalues;
    auto *d_16 = values + 288 * nvalues;

    const auto *q_0_1 = buffer.data(source + 17);
    const auto *q_0_6 = buffer.data(source + 102);
    const auto *q_0_15 = buffer.data(source + 255);
    const auto *q_0_28 = buffer.data(source + 476);
    const auto *q_1_4 = buffer.data(source + 69);
    const auto *q_1_11 = buffer.data(source + 188);
    const auto *q_1_22 = buffer.data(source + 375);
    const auto *q_1_37 = buffer.data(source + 630);
    const auto *q_2_1 = buffer.data(source + 19);
    const auto *q_2_6 = buffer.data(source + 104);
    const auto *q_2_8 = buffer.data(source + 138);
    const auto *q_2_15 = buffer.data(source + 257);
    const auto *q_2_17 = buffer.data(source + 291);
    const auto *q_2_28 = buffer.data(source + 478);
    const auto *q_2_30 = buffer.data(source + 512);
    const auto *q_3_4 = buffer.data(source + 71);
    const auto *q_3_11 = buffer.data(source + 190);
    const auto *q_3_13 = buffer.data(source + 224);
    const auto *q_3_22 = buffer.data(source + 377);
    const auto *q_3_24 = buffer.data(source + 411);
    const auto *q_3_37 = buffer.data(source + 632);
    const auto *q_3_39 = buffer.data(source + 666);
    const auto *q_4_1 = buffer.data(source + 21);
    const auto *q_4_6 = buffer.data(source + 106);
    const auto *q_4_8 = buffer.data(source + 140);
    const auto *q_4_15 = buffer.data(source + 259);
    const auto *q_4_19 = buffer.data(source + 327);
    const auto *q_4_28 = buffer.data(source + 480);
    const auto *q_4_30 = buffer.data(source + 514);
    const auto *q_4_32 = buffer.data(source + 548);
    const auto *q_5_4 = buffer.data(source + 73);
    const auto *q_5_11 = buffer.data(source + 192);
    const auto *q_5_13 = buffer.data(source + 226);
    const auto *q_5_22 = buffer.data(source + 379);
    const auto *q_5_24 = buffer.data(source + 413);
    const auto *q_5_26 = buffer.data(source + 447);
    const auto *q_5_37 = buffer.data(source + 634);
    const auto *q_5_39 = buffer.data(source + 668);
    const auto *q_5_41 = buffer.data(source + 702);
    const auto *q_6_1 = buffer.data(source + 23);
    const auto *q_6_6 = buffer.data(source + 108);
    const auto *q_6_8 = buffer.data(source + 142);
    const auto *q_6_15 = buffer.data(source + 261);
    const auto *q_6_17 = buffer.data(source + 295);
    const auto *q_6_19 = buffer.data(source + 329);
    const auto *q_6_28 = buffer.data(source + 482);
    const auto *q_6_30 = buffer.data(source + 516);
    const auto *q_6_32 = buffer.data(source + 550);
    const auto *q_6_34 = buffer.data(source + 584);
    const auto *q_7_4 = buffer.data(source + 75);
    const auto *q_7_11 = buffer.data(source + 194);
    const auto *q_7_13 = buffer.data(source + 228);
    const auto *q_7_22 = buffer.data(source + 381);
    const auto *q_7_24 = buffer.data(source + 415);
    const auto *q_7_26 = buffer.data(source + 449);
    const auto *q_7_37 = buffer.data(source + 636);
    const auto *q_7_39 = buffer.data(source + 670);
    const auto *q_7_41 = buffer.data(source + 704);
    const auto *q_7_43 = buffer.data(source + 738);
    const auto *q_8_0 = buffer.data(source + 8);
    const auto *q_8_3 = buffer.data(source + 59);
    const auto *q_8_5 = buffer.data(source + 93);
    const auto *q_8_10 = buffer.data(source + 178);
    const auto *q_8_12 = buffer.data(source + 212);
    const auto *q_8_14 = buffer.data(source + 246);
    const auto *q_8_21 = buffer.data(source + 365);
    const auto *q_8_23 = buffer.data(source + 399);
    const auto *q_8_25 = buffer.data(source + 433);
    const auto *q_8_27 = buffer.data(source + 467);
    const auto *q_8_36 = buffer.data(source + 620);
    const auto *q_8_38 = buffer.data(source + 654);
    const auto *q_8_40 = buffer.data(source + 688);
    const auto *q_8_42 = buffer.data(source + 722);
    const auto *q_8_44 = buffer.data(source + 756);
    const auto *q_9_2 = buffer.data(source + 43);
    const auto *q_9_7 = buffer.data(source + 128);
    const auto *q_9_9 = buffer.data(source + 162);
    const auto *q_9_16 = buffer.data(source + 281);
    const auto *q_9_18 = buffer.data(source + 315);
    const auto *q_9_20 = buffer.data(source + 349);
    const auto *q_9_29 = buffer.data(source + 502);
    const auto *q_9_31 = buffer.data(source + 536);
    const auto *q_9_33 = buffer.data(source + 570);
    const auto *q_9_35 = buffer.data(source + 604);
    const auto *q_10_0 = buffer.data(source + 10);
    const auto *q_10_3 = buffer.data(source + 61);
    const auto *q_10_5 = buffer.data(source + 95);
    const auto *q_10_12 = buffer.data(source + 214);
    const auto *q_10_14 = buffer.data(source + 248);
    const auto *q_10_21 = buffer.data(source + 367);
    const auto *q_10_23 = buffer.data(source + 401);
    const auto *q_10_27 = buffer.data(source + 469);
    const auto *q_10_36 = buffer.data(source + 622);
    const auto *q_10_38 = buffer.data(source + 656);
    const auto *q_10_40 = buffer.data(source + 690);
    const auto *q_10_42 = buffer.data(source + 724);
    const auto *q_11_2 = buffer.data(source + 45);
    const auto *q_11_7 = buffer.data(source + 130);
    const auto *q_11_9 = buffer.data(source + 164);
    const auto *q_11_16 = buffer.data(source + 283);
    const auto *q_11_18 = buffer.data(source + 317);
    const auto *q_11_20 = buffer.data(source + 351);
    const auto *q_11_29 = buffer.data(source + 504);
    const auto *q_11_31 = buffer.data(source + 538);
    const auto *q_11_33 = buffer.data(source + 572);
    const auto *q_12_0 = buffer.data(source + 12);
    const auto *q_12_3 = buffer.data(source + 63);
    const auto *q_12_5 = buffer.data(source + 97);
    const auto *q_12_10 = buffer.data(source + 182);
    const auto *q_12_12 = buffer.data(source + 216);
    const auto *q_12_14 = buffer.data(source + 250);
    const auto *q_12_21 = buffer.data(source + 369);
    const auto *q_12_23 = buffer.data(source + 403);
    const auto *q_12_25 = buffer.data(source + 437);
    const auto *q_12_36 = buffer.data(source + 624);
    const auto *q_12_38 = buffer.data(source + 658);
    const auto *q_12_40 = buffer.data(source + 692);
    const auto *q_13_2 = buffer.data(source + 47);
    const auto *q_13_7 = buffer.data(source + 132);
    const auto *q_13_9 = buffer.data(source + 166);
    const auto *q_13_16 = buffer.data(source + 285);
    const auto *q_13_18 = buffer.data(source + 319);
    const auto *q_13_29 = buffer.data(source + 506);
    const auto *q_13_31 = buffer.data(source + 540);
    const auto *q_14_0 = buffer.data(source + 14);
    const auto *q_14_3 = buffer.data(source + 65);
    const auto *q_14_5 = buffer.data(source + 99);
    const auto *q_14_12 = buffer.data(source + 218);
    const auto *q_14_21 = buffer.data(source + 371);
    const auto *q_14_23 = buffer.data(source + 405);
    const auto *q_14_36 = buffer.data(source + 626);
    const auto *q_14_38 = buffer.data(source + 660);
    const auto *q_15_2 = buffer.data(source + 49);
    const auto *q_15_7 = buffer.data(source + 134);
    const auto *q_15_16 = buffer.data(source + 287);
    const auto *q_15_29 = buffer.data(source + 508);
    const auto *q_16_0 = buffer.data(source + 16);
    const auto *q_16_3 = buffer.data(source + 67);
    const auto *q_16_10 = buffer.data(source + 186);
    const auto *q_16_21 = buffer.data(source + 373);
    const auto *q_16_36 = buffer.data(source + 628);

#pragma omp simd aligned(q_0_1, q_0_6, q_0_15, q_0_28, q_1_4, q_1_11, q_1_22, \
                         q_1_37 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        d_0[k] = f_0 * q_0_1[k]
                 - f_1 * q_0_6[k]
                 + f_1 * q_0_15[k]
                 - f_0 * q_0_28[k];

        d_1[k] = f_2 * q_1_4[k]
                 - f_3 * q_1_11[k]
                 + f_4 * q_1_22[k]
                 - f_5 * q_1_37[k];
    }

#pragma omp simd aligned(q_2_1, q_2_6, q_2_8, q_2_15, q_2_17, q_2_28, q_2_30, q_3_4, q_3_11, \
                         q_3_13, q_3_22, q_3_24, q_3_37, q_3_39 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        d_2[k] = -f_6 * q_2_1[k]
                 + f_7 * q_2_6[k]
                 + f_8 * q_2_8[k]
                 + f_7 * q_2_15[k]
                 - f_9 * q_2_17[k]
                 - f_6 * q_2_28[k]
                 + f_8 * q_2_30[k];

        d_3[k] = -f_10 * q_3_4[k]
                 + f_10 * q_3_11[k]
                 + f_11 * q_3_13[k]
                 + f_12 * q_3_22[k]
                 - f_13 * q_3_24[k]
                 - f_14 * q_3_37[k]
                 + f_15 * q_3_39[k];
    }

#pragma omp simd aligned(q_4_1, q_4_6, q_4_8, q_4_15, q_4_19, q_4_28, q_4_30, \
                         q_4_32 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        d_4[k] = f_16 * q_4_1[k]
                 + f_16 * q_4_6[k]
                 - f_17 * q_4_8[k]
                 - f_16 * q_4_15[k]
                 + f_18 * q_4_19[k]
                 - f_16 * q_4_28[k]
                 + f_17 * q_4_30[k]
                 - f_18 * q_4_32[k];
    }

#pragma omp simd aligned(q_5_4, q_5_11, q_5_13, q_5_22, q_5_24, q_5_26, q_5_37, q_5_39, \
                         q_5_41 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        d_5[k] = f_19 * q_5_4[k]
                 + f_20 * q_5_11[k]
                 - f_21 * q_5_13[k]
                 + f_22 * q_5_22[k]
                 - f_23 * q_5_24[k]
                 + f_24 * q_5_26[k]
                 - f_22 * q_5_37[k]
                 + f_25 * q_5_39[k]
                 - f_26 * q_5_41[k];
    }

#pragma omp simd aligned(q_6_1, q_6_6, q_6_8, q_6_15, q_6_17, q_6_19, q_6_28, q_6_30, q_6_32, \
                         q_6_34 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        d_6[k] = -f_27 * q_6_1[k]
                 - f_28 * q_6_6[k]
                 + f_29 * q_6_8[k]
                 - f_28 * q_6_15[k]
                 + f_30 * q_6_17[k]
                 - f_31 * q_6_19[k]
                 - f_27 * q_6_28[k]
                 + f_29 * q_6_30[k]
                 - f_31 * q_6_32[k]
                 + f_32 * q_6_34[k];
    }

#pragma omp simd aligned(q_7_4, q_7_11, q_7_13, q_7_22, q_7_24, q_7_26, q_7_37, q_7_39, \
                         q_7_41, q_7_43 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        d_7[k] = -3.28125 * q_7_4[k]
                 - 9.84375 * q_7_11[k]
                 + 26.25 * q_7_13[k]
                 - 9.84375 * q_7_22[k]
                 + 52.5 * q_7_24[k]
                 - 31.5 * q_7_26[k]
                 - 3.28125 * q_7_37[k]
                 + 26.25 * q_7_39[k]
                 - 31.5 * q_7_41[k]
                 + 6.0 * q_7_43[k];
    }

#pragma omp simd aligned(q_8_0, q_8_3, q_8_5, q_8_10, q_8_12, q_8_14, q_8_21, q_8_23, q_8_25, \
                         q_8_27, q_8_36, q_8_38, q_8_40, q_8_42, \
                         q_8_44 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        d_8[k] = 0.2734375 * q_8_0[k]
                 + 1.09375 * q_8_3[k]
                 - 8.75 * q_8_5[k]
                 + 1.640625 * q_8_10[k]
                 - 26.25 * q_8_12[k]
                 + 26.25 * q_8_14[k]
                 + 1.09375 * q_8_21[k]
                 - 26.25 * q_8_23[k]
                 + 52.5 * q_8_25[k]
                 - 14.0 * q_8_27[k]
                 + 0.2734375 * q_8_36[k]
                 - 8.75 * q_8_38[k]
                 + 26.25 * q_8_40[k]
                 - 14.0 * q_8_42[k]
                 + q_8_44[k];
    }

#pragma omp simd aligned(q_9_2, q_9_7, q_9_9, q_9_16, q_9_18, q_9_20, q_9_29, q_9_31, q_9_33, \
                         q_9_35 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        d_9[k] = -3.28125 * q_9_2[k]
                 - 9.84375 * q_9_7[k]
                 + 26.25 * q_9_9[k]
                 - 9.84375 * q_9_16[k]
                 + 52.5 * q_9_18[k]
                 - 31.5 * q_9_20[k]
                 - 3.28125 * q_9_29[k]
                 + 26.25 * q_9_31[k]
                 - 31.5 * q_9_33[k]
                 + 6.0 * q_9_35[k];
    }

#pragma omp simd aligned(q_10_0, q_10_3, q_10_5, q_10_12, q_10_14, q_10_21, q_10_23, q_10_27, \
                         q_10_36, q_10_38, q_10_40, q_10_42 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        d_10[k] = -f_33 * q_10_0[k]
                  - f_27 * q_10_3[k]
                  + f_34 * q_10_5[k]
                  + f_34 * q_10_12[k]
                  - f_35 * q_10_14[k]
                  + f_27 * q_10_21[k]
                  - f_34 * q_10_23[k]
                  + f_36 * q_10_27[k]
                  + f_33 * q_10_36[k]
                  - f_34 * q_10_38[k]
                  + f_35 * q_10_40[k]
                  - f_36 * q_10_42[k];
    }

#pragma omp simd aligned(q_11_2, q_11_7, q_11_9, q_11_16, q_11_18, q_11_20, q_11_29, q_11_31, \
                         q_11_33 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        d_11[k] = f_22 * q_11_2[k]
                  - f_22 * q_11_7[k]
                  - f_25 * q_11_9[k]
                  - f_20 * q_11_16[k]
                  + f_23 * q_11_18[k]
                  + f_26 * q_11_20[k]
                  - f_19 * q_11_29[k]
                  + f_21 * q_11_31[k]
                  - f_24 * q_11_33[k];
    }

#pragma omp simd aligned(q_12_0, q_12_3, q_12_5, q_12_10, q_12_12, q_12_14, q_12_21, q_12_23, \
                         q_12_25, q_12_36, q_12_38, q_12_40 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        d_12[k] = f_37 * q_12_0[k]
                  - f_16 * q_12_3[k]
                  - f_38 * q_12_5[k]
                  - f_39 * q_12_10[k]
                  + f_40 * q_12_12[k]
                  + f_41 * q_12_14[k]
                  - f_16 * q_12_21[k]
                  + f_40 * q_12_23[k]
                  - f_42 * q_12_25[k]
                  + f_37 * q_12_36[k]
                  - f_38 * q_12_38[k]
                  + f_41 * q_12_40[k];
    }

#pragma omp simd aligned(q_13_2, q_13_7, q_13_9, q_13_16, q_13_18, q_13_29, \
                         q_13_31 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        d_13[k] = -f_14 * q_13_2[k]
                  + f_12 * q_13_7[k]
                  + f_15 * q_13_9[k]
                  + f_10 * q_13_16[k]
                  - f_13 * q_13_18[k]
                  - f_10 * q_13_29[k]
                  + f_11 * q_13_31[k];
    }

#pragma omp simd aligned(q_14_0, q_14_3, q_14_5, q_14_12, q_14_21, q_14_23, q_14_36, q_14_38, \
                         q_15_2, q_15_7, q_15_16, q_15_29 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        d_14[k] = -f_43 * q_14_0[k]
                  + f_7 * q_14_3[k]
                  + f_7 * q_14_5[k]
                  - f_44 * q_14_12[k]
                  - f_7 * q_14_21[k]
                  + f_44 * q_14_23[k]
                  + f_43 * q_14_36[k]
                  - f_7 * q_14_38[k];

        d_15[k] = f_5 * q_15_2[k]
                  - f_4 * q_15_7[k]
                  + f_3 * q_15_16[k]
                  - f_2 * q_15_29[k];
    }

#pragma omp simd aligned(q_16_0, q_16_3, q_16_10, q_16_21, q_16_36 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        d_16[k] = f_45 * q_16_0[k]
                  - f_2 * q_16_3[k]
                  + f_46 * q_16_10[k]
                  - f_2 * q_16_21[k]
                  + f_45 * q_16_36[k];
    }

    for (size_t c = 1; c < 17; c++)
    {
        auto *g_0 = values + (0 + c) * nvalues;
        auto *g_1 = values + (c * 17 + 0) * nvalues;

        const auto *s_1 = buffer.data(source + 17 + c);
        const auto *s_6 = buffer.data(source + 102 + c);
        const auto *s_15 = buffer.data(source + 255 + c);
        const auto *s_28 = buffer.data(source + 476 + c);

#pragma omp simd aligned(s_1, s_6, s_15, s_28 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = f_0 * s_1[k]
                     - f_1 * s_6[k]
                     + f_1 * s_15[k]
                     - f_0 * s_28[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 2; c < 17; c++)
    {
        auto *g_0 = values + (17 + c) * nvalues;
        auto *g_1 = values + (c * 17 + 1) * nvalues;

        const auto *s_4 = buffer.data(source + 68 + c);
        const auto *s_11 = buffer.data(source + 187 + c);
        const auto *s_22 = buffer.data(source + 374 + c);
        const auto *s_37 = buffer.data(source + 629 + c);

#pragma omp simd aligned(s_4, s_11, s_22, s_37 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = f_2 * s_4[k]
                     - f_3 * s_11[k]
                     + f_4 * s_22[k]
                     - f_5 * s_37[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 3; c < 17; c++)
    {
        auto *g_0 = values + (34 + c) * nvalues;
        auto *g_1 = values + (c * 17 + 2) * nvalues;

        const auto *s_1 = buffer.data(source + 17 + c);
        const auto *s_6 = buffer.data(source + 102 + c);
        const auto *s_8 = buffer.data(source + 136 + c);
        const auto *s_15 = buffer.data(source + 255 + c);
        const auto *s_17 = buffer.data(source + 289 + c);
        const auto *s_28 = buffer.data(source + 476 + c);
        const auto *s_30 = buffer.data(source + 510 + c);

#pragma omp simd aligned(s_1, s_6, s_8, s_15, s_17, s_28, s_30 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = -f_6 * s_1[k]
                     + f_7 * s_6[k]
                     + f_8 * s_8[k]
                     + f_7 * s_15[k]
                     - f_9 * s_17[k]
                     - f_6 * s_28[k]
                     + f_8 * s_30[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 4; c < 17; c++)
    {
        auto *g_0 = values + (51 + c) * nvalues;
        auto *g_1 = values + (c * 17 + 3) * nvalues;

        const auto *s_4 = buffer.data(source + 68 + c);
        const auto *s_11 = buffer.data(source + 187 + c);
        const auto *s_13 = buffer.data(source + 221 + c);
        const auto *s_22 = buffer.data(source + 374 + c);
        const auto *s_24 = buffer.data(source + 408 + c);
        const auto *s_37 = buffer.data(source + 629 + c);
        const auto *s_39 = buffer.data(source + 663 + c);

#pragma omp simd aligned(s_4, s_11, s_13, s_22, s_24, s_37, s_39 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = -f_10 * s_4[k]
                     + f_10 * s_11[k]
                     + f_11 * s_13[k]
                     + f_12 * s_22[k]
                     - f_13 * s_24[k]
                     - f_14 * s_37[k]
                     + f_15 * s_39[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 5; c < 17; c++)
    {
        auto *g_0 = values + (68 + c) * nvalues;
        auto *g_1 = values + (c * 17 + 4) * nvalues;

        const auto *s_1 = buffer.data(source + 17 + c);
        const auto *s_6 = buffer.data(source + 102 + c);
        const auto *s_8 = buffer.data(source + 136 + c);
        const auto *s_15 = buffer.data(source + 255 + c);
        const auto *s_19 = buffer.data(source + 323 + c);
        const auto *s_28 = buffer.data(source + 476 + c);
        const auto *s_30 = buffer.data(source + 510 + c);
        const auto *s_32 = buffer.data(source + 544 + c);

#pragma omp simd aligned(s_1, s_6, s_8, s_15, s_19, s_28, s_30, s_32 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = f_16 * s_1[k]
                     + f_16 * s_6[k]
                     - f_17 * s_8[k]
                     - f_16 * s_15[k]
                     + f_18 * s_19[k]
                     - f_16 * s_28[k]
                     + f_17 * s_30[k]
                     - f_18 * s_32[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 6; c < 17; c++)
    {
        auto *g_0 = values + (85 + c) * nvalues;
        auto *g_1 = values + (c * 17 + 5) * nvalues;

        const auto *s_4 = buffer.data(source + 68 + c);
        const auto *s_11 = buffer.data(source + 187 + c);
        const auto *s_13 = buffer.data(source + 221 + c);
        const auto *s_22 = buffer.data(source + 374 + c);
        const auto *s_24 = buffer.data(source + 408 + c);
        const auto *s_26 = buffer.data(source + 442 + c);
        const auto *s_37 = buffer.data(source + 629 + c);
        const auto *s_39 = buffer.data(source + 663 + c);
        const auto *s_41 = buffer.data(source + 697 + c);

#pragma omp simd aligned(s_4, s_11, s_13, s_22, s_24, s_26, s_37, s_39, \
                         s_41 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = f_19 * s_4[k]
                     + f_20 * s_11[k]
                     - f_21 * s_13[k]
                     + f_22 * s_22[k]
                     - f_23 * s_24[k]
                     + f_24 * s_26[k]
                     - f_22 * s_37[k]
                     + f_25 * s_39[k]
                     - f_26 * s_41[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 7; c < 17; c++)
    {
        auto *g_0 = values + (102 + c) * nvalues;
        auto *g_1 = values + (c * 17 + 6) * nvalues;

        const auto *s_1 = buffer.data(source + 17 + c);
        const auto *s_6 = buffer.data(source + 102 + c);
        const auto *s_8 = buffer.data(source + 136 + c);
        const auto *s_15 = buffer.data(source + 255 + c);
        const auto *s_17 = buffer.data(source + 289 + c);
        const auto *s_19 = buffer.data(source + 323 + c);
        const auto *s_28 = buffer.data(source + 476 + c);
        const auto *s_30 = buffer.data(source + 510 + c);
        const auto *s_32 = buffer.data(source + 544 + c);
        const auto *s_34 = buffer.data(source + 578 + c);

#pragma omp simd aligned(s_1, s_6, s_8, s_15, s_17, s_19, s_28, s_30, s_32, \
                         s_34 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = -f_27 * s_1[k]
                     - f_28 * s_6[k]
                     + f_29 * s_8[k]
                     - f_28 * s_15[k]
                     + f_30 * s_17[k]
                     - f_31 * s_19[k]
                     - f_27 * s_28[k]
                     + f_29 * s_30[k]
                     - f_31 * s_32[k]
                     + f_32 * s_34[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 8; c < 17; c++)
    {
        auto *g_0 = values + (119 + c) * nvalues;
        auto *g_1 = values + (c * 17 + 7) * nvalues;

        const auto *s_4 = buffer.data(source + 68 + c);
        const auto *s_11 = buffer.data(source + 187 + c);
        const auto *s_13 = buffer.data(source + 221 + c);
        const auto *s_22 = buffer.data(source + 374 + c);
        const auto *s_24 = buffer.data(source + 408 + c);
        const auto *s_26 = buffer.data(source + 442 + c);
        const auto *s_37 = buffer.data(source + 629 + c);
        const auto *s_39 = buffer.data(source + 663 + c);
        const auto *s_41 = buffer.data(source + 697 + c);
        const auto *s_43 = buffer.data(source + 731 + c);

#pragma omp simd aligned(s_4, s_11, s_13, s_22, s_24, s_26, s_37, s_39, s_41, \
                         s_43 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = -3.28125 * s_4[k]
                     - 9.84375 * s_11[k]
                     + 26.25 * s_13[k]
                     - 9.84375 * s_22[k]
                     + 52.5 * s_24[k]
                     - 31.5 * s_26[k]
                     - 3.28125 * s_37[k]
                     + 26.25 * s_39[k]
                     - 31.5 * s_41[k]
                     + 6.0 * s_43[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 9; c < 17; c++)
    {
        auto *g_0 = values + (136 + c) * nvalues;
        auto *g_1 = values + (c * 17 + 8) * nvalues;

        const auto *s_0 = buffer.data(source + 0 + c);
        const auto *s_3 = buffer.data(source + 51 + c);
        const auto *s_5 = buffer.data(source + 85 + c);
        const auto *s_10 = buffer.data(source + 170 + c);
        const auto *s_12 = buffer.data(source + 204 + c);
        const auto *s_14 = buffer.data(source + 238 + c);
        const auto *s_21 = buffer.data(source + 357 + c);
        const auto *s_23 = buffer.data(source + 391 + c);
        const auto *s_25 = buffer.data(source + 425 + c);
        const auto *s_27 = buffer.data(source + 459 + c);
        const auto *s_36 = buffer.data(source + 612 + c);
        const auto *s_38 = buffer.data(source + 646 + c);
        const auto *s_40 = buffer.data(source + 680 + c);
        const auto *s_42 = buffer.data(source + 714 + c);
        const auto *s_44 = buffer.data(source + 748 + c);

#pragma omp simd aligned(s_0, s_3, s_5, s_10, s_12, s_14, s_21, s_23, s_25, s_27, s_36, s_38, \
                         s_40, s_42, s_44 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = 0.2734375 * s_0[k]
                     + 1.09375 * s_3[k]
                     - 8.75 * s_5[k]
                     + 1.640625 * s_10[k]
                     - 26.25 * s_12[k]
                     + 26.25 * s_14[k]
                     + 1.09375 * s_21[k]
                     - 26.25 * s_23[k]
                     + 52.5 * s_25[k]
                     - 14.0 * s_27[k]
                     + 0.2734375 * s_36[k]
                     - 8.75 * s_38[k]
                     + 26.25 * s_40[k]
                     - 14.0 * s_42[k]
                     + s_44[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 10; c < 17; c++)
    {
        auto *g_0 = values + (153 + c) * nvalues;
        auto *g_1 = values + (c * 17 + 9) * nvalues;

        const auto *s_2 = buffer.data(source + 34 + c);
        const auto *s_7 = buffer.data(source + 119 + c);
        const auto *s_9 = buffer.data(source + 153 + c);
        const auto *s_16 = buffer.data(source + 272 + c);
        const auto *s_18 = buffer.data(source + 306 + c);
        const auto *s_20 = buffer.data(source + 340 + c);
        const auto *s_29 = buffer.data(source + 493 + c);
        const auto *s_31 = buffer.data(source + 527 + c);
        const auto *s_33 = buffer.data(source + 561 + c);
        const auto *s_35 = buffer.data(source + 595 + c);

#pragma omp simd aligned(s_2, s_7, s_9, s_16, s_18, s_20, s_29, s_31, s_33, \
                         s_35 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = -3.28125 * s_2[k]
                     - 9.84375 * s_7[k]
                     + 26.25 * s_9[k]
                     - 9.84375 * s_16[k]
                     + 52.5 * s_18[k]
                     - 31.5 * s_20[k]
                     - 3.28125 * s_29[k]
                     + 26.25 * s_31[k]
                     - 31.5 * s_33[k]
                     + 6.0 * s_35[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 11; c < 17; c++)
    {
        auto *g_0 = values + (170 + c) * nvalues;
        auto *g_1 = values + (c * 17 + 10) * nvalues;

        const auto *s_0 = buffer.data(source + 0 + c);
        const auto *s_3 = buffer.data(source + 51 + c);
        const auto *s_5 = buffer.data(source + 85 + c);
        const auto *s_12 = buffer.data(source + 204 + c);
        const auto *s_14 = buffer.data(source + 238 + c);
        const auto *s_21 = buffer.data(source + 357 + c);
        const auto *s_23 = buffer.data(source + 391 + c);
        const auto *s_27 = buffer.data(source + 459 + c);
        const auto *s_36 = buffer.data(source + 612 + c);
        const auto *s_38 = buffer.data(source + 646 + c);
        const auto *s_40 = buffer.data(source + 680 + c);
        const auto *s_42 = buffer.data(source + 714 + c);

#pragma omp simd aligned(s_0, s_3, s_5, s_12, s_14, s_21, s_23, s_27, s_36, s_38, s_40, \
                         s_42 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = -f_33 * s_0[k]
                     - f_27 * s_3[k]
                     + f_34 * s_5[k]
                     + f_34 * s_12[k]
                     - f_35 * s_14[k]
                     + f_27 * s_21[k]
                     - f_34 * s_23[k]
                     + f_36 * s_27[k]
                     + f_33 * s_36[k]
                     - f_34 * s_38[k]
                     + f_35 * s_40[k]
                     - f_36 * s_42[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 12; c < 17; c++)
    {
        auto *g_0 = values + (187 + c) * nvalues;
        auto *g_1 = values + (c * 17 + 11) * nvalues;

        const auto *s_2 = buffer.data(source + 34 + c);
        const auto *s_7 = buffer.data(source + 119 + c);
        const auto *s_9 = buffer.data(source + 153 + c);
        const auto *s_16 = buffer.data(source + 272 + c);
        const auto *s_18 = buffer.data(source + 306 + c);
        const auto *s_20 = buffer.data(source + 340 + c);
        const auto *s_29 = buffer.data(source + 493 + c);
        const auto *s_31 = buffer.data(source + 527 + c);
        const auto *s_33 = buffer.data(source + 561 + c);

#pragma omp simd aligned(s_2, s_7, s_9, s_16, s_18, s_20, s_29, s_31, \
                         s_33 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = f_22 * s_2[k]
                     - f_22 * s_7[k]
                     - f_25 * s_9[k]
                     - f_20 * s_16[k]
                     + f_23 * s_18[k]
                     + f_26 * s_20[k]
                     - f_19 * s_29[k]
                     + f_21 * s_31[k]
                     - f_24 * s_33[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 13; c < 17; c++)
    {
        auto *g_0 = values + (204 + c) * nvalues;
        auto *g_1 = values + (c * 17 + 12) * nvalues;

        const auto *s_0 = buffer.data(source + 0 + c);
        const auto *s_3 = buffer.data(source + 51 + c);
        const auto *s_5 = buffer.data(source + 85 + c);
        const auto *s_10 = buffer.data(source + 170 + c);
        const auto *s_12 = buffer.data(source + 204 + c);
        const auto *s_14 = buffer.data(source + 238 + c);
        const auto *s_21 = buffer.data(source + 357 + c);
        const auto *s_23 = buffer.data(source + 391 + c);
        const auto *s_25 = buffer.data(source + 425 + c);
        const auto *s_36 = buffer.data(source + 612 + c);
        const auto *s_38 = buffer.data(source + 646 + c);
        const auto *s_40 = buffer.data(source + 680 + c);

#pragma omp simd aligned(s_0, s_3, s_5, s_10, s_12, s_14, s_21, s_23, s_25, s_36, s_38, \
                         s_40 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = f_37 * s_0[k]
                     - f_16 * s_3[k]
                     - f_38 * s_5[k]
                     - f_39 * s_10[k]
                     + f_40 * s_12[k]
                     + f_41 * s_14[k]
                     - f_16 * s_21[k]
                     + f_40 * s_23[k]
                     - f_42 * s_25[k]
                     + f_37 * s_36[k]
                     - f_38 * s_38[k]
                     + f_41 * s_40[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 14; c < 17; c++)
    {
        auto *g_0 = values + (221 + c) * nvalues;
        auto *g_1 = values + (c * 17 + 13) * nvalues;

        const auto *s_2 = buffer.data(source + 34 + c);
        const auto *s_7 = buffer.data(source + 119 + c);
        const auto *s_9 = buffer.data(source + 153 + c);
        const auto *s_16 = buffer.data(source + 272 + c);
        const auto *s_18 = buffer.data(source + 306 + c);
        const auto *s_29 = buffer.data(source + 493 + c);
        const auto *s_31 = buffer.data(source + 527 + c);

#pragma omp simd aligned(s_2, s_7, s_9, s_16, s_18, s_29, s_31 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = -f_14 * s_2[k]
                     + f_12 * s_7[k]
                     + f_15 * s_9[k]
                     + f_10 * s_16[k]
                     - f_13 * s_18[k]
                     - f_10 * s_29[k]
                     + f_11 * s_31[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 15; c < 17; c++)
    {
        auto *g_0 = values + (238 + c) * nvalues;
        auto *g_1 = values + (c * 17 + 14) * nvalues;

        const auto *s_0 = buffer.data(source + 0 + c);
        const auto *s_3 = buffer.data(source + 51 + c);
        const auto *s_5 = buffer.data(source + 85 + c);
        const auto *s_12 = buffer.data(source + 204 + c);
        const auto *s_21 = buffer.data(source + 357 + c);
        const auto *s_23 = buffer.data(source + 391 + c);
        const auto *s_36 = buffer.data(source + 612 + c);
        const auto *s_38 = buffer.data(source + 646 + c);

#pragma omp simd aligned(s_0, s_3, s_5, s_12, s_21, s_23, s_36, s_38 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = -f_43 * s_0[k]
                     + f_7 * s_3[k]
                     + f_7 * s_5[k]
                     - f_44 * s_12[k]
                     - f_7 * s_21[k]
                     + f_44 * s_23[k]
                     + f_43 * s_36[k]
                     - f_7 * s_38[k];
            g_1[k] = g_0[k];
        }
    }

    for (size_t c = 16; c < 17; c++)
    {
        auto *g_0 = values + (255 + c) * nvalues;
        auto *g_1 = values + (c * 17 + 15) * nvalues;

        const auto *s_2 = buffer.data(source + 34 + c);
        const auto *s_7 = buffer.data(source + 119 + c);
        const auto *s_16 = buffer.data(source + 272 + c);
        const auto *s_29 = buffer.data(source + 493 + c);

#pragma omp simd aligned(s_2, s_7, s_16, s_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            g_0[k] = f_5 * s_2[k]
                     - f_4 * s_7[k]
                     + f_3 * s_16[k]
                     - f_2 * s_29[k];
            g_1[k] = g_0[k];
        }
    }
}

}  // namespace simdtrf
