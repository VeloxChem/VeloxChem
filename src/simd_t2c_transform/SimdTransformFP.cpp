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


#include "SimdTransformFP.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_fp(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t fp,
             const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 0.75 * std::sqrt(10.0);
    const auto f_1 = 0.25 * std::sqrt(10.0);
    const auto f_2 = std::sqrt(15.0);
    const auto f_3 = 0.25 * std::sqrt(6.0);
    const auto f_4 = std::sqrt(6.0);
    const auto f_5 = 0.5 * std::sqrt(15.0);

    // NOTE: the rows of the values are not aligned, starting at this combination's
    // offset in the values block, so they are kept out of the clause below.

    auto *g_0 = values + 0 * nvalues;
    auto *g_1 = values + 1 * nvalues;
    auto *g_2 = values + 2 * nvalues;
    auto *g_3 = values + 3 * nvalues;
    auto *g_4 = values + 4 * nvalues;
    auto *g_5 = values + 5 * nvalues;
    auto *g_6 = values + 6 * nvalues;
    auto *g_7 = values + 7 * nvalues;
    auto *g_8 = values + 8 * nvalues;
    auto *g_9 = values + 9 * nvalues;
    auto *g_10 = values + 10 * nvalues;
    auto *g_11 = values + 11 * nvalues;
    auto *g_12 = values + 12 * nvalues;
    auto *g_13 = values + 13 * nvalues;
    auto *g_14 = values + 14 * nvalues;
    auto *g_15 = values + 15 * nvalues;
    auto *g_16 = values + 16 * nvalues;
    auto *g_17 = values + 17 * nvalues;
    auto *g_18 = values + 18 * nvalues;
    auto *g_19 = values + 19 * nvalues;
    auto *g_20 = values + 20 * nvalues;

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_1 = buffer.data(fp + 1);
    const auto *fp_2 = buffer.data(fp + 2);
    const auto *fp_3 = buffer.data(fp + 3);
    const auto *fp_4 = buffer.data(fp + 4);
    const auto *fp_5 = buffer.data(fp + 5);
    const auto *fp_6 = buffer.data(fp + 6);
    const auto *fp_7 = buffer.data(fp + 7);
    const auto *fp_8 = buffer.data(fp + 8);
    const auto *fp_9 = buffer.data(fp + 9);
    const auto *fp_10 = buffer.data(fp + 10);
    const auto *fp_11 = buffer.data(fp + 11);
    const auto *fp_12 = buffer.data(fp + 12);
    const auto *fp_13 = buffer.data(fp + 13);
    const auto *fp_14 = buffer.data(fp + 14);
    const auto *fp_15 = buffer.data(fp + 15);
    const auto *fp_16 = buffer.data(fp + 16);
    const auto *fp_17 = buffer.data(fp + 17);
    const auto *fp_18 = buffer.data(fp + 18);
    const auto *fp_19 = buffer.data(fp + 19);
    const auto *fp_20 = buffer.data(fp + 20);
    const auto *fp_21 = buffer.data(fp + 21);
    const auto *fp_22 = buffer.data(fp + 22);
    const auto *fp_23 = buffer.data(fp + 23);
    const auto *fp_24 = buffer.data(fp + 24);
    const auto *fp_25 = buffer.data(fp + 25);
    const auto *fp_26 = buffer.data(fp + 26);
    const auto *fp_27 = buffer.data(fp + 27);
    const auto *fp_28 = buffer.data(fp + 28);
    const auto *fp_29 = buffer.data(fp + 29);

#pragma omp simd aligned(fp_3, fp_4, fp_5, fp_12, fp_13, fp_14, fp_18, fp_19, \
                         fp_20 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * fp_4[k]
                 - f_1 * fp_19[k];

        g_1[k] = f_0 * fp_5[k]
                 - f_1 * fp_20[k];

        g_2[k] = f_0 * fp_3[k]
                 - f_1 * fp_18[k];

        g_3[k] = f_2 * fp_13[k];

        g_4[k] = f_2 * fp_14[k];

        g_5[k] = f_2 * fp_12[k];
    }

#pragma omp simd aligned(fp_3, fp_4, fp_5, fp_7, fp_18, fp_19, fp_20, fp_22, fp_24, fp_25, \
                         fp_26, fp_28 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = -f_3 * fp_4[k]
                 - f_3 * fp_19[k]
                 + f_4 * fp_25[k];

        g_7[k] = -f_3 * fp_5[k]
                 - f_3 * fp_20[k]
                 + f_4 * fp_26[k];

        g_8[k] = -f_3 * fp_3[k]
                 - f_3 * fp_18[k]
                 + f_4 * fp_24[k];

        g_9[k] = -1.5 * fp_7[k]
                 - 1.5 * fp_22[k]
                 + fp_28[k];
    }

#pragma omp simd aligned(fp_1, fp_2, fp_6, fp_8, fp_10, fp_11, fp_16, fp_17, fp_21, fp_23, \
                         fp_27, fp_29 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = -1.5 * fp_8[k]
                  - 1.5 * fp_23[k]
                  + fp_29[k];

        g_11[k] = -1.5 * fp_6[k]
                  - 1.5 * fp_21[k]
                  + fp_27[k];

        g_12[k] = -f_3 * fp_1[k]
                  - f_3 * fp_10[k]
                  + f_4 * fp_16[k];

        g_13[k] = -f_3 * fp_2[k]
                  - f_3 * fp_11[k]
                  + f_4 * fp_17[k];
    }

#pragma omp simd aligned(fp_0, fp_1, fp_6, fp_7, fp_8, fp_9, fp_10, fp_15, fp_21, fp_22, \
                         fp_23 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = -f_3 * fp_0[k]
                  - f_3 * fp_9[k]
                  + f_4 * fp_15[k];

        g_15[k] = f_5 * fp_7[k]
                  - f_5 * fp_22[k];

        g_16[k] = f_5 * fp_8[k]
                  - f_5 * fp_23[k];

        g_17[k] = f_5 * fp_6[k]
                  - f_5 * fp_21[k];

        g_18[k] = f_1 * fp_1[k]
                  - f_0 * fp_10[k];
    }

#pragma omp simd aligned(fp_0, fp_2, fp_9, fp_11 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = f_1 * fp_2[k]
                  - f_0 * fp_11[k];

        g_20[k] = f_1 * fp_0[k]
                  - f_0 * fp_9[k];
    }
}

}  // namespace simdtrf
