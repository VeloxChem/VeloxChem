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


#include "SimdTransformPF.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_pf(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t pf,
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

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);
    const auto *pf_3 = buffer.data(pf + 3);
    const auto *pf_4 = buffer.data(pf + 4);
    const auto *pf_5 = buffer.data(pf + 5);
    const auto *pf_6 = buffer.data(pf + 6);
    const auto *pf_7 = buffer.data(pf + 7);
    const auto *pf_8 = buffer.data(pf + 8);
    const auto *pf_9 = buffer.data(pf + 9);
    const auto *pf_10 = buffer.data(pf + 10);
    const auto *pf_11 = buffer.data(pf + 11);
    const auto *pf_12 = buffer.data(pf + 12);
    const auto *pf_13 = buffer.data(pf + 13);
    const auto *pf_14 = buffer.data(pf + 14);
    const auto *pf_15 = buffer.data(pf + 15);
    const auto *pf_16 = buffer.data(pf + 16);
    const auto *pf_17 = buffer.data(pf + 17);
    const auto *pf_18 = buffer.data(pf + 18);
    const auto *pf_19 = buffer.data(pf + 19);
    const auto *pf_20 = buffer.data(pf + 20);
    const auto *pf_21 = buffer.data(pf + 21);
    const auto *pf_22 = buffer.data(pf + 22);
    const auto *pf_23 = buffer.data(pf + 23);
    const auto *pf_24 = buffer.data(pf + 24);
    const auto *pf_25 = buffer.data(pf + 25);
    const auto *pf_26 = buffer.data(pf + 26);
    const auto *pf_27 = buffer.data(pf + 27);
    const auto *pf_28 = buffer.data(pf + 28);
    const auto *pf_29 = buffer.data(pf + 29);

#pragma omp simd aligned(pf_10, pf_11, pf_12, pf_13, pf_14, pf_15, pf_16, pf_17, pf_18, \
                         pf_19 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * pf_11[k]
                 - f_1 * pf_16[k];

        g_1[k] = f_2 * pf_14[k];

        g_2[k] = -f_3 * pf_11[k]
                 - f_3 * pf_16[k]
                 + f_4 * pf_18[k];

        g_3[k] = -1.5 * pf_12[k]
                 - 1.5 * pf_17[k]
                 + pf_19[k];

        g_4[k] = -f_3 * pf_10[k]
                 - f_3 * pf_13[k]
                 + f_4 * pf_15[k];

        g_5[k] = f_5 * pf_12[k]
                 - f_5 * pf_17[k];
    }

#pragma omp simd aligned(pf_10, pf_13, pf_21, pf_22, pf_24, pf_26, pf_27, pf_28, \
                         pf_29 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = f_1 * pf_10[k]
                 - f_0 * pf_13[k];

        g_7[k] = f_0 * pf_21[k]
                 - f_1 * pf_26[k];

        g_8[k] = f_2 * pf_24[k];

        g_9[k] = -f_3 * pf_21[k]
                 - f_3 * pf_26[k]
                 + f_4 * pf_28[k];

        g_10[k] = -1.5 * pf_22[k]
                  - 1.5 * pf_27[k]
                  + pf_29[k];
    }

#pragma omp simd aligned(pf_1, pf_4, pf_6, pf_8, pf_20, pf_22, pf_23, pf_25, \
                         pf_27 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_11[k] = -f_3 * pf_20[k]
                  - f_3 * pf_23[k]
                  + f_4 * pf_25[k];

        g_12[k] = f_5 * pf_22[k]
                  - f_5 * pf_27[k];

        g_13[k] = f_1 * pf_20[k]
                  - f_0 * pf_23[k];

        g_14[k] = f_0 * pf_1[k]
                  - f_1 * pf_6[k];

        g_15[k] = f_2 * pf_4[k];

        g_16[k] = -f_3 * pf_1[k]
                  - f_3 * pf_6[k]
                  + f_4 * pf_8[k];
    }

#pragma omp simd aligned(pf_0, pf_2, pf_3, pf_5, pf_7, pf_9 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = -1.5 * pf_2[k]
                  - 1.5 * pf_7[k]
                  + pf_9[k];

        g_18[k] = -f_3 * pf_0[k]
                  - f_3 * pf_3[k]
                  + f_4 * pf_5[k];

        g_19[k] = f_5 * pf_2[k]
                  - f_5 * pf_7[k];

        g_20[k] = f_1 * pf_0[k]
                  - f_0 * pf_3[k];
    }
}

}  // namespace simdtrf
