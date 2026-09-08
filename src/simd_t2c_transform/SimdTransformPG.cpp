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


#include "SimdTransformPG.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_pg(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t pg,
             const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 0.5 * std::sqrt(35.0);
    const auto f_1 = 0.75 * std::sqrt(70.0);
    const auto f_2 = 0.25 * std::sqrt(70.0);
    const auto f_3 = 0.5 * std::sqrt(5.0);
    const auto f_4 = 3.0 * std::sqrt(5.0);
    const auto f_5 = 0.75 * std::sqrt(10.0);
    const auto f_6 = std::sqrt(10.0);
    const auto f_7 = 0.25 * std::sqrt(5.0);
    const auto f_8 = 1.5 * std::sqrt(5.0);
    const auto f_9 = 0.125 * std::sqrt(35.0);
    const auto f_10 = 0.75 * std::sqrt(35.0);

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
    auto *g_21 = values + 21 * nvalues;
    auto *g_22 = values + 22 * nvalues;
    auto *g_23 = values + 23 * nvalues;
    auto *g_24 = values + 24 * nvalues;
    auto *g_25 = values + 25 * nvalues;
    auto *g_26 = values + 26 * nvalues;

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_2 = buffer.data(pg + 2);
    const auto *pg_3 = buffer.data(pg + 3);
    const auto *pg_4 = buffer.data(pg + 4);
    const auto *pg_5 = buffer.data(pg + 5);
    const auto *pg_6 = buffer.data(pg + 6);
    const auto *pg_7 = buffer.data(pg + 7);
    const auto *pg_8 = buffer.data(pg + 8);
    const auto *pg_9 = buffer.data(pg + 9);
    const auto *pg_10 = buffer.data(pg + 10);
    const auto *pg_11 = buffer.data(pg + 11);
    const auto *pg_12 = buffer.data(pg + 12);
    const auto *pg_13 = buffer.data(pg + 13);
    const auto *pg_14 = buffer.data(pg + 14);
    const auto *pg_15 = buffer.data(pg + 15);
    const auto *pg_16 = buffer.data(pg + 16);
    const auto *pg_17 = buffer.data(pg + 17);
    const auto *pg_18 = buffer.data(pg + 18);
    const auto *pg_19 = buffer.data(pg + 19);
    const auto *pg_20 = buffer.data(pg + 20);
    const auto *pg_21 = buffer.data(pg + 21);
    const auto *pg_22 = buffer.data(pg + 22);
    const auto *pg_23 = buffer.data(pg + 23);
    const auto *pg_24 = buffer.data(pg + 24);
    const auto *pg_25 = buffer.data(pg + 25);
    const auto *pg_26 = buffer.data(pg + 26);
    const auto *pg_27 = buffer.data(pg + 27);
    const auto *pg_28 = buffer.data(pg + 28);
    const auto *pg_29 = buffer.data(pg + 29);
    const auto *pg_30 = buffer.data(pg + 30);
    const auto *pg_31 = buffer.data(pg + 31);
    const auto *pg_32 = buffer.data(pg + 32);
    const auto *pg_33 = buffer.data(pg + 33);
    const auto *pg_34 = buffer.data(pg + 34);
    const auto *pg_35 = buffer.data(pg + 35);
    const auto *pg_36 = buffer.data(pg + 36);
    const auto *pg_37 = buffer.data(pg + 37);
    const auto *pg_38 = buffer.data(pg + 38);
    const auto *pg_39 = buffer.data(pg + 39);
    const auto *pg_40 = buffer.data(pg + 40);
    const auto *pg_41 = buffer.data(pg + 41);
    const auto *pg_42 = buffer.data(pg + 42);
    const auto *pg_43 = buffer.data(pg + 43);
    const auto *pg_44 = buffer.data(pg + 44);

#pragma omp simd aligned(pg_16, pg_19, pg_21, pg_23, pg_26, pg_28 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * pg_16[k]
                 - f_0 * pg_21[k];

        g_1[k] = f_1 * pg_19[k]
                 - f_2 * pg_26[k];

        g_2[k] = -f_3 * pg_16[k]
                 - f_3 * pg_21[k]
                 + f_4 * pg_23[k];

        g_3[k] = -f_5 * pg_19[k]
                 - f_5 * pg_26[k]
                 + f_6 * pg_28[k];
    }

#pragma omp simd aligned(pg_15, pg_17, pg_18, pg_20, pg_22, pg_24, pg_25, pg_27, \
                         pg_29 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = 0.375 * pg_15[k]
                 + 0.75 * pg_18[k]
                 - 3.0 * pg_20[k]
                 + 0.375 * pg_25[k]
                 - 3.0 * pg_27[k]
                 + pg_29[k];

        g_5[k] = -f_5 * pg_17[k]
                 - f_5 * pg_22[k]
                 + f_6 * pg_24[k];

        g_6[k] = -f_7 * pg_15[k]
                 + f_8 * pg_20[k]
                 + f_7 * pg_25[k]
                 - f_8 * pg_27[k];

        g_7[k] = f_2 * pg_17[k]
                 - f_1 * pg_22[k];

        g_8[k] = f_9 * pg_15[k]
                 - f_10 * pg_18[k]
                 + f_9 * pg_25[k];
    }

#pragma omp simd aligned(pg_31, pg_34, pg_36, pg_38, pg_41, pg_43 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_9[k] = f_0 * pg_31[k]
                 - f_0 * pg_36[k];

        g_10[k] = f_1 * pg_34[k]
                  - f_2 * pg_41[k];

        g_11[k] = -f_3 * pg_31[k]
                  - f_3 * pg_36[k]
                  + f_4 * pg_38[k];

        g_12[k] = -f_5 * pg_34[k]
                  - f_5 * pg_41[k]
                  + f_6 * pg_43[k];
    }

#pragma omp simd aligned(pg_30, pg_32, pg_33, pg_35, pg_37, pg_39, pg_40, pg_42, \
                         pg_44 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = 0.375 * pg_30[k]
                  + 0.75 * pg_33[k]
                  - 3.0 * pg_35[k]
                  + 0.375 * pg_40[k]
                  - 3.0 * pg_42[k]
                  + pg_44[k];

        g_14[k] = -f_5 * pg_32[k]
                  - f_5 * pg_37[k]
                  + f_6 * pg_39[k];

        g_15[k] = -f_7 * pg_30[k]
                  + f_8 * pg_35[k]
                  + f_7 * pg_40[k]
                  - f_8 * pg_42[k];

        g_16[k] = f_2 * pg_32[k]
                  - f_1 * pg_37[k];

        g_17[k] = f_9 * pg_30[k]
                  - f_10 * pg_33[k]
                  + f_9 * pg_40[k];
    }

#pragma omp simd aligned(pg_1, pg_4, pg_6, pg_8, pg_11, pg_13 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = f_0 * pg_1[k]
                  - f_0 * pg_6[k];

        g_19[k] = f_1 * pg_4[k]
                  - f_2 * pg_11[k];

        g_20[k] = -f_3 * pg_1[k]
                  - f_3 * pg_6[k]
                  + f_4 * pg_8[k];

        g_21[k] = -f_5 * pg_4[k]
                  - f_5 * pg_11[k]
                  + f_6 * pg_13[k];
    }

#pragma omp simd aligned(pg_0, pg_2, pg_3, pg_5, pg_7, pg_9, pg_10, pg_12, \
                         pg_14 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = 0.375 * pg_0[k]
                  + 0.75 * pg_3[k]
                  - 3.0 * pg_5[k]
                  + 0.375 * pg_10[k]
                  - 3.0 * pg_12[k]
                  + pg_14[k];

        g_23[k] = -f_5 * pg_2[k]
                  - f_5 * pg_7[k]
                  + f_6 * pg_9[k];

        g_24[k] = -f_7 * pg_0[k]
                  + f_8 * pg_5[k]
                  + f_7 * pg_10[k]
                  - f_8 * pg_12[k];

        g_25[k] = f_2 * pg_2[k]
                  - f_1 * pg_7[k];

        g_26[k] = f_9 * pg_0[k]
                  - f_10 * pg_3[k]
                  + f_9 * pg_10[k];
    }
}

}  // namespace simdtrf
