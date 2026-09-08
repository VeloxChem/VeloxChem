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


#include "SimdTransformGP.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_gp(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t gp,
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

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_1 = buffer.data(gp + 1);
    const auto *gp_2 = buffer.data(gp + 2);
    const auto *gp_3 = buffer.data(gp + 3);
    const auto *gp_4 = buffer.data(gp + 4);
    const auto *gp_5 = buffer.data(gp + 5);
    const auto *gp_6 = buffer.data(gp + 6);
    const auto *gp_7 = buffer.data(gp + 7);
    const auto *gp_8 = buffer.data(gp + 8);
    const auto *gp_9 = buffer.data(gp + 9);
    const auto *gp_10 = buffer.data(gp + 10);
    const auto *gp_11 = buffer.data(gp + 11);
    const auto *gp_12 = buffer.data(gp + 12);
    const auto *gp_13 = buffer.data(gp + 13);
    const auto *gp_14 = buffer.data(gp + 14);
    const auto *gp_15 = buffer.data(gp + 15);
    const auto *gp_16 = buffer.data(gp + 16);
    const auto *gp_17 = buffer.data(gp + 17);
    const auto *gp_18 = buffer.data(gp + 18);
    const auto *gp_19 = buffer.data(gp + 19);
    const auto *gp_20 = buffer.data(gp + 20);
    const auto *gp_21 = buffer.data(gp + 21);
    const auto *gp_22 = buffer.data(gp + 22);
    const auto *gp_23 = buffer.data(gp + 23);
    const auto *gp_24 = buffer.data(gp + 24);
    const auto *gp_25 = buffer.data(gp + 25);
    const auto *gp_26 = buffer.data(gp + 26);
    const auto *gp_27 = buffer.data(gp + 27);
    const auto *gp_28 = buffer.data(gp + 28);
    const auto *gp_29 = buffer.data(gp + 29);
    const auto *gp_30 = buffer.data(gp + 30);
    const auto *gp_31 = buffer.data(gp + 31);
    const auto *gp_32 = buffer.data(gp + 32);
    const auto *gp_33 = buffer.data(gp + 33);
    const auto *gp_34 = buffer.data(gp + 34);
    const auto *gp_35 = buffer.data(gp + 35);
    const auto *gp_36 = buffer.data(gp + 36);
    const auto *gp_37 = buffer.data(gp + 37);
    const auto *gp_38 = buffer.data(gp + 38);
    const auto *gp_39 = buffer.data(gp + 39);
    const auto *gp_40 = buffer.data(gp + 40);
    const auto *gp_41 = buffer.data(gp + 41);
    const auto *gp_42 = buffer.data(gp + 42);
    const auto *gp_43 = buffer.data(gp + 43);
    const auto *gp_44 = buffer.data(gp + 44);

#pragma omp simd aligned(gp_3, gp_4, gp_5, gp_13, gp_14, gp_18, gp_19, gp_20, gp_34, \
                         gp_35 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * gp_4[k]
                 - f_0 * gp_19[k];

        g_1[k] = f_0 * gp_5[k]
                 - f_0 * gp_20[k];

        g_2[k] = f_0 * gp_3[k]
                 - f_0 * gp_18[k];

        g_3[k] = f_1 * gp_13[k]
                 - f_2 * gp_34[k];

        g_4[k] = f_1 * gp_14[k]
                 - f_2 * gp_35[k];
    }

#pragma omp simd aligned(gp_3, gp_4, gp_5, gp_12, gp_18, gp_19, gp_20, gp_24, gp_25, gp_26, \
                         gp_33 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = f_1 * gp_12[k]
                 - f_2 * gp_33[k];

        g_6[k] = -f_3 * gp_4[k]
                 - f_3 * gp_19[k]
                 + f_4 * gp_25[k];

        g_7[k] = -f_3 * gp_5[k]
                 - f_3 * gp_20[k]
                 + f_4 * gp_26[k];

        g_8[k] = -f_3 * gp_3[k]
                 - f_3 * gp_18[k]
                 + f_4 * gp_24[k];
    }

#pragma omp simd aligned(gp_12, gp_13, gp_14, gp_33, gp_34, gp_35, gp_39, gp_40, \
                         gp_41 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_9[k] = -f_5 * gp_13[k]
                 - f_5 * gp_34[k]
                 + f_6 * gp_40[k];

        g_10[k] = -f_5 * gp_14[k]
                  - f_5 * gp_35[k]
                  + f_6 * gp_41[k];

        g_11[k] = -f_5 * gp_12[k]
                  - f_5 * gp_33[k]
                  + f_6 * gp_39[k];
    }

#pragma omp simd aligned(gp_1, gp_2, gp_10, gp_11, gp_16, gp_17, gp_31, gp_32, gp_37, gp_38, \
                         gp_43, gp_44 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_12[k] = 0.375 * gp_1[k]
                  + 0.75 * gp_10[k]
                  - 3.0 * gp_16[k]
                  + 0.375 * gp_31[k]
                  - 3.0 * gp_37[k]
                  + gp_43[k];

        g_13[k] = 0.375 * gp_2[k]
                  + 0.75 * gp_11[k]
                  - 3.0 * gp_17[k]
                  + 0.375 * gp_32[k]
                  - 3.0 * gp_38[k]
                  + gp_44[k];
    }

#pragma omp simd aligned(gp_0, gp_7, gp_8, gp_9, gp_15, gp_22, gp_23, gp_28, gp_29, gp_30, \
                         gp_36, gp_42 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = 0.375 * gp_0[k]
                  + 0.75 * gp_9[k]
                  - 3.0 * gp_15[k]
                  + 0.375 * gp_30[k]
                  - 3.0 * gp_36[k]
                  + gp_42[k];

        g_15[k] = -f_5 * gp_7[k]
                  - f_5 * gp_22[k]
                  + f_6 * gp_28[k];

        g_16[k] = -f_5 * gp_8[k]
                  - f_5 * gp_23[k]
                  + f_6 * gp_29[k];
    }

#pragma omp simd aligned(gp_1, gp_2, gp_6, gp_16, gp_17, gp_21, gp_27, gp_31, gp_32, gp_37, \
                         gp_38 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = -f_5 * gp_6[k]
                  - f_5 * gp_21[k]
                  + f_6 * gp_27[k];

        g_18[k] = -f_7 * gp_1[k]
                  + f_8 * gp_16[k]
                  + f_7 * gp_31[k]
                  - f_8 * gp_37[k];

        g_19[k] = -f_7 * gp_2[k]
                  + f_8 * gp_17[k]
                  + f_7 * gp_32[k]
                  - f_8 * gp_38[k];
    }

#pragma omp simd aligned(gp_0, gp_6, gp_7, gp_8, gp_15, gp_21, gp_22, gp_23, gp_30, \
                         gp_36 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = -f_7 * gp_0[k]
                  + f_8 * gp_15[k]
                  + f_7 * gp_30[k]
                  - f_8 * gp_36[k];

        g_21[k] = f_2 * gp_7[k]
                  - f_1 * gp_22[k];

        g_22[k] = f_2 * gp_8[k]
                  - f_1 * gp_23[k];

        g_23[k] = f_2 * gp_6[k]
                  - f_1 * gp_21[k];
    }

#pragma omp simd aligned(gp_0, gp_1, gp_2, gp_9, gp_10, gp_11, gp_30, gp_31, \
                         gp_32 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = f_9 * gp_1[k]
                  - f_10 * gp_10[k]
                  + f_9 * gp_31[k];

        g_25[k] = f_9 * gp_2[k]
                  - f_10 * gp_11[k]
                  + f_9 * gp_32[k];

        g_26[k] = f_9 * gp_0[k]
                  - f_10 * gp_9[k]
                  + f_9 * gp_30[k];
    }
}

}  // namespace simdtrf
