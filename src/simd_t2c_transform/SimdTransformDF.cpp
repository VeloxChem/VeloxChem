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


#include "SimdTransformDF.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_df(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t df,
             const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 0.75 * std::sqrt(30.0);
    const auto f_1 = 0.25 * std::sqrt(30.0);
    const auto f_2 = 3.0 * std::sqrt(5.0);
    const auto f_3 = 0.75 * std::sqrt(2.0);
    const auto f_4 = 3.0 * std::sqrt(2.0);
    const auto f_5 = 1.5 * std::sqrt(3.0);
    const auto f_6 = std::sqrt(3.0);
    const auto f_7 = 1.5 * std::sqrt(5.0);
    const auto f_8 = 0.375 * std::sqrt(10.0);
    const auto f_9 = 0.125 * std::sqrt(10.0);
    const auto f_10 = 0.75 * std::sqrt(10.0);
    const auto f_11 = 0.25 * std::sqrt(10.0);
    const auto f_12 = 0.5 * std::sqrt(15.0);
    const auto f_13 = std::sqrt(15.0);
    const auto f_14 = 0.125 * std::sqrt(6.0);
    const auto f_15 = 0.5 * std::sqrt(6.0);
    const auto f_16 = 0.25 * std::sqrt(6.0);
    const auto f_17 = std::sqrt(6.0);
    const auto f_18 = 0.25 * std::sqrt(15.0);
    const auto f_19 = 0.375 * std::sqrt(30.0);
    const auto f_20 = 0.125 * std::sqrt(30.0);
    const auto f_21 = 0.375 * std::sqrt(2.0);
    const auto f_22 = 1.5 * std::sqrt(2.0);
    const auto f_23 = 0.75 * std::sqrt(3.0);
    const auto f_24 = 0.5 * std::sqrt(3.0);
    const auto f_25 = 0.75 * std::sqrt(5.0);

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
    auto *g_27 = values + 27 * nvalues;
    auto *g_28 = values + 28 * nvalues;
    auto *g_29 = values + 29 * nvalues;
    auto *g_30 = values + 30 * nvalues;
    auto *g_31 = values + 31 * nvalues;
    auto *g_32 = values + 32 * nvalues;
    auto *g_33 = values + 33 * nvalues;
    auto *g_34 = values + 34 * nvalues;

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_7 = buffer.data(df + 7);
    const auto *df_8 = buffer.data(df + 8);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_10 = buffer.data(df + 10);
    const auto *df_11 = buffer.data(df + 11);
    const auto *df_12 = buffer.data(df + 12);
    const auto *df_13 = buffer.data(df + 13);
    const auto *df_14 = buffer.data(df + 14);
    const auto *df_15 = buffer.data(df + 15);
    const auto *df_16 = buffer.data(df + 16);
    const auto *df_17 = buffer.data(df + 17);
    const auto *df_18 = buffer.data(df + 18);
    const auto *df_19 = buffer.data(df + 19);
    const auto *df_20 = buffer.data(df + 20);
    const auto *df_21 = buffer.data(df + 21);
    const auto *df_22 = buffer.data(df + 22);
    const auto *df_23 = buffer.data(df + 23);
    const auto *df_24 = buffer.data(df + 24);
    const auto *df_25 = buffer.data(df + 25);
    const auto *df_26 = buffer.data(df + 26);
    const auto *df_27 = buffer.data(df + 27);
    const auto *df_28 = buffer.data(df + 28);
    const auto *df_29 = buffer.data(df + 29);
    const auto *df_30 = buffer.data(df + 30);
    const auto *df_31 = buffer.data(df + 31);
    const auto *df_32 = buffer.data(df + 32);
    const auto *df_33 = buffer.data(df + 33);
    const auto *df_34 = buffer.data(df + 34);
    const auto *df_35 = buffer.data(df + 35);
    const auto *df_36 = buffer.data(df + 36);
    const auto *df_37 = buffer.data(df + 37);
    const auto *df_38 = buffer.data(df + 38);
    const auto *df_39 = buffer.data(df + 39);
    const auto *df_40 = buffer.data(df + 40);
    const auto *df_41 = buffer.data(df + 41);
    const auto *df_42 = buffer.data(df + 42);
    const auto *df_43 = buffer.data(df + 43);
    const auto *df_44 = buffer.data(df + 44);
    const auto *df_45 = buffer.data(df + 45);
    const auto *df_46 = buffer.data(df + 46);
    const auto *df_47 = buffer.data(df + 47);
    const auto *df_48 = buffer.data(df + 48);
    const auto *df_49 = buffer.data(df + 49);
    const auto *df_50 = buffer.data(df + 50);
    const auto *df_51 = buffer.data(df + 51);
    const auto *df_52 = buffer.data(df + 52);
    const auto *df_53 = buffer.data(df + 53);
    const auto *df_54 = buffer.data(df + 54);
    const auto *df_55 = buffer.data(df + 55);
    const auto *df_56 = buffer.data(df + 56);
    const auto *df_57 = buffer.data(df + 57);
    const auto *df_58 = buffer.data(df + 58);
    const auto *df_59 = buffer.data(df + 59);

#pragma omp simd aligned(df_10, df_11, df_12, df_13, df_14, df_15, df_16, df_17, df_18, \
                         df_19 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * df_11[k]
                 - f_1 * df_16[k];

        g_1[k] = f_2 * df_14[k];

        g_2[k] = -f_3 * df_11[k]
                 - f_3 * df_16[k]
                 + f_4 * df_18[k];

        g_3[k] = -f_5 * df_12[k]
                 - f_5 * df_17[k]
                 + f_6 * df_19[k];

        g_4[k] = -f_3 * df_10[k]
                 - f_3 * df_13[k]
                 + f_4 * df_15[k];

        g_5[k] = f_7 * df_12[k]
                 - f_7 * df_17[k];
    }

#pragma omp simd aligned(df_10, df_13, df_41, df_42, df_44, df_46, df_47, df_48, \
                         df_49 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = f_1 * df_10[k]
                 - f_0 * df_13[k];

        g_7[k] = f_0 * df_41[k]
                 - f_1 * df_46[k];

        g_8[k] = f_2 * df_44[k];

        g_9[k] = -f_3 * df_41[k]
                 - f_3 * df_46[k]
                 + f_4 * df_48[k];

        g_10[k] = -f_5 * df_42[k]
                  - f_5 * df_47[k]
                  + f_6 * df_49[k];
    }

#pragma omp simd aligned(df_1, df_6, df_31, df_36, df_40, df_42, df_43, df_45, df_47, df_51, \
                         df_56 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_11[k] = -f_3 * df_40[k]
                  - f_3 * df_43[k]
                  + f_4 * df_45[k];

        g_12[k] = f_7 * df_42[k]
                  - f_7 * df_47[k];

        g_13[k] = f_1 * df_40[k]
                  - f_0 * df_43[k];

        g_14[k] = -f_8 * df_1[k]
                  + f_9 * df_6[k]
                  - f_8 * df_31[k]
                  + f_9 * df_36[k]
                  + f_10 * df_51[k]
                  - f_11 * df_56[k];
    }

#pragma omp simd aligned(df_1, df_4, df_6, df_8, df_31, df_34, df_36, df_38, df_51, df_54, \
                         df_56, df_58 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_15[k] = -f_12 * df_4[k]
                  - f_12 * df_34[k]
                  + f_13 * df_54[k];

        g_16[k] = f_14 * df_1[k]
                  + f_14 * df_6[k]
                  - f_15 * df_8[k]
                  + f_14 * df_31[k]
                  + f_14 * df_36[k]
                  - f_15 * df_38[k]
                  - f_16 * df_51[k]
                  - f_16 * df_56[k]
                  + f_17 * df_58[k];
    }

#pragma omp simd aligned(df_2, df_7, df_9, df_32, df_37, df_39, df_52, df_57, \
                         df_59 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = 0.75 * df_2[k]
                  + 0.75 * df_7[k]
                  - 0.5 * df_9[k]
                  + 0.75 * df_32[k]
                  + 0.75 * df_37[k]
                  - 0.5 * df_39[k]
                  - 1.5 * df_52[k]
                  - 1.5 * df_57[k]
                  + df_59[k];
    }

#pragma omp simd aligned(df_0, df_3, df_5, df_30, df_33, df_35, df_50, df_53, \
                         df_55 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = f_14 * df_0[k]
                  + f_14 * df_3[k]
                  - f_15 * df_5[k]
                  + f_14 * df_30[k]
                  + f_14 * df_33[k]
                  - f_15 * df_35[k]
                  - f_16 * df_50[k]
                  - f_16 * df_53[k]
                  + f_17 * df_55[k];
    }

#pragma omp simd aligned(df_0, df_2, df_3, df_7, df_30, df_32, df_33, df_37, df_50, df_52, \
                         df_53, df_57 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = -f_18 * df_2[k]
                  + f_18 * df_7[k]
                  - f_18 * df_32[k]
                  + f_18 * df_37[k]
                  + f_12 * df_52[k]
                  - f_12 * df_57[k];

        g_20[k] = -f_9 * df_0[k]
                  + f_8 * df_3[k]
                  - f_9 * df_30[k]
                  + f_8 * df_33[k]
                  + f_11 * df_50[k]
                  - f_10 * df_53[k];
    }

#pragma omp simd aligned(df_20, df_21, df_22, df_23, df_24, df_25, df_26, df_27, df_28, \
                         df_29 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = f_0 * df_21[k]
                  - f_1 * df_26[k];

        g_22[k] = f_2 * df_24[k];

        g_23[k] = -f_3 * df_21[k]
                  - f_3 * df_26[k]
                  + f_4 * df_28[k];

        g_24[k] = -f_5 * df_22[k]
                  - f_5 * df_27[k]
                  + f_6 * df_29[k];

        g_25[k] = -f_3 * df_20[k]
                  - f_3 * df_23[k]
                  + f_4 * df_25[k];

        g_26[k] = f_7 * df_22[k]
                  - f_7 * df_27[k];
    }

#pragma omp simd aligned(df_1, df_4, df_6, df_8, df_20, df_23, df_31, df_34, df_36, \
                         df_38 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = f_1 * df_20[k]
                  - f_0 * df_23[k];

        g_28[k] = f_19 * df_1[k]
                  - f_20 * df_6[k]
                  - f_19 * df_31[k]
                  + f_20 * df_36[k];

        g_29[k] = f_7 * df_4[k]
                  - f_7 * df_34[k];

        g_30[k] = -f_21 * df_1[k]
                  - f_21 * df_6[k]
                  + f_22 * df_8[k]
                  + f_21 * df_31[k]
                  + f_21 * df_36[k]
                  - f_22 * df_38[k];
    }

#pragma omp simd aligned(df_0, df_2, df_3, df_5, df_7, df_9, df_30, df_32, df_33, df_35, \
                         df_37, df_39 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = -f_23 * df_2[k]
                  - f_23 * df_7[k]
                  + f_24 * df_9[k]
                  + f_23 * df_32[k]
                  + f_23 * df_37[k]
                  - f_24 * df_39[k];

        g_32[k] = -f_21 * df_0[k]
                  - f_21 * df_3[k]
                  + f_22 * df_5[k]
                  + f_21 * df_30[k]
                  + f_21 * df_33[k]
                  - f_22 * df_35[k];

        g_33[k] = f_25 * df_2[k]
                  - f_25 * df_7[k]
                  - f_25 * df_32[k]
                  + f_25 * df_37[k];

        g_34[k] = f_20 * df_0[k]
                  - f_19 * df_3[k]
                  - f_20 * df_30[k]
                  + f_19 * df_33[k];
    }
}

}  // namespace simdtrf
