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


#include "SimdTransformHP.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_hp(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t hp,
             const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 0.9375 * std::sqrt(14.0);
    const auto f_1 = 1.875 * std::sqrt(14.0);
    const auto f_2 = 0.1875 * std::sqrt(14.0);
    const auto f_3 = 1.5 * std::sqrt(35.0);
    const auto f_4 = 0.1875 * std::sqrt(70.0);
    const auto f_5 = 0.125 * std::sqrt(70.0);
    const auto f_6 = 1.5 * std::sqrt(70.0);
    const auto f_7 = 0.0625 * std::sqrt(70.0);
    const auto f_8 = 0.5 * std::sqrt(70.0);
    const auto f_9 = 0.5 * std::sqrt(105.0);
    const auto f_10 = std::sqrt(105.0);
    const auto f_11 = 0.125 * std::sqrt(15.0);
    const auto f_12 = 0.25 * std::sqrt(15.0);
    const auto f_13 = 1.5 * std::sqrt(15.0);
    const auto f_14 = std::sqrt(15.0);
    const auto f_15 = 0.25 * std::sqrt(105.0);
    const auto f_16 = 0.375 * std::sqrt(35.0);
    const auto f_17 = 2.25 * std::sqrt(35.0);

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

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);
    const auto *hp_3 = buffer.data(hp + 3);
    const auto *hp_4 = buffer.data(hp + 4);
    const auto *hp_5 = buffer.data(hp + 5);
    const auto *hp_6 = buffer.data(hp + 6);
    const auto *hp_7 = buffer.data(hp + 7);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_9 = buffer.data(hp + 9);
    const auto *hp_10 = buffer.data(hp + 10);
    const auto *hp_11 = buffer.data(hp + 11);
    const auto *hp_12 = buffer.data(hp + 12);
    const auto *hp_13 = buffer.data(hp + 13);
    const auto *hp_14 = buffer.data(hp + 14);
    const auto *hp_15 = buffer.data(hp + 15);
    const auto *hp_16 = buffer.data(hp + 16);
    const auto *hp_17 = buffer.data(hp + 17);
    const auto *hp_18 = buffer.data(hp + 18);
    const auto *hp_19 = buffer.data(hp + 19);
    const auto *hp_20 = buffer.data(hp + 20);
    const auto *hp_21 = buffer.data(hp + 21);
    const auto *hp_22 = buffer.data(hp + 22);
    const auto *hp_23 = buffer.data(hp + 23);
    const auto *hp_24 = buffer.data(hp + 24);
    const auto *hp_25 = buffer.data(hp + 25);
    const auto *hp_26 = buffer.data(hp + 26);
    const auto *hp_27 = buffer.data(hp + 27);
    const auto *hp_28 = buffer.data(hp + 28);
    const auto *hp_29 = buffer.data(hp + 29);
    const auto *hp_30 = buffer.data(hp + 30);
    const auto *hp_31 = buffer.data(hp + 31);
    const auto *hp_32 = buffer.data(hp + 32);
    const auto *hp_33 = buffer.data(hp + 33);
    const auto *hp_34 = buffer.data(hp + 34);
    const auto *hp_35 = buffer.data(hp + 35);
    const auto *hp_36 = buffer.data(hp + 36);
    const auto *hp_37 = buffer.data(hp + 37);
    const auto *hp_38 = buffer.data(hp + 38);
    const auto *hp_39 = buffer.data(hp + 39);
    const auto *hp_40 = buffer.data(hp + 40);
    const auto *hp_41 = buffer.data(hp + 41);
    const auto *hp_42 = buffer.data(hp + 42);
    const auto *hp_43 = buffer.data(hp + 43);
    const auto *hp_44 = buffer.data(hp + 44);
    const auto *hp_45 = buffer.data(hp + 45);
    const auto *hp_46 = buffer.data(hp + 46);
    const auto *hp_47 = buffer.data(hp + 47);
    const auto *hp_48 = buffer.data(hp + 48);
    const auto *hp_49 = buffer.data(hp + 49);
    const auto *hp_50 = buffer.data(hp + 50);
    const auto *hp_51 = buffer.data(hp + 51);
    const auto *hp_52 = buffer.data(hp + 52);
    const auto *hp_53 = buffer.data(hp + 53);
    const auto *hp_54 = buffer.data(hp + 54);
    const auto *hp_55 = buffer.data(hp + 55);
    const auto *hp_56 = buffer.data(hp + 56);
    const auto *hp_57 = buffer.data(hp + 57);
    const auto *hp_58 = buffer.data(hp + 58);
    const auto *hp_59 = buffer.data(hp + 59);
    const auto *hp_60 = buffer.data(hp + 60);
    const auto *hp_61 = buffer.data(hp + 61);
    const auto *hp_62 = buffer.data(hp + 62);

#pragma omp simd aligned(hp_3, hp_4, hp_5, hp_13, hp_18, hp_19, hp_20, hp_34, hp_45, hp_46, \
                         hp_47 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * hp_4[k]
                 - f_1 * hp_19[k]
                 + f_2 * hp_46[k];

        g_1[k] = f_0 * hp_5[k]
                 - f_1 * hp_20[k]
                 + f_2 * hp_47[k];

        g_2[k] = f_0 * hp_3[k]
                 - f_1 * hp_18[k]
                 + f_2 * hp_45[k];

        g_3[k] = f_3 * hp_13[k]
                 - f_3 * hp_34[k];
    }

#pragma omp simd aligned(hp_4, hp_12, hp_14, hp_19, hp_25, hp_33, hp_35, hp_46, \
                         hp_52 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = f_3 * hp_14[k]
                 - f_3 * hp_35[k];

        g_5[k] = f_3 * hp_12[k]
                 - f_3 * hp_33[k];

        g_6[k] = -f_4 * hp_4[k]
                 - f_5 * hp_19[k]
                 + f_6 * hp_25[k]
                 + f_7 * hp_46[k]
                 - f_8 * hp_52[k];
    }

#pragma omp simd aligned(hp_3, hp_5, hp_13, hp_18, hp_20, hp_24, hp_26, hp_34, hp_40, hp_45, \
                         hp_47, hp_51, hp_53 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = -f_4 * hp_5[k]
                 - f_5 * hp_20[k]
                 + f_6 * hp_26[k]
                 + f_7 * hp_47[k]
                 - f_8 * hp_53[k];

        g_8[k] = -f_4 * hp_3[k]
                 - f_5 * hp_18[k]
                 + f_6 * hp_24[k]
                 + f_7 * hp_45[k]
                 - f_8 * hp_51[k];

        g_9[k] = -f_9 * hp_13[k]
                 - f_9 * hp_34[k]
                 + f_10 * hp_40[k];
    }

#pragma omp simd aligned(hp_4, hp_12, hp_14, hp_19, hp_25, hp_33, hp_35, hp_39, hp_41, hp_46, \
                         hp_52, hp_58 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = -f_9 * hp_14[k]
                  - f_9 * hp_35[k]
                  + f_10 * hp_41[k];

        g_11[k] = -f_9 * hp_12[k]
                  - f_9 * hp_33[k]
                  + f_10 * hp_39[k];

        g_12[k] = f_11 * hp_4[k]
                  + f_12 * hp_19[k]
                  - f_13 * hp_25[k]
                  + f_11 * hp_46[k]
                  - f_13 * hp_52[k]
                  + f_14 * hp_58[k];
    }

#pragma omp simd aligned(hp_3, hp_5, hp_18, hp_20, hp_24, hp_26, hp_45, hp_47, hp_51, hp_53, \
                         hp_57, hp_59 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = f_11 * hp_5[k]
                  + f_12 * hp_20[k]
                  - f_13 * hp_26[k]
                  + f_11 * hp_47[k]
                  - f_13 * hp_53[k]
                  + f_14 * hp_59[k];

        g_14[k] = f_11 * hp_3[k]
                  + f_12 * hp_18[k]
                  - f_13 * hp_24[k]
                  + f_11 * hp_45[k]
                  - f_13 * hp_51[k]
                  + f_14 * hp_57[k];
    }

#pragma omp simd aligned(hp_7, hp_8, hp_22, hp_23, hp_28, hp_29, hp_49, hp_50, hp_55, hp_56, \
                         hp_61, hp_62 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_15[k] = 1.875 * hp_7[k]
                  + 3.75 * hp_22[k]
                  - 5.0 * hp_28[k]
                  + 1.875 * hp_49[k]
                  - 5.0 * hp_55[k]
                  + hp_61[k];

        g_16[k] = 1.875 * hp_8[k]
                  + 3.75 * hp_23[k]
                  - 5.0 * hp_29[k]
                  + 1.875 * hp_50[k]
                  - 5.0 * hp_56[k]
                  + hp_62[k];
    }

#pragma omp simd aligned(hp_1, hp_6, hp_10, hp_16, hp_21, hp_27, hp_31, hp_37, hp_43, hp_48, \
                         hp_54, hp_60 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = 1.875 * hp_6[k]
                  + 3.75 * hp_21[k]
                  - 5.0 * hp_27[k]
                  + 1.875 * hp_48[k]
                  - 5.0 * hp_54[k]
                  + hp_60[k];

        g_18[k] = f_11 * hp_1[k]
                  + f_12 * hp_10[k]
                  - f_13 * hp_16[k]
                  + f_11 * hp_31[k]
                  - f_13 * hp_37[k]
                  + f_14 * hp_43[k];
    }

#pragma omp simd aligned(hp_0, hp_2, hp_9, hp_11, hp_15, hp_17, hp_30, hp_32, hp_36, hp_38, \
                         hp_42, hp_44 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = f_11 * hp_2[k]
                  + f_12 * hp_11[k]
                  - f_13 * hp_17[k]
                  + f_11 * hp_32[k]
                  - f_13 * hp_38[k]
                  + f_14 * hp_44[k];

        g_20[k] = f_11 * hp_0[k]
                  + f_12 * hp_9[k]
                  - f_13 * hp_15[k]
                  + f_11 * hp_30[k]
                  - f_13 * hp_36[k]
                  + f_14 * hp_42[k];
    }

#pragma omp simd aligned(hp_6, hp_7, hp_8, hp_27, hp_28, hp_29, hp_48, hp_49, hp_50, hp_54, \
                         hp_55, hp_56 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = -f_15 * hp_7[k]
                  + f_9 * hp_28[k]
                  + f_15 * hp_49[k]
                  - f_9 * hp_55[k];

        g_22[k] = -f_15 * hp_8[k]
                  + f_9 * hp_29[k]
                  + f_15 * hp_50[k]
                  - f_9 * hp_56[k];

        g_23[k] = -f_15 * hp_6[k]
                  + f_9 * hp_27[k]
                  + f_15 * hp_48[k]
                  - f_9 * hp_54[k];
    }

#pragma omp simd aligned(hp_1, hp_2, hp_10, hp_11, hp_16, hp_17, hp_31, hp_32, hp_37, \
                         hp_38 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = -f_7 * hp_1[k]
                  + f_5 * hp_10[k]
                  + f_8 * hp_16[k]
                  + f_4 * hp_31[k]
                  - f_6 * hp_37[k];

        g_25[k] = -f_7 * hp_2[k]
                  + f_5 * hp_11[k]
                  + f_8 * hp_17[k]
                  + f_4 * hp_32[k]
                  - f_6 * hp_38[k];
    }

#pragma omp simd aligned(hp_0, hp_7, hp_8, hp_9, hp_15, hp_22, hp_23, hp_30, hp_36, hp_49, \
                         hp_50 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = -f_7 * hp_0[k]
                  + f_5 * hp_9[k]
                  + f_8 * hp_15[k]
                  + f_4 * hp_30[k]
                  - f_6 * hp_36[k];

        g_27[k] = f_16 * hp_7[k]
                  - f_17 * hp_22[k]
                  + f_16 * hp_49[k];

        g_28[k] = f_16 * hp_8[k]
                  - f_17 * hp_23[k]
                  + f_16 * hp_50[k];
    }

#pragma omp simd aligned(hp_0, hp_1, hp_2, hp_6, hp_9, hp_10, hp_11, hp_21, hp_30, hp_31, \
                         hp_32, hp_48 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = f_16 * hp_6[k]
                  - f_17 * hp_21[k]
                  + f_16 * hp_48[k];

        g_30[k] = f_2 * hp_1[k]
                  - f_1 * hp_10[k]
                  + f_0 * hp_31[k];

        g_31[k] = f_2 * hp_2[k]
                  - f_1 * hp_11[k]
                  + f_0 * hp_32[k];

        g_32[k] = f_2 * hp_0[k]
                  - f_1 * hp_9[k]
                  + f_0 * hp_30[k];
    }
}

}  // namespace simdtrf
