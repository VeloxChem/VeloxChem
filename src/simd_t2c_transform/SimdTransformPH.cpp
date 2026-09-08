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


#include "SimdTransformPH.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_ph(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t ph,
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

    const auto *ph_0 = buffer.data(ph + 0);
    const auto *ph_1 = buffer.data(ph + 1);
    const auto *ph_2 = buffer.data(ph + 2);
    const auto *ph_3 = buffer.data(ph + 3);
    const auto *ph_4 = buffer.data(ph + 4);
    const auto *ph_5 = buffer.data(ph + 5);
    const auto *ph_6 = buffer.data(ph + 6);
    const auto *ph_7 = buffer.data(ph + 7);
    const auto *ph_8 = buffer.data(ph + 8);
    const auto *ph_9 = buffer.data(ph + 9);
    const auto *ph_10 = buffer.data(ph + 10);
    const auto *ph_11 = buffer.data(ph + 11);
    const auto *ph_12 = buffer.data(ph + 12);
    const auto *ph_13 = buffer.data(ph + 13);
    const auto *ph_14 = buffer.data(ph + 14);
    const auto *ph_15 = buffer.data(ph + 15);
    const auto *ph_16 = buffer.data(ph + 16);
    const auto *ph_17 = buffer.data(ph + 17);
    const auto *ph_18 = buffer.data(ph + 18);
    const auto *ph_19 = buffer.data(ph + 19);
    const auto *ph_20 = buffer.data(ph + 20);
    const auto *ph_21 = buffer.data(ph + 21);
    const auto *ph_22 = buffer.data(ph + 22);
    const auto *ph_23 = buffer.data(ph + 23);
    const auto *ph_24 = buffer.data(ph + 24);
    const auto *ph_25 = buffer.data(ph + 25);
    const auto *ph_26 = buffer.data(ph + 26);
    const auto *ph_27 = buffer.data(ph + 27);
    const auto *ph_28 = buffer.data(ph + 28);
    const auto *ph_29 = buffer.data(ph + 29);
    const auto *ph_30 = buffer.data(ph + 30);
    const auto *ph_31 = buffer.data(ph + 31);
    const auto *ph_32 = buffer.data(ph + 32);
    const auto *ph_33 = buffer.data(ph + 33);
    const auto *ph_34 = buffer.data(ph + 34);
    const auto *ph_35 = buffer.data(ph + 35);
    const auto *ph_36 = buffer.data(ph + 36);
    const auto *ph_37 = buffer.data(ph + 37);
    const auto *ph_38 = buffer.data(ph + 38);
    const auto *ph_39 = buffer.data(ph + 39);
    const auto *ph_40 = buffer.data(ph + 40);
    const auto *ph_41 = buffer.data(ph + 41);
    const auto *ph_42 = buffer.data(ph + 42);
    const auto *ph_43 = buffer.data(ph + 43);
    const auto *ph_44 = buffer.data(ph + 44);
    const auto *ph_45 = buffer.data(ph + 45);
    const auto *ph_46 = buffer.data(ph + 46);
    const auto *ph_47 = buffer.data(ph + 47);
    const auto *ph_48 = buffer.data(ph + 48);
    const auto *ph_49 = buffer.data(ph + 49);
    const auto *ph_50 = buffer.data(ph + 50);
    const auto *ph_51 = buffer.data(ph + 51);
    const auto *ph_52 = buffer.data(ph + 52);
    const auto *ph_53 = buffer.data(ph + 53);
    const auto *ph_54 = buffer.data(ph + 54);
    const auto *ph_55 = buffer.data(ph + 55);
    const auto *ph_56 = buffer.data(ph + 56);
    const auto *ph_57 = buffer.data(ph + 57);
    const auto *ph_58 = buffer.data(ph + 58);
    const auto *ph_59 = buffer.data(ph + 59);
    const auto *ph_60 = buffer.data(ph + 60);
    const auto *ph_61 = buffer.data(ph + 61);
    const auto *ph_62 = buffer.data(ph + 62);

#pragma omp simd aligned(ph_22, ph_25, ph_27, ph_29, ph_32, ph_34, ph_36, ph_38, \
                         ph_40 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * ph_22[k]
                 - f_1 * ph_27[k]
                 + f_2 * ph_36[k];

        g_1[k] = f_3 * ph_25[k]
                 - f_3 * ph_32[k];

        g_2[k] = -f_4 * ph_22[k]
                 - f_5 * ph_27[k]
                 + f_6 * ph_29[k]
                 + f_7 * ph_36[k]
                 - f_8 * ph_38[k];

        g_3[k] = -f_9 * ph_25[k]
                 - f_9 * ph_32[k]
                 + f_10 * ph_34[k];

        g_4[k] = f_11 * ph_22[k]
                 + f_12 * ph_27[k]
                 - f_13 * ph_29[k]
                 + f_11 * ph_36[k]
                 - f_13 * ph_38[k]
                 + f_14 * ph_40[k];
    }

#pragma omp simd aligned(ph_21, ph_23, ph_24, ph_26, ph_28, ph_30, ph_31, ph_33, ph_35, ph_37, \
                         ph_39, ph_41 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = 1.875 * ph_23[k]
                 + 3.75 * ph_28[k]
                 - 5.0 * ph_30[k]
                 + 1.875 * ph_37[k]
                 - 5.0 * ph_39[k]
                 + ph_41[k];

        g_6[k] = f_11 * ph_21[k]
                 + f_12 * ph_24[k]
                 - f_13 * ph_26[k]
                 + f_11 * ph_31[k]
                 - f_13 * ph_33[k]
                 + f_14 * ph_35[k];

        g_7[k] = -f_15 * ph_23[k]
                 + f_9 * ph_30[k]
                 + f_15 * ph_37[k]
                 - f_9 * ph_39[k];

        g_8[k] = -f_7 * ph_21[k]
                 + f_5 * ph_24[k]
                 + f_8 * ph_26[k]
                 + f_4 * ph_31[k]
                 - f_6 * ph_33[k];
    }

#pragma omp simd aligned(ph_21, ph_23, ph_24, ph_28, ph_31, ph_37, ph_43, ph_46, ph_48, ph_53, \
                         ph_57 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_9[k] = f_16 * ph_23[k]
                 - f_17 * ph_28[k]
                 + f_16 * ph_37[k];

        g_10[k] = f_2 * ph_21[k]
                  - f_1 * ph_24[k]
                  + f_0 * ph_31[k];

        g_11[k] = f_0 * ph_43[k]
                  - f_1 * ph_48[k]
                  + f_2 * ph_57[k];

        g_12[k] = f_3 * ph_46[k]
                  - f_3 * ph_53[k];
    }

#pragma omp simd aligned(ph_43, ph_46, ph_48, ph_50, ph_53, ph_55, ph_57, ph_59, \
                         ph_61 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = -f_4 * ph_43[k]
                  - f_5 * ph_48[k]
                  + f_6 * ph_50[k]
                  + f_7 * ph_57[k]
                  - f_8 * ph_59[k];

        g_14[k] = -f_9 * ph_46[k]
                  - f_9 * ph_53[k]
                  + f_10 * ph_55[k];

        g_15[k] = f_11 * ph_43[k]
                  + f_12 * ph_48[k]
                  - f_13 * ph_50[k]
                  + f_11 * ph_57[k]
                  - f_13 * ph_59[k]
                  + f_14 * ph_61[k];
    }

#pragma omp simd aligned(ph_42, ph_44, ph_45, ph_47, ph_49, ph_51, ph_52, ph_54, ph_56, ph_58, \
                         ph_60, ph_62 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_16[k] = 1.875 * ph_44[k]
                  + 3.75 * ph_49[k]
                  - 5.0 * ph_51[k]
                  + 1.875 * ph_58[k]
                  - 5.0 * ph_60[k]
                  + ph_62[k];

        g_17[k] = f_11 * ph_42[k]
                  + f_12 * ph_45[k]
                  - f_13 * ph_47[k]
                  + f_11 * ph_52[k]
                  - f_13 * ph_54[k]
                  + f_14 * ph_56[k];

        g_18[k] = -f_15 * ph_44[k]
                  + f_9 * ph_51[k]
                  + f_15 * ph_58[k]
                  - f_9 * ph_60[k];

        g_19[k] = -f_7 * ph_42[k]
                  + f_5 * ph_45[k]
                  + f_8 * ph_47[k]
                  + f_4 * ph_52[k]
                  - f_6 * ph_54[k];
    }

#pragma omp simd aligned(ph_1, ph_4, ph_6, ph_11, ph_15, ph_42, ph_44, ph_45, ph_49, ph_52, \
                         ph_58 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = f_16 * ph_44[k]
                  - f_17 * ph_49[k]
                  + f_16 * ph_58[k];

        g_21[k] = f_2 * ph_42[k]
                  - f_1 * ph_45[k]
                  + f_0 * ph_52[k];

        g_22[k] = f_0 * ph_1[k]
                  - f_1 * ph_6[k]
                  + f_2 * ph_15[k];

        g_23[k] = f_3 * ph_4[k]
                  - f_3 * ph_11[k];
    }

#pragma omp simd aligned(ph_1, ph_4, ph_6, ph_8, ph_11, ph_13, ph_15, ph_17, \
                         ph_19 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = -f_4 * ph_1[k]
                  - f_5 * ph_6[k]
                  + f_6 * ph_8[k]
                  + f_7 * ph_15[k]
                  - f_8 * ph_17[k];

        g_25[k] = -f_9 * ph_4[k]
                  - f_9 * ph_11[k]
                  + f_10 * ph_13[k];

        g_26[k] = f_11 * ph_1[k]
                  + f_12 * ph_6[k]
                  - f_13 * ph_8[k]
                  + f_11 * ph_15[k]
                  - f_13 * ph_17[k]
                  + f_14 * ph_19[k];
    }

#pragma omp simd aligned(ph_0, ph_2, ph_3, ph_5, ph_7, ph_9, ph_10, ph_12, ph_14, ph_16, \
                         ph_18, ph_20 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = 1.875 * ph_2[k]
                  + 3.75 * ph_7[k]
                  - 5.0 * ph_9[k]
                  + 1.875 * ph_16[k]
                  - 5.0 * ph_18[k]
                  + ph_20[k];

        g_28[k] = f_11 * ph_0[k]
                  + f_12 * ph_3[k]
                  - f_13 * ph_5[k]
                  + f_11 * ph_10[k]
                  - f_13 * ph_12[k]
                  + f_14 * ph_14[k];

        g_29[k] = -f_15 * ph_2[k]
                  + f_9 * ph_9[k]
                  + f_15 * ph_16[k]
                  - f_9 * ph_18[k];

        g_30[k] = -f_7 * ph_0[k]
                  + f_5 * ph_3[k]
                  + f_8 * ph_5[k]
                  + f_4 * ph_10[k]
                  - f_6 * ph_12[k];
    }

#pragma omp simd aligned(ph_0, ph_2, ph_3, ph_7, ph_10, ph_16 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = f_16 * ph_2[k]
                  - f_17 * ph_7[k]
                  + f_16 * ph_16[k];

        g_32[k] = f_2 * ph_0[k]
                  - f_1 * ph_3[k]
                  + f_0 * ph_10[k];
    }
}

}  // namespace simdtrf
