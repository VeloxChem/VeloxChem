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


#include "SimdTransformPI.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_pi(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t pi,
             const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 0.1875 * std::sqrt(462.0);
    const auto f_1 = 0.625 * std::sqrt(462.0);
    const auto f_2 = 0.9375 * std::sqrt(154.0);
    const auto f_3 = 1.875 * std::sqrt(154.0);
    const auto f_4 = 0.1875 * std::sqrt(154.0);
    const auto f_5 = 0.75 * std::sqrt(7.0);
    const auto f_6 = 7.5 * std::sqrt(7.0);
    const auto f_7 = 0.5625 * std::sqrt(210.0);
    const auto f_8 = 0.375 * std::sqrt(210.0);
    const auto f_9 = 1.5 * std::sqrt(210.0);
    const auto f_10 = 0.1875 * std::sqrt(210.0);
    const auto f_11 = 0.5 * std::sqrt(210.0);
    const auto f_12 = 0.0625 * std::sqrt(210.0);
    const auto f_13 = 0.125 * std::sqrt(210.0);
    const auto f_14 = std::sqrt(210.0);
    const auto f_15 = 0.625 * std::sqrt(21.0);
    const auto f_16 = 1.25 * std::sqrt(21.0);
    const auto f_17 = 2.5 * std::sqrt(21.0);
    const auto f_18 = std::sqrt(21.0);
    const auto f_19 = 0.03125 * std::sqrt(210.0);
    const auto f_20 = 0.1875 * std::sqrt(7.0);
    const auto f_21 = 0.9375 * std::sqrt(7.0);
    const auto f_22 = 1.875 * std::sqrt(7.0);
    const auto f_23 = 11.25 * std::sqrt(7.0);
    const auto f_24 = 0.03125 * std::sqrt(462.0);
    const auto f_25 = 0.46875 * std::sqrt(462.0);

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
    auto *g_35 = values + 35 * nvalues;
    auto *g_36 = values + 36 * nvalues;
    auto *g_37 = values + 37 * nvalues;
    auto *g_38 = values + 38 * nvalues;

    const auto *pi_0 = buffer.data(pi + 0);
    const auto *pi_1 = buffer.data(pi + 1);
    const auto *pi_2 = buffer.data(pi + 2);
    const auto *pi_3 = buffer.data(pi + 3);
    const auto *pi_4 = buffer.data(pi + 4);
    const auto *pi_5 = buffer.data(pi + 5);
    const auto *pi_6 = buffer.data(pi + 6);
    const auto *pi_7 = buffer.data(pi + 7);
    const auto *pi_8 = buffer.data(pi + 8);
    const auto *pi_9 = buffer.data(pi + 9);
    const auto *pi_10 = buffer.data(pi + 10);
    const auto *pi_11 = buffer.data(pi + 11);
    const auto *pi_12 = buffer.data(pi + 12);
    const auto *pi_13 = buffer.data(pi + 13);
    const auto *pi_14 = buffer.data(pi + 14);
    const auto *pi_15 = buffer.data(pi + 15);
    const auto *pi_16 = buffer.data(pi + 16);
    const auto *pi_17 = buffer.data(pi + 17);
    const auto *pi_18 = buffer.data(pi + 18);
    const auto *pi_19 = buffer.data(pi + 19);
    const auto *pi_20 = buffer.data(pi + 20);
    const auto *pi_21 = buffer.data(pi + 21);
    const auto *pi_22 = buffer.data(pi + 22);
    const auto *pi_23 = buffer.data(pi + 23);
    const auto *pi_24 = buffer.data(pi + 24);
    const auto *pi_25 = buffer.data(pi + 25);
    const auto *pi_26 = buffer.data(pi + 26);
    const auto *pi_27 = buffer.data(pi + 27);
    const auto *pi_28 = buffer.data(pi + 28);
    const auto *pi_29 = buffer.data(pi + 29);
    const auto *pi_30 = buffer.data(pi + 30);
    const auto *pi_31 = buffer.data(pi + 31);
    const auto *pi_32 = buffer.data(pi + 32);
    const auto *pi_33 = buffer.data(pi + 33);
    const auto *pi_34 = buffer.data(pi + 34);
    const auto *pi_35 = buffer.data(pi + 35);
    const auto *pi_36 = buffer.data(pi + 36);
    const auto *pi_37 = buffer.data(pi + 37);
    const auto *pi_38 = buffer.data(pi + 38);
    const auto *pi_39 = buffer.data(pi + 39);
    const auto *pi_40 = buffer.data(pi + 40);
    const auto *pi_41 = buffer.data(pi + 41);
    const auto *pi_42 = buffer.data(pi + 42);
    const auto *pi_43 = buffer.data(pi + 43);
    const auto *pi_44 = buffer.data(pi + 44);
    const auto *pi_45 = buffer.data(pi + 45);
    const auto *pi_46 = buffer.data(pi + 46);
    const auto *pi_47 = buffer.data(pi + 47);
    const auto *pi_48 = buffer.data(pi + 48);
    const auto *pi_49 = buffer.data(pi + 49);
    const auto *pi_50 = buffer.data(pi + 50);
    const auto *pi_51 = buffer.data(pi + 51);
    const auto *pi_52 = buffer.data(pi + 52);
    const auto *pi_53 = buffer.data(pi + 53);
    const auto *pi_54 = buffer.data(pi + 54);
    const auto *pi_55 = buffer.data(pi + 55);
    const auto *pi_56 = buffer.data(pi + 56);
    const auto *pi_57 = buffer.data(pi + 57);
    const auto *pi_58 = buffer.data(pi + 58);
    const auto *pi_59 = buffer.data(pi + 59);
    const auto *pi_60 = buffer.data(pi + 60);
    const auto *pi_61 = buffer.data(pi + 61);
    const auto *pi_62 = buffer.data(pi + 62);
    const auto *pi_63 = buffer.data(pi + 63);
    const auto *pi_64 = buffer.data(pi + 64);
    const auto *pi_65 = buffer.data(pi + 65);
    const auto *pi_66 = buffer.data(pi + 66);
    const auto *pi_67 = buffer.data(pi + 67);
    const auto *pi_68 = buffer.data(pi + 68);
    const auto *pi_69 = buffer.data(pi + 69);
    const auto *pi_70 = buffer.data(pi + 70);
    const auto *pi_71 = buffer.data(pi + 71);
    const auto *pi_72 = buffer.data(pi + 72);
    const auto *pi_73 = buffer.data(pi + 73);
    const auto *pi_74 = buffer.data(pi + 74);
    const auto *pi_75 = buffer.data(pi + 75);
    const auto *pi_76 = buffer.data(pi + 76);
    const auto *pi_77 = buffer.data(pi + 77);
    const auto *pi_78 = buffer.data(pi + 78);
    const auto *pi_79 = buffer.data(pi + 79);
    const auto *pi_80 = buffer.data(pi + 80);
    const auto *pi_81 = buffer.data(pi + 81);
    const auto *pi_82 = buffer.data(pi + 82);
    const auto *pi_83 = buffer.data(pi + 83);

#pragma omp simd aligned(pi_29, pi_32, pi_34, pi_36, pi_39, pi_41, pi_43, pi_45, pi_47, pi_50, \
                         pi_52 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * pi_29[k]
                 - f_1 * pi_34[k]
                 + f_0 * pi_43[k];

        g_1[k] = f_2 * pi_32[k]
                 - f_3 * pi_39[k]
                 + f_4 * pi_50[k];

        g_2[k] = -f_5 * pi_29[k]
                 + f_6 * pi_36[k]
                 + f_5 * pi_43[k]
                 - f_6 * pi_45[k];

        g_3[k] = -f_7 * pi_32[k]
                 - f_8 * pi_39[k]
                 + f_9 * pi_41[k]
                 + f_10 * pi_50[k]
                 - f_11 * pi_52[k];

        g_4[k] = f_12 * pi_29[k]
                 + f_13 * pi_34[k]
                 - f_14 * pi_36[k]
                 + f_12 * pi_43[k]
                 - f_14 * pi_45[k]
                 + f_14 * pi_47[k];
    }

#pragma omp simd aligned(pi_32, pi_39, pi_41, pi_50, pi_52, pi_54 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = f_15 * pi_32[k]
                 + f_16 * pi_39[k]
                 - f_17 * pi_41[k]
                 + f_15 * pi_50[k]
                 - f_17 * pi_52[k]
                 + f_18 * pi_54[k];
    }

#pragma omp simd aligned(pi_28, pi_31, pi_33, pi_38, pi_40, pi_42, pi_49, pi_51, pi_53, \
                         pi_55 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = -0.3125 * pi_28[k]
                 - 0.9375 * pi_31[k]
                 + 5.625 * pi_33[k]
                 - 0.9375 * pi_38[k]
                 + 11.25 * pi_40[k]
                 - 7.5 * pi_42[k]
                 - 0.3125 * pi_49[k]
                 + 5.625 * pi_51[k]
                 - 7.5 * pi_53[k]
                 + pi_55[k];
    }

#pragma omp simd aligned(pi_28, pi_30, pi_31, pi_33, pi_35, pi_37, pi_38, pi_42, pi_44, pi_46, \
                         pi_48, pi_49, pi_51, pi_53 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = f_15 * pi_30[k]
                 + f_16 * pi_35[k]
                 - f_17 * pi_37[k]
                 + f_15 * pi_44[k]
                 - f_17 * pi_46[k]
                 + f_18 * pi_48[k];

        g_8[k] = f_19 * pi_28[k]
                 + f_19 * pi_31[k]
                 - f_11 * pi_33[k]
                 - f_19 * pi_38[k]
                 + f_11 * pi_42[k]
                 - f_19 * pi_49[k]
                 + f_11 * pi_51[k]
                 - f_11 * pi_53[k];
    }

#pragma omp simd aligned(pi_28, pi_30, pi_31, pi_33, pi_35, pi_37, pi_38, pi_40, pi_44, pi_46, \
                         pi_49, pi_51 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_9[k] = -f_10 * pi_30[k]
                 + f_8 * pi_35[k]
                 + f_11 * pi_37[k]
                 + f_7 * pi_44[k]
                 - f_9 * pi_46[k];

        g_10[k] = -f_20 * pi_28[k]
                  + f_21 * pi_31[k]
                  + f_22 * pi_33[k]
                  + f_21 * pi_38[k]
                  - f_23 * pi_40[k]
                  - f_20 * pi_49[k]
                  + f_22 * pi_51[k];

        g_11[k] = f_4 * pi_30[k]
                  - f_3 * pi_35[k]
                  + f_2 * pi_44[k];

        g_12[k] = f_24 * pi_28[k]
                  - f_25 * pi_31[k]
                  + f_25 * pi_38[k]
                  - f_24 * pi_49[k];
    }

#pragma omp simd aligned(pi_57, pi_60, pi_62, pi_64, pi_67, pi_69, pi_71, pi_73, pi_75, pi_78, \
                         pi_80 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = f_0 * pi_57[k]
                  - f_1 * pi_62[k]
                  + f_0 * pi_71[k];

        g_14[k] = f_2 * pi_60[k]
                  - f_3 * pi_67[k]
                  + f_4 * pi_78[k];

        g_15[k] = -f_5 * pi_57[k]
                  + f_6 * pi_64[k]
                  + f_5 * pi_71[k]
                  - f_6 * pi_73[k];

        g_16[k] = -f_7 * pi_60[k]
                  - f_8 * pi_67[k]
                  + f_9 * pi_69[k]
                  + f_10 * pi_78[k]
                  - f_11 * pi_80[k];

        g_17[k] = f_12 * pi_57[k]
                  + f_13 * pi_62[k]
                  - f_14 * pi_64[k]
                  + f_12 * pi_71[k]
                  - f_14 * pi_73[k]
                  + f_14 * pi_75[k];
    }

#pragma omp simd aligned(pi_60, pi_67, pi_69, pi_78, pi_80, pi_82 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = f_15 * pi_60[k]
                  + f_16 * pi_67[k]
                  - f_17 * pi_69[k]
                  + f_15 * pi_78[k]
                  - f_17 * pi_80[k]
                  + f_18 * pi_82[k];
    }

#pragma omp simd aligned(pi_56, pi_59, pi_61, pi_66, pi_68, pi_70, pi_77, pi_79, pi_81, \
                         pi_83 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = -0.3125 * pi_56[k]
                  - 0.9375 * pi_59[k]
                  + 5.625 * pi_61[k]
                  - 0.9375 * pi_66[k]
                  + 11.25 * pi_68[k]
                  - 7.5 * pi_70[k]
                  - 0.3125 * pi_77[k]
                  + 5.625 * pi_79[k]
                  - 7.5 * pi_81[k]
                  + pi_83[k];
    }

#pragma omp simd aligned(pi_56, pi_58, pi_59, pi_61, pi_63, pi_65, pi_66, pi_70, pi_72, pi_74, \
                         pi_76, pi_77, pi_79, pi_81 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = f_15 * pi_58[k]
                  + f_16 * pi_63[k]
                  - f_17 * pi_65[k]
                  + f_15 * pi_72[k]
                  - f_17 * pi_74[k]
                  + f_18 * pi_76[k];

        g_21[k] = f_19 * pi_56[k]
                  + f_19 * pi_59[k]
                  - f_11 * pi_61[k]
                  - f_19 * pi_66[k]
                  + f_11 * pi_70[k]
                  - f_19 * pi_77[k]
                  + f_11 * pi_79[k]
                  - f_11 * pi_81[k];
    }

#pragma omp simd aligned(pi_56, pi_58, pi_59, pi_61, pi_63, pi_65, pi_66, pi_68, pi_72, pi_74, \
                         pi_77, pi_79 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = -f_10 * pi_58[k]
                  + f_8 * pi_63[k]
                  + f_11 * pi_65[k]
                  + f_7 * pi_72[k]
                  - f_9 * pi_74[k];

        g_23[k] = -f_20 * pi_56[k]
                  + f_21 * pi_59[k]
                  + f_22 * pi_61[k]
                  + f_21 * pi_66[k]
                  - f_23 * pi_68[k]
                  - f_20 * pi_77[k]
                  + f_22 * pi_79[k];

        g_24[k] = f_4 * pi_58[k]
                  - f_3 * pi_63[k]
                  + f_2 * pi_72[k];

        g_25[k] = f_24 * pi_56[k]
                  - f_25 * pi_59[k]
                  + f_25 * pi_66[k]
                  - f_24 * pi_77[k];
    }

#pragma omp simd aligned(pi_1, pi_4, pi_6, pi_8, pi_11, pi_13, pi_15, pi_17, pi_19, pi_22, \
                         pi_24 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = f_0 * pi_1[k]
                  - f_1 * pi_6[k]
                  + f_0 * pi_15[k];

        g_27[k] = f_2 * pi_4[k]
                  - f_3 * pi_11[k]
                  + f_4 * pi_22[k];

        g_28[k] = -f_5 * pi_1[k]
                  + f_6 * pi_8[k]
                  + f_5 * pi_15[k]
                  - f_6 * pi_17[k];

        g_29[k] = -f_7 * pi_4[k]
                  - f_8 * pi_11[k]
                  + f_9 * pi_13[k]
                  + f_10 * pi_22[k]
                  - f_11 * pi_24[k];

        g_30[k] = f_12 * pi_1[k]
                  + f_13 * pi_6[k]
                  - f_14 * pi_8[k]
                  + f_12 * pi_15[k]
                  - f_14 * pi_17[k]
                  + f_14 * pi_19[k];
    }

#pragma omp simd aligned(pi_4, pi_11, pi_13, pi_22, pi_24, pi_26 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = f_15 * pi_4[k]
                  + f_16 * pi_11[k]
                  - f_17 * pi_13[k]
                  + f_15 * pi_22[k]
                  - f_17 * pi_24[k]
                  + f_18 * pi_26[k];
    }

#pragma omp simd aligned(pi_0, pi_3, pi_5, pi_10, pi_12, pi_14, pi_21, pi_23, pi_25, \
                         pi_27 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = -0.3125 * pi_0[k]
                  - 0.9375 * pi_3[k]
                  + 5.625 * pi_5[k]
                  - 0.9375 * pi_10[k]
                  + 11.25 * pi_12[k]
                  - 7.5 * pi_14[k]
                  - 0.3125 * pi_21[k]
                  + 5.625 * pi_23[k]
                  - 7.5 * pi_25[k]
                  + pi_27[k];
    }

#pragma omp simd aligned(pi_0, pi_2, pi_3, pi_5, pi_7, pi_9, pi_10, pi_14, pi_16, pi_18, \
                         pi_20, pi_21, pi_23, pi_25 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = f_15 * pi_2[k]
                  + f_16 * pi_7[k]
                  - f_17 * pi_9[k]
                  + f_15 * pi_16[k]
                  - f_17 * pi_18[k]
                  + f_18 * pi_20[k];

        g_34[k] = f_19 * pi_0[k]
                  + f_19 * pi_3[k]
                  - f_11 * pi_5[k]
                  - f_19 * pi_10[k]
                  + f_11 * pi_14[k]
                  - f_19 * pi_21[k]
                  + f_11 * pi_23[k]
                  - f_11 * pi_25[k];
    }

#pragma omp simd aligned(pi_0, pi_2, pi_3, pi_5, pi_7, pi_9, pi_10, pi_12, pi_16, pi_18, \
                         pi_21, pi_23 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = -f_10 * pi_2[k]
                  + f_8 * pi_7[k]
                  + f_11 * pi_9[k]
                  + f_7 * pi_16[k]
                  - f_9 * pi_18[k];

        g_36[k] = -f_20 * pi_0[k]
                  + f_21 * pi_3[k]
                  + f_22 * pi_5[k]
                  + f_21 * pi_10[k]
                  - f_23 * pi_12[k]
                  - f_20 * pi_21[k]
                  + f_22 * pi_23[k];

        g_37[k] = f_4 * pi_2[k]
                  - f_3 * pi_7[k]
                  + f_2 * pi_16[k];

        g_38[k] = f_24 * pi_0[k]
                  - f_25 * pi_3[k]
                  + f_25 * pi_10[k]
                  - f_24 * pi_21[k];
    }
}

}  // namespace simdtrf
