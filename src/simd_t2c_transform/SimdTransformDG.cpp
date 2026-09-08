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


#include "SimdTransformDG.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_dg(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t dg,
             const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 0.5 * std::sqrt(105.0);
    const auto f_1 = 0.75 * std::sqrt(210.0);
    const auto f_2 = 0.25 * std::sqrt(210.0);
    const auto f_3 = 0.5 * std::sqrt(15.0);
    const auto f_4 = 3.0 * std::sqrt(15.0);
    const auto f_5 = 0.75 * std::sqrt(30.0);
    const auto f_6 = std::sqrt(30.0);
    const auto f_7 = 0.375 * std::sqrt(3.0);
    const auto f_8 = 0.75 * std::sqrt(3.0);
    const auto f_9 = 3.0 * std::sqrt(3.0);
    const auto f_10 = std::sqrt(3.0);
    const auto f_11 = 0.25 * std::sqrt(15.0);
    const auto f_12 = 1.5 * std::sqrt(15.0);
    const auto f_13 = 0.125 * std::sqrt(105.0);
    const auto f_14 = 0.75 * std::sqrt(105.0);
    const auto f_15 = 0.25 * std::sqrt(35.0);
    const auto f_16 = 0.5 * std::sqrt(35.0);
    const auto f_17 = 0.375 * std::sqrt(70.0);
    const auto f_18 = 0.125 * std::sqrt(70.0);
    const auto f_19 = 0.75 * std::sqrt(70.0);
    const auto f_20 = 0.25 * std::sqrt(70.0);
    const auto f_21 = 0.25 * std::sqrt(5.0);
    const auto f_22 = 1.5 * std::sqrt(5.0);
    const auto f_23 = 0.5 * std::sqrt(5.0);
    const auto f_24 = 3.0 * std::sqrt(5.0);
    const auto f_25 = 0.375 * std::sqrt(10.0);
    const auto f_26 = 0.5 * std::sqrt(10.0);
    const auto f_27 = 0.75 * std::sqrt(10.0);
    const auto f_28 = std::sqrt(10.0);
    const auto f_29 = 0.125 * std::sqrt(5.0);
    const auto f_30 = 0.75 * std::sqrt(5.0);
    const auto f_31 = 0.0625 * std::sqrt(35.0);
    const auto f_32 = 0.375 * std::sqrt(35.0);
    const auto f_33 = 0.125 * std::sqrt(35.0);
    const auto f_34 = 0.75 * std::sqrt(35.0);
    const auto f_35 = 0.25 * std::sqrt(105.0);
    const auto f_36 = 0.375 * std::sqrt(210.0);
    const auto f_37 = 0.125 * std::sqrt(210.0);
    const auto f_38 = 0.375 * std::sqrt(30.0);
    const auto f_39 = 0.5 * std::sqrt(30.0);
    const auto f_40 = 0.1875 * std::sqrt(3.0);
    const auto f_41 = 1.5 * std::sqrt(3.0);
    const auto f_42 = 0.5 * std::sqrt(3.0);
    const auto f_43 = 0.125 * std::sqrt(15.0);
    const auto f_44 = 0.75 * std::sqrt(15.0);
    const auto f_45 = 0.0625 * std::sqrt(105.0);
    const auto f_46 = 0.375 * std::sqrt(105.0);

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
    auto *g_39 = values + 39 * nvalues;
    auto *g_40 = values + 40 * nvalues;
    auto *g_41 = values + 41 * nvalues;
    auto *g_42 = values + 42 * nvalues;
    auto *g_43 = values + 43 * nvalues;
    auto *g_44 = values + 44 * nvalues;

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_11 = buffer.data(dg + 11);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_15 = buffer.data(dg + 15);
    const auto *dg_16 = buffer.data(dg + 16);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_19 = buffer.data(dg + 19);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_21 = buffer.data(dg + 21);
    const auto *dg_22 = buffer.data(dg + 22);
    const auto *dg_23 = buffer.data(dg + 23);
    const auto *dg_24 = buffer.data(dg + 24);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_26 = buffer.data(dg + 26);
    const auto *dg_27 = buffer.data(dg + 27);
    const auto *dg_28 = buffer.data(dg + 28);
    const auto *dg_29 = buffer.data(dg + 29);
    const auto *dg_30 = buffer.data(dg + 30);
    const auto *dg_31 = buffer.data(dg + 31);
    const auto *dg_32 = buffer.data(dg + 32);
    const auto *dg_33 = buffer.data(dg + 33);
    const auto *dg_34 = buffer.data(dg + 34);
    const auto *dg_35 = buffer.data(dg + 35);
    const auto *dg_36 = buffer.data(dg + 36);
    const auto *dg_37 = buffer.data(dg + 37);
    const auto *dg_38 = buffer.data(dg + 38);
    const auto *dg_39 = buffer.data(dg + 39);
    const auto *dg_40 = buffer.data(dg + 40);
    const auto *dg_41 = buffer.data(dg + 41);
    const auto *dg_42 = buffer.data(dg + 42);
    const auto *dg_43 = buffer.data(dg + 43);
    const auto *dg_44 = buffer.data(dg + 44);
    const auto *dg_45 = buffer.data(dg + 45);
    const auto *dg_46 = buffer.data(dg + 46);
    const auto *dg_47 = buffer.data(dg + 47);
    const auto *dg_48 = buffer.data(dg + 48);
    const auto *dg_49 = buffer.data(dg + 49);
    const auto *dg_50 = buffer.data(dg + 50);
    const auto *dg_51 = buffer.data(dg + 51);
    const auto *dg_52 = buffer.data(dg + 52);
    const auto *dg_53 = buffer.data(dg + 53);
    const auto *dg_54 = buffer.data(dg + 54);
    const auto *dg_55 = buffer.data(dg + 55);
    const auto *dg_56 = buffer.data(dg + 56);
    const auto *dg_57 = buffer.data(dg + 57);
    const auto *dg_58 = buffer.data(dg + 58);
    const auto *dg_59 = buffer.data(dg + 59);
    const auto *dg_60 = buffer.data(dg + 60);
    const auto *dg_61 = buffer.data(dg + 61);
    const auto *dg_62 = buffer.data(dg + 62);
    const auto *dg_63 = buffer.data(dg + 63);
    const auto *dg_64 = buffer.data(dg + 64);
    const auto *dg_65 = buffer.data(dg + 65);
    const auto *dg_66 = buffer.data(dg + 66);
    const auto *dg_67 = buffer.data(dg + 67);
    const auto *dg_68 = buffer.data(dg + 68);
    const auto *dg_69 = buffer.data(dg + 69);
    const auto *dg_70 = buffer.data(dg + 70);
    const auto *dg_71 = buffer.data(dg + 71);
    const auto *dg_72 = buffer.data(dg + 72);
    const auto *dg_73 = buffer.data(dg + 73);
    const auto *dg_74 = buffer.data(dg + 74);
    const auto *dg_75 = buffer.data(dg + 75);
    const auto *dg_76 = buffer.data(dg + 76);
    const auto *dg_77 = buffer.data(dg + 77);
    const auto *dg_78 = buffer.data(dg + 78);
    const auto *dg_79 = buffer.data(dg + 79);
    const auto *dg_80 = buffer.data(dg + 80);
    const auto *dg_81 = buffer.data(dg + 81);
    const auto *dg_82 = buffer.data(dg + 82);
    const auto *dg_83 = buffer.data(dg + 83);
    const auto *dg_84 = buffer.data(dg + 84);
    const auto *dg_85 = buffer.data(dg + 85);
    const auto *dg_86 = buffer.data(dg + 86);
    const auto *dg_87 = buffer.data(dg + 87);
    const auto *dg_88 = buffer.data(dg + 88);
    const auto *dg_89 = buffer.data(dg + 89);

#pragma omp simd aligned(dg_16, dg_19, dg_21, dg_23, dg_26, dg_28 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * dg_16[k]
                 - f_0 * dg_21[k];

        g_1[k] = f_1 * dg_19[k]
                 - f_2 * dg_26[k];

        g_2[k] = -f_3 * dg_16[k]
                 - f_3 * dg_21[k]
                 + f_4 * dg_23[k];

        g_3[k] = -f_5 * dg_19[k]
                 - f_5 * dg_26[k]
                 + f_6 * dg_28[k];
    }

#pragma omp simd aligned(dg_15, dg_17, dg_18, dg_20, dg_22, dg_24, dg_25, dg_27, \
                         dg_29 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = f_7 * dg_15[k]
                 + f_8 * dg_18[k]
                 - f_9 * dg_20[k]
                 + f_7 * dg_25[k]
                 - f_9 * dg_27[k]
                 + f_10 * dg_29[k];

        g_5[k] = -f_5 * dg_17[k]
                 - f_5 * dg_22[k]
                 + f_6 * dg_24[k];

        g_6[k] = -f_11 * dg_15[k]
                 + f_12 * dg_20[k]
                 + f_11 * dg_25[k]
                 - f_12 * dg_27[k];

        g_7[k] = f_2 * dg_17[k]
                 - f_1 * dg_22[k];

        g_8[k] = f_13 * dg_15[k]
                 - f_14 * dg_18[k]
                 + f_13 * dg_25[k];
    }

#pragma omp simd aligned(dg_61, dg_64, dg_66, dg_68, dg_71, dg_73 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_9[k] = f_0 * dg_61[k]
                 - f_0 * dg_66[k];

        g_10[k] = f_1 * dg_64[k]
                  - f_2 * dg_71[k];

        g_11[k] = -f_3 * dg_61[k]
                  - f_3 * dg_66[k]
                  + f_4 * dg_68[k];

        g_12[k] = -f_5 * dg_64[k]
                  - f_5 * dg_71[k]
                  + f_6 * dg_73[k];
    }

#pragma omp simd aligned(dg_60, dg_62, dg_63, dg_65, dg_67, dg_69, dg_70, dg_72, \
                         dg_74 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = f_7 * dg_60[k]
                  + f_8 * dg_63[k]
                  - f_9 * dg_65[k]
                  + f_7 * dg_70[k]
                  - f_9 * dg_72[k]
                  + f_10 * dg_74[k];

        g_14[k] = -f_5 * dg_62[k]
                  - f_5 * dg_67[k]
                  + f_6 * dg_69[k];

        g_15[k] = -f_11 * dg_60[k]
                  + f_12 * dg_65[k]
                  + f_11 * dg_70[k]
                  - f_12 * dg_72[k];

        g_16[k] = f_2 * dg_62[k]
                  - f_1 * dg_67[k];

        g_17[k] = f_13 * dg_60[k]
                  - f_14 * dg_63[k]
                  + f_13 * dg_70[k];
    }

#pragma omp simd aligned(dg_1, dg_4, dg_6, dg_11, dg_46, dg_49, dg_51, dg_56, dg_76, dg_79, \
                         dg_81, dg_86 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = -f_15 * dg_1[k]
                  + f_15 * dg_6[k]
                  - f_15 * dg_46[k]
                  + f_15 * dg_51[k]
                  + f_16 * dg_76[k]
                  - f_16 * dg_81[k];

        g_19[k] = -f_17 * dg_4[k]
                  + f_18 * dg_11[k]
                  - f_17 * dg_49[k]
                  + f_18 * dg_56[k]
                  + f_19 * dg_79[k]
                  - f_20 * dg_86[k];
    }

#pragma omp simd aligned(dg_1, dg_6, dg_8, dg_46, dg_51, dg_53, dg_76, dg_81, \
                         dg_83 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = f_21 * dg_1[k]
                  + f_21 * dg_6[k]
                  - f_22 * dg_8[k]
                  + f_21 * dg_46[k]
                  + f_21 * dg_51[k]
                  - f_22 * dg_53[k]
                  - f_23 * dg_76[k]
                  - f_23 * dg_81[k]
                  + f_24 * dg_83[k];
    }

#pragma omp simd aligned(dg_4, dg_11, dg_13, dg_49, dg_56, dg_58, dg_79, dg_86, \
                         dg_88 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = f_25 * dg_4[k]
                  + f_25 * dg_11[k]
                  - f_26 * dg_13[k]
                  + f_25 * dg_49[k]
                  + f_25 * dg_56[k]
                  - f_26 * dg_58[k]
                  - f_27 * dg_79[k]
                  - f_27 * dg_86[k]
                  + f_28 * dg_88[k];
    }

#pragma omp simd aligned(dg_0, dg_3, dg_5, dg_10, dg_12, dg_14, dg_45, dg_48, dg_50, dg_55, \
                         dg_57, dg_59, dg_75, dg_78, dg_80, dg_85, dg_87, \
                         dg_89 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = -0.1875 * dg_0[k]
                  - 0.375 * dg_3[k]
                  + 1.5 * dg_5[k]
                  - 0.1875 * dg_10[k]
                  + 1.5 * dg_12[k]
                  - 0.5 * dg_14[k]
                  - 0.1875 * dg_45[k]
                  - 0.375 * dg_48[k]
                  + 1.5 * dg_50[k]
                  - 0.1875 * dg_55[k]
                  + 1.5 * dg_57[k]
                  - 0.5 * dg_59[k]
                  + 0.375 * dg_75[k]
                  + 0.75 * dg_78[k]
                  - 3.0 * dg_80[k]
                  + 0.375 * dg_85[k]
                  - 3.0 * dg_87[k]
                  + dg_89[k];
    }

#pragma omp simd aligned(dg_2, dg_7, dg_9, dg_47, dg_52, dg_54, dg_77, dg_82, \
                         dg_84 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = f_25 * dg_2[k]
                  + f_25 * dg_7[k]
                  - f_26 * dg_9[k]
                  + f_25 * dg_47[k]
                  + f_25 * dg_52[k]
                  - f_26 * dg_54[k]
                  - f_27 * dg_77[k]
                  - f_27 * dg_82[k]
                  + f_28 * dg_84[k];
    }

#pragma omp simd aligned(dg_0, dg_5, dg_10, dg_12, dg_45, dg_50, dg_55, dg_57, dg_75, dg_80, \
                         dg_85, dg_87 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = f_29 * dg_0[k]
                  - f_30 * dg_5[k]
                  - f_29 * dg_10[k]
                  + f_30 * dg_12[k]
                  + f_29 * dg_45[k]
                  - f_30 * dg_50[k]
                  - f_29 * dg_55[k]
                  + f_30 * dg_57[k]
                  - f_21 * dg_75[k]
                  + f_22 * dg_80[k]
                  + f_21 * dg_85[k]
                  - f_22 * dg_87[k];
    }

#pragma omp simd aligned(dg_2, dg_7, dg_47, dg_52, dg_77, dg_82 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = -f_18 * dg_2[k]
                  + f_17 * dg_7[k]
                  - f_18 * dg_47[k]
                  + f_17 * dg_52[k]
                  + f_20 * dg_77[k]
                  - f_19 * dg_82[k];
    }

#pragma omp simd aligned(dg_0, dg_3, dg_10, dg_31, dg_34, dg_36, dg_41, dg_45, dg_48, dg_55, \
                         dg_75, dg_78, dg_85 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = -f_31 * dg_0[k]
                  + f_32 * dg_3[k]
                  - f_31 * dg_10[k]
                  - f_31 * dg_45[k]
                  + f_32 * dg_48[k]
                  - f_31 * dg_55[k]
                  + f_33 * dg_75[k]
                  - f_34 * dg_78[k]
                  + f_33 * dg_85[k];

        g_27[k] = f_0 * dg_31[k]
                  - f_0 * dg_36[k];

        g_28[k] = f_1 * dg_34[k]
                  - f_2 * dg_41[k];
    }

#pragma omp simd aligned(dg_30, dg_31, dg_33, dg_34, dg_35, dg_36, dg_38, dg_40, dg_41, dg_42, \
                         dg_43, dg_44 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = -f_3 * dg_31[k]
                  - f_3 * dg_36[k]
                  + f_4 * dg_38[k];

        g_30[k] = -f_5 * dg_34[k]
                  - f_5 * dg_41[k]
                  + f_6 * dg_43[k];

        g_31[k] = f_7 * dg_30[k]
                  + f_8 * dg_33[k]
                  - f_9 * dg_35[k]
                  + f_7 * dg_40[k]
                  - f_9 * dg_42[k]
                  + f_10 * dg_44[k];
    }

#pragma omp simd aligned(dg_30, dg_32, dg_33, dg_35, dg_37, dg_39, dg_40, \
                         dg_42 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = -f_5 * dg_32[k]
                  - f_5 * dg_37[k]
                  + f_6 * dg_39[k];

        g_33[k] = -f_11 * dg_30[k]
                  + f_12 * dg_35[k]
                  + f_11 * dg_40[k]
                  - f_12 * dg_42[k];

        g_34[k] = f_2 * dg_32[k]
                  - f_1 * dg_37[k];

        g_35[k] = f_13 * dg_30[k]
                  - f_14 * dg_33[k]
                  + f_13 * dg_40[k];
    }

#pragma omp simd aligned(dg_1, dg_4, dg_6, dg_8, dg_11, dg_13, dg_46, dg_49, dg_51, dg_53, \
                         dg_56, dg_58 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_36[k] = f_35 * dg_1[k]
                  - f_35 * dg_6[k]
                  - f_35 * dg_46[k]
                  + f_35 * dg_51[k];

        g_37[k] = f_36 * dg_4[k]
                  - f_37 * dg_11[k]
                  - f_36 * dg_49[k]
                  + f_37 * dg_56[k];

        g_38[k] = -f_11 * dg_1[k]
                  - f_11 * dg_6[k]
                  + f_12 * dg_8[k]
                  + f_11 * dg_46[k]
                  + f_11 * dg_51[k]
                  - f_12 * dg_53[k];

        g_39[k] = -f_38 * dg_4[k]
                  - f_38 * dg_11[k]
                  + f_39 * dg_13[k]
                  + f_38 * dg_49[k]
                  + f_38 * dg_56[k]
                  - f_39 * dg_58[k];
    }

#pragma omp simd aligned(dg_0, dg_3, dg_5, dg_10, dg_12, dg_14, dg_45, dg_48, dg_50, dg_55, \
                         dg_57, dg_59 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_40[k] = f_40 * dg_0[k]
                  + f_7 * dg_3[k]
                  - f_41 * dg_5[k]
                  + f_40 * dg_10[k]
                  - f_41 * dg_12[k]
                  + f_42 * dg_14[k]
                  - f_40 * dg_45[k]
                  - f_7 * dg_48[k]
                  + f_41 * dg_50[k]
                  - f_40 * dg_55[k]
                  + f_41 * dg_57[k]
                  - f_42 * dg_59[k];
    }

#pragma omp simd aligned(dg_0, dg_2, dg_5, dg_7, dg_9, dg_10, dg_12, dg_45, dg_47, dg_50, \
                         dg_52, dg_54, dg_55, dg_57 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_41[k] = -f_38 * dg_2[k]
                  - f_38 * dg_7[k]
                  + f_39 * dg_9[k]
                  + f_38 * dg_47[k]
                  + f_38 * dg_52[k]
                  - f_39 * dg_54[k];

        g_42[k] = -f_43 * dg_0[k]
                  + f_44 * dg_5[k]
                  + f_43 * dg_10[k]
                  - f_44 * dg_12[k]
                  + f_43 * dg_45[k]
                  - f_44 * dg_50[k]
                  - f_43 * dg_55[k]
                  + f_44 * dg_57[k];
    }

#pragma omp simd aligned(dg_0, dg_2, dg_3, dg_7, dg_10, dg_45, dg_47, dg_48, dg_52, \
                         dg_55 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_43[k] = f_37 * dg_2[k]
                  - f_36 * dg_7[k]
                  - f_37 * dg_47[k]
                  + f_36 * dg_52[k];

        g_44[k] = f_45 * dg_0[k]
                  - f_46 * dg_3[k]
                  + f_45 * dg_10[k]
                  - f_45 * dg_45[k]
                  + f_46 * dg_48[k]
                  - f_45 * dg_55[k];
    }
}

}  // namespace simdtrf
