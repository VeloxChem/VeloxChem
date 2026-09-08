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


#include "SimdTransformGD.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_gd(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t gd,
             const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 0.5 * std::sqrt(105.0);
    const auto f_1 = 0.25 * std::sqrt(35.0);
    const auto f_2 = 0.5 * std::sqrt(35.0);
    const auto f_3 = 0.25 * std::sqrt(105.0);
    const auto f_4 = 0.75 * std::sqrt(210.0);
    const auto f_5 = 0.25 * std::sqrt(210.0);
    const auto f_6 = 0.375 * std::sqrt(70.0);
    const auto f_7 = 0.75 * std::sqrt(70.0);
    const auto f_8 = 0.125 * std::sqrt(70.0);
    const auto f_9 = 0.25 * std::sqrt(70.0);
    const auto f_10 = 0.375 * std::sqrt(210.0);
    const auto f_11 = 0.125 * std::sqrt(210.0);
    const auto f_12 = 0.5 * std::sqrt(15.0);
    const auto f_13 = 3.0 * std::sqrt(15.0);
    const auto f_14 = 0.25 * std::sqrt(5.0);
    const auto f_15 = 0.5 * std::sqrt(5.0);
    const auto f_16 = 1.5 * std::sqrt(5.0);
    const auto f_17 = 3.0 * std::sqrt(5.0);
    const auto f_18 = 0.25 * std::sqrt(15.0);
    const auto f_19 = 1.5 * std::sqrt(15.0);
    const auto f_20 = 0.75 * std::sqrt(30.0);
    const auto f_21 = std::sqrt(30.0);
    const auto f_22 = 0.375 * std::sqrt(10.0);
    const auto f_23 = 0.75 * std::sqrt(10.0);
    const auto f_24 = 0.5 * std::sqrt(10.0);
    const auto f_25 = std::sqrt(10.0);
    const auto f_26 = 0.375 * std::sqrt(30.0);
    const auto f_27 = 0.5 * std::sqrt(30.0);
    const auto f_28 = 0.375 * std::sqrt(3.0);
    const auto f_29 = 0.75 * std::sqrt(3.0);
    const auto f_30 = 3.0 * std::sqrt(3.0);
    const auto f_31 = std::sqrt(3.0);
    const auto f_32 = 0.1875 * std::sqrt(3.0);
    const auto f_33 = 1.5 * std::sqrt(3.0);
    const auto f_34 = 0.5 * std::sqrt(3.0);
    const auto f_35 = 0.125 * std::sqrt(5.0);
    const auto f_36 = 0.75 * std::sqrt(5.0);
    const auto f_37 = 0.125 * std::sqrt(15.0);
    const auto f_38 = 0.75 * std::sqrt(15.0);
    const auto f_39 = 0.125 * std::sqrt(105.0);
    const auto f_40 = 0.75 * std::sqrt(105.0);
    const auto f_41 = 0.0625 * std::sqrt(35.0);
    const auto f_42 = 0.125 * std::sqrt(35.0);
    const auto f_43 = 0.375 * std::sqrt(35.0);
    const auto f_44 = 0.75 * std::sqrt(35.0);
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

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_1 = buffer.data(gd + 1);
    const auto *gd_2 = buffer.data(gd + 2);
    const auto *gd_3 = buffer.data(gd + 3);
    const auto *gd_4 = buffer.data(gd + 4);
    const auto *gd_5 = buffer.data(gd + 5);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_7 = buffer.data(gd + 7);
    const auto *gd_8 = buffer.data(gd + 8);
    const auto *gd_9 = buffer.data(gd + 9);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_13 = buffer.data(gd + 13);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_15 = buffer.data(gd + 15);
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_17 = buffer.data(gd + 17);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_20 = buffer.data(gd + 20);
    const auto *gd_21 = buffer.data(gd + 21);
    const auto *gd_22 = buffer.data(gd + 22);
    const auto *gd_23 = buffer.data(gd + 23);
    const auto *gd_24 = buffer.data(gd + 24);
    const auto *gd_25 = buffer.data(gd + 25);
    const auto *gd_26 = buffer.data(gd + 26);
    const auto *gd_27 = buffer.data(gd + 27);
    const auto *gd_28 = buffer.data(gd + 28);
    const auto *gd_29 = buffer.data(gd + 29);
    const auto *gd_30 = buffer.data(gd + 30);
    const auto *gd_31 = buffer.data(gd + 31);
    const auto *gd_32 = buffer.data(gd + 32);
    const auto *gd_33 = buffer.data(gd + 33);
    const auto *gd_34 = buffer.data(gd + 34);
    const auto *gd_35 = buffer.data(gd + 35);
    const auto *gd_36 = buffer.data(gd + 36);
    const auto *gd_37 = buffer.data(gd + 37);
    const auto *gd_38 = buffer.data(gd + 38);
    const auto *gd_39 = buffer.data(gd + 39);
    const auto *gd_40 = buffer.data(gd + 40);
    const auto *gd_41 = buffer.data(gd + 41);
    const auto *gd_42 = buffer.data(gd + 42);
    const auto *gd_43 = buffer.data(gd + 43);
    const auto *gd_44 = buffer.data(gd + 44);
    const auto *gd_45 = buffer.data(gd + 45);
    const auto *gd_46 = buffer.data(gd + 46);
    const auto *gd_47 = buffer.data(gd + 47);
    const auto *gd_48 = buffer.data(gd + 48);
    const auto *gd_49 = buffer.data(gd + 49);
    const auto *gd_50 = buffer.data(gd + 50);
    const auto *gd_51 = buffer.data(gd + 51);
    const auto *gd_52 = buffer.data(gd + 52);
    const auto *gd_53 = buffer.data(gd + 53);
    const auto *gd_54 = buffer.data(gd + 54);
    const auto *gd_55 = buffer.data(gd + 55);
    const auto *gd_56 = buffer.data(gd + 56);
    const auto *gd_57 = buffer.data(gd + 57);
    const auto *gd_58 = buffer.data(gd + 58);
    const auto *gd_59 = buffer.data(gd + 59);
    const auto *gd_60 = buffer.data(gd + 60);
    const auto *gd_61 = buffer.data(gd + 61);
    const auto *gd_62 = buffer.data(gd + 62);
    const auto *gd_63 = buffer.data(gd + 63);
    const auto *gd_64 = buffer.data(gd + 64);
    const auto *gd_65 = buffer.data(gd + 65);
    const auto *gd_66 = buffer.data(gd + 66);
    const auto *gd_67 = buffer.data(gd + 67);
    const auto *gd_68 = buffer.data(gd + 68);
    const auto *gd_69 = buffer.data(gd + 69);
    const auto *gd_70 = buffer.data(gd + 70);
    const auto *gd_71 = buffer.data(gd + 71);
    const auto *gd_72 = buffer.data(gd + 72);
    const auto *gd_73 = buffer.data(gd + 73);
    const auto *gd_74 = buffer.data(gd + 74);
    const auto *gd_75 = buffer.data(gd + 75);
    const auto *gd_76 = buffer.data(gd + 76);
    const auto *gd_77 = buffer.data(gd + 77);
    const auto *gd_78 = buffer.data(gd + 78);
    const auto *gd_79 = buffer.data(gd + 79);
    const auto *gd_80 = buffer.data(gd + 80);
    const auto *gd_81 = buffer.data(gd + 81);
    const auto *gd_82 = buffer.data(gd + 82);
    const auto *gd_83 = buffer.data(gd + 83);
    const auto *gd_84 = buffer.data(gd + 84);
    const auto *gd_85 = buffer.data(gd + 85);
    const auto *gd_86 = buffer.data(gd + 86);
    const auto *gd_87 = buffer.data(gd + 87);
    const auto *gd_88 = buffer.data(gd + 88);
    const auto *gd_89 = buffer.data(gd + 89);

#pragma omp simd aligned(gd_6, gd_7, gd_8, gd_9, gd_10, gd_11, gd_36, gd_37, gd_38, gd_39, \
                         gd_40, gd_41 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * gd_7[k]
                 - f_0 * gd_37[k];

        g_1[k] = f_0 * gd_10[k]
                 - f_0 * gd_40[k];

        g_2[k] = -f_1 * gd_6[k]
                 - f_1 * gd_9[k]
                 + f_2 * gd_11[k]
                 + f_1 * gd_36[k]
                 + f_1 * gd_39[k]
                 - f_2 * gd_41[k];

        g_3[k] = f_0 * gd_8[k]
                 - f_0 * gd_38[k];
    }

#pragma omp simd aligned(gd_6, gd_9, gd_25, gd_28, gd_36, gd_39, gd_67, \
                         gd_70 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = f_3 * gd_6[k]
                 - f_3 * gd_9[k]
                 - f_3 * gd_36[k]
                 + f_3 * gd_39[k];

        g_5[k] = f_4 * gd_25[k]
                 - f_5 * gd_67[k];

        g_6[k] = f_4 * gd_28[k]
                 - f_5 * gd_70[k];
    }

#pragma omp simd aligned(gd_7, gd_24, gd_26, gd_27, gd_29, gd_37, gd_49, gd_66, gd_68, gd_69, \
                         gd_71 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = -f_6 * gd_24[k]
                 - f_6 * gd_27[k]
                 + f_7 * gd_29[k]
                 + f_8 * gd_66[k]
                 + f_8 * gd_69[k]
                 - f_9 * gd_71[k];

        g_8[k] = f_4 * gd_26[k]
                 - f_5 * gd_68[k];

        g_9[k] = f_10 * gd_24[k]
                 - f_10 * gd_27[k]
                 - f_11 * gd_66[k]
                 + f_11 * gd_69[k];

        g_10[k] = -f_12 * gd_7[k]
                  - f_12 * gd_37[k]
                  + f_13 * gd_49[k];
    }

#pragma omp simd aligned(gd_6, gd_9, gd_10, gd_11, gd_36, gd_39, gd_40, gd_41, gd_48, gd_51, \
                         gd_52, gd_53 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_11[k] = -f_12 * gd_10[k]
                  - f_12 * gd_40[k]
                  + f_13 * gd_52[k];

        g_12[k] = f_14 * gd_6[k]
                  + f_14 * gd_9[k]
                  - f_15 * gd_11[k]
                  + f_14 * gd_36[k]
                  + f_14 * gd_39[k]
                  - f_15 * gd_41[k]
                  - f_16 * gd_48[k]
                  - f_16 * gd_51[k]
                  + f_17 * gd_53[k];
    }

#pragma omp simd aligned(gd_6, gd_8, gd_9, gd_25, gd_36, gd_38, gd_39, gd_48, gd_50, gd_51, \
                         gd_67, gd_79 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = -f_12 * gd_8[k]
                  - f_12 * gd_38[k]
                  + f_13 * gd_50[k];

        g_14[k] = -f_18 * gd_6[k]
                  + f_18 * gd_9[k]
                  - f_18 * gd_36[k]
                  + f_18 * gd_39[k]
                  + f_19 * gd_48[k]
                  - f_19 * gd_51[k];

        g_15[k] = -f_20 * gd_25[k]
                  - f_20 * gd_67[k]
                  + f_21 * gd_79[k];
    }

#pragma omp simd aligned(gd_24, gd_27, gd_28, gd_29, gd_66, gd_69, gd_70, gd_71, gd_78, gd_81, \
                         gd_82, gd_83 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_16[k] = -f_20 * gd_28[k]
                  - f_20 * gd_70[k]
                  + f_21 * gd_82[k];

        g_17[k] = f_22 * gd_24[k]
                  + f_22 * gd_27[k]
                  - f_23 * gd_29[k]
                  + f_22 * gd_66[k]
                  + f_22 * gd_69[k]
                  - f_23 * gd_71[k]
                  - f_24 * gd_78[k]
                  - f_24 * gd_81[k]
                  + f_25 * gd_83[k];
    }

#pragma omp simd aligned(gd_24, gd_26, gd_27, gd_66, gd_68, gd_69, gd_78, gd_80, \
                         gd_81 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = -f_20 * gd_26[k]
                  - f_20 * gd_68[k]
                  + f_21 * gd_80[k];

        g_19[k] = -f_26 * gd_24[k]
                  + f_26 * gd_27[k]
                  - f_26 * gd_66[k]
                  + f_26 * gd_69[k]
                  + f_27 * gd_78[k]
                  - f_27 * gd_81[k];
    }

#pragma omp simd aligned(gd_1, gd_4, gd_19, gd_22, gd_31, gd_34, gd_61, gd_64, gd_73, gd_76, \
                         gd_85, gd_88 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = f_28 * gd_1[k]
                  + f_29 * gd_19[k]
                  - f_30 * gd_31[k]
                  + f_28 * gd_61[k]
                  - f_30 * gd_73[k]
                  + f_31 * gd_85[k];

        g_21[k] = f_28 * gd_4[k]
                  + f_29 * gd_22[k]
                  - f_30 * gd_34[k]
                  + f_28 * gd_64[k]
                  - f_30 * gd_76[k]
                  + f_31 * gd_88[k];
    }

#pragma omp simd aligned(gd_0, gd_3, gd_5, gd_18, gd_21, gd_23, gd_30, gd_33, gd_35, gd_60, \
                         gd_63, gd_65, gd_72, gd_75, gd_77, gd_84, gd_87, \
                         gd_89 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = -0.1875 * gd_0[k]
                  - 0.1875 * gd_3[k]
                  + 0.375 * gd_5[k]
                  - 0.375 * gd_18[k]
                  - 0.375 * gd_21[k]
                  + 0.75 * gd_23[k]
                  + 1.5 * gd_30[k]
                  + 1.5 * gd_33[k]
                  - 3.0 * gd_35[k]
                  - 0.1875 * gd_60[k]
                  - 0.1875 * gd_63[k]
                  + 0.375 * gd_65[k]
                  + 1.5 * gd_72[k]
                  + 1.5 * gd_75[k]
                  - 3.0 * gd_77[k]
                  - 0.5 * gd_84[k]
                  - 0.5 * gd_87[k]
                  + gd_89[k];
    }

#pragma omp simd aligned(gd_2, gd_20, gd_32, gd_62, gd_74, gd_86 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = f_28 * gd_2[k]
                  + f_29 * gd_20[k]
                  - f_30 * gd_32[k]
                  + f_28 * gd_62[k]
                  - f_30 * gd_74[k]
                  + f_31 * gd_86[k];
    }

#pragma omp simd aligned(gd_0, gd_3, gd_18, gd_21, gd_30, gd_33, gd_60, gd_63, gd_72, gd_75, \
                         gd_84, gd_87 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = f_32 * gd_0[k]
                  - f_32 * gd_3[k]
                  + f_28 * gd_18[k]
                  - f_28 * gd_21[k]
                  - f_33 * gd_30[k]
                  + f_33 * gd_33[k]
                  + f_32 * gd_60[k]
                  - f_32 * gd_63[k]
                  - f_33 * gd_72[k]
                  + f_33 * gd_75[k]
                  + f_34 * gd_84[k]
                  - f_34 * gd_87[k];
    }

#pragma omp simd aligned(gd_13, gd_16, gd_43, gd_46, gd_55, gd_58 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = -f_20 * gd_13[k]
                  - f_20 * gd_43[k]
                  + f_21 * gd_55[k];

        g_26[k] = -f_20 * gd_16[k]
                  - f_20 * gd_46[k]
                  + f_21 * gd_58[k];
    }

#pragma omp simd aligned(gd_12, gd_14, gd_15, gd_17, gd_42, gd_44, gd_45, gd_47, gd_54, gd_56, \
                         gd_57, gd_59 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = f_22 * gd_12[k]
                  + f_22 * gd_15[k]
                  - f_23 * gd_17[k]
                  + f_22 * gd_42[k]
                  + f_22 * gd_45[k]
                  - f_23 * gd_47[k]
                  - f_24 * gd_54[k]
                  - f_24 * gd_57[k]
                  + f_25 * gd_59[k];

        g_28[k] = -f_20 * gd_14[k]
                  - f_20 * gd_44[k]
                  + f_21 * gd_56[k];

        g_29[k] = -f_26 * gd_12[k]
                  + f_26 * gd_15[k]
                  - f_26 * gd_42[k]
                  + f_26 * gd_45[k]
                  + f_27 * gd_54[k]
                  - f_27 * gd_57[k];
    }

#pragma omp simd aligned(gd_1, gd_4, gd_31, gd_34, gd_61, gd_64, gd_73, \
                         gd_76 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = -f_18 * gd_1[k]
                  + f_19 * gd_31[k]
                  + f_18 * gd_61[k]
                  - f_19 * gd_73[k];

        g_31[k] = -f_18 * gd_4[k]
                  + f_19 * gd_34[k]
                  + f_18 * gd_64[k]
                  - f_19 * gd_76[k];
    }

#pragma omp simd aligned(gd_0, gd_3, gd_5, gd_30, gd_33, gd_35, gd_60, gd_63, gd_65, gd_72, \
                         gd_75, gd_77 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = f_35 * gd_0[k]
                  + f_35 * gd_3[k]
                  - f_14 * gd_5[k]
                  - f_36 * gd_30[k]
                  - f_36 * gd_33[k]
                  + f_16 * gd_35[k]
                  - f_35 * gd_60[k]
                  - f_35 * gd_63[k]
                  + f_14 * gd_65[k]
                  + f_36 * gd_72[k]
                  + f_36 * gd_75[k]
                  - f_16 * gd_77[k];
    }

#pragma omp simd aligned(gd_0, gd_2, gd_3, gd_30, gd_32, gd_33, gd_60, gd_62, gd_63, gd_72, \
                         gd_74, gd_75 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = -f_18 * gd_2[k]
                  + f_19 * gd_32[k]
                  + f_18 * gd_62[k]
                  - f_19 * gd_74[k];

        g_34[k] = -f_37 * gd_0[k]
                  + f_37 * gd_3[k]
                  + f_38 * gd_30[k]
                  - f_38 * gd_33[k]
                  + f_37 * gd_60[k]
                  - f_37 * gd_63[k]
                  - f_38 * gd_72[k]
                  + f_38 * gd_75[k];
    }

#pragma omp simd aligned(gd_12, gd_13, gd_14, gd_15, gd_16, gd_17, gd_42, gd_43, gd_44, gd_45, \
                         gd_46, gd_47 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = f_5 * gd_13[k]
                  - f_4 * gd_43[k];

        g_36[k] = f_5 * gd_16[k]
                  - f_4 * gd_46[k];

        g_37[k] = -f_8 * gd_12[k]
                  - f_8 * gd_15[k]
                  + f_9 * gd_17[k]
                  + f_6 * gd_42[k]
                  + f_6 * gd_45[k]
                  - f_7 * gd_47[k];

        g_38[k] = f_5 * gd_14[k]
                  - f_4 * gd_44[k];
    }

#pragma omp simd aligned(gd_1, gd_4, gd_12, gd_15, gd_19, gd_22, gd_42, gd_45, gd_61, \
                         gd_64 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_39[k] = f_11 * gd_12[k]
                  - f_11 * gd_15[k]
                  - f_10 * gd_42[k]
                  + f_10 * gd_45[k];

        g_40[k] = f_39 * gd_1[k]
                  - f_40 * gd_19[k]
                  + f_39 * gd_61[k];

        g_41[k] = f_39 * gd_4[k]
                  - f_40 * gd_22[k]
                  + f_39 * gd_64[k];
    }

#pragma omp simd aligned(gd_0, gd_2, gd_3, gd_5, gd_18, gd_20, gd_21, gd_23, gd_60, gd_62, \
                         gd_63, gd_65 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = -f_41 * gd_0[k]
                  - f_41 * gd_3[k]
                  + f_42 * gd_5[k]
                  + f_43 * gd_18[k]
                  + f_43 * gd_21[k]
                  - f_44 * gd_23[k]
                  - f_41 * gd_60[k]
                  - f_41 * gd_63[k]
                  + f_42 * gd_65[k];

        g_43[k] = f_39 * gd_2[k]
                  - f_40 * gd_20[k]
                  + f_39 * gd_62[k];

        g_44[k] = f_45 * gd_0[k]
                  - f_45 * gd_3[k]
                  - f_46 * gd_18[k]
                  + f_46 * gd_21[k]
                  + f_45 * gd_60[k]
                  - f_45 * gd_63[k];
    }
}

}  // namespace simdtrf
