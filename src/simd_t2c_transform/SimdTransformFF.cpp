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


#include "SimdTransformFF.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_ff(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t ff,
             const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 3.75 * std::sqrt(6.0);
    const auto f_1 = 1.25 * std::sqrt(6.0);
    const auto f_2 = 0.375 * std::sqrt(15.0);
    const auto f_3 = 1.5 * std::sqrt(15.0);
    const auto f_4 = 0.125 * std::sqrt(15.0);
    const auto f_5 = 0.5 * std::sqrt(15.0);
    const auto f_6 = 1.125 * std::sqrt(10.0);
    const auto f_7 = 0.75 * std::sqrt(10.0);
    const auto f_8 = 0.375 * std::sqrt(10.0);
    const auto f_9 = 0.25 * std::sqrt(10.0);
    const auto f_10 = 1.875 * std::sqrt(6.0);
    const auto f_11 = 0.625 * std::sqrt(6.0);
    const auto f_12 = 3.0 * std::sqrt(10.0);
    const auto f_13 = std::sqrt(15.0);
    const auto f_14 = 0.375 * std::sqrt(6.0);
    const auto f_15 = 0.25 * std::sqrt(6.0);
    const auto f_16 = 1.5 * std::sqrt(6.0);
    const auto f_17 = std::sqrt(6.0);
    const auto f_18 = 1.5 * std::sqrt(10.0);
    const auto f_19 = 0.75 * std::sqrt(15.0);

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
    auto *g_45 = values + 45 * nvalues;
    auto *g_46 = values + 46 * nvalues;
    auto *g_47 = values + 47 * nvalues;
    auto *g_48 = values + 48 * nvalues;

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_4 = buffer.data(ff + 4);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_6 = buffer.data(ff + 6);
    const auto *ff_7 = buffer.data(ff + 7);
    const auto *ff_8 = buffer.data(ff + 8);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_11 = buffer.data(ff + 11);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_14 = buffer.data(ff + 14);
    const auto *ff_15 = buffer.data(ff + 15);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_17 = buffer.data(ff + 17);
    const auto *ff_18 = buffer.data(ff + 18);
    const auto *ff_19 = buffer.data(ff + 19);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_21 = buffer.data(ff + 21);
    const auto *ff_22 = buffer.data(ff + 22);
    const auto *ff_23 = buffer.data(ff + 23);
    const auto *ff_24 = buffer.data(ff + 24);
    const auto *ff_25 = buffer.data(ff + 25);
    const auto *ff_26 = buffer.data(ff + 26);
    const auto *ff_27 = buffer.data(ff + 27);
    const auto *ff_28 = buffer.data(ff + 28);
    const auto *ff_29 = buffer.data(ff + 29);
    const auto *ff_30 = buffer.data(ff + 30);
    const auto *ff_31 = buffer.data(ff + 31);
    const auto *ff_32 = buffer.data(ff + 32);
    const auto *ff_33 = buffer.data(ff + 33);
    const auto *ff_34 = buffer.data(ff + 34);
    const auto *ff_35 = buffer.data(ff + 35);
    const auto *ff_36 = buffer.data(ff + 36);
    const auto *ff_37 = buffer.data(ff + 37);
    const auto *ff_38 = buffer.data(ff + 38);
    const auto *ff_39 = buffer.data(ff + 39);
    const auto *ff_40 = buffer.data(ff + 40);
    const auto *ff_41 = buffer.data(ff + 41);
    const auto *ff_42 = buffer.data(ff + 42);
    const auto *ff_43 = buffer.data(ff + 43);
    const auto *ff_44 = buffer.data(ff + 44);
    const auto *ff_45 = buffer.data(ff + 45);
    const auto *ff_46 = buffer.data(ff + 46);
    const auto *ff_47 = buffer.data(ff + 47);
    const auto *ff_48 = buffer.data(ff + 48);
    const auto *ff_49 = buffer.data(ff + 49);
    const auto *ff_50 = buffer.data(ff + 50);
    const auto *ff_51 = buffer.data(ff + 51);
    const auto *ff_52 = buffer.data(ff + 52);
    const auto *ff_53 = buffer.data(ff + 53);
    const auto *ff_54 = buffer.data(ff + 54);
    const auto *ff_55 = buffer.data(ff + 55);
    const auto *ff_56 = buffer.data(ff + 56);
    const auto *ff_57 = buffer.data(ff + 57);
    const auto *ff_58 = buffer.data(ff + 58);
    const auto *ff_59 = buffer.data(ff + 59);
    const auto *ff_60 = buffer.data(ff + 60);
    const auto *ff_61 = buffer.data(ff + 61);
    const auto *ff_62 = buffer.data(ff + 62);
    const auto *ff_63 = buffer.data(ff + 63);
    const auto *ff_64 = buffer.data(ff + 64);
    const auto *ff_65 = buffer.data(ff + 65);
    const auto *ff_66 = buffer.data(ff + 66);
    const auto *ff_67 = buffer.data(ff + 67);
    const auto *ff_68 = buffer.data(ff + 68);
    const auto *ff_69 = buffer.data(ff + 69);
    const auto *ff_70 = buffer.data(ff + 70);
    const auto *ff_71 = buffer.data(ff + 71);
    const auto *ff_72 = buffer.data(ff + 72);
    const auto *ff_73 = buffer.data(ff + 73);
    const auto *ff_74 = buffer.data(ff + 74);
    const auto *ff_75 = buffer.data(ff + 75);
    const auto *ff_76 = buffer.data(ff + 76);
    const auto *ff_77 = buffer.data(ff + 77);
    const auto *ff_78 = buffer.data(ff + 78);
    const auto *ff_79 = buffer.data(ff + 79);
    const auto *ff_80 = buffer.data(ff + 80);
    const auto *ff_81 = buffer.data(ff + 81);
    const auto *ff_82 = buffer.data(ff + 82);
    const auto *ff_83 = buffer.data(ff + 83);
    const auto *ff_84 = buffer.data(ff + 84);
    const auto *ff_85 = buffer.data(ff + 85);
    const auto *ff_86 = buffer.data(ff + 86);
    const auto *ff_87 = buffer.data(ff + 87);
    const auto *ff_88 = buffer.data(ff + 88);
    const auto *ff_89 = buffer.data(ff + 89);
    const auto *ff_90 = buffer.data(ff + 90);
    const auto *ff_91 = buffer.data(ff + 91);
    const auto *ff_92 = buffer.data(ff + 92);
    const auto *ff_93 = buffer.data(ff + 93);
    const auto *ff_94 = buffer.data(ff + 94);
    const auto *ff_95 = buffer.data(ff + 95);
    const auto *ff_96 = buffer.data(ff + 96);
    const auto *ff_97 = buffer.data(ff + 97);
    const auto *ff_98 = buffer.data(ff + 98);
    const auto *ff_99 = buffer.data(ff + 99);

#pragma omp simd aligned(ff_11, ff_14, ff_16, ff_18, ff_61, ff_64, ff_66, \
                         ff_68 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = 5.625 * ff_11[k]
                 - 1.875 * ff_16[k]
                 - 1.875 * ff_61[k]
                 + 0.625 * ff_66[k];

        g_1[k] = f_0 * ff_14[k]
                 - f_1 * ff_64[k];

        g_2[k] = -f_2 * ff_11[k]
                 - f_2 * ff_16[k]
                 + f_3 * ff_18[k]
                 + f_4 * ff_61[k]
                 + f_4 * ff_66[k]
                 - f_5 * ff_68[k];
    }

#pragma omp simd aligned(ff_10, ff_12, ff_13, ff_15, ff_17, ff_19, ff_60, ff_62, ff_63, ff_65, \
                         ff_67, ff_69 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = -f_6 * ff_12[k]
                 - f_6 * ff_17[k]
                 + f_7 * ff_19[k]
                 + f_8 * ff_62[k]
                 + f_8 * ff_67[k]
                 - f_9 * ff_69[k];

        g_4[k] = -f_2 * ff_10[k]
                 - f_2 * ff_13[k]
                 + f_3 * ff_15[k]
                 + f_4 * ff_60[k]
                 + f_4 * ff_63[k]
                 - f_5 * ff_65[k];

        g_5[k] = f_10 * ff_12[k]
                 - f_10 * ff_17[k]
                 - f_11 * ff_62[k]
                 + f_11 * ff_67[k];

        g_6[k] = 1.875 * ff_10[k]
                 - 5.625 * ff_13[k]
                 - 0.625 * ff_60[k]
                 + 1.875 * ff_63[k];
    }

#pragma omp simd aligned(ff_40, ff_41, ff_42, ff_43, ff_44, ff_45, ff_46, ff_47, ff_48, \
                         ff_49 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = f_0 * ff_41[k]
                 - f_1 * ff_46[k];

        g_8[k] = 15.0 * ff_44[k];

        g_9[k] = -f_7 * ff_41[k]
                 - f_7 * ff_46[k]
                 + f_12 * ff_48[k];

        g_10[k] = -f_3 * ff_42[k]
                  - f_3 * ff_47[k]
                  + f_13 * ff_49[k];

        g_11[k] = -f_7 * ff_40[k]
                  - f_7 * ff_43[k]
                  + f_12 * ff_45[k];

        g_12[k] = 7.5 * ff_42[k]
                  - 7.5 * ff_47[k];
    }

#pragma omp simd aligned(ff_11, ff_14, ff_16, ff_40, ff_43, ff_61, ff_64, ff_66, ff_81, ff_84, \
                         ff_86 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = f_1 * ff_40[k]
                  - f_0 * ff_43[k];

        g_14[k] = -f_2 * ff_11[k]
                  + f_4 * ff_16[k]
                  - f_2 * ff_61[k]
                  + f_4 * ff_66[k]
                  + f_3 * ff_81[k]
                  - f_5 * ff_86[k];

        g_15[k] = -f_7 * ff_14[k]
                  - f_7 * ff_64[k]
                  + f_12 * ff_84[k];
    }

#pragma omp simd aligned(ff_11, ff_16, ff_18, ff_61, ff_66, ff_68, ff_81, ff_86, \
                         ff_88 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_16[k] = 0.375 * ff_11[k]
                  + 0.375 * ff_16[k]
                  - 1.5 * ff_18[k]
                  + 0.375 * ff_61[k]
                  + 0.375 * ff_66[k]
                  - 1.5 * ff_68[k]
                  - 1.5 * ff_81[k]
                  - 1.5 * ff_86[k]
                  + 6.0 * ff_88[k];
    }

#pragma omp simd aligned(ff_12, ff_17, ff_19, ff_62, ff_67, ff_69, ff_82, ff_87, \
                         ff_89 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = f_14 * ff_12[k]
                  + f_14 * ff_17[k]
                  - f_15 * ff_19[k]
                  + f_14 * ff_62[k]
                  + f_14 * ff_67[k]
                  - f_15 * ff_69[k]
                  - f_16 * ff_82[k]
                  - f_16 * ff_87[k]
                  + f_17 * ff_89[k];
    }

#pragma omp simd aligned(ff_10, ff_13, ff_15, ff_60, ff_63, ff_65, ff_80, ff_83, \
                         ff_85 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = 0.375 * ff_10[k]
                  + 0.375 * ff_13[k]
                  - 1.5 * ff_15[k]
                  + 0.375 * ff_60[k]
                  + 0.375 * ff_63[k]
                  - 1.5 * ff_65[k]
                  - 1.5 * ff_80[k]
                  - 1.5 * ff_83[k]
                  + 6.0 * ff_85[k];
    }

#pragma omp simd aligned(ff_10, ff_12, ff_13, ff_17, ff_60, ff_62, ff_63, ff_67, ff_80, ff_82, \
                         ff_83, ff_87 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = -f_8 * ff_12[k]
                  + f_8 * ff_17[k]
                  - f_8 * ff_62[k]
                  + f_8 * ff_67[k]
                  + f_18 * ff_82[k]
                  - f_18 * ff_87[k];

        g_20[k] = -f_4 * ff_10[k]
                  + f_2 * ff_13[k]
                  - f_4 * ff_60[k]
                  + f_2 * ff_63[k]
                  + f_5 * ff_80[k]
                  - f_3 * ff_83[k];
    }

#pragma omp simd aligned(ff_21, ff_24, ff_26, ff_28, ff_71, ff_74, ff_76, ff_78, ff_91, ff_94, \
                         ff_96, ff_98 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = -f_6 * ff_21[k]
                  + f_8 * ff_26[k]
                  - f_6 * ff_71[k]
                  + f_8 * ff_76[k]
                  + f_7 * ff_91[k]
                  - f_9 * ff_96[k];

        g_22[k] = -f_3 * ff_24[k]
                  - f_3 * ff_74[k]
                  + f_13 * ff_94[k];

        g_23[k] = f_14 * ff_21[k]
                  + f_14 * ff_26[k]
                  - f_16 * ff_28[k]
                  + f_14 * ff_71[k]
                  + f_14 * ff_76[k]
                  - f_16 * ff_78[k]
                  - f_15 * ff_91[k]
                  - f_15 * ff_96[k]
                  + f_17 * ff_98[k];
    }

#pragma omp simd aligned(ff_22, ff_27, ff_29, ff_72, ff_77, ff_79, ff_92, ff_97, \
                         ff_99 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = 2.25 * ff_22[k]
                  + 2.25 * ff_27[k]
                  - 1.5 * ff_29[k]
                  + 2.25 * ff_72[k]
                  + 2.25 * ff_77[k]
                  - 1.5 * ff_79[k]
                  - 1.5 * ff_92[k]
                  - 1.5 * ff_97[k]
                  + ff_99[k];
    }

#pragma omp simd aligned(ff_20, ff_23, ff_25, ff_70, ff_73, ff_75, ff_90, ff_93, \
                         ff_95 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = f_14 * ff_20[k]
                  + f_14 * ff_23[k]
                  - f_16 * ff_25[k]
                  + f_14 * ff_70[k]
                  + f_14 * ff_73[k]
                  - f_16 * ff_75[k]
                  - f_15 * ff_90[k]
                  - f_15 * ff_93[k]
                  + f_17 * ff_95[k];
    }

#pragma omp simd aligned(ff_20, ff_22, ff_23, ff_27, ff_70, ff_72, ff_73, ff_77, ff_90, ff_92, \
                         ff_93, ff_97 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = -f_19 * ff_22[k]
                  + f_19 * ff_27[k]
                  - f_19 * ff_72[k]
                  + f_19 * ff_77[k]
                  + f_5 * ff_92[k]
                  - f_5 * ff_97[k];

        g_27[k] = -f_8 * ff_20[k]
                  + f_6 * ff_23[k]
                  - f_8 * ff_70[k]
                  + f_6 * ff_73[k]
                  + f_9 * ff_90[k]
                  - f_7 * ff_93[k];
    }

#pragma omp simd aligned(ff_1, ff_4, ff_6, ff_8, ff_31, ff_34, ff_36, ff_38, ff_51, ff_54, \
                         ff_56, ff_58 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = -f_2 * ff_1[k]
                  + f_4 * ff_6[k]
                  - f_2 * ff_31[k]
                  + f_4 * ff_36[k]
                  + f_3 * ff_51[k]
                  - f_5 * ff_56[k];

        g_29[k] = -f_7 * ff_4[k]
                  - f_7 * ff_34[k]
                  + f_12 * ff_54[k];

        g_30[k] = 0.375 * ff_1[k]
                  + 0.375 * ff_6[k]
                  - 1.5 * ff_8[k]
                  + 0.375 * ff_31[k]
                  + 0.375 * ff_36[k]
                  - 1.5 * ff_38[k]
                  - 1.5 * ff_51[k]
                  - 1.5 * ff_56[k]
                  + 6.0 * ff_58[k];
    }

#pragma omp simd aligned(ff_2, ff_7, ff_9, ff_32, ff_37, ff_39, ff_52, ff_57, \
                         ff_59 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = f_14 * ff_2[k]
                  + f_14 * ff_7[k]
                  - f_15 * ff_9[k]
                  + f_14 * ff_32[k]
                  + f_14 * ff_37[k]
                  - f_15 * ff_39[k]
                  - f_16 * ff_52[k]
                  - f_16 * ff_57[k]
                  + f_17 * ff_59[k];
    }

#pragma omp simd aligned(ff_0, ff_3, ff_5, ff_30, ff_33, ff_35, ff_50, ff_53, \
                         ff_55 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = 0.375 * ff_0[k]
                  + 0.375 * ff_3[k]
                  - 1.5 * ff_5[k]
                  + 0.375 * ff_30[k]
                  + 0.375 * ff_33[k]
                  - 1.5 * ff_35[k]
                  - 1.5 * ff_50[k]
                  - 1.5 * ff_53[k]
                  + 6.0 * ff_55[k];
    }

#pragma omp simd aligned(ff_0, ff_2, ff_3, ff_7, ff_30, ff_32, ff_33, ff_37, ff_50, ff_52, \
                         ff_53, ff_57 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = -f_8 * ff_2[k]
                  + f_8 * ff_7[k]
                  - f_8 * ff_32[k]
                  + f_8 * ff_37[k]
                  + f_18 * ff_52[k]
                  - f_18 * ff_57[k];

        g_34[k] = -f_4 * ff_0[k]
                  + f_2 * ff_3[k]
                  - f_4 * ff_30[k]
                  + f_2 * ff_33[k]
                  + f_5 * ff_50[k]
                  - f_3 * ff_53[k];
    }

#pragma omp simd aligned(ff_21, ff_24, ff_26, ff_28, ff_71, ff_74, ff_76, \
                         ff_78 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = f_10 * ff_21[k]
                  - f_11 * ff_26[k]
                  - f_10 * ff_71[k]
                  + f_11 * ff_76[k];

        g_36[k] = 7.5 * ff_24[k]
                  - 7.5 * ff_74[k];

        g_37[k] = -f_8 * ff_21[k]
                  - f_8 * ff_26[k]
                  + f_18 * ff_28[k]
                  + f_8 * ff_71[k]
                  + f_8 * ff_76[k]
                  - f_18 * ff_78[k];
    }

#pragma omp simd aligned(ff_20, ff_22, ff_23, ff_25, ff_27, ff_29, ff_70, ff_72, ff_73, ff_75, \
                         ff_77, ff_79 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = -f_19 * ff_22[k]
                  - f_19 * ff_27[k]
                  + f_5 * ff_29[k]
                  + f_19 * ff_72[k]
                  + f_19 * ff_77[k]
                  - f_5 * ff_79[k];

        g_39[k] = -f_8 * ff_20[k]
                  - f_8 * ff_23[k]
                  + f_18 * ff_25[k]
                  + f_8 * ff_70[k]
                  + f_8 * ff_73[k]
                  - f_18 * ff_75[k];

        g_40[k] = 3.75 * ff_22[k]
                  - 3.75 * ff_27[k]
                  - 3.75 * ff_72[k]
                  + 3.75 * ff_77[k];

        g_41[k] = f_11 * ff_20[k]
                  - f_10 * ff_23[k]
                  - f_11 * ff_70[k]
                  + f_10 * ff_73[k];
    }

#pragma omp simd aligned(ff_1, ff_4, ff_6, ff_8, ff_31, ff_34, ff_36, \
                         ff_38 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = 1.875 * ff_1[k]
                  - 0.625 * ff_6[k]
                  - 5.625 * ff_31[k]
                  + 1.875 * ff_36[k];

        g_43[k] = f_1 * ff_4[k]
                  - f_0 * ff_34[k];

        g_44[k] = -f_4 * ff_1[k]
                  - f_4 * ff_6[k]
                  + f_5 * ff_8[k]
                  + f_2 * ff_31[k]
                  + f_2 * ff_36[k]
                  - f_3 * ff_38[k];
    }

#pragma omp simd aligned(ff_0, ff_2, ff_3, ff_5, ff_7, ff_9, ff_30, ff_32, ff_33, ff_35, \
                         ff_37, ff_39 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_45[k] = -f_8 * ff_2[k]
                  - f_8 * ff_7[k]
                  + f_9 * ff_9[k]
                  + f_6 * ff_32[k]
                  + f_6 * ff_37[k]
                  - f_7 * ff_39[k];

        g_46[k] = -f_4 * ff_0[k]
                  - f_4 * ff_3[k]
                  + f_5 * ff_5[k]
                  + f_2 * ff_30[k]
                  + f_2 * ff_33[k]
                  - f_3 * ff_35[k];

        g_47[k] = f_11 * ff_2[k]
                  - f_11 * ff_7[k]
                  - f_10 * ff_32[k]
                  + f_10 * ff_37[k];

        g_48[k] = 0.625 * ff_0[k]
                  - 1.875 * ff_3[k]
                  - 1.875 * ff_30[k]
                  + 5.625 * ff_33[k];
    }
}

auto
transform_ff_tri(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t ff,
                 const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 3.75 * std::sqrt(6.0);
    const auto f_1 = 1.25 * std::sqrt(6.0);
    const auto f_2 = 0.375 * std::sqrt(15.0);
    const auto f_3 = 1.5 * std::sqrt(15.0);
    const auto f_4 = 0.125 * std::sqrt(15.0);
    const auto f_5 = 0.5 * std::sqrt(15.0);
    const auto f_6 = 1.125 * std::sqrt(10.0);
    const auto f_7 = 0.75 * std::sqrt(10.0);
    const auto f_8 = 0.375 * std::sqrt(10.0);
    const auto f_9 = 0.25 * std::sqrt(10.0);
    const auto f_10 = 1.875 * std::sqrt(6.0);
    const auto f_11 = 0.625 * std::sqrt(6.0);
    const auto f_12 = 3.0 * std::sqrt(10.0);
    const auto f_13 = std::sqrt(15.0);
    const auto f_14 = 0.375 * std::sqrt(6.0);
    const auto f_15 = 0.25 * std::sqrt(6.0);
    const auto f_16 = 1.5 * std::sqrt(6.0);
    const auto f_17 = std::sqrt(6.0);
    const auto f_18 = 1.5 * std::sqrt(10.0);
    const auto f_19 = 0.75 * std::sqrt(15.0);

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
    auto *g_45 = values + 45 * nvalues;
    auto *g_46 = values + 46 * nvalues;
    auto *g_47 = values + 47 * nvalues;
    auto *g_48 = values + 48 * nvalues;

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_7 = buffer.data(ff + 7);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_11 = buffer.data(ff + 11);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_14 = buffer.data(ff + 14);
    const auto *ff_15 = buffer.data(ff + 15);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_17 = buffer.data(ff + 17);
    const auto *ff_18 = buffer.data(ff + 18);
    const auto *ff_19 = buffer.data(ff + 19);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_22 = buffer.data(ff + 22);
    const auto *ff_23 = buffer.data(ff + 23);
    const auto *ff_25 = buffer.data(ff + 25);
    const auto *ff_27 = buffer.data(ff + 27);
    const auto *ff_29 = buffer.data(ff + 29);
    const auto *ff_30 = buffer.data(ff + 30);
    const auto *ff_32 = buffer.data(ff + 32);
    const auto *ff_33 = buffer.data(ff + 33);
    const auto *ff_35 = buffer.data(ff + 35);
    const auto *ff_37 = buffer.data(ff + 37);
    const auto *ff_40 = buffer.data(ff + 40);
    const auto *ff_41 = buffer.data(ff + 41);
    const auto *ff_42 = buffer.data(ff + 42);
    const auto *ff_43 = buffer.data(ff + 43);
    const auto *ff_44 = buffer.data(ff + 44);
    const auto *ff_45 = buffer.data(ff + 45);
    const auto *ff_46 = buffer.data(ff + 46);
    const auto *ff_47 = buffer.data(ff + 47);
    const auto *ff_48 = buffer.data(ff + 48);
    const auto *ff_49 = buffer.data(ff + 49);
    const auto *ff_50 = buffer.data(ff + 50);
    const auto *ff_52 = buffer.data(ff + 52);
    const auto *ff_53 = buffer.data(ff + 53);
    const auto *ff_55 = buffer.data(ff + 55);
    const auto *ff_57 = buffer.data(ff + 57);
    const auto *ff_60 = buffer.data(ff + 60);
    const auto *ff_61 = buffer.data(ff + 61);
    const auto *ff_62 = buffer.data(ff + 62);
    const auto *ff_63 = buffer.data(ff + 63);
    const auto *ff_64 = buffer.data(ff + 64);
    const auto *ff_65 = buffer.data(ff + 65);
    const auto *ff_66 = buffer.data(ff + 66);
    const auto *ff_67 = buffer.data(ff + 67);
    const auto *ff_68 = buffer.data(ff + 68);
    const auto *ff_69 = buffer.data(ff + 69);
    const auto *ff_70 = buffer.data(ff + 70);
    const auto *ff_72 = buffer.data(ff + 72);
    const auto *ff_73 = buffer.data(ff + 73);
    const auto *ff_75 = buffer.data(ff + 75);
    const auto *ff_77 = buffer.data(ff + 77);
    const auto *ff_79 = buffer.data(ff + 79);
    const auto *ff_80 = buffer.data(ff + 80);
    const auto *ff_81 = buffer.data(ff + 81);
    const auto *ff_82 = buffer.data(ff + 82);
    const auto *ff_83 = buffer.data(ff + 83);
    const auto *ff_85 = buffer.data(ff + 85);
    const auto *ff_86 = buffer.data(ff + 86);
    const auto *ff_87 = buffer.data(ff + 87);
    const auto *ff_88 = buffer.data(ff + 88);
    const auto *ff_89 = buffer.data(ff + 89);
    const auto *ff_90 = buffer.data(ff + 90);
    const auto *ff_92 = buffer.data(ff + 92);
    const auto *ff_93 = buffer.data(ff + 93);
    const auto *ff_95 = buffer.data(ff + 95);
    const auto *ff_97 = buffer.data(ff + 97);
    const auto *ff_99 = buffer.data(ff + 99);

#pragma omp simd aligned(ff_11, ff_14, ff_16, ff_18, ff_61, ff_64, ff_66, \
                         ff_68 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = 5.625 * ff_11[k]
                 - 1.875 * ff_16[k]
                 - 1.875 * ff_61[k]
                 + 0.625 * ff_66[k];

        g_1[k] = f_0 * ff_14[k]
                 - f_1 * ff_64[k];
        g_7[k] = g_1[k];

        g_2[k] = -f_2 * ff_11[k]
                 - f_2 * ff_16[k]
                 + f_3 * ff_18[k]
                 + f_4 * ff_61[k]
                 + f_4 * ff_66[k]
                 - f_5 * ff_68[k];
        g_14[k] = g_2[k];
    }

#pragma omp simd aligned(ff_10, ff_12, ff_13, ff_15, ff_17, ff_19, ff_60, ff_62, ff_63, ff_65, \
                         ff_67, ff_69 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = -f_6 * ff_12[k]
                 - f_6 * ff_17[k]
                 + f_7 * ff_19[k]
                 + f_8 * ff_62[k]
                 + f_8 * ff_67[k]
                 - f_9 * ff_69[k];
        g_21[k] = g_3[k];

        g_4[k] = -f_2 * ff_10[k]
                 - f_2 * ff_13[k]
                 + f_3 * ff_15[k]
                 + f_4 * ff_60[k]
                 + f_4 * ff_63[k]
                 - f_5 * ff_65[k];
        g_28[k] = g_4[k];

        g_5[k] = f_10 * ff_12[k]
                 - f_10 * ff_17[k]
                 - f_11 * ff_62[k]
                 + f_11 * ff_67[k];
        g_35[k] = g_5[k];

        g_6[k] = 1.875 * ff_10[k]
                 - 5.625 * ff_13[k]
                 - 0.625 * ff_60[k]
                 + 1.875 * ff_63[k];
        g_42[k] = g_6[k];
    }

#pragma omp simd aligned(ff_40, ff_41, ff_42, ff_43, ff_44, ff_45, ff_46, ff_47, ff_48, \
                         ff_49 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_8[k] = 15.0 * ff_44[k];

        g_9[k] = -f_7 * ff_41[k]
                 - f_7 * ff_46[k]
                 + f_12 * ff_48[k];
        g_15[k] = g_9[k];

        g_10[k] = -f_3 * ff_42[k]
                  - f_3 * ff_47[k]
                  + f_13 * ff_49[k];
        g_22[k] = g_10[k];

        g_11[k] = -f_7 * ff_40[k]
                  - f_7 * ff_43[k]
                  + f_12 * ff_45[k];
        g_29[k] = g_11[k];

        g_12[k] = 7.5 * ff_42[k]
                  - 7.5 * ff_47[k];
        g_36[k] = g_12[k];

        g_13[k] = f_1 * ff_40[k]
                  - f_0 * ff_43[k];
        g_43[k] = g_13[k];
    }

#pragma omp simd aligned(ff_11, ff_16, ff_18, ff_61, ff_66, ff_68, ff_81, ff_86, \
                         ff_88 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_16[k] = 0.375 * ff_11[k]
                  + 0.375 * ff_16[k]
                  - 1.5 * ff_18[k]
                  + 0.375 * ff_61[k]
                  + 0.375 * ff_66[k]
                  - 1.5 * ff_68[k]
                  - 1.5 * ff_81[k]
                  - 1.5 * ff_86[k]
                  + 6.0 * ff_88[k];
    }

#pragma omp simd aligned(ff_12, ff_17, ff_19, ff_62, ff_67, ff_69, ff_82, ff_87, \
                         ff_89 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = f_14 * ff_12[k]
                  + f_14 * ff_17[k]
                  - f_15 * ff_19[k]
                  + f_14 * ff_62[k]
                  + f_14 * ff_67[k]
                  - f_15 * ff_69[k]
                  - f_16 * ff_82[k]
                  - f_16 * ff_87[k]
                  + f_17 * ff_89[k];
        g_23[k] = g_17[k];
    }

#pragma omp simd aligned(ff_10, ff_13, ff_15, ff_60, ff_63, ff_65, ff_80, ff_83, \
                         ff_85 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = 0.375 * ff_10[k]
                  + 0.375 * ff_13[k]
                  - 1.5 * ff_15[k]
                  + 0.375 * ff_60[k]
                  + 0.375 * ff_63[k]
                  - 1.5 * ff_65[k]
                  - 1.5 * ff_80[k]
                  - 1.5 * ff_83[k]
                  + 6.0 * ff_85[k];
        g_30[k] = g_18[k];
    }

#pragma omp simd aligned(ff_10, ff_12, ff_13, ff_17, ff_60, ff_62, ff_63, ff_67, ff_80, ff_82, \
                         ff_83, ff_87 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = -f_8 * ff_12[k]
                  + f_8 * ff_17[k]
                  - f_8 * ff_62[k]
                  + f_8 * ff_67[k]
                  + f_18 * ff_82[k]
                  - f_18 * ff_87[k];
        g_37[k] = g_19[k];

        g_20[k] = -f_4 * ff_10[k]
                  + f_2 * ff_13[k]
                  - f_4 * ff_60[k]
                  + f_2 * ff_63[k]
                  + f_5 * ff_80[k]
                  - f_3 * ff_83[k];
        g_44[k] = g_20[k];
    }

#pragma omp simd aligned(ff_22, ff_27, ff_29, ff_72, ff_77, ff_79, ff_92, ff_97, \
                         ff_99 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = 2.25 * ff_22[k]
                  + 2.25 * ff_27[k]
                  - 1.5 * ff_29[k]
                  + 2.25 * ff_72[k]
                  + 2.25 * ff_77[k]
                  - 1.5 * ff_79[k]
                  - 1.5 * ff_92[k]
                  - 1.5 * ff_97[k]
                  + ff_99[k];
    }

#pragma omp simd aligned(ff_20, ff_23, ff_25, ff_70, ff_73, ff_75, ff_90, ff_93, \
                         ff_95 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = f_14 * ff_20[k]
                  + f_14 * ff_23[k]
                  - f_16 * ff_25[k]
                  + f_14 * ff_70[k]
                  + f_14 * ff_73[k]
                  - f_16 * ff_75[k]
                  - f_15 * ff_90[k]
                  - f_15 * ff_93[k]
                  + f_17 * ff_95[k];
        g_31[k] = g_25[k];
    }

#pragma omp simd aligned(ff_20, ff_22, ff_23, ff_27, ff_70, ff_72, ff_73, ff_77, ff_90, ff_92, \
                         ff_93, ff_97 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = -f_19 * ff_22[k]
                  + f_19 * ff_27[k]
                  - f_19 * ff_72[k]
                  + f_19 * ff_77[k]
                  + f_5 * ff_92[k]
                  - f_5 * ff_97[k];
        g_38[k] = g_26[k];

        g_27[k] = -f_8 * ff_20[k]
                  + f_6 * ff_23[k]
                  - f_8 * ff_70[k]
                  + f_6 * ff_73[k]
                  + f_9 * ff_90[k]
                  - f_7 * ff_93[k];
        g_45[k] = g_27[k];
    }

#pragma omp simd aligned(ff_0, ff_3, ff_5, ff_30, ff_33, ff_35, ff_50, ff_53, \
                         ff_55 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = 0.375 * ff_0[k]
                  + 0.375 * ff_3[k]
                  - 1.5 * ff_5[k]
                  + 0.375 * ff_30[k]
                  + 0.375 * ff_33[k]
                  - 1.5 * ff_35[k]
                  - 1.5 * ff_50[k]
                  - 1.5 * ff_53[k]
                  + 6.0 * ff_55[k];
    }

#pragma omp simd aligned(ff_0, ff_2, ff_3, ff_7, ff_30, ff_32, ff_33, ff_37, ff_50, ff_52, \
                         ff_53, ff_57 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = -f_8 * ff_2[k]
                  + f_8 * ff_7[k]
                  - f_8 * ff_32[k]
                  + f_8 * ff_37[k]
                  + f_18 * ff_52[k]
                  - f_18 * ff_57[k];
        g_39[k] = g_33[k];

        g_34[k] = -f_4 * ff_0[k]
                  + f_2 * ff_3[k]
                  - f_4 * ff_30[k]
                  + f_2 * ff_33[k]
                  + f_5 * ff_50[k]
                  - f_3 * ff_53[k];
        g_46[k] = g_34[k];
    }

#pragma omp simd aligned(ff_0, ff_3, ff_20, ff_22, ff_23, ff_27, ff_30, ff_33, ff_70, ff_72, \
                         ff_73, ff_77 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_40[k] = 3.75 * ff_22[k]
                  - 3.75 * ff_27[k]
                  - 3.75 * ff_72[k]
                  + 3.75 * ff_77[k];

        g_41[k] = f_11 * ff_20[k]
                  - f_10 * ff_23[k]
                  - f_11 * ff_70[k]
                  + f_10 * ff_73[k];
        g_47[k] = g_41[k];

        g_48[k] = 0.625 * ff_0[k]
                  - 1.875 * ff_3[k]
                  - 1.875 * ff_30[k]
                  + 5.625 * ff_33[k];
    }
}

}  // namespace simdtrf
