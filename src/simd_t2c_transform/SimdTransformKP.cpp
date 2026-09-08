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


#include "SimdTransformKP.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_kp(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t kp,
             const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 0.21875 * std::sqrt(429.0);
    const auto f_1 = 1.09375 * std::sqrt(429.0);
    const auto f_2 = 0.65625 * std::sqrt(429.0);
    const auto f_3 = 0.03125 * std::sqrt(429.0);
    const auto f_4 = 0.1875 * std::sqrt(6006.0);
    const auto f_5 = 0.625 * std::sqrt(6006.0);
    const auto f_6 = 0.15625 * std::sqrt(231.0);
    const auto f_7 = 1.875 * std::sqrt(231.0);
    const auto f_8 = 0.28125 * std::sqrt(231.0);
    const auto f_9 = 3.75 * std::sqrt(231.0);
    const auto f_10 = 0.03125 * std::sqrt(231.0);
    const auto f_11 = 0.375 * std::sqrt(231.0);
    const auto f_12 = 0.75 * std::sqrt(231.0);
    const auto f_13 = 2.5 * std::sqrt(231.0);
    const auto f_14 = 0.28125 * std::sqrt(21.0);
    const auto f_15 = 0.46875 * std::sqrt(21.0);
    const auto f_16 = 5.625 * std::sqrt(21.0);
    const auto f_17 = 0.09375 * std::sqrt(21.0);
    const auto f_18 = 3.75 * std::sqrt(21.0);
    const auto f_19 = 7.5 * std::sqrt(21.0);
    const auto f_20 = 1.875 * std::sqrt(21.0);
    const auto f_21 = 2.5 * std::sqrt(21.0);
    const auto f_22 = 0.9375 * std::sqrt(42.0);
    const auto f_23 = 1.875 * std::sqrt(42.0);
    const auto f_24 = 5.0 * std::sqrt(42.0);
    const auto f_25 = 3.0 * std::sqrt(42.0);
    const auto f_26 = 0.15625 * std::sqrt(7.0);
    const auto f_27 = 0.46875 * std::sqrt(7.0);
    const auto f_28 = 3.75 * std::sqrt(7.0);
    const auto f_29 = 7.5 * std::sqrt(7.0);
    const auto f_30 = 2.0 * std::sqrt(7.0);
    const auto f_31 = 0.46875 * std::sqrt(42.0);
    const auto f_32 = 2.5 * std::sqrt(42.0);
    const auto f_33 = 1.5 * std::sqrt(42.0);
    const auto f_34 = 0.1875 * std::sqrt(231.0);
    const auto f_35 = 0.9375 * std::sqrt(231.0);
    const auto f_36 = 0.625 * std::sqrt(231.0);
    const auto f_37 = 0.03125 * std::sqrt(6006.0);
    const auto f_38 = 0.46875 * std::sqrt(6006.0);

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

    const auto *kp_0 = buffer.data(kp + 0);
    const auto *kp_1 = buffer.data(kp + 1);
    const auto *kp_2 = buffer.data(kp + 2);
    const auto *kp_3 = buffer.data(kp + 3);
    const auto *kp_4 = buffer.data(kp + 4);
    const auto *kp_5 = buffer.data(kp + 5);
    const auto *kp_6 = buffer.data(kp + 6);
    const auto *kp_7 = buffer.data(kp + 7);
    const auto *kp_8 = buffer.data(kp + 8);
    const auto *kp_9 = buffer.data(kp + 9);
    const auto *kp_10 = buffer.data(kp + 10);
    const auto *kp_11 = buffer.data(kp + 11);
    const auto *kp_12 = buffer.data(kp + 12);
    const auto *kp_13 = buffer.data(kp + 13);
    const auto *kp_14 = buffer.data(kp + 14);
    const auto *kp_15 = buffer.data(kp + 15);
    const auto *kp_16 = buffer.data(kp + 16);
    const auto *kp_17 = buffer.data(kp + 17);
    const auto *kp_18 = buffer.data(kp + 18);
    const auto *kp_19 = buffer.data(kp + 19);
    const auto *kp_20 = buffer.data(kp + 20);
    const auto *kp_21 = buffer.data(kp + 21);
    const auto *kp_22 = buffer.data(kp + 22);
    const auto *kp_23 = buffer.data(kp + 23);
    const auto *kp_24 = buffer.data(kp + 24);
    const auto *kp_25 = buffer.data(kp + 25);
    const auto *kp_26 = buffer.data(kp + 26);
    const auto *kp_27 = buffer.data(kp + 27);
    const auto *kp_28 = buffer.data(kp + 28);
    const auto *kp_29 = buffer.data(kp + 29);
    const auto *kp_30 = buffer.data(kp + 30);
    const auto *kp_31 = buffer.data(kp + 31);
    const auto *kp_32 = buffer.data(kp + 32);
    const auto *kp_33 = buffer.data(kp + 33);
    const auto *kp_34 = buffer.data(kp + 34);
    const auto *kp_35 = buffer.data(kp + 35);
    const auto *kp_36 = buffer.data(kp + 36);
    const auto *kp_37 = buffer.data(kp + 37);
    const auto *kp_38 = buffer.data(kp + 38);
    const auto *kp_39 = buffer.data(kp + 39);
    const auto *kp_40 = buffer.data(kp + 40);
    const auto *kp_41 = buffer.data(kp + 41);
    const auto *kp_42 = buffer.data(kp + 42);
    const auto *kp_43 = buffer.data(kp + 43);
    const auto *kp_44 = buffer.data(kp + 44);
    const auto *kp_45 = buffer.data(kp + 45);
    const auto *kp_46 = buffer.data(kp + 46);
    const auto *kp_47 = buffer.data(kp + 47);
    const auto *kp_48 = buffer.data(kp + 48);
    const auto *kp_49 = buffer.data(kp + 49);
    const auto *kp_50 = buffer.data(kp + 50);
    const auto *kp_51 = buffer.data(kp + 51);
    const auto *kp_52 = buffer.data(kp + 52);
    const auto *kp_53 = buffer.data(kp + 53);
    const auto *kp_54 = buffer.data(kp + 54);
    const auto *kp_55 = buffer.data(kp + 55);
    const auto *kp_56 = buffer.data(kp + 56);
    const auto *kp_57 = buffer.data(kp + 57);
    const auto *kp_58 = buffer.data(kp + 58);
    const auto *kp_59 = buffer.data(kp + 59);
    const auto *kp_60 = buffer.data(kp + 60);
    const auto *kp_61 = buffer.data(kp + 61);
    const auto *kp_62 = buffer.data(kp + 62);
    const auto *kp_63 = buffer.data(kp + 63);
    const auto *kp_64 = buffer.data(kp + 64);
    const auto *kp_65 = buffer.data(kp + 65);
    const auto *kp_66 = buffer.data(kp + 66);
    const auto *kp_67 = buffer.data(kp + 67);
    const auto *kp_68 = buffer.data(kp + 68);
    const auto *kp_69 = buffer.data(kp + 69);
    const auto *kp_70 = buffer.data(kp + 70);
    const auto *kp_71 = buffer.data(kp + 71);
    const auto *kp_72 = buffer.data(kp + 72);
    const auto *kp_73 = buffer.data(kp + 73);
    const auto *kp_74 = buffer.data(kp + 74);
    const auto *kp_75 = buffer.data(kp + 75);
    const auto *kp_76 = buffer.data(kp + 76);
    const auto *kp_77 = buffer.data(kp + 77);
    const auto *kp_78 = buffer.data(kp + 78);
    const auto *kp_79 = buffer.data(kp + 79);
    const auto *kp_80 = buffer.data(kp + 80);
    const auto *kp_81 = buffer.data(kp + 81);
    const auto *kp_82 = buffer.data(kp + 82);
    const auto *kp_83 = buffer.data(kp + 83);
    const auto *kp_84 = buffer.data(kp + 84);
    const auto *kp_85 = buffer.data(kp + 85);
    const auto *kp_86 = buffer.data(kp + 86);
    const auto *kp_87 = buffer.data(kp + 87);
    const auto *kp_88 = buffer.data(kp + 88);
    const auto *kp_89 = buffer.data(kp + 89);
    const auto *kp_90 = buffer.data(kp + 90);
    const auto *kp_91 = buffer.data(kp + 91);
    const auto *kp_92 = buffer.data(kp + 92);
    const auto *kp_93 = buffer.data(kp + 93);
    const auto *kp_94 = buffer.data(kp + 94);
    const auto *kp_95 = buffer.data(kp + 95);
    const auto *kp_96 = buffer.data(kp + 96);
    const auto *kp_97 = buffer.data(kp + 97);
    const auto *kp_98 = buffer.data(kp + 98);
    const auto *kp_99 = buffer.data(kp + 99);
    const auto *kp_100 = buffer.data(kp + 100);
    const auto *kp_101 = buffer.data(kp + 101);
    const auto *kp_102 = buffer.data(kp + 102);
    const auto *kp_103 = buffer.data(kp + 103);
    const auto *kp_104 = buffer.data(kp + 104);
    const auto *kp_105 = buffer.data(kp + 105);
    const auto *kp_106 = buffer.data(kp + 106);
    const auto *kp_107 = buffer.data(kp + 107);

#pragma omp simd aligned(kp_3, kp_4, kp_5, kp_18, kp_19, kp_20, kp_45, kp_46, kp_47, kp_84, \
                         kp_85, kp_86 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * kp_4[k]
                 - f_1 * kp_19[k]
                 + f_2 * kp_46[k]
                 - f_3 * kp_85[k];

        g_1[k] = f_0 * kp_5[k]
                 - f_1 * kp_20[k]
                 + f_2 * kp_47[k]
                 - f_3 * kp_86[k];

        g_2[k] = f_0 * kp_3[k]
                 - f_1 * kp_18[k]
                 + f_2 * kp_45[k]
                 - f_3 * kp_84[k];
    }

#pragma omp simd aligned(kp_12, kp_13, kp_14, kp_33, kp_34, kp_35, kp_66, kp_67, \
                         kp_68 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = f_4 * kp_13[k]
                 - f_5 * kp_34[k]
                 + f_4 * kp_67[k];

        g_4[k] = f_4 * kp_14[k]
                 - f_5 * kp_35[k]
                 + f_4 * kp_68[k];

        g_5[k] = f_4 * kp_12[k]
                 - f_5 * kp_33[k]
                 + f_4 * kp_66[k];
    }

#pragma omp simd aligned(kp_4, kp_5, kp_19, kp_20, kp_25, kp_26, kp_46, kp_47, kp_52, kp_53, \
                         kp_85, kp_86, kp_91, kp_92 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = -f_6 * kp_4[k]
                 + f_6 * kp_19[k]
                 + f_7 * kp_25[k]
                 + f_8 * kp_46[k]
                 - f_9 * kp_52[k]
                 - f_10 * kp_85[k]
                 + f_11 * kp_91[k];

        g_7[k] = -f_6 * kp_5[k]
                 + f_6 * kp_20[k]
                 + f_7 * kp_26[k]
                 + f_8 * kp_47[k]
                 - f_9 * kp_53[k]
                 - f_10 * kp_86[k]
                 + f_11 * kp_92[k];
    }

#pragma omp simd aligned(kp_3, kp_13, kp_18, kp_24, kp_40, kp_45, kp_51, kp_67, kp_73, kp_84, \
                         kp_90 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_8[k] = -f_6 * kp_3[k]
                 + f_6 * kp_18[k]
                 + f_7 * kp_24[k]
                 + f_8 * kp_45[k]
                 - f_9 * kp_51[k]
                 - f_10 * kp_84[k]
                 + f_11 * kp_90[k];

        g_9[k] = -f_12 * kp_13[k]
                 + f_13 * kp_40[k]
                 + f_12 * kp_67[k]
                 - f_13 * kp_73[k];
    }

#pragma omp simd aligned(kp_12, kp_14, kp_39, kp_41, kp_66, kp_68, kp_72, \
                         kp_74 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = -f_12 * kp_14[k]
                  + f_13 * kp_41[k]
                  + f_12 * kp_68[k]
                  - f_13 * kp_74[k];

        g_11[k] = -f_12 * kp_12[k]
                  + f_13 * kp_39[k]
                  + f_12 * kp_66[k]
                  - f_13 * kp_72[k];
    }

#pragma omp simd aligned(kp_4, kp_19, kp_25, kp_46, kp_52, kp_58, kp_85, kp_91, \
                         kp_97 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_12[k] = f_14 * kp_4[k]
                  + f_15 * kp_19[k]
                  - f_16 * kp_25[k]
                  + f_17 * kp_46[k]
                  - f_18 * kp_52[k]
                  + f_19 * kp_58[k]
                  - f_17 * kp_85[k]
                  + f_20 * kp_91[k]
                  - f_21 * kp_97[k];
    }

#pragma omp simd aligned(kp_5, kp_20, kp_26, kp_47, kp_53, kp_59, kp_86, kp_92, \
                         kp_98 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = f_14 * kp_5[k]
                  + f_15 * kp_20[k]
                  - f_16 * kp_26[k]
                  + f_17 * kp_47[k]
                  - f_18 * kp_53[k]
                  + f_19 * kp_59[k]
                  - f_17 * kp_86[k]
                  + f_20 * kp_92[k]
                  - f_21 * kp_98[k];
    }

#pragma omp simd aligned(kp_3, kp_18, kp_24, kp_45, kp_51, kp_57, kp_84, kp_90, \
                         kp_96 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = f_14 * kp_3[k]
                  + f_15 * kp_18[k]
                  - f_16 * kp_24[k]
                  + f_17 * kp_45[k]
                  - f_18 * kp_51[k]
                  + f_19 * kp_57[k]
                  - f_17 * kp_84[k]
                  + f_20 * kp_90[k]
                  - f_21 * kp_96[k];
    }

#pragma omp simd aligned(kp_13, kp_14, kp_34, kp_35, kp_40, kp_41, kp_67, kp_68, kp_73, kp_74, \
                         kp_79, kp_80 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_15[k] = f_22 * kp_13[k]
                  + f_23 * kp_34[k]
                  - f_24 * kp_40[k]
                  + f_22 * kp_67[k]
                  - f_24 * kp_73[k]
                  + f_25 * kp_79[k];

        g_16[k] = f_22 * kp_14[k]
                  + f_23 * kp_35[k]
                  - f_24 * kp_41[k]
                  + f_22 * kp_68[k]
                  - f_24 * kp_74[k]
                  + f_25 * kp_80[k];
    }

#pragma omp simd aligned(kp_12, kp_33, kp_39, kp_66, kp_72, kp_78 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = f_22 * kp_12[k]
                  + f_23 * kp_33[k]
                  - f_24 * kp_39[k]
                  + f_22 * kp_66[k]
                  - f_24 * kp_72[k]
                  + f_25 * kp_78[k];
    }

#pragma omp simd aligned(kp_4, kp_19, kp_25, kp_46, kp_52, kp_58, kp_85, kp_91, kp_97, \
                         kp_103 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = -f_26 * kp_4[k]
                  - f_27 * kp_19[k]
                  + f_28 * kp_25[k]
                  - f_27 * kp_46[k]
                  + f_29 * kp_52[k]
                  - f_29 * kp_58[k]
                  - f_26 * kp_85[k]
                  + f_28 * kp_91[k]
                  - f_29 * kp_97[k]
                  + f_30 * kp_103[k];
    }

#pragma omp simd aligned(kp_5, kp_20, kp_26, kp_47, kp_53, kp_59, kp_86, kp_92, kp_98, \
                         kp_104 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = -f_26 * kp_5[k]
                  - f_27 * kp_20[k]
                  + f_28 * kp_26[k]
                  - f_27 * kp_47[k]
                  + f_29 * kp_53[k]
                  - f_29 * kp_59[k]
                  - f_26 * kp_86[k]
                  + f_28 * kp_92[k]
                  - f_29 * kp_98[k]
                  + f_30 * kp_104[k];
    }

#pragma omp simd aligned(kp_3, kp_18, kp_24, kp_45, kp_51, kp_57, kp_84, kp_90, kp_96, \
                         kp_102 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = -f_26 * kp_3[k]
                  - f_27 * kp_18[k]
                  + f_28 * kp_24[k]
                  - f_27 * kp_45[k]
                  + f_29 * kp_51[k]
                  - f_29 * kp_57[k]
                  - f_26 * kp_84[k]
                  + f_28 * kp_90[k]
                  - f_29 * kp_96[k]
                  + f_30 * kp_102[k];
    }

#pragma omp simd aligned(kp_7, kp_22, kp_28, kp_49, kp_55, kp_61, kp_88, kp_94, kp_100, \
                         kp_106 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = -2.1875 * kp_7[k]
                  - 6.5625 * kp_22[k]
                  + 13.125 * kp_28[k]
                  - 6.5625 * kp_49[k]
                  + 26.25 * kp_55[k]
                  - 10.5 * kp_61[k]
                  - 2.1875 * kp_88[k]
                  + 13.125 * kp_94[k]
                  - 10.5 * kp_100[k]
                  + kp_106[k];
    }

#pragma omp simd aligned(kp_8, kp_23, kp_29, kp_50, kp_56, kp_62, kp_89, kp_95, kp_101, \
                         kp_107 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = -2.1875 * kp_8[k]
                  - 6.5625 * kp_23[k]
                  + 13.125 * kp_29[k]
                  - 6.5625 * kp_50[k]
                  + 26.25 * kp_56[k]
                  - 10.5 * kp_62[k]
                  - 2.1875 * kp_89[k]
                  + 13.125 * kp_95[k]
                  - 10.5 * kp_101[k]
                  + kp_107[k];
    }

#pragma omp simd aligned(kp_6, kp_21, kp_27, kp_48, kp_54, kp_60, kp_87, kp_93, kp_99, \
                         kp_105 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = -2.1875 * kp_6[k]
                  - 6.5625 * kp_21[k]
                  + 13.125 * kp_27[k]
                  - 6.5625 * kp_48[k]
                  + 26.25 * kp_54[k]
                  - 10.5 * kp_60[k]
                  - 2.1875 * kp_87[k]
                  + 13.125 * kp_93[k]
                  - 10.5 * kp_99[k]
                  + kp_105[k];
    }

#pragma omp simd aligned(kp_1, kp_10, kp_16, kp_31, kp_37, kp_43, kp_64, kp_70, kp_76, \
                         kp_82 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = -f_26 * kp_1[k]
                  - f_27 * kp_10[k]
                  + f_28 * kp_16[k]
                  - f_27 * kp_31[k]
                  + f_29 * kp_37[k]
                  - f_29 * kp_43[k]
                  - f_26 * kp_64[k]
                  + f_28 * kp_70[k]
                  - f_29 * kp_76[k]
                  + f_30 * kp_82[k];
    }

#pragma omp simd aligned(kp_2, kp_11, kp_17, kp_32, kp_38, kp_44, kp_65, kp_71, kp_77, \
                         kp_83 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = -f_26 * kp_2[k]
                  - f_27 * kp_11[k]
                  + f_28 * kp_17[k]
                  - f_27 * kp_32[k]
                  + f_29 * kp_38[k]
                  - f_29 * kp_44[k]
                  - f_26 * kp_65[k]
                  + f_28 * kp_71[k]
                  - f_29 * kp_77[k]
                  + f_30 * kp_83[k];
    }

#pragma omp simd aligned(kp_0, kp_9, kp_15, kp_30, kp_36, kp_42, kp_63, kp_69, kp_75, \
                         kp_81 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = -f_26 * kp_0[k]
                  - f_27 * kp_9[k]
                  + f_28 * kp_15[k]
                  - f_27 * kp_30[k]
                  + f_29 * kp_36[k]
                  - f_29 * kp_42[k]
                  - f_26 * kp_63[k]
                  + f_28 * kp_69[k]
                  - f_29 * kp_75[k]
                  + f_30 * kp_81[k];
    }

#pragma omp simd aligned(kp_7, kp_22, kp_28, kp_49, kp_61, kp_88, kp_94, \
                         kp_100 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = f_31 * kp_7[k]
                  + f_31 * kp_22[k]
                  - f_32 * kp_28[k]
                  - f_31 * kp_49[k]
                  + f_33 * kp_61[k]
                  - f_31 * kp_88[k]
                  + f_32 * kp_94[k]
                  - f_33 * kp_100[k];
    }

#pragma omp simd aligned(kp_8, kp_23, kp_29, kp_50, kp_62, kp_89, kp_95, \
                         kp_101 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = f_31 * kp_8[k]
                  + f_31 * kp_23[k]
                  - f_32 * kp_29[k]
                  - f_31 * kp_50[k]
                  + f_33 * kp_62[k]
                  - f_31 * kp_89[k]
                  + f_32 * kp_95[k]
                  - f_33 * kp_101[k];
    }

#pragma omp simd aligned(kp_6, kp_21, kp_27, kp_48, kp_60, kp_87, kp_93, \
                         kp_99 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = f_31 * kp_6[k]
                  + f_31 * kp_21[k]
                  - f_32 * kp_27[k]
                  - f_31 * kp_48[k]
                  + f_33 * kp_60[k]
                  - f_31 * kp_87[k]
                  + f_32 * kp_93[k]
                  - f_33 * kp_99[k];
    }

#pragma omp simd aligned(kp_1, kp_10, kp_16, kp_31, kp_37, kp_43, kp_64, kp_70, \
                         kp_76 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = f_17 * kp_1[k]
                  - f_17 * kp_10[k]
                  - f_20 * kp_16[k]
                  - f_15 * kp_31[k]
                  + f_18 * kp_37[k]
                  + f_21 * kp_43[k]
                  - f_14 * kp_64[k]
                  + f_16 * kp_70[k]
                  - f_19 * kp_76[k];
    }

#pragma omp simd aligned(kp_2, kp_11, kp_17, kp_32, kp_38, kp_44, kp_65, kp_71, \
                         kp_77 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = f_17 * kp_2[k]
                  - f_17 * kp_11[k]
                  - f_20 * kp_17[k]
                  - f_15 * kp_32[k]
                  + f_18 * kp_38[k]
                  + f_21 * kp_44[k]
                  - f_14 * kp_65[k]
                  + f_16 * kp_71[k]
                  - f_19 * kp_77[k];
    }

#pragma omp simd aligned(kp_0, kp_9, kp_15, kp_30, kp_36, kp_42, kp_63, kp_69, \
                         kp_75 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = f_17 * kp_0[k]
                  - f_17 * kp_9[k]
                  - f_20 * kp_15[k]
                  - f_15 * kp_30[k]
                  + f_18 * kp_36[k]
                  + f_21 * kp_42[k]
                  - f_14 * kp_63[k]
                  + f_16 * kp_69[k]
                  - f_19 * kp_75[k];
    }

#pragma omp simd aligned(kp_7, kp_8, kp_22, kp_23, kp_28, kp_29, kp_49, kp_50, kp_55, kp_56, \
                         kp_88, kp_89, kp_94, kp_95 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = -f_34 * kp_7[k]
                  + f_35 * kp_22[k]
                  + f_36 * kp_28[k]
                  + f_35 * kp_49[k]
                  - f_9 * kp_55[k]
                  - f_34 * kp_88[k]
                  + f_36 * kp_94[k];

        g_34[k] = -f_34 * kp_8[k]
                  + f_35 * kp_23[k]
                  + f_36 * kp_29[k]
                  + f_35 * kp_50[k]
                  - f_9 * kp_56[k]
                  - f_34 * kp_89[k]
                  + f_36 * kp_95[k];
    }

#pragma omp simd aligned(kp_1, kp_6, kp_10, kp_16, kp_21, kp_27, kp_31, kp_37, kp_48, kp_54, \
                         kp_64, kp_70, kp_87, kp_93 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = -f_34 * kp_6[k]
                  + f_35 * kp_21[k]
                  + f_36 * kp_27[k]
                  + f_35 * kp_48[k]
                  - f_9 * kp_54[k]
                  - f_34 * kp_87[k]
                  + f_36 * kp_93[k];

        g_36[k] = -f_10 * kp_1[k]
                  + f_8 * kp_10[k]
                  + f_11 * kp_16[k]
                  + f_6 * kp_31[k]
                  - f_9 * kp_37[k]
                  - f_6 * kp_64[k]
                  + f_7 * kp_70[k];
    }

#pragma omp simd aligned(kp_0, kp_2, kp_9, kp_11, kp_15, kp_17, kp_30, kp_32, kp_36, kp_38, \
                         kp_63, kp_65, kp_69, kp_71 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_37[k] = -f_10 * kp_2[k]
                  + f_8 * kp_11[k]
                  + f_11 * kp_17[k]
                  + f_6 * kp_32[k]
                  - f_9 * kp_38[k]
                  - f_6 * kp_65[k]
                  + f_7 * kp_71[k];

        g_38[k] = -f_10 * kp_0[k]
                  + f_8 * kp_9[k]
                  + f_11 * kp_15[k]
                  + f_6 * kp_30[k]
                  - f_9 * kp_36[k]
                  - f_6 * kp_63[k]
                  + f_7 * kp_69[k];
    }

#pragma omp simd aligned(kp_6, kp_7, kp_8, kp_21, kp_22, kp_23, kp_48, kp_49, kp_50, kp_87, \
                         kp_88, kp_89 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_39[k] = f_37 * kp_7[k]
                  - f_38 * kp_22[k]
                  + f_38 * kp_49[k]
                  - f_37 * kp_88[k];

        g_40[k] = f_37 * kp_8[k]
                  - f_38 * kp_23[k]
                  + f_38 * kp_50[k]
                  - f_37 * kp_89[k];

        g_41[k] = f_37 * kp_6[k]
                  - f_38 * kp_21[k]
                  + f_38 * kp_48[k]
                  - f_37 * kp_87[k];
    }

#pragma omp simd aligned(kp_0, kp_1, kp_2, kp_9, kp_10, kp_11, kp_30, kp_31, kp_32, kp_63, \
                         kp_64, kp_65 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = f_3 * kp_1[k]
                  - f_2 * kp_10[k]
                  + f_1 * kp_31[k]
                  - f_0 * kp_64[k];

        g_43[k] = f_3 * kp_2[k]
                  - f_2 * kp_11[k]
                  + f_1 * kp_32[k]
                  - f_0 * kp_65[k];

        g_44[k] = f_3 * kp_0[k]
                  - f_2 * kp_9[k]
                  + f_1 * kp_30[k]
                  - f_0 * kp_63[k];
    }
}

}  // namespace simdtrf
