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


#include "SimdTransformLP.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_lp(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t lp,
             const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 0.1875 * std::sqrt(715.0);
    const auto f_1 = 1.3125 * std::sqrt(715.0);
    const auto f_2 = 0.65625 * std::sqrt(715.0);
    const auto f_3 = 3.28125 * std::sqrt(715.0);
    const auto f_4 = 1.96875 * std::sqrt(715.0);
    const auto f_5 = 0.09375 * std::sqrt(715.0);
    const auto f_6 = 0.09375 * std::sqrt(858.0);
    const auto f_7 = 0.21875 * std::sqrt(858.0);
    const auto f_8 = 1.3125 * std::sqrt(858.0);
    const auto f_9 = 4.375 * std::sqrt(858.0);
    const auto f_10 = 0.46875 * std::sqrt(1001.0);
    const auto f_11 = 1.875 * std::sqrt(1001.0);
    const auto f_12 = 0.84375 * std::sqrt(1001.0);
    const auto f_13 = 3.75 * std::sqrt(1001.0);
    const auto f_14 = 0.09375 * std::sqrt(1001.0);
    const auto f_15 = 0.375 * std::sqrt(1001.0);
    const auto f_16 = 0.1875 * std::sqrt(77.0);
    const auto f_17 = 4.5 * std::sqrt(77.0);
    const auto f_18 = 7.5 * std::sqrt(77.0);
    const auto f_19 = 0.28125 * std::sqrt(1155.0);
    const auto f_20 = 0.46875 * std::sqrt(1155.0);
    const auto f_21 = 1.875 * std::sqrt(1155.0);
    const auto f_22 = 0.09375 * std::sqrt(1155.0);
    const auto f_23 = 1.25 * std::sqrt(1155.0);
    const auto f_24 = 1.5 * std::sqrt(1155.0);
    const auto f_25 = 0.625 * std::sqrt(1155.0);
    const auto f_26 = 0.5 * std::sqrt(1155.0);
    const auto f_27 = 0.09375 * std::sqrt(70.0);
    const auto f_28 = 0.28125 * std::sqrt(70.0);
    const auto f_29 = 2.8125 * std::sqrt(70.0);
    const auto f_30 = 5.625 * std::sqrt(70.0);
    const auto f_31 = 7.5 * std::sqrt(70.0);
    const auto f_32 = 3.0 * std::sqrt(70.0);
    const auto f_33 = 0.046875 * std::sqrt(70.0);
    const auto f_34 = 1.40625 * std::sqrt(70.0);
    const auto f_35 = 3.75 * std::sqrt(70.0);
    const auto f_36 = 1.5 * std::sqrt(70.0);
    const auto f_37 = 0.046875 * std::sqrt(77.0);
    const auto f_38 = 1.125 * std::sqrt(77.0);
    const auto f_39 = 0.46875 * std::sqrt(77.0);
    const auto f_40 = 5.625 * std::sqrt(77.0);
    const auto f_41 = 1.875 * std::sqrt(77.0);
    const auto f_42 = 11.25 * std::sqrt(77.0);
    const auto f_43 = 0.015625 * std::sqrt(858.0);
    const auto f_44 = 3.28125 * std::sqrt(858.0);
    const auto f_45 = 0.0234375 * std::sqrt(715.0);
    const auto f_46 = 1.640625 * std::sqrt(715.0);

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
    auto *g_49 = values + 49 * nvalues;
    auto *g_50 = values + 50 * nvalues;

    const auto *lp_0 = buffer.data(lp + 0);
    const auto *lp_1 = buffer.data(lp + 1);
    const auto *lp_2 = buffer.data(lp + 2);
    const auto *lp_3 = buffer.data(lp + 3);
    const auto *lp_4 = buffer.data(lp + 4);
    const auto *lp_5 = buffer.data(lp + 5);
    const auto *lp_6 = buffer.data(lp + 6);
    const auto *lp_7 = buffer.data(lp + 7);
    const auto *lp_8 = buffer.data(lp + 8);
    const auto *lp_9 = buffer.data(lp + 9);
    const auto *lp_10 = buffer.data(lp + 10);
    const auto *lp_11 = buffer.data(lp + 11);
    const auto *lp_12 = buffer.data(lp + 12);
    const auto *lp_13 = buffer.data(lp + 13);
    const auto *lp_14 = buffer.data(lp + 14);
    const auto *lp_15 = buffer.data(lp + 15);
    const auto *lp_16 = buffer.data(lp + 16);
    const auto *lp_17 = buffer.data(lp + 17);
    const auto *lp_18 = buffer.data(lp + 18);
    const auto *lp_19 = buffer.data(lp + 19);
    const auto *lp_20 = buffer.data(lp + 20);
    const auto *lp_21 = buffer.data(lp + 21);
    const auto *lp_22 = buffer.data(lp + 22);
    const auto *lp_23 = buffer.data(lp + 23);
    const auto *lp_24 = buffer.data(lp + 24);
    const auto *lp_25 = buffer.data(lp + 25);
    const auto *lp_26 = buffer.data(lp + 26);
    const auto *lp_27 = buffer.data(lp + 27);
    const auto *lp_28 = buffer.data(lp + 28);
    const auto *lp_29 = buffer.data(lp + 29);
    const auto *lp_30 = buffer.data(lp + 30);
    const auto *lp_31 = buffer.data(lp + 31);
    const auto *lp_32 = buffer.data(lp + 32);
    const auto *lp_33 = buffer.data(lp + 33);
    const auto *lp_34 = buffer.data(lp + 34);
    const auto *lp_35 = buffer.data(lp + 35);
    const auto *lp_36 = buffer.data(lp + 36);
    const auto *lp_37 = buffer.data(lp + 37);
    const auto *lp_38 = buffer.data(lp + 38);
    const auto *lp_39 = buffer.data(lp + 39);
    const auto *lp_40 = buffer.data(lp + 40);
    const auto *lp_41 = buffer.data(lp + 41);
    const auto *lp_42 = buffer.data(lp + 42);
    const auto *lp_43 = buffer.data(lp + 43);
    const auto *lp_44 = buffer.data(lp + 44);
    const auto *lp_45 = buffer.data(lp + 45);
    const auto *lp_46 = buffer.data(lp + 46);
    const auto *lp_47 = buffer.data(lp + 47);
    const auto *lp_48 = buffer.data(lp + 48);
    const auto *lp_49 = buffer.data(lp + 49);
    const auto *lp_50 = buffer.data(lp + 50);
    const auto *lp_51 = buffer.data(lp + 51);
    const auto *lp_52 = buffer.data(lp + 52);
    const auto *lp_53 = buffer.data(lp + 53);
    const auto *lp_54 = buffer.data(lp + 54);
    const auto *lp_55 = buffer.data(lp + 55);
    const auto *lp_56 = buffer.data(lp + 56);
    const auto *lp_57 = buffer.data(lp + 57);
    const auto *lp_58 = buffer.data(lp + 58);
    const auto *lp_59 = buffer.data(lp + 59);
    const auto *lp_60 = buffer.data(lp + 60);
    const auto *lp_61 = buffer.data(lp + 61);
    const auto *lp_62 = buffer.data(lp + 62);
    const auto *lp_63 = buffer.data(lp + 63);
    const auto *lp_64 = buffer.data(lp + 64);
    const auto *lp_65 = buffer.data(lp + 65);
    const auto *lp_66 = buffer.data(lp + 66);
    const auto *lp_67 = buffer.data(lp + 67);
    const auto *lp_68 = buffer.data(lp + 68);
    const auto *lp_69 = buffer.data(lp + 69);
    const auto *lp_70 = buffer.data(lp + 70);
    const auto *lp_71 = buffer.data(lp + 71);
    const auto *lp_72 = buffer.data(lp + 72);
    const auto *lp_73 = buffer.data(lp + 73);
    const auto *lp_74 = buffer.data(lp + 74);
    const auto *lp_75 = buffer.data(lp + 75);
    const auto *lp_76 = buffer.data(lp + 76);
    const auto *lp_77 = buffer.data(lp + 77);
    const auto *lp_78 = buffer.data(lp + 78);
    const auto *lp_79 = buffer.data(lp + 79);
    const auto *lp_80 = buffer.data(lp + 80);
    const auto *lp_81 = buffer.data(lp + 81);
    const auto *lp_82 = buffer.data(lp + 82);
    const auto *lp_83 = buffer.data(lp + 83);
    const auto *lp_84 = buffer.data(lp + 84);
    const auto *lp_85 = buffer.data(lp + 85);
    const auto *lp_86 = buffer.data(lp + 86);
    const auto *lp_87 = buffer.data(lp + 87);
    const auto *lp_88 = buffer.data(lp + 88);
    const auto *lp_89 = buffer.data(lp + 89);
    const auto *lp_90 = buffer.data(lp + 90);
    const auto *lp_91 = buffer.data(lp + 91);
    const auto *lp_92 = buffer.data(lp + 92);
    const auto *lp_93 = buffer.data(lp + 93);
    const auto *lp_94 = buffer.data(lp + 94);
    const auto *lp_95 = buffer.data(lp + 95);
    const auto *lp_96 = buffer.data(lp + 96);
    const auto *lp_97 = buffer.data(lp + 97);
    const auto *lp_98 = buffer.data(lp + 98);
    const auto *lp_99 = buffer.data(lp + 99);
    const auto *lp_100 = buffer.data(lp + 100);
    const auto *lp_101 = buffer.data(lp + 101);
    const auto *lp_102 = buffer.data(lp + 102);
    const auto *lp_103 = buffer.data(lp + 103);
    const auto *lp_104 = buffer.data(lp + 104);
    const auto *lp_105 = buffer.data(lp + 105);
    const auto *lp_106 = buffer.data(lp + 106);
    const auto *lp_107 = buffer.data(lp + 107);
    const auto *lp_108 = buffer.data(lp + 108);
    const auto *lp_109 = buffer.data(lp + 109);
    const auto *lp_110 = buffer.data(lp + 110);
    const auto *lp_111 = buffer.data(lp + 111);
    const auto *lp_112 = buffer.data(lp + 112);
    const auto *lp_113 = buffer.data(lp + 113);
    const auto *lp_114 = buffer.data(lp + 114);
    const auto *lp_115 = buffer.data(lp + 115);
    const auto *lp_116 = buffer.data(lp + 116);
    const auto *lp_117 = buffer.data(lp + 117);
    const auto *lp_118 = buffer.data(lp + 118);
    const auto *lp_119 = buffer.data(lp + 119);
    const auto *lp_120 = buffer.data(lp + 120);
    const auto *lp_121 = buffer.data(lp + 121);
    const auto *lp_122 = buffer.data(lp + 122);
    const auto *lp_123 = buffer.data(lp + 123);
    const auto *lp_124 = buffer.data(lp + 124);
    const auto *lp_125 = buffer.data(lp + 125);
    const auto *lp_126 = buffer.data(lp + 126);
    const auto *lp_127 = buffer.data(lp + 127);
    const auto *lp_128 = buffer.data(lp + 128);
    const auto *lp_129 = buffer.data(lp + 129);
    const auto *lp_130 = buffer.data(lp + 130);
    const auto *lp_131 = buffer.data(lp + 131);
    const auto *lp_132 = buffer.data(lp + 132);
    const auto *lp_133 = buffer.data(lp + 133);
    const auto *lp_134 = buffer.data(lp + 134);

#pragma omp simd aligned(lp_3, lp_4, lp_5, lp_18, lp_19, lp_20, lp_45, lp_46, lp_47, lp_84, \
                         lp_85, lp_86 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * lp_4[k]
                 - f_1 * lp_19[k]
                 + f_1 * lp_46[k]
                 - f_0 * lp_85[k];

        g_1[k] = f_0 * lp_5[k]
                 - f_1 * lp_20[k]
                 + f_1 * lp_47[k]
                 - f_0 * lp_86[k];

        g_2[k] = f_0 * lp_3[k]
                 - f_1 * lp_18[k]
                 + f_1 * lp_45[k]
                 - f_0 * lp_84[k];
    }

#pragma omp simd aligned(lp_12, lp_13, lp_14, lp_33, lp_34, lp_35, lp_66, lp_67, lp_68, \
                         lp_111, lp_112, lp_113 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = f_2 * lp_13[k]
                 - f_3 * lp_34[k]
                 + f_4 * lp_67[k]
                 - f_5 * lp_112[k];

        g_4[k] = f_2 * lp_14[k]
                 - f_3 * lp_35[k]
                 + f_4 * lp_68[k]
                 - f_5 * lp_113[k];

        g_5[k] = f_2 * lp_12[k]
                 - f_3 * lp_33[k]
                 + f_4 * lp_66[k]
                 - f_5 * lp_111[k];
    }

#pragma omp simd aligned(lp_4, lp_5, lp_19, lp_20, lp_25, lp_26, lp_46, lp_47, lp_52, lp_53, \
                         lp_85, lp_86, lp_91, lp_92 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = -f_6 * lp_4[k]
                 + f_7 * lp_19[k]
                 + f_8 * lp_25[k]
                 + f_7 * lp_46[k]
                 - f_9 * lp_52[k]
                 - f_6 * lp_85[k]
                 + f_8 * lp_91[k];

        g_7[k] = -f_6 * lp_5[k]
                 + f_7 * lp_20[k]
                 + f_8 * lp_26[k]
                 + f_7 * lp_47[k]
                 - f_9 * lp_53[k]
                 - f_6 * lp_86[k]
                 + f_8 * lp_92[k];
    }

#pragma omp simd aligned(lp_3, lp_13, lp_18, lp_24, lp_34, lp_40, lp_45, lp_51, lp_67, lp_73, \
                         lp_84, lp_90, lp_112, lp_118 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_8[k] = -f_6 * lp_3[k]
                 + f_7 * lp_18[k]
                 + f_8 * lp_24[k]
                 + f_7 * lp_45[k]
                 - f_9 * lp_51[k]
                 - f_6 * lp_84[k]
                 + f_8 * lp_90[k];

        g_9[k] = -f_10 * lp_13[k]
                 + f_10 * lp_34[k]
                 + f_11 * lp_40[k]
                 + f_12 * lp_67[k]
                 - f_13 * lp_73[k]
                 - f_14 * lp_112[k]
                 + f_15 * lp_118[k];
    }

#pragma omp simd aligned(lp_12, lp_14, lp_33, lp_35, lp_39, lp_41, lp_66, lp_68, lp_72, lp_74, \
                         lp_111, lp_113, lp_117, lp_119 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = -f_10 * lp_14[k]
                  + f_10 * lp_35[k]
                  + f_11 * lp_41[k]
                  + f_12 * lp_68[k]
                  - f_13 * lp_74[k]
                  - f_14 * lp_113[k]
                  + f_15 * lp_119[k];

        g_11[k] = -f_10 * lp_12[k]
                  + f_10 * lp_33[k]
                  + f_11 * lp_39[k]
                  + f_12 * lp_66[k]
                  - f_13 * lp_72[k]
                  - f_14 * lp_111[k]
                  + f_15 * lp_117[k];
    }

#pragma omp simd aligned(lp_4, lp_19, lp_25, lp_46, lp_58, lp_85, lp_91, \
                         lp_97 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_12[k] = f_16 * lp_4[k]
                  + f_16 * lp_19[k]
                  - f_17 * lp_25[k]
                  - f_16 * lp_46[k]
                  + f_18 * lp_58[k]
                  - f_16 * lp_85[k]
                  + f_17 * lp_91[k]
                  - f_18 * lp_97[k];
    }

#pragma omp simd aligned(lp_5, lp_20, lp_26, lp_47, lp_59, lp_86, lp_92, \
                         lp_98 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = f_16 * lp_5[k]
                  + f_16 * lp_20[k]
                  - f_17 * lp_26[k]
                  - f_16 * lp_47[k]
                  + f_18 * lp_59[k]
                  - f_16 * lp_86[k]
                  + f_17 * lp_92[k]
                  - f_18 * lp_98[k];
    }

#pragma omp simd aligned(lp_3, lp_18, lp_24, lp_45, lp_57, lp_84, lp_90, \
                         lp_96 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = f_16 * lp_3[k]
                  + f_16 * lp_18[k]
                  - f_17 * lp_24[k]
                  - f_16 * lp_45[k]
                  + f_18 * lp_57[k]
                  - f_16 * lp_84[k]
                  + f_17 * lp_90[k]
                  - f_18 * lp_96[k];
    }

#pragma omp simd aligned(lp_13, lp_34, lp_40, lp_67, lp_73, lp_79, lp_112, lp_118, \
                         lp_124 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_15[k] = f_19 * lp_13[k]
                  + f_20 * lp_34[k]
                  - f_21 * lp_40[k]
                  + f_22 * lp_67[k]
                  - f_23 * lp_73[k]
                  + f_24 * lp_79[k]
                  - f_22 * lp_112[k]
                  + f_25 * lp_118[k]
                  - f_26 * lp_124[k];
    }

#pragma omp simd aligned(lp_14, lp_35, lp_41, lp_68, lp_74, lp_80, lp_113, lp_119, \
                         lp_125 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_16[k] = f_19 * lp_14[k]
                  + f_20 * lp_35[k]
                  - f_21 * lp_41[k]
                  + f_22 * lp_68[k]
                  - f_23 * lp_74[k]
                  + f_24 * lp_80[k]
                  - f_22 * lp_113[k]
                  + f_25 * lp_119[k]
                  - f_26 * lp_125[k];
    }

#pragma omp simd aligned(lp_12, lp_33, lp_39, lp_66, lp_72, lp_78, lp_111, lp_117, \
                         lp_123 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = f_19 * lp_12[k]
                  + f_20 * lp_33[k]
                  - f_21 * lp_39[k]
                  + f_22 * lp_66[k]
                  - f_23 * lp_72[k]
                  + f_24 * lp_78[k]
                  - f_22 * lp_111[k]
                  + f_25 * lp_117[k]
                  - f_26 * lp_123[k];
    }

#pragma omp simd aligned(lp_4, lp_19, lp_25, lp_46, lp_52, lp_58, lp_85, lp_91, lp_97, \
                         lp_103 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = -f_27 * lp_4[k]
                  - f_28 * lp_19[k]
                  + f_29 * lp_25[k]
                  - f_28 * lp_46[k]
                  + f_30 * lp_52[k]
                  - f_31 * lp_58[k]
                  - f_27 * lp_85[k]
                  + f_29 * lp_91[k]
                  - f_31 * lp_97[k]
                  + f_32 * lp_103[k];
    }

#pragma omp simd aligned(lp_5, lp_20, lp_26, lp_47, lp_53, lp_59, lp_86, lp_92, lp_98, \
                         lp_104 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = -f_27 * lp_5[k]
                  - f_28 * lp_20[k]
                  + f_29 * lp_26[k]
                  - f_28 * lp_47[k]
                  + f_30 * lp_53[k]
                  - f_31 * lp_59[k]
                  - f_27 * lp_86[k]
                  + f_29 * lp_92[k]
                  - f_31 * lp_98[k]
                  + f_32 * lp_104[k];
    }

#pragma omp simd aligned(lp_3, lp_18, lp_24, lp_45, lp_51, lp_57, lp_84, lp_90, lp_96, \
                         lp_102 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = -f_27 * lp_3[k]
                  - f_28 * lp_18[k]
                  + f_29 * lp_24[k]
                  - f_28 * lp_45[k]
                  + f_30 * lp_51[k]
                  - f_31 * lp_57[k]
                  - f_27 * lp_84[k]
                  + f_29 * lp_90[k]
                  - f_31 * lp_96[k]
                  + f_32 * lp_102[k];
    }

#pragma omp simd aligned(lp_13, lp_34, lp_40, lp_67, lp_73, lp_79, lp_112, lp_118, lp_124, \
                         lp_130 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = -3.28125 * lp_13[k]
                  - 9.84375 * lp_34[k]
                  + 26.25 * lp_40[k]
                  - 9.84375 * lp_67[k]
                  + 52.5 * lp_73[k]
                  - 31.5 * lp_79[k]
                  - 3.28125 * lp_112[k]
                  + 26.25 * lp_118[k]
                  - 31.5 * lp_124[k]
                  + 6.0 * lp_130[k];
    }

#pragma omp simd aligned(lp_14, lp_35, lp_41, lp_68, lp_74, lp_80, lp_113, lp_119, lp_125, \
                         lp_131 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = -3.28125 * lp_14[k]
                  - 9.84375 * lp_35[k]
                  + 26.25 * lp_41[k]
                  - 9.84375 * lp_68[k]
                  + 52.5 * lp_74[k]
                  - 31.5 * lp_80[k]
                  - 3.28125 * lp_113[k]
                  + 26.25 * lp_119[k]
                  - 31.5 * lp_125[k]
                  + 6.0 * lp_131[k];
    }

#pragma omp simd aligned(lp_12, lp_33, lp_39, lp_66, lp_72, lp_78, lp_111, lp_117, lp_123, \
                         lp_129 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = -3.28125 * lp_12[k]
                  - 9.84375 * lp_33[k]
                  + 26.25 * lp_39[k]
                  - 9.84375 * lp_66[k]
                  + 52.5 * lp_72[k]
                  - 31.5 * lp_78[k]
                  - 3.28125 * lp_111[k]
                  + 26.25 * lp_117[k]
                  - 31.5 * lp_123[k]
                  + 6.0 * lp_129[k];
    }

#pragma omp simd aligned(lp_1, lp_10, lp_16, lp_31, lp_37, lp_43, lp_64, lp_70, lp_76, lp_82, \
                         lp_109, lp_115, lp_121, lp_127, lp_133 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = 0.2734375 * lp_1[k]
                  + 1.09375 * lp_10[k]
                  - 8.75 * lp_16[k]
                  + 1.640625 * lp_31[k]
                  - 26.25 * lp_37[k]
                  + 26.25 * lp_43[k]
                  + 1.09375 * lp_64[k]
                  - 26.25 * lp_70[k]
                  + 52.5 * lp_76[k]
                  - 14.0 * lp_82[k]
                  + 0.2734375 * lp_109[k]
                  - 8.75 * lp_115[k]
                  + 26.25 * lp_121[k]
                  - 14.0 * lp_127[k]
                  + lp_133[k];
    }

#pragma omp simd aligned(lp_2, lp_11, lp_17, lp_32, lp_38, lp_44, lp_65, lp_71, lp_77, lp_83, \
                         lp_110, lp_116, lp_122, lp_128, lp_134 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = 0.2734375 * lp_2[k]
                  + 1.09375 * lp_11[k]
                  - 8.75 * lp_17[k]
                  + 1.640625 * lp_32[k]
                  - 26.25 * lp_38[k]
                  + 26.25 * lp_44[k]
                  + 1.09375 * lp_65[k]
                  - 26.25 * lp_71[k]
                  + 52.5 * lp_77[k]
                  - 14.0 * lp_83[k]
                  + 0.2734375 * lp_110[k]
                  - 8.75 * lp_116[k]
                  + 26.25 * lp_122[k]
                  - 14.0 * lp_128[k]
                  + lp_134[k];
    }

#pragma omp simd aligned(lp_0, lp_9, lp_15, lp_30, lp_36, lp_42, lp_63, lp_69, lp_75, lp_81, \
                         lp_108, lp_114, lp_120, lp_126, lp_132 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = 0.2734375 * lp_0[k]
                  + 1.09375 * lp_9[k]
                  - 8.75 * lp_15[k]
                  + 1.640625 * lp_30[k]
                  - 26.25 * lp_36[k]
                  + 26.25 * lp_42[k]
                  + 1.09375 * lp_63[k]
                  - 26.25 * lp_69[k]
                  + 52.5 * lp_75[k]
                  - 14.0 * lp_81[k]
                  + 0.2734375 * lp_108[k]
                  - 8.75 * lp_114[k]
                  + 26.25 * lp_120[k]
                  - 14.0 * lp_126[k]
                  + lp_132[k];
    }

#pragma omp simd aligned(lp_7, lp_22, lp_28, lp_49, lp_55, lp_61, lp_88, lp_94, lp_100, \
                         lp_106 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = -3.28125 * lp_7[k]
                  - 9.84375 * lp_22[k]
                  + 26.25 * lp_28[k]
                  - 9.84375 * lp_49[k]
                  + 52.5 * lp_55[k]
                  - 31.5 * lp_61[k]
                  - 3.28125 * lp_88[k]
                  + 26.25 * lp_94[k]
                  - 31.5 * lp_100[k]
                  + 6.0 * lp_106[k];
    }

#pragma omp simd aligned(lp_8, lp_23, lp_29, lp_50, lp_56, lp_62, lp_89, lp_95, lp_101, \
                         lp_107 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = -3.28125 * lp_8[k]
                  - 9.84375 * lp_23[k]
                  + 26.25 * lp_29[k]
                  - 9.84375 * lp_50[k]
                  + 52.5 * lp_56[k]
                  - 31.5 * lp_62[k]
                  - 3.28125 * lp_89[k]
                  + 26.25 * lp_95[k]
                  - 31.5 * lp_101[k]
                  + 6.0 * lp_107[k];
    }

#pragma omp simd aligned(lp_6, lp_21, lp_27, lp_48, lp_54, lp_60, lp_87, lp_93, lp_99, \
                         lp_105 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = -3.28125 * lp_6[k]
                  - 9.84375 * lp_21[k]
                  + 26.25 * lp_27[k]
                  - 9.84375 * lp_48[k]
                  + 52.5 * lp_54[k]
                  - 31.5 * lp_60[k]
                  - 3.28125 * lp_87[k]
                  + 26.25 * lp_93[k]
                  - 31.5 * lp_99[k]
                  + 6.0 * lp_105[k];
    }

#pragma omp simd aligned(lp_1, lp_10, lp_16, lp_37, lp_43, lp_64, lp_70, lp_82, lp_109, \
                         lp_115, lp_121, lp_127 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = -f_33 * lp_1[k]
                  - f_27 * lp_10[k]
                  + f_34 * lp_16[k]
                  + f_34 * lp_37[k]
                  - f_35 * lp_43[k]
                  + f_27 * lp_64[k]
                  - f_34 * lp_70[k]
                  + f_36 * lp_82[k]
                  + f_33 * lp_109[k]
                  - f_34 * lp_115[k]
                  + f_35 * lp_121[k]
                  - f_36 * lp_127[k];
    }

#pragma omp simd aligned(lp_2, lp_11, lp_17, lp_38, lp_44, lp_65, lp_71, lp_83, lp_110, \
                         lp_116, lp_122, lp_128 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = -f_33 * lp_2[k]
                  - f_27 * lp_11[k]
                  + f_34 * lp_17[k]
                  + f_34 * lp_38[k]
                  - f_35 * lp_44[k]
                  + f_27 * lp_65[k]
                  - f_34 * lp_71[k]
                  + f_36 * lp_83[k]
                  + f_33 * lp_110[k]
                  - f_34 * lp_116[k]
                  + f_35 * lp_122[k]
                  - f_36 * lp_128[k];
    }

#pragma omp simd aligned(lp_0, lp_9, lp_15, lp_36, lp_42, lp_63, lp_69, lp_81, lp_108, lp_114, \
                         lp_120, lp_126 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = -f_33 * lp_0[k]
                  - f_27 * lp_9[k]
                  + f_34 * lp_15[k]
                  + f_34 * lp_36[k]
                  - f_35 * lp_42[k]
                  + f_27 * lp_63[k]
                  - f_34 * lp_69[k]
                  + f_36 * lp_81[k]
                  + f_33 * lp_108[k]
                  - f_34 * lp_114[k]
                  + f_35 * lp_120[k]
                  - f_36 * lp_126[k];
    }

#pragma omp simd aligned(lp_7, lp_22, lp_28, lp_49, lp_55, lp_61, lp_88, lp_94, \
                         lp_100 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = f_22 * lp_7[k]
                  - f_22 * lp_22[k]
                  - f_25 * lp_28[k]
                  - f_20 * lp_49[k]
                  + f_23 * lp_55[k]
                  + f_26 * lp_61[k]
                  - f_19 * lp_88[k]
                  + f_21 * lp_94[k]
                  - f_24 * lp_100[k];
    }

#pragma omp simd aligned(lp_8, lp_23, lp_29, lp_50, lp_56, lp_62, lp_89, lp_95, \
                         lp_101 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = f_22 * lp_8[k]
                  - f_22 * lp_23[k]
                  - f_25 * lp_29[k]
                  - f_20 * lp_50[k]
                  + f_23 * lp_56[k]
                  + f_26 * lp_62[k]
                  - f_19 * lp_89[k]
                  + f_21 * lp_95[k]
                  - f_24 * lp_101[k];
    }

#pragma omp simd aligned(lp_6, lp_21, lp_27, lp_48, lp_54, lp_60, lp_87, lp_93, \
                         lp_99 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = f_22 * lp_6[k]
                  - f_22 * lp_21[k]
                  - f_25 * lp_27[k]
                  - f_20 * lp_48[k]
                  + f_23 * lp_54[k]
                  + f_26 * lp_60[k]
                  - f_19 * lp_87[k]
                  + f_21 * lp_93[k]
                  - f_24 * lp_99[k];
    }

#pragma omp simd aligned(lp_1, lp_10, lp_16, lp_31, lp_37, lp_43, lp_64, lp_70, lp_76, lp_109, \
                         lp_115, lp_121 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_36[k] = f_37 * lp_1[k]
                  - f_16 * lp_10[k]
                  - f_38 * lp_16[k]
                  - f_39 * lp_31[k]
                  + f_40 * lp_37[k]
                  + f_41 * lp_43[k]
                  - f_16 * lp_64[k]
                  + f_40 * lp_70[k]
                  - f_42 * lp_76[k]
                  + f_37 * lp_109[k]
                  - f_38 * lp_115[k]
                  + f_41 * lp_121[k];
    }

#pragma omp simd aligned(lp_2, lp_11, lp_17, lp_32, lp_38, lp_44, lp_65, lp_71, lp_77, lp_110, \
                         lp_116, lp_122 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_37[k] = f_37 * lp_2[k]
                  - f_16 * lp_11[k]
                  - f_38 * lp_17[k]
                  - f_39 * lp_32[k]
                  + f_40 * lp_38[k]
                  + f_41 * lp_44[k]
                  - f_16 * lp_65[k]
                  + f_40 * lp_71[k]
                  - f_42 * lp_77[k]
                  + f_37 * lp_110[k]
                  - f_38 * lp_116[k]
                  + f_41 * lp_122[k];
    }

#pragma omp simd aligned(lp_0, lp_9, lp_15, lp_30, lp_36, lp_42, lp_63, lp_69, lp_75, lp_108, \
                         lp_114, lp_120 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = f_37 * lp_0[k]
                  - f_16 * lp_9[k]
                  - f_38 * lp_15[k]
                  - f_39 * lp_30[k]
                  + f_40 * lp_36[k]
                  + f_41 * lp_42[k]
                  - f_16 * lp_63[k]
                  + f_40 * lp_69[k]
                  - f_42 * lp_75[k]
                  + f_37 * lp_108[k]
                  - f_38 * lp_114[k]
                  + f_41 * lp_120[k];
    }

#pragma omp simd aligned(lp_7, lp_8, lp_22, lp_23, lp_28, lp_29, lp_49, lp_50, lp_55, lp_56, \
                         lp_88, lp_89, lp_94, lp_95 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_39[k] = -f_14 * lp_7[k]
                  + f_12 * lp_22[k]
                  + f_15 * lp_28[k]
                  + f_10 * lp_49[k]
                  - f_13 * lp_55[k]
                  - f_10 * lp_88[k]
                  + f_11 * lp_94[k];

        g_40[k] = -f_14 * lp_8[k]
                  + f_12 * lp_23[k]
                  + f_15 * lp_29[k]
                  + f_10 * lp_50[k]
                  - f_13 * lp_56[k]
                  - f_10 * lp_89[k]
                  + f_11 * lp_95[k];
    }

#pragma omp simd aligned(lp_6, lp_21, lp_27, lp_48, lp_54, lp_87, \
                         lp_93 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_41[k] = -f_14 * lp_6[k]
                  + f_12 * lp_21[k]
                  + f_15 * lp_27[k]
                  + f_10 * lp_48[k]
                  - f_13 * lp_54[k]
                  - f_10 * lp_87[k]
                  + f_11 * lp_93[k];
    }

#pragma omp simd aligned(lp_1, lp_10, lp_16, lp_37, lp_64, lp_70, lp_109, \
                         lp_115 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = -f_43 * lp_1[k]
                  + f_7 * lp_10[k]
                  + f_7 * lp_16[k]
                  - f_44 * lp_37[k]
                  - f_7 * lp_64[k]
                  + f_44 * lp_70[k]
                  + f_43 * lp_109[k]
                  - f_7 * lp_115[k];
    }

#pragma omp simd aligned(lp_2, lp_11, lp_17, lp_38, lp_65, lp_71, lp_110, \
                         lp_116 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_43[k] = -f_43 * lp_2[k]
                  + f_7 * lp_11[k]
                  + f_7 * lp_17[k]
                  - f_44 * lp_38[k]
                  - f_7 * lp_65[k]
                  + f_44 * lp_71[k]
                  + f_43 * lp_110[k]
                  - f_7 * lp_116[k];
    }

#pragma omp simd aligned(lp_0, lp_7, lp_9, lp_15, lp_22, lp_36, lp_49, lp_63, lp_69, lp_88, \
                         lp_108, lp_114 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_44[k] = -f_43 * lp_0[k]
                  + f_7 * lp_9[k]
                  + f_7 * lp_15[k]
                  - f_44 * lp_36[k]
                  - f_7 * lp_63[k]
                  + f_44 * lp_69[k]
                  + f_43 * lp_108[k]
                  - f_7 * lp_114[k];

        g_45[k] = f_5 * lp_7[k]
                  - f_4 * lp_22[k]
                  + f_3 * lp_49[k]
                  - f_2 * lp_88[k];
    }

#pragma omp simd aligned(lp_1, lp_6, lp_8, lp_10, lp_21, lp_23, lp_31, lp_48, lp_50, lp_64, \
                         lp_87, lp_89, lp_109 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_46[k] = f_5 * lp_8[k]
                  - f_4 * lp_23[k]
                  + f_3 * lp_50[k]
                  - f_2 * lp_89[k];

        g_47[k] = f_5 * lp_6[k]
                  - f_4 * lp_21[k]
                  + f_3 * lp_48[k]
                  - f_2 * lp_87[k];

        g_48[k] = f_45 * lp_1[k]
                  - f_2 * lp_10[k]
                  + f_46 * lp_31[k]
                  - f_2 * lp_64[k]
                  + f_45 * lp_109[k];
    }

#pragma omp simd aligned(lp_0, lp_2, lp_9, lp_11, lp_30, lp_32, lp_63, lp_65, lp_108, \
                         lp_110 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_49[k] = f_45 * lp_2[k]
                  - f_2 * lp_11[k]
                  + f_46 * lp_32[k]
                  - f_2 * lp_65[k]
                  + f_45 * lp_110[k];

        g_50[k] = f_45 * lp_0[k]
                  - f_2 * lp_9[k]
                  + f_46 * lp_30[k]
                  - f_2 * lp_63[k]
                  + f_45 * lp_108[k];
    }
}

}  // namespace simdtrf
