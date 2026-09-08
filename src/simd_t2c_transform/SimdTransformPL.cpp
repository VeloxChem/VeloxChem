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


#include "SimdTransformPL.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_pl(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t pl,
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

    const auto *pl_0 = buffer.data(pl + 0);
    const auto *pl_1 = buffer.data(pl + 1);
    const auto *pl_2 = buffer.data(pl + 2);
    const auto *pl_3 = buffer.data(pl + 3);
    const auto *pl_4 = buffer.data(pl + 4);
    const auto *pl_5 = buffer.data(pl + 5);
    const auto *pl_6 = buffer.data(pl + 6);
    const auto *pl_7 = buffer.data(pl + 7);
    const auto *pl_8 = buffer.data(pl + 8);
    const auto *pl_9 = buffer.data(pl + 9);
    const auto *pl_10 = buffer.data(pl + 10);
    const auto *pl_11 = buffer.data(pl + 11);
    const auto *pl_12 = buffer.data(pl + 12);
    const auto *pl_13 = buffer.data(pl + 13);
    const auto *pl_14 = buffer.data(pl + 14);
    const auto *pl_15 = buffer.data(pl + 15);
    const auto *pl_16 = buffer.data(pl + 16);
    const auto *pl_17 = buffer.data(pl + 17);
    const auto *pl_18 = buffer.data(pl + 18);
    const auto *pl_19 = buffer.data(pl + 19);
    const auto *pl_20 = buffer.data(pl + 20);
    const auto *pl_21 = buffer.data(pl + 21);
    const auto *pl_22 = buffer.data(pl + 22);
    const auto *pl_23 = buffer.data(pl + 23);
    const auto *pl_24 = buffer.data(pl + 24);
    const auto *pl_25 = buffer.data(pl + 25);
    const auto *pl_26 = buffer.data(pl + 26);
    const auto *pl_27 = buffer.data(pl + 27);
    const auto *pl_28 = buffer.data(pl + 28);
    const auto *pl_29 = buffer.data(pl + 29);
    const auto *pl_30 = buffer.data(pl + 30);
    const auto *pl_31 = buffer.data(pl + 31);
    const auto *pl_32 = buffer.data(pl + 32);
    const auto *pl_33 = buffer.data(pl + 33);
    const auto *pl_34 = buffer.data(pl + 34);
    const auto *pl_35 = buffer.data(pl + 35);
    const auto *pl_36 = buffer.data(pl + 36);
    const auto *pl_37 = buffer.data(pl + 37);
    const auto *pl_38 = buffer.data(pl + 38);
    const auto *pl_39 = buffer.data(pl + 39);
    const auto *pl_40 = buffer.data(pl + 40);
    const auto *pl_41 = buffer.data(pl + 41);
    const auto *pl_42 = buffer.data(pl + 42);
    const auto *pl_43 = buffer.data(pl + 43);
    const auto *pl_44 = buffer.data(pl + 44);
    const auto *pl_45 = buffer.data(pl + 45);
    const auto *pl_46 = buffer.data(pl + 46);
    const auto *pl_47 = buffer.data(pl + 47);
    const auto *pl_48 = buffer.data(pl + 48);
    const auto *pl_49 = buffer.data(pl + 49);
    const auto *pl_50 = buffer.data(pl + 50);
    const auto *pl_51 = buffer.data(pl + 51);
    const auto *pl_52 = buffer.data(pl + 52);
    const auto *pl_53 = buffer.data(pl + 53);
    const auto *pl_54 = buffer.data(pl + 54);
    const auto *pl_55 = buffer.data(pl + 55);
    const auto *pl_56 = buffer.data(pl + 56);
    const auto *pl_57 = buffer.data(pl + 57);
    const auto *pl_58 = buffer.data(pl + 58);
    const auto *pl_59 = buffer.data(pl + 59);
    const auto *pl_60 = buffer.data(pl + 60);
    const auto *pl_61 = buffer.data(pl + 61);
    const auto *pl_62 = buffer.data(pl + 62);
    const auto *pl_63 = buffer.data(pl + 63);
    const auto *pl_64 = buffer.data(pl + 64);
    const auto *pl_65 = buffer.data(pl + 65);
    const auto *pl_66 = buffer.data(pl + 66);
    const auto *pl_67 = buffer.data(pl + 67);
    const auto *pl_68 = buffer.data(pl + 68);
    const auto *pl_69 = buffer.data(pl + 69);
    const auto *pl_70 = buffer.data(pl + 70);
    const auto *pl_71 = buffer.data(pl + 71);
    const auto *pl_72 = buffer.data(pl + 72);
    const auto *pl_73 = buffer.data(pl + 73);
    const auto *pl_74 = buffer.data(pl + 74);
    const auto *pl_75 = buffer.data(pl + 75);
    const auto *pl_76 = buffer.data(pl + 76);
    const auto *pl_77 = buffer.data(pl + 77);
    const auto *pl_78 = buffer.data(pl + 78);
    const auto *pl_79 = buffer.data(pl + 79);
    const auto *pl_80 = buffer.data(pl + 80);
    const auto *pl_81 = buffer.data(pl + 81);
    const auto *pl_82 = buffer.data(pl + 82);
    const auto *pl_83 = buffer.data(pl + 83);
    const auto *pl_84 = buffer.data(pl + 84);
    const auto *pl_85 = buffer.data(pl + 85);
    const auto *pl_86 = buffer.data(pl + 86);
    const auto *pl_87 = buffer.data(pl + 87);
    const auto *pl_88 = buffer.data(pl + 88);
    const auto *pl_89 = buffer.data(pl + 89);
    const auto *pl_90 = buffer.data(pl + 90);
    const auto *pl_91 = buffer.data(pl + 91);
    const auto *pl_92 = buffer.data(pl + 92);
    const auto *pl_93 = buffer.data(pl + 93);
    const auto *pl_94 = buffer.data(pl + 94);
    const auto *pl_95 = buffer.data(pl + 95);
    const auto *pl_96 = buffer.data(pl + 96);
    const auto *pl_97 = buffer.data(pl + 97);
    const auto *pl_98 = buffer.data(pl + 98);
    const auto *pl_99 = buffer.data(pl + 99);
    const auto *pl_100 = buffer.data(pl + 100);
    const auto *pl_101 = buffer.data(pl + 101);
    const auto *pl_102 = buffer.data(pl + 102);
    const auto *pl_103 = buffer.data(pl + 103);
    const auto *pl_104 = buffer.data(pl + 104);
    const auto *pl_105 = buffer.data(pl + 105);
    const auto *pl_106 = buffer.data(pl + 106);
    const auto *pl_107 = buffer.data(pl + 107);
    const auto *pl_108 = buffer.data(pl + 108);
    const auto *pl_109 = buffer.data(pl + 109);
    const auto *pl_110 = buffer.data(pl + 110);
    const auto *pl_111 = buffer.data(pl + 111);
    const auto *pl_112 = buffer.data(pl + 112);
    const auto *pl_113 = buffer.data(pl + 113);
    const auto *pl_114 = buffer.data(pl + 114);
    const auto *pl_115 = buffer.data(pl + 115);
    const auto *pl_116 = buffer.data(pl + 116);
    const auto *pl_117 = buffer.data(pl + 117);
    const auto *pl_118 = buffer.data(pl + 118);
    const auto *pl_119 = buffer.data(pl + 119);
    const auto *pl_120 = buffer.data(pl + 120);
    const auto *pl_121 = buffer.data(pl + 121);
    const auto *pl_122 = buffer.data(pl + 122);
    const auto *pl_123 = buffer.data(pl + 123);
    const auto *pl_124 = buffer.data(pl + 124);
    const auto *pl_125 = buffer.data(pl + 125);
    const auto *pl_126 = buffer.data(pl + 126);
    const auto *pl_127 = buffer.data(pl + 127);
    const auto *pl_128 = buffer.data(pl + 128);
    const auto *pl_129 = buffer.data(pl + 129);
    const auto *pl_130 = buffer.data(pl + 130);
    const auto *pl_131 = buffer.data(pl + 131);
    const auto *pl_132 = buffer.data(pl + 132);
    const auto *pl_133 = buffer.data(pl + 133);
    const auto *pl_134 = buffer.data(pl + 134);

#pragma omp simd aligned(pl_46, pl_49, pl_51, pl_53, pl_56, pl_60, pl_62, pl_67, pl_73, pl_75, \
                         pl_82 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * pl_46[k]
                 - f_1 * pl_51[k]
                 + f_1 * pl_60[k]
                 - f_0 * pl_73[k];

        g_1[k] = f_2 * pl_49[k]
                 - f_3 * pl_56[k]
                 + f_4 * pl_67[k]
                 - f_5 * pl_82[k];

        g_2[k] = -f_6 * pl_46[k]
                 + f_7 * pl_51[k]
                 + f_8 * pl_53[k]
                 + f_7 * pl_60[k]
                 - f_9 * pl_62[k]
                 - f_6 * pl_73[k]
                 + f_8 * pl_75[k];
    }

#pragma omp simd aligned(pl_49, pl_56, pl_58, pl_67, pl_69, pl_82, \
                         pl_84 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = -f_10 * pl_49[k]
                 + f_10 * pl_56[k]
                 + f_11 * pl_58[k]
                 + f_12 * pl_67[k]
                 - f_13 * pl_69[k]
                 - f_14 * pl_82[k]
                 + f_15 * pl_84[k];
    }

#pragma omp simd aligned(pl_46, pl_51, pl_53, pl_60, pl_64, pl_73, pl_75, \
                         pl_77 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = f_16 * pl_46[k]
                 + f_16 * pl_51[k]
                 - f_17 * pl_53[k]
                 - f_16 * pl_60[k]
                 + f_18 * pl_64[k]
                 - f_16 * pl_73[k]
                 + f_17 * pl_75[k]
                 - f_18 * pl_77[k];
    }

#pragma omp simd aligned(pl_49, pl_56, pl_58, pl_67, pl_69, pl_71, pl_82, pl_84, \
                         pl_86 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = f_19 * pl_49[k]
                 + f_20 * pl_56[k]
                 - f_21 * pl_58[k]
                 + f_22 * pl_67[k]
                 - f_23 * pl_69[k]
                 + f_24 * pl_71[k]
                 - f_22 * pl_82[k]
                 + f_25 * pl_84[k]
                 - f_26 * pl_86[k];
    }

#pragma omp simd aligned(pl_46, pl_51, pl_53, pl_60, pl_62, pl_64, pl_73, pl_75, pl_77, \
                         pl_79 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = -f_27 * pl_46[k]
                 - f_28 * pl_51[k]
                 + f_29 * pl_53[k]
                 - f_28 * pl_60[k]
                 + f_30 * pl_62[k]
                 - f_31 * pl_64[k]
                 - f_27 * pl_73[k]
                 + f_29 * pl_75[k]
                 - f_31 * pl_77[k]
                 + f_32 * pl_79[k];
    }

#pragma omp simd aligned(pl_49, pl_56, pl_58, pl_67, pl_69, pl_71, pl_82, pl_84, pl_86, \
                         pl_88 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = -3.28125 * pl_49[k]
                 - 9.84375 * pl_56[k]
                 + 26.25 * pl_58[k]
                 - 9.84375 * pl_67[k]
                 + 52.5 * pl_69[k]
                 - 31.5 * pl_71[k]
                 - 3.28125 * pl_82[k]
                 + 26.25 * pl_84[k]
                 - 31.5 * pl_86[k]
                 + 6.0 * pl_88[k];
    }

#pragma omp simd aligned(pl_45, pl_48, pl_50, pl_55, pl_57, pl_59, pl_66, pl_68, pl_70, pl_72, \
                         pl_81, pl_83, pl_85, pl_87, pl_89 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_8[k] = 0.2734375 * pl_45[k]
                 + 1.09375 * pl_48[k]
                 - 8.75 * pl_50[k]
                 + 1.640625 * pl_55[k]
                 - 26.25 * pl_57[k]
                 + 26.25 * pl_59[k]
                 + 1.09375 * pl_66[k]
                 - 26.25 * pl_68[k]
                 + 52.5 * pl_70[k]
                 - 14.0 * pl_72[k]
                 + 0.2734375 * pl_81[k]
                 - 8.75 * pl_83[k]
                 + 26.25 * pl_85[k]
                 - 14.0 * pl_87[k]
                 + pl_89[k];
    }

#pragma omp simd aligned(pl_47, pl_52, pl_54, pl_61, pl_63, pl_65, pl_74, pl_76, pl_78, \
                         pl_80 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_9[k] = -3.28125 * pl_47[k]
                 - 9.84375 * pl_52[k]
                 + 26.25 * pl_54[k]
                 - 9.84375 * pl_61[k]
                 + 52.5 * pl_63[k]
                 - 31.5 * pl_65[k]
                 - 3.28125 * pl_74[k]
                 + 26.25 * pl_76[k]
                 - 31.5 * pl_78[k]
                 + 6.0 * pl_80[k];
    }

#pragma omp simd aligned(pl_45, pl_48, pl_50, pl_57, pl_59, pl_66, pl_68, pl_72, pl_81, pl_83, \
                         pl_85, pl_87 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = -f_33 * pl_45[k]
                  - f_27 * pl_48[k]
                  + f_34 * pl_50[k]
                  + f_34 * pl_57[k]
                  - f_35 * pl_59[k]
                  + f_27 * pl_66[k]
                  - f_34 * pl_68[k]
                  + f_36 * pl_72[k]
                  + f_33 * pl_81[k]
                  - f_34 * pl_83[k]
                  + f_35 * pl_85[k]
                  - f_36 * pl_87[k];
    }

#pragma omp simd aligned(pl_47, pl_52, pl_54, pl_61, pl_63, pl_65, pl_74, pl_76, \
                         pl_78 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_11[k] = f_22 * pl_47[k]
                  - f_22 * pl_52[k]
                  - f_25 * pl_54[k]
                  - f_20 * pl_61[k]
                  + f_23 * pl_63[k]
                  + f_26 * pl_65[k]
                  - f_19 * pl_74[k]
                  + f_21 * pl_76[k]
                  - f_24 * pl_78[k];
    }

#pragma omp simd aligned(pl_45, pl_48, pl_50, pl_55, pl_57, pl_59, pl_66, pl_68, pl_70, pl_81, \
                         pl_83, pl_85 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_12[k] = f_37 * pl_45[k]
                  - f_16 * pl_48[k]
                  - f_38 * pl_50[k]
                  - f_39 * pl_55[k]
                  + f_40 * pl_57[k]
                  + f_41 * pl_59[k]
                  - f_16 * pl_66[k]
                  + f_40 * pl_68[k]
                  - f_42 * pl_70[k]
                  + f_37 * pl_81[k]
                  - f_38 * pl_83[k]
                  + f_41 * pl_85[k];
    }

#pragma omp simd aligned(pl_47, pl_52, pl_54, pl_61, pl_63, pl_74, \
                         pl_76 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = -f_14 * pl_47[k]
                  + f_12 * pl_52[k]
                  + f_15 * pl_54[k]
                  + f_10 * pl_61[k]
                  - f_13 * pl_63[k]
                  - f_10 * pl_74[k]
                  + f_11 * pl_76[k];
    }

#pragma omp simd aligned(pl_45, pl_47, pl_48, pl_50, pl_52, pl_55, pl_57, pl_61, pl_66, pl_68, \
                         pl_74, pl_81, pl_83 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = -f_43 * pl_45[k]
                  + f_7 * pl_48[k]
                  + f_7 * pl_50[k]
                  - f_44 * pl_57[k]
                  - f_7 * pl_66[k]
                  + f_44 * pl_68[k]
                  + f_43 * pl_81[k]
                  - f_7 * pl_83[k];

        g_15[k] = f_5 * pl_47[k]
                  - f_4 * pl_52[k]
                  + f_3 * pl_61[k]
                  - f_2 * pl_74[k];

        g_16[k] = f_45 * pl_45[k]
                  - f_2 * pl_48[k]
                  + f_46 * pl_55[k]
                  - f_2 * pl_66[k]
                  + f_45 * pl_81[k];
    }

#pragma omp simd aligned(pl_91, pl_94, pl_96, pl_98, pl_101, pl_105, pl_107, pl_112, pl_118, \
                         pl_120, pl_127 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = f_0 * pl_91[k]
                  - f_1 * pl_96[k]
                  + f_1 * pl_105[k]
                  - f_0 * pl_118[k];

        g_18[k] = f_2 * pl_94[k]
                  - f_3 * pl_101[k]
                  + f_4 * pl_112[k]
                  - f_5 * pl_127[k];

        g_19[k] = -f_6 * pl_91[k]
                  + f_7 * pl_96[k]
                  + f_8 * pl_98[k]
                  + f_7 * pl_105[k]
                  - f_9 * pl_107[k]
                  - f_6 * pl_118[k]
                  + f_8 * pl_120[k];
    }

#pragma omp simd aligned(pl_94, pl_101, pl_103, pl_112, pl_114, pl_127, \
                         pl_129 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = -f_10 * pl_94[k]
                  + f_10 * pl_101[k]
                  + f_11 * pl_103[k]
                  + f_12 * pl_112[k]
                  - f_13 * pl_114[k]
                  - f_14 * pl_127[k]
                  + f_15 * pl_129[k];
    }

#pragma omp simd aligned(pl_91, pl_96, pl_98, pl_105, pl_109, pl_118, pl_120, \
                         pl_122 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = f_16 * pl_91[k]
                  + f_16 * pl_96[k]
                  - f_17 * pl_98[k]
                  - f_16 * pl_105[k]
                  + f_18 * pl_109[k]
                  - f_16 * pl_118[k]
                  + f_17 * pl_120[k]
                  - f_18 * pl_122[k];
    }

#pragma omp simd aligned(pl_94, pl_101, pl_103, pl_112, pl_114, pl_116, pl_127, pl_129, \
                         pl_131 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = f_19 * pl_94[k]
                  + f_20 * pl_101[k]
                  - f_21 * pl_103[k]
                  + f_22 * pl_112[k]
                  - f_23 * pl_114[k]
                  + f_24 * pl_116[k]
                  - f_22 * pl_127[k]
                  + f_25 * pl_129[k]
                  - f_26 * pl_131[k];
    }

#pragma omp simd aligned(pl_91, pl_96, pl_98, pl_105, pl_107, pl_109, pl_118, pl_120, pl_122, \
                         pl_124 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = -f_27 * pl_91[k]
                  - f_28 * pl_96[k]
                  + f_29 * pl_98[k]
                  - f_28 * pl_105[k]
                  + f_30 * pl_107[k]
                  - f_31 * pl_109[k]
                  - f_27 * pl_118[k]
                  + f_29 * pl_120[k]
                  - f_31 * pl_122[k]
                  + f_32 * pl_124[k];
    }

#pragma omp simd aligned(pl_94, pl_101, pl_103, pl_112, pl_114, pl_116, pl_127, pl_129, \
                         pl_131, pl_133 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = -3.28125 * pl_94[k]
                  - 9.84375 * pl_101[k]
                  + 26.25 * pl_103[k]
                  - 9.84375 * pl_112[k]
                  + 52.5 * pl_114[k]
                  - 31.5 * pl_116[k]
                  - 3.28125 * pl_127[k]
                  + 26.25 * pl_129[k]
                  - 31.5 * pl_131[k]
                  + 6.0 * pl_133[k];
    }

#pragma omp simd aligned(pl_90, pl_93, pl_95, pl_100, pl_102, pl_104, pl_111, pl_113, pl_115, \
                         pl_117, pl_126, pl_128, pl_130, pl_132, \
                         pl_134 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = 0.2734375 * pl_90[k]
                  + 1.09375 * pl_93[k]
                  - 8.75 * pl_95[k]
                  + 1.640625 * pl_100[k]
                  - 26.25 * pl_102[k]
                  + 26.25 * pl_104[k]
                  + 1.09375 * pl_111[k]
                  - 26.25 * pl_113[k]
                  + 52.5 * pl_115[k]
                  - 14.0 * pl_117[k]
                  + 0.2734375 * pl_126[k]
                  - 8.75 * pl_128[k]
                  + 26.25 * pl_130[k]
                  - 14.0 * pl_132[k]
                  + pl_134[k];
    }

#pragma omp simd aligned(pl_92, pl_97, pl_99, pl_106, pl_108, pl_110, pl_119, pl_121, pl_123, \
                         pl_125 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = -3.28125 * pl_92[k]
                  - 9.84375 * pl_97[k]
                  + 26.25 * pl_99[k]
                  - 9.84375 * pl_106[k]
                  + 52.5 * pl_108[k]
                  - 31.5 * pl_110[k]
                  - 3.28125 * pl_119[k]
                  + 26.25 * pl_121[k]
                  - 31.5 * pl_123[k]
                  + 6.0 * pl_125[k];
    }

#pragma omp simd aligned(pl_90, pl_93, pl_95, pl_102, pl_104, pl_111, pl_113, pl_117, pl_126, \
                         pl_128, pl_130, pl_132 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = -f_33 * pl_90[k]
                  - f_27 * pl_93[k]
                  + f_34 * pl_95[k]
                  + f_34 * pl_102[k]
                  - f_35 * pl_104[k]
                  + f_27 * pl_111[k]
                  - f_34 * pl_113[k]
                  + f_36 * pl_117[k]
                  + f_33 * pl_126[k]
                  - f_34 * pl_128[k]
                  + f_35 * pl_130[k]
                  - f_36 * pl_132[k];
    }

#pragma omp simd aligned(pl_92, pl_97, pl_99, pl_106, pl_108, pl_110, pl_119, pl_121, \
                         pl_123 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = f_22 * pl_92[k]
                  - f_22 * pl_97[k]
                  - f_25 * pl_99[k]
                  - f_20 * pl_106[k]
                  + f_23 * pl_108[k]
                  + f_26 * pl_110[k]
                  - f_19 * pl_119[k]
                  + f_21 * pl_121[k]
                  - f_24 * pl_123[k];
    }

#pragma omp simd aligned(pl_90, pl_93, pl_95, pl_100, pl_102, pl_104, pl_111, pl_113, pl_115, \
                         pl_126, pl_128, pl_130 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = f_37 * pl_90[k]
                  - f_16 * pl_93[k]
                  - f_38 * pl_95[k]
                  - f_39 * pl_100[k]
                  + f_40 * pl_102[k]
                  + f_41 * pl_104[k]
                  - f_16 * pl_111[k]
                  + f_40 * pl_113[k]
                  - f_42 * pl_115[k]
                  + f_37 * pl_126[k]
                  - f_38 * pl_128[k]
                  + f_41 * pl_130[k];
    }

#pragma omp simd aligned(pl_92, pl_97, pl_99, pl_106, pl_108, pl_119, \
                         pl_121 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = -f_14 * pl_92[k]
                  + f_12 * pl_97[k]
                  + f_15 * pl_99[k]
                  + f_10 * pl_106[k]
                  - f_13 * pl_108[k]
                  - f_10 * pl_119[k]
                  + f_11 * pl_121[k];
    }

#pragma omp simd aligned(pl_90, pl_92, pl_93, pl_95, pl_97, pl_100, pl_102, pl_106, pl_111, \
                         pl_113, pl_119, pl_126, pl_128 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = -f_43 * pl_90[k]
                  + f_7 * pl_93[k]
                  + f_7 * pl_95[k]
                  - f_44 * pl_102[k]
                  - f_7 * pl_111[k]
                  + f_44 * pl_113[k]
                  + f_43 * pl_126[k]
                  - f_7 * pl_128[k];

        g_32[k] = f_5 * pl_92[k]
                  - f_4 * pl_97[k]
                  + f_3 * pl_106[k]
                  - f_2 * pl_119[k];

        g_33[k] = f_45 * pl_90[k]
                  - f_2 * pl_93[k]
                  + f_46 * pl_100[k]
                  - f_2 * pl_111[k]
                  + f_45 * pl_126[k];
    }

#pragma omp simd aligned(pl_1, pl_4, pl_6, pl_8, pl_11, pl_15, pl_17, pl_22, pl_28, pl_30, \
                         pl_37 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = f_0 * pl_1[k]
                  - f_1 * pl_6[k]
                  + f_1 * pl_15[k]
                  - f_0 * pl_28[k];

        g_35[k] = f_2 * pl_4[k]
                  - f_3 * pl_11[k]
                  + f_4 * pl_22[k]
                  - f_5 * pl_37[k];

        g_36[k] = -f_6 * pl_1[k]
                  + f_7 * pl_6[k]
                  + f_8 * pl_8[k]
                  + f_7 * pl_15[k]
                  - f_9 * pl_17[k]
                  - f_6 * pl_28[k]
                  + f_8 * pl_30[k];
    }

#pragma omp simd aligned(pl_4, pl_11, pl_13, pl_22, pl_24, pl_37, \
                         pl_39 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_37[k] = -f_10 * pl_4[k]
                  + f_10 * pl_11[k]
                  + f_11 * pl_13[k]
                  + f_12 * pl_22[k]
                  - f_13 * pl_24[k]
                  - f_14 * pl_37[k]
                  + f_15 * pl_39[k];
    }

#pragma omp simd aligned(pl_1, pl_6, pl_8, pl_15, pl_19, pl_28, pl_30, \
                         pl_32 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = f_16 * pl_1[k]
                  + f_16 * pl_6[k]
                  - f_17 * pl_8[k]
                  - f_16 * pl_15[k]
                  + f_18 * pl_19[k]
                  - f_16 * pl_28[k]
                  + f_17 * pl_30[k]
                  - f_18 * pl_32[k];
    }

#pragma omp simd aligned(pl_4, pl_11, pl_13, pl_22, pl_24, pl_26, pl_37, pl_39, \
                         pl_41 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_39[k] = f_19 * pl_4[k]
                  + f_20 * pl_11[k]
                  - f_21 * pl_13[k]
                  + f_22 * pl_22[k]
                  - f_23 * pl_24[k]
                  + f_24 * pl_26[k]
                  - f_22 * pl_37[k]
                  + f_25 * pl_39[k]
                  - f_26 * pl_41[k];
    }

#pragma omp simd aligned(pl_1, pl_6, pl_8, pl_15, pl_17, pl_19, pl_28, pl_30, pl_32, \
                         pl_34 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_40[k] = -f_27 * pl_1[k]
                  - f_28 * pl_6[k]
                  + f_29 * pl_8[k]
                  - f_28 * pl_15[k]
                  + f_30 * pl_17[k]
                  - f_31 * pl_19[k]
                  - f_27 * pl_28[k]
                  + f_29 * pl_30[k]
                  - f_31 * pl_32[k]
                  + f_32 * pl_34[k];
    }

#pragma omp simd aligned(pl_4, pl_11, pl_13, pl_22, pl_24, pl_26, pl_37, pl_39, pl_41, \
                         pl_43 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_41[k] = -3.28125 * pl_4[k]
                  - 9.84375 * pl_11[k]
                  + 26.25 * pl_13[k]
                  - 9.84375 * pl_22[k]
                  + 52.5 * pl_24[k]
                  - 31.5 * pl_26[k]
                  - 3.28125 * pl_37[k]
                  + 26.25 * pl_39[k]
                  - 31.5 * pl_41[k]
                  + 6.0 * pl_43[k];
    }

#pragma omp simd aligned(pl_0, pl_3, pl_5, pl_10, pl_12, pl_14, pl_21, pl_23, pl_25, pl_27, \
                         pl_36, pl_38, pl_40, pl_42, pl_44 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = 0.2734375 * pl_0[k]
                  + 1.09375 * pl_3[k]
                  - 8.75 * pl_5[k]
                  + 1.640625 * pl_10[k]
                  - 26.25 * pl_12[k]
                  + 26.25 * pl_14[k]
                  + 1.09375 * pl_21[k]
                  - 26.25 * pl_23[k]
                  + 52.5 * pl_25[k]
                  - 14.0 * pl_27[k]
                  + 0.2734375 * pl_36[k]
                  - 8.75 * pl_38[k]
                  + 26.25 * pl_40[k]
                  - 14.0 * pl_42[k]
                  + pl_44[k];
    }

#pragma omp simd aligned(pl_2, pl_7, pl_9, pl_16, pl_18, pl_20, pl_29, pl_31, pl_33, \
                         pl_35 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_43[k] = -3.28125 * pl_2[k]
                  - 9.84375 * pl_7[k]
                  + 26.25 * pl_9[k]
                  - 9.84375 * pl_16[k]
                  + 52.5 * pl_18[k]
                  - 31.5 * pl_20[k]
                  - 3.28125 * pl_29[k]
                  + 26.25 * pl_31[k]
                  - 31.5 * pl_33[k]
                  + 6.0 * pl_35[k];
    }

#pragma omp simd aligned(pl_0, pl_3, pl_5, pl_12, pl_14, pl_21, pl_23, pl_27, pl_36, pl_38, \
                         pl_40, pl_42 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_44[k] = -f_33 * pl_0[k]
                  - f_27 * pl_3[k]
                  + f_34 * pl_5[k]
                  + f_34 * pl_12[k]
                  - f_35 * pl_14[k]
                  + f_27 * pl_21[k]
                  - f_34 * pl_23[k]
                  + f_36 * pl_27[k]
                  + f_33 * pl_36[k]
                  - f_34 * pl_38[k]
                  + f_35 * pl_40[k]
                  - f_36 * pl_42[k];
    }

#pragma omp simd aligned(pl_2, pl_7, pl_9, pl_16, pl_18, pl_20, pl_29, pl_31, \
                         pl_33 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_45[k] = f_22 * pl_2[k]
                  - f_22 * pl_7[k]
                  - f_25 * pl_9[k]
                  - f_20 * pl_16[k]
                  + f_23 * pl_18[k]
                  + f_26 * pl_20[k]
                  - f_19 * pl_29[k]
                  + f_21 * pl_31[k]
                  - f_24 * pl_33[k];
    }

#pragma omp simd aligned(pl_0, pl_3, pl_5, pl_10, pl_12, pl_14, pl_21, pl_23, pl_25, pl_36, \
                         pl_38, pl_40 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_46[k] = f_37 * pl_0[k]
                  - f_16 * pl_3[k]
                  - f_38 * pl_5[k]
                  - f_39 * pl_10[k]
                  + f_40 * pl_12[k]
                  + f_41 * pl_14[k]
                  - f_16 * pl_21[k]
                  + f_40 * pl_23[k]
                  - f_42 * pl_25[k]
                  + f_37 * pl_36[k]
                  - f_38 * pl_38[k]
                  + f_41 * pl_40[k];
    }

#pragma omp simd aligned(pl_2, pl_7, pl_9, pl_16, pl_18, pl_29, pl_31 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_47[k] = -f_14 * pl_2[k]
                  + f_12 * pl_7[k]
                  + f_15 * pl_9[k]
                  + f_10 * pl_16[k]
                  - f_13 * pl_18[k]
                  - f_10 * pl_29[k]
                  + f_11 * pl_31[k];
    }

#pragma omp simd aligned(pl_0, pl_2, pl_3, pl_5, pl_7, pl_10, pl_12, pl_16, pl_21, pl_23, \
                         pl_29, pl_36, pl_38 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_48[k] = -f_43 * pl_0[k]
                  + f_7 * pl_3[k]
                  + f_7 * pl_5[k]
                  - f_44 * pl_12[k]
                  - f_7 * pl_21[k]
                  + f_44 * pl_23[k]
                  + f_43 * pl_36[k]
                  - f_7 * pl_38[k];

        g_49[k] = f_5 * pl_2[k]
                  - f_4 * pl_7[k]
                  + f_3 * pl_16[k]
                  - f_2 * pl_29[k];

        g_50[k] = f_45 * pl_0[k]
                  - f_2 * pl_3[k]
                  + f_46 * pl_10[k]
                  - f_2 * pl_21[k]
                  + f_45 * pl_36[k];
    }
}

}  // namespace simdtrf
