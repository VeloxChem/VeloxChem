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


#include "SimdTransformHD.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_hd(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t hd,
             const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 0.9375 * std::sqrt(42.0);
    const auto f_1 = 1.875 * std::sqrt(42.0);
    const auto f_2 = 0.1875 * std::sqrt(42.0);
    const auto f_3 = 0.46875 * std::sqrt(14.0);
    const auto f_4 = 0.9375 * std::sqrt(14.0);
    const auto f_5 = 1.875 * std::sqrt(14.0);
    const auto f_6 = 0.09375 * std::sqrt(14.0);
    const auto f_7 = 0.1875 * std::sqrt(14.0);
    const auto f_8 = 0.46875 * std::sqrt(42.0);
    const auto f_9 = 0.09375 * std::sqrt(42.0);
    const auto f_10 = 1.5 * std::sqrt(105.0);
    const auto f_11 = 0.75 * std::sqrt(35.0);
    const auto f_12 = 1.5 * std::sqrt(35.0);
    const auto f_13 = 0.75 * std::sqrt(105.0);
    const auto f_14 = 0.1875 * std::sqrt(210.0);
    const auto f_15 = 0.125 * std::sqrt(210.0);
    const auto f_16 = 1.5 * std::sqrt(210.0);
    const auto f_17 = 0.0625 * std::sqrt(210.0);
    const auto f_18 = 0.5 * std::sqrt(210.0);
    const auto f_19 = 0.09375 * std::sqrt(70.0);
    const auto f_20 = 0.1875 * std::sqrt(70.0);
    const auto f_21 = 0.0625 * std::sqrt(70.0);
    const auto f_22 = 0.125 * std::sqrt(70.0);
    const auto f_23 = 0.75 * std::sqrt(70.0);
    const auto f_24 = 1.5 * std::sqrt(70.0);
    const auto f_25 = 0.03125 * std::sqrt(70.0);
    const auto f_26 = 0.25 * std::sqrt(70.0);
    const auto f_27 = 0.5 * std::sqrt(70.0);
    const auto f_28 = 0.09375 * std::sqrt(210.0);
    const auto f_29 = 0.75 * std::sqrt(210.0);
    const auto f_30 = 0.03125 * std::sqrt(210.0);
    const auto f_31 = 0.25 * std::sqrt(210.0);
    const auto f_32 = 3.0 * std::sqrt(35.0);
    const auto f_33 = 0.25 * std::sqrt(105.0);
    const auto f_34 = 0.5 * std::sqrt(105.0);
    const auto f_35 = std::sqrt(105.0);
    const auto f_36 = 0.375 * std::sqrt(5.0);
    const auto f_37 = 0.75 * std::sqrt(5.0);
    const auto f_38 = 4.5 * std::sqrt(5.0);
    const auto f_39 = 3.0 * std::sqrt(5.0);
    const auto f_40 = 0.0625 * std::sqrt(15.0);
    const auto f_41 = 0.125 * std::sqrt(15.0);
    const auto f_42 = 0.25 * std::sqrt(15.0);
    const auto f_43 = 0.75 * std::sqrt(15.0);
    const auto f_44 = 1.5 * std::sqrt(15.0);
    const auto f_45 = 0.5 * std::sqrt(15.0);
    const auto f_46 = std::sqrt(15.0);
    const auto f_47 = 0.1875 * std::sqrt(5.0);
    const auto f_48 = 2.25 * std::sqrt(5.0);
    const auto f_49 = 1.5 * std::sqrt(5.0);
    const auto f_50 = 1.875 * std::sqrt(3.0);
    const auto f_51 = 3.75 * std::sqrt(3.0);
    const auto f_52 = 5.0 * std::sqrt(3.0);
    const auto f_53 = std::sqrt(3.0);
    const auto f_54 = 0.9375 * std::sqrt(3.0);
    const auto f_55 = 2.5 * std::sqrt(3.0);
    const auto f_56 = 0.5 * std::sqrt(3.0);
    const auto f_57 = 0.125 * std::sqrt(105.0);
    const auto f_58 = 0.375 * std::sqrt(35.0);
    const auto f_59 = 0.375 * std::sqrt(105.0);
    const auto f_60 = 2.25 * std::sqrt(105.0);
    const auto f_61 = 0.1875 * std::sqrt(35.0);
    const auto f_62 = 1.125 * std::sqrt(35.0);
    const auto f_63 = 2.25 * std::sqrt(35.0);
    const auto f_64 = 0.1875 * std::sqrt(105.0);
    const auto f_65 = 1.125 * std::sqrt(105.0);

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
    auto *g_51 = values + 51 * nvalues;
    auto *g_52 = values + 52 * nvalues;
    auto *g_53 = values + 53 * nvalues;
    auto *g_54 = values + 54 * nvalues;

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_1 = buffer.data(hd + 1);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_3 = buffer.data(hd + 3);
    const auto *hd_4 = buffer.data(hd + 4);
    const auto *hd_5 = buffer.data(hd + 5);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_7 = buffer.data(hd + 7);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_9 = buffer.data(hd + 9);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_13 = buffer.data(hd + 13);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_15 = buffer.data(hd + 15);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_22 = buffer.data(hd + 22);
    const auto *hd_23 = buffer.data(hd + 23);
    const auto *hd_24 = buffer.data(hd + 24);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_26 = buffer.data(hd + 26);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_30 = buffer.data(hd + 30);
    const auto *hd_31 = buffer.data(hd + 31);
    const auto *hd_32 = buffer.data(hd + 32);
    const auto *hd_33 = buffer.data(hd + 33);
    const auto *hd_34 = buffer.data(hd + 34);
    const auto *hd_35 = buffer.data(hd + 35);
    const auto *hd_36 = buffer.data(hd + 36);
    const auto *hd_37 = buffer.data(hd + 37);
    const auto *hd_38 = buffer.data(hd + 38);
    const auto *hd_39 = buffer.data(hd + 39);
    const auto *hd_40 = buffer.data(hd + 40);
    const auto *hd_41 = buffer.data(hd + 41);
    const auto *hd_42 = buffer.data(hd + 42);
    const auto *hd_43 = buffer.data(hd + 43);
    const auto *hd_44 = buffer.data(hd + 44);
    const auto *hd_45 = buffer.data(hd + 45);
    const auto *hd_46 = buffer.data(hd + 46);
    const auto *hd_47 = buffer.data(hd + 47);
    const auto *hd_48 = buffer.data(hd + 48);
    const auto *hd_49 = buffer.data(hd + 49);
    const auto *hd_50 = buffer.data(hd + 50);
    const auto *hd_51 = buffer.data(hd + 51);
    const auto *hd_52 = buffer.data(hd + 52);
    const auto *hd_53 = buffer.data(hd + 53);
    const auto *hd_54 = buffer.data(hd + 54);
    const auto *hd_55 = buffer.data(hd + 55);
    const auto *hd_56 = buffer.data(hd + 56);
    const auto *hd_57 = buffer.data(hd + 57);
    const auto *hd_58 = buffer.data(hd + 58);
    const auto *hd_59 = buffer.data(hd + 59);
    const auto *hd_60 = buffer.data(hd + 60);
    const auto *hd_61 = buffer.data(hd + 61);
    const auto *hd_62 = buffer.data(hd + 62);
    const auto *hd_63 = buffer.data(hd + 63);
    const auto *hd_64 = buffer.data(hd + 64);
    const auto *hd_65 = buffer.data(hd + 65);
    const auto *hd_66 = buffer.data(hd + 66);
    const auto *hd_67 = buffer.data(hd + 67);
    const auto *hd_68 = buffer.data(hd + 68);
    const auto *hd_69 = buffer.data(hd + 69);
    const auto *hd_70 = buffer.data(hd + 70);
    const auto *hd_71 = buffer.data(hd + 71);
    const auto *hd_72 = buffer.data(hd + 72);
    const auto *hd_73 = buffer.data(hd + 73);
    const auto *hd_74 = buffer.data(hd + 74);
    const auto *hd_75 = buffer.data(hd + 75);
    const auto *hd_76 = buffer.data(hd + 76);
    const auto *hd_77 = buffer.data(hd + 77);
    const auto *hd_78 = buffer.data(hd + 78);
    const auto *hd_79 = buffer.data(hd + 79);
    const auto *hd_80 = buffer.data(hd + 80);
    const auto *hd_81 = buffer.data(hd + 81);
    const auto *hd_82 = buffer.data(hd + 82);
    const auto *hd_83 = buffer.data(hd + 83);
    const auto *hd_84 = buffer.data(hd + 84);
    const auto *hd_85 = buffer.data(hd + 85);
    const auto *hd_86 = buffer.data(hd + 86);
    const auto *hd_87 = buffer.data(hd + 87);
    const auto *hd_88 = buffer.data(hd + 88);
    const auto *hd_89 = buffer.data(hd + 89);
    const auto *hd_90 = buffer.data(hd + 90);
    const auto *hd_91 = buffer.data(hd + 91);
    const auto *hd_92 = buffer.data(hd + 92);
    const auto *hd_93 = buffer.data(hd + 93);
    const auto *hd_94 = buffer.data(hd + 94);
    const auto *hd_95 = buffer.data(hd + 95);
    const auto *hd_96 = buffer.data(hd + 96);
    const auto *hd_97 = buffer.data(hd + 97);
    const auto *hd_98 = buffer.data(hd + 98);
    const auto *hd_99 = buffer.data(hd + 99);
    const auto *hd_100 = buffer.data(hd + 100);
    const auto *hd_101 = buffer.data(hd + 101);
    const auto *hd_102 = buffer.data(hd + 102);
    const auto *hd_103 = buffer.data(hd + 103);
    const auto *hd_104 = buffer.data(hd + 104);
    const auto *hd_105 = buffer.data(hd + 105);
    const auto *hd_106 = buffer.data(hd + 106);
    const auto *hd_107 = buffer.data(hd + 107);
    const auto *hd_108 = buffer.data(hd + 108);
    const auto *hd_109 = buffer.data(hd + 109);
    const auto *hd_110 = buffer.data(hd + 110);
    const auto *hd_111 = buffer.data(hd + 111);
    const auto *hd_112 = buffer.data(hd + 112);
    const auto *hd_113 = buffer.data(hd + 113);
    const auto *hd_114 = buffer.data(hd + 114);
    const auto *hd_115 = buffer.data(hd + 115);
    const auto *hd_116 = buffer.data(hd + 116);
    const auto *hd_117 = buffer.data(hd + 117);
    const auto *hd_118 = buffer.data(hd + 118);
    const auto *hd_119 = buffer.data(hd + 119);
    const auto *hd_120 = buffer.data(hd + 120);
    const auto *hd_121 = buffer.data(hd + 121);
    const auto *hd_122 = buffer.data(hd + 122);
    const auto *hd_123 = buffer.data(hd + 123);
    const auto *hd_124 = buffer.data(hd + 124);
    const auto *hd_125 = buffer.data(hd + 125);

#pragma omp simd aligned(hd_7, hd_10, hd_37, hd_40, hd_91, hd_94 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * hd_7[k]
                 - f_1 * hd_37[k]
                 + f_2 * hd_91[k];

        g_1[k] = f_0 * hd_10[k]
                 - f_1 * hd_40[k]
                 + f_2 * hd_94[k];
    }

#pragma omp simd aligned(hd_6, hd_8, hd_9, hd_11, hd_36, hd_38, hd_39, hd_41, hd_90, hd_92, \
                         hd_93, hd_95 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_2[k] = -f_3 * hd_6[k]
                 - f_3 * hd_9[k]
                 + f_4 * hd_11[k]
                 + f_4 * hd_36[k]
                 + f_4 * hd_39[k]
                 - f_5 * hd_41[k]
                 - f_6 * hd_90[k]
                 - f_6 * hd_93[k]
                 + f_7 * hd_95[k];

        g_3[k] = f_0 * hd_8[k]
                 - f_1 * hd_38[k]
                 + f_2 * hd_92[k];

        g_4[k] = f_8 * hd_6[k]
                 - f_8 * hd_9[k]
                 - f_0 * hd_36[k]
                 + f_0 * hd_39[k]
                 + f_9 * hd_90[k]
                 - f_9 * hd_93[k];
    }

#pragma omp simd aligned(hd_24, hd_25, hd_26, hd_27, hd_28, hd_29, hd_66, hd_67, hd_68, hd_69, \
                         hd_70, hd_71 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = f_10 * hd_25[k]
                 - f_10 * hd_67[k];

        g_6[k] = f_10 * hd_28[k]
                 - f_10 * hd_70[k];

        g_7[k] = -f_11 * hd_24[k]
                 - f_11 * hd_27[k]
                 + f_12 * hd_29[k]
                 + f_11 * hd_66[k]
                 + f_11 * hd_69[k]
                 - f_12 * hd_71[k];

        g_8[k] = f_10 * hd_26[k]
                 - f_10 * hd_68[k];
    }

#pragma omp simd aligned(hd_7, hd_24, hd_27, hd_37, hd_49, hd_66, hd_69, hd_91, \
                         hd_103 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_9[k] = f_13 * hd_24[k]
                 - f_13 * hd_27[k]
                 - f_13 * hd_66[k]
                 + f_13 * hd_69[k];

        g_10[k] = -f_14 * hd_7[k]
                  - f_15 * hd_37[k]
                  + f_16 * hd_49[k]
                  + f_17 * hd_91[k]
                  - f_18 * hd_103[k];
    }

#pragma omp simd aligned(hd_10, hd_40, hd_52, hd_94, hd_106 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_11[k] = -f_14 * hd_10[k]
                  - f_15 * hd_40[k]
                  + f_16 * hd_52[k]
                  + f_17 * hd_94[k]
                  - f_18 * hd_106[k];
    }

#pragma omp simd aligned(hd_6, hd_9, hd_11, hd_36, hd_39, hd_41, hd_48, hd_51, hd_53, hd_90, \
                         hd_93, hd_95, hd_102, hd_105, hd_107 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_12[k] = f_19 * hd_6[k]
                  + f_19 * hd_9[k]
                  - f_20 * hd_11[k]
                  + f_21 * hd_36[k]
                  + f_21 * hd_39[k]
                  - f_22 * hd_41[k]
                  - f_23 * hd_48[k]
                  - f_23 * hd_51[k]
                  + f_24 * hd_53[k]
                  - f_25 * hd_90[k]
                  - f_25 * hd_93[k]
                  + f_21 * hd_95[k]
                  + f_26 * hd_102[k]
                  + f_26 * hd_105[k]
                  - f_27 * hd_107[k];
    }

#pragma omp simd aligned(hd_8, hd_38, hd_50, hd_92, hd_104 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = -f_14 * hd_8[k]
                  - f_15 * hd_38[k]
                  + f_16 * hd_50[k]
                  + f_17 * hd_92[k]
                  - f_18 * hd_104[k];
    }

#pragma omp simd aligned(hd_6, hd_9, hd_25, hd_36, hd_39, hd_48, hd_51, hd_67, hd_79, hd_90, \
                         hd_93, hd_102, hd_105 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = -f_28 * hd_6[k]
                  + f_28 * hd_9[k]
                  - f_17 * hd_36[k]
                  + f_17 * hd_39[k]
                  + f_29 * hd_48[k]
                  - f_29 * hd_51[k]
                  + f_30 * hd_90[k]
                  - f_30 * hd_93[k]
                  - f_31 * hd_102[k]
                  + f_31 * hd_105[k];

        g_15[k] = -f_12 * hd_25[k]
                  - f_12 * hd_67[k]
                  + f_32 * hd_79[k];
    }

#pragma omp simd aligned(hd_24, hd_27, hd_28, hd_29, hd_66, hd_69, hd_70, hd_71, hd_78, hd_81, \
                         hd_82, hd_83 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_16[k] = -f_12 * hd_28[k]
                  - f_12 * hd_70[k]
                  + f_32 * hd_82[k];

        g_17[k] = f_33 * hd_24[k]
                  + f_33 * hd_27[k]
                  - f_34 * hd_29[k]
                  + f_33 * hd_66[k]
                  + f_33 * hd_69[k]
                  - f_34 * hd_71[k]
                  - f_34 * hd_78[k]
                  - f_34 * hd_81[k]
                  + f_35 * hd_83[k];
    }

#pragma omp simd aligned(hd_24, hd_26, hd_27, hd_66, hd_68, hd_69, hd_78, hd_80, \
                         hd_81 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = -f_12 * hd_26[k]
                  - f_12 * hd_68[k]
                  + f_32 * hd_80[k];

        g_19[k] = -f_11 * hd_24[k]
                  + f_11 * hd_27[k]
                  - f_11 * hd_66[k]
                  + f_11 * hd_69[k]
                  + f_12 * hd_78[k]
                  - f_12 * hd_81[k];
    }

#pragma omp simd aligned(hd_7, hd_10, hd_37, hd_40, hd_49, hd_52, hd_91, hd_94, hd_103, \
                         hd_106, hd_115, hd_118 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = f_36 * hd_7[k]
                  + f_37 * hd_37[k]
                  - f_38 * hd_49[k]
                  + f_36 * hd_91[k]
                  - f_38 * hd_103[k]
                  + f_39 * hd_115[k];

        g_21[k] = f_36 * hd_10[k]
                  + f_37 * hd_40[k]
                  - f_38 * hd_52[k]
                  + f_36 * hd_94[k]
                  - f_38 * hd_106[k]
                  + f_39 * hd_118[k];
    }

#pragma omp simd aligned(hd_6, hd_9, hd_11, hd_36, hd_39, hd_41, hd_48, hd_51, hd_53, hd_90, \
                         hd_93, hd_95, hd_102, hd_105, hd_107, hd_114, hd_117, \
                         hd_119 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = -f_40 * hd_6[k]
                  - f_40 * hd_9[k]
                  + f_41 * hd_11[k]
                  - f_41 * hd_36[k]
                  - f_41 * hd_39[k]
                  + f_42 * hd_41[k]
                  + f_43 * hd_48[k]
                  + f_43 * hd_51[k]
                  - f_44 * hd_53[k]
                  - f_40 * hd_90[k]
                  - f_40 * hd_93[k]
                  + f_41 * hd_95[k]
                  + f_43 * hd_102[k]
                  + f_43 * hd_105[k]
                  - f_44 * hd_107[k]
                  - f_45 * hd_114[k]
                  - f_45 * hd_117[k]
                  + f_46 * hd_119[k];
    }

#pragma omp simd aligned(hd_8, hd_38, hd_50, hd_92, hd_104, hd_116 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = f_36 * hd_8[k]
                  + f_37 * hd_38[k]
                  - f_38 * hd_50[k]
                  + f_36 * hd_92[k]
                  - f_38 * hd_104[k]
                  + f_39 * hd_116[k];
    }

#pragma omp simd aligned(hd_6, hd_9, hd_36, hd_39, hd_48, hd_51, hd_90, hd_93, hd_102, hd_105, \
                         hd_114, hd_117 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = f_47 * hd_6[k]
                  - f_47 * hd_9[k]
                  + f_36 * hd_36[k]
                  - f_36 * hd_39[k]
                  - f_48 * hd_48[k]
                  + f_48 * hd_51[k]
                  + f_47 * hd_90[k]
                  - f_47 * hd_93[k]
                  - f_48 * hd_102[k]
                  + f_48 * hd_105[k]
                  + f_49 * hd_114[k]
                  - f_49 * hd_117[k];
    }

#pragma omp simd aligned(hd_13, hd_16, hd_43, hd_46, hd_55, hd_58, hd_97, hd_100, hd_109, \
                         hd_112, hd_121, hd_124 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = f_50 * hd_13[k]
                  + f_51 * hd_43[k]
                  - f_52 * hd_55[k]
                  + f_50 * hd_97[k]
                  - f_52 * hd_109[k]
                  + f_53 * hd_121[k];

        g_26[k] = f_50 * hd_16[k]
                  + f_51 * hd_46[k]
                  - f_52 * hd_58[k]
                  + f_50 * hd_100[k]
                  - f_52 * hd_112[k]
                  + f_53 * hd_124[k];
    }

#pragma omp simd aligned(hd_12, hd_15, hd_17, hd_42, hd_45, hd_47, hd_54, hd_57, hd_59, hd_96, \
                         hd_99, hd_101, hd_108, hd_111, hd_113, hd_120, hd_123, \
                         hd_125 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = -0.9375 * hd_12[k]
                  - 0.9375 * hd_15[k]
                  + 1.875 * hd_17[k]
                  - 1.875 * hd_42[k]
                  - 1.875 * hd_45[k]
                  + 3.75 * hd_47[k]
                  + 2.5 * hd_54[k]
                  + 2.5 * hd_57[k]
                  - 5.0 * hd_59[k]
                  - 0.9375 * hd_96[k]
                  - 0.9375 * hd_99[k]
                  + 1.875 * hd_101[k]
                  + 2.5 * hd_108[k]
                  + 2.5 * hd_111[k]
                  - 5.0 * hd_113[k]
                  - 0.5 * hd_120[k]
                  - 0.5 * hd_123[k]
                  + hd_125[k];
    }

#pragma omp simd aligned(hd_14, hd_44, hd_56, hd_98, hd_110, hd_122 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = f_50 * hd_14[k]
                  + f_51 * hd_44[k]
                  - f_52 * hd_56[k]
                  + f_50 * hd_98[k]
                  - f_52 * hd_110[k]
                  + f_53 * hd_122[k];
    }

#pragma omp simd aligned(hd_12, hd_15, hd_42, hd_45, hd_54, hd_57, hd_96, hd_99, hd_108, \
                         hd_111, hd_120, hd_123 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = f_54 * hd_12[k]
                  - f_54 * hd_15[k]
                  + f_50 * hd_42[k]
                  - f_50 * hd_45[k]
                  - f_55 * hd_54[k]
                  + f_55 * hd_57[k]
                  + f_54 * hd_96[k]
                  - f_54 * hd_99[k]
                  - f_55 * hd_108[k]
                  + f_55 * hd_111[k]
                  + f_56 * hd_120[k]
                  - f_56 * hd_123[k];
    }

#pragma omp simd aligned(hd_1, hd_4, hd_19, hd_22, hd_31, hd_34, hd_61, hd_64, hd_73, hd_76, \
                         hd_85, hd_88 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = f_36 * hd_1[k]
                  + f_37 * hd_19[k]
                  - f_38 * hd_31[k]
                  + f_36 * hd_61[k]
                  - f_38 * hd_73[k]
                  + f_39 * hd_85[k];

        g_31[k] = f_36 * hd_4[k]
                  + f_37 * hd_22[k]
                  - f_38 * hd_34[k]
                  + f_36 * hd_64[k]
                  - f_38 * hd_76[k]
                  + f_39 * hd_88[k];
    }

#pragma omp simd aligned(hd_0, hd_3, hd_5, hd_18, hd_21, hd_23, hd_30, hd_33, hd_35, hd_60, \
                         hd_63, hd_65, hd_72, hd_75, hd_77, hd_84, hd_87, \
                         hd_89 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = -f_40 * hd_0[k]
                  - f_40 * hd_3[k]
                  + f_41 * hd_5[k]
                  - f_41 * hd_18[k]
                  - f_41 * hd_21[k]
                  + f_42 * hd_23[k]
                  + f_43 * hd_30[k]
                  + f_43 * hd_33[k]
                  - f_44 * hd_35[k]
                  - f_40 * hd_60[k]
                  - f_40 * hd_63[k]
                  + f_41 * hd_65[k]
                  + f_43 * hd_72[k]
                  + f_43 * hd_75[k]
                  - f_44 * hd_77[k]
                  - f_45 * hd_84[k]
                  - f_45 * hd_87[k]
                  + f_46 * hd_89[k];
    }

#pragma omp simd aligned(hd_2, hd_20, hd_32, hd_62, hd_74, hd_86 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = f_36 * hd_2[k]
                  + f_37 * hd_20[k]
                  - f_38 * hd_32[k]
                  + f_36 * hd_62[k]
                  - f_38 * hd_74[k]
                  + f_39 * hd_86[k];
    }

#pragma omp simd aligned(hd_0, hd_3, hd_18, hd_21, hd_30, hd_33, hd_60, hd_63, hd_72, hd_75, \
                         hd_84, hd_87 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = f_47 * hd_0[k]
                  - f_47 * hd_3[k]
                  + f_36 * hd_18[k]
                  - f_36 * hd_21[k]
                  - f_48 * hd_30[k]
                  + f_48 * hd_33[k]
                  + f_47 * hd_60[k]
                  - f_47 * hd_63[k]
                  - f_48 * hd_72[k]
                  + f_48 * hd_75[k]
                  + f_49 * hd_84[k]
                  - f_49 * hd_87[k];
    }

#pragma omp simd aligned(hd_13, hd_16, hd_55, hd_58, hd_97, hd_100, hd_109, \
                         hd_112 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = -f_11 * hd_13[k]
                  + f_12 * hd_55[k]
                  + f_11 * hd_97[k]
                  - f_12 * hd_109[k];

        g_36[k] = -f_11 * hd_16[k]
                  + f_12 * hd_58[k]
                  + f_11 * hd_100[k]
                  - f_12 * hd_112[k];
    }

#pragma omp simd aligned(hd_12, hd_15, hd_17, hd_54, hd_57, hd_59, hd_96, hd_99, hd_101, \
                         hd_108, hd_111, hd_113 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_37[k] = f_57 * hd_12[k]
                  + f_57 * hd_15[k]
                  - f_33 * hd_17[k]
                  - f_33 * hd_54[k]
                  - f_33 * hd_57[k]
                  + f_34 * hd_59[k]
                  - f_57 * hd_96[k]
                  - f_57 * hd_99[k]
                  + f_33 * hd_101[k]
                  + f_33 * hd_108[k]
                  + f_33 * hd_111[k]
                  - f_34 * hd_113[k];
    }

#pragma omp simd aligned(hd_12, hd_14, hd_15, hd_54, hd_56, hd_57, hd_96, hd_98, hd_99, \
                         hd_108, hd_110, hd_111 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = -f_11 * hd_14[k]
                  + f_12 * hd_56[k]
                  + f_11 * hd_98[k]
                  - f_12 * hd_110[k];

        g_39[k] = -f_58 * hd_12[k]
                  + f_58 * hd_15[k]
                  + f_11 * hd_54[k]
                  - f_11 * hd_57[k]
                  + f_58 * hd_96[k]
                  - f_58 * hd_99[k]
                  - f_11 * hd_108[k]
                  + f_11 * hd_111[k];
    }

#pragma omp simd aligned(hd_1, hd_4, hd_19, hd_22, hd_31, hd_34, hd_61, hd_64, hd_73, \
                         hd_76 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_40[k] = -f_17 * hd_1[k]
                  + f_15 * hd_19[k]
                  + f_18 * hd_31[k]
                  + f_14 * hd_61[k]
                  - f_16 * hd_73[k];

        g_41[k] = -f_17 * hd_4[k]
                  + f_15 * hd_22[k]
                  + f_18 * hd_34[k]
                  + f_14 * hd_64[k]
                  - f_16 * hd_76[k];
    }

#pragma omp simd aligned(hd_0, hd_3, hd_5, hd_18, hd_21, hd_23, hd_30, hd_33, hd_35, hd_60, \
                         hd_63, hd_65, hd_72, hd_75, hd_77 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = f_25 * hd_0[k]
                  + f_25 * hd_3[k]
                  - f_21 * hd_5[k]
                  - f_21 * hd_18[k]
                  - f_21 * hd_21[k]
                  + f_22 * hd_23[k]
                  - f_26 * hd_30[k]
                  - f_26 * hd_33[k]
                  + f_27 * hd_35[k]
                  - f_19 * hd_60[k]
                  - f_19 * hd_63[k]
                  + f_20 * hd_65[k]
                  + f_23 * hd_72[k]
                  + f_23 * hd_75[k]
                  - f_24 * hd_77[k];
    }

#pragma omp simd aligned(hd_2, hd_20, hd_32, hd_62, hd_74 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_43[k] = -f_17 * hd_2[k]
                  + f_15 * hd_20[k]
                  + f_18 * hd_32[k]
                  + f_14 * hd_62[k]
                  - f_16 * hd_74[k];
    }

#pragma omp simd aligned(hd_0, hd_3, hd_13, hd_18, hd_21, hd_30, hd_33, hd_43, hd_60, hd_63, \
                         hd_72, hd_75, hd_97 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_44[k] = -f_30 * hd_0[k]
                  + f_30 * hd_3[k]
                  + f_17 * hd_18[k]
                  - f_17 * hd_21[k]
                  + f_31 * hd_30[k]
                  - f_31 * hd_33[k]
                  + f_28 * hd_60[k]
                  - f_28 * hd_63[k]
                  - f_29 * hd_72[k]
                  + f_29 * hd_75[k];

        g_45[k] = f_59 * hd_13[k]
                  - f_60 * hd_43[k]
                  + f_59 * hd_97[k];
    }

#pragma omp simd aligned(hd_12, hd_15, hd_16, hd_17, hd_42, hd_45, hd_46, hd_47, hd_96, hd_99, \
                         hd_100, hd_101 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_46[k] = f_59 * hd_16[k]
                  - f_60 * hd_46[k]
                  + f_59 * hd_100[k];

        g_47[k] = -f_61 * hd_12[k]
                  - f_61 * hd_15[k]
                  + f_58 * hd_17[k]
                  + f_62 * hd_42[k]
                  + f_62 * hd_45[k]
                  - f_63 * hd_47[k]
                  - f_61 * hd_96[k]
                  - f_61 * hd_99[k]
                  + f_58 * hd_101[k];
    }

#pragma omp simd aligned(hd_1, hd_12, hd_14, hd_15, hd_19, hd_42, hd_44, hd_45, hd_61, hd_96, \
                         hd_98, hd_99 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_48[k] = f_59 * hd_14[k]
                  - f_60 * hd_44[k]
                  + f_59 * hd_98[k];

        g_49[k] = f_64 * hd_12[k]
                  - f_64 * hd_15[k]
                  - f_65 * hd_42[k]
                  + f_65 * hd_45[k]
                  + f_64 * hd_96[k]
                  - f_64 * hd_99[k];

        g_50[k] = f_2 * hd_1[k]
                  - f_1 * hd_19[k]
                  + f_0 * hd_61[k];
    }

#pragma omp simd aligned(hd_0, hd_3, hd_4, hd_5, hd_18, hd_21, hd_22, hd_23, hd_60, hd_63, \
                         hd_64, hd_65 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_51[k] = f_2 * hd_4[k]
                  - f_1 * hd_22[k]
                  + f_0 * hd_64[k];

        g_52[k] = -f_6 * hd_0[k]
                  - f_6 * hd_3[k]
                  + f_7 * hd_5[k]
                  + f_4 * hd_18[k]
                  + f_4 * hd_21[k]
                  - f_5 * hd_23[k]
                  - f_3 * hd_60[k]
                  - f_3 * hd_63[k]
                  + f_4 * hd_65[k];
    }

#pragma omp simd aligned(hd_0, hd_2, hd_3, hd_18, hd_20, hd_21, hd_60, hd_62, \
                         hd_63 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_53[k] = f_2 * hd_2[k]
                  - f_1 * hd_20[k]
                  + f_0 * hd_62[k];

        g_54[k] = f_9 * hd_0[k]
                  - f_9 * hd_3[k]
                  - f_0 * hd_18[k]
                  + f_0 * hd_21[k]
                  + f_8 * hd_60[k]
                  - f_8 * hd_63[k];
    }
}

}  // namespace simdtrf
