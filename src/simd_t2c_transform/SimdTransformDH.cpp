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


#include "SimdTransformDH.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_dh(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t dh,
             const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 0.9375 * std::sqrt(42.0);
    const auto f_1 = 1.875 * std::sqrt(42.0);
    const auto f_2 = 0.1875 * std::sqrt(42.0);
    const auto f_3 = 1.5 * std::sqrt(105.0);
    const auto f_4 = 0.1875 * std::sqrt(210.0);
    const auto f_5 = 0.125 * std::sqrt(210.0);
    const auto f_6 = 1.5 * std::sqrt(210.0);
    const auto f_7 = 0.0625 * std::sqrt(210.0);
    const auto f_8 = 0.5 * std::sqrt(210.0);
    const auto f_9 = 1.5 * std::sqrt(35.0);
    const auto f_10 = 3.0 * std::sqrt(35.0);
    const auto f_11 = 0.375 * std::sqrt(5.0);
    const auto f_12 = 0.75 * std::sqrt(5.0);
    const auto f_13 = 4.5 * std::sqrt(5.0);
    const auto f_14 = 3.0 * std::sqrt(5.0);
    const auto f_15 = 1.875 * std::sqrt(3.0);
    const auto f_16 = 3.75 * std::sqrt(3.0);
    const auto f_17 = 5.0 * std::sqrt(3.0);
    const auto f_18 = std::sqrt(3.0);
    const auto f_19 = 0.75 * std::sqrt(35.0);
    const auto f_20 = 0.375 * std::sqrt(105.0);
    const auto f_21 = 2.25 * std::sqrt(105.0);
    const auto f_22 = 0.46875 * std::sqrt(14.0);
    const auto f_23 = 0.9375 * std::sqrt(14.0);
    const auto f_24 = 0.09375 * std::sqrt(14.0);
    const auto f_25 = 1.875 * std::sqrt(14.0);
    const auto f_26 = 0.1875 * std::sqrt(14.0);
    const auto f_27 = 0.09375 * std::sqrt(70.0);
    const auto f_28 = 0.0625 * std::sqrt(70.0);
    const auto f_29 = 0.75 * std::sqrt(70.0);
    const auto f_30 = 0.03125 * std::sqrt(70.0);
    const auto f_31 = 0.25 * std::sqrt(70.0);
    const auto f_32 = 0.1875 * std::sqrt(70.0);
    const auto f_33 = 0.125 * std::sqrt(70.0);
    const auto f_34 = 1.5 * std::sqrt(70.0);
    const auto f_35 = 0.5 * std::sqrt(70.0);
    const auto f_36 = 0.25 * std::sqrt(105.0);
    const auto f_37 = 0.5 * std::sqrt(105.0);
    const auto f_38 = std::sqrt(105.0);
    const auto f_39 = 0.0625 * std::sqrt(15.0);
    const auto f_40 = 0.125 * std::sqrt(15.0);
    const auto f_41 = 0.75 * std::sqrt(15.0);
    const auto f_42 = 0.5 * std::sqrt(15.0);
    const auto f_43 = 0.25 * std::sqrt(15.0);
    const auto f_44 = 1.5 * std::sqrt(15.0);
    const auto f_45 = std::sqrt(15.0);
    const auto f_46 = 0.125 * std::sqrt(105.0);
    const auto f_47 = 0.1875 * std::sqrt(35.0);
    const auto f_48 = 1.125 * std::sqrt(35.0);
    const auto f_49 = 0.375 * std::sqrt(35.0);
    const auto f_50 = 2.25 * std::sqrt(35.0);
    const auto f_51 = 0.46875 * std::sqrt(42.0);
    const auto f_52 = 0.09375 * std::sqrt(42.0);
    const auto f_53 = 0.75 * std::sqrt(105.0);
    const auto f_54 = 0.09375 * std::sqrt(210.0);
    const auto f_55 = 0.75 * std::sqrt(210.0);
    const auto f_56 = 0.03125 * std::sqrt(210.0);
    const auto f_57 = 0.25 * std::sqrt(210.0);
    const auto f_58 = 0.1875 * std::sqrt(5.0);
    const auto f_59 = 2.25 * std::sqrt(5.0);
    const auto f_60 = 1.5 * std::sqrt(5.0);
    const auto f_61 = 0.9375 * std::sqrt(3.0);
    const auto f_62 = 2.5 * std::sqrt(3.0);
    const auto f_63 = 0.5 * std::sqrt(3.0);
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

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);
    const auto *dh_3 = buffer.data(dh + 3);
    const auto *dh_4 = buffer.data(dh + 4);
    const auto *dh_5 = buffer.data(dh + 5);
    const auto *dh_6 = buffer.data(dh + 6);
    const auto *dh_7 = buffer.data(dh + 7);
    const auto *dh_8 = buffer.data(dh + 8);
    const auto *dh_9 = buffer.data(dh + 9);
    const auto *dh_10 = buffer.data(dh + 10);
    const auto *dh_11 = buffer.data(dh + 11);
    const auto *dh_12 = buffer.data(dh + 12);
    const auto *dh_13 = buffer.data(dh + 13);
    const auto *dh_14 = buffer.data(dh + 14);
    const auto *dh_15 = buffer.data(dh + 15);
    const auto *dh_16 = buffer.data(dh + 16);
    const auto *dh_17 = buffer.data(dh + 17);
    const auto *dh_18 = buffer.data(dh + 18);
    const auto *dh_19 = buffer.data(dh + 19);
    const auto *dh_20 = buffer.data(dh + 20);
    const auto *dh_21 = buffer.data(dh + 21);
    const auto *dh_22 = buffer.data(dh + 22);
    const auto *dh_23 = buffer.data(dh + 23);
    const auto *dh_24 = buffer.data(dh + 24);
    const auto *dh_25 = buffer.data(dh + 25);
    const auto *dh_26 = buffer.data(dh + 26);
    const auto *dh_27 = buffer.data(dh + 27);
    const auto *dh_28 = buffer.data(dh + 28);
    const auto *dh_29 = buffer.data(dh + 29);
    const auto *dh_30 = buffer.data(dh + 30);
    const auto *dh_31 = buffer.data(dh + 31);
    const auto *dh_32 = buffer.data(dh + 32);
    const auto *dh_33 = buffer.data(dh + 33);
    const auto *dh_34 = buffer.data(dh + 34);
    const auto *dh_35 = buffer.data(dh + 35);
    const auto *dh_36 = buffer.data(dh + 36);
    const auto *dh_37 = buffer.data(dh + 37);
    const auto *dh_38 = buffer.data(dh + 38);
    const auto *dh_39 = buffer.data(dh + 39);
    const auto *dh_40 = buffer.data(dh + 40);
    const auto *dh_41 = buffer.data(dh + 41);
    const auto *dh_42 = buffer.data(dh + 42);
    const auto *dh_43 = buffer.data(dh + 43);
    const auto *dh_44 = buffer.data(dh + 44);
    const auto *dh_45 = buffer.data(dh + 45);
    const auto *dh_46 = buffer.data(dh + 46);
    const auto *dh_47 = buffer.data(dh + 47);
    const auto *dh_48 = buffer.data(dh + 48);
    const auto *dh_49 = buffer.data(dh + 49);
    const auto *dh_50 = buffer.data(dh + 50);
    const auto *dh_51 = buffer.data(dh + 51);
    const auto *dh_52 = buffer.data(dh + 52);
    const auto *dh_53 = buffer.data(dh + 53);
    const auto *dh_54 = buffer.data(dh + 54);
    const auto *dh_55 = buffer.data(dh + 55);
    const auto *dh_56 = buffer.data(dh + 56);
    const auto *dh_57 = buffer.data(dh + 57);
    const auto *dh_58 = buffer.data(dh + 58);
    const auto *dh_59 = buffer.data(dh + 59);
    const auto *dh_60 = buffer.data(dh + 60);
    const auto *dh_61 = buffer.data(dh + 61);
    const auto *dh_62 = buffer.data(dh + 62);
    const auto *dh_63 = buffer.data(dh + 63);
    const auto *dh_64 = buffer.data(dh + 64);
    const auto *dh_65 = buffer.data(dh + 65);
    const auto *dh_66 = buffer.data(dh + 66);
    const auto *dh_67 = buffer.data(dh + 67);
    const auto *dh_68 = buffer.data(dh + 68);
    const auto *dh_69 = buffer.data(dh + 69);
    const auto *dh_70 = buffer.data(dh + 70);
    const auto *dh_71 = buffer.data(dh + 71);
    const auto *dh_72 = buffer.data(dh + 72);
    const auto *dh_73 = buffer.data(dh + 73);
    const auto *dh_74 = buffer.data(dh + 74);
    const auto *dh_75 = buffer.data(dh + 75);
    const auto *dh_76 = buffer.data(dh + 76);
    const auto *dh_77 = buffer.data(dh + 77);
    const auto *dh_78 = buffer.data(dh + 78);
    const auto *dh_79 = buffer.data(dh + 79);
    const auto *dh_80 = buffer.data(dh + 80);
    const auto *dh_81 = buffer.data(dh + 81);
    const auto *dh_82 = buffer.data(dh + 82);
    const auto *dh_83 = buffer.data(dh + 83);
    const auto *dh_84 = buffer.data(dh + 84);
    const auto *dh_85 = buffer.data(dh + 85);
    const auto *dh_86 = buffer.data(dh + 86);
    const auto *dh_87 = buffer.data(dh + 87);
    const auto *dh_88 = buffer.data(dh + 88);
    const auto *dh_89 = buffer.data(dh + 89);
    const auto *dh_90 = buffer.data(dh + 90);
    const auto *dh_91 = buffer.data(dh + 91);
    const auto *dh_92 = buffer.data(dh + 92);
    const auto *dh_93 = buffer.data(dh + 93);
    const auto *dh_94 = buffer.data(dh + 94);
    const auto *dh_95 = buffer.data(dh + 95);
    const auto *dh_96 = buffer.data(dh + 96);
    const auto *dh_97 = buffer.data(dh + 97);
    const auto *dh_98 = buffer.data(dh + 98);
    const auto *dh_99 = buffer.data(dh + 99);
    const auto *dh_100 = buffer.data(dh + 100);
    const auto *dh_101 = buffer.data(dh + 101);
    const auto *dh_102 = buffer.data(dh + 102);
    const auto *dh_103 = buffer.data(dh + 103);
    const auto *dh_104 = buffer.data(dh + 104);
    const auto *dh_105 = buffer.data(dh + 105);
    const auto *dh_106 = buffer.data(dh + 106);
    const auto *dh_107 = buffer.data(dh + 107);
    const auto *dh_108 = buffer.data(dh + 108);
    const auto *dh_109 = buffer.data(dh + 109);
    const auto *dh_110 = buffer.data(dh + 110);
    const auto *dh_111 = buffer.data(dh + 111);
    const auto *dh_112 = buffer.data(dh + 112);
    const auto *dh_113 = buffer.data(dh + 113);
    const auto *dh_114 = buffer.data(dh + 114);
    const auto *dh_115 = buffer.data(dh + 115);
    const auto *dh_116 = buffer.data(dh + 116);
    const auto *dh_117 = buffer.data(dh + 117);
    const auto *dh_118 = buffer.data(dh + 118);
    const auto *dh_119 = buffer.data(dh + 119);
    const auto *dh_120 = buffer.data(dh + 120);
    const auto *dh_121 = buffer.data(dh + 121);
    const auto *dh_122 = buffer.data(dh + 122);
    const auto *dh_123 = buffer.data(dh + 123);
    const auto *dh_124 = buffer.data(dh + 124);
    const auto *dh_125 = buffer.data(dh + 125);

#pragma omp simd aligned(dh_22, dh_25, dh_27, dh_29, dh_32, dh_34, dh_36, dh_38, \
                         dh_40 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * dh_22[k]
                 - f_1 * dh_27[k]
                 + f_2 * dh_36[k];

        g_1[k] = f_3 * dh_25[k]
                 - f_3 * dh_32[k];

        g_2[k] = -f_4 * dh_22[k]
                 - f_5 * dh_27[k]
                 + f_6 * dh_29[k]
                 + f_7 * dh_36[k]
                 - f_8 * dh_38[k];

        g_3[k] = -f_9 * dh_25[k]
                 - f_9 * dh_32[k]
                 + f_10 * dh_34[k];

        g_4[k] = f_11 * dh_22[k]
                 + f_12 * dh_27[k]
                 - f_13 * dh_29[k]
                 + f_11 * dh_36[k]
                 - f_13 * dh_38[k]
                 + f_14 * dh_40[k];
    }

#pragma omp simd aligned(dh_21, dh_23, dh_24, dh_26, dh_28, dh_30, dh_31, dh_33, dh_35, dh_37, \
                         dh_39, dh_41 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = f_15 * dh_23[k]
                 + f_16 * dh_28[k]
                 - f_17 * dh_30[k]
                 + f_15 * dh_37[k]
                 - f_17 * dh_39[k]
                 + f_18 * dh_41[k];

        g_6[k] = f_11 * dh_21[k]
                 + f_12 * dh_24[k]
                 - f_13 * dh_26[k]
                 + f_11 * dh_31[k]
                 - f_13 * dh_33[k]
                 + f_14 * dh_35[k];

        g_7[k] = -f_19 * dh_23[k]
                 + f_9 * dh_30[k]
                 + f_19 * dh_37[k]
                 - f_9 * dh_39[k];

        g_8[k] = -f_7 * dh_21[k]
                 + f_5 * dh_24[k]
                 + f_8 * dh_26[k]
                 + f_4 * dh_31[k]
                 - f_6 * dh_33[k];
    }

#pragma omp simd aligned(dh_21, dh_23, dh_24, dh_28, dh_31, dh_37, dh_85, dh_88, dh_90, dh_95, \
                         dh_99 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_9[k] = f_20 * dh_23[k]
                 - f_21 * dh_28[k]
                 + f_20 * dh_37[k];

        g_10[k] = f_2 * dh_21[k]
                  - f_1 * dh_24[k]
                  + f_0 * dh_31[k];

        g_11[k] = f_0 * dh_85[k]
                  - f_1 * dh_90[k]
                  + f_2 * dh_99[k];

        g_12[k] = f_3 * dh_88[k]
                  - f_3 * dh_95[k];
    }

#pragma omp simd aligned(dh_85, dh_88, dh_90, dh_92, dh_95, dh_97, dh_99, dh_101, \
                         dh_103 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = -f_4 * dh_85[k]
                  - f_5 * dh_90[k]
                  + f_6 * dh_92[k]
                  + f_7 * dh_99[k]
                  - f_8 * dh_101[k];

        g_14[k] = -f_9 * dh_88[k]
                  - f_9 * dh_95[k]
                  + f_10 * dh_97[k];

        g_15[k] = f_11 * dh_85[k]
                  + f_12 * dh_90[k]
                  - f_13 * dh_92[k]
                  + f_11 * dh_99[k]
                  - f_13 * dh_101[k]
                  + f_14 * dh_103[k];
    }

#pragma omp simd aligned(dh_84, dh_86, dh_87, dh_89, dh_91, dh_93, dh_94, dh_96, dh_98, \
                         dh_100, dh_102, dh_104 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_16[k] = f_15 * dh_86[k]
                  + f_16 * dh_91[k]
                  - f_17 * dh_93[k]
                  + f_15 * dh_100[k]
                  - f_17 * dh_102[k]
                  + f_18 * dh_104[k];

        g_17[k] = f_11 * dh_84[k]
                  + f_12 * dh_87[k]
                  - f_13 * dh_89[k]
                  + f_11 * dh_94[k]
                  - f_13 * dh_96[k]
                  + f_14 * dh_98[k];

        g_18[k] = -f_19 * dh_86[k]
                  + f_9 * dh_93[k]
                  + f_19 * dh_100[k]
                  - f_9 * dh_102[k];

        g_19[k] = -f_7 * dh_84[k]
                  + f_5 * dh_87[k]
                  + f_8 * dh_89[k]
                  + f_4 * dh_94[k]
                  - f_6 * dh_96[k];
    }

#pragma omp simd aligned(dh_84, dh_86, dh_87, dh_91, dh_94, dh_100 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = f_20 * dh_86[k]
                  - f_21 * dh_91[k]
                  + f_20 * dh_100[k];

        g_21[k] = f_2 * dh_84[k]
                  - f_1 * dh_87[k]
                  + f_0 * dh_94[k];
    }

#pragma omp simd aligned(dh_1, dh_6, dh_15, dh_64, dh_69, dh_78, dh_106, dh_111, \
                         dh_120 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = -f_22 * dh_1[k]
                  + f_23 * dh_6[k]
                  - f_24 * dh_15[k]
                  - f_22 * dh_64[k]
                  + f_23 * dh_69[k]
                  - f_24 * dh_78[k]
                  + f_23 * dh_106[k]
                  - f_25 * dh_111[k]
                  + f_26 * dh_120[k];
    }

#pragma omp simd aligned(dh_4, dh_11, dh_67, dh_74, dh_109, dh_116 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = -f_19 * dh_4[k]
                  + f_19 * dh_11[k]
                  - f_19 * dh_67[k]
                  + f_19 * dh_74[k]
                  + f_9 * dh_109[k]
                  - f_9 * dh_116[k];
    }

#pragma omp simd aligned(dh_1, dh_6, dh_8, dh_15, dh_17, dh_64, dh_69, dh_71, dh_78, dh_80, \
                         dh_106, dh_111, dh_113, dh_120, dh_122 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = f_27 * dh_1[k]
                  + f_28 * dh_6[k]
                  - f_29 * dh_8[k]
                  - f_30 * dh_15[k]
                  + f_31 * dh_17[k]
                  + f_27 * dh_64[k]
                  + f_28 * dh_69[k]
                  - f_29 * dh_71[k]
                  - f_30 * dh_78[k]
                  + f_31 * dh_80[k]
                  - f_32 * dh_106[k]
                  - f_33 * dh_111[k]
                  + f_34 * dh_113[k]
                  + f_28 * dh_120[k]
                  - f_35 * dh_122[k];
    }

#pragma omp simd aligned(dh_4, dh_11, dh_13, dh_67, dh_74, dh_76, dh_109, dh_116, \
                         dh_118 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = f_36 * dh_4[k]
                  + f_36 * dh_11[k]
                  - f_37 * dh_13[k]
                  + f_36 * dh_67[k]
                  + f_36 * dh_74[k]
                  - f_37 * dh_76[k]
                  - f_37 * dh_109[k]
                  - f_37 * dh_116[k]
                  + f_38 * dh_118[k];
    }

#pragma omp simd aligned(dh_1, dh_6, dh_8, dh_15, dh_17, dh_19, dh_64, dh_69, dh_71, dh_78, \
                         dh_80, dh_82, dh_106, dh_111, dh_113, dh_120, dh_122, \
                         dh_124 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = -f_39 * dh_1[k]
                  - f_40 * dh_6[k]
                  + f_41 * dh_8[k]
                  - f_39 * dh_15[k]
                  + f_41 * dh_17[k]
                  - f_42 * dh_19[k]
                  - f_39 * dh_64[k]
                  - f_40 * dh_69[k]
                  + f_41 * dh_71[k]
                  - f_39 * dh_78[k]
                  + f_41 * dh_80[k]
                  - f_42 * dh_82[k]
                  + f_40 * dh_106[k]
                  + f_43 * dh_111[k]
                  - f_44 * dh_113[k]
                  + f_40 * dh_120[k]
                  - f_44 * dh_122[k]
                  + f_45 * dh_124[k];
    }

#pragma omp simd aligned(dh_2, dh_7, dh_9, dh_16, dh_18, dh_20, dh_65, dh_70, dh_72, dh_79, \
                         dh_81, dh_83, dh_107, dh_112, dh_114, dh_121, dh_123, \
                         dh_125 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = -0.9375 * dh_2[k]
                  - 1.875 * dh_7[k]
                  + 2.5 * dh_9[k]
                  - 0.9375 * dh_16[k]
                  + 2.5 * dh_18[k]
                  - 0.5 * dh_20[k]
                  - 0.9375 * dh_65[k]
                  - 1.875 * dh_70[k]
                  + 2.5 * dh_72[k]
                  - 0.9375 * dh_79[k]
                  + 2.5 * dh_81[k]
                  - 0.5 * dh_83[k]
                  + 1.875 * dh_107[k]
                  + 3.75 * dh_112[k]
                  - 5.0 * dh_114[k]
                  + 1.875 * dh_121[k]
                  - 5.0 * dh_123[k]
                  + dh_125[k];
    }

#pragma omp simd aligned(dh_0, dh_3, dh_5, dh_10, dh_12, dh_14, dh_63, dh_66, dh_68, dh_73, \
                         dh_75, dh_77, dh_105, dh_108, dh_110, dh_115, dh_117, \
                         dh_119 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = -f_39 * dh_0[k]
                  - f_40 * dh_3[k]
                  + f_41 * dh_5[k]
                  - f_39 * dh_10[k]
                  + f_41 * dh_12[k]
                  - f_42 * dh_14[k]
                  - f_39 * dh_63[k]
                  - f_40 * dh_66[k]
                  + f_41 * dh_68[k]
                  - f_39 * dh_73[k]
                  + f_41 * dh_75[k]
                  - f_42 * dh_77[k]
                  + f_40 * dh_105[k]
                  + f_43 * dh_108[k]
                  - f_44 * dh_110[k]
                  + f_40 * dh_115[k]
                  - f_44 * dh_117[k]
                  + f_45 * dh_119[k];
    }

#pragma omp simd aligned(dh_2, dh_9, dh_16, dh_18, dh_65, dh_72, dh_79, dh_81, dh_107, dh_114, \
                         dh_121, dh_123 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = f_46 * dh_2[k]
                  - f_36 * dh_9[k]
                  - f_46 * dh_16[k]
                  + f_36 * dh_18[k]
                  + f_46 * dh_65[k]
                  - f_36 * dh_72[k]
                  - f_46 * dh_79[k]
                  + f_36 * dh_81[k]
                  - f_36 * dh_107[k]
                  + f_37 * dh_114[k]
                  + f_36 * dh_121[k]
                  - f_37 * dh_123[k];
    }

#pragma omp simd aligned(dh_0, dh_3, dh_5, dh_10, dh_12, dh_63, dh_66, dh_68, dh_73, dh_75, \
                         dh_105, dh_108, dh_110, dh_115, dh_117 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = f_30 * dh_0[k]
                  - f_28 * dh_3[k]
                  - f_31 * dh_5[k]
                  - f_27 * dh_10[k]
                  + f_29 * dh_12[k]
                  + f_30 * dh_63[k]
                  - f_28 * dh_66[k]
                  - f_31 * dh_68[k]
                  - f_27 * dh_73[k]
                  + f_29 * dh_75[k]
                  - f_28 * dh_105[k]
                  + f_33 * dh_108[k]
                  + f_35 * dh_110[k]
                  + f_32 * dh_115[k]
                  - f_34 * dh_117[k];
    }

#pragma omp simd aligned(dh_2, dh_7, dh_16, dh_65, dh_70, dh_79, dh_107, dh_112, \
                         dh_121 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = -f_47 * dh_2[k]
                  + f_48 * dh_7[k]
                  - f_47 * dh_16[k]
                  - f_47 * dh_65[k]
                  + f_48 * dh_70[k]
                  - f_47 * dh_79[k]
                  + f_49 * dh_107[k]
                  - f_50 * dh_112[k]
                  + f_49 * dh_121[k];
    }

#pragma omp simd aligned(dh_0, dh_3, dh_10, dh_43, dh_48, dh_57, dh_63, dh_66, dh_73, dh_105, \
                         dh_108, dh_115 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = -f_24 * dh_0[k]
                  + f_23 * dh_3[k]
                  - f_22 * dh_10[k]
                  - f_24 * dh_63[k]
                  + f_23 * dh_66[k]
                  - f_22 * dh_73[k]
                  + f_26 * dh_105[k]
                  - f_25 * dh_108[k]
                  + f_23 * dh_115[k];

        g_33[k] = f_0 * dh_43[k]
                  - f_1 * dh_48[k]
                  + f_2 * dh_57[k];
    }

#pragma omp simd aligned(dh_43, dh_46, dh_48, dh_50, dh_53, dh_55, dh_57, dh_59, \
                         dh_61 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = f_3 * dh_46[k]
                  - f_3 * dh_53[k];

        g_35[k] = -f_4 * dh_43[k]
                  - f_5 * dh_48[k]
                  + f_6 * dh_50[k]
                  + f_7 * dh_57[k]
                  - f_8 * dh_59[k];

        g_36[k] = -f_9 * dh_46[k]
                  - f_9 * dh_53[k]
                  + f_10 * dh_55[k];

        g_37[k] = f_11 * dh_43[k]
                  + f_12 * dh_48[k]
                  - f_13 * dh_50[k]
                  + f_11 * dh_57[k]
                  - f_13 * dh_59[k]
                  + f_14 * dh_61[k];
    }

#pragma omp simd aligned(dh_42, dh_44, dh_45, dh_47, dh_49, dh_51, dh_52, dh_54, dh_56, dh_58, \
                         dh_60, dh_62 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = f_15 * dh_44[k]
                  + f_16 * dh_49[k]
                  - f_17 * dh_51[k]
                  + f_15 * dh_58[k]
                  - f_17 * dh_60[k]
                  + f_18 * dh_62[k];

        g_39[k] = f_11 * dh_42[k]
                  + f_12 * dh_45[k]
                  - f_13 * dh_47[k]
                  + f_11 * dh_52[k]
                  - f_13 * dh_54[k]
                  + f_14 * dh_56[k];

        g_40[k] = -f_19 * dh_44[k]
                  + f_9 * dh_51[k]
                  + f_19 * dh_58[k]
                  - f_9 * dh_60[k];

        g_41[k] = -f_7 * dh_42[k]
                  + f_5 * dh_45[k]
                  + f_8 * dh_47[k]
                  + f_4 * dh_52[k]
                  - f_6 * dh_54[k];
    }

#pragma omp simd aligned(dh_1, dh_6, dh_15, dh_42, dh_44, dh_45, dh_49, dh_52, dh_58, dh_64, \
                         dh_69, dh_78 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = f_20 * dh_44[k]
                  - f_21 * dh_49[k]
                  + f_20 * dh_58[k];

        g_43[k] = f_2 * dh_42[k]
                  - f_1 * dh_45[k]
                  + f_0 * dh_52[k];

        g_44[k] = f_51 * dh_1[k]
                  - f_0 * dh_6[k]
                  + f_52 * dh_15[k]
                  - f_51 * dh_64[k]
                  + f_0 * dh_69[k]
                  - f_52 * dh_78[k];
    }

#pragma omp simd aligned(dh_1, dh_4, dh_6, dh_8, dh_11, dh_15, dh_17, dh_64, dh_67, dh_69, \
                         dh_71, dh_74, dh_78, dh_80 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_45[k] = f_53 * dh_4[k]
                  - f_53 * dh_11[k]
                  - f_53 * dh_67[k]
                  + f_53 * dh_74[k];

        g_46[k] = -f_54 * dh_1[k]
                  - f_7 * dh_6[k]
                  + f_55 * dh_8[k]
                  + f_56 * dh_15[k]
                  - f_57 * dh_17[k]
                  + f_54 * dh_64[k]
                  + f_7 * dh_69[k]
                  - f_55 * dh_71[k]
                  - f_56 * dh_78[k]
                  + f_57 * dh_80[k];
    }

#pragma omp simd aligned(dh_4, dh_11, dh_13, dh_67, dh_74, dh_76 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_47[k] = -f_19 * dh_4[k]
                  - f_19 * dh_11[k]
                  + f_9 * dh_13[k]
                  + f_19 * dh_67[k]
                  + f_19 * dh_74[k]
                  - f_9 * dh_76[k];
    }

#pragma omp simd aligned(dh_1, dh_6, dh_8, dh_15, dh_17, dh_19, dh_64, dh_69, dh_71, dh_78, \
                         dh_80, dh_82 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_48[k] = f_58 * dh_1[k]
                  + f_11 * dh_6[k]
                  - f_59 * dh_8[k]
                  + f_58 * dh_15[k]
                  - f_59 * dh_17[k]
                  + f_60 * dh_19[k]
                  - f_58 * dh_64[k]
                  - f_11 * dh_69[k]
                  + f_59 * dh_71[k]
                  - f_58 * dh_78[k]
                  + f_59 * dh_80[k]
                  - f_60 * dh_82[k];
    }

#pragma omp simd aligned(dh_2, dh_7, dh_9, dh_16, dh_18, dh_20, dh_65, dh_70, dh_72, dh_79, \
                         dh_81, dh_83 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_49[k] = f_61 * dh_2[k]
                  + f_15 * dh_7[k]
                  - f_62 * dh_9[k]
                  + f_61 * dh_16[k]
                  - f_62 * dh_18[k]
                  + f_63 * dh_20[k]
                  - f_61 * dh_65[k]
                  - f_15 * dh_70[k]
                  + f_62 * dh_72[k]
                  - f_61 * dh_79[k]
                  + f_62 * dh_81[k]
                  - f_63 * dh_83[k];
    }

#pragma omp simd aligned(dh_0, dh_3, dh_5, dh_10, dh_12, dh_14, dh_63, dh_66, dh_68, dh_73, \
                         dh_75, dh_77 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_50[k] = f_58 * dh_0[k]
                  + f_11 * dh_3[k]
                  - f_59 * dh_5[k]
                  + f_58 * dh_10[k]
                  - f_59 * dh_12[k]
                  + f_60 * dh_14[k]
                  - f_58 * dh_63[k]
                  - f_11 * dh_66[k]
                  + f_59 * dh_68[k]
                  - f_58 * dh_73[k]
                  + f_59 * dh_75[k]
                  - f_60 * dh_77[k];
    }

#pragma omp simd aligned(dh_2, dh_9, dh_16, dh_18, dh_65, dh_72, dh_79, \
                         dh_81 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_51[k] = -f_49 * dh_2[k]
                  + f_19 * dh_9[k]
                  + f_49 * dh_16[k]
                  - f_19 * dh_18[k]
                  + f_49 * dh_65[k]
                  - f_19 * dh_72[k]
                  - f_49 * dh_79[k]
                  + f_19 * dh_81[k];
    }

#pragma omp simd aligned(dh_0, dh_3, dh_5, dh_10, dh_12, dh_63, dh_66, dh_68, dh_73, \
                         dh_75 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_52[k] = -f_56 * dh_0[k]
                  + f_7 * dh_3[k]
                  + f_57 * dh_5[k]
                  + f_54 * dh_10[k]
                  - f_55 * dh_12[k]
                  + f_56 * dh_63[k]
                  - f_7 * dh_66[k]
                  - f_57 * dh_68[k]
                  - f_54 * dh_73[k]
                  + f_55 * dh_75[k];
    }

#pragma omp simd aligned(dh_0, dh_2, dh_3, dh_7, dh_10, dh_16, dh_63, dh_65, dh_66, dh_70, \
                         dh_73, dh_79 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_53[k] = f_64 * dh_2[k]
                  - f_65 * dh_7[k]
                  + f_64 * dh_16[k]
                  - f_64 * dh_65[k]
                  + f_65 * dh_70[k]
                  - f_64 * dh_79[k];

        g_54[k] = f_52 * dh_0[k]
                  - f_0 * dh_3[k]
                  + f_51 * dh_10[k]
                  - f_52 * dh_63[k]
                  + f_0 * dh_66[k]
                  - f_51 * dh_73[k];
    }
}

}  // namespace simdtrf
