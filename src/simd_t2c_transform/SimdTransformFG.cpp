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


#include "SimdTransformFG.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_fg(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t fg,
             const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 1.875 * std::sqrt(14.0);
    const auto f_1 = 0.625 * std::sqrt(14.0);
    const auto f_2 = 5.625 * std::sqrt(7.0);
    const auto f_3 = 1.875 * std::sqrt(7.0);
    const auto f_4 = 0.625 * std::sqrt(7.0);
    const auto f_5 = 1.875 * std::sqrt(2.0);
    const auto f_6 = 11.25 * std::sqrt(2.0);
    const auto f_7 = 0.625 * std::sqrt(2.0);
    const auto f_8 = 3.75 * std::sqrt(2.0);
    const auto f_9 = 0.28125 * std::sqrt(10.0);
    const auto f_10 = 0.5625 * std::sqrt(10.0);
    const auto f_11 = 2.25 * std::sqrt(10.0);
    const auto f_12 = 0.75 * std::sqrt(10.0);
    const auto f_13 = 0.09375 * std::sqrt(10.0);
    const auto f_14 = 0.1875 * std::sqrt(10.0);
    const auto f_15 = 0.25 * std::sqrt(10.0);
    const auto f_16 = 0.9375 * std::sqrt(2.0);
    const auto f_17 = 5.625 * std::sqrt(2.0);
    const auto f_18 = 0.3125 * std::sqrt(2.0);
    const auto f_19 = 0.46875 * std::sqrt(14.0);
    const auto f_20 = 2.8125 * std::sqrt(14.0);
    const auto f_21 = 0.15625 * std::sqrt(14.0);
    const auto f_22 = 0.9375 * std::sqrt(14.0);
    const auto f_23 = 2.5 * std::sqrt(21.0);
    const auto f_24 = 3.75 * std::sqrt(42.0);
    const auto f_25 = 1.25 * std::sqrt(42.0);
    const auto f_26 = 2.5 * std::sqrt(3.0);
    const auto f_27 = 15.0 * std::sqrt(3.0);
    const auto f_28 = 3.75 * std::sqrt(6.0);
    const auto f_29 = 5.0 * std::sqrt(6.0);
    const auto f_30 = 0.375 * std::sqrt(15.0);
    const auto f_31 = 0.75 * std::sqrt(15.0);
    const auto f_32 = 3.0 * std::sqrt(15.0);
    const auto f_33 = std::sqrt(15.0);
    const auto f_34 = 1.25 * std::sqrt(3.0);
    const auto f_35 = 7.5 * std::sqrt(3.0);
    const auto f_36 = 0.625 * std::sqrt(21.0);
    const auto f_37 = 3.75 * std::sqrt(21.0);
    const auto f_38 = 0.125 * std::sqrt(210.0);
    const auto f_39 = 0.5 * std::sqrt(210.0);
    const auto f_40 = 0.375 * std::sqrt(105.0);
    const auto f_41 = 0.125 * std::sqrt(105.0);
    const auto f_42 = 1.5 * std::sqrt(105.0);
    const auto f_43 = 0.5 * std::sqrt(105.0);
    const auto f_44 = 0.125 * std::sqrt(30.0);
    const auto f_45 = 0.75 * std::sqrt(30.0);
    const auto f_46 = 0.5 * std::sqrt(30.0);
    const auto f_47 = 3.0 * std::sqrt(30.0);
    const auto f_48 = 0.5 * std::sqrt(15.0);
    const auto f_49 = 1.5 * std::sqrt(15.0);
    const auto f_50 = 2.0 * std::sqrt(15.0);
    const auto f_51 = 0.09375 * std::sqrt(6.0);
    const auto f_52 = 0.1875 * std::sqrt(6.0);
    const auto f_53 = 0.75 * std::sqrt(6.0);
    const auto f_54 = 0.25 * std::sqrt(6.0);
    const auto f_55 = 0.375 * std::sqrt(6.0);
    const auto f_56 = 3.0 * std::sqrt(6.0);
    const auto f_57 = std::sqrt(6.0);
    const auto f_58 = 0.0625 * std::sqrt(30.0);
    const auto f_59 = 0.375 * std::sqrt(30.0);
    const auto f_60 = 0.25 * std::sqrt(30.0);
    const auto f_61 = 1.5 * std::sqrt(30.0);
    const auto f_62 = 0.03125 * std::sqrt(210.0);
    const auto f_63 = 0.1875 * std::sqrt(210.0);
    const auto f_64 = 0.75 * std::sqrt(210.0);
    const auto f_65 = 0.75 * std::sqrt(35.0);
    const auto f_66 = 0.5 * std::sqrt(35.0);
    const auto f_67 = 1.125 * std::sqrt(70.0);
    const auto f_68 = 0.375 * std::sqrt(70.0);
    const auto f_69 = 0.75 * std::sqrt(70.0);
    const auto f_70 = 0.25 * std::sqrt(70.0);
    const auto f_71 = 0.75 * std::sqrt(5.0);
    const auto f_72 = 4.5 * std::sqrt(5.0);
    const auto f_73 = 0.5 * std::sqrt(5.0);
    const auto f_74 = 3.0 * std::sqrt(5.0);
    const auto f_75 = 1.125 * std::sqrt(10.0);
    const auto f_76 = 1.5 * std::sqrt(10.0);
    const auto f_77 = std::sqrt(10.0);
    const auto f_78 = 0.375 * std::sqrt(5.0);
    const auto f_79 = 2.25 * std::sqrt(5.0);
    const auto f_80 = 0.25 * std::sqrt(5.0);
    const auto f_81 = 1.5 * std::sqrt(5.0);
    const auto f_82 = 0.1875 * std::sqrt(35.0);
    const auto f_83 = 1.125 * std::sqrt(35.0);
    const auto f_84 = 0.125 * std::sqrt(35.0);
    const auto f_85 = 1.25 * std::sqrt(21.0);
    const auto f_86 = 1.875 * std::sqrt(42.0);
    const auto f_87 = 0.625 * std::sqrt(42.0);
    const auto f_88 = 1.875 * std::sqrt(6.0);
    const auto f_89 = 2.5 * std::sqrt(6.0);
    const auto f_90 = 0.1875 * std::sqrt(15.0);
    const auto f_91 = 0.625 * std::sqrt(3.0);
    const auto f_92 = 3.75 * std::sqrt(3.0);
    const auto f_93 = 0.3125 * std::sqrt(21.0);
    const auto f_94 = 1.875 * std::sqrt(21.0);

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
    auto *g_55 = values + 55 * nvalues;
    auto *g_56 = values + 56 * nvalues;
    auto *g_57 = values + 57 * nvalues;
    auto *g_58 = values + 58 * nvalues;
    auto *g_59 = values + 59 * nvalues;
    auto *g_60 = values + 60 * nvalues;
    auto *g_61 = values + 61 * nvalues;
    auto *g_62 = values + 62 * nvalues;

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_11 = buffer.data(fg + 11);
    const auto *fg_12 = buffer.data(fg + 12);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_14 = buffer.data(fg + 14);
    const auto *fg_15 = buffer.data(fg + 15);
    const auto *fg_16 = buffer.data(fg + 16);
    const auto *fg_17 = buffer.data(fg + 17);
    const auto *fg_18 = buffer.data(fg + 18);
    const auto *fg_19 = buffer.data(fg + 19);
    const auto *fg_20 = buffer.data(fg + 20);
    const auto *fg_21 = buffer.data(fg + 21);
    const auto *fg_22 = buffer.data(fg + 22);
    const auto *fg_23 = buffer.data(fg + 23);
    const auto *fg_24 = buffer.data(fg + 24);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_26 = buffer.data(fg + 26);
    const auto *fg_27 = buffer.data(fg + 27);
    const auto *fg_28 = buffer.data(fg + 28);
    const auto *fg_29 = buffer.data(fg + 29);
    const auto *fg_30 = buffer.data(fg + 30);
    const auto *fg_31 = buffer.data(fg + 31);
    const auto *fg_32 = buffer.data(fg + 32);
    const auto *fg_33 = buffer.data(fg + 33);
    const auto *fg_34 = buffer.data(fg + 34);
    const auto *fg_35 = buffer.data(fg + 35);
    const auto *fg_36 = buffer.data(fg + 36);
    const auto *fg_37 = buffer.data(fg + 37);
    const auto *fg_38 = buffer.data(fg + 38);
    const auto *fg_39 = buffer.data(fg + 39);
    const auto *fg_40 = buffer.data(fg + 40);
    const auto *fg_41 = buffer.data(fg + 41);
    const auto *fg_42 = buffer.data(fg + 42);
    const auto *fg_43 = buffer.data(fg + 43);
    const auto *fg_44 = buffer.data(fg + 44);
    const auto *fg_45 = buffer.data(fg + 45);
    const auto *fg_46 = buffer.data(fg + 46);
    const auto *fg_47 = buffer.data(fg + 47);
    const auto *fg_48 = buffer.data(fg + 48);
    const auto *fg_49 = buffer.data(fg + 49);
    const auto *fg_50 = buffer.data(fg + 50);
    const auto *fg_51 = buffer.data(fg + 51);
    const auto *fg_52 = buffer.data(fg + 52);
    const auto *fg_53 = buffer.data(fg + 53);
    const auto *fg_54 = buffer.data(fg + 54);
    const auto *fg_55 = buffer.data(fg + 55);
    const auto *fg_56 = buffer.data(fg + 56);
    const auto *fg_57 = buffer.data(fg + 57);
    const auto *fg_58 = buffer.data(fg + 58);
    const auto *fg_59 = buffer.data(fg + 59);
    const auto *fg_60 = buffer.data(fg + 60);
    const auto *fg_61 = buffer.data(fg + 61);
    const auto *fg_62 = buffer.data(fg + 62);
    const auto *fg_63 = buffer.data(fg + 63);
    const auto *fg_64 = buffer.data(fg + 64);
    const auto *fg_65 = buffer.data(fg + 65);
    const auto *fg_66 = buffer.data(fg + 66);
    const auto *fg_67 = buffer.data(fg + 67);
    const auto *fg_68 = buffer.data(fg + 68);
    const auto *fg_69 = buffer.data(fg + 69);
    const auto *fg_70 = buffer.data(fg + 70);
    const auto *fg_71 = buffer.data(fg + 71);
    const auto *fg_72 = buffer.data(fg + 72);
    const auto *fg_73 = buffer.data(fg + 73);
    const auto *fg_74 = buffer.data(fg + 74);
    const auto *fg_75 = buffer.data(fg + 75);
    const auto *fg_76 = buffer.data(fg + 76);
    const auto *fg_77 = buffer.data(fg + 77);
    const auto *fg_78 = buffer.data(fg + 78);
    const auto *fg_79 = buffer.data(fg + 79);
    const auto *fg_80 = buffer.data(fg + 80);
    const auto *fg_81 = buffer.data(fg + 81);
    const auto *fg_82 = buffer.data(fg + 82);
    const auto *fg_83 = buffer.data(fg + 83);
    const auto *fg_84 = buffer.data(fg + 84);
    const auto *fg_85 = buffer.data(fg + 85);
    const auto *fg_86 = buffer.data(fg + 86);
    const auto *fg_87 = buffer.data(fg + 87);
    const auto *fg_88 = buffer.data(fg + 88);
    const auto *fg_89 = buffer.data(fg + 89);
    const auto *fg_90 = buffer.data(fg + 90);
    const auto *fg_91 = buffer.data(fg + 91);
    const auto *fg_92 = buffer.data(fg + 92);
    const auto *fg_93 = buffer.data(fg + 93);
    const auto *fg_94 = buffer.data(fg + 94);
    const auto *fg_95 = buffer.data(fg + 95);
    const auto *fg_96 = buffer.data(fg + 96);
    const auto *fg_97 = buffer.data(fg + 97);
    const auto *fg_98 = buffer.data(fg + 98);
    const auto *fg_99 = buffer.data(fg + 99);
    const auto *fg_100 = buffer.data(fg + 100);
    const auto *fg_101 = buffer.data(fg + 101);
    const auto *fg_102 = buffer.data(fg + 102);
    const auto *fg_103 = buffer.data(fg + 103);
    const auto *fg_104 = buffer.data(fg + 104);
    const auto *fg_105 = buffer.data(fg + 105);
    const auto *fg_106 = buffer.data(fg + 106);
    const auto *fg_107 = buffer.data(fg + 107);
    const auto *fg_108 = buffer.data(fg + 108);
    const auto *fg_109 = buffer.data(fg + 109);
    const auto *fg_110 = buffer.data(fg + 110);
    const auto *fg_111 = buffer.data(fg + 111);
    const auto *fg_112 = buffer.data(fg + 112);
    const auto *fg_113 = buffer.data(fg + 113);
    const auto *fg_114 = buffer.data(fg + 114);
    const auto *fg_115 = buffer.data(fg + 115);
    const auto *fg_116 = buffer.data(fg + 116);
    const auto *fg_117 = buffer.data(fg + 117);
    const auto *fg_118 = buffer.data(fg + 118);
    const auto *fg_119 = buffer.data(fg + 119);
    const auto *fg_120 = buffer.data(fg + 120);
    const auto *fg_121 = buffer.data(fg + 121);
    const auto *fg_122 = buffer.data(fg + 122);
    const auto *fg_123 = buffer.data(fg + 123);
    const auto *fg_124 = buffer.data(fg + 124);
    const auto *fg_125 = buffer.data(fg + 125);
    const auto *fg_126 = buffer.data(fg + 126);
    const auto *fg_127 = buffer.data(fg + 127);
    const auto *fg_128 = buffer.data(fg + 128);
    const auto *fg_129 = buffer.data(fg + 129);
    const auto *fg_130 = buffer.data(fg + 130);
    const auto *fg_131 = buffer.data(fg + 131);
    const auto *fg_132 = buffer.data(fg + 132);
    const auto *fg_133 = buffer.data(fg + 133);
    const auto *fg_134 = buffer.data(fg + 134);
    const auto *fg_135 = buffer.data(fg + 135);
    const auto *fg_136 = buffer.data(fg + 136);
    const auto *fg_137 = buffer.data(fg + 137);
    const auto *fg_138 = buffer.data(fg + 138);
    const auto *fg_139 = buffer.data(fg + 139);
    const auto *fg_140 = buffer.data(fg + 140);
    const auto *fg_141 = buffer.data(fg + 141);
    const auto *fg_142 = buffer.data(fg + 142);
    const auto *fg_143 = buffer.data(fg + 143);
    const auto *fg_144 = buffer.data(fg + 144);
    const auto *fg_145 = buffer.data(fg + 145);
    const auto *fg_146 = buffer.data(fg + 146);
    const auto *fg_147 = buffer.data(fg + 147);
    const auto *fg_148 = buffer.data(fg + 148);
    const auto *fg_149 = buffer.data(fg + 149);

#pragma omp simd aligned(fg_16, fg_19, fg_21, fg_23, fg_26, fg_28, fg_91, fg_94, fg_96, fg_98, \
                         fg_101, fg_103 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * fg_16[k]
                 - f_0 * fg_21[k]
                 - f_1 * fg_91[k]
                 + f_1 * fg_96[k];

        g_1[k] = f_2 * fg_19[k]
                 - f_3 * fg_26[k]
                 - f_3 * fg_94[k]
                 + f_4 * fg_101[k];

        g_2[k] = -f_5 * fg_16[k]
                 - f_5 * fg_21[k]
                 + f_6 * fg_23[k]
                 + f_7 * fg_91[k]
                 + f_7 * fg_96[k]
                 - f_8 * fg_98[k];

        g_3[k] = -5.625 * fg_19[k]
                 - 5.625 * fg_26[k]
                 + 7.5 * fg_28[k]
                 + 1.875 * fg_94[k]
                 + 1.875 * fg_101[k]
                 - 2.5 * fg_103[k];
    }

#pragma omp simd aligned(fg_15, fg_18, fg_20, fg_25, fg_27, fg_29, fg_90, fg_93, fg_95, \
                         fg_100, fg_102, fg_104 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = f_9 * fg_15[k]
                 + f_10 * fg_18[k]
                 - f_11 * fg_20[k]
                 + f_9 * fg_25[k]
                 - f_11 * fg_27[k]
                 + f_12 * fg_29[k]
                 - f_13 * fg_90[k]
                 - f_14 * fg_93[k]
                 + f_12 * fg_95[k]
                 - f_13 * fg_100[k]
                 + f_12 * fg_102[k]
                 - f_15 * fg_104[k];
    }

#pragma omp simd aligned(fg_15, fg_17, fg_20, fg_22, fg_24, fg_25, fg_27, fg_90, fg_92, fg_95, \
                         fg_97, fg_99, fg_100, fg_102 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = -5.625 * fg_17[k]
                 - 5.625 * fg_22[k]
                 + 7.5 * fg_24[k]
                 + 1.875 * fg_92[k]
                 + 1.875 * fg_97[k]
                 - 2.5 * fg_99[k];

        g_6[k] = -f_16 * fg_15[k]
                 + f_17 * fg_20[k]
                 + f_16 * fg_25[k]
                 - f_17 * fg_27[k]
                 + f_18 * fg_90[k]
                 - f_5 * fg_95[k]
                 - f_18 * fg_100[k]
                 + f_5 * fg_102[k];
    }

#pragma omp simd aligned(fg_15, fg_17, fg_18, fg_22, fg_25, fg_61, fg_66, fg_90, fg_92, fg_93, \
                         fg_97, fg_100 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = f_3 * fg_17[k]
                 - f_2 * fg_22[k]
                 - f_4 * fg_92[k]
                 + f_3 * fg_97[k];

        g_8[k] = f_19 * fg_15[k]
                 - f_20 * fg_18[k]
                 + f_19 * fg_25[k]
                 - f_21 * fg_90[k]
                 + f_22 * fg_93[k]
                 - f_21 * fg_100[k];

        g_9[k] = f_23 * fg_61[k]
                 - f_23 * fg_66[k];
    }

#pragma omp simd aligned(fg_60, fg_61, fg_63, fg_64, fg_65, fg_66, fg_68, fg_70, fg_71, fg_72, \
                         fg_73, fg_74 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = f_24 * fg_64[k]
                  - f_25 * fg_71[k];

        g_11[k] = -f_26 * fg_61[k]
                  - f_26 * fg_66[k]
                  + f_27 * fg_68[k];

        g_12[k] = -f_28 * fg_64[k]
                  - f_28 * fg_71[k]
                  + f_29 * fg_73[k];

        g_13[k] = f_30 * fg_60[k]
                  + f_31 * fg_63[k]
                  - f_32 * fg_65[k]
                  + f_30 * fg_70[k]
                  - f_32 * fg_72[k]
                  + f_33 * fg_74[k];
    }

#pragma omp simd aligned(fg_60, fg_62, fg_63, fg_65, fg_67, fg_69, fg_70, \
                         fg_72 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = -f_28 * fg_62[k]
                  - f_28 * fg_67[k]
                  + f_29 * fg_69[k];

        g_15[k] = -f_34 * fg_60[k]
                  + f_35 * fg_65[k]
                  + f_34 * fg_70[k]
                  - f_35 * fg_72[k];

        g_16[k] = f_25 * fg_62[k]
                  - f_24 * fg_67[k];

        g_17[k] = f_36 * fg_60[k]
                  - f_37 * fg_63[k]
                  + f_36 * fg_70[k];
    }

#pragma omp simd aligned(fg_16, fg_19, fg_21, fg_26, fg_91, fg_94, fg_96, fg_101, fg_121, \
                         fg_124, fg_126, fg_131 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = -f_38 * fg_16[k]
                  + f_38 * fg_21[k]
                  - f_38 * fg_91[k]
                  + f_38 * fg_96[k]
                  + f_39 * fg_121[k]
                  - f_39 * fg_126[k];

        g_19[k] = -f_40 * fg_19[k]
                  + f_41 * fg_26[k]
                  - f_40 * fg_94[k]
                  + f_41 * fg_101[k]
                  + f_42 * fg_124[k]
                  - f_43 * fg_131[k];
    }

#pragma omp simd aligned(fg_16, fg_21, fg_23, fg_91, fg_96, fg_98, fg_121, fg_126, \
                         fg_128 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = f_44 * fg_16[k]
                  + f_44 * fg_21[k]
                  - f_45 * fg_23[k]
                  + f_44 * fg_91[k]
                  + f_44 * fg_96[k]
                  - f_45 * fg_98[k]
                  - f_46 * fg_121[k]
                  - f_46 * fg_126[k]
                  + f_47 * fg_128[k];
    }

#pragma omp simd aligned(fg_19, fg_26, fg_28, fg_94, fg_101, fg_103, fg_124, fg_131, \
                         fg_133 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = f_30 * fg_19[k]
                  + f_30 * fg_26[k]
                  - f_48 * fg_28[k]
                  + f_30 * fg_94[k]
                  + f_30 * fg_101[k]
                  - f_48 * fg_103[k]
                  - f_49 * fg_124[k]
                  - f_49 * fg_131[k]
                  + f_50 * fg_133[k];
    }

#pragma omp simd aligned(fg_15, fg_18, fg_20, fg_25, fg_27, fg_29, fg_90, fg_93, fg_95, \
                         fg_100, fg_102, fg_104, fg_120, fg_123, fg_125, fg_130, fg_132, \
                         fg_134 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = -f_51 * fg_15[k]
                  - f_52 * fg_18[k]
                  + f_53 * fg_20[k]
                  - f_51 * fg_25[k]
                  + f_53 * fg_27[k]
                  - f_54 * fg_29[k]
                  - f_51 * fg_90[k]
                  - f_52 * fg_93[k]
                  + f_53 * fg_95[k]
                  - f_51 * fg_100[k]
                  + f_53 * fg_102[k]
                  - f_54 * fg_104[k]
                  + f_55 * fg_120[k]
                  + f_53 * fg_123[k]
                  - f_56 * fg_125[k]
                  + f_55 * fg_130[k]
                  - f_56 * fg_132[k]
                  + f_57 * fg_134[k];
    }

#pragma omp simd aligned(fg_17, fg_22, fg_24, fg_92, fg_97, fg_99, fg_122, fg_127, \
                         fg_129 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = f_30 * fg_17[k]
                  + f_30 * fg_22[k]
                  - f_48 * fg_24[k]
                  + f_30 * fg_92[k]
                  + f_30 * fg_97[k]
                  - f_48 * fg_99[k]
                  - f_49 * fg_122[k]
                  - f_49 * fg_127[k]
                  + f_50 * fg_129[k];
    }

#pragma omp simd aligned(fg_15, fg_20, fg_25, fg_27, fg_90, fg_95, fg_100, fg_102, fg_120, \
                         fg_125, fg_130, fg_132 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = f_58 * fg_15[k]
                  - f_59 * fg_20[k]
                  - f_58 * fg_25[k]
                  + f_59 * fg_27[k]
                  + f_58 * fg_90[k]
                  - f_59 * fg_95[k]
                  - f_58 * fg_100[k]
                  + f_59 * fg_102[k]
                  - f_60 * fg_120[k]
                  + f_61 * fg_125[k]
                  + f_60 * fg_130[k]
                  - f_61 * fg_132[k];
    }

#pragma omp simd aligned(fg_17, fg_22, fg_92, fg_97, fg_122, fg_127 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = -f_41 * fg_17[k]
                  + f_40 * fg_22[k]
                  - f_41 * fg_92[k]
                  + f_40 * fg_97[k]
                  + f_43 * fg_122[k]
                  - f_42 * fg_127[k];
    }

#pragma omp simd aligned(fg_15, fg_18, fg_25, fg_90, fg_93, fg_100, fg_120, fg_123, \
                         fg_130 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = -f_62 * fg_15[k]
                  + f_63 * fg_18[k]
                  - f_62 * fg_25[k]
                  - f_62 * fg_90[k]
                  + f_63 * fg_93[k]
                  - f_62 * fg_100[k]
                  + f_38 * fg_120[k]
                  - f_64 * fg_123[k]
                  + f_38 * fg_130[k];
    }

#pragma omp simd aligned(fg_31, fg_34, fg_36, fg_41, fg_106, fg_109, fg_111, fg_116, fg_136, \
                         fg_139, fg_141, fg_146 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = -f_65 * fg_31[k]
                  + f_65 * fg_36[k]
                  - f_65 * fg_106[k]
                  + f_65 * fg_111[k]
                  + f_66 * fg_136[k]
                  - f_66 * fg_141[k];

        g_28[k] = -f_67 * fg_34[k]
                  + f_68 * fg_41[k]
                  - f_67 * fg_109[k]
                  + f_68 * fg_116[k]
                  + f_69 * fg_139[k]
                  - f_70 * fg_146[k];
    }

#pragma omp simd aligned(fg_31, fg_36, fg_38, fg_106, fg_111, fg_113, fg_136, fg_141, \
                         fg_143 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = f_71 * fg_31[k]
                  + f_71 * fg_36[k]
                  - f_72 * fg_38[k]
                  + f_71 * fg_106[k]
                  + f_71 * fg_111[k]
                  - f_72 * fg_113[k]
                  - f_73 * fg_136[k]
                  - f_73 * fg_141[k]
                  + f_74 * fg_143[k];
    }

#pragma omp simd aligned(fg_34, fg_41, fg_43, fg_109, fg_116, fg_118, fg_139, fg_146, \
                         fg_148 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = f_75 * fg_34[k]
                  + f_75 * fg_41[k]
                  - f_76 * fg_43[k]
                  + f_75 * fg_109[k]
                  + f_75 * fg_116[k]
                  - f_76 * fg_118[k]
                  - f_12 * fg_139[k]
                  - f_12 * fg_146[k]
                  + f_77 * fg_148[k];
    }

#pragma omp simd aligned(fg_30, fg_33, fg_35, fg_40, fg_42, fg_44, fg_105, fg_108, fg_110, \
                         fg_115, fg_117, fg_119, fg_135, fg_138, fg_140, fg_145, fg_147, \
                         fg_149 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = -0.5625 * fg_30[k]
                  - 1.125 * fg_33[k]
                  + 4.5 * fg_35[k]
                  - 0.5625 * fg_40[k]
                  + 4.5 * fg_42[k]
                  - 1.5 * fg_44[k]
                  - 0.5625 * fg_105[k]
                  - 1.125 * fg_108[k]
                  + 4.5 * fg_110[k]
                  - 0.5625 * fg_115[k]
                  + 4.5 * fg_117[k]
                  - 1.5 * fg_119[k]
                  + 0.375 * fg_135[k]
                  + 0.75 * fg_138[k]
                  - 3.0 * fg_140[k]
                  + 0.375 * fg_145[k]
                  - 3.0 * fg_147[k]
                  + fg_149[k];
    }

#pragma omp simd aligned(fg_32, fg_37, fg_39, fg_107, fg_112, fg_114, fg_137, fg_142, \
                         fg_144 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = f_75 * fg_32[k]
                  + f_75 * fg_37[k]
                  - f_76 * fg_39[k]
                  + f_75 * fg_107[k]
                  + f_75 * fg_112[k]
                  - f_76 * fg_114[k]
                  - f_12 * fg_137[k]
                  - f_12 * fg_142[k]
                  + f_77 * fg_144[k];
    }

#pragma omp simd aligned(fg_30, fg_35, fg_40, fg_42, fg_105, fg_110, fg_115, fg_117, fg_135, \
                         fg_140, fg_145, fg_147 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = f_78 * fg_30[k]
                  - f_79 * fg_35[k]
                  - f_78 * fg_40[k]
                  + f_79 * fg_42[k]
                  + f_78 * fg_105[k]
                  - f_79 * fg_110[k]
                  - f_78 * fg_115[k]
                  + f_79 * fg_117[k]
                  - f_80 * fg_135[k]
                  + f_81 * fg_140[k]
                  + f_80 * fg_145[k]
                  - f_81 * fg_147[k];
    }

#pragma omp simd aligned(fg_32, fg_37, fg_107, fg_112, fg_137, fg_142 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = -f_68 * fg_32[k]
                  + f_67 * fg_37[k]
                  - f_68 * fg_107[k]
                  + f_67 * fg_112[k]
                  + f_70 * fg_137[k]
                  - f_69 * fg_142[k];
    }

#pragma omp simd aligned(fg_30, fg_33, fg_40, fg_105, fg_108, fg_115, fg_135, fg_138, \
                         fg_145 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = -f_82 * fg_30[k]
                  + f_83 * fg_33[k]
                  - f_82 * fg_40[k]
                  - f_82 * fg_105[k]
                  + f_83 * fg_108[k]
                  - f_82 * fg_115[k]
                  + f_84 * fg_135[k]
                  - f_65 * fg_138[k]
                  + f_84 * fg_145[k];
    }

#pragma omp simd aligned(fg_1, fg_4, fg_6, fg_11, fg_46, fg_49, fg_51, fg_56, fg_76, fg_79, \
                         fg_81, fg_86 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_36[k] = -f_38 * fg_1[k]
                  + f_38 * fg_6[k]
                  - f_38 * fg_46[k]
                  + f_38 * fg_51[k]
                  + f_39 * fg_76[k]
                  - f_39 * fg_81[k];

        g_37[k] = -f_40 * fg_4[k]
                  + f_41 * fg_11[k]
                  - f_40 * fg_49[k]
                  + f_41 * fg_56[k]
                  + f_42 * fg_79[k]
                  - f_43 * fg_86[k];
    }

#pragma omp simd aligned(fg_1, fg_6, fg_8, fg_46, fg_51, fg_53, fg_76, fg_81, \
                         fg_83 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = f_44 * fg_1[k]
                  + f_44 * fg_6[k]
                  - f_45 * fg_8[k]
                  + f_44 * fg_46[k]
                  + f_44 * fg_51[k]
                  - f_45 * fg_53[k]
                  - f_46 * fg_76[k]
                  - f_46 * fg_81[k]
                  + f_47 * fg_83[k];
    }

#pragma omp simd aligned(fg_4, fg_11, fg_13, fg_49, fg_56, fg_58, fg_79, fg_86, \
                         fg_88 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_39[k] = f_30 * fg_4[k]
                  + f_30 * fg_11[k]
                  - f_48 * fg_13[k]
                  + f_30 * fg_49[k]
                  + f_30 * fg_56[k]
                  - f_48 * fg_58[k]
                  - f_49 * fg_79[k]
                  - f_49 * fg_86[k]
                  + f_50 * fg_88[k];
    }

#pragma omp simd aligned(fg_0, fg_3, fg_5, fg_10, fg_12, fg_14, fg_45, fg_48, fg_50, fg_55, \
                         fg_57, fg_59, fg_75, fg_78, fg_80, fg_85, fg_87, \
                         fg_89 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_40[k] = -f_51 * fg_0[k]
                  - f_52 * fg_3[k]
                  + f_53 * fg_5[k]
                  - f_51 * fg_10[k]
                  + f_53 * fg_12[k]
                  - f_54 * fg_14[k]
                  - f_51 * fg_45[k]
                  - f_52 * fg_48[k]
                  + f_53 * fg_50[k]
                  - f_51 * fg_55[k]
                  + f_53 * fg_57[k]
                  - f_54 * fg_59[k]
                  + f_55 * fg_75[k]
                  + f_53 * fg_78[k]
                  - f_56 * fg_80[k]
                  + f_55 * fg_85[k]
                  - f_56 * fg_87[k]
                  + f_57 * fg_89[k];
    }

#pragma omp simd aligned(fg_2, fg_7, fg_9, fg_47, fg_52, fg_54, fg_77, fg_82, \
                         fg_84 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_41[k] = f_30 * fg_2[k]
                  + f_30 * fg_7[k]
                  - f_48 * fg_9[k]
                  + f_30 * fg_47[k]
                  + f_30 * fg_52[k]
                  - f_48 * fg_54[k]
                  - f_49 * fg_77[k]
                  - f_49 * fg_82[k]
                  + f_50 * fg_84[k];
    }

#pragma omp simd aligned(fg_0, fg_5, fg_10, fg_12, fg_45, fg_50, fg_55, fg_57, fg_75, fg_80, \
                         fg_85, fg_87 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = f_58 * fg_0[k]
                  - f_59 * fg_5[k]
                  - f_58 * fg_10[k]
                  + f_59 * fg_12[k]
                  + f_58 * fg_45[k]
                  - f_59 * fg_50[k]
                  - f_58 * fg_55[k]
                  + f_59 * fg_57[k]
                  - f_60 * fg_75[k]
                  + f_61 * fg_80[k]
                  + f_60 * fg_85[k]
                  - f_61 * fg_87[k];
    }

#pragma omp simd aligned(fg_2, fg_7, fg_47, fg_52, fg_77, fg_82 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_43[k] = -f_41 * fg_2[k]
                  + f_40 * fg_7[k]
                  - f_41 * fg_47[k]
                  + f_40 * fg_52[k]
                  + f_43 * fg_77[k]
                  - f_42 * fg_82[k];
    }

#pragma omp simd aligned(fg_0, fg_3, fg_10, fg_31, fg_36, fg_45, fg_48, fg_55, fg_75, fg_78, \
                         fg_85, fg_106, fg_111 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_44[k] = -f_62 * fg_0[k]
                  + f_63 * fg_3[k]
                  - f_62 * fg_10[k]
                  - f_62 * fg_45[k]
                  + f_63 * fg_48[k]
                  - f_62 * fg_55[k]
                  + f_38 * fg_75[k]
                  - f_64 * fg_78[k]
                  + f_38 * fg_85[k];

        g_45[k] = f_85 * fg_31[k]
                  - f_85 * fg_36[k]
                  - f_85 * fg_106[k]
                  + f_85 * fg_111[k];
    }

#pragma omp simd aligned(fg_31, fg_34, fg_36, fg_38, fg_41, fg_43, fg_106, fg_109, fg_111, \
                         fg_113, fg_116, fg_118 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_46[k] = f_86 * fg_34[k]
                  - f_87 * fg_41[k]
                  - f_86 * fg_109[k]
                  + f_87 * fg_116[k];

        g_47[k] = -f_34 * fg_31[k]
                  - f_34 * fg_36[k]
                  + f_35 * fg_38[k]
                  + f_34 * fg_106[k]
                  + f_34 * fg_111[k]
                  - f_35 * fg_113[k];

        g_48[k] = -f_88 * fg_34[k]
                  - f_88 * fg_41[k]
                  + f_89 * fg_43[k]
                  + f_88 * fg_109[k]
                  + f_88 * fg_116[k]
                  - f_89 * fg_118[k];
    }

#pragma omp simd aligned(fg_30, fg_33, fg_35, fg_40, fg_42, fg_44, fg_105, fg_108, fg_110, \
                         fg_115, fg_117, fg_119 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_49[k] = f_90 * fg_30[k]
                  + f_30 * fg_33[k]
                  - f_49 * fg_35[k]
                  + f_90 * fg_40[k]
                  - f_49 * fg_42[k]
                  + f_48 * fg_44[k]
                  - f_90 * fg_105[k]
                  - f_30 * fg_108[k]
                  + f_49 * fg_110[k]
                  - f_90 * fg_115[k]
                  + f_49 * fg_117[k]
                  - f_48 * fg_119[k];
    }

#pragma omp simd aligned(fg_30, fg_32, fg_35, fg_37, fg_39, fg_40, fg_42, fg_105, fg_107, \
                         fg_110, fg_112, fg_114, fg_115, fg_117 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_50[k] = -f_88 * fg_32[k]
                  - f_88 * fg_37[k]
                  + f_89 * fg_39[k]
                  + f_88 * fg_107[k]
                  + f_88 * fg_112[k]
                  - f_89 * fg_114[k];

        g_51[k] = -f_91 * fg_30[k]
                  + f_92 * fg_35[k]
                  + f_91 * fg_40[k]
                  - f_92 * fg_42[k]
                  + f_91 * fg_105[k]
                  - f_92 * fg_110[k]
                  - f_91 * fg_115[k]
                  + f_92 * fg_117[k];
    }

#pragma omp simd aligned(fg_30, fg_32, fg_33, fg_37, fg_40, fg_105, fg_107, fg_108, fg_112, \
                         fg_115 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_52[k] = f_87 * fg_32[k]
                  - f_86 * fg_37[k]
                  - f_87 * fg_107[k]
                  + f_86 * fg_112[k];

        g_53[k] = f_93 * fg_30[k]
                  - f_94 * fg_33[k]
                  + f_93 * fg_40[k]
                  - f_93 * fg_105[k]
                  + f_94 * fg_108[k]
                  - f_93 * fg_115[k];
    }

#pragma omp simd aligned(fg_1, fg_4, fg_6, fg_8, fg_11, fg_13, fg_46, fg_49, fg_51, fg_53, \
                         fg_56, fg_58 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_54[k] = f_1 * fg_1[k]
                  - f_1 * fg_6[k]
                  - f_0 * fg_46[k]
                  + f_0 * fg_51[k];

        g_55[k] = f_3 * fg_4[k]
                  - f_4 * fg_11[k]
                  - f_2 * fg_49[k]
                  + f_3 * fg_56[k];

        g_56[k] = -f_7 * fg_1[k]
                  - f_7 * fg_6[k]
                  + f_8 * fg_8[k]
                  + f_5 * fg_46[k]
                  + f_5 * fg_51[k]
                  - f_6 * fg_53[k];

        g_57[k] = -1.875 * fg_4[k]
                  - 1.875 * fg_11[k]
                  + 2.5 * fg_13[k]
                  + 5.625 * fg_49[k]
                  + 5.625 * fg_56[k]
                  - 7.5 * fg_58[k];
    }

#pragma omp simd aligned(fg_0, fg_3, fg_5, fg_10, fg_12, fg_14, fg_45, fg_48, fg_50, fg_55, \
                         fg_57, fg_59 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_58[k] = f_13 * fg_0[k]
                  + f_14 * fg_3[k]
                  - f_12 * fg_5[k]
                  + f_13 * fg_10[k]
                  - f_12 * fg_12[k]
                  + f_15 * fg_14[k]
                  - f_9 * fg_45[k]
                  - f_10 * fg_48[k]
                  + f_11 * fg_50[k]
                  - f_9 * fg_55[k]
                  + f_11 * fg_57[k]
                  - f_12 * fg_59[k];
    }

#pragma omp simd aligned(fg_0, fg_2, fg_5, fg_7, fg_9, fg_10, fg_12, fg_45, fg_47, fg_50, \
                         fg_52, fg_54, fg_55, fg_57 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_59[k] = -1.875 * fg_2[k]
                  - 1.875 * fg_7[k]
                  + 2.5 * fg_9[k]
                  + 5.625 * fg_47[k]
                  + 5.625 * fg_52[k]
                  - 7.5 * fg_54[k];

        g_60[k] = -f_18 * fg_0[k]
                  + f_5 * fg_5[k]
                  + f_18 * fg_10[k]
                  - f_5 * fg_12[k]
                  + f_16 * fg_45[k]
                  - f_17 * fg_50[k]
                  - f_16 * fg_55[k]
                  + f_17 * fg_57[k];
    }

#pragma omp simd aligned(fg_0, fg_2, fg_3, fg_7, fg_10, fg_45, fg_47, fg_48, fg_52, \
                         fg_55 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_61[k] = f_4 * fg_2[k]
                  - f_3 * fg_7[k]
                  - f_3 * fg_47[k]
                  + f_2 * fg_52[k];

        g_62[k] = f_21 * fg_0[k]
                  - f_22 * fg_3[k]
                  + f_21 * fg_10[k]
                  - f_19 * fg_45[k]
                  + f_20 * fg_48[k]
                  - f_19 * fg_55[k];
    }
}

}  // namespace simdtrf
