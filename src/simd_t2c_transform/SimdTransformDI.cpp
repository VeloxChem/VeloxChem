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


#include "SimdTransformDI.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_di(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t di,
             const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 0.5625 * std::sqrt(154.0);
    const auto f_1 = 1.875 * std::sqrt(154.0);
    const auto f_2 = 0.9375 * std::sqrt(462.0);
    const auto f_3 = 1.875 * std::sqrt(462.0);
    const auto f_4 = 0.1875 * std::sqrt(462.0);
    const auto f_5 = 0.75 * std::sqrt(21.0);
    const auto f_6 = 7.5 * std::sqrt(21.0);
    const auto f_7 = 1.6875 * std::sqrt(70.0);
    const auto f_8 = 1.125 * std::sqrt(70.0);
    const auto f_9 = 4.5 * std::sqrt(70.0);
    const auto f_10 = 0.5625 * std::sqrt(70.0);
    const auto f_11 = 1.5 * std::sqrt(70.0);
    const auto f_12 = 0.1875 * std::sqrt(70.0);
    const auto f_13 = 0.375 * std::sqrt(70.0);
    const auto f_14 = 3.0 * std::sqrt(70.0);
    const auto f_15 = 1.875 * std::sqrt(7.0);
    const auto f_16 = 3.75 * std::sqrt(7.0);
    const auto f_17 = 7.5 * std::sqrt(7.0);
    const auto f_18 = 3.0 * std::sqrt(7.0);
    const auto f_19 = 0.3125 * std::sqrt(3.0);
    const auto f_20 = 0.9375 * std::sqrt(3.0);
    const auto f_21 = 5.625 * std::sqrt(3.0);
    const auto f_22 = 11.25 * std::sqrt(3.0);
    const auto f_23 = 7.5 * std::sqrt(3.0);
    const auto f_24 = std::sqrt(3.0);
    const auto f_25 = 0.09375 * std::sqrt(70.0);
    const auto f_26 = 0.1875 * std::sqrt(21.0);
    const auto f_27 = 0.9375 * std::sqrt(21.0);
    const auto f_28 = 1.875 * std::sqrt(21.0);
    const auto f_29 = 11.25 * std::sqrt(21.0);
    const auto f_30 = 0.09375 * std::sqrt(154.0);
    const auto f_31 = 1.40625 * std::sqrt(154.0);
    const auto f_32 = 0.09375 * std::sqrt(462.0);
    const auto f_33 = 0.3125 * std::sqrt(462.0);
    const auto f_34 = 0.625 * std::sqrt(462.0);
    const auto f_35 = 0.46875 * std::sqrt(154.0);
    const auto f_36 = 0.9375 * std::sqrt(154.0);
    const auto f_37 = 0.1875 * std::sqrt(154.0);
    const auto f_38 = 0.375 * std::sqrt(7.0);
    const auto f_39 = 0.75 * std::sqrt(7.0);
    const auto f_40 = 0.28125 * std::sqrt(210.0);
    const auto f_41 = 0.1875 * std::sqrt(210.0);
    const auto f_42 = 0.75 * std::sqrt(210.0);
    const auto f_43 = 0.09375 * std::sqrt(210.0);
    const auto f_44 = 0.25 * std::sqrt(210.0);
    const auto f_45 = 0.5625 * std::sqrt(210.0);
    const auto f_46 = 0.375 * std::sqrt(210.0);
    const auto f_47 = 1.5 * std::sqrt(210.0);
    const auto f_48 = 0.5 * std::sqrt(210.0);
    const auto f_49 = 0.03125 * std::sqrt(210.0);
    const auto f_50 = 0.0625 * std::sqrt(210.0);
    const auto f_51 = 0.125 * std::sqrt(210.0);
    const auto f_52 = std::sqrt(210.0);
    const auto f_53 = 0.3125 * std::sqrt(21.0);
    const auto f_54 = 0.625 * std::sqrt(21.0);
    const auto f_55 = 1.25 * std::sqrt(21.0);
    const auto f_56 = 0.5 * std::sqrt(21.0);
    const auto f_57 = 2.5 * std::sqrt(21.0);
    const auto f_58 = std::sqrt(21.0);
    const auto f_59 = 0.015625 * std::sqrt(210.0);
    const auto f_60 = 0.09375 * std::sqrt(7.0);
    const auto f_61 = 0.46875 * std::sqrt(7.0);
    const auto f_62 = 0.9375 * std::sqrt(7.0);
    const auto f_63 = 5.625 * std::sqrt(7.0);
    const auto f_64 = 0.1875 * std::sqrt(7.0);
    const auto f_65 = 11.25 * std::sqrt(7.0);
    const auto f_66 = 0.015625 * std::sqrt(462.0);
    const auto f_67 = 0.234375 * std::sqrt(462.0);
    const auto f_68 = 0.03125 * std::sqrt(462.0);
    const auto f_69 = 0.46875 * std::sqrt(462.0);
    const auto f_70 = 0.28125 * std::sqrt(154.0);
    const auto f_71 = 0.375 * std::sqrt(21.0);
    const auto f_72 = 3.75 * std::sqrt(21.0);
    const auto f_73 = 0.84375 * std::sqrt(70.0);
    const auto f_74 = 2.25 * std::sqrt(70.0);
    const auto f_75 = 0.28125 * std::sqrt(70.0);
    const auto f_76 = 0.75 * std::sqrt(70.0);
    const auto f_77 = 1.5 * std::sqrt(7.0);
    const auto f_78 = 0.15625 * std::sqrt(3.0);
    const auto f_79 = 0.46875 * std::sqrt(3.0);
    const auto f_80 = 2.8125 * std::sqrt(3.0);
    const auto f_81 = 3.75 * std::sqrt(3.0);
    const auto f_82 = 0.5 * std::sqrt(3.0);
    const auto f_83 = 0.046875 * std::sqrt(70.0);
    const auto f_84 = 0.09375 * std::sqrt(21.0);
    const auto f_85 = 0.46875 * std::sqrt(21.0);
    const auto f_86 = 5.625 * std::sqrt(21.0);
    const auto f_87 = 0.046875 * std::sqrt(154.0);
    const auto f_88 = 0.703125 * std::sqrt(154.0);

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
    auto *g_63 = values + 63 * nvalues;
    auto *g_64 = values + 64 * nvalues;

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_1 = buffer.data(di + 1);
    const auto *di_2 = buffer.data(di + 2);
    const auto *di_3 = buffer.data(di + 3);
    const auto *di_4 = buffer.data(di + 4);
    const auto *di_5 = buffer.data(di + 5);
    const auto *di_6 = buffer.data(di + 6);
    const auto *di_7 = buffer.data(di + 7);
    const auto *di_8 = buffer.data(di + 8);
    const auto *di_9 = buffer.data(di + 9);
    const auto *di_10 = buffer.data(di + 10);
    const auto *di_11 = buffer.data(di + 11);
    const auto *di_12 = buffer.data(di + 12);
    const auto *di_13 = buffer.data(di + 13);
    const auto *di_14 = buffer.data(di + 14);
    const auto *di_15 = buffer.data(di + 15);
    const auto *di_16 = buffer.data(di + 16);
    const auto *di_17 = buffer.data(di + 17);
    const auto *di_18 = buffer.data(di + 18);
    const auto *di_19 = buffer.data(di + 19);
    const auto *di_20 = buffer.data(di + 20);
    const auto *di_21 = buffer.data(di + 21);
    const auto *di_22 = buffer.data(di + 22);
    const auto *di_23 = buffer.data(di + 23);
    const auto *di_24 = buffer.data(di + 24);
    const auto *di_25 = buffer.data(di + 25);
    const auto *di_26 = buffer.data(di + 26);
    const auto *di_27 = buffer.data(di + 27);
    const auto *di_28 = buffer.data(di + 28);
    const auto *di_29 = buffer.data(di + 29);
    const auto *di_30 = buffer.data(di + 30);
    const auto *di_31 = buffer.data(di + 31);
    const auto *di_32 = buffer.data(di + 32);
    const auto *di_33 = buffer.data(di + 33);
    const auto *di_34 = buffer.data(di + 34);
    const auto *di_35 = buffer.data(di + 35);
    const auto *di_36 = buffer.data(di + 36);
    const auto *di_37 = buffer.data(di + 37);
    const auto *di_38 = buffer.data(di + 38);
    const auto *di_39 = buffer.data(di + 39);
    const auto *di_40 = buffer.data(di + 40);
    const auto *di_41 = buffer.data(di + 41);
    const auto *di_42 = buffer.data(di + 42);
    const auto *di_43 = buffer.data(di + 43);
    const auto *di_44 = buffer.data(di + 44);
    const auto *di_45 = buffer.data(di + 45);
    const auto *di_46 = buffer.data(di + 46);
    const auto *di_47 = buffer.data(di + 47);
    const auto *di_48 = buffer.data(di + 48);
    const auto *di_49 = buffer.data(di + 49);
    const auto *di_50 = buffer.data(di + 50);
    const auto *di_51 = buffer.data(di + 51);
    const auto *di_52 = buffer.data(di + 52);
    const auto *di_53 = buffer.data(di + 53);
    const auto *di_54 = buffer.data(di + 54);
    const auto *di_55 = buffer.data(di + 55);
    const auto *di_56 = buffer.data(di + 56);
    const auto *di_57 = buffer.data(di + 57);
    const auto *di_58 = buffer.data(di + 58);
    const auto *di_59 = buffer.data(di + 59);
    const auto *di_60 = buffer.data(di + 60);
    const auto *di_61 = buffer.data(di + 61);
    const auto *di_62 = buffer.data(di + 62);
    const auto *di_63 = buffer.data(di + 63);
    const auto *di_64 = buffer.data(di + 64);
    const auto *di_65 = buffer.data(di + 65);
    const auto *di_66 = buffer.data(di + 66);
    const auto *di_67 = buffer.data(di + 67);
    const auto *di_68 = buffer.data(di + 68);
    const auto *di_69 = buffer.data(di + 69);
    const auto *di_70 = buffer.data(di + 70);
    const auto *di_71 = buffer.data(di + 71);
    const auto *di_72 = buffer.data(di + 72);
    const auto *di_73 = buffer.data(di + 73);
    const auto *di_74 = buffer.data(di + 74);
    const auto *di_75 = buffer.data(di + 75);
    const auto *di_76 = buffer.data(di + 76);
    const auto *di_77 = buffer.data(di + 77);
    const auto *di_78 = buffer.data(di + 78);
    const auto *di_79 = buffer.data(di + 79);
    const auto *di_80 = buffer.data(di + 80);
    const auto *di_81 = buffer.data(di + 81);
    const auto *di_82 = buffer.data(di + 82);
    const auto *di_83 = buffer.data(di + 83);
    const auto *di_84 = buffer.data(di + 84);
    const auto *di_85 = buffer.data(di + 85);
    const auto *di_86 = buffer.data(di + 86);
    const auto *di_87 = buffer.data(di + 87);
    const auto *di_88 = buffer.data(di + 88);
    const auto *di_89 = buffer.data(di + 89);
    const auto *di_90 = buffer.data(di + 90);
    const auto *di_91 = buffer.data(di + 91);
    const auto *di_92 = buffer.data(di + 92);
    const auto *di_93 = buffer.data(di + 93);
    const auto *di_94 = buffer.data(di + 94);
    const auto *di_95 = buffer.data(di + 95);
    const auto *di_96 = buffer.data(di + 96);
    const auto *di_97 = buffer.data(di + 97);
    const auto *di_98 = buffer.data(di + 98);
    const auto *di_99 = buffer.data(di + 99);
    const auto *di_100 = buffer.data(di + 100);
    const auto *di_101 = buffer.data(di + 101);
    const auto *di_102 = buffer.data(di + 102);
    const auto *di_103 = buffer.data(di + 103);
    const auto *di_104 = buffer.data(di + 104);
    const auto *di_105 = buffer.data(di + 105);
    const auto *di_106 = buffer.data(di + 106);
    const auto *di_107 = buffer.data(di + 107);
    const auto *di_108 = buffer.data(di + 108);
    const auto *di_109 = buffer.data(di + 109);
    const auto *di_110 = buffer.data(di + 110);
    const auto *di_111 = buffer.data(di + 111);
    const auto *di_112 = buffer.data(di + 112);
    const auto *di_113 = buffer.data(di + 113);
    const auto *di_114 = buffer.data(di + 114);
    const auto *di_115 = buffer.data(di + 115);
    const auto *di_116 = buffer.data(di + 116);
    const auto *di_117 = buffer.data(di + 117);
    const auto *di_118 = buffer.data(di + 118);
    const auto *di_119 = buffer.data(di + 119);
    const auto *di_120 = buffer.data(di + 120);
    const auto *di_121 = buffer.data(di + 121);
    const auto *di_122 = buffer.data(di + 122);
    const auto *di_123 = buffer.data(di + 123);
    const auto *di_124 = buffer.data(di + 124);
    const auto *di_125 = buffer.data(di + 125);
    const auto *di_126 = buffer.data(di + 126);
    const auto *di_127 = buffer.data(di + 127);
    const auto *di_128 = buffer.data(di + 128);
    const auto *di_129 = buffer.data(di + 129);
    const auto *di_130 = buffer.data(di + 130);
    const auto *di_131 = buffer.data(di + 131);
    const auto *di_132 = buffer.data(di + 132);
    const auto *di_133 = buffer.data(di + 133);
    const auto *di_134 = buffer.data(di + 134);
    const auto *di_135 = buffer.data(di + 135);
    const auto *di_136 = buffer.data(di + 136);
    const auto *di_137 = buffer.data(di + 137);
    const auto *di_138 = buffer.data(di + 138);
    const auto *di_139 = buffer.data(di + 139);
    const auto *di_140 = buffer.data(di + 140);
    const auto *di_141 = buffer.data(di + 141);
    const auto *di_142 = buffer.data(di + 142);
    const auto *di_143 = buffer.data(di + 143);
    const auto *di_144 = buffer.data(di + 144);
    const auto *di_145 = buffer.data(di + 145);
    const auto *di_146 = buffer.data(di + 146);
    const auto *di_147 = buffer.data(di + 147);
    const auto *di_148 = buffer.data(di + 148);
    const auto *di_149 = buffer.data(di + 149);
    const auto *di_150 = buffer.data(di + 150);
    const auto *di_151 = buffer.data(di + 151);
    const auto *di_152 = buffer.data(di + 152);
    const auto *di_153 = buffer.data(di + 153);
    const auto *di_154 = buffer.data(di + 154);
    const auto *di_155 = buffer.data(di + 155);
    const auto *di_156 = buffer.data(di + 156);
    const auto *di_157 = buffer.data(di + 157);
    const auto *di_158 = buffer.data(di + 158);
    const auto *di_159 = buffer.data(di + 159);
    const auto *di_160 = buffer.data(di + 160);
    const auto *di_161 = buffer.data(di + 161);
    const auto *di_162 = buffer.data(di + 162);
    const auto *di_163 = buffer.data(di + 163);
    const auto *di_164 = buffer.data(di + 164);
    const auto *di_165 = buffer.data(di + 165);
    const auto *di_166 = buffer.data(di + 166);
    const auto *di_167 = buffer.data(di + 167);

#pragma omp simd aligned(di_29, di_32, di_34, di_36, di_39, di_41, di_43, di_45, di_47, di_50, \
                         di_52 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * di_29[k]
                 - f_1 * di_34[k]
                 + f_0 * di_43[k];

        g_1[k] = f_2 * di_32[k]
                 - f_3 * di_39[k]
                 + f_4 * di_50[k];

        g_2[k] = -f_5 * di_29[k]
                 + f_6 * di_36[k]
                 + f_5 * di_43[k]
                 - f_6 * di_45[k];

        g_3[k] = -f_7 * di_32[k]
                 - f_8 * di_39[k]
                 + f_9 * di_41[k]
                 + f_10 * di_50[k]
                 - f_11 * di_52[k];

        g_4[k] = f_12 * di_29[k]
                 + f_13 * di_34[k]
                 - f_14 * di_36[k]
                 + f_12 * di_43[k]
                 - f_14 * di_45[k]
                 + f_14 * di_47[k];
    }

#pragma omp simd aligned(di_32, di_39, di_41, di_50, di_52, di_54 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = f_15 * di_32[k]
                 + f_16 * di_39[k]
                 - f_17 * di_41[k]
                 + f_15 * di_50[k]
                 - f_17 * di_52[k]
                 + f_18 * di_54[k];
    }

#pragma omp simd aligned(di_28, di_31, di_33, di_38, di_40, di_42, di_49, di_51, di_53, \
                         di_55 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = -f_19 * di_28[k]
                 - f_20 * di_31[k]
                 + f_21 * di_33[k]
                 - f_20 * di_38[k]
                 + f_22 * di_40[k]
                 - f_23 * di_42[k]
                 - f_19 * di_49[k]
                 + f_21 * di_51[k]
                 - f_23 * di_53[k]
                 + f_24 * di_55[k];
    }

#pragma omp simd aligned(di_28, di_30, di_31, di_33, di_35, di_37, di_38, di_42, di_44, di_46, \
                         di_48, di_49, di_51, di_53 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = f_15 * di_30[k]
                 + f_16 * di_35[k]
                 - f_17 * di_37[k]
                 + f_15 * di_44[k]
                 - f_17 * di_46[k]
                 + f_18 * di_48[k];

        g_8[k] = f_25 * di_28[k]
                 + f_25 * di_31[k]
                 - f_11 * di_33[k]
                 - f_25 * di_38[k]
                 + f_11 * di_42[k]
                 - f_25 * di_49[k]
                 + f_11 * di_51[k]
                 - f_11 * di_53[k];
    }

#pragma omp simd aligned(di_28, di_30, di_31, di_33, di_35, di_37, di_38, di_40, di_44, di_46, \
                         di_49, di_51 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_9[k] = -f_10 * di_30[k]
                 + f_8 * di_35[k]
                 + f_11 * di_37[k]
                 + f_7 * di_44[k]
                 - f_9 * di_46[k];

        g_10[k] = -f_26 * di_28[k]
                  + f_27 * di_31[k]
                  + f_28 * di_33[k]
                  + f_27 * di_38[k]
                  - f_29 * di_40[k]
                  - f_26 * di_49[k]
                  + f_28 * di_51[k];

        g_11[k] = f_4 * di_30[k]
                  - f_3 * di_35[k]
                  + f_2 * di_44[k];

        g_12[k] = f_30 * di_28[k]
                  - f_31 * di_31[k]
                  + f_31 * di_38[k]
                  - f_30 * di_49[k];
    }

#pragma omp simd aligned(di_113, di_116, di_118, di_120, di_123, di_125, di_127, di_129, \
                         di_131, di_134, di_136 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = f_0 * di_113[k]
                  - f_1 * di_118[k]
                  + f_0 * di_127[k];

        g_14[k] = f_2 * di_116[k]
                  - f_3 * di_123[k]
                  + f_4 * di_134[k];

        g_15[k] = -f_5 * di_113[k]
                  + f_6 * di_120[k]
                  + f_5 * di_127[k]
                  - f_6 * di_129[k];

        g_16[k] = -f_7 * di_116[k]
                  - f_8 * di_123[k]
                  + f_9 * di_125[k]
                  + f_10 * di_134[k]
                  - f_11 * di_136[k];

        g_17[k] = f_12 * di_113[k]
                  + f_13 * di_118[k]
                  - f_14 * di_120[k]
                  + f_12 * di_127[k]
                  - f_14 * di_129[k]
                  + f_14 * di_131[k];
    }

#pragma omp simd aligned(di_116, di_123, di_125, di_134, di_136, \
                         di_138 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = f_15 * di_116[k]
                  + f_16 * di_123[k]
                  - f_17 * di_125[k]
                  + f_15 * di_134[k]
                  - f_17 * di_136[k]
                  + f_18 * di_138[k];
    }

#pragma omp simd aligned(di_112, di_115, di_117, di_122, di_124, di_126, di_133, di_135, \
                         di_137, di_139 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = -f_19 * di_112[k]
                  - f_20 * di_115[k]
                  + f_21 * di_117[k]
                  - f_20 * di_122[k]
                  + f_22 * di_124[k]
                  - f_23 * di_126[k]
                  - f_19 * di_133[k]
                  + f_21 * di_135[k]
                  - f_23 * di_137[k]
                  + f_24 * di_139[k];
    }

#pragma omp simd aligned(di_112, di_114, di_115, di_117, di_119, di_121, di_122, di_126, \
                         di_128, di_130, di_132, di_133, di_135, \
                         di_137 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = f_15 * di_114[k]
                  + f_16 * di_119[k]
                  - f_17 * di_121[k]
                  + f_15 * di_128[k]
                  - f_17 * di_130[k]
                  + f_18 * di_132[k];

        g_21[k] = f_25 * di_112[k]
                  + f_25 * di_115[k]
                  - f_11 * di_117[k]
                  - f_25 * di_122[k]
                  + f_11 * di_126[k]
                  - f_25 * di_133[k]
                  + f_11 * di_135[k]
                  - f_11 * di_137[k];
    }

#pragma omp simd aligned(di_112, di_114, di_115, di_117, di_119, di_121, di_122, di_124, \
                         di_128, di_130, di_133, di_135 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = -f_10 * di_114[k]
                  + f_8 * di_119[k]
                  + f_11 * di_121[k]
                  + f_7 * di_128[k]
                  - f_9 * di_130[k];

        g_23[k] = -f_26 * di_112[k]
                  + f_27 * di_115[k]
                  + f_28 * di_117[k]
                  + f_27 * di_122[k]
                  - f_29 * di_124[k]
                  - f_26 * di_133[k]
                  + f_28 * di_135[k];

        g_24[k] = f_4 * di_114[k]
                  - f_3 * di_119[k]
                  + f_2 * di_128[k];

        g_25[k] = f_30 * di_112[k]
                  - f_31 * di_115[k]
                  + f_31 * di_122[k]
                  - f_30 * di_133[k];
    }

#pragma omp simd aligned(di_1, di_6, di_15, di_85, di_90, di_99, di_141, di_146, \
                         di_155 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = -f_32 * di_1[k]
                  + f_33 * di_6[k]
                  - f_32 * di_15[k]
                  - f_32 * di_85[k]
                  + f_33 * di_90[k]
                  - f_32 * di_99[k]
                  + f_4 * di_141[k]
                  - f_34 * di_146[k]
                  + f_4 * di_155[k];
    }

#pragma omp simd aligned(di_4, di_11, di_22, di_88, di_95, di_106, di_144, di_151, \
                         di_162 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = -f_35 * di_4[k]
                  + f_36 * di_11[k]
                  - f_30 * di_22[k]
                  - f_35 * di_88[k]
                  + f_36 * di_95[k]
                  - f_30 * di_106[k]
                  + f_36 * di_144[k]
                  - f_1 * di_151[k]
                  + f_37 * di_162[k];
    }

#pragma omp simd aligned(di_1, di_8, di_15, di_17, di_85, di_92, di_99, di_101, di_141, \
                         di_148, di_155, di_157 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = f_38 * di_1[k]
                  - f_16 * di_8[k]
                  - f_38 * di_15[k]
                  + f_16 * di_17[k]
                  + f_38 * di_85[k]
                  - f_16 * di_92[k]
                  - f_38 * di_99[k]
                  + f_16 * di_101[k]
                  - f_39 * di_141[k]
                  + f_17 * di_148[k]
                  + f_39 * di_155[k]
                  - f_17 * di_157[k];
    }

#pragma omp simd aligned(di_4, di_11, di_13, di_22, di_24, di_88, di_95, di_97, di_106, \
                         di_108, di_144, di_151, di_153, di_162, \
                         di_164 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = f_40 * di_4[k]
                  + f_41 * di_11[k]
                  - f_42 * di_13[k]
                  - f_43 * di_22[k]
                  + f_44 * di_24[k]
                  + f_40 * di_88[k]
                  + f_41 * di_95[k]
                  - f_42 * di_97[k]
                  - f_43 * di_106[k]
                  + f_44 * di_108[k]
                  - f_45 * di_144[k]
                  - f_46 * di_151[k]
                  + f_47 * di_153[k]
                  + f_41 * di_162[k]
                  - f_48 * di_164[k];
    }

#pragma omp simd aligned(di_1, di_6, di_8, di_15, di_17, di_19, di_85, di_90, di_92, di_99, \
                         di_101, di_103, di_141, di_146, di_148, di_155, di_157, \
                         di_159 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = -f_49 * di_1[k]
                  - f_50 * di_6[k]
                  + f_48 * di_8[k]
                  - f_49 * di_15[k]
                  + f_48 * di_17[k]
                  - f_48 * di_19[k]
                  - f_49 * di_85[k]
                  - f_50 * di_90[k]
                  + f_48 * di_92[k]
                  - f_49 * di_99[k]
                  + f_48 * di_101[k]
                  - f_48 * di_103[k]
                  + f_50 * di_141[k]
                  + f_51 * di_146[k]
                  - f_52 * di_148[k]
                  + f_50 * di_155[k]
                  - f_52 * di_157[k]
                  + f_52 * di_159[k];
    }

#pragma omp simd aligned(di_4, di_11, di_13, di_22, di_24, di_26, di_88, di_95, di_97, di_106, \
                         di_108, di_110, di_144, di_151, di_153, di_162, di_164, \
                         di_166 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = -f_53 * di_4[k]
                  - f_54 * di_11[k]
                  + f_55 * di_13[k]
                  - f_53 * di_22[k]
                  + f_55 * di_24[k]
                  - f_56 * di_26[k]
                  - f_53 * di_88[k]
                  - f_54 * di_95[k]
                  + f_55 * di_97[k]
                  - f_53 * di_106[k]
                  + f_55 * di_108[k]
                  - f_56 * di_110[k]
                  + f_54 * di_144[k]
                  + f_55 * di_151[k]
                  - f_57 * di_153[k]
                  + f_54 * di_162[k]
                  - f_57 * di_164[k]
                  + f_58 * di_166[k];
    }

#pragma omp simd aligned(di_0, di_3, di_5, di_10, di_12, di_14, di_21, di_23, di_25, di_27, \
                         di_84, di_87, di_89, di_94, di_96, di_98, di_105, di_107, di_109, \
                         di_111, di_140, di_143, di_145, di_150, di_152, di_154, di_161, \
                         di_163, di_165, di_167 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = 0.15625 * di_0[k]
                  + 0.46875 * di_3[k]
                  - 2.8125 * di_5[k]
                  + 0.46875 * di_10[k]
                  - 5.625 * di_12[k]
                  + 3.75 * di_14[k]
                  + 0.15625 * di_21[k]
                  - 2.8125 * di_23[k]
                  + 3.75 * di_25[k]
                  - 0.5 * di_27[k]
                  + 0.15625 * di_84[k]
                  + 0.46875 * di_87[k]
                  - 2.8125 * di_89[k]
                  + 0.46875 * di_94[k]
                  - 5.625 * di_96[k]
                  + 3.75 * di_98[k]
                  + 0.15625 * di_105[k]
                  - 2.8125 * di_107[k]
                  + 3.75 * di_109[k]
                  - 0.5 * di_111[k]
                  - 0.3125 * di_140[k]
                  - 0.9375 * di_143[k]
                  + 5.625 * di_145[k]
                  - 0.9375 * di_150[k]
                  + 11.25 * di_152[k]
                  - 7.5 * di_154[k]
                  - 0.3125 * di_161[k]
                  + 5.625 * di_163[k]
                  - 7.5 * di_165[k]
                  + di_167[k];
    }

#pragma omp simd aligned(di_2, di_7, di_9, di_16, di_18, di_20, di_86, di_91, di_93, di_100, \
                         di_102, di_104, di_142, di_147, di_149, di_156, di_158, \
                         di_160 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = -f_53 * di_2[k]
                  - f_54 * di_7[k]
                  + f_55 * di_9[k]
                  - f_53 * di_16[k]
                  + f_55 * di_18[k]
                  - f_56 * di_20[k]
                  - f_53 * di_86[k]
                  - f_54 * di_91[k]
                  + f_55 * di_93[k]
                  - f_53 * di_100[k]
                  + f_55 * di_102[k]
                  - f_56 * di_104[k]
                  + f_54 * di_142[k]
                  + f_55 * di_147[k]
                  - f_57 * di_149[k]
                  + f_54 * di_156[k]
                  - f_57 * di_158[k]
                  + f_58 * di_160[k];
    }

#pragma omp simd aligned(di_0, di_3, di_5, di_10, di_14, di_21, di_23, di_25, di_84, di_87, \
                         di_89, di_94, di_98, di_105, di_107, di_109, di_140, di_143, di_145, \
                         di_150, di_154, di_161, di_163, di_165 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = -f_59 * di_0[k]
                  - f_59 * di_3[k]
                  + f_44 * di_5[k]
                  + f_59 * di_10[k]
                  - f_44 * di_14[k]
                  + f_59 * di_21[k]
                  - f_44 * di_23[k]
                  + f_44 * di_25[k]
                  - f_59 * di_84[k]
                  - f_59 * di_87[k]
                  + f_44 * di_89[k]
                  + f_59 * di_94[k]
                  - f_44 * di_98[k]
                  + f_59 * di_105[k]
                  - f_44 * di_107[k]
                  + f_44 * di_109[k]
                  + f_49 * di_140[k]
                  + f_49 * di_143[k]
                  - f_48 * di_145[k]
                  - f_49 * di_150[k]
                  + f_48 * di_154[k]
                  - f_49 * di_161[k]
                  + f_48 * di_163[k]
                  - f_48 * di_165[k];
    }

#pragma omp simd aligned(di_2, di_7, di_9, di_16, di_18, di_86, di_91, di_93, di_100, di_102, \
                         di_142, di_147, di_149, di_156, di_158 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = f_43 * di_2[k]
                  - f_41 * di_7[k]
                  - f_44 * di_9[k]
                  - f_40 * di_16[k]
                  + f_42 * di_18[k]
                  + f_43 * di_86[k]
                  - f_41 * di_91[k]
                  - f_44 * di_93[k]
                  - f_40 * di_100[k]
                  + f_42 * di_102[k]
                  - f_41 * di_142[k]
                  + f_46 * di_147[k]
                  + f_48 * di_149[k]
                  + f_45 * di_156[k]
                  - f_47 * di_158[k];
    }

#pragma omp simd aligned(di_0, di_3, di_5, di_10, di_12, di_21, di_23, di_84, di_87, di_89, \
                         di_94, di_96, di_105, di_107, di_140, di_143, di_145, di_150, di_152, \
                         di_161, di_163 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_36[k] = f_60 * di_0[k]
                  - f_61 * di_3[k]
                  - f_62 * di_5[k]
                  - f_61 * di_10[k]
                  + f_63 * di_12[k]
                  + f_60 * di_21[k]
                  - f_62 * di_23[k]
                  + f_60 * di_84[k]
                  - f_61 * di_87[k]
                  - f_62 * di_89[k]
                  - f_61 * di_94[k]
                  + f_63 * di_96[k]
                  + f_60 * di_105[k]
                  - f_62 * di_107[k]
                  - f_64 * di_140[k]
                  + f_62 * di_143[k]
                  + f_15 * di_145[k]
                  + f_62 * di_150[k]
                  - f_65 * di_152[k]
                  - f_64 * di_161[k]
                  + f_15 * di_163[k];
    }

#pragma omp simd aligned(di_2, di_7, di_16, di_86, di_91, di_100, di_142, di_147, \
                         di_156 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_37[k] = -f_30 * di_2[k]
                  + f_36 * di_7[k]
                  - f_35 * di_16[k]
                  - f_30 * di_86[k]
                  + f_36 * di_91[k]
                  - f_35 * di_100[k]
                  + f_37 * di_142[k]
                  - f_1 * di_147[k]
                  + f_36 * di_156[k];
    }

#pragma omp simd aligned(di_0, di_3, di_10, di_21, di_84, di_87, di_94, di_105, di_140, \
                         di_143, di_150, di_161 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = -f_66 * di_0[k]
                  + f_67 * di_3[k]
                  - f_67 * di_10[k]
                  + f_66 * di_21[k]
                  - f_66 * di_84[k]
                  + f_67 * di_87[k]
                  - f_67 * di_94[k]
                  + f_66 * di_105[k]
                  + f_68 * di_140[k]
                  - f_69 * di_143[k]
                  + f_69 * di_150[k]
                  - f_68 * di_161[k];
    }

#pragma omp simd aligned(di_57, di_60, di_62, di_64, di_67, di_69, di_71, di_73, di_75, di_78, \
                         di_80 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_39[k] = f_0 * di_57[k]
                  - f_1 * di_62[k]
                  + f_0 * di_71[k];

        g_40[k] = f_2 * di_60[k]
                  - f_3 * di_67[k]
                  + f_4 * di_78[k];

        g_41[k] = -f_5 * di_57[k]
                  + f_6 * di_64[k]
                  + f_5 * di_71[k]
                  - f_6 * di_73[k];

        g_42[k] = -f_7 * di_60[k]
                  - f_8 * di_67[k]
                  + f_9 * di_69[k]
                  + f_10 * di_78[k]
                  - f_11 * di_80[k];

        g_43[k] = f_12 * di_57[k]
                  + f_13 * di_62[k]
                  - f_14 * di_64[k]
                  + f_12 * di_71[k]
                  - f_14 * di_73[k]
                  + f_14 * di_75[k];
    }

#pragma omp simd aligned(di_60, di_67, di_69, di_78, di_80, di_82 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_44[k] = f_15 * di_60[k]
                  + f_16 * di_67[k]
                  - f_17 * di_69[k]
                  + f_15 * di_78[k]
                  - f_17 * di_80[k]
                  + f_18 * di_82[k];
    }

#pragma omp simd aligned(di_56, di_59, di_61, di_66, di_68, di_70, di_77, di_79, di_81, \
                         di_83 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_45[k] = -f_19 * di_56[k]
                  - f_20 * di_59[k]
                  + f_21 * di_61[k]
                  - f_20 * di_66[k]
                  + f_22 * di_68[k]
                  - f_23 * di_70[k]
                  - f_19 * di_77[k]
                  + f_21 * di_79[k]
                  - f_23 * di_81[k]
                  + f_24 * di_83[k];
    }

#pragma omp simd aligned(di_56, di_58, di_59, di_61, di_63, di_65, di_66, di_70, di_72, di_74, \
                         di_76, di_77, di_79, di_81 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_46[k] = f_15 * di_58[k]
                  + f_16 * di_63[k]
                  - f_17 * di_65[k]
                  + f_15 * di_72[k]
                  - f_17 * di_74[k]
                  + f_18 * di_76[k];

        g_47[k] = f_25 * di_56[k]
                  + f_25 * di_59[k]
                  - f_11 * di_61[k]
                  - f_25 * di_66[k]
                  + f_11 * di_70[k]
                  - f_25 * di_77[k]
                  + f_11 * di_79[k]
                  - f_11 * di_81[k];
    }

#pragma omp simd aligned(di_56, di_58, di_59, di_61, di_63, di_65, di_66, di_68, di_72, di_74, \
                         di_77, di_79 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_48[k] = -f_10 * di_58[k]
                  + f_8 * di_63[k]
                  + f_11 * di_65[k]
                  + f_7 * di_72[k]
                  - f_9 * di_74[k];

        g_49[k] = -f_26 * di_56[k]
                  + f_27 * di_59[k]
                  + f_28 * di_61[k]
                  + f_27 * di_66[k]
                  - f_29 * di_68[k]
                  - f_26 * di_77[k]
                  + f_28 * di_79[k];

        g_50[k] = f_4 * di_58[k]
                  - f_3 * di_63[k]
                  + f_2 * di_72[k];

        g_51[k] = f_30 * di_56[k]
                  - f_31 * di_59[k]
                  + f_31 * di_66[k]
                  - f_30 * di_77[k];
    }

#pragma omp simd aligned(di_1, di_4, di_6, di_11, di_15, di_22, di_85, di_88, di_90, di_95, \
                         di_99, di_106 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_52[k] = f_70 * di_1[k]
                  - f_36 * di_6[k]
                  + f_70 * di_15[k]
                  - f_70 * di_85[k]
                  + f_36 * di_90[k]
                  - f_70 * di_99[k];

        g_53[k] = f_69 * di_4[k]
                  - f_2 * di_11[k]
                  + f_32 * di_22[k]
                  - f_69 * di_88[k]
                  + f_2 * di_95[k]
                  - f_32 * di_106[k];
    }

#pragma omp simd aligned(di_1, di_8, di_15, di_17, di_85, di_92, di_99, \
                         di_101 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_54[k] = -f_71 * di_1[k]
                  + f_72 * di_8[k]
                  + f_71 * di_15[k]
                  - f_72 * di_17[k]
                  + f_71 * di_85[k]
                  - f_72 * di_92[k]
                  - f_71 * di_99[k]
                  + f_72 * di_101[k];
    }

#pragma omp simd aligned(di_4, di_11, di_13, di_22, di_24, di_88, di_95, di_97, di_106, \
                         di_108 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_55[k] = -f_73 * di_4[k]
                  - f_10 * di_11[k]
                  + f_74 * di_13[k]
                  + f_75 * di_22[k]
                  - f_76 * di_24[k]
                  + f_73 * di_88[k]
                  + f_10 * di_95[k]
                  - f_74 * di_97[k]
                  - f_75 * di_106[k]
                  + f_76 * di_108[k];
    }

#pragma omp simd aligned(di_1, di_6, di_8, di_15, di_17, di_19, di_85, di_90, di_92, di_99, \
                         di_101, di_103 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_56[k] = f_25 * di_1[k]
                  + f_12 * di_6[k]
                  - f_11 * di_8[k]
                  + f_25 * di_15[k]
                  - f_11 * di_17[k]
                  + f_11 * di_19[k]
                  - f_25 * di_85[k]
                  - f_12 * di_90[k]
                  + f_11 * di_92[k]
                  - f_25 * di_99[k]
                  + f_11 * di_101[k]
                  - f_11 * di_103[k];
    }

#pragma omp simd aligned(di_4, di_11, di_13, di_22, di_24, di_26, di_88, di_95, di_97, di_106, \
                         di_108, di_110 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_57[k] = f_62 * di_4[k]
                  + f_15 * di_11[k]
                  - f_16 * di_13[k]
                  + f_62 * di_22[k]
                  - f_16 * di_24[k]
                  + f_77 * di_26[k]
                  - f_62 * di_88[k]
                  - f_15 * di_95[k]
                  + f_16 * di_97[k]
                  - f_62 * di_106[k]
                  + f_16 * di_108[k]
                  - f_77 * di_110[k];
    }

#pragma omp simd aligned(di_0, di_3, di_5, di_10, di_12, di_14, di_21, di_23, di_25, di_27, \
                         di_84, di_87, di_89, di_94, di_96, di_98, di_105, di_107, di_109, \
                         di_111 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_58[k] = -f_78 * di_0[k]
                  - f_79 * di_3[k]
                  + f_80 * di_5[k]
                  - f_79 * di_10[k]
                  + f_21 * di_12[k]
                  - f_81 * di_14[k]
                  - f_78 * di_21[k]
                  + f_80 * di_23[k]
                  - f_81 * di_25[k]
                  + f_82 * di_27[k]
                  + f_78 * di_84[k]
                  + f_79 * di_87[k]
                  - f_80 * di_89[k]
                  + f_79 * di_94[k]
                  - f_21 * di_96[k]
                  + f_81 * di_98[k]
                  + f_78 * di_105[k]
                  - f_80 * di_107[k]
                  + f_81 * di_109[k]
                  - f_82 * di_111[k];
    }

#pragma omp simd aligned(di_2, di_7, di_9, di_16, di_18, di_20, di_86, di_91, di_93, di_100, \
                         di_102, di_104 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_59[k] = f_62 * di_2[k]
                  + f_15 * di_7[k]
                  - f_16 * di_9[k]
                  + f_62 * di_16[k]
                  - f_16 * di_18[k]
                  + f_77 * di_20[k]
                  - f_62 * di_86[k]
                  - f_15 * di_91[k]
                  + f_16 * di_93[k]
                  - f_62 * di_100[k]
                  + f_16 * di_102[k]
                  - f_77 * di_104[k];
    }

#pragma omp simd aligned(di_0, di_3, di_5, di_10, di_14, di_21, di_23, di_25, di_84, di_87, \
                         di_89, di_94, di_98, di_105, di_107, di_109 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_60[k] = f_83 * di_0[k]
                  + f_83 * di_3[k]
                  - f_76 * di_5[k]
                  - f_83 * di_10[k]
                  + f_76 * di_14[k]
                  - f_83 * di_21[k]
                  + f_76 * di_23[k]
                  - f_76 * di_25[k]
                  - f_83 * di_84[k]
                  - f_83 * di_87[k]
                  + f_76 * di_89[k]
                  + f_83 * di_94[k]
                  - f_76 * di_98[k]
                  + f_83 * di_105[k]
                  - f_76 * di_107[k]
                  + f_76 * di_109[k];
    }

#pragma omp simd aligned(di_2, di_7, di_9, di_16, di_18, di_86, di_91, di_93, di_100, \
                         di_102 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_61[k] = -f_75 * di_2[k]
                  + f_10 * di_7[k]
                  + f_76 * di_9[k]
                  + f_73 * di_16[k]
                  - f_74 * di_18[k]
                  + f_75 * di_86[k]
                  - f_10 * di_91[k]
                  - f_76 * di_93[k]
                  - f_73 * di_100[k]
                  + f_74 * di_102[k];
    }

#pragma omp simd aligned(di_0, di_3, di_5, di_10, di_12, di_21, di_23, di_84, di_87, di_89, \
                         di_94, di_96, di_105, di_107 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_62[k] = -f_84 * di_0[k]
                  + f_85 * di_3[k]
                  + f_27 * di_5[k]
                  + f_85 * di_10[k]
                  - f_86 * di_12[k]
                  - f_84 * di_21[k]
                  + f_27 * di_23[k]
                  + f_84 * di_84[k]
                  - f_85 * di_87[k]
                  - f_27 * di_89[k]
                  - f_85 * di_94[k]
                  + f_86 * di_96[k]
                  + f_84 * di_105[k]
                  - f_27 * di_107[k];
    }

#pragma omp simd aligned(di_0, di_2, di_3, di_7, di_10, di_16, di_21, di_84, di_86, di_87, \
                         di_91, di_94, di_100, di_105 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_63[k] = f_32 * di_2[k]
                  - f_2 * di_7[k]
                  + f_69 * di_16[k]
                  - f_32 * di_86[k]
                  + f_2 * di_91[k]
                  - f_69 * di_100[k];

        g_64[k] = f_87 * di_0[k]
                  - f_88 * di_3[k]
                  + f_88 * di_10[k]
                  - f_87 * di_21[k]
                  - f_87 * di_84[k]
                  + f_88 * di_87[k]
                  - f_88 * di_94[k]
                  + f_87 * di_105[k];
    }
}

}  // namespace simdtrf
