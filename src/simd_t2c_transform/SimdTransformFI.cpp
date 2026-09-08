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


#include "SimdTransformFI.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_fi(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t fi,
             const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 0.28125 * std::sqrt(1155.0);
    const auto f_1 = 0.9375 * std::sqrt(1155.0);
    const auto f_2 = 0.09375 * std::sqrt(1155.0);
    const auto f_3 = 0.3125 * std::sqrt(1155.0);
    const auto f_4 = 1.40625 * std::sqrt(385.0);
    const auto f_5 = 2.8125 * std::sqrt(385.0);
    const auto f_6 = 0.28125 * std::sqrt(385.0);
    const auto f_7 = 0.46875 * std::sqrt(385.0);
    const auto f_8 = 0.9375 * std::sqrt(385.0);
    const auto f_9 = 0.09375 * std::sqrt(385.0);
    const auto f_10 = 0.5625 * std::sqrt(70.0);
    const auto f_11 = 5.625 * std::sqrt(70.0);
    const auto f_12 = 0.1875 * std::sqrt(70.0);
    const auto f_13 = 1.875 * std::sqrt(70.0);
    const auto f_14 = 4.21875 * std::sqrt(21.0);
    const auto f_15 = 2.8125 * std::sqrt(21.0);
    const auto f_16 = 11.25 * std::sqrt(21.0);
    const auto f_17 = 1.40625 * std::sqrt(21.0);
    const auto f_18 = 3.75 * std::sqrt(21.0);
    const auto f_19 = 0.9375 * std::sqrt(21.0);
    const auto f_20 = 0.46875 * std::sqrt(21.0);
    const auto f_21 = 1.25 * std::sqrt(21.0);
    const auto f_22 = 7.5 * std::sqrt(21.0);
    const auto f_23 = 0.15625 * std::sqrt(21.0);
    const auto f_24 = 0.3125 * std::sqrt(21.0);
    const auto f_25 = 2.5 * std::sqrt(21.0);
    const auto f_26 = 0.46875 * std::sqrt(210.0);
    const auto f_27 = 0.9375 * std::sqrt(210.0);
    const auto f_28 = 1.875 * std::sqrt(210.0);
    const auto f_29 = 0.75 * std::sqrt(210.0);
    const auto f_30 = 0.15625 * std::sqrt(210.0);
    const auto f_31 = 0.3125 * std::sqrt(210.0);
    const auto f_32 = 0.625 * std::sqrt(210.0);
    const auto f_33 = 0.25 * std::sqrt(210.0);
    const auto f_34 = 0.234375 * std::sqrt(10.0);
    const auto f_35 = 0.703125 * std::sqrt(10.0);
    const auto f_36 = 4.21875 * std::sqrt(10.0);
    const auto f_37 = 8.4375 * std::sqrt(10.0);
    const auto f_38 = 5.625 * std::sqrt(10.0);
    const auto f_39 = 0.75 * std::sqrt(10.0);
    const auto f_40 = 0.078125 * std::sqrt(10.0);
    const auto f_41 = 1.40625 * std::sqrt(10.0);
    const auto f_42 = 2.8125 * std::sqrt(10.0);
    const auto f_43 = 1.875 * std::sqrt(10.0);
    const auto f_44 = 0.25 * std::sqrt(10.0);
    const auto f_45 = 0.234375 * std::sqrt(21.0);
    const auto f_46 = 0.078125 * std::sqrt(21.0);
    const auto f_47 = 0.140625 * std::sqrt(70.0);
    const auto f_48 = 0.703125 * std::sqrt(70.0);
    const auto f_49 = 1.40625 * std::sqrt(70.0);
    const auto f_50 = 8.4375 * std::sqrt(70.0);
    const auto f_51 = 0.046875 * std::sqrt(70.0);
    const auto f_52 = 0.234375 * std::sqrt(70.0);
    const auto f_53 = 0.46875 * std::sqrt(70.0);
    const auto f_54 = 2.8125 * std::sqrt(70.0);
    const auto f_55 = 0.046875 * std::sqrt(1155.0);
    const auto f_56 = 0.703125 * std::sqrt(1155.0);
    const auto f_57 = 0.015625 * std::sqrt(1155.0);
    const auto f_58 = 0.234375 * std::sqrt(1155.0);
    const auto f_59 = 0.5625 * std::sqrt(770.0);
    const auto f_60 = 1.875 * std::sqrt(770.0);
    const auto f_61 = 0.9375 * std::sqrt(2310.0);
    const auto f_62 = 1.875 * std::sqrt(2310.0);
    const auto f_63 = 0.1875 * std::sqrt(2310.0);
    const auto f_64 = 0.75 * std::sqrt(105.0);
    const auto f_65 = 7.5 * std::sqrt(105.0);
    const auto f_66 = 8.4375 * std::sqrt(14.0);
    const auto f_67 = 5.625 * std::sqrt(14.0);
    const auto f_68 = 22.5 * std::sqrt(14.0);
    const auto f_69 = 2.8125 * std::sqrt(14.0);
    const auto f_70 = 7.5 * std::sqrt(14.0);
    const auto f_71 = 0.9375 * std::sqrt(14.0);
    const auto f_72 = 1.875 * std::sqrt(14.0);
    const auto f_73 = 15.0 * std::sqrt(14.0);
    const auto f_74 = 1.875 * std::sqrt(35.0);
    const auto f_75 = 3.75 * std::sqrt(35.0);
    const auto f_76 = 7.5 * std::sqrt(35.0);
    const auto f_77 = 3.0 * std::sqrt(35.0);
    const auto f_78 = 0.3125 * std::sqrt(15.0);
    const auto f_79 = 0.9375 * std::sqrt(15.0);
    const auto f_80 = 5.625 * std::sqrt(15.0);
    const auto f_81 = 11.25 * std::sqrt(15.0);
    const auto f_82 = 7.5 * std::sqrt(15.0);
    const auto f_83 = std::sqrt(15.0);
    const auto f_84 = 0.46875 * std::sqrt(14.0);
    const auto f_85 = 0.1875 * std::sqrt(105.0);
    const auto f_86 = 0.9375 * std::sqrt(105.0);
    const auto f_87 = 1.875 * std::sqrt(105.0);
    const auto f_88 = 11.25 * std::sqrt(105.0);
    const auto f_89 = 0.09375 * std::sqrt(770.0);
    const auto f_90 = 1.40625 * std::sqrt(770.0);
    const auto f_91 = 0.28125 * std::sqrt(77.0);
    const auto f_92 = 0.9375 * std::sqrt(77.0);
    const auto f_93 = 1.125 * std::sqrt(77.0);
    const auto f_94 = 3.75 * std::sqrt(77.0);
    const auto f_95 = 0.46875 * std::sqrt(231.0);
    const auto f_96 = 0.9375 * std::sqrt(231.0);
    const auto f_97 = 0.09375 * std::sqrt(231.0);
    const auto f_98 = 1.875 * std::sqrt(231.0);
    const auto f_99 = 3.75 * std::sqrt(231.0);
    const auto f_100 = 0.375 * std::sqrt(231.0);
    const auto f_101 = 0.1875 * std::sqrt(42.0);
    const auto f_102 = 1.875 * std::sqrt(42.0);
    const auto f_103 = 0.75 * std::sqrt(42.0);
    const auto f_104 = 7.5 * std::sqrt(42.0);
    const auto f_105 = 0.84375 * std::sqrt(35.0);
    const auto f_106 = 0.5625 * std::sqrt(35.0);
    const auto f_107 = 2.25 * std::sqrt(35.0);
    const auto f_108 = 0.28125 * std::sqrt(35.0);
    const auto f_109 = 0.75 * std::sqrt(35.0);
    const auto f_110 = 3.375 * std::sqrt(35.0);
    const auto f_111 = 9.0 * std::sqrt(35.0);
    const auto f_112 = 1.125 * std::sqrt(35.0);
    const auto f_113 = 0.09375 * std::sqrt(35.0);
    const auto f_114 = 0.1875 * std::sqrt(35.0);
    const auto f_115 = 1.5 * std::sqrt(35.0);
    const auto f_116 = 0.375 * std::sqrt(35.0);
    const auto f_117 = 6.0 * std::sqrt(35.0);
    const auto f_118 = 0.75 * std::sqrt(14.0);
    const auto f_119 = 3.75 * std::sqrt(14.0);
    const auto f_120 = 3.0 * std::sqrt(14.0);
    const auto f_121 = 0.078125 * std::sqrt(6.0);
    const auto f_122 = 0.234375 * std::sqrt(6.0);
    const auto f_123 = 1.40625 * std::sqrt(6.0);
    const auto f_124 = 2.8125 * std::sqrt(6.0);
    const auto f_125 = 1.875 * std::sqrt(6.0);
    const auto f_126 = 0.25 * std::sqrt(6.0);
    const auto f_127 = 0.3125 * std::sqrt(6.0);
    const auto f_128 = 0.9375 * std::sqrt(6.0);
    const auto f_129 = 5.625 * std::sqrt(6.0);
    const auto f_130 = 11.25 * std::sqrt(6.0);
    const auto f_131 = 7.5 * std::sqrt(6.0);
    const auto f_132 = std::sqrt(6.0);
    const auto f_133 = 0.046875 * std::sqrt(35.0);
    const auto f_134 = 0.046875 * std::sqrt(42.0);
    const auto f_135 = 0.234375 * std::sqrt(42.0);
    const auto f_136 = 0.46875 * std::sqrt(42.0);
    const auto f_137 = 2.8125 * std::sqrt(42.0);
    const auto f_138 = 0.9375 * std::sqrt(42.0);
    const auto f_139 = 11.25 * std::sqrt(42.0);
    const auto f_140 = 0.046875 * std::sqrt(77.0);
    const auto f_141 = 0.703125 * std::sqrt(77.0);
    const auto f_142 = 0.1875 * std::sqrt(77.0);
    const auto f_143 = 2.8125 * std::sqrt(77.0);
    const auto f_144 = 0.28125 * std::sqrt(462.0);
    const auto f_145 = 0.9375 * std::sqrt(462.0);
    const auto f_146 = 0.1875 * std::sqrt(462.0);
    const auto f_147 = 0.625 * std::sqrt(462.0);
    const auto f_148 = 1.40625 * std::sqrt(154.0);
    const auto f_149 = 2.8125 * std::sqrt(154.0);
    const auto f_150 = 0.28125 * std::sqrt(154.0);
    const auto f_151 = 0.9375 * std::sqrt(154.0);
    const auto f_152 = 1.875 * std::sqrt(154.0);
    const auto f_153 = 0.1875 * std::sqrt(154.0);
    const auto f_154 = 1.125 * std::sqrt(7.0);
    const auto f_155 = 11.25 * std::sqrt(7.0);
    const auto f_156 = 0.75 * std::sqrt(7.0);
    const auto f_157 = 7.5 * std::sqrt(7.0);
    const auto f_158 = 0.84375 * std::sqrt(210.0);
    const auto f_159 = 0.5625 * std::sqrt(210.0);
    const auto f_160 = 2.25 * std::sqrt(210.0);
    const auto f_161 = 0.28125 * std::sqrt(210.0);
    const auto f_162 = 0.375 * std::sqrt(210.0);
    const auto f_163 = 1.5 * std::sqrt(210.0);
    const auto f_164 = 0.1875 * std::sqrt(210.0);
    const auto f_165 = 0.5 * std::sqrt(210.0);
    const auto f_166 = 0.09375 * std::sqrt(210.0);
    const auto f_167 = 0.0625 * std::sqrt(210.0);
    const auto f_168 = 0.125 * std::sqrt(210.0);
    const auto f_169 = std::sqrt(210.0);
    const auto f_170 = 1.875 * std::sqrt(21.0);
    const auto f_171 = 1.5 * std::sqrt(21.0);
    const auto f_172 = 0.625 * std::sqrt(21.0);
    const auto f_173 = std::sqrt(21.0);
    const auto f_174 = 0.046875 * std::sqrt(210.0);
    const auto f_175 = 0.03125 * std::sqrt(210.0);
    const auto f_176 = 0.28125 * std::sqrt(7.0);
    const auto f_177 = 1.40625 * std::sqrt(7.0);
    const auto f_178 = 2.8125 * std::sqrt(7.0);
    const auto f_179 = 16.875 * std::sqrt(7.0);
    const auto f_180 = 0.1875 * std::sqrt(7.0);
    const auto f_181 = 0.9375 * std::sqrt(7.0);
    const auto f_182 = 1.875 * std::sqrt(7.0);
    const auto f_183 = 0.046875 * std::sqrt(462.0);
    const auto f_184 = 0.703125 * std::sqrt(462.0);
    const auto f_185 = 0.03125 * std::sqrt(462.0);
    const auto f_186 = 0.46875 * std::sqrt(462.0);
    const auto f_187 = 0.28125 * std::sqrt(770.0);
    const auto f_188 = 0.9375 * std::sqrt(770.0);
    const auto f_189 = 0.46875 * std::sqrt(2310.0);
    const auto f_190 = 0.09375 * std::sqrt(2310.0);
    const auto f_191 = 0.375 * std::sqrt(105.0);
    const auto f_192 = 3.75 * std::sqrt(105.0);
    const auto f_193 = 4.21875 * std::sqrt(14.0);
    const auto f_194 = 11.25 * std::sqrt(14.0);
    const auto f_195 = 1.40625 * std::sqrt(14.0);
    const auto f_196 = 0.9375 * std::sqrt(35.0);
    const auto f_197 = 0.15625 * std::sqrt(15.0);
    const auto f_198 = 0.46875 * std::sqrt(15.0);
    const auto f_199 = 2.8125 * std::sqrt(15.0);
    const auto f_200 = 3.75 * std::sqrt(15.0);
    const auto f_201 = 0.5 * std::sqrt(15.0);
    const auto f_202 = 0.234375 * std::sqrt(14.0);
    const auto f_203 = 0.09375 * std::sqrt(105.0);
    const auto f_204 = 0.46875 * std::sqrt(105.0);
    const auto f_205 = 5.625 * std::sqrt(105.0);
    const auto f_206 = 0.046875 * std::sqrt(770.0);
    const auto f_207 = 0.703125 * std::sqrt(770.0);

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
    auto *g_65 = values + 65 * nvalues;
    auto *g_66 = values + 66 * nvalues;
    auto *g_67 = values + 67 * nvalues;
    auto *g_68 = values + 68 * nvalues;
    auto *g_69 = values + 69 * nvalues;
    auto *g_70 = values + 70 * nvalues;
    auto *g_71 = values + 71 * nvalues;
    auto *g_72 = values + 72 * nvalues;
    auto *g_73 = values + 73 * nvalues;
    auto *g_74 = values + 74 * nvalues;
    auto *g_75 = values + 75 * nvalues;
    auto *g_76 = values + 76 * nvalues;
    auto *g_77 = values + 77 * nvalues;
    auto *g_78 = values + 78 * nvalues;
    auto *g_79 = values + 79 * nvalues;
    auto *g_80 = values + 80 * nvalues;
    auto *g_81 = values + 81 * nvalues;
    auto *g_82 = values + 82 * nvalues;
    auto *g_83 = values + 83 * nvalues;
    auto *g_84 = values + 84 * nvalues;
    auto *g_85 = values + 85 * nvalues;
    auto *g_86 = values + 86 * nvalues;
    auto *g_87 = values + 87 * nvalues;
    auto *g_88 = values + 88 * nvalues;
    auto *g_89 = values + 89 * nvalues;
    auto *g_90 = values + 90 * nvalues;

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_1 = buffer.data(fi + 1);
    const auto *fi_2 = buffer.data(fi + 2);
    const auto *fi_3 = buffer.data(fi + 3);
    const auto *fi_4 = buffer.data(fi + 4);
    const auto *fi_5 = buffer.data(fi + 5);
    const auto *fi_6 = buffer.data(fi + 6);
    const auto *fi_7 = buffer.data(fi + 7);
    const auto *fi_8 = buffer.data(fi + 8);
    const auto *fi_9 = buffer.data(fi + 9);
    const auto *fi_10 = buffer.data(fi + 10);
    const auto *fi_11 = buffer.data(fi + 11);
    const auto *fi_12 = buffer.data(fi + 12);
    const auto *fi_13 = buffer.data(fi + 13);
    const auto *fi_14 = buffer.data(fi + 14);
    const auto *fi_15 = buffer.data(fi + 15);
    const auto *fi_16 = buffer.data(fi + 16);
    const auto *fi_17 = buffer.data(fi + 17);
    const auto *fi_18 = buffer.data(fi + 18);
    const auto *fi_19 = buffer.data(fi + 19);
    const auto *fi_20 = buffer.data(fi + 20);
    const auto *fi_21 = buffer.data(fi + 21);
    const auto *fi_22 = buffer.data(fi + 22);
    const auto *fi_23 = buffer.data(fi + 23);
    const auto *fi_24 = buffer.data(fi + 24);
    const auto *fi_25 = buffer.data(fi + 25);
    const auto *fi_26 = buffer.data(fi + 26);
    const auto *fi_27 = buffer.data(fi + 27);
    const auto *fi_28 = buffer.data(fi + 28);
    const auto *fi_29 = buffer.data(fi + 29);
    const auto *fi_30 = buffer.data(fi + 30);
    const auto *fi_31 = buffer.data(fi + 31);
    const auto *fi_32 = buffer.data(fi + 32);
    const auto *fi_33 = buffer.data(fi + 33);
    const auto *fi_34 = buffer.data(fi + 34);
    const auto *fi_35 = buffer.data(fi + 35);
    const auto *fi_36 = buffer.data(fi + 36);
    const auto *fi_37 = buffer.data(fi + 37);
    const auto *fi_38 = buffer.data(fi + 38);
    const auto *fi_39 = buffer.data(fi + 39);
    const auto *fi_40 = buffer.data(fi + 40);
    const auto *fi_41 = buffer.data(fi + 41);
    const auto *fi_42 = buffer.data(fi + 42);
    const auto *fi_43 = buffer.data(fi + 43);
    const auto *fi_44 = buffer.data(fi + 44);
    const auto *fi_45 = buffer.data(fi + 45);
    const auto *fi_46 = buffer.data(fi + 46);
    const auto *fi_47 = buffer.data(fi + 47);
    const auto *fi_48 = buffer.data(fi + 48);
    const auto *fi_49 = buffer.data(fi + 49);
    const auto *fi_50 = buffer.data(fi + 50);
    const auto *fi_51 = buffer.data(fi + 51);
    const auto *fi_52 = buffer.data(fi + 52);
    const auto *fi_53 = buffer.data(fi + 53);
    const auto *fi_54 = buffer.data(fi + 54);
    const auto *fi_55 = buffer.data(fi + 55);
    const auto *fi_56 = buffer.data(fi + 56);
    const auto *fi_57 = buffer.data(fi + 57);
    const auto *fi_58 = buffer.data(fi + 58);
    const auto *fi_59 = buffer.data(fi + 59);
    const auto *fi_60 = buffer.data(fi + 60);
    const auto *fi_61 = buffer.data(fi + 61);
    const auto *fi_62 = buffer.data(fi + 62);
    const auto *fi_63 = buffer.data(fi + 63);
    const auto *fi_64 = buffer.data(fi + 64);
    const auto *fi_65 = buffer.data(fi + 65);
    const auto *fi_66 = buffer.data(fi + 66);
    const auto *fi_67 = buffer.data(fi + 67);
    const auto *fi_68 = buffer.data(fi + 68);
    const auto *fi_69 = buffer.data(fi + 69);
    const auto *fi_70 = buffer.data(fi + 70);
    const auto *fi_71 = buffer.data(fi + 71);
    const auto *fi_72 = buffer.data(fi + 72);
    const auto *fi_73 = buffer.data(fi + 73);
    const auto *fi_74 = buffer.data(fi + 74);
    const auto *fi_75 = buffer.data(fi + 75);
    const auto *fi_76 = buffer.data(fi + 76);
    const auto *fi_77 = buffer.data(fi + 77);
    const auto *fi_78 = buffer.data(fi + 78);
    const auto *fi_79 = buffer.data(fi + 79);
    const auto *fi_80 = buffer.data(fi + 80);
    const auto *fi_81 = buffer.data(fi + 81);
    const auto *fi_82 = buffer.data(fi + 82);
    const auto *fi_83 = buffer.data(fi + 83);
    const auto *fi_84 = buffer.data(fi + 84);
    const auto *fi_85 = buffer.data(fi + 85);
    const auto *fi_86 = buffer.data(fi + 86);
    const auto *fi_87 = buffer.data(fi + 87);
    const auto *fi_88 = buffer.data(fi + 88);
    const auto *fi_89 = buffer.data(fi + 89);
    const auto *fi_90 = buffer.data(fi + 90);
    const auto *fi_91 = buffer.data(fi + 91);
    const auto *fi_92 = buffer.data(fi + 92);
    const auto *fi_93 = buffer.data(fi + 93);
    const auto *fi_94 = buffer.data(fi + 94);
    const auto *fi_95 = buffer.data(fi + 95);
    const auto *fi_96 = buffer.data(fi + 96);
    const auto *fi_97 = buffer.data(fi + 97);
    const auto *fi_98 = buffer.data(fi + 98);
    const auto *fi_99 = buffer.data(fi + 99);
    const auto *fi_100 = buffer.data(fi + 100);
    const auto *fi_101 = buffer.data(fi + 101);
    const auto *fi_102 = buffer.data(fi + 102);
    const auto *fi_103 = buffer.data(fi + 103);
    const auto *fi_104 = buffer.data(fi + 104);
    const auto *fi_105 = buffer.data(fi + 105);
    const auto *fi_106 = buffer.data(fi + 106);
    const auto *fi_107 = buffer.data(fi + 107);
    const auto *fi_108 = buffer.data(fi + 108);
    const auto *fi_109 = buffer.data(fi + 109);
    const auto *fi_110 = buffer.data(fi + 110);
    const auto *fi_111 = buffer.data(fi + 111);
    const auto *fi_112 = buffer.data(fi + 112);
    const auto *fi_113 = buffer.data(fi + 113);
    const auto *fi_114 = buffer.data(fi + 114);
    const auto *fi_115 = buffer.data(fi + 115);
    const auto *fi_116 = buffer.data(fi + 116);
    const auto *fi_117 = buffer.data(fi + 117);
    const auto *fi_118 = buffer.data(fi + 118);
    const auto *fi_119 = buffer.data(fi + 119);
    const auto *fi_120 = buffer.data(fi + 120);
    const auto *fi_121 = buffer.data(fi + 121);
    const auto *fi_122 = buffer.data(fi + 122);
    const auto *fi_123 = buffer.data(fi + 123);
    const auto *fi_124 = buffer.data(fi + 124);
    const auto *fi_125 = buffer.data(fi + 125);
    const auto *fi_126 = buffer.data(fi + 126);
    const auto *fi_127 = buffer.data(fi + 127);
    const auto *fi_128 = buffer.data(fi + 128);
    const auto *fi_129 = buffer.data(fi + 129);
    const auto *fi_130 = buffer.data(fi + 130);
    const auto *fi_131 = buffer.data(fi + 131);
    const auto *fi_132 = buffer.data(fi + 132);
    const auto *fi_133 = buffer.data(fi + 133);
    const auto *fi_134 = buffer.data(fi + 134);
    const auto *fi_135 = buffer.data(fi + 135);
    const auto *fi_136 = buffer.data(fi + 136);
    const auto *fi_137 = buffer.data(fi + 137);
    const auto *fi_138 = buffer.data(fi + 138);
    const auto *fi_139 = buffer.data(fi + 139);
    const auto *fi_140 = buffer.data(fi + 140);
    const auto *fi_141 = buffer.data(fi + 141);
    const auto *fi_142 = buffer.data(fi + 142);
    const auto *fi_143 = buffer.data(fi + 143);
    const auto *fi_144 = buffer.data(fi + 144);
    const auto *fi_145 = buffer.data(fi + 145);
    const auto *fi_146 = buffer.data(fi + 146);
    const auto *fi_147 = buffer.data(fi + 147);
    const auto *fi_148 = buffer.data(fi + 148);
    const auto *fi_149 = buffer.data(fi + 149);
    const auto *fi_150 = buffer.data(fi + 150);
    const auto *fi_151 = buffer.data(fi + 151);
    const auto *fi_152 = buffer.data(fi + 152);
    const auto *fi_153 = buffer.data(fi + 153);
    const auto *fi_154 = buffer.data(fi + 154);
    const auto *fi_155 = buffer.data(fi + 155);
    const auto *fi_156 = buffer.data(fi + 156);
    const auto *fi_157 = buffer.data(fi + 157);
    const auto *fi_158 = buffer.data(fi + 158);
    const auto *fi_159 = buffer.data(fi + 159);
    const auto *fi_160 = buffer.data(fi + 160);
    const auto *fi_161 = buffer.data(fi + 161);
    const auto *fi_162 = buffer.data(fi + 162);
    const auto *fi_163 = buffer.data(fi + 163);
    const auto *fi_164 = buffer.data(fi + 164);
    const auto *fi_165 = buffer.data(fi + 165);
    const auto *fi_166 = buffer.data(fi + 166);
    const auto *fi_167 = buffer.data(fi + 167);
    const auto *fi_168 = buffer.data(fi + 168);
    const auto *fi_169 = buffer.data(fi + 169);
    const auto *fi_170 = buffer.data(fi + 170);
    const auto *fi_171 = buffer.data(fi + 171);
    const auto *fi_172 = buffer.data(fi + 172);
    const auto *fi_173 = buffer.data(fi + 173);
    const auto *fi_174 = buffer.data(fi + 174);
    const auto *fi_175 = buffer.data(fi + 175);
    const auto *fi_176 = buffer.data(fi + 176);
    const auto *fi_177 = buffer.data(fi + 177);
    const auto *fi_178 = buffer.data(fi + 178);
    const auto *fi_179 = buffer.data(fi + 179);
    const auto *fi_180 = buffer.data(fi + 180);
    const auto *fi_181 = buffer.data(fi + 181);
    const auto *fi_182 = buffer.data(fi + 182);
    const auto *fi_183 = buffer.data(fi + 183);
    const auto *fi_184 = buffer.data(fi + 184);
    const auto *fi_185 = buffer.data(fi + 185);
    const auto *fi_186 = buffer.data(fi + 186);
    const auto *fi_187 = buffer.data(fi + 187);
    const auto *fi_188 = buffer.data(fi + 188);
    const auto *fi_189 = buffer.data(fi + 189);
    const auto *fi_190 = buffer.data(fi + 190);
    const auto *fi_191 = buffer.data(fi + 191);
    const auto *fi_192 = buffer.data(fi + 192);
    const auto *fi_193 = buffer.data(fi + 193);
    const auto *fi_194 = buffer.data(fi + 194);
    const auto *fi_195 = buffer.data(fi + 195);
    const auto *fi_196 = buffer.data(fi + 196);
    const auto *fi_197 = buffer.data(fi + 197);
    const auto *fi_198 = buffer.data(fi + 198);
    const auto *fi_199 = buffer.data(fi + 199);
    const auto *fi_200 = buffer.data(fi + 200);
    const auto *fi_201 = buffer.data(fi + 201);
    const auto *fi_202 = buffer.data(fi + 202);
    const auto *fi_203 = buffer.data(fi + 203);
    const auto *fi_204 = buffer.data(fi + 204);
    const auto *fi_205 = buffer.data(fi + 205);
    const auto *fi_206 = buffer.data(fi + 206);
    const auto *fi_207 = buffer.data(fi + 207);
    const auto *fi_208 = buffer.data(fi + 208);
    const auto *fi_209 = buffer.data(fi + 209);
    const auto *fi_210 = buffer.data(fi + 210);
    const auto *fi_211 = buffer.data(fi + 211);
    const auto *fi_212 = buffer.data(fi + 212);
    const auto *fi_213 = buffer.data(fi + 213);
    const auto *fi_214 = buffer.data(fi + 214);
    const auto *fi_215 = buffer.data(fi + 215);
    const auto *fi_216 = buffer.data(fi + 216);
    const auto *fi_217 = buffer.data(fi + 217);
    const auto *fi_218 = buffer.data(fi + 218);
    const auto *fi_219 = buffer.data(fi + 219);
    const auto *fi_220 = buffer.data(fi + 220);
    const auto *fi_221 = buffer.data(fi + 221);
    const auto *fi_222 = buffer.data(fi + 222);
    const auto *fi_223 = buffer.data(fi + 223);
    const auto *fi_224 = buffer.data(fi + 224);
    const auto *fi_225 = buffer.data(fi + 225);
    const auto *fi_226 = buffer.data(fi + 226);
    const auto *fi_227 = buffer.data(fi + 227);
    const auto *fi_228 = buffer.data(fi + 228);
    const auto *fi_229 = buffer.data(fi + 229);
    const auto *fi_230 = buffer.data(fi + 230);
    const auto *fi_231 = buffer.data(fi + 231);
    const auto *fi_232 = buffer.data(fi + 232);
    const auto *fi_233 = buffer.data(fi + 233);
    const auto *fi_234 = buffer.data(fi + 234);
    const auto *fi_235 = buffer.data(fi + 235);
    const auto *fi_236 = buffer.data(fi + 236);
    const auto *fi_237 = buffer.data(fi + 237);
    const auto *fi_238 = buffer.data(fi + 238);
    const auto *fi_239 = buffer.data(fi + 239);
    const auto *fi_240 = buffer.data(fi + 240);
    const auto *fi_241 = buffer.data(fi + 241);
    const auto *fi_242 = buffer.data(fi + 242);
    const auto *fi_243 = buffer.data(fi + 243);
    const auto *fi_244 = buffer.data(fi + 244);
    const auto *fi_245 = buffer.data(fi + 245);
    const auto *fi_246 = buffer.data(fi + 246);
    const auto *fi_247 = buffer.data(fi + 247);
    const auto *fi_248 = buffer.data(fi + 248);
    const auto *fi_249 = buffer.data(fi + 249);
    const auto *fi_250 = buffer.data(fi + 250);
    const auto *fi_251 = buffer.data(fi + 251);
    const auto *fi_252 = buffer.data(fi + 252);
    const auto *fi_253 = buffer.data(fi + 253);
    const auto *fi_254 = buffer.data(fi + 254);
    const auto *fi_255 = buffer.data(fi + 255);
    const auto *fi_256 = buffer.data(fi + 256);
    const auto *fi_257 = buffer.data(fi + 257);
    const auto *fi_258 = buffer.data(fi + 258);
    const auto *fi_259 = buffer.data(fi + 259);
    const auto *fi_260 = buffer.data(fi + 260);
    const auto *fi_261 = buffer.data(fi + 261);
    const auto *fi_262 = buffer.data(fi + 262);
    const auto *fi_263 = buffer.data(fi + 263);
    const auto *fi_264 = buffer.data(fi + 264);
    const auto *fi_265 = buffer.data(fi + 265);
    const auto *fi_266 = buffer.data(fi + 266);
    const auto *fi_267 = buffer.data(fi + 267);
    const auto *fi_268 = buffer.data(fi + 268);
    const auto *fi_269 = buffer.data(fi + 269);
    const auto *fi_270 = buffer.data(fi + 270);
    const auto *fi_271 = buffer.data(fi + 271);
    const auto *fi_272 = buffer.data(fi + 272);
    const auto *fi_273 = buffer.data(fi + 273);
    const auto *fi_274 = buffer.data(fi + 274);
    const auto *fi_275 = buffer.data(fi + 275);
    const auto *fi_276 = buffer.data(fi + 276);
    const auto *fi_277 = buffer.data(fi + 277);
    const auto *fi_278 = buffer.data(fi + 278);
    const auto *fi_279 = buffer.data(fi + 279);

#pragma omp simd aligned(fi_29, fi_32, fi_34, fi_39, fi_43, fi_50, fi_169, fi_172, fi_174, \
                         fi_179, fi_183, fi_190 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * fi_29[k]
                 - f_1 * fi_34[k]
                 + f_0 * fi_43[k]
                 - f_2 * fi_169[k]
                 + f_3 * fi_174[k]
                 - f_2 * fi_183[k];

        g_1[k] = f_4 * fi_32[k]
                 - f_5 * fi_39[k]
                 + f_6 * fi_50[k]
                 - f_7 * fi_172[k]
                 + f_8 * fi_179[k]
                 - f_9 * fi_190[k];
    }

#pragma omp simd aligned(fi_29, fi_36, fi_43, fi_45, fi_169, fi_176, fi_183, \
                         fi_185 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_2[k] = -f_10 * fi_29[k]
                 + f_11 * fi_36[k]
                 + f_10 * fi_43[k]
                 - f_11 * fi_45[k]
                 + f_12 * fi_169[k]
                 - f_13 * fi_176[k]
                 - f_12 * fi_183[k]
                 + f_13 * fi_185[k];
    }

#pragma omp simd aligned(fi_32, fi_39, fi_41, fi_50, fi_52, fi_172, fi_179, fi_181, fi_190, \
                         fi_192 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = -f_14 * fi_32[k]
                 - f_15 * fi_39[k]
                 + f_16 * fi_41[k]
                 + f_17 * fi_50[k]
                 - f_18 * fi_52[k]
                 + f_17 * fi_172[k]
                 + f_19 * fi_179[k]
                 - f_18 * fi_181[k]
                 - f_20 * fi_190[k]
                 + f_21 * fi_192[k];
    }

#pragma omp simd aligned(fi_29, fi_34, fi_36, fi_43, fi_45, fi_47, fi_169, fi_174, fi_176, \
                         fi_183, fi_185, fi_187 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = f_20 * fi_29[k]
                 + f_19 * fi_34[k]
                 - f_22 * fi_36[k]
                 + f_20 * fi_43[k]
                 - f_22 * fi_45[k]
                 + f_22 * fi_47[k]
                 - f_23 * fi_169[k]
                 - f_24 * fi_174[k]
                 + f_25 * fi_176[k]
                 - f_23 * fi_183[k]
                 + f_25 * fi_185[k]
                 - f_25 * fi_187[k];
    }

#pragma omp simd aligned(fi_32, fi_39, fi_41, fi_50, fi_52, fi_54, fi_172, fi_179, fi_181, \
                         fi_190, fi_192, fi_194 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = f_26 * fi_32[k]
                 + f_27 * fi_39[k]
                 - f_28 * fi_41[k]
                 + f_26 * fi_50[k]
                 - f_28 * fi_52[k]
                 + f_29 * fi_54[k]
                 - f_30 * fi_172[k]
                 - f_31 * fi_179[k]
                 + f_32 * fi_181[k]
                 - f_30 * fi_190[k]
                 + f_32 * fi_192[k]
                 - f_33 * fi_194[k];
    }

#pragma omp simd aligned(fi_28, fi_31, fi_33, fi_38, fi_40, fi_42, fi_49, fi_51, fi_53, fi_55, \
                         fi_168, fi_171, fi_173, fi_178, fi_180, fi_182, fi_189, fi_191, \
                         fi_193, fi_195 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = -f_34 * fi_28[k]
                 - f_35 * fi_31[k]
                 + f_36 * fi_33[k]
                 - f_35 * fi_38[k]
                 + f_37 * fi_40[k]
                 - f_38 * fi_42[k]
                 - f_34 * fi_49[k]
                 + f_36 * fi_51[k]
                 - f_38 * fi_53[k]
                 + f_39 * fi_55[k]
                 + f_40 * fi_168[k]
                 + f_34 * fi_171[k]
                 - f_41 * fi_173[k]
                 + f_34 * fi_178[k]
                 - f_42 * fi_180[k]
                 + f_43 * fi_182[k]
                 + f_40 * fi_189[k]
                 - f_41 * fi_191[k]
                 + f_43 * fi_193[k]
                 - f_44 * fi_195[k];
    }

#pragma omp simd aligned(fi_30, fi_35, fi_37, fi_44, fi_46, fi_48, fi_170, fi_175, fi_177, \
                         fi_184, fi_186, fi_188 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = f_26 * fi_30[k]
                 + f_27 * fi_35[k]
                 - f_28 * fi_37[k]
                 + f_26 * fi_44[k]
                 - f_28 * fi_46[k]
                 + f_29 * fi_48[k]
                 - f_30 * fi_170[k]
                 - f_31 * fi_175[k]
                 + f_32 * fi_177[k]
                 - f_30 * fi_184[k]
                 + f_32 * fi_186[k]
                 - f_33 * fi_188[k];
    }

#pragma omp simd aligned(fi_28, fi_31, fi_33, fi_38, fi_42, fi_49, fi_51, fi_53, fi_168, \
                         fi_171, fi_173, fi_178, fi_182, fi_189, fi_191, \
                         fi_193 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_8[k] = f_45 * fi_28[k]
                 + f_45 * fi_31[k]
                 - f_18 * fi_33[k]
                 - f_45 * fi_38[k]
                 + f_18 * fi_42[k]
                 - f_45 * fi_49[k]
                 + f_18 * fi_51[k]
                 - f_18 * fi_53[k]
                 - f_46 * fi_168[k]
                 - f_46 * fi_171[k]
                 + f_21 * fi_173[k]
                 + f_46 * fi_178[k]
                 - f_21 * fi_182[k]
                 + f_46 * fi_189[k]
                 - f_21 * fi_191[k]
                 + f_21 * fi_193[k];
    }

#pragma omp simd aligned(fi_30, fi_35, fi_37, fi_44, fi_46, fi_170, fi_175, fi_177, fi_184, \
                         fi_186 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_9[k] = -f_17 * fi_30[k]
                 + f_15 * fi_35[k]
                 + f_18 * fi_37[k]
                 + f_14 * fi_44[k]
                 - f_16 * fi_46[k]
                 + f_20 * fi_170[k]
                 - f_19 * fi_175[k]
                 - f_21 * fi_177[k]
                 - f_17 * fi_184[k]
                 + f_18 * fi_186[k];
    }

#pragma omp simd aligned(fi_28, fi_31, fi_33, fi_38, fi_40, fi_49, fi_51, fi_168, fi_171, \
                         fi_173, fi_178, fi_180, fi_189, fi_191 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = -f_47 * fi_28[k]
                  + f_48 * fi_31[k]
                  + f_49 * fi_33[k]
                  + f_48 * fi_38[k]
                  - f_50 * fi_40[k]
                  - f_47 * fi_49[k]
                  + f_49 * fi_51[k]
                  + f_51 * fi_168[k]
                  - f_52 * fi_171[k]
                  - f_53 * fi_173[k]
                  - f_52 * fi_178[k]
                  + f_54 * fi_180[k]
                  + f_51 * fi_189[k]
                  - f_53 * fi_191[k];
    }

#pragma omp simd aligned(fi_28, fi_30, fi_31, fi_35, fi_38, fi_44, fi_49, fi_168, fi_170, \
                         fi_171, fi_175, fi_178, fi_184, fi_189 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_11[k] = f_6 * fi_30[k]
                  - f_5 * fi_35[k]
                  + f_4 * fi_44[k]
                  - f_9 * fi_170[k]
                  + f_8 * fi_175[k]
                  - f_7 * fi_184[k];

        g_12[k] = f_55 * fi_28[k]
                  - f_56 * fi_31[k]
                  + f_56 * fi_38[k]
                  - f_55 * fi_49[k]
                  - f_57 * fi_168[k]
                  + f_58 * fi_171[k]
                  - f_58 * fi_178[k]
                  + f_57 * fi_189[k];
    }

#pragma omp simd aligned(fi_113, fi_116, fi_118, fi_120, fi_123, fi_125, fi_127, fi_129, \
                         fi_131, fi_134, fi_136 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = f_59 * fi_113[k]
                  - f_60 * fi_118[k]
                  + f_59 * fi_127[k];

        g_14[k] = f_61 * fi_116[k]
                  - f_62 * fi_123[k]
                  + f_63 * fi_134[k];

        g_15[k] = -f_64 * fi_113[k]
                  + f_65 * fi_120[k]
                  + f_64 * fi_127[k]
                  - f_65 * fi_129[k];

        g_16[k] = -f_66 * fi_116[k]
                  - f_67 * fi_123[k]
                  + f_68 * fi_125[k]
                  + f_69 * fi_134[k]
                  - f_70 * fi_136[k];

        g_17[k] = f_71 * fi_113[k]
                  + f_72 * fi_118[k]
                  - f_73 * fi_120[k]
                  + f_71 * fi_127[k]
                  - f_73 * fi_129[k]
                  + f_73 * fi_131[k];
    }

#pragma omp simd aligned(fi_116, fi_123, fi_125, fi_134, fi_136, \
                         fi_138 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = f_74 * fi_116[k]
                  + f_75 * fi_123[k]
                  - f_76 * fi_125[k]
                  + f_74 * fi_134[k]
                  - f_76 * fi_136[k]
                  + f_77 * fi_138[k];
    }

#pragma omp simd aligned(fi_112, fi_115, fi_117, fi_122, fi_124, fi_126, fi_133, fi_135, \
                         fi_137, fi_139 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = -f_78 * fi_112[k]
                  - f_79 * fi_115[k]
                  + f_80 * fi_117[k]
                  - f_79 * fi_122[k]
                  + f_81 * fi_124[k]
                  - f_82 * fi_126[k]
                  - f_78 * fi_133[k]
                  + f_80 * fi_135[k]
                  - f_82 * fi_137[k]
                  + f_83 * fi_139[k];
    }

#pragma omp simd aligned(fi_112, fi_114, fi_115, fi_117, fi_119, fi_121, fi_122, fi_126, \
                         fi_128, fi_130, fi_132, fi_133, fi_135, \
                         fi_137 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = f_74 * fi_114[k]
                  + f_75 * fi_119[k]
                  - f_76 * fi_121[k]
                  + f_74 * fi_128[k]
                  - f_76 * fi_130[k]
                  + f_77 * fi_132[k];

        g_21[k] = f_84 * fi_112[k]
                  + f_84 * fi_115[k]
                  - f_70 * fi_117[k]
                  - f_84 * fi_122[k]
                  + f_70 * fi_126[k]
                  - f_84 * fi_133[k]
                  + f_70 * fi_135[k]
                  - f_70 * fi_137[k];
    }

#pragma omp simd aligned(fi_112, fi_114, fi_115, fi_117, fi_119, fi_121, fi_122, fi_124, \
                         fi_128, fi_130, fi_133, fi_135 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = -f_69 * fi_114[k]
                  + f_67 * fi_119[k]
                  + f_70 * fi_121[k]
                  + f_66 * fi_128[k]
                  - f_68 * fi_130[k];

        g_23[k] = -f_85 * fi_112[k]
                  + f_86 * fi_115[k]
                  + f_87 * fi_117[k]
                  + f_86 * fi_122[k]
                  - f_88 * fi_124[k]
                  - f_85 * fi_133[k]
                  + f_87 * fi_135[k];

        g_24[k] = f_63 * fi_114[k]
                  - f_62 * fi_119[k]
                  + f_61 * fi_128[k];

        g_25[k] = f_89 * fi_112[k]
                  - f_90 * fi_115[k]
                  + f_90 * fi_122[k]
                  - f_89 * fi_133[k];
    }

#pragma omp simd aligned(fi_29, fi_34, fi_43, fi_169, fi_174, fi_183, fi_225, fi_230, \
                         fi_239 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = -f_91 * fi_29[k]
                  + f_92 * fi_34[k]
                  - f_91 * fi_43[k]
                  - f_91 * fi_169[k]
                  + f_92 * fi_174[k]
                  - f_91 * fi_183[k]
                  + f_93 * fi_225[k]
                  - f_94 * fi_230[k]
                  + f_93 * fi_239[k];
    }

#pragma omp simd aligned(fi_32, fi_39, fi_50, fi_172, fi_179, fi_190, fi_228, fi_235, \
                         fi_246 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = -f_95 * fi_32[k]
                  + f_96 * fi_39[k]
                  - f_97 * fi_50[k]
                  - f_95 * fi_172[k]
                  + f_96 * fi_179[k]
                  - f_97 * fi_190[k]
                  + f_98 * fi_228[k]
                  - f_99 * fi_235[k]
                  + f_100 * fi_246[k];
    }

#pragma omp simd aligned(fi_29, fi_36, fi_43, fi_45, fi_169, fi_176, fi_183, fi_185, fi_225, \
                         fi_232, fi_239, fi_241 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = f_101 * fi_29[k]
                  - f_102 * fi_36[k]
                  - f_101 * fi_43[k]
                  + f_102 * fi_45[k]
                  + f_101 * fi_169[k]
                  - f_102 * fi_176[k]
                  - f_101 * fi_183[k]
                  + f_102 * fi_185[k]
                  - f_103 * fi_225[k]
                  + f_104 * fi_232[k]
                  + f_103 * fi_239[k]
                  - f_104 * fi_241[k];
    }

#pragma omp simd aligned(fi_32, fi_39, fi_41, fi_50, fi_52, fi_172, fi_179, fi_181, fi_190, \
                         fi_192, fi_228, fi_235, fi_237, fi_246, \
                         fi_248 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = f_105 * fi_32[k]
                  + f_106 * fi_39[k]
                  - f_107 * fi_41[k]
                  - f_108 * fi_50[k]
                  + f_109 * fi_52[k]
                  + f_105 * fi_172[k]
                  + f_106 * fi_179[k]
                  - f_107 * fi_181[k]
                  - f_108 * fi_190[k]
                  + f_109 * fi_192[k]
                  - f_110 * fi_228[k]
                  - f_107 * fi_235[k]
                  + f_111 * fi_237[k]
                  + f_112 * fi_246[k]
                  - f_77 * fi_248[k];
    }

#pragma omp simd aligned(fi_29, fi_34, fi_36, fi_43, fi_45, fi_47, fi_169, fi_174, fi_176, \
                         fi_183, fi_185, fi_187, fi_225, fi_230, fi_232, fi_239, fi_241, \
                         fi_243 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = -f_113 * fi_29[k]
                  - f_114 * fi_34[k]
                  + f_115 * fi_36[k]
                  - f_113 * fi_43[k]
                  + f_115 * fi_45[k]
                  - f_115 * fi_47[k]
                  - f_113 * fi_169[k]
                  - f_114 * fi_174[k]
                  + f_115 * fi_176[k]
                  - f_113 * fi_183[k]
                  + f_115 * fi_185[k]
                  - f_115 * fi_187[k]
                  + f_116 * fi_225[k]
                  + f_109 * fi_230[k]
                  - f_117 * fi_232[k]
                  + f_116 * fi_239[k]
                  - f_117 * fi_241[k]
                  + f_117 * fi_243[k];
    }

#pragma omp simd aligned(fi_32, fi_39, fi_41, fi_50, fi_52, fi_54, fi_172, fi_179, fi_181, \
                         fi_190, fi_192, fi_194, fi_228, fi_235, fi_237, fi_246, fi_248, \
                         fi_250 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = -f_84 * fi_32[k]
                  - f_71 * fi_39[k]
                  + f_72 * fi_41[k]
                  - f_84 * fi_50[k]
                  + f_72 * fi_52[k]
                  - f_118 * fi_54[k]
                  - f_84 * fi_172[k]
                  - f_71 * fi_179[k]
                  + f_72 * fi_181[k]
                  - f_84 * fi_190[k]
                  + f_72 * fi_192[k]
                  - f_118 * fi_194[k]
                  + f_72 * fi_228[k]
                  + f_119 * fi_235[k]
                  - f_70 * fi_237[k]
                  + f_72 * fi_246[k]
                  - f_70 * fi_248[k]
                  + f_120 * fi_250[k];
    }

#pragma omp simd aligned(fi_28, fi_31, fi_33, fi_38, fi_40, fi_42, fi_49, fi_51, fi_53, fi_55, \
                         fi_168, fi_171, fi_173, fi_178, fi_180, fi_182, fi_189, fi_191, \
                         fi_193, fi_195, fi_224, fi_227, fi_229, fi_234, fi_236, fi_238, \
                         fi_245, fi_247, fi_249, fi_251 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = f_121 * fi_28[k]
                  + f_122 * fi_31[k]
                  - f_123 * fi_33[k]
                  + f_122 * fi_38[k]
                  - f_124 * fi_40[k]
                  + f_125 * fi_42[k]
                  + f_121 * fi_49[k]
                  - f_123 * fi_51[k]
                  + f_125 * fi_53[k]
                  - f_126 * fi_55[k]
                  + f_121 * fi_168[k]
                  + f_122 * fi_171[k]
                  - f_123 * fi_173[k]
                  + f_122 * fi_178[k]
                  - f_124 * fi_180[k]
                  + f_125 * fi_182[k]
                  + f_121 * fi_189[k]
                  - f_123 * fi_191[k]
                  + f_125 * fi_193[k]
                  - f_126 * fi_195[k]
                  - f_127 * fi_224[k]
                  - f_128 * fi_227[k]
                  + f_129 * fi_229[k]
                  - f_128 * fi_234[k]
                  + f_130 * fi_236[k]
                  - f_131 * fi_238[k]
                  - f_127 * fi_245[k]
                  + f_129 * fi_247[k]
                  - f_131 * fi_249[k]
                  + f_132 * fi_251[k];
    }

#pragma omp simd aligned(fi_30, fi_35, fi_37, fi_44, fi_46, fi_48, fi_170, fi_175, fi_177, \
                         fi_184, fi_186, fi_188, fi_226, fi_231, fi_233, fi_240, fi_242, \
                         fi_244 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = -f_84 * fi_30[k]
                  - f_71 * fi_35[k]
                  + f_72 * fi_37[k]
                  - f_84 * fi_44[k]
                  + f_72 * fi_46[k]
                  - f_118 * fi_48[k]
                  - f_84 * fi_170[k]
                  - f_71 * fi_175[k]
                  + f_72 * fi_177[k]
                  - f_84 * fi_184[k]
                  + f_72 * fi_186[k]
                  - f_118 * fi_188[k]
                  + f_72 * fi_226[k]
                  + f_119 * fi_231[k]
                  - f_70 * fi_233[k]
                  + f_72 * fi_240[k]
                  - f_70 * fi_242[k]
                  + f_120 * fi_244[k];
    }

#pragma omp simd aligned(fi_28, fi_31, fi_33, fi_38, fi_42, fi_49, fi_51, fi_53, fi_168, \
                         fi_171, fi_173, fi_178, fi_182, fi_189, fi_191, fi_193, fi_224, \
                         fi_227, fi_229, fi_234, fi_238, fi_245, fi_247, \
                         fi_249 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = -f_133 * fi_28[k]
                  - f_133 * fi_31[k]
                  + f_109 * fi_33[k]
                  + f_133 * fi_38[k]
                  - f_109 * fi_42[k]
                  + f_133 * fi_49[k]
                  - f_109 * fi_51[k]
                  + f_109 * fi_53[k]
                  - f_133 * fi_168[k]
                  - f_133 * fi_171[k]
                  + f_109 * fi_173[k]
                  + f_133 * fi_178[k]
                  - f_109 * fi_182[k]
                  + f_133 * fi_189[k]
                  - f_109 * fi_191[k]
                  + f_109 * fi_193[k]
                  + f_114 * fi_224[k]
                  + f_114 * fi_227[k]
                  - f_77 * fi_229[k]
                  - f_114 * fi_234[k]
                  + f_77 * fi_238[k]
                  - f_114 * fi_245[k]
                  + f_77 * fi_247[k]
                  - f_77 * fi_249[k];
    }

#pragma omp simd aligned(fi_30, fi_35, fi_37, fi_44, fi_46, fi_170, fi_175, fi_177, fi_184, \
                         fi_186, fi_226, fi_231, fi_233, fi_240, \
                         fi_242 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = f_108 * fi_30[k]
                  - f_106 * fi_35[k]
                  - f_109 * fi_37[k]
                  - f_105 * fi_44[k]
                  + f_107 * fi_46[k]
                  + f_108 * fi_170[k]
                  - f_106 * fi_175[k]
                  - f_109 * fi_177[k]
                  - f_105 * fi_184[k]
                  + f_107 * fi_186[k]
                  - f_112 * fi_226[k]
                  + f_107 * fi_231[k]
                  + f_77 * fi_233[k]
                  + f_110 * fi_240[k]
                  - f_111 * fi_242[k];
    }

#pragma omp simd aligned(fi_28, fi_31, fi_33, fi_38, fi_40, fi_49, fi_51, fi_168, fi_171, \
                         fi_173, fi_178, fi_180, fi_189, fi_191, fi_224, fi_227, fi_229, \
                         fi_234, fi_236, fi_245, fi_247 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_36[k] = f_134 * fi_28[k]
                  - f_135 * fi_31[k]
                  - f_136 * fi_33[k]
                  - f_135 * fi_38[k]
                  + f_137 * fi_40[k]
                  + f_134 * fi_49[k]
                  - f_136 * fi_51[k]
                  + f_134 * fi_168[k]
                  - f_135 * fi_171[k]
                  - f_136 * fi_173[k]
                  - f_135 * fi_178[k]
                  + f_137 * fi_180[k]
                  + f_134 * fi_189[k]
                  - f_136 * fi_191[k]
                  - f_101 * fi_224[k]
                  + f_138 * fi_227[k]
                  + f_102 * fi_229[k]
                  + f_138 * fi_234[k]
                  - f_139 * fi_236[k]
                  - f_101 * fi_245[k]
                  + f_102 * fi_247[k];
    }

#pragma omp simd aligned(fi_30, fi_35, fi_44, fi_170, fi_175, fi_184, fi_226, fi_231, \
                         fi_240 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_37[k] = -f_97 * fi_30[k]
                  + f_96 * fi_35[k]
                  - f_95 * fi_44[k]
                  - f_97 * fi_170[k]
                  + f_96 * fi_175[k]
                  - f_95 * fi_184[k]
                  + f_100 * fi_226[k]
                  - f_99 * fi_231[k]
                  + f_98 * fi_240[k];
    }

#pragma omp simd aligned(fi_28, fi_31, fi_38, fi_49, fi_168, fi_171, fi_178, fi_189, fi_224, \
                         fi_227, fi_234, fi_245 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = -f_140 * fi_28[k]
                  + f_141 * fi_31[k]
                  - f_141 * fi_38[k]
                  + f_140 * fi_49[k]
                  - f_140 * fi_168[k]
                  + f_141 * fi_171[k]
                  - f_141 * fi_178[k]
                  + f_140 * fi_189[k]
                  + f_142 * fi_224[k]
                  - f_143 * fi_227[k]
                  + f_143 * fi_234[k]
                  - f_142 * fi_245[k];
    }

#pragma omp simd aligned(fi_57, fi_62, fi_71, fi_197, fi_202, fi_211, fi_253, fi_258, \
                         fi_267 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_39[k] = -f_144 * fi_57[k]
                  + f_145 * fi_62[k]
                  - f_144 * fi_71[k]
                  - f_144 * fi_197[k]
                  + f_145 * fi_202[k]
                  - f_144 * fi_211[k]
                  + f_146 * fi_253[k]
                  - f_147 * fi_258[k]
                  + f_146 * fi_267[k];
    }

#pragma omp simd aligned(fi_60, fi_67, fi_78, fi_200, fi_207, fi_218, fi_256, fi_263, \
                         fi_274 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_40[k] = -f_148 * fi_60[k]
                  + f_149 * fi_67[k]
                  - f_150 * fi_78[k]
                  - f_148 * fi_200[k]
                  + f_149 * fi_207[k]
                  - f_150 * fi_218[k]
                  + f_151 * fi_256[k]
                  - f_152 * fi_263[k]
                  + f_153 * fi_274[k];
    }

#pragma omp simd aligned(fi_57, fi_64, fi_71, fi_73, fi_197, fi_204, fi_211, fi_213, fi_253, \
                         fi_260, fi_267, fi_269 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_41[k] = f_154 * fi_57[k]
                  - f_155 * fi_64[k]
                  - f_154 * fi_71[k]
                  + f_155 * fi_73[k]
                  + f_154 * fi_197[k]
                  - f_155 * fi_204[k]
                  - f_154 * fi_211[k]
                  + f_155 * fi_213[k]
                  - f_156 * fi_253[k]
                  + f_157 * fi_260[k]
                  + f_156 * fi_267[k]
                  - f_157 * fi_269[k];
    }

#pragma omp simd aligned(fi_60, fi_67, fi_69, fi_78, fi_80, fi_200, fi_207, fi_209, fi_218, \
                         fi_220, fi_256, fi_263, fi_265, fi_274, \
                         fi_276 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = f_158 * fi_60[k]
                  + f_159 * fi_67[k]
                  - f_160 * fi_69[k]
                  - f_161 * fi_78[k]
                  + f_29 * fi_80[k]
                  + f_158 * fi_200[k]
                  + f_159 * fi_207[k]
                  - f_160 * fi_209[k]
                  - f_161 * fi_218[k]
                  + f_29 * fi_220[k]
                  - f_159 * fi_256[k]
                  - f_162 * fi_263[k]
                  + f_163 * fi_265[k]
                  + f_164 * fi_274[k]
                  - f_165 * fi_276[k];
    }

#pragma omp simd aligned(fi_57, fi_62, fi_64, fi_71, fi_73, fi_75, fi_197, fi_202, fi_204, \
                         fi_211, fi_213, fi_215, fi_253, fi_258, fi_260, fi_267, fi_269, \
                         fi_271 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_43[k] = -f_166 * fi_57[k]
                  - f_164 * fi_62[k]
                  + f_163 * fi_64[k]
                  - f_166 * fi_71[k]
                  + f_163 * fi_73[k]
                  - f_163 * fi_75[k]
                  - f_166 * fi_197[k]
                  - f_164 * fi_202[k]
                  + f_163 * fi_204[k]
                  - f_166 * fi_211[k]
                  + f_163 * fi_213[k]
                  - f_163 * fi_215[k]
                  + f_167 * fi_253[k]
                  + f_168 * fi_258[k]
                  - f_169 * fi_260[k]
                  + f_167 * fi_267[k]
                  - f_169 * fi_269[k]
                  + f_169 * fi_271[k];
    }

#pragma omp simd aligned(fi_60, fi_67, fi_69, fi_78, fi_80, fi_82, fi_200, fi_207, fi_209, \
                         fi_218, fi_220, fi_222, fi_256, fi_263, fi_265, fi_274, fi_276, \
                         fi_278 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_44[k] = -f_19 * fi_60[k]
                  - f_170 * fi_67[k]
                  + f_18 * fi_69[k]
                  - f_19 * fi_78[k]
                  + f_18 * fi_80[k]
                  - f_171 * fi_82[k]
                  - f_19 * fi_200[k]
                  - f_170 * fi_207[k]
                  + f_18 * fi_209[k]
                  - f_19 * fi_218[k]
                  + f_18 * fi_220[k]
                  - f_171 * fi_222[k]
                  + f_172 * fi_256[k]
                  + f_21 * fi_263[k]
                  - f_25 * fi_265[k]
                  + f_172 * fi_274[k]
                  - f_25 * fi_276[k]
                  + f_173 * fi_278[k];
    }

#pragma omp simd aligned(fi_56, fi_59, fi_61, fi_66, fi_68, fi_70, fi_77, fi_79, fi_81, fi_83, \
                         fi_196, fi_199, fi_201, fi_206, fi_208, fi_210, fi_217, fi_219, \
                         fi_221, fi_223, fi_252, fi_255, fi_257, fi_262, fi_264, fi_266, \
                         fi_273, fi_275, fi_277, fi_279 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_45[k] = 0.46875 * fi_56[k]
                  + 1.40625 * fi_59[k]
                  - 8.4375 * fi_61[k]
                  + 1.40625 * fi_66[k]
                  - 16.875 * fi_68[k]
                  + 11.25 * fi_70[k]
                  + 0.46875 * fi_77[k]
                  - 8.4375 * fi_79[k]
                  + 11.25 * fi_81[k]
                  - 1.5 * fi_83[k]
                  + 0.46875 * fi_196[k]
                  + 1.40625 * fi_199[k]
                  - 8.4375 * fi_201[k]
                  + 1.40625 * fi_206[k]
                  - 16.875 * fi_208[k]
                  + 11.25 * fi_210[k]
                  + 0.46875 * fi_217[k]
                  - 8.4375 * fi_219[k]
                  + 11.25 * fi_221[k]
                  - 1.5 * fi_223[k]
                  - 0.3125 * fi_252[k]
                  - 0.9375 * fi_255[k]
                  + 5.625 * fi_257[k]
                  - 0.9375 * fi_262[k]
                  + 11.25 * fi_264[k]
                  - 7.5 * fi_266[k]
                  - 0.3125 * fi_273[k]
                  + 5.625 * fi_275[k]
                  - 7.5 * fi_277[k]
                  + fi_279[k];
    }

#pragma omp simd aligned(fi_58, fi_63, fi_65, fi_72, fi_74, fi_76, fi_198, fi_203, fi_205, \
                         fi_212, fi_214, fi_216, fi_254, fi_259, fi_261, fi_268, fi_270, \
                         fi_272 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_46[k] = -f_19 * fi_58[k]
                  - f_170 * fi_63[k]
                  + f_18 * fi_65[k]
                  - f_19 * fi_72[k]
                  + f_18 * fi_74[k]
                  - f_171 * fi_76[k]
                  - f_19 * fi_198[k]
                  - f_170 * fi_203[k]
                  + f_18 * fi_205[k]
                  - f_19 * fi_212[k]
                  + f_18 * fi_214[k]
                  - f_171 * fi_216[k]
                  + f_172 * fi_254[k]
                  + f_21 * fi_259[k]
                  - f_25 * fi_261[k]
                  + f_172 * fi_268[k]
                  - f_25 * fi_270[k]
                  + f_173 * fi_272[k];
    }

#pragma omp simd aligned(fi_56, fi_59, fi_61, fi_66, fi_70, fi_77, fi_79, fi_81, fi_196, \
                         fi_199, fi_201, fi_206, fi_210, fi_217, fi_219, fi_221, fi_252, \
                         fi_255, fi_257, fi_262, fi_266, fi_273, fi_275, \
                         fi_277 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_47[k] = -f_174 * fi_56[k]
                  - f_174 * fi_59[k]
                  + f_29 * fi_61[k]
                  + f_174 * fi_66[k]
                  - f_29 * fi_70[k]
                  + f_174 * fi_77[k]
                  - f_29 * fi_79[k]
                  + f_29 * fi_81[k]
                  - f_174 * fi_196[k]
                  - f_174 * fi_199[k]
                  + f_29 * fi_201[k]
                  + f_174 * fi_206[k]
                  - f_29 * fi_210[k]
                  + f_174 * fi_217[k]
                  - f_29 * fi_219[k]
                  + f_29 * fi_221[k]
                  + f_175 * fi_252[k]
                  + f_175 * fi_255[k]
                  - f_165 * fi_257[k]
                  - f_175 * fi_262[k]
                  + f_165 * fi_266[k]
                  - f_175 * fi_273[k]
                  + f_165 * fi_275[k]
                  - f_165 * fi_277[k];
    }

#pragma omp simd aligned(fi_58, fi_63, fi_65, fi_72, fi_74, fi_198, fi_203, fi_205, fi_212, \
                         fi_214, fi_254, fi_259, fi_261, fi_268, \
                         fi_270 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_48[k] = f_161 * fi_58[k]
                  - f_159 * fi_63[k]
                  - f_29 * fi_65[k]
                  - f_158 * fi_72[k]
                  + f_160 * fi_74[k]
                  + f_161 * fi_198[k]
                  - f_159 * fi_203[k]
                  - f_29 * fi_205[k]
                  - f_158 * fi_212[k]
                  + f_160 * fi_214[k]
                  - f_164 * fi_254[k]
                  + f_162 * fi_259[k]
                  + f_165 * fi_261[k]
                  + f_159 * fi_268[k]
                  - f_163 * fi_270[k];
    }

#pragma omp simd aligned(fi_56, fi_59, fi_61, fi_66, fi_68, fi_77, fi_79, fi_196, fi_199, \
                         fi_201, fi_206, fi_208, fi_217, fi_219, fi_252, fi_255, fi_257, \
                         fi_262, fi_264, fi_273, fi_275 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_49[k] = f_176 * fi_56[k]
                  - f_177 * fi_59[k]
                  - f_178 * fi_61[k]
                  - f_177 * fi_66[k]
                  + f_179 * fi_68[k]
                  + f_176 * fi_77[k]
                  - f_178 * fi_79[k]
                  + f_176 * fi_196[k]
                  - f_177 * fi_199[k]
                  - f_178 * fi_201[k]
                  - f_177 * fi_206[k]
                  + f_179 * fi_208[k]
                  + f_176 * fi_217[k]
                  - f_178 * fi_219[k]
                  - f_180 * fi_252[k]
                  + f_181 * fi_255[k]
                  + f_182 * fi_257[k]
                  + f_181 * fi_262[k]
                  - f_155 * fi_264[k]
                  - f_180 * fi_273[k]
                  + f_182 * fi_275[k];
    }

#pragma omp simd aligned(fi_58, fi_63, fi_72, fi_198, fi_203, fi_212, fi_254, fi_259, \
                         fi_268 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_50[k] = -f_150 * fi_58[k]
                  + f_149 * fi_63[k]
                  - f_148 * fi_72[k]
                  - f_150 * fi_198[k]
                  + f_149 * fi_203[k]
                  - f_148 * fi_212[k]
                  + f_153 * fi_254[k]
                  - f_152 * fi_259[k]
                  + f_151 * fi_268[k];
    }

#pragma omp simd aligned(fi_56, fi_59, fi_66, fi_77, fi_196, fi_199, fi_206, fi_217, fi_252, \
                         fi_255, fi_262, fi_273 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_51[k] = -f_183 * fi_56[k]
                  + f_184 * fi_59[k]
                  - f_184 * fi_66[k]
                  + f_183 * fi_77[k]
                  - f_183 * fi_196[k]
                  + f_184 * fi_199[k]
                  - f_184 * fi_206[k]
                  + f_183 * fi_217[k]
                  + f_185 * fi_252[k]
                  - f_186 * fi_255[k]
                  + f_186 * fi_262[k]
                  - f_185 * fi_273[k];
    }

#pragma omp simd aligned(fi_1, fi_6, fi_15, fi_85, fi_90, fi_99, fi_141, fi_146, \
                         fi_155 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_52[k] = -f_91 * fi_1[k]
                  + f_92 * fi_6[k]
                  - f_91 * fi_15[k]
                  - f_91 * fi_85[k]
                  + f_92 * fi_90[k]
                  - f_91 * fi_99[k]
                  + f_93 * fi_141[k]
                  - f_94 * fi_146[k]
                  + f_93 * fi_155[k];
    }

#pragma omp simd aligned(fi_4, fi_11, fi_22, fi_88, fi_95, fi_106, fi_144, fi_151, \
                         fi_162 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_53[k] = -f_95 * fi_4[k]
                  + f_96 * fi_11[k]
                  - f_97 * fi_22[k]
                  - f_95 * fi_88[k]
                  + f_96 * fi_95[k]
                  - f_97 * fi_106[k]
                  + f_98 * fi_144[k]
                  - f_99 * fi_151[k]
                  + f_100 * fi_162[k];
    }

#pragma omp simd aligned(fi_1, fi_8, fi_15, fi_17, fi_85, fi_92, fi_99, fi_101, fi_141, \
                         fi_148, fi_155, fi_157 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_54[k] = f_101 * fi_1[k]
                  - f_102 * fi_8[k]
                  - f_101 * fi_15[k]
                  + f_102 * fi_17[k]
                  + f_101 * fi_85[k]
                  - f_102 * fi_92[k]
                  - f_101 * fi_99[k]
                  + f_102 * fi_101[k]
                  - f_103 * fi_141[k]
                  + f_104 * fi_148[k]
                  + f_103 * fi_155[k]
                  - f_104 * fi_157[k];
    }

#pragma omp simd aligned(fi_4, fi_11, fi_13, fi_22, fi_24, fi_88, fi_95, fi_97, fi_106, \
                         fi_108, fi_144, fi_151, fi_153, fi_162, \
                         fi_164 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_55[k] = f_105 * fi_4[k]
                  + f_106 * fi_11[k]
                  - f_107 * fi_13[k]
                  - f_108 * fi_22[k]
                  + f_109 * fi_24[k]
                  + f_105 * fi_88[k]
                  + f_106 * fi_95[k]
                  - f_107 * fi_97[k]
                  - f_108 * fi_106[k]
                  + f_109 * fi_108[k]
                  - f_110 * fi_144[k]
                  - f_107 * fi_151[k]
                  + f_111 * fi_153[k]
                  + f_112 * fi_162[k]
                  - f_77 * fi_164[k];
    }

#pragma omp simd aligned(fi_1, fi_6, fi_8, fi_15, fi_17, fi_19, fi_85, fi_90, fi_92, fi_99, \
                         fi_101, fi_103, fi_141, fi_146, fi_148, fi_155, fi_157, \
                         fi_159 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_56[k] = -f_113 * fi_1[k]
                  - f_114 * fi_6[k]
                  + f_115 * fi_8[k]
                  - f_113 * fi_15[k]
                  + f_115 * fi_17[k]
                  - f_115 * fi_19[k]
                  - f_113 * fi_85[k]
                  - f_114 * fi_90[k]
                  + f_115 * fi_92[k]
                  - f_113 * fi_99[k]
                  + f_115 * fi_101[k]
                  - f_115 * fi_103[k]
                  + f_116 * fi_141[k]
                  + f_109 * fi_146[k]
                  - f_117 * fi_148[k]
                  + f_116 * fi_155[k]
                  - f_117 * fi_157[k]
                  + f_117 * fi_159[k];
    }

#pragma omp simd aligned(fi_4, fi_11, fi_13, fi_22, fi_24, fi_26, fi_88, fi_95, fi_97, fi_106, \
                         fi_108, fi_110, fi_144, fi_151, fi_153, fi_162, fi_164, \
                         fi_166 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_57[k] = -f_84 * fi_4[k]
                  - f_71 * fi_11[k]
                  + f_72 * fi_13[k]
                  - f_84 * fi_22[k]
                  + f_72 * fi_24[k]
                  - f_118 * fi_26[k]
                  - f_84 * fi_88[k]
                  - f_71 * fi_95[k]
                  + f_72 * fi_97[k]
                  - f_84 * fi_106[k]
                  + f_72 * fi_108[k]
                  - f_118 * fi_110[k]
                  + f_72 * fi_144[k]
                  + f_119 * fi_151[k]
                  - f_70 * fi_153[k]
                  + f_72 * fi_162[k]
                  - f_70 * fi_164[k]
                  + f_120 * fi_166[k];
    }

#pragma omp simd aligned(fi_0, fi_3, fi_5, fi_10, fi_12, fi_14, fi_21, fi_23, fi_25, fi_27, \
                         fi_84, fi_87, fi_89, fi_94, fi_96, fi_98, fi_105, fi_107, fi_109, \
                         fi_111, fi_140, fi_143, fi_145, fi_150, fi_152, fi_154, fi_161, \
                         fi_163, fi_165, fi_167 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_58[k] = f_121 * fi_0[k]
                  + f_122 * fi_3[k]
                  - f_123 * fi_5[k]
                  + f_122 * fi_10[k]
                  - f_124 * fi_12[k]
                  + f_125 * fi_14[k]
                  + f_121 * fi_21[k]
                  - f_123 * fi_23[k]
                  + f_125 * fi_25[k]
                  - f_126 * fi_27[k]
                  + f_121 * fi_84[k]
                  + f_122 * fi_87[k]
                  - f_123 * fi_89[k]
                  + f_122 * fi_94[k]
                  - f_124 * fi_96[k]
                  + f_125 * fi_98[k]
                  + f_121 * fi_105[k]
                  - f_123 * fi_107[k]
                  + f_125 * fi_109[k]
                  - f_126 * fi_111[k]
                  - f_127 * fi_140[k]
                  - f_128 * fi_143[k]
                  + f_129 * fi_145[k]
                  - f_128 * fi_150[k]
                  + f_130 * fi_152[k]
                  - f_131 * fi_154[k]
                  - f_127 * fi_161[k]
                  + f_129 * fi_163[k]
                  - f_131 * fi_165[k]
                  + f_132 * fi_167[k];
    }

#pragma omp simd aligned(fi_2, fi_7, fi_9, fi_16, fi_18, fi_20, fi_86, fi_91, fi_93, fi_100, \
                         fi_102, fi_104, fi_142, fi_147, fi_149, fi_156, fi_158, \
                         fi_160 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_59[k] = -f_84 * fi_2[k]
                  - f_71 * fi_7[k]
                  + f_72 * fi_9[k]
                  - f_84 * fi_16[k]
                  + f_72 * fi_18[k]
                  - f_118 * fi_20[k]
                  - f_84 * fi_86[k]
                  - f_71 * fi_91[k]
                  + f_72 * fi_93[k]
                  - f_84 * fi_100[k]
                  + f_72 * fi_102[k]
                  - f_118 * fi_104[k]
                  + f_72 * fi_142[k]
                  + f_119 * fi_147[k]
                  - f_70 * fi_149[k]
                  + f_72 * fi_156[k]
                  - f_70 * fi_158[k]
                  + f_120 * fi_160[k];
    }

#pragma omp simd aligned(fi_0, fi_3, fi_5, fi_10, fi_14, fi_21, fi_23, fi_25, fi_84, fi_87, \
                         fi_89, fi_94, fi_98, fi_105, fi_107, fi_109, fi_140, fi_143, fi_145, \
                         fi_150, fi_154, fi_161, fi_163, fi_165 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_60[k] = -f_133 * fi_0[k]
                  - f_133 * fi_3[k]
                  + f_109 * fi_5[k]
                  + f_133 * fi_10[k]
                  - f_109 * fi_14[k]
                  + f_133 * fi_21[k]
                  - f_109 * fi_23[k]
                  + f_109 * fi_25[k]
                  - f_133 * fi_84[k]
                  - f_133 * fi_87[k]
                  + f_109 * fi_89[k]
                  + f_133 * fi_94[k]
                  - f_109 * fi_98[k]
                  + f_133 * fi_105[k]
                  - f_109 * fi_107[k]
                  + f_109 * fi_109[k]
                  + f_114 * fi_140[k]
                  + f_114 * fi_143[k]
                  - f_77 * fi_145[k]
                  - f_114 * fi_150[k]
                  + f_77 * fi_154[k]
                  - f_114 * fi_161[k]
                  + f_77 * fi_163[k]
                  - f_77 * fi_165[k];
    }

#pragma omp simd aligned(fi_2, fi_7, fi_9, fi_16, fi_18, fi_86, fi_91, fi_93, fi_100, fi_102, \
                         fi_142, fi_147, fi_149, fi_156, fi_158 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_61[k] = f_108 * fi_2[k]
                  - f_106 * fi_7[k]
                  - f_109 * fi_9[k]
                  - f_105 * fi_16[k]
                  + f_107 * fi_18[k]
                  + f_108 * fi_86[k]
                  - f_106 * fi_91[k]
                  - f_109 * fi_93[k]
                  - f_105 * fi_100[k]
                  + f_107 * fi_102[k]
                  - f_112 * fi_142[k]
                  + f_107 * fi_147[k]
                  + f_77 * fi_149[k]
                  + f_110 * fi_156[k]
                  - f_111 * fi_158[k];
    }

#pragma omp simd aligned(fi_0, fi_3, fi_5, fi_10, fi_12, fi_21, fi_23, fi_84, fi_87, fi_89, \
                         fi_94, fi_96, fi_105, fi_107, fi_140, fi_143, fi_145, fi_150, fi_152, \
                         fi_161, fi_163 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_62[k] = f_134 * fi_0[k]
                  - f_135 * fi_3[k]
                  - f_136 * fi_5[k]
                  - f_135 * fi_10[k]
                  + f_137 * fi_12[k]
                  + f_134 * fi_21[k]
                  - f_136 * fi_23[k]
                  + f_134 * fi_84[k]
                  - f_135 * fi_87[k]
                  - f_136 * fi_89[k]
                  - f_135 * fi_94[k]
                  + f_137 * fi_96[k]
                  + f_134 * fi_105[k]
                  - f_136 * fi_107[k]
                  - f_101 * fi_140[k]
                  + f_138 * fi_143[k]
                  + f_102 * fi_145[k]
                  + f_138 * fi_150[k]
                  - f_139 * fi_152[k]
                  - f_101 * fi_161[k]
                  + f_102 * fi_163[k];
    }

#pragma omp simd aligned(fi_2, fi_7, fi_16, fi_86, fi_91, fi_100, fi_142, fi_147, \
                         fi_156 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_63[k] = -f_97 * fi_2[k]
                  + f_96 * fi_7[k]
                  - f_95 * fi_16[k]
                  - f_97 * fi_86[k]
                  + f_96 * fi_91[k]
                  - f_95 * fi_100[k]
                  + f_100 * fi_142[k]
                  - f_99 * fi_147[k]
                  + f_98 * fi_156[k];
    }

#pragma omp simd aligned(fi_0, fi_3, fi_10, fi_21, fi_84, fi_87, fi_94, fi_105, fi_140, \
                         fi_143, fi_150, fi_161 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_64[k] = -f_140 * fi_0[k]
                  + f_141 * fi_3[k]
                  - f_141 * fi_10[k]
                  + f_140 * fi_21[k]
                  - f_140 * fi_84[k]
                  + f_141 * fi_87[k]
                  - f_141 * fi_94[k]
                  + f_140 * fi_105[k]
                  + f_142 * fi_140[k]
                  - f_143 * fi_143[k]
                  + f_143 * fi_150[k]
                  - f_142 * fi_161[k];
    }

#pragma omp simd aligned(fi_57, fi_60, fi_62, fi_67, fi_71, fi_78, fi_197, fi_200, fi_202, \
                         fi_207, fi_211, fi_218 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_65[k] = f_187 * fi_57[k]
                  - f_188 * fi_62[k]
                  + f_187 * fi_71[k]
                  - f_187 * fi_197[k]
                  + f_188 * fi_202[k]
                  - f_187 * fi_211[k];

        g_66[k] = f_189 * fi_60[k]
                  - f_61 * fi_67[k]
                  + f_190 * fi_78[k]
                  - f_189 * fi_200[k]
                  + f_61 * fi_207[k]
                  - f_190 * fi_218[k];
    }

#pragma omp simd aligned(fi_57, fi_64, fi_71, fi_73, fi_197, fi_204, fi_211, \
                         fi_213 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_67[k] = -f_191 * fi_57[k]
                  + f_192 * fi_64[k]
                  + f_191 * fi_71[k]
                  - f_192 * fi_73[k]
                  + f_191 * fi_197[k]
                  - f_192 * fi_204[k]
                  - f_191 * fi_211[k]
                  + f_192 * fi_213[k];
    }

#pragma omp simd aligned(fi_60, fi_67, fi_69, fi_78, fi_80, fi_200, fi_207, fi_209, fi_218, \
                         fi_220 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_68[k] = -f_193 * fi_60[k]
                  - f_69 * fi_67[k]
                  + f_194 * fi_69[k]
                  + f_195 * fi_78[k]
                  - f_119 * fi_80[k]
                  + f_193 * fi_200[k]
                  + f_69 * fi_207[k]
                  - f_194 * fi_209[k]
                  - f_195 * fi_218[k]
                  + f_119 * fi_220[k];
    }

#pragma omp simd aligned(fi_57, fi_62, fi_64, fi_71, fi_73, fi_75, fi_197, fi_202, fi_204, \
                         fi_211, fi_213, fi_215 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_69[k] = f_84 * fi_57[k]
                  + f_71 * fi_62[k]
                  - f_70 * fi_64[k]
                  + f_84 * fi_71[k]
                  - f_70 * fi_73[k]
                  + f_70 * fi_75[k]
                  - f_84 * fi_197[k]
                  - f_71 * fi_202[k]
                  + f_70 * fi_204[k]
                  - f_84 * fi_211[k]
                  + f_70 * fi_213[k]
                  - f_70 * fi_215[k];
    }

#pragma omp simd aligned(fi_60, fi_67, fi_69, fi_78, fi_80, fi_82, fi_200, fi_207, fi_209, \
                         fi_218, fi_220, fi_222 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_70[k] = f_196 * fi_60[k]
                  + f_74 * fi_67[k]
                  - f_75 * fi_69[k]
                  + f_196 * fi_78[k]
                  - f_75 * fi_80[k]
                  + f_115 * fi_82[k]
                  - f_196 * fi_200[k]
                  - f_74 * fi_207[k]
                  + f_75 * fi_209[k]
                  - f_196 * fi_218[k]
                  + f_75 * fi_220[k]
                  - f_115 * fi_222[k];
    }

#pragma omp simd aligned(fi_56, fi_59, fi_61, fi_66, fi_68, fi_70, fi_77, fi_79, fi_81, fi_83, \
                         fi_196, fi_199, fi_201, fi_206, fi_208, fi_210, fi_217, fi_219, \
                         fi_221, fi_223 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_71[k] = -f_197 * fi_56[k]
                  - f_198 * fi_59[k]
                  + f_199 * fi_61[k]
                  - f_198 * fi_66[k]
                  + f_80 * fi_68[k]
                  - f_200 * fi_70[k]
                  - f_197 * fi_77[k]
                  + f_199 * fi_79[k]
                  - f_200 * fi_81[k]
                  + f_201 * fi_83[k]
                  + f_197 * fi_196[k]
                  + f_198 * fi_199[k]
                  - f_199 * fi_201[k]
                  + f_198 * fi_206[k]
                  - f_80 * fi_208[k]
                  + f_200 * fi_210[k]
                  + f_197 * fi_217[k]
                  - f_199 * fi_219[k]
                  + f_200 * fi_221[k]
                  - f_201 * fi_223[k];
    }

#pragma omp simd aligned(fi_58, fi_63, fi_65, fi_72, fi_74, fi_76, fi_198, fi_203, fi_205, \
                         fi_212, fi_214, fi_216 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_72[k] = f_196 * fi_58[k]
                  + f_74 * fi_63[k]
                  - f_75 * fi_65[k]
                  + f_196 * fi_72[k]
                  - f_75 * fi_74[k]
                  + f_115 * fi_76[k]
                  - f_196 * fi_198[k]
                  - f_74 * fi_203[k]
                  + f_75 * fi_205[k]
                  - f_196 * fi_212[k]
                  + f_75 * fi_214[k]
                  - f_115 * fi_216[k];
    }

#pragma omp simd aligned(fi_56, fi_59, fi_61, fi_66, fi_70, fi_77, fi_79, fi_81, fi_196, \
                         fi_199, fi_201, fi_206, fi_210, fi_217, fi_219, \
                         fi_221 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_73[k] = f_202 * fi_56[k]
                  + f_202 * fi_59[k]
                  - f_119 * fi_61[k]
                  - f_202 * fi_66[k]
                  + f_119 * fi_70[k]
                  - f_202 * fi_77[k]
                  + f_119 * fi_79[k]
                  - f_119 * fi_81[k]
                  - f_202 * fi_196[k]
                  - f_202 * fi_199[k]
                  + f_119 * fi_201[k]
                  + f_202 * fi_206[k]
                  - f_119 * fi_210[k]
                  + f_202 * fi_217[k]
                  - f_119 * fi_219[k]
                  + f_119 * fi_221[k];
    }

#pragma omp simd aligned(fi_58, fi_63, fi_65, fi_72, fi_74, fi_198, fi_203, fi_205, fi_212, \
                         fi_214 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_74[k] = -f_195 * fi_58[k]
                  + f_69 * fi_63[k]
                  + f_119 * fi_65[k]
                  + f_193 * fi_72[k]
                  - f_194 * fi_74[k]
                  + f_195 * fi_198[k]
                  - f_69 * fi_203[k]
                  - f_119 * fi_205[k]
                  - f_193 * fi_212[k]
                  + f_194 * fi_214[k];
    }

#pragma omp simd aligned(fi_56, fi_59, fi_61, fi_66, fi_68, fi_77, fi_79, fi_196, fi_199, \
                         fi_201, fi_206, fi_208, fi_217, fi_219 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_75[k] = -f_203 * fi_56[k]
                  + f_204 * fi_59[k]
                  + f_86 * fi_61[k]
                  + f_204 * fi_66[k]
                  - f_205 * fi_68[k]
                  - f_203 * fi_77[k]
                  + f_86 * fi_79[k]
                  + f_203 * fi_196[k]
                  - f_204 * fi_199[k]
                  - f_86 * fi_201[k]
                  - f_204 * fi_206[k]
                  + f_205 * fi_208[k]
                  + f_203 * fi_217[k]
                  - f_86 * fi_219[k];
    }

#pragma omp simd aligned(fi_56, fi_58, fi_59, fi_63, fi_66, fi_72, fi_77, fi_196, fi_198, \
                         fi_199, fi_203, fi_206, fi_212, fi_217 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_76[k] = f_190 * fi_58[k]
                  - f_61 * fi_63[k]
                  + f_189 * fi_72[k]
                  - f_190 * fi_198[k]
                  + f_61 * fi_203[k]
                  - f_189 * fi_212[k];

        g_77[k] = f_206 * fi_56[k]
                  - f_207 * fi_59[k]
                  + f_207 * fi_66[k]
                  - f_206 * fi_77[k]
                  - f_206 * fi_196[k]
                  + f_207 * fi_199[k]
                  - f_207 * fi_206[k]
                  + f_206 * fi_217[k];
    }

#pragma omp simd aligned(fi_1, fi_4, fi_6, fi_11, fi_15, fi_22, fi_85, fi_88, fi_90, fi_95, \
                         fi_99, fi_106 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_78[k] = f_2 * fi_1[k]
                  - f_3 * fi_6[k]
                  + f_2 * fi_15[k]
                  - f_0 * fi_85[k]
                  + f_1 * fi_90[k]
                  - f_0 * fi_99[k];

        g_79[k] = f_7 * fi_4[k]
                  - f_8 * fi_11[k]
                  + f_9 * fi_22[k]
                  - f_4 * fi_88[k]
                  + f_5 * fi_95[k]
                  - f_6 * fi_106[k];
    }

#pragma omp simd aligned(fi_1, fi_8, fi_15, fi_17, fi_85, fi_92, fi_99, \
                         fi_101 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_80[k] = -f_12 * fi_1[k]
                  + f_13 * fi_8[k]
                  + f_12 * fi_15[k]
                  - f_13 * fi_17[k]
                  + f_10 * fi_85[k]
                  - f_11 * fi_92[k]
                  - f_10 * fi_99[k]
                  + f_11 * fi_101[k];
    }

#pragma omp simd aligned(fi_4, fi_11, fi_13, fi_22, fi_24, fi_88, fi_95, fi_97, fi_106, \
                         fi_108 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_81[k] = -f_17 * fi_4[k]
                  - f_19 * fi_11[k]
                  + f_18 * fi_13[k]
                  + f_20 * fi_22[k]
                  - f_21 * fi_24[k]
                  + f_14 * fi_88[k]
                  + f_15 * fi_95[k]
                  - f_16 * fi_97[k]
                  - f_17 * fi_106[k]
                  + f_18 * fi_108[k];
    }

#pragma omp simd aligned(fi_1, fi_6, fi_8, fi_15, fi_17, fi_19, fi_85, fi_90, fi_92, fi_99, \
                         fi_101, fi_103 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_82[k] = f_23 * fi_1[k]
                  + f_24 * fi_6[k]
                  - f_25 * fi_8[k]
                  + f_23 * fi_15[k]
                  - f_25 * fi_17[k]
                  + f_25 * fi_19[k]
                  - f_20 * fi_85[k]
                  - f_19 * fi_90[k]
                  + f_22 * fi_92[k]
                  - f_20 * fi_99[k]
                  + f_22 * fi_101[k]
                  - f_22 * fi_103[k];
    }

#pragma omp simd aligned(fi_4, fi_11, fi_13, fi_22, fi_24, fi_26, fi_88, fi_95, fi_97, fi_106, \
                         fi_108, fi_110 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_83[k] = f_30 * fi_4[k]
                  + f_31 * fi_11[k]
                  - f_32 * fi_13[k]
                  + f_30 * fi_22[k]
                  - f_32 * fi_24[k]
                  + f_33 * fi_26[k]
                  - f_26 * fi_88[k]
                  - f_27 * fi_95[k]
                  + f_28 * fi_97[k]
                  - f_26 * fi_106[k]
                  + f_28 * fi_108[k]
                  - f_29 * fi_110[k];
    }

#pragma omp simd aligned(fi_0, fi_3, fi_5, fi_10, fi_12, fi_14, fi_21, fi_23, fi_25, fi_27, \
                         fi_84, fi_87, fi_89, fi_94, fi_96, fi_98, fi_105, fi_107, fi_109, \
                         fi_111 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_84[k] = -f_40 * fi_0[k]
                  - f_34 * fi_3[k]
                  + f_41 * fi_5[k]
                  - f_34 * fi_10[k]
                  + f_42 * fi_12[k]
                  - f_43 * fi_14[k]
                  - f_40 * fi_21[k]
                  + f_41 * fi_23[k]
                  - f_43 * fi_25[k]
                  + f_44 * fi_27[k]
                  + f_34 * fi_84[k]
                  + f_35 * fi_87[k]
                  - f_36 * fi_89[k]
                  + f_35 * fi_94[k]
                  - f_37 * fi_96[k]
                  + f_38 * fi_98[k]
                  + f_34 * fi_105[k]
                  - f_36 * fi_107[k]
                  + f_38 * fi_109[k]
                  - f_39 * fi_111[k];
    }

#pragma omp simd aligned(fi_2, fi_7, fi_9, fi_16, fi_18, fi_20, fi_86, fi_91, fi_93, fi_100, \
                         fi_102, fi_104 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_85[k] = f_30 * fi_2[k]
                  + f_31 * fi_7[k]
                  - f_32 * fi_9[k]
                  + f_30 * fi_16[k]
                  - f_32 * fi_18[k]
                  + f_33 * fi_20[k]
                  - f_26 * fi_86[k]
                  - f_27 * fi_91[k]
                  + f_28 * fi_93[k]
                  - f_26 * fi_100[k]
                  + f_28 * fi_102[k]
                  - f_29 * fi_104[k];
    }

#pragma omp simd aligned(fi_0, fi_3, fi_5, fi_10, fi_14, fi_21, fi_23, fi_25, fi_84, fi_87, \
                         fi_89, fi_94, fi_98, fi_105, fi_107, fi_109 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_86[k] = f_46 * fi_0[k]
                  + f_46 * fi_3[k]
                  - f_21 * fi_5[k]
                  - f_46 * fi_10[k]
                  + f_21 * fi_14[k]
                  - f_46 * fi_21[k]
                  + f_21 * fi_23[k]
                  - f_21 * fi_25[k]
                  - f_45 * fi_84[k]
                  - f_45 * fi_87[k]
                  + f_18 * fi_89[k]
                  + f_45 * fi_94[k]
                  - f_18 * fi_98[k]
                  + f_45 * fi_105[k]
                  - f_18 * fi_107[k]
                  + f_18 * fi_109[k];
    }

#pragma omp simd aligned(fi_2, fi_7, fi_9, fi_16, fi_18, fi_86, fi_91, fi_93, fi_100, \
                         fi_102 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_87[k] = -f_20 * fi_2[k]
                  + f_19 * fi_7[k]
                  + f_21 * fi_9[k]
                  + f_17 * fi_16[k]
                  - f_18 * fi_18[k]
                  + f_17 * fi_86[k]
                  - f_15 * fi_91[k]
                  - f_18 * fi_93[k]
                  - f_14 * fi_100[k]
                  + f_16 * fi_102[k];
    }

#pragma omp simd aligned(fi_0, fi_3, fi_5, fi_10, fi_12, fi_21, fi_23, fi_84, fi_87, fi_89, \
                         fi_94, fi_96, fi_105, fi_107 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_88[k] = -f_51 * fi_0[k]
                  + f_52 * fi_3[k]
                  + f_53 * fi_5[k]
                  + f_52 * fi_10[k]
                  - f_54 * fi_12[k]
                  - f_51 * fi_21[k]
                  + f_53 * fi_23[k]
                  + f_47 * fi_84[k]
                  - f_48 * fi_87[k]
                  - f_49 * fi_89[k]
                  - f_48 * fi_94[k]
                  + f_50 * fi_96[k]
                  + f_47 * fi_105[k]
                  - f_49 * fi_107[k];
    }

#pragma omp simd aligned(fi_0, fi_2, fi_3, fi_7, fi_10, fi_16, fi_21, fi_84, fi_86, fi_87, \
                         fi_91, fi_94, fi_100, fi_105 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_89[k] = f_9 * fi_2[k]
                  - f_8 * fi_7[k]
                  + f_7 * fi_16[k]
                  - f_6 * fi_86[k]
                  + f_5 * fi_91[k]
                  - f_4 * fi_100[k];

        g_90[k] = f_57 * fi_0[k]
                  - f_58 * fi_3[k]
                  + f_58 * fi_10[k]
                  - f_57 * fi_21[k]
                  - f_55 * fi_84[k]
                  + f_56 * fi_87[k]
                  - f_56 * fi_94[k]
                  + f_55 * fi_105[k];
    }
}

}  // namespace simdtrf
