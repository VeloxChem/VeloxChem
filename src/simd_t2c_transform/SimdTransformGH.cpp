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


#include "SimdTransformGH.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_gh(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t gh,
             const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 3.28125 * std::sqrt(10.0);
    const auto f_1 = 6.5625 * std::sqrt(10.0);
    const auto f_2 = 0.65625 * std::sqrt(10.0);
    const auto f_3 = 3.28125 * std::sqrt(2.0);
    const auto f_4 = 2.1875 * std::sqrt(2.0);
    const auto f_5 = 26.25 * std::sqrt(2.0);
    const auto f_6 = 1.09375 * std::sqrt(2.0);
    const auto f_7 = 8.75 * std::sqrt(2.0);
    const auto f_8 = 8.75 * std::sqrt(3.0);
    const auto f_9 = 17.5 * std::sqrt(3.0);
    const auto f_10 = 0.3125 * std::sqrt(21.0);
    const auto f_11 = 0.625 * std::sqrt(21.0);
    const auto f_12 = 3.75 * std::sqrt(21.0);
    const auto f_13 = 2.5 * std::sqrt(21.0);
    const auto f_14 = 0.9375 * std::sqrt(35.0);
    const auto f_15 = 1.875 * std::sqrt(35.0);
    const auto f_16 = 2.5 * std::sqrt(35.0);
    const auto f_17 = 0.5 * std::sqrt(35.0);
    const auto f_18 = 4.375 * std::sqrt(3.0);
    const auto f_19 = 9.84375 * std::sqrt(5.0);
    const auto f_20 = 19.6875 * std::sqrt(5.0);
    const auto f_21 = 1.96875 * std::sqrt(5.0);
    const auto f_22 = 3.28125 * std::sqrt(5.0);
    const auto f_23 = 6.5625 * std::sqrt(5.0);
    const auto f_24 = 0.65625 * std::sqrt(5.0);
    const auto f_25 = 39.375 * std::sqrt(2.0);
    const auto f_26 = 13.125 * std::sqrt(2.0);
    const auto f_27 = 13.125 * std::sqrt(6.0);
    const auto f_28 = 26.25 * std::sqrt(6.0);
    const auto f_29 = 4.375 * std::sqrt(6.0);
    const auto f_30 = 8.75 * std::sqrt(6.0);
    const auto f_31 = 0.46875 * std::sqrt(42.0);
    const auto f_32 = 0.9375 * std::sqrt(42.0);
    const auto f_33 = 5.625 * std::sqrt(42.0);
    const auto f_34 = 3.75 * std::sqrt(42.0);
    const auto f_35 = 0.15625 * std::sqrt(42.0);
    const auto f_36 = 0.3125 * std::sqrt(42.0);
    const auto f_37 = 1.875 * std::sqrt(42.0);
    const auto f_38 = 1.25 * std::sqrt(42.0);
    const auto f_39 = 1.40625 * std::sqrt(70.0);
    const auto f_40 = 2.8125 * std::sqrt(70.0);
    const auto f_41 = 3.75 * std::sqrt(70.0);
    const auto f_42 = 0.75 * std::sqrt(70.0);
    const auto f_43 = 0.46875 * std::sqrt(70.0);
    const auto f_44 = 0.9375 * std::sqrt(70.0);
    const auto f_45 = 1.25 * std::sqrt(70.0);
    const auto f_46 = 0.25 * std::sqrt(70.0);
    const auto f_47 = 6.5625 * std::sqrt(6.0);
    const auto f_48 = 2.1875 * std::sqrt(6.0);
    const auto f_49 = 9.84375 * std::sqrt(2.0);
    const auto f_50 = 59.0625 * std::sqrt(2.0);
    const auto f_51 = 19.6875 * std::sqrt(2.0);
    const auto f_52 = 0.09375 * std::sqrt(70.0);
    const auto f_53 = 5.625 * std::sqrt(70.0);
    const auto f_54 = 0.5625 * std::sqrt(70.0);
    const auto f_55 = 3.75 * std::sqrt(7.0);
    const auto f_56 = 22.5 * std::sqrt(7.0);
    const auto f_57 = 0.46875 * std::sqrt(14.0);
    const auto f_58 = 0.3125 * std::sqrt(14.0);
    const auto f_59 = 3.75 * std::sqrt(14.0);
    const auto f_60 = 0.15625 * std::sqrt(14.0);
    const auto f_61 = 1.25 * std::sqrt(14.0);
    const auto f_62 = 2.8125 * std::sqrt(14.0);
    const auto f_63 = 1.875 * std::sqrt(14.0);
    const auto f_64 = 22.5 * std::sqrt(14.0);
    const auto f_65 = 0.9375 * std::sqrt(14.0);
    const auto f_66 = 7.5 * std::sqrt(14.0);
    const auto f_67 = 1.25 * std::sqrt(21.0);
    const auto f_68 = 7.5 * std::sqrt(21.0);
    const auto f_69 = 15.0 * std::sqrt(21.0);
    const auto f_70 = 0.3125 * std::sqrt(3.0);
    const auto f_71 = 0.625 * std::sqrt(3.0);
    const auto f_72 = 3.75 * std::sqrt(3.0);
    const auto f_73 = 2.5 * std::sqrt(3.0);
    const auto f_74 = 1.875 * std::sqrt(3.0);
    const auto f_75 = 22.5 * std::sqrt(3.0);
    const auto f_76 = 15.0 * std::sqrt(3.0);
    const auto f_77 = 0.9375 * std::sqrt(5.0);
    const auto f_78 = 1.875 * std::sqrt(5.0);
    const auto f_79 = 2.5 * std::sqrt(5.0);
    const auto f_80 = 0.5 * std::sqrt(5.0);
    const auto f_81 = 5.625 * std::sqrt(5.0);
    const auto f_82 = 11.25 * std::sqrt(5.0);
    const auto f_83 = 15.0 * std::sqrt(5.0);
    const auto f_84 = 3.0 * std::sqrt(5.0);
    const auto f_85 = 0.9375 * std::sqrt(7.0);
    const auto f_86 = 5.625 * std::sqrt(7.0);
    const auto f_87 = 33.75 * std::sqrt(7.0);
    const auto f_88 = 1.40625 * std::sqrt(35.0);
    const auto f_89 = 2.8125 * std::sqrt(35.0);
    const auto f_90 = 0.28125 * std::sqrt(35.0);
    const auto f_91 = 3.75 * std::sqrt(35.0);
    const auto f_92 = 0.375 * std::sqrt(35.0);
    const auto f_93 = 5.625 * std::sqrt(14.0);
    const auto f_94 = 1.40625 * std::sqrt(7.0);
    const auto f_95 = 11.25 * std::sqrt(7.0);
    const auto f_96 = 0.46875 * std::sqrt(7.0);
    const auto f_97 = 1.875 * std::sqrt(7.0);
    const auto f_98 = 1.25 * std::sqrt(7.0);
    const auto f_99 = 15.0 * std::sqrt(7.0);
    const auto f_100 = 0.625 * std::sqrt(7.0);
    const auto f_101 = 5.0 * std::sqrt(7.0);
    const auto f_102 = 2.5 * std::sqrt(42.0);
    const auto f_103 = 5.0 * std::sqrt(42.0);
    const auto f_104 = 0.46875 * std::sqrt(6.0);
    const auto f_105 = 0.9375 * std::sqrt(6.0);
    const auto f_106 = 5.625 * std::sqrt(6.0);
    const auto f_107 = 3.75 * std::sqrt(6.0);
    const auto f_108 = 0.625 * std::sqrt(6.0);
    const auto f_109 = 1.25 * std::sqrt(6.0);
    const auto f_110 = 7.5 * std::sqrt(6.0);
    const auto f_111 = 5.0 * std::sqrt(6.0);
    const auto f_112 = 1.40625 * std::sqrt(10.0);
    const auto f_113 = 2.8125 * std::sqrt(10.0);
    const auto f_114 = 3.75 * std::sqrt(10.0);
    const auto f_115 = 0.75 * std::sqrt(10.0);
    const auto f_116 = 1.875 * std::sqrt(10.0);
    const auto f_117 = 5.0 * std::sqrt(10.0);
    const auto f_118 = std::sqrt(10.0);
    const auto f_119 = 1.40625 * std::sqrt(14.0);
    const auto f_120 = 8.4375 * std::sqrt(14.0);
    const auto f_121 = 11.25 * std::sqrt(14.0);
    const auto f_122 = 0.3515625 * std::sqrt(14.0);
    const auto f_123 = 0.703125 * std::sqrt(14.0);
    const auto f_124 = 0.0703125 * std::sqrt(14.0);
    const auto f_125 = 0.140625 * std::sqrt(14.0);
    const auto f_126 = 0.5625 * std::sqrt(14.0);
    const auto f_127 = 0.1875 * std::sqrt(14.0);
    const auto f_128 = 0.5625 * std::sqrt(35.0);
    const auto f_129 = 1.125 * std::sqrt(35.0);
    const auto f_130 = 4.5 * std::sqrt(35.0);
    const auto f_131 = 1.5 * std::sqrt(35.0);
    const auto f_132 = 0.0703125 * std::sqrt(70.0);
    const auto f_133 = 0.046875 * std::sqrt(70.0);
    const auto f_134 = 0.0234375 * std::sqrt(70.0);
    const auto f_135 = 0.1875 * std::sqrt(70.0);
    const auto f_136 = 0.140625 * std::sqrt(70.0);
    const auto f_137 = 1.125 * std::sqrt(70.0);
    const auto f_138 = 0.375 * std::sqrt(70.0);
    const auto f_139 = 4.5 * std::sqrt(70.0);
    const auto f_140 = 1.5 * std::sqrt(70.0);
    const auto f_141 = 0.125 * std::sqrt(70.0);
    const auto f_142 = 0.0625 * std::sqrt(70.0);
    const auto f_143 = 0.5 * std::sqrt(70.0);
    const auto f_144 = 0.1875 * std::sqrt(105.0);
    const auto f_145 = 0.375 * std::sqrt(105.0);
    const auto f_146 = 0.75 * std::sqrt(105.0);
    const auto f_147 = 1.5 * std::sqrt(105.0);
    const auto f_148 = 3.0 * std::sqrt(105.0);
    const auto f_149 = 0.5 * std::sqrt(105.0);
    const auto f_150 = std::sqrt(105.0);
    const auto f_151 = 0.046875 * std::sqrt(15.0);
    const auto f_152 = 0.09375 * std::sqrt(15.0);
    const auto f_153 = 0.5625 * std::sqrt(15.0);
    const auto f_154 = 0.375 * std::sqrt(15.0);
    const auto f_155 = 0.1875 * std::sqrt(15.0);
    const auto f_156 = 1.125 * std::sqrt(15.0);
    const auto f_157 = 0.75 * std::sqrt(15.0);
    const auto f_158 = 4.5 * std::sqrt(15.0);
    const auto f_159 = 3.0 * std::sqrt(15.0);
    const auto f_160 = 0.125 * std::sqrt(15.0);
    const auto f_161 = 0.25 * std::sqrt(15.0);
    const auto f_162 = 1.5 * std::sqrt(15.0);
    const auto f_163 = std::sqrt(15.0);
    const auto f_164 = 0.09375 * std::sqrt(105.0);
    const auto f_165 = 0.25 * std::sqrt(105.0);
    const auto f_166 = 0.140625 * std::sqrt(35.0);
    const auto f_167 = 0.84375 * std::sqrt(35.0);
    const auto f_168 = 1.6875 * std::sqrt(35.0);
    const auto f_169 = 6.75 * std::sqrt(35.0);
    const auto f_170 = 2.25 * std::sqrt(35.0);
    const auto f_171 = 0.234375 * std::sqrt(70.0);
    const auto f_172 = 0.28125 * std::sqrt(70.0);
    const auto f_173 = 0.234375 * std::sqrt(14.0);
    const auto f_174 = 0.078125 * std::sqrt(14.0);
    const auto f_175 = 0.625 * std::sqrt(14.0);
    const auto f_176 = 0.15625 * std::sqrt(3.0);
    const auto f_177 = 1.25 * std::sqrt(3.0);
    const auto f_178 = 0.9375 * std::sqrt(3.0);
    const auto f_179 = 11.25 * std::sqrt(3.0);
    const auto f_180 = 7.5 * std::sqrt(3.0);
    const auto f_181 = 0.46875 * std::sqrt(5.0);
    const auto f_182 = 1.25 * std::sqrt(5.0);
    const auto f_183 = 0.25 * std::sqrt(5.0);
    const auto f_184 = 2.8125 * std::sqrt(5.0);
    const auto f_185 = 7.5 * std::sqrt(5.0);
    const auto f_186 = 1.5 * std::sqrt(5.0);
    const auto f_187 = 1.875 * std::sqrt(21.0);
    const auto f_188 = 2.8125 * std::sqrt(7.0);
    const auto f_189 = 16.875 * std::sqrt(7.0);
    const auto f_190 = 0.8203125 * std::sqrt(10.0);
    const auto f_191 = 1.640625 * std::sqrt(10.0);
    const auto f_192 = 0.1640625 * std::sqrt(10.0);
    const auto f_193 = 4.921875 * std::sqrt(10.0);
    const auto f_194 = 9.84375 * std::sqrt(10.0);
    const auto f_195 = 0.984375 * std::sqrt(10.0);
    const auto f_196 = 0.8203125 * std::sqrt(2.0);
    const auto f_197 = 0.546875 * std::sqrt(2.0);
    const auto f_198 = 6.5625 * std::sqrt(2.0);
    const auto f_199 = 0.2734375 * std::sqrt(2.0);
    const auto f_200 = 4.921875 * std::sqrt(2.0);
    const auto f_201 = 1.640625 * std::sqrt(2.0);
    const auto f_202 = 2.1875 * std::sqrt(3.0);
    const auto f_203 = 13.125 * std::sqrt(3.0);
    const auto f_204 = 26.25 * std::sqrt(3.0);
    const auto f_205 = 0.078125 * std::sqrt(21.0);
    const auto f_206 = 0.15625 * std::sqrt(21.0);
    const auto f_207 = 0.9375 * std::sqrt(21.0);
    const auto f_208 = 0.46875 * std::sqrt(21.0);
    const auto f_209 = 5.625 * std::sqrt(21.0);
    const auto f_210 = 0.234375 * std::sqrt(35.0);
    const auto f_211 = 0.46875 * std::sqrt(35.0);
    const auto f_212 = 0.625 * std::sqrt(35.0);
    const auto f_213 = 0.125 * std::sqrt(35.0);
    const auto f_214 = 0.75 * std::sqrt(35.0);
    const auto f_215 = 1.09375 * std::sqrt(3.0);
    const auto f_216 = 6.5625 * std::sqrt(3.0);

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
    auto *g_91 = values + 91 * nvalues;
    auto *g_92 = values + 92 * nvalues;
    auto *g_93 = values + 93 * nvalues;
    auto *g_94 = values + 94 * nvalues;
    auto *g_95 = values + 95 * nvalues;
    auto *g_96 = values + 96 * nvalues;
    auto *g_97 = values + 97 * nvalues;
    auto *g_98 = values + 98 * nvalues;

    const auto *gh_0 = buffer.data(gh + 0);
    const auto *gh_1 = buffer.data(gh + 1);
    const auto *gh_2 = buffer.data(gh + 2);
    const auto *gh_3 = buffer.data(gh + 3);
    const auto *gh_4 = buffer.data(gh + 4);
    const auto *gh_5 = buffer.data(gh + 5);
    const auto *gh_6 = buffer.data(gh + 6);
    const auto *gh_7 = buffer.data(gh + 7);
    const auto *gh_8 = buffer.data(gh + 8);
    const auto *gh_9 = buffer.data(gh + 9);
    const auto *gh_10 = buffer.data(gh + 10);
    const auto *gh_11 = buffer.data(gh + 11);
    const auto *gh_12 = buffer.data(gh + 12);
    const auto *gh_13 = buffer.data(gh + 13);
    const auto *gh_14 = buffer.data(gh + 14);
    const auto *gh_15 = buffer.data(gh + 15);
    const auto *gh_16 = buffer.data(gh + 16);
    const auto *gh_17 = buffer.data(gh + 17);
    const auto *gh_18 = buffer.data(gh + 18);
    const auto *gh_19 = buffer.data(gh + 19);
    const auto *gh_20 = buffer.data(gh + 20);
    const auto *gh_21 = buffer.data(gh + 21);
    const auto *gh_22 = buffer.data(gh + 22);
    const auto *gh_23 = buffer.data(gh + 23);
    const auto *gh_24 = buffer.data(gh + 24);
    const auto *gh_25 = buffer.data(gh + 25);
    const auto *gh_26 = buffer.data(gh + 26);
    const auto *gh_27 = buffer.data(gh + 27);
    const auto *gh_28 = buffer.data(gh + 28);
    const auto *gh_29 = buffer.data(gh + 29);
    const auto *gh_30 = buffer.data(gh + 30);
    const auto *gh_31 = buffer.data(gh + 31);
    const auto *gh_32 = buffer.data(gh + 32);
    const auto *gh_33 = buffer.data(gh + 33);
    const auto *gh_34 = buffer.data(gh + 34);
    const auto *gh_35 = buffer.data(gh + 35);
    const auto *gh_36 = buffer.data(gh + 36);
    const auto *gh_37 = buffer.data(gh + 37);
    const auto *gh_38 = buffer.data(gh + 38);
    const auto *gh_39 = buffer.data(gh + 39);
    const auto *gh_40 = buffer.data(gh + 40);
    const auto *gh_41 = buffer.data(gh + 41);
    const auto *gh_42 = buffer.data(gh + 42);
    const auto *gh_43 = buffer.data(gh + 43);
    const auto *gh_44 = buffer.data(gh + 44);
    const auto *gh_45 = buffer.data(gh + 45);
    const auto *gh_46 = buffer.data(gh + 46);
    const auto *gh_47 = buffer.data(gh + 47);
    const auto *gh_48 = buffer.data(gh + 48);
    const auto *gh_49 = buffer.data(gh + 49);
    const auto *gh_50 = buffer.data(gh + 50);
    const auto *gh_51 = buffer.data(gh + 51);
    const auto *gh_52 = buffer.data(gh + 52);
    const auto *gh_53 = buffer.data(gh + 53);
    const auto *gh_54 = buffer.data(gh + 54);
    const auto *gh_55 = buffer.data(gh + 55);
    const auto *gh_56 = buffer.data(gh + 56);
    const auto *gh_57 = buffer.data(gh + 57);
    const auto *gh_58 = buffer.data(gh + 58);
    const auto *gh_59 = buffer.data(gh + 59);
    const auto *gh_60 = buffer.data(gh + 60);
    const auto *gh_61 = buffer.data(gh + 61);
    const auto *gh_62 = buffer.data(gh + 62);
    const auto *gh_63 = buffer.data(gh + 63);
    const auto *gh_64 = buffer.data(gh + 64);
    const auto *gh_65 = buffer.data(gh + 65);
    const auto *gh_66 = buffer.data(gh + 66);
    const auto *gh_67 = buffer.data(gh + 67);
    const auto *gh_68 = buffer.data(gh + 68);
    const auto *gh_69 = buffer.data(gh + 69);
    const auto *gh_70 = buffer.data(gh + 70);
    const auto *gh_71 = buffer.data(gh + 71);
    const auto *gh_72 = buffer.data(gh + 72);
    const auto *gh_73 = buffer.data(gh + 73);
    const auto *gh_74 = buffer.data(gh + 74);
    const auto *gh_75 = buffer.data(gh + 75);
    const auto *gh_76 = buffer.data(gh + 76);
    const auto *gh_77 = buffer.data(gh + 77);
    const auto *gh_78 = buffer.data(gh + 78);
    const auto *gh_79 = buffer.data(gh + 79);
    const auto *gh_80 = buffer.data(gh + 80);
    const auto *gh_81 = buffer.data(gh + 81);
    const auto *gh_82 = buffer.data(gh + 82);
    const auto *gh_83 = buffer.data(gh + 83);
    const auto *gh_84 = buffer.data(gh + 84);
    const auto *gh_85 = buffer.data(gh + 85);
    const auto *gh_86 = buffer.data(gh + 86);
    const auto *gh_87 = buffer.data(gh + 87);
    const auto *gh_88 = buffer.data(gh + 88);
    const auto *gh_89 = buffer.data(gh + 89);
    const auto *gh_90 = buffer.data(gh + 90);
    const auto *gh_91 = buffer.data(gh + 91);
    const auto *gh_92 = buffer.data(gh + 92);
    const auto *gh_93 = buffer.data(gh + 93);
    const auto *gh_94 = buffer.data(gh + 94);
    const auto *gh_95 = buffer.data(gh + 95);
    const auto *gh_96 = buffer.data(gh + 96);
    const auto *gh_97 = buffer.data(gh + 97);
    const auto *gh_98 = buffer.data(gh + 98);
    const auto *gh_99 = buffer.data(gh + 99);
    const auto *gh_100 = buffer.data(gh + 100);
    const auto *gh_101 = buffer.data(gh + 101);
    const auto *gh_102 = buffer.data(gh + 102);
    const auto *gh_103 = buffer.data(gh + 103);
    const auto *gh_104 = buffer.data(gh + 104);
    const auto *gh_105 = buffer.data(gh + 105);
    const auto *gh_106 = buffer.data(gh + 106);
    const auto *gh_107 = buffer.data(gh + 107);
    const auto *gh_108 = buffer.data(gh + 108);
    const auto *gh_109 = buffer.data(gh + 109);
    const auto *gh_110 = buffer.data(gh + 110);
    const auto *gh_111 = buffer.data(gh + 111);
    const auto *gh_112 = buffer.data(gh + 112);
    const auto *gh_113 = buffer.data(gh + 113);
    const auto *gh_114 = buffer.data(gh + 114);
    const auto *gh_115 = buffer.data(gh + 115);
    const auto *gh_116 = buffer.data(gh + 116);
    const auto *gh_117 = buffer.data(gh + 117);
    const auto *gh_118 = buffer.data(gh + 118);
    const auto *gh_119 = buffer.data(gh + 119);
    const auto *gh_120 = buffer.data(gh + 120);
    const auto *gh_121 = buffer.data(gh + 121);
    const auto *gh_122 = buffer.data(gh + 122);
    const auto *gh_123 = buffer.data(gh + 123);
    const auto *gh_124 = buffer.data(gh + 124);
    const auto *gh_125 = buffer.data(gh + 125);
    const auto *gh_126 = buffer.data(gh + 126);
    const auto *gh_127 = buffer.data(gh + 127);
    const auto *gh_128 = buffer.data(gh + 128);
    const auto *gh_129 = buffer.data(gh + 129);
    const auto *gh_130 = buffer.data(gh + 130);
    const auto *gh_131 = buffer.data(gh + 131);
    const auto *gh_132 = buffer.data(gh + 132);
    const auto *gh_133 = buffer.data(gh + 133);
    const auto *gh_134 = buffer.data(gh + 134);
    const auto *gh_135 = buffer.data(gh + 135);
    const auto *gh_136 = buffer.data(gh + 136);
    const auto *gh_137 = buffer.data(gh + 137);
    const auto *gh_138 = buffer.data(gh + 138);
    const auto *gh_139 = buffer.data(gh + 139);
    const auto *gh_140 = buffer.data(gh + 140);
    const auto *gh_141 = buffer.data(gh + 141);
    const auto *gh_142 = buffer.data(gh + 142);
    const auto *gh_143 = buffer.data(gh + 143);
    const auto *gh_144 = buffer.data(gh + 144);
    const auto *gh_145 = buffer.data(gh + 145);
    const auto *gh_146 = buffer.data(gh + 146);
    const auto *gh_147 = buffer.data(gh + 147);
    const auto *gh_148 = buffer.data(gh + 148);
    const auto *gh_149 = buffer.data(gh + 149);
    const auto *gh_150 = buffer.data(gh + 150);
    const auto *gh_151 = buffer.data(gh + 151);
    const auto *gh_152 = buffer.data(gh + 152);
    const auto *gh_153 = buffer.data(gh + 153);
    const auto *gh_154 = buffer.data(gh + 154);
    const auto *gh_155 = buffer.data(gh + 155);
    const auto *gh_156 = buffer.data(gh + 156);
    const auto *gh_157 = buffer.data(gh + 157);
    const auto *gh_158 = buffer.data(gh + 158);
    const auto *gh_159 = buffer.data(gh + 159);
    const auto *gh_160 = buffer.data(gh + 160);
    const auto *gh_161 = buffer.data(gh + 161);
    const auto *gh_162 = buffer.data(gh + 162);
    const auto *gh_163 = buffer.data(gh + 163);
    const auto *gh_164 = buffer.data(gh + 164);
    const auto *gh_165 = buffer.data(gh + 165);
    const auto *gh_166 = buffer.data(gh + 166);
    const auto *gh_167 = buffer.data(gh + 167);
    const auto *gh_168 = buffer.data(gh + 168);
    const auto *gh_169 = buffer.data(gh + 169);
    const auto *gh_170 = buffer.data(gh + 170);
    const auto *gh_171 = buffer.data(gh + 171);
    const auto *gh_172 = buffer.data(gh + 172);
    const auto *gh_173 = buffer.data(gh + 173);
    const auto *gh_174 = buffer.data(gh + 174);
    const auto *gh_175 = buffer.data(gh + 175);
    const auto *gh_176 = buffer.data(gh + 176);
    const auto *gh_177 = buffer.data(gh + 177);
    const auto *gh_178 = buffer.data(gh + 178);
    const auto *gh_179 = buffer.data(gh + 179);
    const auto *gh_180 = buffer.data(gh + 180);
    const auto *gh_181 = buffer.data(gh + 181);
    const auto *gh_182 = buffer.data(gh + 182);
    const auto *gh_183 = buffer.data(gh + 183);
    const auto *gh_184 = buffer.data(gh + 184);
    const auto *gh_185 = buffer.data(gh + 185);
    const auto *gh_186 = buffer.data(gh + 186);
    const auto *gh_187 = buffer.data(gh + 187);
    const auto *gh_188 = buffer.data(gh + 188);
    const auto *gh_189 = buffer.data(gh + 189);
    const auto *gh_190 = buffer.data(gh + 190);
    const auto *gh_191 = buffer.data(gh + 191);
    const auto *gh_192 = buffer.data(gh + 192);
    const auto *gh_193 = buffer.data(gh + 193);
    const auto *gh_194 = buffer.data(gh + 194);
    const auto *gh_195 = buffer.data(gh + 195);
    const auto *gh_196 = buffer.data(gh + 196);
    const auto *gh_197 = buffer.data(gh + 197);
    const auto *gh_198 = buffer.data(gh + 198);
    const auto *gh_199 = buffer.data(gh + 199);
    const auto *gh_200 = buffer.data(gh + 200);
    const auto *gh_201 = buffer.data(gh + 201);
    const auto *gh_202 = buffer.data(gh + 202);
    const auto *gh_203 = buffer.data(gh + 203);
    const auto *gh_204 = buffer.data(gh + 204);
    const auto *gh_205 = buffer.data(gh + 205);
    const auto *gh_206 = buffer.data(gh + 206);
    const auto *gh_207 = buffer.data(gh + 207);
    const auto *gh_208 = buffer.data(gh + 208);
    const auto *gh_209 = buffer.data(gh + 209);
    const auto *gh_210 = buffer.data(gh + 210);
    const auto *gh_211 = buffer.data(gh + 211);
    const auto *gh_212 = buffer.data(gh + 212);
    const auto *gh_213 = buffer.data(gh + 213);
    const auto *gh_214 = buffer.data(gh + 214);
    const auto *gh_215 = buffer.data(gh + 215);
    const auto *gh_216 = buffer.data(gh + 216);
    const auto *gh_217 = buffer.data(gh + 217);
    const auto *gh_218 = buffer.data(gh + 218);
    const auto *gh_219 = buffer.data(gh + 219);
    const auto *gh_220 = buffer.data(gh + 220);
    const auto *gh_221 = buffer.data(gh + 221);
    const auto *gh_222 = buffer.data(gh + 222);
    const auto *gh_223 = buffer.data(gh + 223);
    const auto *gh_224 = buffer.data(gh + 224);
    const auto *gh_225 = buffer.data(gh + 225);
    const auto *gh_226 = buffer.data(gh + 226);
    const auto *gh_227 = buffer.data(gh + 227);
    const auto *gh_228 = buffer.data(gh + 228);
    const auto *gh_229 = buffer.data(gh + 229);
    const auto *gh_230 = buffer.data(gh + 230);
    const auto *gh_231 = buffer.data(gh + 231);
    const auto *gh_232 = buffer.data(gh + 232);
    const auto *gh_233 = buffer.data(gh + 233);
    const auto *gh_234 = buffer.data(gh + 234);
    const auto *gh_235 = buffer.data(gh + 235);
    const auto *gh_236 = buffer.data(gh + 236);
    const auto *gh_237 = buffer.data(gh + 237);
    const auto *gh_238 = buffer.data(gh + 238);
    const auto *gh_239 = buffer.data(gh + 239);
    const auto *gh_240 = buffer.data(gh + 240);
    const auto *gh_241 = buffer.data(gh + 241);
    const auto *gh_242 = buffer.data(gh + 242);
    const auto *gh_243 = buffer.data(gh + 243);
    const auto *gh_244 = buffer.data(gh + 244);
    const auto *gh_245 = buffer.data(gh + 245);
    const auto *gh_246 = buffer.data(gh + 246);
    const auto *gh_247 = buffer.data(gh + 247);
    const auto *gh_248 = buffer.data(gh + 248);
    const auto *gh_249 = buffer.data(gh + 249);
    const auto *gh_250 = buffer.data(gh + 250);
    const auto *gh_251 = buffer.data(gh + 251);
    const auto *gh_252 = buffer.data(gh + 252);
    const auto *gh_253 = buffer.data(gh + 253);
    const auto *gh_254 = buffer.data(gh + 254);
    const auto *gh_255 = buffer.data(gh + 255);
    const auto *gh_256 = buffer.data(gh + 256);
    const auto *gh_257 = buffer.data(gh + 257);
    const auto *gh_258 = buffer.data(gh + 258);
    const auto *gh_259 = buffer.data(gh + 259);
    const auto *gh_260 = buffer.data(gh + 260);
    const auto *gh_261 = buffer.data(gh + 261);
    const auto *gh_262 = buffer.data(gh + 262);
    const auto *gh_263 = buffer.data(gh + 263);
    const auto *gh_264 = buffer.data(gh + 264);
    const auto *gh_265 = buffer.data(gh + 265);
    const auto *gh_266 = buffer.data(gh + 266);
    const auto *gh_267 = buffer.data(gh + 267);
    const auto *gh_268 = buffer.data(gh + 268);
    const auto *gh_269 = buffer.data(gh + 269);
    const auto *gh_270 = buffer.data(gh + 270);
    const auto *gh_271 = buffer.data(gh + 271);
    const auto *gh_272 = buffer.data(gh + 272);
    const auto *gh_273 = buffer.data(gh + 273);
    const auto *gh_274 = buffer.data(gh + 274);
    const auto *gh_275 = buffer.data(gh + 275);
    const auto *gh_276 = buffer.data(gh + 276);
    const auto *gh_277 = buffer.data(gh + 277);
    const auto *gh_278 = buffer.data(gh + 278);
    const auto *gh_279 = buffer.data(gh + 279);
    const auto *gh_280 = buffer.data(gh + 280);
    const auto *gh_281 = buffer.data(gh + 281);
    const auto *gh_282 = buffer.data(gh + 282);
    const auto *gh_283 = buffer.data(gh + 283);
    const auto *gh_284 = buffer.data(gh + 284);
    const auto *gh_285 = buffer.data(gh + 285);
    const auto *gh_286 = buffer.data(gh + 286);
    const auto *gh_287 = buffer.data(gh + 287);
    const auto *gh_288 = buffer.data(gh + 288);
    const auto *gh_289 = buffer.data(gh + 289);
    const auto *gh_290 = buffer.data(gh + 290);
    const auto *gh_291 = buffer.data(gh + 291);
    const auto *gh_292 = buffer.data(gh + 292);
    const auto *gh_293 = buffer.data(gh + 293);
    const auto *gh_294 = buffer.data(gh + 294);
    const auto *gh_295 = buffer.data(gh + 295);
    const auto *gh_296 = buffer.data(gh + 296);
    const auto *gh_297 = buffer.data(gh + 297);
    const auto *gh_298 = buffer.data(gh + 298);
    const auto *gh_299 = buffer.data(gh + 299);
    const auto *gh_300 = buffer.data(gh + 300);
    const auto *gh_301 = buffer.data(gh + 301);
    const auto *gh_302 = buffer.data(gh + 302);
    const auto *gh_303 = buffer.data(gh + 303);
    const auto *gh_304 = buffer.data(gh + 304);
    const auto *gh_305 = buffer.data(gh + 305);
    const auto *gh_306 = buffer.data(gh + 306);
    const auto *gh_307 = buffer.data(gh + 307);
    const auto *gh_308 = buffer.data(gh + 308);
    const auto *gh_309 = buffer.data(gh + 309);
    const auto *gh_310 = buffer.data(gh + 310);
    const auto *gh_311 = buffer.data(gh + 311);
    const auto *gh_312 = buffer.data(gh + 312);
    const auto *gh_313 = buffer.data(gh + 313);
    const auto *gh_314 = buffer.data(gh + 314);

#pragma omp simd aligned(gh_22, gh_25, gh_27, gh_32, gh_36, gh_127, gh_130, gh_132, gh_137, \
                         gh_141 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * gh_22[k]
                 - f_1 * gh_27[k]
                 + f_2 * gh_36[k]
                 - f_0 * gh_127[k]
                 + f_1 * gh_132[k]
                 - f_2 * gh_141[k];

        g_1[k] = 26.25 * gh_25[k]
                 - 26.25 * gh_32[k]
                 - 26.25 * gh_130[k]
                 + 26.25 * gh_137[k];
    }

#pragma omp simd aligned(gh_22, gh_27, gh_29, gh_36, gh_38, gh_127, gh_132, gh_134, gh_141, \
                         gh_143 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_2[k] = -f_3 * gh_22[k]
                 - f_4 * gh_27[k]
                 + f_5 * gh_29[k]
                 + f_6 * gh_36[k]
                 - f_7 * gh_38[k]
                 + f_3 * gh_127[k]
                 + f_4 * gh_132[k]
                 - f_5 * gh_134[k]
                 - f_6 * gh_141[k]
                 + f_7 * gh_143[k];
    }

#pragma omp simd aligned(gh_25, gh_32, gh_34, gh_130, gh_137, gh_139 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = -f_8 * gh_25[k]
                 - f_8 * gh_32[k]
                 + f_9 * gh_34[k]
                 + f_8 * gh_130[k]
                 + f_8 * gh_137[k]
                 - f_9 * gh_139[k];
    }

#pragma omp simd aligned(gh_22, gh_27, gh_29, gh_36, gh_38, gh_40, gh_127, gh_132, gh_134, \
                         gh_141, gh_143, gh_145 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = f_10 * gh_22[k]
                 + f_11 * gh_27[k]
                 - f_12 * gh_29[k]
                 + f_10 * gh_36[k]
                 - f_12 * gh_38[k]
                 + f_13 * gh_40[k]
                 - f_10 * gh_127[k]
                 - f_11 * gh_132[k]
                 + f_12 * gh_134[k]
                 - f_10 * gh_141[k]
                 + f_12 * gh_143[k]
                 - f_13 * gh_145[k];
    }

#pragma omp simd aligned(gh_23, gh_28, gh_30, gh_37, gh_39, gh_41, gh_128, gh_133, gh_135, \
                         gh_142, gh_144, gh_146 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = f_14 * gh_23[k]
                 + f_15 * gh_28[k]
                 - f_16 * gh_30[k]
                 + f_14 * gh_37[k]
                 - f_16 * gh_39[k]
                 + f_17 * gh_41[k]
                 - f_14 * gh_128[k]
                 - f_15 * gh_133[k]
                 + f_16 * gh_135[k]
                 - f_14 * gh_142[k]
                 + f_16 * gh_144[k]
                 - f_17 * gh_146[k];
    }

#pragma omp simd aligned(gh_21, gh_24, gh_26, gh_31, gh_33, gh_35, gh_126, gh_129, gh_131, \
                         gh_136, gh_138, gh_140 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = f_10 * gh_21[k]
                 + f_11 * gh_24[k]
                 - f_12 * gh_26[k]
                 + f_10 * gh_31[k]
                 - f_12 * gh_33[k]
                 + f_13 * gh_35[k]
                 - f_10 * gh_126[k]
                 - f_11 * gh_129[k]
                 + f_12 * gh_131[k]
                 - f_10 * gh_136[k]
                 + f_12 * gh_138[k]
                 - f_13 * gh_140[k];
    }

#pragma omp simd aligned(gh_23, gh_30, gh_37, gh_39, gh_128, gh_135, gh_142, \
                         gh_144 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = -f_18 * gh_23[k]
                 + f_8 * gh_30[k]
                 + f_18 * gh_37[k]
                 - f_8 * gh_39[k]
                 + f_18 * gh_128[k]
                 - f_8 * gh_135[k]
                 - f_18 * gh_142[k]
                 + f_8 * gh_144[k];
    }

#pragma omp simd aligned(gh_21, gh_24, gh_26, gh_31, gh_33, gh_126, gh_129, gh_131, gh_136, \
                         gh_138 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_8[k] = -f_6 * gh_21[k]
                 + f_4 * gh_24[k]
                 + f_7 * gh_26[k]
                 + f_3 * gh_31[k]
                 - f_5 * gh_33[k]
                 + f_6 * gh_126[k]
                 - f_4 * gh_129[k]
                 - f_7 * gh_131[k]
                 - f_3 * gh_136[k]
                 + f_5 * gh_138[k];
    }

#pragma omp simd aligned(gh_21, gh_23, gh_24, gh_28, gh_31, gh_37, gh_126, gh_128, gh_129, \
                         gh_133, gh_136, gh_142 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_9[k] = 6.5625 * gh_23[k]
                 - 39.375 * gh_28[k]
                 + 6.5625 * gh_37[k]
                 - 6.5625 * gh_128[k]
                 + 39.375 * gh_133[k]
                 - 6.5625 * gh_142[k];

        g_10[k] = f_2 * gh_21[k]
                  - f_1 * gh_24[k]
                  + f_0 * gh_31[k]
                  - f_2 * gh_126[k]
                  + f_1 * gh_129[k]
                  - f_0 * gh_136[k];
    }

#pragma omp simd aligned(gh_85, gh_88, gh_90, gh_95, gh_99, gh_232, gh_235, gh_237, gh_242, \
                         gh_246 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_11[k] = f_19 * gh_85[k]
                  - f_20 * gh_90[k]
                  + f_21 * gh_99[k]
                  - f_22 * gh_232[k]
                  + f_23 * gh_237[k]
                  - f_24 * gh_246[k];

        g_12[k] = f_25 * gh_88[k]
                  - f_25 * gh_95[k]
                  - f_26 * gh_235[k]
                  + f_26 * gh_242[k];
    }

#pragma omp simd aligned(gh_85, gh_90, gh_92, gh_99, gh_101, gh_232, gh_237, gh_239, gh_246, \
                         gh_248 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = -9.84375 * gh_85[k]
                  - 6.5625 * gh_90[k]
                  + 78.75 * gh_92[k]
                  + 3.28125 * gh_99[k]
                  - 26.25 * gh_101[k]
                  + 3.28125 * gh_232[k]
                  + 2.1875 * gh_237[k]
                  - 26.25 * gh_239[k]
                  - 1.09375 * gh_246[k]
                  + 8.75 * gh_248[k];
    }

#pragma omp simd aligned(gh_88, gh_95, gh_97, gh_235, gh_242, gh_244 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = -f_27 * gh_88[k]
                  - f_27 * gh_95[k]
                  + f_28 * gh_97[k]
                  + f_29 * gh_235[k]
                  + f_29 * gh_242[k]
                  - f_30 * gh_244[k];
    }

#pragma omp simd aligned(gh_85, gh_90, gh_92, gh_99, gh_101, gh_103, gh_232, gh_237, gh_239, \
                         gh_246, gh_248, gh_250 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_15[k] = f_31 * gh_85[k]
                  + f_32 * gh_90[k]
                  - f_33 * gh_92[k]
                  + f_31 * gh_99[k]
                  - f_33 * gh_101[k]
                  + f_34 * gh_103[k]
                  - f_35 * gh_232[k]
                  - f_36 * gh_237[k]
                  + f_37 * gh_239[k]
                  - f_35 * gh_246[k]
                  + f_37 * gh_248[k]
                  - f_38 * gh_250[k];
    }

#pragma omp simd aligned(gh_86, gh_91, gh_93, gh_100, gh_102, gh_104, gh_233, gh_238, gh_240, \
                         gh_247, gh_249, gh_251 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_16[k] = f_39 * gh_86[k]
                  + f_40 * gh_91[k]
                  - f_41 * gh_93[k]
                  + f_39 * gh_100[k]
                  - f_41 * gh_102[k]
                  + f_42 * gh_104[k]
                  - f_43 * gh_233[k]
                  - f_44 * gh_238[k]
                  + f_45 * gh_240[k]
                  - f_43 * gh_247[k]
                  + f_45 * gh_249[k]
                  - f_46 * gh_251[k];
    }

#pragma omp simd aligned(gh_84, gh_87, gh_89, gh_94, gh_96, gh_98, gh_231, gh_234, gh_236, \
                         gh_241, gh_243, gh_245 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = f_31 * gh_84[k]
                  + f_32 * gh_87[k]
                  - f_33 * gh_89[k]
                  + f_31 * gh_94[k]
                  - f_33 * gh_96[k]
                  + f_34 * gh_98[k]
                  - f_35 * gh_231[k]
                  - f_36 * gh_234[k]
                  + f_37 * gh_236[k]
                  - f_35 * gh_241[k]
                  + f_37 * gh_243[k]
                  - f_38 * gh_245[k];
    }

#pragma omp simd aligned(gh_86, gh_93, gh_100, gh_102, gh_233, gh_240, gh_247, \
                         gh_249 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = -f_47 * gh_86[k]
                  + f_27 * gh_93[k]
                  + f_47 * gh_100[k]
                  - f_27 * gh_102[k]
                  + f_48 * gh_233[k]
                  - f_29 * gh_240[k]
                  - f_48 * gh_247[k]
                  + f_29 * gh_249[k];
    }

#pragma omp simd aligned(gh_84, gh_87, gh_89, gh_94, gh_96, gh_231, gh_234, gh_236, gh_241, \
                         gh_243 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = -3.28125 * gh_84[k]
                  + 6.5625 * gh_87[k]
                  + 26.25 * gh_89[k]
                  + 9.84375 * gh_94[k]
                  - 78.75 * gh_96[k]
                  + 1.09375 * gh_231[k]
                  - 2.1875 * gh_234[k]
                  - 8.75 * gh_236[k]
                  - 3.28125 * gh_241[k]
                  + 26.25 * gh_243[k];
    }

#pragma omp simd aligned(gh_84, gh_86, gh_87, gh_91, gh_94, gh_100, gh_231, gh_233, gh_234, \
                         gh_238, gh_241, gh_247 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = f_49 * gh_86[k]
                  - f_50 * gh_91[k]
                  + f_49 * gh_100[k]
                  - f_3 * gh_233[k]
                  + f_51 * gh_238[k]
                  - f_3 * gh_247[k];

        g_21[k] = f_21 * gh_84[k]
                  - f_20 * gh_87[k]
                  + f_19 * gh_94[k]
                  - f_24 * gh_231[k]
                  + f_23 * gh_234[k]
                  - f_22 * gh_241[k];
    }

#pragma omp simd aligned(gh_22, gh_27, gh_36, gh_127, gh_132, gh_141, gh_169, gh_174, \
                         gh_183 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = -f_43 * gh_22[k]
                  + f_44 * gh_27[k]
                  - f_52 * gh_36[k]
                  - f_43 * gh_127[k]
                  + f_44 * gh_132[k]
                  - f_52 * gh_141[k]
                  + f_40 * gh_169[k]
                  - f_53 * gh_174[k]
                  + f_54 * gh_183[k];
    }

#pragma omp simd aligned(gh_25, gh_32, gh_130, gh_137, gh_172, gh_179 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = -f_55 * gh_25[k]
                  + f_55 * gh_32[k]
                  - f_55 * gh_130[k]
                  + f_55 * gh_137[k]
                  + f_56 * gh_172[k]
                  - f_56 * gh_179[k];
    }

#pragma omp simd aligned(gh_22, gh_27, gh_29, gh_36, gh_38, gh_127, gh_132, gh_134, gh_141, \
                         gh_143, gh_169, gh_174, gh_176, gh_183, \
                         gh_185 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = f_57 * gh_22[k]
                  + f_58 * gh_27[k]
                  - f_59 * gh_29[k]
                  - f_60 * gh_36[k]
                  + f_61 * gh_38[k]
                  + f_57 * gh_127[k]
                  + f_58 * gh_132[k]
                  - f_59 * gh_134[k]
                  - f_60 * gh_141[k]
                  + f_61 * gh_143[k]
                  - f_62 * gh_169[k]
                  - f_63 * gh_174[k]
                  + f_64 * gh_176[k]
                  + f_65 * gh_183[k]
                  - f_66 * gh_185[k];
    }

#pragma omp simd aligned(gh_25, gh_32, gh_34, gh_130, gh_137, gh_139, gh_172, gh_179, \
                         gh_181 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = f_67 * gh_25[k]
                  + f_67 * gh_32[k]
                  - f_13 * gh_34[k]
                  + f_67 * gh_130[k]
                  + f_67 * gh_137[k]
                  - f_13 * gh_139[k]
                  - f_68 * gh_172[k]
                  - f_68 * gh_179[k]
                  + f_69 * gh_181[k];
    }

#pragma omp simd aligned(gh_22, gh_27, gh_29, gh_36, gh_38, gh_40, gh_127, gh_132, gh_134, \
                         gh_141, gh_143, gh_145, gh_169, gh_174, gh_176, gh_183, gh_185, \
                         gh_187 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = -f_70 * gh_22[k]
                  - f_71 * gh_27[k]
                  + f_72 * gh_29[k]
                  - f_70 * gh_36[k]
                  + f_72 * gh_38[k]
                  - f_73 * gh_40[k]
                  - f_70 * gh_127[k]
                  - f_71 * gh_132[k]
                  + f_72 * gh_134[k]
                  - f_70 * gh_141[k]
                  + f_72 * gh_143[k]
                  - f_73 * gh_145[k]
                  + f_74 * gh_169[k]
                  + f_72 * gh_174[k]
                  - f_75 * gh_176[k]
                  + f_74 * gh_183[k]
                  - f_75 * gh_185[k]
                  + f_76 * gh_187[k];
    }

#pragma omp simd aligned(gh_23, gh_28, gh_30, gh_37, gh_39, gh_41, gh_128, gh_133, gh_135, \
                         gh_142, gh_144, gh_146, gh_170, gh_175, gh_177, gh_184, gh_186, \
                         gh_188 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = -f_77 * gh_23[k]
                  - f_78 * gh_28[k]
                  + f_79 * gh_30[k]
                  - f_77 * gh_37[k]
                  + f_79 * gh_39[k]
                  - f_80 * gh_41[k]
                  - f_77 * gh_128[k]
                  - f_78 * gh_133[k]
                  + f_79 * gh_135[k]
                  - f_77 * gh_142[k]
                  + f_79 * gh_144[k]
                  - f_80 * gh_146[k]
                  + f_81 * gh_170[k]
                  + f_82 * gh_175[k]
                  - f_83 * gh_177[k]
                  + f_81 * gh_184[k]
                  - f_83 * gh_186[k]
                  + f_84 * gh_188[k];
    }

#pragma omp simd aligned(gh_21, gh_24, gh_26, gh_31, gh_33, gh_35, gh_126, gh_129, gh_131, \
                         gh_136, gh_138, gh_140, gh_168, gh_171, gh_173, gh_178, gh_180, \
                         gh_182 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = -f_70 * gh_21[k]
                  - f_71 * gh_24[k]
                  + f_72 * gh_26[k]
                  - f_70 * gh_31[k]
                  + f_72 * gh_33[k]
                  - f_73 * gh_35[k]
                  - f_70 * gh_126[k]
                  - f_71 * gh_129[k]
                  + f_72 * gh_131[k]
                  - f_70 * gh_136[k]
                  + f_72 * gh_138[k]
                  - f_73 * gh_140[k]
                  + f_74 * gh_168[k]
                  + f_72 * gh_171[k]
                  - f_75 * gh_173[k]
                  + f_74 * gh_178[k]
                  - f_75 * gh_180[k]
                  + f_76 * gh_182[k];
    }

#pragma omp simd aligned(gh_23, gh_30, gh_37, gh_39, gh_128, gh_135, gh_142, gh_144, gh_170, \
                         gh_177, gh_184, gh_186 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = f_11 * gh_23[k]
                  - f_67 * gh_30[k]
                  - f_11 * gh_37[k]
                  + f_67 * gh_39[k]
                  + f_11 * gh_128[k]
                  - f_67 * gh_135[k]
                  - f_11 * gh_142[k]
                  + f_67 * gh_144[k]
                  - f_12 * gh_170[k]
                  + f_68 * gh_177[k]
                  + f_12 * gh_184[k]
                  - f_68 * gh_186[k];
    }

#pragma omp simd aligned(gh_21, gh_24, gh_26, gh_31, gh_33, gh_126, gh_129, gh_131, gh_136, \
                         gh_138, gh_168, gh_171, gh_173, gh_178, \
                         gh_180 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = f_60 * gh_21[k]
                  - f_58 * gh_24[k]
                  - f_61 * gh_26[k]
                  - f_57 * gh_31[k]
                  + f_59 * gh_33[k]
                  + f_60 * gh_126[k]
                  - f_58 * gh_129[k]
                  - f_61 * gh_131[k]
                  - f_57 * gh_136[k]
                  + f_59 * gh_138[k]
                  - f_65 * gh_168[k]
                  + f_63 * gh_171[k]
                  + f_66 * gh_173[k]
                  + f_62 * gh_178[k]
                  - f_64 * gh_180[k];
    }

#pragma omp simd aligned(gh_23, gh_28, gh_37, gh_128, gh_133, gh_142, gh_170, gh_175, \
                         gh_184 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = -f_85 * gh_23[k]
                  + f_86 * gh_28[k]
                  - f_85 * gh_37[k]
                  - f_85 * gh_128[k]
                  + f_86 * gh_133[k]
                  - f_85 * gh_142[k]
                  + f_86 * gh_170[k]
                  - f_87 * gh_175[k]
                  + f_86 * gh_184[k];
    }

#pragma omp simd aligned(gh_21, gh_24, gh_31, gh_126, gh_129, gh_136, gh_168, gh_171, \
                         gh_178 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = -f_52 * gh_21[k]
                  + f_44 * gh_24[k]
                  - f_43 * gh_31[k]
                  - f_52 * gh_126[k]
                  + f_44 * gh_129[k]
                  - f_43 * gh_136[k]
                  + f_54 * gh_168[k]
                  - f_53 * gh_171[k]
                  + f_40 * gh_178[k];
    }

#pragma omp simd aligned(gh_85, gh_90, gh_99, gh_232, gh_237, gh_246, gh_274, gh_279, \
                         gh_288 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = -f_88 * gh_85[k]
                  + f_89 * gh_90[k]
                  - f_90 * gh_99[k]
                  - f_88 * gh_232[k]
                  + f_89 * gh_237[k]
                  - f_90 * gh_246[k]
                  + f_15 * gh_274[k]
                  - f_91 * gh_279[k]
                  + f_92 * gh_288[k];
    }

#pragma omp simd aligned(gh_88, gh_95, gh_235, gh_242, gh_277, gh_284 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = -f_93 * gh_88[k]
                  + f_93 * gh_95[k]
                  - f_93 * gh_235[k]
                  + f_93 * gh_242[k]
                  + f_66 * gh_277[k]
                  - f_66 * gh_284[k];
    }

#pragma omp simd aligned(gh_85, gh_90, gh_92, gh_99, gh_101, gh_232, gh_237, gh_239, gh_246, \
                         gh_248, gh_274, gh_279, gh_281, gh_288, \
                         gh_290 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = f_94 * gh_85[k]
                  + f_85 * gh_90[k]
                  - f_95 * gh_92[k]
                  - f_96 * gh_99[k]
                  + f_55 * gh_101[k]
                  + f_94 * gh_232[k]
                  + f_85 * gh_237[k]
                  - f_95 * gh_239[k]
                  - f_96 * gh_246[k]
                  + f_55 * gh_248[k]
                  - f_97 * gh_274[k]
                  - f_98 * gh_279[k]
                  + f_99 * gh_281[k]
                  + f_100 * gh_288[k]
                  - f_101 * gh_290[k];
    }

#pragma omp simd aligned(gh_88, gh_95, gh_97, gh_235, gh_242, gh_244, gh_277, gh_284, \
                         gh_286 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_36[k] = f_37 * gh_88[k]
                  + f_37 * gh_95[k]
                  - f_34 * gh_97[k]
                  + f_37 * gh_235[k]
                  + f_37 * gh_242[k]
                  - f_34 * gh_244[k]
                  - f_102 * gh_277[k]
                  - f_102 * gh_284[k]
                  + f_103 * gh_286[k];
    }

#pragma omp simd aligned(gh_85, gh_90, gh_92, gh_99, gh_101, gh_103, gh_232, gh_237, gh_239, \
                         gh_246, gh_248, gh_250, gh_274, gh_279, gh_281, gh_288, gh_290, \
                         gh_292 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_37[k] = -f_104 * gh_85[k]
                  - f_105 * gh_90[k]
                  + f_106 * gh_92[k]
                  - f_104 * gh_99[k]
                  + f_106 * gh_101[k]
                  - f_107 * gh_103[k]
                  - f_104 * gh_232[k]
                  - f_105 * gh_237[k]
                  + f_106 * gh_239[k]
                  - f_104 * gh_246[k]
                  + f_106 * gh_248[k]
                  - f_107 * gh_250[k]
                  + f_108 * gh_274[k]
                  + f_109 * gh_279[k]
                  - f_110 * gh_281[k]
                  + f_108 * gh_288[k]
                  - f_110 * gh_290[k]
                  + f_111 * gh_292[k];
    }

#pragma omp simd aligned(gh_86, gh_91, gh_93, gh_100, gh_102, gh_104, gh_233, gh_238, gh_240, \
                         gh_247, gh_249, gh_251, gh_275, gh_280, gh_282, gh_289, gh_291, \
                         gh_293 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = -f_112 * gh_86[k]
                  - f_113 * gh_91[k]
                  + f_114 * gh_93[k]
                  - f_112 * gh_100[k]
                  + f_114 * gh_102[k]
                  - f_115 * gh_104[k]
                  - f_112 * gh_233[k]
                  - f_113 * gh_238[k]
                  + f_114 * gh_240[k]
                  - f_112 * gh_247[k]
                  + f_114 * gh_249[k]
                  - f_115 * gh_251[k]
                  + f_116 * gh_275[k]
                  + f_114 * gh_280[k]
                  - f_117 * gh_282[k]
                  + f_116 * gh_289[k]
                  - f_117 * gh_291[k]
                  + f_118 * gh_293[k];
    }

#pragma omp simd aligned(gh_84, gh_87, gh_89, gh_94, gh_96, gh_98, gh_231, gh_234, gh_236, \
                         gh_241, gh_243, gh_245, gh_273, gh_276, gh_278, gh_283, gh_285, \
                         gh_287 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_39[k] = -f_104 * gh_84[k]
                  - f_105 * gh_87[k]
                  + f_106 * gh_89[k]
                  - f_104 * gh_94[k]
                  + f_106 * gh_96[k]
                  - f_107 * gh_98[k]
                  - f_104 * gh_231[k]
                  - f_105 * gh_234[k]
                  + f_106 * gh_236[k]
                  - f_104 * gh_241[k]
                  + f_106 * gh_243[k]
                  - f_107 * gh_245[k]
                  + f_108 * gh_273[k]
                  + f_109 * gh_276[k]
                  - f_110 * gh_278[k]
                  + f_108 * gh_283[k]
                  - f_110 * gh_285[k]
                  + f_111 * gh_287[k];
    }

#pragma omp simd aligned(gh_86, gh_93, gh_100, gh_102, gh_233, gh_240, gh_247, gh_249, gh_275, \
                         gh_282, gh_289, gh_291 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_40[k] = f_32 * gh_86[k]
                  - f_37 * gh_93[k]
                  - f_32 * gh_100[k]
                  + f_37 * gh_102[k]
                  + f_32 * gh_233[k]
                  - f_37 * gh_240[k]
                  - f_32 * gh_247[k]
                  + f_37 * gh_249[k]
                  - f_38 * gh_275[k]
                  + f_102 * gh_282[k]
                  + f_38 * gh_289[k]
                  - f_102 * gh_291[k];
    }

#pragma omp simd aligned(gh_84, gh_87, gh_89, gh_94, gh_96, gh_231, gh_234, gh_236, gh_241, \
                         gh_243, gh_273, gh_276, gh_278, gh_283, \
                         gh_285 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_41[k] = f_96 * gh_84[k]
                  - f_85 * gh_87[k]
                  - f_55 * gh_89[k]
                  - f_94 * gh_94[k]
                  + f_95 * gh_96[k]
                  + f_96 * gh_231[k]
                  - f_85 * gh_234[k]
                  - f_55 * gh_236[k]
                  - f_94 * gh_241[k]
                  + f_95 * gh_243[k]
                  - f_100 * gh_273[k]
                  + f_98 * gh_276[k]
                  + f_101 * gh_278[k]
                  + f_97 * gh_283[k]
                  - f_99 * gh_285[k];
    }

#pragma omp simd aligned(gh_86, gh_91, gh_100, gh_233, gh_238, gh_247, gh_275, gh_280, \
                         gh_289 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = -f_119 * gh_86[k]
                  + f_120 * gh_91[k]
                  - f_119 * gh_100[k]
                  - f_119 * gh_233[k]
                  + f_120 * gh_238[k]
                  - f_119 * gh_247[k]
                  + f_63 * gh_275[k]
                  - f_121 * gh_280[k]
                  + f_63 * gh_289[k];
    }

#pragma omp simd aligned(gh_84, gh_87, gh_94, gh_231, gh_234, gh_241, gh_273, gh_276, \
                         gh_283 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_43[k] = -f_90 * gh_84[k]
                  + f_89 * gh_87[k]
                  - f_88 * gh_94[k]
                  - f_90 * gh_231[k]
                  + f_89 * gh_234[k]
                  - f_88 * gh_241[k]
                  + f_92 * gh_273[k]
                  - f_91 * gh_276[k]
                  + f_15 * gh_283[k];
    }

#pragma omp simd aligned(gh_1, gh_6, gh_15, gh_64, gh_69, gh_78, gh_106, gh_111, gh_120, \
                         gh_211, gh_216, gh_225, gh_253, gh_258, gh_267, gh_295, gh_300, \
                         gh_309 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_44[k] = f_122 * gh_1[k]
                  - f_123 * gh_6[k]
                  + f_124 * gh_15[k]
                  + f_123 * gh_64[k]
                  - f_119 * gh_69[k]
                  + f_125 * gh_78[k]
                  - f_62 * gh_106[k]
                  + f_93 * gh_111[k]
                  - f_126 * gh_120[k]
                  + f_122 * gh_211[k]
                  - f_123 * gh_216[k]
                  + f_124 * gh_225[k]
                  - f_62 * gh_253[k]
                  + f_93 * gh_258[k]
                  - f_126 * gh_267[k]
                  + f_65 * gh_295[k]
                  - f_63 * gh_300[k]
                  + f_127 * gh_309[k];
    }

#pragma omp simd aligned(gh_4, gh_11, gh_67, gh_74, gh_109, gh_116, gh_214, gh_221, gh_256, \
                         gh_263, gh_298, gh_305 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_45[k] = f_128 * gh_4[k]
                  - f_128 * gh_11[k]
                  + f_129 * gh_67[k]
                  - f_129 * gh_74[k]
                  - f_130 * gh_109[k]
                  + f_130 * gh_116[k]
                  + f_128 * gh_214[k]
                  - f_128 * gh_221[k]
                  - f_130 * gh_256[k]
                  + f_130 * gh_263[k]
                  + f_131 * gh_298[k]
                  - f_131 * gh_305[k];
    }

#pragma omp simd aligned(gh_1, gh_6, gh_8, gh_15, gh_17, gh_64, gh_69, gh_71, gh_78, gh_80, \
                         gh_106, gh_111, gh_113, gh_120, gh_122, gh_211, gh_216, gh_218, \
                         gh_225, gh_227, gh_253, gh_258, gh_260, gh_267, gh_269, gh_295, \
                         gh_300, gh_302, gh_309, gh_311 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_46[k] = -f_132 * gh_1[k]
                  - f_133 * gh_6[k]
                  + f_54 * gh_8[k]
                  + f_134 * gh_15[k]
                  - f_135 * gh_17[k]
                  - f_136 * gh_64[k]
                  - f_52 * gh_69[k]
                  + f_137 * gh_71[k]
                  + f_133 * gh_78[k]
                  - f_138 * gh_80[k]
                  + f_54 * gh_106[k]
                  + f_138 * gh_111[k]
                  - f_139 * gh_113[k]
                  - f_135 * gh_120[k]
                  + f_140 * gh_122[k]
                  - f_132 * gh_211[k]
                  - f_133 * gh_216[k]
                  + f_54 * gh_218[k]
                  + f_134 * gh_225[k]
                  - f_135 * gh_227[k]
                  + f_54 * gh_253[k]
                  + f_138 * gh_258[k]
                  - f_139 * gh_260[k]
                  - f_135 * gh_267[k]
                  + f_140 * gh_269[k]
                  - f_135 * gh_295[k]
                  - f_141 * gh_300[k]
                  + f_140 * gh_302[k]
                  + f_142 * gh_309[k]
                  - f_143 * gh_311[k];
    }

#pragma omp simd aligned(gh_4, gh_11, gh_13, gh_67, gh_74, gh_76, gh_109, gh_116, gh_118, \
                         gh_214, gh_221, gh_223, gh_256, gh_263, gh_265, gh_298, gh_305, \
                         gh_307 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_47[k] = -f_144 * gh_4[k]
                  - f_144 * gh_11[k]
                  + f_145 * gh_13[k]
                  - f_145 * gh_67[k]
                  - f_145 * gh_74[k]
                  + f_146 * gh_76[k]
                  + f_147 * gh_109[k]
                  + f_147 * gh_116[k]
                  - f_148 * gh_118[k]
                  - f_144 * gh_214[k]
                  - f_144 * gh_221[k]
                  + f_145 * gh_223[k]
                  + f_147 * gh_256[k]
                  + f_147 * gh_263[k]
                  - f_148 * gh_265[k]
                  - f_149 * gh_298[k]
                  - f_149 * gh_305[k]
                  + f_150 * gh_307[k];
    }

#pragma omp simd aligned(gh_1, gh_6, gh_8, gh_15, gh_17, gh_19, gh_64, gh_69, gh_71, gh_78, \
                         gh_80, gh_82, gh_106, gh_111, gh_113, gh_120, gh_122, gh_124, gh_211, \
                         gh_216, gh_218, gh_225, gh_227, gh_229, gh_253, gh_258, gh_260, \
                         gh_267, gh_269, gh_271, gh_295, gh_300, gh_302, gh_309, gh_311, \
                         gh_313 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_48[k] = f_151 * gh_1[k]
                  + f_152 * gh_6[k]
                  - f_153 * gh_8[k]
                  + f_151 * gh_15[k]
                  - f_153 * gh_17[k]
                  + f_154 * gh_19[k]
                  + f_152 * gh_64[k]
                  + f_155 * gh_69[k]
                  - f_156 * gh_71[k]
                  + f_152 * gh_78[k]
                  - f_156 * gh_80[k]
                  + f_157 * gh_82[k]
                  - f_154 * gh_106[k]
                  - f_157 * gh_111[k]
                  + f_158 * gh_113[k]
                  - f_154 * gh_120[k]
                  + f_158 * gh_122[k]
                  - f_159 * gh_124[k]
                  + f_151 * gh_211[k]
                  + f_152 * gh_216[k]
                  - f_153 * gh_218[k]
                  + f_151 * gh_225[k]
                  - f_153 * gh_227[k]
                  + f_154 * gh_229[k]
                  - f_154 * gh_253[k]
                  - f_157 * gh_258[k]
                  + f_158 * gh_260[k]
                  - f_154 * gh_267[k]
                  + f_158 * gh_269[k]
                  - f_159 * gh_271[k]
                  + f_160 * gh_295[k]
                  + f_161 * gh_300[k]
                  - f_162 * gh_302[k]
                  + f_160 * gh_309[k]
                  - f_162 * gh_311[k]
                  + f_163 * gh_313[k];
    }

#pragma omp simd aligned(gh_2, gh_7, gh_9, gh_16, gh_18, gh_20, gh_65, gh_70, gh_72, gh_79, \
                         gh_81, gh_83, gh_107, gh_112, gh_114, gh_121, gh_123, gh_125, gh_212, \
                         gh_217, gh_219, gh_226, gh_228, gh_230, gh_254, gh_259, gh_261, \
                         gh_268, gh_270, gh_272, gh_296, gh_301, gh_303, gh_310, gh_312, \
                         gh_314 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_49[k] = 0.703125 * gh_2[k]
                  + 1.40625 * gh_7[k]
                  - 1.875 * gh_9[k]
                  + 0.703125 * gh_16[k]
                  - 1.875 * gh_18[k]
                  + 0.375 * gh_20[k]
                  + 1.40625 * gh_65[k]
                  + 2.8125 * gh_70[k]
                  - 3.75 * gh_72[k]
                  + 1.40625 * gh_79[k]
                  - 3.75 * gh_81[k]
                  + 0.75 * gh_83[k]
                  - 5.625 * gh_107[k]
                  - 11.25 * gh_112[k]
                  + 15.0 * gh_114[k]
                  - 5.625 * gh_121[k]
                  + 15.0 * gh_123[k]
                  - 3.0 * gh_125[k]
                  + 0.703125 * gh_212[k]
                  + 1.40625 * gh_217[k]
                  - 1.875 * gh_219[k]
                  + 0.703125 * gh_226[k]
                  - 1.875 * gh_228[k]
                  + 0.375 * gh_230[k]
                  - 5.625 * gh_254[k]
                  - 11.25 * gh_259[k]
                  + 15.0 * gh_261[k]
                  - 5.625 * gh_268[k]
                  + 15.0 * gh_270[k]
                  - 3.0 * gh_272[k]
                  + 1.875 * gh_296[k]
                  + 3.75 * gh_301[k]
                  - 5.0 * gh_303[k]
                  + 1.875 * gh_310[k]
                  - 5.0 * gh_312[k]
                  + gh_314[k];
    }

#pragma omp simd aligned(gh_0, gh_3, gh_5, gh_10, gh_12, gh_14, gh_63, gh_66, gh_68, gh_73, \
                         gh_75, gh_77, gh_105, gh_108, gh_110, gh_115, gh_117, gh_119, gh_210, \
                         gh_213, gh_215, gh_220, gh_222, gh_224, gh_252, gh_255, gh_257, \
                         gh_262, gh_264, gh_266, gh_294, gh_297, gh_299, gh_304, gh_306, \
                         gh_308 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_50[k] = f_151 * gh_0[k]
                  + f_152 * gh_3[k]
                  - f_153 * gh_5[k]
                  + f_151 * gh_10[k]
                  - f_153 * gh_12[k]
                  + f_154 * gh_14[k]
                  + f_152 * gh_63[k]
                  + f_155 * gh_66[k]
                  - f_156 * gh_68[k]
                  + f_152 * gh_73[k]
                  - f_156 * gh_75[k]
                  + f_157 * gh_77[k]
                  - f_154 * gh_105[k]
                  - f_157 * gh_108[k]
                  + f_158 * gh_110[k]
                  - f_154 * gh_115[k]
                  + f_158 * gh_117[k]
                  - f_159 * gh_119[k]
                  + f_151 * gh_210[k]
                  + f_152 * gh_213[k]
                  - f_153 * gh_215[k]
                  + f_151 * gh_220[k]
                  - f_153 * gh_222[k]
                  + f_154 * gh_224[k]
                  - f_154 * gh_252[k]
                  - f_157 * gh_255[k]
                  + f_158 * gh_257[k]
                  - f_154 * gh_262[k]
                  + f_158 * gh_264[k]
                  - f_159 * gh_266[k]
                  + f_160 * gh_294[k]
                  + f_161 * gh_297[k]
                  - f_162 * gh_299[k]
                  + f_160 * gh_304[k]
                  - f_162 * gh_306[k]
                  + f_163 * gh_308[k];
    }

#pragma omp simd aligned(gh_2, gh_9, gh_16, gh_18, gh_65, gh_72, gh_79, gh_81, gh_107, gh_114, \
                         gh_121, gh_123, gh_212, gh_219, gh_226, gh_228, gh_254, gh_261, \
                         gh_268, gh_270, gh_296, gh_303, gh_310, \
                         gh_312 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_51[k] = -f_164 * gh_2[k]
                  + f_144 * gh_9[k]
                  + f_164 * gh_16[k]
                  - f_144 * gh_18[k]
                  - f_144 * gh_65[k]
                  + f_145 * gh_72[k]
                  + f_144 * gh_79[k]
                  - f_145 * gh_81[k]
                  + f_146 * gh_107[k]
                  - f_147 * gh_114[k]
                  - f_146 * gh_121[k]
                  + f_147 * gh_123[k]
                  - f_164 * gh_212[k]
                  + f_144 * gh_219[k]
                  + f_164 * gh_226[k]
                  - f_144 * gh_228[k]
                  + f_146 * gh_254[k]
                  - f_147 * gh_261[k]
                  - f_146 * gh_268[k]
                  + f_147 * gh_270[k]
                  - f_165 * gh_296[k]
                  + f_149 * gh_303[k]
                  + f_165 * gh_310[k]
                  - f_149 * gh_312[k];
    }

#pragma omp simd aligned(gh_0, gh_3, gh_5, gh_10, gh_12, gh_63, gh_66, gh_68, gh_73, gh_75, \
                         gh_105, gh_108, gh_110, gh_115, gh_117, gh_210, gh_213, gh_215, \
                         gh_220, gh_222, gh_252, gh_255, gh_257, gh_262, gh_264, gh_294, \
                         gh_297, gh_299, gh_304, gh_306 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_52[k] = -f_134 * gh_0[k]
                  + f_133 * gh_3[k]
                  + f_135 * gh_5[k]
                  + f_132 * gh_10[k]
                  - f_54 * gh_12[k]
                  - f_133 * gh_63[k]
                  + f_52 * gh_66[k]
                  + f_138 * gh_68[k]
                  + f_136 * gh_73[k]
                  - f_137 * gh_75[k]
                  + f_135 * gh_105[k]
                  - f_138 * gh_108[k]
                  - f_140 * gh_110[k]
                  - f_54 * gh_115[k]
                  + f_139 * gh_117[k]
                  - f_134 * gh_210[k]
                  + f_133 * gh_213[k]
                  + f_135 * gh_215[k]
                  + f_132 * gh_220[k]
                  - f_54 * gh_222[k]
                  + f_135 * gh_252[k]
                  - f_138 * gh_255[k]
                  - f_140 * gh_257[k]
                  - f_54 * gh_262[k]
                  + f_139 * gh_264[k]
                  - f_142 * gh_294[k]
                  + f_141 * gh_297[k]
                  + f_143 * gh_299[k]
                  + f_135 * gh_304[k]
                  - f_140 * gh_306[k];
    }

#pragma omp simd aligned(gh_2, gh_7, gh_16, gh_65, gh_70, gh_79, gh_107, gh_112, gh_121, \
                         gh_212, gh_217, gh_226, gh_254, gh_259, gh_268, gh_296, gh_301, \
                         gh_310 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_53[k] = f_166 * gh_2[k]
                  - f_167 * gh_7[k]
                  + f_166 * gh_16[k]
                  + f_90 * gh_65[k]
                  - f_168 * gh_70[k]
                  + f_90 * gh_79[k]
                  - f_129 * gh_107[k]
                  + f_169 * gh_112[k]
                  - f_129 * gh_121[k]
                  + f_166 * gh_212[k]
                  - f_167 * gh_217[k]
                  + f_166 * gh_226[k]
                  - f_129 * gh_254[k]
                  + f_169 * gh_259[k]
                  - f_129 * gh_268[k]
                  + f_92 * gh_296[k]
                  - f_170 * gh_301[k]
                  + f_92 * gh_310[k];
    }

#pragma omp simd aligned(gh_0, gh_3, gh_10, gh_63, gh_66, gh_73, gh_105, gh_108, gh_115, \
                         gh_210, gh_213, gh_220, gh_252, gh_255, gh_262, gh_294, gh_297, \
                         gh_304 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_54[k] = f_124 * gh_0[k]
                  - f_123 * gh_3[k]
                  + f_122 * gh_10[k]
                  + f_125 * gh_63[k]
                  - f_119 * gh_66[k]
                  + f_123 * gh_73[k]
                  - f_126 * gh_105[k]
                  + f_93 * gh_108[k]
                  - f_62 * gh_115[k]
                  + f_124 * gh_210[k]
                  - f_123 * gh_213[k]
                  + f_122 * gh_220[k]
                  - f_126 * gh_252[k]
                  + f_93 * gh_255[k]
                  - f_62 * gh_262[k]
                  + f_127 * gh_294[k]
                  - f_63 * gh_297[k]
                  + f_65 * gh_304[k];
    }

#pragma omp simd aligned(gh_43, gh_48, gh_57, gh_148, gh_153, gh_162, gh_190, gh_195, \
                         gh_204 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_55[k] = -f_88 * gh_43[k]
                  + f_89 * gh_48[k]
                  - f_90 * gh_57[k]
                  - f_88 * gh_148[k]
                  + f_89 * gh_153[k]
                  - f_90 * gh_162[k]
                  + f_15 * gh_190[k]
                  - f_91 * gh_195[k]
                  + f_92 * gh_204[k];
    }

#pragma omp simd aligned(gh_46, gh_53, gh_151, gh_158, gh_193, gh_200 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_56[k] = -f_93 * gh_46[k]
                  + f_93 * gh_53[k]
                  - f_93 * gh_151[k]
                  + f_93 * gh_158[k]
                  + f_66 * gh_193[k]
                  - f_66 * gh_200[k];
    }

#pragma omp simd aligned(gh_43, gh_48, gh_50, gh_57, gh_59, gh_148, gh_153, gh_155, gh_162, \
                         gh_164, gh_190, gh_195, gh_197, gh_204, \
                         gh_206 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_57[k] = f_94 * gh_43[k]
                  + f_85 * gh_48[k]
                  - f_95 * gh_50[k]
                  - f_96 * gh_57[k]
                  + f_55 * gh_59[k]
                  + f_94 * gh_148[k]
                  + f_85 * gh_153[k]
                  - f_95 * gh_155[k]
                  - f_96 * gh_162[k]
                  + f_55 * gh_164[k]
                  - f_97 * gh_190[k]
                  - f_98 * gh_195[k]
                  + f_99 * gh_197[k]
                  + f_100 * gh_204[k]
                  - f_101 * gh_206[k];
    }

#pragma omp simd aligned(gh_46, gh_53, gh_55, gh_151, gh_158, gh_160, gh_193, gh_200, \
                         gh_202 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_58[k] = f_37 * gh_46[k]
                  + f_37 * gh_53[k]
                  - f_34 * gh_55[k]
                  + f_37 * gh_151[k]
                  + f_37 * gh_158[k]
                  - f_34 * gh_160[k]
                  - f_102 * gh_193[k]
                  - f_102 * gh_200[k]
                  + f_103 * gh_202[k];
    }

#pragma omp simd aligned(gh_43, gh_48, gh_50, gh_57, gh_59, gh_61, gh_148, gh_153, gh_155, \
                         gh_162, gh_164, gh_166, gh_190, gh_195, gh_197, gh_204, gh_206, \
                         gh_208 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_59[k] = -f_104 * gh_43[k]
                  - f_105 * gh_48[k]
                  + f_106 * gh_50[k]
                  - f_104 * gh_57[k]
                  + f_106 * gh_59[k]
                  - f_107 * gh_61[k]
                  - f_104 * gh_148[k]
                  - f_105 * gh_153[k]
                  + f_106 * gh_155[k]
                  - f_104 * gh_162[k]
                  + f_106 * gh_164[k]
                  - f_107 * gh_166[k]
                  + f_108 * gh_190[k]
                  + f_109 * gh_195[k]
                  - f_110 * gh_197[k]
                  + f_108 * gh_204[k]
                  - f_110 * gh_206[k]
                  + f_111 * gh_208[k];
    }

#pragma omp simd aligned(gh_44, gh_49, gh_51, gh_58, gh_60, gh_62, gh_149, gh_154, gh_156, \
                         gh_163, gh_165, gh_167, gh_191, gh_196, gh_198, gh_205, gh_207, \
                         gh_209 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_60[k] = -f_112 * gh_44[k]
                  - f_113 * gh_49[k]
                  + f_114 * gh_51[k]
                  - f_112 * gh_58[k]
                  + f_114 * gh_60[k]
                  - f_115 * gh_62[k]
                  - f_112 * gh_149[k]
                  - f_113 * gh_154[k]
                  + f_114 * gh_156[k]
                  - f_112 * gh_163[k]
                  + f_114 * gh_165[k]
                  - f_115 * gh_167[k]
                  + f_116 * gh_191[k]
                  + f_114 * gh_196[k]
                  - f_117 * gh_198[k]
                  + f_116 * gh_205[k]
                  - f_117 * gh_207[k]
                  + f_118 * gh_209[k];
    }

#pragma omp simd aligned(gh_42, gh_45, gh_47, gh_52, gh_54, gh_56, gh_147, gh_150, gh_152, \
                         gh_157, gh_159, gh_161, gh_189, gh_192, gh_194, gh_199, gh_201, \
                         gh_203 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_61[k] = -f_104 * gh_42[k]
                  - f_105 * gh_45[k]
                  + f_106 * gh_47[k]
                  - f_104 * gh_52[k]
                  + f_106 * gh_54[k]
                  - f_107 * gh_56[k]
                  - f_104 * gh_147[k]
                  - f_105 * gh_150[k]
                  + f_106 * gh_152[k]
                  - f_104 * gh_157[k]
                  + f_106 * gh_159[k]
                  - f_107 * gh_161[k]
                  + f_108 * gh_189[k]
                  + f_109 * gh_192[k]
                  - f_110 * gh_194[k]
                  + f_108 * gh_199[k]
                  - f_110 * gh_201[k]
                  + f_111 * gh_203[k];
    }

#pragma omp simd aligned(gh_44, gh_51, gh_58, gh_60, gh_149, gh_156, gh_163, gh_165, gh_191, \
                         gh_198, gh_205, gh_207 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_62[k] = f_32 * gh_44[k]
                  - f_37 * gh_51[k]
                  - f_32 * gh_58[k]
                  + f_37 * gh_60[k]
                  + f_32 * gh_149[k]
                  - f_37 * gh_156[k]
                  - f_32 * gh_163[k]
                  + f_37 * gh_165[k]
                  - f_38 * gh_191[k]
                  + f_102 * gh_198[k]
                  + f_38 * gh_205[k]
                  - f_102 * gh_207[k];
    }

#pragma omp simd aligned(gh_42, gh_45, gh_47, gh_52, gh_54, gh_147, gh_150, gh_152, gh_157, \
                         gh_159, gh_189, gh_192, gh_194, gh_199, \
                         gh_201 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_63[k] = f_96 * gh_42[k]
                  - f_85 * gh_45[k]
                  - f_55 * gh_47[k]
                  - f_94 * gh_52[k]
                  + f_95 * gh_54[k]
                  + f_96 * gh_147[k]
                  - f_85 * gh_150[k]
                  - f_55 * gh_152[k]
                  - f_94 * gh_157[k]
                  + f_95 * gh_159[k]
                  - f_100 * gh_189[k]
                  + f_98 * gh_192[k]
                  + f_101 * gh_194[k]
                  + f_97 * gh_199[k]
                  - f_99 * gh_201[k];
    }

#pragma omp simd aligned(gh_44, gh_49, gh_58, gh_149, gh_154, gh_163, gh_191, gh_196, \
                         gh_205 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_64[k] = -f_119 * gh_44[k]
                  + f_120 * gh_49[k]
                  - f_119 * gh_58[k]
                  - f_119 * gh_149[k]
                  + f_120 * gh_154[k]
                  - f_119 * gh_163[k]
                  + f_63 * gh_191[k]
                  - f_121 * gh_196[k]
                  + f_63 * gh_205[k];
    }

#pragma omp simd aligned(gh_42, gh_45, gh_52, gh_147, gh_150, gh_157, gh_189, gh_192, \
                         gh_199 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_65[k] = -f_90 * gh_42[k]
                  + f_89 * gh_45[k]
                  - f_88 * gh_52[k]
                  - f_90 * gh_147[k]
                  + f_89 * gh_150[k]
                  - f_88 * gh_157[k]
                  + f_92 * gh_189[k]
                  - f_91 * gh_192[k]
                  + f_15 * gh_199[k];
    }

#pragma omp simd aligned(gh_1, gh_6, gh_15, gh_106, gh_111, gh_120, gh_211, gh_216, gh_225, \
                         gh_253, gh_258, gh_267 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_66[k] = -f_171 * gh_1[k]
                  + f_43 * gh_6[k]
                  - f_133 * gh_15[k]
                  + f_39 * gh_106[k]
                  - f_40 * gh_111[k]
                  + f_172 * gh_120[k]
                  + f_171 * gh_211[k]
                  - f_43 * gh_216[k]
                  + f_133 * gh_225[k]
                  - f_39 * gh_253[k]
                  + f_40 * gh_258[k]
                  - f_172 * gh_267[k];
    }

#pragma omp simd aligned(gh_4, gh_11, gh_109, gh_116, gh_214, gh_221, gh_256, \
                         gh_263 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_67[k] = -f_97 * gh_4[k]
                  + f_97 * gh_11[k]
                  + f_95 * gh_109[k]
                  - f_95 * gh_116[k]
                  + f_97 * gh_214[k]
                  - f_97 * gh_221[k]
                  - f_95 * gh_256[k]
                  + f_95 * gh_263[k];
    }

#pragma omp simd aligned(gh_1, gh_6, gh_8, gh_15, gh_17, gh_106, gh_111, gh_113, gh_120, \
                         gh_122, gh_211, gh_216, gh_218, gh_225, gh_227, gh_253, gh_258, \
                         gh_260, gh_267, gh_269 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_68[k] = f_173 * gh_1[k]
                  + f_60 * gh_6[k]
                  - f_63 * gh_8[k]
                  - f_174 * gh_15[k]
                  + f_175 * gh_17[k]
                  - f_119 * gh_106[k]
                  - f_65 * gh_111[k]
                  + f_121 * gh_113[k]
                  + f_57 * gh_120[k]
                  - f_59 * gh_122[k]
                  - f_173 * gh_211[k]
                  - f_60 * gh_216[k]
                  + f_63 * gh_218[k]
                  + f_174 * gh_225[k]
                  - f_175 * gh_227[k]
                  + f_119 * gh_253[k]
                  + f_65 * gh_258[k]
                  - f_121 * gh_260[k]
                  - f_57 * gh_267[k]
                  + f_59 * gh_269[k];
    }

#pragma omp simd aligned(gh_4, gh_11, gh_13, gh_109, gh_116, gh_118, gh_214, gh_221, gh_223, \
                         gh_256, gh_263, gh_265 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_69[k] = f_11 * gh_4[k]
                  + f_11 * gh_11[k]
                  - f_67 * gh_13[k]
                  - f_12 * gh_109[k]
                  - f_12 * gh_116[k]
                  + f_68 * gh_118[k]
                  - f_11 * gh_214[k]
                  - f_11 * gh_221[k]
                  + f_67 * gh_223[k]
                  + f_12 * gh_256[k]
                  + f_12 * gh_263[k]
                  - f_68 * gh_265[k];
    }

#pragma omp simd aligned(gh_1, gh_6, gh_8, gh_15, gh_17, gh_19, gh_106, gh_111, gh_113, \
                         gh_120, gh_122, gh_124, gh_211, gh_216, gh_218, gh_225, gh_227, \
                         gh_229, gh_253, gh_258, gh_260, gh_267, gh_269, \
                         gh_271 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_70[k] = -f_176 * gh_1[k]
                  - f_70 * gh_6[k]
                  + f_74 * gh_8[k]
                  - f_176 * gh_15[k]
                  + f_74 * gh_17[k]
                  - f_177 * gh_19[k]
                  + f_178 * gh_106[k]
                  + f_74 * gh_111[k]
                  - f_179 * gh_113[k]
                  + f_178 * gh_120[k]
                  - f_179 * gh_122[k]
                  + f_180 * gh_124[k]
                  + f_176 * gh_211[k]
                  + f_70 * gh_216[k]
                  - f_74 * gh_218[k]
                  + f_176 * gh_225[k]
                  - f_74 * gh_227[k]
                  + f_177 * gh_229[k]
                  - f_178 * gh_253[k]
                  - f_74 * gh_258[k]
                  + f_179 * gh_260[k]
                  - f_178 * gh_267[k]
                  + f_179 * gh_269[k]
                  - f_180 * gh_271[k];
    }

#pragma omp simd aligned(gh_2, gh_7, gh_9, gh_16, gh_18, gh_20, gh_107, gh_112, gh_114, \
                         gh_121, gh_123, gh_125, gh_212, gh_217, gh_219, gh_226, gh_228, \
                         gh_230, gh_254, gh_259, gh_261, gh_268, gh_270, \
                         gh_272 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_71[k] = -f_181 * gh_2[k]
                  - f_77 * gh_7[k]
                  + f_182 * gh_9[k]
                  - f_181 * gh_16[k]
                  + f_182 * gh_18[k]
                  - f_183 * gh_20[k]
                  + f_184 * gh_107[k]
                  + f_81 * gh_112[k]
                  - f_185 * gh_114[k]
                  + f_184 * gh_121[k]
                  - f_185 * gh_123[k]
                  + f_186 * gh_125[k]
                  + f_181 * gh_212[k]
                  + f_77 * gh_217[k]
                  - f_182 * gh_219[k]
                  + f_181 * gh_226[k]
                  - f_182 * gh_228[k]
                  + f_183 * gh_230[k]
                  - f_184 * gh_254[k]
                  - f_81 * gh_259[k]
                  + f_185 * gh_261[k]
                  - f_184 * gh_268[k]
                  + f_185 * gh_270[k]
                  - f_186 * gh_272[k];
    }

#pragma omp simd aligned(gh_0, gh_3, gh_5, gh_10, gh_12, gh_14, gh_105, gh_108, gh_110, \
                         gh_115, gh_117, gh_119, gh_210, gh_213, gh_215, gh_220, gh_222, \
                         gh_224, gh_252, gh_255, gh_257, gh_262, gh_264, \
                         gh_266 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_72[k] = -f_176 * gh_0[k]
                  - f_70 * gh_3[k]
                  + f_74 * gh_5[k]
                  - f_176 * gh_10[k]
                  + f_74 * gh_12[k]
                  - f_177 * gh_14[k]
                  + f_178 * gh_105[k]
                  + f_74 * gh_108[k]
                  - f_179 * gh_110[k]
                  + f_178 * gh_115[k]
                  - f_179 * gh_117[k]
                  + f_180 * gh_119[k]
                  + f_176 * gh_210[k]
                  + f_70 * gh_213[k]
                  - f_74 * gh_215[k]
                  + f_176 * gh_220[k]
                  - f_74 * gh_222[k]
                  + f_177 * gh_224[k]
                  - f_178 * gh_252[k]
                  - f_74 * gh_255[k]
                  + f_179 * gh_257[k]
                  - f_178 * gh_262[k]
                  + f_179 * gh_264[k]
                  - f_180 * gh_266[k];
    }

#pragma omp simd aligned(gh_2, gh_9, gh_16, gh_18, gh_107, gh_114, gh_121, gh_123, gh_212, \
                         gh_219, gh_226, gh_228, gh_254, gh_261, gh_268, \
                         gh_270 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_73[k] = f_10 * gh_2[k]
                  - f_11 * gh_9[k]
                  - f_10 * gh_16[k]
                  + f_11 * gh_18[k]
                  - f_187 * gh_107[k]
                  + f_12 * gh_114[k]
                  + f_187 * gh_121[k]
                  - f_12 * gh_123[k]
                  - f_10 * gh_212[k]
                  + f_11 * gh_219[k]
                  + f_10 * gh_226[k]
                  - f_11 * gh_228[k]
                  + f_187 * gh_254[k]
                  - f_12 * gh_261[k]
                  - f_187 * gh_268[k]
                  + f_12 * gh_270[k];
    }

#pragma omp simd aligned(gh_0, gh_3, gh_5, gh_10, gh_12, gh_105, gh_108, gh_110, gh_115, \
                         gh_117, gh_210, gh_213, gh_215, gh_220, gh_222, gh_252, gh_255, \
                         gh_257, gh_262, gh_264 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_74[k] = f_174 * gh_0[k]
                  - f_60 * gh_3[k]
                  - f_175 * gh_5[k]
                  - f_173 * gh_10[k]
                  + f_63 * gh_12[k]
                  - f_57 * gh_105[k]
                  + f_65 * gh_108[k]
                  + f_59 * gh_110[k]
                  + f_119 * gh_115[k]
                  - f_121 * gh_117[k]
                  - f_174 * gh_210[k]
                  + f_60 * gh_213[k]
                  + f_175 * gh_215[k]
                  + f_173 * gh_220[k]
                  - f_63 * gh_222[k]
                  + f_57 * gh_252[k]
                  - f_65 * gh_255[k]
                  - f_59 * gh_257[k]
                  - f_119 * gh_262[k]
                  + f_121 * gh_264[k];
    }

#pragma omp simd aligned(gh_2, gh_7, gh_16, gh_107, gh_112, gh_121, gh_212, gh_217, gh_226, \
                         gh_254, gh_259, gh_268 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_75[k] = -f_96 * gh_2[k]
                  + f_188 * gh_7[k]
                  - f_96 * gh_16[k]
                  + f_188 * gh_107[k]
                  - f_189 * gh_112[k]
                  + f_188 * gh_121[k]
                  + f_96 * gh_212[k]
                  - f_188 * gh_217[k]
                  + f_96 * gh_226[k]
                  - f_188 * gh_254[k]
                  + f_189 * gh_259[k]
                  - f_188 * gh_268[k];
    }

#pragma omp simd aligned(gh_0, gh_3, gh_10, gh_105, gh_108, gh_115, gh_210, gh_213, gh_220, \
                         gh_252, gh_255, gh_262 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_76[k] = -f_133 * gh_0[k]
                  + f_43 * gh_3[k]
                  - f_171 * gh_10[k]
                  + f_172 * gh_105[k]
                  - f_40 * gh_108[k]
                  + f_39 * gh_115[k]
                  + f_133 * gh_210[k]
                  - f_43 * gh_213[k]
                  + f_171 * gh_220[k]
                  - f_172 * gh_252[k]
                  + f_40 * gh_255[k]
                  - f_39 * gh_262[k];
    }

#pragma omp simd aligned(gh_43, gh_46, gh_48, gh_53, gh_57, gh_148, gh_151, gh_153, gh_158, \
                         gh_162 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_77[k] = f_22 * gh_43[k]
                  - f_23 * gh_48[k]
                  + f_24 * gh_57[k]
                  - f_19 * gh_148[k]
                  + f_20 * gh_153[k]
                  - f_21 * gh_162[k];

        g_78[k] = f_26 * gh_46[k]
                  - f_26 * gh_53[k]
                  - f_25 * gh_151[k]
                  + f_25 * gh_158[k];
    }

#pragma omp simd aligned(gh_43, gh_48, gh_50, gh_57, gh_59, gh_148, gh_153, gh_155, gh_162, \
                         gh_164 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_79[k] = -3.28125 * gh_43[k]
                  - 2.1875 * gh_48[k]
                  + 26.25 * gh_50[k]
                  + 1.09375 * gh_57[k]
                  - 8.75 * gh_59[k]
                  + 9.84375 * gh_148[k]
                  + 6.5625 * gh_153[k]
                  - 78.75 * gh_155[k]
                  - 3.28125 * gh_162[k]
                  + 26.25 * gh_164[k];
    }

#pragma omp simd aligned(gh_46, gh_53, gh_55, gh_151, gh_158, gh_160 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_80[k] = -f_29 * gh_46[k]
                  - f_29 * gh_53[k]
                  + f_30 * gh_55[k]
                  + f_27 * gh_151[k]
                  + f_27 * gh_158[k]
                  - f_28 * gh_160[k];
    }

#pragma omp simd aligned(gh_43, gh_48, gh_50, gh_57, gh_59, gh_61, gh_148, gh_153, gh_155, \
                         gh_162, gh_164, gh_166 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_81[k] = f_35 * gh_43[k]
                  + f_36 * gh_48[k]
                  - f_37 * gh_50[k]
                  + f_35 * gh_57[k]
                  - f_37 * gh_59[k]
                  + f_38 * gh_61[k]
                  - f_31 * gh_148[k]
                  - f_32 * gh_153[k]
                  + f_33 * gh_155[k]
                  - f_31 * gh_162[k]
                  + f_33 * gh_164[k]
                  - f_34 * gh_166[k];
    }

#pragma omp simd aligned(gh_44, gh_49, gh_51, gh_58, gh_60, gh_62, gh_149, gh_154, gh_156, \
                         gh_163, gh_165, gh_167 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_82[k] = f_43 * gh_44[k]
                  + f_44 * gh_49[k]
                  - f_45 * gh_51[k]
                  + f_43 * gh_58[k]
                  - f_45 * gh_60[k]
                  + f_46 * gh_62[k]
                  - f_39 * gh_149[k]
                  - f_40 * gh_154[k]
                  + f_41 * gh_156[k]
                  - f_39 * gh_163[k]
                  + f_41 * gh_165[k]
                  - f_42 * gh_167[k];
    }

#pragma omp simd aligned(gh_42, gh_45, gh_47, gh_52, gh_54, gh_56, gh_147, gh_150, gh_152, \
                         gh_157, gh_159, gh_161 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_83[k] = f_35 * gh_42[k]
                  + f_36 * gh_45[k]
                  - f_37 * gh_47[k]
                  + f_35 * gh_52[k]
                  - f_37 * gh_54[k]
                  + f_38 * gh_56[k]
                  - f_31 * gh_147[k]
                  - f_32 * gh_150[k]
                  + f_33 * gh_152[k]
                  - f_31 * gh_157[k]
                  + f_33 * gh_159[k]
                  - f_34 * gh_161[k];
    }

#pragma omp simd aligned(gh_44, gh_51, gh_58, gh_60, gh_149, gh_156, gh_163, \
                         gh_165 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_84[k] = -f_48 * gh_44[k]
                  + f_29 * gh_51[k]
                  + f_48 * gh_58[k]
                  - f_29 * gh_60[k]
                  + f_47 * gh_149[k]
                  - f_27 * gh_156[k]
                  - f_47 * gh_163[k]
                  + f_27 * gh_165[k];
    }

#pragma omp simd aligned(gh_42, gh_45, gh_47, gh_52, gh_54, gh_147, gh_150, gh_152, gh_157, \
                         gh_159 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_85[k] = -1.09375 * gh_42[k]
                  + 2.1875 * gh_45[k]
                  + 8.75 * gh_47[k]
                  + 3.28125 * gh_52[k]
                  - 26.25 * gh_54[k]
                  + 3.28125 * gh_147[k]
                  - 6.5625 * gh_150[k]
                  - 26.25 * gh_152[k]
                  - 9.84375 * gh_157[k]
                  + 78.75 * gh_159[k];
    }

#pragma omp simd aligned(gh_42, gh_44, gh_45, gh_49, gh_52, gh_58, gh_147, gh_149, gh_150, \
                         gh_154, gh_157, gh_163 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_86[k] = f_3 * gh_44[k]
                  - f_51 * gh_49[k]
                  + f_3 * gh_58[k]
                  - f_49 * gh_149[k]
                  + f_50 * gh_154[k]
                  - f_49 * gh_163[k];

        g_87[k] = f_24 * gh_42[k]
                  - f_23 * gh_45[k]
                  + f_22 * gh_52[k]
                  - f_21 * gh_147[k]
                  + f_20 * gh_150[k]
                  - f_19 * gh_157[k];
    }

#pragma omp simd aligned(gh_1, gh_6, gh_15, gh_64, gh_69, gh_78, gh_211, gh_216, \
                         gh_225 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_88[k] = f_190 * gh_1[k]
                  - f_191 * gh_6[k]
                  + f_192 * gh_15[k]
                  - f_193 * gh_64[k]
                  + f_194 * gh_69[k]
                  - f_195 * gh_78[k]
                  + f_190 * gh_211[k]
                  - f_191 * gh_216[k]
                  + f_192 * gh_225[k];
    }

#pragma omp simd aligned(gh_4, gh_11, gh_67, gh_74, gh_214, gh_221 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_89[k] = 6.5625 * gh_4[k]
                  - 6.5625 * gh_11[k]
                  - 39.375 * gh_67[k]
                  + 39.375 * gh_74[k]
                  + 6.5625 * gh_214[k]
                  - 6.5625 * gh_221[k];
    }

#pragma omp simd aligned(gh_1, gh_6, gh_8, gh_15, gh_17, gh_64, gh_69, gh_71, gh_78, gh_80, \
                         gh_211, gh_216, gh_218, gh_225, gh_227 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_90[k] = -f_196 * gh_1[k]
                  - f_197 * gh_6[k]
                  + f_198 * gh_8[k]
                  + f_199 * gh_15[k]
                  - f_4 * gh_17[k]
                  + f_200 * gh_64[k]
                  + f_3 * gh_69[k]
                  - f_25 * gh_71[k]
                  - f_201 * gh_78[k]
                  + f_26 * gh_80[k]
                  - f_196 * gh_211[k]
                  - f_197 * gh_216[k]
                  + f_198 * gh_218[k]
                  + f_199 * gh_225[k]
                  - f_4 * gh_227[k];
    }

#pragma omp simd aligned(gh_4, gh_11, gh_13, gh_67, gh_74, gh_76, gh_214, gh_221, \
                         gh_223 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_91[k] = -f_202 * gh_4[k]
                  - f_202 * gh_11[k]
                  + f_18 * gh_13[k]
                  + f_203 * gh_67[k]
                  + f_203 * gh_74[k]
                  - f_204 * gh_76[k]
                  - f_202 * gh_214[k]
                  - f_202 * gh_221[k]
                  + f_18 * gh_223[k];
    }

#pragma omp simd aligned(gh_1, gh_6, gh_8, gh_15, gh_17, gh_19, gh_64, gh_69, gh_71, gh_78, \
                         gh_80, gh_82, gh_211, gh_216, gh_218, gh_225, gh_227, \
                         gh_229 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_92[k] = f_205 * gh_1[k]
                  + f_206 * gh_6[k]
                  - f_207 * gh_8[k]
                  + f_205 * gh_15[k]
                  - f_207 * gh_17[k]
                  + f_11 * gh_19[k]
                  - f_208 * gh_64[k]
                  - f_207 * gh_69[k]
                  + f_209 * gh_71[k]
                  - f_208 * gh_78[k]
                  + f_209 * gh_80[k]
                  - f_12 * gh_82[k]
                  + f_205 * gh_211[k]
                  + f_206 * gh_216[k]
                  - f_207 * gh_218[k]
                  + f_205 * gh_225[k]
                  - f_207 * gh_227[k]
                  + f_11 * gh_229[k];
    }

#pragma omp simd aligned(gh_2, gh_7, gh_9, gh_16, gh_18, gh_20, gh_65, gh_70, gh_72, gh_79, \
                         gh_81, gh_83, gh_212, gh_217, gh_219, gh_226, gh_228, \
                         gh_230 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_93[k] = f_210 * gh_2[k]
                  + f_211 * gh_7[k]
                  - f_212 * gh_9[k]
                  + f_210 * gh_16[k]
                  - f_212 * gh_18[k]
                  + f_213 * gh_20[k]
                  - f_88 * gh_65[k]
                  - f_89 * gh_70[k]
                  + f_91 * gh_72[k]
                  - f_88 * gh_79[k]
                  + f_91 * gh_81[k]
                  - f_214 * gh_83[k]
                  + f_210 * gh_212[k]
                  + f_211 * gh_217[k]
                  - f_212 * gh_219[k]
                  + f_210 * gh_226[k]
                  - f_212 * gh_228[k]
                  + f_213 * gh_230[k];
    }

#pragma omp simd aligned(gh_0, gh_3, gh_5, gh_10, gh_12, gh_14, gh_63, gh_66, gh_68, gh_73, \
                         gh_75, gh_77, gh_210, gh_213, gh_215, gh_220, gh_222, \
                         gh_224 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_94[k] = f_205 * gh_0[k]
                  + f_206 * gh_3[k]
                  - f_207 * gh_5[k]
                  + f_205 * gh_10[k]
                  - f_207 * gh_12[k]
                  + f_11 * gh_14[k]
                  - f_208 * gh_63[k]
                  - f_207 * gh_66[k]
                  + f_209 * gh_68[k]
                  - f_208 * gh_73[k]
                  + f_209 * gh_75[k]
                  - f_12 * gh_77[k]
                  + f_205 * gh_210[k]
                  + f_206 * gh_213[k]
                  - f_207 * gh_215[k]
                  + f_205 * gh_220[k]
                  - f_207 * gh_222[k]
                  + f_11 * gh_224[k];
    }

#pragma omp simd aligned(gh_2, gh_9, gh_16, gh_18, gh_65, gh_72, gh_79, gh_81, gh_212, gh_219, \
                         gh_226, gh_228 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_95[k] = -f_215 * gh_2[k]
                  + f_202 * gh_9[k]
                  + f_215 * gh_16[k]
                  - f_202 * gh_18[k]
                  + f_216 * gh_65[k]
                  - f_203 * gh_72[k]
                  - f_216 * gh_79[k]
                  + f_203 * gh_81[k]
                  - f_215 * gh_212[k]
                  + f_202 * gh_219[k]
                  + f_215 * gh_226[k]
                  - f_202 * gh_228[k];
    }

#pragma omp simd aligned(gh_0, gh_3, gh_5, gh_10, gh_12, gh_63, gh_66, gh_68, gh_73, gh_75, \
                         gh_210, gh_213, gh_215, gh_220, gh_222 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_96[k] = -f_199 * gh_0[k]
                  + f_197 * gh_3[k]
                  + f_4 * gh_5[k]
                  + f_196 * gh_10[k]
                  - f_198 * gh_12[k]
                  + f_201 * gh_63[k]
                  - f_3 * gh_66[k]
                  - f_26 * gh_68[k]
                  - f_200 * gh_73[k]
                  + f_25 * gh_75[k]
                  - f_199 * gh_210[k]
                  + f_197 * gh_213[k]
                  + f_4 * gh_215[k]
                  + f_196 * gh_220[k]
                  - f_198 * gh_222[k];
    }

#pragma omp simd aligned(gh_2, gh_7, gh_16, gh_65, gh_70, gh_79, gh_212, gh_217, \
                         gh_226 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_97[k] = 1.640625 * gh_2[k]
                  - 9.84375 * gh_7[k]
                  + 1.640625 * gh_16[k]
                  - 9.84375 * gh_65[k]
                  + 59.0625 * gh_70[k]
                  - 9.84375 * gh_79[k]
                  + 1.640625 * gh_212[k]
                  - 9.84375 * gh_217[k]
                  + 1.640625 * gh_226[k];
    }

#pragma omp simd aligned(gh_0, gh_3, gh_10, gh_63, gh_66, gh_73, gh_210, gh_213, \
                         gh_220 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_98[k] = f_192 * gh_0[k]
                  - f_191 * gh_3[k]
                  + f_190 * gh_10[k]
                  - f_195 * gh_63[k]
                  + f_194 * gh_66[k]
                  - f_193 * gh_73[k]
                  + f_192 * gh_210[k]
                  - f_191 * gh_213[k]
                  + f_190 * gh_220[k];
    }
}

}  // namespace simdtrf
