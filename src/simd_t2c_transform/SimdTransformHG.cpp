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


#include "SimdTransformHG.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_hg(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t hg,
             const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 3.28125 * std::sqrt(10.0);
    const auto f_1 = 6.5625 * std::sqrt(10.0);
    const auto f_2 = 0.65625 * std::sqrt(10.0);
    const auto f_3 = 9.84375 * std::sqrt(5.0);
    const auto f_4 = 3.28125 * std::sqrt(5.0);
    const auto f_5 = 19.6875 * std::sqrt(5.0);
    const auto f_6 = 6.5625 * std::sqrt(5.0);
    const auto f_7 = 1.96875 * std::sqrt(5.0);
    const auto f_8 = 0.65625 * std::sqrt(5.0);
    const auto f_9 = 0.46875 * std::sqrt(70.0);
    const auto f_10 = 2.8125 * std::sqrt(70.0);
    const auto f_11 = 0.9375 * std::sqrt(70.0);
    const auto f_12 = 5.625 * std::sqrt(70.0);
    const auto f_13 = 0.09375 * std::sqrt(70.0);
    const auto f_14 = 0.5625 * std::sqrt(70.0);
    const auto f_15 = 1.40625 * std::sqrt(35.0);
    const auto f_16 = 1.875 * std::sqrt(35.0);
    const auto f_17 = 2.8125 * std::sqrt(35.0);
    const auto f_18 = 3.75 * std::sqrt(35.0);
    const auto f_19 = 0.28125 * std::sqrt(35.0);
    const auto f_20 = 0.375 * std::sqrt(35.0);
    const auto f_21 = 0.3515625 * std::sqrt(14.0);
    const auto f_22 = 0.703125 * std::sqrt(14.0);
    const auto f_23 = 2.8125 * std::sqrt(14.0);
    const auto f_24 = 0.9375 * std::sqrt(14.0);
    const auto f_25 = 1.40625 * std::sqrt(14.0);
    const auto f_26 = 5.625 * std::sqrt(14.0);
    const auto f_27 = 1.875 * std::sqrt(14.0);
    const auto f_28 = 0.0703125 * std::sqrt(14.0);
    const auto f_29 = 0.140625 * std::sqrt(14.0);
    const auto f_30 = 0.5625 * std::sqrt(14.0);
    const auto f_31 = 0.1875 * std::sqrt(14.0);
    const auto f_32 = 0.234375 * std::sqrt(70.0);
    const auto f_33 = 1.40625 * std::sqrt(70.0);
    const auto f_34 = 0.046875 * std::sqrt(70.0);
    const auto f_35 = 0.28125 * std::sqrt(70.0);
    const auto f_36 = 0.8203125 * std::sqrt(10.0);
    const auto f_37 = 4.921875 * std::sqrt(10.0);
    const auto f_38 = 1.640625 * std::sqrt(10.0);
    const auto f_39 = 9.84375 * std::sqrt(10.0);
    const auto f_40 = 0.1640625 * std::sqrt(10.0);
    const auto f_41 = 0.984375 * std::sqrt(10.0);
    const auto f_42 = 39.375 * std::sqrt(2.0);
    const auto f_43 = 13.125 * std::sqrt(2.0);
    const auto f_44 = 3.75 * std::sqrt(7.0);
    const auto f_45 = 22.5 * std::sqrt(7.0);
    const auto f_46 = 7.5 * std::sqrt(14.0);
    const auto f_47 = 0.5625 * std::sqrt(35.0);
    const auto f_48 = 1.125 * std::sqrt(35.0);
    const auto f_49 = 4.5 * std::sqrt(35.0);
    const auto f_50 = 1.5 * std::sqrt(35.0);
    const auto f_51 = 1.875 * std::sqrt(7.0);
    const auto f_52 = 11.25 * std::sqrt(7.0);
    const auto f_53 = 3.28125 * std::sqrt(2.0);
    const auto f_54 = 2.1875 * std::sqrt(2.0);
    const auto f_55 = 26.25 * std::sqrt(2.0);
    const auto f_56 = 1.09375 * std::sqrt(2.0);
    const auto f_57 = 8.75 * std::sqrt(2.0);
    const auto f_58 = 0.46875 * std::sqrt(14.0);
    const auto f_59 = 0.3125 * std::sqrt(14.0);
    const auto f_60 = 3.75 * std::sqrt(14.0);
    const auto f_61 = 22.5 * std::sqrt(14.0);
    const auto f_62 = 0.15625 * std::sqrt(14.0);
    const auto f_63 = 1.25 * std::sqrt(14.0);
    const auto f_64 = 1.40625 * std::sqrt(7.0);
    const auto f_65 = 0.9375 * std::sqrt(7.0);
    const auto f_66 = 1.25 * std::sqrt(7.0);
    const auto f_67 = 15.0 * std::sqrt(7.0);
    const auto f_68 = 0.46875 * std::sqrt(7.0);
    const auto f_69 = 0.625 * std::sqrt(7.0);
    const auto f_70 = 5.0 * std::sqrt(7.0);
    const auto f_71 = 0.0703125 * std::sqrt(70.0);
    const auto f_72 = 0.140625 * std::sqrt(70.0);
    const auto f_73 = 0.1875 * std::sqrt(70.0);
    const auto f_74 = 0.375 * std::sqrt(70.0);
    const auto f_75 = 0.125 * std::sqrt(70.0);
    const auto f_76 = 1.125 * std::sqrt(70.0);
    const auto f_77 = 4.5 * std::sqrt(70.0);
    const auto f_78 = 1.5 * std::sqrt(70.0);
    const auto f_79 = 0.0234375 * std::sqrt(70.0);
    const auto f_80 = 0.0625 * std::sqrt(70.0);
    const auto f_81 = 0.5 * std::sqrt(70.0);
    const auto f_82 = 0.234375 * std::sqrt(14.0);
    const auto f_83 = 11.25 * std::sqrt(14.0);
    const auto f_84 = 0.078125 * std::sqrt(14.0);
    const auto f_85 = 0.625 * std::sqrt(14.0);
    const auto f_86 = 0.8203125 * std::sqrt(2.0);
    const auto f_87 = 4.921875 * std::sqrt(2.0);
    const auto f_88 = 0.546875 * std::sqrt(2.0);
    const auto f_89 = 6.5625 * std::sqrt(2.0);
    const auto f_90 = 0.2734375 * std::sqrt(2.0);
    const auto f_91 = 1.640625 * std::sqrt(2.0);
    const auto f_92 = 8.75 * std::sqrt(3.0);
    const auto f_93 = 17.5 * std::sqrt(3.0);
    const auto f_94 = 13.125 * std::sqrt(6.0);
    const auto f_95 = 4.375 * std::sqrt(6.0);
    const auto f_96 = 26.25 * std::sqrt(6.0);
    const auto f_97 = 8.75 * std::sqrt(6.0);
    const auto f_98 = 1.25 * std::sqrt(21.0);
    const auto f_99 = 7.5 * std::sqrt(21.0);
    const auto f_100 = 2.5 * std::sqrt(21.0);
    const auto f_101 = 15.0 * std::sqrt(21.0);
    const auto f_102 = 1.875 * std::sqrt(42.0);
    const auto f_103 = 2.5 * std::sqrt(42.0);
    const auto f_104 = 3.75 * std::sqrt(42.0);
    const auto f_105 = 5.0 * std::sqrt(42.0);
    const auto f_106 = 0.1875 * std::sqrt(105.0);
    const auto f_107 = 0.375 * std::sqrt(105.0);
    const auto f_108 = 1.5 * std::sqrt(105.0);
    const auto f_109 = 0.5 * std::sqrt(105.0);
    const auto f_110 = 0.75 * std::sqrt(105.0);
    const auto f_111 = 3.0 * std::sqrt(105.0);
    const auto f_112 = std::sqrt(105.0);
    const auto f_113 = 0.625 * std::sqrt(21.0);
    const auto f_114 = 3.75 * std::sqrt(21.0);
    const auto f_115 = 2.1875 * std::sqrt(3.0);
    const auto f_116 = 13.125 * std::sqrt(3.0);
    const auto f_117 = 4.375 * std::sqrt(3.0);
    const auto f_118 = 26.25 * std::sqrt(3.0);
    const auto f_119 = 0.3125 * std::sqrt(21.0);
    const auto f_120 = 0.46875 * std::sqrt(42.0);
    const auto f_121 = 0.15625 * std::sqrt(42.0);
    const auto f_122 = 0.9375 * std::sqrt(42.0);
    const auto f_123 = 0.3125 * std::sqrt(42.0);
    const auto f_124 = 5.625 * std::sqrt(42.0);
    const auto f_125 = 1.25 * std::sqrt(42.0);
    const auto f_126 = 0.3125 * std::sqrt(3.0);
    const auto f_127 = 1.875 * std::sqrt(3.0);
    const auto f_128 = 0.625 * std::sqrt(3.0);
    const auto f_129 = 3.75 * std::sqrt(3.0);
    const auto f_130 = 22.5 * std::sqrt(3.0);
    const auto f_131 = 2.5 * std::sqrt(3.0);
    const auto f_132 = 15.0 * std::sqrt(3.0);
    const auto f_133 = 0.46875 * std::sqrt(6.0);
    const auto f_134 = 0.625 * std::sqrt(6.0);
    const auto f_135 = 0.9375 * std::sqrt(6.0);
    const auto f_136 = 1.25 * std::sqrt(6.0);
    const auto f_137 = 5.625 * std::sqrt(6.0);
    const auto f_138 = 7.5 * std::sqrt(6.0);
    const auto f_139 = 3.75 * std::sqrt(6.0);
    const auto f_140 = 5.0 * std::sqrt(6.0);
    const auto f_141 = 0.046875 * std::sqrt(15.0);
    const auto f_142 = 0.09375 * std::sqrt(15.0);
    const auto f_143 = 0.375 * std::sqrt(15.0);
    const auto f_144 = 0.125 * std::sqrt(15.0);
    const auto f_145 = 0.1875 * std::sqrt(15.0);
    const auto f_146 = 0.75 * std::sqrt(15.0);
    const auto f_147 = 0.25 * std::sqrt(15.0);
    const auto f_148 = 0.5625 * std::sqrt(15.0);
    const auto f_149 = 1.125 * std::sqrt(15.0);
    const auto f_150 = 4.5 * std::sqrt(15.0);
    const auto f_151 = 1.5 * std::sqrt(15.0);
    const auto f_152 = 3.0 * std::sqrt(15.0);
    const auto f_153 = std::sqrt(15.0);
    const auto f_154 = 0.15625 * std::sqrt(3.0);
    const auto f_155 = 0.9375 * std::sqrt(3.0);
    const auto f_156 = 11.25 * std::sqrt(3.0);
    const auto f_157 = 1.25 * std::sqrt(3.0);
    const auto f_158 = 7.5 * std::sqrt(3.0);
    const auto f_159 = 0.078125 * std::sqrt(21.0);
    const auto f_160 = 0.46875 * std::sqrt(21.0);
    const auto f_161 = 0.15625 * std::sqrt(21.0);
    const auto f_162 = 0.9375 * std::sqrt(21.0);
    const auto f_163 = 5.625 * std::sqrt(21.0);
    const auto f_164 = 0.9375 * std::sqrt(35.0);
    const auto f_165 = 2.5 * std::sqrt(35.0);
    const auto f_166 = 0.5 * std::sqrt(35.0);
    const auto f_167 = 3.75 * std::sqrt(70.0);
    const auto f_168 = 1.25 * std::sqrt(70.0);
    const auto f_169 = 0.75 * std::sqrt(70.0);
    const auto f_170 = 0.25 * std::sqrt(70.0);
    const auto f_171 = 0.9375 * std::sqrt(5.0);
    const auto f_172 = 5.625 * std::sqrt(5.0);
    const auto f_173 = 1.875 * std::sqrt(5.0);
    const auto f_174 = 11.25 * std::sqrt(5.0);
    const auto f_175 = 2.5 * std::sqrt(5.0);
    const auto f_176 = 15.0 * std::sqrt(5.0);
    const auto f_177 = 0.5 * std::sqrt(5.0);
    const auto f_178 = 3.0 * std::sqrt(5.0);
    const auto f_179 = 1.40625 * std::sqrt(10.0);
    const auto f_180 = 1.875 * std::sqrt(10.0);
    const auto f_181 = 2.8125 * std::sqrt(10.0);
    const auto f_182 = 3.75 * std::sqrt(10.0);
    const auto f_183 = 5.0 * std::sqrt(10.0);
    const auto f_184 = 0.75 * std::sqrt(10.0);
    const auto f_185 = std::sqrt(10.0);
    const auto f_186 = 0.46875 * std::sqrt(5.0);
    const auto f_187 = 2.8125 * std::sqrt(5.0);
    const auto f_188 = 1.25 * std::sqrt(5.0);
    const auto f_189 = 7.5 * std::sqrt(5.0);
    const auto f_190 = 0.25 * std::sqrt(5.0);
    const auto f_191 = 1.5 * std::sqrt(5.0);
    const auto f_192 = 0.234375 * std::sqrt(35.0);
    const auto f_193 = 0.46875 * std::sqrt(35.0);
    const auto f_194 = 0.625 * std::sqrt(35.0);
    const auto f_195 = 0.125 * std::sqrt(35.0);
    const auto f_196 = 0.75 * std::sqrt(35.0);
    const auto f_197 = 6.5625 * std::sqrt(6.0);
    const auto f_198 = 2.1875 * std::sqrt(6.0);
    const auto f_199 = 0.09375 * std::sqrt(105.0);
    const auto f_200 = 0.25 * std::sqrt(105.0);
    const auto f_201 = 1.875 * std::sqrt(21.0);
    const auto f_202 = 1.09375 * std::sqrt(3.0);
    const auto f_203 = 6.5625 * std::sqrt(3.0);
    const auto f_204 = 9.84375 * std::sqrt(2.0);
    const auto f_205 = 59.0625 * std::sqrt(2.0);
    const auto f_206 = 19.6875 * std::sqrt(2.0);
    const auto f_207 = 5.625 * std::sqrt(7.0);
    const auto f_208 = 33.75 * std::sqrt(7.0);
    const auto f_209 = 8.4375 * std::sqrt(14.0);
    const auto f_210 = 0.140625 * std::sqrt(35.0);
    const auto f_211 = 0.84375 * std::sqrt(35.0);
    const auto f_212 = 1.6875 * std::sqrt(35.0);
    const auto f_213 = 6.75 * std::sqrt(35.0);
    const auto f_214 = 2.25 * std::sqrt(35.0);
    const auto f_215 = 2.8125 * std::sqrt(7.0);
    const auto f_216 = 16.875 * std::sqrt(7.0);

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

    const auto *hg_0 = buffer.data(hg + 0);
    const auto *hg_1 = buffer.data(hg + 1);
    const auto *hg_2 = buffer.data(hg + 2);
    const auto *hg_3 = buffer.data(hg + 3);
    const auto *hg_4 = buffer.data(hg + 4);
    const auto *hg_5 = buffer.data(hg + 5);
    const auto *hg_6 = buffer.data(hg + 6);
    const auto *hg_7 = buffer.data(hg + 7);
    const auto *hg_8 = buffer.data(hg + 8);
    const auto *hg_9 = buffer.data(hg + 9);
    const auto *hg_10 = buffer.data(hg + 10);
    const auto *hg_11 = buffer.data(hg + 11);
    const auto *hg_12 = buffer.data(hg + 12);
    const auto *hg_13 = buffer.data(hg + 13);
    const auto *hg_14 = buffer.data(hg + 14);
    const auto *hg_15 = buffer.data(hg + 15);
    const auto *hg_16 = buffer.data(hg + 16);
    const auto *hg_17 = buffer.data(hg + 17);
    const auto *hg_18 = buffer.data(hg + 18);
    const auto *hg_19 = buffer.data(hg + 19);
    const auto *hg_20 = buffer.data(hg + 20);
    const auto *hg_21 = buffer.data(hg + 21);
    const auto *hg_22 = buffer.data(hg + 22);
    const auto *hg_23 = buffer.data(hg + 23);
    const auto *hg_24 = buffer.data(hg + 24);
    const auto *hg_25 = buffer.data(hg + 25);
    const auto *hg_26 = buffer.data(hg + 26);
    const auto *hg_27 = buffer.data(hg + 27);
    const auto *hg_28 = buffer.data(hg + 28);
    const auto *hg_29 = buffer.data(hg + 29);
    const auto *hg_30 = buffer.data(hg + 30);
    const auto *hg_31 = buffer.data(hg + 31);
    const auto *hg_32 = buffer.data(hg + 32);
    const auto *hg_33 = buffer.data(hg + 33);
    const auto *hg_34 = buffer.data(hg + 34);
    const auto *hg_35 = buffer.data(hg + 35);
    const auto *hg_36 = buffer.data(hg + 36);
    const auto *hg_37 = buffer.data(hg + 37);
    const auto *hg_38 = buffer.data(hg + 38);
    const auto *hg_39 = buffer.data(hg + 39);
    const auto *hg_40 = buffer.data(hg + 40);
    const auto *hg_41 = buffer.data(hg + 41);
    const auto *hg_42 = buffer.data(hg + 42);
    const auto *hg_43 = buffer.data(hg + 43);
    const auto *hg_44 = buffer.data(hg + 44);
    const auto *hg_45 = buffer.data(hg + 45);
    const auto *hg_46 = buffer.data(hg + 46);
    const auto *hg_47 = buffer.data(hg + 47);
    const auto *hg_48 = buffer.data(hg + 48);
    const auto *hg_49 = buffer.data(hg + 49);
    const auto *hg_50 = buffer.data(hg + 50);
    const auto *hg_51 = buffer.data(hg + 51);
    const auto *hg_52 = buffer.data(hg + 52);
    const auto *hg_53 = buffer.data(hg + 53);
    const auto *hg_54 = buffer.data(hg + 54);
    const auto *hg_55 = buffer.data(hg + 55);
    const auto *hg_56 = buffer.data(hg + 56);
    const auto *hg_57 = buffer.data(hg + 57);
    const auto *hg_58 = buffer.data(hg + 58);
    const auto *hg_59 = buffer.data(hg + 59);
    const auto *hg_60 = buffer.data(hg + 60);
    const auto *hg_61 = buffer.data(hg + 61);
    const auto *hg_62 = buffer.data(hg + 62);
    const auto *hg_63 = buffer.data(hg + 63);
    const auto *hg_64 = buffer.data(hg + 64);
    const auto *hg_65 = buffer.data(hg + 65);
    const auto *hg_66 = buffer.data(hg + 66);
    const auto *hg_67 = buffer.data(hg + 67);
    const auto *hg_68 = buffer.data(hg + 68);
    const auto *hg_69 = buffer.data(hg + 69);
    const auto *hg_70 = buffer.data(hg + 70);
    const auto *hg_71 = buffer.data(hg + 71);
    const auto *hg_72 = buffer.data(hg + 72);
    const auto *hg_73 = buffer.data(hg + 73);
    const auto *hg_74 = buffer.data(hg + 74);
    const auto *hg_75 = buffer.data(hg + 75);
    const auto *hg_76 = buffer.data(hg + 76);
    const auto *hg_77 = buffer.data(hg + 77);
    const auto *hg_78 = buffer.data(hg + 78);
    const auto *hg_79 = buffer.data(hg + 79);
    const auto *hg_80 = buffer.data(hg + 80);
    const auto *hg_81 = buffer.data(hg + 81);
    const auto *hg_82 = buffer.data(hg + 82);
    const auto *hg_83 = buffer.data(hg + 83);
    const auto *hg_84 = buffer.data(hg + 84);
    const auto *hg_85 = buffer.data(hg + 85);
    const auto *hg_86 = buffer.data(hg + 86);
    const auto *hg_87 = buffer.data(hg + 87);
    const auto *hg_88 = buffer.data(hg + 88);
    const auto *hg_89 = buffer.data(hg + 89);
    const auto *hg_90 = buffer.data(hg + 90);
    const auto *hg_91 = buffer.data(hg + 91);
    const auto *hg_92 = buffer.data(hg + 92);
    const auto *hg_93 = buffer.data(hg + 93);
    const auto *hg_94 = buffer.data(hg + 94);
    const auto *hg_95 = buffer.data(hg + 95);
    const auto *hg_96 = buffer.data(hg + 96);
    const auto *hg_97 = buffer.data(hg + 97);
    const auto *hg_98 = buffer.data(hg + 98);
    const auto *hg_99 = buffer.data(hg + 99);
    const auto *hg_100 = buffer.data(hg + 100);
    const auto *hg_101 = buffer.data(hg + 101);
    const auto *hg_102 = buffer.data(hg + 102);
    const auto *hg_103 = buffer.data(hg + 103);
    const auto *hg_104 = buffer.data(hg + 104);
    const auto *hg_105 = buffer.data(hg + 105);
    const auto *hg_106 = buffer.data(hg + 106);
    const auto *hg_107 = buffer.data(hg + 107);
    const auto *hg_108 = buffer.data(hg + 108);
    const auto *hg_109 = buffer.data(hg + 109);
    const auto *hg_110 = buffer.data(hg + 110);
    const auto *hg_111 = buffer.data(hg + 111);
    const auto *hg_112 = buffer.data(hg + 112);
    const auto *hg_113 = buffer.data(hg + 113);
    const auto *hg_114 = buffer.data(hg + 114);
    const auto *hg_115 = buffer.data(hg + 115);
    const auto *hg_116 = buffer.data(hg + 116);
    const auto *hg_117 = buffer.data(hg + 117);
    const auto *hg_118 = buffer.data(hg + 118);
    const auto *hg_119 = buffer.data(hg + 119);
    const auto *hg_120 = buffer.data(hg + 120);
    const auto *hg_121 = buffer.data(hg + 121);
    const auto *hg_122 = buffer.data(hg + 122);
    const auto *hg_123 = buffer.data(hg + 123);
    const auto *hg_124 = buffer.data(hg + 124);
    const auto *hg_125 = buffer.data(hg + 125);
    const auto *hg_126 = buffer.data(hg + 126);
    const auto *hg_127 = buffer.data(hg + 127);
    const auto *hg_128 = buffer.data(hg + 128);
    const auto *hg_129 = buffer.data(hg + 129);
    const auto *hg_130 = buffer.data(hg + 130);
    const auto *hg_131 = buffer.data(hg + 131);
    const auto *hg_132 = buffer.data(hg + 132);
    const auto *hg_133 = buffer.data(hg + 133);
    const auto *hg_134 = buffer.data(hg + 134);
    const auto *hg_135 = buffer.data(hg + 135);
    const auto *hg_136 = buffer.data(hg + 136);
    const auto *hg_137 = buffer.data(hg + 137);
    const auto *hg_138 = buffer.data(hg + 138);
    const auto *hg_139 = buffer.data(hg + 139);
    const auto *hg_140 = buffer.data(hg + 140);
    const auto *hg_141 = buffer.data(hg + 141);
    const auto *hg_142 = buffer.data(hg + 142);
    const auto *hg_143 = buffer.data(hg + 143);
    const auto *hg_144 = buffer.data(hg + 144);
    const auto *hg_145 = buffer.data(hg + 145);
    const auto *hg_146 = buffer.data(hg + 146);
    const auto *hg_147 = buffer.data(hg + 147);
    const auto *hg_148 = buffer.data(hg + 148);
    const auto *hg_149 = buffer.data(hg + 149);
    const auto *hg_150 = buffer.data(hg + 150);
    const auto *hg_151 = buffer.data(hg + 151);
    const auto *hg_152 = buffer.data(hg + 152);
    const auto *hg_153 = buffer.data(hg + 153);
    const auto *hg_154 = buffer.data(hg + 154);
    const auto *hg_155 = buffer.data(hg + 155);
    const auto *hg_156 = buffer.data(hg + 156);
    const auto *hg_157 = buffer.data(hg + 157);
    const auto *hg_158 = buffer.data(hg + 158);
    const auto *hg_159 = buffer.data(hg + 159);
    const auto *hg_160 = buffer.data(hg + 160);
    const auto *hg_161 = buffer.data(hg + 161);
    const auto *hg_162 = buffer.data(hg + 162);
    const auto *hg_163 = buffer.data(hg + 163);
    const auto *hg_164 = buffer.data(hg + 164);
    const auto *hg_165 = buffer.data(hg + 165);
    const auto *hg_166 = buffer.data(hg + 166);
    const auto *hg_167 = buffer.data(hg + 167);
    const auto *hg_168 = buffer.data(hg + 168);
    const auto *hg_169 = buffer.data(hg + 169);
    const auto *hg_170 = buffer.data(hg + 170);
    const auto *hg_171 = buffer.data(hg + 171);
    const auto *hg_172 = buffer.data(hg + 172);
    const auto *hg_173 = buffer.data(hg + 173);
    const auto *hg_174 = buffer.data(hg + 174);
    const auto *hg_175 = buffer.data(hg + 175);
    const auto *hg_176 = buffer.data(hg + 176);
    const auto *hg_177 = buffer.data(hg + 177);
    const auto *hg_178 = buffer.data(hg + 178);
    const auto *hg_179 = buffer.data(hg + 179);
    const auto *hg_180 = buffer.data(hg + 180);
    const auto *hg_181 = buffer.data(hg + 181);
    const auto *hg_182 = buffer.data(hg + 182);
    const auto *hg_183 = buffer.data(hg + 183);
    const auto *hg_184 = buffer.data(hg + 184);
    const auto *hg_185 = buffer.data(hg + 185);
    const auto *hg_186 = buffer.data(hg + 186);
    const auto *hg_187 = buffer.data(hg + 187);
    const auto *hg_188 = buffer.data(hg + 188);
    const auto *hg_189 = buffer.data(hg + 189);
    const auto *hg_190 = buffer.data(hg + 190);
    const auto *hg_191 = buffer.data(hg + 191);
    const auto *hg_192 = buffer.data(hg + 192);
    const auto *hg_193 = buffer.data(hg + 193);
    const auto *hg_194 = buffer.data(hg + 194);
    const auto *hg_195 = buffer.data(hg + 195);
    const auto *hg_196 = buffer.data(hg + 196);
    const auto *hg_197 = buffer.data(hg + 197);
    const auto *hg_198 = buffer.data(hg + 198);
    const auto *hg_199 = buffer.data(hg + 199);
    const auto *hg_200 = buffer.data(hg + 200);
    const auto *hg_201 = buffer.data(hg + 201);
    const auto *hg_202 = buffer.data(hg + 202);
    const auto *hg_203 = buffer.data(hg + 203);
    const auto *hg_204 = buffer.data(hg + 204);
    const auto *hg_205 = buffer.data(hg + 205);
    const auto *hg_206 = buffer.data(hg + 206);
    const auto *hg_207 = buffer.data(hg + 207);
    const auto *hg_208 = buffer.data(hg + 208);
    const auto *hg_209 = buffer.data(hg + 209);
    const auto *hg_210 = buffer.data(hg + 210);
    const auto *hg_211 = buffer.data(hg + 211);
    const auto *hg_212 = buffer.data(hg + 212);
    const auto *hg_213 = buffer.data(hg + 213);
    const auto *hg_214 = buffer.data(hg + 214);
    const auto *hg_215 = buffer.data(hg + 215);
    const auto *hg_216 = buffer.data(hg + 216);
    const auto *hg_217 = buffer.data(hg + 217);
    const auto *hg_218 = buffer.data(hg + 218);
    const auto *hg_219 = buffer.data(hg + 219);
    const auto *hg_220 = buffer.data(hg + 220);
    const auto *hg_221 = buffer.data(hg + 221);
    const auto *hg_222 = buffer.data(hg + 222);
    const auto *hg_223 = buffer.data(hg + 223);
    const auto *hg_224 = buffer.data(hg + 224);
    const auto *hg_225 = buffer.data(hg + 225);
    const auto *hg_226 = buffer.data(hg + 226);
    const auto *hg_227 = buffer.data(hg + 227);
    const auto *hg_228 = buffer.data(hg + 228);
    const auto *hg_229 = buffer.data(hg + 229);
    const auto *hg_230 = buffer.data(hg + 230);
    const auto *hg_231 = buffer.data(hg + 231);
    const auto *hg_232 = buffer.data(hg + 232);
    const auto *hg_233 = buffer.data(hg + 233);
    const auto *hg_234 = buffer.data(hg + 234);
    const auto *hg_235 = buffer.data(hg + 235);
    const auto *hg_236 = buffer.data(hg + 236);
    const auto *hg_237 = buffer.data(hg + 237);
    const auto *hg_238 = buffer.data(hg + 238);
    const auto *hg_239 = buffer.data(hg + 239);
    const auto *hg_240 = buffer.data(hg + 240);
    const auto *hg_241 = buffer.data(hg + 241);
    const auto *hg_242 = buffer.data(hg + 242);
    const auto *hg_243 = buffer.data(hg + 243);
    const auto *hg_244 = buffer.data(hg + 244);
    const auto *hg_245 = buffer.data(hg + 245);
    const auto *hg_246 = buffer.data(hg + 246);
    const auto *hg_247 = buffer.data(hg + 247);
    const auto *hg_248 = buffer.data(hg + 248);
    const auto *hg_249 = buffer.data(hg + 249);
    const auto *hg_250 = buffer.data(hg + 250);
    const auto *hg_251 = buffer.data(hg + 251);
    const auto *hg_252 = buffer.data(hg + 252);
    const auto *hg_253 = buffer.data(hg + 253);
    const auto *hg_254 = buffer.data(hg + 254);
    const auto *hg_255 = buffer.data(hg + 255);
    const auto *hg_256 = buffer.data(hg + 256);
    const auto *hg_257 = buffer.data(hg + 257);
    const auto *hg_258 = buffer.data(hg + 258);
    const auto *hg_259 = buffer.data(hg + 259);
    const auto *hg_260 = buffer.data(hg + 260);
    const auto *hg_261 = buffer.data(hg + 261);
    const auto *hg_262 = buffer.data(hg + 262);
    const auto *hg_263 = buffer.data(hg + 263);
    const auto *hg_264 = buffer.data(hg + 264);
    const auto *hg_265 = buffer.data(hg + 265);
    const auto *hg_266 = buffer.data(hg + 266);
    const auto *hg_267 = buffer.data(hg + 267);
    const auto *hg_268 = buffer.data(hg + 268);
    const auto *hg_269 = buffer.data(hg + 269);
    const auto *hg_270 = buffer.data(hg + 270);
    const auto *hg_271 = buffer.data(hg + 271);
    const auto *hg_272 = buffer.data(hg + 272);
    const auto *hg_273 = buffer.data(hg + 273);
    const auto *hg_274 = buffer.data(hg + 274);
    const auto *hg_275 = buffer.data(hg + 275);
    const auto *hg_276 = buffer.data(hg + 276);
    const auto *hg_277 = buffer.data(hg + 277);
    const auto *hg_278 = buffer.data(hg + 278);
    const auto *hg_279 = buffer.data(hg + 279);
    const auto *hg_280 = buffer.data(hg + 280);
    const auto *hg_281 = buffer.data(hg + 281);
    const auto *hg_282 = buffer.data(hg + 282);
    const auto *hg_283 = buffer.data(hg + 283);
    const auto *hg_284 = buffer.data(hg + 284);
    const auto *hg_285 = buffer.data(hg + 285);
    const auto *hg_286 = buffer.data(hg + 286);
    const auto *hg_287 = buffer.data(hg + 287);
    const auto *hg_288 = buffer.data(hg + 288);
    const auto *hg_289 = buffer.data(hg + 289);
    const auto *hg_290 = buffer.data(hg + 290);
    const auto *hg_291 = buffer.data(hg + 291);
    const auto *hg_292 = buffer.data(hg + 292);
    const auto *hg_293 = buffer.data(hg + 293);
    const auto *hg_294 = buffer.data(hg + 294);
    const auto *hg_295 = buffer.data(hg + 295);
    const auto *hg_296 = buffer.data(hg + 296);
    const auto *hg_297 = buffer.data(hg + 297);
    const auto *hg_298 = buffer.data(hg + 298);
    const auto *hg_299 = buffer.data(hg + 299);
    const auto *hg_300 = buffer.data(hg + 300);
    const auto *hg_301 = buffer.data(hg + 301);
    const auto *hg_302 = buffer.data(hg + 302);
    const auto *hg_303 = buffer.data(hg + 303);
    const auto *hg_304 = buffer.data(hg + 304);
    const auto *hg_305 = buffer.data(hg + 305);
    const auto *hg_306 = buffer.data(hg + 306);
    const auto *hg_307 = buffer.data(hg + 307);
    const auto *hg_308 = buffer.data(hg + 308);
    const auto *hg_309 = buffer.data(hg + 309);
    const auto *hg_310 = buffer.data(hg + 310);
    const auto *hg_311 = buffer.data(hg + 311);
    const auto *hg_312 = buffer.data(hg + 312);
    const auto *hg_313 = buffer.data(hg + 313);
    const auto *hg_314 = buffer.data(hg + 314);

#pragma omp simd aligned(hg_16, hg_19, hg_21, hg_26, hg_91, hg_94, hg_96, hg_101, hg_226, \
                         hg_229, hg_231, hg_236 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * hg_16[k]
                 - f_0 * hg_21[k]
                 - f_1 * hg_91[k]
                 + f_1 * hg_96[k]
                 + f_2 * hg_226[k]
                 - f_2 * hg_231[k];

        g_1[k] = f_3 * hg_19[k]
                 - f_4 * hg_26[k]
                 - f_5 * hg_94[k]
                 + f_6 * hg_101[k]
                 + f_7 * hg_229[k]
                 - f_8 * hg_236[k];
    }

#pragma omp simd aligned(hg_16, hg_21, hg_23, hg_91, hg_96, hg_98, hg_226, hg_231, \
                         hg_233 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_2[k] = -f_9 * hg_16[k]
                 - f_9 * hg_21[k]
                 + f_10 * hg_23[k]
                 + f_11 * hg_91[k]
                 + f_11 * hg_96[k]
                 - f_12 * hg_98[k]
                 - f_13 * hg_226[k]
                 - f_13 * hg_231[k]
                 + f_14 * hg_233[k];
    }

#pragma omp simd aligned(hg_19, hg_26, hg_28, hg_94, hg_101, hg_103, hg_229, hg_236, \
                         hg_238 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = -f_15 * hg_19[k]
                 - f_15 * hg_26[k]
                 + f_16 * hg_28[k]
                 + f_17 * hg_94[k]
                 + f_17 * hg_101[k]
                 - f_18 * hg_103[k]
                 - f_19 * hg_229[k]
                 - f_19 * hg_236[k]
                 + f_20 * hg_238[k];
    }

#pragma omp simd aligned(hg_15, hg_18, hg_20, hg_25, hg_27, hg_29, hg_90, hg_93, hg_95, \
                         hg_100, hg_102, hg_104, hg_225, hg_228, hg_230, hg_235, hg_237, \
                         hg_239 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = f_21 * hg_15[k]
                 + f_22 * hg_18[k]
                 - f_23 * hg_20[k]
                 + f_21 * hg_25[k]
                 - f_23 * hg_27[k]
                 + f_24 * hg_29[k]
                 - f_22 * hg_90[k]
                 - f_25 * hg_93[k]
                 + f_26 * hg_95[k]
                 - f_22 * hg_100[k]
                 + f_26 * hg_102[k]
                 - f_27 * hg_104[k]
                 + f_28 * hg_225[k]
                 + f_29 * hg_228[k]
                 - f_30 * hg_230[k]
                 + f_28 * hg_235[k]
                 - f_30 * hg_237[k]
                 + f_31 * hg_239[k];
    }

#pragma omp simd aligned(hg_17, hg_22, hg_24, hg_92, hg_97, hg_99, hg_227, hg_232, \
                         hg_234 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = -f_15 * hg_17[k]
                 - f_15 * hg_22[k]
                 + f_16 * hg_24[k]
                 + f_17 * hg_92[k]
                 + f_17 * hg_97[k]
                 - f_18 * hg_99[k]
                 - f_19 * hg_227[k]
                 - f_19 * hg_232[k]
                 + f_20 * hg_234[k];
    }

#pragma omp simd aligned(hg_15, hg_20, hg_25, hg_27, hg_90, hg_95, hg_100, hg_102, hg_225, \
                         hg_230, hg_235, hg_237 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = -f_32 * hg_15[k]
                 + f_33 * hg_20[k]
                 + f_32 * hg_25[k]
                 - f_33 * hg_27[k]
                 + f_9 * hg_90[k]
                 - f_10 * hg_95[k]
                 - f_9 * hg_100[k]
                 + f_10 * hg_102[k]
                 - f_34 * hg_225[k]
                 + f_35 * hg_230[k]
                 + f_34 * hg_235[k]
                 - f_35 * hg_237[k];
    }

#pragma omp simd aligned(hg_17, hg_22, hg_92, hg_97, hg_227, hg_232 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = f_4 * hg_17[k]
                 - f_3 * hg_22[k]
                 - f_6 * hg_92[k]
                 + f_5 * hg_97[k]
                 + f_8 * hg_227[k]
                 - f_7 * hg_232[k];
    }

#pragma omp simd aligned(hg_15, hg_18, hg_25, hg_61, hg_66, hg_90, hg_93, hg_100, hg_166, \
                         hg_171, hg_225, hg_228, hg_235 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_8[k] = f_36 * hg_15[k]
                 - f_37 * hg_18[k]
                 + f_36 * hg_25[k]
                 - f_38 * hg_90[k]
                 + f_39 * hg_93[k]
                 - f_38 * hg_100[k]
                 + f_40 * hg_225[k]
                 - f_41 * hg_228[k]
                 + f_40 * hg_235[k];

        g_9[k] = 26.25 * hg_61[k]
                 - 26.25 * hg_66[k]
                 - 26.25 * hg_166[k]
                 + 26.25 * hg_171[k];
    }

#pragma omp simd aligned(hg_61, hg_64, hg_66, hg_68, hg_71, hg_73, hg_166, hg_169, hg_171, \
                         hg_173, hg_176, hg_178 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = f_42 * hg_64[k]
                  - f_43 * hg_71[k]
                  - f_42 * hg_169[k]
                  + f_43 * hg_176[k];

        g_11[k] = -f_44 * hg_61[k]
                  - f_44 * hg_66[k]
                  + f_45 * hg_68[k]
                  + f_44 * hg_166[k]
                  + f_44 * hg_171[k]
                  - f_45 * hg_173[k];

        g_12[k] = -f_26 * hg_64[k]
                  - f_26 * hg_71[k]
                  + f_46 * hg_73[k]
                  + f_26 * hg_169[k]
                  + f_26 * hg_176[k]
                  - f_46 * hg_178[k];
    }

#pragma omp simd aligned(hg_60, hg_63, hg_65, hg_70, hg_72, hg_74, hg_165, hg_168, hg_170, \
                         hg_175, hg_177, hg_179 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = f_47 * hg_60[k]
                  + f_48 * hg_63[k]
                  - f_49 * hg_65[k]
                  + f_47 * hg_70[k]
                  - f_49 * hg_72[k]
                  + f_50 * hg_74[k]
                  - f_47 * hg_165[k]
                  - f_48 * hg_168[k]
                  + f_49 * hg_170[k]
                  - f_47 * hg_175[k]
                  + f_49 * hg_177[k]
                  - f_50 * hg_179[k];
    }

#pragma omp simd aligned(hg_60, hg_62, hg_65, hg_67, hg_69, hg_70, hg_72, hg_165, hg_167, \
                         hg_170, hg_172, hg_174, hg_175, hg_177 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = -f_26 * hg_62[k]
                  - f_26 * hg_67[k]
                  + f_46 * hg_69[k]
                  + f_26 * hg_167[k]
                  + f_26 * hg_172[k]
                  - f_46 * hg_174[k];

        g_15[k] = -f_51 * hg_60[k]
                  + f_52 * hg_65[k]
                  + f_51 * hg_70[k]
                  - f_52 * hg_72[k]
                  + f_51 * hg_165[k]
                  - f_52 * hg_170[k]
                  - f_51 * hg_175[k]
                  + f_52 * hg_177[k];
    }

#pragma omp simd aligned(hg_60, hg_62, hg_63, hg_67, hg_70, hg_165, hg_167, hg_168, hg_172, \
                         hg_175 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_16[k] = f_43 * hg_62[k]
                  - f_42 * hg_67[k]
                  - f_43 * hg_167[k]
                  + f_42 * hg_172[k];

        g_17[k] = 6.5625 * hg_60[k]
                  - 39.375 * hg_63[k]
                  + 6.5625 * hg_70[k]
                  - 6.5625 * hg_165[k]
                  + 39.375 * hg_168[k]
                  - 6.5625 * hg_175[k];
    }

#pragma omp simd aligned(hg_16, hg_21, hg_91, hg_96, hg_121, hg_126, hg_226, hg_231, hg_256, \
                         hg_261 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = -f_53 * hg_16[k]
                  + f_53 * hg_21[k]
                  - f_54 * hg_91[k]
                  + f_54 * hg_96[k]
                  + f_55 * hg_121[k]
                  - f_55 * hg_126[k]
                  + f_56 * hg_226[k]
                  - f_56 * hg_231[k]
                  - f_57 * hg_256[k]
                  + f_57 * hg_261[k];
    }

#pragma omp simd aligned(hg_19, hg_26, hg_94, hg_101, hg_124, hg_131, hg_229, hg_236, hg_259, \
                         hg_266 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = -9.84375 * hg_19[k]
                  + 3.28125 * hg_26[k]
                  - 6.5625 * hg_94[k]
                  + 2.1875 * hg_101[k]
                  + 78.75 * hg_124[k]
                  - 26.25 * hg_131[k]
                  + 3.28125 * hg_229[k]
                  - 1.09375 * hg_236[k]
                  - 26.25 * hg_259[k]
                  + 8.75 * hg_266[k];
    }

#pragma omp simd aligned(hg_16, hg_21, hg_23, hg_91, hg_96, hg_98, hg_121, hg_126, hg_128, \
                         hg_226, hg_231, hg_233, hg_256, hg_261, \
                         hg_263 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = f_58 * hg_16[k]
                  + f_58 * hg_21[k]
                  - f_23 * hg_23[k]
                  + f_59 * hg_91[k]
                  + f_59 * hg_96[k]
                  - f_27 * hg_98[k]
                  - f_60 * hg_121[k]
                  - f_60 * hg_126[k]
                  + f_61 * hg_128[k]
                  - f_62 * hg_226[k]
                  - f_62 * hg_231[k]
                  + f_24 * hg_233[k]
                  + f_63 * hg_256[k]
                  + f_63 * hg_261[k]
                  - f_46 * hg_263[k];
    }

#pragma omp simd aligned(hg_19, hg_26, hg_28, hg_94, hg_101, hg_103, hg_124, hg_131, hg_133, \
                         hg_229, hg_236, hg_238, hg_259, hg_266, \
                         hg_268 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = f_64 * hg_19[k]
                  + f_64 * hg_26[k]
                  - f_51 * hg_28[k]
                  + f_65 * hg_94[k]
                  + f_65 * hg_101[k]
                  - f_66 * hg_103[k]
                  - f_52 * hg_124[k]
                  - f_52 * hg_131[k]
                  + f_67 * hg_133[k]
                  - f_68 * hg_229[k]
                  - f_68 * hg_236[k]
                  + f_69 * hg_238[k]
                  + f_44 * hg_259[k]
                  + f_44 * hg_266[k]
                  - f_70 * hg_268[k];
    }

#pragma omp simd aligned(hg_15, hg_18, hg_20, hg_25, hg_27, hg_29, hg_90, hg_93, hg_95, \
                         hg_100, hg_102, hg_104, hg_120, hg_123, hg_125, hg_130, hg_132, \
                         hg_134, hg_225, hg_228, hg_230, hg_235, hg_237, hg_239, hg_255, \
                         hg_258, hg_260, hg_265, hg_267, hg_269 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = -f_71 * hg_15[k]
                  - f_72 * hg_18[k]
                  + f_14 * hg_20[k]
                  - f_71 * hg_25[k]
                  + f_14 * hg_27[k]
                  - f_73 * hg_29[k]
                  - f_34 * hg_90[k]
                  - f_13 * hg_93[k]
                  + f_74 * hg_95[k]
                  - f_34 * hg_100[k]
                  + f_74 * hg_102[k]
                  - f_75 * hg_104[k]
                  + f_14 * hg_120[k]
                  + f_76 * hg_123[k]
                  - f_77 * hg_125[k]
                  + f_14 * hg_130[k]
                  - f_77 * hg_132[k]
                  + f_78 * hg_134[k]
                  + f_79 * hg_225[k]
                  + f_34 * hg_228[k]
                  - f_73 * hg_230[k]
                  + f_79 * hg_235[k]
                  - f_73 * hg_237[k]
                  + f_80 * hg_239[k]
                  - f_73 * hg_255[k]
                  - f_74 * hg_258[k]
                  + f_78 * hg_260[k]
                  - f_73 * hg_265[k]
                  + f_78 * hg_267[k]
                  - f_81 * hg_269[k];
    }

#pragma omp simd aligned(hg_17, hg_22, hg_24, hg_92, hg_97, hg_99, hg_122, hg_127, hg_129, \
                         hg_227, hg_232, hg_234, hg_257, hg_262, \
                         hg_264 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = f_64 * hg_17[k]
                  + f_64 * hg_22[k]
                  - f_51 * hg_24[k]
                  + f_65 * hg_92[k]
                  + f_65 * hg_97[k]
                  - f_66 * hg_99[k]
                  - f_52 * hg_122[k]
                  - f_52 * hg_127[k]
                  + f_67 * hg_129[k]
                  - f_68 * hg_227[k]
                  - f_68 * hg_232[k]
                  + f_69 * hg_234[k]
                  + f_44 * hg_257[k]
                  + f_44 * hg_262[k]
                  - f_70 * hg_264[k];
    }

#pragma omp simd aligned(hg_15, hg_20, hg_25, hg_27, hg_90, hg_95, hg_100, hg_102, hg_120, \
                         hg_125, hg_130, hg_132, hg_225, hg_230, hg_235, hg_237, hg_255, \
                         hg_260, hg_265, hg_267 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = f_82 * hg_15[k]
                  - f_25 * hg_20[k]
                  - f_82 * hg_25[k]
                  + f_25 * hg_27[k]
                  + f_62 * hg_90[k]
                  - f_24 * hg_95[k]
                  - f_62 * hg_100[k]
                  + f_24 * hg_102[k]
                  - f_27 * hg_120[k]
                  + f_83 * hg_125[k]
                  + f_27 * hg_130[k]
                  - f_83 * hg_132[k]
                  - f_84 * hg_225[k]
                  + f_58 * hg_230[k]
                  + f_84 * hg_235[k]
                  - f_58 * hg_237[k]
                  + f_85 * hg_255[k]
                  - f_60 * hg_260[k]
                  - f_85 * hg_265[k]
                  + f_60 * hg_267[k];
    }

#pragma omp simd aligned(hg_17, hg_22, hg_92, hg_97, hg_122, hg_127, hg_227, hg_232, hg_257, \
                         hg_262 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = -3.28125 * hg_17[k]
                  + 9.84375 * hg_22[k]
                  - 2.1875 * hg_92[k]
                  + 6.5625 * hg_97[k]
                  + 26.25 * hg_122[k]
                  - 78.75 * hg_127[k]
                  + 1.09375 * hg_227[k]
                  - 3.28125 * hg_232[k]
                  - 8.75 * hg_257[k]
                  + 26.25 * hg_262[k];
    }

#pragma omp simd aligned(hg_15, hg_18, hg_25, hg_90, hg_93, hg_100, hg_120, hg_123, hg_130, \
                         hg_225, hg_228, hg_235, hg_255, hg_258, \
                         hg_265 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = -f_86 * hg_15[k]
                  + f_87 * hg_18[k]
                  - f_86 * hg_25[k]
                  - f_88 * hg_90[k]
                  + f_53 * hg_93[k]
                  - f_88 * hg_100[k]
                  + f_89 * hg_120[k]
                  - f_42 * hg_123[k]
                  + f_89 * hg_130[k]
                  + f_90 * hg_225[k]
                  - f_91 * hg_228[k]
                  + f_90 * hg_235[k]
                  - f_54 * hg_255[k]
                  + f_43 * hg_258[k]
                  - f_54 * hg_265[k];
    }

#pragma omp simd aligned(hg_61, hg_64, hg_66, hg_71, hg_166, hg_169, hg_171, hg_176, hg_196, \
                         hg_199, hg_201, hg_206 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = -f_92 * hg_61[k]
                  + f_92 * hg_66[k]
                  - f_92 * hg_166[k]
                  + f_92 * hg_171[k]
                  + f_93 * hg_196[k]
                  - f_93 * hg_201[k];

        g_28[k] = -f_94 * hg_64[k]
                  + f_95 * hg_71[k]
                  - f_94 * hg_169[k]
                  + f_95 * hg_176[k]
                  + f_96 * hg_199[k]
                  - f_97 * hg_206[k];
    }

#pragma omp simd aligned(hg_61, hg_66, hg_68, hg_166, hg_171, hg_173, hg_196, hg_201, \
                         hg_203 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = f_98 * hg_61[k]
                  + f_98 * hg_66[k]
                  - f_99 * hg_68[k]
                  + f_98 * hg_166[k]
                  + f_98 * hg_171[k]
                  - f_99 * hg_173[k]
                  - f_100 * hg_196[k]
                  - f_100 * hg_201[k]
                  + f_101 * hg_203[k];
    }

#pragma omp simd aligned(hg_64, hg_71, hg_73, hg_169, hg_176, hg_178, hg_199, hg_206, \
                         hg_208 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = f_102 * hg_64[k]
                  + f_102 * hg_71[k]
                  - f_103 * hg_73[k]
                  + f_102 * hg_169[k]
                  + f_102 * hg_176[k]
                  - f_103 * hg_178[k]
                  - f_104 * hg_199[k]
                  - f_104 * hg_206[k]
                  + f_105 * hg_208[k];
    }

#pragma omp simd aligned(hg_60, hg_63, hg_65, hg_70, hg_72, hg_74, hg_165, hg_168, hg_170, \
                         hg_175, hg_177, hg_179, hg_195, hg_198, hg_200, hg_205, hg_207, \
                         hg_209 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = -f_106 * hg_60[k]
                  - f_107 * hg_63[k]
                  + f_108 * hg_65[k]
                  - f_106 * hg_70[k]
                  + f_108 * hg_72[k]
                  - f_109 * hg_74[k]
                  - f_106 * hg_165[k]
                  - f_107 * hg_168[k]
                  + f_108 * hg_170[k]
                  - f_106 * hg_175[k]
                  + f_108 * hg_177[k]
                  - f_109 * hg_179[k]
                  + f_107 * hg_195[k]
                  + f_110 * hg_198[k]
                  - f_111 * hg_200[k]
                  + f_107 * hg_205[k]
                  - f_111 * hg_207[k]
                  + f_112 * hg_209[k];
    }

#pragma omp simd aligned(hg_62, hg_67, hg_69, hg_167, hg_172, hg_174, hg_197, hg_202, \
                         hg_204 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = f_102 * hg_62[k]
                  + f_102 * hg_67[k]
                  - f_103 * hg_69[k]
                  + f_102 * hg_167[k]
                  + f_102 * hg_172[k]
                  - f_103 * hg_174[k]
                  - f_104 * hg_197[k]
                  - f_104 * hg_202[k]
                  + f_105 * hg_204[k];
    }

#pragma omp simd aligned(hg_60, hg_65, hg_70, hg_72, hg_165, hg_170, hg_175, hg_177, hg_195, \
                         hg_200, hg_205, hg_207 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = f_113 * hg_60[k]
                  - f_114 * hg_65[k]
                  - f_113 * hg_70[k]
                  + f_114 * hg_72[k]
                  + f_113 * hg_165[k]
                  - f_114 * hg_170[k]
                  - f_113 * hg_175[k]
                  + f_114 * hg_177[k]
                  - f_98 * hg_195[k]
                  + f_99 * hg_200[k]
                  + f_98 * hg_205[k]
                  - f_99 * hg_207[k];
    }

#pragma omp simd aligned(hg_62, hg_67, hg_167, hg_172, hg_197, hg_202 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = -f_95 * hg_62[k]
                  + f_94 * hg_67[k]
                  - f_95 * hg_167[k]
                  + f_94 * hg_172[k]
                  + f_97 * hg_197[k]
                  - f_96 * hg_202[k];
    }

#pragma omp simd aligned(hg_60, hg_63, hg_70, hg_165, hg_168, hg_175, hg_195, hg_198, \
                         hg_205 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = -f_115 * hg_60[k]
                  + f_116 * hg_63[k]
                  - f_115 * hg_70[k]
                  - f_115 * hg_165[k]
                  + f_116 * hg_168[k]
                  - f_115 * hg_175[k]
                  + f_117 * hg_195[k]
                  - f_118 * hg_198[k]
                  + f_117 * hg_205[k];
    }

#pragma omp simd aligned(hg_16, hg_21, hg_91, hg_96, hg_121, hg_126, hg_226, hg_231, hg_256, \
                         hg_261, hg_286, hg_291 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_36[k] = f_119 * hg_16[k]
                  - f_119 * hg_21[k]
                  + f_113 * hg_91[k]
                  - f_113 * hg_96[k]
                  - f_114 * hg_121[k]
                  + f_114 * hg_126[k]
                  + f_119 * hg_226[k]
                  - f_119 * hg_231[k]
                  - f_114 * hg_256[k]
                  + f_114 * hg_261[k]
                  + f_100 * hg_286[k]
                  - f_100 * hg_291[k];
    }

#pragma omp simd aligned(hg_19, hg_26, hg_94, hg_101, hg_124, hg_131, hg_229, hg_236, hg_259, \
                         hg_266, hg_289, hg_296 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_37[k] = f_120 * hg_19[k]
                  - f_121 * hg_26[k]
                  + f_122 * hg_94[k]
                  - f_123 * hg_101[k]
                  - f_124 * hg_124[k]
                  + f_102 * hg_131[k]
                  + f_120 * hg_229[k]
                  - f_121 * hg_236[k]
                  - f_124 * hg_259[k]
                  + f_102 * hg_266[k]
                  + f_104 * hg_289[k]
                  - f_125 * hg_296[k];
    }

#pragma omp simd aligned(hg_16, hg_21, hg_23, hg_91, hg_96, hg_98, hg_121, hg_126, hg_128, \
                         hg_226, hg_231, hg_233, hg_256, hg_261, hg_263, hg_286, hg_291, \
                         hg_293 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = -f_126 * hg_16[k]
                  - f_126 * hg_21[k]
                  + f_127 * hg_23[k]
                  - f_128 * hg_91[k]
                  - f_128 * hg_96[k]
                  + f_129 * hg_98[k]
                  + f_129 * hg_121[k]
                  + f_129 * hg_126[k]
                  - f_130 * hg_128[k]
                  - f_126 * hg_226[k]
                  - f_126 * hg_231[k]
                  + f_127 * hg_233[k]
                  + f_129 * hg_256[k]
                  + f_129 * hg_261[k]
                  - f_130 * hg_263[k]
                  - f_131 * hg_286[k]
                  - f_131 * hg_291[k]
                  + f_132 * hg_293[k];
    }

#pragma omp simd aligned(hg_19, hg_26, hg_28, hg_94, hg_101, hg_103, hg_124, hg_131, hg_133, \
                         hg_229, hg_236, hg_238, hg_259, hg_266, hg_268, hg_289, hg_296, \
                         hg_298 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_39[k] = -f_133 * hg_19[k]
                  - f_133 * hg_26[k]
                  + f_134 * hg_28[k]
                  - f_135 * hg_94[k]
                  - f_135 * hg_101[k]
                  + f_136 * hg_103[k]
                  + f_137 * hg_124[k]
                  + f_137 * hg_131[k]
                  - f_138 * hg_133[k]
                  - f_133 * hg_229[k]
                  - f_133 * hg_236[k]
                  + f_134 * hg_238[k]
                  + f_137 * hg_259[k]
                  + f_137 * hg_266[k]
                  - f_138 * hg_268[k]
                  - f_139 * hg_289[k]
                  - f_139 * hg_296[k]
                  + f_140 * hg_298[k];
    }

#pragma omp simd aligned(hg_15, hg_18, hg_20, hg_25, hg_27, hg_29, hg_90, hg_93, hg_95, \
                         hg_100, hg_102, hg_104, hg_120, hg_123, hg_125, hg_130, hg_132, \
                         hg_134, hg_225, hg_228, hg_230, hg_235, hg_237, hg_239, hg_255, \
                         hg_258, hg_260, hg_265, hg_267, hg_269, hg_285, hg_288, hg_290, \
                         hg_295, hg_297, hg_299 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_40[k] = f_141 * hg_15[k]
                  + f_142 * hg_18[k]
                  - f_143 * hg_20[k]
                  + f_141 * hg_25[k]
                  - f_143 * hg_27[k]
                  + f_144 * hg_29[k]
                  + f_142 * hg_90[k]
                  + f_145 * hg_93[k]
                  - f_146 * hg_95[k]
                  + f_142 * hg_100[k]
                  - f_146 * hg_102[k]
                  + f_147 * hg_104[k]
                  - f_148 * hg_120[k]
                  - f_149 * hg_123[k]
                  + f_150 * hg_125[k]
                  - f_148 * hg_130[k]
                  + f_150 * hg_132[k]
                  - f_151 * hg_134[k]
                  + f_141 * hg_225[k]
                  + f_142 * hg_228[k]
                  - f_143 * hg_230[k]
                  + f_141 * hg_235[k]
                  - f_143 * hg_237[k]
                  + f_144 * hg_239[k]
                  - f_148 * hg_255[k]
                  - f_149 * hg_258[k]
                  + f_150 * hg_260[k]
                  - f_148 * hg_265[k]
                  + f_150 * hg_267[k]
                  - f_151 * hg_269[k]
                  + f_143 * hg_285[k]
                  + f_146 * hg_288[k]
                  - f_152 * hg_290[k]
                  + f_143 * hg_295[k]
                  - f_152 * hg_297[k]
                  + f_153 * hg_299[k];
    }

#pragma omp simd aligned(hg_17, hg_22, hg_24, hg_92, hg_97, hg_99, hg_122, hg_127, hg_129, \
                         hg_227, hg_232, hg_234, hg_257, hg_262, hg_264, hg_287, hg_292, \
                         hg_294 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_41[k] = -f_133 * hg_17[k]
                  - f_133 * hg_22[k]
                  + f_134 * hg_24[k]
                  - f_135 * hg_92[k]
                  - f_135 * hg_97[k]
                  + f_136 * hg_99[k]
                  + f_137 * hg_122[k]
                  + f_137 * hg_127[k]
                  - f_138 * hg_129[k]
                  - f_133 * hg_227[k]
                  - f_133 * hg_232[k]
                  + f_134 * hg_234[k]
                  + f_137 * hg_257[k]
                  + f_137 * hg_262[k]
                  - f_138 * hg_264[k]
                  - f_139 * hg_287[k]
                  - f_139 * hg_292[k]
                  + f_140 * hg_294[k];
    }

#pragma omp simd aligned(hg_15, hg_20, hg_25, hg_27, hg_90, hg_95, hg_100, hg_102, hg_120, \
                         hg_125, hg_130, hg_132, hg_225, hg_230, hg_235, hg_237, hg_255, \
                         hg_260, hg_265, hg_267, hg_285, hg_290, hg_295, \
                         hg_297 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = -f_154 * hg_15[k]
                  + f_155 * hg_20[k]
                  + f_154 * hg_25[k]
                  - f_155 * hg_27[k]
                  - f_126 * hg_90[k]
                  + f_127 * hg_95[k]
                  + f_126 * hg_100[k]
                  - f_127 * hg_102[k]
                  + f_127 * hg_120[k]
                  - f_156 * hg_125[k]
                  - f_127 * hg_130[k]
                  + f_156 * hg_132[k]
                  - f_154 * hg_225[k]
                  + f_155 * hg_230[k]
                  + f_154 * hg_235[k]
                  - f_155 * hg_237[k]
                  + f_127 * hg_255[k]
                  - f_156 * hg_260[k]
                  - f_127 * hg_265[k]
                  + f_156 * hg_267[k]
                  - f_157 * hg_285[k]
                  + f_158 * hg_290[k]
                  + f_157 * hg_295[k]
                  - f_158 * hg_297[k];
    }

#pragma omp simd aligned(hg_17, hg_22, hg_92, hg_97, hg_122, hg_127, hg_227, hg_232, hg_257, \
                         hg_262, hg_287, hg_292 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_43[k] = f_121 * hg_17[k]
                  - f_120 * hg_22[k]
                  + f_123 * hg_92[k]
                  - f_122 * hg_97[k]
                  - f_102 * hg_122[k]
                  + f_124 * hg_127[k]
                  + f_121 * hg_227[k]
                  - f_120 * hg_232[k]
                  - f_102 * hg_257[k]
                  + f_124 * hg_262[k]
                  + f_125 * hg_287[k]
                  - f_104 * hg_292[k];
    }

#pragma omp simd aligned(hg_15, hg_18, hg_25, hg_90, hg_93, hg_100, hg_120, hg_123, hg_130, \
                         hg_225, hg_228, hg_235, hg_255, hg_258, hg_265, hg_285, hg_288, \
                         hg_295 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_44[k] = f_159 * hg_15[k]
                  - f_160 * hg_18[k]
                  + f_159 * hg_25[k]
                  + f_161 * hg_90[k]
                  - f_162 * hg_93[k]
                  + f_161 * hg_100[k]
                  - f_162 * hg_120[k]
                  + f_163 * hg_123[k]
                  - f_162 * hg_130[k]
                  + f_159 * hg_225[k]
                  - f_160 * hg_228[k]
                  + f_159 * hg_235[k]
                  - f_162 * hg_255[k]
                  + f_163 * hg_258[k]
                  - f_162 * hg_265[k]
                  + f_113 * hg_285[k]
                  - f_114 * hg_288[k]
                  + f_113 * hg_295[k];
    }

#pragma omp simd aligned(hg_31, hg_36, hg_106, hg_111, hg_136, hg_141, hg_241, hg_246, hg_271, \
                         hg_276, hg_301, hg_306 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_45[k] = f_164 * hg_31[k]
                  - f_164 * hg_36[k]
                  + f_16 * hg_106[k]
                  - f_16 * hg_111[k]
                  - f_165 * hg_136[k]
                  + f_165 * hg_141[k]
                  + f_164 * hg_241[k]
                  - f_164 * hg_246[k]
                  - f_165 * hg_271[k]
                  + f_165 * hg_276[k]
                  + f_166 * hg_301[k]
                  - f_166 * hg_306[k];
    }

#pragma omp simd aligned(hg_34, hg_41, hg_109, hg_116, hg_139, hg_146, hg_244, hg_251, hg_274, \
                         hg_281, hg_304, hg_311 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_46[k] = f_33 * hg_34[k]
                  - f_9 * hg_41[k]
                  + f_10 * hg_109[k]
                  - f_11 * hg_116[k]
                  - f_167 * hg_139[k]
                  + f_168 * hg_146[k]
                  + f_33 * hg_244[k]
                  - f_9 * hg_251[k]
                  - f_167 * hg_274[k]
                  + f_168 * hg_281[k]
                  + f_169 * hg_304[k]
                  - f_170 * hg_311[k];
    }

#pragma omp simd aligned(hg_31, hg_36, hg_38, hg_106, hg_111, hg_113, hg_136, hg_141, hg_143, \
                         hg_241, hg_246, hg_248, hg_271, hg_276, hg_278, hg_301, hg_306, \
                         hg_308 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_47[k] = -f_171 * hg_31[k]
                  - f_171 * hg_36[k]
                  + f_172 * hg_38[k]
                  - f_173 * hg_106[k]
                  - f_173 * hg_111[k]
                  + f_174 * hg_113[k]
                  + f_175 * hg_136[k]
                  + f_175 * hg_141[k]
                  - f_176 * hg_143[k]
                  - f_171 * hg_241[k]
                  - f_171 * hg_246[k]
                  + f_172 * hg_248[k]
                  + f_175 * hg_271[k]
                  + f_175 * hg_276[k]
                  - f_176 * hg_278[k]
                  - f_177 * hg_301[k]
                  - f_177 * hg_306[k]
                  + f_178 * hg_308[k];
    }

#pragma omp simd aligned(hg_34, hg_41, hg_43, hg_109, hg_116, hg_118, hg_139, hg_146, hg_148, \
                         hg_244, hg_251, hg_253, hg_274, hg_281, hg_283, hg_304, hg_311, \
                         hg_313 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_48[k] = -f_179 * hg_34[k]
                  - f_179 * hg_41[k]
                  + f_180 * hg_43[k]
                  - f_181 * hg_109[k]
                  - f_181 * hg_116[k]
                  + f_182 * hg_118[k]
                  + f_182 * hg_139[k]
                  + f_182 * hg_146[k]
                  - f_183 * hg_148[k]
                  - f_179 * hg_244[k]
                  - f_179 * hg_251[k]
                  + f_180 * hg_253[k]
                  + f_182 * hg_274[k]
                  + f_182 * hg_281[k]
                  - f_183 * hg_283[k]
                  - f_184 * hg_304[k]
                  - f_184 * hg_311[k]
                  + f_185 * hg_313[k];
    }

#pragma omp simd aligned(hg_30, hg_33, hg_35, hg_40, hg_42, hg_44, hg_105, hg_108, hg_110, \
                         hg_115, hg_117, hg_119, hg_135, hg_138, hg_140, hg_145, hg_147, \
                         hg_149, hg_240, hg_243, hg_245, hg_250, hg_252, hg_254, hg_270, \
                         hg_273, hg_275, hg_280, hg_282, hg_284, hg_300, hg_303, hg_305, \
                         hg_310, hg_312, hg_314 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_49[k] = 0.703125 * hg_30[k]
                  + 1.40625 * hg_33[k]
                  - 5.625 * hg_35[k]
                  + 0.703125 * hg_40[k]
                  - 5.625 * hg_42[k]
                  + 1.875 * hg_44[k]
                  + 1.40625 * hg_105[k]
                  + 2.8125 * hg_108[k]
                  - 11.25 * hg_110[k]
                  + 1.40625 * hg_115[k]
                  - 11.25 * hg_117[k]
                  + 3.75 * hg_119[k]
                  - 1.875 * hg_135[k]
                  - 3.75 * hg_138[k]
                  + 15.0 * hg_140[k]
                  - 1.875 * hg_145[k]
                  + 15.0 * hg_147[k]
                  - 5.0 * hg_149[k]
                  + 0.703125 * hg_240[k]
                  + 1.40625 * hg_243[k]
                  - 5.625 * hg_245[k]
                  + 0.703125 * hg_250[k]
                  - 5.625 * hg_252[k]
                  + 1.875 * hg_254[k]
                  - 1.875 * hg_270[k]
                  - 3.75 * hg_273[k]
                  + 15.0 * hg_275[k]
                  - 1.875 * hg_280[k]
                  + 15.0 * hg_282[k]
                  - 5.0 * hg_284[k]
                  + 0.375 * hg_300[k]
                  + 0.75 * hg_303[k]
                  - 3.0 * hg_305[k]
                  + 0.375 * hg_310[k]
                  - 3.0 * hg_312[k]
                  + hg_314[k];
    }

#pragma omp simd aligned(hg_32, hg_37, hg_39, hg_107, hg_112, hg_114, hg_137, hg_142, hg_144, \
                         hg_242, hg_247, hg_249, hg_272, hg_277, hg_279, hg_302, hg_307, \
                         hg_309 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_50[k] = -f_179 * hg_32[k]
                  - f_179 * hg_37[k]
                  + f_180 * hg_39[k]
                  - f_181 * hg_107[k]
                  - f_181 * hg_112[k]
                  + f_182 * hg_114[k]
                  + f_182 * hg_137[k]
                  + f_182 * hg_142[k]
                  - f_183 * hg_144[k]
                  - f_179 * hg_242[k]
                  - f_179 * hg_247[k]
                  + f_180 * hg_249[k]
                  + f_182 * hg_272[k]
                  + f_182 * hg_277[k]
                  - f_183 * hg_279[k]
                  - f_184 * hg_302[k]
                  - f_184 * hg_307[k]
                  + f_185 * hg_309[k];
    }

#pragma omp simd aligned(hg_30, hg_35, hg_40, hg_42, hg_105, hg_110, hg_115, hg_117, hg_135, \
                         hg_140, hg_145, hg_147, hg_240, hg_245, hg_250, hg_252, hg_270, \
                         hg_275, hg_280, hg_282, hg_300, hg_305, hg_310, \
                         hg_312 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_51[k] = -f_186 * hg_30[k]
                  + f_187 * hg_35[k]
                  + f_186 * hg_40[k]
                  - f_187 * hg_42[k]
                  - f_171 * hg_105[k]
                  + f_172 * hg_110[k]
                  + f_171 * hg_115[k]
                  - f_172 * hg_117[k]
                  + f_188 * hg_135[k]
                  - f_189 * hg_140[k]
                  - f_188 * hg_145[k]
                  + f_189 * hg_147[k]
                  - f_186 * hg_240[k]
                  + f_187 * hg_245[k]
                  + f_186 * hg_250[k]
                  - f_187 * hg_252[k]
                  + f_188 * hg_270[k]
                  - f_189 * hg_275[k]
                  - f_188 * hg_280[k]
                  + f_189 * hg_282[k]
                  - f_190 * hg_300[k]
                  + f_191 * hg_305[k]
                  + f_190 * hg_310[k]
                  - f_191 * hg_312[k];
    }

#pragma omp simd aligned(hg_32, hg_37, hg_107, hg_112, hg_137, hg_142, hg_242, hg_247, hg_272, \
                         hg_277, hg_302, hg_307 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_52[k] = f_9 * hg_32[k]
                  - f_33 * hg_37[k]
                  + f_11 * hg_107[k]
                  - f_10 * hg_112[k]
                  - f_168 * hg_137[k]
                  + f_167 * hg_142[k]
                  + f_9 * hg_242[k]
                  - f_33 * hg_247[k]
                  - f_168 * hg_272[k]
                  + f_167 * hg_277[k]
                  + f_170 * hg_302[k]
                  - f_169 * hg_307[k];
    }

#pragma omp simd aligned(hg_30, hg_33, hg_40, hg_105, hg_108, hg_115, hg_135, hg_138, hg_145, \
                         hg_240, hg_243, hg_250, hg_270, hg_273, hg_280, hg_300, hg_303, \
                         hg_310 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_53[k] = f_192 * hg_30[k]
                  - f_15 * hg_33[k]
                  + f_192 * hg_40[k]
                  + f_193 * hg_105[k]
                  - f_17 * hg_108[k]
                  + f_193 * hg_115[k]
                  - f_194 * hg_135[k]
                  + f_18 * hg_138[k]
                  - f_194 * hg_145[k]
                  + f_192 * hg_240[k]
                  - f_15 * hg_243[k]
                  + f_192 * hg_250[k]
                  - f_194 * hg_270[k]
                  + f_18 * hg_273[k]
                  - f_194 * hg_280[k]
                  + f_195 * hg_300[k]
                  - f_196 * hg_303[k]
                  + f_195 * hg_310[k];
    }

#pragma omp simd aligned(hg_1, hg_6, hg_46, hg_51, hg_76, hg_81, hg_151, hg_156, hg_181, \
                         hg_186, hg_211, hg_216 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_54[k] = f_119 * hg_1[k]
                  - f_119 * hg_6[k]
                  + f_113 * hg_46[k]
                  - f_113 * hg_51[k]
                  - f_114 * hg_76[k]
                  + f_114 * hg_81[k]
                  + f_119 * hg_151[k]
                  - f_119 * hg_156[k]
                  - f_114 * hg_181[k]
                  + f_114 * hg_186[k]
                  + f_100 * hg_211[k]
                  - f_100 * hg_216[k];
    }

#pragma omp simd aligned(hg_4, hg_11, hg_49, hg_56, hg_79, hg_86, hg_154, hg_161, hg_184, \
                         hg_191, hg_214, hg_221 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_55[k] = f_120 * hg_4[k]
                  - f_121 * hg_11[k]
                  + f_122 * hg_49[k]
                  - f_123 * hg_56[k]
                  - f_124 * hg_79[k]
                  + f_102 * hg_86[k]
                  + f_120 * hg_154[k]
                  - f_121 * hg_161[k]
                  - f_124 * hg_184[k]
                  + f_102 * hg_191[k]
                  + f_104 * hg_214[k]
                  - f_125 * hg_221[k];
    }

#pragma omp simd aligned(hg_1, hg_6, hg_8, hg_46, hg_51, hg_53, hg_76, hg_81, hg_83, hg_151, \
                         hg_156, hg_158, hg_181, hg_186, hg_188, hg_211, hg_216, \
                         hg_218 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_56[k] = -f_126 * hg_1[k]
                  - f_126 * hg_6[k]
                  + f_127 * hg_8[k]
                  - f_128 * hg_46[k]
                  - f_128 * hg_51[k]
                  + f_129 * hg_53[k]
                  + f_129 * hg_76[k]
                  + f_129 * hg_81[k]
                  - f_130 * hg_83[k]
                  - f_126 * hg_151[k]
                  - f_126 * hg_156[k]
                  + f_127 * hg_158[k]
                  + f_129 * hg_181[k]
                  + f_129 * hg_186[k]
                  - f_130 * hg_188[k]
                  - f_131 * hg_211[k]
                  - f_131 * hg_216[k]
                  + f_132 * hg_218[k];
    }

#pragma omp simd aligned(hg_4, hg_11, hg_13, hg_49, hg_56, hg_58, hg_79, hg_86, hg_88, hg_154, \
                         hg_161, hg_163, hg_184, hg_191, hg_193, hg_214, hg_221, \
                         hg_223 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_57[k] = -f_133 * hg_4[k]
                  - f_133 * hg_11[k]
                  + f_134 * hg_13[k]
                  - f_135 * hg_49[k]
                  - f_135 * hg_56[k]
                  + f_136 * hg_58[k]
                  + f_137 * hg_79[k]
                  + f_137 * hg_86[k]
                  - f_138 * hg_88[k]
                  - f_133 * hg_154[k]
                  - f_133 * hg_161[k]
                  + f_134 * hg_163[k]
                  + f_137 * hg_184[k]
                  + f_137 * hg_191[k]
                  - f_138 * hg_193[k]
                  - f_139 * hg_214[k]
                  - f_139 * hg_221[k]
                  + f_140 * hg_223[k];
    }

#pragma omp simd aligned(hg_0, hg_3, hg_5, hg_10, hg_12, hg_14, hg_45, hg_48, hg_50, hg_55, \
                         hg_57, hg_59, hg_75, hg_78, hg_80, hg_85, hg_87, hg_89, hg_150, \
                         hg_153, hg_155, hg_160, hg_162, hg_164, hg_180, hg_183, hg_185, \
                         hg_190, hg_192, hg_194, hg_210, hg_213, hg_215, hg_220, hg_222, \
                         hg_224 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_58[k] = f_141 * hg_0[k]
                  + f_142 * hg_3[k]
                  - f_143 * hg_5[k]
                  + f_141 * hg_10[k]
                  - f_143 * hg_12[k]
                  + f_144 * hg_14[k]
                  + f_142 * hg_45[k]
                  + f_145 * hg_48[k]
                  - f_146 * hg_50[k]
                  + f_142 * hg_55[k]
                  - f_146 * hg_57[k]
                  + f_147 * hg_59[k]
                  - f_148 * hg_75[k]
                  - f_149 * hg_78[k]
                  + f_150 * hg_80[k]
                  - f_148 * hg_85[k]
                  + f_150 * hg_87[k]
                  - f_151 * hg_89[k]
                  + f_141 * hg_150[k]
                  + f_142 * hg_153[k]
                  - f_143 * hg_155[k]
                  + f_141 * hg_160[k]
                  - f_143 * hg_162[k]
                  + f_144 * hg_164[k]
                  - f_148 * hg_180[k]
                  - f_149 * hg_183[k]
                  + f_150 * hg_185[k]
                  - f_148 * hg_190[k]
                  + f_150 * hg_192[k]
                  - f_151 * hg_194[k]
                  + f_143 * hg_210[k]
                  + f_146 * hg_213[k]
                  - f_152 * hg_215[k]
                  + f_143 * hg_220[k]
                  - f_152 * hg_222[k]
                  + f_153 * hg_224[k];
    }

#pragma omp simd aligned(hg_2, hg_7, hg_9, hg_47, hg_52, hg_54, hg_77, hg_82, hg_84, hg_152, \
                         hg_157, hg_159, hg_182, hg_187, hg_189, hg_212, hg_217, \
                         hg_219 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_59[k] = -f_133 * hg_2[k]
                  - f_133 * hg_7[k]
                  + f_134 * hg_9[k]
                  - f_135 * hg_47[k]
                  - f_135 * hg_52[k]
                  + f_136 * hg_54[k]
                  + f_137 * hg_77[k]
                  + f_137 * hg_82[k]
                  - f_138 * hg_84[k]
                  - f_133 * hg_152[k]
                  - f_133 * hg_157[k]
                  + f_134 * hg_159[k]
                  + f_137 * hg_182[k]
                  + f_137 * hg_187[k]
                  - f_138 * hg_189[k]
                  - f_139 * hg_212[k]
                  - f_139 * hg_217[k]
                  + f_140 * hg_219[k];
    }

#pragma omp simd aligned(hg_0, hg_5, hg_10, hg_12, hg_45, hg_50, hg_55, hg_57, hg_75, hg_80, \
                         hg_85, hg_87, hg_150, hg_155, hg_160, hg_162, hg_180, hg_185, hg_190, \
                         hg_192, hg_210, hg_215, hg_220, hg_222 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_60[k] = -f_154 * hg_0[k]
                  + f_155 * hg_5[k]
                  + f_154 * hg_10[k]
                  - f_155 * hg_12[k]
                  - f_126 * hg_45[k]
                  + f_127 * hg_50[k]
                  + f_126 * hg_55[k]
                  - f_127 * hg_57[k]
                  + f_127 * hg_75[k]
                  - f_156 * hg_80[k]
                  - f_127 * hg_85[k]
                  + f_156 * hg_87[k]
                  - f_154 * hg_150[k]
                  + f_155 * hg_155[k]
                  + f_154 * hg_160[k]
                  - f_155 * hg_162[k]
                  + f_127 * hg_180[k]
                  - f_156 * hg_185[k]
                  - f_127 * hg_190[k]
                  + f_156 * hg_192[k]
                  - f_157 * hg_210[k]
                  + f_158 * hg_215[k]
                  + f_157 * hg_220[k]
                  - f_158 * hg_222[k];
    }

#pragma omp simd aligned(hg_2, hg_7, hg_47, hg_52, hg_77, hg_82, hg_152, hg_157, hg_182, \
                         hg_187, hg_212, hg_217 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_61[k] = f_121 * hg_2[k]
                  - f_120 * hg_7[k]
                  + f_123 * hg_47[k]
                  - f_122 * hg_52[k]
                  - f_102 * hg_77[k]
                  + f_124 * hg_82[k]
                  + f_121 * hg_152[k]
                  - f_120 * hg_157[k]
                  - f_102 * hg_182[k]
                  + f_124 * hg_187[k]
                  + f_125 * hg_212[k]
                  - f_104 * hg_217[k];
    }

#pragma omp simd aligned(hg_0, hg_3, hg_10, hg_45, hg_48, hg_55, hg_75, hg_78, hg_85, hg_150, \
                         hg_153, hg_160, hg_180, hg_183, hg_190, hg_210, hg_213, \
                         hg_220 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_62[k] = f_159 * hg_0[k]
                  - f_160 * hg_3[k]
                  + f_159 * hg_10[k]
                  + f_161 * hg_45[k]
                  - f_162 * hg_48[k]
                  + f_161 * hg_55[k]
                  - f_162 * hg_75[k]
                  + f_163 * hg_78[k]
                  - f_162 * hg_85[k]
                  + f_159 * hg_150[k]
                  - f_160 * hg_153[k]
                  + f_159 * hg_160[k]
                  - f_162 * hg_180[k]
                  + f_163 * hg_183[k]
                  - f_162 * hg_190[k]
                  + f_113 * hg_210[k]
                  - f_114 * hg_213[k]
                  + f_113 * hg_220[k];
    }

#pragma omp simd aligned(hg_31, hg_36, hg_136, hg_141, hg_241, hg_246, hg_271, \
                         hg_276 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_63[k] = -f_117 * hg_31[k]
                  + f_117 * hg_36[k]
                  + f_92 * hg_136[k]
                  - f_92 * hg_141[k]
                  + f_117 * hg_241[k]
                  - f_117 * hg_246[k]
                  - f_92 * hg_271[k]
                  + f_92 * hg_276[k];
    }

#pragma omp simd aligned(hg_34, hg_41, hg_139, hg_146, hg_244, hg_251, hg_274, \
                         hg_281 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_64[k] = -f_197 * hg_34[k]
                  + f_198 * hg_41[k]
                  + f_94 * hg_139[k]
                  - f_95 * hg_146[k]
                  + f_197 * hg_244[k]
                  - f_198 * hg_251[k]
                  - f_94 * hg_274[k]
                  + f_95 * hg_281[k];
    }

#pragma omp simd aligned(hg_31, hg_36, hg_38, hg_136, hg_141, hg_143, hg_241, hg_246, hg_248, \
                         hg_271, hg_276, hg_278 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_65[k] = f_113 * hg_31[k]
                  + f_113 * hg_36[k]
                  - f_114 * hg_38[k]
                  - f_98 * hg_136[k]
                  - f_98 * hg_141[k]
                  + f_99 * hg_143[k]
                  - f_113 * hg_241[k]
                  - f_113 * hg_246[k]
                  + f_114 * hg_248[k]
                  + f_98 * hg_271[k]
                  + f_98 * hg_276[k]
                  - f_99 * hg_278[k];
    }

#pragma omp simd aligned(hg_34, hg_41, hg_43, hg_139, hg_146, hg_148, hg_244, hg_251, hg_253, \
                         hg_274, hg_281, hg_283 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_66[k] = f_122 * hg_34[k]
                  + f_122 * hg_41[k]
                  - f_125 * hg_43[k]
                  - f_102 * hg_139[k]
                  - f_102 * hg_146[k]
                  + f_103 * hg_148[k]
                  - f_122 * hg_244[k]
                  - f_122 * hg_251[k]
                  + f_125 * hg_253[k]
                  + f_102 * hg_274[k]
                  + f_102 * hg_281[k]
                  - f_103 * hg_283[k];
    }

#pragma omp simd aligned(hg_30, hg_33, hg_35, hg_40, hg_42, hg_44, hg_135, hg_138, hg_140, \
                         hg_145, hg_147, hg_149, hg_240, hg_243, hg_245, hg_250, hg_252, \
                         hg_254, hg_270, hg_273, hg_275, hg_280, hg_282, \
                         hg_284 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_67[k] = -f_199 * hg_30[k]
                  - f_106 * hg_33[k]
                  + f_110 * hg_35[k]
                  - f_199 * hg_40[k]
                  + f_110 * hg_42[k]
                  - f_200 * hg_44[k]
                  + f_106 * hg_135[k]
                  + f_107 * hg_138[k]
                  - f_108 * hg_140[k]
                  + f_106 * hg_145[k]
                  - f_108 * hg_147[k]
                  + f_109 * hg_149[k]
                  + f_199 * hg_240[k]
                  + f_106 * hg_243[k]
                  - f_110 * hg_245[k]
                  + f_199 * hg_250[k]
                  - f_110 * hg_252[k]
                  + f_200 * hg_254[k]
                  - f_106 * hg_270[k]
                  - f_107 * hg_273[k]
                  + f_108 * hg_275[k]
                  - f_106 * hg_280[k]
                  + f_108 * hg_282[k]
                  - f_109 * hg_284[k];
    }

#pragma omp simd aligned(hg_32, hg_37, hg_39, hg_137, hg_142, hg_144, hg_242, hg_247, hg_249, \
                         hg_272, hg_277, hg_279 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_68[k] = f_122 * hg_32[k]
                  + f_122 * hg_37[k]
                  - f_125 * hg_39[k]
                  - f_102 * hg_137[k]
                  - f_102 * hg_142[k]
                  + f_103 * hg_144[k]
                  - f_122 * hg_242[k]
                  - f_122 * hg_247[k]
                  + f_125 * hg_249[k]
                  + f_102 * hg_272[k]
                  + f_102 * hg_277[k]
                  - f_103 * hg_279[k];
    }

#pragma omp simd aligned(hg_30, hg_35, hg_40, hg_42, hg_135, hg_140, hg_145, hg_147, hg_240, \
                         hg_245, hg_250, hg_252, hg_270, hg_275, hg_280, \
                         hg_282 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_69[k] = f_119 * hg_30[k]
                  - f_201 * hg_35[k]
                  - f_119 * hg_40[k]
                  + f_201 * hg_42[k]
                  - f_113 * hg_135[k]
                  + f_114 * hg_140[k]
                  + f_113 * hg_145[k]
                  - f_114 * hg_147[k]
                  - f_119 * hg_240[k]
                  + f_201 * hg_245[k]
                  + f_119 * hg_250[k]
                  - f_201 * hg_252[k]
                  + f_113 * hg_270[k]
                  - f_114 * hg_275[k]
                  - f_113 * hg_280[k]
                  + f_114 * hg_282[k];
    }

#pragma omp simd aligned(hg_32, hg_37, hg_137, hg_142, hg_242, hg_247, hg_272, \
                         hg_277 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_70[k] = -f_198 * hg_32[k]
                  + f_197 * hg_37[k]
                  + f_95 * hg_137[k]
                  - f_94 * hg_142[k]
                  + f_198 * hg_242[k]
                  - f_197 * hg_247[k]
                  - f_95 * hg_272[k]
                  + f_94 * hg_277[k];
    }

#pragma omp simd aligned(hg_30, hg_33, hg_40, hg_135, hg_138, hg_145, hg_240, hg_243, hg_250, \
                         hg_270, hg_273, hg_280 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_71[k] = -f_202 * hg_30[k]
                  + f_203 * hg_33[k]
                  - f_202 * hg_40[k]
                  + f_115 * hg_135[k]
                  - f_116 * hg_138[k]
                  + f_115 * hg_145[k]
                  + f_202 * hg_240[k]
                  - f_203 * hg_243[k]
                  + f_202 * hg_250[k]
                  - f_115 * hg_270[k]
                  + f_116 * hg_273[k]
                  - f_115 * hg_280[k];
    }

#pragma omp simd aligned(hg_1, hg_6, hg_46, hg_51, hg_76, hg_81, hg_151, hg_156, hg_181, \
                         hg_186 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_72[k] = -f_56 * hg_1[k]
                  + f_56 * hg_6[k]
                  + f_54 * hg_46[k]
                  - f_54 * hg_51[k]
                  + f_57 * hg_76[k]
                  - f_57 * hg_81[k]
                  + f_53 * hg_151[k]
                  - f_53 * hg_156[k]
                  - f_55 * hg_181[k]
                  + f_55 * hg_186[k];
    }

#pragma omp simd aligned(hg_4, hg_11, hg_49, hg_56, hg_79, hg_86, hg_154, hg_161, hg_184, \
                         hg_191 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_73[k] = -3.28125 * hg_4[k]
                  + 1.09375 * hg_11[k]
                  + 6.5625 * hg_49[k]
                  - 2.1875 * hg_56[k]
                  + 26.25 * hg_79[k]
                  - 8.75 * hg_86[k]
                  + 9.84375 * hg_154[k]
                  - 3.28125 * hg_161[k]
                  - 78.75 * hg_184[k]
                  + 26.25 * hg_191[k];
    }

#pragma omp simd aligned(hg_1, hg_6, hg_8, hg_46, hg_51, hg_53, hg_76, hg_81, hg_83, hg_151, \
                         hg_156, hg_158, hg_181, hg_186, hg_188 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_74[k] = f_62 * hg_1[k]
                  + f_62 * hg_6[k]
                  - f_24 * hg_8[k]
                  - f_59 * hg_46[k]
                  - f_59 * hg_51[k]
                  + f_27 * hg_53[k]
                  - f_63 * hg_76[k]
                  - f_63 * hg_81[k]
                  + f_46 * hg_83[k]
                  - f_58 * hg_151[k]
                  - f_58 * hg_156[k]
                  + f_23 * hg_158[k]
                  + f_60 * hg_181[k]
                  + f_60 * hg_186[k]
                  - f_61 * hg_188[k];
    }

#pragma omp simd aligned(hg_4, hg_11, hg_13, hg_49, hg_56, hg_58, hg_79, hg_86, hg_88, hg_154, \
                         hg_161, hg_163, hg_184, hg_191, hg_193 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_75[k] = f_68 * hg_4[k]
                  + f_68 * hg_11[k]
                  - f_69 * hg_13[k]
                  - f_65 * hg_49[k]
                  - f_65 * hg_56[k]
                  + f_66 * hg_58[k]
                  - f_44 * hg_79[k]
                  - f_44 * hg_86[k]
                  + f_70 * hg_88[k]
                  - f_64 * hg_154[k]
                  - f_64 * hg_161[k]
                  + f_51 * hg_163[k]
                  + f_52 * hg_184[k]
                  + f_52 * hg_191[k]
                  - f_67 * hg_193[k];
    }

#pragma omp simd aligned(hg_0, hg_3, hg_5, hg_10, hg_12, hg_14, hg_45, hg_48, hg_50, hg_55, \
                         hg_57, hg_59, hg_75, hg_78, hg_80, hg_85, hg_87, hg_89, hg_150, \
                         hg_153, hg_155, hg_160, hg_162, hg_164, hg_180, hg_183, hg_185, \
                         hg_190, hg_192, hg_194 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_76[k] = -f_79 * hg_0[k]
                  - f_34 * hg_3[k]
                  + f_73 * hg_5[k]
                  - f_79 * hg_10[k]
                  + f_73 * hg_12[k]
                  - f_80 * hg_14[k]
                  + f_34 * hg_45[k]
                  + f_13 * hg_48[k]
                  - f_74 * hg_50[k]
                  + f_34 * hg_55[k]
                  - f_74 * hg_57[k]
                  + f_75 * hg_59[k]
                  + f_73 * hg_75[k]
                  + f_74 * hg_78[k]
                  - f_78 * hg_80[k]
                  + f_73 * hg_85[k]
                  - f_78 * hg_87[k]
                  + f_81 * hg_89[k]
                  + f_71 * hg_150[k]
                  + f_72 * hg_153[k]
                  - f_14 * hg_155[k]
                  + f_71 * hg_160[k]
                  - f_14 * hg_162[k]
                  + f_73 * hg_164[k]
                  - f_14 * hg_180[k]
                  - f_76 * hg_183[k]
                  + f_77 * hg_185[k]
                  - f_14 * hg_190[k]
                  + f_77 * hg_192[k]
                  - f_78 * hg_194[k];
    }

#pragma omp simd aligned(hg_2, hg_7, hg_9, hg_47, hg_52, hg_54, hg_77, hg_82, hg_84, hg_152, \
                         hg_157, hg_159, hg_182, hg_187, hg_189 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_77[k] = f_68 * hg_2[k]
                  + f_68 * hg_7[k]
                  - f_69 * hg_9[k]
                  - f_65 * hg_47[k]
                  - f_65 * hg_52[k]
                  + f_66 * hg_54[k]
                  - f_44 * hg_77[k]
                  - f_44 * hg_82[k]
                  + f_70 * hg_84[k]
                  - f_64 * hg_152[k]
                  - f_64 * hg_157[k]
                  + f_51 * hg_159[k]
                  + f_52 * hg_182[k]
                  + f_52 * hg_187[k]
                  - f_67 * hg_189[k];
    }

#pragma omp simd aligned(hg_0, hg_5, hg_10, hg_12, hg_45, hg_50, hg_55, hg_57, hg_75, hg_80, \
                         hg_85, hg_87, hg_150, hg_155, hg_160, hg_162, hg_180, hg_185, hg_190, \
                         hg_192 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_78[k] = f_84 * hg_0[k]
                  - f_58 * hg_5[k]
                  - f_84 * hg_10[k]
                  + f_58 * hg_12[k]
                  - f_62 * hg_45[k]
                  + f_24 * hg_50[k]
                  + f_62 * hg_55[k]
                  - f_24 * hg_57[k]
                  - f_85 * hg_75[k]
                  + f_60 * hg_80[k]
                  + f_85 * hg_85[k]
                  - f_60 * hg_87[k]
                  - f_82 * hg_150[k]
                  + f_25 * hg_155[k]
                  + f_82 * hg_160[k]
                  - f_25 * hg_162[k]
                  + f_27 * hg_180[k]
                  - f_83 * hg_185[k]
                  - f_27 * hg_190[k]
                  + f_83 * hg_192[k];
    }

#pragma omp simd aligned(hg_2, hg_7, hg_47, hg_52, hg_77, hg_82, hg_152, hg_157, hg_182, \
                         hg_187 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_79[k] = -1.09375 * hg_2[k]
                  + 3.28125 * hg_7[k]
                  + 2.1875 * hg_47[k]
                  - 6.5625 * hg_52[k]
                  + 8.75 * hg_77[k]
                  - 26.25 * hg_82[k]
                  + 3.28125 * hg_152[k]
                  - 9.84375 * hg_157[k]
                  - 26.25 * hg_182[k]
                  + 78.75 * hg_187[k];
    }

#pragma omp simd aligned(hg_0, hg_3, hg_10, hg_45, hg_48, hg_55, hg_75, hg_78, hg_85, hg_150, \
                         hg_153, hg_160, hg_180, hg_183, hg_190 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_80[k] = -f_90 * hg_0[k]
                  + f_91 * hg_3[k]
                  - f_90 * hg_10[k]
                  + f_88 * hg_45[k]
                  - f_53 * hg_48[k]
                  + f_88 * hg_55[k]
                  + f_54 * hg_75[k]
                  - f_43 * hg_78[k]
                  + f_54 * hg_85[k]
                  + f_86 * hg_150[k]
                  - f_87 * hg_153[k]
                  + f_86 * hg_160[k]
                  - f_89 * hg_180[k]
                  + f_42 * hg_183[k]
                  - f_89 * hg_190[k];
    }

#pragma omp simd aligned(hg_31, hg_34, hg_36, hg_41, hg_106, hg_109, hg_111, hg_116, hg_241, \
                         hg_244, hg_246, hg_251 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_81[k] = 6.5625 * hg_31[k]
                  - 6.5625 * hg_36[k]
                  - 39.375 * hg_106[k]
                  + 39.375 * hg_111[k]
                  + 6.5625 * hg_241[k]
                  - 6.5625 * hg_246[k];

        g_82[k] = f_204 * hg_34[k]
                  - f_53 * hg_41[k]
                  - f_205 * hg_109[k]
                  + f_206 * hg_116[k]
                  + f_204 * hg_244[k]
                  - f_53 * hg_251[k];
    }

#pragma omp simd aligned(hg_31, hg_36, hg_38, hg_106, hg_111, hg_113, hg_241, hg_246, \
                         hg_248 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_83[k] = -f_65 * hg_31[k]
                  - f_65 * hg_36[k]
                  + f_207 * hg_38[k]
                  + f_207 * hg_106[k]
                  + f_207 * hg_111[k]
                  - f_208 * hg_113[k]
                  - f_65 * hg_241[k]
                  - f_65 * hg_246[k]
                  + f_207 * hg_248[k];
    }

#pragma omp simd aligned(hg_34, hg_41, hg_43, hg_109, hg_116, hg_118, hg_244, hg_251, \
                         hg_253 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_84[k] = -f_25 * hg_34[k]
                  - f_25 * hg_41[k]
                  + f_27 * hg_43[k]
                  + f_209 * hg_109[k]
                  + f_209 * hg_116[k]
                  - f_83 * hg_118[k]
                  - f_25 * hg_244[k]
                  - f_25 * hg_251[k]
                  + f_27 * hg_253[k];
    }

#pragma omp simd aligned(hg_30, hg_33, hg_35, hg_40, hg_42, hg_44, hg_105, hg_108, hg_110, \
                         hg_115, hg_117, hg_119, hg_240, hg_243, hg_245, hg_250, hg_252, \
                         hg_254 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_85[k] = f_210 * hg_30[k]
                  + f_19 * hg_33[k]
                  - f_48 * hg_35[k]
                  + f_210 * hg_40[k]
                  - f_48 * hg_42[k]
                  + f_20 * hg_44[k]
                  - f_211 * hg_105[k]
                  - f_212 * hg_108[k]
                  + f_213 * hg_110[k]
                  - f_211 * hg_115[k]
                  + f_213 * hg_117[k]
                  - f_214 * hg_119[k]
                  + f_210 * hg_240[k]
                  + f_19 * hg_243[k]
                  - f_48 * hg_245[k]
                  + f_210 * hg_250[k]
                  - f_48 * hg_252[k]
                  + f_20 * hg_254[k];
    }

#pragma omp simd aligned(hg_32, hg_37, hg_39, hg_107, hg_112, hg_114, hg_242, hg_247, \
                         hg_249 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_86[k] = -f_25 * hg_32[k]
                  - f_25 * hg_37[k]
                  + f_27 * hg_39[k]
                  + f_209 * hg_107[k]
                  + f_209 * hg_112[k]
                  - f_83 * hg_114[k]
                  - f_25 * hg_242[k]
                  - f_25 * hg_247[k]
                  + f_27 * hg_249[k];
    }

#pragma omp simd aligned(hg_30, hg_35, hg_40, hg_42, hg_105, hg_110, hg_115, hg_117, hg_240, \
                         hg_245, hg_250, hg_252 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_87[k] = -f_68 * hg_30[k]
                  + f_215 * hg_35[k]
                  + f_68 * hg_40[k]
                  - f_215 * hg_42[k]
                  + f_215 * hg_105[k]
                  - f_216 * hg_110[k]
                  - f_215 * hg_115[k]
                  + f_216 * hg_117[k]
                  - f_68 * hg_240[k]
                  + f_215 * hg_245[k]
                  + f_68 * hg_250[k]
                  - f_215 * hg_252[k];
    }

#pragma omp simd aligned(hg_32, hg_37, hg_107, hg_112, hg_242, hg_247 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_88[k] = f_53 * hg_32[k]
                  - f_204 * hg_37[k]
                  - f_206 * hg_107[k]
                  + f_205 * hg_112[k]
                  + f_53 * hg_242[k]
                  - f_204 * hg_247[k];
    }

#pragma omp simd aligned(hg_30, hg_33, hg_40, hg_105, hg_108, hg_115, hg_240, hg_243, \
                         hg_250 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_89[k] = 1.640625 * hg_30[k]
                  - 9.84375 * hg_33[k]
                  + 1.640625 * hg_40[k]
                  - 9.84375 * hg_105[k]
                  + 59.0625 * hg_108[k]
                  - 9.84375 * hg_115[k]
                  + 1.640625 * hg_240[k]
                  - 9.84375 * hg_243[k]
                  + 1.640625 * hg_250[k];
    }

#pragma omp simd aligned(hg_1, hg_4, hg_6, hg_11, hg_46, hg_49, hg_51, hg_56, hg_151, hg_154, \
                         hg_156, hg_161 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_90[k] = f_2 * hg_1[k]
                  - f_2 * hg_6[k]
                  - f_1 * hg_46[k]
                  + f_1 * hg_51[k]
                  + f_0 * hg_151[k]
                  - f_0 * hg_156[k];

        g_91[k] = f_7 * hg_4[k]
                  - f_8 * hg_11[k]
                  - f_5 * hg_49[k]
                  + f_6 * hg_56[k]
                  + f_3 * hg_154[k]
                  - f_4 * hg_161[k];
    }

#pragma omp simd aligned(hg_1, hg_6, hg_8, hg_46, hg_51, hg_53, hg_151, hg_156, \
                         hg_158 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_92[k] = -f_13 * hg_1[k]
                  - f_13 * hg_6[k]
                  + f_14 * hg_8[k]
                  + f_11 * hg_46[k]
                  + f_11 * hg_51[k]
                  - f_12 * hg_53[k]
                  - f_9 * hg_151[k]
                  - f_9 * hg_156[k]
                  + f_10 * hg_158[k];
    }

#pragma omp simd aligned(hg_4, hg_11, hg_13, hg_49, hg_56, hg_58, hg_154, hg_161, \
                         hg_163 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_93[k] = -f_19 * hg_4[k]
                  - f_19 * hg_11[k]
                  + f_20 * hg_13[k]
                  + f_17 * hg_49[k]
                  + f_17 * hg_56[k]
                  - f_18 * hg_58[k]
                  - f_15 * hg_154[k]
                  - f_15 * hg_161[k]
                  + f_16 * hg_163[k];
    }

#pragma omp simd aligned(hg_0, hg_3, hg_5, hg_10, hg_12, hg_14, hg_45, hg_48, hg_50, hg_55, \
                         hg_57, hg_59, hg_150, hg_153, hg_155, hg_160, hg_162, \
                         hg_164 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_94[k] = f_28 * hg_0[k]
                  + f_29 * hg_3[k]
                  - f_30 * hg_5[k]
                  + f_28 * hg_10[k]
                  - f_30 * hg_12[k]
                  + f_31 * hg_14[k]
                  - f_22 * hg_45[k]
                  - f_25 * hg_48[k]
                  + f_26 * hg_50[k]
                  - f_22 * hg_55[k]
                  + f_26 * hg_57[k]
                  - f_27 * hg_59[k]
                  + f_21 * hg_150[k]
                  + f_22 * hg_153[k]
                  - f_23 * hg_155[k]
                  + f_21 * hg_160[k]
                  - f_23 * hg_162[k]
                  + f_24 * hg_164[k];
    }

#pragma omp simd aligned(hg_2, hg_7, hg_9, hg_47, hg_52, hg_54, hg_152, hg_157, \
                         hg_159 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_95[k] = -f_19 * hg_2[k]
                  - f_19 * hg_7[k]
                  + f_20 * hg_9[k]
                  + f_17 * hg_47[k]
                  + f_17 * hg_52[k]
                  - f_18 * hg_54[k]
                  - f_15 * hg_152[k]
                  - f_15 * hg_157[k]
                  + f_16 * hg_159[k];
    }

#pragma omp simd aligned(hg_0, hg_5, hg_10, hg_12, hg_45, hg_50, hg_55, hg_57, hg_150, hg_155, \
                         hg_160, hg_162 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_96[k] = -f_34 * hg_0[k]
                  + f_35 * hg_5[k]
                  + f_34 * hg_10[k]
                  - f_35 * hg_12[k]
                  + f_9 * hg_45[k]
                  - f_10 * hg_50[k]
                  - f_9 * hg_55[k]
                  + f_10 * hg_57[k]
                  - f_32 * hg_150[k]
                  + f_33 * hg_155[k]
                  + f_32 * hg_160[k]
                  - f_33 * hg_162[k];
    }

#pragma omp simd aligned(hg_2, hg_7, hg_47, hg_52, hg_152, hg_157 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_97[k] = f_8 * hg_2[k]
                  - f_7 * hg_7[k]
                  - f_6 * hg_47[k]
                  + f_5 * hg_52[k]
                  + f_4 * hg_152[k]
                  - f_3 * hg_157[k];
    }

#pragma omp simd aligned(hg_0, hg_3, hg_10, hg_45, hg_48, hg_55, hg_150, hg_153, \
                         hg_160 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_98[k] = f_40 * hg_0[k]
                  - f_41 * hg_3[k]
                  + f_40 * hg_10[k]
                  - f_38 * hg_45[k]
                  + f_39 * hg_48[k]
                  - f_38 * hg_55[k]
                  + f_36 * hg_150[k]
                  - f_37 * hg_153[k]
                  + f_36 * hg_160[k];
    }
}

}  // namespace simdtrf
