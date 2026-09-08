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


#include "SimdTransformIF.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_if(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t if_,
             const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 0.28125 * std::sqrt(1155.0);
    const auto f_1 = 0.09375 * std::sqrt(1155.0);
    const auto f_2 = 0.9375 * std::sqrt(1155.0);
    const auto f_3 = 0.3125 * std::sqrt(1155.0);
    const auto f_4 = 0.5625 * std::sqrt(770.0);
    const auto f_5 = 1.875 * std::sqrt(770.0);
    const auto f_6 = 0.28125 * std::sqrt(77.0);
    const auto f_7 = 1.125 * std::sqrt(77.0);
    const auto f_8 = 0.9375 * std::sqrt(77.0);
    const auto f_9 = 3.75 * std::sqrt(77.0);
    const auto f_10 = 0.28125 * std::sqrt(462.0);
    const auto f_11 = 0.1875 * std::sqrt(462.0);
    const auto f_12 = 0.9375 * std::sqrt(462.0);
    const auto f_13 = 0.625 * std::sqrt(462.0);
    const auto f_14 = 0.28125 * std::sqrt(770.0);
    const auto f_15 = 0.9375 * std::sqrt(770.0);
    const auto f_16 = 1.40625 * std::sqrt(385.0);
    const auto f_17 = 0.46875 * std::sqrt(385.0);
    const auto f_18 = 2.8125 * std::sqrt(385.0);
    const auto f_19 = 0.9375 * std::sqrt(385.0);
    const auto f_20 = 0.28125 * std::sqrt(385.0);
    const auto f_21 = 0.09375 * std::sqrt(385.0);
    const auto f_22 = 0.9375 * std::sqrt(2310.0);
    const auto f_23 = 1.875 * std::sqrt(2310.0);
    const auto f_24 = 0.1875 * std::sqrt(2310.0);
    const auto f_25 = 0.46875 * std::sqrt(231.0);
    const auto f_26 = 1.875 * std::sqrt(231.0);
    const auto f_27 = 0.9375 * std::sqrt(231.0);
    const auto f_28 = 3.75 * std::sqrt(231.0);
    const auto f_29 = 0.09375 * std::sqrt(231.0);
    const auto f_30 = 0.375 * std::sqrt(231.0);
    const auto f_31 = 1.40625 * std::sqrt(154.0);
    const auto f_32 = 0.9375 * std::sqrt(154.0);
    const auto f_33 = 2.8125 * std::sqrt(154.0);
    const auto f_34 = 1.875 * std::sqrt(154.0);
    const auto f_35 = 0.28125 * std::sqrt(154.0);
    const auto f_36 = 0.1875 * std::sqrt(154.0);
    const auto f_37 = 0.46875 * std::sqrt(2310.0);
    const auto f_38 = 0.09375 * std::sqrt(2310.0);
    const auto f_39 = 0.5625 * std::sqrt(70.0);
    const auto f_40 = 0.1875 * std::sqrt(70.0);
    const auto f_41 = 5.625 * std::sqrt(70.0);
    const auto f_42 = 1.875 * std::sqrt(70.0);
    const auto f_43 = 0.75 * std::sqrt(105.0);
    const auto f_44 = 7.5 * std::sqrt(105.0);
    const auto f_45 = 0.1875 * std::sqrt(42.0);
    const auto f_46 = 0.75 * std::sqrt(42.0);
    const auto f_47 = 1.875 * std::sqrt(42.0);
    const auto f_48 = 7.5 * std::sqrt(42.0);
    const auto f_49 = 1.125 * std::sqrt(7.0);
    const auto f_50 = 0.75 * std::sqrt(7.0);
    const auto f_51 = 11.25 * std::sqrt(7.0);
    const auto f_52 = 7.5 * std::sqrt(7.0);
    const auto f_53 = 0.375 * std::sqrt(105.0);
    const auto f_54 = 3.75 * std::sqrt(105.0);
    const auto f_55 = 4.21875 * std::sqrt(21.0);
    const auto f_56 = 1.40625 * std::sqrt(21.0);
    const auto f_57 = 2.8125 * std::sqrt(21.0);
    const auto f_58 = 0.9375 * std::sqrt(21.0);
    const auto f_59 = 11.25 * std::sqrt(21.0);
    const auto f_60 = 3.75 * std::sqrt(21.0);
    const auto f_61 = 0.46875 * std::sqrt(21.0);
    const auto f_62 = 1.25 * std::sqrt(21.0);
    const auto f_63 = 8.4375 * std::sqrt(14.0);
    const auto f_64 = 5.625 * std::sqrt(14.0);
    const auto f_65 = 22.5 * std::sqrt(14.0);
    const auto f_66 = 2.8125 * std::sqrt(14.0);
    const auto f_67 = 7.5 * std::sqrt(14.0);
    const auto f_68 = 0.84375 * std::sqrt(35.0);
    const auto f_69 = 3.375 * std::sqrt(35.0);
    const auto f_70 = 0.5625 * std::sqrt(35.0);
    const auto f_71 = 2.25 * std::sqrt(35.0);
    const auto f_72 = 9.0 * std::sqrt(35.0);
    const auto f_73 = 0.28125 * std::sqrt(35.0);
    const auto f_74 = 1.125 * std::sqrt(35.0);
    const auto f_75 = 0.75 * std::sqrt(35.0);
    const auto f_76 = 3.0 * std::sqrt(35.0);
    const auto f_77 = 0.84375 * std::sqrt(210.0);
    const auto f_78 = 0.5625 * std::sqrt(210.0);
    const auto f_79 = 0.375 * std::sqrt(210.0);
    const auto f_80 = 2.25 * std::sqrt(210.0);
    const auto f_81 = 1.5 * std::sqrt(210.0);
    const auto f_82 = 0.28125 * std::sqrt(210.0);
    const auto f_83 = 0.1875 * std::sqrt(210.0);
    const auto f_84 = 0.75 * std::sqrt(210.0);
    const auto f_85 = 0.5 * std::sqrt(210.0);
    const auto f_86 = 4.21875 * std::sqrt(14.0);
    const auto f_87 = 11.25 * std::sqrt(14.0);
    const auto f_88 = 1.40625 * std::sqrt(14.0);
    const auto f_89 = 3.75 * std::sqrt(14.0);
    const auto f_90 = 0.15625 * std::sqrt(21.0);
    const auto f_91 = 0.3125 * std::sqrt(21.0);
    const auto f_92 = 7.5 * std::sqrt(21.0);
    const auto f_93 = 2.5 * std::sqrt(21.0);
    const auto f_94 = 0.9375 * std::sqrt(14.0);
    const auto f_95 = 1.875 * std::sqrt(14.0);
    const auto f_96 = 15.0 * std::sqrt(14.0);
    const auto f_97 = 0.09375 * std::sqrt(35.0);
    const auto f_98 = 0.375 * std::sqrt(35.0);
    const auto f_99 = 0.1875 * std::sqrt(35.0);
    const auto f_100 = 1.5 * std::sqrt(35.0);
    const auto f_101 = 6.0 * std::sqrt(35.0);
    const auto f_102 = 0.09375 * std::sqrt(210.0);
    const auto f_103 = 0.0625 * std::sqrt(210.0);
    const auto f_104 = 0.125 * std::sqrt(210.0);
    const auto f_105 = std::sqrt(210.0);
    const auto f_106 = 0.46875 * std::sqrt(14.0);
    const auto f_107 = 0.46875 * std::sqrt(210.0);
    const auto f_108 = 0.15625 * std::sqrt(210.0);
    const auto f_109 = 0.9375 * std::sqrt(210.0);
    const auto f_110 = 0.3125 * std::sqrt(210.0);
    const auto f_111 = 1.875 * std::sqrt(210.0);
    const auto f_112 = 0.625 * std::sqrt(210.0);
    const auto f_113 = 0.25 * std::sqrt(210.0);
    const auto f_114 = 1.875 * std::sqrt(35.0);
    const auto f_115 = 3.75 * std::sqrt(35.0);
    const auto f_116 = 7.5 * std::sqrt(35.0);
    const auto f_117 = 0.75 * std::sqrt(14.0);
    const auto f_118 = 3.0 * std::sqrt(14.0);
    const auto f_119 = 0.625 * std::sqrt(21.0);
    const auto f_120 = 1.875 * std::sqrt(21.0);
    const auto f_121 = 1.5 * std::sqrt(21.0);
    const auto f_122 = std::sqrt(21.0);
    const auto f_123 = 0.9375 * std::sqrt(35.0);
    const auto f_124 = 0.234375 * std::sqrt(10.0);
    const auto f_125 = 0.078125 * std::sqrt(10.0);
    const auto f_126 = 0.703125 * std::sqrt(10.0);
    const auto f_127 = 4.21875 * std::sqrt(10.0);
    const auto f_128 = 1.40625 * std::sqrt(10.0);
    const auto f_129 = 8.4375 * std::sqrt(10.0);
    const auto f_130 = 2.8125 * std::sqrt(10.0);
    const auto f_131 = 5.625 * std::sqrt(10.0);
    const auto f_132 = 1.875 * std::sqrt(10.0);
    const auto f_133 = 0.75 * std::sqrt(10.0);
    const auto f_134 = 0.25 * std::sqrt(10.0);
    const auto f_135 = 0.3125 * std::sqrt(15.0);
    const auto f_136 = 0.9375 * std::sqrt(15.0);
    const auto f_137 = 5.625 * std::sqrt(15.0);
    const auto f_138 = 11.25 * std::sqrt(15.0);
    const auto f_139 = 7.5 * std::sqrt(15.0);
    const auto f_140 = std::sqrt(15.0);
    const auto f_141 = 0.078125 * std::sqrt(6.0);
    const auto f_142 = 0.3125 * std::sqrt(6.0);
    const auto f_143 = 0.234375 * std::sqrt(6.0);
    const auto f_144 = 0.9375 * std::sqrt(6.0);
    const auto f_145 = 1.40625 * std::sqrt(6.0);
    const auto f_146 = 5.625 * std::sqrt(6.0);
    const auto f_147 = 2.8125 * std::sqrt(6.0);
    const auto f_148 = 11.25 * std::sqrt(6.0);
    const auto f_149 = 1.875 * std::sqrt(6.0);
    const auto f_150 = 7.5 * std::sqrt(6.0);
    const auto f_151 = 0.25 * std::sqrt(6.0);
    const auto f_152 = std::sqrt(6.0);
    const auto f_153 = 0.15625 * std::sqrt(15.0);
    const auto f_154 = 0.46875 * std::sqrt(15.0);
    const auto f_155 = 2.8125 * std::sqrt(15.0);
    const auto f_156 = 3.75 * std::sqrt(15.0);
    const auto f_157 = 0.5 * std::sqrt(15.0);
    const auto f_158 = 0.234375 * std::sqrt(21.0);
    const auto f_159 = 0.078125 * std::sqrt(21.0);
    const auto f_160 = 0.046875 * std::sqrt(35.0);
    const auto f_161 = 0.046875 * std::sqrt(210.0);
    const auto f_162 = 0.03125 * std::sqrt(210.0);
    const auto f_163 = 0.234375 * std::sqrt(14.0);
    const auto f_164 = 0.140625 * std::sqrt(70.0);
    const auto f_165 = 0.046875 * std::sqrt(70.0);
    const auto f_166 = 0.703125 * std::sqrt(70.0);
    const auto f_167 = 0.234375 * std::sqrt(70.0);
    const auto f_168 = 1.40625 * std::sqrt(70.0);
    const auto f_169 = 0.46875 * std::sqrt(70.0);
    const auto f_170 = 8.4375 * std::sqrt(70.0);
    const auto f_171 = 2.8125 * std::sqrt(70.0);
    const auto f_172 = 0.1875 * std::sqrt(105.0);
    const auto f_173 = 0.9375 * std::sqrt(105.0);
    const auto f_174 = 1.875 * std::sqrt(105.0);
    const auto f_175 = 11.25 * std::sqrt(105.0);
    const auto f_176 = 0.046875 * std::sqrt(42.0);
    const auto f_177 = 0.234375 * std::sqrt(42.0);
    const auto f_178 = 0.9375 * std::sqrt(42.0);
    const auto f_179 = 0.46875 * std::sqrt(42.0);
    const auto f_180 = 2.8125 * std::sqrt(42.0);
    const auto f_181 = 11.25 * std::sqrt(42.0);
    const auto f_182 = 0.28125 * std::sqrt(7.0);
    const auto f_183 = 0.1875 * std::sqrt(7.0);
    const auto f_184 = 1.40625 * std::sqrt(7.0);
    const auto f_185 = 0.9375 * std::sqrt(7.0);
    const auto f_186 = 2.8125 * std::sqrt(7.0);
    const auto f_187 = 1.875 * std::sqrt(7.0);
    const auto f_188 = 16.875 * std::sqrt(7.0);
    const auto f_189 = 0.09375 * std::sqrt(105.0);
    const auto f_190 = 0.46875 * std::sqrt(105.0);
    const auto f_191 = 5.625 * std::sqrt(105.0);
    const auto f_192 = 0.046875 * std::sqrt(1155.0);
    const auto f_193 = 0.015625 * std::sqrt(1155.0);
    const auto f_194 = 0.703125 * std::sqrt(1155.0);
    const auto f_195 = 0.234375 * std::sqrt(1155.0);
    const auto f_196 = 0.09375 * std::sqrt(770.0);
    const auto f_197 = 1.40625 * std::sqrt(770.0);
    const auto f_198 = 0.046875 * std::sqrt(77.0);
    const auto f_199 = 0.1875 * std::sqrt(77.0);
    const auto f_200 = 0.703125 * std::sqrt(77.0);
    const auto f_201 = 2.8125 * std::sqrt(77.0);
    const auto f_202 = 0.046875 * std::sqrt(462.0);
    const auto f_203 = 0.03125 * std::sqrt(462.0);
    const auto f_204 = 0.703125 * std::sqrt(462.0);
    const auto f_205 = 0.46875 * std::sqrt(462.0);
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

    const auto *if__0 = buffer.data(if_ + 0);
    const auto *if__1 = buffer.data(if_ + 1);
    const auto *if__2 = buffer.data(if_ + 2);
    const auto *if__3 = buffer.data(if_ + 3);
    const auto *if__4 = buffer.data(if_ + 4);
    const auto *if__5 = buffer.data(if_ + 5);
    const auto *if__6 = buffer.data(if_ + 6);
    const auto *if__7 = buffer.data(if_ + 7);
    const auto *if__8 = buffer.data(if_ + 8);
    const auto *if__9 = buffer.data(if_ + 9);
    const auto *if__10 = buffer.data(if_ + 10);
    const auto *if__11 = buffer.data(if_ + 11);
    const auto *if__12 = buffer.data(if_ + 12);
    const auto *if__13 = buffer.data(if_ + 13);
    const auto *if__14 = buffer.data(if_ + 14);
    const auto *if__15 = buffer.data(if_ + 15);
    const auto *if__16 = buffer.data(if_ + 16);
    const auto *if__17 = buffer.data(if_ + 17);
    const auto *if__18 = buffer.data(if_ + 18);
    const auto *if__19 = buffer.data(if_ + 19);
    const auto *if__20 = buffer.data(if_ + 20);
    const auto *if__21 = buffer.data(if_ + 21);
    const auto *if__22 = buffer.data(if_ + 22);
    const auto *if__23 = buffer.data(if_ + 23);
    const auto *if__24 = buffer.data(if_ + 24);
    const auto *if__25 = buffer.data(if_ + 25);
    const auto *if__26 = buffer.data(if_ + 26);
    const auto *if__27 = buffer.data(if_ + 27);
    const auto *if__28 = buffer.data(if_ + 28);
    const auto *if__29 = buffer.data(if_ + 29);
    const auto *if__30 = buffer.data(if_ + 30);
    const auto *if__31 = buffer.data(if_ + 31);
    const auto *if__32 = buffer.data(if_ + 32);
    const auto *if__33 = buffer.data(if_ + 33);
    const auto *if__34 = buffer.data(if_ + 34);
    const auto *if__35 = buffer.data(if_ + 35);
    const auto *if__36 = buffer.data(if_ + 36);
    const auto *if__37 = buffer.data(if_ + 37);
    const auto *if__38 = buffer.data(if_ + 38);
    const auto *if__39 = buffer.data(if_ + 39);
    const auto *if__40 = buffer.data(if_ + 40);
    const auto *if__41 = buffer.data(if_ + 41);
    const auto *if__42 = buffer.data(if_ + 42);
    const auto *if__43 = buffer.data(if_ + 43);
    const auto *if__44 = buffer.data(if_ + 44);
    const auto *if__45 = buffer.data(if_ + 45);
    const auto *if__46 = buffer.data(if_ + 46);
    const auto *if__47 = buffer.data(if_ + 47);
    const auto *if__48 = buffer.data(if_ + 48);
    const auto *if__49 = buffer.data(if_ + 49);
    const auto *if__50 = buffer.data(if_ + 50);
    const auto *if__51 = buffer.data(if_ + 51);
    const auto *if__52 = buffer.data(if_ + 52);
    const auto *if__53 = buffer.data(if_ + 53);
    const auto *if__54 = buffer.data(if_ + 54);
    const auto *if__55 = buffer.data(if_ + 55);
    const auto *if__56 = buffer.data(if_ + 56);
    const auto *if__57 = buffer.data(if_ + 57);
    const auto *if__58 = buffer.data(if_ + 58);
    const auto *if__59 = buffer.data(if_ + 59);
    const auto *if__60 = buffer.data(if_ + 60);
    const auto *if__61 = buffer.data(if_ + 61);
    const auto *if__62 = buffer.data(if_ + 62);
    const auto *if__63 = buffer.data(if_ + 63);
    const auto *if__64 = buffer.data(if_ + 64);
    const auto *if__65 = buffer.data(if_ + 65);
    const auto *if__66 = buffer.data(if_ + 66);
    const auto *if__67 = buffer.data(if_ + 67);
    const auto *if__68 = buffer.data(if_ + 68);
    const auto *if__69 = buffer.data(if_ + 69);
    const auto *if__70 = buffer.data(if_ + 70);
    const auto *if__71 = buffer.data(if_ + 71);
    const auto *if__72 = buffer.data(if_ + 72);
    const auto *if__73 = buffer.data(if_ + 73);
    const auto *if__74 = buffer.data(if_ + 74);
    const auto *if__75 = buffer.data(if_ + 75);
    const auto *if__76 = buffer.data(if_ + 76);
    const auto *if__77 = buffer.data(if_ + 77);
    const auto *if__78 = buffer.data(if_ + 78);
    const auto *if__79 = buffer.data(if_ + 79);
    const auto *if__80 = buffer.data(if_ + 80);
    const auto *if__81 = buffer.data(if_ + 81);
    const auto *if__82 = buffer.data(if_ + 82);
    const auto *if__83 = buffer.data(if_ + 83);
    const auto *if__84 = buffer.data(if_ + 84);
    const auto *if__85 = buffer.data(if_ + 85);
    const auto *if__86 = buffer.data(if_ + 86);
    const auto *if__87 = buffer.data(if_ + 87);
    const auto *if__88 = buffer.data(if_ + 88);
    const auto *if__89 = buffer.data(if_ + 89);
    const auto *if__90 = buffer.data(if_ + 90);
    const auto *if__91 = buffer.data(if_ + 91);
    const auto *if__92 = buffer.data(if_ + 92);
    const auto *if__93 = buffer.data(if_ + 93);
    const auto *if__94 = buffer.data(if_ + 94);
    const auto *if__95 = buffer.data(if_ + 95);
    const auto *if__96 = buffer.data(if_ + 96);
    const auto *if__97 = buffer.data(if_ + 97);
    const auto *if__98 = buffer.data(if_ + 98);
    const auto *if__99 = buffer.data(if_ + 99);
    const auto *if__100 = buffer.data(if_ + 100);
    const auto *if__101 = buffer.data(if_ + 101);
    const auto *if__102 = buffer.data(if_ + 102);
    const auto *if__103 = buffer.data(if_ + 103);
    const auto *if__104 = buffer.data(if_ + 104);
    const auto *if__105 = buffer.data(if_ + 105);
    const auto *if__106 = buffer.data(if_ + 106);
    const auto *if__107 = buffer.data(if_ + 107);
    const auto *if__108 = buffer.data(if_ + 108);
    const auto *if__109 = buffer.data(if_ + 109);
    const auto *if__110 = buffer.data(if_ + 110);
    const auto *if__111 = buffer.data(if_ + 111);
    const auto *if__112 = buffer.data(if_ + 112);
    const auto *if__113 = buffer.data(if_ + 113);
    const auto *if__114 = buffer.data(if_ + 114);
    const auto *if__115 = buffer.data(if_ + 115);
    const auto *if__116 = buffer.data(if_ + 116);
    const auto *if__117 = buffer.data(if_ + 117);
    const auto *if__118 = buffer.data(if_ + 118);
    const auto *if__119 = buffer.data(if_ + 119);
    const auto *if__120 = buffer.data(if_ + 120);
    const auto *if__121 = buffer.data(if_ + 121);
    const auto *if__122 = buffer.data(if_ + 122);
    const auto *if__123 = buffer.data(if_ + 123);
    const auto *if__124 = buffer.data(if_ + 124);
    const auto *if__125 = buffer.data(if_ + 125);
    const auto *if__126 = buffer.data(if_ + 126);
    const auto *if__127 = buffer.data(if_ + 127);
    const auto *if__128 = buffer.data(if_ + 128);
    const auto *if__129 = buffer.data(if_ + 129);
    const auto *if__130 = buffer.data(if_ + 130);
    const auto *if__131 = buffer.data(if_ + 131);
    const auto *if__132 = buffer.data(if_ + 132);
    const auto *if__133 = buffer.data(if_ + 133);
    const auto *if__134 = buffer.data(if_ + 134);
    const auto *if__135 = buffer.data(if_ + 135);
    const auto *if__136 = buffer.data(if_ + 136);
    const auto *if__137 = buffer.data(if_ + 137);
    const auto *if__138 = buffer.data(if_ + 138);
    const auto *if__139 = buffer.data(if_ + 139);
    const auto *if__140 = buffer.data(if_ + 140);
    const auto *if__141 = buffer.data(if_ + 141);
    const auto *if__142 = buffer.data(if_ + 142);
    const auto *if__143 = buffer.data(if_ + 143);
    const auto *if__144 = buffer.data(if_ + 144);
    const auto *if__145 = buffer.data(if_ + 145);
    const auto *if__146 = buffer.data(if_ + 146);
    const auto *if__147 = buffer.data(if_ + 147);
    const auto *if__148 = buffer.data(if_ + 148);
    const auto *if__149 = buffer.data(if_ + 149);
    const auto *if__150 = buffer.data(if_ + 150);
    const auto *if__151 = buffer.data(if_ + 151);
    const auto *if__152 = buffer.data(if_ + 152);
    const auto *if__153 = buffer.data(if_ + 153);
    const auto *if__154 = buffer.data(if_ + 154);
    const auto *if__155 = buffer.data(if_ + 155);
    const auto *if__156 = buffer.data(if_ + 156);
    const auto *if__157 = buffer.data(if_ + 157);
    const auto *if__158 = buffer.data(if_ + 158);
    const auto *if__159 = buffer.data(if_ + 159);
    const auto *if__160 = buffer.data(if_ + 160);
    const auto *if__161 = buffer.data(if_ + 161);
    const auto *if__162 = buffer.data(if_ + 162);
    const auto *if__163 = buffer.data(if_ + 163);
    const auto *if__164 = buffer.data(if_ + 164);
    const auto *if__165 = buffer.data(if_ + 165);
    const auto *if__166 = buffer.data(if_ + 166);
    const auto *if__167 = buffer.data(if_ + 167);
    const auto *if__168 = buffer.data(if_ + 168);
    const auto *if__169 = buffer.data(if_ + 169);
    const auto *if__170 = buffer.data(if_ + 170);
    const auto *if__171 = buffer.data(if_ + 171);
    const auto *if__172 = buffer.data(if_ + 172);
    const auto *if__173 = buffer.data(if_ + 173);
    const auto *if__174 = buffer.data(if_ + 174);
    const auto *if__175 = buffer.data(if_ + 175);
    const auto *if__176 = buffer.data(if_ + 176);
    const auto *if__177 = buffer.data(if_ + 177);
    const auto *if__178 = buffer.data(if_ + 178);
    const auto *if__179 = buffer.data(if_ + 179);
    const auto *if__180 = buffer.data(if_ + 180);
    const auto *if__181 = buffer.data(if_ + 181);
    const auto *if__182 = buffer.data(if_ + 182);
    const auto *if__183 = buffer.data(if_ + 183);
    const auto *if__184 = buffer.data(if_ + 184);
    const auto *if__185 = buffer.data(if_ + 185);
    const auto *if__186 = buffer.data(if_ + 186);
    const auto *if__187 = buffer.data(if_ + 187);
    const auto *if__188 = buffer.data(if_ + 188);
    const auto *if__189 = buffer.data(if_ + 189);
    const auto *if__190 = buffer.data(if_ + 190);
    const auto *if__191 = buffer.data(if_ + 191);
    const auto *if__192 = buffer.data(if_ + 192);
    const auto *if__193 = buffer.data(if_ + 193);
    const auto *if__194 = buffer.data(if_ + 194);
    const auto *if__195 = buffer.data(if_ + 195);
    const auto *if__196 = buffer.data(if_ + 196);
    const auto *if__197 = buffer.data(if_ + 197);
    const auto *if__198 = buffer.data(if_ + 198);
    const auto *if__199 = buffer.data(if_ + 199);
    const auto *if__200 = buffer.data(if_ + 200);
    const auto *if__201 = buffer.data(if_ + 201);
    const auto *if__202 = buffer.data(if_ + 202);
    const auto *if__203 = buffer.data(if_ + 203);
    const auto *if__204 = buffer.data(if_ + 204);
    const auto *if__205 = buffer.data(if_ + 205);
    const auto *if__206 = buffer.data(if_ + 206);
    const auto *if__207 = buffer.data(if_ + 207);
    const auto *if__208 = buffer.data(if_ + 208);
    const auto *if__209 = buffer.data(if_ + 209);
    const auto *if__210 = buffer.data(if_ + 210);
    const auto *if__211 = buffer.data(if_ + 211);
    const auto *if__212 = buffer.data(if_ + 212);
    const auto *if__213 = buffer.data(if_ + 213);
    const auto *if__214 = buffer.data(if_ + 214);
    const auto *if__215 = buffer.data(if_ + 215);
    const auto *if__216 = buffer.data(if_ + 216);
    const auto *if__217 = buffer.data(if_ + 217);
    const auto *if__218 = buffer.data(if_ + 218);
    const auto *if__219 = buffer.data(if_ + 219);
    const auto *if__220 = buffer.data(if_ + 220);
    const auto *if__221 = buffer.data(if_ + 221);
    const auto *if__222 = buffer.data(if_ + 222);
    const auto *if__223 = buffer.data(if_ + 223);
    const auto *if__224 = buffer.data(if_ + 224);
    const auto *if__225 = buffer.data(if_ + 225);
    const auto *if__226 = buffer.data(if_ + 226);
    const auto *if__227 = buffer.data(if_ + 227);
    const auto *if__228 = buffer.data(if_ + 228);
    const auto *if__229 = buffer.data(if_ + 229);
    const auto *if__230 = buffer.data(if_ + 230);
    const auto *if__231 = buffer.data(if_ + 231);
    const auto *if__232 = buffer.data(if_ + 232);
    const auto *if__233 = buffer.data(if_ + 233);
    const auto *if__234 = buffer.data(if_ + 234);
    const auto *if__235 = buffer.data(if_ + 235);
    const auto *if__236 = buffer.data(if_ + 236);
    const auto *if__237 = buffer.data(if_ + 237);
    const auto *if__238 = buffer.data(if_ + 238);
    const auto *if__239 = buffer.data(if_ + 239);
    const auto *if__240 = buffer.data(if_ + 240);
    const auto *if__241 = buffer.data(if_ + 241);
    const auto *if__242 = buffer.data(if_ + 242);
    const auto *if__243 = buffer.data(if_ + 243);
    const auto *if__244 = buffer.data(if_ + 244);
    const auto *if__245 = buffer.data(if_ + 245);
    const auto *if__246 = buffer.data(if_ + 246);
    const auto *if__247 = buffer.data(if_ + 247);
    const auto *if__248 = buffer.data(if_ + 248);
    const auto *if__249 = buffer.data(if_ + 249);
    const auto *if__250 = buffer.data(if_ + 250);
    const auto *if__251 = buffer.data(if_ + 251);
    const auto *if__252 = buffer.data(if_ + 252);
    const auto *if__253 = buffer.data(if_ + 253);
    const auto *if__254 = buffer.data(if_ + 254);
    const auto *if__255 = buffer.data(if_ + 255);
    const auto *if__256 = buffer.data(if_ + 256);
    const auto *if__257 = buffer.data(if_ + 257);
    const auto *if__258 = buffer.data(if_ + 258);
    const auto *if__259 = buffer.data(if_ + 259);
    const auto *if__260 = buffer.data(if_ + 260);
    const auto *if__261 = buffer.data(if_ + 261);
    const auto *if__262 = buffer.data(if_ + 262);
    const auto *if__263 = buffer.data(if_ + 263);
    const auto *if__264 = buffer.data(if_ + 264);
    const auto *if__265 = buffer.data(if_ + 265);
    const auto *if__266 = buffer.data(if_ + 266);
    const auto *if__267 = buffer.data(if_ + 267);
    const auto *if__268 = buffer.data(if_ + 268);
    const auto *if__269 = buffer.data(if_ + 269);
    const auto *if__270 = buffer.data(if_ + 270);
    const auto *if__271 = buffer.data(if_ + 271);
    const auto *if__272 = buffer.data(if_ + 272);
    const auto *if__273 = buffer.data(if_ + 273);
    const auto *if__274 = buffer.data(if_ + 274);
    const auto *if__275 = buffer.data(if_ + 275);
    const auto *if__276 = buffer.data(if_ + 276);
    const auto *if__277 = buffer.data(if_ + 277);
    const auto *if__278 = buffer.data(if_ + 278);
    const auto *if__279 = buffer.data(if_ + 279);

#pragma omp simd aligned(if__11, if__14, if__16, if__18, if__61, if__64, if__66, if__68, \
                         if__151, if__154, if__156, if__158 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * if__11[k]
                 - f_1 * if__16[k]
                 - f_2 * if__61[k]
                 + f_3 * if__66[k]
                 + f_0 * if__151[k]
                 - f_1 * if__156[k];

        g_1[k] = f_4 * if__14[k]
                 - f_5 * if__64[k]
                 + f_4 * if__154[k];

        g_2[k] = -f_6 * if__11[k]
                 - f_6 * if__16[k]
                 + f_7 * if__18[k]
                 + f_8 * if__61[k]
                 + f_8 * if__66[k]
                 - f_9 * if__68[k]
                 - f_6 * if__151[k]
                 - f_6 * if__156[k]
                 + f_7 * if__158[k];
    }

#pragma omp simd aligned(if__12, if__17, if__19, if__62, if__67, if__69, if__152, if__157, \
                         if__159 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = -f_10 * if__12[k]
                 - f_10 * if__17[k]
                 + f_11 * if__19[k]
                 + f_12 * if__62[k]
                 + f_12 * if__67[k]
                 - f_13 * if__69[k]
                 - f_10 * if__152[k]
                 - f_10 * if__157[k]
                 + f_11 * if__159[k];
    }

#pragma omp simd aligned(if__10, if__13, if__15, if__60, if__63, if__65, if__150, if__153, \
                         if__155 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = -f_6 * if__10[k]
                 - f_6 * if__13[k]
                 + f_7 * if__15[k]
                 + f_8 * if__60[k]
                 + f_8 * if__63[k]
                 - f_9 * if__65[k]
                 - f_6 * if__150[k]
                 - f_6 * if__153[k]
                 + f_7 * if__155[k];
    }

#pragma omp simd aligned(if__10, if__12, if__13, if__17, if__60, if__62, if__63, if__67, \
                         if__150, if__152, if__153, if__157 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = f_14 * if__12[k]
                 - f_14 * if__17[k]
                 - f_15 * if__62[k]
                 + f_15 * if__67[k]
                 + f_14 * if__152[k]
                 - f_14 * if__157[k];

        g_6[k] = f_1 * if__10[k]
                 - f_0 * if__13[k]
                 - f_3 * if__60[k]
                 + f_2 * if__63[k]
                 + f_1 * if__150[k]
                 - f_0 * if__153[k];
    }

#pragma omp simd aligned(if__41, if__44, if__46, if__48, if__111, if__114, if__116, if__118, \
                         if__221, if__224, if__226, if__228 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = f_16 * if__41[k]
                 - f_17 * if__46[k]
                 - f_18 * if__111[k]
                 + f_19 * if__116[k]
                 + f_20 * if__221[k]
                 - f_21 * if__226[k];

        g_8[k] = f_22 * if__44[k]
                 - f_23 * if__114[k]
                 + f_24 * if__224[k];

        g_9[k] = -f_25 * if__41[k]
                 - f_25 * if__46[k]
                 + f_26 * if__48[k]
                 + f_27 * if__111[k]
                 + f_27 * if__116[k]
                 - f_28 * if__118[k]
                 - f_29 * if__221[k]
                 - f_29 * if__226[k]
                 + f_30 * if__228[k];
    }

#pragma omp simd aligned(if__42, if__47, if__49, if__112, if__117, if__119, if__222, if__227, \
                         if__229 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = -f_31 * if__42[k]
                  - f_31 * if__47[k]
                  + f_32 * if__49[k]
                  + f_33 * if__112[k]
                  + f_33 * if__117[k]
                  - f_34 * if__119[k]
                  - f_35 * if__222[k]
                  - f_35 * if__227[k]
                  + f_36 * if__229[k];
    }

#pragma omp simd aligned(if__40, if__43, if__45, if__110, if__113, if__115, if__220, if__223, \
                         if__225 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_11[k] = -f_25 * if__40[k]
                  - f_25 * if__43[k]
                  + f_26 * if__45[k]
                  + f_27 * if__110[k]
                  + f_27 * if__113[k]
                  - f_28 * if__115[k]
                  - f_29 * if__220[k]
                  - f_29 * if__223[k]
                  + f_30 * if__225[k];
    }

#pragma omp simd aligned(if__40, if__42, if__43, if__47, if__110, if__112, if__113, if__117, \
                         if__220, if__222, if__223, if__227 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_12[k] = f_37 * if__42[k]
                  - f_37 * if__47[k]
                  - f_22 * if__112[k]
                  + f_22 * if__117[k]
                  + f_38 * if__222[k]
                  - f_38 * if__227[k];

        g_13[k] = f_17 * if__40[k]
                  - f_16 * if__43[k]
                  - f_19 * if__110[k]
                  + f_18 * if__113[k]
                  + f_21 * if__220[k]
                  - f_20 * if__223[k];
    }

#pragma omp simd aligned(if__11, if__14, if__16, if__81, if__84, if__86, if__151, if__154, \
                         if__156, if__171, if__174, if__176 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = -f_39 * if__11[k]
                  + f_40 * if__16[k]
                  + f_41 * if__81[k]
                  - f_42 * if__86[k]
                  + f_39 * if__151[k]
                  - f_40 * if__156[k]
                  - f_41 * if__171[k]
                  + f_42 * if__176[k];

        g_15[k] = -f_43 * if__14[k]
                  + f_44 * if__84[k]
                  + f_43 * if__154[k]
                  - f_44 * if__174[k];
    }

#pragma omp simd aligned(if__11, if__16, if__18, if__81, if__86, if__88, if__151, if__156, \
                         if__158, if__171, if__176, if__178 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_16[k] = f_45 * if__11[k]
                  + f_45 * if__16[k]
                  - f_46 * if__18[k]
                  - f_47 * if__81[k]
                  - f_47 * if__86[k]
                  + f_48 * if__88[k]
                  - f_45 * if__151[k]
                  - f_45 * if__156[k]
                  + f_46 * if__158[k]
                  + f_47 * if__171[k]
                  + f_47 * if__176[k]
                  - f_48 * if__178[k];
    }

#pragma omp simd aligned(if__12, if__17, if__19, if__82, if__87, if__89, if__152, if__157, \
                         if__159, if__172, if__177, if__179 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = f_49 * if__12[k]
                  + f_49 * if__17[k]
                  - f_50 * if__19[k]
                  - f_51 * if__82[k]
                  - f_51 * if__87[k]
                  + f_52 * if__89[k]
                  - f_49 * if__152[k]
                  - f_49 * if__157[k]
                  + f_50 * if__159[k]
                  + f_51 * if__172[k]
                  + f_51 * if__177[k]
                  - f_52 * if__179[k];
    }

#pragma omp simd aligned(if__10, if__13, if__15, if__80, if__83, if__85, if__150, if__153, \
                         if__155, if__170, if__173, if__175 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = f_45 * if__10[k]
                  + f_45 * if__13[k]
                  - f_46 * if__15[k]
                  - f_47 * if__80[k]
                  - f_47 * if__83[k]
                  + f_48 * if__85[k]
                  - f_45 * if__150[k]
                  - f_45 * if__153[k]
                  + f_46 * if__155[k]
                  + f_47 * if__170[k]
                  + f_47 * if__173[k]
                  - f_48 * if__175[k];
    }

#pragma omp simd aligned(if__12, if__17, if__82, if__87, if__152, if__157, if__172, \
                         if__177 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = -f_53 * if__12[k]
                  + f_53 * if__17[k]
                  + f_54 * if__82[k]
                  - f_54 * if__87[k]
                  + f_53 * if__152[k]
                  - f_53 * if__157[k]
                  - f_54 * if__172[k]
                  + f_54 * if__177[k];
    }

#pragma omp simd aligned(if__10, if__13, if__80, if__83, if__150, if__153, if__170, \
                         if__173 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = -f_40 * if__10[k]
                  + f_39 * if__13[k]
                  + f_42 * if__80[k]
                  - f_41 * if__83[k]
                  + f_40 * if__150[k]
                  - f_39 * if__153[k]
                  - f_42 * if__170[k]
                  + f_41 * if__173[k];
    }

#pragma omp simd aligned(if__41, if__46, if__111, if__116, if__131, if__136, if__221, if__226, \
                         if__241, if__246 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = -f_55 * if__41[k]
                  + f_56 * if__46[k]
                  - f_57 * if__111[k]
                  + f_58 * if__116[k]
                  + f_59 * if__131[k]
                  - f_60 * if__136[k]
                  + f_56 * if__221[k]
                  - f_61 * if__226[k]
                  - f_60 * if__241[k]
                  + f_62 * if__246[k];
    }

#pragma omp simd aligned(if__44, if__114, if__134, if__224, if__244 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = -f_63 * if__44[k]
                  - f_64 * if__114[k]
                  + f_65 * if__134[k]
                  + f_66 * if__224[k]
                  - f_67 * if__244[k];
    }

#pragma omp simd aligned(if__41, if__46, if__48, if__111, if__116, if__118, if__131, if__136, \
                         if__138, if__221, if__226, if__228, if__241, if__246, \
                         if__248 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = f_68 * if__41[k]
                  + f_68 * if__46[k]
                  - f_69 * if__48[k]
                  + f_70 * if__111[k]
                  + f_70 * if__116[k]
                  - f_71 * if__118[k]
                  - f_71 * if__131[k]
                  - f_71 * if__136[k]
                  + f_72 * if__138[k]
                  - f_73 * if__221[k]
                  - f_73 * if__226[k]
                  + f_74 * if__228[k]
                  + f_75 * if__241[k]
                  + f_75 * if__246[k]
                  - f_76 * if__248[k];
    }

#pragma omp simd aligned(if__42, if__47, if__49, if__112, if__117, if__119, if__132, if__137, \
                         if__139, if__222, if__227, if__229, if__242, if__247, \
                         if__249 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = f_77 * if__42[k]
                  + f_77 * if__47[k]
                  - f_78 * if__49[k]
                  + f_78 * if__112[k]
                  + f_78 * if__117[k]
                  - f_79 * if__119[k]
                  - f_80 * if__132[k]
                  - f_80 * if__137[k]
                  + f_81 * if__139[k]
                  - f_82 * if__222[k]
                  - f_82 * if__227[k]
                  + f_83 * if__229[k]
                  + f_84 * if__242[k]
                  + f_84 * if__247[k]
                  - f_85 * if__249[k];
    }

#pragma omp simd aligned(if__40, if__43, if__45, if__110, if__113, if__115, if__130, if__133, \
                         if__135, if__220, if__223, if__225, if__240, if__243, \
                         if__245 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = f_68 * if__40[k]
                  + f_68 * if__43[k]
                  - f_69 * if__45[k]
                  + f_70 * if__110[k]
                  + f_70 * if__113[k]
                  - f_71 * if__115[k]
                  - f_71 * if__130[k]
                  - f_71 * if__133[k]
                  + f_72 * if__135[k]
                  - f_73 * if__220[k]
                  - f_73 * if__223[k]
                  + f_74 * if__225[k]
                  + f_75 * if__240[k]
                  + f_75 * if__243[k]
                  - f_76 * if__245[k];
    }

#pragma omp simd aligned(if__42, if__47, if__112, if__117, if__132, if__137, if__222, if__227, \
                         if__242, if__247 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = -f_86 * if__42[k]
                  + f_86 * if__47[k]
                  - f_66 * if__112[k]
                  + f_66 * if__117[k]
                  + f_87 * if__132[k]
                  - f_87 * if__137[k]
                  + f_88 * if__222[k]
                  - f_88 * if__227[k]
                  - f_89 * if__242[k]
                  + f_89 * if__247[k];
    }

#pragma omp simd aligned(if__40, if__43, if__110, if__113, if__130, if__133, if__220, if__223, \
                         if__240, if__243 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = -f_56 * if__40[k]
                  + f_55 * if__43[k]
                  - f_58 * if__110[k]
                  + f_57 * if__113[k]
                  + f_60 * if__130[k]
                  - f_59 * if__133[k]
                  + f_61 * if__220[k]
                  - f_56 * if__223[k]
                  - f_62 * if__240[k]
                  + f_60 * if__243[k];
    }

#pragma omp simd aligned(if__11, if__16, if__61, if__66, if__81, if__86, if__151, if__156, \
                         if__171, if__176, if__191, if__196 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = f_61 * if__11[k]
                  - f_90 * if__16[k]
                  + f_58 * if__61[k]
                  - f_91 * if__66[k]
                  - f_92 * if__81[k]
                  + f_93 * if__86[k]
                  + f_61 * if__151[k]
                  - f_90 * if__156[k]
                  - f_92 * if__171[k]
                  + f_93 * if__176[k]
                  + f_92 * if__191[k]
                  - f_93 * if__196[k];
    }

#pragma omp simd aligned(if__14, if__64, if__84, if__154, if__174, \
                         if__194 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = f_94 * if__14[k]
                  + f_95 * if__64[k]
                  - f_96 * if__84[k]
                  + f_94 * if__154[k]
                  - f_96 * if__174[k]
                  + f_96 * if__194[k];
    }

#pragma omp simd aligned(if__11, if__16, if__18, if__61, if__66, if__68, if__81, if__86, \
                         if__88, if__151, if__156, if__158, if__171, if__176, if__178, \
                         if__191, if__196, if__198 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = -f_97 * if__11[k]
                  - f_97 * if__16[k]
                  + f_98 * if__18[k]
                  - f_99 * if__61[k]
                  - f_99 * if__66[k]
                  + f_75 * if__68[k]
                  + f_100 * if__81[k]
                  + f_100 * if__86[k]
                  - f_101 * if__88[k]
                  - f_97 * if__151[k]
                  - f_97 * if__156[k]
                  + f_98 * if__158[k]
                  + f_100 * if__171[k]
                  + f_100 * if__176[k]
                  - f_101 * if__178[k]
                  - f_100 * if__191[k]
                  - f_100 * if__196[k]
                  + f_101 * if__198[k];
    }

#pragma omp simd aligned(if__12, if__17, if__19, if__62, if__67, if__69, if__82, if__87, \
                         if__89, if__152, if__157, if__159, if__172, if__177, if__179, \
                         if__192, if__197, if__199 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = -f_102 * if__12[k]
                  - f_102 * if__17[k]
                  + f_103 * if__19[k]
                  - f_83 * if__62[k]
                  - f_83 * if__67[k]
                  + f_104 * if__69[k]
                  + f_81 * if__82[k]
                  + f_81 * if__87[k]
                  - f_105 * if__89[k]
                  - f_102 * if__152[k]
                  - f_102 * if__157[k]
                  + f_103 * if__159[k]
                  + f_81 * if__172[k]
                  + f_81 * if__177[k]
                  - f_105 * if__179[k]
                  - f_81 * if__192[k]
                  - f_81 * if__197[k]
                  + f_105 * if__199[k];
    }

#pragma omp simd aligned(if__10, if__13, if__15, if__60, if__63, if__65, if__80, if__83, \
                         if__85, if__150, if__153, if__155, if__170, if__173, if__175, \
                         if__190, if__193, if__195 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = -f_97 * if__10[k]
                  - f_97 * if__13[k]
                  + f_98 * if__15[k]
                  - f_99 * if__60[k]
                  - f_99 * if__63[k]
                  + f_75 * if__65[k]
                  + f_100 * if__80[k]
                  + f_100 * if__83[k]
                  - f_101 * if__85[k]
                  - f_97 * if__150[k]
                  - f_97 * if__153[k]
                  + f_98 * if__155[k]
                  + f_100 * if__170[k]
                  + f_100 * if__173[k]
                  - f_101 * if__175[k]
                  - f_100 * if__190[k]
                  - f_100 * if__193[k]
                  + f_101 * if__195[k];
    }

#pragma omp simd aligned(if__12, if__17, if__62, if__67, if__82, if__87, if__152, if__157, \
                         if__172, if__177, if__192, if__197 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = f_106 * if__12[k]
                  - f_106 * if__17[k]
                  + f_94 * if__62[k]
                  - f_94 * if__67[k]
                  - f_67 * if__82[k]
                  + f_67 * if__87[k]
                  + f_106 * if__152[k]
                  - f_106 * if__157[k]
                  - f_67 * if__172[k]
                  + f_67 * if__177[k]
                  + f_67 * if__192[k]
                  - f_67 * if__197[k];
    }

#pragma omp simd aligned(if__10, if__13, if__60, if__63, if__80, if__83, if__150, if__153, \
                         if__170, if__173, if__190, if__193 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = f_90 * if__10[k]
                  - f_61 * if__13[k]
                  + f_91 * if__60[k]
                  - f_58 * if__63[k]
                  - f_93 * if__80[k]
                  + f_92 * if__83[k]
                  + f_90 * if__150[k]
                  - f_61 * if__153[k]
                  - f_93 * if__170[k]
                  + f_92 * if__173[k]
                  + f_93 * if__190[k]
                  - f_92 * if__193[k];
    }

#pragma omp simd aligned(if__41, if__46, if__111, if__116, if__131, if__136, if__221, if__226, \
                         if__241, if__246, if__261, if__266 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = f_107 * if__41[k]
                  - f_108 * if__46[k]
                  + f_109 * if__111[k]
                  - f_110 * if__116[k]
                  - f_111 * if__131[k]
                  + f_112 * if__136[k]
                  + f_107 * if__221[k]
                  - f_108 * if__226[k]
                  - f_111 * if__241[k]
                  + f_112 * if__246[k]
                  + f_84 * if__261[k]
                  - f_113 * if__266[k];
    }

#pragma omp simd aligned(if__44, if__114, if__134, if__224, if__244, \
                         if__264 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_36[k] = f_114 * if__44[k]
                  + f_115 * if__114[k]
                  - f_116 * if__134[k]
                  + f_114 * if__224[k]
                  - f_116 * if__244[k]
                  + f_76 * if__264[k];
    }

#pragma omp simd aligned(if__41, if__46, if__48, if__111, if__116, if__118, if__131, if__136, \
                         if__138, if__221, if__226, if__228, if__241, if__246, if__248, \
                         if__261, if__266, if__268 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_37[k] = -f_106 * if__41[k]
                  - f_106 * if__46[k]
                  + f_95 * if__48[k]
                  - f_94 * if__111[k]
                  - f_94 * if__116[k]
                  + f_89 * if__118[k]
                  + f_95 * if__131[k]
                  + f_95 * if__136[k]
                  - f_67 * if__138[k]
                  - f_106 * if__221[k]
                  - f_106 * if__226[k]
                  + f_95 * if__228[k]
                  + f_95 * if__241[k]
                  + f_95 * if__246[k]
                  - f_67 * if__248[k]
                  - f_117 * if__261[k]
                  - f_117 * if__266[k]
                  + f_118 * if__268[k];
    }

#pragma omp simd aligned(if__42, if__47, if__49, if__112, if__117, if__119, if__132, if__137, \
                         if__139, if__222, if__227, if__229, if__242, if__247, if__249, \
                         if__262, if__267, if__269 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = -f_58 * if__42[k]
                  - f_58 * if__47[k]
                  + f_119 * if__49[k]
                  - f_120 * if__112[k]
                  - f_120 * if__117[k]
                  + f_62 * if__119[k]
                  + f_60 * if__132[k]
                  + f_60 * if__137[k]
                  - f_93 * if__139[k]
                  - f_58 * if__222[k]
                  - f_58 * if__227[k]
                  + f_119 * if__229[k]
                  + f_60 * if__242[k]
                  + f_60 * if__247[k]
                  - f_93 * if__249[k]
                  - f_121 * if__262[k]
                  - f_121 * if__267[k]
                  + f_122 * if__269[k];
    }

#pragma omp simd aligned(if__40, if__43, if__45, if__110, if__113, if__115, if__130, if__133, \
                         if__135, if__220, if__223, if__225, if__240, if__243, if__245, \
                         if__260, if__263, if__265 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_39[k] = -f_106 * if__40[k]
                  - f_106 * if__43[k]
                  + f_95 * if__45[k]
                  - f_94 * if__110[k]
                  - f_94 * if__113[k]
                  + f_89 * if__115[k]
                  + f_95 * if__130[k]
                  + f_95 * if__133[k]
                  - f_67 * if__135[k]
                  - f_106 * if__220[k]
                  - f_106 * if__223[k]
                  + f_95 * if__225[k]
                  + f_95 * if__240[k]
                  + f_95 * if__243[k]
                  - f_67 * if__245[k]
                  - f_117 * if__260[k]
                  - f_117 * if__263[k]
                  + f_118 * if__265[k];
    }

#pragma omp simd aligned(if__42, if__47, if__112, if__117, if__132, if__137, if__222, if__227, \
                         if__242, if__247, if__262, if__267 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_40[k] = f_123 * if__42[k]
                  - f_123 * if__47[k]
                  + f_114 * if__112[k]
                  - f_114 * if__117[k]
                  - f_115 * if__132[k]
                  + f_115 * if__137[k]
                  + f_123 * if__222[k]
                  - f_123 * if__227[k]
                  - f_115 * if__242[k]
                  + f_115 * if__247[k]
                  + f_100 * if__262[k]
                  - f_100 * if__267[k];
    }

#pragma omp simd aligned(if__40, if__43, if__110, if__113, if__130, if__133, if__220, if__223, \
                         if__240, if__243, if__260, if__263 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_41[k] = f_108 * if__40[k]
                  - f_107 * if__43[k]
                  + f_110 * if__110[k]
                  - f_109 * if__113[k]
                  - f_112 * if__130[k]
                  + f_111 * if__133[k]
                  + f_108 * if__220[k]
                  - f_107 * if__223[k]
                  - f_112 * if__240[k]
                  + f_111 * if__243[k]
                  + f_113 * if__260[k]
                  - f_84 * if__263[k];
    }

#pragma omp simd aligned(if__1, if__6, if__31, if__36, if__51, if__56, if__101, if__106, \
                         if__121, if__126, if__141, if__146, if__211, if__216, if__231, \
                         if__236, if__251, if__256, if__271, if__276 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = -f_124 * if__1[k]
                  + f_125 * if__6[k]
                  - f_126 * if__31[k]
                  + f_124 * if__36[k]
                  + f_127 * if__51[k]
                  - f_128 * if__56[k]
                  - f_126 * if__101[k]
                  + f_124 * if__106[k]
                  + f_129 * if__121[k]
                  - f_130 * if__126[k]
                  - f_131 * if__141[k]
                  + f_132 * if__146[k]
                  - f_124 * if__211[k]
                  + f_125 * if__216[k]
                  + f_127 * if__231[k]
                  - f_128 * if__236[k]
                  - f_131 * if__251[k]
                  + f_132 * if__256[k]
                  + f_133 * if__271[k]
                  - f_134 * if__276[k];
    }

#pragma omp simd aligned(if__4, if__34, if__54, if__104, if__124, if__144, if__214, if__234, \
                         if__254, if__274 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_43[k] = -f_135 * if__4[k]
                  - f_136 * if__34[k]
                  + f_137 * if__54[k]
                  - f_136 * if__104[k]
                  + f_138 * if__124[k]
                  - f_139 * if__144[k]
                  - f_135 * if__214[k]
                  + f_137 * if__234[k]
                  - f_139 * if__254[k]
                  + f_140 * if__274[k];
    }

#pragma omp simd aligned(if__1, if__6, if__8, if__31, if__36, if__38, if__51, if__56, if__58, \
                         if__101, if__106, if__108, if__121, if__126, if__128, if__141, \
                         if__146, if__148, if__211, if__216, if__218, if__231, if__236, \
                         if__238, if__251, if__256, if__258, if__271, if__276, \
                         if__278 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_44[k] = f_141 * if__1[k]
                  + f_141 * if__6[k]
                  - f_142 * if__8[k]
                  + f_143 * if__31[k]
                  + f_143 * if__36[k]
                  - f_144 * if__38[k]
                  - f_145 * if__51[k]
                  - f_145 * if__56[k]
                  + f_146 * if__58[k]
                  + f_143 * if__101[k]
                  + f_143 * if__106[k]
                  - f_144 * if__108[k]
                  - f_147 * if__121[k]
                  - f_147 * if__126[k]
                  + f_148 * if__128[k]
                  + f_149 * if__141[k]
                  + f_149 * if__146[k]
                  - f_150 * if__148[k]
                  + f_141 * if__211[k]
                  + f_141 * if__216[k]
                  - f_142 * if__218[k]
                  - f_145 * if__231[k]
                  - f_145 * if__236[k]
                  + f_146 * if__238[k]
                  + f_149 * if__251[k]
                  + f_149 * if__256[k]
                  - f_150 * if__258[k]
                  - f_151 * if__271[k]
                  - f_151 * if__276[k]
                  + f_152 * if__278[k];
    }

#pragma omp simd aligned(if__2, if__7, if__9, if__32, if__37, if__39, if__52, if__57, if__59, \
                         if__102, if__107, if__109, if__122, if__127, if__129, if__142, \
                         if__147, if__149, if__212, if__217, if__219, if__232, if__237, \
                         if__239, if__252, if__257, if__259, if__272, if__277, \
                         if__279 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_45[k] = 0.46875 * if__2[k]
                  + 0.46875 * if__7[k]
                  - 0.3125 * if__9[k]
                  + 1.40625 * if__32[k]
                  + 1.40625 * if__37[k]
                  - 0.9375 * if__39[k]
                  - 8.4375 * if__52[k]
                  - 8.4375 * if__57[k]
                  + 5.625 * if__59[k]
                  + 1.40625 * if__102[k]
                  + 1.40625 * if__107[k]
                  - 0.9375 * if__109[k]
                  - 16.875 * if__122[k]
                  - 16.875 * if__127[k]
                  + 11.25 * if__129[k]
                  + 11.25 * if__142[k]
                  + 11.25 * if__147[k]
                  - 7.5 * if__149[k]
                  + 0.46875 * if__212[k]
                  + 0.46875 * if__217[k]
                  - 0.3125 * if__219[k]
                  - 8.4375 * if__232[k]
                  - 8.4375 * if__237[k]
                  + 5.625 * if__239[k]
                  + 11.25 * if__252[k]
                  + 11.25 * if__257[k]
                  - 7.5 * if__259[k]
                  - 1.5 * if__272[k]
                  - 1.5 * if__277[k]
                  + if__279[k];
    }

#pragma omp simd aligned(if__0, if__3, if__5, if__30, if__33, if__35, if__50, if__53, if__55, \
                         if__100, if__103, if__105, if__120, if__123, if__125, if__140, \
                         if__143, if__145, if__210, if__213, if__215, if__230, if__233, \
                         if__235, if__250, if__253, if__255, if__270, if__273, \
                         if__275 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_46[k] = f_141 * if__0[k]
                  + f_141 * if__3[k]
                  - f_142 * if__5[k]
                  + f_143 * if__30[k]
                  + f_143 * if__33[k]
                  - f_144 * if__35[k]
                  - f_145 * if__50[k]
                  - f_145 * if__53[k]
                  + f_146 * if__55[k]
                  + f_143 * if__100[k]
                  + f_143 * if__103[k]
                  - f_144 * if__105[k]
                  - f_147 * if__120[k]
                  - f_147 * if__123[k]
                  + f_148 * if__125[k]
                  + f_149 * if__140[k]
                  + f_149 * if__143[k]
                  - f_150 * if__145[k]
                  + f_141 * if__210[k]
                  + f_141 * if__213[k]
                  - f_142 * if__215[k]
                  - f_145 * if__230[k]
                  - f_145 * if__233[k]
                  + f_146 * if__235[k]
                  + f_149 * if__250[k]
                  + f_149 * if__253[k]
                  - f_150 * if__255[k]
                  - f_151 * if__270[k]
                  - f_151 * if__273[k]
                  + f_152 * if__275[k];
    }

#pragma omp simd aligned(if__2, if__7, if__32, if__37, if__52, if__57, if__102, if__107, \
                         if__122, if__127, if__142, if__147, if__212, if__217, if__232, \
                         if__237, if__252, if__257, if__272, if__277 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_47[k] = -f_153 * if__2[k]
                  + f_153 * if__7[k]
                  - f_154 * if__32[k]
                  + f_154 * if__37[k]
                  + f_155 * if__52[k]
                  - f_155 * if__57[k]
                  - f_154 * if__102[k]
                  + f_154 * if__107[k]
                  + f_137 * if__122[k]
                  - f_137 * if__127[k]
                  - f_156 * if__142[k]
                  + f_156 * if__147[k]
                  - f_153 * if__212[k]
                  + f_153 * if__217[k]
                  + f_155 * if__232[k]
                  - f_155 * if__237[k]
                  - f_156 * if__252[k]
                  + f_156 * if__257[k]
                  + f_157 * if__272[k]
                  - f_157 * if__277[k];
    }

#pragma omp simd aligned(if__0, if__3, if__30, if__33, if__50, if__53, if__100, if__103, \
                         if__120, if__123, if__140, if__143, if__210, if__213, if__230, \
                         if__233, if__250, if__253, if__270, if__273 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_48[k] = -f_125 * if__0[k]
                  + f_124 * if__3[k]
                  - f_124 * if__30[k]
                  + f_126 * if__33[k]
                  + f_128 * if__50[k]
                  - f_127 * if__53[k]
                  - f_124 * if__100[k]
                  + f_126 * if__103[k]
                  + f_130 * if__120[k]
                  - f_129 * if__123[k]
                  - f_132 * if__140[k]
                  + f_131 * if__143[k]
                  - f_125 * if__210[k]
                  + f_124 * if__213[k]
                  + f_128 * if__230[k]
                  - f_127 * if__233[k]
                  - f_132 * if__250[k]
                  + f_131 * if__253[k]
                  + f_134 * if__270[k]
                  - f_133 * if__273[k];
    }

#pragma omp simd aligned(if__21, if__26, if__71, if__76, if__91, if__96, if__161, if__166, \
                         if__181, if__186, if__201, if__206 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_49[k] = f_107 * if__21[k]
                  - f_108 * if__26[k]
                  + f_109 * if__71[k]
                  - f_110 * if__76[k]
                  - f_111 * if__91[k]
                  + f_112 * if__96[k]
                  + f_107 * if__161[k]
                  - f_108 * if__166[k]
                  - f_111 * if__181[k]
                  + f_112 * if__186[k]
                  + f_84 * if__201[k]
                  - f_113 * if__206[k];
    }

#pragma omp simd aligned(if__24, if__74, if__94, if__164, if__184, \
                         if__204 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_50[k] = f_114 * if__24[k]
                  + f_115 * if__74[k]
                  - f_116 * if__94[k]
                  + f_114 * if__164[k]
                  - f_116 * if__184[k]
                  + f_76 * if__204[k];
    }

#pragma omp simd aligned(if__21, if__26, if__28, if__71, if__76, if__78, if__91, if__96, \
                         if__98, if__161, if__166, if__168, if__181, if__186, if__188, \
                         if__201, if__206, if__208 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_51[k] = -f_106 * if__21[k]
                  - f_106 * if__26[k]
                  + f_95 * if__28[k]
                  - f_94 * if__71[k]
                  - f_94 * if__76[k]
                  + f_89 * if__78[k]
                  + f_95 * if__91[k]
                  + f_95 * if__96[k]
                  - f_67 * if__98[k]
                  - f_106 * if__161[k]
                  - f_106 * if__166[k]
                  + f_95 * if__168[k]
                  + f_95 * if__181[k]
                  + f_95 * if__186[k]
                  - f_67 * if__188[k]
                  - f_117 * if__201[k]
                  - f_117 * if__206[k]
                  + f_118 * if__208[k];
    }

#pragma omp simd aligned(if__22, if__27, if__29, if__72, if__77, if__79, if__92, if__97, \
                         if__99, if__162, if__167, if__169, if__182, if__187, if__189, \
                         if__202, if__207, if__209 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_52[k] = -f_58 * if__22[k]
                  - f_58 * if__27[k]
                  + f_119 * if__29[k]
                  - f_120 * if__72[k]
                  - f_120 * if__77[k]
                  + f_62 * if__79[k]
                  + f_60 * if__92[k]
                  + f_60 * if__97[k]
                  - f_93 * if__99[k]
                  - f_58 * if__162[k]
                  - f_58 * if__167[k]
                  + f_119 * if__169[k]
                  + f_60 * if__182[k]
                  + f_60 * if__187[k]
                  - f_93 * if__189[k]
                  - f_121 * if__202[k]
                  - f_121 * if__207[k]
                  + f_122 * if__209[k];
    }

#pragma omp simd aligned(if__20, if__23, if__25, if__70, if__73, if__75, if__90, if__93, \
                         if__95, if__160, if__163, if__165, if__180, if__183, if__185, \
                         if__200, if__203, if__205 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_53[k] = -f_106 * if__20[k]
                  - f_106 * if__23[k]
                  + f_95 * if__25[k]
                  - f_94 * if__70[k]
                  - f_94 * if__73[k]
                  + f_89 * if__75[k]
                  + f_95 * if__90[k]
                  + f_95 * if__93[k]
                  - f_67 * if__95[k]
                  - f_106 * if__160[k]
                  - f_106 * if__163[k]
                  + f_95 * if__165[k]
                  + f_95 * if__180[k]
                  + f_95 * if__183[k]
                  - f_67 * if__185[k]
                  - f_117 * if__200[k]
                  - f_117 * if__203[k]
                  + f_118 * if__205[k];
    }

#pragma omp simd aligned(if__22, if__27, if__72, if__77, if__92, if__97, if__162, if__167, \
                         if__182, if__187, if__202, if__207 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_54[k] = f_123 * if__22[k]
                  - f_123 * if__27[k]
                  + f_114 * if__72[k]
                  - f_114 * if__77[k]
                  - f_115 * if__92[k]
                  + f_115 * if__97[k]
                  + f_123 * if__162[k]
                  - f_123 * if__167[k]
                  - f_115 * if__182[k]
                  + f_115 * if__187[k]
                  + f_100 * if__202[k]
                  - f_100 * if__207[k];
    }

#pragma omp simd aligned(if__20, if__23, if__70, if__73, if__90, if__93, if__160, if__163, \
                         if__180, if__183, if__200, if__203 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_55[k] = f_108 * if__20[k]
                  - f_107 * if__23[k]
                  + f_110 * if__70[k]
                  - f_109 * if__73[k]
                  - f_112 * if__90[k]
                  + f_111 * if__93[k]
                  + f_108 * if__160[k]
                  - f_107 * if__163[k]
                  - f_112 * if__180[k]
                  + f_111 * if__183[k]
                  + f_113 * if__200[k]
                  - f_84 * if__203[k];
    }

#pragma omp simd aligned(if__1, if__6, if__31, if__36, if__51, if__56, if__101, if__106, \
                         if__141, if__146, if__211, if__216, if__231, if__236, if__251, \
                         if__256 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_56[k] = f_158 * if__1[k]
                  - f_159 * if__6[k]
                  + f_158 * if__31[k]
                  - f_159 * if__36[k]
                  - f_60 * if__51[k]
                  + f_62 * if__56[k]
                  - f_158 * if__101[k]
                  + f_159 * if__106[k]
                  + f_60 * if__141[k]
                  - f_62 * if__146[k]
                  - f_158 * if__211[k]
                  + f_159 * if__216[k]
                  + f_60 * if__231[k]
                  - f_62 * if__236[k]
                  - f_60 * if__251[k]
                  + f_62 * if__256[k];
    }

#pragma omp simd aligned(if__4, if__34, if__54, if__104, if__144, if__214, if__234, \
                         if__254 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_57[k] = f_106 * if__4[k]
                  + f_106 * if__34[k]
                  - f_67 * if__54[k]
                  - f_106 * if__104[k]
                  + f_67 * if__144[k]
                  - f_106 * if__214[k]
                  + f_67 * if__234[k]
                  - f_67 * if__254[k];
    }

#pragma omp simd aligned(if__1, if__6, if__8, if__31, if__36, if__38, if__51, if__56, if__58, \
                         if__101, if__106, if__108, if__141, if__146, if__148, if__211, \
                         if__216, if__218, if__231, if__236, if__238, if__251, if__256, \
                         if__258 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_58[k] = -f_160 * if__1[k]
                  - f_160 * if__6[k]
                  + f_99 * if__8[k]
                  - f_160 * if__31[k]
                  - f_160 * if__36[k]
                  + f_99 * if__38[k]
                  + f_75 * if__51[k]
                  + f_75 * if__56[k]
                  - f_76 * if__58[k]
                  + f_160 * if__101[k]
                  + f_160 * if__106[k]
                  - f_99 * if__108[k]
                  - f_75 * if__141[k]
                  - f_75 * if__146[k]
                  + f_76 * if__148[k]
                  + f_160 * if__211[k]
                  + f_160 * if__216[k]
                  - f_99 * if__218[k]
                  - f_75 * if__231[k]
                  - f_75 * if__236[k]
                  + f_76 * if__238[k]
                  + f_75 * if__251[k]
                  + f_75 * if__256[k]
                  - f_76 * if__258[k];
    }

#pragma omp simd aligned(if__2, if__7, if__9, if__32, if__37, if__39, if__52, if__57, if__59, \
                         if__102, if__107, if__109, if__142, if__147, if__149, if__212, \
                         if__217, if__219, if__232, if__237, if__239, if__252, if__257, \
                         if__259 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_59[k] = -f_161 * if__2[k]
                  - f_161 * if__7[k]
                  + f_162 * if__9[k]
                  - f_161 * if__32[k]
                  - f_161 * if__37[k]
                  + f_162 * if__39[k]
                  + f_84 * if__52[k]
                  + f_84 * if__57[k]
                  - f_85 * if__59[k]
                  + f_161 * if__102[k]
                  + f_161 * if__107[k]
                  - f_162 * if__109[k]
                  - f_84 * if__142[k]
                  - f_84 * if__147[k]
                  + f_85 * if__149[k]
                  + f_161 * if__212[k]
                  + f_161 * if__217[k]
                  - f_162 * if__219[k]
                  - f_84 * if__232[k]
                  - f_84 * if__237[k]
                  + f_85 * if__239[k]
                  + f_84 * if__252[k]
                  + f_84 * if__257[k]
                  - f_85 * if__259[k];
    }

#pragma omp simd aligned(if__0, if__3, if__5, if__30, if__33, if__35, if__50, if__53, if__55, \
                         if__100, if__103, if__105, if__140, if__143, if__145, if__210, \
                         if__213, if__215, if__230, if__233, if__235, if__250, if__253, \
                         if__255 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_60[k] = -f_160 * if__0[k]
                  - f_160 * if__3[k]
                  + f_99 * if__5[k]
                  - f_160 * if__30[k]
                  - f_160 * if__33[k]
                  + f_99 * if__35[k]
                  + f_75 * if__50[k]
                  + f_75 * if__53[k]
                  - f_76 * if__55[k]
                  + f_160 * if__100[k]
                  + f_160 * if__103[k]
                  - f_99 * if__105[k]
                  - f_75 * if__140[k]
                  - f_75 * if__143[k]
                  + f_76 * if__145[k]
                  + f_160 * if__210[k]
                  + f_160 * if__213[k]
                  - f_99 * if__215[k]
                  - f_75 * if__230[k]
                  - f_75 * if__233[k]
                  + f_76 * if__235[k]
                  + f_75 * if__250[k]
                  + f_75 * if__253[k]
                  - f_76 * if__255[k];
    }

#pragma omp simd aligned(if__2, if__7, if__32, if__37, if__52, if__57, if__102, if__107, \
                         if__142, if__147, if__212, if__217, if__232, if__237, if__252, \
                         if__257 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_61[k] = f_163 * if__2[k]
                  - f_163 * if__7[k]
                  + f_163 * if__32[k]
                  - f_163 * if__37[k]
                  - f_89 * if__52[k]
                  + f_89 * if__57[k]
                  - f_163 * if__102[k]
                  + f_163 * if__107[k]
                  + f_89 * if__142[k]
                  - f_89 * if__147[k]
                  - f_163 * if__212[k]
                  + f_163 * if__217[k]
                  + f_89 * if__232[k]
                  - f_89 * if__237[k]
                  - f_89 * if__252[k]
                  + f_89 * if__257[k];
    }

#pragma omp simd aligned(if__0, if__3, if__30, if__33, if__50, if__53, if__100, if__103, \
                         if__140, if__143, if__210, if__213, if__230, if__233, if__250, \
                         if__253 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_62[k] = f_159 * if__0[k]
                  - f_158 * if__3[k]
                  + f_159 * if__30[k]
                  - f_158 * if__33[k]
                  - f_62 * if__50[k]
                  + f_60 * if__53[k]
                  - f_159 * if__100[k]
                  + f_158 * if__103[k]
                  + f_62 * if__140[k]
                  - f_60 * if__143[k]
                  - f_159 * if__210[k]
                  + f_158 * if__213[k]
                  + f_62 * if__230[k]
                  - f_60 * if__233[k]
                  - f_62 * if__250[k]
                  + f_60 * if__253[k];
    }

#pragma omp simd aligned(if__21, if__26, if__71, if__76, if__91, if__96, if__161, if__166, \
                         if__181, if__186 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_63[k] = -f_56 * if__21[k]
                  + f_61 * if__26[k]
                  + f_57 * if__71[k]
                  - f_58 * if__76[k]
                  + f_60 * if__91[k]
                  - f_62 * if__96[k]
                  + f_55 * if__161[k]
                  - f_56 * if__166[k]
                  - f_59 * if__181[k]
                  + f_60 * if__186[k];
    }

#pragma omp simd aligned(if__24, if__74, if__94, if__164, if__184 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_64[k] = -f_66 * if__24[k]
                  + f_64 * if__74[k]
                  + f_67 * if__94[k]
                  + f_63 * if__164[k]
                  - f_65 * if__184[k];
    }

#pragma omp simd aligned(if__21, if__26, if__28, if__71, if__76, if__78, if__91, if__96, \
                         if__98, if__161, if__166, if__168, if__181, if__186, \
                         if__188 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_65[k] = f_73 * if__21[k]
                  + f_73 * if__26[k]
                  - f_74 * if__28[k]
                  - f_70 * if__71[k]
                  - f_70 * if__76[k]
                  + f_71 * if__78[k]
                  - f_75 * if__91[k]
                  - f_75 * if__96[k]
                  + f_76 * if__98[k]
                  - f_68 * if__161[k]
                  - f_68 * if__166[k]
                  + f_69 * if__168[k]
                  + f_71 * if__181[k]
                  + f_71 * if__186[k]
                  - f_72 * if__188[k];
    }

#pragma omp simd aligned(if__22, if__27, if__29, if__72, if__77, if__79, if__92, if__97, \
                         if__99, if__162, if__167, if__169, if__182, if__187, \
                         if__189 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_66[k] = f_82 * if__22[k]
                  + f_82 * if__27[k]
                  - f_83 * if__29[k]
                  - f_78 * if__72[k]
                  - f_78 * if__77[k]
                  + f_79 * if__79[k]
                  - f_84 * if__92[k]
                  - f_84 * if__97[k]
                  + f_85 * if__99[k]
                  - f_77 * if__162[k]
                  - f_77 * if__167[k]
                  + f_78 * if__169[k]
                  + f_80 * if__182[k]
                  + f_80 * if__187[k]
                  - f_81 * if__189[k];
    }

#pragma omp simd aligned(if__20, if__23, if__25, if__70, if__73, if__75, if__90, if__93, \
                         if__95, if__160, if__163, if__165, if__180, if__183, \
                         if__185 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_67[k] = f_73 * if__20[k]
                  + f_73 * if__23[k]
                  - f_74 * if__25[k]
                  - f_70 * if__70[k]
                  - f_70 * if__73[k]
                  + f_71 * if__75[k]
                  - f_75 * if__90[k]
                  - f_75 * if__93[k]
                  + f_76 * if__95[k]
                  - f_68 * if__160[k]
                  - f_68 * if__163[k]
                  + f_69 * if__165[k]
                  + f_71 * if__180[k]
                  + f_71 * if__183[k]
                  - f_72 * if__185[k];
    }

#pragma omp simd aligned(if__22, if__27, if__72, if__77, if__92, if__97, if__162, if__167, \
                         if__182, if__187 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_68[k] = -f_88 * if__22[k]
                  + f_88 * if__27[k]
                  + f_66 * if__72[k]
                  - f_66 * if__77[k]
                  + f_89 * if__92[k]
                  - f_89 * if__97[k]
                  + f_86 * if__162[k]
                  - f_86 * if__167[k]
                  - f_87 * if__182[k]
                  + f_87 * if__187[k];
    }

#pragma omp simd aligned(if__20, if__23, if__70, if__73, if__90, if__93, if__160, if__163, \
                         if__180, if__183 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_69[k] = -f_61 * if__20[k]
                  + f_56 * if__23[k]
                  + f_58 * if__70[k]
                  - f_57 * if__73[k]
                  + f_62 * if__90[k]
                  - f_60 * if__93[k]
                  + f_56 * if__160[k]
                  - f_55 * if__163[k]
                  - f_60 * if__180[k]
                  + f_59 * if__183[k];
    }

#pragma omp simd aligned(if__1, if__6, if__31, if__36, if__51, if__56, if__101, if__106, \
                         if__121, if__126, if__211, if__216, if__231, \
                         if__236 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_70[k] = -f_164 * if__1[k]
                  + f_165 * if__6[k]
                  + f_166 * if__31[k]
                  - f_167 * if__36[k]
                  + f_168 * if__51[k]
                  - f_169 * if__56[k]
                  + f_166 * if__101[k]
                  - f_167 * if__106[k]
                  - f_170 * if__121[k]
                  + f_171 * if__126[k]
                  - f_164 * if__211[k]
                  + f_165 * if__216[k]
                  + f_168 * if__231[k]
                  - f_169 * if__236[k];
    }

#pragma omp simd aligned(if__4, if__34, if__54, if__104, if__124, if__214, \
                         if__234 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_71[k] = -f_172 * if__4[k]
                  + f_173 * if__34[k]
                  + f_174 * if__54[k]
                  + f_173 * if__104[k]
                  - f_175 * if__124[k]
                  - f_172 * if__214[k]
                  + f_174 * if__234[k];
    }

#pragma omp simd aligned(if__1, if__6, if__8, if__31, if__36, if__38, if__51, if__56, if__58, \
                         if__101, if__106, if__108, if__121, if__126, if__128, if__211, \
                         if__216, if__218, if__231, if__236, if__238 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_72[k] = f_176 * if__1[k]
                  + f_176 * if__6[k]
                  - f_45 * if__8[k]
                  - f_177 * if__31[k]
                  - f_177 * if__36[k]
                  + f_178 * if__38[k]
                  - f_179 * if__51[k]
                  - f_179 * if__56[k]
                  + f_47 * if__58[k]
                  - f_177 * if__101[k]
                  - f_177 * if__106[k]
                  + f_178 * if__108[k]
                  + f_180 * if__121[k]
                  + f_180 * if__126[k]
                  - f_181 * if__128[k]
                  + f_176 * if__211[k]
                  + f_176 * if__216[k]
                  - f_45 * if__218[k]
                  - f_179 * if__231[k]
                  - f_179 * if__236[k]
                  + f_47 * if__238[k];
    }

#pragma omp simd aligned(if__2, if__7, if__9, if__32, if__37, if__39, if__52, if__57, if__59, \
                         if__102, if__107, if__109, if__122, if__127, if__129, if__212, \
                         if__217, if__219, if__232, if__237, if__239 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_73[k] = f_182 * if__2[k]
                  + f_182 * if__7[k]
                  - f_183 * if__9[k]
                  - f_184 * if__32[k]
                  - f_184 * if__37[k]
                  + f_185 * if__39[k]
                  - f_186 * if__52[k]
                  - f_186 * if__57[k]
                  + f_187 * if__59[k]
                  - f_184 * if__102[k]
                  - f_184 * if__107[k]
                  + f_185 * if__109[k]
                  + f_188 * if__122[k]
                  + f_188 * if__127[k]
                  - f_51 * if__129[k]
                  + f_182 * if__212[k]
                  + f_182 * if__217[k]
                  - f_183 * if__219[k]
                  - f_186 * if__232[k]
                  - f_186 * if__237[k]
                  + f_187 * if__239[k];
    }

#pragma omp simd aligned(if__0, if__3, if__5, if__30, if__33, if__35, if__50, if__53, if__55, \
                         if__100, if__103, if__105, if__120, if__123, if__125, if__210, \
                         if__213, if__215, if__230, if__233, if__235 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_74[k] = f_176 * if__0[k]
                  + f_176 * if__3[k]
                  - f_45 * if__5[k]
                  - f_177 * if__30[k]
                  - f_177 * if__33[k]
                  + f_178 * if__35[k]
                  - f_179 * if__50[k]
                  - f_179 * if__53[k]
                  + f_47 * if__55[k]
                  - f_177 * if__100[k]
                  - f_177 * if__103[k]
                  + f_178 * if__105[k]
                  + f_180 * if__120[k]
                  + f_180 * if__123[k]
                  - f_181 * if__125[k]
                  + f_176 * if__210[k]
                  + f_176 * if__213[k]
                  - f_45 * if__215[k]
                  - f_179 * if__230[k]
                  - f_179 * if__233[k]
                  + f_47 * if__235[k];
    }

#pragma omp simd aligned(if__2, if__7, if__32, if__37, if__52, if__57, if__102, if__107, \
                         if__122, if__127, if__212, if__217, if__232, \
                         if__237 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_75[k] = -f_189 * if__2[k]
                  + f_189 * if__7[k]
                  + f_190 * if__32[k]
                  - f_190 * if__37[k]
                  + f_173 * if__52[k]
                  - f_173 * if__57[k]
                  + f_190 * if__102[k]
                  - f_190 * if__107[k]
                  - f_191 * if__122[k]
                  + f_191 * if__127[k]
                  - f_189 * if__212[k]
                  + f_189 * if__217[k]
                  + f_173 * if__232[k]
                  - f_173 * if__237[k];
    }

#pragma omp simd aligned(if__0, if__3, if__30, if__33, if__50, if__53, if__100, if__103, \
                         if__120, if__123, if__210, if__213, if__230, \
                         if__233 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_76[k] = -f_165 * if__0[k]
                  + f_164 * if__3[k]
                  + f_167 * if__30[k]
                  - f_166 * if__33[k]
                  + f_169 * if__50[k]
                  - f_168 * if__53[k]
                  + f_167 * if__100[k]
                  - f_166 * if__103[k]
                  - f_171 * if__120[k]
                  + f_170 * if__123[k]
                  - f_165 * if__210[k]
                  + f_164 * if__213[k]
                  + f_169 * if__230[k]
                  - f_168 * if__233[k];
    }

#pragma omp simd aligned(if__21, if__24, if__26, if__28, if__71, if__74, if__76, if__78, \
                         if__161, if__164, if__166, if__168 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_77[k] = f_20 * if__21[k]
                  - f_21 * if__26[k]
                  - f_18 * if__71[k]
                  + f_19 * if__76[k]
                  + f_16 * if__161[k]
                  - f_17 * if__166[k];

        g_78[k] = f_24 * if__24[k]
                  - f_23 * if__74[k]
                  + f_22 * if__164[k];

        g_79[k] = -f_29 * if__21[k]
                  - f_29 * if__26[k]
                  + f_30 * if__28[k]
                  + f_27 * if__71[k]
                  + f_27 * if__76[k]
                  - f_28 * if__78[k]
                  - f_25 * if__161[k]
                  - f_25 * if__166[k]
                  + f_26 * if__168[k];
    }

#pragma omp simd aligned(if__22, if__27, if__29, if__72, if__77, if__79, if__162, if__167, \
                         if__169 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_80[k] = -f_35 * if__22[k]
                  - f_35 * if__27[k]
                  + f_36 * if__29[k]
                  + f_33 * if__72[k]
                  + f_33 * if__77[k]
                  - f_34 * if__79[k]
                  - f_31 * if__162[k]
                  - f_31 * if__167[k]
                  + f_32 * if__169[k];
    }

#pragma omp simd aligned(if__20, if__23, if__25, if__70, if__73, if__75, if__160, if__163, \
                         if__165 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_81[k] = -f_29 * if__20[k]
                  - f_29 * if__23[k]
                  + f_30 * if__25[k]
                  + f_27 * if__70[k]
                  + f_27 * if__73[k]
                  - f_28 * if__75[k]
                  - f_25 * if__160[k]
                  - f_25 * if__163[k]
                  + f_26 * if__165[k];
    }

#pragma omp simd aligned(if__20, if__22, if__23, if__27, if__70, if__72, if__73, if__77, \
                         if__160, if__162, if__163, if__167 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_82[k] = f_38 * if__22[k]
                  - f_38 * if__27[k]
                  - f_22 * if__72[k]
                  + f_22 * if__77[k]
                  + f_37 * if__162[k]
                  - f_37 * if__167[k];

        g_83[k] = f_21 * if__20[k]
                  - f_20 * if__23[k]
                  - f_19 * if__70[k]
                  + f_18 * if__73[k]
                  + f_17 * if__160[k]
                  - f_16 * if__163[k];
    }

#pragma omp simd aligned(if__1, if__4, if__6, if__31, if__34, if__36, if__101, if__104, \
                         if__106, if__211, if__214, if__216 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_84[k] = f_192 * if__1[k]
                  - f_193 * if__6[k]
                  - f_194 * if__31[k]
                  + f_195 * if__36[k]
                  + f_194 * if__101[k]
                  - f_195 * if__106[k]
                  - f_192 * if__211[k]
                  + f_193 * if__216[k];

        g_85[k] = f_196 * if__4[k]
                  - f_197 * if__34[k]
                  + f_197 * if__104[k]
                  - f_196 * if__214[k];
    }

#pragma omp simd aligned(if__1, if__6, if__8, if__31, if__36, if__38, if__101, if__106, \
                         if__108, if__211, if__216, if__218 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_86[k] = -f_198 * if__1[k]
                  - f_198 * if__6[k]
                  + f_199 * if__8[k]
                  + f_200 * if__31[k]
                  + f_200 * if__36[k]
                  - f_201 * if__38[k]
                  - f_200 * if__101[k]
                  - f_200 * if__106[k]
                  + f_201 * if__108[k]
                  + f_198 * if__211[k]
                  + f_198 * if__216[k]
                  - f_199 * if__218[k];
    }

#pragma omp simd aligned(if__2, if__7, if__9, if__32, if__37, if__39, if__102, if__107, \
                         if__109, if__212, if__217, if__219 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_87[k] = -f_202 * if__2[k]
                  - f_202 * if__7[k]
                  + f_203 * if__9[k]
                  + f_204 * if__32[k]
                  + f_204 * if__37[k]
                  - f_205 * if__39[k]
                  - f_204 * if__102[k]
                  - f_204 * if__107[k]
                  + f_205 * if__109[k]
                  + f_202 * if__212[k]
                  + f_202 * if__217[k]
                  - f_203 * if__219[k];
    }

#pragma omp simd aligned(if__0, if__3, if__5, if__30, if__33, if__35, if__100, if__103, \
                         if__105, if__210, if__213, if__215 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_88[k] = -f_198 * if__0[k]
                  - f_198 * if__3[k]
                  + f_199 * if__5[k]
                  + f_200 * if__30[k]
                  + f_200 * if__33[k]
                  - f_201 * if__35[k]
                  - f_200 * if__100[k]
                  - f_200 * if__103[k]
                  + f_201 * if__105[k]
                  + f_198 * if__210[k]
                  + f_198 * if__213[k]
                  - f_199 * if__215[k];
    }

#pragma omp simd aligned(if__2, if__7, if__32, if__37, if__102, if__107, if__212, \
                         if__217 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_89[k] = f_206 * if__2[k]
                  - f_206 * if__7[k]
                  - f_207 * if__32[k]
                  + f_207 * if__37[k]
                  + f_207 * if__102[k]
                  - f_207 * if__107[k]
                  - f_206 * if__212[k]
                  + f_206 * if__217[k];
    }

#pragma omp simd aligned(if__0, if__3, if__30, if__33, if__100, if__103, if__210, \
                         if__213 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_90[k] = f_193 * if__0[k]
                  - f_192 * if__3[k]
                  - f_195 * if__30[k]
                  + f_194 * if__33[k]
                  + f_195 * if__100[k]
                  - f_194 * if__103[k]
                  - f_193 * if__210[k]
                  + f_192 * if__213[k];
    }
}

}  // namespace simdtrf
