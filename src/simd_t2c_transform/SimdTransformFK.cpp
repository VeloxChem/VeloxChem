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


#include "SimdTransformFK.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_fk(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t fk,
             const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 0.1640625 * std::sqrt(4290.0);
    const auto f_1 = 0.8203125 * std::sqrt(4290.0);
    const auto f_2 = 0.4921875 * std::sqrt(4290.0);
    const auto f_3 = 0.0234375 * std::sqrt(4290.0);
    const auto f_4 = 0.0546875 * std::sqrt(4290.0);
    const auto f_5 = 0.2734375 * std::sqrt(4290.0);
    const auto f_6 = 0.0078125 * std::sqrt(4290.0);
    const auto f_7 = 0.28125 * std::sqrt(15015.0);
    const auto f_8 = 0.9375 * std::sqrt(15015.0);
    const auto f_9 = 0.09375 * std::sqrt(15015.0);
    const auto f_10 = 0.3125 * std::sqrt(15015.0);
    const auto f_11 = 0.1171875 * std::sqrt(2310.0);
    const auto f_12 = 1.40625 * std::sqrt(2310.0);
    const auto f_13 = 0.2109375 * std::sqrt(2310.0);
    const auto f_14 = 2.8125 * std::sqrt(2310.0);
    const auto f_15 = 0.0234375 * std::sqrt(2310.0);
    const auto f_16 = 0.28125 * std::sqrt(2310.0);
    const auto f_17 = 0.0390625 * std::sqrt(2310.0);
    const auto f_18 = 0.46875 * std::sqrt(2310.0);
    const auto f_19 = 0.0703125 * std::sqrt(2310.0);
    const auto f_20 = 0.9375 * std::sqrt(2310.0);
    const auto f_21 = 0.0078125 * std::sqrt(2310.0);
    const auto f_22 = 0.09375 * std::sqrt(2310.0);
    const auto f_23 = 0.5625 * std::sqrt(2310.0);
    const auto f_24 = 1.875 * std::sqrt(2310.0);
    const auto f_25 = 0.1875 * std::sqrt(2310.0);
    const auto f_26 = 0.625 * std::sqrt(2310.0);
    const auto f_27 = 0.2109375 * std::sqrt(210.0);
    const auto f_28 = 0.3515625 * std::sqrt(210.0);
    const auto f_29 = 4.21875 * std::sqrt(210.0);
    const auto f_30 = 0.0703125 * std::sqrt(210.0);
    const auto f_31 = 2.8125 * std::sqrt(210.0);
    const auto f_32 = 5.625 * std::sqrt(210.0);
    const auto f_33 = 1.40625 * std::sqrt(210.0);
    const auto f_34 = 1.875 * std::sqrt(210.0);
    const auto f_35 = 0.1171875 * std::sqrt(210.0);
    const auto f_36 = 0.0234375 * std::sqrt(210.0);
    const auto f_37 = 0.9375 * std::sqrt(210.0);
    const auto f_38 = 0.46875 * std::sqrt(210.0);
    const auto f_39 = 0.625 * std::sqrt(210.0);
    const auto f_40 = 1.40625 * std::sqrt(105.0);
    const auto f_41 = 2.8125 * std::sqrt(105.0);
    const auto f_42 = 7.5 * std::sqrt(105.0);
    const auto f_43 = 4.5 * std::sqrt(105.0);
    const auto f_44 = 0.46875 * std::sqrt(105.0);
    const auto f_45 = 0.9375 * std::sqrt(105.0);
    const auto f_46 = 2.5 * std::sqrt(105.0);
    const auto f_47 = 1.5 * std::sqrt(105.0);
    const auto f_48 = 0.1171875 * std::sqrt(70.0);
    const auto f_49 = 0.3515625 * std::sqrt(70.0);
    const auto f_50 = 2.8125 * std::sqrt(70.0);
    const auto f_51 = 5.625 * std::sqrt(70.0);
    const auto f_52 = 1.5 * std::sqrt(70.0);
    const auto f_53 = 0.0390625 * std::sqrt(70.0);
    const auto f_54 = 0.9375 * std::sqrt(70.0);
    const auto f_55 = 1.875 * std::sqrt(70.0);
    const auto f_56 = 0.5 * std::sqrt(70.0);
    const auto f_57 = 1.640625 * std::sqrt(10.0);
    const auto f_58 = 4.921875 * std::sqrt(10.0);
    const auto f_59 = 9.84375 * std::sqrt(10.0);
    const auto f_60 = 19.6875 * std::sqrt(10.0);
    const auto f_61 = 7.875 * std::sqrt(10.0);
    const auto f_62 = 0.75 * std::sqrt(10.0);
    const auto f_63 = 0.546875 * std::sqrt(10.0);
    const auto f_64 = 3.28125 * std::sqrt(10.0);
    const auto f_65 = 6.5625 * std::sqrt(10.0);
    const auto f_66 = 2.625 * std::sqrt(10.0);
    const auto f_67 = 0.25 * std::sqrt(10.0);
    const auto f_68 = 0.703125 * std::sqrt(105.0);
    const auto f_69 = 3.75 * std::sqrt(105.0);
    const auto f_70 = 2.25 * std::sqrt(105.0);
    const auto f_71 = 0.234375 * std::sqrt(105.0);
    const auto f_72 = 1.25 * std::sqrt(105.0);
    const auto f_73 = 0.75 * std::sqrt(105.0);
    const auto f_74 = 0.140625 * std::sqrt(2310.0);
    const auto f_75 = 0.703125 * std::sqrt(2310.0);
    const auto f_76 = 0.046875 * std::sqrt(2310.0);
    const auto f_77 = 0.234375 * std::sqrt(2310.0);
    const auto f_78 = 0.15625 * std::sqrt(2310.0);
    const auto f_79 = 0.046875 * std::sqrt(15015.0);
    const auto f_80 = 0.703125 * std::sqrt(15015.0);
    const auto f_81 = 0.015625 * std::sqrt(15015.0);
    const auto f_82 = 0.234375 * std::sqrt(15015.0);
    const auto f_83 = 0.65625 * std::sqrt(715.0);
    const auto f_84 = 3.28125 * std::sqrt(715.0);
    const auto f_85 = 1.96875 * std::sqrt(715.0);
    const auto f_86 = 0.09375 * std::sqrt(715.0);
    const auto f_87 = 0.5625 * std::sqrt(10010.0);
    const auto f_88 = 1.875 * std::sqrt(10010.0);
    const auto f_89 = 0.46875 * std::sqrt(385.0);
    const auto f_90 = 5.625 * std::sqrt(385.0);
    const auto f_91 = 0.84375 * std::sqrt(385.0);
    const auto f_92 = 11.25 * std::sqrt(385.0);
    const auto f_93 = 0.09375 * std::sqrt(385.0);
    const auto f_94 = 1.125 * std::sqrt(385.0);
    const auto f_95 = 2.25 * std::sqrt(385.0);
    const auto f_96 = 7.5 * std::sqrt(385.0);
    const auto f_97 = 0.84375 * std::sqrt(35.0);
    const auto f_98 = 1.40625 * std::sqrt(35.0);
    const auto f_99 = 16.875 * std::sqrt(35.0);
    const auto f_100 = 0.28125 * std::sqrt(35.0);
    const auto f_101 = 11.25 * std::sqrt(35.0);
    const auto f_102 = 22.5 * std::sqrt(35.0);
    const auto f_103 = 5.625 * std::sqrt(35.0);
    const auto f_104 = 7.5 * std::sqrt(35.0);
    const auto f_105 = 15.0 * std::sqrt(70.0);
    const auto f_106 = 9.0 * std::sqrt(70.0);
    const auto f_107 = 0.15625 * std::sqrt(105.0);
    const auto f_108 = 2.0 * std::sqrt(105.0);
    const auto f_109 = 2.1875 * std::sqrt(15.0);
    const auto f_110 = 6.5625 * std::sqrt(15.0);
    const auto f_111 = 13.125 * std::sqrt(15.0);
    const auto f_112 = 26.25 * std::sqrt(15.0);
    const auto f_113 = 10.5 * std::sqrt(15.0);
    const auto f_114 = std::sqrt(15.0);
    const auto f_115 = 1.40625 * std::sqrt(70.0);
    const auto f_116 = 7.5 * std::sqrt(70.0);
    const auto f_117 = 4.5 * std::sqrt(70.0);
    const auto f_118 = 0.5625 * std::sqrt(385.0);
    const auto f_119 = 2.8125 * std::sqrt(385.0);
    const auto f_120 = 1.875 * std::sqrt(385.0);
    const auto f_121 = 0.09375 * std::sqrt(10010.0);
    const auto f_122 = 1.40625 * std::sqrt(10010.0);
    const auto f_123 = 0.1640625 * std::sqrt(286.0);
    const auto f_124 = 0.8203125 * std::sqrt(286.0);
    const auto f_125 = 0.4921875 * std::sqrt(286.0);
    const auto f_126 = 0.0234375 * std::sqrt(286.0);
    const auto f_127 = 0.65625 * std::sqrt(286.0);
    const auto f_128 = 3.28125 * std::sqrt(286.0);
    const auto f_129 = 1.96875 * std::sqrt(286.0);
    const auto f_130 = 0.09375 * std::sqrt(286.0);
    const auto f_131 = 0.28125 * std::sqrt(1001.0);
    const auto f_132 = 0.9375 * std::sqrt(1001.0);
    const auto f_133 = 1.125 * std::sqrt(1001.0);
    const auto f_134 = 3.75 * std::sqrt(1001.0);
    const auto f_135 = 0.1171875 * std::sqrt(154.0);
    const auto f_136 = 1.40625 * std::sqrt(154.0);
    const auto f_137 = 0.2109375 * std::sqrt(154.0);
    const auto f_138 = 2.8125 * std::sqrt(154.0);
    const auto f_139 = 0.0234375 * std::sqrt(154.0);
    const auto f_140 = 0.28125 * std::sqrt(154.0);
    const auto f_141 = 0.46875 * std::sqrt(154.0);
    const auto f_142 = 5.625 * std::sqrt(154.0);
    const auto f_143 = 0.84375 * std::sqrt(154.0);
    const auto f_144 = 11.25 * std::sqrt(154.0);
    const auto f_145 = 0.09375 * std::sqrt(154.0);
    const auto f_146 = 1.125 * std::sqrt(154.0);
    const auto f_147 = 0.5625 * std::sqrt(154.0);
    const auto f_148 = 1.875 * std::sqrt(154.0);
    const auto f_149 = 2.25 * std::sqrt(154.0);
    const auto f_150 = 7.5 * std::sqrt(154.0);
    const auto f_151 = 0.2109375 * std::sqrt(14.0);
    const auto f_152 = 0.3515625 * std::sqrt(14.0);
    const auto f_153 = 4.21875 * std::sqrt(14.0);
    const auto f_154 = 0.0703125 * std::sqrt(14.0);
    const auto f_155 = 2.8125 * std::sqrt(14.0);
    const auto f_156 = 5.625 * std::sqrt(14.0);
    const auto f_157 = 1.40625 * std::sqrt(14.0);
    const auto f_158 = 1.875 * std::sqrt(14.0);
    const auto f_159 = 0.84375 * std::sqrt(14.0);
    const auto f_160 = 16.875 * std::sqrt(14.0);
    const auto f_161 = 0.28125 * std::sqrt(14.0);
    const auto f_162 = 11.25 * std::sqrt(14.0);
    const auto f_163 = 22.5 * std::sqrt(14.0);
    const auto f_164 = 7.5 * std::sqrt(14.0);
    const auto f_165 = 1.40625 * std::sqrt(7.0);
    const auto f_166 = 2.8125 * std::sqrt(7.0);
    const auto f_167 = 7.5 * std::sqrt(7.0);
    const auto f_168 = 4.5 * std::sqrt(7.0);
    const auto f_169 = 5.625 * std::sqrt(7.0);
    const auto f_170 = 11.25 * std::sqrt(7.0);
    const auto f_171 = 30.0 * std::sqrt(7.0);
    const auto f_172 = 18.0 * std::sqrt(7.0);
    const auto f_173 = 0.0390625 * std::sqrt(42.0);
    const auto f_174 = 0.1171875 * std::sqrt(42.0);
    const auto f_175 = 0.9375 * std::sqrt(42.0);
    const auto f_176 = 1.875 * std::sqrt(42.0);
    const auto f_177 = 0.5 * std::sqrt(42.0);
    const auto f_178 = 0.15625 * std::sqrt(42.0);
    const auto f_179 = 0.46875 * std::sqrt(42.0);
    const auto f_180 = 3.75 * std::sqrt(42.0);
    const auto f_181 = 7.5 * std::sqrt(42.0);
    const auto f_182 = 2.0 * std::sqrt(42.0);
    const auto f_183 = 0.546875 * std::sqrt(6.0);
    const auto f_184 = 1.640625 * std::sqrt(6.0);
    const auto f_185 = 3.28125 * std::sqrt(6.0);
    const auto f_186 = 6.5625 * std::sqrt(6.0);
    const auto f_187 = 2.625 * std::sqrt(6.0);
    const auto f_188 = 0.25 * std::sqrt(6.0);
    const auto f_189 = 2.1875 * std::sqrt(6.0);
    const auto f_190 = 13.125 * std::sqrt(6.0);
    const auto f_191 = 26.25 * std::sqrt(6.0);
    const auto f_192 = 10.5 * std::sqrt(6.0);
    const auto f_193 = std::sqrt(6.0);
    const auto f_194 = 0.703125 * std::sqrt(7.0);
    const auto f_195 = 3.75 * std::sqrt(7.0);
    const auto f_196 = 2.25 * std::sqrt(7.0);
    const auto f_197 = 15.0 * std::sqrt(7.0);
    const auto f_198 = 9.0 * std::sqrt(7.0);
    const auto f_199 = 0.140625 * std::sqrt(154.0);
    const auto f_200 = 0.703125 * std::sqrt(154.0);
    const auto f_201 = 0.046875 * std::sqrt(1001.0);
    const auto f_202 = 0.703125 * std::sqrt(1001.0);
    const auto f_203 = 0.1875 * std::sqrt(1001.0);
    const auto f_204 = 2.8125 * std::sqrt(1001.0);
    const auto f_205 = 0.328125 * std::sqrt(429.0);
    const auto f_206 = 1.640625 * std::sqrt(429.0);
    const auto f_207 = 0.984375 * std::sqrt(429.0);
    const auto f_208 = 0.046875 * std::sqrt(429.0);
    const auto f_209 = 0.21875 * std::sqrt(429.0);
    const auto f_210 = 1.09375 * std::sqrt(429.0);
    const auto f_211 = 0.65625 * std::sqrt(429.0);
    const auto f_212 = 0.03125 * std::sqrt(429.0);
    const auto f_213 = 0.28125 * std::sqrt(6006.0);
    const auto f_214 = 0.9375 * std::sqrt(6006.0);
    const auto f_215 = 0.1875 * std::sqrt(6006.0);
    const auto f_216 = 0.625 * std::sqrt(6006.0);
    const auto f_217 = 0.234375 * std::sqrt(231.0);
    const auto f_218 = 2.8125 * std::sqrt(231.0);
    const auto f_219 = 0.421875 * std::sqrt(231.0);
    const auto f_220 = 5.625 * std::sqrt(231.0);
    const auto f_221 = 0.046875 * std::sqrt(231.0);
    const auto f_222 = 0.5625 * std::sqrt(231.0);
    const auto f_223 = 0.15625 * std::sqrt(231.0);
    const auto f_224 = 1.875 * std::sqrt(231.0);
    const auto f_225 = 0.28125 * std::sqrt(231.0);
    const auto f_226 = 3.75 * std::sqrt(231.0);
    const auto f_227 = 0.03125 * std::sqrt(231.0);
    const auto f_228 = 0.375 * std::sqrt(231.0);
    const auto f_229 = 1.125 * std::sqrt(231.0);
    const auto f_230 = 0.75 * std::sqrt(231.0);
    const auto f_231 = 2.5 * std::sqrt(231.0);
    const auto f_232 = 0.421875 * std::sqrt(21.0);
    const auto f_233 = 0.703125 * std::sqrt(21.0);
    const auto f_234 = 8.4375 * std::sqrt(21.0);
    const auto f_235 = 0.140625 * std::sqrt(21.0);
    const auto f_236 = 5.625 * std::sqrt(21.0);
    const auto f_237 = 11.25 * std::sqrt(21.0);
    const auto f_238 = 2.8125 * std::sqrt(21.0);
    const auto f_239 = 3.75 * std::sqrt(21.0);
    const auto f_240 = 0.28125 * std::sqrt(21.0);
    const auto f_241 = 0.46875 * std::sqrt(21.0);
    const auto f_242 = 0.09375 * std::sqrt(21.0);
    const auto f_243 = 7.5 * std::sqrt(21.0);
    const auto f_244 = 1.875 * std::sqrt(21.0);
    const auto f_245 = 2.5 * std::sqrt(21.0);
    const auto f_246 = 1.40625 * std::sqrt(42.0);
    const auto f_247 = 2.8125 * std::sqrt(42.0);
    const auto f_248 = 4.5 * std::sqrt(42.0);
    const auto f_249 = 5.0 * std::sqrt(42.0);
    const auto f_250 = 3.0 * std::sqrt(42.0);
    const auto f_251 = 0.234375 * std::sqrt(7.0);
    const auto f_252 = 3.0 * std::sqrt(7.0);
    const auto f_253 = 0.15625 * std::sqrt(7.0);
    const auto f_254 = 0.46875 * std::sqrt(7.0);
    const auto f_255 = 2.0 * std::sqrt(7.0);
    const auto f_256 = 0.703125 * std::sqrt(42.0);
    const auto f_257 = 2.25 * std::sqrt(42.0);
    const auto f_258 = 2.5 * std::sqrt(42.0);
    const auto f_259 = 1.5 * std::sqrt(42.0);
    const auto f_260 = 1.40625 * std::sqrt(231.0);
    const auto f_261 = 0.9375 * std::sqrt(231.0);
    const auto f_262 = 0.1875 * std::sqrt(231.0);
    const auto f_263 = 0.625 * std::sqrt(231.0);
    const auto f_264 = 0.046875 * std::sqrt(6006.0);
    const auto f_265 = 0.703125 * std::sqrt(6006.0);
    const auto f_266 = 0.03125 * std::sqrt(6006.0);
    const auto f_267 = 0.46875 * std::sqrt(6006.0);
    const auto f_268 = 0.328125 * std::sqrt(715.0);
    const auto f_269 = 1.640625 * std::sqrt(715.0);
    const auto f_270 = 0.984375 * std::sqrt(715.0);
    const auto f_271 = 0.046875 * std::sqrt(715.0);
    const auto f_272 = 0.28125 * std::sqrt(10010.0);
    const auto f_273 = 0.9375 * std::sqrt(10010.0);
    const auto f_274 = 0.234375 * std::sqrt(385.0);
    const auto f_275 = 0.421875 * std::sqrt(385.0);
    const auto f_276 = 0.046875 * std::sqrt(385.0);
    const auto f_277 = 3.75 * std::sqrt(385.0);
    const auto f_278 = 0.421875 * std::sqrt(35.0);
    const auto f_279 = 0.703125 * std::sqrt(35.0);
    const auto f_280 = 8.4375 * std::sqrt(35.0);
    const auto f_281 = 0.140625 * std::sqrt(35.0);
    const auto f_282 = 2.8125 * std::sqrt(35.0);
    const auto f_283 = 3.75 * std::sqrt(35.0);
    const auto f_284 = 0.078125 * std::sqrt(105.0);
    const auto f_285 = 1.875 * std::sqrt(105.0);
    const auto f_286 = std::sqrt(105.0);
    const auto f_287 = 1.09375 * std::sqrt(15.0);
    const auto f_288 = 3.28125 * std::sqrt(15.0);
    const auto f_289 = 5.25 * std::sqrt(15.0);
    const auto f_290 = 0.5 * std::sqrt(15.0);
    const auto f_291 = 0.703125 * std::sqrt(70.0);
    const auto f_292 = 3.75 * std::sqrt(70.0);
    const auto f_293 = 2.25 * std::sqrt(70.0);
    const auto f_294 = 0.28125 * std::sqrt(385.0);
    const auto f_295 = 1.40625 * std::sqrt(385.0);
    const auto f_296 = 0.9375 * std::sqrt(385.0);
    const auto f_297 = 0.046875 * std::sqrt(10010.0);
    const auto f_298 = 0.703125 * std::sqrt(10010.0);

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
    auto *g_99 = values + 99 * nvalues;
    auto *g_100 = values + 100 * nvalues;
    auto *g_101 = values + 101 * nvalues;
    auto *g_102 = values + 102 * nvalues;
    auto *g_103 = values + 103 * nvalues;
    auto *g_104 = values + 104 * nvalues;

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_1 = buffer.data(fk + 1);
    const auto *fk_2 = buffer.data(fk + 2);
    const auto *fk_3 = buffer.data(fk + 3);
    const auto *fk_4 = buffer.data(fk + 4);
    const auto *fk_5 = buffer.data(fk + 5);
    const auto *fk_6 = buffer.data(fk + 6);
    const auto *fk_7 = buffer.data(fk + 7);
    const auto *fk_8 = buffer.data(fk + 8);
    const auto *fk_9 = buffer.data(fk + 9);
    const auto *fk_10 = buffer.data(fk + 10);
    const auto *fk_11 = buffer.data(fk + 11);
    const auto *fk_12 = buffer.data(fk + 12);
    const auto *fk_13 = buffer.data(fk + 13);
    const auto *fk_14 = buffer.data(fk + 14);
    const auto *fk_15 = buffer.data(fk + 15);
    const auto *fk_16 = buffer.data(fk + 16);
    const auto *fk_17 = buffer.data(fk + 17);
    const auto *fk_18 = buffer.data(fk + 18);
    const auto *fk_19 = buffer.data(fk + 19);
    const auto *fk_20 = buffer.data(fk + 20);
    const auto *fk_21 = buffer.data(fk + 21);
    const auto *fk_22 = buffer.data(fk + 22);
    const auto *fk_23 = buffer.data(fk + 23);
    const auto *fk_24 = buffer.data(fk + 24);
    const auto *fk_25 = buffer.data(fk + 25);
    const auto *fk_26 = buffer.data(fk + 26);
    const auto *fk_27 = buffer.data(fk + 27);
    const auto *fk_28 = buffer.data(fk + 28);
    const auto *fk_29 = buffer.data(fk + 29);
    const auto *fk_30 = buffer.data(fk + 30);
    const auto *fk_31 = buffer.data(fk + 31);
    const auto *fk_32 = buffer.data(fk + 32);
    const auto *fk_33 = buffer.data(fk + 33);
    const auto *fk_34 = buffer.data(fk + 34);
    const auto *fk_35 = buffer.data(fk + 35);
    const auto *fk_36 = buffer.data(fk + 36);
    const auto *fk_37 = buffer.data(fk + 37);
    const auto *fk_38 = buffer.data(fk + 38);
    const auto *fk_39 = buffer.data(fk + 39);
    const auto *fk_40 = buffer.data(fk + 40);
    const auto *fk_41 = buffer.data(fk + 41);
    const auto *fk_42 = buffer.data(fk + 42);
    const auto *fk_43 = buffer.data(fk + 43);
    const auto *fk_44 = buffer.data(fk + 44);
    const auto *fk_45 = buffer.data(fk + 45);
    const auto *fk_46 = buffer.data(fk + 46);
    const auto *fk_47 = buffer.data(fk + 47);
    const auto *fk_48 = buffer.data(fk + 48);
    const auto *fk_49 = buffer.data(fk + 49);
    const auto *fk_50 = buffer.data(fk + 50);
    const auto *fk_51 = buffer.data(fk + 51);
    const auto *fk_52 = buffer.data(fk + 52);
    const auto *fk_53 = buffer.data(fk + 53);
    const auto *fk_54 = buffer.data(fk + 54);
    const auto *fk_55 = buffer.data(fk + 55);
    const auto *fk_56 = buffer.data(fk + 56);
    const auto *fk_57 = buffer.data(fk + 57);
    const auto *fk_58 = buffer.data(fk + 58);
    const auto *fk_59 = buffer.data(fk + 59);
    const auto *fk_60 = buffer.data(fk + 60);
    const auto *fk_61 = buffer.data(fk + 61);
    const auto *fk_62 = buffer.data(fk + 62);
    const auto *fk_63 = buffer.data(fk + 63);
    const auto *fk_64 = buffer.data(fk + 64);
    const auto *fk_65 = buffer.data(fk + 65);
    const auto *fk_66 = buffer.data(fk + 66);
    const auto *fk_67 = buffer.data(fk + 67);
    const auto *fk_68 = buffer.data(fk + 68);
    const auto *fk_69 = buffer.data(fk + 69);
    const auto *fk_70 = buffer.data(fk + 70);
    const auto *fk_71 = buffer.data(fk + 71);
    const auto *fk_72 = buffer.data(fk + 72);
    const auto *fk_73 = buffer.data(fk + 73);
    const auto *fk_74 = buffer.data(fk + 74);
    const auto *fk_75 = buffer.data(fk + 75);
    const auto *fk_76 = buffer.data(fk + 76);
    const auto *fk_77 = buffer.data(fk + 77);
    const auto *fk_78 = buffer.data(fk + 78);
    const auto *fk_79 = buffer.data(fk + 79);
    const auto *fk_80 = buffer.data(fk + 80);
    const auto *fk_81 = buffer.data(fk + 81);
    const auto *fk_82 = buffer.data(fk + 82);
    const auto *fk_83 = buffer.data(fk + 83);
    const auto *fk_84 = buffer.data(fk + 84);
    const auto *fk_85 = buffer.data(fk + 85);
    const auto *fk_86 = buffer.data(fk + 86);
    const auto *fk_87 = buffer.data(fk + 87);
    const auto *fk_88 = buffer.data(fk + 88);
    const auto *fk_89 = buffer.data(fk + 89);
    const auto *fk_90 = buffer.data(fk + 90);
    const auto *fk_91 = buffer.data(fk + 91);
    const auto *fk_92 = buffer.data(fk + 92);
    const auto *fk_93 = buffer.data(fk + 93);
    const auto *fk_94 = buffer.data(fk + 94);
    const auto *fk_95 = buffer.data(fk + 95);
    const auto *fk_96 = buffer.data(fk + 96);
    const auto *fk_97 = buffer.data(fk + 97);
    const auto *fk_98 = buffer.data(fk + 98);
    const auto *fk_99 = buffer.data(fk + 99);
    const auto *fk_100 = buffer.data(fk + 100);
    const auto *fk_101 = buffer.data(fk + 101);
    const auto *fk_102 = buffer.data(fk + 102);
    const auto *fk_103 = buffer.data(fk + 103);
    const auto *fk_104 = buffer.data(fk + 104);
    const auto *fk_105 = buffer.data(fk + 105);
    const auto *fk_106 = buffer.data(fk + 106);
    const auto *fk_107 = buffer.data(fk + 107);
    const auto *fk_108 = buffer.data(fk + 108);
    const auto *fk_109 = buffer.data(fk + 109);
    const auto *fk_110 = buffer.data(fk + 110);
    const auto *fk_111 = buffer.data(fk + 111);
    const auto *fk_112 = buffer.data(fk + 112);
    const auto *fk_113 = buffer.data(fk + 113);
    const auto *fk_114 = buffer.data(fk + 114);
    const auto *fk_115 = buffer.data(fk + 115);
    const auto *fk_116 = buffer.data(fk + 116);
    const auto *fk_117 = buffer.data(fk + 117);
    const auto *fk_118 = buffer.data(fk + 118);
    const auto *fk_119 = buffer.data(fk + 119);
    const auto *fk_120 = buffer.data(fk + 120);
    const auto *fk_121 = buffer.data(fk + 121);
    const auto *fk_122 = buffer.data(fk + 122);
    const auto *fk_123 = buffer.data(fk + 123);
    const auto *fk_124 = buffer.data(fk + 124);
    const auto *fk_125 = buffer.data(fk + 125);
    const auto *fk_126 = buffer.data(fk + 126);
    const auto *fk_127 = buffer.data(fk + 127);
    const auto *fk_128 = buffer.data(fk + 128);
    const auto *fk_129 = buffer.data(fk + 129);
    const auto *fk_130 = buffer.data(fk + 130);
    const auto *fk_131 = buffer.data(fk + 131);
    const auto *fk_132 = buffer.data(fk + 132);
    const auto *fk_133 = buffer.data(fk + 133);
    const auto *fk_134 = buffer.data(fk + 134);
    const auto *fk_135 = buffer.data(fk + 135);
    const auto *fk_136 = buffer.data(fk + 136);
    const auto *fk_137 = buffer.data(fk + 137);
    const auto *fk_138 = buffer.data(fk + 138);
    const auto *fk_139 = buffer.data(fk + 139);
    const auto *fk_140 = buffer.data(fk + 140);
    const auto *fk_141 = buffer.data(fk + 141);
    const auto *fk_142 = buffer.data(fk + 142);
    const auto *fk_143 = buffer.data(fk + 143);
    const auto *fk_144 = buffer.data(fk + 144);
    const auto *fk_145 = buffer.data(fk + 145);
    const auto *fk_146 = buffer.data(fk + 146);
    const auto *fk_147 = buffer.data(fk + 147);
    const auto *fk_148 = buffer.data(fk + 148);
    const auto *fk_149 = buffer.data(fk + 149);
    const auto *fk_150 = buffer.data(fk + 150);
    const auto *fk_151 = buffer.data(fk + 151);
    const auto *fk_152 = buffer.data(fk + 152);
    const auto *fk_153 = buffer.data(fk + 153);
    const auto *fk_154 = buffer.data(fk + 154);
    const auto *fk_155 = buffer.data(fk + 155);
    const auto *fk_156 = buffer.data(fk + 156);
    const auto *fk_157 = buffer.data(fk + 157);
    const auto *fk_158 = buffer.data(fk + 158);
    const auto *fk_159 = buffer.data(fk + 159);
    const auto *fk_160 = buffer.data(fk + 160);
    const auto *fk_161 = buffer.data(fk + 161);
    const auto *fk_162 = buffer.data(fk + 162);
    const auto *fk_163 = buffer.data(fk + 163);
    const auto *fk_164 = buffer.data(fk + 164);
    const auto *fk_165 = buffer.data(fk + 165);
    const auto *fk_166 = buffer.data(fk + 166);
    const auto *fk_167 = buffer.data(fk + 167);
    const auto *fk_168 = buffer.data(fk + 168);
    const auto *fk_169 = buffer.data(fk + 169);
    const auto *fk_170 = buffer.data(fk + 170);
    const auto *fk_171 = buffer.data(fk + 171);
    const auto *fk_172 = buffer.data(fk + 172);
    const auto *fk_173 = buffer.data(fk + 173);
    const auto *fk_174 = buffer.data(fk + 174);
    const auto *fk_175 = buffer.data(fk + 175);
    const auto *fk_176 = buffer.data(fk + 176);
    const auto *fk_177 = buffer.data(fk + 177);
    const auto *fk_178 = buffer.data(fk + 178);
    const auto *fk_179 = buffer.data(fk + 179);
    const auto *fk_180 = buffer.data(fk + 180);
    const auto *fk_181 = buffer.data(fk + 181);
    const auto *fk_182 = buffer.data(fk + 182);
    const auto *fk_183 = buffer.data(fk + 183);
    const auto *fk_184 = buffer.data(fk + 184);
    const auto *fk_185 = buffer.data(fk + 185);
    const auto *fk_186 = buffer.data(fk + 186);
    const auto *fk_187 = buffer.data(fk + 187);
    const auto *fk_188 = buffer.data(fk + 188);
    const auto *fk_189 = buffer.data(fk + 189);
    const auto *fk_190 = buffer.data(fk + 190);
    const auto *fk_191 = buffer.data(fk + 191);
    const auto *fk_192 = buffer.data(fk + 192);
    const auto *fk_193 = buffer.data(fk + 193);
    const auto *fk_194 = buffer.data(fk + 194);
    const auto *fk_195 = buffer.data(fk + 195);
    const auto *fk_196 = buffer.data(fk + 196);
    const auto *fk_197 = buffer.data(fk + 197);
    const auto *fk_198 = buffer.data(fk + 198);
    const auto *fk_199 = buffer.data(fk + 199);
    const auto *fk_200 = buffer.data(fk + 200);
    const auto *fk_201 = buffer.data(fk + 201);
    const auto *fk_202 = buffer.data(fk + 202);
    const auto *fk_203 = buffer.data(fk + 203);
    const auto *fk_204 = buffer.data(fk + 204);
    const auto *fk_205 = buffer.data(fk + 205);
    const auto *fk_206 = buffer.data(fk + 206);
    const auto *fk_207 = buffer.data(fk + 207);
    const auto *fk_208 = buffer.data(fk + 208);
    const auto *fk_209 = buffer.data(fk + 209);
    const auto *fk_210 = buffer.data(fk + 210);
    const auto *fk_211 = buffer.data(fk + 211);
    const auto *fk_212 = buffer.data(fk + 212);
    const auto *fk_213 = buffer.data(fk + 213);
    const auto *fk_214 = buffer.data(fk + 214);
    const auto *fk_215 = buffer.data(fk + 215);
    const auto *fk_216 = buffer.data(fk + 216);
    const auto *fk_217 = buffer.data(fk + 217);
    const auto *fk_218 = buffer.data(fk + 218);
    const auto *fk_219 = buffer.data(fk + 219);
    const auto *fk_220 = buffer.data(fk + 220);
    const auto *fk_221 = buffer.data(fk + 221);
    const auto *fk_222 = buffer.data(fk + 222);
    const auto *fk_223 = buffer.data(fk + 223);
    const auto *fk_224 = buffer.data(fk + 224);
    const auto *fk_225 = buffer.data(fk + 225);
    const auto *fk_226 = buffer.data(fk + 226);
    const auto *fk_227 = buffer.data(fk + 227);
    const auto *fk_228 = buffer.data(fk + 228);
    const auto *fk_229 = buffer.data(fk + 229);
    const auto *fk_230 = buffer.data(fk + 230);
    const auto *fk_231 = buffer.data(fk + 231);
    const auto *fk_232 = buffer.data(fk + 232);
    const auto *fk_233 = buffer.data(fk + 233);
    const auto *fk_234 = buffer.data(fk + 234);
    const auto *fk_235 = buffer.data(fk + 235);
    const auto *fk_236 = buffer.data(fk + 236);
    const auto *fk_237 = buffer.data(fk + 237);
    const auto *fk_238 = buffer.data(fk + 238);
    const auto *fk_239 = buffer.data(fk + 239);
    const auto *fk_240 = buffer.data(fk + 240);
    const auto *fk_241 = buffer.data(fk + 241);
    const auto *fk_242 = buffer.data(fk + 242);
    const auto *fk_243 = buffer.data(fk + 243);
    const auto *fk_244 = buffer.data(fk + 244);
    const auto *fk_245 = buffer.data(fk + 245);
    const auto *fk_246 = buffer.data(fk + 246);
    const auto *fk_247 = buffer.data(fk + 247);
    const auto *fk_248 = buffer.data(fk + 248);
    const auto *fk_249 = buffer.data(fk + 249);
    const auto *fk_250 = buffer.data(fk + 250);
    const auto *fk_251 = buffer.data(fk + 251);
    const auto *fk_252 = buffer.data(fk + 252);
    const auto *fk_253 = buffer.data(fk + 253);
    const auto *fk_254 = buffer.data(fk + 254);
    const auto *fk_255 = buffer.data(fk + 255);
    const auto *fk_256 = buffer.data(fk + 256);
    const auto *fk_257 = buffer.data(fk + 257);
    const auto *fk_258 = buffer.data(fk + 258);
    const auto *fk_259 = buffer.data(fk + 259);
    const auto *fk_260 = buffer.data(fk + 260);
    const auto *fk_261 = buffer.data(fk + 261);
    const auto *fk_262 = buffer.data(fk + 262);
    const auto *fk_263 = buffer.data(fk + 263);
    const auto *fk_264 = buffer.data(fk + 264);
    const auto *fk_265 = buffer.data(fk + 265);
    const auto *fk_266 = buffer.data(fk + 266);
    const auto *fk_267 = buffer.data(fk + 267);
    const auto *fk_268 = buffer.data(fk + 268);
    const auto *fk_269 = buffer.data(fk + 269);
    const auto *fk_270 = buffer.data(fk + 270);
    const auto *fk_271 = buffer.data(fk + 271);
    const auto *fk_272 = buffer.data(fk + 272);
    const auto *fk_273 = buffer.data(fk + 273);
    const auto *fk_274 = buffer.data(fk + 274);
    const auto *fk_275 = buffer.data(fk + 275);
    const auto *fk_276 = buffer.data(fk + 276);
    const auto *fk_277 = buffer.data(fk + 277);
    const auto *fk_278 = buffer.data(fk + 278);
    const auto *fk_279 = buffer.data(fk + 279);
    const auto *fk_280 = buffer.data(fk + 280);
    const auto *fk_281 = buffer.data(fk + 281);
    const auto *fk_282 = buffer.data(fk + 282);
    const auto *fk_283 = buffer.data(fk + 283);
    const auto *fk_284 = buffer.data(fk + 284);
    const auto *fk_285 = buffer.data(fk + 285);
    const auto *fk_286 = buffer.data(fk + 286);
    const auto *fk_287 = buffer.data(fk + 287);
    const auto *fk_288 = buffer.data(fk + 288);
    const auto *fk_289 = buffer.data(fk + 289);
    const auto *fk_290 = buffer.data(fk + 290);
    const auto *fk_291 = buffer.data(fk + 291);
    const auto *fk_292 = buffer.data(fk + 292);
    const auto *fk_293 = buffer.data(fk + 293);
    const auto *fk_294 = buffer.data(fk + 294);
    const auto *fk_295 = buffer.data(fk + 295);
    const auto *fk_296 = buffer.data(fk + 296);
    const auto *fk_297 = buffer.data(fk + 297);
    const auto *fk_298 = buffer.data(fk + 298);
    const auto *fk_299 = buffer.data(fk + 299);
    const auto *fk_300 = buffer.data(fk + 300);
    const auto *fk_301 = buffer.data(fk + 301);
    const auto *fk_302 = buffer.data(fk + 302);
    const auto *fk_303 = buffer.data(fk + 303);
    const auto *fk_304 = buffer.data(fk + 304);
    const auto *fk_305 = buffer.data(fk + 305);
    const auto *fk_306 = buffer.data(fk + 306);
    const auto *fk_307 = buffer.data(fk + 307);
    const auto *fk_308 = buffer.data(fk + 308);
    const auto *fk_309 = buffer.data(fk + 309);
    const auto *fk_310 = buffer.data(fk + 310);
    const auto *fk_311 = buffer.data(fk + 311);
    const auto *fk_312 = buffer.data(fk + 312);
    const auto *fk_313 = buffer.data(fk + 313);
    const auto *fk_314 = buffer.data(fk + 314);
    const auto *fk_315 = buffer.data(fk + 315);
    const auto *fk_316 = buffer.data(fk + 316);
    const auto *fk_317 = buffer.data(fk + 317);
    const auto *fk_318 = buffer.data(fk + 318);
    const auto *fk_319 = buffer.data(fk + 319);
    const auto *fk_320 = buffer.data(fk + 320);
    const auto *fk_321 = buffer.data(fk + 321);
    const auto *fk_322 = buffer.data(fk + 322);
    const auto *fk_323 = buffer.data(fk + 323);
    const auto *fk_324 = buffer.data(fk + 324);
    const auto *fk_325 = buffer.data(fk + 325);
    const auto *fk_326 = buffer.data(fk + 326);
    const auto *fk_327 = buffer.data(fk + 327);
    const auto *fk_328 = buffer.data(fk + 328);
    const auto *fk_329 = buffer.data(fk + 329);
    const auto *fk_330 = buffer.data(fk + 330);
    const auto *fk_331 = buffer.data(fk + 331);
    const auto *fk_332 = buffer.data(fk + 332);
    const auto *fk_333 = buffer.data(fk + 333);
    const auto *fk_334 = buffer.data(fk + 334);
    const auto *fk_335 = buffer.data(fk + 335);
    const auto *fk_336 = buffer.data(fk + 336);
    const auto *fk_337 = buffer.data(fk + 337);
    const auto *fk_338 = buffer.data(fk + 338);
    const auto *fk_339 = buffer.data(fk + 339);
    const auto *fk_340 = buffer.data(fk + 340);
    const auto *fk_341 = buffer.data(fk + 341);
    const auto *fk_342 = buffer.data(fk + 342);
    const auto *fk_343 = buffer.data(fk + 343);
    const auto *fk_344 = buffer.data(fk + 344);
    const auto *fk_345 = buffer.data(fk + 345);
    const auto *fk_346 = buffer.data(fk + 346);
    const auto *fk_347 = buffer.data(fk + 347);
    const auto *fk_348 = buffer.data(fk + 348);
    const auto *fk_349 = buffer.data(fk + 349);
    const auto *fk_350 = buffer.data(fk + 350);
    const auto *fk_351 = buffer.data(fk + 351);
    const auto *fk_352 = buffer.data(fk + 352);
    const auto *fk_353 = buffer.data(fk + 353);
    const auto *fk_354 = buffer.data(fk + 354);
    const auto *fk_355 = buffer.data(fk + 355);
    const auto *fk_356 = buffer.data(fk + 356);
    const auto *fk_357 = buffer.data(fk + 357);
    const auto *fk_358 = buffer.data(fk + 358);
    const auto *fk_359 = buffer.data(fk + 359);

#pragma omp simd aligned(fk_37, fk_40, fk_42, fk_47, fk_51, fk_58, fk_64, fk_217, fk_220, \
                         fk_222, fk_227, fk_231, fk_238, fk_244 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * fk_37[k]
                 - f_1 * fk_42[k]
                 + f_2 * fk_51[k]
                 - f_3 * fk_64[k]
                 - f_4 * fk_217[k]
                 + f_5 * fk_222[k]
                 - f_0 * fk_231[k]
                 + f_6 * fk_244[k];

        g_1[k] = f_7 * fk_40[k]
                 - f_8 * fk_47[k]
                 + f_7 * fk_58[k]
                 - f_9 * fk_220[k]
                 + f_10 * fk_227[k]
                 - f_9 * fk_238[k];
    }

#pragma omp simd aligned(fk_37, fk_42, fk_44, fk_51, fk_53, fk_64, fk_66, fk_217, fk_222, \
                         fk_224, fk_231, fk_233, fk_244, fk_246 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_2[k] = -f_11 * fk_37[k]
                 + f_11 * fk_42[k]
                 + f_12 * fk_44[k]
                 + f_13 * fk_51[k]
                 - f_14 * fk_53[k]
                 - f_15 * fk_64[k]
                 + f_16 * fk_66[k]
                 + f_17 * fk_217[k]
                 - f_17 * fk_222[k]
                 - f_18 * fk_224[k]
                 - f_19 * fk_231[k]
                 + f_20 * fk_233[k]
                 + f_21 * fk_244[k]
                 - f_22 * fk_246[k];
    }

#pragma omp simd aligned(fk_40, fk_49, fk_58, fk_60, fk_220, fk_229, fk_238, \
                         fk_240 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = -f_23 * fk_40[k]
                 + f_24 * fk_49[k]
                 + f_23 * fk_58[k]
                 - f_24 * fk_60[k]
                 + f_25 * fk_220[k]
                 - f_26 * fk_229[k]
                 - f_25 * fk_238[k]
                 + f_26 * fk_240[k];
    }

#pragma omp simd aligned(fk_37, fk_42, fk_44, fk_51, fk_53, fk_55, fk_64, fk_66, fk_68, \
                         fk_217, fk_222, fk_224, fk_231, fk_233, fk_235, fk_244, fk_246, \
                         fk_248 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = f_27 * fk_37[k]
                 + f_28 * fk_42[k]
                 - f_29 * fk_44[k]
                 + f_30 * fk_51[k]
                 - f_31 * fk_53[k]
                 + f_32 * fk_55[k]
                 - f_30 * fk_64[k]
                 + f_33 * fk_66[k]
                 - f_34 * fk_68[k]
                 - f_30 * fk_217[k]
                 - f_35 * fk_222[k]
                 + f_33 * fk_224[k]
                 - f_36 * fk_231[k]
                 + f_37 * fk_233[k]
                 - f_34 * fk_235[k]
                 + f_36 * fk_244[k]
                 - f_38 * fk_246[k]
                 + f_39 * fk_248[k];
    }

#pragma omp simd aligned(fk_40, fk_47, fk_49, fk_58, fk_60, fk_62, fk_220, fk_227, fk_229, \
                         fk_238, fk_240, fk_242 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = f_40 * fk_40[k]
                 + f_41 * fk_47[k]
                 - f_42 * fk_49[k]
                 + f_40 * fk_58[k]
                 - f_42 * fk_60[k]
                 + f_43 * fk_62[k]
                 - f_44 * fk_220[k]
                 - f_45 * fk_227[k]
                 + f_46 * fk_229[k]
                 - f_44 * fk_238[k]
                 + f_46 * fk_240[k]
                 - f_47 * fk_242[k];
    }

#pragma omp simd aligned(fk_37, fk_42, fk_44, fk_51, fk_53, fk_55, fk_64, fk_66, fk_68, fk_70, \
                         fk_217, fk_222, fk_224, fk_231, fk_233, fk_235, fk_244, fk_246, \
                         fk_248, fk_250 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = -f_48 * fk_37[k]
                 - f_49 * fk_42[k]
                 + f_50 * fk_44[k]
                 - f_49 * fk_51[k]
                 + f_51 * fk_53[k]
                 - f_51 * fk_55[k]
                 - f_48 * fk_64[k]
                 + f_50 * fk_66[k]
                 - f_51 * fk_68[k]
                 + f_52 * fk_70[k]
                 + f_53 * fk_217[k]
                 + f_48 * fk_222[k]
                 - f_54 * fk_224[k]
                 + f_48 * fk_231[k]
                 - f_55 * fk_233[k]
                 + f_55 * fk_235[k]
                 + f_53 * fk_244[k]
                 - f_54 * fk_246[k]
                 + f_55 * fk_248[k]
                 - f_56 * fk_250[k];
    }

#pragma omp simd aligned(fk_38, fk_43, fk_45, fk_52, fk_54, fk_56, fk_65, fk_67, fk_69, fk_71, \
                         fk_218, fk_223, fk_225, fk_232, fk_234, fk_236, fk_245, fk_247, \
                         fk_249, fk_251 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = -f_57 * fk_38[k]
                 - f_58 * fk_43[k]
                 + f_59 * fk_45[k]
                 - f_58 * fk_52[k]
                 + f_60 * fk_54[k]
                 - f_61 * fk_56[k]
                 - f_57 * fk_65[k]
                 + f_59 * fk_67[k]
                 - f_61 * fk_69[k]
                 + f_62 * fk_71[k]
                 + f_63 * fk_218[k]
                 + f_57 * fk_223[k]
                 - f_64 * fk_225[k]
                 + f_57 * fk_232[k]
                 - f_65 * fk_234[k]
                 + f_66 * fk_236[k]
                 + f_63 * fk_245[k]
                 - f_64 * fk_247[k]
                 + f_66 * fk_249[k]
                 - f_67 * fk_251[k];
    }

#pragma omp simd aligned(fk_36, fk_39, fk_41, fk_46, fk_48, fk_50, fk_57, fk_59, fk_61, fk_63, \
                         fk_216, fk_219, fk_221, fk_226, fk_228, fk_230, fk_237, fk_239, \
                         fk_241, fk_243 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_8[k] = -f_48 * fk_36[k]
                 - f_49 * fk_39[k]
                 + f_50 * fk_41[k]
                 - f_49 * fk_46[k]
                 + f_51 * fk_48[k]
                 - f_51 * fk_50[k]
                 - f_48 * fk_57[k]
                 + f_50 * fk_59[k]
                 - f_51 * fk_61[k]
                 + f_52 * fk_63[k]
                 + f_53 * fk_216[k]
                 + f_48 * fk_219[k]
                 - f_54 * fk_221[k]
                 + f_48 * fk_226[k]
                 - f_55 * fk_228[k]
                 + f_55 * fk_230[k]
                 + f_53 * fk_237[k]
                 - f_54 * fk_239[k]
                 + f_55 * fk_241[k]
                 - f_56 * fk_243[k];
    }

#pragma omp simd aligned(fk_38, fk_43, fk_45, fk_52, fk_56, fk_65, fk_67, fk_69, fk_218, \
                         fk_223, fk_225, fk_232, fk_236, fk_245, fk_247, \
                         fk_249 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_9[k] = f_68 * fk_38[k]
                 + f_68 * fk_43[k]
                 - f_69 * fk_45[k]
                 - f_68 * fk_52[k]
                 + f_70 * fk_56[k]
                 - f_68 * fk_65[k]
                 + f_69 * fk_67[k]
                 - f_70 * fk_69[k]
                 - f_71 * fk_218[k]
                 - f_71 * fk_223[k]
                 + f_72 * fk_225[k]
                 + f_71 * fk_232[k]
                 - f_73 * fk_236[k]
                 + f_71 * fk_245[k]
                 - f_72 * fk_247[k]
                 + f_73 * fk_249[k];
    }

#pragma omp simd aligned(fk_36, fk_39, fk_41, fk_46, fk_48, fk_50, fk_57, fk_59, fk_61, \
                         fk_216, fk_219, fk_221, fk_226, fk_228, fk_230, fk_237, fk_239, \
                         fk_241 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = f_30 * fk_36[k]
                  - f_30 * fk_39[k]
                  - f_33 * fk_41[k]
                  - f_28 * fk_46[k]
                  + f_31 * fk_48[k]
                  + f_34 * fk_50[k]
                  - f_27 * fk_57[k]
                  + f_29 * fk_59[k]
                  - f_32 * fk_61[k]
                  - f_36 * fk_216[k]
                  + f_36 * fk_219[k]
                  + f_38 * fk_221[k]
                  + f_35 * fk_226[k]
                  - f_37 * fk_228[k]
                  - f_39 * fk_230[k]
                  + f_30 * fk_237[k]
                  - f_33 * fk_239[k]
                  + f_34 * fk_241[k];
    }

#pragma omp simd aligned(fk_38, fk_43, fk_45, fk_52, fk_54, fk_65, fk_67, fk_218, fk_223, \
                         fk_225, fk_232, fk_234, fk_245, fk_247 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_11[k] = -f_74 * fk_38[k]
                  + f_75 * fk_43[k]
                  + f_18 * fk_45[k]
                  + f_75 * fk_52[k]
                  - f_14 * fk_54[k]
                  - f_74 * fk_65[k]
                  + f_18 * fk_67[k]
                  + f_76 * fk_218[k]
                  - f_77 * fk_223[k]
                  - f_78 * fk_225[k]
                  - f_77 * fk_232[k]
                  + f_20 * fk_234[k]
                  + f_76 * fk_245[k]
                  - f_78 * fk_247[k];
    }

#pragma omp simd aligned(fk_36, fk_39, fk_41, fk_46, fk_48, fk_57, fk_59, fk_216, fk_219, \
                         fk_221, fk_226, fk_228, fk_237, fk_239 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_12[k] = -f_15 * fk_36[k]
                  + f_13 * fk_39[k]
                  + f_16 * fk_41[k]
                  + f_11 * fk_46[k]
                  - f_14 * fk_48[k]
                  - f_11 * fk_57[k]
                  + f_12 * fk_59[k]
                  + f_21 * fk_216[k]
                  - f_19 * fk_219[k]
                  - f_22 * fk_221[k]
                  - f_17 * fk_226[k]
                  + f_20 * fk_228[k]
                  + f_17 * fk_237[k]
                  - f_18 * fk_239[k];
    }

#pragma omp simd aligned(fk_38, fk_43, fk_52, fk_65, fk_218, fk_223, fk_232, \
                         fk_245 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = f_79 * fk_38[k]
                  - f_80 * fk_43[k]
                  + f_80 * fk_52[k]
                  - f_79 * fk_65[k]
                  - f_81 * fk_218[k]
                  + f_82 * fk_223[k]
                  - f_82 * fk_232[k]
                  + f_81 * fk_245[k];
    }

#pragma omp simd aligned(fk_36, fk_39, fk_46, fk_57, fk_145, fk_150, fk_159, fk_172, fk_216, \
                         fk_219, fk_226, fk_237 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = f_3 * fk_36[k]
                  - f_2 * fk_39[k]
                  + f_1 * fk_46[k]
                  - f_0 * fk_57[k]
                  - f_6 * fk_216[k]
                  + f_0 * fk_219[k]
                  - f_5 * fk_226[k]
                  + f_4 * fk_237[k];

        g_15[k] = f_83 * fk_145[k]
                  - f_84 * fk_150[k]
                  + f_85 * fk_159[k]
                  - f_86 * fk_172[k];
    }

#pragma omp simd aligned(fk_145, fk_148, fk_150, fk_152, fk_155, fk_157, fk_159, fk_161, \
                         fk_166, fk_168, fk_172, fk_174 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_16[k] = f_87 * fk_148[k]
                  - f_88 * fk_155[k]
                  + f_87 * fk_166[k];

        g_17[k] = -f_89 * fk_145[k]
                  + f_89 * fk_150[k]
                  + f_90 * fk_152[k]
                  + f_91 * fk_159[k]
                  - f_92 * fk_161[k]
                  - f_93 * fk_172[k]
                  + f_94 * fk_174[k];

        g_18[k] = -f_95 * fk_148[k]
                  + f_96 * fk_157[k]
                  + f_95 * fk_166[k]
                  - f_96 * fk_168[k];
    }

#pragma omp simd aligned(fk_145, fk_150, fk_152, fk_159, fk_161, fk_163, fk_172, fk_174, \
                         fk_176 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = f_97 * fk_145[k]
                  + f_98 * fk_150[k]
                  - f_99 * fk_152[k]
                  + f_100 * fk_159[k]
                  - f_101 * fk_161[k]
                  + f_102 * fk_163[k]
                  - f_100 * fk_172[k]
                  + f_103 * fk_174[k]
                  - f_104 * fk_176[k];
    }

#pragma omp simd aligned(fk_148, fk_155, fk_157, fk_166, fk_168, \
                         fk_170 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = f_50 * fk_148[k]
                  + f_51 * fk_155[k]
                  - f_105 * fk_157[k]
                  + f_50 * fk_166[k]
                  - f_105 * fk_168[k]
                  + f_106 * fk_170[k];
    }

#pragma omp simd aligned(fk_145, fk_150, fk_152, fk_159, fk_161, fk_163, fk_172, fk_174, \
                         fk_176, fk_178 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = -f_107 * fk_145[k]
                  - f_44 * fk_150[k]
                  + f_69 * fk_152[k]
                  - f_44 * fk_159[k]
                  + f_42 * fk_161[k]
                  - f_42 * fk_163[k]
                  - f_107 * fk_172[k]
                  + f_69 * fk_174[k]
                  - f_42 * fk_176[k]
                  + f_108 * fk_178[k];
    }

#pragma omp simd aligned(fk_146, fk_151, fk_153, fk_160, fk_162, fk_164, fk_173, fk_175, \
                         fk_177, fk_179 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = -f_109 * fk_146[k]
                  - f_110 * fk_151[k]
                  + f_111 * fk_153[k]
                  - f_110 * fk_160[k]
                  + f_112 * fk_162[k]
                  - f_113 * fk_164[k]
                  - f_109 * fk_173[k]
                  + f_111 * fk_175[k]
                  - f_113 * fk_177[k]
                  + f_114 * fk_179[k];
    }

#pragma omp simd aligned(fk_144, fk_147, fk_149, fk_154, fk_156, fk_158, fk_165, fk_167, \
                         fk_169, fk_171 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = -f_107 * fk_144[k]
                  - f_44 * fk_147[k]
                  + f_69 * fk_149[k]
                  - f_44 * fk_154[k]
                  + f_42 * fk_156[k]
                  - f_42 * fk_158[k]
                  - f_107 * fk_165[k]
                  + f_69 * fk_167[k]
                  - f_42 * fk_169[k]
                  + f_108 * fk_171[k];
    }

#pragma omp simd aligned(fk_146, fk_151, fk_153, fk_160, fk_164, fk_173, fk_175, \
                         fk_177 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = f_115 * fk_146[k]
                  + f_115 * fk_151[k]
                  - f_116 * fk_153[k]
                  - f_115 * fk_160[k]
                  + f_117 * fk_164[k]
                  - f_115 * fk_173[k]
                  + f_116 * fk_175[k]
                  - f_117 * fk_177[k];
    }

#pragma omp simd aligned(fk_144, fk_147, fk_149, fk_154, fk_156, fk_158, fk_165, fk_167, \
                         fk_169 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = f_100 * fk_144[k]
                  - f_100 * fk_147[k]
                  - f_103 * fk_149[k]
                  - f_98 * fk_154[k]
                  + f_101 * fk_156[k]
                  + f_104 * fk_158[k]
                  - f_97 * fk_165[k]
                  + f_99 * fk_167[k]
                  - f_102 * fk_169[k];
    }

#pragma omp simd aligned(fk_144, fk_146, fk_147, fk_149, fk_151, fk_153, fk_154, fk_156, \
                         fk_160, fk_162, fk_165, fk_167, fk_173, \
                         fk_175 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = -f_118 * fk_146[k]
                  + f_119 * fk_151[k]
                  + f_120 * fk_153[k]
                  + f_119 * fk_160[k]
                  - f_92 * fk_162[k]
                  - f_118 * fk_173[k]
                  + f_120 * fk_175[k];

        g_27[k] = -f_93 * fk_144[k]
                  + f_91 * fk_147[k]
                  + f_94 * fk_149[k]
                  + f_89 * fk_154[k]
                  - f_92 * fk_156[k]
                  - f_89 * fk_165[k]
                  + f_90 * fk_167[k];
    }

#pragma omp simd aligned(fk_144, fk_146, fk_147, fk_151, fk_154, fk_160, fk_165, \
                         fk_173 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = f_121 * fk_146[k]
                  - f_122 * fk_151[k]
                  + f_122 * fk_160[k]
                  - f_121 * fk_173[k];

        g_29[k] = f_86 * fk_144[k]
                  - f_85 * fk_147[k]
                  + f_84 * fk_154[k]
                  - f_83 * fk_165[k];
    }

#pragma omp simd aligned(fk_37, fk_42, fk_51, fk_64, fk_217, fk_222, fk_231, fk_244, fk_289, \
                         fk_294, fk_303, fk_316 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = -f_123 * fk_37[k]
                  + f_124 * fk_42[k]
                  - f_125 * fk_51[k]
                  + f_126 * fk_64[k]
                  - f_123 * fk_217[k]
                  + f_124 * fk_222[k]
                  - f_125 * fk_231[k]
                  + f_126 * fk_244[k]
                  + f_127 * fk_289[k]
                  - f_128 * fk_294[k]
                  + f_129 * fk_303[k]
                  - f_130 * fk_316[k];
    }

#pragma omp simd aligned(fk_40, fk_47, fk_58, fk_220, fk_227, fk_238, fk_292, fk_299, \
                         fk_310 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = -f_131 * fk_40[k]
                  + f_132 * fk_47[k]
                  - f_131 * fk_58[k]
                  - f_131 * fk_220[k]
                  + f_132 * fk_227[k]
                  - f_131 * fk_238[k]
                  + f_133 * fk_292[k]
                  - f_134 * fk_299[k]
                  + f_133 * fk_310[k];
    }

#pragma omp simd aligned(fk_37, fk_42, fk_44, fk_51, fk_53, fk_64, fk_66, fk_217, fk_222, \
                         fk_224, fk_231, fk_233, fk_244, fk_246, fk_289, fk_294, fk_296, \
                         fk_303, fk_305, fk_316, fk_318 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = f_135 * fk_37[k]
                  - f_135 * fk_42[k]
                  - f_136 * fk_44[k]
                  - f_137 * fk_51[k]
                  + f_138 * fk_53[k]
                  + f_139 * fk_64[k]
                  - f_140 * fk_66[k]
                  + f_135 * fk_217[k]
                  - f_135 * fk_222[k]
                  - f_136 * fk_224[k]
                  - f_137 * fk_231[k]
                  + f_138 * fk_233[k]
                  + f_139 * fk_244[k]
                  - f_140 * fk_246[k]
                  - f_141 * fk_289[k]
                  + f_141 * fk_294[k]
                  + f_142 * fk_296[k]
                  + f_143 * fk_303[k]
                  - f_144 * fk_305[k]
                  - f_145 * fk_316[k]
                  + f_146 * fk_318[k];
    }

#pragma omp simd aligned(fk_40, fk_49, fk_58, fk_60, fk_220, fk_229, fk_238, fk_240, fk_292, \
                         fk_301, fk_310, fk_312 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = f_147 * fk_40[k]
                  - f_148 * fk_49[k]
                  - f_147 * fk_58[k]
                  + f_148 * fk_60[k]
                  + f_147 * fk_220[k]
                  - f_148 * fk_229[k]
                  - f_147 * fk_238[k]
                  + f_148 * fk_240[k]
                  - f_149 * fk_292[k]
                  + f_150 * fk_301[k]
                  + f_149 * fk_310[k]
                  - f_150 * fk_312[k];
    }

#pragma omp simd aligned(fk_37, fk_42, fk_44, fk_51, fk_53, fk_55, fk_64, fk_66, fk_68, \
                         fk_217, fk_222, fk_224, fk_231, fk_233, fk_235, fk_244, fk_246, \
                         fk_248, fk_289, fk_294, fk_296, fk_303, fk_305, fk_307, fk_316, \
                         fk_318, fk_320 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = -f_151 * fk_37[k]
                  - f_152 * fk_42[k]
                  + f_153 * fk_44[k]
                  - f_154 * fk_51[k]
                  + f_155 * fk_53[k]
                  - f_156 * fk_55[k]
                  + f_154 * fk_64[k]
                  - f_157 * fk_66[k]
                  + f_158 * fk_68[k]
                  - f_151 * fk_217[k]
                  - f_152 * fk_222[k]
                  + f_153 * fk_224[k]
                  - f_154 * fk_231[k]
                  + f_155 * fk_233[k]
                  - f_156 * fk_235[k]
                  + f_154 * fk_244[k]
                  - f_157 * fk_246[k]
                  + f_158 * fk_248[k]
                  + f_159 * fk_289[k]
                  + f_157 * fk_294[k]
                  - f_160 * fk_296[k]
                  + f_161 * fk_303[k]
                  - f_162 * fk_305[k]
                  + f_163 * fk_307[k]
                  - f_161 * fk_316[k]
                  + f_156 * fk_318[k]
                  - f_164 * fk_320[k];
    }

#pragma omp simd aligned(fk_40, fk_47, fk_49, fk_58, fk_60, fk_62, fk_220, fk_227, fk_229, \
                         fk_238, fk_240, fk_242, fk_292, fk_299, fk_301, fk_310, fk_312, \
                         fk_314 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = -f_165 * fk_40[k]
                  - f_166 * fk_47[k]
                  + f_167 * fk_49[k]
                  - f_165 * fk_58[k]
                  + f_167 * fk_60[k]
                  - f_168 * fk_62[k]
                  - f_165 * fk_220[k]
                  - f_166 * fk_227[k]
                  + f_167 * fk_229[k]
                  - f_165 * fk_238[k]
                  + f_167 * fk_240[k]
                  - f_168 * fk_242[k]
                  + f_169 * fk_292[k]
                  + f_170 * fk_299[k]
                  - f_171 * fk_301[k]
                  + f_169 * fk_310[k]
                  - f_171 * fk_312[k]
                  + f_172 * fk_314[k];
    }

#pragma omp simd aligned(fk_37, fk_42, fk_44, fk_51, fk_53, fk_55, fk_64, fk_66, fk_68, fk_70, \
                         fk_217, fk_222, fk_224, fk_231, fk_233, fk_235, fk_244, fk_246, \
                         fk_248, fk_250, fk_289, fk_294, fk_296, fk_303, fk_305, fk_307, \
                         fk_316, fk_318, fk_320, fk_322 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_36[k] = f_173 * fk_37[k]
                  + f_174 * fk_42[k]
                  - f_175 * fk_44[k]
                  + f_174 * fk_51[k]
                  - f_176 * fk_53[k]
                  + f_176 * fk_55[k]
                  + f_173 * fk_64[k]
                  - f_175 * fk_66[k]
                  + f_176 * fk_68[k]
                  - f_177 * fk_70[k]
                  + f_173 * fk_217[k]
                  + f_174 * fk_222[k]
                  - f_175 * fk_224[k]
                  + f_174 * fk_231[k]
                  - f_176 * fk_233[k]
                  + f_176 * fk_235[k]
                  + f_173 * fk_244[k]
                  - f_175 * fk_246[k]
                  + f_176 * fk_248[k]
                  - f_177 * fk_250[k]
                  - f_178 * fk_289[k]
                  - f_179 * fk_294[k]
                  + f_180 * fk_296[k]
                  - f_179 * fk_303[k]
                  + f_181 * fk_305[k]
                  - f_181 * fk_307[k]
                  - f_178 * fk_316[k]
                  + f_180 * fk_318[k]
                  - f_181 * fk_320[k]
                  + f_182 * fk_322[k];
    }

#pragma omp simd aligned(fk_38, fk_43, fk_45, fk_52, fk_54, fk_56, fk_65, fk_67, fk_69, fk_71, \
                         fk_218, fk_223, fk_225, fk_232, fk_234, fk_236, fk_245, fk_247, \
                         fk_249, fk_251, fk_290, fk_295, fk_297, fk_304, fk_306, fk_308, \
                         fk_317, fk_319, fk_321, fk_323 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_37[k] = f_183 * fk_38[k]
                  + f_184 * fk_43[k]
                  - f_185 * fk_45[k]
                  + f_184 * fk_52[k]
                  - f_186 * fk_54[k]
                  + f_187 * fk_56[k]
                  + f_183 * fk_65[k]
                  - f_185 * fk_67[k]
                  + f_187 * fk_69[k]
                  - f_188 * fk_71[k]
                  + f_183 * fk_218[k]
                  + f_184 * fk_223[k]
                  - f_185 * fk_225[k]
                  + f_184 * fk_232[k]
                  - f_186 * fk_234[k]
                  + f_187 * fk_236[k]
                  + f_183 * fk_245[k]
                  - f_185 * fk_247[k]
                  + f_187 * fk_249[k]
                  - f_188 * fk_251[k]
                  - f_189 * fk_290[k]
                  - f_186 * fk_295[k]
                  + f_190 * fk_297[k]
                  - f_186 * fk_304[k]
                  + f_191 * fk_306[k]
                  - f_192 * fk_308[k]
                  - f_189 * fk_317[k]
                  + f_190 * fk_319[k]
                  - f_192 * fk_321[k]
                  + f_193 * fk_323[k];
    }

#pragma omp simd aligned(fk_36, fk_39, fk_41, fk_46, fk_48, fk_50, fk_57, fk_59, fk_61, fk_63, \
                         fk_216, fk_219, fk_221, fk_226, fk_228, fk_230, fk_237, fk_239, \
                         fk_241, fk_243, fk_288, fk_291, fk_293, fk_298, fk_300, fk_302, \
                         fk_309, fk_311, fk_313, fk_315 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = f_173 * fk_36[k]
                  + f_174 * fk_39[k]
                  - f_175 * fk_41[k]
                  + f_174 * fk_46[k]
                  - f_176 * fk_48[k]
                  + f_176 * fk_50[k]
                  + f_173 * fk_57[k]
                  - f_175 * fk_59[k]
                  + f_176 * fk_61[k]
                  - f_177 * fk_63[k]
                  + f_173 * fk_216[k]
                  + f_174 * fk_219[k]
                  - f_175 * fk_221[k]
                  + f_174 * fk_226[k]
                  - f_176 * fk_228[k]
                  + f_176 * fk_230[k]
                  + f_173 * fk_237[k]
                  - f_175 * fk_239[k]
                  + f_176 * fk_241[k]
                  - f_177 * fk_243[k]
                  - f_178 * fk_288[k]
                  - f_179 * fk_291[k]
                  + f_180 * fk_293[k]
                  - f_179 * fk_298[k]
                  + f_181 * fk_300[k]
                  - f_181 * fk_302[k]
                  - f_178 * fk_309[k]
                  + f_180 * fk_311[k]
                  - f_181 * fk_313[k]
                  + f_182 * fk_315[k];
    }

#pragma omp simd aligned(fk_38, fk_43, fk_45, fk_52, fk_56, fk_65, fk_67, fk_69, fk_218, \
                         fk_223, fk_225, fk_232, fk_236, fk_245, fk_247, fk_249, fk_290, \
                         fk_295, fk_297, fk_304, fk_308, fk_317, fk_319, \
                         fk_321 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_39[k] = -f_194 * fk_38[k]
                  - f_194 * fk_43[k]
                  + f_195 * fk_45[k]
                  + f_194 * fk_52[k]
                  - f_196 * fk_56[k]
                  + f_194 * fk_65[k]
                  - f_195 * fk_67[k]
                  + f_196 * fk_69[k]
                  - f_194 * fk_218[k]
                  - f_194 * fk_223[k]
                  + f_195 * fk_225[k]
                  + f_194 * fk_232[k]
                  - f_196 * fk_236[k]
                  + f_194 * fk_245[k]
                  - f_195 * fk_247[k]
                  + f_196 * fk_249[k]
                  + f_166 * fk_290[k]
                  + f_166 * fk_295[k]
                  - f_197 * fk_297[k]
                  - f_166 * fk_304[k]
                  + f_198 * fk_308[k]
                  - f_166 * fk_317[k]
                  + f_197 * fk_319[k]
                  - f_198 * fk_321[k];
    }

#pragma omp simd aligned(fk_36, fk_39, fk_41, fk_46, fk_48, fk_50, fk_57, fk_59, fk_61, \
                         fk_216, fk_219, fk_221, fk_226, fk_228, fk_230, fk_237, fk_239, \
                         fk_241, fk_288, fk_291, fk_293, fk_298, fk_300, fk_302, fk_309, \
                         fk_311, fk_313 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_40[k] = -f_154 * fk_36[k]
                  + f_154 * fk_39[k]
                  + f_157 * fk_41[k]
                  + f_152 * fk_46[k]
                  - f_155 * fk_48[k]
                  - f_158 * fk_50[k]
                  + f_151 * fk_57[k]
                  - f_153 * fk_59[k]
                  + f_156 * fk_61[k]
                  - f_154 * fk_216[k]
                  + f_154 * fk_219[k]
                  + f_157 * fk_221[k]
                  + f_152 * fk_226[k]
                  - f_155 * fk_228[k]
                  - f_158 * fk_230[k]
                  + f_151 * fk_237[k]
                  - f_153 * fk_239[k]
                  + f_156 * fk_241[k]
                  + f_161 * fk_288[k]
                  - f_161 * fk_291[k]
                  - f_156 * fk_293[k]
                  - f_157 * fk_298[k]
                  + f_162 * fk_300[k]
                  + f_164 * fk_302[k]
                  - f_159 * fk_309[k]
                  + f_160 * fk_311[k]
                  - f_163 * fk_313[k];
    }

#pragma omp simd aligned(fk_38, fk_43, fk_45, fk_52, fk_54, fk_65, fk_67, fk_218, fk_223, \
                         fk_225, fk_232, fk_234, fk_245, fk_247, fk_290, fk_295, fk_297, \
                         fk_304, fk_306, fk_317, fk_319 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_41[k] = f_199 * fk_38[k]
                  - f_200 * fk_43[k]
                  - f_141 * fk_45[k]
                  - f_200 * fk_52[k]
                  + f_138 * fk_54[k]
                  + f_199 * fk_65[k]
                  - f_141 * fk_67[k]
                  + f_199 * fk_218[k]
                  - f_200 * fk_223[k]
                  - f_141 * fk_225[k]
                  - f_200 * fk_232[k]
                  + f_138 * fk_234[k]
                  + f_199 * fk_245[k]
                  - f_141 * fk_247[k]
                  - f_147 * fk_290[k]
                  + f_138 * fk_295[k]
                  + f_148 * fk_297[k]
                  + f_138 * fk_304[k]
                  - f_144 * fk_306[k]
                  - f_147 * fk_317[k]
                  + f_148 * fk_319[k];
    }

#pragma omp simd aligned(fk_36, fk_39, fk_41, fk_46, fk_48, fk_57, fk_59, fk_216, fk_219, \
                         fk_221, fk_226, fk_228, fk_237, fk_239, fk_288, fk_291, fk_293, \
                         fk_298, fk_300, fk_309, fk_311 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = f_139 * fk_36[k]
                  - f_137 * fk_39[k]
                  - f_140 * fk_41[k]
                  - f_135 * fk_46[k]
                  + f_138 * fk_48[k]
                  + f_135 * fk_57[k]
                  - f_136 * fk_59[k]
                  + f_139 * fk_216[k]
                  - f_137 * fk_219[k]
                  - f_140 * fk_221[k]
                  - f_135 * fk_226[k]
                  + f_138 * fk_228[k]
                  + f_135 * fk_237[k]
                  - f_136 * fk_239[k]
                  - f_145 * fk_288[k]
                  + f_143 * fk_291[k]
                  + f_146 * fk_293[k]
                  + f_141 * fk_298[k]
                  - f_144 * fk_300[k]
                  - f_141 * fk_309[k]
                  + f_142 * fk_311[k];
    }

#pragma omp simd aligned(fk_38, fk_43, fk_52, fk_65, fk_218, fk_223, fk_232, fk_245, fk_290, \
                         fk_295, fk_304, fk_317 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_43[k] = -f_201 * fk_38[k]
                  + f_202 * fk_43[k]
                  - f_202 * fk_52[k]
                  + f_201 * fk_65[k]
                  - f_201 * fk_218[k]
                  + f_202 * fk_223[k]
                  - f_202 * fk_232[k]
                  + f_201 * fk_245[k]
                  + f_203 * fk_290[k]
                  - f_204 * fk_295[k]
                  + f_204 * fk_304[k]
                  - f_203 * fk_317[k];
    }

#pragma omp simd aligned(fk_36, fk_39, fk_46, fk_57, fk_216, fk_219, fk_226, fk_237, fk_288, \
                         fk_291, fk_298, fk_309 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_44[k] = -f_126 * fk_36[k]
                  + f_125 * fk_39[k]
                  - f_124 * fk_46[k]
                  + f_123 * fk_57[k]
                  - f_126 * fk_216[k]
                  + f_125 * fk_219[k]
                  - f_124 * fk_226[k]
                  + f_123 * fk_237[k]
                  + f_130 * fk_288[k]
                  - f_129 * fk_291[k]
                  + f_128 * fk_298[k]
                  - f_127 * fk_309[k];
    }

#pragma omp simd aligned(fk_73, fk_78, fk_87, fk_100, fk_253, fk_258, fk_267, fk_280, fk_325, \
                         fk_330, fk_339, fk_352 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_45[k] = -f_205 * fk_73[k]
                  + f_206 * fk_78[k]
                  - f_207 * fk_87[k]
                  + f_208 * fk_100[k]
                  - f_205 * fk_253[k]
                  + f_206 * fk_258[k]
                  - f_207 * fk_267[k]
                  + f_208 * fk_280[k]
                  + f_209 * fk_325[k]
                  - f_210 * fk_330[k]
                  + f_211 * fk_339[k]
                  - f_212 * fk_352[k];
    }

#pragma omp simd aligned(fk_76, fk_83, fk_94, fk_256, fk_263, fk_274, fk_328, fk_335, \
                         fk_346 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_46[k] = -f_213 * fk_76[k]
                  + f_214 * fk_83[k]
                  - f_213 * fk_94[k]
                  - f_213 * fk_256[k]
                  + f_214 * fk_263[k]
                  - f_213 * fk_274[k]
                  + f_215 * fk_328[k]
                  - f_216 * fk_335[k]
                  + f_215 * fk_346[k];
    }

#pragma omp simd aligned(fk_73, fk_78, fk_80, fk_87, fk_89, fk_100, fk_102, fk_253, fk_258, \
                         fk_260, fk_267, fk_269, fk_280, fk_282, fk_325, fk_330, fk_332, \
                         fk_339, fk_341, fk_352, fk_354 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_47[k] = f_217 * fk_73[k]
                  - f_217 * fk_78[k]
                  - f_218 * fk_80[k]
                  - f_219 * fk_87[k]
                  + f_220 * fk_89[k]
                  + f_221 * fk_100[k]
                  - f_222 * fk_102[k]
                  + f_217 * fk_253[k]
                  - f_217 * fk_258[k]
                  - f_218 * fk_260[k]
                  - f_219 * fk_267[k]
                  + f_220 * fk_269[k]
                  + f_221 * fk_280[k]
                  - f_222 * fk_282[k]
                  - f_223 * fk_325[k]
                  + f_223 * fk_330[k]
                  + f_224 * fk_332[k]
                  + f_225 * fk_339[k]
                  - f_226 * fk_341[k]
                  - f_227 * fk_352[k]
                  + f_228 * fk_354[k];
    }

#pragma omp simd aligned(fk_76, fk_85, fk_94, fk_96, fk_256, fk_265, fk_274, fk_276, fk_328, \
                         fk_337, fk_346, fk_348 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_48[k] = f_229 * fk_76[k]
                  - f_226 * fk_85[k]
                  - f_229 * fk_94[k]
                  + f_226 * fk_96[k]
                  + f_229 * fk_256[k]
                  - f_226 * fk_265[k]
                  - f_229 * fk_274[k]
                  + f_226 * fk_276[k]
                  - f_230 * fk_328[k]
                  + f_231 * fk_337[k]
                  + f_230 * fk_346[k]
                  - f_231 * fk_348[k];
    }

#pragma omp simd aligned(fk_73, fk_78, fk_80, fk_87, fk_89, fk_91, fk_100, fk_102, fk_104, \
                         fk_253, fk_258, fk_260, fk_267, fk_269, fk_271, fk_280, fk_282, \
                         fk_284, fk_325, fk_330, fk_332, fk_339, fk_341, fk_343, fk_352, \
                         fk_354, fk_356 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_49[k] = -f_232 * fk_73[k]
                  - f_233 * fk_78[k]
                  + f_234 * fk_80[k]
                  - f_235 * fk_87[k]
                  + f_236 * fk_89[k]
                  - f_237 * fk_91[k]
                  + f_235 * fk_100[k]
                  - f_238 * fk_102[k]
                  + f_239 * fk_104[k]
                  - f_232 * fk_253[k]
                  - f_233 * fk_258[k]
                  + f_234 * fk_260[k]
                  - f_235 * fk_267[k]
                  + f_236 * fk_269[k]
                  - f_237 * fk_271[k]
                  + f_235 * fk_280[k]
                  - f_238 * fk_282[k]
                  + f_239 * fk_284[k]
                  + f_240 * fk_325[k]
                  + f_241 * fk_330[k]
                  - f_236 * fk_332[k]
                  + f_242 * fk_339[k]
                  - f_239 * fk_341[k]
                  + f_243 * fk_343[k]
                  - f_242 * fk_352[k]
                  + f_244 * fk_354[k]
                  - f_245 * fk_356[k];
    }

#pragma omp simd aligned(fk_76, fk_83, fk_85, fk_94, fk_96, fk_98, fk_256, fk_263, fk_265, \
                         fk_274, fk_276, fk_278, fk_328, fk_335, fk_337, fk_346, fk_348, \
                         fk_350 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_50[k] = -f_246 * fk_76[k]
                  - f_247 * fk_83[k]
                  + f_181 * fk_85[k]
                  - f_246 * fk_94[k]
                  + f_181 * fk_96[k]
                  - f_248 * fk_98[k]
                  - f_246 * fk_256[k]
                  - f_247 * fk_263[k]
                  + f_181 * fk_265[k]
                  - f_246 * fk_274[k]
                  + f_181 * fk_276[k]
                  - f_248 * fk_278[k]
                  + f_175 * fk_328[k]
                  + f_176 * fk_335[k]
                  - f_249 * fk_337[k]
                  + f_175 * fk_346[k]
                  - f_249 * fk_348[k]
                  + f_250 * fk_350[k];
    }

#pragma omp simd aligned(fk_73, fk_78, fk_80, fk_87, fk_89, fk_91, fk_100, fk_102, fk_104, \
                         fk_106, fk_253, fk_258, fk_260, fk_267, fk_269, fk_271, fk_280, \
                         fk_282, fk_284, fk_286, fk_325, fk_330, fk_332, fk_339, fk_341, \
                         fk_343, fk_352, fk_354, fk_356, fk_358 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_51[k] = f_251 * fk_73[k]
                  + f_194 * fk_78[k]
                  - f_169 * fk_80[k]
                  + f_194 * fk_87[k]
                  - f_170 * fk_89[k]
                  + f_170 * fk_91[k]
                  + f_251 * fk_100[k]
                  - f_169 * fk_102[k]
                  + f_170 * fk_104[k]
                  - f_252 * fk_106[k]
                  + f_251 * fk_253[k]
                  + f_194 * fk_258[k]
                  - f_169 * fk_260[k]
                  + f_194 * fk_267[k]
                  - f_170 * fk_269[k]
                  + f_170 * fk_271[k]
                  + f_251 * fk_280[k]
                  - f_169 * fk_282[k]
                  + f_170 * fk_284[k]
                  - f_252 * fk_286[k]
                  - f_253 * fk_325[k]
                  - f_254 * fk_330[k]
                  + f_195 * fk_332[k]
                  - f_254 * fk_339[k]
                  + f_167 * fk_341[k]
                  - f_167 * fk_343[k]
                  - f_253 * fk_352[k]
                  + f_195 * fk_354[k]
                  - f_167 * fk_356[k]
                  + f_255 * fk_358[k];
    }

#pragma omp simd aligned(fk_74, fk_79, fk_81, fk_88, fk_90, fk_92, fk_101, fk_103, fk_105, \
                         fk_107, fk_254, fk_259, fk_261, fk_268, fk_270, fk_272, fk_281, \
                         fk_283, fk_285, fk_287, fk_326, fk_331, fk_333, fk_340, fk_342, \
                         fk_344, fk_353, fk_355, fk_357, fk_359 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_52[k] = 3.28125 * fk_74[k]
                  + 9.84375 * fk_79[k]
                  - 19.6875 * fk_81[k]
                  + 9.84375 * fk_88[k]
                  - 39.375 * fk_90[k]
                  + 15.75 * fk_92[k]
                  + 3.28125 * fk_101[k]
                  - 19.6875 * fk_103[k]
                  + 15.75 * fk_105[k]
                  - 1.5 * fk_107[k]
                  + 3.28125 * fk_254[k]
                  + 9.84375 * fk_259[k]
                  - 19.6875 * fk_261[k]
                  + 9.84375 * fk_268[k]
                  - 39.375 * fk_270[k]
                  + 15.75 * fk_272[k]
                  + 3.28125 * fk_281[k]
                  - 19.6875 * fk_283[k]
                  + 15.75 * fk_285[k]
                  - 1.5 * fk_287[k]
                  - 2.1875 * fk_326[k]
                  - 6.5625 * fk_331[k]
                  + 13.125 * fk_333[k]
                  - 6.5625 * fk_340[k]
                  + 26.25 * fk_342[k]
                  - 10.5 * fk_344[k]
                  - 2.1875 * fk_353[k]
                  + 13.125 * fk_355[k]
                  - 10.5 * fk_357[k]
                  + fk_359[k];
    }

#pragma omp simd aligned(fk_72, fk_75, fk_77, fk_82, fk_84, fk_86, fk_93, fk_95, fk_97, fk_99, \
                         fk_252, fk_255, fk_257, fk_262, fk_264, fk_266, fk_273, fk_275, \
                         fk_277, fk_279, fk_324, fk_327, fk_329, fk_334, fk_336, fk_338, \
                         fk_345, fk_347, fk_349, fk_351 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_53[k] = f_251 * fk_72[k]
                  + f_194 * fk_75[k]
                  - f_169 * fk_77[k]
                  + f_194 * fk_82[k]
                  - f_170 * fk_84[k]
                  + f_170 * fk_86[k]
                  + f_251 * fk_93[k]
                  - f_169 * fk_95[k]
                  + f_170 * fk_97[k]
                  - f_252 * fk_99[k]
                  + f_251 * fk_252[k]
                  + f_194 * fk_255[k]
                  - f_169 * fk_257[k]
                  + f_194 * fk_262[k]
                  - f_170 * fk_264[k]
                  + f_170 * fk_266[k]
                  + f_251 * fk_273[k]
                  - f_169 * fk_275[k]
                  + f_170 * fk_277[k]
                  - f_252 * fk_279[k]
                  - f_253 * fk_324[k]
                  - f_254 * fk_327[k]
                  + f_195 * fk_329[k]
                  - f_254 * fk_334[k]
                  + f_167 * fk_336[k]
                  - f_167 * fk_338[k]
                  - f_253 * fk_345[k]
                  + f_195 * fk_347[k]
                  - f_167 * fk_349[k]
                  + f_255 * fk_351[k];
    }

#pragma omp simd aligned(fk_74, fk_79, fk_81, fk_88, fk_92, fk_101, fk_103, fk_105, fk_254, \
                         fk_259, fk_261, fk_268, fk_272, fk_281, fk_283, fk_285, fk_326, \
                         fk_331, fk_333, fk_340, fk_344, fk_353, fk_355, \
                         fk_357 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_54[k] = -f_256 * fk_74[k]
                  - f_256 * fk_79[k]
                  + f_180 * fk_81[k]
                  + f_256 * fk_88[k]
                  - f_257 * fk_92[k]
                  + f_256 * fk_101[k]
                  - f_180 * fk_103[k]
                  + f_257 * fk_105[k]
                  - f_256 * fk_254[k]
                  - f_256 * fk_259[k]
                  + f_180 * fk_261[k]
                  + f_256 * fk_268[k]
                  - f_257 * fk_272[k]
                  + f_256 * fk_281[k]
                  - f_180 * fk_283[k]
                  + f_257 * fk_285[k]
                  + f_179 * fk_326[k]
                  + f_179 * fk_331[k]
                  - f_258 * fk_333[k]
                  - f_179 * fk_340[k]
                  + f_259 * fk_344[k]
                  - f_179 * fk_353[k]
                  + f_258 * fk_355[k]
                  - f_259 * fk_357[k];
    }

#pragma omp simd aligned(fk_72, fk_75, fk_77, fk_82, fk_84, fk_86, fk_93, fk_95, fk_97, \
                         fk_252, fk_255, fk_257, fk_262, fk_264, fk_266, fk_273, fk_275, \
                         fk_277, fk_324, fk_327, fk_329, fk_334, fk_336, fk_338, fk_345, \
                         fk_347, fk_349 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_55[k] = -f_235 * fk_72[k]
                  + f_235 * fk_75[k]
                  + f_238 * fk_77[k]
                  + f_233 * fk_82[k]
                  - f_236 * fk_84[k]
                  - f_239 * fk_86[k]
                  + f_232 * fk_93[k]
                  - f_234 * fk_95[k]
                  + f_237 * fk_97[k]
                  - f_235 * fk_252[k]
                  + f_235 * fk_255[k]
                  + f_238 * fk_257[k]
                  + f_233 * fk_262[k]
                  - f_236 * fk_264[k]
                  - f_239 * fk_266[k]
                  + f_232 * fk_273[k]
                  - f_234 * fk_275[k]
                  + f_237 * fk_277[k]
                  + f_242 * fk_324[k]
                  - f_242 * fk_327[k]
                  - f_244 * fk_329[k]
                  - f_241 * fk_334[k]
                  + f_239 * fk_336[k]
                  + f_245 * fk_338[k]
                  - f_240 * fk_345[k]
                  + f_236 * fk_347[k]
                  - f_243 * fk_349[k];
    }

#pragma omp simd aligned(fk_74, fk_79, fk_81, fk_88, fk_90, fk_101, fk_103, fk_254, fk_259, \
                         fk_261, fk_268, fk_270, fk_281, fk_283, fk_326, fk_331, fk_333, \
                         fk_340, fk_342, fk_353, fk_355 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_56[k] = f_225 * fk_74[k]
                  - f_260 * fk_79[k]
                  - f_261 * fk_81[k]
                  - f_260 * fk_88[k]
                  + f_220 * fk_90[k]
                  + f_225 * fk_101[k]
                  - f_261 * fk_103[k]
                  + f_225 * fk_254[k]
                  - f_260 * fk_259[k]
                  - f_261 * fk_261[k]
                  - f_260 * fk_268[k]
                  + f_220 * fk_270[k]
                  + f_225 * fk_281[k]
                  - f_261 * fk_283[k]
                  - f_262 * fk_326[k]
                  + f_261 * fk_331[k]
                  + f_263 * fk_333[k]
                  + f_261 * fk_340[k]
                  - f_226 * fk_342[k]
                  - f_262 * fk_353[k]
                  + f_263 * fk_355[k];
    }

#pragma omp simd aligned(fk_72, fk_75, fk_77, fk_82, fk_84, fk_93, fk_95, fk_252, fk_255, \
                         fk_257, fk_262, fk_264, fk_273, fk_275, fk_324, fk_327, fk_329, \
                         fk_334, fk_336, fk_345, fk_347 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_57[k] = f_221 * fk_72[k]
                  - f_219 * fk_75[k]
                  - f_222 * fk_77[k]
                  - f_217 * fk_82[k]
                  + f_220 * fk_84[k]
                  + f_217 * fk_93[k]
                  - f_218 * fk_95[k]
                  + f_221 * fk_252[k]
                  - f_219 * fk_255[k]
                  - f_222 * fk_257[k]
                  - f_217 * fk_262[k]
                  + f_220 * fk_264[k]
                  + f_217 * fk_273[k]
                  - f_218 * fk_275[k]
                  - f_227 * fk_324[k]
                  + f_225 * fk_327[k]
                  + f_228 * fk_329[k]
                  + f_223 * fk_334[k]
                  - f_226 * fk_336[k]
                  - f_223 * fk_345[k]
                  + f_224 * fk_347[k];
    }

#pragma omp simd aligned(fk_74, fk_79, fk_88, fk_101, fk_254, fk_259, fk_268, fk_281, fk_326, \
                         fk_331, fk_340, fk_353 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_58[k] = -f_264 * fk_74[k]
                  + f_265 * fk_79[k]
                  - f_265 * fk_88[k]
                  + f_264 * fk_101[k]
                  - f_264 * fk_254[k]
                  + f_265 * fk_259[k]
                  - f_265 * fk_268[k]
                  + f_264 * fk_281[k]
                  + f_266 * fk_326[k]
                  - f_267 * fk_331[k]
                  + f_267 * fk_340[k]
                  - f_266 * fk_353[k];
    }

#pragma omp simd aligned(fk_72, fk_75, fk_82, fk_93, fk_252, fk_255, fk_262, fk_273, fk_324, \
                         fk_327, fk_334, fk_345 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_59[k] = -f_208 * fk_72[k]
                  + f_207 * fk_75[k]
                  - f_206 * fk_82[k]
                  + f_205 * fk_93[k]
                  - f_208 * fk_252[k]
                  + f_207 * fk_255[k]
                  - f_206 * fk_262[k]
                  + f_205 * fk_273[k]
                  + f_212 * fk_324[k]
                  - f_211 * fk_327[k]
                  + f_210 * fk_334[k]
                  - f_209 * fk_345[k];
    }

#pragma omp simd aligned(fk_1, fk_6, fk_15, fk_28, fk_109, fk_114, fk_123, fk_136, fk_181, \
                         fk_186, fk_195, fk_208 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_60[k] = -f_123 * fk_1[k]
                  + f_124 * fk_6[k]
                  - f_125 * fk_15[k]
                  + f_126 * fk_28[k]
                  - f_123 * fk_109[k]
                  + f_124 * fk_114[k]
                  - f_125 * fk_123[k]
                  + f_126 * fk_136[k]
                  + f_127 * fk_181[k]
                  - f_128 * fk_186[k]
                  + f_129 * fk_195[k]
                  - f_130 * fk_208[k];
    }

#pragma omp simd aligned(fk_4, fk_11, fk_22, fk_112, fk_119, fk_130, fk_184, fk_191, \
                         fk_202 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_61[k] = -f_131 * fk_4[k]
                  + f_132 * fk_11[k]
                  - f_131 * fk_22[k]
                  - f_131 * fk_112[k]
                  + f_132 * fk_119[k]
                  - f_131 * fk_130[k]
                  + f_133 * fk_184[k]
                  - f_134 * fk_191[k]
                  + f_133 * fk_202[k];
    }

#pragma omp simd aligned(fk_1, fk_6, fk_8, fk_15, fk_17, fk_28, fk_30, fk_109, fk_114, fk_116, \
                         fk_123, fk_125, fk_136, fk_138, fk_181, fk_186, fk_188, fk_195, \
                         fk_197, fk_208, fk_210 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_62[k] = f_135 * fk_1[k]
                  - f_135 * fk_6[k]
                  - f_136 * fk_8[k]
                  - f_137 * fk_15[k]
                  + f_138 * fk_17[k]
                  + f_139 * fk_28[k]
                  - f_140 * fk_30[k]
                  + f_135 * fk_109[k]
                  - f_135 * fk_114[k]
                  - f_136 * fk_116[k]
                  - f_137 * fk_123[k]
                  + f_138 * fk_125[k]
                  + f_139 * fk_136[k]
                  - f_140 * fk_138[k]
                  - f_141 * fk_181[k]
                  + f_141 * fk_186[k]
                  + f_142 * fk_188[k]
                  + f_143 * fk_195[k]
                  - f_144 * fk_197[k]
                  - f_145 * fk_208[k]
                  + f_146 * fk_210[k];
    }

#pragma omp simd aligned(fk_4, fk_13, fk_22, fk_24, fk_112, fk_121, fk_130, fk_132, fk_184, \
                         fk_193, fk_202, fk_204 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_63[k] = f_147 * fk_4[k]
                  - f_148 * fk_13[k]
                  - f_147 * fk_22[k]
                  + f_148 * fk_24[k]
                  + f_147 * fk_112[k]
                  - f_148 * fk_121[k]
                  - f_147 * fk_130[k]
                  + f_148 * fk_132[k]
                  - f_149 * fk_184[k]
                  + f_150 * fk_193[k]
                  + f_149 * fk_202[k]
                  - f_150 * fk_204[k];
    }

#pragma omp simd aligned(fk_1, fk_6, fk_8, fk_15, fk_17, fk_19, fk_28, fk_30, fk_32, fk_109, \
                         fk_114, fk_116, fk_123, fk_125, fk_127, fk_136, fk_138, fk_140, \
                         fk_181, fk_186, fk_188, fk_195, fk_197, fk_199, fk_208, fk_210, \
                         fk_212 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_64[k] = -f_151 * fk_1[k]
                  - f_152 * fk_6[k]
                  + f_153 * fk_8[k]
                  - f_154 * fk_15[k]
                  + f_155 * fk_17[k]
                  - f_156 * fk_19[k]
                  + f_154 * fk_28[k]
                  - f_157 * fk_30[k]
                  + f_158 * fk_32[k]
                  - f_151 * fk_109[k]
                  - f_152 * fk_114[k]
                  + f_153 * fk_116[k]
                  - f_154 * fk_123[k]
                  + f_155 * fk_125[k]
                  - f_156 * fk_127[k]
                  + f_154 * fk_136[k]
                  - f_157 * fk_138[k]
                  + f_158 * fk_140[k]
                  + f_159 * fk_181[k]
                  + f_157 * fk_186[k]
                  - f_160 * fk_188[k]
                  + f_161 * fk_195[k]
                  - f_162 * fk_197[k]
                  + f_163 * fk_199[k]
                  - f_161 * fk_208[k]
                  + f_156 * fk_210[k]
                  - f_164 * fk_212[k];
    }

#pragma omp simd aligned(fk_4, fk_11, fk_13, fk_22, fk_24, fk_26, fk_112, fk_119, fk_121, \
                         fk_130, fk_132, fk_134, fk_184, fk_191, fk_193, fk_202, fk_204, \
                         fk_206 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_65[k] = -f_165 * fk_4[k]
                  - f_166 * fk_11[k]
                  + f_167 * fk_13[k]
                  - f_165 * fk_22[k]
                  + f_167 * fk_24[k]
                  - f_168 * fk_26[k]
                  - f_165 * fk_112[k]
                  - f_166 * fk_119[k]
                  + f_167 * fk_121[k]
                  - f_165 * fk_130[k]
                  + f_167 * fk_132[k]
                  - f_168 * fk_134[k]
                  + f_169 * fk_184[k]
                  + f_170 * fk_191[k]
                  - f_171 * fk_193[k]
                  + f_169 * fk_202[k]
                  - f_171 * fk_204[k]
                  + f_172 * fk_206[k];
    }

#pragma omp simd aligned(fk_1, fk_6, fk_8, fk_15, fk_17, fk_19, fk_28, fk_30, fk_32, fk_34, \
                         fk_109, fk_114, fk_116, fk_123, fk_125, fk_127, fk_136, fk_138, \
                         fk_140, fk_142, fk_181, fk_186, fk_188, fk_195, fk_197, fk_199, \
                         fk_208, fk_210, fk_212, fk_214 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_66[k] = f_173 * fk_1[k]
                  + f_174 * fk_6[k]
                  - f_175 * fk_8[k]
                  + f_174 * fk_15[k]
                  - f_176 * fk_17[k]
                  + f_176 * fk_19[k]
                  + f_173 * fk_28[k]
                  - f_175 * fk_30[k]
                  + f_176 * fk_32[k]
                  - f_177 * fk_34[k]
                  + f_173 * fk_109[k]
                  + f_174 * fk_114[k]
                  - f_175 * fk_116[k]
                  + f_174 * fk_123[k]
                  - f_176 * fk_125[k]
                  + f_176 * fk_127[k]
                  + f_173 * fk_136[k]
                  - f_175 * fk_138[k]
                  + f_176 * fk_140[k]
                  - f_177 * fk_142[k]
                  - f_178 * fk_181[k]
                  - f_179 * fk_186[k]
                  + f_180 * fk_188[k]
                  - f_179 * fk_195[k]
                  + f_181 * fk_197[k]
                  - f_181 * fk_199[k]
                  - f_178 * fk_208[k]
                  + f_180 * fk_210[k]
                  - f_181 * fk_212[k]
                  + f_182 * fk_214[k];
    }

#pragma omp simd aligned(fk_2, fk_7, fk_9, fk_16, fk_18, fk_20, fk_29, fk_31, fk_33, fk_35, \
                         fk_110, fk_115, fk_117, fk_124, fk_126, fk_128, fk_137, fk_139, \
                         fk_141, fk_143, fk_182, fk_187, fk_189, fk_196, fk_198, fk_200, \
                         fk_209, fk_211, fk_213, fk_215 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_67[k] = f_183 * fk_2[k]
                  + f_184 * fk_7[k]
                  - f_185 * fk_9[k]
                  + f_184 * fk_16[k]
                  - f_186 * fk_18[k]
                  + f_187 * fk_20[k]
                  + f_183 * fk_29[k]
                  - f_185 * fk_31[k]
                  + f_187 * fk_33[k]
                  - f_188 * fk_35[k]
                  + f_183 * fk_110[k]
                  + f_184 * fk_115[k]
                  - f_185 * fk_117[k]
                  + f_184 * fk_124[k]
                  - f_186 * fk_126[k]
                  + f_187 * fk_128[k]
                  + f_183 * fk_137[k]
                  - f_185 * fk_139[k]
                  + f_187 * fk_141[k]
                  - f_188 * fk_143[k]
                  - f_189 * fk_182[k]
                  - f_186 * fk_187[k]
                  + f_190 * fk_189[k]
                  - f_186 * fk_196[k]
                  + f_191 * fk_198[k]
                  - f_192 * fk_200[k]
                  - f_189 * fk_209[k]
                  + f_190 * fk_211[k]
                  - f_192 * fk_213[k]
                  + f_193 * fk_215[k];
    }

#pragma omp simd aligned(fk_0, fk_3, fk_5, fk_10, fk_12, fk_14, fk_21, fk_23, fk_25, fk_27, \
                         fk_108, fk_111, fk_113, fk_118, fk_120, fk_122, fk_129, fk_131, \
                         fk_133, fk_135, fk_180, fk_183, fk_185, fk_190, fk_192, fk_194, \
                         fk_201, fk_203, fk_205, fk_207 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_68[k] = f_173 * fk_0[k]
                  + f_174 * fk_3[k]
                  - f_175 * fk_5[k]
                  + f_174 * fk_10[k]
                  - f_176 * fk_12[k]
                  + f_176 * fk_14[k]
                  + f_173 * fk_21[k]
                  - f_175 * fk_23[k]
                  + f_176 * fk_25[k]
                  - f_177 * fk_27[k]
                  + f_173 * fk_108[k]
                  + f_174 * fk_111[k]
                  - f_175 * fk_113[k]
                  + f_174 * fk_118[k]
                  - f_176 * fk_120[k]
                  + f_176 * fk_122[k]
                  + f_173 * fk_129[k]
                  - f_175 * fk_131[k]
                  + f_176 * fk_133[k]
                  - f_177 * fk_135[k]
                  - f_178 * fk_180[k]
                  - f_179 * fk_183[k]
                  + f_180 * fk_185[k]
                  - f_179 * fk_190[k]
                  + f_181 * fk_192[k]
                  - f_181 * fk_194[k]
                  - f_178 * fk_201[k]
                  + f_180 * fk_203[k]
                  - f_181 * fk_205[k]
                  + f_182 * fk_207[k];
    }

#pragma omp simd aligned(fk_2, fk_7, fk_9, fk_16, fk_20, fk_29, fk_31, fk_33, fk_110, fk_115, \
                         fk_117, fk_124, fk_128, fk_137, fk_139, fk_141, fk_182, fk_187, \
                         fk_189, fk_196, fk_200, fk_209, fk_211, \
                         fk_213 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_69[k] = -f_194 * fk_2[k]
                  - f_194 * fk_7[k]
                  + f_195 * fk_9[k]
                  + f_194 * fk_16[k]
                  - f_196 * fk_20[k]
                  + f_194 * fk_29[k]
                  - f_195 * fk_31[k]
                  + f_196 * fk_33[k]
                  - f_194 * fk_110[k]
                  - f_194 * fk_115[k]
                  + f_195 * fk_117[k]
                  + f_194 * fk_124[k]
                  - f_196 * fk_128[k]
                  + f_194 * fk_137[k]
                  - f_195 * fk_139[k]
                  + f_196 * fk_141[k]
                  + f_166 * fk_182[k]
                  + f_166 * fk_187[k]
                  - f_197 * fk_189[k]
                  - f_166 * fk_196[k]
                  + f_198 * fk_200[k]
                  - f_166 * fk_209[k]
                  + f_197 * fk_211[k]
                  - f_198 * fk_213[k];
    }

#pragma omp simd aligned(fk_0, fk_3, fk_5, fk_10, fk_12, fk_14, fk_21, fk_23, fk_25, fk_108, \
                         fk_111, fk_113, fk_118, fk_120, fk_122, fk_129, fk_131, fk_133, \
                         fk_180, fk_183, fk_185, fk_190, fk_192, fk_194, fk_201, fk_203, \
                         fk_205 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_70[k] = -f_154 * fk_0[k]
                  + f_154 * fk_3[k]
                  + f_157 * fk_5[k]
                  + f_152 * fk_10[k]
                  - f_155 * fk_12[k]
                  - f_158 * fk_14[k]
                  + f_151 * fk_21[k]
                  - f_153 * fk_23[k]
                  + f_156 * fk_25[k]
                  - f_154 * fk_108[k]
                  + f_154 * fk_111[k]
                  + f_157 * fk_113[k]
                  + f_152 * fk_118[k]
                  - f_155 * fk_120[k]
                  - f_158 * fk_122[k]
                  + f_151 * fk_129[k]
                  - f_153 * fk_131[k]
                  + f_156 * fk_133[k]
                  + f_161 * fk_180[k]
                  - f_161 * fk_183[k]
                  - f_156 * fk_185[k]
                  - f_157 * fk_190[k]
                  + f_162 * fk_192[k]
                  + f_164 * fk_194[k]
                  - f_159 * fk_201[k]
                  + f_160 * fk_203[k]
                  - f_163 * fk_205[k];
    }

#pragma omp simd aligned(fk_2, fk_7, fk_9, fk_16, fk_18, fk_29, fk_31, fk_110, fk_115, fk_117, \
                         fk_124, fk_126, fk_137, fk_139, fk_182, fk_187, fk_189, fk_196, \
                         fk_198, fk_209, fk_211 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_71[k] = f_199 * fk_2[k]
                  - f_200 * fk_7[k]
                  - f_141 * fk_9[k]
                  - f_200 * fk_16[k]
                  + f_138 * fk_18[k]
                  + f_199 * fk_29[k]
                  - f_141 * fk_31[k]
                  + f_199 * fk_110[k]
                  - f_200 * fk_115[k]
                  - f_141 * fk_117[k]
                  - f_200 * fk_124[k]
                  + f_138 * fk_126[k]
                  + f_199 * fk_137[k]
                  - f_141 * fk_139[k]
                  - f_147 * fk_182[k]
                  + f_138 * fk_187[k]
                  + f_148 * fk_189[k]
                  + f_138 * fk_196[k]
                  - f_144 * fk_198[k]
                  - f_147 * fk_209[k]
                  + f_148 * fk_211[k];
    }

#pragma omp simd aligned(fk_0, fk_3, fk_5, fk_10, fk_12, fk_21, fk_23, fk_108, fk_111, fk_113, \
                         fk_118, fk_120, fk_129, fk_131, fk_180, fk_183, fk_185, fk_190, \
                         fk_192, fk_201, fk_203 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_72[k] = f_139 * fk_0[k]
                  - f_137 * fk_3[k]
                  - f_140 * fk_5[k]
                  - f_135 * fk_10[k]
                  + f_138 * fk_12[k]
                  + f_135 * fk_21[k]
                  - f_136 * fk_23[k]
                  + f_139 * fk_108[k]
                  - f_137 * fk_111[k]
                  - f_140 * fk_113[k]
                  - f_135 * fk_118[k]
                  + f_138 * fk_120[k]
                  + f_135 * fk_129[k]
                  - f_136 * fk_131[k]
                  - f_145 * fk_180[k]
                  + f_143 * fk_183[k]
                  + f_146 * fk_185[k]
                  + f_141 * fk_190[k]
                  - f_144 * fk_192[k]
                  - f_141 * fk_201[k]
                  + f_142 * fk_203[k];
    }

#pragma omp simd aligned(fk_2, fk_7, fk_16, fk_29, fk_110, fk_115, fk_124, fk_137, fk_182, \
                         fk_187, fk_196, fk_209 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_73[k] = -f_201 * fk_2[k]
                  + f_202 * fk_7[k]
                  - f_202 * fk_16[k]
                  + f_201 * fk_29[k]
                  - f_201 * fk_110[k]
                  + f_202 * fk_115[k]
                  - f_202 * fk_124[k]
                  + f_201 * fk_137[k]
                  + f_203 * fk_182[k]
                  - f_204 * fk_187[k]
                  + f_204 * fk_196[k]
                  - f_203 * fk_209[k];
    }

#pragma omp simd aligned(fk_0, fk_3, fk_10, fk_21, fk_108, fk_111, fk_118, fk_129, fk_180, \
                         fk_183, fk_190, fk_201 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_74[k] = -f_126 * fk_0[k]
                  + f_125 * fk_3[k]
                  - f_124 * fk_10[k]
                  + f_123 * fk_21[k]
                  - f_126 * fk_108[k]
                  + f_125 * fk_111[k]
                  - f_124 * fk_118[k]
                  + f_123 * fk_129[k]
                  + f_130 * fk_180[k]
                  - f_129 * fk_183[k]
                  + f_128 * fk_190[k]
                  - f_127 * fk_201[k];
    }

#pragma omp simd aligned(fk_73, fk_76, fk_78, fk_83, fk_87, fk_94, fk_100, fk_253, fk_256, \
                         fk_258, fk_263, fk_267, fk_274, fk_280 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_75[k] = f_268 * fk_73[k]
                  - f_269 * fk_78[k]
                  + f_270 * fk_87[k]
                  - f_271 * fk_100[k]
                  - f_268 * fk_253[k]
                  + f_269 * fk_258[k]
                  - f_270 * fk_267[k]
                  + f_271 * fk_280[k];

        g_76[k] = f_272 * fk_76[k]
                  - f_273 * fk_83[k]
                  + f_272 * fk_94[k]
                  - f_272 * fk_256[k]
                  + f_273 * fk_263[k]
                  - f_272 * fk_274[k];
    }

#pragma omp simd aligned(fk_73, fk_78, fk_80, fk_87, fk_89, fk_100, fk_102, fk_253, fk_258, \
                         fk_260, fk_267, fk_269, fk_280, fk_282 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_77[k] = -f_274 * fk_73[k]
                  + f_274 * fk_78[k]
                  + f_119 * fk_80[k]
                  + f_275 * fk_87[k]
                  - f_90 * fk_89[k]
                  - f_276 * fk_100[k]
                  + f_118 * fk_102[k]
                  + f_274 * fk_253[k]
                  - f_274 * fk_258[k]
                  - f_119 * fk_260[k]
                  - f_275 * fk_267[k]
                  + f_90 * fk_269[k]
                  + f_276 * fk_280[k]
                  - f_118 * fk_282[k];
    }

#pragma omp simd aligned(fk_76, fk_85, fk_94, fk_96, fk_256, fk_265, fk_274, \
                         fk_276 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_78[k] = -f_94 * fk_76[k]
                  + f_277 * fk_85[k]
                  + f_94 * fk_94[k]
                  - f_277 * fk_96[k]
                  + f_94 * fk_256[k]
                  - f_277 * fk_265[k]
                  - f_94 * fk_274[k]
                  + f_277 * fk_276[k];
    }

#pragma omp simd aligned(fk_73, fk_78, fk_80, fk_87, fk_89, fk_91, fk_100, fk_102, fk_104, \
                         fk_253, fk_258, fk_260, fk_267, fk_269, fk_271, fk_280, fk_282, \
                         fk_284 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_79[k] = f_278 * fk_73[k]
                  + f_279 * fk_78[k]
                  - f_280 * fk_80[k]
                  + f_281 * fk_87[k]
                  - f_103 * fk_89[k]
                  + f_101 * fk_91[k]
                  - f_281 * fk_100[k]
                  + f_282 * fk_102[k]
                  - f_283 * fk_104[k]
                  - f_278 * fk_253[k]
                  - f_279 * fk_258[k]
                  + f_280 * fk_260[k]
                  - f_281 * fk_267[k]
                  + f_103 * fk_269[k]
                  - f_101 * fk_271[k]
                  + f_281 * fk_280[k]
                  - f_282 * fk_282[k]
                  + f_283 * fk_284[k];
    }

#pragma omp simd aligned(fk_76, fk_83, fk_85, fk_94, fk_96, fk_98, fk_256, fk_263, fk_265, \
                         fk_274, fk_276, fk_278 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_80[k] = f_115 * fk_76[k]
                  + f_50 * fk_83[k]
                  - f_116 * fk_85[k]
                  + f_115 * fk_94[k]
                  - f_116 * fk_96[k]
                  + f_117 * fk_98[k]
                  - f_115 * fk_256[k]
                  - f_50 * fk_263[k]
                  + f_116 * fk_265[k]
                  - f_115 * fk_274[k]
                  + f_116 * fk_276[k]
                  - f_117 * fk_278[k];
    }

#pragma omp simd aligned(fk_73, fk_78, fk_80, fk_87, fk_89, fk_91, fk_100, fk_102, fk_104, \
                         fk_106, fk_253, fk_258, fk_260, fk_267, fk_269, fk_271, fk_280, \
                         fk_282, fk_284, fk_286 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_81[k] = -f_284 * fk_73[k]
                  - f_71 * fk_78[k]
                  + f_285 * fk_80[k]
                  - f_71 * fk_87[k]
                  + f_69 * fk_89[k]
                  - f_69 * fk_91[k]
                  - f_284 * fk_100[k]
                  + f_285 * fk_102[k]
                  - f_69 * fk_104[k]
                  + f_286 * fk_106[k]
                  + f_284 * fk_253[k]
                  + f_71 * fk_258[k]
                  - f_285 * fk_260[k]
                  + f_71 * fk_267[k]
                  - f_69 * fk_269[k]
                  + f_69 * fk_271[k]
                  + f_284 * fk_280[k]
                  - f_285 * fk_282[k]
                  + f_69 * fk_284[k]
                  - f_286 * fk_286[k];
    }

#pragma omp simd aligned(fk_74, fk_79, fk_81, fk_88, fk_90, fk_92, fk_101, fk_103, fk_105, \
                         fk_107, fk_254, fk_259, fk_261, fk_268, fk_270, fk_272, fk_281, \
                         fk_283, fk_285, fk_287 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_82[k] = -f_287 * fk_74[k]
                  - f_288 * fk_79[k]
                  + f_110 * fk_81[k]
                  - f_288 * fk_88[k]
                  + f_111 * fk_90[k]
                  - f_289 * fk_92[k]
                  - f_287 * fk_101[k]
                  + f_110 * fk_103[k]
                  - f_289 * fk_105[k]
                  + f_290 * fk_107[k]
                  + f_287 * fk_254[k]
                  + f_288 * fk_259[k]
                  - f_110 * fk_261[k]
                  + f_288 * fk_268[k]
                  - f_111 * fk_270[k]
                  + f_289 * fk_272[k]
                  + f_287 * fk_281[k]
                  - f_110 * fk_283[k]
                  + f_289 * fk_285[k]
                  - f_290 * fk_287[k];
    }

#pragma omp simd aligned(fk_72, fk_75, fk_77, fk_82, fk_84, fk_86, fk_93, fk_95, fk_97, fk_99, \
                         fk_252, fk_255, fk_257, fk_262, fk_264, fk_266, fk_273, fk_275, \
                         fk_277, fk_279 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_83[k] = -f_284 * fk_72[k]
                  - f_71 * fk_75[k]
                  + f_285 * fk_77[k]
                  - f_71 * fk_82[k]
                  + f_69 * fk_84[k]
                  - f_69 * fk_86[k]
                  - f_284 * fk_93[k]
                  + f_285 * fk_95[k]
                  - f_69 * fk_97[k]
                  + f_286 * fk_99[k]
                  + f_284 * fk_252[k]
                  + f_71 * fk_255[k]
                  - f_285 * fk_257[k]
                  + f_71 * fk_262[k]
                  - f_69 * fk_264[k]
                  + f_69 * fk_266[k]
                  + f_284 * fk_273[k]
                  - f_285 * fk_275[k]
                  + f_69 * fk_277[k]
                  - f_286 * fk_279[k];
    }

#pragma omp simd aligned(fk_74, fk_79, fk_81, fk_88, fk_92, fk_101, fk_103, fk_105, fk_254, \
                         fk_259, fk_261, fk_268, fk_272, fk_281, fk_283, \
                         fk_285 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_84[k] = f_291 * fk_74[k]
                  + f_291 * fk_79[k]
                  - f_292 * fk_81[k]
                  - f_291 * fk_88[k]
                  + f_293 * fk_92[k]
                  - f_291 * fk_101[k]
                  + f_292 * fk_103[k]
                  - f_293 * fk_105[k]
                  - f_291 * fk_254[k]
                  - f_291 * fk_259[k]
                  + f_292 * fk_261[k]
                  + f_291 * fk_268[k]
                  - f_293 * fk_272[k]
                  + f_291 * fk_281[k]
                  - f_292 * fk_283[k]
                  + f_293 * fk_285[k];
    }

#pragma omp simd aligned(fk_72, fk_75, fk_77, fk_82, fk_84, fk_86, fk_93, fk_95, fk_97, \
                         fk_252, fk_255, fk_257, fk_262, fk_264, fk_266, fk_273, fk_275, \
                         fk_277 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_85[k] = f_281 * fk_72[k]
                  - f_281 * fk_75[k]
                  - f_282 * fk_77[k]
                  - f_279 * fk_82[k]
                  + f_103 * fk_84[k]
                  + f_283 * fk_86[k]
                  - f_278 * fk_93[k]
                  + f_280 * fk_95[k]
                  - f_101 * fk_97[k]
                  - f_281 * fk_252[k]
                  + f_281 * fk_255[k]
                  + f_282 * fk_257[k]
                  + f_279 * fk_262[k]
                  - f_103 * fk_264[k]
                  - f_283 * fk_266[k]
                  + f_278 * fk_273[k]
                  - f_280 * fk_275[k]
                  + f_101 * fk_277[k];
    }

#pragma omp simd aligned(fk_74, fk_79, fk_81, fk_88, fk_90, fk_101, fk_103, fk_254, fk_259, \
                         fk_261, fk_268, fk_270, fk_281, fk_283 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_86[k] = -f_294 * fk_74[k]
                  + f_295 * fk_79[k]
                  + f_296 * fk_81[k]
                  + f_295 * fk_88[k]
                  - f_90 * fk_90[k]
                  - f_294 * fk_101[k]
                  + f_296 * fk_103[k]
                  + f_294 * fk_254[k]
                  - f_295 * fk_259[k]
                  - f_296 * fk_261[k]
                  - f_295 * fk_268[k]
                  + f_90 * fk_270[k]
                  + f_294 * fk_281[k]
                  - f_296 * fk_283[k];
    }

#pragma omp simd aligned(fk_72, fk_75, fk_77, fk_82, fk_84, fk_93, fk_95, fk_252, fk_255, \
                         fk_257, fk_262, fk_264, fk_273, fk_275 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_87[k] = -f_276 * fk_72[k]
                  + f_275 * fk_75[k]
                  + f_118 * fk_77[k]
                  + f_274 * fk_82[k]
                  - f_90 * fk_84[k]
                  - f_274 * fk_93[k]
                  + f_119 * fk_95[k]
                  + f_276 * fk_252[k]
                  - f_275 * fk_255[k]
                  - f_118 * fk_257[k]
                  - f_274 * fk_262[k]
                  + f_90 * fk_264[k]
                  + f_274 * fk_273[k]
                  - f_119 * fk_275[k];
    }

#pragma omp simd aligned(fk_74, fk_79, fk_88, fk_101, fk_254, fk_259, fk_268, \
                         fk_281 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_88[k] = f_297 * fk_74[k]
                  - f_298 * fk_79[k]
                  + f_298 * fk_88[k]
                  - f_297 * fk_101[k]
                  - f_297 * fk_254[k]
                  + f_298 * fk_259[k]
                  - f_298 * fk_268[k]
                  + f_297 * fk_281[k];
    }

#pragma omp simd aligned(fk_72, fk_75, fk_82, fk_93, fk_252, fk_255, fk_262, \
                         fk_273 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_89[k] = f_271 * fk_72[k]
                  - f_270 * fk_75[k]
                  + f_269 * fk_82[k]
                  - f_268 * fk_93[k]
                  - f_271 * fk_252[k]
                  + f_270 * fk_255[k]
                  - f_269 * fk_262[k]
                  + f_268 * fk_273[k];
    }

#pragma omp simd aligned(fk_1, fk_4, fk_6, fk_11, fk_15, fk_22, fk_28, fk_109, fk_112, fk_114, \
                         fk_119, fk_123, fk_130, fk_136 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_90[k] = f_4 * fk_1[k]
                  - f_5 * fk_6[k]
                  + f_0 * fk_15[k]
                  - f_6 * fk_28[k]
                  - f_0 * fk_109[k]
                  + f_1 * fk_114[k]
                  - f_2 * fk_123[k]
                  + f_3 * fk_136[k];

        g_91[k] = f_9 * fk_4[k]
                  - f_10 * fk_11[k]
                  + f_9 * fk_22[k]
                  - f_7 * fk_112[k]
                  + f_8 * fk_119[k]
                  - f_7 * fk_130[k];
    }

#pragma omp simd aligned(fk_1, fk_6, fk_8, fk_15, fk_17, fk_28, fk_30, fk_109, fk_114, fk_116, \
                         fk_123, fk_125, fk_136, fk_138 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_92[k] = -f_17 * fk_1[k]
                  + f_17 * fk_6[k]
                  + f_18 * fk_8[k]
                  + f_19 * fk_15[k]
                  - f_20 * fk_17[k]
                  - f_21 * fk_28[k]
                  + f_22 * fk_30[k]
                  + f_11 * fk_109[k]
                  - f_11 * fk_114[k]
                  - f_12 * fk_116[k]
                  - f_13 * fk_123[k]
                  + f_14 * fk_125[k]
                  + f_15 * fk_136[k]
                  - f_16 * fk_138[k];
    }

#pragma omp simd aligned(fk_4, fk_13, fk_22, fk_24, fk_112, fk_121, fk_130, \
                         fk_132 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_93[k] = -f_25 * fk_4[k]
                  + f_26 * fk_13[k]
                  + f_25 * fk_22[k]
                  - f_26 * fk_24[k]
                  + f_23 * fk_112[k]
                  - f_24 * fk_121[k]
                  - f_23 * fk_130[k]
                  + f_24 * fk_132[k];
    }

#pragma omp simd aligned(fk_1, fk_6, fk_8, fk_15, fk_17, fk_19, fk_28, fk_30, fk_32, fk_109, \
                         fk_114, fk_116, fk_123, fk_125, fk_127, fk_136, fk_138, \
                         fk_140 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_94[k] = f_30 * fk_1[k]
                  + f_35 * fk_6[k]
                  - f_33 * fk_8[k]
                  + f_36 * fk_15[k]
                  - f_37 * fk_17[k]
                  + f_34 * fk_19[k]
                  - f_36 * fk_28[k]
                  + f_38 * fk_30[k]
                  - f_39 * fk_32[k]
                  - f_27 * fk_109[k]
                  - f_28 * fk_114[k]
                  + f_29 * fk_116[k]
                  - f_30 * fk_123[k]
                  + f_31 * fk_125[k]
                  - f_32 * fk_127[k]
                  + f_30 * fk_136[k]
                  - f_33 * fk_138[k]
                  + f_34 * fk_140[k];
    }

#pragma omp simd aligned(fk_4, fk_11, fk_13, fk_22, fk_24, fk_26, fk_112, fk_119, fk_121, \
                         fk_130, fk_132, fk_134 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_95[k] = f_44 * fk_4[k]
                  + f_45 * fk_11[k]
                  - f_46 * fk_13[k]
                  + f_44 * fk_22[k]
                  - f_46 * fk_24[k]
                  + f_47 * fk_26[k]
                  - f_40 * fk_112[k]
                  - f_41 * fk_119[k]
                  + f_42 * fk_121[k]
                  - f_40 * fk_130[k]
                  + f_42 * fk_132[k]
                  - f_43 * fk_134[k];
    }

#pragma omp simd aligned(fk_1, fk_6, fk_8, fk_15, fk_17, fk_19, fk_28, fk_30, fk_32, fk_34, \
                         fk_109, fk_114, fk_116, fk_123, fk_125, fk_127, fk_136, fk_138, \
                         fk_140, fk_142 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_96[k] = -f_53 * fk_1[k]
                  - f_48 * fk_6[k]
                  + f_54 * fk_8[k]
                  - f_48 * fk_15[k]
                  + f_55 * fk_17[k]
                  - f_55 * fk_19[k]
                  - f_53 * fk_28[k]
                  + f_54 * fk_30[k]
                  - f_55 * fk_32[k]
                  + f_56 * fk_34[k]
                  + f_48 * fk_109[k]
                  + f_49 * fk_114[k]
                  - f_50 * fk_116[k]
                  + f_49 * fk_123[k]
                  - f_51 * fk_125[k]
                  + f_51 * fk_127[k]
                  + f_48 * fk_136[k]
                  - f_50 * fk_138[k]
                  + f_51 * fk_140[k]
                  - f_52 * fk_142[k];
    }

#pragma omp simd aligned(fk_2, fk_7, fk_9, fk_16, fk_18, fk_20, fk_29, fk_31, fk_33, fk_35, \
                         fk_110, fk_115, fk_117, fk_124, fk_126, fk_128, fk_137, fk_139, \
                         fk_141, fk_143 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_97[k] = -f_63 * fk_2[k]
                  - f_57 * fk_7[k]
                  + f_64 * fk_9[k]
                  - f_57 * fk_16[k]
                  + f_65 * fk_18[k]
                  - f_66 * fk_20[k]
                  - f_63 * fk_29[k]
                  + f_64 * fk_31[k]
                  - f_66 * fk_33[k]
                  + f_67 * fk_35[k]
                  + f_57 * fk_110[k]
                  + f_58 * fk_115[k]
                  - f_59 * fk_117[k]
                  + f_58 * fk_124[k]
                  - f_60 * fk_126[k]
                  + f_61 * fk_128[k]
                  + f_57 * fk_137[k]
                  - f_59 * fk_139[k]
                  + f_61 * fk_141[k]
                  - f_62 * fk_143[k];
    }

#pragma omp simd aligned(fk_0, fk_3, fk_5, fk_10, fk_12, fk_14, fk_21, fk_23, fk_25, fk_27, \
                         fk_108, fk_111, fk_113, fk_118, fk_120, fk_122, fk_129, fk_131, \
                         fk_133, fk_135 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_98[k] = -f_53 * fk_0[k]
                  - f_48 * fk_3[k]
                  + f_54 * fk_5[k]
                  - f_48 * fk_10[k]
                  + f_55 * fk_12[k]
                  - f_55 * fk_14[k]
                  - f_53 * fk_21[k]
                  + f_54 * fk_23[k]
                  - f_55 * fk_25[k]
                  + f_56 * fk_27[k]
                  + f_48 * fk_108[k]
                  + f_49 * fk_111[k]
                  - f_50 * fk_113[k]
                  + f_49 * fk_118[k]
                  - f_51 * fk_120[k]
                  + f_51 * fk_122[k]
                  + f_48 * fk_129[k]
                  - f_50 * fk_131[k]
                  + f_51 * fk_133[k]
                  - f_52 * fk_135[k];
    }

#pragma omp simd aligned(fk_2, fk_7, fk_9, fk_16, fk_20, fk_29, fk_31, fk_33, fk_110, fk_115, \
                         fk_117, fk_124, fk_128, fk_137, fk_139, \
                         fk_141 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_99[k] = f_71 * fk_2[k]
                  + f_71 * fk_7[k]
                  - f_72 * fk_9[k]
                  - f_71 * fk_16[k]
                  + f_73 * fk_20[k]
                  - f_71 * fk_29[k]
                  + f_72 * fk_31[k]
                  - f_73 * fk_33[k]
                  - f_68 * fk_110[k]
                  - f_68 * fk_115[k]
                  + f_69 * fk_117[k]
                  + f_68 * fk_124[k]
                  - f_70 * fk_128[k]
                  + f_68 * fk_137[k]
                  - f_69 * fk_139[k]
                  + f_70 * fk_141[k];
    }

#pragma omp simd aligned(fk_0, fk_3, fk_5, fk_10, fk_12, fk_14, fk_21, fk_23, fk_25, fk_108, \
                         fk_111, fk_113, fk_118, fk_120, fk_122, fk_129, fk_131, \
                         fk_133 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_100[k] = f_36 * fk_0[k]
                   - f_36 * fk_3[k]
                   - f_38 * fk_5[k]
                   - f_35 * fk_10[k]
                   + f_37 * fk_12[k]
                   + f_39 * fk_14[k]
                   - f_30 * fk_21[k]
                   + f_33 * fk_23[k]
                   - f_34 * fk_25[k]
                   - f_30 * fk_108[k]
                   + f_30 * fk_111[k]
                   + f_33 * fk_113[k]
                   + f_28 * fk_118[k]
                   - f_31 * fk_120[k]
                   - f_34 * fk_122[k]
                   + f_27 * fk_129[k]
                   - f_29 * fk_131[k]
                   + f_32 * fk_133[k];
    }

#pragma omp simd aligned(fk_2, fk_7, fk_9, fk_16, fk_18, fk_29, fk_31, fk_110, fk_115, fk_117, \
                         fk_124, fk_126, fk_137, fk_139 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_101[k] = -f_76 * fk_2[k]
                   + f_77 * fk_7[k]
                   + f_78 * fk_9[k]
                   + f_77 * fk_16[k]
                   - f_20 * fk_18[k]
                   - f_76 * fk_29[k]
                   + f_78 * fk_31[k]
                   + f_74 * fk_110[k]
                   - f_75 * fk_115[k]
                   - f_18 * fk_117[k]
                   - f_75 * fk_124[k]
                   + f_14 * fk_126[k]
                   + f_74 * fk_137[k]
                   - f_18 * fk_139[k];
    }

#pragma omp simd aligned(fk_0, fk_3, fk_5, fk_10, fk_12, fk_21, fk_23, fk_108, fk_111, fk_113, \
                         fk_118, fk_120, fk_129, fk_131 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_102[k] = -f_21 * fk_0[k]
                   + f_19 * fk_3[k]
                   + f_22 * fk_5[k]
                   + f_17 * fk_10[k]
                   - f_20 * fk_12[k]
                   - f_17 * fk_21[k]
                   + f_18 * fk_23[k]
                   + f_15 * fk_108[k]
                   - f_13 * fk_111[k]
                   - f_16 * fk_113[k]
                   - f_11 * fk_118[k]
                   + f_14 * fk_120[k]
                   + f_11 * fk_129[k]
                   - f_12 * fk_131[k];
    }

#pragma omp simd aligned(fk_2, fk_7, fk_16, fk_29, fk_110, fk_115, fk_124, \
                         fk_137 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_103[k] = f_81 * fk_2[k]
                   - f_82 * fk_7[k]
                   + f_82 * fk_16[k]
                   - f_81 * fk_29[k]
                   - f_79 * fk_110[k]
                   + f_80 * fk_115[k]
                   - f_80 * fk_124[k]
                   + f_79 * fk_137[k];
    }

#pragma omp simd aligned(fk_0, fk_3, fk_10, fk_21, fk_108, fk_111, fk_118, \
                         fk_129 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_104[k] = f_6 * fk_0[k]
                   - f_0 * fk_3[k]
                   + f_5 * fk_10[k]
                   - f_4 * fk_21[k]
                   - f_3 * fk_108[k]
                   + f_2 * fk_111[k]
                   - f_1 * fk_118[k]
                   + f_0 * fk_129[k];
    }
}

}  // namespace simdtrf
