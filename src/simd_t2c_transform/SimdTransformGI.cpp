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


#include "SimdTransformGI.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_gi(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t gi,
             const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 0.65625 * std::sqrt(330.0);
    const auto f_1 = 2.1875 * std::sqrt(330.0);
    const auto f_2 = 3.28125 * std::sqrt(110.0);
    const auto f_3 = 6.5625 * std::sqrt(110.0);
    const auto f_4 = 0.65625 * std::sqrt(110.0);
    const auto f_5 = 2.625 * std::sqrt(5.0);
    const auto f_6 = 26.25 * std::sqrt(5.0);
    const auto f_7 = 9.84375 * std::sqrt(6.0);
    const auto f_8 = 6.5625 * std::sqrt(6.0);
    const auto f_9 = 26.25 * std::sqrt(6.0);
    const auto f_10 = 3.28125 * std::sqrt(6.0);
    const auto f_11 = 8.75 * std::sqrt(6.0);
    const auto f_12 = 1.09375 * std::sqrt(6.0);
    const auto f_13 = 2.1875 * std::sqrt(6.0);
    const auto f_14 = 17.5 * std::sqrt(6.0);
    const auto f_15 = 2.1875 * std::sqrt(15.0);
    const auto f_16 = 4.375 * std::sqrt(15.0);
    const auto f_17 = 8.75 * std::sqrt(15.0);
    const auto f_18 = 3.5 * std::sqrt(15.0);
    const auto f_19 = 0.15625 * std::sqrt(35.0);
    const auto f_20 = 0.46875 * std::sqrt(35.0);
    const auto f_21 = 2.8125 * std::sqrt(35.0);
    const auto f_22 = 5.625 * std::sqrt(35.0);
    const auto f_23 = 3.75 * std::sqrt(35.0);
    const auto f_24 = 0.5 * std::sqrt(35.0);
    const auto f_25 = 0.546875 * std::sqrt(6.0);
    const auto f_26 = 0.65625 * std::sqrt(5.0);
    const auto f_27 = 3.28125 * std::sqrt(5.0);
    const auto f_28 = 6.5625 * std::sqrt(5.0);
    const auto f_29 = 39.375 * std::sqrt(5.0);
    const auto f_30 = 0.109375 * std::sqrt(330.0);
    const auto f_31 = 1.640625 * std::sqrt(330.0);
    const auto f_32 = 1.96875 * std::sqrt(165.0);
    const auto f_33 = 6.5625 * std::sqrt(165.0);
    const auto f_34 = 0.65625 * std::sqrt(165.0);
    const auto f_35 = 2.1875 * std::sqrt(165.0);
    const auto f_36 = 9.84375 * std::sqrt(55.0);
    const auto f_37 = 19.6875 * std::sqrt(55.0);
    const auto f_38 = 1.96875 * std::sqrt(55.0);
    const auto f_39 = 3.28125 * std::sqrt(55.0);
    const auto f_40 = 6.5625 * std::sqrt(55.0);
    const auto f_41 = 0.65625 * std::sqrt(55.0);
    const auto f_42 = 3.9375 * std::sqrt(10.0);
    const auto f_43 = 39.375 * std::sqrt(10.0);
    const auto f_44 = 1.3125 * std::sqrt(10.0);
    const auto f_45 = 13.125 * std::sqrt(10.0);
    const auto f_46 = 29.53125 * std::sqrt(3.0);
    const auto f_47 = 19.6875 * std::sqrt(3.0);
    const auto f_48 = 78.75 * std::sqrt(3.0);
    const auto f_49 = 9.84375 * std::sqrt(3.0);
    const auto f_50 = 26.25 * std::sqrt(3.0);
    const auto f_51 = 6.5625 * std::sqrt(3.0);
    const auto f_52 = 3.28125 * std::sqrt(3.0);
    const auto f_53 = 8.75 * std::sqrt(3.0);
    const auto f_54 = 52.5 * std::sqrt(3.0);
    const auto f_55 = 1.09375 * std::sqrt(3.0);
    const auto f_56 = 2.1875 * std::sqrt(3.0);
    const auto f_57 = 17.5 * std::sqrt(3.0);
    const auto f_58 = 3.28125 * std::sqrt(30.0);
    const auto f_59 = 6.5625 * std::sqrt(30.0);
    const auto f_60 = 13.125 * std::sqrt(30.0);
    const auto f_61 = 5.25 * std::sqrt(30.0);
    const auto f_62 = 1.09375 * std::sqrt(30.0);
    const auto f_63 = 2.1875 * std::sqrt(30.0);
    const auto f_64 = 4.375 * std::sqrt(30.0);
    const auto f_65 = 1.75 * std::sqrt(30.0);
    const auto f_66 = 0.234375 * std::sqrt(70.0);
    const auto f_67 = 0.703125 * std::sqrt(70.0);
    const auto f_68 = 4.21875 * std::sqrt(70.0);
    const auto f_69 = 8.4375 * std::sqrt(70.0);
    const auto f_70 = 5.625 * std::sqrt(70.0);
    const auto f_71 = 0.75 * std::sqrt(70.0);
    const auto f_72 = 0.078125 * std::sqrt(70.0);
    const auto f_73 = 1.40625 * std::sqrt(70.0);
    const auto f_74 = 2.8125 * std::sqrt(70.0);
    const auto f_75 = 1.875 * std::sqrt(70.0);
    const auto f_76 = 0.25 * std::sqrt(70.0);
    const auto f_77 = 1.640625 * std::sqrt(3.0);
    const auto f_78 = 0.546875 * std::sqrt(3.0);
    const auto f_79 = 0.984375 * std::sqrt(10.0);
    const auto f_80 = 4.921875 * std::sqrt(10.0);
    const auto f_81 = 9.84375 * std::sqrt(10.0);
    const auto f_82 = 59.0625 * std::sqrt(10.0);
    const auto f_83 = 0.328125 * std::sqrt(10.0);
    const auto f_84 = 1.640625 * std::sqrt(10.0);
    const auto f_85 = 3.28125 * std::sqrt(10.0);
    const auto f_86 = 19.6875 * std::sqrt(10.0);
    const auto f_87 = 0.328125 * std::sqrt(165.0);
    const auto f_88 = 4.921875 * std::sqrt(165.0);
    const auto f_89 = 0.109375 * std::sqrt(165.0);
    const auto f_90 = 1.640625 * std::sqrt(165.0);
    const auto f_91 = 0.09375 * std::sqrt(2310.0);
    const auto f_92 = 0.3125 * std::sqrt(2310.0);
    const auto f_93 = 0.5625 * std::sqrt(2310.0);
    const auto f_94 = 1.875 * std::sqrt(2310.0);
    const auto f_95 = 0.46875 * std::sqrt(770.0);
    const auto f_96 = 0.9375 * std::sqrt(770.0);
    const auto f_97 = 0.09375 * std::sqrt(770.0);
    const auto f_98 = 2.8125 * std::sqrt(770.0);
    const auto f_99 = 5.625 * std::sqrt(770.0);
    const auto f_100 = 0.5625 * std::sqrt(770.0);
    const auto f_101 = 0.375 * std::sqrt(35.0);
    const auto f_102 = 2.25 * std::sqrt(35.0);
    const auto f_103 = 22.5 * std::sqrt(35.0);
    const auto f_104 = 1.40625 * std::sqrt(42.0);
    const auto f_105 = 0.9375 * std::sqrt(42.0);
    const auto f_106 = 3.75 * std::sqrt(42.0);
    const auto f_107 = 0.46875 * std::sqrt(42.0);
    const auto f_108 = 1.25 * std::sqrt(42.0);
    const auto f_109 = 8.4375 * std::sqrt(42.0);
    const auto f_110 = 5.625 * std::sqrt(42.0);
    const auto f_111 = 22.5 * std::sqrt(42.0);
    const auto f_112 = 2.8125 * std::sqrt(42.0);
    const auto f_113 = 7.5 * std::sqrt(42.0);
    const auto f_114 = 0.15625 * std::sqrt(42.0);
    const auto f_115 = 0.3125 * std::sqrt(42.0);
    const auto f_116 = 2.5 * std::sqrt(42.0);
    const auto f_117 = 1.875 * std::sqrt(42.0);
    const auto f_118 = 15.0 * std::sqrt(42.0);
    const auto f_119 = 0.3125 * std::sqrt(105.0);
    const auto f_120 = 0.625 * std::sqrt(105.0);
    const auto f_121 = 1.25 * std::sqrt(105.0);
    const auto f_122 = 0.5 * std::sqrt(105.0);
    const auto f_123 = 1.875 * std::sqrt(105.0);
    const auto f_124 = 3.75 * std::sqrt(105.0);
    const auto f_125 = 7.5 * std::sqrt(105.0);
    const auto f_126 = 3.0 * std::sqrt(105.0);
    const auto f_127 = 0.15625 * std::sqrt(5.0);
    const auto f_128 = 0.46875 * std::sqrt(5.0);
    const auto f_129 = 2.8125 * std::sqrt(5.0);
    const auto f_130 = 5.625 * std::sqrt(5.0);
    const auto f_131 = 3.75 * std::sqrt(5.0);
    const auto f_132 = 0.5 * std::sqrt(5.0);
    const auto f_133 = 0.9375 * std::sqrt(5.0);
    const auto f_134 = 16.875 * std::sqrt(5.0);
    const auto f_135 = 33.75 * std::sqrt(5.0);
    const auto f_136 = 22.5 * std::sqrt(5.0);
    const auto f_137 = 3.0 * std::sqrt(5.0);
    const auto f_138 = 0.078125 * std::sqrt(42.0);
    const auto f_139 = 0.09375 * std::sqrt(35.0);
    const auto f_140 = 0.9375 * std::sqrt(35.0);
    const auto f_141 = 0.5625 * std::sqrt(35.0);
    const auto f_142 = 33.75 * std::sqrt(35.0);
    const auto f_143 = 0.015625 * std::sqrt(2310.0);
    const auto f_144 = 0.234375 * std::sqrt(2310.0);
    const auto f_145 = 1.40625 * std::sqrt(2310.0);
    const auto f_146 = 0.28125 * std::sqrt(1155.0);
    const auto f_147 = 0.9375 * std::sqrt(1155.0);
    const auto f_148 = 0.375 * std::sqrt(1155.0);
    const auto f_149 = 1.25 * std::sqrt(1155.0);
    const auto f_150 = 1.40625 * std::sqrt(385.0);
    const auto f_151 = 2.8125 * std::sqrt(385.0);
    const auto f_152 = 0.28125 * std::sqrt(385.0);
    const auto f_153 = 1.875 * std::sqrt(385.0);
    const auto f_154 = 3.75 * std::sqrt(385.0);
    const auto f_155 = 0.375 * std::sqrt(385.0);
    const auto f_156 = 0.5625 * std::sqrt(70.0);
    const auto f_157 = 7.5 * std::sqrt(70.0);
    const auto f_158 = 4.21875 * std::sqrt(21.0);
    const auto f_159 = 2.8125 * std::sqrt(21.0);
    const auto f_160 = 11.25 * std::sqrt(21.0);
    const auto f_161 = 1.40625 * std::sqrt(21.0);
    const auto f_162 = 3.75 * std::sqrt(21.0);
    const auto f_163 = 5.625 * std::sqrt(21.0);
    const auto f_164 = 15.0 * std::sqrt(21.0);
    const auto f_165 = 1.875 * std::sqrt(21.0);
    const auto f_166 = 5.0 * std::sqrt(21.0);
    const auto f_167 = 0.46875 * std::sqrt(21.0);
    const auto f_168 = 0.9375 * std::sqrt(21.0);
    const auto f_169 = 7.5 * std::sqrt(21.0);
    const auto f_170 = 0.625 * std::sqrt(21.0);
    const auto f_171 = 1.25 * std::sqrt(21.0);
    const auto f_172 = 10.0 * std::sqrt(21.0);
    const auto f_173 = 0.46875 * std::sqrt(210.0);
    const auto f_174 = 0.9375 * std::sqrt(210.0);
    const auto f_175 = 1.875 * std::sqrt(210.0);
    const auto f_176 = 0.75 * std::sqrt(210.0);
    const auto f_177 = 0.625 * std::sqrt(210.0);
    const auto f_178 = 1.25 * std::sqrt(210.0);
    const auto f_179 = 2.5 * std::sqrt(210.0);
    const auto f_180 = std::sqrt(210.0);
    const auto f_181 = 0.234375 * std::sqrt(10.0);
    const auto f_182 = 0.703125 * std::sqrt(10.0);
    const auto f_183 = 4.21875 * std::sqrt(10.0);
    const auto f_184 = 8.4375 * std::sqrt(10.0);
    const auto f_185 = 5.625 * std::sqrt(10.0);
    const auto f_186 = 0.75 * std::sqrt(10.0);
    const auto f_187 = 0.3125 * std::sqrt(10.0);
    const auto f_188 = 0.9375 * std::sqrt(10.0);
    const auto f_189 = 11.25 * std::sqrt(10.0);
    const auto f_190 = 7.5 * std::sqrt(10.0);
    const auto f_191 = std::sqrt(10.0);
    const auto f_192 = 0.234375 * std::sqrt(21.0);
    const auto f_193 = 0.3125 * std::sqrt(21.0);
    const auto f_194 = 0.140625 * std::sqrt(70.0);
    const auto f_195 = 0.1875 * std::sqrt(70.0);
    const auto f_196 = 0.9375 * std::sqrt(70.0);
    const auto f_197 = 11.25 * std::sqrt(70.0);
    const auto f_198 = 0.046875 * std::sqrt(1155.0);
    const auto f_199 = 0.703125 * std::sqrt(1155.0);
    const auto f_200 = 0.0625 * std::sqrt(1155.0);
    const auto f_201 = 0.0703125 * std::sqrt(462.0);
    const auto f_202 = 0.234375 * std::sqrt(462.0);
    const auto f_203 = 0.140625 * std::sqrt(462.0);
    const auto f_204 = 0.46875 * std::sqrt(462.0);
    const auto f_205 = 0.5625 * std::sqrt(462.0);
    const auto f_206 = 1.875 * std::sqrt(462.0);
    const auto f_207 = 0.1875 * std::sqrt(462.0);
    const auto f_208 = 0.625 * std::sqrt(462.0);
    const auto f_209 = 0.3515625 * std::sqrt(154.0);
    const auto f_210 = 0.703125 * std::sqrt(154.0);
    const auto f_211 = 0.0703125 * std::sqrt(154.0);
    const auto f_212 = 1.40625 * std::sqrt(154.0);
    const auto f_213 = 0.140625 * std::sqrt(154.0);
    const auto f_214 = 2.8125 * std::sqrt(154.0);
    const auto f_215 = 5.625 * std::sqrt(154.0);
    const auto f_216 = 0.5625 * std::sqrt(154.0);
    const auto f_217 = 0.9375 * std::sqrt(154.0);
    const auto f_218 = 1.875 * std::sqrt(154.0);
    const auto f_219 = 0.1875 * std::sqrt(154.0);
    const auto f_220 = 0.28125 * std::sqrt(7.0);
    const auto f_221 = 2.8125 * std::sqrt(7.0);
    const auto f_222 = 0.5625 * std::sqrt(7.0);
    const auto f_223 = 5.625 * std::sqrt(7.0);
    const auto f_224 = 2.25 * std::sqrt(7.0);
    const auto f_225 = 22.5 * std::sqrt(7.0);
    const auto f_226 = 0.75 * std::sqrt(7.0);
    const auto f_227 = 7.5 * std::sqrt(7.0);
    const auto f_228 = 0.2109375 * std::sqrt(210.0);
    const auto f_229 = 0.140625 * std::sqrt(210.0);
    const auto f_230 = 0.5625 * std::sqrt(210.0);
    const auto f_231 = 0.0703125 * std::sqrt(210.0);
    const auto f_232 = 0.1875 * std::sqrt(210.0);
    const auto f_233 = 0.421875 * std::sqrt(210.0);
    const auto f_234 = 0.28125 * std::sqrt(210.0);
    const auto f_235 = 1.125 * std::sqrt(210.0);
    const auto f_236 = 0.375 * std::sqrt(210.0);
    const auto f_237 = 1.6875 * std::sqrt(210.0);
    const auto f_238 = 4.5 * std::sqrt(210.0);
    const auto f_239 = 1.5 * std::sqrt(210.0);
    const auto f_240 = 0.5 * std::sqrt(210.0);
    const auto f_241 = 0.0234375 * std::sqrt(210.0);
    const auto f_242 = 0.046875 * std::sqrt(210.0);
    const auto f_243 = 0.09375 * std::sqrt(210.0);
    const auto f_244 = 3.0 * std::sqrt(210.0);
    const auto f_245 = 0.0625 * std::sqrt(210.0);
    const auto f_246 = 0.125 * std::sqrt(210.0);
    const auto f_247 = 0.375 * std::sqrt(21.0);
    const auto f_248 = 0.75 * std::sqrt(21.0);
    const auto f_249 = 3.0 * std::sqrt(21.0);
    const auto f_250 = 2.5 * std::sqrt(21.0);
    const auto f_251 = std::sqrt(21.0);
    const auto f_252 = 0.01171875 * std::sqrt(210.0);
    const auto f_253 = 0.03125 * std::sqrt(210.0);
    const auto f_254 = 0.0703125 * std::sqrt(7.0);
    const auto f_255 = 0.3515625 * std::sqrt(7.0);
    const auto f_256 = 0.703125 * std::sqrt(7.0);
    const auto f_257 = 4.21875 * std::sqrt(7.0);
    const auto f_258 = 0.140625 * std::sqrt(7.0);
    const auto f_259 = 1.40625 * std::sqrt(7.0);
    const auto f_260 = 8.4375 * std::sqrt(7.0);
    const auto f_261 = 33.75 * std::sqrt(7.0);
    const auto f_262 = 0.1875 * std::sqrt(7.0);
    const auto f_263 = 0.9375 * std::sqrt(7.0);
    const auto f_264 = 1.875 * std::sqrt(7.0);
    const auto f_265 = 11.25 * std::sqrt(7.0);
    const auto f_266 = 0.01171875 * std::sqrt(462.0);
    const auto f_267 = 0.17578125 * std::sqrt(462.0);
    const auto f_268 = 0.0234375 * std::sqrt(462.0);
    const auto f_269 = 0.3515625 * std::sqrt(462.0);
    const auto f_270 = 0.09375 * std::sqrt(462.0);
    const auto f_271 = 1.40625 * std::sqrt(462.0);
    const auto f_272 = 0.03125 * std::sqrt(462.0);
    const auto f_273 = 0.046875 * std::sqrt(2310.0);
    const auto f_274 = 0.15625 * std::sqrt(2310.0);
    const auto f_275 = 0.28125 * std::sqrt(2310.0);
    const auto f_276 = 0.9375 * std::sqrt(2310.0);
    const auto f_277 = 0.234375 * std::sqrt(770.0);
    const auto f_278 = 0.046875 * std::sqrt(770.0);
    const auto f_279 = 1.40625 * std::sqrt(770.0);
    const auto f_280 = 0.28125 * std::sqrt(770.0);
    const auto f_281 = 0.1875 * std::sqrt(35.0);
    const auto f_282 = 1.875 * std::sqrt(35.0);
    const auto f_283 = 1.125 * std::sqrt(35.0);
    const auto f_284 = 11.25 * std::sqrt(35.0);
    const auto f_285 = 0.703125 * std::sqrt(42.0);
    const auto f_286 = 0.234375 * std::sqrt(42.0);
    const auto f_287 = 0.625 * std::sqrt(42.0);
    const auto f_288 = 4.21875 * std::sqrt(42.0);
    const auto f_289 = 11.25 * std::sqrt(42.0);
    const auto f_290 = 0.15625 * std::sqrt(105.0);
    const auto f_291 = 0.25 * std::sqrt(105.0);
    const auto f_292 = 0.9375 * std::sqrt(105.0);
    const auto f_293 = 1.5 * std::sqrt(105.0);
    const auto f_294 = 0.078125 * std::sqrt(5.0);
    const auto f_295 = 0.234375 * std::sqrt(5.0);
    const auto f_296 = 1.40625 * std::sqrt(5.0);
    const auto f_297 = 1.875 * std::sqrt(5.0);
    const auto f_298 = 0.25 * std::sqrt(5.0);
    const auto f_299 = 8.4375 * std::sqrt(5.0);
    const auto f_300 = 11.25 * std::sqrt(5.0);
    const auto f_301 = 1.5 * std::sqrt(5.0);
    const auto f_302 = 0.0390625 * std::sqrt(42.0);
    const auto f_303 = 0.046875 * std::sqrt(35.0);
    const auto f_304 = 0.234375 * std::sqrt(35.0);
    const auto f_305 = 0.28125 * std::sqrt(35.0);
    const auto f_306 = 1.40625 * std::sqrt(35.0);
    const auto f_307 = 16.875 * std::sqrt(35.0);
    const auto f_308 = 0.0078125 * std::sqrt(2310.0);
    const auto f_309 = 0.1171875 * std::sqrt(2310.0);
    const auto f_310 = 0.703125 * std::sqrt(2310.0);
    const auto f_311 = 0.1640625 * std::sqrt(330.0);
    const auto f_312 = 0.546875 * std::sqrt(330.0);
    const auto f_313 = 0.984375 * std::sqrt(330.0);
    const auto f_314 = 3.28125 * std::sqrt(330.0);
    const auto f_315 = 0.8203125 * std::sqrt(110.0);
    const auto f_316 = 1.640625 * std::sqrt(110.0);
    const auto f_317 = 0.1640625 * std::sqrt(110.0);
    const auto f_318 = 4.921875 * std::sqrt(110.0);
    const auto f_319 = 9.84375 * std::sqrt(110.0);
    const auto f_320 = 0.984375 * std::sqrt(110.0);
    const auto f_321 = 3.9375 * std::sqrt(5.0);
    const auto f_322 = 2.4609375 * std::sqrt(6.0);
    const auto f_323 = 1.640625 * std::sqrt(6.0);
    const auto f_324 = 0.8203125 * std::sqrt(6.0);
    const auto f_325 = 14.765625 * std::sqrt(6.0);
    const auto f_326 = 39.375 * std::sqrt(6.0);
    const auto f_327 = 4.921875 * std::sqrt(6.0);
    const auto f_328 = 13.125 * std::sqrt(6.0);
    const auto f_329 = 0.2734375 * std::sqrt(6.0);
    const auto f_330 = 4.375 * std::sqrt(6.0);
    const auto f_331 = 0.546875 * std::sqrt(15.0);
    const auto f_332 = 1.09375 * std::sqrt(15.0);
    const auto f_333 = 0.875 * std::sqrt(15.0);
    const auto f_334 = 3.28125 * std::sqrt(15.0);
    const auto f_335 = 6.5625 * std::sqrt(15.0);
    const auto f_336 = 13.125 * std::sqrt(15.0);
    const auto f_337 = 5.25 * std::sqrt(15.0);
    const auto f_338 = 0.0390625 * std::sqrt(35.0);
    const auto f_339 = 0.1171875 * std::sqrt(35.0);
    const auto f_340 = 0.703125 * std::sqrt(35.0);
    const auto f_341 = 0.125 * std::sqrt(35.0);
    const auto f_342 = 4.21875 * std::sqrt(35.0);
    const auto f_343 = 8.4375 * std::sqrt(35.0);
    const auto f_344 = 0.75 * std::sqrt(35.0);
    const auto f_345 = 0.13671875 * std::sqrt(6.0);
    const auto f_346 = 0.1640625 * std::sqrt(5.0);
    const auto f_347 = 0.8203125 * std::sqrt(5.0);
    const auto f_348 = 1.640625 * std::sqrt(5.0);
    const auto f_349 = 9.84375 * std::sqrt(5.0);
    const auto f_350 = 0.984375 * std::sqrt(5.0);
    const auto f_351 = 4.921875 * std::sqrt(5.0);
    const auto f_352 = 59.0625 * std::sqrt(5.0);
    const auto f_353 = 0.02734375 * std::sqrt(330.0);
    const auto f_354 = 0.41015625 * std::sqrt(330.0);
    const auto f_355 = 2.4609375 * std::sqrt(330.0);

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
    auto *g_105 = values + 105 * nvalues;
    auto *g_106 = values + 106 * nvalues;
    auto *g_107 = values + 107 * nvalues;
    auto *g_108 = values + 108 * nvalues;
    auto *g_109 = values + 109 * nvalues;
    auto *g_110 = values + 110 * nvalues;
    auto *g_111 = values + 111 * nvalues;
    auto *g_112 = values + 112 * nvalues;
    auto *g_113 = values + 113 * nvalues;
    auto *g_114 = values + 114 * nvalues;
    auto *g_115 = values + 115 * nvalues;
    auto *g_116 = values + 116 * nvalues;

    const auto *gi_0 = buffer.data(gi + 0);
    const auto *gi_1 = buffer.data(gi + 1);
    const auto *gi_2 = buffer.data(gi + 2);
    const auto *gi_3 = buffer.data(gi + 3);
    const auto *gi_4 = buffer.data(gi + 4);
    const auto *gi_5 = buffer.data(gi + 5);
    const auto *gi_6 = buffer.data(gi + 6);
    const auto *gi_7 = buffer.data(gi + 7);
    const auto *gi_8 = buffer.data(gi + 8);
    const auto *gi_9 = buffer.data(gi + 9);
    const auto *gi_10 = buffer.data(gi + 10);
    const auto *gi_11 = buffer.data(gi + 11);
    const auto *gi_12 = buffer.data(gi + 12);
    const auto *gi_13 = buffer.data(gi + 13);
    const auto *gi_14 = buffer.data(gi + 14);
    const auto *gi_15 = buffer.data(gi + 15);
    const auto *gi_16 = buffer.data(gi + 16);
    const auto *gi_17 = buffer.data(gi + 17);
    const auto *gi_18 = buffer.data(gi + 18);
    const auto *gi_19 = buffer.data(gi + 19);
    const auto *gi_20 = buffer.data(gi + 20);
    const auto *gi_21 = buffer.data(gi + 21);
    const auto *gi_22 = buffer.data(gi + 22);
    const auto *gi_23 = buffer.data(gi + 23);
    const auto *gi_24 = buffer.data(gi + 24);
    const auto *gi_25 = buffer.data(gi + 25);
    const auto *gi_26 = buffer.data(gi + 26);
    const auto *gi_27 = buffer.data(gi + 27);
    const auto *gi_28 = buffer.data(gi + 28);
    const auto *gi_29 = buffer.data(gi + 29);
    const auto *gi_30 = buffer.data(gi + 30);
    const auto *gi_31 = buffer.data(gi + 31);
    const auto *gi_32 = buffer.data(gi + 32);
    const auto *gi_33 = buffer.data(gi + 33);
    const auto *gi_34 = buffer.data(gi + 34);
    const auto *gi_35 = buffer.data(gi + 35);
    const auto *gi_36 = buffer.data(gi + 36);
    const auto *gi_37 = buffer.data(gi + 37);
    const auto *gi_38 = buffer.data(gi + 38);
    const auto *gi_39 = buffer.data(gi + 39);
    const auto *gi_40 = buffer.data(gi + 40);
    const auto *gi_41 = buffer.data(gi + 41);
    const auto *gi_42 = buffer.data(gi + 42);
    const auto *gi_43 = buffer.data(gi + 43);
    const auto *gi_44 = buffer.data(gi + 44);
    const auto *gi_45 = buffer.data(gi + 45);
    const auto *gi_46 = buffer.data(gi + 46);
    const auto *gi_47 = buffer.data(gi + 47);
    const auto *gi_48 = buffer.data(gi + 48);
    const auto *gi_49 = buffer.data(gi + 49);
    const auto *gi_50 = buffer.data(gi + 50);
    const auto *gi_51 = buffer.data(gi + 51);
    const auto *gi_52 = buffer.data(gi + 52);
    const auto *gi_53 = buffer.data(gi + 53);
    const auto *gi_54 = buffer.data(gi + 54);
    const auto *gi_55 = buffer.data(gi + 55);
    const auto *gi_56 = buffer.data(gi + 56);
    const auto *gi_57 = buffer.data(gi + 57);
    const auto *gi_58 = buffer.data(gi + 58);
    const auto *gi_59 = buffer.data(gi + 59);
    const auto *gi_60 = buffer.data(gi + 60);
    const auto *gi_61 = buffer.data(gi + 61);
    const auto *gi_62 = buffer.data(gi + 62);
    const auto *gi_63 = buffer.data(gi + 63);
    const auto *gi_64 = buffer.data(gi + 64);
    const auto *gi_65 = buffer.data(gi + 65);
    const auto *gi_66 = buffer.data(gi + 66);
    const auto *gi_67 = buffer.data(gi + 67);
    const auto *gi_68 = buffer.data(gi + 68);
    const auto *gi_69 = buffer.data(gi + 69);
    const auto *gi_70 = buffer.data(gi + 70);
    const auto *gi_71 = buffer.data(gi + 71);
    const auto *gi_72 = buffer.data(gi + 72);
    const auto *gi_73 = buffer.data(gi + 73);
    const auto *gi_74 = buffer.data(gi + 74);
    const auto *gi_75 = buffer.data(gi + 75);
    const auto *gi_76 = buffer.data(gi + 76);
    const auto *gi_77 = buffer.data(gi + 77);
    const auto *gi_78 = buffer.data(gi + 78);
    const auto *gi_79 = buffer.data(gi + 79);
    const auto *gi_80 = buffer.data(gi + 80);
    const auto *gi_81 = buffer.data(gi + 81);
    const auto *gi_82 = buffer.data(gi + 82);
    const auto *gi_83 = buffer.data(gi + 83);
    const auto *gi_84 = buffer.data(gi + 84);
    const auto *gi_85 = buffer.data(gi + 85);
    const auto *gi_86 = buffer.data(gi + 86);
    const auto *gi_87 = buffer.data(gi + 87);
    const auto *gi_88 = buffer.data(gi + 88);
    const auto *gi_89 = buffer.data(gi + 89);
    const auto *gi_90 = buffer.data(gi + 90);
    const auto *gi_91 = buffer.data(gi + 91);
    const auto *gi_92 = buffer.data(gi + 92);
    const auto *gi_93 = buffer.data(gi + 93);
    const auto *gi_94 = buffer.data(gi + 94);
    const auto *gi_95 = buffer.data(gi + 95);
    const auto *gi_96 = buffer.data(gi + 96);
    const auto *gi_97 = buffer.data(gi + 97);
    const auto *gi_98 = buffer.data(gi + 98);
    const auto *gi_99 = buffer.data(gi + 99);
    const auto *gi_100 = buffer.data(gi + 100);
    const auto *gi_101 = buffer.data(gi + 101);
    const auto *gi_102 = buffer.data(gi + 102);
    const auto *gi_103 = buffer.data(gi + 103);
    const auto *gi_104 = buffer.data(gi + 104);
    const auto *gi_105 = buffer.data(gi + 105);
    const auto *gi_106 = buffer.data(gi + 106);
    const auto *gi_107 = buffer.data(gi + 107);
    const auto *gi_108 = buffer.data(gi + 108);
    const auto *gi_109 = buffer.data(gi + 109);
    const auto *gi_110 = buffer.data(gi + 110);
    const auto *gi_111 = buffer.data(gi + 111);
    const auto *gi_112 = buffer.data(gi + 112);
    const auto *gi_113 = buffer.data(gi + 113);
    const auto *gi_114 = buffer.data(gi + 114);
    const auto *gi_115 = buffer.data(gi + 115);
    const auto *gi_116 = buffer.data(gi + 116);
    const auto *gi_117 = buffer.data(gi + 117);
    const auto *gi_118 = buffer.data(gi + 118);
    const auto *gi_119 = buffer.data(gi + 119);
    const auto *gi_120 = buffer.data(gi + 120);
    const auto *gi_121 = buffer.data(gi + 121);
    const auto *gi_122 = buffer.data(gi + 122);
    const auto *gi_123 = buffer.data(gi + 123);
    const auto *gi_124 = buffer.data(gi + 124);
    const auto *gi_125 = buffer.data(gi + 125);
    const auto *gi_126 = buffer.data(gi + 126);
    const auto *gi_127 = buffer.data(gi + 127);
    const auto *gi_128 = buffer.data(gi + 128);
    const auto *gi_129 = buffer.data(gi + 129);
    const auto *gi_130 = buffer.data(gi + 130);
    const auto *gi_131 = buffer.data(gi + 131);
    const auto *gi_132 = buffer.data(gi + 132);
    const auto *gi_133 = buffer.data(gi + 133);
    const auto *gi_134 = buffer.data(gi + 134);
    const auto *gi_135 = buffer.data(gi + 135);
    const auto *gi_136 = buffer.data(gi + 136);
    const auto *gi_137 = buffer.data(gi + 137);
    const auto *gi_138 = buffer.data(gi + 138);
    const auto *gi_139 = buffer.data(gi + 139);
    const auto *gi_140 = buffer.data(gi + 140);
    const auto *gi_141 = buffer.data(gi + 141);
    const auto *gi_142 = buffer.data(gi + 142);
    const auto *gi_143 = buffer.data(gi + 143);
    const auto *gi_144 = buffer.data(gi + 144);
    const auto *gi_145 = buffer.data(gi + 145);
    const auto *gi_146 = buffer.data(gi + 146);
    const auto *gi_147 = buffer.data(gi + 147);
    const auto *gi_148 = buffer.data(gi + 148);
    const auto *gi_149 = buffer.data(gi + 149);
    const auto *gi_150 = buffer.data(gi + 150);
    const auto *gi_151 = buffer.data(gi + 151);
    const auto *gi_152 = buffer.data(gi + 152);
    const auto *gi_153 = buffer.data(gi + 153);
    const auto *gi_154 = buffer.data(gi + 154);
    const auto *gi_155 = buffer.data(gi + 155);
    const auto *gi_156 = buffer.data(gi + 156);
    const auto *gi_157 = buffer.data(gi + 157);
    const auto *gi_158 = buffer.data(gi + 158);
    const auto *gi_159 = buffer.data(gi + 159);
    const auto *gi_160 = buffer.data(gi + 160);
    const auto *gi_161 = buffer.data(gi + 161);
    const auto *gi_162 = buffer.data(gi + 162);
    const auto *gi_163 = buffer.data(gi + 163);
    const auto *gi_164 = buffer.data(gi + 164);
    const auto *gi_165 = buffer.data(gi + 165);
    const auto *gi_166 = buffer.data(gi + 166);
    const auto *gi_167 = buffer.data(gi + 167);
    const auto *gi_168 = buffer.data(gi + 168);
    const auto *gi_169 = buffer.data(gi + 169);
    const auto *gi_170 = buffer.data(gi + 170);
    const auto *gi_171 = buffer.data(gi + 171);
    const auto *gi_172 = buffer.data(gi + 172);
    const auto *gi_173 = buffer.data(gi + 173);
    const auto *gi_174 = buffer.data(gi + 174);
    const auto *gi_175 = buffer.data(gi + 175);
    const auto *gi_176 = buffer.data(gi + 176);
    const auto *gi_177 = buffer.data(gi + 177);
    const auto *gi_178 = buffer.data(gi + 178);
    const auto *gi_179 = buffer.data(gi + 179);
    const auto *gi_180 = buffer.data(gi + 180);
    const auto *gi_181 = buffer.data(gi + 181);
    const auto *gi_182 = buffer.data(gi + 182);
    const auto *gi_183 = buffer.data(gi + 183);
    const auto *gi_184 = buffer.data(gi + 184);
    const auto *gi_185 = buffer.data(gi + 185);
    const auto *gi_186 = buffer.data(gi + 186);
    const auto *gi_187 = buffer.data(gi + 187);
    const auto *gi_188 = buffer.data(gi + 188);
    const auto *gi_189 = buffer.data(gi + 189);
    const auto *gi_190 = buffer.data(gi + 190);
    const auto *gi_191 = buffer.data(gi + 191);
    const auto *gi_192 = buffer.data(gi + 192);
    const auto *gi_193 = buffer.data(gi + 193);
    const auto *gi_194 = buffer.data(gi + 194);
    const auto *gi_195 = buffer.data(gi + 195);
    const auto *gi_196 = buffer.data(gi + 196);
    const auto *gi_197 = buffer.data(gi + 197);
    const auto *gi_198 = buffer.data(gi + 198);
    const auto *gi_199 = buffer.data(gi + 199);
    const auto *gi_200 = buffer.data(gi + 200);
    const auto *gi_201 = buffer.data(gi + 201);
    const auto *gi_202 = buffer.data(gi + 202);
    const auto *gi_203 = buffer.data(gi + 203);
    const auto *gi_204 = buffer.data(gi + 204);
    const auto *gi_205 = buffer.data(gi + 205);
    const auto *gi_206 = buffer.data(gi + 206);
    const auto *gi_207 = buffer.data(gi + 207);
    const auto *gi_208 = buffer.data(gi + 208);
    const auto *gi_209 = buffer.data(gi + 209);
    const auto *gi_210 = buffer.data(gi + 210);
    const auto *gi_211 = buffer.data(gi + 211);
    const auto *gi_212 = buffer.data(gi + 212);
    const auto *gi_213 = buffer.data(gi + 213);
    const auto *gi_214 = buffer.data(gi + 214);
    const auto *gi_215 = buffer.data(gi + 215);
    const auto *gi_216 = buffer.data(gi + 216);
    const auto *gi_217 = buffer.data(gi + 217);
    const auto *gi_218 = buffer.data(gi + 218);
    const auto *gi_219 = buffer.data(gi + 219);
    const auto *gi_220 = buffer.data(gi + 220);
    const auto *gi_221 = buffer.data(gi + 221);
    const auto *gi_222 = buffer.data(gi + 222);
    const auto *gi_223 = buffer.data(gi + 223);
    const auto *gi_224 = buffer.data(gi + 224);
    const auto *gi_225 = buffer.data(gi + 225);
    const auto *gi_226 = buffer.data(gi + 226);
    const auto *gi_227 = buffer.data(gi + 227);
    const auto *gi_228 = buffer.data(gi + 228);
    const auto *gi_229 = buffer.data(gi + 229);
    const auto *gi_230 = buffer.data(gi + 230);
    const auto *gi_231 = buffer.data(gi + 231);
    const auto *gi_232 = buffer.data(gi + 232);
    const auto *gi_233 = buffer.data(gi + 233);
    const auto *gi_234 = buffer.data(gi + 234);
    const auto *gi_235 = buffer.data(gi + 235);
    const auto *gi_236 = buffer.data(gi + 236);
    const auto *gi_237 = buffer.data(gi + 237);
    const auto *gi_238 = buffer.data(gi + 238);
    const auto *gi_239 = buffer.data(gi + 239);
    const auto *gi_240 = buffer.data(gi + 240);
    const auto *gi_241 = buffer.data(gi + 241);
    const auto *gi_242 = buffer.data(gi + 242);
    const auto *gi_243 = buffer.data(gi + 243);
    const auto *gi_244 = buffer.data(gi + 244);
    const auto *gi_245 = buffer.data(gi + 245);
    const auto *gi_246 = buffer.data(gi + 246);
    const auto *gi_247 = buffer.data(gi + 247);
    const auto *gi_248 = buffer.data(gi + 248);
    const auto *gi_249 = buffer.data(gi + 249);
    const auto *gi_250 = buffer.data(gi + 250);
    const auto *gi_251 = buffer.data(gi + 251);
    const auto *gi_252 = buffer.data(gi + 252);
    const auto *gi_253 = buffer.data(gi + 253);
    const auto *gi_254 = buffer.data(gi + 254);
    const auto *gi_255 = buffer.data(gi + 255);
    const auto *gi_256 = buffer.data(gi + 256);
    const auto *gi_257 = buffer.data(gi + 257);
    const auto *gi_258 = buffer.data(gi + 258);
    const auto *gi_259 = buffer.data(gi + 259);
    const auto *gi_260 = buffer.data(gi + 260);
    const auto *gi_261 = buffer.data(gi + 261);
    const auto *gi_262 = buffer.data(gi + 262);
    const auto *gi_263 = buffer.data(gi + 263);
    const auto *gi_264 = buffer.data(gi + 264);
    const auto *gi_265 = buffer.data(gi + 265);
    const auto *gi_266 = buffer.data(gi + 266);
    const auto *gi_267 = buffer.data(gi + 267);
    const auto *gi_268 = buffer.data(gi + 268);
    const auto *gi_269 = buffer.data(gi + 269);
    const auto *gi_270 = buffer.data(gi + 270);
    const auto *gi_271 = buffer.data(gi + 271);
    const auto *gi_272 = buffer.data(gi + 272);
    const auto *gi_273 = buffer.data(gi + 273);
    const auto *gi_274 = buffer.data(gi + 274);
    const auto *gi_275 = buffer.data(gi + 275);
    const auto *gi_276 = buffer.data(gi + 276);
    const auto *gi_277 = buffer.data(gi + 277);
    const auto *gi_278 = buffer.data(gi + 278);
    const auto *gi_279 = buffer.data(gi + 279);
    const auto *gi_280 = buffer.data(gi + 280);
    const auto *gi_281 = buffer.data(gi + 281);
    const auto *gi_282 = buffer.data(gi + 282);
    const auto *gi_283 = buffer.data(gi + 283);
    const auto *gi_284 = buffer.data(gi + 284);
    const auto *gi_285 = buffer.data(gi + 285);
    const auto *gi_286 = buffer.data(gi + 286);
    const auto *gi_287 = buffer.data(gi + 287);
    const auto *gi_288 = buffer.data(gi + 288);
    const auto *gi_289 = buffer.data(gi + 289);
    const auto *gi_290 = buffer.data(gi + 290);
    const auto *gi_291 = buffer.data(gi + 291);
    const auto *gi_292 = buffer.data(gi + 292);
    const auto *gi_293 = buffer.data(gi + 293);
    const auto *gi_294 = buffer.data(gi + 294);
    const auto *gi_295 = buffer.data(gi + 295);
    const auto *gi_296 = buffer.data(gi + 296);
    const auto *gi_297 = buffer.data(gi + 297);
    const auto *gi_298 = buffer.data(gi + 298);
    const auto *gi_299 = buffer.data(gi + 299);
    const auto *gi_300 = buffer.data(gi + 300);
    const auto *gi_301 = buffer.data(gi + 301);
    const auto *gi_302 = buffer.data(gi + 302);
    const auto *gi_303 = buffer.data(gi + 303);
    const auto *gi_304 = buffer.data(gi + 304);
    const auto *gi_305 = buffer.data(gi + 305);
    const auto *gi_306 = buffer.data(gi + 306);
    const auto *gi_307 = buffer.data(gi + 307);
    const auto *gi_308 = buffer.data(gi + 308);
    const auto *gi_309 = buffer.data(gi + 309);
    const auto *gi_310 = buffer.data(gi + 310);
    const auto *gi_311 = buffer.data(gi + 311);
    const auto *gi_312 = buffer.data(gi + 312);
    const auto *gi_313 = buffer.data(gi + 313);
    const auto *gi_314 = buffer.data(gi + 314);
    const auto *gi_315 = buffer.data(gi + 315);
    const auto *gi_316 = buffer.data(gi + 316);
    const auto *gi_317 = buffer.data(gi + 317);
    const auto *gi_318 = buffer.data(gi + 318);
    const auto *gi_319 = buffer.data(gi + 319);
    const auto *gi_320 = buffer.data(gi + 320);
    const auto *gi_321 = buffer.data(gi + 321);
    const auto *gi_322 = buffer.data(gi + 322);
    const auto *gi_323 = buffer.data(gi + 323);
    const auto *gi_324 = buffer.data(gi + 324);
    const auto *gi_325 = buffer.data(gi + 325);
    const auto *gi_326 = buffer.data(gi + 326);
    const auto *gi_327 = buffer.data(gi + 327);
    const auto *gi_328 = buffer.data(gi + 328);
    const auto *gi_329 = buffer.data(gi + 329);
    const auto *gi_330 = buffer.data(gi + 330);
    const auto *gi_331 = buffer.data(gi + 331);
    const auto *gi_332 = buffer.data(gi + 332);
    const auto *gi_333 = buffer.data(gi + 333);
    const auto *gi_334 = buffer.data(gi + 334);
    const auto *gi_335 = buffer.data(gi + 335);
    const auto *gi_336 = buffer.data(gi + 336);
    const auto *gi_337 = buffer.data(gi + 337);
    const auto *gi_338 = buffer.data(gi + 338);
    const auto *gi_339 = buffer.data(gi + 339);
    const auto *gi_340 = buffer.data(gi + 340);
    const auto *gi_341 = buffer.data(gi + 341);
    const auto *gi_342 = buffer.data(gi + 342);
    const auto *gi_343 = buffer.data(gi + 343);
    const auto *gi_344 = buffer.data(gi + 344);
    const auto *gi_345 = buffer.data(gi + 345);
    const auto *gi_346 = buffer.data(gi + 346);
    const auto *gi_347 = buffer.data(gi + 347);
    const auto *gi_348 = buffer.data(gi + 348);
    const auto *gi_349 = buffer.data(gi + 349);
    const auto *gi_350 = buffer.data(gi + 350);
    const auto *gi_351 = buffer.data(gi + 351);
    const auto *gi_352 = buffer.data(gi + 352);
    const auto *gi_353 = buffer.data(gi + 353);
    const auto *gi_354 = buffer.data(gi + 354);
    const auto *gi_355 = buffer.data(gi + 355);
    const auto *gi_356 = buffer.data(gi + 356);
    const auto *gi_357 = buffer.data(gi + 357);
    const auto *gi_358 = buffer.data(gi + 358);
    const auto *gi_359 = buffer.data(gi + 359);
    const auto *gi_360 = buffer.data(gi + 360);
    const auto *gi_361 = buffer.data(gi + 361);
    const auto *gi_362 = buffer.data(gi + 362);
    const auto *gi_363 = buffer.data(gi + 363);
    const auto *gi_364 = buffer.data(gi + 364);
    const auto *gi_365 = buffer.data(gi + 365);
    const auto *gi_366 = buffer.data(gi + 366);
    const auto *gi_367 = buffer.data(gi + 367);
    const auto *gi_368 = buffer.data(gi + 368);
    const auto *gi_369 = buffer.data(gi + 369);
    const auto *gi_370 = buffer.data(gi + 370);
    const auto *gi_371 = buffer.data(gi + 371);
    const auto *gi_372 = buffer.data(gi + 372);
    const auto *gi_373 = buffer.data(gi + 373);
    const auto *gi_374 = buffer.data(gi + 374);
    const auto *gi_375 = buffer.data(gi + 375);
    const auto *gi_376 = buffer.data(gi + 376);
    const auto *gi_377 = buffer.data(gi + 377);
    const auto *gi_378 = buffer.data(gi + 378);
    const auto *gi_379 = buffer.data(gi + 379);
    const auto *gi_380 = buffer.data(gi + 380);
    const auto *gi_381 = buffer.data(gi + 381);
    const auto *gi_382 = buffer.data(gi + 382);
    const auto *gi_383 = buffer.data(gi + 383);
    const auto *gi_384 = buffer.data(gi + 384);
    const auto *gi_385 = buffer.data(gi + 385);
    const auto *gi_386 = buffer.data(gi + 386);
    const auto *gi_387 = buffer.data(gi + 387);
    const auto *gi_388 = buffer.data(gi + 388);
    const auto *gi_389 = buffer.data(gi + 389);
    const auto *gi_390 = buffer.data(gi + 390);
    const auto *gi_391 = buffer.data(gi + 391);
    const auto *gi_392 = buffer.data(gi + 392);
    const auto *gi_393 = buffer.data(gi + 393);
    const auto *gi_394 = buffer.data(gi + 394);
    const auto *gi_395 = buffer.data(gi + 395);
    const auto *gi_396 = buffer.data(gi + 396);
    const auto *gi_397 = buffer.data(gi + 397);
    const auto *gi_398 = buffer.data(gi + 398);
    const auto *gi_399 = buffer.data(gi + 399);
    const auto *gi_400 = buffer.data(gi + 400);
    const auto *gi_401 = buffer.data(gi + 401);
    const auto *gi_402 = buffer.data(gi + 402);
    const auto *gi_403 = buffer.data(gi + 403);
    const auto *gi_404 = buffer.data(gi + 404);
    const auto *gi_405 = buffer.data(gi + 405);
    const auto *gi_406 = buffer.data(gi + 406);
    const auto *gi_407 = buffer.data(gi + 407);
    const auto *gi_408 = buffer.data(gi + 408);
    const auto *gi_409 = buffer.data(gi + 409);
    const auto *gi_410 = buffer.data(gi + 410);
    const auto *gi_411 = buffer.data(gi + 411);
    const auto *gi_412 = buffer.data(gi + 412);
    const auto *gi_413 = buffer.data(gi + 413);
    const auto *gi_414 = buffer.data(gi + 414);
    const auto *gi_415 = buffer.data(gi + 415);
    const auto *gi_416 = buffer.data(gi + 416);
    const auto *gi_417 = buffer.data(gi + 417);
    const auto *gi_418 = buffer.data(gi + 418);
    const auto *gi_419 = buffer.data(gi + 419);

#pragma omp simd aligned(gi_29, gi_32, gi_34, gi_39, gi_43, gi_50, gi_169, gi_172, gi_174, \
                         gi_179, gi_183, gi_190 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * gi_29[k]
                 - f_1 * gi_34[k]
                 + f_0 * gi_43[k]
                 - f_0 * gi_169[k]
                 + f_1 * gi_174[k]
                 - f_0 * gi_183[k];

        g_1[k] = f_2 * gi_32[k]
                 - f_3 * gi_39[k]
                 + f_4 * gi_50[k]
                 - f_2 * gi_172[k]
                 + f_3 * gi_179[k]
                 - f_4 * gi_190[k];
    }

#pragma omp simd aligned(gi_29, gi_36, gi_43, gi_45, gi_169, gi_176, gi_183, \
                         gi_185 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_2[k] = -f_5 * gi_29[k]
                 + f_6 * gi_36[k]
                 + f_5 * gi_43[k]
                 - f_6 * gi_45[k]
                 + f_5 * gi_169[k]
                 - f_6 * gi_176[k]
                 - f_5 * gi_183[k]
                 + f_6 * gi_185[k];
    }

#pragma omp simd aligned(gi_32, gi_39, gi_41, gi_50, gi_52, gi_172, gi_179, gi_181, gi_190, \
                         gi_192 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = -f_7 * gi_32[k]
                 - f_8 * gi_39[k]
                 + f_9 * gi_41[k]
                 + f_10 * gi_50[k]
                 - f_11 * gi_52[k]
                 + f_7 * gi_172[k]
                 + f_8 * gi_179[k]
                 - f_9 * gi_181[k]
                 - f_10 * gi_190[k]
                 + f_11 * gi_192[k];
    }

#pragma omp simd aligned(gi_29, gi_34, gi_36, gi_43, gi_45, gi_47, gi_169, gi_174, gi_176, \
                         gi_183, gi_185, gi_187 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = f_12 * gi_29[k]
                 + f_13 * gi_34[k]
                 - f_14 * gi_36[k]
                 + f_12 * gi_43[k]
                 - f_14 * gi_45[k]
                 + f_14 * gi_47[k]
                 - f_12 * gi_169[k]
                 - f_13 * gi_174[k]
                 + f_14 * gi_176[k]
                 - f_12 * gi_183[k]
                 + f_14 * gi_185[k]
                 - f_14 * gi_187[k];
    }

#pragma omp simd aligned(gi_32, gi_39, gi_41, gi_50, gi_52, gi_54, gi_172, gi_179, gi_181, \
                         gi_190, gi_192, gi_194 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = f_15 * gi_32[k]
                 + f_16 * gi_39[k]
                 - f_17 * gi_41[k]
                 + f_15 * gi_50[k]
                 - f_17 * gi_52[k]
                 + f_18 * gi_54[k]
                 - f_15 * gi_172[k]
                 - f_16 * gi_179[k]
                 + f_17 * gi_181[k]
                 - f_15 * gi_190[k]
                 + f_17 * gi_192[k]
                 - f_18 * gi_194[k];
    }

#pragma omp simd aligned(gi_28, gi_31, gi_33, gi_38, gi_40, gi_42, gi_49, gi_51, gi_53, gi_55, \
                         gi_168, gi_171, gi_173, gi_178, gi_180, gi_182, gi_189, gi_191, \
                         gi_193, gi_195 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = -f_19 * gi_28[k]
                 - f_20 * gi_31[k]
                 + f_21 * gi_33[k]
                 - f_20 * gi_38[k]
                 + f_22 * gi_40[k]
                 - f_23 * gi_42[k]
                 - f_19 * gi_49[k]
                 + f_21 * gi_51[k]
                 - f_23 * gi_53[k]
                 + f_24 * gi_55[k]
                 + f_19 * gi_168[k]
                 + f_20 * gi_171[k]
                 - f_21 * gi_173[k]
                 + f_20 * gi_178[k]
                 - f_22 * gi_180[k]
                 + f_23 * gi_182[k]
                 + f_19 * gi_189[k]
                 - f_21 * gi_191[k]
                 + f_23 * gi_193[k]
                 - f_24 * gi_195[k];
    }

#pragma omp simd aligned(gi_30, gi_35, gi_37, gi_44, gi_46, gi_48, gi_170, gi_175, gi_177, \
                         gi_184, gi_186, gi_188 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = f_15 * gi_30[k]
                 + f_16 * gi_35[k]
                 - f_17 * gi_37[k]
                 + f_15 * gi_44[k]
                 - f_17 * gi_46[k]
                 + f_18 * gi_48[k]
                 - f_15 * gi_170[k]
                 - f_16 * gi_175[k]
                 + f_17 * gi_177[k]
                 - f_15 * gi_184[k]
                 + f_17 * gi_186[k]
                 - f_18 * gi_188[k];
    }

#pragma omp simd aligned(gi_28, gi_31, gi_33, gi_38, gi_42, gi_49, gi_51, gi_53, gi_168, \
                         gi_171, gi_173, gi_178, gi_182, gi_189, gi_191, \
                         gi_193 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_8[k] = f_25 * gi_28[k]
                 + f_25 * gi_31[k]
                 - f_11 * gi_33[k]
                 - f_25 * gi_38[k]
                 + f_11 * gi_42[k]
                 - f_25 * gi_49[k]
                 + f_11 * gi_51[k]
                 - f_11 * gi_53[k]
                 - f_25 * gi_168[k]
                 - f_25 * gi_171[k]
                 + f_11 * gi_173[k]
                 + f_25 * gi_178[k]
                 - f_11 * gi_182[k]
                 + f_25 * gi_189[k]
                 - f_11 * gi_191[k]
                 + f_11 * gi_193[k];
    }

#pragma omp simd aligned(gi_30, gi_35, gi_37, gi_44, gi_46, gi_170, gi_175, gi_177, gi_184, \
                         gi_186 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_9[k] = -f_10 * gi_30[k]
                 + f_8 * gi_35[k]
                 + f_11 * gi_37[k]
                 + f_7 * gi_44[k]
                 - f_9 * gi_46[k]
                 + f_10 * gi_170[k]
                 - f_8 * gi_175[k]
                 - f_11 * gi_177[k]
                 - f_7 * gi_184[k]
                 + f_9 * gi_186[k];
    }

#pragma omp simd aligned(gi_28, gi_31, gi_33, gi_38, gi_40, gi_49, gi_51, gi_168, gi_171, \
                         gi_173, gi_178, gi_180, gi_189, gi_191 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = -f_26 * gi_28[k]
                  + f_27 * gi_31[k]
                  + f_28 * gi_33[k]
                  + f_27 * gi_38[k]
                  - f_29 * gi_40[k]
                  - f_26 * gi_49[k]
                  + f_28 * gi_51[k]
                  + f_26 * gi_168[k]
                  - f_27 * gi_171[k]
                  - f_28 * gi_173[k]
                  - f_27 * gi_178[k]
                  + f_29 * gi_180[k]
                  + f_26 * gi_189[k]
                  - f_28 * gi_191[k];
    }

#pragma omp simd aligned(gi_28, gi_30, gi_31, gi_35, gi_38, gi_44, gi_49, gi_168, gi_170, \
                         gi_171, gi_175, gi_178, gi_184, gi_189 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_11[k] = f_4 * gi_30[k]
                  - f_3 * gi_35[k]
                  + f_2 * gi_44[k]
                  - f_4 * gi_170[k]
                  + f_3 * gi_175[k]
                  - f_2 * gi_184[k];

        g_12[k] = f_30 * gi_28[k]
                  - f_31 * gi_31[k]
                  + f_31 * gi_38[k]
                  - f_30 * gi_49[k]
                  - f_30 * gi_168[k]
                  + f_31 * gi_171[k]
                  - f_31 * gi_178[k]
                  + f_30 * gi_189[k];
    }

#pragma omp simd aligned(gi_113, gi_116, gi_118, gi_123, gi_127, gi_134, gi_309, gi_312, \
                         gi_314, gi_319, gi_323, gi_330 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = f_32 * gi_113[k]
                  - f_33 * gi_118[k]
                  + f_32 * gi_127[k]
                  - f_34 * gi_309[k]
                  + f_35 * gi_314[k]
                  - f_34 * gi_323[k];

        g_14[k] = f_36 * gi_116[k]
                  - f_37 * gi_123[k]
                  + f_38 * gi_134[k]
                  - f_39 * gi_312[k]
                  + f_40 * gi_319[k]
                  - f_41 * gi_330[k];
    }

#pragma omp simd aligned(gi_113, gi_120, gi_127, gi_129, gi_309, gi_316, gi_323, \
                         gi_325 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_15[k] = -f_42 * gi_113[k]
                  + f_43 * gi_120[k]
                  + f_42 * gi_127[k]
                  - f_43 * gi_129[k]
                  + f_44 * gi_309[k]
                  - f_45 * gi_316[k]
                  - f_44 * gi_323[k]
                  + f_45 * gi_325[k];
    }

#pragma omp simd aligned(gi_116, gi_123, gi_125, gi_134, gi_136, gi_312, gi_319, gi_321, \
                         gi_330, gi_332 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_16[k] = -f_46 * gi_116[k]
                  - f_47 * gi_123[k]
                  + f_48 * gi_125[k]
                  + f_49 * gi_134[k]
                  - f_50 * gi_136[k]
                  + f_49 * gi_312[k]
                  + f_51 * gi_319[k]
                  - f_50 * gi_321[k]
                  - f_52 * gi_330[k]
                  + f_53 * gi_332[k];
    }

#pragma omp simd aligned(gi_113, gi_118, gi_120, gi_127, gi_129, gi_131, gi_309, gi_314, \
                         gi_316, gi_323, gi_325, gi_327 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = f_52 * gi_113[k]
                  + f_51 * gi_118[k]
                  - f_54 * gi_120[k]
                  + f_52 * gi_127[k]
                  - f_54 * gi_129[k]
                  + f_54 * gi_131[k]
                  - f_55 * gi_309[k]
                  - f_56 * gi_314[k]
                  + f_57 * gi_316[k]
                  - f_55 * gi_323[k]
                  + f_57 * gi_325[k]
                  - f_57 * gi_327[k];
    }

#pragma omp simd aligned(gi_116, gi_123, gi_125, gi_134, gi_136, gi_138, gi_312, gi_319, \
                         gi_321, gi_330, gi_332, gi_334 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = f_58 * gi_116[k]
                  + f_59 * gi_123[k]
                  - f_60 * gi_125[k]
                  + f_58 * gi_134[k]
                  - f_60 * gi_136[k]
                  + f_61 * gi_138[k]
                  - f_62 * gi_312[k]
                  - f_63 * gi_319[k]
                  + f_64 * gi_321[k]
                  - f_62 * gi_330[k]
                  + f_64 * gi_332[k]
                  - f_65 * gi_334[k];
    }

#pragma omp simd aligned(gi_112, gi_115, gi_117, gi_122, gi_124, gi_126, gi_133, gi_135, \
                         gi_137, gi_139, gi_308, gi_311, gi_313, gi_318, gi_320, gi_322, \
                         gi_329, gi_331, gi_333, gi_335 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = -f_66 * gi_112[k]
                  - f_67 * gi_115[k]
                  + f_68 * gi_117[k]
                  - f_67 * gi_122[k]
                  + f_69 * gi_124[k]
                  - f_70 * gi_126[k]
                  - f_66 * gi_133[k]
                  + f_68 * gi_135[k]
                  - f_70 * gi_137[k]
                  + f_71 * gi_139[k]
                  + f_72 * gi_308[k]
                  + f_66 * gi_311[k]
                  - f_73 * gi_313[k]
                  + f_66 * gi_318[k]
                  - f_74 * gi_320[k]
                  + f_75 * gi_322[k]
                  + f_72 * gi_329[k]
                  - f_73 * gi_331[k]
                  + f_75 * gi_333[k]
                  - f_76 * gi_335[k];
    }

#pragma omp simd aligned(gi_114, gi_119, gi_121, gi_128, gi_130, gi_132, gi_310, gi_315, \
                         gi_317, gi_324, gi_326, gi_328 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = f_58 * gi_114[k]
                  + f_59 * gi_119[k]
                  - f_60 * gi_121[k]
                  + f_58 * gi_128[k]
                  - f_60 * gi_130[k]
                  + f_61 * gi_132[k]
                  - f_62 * gi_310[k]
                  - f_63 * gi_315[k]
                  + f_64 * gi_317[k]
                  - f_62 * gi_324[k]
                  + f_64 * gi_326[k]
                  - f_65 * gi_328[k];
    }

#pragma omp simd aligned(gi_112, gi_115, gi_117, gi_122, gi_126, gi_133, gi_135, gi_137, \
                         gi_308, gi_311, gi_313, gi_318, gi_322, gi_329, gi_331, \
                         gi_333 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = f_77 * gi_112[k]
                  + f_77 * gi_115[k]
                  - f_50 * gi_117[k]
                  - f_77 * gi_122[k]
                  + f_50 * gi_126[k]
                  - f_77 * gi_133[k]
                  + f_50 * gi_135[k]
                  - f_50 * gi_137[k]
                  - f_78 * gi_308[k]
                  - f_78 * gi_311[k]
                  + f_53 * gi_313[k]
                  + f_78 * gi_318[k]
                  - f_53 * gi_322[k]
                  + f_78 * gi_329[k]
                  - f_53 * gi_331[k]
                  + f_53 * gi_333[k];
    }

#pragma omp simd aligned(gi_114, gi_119, gi_121, gi_128, gi_130, gi_310, gi_315, gi_317, \
                         gi_324, gi_326 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = -f_49 * gi_114[k]
                  + f_47 * gi_119[k]
                  + f_50 * gi_121[k]
                  + f_46 * gi_128[k]
                  - f_48 * gi_130[k]
                  + f_52 * gi_310[k]
                  - f_51 * gi_315[k]
                  - f_53 * gi_317[k]
                  - f_49 * gi_324[k]
                  + f_50 * gi_326[k];
    }

#pragma omp simd aligned(gi_112, gi_115, gi_117, gi_122, gi_124, gi_133, gi_135, gi_308, \
                         gi_311, gi_313, gi_318, gi_320, gi_329, \
                         gi_331 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = -f_79 * gi_112[k]
                  + f_80 * gi_115[k]
                  + f_81 * gi_117[k]
                  + f_80 * gi_122[k]
                  - f_82 * gi_124[k]
                  - f_79 * gi_133[k]
                  + f_81 * gi_135[k]
                  + f_83 * gi_308[k]
                  - f_84 * gi_311[k]
                  - f_85 * gi_313[k]
                  - f_84 * gi_318[k]
                  + f_86 * gi_320[k]
                  + f_83 * gi_329[k]
                  - f_85 * gi_331[k];
    }

#pragma omp simd aligned(gi_112, gi_114, gi_115, gi_119, gi_122, gi_128, gi_133, gi_308, \
                         gi_310, gi_311, gi_315, gi_318, gi_324, \
                         gi_329 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = f_38 * gi_114[k]
                  - f_37 * gi_119[k]
                  + f_36 * gi_128[k]
                  - f_41 * gi_310[k]
                  + f_40 * gi_315[k]
                  - f_39 * gi_324[k];

        g_25[k] = f_87 * gi_112[k]
                  - f_88 * gi_115[k]
                  + f_88 * gi_122[k]
                  - f_87 * gi_133[k]
                  - f_89 * gi_308[k]
                  + f_90 * gi_311[k]
                  - f_90 * gi_318[k]
                  + f_89 * gi_329[k];
    }

#pragma omp simd aligned(gi_29, gi_34, gi_43, gi_169, gi_174, gi_183, gi_225, gi_230, \
                         gi_239 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = -f_91 * gi_29[k]
                  + f_92 * gi_34[k]
                  - f_91 * gi_43[k]
                  - f_91 * gi_169[k]
                  + f_92 * gi_174[k]
                  - f_91 * gi_183[k]
                  + f_93 * gi_225[k]
                  - f_94 * gi_230[k]
                  + f_93 * gi_239[k];
    }

#pragma omp simd aligned(gi_32, gi_39, gi_50, gi_172, gi_179, gi_190, gi_228, gi_235, \
                         gi_246 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = -f_95 * gi_32[k]
                  + f_96 * gi_39[k]
                  - f_97 * gi_50[k]
                  - f_95 * gi_172[k]
                  + f_96 * gi_179[k]
                  - f_97 * gi_190[k]
                  + f_98 * gi_228[k]
                  - f_99 * gi_235[k]
                  + f_100 * gi_246[k];
    }

#pragma omp simd aligned(gi_29, gi_36, gi_43, gi_45, gi_169, gi_176, gi_183, gi_185, gi_225, \
                         gi_232, gi_239, gi_241 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = f_101 * gi_29[k]
                  - f_23 * gi_36[k]
                  - f_101 * gi_43[k]
                  + f_23 * gi_45[k]
                  + f_101 * gi_169[k]
                  - f_23 * gi_176[k]
                  - f_101 * gi_183[k]
                  + f_23 * gi_185[k]
                  - f_102 * gi_225[k]
                  + f_103 * gi_232[k]
                  + f_102 * gi_239[k]
                  - f_103 * gi_241[k];
    }

#pragma omp simd aligned(gi_32, gi_39, gi_41, gi_50, gi_52, gi_172, gi_179, gi_181, gi_190, \
                         gi_192, gi_228, gi_235, gi_237, gi_246, \
                         gi_248 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = f_104 * gi_32[k]
                  + f_105 * gi_39[k]
                  - f_106 * gi_41[k]
                  - f_107 * gi_50[k]
                  + f_108 * gi_52[k]
                  + f_104 * gi_172[k]
                  + f_105 * gi_179[k]
                  - f_106 * gi_181[k]
                  - f_107 * gi_190[k]
                  + f_108 * gi_192[k]
                  - f_109 * gi_228[k]
                  - f_110 * gi_235[k]
                  + f_111 * gi_237[k]
                  + f_112 * gi_246[k]
                  - f_113 * gi_248[k];
    }

#pragma omp simd aligned(gi_29, gi_34, gi_36, gi_43, gi_45, gi_47, gi_169, gi_174, gi_176, \
                         gi_183, gi_185, gi_187, gi_225, gi_230, gi_232, gi_239, gi_241, \
                         gi_243 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = -f_114 * gi_29[k]
                  - f_115 * gi_34[k]
                  + f_116 * gi_36[k]
                  - f_114 * gi_43[k]
                  + f_116 * gi_45[k]
                  - f_116 * gi_47[k]
                  - f_114 * gi_169[k]
                  - f_115 * gi_174[k]
                  + f_116 * gi_176[k]
                  - f_114 * gi_183[k]
                  + f_116 * gi_185[k]
                  - f_116 * gi_187[k]
                  + f_105 * gi_225[k]
                  + f_117 * gi_230[k]
                  - f_118 * gi_232[k]
                  + f_105 * gi_239[k]
                  - f_118 * gi_241[k]
                  + f_118 * gi_243[k];
    }

#pragma omp simd aligned(gi_32, gi_39, gi_41, gi_50, gi_52, gi_54, gi_172, gi_179, gi_181, \
                         gi_190, gi_192, gi_194, gi_228, gi_235, gi_237, gi_246, gi_248, \
                         gi_250 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = -f_119 * gi_32[k]
                  - f_120 * gi_39[k]
                  + f_121 * gi_41[k]
                  - f_119 * gi_50[k]
                  + f_121 * gi_52[k]
                  - f_122 * gi_54[k]
                  - f_119 * gi_172[k]
                  - f_120 * gi_179[k]
                  + f_121 * gi_181[k]
                  - f_119 * gi_190[k]
                  + f_121 * gi_192[k]
                  - f_122 * gi_194[k]
                  + f_123 * gi_228[k]
                  + f_124 * gi_235[k]
                  - f_125 * gi_237[k]
                  + f_123 * gi_246[k]
                  - f_125 * gi_248[k]
                  + f_126 * gi_250[k];
    }

#pragma omp simd aligned(gi_28, gi_31, gi_33, gi_38, gi_40, gi_42, gi_49, gi_51, gi_53, gi_55, \
                         gi_168, gi_171, gi_173, gi_178, gi_180, gi_182, gi_189, gi_191, \
                         gi_193, gi_195, gi_224, gi_227, gi_229, gi_234, gi_236, gi_238, \
                         gi_245, gi_247, gi_249, gi_251 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = f_127 * gi_28[k]
                  + f_128 * gi_31[k]
                  - f_129 * gi_33[k]
                  + f_128 * gi_38[k]
                  - f_130 * gi_40[k]
                  + f_131 * gi_42[k]
                  + f_127 * gi_49[k]
                  - f_129 * gi_51[k]
                  + f_131 * gi_53[k]
                  - f_132 * gi_55[k]
                  + f_127 * gi_168[k]
                  + f_128 * gi_171[k]
                  - f_129 * gi_173[k]
                  + f_128 * gi_178[k]
                  - f_130 * gi_180[k]
                  + f_131 * gi_182[k]
                  + f_127 * gi_189[k]
                  - f_129 * gi_191[k]
                  + f_131 * gi_193[k]
                  - f_132 * gi_195[k]
                  - f_133 * gi_224[k]
                  - f_129 * gi_227[k]
                  + f_134 * gi_229[k]
                  - f_129 * gi_234[k]
                  + f_135 * gi_236[k]
                  - f_136 * gi_238[k]
                  - f_133 * gi_245[k]
                  + f_134 * gi_247[k]
                  - f_136 * gi_249[k]
                  + f_137 * gi_251[k];
    }

#pragma omp simd aligned(gi_30, gi_35, gi_37, gi_44, gi_46, gi_48, gi_170, gi_175, gi_177, \
                         gi_184, gi_186, gi_188, gi_226, gi_231, gi_233, gi_240, gi_242, \
                         gi_244 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = -f_119 * gi_30[k]
                  - f_120 * gi_35[k]
                  + f_121 * gi_37[k]
                  - f_119 * gi_44[k]
                  + f_121 * gi_46[k]
                  - f_122 * gi_48[k]
                  - f_119 * gi_170[k]
                  - f_120 * gi_175[k]
                  + f_121 * gi_177[k]
                  - f_119 * gi_184[k]
                  + f_121 * gi_186[k]
                  - f_122 * gi_188[k]
                  + f_123 * gi_226[k]
                  + f_124 * gi_231[k]
                  - f_125 * gi_233[k]
                  + f_123 * gi_240[k]
                  - f_125 * gi_242[k]
                  + f_126 * gi_244[k];
    }

#pragma omp simd aligned(gi_28, gi_31, gi_33, gi_38, gi_42, gi_49, gi_51, gi_53, gi_168, \
                         gi_171, gi_173, gi_178, gi_182, gi_189, gi_191, gi_193, gi_224, \
                         gi_227, gi_229, gi_234, gi_238, gi_245, gi_247, \
                         gi_249 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = -f_138 * gi_28[k]
                  - f_138 * gi_31[k]
                  + f_108 * gi_33[k]
                  + f_138 * gi_38[k]
                  - f_108 * gi_42[k]
                  + f_138 * gi_49[k]
                  - f_108 * gi_51[k]
                  + f_108 * gi_53[k]
                  - f_138 * gi_168[k]
                  - f_138 * gi_171[k]
                  + f_108 * gi_173[k]
                  + f_138 * gi_178[k]
                  - f_108 * gi_182[k]
                  + f_138 * gi_189[k]
                  - f_108 * gi_191[k]
                  + f_108 * gi_193[k]
                  + f_107 * gi_224[k]
                  + f_107 * gi_227[k]
                  - f_113 * gi_229[k]
                  - f_107 * gi_234[k]
                  + f_113 * gi_238[k]
                  - f_107 * gi_245[k]
                  + f_113 * gi_247[k]
                  - f_113 * gi_249[k];
    }

#pragma omp simd aligned(gi_30, gi_35, gi_37, gi_44, gi_46, gi_170, gi_175, gi_177, gi_184, \
                         gi_186, gi_226, gi_231, gi_233, gi_240, \
                         gi_242 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = f_107 * gi_30[k]
                  - f_105 * gi_35[k]
                  - f_108 * gi_37[k]
                  - f_104 * gi_44[k]
                  + f_106 * gi_46[k]
                  + f_107 * gi_170[k]
                  - f_105 * gi_175[k]
                  - f_108 * gi_177[k]
                  - f_104 * gi_184[k]
                  + f_106 * gi_186[k]
                  - f_112 * gi_226[k]
                  + f_110 * gi_231[k]
                  + f_113 * gi_233[k]
                  + f_109 * gi_240[k]
                  - f_111 * gi_242[k];
    }

#pragma omp simd aligned(gi_28, gi_31, gi_33, gi_38, gi_40, gi_49, gi_51, gi_168, gi_171, \
                         gi_173, gi_178, gi_180, gi_189, gi_191, gi_224, gi_227, gi_229, \
                         gi_234, gi_236, gi_245, gi_247 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_36[k] = f_139 * gi_28[k]
                  - f_20 * gi_31[k]
                  - f_140 * gi_33[k]
                  - f_20 * gi_38[k]
                  + f_22 * gi_40[k]
                  + f_139 * gi_49[k]
                  - f_140 * gi_51[k]
                  + f_139 * gi_168[k]
                  - f_20 * gi_171[k]
                  - f_140 * gi_173[k]
                  - f_20 * gi_178[k]
                  + f_22 * gi_180[k]
                  + f_139 * gi_189[k]
                  - f_140 * gi_191[k]
                  - f_141 * gi_224[k]
                  + f_21 * gi_227[k]
                  + f_22 * gi_229[k]
                  + f_21 * gi_234[k]
                  - f_142 * gi_236[k]
                  - f_141 * gi_245[k]
                  + f_22 * gi_247[k];
    }

#pragma omp simd aligned(gi_30, gi_35, gi_44, gi_170, gi_175, gi_184, gi_226, gi_231, \
                         gi_240 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_37[k] = -f_97 * gi_30[k]
                  + f_96 * gi_35[k]
                  - f_95 * gi_44[k]
                  - f_97 * gi_170[k]
                  + f_96 * gi_175[k]
                  - f_95 * gi_184[k]
                  + f_100 * gi_226[k]
                  - f_99 * gi_231[k]
                  + f_98 * gi_240[k];
    }

#pragma omp simd aligned(gi_28, gi_31, gi_38, gi_49, gi_168, gi_171, gi_178, gi_189, gi_224, \
                         gi_227, gi_234, gi_245 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = -f_143 * gi_28[k]
                  + f_144 * gi_31[k]
                  - f_144 * gi_38[k]
                  + f_143 * gi_49[k]
                  - f_143 * gi_168[k]
                  + f_144 * gi_171[k]
                  - f_144 * gi_178[k]
                  + f_143 * gi_189[k]
                  + f_91 * gi_224[k]
                  - f_145 * gi_227[k]
                  + f_145 * gi_234[k]
                  - f_91 * gi_245[k];
    }

#pragma omp simd aligned(gi_113, gi_118, gi_127, gi_309, gi_314, gi_323, gi_365, gi_370, \
                         gi_379 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_39[k] = -f_146 * gi_113[k]
                  + f_147 * gi_118[k]
                  - f_146 * gi_127[k]
                  - f_146 * gi_309[k]
                  + f_147 * gi_314[k]
                  - f_146 * gi_323[k]
                  + f_148 * gi_365[k]
                  - f_149 * gi_370[k]
                  + f_148 * gi_379[k];
    }

#pragma omp simd aligned(gi_116, gi_123, gi_134, gi_312, gi_319, gi_330, gi_368, gi_375, \
                         gi_386 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_40[k] = -f_150 * gi_116[k]
                  + f_151 * gi_123[k]
                  - f_152 * gi_134[k]
                  - f_150 * gi_312[k]
                  + f_151 * gi_319[k]
                  - f_152 * gi_330[k]
                  + f_153 * gi_368[k]
                  - f_154 * gi_375[k]
                  + f_155 * gi_386[k];
    }

#pragma omp simd aligned(gi_113, gi_120, gi_127, gi_129, gi_309, gi_316, gi_323, gi_325, \
                         gi_365, gi_372, gi_379, gi_381 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_41[k] = f_156 * gi_113[k]
                  - f_70 * gi_120[k]
                  - f_156 * gi_127[k]
                  + f_70 * gi_129[k]
                  + f_156 * gi_309[k]
                  - f_70 * gi_316[k]
                  - f_156 * gi_323[k]
                  + f_70 * gi_325[k]
                  - f_71 * gi_365[k]
                  + f_157 * gi_372[k]
                  + f_71 * gi_379[k]
                  - f_157 * gi_381[k];
    }

#pragma omp simd aligned(gi_116, gi_123, gi_125, gi_134, gi_136, gi_312, gi_319, gi_321, \
                         gi_330, gi_332, gi_368, gi_375, gi_377, gi_386, \
                         gi_388 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = f_158 * gi_116[k]
                  + f_159 * gi_123[k]
                  - f_160 * gi_125[k]
                  - f_161 * gi_134[k]
                  + f_162 * gi_136[k]
                  + f_158 * gi_312[k]
                  + f_159 * gi_319[k]
                  - f_160 * gi_321[k]
                  - f_161 * gi_330[k]
                  + f_162 * gi_332[k]
                  - f_163 * gi_368[k]
                  - f_162 * gi_375[k]
                  + f_164 * gi_377[k]
                  + f_165 * gi_386[k]
                  - f_166 * gi_388[k];
    }

#pragma omp simd aligned(gi_113, gi_118, gi_120, gi_127, gi_129, gi_131, gi_309, gi_314, \
                         gi_316, gi_323, gi_325, gi_327, gi_365, gi_370, gi_372, gi_379, \
                         gi_381, gi_383 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_43[k] = -f_167 * gi_113[k]
                  - f_168 * gi_118[k]
                  + f_169 * gi_120[k]
                  - f_167 * gi_127[k]
                  + f_169 * gi_129[k]
                  - f_169 * gi_131[k]
                  - f_167 * gi_309[k]
                  - f_168 * gi_314[k]
                  + f_169 * gi_316[k]
                  - f_167 * gi_323[k]
                  + f_169 * gi_325[k]
                  - f_169 * gi_327[k]
                  + f_170 * gi_365[k]
                  + f_171 * gi_370[k]
                  - f_172 * gi_372[k]
                  + f_170 * gi_379[k]
                  - f_172 * gi_381[k]
                  + f_172 * gi_383[k];
    }

#pragma omp simd aligned(gi_116, gi_123, gi_125, gi_134, gi_136, gi_138, gi_312, gi_319, \
                         gi_321, gi_330, gi_332, gi_334, gi_368, gi_375, gi_377, gi_386, \
                         gi_388, gi_390 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_44[k] = -f_173 * gi_116[k]
                  - f_174 * gi_123[k]
                  + f_175 * gi_125[k]
                  - f_173 * gi_134[k]
                  + f_175 * gi_136[k]
                  - f_176 * gi_138[k]
                  - f_173 * gi_312[k]
                  - f_174 * gi_319[k]
                  + f_175 * gi_321[k]
                  - f_173 * gi_330[k]
                  + f_175 * gi_332[k]
                  - f_176 * gi_334[k]
                  + f_177 * gi_368[k]
                  + f_178 * gi_375[k]
                  - f_179 * gi_377[k]
                  + f_177 * gi_386[k]
                  - f_179 * gi_388[k]
                  + f_180 * gi_390[k];
    }

#pragma omp simd aligned(gi_112, gi_115, gi_117, gi_122, gi_124, gi_126, gi_133, gi_135, \
                         gi_137, gi_139, gi_308, gi_311, gi_313, gi_318, gi_320, gi_322, \
                         gi_329, gi_331, gi_333, gi_335, gi_364, gi_367, gi_369, gi_374, \
                         gi_376, gi_378, gi_385, gi_387, gi_389, \
                         gi_391 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_45[k] = f_181 * gi_112[k]
                  + f_182 * gi_115[k]
                  - f_183 * gi_117[k]
                  + f_182 * gi_122[k]
                  - f_184 * gi_124[k]
                  + f_185 * gi_126[k]
                  + f_181 * gi_133[k]
                  - f_183 * gi_135[k]
                  + f_185 * gi_137[k]
                  - f_186 * gi_139[k]
                  + f_181 * gi_308[k]
                  + f_182 * gi_311[k]
                  - f_183 * gi_313[k]
                  + f_182 * gi_318[k]
                  - f_184 * gi_320[k]
                  + f_185 * gi_322[k]
                  + f_181 * gi_329[k]
                  - f_183 * gi_331[k]
                  + f_185 * gi_333[k]
                  - f_186 * gi_335[k]
                  - f_187 * gi_364[k]
                  - f_188 * gi_367[k]
                  + f_185 * gi_369[k]
                  - f_188 * gi_374[k]
                  + f_189 * gi_376[k]
                  - f_190 * gi_378[k]
                  - f_187 * gi_385[k]
                  + f_185 * gi_387[k]
                  - f_190 * gi_389[k]
                  + f_191 * gi_391[k];
    }

#pragma omp simd aligned(gi_114, gi_119, gi_121, gi_128, gi_130, gi_132, gi_310, gi_315, \
                         gi_317, gi_324, gi_326, gi_328, gi_366, gi_371, gi_373, gi_380, \
                         gi_382, gi_384 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_46[k] = -f_173 * gi_114[k]
                  - f_174 * gi_119[k]
                  + f_175 * gi_121[k]
                  - f_173 * gi_128[k]
                  + f_175 * gi_130[k]
                  - f_176 * gi_132[k]
                  - f_173 * gi_310[k]
                  - f_174 * gi_315[k]
                  + f_175 * gi_317[k]
                  - f_173 * gi_324[k]
                  + f_175 * gi_326[k]
                  - f_176 * gi_328[k]
                  + f_177 * gi_366[k]
                  + f_178 * gi_371[k]
                  - f_179 * gi_373[k]
                  + f_177 * gi_380[k]
                  - f_179 * gi_382[k]
                  + f_180 * gi_384[k];
    }

#pragma omp simd aligned(gi_112, gi_115, gi_117, gi_122, gi_126, gi_133, gi_135, gi_137, \
                         gi_308, gi_311, gi_313, gi_318, gi_322, gi_329, gi_331, gi_333, \
                         gi_364, gi_367, gi_369, gi_374, gi_378, gi_385, gi_387, \
                         gi_389 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_47[k] = -f_192 * gi_112[k]
                  - f_192 * gi_115[k]
                  + f_162 * gi_117[k]
                  + f_192 * gi_122[k]
                  - f_162 * gi_126[k]
                  + f_192 * gi_133[k]
                  - f_162 * gi_135[k]
                  + f_162 * gi_137[k]
                  - f_192 * gi_308[k]
                  - f_192 * gi_311[k]
                  + f_162 * gi_313[k]
                  + f_192 * gi_318[k]
                  - f_162 * gi_322[k]
                  + f_192 * gi_329[k]
                  - f_162 * gi_331[k]
                  + f_162 * gi_333[k]
                  + f_193 * gi_364[k]
                  + f_193 * gi_367[k]
                  - f_166 * gi_369[k]
                  - f_193 * gi_374[k]
                  + f_166 * gi_378[k]
                  - f_193 * gi_385[k]
                  + f_166 * gi_387[k]
                  - f_166 * gi_389[k];
    }

#pragma omp simd aligned(gi_114, gi_119, gi_121, gi_128, gi_130, gi_310, gi_315, gi_317, \
                         gi_324, gi_326, gi_366, gi_371, gi_373, gi_380, \
                         gi_382 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_48[k] = f_161 * gi_114[k]
                  - f_159 * gi_119[k]
                  - f_162 * gi_121[k]
                  - f_158 * gi_128[k]
                  + f_160 * gi_130[k]
                  + f_161 * gi_310[k]
                  - f_159 * gi_315[k]
                  - f_162 * gi_317[k]
                  - f_158 * gi_324[k]
                  + f_160 * gi_326[k]
                  - f_165 * gi_366[k]
                  + f_162 * gi_371[k]
                  + f_166 * gi_373[k]
                  + f_163 * gi_380[k]
                  - f_164 * gi_382[k];
    }

#pragma omp simd aligned(gi_112, gi_115, gi_117, gi_122, gi_124, gi_133, gi_135, gi_308, \
                         gi_311, gi_313, gi_318, gi_320, gi_329, gi_331, gi_364, gi_367, \
                         gi_369, gi_374, gi_376, gi_385, gi_387 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_49[k] = f_194 * gi_112[k]
                  - f_67 * gi_115[k]
                  - f_73 * gi_117[k]
                  - f_67 * gi_122[k]
                  + f_69 * gi_124[k]
                  + f_194 * gi_133[k]
                  - f_73 * gi_135[k]
                  + f_194 * gi_308[k]
                  - f_67 * gi_311[k]
                  - f_73 * gi_313[k]
                  - f_67 * gi_318[k]
                  + f_69 * gi_320[k]
                  + f_194 * gi_329[k]
                  - f_73 * gi_331[k]
                  - f_195 * gi_364[k]
                  + f_196 * gi_367[k]
                  + f_75 * gi_369[k]
                  + f_196 * gi_374[k]
                  - f_197 * gi_376[k]
                  - f_195 * gi_385[k]
                  + f_75 * gi_387[k];
    }

#pragma omp simd aligned(gi_114, gi_119, gi_128, gi_310, gi_315, gi_324, gi_366, gi_371, \
                         gi_380 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_50[k] = -f_152 * gi_114[k]
                  + f_151 * gi_119[k]
                  - f_150 * gi_128[k]
                  - f_152 * gi_310[k]
                  + f_151 * gi_315[k]
                  - f_150 * gi_324[k]
                  + f_155 * gi_366[k]
                  - f_154 * gi_371[k]
                  + f_153 * gi_380[k];
    }

#pragma omp simd aligned(gi_112, gi_115, gi_122, gi_133, gi_308, gi_311, gi_318, gi_329, \
                         gi_364, gi_367, gi_374, gi_385 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_51[k] = -f_198 * gi_112[k]
                  + f_199 * gi_115[k]
                  - f_199 * gi_122[k]
                  + f_198 * gi_133[k]
                  - f_198 * gi_308[k]
                  + f_199 * gi_311[k]
                  - f_199 * gi_318[k]
                  + f_198 * gi_329[k]
                  + f_200 * gi_364[k]
                  - f_147 * gi_367[k]
                  + f_147 * gi_374[k]
                  - f_200 * gi_385[k];
    }

#pragma omp simd aligned(gi_1, gi_6, gi_15, gi_85, gi_90, gi_99, gi_141, gi_146, gi_155, \
                         gi_281, gi_286, gi_295, gi_337, gi_342, gi_351, gi_393, gi_398, \
                         gi_407 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_52[k] = f_201 * gi_1[k]
                  - f_202 * gi_6[k]
                  + f_201 * gi_15[k]
                  + f_203 * gi_85[k]
                  - f_204 * gi_90[k]
                  + f_203 * gi_99[k]
                  - f_205 * gi_141[k]
                  + f_206 * gi_146[k]
                  - f_205 * gi_155[k]
                  + f_201 * gi_281[k]
                  - f_202 * gi_286[k]
                  + f_201 * gi_295[k]
                  - f_205 * gi_337[k]
                  + f_206 * gi_342[k]
                  - f_205 * gi_351[k]
                  + f_207 * gi_393[k]
                  - f_208 * gi_398[k]
                  + f_207 * gi_407[k];
    }

#pragma omp simd aligned(gi_4, gi_11, gi_22, gi_88, gi_95, gi_106, gi_144, gi_151, gi_162, \
                         gi_284, gi_291, gi_302, gi_340, gi_347, gi_358, gi_396, gi_403, \
                         gi_414 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_53[k] = f_209 * gi_4[k]
                  - f_210 * gi_11[k]
                  + f_211 * gi_22[k]
                  + f_210 * gi_88[k]
                  - f_212 * gi_95[k]
                  + f_213 * gi_106[k]
                  - f_214 * gi_144[k]
                  + f_215 * gi_151[k]
                  - f_216 * gi_162[k]
                  + f_209 * gi_284[k]
                  - f_210 * gi_291[k]
                  + f_211 * gi_302[k]
                  - f_214 * gi_340[k]
                  + f_215 * gi_347[k]
                  - f_216 * gi_358[k]
                  + f_217 * gi_396[k]
                  - f_218 * gi_403[k]
                  + f_219 * gi_414[k];
    }

#pragma omp simd aligned(gi_1, gi_8, gi_15, gi_17, gi_85, gi_92, gi_99, gi_101, gi_141, \
                         gi_148, gi_155, gi_157, gi_281, gi_288, gi_295, gi_297, gi_337, \
                         gi_344, gi_351, gi_353, gi_393, gi_400, gi_407, \
                         gi_409 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_54[k] = -f_220 * gi_1[k]
                  + f_221 * gi_8[k]
                  + f_220 * gi_15[k]
                  - f_221 * gi_17[k]
                  - f_222 * gi_85[k]
                  + f_223 * gi_92[k]
                  + f_222 * gi_99[k]
                  - f_223 * gi_101[k]
                  + f_224 * gi_141[k]
                  - f_225 * gi_148[k]
                  - f_224 * gi_155[k]
                  + f_225 * gi_157[k]
                  - f_220 * gi_281[k]
                  + f_221 * gi_288[k]
                  + f_220 * gi_295[k]
                  - f_221 * gi_297[k]
                  + f_224 * gi_337[k]
                  - f_225 * gi_344[k]
                  - f_224 * gi_351[k]
                  + f_225 * gi_353[k]
                  - f_226 * gi_393[k]
                  + f_227 * gi_400[k]
                  + f_226 * gi_407[k]
                  - f_227 * gi_409[k];
    }

#pragma omp simd aligned(gi_4, gi_11, gi_13, gi_22, gi_24, gi_88, gi_95, gi_97, gi_106, \
                         gi_108, gi_144, gi_151, gi_153, gi_162, gi_164, gi_284, gi_291, \
                         gi_293, gi_302, gi_304, gi_340, gi_347, gi_349, gi_358, gi_360, \
                         gi_396, gi_403, gi_405, gi_414, gi_416 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_55[k] = -f_228 * gi_4[k]
                  - f_229 * gi_11[k]
                  + f_230 * gi_13[k]
                  + f_231 * gi_22[k]
                  - f_232 * gi_24[k]
                  - f_233 * gi_88[k]
                  - f_234 * gi_95[k]
                  + f_235 * gi_97[k]
                  + f_229 * gi_106[k]
                  - f_236 * gi_108[k]
                  + f_237 * gi_144[k]
                  + f_235 * gi_151[k]
                  - f_238 * gi_153[k]
                  - f_230 * gi_162[k]
                  + f_239 * gi_164[k]
                  - f_228 * gi_284[k]
                  - f_229 * gi_291[k]
                  + f_230 * gi_293[k]
                  + f_231 * gi_302[k]
                  - f_232 * gi_304[k]
                  + f_237 * gi_340[k]
                  + f_235 * gi_347[k]
                  - f_238 * gi_349[k]
                  - f_230 * gi_358[k]
                  + f_239 * gi_360[k]
                  - f_230 * gi_396[k]
                  - f_236 * gi_403[k]
                  + f_239 * gi_405[k]
                  + f_232 * gi_414[k]
                  - f_240 * gi_416[k];
    }

#pragma omp simd aligned(gi_1, gi_6, gi_8, gi_15, gi_17, gi_19, gi_85, gi_90, gi_92, gi_99, \
                         gi_101, gi_103, gi_141, gi_146, gi_148, gi_155, gi_157, gi_159, \
                         gi_281, gi_286, gi_288, gi_295, gi_297, gi_299, gi_337, gi_342, \
                         gi_344, gi_351, gi_353, gi_355, gi_393, gi_398, gi_400, gi_407, \
                         gi_409, gi_411 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_56[k] = f_241 * gi_1[k]
                  + f_242 * gi_6[k]
                  - f_236 * gi_8[k]
                  + f_241 * gi_15[k]
                  - f_236 * gi_17[k]
                  + f_236 * gi_19[k]
                  + f_242 * gi_85[k]
                  + f_243 * gi_90[k]
                  - f_176 * gi_92[k]
                  + f_242 * gi_99[k]
                  - f_176 * gi_101[k]
                  + f_176 * gi_103[k]
                  - f_232 * gi_141[k]
                  - f_236 * gi_146[k]
                  + f_244 * gi_148[k]
                  - f_232 * gi_155[k]
                  + f_244 * gi_157[k]
                  - f_244 * gi_159[k]
                  + f_241 * gi_281[k]
                  + f_242 * gi_286[k]
                  - f_236 * gi_288[k]
                  + f_241 * gi_295[k]
                  - f_236 * gi_297[k]
                  + f_236 * gi_299[k]
                  - f_232 * gi_337[k]
                  - f_236 * gi_342[k]
                  + f_244 * gi_344[k]
                  - f_232 * gi_351[k]
                  + f_244 * gi_353[k]
                  - f_244 * gi_355[k]
                  + f_245 * gi_393[k]
                  + f_246 * gi_398[k]
                  - f_180 * gi_400[k]
                  + f_245 * gi_407[k]
                  - f_180 * gi_409[k]
                  + f_180 * gi_411[k];
    }

#pragma omp simd aligned(gi_4, gi_11, gi_13, gi_22, gi_24, gi_26, gi_88, gi_95, gi_97, gi_106, \
                         gi_108, gi_110, gi_144, gi_151, gi_153, gi_162, gi_164, gi_166, \
                         gi_284, gi_291, gi_293, gi_302, gi_304, gi_306, gi_340, gi_347, \
                         gi_349, gi_358, gi_360, gi_362, gi_396, gi_403, gi_405, gi_414, \
                         gi_416, gi_418 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_57[k] = f_192 * gi_4[k]
                  + f_167 * gi_11[k]
                  - f_168 * gi_13[k]
                  + f_192 * gi_22[k]
                  - f_168 * gi_24[k]
                  + f_247 * gi_26[k]
                  + f_167 * gi_88[k]
                  + f_168 * gi_95[k]
                  - f_165 * gi_97[k]
                  + f_167 * gi_106[k]
                  - f_165 * gi_108[k]
                  + f_248 * gi_110[k]
                  - f_165 * gi_144[k]
                  - f_162 * gi_151[k]
                  + f_169 * gi_153[k]
                  - f_165 * gi_162[k]
                  + f_169 * gi_164[k]
                  - f_249 * gi_166[k]
                  + f_192 * gi_284[k]
                  + f_167 * gi_291[k]
                  - f_168 * gi_293[k]
                  + f_192 * gi_302[k]
                  - f_168 * gi_304[k]
                  + f_247 * gi_306[k]
                  - f_165 * gi_340[k]
                  - f_162 * gi_347[k]
                  + f_169 * gi_349[k]
                  - f_165 * gi_358[k]
                  + f_169 * gi_360[k]
                  - f_249 * gi_362[k]
                  + f_170 * gi_396[k]
                  + f_171 * gi_403[k]
                  - f_250 * gi_405[k]
                  + f_170 * gi_414[k]
                  - f_250 * gi_416[k]
                  + f_251 * gi_418[k];
    }

#pragma omp simd aligned(gi_0, gi_3, gi_5, gi_10, gi_12, gi_14, gi_21, gi_23, gi_25, gi_27, \
                         gi_84, gi_87, gi_89, gi_94, gi_96, gi_98, gi_105, gi_107, gi_109, \
                         gi_111, gi_140, gi_143, gi_145, gi_150, gi_152, gi_154, gi_161, \
                         gi_163, gi_165, gi_167, gi_280, gi_283, gi_285, gi_290, gi_292, \
                         gi_294, gi_301, gi_303, gi_305, gi_307, gi_336, gi_339, gi_341, \
                         gi_346, gi_348, gi_350, gi_357, gi_359, gi_361, gi_363, gi_392, \
                         gi_395, gi_397, gi_402, gi_404, gi_406, gi_413, gi_415, gi_417, \
                         gi_419 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_58[k] = -0.1171875 * gi_0[k]
                  - 0.3515625 * gi_3[k]
                  + 2.109375 * gi_5[k]
                  - 0.3515625 * gi_10[k]
                  + 4.21875 * gi_12[k]
                  - 2.8125 * gi_14[k]
                  - 0.1171875 * gi_21[k]
                  + 2.109375 * gi_23[k]
                  - 2.8125 * gi_25[k]
                  + 0.375 * gi_27[k]
                  - 0.234375 * gi_84[k]
                  - 0.703125 * gi_87[k]
                  + 4.21875 * gi_89[k]
                  - 0.703125 * gi_94[k]
                  + 8.4375 * gi_96[k]
                  - 5.625 * gi_98[k]
                  - 0.234375 * gi_105[k]
                  + 4.21875 * gi_107[k]
                  - 5.625 * gi_109[k]
                  + 0.75 * gi_111[k]
                  + 0.9375 * gi_140[k]
                  + 2.8125 * gi_143[k]
                  - 16.875 * gi_145[k]
                  + 2.8125 * gi_150[k]
                  - 33.75 * gi_152[k]
                  + 22.5 * gi_154[k]
                  + 0.9375 * gi_161[k]
                  - 16.875 * gi_163[k]
                  + 22.5 * gi_165[k]
                  - 3.0 * gi_167[k]
                  - 0.1171875 * gi_280[k]
                  - 0.3515625 * gi_283[k]
                  + 2.109375 * gi_285[k]
                  - 0.3515625 * gi_290[k]
                  + 4.21875 * gi_292[k]
                  - 2.8125 * gi_294[k]
                  - 0.1171875 * gi_301[k]
                  + 2.109375 * gi_303[k]
                  - 2.8125 * gi_305[k]
                  + 0.375 * gi_307[k]
                  + 0.9375 * gi_336[k]
                  + 2.8125 * gi_339[k]
                  - 16.875 * gi_341[k]
                  + 2.8125 * gi_346[k]
                  - 33.75 * gi_348[k]
                  + 22.5 * gi_350[k]
                  + 0.9375 * gi_357[k]
                  - 16.875 * gi_359[k]
                  + 22.5 * gi_361[k]
                  - 3.0 * gi_363[k]
                  - 0.3125 * gi_392[k]
                  - 0.9375 * gi_395[k]
                  + 5.625 * gi_397[k]
                  - 0.9375 * gi_402[k]
                  + 11.25 * gi_404[k]
                  - 7.5 * gi_406[k]
                  - 0.3125 * gi_413[k]
                  + 5.625 * gi_415[k]
                  - 7.5 * gi_417[k]
                  + gi_419[k];
    }

#pragma omp simd aligned(gi_2, gi_7, gi_9, gi_16, gi_18, gi_20, gi_86, gi_91, gi_93, gi_100, \
                         gi_102, gi_104, gi_142, gi_147, gi_149, gi_156, gi_158, gi_160, \
                         gi_282, gi_287, gi_289, gi_296, gi_298, gi_300, gi_338, gi_343, \
                         gi_345, gi_352, gi_354, gi_356, gi_394, gi_399, gi_401, gi_408, \
                         gi_410, gi_412 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_59[k] = f_192 * gi_2[k]
                  + f_167 * gi_7[k]
                  - f_168 * gi_9[k]
                  + f_192 * gi_16[k]
                  - f_168 * gi_18[k]
                  + f_247 * gi_20[k]
                  + f_167 * gi_86[k]
                  + f_168 * gi_91[k]
                  - f_165 * gi_93[k]
                  + f_167 * gi_100[k]
                  - f_165 * gi_102[k]
                  + f_248 * gi_104[k]
                  - f_165 * gi_142[k]
                  - f_162 * gi_147[k]
                  + f_169 * gi_149[k]
                  - f_165 * gi_156[k]
                  + f_169 * gi_158[k]
                  - f_249 * gi_160[k]
                  + f_192 * gi_282[k]
                  + f_167 * gi_287[k]
                  - f_168 * gi_289[k]
                  + f_192 * gi_296[k]
                  - f_168 * gi_298[k]
                  + f_247 * gi_300[k]
                  - f_165 * gi_338[k]
                  - f_162 * gi_343[k]
                  + f_169 * gi_345[k]
                  - f_165 * gi_352[k]
                  + f_169 * gi_354[k]
                  - f_249 * gi_356[k]
                  + f_170 * gi_394[k]
                  + f_171 * gi_399[k]
                  - f_250 * gi_401[k]
                  + f_170 * gi_408[k]
                  - f_250 * gi_410[k]
                  + f_251 * gi_412[k];
    }

#pragma omp simd aligned(gi_0, gi_3, gi_5, gi_10, gi_14, gi_21, gi_23, gi_25, gi_84, gi_87, \
                         gi_89, gi_94, gi_98, gi_105, gi_107, gi_109, gi_140, gi_143, gi_145, \
                         gi_150, gi_154, gi_161, gi_163, gi_165, gi_280, gi_283, gi_285, \
                         gi_290, gi_294, gi_301, gi_303, gi_305, gi_336, gi_339, gi_341, \
                         gi_346, gi_350, gi_357, gi_359, gi_361, gi_392, gi_395, gi_397, \
                         gi_402, gi_406, gi_413, gi_415, gi_417 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_60[k] = f_252 * gi_0[k]
                  + f_252 * gi_3[k]
                  - f_232 * gi_5[k]
                  - f_252 * gi_10[k]
                  + f_232 * gi_14[k]
                  - f_252 * gi_21[k]
                  + f_232 * gi_23[k]
                  - f_232 * gi_25[k]
                  + f_241 * gi_84[k]
                  + f_241 * gi_87[k]
                  - f_236 * gi_89[k]
                  - f_241 * gi_94[k]
                  + f_236 * gi_98[k]
                  - f_241 * gi_105[k]
                  + f_236 * gi_107[k]
                  - f_236 * gi_109[k]
                  - f_243 * gi_140[k]
                  - f_243 * gi_143[k]
                  + f_239 * gi_145[k]
                  + f_243 * gi_150[k]
                  - f_239 * gi_154[k]
                  + f_243 * gi_161[k]
                  - f_239 * gi_163[k]
                  + f_239 * gi_165[k]
                  + f_252 * gi_280[k]
                  + f_252 * gi_283[k]
                  - f_232 * gi_285[k]
                  - f_252 * gi_290[k]
                  + f_232 * gi_294[k]
                  - f_252 * gi_301[k]
                  + f_232 * gi_303[k]
                  - f_232 * gi_305[k]
                  - f_243 * gi_336[k]
                  - f_243 * gi_339[k]
                  + f_239 * gi_341[k]
                  + f_243 * gi_346[k]
                  - f_239 * gi_350[k]
                  + f_243 * gi_357[k]
                  - f_239 * gi_359[k]
                  + f_239 * gi_361[k]
                  + f_253 * gi_392[k]
                  + f_253 * gi_395[k]
                  - f_240 * gi_397[k]
                  - f_253 * gi_402[k]
                  + f_240 * gi_406[k]
                  - f_253 * gi_413[k]
                  + f_240 * gi_415[k]
                  - f_240 * gi_417[k];
    }

#pragma omp simd aligned(gi_2, gi_7, gi_9, gi_16, gi_18, gi_86, gi_91, gi_93, gi_100, gi_102, \
                         gi_142, gi_147, gi_149, gi_156, gi_158, gi_282, gi_287, gi_289, \
                         gi_296, gi_298, gi_338, gi_343, gi_345, gi_352, gi_354, gi_394, \
                         gi_399, gi_401, gi_408, gi_410 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_61[k] = -f_231 * gi_2[k]
                  + f_229 * gi_7[k]
                  + f_232 * gi_9[k]
                  + f_228 * gi_16[k]
                  - f_230 * gi_18[k]
                  - f_229 * gi_86[k]
                  + f_234 * gi_91[k]
                  + f_236 * gi_93[k]
                  + f_233 * gi_100[k]
                  - f_235 * gi_102[k]
                  + f_230 * gi_142[k]
                  - f_235 * gi_147[k]
                  - f_239 * gi_149[k]
                  - f_237 * gi_156[k]
                  + f_238 * gi_158[k]
                  - f_231 * gi_282[k]
                  + f_229 * gi_287[k]
                  + f_232 * gi_289[k]
                  + f_228 * gi_296[k]
                  - f_230 * gi_298[k]
                  + f_230 * gi_338[k]
                  - f_235 * gi_343[k]
                  - f_239 * gi_345[k]
                  - f_237 * gi_352[k]
                  + f_238 * gi_354[k]
                  - f_232 * gi_394[k]
                  + f_236 * gi_399[k]
                  + f_240 * gi_401[k]
                  + f_230 * gi_408[k]
                  - f_239 * gi_410[k];
    }

#pragma omp simd aligned(gi_0, gi_3, gi_5, gi_10, gi_12, gi_21, gi_23, gi_84, gi_87, gi_89, \
                         gi_94, gi_96, gi_105, gi_107, gi_140, gi_143, gi_145, gi_150, gi_152, \
                         gi_161, gi_163, gi_280, gi_283, gi_285, gi_290, gi_292, gi_301, \
                         gi_303, gi_336, gi_339, gi_341, gi_346, gi_348, gi_357, gi_359, \
                         gi_392, gi_395, gi_397, gi_402, gi_404, gi_413, \
                         gi_415 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_62[k] = -f_254 * gi_0[k]
                  + f_255 * gi_3[k]
                  + f_256 * gi_5[k]
                  + f_255 * gi_10[k]
                  - f_257 * gi_12[k]
                  - f_254 * gi_21[k]
                  + f_256 * gi_23[k]
                  - f_258 * gi_84[k]
                  + f_256 * gi_87[k]
                  + f_259 * gi_89[k]
                  + f_256 * gi_94[k]
                  - f_260 * gi_96[k]
                  - f_258 * gi_105[k]
                  + f_259 * gi_107[k]
                  + f_222 * gi_140[k]
                  - f_221 * gi_143[k]
                  - f_223 * gi_145[k]
                  - f_221 * gi_150[k]
                  + f_261 * gi_152[k]
                  + f_222 * gi_161[k]
                  - f_223 * gi_163[k]
                  - f_254 * gi_280[k]
                  + f_255 * gi_283[k]
                  + f_256 * gi_285[k]
                  + f_255 * gi_290[k]
                  - f_257 * gi_292[k]
                  - f_254 * gi_301[k]
                  + f_256 * gi_303[k]
                  + f_222 * gi_336[k]
                  - f_221 * gi_339[k]
                  - f_223 * gi_341[k]
                  - f_221 * gi_346[k]
                  + f_261 * gi_348[k]
                  + f_222 * gi_357[k]
                  - f_223 * gi_359[k]
                  - f_262 * gi_392[k]
                  + f_263 * gi_395[k]
                  + f_264 * gi_397[k]
                  + f_263 * gi_402[k]
                  - f_265 * gi_404[k]
                  - f_262 * gi_413[k]
                  + f_264 * gi_415[k];
    }

#pragma omp simd aligned(gi_2, gi_7, gi_16, gi_86, gi_91, gi_100, gi_142, gi_147, gi_156, \
                         gi_282, gi_287, gi_296, gi_338, gi_343, gi_352, gi_394, gi_399, \
                         gi_408 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_63[k] = f_211 * gi_2[k]
                  - f_210 * gi_7[k]
                  + f_209 * gi_16[k]
                  + f_213 * gi_86[k]
                  - f_212 * gi_91[k]
                  + f_210 * gi_100[k]
                  - f_216 * gi_142[k]
                  + f_215 * gi_147[k]
                  - f_214 * gi_156[k]
                  + f_211 * gi_282[k]
                  - f_210 * gi_287[k]
                  + f_209 * gi_296[k]
                  - f_216 * gi_338[k]
                  + f_215 * gi_343[k]
                  - f_214 * gi_352[k]
                  + f_219 * gi_394[k]
                  - f_218 * gi_399[k]
                  + f_217 * gi_408[k];
    }

#pragma omp simd aligned(gi_0, gi_3, gi_10, gi_21, gi_84, gi_87, gi_94, gi_105, gi_140, \
                         gi_143, gi_150, gi_161, gi_280, gi_283, gi_290, gi_301, gi_336, \
                         gi_339, gi_346, gi_357, gi_392, gi_395, gi_402, \
                         gi_413 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_64[k] = f_266 * gi_0[k]
                  - f_267 * gi_3[k]
                  + f_267 * gi_10[k]
                  - f_266 * gi_21[k]
                  + f_268 * gi_84[k]
                  - f_269 * gi_87[k]
                  + f_269 * gi_94[k]
                  - f_268 * gi_105[k]
                  - f_270 * gi_140[k]
                  + f_271 * gi_143[k]
                  - f_271 * gi_150[k]
                  + f_270 * gi_161[k]
                  + f_266 * gi_280[k]
                  - f_267 * gi_283[k]
                  + f_267 * gi_290[k]
                  - f_266 * gi_301[k]
                  - f_270 * gi_336[k]
                  + f_271 * gi_339[k]
                  - f_271 * gi_346[k]
                  + f_270 * gi_357[k]
                  + f_272 * gi_392[k]
                  - f_204 * gi_395[k]
                  + f_204 * gi_402[k]
                  - f_272 * gi_413[k];
    }

#pragma omp simd aligned(gi_57, gi_62, gi_71, gi_197, gi_202, gi_211, gi_253, gi_258, \
                         gi_267 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_65[k] = -f_146 * gi_57[k]
                  + f_147 * gi_62[k]
                  - f_146 * gi_71[k]
                  - f_146 * gi_197[k]
                  + f_147 * gi_202[k]
                  - f_146 * gi_211[k]
                  + f_148 * gi_253[k]
                  - f_149 * gi_258[k]
                  + f_148 * gi_267[k];
    }

#pragma omp simd aligned(gi_60, gi_67, gi_78, gi_200, gi_207, gi_218, gi_256, gi_263, \
                         gi_274 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_66[k] = -f_150 * gi_60[k]
                  + f_151 * gi_67[k]
                  - f_152 * gi_78[k]
                  - f_150 * gi_200[k]
                  + f_151 * gi_207[k]
                  - f_152 * gi_218[k]
                  + f_153 * gi_256[k]
                  - f_154 * gi_263[k]
                  + f_155 * gi_274[k];
    }

#pragma omp simd aligned(gi_57, gi_64, gi_71, gi_73, gi_197, gi_204, gi_211, gi_213, gi_253, \
                         gi_260, gi_267, gi_269 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_67[k] = f_156 * gi_57[k]
                  - f_70 * gi_64[k]
                  - f_156 * gi_71[k]
                  + f_70 * gi_73[k]
                  + f_156 * gi_197[k]
                  - f_70 * gi_204[k]
                  - f_156 * gi_211[k]
                  + f_70 * gi_213[k]
                  - f_71 * gi_253[k]
                  + f_157 * gi_260[k]
                  + f_71 * gi_267[k]
                  - f_157 * gi_269[k];
    }

#pragma omp simd aligned(gi_60, gi_67, gi_69, gi_78, gi_80, gi_200, gi_207, gi_209, gi_218, \
                         gi_220, gi_256, gi_263, gi_265, gi_274, \
                         gi_276 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_68[k] = f_158 * gi_60[k]
                  + f_159 * gi_67[k]
                  - f_160 * gi_69[k]
                  - f_161 * gi_78[k]
                  + f_162 * gi_80[k]
                  + f_158 * gi_200[k]
                  + f_159 * gi_207[k]
                  - f_160 * gi_209[k]
                  - f_161 * gi_218[k]
                  + f_162 * gi_220[k]
                  - f_163 * gi_256[k]
                  - f_162 * gi_263[k]
                  + f_164 * gi_265[k]
                  + f_165 * gi_274[k]
                  - f_166 * gi_276[k];
    }

#pragma omp simd aligned(gi_57, gi_62, gi_64, gi_71, gi_73, gi_75, gi_197, gi_202, gi_204, \
                         gi_211, gi_213, gi_215, gi_253, gi_258, gi_260, gi_267, gi_269, \
                         gi_271 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_69[k] = -f_167 * gi_57[k]
                  - f_168 * gi_62[k]
                  + f_169 * gi_64[k]
                  - f_167 * gi_71[k]
                  + f_169 * gi_73[k]
                  - f_169 * gi_75[k]
                  - f_167 * gi_197[k]
                  - f_168 * gi_202[k]
                  + f_169 * gi_204[k]
                  - f_167 * gi_211[k]
                  + f_169 * gi_213[k]
                  - f_169 * gi_215[k]
                  + f_170 * gi_253[k]
                  + f_171 * gi_258[k]
                  - f_172 * gi_260[k]
                  + f_170 * gi_267[k]
                  - f_172 * gi_269[k]
                  + f_172 * gi_271[k];
    }

#pragma omp simd aligned(gi_60, gi_67, gi_69, gi_78, gi_80, gi_82, gi_200, gi_207, gi_209, \
                         gi_218, gi_220, gi_222, gi_256, gi_263, gi_265, gi_274, gi_276, \
                         gi_278 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_70[k] = -f_173 * gi_60[k]
                  - f_174 * gi_67[k]
                  + f_175 * gi_69[k]
                  - f_173 * gi_78[k]
                  + f_175 * gi_80[k]
                  - f_176 * gi_82[k]
                  - f_173 * gi_200[k]
                  - f_174 * gi_207[k]
                  + f_175 * gi_209[k]
                  - f_173 * gi_218[k]
                  + f_175 * gi_220[k]
                  - f_176 * gi_222[k]
                  + f_177 * gi_256[k]
                  + f_178 * gi_263[k]
                  - f_179 * gi_265[k]
                  + f_177 * gi_274[k]
                  - f_179 * gi_276[k]
                  + f_180 * gi_278[k];
    }

#pragma omp simd aligned(gi_56, gi_59, gi_61, gi_66, gi_68, gi_70, gi_77, gi_79, gi_81, gi_83, \
                         gi_196, gi_199, gi_201, gi_206, gi_208, gi_210, gi_217, gi_219, \
                         gi_221, gi_223, gi_252, gi_255, gi_257, gi_262, gi_264, gi_266, \
                         gi_273, gi_275, gi_277, gi_279 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_71[k] = f_181 * gi_56[k]
                  + f_182 * gi_59[k]
                  - f_183 * gi_61[k]
                  + f_182 * gi_66[k]
                  - f_184 * gi_68[k]
                  + f_185 * gi_70[k]
                  + f_181 * gi_77[k]
                  - f_183 * gi_79[k]
                  + f_185 * gi_81[k]
                  - f_186 * gi_83[k]
                  + f_181 * gi_196[k]
                  + f_182 * gi_199[k]
                  - f_183 * gi_201[k]
                  + f_182 * gi_206[k]
                  - f_184 * gi_208[k]
                  + f_185 * gi_210[k]
                  + f_181 * gi_217[k]
                  - f_183 * gi_219[k]
                  + f_185 * gi_221[k]
                  - f_186 * gi_223[k]
                  - f_187 * gi_252[k]
                  - f_188 * gi_255[k]
                  + f_185 * gi_257[k]
                  - f_188 * gi_262[k]
                  + f_189 * gi_264[k]
                  - f_190 * gi_266[k]
                  - f_187 * gi_273[k]
                  + f_185 * gi_275[k]
                  - f_190 * gi_277[k]
                  + f_191 * gi_279[k];
    }

#pragma omp simd aligned(gi_58, gi_63, gi_65, gi_72, gi_74, gi_76, gi_198, gi_203, gi_205, \
                         gi_212, gi_214, gi_216, gi_254, gi_259, gi_261, gi_268, gi_270, \
                         gi_272 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_72[k] = -f_173 * gi_58[k]
                  - f_174 * gi_63[k]
                  + f_175 * gi_65[k]
                  - f_173 * gi_72[k]
                  + f_175 * gi_74[k]
                  - f_176 * gi_76[k]
                  - f_173 * gi_198[k]
                  - f_174 * gi_203[k]
                  + f_175 * gi_205[k]
                  - f_173 * gi_212[k]
                  + f_175 * gi_214[k]
                  - f_176 * gi_216[k]
                  + f_177 * gi_254[k]
                  + f_178 * gi_259[k]
                  - f_179 * gi_261[k]
                  + f_177 * gi_268[k]
                  - f_179 * gi_270[k]
                  + f_180 * gi_272[k];
    }

#pragma omp simd aligned(gi_56, gi_59, gi_61, gi_66, gi_70, gi_77, gi_79, gi_81, gi_196, \
                         gi_199, gi_201, gi_206, gi_210, gi_217, gi_219, gi_221, gi_252, \
                         gi_255, gi_257, gi_262, gi_266, gi_273, gi_275, \
                         gi_277 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_73[k] = -f_192 * gi_56[k]
                  - f_192 * gi_59[k]
                  + f_162 * gi_61[k]
                  + f_192 * gi_66[k]
                  - f_162 * gi_70[k]
                  + f_192 * gi_77[k]
                  - f_162 * gi_79[k]
                  + f_162 * gi_81[k]
                  - f_192 * gi_196[k]
                  - f_192 * gi_199[k]
                  + f_162 * gi_201[k]
                  + f_192 * gi_206[k]
                  - f_162 * gi_210[k]
                  + f_192 * gi_217[k]
                  - f_162 * gi_219[k]
                  + f_162 * gi_221[k]
                  + f_193 * gi_252[k]
                  + f_193 * gi_255[k]
                  - f_166 * gi_257[k]
                  - f_193 * gi_262[k]
                  + f_166 * gi_266[k]
                  - f_193 * gi_273[k]
                  + f_166 * gi_275[k]
                  - f_166 * gi_277[k];
    }

#pragma omp simd aligned(gi_58, gi_63, gi_65, gi_72, gi_74, gi_198, gi_203, gi_205, gi_212, \
                         gi_214, gi_254, gi_259, gi_261, gi_268, \
                         gi_270 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_74[k] = f_161 * gi_58[k]
                  - f_159 * gi_63[k]
                  - f_162 * gi_65[k]
                  - f_158 * gi_72[k]
                  + f_160 * gi_74[k]
                  + f_161 * gi_198[k]
                  - f_159 * gi_203[k]
                  - f_162 * gi_205[k]
                  - f_158 * gi_212[k]
                  + f_160 * gi_214[k]
                  - f_165 * gi_254[k]
                  + f_162 * gi_259[k]
                  + f_166 * gi_261[k]
                  + f_163 * gi_268[k]
                  - f_164 * gi_270[k];
    }

#pragma omp simd aligned(gi_56, gi_59, gi_61, gi_66, gi_68, gi_77, gi_79, gi_196, gi_199, \
                         gi_201, gi_206, gi_208, gi_217, gi_219, gi_252, gi_255, gi_257, \
                         gi_262, gi_264, gi_273, gi_275 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_75[k] = f_194 * gi_56[k]
                  - f_67 * gi_59[k]
                  - f_73 * gi_61[k]
                  - f_67 * gi_66[k]
                  + f_69 * gi_68[k]
                  + f_194 * gi_77[k]
                  - f_73 * gi_79[k]
                  + f_194 * gi_196[k]
                  - f_67 * gi_199[k]
                  - f_73 * gi_201[k]
                  - f_67 * gi_206[k]
                  + f_69 * gi_208[k]
                  + f_194 * gi_217[k]
                  - f_73 * gi_219[k]
                  - f_195 * gi_252[k]
                  + f_196 * gi_255[k]
                  + f_75 * gi_257[k]
                  + f_196 * gi_262[k]
                  - f_197 * gi_264[k]
                  - f_195 * gi_273[k]
                  + f_75 * gi_275[k];
    }

#pragma omp simd aligned(gi_58, gi_63, gi_72, gi_198, gi_203, gi_212, gi_254, gi_259, \
                         gi_268 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_76[k] = -f_152 * gi_58[k]
                  + f_151 * gi_63[k]
                  - f_150 * gi_72[k]
                  - f_152 * gi_198[k]
                  + f_151 * gi_203[k]
                  - f_150 * gi_212[k]
                  + f_155 * gi_254[k]
                  - f_154 * gi_259[k]
                  + f_153 * gi_268[k];
    }

#pragma omp simd aligned(gi_56, gi_59, gi_66, gi_77, gi_196, gi_199, gi_206, gi_217, gi_252, \
                         gi_255, gi_262, gi_273 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_77[k] = -f_198 * gi_56[k]
                  + f_199 * gi_59[k]
                  - f_199 * gi_66[k]
                  + f_198 * gi_77[k]
                  - f_198 * gi_196[k]
                  + f_199 * gi_199[k]
                  - f_199 * gi_206[k]
                  + f_198 * gi_217[k]
                  + f_200 * gi_252[k]
                  - f_147 * gi_255[k]
                  + f_147 * gi_262[k]
                  - f_200 * gi_273[k];
    }

#pragma omp simd aligned(gi_1, gi_6, gi_15, gi_141, gi_146, gi_155, gi_281, gi_286, gi_295, \
                         gi_337, gi_342, gi_351 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_78[k] = -f_273 * gi_1[k]
                  + f_274 * gi_6[k]
                  - f_273 * gi_15[k]
                  + f_275 * gi_141[k]
                  - f_276 * gi_146[k]
                  + f_275 * gi_155[k]
                  + f_273 * gi_281[k]
                  - f_274 * gi_286[k]
                  + f_273 * gi_295[k]
                  - f_275 * gi_337[k]
                  + f_276 * gi_342[k]
                  - f_275 * gi_351[k];
    }

#pragma omp simd aligned(gi_4, gi_11, gi_22, gi_144, gi_151, gi_162, gi_284, gi_291, gi_302, \
                         gi_340, gi_347, gi_358 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_79[k] = -f_277 * gi_4[k]
                  + f_95 * gi_11[k]
                  - f_278 * gi_22[k]
                  + f_279 * gi_144[k]
                  - f_98 * gi_151[k]
                  + f_280 * gi_162[k]
                  + f_277 * gi_284[k]
                  - f_95 * gi_291[k]
                  + f_278 * gi_302[k]
                  - f_279 * gi_340[k]
                  + f_98 * gi_347[k]
                  - f_280 * gi_358[k];
    }

#pragma omp simd aligned(gi_1, gi_8, gi_15, gi_17, gi_141, gi_148, gi_155, gi_157, gi_281, \
                         gi_288, gi_295, gi_297, gi_337, gi_344, gi_351, \
                         gi_353 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_80[k] = f_281 * gi_1[k]
                  - f_282 * gi_8[k]
                  - f_281 * gi_15[k]
                  + f_282 * gi_17[k]
                  - f_283 * gi_141[k]
                  + f_284 * gi_148[k]
                  + f_283 * gi_155[k]
                  - f_284 * gi_157[k]
                  - f_281 * gi_281[k]
                  + f_282 * gi_288[k]
                  + f_281 * gi_295[k]
                  - f_282 * gi_297[k]
                  + f_283 * gi_337[k]
                  - f_284 * gi_344[k]
                  - f_283 * gi_351[k]
                  + f_284 * gi_353[k];
    }

#pragma omp simd aligned(gi_4, gi_11, gi_13, gi_22, gi_24, gi_144, gi_151, gi_153, gi_162, \
                         gi_164, gi_284, gi_291, gi_293, gi_302, gi_304, gi_340, gi_347, \
                         gi_349, gi_358, gi_360 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_81[k] = f_285 * gi_4[k]
                  + f_107 * gi_11[k]
                  - f_117 * gi_13[k]
                  - f_286 * gi_22[k]
                  + f_287 * gi_24[k]
                  - f_288 * gi_144[k]
                  - f_112 * gi_151[k]
                  + f_289 * gi_153[k]
                  + f_104 * gi_162[k]
                  - f_106 * gi_164[k]
                  - f_285 * gi_284[k]
                  - f_107 * gi_291[k]
                  + f_117 * gi_293[k]
                  + f_286 * gi_302[k]
                  - f_287 * gi_304[k]
                  + f_288 * gi_340[k]
                  + f_112 * gi_347[k]
                  - f_289 * gi_349[k]
                  - f_104 * gi_358[k]
                  + f_106 * gi_360[k];
    }

#pragma omp simd aligned(gi_1, gi_6, gi_8, gi_15, gi_17, gi_19, gi_141, gi_146, gi_148, \
                         gi_155, gi_157, gi_159, gi_281, gi_286, gi_288, gi_295, gi_297, \
                         gi_299, gi_337, gi_342, gi_344, gi_351, gi_353, \
                         gi_355 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_82[k] = -f_138 * gi_1[k]
                  - f_114 * gi_6[k]
                  + f_108 * gi_8[k]
                  - f_138 * gi_15[k]
                  + f_108 * gi_17[k]
                  - f_108 * gi_19[k]
                  + f_107 * gi_141[k]
                  + f_105 * gi_146[k]
                  - f_113 * gi_148[k]
                  + f_107 * gi_155[k]
                  - f_113 * gi_157[k]
                  + f_113 * gi_159[k]
                  + f_138 * gi_281[k]
                  + f_114 * gi_286[k]
                  - f_108 * gi_288[k]
                  + f_138 * gi_295[k]
                  - f_108 * gi_297[k]
                  + f_108 * gi_299[k]
                  - f_107 * gi_337[k]
                  - f_105 * gi_342[k]
                  + f_113 * gi_344[k]
                  - f_107 * gi_351[k]
                  + f_113 * gi_353[k]
                  - f_113 * gi_355[k];
    }

#pragma omp simd aligned(gi_4, gi_11, gi_13, gi_22, gi_24, gi_26, gi_144, gi_151, gi_153, \
                         gi_162, gi_164, gi_166, gi_284, gi_291, gi_293, gi_302, gi_304, \
                         gi_306, gi_340, gi_347, gi_349, gi_358, gi_360, \
                         gi_362 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_83[k] = -f_290 * gi_4[k]
                  - f_119 * gi_11[k]
                  + f_120 * gi_13[k]
                  - f_290 * gi_22[k]
                  + f_120 * gi_24[k]
                  - f_291 * gi_26[k]
                  + f_292 * gi_144[k]
                  + f_123 * gi_151[k]
                  - f_124 * gi_153[k]
                  + f_292 * gi_162[k]
                  - f_124 * gi_164[k]
                  + f_293 * gi_166[k]
                  + f_290 * gi_284[k]
                  + f_119 * gi_291[k]
                  - f_120 * gi_293[k]
                  + f_290 * gi_302[k]
                  - f_120 * gi_304[k]
                  + f_291 * gi_306[k]
                  - f_292 * gi_340[k]
                  - f_123 * gi_347[k]
                  + f_124 * gi_349[k]
                  - f_292 * gi_358[k]
                  + f_124 * gi_360[k]
                  - f_293 * gi_362[k];
    }

#pragma omp simd aligned(gi_0, gi_3, gi_5, gi_10, gi_12, gi_14, gi_21, gi_23, gi_25, gi_27, \
                         gi_140, gi_143, gi_145, gi_150, gi_152, gi_154, gi_161, gi_163, \
                         gi_165, gi_167, gi_280, gi_283, gi_285, gi_290, gi_292, gi_294, \
                         gi_301, gi_303, gi_305, gi_307, gi_336, gi_339, gi_341, gi_346, \
                         gi_348, gi_350, gi_357, gi_359, gi_361, \
                         gi_363 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_84[k] = f_294 * gi_0[k]
                  + f_295 * gi_3[k]
                  - f_296 * gi_5[k]
                  + f_295 * gi_10[k]
                  - f_129 * gi_12[k]
                  + f_297 * gi_14[k]
                  + f_294 * gi_21[k]
                  - f_296 * gi_23[k]
                  + f_297 * gi_25[k]
                  - f_298 * gi_27[k]
                  - f_128 * gi_140[k]
                  - f_296 * gi_143[k]
                  + f_299 * gi_145[k]
                  - f_296 * gi_150[k]
                  + f_134 * gi_152[k]
                  - f_300 * gi_154[k]
                  - f_128 * gi_161[k]
                  + f_299 * gi_163[k]
                  - f_300 * gi_165[k]
                  + f_301 * gi_167[k]
                  - f_294 * gi_280[k]
                  - f_295 * gi_283[k]
                  + f_296 * gi_285[k]
                  - f_295 * gi_290[k]
                  + f_129 * gi_292[k]
                  - f_297 * gi_294[k]
                  - f_294 * gi_301[k]
                  + f_296 * gi_303[k]
                  - f_297 * gi_305[k]
                  + f_298 * gi_307[k]
                  + f_128 * gi_336[k]
                  + f_296 * gi_339[k]
                  - f_299 * gi_341[k]
                  + f_296 * gi_346[k]
                  - f_134 * gi_348[k]
                  + f_300 * gi_350[k]
                  + f_128 * gi_357[k]
                  - f_299 * gi_359[k]
                  + f_300 * gi_361[k]
                  - f_301 * gi_363[k];
    }

#pragma omp simd aligned(gi_2, gi_7, gi_9, gi_16, gi_18, gi_20, gi_142, gi_147, gi_149, \
                         gi_156, gi_158, gi_160, gi_282, gi_287, gi_289, gi_296, gi_298, \
                         gi_300, gi_338, gi_343, gi_345, gi_352, gi_354, \
                         gi_356 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_85[k] = -f_290 * gi_2[k]
                  - f_119 * gi_7[k]
                  + f_120 * gi_9[k]
                  - f_290 * gi_16[k]
                  + f_120 * gi_18[k]
                  - f_291 * gi_20[k]
                  + f_292 * gi_142[k]
                  + f_123 * gi_147[k]
                  - f_124 * gi_149[k]
                  + f_292 * gi_156[k]
                  - f_124 * gi_158[k]
                  + f_293 * gi_160[k]
                  + f_290 * gi_282[k]
                  + f_119 * gi_287[k]
                  - f_120 * gi_289[k]
                  + f_290 * gi_296[k]
                  - f_120 * gi_298[k]
                  + f_291 * gi_300[k]
                  - f_292 * gi_338[k]
                  - f_123 * gi_343[k]
                  + f_124 * gi_345[k]
                  - f_292 * gi_352[k]
                  + f_124 * gi_354[k]
                  - f_293 * gi_356[k];
    }

#pragma omp simd aligned(gi_0, gi_3, gi_5, gi_10, gi_14, gi_21, gi_23, gi_25, gi_140, gi_143, \
                         gi_145, gi_150, gi_154, gi_161, gi_163, gi_165, gi_280, gi_283, \
                         gi_285, gi_290, gi_294, gi_301, gi_303, gi_305, gi_336, gi_339, \
                         gi_341, gi_346, gi_350, gi_357, gi_359, \
                         gi_361 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_86[k] = -f_302 * gi_0[k]
                  - f_302 * gi_3[k]
                  + f_287 * gi_5[k]
                  + f_302 * gi_10[k]
                  - f_287 * gi_14[k]
                  + f_302 * gi_21[k]
                  - f_287 * gi_23[k]
                  + f_287 * gi_25[k]
                  + f_286 * gi_140[k]
                  + f_286 * gi_143[k]
                  - f_106 * gi_145[k]
                  - f_286 * gi_150[k]
                  + f_106 * gi_154[k]
                  - f_286 * gi_161[k]
                  + f_106 * gi_163[k]
                  - f_106 * gi_165[k]
                  + f_302 * gi_280[k]
                  + f_302 * gi_283[k]
                  - f_287 * gi_285[k]
                  - f_302 * gi_290[k]
                  + f_287 * gi_294[k]
                  - f_302 * gi_301[k]
                  + f_287 * gi_303[k]
                  - f_287 * gi_305[k]
                  - f_286 * gi_336[k]
                  - f_286 * gi_339[k]
                  + f_106 * gi_341[k]
                  + f_286 * gi_346[k]
                  - f_106 * gi_350[k]
                  + f_286 * gi_357[k]
                  - f_106 * gi_359[k]
                  + f_106 * gi_361[k];
    }

#pragma omp simd aligned(gi_2, gi_7, gi_9, gi_16, gi_18, gi_142, gi_147, gi_149, gi_156, \
                         gi_158, gi_282, gi_287, gi_289, gi_296, gi_298, gi_338, gi_343, \
                         gi_345, gi_352, gi_354 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_87[k] = f_286 * gi_2[k]
                  - f_107 * gi_7[k]
                  - f_287 * gi_9[k]
                  - f_285 * gi_16[k]
                  + f_117 * gi_18[k]
                  - f_104 * gi_142[k]
                  + f_112 * gi_147[k]
                  + f_106 * gi_149[k]
                  + f_288 * gi_156[k]
                  - f_289 * gi_158[k]
                  - f_286 * gi_282[k]
                  + f_107 * gi_287[k]
                  + f_287 * gi_289[k]
                  + f_285 * gi_296[k]
                  - f_117 * gi_298[k]
                  + f_104 * gi_338[k]
                  - f_112 * gi_343[k]
                  - f_106 * gi_345[k]
                  - f_288 * gi_352[k]
                  + f_289 * gi_354[k];
    }

#pragma omp simd aligned(gi_0, gi_3, gi_5, gi_10, gi_12, gi_21, gi_23, gi_140, gi_143, gi_145, \
                         gi_150, gi_152, gi_161, gi_163, gi_280, gi_283, gi_285, gi_290, \
                         gi_292, gi_301, gi_303, gi_336, gi_339, gi_341, gi_346, gi_348, \
                         gi_357, gi_359 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_88[k] = f_303 * gi_0[k]
                  - f_304 * gi_3[k]
                  - f_20 * gi_5[k]
                  - f_304 * gi_10[k]
                  + f_21 * gi_12[k]
                  + f_303 * gi_21[k]
                  - f_20 * gi_23[k]
                  - f_305 * gi_140[k]
                  + f_306 * gi_143[k]
                  + f_21 * gi_145[k]
                  + f_306 * gi_150[k]
                  - f_307 * gi_152[k]
                  - f_305 * gi_161[k]
                  + f_21 * gi_163[k]
                  - f_303 * gi_280[k]
                  + f_304 * gi_283[k]
                  + f_20 * gi_285[k]
                  + f_304 * gi_290[k]
                  - f_21 * gi_292[k]
                  - f_303 * gi_301[k]
                  + f_20 * gi_303[k]
                  + f_305 * gi_336[k]
                  - f_306 * gi_339[k]
                  - f_21 * gi_341[k]
                  - f_306 * gi_346[k]
                  + f_307 * gi_348[k]
                  + f_305 * gi_357[k]
                  - f_21 * gi_359[k];
    }

#pragma omp simd aligned(gi_2, gi_7, gi_16, gi_142, gi_147, gi_156, gi_282, gi_287, gi_296, \
                         gi_338, gi_343, gi_352 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_89[k] = -f_278 * gi_2[k]
                  + f_95 * gi_7[k]
                  - f_277 * gi_16[k]
                  + f_280 * gi_142[k]
                  - f_98 * gi_147[k]
                  + f_279 * gi_156[k]
                  + f_278 * gi_282[k]
                  - f_95 * gi_287[k]
                  + f_277 * gi_296[k]
                  - f_280 * gi_338[k]
                  + f_98 * gi_343[k]
                  - f_279 * gi_352[k];
    }

#pragma omp simd aligned(gi_0, gi_3, gi_10, gi_21, gi_140, gi_143, gi_150, gi_161, gi_280, \
                         gi_283, gi_290, gi_301, gi_336, gi_339, gi_346, \
                         gi_357 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_90[k] = -f_308 * gi_0[k]
                  + f_309 * gi_3[k]
                  - f_309 * gi_10[k]
                  + f_308 * gi_21[k]
                  + f_273 * gi_140[k]
                  - f_310 * gi_143[k]
                  + f_310 * gi_150[k]
                  - f_273 * gi_161[k]
                  + f_308 * gi_280[k]
                  - f_309 * gi_283[k]
                  + f_309 * gi_290[k]
                  - f_308 * gi_301[k]
                  - f_273 * gi_336[k]
                  + f_310 * gi_339[k]
                  - f_310 * gi_346[k]
                  + f_273 * gi_357[k];
    }

#pragma omp simd aligned(gi_57, gi_60, gi_62, gi_67, gi_71, gi_78, gi_197, gi_200, gi_202, \
                         gi_207, gi_211, gi_218 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_91[k] = f_34 * gi_57[k]
                  - f_35 * gi_62[k]
                  + f_34 * gi_71[k]
                  - f_32 * gi_197[k]
                  + f_33 * gi_202[k]
                  - f_32 * gi_211[k];

        g_92[k] = f_39 * gi_60[k]
                  - f_40 * gi_67[k]
                  + f_41 * gi_78[k]
                  - f_36 * gi_200[k]
                  + f_37 * gi_207[k]
                  - f_38 * gi_218[k];
    }

#pragma omp simd aligned(gi_57, gi_64, gi_71, gi_73, gi_197, gi_204, gi_211, \
                         gi_213 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_93[k] = -f_44 * gi_57[k]
                  + f_45 * gi_64[k]
                  + f_44 * gi_71[k]
                  - f_45 * gi_73[k]
                  + f_42 * gi_197[k]
                  - f_43 * gi_204[k]
                  - f_42 * gi_211[k]
                  + f_43 * gi_213[k];
    }

#pragma omp simd aligned(gi_60, gi_67, gi_69, gi_78, gi_80, gi_200, gi_207, gi_209, gi_218, \
                         gi_220 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_94[k] = -f_49 * gi_60[k]
                  - f_51 * gi_67[k]
                  + f_50 * gi_69[k]
                  + f_52 * gi_78[k]
                  - f_53 * gi_80[k]
                  + f_46 * gi_200[k]
                  + f_47 * gi_207[k]
                  - f_48 * gi_209[k]
                  - f_49 * gi_218[k]
                  + f_50 * gi_220[k];
    }

#pragma omp simd aligned(gi_57, gi_62, gi_64, gi_71, gi_73, gi_75, gi_197, gi_202, gi_204, \
                         gi_211, gi_213, gi_215 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_95[k] = f_55 * gi_57[k]
                  + f_56 * gi_62[k]
                  - f_57 * gi_64[k]
                  + f_55 * gi_71[k]
                  - f_57 * gi_73[k]
                  + f_57 * gi_75[k]
                  - f_52 * gi_197[k]
                  - f_51 * gi_202[k]
                  + f_54 * gi_204[k]
                  - f_52 * gi_211[k]
                  + f_54 * gi_213[k]
                  - f_54 * gi_215[k];
    }

#pragma omp simd aligned(gi_60, gi_67, gi_69, gi_78, gi_80, gi_82, gi_200, gi_207, gi_209, \
                         gi_218, gi_220, gi_222 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_96[k] = f_62 * gi_60[k]
                  + f_63 * gi_67[k]
                  - f_64 * gi_69[k]
                  + f_62 * gi_78[k]
                  - f_64 * gi_80[k]
                  + f_65 * gi_82[k]
                  - f_58 * gi_200[k]
                  - f_59 * gi_207[k]
                  + f_60 * gi_209[k]
                  - f_58 * gi_218[k]
                  + f_60 * gi_220[k]
                  - f_61 * gi_222[k];
    }

#pragma omp simd aligned(gi_56, gi_59, gi_61, gi_66, gi_68, gi_70, gi_77, gi_79, gi_81, gi_83, \
                         gi_196, gi_199, gi_201, gi_206, gi_208, gi_210, gi_217, gi_219, \
                         gi_221, gi_223 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_97[k] = -f_72 * gi_56[k]
                  - f_66 * gi_59[k]
                  + f_73 * gi_61[k]
                  - f_66 * gi_66[k]
                  + f_74 * gi_68[k]
                  - f_75 * gi_70[k]
                  - f_72 * gi_77[k]
                  + f_73 * gi_79[k]
                  - f_75 * gi_81[k]
                  + f_76 * gi_83[k]
                  + f_66 * gi_196[k]
                  + f_67 * gi_199[k]
                  - f_68 * gi_201[k]
                  + f_67 * gi_206[k]
                  - f_69 * gi_208[k]
                  + f_70 * gi_210[k]
                  + f_66 * gi_217[k]
                  - f_68 * gi_219[k]
                  + f_70 * gi_221[k]
                  - f_71 * gi_223[k];
    }

#pragma omp simd aligned(gi_58, gi_63, gi_65, gi_72, gi_74, gi_76, gi_198, gi_203, gi_205, \
                         gi_212, gi_214, gi_216 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_98[k] = f_62 * gi_58[k]
                  + f_63 * gi_63[k]
                  - f_64 * gi_65[k]
                  + f_62 * gi_72[k]
                  - f_64 * gi_74[k]
                  + f_65 * gi_76[k]
                  - f_58 * gi_198[k]
                  - f_59 * gi_203[k]
                  + f_60 * gi_205[k]
                  - f_58 * gi_212[k]
                  + f_60 * gi_214[k]
                  - f_61 * gi_216[k];
    }

#pragma omp simd aligned(gi_56, gi_59, gi_61, gi_66, gi_70, gi_77, gi_79, gi_81, gi_196, \
                         gi_199, gi_201, gi_206, gi_210, gi_217, gi_219, \
                         gi_221 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_99[k] = f_78 * gi_56[k]
                  + f_78 * gi_59[k]
                  - f_53 * gi_61[k]
                  - f_78 * gi_66[k]
                  + f_53 * gi_70[k]
                  - f_78 * gi_77[k]
                  + f_53 * gi_79[k]
                  - f_53 * gi_81[k]
                  - f_77 * gi_196[k]
                  - f_77 * gi_199[k]
                  + f_50 * gi_201[k]
                  + f_77 * gi_206[k]
                  - f_50 * gi_210[k]
                  + f_77 * gi_217[k]
                  - f_50 * gi_219[k]
                  + f_50 * gi_221[k];
    }

#pragma omp simd aligned(gi_58, gi_63, gi_65, gi_72, gi_74, gi_198, gi_203, gi_205, gi_212, \
                         gi_214 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_100[k] = -f_52 * gi_58[k]
                   + f_51 * gi_63[k]
                   + f_53 * gi_65[k]
                   + f_49 * gi_72[k]
                   - f_50 * gi_74[k]
                   + f_49 * gi_198[k]
                   - f_47 * gi_203[k]
                   - f_50 * gi_205[k]
                   - f_46 * gi_212[k]
                   + f_48 * gi_214[k];
    }

#pragma omp simd aligned(gi_56, gi_59, gi_61, gi_66, gi_68, gi_77, gi_79, gi_196, gi_199, \
                         gi_201, gi_206, gi_208, gi_217, gi_219 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_101[k] = -f_83 * gi_56[k]
                   + f_84 * gi_59[k]
                   + f_85 * gi_61[k]
                   + f_84 * gi_66[k]
                   - f_86 * gi_68[k]
                   - f_83 * gi_77[k]
                   + f_85 * gi_79[k]
                   + f_79 * gi_196[k]
                   - f_80 * gi_199[k]
                   - f_81 * gi_201[k]
                   - f_80 * gi_206[k]
                   + f_82 * gi_208[k]
                   + f_79 * gi_217[k]
                   - f_81 * gi_219[k];
    }

#pragma omp simd aligned(gi_56, gi_58, gi_59, gi_63, gi_66, gi_72, gi_77, gi_196, gi_198, \
                         gi_199, gi_203, gi_206, gi_212, gi_217 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_102[k] = f_41 * gi_58[k]
                   - f_40 * gi_63[k]
                   + f_39 * gi_72[k]
                   - f_38 * gi_198[k]
                   + f_37 * gi_203[k]
                   - f_36 * gi_212[k];

        g_103[k] = f_89 * gi_56[k]
                   - f_90 * gi_59[k]
                   + f_90 * gi_66[k]
                   - f_89 * gi_77[k]
                   - f_87 * gi_196[k]
                   + f_88 * gi_199[k]
                   - f_88 * gi_206[k]
                   + f_87 * gi_217[k];
    }

#pragma omp simd aligned(gi_1, gi_6, gi_15, gi_85, gi_90, gi_99, gi_281, gi_286, \
                         gi_295 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_104[k] = f_311 * gi_1[k]
                   - f_312 * gi_6[k]
                   + f_311 * gi_15[k]
                   - f_313 * gi_85[k]
                   + f_314 * gi_90[k]
                   - f_313 * gi_99[k]
                   + f_311 * gi_281[k]
                   - f_312 * gi_286[k]
                   + f_311 * gi_295[k];
    }

#pragma omp simd aligned(gi_4, gi_11, gi_22, gi_88, gi_95, gi_106, gi_284, gi_291, \
                         gi_302 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_105[k] = f_315 * gi_4[k]
                   - f_316 * gi_11[k]
                   + f_317 * gi_22[k]
                   - f_318 * gi_88[k]
                   + f_319 * gi_95[k]
                   - f_320 * gi_106[k]
                   + f_315 * gi_284[k]
                   - f_316 * gi_291[k]
                   + f_317 * gi_302[k];
    }

#pragma omp simd aligned(gi_1, gi_8, gi_15, gi_17, gi_85, gi_92, gi_99, gi_101, gi_281, \
                         gi_288, gi_295, gi_297 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_106[k] = -f_26 * gi_1[k]
                   + f_28 * gi_8[k]
                   + f_26 * gi_15[k]
                   - f_28 * gi_17[k]
                   + f_321 * gi_85[k]
                   - f_29 * gi_92[k]
                   - f_321 * gi_99[k]
                   + f_29 * gi_101[k]
                   - f_26 * gi_281[k]
                   + f_28 * gi_288[k]
                   + f_26 * gi_295[k]
                   - f_28 * gi_297[k];
    }

#pragma omp simd aligned(gi_4, gi_11, gi_13, gi_22, gi_24, gi_88, gi_95, gi_97, gi_106, \
                         gi_108, gi_284, gi_291, gi_293, gi_302, \
                         gi_304 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_107[k] = -f_322 * gi_4[k]
                   - f_323 * gi_11[k]
                   + f_8 * gi_13[k]
                   + f_324 * gi_22[k]
                   - f_13 * gi_24[k]
                   + f_325 * gi_88[k]
                   + f_7 * gi_95[k]
                   - f_326 * gi_97[k]
                   - f_327 * gi_106[k]
                   + f_328 * gi_108[k]
                   - f_322 * gi_284[k]
                   - f_323 * gi_291[k]
                   + f_8 * gi_293[k]
                   + f_324 * gi_302[k]
                   - f_13 * gi_304[k];
    }

#pragma omp simd aligned(gi_1, gi_6, gi_8, gi_15, gi_17, gi_19, gi_85, gi_90, gi_92, gi_99, \
                         gi_101, gi_103, gi_281, gi_286, gi_288, gi_295, gi_297, \
                         gi_299 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_108[k] = f_329 * gi_1[k]
                   + f_25 * gi_6[k]
                   - f_330 * gi_8[k]
                   + f_329 * gi_15[k]
                   - f_330 * gi_17[k]
                   + f_330 * gi_19[k]
                   - f_323 * gi_85[k]
                   - f_10 * gi_90[k]
                   + f_9 * gi_92[k]
                   - f_323 * gi_99[k]
                   + f_9 * gi_101[k]
                   - f_9 * gi_103[k]
                   + f_329 * gi_281[k]
                   + f_25 * gi_286[k]
                   - f_330 * gi_288[k]
                   + f_329 * gi_295[k]
                   - f_330 * gi_297[k]
                   + f_330 * gi_299[k];
    }

#pragma omp simd aligned(gi_4, gi_11, gi_13, gi_22, gi_24, gi_26, gi_88, gi_95, gi_97, gi_106, \
                         gi_108, gi_110, gi_284, gi_291, gi_293, gi_302, gi_304, \
                         gi_306 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_109[k] = f_331 * gi_4[k]
                   + f_332 * gi_11[k]
                   - f_15 * gi_13[k]
                   + f_331 * gi_22[k]
                   - f_15 * gi_24[k]
                   + f_333 * gi_26[k]
                   - f_334 * gi_88[k]
                   - f_335 * gi_95[k]
                   + f_336 * gi_97[k]
                   - f_334 * gi_106[k]
                   + f_336 * gi_108[k]
                   - f_337 * gi_110[k]
                   + f_331 * gi_284[k]
                   + f_332 * gi_291[k]
                   - f_15 * gi_293[k]
                   + f_331 * gi_302[k]
                   - f_15 * gi_304[k]
                   + f_333 * gi_306[k];
    }

#pragma omp simd aligned(gi_0, gi_3, gi_5, gi_10, gi_12, gi_14, gi_21, gi_23, gi_25, gi_27, \
                         gi_84, gi_87, gi_89, gi_94, gi_96, gi_98, gi_105, gi_107, gi_109, \
                         gi_111, gi_280, gi_283, gi_285, gi_290, gi_292, gi_294, gi_301, \
                         gi_303, gi_305, gi_307 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_110[k] = -f_338 * gi_0[k]
                   - f_339 * gi_3[k]
                   + f_340 * gi_5[k]
                   - f_339 * gi_10[k]
                   + f_306 * gi_12[k]
                   - f_140 * gi_14[k]
                   - f_338 * gi_21[k]
                   + f_340 * gi_23[k]
                   - f_140 * gi_25[k]
                   + f_341 * gi_27[k]
                   + f_304 * gi_84[k]
                   + f_340 * gi_87[k]
                   - f_342 * gi_89[k]
                   + f_340 * gi_94[k]
                   - f_343 * gi_96[k]
                   + f_22 * gi_98[k]
                   + f_304 * gi_105[k]
                   - f_342 * gi_107[k]
                   + f_22 * gi_109[k]
                   - f_344 * gi_111[k]
                   - f_338 * gi_280[k]
                   - f_339 * gi_283[k]
                   + f_340 * gi_285[k]
                   - f_339 * gi_290[k]
                   + f_306 * gi_292[k]
                   - f_140 * gi_294[k]
                   - f_338 * gi_301[k]
                   + f_340 * gi_303[k]
                   - f_140 * gi_305[k]
                   + f_341 * gi_307[k];
    }

#pragma omp simd aligned(gi_2, gi_7, gi_9, gi_16, gi_18, gi_20, gi_86, gi_91, gi_93, gi_100, \
                         gi_102, gi_104, gi_282, gi_287, gi_289, gi_296, gi_298, \
                         gi_300 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_111[k] = f_331 * gi_2[k]
                   + f_332 * gi_7[k]
                   - f_15 * gi_9[k]
                   + f_331 * gi_16[k]
                   - f_15 * gi_18[k]
                   + f_333 * gi_20[k]
                   - f_334 * gi_86[k]
                   - f_335 * gi_91[k]
                   + f_336 * gi_93[k]
                   - f_334 * gi_100[k]
                   + f_336 * gi_102[k]
                   - f_337 * gi_104[k]
                   + f_331 * gi_282[k]
                   + f_332 * gi_287[k]
                   - f_15 * gi_289[k]
                   + f_331 * gi_296[k]
                   - f_15 * gi_298[k]
                   + f_333 * gi_300[k];
    }

#pragma omp simd aligned(gi_0, gi_3, gi_5, gi_10, gi_14, gi_21, gi_23, gi_25, gi_84, gi_87, \
                         gi_89, gi_94, gi_98, gi_105, gi_107, gi_109, gi_280, gi_283, gi_285, \
                         gi_290, gi_294, gi_301, gi_303, gi_305 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_112[k] = f_345 * gi_0[k]
                   + f_345 * gi_3[k]
                   - f_13 * gi_5[k]
                   - f_345 * gi_10[k]
                   + f_13 * gi_14[k]
                   - f_345 * gi_21[k]
                   + f_13 * gi_23[k]
                   - f_13 * gi_25[k]
                   - f_324 * gi_84[k]
                   - f_324 * gi_87[k]
                   + f_328 * gi_89[k]
                   + f_324 * gi_94[k]
                   - f_328 * gi_98[k]
                   + f_324 * gi_105[k]
                   - f_328 * gi_107[k]
                   + f_328 * gi_109[k]
                   + f_345 * gi_280[k]
                   + f_345 * gi_283[k]
                   - f_13 * gi_285[k]
                   - f_345 * gi_290[k]
                   + f_13 * gi_294[k]
                   - f_345 * gi_301[k]
                   + f_13 * gi_303[k]
                   - f_13 * gi_305[k];
    }

#pragma omp simd aligned(gi_2, gi_7, gi_9, gi_16, gi_18, gi_86, gi_91, gi_93, gi_100, gi_102, \
                         gi_282, gi_287, gi_289, gi_296, gi_298 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_113[k] = -f_324 * gi_2[k]
                   + f_323 * gi_7[k]
                   + f_13 * gi_9[k]
                   + f_322 * gi_16[k]
                   - f_8 * gi_18[k]
                   + f_327 * gi_86[k]
                   - f_7 * gi_91[k]
                   - f_328 * gi_93[k]
                   - f_325 * gi_100[k]
                   + f_326 * gi_102[k]
                   - f_324 * gi_282[k]
                   + f_323 * gi_287[k]
                   + f_13 * gi_289[k]
                   + f_322 * gi_296[k]
                   - f_8 * gi_298[k];
    }

#pragma omp simd aligned(gi_0, gi_3, gi_5, gi_10, gi_12, gi_21, gi_23, gi_84, gi_87, gi_89, \
                         gi_94, gi_96, gi_105, gi_107, gi_280, gi_283, gi_285, gi_290, gi_292, \
                         gi_301, gi_303 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_114[k] = -f_346 * gi_0[k]
                   + f_347 * gi_3[k]
                   + f_348 * gi_5[k]
                   + f_347 * gi_10[k]
                   - f_349 * gi_12[k]
                   - f_346 * gi_21[k]
                   + f_348 * gi_23[k]
                   + f_350 * gi_84[k]
                   - f_351 * gi_87[k]
                   - f_349 * gi_89[k]
                   - f_351 * gi_94[k]
                   + f_352 * gi_96[k]
                   + f_350 * gi_105[k]
                   - f_349 * gi_107[k]
                   - f_346 * gi_280[k]
                   + f_347 * gi_283[k]
                   + f_348 * gi_285[k]
                   + f_347 * gi_290[k]
                   - f_349 * gi_292[k]
                   - f_346 * gi_301[k]
                   + f_348 * gi_303[k];
    }

#pragma omp simd aligned(gi_2, gi_7, gi_16, gi_86, gi_91, gi_100, gi_282, gi_287, \
                         gi_296 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_115[k] = f_317 * gi_2[k]
                   - f_316 * gi_7[k]
                   + f_315 * gi_16[k]
                   - f_320 * gi_86[k]
                   + f_319 * gi_91[k]
                   - f_318 * gi_100[k]
                   + f_317 * gi_282[k]
                   - f_316 * gi_287[k]
                   + f_315 * gi_296[k];
    }

#pragma omp simd aligned(gi_0, gi_3, gi_10, gi_21, gi_84, gi_87, gi_94, gi_105, gi_280, \
                         gi_283, gi_290, gi_301 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_116[k] = f_353 * gi_0[k]
                   - f_354 * gi_3[k]
                   + f_354 * gi_10[k]
                   - f_353 * gi_21[k]
                   - f_311 * gi_84[k]
                   + f_355 * gi_87[k]
                   - f_355 * gi_94[k]
                   + f_311 * gi_105[k]
                   + f_353 * gi_280[k]
                   - f_354 * gi_283[k]
                   + f_354 * gi_290[k]
                   - f_353 * gi_301[k];
    }
}

}  // namespace simdtrf
