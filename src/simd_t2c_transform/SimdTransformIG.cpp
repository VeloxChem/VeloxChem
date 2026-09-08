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


#include "SimdTransformIG.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_ig(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t ig,
             const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 0.65625 * std::sqrt(330.0);
    const auto f_1 = 2.1875 * std::sqrt(330.0);
    const auto f_2 = 1.96875 * std::sqrt(165.0);
    const auto f_3 = 0.65625 * std::sqrt(165.0);
    const auto f_4 = 6.5625 * std::sqrt(165.0);
    const auto f_5 = 2.1875 * std::sqrt(165.0);
    const auto f_6 = 0.09375 * std::sqrt(2310.0);
    const auto f_7 = 0.5625 * std::sqrt(2310.0);
    const auto f_8 = 0.3125 * std::sqrt(2310.0);
    const auto f_9 = 1.875 * std::sqrt(2310.0);
    const auto f_10 = 0.28125 * std::sqrt(1155.0);
    const auto f_11 = 0.375 * std::sqrt(1155.0);
    const auto f_12 = 0.9375 * std::sqrt(1155.0);
    const auto f_13 = 1.25 * std::sqrt(1155.0);
    const auto f_14 = 0.0703125 * std::sqrt(462.0);
    const auto f_15 = 0.140625 * std::sqrt(462.0);
    const auto f_16 = 0.5625 * std::sqrt(462.0);
    const auto f_17 = 0.1875 * std::sqrt(462.0);
    const auto f_18 = 0.234375 * std::sqrt(462.0);
    const auto f_19 = 0.46875 * std::sqrt(462.0);
    const auto f_20 = 1.875 * std::sqrt(462.0);
    const auto f_21 = 0.625 * std::sqrt(462.0);
    const auto f_22 = 0.046875 * std::sqrt(2310.0);
    const auto f_23 = 0.28125 * std::sqrt(2310.0);
    const auto f_24 = 0.15625 * std::sqrt(2310.0);
    const auto f_25 = 0.9375 * std::sqrt(2310.0);
    const auto f_26 = 0.1640625 * std::sqrt(330.0);
    const auto f_27 = 0.984375 * std::sqrt(330.0);
    const auto f_28 = 0.546875 * std::sqrt(330.0);
    const auto f_29 = 3.28125 * std::sqrt(330.0);
    const auto f_30 = 3.28125 * std::sqrt(110.0);
    const auto f_31 = 6.5625 * std::sqrt(110.0);
    const auto f_32 = 0.65625 * std::sqrt(110.0);
    const auto f_33 = 9.84375 * std::sqrt(55.0);
    const auto f_34 = 3.28125 * std::sqrt(55.0);
    const auto f_35 = 19.6875 * std::sqrt(55.0);
    const auto f_36 = 6.5625 * std::sqrt(55.0);
    const auto f_37 = 1.96875 * std::sqrt(55.0);
    const auto f_38 = 0.65625 * std::sqrt(55.0);
    const auto f_39 = 0.46875 * std::sqrt(770.0);
    const auto f_40 = 2.8125 * std::sqrt(770.0);
    const auto f_41 = 0.9375 * std::sqrt(770.0);
    const auto f_42 = 5.625 * std::sqrt(770.0);
    const auto f_43 = 0.09375 * std::sqrt(770.0);
    const auto f_44 = 0.5625 * std::sqrt(770.0);
    const auto f_45 = 1.40625 * std::sqrt(385.0);
    const auto f_46 = 1.875 * std::sqrt(385.0);
    const auto f_47 = 2.8125 * std::sqrt(385.0);
    const auto f_48 = 3.75 * std::sqrt(385.0);
    const auto f_49 = 0.28125 * std::sqrt(385.0);
    const auto f_50 = 0.375 * std::sqrt(385.0);
    const auto f_51 = 0.3515625 * std::sqrt(154.0);
    const auto f_52 = 0.703125 * std::sqrt(154.0);
    const auto f_53 = 2.8125 * std::sqrt(154.0);
    const auto f_54 = 0.9375 * std::sqrt(154.0);
    const auto f_55 = 1.40625 * std::sqrt(154.0);
    const auto f_56 = 5.625 * std::sqrt(154.0);
    const auto f_57 = 1.875 * std::sqrt(154.0);
    const auto f_58 = 0.0703125 * std::sqrt(154.0);
    const auto f_59 = 0.140625 * std::sqrt(154.0);
    const auto f_60 = 0.5625 * std::sqrt(154.0);
    const auto f_61 = 0.1875 * std::sqrt(154.0);
    const auto f_62 = 0.234375 * std::sqrt(770.0);
    const auto f_63 = 1.40625 * std::sqrt(770.0);
    const auto f_64 = 0.046875 * std::sqrt(770.0);
    const auto f_65 = 0.28125 * std::sqrt(770.0);
    const auto f_66 = 0.8203125 * std::sqrt(110.0);
    const auto f_67 = 4.921875 * std::sqrt(110.0);
    const auto f_68 = 1.640625 * std::sqrt(110.0);
    const auto f_69 = 9.84375 * std::sqrt(110.0);
    const auto f_70 = 0.1640625 * std::sqrt(110.0);
    const auto f_71 = 0.984375 * std::sqrt(110.0);
    const auto f_72 = 2.625 * std::sqrt(5.0);
    const auto f_73 = 26.25 * std::sqrt(5.0);
    const auto f_74 = 3.9375 * std::sqrt(10.0);
    const auto f_75 = 1.3125 * std::sqrt(10.0);
    const auto f_76 = 39.375 * std::sqrt(10.0);
    const auto f_77 = 13.125 * std::sqrt(10.0);
    const auto f_78 = 0.375 * std::sqrt(35.0);
    const auto f_79 = 2.25 * std::sqrt(35.0);
    const auto f_80 = 3.75 * std::sqrt(35.0);
    const auto f_81 = 22.5 * std::sqrt(35.0);
    const auto f_82 = 0.5625 * std::sqrt(70.0);
    const auto f_83 = 0.75 * std::sqrt(70.0);
    const auto f_84 = 5.625 * std::sqrt(70.0);
    const auto f_85 = 7.5 * std::sqrt(70.0);
    const auto f_86 = 0.28125 * std::sqrt(7.0);
    const auto f_87 = 0.5625 * std::sqrt(7.0);
    const auto f_88 = 2.25 * std::sqrt(7.0);
    const auto f_89 = 0.75 * std::sqrt(7.0);
    const auto f_90 = 2.8125 * std::sqrt(7.0);
    const auto f_91 = 5.625 * std::sqrt(7.0);
    const auto f_92 = 22.5 * std::sqrt(7.0);
    const auto f_93 = 7.5 * std::sqrt(7.0);
    const auto f_94 = 0.1875 * std::sqrt(35.0);
    const auto f_95 = 1.125 * std::sqrt(35.0);
    const auto f_96 = 1.875 * std::sqrt(35.0);
    const auto f_97 = 11.25 * std::sqrt(35.0);
    const auto f_98 = 0.65625 * std::sqrt(5.0);
    const auto f_99 = 3.9375 * std::sqrt(5.0);
    const auto f_100 = 6.5625 * std::sqrt(5.0);
    const auto f_101 = 39.375 * std::sqrt(5.0);
    const auto f_102 = 9.84375 * std::sqrt(6.0);
    const auto f_103 = 6.5625 * std::sqrt(6.0);
    const auto f_104 = 26.25 * std::sqrt(6.0);
    const auto f_105 = 3.28125 * std::sqrt(6.0);
    const auto f_106 = 8.75 * std::sqrt(6.0);
    const auto f_107 = 29.53125 * std::sqrt(3.0);
    const auto f_108 = 9.84375 * std::sqrt(3.0);
    const auto f_109 = 19.6875 * std::sqrt(3.0);
    const auto f_110 = 6.5625 * std::sqrt(3.0);
    const auto f_111 = 78.75 * std::sqrt(3.0);
    const auto f_112 = 26.25 * std::sqrt(3.0);
    const auto f_113 = 3.28125 * std::sqrt(3.0);
    const auto f_114 = 8.75 * std::sqrt(3.0);
    const auto f_115 = 1.40625 * std::sqrt(42.0);
    const auto f_116 = 8.4375 * std::sqrt(42.0);
    const auto f_117 = 0.9375 * std::sqrt(42.0);
    const auto f_118 = 5.625 * std::sqrt(42.0);
    const auto f_119 = 3.75 * std::sqrt(42.0);
    const auto f_120 = 22.5 * std::sqrt(42.0);
    const auto f_121 = 0.46875 * std::sqrt(42.0);
    const auto f_122 = 2.8125 * std::sqrt(42.0);
    const auto f_123 = 1.25 * std::sqrt(42.0);
    const auto f_124 = 7.5 * std::sqrt(42.0);
    const auto f_125 = 4.21875 * std::sqrt(21.0);
    const auto f_126 = 5.625 * std::sqrt(21.0);
    const auto f_127 = 2.8125 * std::sqrt(21.0);
    const auto f_128 = 3.75 * std::sqrt(21.0);
    const auto f_129 = 11.25 * std::sqrt(21.0);
    const auto f_130 = 15.0 * std::sqrt(21.0);
    const auto f_131 = 1.40625 * std::sqrt(21.0);
    const auto f_132 = 1.875 * std::sqrt(21.0);
    const auto f_133 = 5.0 * std::sqrt(21.0);
    const auto f_134 = 0.2109375 * std::sqrt(210.0);
    const auto f_135 = 0.421875 * std::sqrt(210.0);
    const auto f_136 = 1.6875 * std::sqrt(210.0);
    const auto f_137 = 0.5625 * std::sqrt(210.0);
    const auto f_138 = 0.140625 * std::sqrt(210.0);
    const auto f_139 = 0.28125 * std::sqrt(210.0);
    const auto f_140 = 1.125 * std::sqrt(210.0);
    const auto f_141 = 0.375 * std::sqrt(210.0);
    const auto f_142 = 4.5 * std::sqrt(210.0);
    const auto f_143 = 1.5 * std::sqrt(210.0);
    const auto f_144 = 0.0703125 * std::sqrt(210.0);
    const auto f_145 = 0.1875 * std::sqrt(210.0);
    const auto f_146 = 0.5 * std::sqrt(210.0);
    const auto f_147 = 0.703125 * std::sqrt(42.0);
    const auto f_148 = 4.21875 * std::sqrt(42.0);
    const auto f_149 = 1.875 * std::sqrt(42.0);
    const auto f_150 = 11.25 * std::sqrt(42.0);
    const auto f_151 = 0.234375 * std::sqrt(42.0);
    const auto f_152 = 0.625 * std::sqrt(42.0);
    const auto f_153 = 2.4609375 * std::sqrt(6.0);
    const auto f_154 = 14.765625 * std::sqrt(6.0);
    const auto f_155 = 1.640625 * std::sqrt(6.0);
    const auto f_156 = 39.375 * std::sqrt(6.0);
    const auto f_157 = 0.8203125 * std::sqrt(6.0);
    const auto f_158 = 4.921875 * std::sqrt(6.0);
    const auto f_159 = 2.1875 * std::sqrt(6.0);
    const auto f_160 = 13.125 * std::sqrt(6.0);
    const auto f_161 = 1.09375 * std::sqrt(6.0);
    const auto f_162 = 17.5 * std::sqrt(6.0);
    const auto f_163 = 1.09375 * std::sqrt(3.0);
    const auto f_164 = 2.1875 * std::sqrt(3.0);
    const auto f_165 = 52.5 * std::sqrt(3.0);
    const auto f_166 = 17.5 * std::sqrt(3.0);
    const auto f_167 = 0.15625 * std::sqrt(42.0);
    const auto f_168 = 0.3125 * std::sqrt(42.0);
    const auto f_169 = 2.5 * std::sqrt(42.0);
    const auto f_170 = 15.0 * std::sqrt(42.0);
    const auto f_171 = 0.46875 * std::sqrt(21.0);
    const auto f_172 = 0.625 * std::sqrt(21.0);
    const auto f_173 = 0.9375 * std::sqrt(21.0);
    const auto f_174 = 1.25 * std::sqrt(21.0);
    const auto f_175 = 7.5 * std::sqrt(21.0);
    const auto f_176 = 10.0 * std::sqrt(21.0);
    const auto f_177 = 0.0234375 * std::sqrt(210.0);
    const auto f_178 = 0.046875 * std::sqrt(210.0);
    const auto f_179 = 0.0625 * std::sqrt(210.0);
    const auto f_180 = 0.09375 * std::sqrt(210.0);
    const auto f_181 = 0.125 * std::sqrt(210.0);
    const auto f_182 = 0.75 * std::sqrt(210.0);
    const auto f_183 = 3.0 * std::sqrt(210.0);
    const auto f_184 = std::sqrt(210.0);
    const auto f_185 = 0.078125 * std::sqrt(42.0);
    const auto f_186 = 0.2734375 * std::sqrt(6.0);
    const auto f_187 = 0.546875 * std::sqrt(6.0);
    const auto f_188 = 4.375 * std::sqrt(6.0);
    const auto f_189 = 2.1875 * std::sqrt(15.0);
    const auto f_190 = 4.375 * std::sqrt(15.0);
    const auto f_191 = 8.75 * std::sqrt(15.0);
    const auto f_192 = 3.5 * std::sqrt(15.0);
    const auto f_193 = 3.28125 * std::sqrt(30.0);
    const auto f_194 = 1.09375 * std::sqrt(30.0);
    const auto f_195 = 6.5625 * std::sqrt(30.0);
    const auto f_196 = 2.1875 * std::sqrt(30.0);
    const auto f_197 = 13.125 * std::sqrt(30.0);
    const auto f_198 = 4.375 * std::sqrt(30.0);
    const auto f_199 = 5.25 * std::sqrt(30.0);
    const auto f_200 = 1.75 * std::sqrt(30.0);
    const auto f_201 = 0.3125 * std::sqrt(105.0);
    const auto f_202 = 1.875 * std::sqrt(105.0);
    const auto f_203 = 0.625 * std::sqrt(105.0);
    const auto f_204 = 3.75 * std::sqrt(105.0);
    const auto f_205 = 1.25 * std::sqrt(105.0);
    const auto f_206 = 7.5 * std::sqrt(105.0);
    const auto f_207 = 0.5 * std::sqrt(105.0);
    const auto f_208 = 3.0 * std::sqrt(105.0);
    const auto f_209 = 0.46875 * std::sqrt(210.0);
    const auto f_210 = 0.625 * std::sqrt(210.0);
    const auto f_211 = 0.9375 * std::sqrt(210.0);
    const auto f_212 = 1.25 * std::sqrt(210.0);
    const auto f_213 = 1.875 * std::sqrt(210.0);
    const auto f_214 = 2.5 * std::sqrt(210.0);
    const auto f_215 = 0.234375 * std::sqrt(21.0);
    const auto f_216 = 2.5 * std::sqrt(21.0);
    const auto f_217 = 0.375 * std::sqrt(21.0);
    const auto f_218 = 0.75 * std::sqrt(21.0);
    const auto f_219 = 3.0 * std::sqrt(21.0);
    const auto f_220 = std::sqrt(21.0);
    const auto f_221 = 0.15625 * std::sqrt(105.0);
    const auto f_222 = 0.9375 * std::sqrt(105.0);
    const auto f_223 = 0.25 * std::sqrt(105.0);
    const auto f_224 = 1.5 * std::sqrt(105.0);
    const auto f_225 = 0.546875 * std::sqrt(15.0);
    const auto f_226 = 3.28125 * std::sqrt(15.0);
    const auto f_227 = 1.09375 * std::sqrt(15.0);
    const auto f_228 = 6.5625 * std::sqrt(15.0);
    const auto f_229 = 13.125 * std::sqrt(15.0);
    const auto f_230 = 0.875 * std::sqrt(15.0);
    const auto f_231 = 5.25 * std::sqrt(15.0);
    const auto f_232 = 0.15625 * std::sqrt(35.0);
    const auto f_233 = 0.46875 * std::sqrt(35.0);
    const auto f_234 = 2.8125 * std::sqrt(35.0);
    const auto f_235 = 5.625 * std::sqrt(35.0);
    const auto f_236 = 0.5 * std::sqrt(35.0);
    const auto f_237 = 0.234375 * std::sqrt(70.0);
    const auto f_238 = 0.078125 * std::sqrt(70.0);
    const auto f_239 = 0.703125 * std::sqrt(70.0);
    const auto f_240 = 4.21875 * std::sqrt(70.0);
    const auto f_241 = 1.40625 * std::sqrt(70.0);
    const auto f_242 = 8.4375 * std::sqrt(70.0);
    const auto f_243 = 2.8125 * std::sqrt(70.0);
    const auto f_244 = 1.875 * std::sqrt(70.0);
    const auto f_245 = 0.25 * std::sqrt(70.0);
    const auto f_246 = 0.15625 * std::sqrt(5.0);
    const auto f_247 = 0.9375 * std::sqrt(5.0);
    const auto f_248 = 0.46875 * std::sqrt(5.0);
    const auto f_249 = 2.8125 * std::sqrt(5.0);
    const auto f_250 = 16.875 * std::sqrt(5.0);
    const auto f_251 = 5.625 * std::sqrt(5.0);
    const auto f_252 = 33.75 * std::sqrt(5.0);
    const auto f_253 = 3.75 * std::sqrt(5.0);
    const auto f_254 = 22.5 * std::sqrt(5.0);
    const auto f_255 = 0.5 * std::sqrt(5.0);
    const auto f_256 = 3.0 * std::sqrt(5.0);
    const auto f_257 = 0.234375 * std::sqrt(10.0);
    const auto f_258 = 0.3125 * std::sqrt(10.0);
    const auto f_259 = 0.703125 * std::sqrt(10.0);
    const auto f_260 = 0.9375 * std::sqrt(10.0);
    const auto f_261 = 4.21875 * std::sqrt(10.0);
    const auto f_262 = 5.625 * std::sqrt(10.0);
    const auto f_263 = 8.4375 * std::sqrt(10.0);
    const auto f_264 = 11.25 * std::sqrt(10.0);
    const auto f_265 = 7.5 * std::sqrt(10.0);
    const auto f_266 = 0.75 * std::sqrt(10.0);
    const auto f_267 = std::sqrt(10.0);
    const auto f_268 = 0.078125 * std::sqrt(5.0);
    const auto f_269 = 0.234375 * std::sqrt(5.0);
    const auto f_270 = 1.40625 * std::sqrt(5.0);
    const auto f_271 = 8.4375 * std::sqrt(5.0);
    const auto f_272 = 1.875 * std::sqrt(5.0);
    const auto f_273 = 11.25 * std::sqrt(5.0);
    const auto f_274 = 0.25 * std::sqrt(5.0);
    const auto f_275 = 1.5 * std::sqrt(5.0);
    const auto f_276 = 0.0390625 * std::sqrt(35.0);
    const auto f_277 = 0.234375 * std::sqrt(35.0);
    const auto f_278 = 0.1171875 * std::sqrt(35.0);
    const auto f_279 = 0.703125 * std::sqrt(35.0);
    const auto f_280 = 4.21875 * std::sqrt(35.0);
    const auto f_281 = 1.40625 * std::sqrt(35.0);
    const auto f_282 = 8.4375 * std::sqrt(35.0);
    const auto f_283 = 0.9375 * std::sqrt(35.0);
    const auto f_284 = 0.125 * std::sqrt(35.0);
    const auto f_285 = 0.75 * std::sqrt(35.0);
    const auto f_286 = 1.640625 * std::sqrt(3.0);
    const auto f_287 = 0.546875 * std::sqrt(3.0);
    const auto f_288 = 0.3125 * std::sqrt(21.0);
    const auto f_289 = 0.01171875 * std::sqrt(210.0);
    const auto f_290 = 0.03125 * std::sqrt(210.0);
    const auto f_291 = 0.0390625 * std::sqrt(42.0);
    const auto f_292 = 0.13671875 * std::sqrt(6.0);
    const auto f_293 = 3.28125 * std::sqrt(5.0);
    const auto f_294 = 0.984375 * std::sqrt(10.0);
    const auto f_295 = 0.328125 * std::sqrt(10.0);
    const auto f_296 = 4.921875 * std::sqrt(10.0);
    const auto f_297 = 1.640625 * std::sqrt(10.0);
    const auto f_298 = 9.84375 * std::sqrt(10.0);
    const auto f_299 = 3.28125 * std::sqrt(10.0);
    const auto f_300 = 59.0625 * std::sqrt(10.0);
    const auto f_301 = 19.6875 * std::sqrt(10.0);
    const auto f_302 = 0.09375 * std::sqrt(35.0);
    const auto f_303 = 0.5625 * std::sqrt(35.0);
    const auto f_304 = 33.75 * std::sqrt(35.0);
    const auto f_305 = 0.140625 * std::sqrt(70.0);
    const auto f_306 = 0.1875 * std::sqrt(70.0);
    const auto f_307 = 0.9375 * std::sqrt(70.0);
    const auto f_308 = 11.25 * std::sqrt(70.0);
    const auto f_309 = 0.0703125 * std::sqrt(7.0);
    const auto f_310 = 0.140625 * std::sqrt(7.0);
    const auto f_311 = 0.1875 * std::sqrt(7.0);
    const auto f_312 = 0.3515625 * std::sqrt(7.0);
    const auto f_313 = 0.703125 * std::sqrt(7.0);
    const auto f_314 = 0.9375 * std::sqrt(7.0);
    const auto f_315 = 1.40625 * std::sqrt(7.0);
    const auto f_316 = 1.875 * std::sqrt(7.0);
    const auto f_317 = 4.21875 * std::sqrt(7.0);
    const auto f_318 = 8.4375 * std::sqrt(7.0);
    const auto f_319 = 33.75 * std::sqrt(7.0);
    const auto f_320 = 11.25 * std::sqrt(7.0);
    const auto f_321 = 0.046875 * std::sqrt(35.0);
    const auto f_322 = 0.28125 * std::sqrt(35.0);
    const auto f_323 = 16.875 * std::sqrt(35.0);
    const auto f_324 = 0.1640625 * std::sqrt(5.0);
    const auto f_325 = 0.984375 * std::sqrt(5.0);
    const auto f_326 = 0.8203125 * std::sqrt(5.0);
    const auto f_327 = 4.921875 * std::sqrt(5.0);
    const auto f_328 = 1.640625 * std::sqrt(5.0);
    const auto f_329 = 9.84375 * std::sqrt(5.0);
    const auto f_330 = 59.0625 * std::sqrt(5.0);
    const auto f_331 = 0.109375 * std::sqrt(330.0);
    const auto f_332 = 1.640625 * std::sqrt(330.0);
    const auto f_333 = 0.328125 * std::sqrt(165.0);
    const auto f_334 = 0.109375 * std::sqrt(165.0);
    const auto f_335 = 4.921875 * std::sqrt(165.0);
    const auto f_336 = 1.640625 * std::sqrt(165.0);
    const auto f_337 = 0.015625 * std::sqrt(2310.0);
    const auto f_338 = 0.234375 * std::sqrt(2310.0);
    const auto f_339 = 1.40625 * std::sqrt(2310.0);
    const auto f_340 = 0.046875 * std::sqrt(1155.0);
    const auto f_341 = 0.0625 * std::sqrt(1155.0);
    const auto f_342 = 0.703125 * std::sqrt(1155.0);
    const auto f_343 = 0.01171875 * std::sqrt(462.0);
    const auto f_344 = 0.0234375 * std::sqrt(462.0);
    const auto f_345 = 0.09375 * std::sqrt(462.0);
    const auto f_346 = 0.03125 * std::sqrt(462.0);
    const auto f_347 = 0.17578125 * std::sqrt(462.0);
    const auto f_348 = 0.3515625 * std::sqrt(462.0);
    const auto f_349 = 1.40625 * std::sqrt(462.0);
    const auto f_350 = 0.0078125 * std::sqrt(2310.0);
    const auto f_351 = 0.1171875 * std::sqrt(2310.0);
    const auto f_352 = 0.703125 * std::sqrt(2310.0);
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

    const auto *ig_0 = buffer.data(ig + 0);
    const auto *ig_1 = buffer.data(ig + 1);
    const auto *ig_2 = buffer.data(ig + 2);
    const auto *ig_3 = buffer.data(ig + 3);
    const auto *ig_4 = buffer.data(ig + 4);
    const auto *ig_5 = buffer.data(ig + 5);
    const auto *ig_6 = buffer.data(ig + 6);
    const auto *ig_7 = buffer.data(ig + 7);
    const auto *ig_8 = buffer.data(ig + 8);
    const auto *ig_9 = buffer.data(ig + 9);
    const auto *ig_10 = buffer.data(ig + 10);
    const auto *ig_11 = buffer.data(ig + 11);
    const auto *ig_12 = buffer.data(ig + 12);
    const auto *ig_13 = buffer.data(ig + 13);
    const auto *ig_14 = buffer.data(ig + 14);
    const auto *ig_15 = buffer.data(ig + 15);
    const auto *ig_16 = buffer.data(ig + 16);
    const auto *ig_17 = buffer.data(ig + 17);
    const auto *ig_18 = buffer.data(ig + 18);
    const auto *ig_19 = buffer.data(ig + 19);
    const auto *ig_20 = buffer.data(ig + 20);
    const auto *ig_21 = buffer.data(ig + 21);
    const auto *ig_22 = buffer.data(ig + 22);
    const auto *ig_23 = buffer.data(ig + 23);
    const auto *ig_24 = buffer.data(ig + 24);
    const auto *ig_25 = buffer.data(ig + 25);
    const auto *ig_26 = buffer.data(ig + 26);
    const auto *ig_27 = buffer.data(ig + 27);
    const auto *ig_28 = buffer.data(ig + 28);
    const auto *ig_29 = buffer.data(ig + 29);
    const auto *ig_30 = buffer.data(ig + 30);
    const auto *ig_31 = buffer.data(ig + 31);
    const auto *ig_32 = buffer.data(ig + 32);
    const auto *ig_33 = buffer.data(ig + 33);
    const auto *ig_34 = buffer.data(ig + 34);
    const auto *ig_35 = buffer.data(ig + 35);
    const auto *ig_36 = buffer.data(ig + 36);
    const auto *ig_37 = buffer.data(ig + 37);
    const auto *ig_38 = buffer.data(ig + 38);
    const auto *ig_39 = buffer.data(ig + 39);
    const auto *ig_40 = buffer.data(ig + 40);
    const auto *ig_41 = buffer.data(ig + 41);
    const auto *ig_42 = buffer.data(ig + 42);
    const auto *ig_43 = buffer.data(ig + 43);
    const auto *ig_44 = buffer.data(ig + 44);
    const auto *ig_45 = buffer.data(ig + 45);
    const auto *ig_46 = buffer.data(ig + 46);
    const auto *ig_47 = buffer.data(ig + 47);
    const auto *ig_48 = buffer.data(ig + 48);
    const auto *ig_49 = buffer.data(ig + 49);
    const auto *ig_50 = buffer.data(ig + 50);
    const auto *ig_51 = buffer.data(ig + 51);
    const auto *ig_52 = buffer.data(ig + 52);
    const auto *ig_53 = buffer.data(ig + 53);
    const auto *ig_54 = buffer.data(ig + 54);
    const auto *ig_55 = buffer.data(ig + 55);
    const auto *ig_56 = buffer.data(ig + 56);
    const auto *ig_57 = buffer.data(ig + 57);
    const auto *ig_58 = buffer.data(ig + 58);
    const auto *ig_59 = buffer.data(ig + 59);
    const auto *ig_60 = buffer.data(ig + 60);
    const auto *ig_61 = buffer.data(ig + 61);
    const auto *ig_62 = buffer.data(ig + 62);
    const auto *ig_63 = buffer.data(ig + 63);
    const auto *ig_64 = buffer.data(ig + 64);
    const auto *ig_65 = buffer.data(ig + 65);
    const auto *ig_66 = buffer.data(ig + 66);
    const auto *ig_67 = buffer.data(ig + 67);
    const auto *ig_68 = buffer.data(ig + 68);
    const auto *ig_69 = buffer.data(ig + 69);
    const auto *ig_70 = buffer.data(ig + 70);
    const auto *ig_71 = buffer.data(ig + 71);
    const auto *ig_72 = buffer.data(ig + 72);
    const auto *ig_73 = buffer.data(ig + 73);
    const auto *ig_74 = buffer.data(ig + 74);
    const auto *ig_75 = buffer.data(ig + 75);
    const auto *ig_76 = buffer.data(ig + 76);
    const auto *ig_77 = buffer.data(ig + 77);
    const auto *ig_78 = buffer.data(ig + 78);
    const auto *ig_79 = buffer.data(ig + 79);
    const auto *ig_80 = buffer.data(ig + 80);
    const auto *ig_81 = buffer.data(ig + 81);
    const auto *ig_82 = buffer.data(ig + 82);
    const auto *ig_83 = buffer.data(ig + 83);
    const auto *ig_84 = buffer.data(ig + 84);
    const auto *ig_85 = buffer.data(ig + 85);
    const auto *ig_86 = buffer.data(ig + 86);
    const auto *ig_87 = buffer.data(ig + 87);
    const auto *ig_88 = buffer.data(ig + 88);
    const auto *ig_89 = buffer.data(ig + 89);
    const auto *ig_90 = buffer.data(ig + 90);
    const auto *ig_91 = buffer.data(ig + 91);
    const auto *ig_92 = buffer.data(ig + 92);
    const auto *ig_93 = buffer.data(ig + 93);
    const auto *ig_94 = buffer.data(ig + 94);
    const auto *ig_95 = buffer.data(ig + 95);
    const auto *ig_96 = buffer.data(ig + 96);
    const auto *ig_97 = buffer.data(ig + 97);
    const auto *ig_98 = buffer.data(ig + 98);
    const auto *ig_99 = buffer.data(ig + 99);
    const auto *ig_100 = buffer.data(ig + 100);
    const auto *ig_101 = buffer.data(ig + 101);
    const auto *ig_102 = buffer.data(ig + 102);
    const auto *ig_103 = buffer.data(ig + 103);
    const auto *ig_104 = buffer.data(ig + 104);
    const auto *ig_105 = buffer.data(ig + 105);
    const auto *ig_106 = buffer.data(ig + 106);
    const auto *ig_107 = buffer.data(ig + 107);
    const auto *ig_108 = buffer.data(ig + 108);
    const auto *ig_109 = buffer.data(ig + 109);
    const auto *ig_110 = buffer.data(ig + 110);
    const auto *ig_111 = buffer.data(ig + 111);
    const auto *ig_112 = buffer.data(ig + 112);
    const auto *ig_113 = buffer.data(ig + 113);
    const auto *ig_114 = buffer.data(ig + 114);
    const auto *ig_115 = buffer.data(ig + 115);
    const auto *ig_116 = buffer.data(ig + 116);
    const auto *ig_117 = buffer.data(ig + 117);
    const auto *ig_118 = buffer.data(ig + 118);
    const auto *ig_119 = buffer.data(ig + 119);
    const auto *ig_120 = buffer.data(ig + 120);
    const auto *ig_121 = buffer.data(ig + 121);
    const auto *ig_122 = buffer.data(ig + 122);
    const auto *ig_123 = buffer.data(ig + 123);
    const auto *ig_124 = buffer.data(ig + 124);
    const auto *ig_125 = buffer.data(ig + 125);
    const auto *ig_126 = buffer.data(ig + 126);
    const auto *ig_127 = buffer.data(ig + 127);
    const auto *ig_128 = buffer.data(ig + 128);
    const auto *ig_129 = buffer.data(ig + 129);
    const auto *ig_130 = buffer.data(ig + 130);
    const auto *ig_131 = buffer.data(ig + 131);
    const auto *ig_132 = buffer.data(ig + 132);
    const auto *ig_133 = buffer.data(ig + 133);
    const auto *ig_134 = buffer.data(ig + 134);
    const auto *ig_135 = buffer.data(ig + 135);
    const auto *ig_136 = buffer.data(ig + 136);
    const auto *ig_137 = buffer.data(ig + 137);
    const auto *ig_138 = buffer.data(ig + 138);
    const auto *ig_139 = buffer.data(ig + 139);
    const auto *ig_140 = buffer.data(ig + 140);
    const auto *ig_141 = buffer.data(ig + 141);
    const auto *ig_142 = buffer.data(ig + 142);
    const auto *ig_143 = buffer.data(ig + 143);
    const auto *ig_144 = buffer.data(ig + 144);
    const auto *ig_145 = buffer.data(ig + 145);
    const auto *ig_146 = buffer.data(ig + 146);
    const auto *ig_147 = buffer.data(ig + 147);
    const auto *ig_148 = buffer.data(ig + 148);
    const auto *ig_149 = buffer.data(ig + 149);
    const auto *ig_150 = buffer.data(ig + 150);
    const auto *ig_151 = buffer.data(ig + 151);
    const auto *ig_152 = buffer.data(ig + 152);
    const auto *ig_153 = buffer.data(ig + 153);
    const auto *ig_154 = buffer.data(ig + 154);
    const auto *ig_155 = buffer.data(ig + 155);
    const auto *ig_156 = buffer.data(ig + 156);
    const auto *ig_157 = buffer.data(ig + 157);
    const auto *ig_158 = buffer.data(ig + 158);
    const auto *ig_159 = buffer.data(ig + 159);
    const auto *ig_160 = buffer.data(ig + 160);
    const auto *ig_161 = buffer.data(ig + 161);
    const auto *ig_162 = buffer.data(ig + 162);
    const auto *ig_163 = buffer.data(ig + 163);
    const auto *ig_164 = buffer.data(ig + 164);
    const auto *ig_165 = buffer.data(ig + 165);
    const auto *ig_166 = buffer.data(ig + 166);
    const auto *ig_167 = buffer.data(ig + 167);
    const auto *ig_168 = buffer.data(ig + 168);
    const auto *ig_169 = buffer.data(ig + 169);
    const auto *ig_170 = buffer.data(ig + 170);
    const auto *ig_171 = buffer.data(ig + 171);
    const auto *ig_172 = buffer.data(ig + 172);
    const auto *ig_173 = buffer.data(ig + 173);
    const auto *ig_174 = buffer.data(ig + 174);
    const auto *ig_175 = buffer.data(ig + 175);
    const auto *ig_176 = buffer.data(ig + 176);
    const auto *ig_177 = buffer.data(ig + 177);
    const auto *ig_178 = buffer.data(ig + 178);
    const auto *ig_179 = buffer.data(ig + 179);
    const auto *ig_180 = buffer.data(ig + 180);
    const auto *ig_181 = buffer.data(ig + 181);
    const auto *ig_182 = buffer.data(ig + 182);
    const auto *ig_183 = buffer.data(ig + 183);
    const auto *ig_184 = buffer.data(ig + 184);
    const auto *ig_185 = buffer.data(ig + 185);
    const auto *ig_186 = buffer.data(ig + 186);
    const auto *ig_187 = buffer.data(ig + 187);
    const auto *ig_188 = buffer.data(ig + 188);
    const auto *ig_189 = buffer.data(ig + 189);
    const auto *ig_190 = buffer.data(ig + 190);
    const auto *ig_191 = buffer.data(ig + 191);
    const auto *ig_192 = buffer.data(ig + 192);
    const auto *ig_193 = buffer.data(ig + 193);
    const auto *ig_194 = buffer.data(ig + 194);
    const auto *ig_195 = buffer.data(ig + 195);
    const auto *ig_196 = buffer.data(ig + 196);
    const auto *ig_197 = buffer.data(ig + 197);
    const auto *ig_198 = buffer.data(ig + 198);
    const auto *ig_199 = buffer.data(ig + 199);
    const auto *ig_200 = buffer.data(ig + 200);
    const auto *ig_201 = buffer.data(ig + 201);
    const auto *ig_202 = buffer.data(ig + 202);
    const auto *ig_203 = buffer.data(ig + 203);
    const auto *ig_204 = buffer.data(ig + 204);
    const auto *ig_205 = buffer.data(ig + 205);
    const auto *ig_206 = buffer.data(ig + 206);
    const auto *ig_207 = buffer.data(ig + 207);
    const auto *ig_208 = buffer.data(ig + 208);
    const auto *ig_209 = buffer.data(ig + 209);
    const auto *ig_210 = buffer.data(ig + 210);
    const auto *ig_211 = buffer.data(ig + 211);
    const auto *ig_212 = buffer.data(ig + 212);
    const auto *ig_213 = buffer.data(ig + 213);
    const auto *ig_214 = buffer.data(ig + 214);
    const auto *ig_215 = buffer.data(ig + 215);
    const auto *ig_216 = buffer.data(ig + 216);
    const auto *ig_217 = buffer.data(ig + 217);
    const auto *ig_218 = buffer.data(ig + 218);
    const auto *ig_219 = buffer.data(ig + 219);
    const auto *ig_220 = buffer.data(ig + 220);
    const auto *ig_221 = buffer.data(ig + 221);
    const auto *ig_222 = buffer.data(ig + 222);
    const auto *ig_223 = buffer.data(ig + 223);
    const auto *ig_224 = buffer.data(ig + 224);
    const auto *ig_225 = buffer.data(ig + 225);
    const auto *ig_226 = buffer.data(ig + 226);
    const auto *ig_227 = buffer.data(ig + 227);
    const auto *ig_228 = buffer.data(ig + 228);
    const auto *ig_229 = buffer.data(ig + 229);
    const auto *ig_230 = buffer.data(ig + 230);
    const auto *ig_231 = buffer.data(ig + 231);
    const auto *ig_232 = buffer.data(ig + 232);
    const auto *ig_233 = buffer.data(ig + 233);
    const auto *ig_234 = buffer.data(ig + 234);
    const auto *ig_235 = buffer.data(ig + 235);
    const auto *ig_236 = buffer.data(ig + 236);
    const auto *ig_237 = buffer.data(ig + 237);
    const auto *ig_238 = buffer.data(ig + 238);
    const auto *ig_239 = buffer.data(ig + 239);
    const auto *ig_240 = buffer.data(ig + 240);
    const auto *ig_241 = buffer.data(ig + 241);
    const auto *ig_242 = buffer.data(ig + 242);
    const auto *ig_243 = buffer.data(ig + 243);
    const auto *ig_244 = buffer.data(ig + 244);
    const auto *ig_245 = buffer.data(ig + 245);
    const auto *ig_246 = buffer.data(ig + 246);
    const auto *ig_247 = buffer.data(ig + 247);
    const auto *ig_248 = buffer.data(ig + 248);
    const auto *ig_249 = buffer.data(ig + 249);
    const auto *ig_250 = buffer.data(ig + 250);
    const auto *ig_251 = buffer.data(ig + 251);
    const auto *ig_252 = buffer.data(ig + 252);
    const auto *ig_253 = buffer.data(ig + 253);
    const auto *ig_254 = buffer.data(ig + 254);
    const auto *ig_255 = buffer.data(ig + 255);
    const auto *ig_256 = buffer.data(ig + 256);
    const auto *ig_257 = buffer.data(ig + 257);
    const auto *ig_258 = buffer.data(ig + 258);
    const auto *ig_259 = buffer.data(ig + 259);
    const auto *ig_260 = buffer.data(ig + 260);
    const auto *ig_261 = buffer.data(ig + 261);
    const auto *ig_262 = buffer.data(ig + 262);
    const auto *ig_263 = buffer.data(ig + 263);
    const auto *ig_264 = buffer.data(ig + 264);
    const auto *ig_265 = buffer.data(ig + 265);
    const auto *ig_266 = buffer.data(ig + 266);
    const auto *ig_267 = buffer.data(ig + 267);
    const auto *ig_268 = buffer.data(ig + 268);
    const auto *ig_269 = buffer.data(ig + 269);
    const auto *ig_270 = buffer.data(ig + 270);
    const auto *ig_271 = buffer.data(ig + 271);
    const auto *ig_272 = buffer.data(ig + 272);
    const auto *ig_273 = buffer.data(ig + 273);
    const auto *ig_274 = buffer.data(ig + 274);
    const auto *ig_275 = buffer.data(ig + 275);
    const auto *ig_276 = buffer.data(ig + 276);
    const auto *ig_277 = buffer.data(ig + 277);
    const auto *ig_278 = buffer.data(ig + 278);
    const auto *ig_279 = buffer.data(ig + 279);
    const auto *ig_280 = buffer.data(ig + 280);
    const auto *ig_281 = buffer.data(ig + 281);
    const auto *ig_282 = buffer.data(ig + 282);
    const auto *ig_283 = buffer.data(ig + 283);
    const auto *ig_284 = buffer.data(ig + 284);
    const auto *ig_285 = buffer.data(ig + 285);
    const auto *ig_286 = buffer.data(ig + 286);
    const auto *ig_287 = buffer.data(ig + 287);
    const auto *ig_288 = buffer.data(ig + 288);
    const auto *ig_289 = buffer.data(ig + 289);
    const auto *ig_290 = buffer.data(ig + 290);
    const auto *ig_291 = buffer.data(ig + 291);
    const auto *ig_292 = buffer.data(ig + 292);
    const auto *ig_293 = buffer.data(ig + 293);
    const auto *ig_294 = buffer.data(ig + 294);
    const auto *ig_295 = buffer.data(ig + 295);
    const auto *ig_296 = buffer.data(ig + 296);
    const auto *ig_297 = buffer.data(ig + 297);
    const auto *ig_298 = buffer.data(ig + 298);
    const auto *ig_299 = buffer.data(ig + 299);
    const auto *ig_300 = buffer.data(ig + 300);
    const auto *ig_301 = buffer.data(ig + 301);
    const auto *ig_302 = buffer.data(ig + 302);
    const auto *ig_303 = buffer.data(ig + 303);
    const auto *ig_304 = buffer.data(ig + 304);
    const auto *ig_305 = buffer.data(ig + 305);
    const auto *ig_306 = buffer.data(ig + 306);
    const auto *ig_307 = buffer.data(ig + 307);
    const auto *ig_308 = buffer.data(ig + 308);
    const auto *ig_309 = buffer.data(ig + 309);
    const auto *ig_310 = buffer.data(ig + 310);
    const auto *ig_311 = buffer.data(ig + 311);
    const auto *ig_312 = buffer.data(ig + 312);
    const auto *ig_313 = buffer.data(ig + 313);
    const auto *ig_314 = buffer.data(ig + 314);
    const auto *ig_315 = buffer.data(ig + 315);
    const auto *ig_316 = buffer.data(ig + 316);
    const auto *ig_317 = buffer.data(ig + 317);
    const auto *ig_318 = buffer.data(ig + 318);
    const auto *ig_319 = buffer.data(ig + 319);
    const auto *ig_320 = buffer.data(ig + 320);
    const auto *ig_321 = buffer.data(ig + 321);
    const auto *ig_322 = buffer.data(ig + 322);
    const auto *ig_323 = buffer.data(ig + 323);
    const auto *ig_324 = buffer.data(ig + 324);
    const auto *ig_325 = buffer.data(ig + 325);
    const auto *ig_326 = buffer.data(ig + 326);
    const auto *ig_327 = buffer.data(ig + 327);
    const auto *ig_328 = buffer.data(ig + 328);
    const auto *ig_329 = buffer.data(ig + 329);
    const auto *ig_330 = buffer.data(ig + 330);
    const auto *ig_331 = buffer.data(ig + 331);
    const auto *ig_332 = buffer.data(ig + 332);
    const auto *ig_333 = buffer.data(ig + 333);
    const auto *ig_334 = buffer.data(ig + 334);
    const auto *ig_335 = buffer.data(ig + 335);
    const auto *ig_336 = buffer.data(ig + 336);
    const auto *ig_337 = buffer.data(ig + 337);
    const auto *ig_338 = buffer.data(ig + 338);
    const auto *ig_339 = buffer.data(ig + 339);
    const auto *ig_340 = buffer.data(ig + 340);
    const auto *ig_341 = buffer.data(ig + 341);
    const auto *ig_342 = buffer.data(ig + 342);
    const auto *ig_343 = buffer.data(ig + 343);
    const auto *ig_344 = buffer.data(ig + 344);
    const auto *ig_345 = buffer.data(ig + 345);
    const auto *ig_346 = buffer.data(ig + 346);
    const auto *ig_347 = buffer.data(ig + 347);
    const auto *ig_348 = buffer.data(ig + 348);
    const auto *ig_349 = buffer.data(ig + 349);
    const auto *ig_350 = buffer.data(ig + 350);
    const auto *ig_351 = buffer.data(ig + 351);
    const auto *ig_352 = buffer.data(ig + 352);
    const auto *ig_353 = buffer.data(ig + 353);
    const auto *ig_354 = buffer.data(ig + 354);
    const auto *ig_355 = buffer.data(ig + 355);
    const auto *ig_356 = buffer.data(ig + 356);
    const auto *ig_357 = buffer.data(ig + 357);
    const auto *ig_358 = buffer.data(ig + 358);
    const auto *ig_359 = buffer.data(ig + 359);
    const auto *ig_360 = buffer.data(ig + 360);
    const auto *ig_361 = buffer.data(ig + 361);
    const auto *ig_362 = buffer.data(ig + 362);
    const auto *ig_363 = buffer.data(ig + 363);
    const auto *ig_364 = buffer.data(ig + 364);
    const auto *ig_365 = buffer.data(ig + 365);
    const auto *ig_366 = buffer.data(ig + 366);
    const auto *ig_367 = buffer.data(ig + 367);
    const auto *ig_368 = buffer.data(ig + 368);
    const auto *ig_369 = buffer.data(ig + 369);
    const auto *ig_370 = buffer.data(ig + 370);
    const auto *ig_371 = buffer.data(ig + 371);
    const auto *ig_372 = buffer.data(ig + 372);
    const auto *ig_373 = buffer.data(ig + 373);
    const auto *ig_374 = buffer.data(ig + 374);
    const auto *ig_375 = buffer.data(ig + 375);
    const auto *ig_376 = buffer.data(ig + 376);
    const auto *ig_377 = buffer.data(ig + 377);
    const auto *ig_378 = buffer.data(ig + 378);
    const auto *ig_379 = buffer.data(ig + 379);
    const auto *ig_380 = buffer.data(ig + 380);
    const auto *ig_381 = buffer.data(ig + 381);
    const auto *ig_382 = buffer.data(ig + 382);
    const auto *ig_383 = buffer.data(ig + 383);
    const auto *ig_384 = buffer.data(ig + 384);
    const auto *ig_385 = buffer.data(ig + 385);
    const auto *ig_386 = buffer.data(ig + 386);
    const auto *ig_387 = buffer.data(ig + 387);
    const auto *ig_388 = buffer.data(ig + 388);
    const auto *ig_389 = buffer.data(ig + 389);
    const auto *ig_390 = buffer.data(ig + 390);
    const auto *ig_391 = buffer.data(ig + 391);
    const auto *ig_392 = buffer.data(ig + 392);
    const auto *ig_393 = buffer.data(ig + 393);
    const auto *ig_394 = buffer.data(ig + 394);
    const auto *ig_395 = buffer.data(ig + 395);
    const auto *ig_396 = buffer.data(ig + 396);
    const auto *ig_397 = buffer.data(ig + 397);
    const auto *ig_398 = buffer.data(ig + 398);
    const auto *ig_399 = buffer.data(ig + 399);
    const auto *ig_400 = buffer.data(ig + 400);
    const auto *ig_401 = buffer.data(ig + 401);
    const auto *ig_402 = buffer.data(ig + 402);
    const auto *ig_403 = buffer.data(ig + 403);
    const auto *ig_404 = buffer.data(ig + 404);
    const auto *ig_405 = buffer.data(ig + 405);
    const auto *ig_406 = buffer.data(ig + 406);
    const auto *ig_407 = buffer.data(ig + 407);
    const auto *ig_408 = buffer.data(ig + 408);
    const auto *ig_409 = buffer.data(ig + 409);
    const auto *ig_410 = buffer.data(ig + 410);
    const auto *ig_411 = buffer.data(ig + 411);
    const auto *ig_412 = buffer.data(ig + 412);
    const auto *ig_413 = buffer.data(ig + 413);
    const auto *ig_414 = buffer.data(ig + 414);
    const auto *ig_415 = buffer.data(ig + 415);
    const auto *ig_416 = buffer.data(ig + 416);
    const auto *ig_417 = buffer.data(ig + 417);
    const auto *ig_418 = buffer.data(ig + 418);
    const auto *ig_419 = buffer.data(ig + 419);

#pragma omp simd aligned(ig_16, ig_19, ig_21, ig_26, ig_91, ig_94, ig_96, ig_101, ig_226, \
                         ig_229, ig_231, ig_236 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * ig_16[k]
                 - f_0 * ig_21[k]
                 - f_1 * ig_91[k]
                 + f_1 * ig_96[k]
                 + f_0 * ig_226[k]
                 - f_0 * ig_231[k];

        g_1[k] = f_2 * ig_19[k]
                 - f_3 * ig_26[k]
                 - f_4 * ig_94[k]
                 + f_5 * ig_101[k]
                 + f_2 * ig_229[k]
                 - f_3 * ig_236[k];
    }

#pragma omp simd aligned(ig_16, ig_21, ig_23, ig_91, ig_96, ig_98, ig_226, ig_231, \
                         ig_233 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_2[k] = -f_6 * ig_16[k]
                 - f_6 * ig_21[k]
                 + f_7 * ig_23[k]
                 + f_8 * ig_91[k]
                 + f_8 * ig_96[k]
                 - f_9 * ig_98[k]
                 - f_6 * ig_226[k]
                 - f_6 * ig_231[k]
                 + f_7 * ig_233[k];
    }

#pragma omp simd aligned(ig_19, ig_26, ig_28, ig_94, ig_101, ig_103, ig_229, ig_236, \
                         ig_238 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = -f_10 * ig_19[k]
                 - f_10 * ig_26[k]
                 + f_11 * ig_28[k]
                 + f_12 * ig_94[k]
                 + f_12 * ig_101[k]
                 - f_13 * ig_103[k]
                 - f_10 * ig_229[k]
                 - f_10 * ig_236[k]
                 + f_11 * ig_238[k];
    }

#pragma omp simd aligned(ig_15, ig_18, ig_20, ig_25, ig_27, ig_29, ig_90, ig_93, ig_95, \
                         ig_100, ig_102, ig_104, ig_225, ig_228, ig_230, ig_235, ig_237, \
                         ig_239 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = f_14 * ig_15[k]
                 + f_15 * ig_18[k]
                 - f_16 * ig_20[k]
                 + f_14 * ig_25[k]
                 - f_16 * ig_27[k]
                 + f_17 * ig_29[k]
                 - f_18 * ig_90[k]
                 - f_19 * ig_93[k]
                 + f_20 * ig_95[k]
                 - f_18 * ig_100[k]
                 + f_20 * ig_102[k]
                 - f_21 * ig_104[k]
                 + f_14 * ig_225[k]
                 + f_15 * ig_228[k]
                 - f_16 * ig_230[k]
                 + f_14 * ig_235[k]
                 - f_16 * ig_237[k]
                 + f_17 * ig_239[k];
    }

#pragma omp simd aligned(ig_17, ig_22, ig_24, ig_92, ig_97, ig_99, ig_227, ig_232, \
                         ig_234 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = -f_10 * ig_17[k]
                 - f_10 * ig_22[k]
                 + f_11 * ig_24[k]
                 + f_12 * ig_92[k]
                 + f_12 * ig_97[k]
                 - f_13 * ig_99[k]
                 - f_10 * ig_227[k]
                 - f_10 * ig_232[k]
                 + f_11 * ig_234[k];
    }

#pragma omp simd aligned(ig_15, ig_20, ig_25, ig_27, ig_90, ig_95, ig_100, ig_102, ig_225, \
                         ig_230, ig_235, ig_237 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = -f_22 * ig_15[k]
                 + f_23 * ig_20[k]
                 + f_22 * ig_25[k]
                 - f_23 * ig_27[k]
                 + f_24 * ig_90[k]
                 - f_25 * ig_95[k]
                 - f_24 * ig_100[k]
                 + f_25 * ig_102[k]
                 - f_22 * ig_225[k]
                 + f_23 * ig_230[k]
                 + f_22 * ig_235[k]
                 - f_23 * ig_237[k];
    }

#pragma omp simd aligned(ig_17, ig_22, ig_92, ig_97, ig_227, ig_232 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = f_3 * ig_17[k]
                 - f_2 * ig_22[k]
                 - f_5 * ig_92[k]
                 + f_4 * ig_97[k]
                 + f_3 * ig_227[k]
                 - f_2 * ig_232[k];
    }

#pragma omp simd aligned(ig_15, ig_18, ig_25, ig_90, ig_93, ig_100, ig_225, ig_228, \
                         ig_235 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_8[k] = f_26 * ig_15[k]
                 - f_27 * ig_18[k]
                 + f_26 * ig_25[k]
                 - f_28 * ig_90[k]
                 + f_29 * ig_93[k]
                 - f_28 * ig_100[k]
                 + f_26 * ig_225[k]
                 - f_27 * ig_228[k]
                 + f_26 * ig_235[k];
    }

#pragma omp simd aligned(ig_61, ig_64, ig_66, ig_71, ig_166, ig_169, ig_171, ig_176, ig_331, \
                         ig_334, ig_336, ig_341 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_9[k] = f_30 * ig_61[k]
                 - f_30 * ig_66[k]
                 - f_31 * ig_166[k]
                 + f_31 * ig_171[k]
                 + f_32 * ig_331[k]
                 - f_32 * ig_336[k];

        g_10[k] = f_33 * ig_64[k]
                  - f_34 * ig_71[k]
                  - f_35 * ig_169[k]
                  + f_36 * ig_176[k]
                  + f_37 * ig_334[k]
                  - f_38 * ig_341[k];
    }

#pragma omp simd aligned(ig_61, ig_66, ig_68, ig_166, ig_171, ig_173, ig_331, ig_336, \
                         ig_338 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_11[k] = -f_39 * ig_61[k]
                  - f_39 * ig_66[k]
                  + f_40 * ig_68[k]
                  + f_41 * ig_166[k]
                  + f_41 * ig_171[k]
                  - f_42 * ig_173[k]
                  - f_43 * ig_331[k]
                  - f_43 * ig_336[k]
                  + f_44 * ig_338[k];
    }

#pragma omp simd aligned(ig_64, ig_71, ig_73, ig_169, ig_176, ig_178, ig_334, ig_341, \
                         ig_343 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_12[k] = -f_45 * ig_64[k]
                  - f_45 * ig_71[k]
                  + f_46 * ig_73[k]
                  + f_47 * ig_169[k]
                  + f_47 * ig_176[k]
                  - f_48 * ig_178[k]
                  - f_49 * ig_334[k]
                  - f_49 * ig_341[k]
                  + f_50 * ig_343[k];
    }

#pragma omp simd aligned(ig_60, ig_63, ig_65, ig_70, ig_72, ig_74, ig_165, ig_168, ig_170, \
                         ig_175, ig_177, ig_179, ig_330, ig_333, ig_335, ig_340, ig_342, \
                         ig_344 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = f_51 * ig_60[k]
                  + f_52 * ig_63[k]
                  - f_53 * ig_65[k]
                  + f_51 * ig_70[k]
                  - f_53 * ig_72[k]
                  + f_54 * ig_74[k]
                  - f_52 * ig_165[k]
                  - f_55 * ig_168[k]
                  + f_56 * ig_170[k]
                  - f_52 * ig_175[k]
                  + f_56 * ig_177[k]
                  - f_57 * ig_179[k]
                  + f_58 * ig_330[k]
                  + f_59 * ig_333[k]
                  - f_60 * ig_335[k]
                  + f_58 * ig_340[k]
                  - f_60 * ig_342[k]
                  + f_61 * ig_344[k];
    }

#pragma omp simd aligned(ig_62, ig_67, ig_69, ig_167, ig_172, ig_174, ig_332, ig_337, \
                         ig_339 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = -f_45 * ig_62[k]
                  - f_45 * ig_67[k]
                  + f_46 * ig_69[k]
                  + f_47 * ig_167[k]
                  + f_47 * ig_172[k]
                  - f_48 * ig_174[k]
                  - f_49 * ig_332[k]
                  - f_49 * ig_337[k]
                  + f_50 * ig_339[k];
    }

#pragma omp simd aligned(ig_60, ig_65, ig_70, ig_72, ig_165, ig_170, ig_175, ig_177, ig_330, \
                         ig_335, ig_340, ig_342 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_15[k] = -f_62 * ig_60[k]
                  + f_63 * ig_65[k]
                  + f_62 * ig_70[k]
                  - f_63 * ig_72[k]
                  + f_39 * ig_165[k]
                  - f_40 * ig_170[k]
                  - f_39 * ig_175[k]
                  + f_40 * ig_177[k]
                  - f_64 * ig_330[k]
                  + f_65 * ig_335[k]
                  + f_64 * ig_340[k]
                  - f_65 * ig_342[k];
    }

#pragma omp simd aligned(ig_62, ig_67, ig_167, ig_172, ig_332, ig_337 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_16[k] = f_34 * ig_62[k]
                  - f_33 * ig_67[k]
                  - f_36 * ig_167[k]
                  + f_35 * ig_172[k]
                  + f_38 * ig_332[k]
                  - f_37 * ig_337[k];
    }

#pragma omp simd aligned(ig_60, ig_63, ig_70, ig_165, ig_168, ig_175, ig_330, ig_333, \
                         ig_340 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = f_66 * ig_60[k]
                  - f_67 * ig_63[k]
                  + f_66 * ig_70[k]
                  - f_68 * ig_165[k]
                  + f_69 * ig_168[k]
                  - f_68 * ig_175[k]
                  + f_70 * ig_330[k]
                  - f_71 * ig_333[k]
                  + f_70 * ig_340[k];
    }

#pragma omp simd aligned(ig_16, ig_21, ig_121, ig_126, ig_226, ig_231, ig_256, \
                         ig_261 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = -f_72 * ig_16[k]
                  + f_72 * ig_21[k]
                  + f_73 * ig_121[k]
                  - f_73 * ig_126[k]
                  + f_72 * ig_226[k]
                  - f_72 * ig_231[k]
                  - f_73 * ig_256[k]
                  + f_73 * ig_261[k];
    }

#pragma omp simd aligned(ig_19, ig_26, ig_124, ig_131, ig_229, ig_236, ig_259, \
                         ig_266 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = -f_74 * ig_19[k]
                  + f_75 * ig_26[k]
                  + f_76 * ig_124[k]
                  - f_77 * ig_131[k]
                  + f_74 * ig_229[k]
                  - f_75 * ig_236[k]
                  - f_76 * ig_259[k]
                  + f_77 * ig_266[k];
    }

#pragma omp simd aligned(ig_16, ig_21, ig_23, ig_121, ig_126, ig_128, ig_226, ig_231, ig_233, \
                         ig_256, ig_261, ig_263 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = f_78 * ig_16[k]
                  + f_78 * ig_21[k]
                  - f_79 * ig_23[k]
                  - f_80 * ig_121[k]
                  - f_80 * ig_126[k]
                  + f_81 * ig_128[k]
                  - f_78 * ig_226[k]
                  - f_78 * ig_231[k]
                  + f_79 * ig_233[k]
                  + f_80 * ig_256[k]
                  + f_80 * ig_261[k]
                  - f_81 * ig_263[k];
    }

#pragma omp simd aligned(ig_19, ig_26, ig_28, ig_124, ig_131, ig_133, ig_229, ig_236, ig_238, \
                         ig_259, ig_266, ig_268 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = f_82 * ig_19[k]
                  + f_82 * ig_26[k]
                  - f_83 * ig_28[k]
                  - f_84 * ig_124[k]
                  - f_84 * ig_131[k]
                  + f_85 * ig_133[k]
                  - f_82 * ig_229[k]
                  - f_82 * ig_236[k]
                  + f_83 * ig_238[k]
                  + f_84 * ig_259[k]
                  + f_84 * ig_266[k]
                  - f_85 * ig_268[k];
    }

#pragma omp simd aligned(ig_15, ig_18, ig_20, ig_25, ig_27, ig_29, ig_120, ig_123, ig_125, \
                         ig_130, ig_132, ig_134, ig_225, ig_228, ig_230, ig_235, ig_237, \
                         ig_239, ig_255, ig_258, ig_260, ig_265, ig_267, \
                         ig_269 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = -f_86 * ig_15[k]
                  - f_87 * ig_18[k]
                  + f_88 * ig_20[k]
                  - f_86 * ig_25[k]
                  + f_88 * ig_27[k]
                  - f_89 * ig_29[k]
                  + f_90 * ig_120[k]
                  + f_91 * ig_123[k]
                  - f_92 * ig_125[k]
                  + f_90 * ig_130[k]
                  - f_92 * ig_132[k]
                  + f_93 * ig_134[k]
                  + f_86 * ig_225[k]
                  + f_87 * ig_228[k]
                  - f_88 * ig_230[k]
                  + f_86 * ig_235[k]
                  - f_88 * ig_237[k]
                  + f_89 * ig_239[k]
                  - f_90 * ig_255[k]
                  - f_91 * ig_258[k]
                  + f_92 * ig_260[k]
                  - f_90 * ig_265[k]
                  + f_92 * ig_267[k]
                  - f_93 * ig_269[k];
    }

#pragma omp simd aligned(ig_17, ig_22, ig_24, ig_122, ig_127, ig_129, ig_227, ig_232, ig_234, \
                         ig_257, ig_262, ig_264 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = f_82 * ig_17[k]
                  + f_82 * ig_22[k]
                  - f_83 * ig_24[k]
                  - f_84 * ig_122[k]
                  - f_84 * ig_127[k]
                  + f_85 * ig_129[k]
                  - f_82 * ig_227[k]
                  - f_82 * ig_232[k]
                  + f_83 * ig_234[k]
                  + f_84 * ig_257[k]
                  + f_84 * ig_262[k]
                  - f_85 * ig_264[k];
    }

#pragma omp simd aligned(ig_15, ig_20, ig_25, ig_27, ig_120, ig_125, ig_130, ig_132, ig_225, \
                         ig_230, ig_235, ig_237, ig_255, ig_260, ig_265, \
                         ig_267 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = f_94 * ig_15[k]
                  - f_95 * ig_20[k]
                  - f_94 * ig_25[k]
                  + f_95 * ig_27[k]
                  - f_96 * ig_120[k]
                  + f_97 * ig_125[k]
                  + f_96 * ig_130[k]
                  - f_97 * ig_132[k]
                  - f_94 * ig_225[k]
                  + f_95 * ig_230[k]
                  + f_94 * ig_235[k]
                  - f_95 * ig_237[k]
                  + f_96 * ig_255[k]
                  - f_97 * ig_260[k]
                  - f_96 * ig_265[k]
                  + f_97 * ig_267[k];
    }

#pragma omp simd aligned(ig_17, ig_22, ig_122, ig_127, ig_227, ig_232, ig_257, \
                         ig_262 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = -f_75 * ig_17[k]
                  + f_74 * ig_22[k]
                  + f_77 * ig_122[k]
                  - f_76 * ig_127[k]
                  + f_75 * ig_227[k]
                  - f_74 * ig_232[k]
                  - f_77 * ig_257[k]
                  + f_76 * ig_262[k];
    }

#pragma omp simd aligned(ig_15, ig_18, ig_25, ig_120, ig_123, ig_130, ig_225, ig_228, ig_235, \
                         ig_255, ig_258, ig_265 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = -f_98 * ig_15[k]
                  + f_99 * ig_18[k]
                  - f_98 * ig_25[k]
                  + f_100 * ig_120[k]
                  - f_101 * ig_123[k]
                  + f_100 * ig_130[k]
                  + f_98 * ig_225[k]
                  - f_99 * ig_228[k]
                  + f_98 * ig_235[k]
                  - f_100 * ig_255[k]
                  + f_101 * ig_258[k]
                  - f_100 * ig_265[k];
    }

#pragma omp simd aligned(ig_61, ig_66, ig_166, ig_171, ig_196, ig_201, ig_331, ig_336, ig_361, \
                         ig_366 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = -f_102 * ig_61[k]
                  + f_102 * ig_66[k]
                  - f_103 * ig_166[k]
                  + f_103 * ig_171[k]
                  + f_104 * ig_196[k]
                  - f_104 * ig_201[k]
                  + f_105 * ig_331[k]
                  - f_105 * ig_336[k]
                  - f_106 * ig_361[k]
                  + f_106 * ig_366[k];
    }

#pragma omp simd aligned(ig_64, ig_71, ig_169, ig_176, ig_199, ig_206, ig_334, ig_341, ig_364, \
                         ig_371 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = -f_107 * ig_64[k]
                  + f_108 * ig_71[k]
                  - f_109 * ig_169[k]
                  + f_110 * ig_176[k]
                  + f_111 * ig_199[k]
                  - f_112 * ig_206[k]
                  + f_108 * ig_334[k]
                  - f_113 * ig_341[k]
                  - f_112 * ig_364[k]
                  + f_114 * ig_371[k];
    }

#pragma omp simd aligned(ig_61, ig_66, ig_68, ig_166, ig_171, ig_173, ig_196, ig_201, ig_203, \
                         ig_331, ig_336, ig_338, ig_361, ig_366, \
                         ig_368 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = f_115 * ig_61[k]
                  + f_115 * ig_66[k]
                  - f_116 * ig_68[k]
                  + f_117 * ig_166[k]
                  + f_117 * ig_171[k]
                  - f_118 * ig_173[k]
                  - f_119 * ig_196[k]
                  - f_119 * ig_201[k]
                  + f_120 * ig_203[k]
                  - f_121 * ig_331[k]
                  - f_121 * ig_336[k]
                  + f_122 * ig_338[k]
                  + f_123 * ig_361[k]
                  + f_123 * ig_366[k]
                  - f_124 * ig_368[k];
    }

#pragma omp simd aligned(ig_64, ig_71, ig_73, ig_169, ig_176, ig_178, ig_199, ig_206, ig_208, \
                         ig_334, ig_341, ig_343, ig_364, ig_371, \
                         ig_373 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = f_125 * ig_64[k]
                  + f_125 * ig_71[k]
                  - f_126 * ig_73[k]
                  + f_127 * ig_169[k]
                  + f_127 * ig_176[k]
                  - f_128 * ig_178[k]
                  - f_129 * ig_199[k]
                  - f_129 * ig_206[k]
                  + f_130 * ig_208[k]
                  - f_131 * ig_334[k]
                  - f_131 * ig_341[k]
                  + f_132 * ig_343[k]
                  + f_128 * ig_364[k]
                  + f_128 * ig_371[k]
                  - f_133 * ig_373[k];
    }

#pragma omp simd aligned(ig_60, ig_63, ig_65, ig_70, ig_72, ig_74, ig_165, ig_168, ig_170, \
                         ig_175, ig_177, ig_179, ig_195, ig_198, ig_200, ig_205, ig_207, \
                         ig_209, ig_330, ig_333, ig_335, ig_340, ig_342, ig_344, ig_360, \
                         ig_363, ig_365, ig_370, ig_372, ig_374 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = -f_134 * ig_60[k]
                  - f_135 * ig_63[k]
                  + f_136 * ig_65[k]
                  - f_134 * ig_70[k]
                  + f_136 * ig_72[k]
                  - f_137 * ig_74[k]
                  - f_138 * ig_165[k]
                  - f_139 * ig_168[k]
                  + f_140 * ig_170[k]
                  - f_138 * ig_175[k]
                  + f_140 * ig_177[k]
                  - f_141 * ig_179[k]
                  + f_137 * ig_195[k]
                  + f_140 * ig_198[k]
                  - f_142 * ig_200[k]
                  + f_137 * ig_205[k]
                  - f_142 * ig_207[k]
                  + f_143 * ig_209[k]
                  + f_144 * ig_330[k]
                  + f_138 * ig_333[k]
                  - f_137 * ig_335[k]
                  + f_144 * ig_340[k]
                  - f_137 * ig_342[k]
                  + f_145 * ig_344[k]
                  - f_145 * ig_360[k]
                  - f_141 * ig_363[k]
                  + f_143 * ig_365[k]
                  - f_145 * ig_370[k]
                  + f_143 * ig_372[k]
                  - f_146 * ig_374[k];
    }

#pragma omp simd aligned(ig_62, ig_67, ig_69, ig_167, ig_172, ig_174, ig_197, ig_202, ig_204, \
                         ig_332, ig_337, ig_339, ig_362, ig_367, \
                         ig_369 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = f_125 * ig_62[k]
                  + f_125 * ig_67[k]
                  - f_126 * ig_69[k]
                  + f_127 * ig_167[k]
                  + f_127 * ig_172[k]
                  - f_128 * ig_174[k]
                  - f_129 * ig_197[k]
                  - f_129 * ig_202[k]
                  + f_130 * ig_204[k]
                  - f_131 * ig_332[k]
                  - f_131 * ig_337[k]
                  + f_132 * ig_339[k]
                  + f_128 * ig_362[k]
                  + f_128 * ig_367[k]
                  - f_133 * ig_369[k];
    }

#pragma omp simd aligned(ig_60, ig_65, ig_70, ig_72, ig_165, ig_170, ig_175, ig_177, ig_195, \
                         ig_200, ig_205, ig_207, ig_330, ig_335, ig_340, ig_342, ig_360, \
                         ig_365, ig_370, ig_372 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = f_147 * ig_60[k]
                  - f_148 * ig_65[k]
                  - f_147 * ig_70[k]
                  + f_148 * ig_72[k]
                  + f_121 * ig_165[k]
                  - f_122 * ig_170[k]
                  - f_121 * ig_175[k]
                  + f_122 * ig_177[k]
                  - f_149 * ig_195[k]
                  + f_150 * ig_200[k]
                  + f_149 * ig_205[k]
                  - f_150 * ig_207[k]
                  - f_151 * ig_330[k]
                  + f_115 * ig_335[k]
                  + f_151 * ig_340[k]
                  - f_115 * ig_342[k]
                  + f_152 * ig_360[k]
                  - f_119 * ig_365[k]
                  - f_152 * ig_370[k]
                  + f_119 * ig_372[k];
    }

#pragma omp simd aligned(ig_62, ig_67, ig_167, ig_172, ig_197, ig_202, ig_332, ig_337, ig_362, \
                         ig_367 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = -f_108 * ig_62[k]
                  + f_107 * ig_67[k]
                  - f_110 * ig_167[k]
                  + f_109 * ig_172[k]
                  + f_112 * ig_197[k]
                  - f_111 * ig_202[k]
                  + f_113 * ig_332[k]
                  - f_108 * ig_337[k]
                  - f_114 * ig_362[k]
                  + f_112 * ig_367[k];
    }

#pragma omp simd aligned(ig_60, ig_63, ig_70, ig_165, ig_168, ig_175, ig_195, ig_198, ig_205, \
                         ig_330, ig_333, ig_340, ig_360, ig_363, \
                         ig_370 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = -f_153 * ig_60[k]
                  + f_154 * ig_63[k]
                  - f_153 * ig_70[k]
                  - f_155 * ig_165[k]
                  + f_102 * ig_168[k]
                  - f_155 * ig_175[k]
                  + f_103 * ig_195[k]
                  - f_156 * ig_198[k]
                  + f_103 * ig_205[k]
                  + f_157 * ig_330[k]
                  - f_158 * ig_333[k]
                  + f_157 * ig_340[k]
                  - f_159 * ig_360[k]
                  + f_160 * ig_363[k]
                  - f_159 * ig_370[k];
    }

#pragma omp simd aligned(ig_16, ig_21, ig_91, ig_96, ig_121, ig_126, ig_226, ig_231, ig_256, \
                         ig_261, ig_286, ig_291 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_36[k] = f_161 * ig_16[k]
                  - f_161 * ig_21[k]
                  + f_159 * ig_91[k]
                  - f_159 * ig_96[k]
                  - f_162 * ig_121[k]
                  + f_162 * ig_126[k]
                  + f_161 * ig_226[k]
                  - f_161 * ig_231[k]
                  - f_162 * ig_256[k]
                  + f_162 * ig_261[k]
                  + f_162 * ig_286[k]
                  - f_162 * ig_291[k];
    }

#pragma omp simd aligned(ig_19, ig_26, ig_94, ig_101, ig_124, ig_131, ig_229, ig_236, ig_259, \
                         ig_266, ig_289, ig_296 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_37[k] = f_113 * ig_19[k]
                  - f_163 * ig_26[k]
                  + f_110 * ig_94[k]
                  - f_164 * ig_101[k]
                  - f_165 * ig_124[k]
                  + f_166 * ig_131[k]
                  + f_113 * ig_229[k]
                  - f_163 * ig_236[k]
                  - f_165 * ig_259[k]
                  + f_166 * ig_266[k]
                  + f_165 * ig_289[k]
                  - f_166 * ig_296[k];
    }

#pragma omp simd aligned(ig_16, ig_21, ig_23, ig_91, ig_96, ig_98, ig_121, ig_126, ig_128, \
                         ig_226, ig_231, ig_233, ig_256, ig_261, ig_263, ig_286, ig_291, \
                         ig_293 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = -f_167 * ig_16[k]
                  - f_167 * ig_21[k]
                  + f_117 * ig_23[k]
                  - f_168 * ig_91[k]
                  - f_168 * ig_96[k]
                  + f_149 * ig_98[k]
                  + f_169 * ig_121[k]
                  + f_169 * ig_126[k]
                  - f_170 * ig_128[k]
                  - f_167 * ig_226[k]
                  - f_167 * ig_231[k]
                  + f_117 * ig_233[k]
                  + f_169 * ig_256[k]
                  + f_169 * ig_261[k]
                  - f_170 * ig_263[k]
                  - f_169 * ig_286[k]
                  - f_169 * ig_291[k]
                  + f_170 * ig_293[k];
    }

#pragma omp simd aligned(ig_19, ig_26, ig_28, ig_94, ig_101, ig_103, ig_124, ig_131, ig_133, \
                         ig_229, ig_236, ig_238, ig_259, ig_266, ig_268, ig_289, ig_296, \
                         ig_298 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_39[k] = -f_171 * ig_19[k]
                  - f_171 * ig_26[k]
                  + f_172 * ig_28[k]
                  - f_173 * ig_94[k]
                  - f_173 * ig_101[k]
                  + f_174 * ig_103[k]
                  + f_175 * ig_124[k]
                  + f_175 * ig_131[k]
                  - f_176 * ig_133[k]
                  - f_171 * ig_229[k]
                  - f_171 * ig_236[k]
                  + f_172 * ig_238[k]
                  + f_175 * ig_259[k]
                  + f_175 * ig_266[k]
                  - f_176 * ig_268[k]
                  - f_175 * ig_289[k]
                  - f_175 * ig_296[k]
                  + f_176 * ig_298[k];
    }

#pragma omp simd aligned(ig_15, ig_18, ig_20, ig_25, ig_27, ig_29, ig_90, ig_93, ig_95, \
                         ig_100, ig_102, ig_104, ig_120, ig_123, ig_125, ig_130, ig_132, \
                         ig_134, ig_225, ig_228, ig_230, ig_235, ig_237, ig_239, ig_255, \
                         ig_258, ig_260, ig_265, ig_267, ig_269, ig_285, ig_288, ig_290, \
                         ig_295, ig_297, ig_299 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_40[k] = f_177 * ig_15[k]
                  + f_178 * ig_18[k]
                  - f_145 * ig_20[k]
                  + f_177 * ig_25[k]
                  - f_145 * ig_27[k]
                  + f_179 * ig_29[k]
                  + f_178 * ig_90[k]
                  + f_180 * ig_93[k]
                  - f_141 * ig_95[k]
                  + f_178 * ig_100[k]
                  - f_141 * ig_102[k]
                  + f_181 * ig_104[k]
                  - f_141 * ig_120[k]
                  - f_182 * ig_123[k]
                  + f_183 * ig_125[k]
                  - f_141 * ig_130[k]
                  + f_183 * ig_132[k]
                  - f_184 * ig_134[k]
                  + f_177 * ig_225[k]
                  + f_178 * ig_228[k]
                  - f_145 * ig_230[k]
                  + f_177 * ig_235[k]
                  - f_145 * ig_237[k]
                  + f_179 * ig_239[k]
                  - f_141 * ig_255[k]
                  - f_182 * ig_258[k]
                  + f_183 * ig_260[k]
                  - f_141 * ig_265[k]
                  + f_183 * ig_267[k]
                  - f_184 * ig_269[k]
                  + f_141 * ig_285[k]
                  + f_182 * ig_288[k]
                  - f_183 * ig_290[k]
                  + f_141 * ig_295[k]
                  - f_183 * ig_297[k]
                  + f_184 * ig_299[k];
    }

#pragma omp simd aligned(ig_17, ig_22, ig_24, ig_92, ig_97, ig_99, ig_122, ig_127, ig_129, \
                         ig_227, ig_232, ig_234, ig_257, ig_262, ig_264, ig_287, ig_292, \
                         ig_294 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_41[k] = -f_171 * ig_17[k]
                  - f_171 * ig_22[k]
                  + f_172 * ig_24[k]
                  - f_173 * ig_92[k]
                  - f_173 * ig_97[k]
                  + f_174 * ig_99[k]
                  + f_175 * ig_122[k]
                  + f_175 * ig_127[k]
                  - f_176 * ig_129[k]
                  - f_171 * ig_227[k]
                  - f_171 * ig_232[k]
                  + f_172 * ig_234[k]
                  + f_175 * ig_257[k]
                  + f_175 * ig_262[k]
                  - f_176 * ig_264[k]
                  - f_175 * ig_287[k]
                  - f_175 * ig_292[k]
                  + f_176 * ig_294[k];
    }

#pragma omp simd aligned(ig_15, ig_20, ig_25, ig_27, ig_90, ig_95, ig_100, ig_102, ig_120, \
                         ig_125, ig_130, ig_132, ig_225, ig_230, ig_235, ig_237, ig_255, \
                         ig_260, ig_265, ig_267, ig_285, ig_290, ig_295, \
                         ig_297 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = -f_185 * ig_15[k]
                  + f_121 * ig_20[k]
                  + f_185 * ig_25[k]
                  - f_121 * ig_27[k]
                  - f_167 * ig_90[k]
                  + f_117 * ig_95[k]
                  + f_167 * ig_100[k]
                  - f_117 * ig_102[k]
                  + f_123 * ig_120[k]
                  - f_124 * ig_125[k]
                  - f_123 * ig_130[k]
                  + f_124 * ig_132[k]
                  - f_185 * ig_225[k]
                  + f_121 * ig_230[k]
                  + f_185 * ig_235[k]
                  - f_121 * ig_237[k]
                  + f_123 * ig_255[k]
                  - f_124 * ig_260[k]
                  - f_123 * ig_265[k]
                  + f_124 * ig_267[k]
                  - f_123 * ig_285[k]
                  + f_124 * ig_290[k]
                  + f_123 * ig_295[k]
                  - f_124 * ig_297[k];
    }

#pragma omp simd aligned(ig_17, ig_22, ig_92, ig_97, ig_122, ig_127, ig_227, ig_232, ig_257, \
                         ig_262, ig_287, ig_292 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_43[k] = f_163 * ig_17[k]
                  - f_113 * ig_22[k]
                  + f_164 * ig_92[k]
                  - f_110 * ig_97[k]
                  - f_166 * ig_122[k]
                  + f_165 * ig_127[k]
                  + f_163 * ig_227[k]
                  - f_113 * ig_232[k]
                  - f_166 * ig_257[k]
                  + f_165 * ig_262[k]
                  + f_166 * ig_287[k]
                  - f_165 * ig_292[k];
    }

#pragma omp simd aligned(ig_15, ig_18, ig_25, ig_90, ig_93, ig_100, ig_120, ig_123, ig_130, \
                         ig_225, ig_228, ig_235, ig_255, ig_258, ig_265, ig_285, ig_288, \
                         ig_295 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_44[k] = f_186 * ig_15[k]
                  - f_155 * ig_18[k]
                  + f_186 * ig_25[k]
                  + f_187 * ig_90[k]
                  - f_105 * ig_93[k]
                  + f_187 * ig_100[k]
                  - f_188 * ig_120[k]
                  + f_104 * ig_123[k]
                  - f_188 * ig_130[k]
                  + f_186 * ig_225[k]
                  - f_155 * ig_228[k]
                  + f_186 * ig_235[k]
                  - f_188 * ig_255[k]
                  + f_104 * ig_258[k]
                  - f_188 * ig_265[k]
                  + f_188 * ig_285[k]
                  - f_104 * ig_288[k]
                  + f_188 * ig_295[k];
    }

#pragma omp simd aligned(ig_61, ig_66, ig_166, ig_171, ig_196, ig_201, ig_331, ig_336, ig_361, \
                         ig_366, ig_391, ig_396 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_45[k] = f_189 * ig_61[k]
                  - f_189 * ig_66[k]
                  + f_190 * ig_166[k]
                  - f_190 * ig_171[k]
                  - f_191 * ig_196[k]
                  + f_191 * ig_201[k]
                  + f_189 * ig_331[k]
                  - f_189 * ig_336[k]
                  - f_191 * ig_361[k]
                  + f_191 * ig_366[k]
                  + f_192 * ig_391[k]
                  - f_192 * ig_396[k];
    }

#pragma omp simd aligned(ig_64, ig_71, ig_169, ig_176, ig_199, ig_206, ig_334, ig_341, ig_364, \
                         ig_371, ig_394, ig_401 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_46[k] = f_193 * ig_64[k]
                  - f_194 * ig_71[k]
                  + f_195 * ig_169[k]
                  - f_196 * ig_176[k]
                  - f_197 * ig_199[k]
                  + f_198 * ig_206[k]
                  + f_193 * ig_334[k]
                  - f_194 * ig_341[k]
                  - f_197 * ig_364[k]
                  + f_198 * ig_371[k]
                  + f_199 * ig_394[k]
                  - f_200 * ig_401[k];
    }

#pragma omp simd aligned(ig_61, ig_66, ig_68, ig_166, ig_171, ig_173, ig_196, ig_201, ig_203, \
                         ig_331, ig_336, ig_338, ig_361, ig_366, ig_368, ig_391, ig_396, \
                         ig_398 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_47[k] = -f_201 * ig_61[k]
                  - f_201 * ig_66[k]
                  + f_202 * ig_68[k]
                  - f_203 * ig_166[k]
                  - f_203 * ig_171[k]
                  + f_204 * ig_173[k]
                  + f_205 * ig_196[k]
                  + f_205 * ig_201[k]
                  - f_206 * ig_203[k]
                  - f_201 * ig_331[k]
                  - f_201 * ig_336[k]
                  + f_202 * ig_338[k]
                  + f_205 * ig_361[k]
                  + f_205 * ig_366[k]
                  - f_206 * ig_368[k]
                  - f_207 * ig_391[k]
                  - f_207 * ig_396[k]
                  + f_208 * ig_398[k];
    }

#pragma omp simd aligned(ig_64, ig_71, ig_73, ig_169, ig_176, ig_178, ig_199, ig_206, ig_208, \
                         ig_334, ig_341, ig_343, ig_364, ig_371, ig_373, ig_394, ig_401, \
                         ig_403 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_48[k] = -f_209 * ig_64[k]
                  - f_209 * ig_71[k]
                  + f_210 * ig_73[k]
                  - f_211 * ig_169[k]
                  - f_211 * ig_176[k]
                  + f_212 * ig_178[k]
                  + f_213 * ig_199[k]
                  + f_213 * ig_206[k]
                  - f_214 * ig_208[k]
                  - f_209 * ig_334[k]
                  - f_209 * ig_341[k]
                  + f_210 * ig_343[k]
                  + f_213 * ig_364[k]
                  + f_213 * ig_371[k]
                  - f_214 * ig_373[k]
                  - f_182 * ig_394[k]
                  - f_182 * ig_401[k]
                  + f_184 * ig_403[k];
    }

#pragma omp simd aligned(ig_60, ig_63, ig_65, ig_70, ig_72, ig_74, ig_165, ig_168, ig_170, \
                         ig_175, ig_177, ig_179, ig_195, ig_198, ig_200, ig_205, ig_207, \
                         ig_209, ig_330, ig_333, ig_335, ig_340, ig_342, ig_344, ig_360, \
                         ig_363, ig_365, ig_370, ig_372, ig_374, ig_390, ig_393, ig_395, \
                         ig_400, ig_402, ig_404 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_49[k] = f_215 * ig_60[k]
                  + f_171 * ig_63[k]
                  - f_132 * ig_65[k]
                  + f_215 * ig_70[k]
                  - f_132 * ig_72[k]
                  + f_172 * ig_74[k]
                  + f_171 * ig_165[k]
                  + f_173 * ig_168[k]
                  - f_128 * ig_170[k]
                  + f_171 * ig_175[k]
                  - f_128 * ig_177[k]
                  + f_174 * ig_179[k]
                  - f_173 * ig_195[k]
                  - f_132 * ig_198[k]
                  + f_175 * ig_200[k]
                  - f_173 * ig_205[k]
                  + f_175 * ig_207[k]
                  - f_216 * ig_209[k]
                  + f_215 * ig_330[k]
                  + f_171 * ig_333[k]
                  - f_132 * ig_335[k]
                  + f_215 * ig_340[k]
                  - f_132 * ig_342[k]
                  + f_172 * ig_344[k]
                  - f_173 * ig_360[k]
                  - f_132 * ig_363[k]
                  + f_175 * ig_365[k]
                  - f_173 * ig_370[k]
                  + f_175 * ig_372[k]
                  - f_216 * ig_374[k]
                  + f_217 * ig_390[k]
                  + f_218 * ig_393[k]
                  - f_219 * ig_395[k]
                  + f_217 * ig_400[k]
                  - f_219 * ig_402[k]
                  + f_220 * ig_404[k];
    }

#pragma omp simd aligned(ig_62, ig_67, ig_69, ig_167, ig_172, ig_174, ig_197, ig_202, ig_204, \
                         ig_332, ig_337, ig_339, ig_362, ig_367, ig_369, ig_392, ig_397, \
                         ig_399 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_50[k] = -f_209 * ig_62[k]
                  - f_209 * ig_67[k]
                  + f_210 * ig_69[k]
                  - f_211 * ig_167[k]
                  - f_211 * ig_172[k]
                  + f_212 * ig_174[k]
                  + f_213 * ig_197[k]
                  + f_213 * ig_202[k]
                  - f_214 * ig_204[k]
                  - f_209 * ig_332[k]
                  - f_209 * ig_337[k]
                  + f_210 * ig_339[k]
                  + f_213 * ig_362[k]
                  + f_213 * ig_367[k]
                  - f_214 * ig_369[k]
                  - f_182 * ig_392[k]
                  - f_182 * ig_397[k]
                  + f_184 * ig_399[k];
    }

#pragma omp simd aligned(ig_60, ig_65, ig_70, ig_72, ig_165, ig_170, ig_175, ig_177, ig_195, \
                         ig_200, ig_205, ig_207, ig_330, ig_335, ig_340, ig_342, ig_360, \
                         ig_365, ig_370, ig_372, ig_390, ig_395, ig_400, \
                         ig_402 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_51[k] = -f_221 * ig_60[k]
                  + f_222 * ig_65[k]
                  + f_221 * ig_70[k]
                  - f_222 * ig_72[k]
                  - f_201 * ig_165[k]
                  + f_202 * ig_170[k]
                  + f_201 * ig_175[k]
                  - f_202 * ig_177[k]
                  + f_203 * ig_195[k]
                  - f_204 * ig_200[k]
                  - f_203 * ig_205[k]
                  + f_204 * ig_207[k]
                  - f_221 * ig_330[k]
                  + f_222 * ig_335[k]
                  + f_221 * ig_340[k]
                  - f_222 * ig_342[k]
                  + f_203 * ig_360[k]
                  - f_204 * ig_365[k]
                  - f_203 * ig_370[k]
                  + f_204 * ig_372[k]
                  - f_223 * ig_390[k]
                  + f_224 * ig_395[k]
                  + f_223 * ig_400[k]
                  - f_224 * ig_402[k];
    }

#pragma omp simd aligned(ig_62, ig_67, ig_167, ig_172, ig_197, ig_202, ig_332, ig_337, ig_362, \
                         ig_367, ig_392, ig_397 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_52[k] = f_194 * ig_62[k]
                  - f_193 * ig_67[k]
                  + f_196 * ig_167[k]
                  - f_195 * ig_172[k]
                  - f_198 * ig_197[k]
                  + f_197 * ig_202[k]
                  + f_194 * ig_332[k]
                  - f_193 * ig_337[k]
                  - f_198 * ig_362[k]
                  + f_197 * ig_367[k]
                  + f_200 * ig_392[k]
                  - f_199 * ig_397[k];
    }

#pragma omp simd aligned(ig_60, ig_63, ig_70, ig_165, ig_168, ig_175, ig_195, ig_198, ig_205, \
                         ig_330, ig_333, ig_340, ig_360, ig_363, ig_370, ig_390, ig_393, \
                         ig_400 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_53[k] = f_225 * ig_60[k]
                  - f_226 * ig_63[k]
                  + f_225 * ig_70[k]
                  + f_227 * ig_165[k]
                  - f_228 * ig_168[k]
                  + f_227 * ig_175[k]
                  - f_189 * ig_195[k]
                  + f_229 * ig_198[k]
                  - f_189 * ig_205[k]
                  + f_225 * ig_330[k]
                  - f_226 * ig_333[k]
                  + f_225 * ig_340[k]
                  - f_189 * ig_360[k]
                  + f_229 * ig_363[k]
                  - f_189 * ig_370[k]
                  + f_230 * ig_390[k]
                  - f_231 * ig_393[k]
                  + f_230 * ig_400[k];
    }

#pragma omp simd aligned(ig_1, ig_6, ig_46, ig_51, ig_76, ig_81, ig_151, ig_156, ig_181, \
                         ig_186, ig_211, ig_216, ig_316, ig_321, ig_346, ig_351, ig_376, \
                         ig_381, ig_406, ig_411 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_54[k] = -f_232 * ig_1[k]
                  + f_232 * ig_6[k]
                  - f_233 * ig_46[k]
                  + f_233 * ig_51[k]
                  + f_234 * ig_76[k]
                  - f_234 * ig_81[k]
                  - f_233 * ig_151[k]
                  + f_233 * ig_156[k]
                  + f_235 * ig_181[k]
                  - f_235 * ig_186[k]
                  - f_80 * ig_211[k]
                  + f_80 * ig_216[k]
                  - f_232 * ig_316[k]
                  + f_232 * ig_321[k]
                  + f_234 * ig_346[k]
                  - f_234 * ig_351[k]
                  - f_80 * ig_376[k]
                  + f_80 * ig_381[k]
                  + f_236 * ig_406[k]
                  - f_236 * ig_411[k];
    }

#pragma omp simd aligned(ig_4, ig_11, ig_49, ig_56, ig_79, ig_86, ig_154, ig_161, ig_184, \
                         ig_191, ig_214, ig_221, ig_319, ig_326, ig_349, ig_356, ig_379, \
                         ig_386, ig_409, ig_416 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_55[k] = -f_237 * ig_4[k]
                  + f_238 * ig_11[k]
                  - f_239 * ig_49[k]
                  + f_237 * ig_56[k]
                  + f_240 * ig_79[k]
                  - f_241 * ig_86[k]
                  - f_239 * ig_154[k]
                  + f_237 * ig_161[k]
                  + f_242 * ig_184[k]
                  - f_243 * ig_191[k]
                  - f_84 * ig_214[k]
                  + f_244 * ig_221[k]
                  - f_237 * ig_319[k]
                  + f_238 * ig_326[k]
                  + f_240 * ig_349[k]
                  - f_241 * ig_356[k]
                  - f_84 * ig_379[k]
                  + f_244 * ig_386[k]
                  + f_83 * ig_409[k]
                  - f_245 * ig_416[k];
    }

#pragma omp simd aligned(ig_1, ig_6, ig_8, ig_46, ig_51, ig_53, ig_76, ig_81, ig_83, ig_151, \
                         ig_156, ig_158, ig_181, ig_186, ig_188, ig_211, ig_216, ig_218, \
                         ig_316, ig_321, ig_323, ig_346, ig_351, ig_353, ig_376, ig_381, \
                         ig_383, ig_406, ig_411, ig_413 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_56[k] = f_246 * ig_1[k]
                  + f_246 * ig_6[k]
                  - f_247 * ig_8[k]
                  + f_248 * ig_46[k]
                  + f_248 * ig_51[k]
                  - f_249 * ig_53[k]
                  - f_249 * ig_76[k]
                  - f_249 * ig_81[k]
                  + f_250 * ig_83[k]
                  + f_248 * ig_151[k]
                  + f_248 * ig_156[k]
                  - f_249 * ig_158[k]
                  - f_251 * ig_181[k]
                  - f_251 * ig_186[k]
                  + f_252 * ig_188[k]
                  + f_253 * ig_211[k]
                  + f_253 * ig_216[k]
                  - f_254 * ig_218[k]
                  + f_246 * ig_316[k]
                  + f_246 * ig_321[k]
                  - f_247 * ig_323[k]
                  - f_249 * ig_346[k]
                  - f_249 * ig_351[k]
                  + f_250 * ig_353[k]
                  + f_253 * ig_376[k]
                  + f_253 * ig_381[k]
                  - f_254 * ig_383[k]
                  - f_255 * ig_406[k]
                  - f_255 * ig_411[k]
                  + f_256 * ig_413[k];
    }

#pragma omp simd aligned(ig_4, ig_11, ig_13, ig_49, ig_56, ig_58, ig_79, ig_86, ig_88, ig_154, \
                         ig_161, ig_163, ig_184, ig_191, ig_193, ig_214, ig_221, ig_223, \
                         ig_319, ig_326, ig_328, ig_349, ig_356, ig_358, ig_379, ig_386, \
                         ig_388, ig_409, ig_416, ig_418 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_57[k] = f_257 * ig_4[k]
                  + f_257 * ig_11[k]
                  - f_258 * ig_13[k]
                  + f_259 * ig_49[k]
                  + f_259 * ig_56[k]
                  - f_260 * ig_58[k]
                  - f_261 * ig_79[k]
                  - f_261 * ig_86[k]
                  + f_262 * ig_88[k]
                  + f_259 * ig_154[k]
                  + f_259 * ig_161[k]
                  - f_260 * ig_163[k]
                  - f_263 * ig_184[k]
                  - f_263 * ig_191[k]
                  + f_264 * ig_193[k]
                  + f_262 * ig_214[k]
                  + f_262 * ig_221[k]
                  - f_265 * ig_223[k]
                  + f_257 * ig_319[k]
                  + f_257 * ig_326[k]
                  - f_258 * ig_328[k]
                  - f_261 * ig_349[k]
                  - f_261 * ig_356[k]
                  + f_262 * ig_358[k]
                  + f_262 * ig_379[k]
                  + f_262 * ig_386[k]
                  - f_265 * ig_388[k]
                  - f_266 * ig_409[k]
                  - f_266 * ig_416[k]
                  + f_267 * ig_418[k];
    }

#pragma omp simd aligned(ig_0, ig_3, ig_5, ig_10, ig_12, ig_14, ig_45, ig_48, ig_50, ig_55, \
                         ig_57, ig_59, ig_75, ig_78, ig_80, ig_85, ig_87, ig_89, ig_150, \
                         ig_153, ig_155, ig_160, ig_162, ig_164, ig_180, ig_183, ig_185, \
                         ig_190, ig_192, ig_194, ig_210, ig_213, ig_215, ig_220, ig_222, \
                         ig_224, ig_315, ig_318, ig_320, ig_325, ig_327, ig_329, ig_345, \
                         ig_348, ig_350, ig_355, ig_357, ig_359, ig_375, ig_378, ig_380, \
                         ig_385, ig_387, ig_389, ig_405, ig_408, ig_410, ig_415, ig_417, \
                         ig_419 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_58[k] = -0.1171875 * ig_0[k]
                  - 0.234375 * ig_3[k]
                  + 0.9375 * ig_5[k]
                  - 0.1171875 * ig_10[k]
                  + 0.9375 * ig_12[k]
                  - 0.3125 * ig_14[k]
                  - 0.3515625 * ig_45[k]
                  - 0.703125 * ig_48[k]
                  + 2.8125 * ig_50[k]
                  - 0.3515625 * ig_55[k]
                  + 2.8125 * ig_57[k]
                  - 0.9375 * ig_59[k]
                  + 2.109375 * ig_75[k]
                  + 4.21875 * ig_78[k]
                  - 16.875 * ig_80[k]
                  + 2.109375 * ig_85[k]
                  - 16.875 * ig_87[k]
                  + 5.625 * ig_89[k]
                  - 0.3515625 * ig_150[k]
                  - 0.703125 * ig_153[k]
                  + 2.8125 * ig_155[k]
                  - 0.3515625 * ig_160[k]
                  + 2.8125 * ig_162[k]
                  - 0.9375 * ig_164[k]
                  + 4.21875 * ig_180[k]
                  + 8.4375 * ig_183[k]
                  - 33.75 * ig_185[k]
                  + 4.21875 * ig_190[k]
                  - 33.75 * ig_192[k]
                  + 11.25 * ig_194[k]
                  - 2.8125 * ig_210[k]
                  - 5.625 * ig_213[k]
                  + 22.5 * ig_215[k]
                  - 2.8125 * ig_220[k]
                  + 22.5 * ig_222[k]
                  - 7.5 * ig_224[k]
                  - 0.1171875 * ig_315[k]
                  - 0.234375 * ig_318[k]
                  + 0.9375 * ig_320[k]
                  - 0.1171875 * ig_325[k]
                  + 0.9375 * ig_327[k]
                  - 0.3125 * ig_329[k]
                  + 2.109375 * ig_345[k]
                  + 4.21875 * ig_348[k]
                  - 16.875 * ig_350[k]
                  + 2.109375 * ig_355[k]
                  - 16.875 * ig_357[k]
                  + 5.625 * ig_359[k]
                  - 2.8125 * ig_375[k]
                  - 5.625 * ig_378[k]
                  + 22.5 * ig_380[k]
                  - 2.8125 * ig_385[k]
                  + 22.5 * ig_387[k]
                  - 7.5 * ig_389[k]
                  + 0.375 * ig_405[k]
                  + 0.75 * ig_408[k]
                  - 3.0 * ig_410[k]
                  + 0.375 * ig_415[k]
                  - 3.0 * ig_417[k]
                  + ig_419[k];
    }

#pragma omp simd aligned(ig_2, ig_7, ig_9, ig_47, ig_52, ig_54, ig_77, ig_82, ig_84, ig_152, \
                         ig_157, ig_159, ig_182, ig_187, ig_189, ig_212, ig_217, ig_219, \
                         ig_317, ig_322, ig_324, ig_347, ig_352, ig_354, ig_377, ig_382, \
                         ig_384, ig_407, ig_412, ig_414 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_59[k] = f_257 * ig_2[k]
                  + f_257 * ig_7[k]
                  - f_258 * ig_9[k]
                  + f_259 * ig_47[k]
                  + f_259 * ig_52[k]
                  - f_260 * ig_54[k]
                  - f_261 * ig_77[k]
                  - f_261 * ig_82[k]
                  + f_262 * ig_84[k]
                  + f_259 * ig_152[k]
                  + f_259 * ig_157[k]
                  - f_260 * ig_159[k]
                  - f_263 * ig_182[k]
                  - f_263 * ig_187[k]
                  + f_264 * ig_189[k]
                  + f_262 * ig_212[k]
                  + f_262 * ig_217[k]
                  - f_265 * ig_219[k]
                  + f_257 * ig_317[k]
                  + f_257 * ig_322[k]
                  - f_258 * ig_324[k]
                  - f_261 * ig_347[k]
                  - f_261 * ig_352[k]
                  + f_262 * ig_354[k]
                  + f_262 * ig_377[k]
                  + f_262 * ig_382[k]
                  - f_265 * ig_384[k]
                  - f_266 * ig_407[k]
                  - f_266 * ig_412[k]
                  + f_267 * ig_414[k];
    }

#pragma omp simd aligned(ig_0, ig_5, ig_10, ig_12, ig_45, ig_50, ig_55, ig_57, ig_75, ig_80, \
                         ig_85, ig_87, ig_150, ig_155, ig_160, ig_162, ig_180, ig_185, ig_190, \
                         ig_192, ig_210, ig_215, ig_220, ig_222, ig_315, ig_320, ig_325, \
                         ig_327, ig_345, ig_350, ig_355, ig_357, ig_375, ig_380, ig_385, \
                         ig_387, ig_405, ig_410, ig_415, ig_417 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_60[k] = f_268 * ig_0[k]
                  - f_248 * ig_5[k]
                  - f_268 * ig_10[k]
                  + f_248 * ig_12[k]
                  + f_269 * ig_45[k]
                  - f_270 * ig_50[k]
                  - f_269 * ig_55[k]
                  + f_270 * ig_57[k]
                  - f_270 * ig_75[k]
                  + f_271 * ig_80[k]
                  + f_270 * ig_85[k]
                  - f_271 * ig_87[k]
                  + f_269 * ig_150[k]
                  - f_270 * ig_155[k]
                  - f_269 * ig_160[k]
                  + f_270 * ig_162[k]
                  - f_249 * ig_180[k]
                  + f_250 * ig_185[k]
                  + f_249 * ig_190[k]
                  - f_250 * ig_192[k]
                  + f_272 * ig_210[k]
                  - f_273 * ig_215[k]
                  - f_272 * ig_220[k]
                  + f_273 * ig_222[k]
                  + f_268 * ig_315[k]
                  - f_248 * ig_320[k]
                  - f_268 * ig_325[k]
                  + f_248 * ig_327[k]
                  - f_270 * ig_345[k]
                  + f_271 * ig_350[k]
                  + f_270 * ig_355[k]
                  - f_271 * ig_357[k]
                  + f_272 * ig_375[k]
                  - f_273 * ig_380[k]
                  - f_272 * ig_385[k]
                  + f_273 * ig_387[k]
                  - f_274 * ig_405[k]
                  + f_275 * ig_410[k]
                  + f_274 * ig_415[k]
                  - f_275 * ig_417[k];
    }

#pragma omp simd aligned(ig_2, ig_7, ig_47, ig_52, ig_77, ig_82, ig_152, ig_157, ig_182, \
                         ig_187, ig_212, ig_217, ig_317, ig_322, ig_347, ig_352, ig_377, \
                         ig_382, ig_407, ig_412 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_61[k] = -f_238 * ig_2[k]
                  + f_237 * ig_7[k]
                  - f_237 * ig_47[k]
                  + f_239 * ig_52[k]
                  + f_241 * ig_77[k]
                  - f_240 * ig_82[k]
                  - f_237 * ig_152[k]
                  + f_239 * ig_157[k]
                  + f_243 * ig_182[k]
                  - f_242 * ig_187[k]
                  - f_244 * ig_212[k]
                  + f_84 * ig_217[k]
                  - f_238 * ig_317[k]
                  + f_237 * ig_322[k]
                  + f_241 * ig_347[k]
                  - f_240 * ig_352[k]
                  - f_244 * ig_377[k]
                  + f_84 * ig_382[k]
                  + f_245 * ig_407[k]
                  - f_83 * ig_412[k];
    }

#pragma omp simd aligned(ig_0, ig_3, ig_10, ig_45, ig_48, ig_55, ig_75, ig_78, ig_85, ig_150, \
                         ig_153, ig_160, ig_180, ig_183, ig_190, ig_210, ig_213, ig_220, \
                         ig_315, ig_318, ig_325, ig_345, ig_348, ig_355, ig_375, ig_378, \
                         ig_385, ig_405, ig_408, ig_415 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_62[k] = -f_276 * ig_0[k]
                  + f_277 * ig_3[k]
                  - f_276 * ig_10[k]
                  - f_278 * ig_45[k]
                  + f_279 * ig_48[k]
                  - f_278 * ig_55[k]
                  + f_279 * ig_75[k]
                  - f_280 * ig_78[k]
                  + f_279 * ig_85[k]
                  - f_278 * ig_150[k]
                  + f_279 * ig_153[k]
                  - f_278 * ig_160[k]
                  + f_281 * ig_180[k]
                  - f_282 * ig_183[k]
                  + f_281 * ig_190[k]
                  - f_283 * ig_210[k]
                  + f_235 * ig_213[k]
                  - f_283 * ig_220[k]
                  - f_276 * ig_315[k]
                  + f_277 * ig_318[k]
                  - f_276 * ig_325[k]
                  + f_279 * ig_345[k]
                  - f_280 * ig_348[k]
                  + f_279 * ig_355[k]
                  - f_283 * ig_375[k]
                  + f_235 * ig_378[k]
                  - f_283 * ig_385[k]
                  + f_284 * ig_405[k]
                  - f_285 * ig_408[k]
                  + f_284 * ig_415[k];
    }

#pragma omp simd aligned(ig_31, ig_36, ig_106, ig_111, ig_136, ig_141, ig_241, ig_246, ig_271, \
                         ig_276, ig_301, ig_306 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_63[k] = f_189 * ig_31[k]
                  - f_189 * ig_36[k]
                  + f_190 * ig_106[k]
                  - f_190 * ig_111[k]
                  - f_191 * ig_136[k]
                  + f_191 * ig_141[k]
                  + f_189 * ig_241[k]
                  - f_189 * ig_246[k]
                  - f_191 * ig_271[k]
                  + f_191 * ig_276[k]
                  + f_192 * ig_301[k]
                  - f_192 * ig_306[k];
    }

#pragma omp simd aligned(ig_34, ig_41, ig_109, ig_116, ig_139, ig_146, ig_244, ig_251, ig_274, \
                         ig_281, ig_304, ig_311 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_64[k] = f_193 * ig_34[k]
                  - f_194 * ig_41[k]
                  + f_195 * ig_109[k]
                  - f_196 * ig_116[k]
                  - f_197 * ig_139[k]
                  + f_198 * ig_146[k]
                  + f_193 * ig_244[k]
                  - f_194 * ig_251[k]
                  - f_197 * ig_274[k]
                  + f_198 * ig_281[k]
                  + f_199 * ig_304[k]
                  - f_200 * ig_311[k];
    }

#pragma omp simd aligned(ig_31, ig_36, ig_38, ig_106, ig_111, ig_113, ig_136, ig_141, ig_143, \
                         ig_241, ig_246, ig_248, ig_271, ig_276, ig_278, ig_301, ig_306, \
                         ig_308 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_65[k] = -f_201 * ig_31[k]
                  - f_201 * ig_36[k]
                  + f_202 * ig_38[k]
                  - f_203 * ig_106[k]
                  - f_203 * ig_111[k]
                  + f_204 * ig_113[k]
                  + f_205 * ig_136[k]
                  + f_205 * ig_141[k]
                  - f_206 * ig_143[k]
                  - f_201 * ig_241[k]
                  - f_201 * ig_246[k]
                  + f_202 * ig_248[k]
                  + f_205 * ig_271[k]
                  + f_205 * ig_276[k]
                  - f_206 * ig_278[k]
                  - f_207 * ig_301[k]
                  - f_207 * ig_306[k]
                  + f_208 * ig_308[k];
    }

#pragma omp simd aligned(ig_34, ig_41, ig_43, ig_109, ig_116, ig_118, ig_139, ig_146, ig_148, \
                         ig_244, ig_251, ig_253, ig_274, ig_281, ig_283, ig_304, ig_311, \
                         ig_313 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_66[k] = -f_209 * ig_34[k]
                  - f_209 * ig_41[k]
                  + f_210 * ig_43[k]
                  - f_211 * ig_109[k]
                  - f_211 * ig_116[k]
                  + f_212 * ig_118[k]
                  + f_213 * ig_139[k]
                  + f_213 * ig_146[k]
                  - f_214 * ig_148[k]
                  - f_209 * ig_244[k]
                  - f_209 * ig_251[k]
                  + f_210 * ig_253[k]
                  + f_213 * ig_274[k]
                  + f_213 * ig_281[k]
                  - f_214 * ig_283[k]
                  - f_182 * ig_304[k]
                  - f_182 * ig_311[k]
                  + f_184 * ig_313[k];
    }

#pragma omp simd aligned(ig_30, ig_33, ig_35, ig_40, ig_42, ig_44, ig_105, ig_108, ig_110, \
                         ig_115, ig_117, ig_119, ig_135, ig_138, ig_140, ig_145, ig_147, \
                         ig_149, ig_240, ig_243, ig_245, ig_250, ig_252, ig_254, ig_270, \
                         ig_273, ig_275, ig_280, ig_282, ig_284, ig_300, ig_303, ig_305, \
                         ig_310, ig_312, ig_314 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_67[k] = f_215 * ig_30[k]
                  + f_171 * ig_33[k]
                  - f_132 * ig_35[k]
                  + f_215 * ig_40[k]
                  - f_132 * ig_42[k]
                  + f_172 * ig_44[k]
                  + f_171 * ig_105[k]
                  + f_173 * ig_108[k]
                  - f_128 * ig_110[k]
                  + f_171 * ig_115[k]
                  - f_128 * ig_117[k]
                  + f_174 * ig_119[k]
                  - f_173 * ig_135[k]
                  - f_132 * ig_138[k]
                  + f_175 * ig_140[k]
                  - f_173 * ig_145[k]
                  + f_175 * ig_147[k]
                  - f_216 * ig_149[k]
                  + f_215 * ig_240[k]
                  + f_171 * ig_243[k]
                  - f_132 * ig_245[k]
                  + f_215 * ig_250[k]
                  - f_132 * ig_252[k]
                  + f_172 * ig_254[k]
                  - f_173 * ig_270[k]
                  - f_132 * ig_273[k]
                  + f_175 * ig_275[k]
                  - f_173 * ig_280[k]
                  + f_175 * ig_282[k]
                  - f_216 * ig_284[k]
                  + f_217 * ig_300[k]
                  + f_218 * ig_303[k]
                  - f_219 * ig_305[k]
                  + f_217 * ig_310[k]
                  - f_219 * ig_312[k]
                  + f_220 * ig_314[k];
    }

#pragma omp simd aligned(ig_32, ig_37, ig_39, ig_107, ig_112, ig_114, ig_137, ig_142, ig_144, \
                         ig_242, ig_247, ig_249, ig_272, ig_277, ig_279, ig_302, ig_307, \
                         ig_309 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_68[k] = -f_209 * ig_32[k]
                  - f_209 * ig_37[k]
                  + f_210 * ig_39[k]
                  - f_211 * ig_107[k]
                  - f_211 * ig_112[k]
                  + f_212 * ig_114[k]
                  + f_213 * ig_137[k]
                  + f_213 * ig_142[k]
                  - f_214 * ig_144[k]
                  - f_209 * ig_242[k]
                  - f_209 * ig_247[k]
                  + f_210 * ig_249[k]
                  + f_213 * ig_272[k]
                  + f_213 * ig_277[k]
                  - f_214 * ig_279[k]
                  - f_182 * ig_302[k]
                  - f_182 * ig_307[k]
                  + f_184 * ig_309[k];
    }

#pragma omp simd aligned(ig_30, ig_35, ig_40, ig_42, ig_105, ig_110, ig_115, ig_117, ig_135, \
                         ig_140, ig_145, ig_147, ig_240, ig_245, ig_250, ig_252, ig_270, \
                         ig_275, ig_280, ig_282, ig_300, ig_305, ig_310, \
                         ig_312 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_69[k] = -f_221 * ig_30[k]
                  + f_222 * ig_35[k]
                  + f_221 * ig_40[k]
                  - f_222 * ig_42[k]
                  - f_201 * ig_105[k]
                  + f_202 * ig_110[k]
                  + f_201 * ig_115[k]
                  - f_202 * ig_117[k]
                  + f_203 * ig_135[k]
                  - f_204 * ig_140[k]
                  - f_203 * ig_145[k]
                  + f_204 * ig_147[k]
                  - f_221 * ig_240[k]
                  + f_222 * ig_245[k]
                  + f_221 * ig_250[k]
                  - f_222 * ig_252[k]
                  + f_203 * ig_270[k]
                  - f_204 * ig_275[k]
                  - f_203 * ig_280[k]
                  + f_204 * ig_282[k]
                  - f_223 * ig_300[k]
                  + f_224 * ig_305[k]
                  + f_223 * ig_310[k]
                  - f_224 * ig_312[k];
    }

#pragma omp simd aligned(ig_32, ig_37, ig_107, ig_112, ig_137, ig_142, ig_242, ig_247, ig_272, \
                         ig_277, ig_302, ig_307 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_70[k] = f_194 * ig_32[k]
                  - f_193 * ig_37[k]
                  + f_196 * ig_107[k]
                  - f_195 * ig_112[k]
                  - f_198 * ig_137[k]
                  + f_197 * ig_142[k]
                  + f_194 * ig_242[k]
                  - f_193 * ig_247[k]
                  - f_198 * ig_272[k]
                  + f_197 * ig_277[k]
                  + f_200 * ig_302[k]
                  - f_199 * ig_307[k];
    }

#pragma omp simd aligned(ig_30, ig_33, ig_40, ig_105, ig_108, ig_115, ig_135, ig_138, ig_145, \
                         ig_240, ig_243, ig_250, ig_270, ig_273, ig_280, ig_300, ig_303, \
                         ig_310 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_71[k] = f_225 * ig_30[k]
                  - f_226 * ig_33[k]
                  + f_225 * ig_40[k]
                  + f_227 * ig_105[k]
                  - f_228 * ig_108[k]
                  + f_227 * ig_115[k]
                  - f_189 * ig_135[k]
                  + f_229 * ig_138[k]
                  - f_189 * ig_145[k]
                  + f_225 * ig_240[k]
                  - f_226 * ig_243[k]
                  + f_225 * ig_250[k]
                  - f_189 * ig_270[k]
                  + f_229 * ig_273[k]
                  - f_189 * ig_280[k]
                  + f_230 * ig_300[k]
                  - f_231 * ig_303[k]
                  + f_230 * ig_310[k];
    }

#pragma omp simd aligned(ig_1, ig_6, ig_46, ig_51, ig_76, ig_81, ig_151, ig_156, ig_211, \
                         ig_216, ig_316, ig_321, ig_346, ig_351, ig_376, \
                         ig_381 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_72[k] = f_187 * ig_1[k]
                  - f_187 * ig_6[k]
                  + f_187 * ig_46[k]
                  - f_187 * ig_51[k]
                  - f_106 * ig_76[k]
                  + f_106 * ig_81[k]
                  - f_187 * ig_151[k]
                  + f_187 * ig_156[k]
                  + f_106 * ig_211[k]
                  - f_106 * ig_216[k]
                  - f_187 * ig_316[k]
                  + f_187 * ig_321[k]
                  + f_106 * ig_346[k]
                  - f_106 * ig_351[k]
                  - f_106 * ig_376[k]
                  + f_106 * ig_381[k];
    }

#pragma omp simd aligned(ig_4, ig_11, ig_49, ig_56, ig_79, ig_86, ig_154, ig_161, ig_214, \
                         ig_221, ig_319, ig_326, ig_349, ig_356, ig_379, \
                         ig_386 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_73[k] = f_286 * ig_4[k]
                  - f_287 * ig_11[k]
                  + f_286 * ig_49[k]
                  - f_287 * ig_56[k]
                  - f_112 * ig_79[k]
                  + f_114 * ig_86[k]
                  - f_286 * ig_154[k]
                  + f_287 * ig_161[k]
                  + f_112 * ig_214[k]
                  - f_114 * ig_221[k]
                  - f_286 * ig_319[k]
                  + f_287 * ig_326[k]
                  + f_112 * ig_349[k]
                  - f_114 * ig_356[k]
                  - f_112 * ig_379[k]
                  + f_114 * ig_386[k];
    }

#pragma omp simd aligned(ig_1, ig_6, ig_8, ig_46, ig_51, ig_53, ig_76, ig_81, ig_83, ig_151, \
                         ig_156, ig_158, ig_211, ig_216, ig_218, ig_316, ig_321, ig_323, \
                         ig_346, ig_351, ig_353, ig_376, ig_381, \
                         ig_383 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_74[k] = -f_185 * ig_1[k]
                  - f_185 * ig_6[k]
                  + f_121 * ig_8[k]
                  - f_185 * ig_46[k]
                  - f_185 * ig_51[k]
                  + f_121 * ig_53[k]
                  + f_123 * ig_76[k]
                  + f_123 * ig_81[k]
                  - f_124 * ig_83[k]
                  + f_185 * ig_151[k]
                  + f_185 * ig_156[k]
                  - f_121 * ig_158[k]
                  - f_123 * ig_211[k]
                  - f_123 * ig_216[k]
                  + f_124 * ig_218[k]
                  + f_185 * ig_316[k]
                  + f_185 * ig_321[k]
                  - f_121 * ig_323[k]
                  - f_123 * ig_346[k]
                  - f_123 * ig_351[k]
                  + f_124 * ig_353[k]
                  + f_123 * ig_376[k]
                  + f_123 * ig_381[k]
                  - f_124 * ig_383[k];
    }

#pragma omp simd aligned(ig_4, ig_11, ig_13, ig_49, ig_56, ig_58, ig_79, ig_86, ig_88, ig_154, \
                         ig_161, ig_163, ig_214, ig_221, ig_223, ig_319, ig_326, ig_328, \
                         ig_349, ig_356, ig_358, ig_379, ig_386, \
                         ig_388 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_75[k] = -f_215 * ig_4[k]
                  - f_215 * ig_11[k]
                  + f_288 * ig_13[k]
                  - f_215 * ig_49[k]
                  - f_215 * ig_56[k]
                  + f_288 * ig_58[k]
                  + f_128 * ig_79[k]
                  + f_128 * ig_86[k]
                  - f_133 * ig_88[k]
                  + f_215 * ig_154[k]
                  + f_215 * ig_161[k]
                  - f_288 * ig_163[k]
                  - f_128 * ig_214[k]
                  - f_128 * ig_221[k]
                  + f_133 * ig_223[k]
                  + f_215 * ig_319[k]
                  + f_215 * ig_326[k]
                  - f_288 * ig_328[k]
                  - f_128 * ig_349[k]
                  - f_128 * ig_356[k]
                  + f_133 * ig_358[k]
                  + f_128 * ig_379[k]
                  + f_128 * ig_386[k]
                  - f_133 * ig_388[k];
    }

#pragma omp simd aligned(ig_0, ig_3, ig_5, ig_10, ig_12, ig_14, ig_45, ig_48, ig_50, ig_55, \
                         ig_57, ig_59, ig_75, ig_78, ig_80, ig_85, ig_87, ig_89, ig_150, \
                         ig_153, ig_155, ig_160, ig_162, ig_164, ig_210, ig_213, ig_215, \
                         ig_220, ig_222, ig_224, ig_315, ig_318, ig_320, ig_325, ig_327, \
                         ig_329, ig_345, ig_348, ig_350, ig_355, ig_357, ig_359, ig_375, \
                         ig_378, ig_380, ig_385, ig_387, ig_389 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_76[k] = f_289 * ig_0[k]
                  + f_177 * ig_3[k]
                  - f_180 * ig_5[k]
                  + f_289 * ig_10[k]
                  - f_180 * ig_12[k]
                  + f_290 * ig_14[k]
                  + f_289 * ig_45[k]
                  + f_177 * ig_48[k]
                  - f_180 * ig_50[k]
                  + f_289 * ig_55[k]
                  - f_180 * ig_57[k]
                  + f_290 * ig_59[k]
                  - f_145 * ig_75[k]
                  - f_141 * ig_78[k]
                  + f_143 * ig_80[k]
                  - f_145 * ig_85[k]
                  + f_143 * ig_87[k]
                  - f_146 * ig_89[k]
                  - f_289 * ig_150[k]
                  - f_177 * ig_153[k]
                  + f_180 * ig_155[k]
                  - f_289 * ig_160[k]
                  + f_180 * ig_162[k]
                  - f_290 * ig_164[k]
                  + f_145 * ig_210[k]
                  + f_141 * ig_213[k]
                  - f_143 * ig_215[k]
                  + f_145 * ig_220[k]
                  - f_143 * ig_222[k]
                  + f_146 * ig_224[k]
                  - f_289 * ig_315[k]
                  - f_177 * ig_318[k]
                  + f_180 * ig_320[k]
                  - f_289 * ig_325[k]
                  + f_180 * ig_327[k]
                  - f_290 * ig_329[k]
                  + f_145 * ig_345[k]
                  + f_141 * ig_348[k]
                  - f_143 * ig_350[k]
                  + f_145 * ig_355[k]
                  - f_143 * ig_357[k]
                  + f_146 * ig_359[k]
                  - f_145 * ig_375[k]
                  - f_141 * ig_378[k]
                  + f_143 * ig_380[k]
                  - f_145 * ig_385[k]
                  + f_143 * ig_387[k]
                  - f_146 * ig_389[k];
    }

#pragma omp simd aligned(ig_2, ig_7, ig_9, ig_47, ig_52, ig_54, ig_77, ig_82, ig_84, ig_152, \
                         ig_157, ig_159, ig_212, ig_217, ig_219, ig_317, ig_322, ig_324, \
                         ig_347, ig_352, ig_354, ig_377, ig_382, \
                         ig_384 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_77[k] = -f_215 * ig_2[k]
                  - f_215 * ig_7[k]
                  + f_288 * ig_9[k]
                  - f_215 * ig_47[k]
                  - f_215 * ig_52[k]
                  + f_288 * ig_54[k]
                  + f_128 * ig_77[k]
                  + f_128 * ig_82[k]
                  - f_133 * ig_84[k]
                  + f_215 * ig_152[k]
                  + f_215 * ig_157[k]
                  - f_288 * ig_159[k]
                  - f_128 * ig_212[k]
                  - f_128 * ig_217[k]
                  + f_133 * ig_219[k]
                  + f_215 * ig_317[k]
                  + f_215 * ig_322[k]
                  - f_288 * ig_324[k]
                  - f_128 * ig_347[k]
                  - f_128 * ig_352[k]
                  + f_133 * ig_354[k]
                  + f_128 * ig_377[k]
                  + f_128 * ig_382[k]
                  - f_133 * ig_384[k];
    }

#pragma omp simd aligned(ig_0, ig_5, ig_10, ig_12, ig_45, ig_50, ig_55, ig_57, ig_75, ig_80, \
                         ig_85, ig_87, ig_150, ig_155, ig_160, ig_162, ig_210, ig_215, ig_220, \
                         ig_222, ig_315, ig_320, ig_325, ig_327, ig_345, ig_350, ig_355, \
                         ig_357, ig_375, ig_380, ig_385, ig_387 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_78[k] = -f_291 * ig_0[k]
                  + f_151 * ig_5[k]
                  + f_291 * ig_10[k]
                  - f_151 * ig_12[k]
                  - f_291 * ig_45[k]
                  + f_151 * ig_50[k]
                  + f_291 * ig_55[k]
                  - f_151 * ig_57[k]
                  + f_152 * ig_75[k]
                  - f_119 * ig_80[k]
                  - f_152 * ig_85[k]
                  + f_119 * ig_87[k]
                  + f_291 * ig_150[k]
                  - f_151 * ig_155[k]
                  - f_291 * ig_160[k]
                  + f_151 * ig_162[k]
                  - f_152 * ig_210[k]
                  + f_119 * ig_215[k]
                  + f_152 * ig_220[k]
                  - f_119 * ig_222[k]
                  + f_291 * ig_315[k]
                  - f_151 * ig_320[k]
                  - f_291 * ig_325[k]
                  + f_151 * ig_327[k]
                  - f_152 * ig_345[k]
                  + f_119 * ig_350[k]
                  + f_152 * ig_355[k]
                  - f_119 * ig_357[k]
                  + f_152 * ig_375[k]
                  - f_119 * ig_380[k]
                  - f_152 * ig_385[k]
                  + f_119 * ig_387[k];
    }

#pragma omp simd aligned(ig_2, ig_7, ig_47, ig_52, ig_77, ig_82, ig_152, ig_157, ig_212, \
                         ig_217, ig_317, ig_322, ig_347, ig_352, ig_377, \
                         ig_382 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_79[k] = f_287 * ig_2[k]
                  - f_286 * ig_7[k]
                  + f_287 * ig_47[k]
                  - f_286 * ig_52[k]
                  - f_114 * ig_77[k]
                  + f_112 * ig_82[k]
                  - f_287 * ig_152[k]
                  + f_286 * ig_157[k]
                  + f_114 * ig_212[k]
                  - f_112 * ig_217[k]
                  - f_287 * ig_317[k]
                  + f_286 * ig_322[k]
                  + f_114 * ig_347[k]
                  - f_112 * ig_352[k]
                  - f_114 * ig_377[k]
                  + f_112 * ig_382[k];
    }

#pragma omp simd aligned(ig_0, ig_3, ig_10, ig_45, ig_48, ig_55, ig_75, ig_78, ig_85, ig_150, \
                         ig_153, ig_160, ig_210, ig_213, ig_220, ig_315, ig_318, ig_325, \
                         ig_345, ig_348, ig_355, ig_375, ig_378, \
                         ig_385 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_80[k] = f_292 * ig_0[k]
                  - f_157 * ig_3[k]
                  + f_292 * ig_10[k]
                  + f_292 * ig_45[k]
                  - f_157 * ig_48[k]
                  + f_292 * ig_55[k]
                  - f_159 * ig_75[k]
                  + f_160 * ig_78[k]
                  - f_159 * ig_85[k]
                  - f_292 * ig_150[k]
                  + f_157 * ig_153[k]
                  - f_292 * ig_160[k]
                  + f_159 * ig_210[k]
                  - f_160 * ig_213[k]
                  + f_159 * ig_220[k]
                  - f_292 * ig_315[k]
                  + f_157 * ig_318[k]
                  - f_292 * ig_325[k]
                  + f_159 * ig_345[k]
                  - f_160 * ig_348[k]
                  + f_159 * ig_355[k]
                  - f_159 * ig_375[k]
                  + f_160 * ig_378[k]
                  - f_159 * ig_385[k];
    }

#pragma omp simd aligned(ig_31, ig_36, ig_106, ig_111, ig_136, ig_141, ig_241, ig_246, ig_271, \
                         ig_276 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_81[k] = -f_105 * ig_31[k]
                  + f_105 * ig_36[k]
                  + f_103 * ig_106[k]
                  - f_103 * ig_111[k]
                  + f_106 * ig_136[k]
                  - f_106 * ig_141[k]
                  + f_102 * ig_241[k]
                  - f_102 * ig_246[k]
                  - f_104 * ig_271[k]
                  + f_104 * ig_276[k];
    }

#pragma omp simd aligned(ig_34, ig_41, ig_109, ig_116, ig_139, ig_146, ig_244, ig_251, ig_274, \
                         ig_281 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_82[k] = -f_108 * ig_34[k]
                  + f_113 * ig_41[k]
                  + f_109 * ig_109[k]
                  - f_110 * ig_116[k]
                  + f_112 * ig_139[k]
                  - f_114 * ig_146[k]
                  + f_107 * ig_244[k]
                  - f_108 * ig_251[k]
                  - f_111 * ig_274[k]
                  + f_112 * ig_281[k];
    }

#pragma omp simd aligned(ig_31, ig_36, ig_38, ig_106, ig_111, ig_113, ig_136, ig_141, ig_143, \
                         ig_241, ig_246, ig_248, ig_271, ig_276, \
                         ig_278 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_83[k] = f_121 * ig_31[k]
                  + f_121 * ig_36[k]
                  - f_122 * ig_38[k]
                  - f_117 * ig_106[k]
                  - f_117 * ig_111[k]
                  + f_118 * ig_113[k]
                  - f_123 * ig_136[k]
                  - f_123 * ig_141[k]
                  + f_124 * ig_143[k]
                  - f_115 * ig_241[k]
                  - f_115 * ig_246[k]
                  + f_116 * ig_248[k]
                  + f_119 * ig_271[k]
                  + f_119 * ig_276[k]
                  - f_120 * ig_278[k];
    }

#pragma omp simd aligned(ig_34, ig_41, ig_43, ig_109, ig_116, ig_118, ig_139, ig_146, ig_148, \
                         ig_244, ig_251, ig_253, ig_274, ig_281, \
                         ig_283 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_84[k] = f_131 * ig_34[k]
                  + f_131 * ig_41[k]
                  - f_132 * ig_43[k]
                  - f_127 * ig_109[k]
                  - f_127 * ig_116[k]
                  + f_128 * ig_118[k]
                  - f_128 * ig_139[k]
                  - f_128 * ig_146[k]
                  + f_133 * ig_148[k]
                  - f_125 * ig_244[k]
                  - f_125 * ig_251[k]
                  + f_126 * ig_253[k]
                  + f_129 * ig_274[k]
                  + f_129 * ig_281[k]
                  - f_130 * ig_283[k];
    }

#pragma omp simd aligned(ig_30, ig_33, ig_35, ig_40, ig_42, ig_44, ig_105, ig_108, ig_110, \
                         ig_115, ig_117, ig_119, ig_135, ig_138, ig_140, ig_145, ig_147, \
                         ig_149, ig_240, ig_243, ig_245, ig_250, ig_252, ig_254, ig_270, \
                         ig_273, ig_275, ig_280, ig_282, ig_284 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_85[k] = -f_144 * ig_30[k]
                  - f_138 * ig_33[k]
                  + f_137 * ig_35[k]
                  - f_144 * ig_40[k]
                  + f_137 * ig_42[k]
                  - f_145 * ig_44[k]
                  + f_138 * ig_105[k]
                  + f_139 * ig_108[k]
                  - f_140 * ig_110[k]
                  + f_138 * ig_115[k]
                  - f_140 * ig_117[k]
                  + f_141 * ig_119[k]
                  + f_145 * ig_135[k]
                  + f_141 * ig_138[k]
                  - f_143 * ig_140[k]
                  + f_145 * ig_145[k]
                  - f_143 * ig_147[k]
                  + f_146 * ig_149[k]
                  + f_134 * ig_240[k]
                  + f_135 * ig_243[k]
                  - f_136 * ig_245[k]
                  + f_134 * ig_250[k]
                  - f_136 * ig_252[k]
                  + f_137 * ig_254[k]
                  - f_137 * ig_270[k]
                  - f_140 * ig_273[k]
                  + f_142 * ig_275[k]
                  - f_137 * ig_280[k]
                  + f_142 * ig_282[k]
                  - f_143 * ig_284[k];
    }

#pragma omp simd aligned(ig_32, ig_37, ig_39, ig_107, ig_112, ig_114, ig_137, ig_142, ig_144, \
                         ig_242, ig_247, ig_249, ig_272, ig_277, \
                         ig_279 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_86[k] = f_131 * ig_32[k]
                  + f_131 * ig_37[k]
                  - f_132 * ig_39[k]
                  - f_127 * ig_107[k]
                  - f_127 * ig_112[k]
                  + f_128 * ig_114[k]
                  - f_128 * ig_137[k]
                  - f_128 * ig_142[k]
                  + f_133 * ig_144[k]
                  - f_125 * ig_242[k]
                  - f_125 * ig_247[k]
                  + f_126 * ig_249[k]
                  + f_129 * ig_272[k]
                  + f_129 * ig_277[k]
                  - f_130 * ig_279[k];
    }

#pragma omp simd aligned(ig_30, ig_35, ig_40, ig_42, ig_105, ig_110, ig_115, ig_117, ig_135, \
                         ig_140, ig_145, ig_147, ig_240, ig_245, ig_250, ig_252, ig_270, \
                         ig_275, ig_280, ig_282 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_87[k] = f_151 * ig_30[k]
                  - f_115 * ig_35[k]
                  - f_151 * ig_40[k]
                  + f_115 * ig_42[k]
                  - f_121 * ig_105[k]
                  + f_122 * ig_110[k]
                  + f_121 * ig_115[k]
                  - f_122 * ig_117[k]
                  - f_152 * ig_135[k]
                  + f_119 * ig_140[k]
                  + f_152 * ig_145[k]
                  - f_119 * ig_147[k]
                  - f_147 * ig_240[k]
                  + f_148 * ig_245[k]
                  + f_147 * ig_250[k]
                  - f_148 * ig_252[k]
                  + f_149 * ig_270[k]
                  - f_150 * ig_275[k]
                  - f_149 * ig_280[k]
                  + f_150 * ig_282[k];
    }

#pragma omp simd aligned(ig_32, ig_37, ig_107, ig_112, ig_137, ig_142, ig_242, ig_247, ig_272, \
                         ig_277 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_88[k] = -f_113 * ig_32[k]
                  + f_108 * ig_37[k]
                  + f_110 * ig_107[k]
                  - f_109 * ig_112[k]
                  + f_114 * ig_137[k]
                  - f_112 * ig_142[k]
                  + f_108 * ig_242[k]
                  - f_107 * ig_247[k]
                  - f_112 * ig_272[k]
                  + f_111 * ig_277[k];
    }

#pragma omp simd aligned(ig_30, ig_33, ig_40, ig_105, ig_108, ig_115, ig_135, ig_138, ig_145, \
                         ig_240, ig_243, ig_250, ig_270, ig_273, \
                         ig_280 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_89[k] = -f_157 * ig_30[k]
                  + f_158 * ig_33[k]
                  - f_157 * ig_40[k]
                  + f_155 * ig_105[k]
                  - f_102 * ig_108[k]
                  + f_155 * ig_115[k]
                  + f_159 * ig_135[k]
                  - f_160 * ig_138[k]
                  + f_159 * ig_145[k]
                  + f_153 * ig_240[k]
                  - f_154 * ig_243[k]
                  + f_153 * ig_250[k]
                  - f_103 * ig_270[k]
                  + f_156 * ig_273[k]
                  - f_103 * ig_280[k];
    }

#pragma omp simd aligned(ig_1, ig_6, ig_46, ig_51, ig_76, ig_81, ig_151, ig_156, ig_181, \
                         ig_186, ig_316, ig_321, ig_346, ig_351 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_90[k] = -f_98 * ig_1[k]
                  + f_98 * ig_6[k]
                  + f_293 * ig_46[k]
                  - f_293 * ig_51[k]
                  + f_100 * ig_76[k]
                  - f_100 * ig_81[k]
                  + f_293 * ig_151[k]
                  - f_293 * ig_156[k]
                  - f_101 * ig_181[k]
                  + f_101 * ig_186[k]
                  - f_98 * ig_316[k]
                  + f_98 * ig_321[k]
                  + f_100 * ig_346[k]
                  - f_100 * ig_351[k];
    }

#pragma omp simd aligned(ig_4, ig_11, ig_49, ig_56, ig_79, ig_86, ig_154, ig_161, ig_184, \
                         ig_191, ig_319, ig_326, ig_349, ig_356 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_91[k] = -f_294 * ig_4[k]
                  + f_295 * ig_11[k]
                  + f_296 * ig_49[k]
                  - f_297 * ig_56[k]
                  + f_298 * ig_79[k]
                  - f_299 * ig_86[k]
                  + f_296 * ig_154[k]
                  - f_297 * ig_161[k]
                  - f_300 * ig_184[k]
                  + f_301 * ig_191[k]
                  - f_294 * ig_319[k]
                  + f_295 * ig_326[k]
                  + f_298 * ig_349[k]
                  - f_299 * ig_356[k];
    }

#pragma omp simd aligned(ig_1, ig_6, ig_8, ig_46, ig_51, ig_53, ig_76, ig_81, ig_83, ig_151, \
                         ig_156, ig_158, ig_181, ig_186, ig_188, ig_316, ig_321, ig_323, \
                         ig_346, ig_351, ig_353 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_92[k] = f_302 * ig_1[k]
                  + f_302 * ig_6[k]
                  - f_303 * ig_8[k]
                  - f_233 * ig_46[k]
                  - f_233 * ig_51[k]
                  + f_234 * ig_53[k]
                  - f_283 * ig_76[k]
                  - f_283 * ig_81[k]
                  + f_235 * ig_83[k]
                  - f_233 * ig_151[k]
                  - f_233 * ig_156[k]
                  + f_234 * ig_158[k]
                  + f_235 * ig_181[k]
                  + f_235 * ig_186[k]
                  - f_304 * ig_188[k]
                  + f_302 * ig_316[k]
                  + f_302 * ig_321[k]
                  - f_303 * ig_323[k]
                  - f_283 * ig_346[k]
                  - f_283 * ig_351[k]
                  + f_235 * ig_353[k];
    }

#pragma omp simd aligned(ig_4, ig_11, ig_13, ig_49, ig_56, ig_58, ig_79, ig_86, ig_88, ig_154, \
                         ig_161, ig_163, ig_184, ig_191, ig_193, ig_319, ig_326, ig_328, \
                         ig_349, ig_356, ig_358 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_93[k] = f_305 * ig_4[k]
                  + f_305 * ig_11[k]
                  - f_306 * ig_13[k]
                  - f_239 * ig_49[k]
                  - f_239 * ig_56[k]
                  + f_307 * ig_58[k]
                  - f_241 * ig_79[k]
                  - f_241 * ig_86[k]
                  + f_244 * ig_88[k]
                  - f_239 * ig_154[k]
                  - f_239 * ig_161[k]
                  + f_307 * ig_163[k]
                  + f_242 * ig_184[k]
                  + f_242 * ig_191[k]
                  - f_308 * ig_193[k]
                  + f_305 * ig_319[k]
                  + f_305 * ig_326[k]
                  - f_306 * ig_328[k]
                  - f_241 * ig_349[k]
                  - f_241 * ig_356[k]
                  + f_244 * ig_358[k];
    }

#pragma omp simd aligned(ig_0, ig_3, ig_5, ig_10, ig_12, ig_14, ig_45, ig_48, ig_50, ig_55, \
                         ig_57, ig_59, ig_75, ig_78, ig_80, ig_85, ig_87, ig_89, ig_150, \
                         ig_153, ig_155, ig_160, ig_162, ig_164, ig_180, ig_183, ig_185, \
                         ig_190, ig_192, ig_194, ig_315, ig_318, ig_320, ig_325, ig_327, \
                         ig_329, ig_345, ig_348, ig_350, ig_355, ig_357, \
                         ig_359 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_94[k] = -f_309 * ig_0[k]
                  - f_310 * ig_3[k]
                  + f_87 * ig_5[k]
                  - f_309 * ig_10[k]
                  + f_87 * ig_12[k]
                  - f_311 * ig_14[k]
                  + f_312 * ig_45[k]
                  + f_313 * ig_48[k]
                  - f_90 * ig_50[k]
                  + f_312 * ig_55[k]
                  - f_90 * ig_57[k]
                  + f_314 * ig_59[k]
                  + f_313 * ig_75[k]
                  + f_315 * ig_78[k]
                  - f_91 * ig_80[k]
                  + f_313 * ig_85[k]
                  - f_91 * ig_87[k]
                  + f_316 * ig_89[k]
                  + f_312 * ig_150[k]
                  + f_313 * ig_153[k]
                  - f_90 * ig_155[k]
                  + f_312 * ig_160[k]
                  - f_90 * ig_162[k]
                  + f_314 * ig_164[k]
                  - f_317 * ig_180[k]
                  - f_318 * ig_183[k]
                  + f_319 * ig_185[k]
                  - f_317 * ig_190[k]
                  + f_319 * ig_192[k]
                  - f_320 * ig_194[k]
                  - f_309 * ig_315[k]
                  - f_310 * ig_318[k]
                  + f_87 * ig_320[k]
                  - f_309 * ig_325[k]
                  + f_87 * ig_327[k]
                  - f_311 * ig_329[k]
                  + f_313 * ig_345[k]
                  + f_315 * ig_348[k]
                  - f_91 * ig_350[k]
                  + f_313 * ig_355[k]
                  - f_91 * ig_357[k]
                  + f_316 * ig_359[k];
    }

#pragma omp simd aligned(ig_2, ig_7, ig_9, ig_47, ig_52, ig_54, ig_77, ig_82, ig_84, ig_152, \
                         ig_157, ig_159, ig_182, ig_187, ig_189, ig_317, ig_322, ig_324, \
                         ig_347, ig_352, ig_354 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_95[k] = f_305 * ig_2[k]
                  + f_305 * ig_7[k]
                  - f_306 * ig_9[k]
                  - f_239 * ig_47[k]
                  - f_239 * ig_52[k]
                  + f_307 * ig_54[k]
                  - f_241 * ig_77[k]
                  - f_241 * ig_82[k]
                  + f_244 * ig_84[k]
                  - f_239 * ig_152[k]
                  - f_239 * ig_157[k]
                  + f_307 * ig_159[k]
                  + f_242 * ig_182[k]
                  + f_242 * ig_187[k]
                  - f_308 * ig_189[k]
                  + f_305 * ig_317[k]
                  + f_305 * ig_322[k]
                  - f_306 * ig_324[k]
                  - f_241 * ig_347[k]
                  - f_241 * ig_352[k]
                  + f_244 * ig_354[k];
    }

#pragma omp simd aligned(ig_0, ig_5, ig_10, ig_12, ig_45, ig_50, ig_55, ig_57, ig_75, ig_80, \
                         ig_85, ig_87, ig_150, ig_155, ig_160, ig_162, ig_180, ig_185, ig_190, \
                         ig_192, ig_315, ig_320, ig_325, ig_327, ig_345, ig_350, ig_355, \
                         ig_357 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_96[k] = f_321 * ig_0[k]
                  - f_322 * ig_5[k]
                  - f_321 * ig_10[k]
                  + f_322 * ig_12[k]
                  - f_277 * ig_45[k]
                  + f_281 * ig_50[k]
                  + f_277 * ig_55[k]
                  - f_281 * ig_57[k]
                  - f_233 * ig_75[k]
                  + f_234 * ig_80[k]
                  + f_233 * ig_85[k]
                  - f_234 * ig_87[k]
                  - f_277 * ig_150[k]
                  + f_281 * ig_155[k]
                  + f_277 * ig_160[k]
                  - f_281 * ig_162[k]
                  + f_234 * ig_180[k]
                  - f_323 * ig_185[k]
                  - f_234 * ig_190[k]
                  + f_323 * ig_192[k]
                  + f_321 * ig_315[k]
                  - f_322 * ig_320[k]
                  - f_321 * ig_325[k]
                  + f_322 * ig_327[k]
                  - f_233 * ig_345[k]
                  + f_234 * ig_350[k]
                  + f_233 * ig_355[k]
                  - f_234 * ig_357[k];
    }

#pragma omp simd aligned(ig_2, ig_7, ig_47, ig_52, ig_77, ig_82, ig_152, ig_157, ig_182, \
                         ig_187, ig_317, ig_322, ig_347, ig_352 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_97[k] = -f_295 * ig_2[k]
                  + f_294 * ig_7[k]
                  + f_297 * ig_47[k]
                  - f_296 * ig_52[k]
                  + f_299 * ig_77[k]
                  - f_298 * ig_82[k]
                  + f_297 * ig_152[k]
                  - f_296 * ig_157[k]
                  - f_301 * ig_182[k]
                  + f_300 * ig_187[k]
                  - f_295 * ig_317[k]
                  + f_294 * ig_322[k]
                  + f_299 * ig_347[k]
                  - f_298 * ig_352[k];
    }

#pragma omp simd aligned(ig_0, ig_3, ig_10, ig_45, ig_48, ig_55, ig_75, ig_78, ig_85, ig_150, \
                         ig_153, ig_160, ig_180, ig_183, ig_190, ig_315, ig_318, ig_325, \
                         ig_345, ig_348, ig_355 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_98[k] = -f_324 * ig_0[k]
                  + f_325 * ig_3[k]
                  - f_324 * ig_10[k]
                  + f_326 * ig_45[k]
                  - f_327 * ig_48[k]
                  + f_326 * ig_55[k]
                  + f_328 * ig_75[k]
                  - f_329 * ig_78[k]
                  + f_328 * ig_85[k]
                  + f_326 * ig_150[k]
                  - f_327 * ig_153[k]
                  + f_326 * ig_160[k]
                  - f_329 * ig_180[k]
                  + f_330 * ig_183[k]
                  - f_329 * ig_190[k]
                  - f_324 * ig_315[k]
                  + f_325 * ig_318[k]
                  - f_324 * ig_325[k]
                  + f_328 * ig_345[k]
                  - f_329 * ig_348[k]
                  + f_328 * ig_355[k];
    }

#pragma omp simd aligned(ig_31, ig_34, ig_36, ig_41, ig_106, ig_109, ig_111, ig_116, ig_241, \
                         ig_244, ig_246, ig_251 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_99[k] = f_32 * ig_31[k]
                  - f_32 * ig_36[k]
                  - f_31 * ig_106[k]
                  + f_31 * ig_111[k]
                  + f_30 * ig_241[k]
                  - f_30 * ig_246[k];

        g_100[k] = f_37 * ig_34[k]
                   - f_38 * ig_41[k]
                   - f_35 * ig_109[k]
                   + f_36 * ig_116[k]
                   + f_33 * ig_244[k]
                   - f_34 * ig_251[k];
    }

#pragma omp simd aligned(ig_31, ig_36, ig_38, ig_106, ig_111, ig_113, ig_241, ig_246, \
                         ig_248 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_101[k] = -f_43 * ig_31[k]
                   - f_43 * ig_36[k]
                   + f_44 * ig_38[k]
                   + f_41 * ig_106[k]
                   + f_41 * ig_111[k]
                   - f_42 * ig_113[k]
                   - f_39 * ig_241[k]
                   - f_39 * ig_246[k]
                   + f_40 * ig_248[k];
    }

#pragma omp simd aligned(ig_34, ig_41, ig_43, ig_109, ig_116, ig_118, ig_244, ig_251, \
                         ig_253 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_102[k] = -f_49 * ig_34[k]
                   - f_49 * ig_41[k]
                   + f_50 * ig_43[k]
                   + f_47 * ig_109[k]
                   + f_47 * ig_116[k]
                   - f_48 * ig_118[k]
                   - f_45 * ig_244[k]
                   - f_45 * ig_251[k]
                   + f_46 * ig_253[k];
    }

#pragma omp simd aligned(ig_30, ig_33, ig_35, ig_40, ig_42, ig_44, ig_105, ig_108, ig_110, \
                         ig_115, ig_117, ig_119, ig_240, ig_243, ig_245, ig_250, ig_252, \
                         ig_254 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_103[k] = f_58 * ig_30[k]
                   + f_59 * ig_33[k]
                   - f_60 * ig_35[k]
                   + f_58 * ig_40[k]
                   - f_60 * ig_42[k]
                   + f_61 * ig_44[k]
                   - f_52 * ig_105[k]
                   - f_55 * ig_108[k]
                   + f_56 * ig_110[k]
                   - f_52 * ig_115[k]
                   + f_56 * ig_117[k]
                   - f_57 * ig_119[k]
                   + f_51 * ig_240[k]
                   + f_52 * ig_243[k]
                   - f_53 * ig_245[k]
                   + f_51 * ig_250[k]
                   - f_53 * ig_252[k]
                   + f_54 * ig_254[k];
    }

#pragma omp simd aligned(ig_32, ig_37, ig_39, ig_107, ig_112, ig_114, ig_242, ig_247, \
                         ig_249 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_104[k] = -f_49 * ig_32[k]
                   - f_49 * ig_37[k]
                   + f_50 * ig_39[k]
                   + f_47 * ig_107[k]
                   + f_47 * ig_112[k]
                   - f_48 * ig_114[k]
                   - f_45 * ig_242[k]
                   - f_45 * ig_247[k]
                   + f_46 * ig_249[k];
    }

#pragma omp simd aligned(ig_30, ig_35, ig_40, ig_42, ig_105, ig_110, ig_115, ig_117, ig_240, \
                         ig_245, ig_250, ig_252 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_105[k] = -f_64 * ig_30[k]
                   + f_65 * ig_35[k]
                   + f_64 * ig_40[k]
                   - f_65 * ig_42[k]
                   + f_39 * ig_105[k]
                   - f_40 * ig_110[k]
                   - f_39 * ig_115[k]
                   + f_40 * ig_117[k]
                   - f_62 * ig_240[k]
                   + f_63 * ig_245[k]
                   + f_62 * ig_250[k]
                   - f_63 * ig_252[k];
    }

#pragma omp simd aligned(ig_32, ig_37, ig_107, ig_112, ig_242, ig_247 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_106[k] = f_38 * ig_32[k]
                   - f_37 * ig_37[k]
                   - f_36 * ig_107[k]
                   + f_35 * ig_112[k]
                   + f_34 * ig_242[k]
                   - f_33 * ig_247[k];
    }

#pragma omp simd aligned(ig_30, ig_33, ig_40, ig_105, ig_108, ig_115, ig_240, ig_243, \
                         ig_250 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_107[k] = f_70 * ig_30[k]
                   - f_71 * ig_33[k]
                   + f_70 * ig_40[k]
                   - f_68 * ig_105[k]
                   + f_69 * ig_108[k]
                   - f_68 * ig_115[k]
                   + f_66 * ig_240[k]
                   - f_67 * ig_243[k]
                   + f_66 * ig_250[k];
    }

#pragma omp simd aligned(ig_1, ig_6, ig_46, ig_51, ig_151, ig_156, ig_316, \
                         ig_321 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_108[k] = f_331 * ig_1[k]
                   - f_331 * ig_6[k]
                   - f_332 * ig_46[k]
                   + f_332 * ig_51[k]
                   + f_332 * ig_151[k]
                   - f_332 * ig_156[k]
                   - f_331 * ig_316[k]
                   + f_331 * ig_321[k];
    }

#pragma omp simd aligned(ig_4, ig_11, ig_49, ig_56, ig_154, ig_161, ig_319, \
                         ig_326 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_109[k] = f_333 * ig_4[k]
                   - f_334 * ig_11[k]
                   - f_335 * ig_49[k]
                   + f_336 * ig_56[k]
                   + f_335 * ig_154[k]
                   - f_336 * ig_161[k]
                   - f_333 * ig_319[k]
                   + f_334 * ig_326[k];
    }

#pragma omp simd aligned(ig_1, ig_6, ig_8, ig_46, ig_51, ig_53, ig_151, ig_156, ig_158, \
                         ig_316, ig_321, ig_323 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_110[k] = -f_337 * ig_1[k]
                   - f_337 * ig_6[k]
                   + f_6 * ig_8[k]
                   + f_338 * ig_46[k]
                   + f_338 * ig_51[k]
                   - f_339 * ig_53[k]
                   - f_338 * ig_151[k]
                   - f_338 * ig_156[k]
                   + f_339 * ig_158[k]
                   + f_337 * ig_316[k]
                   + f_337 * ig_321[k]
                   - f_6 * ig_323[k];
    }

#pragma omp simd aligned(ig_4, ig_11, ig_13, ig_49, ig_56, ig_58, ig_154, ig_161, ig_163, \
                         ig_319, ig_326, ig_328 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_111[k] = -f_340 * ig_4[k]
                   - f_340 * ig_11[k]
                   + f_341 * ig_13[k]
                   + f_342 * ig_49[k]
                   + f_342 * ig_56[k]
                   - f_12 * ig_58[k]
                   - f_342 * ig_154[k]
                   - f_342 * ig_161[k]
                   + f_12 * ig_163[k]
                   + f_340 * ig_319[k]
                   + f_340 * ig_326[k]
                   - f_341 * ig_328[k];
    }

#pragma omp simd aligned(ig_0, ig_3, ig_5, ig_10, ig_12, ig_14, ig_45, ig_48, ig_50, ig_55, \
                         ig_57, ig_59, ig_150, ig_153, ig_155, ig_160, ig_162, ig_164, ig_315, \
                         ig_318, ig_320, ig_325, ig_327, ig_329 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_112[k] = f_343 * ig_0[k]
                   + f_344 * ig_3[k]
                   - f_345 * ig_5[k]
                   + f_343 * ig_10[k]
                   - f_345 * ig_12[k]
                   + f_346 * ig_14[k]
                   - f_347 * ig_45[k]
                   - f_348 * ig_48[k]
                   + f_349 * ig_50[k]
                   - f_347 * ig_55[k]
                   + f_349 * ig_57[k]
                   - f_19 * ig_59[k]
                   + f_347 * ig_150[k]
                   + f_348 * ig_153[k]
                   - f_349 * ig_155[k]
                   + f_347 * ig_160[k]
                   - f_349 * ig_162[k]
                   + f_19 * ig_164[k]
                   - f_343 * ig_315[k]
                   - f_344 * ig_318[k]
                   + f_345 * ig_320[k]
                   - f_343 * ig_325[k]
                   + f_345 * ig_327[k]
                   - f_346 * ig_329[k];
    }

#pragma omp simd aligned(ig_2, ig_7, ig_9, ig_47, ig_52, ig_54, ig_152, ig_157, ig_159, \
                         ig_317, ig_322, ig_324 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_113[k] = -f_340 * ig_2[k]
                   - f_340 * ig_7[k]
                   + f_341 * ig_9[k]
                   + f_342 * ig_47[k]
                   + f_342 * ig_52[k]
                   - f_12 * ig_54[k]
                   - f_342 * ig_152[k]
                   - f_342 * ig_157[k]
                   + f_12 * ig_159[k]
                   + f_340 * ig_317[k]
                   + f_340 * ig_322[k]
                   - f_341 * ig_324[k];
    }

#pragma omp simd aligned(ig_0, ig_5, ig_10, ig_12, ig_45, ig_50, ig_55, ig_57, ig_150, ig_155, \
                         ig_160, ig_162, ig_315, ig_320, ig_325, \
                         ig_327 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_114[k] = -f_350 * ig_0[k]
                   + f_22 * ig_5[k]
                   + f_350 * ig_10[k]
                   - f_22 * ig_12[k]
                   + f_351 * ig_45[k]
                   - f_352 * ig_50[k]
                   - f_351 * ig_55[k]
                   + f_352 * ig_57[k]
                   - f_351 * ig_150[k]
                   + f_352 * ig_155[k]
                   + f_351 * ig_160[k]
                   - f_352 * ig_162[k]
                   + f_350 * ig_315[k]
                   - f_22 * ig_320[k]
                   - f_350 * ig_325[k]
                   + f_22 * ig_327[k];
    }

#pragma omp simd aligned(ig_2, ig_7, ig_47, ig_52, ig_152, ig_157, ig_317, \
                         ig_322 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_115[k] = f_334 * ig_2[k]
                   - f_333 * ig_7[k]
                   - f_336 * ig_47[k]
                   + f_335 * ig_52[k]
                   + f_336 * ig_152[k]
                   - f_335 * ig_157[k]
                   - f_334 * ig_317[k]
                   + f_333 * ig_322[k];
    }

#pragma omp simd aligned(ig_0, ig_3, ig_10, ig_45, ig_48, ig_55, ig_150, ig_153, ig_160, \
                         ig_315, ig_318, ig_325 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_116[k] = f_353 * ig_0[k]
                   - f_26 * ig_3[k]
                   + f_353 * ig_10[k]
                   - f_354 * ig_45[k]
                   + f_355 * ig_48[k]
                   - f_354 * ig_55[k]
                   + f_354 * ig_150[k]
                   - f_355 * ig_153[k]
                   + f_354 * ig_160[k]
                   - f_353 * ig_315[k]
                   + f_26 * ig_318[k]
                   - f_353 * ig_325[k];
    }
}

}  // namespace simdtrf
