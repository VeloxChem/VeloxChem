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


#include "SimdTransformGL.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_gl(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t gl,
             const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 0.46875 * std::sqrt(1001.0);
    const auto f_1 = 3.28125 * std::sqrt(1001.0);
    const auto f_2 = 1.640625 * std::sqrt(1001.0);
    const auto f_3 = 8.203125 * std::sqrt(1001.0);
    const auto f_4 = 4.921875 * std::sqrt(1001.0);
    const auto f_5 = 0.234375 * std::sqrt(1001.0);
    const auto f_6 = 0.046875 * std::sqrt(30030.0);
    const auto f_7 = 0.109375 * std::sqrt(30030.0);
    const auto f_8 = 0.65625 * std::sqrt(30030.0);
    const auto f_9 = 2.1875 * std::sqrt(30030.0);
    const auto f_10 = 1.640625 * std::sqrt(715.0);
    const auto f_11 = 6.5625 * std::sqrt(715.0);
    const auto f_12 = 2.953125 * std::sqrt(715.0);
    const auto f_13 = 13.125 * std::sqrt(715.0);
    const auto f_14 = 0.328125 * std::sqrt(715.0);
    const auto f_15 = 1.3125 * std::sqrt(715.0);
    const auto f_16 = 0.65625 * std::sqrt(55.0);
    const auto f_17 = 15.75 * std::sqrt(55.0);
    const auto f_18 = 26.25 * std::sqrt(55.0);
    const auto f_19 = 4.921875 * std::sqrt(33.0);
    const auto f_20 = 8.203125 * std::sqrt(33.0);
    const auto f_21 = 32.8125 * std::sqrt(33.0);
    const auto f_22 = 1.640625 * std::sqrt(33.0);
    const auto f_23 = 21.875 * std::sqrt(33.0);
    const auto f_24 = 26.25 * std::sqrt(33.0);
    const auto f_25 = 10.9375 * std::sqrt(33.0);
    const auto f_26 = 8.75 * std::sqrt(33.0);
    const auto f_27 = 1.640625 * std::sqrt(2.0);
    const auto f_28 = 4.921875 * std::sqrt(2.0);
    const auto f_29 = 49.21875 * std::sqrt(2.0);
    const auto f_30 = 98.4375 * std::sqrt(2.0);
    const auto f_31 = 131.25 * std::sqrt(2.0);
    const auto f_32 = 52.5 * std::sqrt(2.0);
    const auto f_33 = 1.640625 * std::sqrt(35.0);
    const auto f_34 = 4.921875 * std::sqrt(35.0);
    const auto f_35 = 13.125 * std::sqrt(35.0);
    const auto f_36 = 26.25 * std::sqrt(35.0);
    const auto f_37 = 15.75 * std::sqrt(35.0);
    const auto f_38 = 3.0 * std::sqrt(35.0);
    const auto f_39 = 0.13671875 * std::sqrt(35.0);
    const auto f_40 = 0.546875 * std::sqrt(35.0);
    const auto f_41 = 4.375 * std::sqrt(35.0);
    const auto f_42 = 0.8203125 * std::sqrt(35.0);
    const auto f_43 = 7.0 * std::sqrt(35.0);
    const auto f_44 = 0.5 * std::sqrt(35.0);
    const auto f_45 = 0.8203125 * std::sqrt(2.0);
    const auto f_46 = 24.609375 * std::sqrt(2.0);
    const auto f_47 = 65.625 * std::sqrt(2.0);
    const auto f_48 = 26.25 * std::sqrt(2.0);
    const auto f_49 = 0.1640625 * std::sqrt(55.0);
    const auto f_50 = 3.9375 * std::sqrt(55.0);
    const auto f_51 = 1.640625 * std::sqrt(55.0);
    const auto f_52 = 19.6875 * std::sqrt(55.0);
    const auto f_53 = 6.5625 * std::sqrt(55.0);
    const auto f_54 = 39.375 * std::sqrt(55.0);
    const auto f_55 = 0.0078125 * std::sqrt(30030.0);
    const auto f_56 = 1.640625 * std::sqrt(30030.0);
    const auto f_57 = 0.05859375 * std::sqrt(1001.0);
    const auto f_58 = 4.1015625 * std::sqrt(1001.0);
    const auto f_59 = 0.703125 * std::sqrt(2002.0);
    const auto f_60 = 4.921875 * std::sqrt(2002.0);
    const auto f_61 = 0.234375 * std::sqrt(2002.0);
    const auto f_62 = 1.640625 * std::sqrt(2002.0);
    const auto f_63 = 2.4609375 * std::sqrt(2002.0);
    const auto f_64 = 12.3046875 * std::sqrt(2002.0);
    const auto f_65 = 7.3828125 * std::sqrt(2002.0);
    const auto f_66 = 0.3515625 * std::sqrt(2002.0);
    const auto f_67 = 0.8203125 * std::sqrt(2002.0);
    const auto f_68 = 4.1015625 * std::sqrt(2002.0);
    const auto f_69 = 0.1171875 * std::sqrt(2002.0);
    const auto f_70 = 0.140625 * std::sqrt(15015.0);
    const auto f_71 = 0.328125 * std::sqrt(15015.0);
    const auto f_72 = 1.96875 * std::sqrt(15015.0);
    const auto f_73 = 6.5625 * std::sqrt(15015.0);
    const auto f_74 = 0.046875 * std::sqrt(15015.0);
    const auto f_75 = 0.109375 * std::sqrt(15015.0);
    const auto f_76 = 0.65625 * std::sqrt(15015.0);
    const auto f_77 = 2.1875 * std::sqrt(15015.0);
    const auto f_78 = 2.4609375 * std::sqrt(1430.0);
    const auto f_79 = 9.84375 * std::sqrt(1430.0);
    const auto f_80 = 4.4296875 * std::sqrt(1430.0);
    const auto f_81 = 19.6875 * std::sqrt(1430.0);
    const auto f_82 = 0.4921875 * std::sqrt(1430.0);
    const auto f_83 = 1.96875 * std::sqrt(1430.0);
    const auto f_84 = 0.8203125 * std::sqrt(1430.0);
    const auto f_85 = 3.28125 * std::sqrt(1430.0);
    const auto f_86 = 1.4765625 * std::sqrt(1430.0);
    const auto f_87 = 6.5625 * std::sqrt(1430.0);
    const auto f_88 = 0.1640625 * std::sqrt(1430.0);
    const auto f_89 = 0.65625 * std::sqrt(1430.0);
    const auto f_90 = 0.984375 * std::sqrt(110.0);
    const auto f_91 = 23.625 * std::sqrt(110.0);
    const auto f_92 = 39.375 * std::sqrt(110.0);
    const auto f_93 = 0.328125 * std::sqrt(110.0);
    const auto f_94 = 7.875 * std::sqrt(110.0);
    const auto f_95 = 13.125 * std::sqrt(110.0);
    const auto f_96 = 7.3828125 * std::sqrt(66.0);
    const auto f_97 = 12.3046875 * std::sqrt(66.0);
    const auto f_98 = 49.21875 * std::sqrt(66.0);
    const auto f_99 = 2.4609375 * std::sqrt(66.0);
    const auto f_100 = 32.8125 * std::sqrt(66.0);
    const auto f_101 = 39.375 * std::sqrt(66.0);
    const auto f_102 = 16.40625 * std::sqrt(66.0);
    const auto f_103 = 13.125 * std::sqrt(66.0);
    const auto f_104 = 4.1015625 * std::sqrt(66.0);
    const auto f_105 = 0.8203125 * std::sqrt(66.0);
    const auto f_106 = 10.9375 * std::sqrt(66.0);
    const auto f_107 = 5.46875 * std::sqrt(66.0);
    const auto f_108 = 4.375 * std::sqrt(66.0);
    const auto f_109 = 2.4609375 * std::sqrt(70.0);
    const auto f_110 = 7.3828125 * std::sqrt(70.0);
    const auto f_111 = 19.6875 * std::sqrt(70.0);
    const auto f_112 = 39.375 * std::sqrt(70.0);
    const auto f_113 = 23.625 * std::sqrt(70.0);
    const auto f_114 = 4.5 * std::sqrt(70.0);
    const auto f_115 = 0.8203125 * std::sqrt(70.0);
    const auto f_116 = 6.5625 * std::sqrt(70.0);
    const auto f_117 = 13.125 * std::sqrt(70.0);
    const auto f_118 = 7.875 * std::sqrt(70.0);
    const auto f_119 = 1.5 * std::sqrt(70.0);
    const auto f_120 = 0.205078125 * std::sqrt(70.0);
    const auto f_121 = 1.23046875 * std::sqrt(70.0);
    const auto f_122 = 10.5 * std::sqrt(70.0);
    const auto f_123 = 0.75 * std::sqrt(70.0);
    const auto f_124 = 0.068359375 * std::sqrt(70.0);
    const auto f_125 = 0.2734375 * std::sqrt(70.0);
    const auto f_126 = 2.1875 * std::sqrt(70.0);
    const auto f_127 = 0.41015625 * std::sqrt(70.0);
    const auto f_128 = 3.5 * std::sqrt(70.0);
    const auto f_129 = 0.25 * std::sqrt(70.0);
    const auto f_130 = 0.24609375 * std::sqrt(110.0);
    const auto f_131 = 5.90625 * std::sqrt(110.0);
    const auto f_132 = 2.4609375 * std::sqrt(110.0);
    const auto f_133 = 29.53125 * std::sqrt(110.0);
    const auto f_134 = 9.84375 * std::sqrt(110.0);
    const auto f_135 = 59.0625 * std::sqrt(110.0);
    const auto f_136 = 0.08203125 * std::sqrt(110.0);
    const auto f_137 = 1.96875 * std::sqrt(110.0);
    const auto f_138 = 0.8203125 * std::sqrt(110.0);
    const auto f_139 = 3.28125 * std::sqrt(110.0);
    const auto f_140 = 19.6875 * std::sqrt(110.0);
    const auto f_141 = 0.0234375 * std::sqrt(15015.0);
    const auto f_142 = 4.921875 * std::sqrt(15015.0);
    const auto f_143 = 0.0078125 * std::sqrt(15015.0);
    const auto f_144 = 1.640625 * std::sqrt(15015.0);
    const auto f_145 = 0.087890625 * std::sqrt(2002.0);
    const auto f_146 = 6.15234375 * std::sqrt(2002.0);
    const auto f_147 = 0.029296875 * std::sqrt(2002.0);
    const auto f_148 = 2.05078125 * std::sqrt(2002.0);
    const auto f_149 = 0.46875 * std::sqrt(143.0);
    const auto f_150 = 3.28125 * std::sqrt(143.0);
    const auto f_151 = 2.8125 * std::sqrt(143.0);
    const auto f_152 = 19.6875 * std::sqrt(143.0);
    const auto f_153 = 1.640625 * std::sqrt(143.0);
    const auto f_154 = 8.203125 * std::sqrt(143.0);
    const auto f_155 = 4.921875 * std::sqrt(143.0);
    const auto f_156 = 0.234375 * std::sqrt(143.0);
    const auto f_157 = 9.84375 * std::sqrt(143.0);
    const auto f_158 = 49.21875 * std::sqrt(143.0);
    const auto f_159 = 29.53125 * std::sqrt(143.0);
    const auto f_160 = 1.40625 * std::sqrt(143.0);
    const auto f_161 = 0.046875 * std::sqrt(4290.0);
    const auto f_162 = 0.109375 * std::sqrt(4290.0);
    const auto f_163 = 0.65625 * std::sqrt(4290.0);
    const auto f_164 = 2.1875 * std::sqrt(4290.0);
    const auto f_165 = 0.28125 * std::sqrt(4290.0);
    const auto f_166 = 3.9375 * std::sqrt(4290.0);
    const auto f_167 = 13.125 * std::sqrt(4290.0);
    const auto f_168 = 0.234375 * std::sqrt(5005.0);
    const auto f_169 = 0.9375 * std::sqrt(5005.0);
    const auto f_170 = 0.421875 * std::sqrt(5005.0);
    const auto f_171 = 1.875 * std::sqrt(5005.0);
    const auto f_172 = 0.046875 * std::sqrt(5005.0);
    const auto f_173 = 0.1875 * std::sqrt(5005.0);
    const auto f_174 = 1.40625 * std::sqrt(5005.0);
    const auto f_175 = 5.625 * std::sqrt(5005.0);
    const auto f_176 = 2.53125 * std::sqrt(5005.0);
    const auto f_177 = 11.25 * std::sqrt(5005.0);
    const auto f_178 = 0.28125 * std::sqrt(5005.0);
    const auto f_179 = 1.125 * std::sqrt(5005.0);
    const auto f_180 = 0.09375 * std::sqrt(385.0);
    const auto f_181 = 2.25 * std::sqrt(385.0);
    const auto f_182 = 3.75 * std::sqrt(385.0);
    const auto f_183 = 0.5625 * std::sqrt(385.0);
    const auto f_184 = 13.5 * std::sqrt(385.0);
    const auto f_185 = 22.5 * std::sqrt(385.0);
    const auto f_186 = 0.703125 * std::sqrt(231.0);
    const auto f_187 = 1.171875 * std::sqrt(231.0);
    const auto f_188 = 4.6875 * std::sqrt(231.0);
    const auto f_189 = 0.234375 * std::sqrt(231.0);
    const auto f_190 = 3.125 * std::sqrt(231.0);
    const auto f_191 = 3.75 * std::sqrt(231.0);
    const auto f_192 = 1.5625 * std::sqrt(231.0);
    const auto f_193 = 1.25 * std::sqrt(231.0);
    const auto f_194 = 4.21875 * std::sqrt(231.0);
    const auto f_195 = 7.03125 * std::sqrt(231.0);
    const auto f_196 = 28.125 * std::sqrt(231.0);
    const auto f_197 = 1.40625 * std::sqrt(231.0);
    const auto f_198 = 18.75 * std::sqrt(231.0);
    const auto f_199 = 22.5 * std::sqrt(231.0);
    const auto f_200 = 9.375 * std::sqrt(231.0);
    const auto f_201 = 7.5 * std::sqrt(231.0);
    const auto f_202 = 0.234375 * std::sqrt(14.0);
    const auto f_203 = 0.703125 * std::sqrt(14.0);
    const auto f_204 = 7.03125 * std::sqrt(14.0);
    const auto f_205 = 14.0625 * std::sqrt(14.0);
    const auto f_206 = 18.75 * std::sqrt(14.0);
    const auto f_207 = 7.5 * std::sqrt(14.0);
    const auto f_208 = 1.40625 * std::sqrt(14.0);
    const auto f_209 = 4.21875 * std::sqrt(14.0);
    const auto f_210 = 42.1875 * std::sqrt(14.0);
    const auto f_211 = 84.375 * std::sqrt(14.0);
    const auto f_212 = 112.5 * std::sqrt(14.0);
    const auto f_213 = 45.0 * std::sqrt(14.0);
    const auto f_214 = 1.640625 * std::sqrt(5.0);
    const auto f_215 = 4.921875 * std::sqrt(5.0);
    const auto f_216 = 13.125 * std::sqrt(5.0);
    const auto f_217 = 26.25 * std::sqrt(5.0);
    const auto f_218 = 15.75 * std::sqrt(5.0);
    const auto f_219 = 3.0 * std::sqrt(5.0);
    const auto f_220 = 9.84375 * std::sqrt(5.0);
    const auto f_221 = 29.53125 * std::sqrt(5.0);
    const auto f_222 = 78.75 * std::sqrt(5.0);
    const auto f_223 = 157.5 * std::sqrt(5.0);
    const auto f_224 = 94.5 * std::sqrt(5.0);
    const auto f_225 = 18.0 * std::sqrt(5.0);
    const auto f_226 = 0.13671875 * std::sqrt(5.0);
    const auto f_227 = 0.546875 * std::sqrt(5.0);
    const auto f_228 = 4.375 * std::sqrt(5.0);
    const auto f_229 = 0.8203125 * std::sqrt(5.0);
    const auto f_230 = 7.0 * std::sqrt(5.0);
    const auto f_231 = 0.5 * std::sqrt(5.0);
    const auto f_232 = 3.28125 * std::sqrt(5.0);
    const auto f_233 = 42.0 * std::sqrt(5.0);
    const auto f_234 = 0.1171875 * std::sqrt(14.0);
    const auto f_235 = 3.515625 * std::sqrt(14.0);
    const auto f_236 = 9.375 * std::sqrt(14.0);
    const auto f_237 = 3.75 * std::sqrt(14.0);
    const auto f_238 = 21.09375 * std::sqrt(14.0);
    const auto f_239 = 56.25 * std::sqrt(14.0);
    const auto f_240 = 22.5 * std::sqrt(14.0);
    const auto f_241 = 0.0234375 * std::sqrt(385.0);
    const auto f_242 = 0.234375 * std::sqrt(385.0);
    const auto f_243 = 2.8125 * std::sqrt(385.0);
    const auto f_244 = 0.9375 * std::sqrt(385.0);
    const auto f_245 = 5.625 * std::sqrt(385.0);
    const auto f_246 = 0.140625 * std::sqrt(385.0);
    const auto f_247 = 3.375 * std::sqrt(385.0);
    const auto f_248 = 1.40625 * std::sqrt(385.0);
    const auto f_249 = 16.875 * std::sqrt(385.0);
    const auto f_250 = 33.75 * std::sqrt(385.0);
    const auto f_251 = 0.0078125 * std::sqrt(4290.0);
    const auto f_252 = 1.640625 * std::sqrt(4290.0);
    const auto f_253 = 9.84375 * std::sqrt(4290.0);
    const auto f_254 = 0.05859375 * std::sqrt(143.0);
    const auto f_255 = 4.1015625 * std::sqrt(143.0);
    const auto f_256 = 0.3515625 * std::sqrt(143.0);
    const auto f_257 = 24.609375 * std::sqrt(143.0);
    const auto f_258 = 0.703125 * std::sqrt(286.0);
    const auto f_259 = 4.921875 * std::sqrt(286.0);
    const auto f_260 = 0.9375 * std::sqrt(286.0);
    const auto f_261 = 6.5625 * std::sqrt(286.0);
    const auto f_262 = 2.4609375 * std::sqrt(286.0);
    const auto f_263 = 12.3046875 * std::sqrt(286.0);
    const auto f_264 = 7.3828125 * std::sqrt(286.0);
    const auto f_265 = 0.3515625 * std::sqrt(286.0);
    const auto f_266 = 3.28125 * std::sqrt(286.0);
    const auto f_267 = 16.40625 * std::sqrt(286.0);
    const auto f_268 = 9.84375 * std::sqrt(286.0);
    const auto f_269 = 0.46875 * std::sqrt(286.0);
    const auto f_270 = 0.140625 * std::sqrt(2145.0);
    const auto f_271 = 0.328125 * std::sqrt(2145.0);
    const auto f_272 = 1.96875 * std::sqrt(2145.0);
    const auto f_273 = 6.5625 * std::sqrt(2145.0);
    const auto f_274 = 0.1875 * std::sqrt(2145.0);
    const auto f_275 = 0.4375 * std::sqrt(2145.0);
    const auto f_276 = 2.625 * std::sqrt(2145.0);
    const auto f_277 = 8.75 * std::sqrt(2145.0);
    const auto f_278 = 0.3515625 * std::sqrt(10010.0);
    const auto f_279 = 1.40625 * std::sqrt(10010.0);
    const auto f_280 = 0.6328125 * std::sqrt(10010.0);
    const auto f_281 = 2.8125 * std::sqrt(10010.0);
    const auto f_282 = 0.0703125 * std::sqrt(10010.0);
    const auto f_283 = 0.28125 * std::sqrt(10010.0);
    const auto f_284 = 0.46875 * std::sqrt(10010.0);
    const auto f_285 = 1.875 * std::sqrt(10010.0);
    const auto f_286 = 0.84375 * std::sqrt(10010.0);
    const auto f_287 = 3.75 * std::sqrt(10010.0);
    const auto f_288 = 0.09375 * std::sqrt(10010.0);
    const auto f_289 = 0.375 * std::sqrt(10010.0);
    const auto f_290 = 0.140625 * std::sqrt(770.0);
    const auto f_291 = 3.375 * std::sqrt(770.0);
    const auto f_292 = 5.625 * std::sqrt(770.0);
    const auto f_293 = 0.1875 * std::sqrt(770.0);
    const auto f_294 = 4.5 * std::sqrt(770.0);
    const auto f_295 = 7.5 * std::sqrt(770.0);
    const auto f_296 = 1.0546875 * std::sqrt(462.0);
    const auto f_297 = 1.7578125 * std::sqrt(462.0);
    const auto f_298 = 7.03125 * std::sqrt(462.0);
    const auto f_299 = 0.3515625 * std::sqrt(462.0);
    const auto f_300 = 4.6875 * std::sqrt(462.0);
    const auto f_301 = 5.625 * std::sqrt(462.0);
    const auto f_302 = 2.34375 * std::sqrt(462.0);
    const auto f_303 = 1.875 * std::sqrt(462.0);
    const auto f_304 = 1.40625 * std::sqrt(462.0);
    const auto f_305 = 9.375 * std::sqrt(462.0);
    const auto f_306 = 0.46875 * std::sqrt(462.0);
    const auto f_307 = 6.25 * std::sqrt(462.0);
    const auto f_308 = 7.5 * std::sqrt(462.0);
    const auto f_309 = 3.125 * std::sqrt(462.0);
    const auto f_310 = 2.5 * std::sqrt(462.0);
    const auto f_311 = 0.703125 * std::sqrt(7.0);
    const auto f_312 = 2.109375 * std::sqrt(7.0);
    const auto f_313 = 21.09375 * std::sqrt(7.0);
    const auto f_314 = 42.1875 * std::sqrt(7.0);
    const auto f_315 = 56.25 * std::sqrt(7.0);
    const auto f_316 = 22.5 * std::sqrt(7.0);
    const auto f_317 = 0.9375 * std::sqrt(7.0);
    const auto f_318 = 2.8125 * std::sqrt(7.0);
    const auto f_319 = 28.125 * std::sqrt(7.0);
    const auto f_320 = 75.0 * std::sqrt(7.0);
    const auto f_321 = 30.0 * std::sqrt(7.0);
    const auto f_322 = 2.4609375 * std::sqrt(10.0);
    const auto f_323 = 7.3828125 * std::sqrt(10.0);
    const auto f_324 = 19.6875 * std::sqrt(10.0);
    const auto f_325 = 39.375 * std::sqrt(10.0);
    const auto f_326 = 23.625 * std::sqrt(10.0);
    const auto f_327 = 4.5 * std::sqrt(10.0);
    const auto f_328 = 3.28125 * std::sqrt(10.0);
    const auto f_329 = 9.84375 * std::sqrt(10.0);
    const auto f_330 = 26.25 * std::sqrt(10.0);
    const auto f_331 = 52.5 * std::sqrt(10.0);
    const auto f_332 = 31.5 * std::sqrt(10.0);
    const auto f_333 = 6.0 * std::sqrt(10.0);
    const auto f_334 = 0.205078125 * std::sqrt(10.0);
    const auto f_335 = 0.8203125 * std::sqrt(10.0);
    const auto f_336 = 6.5625 * std::sqrt(10.0);
    const auto f_337 = 1.23046875 * std::sqrt(10.0);
    const auto f_338 = 10.5 * std::sqrt(10.0);
    const auto f_339 = 0.75 * std::sqrt(10.0);
    const auto f_340 = 0.2734375 * std::sqrt(10.0);
    const auto f_341 = 1.09375 * std::sqrt(10.0);
    const auto f_342 = 8.75 * std::sqrt(10.0);
    const auto f_343 = 1.640625 * std::sqrt(10.0);
    const auto f_344 = 14.0 * std::sqrt(10.0);
    const auto f_345 = std::sqrt(10.0);
    const auto f_346 = 0.3515625 * std::sqrt(7.0);
    const auto f_347 = 10.546875 * std::sqrt(7.0);
    const auto f_348 = 11.25 * std::sqrt(7.0);
    const auto f_349 = 0.46875 * std::sqrt(7.0);
    const auto f_350 = 14.0625 * std::sqrt(7.0);
    const auto f_351 = 37.5 * std::sqrt(7.0);
    const auto f_352 = 15.0 * std::sqrt(7.0);
    const auto f_353 = 0.03515625 * std::sqrt(770.0);
    const auto f_354 = 0.84375 * std::sqrt(770.0);
    const auto f_355 = 0.3515625 * std::sqrt(770.0);
    const auto f_356 = 4.21875 * std::sqrt(770.0);
    const auto f_357 = 1.40625 * std::sqrt(770.0);
    const auto f_358 = 8.4375 * std::sqrt(770.0);
    const auto f_359 = 0.046875 * std::sqrt(770.0);
    const auto f_360 = 1.125 * std::sqrt(770.0);
    const auto f_361 = 0.46875 * std::sqrt(770.0);
    const auto f_362 = 1.875 * std::sqrt(770.0);
    const auto f_363 = 11.25 * std::sqrt(770.0);
    const auto f_364 = 0.0234375 * std::sqrt(2145.0);
    const auto f_365 = 4.921875 * std::sqrt(2145.0);
    const auto f_366 = 0.03125 * std::sqrt(2145.0);
    const auto f_367 = 0.087890625 * std::sqrt(286.0);
    const auto f_368 = 6.15234375 * std::sqrt(286.0);
    const auto f_369 = 0.1171875 * std::sqrt(286.0);
    const auto f_370 = 8.203125 * std::sqrt(286.0);
    const auto f_371 = 0.0703125 * std::sqrt(715.0);
    const auto f_372 = 0.4921875 * std::sqrt(715.0);
    const auto f_373 = 0.140625 * std::sqrt(715.0);
    const auto f_374 = 0.984375 * std::sqrt(715.0);
    const auto f_375 = 0.5625 * std::sqrt(715.0);
    const auto f_376 = 3.9375 * std::sqrt(715.0);
    const auto f_377 = 0.1875 * std::sqrt(715.0);
    const auto f_378 = 0.24609375 * std::sqrt(715.0);
    const auto f_379 = 1.23046875 * std::sqrt(715.0);
    const auto f_380 = 0.73828125 * std::sqrt(715.0);
    const auto f_381 = 0.03515625 * std::sqrt(715.0);
    const auto f_382 = 2.4609375 * std::sqrt(715.0);
    const auto f_383 = 1.4765625 * std::sqrt(715.0);
    const auto f_384 = 1.96875 * std::sqrt(715.0);
    const auto f_385 = 9.84375 * std::sqrt(715.0);
    const auto f_386 = 5.90625 * std::sqrt(715.0);
    const auto f_387 = 0.28125 * std::sqrt(715.0);
    const auto f_388 = 0.65625 * std::sqrt(715.0);
    const auto f_389 = 3.28125 * std::sqrt(715.0);
    const auto f_390 = 0.09375 * std::sqrt(715.0);
    const auto f_391 = 0.03515625 * std::sqrt(858.0);
    const auto f_392 = 0.08203125 * std::sqrt(858.0);
    const auto f_393 = 0.4921875 * std::sqrt(858.0);
    const auto f_394 = 1.640625 * std::sqrt(858.0);
    const auto f_395 = 0.0703125 * std::sqrt(858.0);
    const auto f_396 = 0.1640625 * std::sqrt(858.0);
    const auto f_397 = 0.984375 * std::sqrt(858.0);
    const auto f_398 = 3.28125 * std::sqrt(858.0);
    const auto f_399 = 0.28125 * std::sqrt(858.0);
    const auto f_400 = 0.65625 * std::sqrt(858.0);
    const auto f_401 = 3.9375 * std::sqrt(858.0);
    const auto f_402 = 13.125 * std::sqrt(858.0);
    const auto f_403 = 0.09375 * std::sqrt(858.0);
    const auto f_404 = 0.21875 * std::sqrt(858.0);
    const auto f_405 = 1.3125 * std::sqrt(858.0);
    const auto f_406 = 4.375 * std::sqrt(858.0);
    const auto f_407 = 0.17578125 * std::sqrt(1001.0);
    const auto f_408 = 0.703125 * std::sqrt(1001.0);
    const auto f_409 = 0.31640625 * std::sqrt(1001.0);
    const auto f_410 = 1.40625 * std::sqrt(1001.0);
    const auto f_411 = 0.03515625 * std::sqrt(1001.0);
    const auto f_412 = 0.140625 * std::sqrt(1001.0);
    const auto f_413 = 0.3515625 * std::sqrt(1001.0);
    const auto f_414 = 0.6328125 * std::sqrt(1001.0);
    const auto f_415 = 2.8125 * std::sqrt(1001.0);
    const auto f_416 = 0.0703125 * std::sqrt(1001.0);
    const auto f_417 = 0.28125 * std::sqrt(1001.0);
    const auto f_418 = 5.625 * std::sqrt(1001.0);
    const auto f_419 = 2.53125 * std::sqrt(1001.0);
    const auto f_420 = 11.25 * std::sqrt(1001.0);
    const auto f_421 = 1.125 * std::sqrt(1001.0);
    const auto f_422 = 1.875 * std::sqrt(1001.0);
    const auto f_423 = 0.84375 * std::sqrt(1001.0);
    const auto f_424 = 3.75 * std::sqrt(1001.0);
    const auto f_425 = 0.09375 * std::sqrt(1001.0);
    const auto f_426 = 0.375 * std::sqrt(1001.0);
    const auto f_427 = 0.0703125 * std::sqrt(77.0);
    const auto f_428 = 1.6875 * std::sqrt(77.0);
    const auto f_429 = 2.8125 * std::sqrt(77.0);
    const auto f_430 = 0.140625 * std::sqrt(77.0);
    const auto f_431 = 3.375 * std::sqrt(77.0);
    const auto f_432 = 5.625 * std::sqrt(77.0);
    const auto f_433 = 0.5625 * std::sqrt(77.0);
    const auto f_434 = 13.5 * std::sqrt(77.0);
    const auto f_435 = 22.5 * std::sqrt(77.0);
    const auto f_436 = 0.1875 * std::sqrt(77.0);
    const auto f_437 = 4.5 * std::sqrt(77.0);
    const auto f_438 = 7.5 * std::sqrt(77.0);
    const auto f_439 = 0.10546875 * std::sqrt(1155.0);
    const auto f_440 = 0.17578125 * std::sqrt(1155.0);
    const auto f_441 = 0.703125 * std::sqrt(1155.0);
    const auto f_442 = 0.03515625 * std::sqrt(1155.0);
    const auto f_443 = 0.46875 * std::sqrt(1155.0);
    const auto f_444 = 0.5625 * std::sqrt(1155.0);
    const auto f_445 = 0.234375 * std::sqrt(1155.0);
    const auto f_446 = 0.1875 * std::sqrt(1155.0);
    const auto f_447 = 0.2109375 * std::sqrt(1155.0);
    const auto f_448 = 0.3515625 * std::sqrt(1155.0);
    const auto f_449 = 1.40625 * std::sqrt(1155.0);
    const auto f_450 = 0.0703125 * std::sqrt(1155.0);
    const auto f_451 = 0.9375 * std::sqrt(1155.0);
    const auto f_452 = 1.125 * std::sqrt(1155.0);
    const auto f_453 = 0.375 * std::sqrt(1155.0);
    const auto f_454 = 0.84375 * std::sqrt(1155.0);
    const auto f_455 = 5.625 * std::sqrt(1155.0);
    const auto f_456 = 0.28125 * std::sqrt(1155.0);
    const auto f_457 = 3.75 * std::sqrt(1155.0);
    const auto f_458 = 4.5 * std::sqrt(1155.0);
    const auto f_459 = 1.875 * std::sqrt(1155.0);
    const auto f_460 = 1.5 * std::sqrt(1155.0);
    const auto f_461 = 0.09375 * std::sqrt(1155.0);
    const auto f_462 = 1.25 * std::sqrt(1155.0);
    const auto f_463 = 0.625 * std::sqrt(1155.0);
    const auto f_464 = 0.5 * std::sqrt(1155.0);
    const auto f_465 = 0.03515625 * std::sqrt(70.0);
    const auto f_466 = 0.10546875 * std::sqrt(70.0);
    const auto f_467 = 1.0546875 * std::sqrt(70.0);
    const auto f_468 = 2.109375 * std::sqrt(70.0);
    const auto f_469 = 2.8125 * std::sqrt(70.0);
    const auto f_470 = 1.125 * std::sqrt(70.0);
    const auto f_471 = 0.0703125 * std::sqrt(70.0);
    const auto f_472 = 0.2109375 * std::sqrt(70.0);
    const auto f_473 = 4.21875 * std::sqrt(70.0);
    const auto f_474 = 5.625 * std::sqrt(70.0);
    const auto f_475 = 2.25 * std::sqrt(70.0);
    const auto f_476 = 0.28125 * std::sqrt(70.0);
    const auto f_477 = 0.84375 * std::sqrt(70.0);
    const auto f_478 = 8.4375 * std::sqrt(70.0);
    const auto f_479 = 16.875 * std::sqrt(70.0);
    const auto f_480 = 22.5 * std::sqrt(70.0);
    const auto f_481 = 9.0 * std::sqrt(70.0);
    const auto f_482 = 0.09375 * std::sqrt(70.0);
    const auto f_483 = 7.5 * std::sqrt(70.0);
    const auto f_484 = 3.0 * std::sqrt(70.0);
    const auto f_485 = 0.017578125 * std::sqrt(70.0);
    const auto f_486 = 0.52734375 * std::sqrt(70.0);
    const auto f_487 = 1.40625 * std::sqrt(70.0);
    const auto f_488 = 0.5625 * std::sqrt(70.0);
    const auto f_489 = 0.140625 * std::sqrt(70.0);
    const auto f_490 = 11.25 * std::sqrt(70.0);
    const auto f_491 = 0.046875 * std::sqrt(70.0);
    const auto f_492 = 3.75 * std::sqrt(70.0);
    const auto f_493 = 0.017578125 * std::sqrt(77.0);
    const auto f_494 = 0.421875 * std::sqrt(77.0);
    const auto f_495 = 0.17578125 * std::sqrt(77.0);
    const auto f_496 = 2.109375 * std::sqrt(77.0);
    const auto f_497 = 0.703125 * std::sqrt(77.0);
    const auto f_498 = 4.21875 * std::sqrt(77.0);
    const auto f_499 = 0.03515625 * std::sqrt(77.0);
    const auto f_500 = 0.84375 * std::sqrt(77.0);
    const auto f_501 = 0.3515625 * std::sqrt(77.0);
    const auto f_502 = 1.40625 * std::sqrt(77.0);
    const auto f_503 = 8.4375 * std::sqrt(77.0);
    const auto f_504 = 16.875 * std::sqrt(77.0);
    const auto f_505 = 33.75 * std::sqrt(77.0);
    const auto f_506 = 0.046875 * std::sqrt(77.0);
    const auto f_507 = 1.125 * std::sqrt(77.0);
    const auto f_508 = 0.46875 * std::sqrt(77.0);
    const auto f_509 = 1.875 * std::sqrt(77.0);
    const auto f_510 = 11.25 * std::sqrt(77.0);
    const auto f_511 = 0.005859375 * std::sqrt(858.0);
    const auto f_512 = 1.23046875 * std::sqrt(858.0);
    const auto f_513 = 0.01171875 * std::sqrt(858.0);
    const auto f_514 = 2.4609375 * std::sqrt(858.0);
    const auto f_515 = 0.046875 * std::sqrt(858.0);
    const auto f_516 = 9.84375 * std::sqrt(858.0);
    const auto f_517 = 0.015625 * std::sqrt(858.0);
    const auto f_518 = 0.0087890625 * std::sqrt(715.0);
    const auto f_519 = 0.615234375 * std::sqrt(715.0);
    const auto f_520 = 0.017578125 * std::sqrt(715.0);
    const auto f_521 = 4.921875 * std::sqrt(715.0);
    const auto f_522 = 0.0234375 * std::sqrt(715.0);
    const auto f_523 = 0.8203125 * std::sqrt(143.0);
    const auto f_524 = 2.4609375 * std::sqrt(143.0);
    const auto f_525 = 0.1171875 * std::sqrt(143.0);
    const auto f_526 = 14.765625 * std::sqrt(143.0);
    const auto f_527 = 0.703125 * std::sqrt(143.0);
    const auto f_528 = 0.0234375 * std::sqrt(4290.0);
    const auto f_529 = 0.0546875 * std::sqrt(4290.0);
    const auto f_530 = 0.328125 * std::sqrt(4290.0);
    const auto f_531 = 1.09375 * std::sqrt(4290.0);
    const auto f_532 = 0.140625 * std::sqrt(4290.0);
    const auto f_533 = 1.96875 * std::sqrt(4290.0);
    const auto f_534 = 6.5625 * std::sqrt(4290.0);
    const auto f_535 = 0.1171875 * std::sqrt(5005.0);
    const auto f_536 = 0.46875 * std::sqrt(5005.0);
    const auto f_537 = 0.2109375 * std::sqrt(5005.0);
    const auto f_538 = 0.0234375 * std::sqrt(5005.0);
    const auto f_539 = 0.09375 * std::sqrt(5005.0);
    const auto f_540 = 0.703125 * std::sqrt(5005.0);
    const auto f_541 = 2.8125 * std::sqrt(5005.0);
    const auto f_542 = 1.265625 * std::sqrt(5005.0);
    const auto f_543 = 0.140625 * std::sqrt(5005.0);
    const auto f_544 = 0.5625 * std::sqrt(5005.0);
    const auto f_545 = 0.046875 * std::sqrt(385.0);
    const auto f_546 = 1.125 * std::sqrt(385.0);
    const auto f_547 = 1.875 * std::sqrt(385.0);
    const auto f_548 = 0.28125 * std::sqrt(385.0);
    const auto f_549 = 6.75 * std::sqrt(385.0);
    const auto f_550 = 11.25 * std::sqrt(385.0);
    const auto f_551 = 0.3515625 * std::sqrt(231.0);
    const auto f_552 = 0.5859375 * std::sqrt(231.0);
    const auto f_553 = 2.34375 * std::sqrt(231.0);
    const auto f_554 = 0.1171875 * std::sqrt(231.0);
    const auto f_555 = 1.875 * std::sqrt(231.0);
    const auto f_556 = 0.78125 * std::sqrt(231.0);
    const auto f_557 = 0.625 * std::sqrt(231.0);
    const auto f_558 = 2.109375 * std::sqrt(231.0);
    const auto f_559 = 3.515625 * std::sqrt(231.0);
    const auto f_560 = 14.0625 * std::sqrt(231.0);
    const auto f_561 = 11.25 * std::sqrt(231.0);
    const auto f_562 = 0.3515625 * std::sqrt(14.0);
    const auto f_563 = 2.109375 * std::sqrt(14.0);
    const auto f_564 = 2.4609375 * std::sqrt(5.0);
    const auto f_565 = 6.5625 * std::sqrt(5.0);
    const auto f_566 = 7.875 * std::sqrt(5.0);
    const auto f_567 = 1.5 * std::sqrt(5.0);
    const auto f_568 = 14.765625 * std::sqrt(5.0);
    const auto f_569 = 39.375 * std::sqrt(5.0);
    const auto f_570 = 47.25 * std::sqrt(5.0);
    const auto f_571 = 9.0 * std::sqrt(5.0);
    const auto f_572 = 0.068359375 * std::sqrt(5.0);
    const auto f_573 = 0.2734375 * std::sqrt(5.0);
    const auto f_574 = 2.1875 * std::sqrt(5.0);
    const auto f_575 = 0.41015625 * std::sqrt(5.0);
    const auto f_576 = 3.5 * std::sqrt(5.0);
    const auto f_577 = 0.25 * std::sqrt(5.0);
    const auto f_578 = 21.0 * std::sqrt(5.0);
    const auto f_579 = 0.05859375 * std::sqrt(14.0);
    const auto f_580 = 1.7578125 * std::sqrt(14.0);
    const auto f_581 = 4.6875 * std::sqrt(14.0);
    const auto f_582 = 1.875 * std::sqrt(14.0);
    const auto f_583 = 10.546875 * std::sqrt(14.0);
    const auto f_584 = 28.125 * std::sqrt(14.0);
    const auto f_585 = 11.25 * std::sqrt(14.0);
    const auto f_586 = 0.01171875 * std::sqrt(385.0);
    const auto f_587 = 0.1171875 * std::sqrt(385.0);
    const auto f_588 = 0.46875 * std::sqrt(385.0);
    const auto f_589 = 0.0703125 * std::sqrt(385.0);
    const auto f_590 = 1.6875 * std::sqrt(385.0);
    const auto f_591 = 0.703125 * std::sqrt(385.0);
    const auto f_592 = 8.4375 * std::sqrt(385.0);
    const auto f_593 = 0.00390625 * std::sqrt(4290.0);
    const auto f_594 = 0.8203125 * std::sqrt(4290.0);
    const auto f_595 = 4.921875 * std::sqrt(4290.0);
    const auto f_596 = 0.029296875 * std::sqrt(143.0);
    const auto f_597 = 2.05078125 * std::sqrt(143.0);
    const auto f_598 = 0.17578125 * std::sqrt(143.0);
    const auto f_599 = 12.3046875 * std::sqrt(143.0);
    const auto f_600 = 0.1171875 * std::sqrt(1001.0);
    const auto f_601 = 0.8203125 * std::sqrt(1001.0);
    const auto f_602 = 0.41015625 * std::sqrt(1001.0);
    const auto f_603 = 2.05078125 * std::sqrt(1001.0);
    const auto f_604 = 1.23046875 * std::sqrt(1001.0);
    const auto f_605 = 2.4609375 * std::sqrt(1001.0);
    const auto f_606 = 12.3046875 * std::sqrt(1001.0);
    const auto f_607 = 7.3828125 * std::sqrt(1001.0);
    const auto f_608 = 0.01171875 * std::sqrt(30030.0);
    const auto f_609 = 0.02734375 * std::sqrt(30030.0);
    const auto f_610 = 0.1640625 * std::sqrt(30030.0);
    const auto f_611 = 0.546875 * std::sqrt(30030.0);
    const auto f_612 = 0.0703125 * std::sqrt(30030.0);
    const auto f_613 = 0.984375 * std::sqrt(30030.0);
    const auto f_614 = 3.28125 * std::sqrt(30030.0);
    const auto f_615 = 0.41015625 * std::sqrt(715.0);
    const auto f_616 = 0.08203125 * std::sqrt(715.0);
    const auto f_617 = 4.4296875 * std::sqrt(715.0);
    const auto f_618 = 19.6875 * std::sqrt(715.0);
    const auto f_619 = 0.984375 * std::sqrt(55.0);
    const auto f_620 = 23.625 * std::sqrt(55.0);
    const auto f_621 = 1.23046875 * std::sqrt(33.0);
    const auto f_622 = 2.05078125 * std::sqrt(33.0);
    const auto f_623 = 0.41015625 * std::sqrt(33.0);
    const auto f_624 = 5.46875 * std::sqrt(33.0);
    const auto f_625 = 6.5625 * std::sqrt(33.0);
    const auto f_626 = 2.734375 * std::sqrt(33.0);
    const auto f_627 = 2.1875 * std::sqrt(33.0);
    const auto f_628 = 7.3828125 * std::sqrt(33.0);
    const auto f_629 = 12.3046875 * std::sqrt(33.0);
    const auto f_630 = 49.21875 * std::sqrt(33.0);
    const auto f_631 = 2.4609375 * std::sqrt(33.0);
    const auto f_632 = 39.375 * std::sqrt(33.0);
    const auto f_633 = 16.40625 * std::sqrt(33.0);
    const auto f_634 = 13.125 * std::sqrt(33.0);
    const auto f_635 = 0.41015625 * std::sqrt(2.0);
    const auto f_636 = 1.23046875 * std::sqrt(2.0);
    const auto f_637 = 12.3046875 * std::sqrt(2.0);
    const auto f_638 = 32.8125 * std::sqrt(2.0);
    const auto f_639 = 13.125 * std::sqrt(2.0);
    const auto f_640 = 2.4609375 * std::sqrt(2.0);
    const auto f_641 = 7.3828125 * std::sqrt(2.0);
    const auto f_642 = 73.828125 * std::sqrt(2.0);
    const auto f_643 = 147.65625 * std::sqrt(2.0);
    const auto f_644 = 196.875 * std::sqrt(2.0);
    const auto f_645 = 78.75 * std::sqrt(2.0);
    const auto f_646 = 0.41015625 * std::sqrt(35.0);
    const auto f_647 = 1.23046875 * std::sqrt(35.0);
    const auto f_648 = 3.28125 * std::sqrt(35.0);
    const auto f_649 = 6.5625 * std::sqrt(35.0);
    const auto f_650 = 3.9375 * std::sqrt(35.0);
    const auto f_651 = 0.75 * std::sqrt(35.0);
    const auto f_652 = 2.4609375 * std::sqrt(35.0);
    const auto f_653 = 7.3828125 * std::sqrt(35.0);
    const auto f_654 = 19.6875 * std::sqrt(35.0);
    const auto f_655 = 39.375 * std::sqrt(35.0);
    const auto f_656 = 23.625 * std::sqrt(35.0);
    const auto f_657 = 4.5 * std::sqrt(35.0);
    const auto f_658 = 0.0341796875 * std::sqrt(35.0);
    const auto f_659 = 1.09375 * std::sqrt(35.0);
    const auto f_660 = 0.205078125 * std::sqrt(35.0);
    const auto f_661 = 1.75 * std::sqrt(35.0);
    const auto f_662 = 0.125 * std::sqrt(35.0);
    const auto f_663 = 10.5 * std::sqrt(35.0);
    const auto f_664 = 0.205078125 * std::sqrt(2.0);
    const auto f_665 = 6.15234375 * std::sqrt(2.0);
    const auto f_666 = 16.40625 * std::sqrt(2.0);
    const auto f_667 = 6.5625 * std::sqrt(2.0);
    const auto f_668 = 36.9140625 * std::sqrt(2.0);
    const auto f_669 = 39.375 * std::sqrt(2.0);
    const auto f_670 = 0.041015625 * std::sqrt(55.0);
    const auto f_671 = 0.41015625 * std::sqrt(55.0);
    const auto f_672 = 4.921875 * std::sqrt(55.0);
    const auto f_673 = 9.84375 * std::sqrt(55.0);
    const auto f_674 = 0.24609375 * std::sqrt(55.0);
    const auto f_675 = 5.90625 * std::sqrt(55.0);
    const auto f_676 = 2.4609375 * std::sqrt(55.0);
    const auto f_677 = 29.53125 * std::sqrt(55.0);
    const auto f_678 = 59.0625 * std::sqrt(55.0);
    const auto f_679 = 0.001953125 * std::sqrt(30030.0);
    const auto f_680 = 0.41015625 * std::sqrt(30030.0);
    const auto f_681 = 2.4609375 * std::sqrt(30030.0);
    const auto f_682 = 0.0146484375 * std::sqrt(1001.0);
    const auto f_683 = 1.025390625 * std::sqrt(1001.0);
    const auto f_684 = 0.087890625 * std::sqrt(1001.0);
    const auto f_685 = 6.15234375 * std::sqrt(1001.0);

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
    auto *g_117 = values + 117 * nvalues;
    auto *g_118 = values + 118 * nvalues;
    auto *g_119 = values + 119 * nvalues;
    auto *g_120 = values + 120 * nvalues;
    auto *g_121 = values + 121 * nvalues;
    auto *g_122 = values + 122 * nvalues;
    auto *g_123 = values + 123 * nvalues;
    auto *g_124 = values + 124 * nvalues;
    auto *g_125 = values + 125 * nvalues;
    auto *g_126 = values + 126 * nvalues;
    auto *g_127 = values + 127 * nvalues;
    auto *g_128 = values + 128 * nvalues;
    auto *g_129 = values + 129 * nvalues;
    auto *g_130 = values + 130 * nvalues;
    auto *g_131 = values + 131 * nvalues;
    auto *g_132 = values + 132 * nvalues;
    auto *g_133 = values + 133 * nvalues;
    auto *g_134 = values + 134 * nvalues;
    auto *g_135 = values + 135 * nvalues;
    auto *g_136 = values + 136 * nvalues;
    auto *g_137 = values + 137 * nvalues;
    auto *g_138 = values + 138 * nvalues;
    auto *g_139 = values + 139 * nvalues;
    auto *g_140 = values + 140 * nvalues;
    auto *g_141 = values + 141 * nvalues;
    auto *g_142 = values + 142 * nvalues;
    auto *g_143 = values + 143 * nvalues;
    auto *g_144 = values + 144 * nvalues;
    auto *g_145 = values + 145 * nvalues;
    auto *g_146 = values + 146 * nvalues;
    auto *g_147 = values + 147 * nvalues;
    auto *g_148 = values + 148 * nvalues;
    auto *g_149 = values + 149 * nvalues;
    auto *g_150 = values + 150 * nvalues;
    auto *g_151 = values + 151 * nvalues;
    auto *g_152 = values + 152 * nvalues;

    const auto *gl_0 = buffer.data(gl + 0);
    const auto *gl_1 = buffer.data(gl + 1);
    const auto *gl_2 = buffer.data(gl + 2);
    const auto *gl_3 = buffer.data(gl + 3);
    const auto *gl_4 = buffer.data(gl + 4);
    const auto *gl_5 = buffer.data(gl + 5);
    const auto *gl_6 = buffer.data(gl + 6);
    const auto *gl_7 = buffer.data(gl + 7);
    const auto *gl_8 = buffer.data(gl + 8);
    const auto *gl_9 = buffer.data(gl + 9);
    const auto *gl_10 = buffer.data(gl + 10);
    const auto *gl_11 = buffer.data(gl + 11);
    const auto *gl_12 = buffer.data(gl + 12);
    const auto *gl_13 = buffer.data(gl + 13);
    const auto *gl_14 = buffer.data(gl + 14);
    const auto *gl_15 = buffer.data(gl + 15);
    const auto *gl_16 = buffer.data(gl + 16);
    const auto *gl_17 = buffer.data(gl + 17);
    const auto *gl_18 = buffer.data(gl + 18);
    const auto *gl_19 = buffer.data(gl + 19);
    const auto *gl_20 = buffer.data(gl + 20);
    const auto *gl_21 = buffer.data(gl + 21);
    const auto *gl_22 = buffer.data(gl + 22);
    const auto *gl_23 = buffer.data(gl + 23);
    const auto *gl_24 = buffer.data(gl + 24);
    const auto *gl_25 = buffer.data(gl + 25);
    const auto *gl_26 = buffer.data(gl + 26);
    const auto *gl_27 = buffer.data(gl + 27);
    const auto *gl_28 = buffer.data(gl + 28);
    const auto *gl_29 = buffer.data(gl + 29);
    const auto *gl_30 = buffer.data(gl + 30);
    const auto *gl_31 = buffer.data(gl + 31);
    const auto *gl_32 = buffer.data(gl + 32);
    const auto *gl_33 = buffer.data(gl + 33);
    const auto *gl_34 = buffer.data(gl + 34);
    const auto *gl_35 = buffer.data(gl + 35);
    const auto *gl_36 = buffer.data(gl + 36);
    const auto *gl_37 = buffer.data(gl + 37);
    const auto *gl_38 = buffer.data(gl + 38);
    const auto *gl_39 = buffer.data(gl + 39);
    const auto *gl_40 = buffer.data(gl + 40);
    const auto *gl_41 = buffer.data(gl + 41);
    const auto *gl_42 = buffer.data(gl + 42);
    const auto *gl_43 = buffer.data(gl + 43);
    const auto *gl_44 = buffer.data(gl + 44);
    const auto *gl_45 = buffer.data(gl + 45);
    const auto *gl_46 = buffer.data(gl + 46);
    const auto *gl_47 = buffer.data(gl + 47);
    const auto *gl_48 = buffer.data(gl + 48);
    const auto *gl_49 = buffer.data(gl + 49);
    const auto *gl_50 = buffer.data(gl + 50);
    const auto *gl_51 = buffer.data(gl + 51);
    const auto *gl_52 = buffer.data(gl + 52);
    const auto *gl_53 = buffer.data(gl + 53);
    const auto *gl_54 = buffer.data(gl + 54);
    const auto *gl_55 = buffer.data(gl + 55);
    const auto *gl_56 = buffer.data(gl + 56);
    const auto *gl_57 = buffer.data(gl + 57);
    const auto *gl_58 = buffer.data(gl + 58);
    const auto *gl_59 = buffer.data(gl + 59);
    const auto *gl_60 = buffer.data(gl + 60);
    const auto *gl_61 = buffer.data(gl + 61);
    const auto *gl_62 = buffer.data(gl + 62);
    const auto *gl_63 = buffer.data(gl + 63);
    const auto *gl_64 = buffer.data(gl + 64);
    const auto *gl_65 = buffer.data(gl + 65);
    const auto *gl_66 = buffer.data(gl + 66);
    const auto *gl_67 = buffer.data(gl + 67);
    const auto *gl_68 = buffer.data(gl + 68);
    const auto *gl_69 = buffer.data(gl + 69);
    const auto *gl_70 = buffer.data(gl + 70);
    const auto *gl_71 = buffer.data(gl + 71);
    const auto *gl_72 = buffer.data(gl + 72);
    const auto *gl_73 = buffer.data(gl + 73);
    const auto *gl_74 = buffer.data(gl + 74);
    const auto *gl_75 = buffer.data(gl + 75);
    const auto *gl_76 = buffer.data(gl + 76);
    const auto *gl_77 = buffer.data(gl + 77);
    const auto *gl_78 = buffer.data(gl + 78);
    const auto *gl_79 = buffer.data(gl + 79);
    const auto *gl_80 = buffer.data(gl + 80);
    const auto *gl_81 = buffer.data(gl + 81);
    const auto *gl_82 = buffer.data(gl + 82);
    const auto *gl_83 = buffer.data(gl + 83);
    const auto *gl_84 = buffer.data(gl + 84);
    const auto *gl_85 = buffer.data(gl + 85);
    const auto *gl_86 = buffer.data(gl + 86);
    const auto *gl_87 = buffer.data(gl + 87);
    const auto *gl_88 = buffer.data(gl + 88);
    const auto *gl_89 = buffer.data(gl + 89);
    const auto *gl_90 = buffer.data(gl + 90);
    const auto *gl_91 = buffer.data(gl + 91);
    const auto *gl_92 = buffer.data(gl + 92);
    const auto *gl_93 = buffer.data(gl + 93);
    const auto *gl_94 = buffer.data(gl + 94);
    const auto *gl_95 = buffer.data(gl + 95);
    const auto *gl_96 = buffer.data(gl + 96);
    const auto *gl_97 = buffer.data(gl + 97);
    const auto *gl_98 = buffer.data(gl + 98);
    const auto *gl_99 = buffer.data(gl + 99);
    const auto *gl_100 = buffer.data(gl + 100);
    const auto *gl_101 = buffer.data(gl + 101);
    const auto *gl_102 = buffer.data(gl + 102);
    const auto *gl_103 = buffer.data(gl + 103);
    const auto *gl_104 = buffer.data(gl + 104);
    const auto *gl_105 = buffer.data(gl + 105);
    const auto *gl_106 = buffer.data(gl + 106);
    const auto *gl_107 = buffer.data(gl + 107);
    const auto *gl_108 = buffer.data(gl + 108);
    const auto *gl_109 = buffer.data(gl + 109);
    const auto *gl_110 = buffer.data(gl + 110);
    const auto *gl_111 = buffer.data(gl + 111);
    const auto *gl_112 = buffer.data(gl + 112);
    const auto *gl_113 = buffer.data(gl + 113);
    const auto *gl_114 = buffer.data(gl + 114);
    const auto *gl_115 = buffer.data(gl + 115);
    const auto *gl_116 = buffer.data(gl + 116);
    const auto *gl_117 = buffer.data(gl + 117);
    const auto *gl_118 = buffer.data(gl + 118);
    const auto *gl_119 = buffer.data(gl + 119);
    const auto *gl_120 = buffer.data(gl + 120);
    const auto *gl_121 = buffer.data(gl + 121);
    const auto *gl_122 = buffer.data(gl + 122);
    const auto *gl_123 = buffer.data(gl + 123);
    const auto *gl_124 = buffer.data(gl + 124);
    const auto *gl_125 = buffer.data(gl + 125);
    const auto *gl_126 = buffer.data(gl + 126);
    const auto *gl_127 = buffer.data(gl + 127);
    const auto *gl_128 = buffer.data(gl + 128);
    const auto *gl_129 = buffer.data(gl + 129);
    const auto *gl_130 = buffer.data(gl + 130);
    const auto *gl_131 = buffer.data(gl + 131);
    const auto *gl_132 = buffer.data(gl + 132);
    const auto *gl_133 = buffer.data(gl + 133);
    const auto *gl_134 = buffer.data(gl + 134);
    const auto *gl_135 = buffer.data(gl + 135);
    const auto *gl_136 = buffer.data(gl + 136);
    const auto *gl_137 = buffer.data(gl + 137);
    const auto *gl_138 = buffer.data(gl + 138);
    const auto *gl_139 = buffer.data(gl + 139);
    const auto *gl_140 = buffer.data(gl + 140);
    const auto *gl_141 = buffer.data(gl + 141);
    const auto *gl_142 = buffer.data(gl + 142);
    const auto *gl_143 = buffer.data(gl + 143);
    const auto *gl_144 = buffer.data(gl + 144);
    const auto *gl_145 = buffer.data(gl + 145);
    const auto *gl_146 = buffer.data(gl + 146);
    const auto *gl_147 = buffer.data(gl + 147);
    const auto *gl_148 = buffer.data(gl + 148);
    const auto *gl_149 = buffer.data(gl + 149);
    const auto *gl_150 = buffer.data(gl + 150);
    const auto *gl_151 = buffer.data(gl + 151);
    const auto *gl_152 = buffer.data(gl + 152);
    const auto *gl_153 = buffer.data(gl + 153);
    const auto *gl_154 = buffer.data(gl + 154);
    const auto *gl_155 = buffer.data(gl + 155);
    const auto *gl_156 = buffer.data(gl + 156);
    const auto *gl_157 = buffer.data(gl + 157);
    const auto *gl_158 = buffer.data(gl + 158);
    const auto *gl_159 = buffer.data(gl + 159);
    const auto *gl_160 = buffer.data(gl + 160);
    const auto *gl_161 = buffer.data(gl + 161);
    const auto *gl_162 = buffer.data(gl + 162);
    const auto *gl_163 = buffer.data(gl + 163);
    const auto *gl_164 = buffer.data(gl + 164);
    const auto *gl_165 = buffer.data(gl + 165);
    const auto *gl_166 = buffer.data(gl + 166);
    const auto *gl_167 = buffer.data(gl + 167);
    const auto *gl_168 = buffer.data(gl + 168);
    const auto *gl_169 = buffer.data(gl + 169);
    const auto *gl_170 = buffer.data(gl + 170);
    const auto *gl_171 = buffer.data(gl + 171);
    const auto *gl_172 = buffer.data(gl + 172);
    const auto *gl_173 = buffer.data(gl + 173);
    const auto *gl_174 = buffer.data(gl + 174);
    const auto *gl_175 = buffer.data(gl + 175);
    const auto *gl_176 = buffer.data(gl + 176);
    const auto *gl_177 = buffer.data(gl + 177);
    const auto *gl_178 = buffer.data(gl + 178);
    const auto *gl_179 = buffer.data(gl + 179);
    const auto *gl_180 = buffer.data(gl + 180);
    const auto *gl_181 = buffer.data(gl + 181);
    const auto *gl_182 = buffer.data(gl + 182);
    const auto *gl_183 = buffer.data(gl + 183);
    const auto *gl_184 = buffer.data(gl + 184);
    const auto *gl_185 = buffer.data(gl + 185);
    const auto *gl_186 = buffer.data(gl + 186);
    const auto *gl_187 = buffer.data(gl + 187);
    const auto *gl_188 = buffer.data(gl + 188);
    const auto *gl_189 = buffer.data(gl + 189);
    const auto *gl_190 = buffer.data(gl + 190);
    const auto *gl_191 = buffer.data(gl + 191);
    const auto *gl_192 = buffer.data(gl + 192);
    const auto *gl_193 = buffer.data(gl + 193);
    const auto *gl_194 = buffer.data(gl + 194);
    const auto *gl_195 = buffer.data(gl + 195);
    const auto *gl_196 = buffer.data(gl + 196);
    const auto *gl_197 = buffer.data(gl + 197);
    const auto *gl_198 = buffer.data(gl + 198);
    const auto *gl_199 = buffer.data(gl + 199);
    const auto *gl_200 = buffer.data(gl + 200);
    const auto *gl_201 = buffer.data(gl + 201);
    const auto *gl_202 = buffer.data(gl + 202);
    const auto *gl_203 = buffer.data(gl + 203);
    const auto *gl_204 = buffer.data(gl + 204);
    const auto *gl_205 = buffer.data(gl + 205);
    const auto *gl_206 = buffer.data(gl + 206);
    const auto *gl_207 = buffer.data(gl + 207);
    const auto *gl_208 = buffer.data(gl + 208);
    const auto *gl_209 = buffer.data(gl + 209);
    const auto *gl_210 = buffer.data(gl + 210);
    const auto *gl_211 = buffer.data(gl + 211);
    const auto *gl_212 = buffer.data(gl + 212);
    const auto *gl_213 = buffer.data(gl + 213);
    const auto *gl_214 = buffer.data(gl + 214);
    const auto *gl_215 = buffer.data(gl + 215);
    const auto *gl_216 = buffer.data(gl + 216);
    const auto *gl_217 = buffer.data(gl + 217);
    const auto *gl_218 = buffer.data(gl + 218);
    const auto *gl_219 = buffer.data(gl + 219);
    const auto *gl_220 = buffer.data(gl + 220);
    const auto *gl_221 = buffer.data(gl + 221);
    const auto *gl_222 = buffer.data(gl + 222);
    const auto *gl_223 = buffer.data(gl + 223);
    const auto *gl_224 = buffer.data(gl + 224);
    const auto *gl_225 = buffer.data(gl + 225);
    const auto *gl_226 = buffer.data(gl + 226);
    const auto *gl_227 = buffer.data(gl + 227);
    const auto *gl_228 = buffer.data(gl + 228);
    const auto *gl_229 = buffer.data(gl + 229);
    const auto *gl_230 = buffer.data(gl + 230);
    const auto *gl_231 = buffer.data(gl + 231);
    const auto *gl_232 = buffer.data(gl + 232);
    const auto *gl_233 = buffer.data(gl + 233);
    const auto *gl_234 = buffer.data(gl + 234);
    const auto *gl_235 = buffer.data(gl + 235);
    const auto *gl_236 = buffer.data(gl + 236);
    const auto *gl_237 = buffer.data(gl + 237);
    const auto *gl_238 = buffer.data(gl + 238);
    const auto *gl_239 = buffer.data(gl + 239);
    const auto *gl_240 = buffer.data(gl + 240);
    const auto *gl_241 = buffer.data(gl + 241);
    const auto *gl_242 = buffer.data(gl + 242);
    const auto *gl_243 = buffer.data(gl + 243);
    const auto *gl_244 = buffer.data(gl + 244);
    const auto *gl_245 = buffer.data(gl + 245);
    const auto *gl_246 = buffer.data(gl + 246);
    const auto *gl_247 = buffer.data(gl + 247);
    const auto *gl_248 = buffer.data(gl + 248);
    const auto *gl_249 = buffer.data(gl + 249);
    const auto *gl_250 = buffer.data(gl + 250);
    const auto *gl_251 = buffer.data(gl + 251);
    const auto *gl_252 = buffer.data(gl + 252);
    const auto *gl_253 = buffer.data(gl + 253);
    const auto *gl_254 = buffer.data(gl + 254);
    const auto *gl_255 = buffer.data(gl + 255);
    const auto *gl_256 = buffer.data(gl + 256);
    const auto *gl_257 = buffer.data(gl + 257);
    const auto *gl_258 = buffer.data(gl + 258);
    const auto *gl_259 = buffer.data(gl + 259);
    const auto *gl_260 = buffer.data(gl + 260);
    const auto *gl_261 = buffer.data(gl + 261);
    const auto *gl_262 = buffer.data(gl + 262);
    const auto *gl_263 = buffer.data(gl + 263);
    const auto *gl_264 = buffer.data(gl + 264);
    const auto *gl_265 = buffer.data(gl + 265);
    const auto *gl_266 = buffer.data(gl + 266);
    const auto *gl_267 = buffer.data(gl + 267);
    const auto *gl_268 = buffer.data(gl + 268);
    const auto *gl_269 = buffer.data(gl + 269);
    const auto *gl_270 = buffer.data(gl + 270);
    const auto *gl_271 = buffer.data(gl + 271);
    const auto *gl_272 = buffer.data(gl + 272);
    const auto *gl_273 = buffer.data(gl + 273);
    const auto *gl_274 = buffer.data(gl + 274);
    const auto *gl_275 = buffer.data(gl + 275);
    const auto *gl_276 = buffer.data(gl + 276);
    const auto *gl_277 = buffer.data(gl + 277);
    const auto *gl_278 = buffer.data(gl + 278);
    const auto *gl_279 = buffer.data(gl + 279);
    const auto *gl_280 = buffer.data(gl + 280);
    const auto *gl_281 = buffer.data(gl + 281);
    const auto *gl_282 = buffer.data(gl + 282);
    const auto *gl_283 = buffer.data(gl + 283);
    const auto *gl_284 = buffer.data(gl + 284);
    const auto *gl_285 = buffer.data(gl + 285);
    const auto *gl_286 = buffer.data(gl + 286);
    const auto *gl_287 = buffer.data(gl + 287);
    const auto *gl_288 = buffer.data(gl + 288);
    const auto *gl_289 = buffer.data(gl + 289);
    const auto *gl_290 = buffer.data(gl + 290);
    const auto *gl_291 = buffer.data(gl + 291);
    const auto *gl_292 = buffer.data(gl + 292);
    const auto *gl_293 = buffer.data(gl + 293);
    const auto *gl_294 = buffer.data(gl + 294);
    const auto *gl_295 = buffer.data(gl + 295);
    const auto *gl_296 = buffer.data(gl + 296);
    const auto *gl_297 = buffer.data(gl + 297);
    const auto *gl_298 = buffer.data(gl + 298);
    const auto *gl_299 = buffer.data(gl + 299);
    const auto *gl_300 = buffer.data(gl + 300);
    const auto *gl_301 = buffer.data(gl + 301);
    const auto *gl_302 = buffer.data(gl + 302);
    const auto *gl_303 = buffer.data(gl + 303);
    const auto *gl_304 = buffer.data(gl + 304);
    const auto *gl_305 = buffer.data(gl + 305);
    const auto *gl_306 = buffer.data(gl + 306);
    const auto *gl_307 = buffer.data(gl + 307);
    const auto *gl_308 = buffer.data(gl + 308);
    const auto *gl_309 = buffer.data(gl + 309);
    const auto *gl_310 = buffer.data(gl + 310);
    const auto *gl_311 = buffer.data(gl + 311);
    const auto *gl_312 = buffer.data(gl + 312);
    const auto *gl_313 = buffer.data(gl + 313);
    const auto *gl_314 = buffer.data(gl + 314);
    const auto *gl_315 = buffer.data(gl + 315);
    const auto *gl_316 = buffer.data(gl + 316);
    const auto *gl_317 = buffer.data(gl + 317);
    const auto *gl_318 = buffer.data(gl + 318);
    const auto *gl_319 = buffer.data(gl + 319);
    const auto *gl_320 = buffer.data(gl + 320);
    const auto *gl_321 = buffer.data(gl + 321);
    const auto *gl_322 = buffer.data(gl + 322);
    const auto *gl_323 = buffer.data(gl + 323);
    const auto *gl_324 = buffer.data(gl + 324);
    const auto *gl_325 = buffer.data(gl + 325);
    const auto *gl_326 = buffer.data(gl + 326);
    const auto *gl_327 = buffer.data(gl + 327);
    const auto *gl_328 = buffer.data(gl + 328);
    const auto *gl_329 = buffer.data(gl + 329);
    const auto *gl_330 = buffer.data(gl + 330);
    const auto *gl_331 = buffer.data(gl + 331);
    const auto *gl_332 = buffer.data(gl + 332);
    const auto *gl_333 = buffer.data(gl + 333);
    const auto *gl_334 = buffer.data(gl + 334);
    const auto *gl_335 = buffer.data(gl + 335);
    const auto *gl_336 = buffer.data(gl + 336);
    const auto *gl_337 = buffer.data(gl + 337);
    const auto *gl_338 = buffer.data(gl + 338);
    const auto *gl_339 = buffer.data(gl + 339);
    const auto *gl_340 = buffer.data(gl + 340);
    const auto *gl_341 = buffer.data(gl + 341);
    const auto *gl_342 = buffer.data(gl + 342);
    const auto *gl_343 = buffer.data(gl + 343);
    const auto *gl_344 = buffer.data(gl + 344);
    const auto *gl_345 = buffer.data(gl + 345);
    const auto *gl_346 = buffer.data(gl + 346);
    const auto *gl_347 = buffer.data(gl + 347);
    const auto *gl_348 = buffer.data(gl + 348);
    const auto *gl_349 = buffer.data(gl + 349);
    const auto *gl_350 = buffer.data(gl + 350);
    const auto *gl_351 = buffer.data(gl + 351);
    const auto *gl_352 = buffer.data(gl + 352);
    const auto *gl_353 = buffer.data(gl + 353);
    const auto *gl_354 = buffer.data(gl + 354);
    const auto *gl_355 = buffer.data(gl + 355);
    const auto *gl_356 = buffer.data(gl + 356);
    const auto *gl_357 = buffer.data(gl + 357);
    const auto *gl_358 = buffer.data(gl + 358);
    const auto *gl_359 = buffer.data(gl + 359);
    const auto *gl_360 = buffer.data(gl + 360);
    const auto *gl_361 = buffer.data(gl + 361);
    const auto *gl_362 = buffer.data(gl + 362);
    const auto *gl_363 = buffer.data(gl + 363);
    const auto *gl_364 = buffer.data(gl + 364);
    const auto *gl_365 = buffer.data(gl + 365);
    const auto *gl_366 = buffer.data(gl + 366);
    const auto *gl_367 = buffer.data(gl + 367);
    const auto *gl_368 = buffer.data(gl + 368);
    const auto *gl_369 = buffer.data(gl + 369);
    const auto *gl_370 = buffer.data(gl + 370);
    const auto *gl_371 = buffer.data(gl + 371);
    const auto *gl_372 = buffer.data(gl + 372);
    const auto *gl_373 = buffer.data(gl + 373);
    const auto *gl_374 = buffer.data(gl + 374);
    const auto *gl_375 = buffer.data(gl + 375);
    const auto *gl_376 = buffer.data(gl + 376);
    const auto *gl_377 = buffer.data(gl + 377);
    const auto *gl_378 = buffer.data(gl + 378);
    const auto *gl_379 = buffer.data(gl + 379);
    const auto *gl_380 = buffer.data(gl + 380);
    const auto *gl_381 = buffer.data(gl + 381);
    const auto *gl_382 = buffer.data(gl + 382);
    const auto *gl_383 = buffer.data(gl + 383);
    const auto *gl_384 = buffer.data(gl + 384);
    const auto *gl_385 = buffer.data(gl + 385);
    const auto *gl_386 = buffer.data(gl + 386);
    const auto *gl_387 = buffer.data(gl + 387);
    const auto *gl_388 = buffer.data(gl + 388);
    const auto *gl_389 = buffer.data(gl + 389);
    const auto *gl_390 = buffer.data(gl + 390);
    const auto *gl_391 = buffer.data(gl + 391);
    const auto *gl_392 = buffer.data(gl + 392);
    const auto *gl_393 = buffer.data(gl + 393);
    const auto *gl_394 = buffer.data(gl + 394);
    const auto *gl_395 = buffer.data(gl + 395);
    const auto *gl_396 = buffer.data(gl + 396);
    const auto *gl_397 = buffer.data(gl + 397);
    const auto *gl_398 = buffer.data(gl + 398);
    const auto *gl_399 = buffer.data(gl + 399);
    const auto *gl_400 = buffer.data(gl + 400);
    const auto *gl_401 = buffer.data(gl + 401);
    const auto *gl_402 = buffer.data(gl + 402);
    const auto *gl_403 = buffer.data(gl + 403);
    const auto *gl_404 = buffer.data(gl + 404);
    const auto *gl_405 = buffer.data(gl + 405);
    const auto *gl_406 = buffer.data(gl + 406);
    const auto *gl_407 = buffer.data(gl + 407);
    const auto *gl_408 = buffer.data(gl + 408);
    const auto *gl_409 = buffer.data(gl + 409);
    const auto *gl_410 = buffer.data(gl + 410);
    const auto *gl_411 = buffer.data(gl + 411);
    const auto *gl_412 = buffer.data(gl + 412);
    const auto *gl_413 = buffer.data(gl + 413);
    const auto *gl_414 = buffer.data(gl + 414);
    const auto *gl_415 = buffer.data(gl + 415);
    const auto *gl_416 = buffer.data(gl + 416);
    const auto *gl_417 = buffer.data(gl + 417);
    const auto *gl_418 = buffer.data(gl + 418);
    const auto *gl_419 = buffer.data(gl + 419);
    const auto *gl_420 = buffer.data(gl + 420);
    const auto *gl_421 = buffer.data(gl + 421);
    const auto *gl_422 = buffer.data(gl + 422);
    const auto *gl_423 = buffer.data(gl + 423);
    const auto *gl_424 = buffer.data(gl + 424);
    const auto *gl_425 = buffer.data(gl + 425);
    const auto *gl_426 = buffer.data(gl + 426);
    const auto *gl_427 = buffer.data(gl + 427);
    const auto *gl_428 = buffer.data(gl + 428);
    const auto *gl_429 = buffer.data(gl + 429);
    const auto *gl_430 = buffer.data(gl + 430);
    const auto *gl_431 = buffer.data(gl + 431);
    const auto *gl_432 = buffer.data(gl + 432);
    const auto *gl_433 = buffer.data(gl + 433);
    const auto *gl_434 = buffer.data(gl + 434);
    const auto *gl_435 = buffer.data(gl + 435);
    const auto *gl_436 = buffer.data(gl + 436);
    const auto *gl_437 = buffer.data(gl + 437);
    const auto *gl_438 = buffer.data(gl + 438);
    const auto *gl_439 = buffer.data(gl + 439);
    const auto *gl_440 = buffer.data(gl + 440);
    const auto *gl_441 = buffer.data(gl + 441);
    const auto *gl_442 = buffer.data(gl + 442);
    const auto *gl_443 = buffer.data(gl + 443);
    const auto *gl_444 = buffer.data(gl + 444);
    const auto *gl_445 = buffer.data(gl + 445);
    const auto *gl_446 = buffer.data(gl + 446);
    const auto *gl_447 = buffer.data(gl + 447);
    const auto *gl_448 = buffer.data(gl + 448);
    const auto *gl_449 = buffer.data(gl + 449);
    const auto *gl_450 = buffer.data(gl + 450);
    const auto *gl_451 = buffer.data(gl + 451);
    const auto *gl_452 = buffer.data(gl + 452);
    const auto *gl_453 = buffer.data(gl + 453);
    const auto *gl_454 = buffer.data(gl + 454);
    const auto *gl_455 = buffer.data(gl + 455);
    const auto *gl_456 = buffer.data(gl + 456);
    const auto *gl_457 = buffer.data(gl + 457);
    const auto *gl_458 = buffer.data(gl + 458);
    const auto *gl_459 = buffer.data(gl + 459);
    const auto *gl_460 = buffer.data(gl + 460);
    const auto *gl_461 = buffer.data(gl + 461);
    const auto *gl_462 = buffer.data(gl + 462);
    const auto *gl_463 = buffer.data(gl + 463);
    const auto *gl_464 = buffer.data(gl + 464);
    const auto *gl_465 = buffer.data(gl + 465);
    const auto *gl_466 = buffer.data(gl + 466);
    const auto *gl_467 = buffer.data(gl + 467);
    const auto *gl_468 = buffer.data(gl + 468);
    const auto *gl_469 = buffer.data(gl + 469);
    const auto *gl_470 = buffer.data(gl + 470);
    const auto *gl_471 = buffer.data(gl + 471);
    const auto *gl_472 = buffer.data(gl + 472);
    const auto *gl_473 = buffer.data(gl + 473);
    const auto *gl_474 = buffer.data(gl + 474);
    const auto *gl_475 = buffer.data(gl + 475);
    const auto *gl_476 = buffer.data(gl + 476);
    const auto *gl_477 = buffer.data(gl + 477);
    const auto *gl_478 = buffer.data(gl + 478);
    const auto *gl_479 = buffer.data(gl + 479);
    const auto *gl_480 = buffer.data(gl + 480);
    const auto *gl_481 = buffer.data(gl + 481);
    const auto *gl_482 = buffer.data(gl + 482);
    const auto *gl_483 = buffer.data(gl + 483);
    const auto *gl_484 = buffer.data(gl + 484);
    const auto *gl_485 = buffer.data(gl + 485);
    const auto *gl_486 = buffer.data(gl + 486);
    const auto *gl_487 = buffer.data(gl + 487);
    const auto *gl_488 = buffer.data(gl + 488);
    const auto *gl_489 = buffer.data(gl + 489);
    const auto *gl_490 = buffer.data(gl + 490);
    const auto *gl_491 = buffer.data(gl + 491);
    const auto *gl_492 = buffer.data(gl + 492);
    const auto *gl_493 = buffer.data(gl + 493);
    const auto *gl_494 = buffer.data(gl + 494);
    const auto *gl_495 = buffer.data(gl + 495);
    const auto *gl_496 = buffer.data(gl + 496);
    const auto *gl_497 = buffer.data(gl + 497);
    const auto *gl_498 = buffer.data(gl + 498);
    const auto *gl_499 = buffer.data(gl + 499);
    const auto *gl_500 = buffer.data(gl + 500);
    const auto *gl_501 = buffer.data(gl + 501);
    const auto *gl_502 = buffer.data(gl + 502);
    const auto *gl_503 = buffer.data(gl + 503);
    const auto *gl_504 = buffer.data(gl + 504);
    const auto *gl_505 = buffer.data(gl + 505);
    const auto *gl_506 = buffer.data(gl + 506);
    const auto *gl_507 = buffer.data(gl + 507);
    const auto *gl_508 = buffer.data(gl + 508);
    const auto *gl_509 = buffer.data(gl + 509);
    const auto *gl_510 = buffer.data(gl + 510);
    const auto *gl_511 = buffer.data(gl + 511);
    const auto *gl_512 = buffer.data(gl + 512);
    const auto *gl_513 = buffer.data(gl + 513);
    const auto *gl_514 = buffer.data(gl + 514);
    const auto *gl_515 = buffer.data(gl + 515);
    const auto *gl_516 = buffer.data(gl + 516);
    const auto *gl_517 = buffer.data(gl + 517);
    const auto *gl_518 = buffer.data(gl + 518);
    const auto *gl_519 = buffer.data(gl + 519);
    const auto *gl_520 = buffer.data(gl + 520);
    const auto *gl_521 = buffer.data(gl + 521);
    const auto *gl_522 = buffer.data(gl + 522);
    const auto *gl_523 = buffer.data(gl + 523);
    const auto *gl_524 = buffer.data(gl + 524);
    const auto *gl_525 = buffer.data(gl + 525);
    const auto *gl_526 = buffer.data(gl + 526);
    const auto *gl_527 = buffer.data(gl + 527);
    const auto *gl_528 = buffer.data(gl + 528);
    const auto *gl_529 = buffer.data(gl + 529);
    const auto *gl_530 = buffer.data(gl + 530);
    const auto *gl_531 = buffer.data(gl + 531);
    const auto *gl_532 = buffer.data(gl + 532);
    const auto *gl_533 = buffer.data(gl + 533);
    const auto *gl_534 = buffer.data(gl + 534);
    const auto *gl_535 = buffer.data(gl + 535);
    const auto *gl_536 = buffer.data(gl + 536);
    const auto *gl_537 = buffer.data(gl + 537);
    const auto *gl_538 = buffer.data(gl + 538);
    const auto *gl_539 = buffer.data(gl + 539);
    const auto *gl_540 = buffer.data(gl + 540);
    const auto *gl_541 = buffer.data(gl + 541);
    const auto *gl_542 = buffer.data(gl + 542);
    const auto *gl_543 = buffer.data(gl + 543);
    const auto *gl_544 = buffer.data(gl + 544);
    const auto *gl_545 = buffer.data(gl + 545);
    const auto *gl_546 = buffer.data(gl + 546);
    const auto *gl_547 = buffer.data(gl + 547);
    const auto *gl_548 = buffer.data(gl + 548);
    const auto *gl_549 = buffer.data(gl + 549);
    const auto *gl_550 = buffer.data(gl + 550);
    const auto *gl_551 = buffer.data(gl + 551);
    const auto *gl_552 = buffer.data(gl + 552);
    const auto *gl_553 = buffer.data(gl + 553);
    const auto *gl_554 = buffer.data(gl + 554);
    const auto *gl_555 = buffer.data(gl + 555);
    const auto *gl_556 = buffer.data(gl + 556);
    const auto *gl_557 = buffer.data(gl + 557);
    const auto *gl_558 = buffer.data(gl + 558);
    const auto *gl_559 = buffer.data(gl + 559);
    const auto *gl_560 = buffer.data(gl + 560);
    const auto *gl_561 = buffer.data(gl + 561);
    const auto *gl_562 = buffer.data(gl + 562);
    const auto *gl_563 = buffer.data(gl + 563);
    const auto *gl_564 = buffer.data(gl + 564);
    const auto *gl_565 = buffer.data(gl + 565);
    const auto *gl_566 = buffer.data(gl + 566);
    const auto *gl_567 = buffer.data(gl + 567);
    const auto *gl_568 = buffer.data(gl + 568);
    const auto *gl_569 = buffer.data(gl + 569);
    const auto *gl_570 = buffer.data(gl + 570);
    const auto *gl_571 = buffer.data(gl + 571);
    const auto *gl_572 = buffer.data(gl + 572);
    const auto *gl_573 = buffer.data(gl + 573);
    const auto *gl_574 = buffer.data(gl + 574);
    const auto *gl_575 = buffer.data(gl + 575);
    const auto *gl_576 = buffer.data(gl + 576);
    const auto *gl_577 = buffer.data(gl + 577);
    const auto *gl_578 = buffer.data(gl + 578);
    const auto *gl_579 = buffer.data(gl + 579);
    const auto *gl_580 = buffer.data(gl + 580);
    const auto *gl_581 = buffer.data(gl + 581);
    const auto *gl_582 = buffer.data(gl + 582);
    const auto *gl_583 = buffer.data(gl + 583);
    const auto *gl_584 = buffer.data(gl + 584);
    const auto *gl_585 = buffer.data(gl + 585);
    const auto *gl_586 = buffer.data(gl + 586);
    const auto *gl_587 = buffer.data(gl + 587);
    const auto *gl_588 = buffer.data(gl + 588);
    const auto *gl_589 = buffer.data(gl + 589);
    const auto *gl_590 = buffer.data(gl + 590);
    const auto *gl_591 = buffer.data(gl + 591);
    const auto *gl_592 = buffer.data(gl + 592);
    const auto *gl_593 = buffer.data(gl + 593);
    const auto *gl_594 = buffer.data(gl + 594);
    const auto *gl_595 = buffer.data(gl + 595);
    const auto *gl_596 = buffer.data(gl + 596);
    const auto *gl_597 = buffer.data(gl + 597);
    const auto *gl_598 = buffer.data(gl + 598);
    const auto *gl_599 = buffer.data(gl + 599);
    const auto *gl_600 = buffer.data(gl + 600);
    const auto *gl_601 = buffer.data(gl + 601);
    const auto *gl_602 = buffer.data(gl + 602);
    const auto *gl_603 = buffer.data(gl + 603);
    const auto *gl_604 = buffer.data(gl + 604);
    const auto *gl_605 = buffer.data(gl + 605);
    const auto *gl_606 = buffer.data(gl + 606);
    const auto *gl_607 = buffer.data(gl + 607);
    const auto *gl_608 = buffer.data(gl + 608);
    const auto *gl_609 = buffer.data(gl + 609);
    const auto *gl_610 = buffer.data(gl + 610);
    const auto *gl_611 = buffer.data(gl + 611);
    const auto *gl_612 = buffer.data(gl + 612);
    const auto *gl_613 = buffer.data(gl + 613);
    const auto *gl_614 = buffer.data(gl + 614);
    const auto *gl_615 = buffer.data(gl + 615);
    const auto *gl_616 = buffer.data(gl + 616);
    const auto *gl_617 = buffer.data(gl + 617);
    const auto *gl_618 = buffer.data(gl + 618);
    const auto *gl_619 = buffer.data(gl + 619);
    const auto *gl_620 = buffer.data(gl + 620);
    const auto *gl_621 = buffer.data(gl + 621);
    const auto *gl_622 = buffer.data(gl + 622);
    const auto *gl_623 = buffer.data(gl + 623);
    const auto *gl_624 = buffer.data(gl + 624);
    const auto *gl_625 = buffer.data(gl + 625);
    const auto *gl_626 = buffer.data(gl + 626);
    const auto *gl_627 = buffer.data(gl + 627);
    const auto *gl_628 = buffer.data(gl + 628);
    const auto *gl_629 = buffer.data(gl + 629);
    const auto *gl_630 = buffer.data(gl + 630);
    const auto *gl_631 = buffer.data(gl + 631);
    const auto *gl_632 = buffer.data(gl + 632);
    const auto *gl_633 = buffer.data(gl + 633);
    const auto *gl_634 = buffer.data(gl + 634);
    const auto *gl_635 = buffer.data(gl + 635);
    const auto *gl_636 = buffer.data(gl + 636);
    const auto *gl_637 = buffer.data(gl + 637);
    const auto *gl_638 = buffer.data(gl + 638);
    const auto *gl_639 = buffer.data(gl + 639);
    const auto *gl_640 = buffer.data(gl + 640);
    const auto *gl_641 = buffer.data(gl + 641);
    const auto *gl_642 = buffer.data(gl + 642);
    const auto *gl_643 = buffer.data(gl + 643);
    const auto *gl_644 = buffer.data(gl + 644);
    const auto *gl_645 = buffer.data(gl + 645);
    const auto *gl_646 = buffer.data(gl + 646);
    const auto *gl_647 = buffer.data(gl + 647);
    const auto *gl_648 = buffer.data(gl + 648);
    const auto *gl_649 = buffer.data(gl + 649);
    const auto *gl_650 = buffer.data(gl + 650);
    const auto *gl_651 = buffer.data(gl + 651);
    const auto *gl_652 = buffer.data(gl + 652);
    const auto *gl_653 = buffer.data(gl + 653);
    const auto *gl_654 = buffer.data(gl + 654);
    const auto *gl_655 = buffer.data(gl + 655);
    const auto *gl_656 = buffer.data(gl + 656);
    const auto *gl_657 = buffer.data(gl + 657);
    const auto *gl_658 = buffer.data(gl + 658);
    const auto *gl_659 = buffer.data(gl + 659);
    const auto *gl_660 = buffer.data(gl + 660);
    const auto *gl_661 = buffer.data(gl + 661);
    const auto *gl_662 = buffer.data(gl + 662);
    const auto *gl_663 = buffer.data(gl + 663);
    const auto *gl_664 = buffer.data(gl + 664);
    const auto *gl_665 = buffer.data(gl + 665);
    const auto *gl_666 = buffer.data(gl + 666);
    const auto *gl_667 = buffer.data(gl + 667);
    const auto *gl_668 = buffer.data(gl + 668);
    const auto *gl_669 = buffer.data(gl + 669);
    const auto *gl_670 = buffer.data(gl + 670);
    const auto *gl_671 = buffer.data(gl + 671);
    const auto *gl_672 = buffer.data(gl + 672);
    const auto *gl_673 = buffer.data(gl + 673);
    const auto *gl_674 = buffer.data(gl + 674);

#pragma omp simd aligned(gl_46, gl_51, gl_60, gl_73, gl_271, gl_276, gl_285, \
                         gl_298 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * gl_46[k]
                 - f_1 * gl_51[k]
                 + f_1 * gl_60[k]
                 - f_0 * gl_73[k]
                 - f_0 * gl_271[k]
                 + f_1 * gl_276[k]
                 - f_1 * gl_285[k]
                 + f_0 * gl_298[k];
    }

#pragma omp simd aligned(gl_49, gl_56, gl_67, gl_82, gl_274, gl_281, gl_292, \
                         gl_307 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_1[k] = f_2 * gl_49[k]
                 - f_3 * gl_56[k]
                 + f_4 * gl_67[k]
                 - f_5 * gl_82[k]
                 - f_2 * gl_274[k]
                 + f_3 * gl_281[k]
                 - f_4 * gl_292[k]
                 + f_5 * gl_307[k];
    }

#pragma omp simd aligned(gl_46, gl_51, gl_53, gl_60, gl_62, gl_73, gl_75, gl_271, gl_276, \
                         gl_278, gl_285, gl_287, gl_298, gl_300 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_2[k] = -f_6 * gl_46[k]
                 + f_7 * gl_51[k]
                 + f_8 * gl_53[k]
                 + f_7 * gl_60[k]
                 - f_9 * gl_62[k]
                 - f_6 * gl_73[k]
                 + f_8 * gl_75[k]
                 + f_6 * gl_271[k]
                 - f_7 * gl_276[k]
                 - f_8 * gl_278[k]
                 - f_7 * gl_285[k]
                 + f_9 * gl_287[k]
                 + f_6 * gl_298[k]
                 - f_8 * gl_300[k];
    }

#pragma omp simd aligned(gl_49, gl_56, gl_58, gl_67, gl_69, gl_82, gl_84, gl_274, gl_281, \
                         gl_283, gl_292, gl_294, gl_307, gl_309 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = -f_10 * gl_49[k]
                 + f_10 * gl_56[k]
                 + f_11 * gl_58[k]
                 + f_12 * gl_67[k]
                 - f_13 * gl_69[k]
                 - f_14 * gl_82[k]
                 + f_15 * gl_84[k]
                 + f_10 * gl_274[k]
                 - f_10 * gl_281[k]
                 - f_11 * gl_283[k]
                 - f_12 * gl_292[k]
                 + f_13 * gl_294[k]
                 + f_14 * gl_307[k]
                 - f_15 * gl_309[k];
    }

#pragma omp simd aligned(gl_46, gl_51, gl_53, gl_60, gl_64, gl_73, gl_75, gl_77, gl_271, \
                         gl_276, gl_278, gl_285, gl_289, gl_298, gl_300, \
                         gl_302 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = f_16 * gl_46[k]
                 + f_16 * gl_51[k]
                 - f_17 * gl_53[k]
                 - f_16 * gl_60[k]
                 + f_18 * gl_64[k]
                 - f_16 * gl_73[k]
                 + f_17 * gl_75[k]
                 - f_18 * gl_77[k]
                 - f_16 * gl_271[k]
                 - f_16 * gl_276[k]
                 + f_17 * gl_278[k]
                 + f_16 * gl_285[k]
                 - f_18 * gl_289[k]
                 + f_16 * gl_298[k]
                 - f_17 * gl_300[k]
                 + f_18 * gl_302[k];
    }

#pragma omp simd aligned(gl_49, gl_56, gl_58, gl_67, gl_69, gl_71, gl_82, gl_84, gl_86, \
                         gl_274, gl_281, gl_283, gl_292, gl_294, gl_296, gl_307, gl_309, \
                         gl_311 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = f_19 * gl_49[k]
                 + f_20 * gl_56[k]
                 - f_21 * gl_58[k]
                 + f_22 * gl_67[k]
                 - f_23 * gl_69[k]
                 + f_24 * gl_71[k]
                 - f_22 * gl_82[k]
                 + f_25 * gl_84[k]
                 - f_26 * gl_86[k]
                 - f_19 * gl_274[k]
                 - f_20 * gl_281[k]
                 + f_21 * gl_283[k]
                 - f_22 * gl_292[k]
                 + f_23 * gl_294[k]
                 - f_24 * gl_296[k]
                 + f_22 * gl_307[k]
                 - f_25 * gl_309[k]
                 + f_26 * gl_311[k];
    }

#pragma omp simd aligned(gl_46, gl_51, gl_53, gl_60, gl_62, gl_64, gl_73, gl_75, gl_77, gl_79, \
                         gl_271, gl_276, gl_278, gl_285, gl_287, gl_289, gl_298, gl_300, \
                         gl_302, gl_304 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = -f_27 * gl_46[k]
                 - f_28 * gl_51[k]
                 + f_29 * gl_53[k]
                 - f_28 * gl_60[k]
                 + f_30 * gl_62[k]
                 - f_31 * gl_64[k]
                 - f_27 * gl_73[k]
                 + f_29 * gl_75[k]
                 - f_31 * gl_77[k]
                 + f_32 * gl_79[k]
                 + f_27 * gl_271[k]
                 + f_28 * gl_276[k]
                 - f_29 * gl_278[k]
                 + f_28 * gl_285[k]
                 - f_30 * gl_287[k]
                 + f_31 * gl_289[k]
                 + f_27 * gl_298[k]
                 - f_29 * gl_300[k]
                 + f_31 * gl_302[k]
                 - f_32 * gl_304[k];
    }

#pragma omp simd aligned(gl_49, gl_56, gl_58, gl_67, gl_69, gl_71, gl_82, gl_84, gl_86, gl_88, \
                         gl_274, gl_281, gl_283, gl_292, gl_294, gl_296, gl_307, gl_309, \
                         gl_311, gl_313 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = -f_33 * gl_49[k]
                 - f_34 * gl_56[k]
                 + f_35 * gl_58[k]
                 - f_34 * gl_67[k]
                 + f_36 * gl_69[k]
                 - f_37 * gl_71[k]
                 - f_33 * gl_82[k]
                 + f_35 * gl_84[k]
                 - f_37 * gl_86[k]
                 + f_38 * gl_88[k]
                 + f_33 * gl_274[k]
                 + f_34 * gl_281[k]
                 - f_35 * gl_283[k]
                 + f_34 * gl_292[k]
                 - f_36 * gl_294[k]
                 + f_37 * gl_296[k]
                 + f_33 * gl_307[k]
                 - f_35 * gl_309[k]
                 + f_37 * gl_311[k]
                 - f_38 * gl_313[k];
    }

#pragma omp simd aligned(gl_45, gl_48, gl_50, gl_55, gl_57, gl_59, gl_66, gl_68, gl_70, gl_72, \
                         gl_81, gl_83, gl_85, gl_87, gl_89, gl_270, gl_273, gl_275, gl_280, \
                         gl_282, gl_284, gl_291, gl_293, gl_295, gl_297, gl_306, gl_308, \
                         gl_310, gl_312, gl_314 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_8[k] = f_39 * gl_45[k]
                 + f_40 * gl_48[k]
                 - f_41 * gl_50[k]
                 + f_42 * gl_55[k]
                 - f_35 * gl_57[k]
                 + f_35 * gl_59[k]
                 + f_40 * gl_66[k]
                 - f_35 * gl_68[k]
                 + f_36 * gl_70[k]
                 - f_43 * gl_72[k]
                 + f_39 * gl_81[k]
                 - f_41 * gl_83[k]
                 + f_35 * gl_85[k]
                 - f_43 * gl_87[k]
                 + f_44 * gl_89[k]
                 - f_39 * gl_270[k]
                 - f_40 * gl_273[k]
                 + f_41 * gl_275[k]
                 - f_42 * gl_280[k]
                 + f_35 * gl_282[k]
                 - f_35 * gl_284[k]
                 - f_40 * gl_291[k]
                 + f_35 * gl_293[k]
                 - f_36 * gl_295[k]
                 + f_43 * gl_297[k]
                 - f_39 * gl_306[k]
                 + f_41 * gl_308[k]
                 - f_35 * gl_310[k]
                 + f_43 * gl_312[k]
                 - f_44 * gl_314[k];
    }

#pragma omp simd aligned(gl_47, gl_52, gl_54, gl_61, gl_63, gl_65, gl_74, gl_76, gl_78, gl_80, \
                         gl_272, gl_277, gl_279, gl_286, gl_288, gl_290, gl_299, gl_301, \
                         gl_303, gl_305 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_9[k] = -f_33 * gl_47[k]
                 - f_34 * gl_52[k]
                 + f_35 * gl_54[k]
                 - f_34 * gl_61[k]
                 + f_36 * gl_63[k]
                 - f_37 * gl_65[k]
                 - f_33 * gl_74[k]
                 + f_35 * gl_76[k]
                 - f_37 * gl_78[k]
                 + f_38 * gl_80[k]
                 + f_33 * gl_272[k]
                 + f_34 * gl_277[k]
                 - f_35 * gl_279[k]
                 + f_34 * gl_286[k]
                 - f_36 * gl_288[k]
                 + f_37 * gl_290[k]
                 + f_33 * gl_299[k]
                 - f_35 * gl_301[k]
                 + f_37 * gl_303[k]
                 - f_38 * gl_305[k];
    }

#pragma omp simd aligned(gl_45, gl_48, gl_50, gl_57, gl_59, gl_66, gl_68, gl_72, gl_81, gl_83, \
                         gl_85, gl_87, gl_270, gl_273, gl_275, gl_282, gl_284, gl_291, gl_293, \
                         gl_297, gl_306, gl_308, gl_310, gl_312 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = -f_45 * gl_45[k]
                  - f_27 * gl_48[k]
                  + f_46 * gl_50[k]
                  + f_46 * gl_57[k]
                  - f_47 * gl_59[k]
                  + f_27 * gl_66[k]
                  - f_46 * gl_68[k]
                  + f_48 * gl_72[k]
                  + f_45 * gl_81[k]
                  - f_46 * gl_83[k]
                  + f_47 * gl_85[k]
                  - f_48 * gl_87[k]
                  + f_45 * gl_270[k]
                  + f_27 * gl_273[k]
                  - f_46 * gl_275[k]
                  - f_46 * gl_282[k]
                  + f_47 * gl_284[k]
                  - f_27 * gl_291[k]
                  + f_46 * gl_293[k]
                  - f_48 * gl_297[k]
                  - f_45 * gl_306[k]
                  + f_46 * gl_308[k]
                  - f_47 * gl_310[k]
                  + f_48 * gl_312[k];
    }

#pragma omp simd aligned(gl_47, gl_52, gl_54, gl_61, gl_63, gl_65, gl_74, gl_76, gl_78, \
                         gl_272, gl_277, gl_279, gl_286, gl_288, gl_290, gl_299, gl_301, \
                         gl_303 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_11[k] = f_22 * gl_47[k]
                  - f_22 * gl_52[k]
                  - f_25 * gl_54[k]
                  - f_20 * gl_61[k]
                  + f_23 * gl_63[k]
                  + f_26 * gl_65[k]
                  - f_19 * gl_74[k]
                  + f_21 * gl_76[k]
                  - f_24 * gl_78[k]
                  - f_22 * gl_272[k]
                  + f_22 * gl_277[k]
                  + f_25 * gl_279[k]
                  + f_20 * gl_286[k]
                  - f_23 * gl_288[k]
                  - f_26 * gl_290[k]
                  + f_19 * gl_299[k]
                  - f_21 * gl_301[k]
                  + f_24 * gl_303[k];
    }

#pragma omp simd aligned(gl_45, gl_48, gl_50, gl_55, gl_57, gl_59, gl_66, gl_68, gl_70, gl_81, \
                         gl_83, gl_85, gl_270, gl_273, gl_275, gl_280, gl_282, gl_284, gl_291, \
                         gl_293, gl_295, gl_306, gl_308, gl_310 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_12[k] = f_49 * gl_45[k]
                  - f_16 * gl_48[k]
                  - f_50 * gl_50[k]
                  - f_51 * gl_55[k]
                  + f_52 * gl_57[k]
                  + f_53 * gl_59[k]
                  - f_16 * gl_66[k]
                  + f_52 * gl_68[k]
                  - f_54 * gl_70[k]
                  + f_49 * gl_81[k]
                  - f_50 * gl_83[k]
                  + f_53 * gl_85[k]
                  - f_49 * gl_270[k]
                  + f_16 * gl_273[k]
                  + f_50 * gl_275[k]
                  + f_51 * gl_280[k]
                  - f_52 * gl_282[k]
                  - f_53 * gl_284[k]
                  + f_16 * gl_291[k]
                  - f_52 * gl_293[k]
                  + f_54 * gl_295[k]
                  - f_49 * gl_306[k]
                  + f_50 * gl_308[k]
                  - f_53 * gl_310[k];
    }

#pragma omp simd aligned(gl_47, gl_52, gl_54, gl_61, gl_63, gl_74, gl_76, gl_272, gl_277, \
                         gl_279, gl_286, gl_288, gl_299, gl_301 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = -f_14 * gl_47[k]
                  + f_12 * gl_52[k]
                  + f_15 * gl_54[k]
                  + f_10 * gl_61[k]
                  - f_13 * gl_63[k]
                  - f_10 * gl_74[k]
                  + f_11 * gl_76[k]
                  + f_14 * gl_272[k]
                  - f_12 * gl_277[k]
                  - f_15 * gl_279[k]
                  - f_10 * gl_286[k]
                  + f_13 * gl_288[k]
                  + f_10 * gl_299[k]
                  - f_11 * gl_301[k];
    }

#pragma omp simd aligned(gl_45, gl_48, gl_50, gl_57, gl_66, gl_68, gl_81, gl_83, gl_270, \
                         gl_273, gl_275, gl_282, gl_291, gl_293, gl_306, \
                         gl_308 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = -f_55 * gl_45[k]
                  + f_7 * gl_48[k]
                  + f_7 * gl_50[k]
                  - f_56 * gl_57[k]
                  - f_7 * gl_66[k]
                  + f_56 * gl_68[k]
                  + f_55 * gl_81[k]
                  - f_7 * gl_83[k]
                  + f_55 * gl_270[k]
                  - f_7 * gl_273[k]
                  - f_7 * gl_275[k]
                  + f_56 * gl_282[k]
                  + f_7 * gl_291[k]
                  - f_56 * gl_293[k]
                  - f_55 * gl_306[k]
                  + f_7 * gl_308[k];
    }

#pragma omp simd aligned(gl_47, gl_52, gl_61, gl_74, gl_272, gl_277, gl_286, \
                         gl_299 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_15[k] = f_5 * gl_47[k]
                  - f_4 * gl_52[k]
                  + f_3 * gl_61[k]
                  - f_2 * gl_74[k]
                  - f_5 * gl_272[k]
                  + f_4 * gl_277[k]
                  - f_3 * gl_286[k]
                  + f_2 * gl_299[k];
    }

#pragma omp simd aligned(gl_45, gl_48, gl_55, gl_66, gl_81, gl_270, gl_273, gl_280, gl_291, \
                         gl_306 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_16[k] = f_57 * gl_45[k]
                  - f_2 * gl_48[k]
                  + f_58 * gl_55[k]
                  - f_2 * gl_66[k]
                  + f_57 * gl_81[k]
                  - f_57 * gl_270[k]
                  + f_2 * gl_273[k]
                  - f_58 * gl_280[k]
                  + f_2 * gl_291[k]
                  - f_57 * gl_306[k];
    }

#pragma omp simd aligned(gl_181, gl_186, gl_195, gl_208, gl_496, gl_501, gl_510, \
                         gl_523 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = f_59 * gl_181[k]
                  - f_60 * gl_186[k]
                  + f_60 * gl_195[k]
                  - f_59 * gl_208[k]
                  - f_61 * gl_496[k]
                  + f_62 * gl_501[k]
                  - f_62 * gl_510[k]
                  + f_61 * gl_523[k];
    }

#pragma omp simd aligned(gl_184, gl_191, gl_202, gl_217, gl_499, gl_506, gl_517, \
                         gl_532 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = f_63 * gl_184[k]
                  - f_64 * gl_191[k]
                  + f_65 * gl_202[k]
                  - f_66 * gl_217[k]
                  - f_67 * gl_499[k]
                  + f_68 * gl_506[k]
                  - f_63 * gl_517[k]
                  + f_69 * gl_532[k];
    }

#pragma omp simd aligned(gl_181, gl_186, gl_188, gl_195, gl_197, gl_208, gl_210, gl_496, \
                         gl_501, gl_503, gl_510, gl_512, gl_523, \
                         gl_525 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = -f_70 * gl_181[k]
                  + f_71 * gl_186[k]
                  + f_72 * gl_188[k]
                  + f_71 * gl_195[k]
                  - f_73 * gl_197[k]
                  - f_70 * gl_208[k]
                  + f_72 * gl_210[k]
                  + f_74 * gl_496[k]
                  - f_75 * gl_501[k]
                  - f_76 * gl_503[k]
                  - f_75 * gl_510[k]
                  + f_77 * gl_512[k]
                  + f_74 * gl_523[k]
                  - f_76 * gl_525[k];
    }

#pragma omp simd aligned(gl_184, gl_191, gl_193, gl_202, gl_204, gl_217, gl_219, gl_499, \
                         gl_506, gl_508, gl_517, gl_519, gl_532, \
                         gl_534 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = -f_78 * gl_184[k]
                  + f_78 * gl_191[k]
                  + f_79 * gl_193[k]
                  + f_80 * gl_202[k]
                  - f_81 * gl_204[k]
                  - f_82 * gl_217[k]
                  + f_83 * gl_219[k]
                  + f_84 * gl_499[k]
                  - f_84 * gl_506[k]
                  - f_85 * gl_508[k]
                  - f_86 * gl_517[k]
                  + f_87 * gl_519[k]
                  + f_88 * gl_532[k]
                  - f_89 * gl_534[k];
    }

#pragma omp simd aligned(gl_181, gl_186, gl_188, gl_195, gl_199, gl_208, gl_210, gl_212, \
                         gl_496, gl_501, gl_503, gl_510, gl_514, gl_523, gl_525, \
                         gl_527 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = f_90 * gl_181[k]
                  + f_90 * gl_186[k]
                  - f_91 * gl_188[k]
                  - f_90 * gl_195[k]
                  + f_92 * gl_199[k]
                  - f_90 * gl_208[k]
                  + f_91 * gl_210[k]
                  - f_92 * gl_212[k]
                  - f_93 * gl_496[k]
                  - f_93 * gl_501[k]
                  + f_94 * gl_503[k]
                  + f_93 * gl_510[k]
                  - f_95 * gl_514[k]
                  + f_93 * gl_523[k]
                  - f_94 * gl_525[k]
                  + f_95 * gl_527[k];
    }

#pragma omp simd aligned(gl_184, gl_191, gl_193, gl_202, gl_204, gl_206, gl_217, gl_219, \
                         gl_221, gl_499, gl_506, gl_508, gl_517, gl_519, gl_521, gl_532, \
                         gl_534, gl_536 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = f_96 * gl_184[k]
                  + f_97 * gl_191[k]
                  - f_98 * gl_193[k]
                  + f_99 * gl_202[k]
                  - f_100 * gl_204[k]
                  + f_101 * gl_206[k]
                  - f_99 * gl_217[k]
                  + f_102 * gl_219[k]
                  - f_103 * gl_221[k]
                  - f_99 * gl_499[k]
                  - f_104 * gl_506[k]
                  + f_102 * gl_508[k]
                  - f_105 * gl_517[k]
                  + f_106 * gl_519[k]
                  - f_103 * gl_521[k]
                  + f_105 * gl_532[k]
                  - f_107 * gl_534[k]
                  + f_108 * gl_536[k];
    }

#pragma omp simd aligned(gl_181, gl_186, gl_188, gl_195, gl_197, gl_199, gl_208, gl_210, \
                         gl_212, gl_214, gl_496, gl_501, gl_503, gl_510, gl_512, gl_514, \
                         gl_523, gl_525, gl_527, gl_529 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = -4.921875 * gl_181[k]
                  - 14.765625 * gl_186[k]
                  + 147.65625 * gl_188[k]
                  - 14.765625 * gl_195[k]
                  + 295.3125 * gl_197[k]
                  - 393.75 * gl_199[k]
                  - 4.921875 * gl_208[k]
                  + 147.65625 * gl_210[k]
                  - 393.75 * gl_212[k]
                  + 157.5 * gl_214[k]
                  + 1.640625 * gl_496[k]
                  + 4.921875 * gl_501[k]
                  - 49.21875 * gl_503[k]
                  + 4.921875 * gl_510[k]
                  - 98.4375 * gl_512[k]
                  + 131.25 * gl_514[k]
                  + 1.640625 * gl_523[k]
                  - 49.21875 * gl_525[k]
                  + 131.25 * gl_527[k]
                  - 52.5 * gl_529[k];
    }

#pragma omp simd aligned(gl_184, gl_191, gl_193, gl_202, gl_204, gl_206, gl_217, gl_219, \
                         gl_221, gl_223, gl_499, gl_506, gl_508, gl_517, gl_519, gl_521, \
                         gl_532, gl_534, gl_536, gl_538 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = -f_109 * gl_184[k]
                  - f_110 * gl_191[k]
                  + f_111 * gl_193[k]
                  - f_110 * gl_202[k]
                  + f_112 * gl_204[k]
                  - f_113 * gl_206[k]
                  - f_109 * gl_217[k]
                  + f_111 * gl_219[k]
                  - f_113 * gl_221[k]
                  + f_114 * gl_223[k]
                  + f_115 * gl_499[k]
                  + f_109 * gl_506[k]
                  - f_116 * gl_508[k]
                  + f_109 * gl_517[k]
                  - f_117 * gl_519[k]
                  + f_118 * gl_521[k]
                  + f_115 * gl_532[k]
                  - f_116 * gl_534[k]
                  + f_118 * gl_536[k]
                  - f_119 * gl_538[k];
    }

#pragma omp simd aligned(gl_180, gl_183, gl_185, gl_190, gl_192, gl_194, gl_201, gl_203, \
                         gl_205, gl_207, gl_216, gl_218, gl_220, gl_222, gl_224, gl_495, \
                         gl_498, gl_500, gl_505, gl_507, gl_509, gl_516, gl_518, gl_520, \
                         gl_522, gl_531, gl_533, gl_535, gl_537, \
                         gl_539 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = f_120 * gl_180[k]
                  + f_115 * gl_183[k]
                  - f_116 * gl_185[k]
                  + f_121 * gl_190[k]
                  - f_111 * gl_192[k]
                  + f_111 * gl_194[k]
                  + f_115 * gl_201[k]
                  - f_111 * gl_203[k]
                  + f_112 * gl_205[k]
                  - f_122 * gl_207[k]
                  + f_120 * gl_216[k]
                  - f_116 * gl_218[k]
                  + f_111 * gl_220[k]
                  - f_122 * gl_222[k]
                  + f_123 * gl_224[k]
                  - f_124 * gl_495[k]
                  - f_125 * gl_498[k]
                  + f_126 * gl_500[k]
                  - f_127 * gl_505[k]
                  + f_116 * gl_507[k]
                  - f_116 * gl_509[k]
                  - f_125 * gl_516[k]
                  + f_116 * gl_518[k]
                  - f_117 * gl_520[k]
                  + f_128 * gl_522[k]
                  - f_124 * gl_531[k]
                  + f_126 * gl_533[k]
                  - f_116 * gl_535[k]
                  + f_128 * gl_537[k]
                  - f_129 * gl_539[k];
    }

#pragma omp simd aligned(gl_182, gl_187, gl_189, gl_196, gl_198, gl_200, gl_209, gl_211, \
                         gl_213, gl_215, gl_497, gl_502, gl_504, gl_511, gl_513, gl_515, \
                         gl_524, gl_526, gl_528, gl_530 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = -f_109 * gl_182[k]
                  - f_110 * gl_187[k]
                  + f_111 * gl_189[k]
                  - f_110 * gl_196[k]
                  + f_112 * gl_198[k]
                  - f_113 * gl_200[k]
                  - f_109 * gl_209[k]
                  + f_111 * gl_211[k]
                  - f_113 * gl_213[k]
                  + f_114 * gl_215[k]
                  + f_115 * gl_497[k]
                  + f_109 * gl_502[k]
                  - f_116 * gl_504[k]
                  + f_109 * gl_511[k]
                  - f_117 * gl_513[k]
                  + f_118 * gl_515[k]
                  + f_115 * gl_524[k]
                  - f_116 * gl_526[k]
                  + f_118 * gl_528[k]
                  - f_119 * gl_530[k];
    }

#pragma omp simd aligned(gl_180, gl_183, gl_185, gl_192, gl_194, gl_201, gl_203, gl_207, \
                         gl_216, gl_218, gl_220, gl_222, gl_495, gl_498, gl_500, gl_507, \
                         gl_509, gl_516, gl_518, gl_522, gl_531, gl_533, gl_535, \
                         gl_537 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = -2.4609375 * gl_180[k]
                  - 4.921875 * gl_183[k]
                  + 73.828125 * gl_185[k]
                  + 73.828125 * gl_192[k]
                  - 196.875 * gl_194[k]
                  + 4.921875 * gl_201[k]
                  - 73.828125 * gl_203[k]
                  + 78.75 * gl_207[k]
                  + 2.4609375 * gl_216[k]
                  - 73.828125 * gl_218[k]
                  + 196.875 * gl_220[k]
                  - 78.75 * gl_222[k]
                  + 0.8203125 * gl_495[k]
                  + 1.640625 * gl_498[k]
                  - 24.609375 * gl_500[k]
                  - 24.609375 * gl_507[k]
                  + 65.625 * gl_509[k]
                  - 1.640625 * gl_516[k]
                  + 24.609375 * gl_518[k]
                  - 26.25 * gl_522[k]
                  - 0.8203125 * gl_531[k]
                  + 24.609375 * gl_533[k]
                  - 65.625 * gl_535[k]
                  + 26.25 * gl_537[k];
    }

#pragma omp simd aligned(gl_182, gl_187, gl_189, gl_196, gl_198, gl_200, gl_209, gl_211, \
                         gl_213, gl_497, gl_502, gl_504, gl_511, gl_513, gl_515, gl_524, \
                         gl_526, gl_528 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = f_99 * gl_182[k]
                  - f_99 * gl_187[k]
                  - f_102 * gl_189[k]
                  - f_97 * gl_196[k]
                  + f_100 * gl_198[k]
                  + f_103 * gl_200[k]
                  - f_96 * gl_209[k]
                  + f_98 * gl_211[k]
                  - f_101 * gl_213[k]
                  - f_105 * gl_497[k]
                  + f_105 * gl_502[k]
                  + f_107 * gl_504[k]
                  + f_104 * gl_511[k]
                  - f_106 * gl_513[k]
                  - f_108 * gl_515[k]
                  + f_99 * gl_524[k]
                  - f_102 * gl_526[k]
                  + f_103 * gl_528[k];
    }

#pragma omp simd aligned(gl_180, gl_183, gl_185, gl_190, gl_192, gl_194, gl_201, gl_203, \
                         gl_205, gl_216, gl_218, gl_220, gl_495, gl_498, gl_500, gl_505, \
                         gl_507, gl_509, gl_516, gl_518, gl_520, gl_531, gl_533, \
                         gl_535 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = f_130 * gl_180[k]
                  - f_90 * gl_183[k]
                  - f_131 * gl_185[k]
                  - f_132 * gl_190[k]
                  + f_133 * gl_192[k]
                  + f_134 * gl_194[k]
                  - f_90 * gl_201[k]
                  + f_133 * gl_203[k]
                  - f_135 * gl_205[k]
                  + f_130 * gl_216[k]
                  - f_131 * gl_218[k]
                  + f_134 * gl_220[k]
                  - f_136 * gl_495[k]
                  + f_93 * gl_498[k]
                  + f_137 * gl_500[k]
                  + f_138 * gl_505[k]
                  - f_134 * gl_507[k]
                  - f_139 * gl_509[k]
                  + f_93 * gl_516[k]
                  - f_134 * gl_518[k]
                  + f_140 * gl_520[k]
                  - f_136 * gl_531[k]
                  + f_137 * gl_533[k]
                  - f_139 * gl_535[k];
    }

#pragma omp simd aligned(gl_182, gl_187, gl_189, gl_196, gl_198, gl_209, gl_211, gl_497, \
                         gl_502, gl_504, gl_511, gl_513, gl_524, \
                         gl_526 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = -f_82 * gl_182[k]
                  + f_80 * gl_187[k]
                  + f_83 * gl_189[k]
                  + f_78 * gl_196[k]
                  - f_81 * gl_198[k]
                  - f_78 * gl_209[k]
                  + f_79 * gl_211[k]
                  + f_88 * gl_497[k]
                  - f_86 * gl_502[k]
                  - f_89 * gl_504[k]
                  - f_84 * gl_511[k]
                  + f_87 * gl_513[k]
                  + f_84 * gl_524[k]
                  - f_85 * gl_526[k];
    }

#pragma omp simd aligned(gl_180, gl_183, gl_185, gl_192, gl_201, gl_203, gl_216, gl_218, \
                         gl_495, gl_498, gl_500, gl_507, gl_516, gl_518, gl_531, \
                         gl_533 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = -f_141 * gl_180[k]
                  + f_71 * gl_183[k]
                  + f_71 * gl_185[k]
                  - f_142 * gl_192[k]
                  - f_71 * gl_201[k]
                  + f_142 * gl_203[k]
                  + f_141 * gl_216[k]
                  - f_71 * gl_218[k]
                  + f_143 * gl_495[k]
                  - f_75 * gl_498[k]
                  - f_75 * gl_500[k]
                  + f_144 * gl_507[k]
                  + f_75 * gl_516[k]
                  - f_144 * gl_518[k]
                  - f_143 * gl_531[k]
                  + f_75 * gl_533[k];
    }

#pragma omp simd aligned(gl_182, gl_187, gl_196, gl_209, gl_497, gl_502, gl_511, \
                         gl_524 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = f_66 * gl_182[k]
                  - f_65 * gl_187[k]
                  + f_64 * gl_196[k]
                  - f_63 * gl_209[k]
                  - f_69 * gl_497[k]
                  + f_63 * gl_502[k]
                  - f_68 * gl_511[k]
                  + f_67 * gl_524[k];
    }

#pragma omp simd aligned(gl_180, gl_183, gl_190, gl_201, gl_216, gl_495, gl_498, gl_505, \
                         gl_516, gl_531 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = f_145 * gl_180[k]
                  - f_63 * gl_183[k]
                  + f_146 * gl_190[k]
                  - f_63 * gl_201[k]
                  + f_145 * gl_216[k]
                  - f_147 * gl_495[k]
                  + f_67 * gl_498[k]
                  - f_148 * gl_505[k]
                  + f_67 * gl_516[k]
                  - f_147 * gl_531[k];
    }

#pragma omp simd aligned(gl_46, gl_51, gl_60, gl_73, gl_271, gl_276, gl_285, gl_298, gl_361, \
                         gl_366, gl_375, gl_388 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = -f_149 * gl_46[k]
                  + f_150 * gl_51[k]
                  - f_150 * gl_60[k]
                  + f_149 * gl_73[k]
                  - f_149 * gl_271[k]
                  + f_150 * gl_276[k]
                  - f_150 * gl_285[k]
                  + f_149 * gl_298[k]
                  + f_151 * gl_361[k]
                  - f_152 * gl_366[k]
                  + f_152 * gl_375[k]
                  - f_151 * gl_388[k];
    }

#pragma omp simd aligned(gl_49, gl_56, gl_67, gl_82, gl_274, gl_281, gl_292, gl_307, gl_364, \
                         gl_371, gl_382, gl_397 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = -f_153 * gl_49[k]
                  + f_154 * gl_56[k]
                  - f_155 * gl_67[k]
                  + f_156 * gl_82[k]
                  - f_153 * gl_274[k]
                  + f_154 * gl_281[k]
                  - f_155 * gl_292[k]
                  + f_156 * gl_307[k]
                  + f_157 * gl_364[k]
                  - f_158 * gl_371[k]
                  + f_159 * gl_382[k]
                  - f_160 * gl_397[k];
    }

#pragma omp simd aligned(gl_46, gl_51, gl_53, gl_60, gl_62, gl_73, gl_75, gl_271, gl_276, \
                         gl_278, gl_285, gl_287, gl_298, gl_300, gl_361, gl_366, gl_368, \
                         gl_375, gl_377, gl_388, gl_390 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_36[k] = f_161 * gl_46[k]
                  - f_162 * gl_51[k]
                  - f_163 * gl_53[k]
                  - f_162 * gl_60[k]
                  + f_164 * gl_62[k]
                  + f_161 * gl_73[k]
                  - f_163 * gl_75[k]
                  + f_161 * gl_271[k]
                  - f_162 * gl_276[k]
                  - f_163 * gl_278[k]
                  - f_162 * gl_285[k]
                  + f_164 * gl_287[k]
                  + f_161 * gl_298[k]
                  - f_163 * gl_300[k]
                  - f_165 * gl_361[k]
                  + f_163 * gl_366[k]
                  + f_166 * gl_368[k]
                  + f_163 * gl_375[k]
                  - f_167 * gl_377[k]
                  - f_165 * gl_388[k]
                  + f_166 * gl_390[k];
    }

#pragma omp simd aligned(gl_49, gl_56, gl_58, gl_67, gl_69, gl_82, gl_84, gl_274, gl_281, \
                         gl_283, gl_292, gl_294, gl_307, gl_309, gl_364, gl_371, gl_373, \
                         gl_382, gl_384, gl_397, gl_399 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_37[k] = f_168 * gl_49[k]
                  - f_168 * gl_56[k]
                  - f_169 * gl_58[k]
                  - f_170 * gl_67[k]
                  + f_171 * gl_69[k]
                  + f_172 * gl_82[k]
                  - f_173 * gl_84[k]
                  + f_168 * gl_274[k]
                  - f_168 * gl_281[k]
                  - f_169 * gl_283[k]
                  - f_170 * gl_292[k]
                  + f_171 * gl_294[k]
                  + f_172 * gl_307[k]
                  - f_173 * gl_309[k]
                  - f_174 * gl_364[k]
                  + f_174 * gl_371[k]
                  + f_175 * gl_373[k]
                  + f_176 * gl_382[k]
                  - f_177 * gl_384[k]
                  - f_178 * gl_397[k]
                  + f_179 * gl_399[k];
    }

#pragma omp simd aligned(gl_46, gl_51, gl_53, gl_60, gl_64, gl_73, gl_75, gl_77, gl_271, \
                         gl_276, gl_278, gl_285, gl_289, gl_298, gl_300, gl_302, gl_361, \
                         gl_366, gl_368, gl_375, gl_379, gl_388, gl_390, \
                         gl_392 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = -f_180 * gl_46[k]
                  - f_180 * gl_51[k]
                  + f_181 * gl_53[k]
                  + f_180 * gl_60[k]
                  - f_182 * gl_64[k]
                  + f_180 * gl_73[k]
                  - f_181 * gl_75[k]
                  + f_182 * gl_77[k]
                  - f_180 * gl_271[k]
                  - f_180 * gl_276[k]
                  + f_181 * gl_278[k]
                  + f_180 * gl_285[k]
                  - f_182 * gl_289[k]
                  + f_180 * gl_298[k]
                  - f_181 * gl_300[k]
                  + f_182 * gl_302[k]
                  + f_183 * gl_361[k]
                  + f_183 * gl_366[k]
                  - f_184 * gl_368[k]
                  - f_183 * gl_375[k]
                  + f_185 * gl_379[k]
                  - f_183 * gl_388[k]
                  + f_184 * gl_390[k]
                  - f_185 * gl_392[k];
    }

#pragma omp simd aligned(gl_49, gl_56, gl_58, gl_67, gl_69, gl_71, gl_82, gl_84, gl_86, \
                         gl_274, gl_281, gl_283, gl_292, gl_294, gl_296, gl_307, gl_309, \
                         gl_311, gl_364, gl_371, gl_373, gl_382, gl_384, gl_386, gl_397, \
                         gl_399, gl_401 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_39[k] = -f_186 * gl_49[k]
                  - f_187 * gl_56[k]
                  + f_188 * gl_58[k]
                  - f_189 * gl_67[k]
                  + f_190 * gl_69[k]
                  - f_191 * gl_71[k]
                  + f_189 * gl_82[k]
                  - f_192 * gl_84[k]
                  + f_193 * gl_86[k]
                  - f_186 * gl_274[k]
                  - f_187 * gl_281[k]
                  + f_188 * gl_283[k]
                  - f_189 * gl_292[k]
                  + f_190 * gl_294[k]
                  - f_191 * gl_296[k]
                  + f_189 * gl_307[k]
                  - f_192 * gl_309[k]
                  + f_193 * gl_311[k]
                  + f_194 * gl_364[k]
                  + f_195 * gl_371[k]
                  - f_196 * gl_373[k]
                  + f_197 * gl_382[k]
                  - f_198 * gl_384[k]
                  + f_199 * gl_386[k]
                  - f_197 * gl_397[k]
                  + f_200 * gl_399[k]
                  - f_201 * gl_401[k];
    }

#pragma omp simd aligned(gl_46, gl_51, gl_53, gl_60, gl_62, gl_64, gl_73, gl_75, gl_77, gl_79, \
                         gl_271, gl_276, gl_278, gl_285, gl_287, gl_289, gl_298, gl_300, \
                         gl_302, gl_304, gl_361, gl_366, gl_368, gl_375, gl_377, gl_379, \
                         gl_388, gl_390, gl_392, gl_394 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_40[k] = f_202 * gl_46[k]
                  + f_203 * gl_51[k]
                  - f_204 * gl_53[k]
                  + f_203 * gl_60[k]
                  - f_205 * gl_62[k]
                  + f_206 * gl_64[k]
                  + f_202 * gl_73[k]
                  - f_204 * gl_75[k]
                  + f_206 * gl_77[k]
                  - f_207 * gl_79[k]
                  + f_202 * gl_271[k]
                  + f_203 * gl_276[k]
                  - f_204 * gl_278[k]
                  + f_203 * gl_285[k]
                  - f_205 * gl_287[k]
                  + f_206 * gl_289[k]
                  + f_202 * gl_298[k]
                  - f_204 * gl_300[k]
                  + f_206 * gl_302[k]
                  - f_207 * gl_304[k]
                  - f_208 * gl_361[k]
                  - f_209 * gl_366[k]
                  + f_210 * gl_368[k]
                  - f_209 * gl_375[k]
                  + f_211 * gl_377[k]
                  - f_212 * gl_379[k]
                  - f_208 * gl_388[k]
                  + f_210 * gl_390[k]
                  - f_212 * gl_392[k]
                  + f_213 * gl_394[k];
    }

#pragma omp simd aligned(gl_49, gl_56, gl_58, gl_67, gl_69, gl_71, gl_82, gl_84, gl_86, gl_88, \
                         gl_274, gl_281, gl_283, gl_292, gl_294, gl_296, gl_307, gl_309, \
                         gl_311, gl_313, gl_364, gl_371, gl_373, gl_382, gl_384, gl_386, \
                         gl_397, gl_399, gl_401, gl_403 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_41[k] = f_214 * gl_49[k]
                  + f_215 * gl_56[k]
                  - f_216 * gl_58[k]
                  + f_215 * gl_67[k]
                  - f_217 * gl_69[k]
                  + f_218 * gl_71[k]
                  + f_214 * gl_82[k]
                  - f_216 * gl_84[k]
                  + f_218 * gl_86[k]
                  - f_219 * gl_88[k]
                  + f_214 * gl_274[k]
                  + f_215 * gl_281[k]
                  - f_216 * gl_283[k]
                  + f_215 * gl_292[k]
                  - f_217 * gl_294[k]
                  + f_218 * gl_296[k]
                  + f_214 * gl_307[k]
                  - f_216 * gl_309[k]
                  + f_218 * gl_311[k]
                  - f_219 * gl_313[k]
                  - f_220 * gl_364[k]
                  - f_221 * gl_371[k]
                  + f_222 * gl_373[k]
                  - f_221 * gl_382[k]
                  + f_223 * gl_384[k]
                  - f_224 * gl_386[k]
                  - f_220 * gl_397[k]
                  + f_222 * gl_399[k]
                  - f_224 * gl_401[k]
                  + f_225 * gl_403[k];
    }

#pragma omp simd aligned(gl_45, gl_48, gl_50, gl_55, gl_57, gl_59, gl_66, gl_68, gl_70, gl_72, \
                         gl_81, gl_83, gl_85, gl_87, gl_89, gl_270, gl_273, gl_275, gl_280, \
                         gl_282, gl_284, gl_291, gl_293, gl_295, gl_297, gl_306, gl_308, \
                         gl_310, gl_312, gl_314, gl_360, gl_363, gl_365, gl_370, gl_372, \
                         gl_374, gl_381, gl_383, gl_385, gl_387, gl_396, gl_398, gl_400, \
                         gl_402, gl_404 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = -f_226 * gl_45[k]
                  - f_227 * gl_48[k]
                  + f_228 * gl_50[k]
                  - f_229 * gl_55[k]
                  + f_216 * gl_57[k]
                  - f_216 * gl_59[k]
                  - f_227 * gl_66[k]
                  + f_216 * gl_68[k]
                  - f_217 * gl_70[k]
                  + f_230 * gl_72[k]
                  - f_226 * gl_81[k]
                  + f_228 * gl_83[k]
                  - f_216 * gl_85[k]
                  + f_230 * gl_87[k]
                  - f_231 * gl_89[k]
                  - f_226 * gl_270[k]
                  - f_227 * gl_273[k]
                  + f_228 * gl_275[k]
                  - f_229 * gl_280[k]
                  + f_216 * gl_282[k]
                  - f_216 * gl_284[k]
                  - f_227 * gl_291[k]
                  + f_216 * gl_293[k]
                  - f_217 * gl_295[k]
                  + f_230 * gl_297[k]
                  - f_226 * gl_306[k]
                  + f_228 * gl_308[k]
                  - f_216 * gl_310[k]
                  + f_230 * gl_312[k]
                  - f_231 * gl_314[k]
                  + f_229 * gl_360[k]
                  + f_232 * gl_363[k]
                  - f_217 * gl_365[k]
                  + f_215 * gl_370[k]
                  - f_222 * gl_372[k]
                  + f_222 * gl_374[k]
                  + f_232 * gl_381[k]
                  - f_222 * gl_383[k]
                  + f_223 * gl_385[k]
                  - f_233 * gl_387[k]
                  + f_229 * gl_396[k]
                  - f_217 * gl_398[k]
                  + f_222 * gl_400[k]
                  - f_233 * gl_402[k]
                  + f_219 * gl_404[k];
    }

#pragma omp simd aligned(gl_47, gl_52, gl_54, gl_61, gl_63, gl_65, gl_74, gl_76, gl_78, gl_80, \
                         gl_272, gl_277, gl_279, gl_286, gl_288, gl_290, gl_299, gl_301, \
                         gl_303, gl_305, gl_362, gl_367, gl_369, gl_376, gl_378, gl_380, \
                         gl_389, gl_391, gl_393, gl_395 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_43[k] = f_214 * gl_47[k]
                  + f_215 * gl_52[k]
                  - f_216 * gl_54[k]
                  + f_215 * gl_61[k]
                  - f_217 * gl_63[k]
                  + f_218 * gl_65[k]
                  + f_214 * gl_74[k]
                  - f_216 * gl_76[k]
                  + f_218 * gl_78[k]
                  - f_219 * gl_80[k]
                  + f_214 * gl_272[k]
                  + f_215 * gl_277[k]
                  - f_216 * gl_279[k]
                  + f_215 * gl_286[k]
                  - f_217 * gl_288[k]
                  + f_218 * gl_290[k]
                  + f_214 * gl_299[k]
                  - f_216 * gl_301[k]
                  + f_218 * gl_303[k]
                  - f_219 * gl_305[k]
                  - f_220 * gl_362[k]
                  - f_221 * gl_367[k]
                  + f_222 * gl_369[k]
                  - f_221 * gl_376[k]
                  + f_223 * gl_378[k]
                  - f_224 * gl_380[k]
                  - f_220 * gl_389[k]
                  + f_222 * gl_391[k]
                  - f_224 * gl_393[k]
                  + f_225 * gl_395[k];
    }

#pragma omp simd aligned(gl_45, gl_48, gl_50, gl_57, gl_59, gl_66, gl_68, gl_72, gl_81, gl_83, \
                         gl_85, gl_87, gl_270, gl_273, gl_275, gl_282, gl_284, gl_291, gl_293, \
                         gl_297, gl_306, gl_308, gl_310, gl_312, gl_360, gl_363, gl_365, \
                         gl_372, gl_374, gl_381, gl_383, gl_387, gl_396, gl_398, gl_400, \
                         gl_402 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_44[k] = f_234 * gl_45[k]
                  + f_202 * gl_48[k]
                  - f_235 * gl_50[k]
                  - f_235 * gl_57[k]
                  + f_236 * gl_59[k]
                  - f_202 * gl_66[k]
                  + f_235 * gl_68[k]
                  - f_237 * gl_72[k]
                  - f_234 * gl_81[k]
                  + f_235 * gl_83[k]
                  - f_236 * gl_85[k]
                  + f_237 * gl_87[k]
                  + f_234 * gl_270[k]
                  + f_202 * gl_273[k]
                  - f_235 * gl_275[k]
                  - f_235 * gl_282[k]
                  + f_236 * gl_284[k]
                  - f_202 * gl_291[k]
                  + f_235 * gl_293[k]
                  - f_237 * gl_297[k]
                  - f_234 * gl_306[k]
                  + f_235 * gl_308[k]
                  - f_236 * gl_310[k]
                  + f_237 * gl_312[k]
                  - f_203 * gl_360[k]
                  - f_208 * gl_363[k]
                  + f_238 * gl_365[k]
                  + f_238 * gl_372[k]
                  - f_239 * gl_374[k]
                  + f_208 * gl_381[k]
                  - f_238 * gl_383[k]
                  + f_240 * gl_387[k]
                  + f_203 * gl_396[k]
                  - f_238 * gl_398[k]
                  + f_239 * gl_400[k]
                  - f_240 * gl_402[k];
    }

#pragma omp simd aligned(gl_47, gl_52, gl_54, gl_61, gl_63, gl_65, gl_74, gl_76, gl_78, \
                         gl_272, gl_277, gl_279, gl_286, gl_288, gl_290, gl_299, gl_301, \
                         gl_303, gl_362, gl_367, gl_369, gl_376, gl_378, gl_380, gl_389, \
                         gl_391, gl_393 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_45[k] = -f_189 * gl_47[k]
                  + f_189 * gl_52[k]
                  + f_192 * gl_54[k]
                  + f_187 * gl_61[k]
                  - f_190 * gl_63[k]
                  - f_193 * gl_65[k]
                  + f_186 * gl_74[k]
                  - f_188 * gl_76[k]
                  + f_191 * gl_78[k]
                  - f_189 * gl_272[k]
                  + f_189 * gl_277[k]
                  + f_192 * gl_279[k]
                  + f_187 * gl_286[k]
                  - f_190 * gl_288[k]
                  - f_193 * gl_290[k]
                  + f_186 * gl_299[k]
                  - f_188 * gl_301[k]
                  + f_191 * gl_303[k]
                  + f_197 * gl_362[k]
                  - f_197 * gl_367[k]
                  - f_200 * gl_369[k]
                  - f_195 * gl_376[k]
                  + f_198 * gl_378[k]
                  + f_201 * gl_380[k]
                  - f_194 * gl_389[k]
                  + f_196 * gl_391[k]
                  - f_199 * gl_393[k];
    }

#pragma omp simd aligned(gl_45, gl_48, gl_50, gl_55, gl_57, gl_59, gl_66, gl_68, gl_70, gl_81, \
                         gl_83, gl_85, gl_270, gl_273, gl_275, gl_280, gl_282, gl_284, gl_291, \
                         gl_293, gl_295, gl_306, gl_308, gl_310, gl_360, gl_363, gl_365, \
                         gl_370, gl_372, gl_374, gl_381, gl_383, gl_385, gl_396, gl_398, \
                         gl_400 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_46[k] = -f_241 * gl_45[k]
                  + f_180 * gl_48[k]
                  + f_183 * gl_50[k]
                  + f_242 * gl_55[k]
                  - f_243 * gl_57[k]
                  - f_244 * gl_59[k]
                  + f_180 * gl_66[k]
                  - f_243 * gl_68[k]
                  + f_245 * gl_70[k]
                  - f_241 * gl_81[k]
                  + f_183 * gl_83[k]
                  - f_244 * gl_85[k]
                  - f_241 * gl_270[k]
                  + f_180 * gl_273[k]
                  + f_183 * gl_275[k]
                  + f_242 * gl_280[k]
                  - f_243 * gl_282[k]
                  - f_244 * gl_284[k]
                  + f_180 * gl_291[k]
                  - f_243 * gl_293[k]
                  + f_245 * gl_295[k]
                  - f_241 * gl_306[k]
                  + f_183 * gl_308[k]
                  - f_244 * gl_310[k]
                  + f_246 * gl_360[k]
                  - f_183 * gl_363[k]
                  - f_247 * gl_365[k]
                  - f_248 * gl_370[k]
                  + f_249 * gl_372[k]
                  + f_245 * gl_374[k]
                  - f_183 * gl_381[k]
                  + f_249 * gl_383[k]
                  - f_250 * gl_385[k]
                  + f_246 * gl_396[k]
                  - f_247 * gl_398[k]
                  + f_245 * gl_400[k];
    }

#pragma omp simd aligned(gl_47, gl_52, gl_54, gl_61, gl_63, gl_74, gl_76, gl_272, gl_277, \
                         gl_279, gl_286, gl_288, gl_299, gl_301, gl_362, gl_367, gl_369, \
                         gl_376, gl_378, gl_389, gl_391 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_47[k] = f_172 * gl_47[k]
                  - f_170 * gl_52[k]
                  - f_173 * gl_54[k]
                  - f_168 * gl_61[k]
                  + f_171 * gl_63[k]
                  + f_168 * gl_74[k]
                  - f_169 * gl_76[k]
                  + f_172 * gl_272[k]
                  - f_170 * gl_277[k]
                  - f_173 * gl_279[k]
                  - f_168 * gl_286[k]
                  + f_171 * gl_288[k]
                  + f_168 * gl_299[k]
                  - f_169 * gl_301[k]
                  - f_178 * gl_362[k]
                  + f_176 * gl_367[k]
                  + f_179 * gl_369[k]
                  + f_174 * gl_376[k]
                  - f_177 * gl_378[k]
                  - f_174 * gl_389[k]
                  + f_175 * gl_391[k];
    }

#pragma omp simd aligned(gl_45, gl_48, gl_50, gl_57, gl_66, gl_68, gl_81, gl_83, gl_270, \
                         gl_273, gl_275, gl_282, gl_291, gl_293, gl_306, gl_308, gl_360, \
                         gl_363, gl_365, gl_372, gl_381, gl_383, gl_396, \
                         gl_398 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_48[k] = f_251 * gl_45[k]
                  - f_162 * gl_48[k]
                  - f_162 * gl_50[k]
                  + f_252 * gl_57[k]
                  + f_162 * gl_66[k]
                  - f_252 * gl_68[k]
                  - f_251 * gl_81[k]
                  + f_162 * gl_83[k]
                  + f_251 * gl_270[k]
                  - f_162 * gl_273[k]
                  - f_162 * gl_275[k]
                  + f_252 * gl_282[k]
                  + f_162 * gl_291[k]
                  - f_252 * gl_293[k]
                  - f_251 * gl_306[k]
                  + f_162 * gl_308[k]
                  - f_161 * gl_360[k]
                  + f_163 * gl_363[k]
                  + f_163 * gl_365[k]
                  - f_253 * gl_372[k]
                  - f_163 * gl_381[k]
                  + f_253 * gl_383[k]
                  + f_161 * gl_396[k]
                  - f_163 * gl_398[k];
    }

#pragma omp simd aligned(gl_47, gl_52, gl_61, gl_74, gl_272, gl_277, gl_286, gl_299, gl_362, \
                         gl_367, gl_376, gl_389 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_49[k] = -f_156 * gl_47[k]
                  + f_155 * gl_52[k]
                  - f_154 * gl_61[k]
                  + f_153 * gl_74[k]
                  - f_156 * gl_272[k]
                  + f_155 * gl_277[k]
                  - f_154 * gl_286[k]
                  + f_153 * gl_299[k]
                  + f_160 * gl_362[k]
                  - f_159 * gl_367[k]
                  + f_158 * gl_376[k]
                  - f_157 * gl_389[k];
    }

#pragma omp simd aligned(gl_45, gl_48, gl_55, gl_66, gl_81, gl_270, gl_273, gl_280, gl_291, \
                         gl_306, gl_360, gl_363, gl_370, gl_381, \
                         gl_396 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_50[k] = -f_254 * gl_45[k]
                  + f_153 * gl_48[k]
                  - f_255 * gl_55[k]
                  + f_153 * gl_66[k]
                  - f_254 * gl_81[k]
                  - f_254 * gl_270[k]
                  + f_153 * gl_273[k]
                  - f_255 * gl_280[k]
                  + f_153 * gl_291[k]
                  - f_254 * gl_306[k]
                  + f_256 * gl_360[k]
                  - f_157 * gl_363[k]
                  + f_257 * gl_370[k]
                  - f_157 * gl_381[k]
                  + f_256 * gl_396[k];
    }

#pragma omp simd aligned(gl_181, gl_186, gl_195, gl_208, gl_496, gl_501, gl_510, gl_523, \
                         gl_586, gl_591, gl_600, gl_613 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_51[k] = -f_258 * gl_181[k]
                  + f_259 * gl_186[k]
                  - f_259 * gl_195[k]
                  + f_258 * gl_208[k]
                  - f_258 * gl_496[k]
                  + f_259 * gl_501[k]
                  - f_259 * gl_510[k]
                  + f_258 * gl_523[k]
                  + f_260 * gl_586[k]
                  - f_261 * gl_591[k]
                  + f_261 * gl_600[k]
                  - f_260 * gl_613[k];
    }

#pragma omp simd aligned(gl_184, gl_191, gl_202, gl_217, gl_499, gl_506, gl_517, gl_532, \
                         gl_589, gl_596, gl_607, gl_622 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_52[k] = -f_262 * gl_184[k]
                  + f_263 * gl_191[k]
                  - f_264 * gl_202[k]
                  + f_265 * gl_217[k]
                  - f_262 * gl_499[k]
                  + f_263 * gl_506[k]
                  - f_264 * gl_517[k]
                  + f_265 * gl_532[k]
                  + f_266 * gl_589[k]
                  - f_267 * gl_596[k]
                  + f_268 * gl_607[k]
                  - f_269 * gl_622[k];
    }

#pragma omp simd aligned(gl_181, gl_186, gl_188, gl_195, gl_197, gl_208, gl_210, gl_496, \
                         gl_501, gl_503, gl_510, gl_512, gl_523, gl_525, gl_586, gl_591, \
                         gl_593, gl_600, gl_602, gl_613, gl_615 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_53[k] = f_270 * gl_181[k]
                  - f_271 * gl_186[k]
                  - f_272 * gl_188[k]
                  - f_271 * gl_195[k]
                  + f_273 * gl_197[k]
                  + f_270 * gl_208[k]
                  - f_272 * gl_210[k]
                  + f_270 * gl_496[k]
                  - f_271 * gl_501[k]
                  - f_272 * gl_503[k]
                  - f_271 * gl_510[k]
                  + f_273 * gl_512[k]
                  + f_270 * gl_523[k]
                  - f_272 * gl_525[k]
                  - f_274 * gl_586[k]
                  + f_275 * gl_591[k]
                  + f_276 * gl_593[k]
                  + f_275 * gl_600[k]
                  - f_277 * gl_602[k]
                  - f_274 * gl_613[k]
                  + f_276 * gl_615[k];
    }

#pragma omp simd aligned(gl_184, gl_191, gl_193, gl_202, gl_204, gl_217, gl_219, gl_499, \
                         gl_506, gl_508, gl_517, gl_519, gl_532, gl_534, gl_589, gl_596, \
                         gl_598, gl_607, gl_609, gl_622, gl_624 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_54[k] = f_278 * gl_184[k]
                  - f_278 * gl_191[k]
                  - f_279 * gl_193[k]
                  - f_280 * gl_202[k]
                  + f_281 * gl_204[k]
                  + f_282 * gl_217[k]
                  - f_283 * gl_219[k]
                  + f_278 * gl_499[k]
                  - f_278 * gl_506[k]
                  - f_279 * gl_508[k]
                  - f_280 * gl_517[k]
                  + f_281 * gl_519[k]
                  + f_282 * gl_532[k]
                  - f_283 * gl_534[k]
                  - f_284 * gl_589[k]
                  + f_284 * gl_596[k]
                  + f_285 * gl_598[k]
                  + f_286 * gl_607[k]
                  - f_287 * gl_609[k]
                  - f_288 * gl_622[k]
                  + f_289 * gl_624[k];
    }

#pragma omp simd aligned(gl_181, gl_186, gl_188, gl_195, gl_199, gl_208, gl_210, gl_212, \
                         gl_496, gl_501, gl_503, gl_510, gl_514, gl_523, gl_525, gl_527, \
                         gl_586, gl_591, gl_593, gl_600, gl_604, gl_613, gl_615, \
                         gl_617 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_55[k] = -f_290 * gl_181[k]
                  - f_290 * gl_186[k]
                  + f_291 * gl_188[k]
                  + f_290 * gl_195[k]
                  - f_292 * gl_199[k]
                  + f_290 * gl_208[k]
                  - f_291 * gl_210[k]
                  + f_292 * gl_212[k]
                  - f_290 * gl_496[k]
                  - f_290 * gl_501[k]
                  + f_291 * gl_503[k]
                  + f_290 * gl_510[k]
                  - f_292 * gl_514[k]
                  + f_290 * gl_523[k]
                  - f_291 * gl_525[k]
                  + f_292 * gl_527[k]
                  + f_293 * gl_586[k]
                  + f_293 * gl_591[k]
                  - f_294 * gl_593[k]
                  - f_293 * gl_600[k]
                  + f_295 * gl_604[k]
                  - f_293 * gl_613[k]
                  + f_294 * gl_615[k]
                  - f_295 * gl_617[k];
    }

#pragma omp simd aligned(gl_184, gl_191, gl_193, gl_202, gl_204, gl_206, gl_217, gl_219, \
                         gl_221, gl_499, gl_506, gl_508, gl_517, gl_519, gl_521, gl_532, \
                         gl_534, gl_536, gl_589, gl_596, gl_598, gl_607, gl_609, gl_611, \
                         gl_622, gl_624, gl_626 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_56[k] = -f_296 * gl_184[k]
                  - f_297 * gl_191[k]
                  + f_298 * gl_193[k]
                  - f_299 * gl_202[k]
                  + f_300 * gl_204[k]
                  - f_301 * gl_206[k]
                  + f_299 * gl_217[k]
                  - f_302 * gl_219[k]
                  + f_303 * gl_221[k]
                  - f_296 * gl_499[k]
                  - f_297 * gl_506[k]
                  + f_298 * gl_508[k]
                  - f_299 * gl_517[k]
                  + f_300 * gl_519[k]
                  - f_301 * gl_521[k]
                  + f_299 * gl_532[k]
                  - f_302 * gl_534[k]
                  + f_303 * gl_536[k]
                  + f_304 * gl_589[k]
                  + f_302 * gl_596[k]
                  - f_305 * gl_598[k]
                  + f_306 * gl_607[k]
                  - f_307 * gl_609[k]
                  + f_308 * gl_611[k]
                  - f_306 * gl_622[k]
                  + f_309 * gl_624[k]
                  - f_310 * gl_626[k];
    }

#pragma omp simd aligned(gl_181, gl_186, gl_188, gl_195, gl_197, gl_199, gl_208, gl_210, \
                         gl_212, gl_214, gl_496, gl_501, gl_503, gl_510, gl_512, gl_514, \
                         gl_523, gl_525, gl_527, gl_529, gl_586, gl_591, gl_593, gl_600, \
                         gl_602, gl_604, gl_613, gl_615, gl_617, \
                         gl_619 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_57[k] = f_311 * gl_181[k]
                  + f_312 * gl_186[k]
                  - f_313 * gl_188[k]
                  + f_312 * gl_195[k]
                  - f_314 * gl_197[k]
                  + f_315 * gl_199[k]
                  + f_311 * gl_208[k]
                  - f_313 * gl_210[k]
                  + f_315 * gl_212[k]
                  - f_316 * gl_214[k]
                  + f_311 * gl_496[k]
                  + f_312 * gl_501[k]
                  - f_313 * gl_503[k]
                  + f_312 * gl_510[k]
                  - f_314 * gl_512[k]
                  + f_315 * gl_514[k]
                  + f_311 * gl_523[k]
                  - f_313 * gl_525[k]
                  + f_315 * gl_527[k]
                  - f_316 * gl_529[k]
                  - f_317 * gl_586[k]
                  - f_318 * gl_591[k]
                  + f_319 * gl_593[k]
                  - f_318 * gl_600[k]
                  + f_315 * gl_602[k]
                  - f_320 * gl_604[k]
                  - f_317 * gl_613[k]
                  + f_319 * gl_615[k]
                  - f_320 * gl_617[k]
                  + f_321 * gl_619[k];
    }

#pragma omp simd aligned(gl_184, gl_191, gl_193, gl_202, gl_204, gl_206, gl_217, gl_219, \
                         gl_221, gl_223, gl_499, gl_506, gl_508, gl_517, gl_519, gl_521, \
                         gl_532, gl_534, gl_536, gl_538, gl_589, gl_596, gl_598, gl_607, \
                         gl_609, gl_611, gl_622, gl_624, gl_626, \
                         gl_628 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_58[k] = f_322 * gl_184[k]
                  + f_323 * gl_191[k]
                  - f_324 * gl_193[k]
                  + f_323 * gl_202[k]
                  - f_325 * gl_204[k]
                  + f_326 * gl_206[k]
                  + f_322 * gl_217[k]
                  - f_324 * gl_219[k]
                  + f_326 * gl_221[k]
                  - f_327 * gl_223[k]
                  + f_322 * gl_499[k]
                  + f_323 * gl_506[k]
                  - f_324 * gl_508[k]
                  + f_323 * gl_517[k]
                  - f_325 * gl_519[k]
                  + f_326 * gl_521[k]
                  + f_322 * gl_532[k]
                  - f_324 * gl_534[k]
                  + f_326 * gl_536[k]
                  - f_327 * gl_538[k]
                  - f_328 * gl_589[k]
                  - f_329 * gl_596[k]
                  + f_330 * gl_598[k]
                  - f_329 * gl_607[k]
                  + f_331 * gl_609[k]
                  - f_332 * gl_611[k]
                  - f_328 * gl_622[k]
                  + f_330 * gl_624[k]
                  - f_332 * gl_626[k]
                  + f_333 * gl_628[k];
    }

#pragma omp simd aligned(gl_180, gl_183, gl_185, gl_190, gl_192, gl_194, gl_201, gl_203, \
                         gl_205, gl_207, gl_216, gl_218, gl_220, gl_222, gl_224, gl_495, \
                         gl_498, gl_500, gl_505, gl_507, gl_509, gl_516, gl_518, gl_520, \
                         gl_522, gl_531, gl_533, gl_535, gl_537, gl_539, gl_585, gl_588, \
                         gl_590, gl_595, gl_597, gl_599, gl_606, gl_608, gl_610, gl_612, \
                         gl_621, gl_623, gl_625, gl_627, gl_629 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_59[k] = -f_334 * gl_180[k]
                  - f_335 * gl_183[k]
                  + f_336 * gl_185[k]
                  - f_337 * gl_190[k]
                  + f_324 * gl_192[k]
                  - f_324 * gl_194[k]
                  - f_335 * gl_201[k]
                  + f_324 * gl_203[k]
                  - f_325 * gl_205[k]
                  + f_338 * gl_207[k]
                  - f_334 * gl_216[k]
                  + f_336 * gl_218[k]
                  - f_324 * gl_220[k]
                  + f_338 * gl_222[k]
                  - f_339 * gl_224[k]
                  - f_334 * gl_495[k]
                  - f_335 * gl_498[k]
                  + f_336 * gl_500[k]
                  - f_337 * gl_505[k]
                  + f_324 * gl_507[k]
                  - f_324 * gl_509[k]
                  - f_335 * gl_516[k]
                  + f_324 * gl_518[k]
                  - f_325 * gl_520[k]
                  + f_338 * gl_522[k]
                  - f_334 * gl_531[k]
                  + f_336 * gl_533[k]
                  - f_324 * gl_535[k]
                  + f_338 * gl_537[k]
                  - f_339 * gl_539[k]
                  + f_340 * gl_585[k]
                  + f_341 * gl_588[k]
                  - f_342 * gl_590[k]
                  + f_343 * gl_595[k]
                  - f_330 * gl_597[k]
                  + f_330 * gl_599[k]
                  + f_341 * gl_606[k]
                  - f_330 * gl_608[k]
                  + f_331 * gl_610[k]
                  - f_344 * gl_612[k]
                  + f_340 * gl_621[k]
                  - f_342 * gl_623[k]
                  + f_330 * gl_625[k]
                  - f_344 * gl_627[k]
                  + f_345 * gl_629[k];
    }

#pragma omp simd aligned(gl_182, gl_187, gl_189, gl_196, gl_198, gl_200, gl_209, gl_211, \
                         gl_213, gl_215, gl_497, gl_502, gl_504, gl_511, gl_513, gl_515, \
                         gl_524, gl_526, gl_528, gl_530, gl_587, gl_592, gl_594, gl_601, \
                         gl_603, gl_605, gl_614, gl_616, gl_618, \
                         gl_620 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_60[k] = f_322 * gl_182[k]
                  + f_323 * gl_187[k]
                  - f_324 * gl_189[k]
                  + f_323 * gl_196[k]
                  - f_325 * gl_198[k]
                  + f_326 * gl_200[k]
                  + f_322 * gl_209[k]
                  - f_324 * gl_211[k]
                  + f_326 * gl_213[k]
                  - f_327 * gl_215[k]
                  + f_322 * gl_497[k]
                  + f_323 * gl_502[k]
                  - f_324 * gl_504[k]
                  + f_323 * gl_511[k]
                  - f_325 * gl_513[k]
                  + f_326 * gl_515[k]
                  + f_322 * gl_524[k]
                  - f_324 * gl_526[k]
                  + f_326 * gl_528[k]
                  - f_327 * gl_530[k]
                  - f_328 * gl_587[k]
                  - f_329 * gl_592[k]
                  + f_330 * gl_594[k]
                  - f_329 * gl_601[k]
                  + f_331 * gl_603[k]
                  - f_332 * gl_605[k]
                  - f_328 * gl_614[k]
                  + f_330 * gl_616[k]
                  - f_332 * gl_618[k]
                  + f_333 * gl_620[k];
    }

#pragma omp simd aligned(gl_180, gl_183, gl_185, gl_192, gl_194, gl_201, gl_203, gl_207, \
                         gl_216, gl_218, gl_220, gl_222, gl_495, gl_498, gl_500, gl_507, \
                         gl_509, gl_516, gl_518, gl_522, gl_531, gl_533, gl_535, gl_537, \
                         gl_585, gl_588, gl_590, gl_597, gl_599, gl_606, gl_608, gl_612, \
                         gl_621, gl_623, gl_625, gl_627 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_61[k] = f_346 * gl_180[k]
                  + f_311 * gl_183[k]
                  - f_347 * gl_185[k]
                  - f_347 * gl_192[k]
                  + f_319 * gl_194[k]
                  - f_311 * gl_201[k]
                  + f_347 * gl_203[k]
                  - f_348 * gl_207[k]
                  - f_346 * gl_216[k]
                  + f_347 * gl_218[k]
                  - f_319 * gl_220[k]
                  + f_348 * gl_222[k]
                  + f_346 * gl_495[k]
                  + f_311 * gl_498[k]
                  - f_347 * gl_500[k]
                  - f_347 * gl_507[k]
                  + f_319 * gl_509[k]
                  - f_311 * gl_516[k]
                  + f_347 * gl_518[k]
                  - f_348 * gl_522[k]
                  - f_346 * gl_531[k]
                  + f_347 * gl_533[k]
                  - f_319 * gl_535[k]
                  + f_348 * gl_537[k]
                  - f_349 * gl_585[k]
                  - f_317 * gl_588[k]
                  + f_350 * gl_590[k]
                  + f_350 * gl_597[k]
                  - f_351 * gl_599[k]
                  + f_317 * gl_606[k]
                  - f_350 * gl_608[k]
                  + f_352 * gl_612[k]
                  + f_349 * gl_621[k]
                  - f_350 * gl_623[k]
                  + f_351 * gl_625[k]
                  - f_352 * gl_627[k];
    }

#pragma omp simd aligned(gl_182, gl_187, gl_189, gl_196, gl_198, gl_200, gl_209, gl_211, \
                         gl_213, gl_497, gl_502, gl_504, gl_511, gl_513, gl_515, gl_524, \
                         gl_526, gl_528, gl_587, gl_592, gl_594, gl_601, gl_603, gl_605, \
                         gl_614, gl_616, gl_618 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_62[k] = -f_299 * gl_182[k]
                  + f_299 * gl_187[k]
                  + f_302 * gl_189[k]
                  + f_297 * gl_196[k]
                  - f_300 * gl_198[k]
                  - f_303 * gl_200[k]
                  + f_296 * gl_209[k]
                  - f_298 * gl_211[k]
                  + f_301 * gl_213[k]
                  - f_299 * gl_497[k]
                  + f_299 * gl_502[k]
                  + f_302 * gl_504[k]
                  + f_297 * gl_511[k]
                  - f_300 * gl_513[k]
                  - f_303 * gl_515[k]
                  + f_296 * gl_524[k]
                  - f_298 * gl_526[k]
                  + f_301 * gl_528[k]
                  + f_306 * gl_587[k]
                  - f_306 * gl_592[k]
                  - f_309 * gl_594[k]
                  - f_302 * gl_601[k]
                  + f_307 * gl_603[k]
                  + f_310 * gl_605[k]
                  - f_304 * gl_614[k]
                  + f_305 * gl_616[k]
                  - f_308 * gl_618[k];
    }

#pragma omp simd aligned(gl_180, gl_183, gl_185, gl_190, gl_192, gl_194, gl_201, gl_203, \
                         gl_205, gl_216, gl_218, gl_220, gl_495, gl_498, gl_500, gl_505, \
                         gl_507, gl_509, gl_516, gl_518, gl_520, gl_531, gl_533, gl_535, \
                         gl_585, gl_588, gl_590, gl_595, gl_597, gl_599, gl_606, gl_608, \
                         gl_610, gl_621, gl_623, gl_625 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_63[k] = -f_353 * gl_180[k]
                  + f_290 * gl_183[k]
                  + f_354 * gl_185[k]
                  + f_355 * gl_190[k]
                  - f_356 * gl_192[k]
                  - f_357 * gl_194[k]
                  + f_290 * gl_201[k]
                  - f_356 * gl_203[k]
                  + f_358 * gl_205[k]
                  - f_353 * gl_216[k]
                  + f_354 * gl_218[k]
                  - f_357 * gl_220[k]
                  - f_353 * gl_495[k]
                  + f_290 * gl_498[k]
                  + f_354 * gl_500[k]
                  + f_355 * gl_505[k]
                  - f_356 * gl_507[k]
                  - f_357 * gl_509[k]
                  + f_290 * gl_516[k]
                  - f_356 * gl_518[k]
                  + f_358 * gl_520[k]
                  - f_353 * gl_531[k]
                  + f_354 * gl_533[k]
                  - f_357 * gl_535[k]
                  + f_359 * gl_585[k]
                  - f_293 * gl_588[k]
                  - f_360 * gl_590[k]
                  - f_361 * gl_595[k]
                  + f_292 * gl_597[k]
                  + f_362 * gl_599[k]
                  - f_293 * gl_606[k]
                  + f_292 * gl_608[k]
                  - f_363 * gl_610[k]
                  + f_359 * gl_621[k]
                  - f_360 * gl_623[k]
                  + f_362 * gl_625[k];
    }

#pragma omp simd aligned(gl_182, gl_187, gl_189, gl_196, gl_198, gl_209, gl_211, gl_497, \
                         gl_502, gl_504, gl_511, gl_513, gl_524, gl_526, gl_587, gl_592, \
                         gl_594, gl_601, gl_603, gl_614, gl_616 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_64[k] = f_282 * gl_182[k]
                  - f_280 * gl_187[k]
                  - f_283 * gl_189[k]
                  - f_278 * gl_196[k]
                  + f_281 * gl_198[k]
                  + f_278 * gl_209[k]
                  - f_279 * gl_211[k]
                  + f_282 * gl_497[k]
                  - f_280 * gl_502[k]
                  - f_283 * gl_504[k]
                  - f_278 * gl_511[k]
                  + f_281 * gl_513[k]
                  + f_278 * gl_524[k]
                  - f_279 * gl_526[k]
                  - f_288 * gl_587[k]
                  + f_286 * gl_592[k]
                  + f_289 * gl_594[k]
                  + f_284 * gl_601[k]
                  - f_287 * gl_603[k]
                  - f_284 * gl_614[k]
                  + f_285 * gl_616[k];
    }

#pragma omp simd aligned(gl_180, gl_183, gl_185, gl_192, gl_201, gl_203, gl_216, gl_218, \
                         gl_495, gl_498, gl_500, gl_507, gl_516, gl_518, gl_531, gl_533, \
                         gl_585, gl_588, gl_590, gl_597, gl_606, gl_608, gl_621, \
                         gl_623 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_65[k] = f_364 * gl_180[k]
                  - f_271 * gl_183[k]
                  - f_271 * gl_185[k]
                  + f_365 * gl_192[k]
                  + f_271 * gl_201[k]
                  - f_365 * gl_203[k]
                  - f_364 * gl_216[k]
                  + f_271 * gl_218[k]
                  + f_364 * gl_495[k]
                  - f_271 * gl_498[k]
                  - f_271 * gl_500[k]
                  + f_365 * gl_507[k]
                  + f_271 * gl_516[k]
                  - f_365 * gl_518[k]
                  - f_364 * gl_531[k]
                  + f_271 * gl_533[k]
                  - f_366 * gl_585[k]
                  + f_275 * gl_588[k]
                  + f_275 * gl_590[k]
                  - f_273 * gl_597[k]
                  - f_275 * gl_606[k]
                  + f_273 * gl_608[k]
                  + f_366 * gl_621[k]
                  - f_275 * gl_623[k];
    }

#pragma omp simd aligned(gl_182, gl_187, gl_196, gl_209, gl_497, gl_502, gl_511, gl_524, \
                         gl_587, gl_592, gl_601, gl_614 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_66[k] = -f_265 * gl_182[k]
                  + f_264 * gl_187[k]
                  - f_263 * gl_196[k]
                  + f_262 * gl_209[k]
                  - f_265 * gl_497[k]
                  + f_264 * gl_502[k]
                  - f_263 * gl_511[k]
                  + f_262 * gl_524[k]
                  + f_269 * gl_587[k]
                  - f_268 * gl_592[k]
                  + f_267 * gl_601[k]
                  - f_266 * gl_614[k];
    }

#pragma omp simd aligned(gl_180, gl_183, gl_190, gl_201, gl_216, gl_495, gl_498, gl_505, \
                         gl_516, gl_531, gl_585, gl_588, gl_595, gl_606, \
                         gl_621 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_67[k] = -f_367 * gl_180[k]
                  + f_262 * gl_183[k]
                  - f_368 * gl_190[k]
                  + f_262 * gl_201[k]
                  - f_367 * gl_216[k]
                  - f_367 * gl_495[k]
                  + f_262 * gl_498[k]
                  - f_368 * gl_505[k]
                  + f_262 * gl_516[k]
                  - f_367 * gl_531[k]
                  + f_369 * gl_585[k]
                  - f_266 * gl_588[k]
                  + f_370 * gl_595[k]
                  - f_266 * gl_606[k]
                  + f_369 * gl_621[k];
    }

#pragma omp simd aligned(gl_1, gl_6, gl_15, gl_28, gl_136, gl_141, gl_150, gl_163, gl_226, \
                         gl_231, gl_240, gl_253, gl_451, gl_456, gl_465, gl_478, gl_541, \
                         gl_546, gl_555, gl_568, gl_631, gl_636, gl_645, \
                         gl_658 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_68[k] = f_371 * gl_1[k]
                  - f_372 * gl_6[k]
                  + f_372 * gl_15[k]
                  - f_371 * gl_28[k]
                  + f_373 * gl_136[k]
                  - f_374 * gl_141[k]
                  + f_374 * gl_150[k]
                  - f_373 * gl_163[k]
                  - f_375 * gl_226[k]
                  + f_376 * gl_231[k]
                  - f_376 * gl_240[k]
                  + f_375 * gl_253[k]
                  + f_371 * gl_451[k]
                  - f_372 * gl_456[k]
                  + f_372 * gl_465[k]
                  - f_371 * gl_478[k]
                  - f_375 * gl_541[k]
                  + f_376 * gl_546[k]
                  - f_376 * gl_555[k]
                  + f_375 * gl_568[k]
                  + f_377 * gl_631[k]
                  - f_15 * gl_636[k]
                  + f_15 * gl_645[k]
                  - f_377 * gl_658[k];
    }

#pragma omp simd aligned(gl_4, gl_11, gl_22, gl_37, gl_139, gl_146, gl_157, gl_172, gl_229, \
                         gl_236, gl_247, gl_262, gl_454, gl_461, gl_472, gl_487, gl_544, \
                         gl_551, gl_562, gl_577, gl_634, gl_641, gl_652, \
                         gl_667 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_69[k] = f_378 * gl_4[k]
                  - f_379 * gl_11[k]
                  + f_380 * gl_22[k]
                  - f_381 * gl_37[k]
                  + f_372 * gl_139[k]
                  - f_382 * gl_146[k]
                  + f_383 * gl_157[k]
                  - f_371 * gl_172[k]
                  - f_384 * gl_229[k]
                  + f_385 * gl_236[k]
                  - f_386 * gl_247[k]
                  + f_387 * gl_262[k]
                  + f_378 * gl_454[k]
                  - f_379 * gl_461[k]
                  + f_380 * gl_472[k]
                  - f_381 * gl_487[k]
                  - f_384 * gl_544[k]
                  + f_385 * gl_551[k]
                  - f_386 * gl_562[k]
                  + f_387 * gl_577[k]
                  + f_388 * gl_634[k]
                  - f_389 * gl_641[k]
                  + f_384 * gl_652[k]
                  - f_390 * gl_667[k];
    }

#pragma omp simd aligned(gl_1, gl_6, gl_8, gl_15, gl_17, gl_28, gl_30, gl_136, gl_141, gl_143, \
                         gl_150, gl_152, gl_163, gl_165, gl_226, gl_231, gl_233, gl_240, \
                         gl_242, gl_253, gl_255, gl_451, gl_456, gl_458, gl_465, gl_467, \
                         gl_478, gl_480, gl_541, gl_546, gl_548, gl_555, gl_557, gl_568, \
                         gl_570, gl_631, gl_636, gl_638, gl_645, gl_647, gl_658, \
                         gl_660 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_70[k] = -f_391 * gl_1[k]
                  + f_392 * gl_6[k]
                  + f_393 * gl_8[k]
                  + f_392 * gl_15[k]
                  - f_394 * gl_17[k]
                  - f_391 * gl_28[k]
                  + f_393 * gl_30[k]
                  - f_395 * gl_136[k]
                  + f_396 * gl_141[k]
                  + f_397 * gl_143[k]
                  + f_396 * gl_150[k]
                  - f_398 * gl_152[k]
                  - f_395 * gl_163[k]
                  + f_397 * gl_165[k]
                  + f_399 * gl_226[k]
                  - f_400 * gl_231[k]
                  - f_401 * gl_233[k]
                  - f_400 * gl_240[k]
                  + f_402 * gl_242[k]
                  + f_399 * gl_253[k]
                  - f_401 * gl_255[k]
                  - f_391 * gl_451[k]
                  + f_392 * gl_456[k]
                  + f_393 * gl_458[k]
                  + f_392 * gl_465[k]
                  - f_394 * gl_467[k]
                  - f_391 * gl_478[k]
                  + f_393 * gl_480[k]
                  + f_399 * gl_541[k]
                  - f_400 * gl_546[k]
                  - f_401 * gl_548[k]
                  - f_400 * gl_555[k]
                  + f_402 * gl_557[k]
                  + f_399 * gl_568[k]
                  - f_401 * gl_570[k]
                  - f_403 * gl_631[k]
                  + f_404 * gl_636[k]
                  + f_405 * gl_638[k]
                  + f_404 * gl_645[k]
                  - f_406 * gl_647[k]
                  - f_403 * gl_658[k]
                  + f_405 * gl_660[k];
    }

#pragma omp simd aligned(gl_4, gl_11, gl_13, gl_22, gl_24, gl_37, gl_39, gl_139, gl_146, \
                         gl_148, gl_157, gl_159, gl_172, gl_174, gl_229, gl_236, gl_238, \
                         gl_247, gl_249, gl_262, gl_264, gl_454, gl_461, gl_463, gl_472, \
                         gl_474, gl_487, gl_489, gl_544, gl_551, gl_553, gl_562, gl_564, \
                         gl_577, gl_579, gl_634, gl_641, gl_643, gl_652, gl_654, gl_667, \
                         gl_669 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_71[k] = -f_407 * gl_4[k]
                  + f_407 * gl_11[k]
                  + f_408 * gl_13[k]
                  + f_409 * gl_22[k]
                  - f_410 * gl_24[k]
                  - f_411 * gl_37[k]
                  + f_412 * gl_39[k]
                  - f_413 * gl_139[k]
                  + f_413 * gl_146[k]
                  + f_410 * gl_148[k]
                  + f_414 * gl_157[k]
                  - f_415 * gl_159[k]
                  - f_416 * gl_172[k]
                  + f_417 * gl_174[k]
                  + f_410 * gl_229[k]
                  - f_410 * gl_236[k]
                  - f_418 * gl_238[k]
                  - f_419 * gl_247[k]
                  + f_420 * gl_249[k]
                  + f_417 * gl_262[k]
                  - f_421 * gl_264[k]
                  - f_407 * gl_454[k]
                  + f_407 * gl_461[k]
                  + f_408 * gl_463[k]
                  + f_409 * gl_472[k]
                  - f_410 * gl_474[k]
                  - f_411 * gl_487[k]
                  + f_412 * gl_489[k]
                  + f_410 * gl_544[k]
                  - f_410 * gl_551[k]
                  - f_418 * gl_553[k]
                  - f_419 * gl_562[k]
                  + f_420 * gl_564[k]
                  + f_417 * gl_577[k]
                  - f_421 * gl_579[k]
                  - f_0 * gl_634[k]
                  + f_0 * gl_641[k]
                  + f_422 * gl_643[k]
                  + f_423 * gl_652[k]
                  - f_424 * gl_654[k]
                  - f_425 * gl_667[k]
                  + f_426 * gl_669[k];
    }

#pragma omp simd aligned(gl_1, gl_6, gl_8, gl_15, gl_19, gl_28, gl_30, gl_32, gl_136, gl_141, \
                         gl_143, gl_150, gl_154, gl_163, gl_165, gl_167, gl_226, gl_231, \
                         gl_233, gl_240, gl_244, gl_253, gl_255, gl_257, gl_451, gl_456, \
                         gl_458, gl_465, gl_469, gl_478, gl_480, gl_482, gl_541, gl_546, \
                         gl_548, gl_555, gl_559, gl_568, gl_570, gl_572, gl_631, gl_636, \
                         gl_638, gl_645, gl_649, gl_658, gl_660, \
                         gl_662 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_72[k] = f_427 * gl_1[k]
                  + f_427 * gl_6[k]
                  - f_428 * gl_8[k]
                  - f_427 * gl_15[k]
                  + f_429 * gl_19[k]
                  - f_427 * gl_28[k]
                  + f_428 * gl_30[k]
                  - f_429 * gl_32[k]
                  + f_430 * gl_136[k]
                  + f_430 * gl_141[k]
                  - f_431 * gl_143[k]
                  - f_430 * gl_150[k]
                  + f_432 * gl_154[k]
                  - f_430 * gl_163[k]
                  + f_431 * gl_165[k]
                  - f_432 * gl_167[k]
                  - f_433 * gl_226[k]
                  - f_433 * gl_231[k]
                  + f_434 * gl_233[k]
                  + f_433 * gl_240[k]
                  - f_435 * gl_244[k]
                  + f_433 * gl_253[k]
                  - f_434 * gl_255[k]
                  + f_435 * gl_257[k]
                  + f_427 * gl_451[k]
                  + f_427 * gl_456[k]
                  - f_428 * gl_458[k]
                  - f_427 * gl_465[k]
                  + f_429 * gl_469[k]
                  - f_427 * gl_478[k]
                  + f_428 * gl_480[k]
                  - f_429 * gl_482[k]
                  - f_433 * gl_541[k]
                  - f_433 * gl_546[k]
                  + f_434 * gl_548[k]
                  + f_433 * gl_555[k]
                  - f_435 * gl_559[k]
                  + f_433 * gl_568[k]
                  - f_434 * gl_570[k]
                  + f_435 * gl_572[k]
                  + f_436 * gl_631[k]
                  + f_436 * gl_636[k]
                  - f_437 * gl_638[k]
                  - f_436 * gl_645[k]
                  + f_438 * gl_649[k]
                  - f_436 * gl_658[k]
                  + f_437 * gl_660[k]
                  - f_438 * gl_662[k];
    }

#pragma omp simd aligned(gl_4, gl_11, gl_13, gl_22, gl_24, gl_26, gl_37, gl_39, gl_41, gl_139, \
                         gl_146, gl_148, gl_157, gl_159, gl_161, gl_172, gl_174, gl_176, \
                         gl_229, gl_236, gl_238, gl_247, gl_249, gl_251, gl_262, gl_264, \
                         gl_266, gl_454, gl_461, gl_463, gl_472, gl_474, gl_476, gl_487, \
                         gl_489, gl_491, gl_544, gl_551, gl_553, gl_562, gl_564, gl_566, \
                         gl_577, gl_579, gl_581, gl_634, gl_641, gl_643, gl_652, gl_654, \
                         gl_656, gl_667, gl_669, gl_671 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_73[k] = f_439 * gl_4[k]
                  + f_440 * gl_11[k]
                  - f_441 * gl_13[k]
                  + f_442 * gl_22[k]
                  - f_443 * gl_24[k]
                  + f_444 * gl_26[k]
                  - f_442 * gl_37[k]
                  + f_445 * gl_39[k]
                  - f_446 * gl_41[k]
                  + f_447 * gl_139[k]
                  + f_448 * gl_146[k]
                  - f_449 * gl_148[k]
                  + f_450 * gl_157[k]
                  - f_451 * gl_159[k]
                  + f_452 * gl_161[k]
                  - f_450 * gl_172[k]
                  + f_443 * gl_174[k]
                  - f_453 * gl_176[k]
                  - f_454 * gl_229[k]
                  - f_449 * gl_236[k]
                  + f_455 * gl_238[k]
                  - f_456 * gl_247[k]
                  + f_457 * gl_249[k]
                  - f_458 * gl_251[k]
                  + f_456 * gl_262[k]
                  - f_459 * gl_264[k]
                  + f_460 * gl_266[k]
                  + f_439 * gl_454[k]
                  + f_440 * gl_461[k]
                  - f_441 * gl_463[k]
                  + f_442 * gl_472[k]
                  - f_443 * gl_474[k]
                  + f_444 * gl_476[k]
                  - f_442 * gl_487[k]
                  + f_445 * gl_489[k]
                  - f_446 * gl_491[k]
                  - f_454 * gl_544[k]
                  - f_449 * gl_551[k]
                  + f_455 * gl_553[k]
                  - f_456 * gl_562[k]
                  + f_457 * gl_564[k]
                  - f_458 * gl_566[k]
                  + f_456 * gl_577[k]
                  - f_459 * gl_579[k]
                  + f_460 * gl_581[k]
                  + f_456 * gl_634[k]
                  + f_443 * gl_641[k]
                  - f_459 * gl_643[k]
                  + f_461 * gl_652[k]
                  - f_462 * gl_654[k]
                  + f_460 * gl_656[k]
                  - f_461 * gl_667[k]
                  + f_463 * gl_669[k]
                  - f_464 * gl_671[k];
    }

#pragma omp simd aligned(gl_1, gl_6, gl_8, gl_15, gl_17, gl_19, gl_28, gl_30, gl_32, gl_34, \
                         gl_136, gl_141, gl_143, gl_150, gl_152, gl_154, gl_163, gl_165, \
                         gl_167, gl_169, gl_226, gl_231, gl_233, gl_240, gl_242, gl_244, \
                         gl_253, gl_255, gl_257, gl_259, gl_451, gl_456, gl_458, gl_465, \
                         gl_467, gl_469, gl_478, gl_480, gl_482, gl_484, gl_541, gl_546, \
                         gl_548, gl_555, gl_557, gl_559, gl_568, gl_570, gl_572, gl_574, \
                         gl_631, gl_636, gl_638, gl_645, gl_647, gl_649, gl_658, gl_660, \
                         gl_662, gl_664 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_74[k] = -f_465 * gl_1[k]
                  - f_466 * gl_6[k]
                  + f_467 * gl_8[k]
                  - f_466 * gl_15[k]
                  + f_468 * gl_17[k]
                  - f_469 * gl_19[k]
                  - f_465 * gl_28[k]
                  + f_467 * gl_30[k]
                  - f_469 * gl_32[k]
                  + f_470 * gl_34[k]
                  - f_471 * gl_136[k]
                  - f_472 * gl_141[k]
                  + f_468 * gl_143[k]
                  - f_472 * gl_150[k]
                  + f_473 * gl_152[k]
                  - f_474 * gl_154[k]
                  - f_471 * gl_163[k]
                  + f_468 * gl_165[k]
                  - f_474 * gl_167[k]
                  + f_475 * gl_169[k]
                  + f_476 * gl_226[k]
                  + f_477 * gl_231[k]
                  - f_478 * gl_233[k]
                  + f_477 * gl_240[k]
                  - f_479 * gl_242[k]
                  + f_480 * gl_244[k]
                  + f_476 * gl_253[k]
                  - f_478 * gl_255[k]
                  + f_480 * gl_257[k]
                  - f_481 * gl_259[k]
                  - f_465 * gl_451[k]
                  - f_466 * gl_456[k]
                  + f_467 * gl_458[k]
                  - f_466 * gl_465[k]
                  + f_468 * gl_467[k]
                  - f_469 * gl_469[k]
                  - f_465 * gl_478[k]
                  + f_467 * gl_480[k]
                  - f_469 * gl_482[k]
                  + f_470 * gl_484[k]
                  + f_476 * gl_541[k]
                  + f_477 * gl_546[k]
                  - f_478 * gl_548[k]
                  + f_477 * gl_555[k]
                  - f_479 * gl_557[k]
                  + f_480 * gl_559[k]
                  + f_476 * gl_568[k]
                  - f_478 * gl_570[k]
                  + f_480 * gl_572[k]
                  - f_481 * gl_574[k]
                  - f_482 * gl_631[k]
                  - f_476 * gl_636[k]
                  + f_469 * gl_638[k]
                  - f_476 * gl_645[k]
                  + f_474 * gl_647[k]
                  - f_483 * gl_649[k]
                  - f_482 * gl_658[k]
                  + f_469 * gl_660[k]
                  - f_483 * gl_662[k]
                  + f_484 * gl_664[k];
    }

#pragma omp simd aligned(gl_4, gl_11, gl_13, gl_22, gl_24, gl_26, gl_37, gl_39, gl_41, gl_43, \
                         gl_139, gl_146, gl_148, gl_157, gl_159, gl_161, gl_172, gl_174, \
                         gl_176, gl_178, gl_229, gl_236, gl_238, gl_247, gl_249, gl_251, \
                         gl_262, gl_264, gl_266, gl_268, gl_454, gl_461, gl_463, gl_472, \
                         gl_474, gl_476, gl_487, gl_489, gl_491, gl_493, gl_544, gl_551, \
                         gl_553, gl_562, gl_564, gl_566, gl_577, gl_579, gl_581, gl_583, \
                         gl_634, gl_641, gl_643, gl_652, gl_654, gl_656, gl_667, gl_669, \
                         gl_671, gl_673 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_75[k] = -1.23046875 * gl_4[k]
                  - 3.69140625 * gl_11[k]
                  + 9.84375 * gl_13[k]
                  - 3.69140625 * gl_22[k]
                  + 19.6875 * gl_24[k]
                  - 11.8125 * gl_26[k]
                  - 1.23046875 * gl_37[k]
                  + 9.84375 * gl_39[k]
                  - 11.8125 * gl_41[k]
                  + 2.25 * gl_43[k]
                  - 2.4609375 * gl_139[k]
                  - 7.3828125 * gl_146[k]
                  + 19.6875 * gl_148[k]
                  - 7.3828125 * gl_157[k]
                  + 39.375 * gl_159[k]
                  - 23.625 * gl_161[k]
                  - 2.4609375 * gl_172[k]
                  + 19.6875 * gl_174[k]
                  - 23.625 * gl_176[k]
                  + 4.5 * gl_178[k]
                  + 9.84375 * gl_229[k]
                  + 29.53125 * gl_236[k]
                  - 78.75 * gl_238[k]
                  + 29.53125 * gl_247[k]
                  - 157.5 * gl_249[k]
                  + 94.5 * gl_251[k]
                  + 9.84375 * gl_262[k]
                  - 78.75 * gl_264[k]
                  + 94.5 * gl_266[k]
                  - 18.0 * gl_268[k]
                  - 1.23046875 * gl_454[k]
                  - 3.69140625 * gl_461[k]
                  + 9.84375 * gl_463[k]
                  - 3.69140625 * gl_472[k]
                  + 19.6875 * gl_474[k]
                  - 11.8125 * gl_476[k]
                  - 1.23046875 * gl_487[k]
                  + 9.84375 * gl_489[k]
                  - 11.8125 * gl_491[k]
                  + 2.25 * gl_493[k]
                  + 9.84375 * gl_544[k]
                  + 29.53125 * gl_551[k]
                  - 78.75 * gl_553[k]
                  + 29.53125 * gl_562[k]
                  - 157.5 * gl_564[k]
                  + 94.5 * gl_566[k]
                  + 9.84375 * gl_577[k]
                  - 78.75 * gl_579[k]
                  + 94.5 * gl_581[k]
                  - 18.0 * gl_583[k]
                  - 3.28125 * gl_634[k]
                  - 9.84375 * gl_641[k]
                  + 26.25 * gl_643[k]
                  - 9.84375 * gl_652[k]
                  + 52.5 * gl_654[k]
                  - 31.5 * gl_656[k]
                  - 3.28125 * gl_667[k]
                  + 26.25 * gl_669[k]
                  - 31.5 * gl_671[k]
                  + 6.0 * gl_673[k];
    }

#pragma omp simd aligned(gl_0, gl_3, gl_5, gl_10, gl_12, gl_14, gl_21, gl_23, gl_25, gl_27, \
                         gl_36, gl_38, gl_40, gl_42, gl_44, gl_135, gl_138, gl_140, gl_145, \
                         gl_147, gl_149, gl_156, gl_158, gl_160, gl_162, gl_171, gl_173, \
                         gl_175, gl_177, gl_179, gl_225, gl_228, gl_230, gl_235, gl_237, \
                         gl_239, gl_246, gl_248, gl_250, gl_252, gl_261, gl_263, gl_265, \
                         gl_267, gl_269, gl_450, gl_453, gl_455, gl_460, gl_462, gl_464, \
                         gl_471, gl_473, gl_475, gl_477, gl_486, gl_488, gl_490, gl_492, \
                         gl_494, gl_540, gl_543, gl_545, gl_550, gl_552, gl_554, gl_561, \
                         gl_563, gl_565, gl_567, gl_576, gl_578, gl_580, gl_582, gl_584, \
                         gl_630, gl_633, gl_635, gl_640, gl_642, gl_644, gl_651, gl_653, \
                         gl_655, gl_657, gl_666, gl_668, gl_670, gl_672, \
                         gl_674 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_76[k] = 0.1025390625 * gl_0[k]
                  + 0.41015625 * gl_3[k]
                  - 3.28125 * gl_5[k]
                  + 0.615234375 * gl_10[k]
                  - 9.84375 * gl_12[k]
                  + 9.84375 * gl_14[k]
                  + 0.41015625 * gl_21[k]
                  - 9.84375 * gl_23[k]
                  + 19.6875 * gl_25[k]
                  - 5.25 * gl_27[k]
                  + 0.1025390625 * gl_36[k]
                  - 3.28125 * gl_38[k]
                  + 9.84375 * gl_40[k]
                  - 5.25 * gl_42[k]
                  + 0.375 * gl_44[k]
                  + 0.205078125 * gl_135[k]
                  + 0.8203125 * gl_138[k]
                  - 6.5625 * gl_140[k]
                  + 1.23046875 * gl_145[k]
                  - 19.6875 * gl_147[k]
                  + 19.6875 * gl_149[k]
                  + 0.8203125 * gl_156[k]
                  - 19.6875 * gl_158[k]
                  + 39.375 * gl_160[k]
                  - 10.5 * gl_162[k]
                  + 0.205078125 * gl_171[k]
                  - 6.5625 * gl_173[k]
                  + 19.6875 * gl_175[k]
                  - 10.5 * gl_177[k]
                  + 0.75 * gl_179[k]
                  - 0.8203125 * gl_225[k]
                  - 3.28125 * gl_228[k]
                  + 26.25 * gl_230[k]
                  - 4.921875 * gl_235[k]
                  + 78.75 * gl_237[k]
                  - 78.75 * gl_239[k]
                  - 3.28125 * gl_246[k]
                  + 78.75 * gl_248[k]
                  - 157.5 * gl_250[k]
                  + 42.0 * gl_252[k]
                  - 0.8203125 * gl_261[k]
                  + 26.25 * gl_263[k]
                  - 78.75 * gl_265[k]
                  + 42.0 * gl_267[k]
                  - 3.0 * gl_269[k]
                  + 0.1025390625 * gl_450[k]
                  + 0.41015625 * gl_453[k]
                  - 3.28125 * gl_455[k]
                  + 0.615234375 * gl_460[k]
                  - 9.84375 * gl_462[k]
                  + 9.84375 * gl_464[k]
                  + 0.41015625 * gl_471[k]
                  - 9.84375 * gl_473[k]
                  + 19.6875 * gl_475[k]
                  - 5.25 * gl_477[k]
                  + 0.1025390625 * gl_486[k]
                  - 3.28125 * gl_488[k]
                  + 9.84375 * gl_490[k]
                  - 5.25 * gl_492[k]
                  + 0.375 * gl_494[k]
                  - 0.8203125 * gl_540[k]
                  - 3.28125 * gl_543[k]
                  + 26.25 * gl_545[k]
                  - 4.921875 * gl_550[k]
                  + 78.75 * gl_552[k]
                  - 78.75 * gl_554[k]
                  - 3.28125 * gl_561[k]
                  + 78.75 * gl_563[k]
                  - 157.5 * gl_565[k]
                  + 42.0 * gl_567[k]
                  - 0.8203125 * gl_576[k]
                  + 26.25 * gl_578[k]
                  - 78.75 * gl_580[k]
                  + 42.0 * gl_582[k]
                  - 3.0 * gl_584[k]
                  + 0.2734375 * gl_630[k]
                  + 1.09375 * gl_633[k]
                  - 8.75 * gl_635[k]
                  + 1.640625 * gl_640[k]
                  - 26.25 * gl_642[k]
                  + 26.25 * gl_644[k]
                  + 1.09375 * gl_651[k]
                  - 26.25 * gl_653[k]
                  + 52.5 * gl_655[k]
                  - 14.0 * gl_657[k]
                  + 0.2734375 * gl_666[k]
                  - 8.75 * gl_668[k]
                  + 26.25 * gl_670[k]
                  - 14.0 * gl_672[k]
                  + gl_674[k];
    }

#pragma omp simd aligned(gl_2, gl_7, gl_9, gl_16, gl_18, gl_20, gl_29, gl_31, gl_33, gl_35, \
                         gl_137, gl_142, gl_144, gl_151, gl_153, gl_155, gl_164, gl_166, \
                         gl_168, gl_170, gl_227, gl_232, gl_234, gl_241, gl_243, gl_245, \
                         gl_254, gl_256, gl_258, gl_260, gl_452, gl_457, gl_459, gl_466, \
                         gl_468, gl_470, gl_479, gl_481, gl_483, gl_485, gl_542, gl_547, \
                         gl_549, gl_556, gl_558, gl_560, gl_569, gl_571, gl_573, gl_575, \
                         gl_632, gl_637, gl_639, gl_646, gl_648, gl_650, gl_659, gl_661, \
                         gl_663, gl_665 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_77[k] = -1.23046875 * gl_2[k]
                  - 3.69140625 * gl_7[k]
                  + 9.84375 * gl_9[k]
                  - 3.69140625 * gl_16[k]
                  + 19.6875 * gl_18[k]
                  - 11.8125 * gl_20[k]
                  - 1.23046875 * gl_29[k]
                  + 9.84375 * gl_31[k]
                  - 11.8125 * gl_33[k]
                  + 2.25 * gl_35[k]
                  - 2.4609375 * gl_137[k]
                  - 7.3828125 * gl_142[k]
                  + 19.6875 * gl_144[k]
                  - 7.3828125 * gl_151[k]
                  + 39.375 * gl_153[k]
                  - 23.625 * gl_155[k]
                  - 2.4609375 * gl_164[k]
                  + 19.6875 * gl_166[k]
                  - 23.625 * gl_168[k]
                  + 4.5 * gl_170[k]
                  + 9.84375 * gl_227[k]
                  + 29.53125 * gl_232[k]
                  - 78.75 * gl_234[k]
                  + 29.53125 * gl_241[k]
                  - 157.5 * gl_243[k]
                  + 94.5 * gl_245[k]
                  + 9.84375 * gl_254[k]
                  - 78.75 * gl_256[k]
                  + 94.5 * gl_258[k]
                  - 18.0 * gl_260[k]
                  - 1.23046875 * gl_452[k]
                  - 3.69140625 * gl_457[k]
                  + 9.84375 * gl_459[k]
                  - 3.69140625 * gl_466[k]
                  + 19.6875 * gl_468[k]
                  - 11.8125 * gl_470[k]
                  - 1.23046875 * gl_479[k]
                  + 9.84375 * gl_481[k]
                  - 11.8125 * gl_483[k]
                  + 2.25 * gl_485[k]
                  + 9.84375 * gl_542[k]
                  + 29.53125 * gl_547[k]
                  - 78.75 * gl_549[k]
                  + 29.53125 * gl_556[k]
                  - 157.5 * gl_558[k]
                  + 94.5 * gl_560[k]
                  + 9.84375 * gl_569[k]
                  - 78.75 * gl_571[k]
                  + 94.5 * gl_573[k]
                  - 18.0 * gl_575[k]
                  - 3.28125 * gl_632[k]
                  - 9.84375 * gl_637[k]
                  + 26.25 * gl_639[k]
                  - 9.84375 * gl_646[k]
                  + 52.5 * gl_648[k]
                  - 31.5 * gl_650[k]
                  - 3.28125 * gl_659[k]
                  + 26.25 * gl_661[k]
                  - 31.5 * gl_663[k]
                  + 6.0 * gl_665[k];
    }

#pragma omp simd aligned(gl_0, gl_3, gl_5, gl_12, gl_14, gl_21, gl_23, gl_27, gl_36, gl_38, \
                         gl_40, gl_42, gl_135, gl_138, gl_140, gl_147, gl_149, gl_156, gl_158, \
                         gl_162, gl_171, gl_173, gl_175, gl_177, gl_225, gl_228, gl_230, \
                         gl_237, gl_239, gl_246, gl_248, gl_252, gl_261, gl_263, gl_265, \
                         gl_267, gl_450, gl_453, gl_455, gl_462, gl_464, gl_471, gl_473, \
                         gl_477, gl_486, gl_488, gl_490, gl_492, gl_540, gl_543, gl_545, \
                         gl_552, gl_554, gl_561, gl_563, gl_567, gl_576, gl_578, gl_580, \
                         gl_582, gl_630, gl_633, gl_635, gl_642, gl_644, gl_651, gl_653, \
                         gl_657, gl_666, gl_668, gl_670, gl_672 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_78[k] = -f_485 * gl_0[k]
                  - f_465 * gl_3[k]
                  + f_486 * gl_5[k]
                  + f_486 * gl_12[k]
                  - f_487 * gl_14[k]
                  + f_465 * gl_21[k]
                  - f_486 * gl_23[k]
                  + f_488 * gl_27[k]
                  + f_485 * gl_36[k]
                  - f_486 * gl_38[k]
                  + f_487 * gl_40[k]
                  - f_488 * gl_42[k]
                  - f_465 * gl_135[k]
                  - f_471 * gl_138[k]
                  + f_467 * gl_140[k]
                  + f_467 * gl_147[k]
                  - f_469 * gl_149[k]
                  + f_471 * gl_156[k]
                  - f_467 * gl_158[k]
                  + f_470 * gl_162[k]
                  + f_465 * gl_171[k]
                  - f_467 * gl_173[k]
                  + f_469 * gl_175[k]
                  - f_470 * gl_177[k]
                  + f_489 * gl_225[k]
                  + f_476 * gl_228[k]
                  - f_473 * gl_230[k]
                  - f_473 * gl_237[k]
                  + f_490 * gl_239[k]
                  - f_476 * gl_246[k]
                  + f_473 * gl_248[k]
                  - f_114 * gl_252[k]
                  - f_489 * gl_261[k]
                  + f_473 * gl_263[k]
                  - f_490 * gl_265[k]
                  + f_114 * gl_267[k]
                  - f_485 * gl_450[k]
                  - f_465 * gl_453[k]
                  + f_486 * gl_455[k]
                  + f_486 * gl_462[k]
                  - f_487 * gl_464[k]
                  + f_465 * gl_471[k]
                  - f_486 * gl_473[k]
                  + f_488 * gl_477[k]
                  + f_485 * gl_486[k]
                  - f_486 * gl_488[k]
                  + f_487 * gl_490[k]
                  - f_488 * gl_492[k]
                  + f_489 * gl_540[k]
                  + f_476 * gl_543[k]
                  - f_473 * gl_545[k]
                  - f_473 * gl_552[k]
                  + f_490 * gl_554[k]
                  - f_476 * gl_561[k]
                  + f_473 * gl_563[k]
                  - f_114 * gl_567[k]
                  - f_489 * gl_576[k]
                  + f_473 * gl_578[k]
                  - f_490 * gl_580[k]
                  + f_114 * gl_582[k]
                  - f_491 * gl_630[k]
                  - f_482 * gl_633[k]
                  + f_487 * gl_635[k]
                  + f_487 * gl_642[k]
                  - f_492 * gl_644[k]
                  + f_482 * gl_651[k]
                  - f_487 * gl_653[k]
                  + f_119 * gl_657[k]
                  + f_491 * gl_666[k]
                  - f_487 * gl_668[k]
                  + f_492 * gl_670[k]
                  - f_119 * gl_672[k];
    }

#pragma omp simd aligned(gl_2, gl_7, gl_9, gl_16, gl_18, gl_20, gl_29, gl_31, gl_33, gl_137, \
                         gl_142, gl_144, gl_151, gl_153, gl_155, gl_164, gl_166, gl_168, \
                         gl_227, gl_232, gl_234, gl_241, gl_243, gl_245, gl_254, gl_256, \
                         gl_258, gl_452, gl_457, gl_459, gl_466, gl_468, gl_470, gl_479, \
                         gl_481, gl_483, gl_542, gl_547, gl_549, gl_556, gl_558, gl_560, \
                         gl_569, gl_571, gl_573, gl_632, gl_637, gl_639, gl_646, gl_648, \
                         gl_650, gl_659, gl_661, gl_663 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_79[k] = f_442 * gl_2[k]
                  - f_442 * gl_7[k]
                  - f_445 * gl_9[k]
                  - f_440 * gl_16[k]
                  + f_443 * gl_18[k]
                  + f_446 * gl_20[k]
                  - f_439 * gl_29[k]
                  + f_441 * gl_31[k]
                  - f_444 * gl_33[k]
                  + f_450 * gl_137[k]
                  - f_450 * gl_142[k]
                  - f_443 * gl_144[k]
                  - f_448 * gl_151[k]
                  + f_451 * gl_153[k]
                  + f_453 * gl_155[k]
                  - f_447 * gl_164[k]
                  + f_449 * gl_166[k]
                  - f_452 * gl_168[k]
                  - f_456 * gl_227[k]
                  + f_456 * gl_232[k]
                  + f_459 * gl_234[k]
                  + f_449 * gl_241[k]
                  - f_457 * gl_243[k]
                  - f_460 * gl_245[k]
                  + f_454 * gl_254[k]
                  - f_455 * gl_256[k]
                  + f_458 * gl_258[k]
                  + f_442 * gl_452[k]
                  - f_442 * gl_457[k]
                  - f_445 * gl_459[k]
                  - f_440 * gl_466[k]
                  + f_443 * gl_468[k]
                  + f_446 * gl_470[k]
                  - f_439 * gl_479[k]
                  + f_441 * gl_481[k]
                  - f_444 * gl_483[k]
                  - f_456 * gl_542[k]
                  + f_456 * gl_547[k]
                  + f_459 * gl_549[k]
                  + f_449 * gl_556[k]
                  - f_457 * gl_558[k]
                  - f_460 * gl_560[k]
                  + f_454 * gl_569[k]
                  - f_455 * gl_571[k]
                  + f_458 * gl_573[k]
                  + f_461 * gl_632[k]
                  - f_461 * gl_637[k]
                  - f_463 * gl_639[k]
                  - f_443 * gl_646[k]
                  + f_462 * gl_648[k]
                  + f_464 * gl_650[k]
                  - f_456 * gl_659[k]
                  + f_459 * gl_661[k]
                  - f_460 * gl_663[k];
    }

#pragma omp simd aligned(gl_0, gl_3, gl_5, gl_10, gl_12, gl_14, gl_21, gl_23, gl_25, gl_36, \
                         gl_38, gl_40, gl_135, gl_138, gl_140, gl_145, gl_147, gl_149, gl_156, \
                         gl_158, gl_160, gl_171, gl_173, gl_175, gl_225, gl_228, gl_230, \
                         gl_235, gl_237, gl_239, gl_246, gl_248, gl_250, gl_261, gl_263, \
                         gl_265, gl_450, gl_453, gl_455, gl_460, gl_462, gl_464, gl_471, \
                         gl_473, gl_475, gl_486, gl_488, gl_490, gl_540, gl_543, gl_545, \
                         gl_550, gl_552, gl_554, gl_561, gl_563, gl_565, gl_576, gl_578, \
                         gl_580, gl_630, gl_633, gl_635, gl_640, gl_642, gl_644, gl_651, \
                         gl_653, gl_655, gl_666, gl_668, gl_670 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_80[k] = f_493 * gl_0[k]
                  - f_427 * gl_3[k]
                  - f_494 * gl_5[k]
                  - f_495 * gl_10[k]
                  + f_496 * gl_12[k]
                  + f_497 * gl_14[k]
                  - f_427 * gl_21[k]
                  + f_496 * gl_23[k]
                  - f_498 * gl_25[k]
                  + f_493 * gl_36[k]
                  - f_494 * gl_38[k]
                  + f_497 * gl_40[k]
                  + f_499 * gl_135[k]
                  - f_430 * gl_138[k]
                  - f_500 * gl_140[k]
                  - f_501 * gl_145[k]
                  + f_498 * gl_147[k]
                  + f_502 * gl_149[k]
                  - f_430 * gl_156[k]
                  + f_498 * gl_158[k]
                  - f_503 * gl_160[k]
                  + f_499 * gl_171[k]
                  - f_500 * gl_173[k]
                  + f_502 * gl_175[k]
                  - f_430 * gl_225[k]
                  + f_433 * gl_228[k]
                  + f_431 * gl_230[k]
                  + f_502 * gl_235[k]
                  - f_504 * gl_237[k]
                  - f_432 * gl_239[k]
                  + f_433 * gl_246[k]
                  - f_504 * gl_248[k]
                  + f_505 * gl_250[k]
                  - f_430 * gl_261[k]
                  + f_431 * gl_263[k]
                  - f_432 * gl_265[k]
                  + f_493 * gl_450[k]
                  - f_427 * gl_453[k]
                  - f_494 * gl_455[k]
                  - f_495 * gl_460[k]
                  + f_496 * gl_462[k]
                  + f_497 * gl_464[k]
                  - f_427 * gl_471[k]
                  + f_496 * gl_473[k]
                  - f_498 * gl_475[k]
                  + f_493 * gl_486[k]
                  - f_494 * gl_488[k]
                  + f_497 * gl_490[k]
                  - f_430 * gl_540[k]
                  + f_433 * gl_543[k]
                  + f_431 * gl_545[k]
                  + f_502 * gl_550[k]
                  - f_504 * gl_552[k]
                  - f_432 * gl_554[k]
                  + f_433 * gl_561[k]
                  - f_504 * gl_563[k]
                  + f_505 * gl_565[k]
                  - f_430 * gl_576[k]
                  + f_431 * gl_578[k]
                  - f_432 * gl_580[k]
                  + f_506 * gl_630[k]
                  - f_436 * gl_633[k]
                  - f_507 * gl_635[k]
                  - f_508 * gl_640[k]
                  + f_432 * gl_642[k]
                  + f_509 * gl_644[k]
                  - f_436 * gl_651[k]
                  + f_432 * gl_653[k]
                  - f_510 * gl_655[k]
                  + f_506 * gl_666[k]
                  - f_507 * gl_668[k]
                  + f_509 * gl_670[k];
    }

#pragma omp simd aligned(gl_2, gl_7, gl_9, gl_16, gl_18, gl_29, gl_31, gl_137, gl_142, gl_144, \
                         gl_151, gl_153, gl_164, gl_166, gl_227, gl_232, gl_234, gl_241, \
                         gl_243, gl_254, gl_256, gl_452, gl_457, gl_459, gl_466, gl_468, \
                         gl_479, gl_481, gl_542, gl_547, gl_549, gl_556, gl_558, gl_569, \
                         gl_571, gl_632, gl_637, gl_639, gl_646, gl_648, gl_659, \
                         gl_661 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_81[k] = -f_411 * gl_2[k]
                  + f_409 * gl_7[k]
                  + f_412 * gl_9[k]
                  + f_407 * gl_16[k]
                  - f_410 * gl_18[k]
                  - f_407 * gl_29[k]
                  + f_408 * gl_31[k]
                  - f_416 * gl_137[k]
                  + f_414 * gl_142[k]
                  + f_417 * gl_144[k]
                  + f_413 * gl_151[k]
                  - f_415 * gl_153[k]
                  - f_413 * gl_164[k]
                  + f_410 * gl_166[k]
                  + f_417 * gl_227[k]
                  - f_419 * gl_232[k]
                  - f_421 * gl_234[k]
                  - f_410 * gl_241[k]
                  + f_420 * gl_243[k]
                  + f_410 * gl_254[k]
                  - f_418 * gl_256[k]
                  - f_411 * gl_452[k]
                  + f_409 * gl_457[k]
                  + f_412 * gl_459[k]
                  + f_407 * gl_466[k]
                  - f_410 * gl_468[k]
                  - f_407 * gl_479[k]
                  + f_408 * gl_481[k]
                  + f_417 * gl_542[k]
                  - f_419 * gl_547[k]
                  - f_421 * gl_549[k]
                  - f_410 * gl_556[k]
                  + f_420 * gl_558[k]
                  + f_410 * gl_569[k]
                  - f_418 * gl_571[k]
                  - f_425 * gl_632[k]
                  + f_423 * gl_637[k]
                  + f_426 * gl_639[k]
                  + f_0 * gl_646[k]
                  - f_424 * gl_648[k]
                  - f_0 * gl_659[k]
                  + f_422 * gl_661[k];
    }

#pragma omp simd aligned(gl_0, gl_3, gl_5, gl_12, gl_21, gl_23, gl_36, gl_38, gl_135, gl_138, \
                         gl_140, gl_147, gl_156, gl_158, gl_171, gl_173, gl_225, gl_228, \
                         gl_230, gl_237, gl_246, gl_248, gl_261, gl_263, gl_450, gl_453, \
                         gl_455, gl_462, gl_471, gl_473, gl_486, gl_488, gl_540, gl_543, \
                         gl_545, gl_552, gl_561, gl_563, gl_576, gl_578, gl_630, gl_633, \
                         gl_635, gl_642, gl_651, gl_653, gl_666, \
                         gl_668 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_82[k] = -f_511 * gl_0[k]
                  + f_392 * gl_3[k]
                  + f_392 * gl_5[k]
                  - f_512 * gl_12[k]
                  - f_392 * gl_21[k]
                  + f_512 * gl_23[k]
                  + f_511 * gl_36[k]
                  - f_392 * gl_38[k]
                  - f_513 * gl_135[k]
                  + f_396 * gl_138[k]
                  + f_396 * gl_140[k]
                  - f_514 * gl_147[k]
                  - f_396 * gl_156[k]
                  + f_514 * gl_158[k]
                  + f_513 * gl_171[k]
                  - f_396 * gl_173[k]
                  + f_515 * gl_225[k]
                  - f_400 * gl_228[k]
                  - f_400 * gl_230[k]
                  + f_516 * gl_237[k]
                  + f_400 * gl_246[k]
                  - f_516 * gl_248[k]
                  - f_515 * gl_261[k]
                  + f_400 * gl_263[k]
                  - f_511 * gl_450[k]
                  + f_392 * gl_453[k]
                  + f_392 * gl_455[k]
                  - f_512 * gl_462[k]
                  - f_392 * gl_471[k]
                  + f_512 * gl_473[k]
                  + f_511 * gl_486[k]
                  - f_392 * gl_488[k]
                  + f_515 * gl_540[k]
                  - f_400 * gl_543[k]
                  - f_400 * gl_545[k]
                  + f_516 * gl_552[k]
                  + f_400 * gl_561[k]
                  - f_516 * gl_563[k]
                  - f_515 * gl_576[k]
                  + f_400 * gl_578[k]
                  - f_517 * gl_630[k]
                  + f_404 * gl_633[k]
                  + f_404 * gl_635[k]
                  - f_398 * gl_642[k]
                  - f_404 * gl_651[k]
                  + f_398 * gl_653[k]
                  + f_517 * gl_666[k]
                  - f_404 * gl_668[k];
    }

#pragma omp simd aligned(gl_2, gl_7, gl_16, gl_29, gl_137, gl_142, gl_151, gl_164, gl_227, \
                         gl_232, gl_241, gl_254, gl_452, gl_457, gl_466, gl_479, gl_542, \
                         gl_547, gl_556, gl_569, gl_632, gl_637, gl_646, \
                         gl_659 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_83[k] = f_381 * gl_2[k]
                  - f_380 * gl_7[k]
                  + f_379 * gl_16[k]
                  - f_378 * gl_29[k]
                  + f_371 * gl_137[k]
                  - f_383 * gl_142[k]
                  + f_382 * gl_151[k]
                  - f_372 * gl_164[k]
                  - f_387 * gl_227[k]
                  + f_386 * gl_232[k]
                  - f_385 * gl_241[k]
                  + f_384 * gl_254[k]
                  + f_381 * gl_452[k]
                  - f_380 * gl_457[k]
                  + f_379 * gl_466[k]
                  - f_378 * gl_479[k]
                  - f_387 * gl_542[k]
                  + f_386 * gl_547[k]
                  - f_385 * gl_556[k]
                  + f_384 * gl_569[k]
                  + f_390 * gl_632[k]
                  - f_384 * gl_637[k]
                  + f_389 * gl_646[k]
                  - f_388 * gl_659[k];
    }

#pragma omp simd aligned(gl_0, gl_3, gl_10, gl_21, gl_36, gl_135, gl_138, gl_145, gl_156, \
                         gl_171, gl_225, gl_228, gl_235, gl_246, gl_261, gl_450, gl_453, \
                         gl_460, gl_471, gl_486, gl_540, gl_543, gl_550, gl_561, gl_576, \
                         gl_630, gl_633, gl_640, gl_651, gl_666 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_84[k] = f_518 * gl_0[k]
                  - f_378 * gl_3[k]
                  + f_519 * gl_10[k]
                  - f_378 * gl_21[k]
                  + f_518 * gl_36[k]
                  + f_520 * gl_135[k]
                  - f_372 * gl_138[k]
                  + f_379 * gl_145[k]
                  - f_372 * gl_156[k]
                  + f_520 * gl_171[k]
                  - f_371 * gl_225[k]
                  + f_384 * gl_228[k]
                  - f_521 * gl_235[k]
                  + f_384 * gl_246[k]
                  - f_371 * gl_261[k]
                  + f_518 * gl_450[k]
                  - f_378 * gl_453[k]
                  + f_519 * gl_460[k]
                  - f_378 * gl_471[k]
                  + f_518 * gl_486[k]
                  - f_371 * gl_540[k]
                  + f_384 * gl_543[k]
                  - f_521 * gl_550[k]
                  + f_384 * gl_561[k]
                  - f_371 * gl_576[k]
                  + f_522 * gl_630[k]
                  - f_388 * gl_633[k]
                  + f_10 * gl_640[k]
                  - f_388 * gl_651[k]
                  + f_522 * gl_666[k];
    }

#pragma omp simd aligned(gl_91, gl_96, gl_105, gl_118, gl_316, gl_321, gl_330, gl_343, gl_406, \
                         gl_411, gl_420, gl_433 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_85[k] = -f_258 * gl_91[k]
                  + f_259 * gl_96[k]
                  - f_259 * gl_105[k]
                  + f_258 * gl_118[k]
                  - f_258 * gl_316[k]
                  + f_259 * gl_321[k]
                  - f_259 * gl_330[k]
                  + f_258 * gl_343[k]
                  + f_260 * gl_406[k]
                  - f_261 * gl_411[k]
                  + f_261 * gl_420[k]
                  - f_260 * gl_433[k];
    }

#pragma omp simd aligned(gl_94, gl_101, gl_112, gl_127, gl_319, gl_326, gl_337, gl_352, \
                         gl_409, gl_416, gl_427, gl_442 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_86[k] = -f_262 * gl_94[k]
                  + f_263 * gl_101[k]
                  - f_264 * gl_112[k]
                  + f_265 * gl_127[k]
                  - f_262 * gl_319[k]
                  + f_263 * gl_326[k]
                  - f_264 * gl_337[k]
                  + f_265 * gl_352[k]
                  + f_266 * gl_409[k]
                  - f_267 * gl_416[k]
                  + f_268 * gl_427[k]
                  - f_269 * gl_442[k];
    }

#pragma omp simd aligned(gl_91, gl_96, gl_98, gl_105, gl_107, gl_118, gl_120, gl_316, gl_321, \
                         gl_323, gl_330, gl_332, gl_343, gl_345, gl_406, gl_411, gl_413, \
                         gl_420, gl_422, gl_433, gl_435 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_87[k] = f_270 * gl_91[k]
                  - f_271 * gl_96[k]
                  - f_272 * gl_98[k]
                  - f_271 * gl_105[k]
                  + f_273 * gl_107[k]
                  + f_270 * gl_118[k]
                  - f_272 * gl_120[k]
                  + f_270 * gl_316[k]
                  - f_271 * gl_321[k]
                  - f_272 * gl_323[k]
                  - f_271 * gl_330[k]
                  + f_273 * gl_332[k]
                  + f_270 * gl_343[k]
                  - f_272 * gl_345[k]
                  - f_274 * gl_406[k]
                  + f_275 * gl_411[k]
                  + f_276 * gl_413[k]
                  + f_275 * gl_420[k]
                  - f_277 * gl_422[k]
                  - f_274 * gl_433[k]
                  + f_276 * gl_435[k];
    }

#pragma omp simd aligned(gl_94, gl_101, gl_103, gl_112, gl_114, gl_127, gl_129, gl_319, \
                         gl_326, gl_328, gl_337, gl_339, gl_352, gl_354, gl_409, gl_416, \
                         gl_418, gl_427, gl_429, gl_442, gl_444 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_88[k] = f_278 * gl_94[k]
                  - f_278 * gl_101[k]
                  - f_279 * gl_103[k]
                  - f_280 * gl_112[k]
                  + f_281 * gl_114[k]
                  + f_282 * gl_127[k]
                  - f_283 * gl_129[k]
                  + f_278 * gl_319[k]
                  - f_278 * gl_326[k]
                  - f_279 * gl_328[k]
                  - f_280 * gl_337[k]
                  + f_281 * gl_339[k]
                  + f_282 * gl_352[k]
                  - f_283 * gl_354[k]
                  - f_284 * gl_409[k]
                  + f_284 * gl_416[k]
                  + f_285 * gl_418[k]
                  + f_286 * gl_427[k]
                  - f_287 * gl_429[k]
                  - f_288 * gl_442[k]
                  + f_289 * gl_444[k];
    }

#pragma omp simd aligned(gl_91, gl_96, gl_98, gl_105, gl_109, gl_118, gl_120, gl_122, gl_316, \
                         gl_321, gl_323, gl_330, gl_334, gl_343, gl_345, gl_347, gl_406, \
                         gl_411, gl_413, gl_420, gl_424, gl_433, gl_435, \
                         gl_437 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_89[k] = -f_290 * gl_91[k]
                  - f_290 * gl_96[k]
                  + f_291 * gl_98[k]
                  + f_290 * gl_105[k]
                  - f_292 * gl_109[k]
                  + f_290 * gl_118[k]
                  - f_291 * gl_120[k]
                  + f_292 * gl_122[k]
                  - f_290 * gl_316[k]
                  - f_290 * gl_321[k]
                  + f_291 * gl_323[k]
                  + f_290 * gl_330[k]
                  - f_292 * gl_334[k]
                  + f_290 * gl_343[k]
                  - f_291 * gl_345[k]
                  + f_292 * gl_347[k]
                  + f_293 * gl_406[k]
                  + f_293 * gl_411[k]
                  - f_294 * gl_413[k]
                  - f_293 * gl_420[k]
                  + f_295 * gl_424[k]
                  - f_293 * gl_433[k]
                  + f_294 * gl_435[k]
                  - f_295 * gl_437[k];
    }

#pragma omp simd aligned(gl_94, gl_101, gl_103, gl_112, gl_114, gl_116, gl_127, gl_129, \
                         gl_131, gl_319, gl_326, gl_328, gl_337, gl_339, gl_341, gl_352, \
                         gl_354, gl_356, gl_409, gl_416, gl_418, gl_427, gl_429, gl_431, \
                         gl_442, gl_444, gl_446 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_90[k] = -f_296 * gl_94[k]
                  - f_297 * gl_101[k]
                  + f_298 * gl_103[k]
                  - f_299 * gl_112[k]
                  + f_300 * gl_114[k]
                  - f_301 * gl_116[k]
                  + f_299 * gl_127[k]
                  - f_302 * gl_129[k]
                  + f_303 * gl_131[k]
                  - f_296 * gl_319[k]
                  - f_297 * gl_326[k]
                  + f_298 * gl_328[k]
                  - f_299 * gl_337[k]
                  + f_300 * gl_339[k]
                  - f_301 * gl_341[k]
                  + f_299 * gl_352[k]
                  - f_302 * gl_354[k]
                  + f_303 * gl_356[k]
                  + f_304 * gl_409[k]
                  + f_302 * gl_416[k]
                  - f_305 * gl_418[k]
                  + f_306 * gl_427[k]
                  - f_307 * gl_429[k]
                  + f_308 * gl_431[k]
                  - f_306 * gl_442[k]
                  + f_309 * gl_444[k]
                  - f_310 * gl_446[k];
    }

#pragma omp simd aligned(gl_91, gl_96, gl_98, gl_105, gl_107, gl_109, gl_118, gl_120, gl_122, \
                         gl_124, gl_316, gl_321, gl_323, gl_330, gl_332, gl_334, gl_343, \
                         gl_345, gl_347, gl_349, gl_406, gl_411, gl_413, gl_420, gl_422, \
                         gl_424, gl_433, gl_435, gl_437, gl_439 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_91[k] = f_311 * gl_91[k]
                  + f_312 * gl_96[k]
                  - f_313 * gl_98[k]
                  + f_312 * gl_105[k]
                  - f_314 * gl_107[k]
                  + f_315 * gl_109[k]
                  + f_311 * gl_118[k]
                  - f_313 * gl_120[k]
                  + f_315 * gl_122[k]
                  - f_316 * gl_124[k]
                  + f_311 * gl_316[k]
                  + f_312 * gl_321[k]
                  - f_313 * gl_323[k]
                  + f_312 * gl_330[k]
                  - f_314 * gl_332[k]
                  + f_315 * gl_334[k]
                  + f_311 * gl_343[k]
                  - f_313 * gl_345[k]
                  + f_315 * gl_347[k]
                  - f_316 * gl_349[k]
                  - f_317 * gl_406[k]
                  - f_318 * gl_411[k]
                  + f_319 * gl_413[k]
                  - f_318 * gl_420[k]
                  + f_315 * gl_422[k]
                  - f_320 * gl_424[k]
                  - f_317 * gl_433[k]
                  + f_319 * gl_435[k]
                  - f_320 * gl_437[k]
                  + f_321 * gl_439[k];
    }

#pragma omp simd aligned(gl_94, gl_101, gl_103, gl_112, gl_114, gl_116, gl_127, gl_129, \
                         gl_131, gl_133, gl_319, gl_326, gl_328, gl_337, gl_339, gl_341, \
                         gl_352, gl_354, gl_356, gl_358, gl_409, gl_416, gl_418, gl_427, \
                         gl_429, gl_431, gl_442, gl_444, gl_446, \
                         gl_448 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_92[k] = f_322 * gl_94[k]
                  + f_323 * gl_101[k]
                  - f_324 * gl_103[k]
                  + f_323 * gl_112[k]
                  - f_325 * gl_114[k]
                  + f_326 * gl_116[k]
                  + f_322 * gl_127[k]
                  - f_324 * gl_129[k]
                  + f_326 * gl_131[k]
                  - f_327 * gl_133[k]
                  + f_322 * gl_319[k]
                  + f_323 * gl_326[k]
                  - f_324 * gl_328[k]
                  + f_323 * gl_337[k]
                  - f_325 * gl_339[k]
                  + f_326 * gl_341[k]
                  + f_322 * gl_352[k]
                  - f_324 * gl_354[k]
                  + f_326 * gl_356[k]
                  - f_327 * gl_358[k]
                  - f_328 * gl_409[k]
                  - f_329 * gl_416[k]
                  + f_330 * gl_418[k]
                  - f_329 * gl_427[k]
                  + f_331 * gl_429[k]
                  - f_332 * gl_431[k]
                  - f_328 * gl_442[k]
                  + f_330 * gl_444[k]
                  - f_332 * gl_446[k]
                  + f_333 * gl_448[k];
    }

#pragma omp simd aligned(gl_90, gl_93, gl_95, gl_100, gl_102, gl_104, gl_111, gl_113, gl_115, \
                         gl_117, gl_126, gl_128, gl_130, gl_132, gl_134, gl_315, gl_318, \
                         gl_320, gl_325, gl_327, gl_329, gl_336, gl_338, gl_340, gl_342, \
                         gl_351, gl_353, gl_355, gl_357, gl_359, gl_405, gl_408, gl_410, \
                         gl_415, gl_417, gl_419, gl_426, gl_428, gl_430, gl_432, gl_441, \
                         gl_443, gl_445, gl_447, gl_449 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_93[k] = -f_334 * gl_90[k]
                  - f_335 * gl_93[k]
                  + f_336 * gl_95[k]
                  - f_337 * gl_100[k]
                  + f_324 * gl_102[k]
                  - f_324 * gl_104[k]
                  - f_335 * gl_111[k]
                  + f_324 * gl_113[k]
                  - f_325 * gl_115[k]
                  + f_338 * gl_117[k]
                  - f_334 * gl_126[k]
                  + f_336 * gl_128[k]
                  - f_324 * gl_130[k]
                  + f_338 * gl_132[k]
                  - f_339 * gl_134[k]
                  - f_334 * gl_315[k]
                  - f_335 * gl_318[k]
                  + f_336 * gl_320[k]
                  - f_337 * gl_325[k]
                  + f_324 * gl_327[k]
                  - f_324 * gl_329[k]
                  - f_335 * gl_336[k]
                  + f_324 * gl_338[k]
                  - f_325 * gl_340[k]
                  + f_338 * gl_342[k]
                  - f_334 * gl_351[k]
                  + f_336 * gl_353[k]
                  - f_324 * gl_355[k]
                  + f_338 * gl_357[k]
                  - f_339 * gl_359[k]
                  + f_340 * gl_405[k]
                  + f_341 * gl_408[k]
                  - f_342 * gl_410[k]
                  + f_343 * gl_415[k]
                  - f_330 * gl_417[k]
                  + f_330 * gl_419[k]
                  + f_341 * gl_426[k]
                  - f_330 * gl_428[k]
                  + f_331 * gl_430[k]
                  - f_344 * gl_432[k]
                  + f_340 * gl_441[k]
                  - f_342 * gl_443[k]
                  + f_330 * gl_445[k]
                  - f_344 * gl_447[k]
                  + f_345 * gl_449[k];
    }

#pragma omp simd aligned(gl_92, gl_97, gl_99, gl_106, gl_108, gl_110, gl_119, gl_121, gl_123, \
                         gl_125, gl_317, gl_322, gl_324, gl_331, gl_333, gl_335, gl_344, \
                         gl_346, gl_348, gl_350, gl_407, gl_412, gl_414, gl_421, gl_423, \
                         gl_425, gl_434, gl_436, gl_438, gl_440 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_94[k] = f_322 * gl_92[k]
                  + f_323 * gl_97[k]
                  - f_324 * gl_99[k]
                  + f_323 * gl_106[k]
                  - f_325 * gl_108[k]
                  + f_326 * gl_110[k]
                  + f_322 * gl_119[k]
                  - f_324 * gl_121[k]
                  + f_326 * gl_123[k]
                  - f_327 * gl_125[k]
                  + f_322 * gl_317[k]
                  + f_323 * gl_322[k]
                  - f_324 * gl_324[k]
                  + f_323 * gl_331[k]
                  - f_325 * gl_333[k]
                  + f_326 * gl_335[k]
                  + f_322 * gl_344[k]
                  - f_324 * gl_346[k]
                  + f_326 * gl_348[k]
                  - f_327 * gl_350[k]
                  - f_328 * gl_407[k]
                  - f_329 * gl_412[k]
                  + f_330 * gl_414[k]
                  - f_329 * gl_421[k]
                  + f_331 * gl_423[k]
                  - f_332 * gl_425[k]
                  - f_328 * gl_434[k]
                  + f_330 * gl_436[k]
                  - f_332 * gl_438[k]
                  + f_333 * gl_440[k];
    }

#pragma omp simd aligned(gl_90, gl_93, gl_95, gl_102, gl_104, gl_111, gl_113, gl_117, gl_126, \
                         gl_128, gl_130, gl_132, gl_315, gl_318, gl_320, gl_327, gl_329, \
                         gl_336, gl_338, gl_342, gl_351, gl_353, gl_355, gl_357, gl_405, \
                         gl_408, gl_410, gl_417, gl_419, gl_426, gl_428, gl_432, gl_441, \
                         gl_443, gl_445, gl_447 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_95[k] = f_346 * gl_90[k]
                  + f_311 * gl_93[k]
                  - f_347 * gl_95[k]
                  - f_347 * gl_102[k]
                  + f_319 * gl_104[k]
                  - f_311 * gl_111[k]
                  + f_347 * gl_113[k]
                  - f_348 * gl_117[k]
                  - f_346 * gl_126[k]
                  + f_347 * gl_128[k]
                  - f_319 * gl_130[k]
                  + f_348 * gl_132[k]
                  + f_346 * gl_315[k]
                  + f_311 * gl_318[k]
                  - f_347 * gl_320[k]
                  - f_347 * gl_327[k]
                  + f_319 * gl_329[k]
                  - f_311 * gl_336[k]
                  + f_347 * gl_338[k]
                  - f_348 * gl_342[k]
                  - f_346 * gl_351[k]
                  + f_347 * gl_353[k]
                  - f_319 * gl_355[k]
                  + f_348 * gl_357[k]
                  - f_349 * gl_405[k]
                  - f_317 * gl_408[k]
                  + f_350 * gl_410[k]
                  + f_350 * gl_417[k]
                  - f_351 * gl_419[k]
                  + f_317 * gl_426[k]
                  - f_350 * gl_428[k]
                  + f_352 * gl_432[k]
                  + f_349 * gl_441[k]
                  - f_350 * gl_443[k]
                  + f_351 * gl_445[k]
                  - f_352 * gl_447[k];
    }

#pragma omp simd aligned(gl_92, gl_97, gl_99, gl_106, gl_108, gl_110, gl_119, gl_121, gl_123, \
                         gl_317, gl_322, gl_324, gl_331, gl_333, gl_335, gl_344, gl_346, \
                         gl_348, gl_407, gl_412, gl_414, gl_421, gl_423, gl_425, gl_434, \
                         gl_436, gl_438 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_96[k] = -f_299 * gl_92[k]
                  + f_299 * gl_97[k]
                  + f_302 * gl_99[k]
                  + f_297 * gl_106[k]
                  - f_300 * gl_108[k]
                  - f_303 * gl_110[k]
                  + f_296 * gl_119[k]
                  - f_298 * gl_121[k]
                  + f_301 * gl_123[k]
                  - f_299 * gl_317[k]
                  + f_299 * gl_322[k]
                  + f_302 * gl_324[k]
                  + f_297 * gl_331[k]
                  - f_300 * gl_333[k]
                  - f_303 * gl_335[k]
                  + f_296 * gl_344[k]
                  - f_298 * gl_346[k]
                  + f_301 * gl_348[k]
                  + f_306 * gl_407[k]
                  - f_306 * gl_412[k]
                  - f_309 * gl_414[k]
                  - f_302 * gl_421[k]
                  + f_307 * gl_423[k]
                  + f_310 * gl_425[k]
                  - f_304 * gl_434[k]
                  + f_305 * gl_436[k]
                  - f_308 * gl_438[k];
    }

#pragma omp simd aligned(gl_90, gl_93, gl_95, gl_100, gl_102, gl_104, gl_111, gl_113, gl_115, \
                         gl_126, gl_128, gl_130, gl_315, gl_318, gl_320, gl_325, gl_327, \
                         gl_329, gl_336, gl_338, gl_340, gl_351, gl_353, gl_355, gl_405, \
                         gl_408, gl_410, gl_415, gl_417, gl_419, gl_426, gl_428, gl_430, \
                         gl_441, gl_443, gl_445 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_97[k] = -f_353 * gl_90[k]
                  + f_290 * gl_93[k]
                  + f_354 * gl_95[k]
                  + f_355 * gl_100[k]
                  - f_356 * gl_102[k]
                  - f_357 * gl_104[k]
                  + f_290 * gl_111[k]
                  - f_356 * gl_113[k]
                  + f_358 * gl_115[k]
                  - f_353 * gl_126[k]
                  + f_354 * gl_128[k]
                  - f_357 * gl_130[k]
                  - f_353 * gl_315[k]
                  + f_290 * gl_318[k]
                  + f_354 * gl_320[k]
                  + f_355 * gl_325[k]
                  - f_356 * gl_327[k]
                  - f_357 * gl_329[k]
                  + f_290 * gl_336[k]
                  - f_356 * gl_338[k]
                  + f_358 * gl_340[k]
                  - f_353 * gl_351[k]
                  + f_354 * gl_353[k]
                  - f_357 * gl_355[k]
                  + f_359 * gl_405[k]
                  - f_293 * gl_408[k]
                  - f_360 * gl_410[k]
                  - f_361 * gl_415[k]
                  + f_292 * gl_417[k]
                  + f_362 * gl_419[k]
                  - f_293 * gl_426[k]
                  + f_292 * gl_428[k]
                  - f_363 * gl_430[k]
                  + f_359 * gl_441[k]
                  - f_360 * gl_443[k]
                  + f_362 * gl_445[k];
    }

#pragma omp simd aligned(gl_92, gl_97, gl_99, gl_106, gl_108, gl_119, gl_121, gl_317, gl_322, \
                         gl_324, gl_331, gl_333, gl_344, gl_346, gl_407, gl_412, gl_414, \
                         gl_421, gl_423, gl_434, gl_436 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_98[k] = f_282 * gl_92[k]
                  - f_280 * gl_97[k]
                  - f_283 * gl_99[k]
                  - f_278 * gl_106[k]
                  + f_281 * gl_108[k]
                  + f_278 * gl_119[k]
                  - f_279 * gl_121[k]
                  + f_282 * gl_317[k]
                  - f_280 * gl_322[k]
                  - f_283 * gl_324[k]
                  - f_278 * gl_331[k]
                  + f_281 * gl_333[k]
                  + f_278 * gl_344[k]
                  - f_279 * gl_346[k]
                  - f_288 * gl_407[k]
                  + f_286 * gl_412[k]
                  + f_289 * gl_414[k]
                  + f_284 * gl_421[k]
                  - f_287 * gl_423[k]
                  - f_284 * gl_434[k]
                  + f_285 * gl_436[k];
    }

#pragma omp simd aligned(gl_90, gl_93, gl_95, gl_102, gl_111, gl_113, gl_126, gl_128, gl_315, \
                         gl_318, gl_320, gl_327, gl_336, gl_338, gl_351, gl_353, gl_405, \
                         gl_408, gl_410, gl_417, gl_426, gl_428, gl_441, \
                         gl_443 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_99[k] = f_364 * gl_90[k]
                  - f_271 * gl_93[k]
                  - f_271 * gl_95[k]
                  + f_365 * gl_102[k]
                  + f_271 * gl_111[k]
                  - f_365 * gl_113[k]
                  - f_364 * gl_126[k]
                  + f_271 * gl_128[k]
                  + f_364 * gl_315[k]
                  - f_271 * gl_318[k]
                  - f_271 * gl_320[k]
                  + f_365 * gl_327[k]
                  + f_271 * gl_336[k]
                  - f_365 * gl_338[k]
                  - f_364 * gl_351[k]
                  + f_271 * gl_353[k]
                  - f_366 * gl_405[k]
                  + f_275 * gl_408[k]
                  + f_275 * gl_410[k]
                  - f_273 * gl_417[k]
                  - f_275 * gl_426[k]
                  + f_273 * gl_428[k]
                  + f_366 * gl_441[k]
                  - f_275 * gl_443[k];
    }

#pragma omp simd aligned(gl_92, gl_97, gl_106, gl_119, gl_317, gl_322, gl_331, gl_344, gl_407, \
                         gl_412, gl_421, gl_434 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_100[k] = -f_265 * gl_92[k]
                   + f_264 * gl_97[k]
                   - f_263 * gl_106[k]
                   + f_262 * gl_119[k]
                   - f_265 * gl_317[k]
                   + f_264 * gl_322[k]
                   - f_263 * gl_331[k]
                   + f_262 * gl_344[k]
                   + f_269 * gl_407[k]
                   - f_268 * gl_412[k]
                   + f_267 * gl_421[k]
                   - f_266 * gl_434[k];
    }

#pragma omp simd aligned(gl_90, gl_93, gl_100, gl_111, gl_126, gl_315, gl_318, gl_325, gl_336, \
                         gl_351, gl_405, gl_408, gl_415, gl_426, \
                         gl_441 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_101[k] = -f_367 * gl_90[k]
                   + f_262 * gl_93[k]
                   - f_368 * gl_100[k]
                   + f_262 * gl_111[k]
                   - f_367 * gl_126[k]
                   - f_367 * gl_315[k]
                   + f_262 * gl_318[k]
                   - f_368 * gl_325[k]
                   + f_262 * gl_336[k]
                   - f_367 * gl_351[k]
                   + f_369 * gl_405[k]
                   - f_266 * gl_408[k]
                   + f_370 * gl_415[k]
                   - f_266 * gl_426[k]
                   + f_369 * gl_441[k];
    }

#pragma omp simd aligned(gl_1, gl_6, gl_15, gl_28, gl_226, gl_231, gl_240, gl_253, gl_451, \
                         gl_456, gl_465, gl_478, gl_541, gl_546, gl_555, \
                         gl_568 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_102[k] = -f_156 * gl_1[k]
                   + f_153 * gl_6[k]
                   - f_153 * gl_15[k]
                   + f_156 * gl_28[k]
                   + f_160 * gl_226[k]
                   - f_157 * gl_231[k]
                   + f_157 * gl_240[k]
                   - f_160 * gl_253[k]
                   + f_156 * gl_451[k]
                   - f_153 * gl_456[k]
                   + f_153 * gl_465[k]
                   - f_156 * gl_478[k]
                   - f_160 * gl_541[k]
                   + f_157 * gl_546[k]
                   - f_157 * gl_555[k]
                   + f_160 * gl_568[k];
    }

#pragma omp simd aligned(gl_4, gl_11, gl_22, gl_37, gl_229, gl_236, gl_247, gl_262, gl_454, \
                         gl_461, gl_472, gl_487, gl_544, gl_551, gl_562, \
                         gl_577 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_103[k] = -f_523 * gl_4[k]
                   + f_255 * gl_11[k]
                   - f_524 * gl_22[k]
                   + f_525 * gl_37[k]
                   + f_155 * gl_229[k]
                   - f_257 * gl_236[k]
                   + f_526 * gl_247[k]
                   - f_527 * gl_262[k]
                   + f_523 * gl_454[k]
                   - f_255 * gl_461[k]
                   + f_524 * gl_472[k]
                   - f_525 * gl_487[k]
                   - f_155 * gl_544[k]
                   + f_257 * gl_551[k]
                   - f_526 * gl_562[k]
                   + f_527 * gl_577[k];
    }

#pragma omp simd aligned(gl_1, gl_6, gl_8, gl_15, gl_17, gl_28, gl_30, gl_226, gl_231, gl_233, \
                         gl_240, gl_242, gl_253, gl_255, gl_451, gl_456, gl_458, gl_465, \
                         gl_467, gl_478, gl_480, gl_541, gl_546, gl_548, gl_555, gl_557, \
                         gl_568, gl_570 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_104[k] = f_528 * gl_1[k]
                   - f_529 * gl_6[k]
                   - f_530 * gl_8[k]
                   - f_529 * gl_15[k]
                   + f_531 * gl_17[k]
                   + f_528 * gl_28[k]
                   - f_530 * gl_30[k]
                   - f_532 * gl_226[k]
                   + f_530 * gl_231[k]
                   + f_533 * gl_233[k]
                   + f_530 * gl_240[k]
                   - f_534 * gl_242[k]
                   - f_532 * gl_253[k]
                   + f_533 * gl_255[k]
                   - f_528 * gl_451[k]
                   + f_529 * gl_456[k]
                   + f_530 * gl_458[k]
                   + f_529 * gl_465[k]
                   - f_531 * gl_467[k]
                   - f_528 * gl_478[k]
                   + f_530 * gl_480[k]
                   + f_532 * gl_541[k]
                   - f_530 * gl_546[k]
                   - f_533 * gl_548[k]
                   - f_530 * gl_555[k]
                   + f_534 * gl_557[k]
                   + f_532 * gl_568[k]
                   - f_533 * gl_570[k];
    }

#pragma omp simd aligned(gl_4, gl_11, gl_13, gl_22, gl_24, gl_37, gl_39, gl_229, gl_236, \
                         gl_238, gl_247, gl_249, gl_262, gl_264, gl_454, gl_461, gl_463, \
                         gl_472, gl_474, gl_487, gl_489, gl_544, gl_551, gl_553, gl_562, \
                         gl_564, gl_577, gl_579 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_105[k] = f_535 * gl_4[k]
                   - f_535 * gl_11[k]
                   - f_536 * gl_13[k]
                   - f_537 * gl_22[k]
                   + f_169 * gl_24[k]
                   + f_538 * gl_37[k]
                   - f_539 * gl_39[k]
                   - f_540 * gl_229[k]
                   + f_540 * gl_236[k]
                   + f_541 * gl_238[k]
                   + f_542 * gl_247[k]
                   - f_175 * gl_249[k]
                   - f_543 * gl_262[k]
                   + f_544 * gl_264[k]
                   - f_535 * gl_454[k]
                   + f_535 * gl_461[k]
                   + f_536 * gl_463[k]
                   + f_537 * gl_472[k]
                   - f_169 * gl_474[k]
                   - f_538 * gl_487[k]
                   + f_539 * gl_489[k]
                   + f_540 * gl_544[k]
                   - f_540 * gl_551[k]
                   - f_541 * gl_553[k]
                   - f_542 * gl_562[k]
                   + f_175 * gl_564[k]
                   + f_543 * gl_577[k]
                   - f_544 * gl_579[k];
    }

#pragma omp simd aligned(gl_1, gl_6, gl_8, gl_15, gl_19, gl_28, gl_30, gl_32, gl_226, gl_231, \
                         gl_233, gl_240, gl_244, gl_253, gl_255, gl_257, gl_451, gl_456, \
                         gl_458, gl_465, gl_469, gl_478, gl_480, gl_482, gl_541, gl_546, \
                         gl_548, gl_555, gl_559, gl_568, gl_570, \
                         gl_572 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_106[k] = -f_545 * gl_1[k]
                   - f_545 * gl_6[k]
                   + f_546 * gl_8[k]
                   + f_545 * gl_15[k]
                   - f_547 * gl_19[k]
                   + f_545 * gl_28[k]
                   - f_546 * gl_30[k]
                   + f_547 * gl_32[k]
                   + f_548 * gl_226[k]
                   + f_548 * gl_231[k]
                   - f_549 * gl_233[k]
                   - f_548 * gl_240[k]
                   + f_550 * gl_244[k]
                   - f_548 * gl_253[k]
                   + f_549 * gl_255[k]
                   - f_550 * gl_257[k]
                   + f_545 * gl_451[k]
                   + f_545 * gl_456[k]
                   - f_546 * gl_458[k]
                   - f_545 * gl_465[k]
                   + f_547 * gl_469[k]
                   - f_545 * gl_478[k]
                   + f_546 * gl_480[k]
                   - f_547 * gl_482[k]
                   - f_548 * gl_541[k]
                   - f_548 * gl_546[k]
                   + f_549 * gl_548[k]
                   + f_548 * gl_555[k]
                   - f_550 * gl_559[k]
                   + f_548 * gl_568[k]
                   - f_549 * gl_570[k]
                   + f_550 * gl_572[k];
    }

#pragma omp simd aligned(gl_4, gl_11, gl_13, gl_22, gl_24, gl_26, gl_37, gl_39, gl_41, gl_229, \
                         gl_236, gl_238, gl_247, gl_249, gl_251, gl_262, gl_264, gl_266, \
                         gl_454, gl_461, gl_463, gl_472, gl_474, gl_476, gl_487, gl_489, \
                         gl_491, gl_544, gl_551, gl_553, gl_562, gl_564, gl_566, gl_577, \
                         gl_579, gl_581 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_107[k] = -f_551 * gl_4[k]
                   - f_552 * gl_11[k]
                   + f_553 * gl_13[k]
                   - f_554 * gl_22[k]
                   + f_192 * gl_24[k]
                   - f_555 * gl_26[k]
                   + f_554 * gl_37[k]
                   - f_556 * gl_39[k]
                   + f_557 * gl_41[k]
                   + f_558 * gl_229[k]
                   + f_559 * gl_236[k]
                   - f_560 * gl_238[k]
                   + f_186 * gl_247[k]
                   - f_200 * gl_249[k]
                   + f_561 * gl_251[k]
                   - f_186 * gl_262[k]
                   + f_188 * gl_264[k]
                   - f_191 * gl_266[k]
                   + f_551 * gl_454[k]
                   + f_552 * gl_461[k]
                   - f_553 * gl_463[k]
                   + f_554 * gl_472[k]
                   - f_192 * gl_474[k]
                   + f_555 * gl_476[k]
                   - f_554 * gl_487[k]
                   + f_556 * gl_489[k]
                   - f_557 * gl_491[k]
                   - f_558 * gl_544[k]
                   - f_559 * gl_551[k]
                   + f_560 * gl_553[k]
                   - f_186 * gl_562[k]
                   + f_200 * gl_564[k]
                   - f_561 * gl_566[k]
                   + f_186 * gl_577[k]
                   - f_188 * gl_579[k]
                   + f_191 * gl_581[k];
    }

#pragma omp simd aligned(gl_1, gl_6, gl_8, gl_15, gl_17, gl_19, gl_28, gl_30, gl_32, gl_34, \
                         gl_226, gl_231, gl_233, gl_240, gl_242, gl_244, gl_253, gl_255, \
                         gl_257, gl_259, gl_451, gl_456, gl_458, gl_465, gl_467, gl_469, \
                         gl_478, gl_480, gl_482, gl_484, gl_541, gl_546, gl_548, gl_555, \
                         gl_557, gl_559, gl_568, gl_570, gl_572, \
                         gl_574 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_108[k] = f_234 * gl_1[k]
                   + f_562 * gl_6[k]
                   - f_235 * gl_8[k]
                   + f_562 * gl_15[k]
                   - f_204 * gl_17[k]
                   + f_236 * gl_19[k]
                   + f_234 * gl_28[k]
                   - f_235 * gl_30[k]
                   + f_236 * gl_32[k]
                   - f_237 * gl_34[k]
                   - f_203 * gl_226[k]
                   - f_563 * gl_231[k]
                   + f_238 * gl_233[k]
                   - f_563 * gl_240[k]
                   + f_210 * gl_242[k]
                   - f_239 * gl_244[k]
                   - f_203 * gl_253[k]
                   + f_238 * gl_255[k]
                   - f_239 * gl_257[k]
                   + f_240 * gl_259[k]
                   - f_234 * gl_451[k]
                   - f_562 * gl_456[k]
                   + f_235 * gl_458[k]
                   - f_562 * gl_465[k]
                   + f_204 * gl_467[k]
                   - f_236 * gl_469[k]
                   - f_234 * gl_478[k]
                   + f_235 * gl_480[k]
                   - f_236 * gl_482[k]
                   + f_237 * gl_484[k]
                   + f_203 * gl_541[k]
                   + f_563 * gl_546[k]
                   - f_238 * gl_548[k]
                   + f_563 * gl_555[k]
                   - f_210 * gl_557[k]
                   + f_239 * gl_559[k]
                   + f_203 * gl_568[k]
                   - f_238 * gl_570[k]
                   + f_239 * gl_572[k]
                   - f_240 * gl_574[k];
    }

#pragma omp simd aligned(gl_4, gl_11, gl_13, gl_22, gl_24, gl_26, gl_37, gl_39, gl_41, gl_43, \
                         gl_229, gl_236, gl_238, gl_247, gl_249, gl_251, gl_262, gl_264, \
                         gl_266, gl_268, gl_454, gl_461, gl_463, gl_472, gl_474, gl_476, \
                         gl_487, gl_489, gl_491, gl_493, gl_544, gl_551, gl_553, gl_562, \
                         gl_564, gl_566, gl_577, gl_579, gl_581, \
                         gl_583 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_109[k] = f_229 * gl_4[k]
                   + f_564 * gl_11[k]
                   - f_565 * gl_13[k]
                   + f_564 * gl_22[k]
                   - f_216 * gl_24[k]
                   + f_566 * gl_26[k]
                   + f_229 * gl_37[k]
                   - f_565 * gl_39[k]
                   + f_566 * gl_41[k]
                   - f_567 * gl_43[k]
                   - f_215 * gl_229[k]
                   - f_568 * gl_236[k]
                   + f_569 * gl_238[k]
                   - f_568 * gl_247[k]
                   + f_222 * gl_249[k]
                   - f_570 * gl_251[k]
                   - f_215 * gl_262[k]
                   + f_569 * gl_264[k]
                   - f_570 * gl_266[k]
                   + f_571 * gl_268[k]
                   - f_229 * gl_454[k]
                   - f_564 * gl_461[k]
                   + f_565 * gl_463[k]
                   - f_564 * gl_472[k]
                   + f_216 * gl_474[k]
                   - f_566 * gl_476[k]
                   - f_229 * gl_487[k]
                   + f_565 * gl_489[k]
                   - f_566 * gl_491[k]
                   + f_567 * gl_493[k]
                   + f_215 * gl_544[k]
                   + f_568 * gl_551[k]
                   - f_569 * gl_553[k]
                   + f_568 * gl_562[k]
                   - f_222 * gl_564[k]
                   + f_570 * gl_566[k]
                   + f_215 * gl_577[k]
                   - f_569 * gl_579[k]
                   + f_570 * gl_581[k]
                   - f_571 * gl_583[k];
    }

#pragma omp simd aligned(gl_0, gl_3, gl_5, gl_10, gl_12, gl_14, gl_21, gl_23, gl_25, gl_27, \
                         gl_36, gl_38, gl_40, gl_42, gl_44, gl_225, gl_228, gl_230, gl_235, \
                         gl_237, gl_239, gl_246, gl_248, gl_250, gl_252, gl_261, gl_263, \
                         gl_265, gl_267, gl_269, gl_450, gl_453, gl_455, gl_460, gl_462, \
                         gl_464, gl_471, gl_473, gl_475, gl_477, gl_486, gl_488, gl_490, \
                         gl_492, gl_494, gl_540, gl_543, gl_545, gl_550, gl_552, gl_554, \
                         gl_561, gl_563, gl_565, gl_567, gl_576, gl_578, gl_580, gl_582, \
                         gl_584 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_110[k] = -f_572 * gl_0[k]
                   - f_573 * gl_3[k]
                   + f_574 * gl_5[k]
                   - f_575 * gl_10[k]
                   + f_565 * gl_12[k]
                   - f_565 * gl_14[k]
                   - f_573 * gl_21[k]
                   + f_565 * gl_23[k]
                   - f_216 * gl_25[k]
                   + f_576 * gl_27[k]
                   - f_572 * gl_36[k]
                   + f_574 * gl_38[k]
                   - f_565 * gl_40[k]
                   + f_576 * gl_42[k]
                   - f_577 * gl_44[k]
                   + f_575 * gl_225[k]
                   + f_214 * gl_228[k]
                   - f_216 * gl_230[k]
                   + f_564 * gl_235[k]
                   - f_569 * gl_237[k]
                   + f_569 * gl_239[k]
                   + f_214 * gl_246[k]
                   - f_569 * gl_248[k]
                   + f_222 * gl_250[k]
                   - f_578 * gl_252[k]
                   + f_575 * gl_261[k]
                   - f_216 * gl_263[k]
                   + f_569 * gl_265[k]
                   - f_578 * gl_267[k]
                   + f_567 * gl_269[k]
                   + f_572 * gl_450[k]
                   + f_573 * gl_453[k]
                   - f_574 * gl_455[k]
                   + f_575 * gl_460[k]
                   - f_565 * gl_462[k]
                   + f_565 * gl_464[k]
                   + f_573 * gl_471[k]
                   - f_565 * gl_473[k]
                   + f_216 * gl_475[k]
                   - f_576 * gl_477[k]
                   + f_572 * gl_486[k]
                   - f_574 * gl_488[k]
                   + f_565 * gl_490[k]
                   - f_576 * gl_492[k]
                   + f_577 * gl_494[k]
                   - f_575 * gl_540[k]
                   - f_214 * gl_543[k]
                   + f_216 * gl_545[k]
                   - f_564 * gl_550[k]
                   + f_569 * gl_552[k]
                   - f_569 * gl_554[k]
                   - f_214 * gl_561[k]
                   + f_569 * gl_563[k]
                   - f_222 * gl_565[k]
                   + f_578 * gl_567[k]
                   - f_575 * gl_576[k]
                   + f_216 * gl_578[k]
                   - f_569 * gl_580[k]
                   + f_578 * gl_582[k]
                   - f_567 * gl_584[k];
    }

#pragma omp simd aligned(gl_2, gl_7, gl_9, gl_16, gl_18, gl_20, gl_29, gl_31, gl_33, gl_35, \
                         gl_227, gl_232, gl_234, gl_241, gl_243, gl_245, gl_254, gl_256, \
                         gl_258, gl_260, gl_452, gl_457, gl_459, gl_466, gl_468, gl_470, \
                         gl_479, gl_481, gl_483, gl_485, gl_542, gl_547, gl_549, gl_556, \
                         gl_558, gl_560, gl_569, gl_571, gl_573, \
                         gl_575 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_111[k] = f_229 * gl_2[k]
                   + f_564 * gl_7[k]
                   - f_565 * gl_9[k]
                   + f_564 * gl_16[k]
                   - f_216 * gl_18[k]
                   + f_566 * gl_20[k]
                   + f_229 * gl_29[k]
                   - f_565 * gl_31[k]
                   + f_566 * gl_33[k]
                   - f_567 * gl_35[k]
                   - f_215 * gl_227[k]
                   - f_568 * gl_232[k]
                   + f_569 * gl_234[k]
                   - f_568 * gl_241[k]
                   + f_222 * gl_243[k]
                   - f_570 * gl_245[k]
                   - f_215 * gl_254[k]
                   + f_569 * gl_256[k]
                   - f_570 * gl_258[k]
                   + f_571 * gl_260[k]
                   - f_229 * gl_452[k]
                   - f_564 * gl_457[k]
                   + f_565 * gl_459[k]
                   - f_564 * gl_466[k]
                   + f_216 * gl_468[k]
                   - f_566 * gl_470[k]
                   - f_229 * gl_479[k]
                   + f_565 * gl_481[k]
                   - f_566 * gl_483[k]
                   + f_567 * gl_485[k]
                   + f_215 * gl_542[k]
                   + f_568 * gl_547[k]
                   - f_569 * gl_549[k]
                   + f_568 * gl_556[k]
                   - f_222 * gl_558[k]
                   + f_570 * gl_560[k]
                   + f_215 * gl_569[k]
                   - f_569 * gl_571[k]
                   + f_570 * gl_573[k]
                   - f_571 * gl_575[k];
    }

#pragma omp simd aligned(gl_0, gl_3, gl_5, gl_12, gl_14, gl_21, gl_23, gl_27, gl_36, gl_38, \
                         gl_40, gl_42, gl_225, gl_228, gl_230, gl_237, gl_239, gl_246, gl_248, \
                         gl_252, gl_261, gl_263, gl_265, gl_267, gl_450, gl_453, gl_455, \
                         gl_462, gl_464, gl_471, gl_473, gl_477, gl_486, gl_488, gl_490, \
                         gl_492, gl_540, gl_543, gl_545, gl_552, gl_554, gl_561, gl_563, \
                         gl_567, gl_576, gl_578, gl_580, gl_582 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_112[k] = f_579 * gl_0[k]
                   + f_234 * gl_3[k]
                   - f_580 * gl_5[k]
                   - f_580 * gl_12[k]
                   + f_581 * gl_14[k]
                   - f_234 * gl_21[k]
                   + f_580 * gl_23[k]
                   - f_582 * gl_27[k]
                   - f_579 * gl_36[k]
                   + f_580 * gl_38[k]
                   - f_581 * gl_40[k]
                   + f_582 * gl_42[k]
                   - f_562 * gl_225[k]
                   - f_203 * gl_228[k]
                   + f_583 * gl_230[k]
                   + f_583 * gl_237[k]
                   - f_584 * gl_239[k]
                   + f_203 * gl_246[k]
                   - f_583 * gl_248[k]
                   + f_585 * gl_252[k]
                   + f_562 * gl_261[k]
                   - f_583 * gl_263[k]
                   + f_584 * gl_265[k]
                   - f_585 * gl_267[k]
                   - f_579 * gl_450[k]
                   - f_234 * gl_453[k]
                   + f_580 * gl_455[k]
                   + f_580 * gl_462[k]
                   - f_581 * gl_464[k]
                   + f_234 * gl_471[k]
                   - f_580 * gl_473[k]
                   + f_582 * gl_477[k]
                   + f_579 * gl_486[k]
                   - f_580 * gl_488[k]
                   + f_581 * gl_490[k]
                   - f_582 * gl_492[k]
                   + f_562 * gl_540[k]
                   + f_203 * gl_543[k]
                   - f_583 * gl_545[k]
                   - f_583 * gl_552[k]
                   + f_584 * gl_554[k]
                   - f_203 * gl_561[k]
                   + f_583 * gl_563[k]
                   - f_585 * gl_567[k]
                   - f_562 * gl_576[k]
                   + f_583 * gl_578[k]
                   - f_584 * gl_580[k]
                   + f_585 * gl_582[k];
    }

#pragma omp simd aligned(gl_2, gl_7, gl_9, gl_16, gl_18, gl_20, gl_29, gl_31, gl_33, gl_227, \
                         gl_232, gl_234, gl_241, gl_243, gl_245, gl_254, gl_256, gl_258, \
                         gl_452, gl_457, gl_459, gl_466, gl_468, gl_470, gl_479, gl_481, \
                         gl_483, gl_542, gl_547, gl_549, gl_556, gl_558, gl_560, gl_569, \
                         gl_571, gl_573 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_113[k] = -f_554 * gl_2[k]
                   + f_554 * gl_7[k]
                   + f_556 * gl_9[k]
                   + f_552 * gl_16[k]
                   - f_192 * gl_18[k]
                   - f_557 * gl_20[k]
                   + f_551 * gl_29[k]
                   - f_553 * gl_31[k]
                   + f_555 * gl_33[k]
                   + f_186 * gl_227[k]
                   - f_186 * gl_232[k]
                   - f_188 * gl_234[k]
                   - f_559 * gl_241[k]
                   + f_200 * gl_243[k]
                   + f_191 * gl_245[k]
                   - f_558 * gl_254[k]
                   + f_560 * gl_256[k]
                   - f_561 * gl_258[k]
                   + f_554 * gl_452[k]
                   - f_554 * gl_457[k]
                   - f_556 * gl_459[k]
                   - f_552 * gl_466[k]
                   + f_192 * gl_468[k]
                   + f_557 * gl_470[k]
                   - f_551 * gl_479[k]
                   + f_553 * gl_481[k]
                   - f_555 * gl_483[k]
                   - f_186 * gl_542[k]
                   + f_186 * gl_547[k]
                   + f_188 * gl_549[k]
                   + f_559 * gl_556[k]
                   - f_200 * gl_558[k]
                   - f_191 * gl_560[k]
                   + f_558 * gl_569[k]
                   - f_560 * gl_571[k]
                   + f_561 * gl_573[k];
    }

#pragma omp simd aligned(gl_0, gl_3, gl_5, gl_10, gl_12, gl_14, gl_21, gl_23, gl_25, gl_36, \
                         gl_38, gl_40, gl_225, gl_228, gl_230, gl_235, gl_237, gl_239, gl_246, \
                         gl_248, gl_250, gl_261, gl_263, gl_265, gl_450, gl_453, gl_455, \
                         gl_460, gl_462, gl_464, gl_471, gl_473, gl_475, gl_486, gl_488, \
                         gl_490, gl_540, gl_543, gl_545, gl_550, gl_552, gl_554, gl_561, \
                         gl_563, gl_565, gl_576, gl_578, gl_580 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_114[k] = -f_586 * gl_0[k]
                   + f_545 * gl_3[k]
                   + f_548 * gl_5[k]
                   + f_587 * gl_10[k]
                   - f_248 * gl_12[k]
                   - f_588 * gl_14[k]
                   + f_545 * gl_21[k]
                   - f_248 * gl_23[k]
                   + f_243 * gl_25[k]
                   - f_586 * gl_36[k]
                   + f_548 * gl_38[k]
                   - f_588 * gl_40[k]
                   + f_589 * gl_225[k]
                   - f_548 * gl_228[k]
                   - f_590 * gl_230[k]
                   - f_591 * gl_235[k]
                   + f_592 * gl_237[k]
                   + f_243 * gl_239[k]
                   - f_548 * gl_246[k]
                   + f_592 * gl_248[k]
                   - f_249 * gl_250[k]
                   + f_589 * gl_261[k]
                   - f_590 * gl_263[k]
                   + f_243 * gl_265[k]
                   + f_586 * gl_450[k]
                   - f_545 * gl_453[k]
                   - f_548 * gl_455[k]
                   - f_587 * gl_460[k]
                   + f_248 * gl_462[k]
                   + f_588 * gl_464[k]
                   - f_545 * gl_471[k]
                   + f_248 * gl_473[k]
                   - f_243 * gl_475[k]
                   + f_586 * gl_486[k]
                   - f_548 * gl_488[k]
                   + f_588 * gl_490[k]
                   - f_589 * gl_540[k]
                   + f_548 * gl_543[k]
                   + f_590 * gl_545[k]
                   + f_591 * gl_550[k]
                   - f_592 * gl_552[k]
                   - f_243 * gl_554[k]
                   + f_548 * gl_561[k]
                   - f_592 * gl_563[k]
                   + f_249 * gl_565[k]
                   - f_589 * gl_576[k]
                   + f_590 * gl_578[k]
                   - f_243 * gl_580[k];
    }

#pragma omp simd aligned(gl_2, gl_7, gl_9, gl_16, gl_18, gl_29, gl_31, gl_227, gl_232, gl_234, \
                         gl_241, gl_243, gl_254, gl_256, gl_452, gl_457, gl_459, gl_466, \
                         gl_468, gl_479, gl_481, gl_542, gl_547, gl_549, gl_556, gl_558, \
                         gl_569, gl_571 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_115[k] = f_538 * gl_2[k]
                   - f_537 * gl_7[k]
                   - f_539 * gl_9[k]
                   - f_535 * gl_16[k]
                   + f_169 * gl_18[k]
                   + f_535 * gl_29[k]
                   - f_536 * gl_31[k]
                   - f_543 * gl_227[k]
                   + f_542 * gl_232[k]
                   + f_544 * gl_234[k]
                   + f_540 * gl_241[k]
                   - f_175 * gl_243[k]
                   - f_540 * gl_254[k]
                   + f_541 * gl_256[k]
                   - f_538 * gl_452[k]
                   + f_537 * gl_457[k]
                   + f_539 * gl_459[k]
                   + f_535 * gl_466[k]
                   - f_169 * gl_468[k]
                   - f_535 * gl_479[k]
                   + f_536 * gl_481[k]
                   + f_543 * gl_542[k]
                   - f_542 * gl_547[k]
                   - f_544 * gl_549[k]
                   - f_540 * gl_556[k]
                   + f_175 * gl_558[k]
                   + f_540 * gl_569[k]
                   - f_541 * gl_571[k];
    }

#pragma omp simd aligned(gl_0, gl_3, gl_5, gl_12, gl_21, gl_23, gl_36, gl_38, gl_225, gl_228, \
                         gl_230, gl_237, gl_246, gl_248, gl_261, gl_263, gl_450, gl_453, \
                         gl_455, gl_462, gl_471, gl_473, gl_486, gl_488, gl_540, gl_543, \
                         gl_545, gl_552, gl_561, gl_563, gl_576, \
                         gl_578 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_116[k] = f_593 * gl_0[k]
                   - f_529 * gl_3[k]
                   - f_529 * gl_5[k]
                   + f_594 * gl_12[k]
                   + f_529 * gl_21[k]
                   - f_594 * gl_23[k]
                   - f_593 * gl_36[k]
                   + f_529 * gl_38[k]
                   - f_528 * gl_225[k]
                   + f_530 * gl_228[k]
                   + f_530 * gl_230[k]
                   - f_595 * gl_237[k]
                   - f_530 * gl_246[k]
                   + f_595 * gl_248[k]
                   + f_528 * gl_261[k]
                   - f_530 * gl_263[k]
                   - f_593 * gl_450[k]
                   + f_529 * gl_453[k]
                   + f_529 * gl_455[k]
                   - f_594 * gl_462[k]
                   - f_529 * gl_471[k]
                   + f_594 * gl_473[k]
                   + f_593 * gl_486[k]
                   - f_529 * gl_488[k]
                   + f_528 * gl_540[k]
                   - f_530 * gl_543[k]
                   - f_530 * gl_545[k]
                   + f_595 * gl_552[k]
                   + f_530 * gl_561[k]
                   - f_595 * gl_563[k]
                   - f_528 * gl_576[k]
                   + f_530 * gl_578[k];
    }

#pragma omp simd aligned(gl_2, gl_7, gl_16, gl_29, gl_227, gl_232, gl_241, gl_254, gl_452, \
                         gl_457, gl_466, gl_479, gl_542, gl_547, gl_556, \
                         gl_569 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_117[k] = -f_525 * gl_2[k]
                   + f_524 * gl_7[k]
                   - f_255 * gl_16[k]
                   + f_523 * gl_29[k]
                   + f_527 * gl_227[k]
                   - f_526 * gl_232[k]
                   + f_257 * gl_241[k]
                   - f_155 * gl_254[k]
                   + f_525 * gl_452[k]
                   - f_524 * gl_457[k]
                   + f_255 * gl_466[k]
                   - f_523 * gl_479[k]
                   - f_527 * gl_542[k]
                   + f_526 * gl_547[k]
                   - f_257 * gl_556[k]
                   + f_155 * gl_569[k];
    }

#pragma omp simd aligned(gl_0, gl_3, gl_10, gl_21, gl_36, gl_225, gl_228, gl_235, gl_246, \
                         gl_261, gl_450, gl_453, gl_460, gl_471, gl_486, gl_540, gl_543, \
                         gl_550, gl_561, gl_576 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_118[k] = -f_596 * gl_0[k]
                   + f_523 * gl_3[k]
                   - f_597 * gl_10[k]
                   + f_523 * gl_21[k]
                   - f_596 * gl_36[k]
                   + f_598 * gl_225[k]
                   - f_155 * gl_228[k]
                   + f_599 * gl_235[k]
                   - f_155 * gl_246[k]
                   + f_598 * gl_261[k]
                   + f_596 * gl_450[k]
                   - f_523 * gl_453[k]
                   + f_597 * gl_460[k]
                   - f_523 * gl_471[k]
                   + f_596 * gl_486[k]
                   - f_598 * gl_540[k]
                   + f_155 * gl_543[k]
                   - f_599 * gl_550[k]
                   + f_155 * gl_561[k]
                   - f_598 * gl_576[k];
    }

#pragma omp simd aligned(gl_91, gl_96, gl_105, gl_118, gl_316, gl_321, gl_330, \
                         gl_343 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_119[k] = f_61 * gl_91[k]
                   - f_62 * gl_96[k]
                   + f_62 * gl_105[k]
                   - f_61 * gl_118[k]
                   - f_59 * gl_316[k]
                   + f_60 * gl_321[k]
                   - f_60 * gl_330[k]
                   + f_59 * gl_343[k];
    }

#pragma omp simd aligned(gl_94, gl_101, gl_112, gl_127, gl_319, gl_326, gl_337, \
                         gl_352 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_120[k] = f_67 * gl_94[k]
                   - f_68 * gl_101[k]
                   + f_63 * gl_112[k]
                   - f_69 * gl_127[k]
                   - f_63 * gl_319[k]
                   + f_64 * gl_326[k]
                   - f_65 * gl_337[k]
                   + f_66 * gl_352[k];
    }

#pragma omp simd aligned(gl_91, gl_96, gl_98, gl_105, gl_107, gl_118, gl_120, gl_316, gl_321, \
                         gl_323, gl_330, gl_332, gl_343, gl_345 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_121[k] = -f_74 * gl_91[k]
                   + f_75 * gl_96[k]
                   + f_76 * gl_98[k]
                   + f_75 * gl_105[k]
                   - f_77 * gl_107[k]
                   - f_74 * gl_118[k]
                   + f_76 * gl_120[k]
                   + f_70 * gl_316[k]
                   - f_71 * gl_321[k]
                   - f_72 * gl_323[k]
                   - f_71 * gl_330[k]
                   + f_73 * gl_332[k]
                   + f_70 * gl_343[k]
                   - f_72 * gl_345[k];
    }

#pragma omp simd aligned(gl_94, gl_101, gl_103, gl_112, gl_114, gl_127, gl_129, gl_319, \
                         gl_326, gl_328, gl_337, gl_339, gl_352, \
                         gl_354 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_122[k] = -f_84 * gl_94[k]
                   + f_84 * gl_101[k]
                   + f_85 * gl_103[k]
                   + f_86 * gl_112[k]
                   - f_87 * gl_114[k]
                   - f_88 * gl_127[k]
                   + f_89 * gl_129[k]
                   + f_78 * gl_319[k]
                   - f_78 * gl_326[k]
                   - f_79 * gl_328[k]
                   - f_80 * gl_337[k]
                   + f_81 * gl_339[k]
                   + f_82 * gl_352[k]
                   - f_83 * gl_354[k];
    }

#pragma omp simd aligned(gl_91, gl_96, gl_98, gl_105, gl_109, gl_118, gl_120, gl_122, gl_316, \
                         gl_321, gl_323, gl_330, gl_334, gl_343, gl_345, \
                         gl_347 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_123[k] = f_93 * gl_91[k]
                   + f_93 * gl_96[k]
                   - f_94 * gl_98[k]
                   - f_93 * gl_105[k]
                   + f_95 * gl_109[k]
                   - f_93 * gl_118[k]
                   + f_94 * gl_120[k]
                   - f_95 * gl_122[k]
                   - f_90 * gl_316[k]
                   - f_90 * gl_321[k]
                   + f_91 * gl_323[k]
                   + f_90 * gl_330[k]
                   - f_92 * gl_334[k]
                   + f_90 * gl_343[k]
                   - f_91 * gl_345[k]
                   + f_92 * gl_347[k];
    }

#pragma omp simd aligned(gl_94, gl_101, gl_103, gl_112, gl_114, gl_116, gl_127, gl_129, \
                         gl_131, gl_319, gl_326, gl_328, gl_337, gl_339, gl_341, gl_352, \
                         gl_354, gl_356 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_124[k] = f_99 * gl_94[k]
                   + f_104 * gl_101[k]
                   - f_102 * gl_103[k]
                   + f_105 * gl_112[k]
                   - f_106 * gl_114[k]
                   + f_103 * gl_116[k]
                   - f_105 * gl_127[k]
                   + f_107 * gl_129[k]
                   - f_108 * gl_131[k]
                   - f_96 * gl_319[k]
                   - f_97 * gl_326[k]
                   + f_98 * gl_328[k]
                   - f_99 * gl_337[k]
                   + f_100 * gl_339[k]
                   - f_101 * gl_341[k]
                   + f_99 * gl_352[k]
                   - f_102 * gl_354[k]
                   + f_103 * gl_356[k];
    }

#pragma omp simd aligned(gl_91, gl_96, gl_98, gl_105, gl_107, gl_109, gl_118, gl_120, gl_122, \
                         gl_124, gl_316, gl_321, gl_323, gl_330, gl_332, gl_334, gl_343, \
                         gl_345, gl_347, gl_349 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_125[k] = -1.640625 * gl_91[k]
                   - 4.921875 * gl_96[k]
                   + 49.21875 * gl_98[k]
                   - 4.921875 * gl_105[k]
                   + 98.4375 * gl_107[k]
                   - 131.25 * gl_109[k]
                   - 1.640625 * gl_118[k]
                   + 49.21875 * gl_120[k]
                   - 131.25 * gl_122[k]
                   + 52.5 * gl_124[k]
                   + 4.921875 * gl_316[k]
                   + 14.765625 * gl_321[k]
                   - 147.65625 * gl_323[k]
                   + 14.765625 * gl_330[k]
                   - 295.3125 * gl_332[k]
                   + 393.75 * gl_334[k]
                   + 4.921875 * gl_343[k]
                   - 147.65625 * gl_345[k]
                   + 393.75 * gl_347[k]
                   - 157.5 * gl_349[k];
    }

#pragma omp simd aligned(gl_94, gl_101, gl_103, gl_112, gl_114, gl_116, gl_127, gl_129, \
                         gl_131, gl_133, gl_319, gl_326, gl_328, gl_337, gl_339, gl_341, \
                         gl_352, gl_354, gl_356, gl_358 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_126[k] = -f_115 * gl_94[k]
                   - f_109 * gl_101[k]
                   + f_116 * gl_103[k]
                   - f_109 * gl_112[k]
                   + f_117 * gl_114[k]
                   - f_118 * gl_116[k]
                   - f_115 * gl_127[k]
                   + f_116 * gl_129[k]
                   - f_118 * gl_131[k]
                   + f_119 * gl_133[k]
                   + f_109 * gl_319[k]
                   + f_110 * gl_326[k]
                   - f_111 * gl_328[k]
                   + f_110 * gl_337[k]
                   - f_112 * gl_339[k]
                   + f_113 * gl_341[k]
                   + f_109 * gl_352[k]
                   - f_111 * gl_354[k]
                   + f_113 * gl_356[k]
                   - f_114 * gl_358[k];
    }

#pragma omp simd aligned(gl_90, gl_93, gl_95, gl_100, gl_102, gl_104, gl_111, gl_113, gl_115, \
                         gl_117, gl_126, gl_128, gl_130, gl_132, gl_134, gl_315, gl_318, \
                         gl_320, gl_325, gl_327, gl_329, gl_336, gl_338, gl_340, gl_342, \
                         gl_351, gl_353, gl_355, gl_357, gl_359 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_127[k] = f_124 * gl_90[k]
                   + f_125 * gl_93[k]
                   - f_126 * gl_95[k]
                   + f_127 * gl_100[k]
                   - f_116 * gl_102[k]
                   + f_116 * gl_104[k]
                   + f_125 * gl_111[k]
                   - f_116 * gl_113[k]
                   + f_117 * gl_115[k]
                   - f_128 * gl_117[k]
                   + f_124 * gl_126[k]
                   - f_126 * gl_128[k]
                   + f_116 * gl_130[k]
                   - f_128 * gl_132[k]
                   + f_129 * gl_134[k]
                   - f_120 * gl_315[k]
                   - f_115 * gl_318[k]
                   + f_116 * gl_320[k]
                   - f_121 * gl_325[k]
                   + f_111 * gl_327[k]
                   - f_111 * gl_329[k]
                   - f_115 * gl_336[k]
                   + f_111 * gl_338[k]
                   - f_112 * gl_340[k]
                   + f_122 * gl_342[k]
                   - f_120 * gl_351[k]
                   + f_116 * gl_353[k]
                   - f_111 * gl_355[k]
                   + f_122 * gl_357[k]
                   - f_123 * gl_359[k];
    }

#pragma omp simd aligned(gl_92, gl_97, gl_99, gl_106, gl_108, gl_110, gl_119, gl_121, gl_123, \
                         gl_125, gl_317, gl_322, gl_324, gl_331, gl_333, gl_335, gl_344, \
                         gl_346, gl_348, gl_350 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_128[k] = -f_115 * gl_92[k]
                   - f_109 * gl_97[k]
                   + f_116 * gl_99[k]
                   - f_109 * gl_106[k]
                   + f_117 * gl_108[k]
                   - f_118 * gl_110[k]
                   - f_115 * gl_119[k]
                   + f_116 * gl_121[k]
                   - f_118 * gl_123[k]
                   + f_119 * gl_125[k]
                   + f_109 * gl_317[k]
                   + f_110 * gl_322[k]
                   - f_111 * gl_324[k]
                   + f_110 * gl_331[k]
                   - f_112 * gl_333[k]
                   + f_113 * gl_335[k]
                   + f_109 * gl_344[k]
                   - f_111 * gl_346[k]
                   + f_113 * gl_348[k]
                   - f_114 * gl_350[k];
    }

#pragma omp simd aligned(gl_90, gl_93, gl_95, gl_102, gl_104, gl_111, gl_113, gl_117, gl_126, \
                         gl_128, gl_130, gl_132, gl_315, gl_318, gl_320, gl_327, gl_329, \
                         gl_336, gl_338, gl_342, gl_351, gl_353, gl_355, \
                         gl_357 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_129[k] = -0.8203125 * gl_90[k]
                   - 1.640625 * gl_93[k]
                   + 24.609375 * gl_95[k]
                   + 24.609375 * gl_102[k]
                   - 65.625 * gl_104[k]
                   + 1.640625 * gl_111[k]
                   - 24.609375 * gl_113[k]
                   + 26.25 * gl_117[k]
                   + 0.8203125 * gl_126[k]
                   - 24.609375 * gl_128[k]
                   + 65.625 * gl_130[k]
                   - 26.25 * gl_132[k]
                   + 2.4609375 * gl_315[k]
                   + 4.921875 * gl_318[k]
                   - 73.828125 * gl_320[k]
                   - 73.828125 * gl_327[k]
                   + 196.875 * gl_329[k]
                   - 4.921875 * gl_336[k]
                   + 73.828125 * gl_338[k]
                   - 78.75 * gl_342[k]
                   - 2.4609375 * gl_351[k]
                   + 73.828125 * gl_353[k]
                   - 196.875 * gl_355[k]
                   + 78.75 * gl_357[k];
    }

#pragma omp simd aligned(gl_92, gl_97, gl_99, gl_106, gl_108, gl_110, gl_119, gl_121, gl_123, \
                         gl_317, gl_322, gl_324, gl_331, gl_333, gl_335, gl_344, gl_346, \
                         gl_348 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_130[k] = f_105 * gl_92[k]
                   - f_105 * gl_97[k]
                   - f_107 * gl_99[k]
                   - f_104 * gl_106[k]
                   + f_106 * gl_108[k]
                   + f_108 * gl_110[k]
                   - f_99 * gl_119[k]
                   + f_102 * gl_121[k]
                   - f_103 * gl_123[k]
                   - f_99 * gl_317[k]
                   + f_99 * gl_322[k]
                   + f_102 * gl_324[k]
                   + f_97 * gl_331[k]
                   - f_100 * gl_333[k]
                   - f_103 * gl_335[k]
                   + f_96 * gl_344[k]
                   - f_98 * gl_346[k]
                   + f_101 * gl_348[k];
    }

#pragma omp simd aligned(gl_90, gl_93, gl_95, gl_100, gl_102, gl_104, gl_111, gl_113, gl_115, \
                         gl_126, gl_128, gl_130, gl_315, gl_318, gl_320, gl_325, gl_327, \
                         gl_329, gl_336, gl_338, gl_340, gl_351, gl_353, \
                         gl_355 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_131[k] = f_136 * gl_90[k]
                   - f_93 * gl_93[k]
                   - f_137 * gl_95[k]
                   - f_138 * gl_100[k]
                   + f_134 * gl_102[k]
                   + f_139 * gl_104[k]
                   - f_93 * gl_111[k]
                   + f_134 * gl_113[k]
                   - f_140 * gl_115[k]
                   + f_136 * gl_126[k]
                   - f_137 * gl_128[k]
                   + f_139 * gl_130[k]
                   - f_130 * gl_315[k]
                   + f_90 * gl_318[k]
                   + f_131 * gl_320[k]
                   + f_132 * gl_325[k]
                   - f_133 * gl_327[k]
                   - f_134 * gl_329[k]
                   + f_90 * gl_336[k]
                   - f_133 * gl_338[k]
                   + f_135 * gl_340[k]
                   - f_130 * gl_351[k]
                   + f_131 * gl_353[k]
                   - f_134 * gl_355[k];
    }

#pragma omp simd aligned(gl_92, gl_97, gl_99, gl_106, gl_108, gl_119, gl_121, gl_317, gl_322, \
                         gl_324, gl_331, gl_333, gl_344, gl_346 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_132[k] = -f_88 * gl_92[k]
                   + f_86 * gl_97[k]
                   + f_89 * gl_99[k]
                   + f_84 * gl_106[k]
                   - f_87 * gl_108[k]
                   - f_84 * gl_119[k]
                   + f_85 * gl_121[k]
                   + f_82 * gl_317[k]
                   - f_80 * gl_322[k]
                   - f_83 * gl_324[k]
                   - f_78 * gl_331[k]
                   + f_81 * gl_333[k]
                   + f_78 * gl_344[k]
                   - f_79 * gl_346[k];
    }

#pragma omp simd aligned(gl_90, gl_93, gl_95, gl_102, gl_111, gl_113, gl_126, gl_128, gl_315, \
                         gl_318, gl_320, gl_327, gl_336, gl_338, gl_351, \
                         gl_353 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_133[k] = -f_143 * gl_90[k]
                   + f_75 * gl_93[k]
                   + f_75 * gl_95[k]
                   - f_144 * gl_102[k]
                   - f_75 * gl_111[k]
                   + f_144 * gl_113[k]
                   + f_143 * gl_126[k]
                   - f_75 * gl_128[k]
                   + f_141 * gl_315[k]
                   - f_71 * gl_318[k]
                   - f_71 * gl_320[k]
                   + f_142 * gl_327[k]
                   + f_71 * gl_336[k]
                   - f_142 * gl_338[k]
                   - f_141 * gl_351[k]
                   + f_71 * gl_353[k];
    }

#pragma omp simd aligned(gl_92, gl_97, gl_106, gl_119, gl_317, gl_322, gl_331, \
                         gl_344 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_134[k] = f_69 * gl_92[k]
                   - f_63 * gl_97[k]
                   + f_68 * gl_106[k]
                   - f_67 * gl_119[k]
                   - f_66 * gl_317[k]
                   + f_65 * gl_322[k]
                   - f_64 * gl_331[k]
                   + f_63 * gl_344[k];
    }

#pragma omp simd aligned(gl_90, gl_93, gl_100, gl_111, gl_126, gl_315, gl_318, gl_325, gl_336, \
                         gl_351 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_135[k] = f_147 * gl_90[k]
                   - f_67 * gl_93[k]
                   + f_148 * gl_100[k]
                   - f_67 * gl_111[k]
                   + f_147 * gl_126[k]
                   - f_145 * gl_315[k]
                   + f_63 * gl_318[k]
                   - f_146 * gl_325[k]
                   + f_63 * gl_336[k]
                   - f_145 * gl_351[k];
    }

#pragma omp simd aligned(gl_1, gl_6, gl_15, gl_28, gl_136, gl_141, gl_150, gl_163, gl_451, \
                         gl_456, gl_465, gl_478 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_136[k] = f_600 * gl_1[k]
                   - f_601 * gl_6[k]
                   + f_601 * gl_15[k]
                   - f_600 * gl_28[k]
                   - f_408 * gl_136[k]
                   + f_4 * gl_141[k]
                   - f_4 * gl_150[k]
                   + f_408 * gl_163[k]
                   + f_600 * gl_451[k]
                   - f_601 * gl_456[k]
                   + f_601 * gl_465[k]
                   - f_600 * gl_478[k];
    }

#pragma omp simd aligned(gl_4, gl_11, gl_22, gl_37, gl_139, gl_146, gl_157, gl_172, gl_454, \
                         gl_461, gl_472, gl_487 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_137[k] = f_602 * gl_4[k]
                   - f_603 * gl_11[k]
                   + f_604 * gl_22[k]
                   - f_57 * gl_37[k]
                   - f_605 * gl_139[k]
                   + f_606 * gl_146[k]
                   - f_607 * gl_157[k]
                   + f_413 * gl_172[k]
                   + f_602 * gl_454[k]
                   - f_603 * gl_461[k]
                   + f_604 * gl_472[k]
                   - f_57 * gl_487[k];
    }

#pragma omp simd aligned(gl_1, gl_6, gl_8, gl_15, gl_17, gl_28, gl_30, gl_136, gl_141, gl_143, \
                         gl_150, gl_152, gl_163, gl_165, gl_451, gl_456, gl_458, gl_465, \
                         gl_467, gl_478, gl_480 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_138[k] = -f_608 * gl_1[k]
                   + f_609 * gl_6[k]
                   + f_610 * gl_8[k]
                   + f_609 * gl_15[k]
                   - f_611 * gl_17[k]
                   - f_608 * gl_28[k]
                   + f_610 * gl_30[k]
                   + f_612 * gl_136[k]
                   - f_610 * gl_141[k]
                   - f_613 * gl_143[k]
                   - f_610 * gl_150[k]
                   + f_614 * gl_152[k]
                   + f_612 * gl_163[k]
                   - f_613 * gl_165[k]
                   - f_608 * gl_451[k]
                   + f_609 * gl_456[k]
                   + f_610 * gl_458[k]
                   + f_609 * gl_465[k]
                   - f_611 * gl_467[k]
                   - f_608 * gl_478[k]
                   + f_610 * gl_480[k];
    }

#pragma omp simd aligned(gl_4, gl_11, gl_13, gl_22, gl_24, gl_37, gl_39, gl_139, gl_146, \
                         gl_148, gl_157, gl_159, gl_172, gl_174, gl_454, gl_461, gl_463, \
                         gl_472, gl_474, gl_487, gl_489 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_139[k] = -f_615 * gl_4[k]
                   + f_615 * gl_11[k]
                   + f_10 * gl_13[k]
                   + f_380 * gl_22[k]
                   - f_389 * gl_24[k]
                   - f_616 * gl_37[k]
                   + f_14 * gl_39[k]
                   + f_382 * gl_139[k]
                   - f_382 * gl_146[k]
                   - f_385 * gl_148[k]
                   - f_617 * gl_157[k]
                   + f_618 * gl_159[k]
                   + f_372 * gl_172[k]
                   - f_384 * gl_174[k]
                   - f_615 * gl_454[k]
                   + f_615 * gl_461[k]
                   + f_10 * gl_463[k]
                   + f_380 * gl_472[k]
                   - f_389 * gl_474[k]
                   - f_616 * gl_487[k]
                   + f_14 * gl_489[k];
    }

#pragma omp simd aligned(gl_1, gl_6, gl_8, gl_15, gl_19, gl_28, gl_30, gl_32, gl_136, gl_141, \
                         gl_143, gl_150, gl_154, gl_163, gl_165, gl_167, gl_451, gl_456, \
                         gl_458, gl_465, gl_469, gl_478, gl_480, \
                         gl_482 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_140[k] = f_49 * gl_1[k]
                   + f_49 * gl_6[k]
                   - f_50 * gl_8[k]
                   - f_49 * gl_15[k]
                   + f_53 * gl_19[k]
                   - f_49 * gl_28[k]
                   + f_50 * gl_30[k]
                   - f_53 * gl_32[k]
                   - f_619 * gl_136[k]
                   - f_619 * gl_141[k]
                   + f_620 * gl_143[k]
                   + f_619 * gl_150[k]
                   - f_54 * gl_154[k]
                   + f_619 * gl_163[k]
                   - f_620 * gl_165[k]
                   + f_54 * gl_167[k]
                   + f_49 * gl_451[k]
                   + f_49 * gl_456[k]
                   - f_50 * gl_458[k]
                   - f_49 * gl_465[k]
                   + f_53 * gl_469[k]
                   - f_49 * gl_478[k]
                   + f_50 * gl_480[k]
                   - f_53 * gl_482[k];
    }

#pragma omp simd aligned(gl_4, gl_11, gl_13, gl_22, gl_24, gl_26, gl_37, gl_39, gl_41, gl_139, \
                         gl_146, gl_148, gl_157, gl_159, gl_161, gl_172, gl_174, gl_176, \
                         gl_454, gl_461, gl_463, gl_472, gl_474, gl_476, gl_487, gl_489, \
                         gl_491 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_141[k] = f_621 * gl_4[k]
                   + f_622 * gl_11[k]
                   - f_20 * gl_13[k]
                   + f_623 * gl_22[k]
                   - f_624 * gl_24[k]
                   + f_625 * gl_26[k]
                   - f_623 * gl_37[k]
                   + f_626 * gl_39[k]
                   - f_627 * gl_41[k]
                   - f_628 * gl_139[k]
                   - f_629 * gl_146[k]
                   + f_630 * gl_148[k]
                   - f_631 * gl_157[k]
                   + f_21 * gl_159[k]
                   - f_632 * gl_161[k]
                   + f_631 * gl_172[k]
                   - f_633 * gl_174[k]
                   + f_634 * gl_176[k]
                   + f_621 * gl_454[k]
                   + f_622 * gl_461[k]
                   - f_20 * gl_463[k]
                   + f_623 * gl_472[k]
                   - f_624 * gl_474[k]
                   + f_625 * gl_476[k]
                   - f_623 * gl_487[k]
                   + f_626 * gl_489[k]
                   - f_627 * gl_491[k];
    }

#pragma omp simd aligned(gl_1, gl_6, gl_8, gl_15, gl_17, gl_19, gl_28, gl_30, gl_32, gl_34, \
                         gl_136, gl_141, gl_143, gl_150, gl_152, gl_154, gl_163, gl_165, \
                         gl_167, gl_169, gl_451, gl_456, gl_458, gl_465, gl_467, gl_469, \
                         gl_478, gl_480, gl_482, gl_484 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_142[k] = -f_635 * gl_1[k]
                   - f_636 * gl_6[k]
                   + f_637 * gl_8[k]
                   - f_636 * gl_15[k]
                   + f_46 * gl_17[k]
                   - f_638 * gl_19[k]
                   - f_635 * gl_28[k]
                   + f_637 * gl_30[k]
                   - f_638 * gl_32[k]
                   + f_639 * gl_34[k]
                   + f_640 * gl_136[k]
                   + f_641 * gl_141[k]
                   - f_642 * gl_143[k]
                   + f_641 * gl_150[k]
                   - f_643 * gl_152[k]
                   + f_644 * gl_154[k]
                   + f_640 * gl_163[k]
                   - f_642 * gl_165[k]
                   + f_644 * gl_167[k]
                   - f_645 * gl_169[k]
                   - f_635 * gl_451[k]
                   - f_636 * gl_456[k]
                   + f_637 * gl_458[k]
                   - f_636 * gl_465[k]
                   + f_46 * gl_467[k]
                   - f_638 * gl_469[k]
                   - f_635 * gl_478[k]
                   + f_637 * gl_480[k]
                   - f_638 * gl_482[k]
                   + f_639 * gl_484[k];
    }

#pragma omp simd aligned(gl_4, gl_11, gl_13, gl_22, gl_24, gl_26, gl_37, gl_39, gl_41, gl_43, \
                         gl_139, gl_146, gl_148, gl_157, gl_159, gl_161, gl_172, gl_174, \
                         gl_176, gl_178, gl_454, gl_461, gl_463, gl_472, gl_474, gl_476, \
                         gl_487, gl_489, gl_491, gl_493 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_143[k] = -f_646 * gl_4[k]
                   - f_647 * gl_11[k]
                   + f_648 * gl_13[k]
                   - f_647 * gl_22[k]
                   + f_649 * gl_24[k]
                   - f_650 * gl_26[k]
                   - f_646 * gl_37[k]
                   + f_648 * gl_39[k]
                   - f_650 * gl_41[k]
                   + f_651 * gl_43[k]
                   + f_652 * gl_139[k]
                   + f_653 * gl_146[k]
                   - f_654 * gl_148[k]
                   + f_653 * gl_157[k]
                   - f_655 * gl_159[k]
                   + f_656 * gl_161[k]
                   + f_652 * gl_172[k]
                   - f_654 * gl_174[k]
                   + f_656 * gl_176[k]
                   - f_657 * gl_178[k]
                   - f_646 * gl_454[k]
                   - f_647 * gl_461[k]
                   + f_648 * gl_463[k]
                   - f_647 * gl_472[k]
                   + f_649 * gl_474[k]
                   - f_650 * gl_476[k]
                   - f_646 * gl_487[k]
                   + f_648 * gl_489[k]
                   - f_650 * gl_491[k]
                   + f_651 * gl_493[k];
    }

#pragma omp simd aligned(gl_0, gl_3, gl_5, gl_10, gl_12, gl_14, gl_21, gl_23, gl_25, gl_27, \
                         gl_36, gl_38, gl_40, gl_42, gl_44, gl_135, gl_138, gl_140, gl_145, \
                         gl_147, gl_149, gl_156, gl_158, gl_160, gl_162, gl_171, gl_173, \
                         gl_175, gl_177, gl_179, gl_450, gl_453, gl_455, gl_460, gl_462, \
                         gl_464, gl_471, gl_473, gl_475, gl_477, gl_486, gl_488, gl_490, \
                         gl_492, gl_494 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_144[k] = f_658 * gl_0[k]
                   + f_39 * gl_3[k]
                   - f_659 * gl_5[k]
                   + f_660 * gl_10[k]
                   - f_648 * gl_12[k]
                   + f_648 * gl_14[k]
                   + f_39 * gl_21[k]
                   - f_648 * gl_23[k]
                   + f_649 * gl_25[k]
                   - f_661 * gl_27[k]
                   + f_658 * gl_36[k]
                   - f_659 * gl_38[k]
                   + f_648 * gl_40[k]
                   - f_661 * gl_42[k]
                   + f_662 * gl_44[k]
                   - f_660 * gl_135[k]
                   - f_42 * gl_138[k]
                   + f_649 * gl_140[k]
                   - f_647 * gl_145[k]
                   + f_654 * gl_147[k]
                   - f_654 * gl_149[k]
                   - f_42 * gl_156[k]
                   + f_654 * gl_158[k]
                   - f_655 * gl_160[k]
                   + f_663 * gl_162[k]
                   - f_660 * gl_171[k]
                   + f_649 * gl_173[k]
                   - f_654 * gl_175[k]
                   + f_663 * gl_177[k]
                   - f_651 * gl_179[k]
                   + f_658 * gl_450[k]
                   + f_39 * gl_453[k]
                   - f_659 * gl_455[k]
                   + f_660 * gl_460[k]
                   - f_648 * gl_462[k]
                   + f_648 * gl_464[k]
                   + f_39 * gl_471[k]
                   - f_648 * gl_473[k]
                   + f_649 * gl_475[k]
                   - f_661 * gl_477[k]
                   + f_658 * gl_486[k]
                   - f_659 * gl_488[k]
                   + f_648 * gl_490[k]
                   - f_661 * gl_492[k]
                   + f_662 * gl_494[k];
    }

#pragma omp simd aligned(gl_2, gl_7, gl_9, gl_16, gl_18, gl_20, gl_29, gl_31, gl_33, gl_35, \
                         gl_137, gl_142, gl_144, gl_151, gl_153, gl_155, gl_164, gl_166, \
                         gl_168, gl_170, gl_452, gl_457, gl_459, gl_466, gl_468, gl_470, \
                         gl_479, gl_481, gl_483, gl_485 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_145[k] = -f_646 * gl_2[k]
                   - f_647 * gl_7[k]
                   + f_648 * gl_9[k]
                   - f_647 * gl_16[k]
                   + f_649 * gl_18[k]
                   - f_650 * gl_20[k]
                   - f_646 * gl_29[k]
                   + f_648 * gl_31[k]
                   - f_650 * gl_33[k]
                   + f_651 * gl_35[k]
                   + f_652 * gl_137[k]
                   + f_653 * gl_142[k]
                   - f_654 * gl_144[k]
                   + f_653 * gl_151[k]
                   - f_655 * gl_153[k]
                   + f_656 * gl_155[k]
                   + f_652 * gl_164[k]
                   - f_654 * gl_166[k]
                   + f_656 * gl_168[k]
                   - f_657 * gl_170[k]
                   - f_646 * gl_452[k]
                   - f_647 * gl_457[k]
                   + f_648 * gl_459[k]
                   - f_647 * gl_466[k]
                   + f_649 * gl_468[k]
                   - f_650 * gl_470[k]
                   - f_646 * gl_479[k]
                   + f_648 * gl_481[k]
                   - f_650 * gl_483[k]
                   + f_651 * gl_485[k];
    }

#pragma omp simd aligned(gl_0, gl_3, gl_5, gl_12, gl_14, gl_21, gl_23, gl_27, gl_36, gl_38, \
                         gl_40, gl_42, gl_135, gl_138, gl_140, gl_147, gl_149, gl_156, gl_158, \
                         gl_162, gl_171, gl_173, gl_175, gl_177, gl_450, gl_453, gl_455, \
                         gl_462, gl_464, gl_471, gl_473, gl_477, gl_486, gl_488, gl_490, \
                         gl_492 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_146[k] = -f_664 * gl_0[k]
                   - f_635 * gl_3[k]
                   + f_665 * gl_5[k]
                   + f_665 * gl_12[k]
                   - f_666 * gl_14[k]
                   + f_635 * gl_21[k]
                   - f_665 * gl_23[k]
                   + f_667 * gl_27[k]
                   + f_664 * gl_36[k]
                   - f_665 * gl_38[k]
                   + f_666 * gl_40[k]
                   - f_667 * gl_42[k]
                   + f_636 * gl_135[k]
                   + f_640 * gl_138[k]
                   - f_668 * gl_140[k]
                   - f_668 * gl_147[k]
                   + f_30 * gl_149[k]
                   - f_640 * gl_156[k]
                   + f_668 * gl_158[k]
                   - f_669 * gl_162[k]
                   - f_636 * gl_171[k]
                   + f_668 * gl_173[k]
                   - f_30 * gl_175[k]
                   + f_669 * gl_177[k]
                   - f_664 * gl_450[k]
                   - f_635 * gl_453[k]
                   + f_665 * gl_455[k]
                   + f_665 * gl_462[k]
                   - f_666 * gl_464[k]
                   + f_635 * gl_471[k]
                   - f_665 * gl_473[k]
                   + f_667 * gl_477[k]
                   + f_664 * gl_486[k]
                   - f_665 * gl_488[k]
                   + f_666 * gl_490[k]
                   - f_667 * gl_492[k];
    }

#pragma omp simd aligned(gl_2, gl_7, gl_9, gl_16, gl_18, gl_20, gl_29, gl_31, gl_33, gl_137, \
                         gl_142, gl_144, gl_151, gl_153, gl_155, gl_164, gl_166, gl_168, \
                         gl_452, gl_457, gl_459, gl_466, gl_468, gl_470, gl_479, gl_481, \
                         gl_483 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_147[k] = f_623 * gl_2[k]
                   - f_623 * gl_7[k]
                   - f_626 * gl_9[k]
                   - f_622 * gl_16[k]
                   + f_624 * gl_18[k]
                   + f_627 * gl_20[k]
                   - f_621 * gl_29[k]
                   + f_20 * gl_31[k]
                   - f_625 * gl_33[k]
                   - f_631 * gl_137[k]
                   + f_631 * gl_142[k]
                   + f_633 * gl_144[k]
                   + f_629 * gl_151[k]
                   - f_21 * gl_153[k]
                   - f_634 * gl_155[k]
                   + f_628 * gl_164[k]
                   - f_630 * gl_166[k]
                   + f_632 * gl_168[k]
                   + f_623 * gl_452[k]
                   - f_623 * gl_457[k]
                   - f_626 * gl_459[k]
                   - f_622 * gl_466[k]
                   + f_624 * gl_468[k]
                   + f_627 * gl_470[k]
                   - f_621 * gl_479[k]
                   + f_20 * gl_481[k]
                   - f_625 * gl_483[k];
    }

#pragma omp simd aligned(gl_0, gl_3, gl_5, gl_10, gl_12, gl_14, gl_21, gl_23, gl_25, gl_36, \
                         gl_38, gl_40, gl_135, gl_138, gl_140, gl_145, gl_147, gl_149, gl_156, \
                         gl_158, gl_160, gl_171, gl_173, gl_175, gl_450, gl_453, gl_455, \
                         gl_460, gl_462, gl_464, gl_471, gl_473, gl_475, gl_486, gl_488, \
                         gl_490 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_148[k] = f_670 * gl_0[k]
                   - f_49 * gl_3[k]
                   - f_619 * gl_5[k]
                   - f_671 * gl_10[k]
                   + f_672 * gl_12[k]
                   + f_51 * gl_14[k]
                   - f_49 * gl_21[k]
                   + f_672 * gl_23[k]
                   - f_673 * gl_25[k]
                   + f_670 * gl_36[k]
                   - f_619 * gl_38[k]
                   + f_51 * gl_40[k]
                   - f_674 * gl_135[k]
                   + f_619 * gl_138[k]
                   + f_675 * gl_140[k]
                   + f_676 * gl_145[k]
                   - f_677 * gl_147[k]
                   - f_673 * gl_149[k]
                   + f_619 * gl_156[k]
                   - f_677 * gl_158[k]
                   + f_678 * gl_160[k]
                   - f_674 * gl_171[k]
                   + f_675 * gl_173[k]
                   - f_673 * gl_175[k]
                   + f_670 * gl_450[k]
                   - f_49 * gl_453[k]
                   - f_619 * gl_455[k]
                   - f_671 * gl_460[k]
                   + f_672 * gl_462[k]
                   + f_51 * gl_464[k]
                   - f_49 * gl_471[k]
                   + f_672 * gl_473[k]
                   - f_673 * gl_475[k]
                   + f_670 * gl_486[k]
                   - f_619 * gl_488[k]
                   + f_51 * gl_490[k];
    }

#pragma omp simd aligned(gl_2, gl_7, gl_9, gl_16, gl_18, gl_29, gl_31, gl_137, gl_142, gl_144, \
                         gl_151, gl_153, gl_164, gl_166, gl_452, gl_457, gl_459, gl_466, \
                         gl_468, gl_479, gl_481 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_149[k] = -f_616 * gl_2[k]
                   + f_380 * gl_7[k]
                   + f_14 * gl_9[k]
                   + f_615 * gl_16[k]
                   - f_389 * gl_18[k]
                   - f_615 * gl_29[k]
                   + f_10 * gl_31[k]
                   + f_372 * gl_137[k]
                   - f_617 * gl_142[k]
                   - f_384 * gl_144[k]
                   - f_382 * gl_151[k]
                   + f_618 * gl_153[k]
                   + f_382 * gl_164[k]
                   - f_385 * gl_166[k]
                   - f_616 * gl_452[k]
                   + f_380 * gl_457[k]
                   + f_14 * gl_459[k]
                   + f_615 * gl_466[k]
                   - f_389 * gl_468[k]
                   - f_615 * gl_479[k]
                   + f_10 * gl_481[k];
    }

#pragma omp simd aligned(gl_0, gl_3, gl_5, gl_12, gl_21, gl_23, gl_36, gl_38, gl_135, gl_138, \
                         gl_140, gl_147, gl_156, gl_158, gl_171, gl_173, gl_450, gl_453, \
                         gl_455, gl_462, gl_471, gl_473, gl_486, \
                         gl_488 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_150[k] = -f_679 * gl_0[k]
                   + f_609 * gl_3[k]
                   + f_609 * gl_5[k]
                   - f_680 * gl_12[k]
                   - f_609 * gl_21[k]
                   + f_680 * gl_23[k]
                   + f_679 * gl_36[k]
                   - f_609 * gl_38[k]
                   + f_608 * gl_135[k]
                   - f_610 * gl_138[k]
                   - f_610 * gl_140[k]
                   + f_681 * gl_147[k]
                   + f_610 * gl_156[k]
                   - f_681 * gl_158[k]
                   - f_608 * gl_171[k]
                   + f_610 * gl_173[k]
                   - f_679 * gl_450[k]
                   + f_609 * gl_453[k]
                   + f_609 * gl_455[k]
                   - f_680 * gl_462[k]
                   - f_609 * gl_471[k]
                   + f_680 * gl_473[k]
                   + f_679 * gl_486[k]
                   - f_609 * gl_488[k];
    }

#pragma omp simd aligned(gl_2, gl_7, gl_16, gl_29, gl_137, gl_142, gl_151, gl_164, gl_452, \
                         gl_457, gl_466, gl_479 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_151[k] = f_57 * gl_2[k]
                   - f_604 * gl_7[k]
                   + f_603 * gl_16[k]
                   - f_602 * gl_29[k]
                   - f_413 * gl_137[k]
                   + f_607 * gl_142[k]
                   - f_606 * gl_151[k]
                   + f_605 * gl_164[k]
                   + f_57 * gl_452[k]
                   - f_604 * gl_457[k]
                   + f_603 * gl_466[k]
                   - f_602 * gl_479[k];
    }

#pragma omp simd aligned(gl_0, gl_3, gl_10, gl_21, gl_36, gl_135, gl_138, gl_145, gl_156, \
                         gl_171, gl_450, gl_453, gl_460, gl_471, \
                         gl_486 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_152[k] = f_682 * gl_0[k]
                   - f_602 * gl_3[k]
                   + f_683 * gl_10[k]
                   - f_602 * gl_21[k]
                   + f_682 * gl_36[k]
                   - f_684 * gl_135[k]
                   + f_605 * gl_138[k]
                   - f_685 * gl_145[k]
                   + f_605 * gl_156[k]
                   - f_684 * gl_171[k]
                   + f_682 * gl_450[k]
                   - f_602 * gl_453[k]
                   + f_683 * gl_460[k]
                   - f_602 * gl_471[k]
                   + f_682 * gl_486[k];
    }
}

}  // namespace simdtrf
