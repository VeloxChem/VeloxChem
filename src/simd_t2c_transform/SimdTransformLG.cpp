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


#include "SimdTransformLG.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_lg(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t lg,
             const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 0.46875 * std::sqrt(1001.0);
    const auto f_1 = 3.28125 * std::sqrt(1001.0);
    const auto f_2 = 0.703125 * std::sqrt(2002.0);
    const auto f_3 = 0.234375 * std::sqrt(2002.0);
    const auto f_4 = 4.921875 * std::sqrt(2002.0);
    const auto f_5 = 1.640625 * std::sqrt(2002.0);
    const auto f_6 = 0.46875 * std::sqrt(143.0);
    const auto f_7 = 2.8125 * std::sqrt(143.0);
    const auto f_8 = 3.28125 * std::sqrt(143.0);
    const auto f_9 = 19.6875 * std::sqrt(143.0);
    const auto f_10 = 0.703125 * std::sqrt(286.0);
    const auto f_11 = 0.9375 * std::sqrt(286.0);
    const auto f_12 = 4.921875 * std::sqrt(286.0);
    const auto f_13 = 6.5625 * std::sqrt(286.0);
    const auto f_14 = 0.0703125 * std::sqrt(715.0);
    const auto f_15 = 0.140625 * std::sqrt(715.0);
    const auto f_16 = 0.5625 * std::sqrt(715.0);
    const auto f_17 = 0.1875 * std::sqrt(715.0);
    const auto f_18 = 0.4921875 * std::sqrt(715.0);
    const auto f_19 = 0.984375 * std::sqrt(715.0);
    const auto f_20 = 3.9375 * std::sqrt(715.0);
    const auto f_21 = 1.3125 * std::sqrt(715.0);
    const auto f_22 = 0.234375 * std::sqrt(143.0);
    const auto f_23 = 1.40625 * std::sqrt(143.0);
    const auto f_24 = 1.640625 * std::sqrt(143.0);
    const auto f_25 = 9.84375 * std::sqrt(143.0);
    const auto f_26 = 0.1171875 * std::sqrt(1001.0);
    const auto f_27 = 0.703125 * std::sqrt(1001.0);
    const auto f_28 = 0.8203125 * std::sqrt(1001.0);
    const auto f_29 = 4.921875 * std::sqrt(1001.0);
    const auto f_30 = 1.640625 * std::sqrt(1001.0);
    const auto f_31 = 8.203125 * std::sqrt(1001.0);
    const auto f_32 = 0.234375 * std::sqrt(1001.0);
    const auto f_33 = 2.4609375 * std::sqrt(2002.0);
    const auto f_34 = 0.8203125 * std::sqrt(2002.0);
    const auto f_35 = 12.3046875 * std::sqrt(2002.0);
    const auto f_36 = 4.1015625 * std::sqrt(2002.0);
    const auto f_37 = 7.3828125 * std::sqrt(2002.0);
    const auto f_38 = 0.3515625 * std::sqrt(2002.0);
    const auto f_39 = 0.1171875 * std::sqrt(2002.0);
    const auto f_40 = 8.203125 * std::sqrt(143.0);
    const auto f_41 = 49.21875 * std::sqrt(143.0);
    const auto f_42 = 4.921875 * std::sqrt(143.0);
    const auto f_43 = 29.53125 * std::sqrt(143.0);
    const auto f_44 = 2.4609375 * std::sqrt(286.0);
    const auto f_45 = 3.28125 * std::sqrt(286.0);
    const auto f_46 = 12.3046875 * std::sqrt(286.0);
    const auto f_47 = 16.40625 * std::sqrt(286.0);
    const auto f_48 = 7.3828125 * std::sqrt(286.0);
    const auto f_49 = 9.84375 * std::sqrt(286.0);
    const auto f_50 = 0.3515625 * std::sqrt(286.0);
    const auto f_51 = 0.46875 * std::sqrt(286.0);
    const auto f_52 = 0.24609375 * std::sqrt(715.0);
    const auto f_53 = 1.96875 * std::sqrt(715.0);
    const auto f_54 = 0.65625 * std::sqrt(715.0);
    const auto f_55 = 1.23046875 * std::sqrt(715.0);
    const auto f_56 = 2.4609375 * std::sqrt(715.0);
    const auto f_57 = 9.84375 * std::sqrt(715.0);
    const auto f_58 = 3.28125 * std::sqrt(715.0);
    const auto f_59 = 0.73828125 * std::sqrt(715.0);
    const auto f_60 = 1.4765625 * std::sqrt(715.0);
    const auto f_61 = 5.90625 * std::sqrt(715.0);
    const auto f_62 = 0.03515625 * std::sqrt(715.0);
    const auto f_63 = 0.28125 * std::sqrt(715.0);
    const auto f_64 = 0.09375 * std::sqrt(715.0);
    const auto f_65 = 0.8203125 * std::sqrt(143.0);
    const auto f_66 = 4.1015625 * std::sqrt(143.0);
    const auto f_67 = 24.609375 * std::sqrt(143.0);
    const auto f_68 = 2.4609375 * std::sqrt(143.0);
    const auto f_69 = 14.765625 * std::sqrt(143.0);
    const auto f_70 = 0.1171875 * std::sqrt(143.0);
    const auto f_71 = 0.703125 * std::sqrt(143.0);
    const auto f_72 = 0.41015625 * std::sqrt(1001.0);
    const auto f_73 = 2.4609375 * std::sqrt(1001.0);
    const auto f_74 = 2.05078125 * std::sqrt(1001.0);
    const auto f_75 = 12.3046875 * std::sqrt(1001.0);
    const auto f_76 = 1.23046875 * std::sqrt(1001.0);
    const auto f_77 = 7.3828125 * std::sqrt(1001.0);
    const auto f_78 = 0.05859375 * std::sqrt(1001.0);
    const auto f_79 = 0.3515625 * std::sqrt(1001.0);
    const auto f_80 = 0.046875 * std::sqrt(30030.0);
    const auto f_81 = 0.109375 * std::sqrt(30030.0);
    const auto f_82 = 0.65625 * std::sqrt(30030.0);
    const auto f_83 = 2.1875 * std::sqrt(30030.0);
    const auto f_84 = 0.140625 * std::sqrt(15015.0);
    const auto f_85 = 0.046875 * std::sqrt(15015.0);
    const auto f_86 = 0.328125 * std::sqrt(15015.0);
    const auto f_87 = 0.109375 * std::sqrt(15015.0);
    const auto f_88 = 1.96875 * std::sqrt(15015.0);
    const auto f_89 = 0.65625 * std::sqrt(15015.0);
    const auto f_90 = 6.5625 * std::sqrt(15015.0);
    const auto f_91 = 2.1875 * std::sqrt(15015.0);
    const auto f_92 = 0.046875 * std::sqrt(4290.0);
    const auto f_93 = 0.28125 * std::sqrt(4290.0);
    const auto f_94 = 0.109375 * std::sqrt(4290.0);
    const auto f_95 = 0.65625 * std::sqrt(4290.0);
    const auto f_96 = 3.9375 * std::sqrt(4290.0);
    const auto f_97 = 2.1875 * std::sqrt(4290.0);
    const auto f_98 = 13.125 * std::sqrt(4290.0);
    const auto f_99 = 0.140625 * std::sqrt(2145.0);
    const auto f_100 = 0.1875 * std::sqrt(2145.0);
    const auto f_101 = 0.328125 * std::sqrt(2145.0);
    const auto f_102 = 0.4375 * std::sqrt(2145.0);
    const auto f_103 = 1.96875 * std::sqrt(2145.0);
    const auto f_104 = 2.625 * std::sqrt(2145.0);
    const auto f_105 = 6.5625 * std::sqrt(2145.0);
    const auto f_106 = 8.75 * std::sqrt(2145.0);
    const auto f_107 = 0.03515625 * std::sqrt(858.0);
    const auto f_108 = 0.0703125 * std::sqrt(858.0);
    const auto f_109 = 0.28125 * std::sqrt(858.0);
    const auto f_110 = 0.09375 * std::sqrt(858.0);
    const auto f_111 = 0.08203125 * std::sqrt(858.0);
    const auto f_112 = 0.1640625 * std::sqrt(858.0);
    const auto f_113 = 0.65625 * std::sqrt(858.0);
    const auto f_114 = 0.21875 * std::sqrt(858.0);
    const auto f_115 = 0.4921875 * std::sqrt(858.0);
    const auto f_116 = 0.984375 * std::sqrt(858.0);
    const auto f_117 = 3.9375 * std::sqrt(858.0);
    const auto f_118 = 1.3125 * std::sqrt(858.0);
    const auto f_119 = 1.640625 * std::sqrt(858.0);
    const auto f_120 = 3.28125 * std::sqrt(858.0);
    const auto f_121 = 13.125 * std::sqrt(858.0);
    const auto f_122 = 4.375 * std::sqrt(858.0);
    const auto f_123 = 0.0234375 * std::sqrt(4290.0);
    const auto f_124 = 0.140625 * std::sqrt(4290.0);
    const auto f_125 = 0.0546875 * std::sqrt(4290.0);
    const auto f_126 = 0.328125 * std::sqrt(4290.0);
    const auto f_127 = 1.96875 * std::sqrt(4290.0);
    const auto f_128 = 1.09375 * std::sqrt(4290.0);
    const auto f_129 = 6.5625 * std::sqrt(4290.0);
    const auto f_130 = 0.01171875 * std::sqrt(30030.0);
    const auto f_131 = 0.0703125 * std::sqrt(30030.0);
    const auto f_132 = 0.02734375 * std::sqrt(30030.0);
    const auto f_133 = 0.1640625 * std::sqrt(30030.0);
    const auto f_134 = 0.984375 * std::sqrt(30030.0);
    const auto f_135 = 0.546875 * std::sqrt(30030.0);
    const auto f_136 = 3.28125 * std::sqrt(30030.0);
    const auto f_137 = 1.640625 * std::sqrt(715.0);
    const auto f_138 = 6.5625 * std::sqrt(715.0);
    const auto f_139 = 2.953125 * std::sqrt(715.0);
    const auto f_140 = 13.125 * std::sqrt(715.0);
    const auto f_141 = 0.328125 * std::sqrt(715.0);
    const auto f_142 = 2.4609375 * std::sqrt(1430.0);
    const auto f_143 = 0.8203125 * std::sqrt(1430.0);
    const auto f_144 = 9.84375 * std::sqrt(1430.0);
    const auto f_145 = 3.28125 * std::sqrt(1430.0);
    const auto f_146 = 4.4296875 * std::sqrt(1430.0);
    const auto f_147 = 1.4765625 * std::sqrt(1430.0);
    const auto f_148 = 19.6875 * std::sqrt(1430.0);
    const auto f_149 = 6.5625 * std::sqrt(1430.0);
    const auto f_150 = 0.4921875 * std::sqrt(1430.0);
    const auto f_151 = 0.1640625 * std::sqrt(1430.0);
    const auto f_152 = 1.96875 * std::sqrt(1430.0);
    const auto f_153 = 0.65625 * std::sqrt(1430.0);
    const auto f_154 = 0.234375 * std::sqrt(5005.0);
    const auto f_155 = 1.40625 * std::sqrt(5005.0);
    const auto f_156 = 0.9375 * std::sqrt(5005.0);
    const auto f_157 = 5.625 * std::sqrt(5005.0);
    const auto f_158 = 0.421875 * std::sqrt(5005.0);
    const auto f_159 = 2.53125 * std::sqrt(5005.0);
    const auto f_160 = 1.875 * std::sqrt(5005.0);
    const auto f_161 = 11.25 * std::sqrt(5005.0);
    const auto f_162 = 0.046875 * std::sqrt(5005.0);
    const auto f_163 = 0.28125 * std::sqrt(5005.0);
    const auto f_164 = 0.1875 * std::sqrt(5005.0);
    const auto f_165 = 1.125 * std::sqrt(5005.0);
    const auto f_166 = 0.3515625 * std::sqrt(10010.0);
    const auto f_167 = 0.46875 * std::sqrt(10010.0);
    const auto f_168 = 1.40625 * std::sqrt(10010.0);
    const auto f_169 = 1.875 * std::sqrt(10010.0);
    const auto f_170 = 0.6328125 * std::sqrt(10010.0);
    const auto f_171 = 0.84375 * std::sqrt(10010.0);
    const auto f_172 = 2.8125 * std::sqrt(10010.0);
    const auto f_173 = 3.75 * std::sqrt(10010.0);
    const auto f_174 = 0.0703125 * std::sqrt(10010.0);
    const auto f_175 = 0.09375 * std::sqrt(10010.0);
    const auto f_176 = 0.28125 * std::sqrt(10010.0);
    const auto f_177 = 0.375 * std::sqrt(10010.0);
    const auto f_178 = 0.17578125 * std::sqrt(1001.0);
    const auto f_179 = 1.40625 * std::sqrt(1001.0);
    const auto f_180 = 5.625 * std::sqrt(1001.0);
    const auto f_181 = 1.875 * std::sqrt(1001.0);
    const auto f_182 = 0.31640625 * std::sqrt(1001.0);
    const auto f_183 = 0.6328125 * std::sqrt(1001.0);
    const auto f_184 = 2.53125 * std::sqrt(1001.0);
    const auto f_185 = 0.84375 * std::sqrt(1001.0);
    const auto f_186 = 2.8125 * std::sqrt(1001.0);
    const auto f_187 = 11.25 * std::sqrt(1001.0);
    const auto f_188 = 3.75 * std::sqrt(1001.0);
    const auto f_189 = 0.03515625 * std::sqrt(1001.0);
    const auto f_190 = 0.0703125 * std::sqrt(1001.0);
    const auto f_191 = 0.28125 * std::sqrt(1001.0);
    const auto f_192 = 0.09375 * std::sqrt(1001.0);
    const auto f_193 = 0.140625 * std::sqrt(1001.0);
    const auto f_194 = 1.125 * std::sqrt(1001.0);
    const auto f_195 = 0.375 * std::sqrt(1001.0);
    const auto f_196 = 0.1171875 * std::sqrt(5005.0);
    const auto f_197 = 0.703125 * std::sqrt(5005.0);
    const auto f_198 = 0.46875 * std::sqrt(5005.0);
    const auto f_199 = 2.8125 * std::sqrt(5005.0);
    const auto f_200 = 0.2109375 * std::sqrt(5005.0);
    const auto f_201 = 1.265625 * std::sqrt(5005.0);
    const auto f_202 = 0.0234375 * std::sqrt(5005.0);
    const auto f_203 = 0.140625 * std::sqrt(5005.0);
    const auto f_204 = 0.09375 * std::sqrt(5005.0);
    const auto f_205 = 0.5625 * std::sqrt(5005.0);
    const auto f_206 = 0.41015625 * std::sqrt(715.0);
    const auto f_207 = 4.4296875 * std::sqrt(715.0);
    const auto f_208 = 19.6875 * std::sqrt(715.0);
    const auto f_209 = 0.08203125 * std::sqrt(715.0);
    const auto f_210 = 0.65625 * std::sqrt(55.0);
    const auto f_211 = 15.75 * std::sqrt(55.0);
    const auto f_212 = 26.25 * std::sqrt(55.0);
    const auto f_213 = 0.984375 * std::sqrt(110.0);
    const auto f_214 = 0.328125 * std::sqrt(110.0);
    const auto f_215 = 23.625 * std::sqrt(110.0);
    const auto f_216 = 7.875 * std::sqrt(110.0);
    const auto f_217 = 39.375 * std::sqrt(110.0);
    const auto f_218 = 13.125 * std::sqrt(110.0);
    const auto f_219 = 0.09375 * std::sqrt(385.0);
    const auto f_220 = 0.5625 * std::sqrt(385.0);
    const auto f_221 = 2.25 * std::sqrt(385.0);
    const auto f_222 = 13.5 * std::sqrt(385.0);
    const auto f_223 = 3.75 * std::sqrt(385.0);
    const auto f_224 = 22.5 * std::sqrt(385.0);
    const auto f_225 = 0.140625 * std::sqrt(770.0);
    const auto f_226 = 0.1875 * std::sqrt(770.0);
    const auto f_227 = 3.375 * std::sqrt(770.0);
    const auto f_228 = 4.5 * std::sqrt(770.0);
    const auto f_229 = 5.625 * std::sqrt(770.0);
    const auto f_230 = 7.5 * std::sqrt(770.0);
    const auto f_231 = 0.0703125 * std::sqrt(77.0);
    const auto f_232 = 0.140625 * std::sqrt(77.0);
    const auto f_233 = 0.5625 * std::sqrt(77.0);
    const auto f_234 = 0.1875 * std::sqrt(77.0);
    const auto f_235 = 1.6875 * std::sqrt(77.0);
    const auto f_236 = 3.375 * std::sqrt(77.0);
    const auto f_237 = 13.5 * std::sqrt(77.0);
    const auto f_238 = 4.5 * std::sqrt(77.0);
    const auto f_239 = 2.8125 * std::sqrt(77.0);
    const auto f_240 = 5.625 * std::sqrt(77.0);
    const auto f_241 = 22.5 * std::sqrt(77.0);
    const auto f_242 = 7.5 * std::sqrt(77.0);
    const auto f_243 = 0.046875 * std::sqrt(385.0);
    const auto f_244 = 0.28125 * std::sqrt(385.0);
    const auto f_245 = 1.125 * std::sqrt(385.0);
    const auto f_246 = 6.75 * std::sqrt(385.0);
    const auto f_247 = 1.875 * std::sqrt(385.0);
    const auto f_248 = 11.25 * std::sqrt(385.0);
    const auto f_249 = 0.1640625 * std::sqrt(55.0);
    const auto f_250 = 0.984375 * std::sqrt(55.0);
    const auto f_251 = 3.9375 * std::sqrt(55.0);
    const auto f_252 = 23.625 * std::sqrt(55.0);
    const auto f_253 = 6.5625 * std::sqrt(55.0);
    const auto f_254 = 39.375 * std::sqrt(55.0);
    const auto f_255 = 4.921875 * std::sqrt(33.0);
    const auto f_256 = 8.203125 * std::sqrt(33.0);
    const auto f_257 = 32.8125 * std::sqrt(33.0);
    const auto f_258 = 1.640625 * std::sqrt(33.0);
    const auto f_259 = 21.875 * std::sqrt(33.0);
    const auto f_260 = 26.25 * std::sqrt(33.0);
    const auto f_261 = 10.9375 * std::sqrt(33.0);
    const auto f_262 = 8.75 * std::sqrt(33.0);
    const auto f_263 = 7.3828125 * std::sqrt(66.0);
    const auto f_264 = 2.4609375 * std::sqrt(66.0);
    const auto f_265 = 12.3046875 * std::sqrt(66.0);
    const auto f_266 = 4.1015625 * std::sqrt(66.0);
    const auto f_267 = 49.21875 * std::sqrt(66.0);
    const auto f_268 = 16.40625 * std::sqrt(66.0);
    const auto f_269 = 0.8203125 * std::sqrt(66.0);
    const auto f_270 = 32.8125 * std::sqrt(66.0);
    const auto f_271 = 10.9375 * std::sqrt(66.0);
    const auto f_272 = 39.375 * std::sqrt(66.0);
    const auto f_273 = 13.125 * std::sqrt(66.0);
    const auto f_274 = 5.46875 * std::sqrt(66.0);
    const auto f_275 = 4.375 * std::sqrt(66.0);
    const auto f_276 = 0.703125 * std::sqrt(231.0);
    const auto f_277 = 4.21875 * std::sqrt(231.0);
    const auto f_278 = 1.171875 * std::sqrt(231.0);
    const auto f_279 = 7.03125 * std::sqrt(231.0);
    const auto f_280 = 4.6875 * std::sqrt(231.0);
    const auto f_281 = 28.125 * std::sqrt(231.0);
    const auto f_282 = 0.234375 * std::sqrt(231.0);
    const auto f_283 = 1.40625 * std::sqrt(231.0);
    const auto f_284 = 3.125 * std::sqrt(231.0);
    const auto f_285 = 18.75 * std::sqrt(231.0);
    const auto f_286 = 3.75 * std::sqrt(231.0);
    const auto f_287 = 22.5 * std::sqrt(231.0);
    const auto f_288 = 1.5625 * std::sqrt(231.0);
    const auto f_289 = 9.375 * std::sqrt(231.0);
    const auto f_290 = 1.25 * std::sqrt(231.0);
    const auto f_291 = 7.5 * std::sqrt(231.0);
    const auto f_292 = 1.0546875 * std::sqrt(462.0);
    const auto f_293 = 1.40625 * std::sqrt(462.0);
    const auto f_294 = 1.7578125 * std::sqrt(462.0);
    const auto f_295 = 2.34375 * std::sqrt(462.0);
    const auto f_296 = 7.03125 * std::sqrt(462.0);
    const auto f_297 = 9.375 * std::sqrt(462.0);
    const auto f_298 = 0.3515625 * std::sqrt(462.0);
    const auto f_299 = 0.46875 * std::sqrt(462.0);
    const auto f_300 = 4.6875 * std::sqrt(462.0);
    const auto f_301 = 6.25 * std::sqrt(462.0);
    const auto f_302 = 5.625 * std::sqrt(462.0);
    const auto f_303 = 7.5 * std::sqrt(462.0);
    const auto f_304 = 3.125 * std::sqrt(462.0);
    const auto f_305 = 1.875 * std::sqrt(462.0);
    const auto f_306 = 2.5 * std::sqrt(462.0);
    const auto f_307 = 0.10546875 * std::sqrt(1155.0);
    const auto f_308 = 0.2109375 * std::sqrt(1155.0);
    const auto f_309 = 0.84375 * std::sqrt(1155.0);
    const auto f_310 = 0.28125 * std::sqrt(1155.0);
    const auto f_311 = 0.17578125 * std::sqrt(1155.0);
    const auto f_312 = 0.3515625 * std::sqrt(1155.0);
    const auto f_313 = 1.40625 * std::sqrt(1155.0);
    const auto f_314 = 0.46875 * std::sqrt(1155.0);
    const auto f_315 = 0.703125 * std::sqrt(1155.0);
    const auto f_316 = 5.625 * std::sqrt(1155.0);
    const auto f_317 = 1.875 * std::sqrt(1155.0);
    const auto f_318 = 0.03515625 * std::sqrt(1155.0);
    const auto f_319 = 0.0703125 * std::sqrt(1155.0);
    const auto f_320 = 0.09375 * std::sqrt(1155.0);
    const auto f_321 = 0.9375 * std::sqrt(1155.0);
    const auto f_322 = 3.75 * std::sqrt(1155.0);
    const auto f_323 = 1.25 * std::sqrt(1155.0);
    const auto f_324 = 0.5625 * std::sqrt(1155.0);
    const auto f_325 = 1.125 * std::sqrt(1155.0);
    const auto f_326 = 4.5 * std::sqrt(1155.0);
    const auto f_327 = 1.5 * std::sqrt(1155.0);
    const auto f_328 = 0.234375 * std::sqrt(1155.0);
    const auto f_329 = 0.625 * std::sqrt(1155.0);
    const auto f_330 = 0.1875 * std::sqrt(1155.0);
    const auto f_331 = 0.375 * std::sqrt(1155.0);
    const auto f_332 = 0.5 * std::sqrt(1155.0);
    const auto f_333 = 0.3515625 * std::sqrt(231.0);
    const auto f_334 = 2.109375 * std::sqrt(231.0);
    const auto f_335 = 0.5859375 * std::sqrt(231.0);
    const auto f_336 = 3.515625 * std::sqrt(231.0);
    const auto f_337 = 2.34375 * std::sqrt(231.0);
    const auto f_338 = 14.0625 * std::sqrt(231.0);
    const auto f_339 = 0.1171875 * std::sqrt(231.0);
    const auto f_340 = 1.875 * std::sqrt(231.0);
    const auto f_341 = 11.25 * std::sqrt(231.0);
    const auto f_342 = 0.78125 * std::sqrt(231.0);
    const auto f_343 = 0.625 * std::sqrt(231.0);
    const auto f_344 = 1.23046875 * std::sqrt(33.0);
    const auto f_345 = 7.3828125 * std::sqrt(33.0);
    const auto f_346 = 2.05078125 * std::sqrt(33.0);
    const auto f_347 = 12.3046875 * std::sqrt(33.0);
    const auto f_348 = 49.21875 * std::sqrt(33.0);
    const auto f_349 = 0.41015625 * std::sqrt(33.0);
    const auto f_350 = 2.4609375 * std::sqrt(33.0);
    const auto f_351 = 5.46875 * std::sqrt(33.0);
    const auto f_352 = 6.5625 * std::sqrt(33.0);
    const auto f_353 = 39.375 * std::sqrt(33.0);
    const auto f_354 = 2.734375 * std::sqrt(33.0);
    const auto f_355 = 16.40625 * std::sqrt(33.0);
    const auto f_356 = 2.1875 * std::sqrt(33.0);
    const auto f_357 = 13.125 * std::sqrt(33.0);
    const auto f_358 = 1.640625 * std::sqrt(2.0);
    const auto f_359 = 4.921875 * std::sqrt(2.0);
    const auto f_360 = 49.21875 * std::sqrt(2.0);
    const auto f_361 = 98.4375 * std::sqrt(2.0);
    const auto f_362 = 131.25 * std::sqrt(2.0);
    const auto f_363 = 52.5 * std::sqrt(2.0);
    const auto f_364 = 0.234375 * std::sqrt(14.0);
    const auto f_365 = 1.40625 * std::sqrt(14.0);
    const auto f_366 = 0.703125 * std::sqrt(14.0);
    const auto f_367 = 4.21875 * std::sqrt(14.0);
    const auto f_368 = 7.03125 * std::sqrt(14.0);
    const auto f_369 = 42.1875 * std::sqrt(14.0);
    const auto f_370 = 14.0625 * std::sqrt(14.0);
    const auto f_371 = 84.375 * std::sqrt(14.0);
    const auto f_372 = 18.75 * std::sqrt(14.0);
    const auto f_373 = 112.5 * std::sqrt(14.0);
    const auto f_374 = 7.5 * std::sqrt(14.0);
    const auto f_375 = 45.0 * std::sqrt(14.0);
    const auto f_376 = 0.703125 * std::sqrt(7.0);
    const auto f_377 = 0.9375 * std::sqrt(7.0);
    const auto f_378 = 2.109375 * std::sqrt(7.0);
    const auto f_379 = 2.8125 * std::sqrt(7.0);
    const auto f_380 = 21.09375 * std::sqrt(7.0);
    const auto f_381 = 28.125 * std::sqrt(7.0);
    const auto f_382 = 42.1875 * std::sqrt(7.0);
    const auto f_383 = 56.25 * std::sqrt(7.0);
    const auto f_384 = 75.0 * std::sqrt(7.0);
    const auto f_385 = 22.5 * std::sqrt(7.0);
    const auto f_386 = 30.0 * std::sqrt(7.0);
    const auto f_387 = 0.03515625 * std::sqrt(70.0);
    const auto f_388 = 0.0703125 * std::sqrt(70.0);
    const auto f_389 = 0.28125 * std::sqrt(70.0);
    const auto f_390 = 0.09375 * std::sqrt(70.0);
    const auto f_391 = 0.10546875 * std::sqrt(70.0);
    const auto f_392 = 0.2109375 * std::sqrt(70.0);
    const auto f_393 = 0.84375 * std::sqrt(70.0);
    const auto f_394 = 1.0546875 * std::sqrt(70.0);
    const auto f_395 = 2.109375 * std::sqrt(70.0);
    const auto f_396 = 8.4375 * std::sqrt(70.0);
    const auto f_397 = 2.8125 * std::sqrt(70.0);
    const auto f_398 = 4.21875 * std::sqrt(70.0);
    const auto f_399 = 16.875 * std::sqrt(70.0);
    const auto f_400 = 5.625 * std::sqrt(70.0);
    const auto f_401 = 22.5 * std::sqrt(70.0);
    const auto f_402 = 7.5 * std::sqrt(70.0);
    const auto f_403 = 1.125 * std::sqrt(70.0);
    const auto f_404 = 2.25 * std::sqrt(70.0);
    const auto f_405 = 9.0 * std::sqrt(70.0);
    const auto f_406 = 3.0 * std::sqrt(70.0);
    const auto f_407 = 0.1171875 * std::sqrt(14.0);
    const auto f_408 = 0.3515625 * std::sqrt(14.0);
    const auto f_409 = 2.109375 * std::sqrt(14.0);
    const auto f_410 = 3.515625 * std::sqrt(14.0);
    const auto f_411 = 21.09375 * std::sqrt(14.0);
    const auto f_412 = 9.375 * std::sqrt(14.0);
    const auto f_413 = 56.25 * std::sqrt(14.0);
    const auto f_414 = 3.75 * std::sqrt(14.0);
    const auto f_415 = 22.5 * std::sqrt(14.0);
    const auto f_416 = 0.41015625 * std::sqrt(2.0);
    const auto f_417 = 2.4609375 * std::sqrt(2.0);
    const auto f_418 = 1.23046875 * std::sqrt(2.0);
    const auto f_419 = 7.3828125 * std::sqrt(2.0);
    const auto f_420 = 12.3046875 * std::sqrt(2.0);
    const auto f_421 = 73.828125 * std::sqrt(2.0);
    const auto f_422 = 24.609375 * std::sqrt(2.0);
    const auto f_423 = 147.65625 * std::sqrt(2.0);
    const auto f_424 = 32.8125 * std::sqrt(2.0);
    const auto f_425 = 196.875 * std::sqrt(2.0);
    const auto f_426 = 13.125 * std::sqrt(2.0);
    const auto f_427 = 78.75 * std::sqrt(2.0);
    const auto f_428 = 1.640625 * std::sqrt(35.0);
    const auto f_429 = 4.921875 * std::sqrt(35.0);
    const auto f_430 = 13.125 * std::sqrt(35.0);
    const auto f_431 = 26.25 * std::sqrt(35.0);
    const auto f_432 = 15.75 * std::sqrt(35.0);
    const auto f_433 = 3.0 * std::sqrt(35.0);
    const auto f_434 = 2.4609375 * std::sqrt(70.0);
    const auto f_435 = 0.8203125 * std::sqrt(70.0);
    const auto f_436 = 7.3828125 * std::sqrt(70.0);
    const auto f_437 = 19.6875 * std::sqrt(70.0);
    const auto f_438 = 6.5625 * std::sqrt(70.0);
    const auto f_439 = 39.375 * std::sqrt(70.0);
    const auto f_440 = 13.125 * std::sqrt(70.0);
    const auto f_441 = 23.625 * std::sqrt(70.0);
    const auto f_442 = 7.875 * std::sqrt(70.0);
    const auto f_443 = 4.5 * std::sqrt(70.0);
    const auto f_444 = 1.5 * std::sqrt(70.0);
    const auto f_445 = 1.640625 * std::sqrt(5.0);
    const auto f_446 = 9.84375 * std::sqrt(5.0);
    const auto f_447 = 4.921875 * std::sqrt(5.0);
    const auto f_448 = 29.53125 * std::sqrt(5.0);
    const auto f_449 = 13.125 * std::sqrt(5.0);
    const auto f_450 = 78.75 * std::sqrt(5.0);
    const auto f_451 = 26.25 * std::sqrt(5.0);
    const auto f_452 = 157.5 * std::sqrt(5.0);
    const auto f_453 = 15.75 * std::sqrt(5.0);
    const auto f_454 = 94.5 * std::sqrt(5.0);
    const auto f_455 = 3.0 * std::sqrt(5.0);
    const auto f_456 = 18.0 * std::sqrt(5.0);
    const auto f_457 = 2.4609375 * std::sqrt(10.0);
    const auto f_458 = 3.28125 * std::sqrt(10.0);
    const auto f_459 = 7.3828125 * std::sqrt(10.0);
    const auto f_460 = 9.84375 * std::sqrt(10.0);
    const auto f_461 = 19.6875 * std::sqrt(10.0);
    const auto f_462 = 26.25 * std::sqrt(10.0);
    const auto f_463 = 39.375 * std::sqrt(10.0);
    const auto f_464 = 52.5 * std::sqrt(10.0);
    const auto f_465 = 23.625 * std::sqrt(10.0);
    const auto f_466 = 31.5 * std::sqrt(10.0);
    const auto f_467 = 4.5 * std::sqrt(10.0);
    const auto f_468 = 6.0 * std::sqrt(10.0);
    const auto f_469 = 0.8203125 * std::sqrt(5.0);
    const auto f_470 = 2.4609375 * std::sqrt(5.0);
    const auto f_471 = 14.765625 * std::sqrt(5.0);
    const auto f_472 = 6.5625 * std::sqrt(5.0);
    const auto f_473 = 39.375 * std::sqrt(5.0);
    const auto f_474 = 7.875 * std::sqrt(5.0);
    const auto f_475 = 47.25 * std::sqrt(5.0);
    const auto f_476 = 1.5 * std::sqrt(5.0);
    const auto f_477 = 9.0 * std::sqrt(5.0);
    const auto f_478 = 0.41015625 * std::sqrt(35.0);
    const auto f_479 = 2.4609375 * std::sqrt(35.0);
    const auto f_480 = 1.23046875 * std::sqrt(35.0);
    const auto f_481 = 7.3828125 * std::sqrt(35.0);
    const auto f_482 = 3.28125 * std::sqrt(35.0);
    const auto f_483 = 19.6875 * std::sqrt(35.0);
    const auto f_484 = 6.5625 * std::sqrt(35.0);
    const auto f_485 = 39.375 * std::sqrt(35.0);
    const auto f_486 = 3.9375 * std::sqrt(35.0);
    const auto f_487 = 23.625 * std::sqrt(35.0);
    const auto f_488 = 0.75 * std::sqrt(35.0);
    const auto f_489 = 4.5 * std::sqrt(35.0);
    const auto f_490 = 0.13671875 * std::sqrt(35.0);
    const auto f_491 = 0.546875 * std::sqrt(35.0);
    const auto f_492 = 4.375 * std::sqrt(35.0);
    const auto f_493 = 0.8203125 * std::sqrt(35.0);
    const auto f_494 = 7.0 * std::sqrt(35.0);
    const auto f_495 = 0.5 * std::sqrt(35.0);
    const auto f_496 = 0.205078125 * std::sqrt(70.0);
    const auto f_497 = 0.068359375 * std::sqrt(70.0);
    const auto f_498 = 0.2734375 * std::sqrt(70.0);
    const auto f_499 = 2.1875 * std::sqrt(70.0);
    const auto f_500 = 1.23046875 * std::sqrt(70.0);
    const auto f_501 = 0.41015625 * std::sqrt(70.0);
    const auto f_502 = 10.5 * std::sqrt(70.0);
    const auto f_503 = 3.5 * std::sqrt(70.0);
    const auto f_504 = 0.75 * std::sqrt(70.0);
    const auto f_505 = 0.25 * std::sqrt(70.0);
    const auto f_506 = 0.13671875 * std::sqrt(5.0);
    const auto f_507 = 0.546875 * std::sqrt(5.0);
    const auto f_508 = 3.28125 * std::sqrt(5.0);
    const auto f_509 = 4.375 * std::sqrt(5.0);
    const auto f_510 = 7.0 * std::sqrt(5.0);
    const auto f_511 = 42.0 * std::sqrt(5.0);
    const auto f_512 = 0.5 * std::sqrt(5.0);
    const auto f_513 = 0.205078125 * std::sqrt(10.0);
    const auto f_514 = 0.2734375 * std::sqrt(10.0);
    const auto f_515 = 0.8203125 * std::sqrt(10.0);
    const auto f_516 = 1.09375 * std::sqrt(10.0);
    const auto f_517 = 6.5625 * std::sqrt(10.0);
    const auto f_518 = 8.75 * std::sqrt(10.0);
    const auto f_519 = 1.23046875 * std::sqrt(10.0);
    const auto f_520 = 1.640625 * std::sqrt(10.0);
    const auto f_521 = 10.5 * std::sqrt(10.0);
    const auto f_522 = 14.0 * std::sqrt(10.0);
    const auto f_523 = 0.75 * std::sqrt(10.0);
    const auto f_524 = std::sqrt(10.0);
    const auto f_525 = 0.068359375 * std::sqrt(5.0);
    const auto f_526 = 0.41015625 * std::sqrt(5.0);
    const auto f_527 = 0.2734375 * std::sqrt(5.0);
    const auto f_528 = 2.1875 * std::sqrt(5.0);
    const auto f_529 = 3.5 * std::sqrt(5.0);
    const auto f_530 = 21.0 * std::sqrt(5.0);
    const auto f_531 = 0.25 * std::sqrt(5.0);
    const auto f_532 = 0.0341796875 * std::sqrt(35.0);
    const auto f_533 = 0.205078125 * std::sqrt(35.0);
    const auto f_534 = 1.09375 * std::sqrt(35.0);
    const auto f_535 = 1.75 * std::sqrt(35.0);
    const auto f_536 = 10.5 * std::sqrt(35.0);
    const auto f_537 = 0.125 * std::sqrt(35.0);
    const auto f_538 = 0.8203125 * std::sqrt(2.0);
    const auto f_539 = 65.625 * std::sqrt(2.0);
    const auto f_540 = 26.25 * std::sqrt(2.0);
    const auto f_541 = 0.3515625 * std::sqrt(7.0);
    const auto f_542 = 0.46875 * std::sqrt(7.0);
    const auto f_543 = 10.546875 * std::sqrt(7.0);
    const auto f_544 = 14.0625 * std::sqrt(7.0);
    const auto f_545 = 37.5 * std::sqrt(7.0);
    const auto f_546 = 11.25 * std::sqrt(7.0);
    const auto f_547 = 15.0 * std::sqrt(7.0);
    const auto f_548 = 0.017578125 * std::sqrt(70.0);
    const auto f_549 = 0.140625 * std::sqrt(70.0);
    const auto f_550 = 0.046875 * std::sqrt(70.0);
    const auto f_551 = 0.52734375 * std::sqrt(70.0);
    const auto f_552 = 1.40625 * std::sqrt(70.0);
    const auto f_553 = 11.25 * std::sqrt(70.0);
    const auto f_554 = 3.75 * std::sqrt(70.0);
    const auto f_555 = 0.5625 * std::sqrt(70.0);
    const auto f_556 = 0.05859375 * std::sqrt(14.0);
    const auto f_557 = 1.7578125 * std::sqrt(14.0);
    const auto f_558 = 10.546875 * std::sqrt(14.0);
    const auto f_559 = 4.6875 * std::sqrt(14.0);
    const auto f_560 = 28.125 * std::sqrt(14.0);
    const auto f_561 = 1.875 * std::sqrt(14.0);
    const auto f_562 = 11.25 * std::sqrt(14.0);
    const auto f_563 = 0.205078125 * std::sqrt(2.0);
    const auto f_564 = 6.15234375 * std::sqrt(2.0);
    const auto f_565 = 36.9140625 * std::sqrt(2.0);
    const auto f_566 = 16.40625 * std::sqrt(2.0);
    const auto f_567 = 6.5625 * std::sqrt(2.0);
    const auto f_568 = 39.375 * std::sqrt(2.0);
    const auto f_569 = 1.640625 * std::sqrt(55.0);
    const auto f_570 = 19.6875 * std::sqrt(55.0);
    const auto f_571 = 0.24609375 * std::sqrt(110.0);
    const auto f_572 = 0.08203125 * std::sqrt(110.0);
    const auto f_573 = 5.90625 * std::sqrt(110.0);
    const auto f_574 = 1.96875 * std::sqrt(110.0);
    const auto f_575 = 2.4609375 * std::sqrt(110.0);
    const auto f_576 = 0.8203125 * std::sqrt(110.0);
    const auto f_577 = 29.53125 * std::sqrt(110.0);
    const auto f_578 = 9.84375 * std::sqrt(110.0);
    const auto f_579 = 3.28125 * std::sqrt(110.0);
    const auto f_580 = 59.0625 * std::sqrt(110.0);
    const auto f_581 = 19.6875 * std::sqrt(110.0);
    const auto f_582 = 0.0234375 * std::sqrt(385.0);
    const auto f_583 = 0.140625 * std::sqrt(385.0);
    const auto f_584 = 3.375 * std::sqrt(385.0);
    const auto f_585 = 0.234375 * std::sqrt(385.0);
    const auto f_586 = 1.40625 * std::sqrt(385.0);
    const auto f_587 = 2.8125 * std::sqrt(385.0);
    const auto f_588 = 16.875 * std::sqrt(385.0);
    const auto f_589 = 0.9375 * std::sqrt(385.0);
    const auto f_590 = 5.625 * std::sqrt(385.0);
    const auto f_591 = 33.75 * std::sqrt(385.0);
    const auto f_592 = 0.03515625 * std::sqrt(770.0);
    const auto f_593 = 0.046875 * std::sqrt(770.0);
    const auto f_594 = 0.84375 * std::sqrt(770.0);
    const auto f_595 = 1.125 * std::sqrt(770.0);
    const auto f_596 = 0.3515625 * std::sqrt(770.0);
    const auto f_597 = 0.46875 * std::sqrt(770.0);
    const auto f_598 = 4.21875 * std::sqrt(770.0);
    const auto f_599 = 1.40625 * std::sqrt(770.0);
    const auto f_600 = 1.875 * std::sqrt(770.0);
    const auto f_601 = 8.4375 * std::sqrt(770.0);
    const auto f_602 = 11.25 * std::sqrt(770.0);
    const auto f_603 = 0.017578125 * std::sqrt(77.0);
    const auto f_604 = 0.03515625 * std::sqrt(77.0);
    const auto f_605 = 0.046875 * std::sqrt(77.0);
    const auto f_606 = 0.421875 * std::sqrt(77.0);
    const auto f_607 = 0.84375 * std::sqrt(77.0);
    const auto f_608 = 1.125 * std::sqrt(77.0);
    const auto f_609 = 0.17578125 * std::sqrt(77.0);
    const auto f_610 = 0.3515625 * std::sqrt(77.0);
    const auto f_611 = 1.40625 * std::sqrt(77.0);
    const auto f_612 = 0.46875 * std::sqrt(77.0);
    const auto f_613 = 2.109375 * std::sqrt(77.0);
    const auto f_614 = 4.21875 * std::sqrt(77.0);
    const auto f_615 = 16.875 * std::sqrt(77.0);
    const auto f_616 = 0.703125 * std::sqrt(77.0);
    const auto f_617 = 1.875 * std::sqrt(77.0);
    const auto f_618 = 8.4375 * std::sqrt(77.0);
    const auto f_619 = 33.75 * std::sqrt(77.0);
    const auto f_620 = 11.25 * std::sqrt(77.0);
    const auto f_621 = 0.01171875 * std::sqrt(385.0);
    const auto f_622 = 0.0703125 * std::sqrt(385.0);
    const auto f_623 = 1.6875 * std::sqrt(385.0);
    const auto f_624 = 0.1171875 * std::sqrt(385.0);
    const auto f_625 = 0.703125 * std::sqrt(385.0);
    const auto f_626 = 8.4375 * std::sqrt(385.0);
    const auto f_627 = 0.46875 * std::sqrt(385.0);
    const auto f_628 = 0.041015625 * std::sqrt(55.0);
    const auto f_629 = 0.24609375 * std::sqrt(55.0);
    const auto f_630 = 5.90625 * std::sqrt(55.0);
    const auto f_631 = 0.41015625 * std::sqrt(55.0);
    const auto f_632 = 2.4609375 * std::sqrt(55.0);
    const auto f_633 = 4.921875 * std::sqrt(55.0);
    const auto f_634 = 29.53125 * std::sqrt(55.0);
    const auto f_635 = 9.84375 * std::sqrt(55.0);
    const auto f_636 = 59.0625 * std::sqrt(55.0);
    const auto f_637 = 0.0078125 * std::sqrt(30030.0);
    const auto f_638 = 1.640625 * std::sqrt(30030.0);
    const auto f_639 = 0.0234375 * std::sqrt(15015.0);
    const auto f_640 = 0.0078125 * std::sqrt(15015.0);
    const auto f_641 = 4.921875 * std::sqrt(15015.0);
    const auto f_642 = 1.640625 * std::sqrt(15015.0);
    const auto f_643 = 0.0078125 * std::sqrt(4290.0);
    const auto f_644 = 1.640625 * std::sqrt(4290.0);
    const auto f_645 = 9.84375 * std::sqrt(4290.0);
    const auto f_646 = 0.0234375 * std::sqrt(2145.0);
    const auto f_647 = 0.03125 * std::sqrt(2145.0);
    const auto f_648 = 4.921875 * std::sqrt(2145.0);
    const auto f_649 = 0.005859375 * std::sqrt(858.0);
    const auto f_650 = 0.01171875 * std::sqrt(858.0);
    const auto f_651 = 0.046875 * std::sqrt(858.0);
    const auto f_652 = 0.015625 * std::sqrt(858.0);
    const auto f_653 = 1.23046875 * std::sqrt(858.0);
    const auto f_654 = 2.4609375 * std::sqrt(858.0);
    const auto f_655 = 9.84375 * std::sqrt(858.0);
    const auto f_656 = 0.00390625 * std::sqrt(4290.0);
    const auto f_657 = 0.8203125 * std::sqrt(4290.0);
    const auto f_658 = 4.921875 * std::sqrt(4290.0);
    const auto f_659 = 0.001953125 * std::sqrt(30030.0);
    const auto f_660 = 0.41015625 * std::sqrt(30030.0);
    const auto f_661 = 2.4609375 * std::sqrt(30030.0);
    const auto f_662 = 4.1015625 * std::sqrt(1001.0);
    const auto f_663 = 0.087890625 * std::sqrt(2002.0);
    const auto f_664 = 0.029296875 * std::sqrt(2002.0);
    const auto f_665 = 6.15234375 * std::sqrt(2002.0);
    const auto f_666 = 2.05078125 * std::sqrt(2002.0);
    const auto f_667 = 0.05859375 * std::sqrt(143.0);
    const auto f_668 = 0.3515625 * std::sqrt(143.0);
    const auto f_669 = 0.087890625 * std::sqrt(286.0);
    const auto f_670 = 0.1171875 * std::sqrt(286.0);
    const auto f_671 = 6.15234375 * std::sqrt(286.0);
    const auto f_672 = 8.203125 * std::sqrt(286.0);
    const auto f_673 = 0.0087890625 * std::sqrt(715.0);
    const auto f_674 = 0.017578125 * std::sqrt(715.0);
    const auto f_675 = 0.0234375 * std::sqrt(715.0);
    const auto f_676 = 0.615234375 * std::sqrt(715.0);
    const auto f_677 = 4.921875 * std::sqrt(715.0);
    const auto f_678 = 0.029296875 * std::sqrt(143.0);
    const auto f_679 = 0.17578125 * std::sqrt(143.0);
    const auto f_680 = 2.05078125 * std::sqrt(143.0);
    const auto f_681 = 12.3046875 * std::sqrt(143.0);
    const auto f_682 = 0.0146484375 * std::sqrt(1001.0);
    const auto f_683 = 0.087890625 * std::sqrt(1001.0);
    const auto f_684 = 1.025390625 * std::sqrt(1001.0);
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

    const auto *lg_0 = buffer.data(lg + 0);
    const auto *lg_1 = buffer.data(lg + 1);
    const auto *lg_2 = buffer.data(lg + 2);
    const auto *lg_3 = buffer.data(lg + 3);
    const auto *lg_4 = buffer.data(lg + 4);
    const auto *lg_5 = buffer.data(lg + 5);
    const auto *lg_6 = buffer.data(lg + 6);
    const auto *lg_7 = buffer.data(lg + 7);
    const auto *lg_8 = buffer.data(lg + 8);
    const auto *lg_9 = buffer.data(lg + 9);
    const auto *lg_10 = buffer.data(lg + 10);
    const auto *lg_11 = buffer.data(lg + 11);
    const auto *lg_12 = buffer.data(lg + 12);
    const auto *lg_13 = buffer.data(lg + 13);
    const auto *lg_14 = buffer.data(lg + 14);
    const auto *lg_15 = buffer.data(lg + 15);
    const auto *lg_16 = buffer.data(lg + 16);
    const auto *lg_17 = buffer.data(lg + 17);
    const auto *lg_18 = buffer.data(lg + 18);
    const auto *lg_19 = buffer.data(lg + 19);
    const auto *lg_20 = buffer.data(lg + 20);
    const auto *lg_21 = buffer.data(lg + 21);
    const auto *lg_22 = buffer.data(lg + 22);
    const auto *lg_23 = buffer.data(lg + 23);
    const auto *lg_24 = buffer.data(lg + 24);
    const auto *lg_25 = buffer.data(lg + 25);
    const auto *lg_26 = buffer.data(lg + 26);
    const auto *lg_27 = buffer.data(lg + 27);
    const auto *lg_28 = buffer.data(lg + 28);
    const auto *lg_29 = buffer.data(lg + 29);
    const auto *lg_30 = buffer.data(lg + 30);
    const auto *lg_31 = buffer.data(lg + 31);
    const auto *lg_32 = buffer.data(lg + 32);
    const auto *lg_33 = buffer.data(lg + 33);
    const auto *lg_34 = buffer.data(lg + 34);
    const auto *lg_35 = buffer.data(lg + 35);
    const auto *lg_36 = buffer.data(lg + 36);
    const auto *lg_37 = buffer.data(lg + 37);
    const auto *lg_38 = buffer.data(lg + 38);
    const auto *lg_39 = buffer.data(lg + 39);
    const auto *lg_40 = buffer.data(lg + 40);
    const auto *lg_41 = buffer.data(lg + 41);
    const auto *lg_42 = buffer.data(lg + 42);
    const auto *lg_43 = buffer.data(lg + 43);
    const auto *lg_44 = buffer.data(lg + 44);
    const auto *lg_45 = buffer.data(lg + 45);
    const auto *lg_46 = buffer.data(lg + 46);
    const auto *lg_47 = buffer.data(lg + 47);
    const auto *lg_48 = buffer.data(lg + 48);
    const auto *lg_49 = buffer.data(lg + 49);
    const auto *lg_50 = buffer.data(lg + 50);
    const auto *lg_51 = buffer.data(lg + 51);
    const auto *lg_52 = buffer.data(lg + 52);
    const auto *lg_53 = buffer.data(lg + 53);
    const auto *lg_54 = buffer.data(lg + 54);
    const auto *lg_55 = buffer.data(lg + 55);
    const auto *lg_56 = buffer.data(lg + 56);
    const auto *lg_57 = buffer.data(lg + 57);
    const auto *lg_58 = buffer.data(lg + 58);
    const auto *lg_59 = buffer.data(lg + 59);
    const auto *lg_60 = buffer.data(lg + 60);
    const auto *lg_61 = buffer.data(lg + 61);
    const auto *lg_62 = buffer.data(lg + 62);
    const auto *lg_63 = buffer.data(lg + 63);
    const auto *lg_64 = buffer.data(lg + 64);
    const auto *lg_65 = buffer.data(lg + 65);
    const auto *lg_66 = buffer.data(lg + 66);
    const auto *lg_67 = buffer.data(lg + 67);
    const auto *lg_68 = buffer.data(lg + 68);
    const auto *lg_69 = buffer.data(lg + 69);
    const auto *lg_70 = buffer.data(lg + 70);
    const auto *lg_71 = buffer.data(lg + 71);
    const auto *lg_72 = buffer.data(lg + 72);
    const auto *lg_73 = buffer.data(lg + 73);
    const auto *lg_74 = buffer.data(lg + 74);
    const auto *lg_75 = buffer.data(lg + 75);
    const auto *lg_76 = buffer.data(lg + 76);
    const auto *lg_77 = buffer.data(lg + 77);
    const auto *lg_78 = buffer.data(lg + 78);
    const auto *lg_79 = buffer.data(lg + 79);
    const auto *lg_80 = buffer.data(lg + 80);
    const auto *lg_81 = buffer.data(lg + 81);
    const auto *lg_82 = buffer.data(lg + 82);
    const auto *lg_83 = buffer.data(lg + 83);
    const auto *lg_84 = buffer.data(lg + 84);
    const auto *lg_85 = buffer.data(lg + 85);
    const auto *lg_86 = buffer.data(lg + 86);
    const auto *lg_87 = buffer.data(lg + 87);
    const auto *lg_88 = buffer.data(lg + 88);
    const auto *lg_89 = buffer.data(lg + 89);
    const auto *lg_90 = buffer.data(lg + 90);
    const auto *lg_91 = buffer.data(lg + 91);
    const auto *lg_92 = buffer.data(lg + 92);
    const auto *lg_93 = buffer.data(lg + 93);
    const auto *lg_94 = buffer.data(lg + 94);
    const auto *lg_95 = buffer.data(lg + 95);
    const auto *lg_96 = buffer.data(lg + 96);
    const auto *lg_97 = buffer.data(lg + 97);
    const auto *lg_98 = buffer.data(lg + 98);
    const auto *lg_99 = buffer.data(lg + 99);
    const auto *lg_100 = buffer.data(lg + 100);
    const auto *lg_101 = buffer.data(lg + 101);
    const auto *lg_102 = buffer.data(lg + 102);
    const auto *lg_103 = buffer.data(lg + 103);
    const auto *lg_104 = buffer.data(lg + 104);
    const auto *lg_105 = buffer.data(lg + 105);
    const auto *lg_106 = buffer.data(lg + 106);
    const auto *lg_107 = buffer.data(lg + 107);
    const auto *lg_108 = buffer.data(lg + 108);
    const auto *lg_109 = buffer.data(lg + 109);
    const auto *lg_110 = buffer.data(lg + 110);
    const auto *lg_111 = buffer.data(lg + 111);
    const auto *lg_112 = buffer.data(lg + 112);
    const auto *lg_113 = buffer.data(lg + 113);
    const auto *lg_114 = buffer.data(lg + 114);
    const auto *lg_115 = buffer.data(lg + 115);
    const auto *lg_116 = buffer.data(lg + 116);
    const auto *lg_117 = buffer.data(lg + 117);
    const auto *lg_118 = buffer.data(lg + 118);
    const auto *lg_119 = buffer.data(lg + 119);
    const auto *lg_120 = buffer.data(lg + 120);
    const auto *lg_121 = buffer.data(lg + 121);
    const auto *lg_122 = buffer.data(lg + 122);
    const auto *lg_123 = buffer.data(lg + 123);
    const auto *lg_124 = buffer.data(lg + 124);
    const auto *lg_125 = buffer.data(lg + 125);
    const auto *lg_126 = buffer.data(lg + 126);
    const auto *lg_127 = buffer.data(lg + 127);
    const auto *lg_128 = buffer.data(lg + 128);
    const auto *lg_129 = buffer.data(lg + 129);
    const auto *lg_130 = buffer.data(lg + 130);
    const auto *lg_131 = buffer.data(lg + 131);
    const auto *lg_132 = buffer.data(lg + 132);
    const auto *lg_133 = buffer.data(lg + 133);
    const auto *lg_134 = buffer.data(lg + 134);
    const auto *lg_135 = buffer.data(lg + 135);
    const auto *lg_136 = buffer.data(lg + 136);
    const auto *lg_137 = buffer.data(lg + 137);
    const auto *lg_138 = buffer.data(lg + 138);
    const auto *lg_139 = buffer.data(lg + 139);
    const auto *lg_140 = buffer.data(lg + 140);
    const auto *lg_141 = buffer.data(lg + 141);
    const auto *lg_142 = buffer.data(lg + 142);
    const auto *lg_143 = buffer.data(lg + 143);
    const auto *lg_144 = buffer.data(lg + 144);
    const auto *lg_145 = buffer.data(lg + 145);
    const auto *lg_146 = buffer.data(lg + 146);
    const auto *lg_147 = buffer.data(lg + 147);
    const auto *lg_148 = buffer.data(lg + 148);
    const auto *lg_149 = buffer.data(lg + 149);
    const auto *lg_150 = buffer.data(lg + 150);
    const auto *lg_151 = buffer.data(lg + 151);
    const auto *lg_152 = buffer.data(lg + 152);
    const auto *lg_153 = buffer.data(lg + 153);
    const auto *lg_154 = buffer.data(lg + 154);
    const auto *lg_155 = buffer.data(lg + 155);
    const auto *lg_156 = buffer.data(lg + 156);
    const auto *lg_157 = buffer.data(lg + 157);
    const auto *lg_158 = buffer.data(lg + 158);
    const auto *lg_159 = buffer.data(lg + 159);
    const auto *lg_160 = buffer.data(lg + 160);
    const auto *lg_161 = buffer.data(lg + 161);
    const auto *lg_162 = buffer.data(lg + 162);
    const auto *lg_163 = buffer.data(lg + 163);
    const auto *lg_164 = buffer.data(lg + 164);
    const auto *lg_165 = buffer.data(lg + 165);
    const auto *lg_166 = buffer.data(lg + 166);
    const auto *lg_167 = buffer.data(lg + 167);
    const auto *lg_168 = buffer.data(lg + 168);
    const auto *lg_169 = buffer.data(lg + 169);
    const auto *lg_170 = buffer.data(lg + 170);
    const auto *lg_171 = buffer.data(lg + 171);
    const auto *lg_172 = buffer.data(lg + 172);
    const auto *lg_173 = buffer.data(lg + 173);
    const auto *lg_174 = buffer.data(lg + 174);
    const auto *lg_175 = buffer.data(lg + 175);
    const auto *lg_176 = buffer.data(lg + 176);
    const auto *lg_177 = buffer.data(lg + 177);
    const auto *lg_178 = buffer.data(lg + 178);
    const auto *lg_179 = buffer.data(lg + 179);
    const auto *lg_180 = buffer.data(lg + 180);
    const auto *lg_181 = buffer.data(lg + 181);
    const auto *lg_182 = buffer.data(lg + 182);
    const auto *lg_183 = buffer.data(lg + 183);
    const auto *lg_184 = buffer.data(lg + 184);
    const auto *lg_185 = buffer.data(lg + 185);
    const auto *lg_186 = buffer.data(lg + 186);
    const auto *lg_187 = buffer.data(lg + 187);
    const auto *lg_188 = buffer.data(lg + 188);
    const auto *lg_189 = buffer.data(lg + 189);
    const auto *lg_190 = buffer.data(lg + 190);
    const auto *lg_191 = buffer.data(lg + 191);
    const auto *lg_192 = buffer.data(lg + 192);
    const auto *lg_193 = buffer.data(lg + 193);
    const auto *lg_194 = buffer.data(lg + 194);
    const auto *lg_195 = buffer.data(lg + 195);
    const auto *lg_196 = buffer.data(lg + 196);
    const auto *lg_197 = buffer.data(lg + 197);
    const auto *lg_198 = buffer.data(lg + 198);
    const auto *lg_199 = buffer.data(lg + 199);
    const auto *lg_200 = buffer.data(lg + 200);
    const auto *lg_201 = buffer.data(lg + 201);
    const auto *lg_202 = buffer.data(lg + 202);
    const auto *lg_203 = buffer.data(lg + 203);
    const auto *lg_204 = buffer.data(lg + 204);
    const auto *lg_205 = buffer.data(lg + 205);
    const auto *lg_206 = buffer.data(lg + 206);
    const auto *lg_207 = buffer.data(lg + 207);
    const auto *lg_208 = buffer.data(lg + 208);
    const auto *lg_209 = buffer.data(lg + 209);
    const auto *lg_210 = buffer.data(lg + 210);
    const auto *lg_211 = buffer.data(lg + 211);
    const auto *lg_212 = buffer.data(lg + 212);
    const auto *lg_213 = buffer.data(lg + 213);
    const auto *lg_214 = buffer.data(lg + 214);
    const auto *lg_215 = buffer.data(lg + 215);
    const auto *lg_216 = buffer.data(lg + 216);
    const auto *lg_217 = buffer.data(lg + 217);
    const auto *lg_218 = buffer.data(lg + 218);
    const auto *lg_219 = buffer.data(lg + 219);
    const auto *lg_220 = buffer.data(lg + 220);
    const auto *lg_221 = buffer.data(lg + 221);
    const auto *lg_222 = buffer.data(lg + 222);
    const auto *lg_223 = buffer.data(lg + 223);
    const auto *lg_224 = buffer.data(lg + 224);
    const auto *lg_225 = buffer.data(lg + 225);
    const auto *lg_226 = buffer.data(lg + 226);
    const auto *lg_227 = buffer.data(lg + 227);
    const auto *lg_228 = buffer.data(lg + 228);
    const auto *lg_229 = buffer.data(lg + 229);
    const auto *lg_230 = buffer.data(lg + 230);
    const auto *lg_231 = buffer.data(lg + 231);
    const auto *lg_232 = buffer.data(lg + 232);
    const auto *lg_233 = buffer.data(lg + 233);
    const auto *lg_234 = buffer.data(lg + 234);
    const auto *lg_235 = buffer.data(lg + 235);
    const auto *lg_236 = buffer.data(lg + 236);
    const auto *lg_237 = buffer.data(lg + 237);
    const auto *lg_238 = buffer.data(lg + 238);
    const auto *lg_239 = buffer.data(lg + 239);
    const auto *lg_240 = buffer.data(lg + 240);
    const auto *lg_241 = buffer.data(lg + 241);
    const auto *lg_242 = buffer.data(lg + 242);
    const auto *lg_243 = buffer.data(lg + 243);
    const auto *lg_244 = buffer.data(lg + 244);
    const auto *lg_245 = buffer.data(lg + 245);
    const auto *lg_246 = buffer.data(lg + 246);
    const auto *lg_247 = buffer.data(lg + 247);
    const auto *lg_248 = buffer.data(lg + 248);
    const auto *lg_249 = buffer.data(lg + 249);
    const auto *lg_250 = buffer.data(lg + 250);
    const auto *lg_251 = buffer.data(lg + 251);
    const auto *lg_252 = buffer.data(lg + 252);
    const auto *lg_253 = buffer.data(lg + 253);
    const auto *lg_254 = buffer.data(lg + 254);
    const auto *lg_255 = buffer.data(lg + 255);
    const auto *lg_256 = buffer.data(lg + 256);
    const auto *lg_257 = buffer.data(lg + 257);
    const auto *lg_258 = buffer.data(lg + 258);
    const auto *lg_259 = buffer.data(lg + 259);
    const auto *lg_260 = buffer.data(lg + 260);
    const auto *lg_261 = buffer.data(lg + 261);
    const auto *lg_262 = buffer.data(lg + 262);
    const auto *lg_263 = buffer.data(lg + 263);
    const auto *lg_264 = buffer.data(lg + 264);
    const auto *lg_265 = buffer.data(lg + 265);
    const auto *lg_266 = buffer.data(lg + 266);
    const auto *lg_267 = buffer.data(lg + 267);
    const auto *lg_268 = buffer.data(lg + 268);
    const auto *lg_269 = buffer.data(lg + 269);
    const auto *lg_270 = buffer.data(lg + 270);
    const auto *lg_271 = buffer.data(lg + 271);
    const auto *lg_272 = buffer.data(lg + 272);
    const auto *lg_273 = buffer.data(lg + 273);
    const auto *lg_274 = buffer.data(lg + 274);
    const auto *lg_275 = buffer.data(lg + 275);
    const auto *lg_276 = buffer.data(lg + 276);
    const auto *lg_277 = buffer.data(lg + 277);
    const auto *lg_278 = buffer.data(lg + 278);
    const auto *lg_279 = buffer.data(lg + 279);
    const auto *lg_280 = buffer.data(lg + 280);
    const auto *lg_281 = buffer.data(lg + 281);
    const auto *lg_282 = buffer.data(lg + 282);
    const auto *lg_283 = buffer.data(lg + 283);
    const auto *lg_284 = buffer.data(lg + 284);
    const auto *lg_285 = buffer.data(lg + 285);
    const auto *lg_286 = buffer.data(lg + 286);
    const auto *lg_287 = buffer.data(lg + 287);
    const auto *lg_288 = buffer.data(lg + 288);
    const auto *lg_289 = buffer.data(lg + 289);
    const auto *lg_290 = buffer.data(lg + 290);
    const auto *lg_291 = buffer.data(lg + 291);
    const auto *lg_292 = buffer.data(lg + 292);
    const auto *lg_293 = buffer.data(lg + 293);
    const auto *lg_294 = buffer.data(lg + 294);
    const auto *lg_295 = buffer.data(lg + 295);
    const auto *lg_296 = buffer.data(lg + 296);
    const auto *lg_297 = buffer.data(lg + 297);
    const auto *lg_298 = buffer.data(lg + 298);
    const auto *lg_299 = buffer.data(lg + 299);
    const auto *lg_300 = buffer.data(lg + 300);
    const auto *lg_301 = buffer.data(lg + 301);
    const auto *lg_302 = buffer.data(lg + 302);
    const auto *lg_303 = buffer.data(lg + 303);
    const auto *lg_304 = buffer.data(lg + 304);
    const auto *lg_305 = buffer.data(lg + 305);
    const auto *lg_306 = buffer.data(lg + 306);
    const auto *lg_307 = buffer.data(lg + 307);
    const auto *lg_308 = buffer.data(lg + 308);
    const auto *lg_309 = buffer.data(lg + 309);
    const auto *lg_310 = buffer.data(lg + 310);
    const auto *lg_311 = buffer.data(lg + 311);
    const auto *lg_312 = buffer.data(lg + 312);
    const auto *lg_313 = buffer.data(lg + 313);
    const auto *lg_314 = buffer.data(lg + 314);
    const auto *lg_315 = buffer.data(lg + 315);
    const auto *lg_316 = buffer.data(lg + 316);
    const auto *lg_317 = buffer.data(lg + 317);
    const auto *lg_318 = buffer.data(lg + 318);
    const auto *lg_319 = buffer.data(lg + 319);
    const auto *lg_320 = buffer.data(lg + 320);
    const auto *lg_321 = buffer.data(lg + 321);
    const auto *lg_322 = buffer.data(lg + 322);
    const auto *lg_323 = buffer.data(lg + 323);
    const auto *lg_324 = buffer.data(lg + 324);
    const auto *lg_325 = buffer.data(lg + 325);
    const auto *lg_326 = buffer.data(lg + 326);
    const auto *lg_327 = buffer.data(lg + 327);
    const auto *lg_328 = buffer.data(lg + 328);
    const auto *lg_329 = buffer.data(lg + 329);
    const auto *lg_330 = buffer.data(lg + 330);
    const auto *lg_331 = buffer.data(lg + 331);
    const auto *lg_332 = buffer.data(lg + 332);
    const auto *lg_333 = buffer.data(lg + 333);
    const auto *lg_334 = buffer.data(lg + 334);
    const auto *lg_335 = buffer.data(lg + 335);
    const auto *lg_336 = buffer.data(lg + 336);
    const auto *lg_337 = buffer.data(lg + 337);
    const auto *lg_338 = buffer.data(lg + 338);
    const auto *lg_339 = buffer.data(lg + 339);
    const auto *lg_340 = buffer.data(lg + 340);
    const auto *lg_341 = buffer.data(lg + 341);
    const auto *lg_342 = buffer.data(lg + 342);
    const auto *lg_343 = buffer.data(lg + 343);
    const auto *lg_344 = buffer.data(lg + 344);
    const auto *lg_345 = buffer.data(lg + 345);
    const auto *lg_346 = buffer.data(lg + 346);
    const auto *lg_347 = buffer.data(lg + 347);
    const auto *lg_348 = buffer.data(lg + 348);
    const auto *lg_349 = buffer.data(lg + 349);
    const auto *lg_350 = buffer.data(lg + 350);
    const auto *lg_351 = buffer.data(lg + 351);
    const auto *lg_352 = buffer.data(lg + 352);
    const auto *lg_353 = buffer.data(lg + 353);
    const auto *lg_354 = buffer.data(lg + 354);
    const auto *lg_355 = buffer.data(lg + 355);
    const auto *lg_356 = buffer.data(lg + 356);
    const auto *lg_357 = buffer.data(lg + 357);
    const auto *lg_358 = buffer.data(lg + 358);
    const auto *lg_359 = buffer.data(lg + 359);
    const auto *lg_360 = buffer.data(lg + 360);
    const auto *lg_361 = buffer.data(lg + 361);
    const auto *lg_362 = buffer.data(lg + 362);
    const auto *lg_363 = buffer.data(lg + 363);
    const auto *lg_364 = buffer.data(lg + 364);
    const auto *lg_365 = buffer.data(lg + 365);
    const auto *lg_366 = buffer.data(lg + 366);
    const auto *lg_367 = buffer.data(lg + 367);
    const auto *lg_368 = buffer.data(lg + 368);
    const auto *lg_369 = buffer.data(lg + 369);
    const auto *lg_370 = buffer.data(lg + 370);
    const auto *lg_371 = buffer.data(lg + 371);
    const auto *lg_372 = buffer.data(lg + 372);
    const auto *lg_373 = buffer.data(lg + 373);
    const auto *lg_374 = buffer.data(lg + 374);
    const auto *lg_375 = buffer.data(lg + 375);
    const auto *lg_376 = buffer.data(lg + 376);
    const auto *lg_377 = buffer.data(lg + 377);
    const auto *lg_378 = buffer.data(lg + 378);
    const auto *lg_379 = buffer.data(lg + 379);
    const auto *lg_380 = buffer.data(lg + 380);
    const auto *lg_381 = buffer.data(lg + 381);
    const auto *lg_382 = buffer.data(lg + 382);
    const auto *lg_383 = buffer.data(lg + 383);
    const auto *lg_384 = buffer.data(lg + 384);
    const auto *lg_385 = buffer.data(lg + 385);
    const auto *lg_386 = buffer.data(lg + 386);
    const auto *lg_387 = buffer.data(lg + 387);
    const auto *lg_388 = buffer.data(lg + 388);
    const auto *lg_389 = buffer.data(lg + 389);
    const auto *lg_390 = buffer.data(lg + 390);
    const auto *lg_391 = buffer.data(lg + 391);
    const auto *lg_392 = buffer.data(lg + 392);
    const auto *lg_393 = buffer.data(lg + 393);
    const auto *lg_394 = buffer.data(lg + 394);
    const auto *lg_395 = buffer.data(lg + 395);
    const auto *lg_396 = buffer.data(lg + 396);
    const auto *lg_397 = buffer.data(lg + 397);
    const auto *lg_398 = buffer.data(lg + 398);
    const auto *lg_399 = buffer.data(lg + 399);
    const auto *lg_400 = buffer.data(lg + 400);
    const auto *lg_401 = buffer.data(lg + 401);
    const auto *lg_402 = buffer.data(lg + 402);
    const auto *lg_403 = buffer.data(lg + 403);
    const auto *lg_404 = buffer.data(lg + 404);
    const auto *lg_405 = buffer.data(lg + 405);
    const auto *lg_406 = buffer.data(lg + 406);
    const auto *lg_407 = buffer.data(lg + 407);
    const auto *lg_408 = buffer.data(lg + 408);
    const auto *lg_409 = buffer.data(lg + 409);
    const auto *lg_410 = buffer.data(lg + 410);
    const auto *lg_411 = buffer.data(lg + 411);
    const auto *lg_412 = buffer.data(lg + 412);
    const auto *lg_413 = buffer.data(lg + 413);
    const auto *lg_414 = buffer.data(lg + 414);
    const auto *lg_415 = buffer.data(lg + 415);
    const auto *lg_416 = buffer.data(lg + 416);
    const auto *lg_417 = buffer.data(lg + 417);
    const auto *lg_418 = buffer.data(lg + 418);
    const auto *lg_419 = buffer.data(lg + 419);
    const auto *lg_420 = buffer.data(lg + 420);
    const auto *lg_421 = buffer.data(lg + 421);
    const auto *lg_422 = buffer.data(lg + 422);
    const auto *lg_423 = buffer.data(lg + 423);
    const auto *lg_424 = buffer.data(lg + 424);
    const auto *lg_425 = buffer.data(lg + 425);
    const auto *lg_426 = buffer.data(lg + 426);
    const auto *lg_427 = buffer.data(lg + 427);
    const auto *lg_428 = buffer.data(lg + 428);
    const auto *lg_429 = buffer.data(lg + 429);
    const auto *lg_430 = buffer.data(lg + 430);
    const auto *lg_431 = buffer.data(lg + 431);
    const auto *lg_432 = buffer.data(lg + 432);
    const auto *lg_433 = buffer.data(lg + 433);
    const auto *lg_434 = buffer.data(lg + 434);
    const auto *lg_435 = buffer.data(lg + 435);
    const auto *lg_436 = buffer.data(lg + 436);
    const auto *lg_437 = buffer.data(lg + 437);
    const auto *lg_438 = buffer.data(lg + 438);
    const auto *lg_439 = buffer.data(lg + 439);
    const auto *lg_440 = buffer.data(lg + 440);
    const auto *lg_441 = buffer.data(lg + 441);
    const auto *lg_442 = buffer.data(lg + 442);
    const auto *lg_443 = buffer.data(lg + 443);
    const auto *lg_444 = buffer.data(lg + 444);
    const auto *lg_445 = buffer.data(lg + 445);
    const auto *lg_446 = buffer.data(lg + 446);
    const auto *lg_447 = buffer.data(lg + 447);
    const auto *lg_448 = buffer.data(lg + 448);
    const auto *lg_449 = buffer.data(lg + 449);
    const auto *lg_450 = buffer.data(lg + 450);
    const auto *lg_451 = buffer.data(lg + 451);
    const auto *lg_452 = buffer.data(lg + 452);
    const auto *lg_453 = buffer.data(lg + 453);
    const auto *lg_454 = buffer.data(lg + 454);
    const auto *lg_455 = buffer.data(lg + 455);
    const auto *lg_456 = buffer.data(lg + 456);
    const auto *lg_457 = buffer.data(lg + 457);
    const auto *lg_458 = buffer.data(lg + 458);
    const auto *lg_459 = buffer.data(lg + 459);
    const auto *lg_460 = buffer.data(lg + 460);
    const auto *lg_461 = buffer.data(lg + 461);
    const auto *lg_462 = buffer.data(lg + 462);
    const auto *lg_463 = buffer.data(lg + 463);
    const auto *lg_464 = buffer.data(lg + 464);
    const auto *lg_465 = buffer.data(lg + 465);
    const auto *lg_466 = buffer.data(lg + 466);
    const auto *lg_467 = buffer.data(lg + 467);
    const auto *lg_468 = buffer.data(lg + 468);
    const auto *lg_469 = buffer.data(lg + 469);
    const auto *lg_470 = buffer.data(lg + 470);
    const auto *lg_471 = buffer.data(lg + 471);
    const auto *lg_472 = buffer.data(lg + 472);
    const auto *lg_473 = buffer.data(lg + 473);
    const auto *lg_474 = buffer.data(lg + 474);
    const auto *lg_475 = buffer.data(lg + 475);
    const auto *lg_476 = buffer.data(lg + 476);
    const auto *lg_477 = buffer.data(lg + 477);
    const auto *lg_478 = buffer.data(lg + 478);
    const auto *lg_479 = buffer.data(lg + 479);
    const auto *lg_480 = buffer.data(lg + 480);
    const auto *lg_481 = buffer.data(lg + 481);
    const auto *lg_482 = buffer.data(lg + 482);
    const auto *lg_483 = buffer.data(lg + 483);
    const auto *lg_484 = buffer.data(lg + 484);
    const auto *lg_485 = buffer.data(lg + 485);
    const auto *lg_486 = buffer.data(lg + 486);
    const auto *lg_487 = buffer.data(lg + 487);
    const auto *lg_488 = buffer.data(lg + 488);
    const auto *lg_489 = buffer.data(lg + 489);
    const auto *lg_490 = buffer.data(lg + 490);
    const auto *lg_491 = buffer.data(lg + 491);
    const auto *lg_492 = buffer.data(lg + 492);
    const auto *lg_493 = buffer.data(lg + 493);
    const auto *lg_494 = buffer.data(lg + 494);
    const auto *lg_495 = buffer.data(lg + 495);
    const auto *lg_496 = buffer.data(lg + 496);
    const auto *lg_497 = buffer.data(lg + 497);
    const auto *lg_498 = buffer.data(lg + 498);
    const auto *lg_499 = buffer.data(lg + 499);
    const auto *lg_500 = buffer.data(lg + 500);
    const auto *lg_501 = buffer.data(lg + 501);
    const auto *lg_502 = buffer.data(lg + 502);
    const auto *lg_503 = buffer.data(lg + 503);
    const auto *lg_504 = buffer.data(lg + 504);
    const auto *lg_505 = buffer.data(lg + 505);
    const auto *lg_506 = buffer.data(lg + 506);
    const auto *lg_507 = buffer.data(lg + 507);
    const auto *lg_508 = buffer.data(lg + 508);
    const auto *lg_509 = buffer.data(lg + 509);
    const auto *lg_510 = buffer.data(lg + 510);
    const auto *lg_511 = buffer.data(lg + 511);
    const auto *lg_512 = buffer.data(lg + 512);
    const auto *lg_513 = buffer.data(lg + 513);
    const auto *lg_514 = buffer.data(lg + 514);
    const auto *lg_515 = buffer.data(lg + 515);
    const auto *lg_516 = buffer.data(lg + 516);
    const auto *lg_517 = buffer.data(lg + 517);
    const auto *lg_518 = buffer.data(lg + 518);
    const auto *lg_519 = buffer.data(lg + 519);
    const auto *lg_520 = buffer.data(lg + 520);
    const auto *lg_521 = buffer.data(lg + 521);
    const auto *lg_522 = buffer.data(lg + 522);
    const auto *lg_523 = buffer.data(lg + 523);
    const auto *lg_524 = buffer.data(lg + 524);
    const auto *lg_525 = buffer.data(lg + 525);
    const auto *lg_526 = buffer.data(lg + 526);
    const auto *lg_527 = buffer.data(lg + 527);
    const auto *lg_528 = buffer.data(lg + 528);
    const auto *lg_529 = buffer.data(lg + 529);
    const auto *lg_530 = buffer.data(lg + 530);
    const auto *lg_531 = buffer.data(lg + 531);
    const auto *lg_532 = buffer.data(lg + 532);
    const auto *lg_533 = buffer.data(lg + 533);
    const auto *lg_534 = buffer.data(lg + 534);
    const auto *lg_535 = buffer.data(lg + 535);
    const auto *lg_536 = buffer.data(lg + 536);
    const auto *lg_537 = buffer.data(lg + 537);
    const auto *lg_538 = buffer.data(lg + 538);
    const auto *lg_539 = buffer.data(lg + 539);
    const auto *lg_540 = buffer.data(lg + 540);
    const auto *lg_541 = buffer.data(lg + 541);
    const auto *lg_542 = buffer.data(lg + 542);
    const auto *lg_543 = buffer.data(lg + 543);
    const auto *lg_544 = buffer.data(lg + 544);
    const auto *lg_545 = buffer.data(lg + 545);
    const auto *lg_546 = buffer.data(lg + 546);
    const auto *lg_547 = buffer.data(lg + 547);
    const auto *lg_548 = buffer.data(lg + 548);
    const auto *lg_549 = buffer.data(lg + 549);
    const auto *lg_550 = buffer.data(lg + 550);
    const auto *lg_551 = buffer.data(lg + 551);
    const auto *lg_552 = buffer.data(lg + 552);
    const auto *lg_553 = buffer.data(lg + 553);
    const auto *lg_554 = buffer.data(lg + 554);
    const auto *lg_555 = buffer.data(lg + 555);
    const auto *lg_556 = buffer.data(lg + 556);
    const auto *lg_557 = buffer.data(lg + 557);
    const auto *lg_558 = buffer.data(lg + 558);
    const auto *lg_559 = buffer.data(lg + 559);
    const auto *lg_560 = buffer.data(lg + 560);
    const auto *lg_561 = buffer.data(lg + 561);
    const auto *lg_562 = buffer.data(lg + 562);
    const auto *lg_563 = buffer.data(lg + 563);
    const auto *lg_564 = buffer.data(lg + 564);
    const auto *lg_565 = buffer.data(lg + 565);
    const auto *lg_566 = buffer.data(lg + 566);
    const auto *lg_567 = buffer.data(lg + 567);
    const auto *lg_568 = buffer.data(lg + 568);
    const auto *lg_569 = buffer.data(lg + 569);
    const auto *lg_570 = buffer.data(lg + 570);
    const auto *lg_571 = buffer.data(lg + 571);
    const auto *lg_572 = buffer.data(lg + 572);
    const auto *lg_573 = buffer.data(lg + 573);
    const auto *lg_574 = buffer.data(lg + 574);
    const auto *lg_575 = buffer.data(lg + 575);
    const auto *lg_576 = buffer.data(lg + 576);
    const auto *lg_577 = buffer.data(lg + 577);
    const auto *lg_578 = buffer.data(lg + 578);
    const auto *lg_579 = buffer.data(lg + 579);
    const auto *lg_580 = buffer.data(lg + 580);
    const auto *lg_581 = buffer.data(lg + 581);
    const auto *lg_582 = buffer.data(lg + 582);
    const auto *lg_583 = buffer.data(lg + 583);
    const auto *lg_584 = buffer.data(lg + 584);
    const auto *lg_585 = buffer.data(lg + 585);
    const auto *lg_586 = buffer.data(lg + 586);
    const auto *lg_587 = buffer.data(lg + 587);
    const auto *lg_588 = buffer.data(lg + 588);
    const auto *lg_589 = buffer.data(lg + 589);
    const auto *lg_590 = buffer.data(lg + 590);
    const auto *lg_591 = buffer.data(lg + 591);
    const auto *lg_592 = buffer.data(lg + 592);
    const auto *lg_593 = buffer.data(lg + 593);
    const auto *lg_594 = buffer.data(lg + 594);
    const auto *lg_595 = buffer.data(lg + 595);
    const auto *lg_596 = buffer.data(lg + 596);
    const auto *lg_597 = buffer.data(lg + 597);
    const auto *lg_598 = buffer.data(lg + 598);
    const auto *lg_599 = buffer.data(lg + 599);
    const auto *lg_600 = buffer.data(lg + 600);
    const auto *lg_601 = buffer.data(lg + 601);
    const auto *lg_602 = buffer.data(lg + 602);
    const auto *lg_603 = buffer.data(lg + 603);
    const auto *lg_604 = buffer.data(lg + 604);
    const auto *lg_605 = buffer.data(lg + 605);
    const auto *lg_606 = buffer.data(lg + 606);
    const auto *lg_607 = buffer.data(lg + 607);
    const auto *lg_608 = buffer.data(lg + 608);
    const auto *lg_609 = buffer.data(lg + 609);
    const auto *lg_610 = buffer.data(lg + 610);
    const auto *lg_611 = buffer.data(lg + 611);
    const auto *lg_612 = buffer.data(lg + 612);
    const auto *lg_613 = buffer.data(lg + 613);
    const auto *lg_614 = buffer.data(lg + 614);
    const auto *lg_615 = buffer.data(lg + 615);
    const auto *lg_616 = buffer.data(lg + 616);
    const auto *lg_617 = buffer.data(lg + 617);
    const auto *lg_618 = buffer.data(lg + 618);
    const auto *lg_619 = buffer.data(lg + 619);
    const auto *lg_620 = buffer.data(lg + 620);
    const auto *lg_621 = buffer.data(lg + 621);
    const auto *lg_622 = buffer.data(lg + 622);
    const auto *lg_623 = buffer.data(lg + 623);
    const auto *lg_624 = buffer.data(lg + 624);
    const auto *lg_625 = buffer.data(lg + 625);
    const auto *lg_626 = buffer.data(lg + 626);
    const auto *lg_627 = buffer.data(lg + 627);
    const auto *lg_628 = buffer.data(lg + 628);
    const auto *lg_629 = buffer.data(lg + 629);
    const auto *lg_630 = buffer.data(lg + 630);
    const auto *lg_631 = buffer.data(lg + 631);
    const auto *lg_632 = buffer.data(lg + 632);
    const auto *lg_633 = buffer.data(lg + 633);
    const auto *lg_634 = buffer.data(lg + 634);
    const auto *lg_635 = buffer.data(lg + 635);
    const auto *lg_636 = buffer.data(lg + 636);
    const auto *lg_637 = buffer.data(lg + 637);
    const auto *lg_638 = buffer.data(lg + 638);
    const auto *lg_639 = buffer.data(lg + 639);
    const auto *lg_640 = buffer.data(lg + 640);
    const auto *lg_641 = buffer.data(lg + 641);
    const auto *lg_642 = buffer.data(lg + 642);
    const auto *lg_643 = buffer.data(lg + 643);
    const auto *lg_644 = buffer.data(lg + 644);
    const auto *lg_645 = buffer.data(lg + 645);
    const auto *lg_646 = buffer.data(lg + 646);
    const auto *lg_647 = buffer.data(lg + 647);
    const auto *lg_648 = buffer.data(lg + 648);
    const auto *lg_649 = buffer.data(lg + 649);
    const auto *lg_650 = buffer.data(lg + 650);
    const auto *lg_651 = buffer.data(lg + 651);
    const auto *lg_652 = buffer.data(lg + 652);
    const auto *lg_653 = buffer.data(lg + 653);
    const auto *lg_654 = buffer.data(lg + 654);
    const auto *lg_655 = buffer.data(lg + 655);
    const auto *lg_656 = buffer.data(lg + 656);
    const auto *lg_657 = buffer.data(lg + 657);
    const auto *lg_658 = buffer.data(lg + 658);
    const auto *lg_659 = buffer.data(lg + 659);
    const auto *lg_660 = buffer.data(lg + 660);
    const auto *lg_661 = buffer.data(lg + 661);
    const auto *lg_662 = buffer.data(lg + 662);
    const auto *lg_663 = buffer.data(lg + 663);
    const auto *lg_664 = buffer.data(lg + 664);
    const auto *lg_665 = buffer.data(lg + 665);
    const auto *lg_666 = buffer.data(lg + 666);
    const auto *lg_667 = buffer.data(lg + 667);
    const auto *lg_668 = buffer.data(lg + 668);
    const auto *lg_669 = buffer.data(lg + 669);
    const auto *lg_670 = buffer.data(lg + 670);
    const auto *lg_671 = buffer.data(lg + 671);
    const auto *lg_672 = buffer.data(lg + 672);
    const auto *lg_673 = buffer.data(lg + 673);
    const auto *lg_674 = buffer.data(lg + 674);

#pragma omp simd aligned(lg_16, lg_21, lg_91, lg_96, lg_226, lg_231, lg_421, \
                         lg_426 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * lg_16[k]
                 - f_0 * lg_21[k]
                 - f_1 * lg_91[k]
                 + f_1 * lg_96[k]
                 + f_1 * lg_226[k]
                 - f_1 * lg_231[k]
                 - f_0 * lg_421[k]
                 + f_0 * lg_426[k];
    }

#pragma omp simd aligned(lg_19, lg_26, lg_94, lg_101, lg_229, lg_236, lg_424, \
                         lg_431 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_1[k] = f_2 * lg_19[k]
                 - f_3 * lg_26[k]
                 - f_4 * lg_94[k]
                 + f_5 * lg_101[k]
                 + f_4 * lg_229[k]
                 - f_5 * lg_236[k]
                 - f_2 * lg_424[k]
                 + f_3 * lg_431[k];
    }

#pragma omp simd aligned(lg_16, lg_21, lg_23, lg_91, lg_96, lg_98, lg_226, lg_231, lg_233, \
                         lg_421, lg_426, lg_428 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_2[k] = -f_6 * lg_16[k]
                 - f_6 * lg_21[k]
                 + f_7 * lg_23[k]
                 + f_8 * lg_91[k]
                 + f_8 * lg_96[k]
                 - f_9 * lg_98[k]
                 - f_8 * lg_226[k]
                 - f_8 * lg_231[k]
                 + f_9 * lg_233[k]
                 + f_6 * lg_421[k]
                 + f_6 * lg_426[k]
                 - f_7 * lg_428[k];
    }

#pragma omp simd aligned(lg_19, lg_26, lg_28, lg_94, lg_101, lg_103, lg_229, lg_236, lg_238, \
                         lg_424, lg_431, lg_433 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = -f_10 * lg_19[k]
                 - f_10 * lg_26[k]
                 + f_11 * lg_28[k]
                 + f_12 * lg_94[k]
                 + f_12 * lg_101[k]
                 - f_13 * lg_103[k]
                 - f_12 * lg_229[k]
                 - f_12 * lg_236[k]
                 + f_13 * lg_238[k]
                 + f_10 * lg_424[k]
                 + f_10 * lg_431[k]
                 - f_11 * lg_433[k];
    }

#pragma omp simd aligned(lg_15, lg_18, lg_20, lg_25, lg_27, lg_29, lg_90, lg_93, lg_95, \
                         lg_100, lg_102, lg_104, lg_225, lg_228, lg_230, lg_235, lg_237, \
                         lg_239, lg_420, lg_423, lg_425, lg_430, lg_432, \
                         lg_434 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = f_14 * lg_15[k]
                 + f_15 * lg_18[k]
                 - f_16 * lg_20[k]
                 + f_14 * lg_25[k]
                 - f_16 * lg_27[k]
                 + f_17 * lg_29[k]
                 - f_18 * lg_90[k]
                 - f_19 * lg_93[k]
                 + f_20 * lg_95[k]
                 - f_18 * lg_100[k]
                 + f_20 * lg_102[k]
                 - f_21 * lg_104[k]
                 + f_18 * lg_225[k]
                 + f_19 * lg_228[k]
                 - f_20 * lg_230[k]
                 + f_18 * lg_235[k]
                 - f_20 * lg_237[k]
                 + f_21 * lg_239[k]
                 - f_14 * lg_420[k]
                 - f_15 * lg_423[k]
                 + f_16 * lg_425[k]
                 - f_14 * lg_430[k]
                 + f_16 * lg_432[k]
                 - f_17 * lg_434[k];
    }

#pragma omp simd aligned(lg_17, lg_22, lg_24, lg_92, lg_97, lg_99, lg_227, lg_232, lg_234, \
                         lg_422, lg_427, lg_429 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = -f_10 * lg_17[k]
                 - f_10 * lg_22[k]
                 + f_11 * lg_24[k]
                 + f_12 * lg_92[k]
                 + f_12 * lg_97[k]
                 - f_13 * lg_99[k]
                 - f_12 * lg_227[k]
                 - f_12 * lg_232[k]
                 + f_13 * lg_234[k]
                 + f_10 * lg_422[k]
                 + f_10 * lg_427[k]
                 - f_11 * lg_429[k];
    }

#pragma omp simd aligned(lg_15, lg_20, lg_25, lg_27, lg_90, lg_95, lg_100, lg_102, lg_225, \
                         lg_230, lg_235, lg_237, lg_420, lg_425, lg_430, \
                         lg_432 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = -f_22 * lg_15[k]
                 + f_23 * lg_20[k]
                 + f_22 * lg_25[k]
                 - f_23 * lg_27[k]
                 + f_24 * lg_90[k]
                 - f_25 * lg_95[k]
                 - f_24 * lg_100[k]
                 + f_25 * lg_102[k]
                 - f_24 * lg_225[k]
                 + f_25 * lg_230[k]
                 + f_24 * lg_235[k]
                 - f_25 * lg_237[k]
                 + f_22 * lg_420[k]
                 - f_23 * lg_425[k]
                 - f_22 * lg_430[k]
                 + f_23 * lg_432[k];
    }

#pragma omp simd aligned(lg_17, lg_22, lg_92, lg_97, lg_227, lg_232, lg_422, \
                         lg_427 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = f_3 * lg_17[k]
                 - f_2 * lg_22[k]
                 - f_5 * lg_92[k]
                 + f_4 * lg_97[k]
                 + f_5 * lg_227[k]
                 - f_4 * lg_232[k]
                 - f_3 * lg_422[k]
                 + f_2 * lg_427[k];
    }

#pragma omp simd aligned(lg_15, lg_18, lg_25, lg_90, lg_93, lg_100, lg_225, lg_228, lg_235, \
                         lg_420, lg_423, lg_430 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_8[k] = f_26 * lg_15[k]
                 - f_27 * lg_18[k]
                 + f_26 * lg_25[k]
                 - f_28 * lg_90[k]
                 + f_29 * lg_93[k]
                 - f_28 * lg_100[k]
                 + f_28 * lg_225[k]
                 - f_29 * lg_228[k]
                 + f_28 * lg_235[k]
                 - f_26 * lg_420[k]
                 + f_27 * lg_423[k]
                 - f_26 * lg_430[k];
    }

#pragma omp simd aligned(lg_61, lg_66, lg_166, lg_171, lg_331, lg_336, lg_556, \
                         lg_561 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_9[k] = f_30 * lg_61[k]
                 - f_30 * lg_66[k]
                 - f_31 * lg_166[k]
                 + f_31 * lg_171[k]
                 + f_29 * lg_331[k]
                 - f_29 * lg_336[k]
                 - f_32 * lg_556[k]
                 + f_32 * lg_561[k];
    }

#pragma omp simd aligned(lg_64, lg_71, lg_169, lg_176, lg_334, lg_341, lg_559, \
                         lg_566 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = f_33 * lg_64[k]
                  - f_34 * lg_71[k]
                  - f_35 * lg_169[k]
                  + f_36 * lg_176[k]
                  + f_37 * lg_334[k]
                  - f_33 * lg_341[k]
                  - f_38 * lg_559[k]
                  + f_39 * lg_566[k];
    }

#pragma omp simd aligned(lg_61, lg_66, lg_68, lg_166, lg_171, lg_173, lg_331, lg_336, lg_338, \
                         lg_556, lg_561, lg_563 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_11[k] = -f_24 * lg_61[k]
                  - f_24 * lg_66[k]
                  + f_25 * lg_68[k]
                  + f_40 * lg_166[k]
                  + f_40 * lg_171[k]
                  - f_41 * lg_173[k]
                  - f_42 * lg_331[k]
                  - f_42 * lg_336[k]
                  + f_43 * lg_338[k]
                  + f_22 * lg_556[k]
                  + f_22 * lg_561[k]
                  - f_23 * lg_563[k];
    }

#pragma omp simd aligned(lg_64, lg_71, lg_73, lg_169, lg_176, lg_178, lg_334, lg_341, lg_343, \
                         lg_559, lg_566, lg_568 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_12[k] = -f_44 * lg_64[k]
                  - f_44 * lg_71[k]
                  + f_45 * lg_73[k]
                  + f_46 * lg_169[k]
                  + f_46 * lg_176[k]
                  - f_47 * lg_178[k]
                  - f_48 * lg_334[k]
                  - f_48 * lg_341[k]
                  + f_49 * lg_343[k]
                  + f_50 * lg_559[k]
                  + f_50 * lg_566[k]
                  - f_51 * lg_568[k];
    }

#pragma omp simd aligned(lg_60, lg_63, lg_65, lg_70, lg_72, lg_74, lg_165, lg_168, lg_170, \
                         lg_175, lg_177, lg_179, lg_330, lg_333, lg_335, lg_340, lg_342, \
                         lg_344, lg_555, lg_558, lg_560, lg_565, lg_567, \
                         lg_569 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = f_52 * lg_60[k]
                  + f_18 * lg_63[k]
                  - f_53 * lg_65[k]
                  + f_52 * lg_70[k]
                  - f_53 * lg_72[k]
                  + f_54 * lg_74[k]
                  - f_55 * lg_165[k]
                  - f_56 * lg_168[k]
                  + f_57 * lg_170[k]
                  - f_55 * lg_175[k]
                  + f_57 * lg_177[k]
                  - f_58 * lg_179[k]
                  + f_59 * lg_330[k]
                  + f_60 * lg_333[k]
                  - f_61 * lg_335[k]
                  + f_59 * lg_340[k]
                  - f_61 * lg_342[k]
                  + f_53 * lg_344[k]
                  - f_62 * lg_555[k]
                  - f_14 * lg_558[k]
                  + f_63 * lg_560[k]
                  - f_62 * lg_565[k]
                  + f_63 * lg_567[k]
                  - f_64 * lg_569[k];
    }

#pragma omp simd aligned(lg_62, lg_67, lg_69, lg_167, lg_172, lg_174, lg_332, lg_337, lg_339, \
                         lg_557, lg_562, lg_564 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = -f_44 * lg_62[k]
                  - f_44 * lg_67[k]
                  + f_45 * lg_69[k]
                  + f_46 * lg_167[k]
                  + f_46 * lg_172[k]
                  - f_47 * lg_174[k]
                  - f_48 * lg_332[k]
                  - f_48 * lg_337[k]
                  + f_49 * lg_339[k]
                  + f_50 * lg_557[k]
                  + f_50 * lg_562[k]
                  - f_51 * lg_564[k];
    }

#pragma omp simd aligned(lg_60, lg_65, lg_70, lg_72, lg_165, lg_170, lg_175, lg_177, lg_330, \
                         lg_335, lg_340, lg_342, lg_555, lg_560, lg_565, \
                         lg_567 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_15[k] = -f_65 * lg_60[k]
                  + f_42 * lg_65[k]
                  + f_65 * lg_70[k]
                  - f_42 * lg_72[k]
                  + f_66 * lg_165[k]
                  - f_67 * lg_170[k]
                  - f_66 * lg_175[k]
                  + f_67 * lg_177[k]
                  - f_68 * lg_330[k]
                  + f_69 * lg_335[k]
                  + f_68 * lg_340[k]
                  - f_69 * lg_342[k]
                  + f_70 * lg_555[k]
                  - f_71 * lg_560[k]
                  - f_70 * lg_565[k]
                  + f_71 * lg_567[k];
    }

#pragma omp simd aligned(lg_62, lg_67, lg_167, lg_172, lg_332, lg_337, lg_557, \
                         lg_562 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_16[k] = f_34 * lg_62[k]
                  - f_33 * lg_67[k]
                  - f_36 * lg_167[k]
                  + f_35 * lg_172[k]
                  + f_33 * lg_332[k]
                  - f_37 * lg_337[k]
                  - f_39 * lg_557[k]
                  + f_38 * lg_562[k];
    }

#pragma omp simd aligned(lg_60, lg_63, lg_70, lg_165, lg_168, lg_175, lg_330, lg_333, lg_340, \
                         lg_555, lg_558, lg_565 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = f_72 * lg_60[k]
                  - f_73 * lg_63[k]
                  + f_72 * lg_70[k]
                  - f_74 * lg_165[k]
                  + f_75 * lg_168[k]
                  - f_74 * lg_175[k]
                  + f_76 * lg_330[k]
                  - f_77 * lg_333[k]
                  + f_76 * lg_340[k]
                  - f_78 * lg_555[k]
                  + f_79 * lg_558[k]
                  - f_78 * lg_565[k];
    }

#pragma omp simd aligned(lg_16, lg_21, lg_91, lg_96, lg_121, lg_126, lg_226, lg_231, lg_256, \
                         lg_261, lg_421, lg_426, lg_451, lg_456 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = -f_80 * lg_16[k]
                  + f_80 * lg_21[k]
                  + f_81 * lg_91[k]
                  - f_81 * lg_96[k]
                  + f_82 * lg_121[k]
                  - f_82 * lg_126[k]
                  + f_81 * lg_226[k]
                  - f_81 * lg_231[k]
                  - f_83 * lg_256[k]
                  + f_83 * lg_261[k]
                  - f_80 * lg_421[k]
                  + f_80 * lg_426[k]
                  + f_82 * lg_451[k]
                  - f_82 * lg_456[k];
    }

#pragma omp simd aligned(lg_19, lg_26, lg_94, lg_101, lg_124, lg_131, lg_229, lg_236, lg_259, \
                         lg_266, lg_424, lg_431, lg_454, lg_461 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = -f_84 * lg_19[k]
                  + f_85 * lg_26[k]
                  + f_86 * lg_94[k]
                  - f_87 * lg_101[k]
                  + f_88 * lg_124[k]
                  - f_89 * lg_131[k]
                  + f_86 * lg_229[k]
                  - f_87 * lg_236[k]
                  - f_90 * lg_259[k]
                  + f_91 * lg_266[k]
                  - f_84 * lg_424[k]
                  + f_85 * lg_431[k]
                  + f_88 * lg_454[k]
                  - f_89 * lg_461[k];
    }

#pragma omp simd aligned(lg_16, lg_21, lg_23, lg_91, lg_96, lg_98, lg_121, lg_126, lg_128, \
                         lg_226, lg_231, lg_233, lg_256, lg_261, lg_263, lg_421, lg_426, \
                         lg_428, lg_451, lg_456, lg_458 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = f_92 * lg_16[k]
                  + f_92 * lg_21[k]
                  - f_93 * lg_23[k]
                  - f_94 * lg_91[k]
                  - f_94 * lg_96[k]
                  + f_95 * lg_98[k]
                  - f_95 * lg_121[k]
                  - f_95 * lg_126[k]
                  + f_96 * lg_128[k]
                  - f_94 * lg_226[k]
                  - f_94 * lg_231[k]
                  + f_95 * lg_233[k]
                  + f_97 * lg_256[k]
                  + f_97 * lg_261[k]
                  - f_98 * lg_263[k]
                  + f_92 * lg_421[k]
                  + f_92 * lg_426[k]
                  - f_93 * lg_428[k]
                  - f_95 * lg_451[k]
                  - f_95 * lg_456[k]
                  + f_96 * lg_458[k];
    }

#pragma omp simd aligned(lg_19, lg_26, lg_28, lg_94, lg_101, lg_103, lg_124, lg_131, lg_133, \
                         lg_229, lg_236, lg_238, lg_259, lg_266, lg_268, lg_424, lg_431, \
                         lg_433, lg_454, lg_461, lg_463 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = f_99 * lg_19[k]
                  + f_99 * lg_26[k]
                  - f_100 * lg_28[k]
                  - f_101 * lg_94[k]
                  - f_101 * lg_101[k]
                  + f_102 * lg_103[k]
                  - f_103 * lg_124[k]
                  - f_103 * lg_131[k]
                  + f_104 * lg_133[k]
                  - f_101 * lg_229[k]
                  - f_101 * lg_236[k]
                  + f_102 * lg_238[k]
                  + f_105 * lg_259[k]
                  + f_105 * lg_266[k]
                  - f_106 * lg_268[k]
                  + f_99 * lg_424[k]
                  + f_99 * lg_431[k]
                  - f_100 * lg_433[k]
                  - f_103 * lg_454[k]
                  - f_103 * lg_461[k]
                  + f_104 * lg_463[k];
    }

#pragma omp simd aligned(lg_15, lg_18, lg_20, lg_25, lg_27, lg_29, lg_90, lg_93, lg_95, \
                         lg_100, lg_102, lg_104, lg_120, lg_123, lg_125, lg_130, lg_132, \
                         lg_134, lg_225, lg_228, lg_230, lg_235, lg_237, lg_239, lg_255, \
                         lg_258, lg_260, lg_265, lg_267, lg_269, lg_420, lg_423, lg_425, \
                         lg_430, lg_432, lg_434, lg_450, lg_453, lg_455, lg_460, lg_462, \
                         lg_464 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = -f_107 * lg_15[k]
                  - f_108 * lg_18[k]
                  + f_109 * lg_20[k]
                  - f_107 * lg_25[k]
                  + f_109 * lg_27[k]
                  - f_110 * lg_29[k]
                  + f_111 * lg_90[k]
                  + f_112 * lg_93[k]
                  - f_113 * lg_95[k]
                  + f_111 * lg_100[k]
                  - f_113 * lg_102[k]
                  + f_114 * lg_104[k]
                  + f_115 * lg_120[k]
                  + f_116 * lg_123[k]
                  - f_117 * lg_125[k]
                  + f_115 * lg_130[k]
                  - f_117 * lg_132[k]
                  + f_118 * lg_134[k]
                  + f_111 * lg_225[k]
                  + f_112 * lg_228[k]
                  - f_113 * lg_230[k]
                  + f_111 * lg_235[k]
                  - f_113 * lg_237[k]
                  + f_114 * lg_239[k]
                  - f_119 * lg_255[k]
                  - f_120 * lg_258[k]
                  + f_121 * lg_260[k]
                  - f_119 * lg_265[k]
                  + f_121 * lg_267[k]
                  - f_122 * lg_269[k]
                  - f_107 * lg_420[k]
                  - f_108 * lg_423[k]
                  + f_109 * lg_425[k]
                  - f_107 * lg_430[k]
                  + f_109 * lg_432[k]
                  - f_110 * lg_434[k]
                  + f_115 * lg_450[k]
                  + f_116 * lg_453[k]
                  - f_117 * lg_455[k]
                  + f_115 * lg_460[k]
                  - f_117 * lg_462[k]
                  + f_118 * lg_464[k];
    }

#pragma omp simd aligned(lg_17, lg_22, lg_24, lg_92, lg_97, lg_99, lg_122, lg_127, lg_129, \
                         lg_227, lg_232, lg_234, lg_257, lg_262, lg_264, lg_422, lg_427, \
                         lg_429, lg_452, lg_457, lg_459 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = f_99 * lg_17[k]
                  + f_99 * lg_22[k]
                  - f_100 * lg_24[k]
                  - f_101 * lg_92[k]
                  - f_101 * lg_97[k]
                  + f_102 * lg_99[k]
                  - f_103 * lg_122[k]
                  - f_103 * lg_127[k]
                  + f_104 * lg_129[k]
                  - f_101 * lg_227[k]
                  - f_101 * lg_232[k]
                  + f_102 * lg_234[k]
                  + f_105 * lg_257[k]
                  + f_105 * lg_262[k]
                  - f_106 * lg_264[k]
                  + f_99 * lg_422[k]
                  + f_99 * lg_427[k]
                  - f_100 * lg_429[k]
                  - f_103 * lg_452[k]
                  - f_103 * lg_457[k]
                  + f_104 * lg_459[k];
    }

#pragma omp simd aligned(lg_15, lg_20, lg_25, lg_27, lg_90, lg_95, lg_100, lg_102, lg_120, \
                         lg_125, lg_130, lg_132, lg_225, lg_230, lg_235, lg_237, lg_255, \
                         lg_260, lg_265, lg_267, lg_420, lg_425, lg_430, lg_432, lg_450, \
                         lg_455, lg_460, lg_462 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = f_123 * lg_15[k]
                  - f_124 * lg_20[k]
                  - f_123 * lg_25[k]
                  + f_124 * lg_27[k]
                  - f_125 * lg_90[k]
                  + f_126 * lg_95[k]
                  + f_125 * lg_100[k]
                  - f_126 * lg_102[k]
                  - f_126 * lg_120[k]
                  + f_127 * lg_125[k]
                  + f_126 * lg_130[k]
                  - f_127 * lg_132[k]
                  - f_125 * lg_225[k]
                  + f_126 * lg_230[k]
                  + f_125 * lg_235[k]
                  - f_126 * lg_237[k]
                  + f_128 * lg_255[k]
                  - f_129 * lg_260[k]
                  - f_128 * lg_265[k]
                  + f_129 * lg_267[k]
                  + f_123 * lg_420[k]
                  - f_124 * lg_425[k]
                  - f_123 * lg_430[k]
                  + f_124 * lg_432[k]
                  - f_126 * lg_450[k]
                  + f_127 * lg_455[k]
                  + f_126 * lg_460[k]
                  - f_127 * lg_462[k];
    }

#pragma omp simd aligned(lg_17, lg_22, lg_92, lg_97, lg_122, lg_127, lg_227, lg_232, lg_257, \
                         lg_262, lg_422, lg_427, lg_452, lg_457 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = -f_85 * lg_17[k]
                  + f_84 * lg_22[k]
                  + f_87 * lg_92[k]
                  - f_86 * lg_97[k]
                  + f_89 * lg_122[k]
                  - f_88 * lg_127[k]
                  + f_87 * lg_227[k]
                  - f_86 * lg_232[k]
                  - f_91 * lg_257[k]
                  + f_90 * lg_262[k]
                  - f_85 * lg_422[k]
                  + f_84 * lg_427[k]
                  + f_89 * lg_452[k]
                  - f_88 * lg_457[k];
    }

#pragma omp simd aligned(lg_15, lg_18, lg_25, lg_90, lg_93, lg_100, lg_120, lg_123, lg_130, \
                         lg_225, lg_228, lg_235, lg_255, lg_258, lg_265, lg_420, lg_423, \
                         lg_430, lg_450, lg_453, lg_460 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = -f_130 * lg_15[k]
                  + f_131 * lg_18[k]
                  - f_130 * lg_25[k]
                  + f_132 * lg_90[k]
                  - f_133 * lg_93[k]
                  + f_132 * lg_100[k]
                  + f_133 * lg_120[k]
                  - f_134 * lg_123[k]
                  + f_133 * lg_130[k]
                  + f_132 * lg_225[k]
                  - f_133 * lg_228[k]
                  + f_132 * lg_235[k]
                  - f_135 * lg_255[k]
                  + f_136 * lg_258[k]
                  - f_135 * lg_265[k]
                  - f_130 * lg_420[k]
                  + f_131 * lg_423[k]
                  - f_130 * lg_430[k]
                  + f_133 * lg_450[k]
                  - f_134 * lg_453[k]
                  + f_133 * lg_460[k];
    }

#pragma omp simd aligned(lg_61, lg_66, lg_166, lg_171, lg_196, lg_201, lg_331, lg_336, lg_361, \
                         lg_366, lg_556, lg_561, lg_586, lg_591 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = -f_137 * lg_61[k]
                  + f_137 * lg_66[k]
                  + f_137 * lg_166[k]
                  - f_137 * lg_171[k]
                  + f_138 * lg_196[k]
                  - f_138 * lg_201[k]
                  + f_139 * lg_331[k]
                  - f_139 * lg_336[k]
                  - f_140 * lg_361[k]
                  + f_140 * lg_366[k]
                  - f_141 * lg_556[k]
                  + f_141 * lg_561[k]
                  + f_21 * lg_586[k]
                  - f_21 * lg_591[k];
    }

#pragma omp simd aligned(lg_64, lg_71, lg_169, lg_176, lg_199, lg_206, lg_334, lg_341, lg_364, \
                         lg_371, lg_559, lg_566, lg_589, lg_596 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = -f_142 * lg_64[k]
                  + f_143 * lg_71[k]
                  + f_142 * lg_169[k]
                  - f_143 * lg_176[k]
                  + f_144 * lg_199[k]
                  - f_145 * lg_206[k]
                  + f_146 * lg_334[k]
                  - f_147 * lg_341[k]
                  - f_148 * lg_364[k]
                  + f_149 * lg_371[k]
                  - f_150 * lg_559[k]
                  + f_151 * lg_566[k]
                  + f_152 * lg_589[k]
                  - f_153 * lg_596[k];
    }

#pragma omp simd aligned(lg_61, lg_66, lg_68, lg_166, lg_171, lg_173, lg_196, lg_201, lg_203, \
                         lg_331, lg_336, lg_338, lg_361, lg_366, lg_368, lg_556, lg_561, \
                         lg_563, lg_586, lg_591, lg_593 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = f_154 * lg_61[k]
                  + f_154 * lg_66[k]
                  - f_155 * lg_68[k]
                  - f_154 * lg_166[k]
                  - f_154 * lg_171[k]
                  + f_155 * lg_173[k]
                  - f_156 * lg_196[k]
                  - f_156 * lg_201[k]
                  + f_157 * lg_203[k]
                  - f_158 * lg_331[k]
                  - f_158 * lg_336[k]
                  + f_159 * lg_338[k]
                  + f_160 * lg_361[k]
                  + f_160 * lg_366[k]
                  - f_161 * lg_368[k]
                  + f_162 * lg_556[k]
                  + f_162 * lg_561[k]
                  - f_163 * lg_563[k]
                  - f_164 * lg_586[k]
                  - f_164 * lg_591[k]
                  + f_165 * lg_593[k];
    }

#pragma omp simd aligned(lg_64, lg_71, lg_73, lg_169, lg_176, lg_178, lg_199, lg_206, lg_208, \
                         lg_334, lg_341, lg_343, lg_364, lg_371, lg_373, lg_559, lg_566, \
                         lg_568, lg_589, lg_596, lg_598 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = f_166 * lg_64[k]
                  + f_166 * lg_71[k]
                  - f_167 * lg_73[k]
                  - f_166 * lg_169[k]
                  - f_166 * lg_176[k]
                  + f_167 * lg_178[k]
                  - f_168 * lg_199[k]
                  - f_168 * lg_206[k]
                  + f_169 * lg_208[k]
                  - f_170 * lg_334[k]
                  - f_170 * lg_341[k]
                  + f_171 * lg_343[k]
                  + f_172 * lg_364[k]
                  + f_172 * lg_371[k]
                  - f_173 * lg_373[k]
                  + f_174 * lg_559[k]
                  + f_174 * lg_566[k]
                  - f_175 * lg_568[k]
                  - f_176 * lg_589[k]
                  - f_176 * lg_596[k]
                  + f_177 * lg_598[k];
    }

#pragma omp simd aligned(lg_60, lg_63, lg_65, lg_70, lg_72, lg_74, lg_165, lg_168, lg_170, \
                         lg_175, lg_177, lg_179, lg_195, lg_198, lg_200, lg_205, lg_207, \
                         lg_209, lg_330, lg_333, lg_335, lg_340, lg_342, lg_344, lg_360, \
                         lg_363, lg_365, lg_370, lg_372, lg_374, lg_555, lg_558, lg_560, \
                         lg_565, lg_567, lg_569, lg_585, lg_588, lg_590, lg_595, lg_597, \
                         lg_599 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = -f_178 * lg_60[k]
                  - f_79 * lg_63[k]
                  + f_179 * lg_65[k]
                  - f_178 * lg_70[k]
                  + f_179 * lg_72[k]
                  - f_0 * lg_74[k]
                  + f_178 * lg_165[k]
                  + f_79 * lg_168[k]
                  - f_179 * lg_170[k]
                  + f_178 * lg_175[k]
                  - f_179 * lg_177[k]
                  + f_0 * lg_179[k]
                  + f_27 * lg_195[k]
                  + f_179 * lg_198[k]
                  - f_180 * lg_200[k]
                  + f_27 * lg_205[k]
                  - f_180 * lg_207[k]
                  + f_181 * lg_209[k]
                  + f_182 * lg_330[k]
                  + f_183 * lg_333[k]
                  - f_184 * lg_335[k]
                  + f_182 * lg_340[k]
                  - f_184 * lg_342[k]
                  + f_185 * lg_344[k]
                  - f_179 * lg_360[k]
                  - f_186 * lg_363[k]
                  + f_187 * lg_365[k]
                  - f_179 * lg_370[k]
                  + f_187 * lg_372[k]
                  - f_188 * lg_374[k]
                  - f_189 * lg_555[k]
                  - f_190 * lg_558[k]
                  + f_191 * lg_560[k]
                  - f_189 * lg_565[k]
                  + f_191 * lg_567[k]
                  - f_192 * lg_569[k]
                  + f_193 * lg_585[k]
                  + f_191 * lg_588[k]
                  - f_194 * lg_590[k]
                  + f_193 * lg_595[k]
                  - f_194 * lg_597[k]
                  + f_195 * lg_599[k];
    }

#pragma omp simd aligned(lg_62, lg_67, lg_69, lg_167, lg_172, lg_174, lg_197, lg_202, lg_204, \
                         lg_332, lg_337, lg_339, lg_362, lg_367, lg_369, lg_557, lg_562, \
                         lg_564, lg_587, lg_592, lg_594 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = f_166 * lg_62[k]
                  + f_166 * lg_67[k]
                  - f_167 * lg_69[k]
                  - f_166 * lg_167[k]
                  - f_166 * lg_172[k]
                  + f_167 * lg_174[k]
                  - f_168 * lg_197[k]
                  - f_168 * lg_202[k]
                  + f_169 * lg_204[k]
                  - f_170 * lg_332[k]
                  - f_170 * lg_337[k]
                  + f_171 * lg_339[k]
                  + f_172 * lg_362[k]
                  + f_172 * lg_367[k]
                  - f_173 * lg_369[k]
                  + f_174 * lg_557[k]
                  + f_174 * lg_562[k]
                  - f_175 * lg_564[k]
                  - f_176 * lg_587[k]
                  - f_176 * lg_592[k]
                  + f_177 * lg_594[k];
    }

#pragma omp simd aligned(lg_60, lg_65, lg_70, lg_72, lg_165, lg_170, lg_175, lg_177, lg_195, \
                         lg_200, lg_205, lg_207, lg_330, lg_335, lg_340, lg_342, lg_360, \
                         lg_365, lg_370, lg_372, lg_555, lg_560, lg_565, lg_567, lg_585, \
                         lg_590, lg_595, lg_597 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = f_196 * lg_60[k]
                  - f_197 * lg_65[k]
                  - f_196 * lg_70[k]
                  + f_197 * lg_72[k]
                  - f_196 * lg_165[k]
                  + f_197 * lg_170[k]
                  + f_196 * lg_175[k]
                  - f_197 * lg_177[k]
                  - f_198 * lg_195[k]
                  + f_199 * lg_200[k]
                  + f_198 * lg_205[k]
                  - f_199 * lg_207[k]
                  - f_200 * lg_330[k]
                  + f_201 * lg_335[k]
                  + f_200 * lg_340[k]
                  - f_201 * lg_342[k]
                  + f_156 * lg_360[k]
                  - f_157 * lg_365[k]
                  - f_156 * lg_370[k]
                  + f_157 * lg_372[k]
                  + f_202 * lg_555[k]
                  - f_203 * lg_560[k]
                  - f_202 * lg_565[k]
                  + f_203 * lg_567[k]
                  - f_204 * lg_585[k]
                  + f_205 * lg_590[k]
                  + f_204 * lg_595[k]
                  - f_205 * lg_597[k];
    }

#pragma omp simd aligned(lg_62, lg_67, lg_167, lg_172, lg_197, lg_202, lg_332, lg_337, lg_362, \
                         lg_367, lg_557, lg_562, lg_587, lg_592 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = -f_143 * lg_62[k]
                  + f_142 * lg_67[k]
                  + f_143 * lg_167[k]
                  - f_142 * lg_172[k]
                  + f_145 * lg_197[k]
                  - f_144 * lg_202[k]
                  + f_147 * lg_332[k]
                  - f_146 * lg_337[k]
                  - f_149 * lg_362[k]
                  + f_148 * lg_367[k]
                  - f_151 * lg_557[k]
                  + f_150 * lg_562[k]
                  + f_153 * lg_587[k]
                  - f_152 * lg_592[k];
    }

#pragma omp simd aligned(lg_60, lg_63, lg_70, lg_165, lg_168, lg_175, lg_195, lg_198, lg_205, \
                         lg_330, lg_333, lg_340, lg_360, lg_363, lg_370, lg_555, lg_558, \
                         lg_565, lg_585, lg_588, lg_595 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = -f_206 * lg_60[k]
                  + f_56 * lg_63[k]
                  - f_206 * lg_70[k]
                  + f_206 * lg_165[k]
                  - f_56 * lg_168[k]
                  + f_206 * lg_175[k]
                  + f_137 * lg_195[k]
                  - f_57 * lg_198[k]
                  + f_137 * lg_205[k]
                  + f_59 * lg_330[k]
                  - f_207 * lg_333[k]
                  + f_59 * lg_340[k]
                  - f_58 * lg_360[k]
                  + f_208 * lg_363[k]
                  - f_58 * lg_370[k]
                  - f_209 * lg_555[k]
                  + f_18 * lg_558[k]
                  - f_209 * lg_565[k]
                  + f_141 * lg_585[k]
                  - f_53 * lg_588[k]
                  + f_141 * lg_595[k];
    }

#pragma omp simd aligned(lg_16, lg_21, lg_91, lg_96, lg_121, lg_126, lg_226, lg_231, lg_286, \
                         lg_291, lg_421, lg_426, lg_451, lg_456, lg_481, \
                         lg_486 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_36[k] = f_210 * lg_16[k]
                  - f_210 * lg_21[k]
                  + f_210 * lg_91[k]
                  - f_210 * lg_96[k]
                  - f_211 * lg_121[k]
                  + f_211 * lg_126[k]
                  - f_210 * lg_226[k]
                  + f_210 * lg_231[k]
                  + f_212 * lg_286[k]
                  - f_212 * lg_291[k]
                  - f_210 * lg_421[k]
                  + f_210 * lg_426[k]
                  + f_211 * lg_451[k]
                  - f_211 * lg_456[k]
                  - f_212 * lg_481[k]
                  + f_212 * lg_486[k];
    }

#pragma omp simd aligned(lg_19, lg_26, lg_94, lg_101, lg_124, lg_131, lg_229, lg_236, lg_289, \
                         lg_296, lg_424, lg_431, lg_454, lg_461, lg_484, \
                         lg_491 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_37[k] = f_213 * lg_19[k]
                  - f_214 * lg_26[k]
                  + f_213 * lg_94[k]
                  - f_214 * lg_101[k]
                  - f_215 * lg_124[k]
                  + f_216 * lg_131[k]
                  - f_213 * lg_229[k]
                  + f_214 * lg_236[k]
                  + f_217 * lg_289[k]
                  - f_218 * lg_296[k]
                  - f_213 * lg_424[k]
                  + f_214 * lg_431[k]
                  + f_215 * lg_454[k]
                  - f_216 * lg_461[k]
                  - f_217 * lg_484[k]
                  + f_218 * lg_491[k];
    }

#pragma omp simd aligned(lg_16, lg_21, lg_23, lg_91, lg_96, lg_98, lg_121, lg_126, lg_128, \
                         lg_226, lg_231, lg_233, lg_286, lg_291, lg_293, lg_421, lg_426, \
                         lg_428, lg_451, lg_456, lg_458, lg_481, lg_486, \
                         lg_488 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = -f_219 * lg_16[k]
                  - f_219 * lg_21[k]
                  + f_220 * lg_23[k]
                  - f_219 * lg_91[k]
                  - f_219 * lg_96[k]
                  + f_220 * lg_98[k]
                  + f_221 * lg_121[k]
                  + f_221 * lg_126[k]
                  - f_222 * lg_128[k]
                  + f_219 * lg_226[k]
                  + f_219 * lg_231[k]
                  - f_220 * lg_233[k]
                  - f_223 * lg_286[k]
                  - f_223 * lg_291[k]
                  + f_224 * lg_293[k]
                  + f_219 * lg_421[k]
                  + f_219 * lg_426[k]
                  - f_220 * lg_428[k]
                  - f_221 * lg_451[k]
                  - f_221 * lg_456[k]
                  + f_222 * lg_458[k]
                  + f_223 * lg_481[k]
                  + f_223 * lg_486[k]
                  - f_224 * lg_488[k];
    }

#pragma omp simd aligned(lg_19, lg_26, lg_28, lg_94, lg_101, lg_103, lg_124, lg_131, lg_133, \
                         lg_229, lg_236, lg_238, lg_289, lg_296, lg_298, lg_424, lg_431, \
                         lg_433, lg_454, lg_461, lg_463, lg_484, lg_491, \
                         lg_493 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_39[k] = -f_225 * lg_19[k]
                  - f_225 * lg_26[k]
                  + f_226 * lg_28[k]
                  - f_225 * lg_94[k]
                  - f_225 * lg_101[k]
                  + f_226 * lg_103[k]
                  + f_227 * lg_124[k]
                  + f_227 * lg_131[k]
                  - f_228 * lg_133[k]
                  + f_225 * lg_229[k]
                  + f_225 * lg_236[k]
                  - f_226 * lg_238[k]
                  - f_229 * lg_289[k]
                  - f_229 * lg_296[k]
                  + f_230 * lg_298[k]
                  + f_225 * lg_424[k]
                  + f_225 * lg_431[k]
                  - f_226 * lg_433[k]
                  - f_227 * lg_454[k]
                  - f_227 * lg_461[k]
                  + f_228 * lg_463[k]
                  + f_229 * lg_484[k]
                  + f_229 * lg_491[k]
                  - f_230 * lg_493[k];
    }

#pragma omp simd aligned(lg_15, lg_18, lg_20, lg_25, lg_27, lg_29, lg_90, lg_93, lg_95, \
                         lg_100, lg_102, lg_104, lg_120, lg_123, lg_125, lg_130, lg_132, \
                         lg_134, lg_225, lg_228, lg_230, lg_235, lg_237, lg_239, lg_285, \
                         lg_288, lg_290, lg_295, lg_297, lg_299, lg_420, lg_423, lg_425, \
                         lg_430, lg_432, lg_434, lg_450, lg_453, lg_455, lg_460, lg_462, \
                         lg_464, lg_480, lg_483, lg_485, lg_490, lg_492, \
                         lg_494 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_40[k] = f_231 * lg_15[k]
                  + f_232 * lg_18[k]
                  - f_233 * lg_20[k]
                  + f_231 * lg_25[k]
                  - f_233 * lg_27[k]
                  + f_234 * lg_29[k]
                  + f_231 * lg_90[k]
                  + f_232 * lg_93[k]
                  - f_233 * lg_95[k]
                  + f_231 * lg_100[k]
                  - f_233 * lg_102[k]
                  + f_234 * lg_104[k]
                  - f_235 * lg_120[k]
                  - f_236 * lg_123[k]
                  + f_237 * lg_125[k]
                  - f_235 * lg_130[k]
                  + f_237 * lg_132[k]
                  - f_238 * lg_134[k]
                  - f_231 * lg_225[k]
                  - f_232 * lg_228[k]
                  + f_233 * lg_230[k]
                  - f_231 * lg_235[k]
                  + f_233 * lg_237[k]
                  - f_234 * lg_239[k]
                  + f_239 * lg_285[k]
                  + f_240 * lg_288[k]
                  - f_241 * lg_290[k]
                  + f_239 * lg_295[k]
                  - f_241 * lg_297[k]
                  + f_242 * lg_299[k]
                  - f_231 * lg_420[k]
                  - f_232 * lg_423[k]
                  + f_233 * lg_425[k]
                  - f_231 * lg_430[k]
                  + f_233 * lg_432[k]
                  - f_234 * lg_434[k]
                  + f_235 * lg_450[k]
                  + f_236 * lg_453[k]
                  - f_237 * lg_455[k]
                  + f_235 * lg_460[k]
                  - f_237 * lg_462[k]
                  + f_238 * lg_464[k]
                  - f_239 * lg_480[k]
                  - f_240 * lg_483[k]
                  + f_241 * lg_485[k]
                  - f_239 * lg_490[k]
                  + f_241 * lg_492[k]
                  - f_242 * lg_494[k];
    }

#pragma omp simd aligned(lg_17, lg_22, lg_24, lg_92, lg_97, lg_99, lg_122, lg_127, lg_129, \
                         lg_227, lg_232, lg_234, lg_287, lg_292, lg_294, lg_422, lg_427, \
                         lg_429, lg_452, lg_457, lg_459, lg_482, lg_487, \
                         lg_489 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_41[k] = -f_225 * lg_17[k]
                  - f_225 * lg_22[k]
                  + f_226 * lg_24[k]
                  - f_225 * lg_92[k]
                  - f_225 * lg_97[k]
                  + f_226 * lg_99[k]
                  + f_227 * lg_122[k]
                  + f_227 * lg_127[k]
                  - f_228 * lg_129[k]
                  + f_225 * lg_227[k]
                  + f_225 * lg_232[k]
                  - f_226 * lg_234[k]
                  - f_229 * lg_287[k]
                  - f_229 * lg_292[k]
                  + f_230 * lg_294[k]
                  + f_225 * lg_422[k]
                  + f_225 * lg_427[k]
                  - f_226 * lg_429[k]
                  - f_227 * lg_452[k]
                  - f_227 * lg_457[k]
                  + f_228 * lg_459[k]
                  + f_229 * lg_482[k]
                  + f_229 * lg_487[k]
                  - f_230 * lg_489[k];
    }

#pragma omp simd aligned(lg_15, lg_20, lg_25, lg_27, lg_90, lg_95, lg_100, lg_102, lg_120, \
                         lg_125, lg_130, lg_132, lg_225, lg_230, lg_235, lg_237, lg_285, \
                         lg_290, lg_295, lg_297, lg_420, lg_425, lg_430, lg_432, lg_450, \
                         lg_455, lg_460, lg_462, lg_480, lg_485, lg_490, \
                         lg_492 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = -f_243 * lg_15[k]
                  + f_244 * lg_20[k]
                  + f_243 * lg_25[k]
                  - f_244 * lg_27[k]
                  - f_243 * lg_90[k]
                  + f_244 * lg_95[k]
                  + f_243 * lg_100[k]
                  - f_244 * lg_102[k]
                  + f_245 * lg_120[k]
                  - f_246 * lg_125[k]
                  - f_245 * lg_130[k]
                  + f_246 * lg_132[k]
                  + f_243 * lg_225[k]
                  - f_244 * lg_230[k]
                  - f_243 * lg_235[k]
                  + f_244 * lg_237[k]
                  - f_247 * lg_285[k]
                  + f_248 * lg_290[k]
                  + f_247 * lg_295[k]
                  - f_248 * lg_297[k]
                  + f_243 * lg_420[k]
                  - f_244 * lg_425[k]
                  - f_243 * lg_430[k]
                  + f_244 * lg_432[k]
                  - f_245 * lg_450[k]
                  + f_246 * lg_455[k]
                  + f_245 * lg_460[k]
                  - f_246 * lg_462[k]
                  + f_247 * lg_480[k]
                  - f_248 * lg_485[k]
                  - f_247 * lg_490[k]
                  + f_248 * lg_492[k];
    }

#pragma omp simd aligned(lg_17, lg_22, lg_92, lg_97, lg_122, lg_127, lg_227, lg_232, lg_287, \
                         lg_292, lg_422, lg_427, lg_452, lg_457, lg_482, \
                         lg_487 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_43[k] = f_214 * lg_17[k]
                  - f_213 * lg_22[k]
                  + f_214 * lg_92[k]
                  - f_213 * lg_97[k]
                  - f_216 * lg_122[k]
                  + f_215 * lg_127[k]
                  - f_214 * lg_227[k]
                  + f_213 * lg_232[k]
                  + f_218 * lg_287[k]
                  - f_217 * lg_292[k]
                  - f_214 * lg_422[k]
                  + f_213 * lg_427[k]
                  + f_216 * lg_452[k]
                  - f_215 * lg_457[k]
                  - f_218 * lg_482[k]
                  + f_217 * lg_487[k];
    }

#pragma omp simd aligned(lg_15, lg_18, lg_25, lg_90, lg_93, lg_100, lg_120, lg_123, lg_130, \
                         lg_225, lg_228, lg_235, lg_285, lg_288, lg_295, lg_420, lg_423, \
                         lg_430, lg_450, lg_453, lg_460, lg_480, lg_483, \
                         lg_490 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_44[k] = f_249 * lg_15[k]
                  - f_250 * lg_18[k]
                  + f_249 * lg_25[k]
                  + f_249 * lg_90[k]
                  - f_250 * lg_93[k]
                  + f_249 * lg_100[k]
                  - f_251 * lg_120[k]
                  + f_252 * lg_123[k]
                  - f_251 * lg_130[k]
                  - f_249 * lg_225[k]
                  + f_250 * lg_228[k]
                  - f_249 * lg_235[k]
                  + f_253 * lg_285[k]
                  - f_254 * lg_288[k]
                  + f_253 * lg_295[k]
                  - f_249 * lg_420[k]
                  + f_250 * lg_423[k]
                  - f_249 * lg_430[k]
                  + f_251 * lg_450[k]
                  - f_252 * lg_453[k]
                  + f_251 * lg_460[k]
                  - f_253 * lg_480[k]
                  + f_254 * lg_483[k]
                  - f_253 * lg_490[k];
    }

#pragma omp simd aligned(lg_61, lg_66, lg_166, lg_171, lg_196, lg_201, lg_331, lg_336, lg_361, \
                         lg_366, lg_391, lg_396, lg_556, lg_561, lg_586, lg_591, lg_616, \
                         lg_621 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_45[k] = f_255 * lg_61[k]
                  - f_255 * lg_66[k]
                  + f_256 * lg_166[k]
                  - f_256 * lg_171[k]
                  - f_257 * lg_196[k]
                  + f_257 * lg_201[k]
                  + f_258 * lg_331[k]
                  - f_258 * lg_336[k]
                  - f_259 * lg_361[k]
                  + f_259 * lg_366[k]
                  + f_260 * lg_391[k]
                  - f_260 * lg_396[k]
                  - f_258 * lg_556[k]
                  + f_258 * lg_561[k]
                  + f_261 * lg_586[k]
                  - f_261 * lg_591[k]
                  - f_262 * lg_616[k]
                  + f_262 * lg_621[k];
    }

#pragma omp simd aligned(lg_64, lg_71, lg_169, lg_176, lg_199, lg_206, lg_334, lg_341, lg_364, \
                         lg_371, lg_394, lg_401, lg_559, lg_566, lg_589, lg_596, lg_619, \
                         lg_626 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_46[k] = f_263 * lg_64[k]
                  - f_264 * lg_71[k]
                  + f_265 * lg_169[k]
                  - f_266 * lg_176[k]
                  - f_267 * lg_199[k]
                  + f_268 * lg_206[k]
                  + f_264 * lg_334[k]
                  - f_269 * lg_341[k]
                  - f_270 * lg_364[k]
                  + f_271 * lg_371[k]
                  + f_272 * lg_394[k]
                  - f_273 * lg_401[k]
                  - f_264 * lg_559[k]
                  + f_269 * lg_566[k]
                  + f_268 * lg_589[k]
                  - f_274 * lg_596[k]
                  - f_273 * lg_619[k]
                  + f_275 * lg_626[k];
    }

#pragma omp simd aligned(lg_61, lg_66, lg_68, lg_166, lg_171, lg_173, lg_196, lg_201, lg_203, \
                         lg_331, lg_336, lg_338, lg_361, lg_366, lg_368, lg_391, lg_396, \
                         lg_398, lg_556, lg_561, lg_563, lg_586, lg_591, lg_593, lg_616, \
                         lg_621, lg_623 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_47[k] = -f_276 * lg_61[k]
                  - f_276 * lg_66[k]
                  + f_277 * lg_68[k]
                  - f_278 * lg_166[k]
                  - f_278 * lg_171[k]
                  + f_279 * lg_173[k]
                  + f_280 * lg_196[k]
                  + f_280 * lg_201[k]
                  - f_281 * lg_203[k]
                  - f_282 * lg_331[k]
                  - f_282 * lg_336[k]
                  + f_283 * lg_338[k]
                  + f_284 * lg_361[k]
                  + f_284 * lg_366[k]
                  - f_285 * lg_368[k]
                  - f_286 * lg_391[k]
                  - f_286 * lg_396[k]
                  + f_287 * lg_398[k]
                  + f_282 * lg_556[k]
                  + f_282 * lg_561[k]
                  - f_283 * lg_563[k]
                  - f_288 * lg_586[k]
                  - f_288 * lg_591[k]
                  + f_289 * lg_593[k]
                  + f_290 * lg_616[k]
                  + f_290 * lg_621[k]
                  - f_291 * lg_623[k];
    }

#pragma omp simd aligned(lg_64, lg_71, lg_73, lg_169, lg_176, lg_178, lg_199, lg_206, lg_208, \
                         lg_334, lg_341, lg_343, lg_364, lg_371, lg_373, lg_394, lg_401, \
                         lg_403, lg_559, lg_566, lg_568, lg_589, lg_596, lg_598, lg_619, \
                         lg_626, lg_628 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_48[k] = -f_292 * lg_64[k]
                  - f_292 * lg_71[k]
                  + f_293 * lg_73[k]
                  - f_294 * lg_169[k]
                  - f_294 * lg_176[k]
                  + f_295 * lg_178[k]
                  + f_296 * lg_199[k]
                  + f_296 * lg_206[k]
                  - f_297 * lg_208[k]
                  - f_298 * lg_334[k]
                  - f_298 * lg_341[k]
                  + f_299 * lg_343[k]
                  + f_300 * lg_364[k]
                  + f_300 * lg_371[k]
                  - f_301 * lg_373[k]
                  - f_302 * lg_394[k]
                  - f_302 * lg_401[k]
                  + f_303 * lg_403[k]
                  + f_298 * lg_559[k]
                  + f_298 * lg_566[k]
                  - f_299 * lg_568[k]
                  - f_295 * lg_589[k]
                  - f_295 * lg_596[k]
                  + f_304 * lg_598[k]
                  + f_305 * lg_619[k]
                  + f_305 * lg_626[k]
                  - f_306 * lg_628[k];
    }

#pragma omp simd aligned(lg_60, lg_63, lg_65, lg_70, lg_72, lg_74, lg_165, lg_168, lg_170, \
                         lg_175, lg_177, lg_179, lg_195, lg_198, lg_200, lg_205, lg_207, \
                         lg_209, lg_330, lg_333, lg_335, lg_340, lg_342, lg_344, lg_360, \
                         lg_363, lg_365, lg_370, lg_372, lg_374, lg_390, lg_393, lg_395, \
                         lg_400, lg_402, lg_404, lg_555, lg_558, lg_560, lg_565, lg_567, \
                         lg_569, lg_585, lg_588, lg_590, lg_595, lg_597, lg_599, lg_615, \
                         lg_618, lg_620, lg_625, lg_627, lg_629 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_49[k] = f_307 * lg_60[k]
                  + f_308 * lg_63[k]
                  - f_309 * lg_65[k]
                  + f_307 * lg_70[k]
                  - f_309 * lg_72[k]
                  + f_310 * lg_74[k]
                  + f_311 * lg_165[k]
                  + f_312 * lg_168[k]
                  - f_313 * lg_170[k]
                  + f_311 * lg_175[k]
                  - f_313 * lg_177[k]
                  + f_314 * lg_179[k]
                  - f_315 * lg_195[k]
                  - f_313 * lg_198[k]
                  + f_316 * lg_200[k]
                  - f_315 * lg_205[k]
                  + f_316 * lg_207[k]
                  - f_317 * lg_209[k]
                  + f_318 * lg_330[k]
                  + f_319 * lg_333[k]
                  - f_310 * lg_335[k]
                  + f_318 * lg_340[k]
                  - f_310 * lg_342[k]
                  + f_320 * lg_344[k]
                  - f_314 * lg_360[k]
                  - f_321 * lg_363[k]
                  + f_322 * lg_365[k]
                  - f_314 * lg_370[k]
                  + f_322 * lg_372[k]
                  - f_323 * lg_374[k]
                  + f_324 * lg_390[k]
                  + f_325 * lg_393[k]
                  - f_326 * lg_395[k]
                  + f_324 * lg_400[k]
                  - f_326 * lg_402[k]
                  + f_327 * lg_404[k]
                  - f_318 * lg_555[k]
                  - f_319 * lg_558[k]
                  + f_310 * lg_560[k]
                  - f_318 * lg_565[k]
                  + f_310 * lg_567[k]
                  - f_320 * lg_569[k]
                  + f_328 * lg_585[k]
                  + f_314 * lg_588[k]
                  - f_317 * lg_590[k]
                  + f_328 * lg_595[k]
                  - f_317 * lg_597[k]
                  + f_329 * lg_599[k]
                  - f_330 * lg_615[k]
                  - f_331 * lg_618[k]
                  + f_327 * lg_620[k]
                  - f_330 * lg_625[k]
                  + f_327 * lg_627[k]
                  - f_332 * lg_629[k];
    }

#pragma omp simd aligned(lg_62, lg_67, lg_69, lg_167, lg_172, lg_174, lg_197, lg_202, lg_204, \
                         lg_332, lg_337, lg_339, lg_362, lg_367, lg_369, lg_392, lg_397, \
                         lg_399, lg_557, lg_562, lg_564, lg_587, lg_592, lg_594, lg_617, \
                         lg_622, lg_624 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_50[k] = -f_292 * lg_62[k]
                  - f_292 * lg_67[k]
                  + f_293 * lg_69[k]
                  - f_294 * lg_167[k]
                  - f_294 * lg_172[k]
                  + f_295 * lg_174[k]
                  + f_296 * lg_197[k]
                  + f_296 * lg_202[k]
                  - f_297 * lg_204[k]
                  - f_298 * lg_332[k]
                  - f_298 * lg_337[k]
                  + f_299 * lg_339[k]
                  + f_300 * lg_362[k]
                  + f_300 * lg_367[k]
                  - f_301 * lg_369[k]
                  - f_302 * lg_392[k]
                  - f_302 * lg_397[k]
                  + f_303 * lg_399[k]
                  + f_298 * lg_557[k]
                  + f_298 * lg_562[k]
                  - f_299 * lg_564[k]
                  - f_295 * lg_587[k]
                  - f_295 * lg_592[k]
                  + f_304 * lg_594[k]
                  + f_305 * lg_617[k]
                  + f_305 * lg_622[k]
                  - f_306 * lg_624[k];
    }

#pragma omp simd aligned(lg_60, lg_65, lg_70, lg_72, lg_165, lg_170, lg_175, lg_177, lg_195, \
                         lg_200, lg_205, lg_207, lg_330, lg_335, lg_340, lg_342, lg_360, \
                         lg_365, lg_370, lg_372, lg_390, lg_395, lg_400, lg_402, lg_555, \
                         lg_560, lg_565, lg_567, lg_585, lg_590, lg_595, lg_597, lg_615, \
                         lg_620, lg_625, lg_627 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_51[k] = -f_333 * lg_60[k]
                  + f_334 * lg_65[k]
                  + f_333 * lg_70[k]
                  - f_334 * lg_72[k]
                  - f_335 * lg_165[k]
                  + f_336 * lg_170[k]
                  + f_335 * lg_175[k]
                  - f_336 * lg_177[k]
                  + f_337 * lg_195[k]
                  - f_338 * lg_200[k]
                  - f_337 * lg_205[k]
                  + f_338 * lg_207[k]
                  - f_339 * lg_330[k]
                  + f_276 * lg_335[k]
                  + f_339 * lg_340[k]
                  - f_276 * lg_342[k]
                  + f_288 * lg_360[k]
                  - f_289 * lg_365[k]
                  - f_288 * lg_370[k]
                  + f_289 * lg_372[k]
                  - f_340 * lg_390[k]
                  + f_341 * lg_395[k]
                  + f_340 * lg_400[k]
                  - f_341 * lg_402[k]
                  + f_339 * lg_555[k]
                  - f_276 * lg_560[k]
                  - f_339 * lg_565[k]
                  + f_276 * lg_567[k]
                  - f_342 * lg_585[k]
                  + f_280 * lg_590[k]
                  + f_342 * lg_595[k]
                  - f_280 * lg_597[k]
                  + f_343 * lg_615[k]
                  - f_286 * lg_620[k]
                  - f_343 * lg_625[k]
                  + f_286 * lg_627[k];
    }

#pragma omp simd aligned(lg_62, lg_67, lg_167, lg_172, lg_197, lg_202, lg_332, lg_337, lg_362, \
                         lg_367, lg_392, lg_397, lg_557, lg_562, lg_587, lg_592, lg_617, \
                         lg_622 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_52[k] = f_264 * lg_62[k]
                  - f_263 * lg_67[k]
                  + f_266 * lg_167[k]
                  - f_265 * lg_172[k]
                  - f_268 * lg_197[k]
                  + f_267 * lg_202[k]
                  + f_269 * lg_332[k]
                  - f_264 * lg_337[k]
                  - f_271 * lg_362[k]
                  + f_270 * lg_367[k]
                  + f_273 * lg_392[k]
                  - f_272 * lg_397[k]
                  - f_269 * lg_557[k]
                  + f_264 * lg_562[k]
                  + f_274 * lg_587[k]
                  - f_268 * lg_592[k]
                  - f_275 * lg_617[k]
                  + f_273 * lg_622[k];
    }

#pragma omp simd aligned(lg_60, lg_63, lg_70, lg_165, lg_168, lg_175, lg_195, lg_198, lg_205, \
                         lg_330, lg_333, lg_340, lg_360, lg_363, lg_370, lg_390, lg_393, \
                         lg_400, lg_555, lg_558, lg_565, lg_585, lg_588, lg_595, lg_615, \
                         lg_618, lg_625 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_53[k] = f_344 * lg_60[k]
                  - f_345 * lg_63[k]
                  + f_344 * lg_70[k]
                  + f_346 * lg_165[k]
                  - f_347 * lg_168[k]
                  + f_346 * lg_175[k]
                  - f_256 * lg_195[k]
                  + f_348 * lg_198[k]
                  - f_256 * lg_205[k]
                  + f_349 * lg_330[k]
                  - f_350 * lg_333[k]
                  + f_349 * lg_340[k]
                  - f_351 * lg_360[k]
                  + f_257 * lg_363[k]
                  - f_351 * lg_370[k]
                  + f_352 * lg_390[k]
                  - f_353 * lg_393[k]
                  + f_352 * lg_400[k]
                  - f_349 * lg_555[k]
                  + f_350 * lg_558[k]
                  - f_349 * lg_565[k]
                  + f_354 * lg_585[k]
                  - f_355 * lg_588[k]
                  + f_354 * lg_595[k]
                  - f_356 * lg_615[k]
                  + f_357 * lg_618[k]
                  - f_356 * lg_625[k];
    }

#pragma omp simd aligned(lg_16, lg_21, lg_91, lg_96, lg_121, lg_126, lg_226, lg_231, lg_256, \
                         lg_261, lg_286, lg_291, lg_421, lg_426, lg_451, lg_456, lg_481, \
                         lg_486, lg_511, lg_516 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_54[k] = -f_358 * lg_16[k]
                  + f_358 * lg_21[k]
                  - f_359 * lg_91[k]
                  + f_359 * lg_96[k]
                  + f_360 * lg_121[k]
                  - f_360 * lg_126[k]
                  - f_359 * lg_226[k]
                  + f_359 * lg_231[k]
                  + f_361 * lg_256[k]
                  - f_361 * lg_261[k]
                  - f_362 * lg_286[k]
                  + f_362 * lg_291[k]
                  - f_358 * lg_421[k]
                  + f_358 * lg_426[k]
                  + f_360 * lg_451[k]
                  - f_360 * lg_456[k]
                  - f_362 * lg_481[k]
                  + f_362 * lg_486[k]
                  + f_363 * lg_511[k]
                  - f_363 * lg_516[k];
    }

#pragma omp simd aligned(lg_19, lg_26, lg_94, lg_101, lg_124, lg_131, lg_229, lg_236, lg_259, \
                         lg_266, lg_289, lg_296, lg_424, lg_431, lg_454, lg_461, lg_484, \
                         lg_491, lg_514, lg_521 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_55[k] = -4.921875 * lg_19[k]
                  + 1.640625 * lg_26[k]
                  - 14.765625 * lg_94[k]
                  + 4.921875 * lg_101[k]
                  + 147.65625 * lg_124[k]
                  - 49.21875 * lg_131[k]
                  - 14.765625 * lg_229[k]
                  + 4.921875 * lg_236[k]
                  + 295.3125 * lg_259[k]
                  - 98.4375 * lg_266[k]
                  - 393.75 * lg_289[k]
                  + 131.25 * lg_296[k]
                  - 4.921875 * lg_424[k]
                  + 1.640625 * lg_431[k]
                  + 147.65625 * lg_454[k]
                  - 49.21875 * lg_461[k]
                  - 393.75 * lg_484[k]
                  + 131.25 * lg_491[k]
                  + 157.5 * lg_514[k]
                  - 52.5 * lg_521[k];
    }

#pragma omp simd aligned(lg_16, lg_21, lg_23, lg_91, lg_96, lg_98, lg_121, lg_126, lg_128, \
                         lg_226, lg_231, lg_233, lg_256, lg_261, lg_263, lg_286, lg_291, \
                         lg_293, lg_421, lg_426, lg_428, lg_451, lg_456, lg_458, lg_481, \
                         lg_486, lg_488, lg_511, lg_516, lg_518 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_56[k] = f_364 * lg_16[k]
                  + f_364 * lg_21[k]
                  - f_365 * lg_23[k]
                  + f_366 * lg_91[k]
                  + f_366 * lg_96[k]
                  - f_367 * lg_98[k]
                  - f_368 * lg_121[k]
                  - f_368 * lg_126[k]
                  + f_369 * lg_128[k]
                  + f_366 * lg_226[k]
                  + f_366 * lg_231[k]
                  - f_367 * lg_233[k]
                  - f_370 * lg_256[k]
                  - f_370 * lg_261[k]
                  + f_371 * lg_263[k]
                  + f_372 * lg_286[k]
                  + f_372 * lg_291[k]
                  - f_373 * lg_293[k]
                  + f_364 * lg_421[k]
                  + f_364 * lg_426[k]
                  - f_365 * lg_428[k]
                  - f_368 * lg_451[k]
                  - f_368 * lg_456[k]
                  + f_369 * lg_458[k]
                  + f_372 * lg_481[k]
                  + f_372 * lg_486[k]
                  - f_373 * lg_488[k]
                  - f_374 * lg_511[k]
                  - f_374 * lg_516[k]
                  + f_375 * lg_518[k];
    }

#pragma omp simd aligned(lg_19, lg_26, lg_28, lg_94, lg_101, lg_103, lg_124, lg_131, lg_133, \
                         lg_229, lg_236, lg_238, lg_259, lg_266, lg_268, lg_289, lg_296, \
                         lg_298, lg_424, lg_431, lg_433, lg_454, lg_461, lg_463, lg_484, \
                         lg_491, lg_493, lg_514, lg_521, lg_523 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_57[k] = f_376 * lg_19[k]
                  + f_376 * lg_26[k]
                  - f_377 * lg_28[k]
                  + f_378 * lg_94[k]
                  + f_378 * lg_101[k]
                  - f_379 * lg_103[k]
                  - f_380 * lg_124[k]
                  - f_380 * lg_131[k]
                  + f_381 * lg_133[k]
                  + f_378 * lg_229[k]
                  + f_378 * lg_236[k]
                  - f_379 * lg_238[k]
                  - f_382 * lg_259[k]
                  - f_382 * lg_266[k]
                  + f_383 * lg_268[k]
                  + f_383 * lg_289[k]
                  + f_383 * lg_296[k]
                  - f_384 * lg_298[k]
                  + f_376 * lg_424[k]
                  + f_376 * lg_431[k]
                  - f_377 * lg_433[k]
                  - f_380 * lg_454[k]
                  - f_380 * lg_461[k]
                  + f_381 * lg_463[k]
                  + f_383 * lg_484[k]
                  + f_383 * lg_491[k]
                  - f_384 * lg_493[k]
                  - f_385 * lg_514[k]
                  - f_385 * lg_521[k]
                  + f_386 * lg_523[k];
    }

#pragma omp simd aligned(lg_15, lg_18, lg_20, lg_25, lg_27, lg_29, lg_90, lg_93, lg_95, \
                         lg_100, lg_102, lg_104, lg_120, lg_123, lg_125, lg_130, lg_132, \
                         lg_134, lg_225, lg_228, lg_230, lg_235, lg_237, lg_239, lg_255, \
                         lg_258, lg_260, lg_265, lg_267, lg_269, lg_285, lg_288, lg_290, \
                         lg_295, lg_297, lg_299, lg_420, lg_423, lg_425, lg_430, lg_432, \
                         lg_434, lg_450, lg_453, lg_455, lg_460, lg_462, lg_464, lg_480, \
                         lg_483, lg_485, lg_490, lg_492, lg_494, lg_510, lg_513, lg_515, \
                         lg_520, lg_522, lg_524 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_58[k] = -f_387 * lg_15[k]
                  - f_388 * lg_18[k]
                  + f_389 * lg_20[k]
                  - f_387 * lg_25[k]
                  + f_389 * lg_27[k]
                  - f_390 * lg_29[k]
                  - f_391 * lg_90[k]
                  - f_392 * lg_93[k]
                  + f_393 * lg_95[k]
                  - f_391 * lg_100[k]
                  + f_393 * lg_102[k]
                  - f_389 * lg_104[k]
                  + f_394 * lg_120[k]
                  + f_395 * lg_123[k]
                  - f_396 * lg_125[k]
                  + f_394 * lg_130[k]
                  - f_396 * lg_132[k]
                  + f_397 * lg_134[k]
                  - f_391 * lg_225[k]
                  - f_392 * lg_228[k]
                  + f_393 * lg_230[k]
                  - f_391 * lg_235[k]
                  + f_393 * lg_237[k]
                  - f_389 * lg_239[k]
                  + f_395 * lg_255[k]
                  + f_398 * lg_258[k]
                  - f_399 * lg_260[k]
                  + f_395 * lg_265[k]
                  - f_399 * lg_267[k]
                  + f_400 * lg_269[k]
                  - f_397 * lg_285[k]
                  - f_400 * lg_288[k]
                  + f_401 * lg_290[k]
                  - f_397 * lg_295[k]
                  + f_401 * lg_297[k]
                  - f_402 * lg_299[k]
                  - f_387 * lg_420[k]
                  - f_388 * lg_423[k]
                  + f_389 * lg_425[k]
                  - f_387 * lg_430[k]
                  + f_389 * lg_432[k]
                  - f_390 * lg_434[k]
                  + f_394 * lg_450[k]
                  + f_395 * lg_453[k]
                  - f_396 * lg_455[k]
                  + f_394 * lg_460[k]
                  - f_396 * lg_462[k]
                  + f_397 * lg_464[k]
                  - f_397 * lg_480[k]
                  - f_400 * lg_483[k]
                  + f_401 * lg_485[k]
                  - f_397 * lg_490[k]
                  + f_401 * lg_492[k]
                  - f_402 * lg_494[k]
                  + f_403 * lg_510[k]
                  + f_404 * lg_513[k]
                  - f_405 * lg_515[k]
                  + f_403 * lg_520[k]
                  - f_405 * lg_522[k]
                  + f_406 * lg_524[k];
    }

#pragma omp simd aligned(lg_17, lg_22, lg_24, lg_92, lg_97, lg_99, lg_122, lg_127, lg_129, \
                         lg_227, lg_232, lg_234, lg_257, lg_262, lg_264, lg_287, lg_292, \
                         lg_294, lg_422, lg_427, lg_429, lg_452, lg_457, lg_459, lg_482, \
                         lg_487, lg_489, lg_512, lg_517, lg_519 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_59[k] = f_376 * lg_17[k]
                  + f_376 * lg_22[k]
                  - f_377 * lg_24[k]
                  + f_378 * lg_92[k]
                  + f_378 * lg_97[k]
                  - f_379 * lg_99[k]
                  - f_380 * lg_122[k]
                  - f_380 * lg_127[k]
                  + f_381 * lg_129[k]
                  + f_378 * lg_227[k]
                  + f_378 * lg_232[k]
                  - f_379 * lg_234[k]
                  - f_382 * lg_257[k]
                  - f_382 * lg_262[k]
                  + f_383 * lg_264[k]
                  + f_383 * lg_287[k]
                  + f_383 * lg_292[k]
                  - f_384 * lg_294[k]
                  + f_376 * lg_422[k]
                  + f_376 * lg_427[k]
                  - f_377 * lg_429[k]
                  - f_380 * lg_452[k]
                  - f_380 * lg_457[k]
                  + f_381 * lg_459[k]
                  + f_383 * lg_482[k]
                  + f_383 * lg_487[k]
                  - f_384 * lg_489[k]
                  - f_385 * lg_512[k]
                  - f_385 * lg_517[k]
                  + f_386 * lg_519[k];
    }

#pragma omp simd aligned(lg_15, lg_20, lg_25, lg_27, lg_90, lg_95, lg_100, lg_102, lg_120, \
                         lg_125, lg_130, lg_132, lg_225, lg_230, lg_235, lg_237, lg_255, \
                         lg_260, lg_265, lg_267, lg_285, lg_290, lg_295, lg_297, lg_420, \
                         lg_425, lg_430, lg_432, lg_450, lg_455, lg_460, lg_462, lg_480, \
                         lg_485, lg_490, lg_492, lg_510, lg_515, lg_520, \
                         lg_522 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_60[k] = f_407 * lg_15[k]
                  - f_366 * lg_20[k]
                  - f_407 * lg_25[k]
                  + f_366 * lg_27[k]
                  + f_408 * lg_90[k]
                  - f_409 * lg_95[k]
                  - f_408 * lg_100[k]
                  + f_409 * lg_102[k]
                  - f_410 * lg_120[k]
                  + f_411 * lg_125[k]
                  + f_410 * lg_130[k]
                  - f_411 * lg_132[k]
                  + f_408 * lg_225[k]
                  - f_409 * lg_230[k]
                  - f_408 * lg_235[k]
                  + f_409 * lg_237[k]
                  - f_368 * lg_255[k]
                  + f_369 * lg_260[k]
                  + f_368 * lg_265[k]
                  - f_369 * lg_267[k]
                  + f_412 * lg_285[k]
                  - f_413 * lg_290[k]
                  - f_412 * lg_295[k]
                  + f_413 * lg_297[k]
                  + f_407 * lg_420[k]
                  - f_366 * lg_425[k]
                  - f_407 * lg_430[k]
                  + f_366 * lg_432[k]
                  - f_410 * lg_450[k]
                  + f_411 * lg_455[k]
                  + f_410 * lg_460[k]
                  - f_411 * lg_462[k]
                  + f_412 * lg_480[k]
                  - f_413 * lg_485[k]
                  - f_412 * lg_490[k]
                  + f_413 * lg_492[k]
                  - f_414 * lg_510[k]
                  + f_415 * lg_515[k]
                  + f_414 * lg_520[k]
                  - f_415 * lg_522[k];
    }

#pragma omp simd aligned(lg_17, lg_22, lg_92, lg_97, lg_122, lg_127, lg_227, lg_232, lg_257, \
                         lg_262, lg_287, lg_292, lg_422, lg_427, lg_452, lg_457, lg_482, \
                         lg_487, lg_512, lg_517 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_61[k] = -1.640625 * lg_17[k]
                  + 4.921875 * lg_22[k]
                  - 4.921875 * lg_92[k]
                  + 14.765625 * lg_97[k]
                  + 49.21875 * lg_122[k]
                  - 147.65625 * lg_127[k]
                  - 4.921875 * lg_227[k]
                  + 14.765625 * lg_232[k]
                  + 98.4375 * lg_257[k]
                  - 295.3125 * lg_262[k]
                  - 131.25 * lg_287[k]
                  + 393.75 * lg_292[k]
                  - 1.640625 * lg_422[k]
                  + 4.921875 * lg_427[k]
                  + 49.21875 * lg_452[k]
                  - 147.65625 * lg_457[k]
                  - 131.25 * lg_482[k]
                  + 393.75 * lg_487[k]
                  + 52.5 * lg_512[k]
                  - 157.5 * lg_517[k];
    }

#pragma omp simd aligned(lg_15, lg_18, lg_25, lg_90, lg_93, lg_100, lg_120, lg_123, lg_130, \
                         lg_225, lg_228, lg_235, lg_255, lg_258, lg_265, lg_285, lg_288, \
                         lg_295, lg_420, lg_423, lg_430, lg_450, lg_453, lg_460, lg_480, \
                         lg_483, lg_490, lg_510, lg_513, lg_520 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_62[k] = -f_416 * lg_15[k]
                  + f_417 * lg_18[k]
                  - f_416 * lg_25[k]
                  - f_418 * lg_90[k]
                  + f_419 * lg_93[k]
                  - f_418 * lg_100[k]
                  + f_420 * lg_120[k]
                  - f_421 * lg_123[k]
                  + f_420 * lg_130[k]
                  - f_418 * lg_225[k]
                  + f_419 * lg_228[k]
                  - f_418 * lg_235[k]
                  + f_422 * lg_255[k]
                  - f_423 * lg_258[k]
                  + f_422 * lg_265[k]
                  - f_424 * lg_285[k]
                  + f_425 * lg_288[k]
                  - f_424 * lg_295[k]
                  - f_416 * lg_420[k]
                  + f_417 * lg_423[k]
                  - f_416 * lg_430[k]
                  + f_420 * lg_450[k]
                  - f_421 * lg_453[k]
                  + f_420 * lg_460[k]
                  - f_424 * lg_480[k]
                  + f_425 * lg_483[k]
                  - f_424 * lg_490[k]
                  + f_426 * lg_510[k]
                  - f_427 * lg_513[k]
                  + f_426 * lg_520[k];
    }

#pragma omp simd aligned(lg_61, lg_66, lg_166, lg_171, lg_196, lg_201, lg_331, lg_336, lg_361, \
                         lg_366, lg_391, lg_396, lg_556, lg_561, lg_586, lg_591, lg_616, \
                         lg_621, lg_646, lg_651 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_63[k] = -f_428 * lg_61[k]
                  + f_428 * lg_66[k]
                  - f_429 * lg_166[k]
                  + f_429 * lg_171[k]
                  + f_430 * lg_196[k]
                  - f_430 * lg_201[k]
                  - f_429 * lg_331[k]
                  + f_429 * lg_336[k]
                  + f_431 * lg_361[k]
                  - f_431 * lg_366[k]
                  - f_432 * lg_391[k]
                  + f_432 * lg_396[k]
                  - f_428 * lg_556[k]
                  + f_428 * lg_561[k]
                  + f_430 * lg_586[k]
                  - f_430 * lg_591[k]
                  - f_432 * lg_616[k]
                  + f_432 * lg_621[k]
                  + f_433 * lg_646[k]
                  - f_433 * lg_651[k];
    }

#pragma omp simd aligned(lg_64, lg_71, lg_169, lg_176, lg_199, lg_206, lg_334, lg_341, lg_364, \
                         lg_371, lg_394, lg_401, lg_559, lg_566, lg_589, lg_596, lg_619, \
                         lg_626, lg_649, lg_656 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_64[k] = -f_434 * lg_64[k]
                  + f_435 * lg_71[k]
                  - f_436 * lg_169[k]
                  + f_434 * lg_176[k]
                  + f_437 * lg_199[k]
                  - f_438 * lg_206[k]
                  - f_436 * lg_334[k]
                  + f_434 * lg_341[k]
                  + f_439 * lg_364[k]
                  - f_440 * lg_371[k]
                  - f_441 * lg_394[k]
                  + f_442 * lg_401[k]
                  - f_434 * lg_559[k]
                  + f_435 * lg_566[k]
                  + f_437 * lg_589[k]
                  - f_438 * lg_596[k]
                  - f_441 * lg_619[k]
                  + f_442 * lg_626[k]
                  + f_443 * lg_649[k]
                  - f_444 * lg_656[k];
    }

#pragma omp simd aligned(lg_61, lg_66, lg_68, lg_166, lg_171, lg_173, lg_196, lg_201, lg_203, \
                         lg_331, lg_336, lg_338, lg_361, lg_366, lg_368, lg_391, lg_396, \
                         lg_398, lg_556, lg_561, lg_563, lg_586, lg_591, lg_593, lg_616, \
                         lg_621, lg_623, lg_646, lg_651, lg_653 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_65[k] = f_445 * lg_61[k]
                  + f_445 * lg_66[k]
                  - f_446 * lg_68[k]
                  + f_447 * lg_166[k]
                  + f_447 * lg_171[k]
                  - f_448 * lg_173[k]
                  - f_449 * lg_196[k]
                  - f_449 * lg_201[k]
                  + f_450 * lg_203[k]
                  + f_447 * lg_331[k]
                  + f_447 * lg_336[k]
                  - f_448 * lg_338[k]
                  - f_451 * lg_361[k]
                  - f_451 * lg_366[k]
                  + f_452 * lg_368[k]
                  + f_453 * lg_391[k]
                  + f_453 * lg_396[k]
                  - f_454 * lg_398[k]
                  + f_445 * lg_556[k]
                  + f_445 * lg_561[k]
                  - f_446 * lg_563[k]
                  - f_449 * lg_586[k]
                  - f_449 * lg_591[k]
                  + f_450 * lg_593[k]
                  + f_453 * lg_616[k]
                  + f_453 * lg_621[k]
                  - f_454 * lg_623[k]
                  - f_455 * lg_646[k]
                  - f_455 * lg_651[k]
                  + f_456 * lg_653[k];
    }

#pragma omp simd aligned(lg_64, lg_71, lg_73, lg_169, lg_176, lg_178, lg_199, lg_206, lg_208, \
                         lg_334, lg_341, lg_343, lg_364, lg_371, lg_373, lg_394, lg_401, \
                         lg_403, lg_559, lg_566, lg_568, lg_589, lg_596, lg_598, lg_619, \
                         lg_626, lg_628, lg_649, lg_656, lg_658 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_66[k] = f_457 * lg_64[k]
                  + f_457 * lg_71[k]
                  - f_458 * lg_73[k]
                  + f_459 * lg_169[k]
                  + f_459 * lg_176[k]
                  - f_460 * lg_178[k]
                  - f_461 * lg_199[k]
                  - f_461 * lg_206[k]
                  + f_462 * lg_208[k]
                  + f_459 * lg_334[k]
                  + f_459 * lg_341[k]
                  - f_460 * lg_343[k]
                  - f_463 * lg_364[k]
                  - f_463 * lg_371[k]
                  + f_464 * lg_373[k]
                  + f_465 * lg_394[k]
                  + f_465 * lg_401[k]
                  - f_466 * lg_403[k]
                  + f_457 * lg_559[k]
                  + f_457 * lg_566[k]
                  - f_458 * lg_568[k]
                  - f_461 * lg_589[k]
                  - f_461 * lg_596[k]
                  + f_462 * lg_598[k]
                  + f_465 * lg_619[k]
                  + f_465 * lg_626[k]
                  - f_466 * lg_628[k]
                  - f_467 * lg_649[k]
                  - f_467 * lg_656[k]
                  + f_468 * lg_658[k];
    }

#pragma omp simd aligned(lg_60, lg_63, lg_65, lg_70, lg_72, lg_74, lg_165, lg_168, lg_170, \
                         lg_175, lg_177, lg_179, lg_195, lg_198, lg_200, lg_205, lg_207, \
                         lg_209, lg_330, lg_333, lg_335, lg_340, lg_342, lg_344, lg_360, \
                         lg_363, lg_365, lg_370, lg_372, lg_374, lg_390, lg_393, lg_395, \
                         lg_400, lg_402, lg_404, lg_555, lg_558, lg_560, lg_565, lg_567, \
                         lg_569, lg_585, lg_588, lg_590, lg_595, lg_597, lg_599, lg_615, \
                         lg_618, lg_620, lg_625, lg_627, lg_629, lg_645, lg_648, lg_650, \
                         lg_655, lg_657, lg_659 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_67[k] = -1.23046875 * lg_60[k]
                  - 2.4609375 * lg_63[k]
                  + 9.84375 * lg_65[k]
                  - 1.23046875 * lg_70[k]
                  + 9.84375 * lg_72[k]
                  - 3.28125 * lg_74[k]
                  - 3.69140625 * lg_165[k]
                  - 7.3828125 * lg_168[k]
                  + 29.53125 * lg_170[k]
                  - 3.69140625 * lg_175[k]
                  + 29.53125 * lg_177[k]
                  - 9.84375 * lg_179[k]
                  + 9.84375 * lg_195[k]
                  + 19.6875 * lg_198[k]
                  - 78.75 * lg_200[k]
                  + 9.84375 * lg_205[k]
                  - 78.75 * lg_207[k]
                  + 26.25 * lg_209[k]
                  - 3.69140625 * lg_330[k]
                  - 7.3828125 * lg_333[k]
                  + 29.53125 * lg_335[k]
                  - 3.69140625 * lg_340[k]
                  + 29.53125 * lg_342[k]
                  - 9.84375 * lg_344[k]
                  + 19.6875 * lg_360[k]
                  + 39.375 * lg_363[k]
                  - 157.5 * lg_365[k]
                  + 19.6875 * lg_370[k]
                  - 157.5 * lg_372[k]
                  + 52.5 * lg_374[k]
                  - 11.8125 * lg_390[k]
                  - 23.625 * lg_393[k]
                  + 94.5 * lg_395[k]
                  - 11.8125 * lg_400[k]
                  + 94.5 * lg_402[k]
                  - 31.5 * lg_404[k]
                  - 1.23046875 * lg_555[k]
                  - 2.4609375 * lg_558[k]
                  + 9.84375 * lg_560[k]
                  - 1.23046875 * lg_565[k]
                  + 9.84375 * lg_567[k]
                  - 3.28125 * lg_569[k]
                  + 9.84375 * lg_585[k]
                  + 19.6875 * lg_588[k]
                  - 78.75 * lg_590[k]
                  + 9.84375 * lg_595[k]
                  - 78.75 * lg_597[k]
                  + 26.25 * lg_599[k]
                  - 11.8125 * lg_615[k]
                  - 23.625 * lg_618[k]
                  + 94.5 * lg_620[k]
                  - 11.8125 * lg_625[k]
                  + 94.5 * lg_627[k]
                  - 31.5 * lg_629[k]
                  + 2.25 * lg_645[k]
                  + 4.5 * lg_648[k]
                  - 18.0 * lg_650[k]
                  + 2.25 * lg_655[k]
                  - 18.0 * lg_657[k]
                  + 6.0 * lg_659[k];
    }

#pragma omp simd aligned(lg_62, lg_67, lg_69, lg_167, lg_172, lg_174, lg_197, lg_202, lg_204, \
                         lg_332, lg_337, lg_339, lg_362, lg_367, lg_369, lg_392, lg_397, \
                         lg_399, lg_557, lg_562, lg_564, lg_587, lg_592, lg_594, lg_617, \
                         lg_622, lg_624, lg_647, lg_652, lg_654 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_68[k] = f_457 * lg_62[k]
                  + f_457 * lg_67[k]
                  - f_458 * lg_69[k]
                  + f_459 * lg_167[k]
                  + f_459 * lg_172[k]
                  - f_460 * lg_174[k]
                  - f_461 * lg_197[k]
                  - f_461 * lg_202[k]
                  + f_462 * lg_204[k]
                  + f_459 * lg_332[k]
                  + f_459 * lg_337[k]
                  - f_460 * lg_339[k]
                  - f_463 * lg_362[k]
                  - f_463 * lg_367[k]
                  + f_464 * lg_369[k]
                  + f_465 * lg_392[k]
                  + f_465 * lg_397[k]
                  - f_466 * lg_399[k]
                  + f_457 * lg_557[k]
                  + f_457 * lg_562[k]
                  - f_458 * lg_564[k]
                  - f_461 * lg_587[k]
                  - f_461 * lg_592[k]
                  + f_462 * lg_594[k]
                  + f_465 * lg_617[k]
                  + f_465 * lg_622[k]
                  - f_466 * lg_624[k]
                  - f_467 * lg_647[k]
                  - f_467 * lg_652[k]
                  + f_468 * lg_654[k];
    }

#pragma omp simd aligned(lg_60, lg_65, lg_70, lg_72, lg_165, lg_170, lg_175, lg_177, lg_195, \
                         lg_200, lg_205, lg_207, lg_330, lg_335, lg_340, lg_342, lg_360, \
                         lg_365, lg_370, lg_372, lg_390, lg_395, lg_400, lg_402, lg_555, \
                         lg_560, lg_565, lg_567, lg_585, lg_590, lg_595, lg_597, lg_615, \
                         lg_620, lg_625, lg_627, lg_645, lg_650, lg_655, \
                         lg_657 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_69[k] = f_469 * lg_60[k]
                  - f_447 * lg_65[k]
                  - f_469 * lg_70[k]
                  + f_447 * lg_72[k]
                  + f_470 * lg_165[k]
                  - f_471 * lg_170[k]
                  - f_470 * lg_175[k]
                  + f_471 * lg_177[k]
                  - f_472 * lg_195[k]
                  + f_473 * lg_200[k]
                  + f_472 * lg_205[k]
                  - f_473 * lg_207[k]
                  + f_470 * lg_330[k]
                  - f_471 * lg_335[k]
                  - f_470 * lg_340[k]
                  + f_471 * lg_342[k]
                  - f_449 * lg_360[k]
                  + f_450 * lg_365[k]
                  + f_449 * lg_370[k]
                  - f_450 * lg_372[k]
                  + f_474 * lg_390[k]
                  - f_475 * lg_395[k]
                  - f_474 * lg_400[k]
                  + f_475 * lg_402[k]
                  + f_469 * lg_555[k]
                  - f_447 * lg_560[k]
                  - f_469 * lg_565[k]
                  + f_447 * lg_567[k]
                  - f_472 * lg_585[k]
                  + f_473 * lg_590[k]
                  + f_472 * lg_595[k]
                  - f_473 * lg_597[k]
                  + f_474 * lg_615[k]
                  - f_475 * lg_620[k]
                  - f_474 * lg_625[k]
                  + f_475 * lg_627[k]
                  - f_476 * lg_645[k]
                  + f_477 * lg_650[k]
                  + f_476 * lg_655[k]
                  - f_477 * lg_657[k];
    }

#pragma omp simd aligned(lg_62, lg_67, lg_167, lg_172, lg_197, lg_202, lg_332, lg_337, lg_362, \
                         lg_367, lg_392, lg_397, lg_557, lg_562, lg_587, lg_592, lg_617, \
                         lg_622, lg_647, lg_652 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_70[k] = -f_435 * lg_62[k]
                  + f_434 * lg_67[k]
                  - f_434 * lg_167[k]
                  + f_436 * lg_172[k]
                  + f_438 * lg_197[k]
                  - f_437 * lg_202[k]
                  - f_434 * lg_332[k]
                  + f_436 * lg_337[k]
                  + f_440 * lg_362[k]
                  - f_439 * lg_367[k]
                  - f_442 * lg_392[k]
                  + f_441 * lg_397[k]
                  - f_435 * lg_557[k]
                  + f_434 * lg_562[k]
                  + f_438 * lg_587[k]
                  - f_437 * lg_592[k]
                  - f_442 * lg_617[k]
                  + f_441 * lg_622[k]
                  + f_444 * lg_647[k]
                  - f_443 * lg_652[k];
    }

#pragma omp simd aligned(lg_60, lg_63, lg_70, lg_165, lg_168, lg_175, lg_195, lg_198, lg_205, \
                         lg_330, lg_333, lg_340, lg_360, lg_363, lg_370, lg_390, lg_393, \
                         lg_400, lg_555, lg_558, lg_565, lg_585, lg_588, lg_595, lg_615, \
                         lg_618, lg_625, lg_645, lg_648, lg_655 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_71[k] = -f_478 * lg_60[k]
                  + f_479 * lg_63[k]
                  - f_478 * lg_70[k]
                  - f_480 * lg_165[k]
                  + f_481 * lg_168[k]
                  - f_480 * lg_175[k]
                  + f_482 * lg_195[k]
                  - f_483 * lg_198[k]
                  + f_482 * lg_205[k]
                  - f_480 * lg_330[k]
                  + f_481 * lg_333[k]
                  - f_480 * lg_340[k]
                  + f_484 * lg_360[k]
                  - f_485 * lg_363[k]
                  + f_484 * lg_370[k]
                  - f_486 * lg_390[k]
                  + f_487 * lg_393[k]
                  - f_486 * lg_400[k]
                  - f_478 * lg_555[k]
                  + f_479 * lg_558[k]
                  - f_478 * lg_565[k]
                  + f_482 * lg_585[k]
                  - f_483 * lg_588[k]
                  + f_482 * lg_595[k]
                  - f_486 * lg_615[k]
                  + f_487 * lg_618[k]
                  - f_486 * lg_625[k]
                  + f_488 * lg_645[k]
                  - f_489 * lg_648[k]
                  + f_488 * lg_655[k];
    }

#pragma omp simd aligned(lg_1, lg_6, lg_46, lg_51, lg_76, lg_81, lg_151, lg_156, lg_181, \
                         lg_186, lg_211, lg_216, lg_316, lg_321, lg_346, lg_351, lg_376, \
                         lg_381, lg_406, lg_411, lg_541, lg_546, lg_571, lg_576, lg_601, \
                         lg_606, lg_631, lg_636, lg_661, lg_666 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_72[k] = f_490 * lg_1[k]
                  - f_490 * lg_6[k]
                  + f_491 * lg_46[k]
                  - f_491 * lg_51[k]
                  - f_492 * lg_76[k]
                  + f_492 * lg_81[k]
                  + f_493 * lg_151[k]
                  - f_493 * lg_156[k]
                  - f_430 * lg_181[k]
                  + f_430 * lg_186[k]
                  + f_430 * lg_211[k]
                  - f_430 * lg_216[k]
                  + f_491 * lg_316[k]
                  - f_491 * lg_321[k]
                  - f_430 * lg_346[k]
                  + f_430 * lg_351[k]
                  + f_431 * lg_376[k]
                  - f_431 * lg_381[k]
                  - f_494 * lg_406[k]
                  + f_494 * lg_411[k]
                  + f_490 * lg_541[k]
                  - f_490 * lg_546[k]
                  - f_492 * lg_571[k]
                  + f_492 * lg_576[k]
                  + f_430 * lg_601[k]
                  - f_430 * lg_606[k]
                  - f_494 * lg_631[k]
                  + f_494 * lg_636[k]
                  + f_495 * lg_661[k]
                  - f_495 * lg_666[k];
    }

#pragma omp simd aligned(lg_4, lg_11, lg_49, lg_56, lg_79, lg_86, lg_154, lg_161, lg_184, \
                         lg_191, lg_214, lg_221, lg_319, lg_326, lg_349, lg_356, lg_379, \
                         lg_386, lg_409, lg_416, lg_544, lg_551, lg_574, lg_581, lg_604, \
                         lg_611, lg_634, lg_641, lg_664, lg_671 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_73[k] = f_496 * lg_4[k]
                  - f_497 * lg_11[k]
                  + f_435 * lg_49[k]
                  - f_498 * lg_56[k]
                  - f_438 * lg_79[k]
                  + f_499 * lg_86[k]
                  + f_500 * lg_154[k]
                  - f_501 * lg_161[k]
                  - f_437 * lg_184[k]
                  + f_438 * lg_191[k]
                  + f_437 * lg_214[k]
                  - f_438 * lg_221[k]
                  + f_435 * lg_319[k]
                  - f_498 * lg_326[k]
                  - f_437 * lg_349[k]
                  + f_438 * lg_356[k]
                  + f_439 * lg_379[k]
                  - f_440 * lg_386[k]
                  - f_502 * lg_409[k]
                  + f_503 * lg_416[k]
                  + f_496 * lg_544[k]
                  - f_497 * lg_551[k]
                  - f_438 * lg_574[k]
                  + f_499 * lg_581[k]
                  + f_437 * lg_604[k]
                  - f_438 * lg_611[k]
                  - f_502 * lg_634[k]
                  + f_503 * lg_641[k]
                  + f_504 * lg_664[k]
                  - f_505 * lg_671[k];
    }

#pragma omp simd aligned(lg_1, lg_6, lg_8, lg_46, lg_51, lg_53, lg_76, lg_81, lg_83, lg_151, \
                         lg_156, lg_158, lg_181, lg_186, lg_188, lg_211, lg_216, lg_218, \
                         lg_316, lg_321, lg_323, lg_346, lg_351, lg_353, lg_376, lg_381, \
                         lg_383, lg_406, lg_411, lg_413, lg_541, lg_546, lg_548, lg_571, \
                         lg_576, lg_578, lg_601, lg_606, lg_608, lg_631, lg_636, lg_638, \
                         lg_661, lg_666, lg_668 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_74[k] = -f_506 * lg_1[k]
                  - f_506 * lg_6[k]
                  + f_469 * lg_8[k]
                  - f_507 * lg_46[k]
                  - f_507 * lg_51[k]
                  + f_508 * lg_53[k]
                  + f_509 * lg_76[k]
                  + f_509 * lg_81[k]
                  - f_451 * lg_83[k]
                  - f_469 * lg_151[k]
                  - f_469 * lg_156[k]
                  + f_447 * lg_158[k]
                  + f_449 * lg_181[k]
                  + f_449 * lg_186[k]
                  - f_450 * lg_188[k]
                  - f_449 * lg_211[k]
                  - f_449 * lg_216[k]
                  + f_450 * lg_218[k]
                  - f_507 * lg_316[k]
                  - f_507 * lg_321[k]
                  + f_508 * lg_323[k]
                  + f_449 * lg_346[k]
                  + f_449 * lg_351[k]
                  - f_450 * lg_353[k]
                  - f_451 * lg_376[k]
                  - f_451 * lg_381[k]
                  + f_452 * lg_383[k]
                  + f_510 * lg_406[k]
                  + f_510 * lg_411[k]
                  - f_511 * lg_413[k]
                  - f_506 * lg_541[k]
                  - f_506 * lg_546[k]
                  + f_469 * lg_548[k]
                  + f_509 * lg_571[k]
                  + f_509 * lg_576[k]
                  - f_451 * lg_578[k]
                  - f_449 * lg_601[k]
                  - f_449 * lg_606[k]
                  + f_450 * lg_608[k]
                  + f_510 * lg_631[k]
                  + f_510 * lg_636[k]
                  - f_511 * lg_638[k]
                  - f_512 * lg_661[k]
                  - f_512 * lg_666[k]
                  + f_455 * lg_668[k];
    }

#pragma omp simd aligned(lg_4, lg_11, lg_13, lg_49, lg_56, lg_58, lg_79, lg_86, lg_88, lg_154, \
                         lg_161, lg_163, lg_184, lg_191, lg_193, lg_214, lg_221, lg_223, \
                         lg_319, lg_326, lg_328, lg_349, lg_356, lg_358, lg_379, lg_386, \
                         lg_388, lg_409, lg_416, lg_418, lg_544, lg_551, lg_553, lg_574, \
                         lg_581, lg_583, lg_604, lg_611, lg_613, lg_634, lg_641, lg_643, \
                         lg_664, lg_671, lg_673 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_75[k] = -f_513 * lg_4[k]
                  - f_513 * lg_11[k]
                  + f_514 * lg_13[k]
                  - f_515 * lg_49[k]
                  - f_515 * lg_56[k]
                  + f_516 * lg_58[k]
                  + f_517 * lg_79[k]
                  + f_517 * lg_86[k]
                  - f_518 * lg_88[k]
                  - f_519 * lg_154[k]
                  - f_519 * lg_161[k]
                  + f_520 * lg_163[k]
                  + f_461 * lg_184[k]
                  + f_461 * lg_191[k]
                  - f_462 * lg_193[k]
                  - f_461 * lg_214[k]
                  - f_461 * lg_221[k]
                  + f_462 * lg_223[k]
                  - f_515 * lg_319[k]
                  - f_515 * lg_326[k]
                  + f_516 * lg_328[k]
                  + f_461 * lg_349[k]
                  + f_461 * lg_356[k]
                  - f_462 * lg_358[k]
                  - f_463 * lg_379[k]
                  - f_463 * lg_386[k]
                  + f_464 * lg_388[k]
                  + f_521 * lg_409[k]
                  + f_521 * lg_416[k]
                  - f_522 * lg_418[k]
                  - f_513 * lg_544[k]
                  - f_513 * lg_551[k]
                  + f_514 * lg_553[k]
                  + f_517 * lg_574[k]
                  + f_517 * lg_581[k]
                  - f_518 * lg_583[k]
                  - f_461 * lg_604[k]
                  - f_461 * lg_611[k]
                  + f_462 * lg_613[k]
                  + f_521 * lg_634[k]
                  + f_521 * lg_641[k]
                  - f_522 * lg_643[k]
                  - f_523 * lg_664[k]
                  - f_523 * lg_671[k]
                  + f_524 * lg_673[k];
    }

#pragma omp simd aligned(lg_0, lg_3, lg_5, lg_10, lg_12, lg_14, lg_45, lg_48, lg_50, lg_55, \
                         lg_57, lg_59, lg_75, lg_78, lg_80, lg_85, lg_87, lg_89, lg_150, \
                         lg_153, lg_155, lg_160, lg_162, lg_164, lg_180, lg_183, lg_185, \
                         lg_190, lg_192, lg_194, lg_210, lg_213, lg_215, lg_220, lg_222, \
                         lg_224, lg_315, lg_318, lg_320, lg_325, lg_327, lg_329, lg_345, \
                         lg_348, lg_350, lg_355, lg_357, lg_359, lg_375, lg_378, lg_380, \
                         lg_385, lg_387, lg_389, lg_405, lg_408, lg_410, lg_415, lg_417, \
                         lg_419, lg_540, lg_543, lg_545, lg_550, lg_552, lg_554, lg_570, \
                         lg_573, lg_575, lg_580, lg_582, lg_584, lg_600, lg_603, lg_605, \
                         lg_610, lg_612, lg_614, lg_630, lg_633, lg_635, lg_640, lg_642, \
                         lg_644, lg_660, lg_663, lg_665, lg_670, lg_672, \
                         lg_674 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_76[k] = 0.1025390625 * lg_0[k]
                  + 0.205078125 * lg_3[k]
                  - 0.8203125 * lg_5[k]
                  + 0.1025390625 * lg_10[k]
                  - 0.8203125 * lg_12[k]
                  + 0.2734375 * lg_14[k]
                  + 0.41015625 * lg_45[k]
                  + 0.8203125 * lg_48[k]
                  - 3.28125 * lg_50[k]
                  + 0.41015625 * lg_55[k]
                  - 3.28125 * lg_57[k]
                  + 1.09375 * lg_59[k]
                  - 3.28125 * lg_75[k]
                  - 6.5625 * lg_78[k]
                  + 26.25 * lg_80[k]
                  - 3.28125 * lg_85[k]
                  + 26.25 * lg_87[k]
                  - 8.75 * lg_89[k]
                  + 0.615234375 * lg_150[k]
                  + 1.23046875 * lg_153[k]
                  - 4.921875 * lg_155[k]
                  + 0.615234375 * lg_160[k]
                  - 4.921875 * lg_162[k]
                  + 1.640625 * lg_164[k]
                  - 9.84375 * lg_180[k]
                  - 19.6875 * lg_183[k]
                  + 78.75 * lg_185[k]
                  - 9.84375 * lg_190[k]
                  + 78.75 * lg_192[k]
                  - 26.25 * lg_194[k]
                  + 9.84375 * lg_210[k]
                  + 19.6875 * lg_213[k]
                  - 78.75 * lg_215[k]
                  + 9.84375 * lg_220[k]
                  - 78.75 * lg_222[k]
                  + 26.25 * lg_224[k]
                  + 0.41015625 * lg_315[k]
                  + 0.8203125 * lg_318[k]
                  - 3.28125 * lg_320[k]
                  + 0.41015625 * lg_325[k]
                  - 3.28125 * lg_327[k]
                  + 1.09375 * lg_329[k]
                  - 9.84375 * lg_345[k]
                  - 19.6875 * lg_348[k]
                  + 78.75 * lg_350[k]
                  - 9.84375 * lg_355[k]
                  + 78.75 * lg_357[k]
                  - 26.25 * lg_359[k]
                  + 19.6875 * lg_375[k]
                  + 39.375 * lg_378[k]
                  - 157.5 * lg_380[k]
                  + 19.6875 * lg_385[k]
                  - 157.5 * lg_387[k]
                  + 52.5 * lg_389[k]
                  - 5.25 * lg_405[k]
                  - 10.5 * lg_408[k]
                  + 42.0 * lg_410[k]
                  - 5.25 * lg_415[k]
                  + 42.0 * lg_417[k]
                  - 14.0 * lg_419[k]
                  + 0.1025390625 * lg_540[k]
                  + 0.205078125 * lg_543[k]
                  - 0.8203125 * lg_545[k]
                  + 0.1025390625 * lg_550[k]
                  - 0.8203125 * lg_552[k]
                  + 0.2734375 * lg_554[k]
                  - 3.28125 * lg_570[k]
                  - 6.5625 * lg_573[k]
                  + 26.25 * lg_575[k]
                  - 3.28125 * lg_580[k]
                  + 26.25 * lg_582[k]
                  - 8.75 * lg_584[k]
                  + 9.84375 * lg_600[k]
                  + 19.6875 * lg_603[k]
                  - 78.75 * lg_605[k]
                  + 9.84375 * lg_610[k]
                  - 78.75 * lg_612[k]
                  + 26.25 * lg_614[k]
                  - 5.25 * lg_630[k]
                  - 10.5 * lg_633[k]
                  + 42.0 * lg_635[k]
                  - 5.25 * lg_640[k]
                  + 42.0 * lg_642[k]
                  - 14.0 * lg_644[k]
                  + 0.375 * lg_660[k]
                  + 0.75 * lg_663[k]
                  - 3.0 * lg_665[k]
                  + 0.375 * lg_670[k]
                  - 3.0 * lg_672[k]
                  + lg_674[k];
    }

#pragma omp simd aligned(lg_2, lg_7, lg_9, lg_47, lg_52, lg_54, lg_77, lg_82, lg_84, lg_152, \
                         lg_157, lg_159, lg_182, lg_187, lg_189, lg_212, lg_217, lg_219, \
                         lg_317, lg_322, lg_324, lg_347, lg_352, lg_354, lg_377, lg_382, \
                         lg_384, lg_407, lg_412, lg_414, lg_542, lg_547, lg_549, lg_572, \
                         lg_577, lg_579, lg_602, lg_607, lg_609, lg_632, lg_637, lg_639, \
                         lg_662, lg_667, lg_669 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_77[k] = -f_513 * lg_2[k]
                  - f_513 * lg_7[k]
                  + f_514 * lg_9[k]
                  - f_515 * lg_47[k]
                  - f_515 * lg_52[k]
                  + f_516 * lg_54[k]
                  + f_517 * lg_77[k]
                  + f_517 * lg_82[k]
                  - f_518 * lg_84[k]
                  - f_519 * lg_152[k]
                  - f_519 * lg_157[k]
                  + f_520 * lg_159[k]
                  + f_461 * lg_182[k]
                  + f_461 * lg_187[k]
                  - f_462 * lg_189[k]
                  - f_461 * lg_212[k]
                  - f_461 * lg_217[k]
                  + f_462 * lg_219[k]
                  - f_515 * lg_317[k]
                  - f_515 * lg_322[k]
                  + f_516 * lg_324[k]
                  + f_461 * lg_347[k]
                  + f_461 * lg_352[k]
                  - f_462 * lg_354[k]
                  - f_463 * lg_377[k]
                  - f_463 * lg_382[k]
                  + f_464 * lg_384[k]
                  + f_521 * lg_407[k]
                  + f_521 * lg_412[k]
                  - f_522 * lg_414[k]
                  - f_513 * lg_542[k]
                  - f_513 * lg_547[k]
                  + f_514 * lg_549[k]
                  + f_517 * lg_572[k]
                  + f_517 * lg_577[k]
                  - f_518 * lg_579[k]
                  - f_461 * lg_602[k]
                  - f_461 * lg_607[k]
                  + f_462 * lg_609[k]
                  + f_521 * lg_632[k]
                  + f_521 * lg_637[k]
                  - f_522 * lg_639[k]
                  - f_523 * lg_662[k]
                  - f_523 * lg_667[k]
                  + f_524 * lg_669[k];
    }

#pragma omp simd aligned(lg_0, lg_5, lg_10, lg_12, lg_45, lg_50, lg_55, lg_57, lg_75, lg_80, \
                         lg_85, lg_87, lg_150, lg_155, lg_160, lg_162, lg_180, lg_185, lg_190, \
                         lg_192, lg_210, lg_215, lg_220, lg_222, lg_315, lg_320, lg_325, \
                         lg_327, lg_345, lg_350, lg_355, lg_357, lg_375, lg_380, lg_385, \
                         lg_387, lg_405, lg_410, lg_415, lg_417, lg_540, lg_545, lg_550, \
                         lg_552, lg_570, lg_575, lg_580, lg_582, lg_600, lg_605, lg_610, \
                         lg_612, lg_630, lg_635, lg_640, lg_642, lg_660, lg_665, lg_670, \
                         lg_672 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_78[k] = -f_525 * lg_0[k]
                  + f_526 * lg_5[k]
                  + f_525 * lg_10[k]
                  - f_526 * lg_12[k]
                  - f_527 * lg_45[k]
                  + f_445 * lg_50[k]
                  + f_527 * lg_55[k]
                  - f_445 * lg_57[k]
                  + f_528 * lg_75[k]
                  - f_449 * lg_80[k]
                  - f_528 * lg_85[k]
                  + f_449 * lg_87[k]
                  - f_526 * lg_150[k]
                  + f_470 * lg_155[k]
                  + f_526 * lg_160[k]
                  - f_470 * lg_162[k]
                  + f_472 * lg_180[k]
                  - f_473 * lg_185[k]
                  - f_472 * lg_190[k]
                  + f_473 * lg_192[k]
                  - f_472 * lg_210[k]
                  + f_473 * lg_215[k]
                  + f_472 * lg_220[k]
                  - f_473 * lg_222[k]
                  - f_527 * lg_315[k]
                  + f_445 * lg_320[k]
                  + f_527 * lg_325[k]
                  - f_445 * lg_327[k]
                  + f_472 * lg_345[k]
                  - f_473 * lg_350[k]
                  - f_472 * lg_355[k]
                  + f_473 * lg_357[k]
                  - f_449 * lg_375[k]
                  + f_450 * lg_380[k]
                  + f_449 * lg_385[k]
                  - f_450 * lg_387[k]
                  + f_529 * lg_405[k]
                  - f_530 * lg_410[k]
                  - f_529 * lg_415[k]
                  + f_530 * lg_417[k]
                  - f_525 * lg_540[k]
                  + f_526 * lg_545[k]
                  + f_525 * lg_550[k]
                  - f_526 * lg_552[k]
                  + f_528 * lg_570[k]
                  - f_449 * lg_575[k]
                  - f_528 * lg_580[k]
                  + f_449 * lg_582[k]
                  - f_472 * lg_600[k]
                  + f_473 * lg_605[k]
                  + f_472 * lg_610[k]
                  - f_473 * lg_612[k]
                  + f_529 * lg_630[k]
                  - f_530 * lg_635[k]
                  - f_529 * lg_640[k]
                  + f_530 * lg_642[k]
                  - f_531 * lg_660[k]
                  + f_476 * lg_665[k]
                  + f_531 * lg_670[k]
                  - f_476 * lg_672[k];
    }

#pragma omp simd aligned(lg_2, lg_7, lg_47, lg_52, lg_77, lg_82, lg_152, lg_157, lg_182, \
                         lg_187, lg_212, lg_217, lg_317, lg_322, lg_347, lg_352, lg_377, \
                         lg_382, lg_407, lg_412, lg_542, lg_547, lg_572, lg_577, lg_602, \
                         lg_607, lg_632, lg_637, lg_662, lg_667 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_79[k] = f_497 * lg_2[k]
                  - f_496 * lg_7[k]
                  + f_498 * lg_47[k]
                  - f_435 * lg_52[k]
                  - f_499 * lg_77[k]
                  + f_438 * lg_82[k]
                  + f_501 * lg_152[k]
                  - f_500 * lg_157[k]
                  - f_438 * lg_182[k]
                  + f_437 * lg_187[k]
                  + f_438 * lg_212[k]
                  - f_437 * lg_217[k]
                  + f_498 * lg_317[k]
                  - f_435 * lg_322[k]
                  - f_438 * lg_347[k]
                  + f_437 * lg_352[k]
                  + f_440 * lg_377[k]
                  - f_439 * lg_382[k]
                  - f_503 * lg_407[k]
                  + f_502 * lg_412[k]
                  + f_497 * lg_542[k]
                  - f_496 * lg_547[k]
                  - f_499 * lg_572[k]
                  + f_438 * lg_577[k]
                  + f_438 * lg_602[k]
                  - f_437 * lg_607[k]
                  - f_503 * lg_632[k]
                  + f_502 * lg_637[k]
                  + f_505 * lg_662[k]
                  - f_504 * lg_667[k];
    }

#pragma omp simd aligned(lg_0, lg_3, lg_10, lg_45, lg_48, lg_55, lg_75, lg_78, lg_85, lg_150, \
                         lg_153, lg_160, lg_180, lg_183, lg_190, lg_210, lg_213, lg_220, \
                         lg_315, lg_318, lg_325, lg_345, lg_348, lg_355, lg_375, lg_378, \
                         lg_385, lg_405, lg_408, lg_415, lg_540, lg_543, lg_550, lg_570, \
                         lg_573, lg_580, lg_600, lg_603, lg_610, lg_630, lg_633, lg_640, \
                         lg_660, lg_663, lg_670 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_80[k] = f_532 * lg_0[k]
                  - f_533 * lg_3[k]
                  + f_532 * lg_10[k]
                  + f_490 * lg_45[k]
                  - f_493 * lg_48[k]
                  + f_490 * lg_55[k]
                  - f_534 * lg_75[k]
                  + f_484 * lg_78[k]
                  - f_534 * lg_85[k]
                  + f_533 * lg_150[k]
                  - f_480 * lg_153[k]
                  + f_533 * lg_160[k]
                  - f_482 * lg_180[k]
                  + f_483 * lg_183[k]
                  - f_482 * lg_190[k]
                  + f_482 * lg_210[k]
                  - f_483 * lg_213[k]
                  + f_482 * lg_220[k]
                  + f_490 * lg_315[k]
                  - f_493 * lg_318[k]
                  + f_490 * lg_325[k]
                  - f_482 * lg_345[k]
                  + f_483 * lg_348[k]
                  - f_482 * lg_355[k]
                  + f_484 * lg_375[k]
                  - f_485 * lg_378[k]
                  + f_484 * lg_385[k]
                  - f_535 * lg_405[k]
                  + f_536 * lg_408[k]
                  - f_535 * lg_415[k]
                  + f_532 * lg_540[k]
                  - f_533 * lg_543[k]
                  + f_532 * lg_550[k]
                  - f_534 * lg_570[k]
                  + f_484 * lg_573[k]
                  - f_534 * lg_580[k]
                  + f_482 * lg_600[k]
                  - f_483 * lg_603[k]
                  + f_482 * lg_610[k]
                  - f_535 * lg_630[k]
                  + f_536 * lg_633[k]
                  - f_535 * lg_640[k]
                  + f_537 * lg_660[k]
                  - f_488 * lg_663[k]
                  + f_537 * lg_670[k];
    }

#pragma omp simd aligned(lg_31, lg_36, lg_106, lg_111, lg_136, lg_141, lg_241, lg_246, lg_271, \
                         lg_276, lg_301, lg_306, lg_436, lg_441, lg_466, lg_471, lg_496, \
                         lg_501, lg_526, lg_531 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_81[k] = -f_428 * lg_31[k]
                  + f_428 * lg_36[k]
                  - f_429 * lg_106[k]
                  + f_429 * lg_111[k]
                  + f_430 * lg_136[k]
                  - f_430 * lg_141[k]
                  - f_429 * lg_241[k]
                  + f_429 * lg_246[k]
                  + f_431 * lg_271[k]
                  - f_431 * lg_276[k]
                  - f_432 * lg_301[k]
                  + f_432 * lg_306[k]
                  - f_428 * lg_436[k]
                  + f_428 * lg_441[k]
                  + f_430 * lg_466[k]
                  - f_430 * lg_471[k]
                  - f_432 * lg_496[k]
                  + f_432 * lg_501[k]
                  + f_433 * lg_526[k]
                  - f_433 * lg_531[k];
    }

#pragma omp simd aligned(lg_34, lg_41, lg_109, lg_116, lg_139, lg_146, lg_244, lg_251, lg_274, \
                         lg_281, lg_304, lg_311, lg_439, lg_446, lg_469, lg_476, lg_499, \
                         lg_506, lg_529, lg_536 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_82[k] = -f_434 * lg_34[k]
                  + f_435 * lg_41[k]
                  - f_436 * lg_109[k]
                  + f_434 * lg_116[k]
                  + f_437 * lg_139[k]
                  - f_438 * lg_146[k]
                  - f_436 * lg_244[k]
                  + f_434 * lg_251[k]
                  + f_439 * lg_274[k]
                  - f_440 * lg_281[k]
                  - f_441 * lg_304[k]
                  + f_442 * lg_311[k]
                  - f_434 * lg_439[k]
                  + f_435 * lg_446[k]
                  + f_437 * lg_469[k]
                  - f_438 * lg_476[k]
                  - f_441 * lg_499[k]
                  + f_442 * lg_506[k]
                  + f_443 * lg_529[k]
                  - f_444 * lg_536[k];
    }

#pragma omp simd aligned(lg_31, lg_36, lg_38, lg_106, lg_111, lg_113, lg_136, lg_141, lg_143, \
                         lg_241, lg_246, lg_248, lg_271, lg_276, lg_278, lg_301, lg_306, \
                         lg_308, lg_436, lg_441, lg_443, lg_466, lg_471, lg_473, lg_496, \
                         lg_501, lg_503, lg_526, lg_531, lg_533 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_83[k] = f_445 * lg_31[k]
                  + f_445 * lg_36[k]
                  - f_446 * lg_38[k]
                  + f_447 * lg_106[k]
                  + f_447 * lg_111[k]
                  - f_448 * lg_113[k]
                  - f_449 * lg_136[k]
                  - f_449 * lg_141[k]
                  + f_450 * lg_143[k]
                  + f_447 * lg_241[k]
                  + f_447 * lg_246[k]
                  - f_448 * lg_248[k]
                  - f_451 * lg_271[k]
                  - f_451 * lg_276[k]
                  + f_452 * lg_278[k]
                  + f_453 * lg_301[k]
                  + f_453 * lg_306[k]
                  - f_454 * lg_308[k]
                  + f_445 * lg_436[k]
                  + f_445 * lg_441[k]
                  - f_446 * lg_443[k]
                  - f_449 * lg_466[k]
                  - f_449 * lg_471[k]
                  + f_450 * lg_473[k]
                  + f_453 * lg_496[k]
                  + f_453 * lg_501[k]
                  - f_454 * lg_503[k]
                  - f_455 * lg_526[k]
                  - f_455 * lg_531[k]
                  + f_456 * lg_533[k];
    }

#pragma omp simd aligned(lg_34, lg_41, lg_43, lg_109, lg_116, lg_118, lg_139, lg_146, lg_148, \
                         lg_244, lg_251, lg_253, lg_274, lg_281, lg_283, lg_304, lg_311, \
                         lg_313, lg_439, lg_446, lg_448, lg_469, lg_476, lg_478, lg_499, \
                         lg_506, lg_508, lg_529, lg_536, lg_538 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_84[k] = f_457 * lg_34[k]
                  + f_457 * lg_41[k]
                  - f_458 * lg_43[k]
                  + f_459 * lg_109[k]
                  + f_459 * lg_116[k]
                  - f_460 * lg_118[k]
                  - f_461 * lg_139[k]
                  - f_461 * lg_146[k]
                  + f_462 * lg_148[k]
                  + f_459 * lg_244[k]
                  + f_459 * lg_251[k]
                  - f_460 * lg_253[k]
                  - f_463 * lg_274[k]
                  - f_463 * lg_281[k]
                  + f_464 * lg_283[k]
                  + f_465 * lg_304[k]
                  + f_465 * lg_311[k]
                  - f_466 * lg_313[k]
                  + f_457 * lg_439[k]
                  + f_457 * lg_446[k]
                  - f_458 * lg_448[k]
                  - f_461 * lg_469[k]
                  - f_461 * lg_476[k]
                  + f_462 * lg_478[k]
                  + f_465 * lg_499[k]
                  + f_465 * lg_506[k]
                  - f_466 * lg_508[k]
                  - f_467 * lg_529[k]
                  - f_467 * lg_536[k]
                  + f_468 * lg_538[k];
    }

#pragma omp simd aligned(lg_30, lg_33, lg_35, lg_40, lg_42, lg_44, lg_105, lg_108, lg_110, \
                         lg_115, lg_117, lg_119, lg_135, lg_138, lg_140, lg_145, lg_147, \
                         lg_149, lg_240, lg_243, lg_245, lg_250, lg_252, lg_254, lg_270, \
                         lg_273, lg_275, lg_280, lg_282, lg_284, lg_300, lg_303, lg_305, \
                         lg_310, lg_312, lg_314, lg_435, lg_438, lg_440, lg_445, lg_447, \
                         lg_449, lg_465, lg_468, lg_470, lg_475, lg_477, lg_479, lg_495, \
                         lg_498, lg_500, lg_505, lg_507, lg_509, lg_525, lg_528, lg_530, \
                         lg_535, lg_537, lg_539 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_85[k] = -1.23046875 * lg_30[k]
                  - 2.4609375 * lg_33[k]
                  + 9.84375 * lg_35[k]
                  - 1.23046875 * lg_40[k]
                  + 9.84375 * lg_42[k]
                  - 3.28125 * lg_44[k]
                  - 3.69140625 * lg_105[k]
                  - 7.3828125 * lg_108[k]
                  + 29.53125 * lg_110[k]
                  - 3.69140625 * lg_115[k]
                  + 29.53125 * lg_117[k]
                  - 9.84375 * lg_119[k]
                  + 9.84375 * lg_135[k]
                  + 19.6875 * lg_138[k]
                  - 78.75 * lg_140[k]
                  + 9.84375 * lg_145[k]
                  - 78.75 * lg_147[k]
                  + 26.25 * lg_149[k]
                  - 3.69140625 * lg_240[k]
                  - 7.3828125 * lg_243[k]
                  + 29.53125 * lg_245[k]
                  - 3.69140625 * lg_250[k]
                  + 29.53125 * lg_252[k]
                  - 9.84375 * lg_254[k]
                  + 19.6875 * lg_270[k]
                  + 39.375 * lg_273[k]
                  - 157.5 * lg_275[k]
                  + 19.6875 * lg_280[k]
                  - 157.5 * lg_282[k]
                  + 52.5 * lg_284[k]
                  - 11.8125 * lg_300[k]
                  - 23.625 * lg_303[k]
                  + 94.5 * lg_305[k]
                  - 11.8125 * lg_310[k]
                  + 94.5 * lg_312[k]
                  - 31.5 * lg_314[k]
                  - 1.23046875 * lg_435[k]
                  - 2.4609375 * lg_438[k]
                  + 9.84375 * lg_440[k]
                  - 1.23046875 * lg_445[k]
                  + 9.84375 * lg_447[k]
                  - 3.28125 * lg_449[k]
                  + 9.84375 * lg_465[k]
                  + 19.6875 * lg_468[k]
                  - 78.75 * lg_470[k]
                  + 9.84375 * lg_475[k]
                  - 78.75 * lg_477[k]
                  + 26.25 * lg_479[k]
                  - 11.8125 * lg_495[k]
                  - 23.625 * lg_498[k]
                  + 94.5 * lg_500[k]
                  - 11.8125 * lg_505[k]
                  + 94.5 * lg_507[k]
                  - 31.5 * lg_509[k]
                  + 2.25 * lg_525[k]
                  + 4.5 * lg_528[k]
                  - 18.0 * lg_530[k]
                  + 2.25 * lg_535[k]
                  - 18.0 * lg_537[k]
                  + 6.0 * lg_539[k];
    }

#pragma omp simd aligned(lg_32, lg_37, lg_39, lg_107, lg_112, lg_114, lg_137, lg_142, lg_144, \
                         lg_242, lg_247, lg_249, lg_272, lg_277, lg_279, lg_302, lg_307, \
                         lg_309, lg_437, lg_442, lg_444, lg_467, lg_472, lg_474, lg_497, \
                         lg_502, lg_504, lg_527, lg_532, lg_534 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_86[k] = f_457 * lg_32[k]
                  + f_457 * lg_37[k]
                  - f_458 * lg_39[k]
                  + f_459 * lg_107[k]
                  + f_459 * lg_112[k]
                  - f_460 * lg_114[k]
                  - f_461 * lg_137[k]
                  - f_461 * lg_142[k]
                  + f_462 * lg_144[k]
                  + f_459 * lg_242[k]
                  + f_459 * lg_247[k]
                  - f_460 * lg_249[k]
                  - f_463 * lg_272[k]
                  - f_463 * lg_277[k]
                  + f_464 * lg_279[k]
                  + f_465 * lg_302[k]
                  + f_465 * lg_307[k]
                  - f_466 * lg_309[k]
                  + f_457 * lg_437[k]
                  + f_457 * lg_442[k]
                  - f_458 * lg_444[k]
                  - f_461 * lg_467[k]
                  - f_461 * lg_472[k]
                  + f_462 * lg_474[k]
                  + f_465 * lg_497[k]
                  + f_465 * lg_502[k]
                  - f_466 * lg_504[k]
                  - f_467 * lg_527[k]
                  - f_467 * lg_532[k]
                  + f_468 * lg_534[k];
    }

#pragma omp simd aligned(lg_30, lg_35, lg_40, lg_42, lg_105, lg_110, lg_115, lg_117, lg_135, \
                         lg_140, lg_145, lg_147, lg_240, lg_245, lg_250, lg_252, lg_270, \
                         lg_275, lg_280, lg_282, lg_300, lg_305, lg_310, lg_312, lg_435, \
                         lg_440, lg_445, lg_447, lg_465, lg_470, lg_475, lg_477, lg_495, \
                         lg_500, lg_505, lg_507, lg_525, lg_530, lg_535, \
                         lg_537 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_87[k] = f_469 * lg_30[k]
                  - f_447 * lg_35[k]
                  - f_469 * lg_40[k]
                  + f_447 * lg_42[k]
                  + f_470 * lg_105[k]
                  - f_471 * lg_110[k]
                  - f_470 * lg_115[k]
                  + f_471 * lg_117[k]
                  - f_472 * lg_135[k]
                  + f_473 * lg_140[k]
                  + f_472 * lg_145[k]
                  - f_473 * lg_147[k]
                  + f_470 * lg_240[k]
                  - f_471 * lg_245[k]
                  - f_470 * lg_250[k]
                  + f_471 * lg_252[k]
                  - f_449 * lg_270[k]
                  + f_450 * lg_275[k]
                  + f_449 * lg_280[k]
                  - f_450 * lg_282[k]
                  + f_474 * lg_300[k]
                  - f_475 * lg_305[k]
                  - f_474 * lg_310[k]
                  + f_475 * lg_312[k]
                  + f_469 * lg_435[k]
                  - f_447 * lg_440[k]
                  - f_469 * lg_445[k]
                  + f_447 * lg_447[k]
                  - f_472 * lg_465[k]
                  + f_473 * lg_470[k]
                  + f_472 * lg_475[k]
                  - f_473 * lg_477[k]
                  + f_474 * lg_495[k]
                  - f_475 * lg_500[k]
                  - f_474 * lg_505[k]
                  + f_475 * lg_507[k]
                  - f_476 * lg_525[k]
                  + f_477 * lg_530[k]
                  + f_476 * lg_535[k]
                  - f_477 * lg_537[k];
    }

#pragma omp simd aligned(lg_32, lg_37, lg_107, lg_112, lg_137, lg_142, lg_242, lg_247, lg_272, \
                         lg_277, lg_302, lg_307, lg_437, lg_442, lg_467, lg_472, lg_497, \
                         lg_502, lg_527, lg_532 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_88[k] = -f_435 * lg_32[k]
                  + f_434 * lg_37[k]
                  - f_434 * lg_107[k]
                  + f_436 * lg_112[k]
                  + f_438 * lg_137[k]
                  - f_437 * lg_142[k]
                  - f_434 * lg_242[k]
                  + f_436 * lg_247[k]
                  + f_440 * lg_272[k]
                  - f_439 * lg_277[k]
                  - f_442 * lg_302[k]
                  + f_441 * lg_307[k]
                  - f_435 * lg_437[k]
                  + f_434 * lg_442[k]
                  + f_438 * lg_467[k]
                  - f_437 * lg_472[k]
                  - f_442 * lg_497[k]
                  + f_441 * lg_502[k]
                  + f_444 * lg_527[k]
                  - f_443 * lg_532[k];
    }

#pragma omp simd aligned(lg_30, lg_33, lg_40, lg_105, lg_108, lg_115, lg_135, lg_138, lg_145, \
                         lg_240, lg_243, lg_250, lg_270, lg_273, lg_280, lg_300, lg_303, \
                         lg_310, lg_435, lg_438, lg_445, lg_465, lg_468, lg_475, lg_495, \
                         lg_498, lg_505, lg_525, lg_528, lg_535 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_89[k] = -f_478 * lg_30[k]
                  + f_479 * lg_33[k]
                  - f_478 * lg_40[k]
                  - f_480 * lg_105[k]
                  + f_481 * lg_108[k]
                  - f_480 * lg_115[k]
                  + f_482 * lg_135[k]
                  - f_483 * lg_138[k]
                  + f_482 * lg_145[k]
                  - f_480 * lg_240[k]
                  + f_481 * lg_243[k]
                  - f_480 * lg_250[k]
                  + f_484 * lg_270[k]
                  - f_485 * lg_273[k]
                  + f_484 * lg_280[k]
                  - f_486 * lg_300[k]
                  + f_487 * lg_303[k]
                  - f_486 * lg_310[k]
                  - f_478 * lg_435[k]
                  + f_479 * lg_438[k]
                  - f_478 * lg_445[k]
                  + f_482 * lg_465[k]
                  - f_483 * lg_468[k]
                  + f_482 * lg_475[k]
                  - f_486 * lg_495[k]
                  + f_487 * lg_498[k]
                  - f_486 * lg_505[k]
                  + f_488 * lg_525[k]
                  - f_489 * lg_528[k]
                  + f_488 * lg_535[k];
    }

#pragma omp simd aligned(lg_1, lg_6, lg_46, lg_51, lg_76, lg_81, lg_181, lg_186, lg_211, \
                         lg_216, lg_316, lg_321, lg_346, lg_351, lg_406, lg_411, lg_541, \
                         lg_546, lg_571, lg_576, lg_601, lg_606, lg_631, \
                         lg_636 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_90[k] = -f_538 * lg_1[k]
                  + f_538 * lg_6[k]
                  - f_358 * lg_46[k]
                  + f_358 * lg_51[k]
                  + f_422 * lg_76[k]
                  - f_422 * lg_81[k]
                  + f_422 * lg_181[k]
                  - f_422 * lg_186[k]
                  - f_539 * lg_211[k]
                  + f_539 * lg_216[k]
                  + f_358 * lg_316[k]
                  - f_358 * lg_321[k]
                  - f_422 * lg_346[k]
                  + f_422 * lg_351[k]
                  + f_540 * lg_406[k]
                  - f_540 * lg_411[k]
                  + f_538 * lg_541[k]
                  - f_538 * lg_546[k]
                  - f_422 * lg_571[k]
                  + f_422 * lg_576[k]
                  + f_539 * lg_601[k]
                  - f_539 * lg_606[k]
                  - f_540 * lg_631[k]
                  + f_540 * lg_636[k];
    }

#pragma omp simd aligned(lg_4, lg_11, lg_49, lg_56, lg_79, lg_86, lg_184, lg_191, lg_214, \
                         lg_221, lg_319, lg_326, lg_349, lg_356, lg_409, lg_416, lg_544, \
                         lg_551, lg_574, lg_581, lg_604, lg_611, lg_634, \
                         lg_641 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_91[k] = -2.4609375 * lg_4[k]
                  + 0.8203125 * lg_11[k]
                  - 4.921875 * lg_49[k]
                  + 1.640625 * lg_56[k]
                  + 73.828125 * lg_79[k]
                  - 24.609375 * lg_86[k]
                  + 73.828125 * lg_184[k]
                  - 24.609375 * lg_191[k]
                  - 196.875 * lg_214[k]
                  + 65.625 * lg_221[k]
                  + 4.921875 * lg_319[k]
                  - 1.640625 * lg_326[k]
                  - 73.828125 * lg_349[k]
                  + 24.609375 * lg_356[k]
                  + 78.75 * lg_409[k]
                  - 26.25 * lg_416[k]
                  + 2.4609375 * lg_544[k]
                  - 0.8203125 * lg_551[k]
                  - 73.828125 * lg_574[k]
                  + 24.609375 * lg_581[k]
                  + 196.875 * lg_604[k]
                  - 65.625 * lg_611[k]
                  - 78.75 * lg_634[k]
                  + 26.25 * lg_641[k];
    }

#pragma omp simd aligned(lg_1, lg_6, lg_8, lg_46, lg_51, lg_53, lg_76, lg_81, lg_83, lg_181, \
                         lg_186, lg_188, lg_211, lg_216, lg_218, lg_316, lg_321, lg_323, \
                         lg_346, lg_351, lg_353, lg_406, lg_411, lg_413, lg_541, lg_546, \
                         lg_548, lg_571, lg_576, lg_578, lg_601, lg_606, lg_608, lg_631, \
                         lg_636, lg_638 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_92[k] = f_407 * lg_1[k]
                  + f_407 * lg_6[k]
                  - f_366 * lg_8[k]
                  + f_364 * lg_46[k]
                  + f_364 * lg_51[k]
                  - f_365 * lg_53[k]
                  - f_410 * lg_76[k]
                  - f_410 * lg_81[k]
                  + f_411 * lg_83[k]
                  - f_410 * lg_181[k]
                  - f_410 * lg_186[k]
                  + f_411 * lg_188[k]
                  + f_412 * lg_211[k]
                  + f_412 * lg_216[k]
                  - f_413 * lg_218[k]
                  - f_364 * lg_316[k]
                  - f_364 * lg_321[k]
                  + f_365 * lg_323[k]
                  + f_410 * lg_346[k]
                  + f_410 * lg_351[k]
                  - f_411 * lg_353[k]
                  - f_414 * lg_406[k]
                  - f_414 * lg_411[k]
                  + f_415 * lg_413[k]
                  - f_407 * lg_541[k]
                  - f_407 * lg_546[k]
                  + f_366 * lg_548[k]
                  + f_410 * lg_571[k]
                  + f_410 * lg_576[k]
                  - f_411 * lg_578[k]
                  - f_412 * lg_601[k]
                  - f_412 * lg_606[k]
                  + f_413 * lg_608[k]
                  + f_414 * lg_631[k]
                  + f_414 * lg_636[k]
                  - f_415 * lg_638[k];
    }

#pragma omp simd aligned(lg_4, lg_11, lg_13, lg_49, lg_56, lg_58, lg_79, lg_86, lg_88, lg_184, \
                         lg_191, lg_193, lg_214, lg_221, lg_223, lg_319, lg_326, lg_328, \
                         lg_349, lg_356, lg_358, lg_409, lg_416, lg_418, lg_544, lg_551, \
                         lg_553, lg_574, lg_581, lg_583, lg_604, lg_611, lg_613, lg_634, \
                         lg_641, lg_643 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_93[k] = f_541 * lg_4[k]
                  + f_541 * lg_11[k]
                  - f_542 * lg_13[k]
                  + f_376 * lg_49[k]
                  + f_376 * lg_56[k]
                  - f_377 * lg_58[k]
                  - f_543 * lg_79[k]
                  - f_543 * lg_86[k]
                  + f_544 * lg_88[k]
                  - f_543 * lg_184[k]
                  - f_543 * lg_191[k]
                  + f_544 * lg_193[k]
                  + f_381 * lg_214[k]
                  + f_381 * lg_221[k]
                  - f_545 * lg_223[k]
                  - f_376 * lg_319[k]
                  - f_376 * lg_326[k]
                  + f_377 * lg_328[k]
                  + f_543 * lg_349[k]
                  + f_543 * lg_356[k]
                  - f_544 * lg_358[k]
                  - f_546 * lg_409[k]
                  - f_546 * lg_416[k]
                  + f_547 * lg_418[k]
                  - f_541 * lg_544[k]
                  - f_541 * lg_551[k]
                  + f_542 * lg_553[k]
                  + f_543 * lg_574[k]
                  + f_543 * lg_581[k]
                  - f_544 * lg_583[k]
                  - f_381 * lg_604[k]
                  - f_381 * lg_611[k]
                  + f_545 * lg_613[k]
                  + f_546 * lg_634[k]
                  + f_546 * lg_641[k]
                  - f_547 * lg_643[k];
    }

#pragma omp simd aligned(lg_0, lg_3, lg_5, lg_10, lg_12, lg_14, lg_45, lg_48, lg_50, lg_55, \
                         lg_57, lg_59, lg_75, lg_78, lg_80, lg_85, lg_87, lg_89, lg_180, \
                         lg_183, lg_185, lg_190, lg_192, lg_194, lg_210, lg_213, lg_215, \
                         lg_220, lg_222, lg_224, lg_315, lg_318, lg_320, lg_325, lg_327, \
                         lg_329, lg_345, lg_348, lg_350, lg_355, lg_357, lg_359, lg_405, \
                         lg_408, lg_410, lg_415, lg_417, lg_419, lg_540, lg_543, lg_545, \
                         lg_550, lg_552, lg_554, lg_570, lg_573, lg_575, lg_580, lg_582, \
                         lg_584, lg_600, lg_603, lg_605, lg_610, lg_612, lg_614, lg_630, \
                         lg_633, lg_635, lg_640, lg_642, lg_644 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_94[k] = -f_548 * lg_0[k]
                  - f_387 * lg_3[k]
                  + f_549 * lg_5[k]
                  - f_548 * lg_10[k]
                  + f_549 * lg_12[k]
                  - f_550 * lg_14[k]
                  - f_387 * lg_45[k]
                  - f_388 * lg_48[k]
                  + f_389 * lg_50[k]
                  - f_387 * lg_55[k]
                  + f_389 * lg_57[k]
                  - f_390 * lg_59[k]
                  + f_551 * lg_75[k]
                  + f_394 * lg_78[k]
                  - f_398 * lg_80[k]
                  + f_551 * lg_85[k]
                  - f_398 * lg_87[k]
                  + f_552 * lg_89[k]
                  + f_551 * lg_180[k]
                  + f_394 * lg_183[k]
                  - f_398 * lg_185[k]
                  + f_551 * lg_190[k]
                  - f_398 * lg_192[k]
                  + f_552 * lg_194[k]
                  - f_552 * lg_210[k]
                  - f_397 * lg_213[k]
                  + f_553 * lg_215[k]
                  - f_552 * lg_220[k]
                  + f_553 * lg_222[k]
                  - f_554 * lg_224[k]
                  + f_387 * lg_315[k]
                  + f_388 * lg_318[k]
                  - f_389 * lg_320[k]
                  + f_387 * lg_325[k]
                  - f_389 * lg_327[k]
                  + f_390 * lg_329[k]
                  - f_551 * lg_345[k]
                  - f_394 * lg_348[k]
                  + f_398 * lg_350[k]
                  - f_551 * lg_355[k]
                  + f_398 * lg_357[k]
                  - f_552 * lg_359[k]
                  + f_555 * lg_405[k]
                  + f_403 * lg_408[k]
                  - f_443 * lg_410[k]
                  + f_555 * lg_415[k]
                  - f_443 * lg_417[k]
                  + f_444 * lg_419[k]
                  + f_548 * lg_540[k]
                  + f_387 * lg_543[k]
                  - f_549 * lg_545[k]
                  + f_548 * lg_550[k]
                  - f_549 * lg_552[k]
                  + f_550 * lg_554[k]
                  - f_551 * lg_570[k]
                  - f_394 * lg_573[k]
                  + f_398 * lg_575[k]
                  - f_551 * lg_580[k]
                  + f_398 * lg_582[k]
                  - f_552 * lg_584[k]
                  + f_552 * lg_600[k]
                  + f_397 * lg_603[k]
                  - f_553 * lg_605[k]
                  + f_552 * lg_610[k]
                  - f_553 * lg_612[k]
                  + f_554 * lg_614[k]
                  - f_555 * lg_630[k]
                  - f_403 * lg_633[k]
                  + f_443 * lg_635[k]
                  - f_555 * lg_640[k]
                  + f_443 * lg_642[k]
                  - f_444 * lg_644[k];
    }

#pragma omp simd aligned(lg_2, lg_7, lg_9, lg_47, lg_52, lg_54, lg_77, lg_82, lg_84, lg_182, \
                         lg_187, lg_189, lg_212, lg_217, lg_219, lg_317, lg_322, lg_324, \
                         lg_347, lg_352, lg_354, lg_407, lg_412, lg_414, lg_542, lg_547, \
                         lg_549, lg_572, lg_577, lg_579, lg_602, lg_607, lg_609, lg_632, \
                         lg_637, lg_639 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_95[k] = f_541 * lg_2[k]
                  + f_541 * lg_7[k]
                  - f_542 * lg_9[k]
                  + f_376 * lg_47[k]
                  + f_376 * lg_52[k]
                  - f_377 * lg_54[k]
                  - f_543 * lg_77[k]
                  - f_543 * lg_82[k]
                  + f_544 * lg_84[k]
                  - f_543 * lg_182[k]
                  - f_543 * lg_187[k]
                  + f_544 * lg_189[k]
                  + f_381 * lg_212[k]
                  + f_381 * lg_217[k]
                  - f_545 * lg_219[k]
                  - f_376 * lg_317[k]
                  - f_376 * lg_322[k]
                  + f_377 * lg_324[k]
                  + f_543 * lg_347[k]
                  + f_543 * lg_352[k]
                  - f_544 * lg_354[k]
                  - f_546 * lg_407[k]
                  - f_546 * lg_412[k]
                  + f_547 * lg_414[k]
                  - f_541 * lg_542[k]
                  - f_541 * lg_547[k]
                  + f_542 * lg_549[k]
                  + f_543 * lg_572[k]
                  + f_543 * lg_577[k]
                  - f_544 * lg_579[k]
                  - f_381 * lg_602[k]
                  - f_381 * lg_607[k]
                  + f_545 * lg_609[k]
                  + f_546 * lg_632[k]
                  + f_546 * lg_637[k]
                  - f_547 * lg_639[k];
    }

#pragma omp simd aligned(lg_0, lg_5, lg_10, lg_12, lg_45, lg_50, lg_55, lg_57, lg_75, lg_80, \
                         lg_85, lg_87, lg_180, lg_185, lg_190, lg_192, lg_210, lg_215, lg_220, \
                         lg_222, lg_315, lg_320, lg_325, lg_327, lg_345, lg_350, lg_355, \
                         lg_357, lg_405, lg_410, lg_415, lg_417, lg_540, lg_545, lg_550, \
                         lg_552, lg_570, lg_575, lg_580, lg_582, lg_600, lg_605, lg_610, \
                         lg_612, lg_630, lg_635, lg_640, lg_642 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_96[k] = f_556 * lg_0[k]
                  - f_408 * lg_5[k]
                  - f_556 * lg_10[k]
                  + f_408 * lg_12[k]
                  + f_407 * lg_45[k]
                  - f_366 * lg_50[k]
                  - f_407 * lg_55[k]
                  + f_366 * lg_57[k]
                  - f_557 * lg_75[k]
                  + f_558 * lg_80[k]
                  + f_557 * lg_85[k]
                  - f_558 * lg_87[k]
                  - f_557 * lg_180[k]
                  + f_558 * lg_185[k]
                  + f_557 * lg_190[k]
                  - f_558 * lg_192[k]
                  + f_559 * lg_210[k]
                  - f_560 * lg_215[k]
                  - f_559 * lg_220[k]
                  + f_560 * lg_222[k]
                  - f_407 * lg_315[k]
                  + f_366 * lg_320[k]
                  + f_407 * lg_325[k]
                  - f_366 * lg_327[k]
                  + f_557 * lg_345[k]
                  - f_558 * lg_350[k]
                  - f_557 * lg_355[k]
                  + f_558 * lg_357[k]
                  - f_561 * lg_405[k]
                  + f_562 * lg_410[k]
                  + f_561 * lg_415[k]
                  - f_562 * lg_417[k]
                  - f_556 * lg_540[k]
                  + f_408 * lg_545[k]
                  + f_556 * lg_550[k]
                  - f_408 * lg_552[k]
                  + f_557 * lg_570[k]
                  - f_558 * lg_575[k]
                  - f_557 * lg_580[k]
                  + f_558 * lg_582[k]
                  - f_559 * lg_600[k]
                  + f_560 * lg_605[k]
                  + f_559 * lg_610[k]
                  - f_560 * lg_612[k]
                  + f_561 * lg_630[k]
                  - f_562 * lg_635[k]
                  - f_561 * lg_640[k]
                  + f_562 * lg_642[k];
    }

#pragma omp simd aligned(lg_2, lg_7, lg_47, lg_52, lg_77, lg_82, lg_182, lg_187, lg_212, \
                         lg_217, lg_317, lg_322, lg_347, lg_352, lg_407, lg_412, lg_542, \
                         lg_547, lg_572, lg_577, lg_602, lg_607, lg_632, \
                         lg_637 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_97[k] = -0.8203125 * lg_2[k]
                  + 2.4609375 * lg_7[k]
                  - 1.640625 * lg_47[k]
                  + 4.921875 * lg_52[k]
                  + 24.609375 * lg_77[k]
                  - 73.828125 * lg_82[k]
                  + 24.609375 * lg_182[k]
                  - 73.828125 * lg_187[k]
                  - 65.625 * lg_212[k]
                  + 196.875 * lg_217[k]
                  + 1.640625 * lg_317[k]
                  - 4.921875 * lg_322[k]
                  - 24.609375 * lg_347[k]
                  + 73.828125 * lg_352[k]
                  + 26.25 * lg_407[k]
                  - 78.75 * lg_412[k]
                  + 0.8203125 * lg_542[k]
                  - 2.4609375 * lg_547[k]
                  - 24.609375 * lg_572[k]
                  + 73.828125 * lg_577[k]
                  + 65.625 * lg_602[k]
                  - 196.875 * lg_607[k]
                  - 26.25 * lg_632[k]
                  + 78.75 * lg_637[k];
    }

#pragma omp simd aligned(lg_0, lg_3, lg_10, lg_45, lg_48, lg_55, lg_75, lg_78, lg_85, lg_180, \
                         lg_183, lg_190, lg_210, lg_213, lg_220, lg_315, lg_318, lg_325, \
                         lg_345, lg_348, lg_355, lg_405, lg_408, lg_415, lg_540, lg_543, \
                         lg_550, lg_570, lg_573, lg_580, lg_600, lg_603, lg_610, lg_630, \
                         lg_633, lg_640 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_98[k] = -f_563 * lg_0[k]
                  + f_418 * lg_3[k]
                  - f_563 * lg_10[k]
                  - f_416 * lg_45[k]
                  + f_417 * lg_48[k]
                  - f_416 * lg_55[k]
                  + f_564 * lg_75[k]
                  - f_565 * lg_78[k]
                  + f_564 * lg_85[k]
                  + f_564 * lg_180[k]
                  - f_565 * lg_183[k]
                  + f_564 * lg_190[k]
                  - f_566 * lg_210[k]
                  + f_361 * lg_213[k]
                  - f_566 * lg_220[k]
                  + f_416 * lg_315[k]
                  - f_417 * lg_318[k]
                  + f_416 * lg_325[k]
                  - f_564 * lg_345[k]
                  + f_565 * lg_348[k]
                  - f_564 * lg_355[k]
                  + f_567 * lg_405[k]
                  - f_568 * lg_408[k]
                  + f_567 * lg_415[k]
                  + f_563 * lg_540[k]
                  - f_418 * lg_543[k]
                  + f_563 * lg_550[k]
                  - f_564 * lg_570[k]
                  + f_565 * lg_573[k]
                  - f_564 * lg_580[k]
                  + f_566 * lg_600[k]
                  - f_361 * lg_603[k]
                  + f_566 * lg_610[k]
                  - f_567 * lg_630[k]
                  + f_568 * lg_633[k]
                  - f_567 * lg_640[k];
    }

#pragma omp simd aligned(lg_31, lg_36, lg_106, lg_111, lg_136, lg_141, lg_241, lg_246, lg_271, \
                         lg_276, lg_301, lg_306, lg_436, lg_441, lg_466, lg_471, lg_496, \
                         lg_501 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_99[k] = f_258 * lg_31[k]
                  - f_258 * lg_36[k]
                  - f_258 * lg_106[k]
                  + f_258 * lg_111[k]
                  - f_261 * lg_136[k]
                  + f_261 * lg_141[k]
                  - f_256 * lg_241[k]
                  + f_256 * lg_246[k]
                  + f_259 * lg_271[k]
                  - f_259 * lg_276[k]
                  + f_262 * lg_301[k]
                  - f_262 * lg_306[k]
                  - f_255 * lg_436[k]
                  + f_255 * lg_441[k]
                  + f_257 * lg_466[k]
                  - f_257 * lg_471[k]
                  - f_260 * lg_496[k]
                  + f_260 * lg_501[k];
    }

#pragma omp simd aligned(lg_34, lg_41, lg_109, lg_116, lg_139, lg_146, lg_244, lg_251, lg_274, \
                         lg_281, lg_304, lg_311, lg_439, lg_446, lg_469, lg_476, lg_499, \
                         lg_506 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_100[k] = f_264 * lg_34[k]
                   - f_269 * lg_41[k]
                   - f_264 * lg_109[k]
                   + f_269 * lg_116[k]
                   - f_268 * lg_139[k]
                   + f_274 * lg_146[k]
                   - f_265 * lg_244[k]
                   + f_266 * lg_251[k]
                   + f_270 * lg_274[k]
                   - f_271 * lg_281[k]
                   + f_273 * lg_304[k]
                   - f_275 * lg_311[k]
                   - f_263 * lg_439[k]
                   + f_264 * lg_446[k]
                   + f_267 * lg_469[k]
                   - f_268 * lg_476[k]
                   - f_272 * lg_499[k]
                   + f_273 * lg_506[k];
    }

#pragma omp simd aligned(lg_31, lg_36, lg_38, lg_106, lg_111, lg_113, lg_136, lg_141, lg_143, \
                         lg_241, lg_246, lg_248, lg_271, lg_276, lg_278, lg_301, lg_306, \
                         lg_308, lg_436, lg_441, lg_443, lg_466, lg_471, lg_473, lg_496, \
                         lg_501, lg_503 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_101[k] = -f_282 * lg_31[k]
                   - f_282 * lg_36[k]
                   + f_283 * lg_38[k]
                   + f_282 * lg_106[k]
                   + f_282 * lg_111[k]
                   - f_283 * lg_113[k]
                   + f_288 * lg_136[k]
                   + f_288 * lg_141[k]
                   - f_289 * lg_143[k]
                   + f_278 * lg_241[k]
                   + f_278 * lg_246[k]
                   - f_279 * lg_248[k]
                   - f_284 * lg_271[k]
                   - f_284 * lg_276[k]
                   + f_285 * lg_278[k]
                   - f_290 * lg_301[k]
                   - f_290 * lg_306[k]
                   + f_291 * lg_308[k]
                   + f_276 * lg_436[k]
                   + f_276 * lg_441[k]
                   - f_277 * lg_443[k]
                   - f_280 * lg_466[k]
                   - f_280 * lg_471[k]
                   + f_281 * lg_473[k]
                   + f_286 * lg_496[k]
                   + f_286 * lg_501[k]
                   - f_287 * lg_503[k];
    }

#pragma omp simd aligned(lg_34, lg_41, lg_43, lg_109, lg_116, lg_118, lg_139, lg_146, lg_148, \
                         lg_244, lg_251, lg_253, lg_274, lg_281, lg_283, lg_304, lg_311, \
                         lg_313, lg_439, lg_446, lg_448, lg_469, lg_476, lg_478, lg_499, \
                         lg_506, lg_508 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_102[k] = -f_298 * lg_34[k]
                   - f_298 * lg_41[k]
                   + f_299 * lg_43[k]
                   + f_298 * lg_109[k]
                   + f_298 * lg_116[k]
                   - f_299 * lg_118[k]
                   + f_295 * lg_139[k]
                   + f_295 * lg_146[k]
                   - f_304 * lg_148[k]
                   + f_294 * lg_244[k]
                   + f_294 * lg_251[k]
                   - f_295 * lg_253[k]
                   - f_300 * lg_274[k]
                   - f_300 * lg_281[k]
                   + f_301 * lg_283[k]
                   - f_305 * lg_304[k]
                   - f_305 * lg_311[k]
                   + f_306 * lg_313[k]
                   + f_292 * lg_439[k]
                   + f_292 * lg_446[k]
                   - f_293 * lg_448[k]
                   - f_296 * lg_469[k]
                   - f_296 * lg_476[k]
                   + f_297 * lg_478[k]
                   + f_302 * lg_499[k]
                   + f_302 * lg_506[k]
                   - f_303 * lg_508[k];
    }

#pragma omp simd aligned(lg_30, lg_33, lg_35, lg_40, lg_42, lg_44, lg_105, lg_108, lg_110, \
                         lg_115, lg_117, lg_119, lg_135, lg_138, lg_140, lg_145, lg_147, \
                         lg_149, lg_240, lg_243, lg_245, lg_250, lg_252, lg_254, lg_270, \
                         lg_273, lg_275, lg_280, lg_282, lg_284, lg_300, lg_303, lg_305, \
                         lg_310, lg_312, lg_314, lg_435, lg_438, lg_440, lg_445, lg_447, \
                         lg_449, lg_465, lg_468, lg_470, lg_475, lg_477, lg_479, lg_495, \
                         lg_498, lg_500, lg_505, lg_507, lg_509 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_103[k] = f_318 * lg_30[k]
                   + f_319 * lg_33[k]
                   - f_310 * lg_35[k]
                   + f_318 * lg_40[k]
                   - f_310 * lg_42[k]
                   + f_320 * lg_44[k]
                   - f_318 * lg_105[k]
                   - f_319 * lg_108[k]
                   + f_310 * lg_110[k]
                   - f_318 * lg_115[k]
                   + f_310 * lg_117[k]
                   - f_320 * lg_119[k]
                   - f_328 * lg_135[k]
                   - f_314 * lg_138[k]
                   + f_317 * lg_140[k]
                   - f_328 * lg_145[k]
                   + f_317 * lg_147[k]
                   - f_329 * lg_149[k]
                   - f_311 * lg_240[k]
                   - f_312 * lg_243[k]
                   + f_313 * lg_245[k]
                   - f_311 * lg_250[k]
                   + f_313 * lg_252[k]
                   - f_314 * lg_254[k]
                   + f_314 * lg_270[k]
                   + f_321 * lg_273[k]
                   - f_322 * lg_275[k]
                   + f_314 * lg_280[k]
                   - f_322 * lg_282[k]
                   + f_323 * lg_284[k]
                   + f_330 * lg_300[k]
                   + f_331 * lg_303[k]
                   - f_327 * lg_305[k]
                   + f_330 * lg_310[k]
                   - f_327 * lg_312[k]
                   + f_332 * lg_314[k]
                   - f_307 * lg_435[k]
                   - f_308 * lg_438[k]
                   + f_309 * lg_440[k]
                   - f_307 * lg_445[k]
                   + f_309 * lg_447[k]
                   - f_310 * lg_449[k]
                   + f_315 * lg_465[k]
                   + f_313 * lg_468[k]
                   - f_316 * lg_470[k]
                   + f_315 * lg_475[k]
                   - f_316 * lg_477[k]
                   + f_317 * lg_479[k]
                   - f_324 * lg_495[k]
                   - f_325 * lg_498[k]
                   + f_326 * lg_500[k]
                   - f_324 * lg_505[k]
                   + f_326 * lg_507[k]
                   - f_327 * lg_509[k];
    }

#pragma omp simd aligned(lg_32, lg_37, lg_39, lg_107, lg_112, lg_114, lg_137, lg_142, lg_144, \
                         lg_242, lg_247, lg_249, lg_272, lg_277, lg_279, lg_302, lg_307, \
                         lg_309, lg_437, lg_442, lg_444, lg_467, lg_472, lg_474, lg_497, \
                         lg_502, lg_504 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_104[k] = -f_298 * lg_32[k]
                   - f_298 * lg_37[k]
                   + f_299 * lg_39[k]
                   + f_298 * lg_107[k]
                   + f_298 * lg_112[k]
                   - f_299 * lg_114[k]
                   + f_295 * lg_137[k]
                   + f_295 * lg_142[k]
                   - f_304 * lg_144[k]
                   + f_294 * lg_242[k]
                   + f_294 * lg_247[k]
                   - f_295 * lg_249[k]
                   - f_300 * lg_272[k]
                   - f_300 * lg_277[k]
                   + f_301 * lg_279[k]
                   - f_305 * lg_302[k]
                   - f_305 * lg_307[k]
                   + f_306 * lg_309[k]
                   + f_292 * lg_437[k]
                   + f_292 * lg_442[k]
                   - f_293 * lg_444[k]
                   - f_296 * lg_467[k]
                   - f_296 * lg_472[k]
                   + f_297 * lg_474[k]
                   + f_302 * lg_497[k]
                   + f_302 * lg_502[k]
                   - f_303 * lg_504[k];
    }

#pragma omp simd aligned(lg_30, lg_35, lg_40, lg_42, lg_105, lg_110, lg_115, lg_117, lg_135, \
                         lg_140, lg_145, lg_147, lg_240, lg_245, lg_250, lg_252, lg_270, \
                         lg_275, lg_280, lg_282, lg_300, lg_305, lg_310, lg_312, lg_435, \
                         lg_440, lg_445, lg_447, lg_465, lg_470, lg_475, lg_477, lg_495, \
                         lg_500, lg_505, lg_507 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_105[k] = -f_339 * lg_30[k]
                   + f_276 * lg_35[k]
                   + f_339 * lg_40[k]
                   - f_276 * lg_42[k]
                   + f_339 * lg_105[k]
                   - f_276 * lg_110[k]
                   - f_339 * lg_115[k]
                   + f_276 * lg_117[k]
                   + f_342 * lg_135[k]
                   - f_280 * lg_140[k]
                   - f_342 * lg_145[k]
                   + f_280 * lg_147[k]
                   + f_335 * lg_240[k]
                   - f_336 * lg_245[k]
                   - f_335 * lg_250[k]
                   + f_336 * lg_252[k]
                   - f_288 * lg_270[k]
                   + f_289 * lg_275[k]
                   + f_288 * lg_280[k]
                   - f_289 * lg_282[k]
                   - f_343 * lg_300[k]
                   + f_286 * lg_305[k]
                   + f_343 * lg_310[k]
                   - f_286 * lg_312[k]
                   + f_333 * lg_435[k]
                   - f_334 * lg_440[k]
                   - f_333 * lg_445[k]
                   + f_334 * lg_447[k]
                   - f_337 * lg_465[k]
                   + f_338 * lg_470[k]
                   + f_337 * lg_475[k]
                   - f_338 * lg_477[k]
                   + f_340 * lg_495[k]
                   - f_341 * lg_500[k]
                   - f_340 * lg_505[k]
                   + f_341 * lg_507[k];
    }

#pragma omp simd aligned(lg_32, lg_37, lg_107, lg_112, lg_137, lg_142, lg_242, lg_247, lg_272, \
                         lg_277, lg_302, lg_307, lg_437, lg_442, lg_467, lg_472, lg_497, \
                         lg_502 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_106[k] = f_269 * lg_32[k]
                   - f_264 * lg_37[k]
                   - f_269 * lg_107[k]
                   + f_264 * lg_112[k]
                   - f_274 * lg_137[k]
                   + f_268 * lg_142[k]
                   - f_266 * lg_242[k]
                   + f_265 * lg_247[k]
                   + f_271 * lg_272[k]
                   - f_270 * lg_277[k]
                   + f_275 * lg_302[k]
                   - f_273 * lg_307[k]
                   - f_264 * lg_437[k]
                   + f_263 * lg_442[k]
                   + f_268 * lg_467[k]
                   - f_267 * lg_472[k]
                   - f_273 * lg_497[k]
                   + f_272 * lg_502[k];
    }

#pragma omp simd aligned(lg_30, lg_33, lg_40, lg_105, lg_108, lg_115, lg_135, lg_138, lg_145, \
                         lg_240, lg_243, lg_250, lg_270, lg_273, lg_280, lg_300, lg_303, \
                         lg_310, lg_435, lg_438, lg_445, lg_465, lg_468, lg_475, lg_495, \
                         lg_498, lg_505 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_107[k] = f_349 * lg_30[k]
                   - f_350 * lg_33[k]
                   + f_349 * lg_40[k]
                   - f_349 * lg_105[k]
                   + f_350 * lg_108[k]
                   - f_349 * lg_115[k]
                   - f_354 * lg_135[k]
                   + f_355 * lg_138[k]
                   - f_354 * lg_145[k]
                   - f_346 * lg_240[k]
                   + f_347 * lg_243[k]
                   - f_346 * lg_250[k]
                   + f_351 * lg_270[k]
                   - f_257 * lg_273[k]
                   + f_351 * lg_280[k]
                   + f_356 * lg_300[k]
                   - f_357 * lg_303[k]
                   + f_356 * lg_310[k]
                   - f_344 * lg_435[k]
                   + f_345 * lg_438[k]
                   - f_344 * lg_445[k]
                   + f_256 * lg_465[k]
                   - f_348 * lg_468[k]
                   + f_256 * lg_475[k]
                   - f_352 * lg_495[k]
                   + f_353 * lg_498[k]
                   - f_352 * lg_505[k];
    }

#pragma omp simd aligned(lg_1, lg_6, lg_46, lg_51, lg_76, lg_81, lg_151, lg_156, lg_181, \
                         lg_186, lg_211, lg_216, lg_316, lg_321, lg_346, lg_351, lg_376, \
                         lg_381, lg_541, lg_546, lg_571, lg_576, lg_601, \
                         lg_606 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_108[k] = f_249 * lg_1[k]
                   - f_249 * lg_6[k]
                   - f_210 * lg_46[k]
                   + f_210 * lg_51[k]
                   - f_251 * lg_76[k]
                   + f_251 * lg_81[k]
                   - f_569 * lg_151[k]
                   + f_569 * lg_156[k]
                   + f_570 * lg_181[k]
                   - f_570 * lg_186[k]
                   + f_253 * lg_211[k]
                   - f_253 * lg_216[k]
                   - f_210 * lg_316[k]
                   + f_210 * lg_321[k]
                   + f_570 * lg_346[k]
                   - f_570 * lg_351[k]
                   - f_254 * lg_376[k]
                   + f_254 * lg_381[k]
                   + f_249 * lg_541[k]
                   - f_249 * lg_546[k]
                   - f_251 * lg_571[k]
                   + f_251 * lg_576[k]
                   + f_253 * lg_601[k]
                   - f_253 * lg_606[k];
    }

#pragma omp simd aligned(lg_4, lg_11, lg_49, lg_56, lg_79, lg_86, lg_154, lg_161, lg_184, \
                         lg_191, lg_214, lg_221, lg_319, lg_326, lg_349, lg_356, lg_379, \
                         lg_386, lg_544, lg_551, lg_574, lg_581, lg_604, \
                         lg_611 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_109[k] = f_571 * lg_4[k]
                   - f_572 * lg_11[k]
                   - f_213 * lg_49[k]
                   + f_214 * lg_56[k]
                   - f_573 * lg_79[k]
                   + f_574 * lg_86[k]
                   - f_575 * lg_154[k]
                   + f_576 * lg_161[k]
                   + f_577 * lg_184[k]
                   - f_578 * lg_191[k]
                   + f_578 * lg_214[k]
                   - f_579 * lg_221[k]
                   - f_213 * lg_319[k]
                   + f_214 * lg_326[k]
                   + f_577 * lg_349[k]
                   - f_578 * lg_356[k]
                   - f_580 * lg_379[k]
                   + f_581 * lg_386[k]
                   + f_571 * lg_544[k]
                   - f_572 * lg_551[k]
                   - f_573 * lg_574[k]
                   + f_574 * lg_581[k]
                   + f_578 * lg_604[k]
                   - f_579 * lg_611[k];
    }

#pragma omp simd aligned(lg_1, lg_6, lg_8, lg_46, lg_51, lg_53, lg_76, lg_81, lg_83, lg_151, \
                         lg_156, lg_158, lg_181, lg_186, lg_188, lg_211, lg_216, lg_218, \
                         lg_316, lg_321, lg_323, lg_346, lg_351, lg_353, lg_376, lg_381, \
                         lg_383, lg_541, lg_546, lg_548, lg_571, lg_576, lg_578, lg_601, \
                         lg_606, lg_608 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_110[k] = -f_582 * lg_1[k]
                   - f_582 * lg_6[k]
                   + f_583 * lg_8[k]
                   + f_219 * lg_46[k]
                   + f_219 * lg_51[k]
                   - f_220 * lg_53[k]
                   + f_220 * lg_76[k]
                   + f_220 * lg_81[k]
                   - f_584 * lg_83[k]
                   + f_585 * lg_151[k]
                   + f_585 * lg_156[k]
                   - f_586 * lg_158[k]
                   - f_587 * lg_181[k]
                   - f_587 * lg_186[k]
                   + f_588 * lg_188[k]
                   - f_589 * lg_211[k]
                   - f_589 * lg_216[k]
                   + f_590 * lg_218[k]
                   + f_219 * lg_316[k]
                   + f_219 * lg_321[k]
                   - f_220 * lg_323[k]
                   - f_587 * lg_346[k]
                   - f_587 * lg_351[k]
                   + f_588 * lg_353[k]
                   + f_590 * lg_376[k]
                   + f_590 * lg_381[k]
                   - f_591 * lg_383[k]
                   - f_582 * lg_541[k]
                   - f_582 * lg_546[k]
                   + f_583 * lg_548[k]
                   + f_220 * lg_571[k]
                   + f_220 * lg_576[k]
                   - f_584 * lg_578[k]
                   - f_589 * lg_601[k]
                   - f_589 * lg_606[k]
                   + f_590 * lg_608[k];
    }

#pragma omp simd aligned(lg_4, lg_11, lg_13, lg_49, lg_56, lg_58, lg_79, lg_86, lg_88, lg_154, \
                         lg_161, lg_163, lg_184, lg_191, lg_193, lg_214, lg_221, lg_223, \
                         lg_319, lg_326, lg_328, lg_349, lg_356, lg_358, lg_379, lg_386, \
                         lg_388, lg_544, lg_551, lg_553, lg_574, lg_581, lg_583, lg_604, \
                         lg_611, lg_613 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_111[k] = -f_592 * lg_4[k]
                   - f_592 * lg_11[k]
                   + f_593 * lg_13[k]
                   + f_225 * lg_49[k]
                   + f_225 * lg_56[k]
                   - f_226 * lg_58[k]
                   + f_594 * lg_79[k]
                   + f_594 * lg_86[k]
                   - f_595 * lg_88[k]
                   + f_596 * lg_154[k]
                   + f_596 * lg_161[k]
                   - f_597 * lg_163[k]
                   - f_598 * lg_184[k]
                   - f_598 * lg_191[k]
                   + f_229 * lg_193[k]
                   - f_599 * lg_214[k]
                   - f_599 * lg_221[k]
                   + f_600 * lg_223[k]
                   + f_225 * lg_319[k]
                   + f_225 * lg_326[k]
                   - f_226 * lg_328[k]
                   - f_598 * lg_349[k]
                   - f_598 * lg_356[k]
                   + f_229 * lg_358[k]
                   + f_601 * lg_379[k]
                   + f_601 * lg_386[k]
                   - f_602 * lg_388[k]
                   - f_592 * lg_544[k]
                   - f_592 * lg_551[k]
                   + f_593 * lg_553[k]
                   + f_594 * lg_574[k]
                   + f_594 * lg_581[k]
                   - f_595 * lg_583[k]
                   - f_599 * lg_604[k]
                   - f_599 * lg_611[k]
                   + f_600 * lg_613[k];
    }

#pragma omp simd aligned(lg_0, lg_3, lg_5, lg_10, lg_12, lg_14, lg_45, lg_48, lg_50, lg_55, \
                         lg_57, lg_59, lg_75, lg_78, lg_80, lg_85, lg_87, lg_89, lg_150, \
                         lg_153, lg_155, lg_160, lg_162, lg_164, lg_180, lg_183, lg_185, \
                         lg_190, lg_192, lg_194, lg_210, lg_213, lg_215, lg_220, lg_222, \
                         lg_224, lg_315, lg_318, lg_320, lg_325, lg_327, lg_329, lg_345, \
                         lg_348, lg_350, lg_355, lg_357, lg_359, lg_375, lg_378, lg_380, \
                         lg_385, lg_387, lg_389, lg_540, lg_543, lg_545, lg_550, lg_552, \
                         lg_554, lg_570, lg_573, lg_575, lg_580, lg_582, lg_584, lg_600, \
                         lg_603, lg_605, lg_610, lg_612, lg_614 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_112[k] = f_603 * lg_0[k]
                   + f_604 * lg_3[k]
                   - f_232 * lg_5[k]
                   + f_603 * lg_10[k]
                   - f_232 * lg_12[k]
                   + f_605 * lg_14[k]
                   - f_231 * lg_45[k]
                   - f_232 * lg_48[k]
                   + f_233 * lg_50[k]
                   - f_231 * lg_55[k]
                   + f_233 * lg_57[k]
                   - f_234 * lg_59[k]
                   - f_606 * lg_75[k]
                   - f_607 * lg_78[k]
                   + f_236 * lg_80[k]
                   - f_606 * lg_85[k]
                   + f_236 * lg_87[k]
                   - f_608 * lg_89[k]
                   - f_609 * lg_150[k]
                   - f_610 * lg_153[k]
                   + f_611 * lg_155[k]
                   - f_609 * lg_160[k]
                   + f_611 * lg_162[k]
                   - f_612 * lg_164[k]
                   + f_613 * lg_180[k]
                   + f_614 * lg_183[k]
                   - f_615 * lg_185[k]
                   + f_613 * lg_190[k]
                   - f_615 * lg_192[k]
                   + f_240 * lg_194[k]
                   + f_616 * lg_210[k]
                   + f_611 * lg_213[k]
                   - f_240 * lg_215[k]
                   + f_616 * lg_220[k]
                   - f_240 * lg_222[k]
                   + f_617 * lg_224[k]
                   - f_231 * lg_315[k]
                   - f_232 * lg_318[k]
                   + f_233 * lg_320[k]
                   - f_231 * lg_325[k]
                   + f_233 * lg_327[k]
                   - f_234 * lg_329[k]
                   + f_613 * lg_345[k]
                   + f_614 * lg_348[k]
                   - f_615 * lg_350[k]
                   + f_613 * lg_355[k]
                   - f_615 * lg_357[k]
                   + f_240 * lg_359[k]
                   - f_614 * lg_375[k]
                   - f_618 * lg_378[k]
                   + f_619 * lg_380[k]
                   - f_614 * lg_385[k]
                   + f_619 * lg_387[k]
                   - f_620 * lg_389[k]
                   + f_603 * lg_540[k]
                   + f_604 * lg_543[k]
                   - f_232 * lg_545[k]
                   + f_603 * lg_550[k]
                   - f_232 * lg_552[k]
                   + f_605 * lg_554[k]
                   - f_606 * lg_570[k]
                   - f_607 * lg_573[k]
                   + f_236 * lg_575[k]
                   - f_606 * lg_580[k]
                   + f_236 * lg_582[k]
                   - f_608 * lg_584[k]
                   + f_616 * lg_600[k]
                   + f_611 * lg_603[k]
                   - f_240 * lg_605[k]
                   + f_616 * lg_610[k]
                   - f_240 * lg_612[k]
                   + f_617 * lg_614[k];
    }

#pragma omp simd aligned(lg_2, lg_7, lg_9, lg_47, lg_52, lg_54, lg_77, lg_82, lg_84, lg_152, \
                         lg_157, lg_159, lg_182, lg_187, lg_189, lg_212, lg_217, lg_219, \
                         lg_317, lg_322, lg_324, lg_347, lg_352, lg_354, lg_377, lg_382, \
                         lg_384, lg_542, lg_547, lg_549, lg_572, lg_577, lg_579, lg_602, \
                         lg_607, lg_609 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_113[k] = -f_592 * lg_2[k]
                   - f_592 * lg_7[k]
                   + f_593 * lg_9[k]
                   + f_225 * lg_47[k]
                   + f_225 * lg_52[k]
                   - f_226 * lg_54[k]
                   + f_594 * lg_77[k]
                   + f_594 * lg_82[k]
                   - f_595 * lg_84[k]
                   + f_596 * lg_152[k]
                   + f_596 * lg_157[k]
                   - f_597 * lg_159[k]
                   - f_598 * lg_182[k]
                   - f_598 * lg_187[k]
                   + f_229 * lg_189[k]
                   - f_599 * lg_212[k]
                   - f_599 * lg_217[k]
                   + f_600 * lg_219[k]
                   + f_225 * lg_317[k]
                   + f_225 * lg_322[k]
                   - f_226 * lg_324[k]
                   - f_598 * lg_347[k]
                   - f_598 * lg_352[k]
                   + f_229 * lg_354[k]
                   + f_601 * lg_377[k]
                   + f_601 * lg_382[k]
                   - f_602 * lg_384[k]
                   - f_592 * lg_542[k]
                   - f_592 * lg_547[k]
                   + f_593 * lg_549[k]
                   + f_594 * lg_572[k]
                   + f_594 * lg_577[k]
                   - f_595 * lg_579[k]
                   - f_599 * lg_602[k]
                   - f_599 * lg_607[k]
                   + f_600 * lg_609[k];
    }

#pragma omp simd aligned(lg_0, lg_5, lg_10, lg_12, lg_45, lg_50, lg_55, lg_57, lg_75, lg_80, \
                         lg_85, lg_87, lg_150, lg_155, lg_160, lg_162, lg_180, lg_185, lg_190, \
                         lg_192, lg_210, lg_215, lg_220, lg_222, lg_315, lg_320, lg_325, \
                         lg_327, lg_345, lg_350, lg_355, lg_357, lg_375, lg_380, lg_385, \
                         lg_387, lg_540, lg_545, lg_550, lg_552, lg_570, lg_575, lg_580, \
                         lg_582, lg_600, lg_605, lg_610, lg_612 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_114[k] = -f_621 * lg_0[k]
                   + f_622 * lg_5[k]
                   + f_621 * lg_10[k]
                   - f_622 * lg_12[k]
                   + f_243 * lg_45[k]
                   - f_244 * lg_50[k]
                   - f_243 * lg_55[k]
                   + f_244 * lg_57[k]
                   + f_244 * lg_75[k]
                   - f_623 * lg_80[k]
                   - f_244 * lg_85[k]
                   + f_623 * lg_87[k]
                   + f_624 * lg_150[k]
                   - f_625 * lg_155[k]
                   - f_624 * lg_160[k]
                   + f_625 * lg_162[k]
                   - f_586 * lg_180[k]
                   + f_626 * lg_185[k]
                   + f_586 * lg_190[k]
                   - f_626 * lg_192[k]
                   - f_627 * lg_210[k]
                   + f_587 * lg_215[k]
                   + f_627 * lg_220[k]
                   - f_587 * lg_222[k]
                   + f_243 * lg_315[k]
                   - f_244 * lg_320[k]
                   - f_243 * lg_325[k]
                   + f_244 * lg_327[k]
                   - f_586 * lg_345[k]
                   + f_626 * lg_350[k]
                   + f_586 * lg_355[k]
                   - f_626 * lg_357[k]
                   + f_587 * lg_375[k]
                   - f_588 * lg_380[k]
                   - f_587 * lg_385[k]
                   + f_588 * lg_387[k]
                   - f_621 * lg_540[k]
                   + f_622 * lg_545[k]
                   + f_621 * lg_550[k]
                   - f_622 * lg_552[k]
                   + f_244 * lg_570[k]
                   - f_623 * lg_575[k]
                   - f_244 * lg_580[k]
                   + f_623 * lg_582[k]
                   - f_627 * lg_600[k]
                   + f_587 * lg_605[k]
                   + f_627 * lg_610[k]
                   - f_587 * lg_612[k];
    }

#pragma omp simd aligned(lg_2, lg_7, lg_47, lg_52, lg_77, lg_82, lg_152, lg_157, lg_182, \
                         lg_187, lg_212, lg_217, lg_317, lg_322, lg_347, lg_352, lg_377, \
                         lg_382, lg_542, lg_547, lg_572, lg_577, lg_602, \
                         lg_607 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_115[k] = f_572 * lg_2[k]
                   - f_571 * lg_7[k]
                   - f_214 * lg_47[k]
                   + f_213 * lg_52[k]
                   - f_574 * lg_77[k]
                   + f_573 * lg_82[k]
                   - f_576 * lg_152[k]
                   + f_575 * lg_157[k]
                   + f_578 * lg_182[k]
                   - f_577 * lg_187[k]
                   + f_579 * lg_212[k]
                   - f_578 * lg_217[k]
                   - f_214 * lg_317[k]
                   + f_213 * lg_322[k]
                   + f_578 * lg_347[k]
                   - f_577 * lg_352[k]
                   - f_581 * lg_377[k]
                   + f_580 * lg_382[k]
                   + f_572 * lg_542[k]
                   - f_571 * lg_547[k]
                   - f_574 * lg_572[k]
                   + f_573 * lg_577[k]
                   + f_579 * lg_602[k]
                   - f_578 * lg_607[k];
    }

#pragma omp simd aligned(lg_0, lg_3, lg_10, lg_45, lg_48, lg_55, lg_75, lg_78, lg_85, lg_150, \
                         lg_153, lg_160, lg_180, lg_183, lg_190, lg_210, lg_213, lg_220, \
                         lg_315, lg_318, lg_325, lg_345, lg_348, lg_355, lg_375, lg_378, \
                         lg_385, lg_540, lg_543, lg_550, lg_570, lg_573, lg_580, lg_600, \
                         lg_603, lg_610 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_116[k] = f_628 * lg_0[k]
                   - f_629 * lg_3[k]
                   + f_628 * lg_10[k]
                   - f_249 * lg_45[k]
                   + f_250 * lg_48[k]
                   - f_249 * lg_55[k]
                   - f_250 * lg_75[k]
                   + f_630 * lg_78[k]
                   - f_250 * lg_85[k]
                   - f_631 * lg_150[k]
                   + f_632 * lg_153[k]
                   - f_631 * lg_160[k]
                   + f_633 * lg_180[k]
                   - f_634 * lg_183[k]
                   + f_633 * lg_190[k]
                   + f_569 * lg_210[k]
                   - f_635 * lg_213[k]
                   + f_569 * lg_220[k]
                   - f_249 * lg_315[k]
                   + f_250 * lg_318[k]
                   - f_249 * lg_325[k]
                   + f_633 * lg_345[k]
                   - f_634 * lg_348[k]
                   + f_633 * lg_355[k]
                   - f_635 * lg_375[k]
                   + f_636 * lg_378[k]
                   - f_635 * lg_385[k]
                   + f_628 * lg_540[k]
                   - f_629 * lg_543[k]
                   + f_628 * lg_550[k]
                   - f_250 * lg_570[k]
                   + f_630 * lg_573[k]
                   - f_250 * lg_580[k]
                   + f_569 * lg_600[k]
                   - f_635 * lg_603[k]
                   + f_569 * lg_610[k];
    }

#pragma omp simd aligned(lg_31, lg_36, lg_106, lg_111, lg_136, lg_141, lg_241, lg_246, lg_271, \
                         lg_276, lg_436, lg_441, lg_466, lg_471 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_117[k] = -f_141 * lg_31[k]
                   + f_141 * lg_36[k]
                   + f_139 * lg_106[k]
                   - f_139 * lg_111[k]
                   + f_21 * lg_136[k]
                   - f_21 * lg_141[k]
                   + f_137 * lg_241[k]
                   - f_137 * lg_246[k]
                   - f_140 * lg_271[k]
                   + f_140 * lg_276[k]
                   - f_137 * lg_436[k]
                   + f_137 * lg_441[k]
                   + f_138 * lg_466[k]
                   - f_138 * lg_471[k];
    }

#pragma omp simd aligned(lg_34, lg_41, lg_109, lg_116, lg_139, lg_146, lg_244, lg_251, lg_274, \
                         lg_281, lg_439, lg_446, lg_469, lg_476 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_118[k] = -f_150 * lg_34[k]
                   + f_151 * lg_41[k]
                   + f_146 * lg_109[k]
                   - f_147 * lg_116[k]
                   + f_152 * lg_139[k]
                   - f_153 * lg_146[k]
                   + f_142 * lg_244[k]
                   - f_143 * lg_251[k]
                   - f_148 * lg_274[k]
                   + f_149 * lg_281[k]
                   - f_142 * lg_439[k]
                   + f_143 * lg_446[k]
                   + f_144 * lg_469[k]
                   - f_145 * lg_476[k];
    }

#pragma omp simd aligned(lg_31, lg_36, lg_38, lg_106, lg_111, lg_113, lg_136, lg_141, lg_143, \
                         lg_241, lg_246, lg_248, lg_271, lg_276, lg_278, lg_436, lg_441, \
                         lg_443, lg_466, lg_471, lg_473 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_119[k] = f_162 * lg_31[k]
                   + f_162 * lg_36[k]
                   - f_163 * lg_38[k]
                   - f_158 * lg_106[k]
                   - f_158 * lg_111[k]
                   + f_159 * lg_113[k]
                   - f_164 * lg_136[k]
                   - f_164 * lg_141[k]
                   + f_165 * lg_143[k]
                   - f_154 * lg_241[k]
                   - f_154 * lg_246[k]
                   + f_155 * lg_248[k]
                   + f_160 * lg_271[k]
                   + f_160 * lg_276[k]
                   - f_161 * lg_278[k]
                   + f_154 * lg_436[k]
                   + f_154 * lg_441[k]
                   - f_155 * lg_443[k]
                   - f_156 * lg_466[k]
                   - f_156 * lg_471[k]
                   + f_157 * lg_473[k];
    }

#pragma omp simd aligned(lg_34, lg_41, lg_43, lg_109, lg_116, lg_118, lg_139, lg_146, lg_148, \
                         lg_244, lg_251, lg_253, lg_274, lg_281, lg_283, lg_439, lg_446, \
                         lg_448, lg_469, lg_476, lg_478 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_120[k] = f_174 * lg_34[k]
                   + f_174 * lg_41[k]
                   - f_175 * lg_43[k]
                   - f_170 * lg_109[k]
                   - f_170 * lg_116[k]
                   + f_171 * lg_118[k]
                   - f_176 * lg_139[k]
                   - f_176 * lg_146[k]
                   + f_177 * lg_148[k]
                   - f_166 * lg_244[k]
                   - f_166 * lg_251[k]
                   + f_167 * lg_253[k]
                   + f_172 * lg_274[k]
                   + f_172 * lg_281[k]
                   - f_173 * lg_283[k]
                   + f_166 * lg_439[k]
                   + f_166 * lg_446[k]
                   - f_167 * lg_448[k]
                   - f_168 * lg_469[k]
                   - f_168 * lg_476[k]
                   + f_169 * lg_478[k];
    }

#pragma omp simd aligned(lg_30, lg_33, lg_35, lg_40, lg_42, lg_44, lg_105, lg_108, lg_110, \
                         lg_115, lg_117, lg_119, lg_135, lg_138, lg_140, lg_145, lg_147, \
                         lg_149, lg_240, lg_243, lg_245, lg_250, lg_252, lg_254, lg_270, \
                         lg_273, lg_275, lg_280, lg_282, lg_284, lg_435, lg_438, lg_440, \
                         lg_445, lg_447, lg_449, lg_465, lg_468, lg_470, lg_475, lg_477, \
                         lg_479 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_121[k] = -f_189 * lg_30[k]
                   - f_190 * lg_33[k]
                   + f_191 * lg_35[k]
                   - f_189 * lg_40[k]
                   + f_191 * lg_42[k]
                   - f_192 * lg_44[k]
                   + f_182 * lg_105[k]
                   + f_183 * lg_108[k]
                   - f_184 * lg_110[k]
                   + f_182 * lg_115[k]
                   - f_184 * lg_117[k]
                   + f_185 * lg_119[k]
                   + f_193 * lg_135[k]
                   + f_191 * lg_138[k]
                   - f_194 * lg_140[k]
                   + f_193 * lg_145[k]
                   - f_194 * lg_147[k]
                   + f_195 * lg_149[k]
                   + f_178 * lg_240[k]
                   + f_79 * lg_243[k]
                   - f_179 * lg_245[k]
                   + f_178 * lg_250[k]
                   - f_179 * lg_252[k]
                   + f_0 * lg_254[k]
                   - f_179 * lg_270[k]
                   - f_186 * lg_273[k]
                   + f_187 * lg_275[k]
                   - f_179 * lg_280[k]
                   + f_187 * lg_282[k]
                   - f_188 * lg_284[k]
                   - f_178 * lg_435[k]
                   - f_79 * lg_438[k]
                   + f_179 * lg_440[k]
                   - f_178 * lg_445[k]
                   + f_179 * lg_447[k]
                   - f_0 * lg_449[k]
                   + f_27 * lg_465[k]
                   + f_179 * lg_468[k]
                   - f_180 * lg_470[k]
                   + f_27 * lg_475[k]
                   - f_180 * lg_477[k]
                   + f_181 * lg_479[k];
    }

#pragma omp simd aligned(lg_32, lg_37, lg_39, lg_107, lg_112, lg_114, lg_137, lg_142, lg_144, \
                         lg_242, lg_247, lg_249, lg_272, lg_277, lg_279, lg_437, lg_442, \
                         lg_444, lg_467, lg_472, lg_474 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_122[k] = f_174 * lg_32[k]
                   + f_174 * lg_37[k]
                   - f_175 * lg_39[k]
                   - f_170 * lg_107[k]
                   - f_170 * lg_112[k]
                   + f_171 * lg_114[k]
                   - f_176 * lg_137[k]
                   - f_176 * lg_142[k]
                   + f_177 * lg_144[k]
                   - f_166 * lg_242[k]
                   - f_166 * lg_247[k]
                   + f_167 * lg_249[k]
                   + f_172 * lg_272[k]
                   + f_172 * lg_277[k]
                   - f_173 * lg_279[k]
                   + f_166 * lg_437[k]
                   + f_166 * lg_442[k]
                   - f_167 * lg_444[k]
                   - f_168 * lg_467[k]
                   - f_168 * lg_472[k]
                   + f_169 * lg_474[k];
    }

#pragma omp simd aligned(lg_30, lg_35, lg_40, lg_42, lg_105, lg_110, lg_115, lg_117, lg_135, \
                         lg_140, lg_145, lg_147, lg_240, lg_245, lg_250, lg_252, lg_270, \
                         lg_275, lg_280, lg_282, lg_435, lg_440, lg_445, lg_447, lg_465, \
                         lg_470, lg_475, lg_477 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_123[k] = f_202 * lg_30[k]
                   - f_203 * lg_35[k]
                   - f_202 * lg_40[k]
                   + f_203 * lg_42[k]
                   - f_200 * lg_105[k]
                   + f_201 * lg_110[k]
                   + f_200 * lg_115[k]
                   - f_201 * lg_117[k]
                   - f_204 * lg_135[k]
                   + f_205 * lg_140[k]
                   + f_204 * lg_145[k]
                   - f_205 * lg_147[k]
                   - f_196 * lg_240[k]
                   + f_197 * lg_245[k]
                   + f_196 * lg_250[k]
                   - f_197 * lg_252[k]
                   + f_156 * lg_270[k]
                   - f_157 * lg_275[k]
                   - f_156 * lg_280[k]
                   + f_157 * lg_282[k]
                   + f_196 * lg_435[k]
                   - f_197 * lg_440[k]
                   - f_196 * lg_445[k]
                   + f_197 * lg_447[k]
                   - f_198 * lg_465[k]
                   + f_199 * lg_470[k]
                   + f_198 * lg_475[k]
                   - f_199 * lg_477[k];
    }

#pragma omp simd aligned(lg_32, lg_37, lg_107, lg_112, lg_137, lg_142, lg_242, lg_247, lg_272, \
                         lg_277, lg_437, lg_442, lg_467, lg_472 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_124[k] = -f_151 * lg_32[k]
                   + f_150 * lg_37[k]
                   + f_147 * lg_107[k]
                   - f_146 * lg_112[k]
                   + f_153 * lg_137[k]
                   - f_152 * lg_142[k]
                   + f_143 * lg_242[k]
                   - f_142 * lg_247[k]
                   - f_149 * lg_272[k]
                   + f_148 * lg_277[k]
                   - f_143 * lg_437[k]
                   + f_142 * lg_442[k]
                   + f_145 * lg_467[k]
                   - f_144 * lg_472[k];
    }

#pragma omp simd aligned(lg_30, lg_33, lg_40, lg_105, lg_108, lg_115, lg_135, lg_138, lg_145, \
                         lg_240, lg_243, lg_250, lg_270, lg_273, lg_280, lg_435, lg_438, \
                         lg_445, lg_465, lg_468, lg_475 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_125[k] = -f_209 * lg_30[k]
                   + f_18 * lg_33[k]
                   - f_209 * lg_40[k]
                   + f_59 * lg_105[k]
                   - f_207 * lg_108[k]
                   + f_59 * lg_115[k]
                   + f_141 * lg_135[k]
                   - f_53 * lg_138[k]
                   + f_141 * lg_145[k]
                   + f_206 * lg_240[k]
                   - f_56 * lg_243[k]
                   + f_206 * lg_250[k]
                   - f_58 * lg_270[k]
                   + f_208 * lg_273[k]
                   - f_58 * lg_280[k]
                   - f_206 * lg_435[k]
                   + f_56 * lg_438[k]
                   - f_206 * lg_445[k]
                   + f_137 * lg_465[k]
                   - f_57 * lg_468[k]
                   + f_137 * lg_475[k];
    }

#pragma omp simd aligned(lg_1, lg_6, lg_46, lg_51, lg_76, lg_81, lg_181, lg_186, lg_316, \
                         lg_321, lg_346, lg_351, lg_541, lg_546, lg_571, \
                         lg_576 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_126[k] = -f_637 * lg_1[k]
                   + f_637 * lg_6[k]
                   + f_81 * lg_46[k]
                   - f_81 * lg_51[k]
                   + f_81 * lg_76[k]
                   - f_81 * lg_81[k]
                   - f_638 * lg_181[k]
                   + f_638 * lg_186[k]
                   - f_81 * lg_316[k]
                   + f_81 * lg_321[k]
                   + f_638 * lg_346[k]
                   - f_638 * lg_351[k]
                   + f_637 * lg_541[k]
                   - f_637 * lg_546[k]
                   - f_81 * lg_571[k]
                   + f_81 * lg_576[k];
    }

#pragma omp simd aligned(lg_4, lg_11, lg_49, lg_56, lg_79, lg_86, lg_184, lg_191, lg_319, \
                         lg_326, lg_349, lg_356, lg_544, lg_551, lg_574, \
                         lg_581 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_127[k] = -f_639 * lg_4[k]
                   + f_640 * lg_11[k]
                   + f_86 * lg_49[k]
                   - f_87 * lg_56[k]
                   + f_86 * lg_79[k]
                   - f_87 * lg_86[k]
                   - f_641 * lg_184[k]
                   + f_642 * lg_191[k]
                   - f_86 * lg_319[k]
                   + f_87 * lg_326[k]
                   + f_641 * lg_349[k]
                   - f_642 * lg_356[k]
                   + f_639 * lg_544[k]
                   - f_640 * lg_551[k]
                   - f_86 * lg_574[k]
                   + f_87 * lg_581[k];
    }

#pragma omp simd aligned(lg_1, lg_6, lg_8, lg_46, lg_51, lg_53, lg_76, lg_81, lg_83, lg_181, \
                         lg_186, lg_188, lg_316, lg_321, lg_323, lg_346, lg_351, lg_353, \
                         lg_541, lg_546, lg_548, lg_571, lg_576, \
                         lg_578 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_128[k] = f_643 * lg_1[k]
                   + f_643 * lg_6[k]
                   - f_92 * lg_8[k]
                   - f_94 * lg_46[k]
                   - f_94 * lg_51[k]
                   + f_95 * lg_53[k]
                   - f_94 * lg_76[k]
                   - f_94 * lg_81[k]
                   + f_95 * lg_83[k]
                   + f_644 * lg_181[k]
                   + f_644 * lg_186[k]
                   - f_645 * lg_188[k]
                   + f_94 * lg_316[k]
                   + f_94 * lg_321[k]
                   - f_95 * lg_323[k]
                   - f_644 * lg_346[k]
                   - f_644 * lg_351[k]
                   + f_645 * lg_353[k]
                   - f_643 * lg_541[k]
                   - f_643 * lg_546[k]
                   + f_92 * lg_548[k]
                   + f_94 * lg_571[k]
                   + f_94 * lg_576[k]
                   - f_95 * lg_578[k];
    }

#pragma omp simd aligned(lg_4, lg_11, lg_13, lg_49, lg_56, lg_58, lg_79, lg_86, lg_88, lg_184, \
                         lg_191, lg_193, lg_319, lg_326, lg_328, lg_349, lg_356, lg_358, \
                         lg_544, lg_551, lg_553, lg_574, lg_581, \
                         lg_583 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_129[k] = f_646 * lg_4[k]
                   + f_646 * lg_11[k]
                   - f_647 * lg_13[k]
                   - f_101 * lg_49[k]
                   - f_101 * lg_56[k]
                   + f_102 * lg_58[k]
                   - f_101 * lg_79[k]
                   - f_101 * lg_86[k]
                   + f_102 * lg_88[k]
                   + f_648 * lg_184[k]
                   + f_648 * lg_191[k]
                   - f_105 * lg_193[k]
                   + f_101 * lg_319[k]
                   + f_101 * lg_326[k]
                   - f_102 * lg_328[k]
                   - f_648 * lg_349[k]
                   - f_648 * lg_356[k]
                   + f_105 * lg_358[k]
                   - f_646 * lg_544[k]
                   - f_646 * lg_551[k]
                   + f_647 * lg_553[k]
                   + f_101 * lg_574[k]
                   + f_101 * lg_581[k]
                   - f_102 * lg_583[k];
    }

#pragma omp simd aligned(lg_0, lg_3, lg_5, lg_10, lg_12, lg_14, lg_45, lg_48, lg_50, lg_55, \
                         lg_57, lg_59, lg_75, lg_78, lg_80, lg_85, lg_87, lg_89, lg_180, \
                         lg_183, lg_185, lg_190, lg_192, lg_194, lg_315, lg_318, lg_320, \
                         lg_325, lg_327, lg_329, lg_345, lg_348, lg_350, lg_355, lg_357, \
                         lg_359, lg_540, lg_543, lg_545, lg_550, lg_552, lg_554, lg_570, \
                         lg_573, lg_575, lg_580, lg_582, lg_584 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_130[k] = -f_649 * lg_0[k]
                   - f_650 * lg_3[k]
                   + f_651 * lg_5[k]
                   - f_649 * lg_10[k]
                   + f_651 * lg_12[k]
                   - f_652 * lg_14[k]
                   + f_111 * lg_45[k]
                   + f_112 * lg_48[k]
                   - f_113 * lg_50[k]
                   + f_111 * lg_55[k]
                   - f_113 * lg_57[k]
                   + f_114 * lg_59[k]
                   + f_111 * lg_75[k]
                   + f_112 * lg_78[k]
                   - f_113 * lg_80[k]
                   + f_111 * lg_85[k]
                   - f_113 * lg_87[k]
                   + f_114 * lg_89[k]
                   - f_653 * lg_180[k]
                   - f_654 * lg_183[k]
                   + f_655 * lg_185[k]
                   - f_653 * lg_190[k]
                   + f_655 * lg_192[k]
                   - f_120 * lg_194[k]
                   - f_111 * lg_315[k]
                   - f_112 * lg_318[k]
                   + f_113 * lg_320[k]
                   - f_111 * lg_325[k]
                   + f_113 * lg_327[k]
                   - f_114 * lg_329[k]
                   + f_653 * lg_345[k]
                   + f_654 * lg_348[k]
                   - f_655 * lg_350[k]
                   + f_653 * lg_355[k]
                   - f_655 * lg_357[k]
                   + f_120 * lg_359[k]
                   + f_649 * lg_540[k]
                   + f_650 * lg_543[k]
                   - f_651 * lg_545[k]
                   + f_649 * lg_550[k]
                   - f_651 * lg_552[k]
                   + f_652 * lg_554[k]
                   - f_111 * lg_570[k]
                   - f_112 * lg_573[k]
                   + f_113 * lg_575[k]
                   - f_111 * lg_580[k]
                   + f_113 * lg_582[k]
                   - f_114 * lg_584[k];
    }

#pragma omp simd aligned(lg_2, lg_7, lg_9, lg_47, lg_52, lg_54, lg_77, lg_82, lg_84, lg_182, \
                         lg_187, lg_189, lg_317, lg_322, lg_324, lg_347, lg_352, lg_354, \
                         lg_542, lg_547, lg_549, lg_572, lg_577, \
                         lg_579 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_131[k] = f_646 * lg_2[k]
                   + f_646 * lg_7[k]
                   - f_647 * lg_9[k]
                   - f_101 * lg_47[k]
                   - f_101 * lg_52[k]
                   + f_102 * lg_54[k]
                   - f_101 * lg_77[k]
                   - f_101 * lg_82[k]
                   + f_102 * lg_84[k]
                   + f_648 * lg_182[k]
                   + f_648 * lg_187[k]
                   - f_105 * lg_189[k]
                   + f_101 * lg_317[k]
                   + f_101 * lg_322[k]
                   - f_102 * lg_324[k]
                   - f_648 * lg_347[k]
                   - f_648 * lg_352[k]
                   + f_105 * lg_354[k]
                   - f_646 * lg_542[k]
                   - f_646 * lg_547[k]
                   + f_647 * lg_549[k]
                   + f_101 * lg_572[k]
                   + f_101 * lg_577[k]
                   - f_102 * lg_579[k];
    }

#pragma omp simd aligned(lg_0, lg_5, lg_10, lg_12, lg_45, lg_50, lg_55, lg_57, lg_75, lg_80, \
                         lg_85, lg_87, lg_180, lg_185, lg_190, lg_192, lg_315, lg_320, lg_325, \
                         lg_327, lg_345, lg_350, lg_355, lg_357, lg_540, lg_545, lg_550, \
                         lg_552, lg_570, lg_575, lg_580, lg_582 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_132[k] = f_656 * lg_0[k]
                   - f_123 * lg_5[k]
                   - f_656 * lg_10[k]
                   + f_123 * lg_12[k]
                   - f_125 * lg_45[k]
                   + f_126 * lg_50[k]
                   + f_125 * lg_55[k]
                   - f_126 * lg_57[k]
                   - f_125 * lg_75[k]
                   + f_126 * lg_80[k]
                   + f_125 * lg_85[k]
                   - f_126 * lg_87[k]
                   + f_657 * lg_180[k]
                   - f_658 * lg_185[k]
                   - f_657 * lg_190[k]
                   + f_658 * lg_192[k]
                   + f_125 * lg_315[k]
                   - f_126 * lg_320[k]
                   - f_125 * lg_325[k]
                   + f_126 * lg_327[k]
                   - f_657 * lg_345[k]
                   + f_658 * lg_350[k]
                   + f_657 * lg_355[k]
                   - f_658 * lg_357[k]
                   - f_656 * lg_540[k]
                   + f_123 * lg_545[k]
                   + f_656 * lg_550[k]
                   - f_123 * lg_552[k]
                   + f_125 * lg_570[k]
                   - f_126 * lg_575[k]
                   - f_125 * lg_580[k]
                   + f_126 * lg_582[k];
    }

#pragma omp simd aligned(lg_2, lg_7, lg_47, lg_52, lg_77, lg_82, lg_182, lg_187, lg_317, \
                         lg_322, lg_347, lg_352, lg_542, lg_547, lg_572, \
                         lg_577 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_133[k] = -f_640 * lg_2[k]
                   + f_639 * lg_7[k]
                   + f_87 * lg_47[k]
                   - f_86 * lg_52[k]
                   + f_87 * lg_77[k]
                   - f_86 * lg_82[k]
                   - f_642 * lg_182[k]
                   + f_641 * lg_187[k]
                   - f_87 * lg_317[k]
                   + f_86 * lg_322[k]
                   + f_642 * lg_347[k]
                   - f_641 * lg_352[k]
                   + f_640 * lg_542[k]
                   - f_639 * lg_547[k]
                   - f_87 * lg_572[k]
                   + f_86 * lg_577[k];
    }

#pragma omp simd aligned(lg_0, lg_3, lg_10, lg_45, lg_48, lg_55, lg_75, lg_78, lg_85, lg_180, \
                         lg_183, lg_190, lg_315, lg_318, lg_325, lg_345, lg_348, lg_355, \
                         lg_540, lg_543, lg_550, lg_570, lg_573, \
                         lg_580 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_134[k] = -f_659 * lg_0[k]
                   + f_130 * lg_3[k]
                   - f_659 * lg_10[k]
                   + f_132 * lg_45[k]
                   - f_133 * lg_48[k]
                   + f_132 * lg_55[k]
                   + f_132 * lg_75[k]
                   - f_133 * lg_78[k]
                   + f_132 * lg_85[k]
                   - f_660 * lg_180[k]
                   + f_661 * lg_183[k]
                   - f_660 * lg_190[k]
                   - f_132 * lg_315[k]
                   + f_133 * lg_318[k]
                   - f_132 * lg_325[k]
                   + f_660 * lg_345[k]
                   - f_661 * lg_348[k]
                   + f_660 * lg_355[k]
                   + f_659 * lg_540[k]
                   - f_130 * lg_543[k]
                   + f_659 * lg_550[k]
                   - f_132 * lg_570[k]
                   + f_133 * lg_573[k]
                   - f_132 * lg_580[k];
    }

#pragma omp simd aligned(lg_31, lg_36, lg_106, lg_111, lg_241, lg_246, lg_436, \
                         lg_441 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_135[k] = f_32 * lg_31[k]
                   - f_32 * lg_36[k]
                   - f_29 * lg_106[k]
                   + f_29 * lg_111[k]
                   + f_31 * lg_241[k]
                   - f_31 * lg_246[k]
                   - f_30 * lg_436[k]
                   + f_30 * lg_441[k];
    }

#pragma omp simd aligned(lg_34, lg_41, lg_109, lg_116, lg_244, lg_251, lg_439, \
                         lg_446 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_136[k] = f_38 * lg_34[k]
                   - f_39 * lg_41[k]
                   - f_37 * lg_109[k]
                   + f_33 * lg_116[k]
                   + f_35 * lg_244[k]
                   - f_36 * lg_251[k]
                   - f_33 * lg_439[k]
                   + f_34 * lg_446[k];
    }

#pragma omp simd aligned(lg_31, lg_36, lg_38, lg_106, lg_111, lg_113, lg_241, lg_246, lg_248, \
                         lg_436, lg_441, lg_443 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_137[k] = -f_22 * lg_31[k]
                   - f_22 * lg_36[k]
                   + f_23 * lg_38[k]
                   + f_42 * lg_106[k]
                   + f_42 * lg_111[k]
                   - f_43 * lg_113[k]
                   - f_40 * lg_241[k]
                   - f_40 * lg_246[k]
                   + f_41 * lg_248[k]
                   + f_24 * lg_436[k]
                   + f_24 * lg_441[k]
                   - f_25 * lg_443[k];
    }

#pragma omp simd aligned(lg_34, lg_41, lg_43, lg_109, lg_116, lg_118, lg_244, lg_251, lg_253, \
                         lg_439, lg_446, lg_448 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_138[k] = -f_50 * lg_34[k]
                   - f_50 * lg_41[k]
                   + f_51 * lg_43[k]
                   + f_48 * lg_109[k]
                   + f_48 * lg_116[k]
                   - f_49 * lg_118[k]
                   - f_46 * lg_244[k]
                   - f_46 * lg_251[k]
                   + f_47 * lg_253[k]
                   + f_44 * lg_439[k]
                   + f_44 * lg_446[k]
                   - f_45 * lg_448[k];
    }

#pragma omp simd aligned(lg_30, lg_33, lg_35, lg_40, lg_42, lg_44, lg_105, lg_108, lg_110, \
                         lg_115, lg_117, lg_119, lg_240, lg_243, lg_245, lg_250, lg_252, \
                         lg_254, lg_435, lg_438, lg_440, lg_445, lg_447, \
                         lg_449 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_139[k] = f_62 * lg_30[k]
                   + f_14 * lg_33[k]
                   - f_63 * lg_35[k]
                   + f_62 * lg_40[k]
                   - f_63 * lg_42[k]
                   + f_64 * lg_44[k]
                   - f_59 * lg_105[k]
                   - f_60 * lg_108[k]
                   + f_61 * lg_110[k]
                   - f_59 * lg_115[k]
                   + f_61 * lg_117[k]
                   - f_53 * lg_119[k]
                   + f_55 * lg_240[k]
                   + f_56 * lg_243[k]
                   - f_57 * lg_245[k]
                   + f_55 * lg_250[k]
                   - f_57 * lg_252[k]
                   + f_58 * lg_254[k]
                   - f_52 * lg_435[k]
                   - f_18 * lg_438[k]
                   + f_53 * lg_440[k]
                   - f_52 * lg_445[k]
                   + f_53 * lg_447[k]
                   - f_54 * lg_449[k];
    }

#pragma omp simd aligned(lg_32, lg_37, lg_39, lg_107, lg_112, lg_114, lg_242, lg_247, lg_249, \
                         lg_437, lg_442, lg_444 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_140[k] = -f_50 * lg_32[k]
                   - f_50 * lg_37[k]
                   + f_51 * lg_39[k]
                   + f_48 * lg_107[k]
                   + f_48 * lg_112[k]
                   - f_49 * lg_114[k]
                   - f_46 * lg_242[k]
                   - f_46 * lg_247[k]
                   + f_47 * lg_249[k]
                   + f_44 * lg_437[k]
                   + f_44 * lg_442[k]
                   - f_45 * lg_444[k];
    }

#pragma omp simd aligned(lg_30, lg_35, lg_40, lg_42, lg_105, lg_110, lg_115, lg_117, lg_240, \
                         lg_245, lg_250, lg_252, lg_435, lg_440, lg_445, \
                         lg_447 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_141[k] = -f_70 * lg_30[k]
                   + f_71 * lg_35[k]
                   + f_70 * lg_40[k]
                   - f_71 * lg_42[k]
                   + f_68 * lg_105[k]
                   - f_69 * lg_110[k]
                   - f_68 * lg_115[k]
                   + f_69 * lg_117[k]
                   - f_66 * lg_240[k]
                   + f_67 * lg_245[k]
                   + f_66 * lg_250[k]
                   - f_67 * lg_252[k]
                   + f_65 * lg_435[k]
                   - f_42 * lg_440[k]
                   - f_65 * lg_445[k]
                   + f_42 * lg_447[k];
    }

#pragma omp simd aligned(lg_32, lg_37, lg_107, lg_112, lg_242, lg_247, lg_437, \
                         lg_442 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_142[k] = f_39 * lg_32[k]
                   - f_38 * lg_37[k]
                   - f_33 * lg_107[k]
                   + f_37 * lg_112[k]
                   + f_36 * lg_242[k]
                   - f_35 * lg_247[k]
                   - f_34 * lg_437[k]
                   + f_33 * lg_442[k];
    }

#pragma omp simd aligned(lg_30, lg_33, lg_40, lg_105, lg_108, lg_115, lg_240, lg_243, lg_250, \
                         lg_435, lg_438, lg_445 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_143[k] = f_78 * lg_30[k]
                   - f_79 * lg_33[k]
                   + f_78 * lg_40[k]
                   - f_76 * lg_105[k]
                   + f_77 * lg_108[k]
                   - f_76 * lg_115[k]
                   + f_74 * lg_240[k]
                   - f_75 * lg_243[k]
                   + f_74 * lg_250[k]
                   - f_72 * lg_435[k]
                   + f_73 * lg_438[k]
                   - f_72 * lg_445[k];
    }

#pragma omp simd aligned(lg_1, lg_6, lg_46, lg_51, lg_151, lg_156, lg_316, lg_321, lg_541, \
                         lg_546 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_144[k] = f_78 * lg_1[k]
                   - f_78 * lg_6[k]
                   - f_30 * lg_46[k]
                   + f_30 * lg_51[k]
                   + f_662 * lg_151[k]
                   - f_662 * lg_156[k]
                   - f_30 * lg_316[k]
                   + f_30 * lg_321[k]
                   + f_78 * lg_541[k]
                   - f_78 * lg_546[k];
    }

#pragma omp simd aligned(lg_4, lg_11, lg_49, lg_56, lg_154, lg_161, lg_319, lg_326, lg_544, \
                         lg_551 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_145[k] = f_663 * lg_4[k]
                   - f_664 * lg_11[k]
                   - f_33 * lg_49[k]
                   + f_34 * lg_56[k]
                   + f_665 * lg_154[k]
                   - f_666 * lg_161[k]
                   - f_33 * lg_319[k]
                   + f_34 * lg_326[k]
                   + f_663 * lg_544[k]
                   - f_664 * lg_551[k];
    }

#pragma omp simd aligned(lg_1, lg_6, lg_8, lg_46, lg_51, lg_53, lg_151, lg_156, lg_158, \
                         lg_316, lg_321, lg_323, lg_541, lg_546, \
                         lg_548 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_146[k] = -f_667 * lg_1[k]
                   - f_667 * lg_6[k]
                   + f_668 * lg_8[k]
                   + f_24 * lg_46[k]
                   + f_24 * lg_51[k]
                   - f_25 * lg_53[k]
                   - f_66 * lg_151[k]
                   - f_66 * lg_156[k]
                   + f_67 * lg_158[k]
                   + f_24 * lg_316[k]
                   + f_24 * lg_321[k]
                   - f_25 * lg_323[k]
                   - f_667 * lg_541[k]
                   - f_667 * lg_546[k]
                   + f_668 * lg_548[k];
    }

#pragma omp simd aligned(lg_4, lg_11, lg_13, lg_49, lg_56, lg_58, lg_154, lg_161, lg_163, \
                         lg_319, lg_326, lg_328, lg_544, lg_551, \
                         lg_553 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_147[k] = -f_669 * lg_4[k]
                   - f_669 * lg_11[k]
                   + f_670 * lg_13[k]
                   + f_44 * lg_49[k]
                   + f_44 * lg_56[k]
                   - f_45 * lg_58[k]
                   - f_671 * lg_154[k]
                   - f_671 * lg_161[k]
                   + f_672 * lg_163[k]
                   + f_44 * lg_319[k]
                   + f_44 * lg_326[k]
                   - f_45 * lg_328[k]
                   - f_669 * lg_544[k]
                   - f_669 * lg_551[k]
                   + f_670 * lg_553[k];
    }

#pragma omp simd aligned(lg_0, lg_3, lg_5, lg_10, lg_12, lg_14, lg_45, lg_48, lg_50, lg_55, \
                         lg_57, lg_59, lg_150, lg_153, lg_155, lg_160, lg_162, lg_164, lg_315, \
                         lg_318, lg_320, lg_325, lg_327, lg_329, lg_540, lg_543, lg_545, \
                         lg_550, lg_552, lg_554 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_148[k] = f_673 * lg_0[k]
                   + f_674 * lg_3[k]
                   - f_14 * lg_5[k]
                   + f_673 * lg_10[k]
                   - f_14 * lg_12[k]
                   + f_675 * lg_14[k]
                   - f_52 * lg_45[k]
                   - f_18 * lg_48[k]
                   + f_53 * lg_50[k]
                   - f_52 * lg_55[k]
                   + f_53 * lg_57[k]
                   - f_54 * lg_59[k]
                   + f_676 * lg_150[k]
                   + f_55 * lg_153[k]
                   - f_677 * lg_155[k]
                   + f_676 * lg_160[k]
                   - f_677 * lg_162[k]
                   + f_137 * lg_164[k]
                   - f_52 * lg_315[k]
                   - f_18 * lg_318[k]
                   + f_53 * lg_320[k]
                   - f_52 * lg_325[k]
                   + f_53 * lg_327[k]
                   - f_54 * lg_329[k]
                   + f_673 * lg_540[k]
                   + f_674 * lg_543[k]
                   - f_14 * lg_545[k]
                   + f_673 * lg_550[k]
                   - f_14 * lg_552[k]
                   + f_675 * lg_554[k];
    }

#pragma omp simd aligned(lg_2, lg_7, lg_9, lg_47, lg_52, lg_54, lg_152, lg_157, lg_159, \
                         lg_317, lg_322, lg_324, lg_542, lg_547, \
                         lg_549 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_149[k] = -f_669 * lg_2[k]
                   - f_669 * lg_7[k]
                   + f_670 * lg_9[k]
                   + f_44 * lg_47[k]
                   + f_44 * lg_52[k]
                   - f_45 * lg_54[k]
                   - f_671 * lg_152[k]
                   - f_671 * lg_157[k]
                   + f_672 * lg_159[k]
                   + f_44 * lg_317[k]
                   + f_44 * lg_322[k]
                   - f_45 * lg_324[k]
                   - f_669 * lg_542[k]
                   - f_669 * lg_547[k]
                   + f_670 * lg_549[k];
    }

#pragma omp simd aligned(lg_0, lg_5, lg_10, lg_12, lg_45, lg_50, lg_55, lg_57, lg_150, lg_155, \
                         lg_160, lg_162, lg_315, lg_320, lg_325, lg_327, lg_540, lg_545, \
                         lg_550, lg_552 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_150[k] = -f_678 * lg_0[k]
                   + f_679 * lg_5[k]
                   + f_678 * lg_10[k]
                   - f_679 * lg_12[k]
                   + f_65 * lg_45[k]
                   - f_42 * lg_50[k]
                   - f_65 * lg_55[k]
                   + f_42 * lg_57[k]
                   - f_680 * lg_150[k]
                   + f_681 * lg_155[k]
                   + f_680 * lg_160[k]
                   - f_681 * lg_162[k]
                   + f_65 * lg_315[k]
                   - f_42 * lg_320[k]
                   - f_65 * lg_325[k]
                   + f_42 * lg_327[k]
                   - f_678 * lg_540[k]
                   + f_679 * lg_545[k]
                   + f_678 * lg_550[k]
                   - f_679 * lg_552[k];
    }

#pragma omp simd aligned(lg_2, lg_7, lg_47, lg_52, lg_152, lg_157, lg_317, lg_322, lg_542, \
                         lg_547 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_151[k] = f_664 * lg_2[k]
                   - f_663 * lg_7[k]
                   - f_34 * lg_47[k]
                   + f_33 * lg_52[k]
                   + f_666 * lg_152[k]
                   - f_665 * lg_157[k]
                   - f_34 * lg_317[k]
                   + f_33 * lg_322[k]
                   + f_664 * lg_542[k]
                   - f_663 * lg_547[k];
    }

#pragma omp simd aligned(lg_0, lg_3, lg_10, lg_45, lg_48, lg_55, lg_150, lg_153, lg_160, \
                         lg_315, lg_318, lg_325, lg_540, lg_543, \
                         lg_550 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_152[k] = f_682 * lg_0[k]
                   - f_683 * lg_3[k]
                   + f_682 * lg_10[k]
                   - f_72 * lg_45[k]
                   + f_73 * lg_48[k]
                   - f_72 * lg_55[k]
                   + f_684 * lg_150[k]
                   - f_685 * lg_153[k]
                   + f_684 * lg_160[k]
                   - f_72 * lg_315[k]
                   + f_73 * lg_318[k]
                   - f_72 * lg_325[k]
                   + f_682 * lg_540[k]
                   - f_683 * lg_543[k]
                   + f_682 * lg_550[k];
    }
}

}  // namespace simdtrf
