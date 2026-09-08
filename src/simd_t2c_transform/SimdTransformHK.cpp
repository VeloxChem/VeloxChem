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


#include "SimdTransformHK.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_hk(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t hk,
             const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 0.205078125 * std::sqrt(6006.0);
    const auto f_1 = 1.025390625 * std::sqrt(6006.0);
    const auto f_2 = 0.615234375 * std::sqrt(6006.0);
    const auto f_3 = 0.029296875 * std::sqrt(6006.0);
    const auto f_4 = 0.41015625 * std::sqrt(6006.0);
    const auto f_5 = 2.05078125 * std::sqrt(6006.0);
    const auto f_6 = 1.23046875 * std::sqrt(6006.0);
    const auto f_7 = 0.05859375 * std::sqrt(6006.0);
    const auto f_8 = 0.041015625 * std::sqrt(6006.0);
    const auto f_9 = 0.123046875 * std::sqrt(6006.0);
    const auto f_10 = 0.005859375 * std::sqrt(6006.0);
    const auto f_11 = 2.4609375 * std::sqrt(429.0);
    const auto f_12 = 8.203125 * std::sqrt(429.0);
    const auto f_13 = 4.921875 * std::sqrt(429.0);
    const auto f_14 = 16.40625 * std::sqrt(429.0);
    const auto f_15 = 0.4921875 * std::sqrt(429.0);
    const auto f_16 = 1.640625 * std::sqrt(429.0);
    const auto f_17 = 1.025390625 * std::sqrt(66.0);
    const auto f_18 = 12.3046875 * std::sqrt(66.0);
    const auto f_19 = 1.845703125 * std::sqrt(66.0);
    const auto f_20 = 24.609375 * std::sqrt(66.0);
    const auto f_21 = 0.205078125 * std::sqrt(66.0);
    const auto f_22 = 2.4609375 * std::sqrt(66.0);
    const auto f_23 = 2.05078125 * std::sqrt(66.0);
    const auto f_24 = 3.69140625 * std::sqrt(66.0);
    const auto f_25 = 49.21875 * std::sqrt(66.0);
    const auto f_26 = 0.41015625 * std::sqrt(66.0);
    const auto f_27 = 4.921875 * std::sqrt(66.0);
    const auto f_28 = 0.369140625 * std::sqrt(66.0);
    const auto f_29 = 0.041015625 * std::sqrt(66.0);
    const auto f_30 = 0.4921875 * std::sqrt(66.0);
    const auto f_31 = 16.40625 * std::sqrt(66.0);
    const auto f_32 = 9.84375 * std::sqrt(66.0);
    const auto f_33 = 32.8125 * std::sqrt(66.0);
    const auto f_34 = 0.984375 * std::sqrt(66.0);
    const auto f_35 = 3.28125 * std::sqrt(66.0);
    const auto f_36 = 1.845703125 * std::sqrt(6.0);
    const auto f_37 = 3.076171875 * std::sqrt(6.0);
    const auto f_38 = 36.9140625 * std::sqrt(6.0);
    const auto f_39 = 0.615234375 * std::sqrt(6.0);
    const auto f_40 = 24.609375 * std::sqrt(6.0);
    const auto f_41 = 49.21875 * std::sqrt(6.0);
    const auto f_42 = 12.3046875 * std::sqrt(6.0);
    const auto f_43 = 16.40625 * std::sqrt(6.0);
    const auto f_44 = 3.69140625 * std::sqrt(6.0);
    const auto f_45 = 6.15234375 * std::sqrt(6.0);
    const auto f_46 = 73.828125 * std::sqrt(6.0);
    const auto f_47 = 1.23046875 * std::sqrt(6.0);
    const auto f_48 = 98.4375 * std::sqrt(6.0);
    const auto f_49 = 32.8125 * std::sqrt(6.0);
    const auto f_50 = 0.369140625 * std::sqrt(6.0);
    const auto f_51 = 7.3828125 * std::sqrt(6.0);
    const auto f_52 = 0.123046875 * std::sqrt(6.0);
    const auto f_53 = 4.921875 * std::sqrt(6.0);
    const auto f_54 = 9.84375 * std::sqrt(6.0);
    const auto f_55 = 2.4609375 * std::sqrt(6.0);
    const auto f_56 = 3.28125 * std::sqrt(6.0);
    const auto f_57 = 12.3046875 * std::sqrt(3.0);
    const auto f_58 = 24.609375 * std::sqrt(3.0);
    const auto f_59 = 65.625 * std::sqrt(3.0);
    const auto f_60 = 39.375 * std::sqrt(3.0);
    const auto f_61 = 49.21875 * std::sqrt(3.0);
    const auto f_62 = 131.25 * std::sqrt(3.0);
    const auto f_63 = 78.75 * std::sqrt(3.0);
    const auto f_64 = 2.4609375 * std::sqrt(3.0);
    const auto f_65 = 4.921875 * std::sqrt(3.0);
    const auto f_66 = 13.125 * std::sqrt(3.0);
    const auto f_67 = 7.875 * std::sqrt(3.0);
    const auto f_68 = 1.025390625 * std::sqrt(2.0);
    const auto f_69 = 3.076171875 * std::sqrt(2.0);
    const auto f_70 = 24.609375 * std::sqrt(2.0);
    const auto f_71 = 49.21875 * std::sqrt(2.0);
    const auto f_72 = 13.125 * std::sqrt(2.0);
    const auto f_73 = 2.05078125 * std::sqrt(2.0);
    const auto f_74 = 6.15234375 * std::sqrt(2.0);
    const auto f_75 = 98.4375 * std::sqrt(2.0);
    const auto f_76 = 26.25 * std::sqrt(2.0);
    const auto f_77 = 0.205078125 * std::sqrt(2.0);
    const auto f_78 = 0.615234375 * std::sqrt(2.0);
    const auto f_79 = 4.921875 * std::sqrt(2.0);
    const auto f_80 = 9.84375 * std::sqrt(2.0);
    const auto f_81 = 2.625 * std::sqrt(2.0);
    const auto f_82 = 2.05078125 * std::sqrt(14.0);
    const auto f_83 = 6.15234375 * std::sqrt(14.0);
    const auto f_84 = 12.3046875 * std::sqrt(14.0);
    const auto f_85 = 24.609375 * std::sqrt(14.0);
    const auto f_86 = 9.84375 * std::sqrt(14.0);
    const auto f_87 = 0.9375 * std::sqrt(14.0);
    const auto f_88 = 4.1015625 * std::sqrt(14.0);
    const auto f_89 = 49.21875 * std::sqrt(14.0);
    const auto f_90 = 19.6875 * std::sqrt(14.0);
    const auto f_91 = 1.875 * std::sqrt(14.0);
    const auto f_92 = 0.41015625 * std::sqrt(14.0);
    const auto f_93 = 1.23046875 * std::sqrt(14.0);
    const auto f_94 = 2.4609375 * std::sqrt(14.0);
    const auto f_95 = 4.921875 * std::sqrt(14.0);
    const auto f_96 = 1.96875 * std::sqrt(14.0);
    const auto f_97 = 0.1875 * std::sqrt(14.0);
    const auto f_98 = 6.15234375 * std::sqrt(3.0);
    const auto f_99 = 32.8125 * std::sqrt(3.0);
    const auto f_100 = 19.6875 * std::sqrt(3.0);
    const auto f_101 = 1.23046875 * std::sqrt(3.0);
    const auto f_102 = 6.5625 * std::sqrt(3.0);
    const auto f_103 = 3.9375 * std::sqrt(3.0);
    const auto f_104 = 1.23046875 * std::sqrt(66.0);
    const auto f_105 = 6.15234375 * std::sqrt(66.0);
    const auto f_106 = 4.1015625 * std::sqrt(66.0);
    const auto f_107 = 8.203125 * std::sqrt(66.0);
    const auto f_108 = 0.24609375 * std::sqrt(66.0);
    const auto f_109 = 0.8203125 * std::sqrt(66.0);
    const auto f_110 = 0.41015625 * std::sqrt(429.0);
    const auto f_111 = 6.15234375 * std::sqrt(429.0);
    const auto f_112 = 0.8203125 * std::sqrt(429.0);
    const auto f_113 = 12.3046875 * std::sqrt(429.0);
    const auto f_114 = 0.08203125 * std::sqrt(429.0);
    const auto f_115 = 1.23046875 * std::sqrt(429.0);
    const auto f_116 = 0.328125 * std::sqrt(15015.0);
    const auto f_117 = 1.640625 * std::sqrt(15015.0);
    const auto f_118 = 0.984375 * std::sqrt(15015.0);
    const auto f_119 = 0.046875 * std::sqrt(15015.0);
    const auto f_120 = 1.96875 * std::sqrt(4290.0);
    const auto f_121 = 6.5625 * std::sqrt(4290.0);
    const auto f_122 = 1.640625 * std::sqrt(165.0);
    const auto f_123 = 19.6875 * std::sqrt(165.0);
    const auto f_124 = 2.953125 * std::sqrt(165.0);
    const auto f_125 = 39.375 * std::sqrt(165.0);
    const auto f_126 = 0.328125 * std::sqrt(165.0);
    const auto f_127 = 3.9375 * std::sqrt(165.0);
    const auto f_128 = 7.875 * std::sqrt(165.0);
    const auto f_129 = 26.25 * std::sqrt(165.0);
    const auto f_130 = 2.953125 * std::sqrt(15.0);
    const auto f_131 = 4.921875 * std::sqrt(15.0);
    const auto f_132 = 59.0625 * std::sqrt(15.0);
    const auto f_133 = 0.984375 * std::sqrt(15.0);
    const auto f_134 = 39.375 * std::sqrt(15.0);
    const auto f_135 = 78.75 * std::sqrt(15.0);
    const auto f_136 = 19.6875 * std::sqrt(15.0);
    const auto f_137 = 26.25 * std::sqrt(15.0);
    const auto f_138 = 9.84375 * std::sqrt(30.0);
    const auto f_139 = 19.6875 * std::sqrt(30.0);
    const auto f_140 = 52.5 * std::sqrt(30.0);
    const auto f_141 = 31.5 * std::sqrt(30.0);
    const auto f_142 = 1.640625 * std::sqrt(5.0);
    const auto f_143 = 4.921875 * std::sqrt(5.0);
    const auto f_144 = 39.375 * std::sqrt(5.0);
    const auto f_145 = 78.75 * std::sqrt(5.0);
    const auto f_146 = 21.0 * std::sqrt(5.0);
    const auto f_147 = 3.28125 * std::sqrt(35.0);
    const auto f_148 = 9.84375 * std::sqrt(35.0);
    const auto f_149 = 19.6875 * std::sqrt(35.0);
    const auto f_150 = 39.375 * std::sqrt(35.0);
    const auto f_151 = 15.75 * std::sqrt(35.0);
    const auto f_152 = 1.5 * std::sqrt(35.0);
    const auto f_153 = 4.921875 * std::sqrt(30.0);
    const auto f_154 = 26.25 * std::sqrt(30.0);
    const auto f_155 = 15.75 * std::sqrt(30.0);
    const auto f_156 = 1.96875 * std::sqrt(165.0);
    const auto f_157 = 9.84375 * std::sqrt(165.0);
    const auto f_158 = 6.5625 * std::sqrt(165.0);
    const auto f_159 = 0.328125 * std::sqrt(4290.0);
    const auto f_160 = 4.921875 * std::sqrt(4290.0);
    const auto f_161 = 0.041015625 * std::sqrt(30030.0);
    const auto f_162 = 0.205078125 * std::sqrt(30030.0);
    const auto f_163 = 0.123046875 * std::sqrt(30030.0);
    const auto f_164 = 0.005859375 * std::sqrt(30030.0);
    const auto f_165 = 0.02734375 * std::sqrt(30030.0);
    const auto f_166 = 0.13671875 * std::sqrt(30030.0);
    const auto f_167 = 0.08203125 * std::sqrt(30030.0);
    const auto f_168 = 0.00390625 * std::sqrt(30030.0);
    const auto f_169 = 0.328125 * std::sqrt(30030.0);
    const auto f_170 = 1.640625 * std::sqrt(30030.0);
    const auto f_171 = 0.984375 * std::sqrt(30030.0);
    const auto f_172 = 0.046875 * std::sqrt(30030.0);
    const auto f_173 = 0.013671875 * std::sqrt(30030.0);
    const auto f_174 = 0.068359375 * std::sqrt(30030.0);
    const auto f_175 = 0.001953125 * std::sqrt(30030.0);
    const auto f_176 = 0.109375 * std::sqrt(30030.0);
    const auto f_177 = 0.546875 * std::sqrt(30030.0);
    const auto f_178 = 0.015625 * std::sqrt(30030.0);
    const auto f_179 = 0.4921875 * std::sqrt(2145.0);
    const auto f_180 = 1.640625 * std::sqrt(2145.0);
    const auto f_181 = 0.328125 * std::sqrt(2145.0);
    const auto f_182 = 1.09375 * std::sqrt(2145.0);
    const auto f_183 = 3.9375 * std::sqrt(2145.0);
    const auto f_184 = 13.125 * std::sqrt(2145.0);
    const auto f_185 = 0.1640625 * std::sqrt(2145.0);
    const auto f_186 = 0.546875 * std::sqrt(2145.0);
    const auto f_187 = 1.3125 * std::sqrt(2145.0);
    const auto f_188 = 4.375 * std::sqrt(2145.0);
    const auto f_189 = 0.205078125 * std::sqrt(330.0);
    const auto f_190 = 2.4609375 * std::sqrt(330.0);
    const auto f_191 = 0.369140625 * std::sqrt(330.0);
    const auto f_192 = 4.921875 * std::sqrt(330.0);
    const auto f_193 = 0.041015625 * std::sqrt(330.0);
    const auto f_194 = 0.4921875 * std::sqrt(330.0);
    const auto f_195 = 0.13671875 * std::sqrt(330.0);
    const auto f_196 = 1.640625 * std::sqrt(330.0);
    const auto f_197 = 0.24609375 * std::sqrt(330.0);
    const auto f_198 = 3.28125 * std::sqrt(330.0);
    const auto f_199 = 0.02734375 * std::sqrt(330.0);
    const auto f_200 = 0.328125 * std::sqrt(330.0);
    const auto f_201 = 19.6875 * std::sqrt(330.0);
    const auto f_202 = 2.953125 * std::sqrt(330.0);
    const auto f_203 = 39.375 * std::sqrt(330.0);
    const auto f_204 = 3.9375 * std::sqrt(330.0);
    const auto f_205 = 0.068359375 * std::sqrt(330.0);
    const auto f_206 = 0.8203125 * std::sqrt(330.0);
    const auto f_207 = 0.123046875 * std::sqrt(330.0);
    const auto f_208 = 0.013671875 * std::sqrt(330.0);
    const auto f_209 = 0.1640625 * std::sqrt(330.0);
    const auto f_210 = 0.546875 * std::sqrt(330.0);
    const auto f_211 = 6.5625 * std::sqrt(330.0);
    const auto f_212 = 0.984375 * std::sqrt(330.0);
    const auto f_213 = 13.125 * std::sqrt(330.0);
    const auto f_214 = 0.109375 * std::sqrt(330.0);
    const auto f_215 = 1.3125 * std::sqrt(330.0);
    const auto f_216 = 0.65625 * std::sqrt(330.0);
    const auto f_217 = 2.1875 * std::sqrt(330.0);
    const auto f_218 = 7.875 * std::sqrt(330.0);
    const auto f_219 = 26.25 * std::sqrt(330.0);
    const auto f_220 = 1.09375 * std::sqrt(330.0);
    const auto f_221 = 2.625 * std::sqrt(330.0);
    const auto f_222 = 8.75 * std::sqrt(330.0);
    const auto f_223 = 0.369140625 * std::sqrt(30.0);
    const auto f_224 = 0.615234375 * std::sqrt(30.0);
    const auto f_225 = 7.3828125 * std::sqrt(30.0);
    const auto f_226 = 0.123046875 * std::sqrt(30.0);
    const auto f_227 = 2.4609375 * std::sqrt(30.0);
    const auto f_228 = 3.28125 * std::sqrt(30.0);
    const auto f_229 = 0.24609375 * std::sqrt(30.0);
    const auto f_230 = 0.41015625 * std::sqrt(30.0);
    const auto f_231 = 0.08203125 * std::sqrt(30.0);
    const auto f_232 = 6.5625 * std::sqrt(30.0);
    const auto f_233 = 1.640625 * std::sqrt(30.0);
    const auto f_234 = 2.1875 * std::sqrt(30.0);
    const auto f_235 = 2.953125 * std::sqrt(30.0);
    const auto f_236 = 59.0625 * std::sqrt(30.0);
    const auto f_237 = 0.984375 * std::sqrt(30.0);
    const auto f_238 = 39.375 * std::sqrt(30.0);
    const auto f_239 = 78.75 * std::sqrt(30.0);
    const auto f_240 = 0.205078125 * std::sqrt(30.0);
    const auto f_241 = 0.041015625 * std::sqrt(30.0);
    const auto f_242 = 0.8203125 * std::sqrt(30.0);
    const auto f_243 = 1.09375 * std::sqrt(30.0);
    const auto f_244 = 0.328125 * std::sqrt(30.0);
    const auto f_245 = 13.125 * std::sqrt(30.0);
    const auto f_246 = 8.75 * std::sqrt(30.0);
    const auto f_247 = 2.4609375 * std::sqrt(15.0);
    const auto f_248 = 13.125 * std::sqrt(15.0);
    const auto f_249 = 7.875 * std::sqrt(15.0);
    const auto f_250 = 1.640625 * std::sqrt(15.0);
    const auto f_251 = 3.28125 * std::sqrt(15.0);
    const auto f_252 = 8.75 * std::sqrt(15.0);
    const auto f_253 = 5.25 * std::sqrt(15.0);
    const auto f_254 = 105.0 * std::sqrt(15.0);
    const auto f_255 = 63.0 * std::sqrt(15.0);
    const auto f_256 = 0.8203125 * std::sqrt(15.0);
    const auto f_257 = 4.375 * std::sqrt(15.0);
    const auto f_258 = 2.625 * std::sqrt(15.0);
    const auto f_259 = 6.5625 * std::sqrt(15.0);
    const auto f_260 = 35.0 * std::sqrt(15.0);
    const auto f_261 = 21.0 * std::sqrt(15.0);
    const auto f_262 = 0.205078125 * std::sqrt(10.0);
    const auto f_263 = 0.615234375 * std::sqrt(10.0);
    const auto f_264 = 4.921875 * std::sqrt(10.0);
    const auto f_265 = 9.84375 * std::sqrt(10.0);
    const auto f_266 = 2.625 * std::sqrt(10.0);
    const auto f_267 = 0.13671875 * std::sqrt(10.0);
    const auto f_268 = 0.41015625 * std::sqrt(10.0);
    const auto f_269 = 3.28125 * std::sqrt(10.0);
    const auto f_270 = 6.5625 * std::sqrt(10.0);
    const auto f_271 = 1.75 * std::sqrt(10.0);
    const auto f_272 = 1.640625 * std::sqrt(10.0);
    const auto f_273 = 39.375 * std::sqrt(10.0);
    const auto f_274 = 78.75 * std::sqrt(10.0);
    const auto f_275 = 21.0 * std::sqrt(10.0);
    const auto f_276 = 0.068359375 * std::sqrt(10.0);
    const auto f_277 = 0.875 * std::sqrt(10.0);
    const auto f_278 = 0.546875 * std::sqrt(10.0);
    const auto f_279 = 13.125 * std::sqrt(10.0);
    const auto f_280 = 26.25 * std::sqrt(10.0);
    const auto f_281 = 7.0 * std::sqrt(10.0);
    const auto f_282 = 0.41015625 * std::sqrt(70.0);
    const auto f_283 = 1.23046875 * std::sqrt(70.0);
    const auto f_284 = 2.4609375 * std::sqrt(70.0);
    const auto f_285 = 4.921875 * std::sqrt(70.0);
    const auto f_286 = 1.96875 * std::sqrt(70.0);
    const auto f_287 = 0.1875 * std::sqrt(70.0);
    const auto f_288 = 0.2734375 * std::sqrt(70.0);
    const auto f_289 = 0.8203125 * std::sqrt(70.0);
    const auto f_290 = 1.640625 * std::sqrt(70.0);
    const auto f_291 = 3.28125 * std::sqrt(70.0);
    const auto f_292 = 1.3125 * std::sqrt(70.0);
    const auto f_293 = 0.125 * std::sqrt(70.0);
    const auto f_294 = 9.84375 * std::sqrt(70.0);
    const auto f_295 = 19.6875 * std::sqrt(70.0);
    const auto f_296 = 39.375 * std::sqrt(70.0);
    const auto f_297 = 15.75 * std::sqrt(70.0);
    const auto f_298 = 1.5 * std::sqrt(70.0);
    const auto f_299 = 0.13671875 * std::sqrt(70.0);
    const auto f_300 = 0.65625 * std::sqrt(70.0);
    const auto f_301 = 0.0625 * std::sqrt(70.0);
    const auto f_302 = 1.09375 * std::sqrt(70.0);
    const auto f_303 = 6.5625 * std::sqrt(70.0);
    const auto f_304 = 13.125 * std::sqrt(70.0);
    const auto f_305 = 5.25 * std::sqrt(70.0);
    const auto f_306 = 0.5 * std::sqrt(70.0);
    const auto f_307 = 1.23046875 * std::sqrt(15.0);
    const auto f_308 = 3.9375 * std::sqrt(15.0);
    const auto f_309 = 9.84375 * std::sqrt(15.0);
    const auto f_310 = 52.5 * std::sqrt(15.0);
    const auto f_311 = 31.5 * std::sqrt(15.0);
    const auto f_312 = 0.41015625 * std::sqrt(15.0);
    const auto f_313 = 2.1875 * std::sqrt(15.0);
    const auto f_314 = 1.3125 * std::sqrt(15.0);
    const auto f_315 = 17.5 * std::sqrt(15.0);
    const auto f_316 = 10.5 * std::sqrt(15.0);
    const auto f_317 = 1.23046875 * std::sqrt(330.0);
    const auto f_318 = 1.96875 * std::sqrt(330.0);
    const auto f_319 = 9.84375 * std::sqrt(330.0);
    const auto f_320 = 0.08203125 * std::sqrt(330.0);
    const auto f_321 = 0.41015625 * std::sqrt(330.0);
    const auto f_322 = 0.2734375 * std::sqrt(330.0);
    const auto f_323 = 0.08203125 * std::sqrt(2145.0);
    const auto f_324 = 1.23046875 * std::sqrt(2145.0);
    const auto f_325 = 0.0546875 * std::sqrt(2145.0);
    const auto f_326 = 0.8203125 * std::sqrt(2145.0);
    const auto f_327 = 0.65625 * std::sqrt(2145.0);
    const auto f_328 = 9.84375 * std::sqrt(2145.0);
    const auto f_329 = 0.02734375 * std::sqrt(2145.0);
    const auto f_330 = 0.41015625 * std::sqrt(2145.0);
    const auto f_331 = 0.21875 * std::sqrt(2145.0);
    const auto f_332 = 3.28125 * std::sqrt(2145.0);
    const auto f_333 = 0.328125 * std::sqrt(5005.0);
    const auto f_334 = 1.640625 * std::sqrt(5005.0);
    const auto f_335 = 0.984375 * std::sqrt(5005.0);
    const auto f_336 = 0.046875 * std::sqrt(5005.0);
    const auto f_337 = 0.65625 * std::sqrt(5005.0);
    const auto f_338 = 3.28125 * std::sqrt(5005.0);
    const auto f_339 = 1.96875 * std::sqrt(5005.0);
    const auto f_340 = 0.09375 * std::sqrt(5005.0);
    const auto f_341 = 1.96875 * std::sqrt(1430.0);
    const auto f_342 = 6.5625 * std::sqrt(1430.0);
    const auto f_343 = 3.9375 * std::sqrt(1430.0);
    const auto f_344 = 13.125 * std::sqrt(1430.0);
    const auto f_345 = 1.640625 * std::sqrt(55.0);
    const auto f_346 = 19.6875 * std::sqrt(55.0);
    const auto f_347 = 2.953125 * std::sqrt(55.0);
    const auto f_348 = 39.375 * std::sqrt(55.0);
    const auto f_349 = 0.328125 * std::sqrt(55.0);
    const auto f_350 = 3.9375 * std::sqrt(55.0);
    const auto f_351 = 3.28125 * std::sqrt(55.0);
    const auto f_352 = 5.90625 * std::sqrt(55.0);
    const auto f_353 = 78.75 * std::sqrt(55.0);
    const auto f_354 = 0.65625 * std::sqrt(55.0);
    const auto f_355 = 7.875 * std::sqrt(55.0);
    const auto f_356 = 26.25 * std::sqrt(55.0);
    const auto f_357 = 15.75 * std::sqrt(55.0);
    const auto f_358 = 52.5 * std::sqrt(55.0);
    const auto f_359 = 2.953125 * std::sqrt(5.0);
    const auto f_360 = 59.0625 * std::sqrt(5.0);
    const auto f_361 = 0.984375 * std::sqrt(5.0);
    const auto f_362 = 19.6875 * std::sqrt(5.0);
    const auto f_363 = 26.25 * std::sqrt(5.0);
    const auto f_364 = 5.90625 * std::sqrt(5.0);
    const auto f_365 = 9.84375 * std::sqrt(5.0);
    const auto f_366 = 118.125 * std::sqrt(5.0);
    const auto f_367 = 1.96875 * std::sqrt(5.0);
    const auto f_368 = 157.5 * std::sqrt(5.0);
    const auto f_369 = 52.5 * std::sqrt(5.0);
    const auto f_370 = 19.6875 * std::sqrt(10.0);
    const auto f_371 = 52.5 * std::sqrt(10.0);
    const auto f_372 = 31.5 * std::sqrt(10.0);
    const auto f_373 = 105.0 * std::sqrt(10.0);
    const auto f_374 = 63.0 * std::sqrt(10.0);
    const auto f_375 = 0.546875 * std::sqrt(15.0);
    const auto f_376 = 7.0 * std::sqrt(15.0);
    const auto f_377 = 1.09375 * std::sqrt(15.0);
    const auto f_378 = 14.0 * std::sqrt(15.0);
    const auto f_379 = 1.09375 * std::sqrt(105.0);
    const auto f_380 = 3.28125 * std::sqrt(105.0);
    const auto f_381 = 6.5625 * std::sqrt(105.0);
    const auto f_382 = 13.125 * std::sqrt(105.0);
    const auto f_383 = 5.25 * std::sqrt(105.0);
    const auto f_384 = 0.5 * std::sqrt(105.0);
    const auto f_385 = 2.1875 * std::sqrt(105.0);
    const auto f_386 = 26.25 * std::sqrt(105.0);
    const auto f_387 = 10.5 * std::sqrt(105.0);
    const auto f_388 = std::sqrt(105.0);
    const auto f_389 = 15.75 * std::sqrt(10.0);
    const auto f_390 = 1.96875 * std::sqrt(55.0);
    const auto f_391 = 9.84375 * std::sqrt(55.0);
    const auto f_392 = 6.5625 * std::sqrt(55.0);
    const auto f_393 = 13.125 * std::sqrt(55.0);
    const auto f_394 = 0.328125 * std::sqrt(1430.0);
    const auto f_395 = 4.921875 * std::sqrt(1430.0);
    const auto f_396 = 0.65625 * std::sqrt(1430.0);
    const auto f_397 = 9.84375 * std::sqrt(1430.0);
    const auto f_398 = 0.08203125 * std::sqrt(715.0);
    const auto f_399 = 0.41015625 * std::sqrt(715.0);
    const auto f_400 = 0.24609375 * std::sqrt(715.0);
    const auto f_401 = 0.01171875 * std::sqrt(715.0);
    const auto f_402 = 0.1640625 * std::sqrt(715.0);
    const auto f_403 = 0.8203125 * std::sqrt(715.0);
    const auto f_404 = 0.4921875 * std::sqrt(715.0);
    const auto f_405 = 0.0234375 * std::sqrt(715.0);
    const auto f_406 = 0.984375 * std::sqrt(715.0);
    const auto f_407 = 4.921875 * std::sqrt(715.0);
    const auto f_408 = 2.953125 * std::sqrt(715.0);
    const auto f_409 = 0.140625 * std::sqrt(715.0);
    const auto f_410 = 0.65625 * std::sqrt(715.0);
    const auto f_411 = 3.28125 * std::sqrt(715.0);
    const auto f_412 = 1.96875 * std::sqrt(715.0);
    const auto f_413 = 0.09375 * std::sqrt(715.0);
    const auto f_414 = 0.0703125 * std::sqrt(10010.0);
    const auto f_415 = 0.234375 * std::sqrt(10010.0);
    const auto f_416 = 0.140625 * std::sqrt(10010.0);
    const auto f_417 = 0.46875 * std::sqrt(10010.0);
    const auto f_418 = 0.84375 * std::sqrt(10010.0);
    const auto f_419 = 2.8125 * std::sqrt(10010.0);
    const auto f_420 = 0.5625 * std::sqrt(10010.0);
    const auto f_421 = 1.875 * std::sqrt(10010.0);
    const auto f_422 = 0.05859375 * std::sqrt(385.0);
    const auto f_423 = 0.703125 * std::sqrt(385.0);
    const auto f_424 = 0.10546875 * std::sqrt(385.0);
    const auto f_425 = 1.40625 * std::sqrt(385.0);
    const auto f_426 = 0.01171875 * std::sqrt(385.0);
    const auto f_427 = 0.140625 * std::sqrt(385.0);
    const auto f_428 = 0.1171875 * std::sqrt(385.0);
    const auto f_429 = 0.2109375 * std::sqrt(385.0);
    const auto f_430 = 2.8125 * std::sqrt(385.0);
    const auto f_431 = 0.0234375 * std::sqrt(385.0);
    const auto f_432 = 0.28125 * std::sqrt(385.0);
    const auto f_433 = 8.4375 * std::sqrt(385.0);
    const auto f_434 = 1.265625 * std::sqrt(385.0);
    const auto f_435 = 16.875 * std::sqrt(385.0);
    const auto f_436 = 1.6875 * std::sqrt(385.0);
    const auto f_437 = 0.46875 * std::sqrt(385.0);
    const auto f_438 = 5.625 * std::sqrt(385.0);
    const auto f_439 = 0.84375 * std::sqrt(385.0);
    const auto f_440 = 11.25 * std::sqrt(385.0);
    const auto f_441 = 0.09375 * std::sqrt(385.0);
    const auto f_442 = 1.125 * std::sqrt(385.0);
    const auto f_443 = 0.9375 * std::sqrt(385.0);
    const auto f_444 = 0.5625 * std::sqrt(385.0);
    const auto f_445 = 1.875 * std::sqrt(385.0);
    const auto f_446 = 3.375 * std::sqrt(385.0);
    const auto f_447 = 2.25 * std::sqrt(385.0);
    const auto f_448 = 7.5 * std::sqrt(385.0);
    const auto f_449 = 0.10546875 * std::sqrt(35.0);
    const auto f_450 = 0.17578125 * std::sqrt(35.0);
    const auto f_451 = 2.109375 * std::sqrt(35.0);
    const auto f_452 = 0.03515625 * std::sqrt(35.0);
    const auto f_453 = 1.40625 * std::sqrt(35.0);
    const auto f_454 = 2.8125 * std::sqrt(35.0);
    const auto f_455 = 0.703125 * std::sqrt(35.0);
    const auto f_456 = 0.9375 * std::sqrt(35.0);
    const auto f_457 = 0.2109375 * std::sqrt(35.0);
    const auto f_458 = 0.3515625 * std::sqrt(35.0);
    const auto f_459 = 4.21875 * std::sqrt(35.0);
    const auto f_460 = 0.0703125 * std::sqrt(35.0);
    const auto f_461 = 5.625 * std::sqrt(35.0);
    const auto f_462 = 1.875 * std::sqrt(35.0);
    const auto f_463 = 1.265625 * std::sqrt(35.0);
    const auto f_464 = 25.3125 * std::sqrt(35.0);
    const auto f_465 = 0.421875 * std::sqrt(35.0);
    const auto f_466 = 16.875 * std::sqrt(35.0);
    const auto f_467 = 33.75 * std::sqrt(35.0);
    const auto f_468 = 8.4375 * std::sqrt(35.0);
    const auto f_469 = 11.25 * std::sqrt(35.0);
    const auto f_470 = 0.84375 * std::sqrt(35.0);
    const auto f_471 = 0.28125 * std::sqrt(35.0);
    const auto f_472 = 22.5 * std::sqrt(35.0);
    const auto f_473 = 7.5 * std::sqrt(35.0);
    const auto f_474 = 0.3515625 * std::sqrt(70.0);
    const auto f_475 = 0.703125 * std::sqrt(70.0);
    const auto f_476 = 1.875 * std::sqrt(70.0);
    const auto f_477 = 1.125 * std::sqrt(70.0);
    const auto f_478 = 1.40625 * std::sqrt(70.0);
    const auto f_479 = 3.75 * std::sqrt(70.0);
    const auto f_480 = 2.25 * std::sqrt(70.0);
    const auto f_481 = 4.21875 * std::sqrt(70.0);
    const auto f_482 = 8.4375 * std::sqrt(70.0);
    const auto f_483 = 22.5 * std::sqrt(70.0);
    const auto f_484 = 13.5 * std::sqrt(70.0);
    const auto f_485 = 2.8125 * std::sqrt(70.0);
    const auto f_486 = 5.625 * std::sqrt(70.0);
    const auto f_487 = 15.0 * std::sqrt(70.0);
    const auto f_488 = 9.0 * std::sqrt(70.0);
    const auto f_489 = 0.01953125 * std::sqrt(105.0);
    const auto f_490 = 0.05859375 * std::sqrt(105.0);
    const auto f_491 = 0.46875 * std::sqrt(105.0);
    const auto f_492 = 0.9375 * std::sqrt(105.0);
    const auto f_493 = 0.25 * std::sqrt(105.0);
    const auto f_494 = 0.0390625 * std::sqrt(105.0);
    const auto f_495 = 0.1171875 * std::sqrt(105.0);
    const auto f_496 = 1.875 * std::sqrt(105.0);
    const auto f_497 = 0.234375 * std::sqrt(105.0);
    const auto f_498 = 0.703125 * std::sqrt(105.0);
    const auto f_499 = 5.625 * std::sqrt(105.0);
    const auto f_500 = 11.25 * std::sqrt(105.0);
    const auto f_501 = 3.0 * std::sqrt(105.0);
    const auto f_502 = 0.15625 * std::sqrt(105.0);
    const auto f_503 = 3.75 * std::sqrt(105.0);
    const auto f_504 = 7.5 * std::sqrt(105.0);
    const auto f_505 = 2.0 * std::sqrt(105.0);
    const auto f_506 = 0.2734375 * std::sqrt(15.0);
    const auto f_507 = 0.125 * std::sqrt(15.0);
    const auto f_508 = 0.25 * std::sqrt(15.0);
    const auto f_509 = 15.75 * std::sqrt(15.0);
    const auto f_510 = 1.5 * std::sqrt(15.0);
    const auto f_511 = std::sqrt(15.0);
    const auto f_512 = 0.17578125 * std::sqrt(70.0);
    const auto f_513 = 0.9375 * std::sqrt(70.0);
    const auto f_514 = 0.5625 * std::sqrt(70.0);
    const auto f_515 = 2.109375 * std::sqrt(70.0);
    const auto f_516 = 11.25 * std::sqrt(70.0);
    const auto f_517 = 6.75 * std::sqrt(70.0);
    const auto f_518 = 7.5 * std::sqrt(70.0);
    const auto f_519 = 4.5 * std::sqrt(70.0);
    const auto f_520 = 0.0703125 * std::sqrt(385.0);
    const auto f_521 = 0.3515625 * std::sqrt(385.0);
    const auto f_522 = 0.234375 * std::sqrt(385.0);
    const auto f_523 = 4.21875 * std::sqrt(385.0);
    const auto f_524 = 0.01171875 * std::sqrt(10010.0);
    const auto f_525 = 0.17578125 * std::sqrt(10010.0);
    const auto f_526 = 0.0234375 * std::sqrt(10010.0);
    const auto f_527 = 0.3515625 * std::sqrt(10010.0);
    const auto f_528 = 2.109375 * std::sqrt(10010.0);
    const auto f_529 = 0.09375 * std::sqrt(10010.0);
    const auto f_530 = 1.40625 * std::sqrt(10010.0);
    const auto f_531 = 2.05078125 * std::sqrt(429.0);
    const auto f_532 = 0.05859375 * std::sqrt(429.0);
    const auto f_533 = 4.1015625 * std::sqrt(429.0);
    const auto f_534 = 0.1171875 * std::sqrt(429.0);
    const auto f_535 = 1.09375 * std::sqrt(429.0);
    const auto f_536 = 5.46875 * std::sqrt(429.0);
    const auto f_537 = 3.28125 * std::sqrt(429.0);
    const auto f_538 = 0.15625 * std::sqrt(429.0);
    const auto f_539 = 0.21875 * std::sqrt(429.0);
    const auto f_540 = 0.65625 * std::sqrt(429.0);
    const auto f_541 = 0.03125 * std::sqrt(429.0);
    const auto f_542 = 0.3515625 * std::sqrt(6006.0);
    const auto f_543 = 1.171875 * std::sqrt(6006.0);
    const auto f_544 = 0.703125 * std::sqrt(6006.0);
    const auto f_545 = 2.34375 * std::sqrt(6006.0);
    const auto f_546 = 0.9375 * std::sqrt(6006.0);
    const auto f_547 = 3.125 * std::sqrt(6006.0);
    const auto f_548 = 0.1875 * std::sqrt(6006.0);
    const auto f_549 = 0.625 * std::sqrt(6006.0);
    const auto f_550 = 0.29296875 * std::sqrt(231.0);
    const auto f_551 = 3.515625 * std::sqrt(231.0);
    const auto f_552 = 0.52734375 * std::sqrt(231.0);
    const auto f_553 = 7.03125 * std::sqrt(231.0);
    const auto f_554 = 0.05859375 * std::sqrt(231.0);
    const auto f_555 = 0.703125 * std::sqrt(231.0);
    const auto f_556 = 0.5859375 * std::sqrt(231.0);
    const auto f_557 = 1.0546875 * std::sqrt(231.0);
    const auto f_558 = 14.0625 * std::sqrt(231.0);
    const auto f_559 = 0.1171875 * std::sqrt(231.0);
    const auto f_560 = 1.40625 * std::sqrt(231.0);
    const auto f_561 = 0.78125 * std::sqrt(231.0);
    const auto f_562 = 9.375 * std::sqrt(231.0);
    const auto f_563 = 18.75 * std::sqrt(231.0);
    const auto f_564 = 0.15625 * std::sqrt(231.0);
    const auto f_565 = 1.875 * std::sqrt(231.0);
    const auto f_566 = 0.28125 * std::sqrt(231.0);
    const auto f_567 = 3.75 * std::sqrt(231.0);
    const auto f_568 = 0.03125 * std::sqrt(231.0);
    const auto f_569 = 0.375 * std::sqrt(231.0);
    const auto f_570 = 4.6875 * std::sqrt(231.0);
    const auto f_571 = 2.8125 * std::sqrt(231.0);
    const auto f_572 = 12.5 * std::sqrt(231.0);
    const auto f_573 = 0.75 * std::sqrt(231.0);
    const auto f_574 = 2.5 * std::sqrt(231.0);
    const auto f_575 = 0.52734375 * std::sqrt(21.0);
    const auto f_576 = 0.87890625 * std::sqrt(21.0);
    const auto f_577 = 10.546875 * std::sqrt(21.0);
    const auto f_578 = 0.17578125 * std::sqrt(21.0);
    const auto f_579 = 7.03125 * std::sqrt(21.0);
    const auto f_580 = 14.0625 * std::sqrt(21.0);
    const auto f_581 = 3.515625 * std::sqrt(21.0);
    const auto f_582 = 4.6875 * std::sqrt(21.0);
    const auto f_583 = 1.0546875 * std::sqrt(21.0);
    const auto f_584 = 1.7578125 * std::sqrt(21.0);
    const auto f_585 = 21.09375 * std::sqrt(21.0);
    const auto f_586 = 0.3515625 * std::sqrt(21.0);
    const auto f_587 = 28.125 * std::sqrt(21.0);
    const auto f_588 = 9.375 * std::sqrt(21.0);
    const auto f_589 = 1.40625 * std::sqrt(21.0);
    const auto f_590 = 2.34375 * std::sqrt(21.0);
    const auto f_591 = 0.46875 * std::sqrt(21.0);
    const auto f_592 = 18.75 * std::sqrt(21.0);
    const auto f_593 = 37.5 * std::sqrt(21.0);
    const auto f_594 = 12.5 * std::sqrt(21.0);
    const auto f_595 = 0.28125 * std::sqrt(21.0);
    const auto f_596 = 5.625 * std::sqrt(21.0);
    const auto f_597 = 0.09375 * std::sqrt(21.0);
    const auto f_598 = 3.75 * std::sqrt(21.0);
    const auto f_599 = 7.5 * std::sqrt(21.0);
    const auto f_600 = 1.875 * std::sqrt(21.0);
    const auto f_601 = 2.5 * std::sqrt(21.0);
    const auto f_602 = 1.7578125 * std::sqrt(42.0);
    const auto f_603 = 3.515625 * std::sqrt(42.0);
    const auto f_604 = 9.375 * std::sqrt(42.0);
    const auto f_605 = 5.625 * std::sqrt(42.0);
    const auto f_606 = 7.03125 * std::sqrt(42.0);
    const auto f_607 = 18.75 * std::sqrt(42.0);
    const auto f_608 = 11.25 * std::sqrt(42.0);
    const auto f_609 = 4.6875 * std::sqrt(42.0);
    const auto f_610 = 25.0 * std::sqrt(42.0);
    const auto f_611 = 15.0 * std::sqrt(42.0);
    const auto f_612 = 0.9375 * std::sqrt(42.0);
    const auto f_613 = 1.875 * std::sqrt(42.0);
    const auto f_614 = 5.0 * std::sqrt(42.0);
    const auto f_615 = 3.0 * std::sqrt(42.0);
    const auto f_616 = 0.29296875 * std::sqrt(7.0);
    const auto f_617 = 0.87890625 * std::sqrt(7.0);
    const auto f_618 = 7.03125 * std::sqrt(7.0);
    const auto f_619 = 14.0625 * std::sqrt(7.0);
    const auto f_620 = 3.75 * std::sqrt(7.0);
    const auto f_621 = 0.5859375 * std::sqrt(7.0);
    const auto f_622 = 1.7578125 * std::sqrt(7.0);
    const auto f_623 = 28.125 * std::sqrt(7.0);
    const auto f_624 = 7.5 * std::sqrt(7.0);
    const auto f_625 = 0.78125 * std::sqrt(7.0);
    const auto f_626 = 2.34375 * std::sqrt(7.0);
    const auto f_627 = 18.75 * std::sqrt(7.0);
    const auto f_628 = 37.5 * std::sqrt(7.0);
    const auto f_629 = 10.0 * std::sqrt(7.0);
    const auto f_630 = 0.15625 * std::sqrt(7.0);
    const auto f_631 = 0.46875 * std::sqrt(7.0);
    const auto f_632 = 2.0 * std::sqrt(7.0);
    const auto f_633 = 0.87890625 * std::sqrt(42.0);
    const auto f_634 = 2.8125 * std::sqrt(42.0);
    const auto f_635 = 2.34375 * std::sqrt(42.0);
    const auto f_636 = 12.5 * std::sqrt(42.0);
    const auto f_637 = 7.5 * std::sqrt(42.0);
    const auto f_638 = 0.46875 * std::sqrt(42.0);
    const auto f_639 = 2.5 * std::sqrt(42.0);
    const auto f_640 = 1.5 * std::sqrt(42.0);
    const auto f_641 = 0.3515625 * std::sqrt(231.0);
    const auto f_642 = 1.7578125 * std::sqrt(231.0);
    const auto f_643 = 1.171875 * std::sqrt(231.0);
    const auto f_644 = 2.34375 * std::sqrt(231.0);
    const auto f_645 = 0.9375 * std::sqrt(231.0);
    const auto f_646 = 3.125 * std::sqrt(231.0);
    const auto f_647 = 0.1875 * std::sqrt(231.0);
    const auto f_648 = 0.625 * std::sqrt(231.0);
    const auto f_649 = 0.87890625 * std::sqrt(6006.0);
    const auto f_650 = 0.1171875 * std::sqrt(6006.0);
    const auto f_651 = 1.7578125 * std::sqrt(6006.0);
    const auto f_652 = 0.15625 * std::sqrt(6006.0);
    const auto f_653 = 0.03125 * std::sqrt(6006.0);
    const auto f_654 = 0.46875 * std::sqrt(6006.0);
    const auto f_655 = 0.1640625 * std::sqrt(5005.0);
    const auto f_656 = 0.8203125 * std::sqrt(5005.0);
    const auto f_657 = 0.4921875 * std::sqrt(5005.0);
    const auto f_658 = 0.0234375 * std::sqrt(5005.0);
    const auto f_659 = 0.984375 * std::sqrt(1430.0);
    const auto f_660 = 3.28125 * std::sqrt(1430.0);
    const auto f_661 = 0.8203125 * std::sqrt(55.0);
    const auto f_662 = 1.4765625 * std::sqrt(55.0);
    const auto f_663 = 0.1640625 * std::sqrt(55.0);
    const auto f_664 = 1.4765625 * std::sqrt(5.0);
    const auto f_665 = 2.4609375 * std::sqrt(5.0);
    const auto f_666 = 29.53125 * std::sqrt(5.0);
    const auto f_667 = 0.4921875 * std::sqrt(5.0);
    const auto f_668 = 13.125 * std::sqrt(5.0);
    const auto f_669 = 3.5 * std::sqrt(15.0);
    const auto f_670 = 0.546875 * std::sqrt(105.0);
    const auto f_671 = 1.640625 * std::sqrt(105.0);
    const auto f_672 = 2.625 * std::sqrt(105.0);
    const auto f_673 = 2.4609375 * std::sqrt(10.0);
    const auto f_674 = 7.875 * std::sqrt(10.0);
    const auto f_675 = 0.984375 * std::sqrt(55.0);
    const auto f_676 = 4.921875 * std::sqrt(55.0);
    const auto f_677 = 0.1640625 * std::sqrt(1430.0);
    const auto f_678 = 2.4609375 * std::sqrt(1430.0);
    const auto f_679 = 0.08203125 * std::sqrt(15015.0);
    const auto f_680 = 0.41015625 * std::sqrt(15015.0);
    const auto f_681 = 0.24609375 * std::sqrt(15015.0);
    const auto f_682 = 0.01171875 * std::sqrt(15015.0);
    const auto f_683 = 0.4921875 * std::sqrt(15015.0);
    const auto f_684 = 2.4609375 * std::sqrt(15015.0);
    const auto f_685 = 1.4765625 * std::sqrt(15015.0);
    const auto f_686 = 0.0703125 * std::sqrt(15015.0);
    const auto f_687 = 0.4921875 * std::sqrt(4290.0);
    const auto f_688 = 1.640625 * std::sqrt(4290.0);
    const auto f_689 = 2.953125 * std::sqrt(4290.0);
    const auto f_690 = 9.84375 * std::sqrt(4290.0);
    const auto f_691 = 0.41015625 * std::sqrt(165.0);
    const auto f_692 = 4.921875 * std::sqrt(165.0);
    const auto f_693 = 0.73828125 * std::sqrt(165.0);
    const auto f_694 = 0.08203125 * std::sqrt(165.0);
    const auto f_695 = 0.984375 * std::sqrt(165.0);
    const auto f_696 = 2.4609375 * std::sqrt(165.0);
    const auto f_697 = 29.53125 * std::sqrt(165.0);
    const auto f_698 = 4.4296875 * std::sqrt(165.0);
    const auto f_699 = 59.0625 * std::sqrt(165.0);
    const auto f_700 = 0.4921875 * std::sqrt(165.0);
    const auto f_701 = 5.90625 * std::sqrt(165.0);
    const auto f_702 = 11.8125 * std::sqrt(165.0);
    const auto f_703 = 0.73828125 * std::sqrt(15.0);
    const auto f_704 = 14.765625 * std::sqrt(15.0);
    const auto f_705 = 0.24609375 * std::sqrt(15.0);
    const auto f_706 = 4.4296875 * std::sqrt(15.0);
    const auto f_707 = 7.3828125 * std::sqrt(15.0);
    const auto f_708 = 88.59375 * std::sqrt(15.0);
    const auto f_709 = 1.4765625 * std::sqrt(15.0);
    const auto f_710 = 118.125 * std::sqrt(15.0);
    const auto f_711 = 29.53125 * std::sqrt(15.0);
    const auto f_712 = 7.875 * std::sqrt(30.0);
    const auto f_713 = 14.765625 * std::sqrt(30.0);
    const auto f_714 = 29.53125 * std::sqrt(30.0);
    const auto f_715 = 47.25 * std::sqrt(30.0);
    const auto f_716 = 0.41015625 * std::sqrt(5.0);
    const auto f_717 = 1.23046875 * std::sqrt(5.0);
    const auto f_718 = 5.25 * std::sqrt(5.0);
    const auto f_719 = 7.3828125 * std::sqrt(5.0);
    const auto f_720 = 31.5 * std::sqrt(5.0);
    const auto f_721 = 0.8203125 * std::sqrt(35.0);
    const auto f_722 = 2.4609375 * std::sqrt(35.0);
    const auto f_723 = 4.921875 * std::sqrt(35.0);
    const auto f_724 = 3.9375 * std::sqrt(35.0);
    const auto f_725 = 0.375 * std::sqrt(35.0);
    const auto f_726 = 14.765625 * std::sqrt(35.0);
    const auto f_727 = 29.53125 * std::sqrt(35.0);
    const auto f_728 = 59.0625 * std::sqrt(35.0);
    const auto f_729 = 23.625 * std::sqrt(35.0);
    const auto f_730 = 2.25 * std::sqrt(35.0);
    const auto f_731 = 1.23046875 * std::sqrt(30.0);
    const auto f_732 = 3.9375 * std::sqrt(30.0);
    const auto f_733 = 23.625 * std::sqrt(30.0);
    const auto f_734 = 14.765625 * std::sqrt(165.0);
    const auto f_735 = 0.08203125 * std::sqrt(4290.0);
    const auto f_736 = 1.23046875 * std::sqrt(4290.0);
    const auto f_737 = 7.3828125 * std::sqrt(4290.0);

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
    auto *g_153 = values + 153 * nvalues;
    auto *g_154 = values + 154 * nvalues;
    auto *g_155 = values + 155 * nvalues;
    auto *g_156 = values + 156 * nvalues;
    auto *g_157 = values + 157 * nvalues;
    auto *g_158 = values + 158 * nvalues;
    auto *g_159 = values + 159 * nvalues;
    auto *g_160 = values + 160 * nvalues;
    auto *g_161 = values + 161 * nvalues;
    auto *g_162 = values + 162 * nvalues;
    auto *g_163 = values + 163 * nvalues;
    auto *g_164 = values + 164 * nvalues;

    const auto *hk_0 = buffer.data(hk + 0);
    const auto *hk_1 = buffer.data(hk + 1);
    const auto *hk_2 = buffer.data(hk + 2);
    const auto *hk_3 = buffer.data(hk + 3);
    const auto *hk_4 = buffer.data(hk + 4);
    const auto *hk_5 = buffer.data(hk + 5);
    const auto *hk_6 = buffer.data(hk + 6);
    const auto *hk_7 = buffer.data(hk + 7);
    const auto *hk_8 = buffer.data(hk + 8);
    const auto *hk_9 = buffer.data(hk + 9);
    const auto *hk_10 = buffer.data(hk + 10);
    const auto *hk_11 = buffer.data(hk + 11);
    const auto *hk_12 = buffer.data(hk + 12);
    const auto *hk_13 = buffer.data(hk + 13);
    const auto *hk_14 = buffer.data(hk + 14);
    const auto *hk_15 = buffer.data(hk + 15);
    const auto *hk_16 = buffer.data(hk + 16);
    const auto *hk_17 = buffer.data(hk + 17);
    const auto *hk_18 = buffer.data(hk + 18);
    const auto *hk_19 = buffer.data(hk + 19);
    const auto *hk_20 = buffer.data(hk + 20);
    const auto *hk_21 = buffer.data(hk + 21);
    const auto *hk_22 = buffer.data(hk + 22);
    const auto *hk_23 = buffer.data(hk + 23);
    const auto *hk_24 = buffer.data(hk + 24);
    const auto *hk_25 = buffer.data(hk + 25);
    const auto *hk_26 = buffer.data(hk + 26);
    const auto *hk_27 = buffer.data(hk + 27);
    const auto *hk_28 = buffer.data(hk + 28);
    const auto *hk_29 = buffer.data(hk + 29);
    const auto *hk_30 = buffer.data(hk + 30);
    const auto *hk_31 = buffer.data(hk + 31);
    const auto *hk_32 = buffer.data(hk + 32);
    const auto *hk_33 = buffer.data(hk + 33);
    const auto *hk_34 = buffer.data(hk + 34);
    const auto *hk_35 = buffer.data(hk + 35);
    const auto *hk_36 = buffer.data(hk + 36);
    const auto *hk_37 = buffer.data(hk + 37);
    const auto *hk_38 = buffer.data(hk + 38);
    const auto *hk_39 = buffer.data(hk + 39);
    const auto *hk_40 = buffer.data(hk + 40);
    const auto *hk_41 = buffer.data(hk + 41);
    const auto *hk_42 = buffer.data(hk + 42);
    const auto *hk_43 = buffer.data(hk + 43);
    const auto *hk_44 = buffer.data(hk + 44);
    const auto *hk_45 = buffer.data(hk + 45);
    const auto *hk_46 = buffer.data(hk + 46);
    const auto *hk_47 = buffer.data(hk + 47);
    const auto *hk_48 = buffer.data(hk + 48);
    const auto *hk_49 = buffer.data(hk + 49);
    const auto *hk_50 = buffer.data(hk + 50);
    const auto *hk_51 = buffer.data(hk + 51);
    const auto *hk_52 = buffer.data(hk + 52);
    const auto *hk_53 = buffer.data(hk + 53);
    const auto *hk_54 = buffer.data(hk + 54);
    const auto *hk_55 = buffer.data(hk + 55);
    const auto *hk_56 = buffer.data(hk + 56);
    const auto *hk_57 = buffer.data(hk + 57);
    const auto *hk_58 = buffer.data(hk + 58);
    const auto *hk_59 = buffer.data(hk + 59);
    const auto *hk_60 = buffer.data(hk + 60);
    const auto *hk_61 = buffer.data(hk + 61);
    const auto *hk_62 = buffer.data(hk + 62);
    const auto *hk_63 = buffer.data(hk + 63);
    const auto *hk_64 = buffer.data(hk + 64);
    const auto *hk_65 = buffer.data(hk + 65);
    const auto *hk_66 = buffer.data(hk + 66);
    const auto *hk_67 = buffer.data(hk + 67);
    const auto *hk_68 = buffer.data(hk + 68);
    const auto *hk_69 = buffer.data(hk + 69);
    const auto *hk_70 = buffer.data(hk + 70);
    const auto *hk_71 = buffer.data(hk + 71);
    const auto *hk_72 = buffer.data(hk + 72);
    const auto *hk_73 = buffer.data(hk + 73);
    const auto *hk_74 = buffer.data(hk + 74);
    const auto *hk_75 = buffer.data(hk + 75);
    const auto *hk_76 = buffer.data(hk + 76);
    const auto *hk_77 = buffer.data(hk + 77);
    const auto *hk_78 = buffer.data(hk + 78);
    const auto *hk_79 = buffer.data(hk + 79);
    const auto *hk_80 = buffer.data(hk + 80);
    const auto *hk_81 = buffer.data(hk + 81);
    const auto *hk_82 = buffer.data(hk + 82);
    const auto *hk_83 = buffer.data(hk + 83);
    const auto *hk_84 = buffer.data(hk + 84);
    const auto *hk_85 = buffer.data(hk + 85);
    const auto *hk_86 = buffer.data(hk + 86);
    const auto *hk_87 = buffer.data(hk + 87);
    const auto *hk_88 = buffer.data(hk + 88);
    const auto *hk_89 = buffer.data(hk + 89);
    const auto *hk_90 = buffer.data(hk + 90);
    const auto *hk_91 = buffer.data(hk + 91);
    const auto *hk_92 = buffer.data(hk + 92);
    const auto *hk_93 = buffer.data(hk + 93);
    const auto *hk_94 = buffer.data(hk + 94);
    const auto *hk_95 = buffer.data(hk + 95);
    const auto *hk_96 = buffer.data(hk + 96);
    const auto *hk_97 = buffer.data(hk + 97);
    const auto *hk_98 = buffer.data(hk + 98);
    const auto *hk_99 = buffer.data(hk + 99);
    const auto *hk_100 = buffer.data(hk + 100);
    const auto *hk_101 = buffer.data(hk + 101);
    const auto *hk_102 = buffer.data(hk + 102);
    const auto *hk_103 = buffer.data(hk + 103);
    const auto *hk_104 = buffer.data(hk + 104);
    const auto *hk_105 = buffer.data(hk + 105);
    const auto *hk_106 = buffer.data(hk + 106);
    const auto *hk_107 = buffer.data(hk + 107);
    const auto *hk_108 = buffer.data(hk + 108);
    const auto *hk_109 = buffer.data(hk + 109);
    const auto *hk_110 = buffer.data(hk + 110);
    const auto *hk_111 = buffer.data(hk + 111);
    const auto *hk_112 = buffer.data(hk + 112);
    const auto *hk_113 = buffer.data(hk + 113);
    const auto *hk_114 = buffer.data(hk + 114);
    const auto *hk_115 = buffer.data(hk + 115);
    const auto *hk_116 = buffer.data(hk + 116);
    const auto *hk_117 = buffer.data(hk + 117);
    const auto *hk_118 = buffer.data(hk + 118);
    const auto *hk_119 = buffer.data(hk + 119);
    const auto *hk_120 = buffer.data(hk + 120);
    const auto *hk_121 = buffer.data(hk + 121);
    const auto *hk_122 = buffer.data(hk + 122);
    const auto *hk_123 = buffer.data(hk + 123);
    const auto *hk_124 = buffer.data(hk + 124);
    const auto *hk_125 = buffer.data(hk + 125);
    const auto *hk_126 = buffer.data(hk + 126);
    const auto *hk_127 = buffer.data(hk + 127);
    const auto *hk_128 = buffer.data(hk + 128);
    const auto *hk_129 = buffer.data(hk + 129);
    const auto *hk_130 = buffer.data(hk + 130);
    const auto *hk_131 = buffer.data(hk + 131);
    const auto *hk_132 = buffer.data(hk + 132);
    const auto *hk_133 = buffer.data(hk + 133);
    const auto *hk_134 = buffer.data(hk + 134);
    const auto *hk_135 = buffer.data(hk + 135);
    const auto *hk_136 = buffer.data(hk + 136);
    const auto *hk_137 = buffer.data(hk + 137);
    const auto *hk_138 = buffer.data(hk + 138);
    const auto *hk_139 = buffer.data(hk + 139);
    const auto *hk_140 = buffer.data(hk + 140);
    const auto *hk_141 = buffer.data(hk + 141);
    const auto *hk_142 = buffer.data(hk + 142);
    const auto *hk_143 = buffer.data(hk + 143);
    const auto *hk_144 = buffer.data(hk + 144);
    const auto *hk_145 = buffer.data(hk + 145);
    const auto *hk_146 = buffer.data(hk + 146);
    const auto *hk_147 = buffer.data(hk + 147);
    const auto *hk_148 = buffer.data(hk + 148);
    const auto *hk_149 = buffer.data(hk + 149);
    const auto *hk_150 = buffer.data(hk + 150);
    const auto *hk_151 = buffer.data(hk + 151);
    const auto *hk_152 = buffer.data(hk + 152);
    const auto *hk_153 = buffer.data(hk + 153);
    const auto *hk_154 = buffer.data(hk + 154);
    const auto *hk_155 = buffer.data(hk + 155);
    const auto *hk_156 = buffer.data(hk + 156);
    const auto *hk_157 = buffer.data(hk + 157);
    const auto *hk_158 = buffer.data(hk + 158);
    const auto *hk_159 = buffer.data(hk + 159);
    const auto *hk_160 = buffer.data(hk + 160);
    const auto *hk_161 = buffer.data(hk + 161);
    const auto *hk_162 = buffer.data(hk + 162);
    const auto *hk_163 = buffer.data(hk + 163);
    const auto *hk_164 = buffer.data(hk + 164);
    const auto *hk_165 = buffer.data(hk + 165);
    const auto *hk_166 = buffer.data(hk + 166);
    const auto *hk_167 = buffer.data(hk + 167);
    const auto *hk_168 = buffer.data(hk + 168);
    const auto *hk_169 = buffer.data(hk + 169);
    const auto *hk_170 = buffer.data(hk + 170);
    const auto *hk_171 = buffer.data(hk + 171);
    const auto *hk_172 = buffer.data(hk + 172);
    const auto *hk_173 = buffer.data(hk + 173);
    const auto *hk_174 = buffer.data(hk + 174);
    const auto *hk_175 = buffer.data(hk + 175);
    const auto *hk_176 = buffer.data(hk + 176);
    const auto *hk_177 = buffer.data(hk + 177);
    const auto *hk_178 = buffer.data(hk + 178);
    const auto *hk_179 = buffer.data(hk + 179);
    const auto *hk_180 = buffer.data(hk + 180);
    const auto *hk_181 = buffer.data(hk + 181);
    const auto *hk_182 = buffer.data(hk + 182);
    const auto *hk_183 = buffer.data(hk + 183);
    const auto *hk_184 = buffer.data(hk + 184);
    const auto *hk_185 = buffer.data(hk + 185);
    const auto *hk_186 = buffer.data(hk + 186);
    const auto *hk_187 = buffer.data(hk + 187);
    const auto *hk_188 = buffer.data(hk + 188);
    const auto *hk_189 = buffer.data(hk + 189);
    const auto *hk_190 = buffer.data(hk + 190);
    const auto *hk_191 = buffer.data(hk + 191);
    const auto *hk_192 = buffer.data(hk + 192);
    const auto *hk_193 = buffer.data(hk + 193);
    const auto *hk_194 = buffer.data(hk + 194);
    const auto *hk_195 = buffer.data(hk + 195);
    const auto *hk_196 = buffer.data(hk + 196);
    const auto *hk_197 = buffer.data(hk + 197);
    const auto *hk_198 = buffer.data(hk + 198);
    const auto *hk_199 = buffer.data(hk + 199);
    const auto *hk_200 = buffer.data(hk + 200);
    const auto *hk_201 = buffer.data(hk + 201);
    const auto *hk_202 = buffer.data(hk + 202);
    const auto *hk_203 = buffer.data(hk + 203);
    const auto *hk_204 = buffer.data(hk + 204);
    const auto *hk_205 = buffer.data(hk + 205);
    const auto *hk_206 = buffer.data(hk + 206);
    const auto *hk_207 = buffer.data(hk + 207);
    const auto *hk_208 = buffer.data(hk + 208);
    const auto *hk_209 = buffer.data(hk + 209);
    const auto *hk_210 = buffer.data(hk + 210);
    const auto *hk_211 = buffer.data(hk + 211);
    const auto *hk_212 = buffer.data(hk + 212);
    const auto *hk_213 = buffer.data(hk + 213);
    const auto *hk_214 = buffer.data(hk + 214);
    const auto *hk_215 = buffer.data(hk + 215);
    const auto *hk_216 = buffer.data(hk + 216);
    const auto *hk_217 = buffer.data(hk + 217);
    const auto *hk_218 = buffer.data(hk + 218);
    const auto *hk_219 = buffer.data(hk + 219);
    const auto *hk_220 = buffer.data(hk + 220);
    const auto *hk_221 = buffer.data(hk + 221);
    const auto *hk_222 = buffer.data(hk + 222);
    const auto *hk_223 = buffer.data(hk + 223);
    const auto *hk_224 = buffer.data(hk + 224);
    const auto *hk_225 = buffer.data(hk + 225);
    const auto *hk_226 = buffer.data(hk + 226);
    const auto *hk_227 = buffer.data(hk + 227);
    const auto *hk_228 = buffer.data(hk + 228);
    const auto *hk_229 = buffer.data(hk + 229);
    const auto *hk_230 = buffer.data(hk + 230);
    const auto *hk_231 = buffer.data(hk + 231);
    const auto *hk_232 = buffer.data(hk + 232);
    const auto *hk_233 = buffer.data(hk + 233);
    const auto *hk_234 = buffer.data(hk + 234);
    const auto *hk_235 = buffer.data(hk + 235);
    const auto *hk_236 = buffer.data(hk + 236);
    const auto *hk_237 = buffer.data(hk + 237);
    const auto *hk_238 = buffer.data(hk + 238);
    const auto *hk_239 = buffer.data(hk + 239);
    const auto *hk_240 = buffer.data(hk + 240);
    const auto *hk_241 = buffer.data(hk + 241);
    const auto *hk_242 = buffer.data(hk + 242);
    const auto *hk_243 = buffer.data(hk + 243);
    const auto *hk_244 = buffer.data(hk + 244);
    const auto *hk_245 = buffer.data(hk + 245);
    const auto *hk_246 = buffer.data(hk + 246);
    const auto *hk_247 = buffer.data(hk + 247);
    const auto *hk_248 = buffer.data(hk + 248);
    const auto *hk_249 = buffer.data(hk + 249);
    const auto *hk_250 = buffer.data(hk + 250);
    const auto *hk_251 = buffer.data(hk + 251);
    const auto *hk_252 = buffer.data(hk + 252);
    const auto *hk_253 = buffer.data(hk + 253);
    const auto *hk_254 = buffer.data(hk + 254);
    const auto *hk_255 = buffer.data(hk + 255);
    const auto *hk_256 = buffer.data(hk + 256);
    const auto *hk_257 = buffer.data(hk + 257);
    const auto *hk_258 = buffer.data(hk + 258);
    const auto *hk_259 = buffer.data(hk + 259);
    const auto *hk_260 = buffer.data(hk + 260);
    const auto *hk_261 = buffer.data(hk + 261);
    const auto *hk_262 = buffer.data(hk + 262);
    const auto *hk_263 = buffer.data(hk + 263);
    const auto *hk_264 = buffer.data(hk + 264);
    const auto *hk_265 = buffer.data(hk + 265);
    const auto *hk_266 = buffer.data(hk + 266);
    const auto *hk_267 = buffer.data(hk + 267);
    const auto *hk_268 = buffer.data(hk + 268);
    const auto *hk_269 = buffer.data(hk + 269);
    const auto *hk_270 = buffer.data(hk + 270);
    const auto *hk_271 = buffer.data(hk + 271);
    const auto *hk_272 = buffer.data(hk + 272);
    const auto *hk_273 = buffer.data(hk + 273);
    const auto *hk_274 = buffer.data(hk + 274);
    const auto *hk_275 = buffer.data(hk + 275);
    const auto *hk_276 = buffer.data(hk + 276);
    const auto *hk_277 = buffer.data(hk + 277);
    const auto *hk_278 = buffer.data(hk + 278);
    const auto *hk_279 = buffer.data(hk + 279);
    const auto *hk_280 = buffer.data(hk + 280);
    const auto *hk_281 = buffer.data(hk + 281);
    const auto *hk_282 = buffer.data(hk + 282);
    const auto *hk_283 = buffer.data(hk + 283);
    const auto *hk_284 = buffer.data(hk + 284);
    const auto *hk_285 = buffer.data(hk + 285);
    const auto *hk_286 = buffer.data(hk + 286);
    const auto *hk_287 = buffer.data(hk + 287);
    const auto *hk_288 = buffer.data(hk + 288);
    const auto *hk_289 = buffer.data(hk + 289);
    const auto *hk_290 = buffer.data(hk + 290);
    const auto *hk_291 = buffer.data(hk + 291);
    const auto *hk_292 = buffer.data(hk + 292);
    const auto *hk_293 = buffer.data(hk + 293);
    const auto *hk_294 = buffer.data(hk + 294);
    const auto *hk_295 = buffer.data(hk + 295);
    const auto *hk_296 = buffer.data(hk + 296);
    const auto *hk_297 = buffer.data(hk + 297);
    const auto *hk_298 = buffer.data(hk + 298);
    const auto *hk_299 = buffer.data(hk + 299);
    const auto *hk_300 = buffer.data(hk + 300);
    const auto *hk_301 = buffer.data(hk + 301);
    const auto *hk_302 = buffer.data(hk + 302);
    const auto *hk_303 = buffer.data(hk + 303);
    const auto *hk_304 = buffer.data(hk + 304);
    const auto *hk_305 = buffer.data(hk + 305);
    const auto *hk_306 = buffer.data(hk + 306);
    const auto *hk_307 = buffer.data(hk + 307);
    const auto *hk_308 = buffer.data(hk + 308);
    const auto *hk_309 = buffer.data(hk + 309);
    const auto *hk_310 = buffer.data(hk + 310);
    const auto *hk_311 = buffer.data(hk + 311);
    const auto *hk_312 = buffer.data(hk + 312);
    const auto *hk_313 = buffer.data(hk + 313);
    const auto *hk_314 = buffer.data(hk + 314);
    const auto *hk_315 = buffer.data(hk + 315);
    const auto *hk_316 = buffer.data(hk + 316);
    const auto *hk_317 = buffer.data(hk + 317);
    const auto *hk_318 = buffer.data(hk + 318);
    const auto *hk_319 = buffer.data(hk + 319);
    const auto *hk_320 = buffer.data(hk + 320);
    const auto *hk_321 = buffer.data(hk + 321);
    const auto *hk_322 = buffer.data(hk + 322);
    const auto *hk_323 = buffer.data(hk + 323);
    const auto *hk_324 = buffer.data(hk + 324);
    const auto *hk_325 = buffer.data(hk + 325);
    const auto *hk_326 = buffer.data(hk + 326);
    const auto *hk_327 = buffer.data(hk + 327);
    const auto *hk_328 = buffer.data(hk + 328);
    const auto *hk_329 = buffer.data(hk + 329);
    const auto *hk_330 = buffer.data(hk + 330);
    const auto *hk_331 = buffer.data(hk + 331);
    const auto *hk_332 = buffer.data(hk + 332);
    const auto *hk_333 = buffer.data(hk + 333);
    const auto *hk_334 = buffer.data(hk + 334);
    const auto *hk_335 = buffer.data(hk + 335);
    const auto *hk_336 = buffer.data(hk + 336);
    const auto *hk_337 = buffer.data(hk + 337);
    const auto *hk_338 = buffer.data(hk + 338);
    const auto *hk_339 = buffer.data(hk + 339);
    const auto *hk_340 = buffer.data(hk + 340);
    const auto *hk_341 = buffer.data(hk + 341);
    const auto *hk_342 = buffer.data(hk + 342);
    const auto *hk_343 = buffer.data(hk + 343);
    const auto *hk_344 = buffer.data(hk + 344);
    const auto *hk_345 = buffer.data(hk + 345);
    const auto *hk_346 = buffer.data(hk + 346);
    const auto *hk_347 = buffer.data(hk + 347);
    const auto *hk_348 = buffer.data(hk + 348);
    const auto *hk_349 = buffer.data(hk + 349);
    const auto *hk_350 = buffer.data(hk + 350);
    const auto *hk_351 = buffer.data(hk + 351);
    const auto *hk_352 = buffer.data(hk + 352);
    const auto *hk_353 = buffer.data(hk + 353);
    const auto *hk_354 = buffer.data(hk + 354);
    const auto *hk_355 = buffer.data(hk + 355);
    const auto *hk_356 = buffer.data(hk + 356);
    const auto *hk_357 = buffer.data(hk + 357);
    const auto *hk_358 = buffer.data(hk + 358);
    const auto *hk_359 = buffer.data(hk + 359);
    const auto *hk_360 = buffer.data(hk + 360);
    const auto *hk_361 = buffer.data(hk + 361);
    const auto *hk_362 = buffer.data(hk + 362);
    const auto *hk_363 = buffer.data(hk + 363);
    const auto *hk_364 = buffer.data(hk + 364);
    const auto *hk_365 = buffer.data(hk + 365);
    const auto *hk_366 = buffer.data(hk + 366);
    const auto *hk_367 = buffer.data(hk + 367);
    const auto *hk_368 = buffer.data(hk + 368);
    const auto *hk_369 = buffer.data(hk + 369);
    const auto *hk_370 = buffer.data(hk + 370);
    const auto *hk_371 = buffer.data(hk + 371);
    const auto *hk_372 = buffer.data(hk + 372);
    const auto *hk_373 = buffer.data(hk + 373);
    const auto *hk_374 = buffer.data(hk + 374);
    const auto *hk_375 = buffer.data(hk + 375);
    const auto *hk_376 = buffer.data(hk + 376);
    const auto *hk_377 = buffer.data(hk + 377);
    const auto *hk_378 = buffer.data(hk + 378);
    const auto *hk_379 = buffer.data(hk + 379);
    const auto *hk_380 = buffer.data(hk + 380);
    const auto *hk_381 = buffer.data(hk + 381);
    const auto *hk_382 = buffer.data(hk + 382);
    const auto *hk_383 = buffer.data(hk + 383);
    const auto *hk_384 = buffer.data(hk + 384);
    const auto *hk_385 = buffer.data(hk + 385);
    const auto *hk_386 = buffer.data(hk + 386);
    const auto *hk_387 = buffer.data(hk + 387);
    const auto *hk_388 = buffer.data(hk + 388);
    const auto *hk_389 = buffer.data(hk + 389);
    const auto *hk_390 = buffer.data(hk + 390);
    const auto *hk_391 = buffer.data(hk + 391);
    const auto *hk_392 = buffer.data(hk + 392);
    const auto *hk_393 = buffer.data(hk + 393);
    const auto *hk_394 = buffer.data(hk + 394);
    const auto *hk_395 = buffer.data(hk + 395);
    const auto *hk_396 = buffer.data(hk + 396);
    const auto *hk_397 = buffer.data(hk + 397);
    const auto *hk_398 = buffer.data(hk + 398);
    const auto *hk_399 = buffer.data(hk + 399);
    const auto *hk_400 = buffer.data(hk + 400);
    const auto *hk_401 = buffer.data(hk + 401);
    const auto *hk_402 = buffer.data(hk + 402);
    const auto *hk_403 = buffer.data(hk + 403);
    const auto *hk_404 = buffer.data(hk + 404);
    const auto *hk_405 = buffer.data(hk + 405);
    const auto *hk_406 = buffer.data(hk + 406);
    const auto *hk_407 = buffer.data(hk + 407);
    const auto *hk_408 = buffer.data(hk + 408);
    const auto *hk_409 = buffer.data(hk + 409);
    const auto *hk_410 = buffer.data(hk + 410);
    const auto *hk_411 = buffer.data(hk + 411);
    const auto *hk_412 = buffer.data(hk + 412);
    const auto *hk_413 = buffer.data(hk + 413);
    const auto *hk_414 = buffer.data(hk + 414);
    const auto *hk_415 = buffer.data(hk + 415);
    const auto *hk_416 = buffer.data(hk + 416);
    const auto *hk_417 = buffer.data(hk + 417);
    const auto *hk_418 = buffer.data(hk + 418);
    const auto *hk_419 = buffer.data(hk + 419);
    const auto *hk_420 = buffer.data(hk + 420);
    const auto *hk_421 = buffer.data(hk + 421);
    const auto *hk_422 = buffer.data(hk + 422);
    const auto *hk_423 = buffer.data(hk + 423);
    const auto *hk_424 = buffer.data(hk + 424);
    const auto *hk_425 = buffer.data(hk + 425);
    const auto *hk_426 = buffer.data(hk + 426);
    const auto *hk_427 = buffer.data(hk + 427);
    const auto *hk_428 = buffer.data(hk + 428);
    const auto *hk_429 = buffer.data(hk + 429);
    const auto *hk_430 = buffer.data(hk + 430);
    const auto *hk_431 = buffer.data(hk + 431);
    const auto *hk_432 = buffer.data(hk + 432);
    const auto *hk_433 = buffer.data(hk + 433);
    const auto *hk_434 = buffer.data(hk + 434);
    const auto *hk_435 = buffer.data(hk + 435);
    const auto *hk_436 = buffer.data(hk + 436);
    const auto *hk_437 = buffer.data(hk + 437);
    const auto *hk_438 = buffer.data(hk + 438);
    const auto *hk_439 = buffer.data(hk + 439);
    const auto *hk_440 = buffer.data(hk + 440);
    const auto *hk_441 = buffer.data(hk + 441);
    const auto *hk_442 = buffer.data(hk + 442);
    const auto *hk_443 = buffer.data(hk + 443);
    const auto *hk_444 = buffer.data(hk + 444);
    const auto *hk_445 = buffer.data(hk + 445);
    const auto *hk_446 = buffer.data(hk + 446);
    const auto *hk_447 = buffer.data(hk + 447);
    const auto *hk_448 = buffer.data(hk + 448);
    const auto *hk_449 = buffer.data(hk + 449);
    const auto *hk_450 = buffer.data(hk + 450);
    const auto *hk_451 = buffer.data(hk + 451);
    const auto *hk_452 = buffer.data(hk + 452);
    const auto *hk_453 = buffer.data(hk + 453);
    const auto *hk_454 = buffer.data(hk + 454);
    const auto *hk_455 = buffer.data(hk + 455);
    const auto *hk_456 = buffer.data(hk + 456);
    const auto *hk_457 = buffer.data(hk + 457);
    const auto *hk_458 = buffer.data(hk + 458);
    const auto *hk_459 = buffer.data(hk + 459);
    const auto *hk_460 = buffer.data(hk + 460);
    const auto *hk_461 = buffer.data(hk + 461);
    const auto *hk_462 = buffer.data(hk + 462);
    const auto *hk_463 = buffer.data(hk + 463);
    const auto *hk_464 = buffer.data(hk + 464);
    const auto *hk_465 = buffer.data(hk + 465);
    const auto *hk_466 = buffer.data(hk + 466);
    const auto *hk_467 = buffer.data(hk + 467);
    const auto *hk_468 = buffer.data(hk + 468);
    const auto *hk_469 = buffer.data(hk + 469);
    const auto *hk_470 = buffer.data(hk + 470);
    const auto *hk_471 = buffer.data(hk + 471);
    const auto *hk_472 = buffer.data(hk + 472);
    const auto *hk_473 = buffer.data(hk + 473);
    const auto *hk_474 = buffer.data(hk + 474);
    const auto *hk_475 = buffer.data(hk + 475);
    const auto *hk_476 = buffer.data(hk + 476);
    const auto *hk_477 = buffer.data(hk + 477);
    const auto *hk_478 = buffer.data(hk + 478);
    const auto *hk_479 = buffer.data(hk + 479);
    const auto *hk_480 = buffer.data(hk + 480);
    const auto *hk_481 = buffer.data(hk + 481);
    const auto *hk_482 = buffer.data(hk + 482);
    const auto *hk_483 = buffer.data(hk + 483);
    const auto *hk_484 = buffer.data(hk + 484);
    const auto *hk_485 = buffer.data(hk + 485);
    const auto *hk_486 = buffer.data(hk + 486);
    const auto *hk_487 = buffer.data(hk + 487);
    const auto *hk_488 = buffer.data(hk + 488);
    const auto *hk_489 = buffer.data(hk + 489);
    const auto *hk_490 = buffer.data(hk + 490);
    const auto *hk_491 = buffer.data(hk + 491);
    const auto *hk_492 = buffer.data(hk + 492);
    const auto *hk_493 = buffer.data(hk + 493);
    const auto *hk_494 = buffer.data(hk + 494);
    const auto *hk_495 = buffer.data(hk + 495);
    const auto *hk_496 = buffer.data(hk + 496);
    const auto *hk_497 = buffer.data(hk + 497);
    const auto *hk_498 = buffer.data(hk + 498);
    const auto *hk_499 = buffer.data(hk + 499);
    const auto *hk_500 = buffer.data(hk + 500);
    const auto *hk_501 = buffer.data(hk + 501);
    const auto *hk_502 = buffer.data(hk + 502);
    const auto *hk_503 = buffer.data(hk + 503);
    const auto *hk_504 = buffer.data(hk + 504);
    const auto *hk_505 = buffer.data(hk + 505);
    const auto *hk_506 = buffer.data(hk + 506);
    const auto *hk_507 = buffer.data(hk + 507);
    const auto *hk_508 = buffer.data(hk + 508);
    const auto *hk_509 = buffer.data(hk + 509);
    const auto *hk_510 = buffer.data(hk + 510);
    const auto *hk_511 = buffer.data(hk + 511);
    const auto *hk_512 = buffer.data(hk + 512);
    const auto *hk_513 = buffer.data(hk + 513);
    const auto *hk_514 = buffer.data(hk + 514);
    const auto *hk_515 = buffer.data(hk + 515);
    const auto *hk_516 = buffer.data(hk + 516);
    const auto *hk_517 = buffer.data(hk + 517);
    const auto *hk_518 = buffer.data(hk + 518);
    const auto *hk_519 = buffer.data(hk + 519);
    const auto *hk_520 = buffer.data(hk + 520);
    const auto *hk_521 = buffer.data(hk + 521);
    const auto *hk_522 = buffer.data(hk + 522);
    const auto *hk_523 = buffer.data(hk + 523);
    const auto *hk_524 = buffer.data(hk + 524);
    const auto *hk_525 = buffer.data(hk + 525);
    const auto *hk_526 = buffer.data(hk + 526);
    const auto *hk_527 = buffer.data(hk + 527);
    const auto *hk_528 = buffer.data(hk + 528);
    const auto *hk_529 = buffer.data(hk + 529);
    const auto *hk_530 = buffer.data(hk + 530);
    const auto *hk_531 = buffer.data(hk + 531);
    const auto *hk_532 = buffer.data(hk + 532);
    const auto *hk_533 = buffer.data(hk + 533);
    const auto *hk_534 = buffer.data(hk + 534);
    const auto *hk_535 = buffer.data(hk + 535);
    const auto *hk_536 = buffer.data(hk + 536);
    const auto *hk_537 = buffer.data(hk + 537);
    const auto *hk_538 = buffer.data(hk + 538);
    const auto *hk_539 = buffer.data(hk + 539);
    const auto *hk_540 = buffer.data(hk + 540);
    const auto *hk_541 = buffer.data(hk + 541);
    const auto *hk_542 = buffer.data(hk + 542);
    const auto *hk_543 = buffer.data(hk + 543);
    const auto *hk_544 = buffer.data(hk + 544);
    const auto *hk_545 = buffer.data(hk + 545);
    const auto *hk_546 = buffer.data(hk + 546);
    const auto *hk_547 = buffer.data(hk + 547);
    const auto *hk_548 = buffer.data(hk + 548);
    const auto *hk_549 = buffer.data(hk + 549);
    const auto *hk_550 = buffer.data(hk + 550);
    const auto *hk_551 = buffer.data(hk + 551);
    const auto *hk_552 = buffer.data(hk + 552);
    const auto *hk_553 = buffer.data(hk + 553);
    const auto *hk_554 = buffer.data(hk + 554);
    const auto *hk_555 = buffer.data(hk + 555);
    const auto *hk_556 = buffer.data(hk + 556);
    const auto *hk_557 = buffer.data(hk + 557);
    const auto *hk_558 = buffer.data(hk + 558);
    const auto *hk_559 = buffer.data(hk + 559);
    const auto *hk_560 = buffer.data(hk + 560);
    const auto *hk_561 = buffer.data(hk + 561);
    const auto *hk_562 = buffer.data(hk + 562);
    const auto *hk_563 = buffer.data(hk + 563);
    const auto *hk_564 = buffer.data(hk + 564);
    const auto *hk_565 = buffer.data(hk + 565);
    const auto *hk_566 = buffer.data(hk + 566);
    const auto *hk_567 = buffer.data(hk + 567);
    const auto *hk_568 = buffer.data(hk + 568);
    const auto *hk_569 = buffer.data(hk + 569);
    const auto *hk_570 = buffer.data(hk + 570);
    const auto *hk_571 = buffer.data(hk + 571);
    const auto *hk_572 = buffer.data(hk + 572);
    const auto *hk_573 = buffer.data(hk + 573);
    const auto *hk_574 = buffer.data(hk + 574);
    const auto *hk_575 = buffer.data(hk + 575);
    const auto *hk_576 = buffer.data(hk + 576);
    const auto *hk_577 = buffer.data(hk + 577);
    const auto *hk_578 = buffer.data(hk + 578);
    const auto *hk_579 = buffer.data(hk + 579);
    const auto *hk_580 = buffer.data(hk + 580);
    const auto *hk_581 = buffer.data(hk + 581);
    const auto *hk_582 = buffer.data(hk + 582);
    const auto *hk_583 = buffer.data(hk + 583);
    const auto *hk_584 = buffer.data(hk + 584);
    const auto *hk_585 = buffer.data(hk + 585);
    const auto *hk_586 = buffer.data(hk + 586);
    const auto *hk_587 = buffer.data(hk + 587);
    const auto *hk_588 = buffer.data(hk + 588);
    const auto *hk_589 = buffer.data(hk + 589);
    const auto *hk_590 = buffer.data(hk + 590);
    const auto *hk_591 = buffer.data(hk + 591);
    const auto *hk_592 = buffer.data(hk + 592);
    const auto *hk_593 = buffer.data(hk + 593);
    const auto *hk_594 = buffer.data(hk + 594);
    const auto *hk_595 = buffer.data(hk + 595);
    const auto *hk_596 = buffer.data(hk + 596);
    const auto *hk_597 = buffer.data(hk + 597);
    const auto *hk_598 = buffer.data(hk + 598);
    const auto *hk_599 = buffer.data(hk + 599);
    const auto *hk_600 = buffer.data(hk + 600);
    const auto *hk_601 = buffer.data(hk + 601);
    const auto *hk_602 = buffer.data(hk + 602);
    const auto *hk_603 = buffer.data(hk + 603);
    const auto *hk_604 = buffer.data(hk + 604);
    const auto *hk_605 = buffer.data(hk + 605);
    const auto *hk_606 = buffer.data(hk + 606);
    const auto *hk_607 = buffer.data(hk + 607);
    const auto *hk_608 = buffer.data(hk + 608);
    const auto *hk_609 = buffer.data(hk + 609);
    const auto *hk_610 = buffer.data(hk + 610);
    const auto *hk_611 = buffer.data(hk + 611);
    const auto *hk_612 = buffer.data(hk + 612);
    const auto *hk_613 = buffer.data(hk + 613);
    const auto *hk_614 = buffer.data(hk + 614);
    const auto *hk_615 = buffer.data(hk + 615);
    const auto *hk_616 = buffer.data(hk + 616);
    const auto *hk_617 = buffer.data(hk + 617);
    const auto *hk_618 = buffer.data(hk + 618);
    const auto *hk_619 = buffer.data(hk + 619);
    const auto *hk_620 = buffer.data(hk + 620);
    const auto *hk_621 = buffer.data(hk + 621);
    const auto *hk_622 = buffer.data(hk + 622);
    const auto *hk_623 = buffer.data(hk + 623);
    const auto *hk_624 = buffer.data(hk + 624);
    const auto *hk_625 = buffer.data(hk + 625);
    const auto *hk_626 = buffer.data(hk + 626);
    const auto *hk_627 = buffer.data(hk + 627);
    const auto *hk_628 = buffer.data(hk + 628);
    const auto *hk_629 = buffer.data(hk + 629);
    const auto *hk_630 = buffer.data(hk + 630);
    const auto *hk_631 = buffer.data(hk + 631);
    const auto *hk_632 = buffer.data(hk + 632);
    const auto *hk_633 = buffer.data(hk + 633);
    const auto *hk_634 = buffer.data(hk + 634);
    const auto *hk_635 = buffer.data(hk + 635);
    const auto *hk_636 = buffer.data(hk + 636);
    const auto *hk_637 = buffer.data(hk + 637);
    const auto *hk_638 = buffer.data(hk + 638);
    const auto *hk_639 = buffer.data(hk + 639);
    const auto *hk_640 = buffer.data(hk + 640);
    const auto *hk_641 = buffer.data(hk + 641);
    const auto *hk_642 = buffer.data(hk + 642);
    const auto *hk_643 = buffer.data(hk + 643);
    const auto *hk_644 = buffer.data(hk + 644);
    const auto *hk_645 = buffer.data(hk + 645);
    const auto *hk_646 = buffer.data(hk + 646);
    const auto *hk_647 = buffer.data(hk + 647);
    const auto *hk_648 = buffer.data(hk + 648);
    const auto *hk_649 = buffer.data(hk + 649);
    const auto *hk_650 = buffer.data(hk + 650);
    const auto *hk_651 = buffer.data(hk + 651);
    const auto *hk_652 = buffer.data(hk + 652);
    const auto *hk_653 = buffer.data(hk + 653);
    const auto *hk_654 = buffer.data(hk + 654);
    const auto *hk_655 = buffer.data(hk + 655);
    const auto *hk_656 = buffer.data(hk + 656);
    const auto *hk_657 = buffer.data(hk + 657);
    const auto *hk_658 = buffer.data(hk + 658);
    const auto *hk_659 = buffer.data(hk + 659);
    const auto *hk_660 = buffer.data(hk + 660);
    const auto *hk_661 = buffer.data(hk + 661);
    const auto *hk_662 = buffer.data(hk + 662);
    const auto *hk_663 = buffer.data(hk + 663);
    const auto *hk_664 = buffer.data(hk + 664);
    const auto *hk_665 = buffer.data(hk + 665);
    const auto *hk_666 = buffer.data(hk + 666);
    const auto *hk_667 = buffer.data(hk + 667);
    const auto *hk_668 = buffer.data(hk + 668);
    const auto *hk_669 = buffer.data(hk + 669);
    const auto *hk_670 = buffer.data(hk + 670);
    const auto *hk_671 = buffer.data(hk + 671);
    const auto *hk_672 = buffer.data(hk + 672);
    const auto *hk_673 = buffer.data(hk + 673);
    const auto *hk_674 = buffer.data(hk + 674);
    const auto *hk_675 = buffer.data(hk + 675);
    const auto *hk_676 = buffer.data(hk + 676);
    const auto *hk_677 = buffer.data(hk + 677);
    const auto *hk_678 = buffer.data(hk + 678);
    const auto *hk_679 = buffer.data(hk + 679);
    const auto *hk_680 = buffer.data(hk + 680);
    const auto *hk_681 = buffer.data(hk + 681);
    const auto *hk_682 = buffer.data(hk + 682);
    const auto *hk_683 = buffer.data(hk + 683);
    const auto *hk_684 = buffer.data(hk + 684);
    const auto *hk_685 = buffer.data(hk + 685);
    const auto *hk_686 = buffer.data(hk + 686);
    const auto *hk_687 = buffer.data(hk + 687);
    const auto *hk_688 = buffer.data(hk + 688);
    const auto *hk_689 = buffer.data(hk + 689);
    const auto *hk_690 = buffer.data(hk + 690);
    const auto *hk_691 = buffer.data(hk + 691);
    const auto *hk_692 = buffer.data(hk + 692);
    const auto *hk_693 = buffer.data(hk + 693);
    const auto *hk_694 = buffer.data(hk + 694);
    const auto *hk_695 = buffer.data(hk + 695);
    const auto *hk_696 = buffer.data(hk + 696);
    const auto *hk_697 = buffer.data(hk + 697);
    const auto *hk_698 = buffer.data(hk + 698);
    const auto *hk_699 = buffer.data(hk + 699);
    const auto *hk_700 = buffer.data(hk + 700);
    const auto *hk_701 = buffer.data(hk + 701);
    const auto *hk_702 = buffer.data(hk + 702);
    const auto *hk_703 = buffer.data(hk + 703);
    const auto *hk_704 = buffer.data(hk + 704);
    const auto *hk_705 = buffer.data(hk + 705);
    const auto *hk_706 = buffer.data(hk + 706);
    const auto *hk_707 = buffer.data(hk + 707);
    const auto *hk_708 = buffer.data(hk + 708);
    const auto *hk_709 = buffer.data(hk + 709);
    const auto *hk_710 = buffer.data(hk + 710);
    const auto *hk_711 = buffer.data(hk + 711);
    const auto *hk_712 = buffer.data(hk + 712);
    const auto *hk_713 = buffer.data(hk + 713);
    const auto *hk_714 = buffer.data(hk + 714);
    const auto *hk_715 = buffer.data(hk + 715);
    const auto *hk_716 = buffer.data(hk + 716);
    const auto *hk_717 = buffer.data(hk + 717);
    const auto *hk_718 = buffer.data(hk + 718);
    const auto *hk_719 = buffer.data(hk + 719);
    const auto *hk_720 = buffer.data(hk + 720);
    const auto *hk_721 = buffer.data(hk + 721);
    const auto *hk_722 = buffer.data(hk + 722);
    const auto *hk_723 = buffer.data(hk + 723);
    const auto *hk_724 = buffer.data(hk + 724);
    const auto *hk_725 = buffer.data(hk + 725);
    const auto *hk_726 = buffer.data(hk + 726);
    const auto *hk_727 = buffer.data(hk + 727);
    const auto *hk_728 = buffer.data(hk + 728);
    const auto *hk_729 = buffer.data(hk + 729);
    const auto *hk_730 = buffer.data(hk + 730);
    const auto *hk_731 = buffer.data(hk + 731);
    const auto *hk_732 = buffer.data(hk + 732);
    const auto *hk_733 = buffer.data(hk + 733);
    const auto *hk_734 = buffer.data(hk + 734);
    const auto *hk_735 = buffer.data(hk + 735);
    const auto *hk_736 = buffer.data(hk + 736);
    const auto *hk_737 = buffer.data(hk + 737);
    const auto *hk_738 = buffer.data(hk + 738);
    const auto *hk_739 = buffer.data(hk + 739);
    const auto *hk_740 = buffer.data(hk + 740);
    const auto *hk_741 = buffer.data(hk + 741);
    const auto *hk_742 = buffer.data(hk + 742);
    const auto *hk_743 = buffer.data(hk + 743);
    const auto *hk_744 = buffer.data(hk + 744);
    const auto *hk_745 = buffer.data(hk + 745);
    const auto *hk_746 = buffer.data(hk + 746);
    const auto *hk_747 = buffer.data(hk + 747);
    const auto *hk_748 = buffer.data(hk + 748);
    const auto *hk_749 = buffer.data(hk + 749);
    const auto *hk_750 = buffer.data(hk + 750);
    const auto *hk_751 = buffer.data(hk + 751);
    const auto *hk_752 = buffer.data(hk + 752);
    const auto *hk_753 = buffer.data(hk + 753);
    const auto *hk_754 = buffer.data(hk + 754);
    const auto *hk_755 = buffer.data(hk + 755);

#pragma omp simd aligned(hk_37, hk_42, hk_51, hk_64, hk_217, hk_222, hk_231, hk_244, hk_541, \
                         hk_546, hk_555, hk_568 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * hk_37[k]
                 - f_1 * hk_42[k]
                 + f_2 * hk_51[k]
                 - f_3 * hk_64[k]
                 - f_4 * hk_217[k]
                 + f_5 * hk_222[k]
                 - f_6 * hk_231[k]
                 + f_7 * hk_244[k]
                 + f_8 * hk_541[k]
                 - f_0 * hk_546[k]
                 + f_9 * hk_555[k]
                 - f_10 * hk_568[k];
    }

#pragma omp simd aligned(hk_40, hk_47, hk_58, hk_220, hk_227, hk_238, hk_544, hk_551, \
                         hk_562 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_1[k] = f_11 * hk_40[k]
                 - f_12 * hk_47[k]
                 + f_11 * hk_58[k]
                 - f_13 * hk_220[k]
                 + f_14 * hk_227[k]
                 - f_13 * hk_238[k]
                 + f_15 * hk_544[k]
                 - f_16 * hk_551[k]
                 + f_15 * hk_562[k];
    }

#pragma omp simd aligned(hk_37, hk_42, hk_44, hk_51, hk_53, hk_64, hk_66, hk_217, hk_222, \
                         hk_224, hk_231, hk_233, hk_244, hk_246, hk_541, hk_546, hk_548, \
                         hk_555, hk_557, hk_568, hk_570 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_2[k] = -f_17 * hk_37[k]
                 + f_17 * hk_42[k]
                 + f_18 * hk_44[k]
                 + f_19 * hk_51[k]
                 - f_20 * hk_53[k]
                 - f_21 * hk_64[k]
                 + f_22 * hk_66[k]
                 + f_23 * hk_217[k]
                 - f_23 * hk_222[k]
                 - f_20 * hk_224[k]
                 - f_24 * hk_231[k]
                 + f_25 * hk_233[k]
                 + f_26 * hk_244[k]
                 - f_27 * hk_246[k]
                 - f_21 * hk_541[k]
                 + f_21 * hk_546[k]
                 + f_22 * hk_548[k]
                 + f_28 * hk_555[k]
                 - f_27 * hk_557[k]
                 - f_29 * hk_568[k]
                 + f_30 * hk_570[k];
    }

#pragma omp simd aligned(hk_40, hk_49, hk_58, hk_60, hk_220, hk_229, hk_238, hk_240, hk_544, \
                         hk_553, hk_562, hk_564 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = -f_27 * hk_40[k]
                 + f_31 * hk_49[k]
                 + f_27 * hk_58[k]
                 - f_31 * hk_60[k]
                 + f_32 * hk_220[k]
                 - f_33 * hk_229[k]
                 - f_32 * hk_238[k]
                 + f_33 * hk_240[k]
                 - f_34 * hk_544[k]
                 + f_35 * hk_553[k]
                 + f_34 * hk_562[k]
                 - f_35 * hk_564[k];
    }

#pragma omp simd aligned(hk_37, hk_42, hk_44, hk_51, hk_53, hk_55, hk_64, hk_66, hk_68, \
                         hk_217, hk_222, hk_224, hk_231, hk_233, hk_235, hk_244, hk_246, \
                         hk_248, hk_541, hk_546, hk_548, hk_555, hk_557, hk_559, hk_568, \
                         hk_570, hk_572 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = f_36 * hk_37[k]
                 + f_37 * hk_42[k]
                 - f_38 * hk_44[k]
                 + f_39 * hk_51[k]
                 - f_40 * hk_53[k]
                 + f_41 * hk_55[k]
                 - f_39 * hk_64[k]
                 + f_42 * hk_66[k]
                 - f_43 * hk_68[k]
                 - f_44 * hk_217[k]
                 - f_45 * hk_222[k]
                 + f_46 * hk_224[k]
                 - f_47 * hk_231[k]
                 + f_41 * hk_233[k]
                 - f_48 * hk_235[k]
                 + f_47 * hk_244[k]
                 - f_40 * hk_246[k]
                 + f_49 * hk_248[k]
                 + f_50 * hk_541[k]
                 + f_39 * hk_546[k]
                 - f_51 * hk_548[k]
                 + f_52 * hk_555[k]
                 - f_53 * hk_557[k]
                 + f_54 * hk_559[k]
                 - f_52 * hk_568[k]
                 + f_55 * hk_570[k]
                 - f_56 * hk_572[k];
    }

#pragma omp simd aligned(hk_40, hk_47, hk_49, hk_58, hk_60, hk_62, hk_220, hk_227, hk_229, \
                         hk_238, hk_240, hk_242, hk_544, hk_551, hk_553, hk_562, hk_564, \
                         hk_566 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = f_57 * hk_40[k]
                 + f_58 * hk_47[k]
                 - f_59 * hk_49[k]
                 + f_57 * hk_58[k]
                 - f_59 * hk_60[k]
                 + f_60 * hk_62[k]
                 - f_58 * hk_220[k]
                 - f_61 * hk_227[k]
                 + f_62 * hk_229[k]
                 - f_58 * hk_238[k]
                 + f_62 * hk_240[k]
                 - f_63 * hk_242[k]
                 + f_64 * hk_544[k]
                 + f_65 * hk_551[k]
                 - f_66 * hk_553[k]
                 + f_64 * hk_562[k]
                 - f_66 * hk_564[k]
                 + f_67 * hk_566[k];
    }

#pragma omp simd aligned(hk_37, hk_42, hk_44, hk_51, hk_53, hk_55, hk_64, hk_66, hk_68, hk_70, \
                         hk_217, hk_222, hk_224, hk_231, hk_233, hk_235, hk_244, hk_246, \
                         hk_248, hk_250, hk_541, hk_546, hk_548, hk_555, hk_557, hk_559, \
                         hk_568, hk_570, hk_572, hk_574 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = -f_68 * hk_37[k]
                 - f_69 * hk_42[k]
                 + f_70 * hk_44[k]
                 - f_69 * hk_51[k]
                 + f_71 * hk_53[k]
                 - f_71 * hk_55[k]
                 - f_68 * hk_64[k]
                 + f_70 * hk_66[k]
                 - f_71 * hk_68[k]
                 + f_72 * hk_70[k]
                 + f_73 * hk_217[k]
                 + f_74 * hk_222[k]
                 - f_71 * hk_224[k]
                 + f_74 * hk_231[k]
                 - f_75 * hk_233[k]
                 + f_75 * hk_235[k]
                 + f_73 * hk_244[k]
                 - f_71 * hk_246[k]
                 + f_75 * hk_248[k]
                 - f_76 * hk_250[k]
                 - f_77 * hk_541[k]
                 - f_78 * hk_546[k]
                 + f_79 * hk_548[k]
                 - f_78 * hk_555[k]
                 + f_80 * hk_557[k]
                 - f_80 * hk_559[k]
                 - f_77 * hk_568[k]
                 + f_79 * hk_570[k]
                 - f_80 * hk_572[k]
                 + f_81 * hk_574[k];
    }

#pragma omp simd aligned(hk_38, hk_43, hk_45, hk_52, hk_54, hk_56, hk_65, hk_67, hk_69, hk_71, \
                         hk_218, hk_223, hk_225, hk_232, hk_234, hk_236, hk_245, hk_247, \
                         hk_249, hk_251, hk_542, hk_547, hk_549, hk_556, hk_558, hk_560, \
                         hk_569, hk_571, hk_573, hk_575 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = -f_82 * hk_38[k]
                 - f_83 * hk_43[k]
                 + f_84 * hk_45[k]
                 - f_83 * hk_52[k]
                 + f_85 * hk_54[k]
                 - f_86 * hk_56[k]
                 - f_82 * hk_65[k]
                 + f_84 * hk_67[k]
                 - f_86 * hk_69[k]
                 + f_87 * hk_71[k]
                 + f_88 * hk_218[k]
                 + f_84 * hk_223[k]
                 - f_85 * hk_225[k]
                 + f_84 * hk_232[k]
                 - f_89 * hk_234[k]
                 + f_90 * hk_236[k]
                 + f_88 * hk_245[k]
                 - f_85 * hk_247[k]
                 + f_90 * hk_249[k]
                 - f_91 * hk_251[k]
                 - f_92 * hk_542[k]
                 - f_93 * hk_547[k]
                 + f_94 * hk_549[k]
                 - f_93 * hk_556[k]
                 + f_95 * hk_558[k]
                 - f_96 * hk_560[k]
                 - f_92 * hk_569[k]
                 + f_94 * hk_571[k]
                 - f_96 * hk_573[k]
                 + f_97 * hk_575[k];
    }

#pragma omp simd aligned(hk_36, hk_39, hk_41, hk_46, hk_48, hk_50, hk_57, hk_59, hk_61, hk_63, \
                         hk_216, hk_219, hk_221, hk_226, hk_228, hk_230, hk_237, hk_239, \
                         hk_241, hk_243, hk_540, hk_543, hk_545, hk_550, hk_552, hk_554, \
                         hk_561, hk_563, hk_565, hk_567 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_8[k] = -f_68 * hk_36[k]
                 - f_69 * hk_39[k]
                 + f_70 * hk_41[k]
                 - f_69 * hk_46[k]
                 + f_71 * hk_48[k]
                 - f_71 * hk_50[k]
                 - f_68 * hk_57[k]
                 + f_70 * hk_59[k]
                 - f_71 * hk_61[k]
                 + f_72 * hk_63[k]
                 + f_73 * hk_216[k]
                 + f_74 * hk_219[k]
                 - f_71 * hk_221[k]
                 + f_74 * hk_226[k]
                 - f_75 * hk_228[k]
                 + f_75 * hk_230[k]
                 + f_73 * hk_237[k]
                 - f_71 * hk_239[k]
                 + f_75 * hk_241[k]
                 - f_76 * hk_243[k]
                 - f_77 * hk_540[k]
                 - f_78 * hk_543[k]
                 + f_79 * hk_545[k]
                 - f_78 * hk_550[k]
                 + f_80 * hk_552[k]
                 - f_80 * hk_554[k]
                 - f_77 * hk_561[k]
                 + f_79 * hk_563[k]
                 - f_80 * hk_565[k]
                 + f_81 * hk_567[k];
    }

#pragma omp simd aligned(hk_38, hk_43, hk_45, hk_52, hk_56, hk_65, hk_67, hk_69, hk_218, \
                         hk_223, hk_225, hk_232, hk_236, hk_245, hk_247, hk_249, hk_542, \
                         hk_547, hk_549, hk_556, hk_560, hk_569, hk_571, \
                         hk_573 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_9[k] = f_98 * hk_38[k]
                 + f_98 * hk_43[k]
                 - f_99 * hk_45[k]
                 - f_98 * hk_52[k]
                 + f_100 * hk_56[k]
                 - f_98 * hk_65[k]
                 + f_99 * hk_67[k]
                 - f_100 * hk_69[k]
                 - f_57 * hk_218[k]
                 - f_57 * hk_223[k]
                 + f_59 * hk_225[k]
                 + f_57 * hk_232[k]
                 - f_60 * hk_236[k]
                 + f_57 * hk_245[k]
                 - f_59 * hk_247[k]
                 + f_60 * hk_249[k]
                 + f_101 * hk_542[k]
                 + f_101 * hk_547[k]
                 - f_102 * hk_549[k]
                 - f_101 * hk_556[k]
                 + f_103 * hk_560[k]
                 - f_101 * hk_569[k]
                 + f_102 * hk_571[k]
                 - f_103 * hk_573[k];
    }

#pragma omp simd aligned(hk_36, hk_39, hk_41, hk_46, hk_48, hk_50, hk_57, hk_59, hk_61, \
                         hk_216, hk_219, hk_221, hk_226, hk_228, hk_230, hk_237, hk_239, \
                         hk_241, hk_540, hk_543, hk_545, hk_550, hk_552, hk_554, hk_561, \
                         hk_563, hk_565 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = f_39 * hk_36[k]
                  - f_39 * hk_39[k]
                  - f_42 * hk_41[k]
                  - f_37 * hk_46[k]
                  + f_40 * hk_48[k]
                  + f_43 * hk_50[k]
                  - f_36 * hk_57[k]
                  + f_38 * hk_59[k]
                  - f_41 * hk_61[k]
                  - f_47 * hk_216[k]
                  + f_47 * hk_219[k]
                  + f_40 * hk_221[k]
                  + f_45 * hk_226[k]
                  - f_41 * hk_228[k]
                  - f_49 * hk_230[k]
                  + f_44 * hk_237[k]
                  - f_46 * hk_239[k]
                  + f_48 * hk_241[k]
                  + f_52 * hk_540[k]
                  - f_52 * hk_543[k]
                  - f_55 * hk_545[k]
                  - f_39 * hk_550[k]
                  + f_53 * hk_552[k]
                  + f_56 * hk_554[k]
                  - f_50 * hk_561[k]
                  + f_51 * hk_563[k]
                  - f_54 * hk_565[k];
    }

#pragma omp simd aligned(hk_38, hk_43, hk_45, hk_52, hk_54, hk_65, hk_67, hk_218, hk_223, \
                         hk_225, hk_232, hk_234, hk_245, hk_247, hk_542, hk_547, hk_549, \
                         hk_556, hk_558, hk_569, hk_571 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_11[k] = -f_104 * hk_38[k]
                  + f_105 * hk_43[k]
                  + f_106 * hk_45[k]
                  + f_105 * hk_52[k]
                  - f_20 * hk_54[k]
                  - f_104 * hk_65[k]
                  + f_106 * hk_67[k]
                  + f_22 * hk_218[k]
                  - f_18 * hk_223[k]
                  - f_107 * hk_225[k]
                  - f_18 * hk_232[k]
                  + f_25 * hk_234[k]
                  + f_22 * hk_245[k]
                  - f_107 * hk_247[k]
                  - f_108 * hk_542[k]
                  + f_104 * hk_547[k]
                  + f_109 * hk_549[k]
                  + f_104 * hk_556[k]
                  - f_27 * hk_558[k]
                  - f_108 * hk_569[k]
                  + f_109 * hk_571[k];
    }

#pragma omp simd aligned(hk_36, hk_39, hk_41, hk_46, hk_48, hk_57, hk_59, hk_216, hk_219, \
                         hk_221, hk_226, hk_228, hk_237, hk_239, hk_540, hk_543, hk_545, \
                         hk_550, hk_552, hk_561, hk_563 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_12[k] = -f_21 * hk_36[k]
                  + f_19 * hk_39[k]
                  + f_22 * hk_41[k]
                  + f_17 * hk_46[k]
                  - f_20 * hk_48[k]
                  - f_17 * hk_57[k]
                  + f_18 * hk_59[k]
                  + f_26 * hk_216[k]
                  - f_24 * hk_219[k]
                  - f_27 * hk_221[k]
                  - f_23 * hk_226[k]
                  + f_25 * hk_228[k]
                  + f_23 * hk_237[k]
                  - f_20 * hk_239[k]
                  - f_29 * hk_540[k]
                  + f_28 * hk_543[k]
                  + f_30 * hk_545[k]
                  + f_21 * hk_550[k]
                  - f_27 * hk_552[k]
                  - f_21 * hk_561[k]
                  + f_22 * hk_563[k];
    }

#pragma omp simd aligned(hk_38, hk_43, hk_52, hk_65, hk_218, hk_223, hk_232, hk_245, hk_542, \
                         hk_547, hk_556, hk_569 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = f_110 * hk_38[k]
                  - f_111 * hk_43[k]
                  + f_111 * hk_52[k]
                  - f_110 * hk_65[k]
                  - f_112 * hk_218[k]
                  + f_113 * hk_223[k]
                  - f_113 * hk_232[k]
                  + f_112 * hk_245[k]
                  + f_114 * hk_542[k]
                  - f_115 * hk_547[k]
                  + f_115 * hk_556[k]
                  - f_114 * hk_569[k];
    }

#pragma omp simd aligned(hk_36, hk_39, hk_46, hk_57, hk_216, hk_219, hk_226, hk_237, hk_540, \
                         hk_543, hk_550, hk_561 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = f_3 * hk_36[k]
                  - f_2 * hk_39[k]
                  + f_1 * hk_46[k]
                  - f_0 * hk_57[k]
                  - f_7 * hk_216[k]
                  + f_6 * hk_219[k]
                  - f_5 * hk_226[k]
                  + f_4 * hk_237[k]
                  + f_10 * hk_540[k]
                  - f_9 * hk_543[k]
                  + f_0 * hk_550[k]
                  - f_8 * hk_561[k];
    }

#pragma omp simd aligned(hk_145, hk_148, hk_150, hk_155, hk_159, hk_166, hk_172, hk_397, \
                         hk_400, hk_402, hk_407, hk_411, hk_418, \
                         hk_424 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_15[k] = f_116 * hk_145[k]
                  - f_117 * hk_150[k]
                  + f_118 * hk_159[k]
                  - f_119 * hk_172[k]
                  - f_116 * hk_397[k]
                  + f_117 * hk_402[k]
                  - f_118 * hk_411[k]
                  + f_119 * hk_424[k];

        g_16[k] = f_120 * hk_148[k]
                  - f_121 * hk_155[k]
                  + f_120 * hk_166[k]
                  - f_120 * hk_400[k]
                  + f_121 * hk_407[k]
                  - f_120 * hk_418[k];
    }

#pragma omp simd aligned(hk_145, hk_150, hk_152, hk_159, hk_161, hk_172, hk_174, hk_397, \
                         hk_402, hk_404, hk_411, hk_413, hk_424, \
                         hk_426 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = -f_122 * hk_145[k]
                  + f_122 * hk_150[k]
                  + f_123 * hk_152[k]
                  + f_124 * hk_159[k]
                  - f_125 * hk_161[k]
                  - f_126 * hk_172[k]
                  + f_127 * hk_174[k]
                  + f_122 * hk_397[k]
                  - f_122 * hk_402[k]
                  - f_123 * hk_404[k]
                  - f_124 * hk_411[k]
                  + f_125 * hk_413[k]
                  + f_126 * hk_424[k]
                  - f_127 * hk_426[k];
    }

#pragma omp simd aligned(hk_148, hk_157, hk_166, hk_168, hk_400, hk_409, hk_418, \
                         hk_420 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = -f_128 * hk_148[k]
                  + f_129 * hk_157[k]
                  + f_128 * hk_166[k]
                  - f_129 * hk_168[k]
                  + f_128 * hk_400[k]
                  - f_129 * hk_409[k]
                  - f_128 * hk_418[k]
                  + f_129 * hk_420[k];
    }

#pragma omp simd aligned(hk_145, hk_150, hk_152, hk_159, hk_161, hk_163, hk_172, hk_174, \
                         hk_176, hk_397, hk_402, hk_404, hk_411, hk_413, hk_415, hk_424, \
                         hk_426, hk_428 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = f_130 * hk_145[k]
                  + f_131 * hk_150[k]
                  - f_132 * hk_152[k]
                  + f_133 * hk_159[k]
                  - f_134 * hk_161[k]
                  + f_135 * hk_163[k]
                  - f_133 * hk_172[k]
                  + f_136 * hk_174[k]
                  - f_137 * hk_176[k]
                  - f_130 * hk_397[k]
                  - f_131 * hk_402[k]
                  + f_132 * hk_404[k]
                  - f_133 * hk_411[k]
                  + f_134 * hk_413[k]
                  - f_135 * hk_415[k]
                  + f_133 * hk_424[k]
                  - f_136 * hk_426[k]
                  + f_137 * hk_428[k];
    }

#pragma omp simd aligned(hk_148, hk_155, hk_157, hk_166, hk_168, hk_170, hk_400, hk_407, \
                         hk_409, hk_418, hk_420, hk_422 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = f_138 * hk_148[k]
                  + f_139 * hk_155[k]
                  - f_140 * hk_157[k]
                  + f_138 * hk_166[k]
                  - f_140 * hk_168[k]
                  + f_141 * hk_170[k]
                  - f_138 * hk_400[k]
                  - f_139 * hk_407[k]
                  + f_140 * hk_409[k]
                  - f_138 * hk_418[k]
                  + f_140 * hk_420[k]
                  - f_141 * hk_422[k];
    }

#pragma omp simd aligned(hk_145, hk_150, hk_152, hk_159, hk_161, hk_163, hk_172, hk_174, \
                         hk_176, hk_178, hk_397, hk_402, hk_404, hk_411, hk_413, hk_415, \
                         hk_424, hk_426, hk_428, hk_430 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = -f_142 * hk_145[k]
                  - f_143 * hk_150[k]
                  + f_144 * hk_152[k]
                  - f_143 * hk_159[k]
                  + f_145 * hk_161[k]
                  - f_145 * hk_163[k]
                  - f_142 * hk_172[k]
                  + f_144 * hk_174[k]
                  - f_145 * hk_176[k]
                  + f_146 * hk_178[k]
                  + f_142 * hk_397[k]
                  + f_143 * hk_402[k]
                  - f_144 * hk_404[k]
                  + f_143 * hk_411[k]
                  - f_145 * hk_413[k]
                  + f_145 * hk_415[k]
                  + f_142 * hk_424[k]
                  - f_144 * hk_426[k]
                  + f_145 * hk_428[k]
                  - f_146 * hk_430[k];
    }

#pragma omp simd aligned(hk_146, hk_151, hk_153, hk_160, hk_162, hk_164, hk_173, hk_175, \
                         hk_177, hk_179, hk_398, hk_403, hk_405, hk_412, hk_414, hk_416, \
                         hk_425, hk_427, hk_429, hk_431 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = -f_147 * hk_146[k]
                  - f_148 * hk_151[k]
                  + f_149 * hk_153[k]
                  - f_148 * hk_160[k]
                  + f_150 * hk_162[k]
                  - f_151 * hk_164[k]
                  - f_147 * hk_173[k]
                  + f_149 * hk_175[k]
                  - f_151 * hk_177[k]
                  + f_152 * hk_179[k]
                  + f_147 * hk_398[k]
                  + f_148 * hk_403[k]
                  - f_149 * hk_405[k]
                  + f_148 * hk_412[k]
                  - f_150 * hk_414[k]
                  + f_151 * hk_416[k]
                  + f_147 * hk_425[k]
                  - f_149 * hk_427[k]
                  + f_151 * hk_429[k]
                  - f_152 * hk_431[k];
    }

#pragma omp simd aligned(hk_144, hk_147, hk_149, hk_154, hk_156, hk_158, hk_165, hk_167, \
                         hk_169, hk_171, hk_396, hk_399, hk_401, hk_406, hk_408, hk_410, \
                         hk_417, hk_419, hk_421, hk_423 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = -f_142 * hk_144[k]
                  - f_143 * hk_147[k]
                  + f_144 * hk_149[k]
                  - f_143 * hk_154[k]
                  + f_145 * hk_156[k]
                  - f_145 * hk_158[k]
                  - f_142 * hk_165[k]
                  + f_144 * hk_167[k]
                  - f_145 * hk_169[k]
                  + f_146 * hk_171[k]
                  + f_142 * hk_396[k]
                  + f_143 * hk_399[k]
                  - f_144 * hk_401[k]
                  + f_143 * hk_406[k]
                  - f_145 * hk_408[k]
                  + f_145 * hk_410[k]
                  + f_142 * hk_417[k]
                  - f_144 * hk_419[k]
                  + f_145 * hk_421[k]
                  - f_146 * hk_423[k];
    }

#pragma omp simd aligned(hk_146, hk_151, hk_153, hk_160, hk_164, hk_173, hk_175, hk_177, \
                         hk_398, hk_403, hk_405, hk_412, hk_416, hk_425, hk_427, \
                         hk_429 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = f_153 * hk_146[k]
                  + f_153 * hk_151[k]
                  - f_154 * hk_153[k]
                  - f_153 * hk_160[k]
                  + f_155 * hk_164[k]
                  - f_153 * hk_173[k]
                  + f_154 * hk_175[k]
                  - f_155 * hk_177[k]
                  - f_153 * hk_398[k]
                  - f_153 * hk_403[k]
                  + f_154 * hk_405[k]
                  + f_153 * hk_412[k]
                  - f_155 * hk_416[k]
                  + f_153 * hk_425[k]
                  - f_154 * hk_427[k]
                  + f_155 * hk_429[k];
    }

#pragma omp simd aligned(hk_144, hk_147, hk_149, hk_154, hk_156, hk_158, hk_165, hk_167, \
                         hk_169, hk_396, hk_399, hk_401, hk_406, hk_408, hk_410, hk_417, \
                         hk_419, hk_421 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = f_133 * hk_144[k]
                  - f_133 * hk_147[k]
                  - f_136 * hk_149[k]
                  - f_131 * hk_154[k]
                  + f_134 * hk_156[k]
                  + f_137 * hk_158[k]
                  - f_130 * hk_165[k]
                  + f_132 * hk_167[k]
                  - f_135 * hk_169[k]
                  - f_133 * hk_396[k]
                  + f_133 * hk_399[k]
                  + f_136 * hk_401[k]
                  + f_131 * hk_406[k]
                  - f_134 * hk_408[k]
                  - f_137 * hk_410[k]
                  + f_130 * hk_417[k]
                  - f_132 * hk_419[k]
                  + f_135 * hk_421[k];
    }

#pragma omp simd aligned(hk_146, hk_151, hk_153, hk_160, hk_162, hk_173, hk_175, hk_398, \
                         hk_403, hk_405, hk_412, hk_414, hk_425, \
                         hk_427 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = -f_156 * hk_146[k]
                  + f_157 * hk_151[k]
                  + f_158 * hk_153[k]
                  + f_157 * hk_160[k]
                  - f_125 * hk_162[k]
                  - f_156 * hk_173[k]
                  + f_158 * hk_175[k]
                  + f_156 * hk_398[k]
                  - f_157 * hk_403[k]
                  - f_158 * hk_405[k]
                  - f_157 * hk_412[k]
                  + f_125 * hk_414[k]
                  + f_156 * hk_425[k]
                  - f_158 * hk_427[k];
    }

#pragma omp simd aligned(hk_144, hk_147, hk_149, hk_154, hk_156, hk_165, hk_167, hk_396, \
                         hk_399, hk_401, hk_406, hk_408, hk_417, \
                         hk_419 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = -f_126 * hk_144[k]
                  + f_124 * hk_147[k]
                  + f_127 * hk_149[k]
                  + f_122 * hk_154[k]
                  - f_125 * hk_156[k]
                  - f_122 * hk_165[k]
                  + f_123 * hk_167[k]
                  + f_126 * hk_396[k]
                  - f_124 * hk_399[k]
                  - f_127 * hk_401[k]
                  - f_122 * hk_406[k]
                  + f_125 * hk_408[k]
                  + f_122 * hk_417[k]
                  - f_123 * hk_419[k];
    }

#pragma omp simd aligned(hk_146, hk_151, hk_160, hk_173, hk_398, hk_403, hk_412, \
                         hk_425 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = f_159 * hk_146[k]
                  - f_160 * hk_151[k]
                  + f_160 * hk_160[k]
                  - f_159 * hk_173[k]
                  - f_159 * hk_398[k]
                  + f_160 * hk_403[k]
                  - f_160 * hk_412[k]
                  + f_159 * hk_425[k];
    }

#pragma omp simd aligned(hk_144, hk_147, hk_154, hk_165, hk_396, hk_399, hk_406, \
                         hk_417 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = f_119 * hk_144[k]
                  - f_118 * hk_147[k]
                  + f_117 * hk_154[k]
                  - f_116 * hk_165[k]
                  - f_119 * hk_396[k]
                  + f_118 * hk_399[k]
                  - f_117 * hk_406[k]
                  + f_116 * hk_417[k];
    }

#pragma omp simd aligned(hk_37, hk_42, hk_51, hk_64, hk_217, hk_222, hk_231, hk_244, hk_289, \
                         hk_294, hk_303, hk_316, hk_541, hk_546, hk_555, hk_568, hk_613, \
                         hk_618, hk_627, hk_640 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = -f_161 * hk_37[k]
                  + f_162 * hk_42[k]
                  - f_163 * hk_51[k]
                  + f_164 * hk_64[k]
                  - f_165 * hk_217[k]
                  + f_166 * hk_222[k]
                  - f_167 * hk_231[k]
                  + f_168 * hk_244[k]
                  + f_169 * hk_289[k]
                  - f_170 * hk_294[k]
                  + f_171 * hk_303[k]
                  - f_172 * hk_316[k]
                  + f_173 * hk_541[k]
                  - f_174 * hk_546[k]
                  + f_161 * hk_555[k]
                  - f_175 * hk_568[k]
                  - f_176 * hk_613[k]
                  + f_177 * hk_618[k]
                  - f_169 * hk_627[k]
                  + f_178 * hk_640[k];
    }

#pragma omp simd aligned(hk_40, hk_47, hk_58, hk_220, hk_227, hk_238, hk_292, hk_299, hk_310, \
                         hk_544, hk_551, hk_562, hk_616, hk_623, \
                         hk_634 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = -f_179 * hk_40[k]
                  + f_180 * hk_47[k]
                  - f_179 * hk_58[k]
                  - f_181 * hk_220[k]
                  + f_182 * hk_227[k]
                  - f_181 * hk_238[k]
                  + f_183 * hk_292[k]
                  - f_184 * hk_299[k]
                  + f_183 * hk_310[k]
                  + f_185 * hk_544[k]
                  - f_186 * hk_551[k]
                  + f_185 * hk_562[k]
                  - f_187 * hk_616[k]
                  + f_188 * hk_623[k]
                  - f_187 * hk_634[k];
    }

#pragma omp simd aligned(hk_37, hk_42, hk_44, hk_51, hk_53, hk_64, hk_66, hk_217, hk_222, \
                         hk_224, hk_231, hk_233, hk_244, hk_246, hk_289, hk_294, hk_296, \
                         hk_303, hk_305, hk_316, hk_318, hk_541, hk_546, hk_548, hk_555, \
                         hk_557, hk_568, hk_570, hk_613, hk_618, hk_620, hk_627, hk_629, \
                         hk_640, hk_642 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = f_189 * hk_37[k]
                  - f_189 * hk_42[k]
                  - f_190 * hk_44[k]
                  - f_191 * hk_51[k]
                  + f_192 * hk_53[k]
                  + f_193 * hk_64[k]
                  - f_194 * hk_66[k]
                  + f_195 * hk_217[k]
                  - f_195 * hk_222[k]
                  - f_196 * hk_224[k]
                  - f_197 * hk_231[k]
                  + f_198 * hk_233[k]
                  + f_199 * hk_244[k]
                  - f_200 * hk_246[k]
                  - f_196 * hk_289[k]
                  + f_196 * hk_294[k]
                  + f_201 * hk_296[k]
                  + f_202 * hk_303[k]
                  - f_203 * hk_305[k]
                  - f_200 * hk_316[k]
                  + f_204 * hk_318[k]
                  - f_205 * hk_541[k]
                  + f_205 * hk_546[k]
                  + f_206 * hk_548[k]
                  + f_207 * hk_555[k]
                  - f_196 * hk_557[k]
                  - f_208 * hk_568[k]
                  + f_209 * hk_570[k]
                  + f_210 * hk_613[k]
                  - f_210 * hk_618[k]
                  - f_211 * hk_620[k]
                  - f_212 * hk_627[k]
                  + f_213 * hk_629[k]
                  + f_214 * hk_640[k]
                  - f_215 * hk_642[k];
    }

#pragma omp simd aligned(hk_40, hk_49, hk_58, hk_60, hk_220, hk_229, hk_238, hk_240, hk_292, \
                         hk_301, hk_310, hk_312, hk_544, hk_553, hk_562, hk_564, hk_616, \
                         hk_625, hk_634, hk_636 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = f_212 * hk_40[k]
                  - f_198 * hk_49[k]
                  - f_212 * hk_58[k]
                  + f_198 * hk_60[k]
                  + f_216 * hk_220[k]
                  - f_217 * hk_229[k]
                  - f_216 * hk_238[k]
                  + f_217 * hk_240[k]
                  - f_218 * hk_292[k]
                  + f_219 * hk_301[k]
                  + f_218 * hk_310[k]
                  - f_219 * hk_312[k]
                  - f_200 * hk_544[k]
                  + f_220 * hk_553[k]
                  + f_200 * hk_562[k]
                  - f_220 * hk_564[k]
                  + f_221 * hk_616[k]
                  - f_222 * hk_625[k]
                  - f_221 * hk_634[k]
                  + f_222 * hk_636[k];
    }

#pragma omp simd aligned(hk_37, hk_42, hk_44, hk_51, hk_53, hk_55, hk_64, hk_66, hk_68, \
                         hk_217, hk_222, hk_224, hk_231, hk_233, hk_235, hk_244, hk_246, \
                         hk_248, hk_289, hk_294, hk_296, hk_303, hk_305, hk_307, hk_316, \
                         hk_318, hk_320, hk_541, hk_546, hk_548, hk_555, hk_557, hk_559, \
                         hk_568, hk_570, hk_572, hk_613, hk_618, hk_620, hk_627, hk_629, \
                         hk_631, hk_640, hk_642, hk_644 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = -f_223 * hk_37[k]
                  - f_224 * hk_42[k]
                  + f_225 * hk_44[k]
                  - f_226 * hk_51[k]
                  + f_153 * hk_53[k]
                  - f_138 * hk_55[k]
                  + f_226 * hk_64[k]
                  - f_227 * hk_66[k]
                  + f_228 * hk_68[k]
                  - f_229 * hk_217[k]
                  - f_230 * hk_222[k]
                  + f_153 * hk_224[k]
                  - f_231 * hk_231[k]
                  + f_228 * hk_233[k]
                  - f_232 * hk_235[k]
                  + f_231 * hk_244[k]
                  - f_233 * hk_246[k]
                  + f_234 * hk_248[k]
                  + f_235 * hk_289[k]
                  + f_153 * hk_294[k]
                  - f_236 * hk_296[k]
                  + f_237 * hk_303[k]
                  - f_238 * hk_305[k]
                  + f_239 * hk_307[k]
                  - f_237 * hk_316[k]
                  + f_139 * hk_318[k]
                  - f_154 * hk_320[k]
                  + f_226 * hk_541[k]
                  + f_240 * hk_546[k]
                  - f_227 * hk_548[k]
                  + f_241 * hk_555[k]
                  - f_233 * hk_557[k]
                  + f_228 * hk_559[k]
                  - f_241 * hk_568[k]
                  + f_242 * hk_570[k]
                  - f_243 * hk_572[k]
                  - f_237 * hk_613[k]
                  - f_233 * hk_618[k]
                  + f_139 * hk_620[k]
                  - f_244 * hk_627[k]
                  + f_245 * hk_629[k]
                  - f_154 * hk_631[k]
                  + f_244 * hk_640[k]
                  - f_232 * hk_642[k]
                  + f_246 * hk_644[k];
    }

#pragma omp simd aligned(hk_40, hk_47, hk_49, hk_58, hk_60, hk_62, hk_220, hk_227, hk_229, \
                         hk_238, hk_240, hk_242, hk_292, hk_299, hk_301, hk_310, hk_312, \
                         hk_314, hk_544, hk_551, hk_553, hk_562, hk_564, hk_566, hk_616, \
                         hk_623, hk_625, hk_634, hk_636, hk_638 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = -f_247 * hk_40[k]
                  - f_131 * hk_47[k]
                  + f_248 * hk_49[k]
                  - f_247 * hk_58[k]
                  + f_248 * hk_60[k]
                  - f_249 * hk_62[k]
                  - f_250 * hk_220[k]
                  - f_251 * hk_227[k]
                  + f_252 * hk_229[k]
                  - f_250 * hk_238[k]
                  + f_252 * hk_240[k]
                  - f_253 * hk_242[k]
                  + f_136 * hk_292[k]
                  + f_134 * hk_299[k]
                  - f_254 * hk_301[k]
                  + f_136 * hk_310[k]
                  - f_254 * hk_312[k]
                  + f_255 * hk_314[k]
                  + f_256 * hk_544[k]
                  + f_250 * hk_551[k]
                  - f_257 * hk_553[k]
                  + f_256 * hk_562[k]
                  - f_257 * hk_564[k]
                  + f_258 * hk_566[k]
                  - f_259 * hk_616[k]
                  - f_248 * hk_623[k]
                  + f_260 * hk_625[k]
                  - f_259 * hk_634[k]
                  + f_260 * hk_636[k]
                  - f_261 * hk_638[k];
    }

#pragma omp simd aligned(hk_37, hk_42, hk_44, hk_51, hk_53, hk_55, hk_64, hk_66, hk_68, hk_70, \
                         hk_217, hk_222, hk_224, hk_231, hk_233, hk_235, hk_244, hk_246, \
                         hk_248, hk_250, hk_289, hk_294, hk_296, hk_303, hk_305, hk_307, \
                         hk_316, hk_318, hk_320, hk_322, hk_541, hk_546, hk_548, hk_555, \
                         hk_557, hk_559, hk_568, hk_570, hk_572, hk_574, hk_613, hk_618, \
                         hk_620, hk_627, hk_629, hk_631, hk_640, hk_642, hk_644, \
                         hk_646 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_36[k] = f_262 * hk_37[k]
                  + f_263 * hk_42[k]
                  - f_264 * hk_44[k]
                  + f_263 * hk_51[k]
                  - f_265 * hk_53[k]
                  + f_265 * hk_55[k]
                  + f_262 * hk_64[k]
                  - f_264 * hk_66[k]
                  + f_265 * hk_68[k]
                  - f_266 * hk_70[k]
                  + f_267 * hk_217[k]
                  + f_268 * hk_222[k]
                  - f_269 * hk_224[k]
                  + f_268 * hk_231[k]
                  - f_270 * hk_233[k]
                  + f_270 * hk_235[k]
                  + f_267 * hk_244[k]
                  - f_269 * hk_246[k]
                  + f_270 * hk_248[k]
                  - f_271 * hk_250[k]
                  - f_272 * hk_289[k]
                  - f_264 * hk_294[k]
                  + f_273 * hk_296[k]
                  - f_264 * hk_303[k]
                  + f_274 * hk_305[k]
                  - f_274 * hk_307[k]
                  - f_272 * hk_316[k]
                  + f_273 * hk_318[k]
                  - f_274 * hk_320[k]
                  + f_275 * hk_322[k]
                  - f_276 * hk_541[k]
                  - f_262 * hk_546[k]
                  + f_272 * hk_548[k]
                  - f_262 * hk_555[k]
                  + f_269 * hk_557[k]
                  - f_269 * hk_559[k]
                  - f_276 * hk_568[k]
                  + f_272 * hk_570[k]
                  - f_269 * hk_572[k]
                  + f_277 * hk_574[k]
                  + f_278 * hk_613[k]
                  + f_272 * hk_618[k]
                  - f_279 * hk_620[k]
                  + f_272 * hk_627[k]
                  - f_280 * hk_629[k]
                  + f_280 * hk_631[k]
                  + f_278 * hk_640[k]
                  - f_279 * hk_642[k]
                  + f_280 * hk_644[k]
                  - f_281 * hk_646[k];
    }

#pragma omp simd aligned(hk_38, hk_43, hk_45, hk_52, hk_54, hk_56, hk_65, hk_67, hk_69, hk_71, \
                         hk_218, hk_223, hk_225, hk_232, hk_234, hk_236, hk_245, hk_247, \
                         hk_249, hk_251, hk_290, hk_295, hk_297, hk_304, hk_306, hk_308, \
                         hk_317, hk_319, hk_321, hk_323, hk_542, hk_547, hk_549, hk_556, \
                         hk_558, hk_560, hk_569, hk_571, hk_573, hk_575, hk_614, hk_619, \
                         hk_621, hk_628, hk_630, hk_632, hk_641, hk_643, hk_645, \
                         hk_647 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_37[k] = f_282 * hk_38[k]
                  + f_283 * hk_43[k]
                  - f_284 * hk_45[k]
                  + f_283 * hk_52[k]
                  - f_285 * hk_54[k]
                  + f_286 * hk_56[k]
                  + f_282 * hk_65[k]
                  - f_284 * hk_67[k]
                  + f_286 * hk_69[k]
                  - f_287 * hk_71[k]
                  + f_288 * hk_218[k]
                  + f_289 * hk_223[k]
                  - f_290 * hk_225[k]
                  + f_289 * hk_232[k]
                  - f_291 * hk_234[k]
                  + f_292 * hk_236[k]
                  + f_288 * hk_245[k]
                  - f_290 * hk_247[k]
                  + f_292 * hk_249[k]
                  - f_293 * hk_251[k]
                  - f_291 * hk_290[k]
                  - f_294 * hk_295[k]
                  + f_295 * hk_297[k]
                  - f_294 * hk_304[k]
                  + f_296 * hk_306[k]
                  - f_297 * hk_308[k]
                  - f_291 * hk_317[k]
                  + f_295 * hk_319[k]
                  - f_297 * hk_321[k]
                  + f_298 * hk_323[k]
                  - f_299 * hk_542[k]
                  - f_282 * hk_547[k]
                  + f_289 * hk_549[k]
                  - f_282 * hk_556[k]
                  + f_290 * hk_558[k]
                  - f_300 * hk_560[k]
                  - f_299 * hk_569[k]
                  + f_289 * hk_571[k]
                  - f_300 * hk_573[k]
                  + f_301 * hk_575[k]
                  + f_302 * hk_614[k]
                  + f_291 * hk_619[k]
                  - f_303 * hk_621[k]
                  + f_291 * hk_628[k]
                  - f_304 * hk_630[k]
                  + f_305 * hk_632[k]
                  + f_302 * hk_641[k]
                  - f_303 * hk_643[k]
                  + f_305 * hk_645[k]
                  - f_306 * hk_647[k];
    }

#pragma omp simd aligned(hk_36, hk_39, hk_41, hk_46, hk_48, hk_50, hk_57, hk_59, hk_61, hk_63, \
                         hk_216, hk_219, hk_221, hk_226, hk_228, hk_230, hk_237, hk_239, \
                         hk_241, hk_243, hk_288, hk_291, hk_293, hk_298, hk_300, hk_302, \
                         hk_309, hk_311, hk_313, hk_315, hk_540, hk_543, hk_545, hk_550, \
                         hk_552, hk_554, hk_561, hk_563, hk_565, hk_567, hk_612, hk_615, \
                         hk_617, hk_622, hk_624, hk_626, hk_633, hk_635, hk_637, \
                         hk_639 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = f_262 * hk_36[k]
                  + f_263 * hk_39[k]
                  - f_264 * hk_41[k]
                  + f_263 * hk_46[k]
                  - f_265 * hk_48[k]
                  + f_265 * hk_50[k]
                  + f_262 * hk_57[k]
                  - f_264 * hk_59[k]
                  + f_265 * hk_61[k]
                  - f_266 * hk_63[k]
                  + f_267 * hk_216[k]
                  + f_268 * hk_219[k]
                  - f_269 * hk_221[k]
                  + f_268 * hk_226[k]
                  - f_270 * hk_228[k]
                  + f_270 * hk_230[k]
                  + f_267 * hk_237[k]
                  - f_269 * hk_239[k]
                  + f_270 * hk_241[k]
                  - f_271 * hk_243[k]
                  - f_272 * hk_288[k]
                  - f_264 * hk_291[k]
                  + f_273 * hk_293[k]
                  - f_264 * hk_298[k]
                  + f_274 * hk_300[k]
                  - f_274 * hk_302[k]
                  - f_272 * hk_309[k]
                  + f_273 * hk_311[k]
                  - f_274 * hk_313[k]
                  + f_275 * hk_315[k]
                  - f_276 * hk_540[k]
                  - f_262 * hk_543[k]
                  + f_272 * hk_545[k]
                  - f_262 * hk_550[k]
                  + f_269 * hk_552[k]
                  - f_269 * hk_554[k]
                  - f_276 * hk_561[k]
                  + f_272 * hk_563[k]
                  - f_269 * hk_565[k]
                  + f_277 * hk_567[k]
                  + f_278 * hk_612[k]
                  + f_272 * hk_615[k]
                  - f_279 * hk_617[k]
                  + f_272 * hk_622[k]
                  - f_280 * hk_624[k]
                  + f_280 * hk_626[k]
                  + f_278 * hk_633[k]
                  - f_279 * hk_635[k]
                  + f_280 * hk_637[k]
                  - f_281 * hk_639[k];
    }

#pragma omp simd aligned(hk_38, hk_43, hk_45, hk_52, hk_56, hk_65, hk_67, hk_69, hk_218, \
                         hk_223, hk_225, hk_232, hk_236, hk_245, hk_247, hk_249, hk_290, \
                         hk_295, hk_297, hk_304, hk_308, hk_317, hk_319, hk_321, hk_542, \
                         hk_547, hk_549, hk_556, hk_560, hk_569, hk_571, hk_573, hk_614, \
                         hk_619, hk_621, hk_628, hk_632, hk_641, hk_643, \
                         hk_645 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_39[k] = -f_307 * hk_38[k]
                  - f_307 * hk_43[k]
                  + f_259 * hk_45[k]
                  + f_307 * hk_52[k]
                  - f_308 * hk_56[k]
                  + f_307 * hk_65[k]
                  - f_259 * hk_67[k]
                  + f_308 * hk_69[k]
                  - f_256 * hk_218[k]
                  - f_256 * hk_223[k]
                  + f_257 * hk_225[k]
                  + f_256 * hk_232[k]
                  - f_258 * hk_236[k]
                  + f_256 * hk_245[k]
                  - f_257 * hk_247[k]
                  + f_258 * hk_249[k]
                  + f_309 * hk_290[k]
                  + f_309 * hk_295[k]
                  - f_310 * hk_297[k]
                  - f_309 * hk_304[k]
                  + f_311 * hk_308[k]
                  - f_309 * hk_317[k]
                  + f_310 * hk_319[k]
                  - f_311 * hk_321[k]
                  + f_312 * hk_542[k]
                  + f_312 * hk_547[k]
                  - f_313 * hk_549[k]
                  - f_312 * hk_556[k]
                  + f_314 * hk_560[k]
                  - f_312 * hk_569[k]
                  + f_313 * hk_571[k]
                  - f_314 * hk_573[k]
                  - f_251 * hk_614[k]
                  - f_251 * hk_619[k]
                  + f_315 * hk_621[k]
                  + f_251 * hk_628[k]
                  - f_316 * hk_632[k]
                  + f_251 * hk_641[k]
                  - f_315 * hk_643[k]
                  + f_316 * hk_645[k];
    }

#pragma omp simd aligned(hk_36, hk_39, hk_41, hk_46, hk_48, hk_50, hk_57, hk_59, hk_61, \
                         hk_216, hk_219, hk_221, hk_226, hk_228, hk_230, hk_237, hk_239, \
                         hk_241, hk_288, hk_291, hk_293, hk_298, hk_300, hk_302, hk_309, \
                         hk_311, hk_313, hk_540, hk_543, hk_545, hk_550, hk_552, hk_554, \
                         hk_561, hk_563, hk_565, hk_612, hk_615, hk_617, hk_622, hk_624, \
                         hk_626, hk_633, hk_635, hk_637 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_40[k] = -f_226 * hk_36[k]
                  + f_226 * hk_39[k]
                  + f_227 * hk_41[k]
                  + f_224 * hk_46[k]
                  - f_153 * hk_48[k]
                  - f_228 * hk_50[k]
                  + f_223 * hk_57[k]
                  - f_225 * hk_59[k]
                  + f_138 * hk_61[k]
                  - f_231 * hk_216[k]
                  + f_231 * hk_219[k]
                  + f_233 * hk_221[k]
                  + f_230 * hk_226[k]
                  - f_228 * hk_228[k]
                  - f_234 * hk_230[k]
                  + f_229 * hk_237[k]
                  - f_153 * hk_239[k]
                  + f_232 * hk_241[k]
                  + f_237 * hk_288[k]
                  - f_237 * hk_291[k]
                  - f_139 * hk_293[k]
                  - f_153 * hk_298[k]
                  + f_238 * hk_300[k]
                  + f_154 * hk_302[k]
                  - f_235 * hk_309[k]
                  + f_236 * hk_311[k]
                  - f_239 * hk_313[k]
                  + f_241 * hk_540[k]
                  - f_241 * hk_543[k]
                  - f_242 * hk_545[k]
                  - f_240 * hk_550[k]
                  + f_233 * hk_552[k]
                  + f_243 * hk_554[k]
                  - f_226 * hk_561[k]
                  + f_227 * hk_563[k]
                  - f_228 * hk_565[k]
                  - f_244 * hk_612[k]
                  + f_244 * hk_615[k]
                  + f_232 * hk_617[k]
                  + f_233 * hk_622[k]
                  - f_245 * hk_624[k]
                  - f_246 * hk_626[k]
                  + f_237 * hk_633[k]
                  - f_139 * hk_635[k]
                  + f_154 * hk_637[k];
    }

#pragma omp simd aligned(hk_38, hk_43, hk_45, hk_52, hk_54, hk_65, hk_67, hk_218, hk_223, \
                         hk_225, hk_232, hk_234, hk_245, hk_247, hk_290, hk_295, hk_297, \
                         hk_304, hk_306, hk_317, hk_319, hk_542, hk_547, hk_549, hk_556, \
                         hk_558, hk_569, hk_571, hk_614, hk_619, hk_621, hk_628, hk_630, \
                         hk_641, hk_643 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_41[k] = f_197 * hk_38[k]
                  - f_317 * hk_43[k]
                  - f_206 * hk_45[k]
                  - f_317 * hk_52[k]
                  + f_192 * hk_54[k]
                  + f_197 * hk_65[k]
                  - f_206 * hk_67[k]
                  + f_209 * hk_218[k]
                  - f_206 * hk_223[k]
                  - f_210 * hk_225[k]
                  - f_206 * hk_232[k]
                  + f_198 * hk_234[k]
                  + f_209 * hk_245[k]
                  - f_210 * hk_247[k]
                  - f_318 * hk_290[k]
                  + f_319 * hk_295[k]
                  + f_211 * hk_297[k]
                  + f_319 * hk_304[k]
                  - f_203 * hk_306[k]
                  - f_318 * hk_317[k]
                  + f_211 * hk_319[k]
                  - f_320 * hk_542[k]
                  + f_321 * hk_547[k]
                  + f_322 * hk_549[k]
                  + f_321 * hk_556[k]
                  - f_196 * hk_558[k]
                  - f_320 * hk_569[k]
                  + f_322 * hk_571[k]
                  + f_216 * hk_614[k]
                  - f_198 * hk_619[k]
                  - f_217 * hk_621[k]
                  - f_198 * hk_628[k]
                  + f_213 * hk_630[k]
                  + f_216 * hk_641[k]
                  - f_217 * hk_643[k];
    }

#pragma omp simd aligned(hk_36, hk_39, hk_41, hk_46, hk_48, hk_57, hk_59, hk_216, hk_219, \
                         hk_221, hk_226, hk_228, hk_237, hk_239, hk_288, hk_291, hk_293, \
                         hk_298, hk_300, hk_309, hk_311, hk_540, hk_543, hk_545, hk_550, \
                         hk_552, hk_561, hk_563, hk_612, hk_615, hk_617, hk_622, hk_624, \
                         hk_633, hk_635 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = f_193 * hk_36[k]
                  - f_191 * hk_39[k]
                  - f_194 * hk_41[k]
                  - f_189 * hk_46[k]
                  + f_192 * hk_48[k]
                  + f_189 * hk_57[k]
                  - f_190 * hk_59[k]
                  + f_199 * hk_216[k]
                  - f_197 * hk_219[k]
                  - f_200 * hk_221[k]
                  - f_195 * hk_226[k]
                  + f_198 * hk_228[k]
                  + f_195 * hk_237[k]
                  - f_196 * hk_239[k]
                  - f_200 * hk_288[k]
                  + f_202 * hk_291[k]
                  + f_204 * hk_293[k]
                  + f_196 * hk_298[k]
                  - f_203 * hk_300[k]
                  - f_196 * hk_309[k]
                  + f_201 * hk_311[k]
                  - f_208 * hk_540[k]
                  + f_207 * hk_543[k]
                  + f_209 * hk_545[k]
                  + f_205 * hk_550[k]
                  - f_196 * hk_552[k]
                  - f_205 * hk_561[k]
                  + f_206 * hk_563[k]
                  + f_214 * hk_612[k]
                  - f_212 * hk_615[k]
                  - f_215 * hk_617[k]
                  - f_210 * hk_622[k]
                  + f_213 * hk_624[k]
                  + f_210 * hk_633[k]
                  - f_211 * hk_635[k];
    }

#pragma omp simd aligned(hk_38, hk_43, hk_52, hk_65, hk_218, hk_223, hk_232, hk_245, hk_290, \
                         hk_295, hk_304, hk_317, hk_542, hk_547, hk_556, hk_569, hk_614, \
                         hk_619, hk_628, hk_641 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_43[k] = -f_323 * hk_38[k]
                  + f_324 * hk_43[k]
                  - f_324 * hk_52[k]
                  + f_323 * hk_65[k]
                  - f_325 * hk_218[k]
                  + f_326 * hk_223[k]
                  - f_326 * hk_232[k]
                  + f_325 * hk_245[k]
                  + f_327 * hk_290[k]
                  - f_328 * hk_295[k]
                  + f_328 * hk_304[k]
                  - f_327 * hk_317[k]
                  + f_329 * hk_542[k]
                  - f_330 * hk_547[k]
                  + f_330 * hk_556[k]
                  - f_329 * hk_569[k]
                  - f_331 * hk_614[k]
                  + f_332 * hk_619[k]
                  - f_332 * hk_628[k]
                  + f_331 * hk_641[k];
    }

#pragma omp simd aligned(hk_36, hk_39, hk_46, hk_57, hk_216, hk_219, hk_226, hk_237, hk_288, \
                         hk_291, hk_298, hk_309, hk_540, hk_543, hk_550, hk_561, hk_612, \
                         hk_615, hk_622, hk_633 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_44[k] = -f_164 * hk_36[k]
                  + f_163 * hk_39[k]
                  - f_162 * hk_46[k]
                  + f_161 * hk_57[k]
                  - f_168 * hk_216[k]
                  + f_167 * hk_219[k]
                  - f_166 * hk_226[k]
                  + f_165 * hk_237[k]
                  + f_172 * hk_288[k]
                  - f_171 * hk_291[k]
                  + f_170 * hk_298[k]
                  - f_169 * hk_309[k]
                  + f_175 * hk_540[k]
                  - f_161 * hk_543[k]
                  + f_174 * hk_550[k]
                  - f_173 * hk_561[k]
                  - f_178 * hk_612[k]
                  + f_169 * hk_615[k]
                  - f_177 * hk_622[k]
                  + f_176 * hk_633[k];
    }

#pragma omp simd aligned(hk_145, hk_150, hk_159, hk_172, hk_397, hk_402, hk_411, hk_424, \
                         hk_469, hk_474, hk_483, hk_496 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_45[k] = -f_333 * hk_145[k]
                  + f_334 * hk_150[k]
                  - f_335 * hk_159[k]
                  + f_336 * hk_172[k]
                  - f_333 * hk_397[k]
                  + f_334 * hk_402[k]
                  - f_335 * hk_411[k]
                  + f_336 * hk_424[k]
                  + f_337 * hk_469[k]
                  - f_338 * hk_474[k]
                  + f_339 * hk_483[k]
                  - f_340 * hk_496[k];
    }

#pragma omp simd aligned(hk_148, hk_155, hk_166, hk_400, hk_407, hk_418, hk_472, hk_479, \
                         hk_490 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_46[k] = -f_341 * hk_148[k]
                  + f_342 * hk_155[k]
                  - f_341 * hk_166[k]
                  - f_341 * hk_400[k]
                  + f_342 * hk_407[k]
                  - f_341 * hk_418[k]
                  + f_343 * hk_472[k]
                  - f_344 * hk_479[k]
                  + f_343 * hk_490[k];
    }

#pragma omp simd aligned(hk_145, hk_150, hk_152, hk_159, hk_161, hk_172, hk_174, hk_397, \
                         hk_402, hk_404, hk_411, hk_413, hk_424, hk_426, hk_469, hk_474, \
                         hk_476, hk_483, hk_485, hk_496, hk_498 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_47[k] = f_345 * hk_145[k]
                  - f_345 * hk_150[k]
                  - f_346 * hk_152[k]
                  - f_347 * hk_159[k]
                  + f_348 * hk_161[k]
                  + f_349 * hk_172[k]
                  - f_350 * hk_174[k]
                  + f_345 * hk_397[k]
                  - f_345 * hk_402[k]
                  - f_346 * hk_404[k]
                  - f_347 * hk_411[k]
                  + f_348 * hk_413[k]
                  + f_349 * hk_424[k]
                  - f_350 * hk_426[k]
                  - f_351 * hk_469[k]
                  + f_351 * hk_474[k]
                  + f_348 * hk_476[k]
                  + f_352 * hk_483[k]
                  - f_353 * hk_485[k]
                  - f_354 * hk_496[k]
                  + f_355 * hk_498[k];
    }

#pragma omp simd aligned(hk_148, hk_157, hk_166, hk_168, hk_400, hk_409, hk_418, hk_420, \
                         hk_472, hk_481, hk_490, hk_492 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_48[k] = f_355 * hk_148[k]
                  - f_356 * hk_157[k]
                  - f_355 * hk_166[k]
                  + f_356 * hk_168[k]
                  + f_355 * hk_400[k]
                  - f_356 * hk_409[k]
                  - f_355 * hk_418[k]
                  + f_356 * hk_420[k]
                  - f_357 * hk_472[k]
                  + f_358 * hk_481[k]
                  + f_357 * hk_490[k]
                  - f_358 * hk_492[k];
    }

#pragma omp simd aligned(hk_145, hk_150, hk_152, hk_159, hk_161, hk_163, hk_172, hk_174, \
                         hk_176, hk_397, hk_402, hk_404, hk_411, hk_413, hk_415, hk_424, \
                         hk_426, hk_428, hk_469, hk_474, hk_476, hk_483, hk_485, hk_487, \
                         hk_496, hk_498, hk_500 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_49[k] = -f_359 * hk_145[k]
                  - f_143 * hk_150[k]
                  + f_360 * hk_152[k]
                  - f_361 * hk_159[k]
                  + f_144 * hk_161[k]
                  - f_145 * hk_163[k]
                  + f_361 * hk_172[k]
                  - f_362 * hk_174[k]
                  + f_363 * hk_176[k]
                  - f_359 * hk_397[k]
                  - f_143 * hk_402[k]
                  + f_360 * hk_404[k]
                  - f_361 * hk_411[k]
                  + f_144 * hk_413[k]
                  - f_145 * hk_415[k]
                  + f_361 * hk_424[k]
                  - f_362 * hk_426[k]
                  + f_363 * hk_428[k]
                  + f_364 * hk_469[k]
                  + f_365 * hk_474[k]
                  - f_366 * hk_476[k]
                  + f_367 * hk_483[k]
                  - f_145 * hk_485[k]
                  + f_368 * hk_487[k]
                  - f_367 * hk_496[k]
                  + f_144 * hk_498[k]
                  - f_369 * hk_500[k];
    }

#pragma omp simd aligned(hk_148, hk_155, hk_157, hk_166, hk_168, hk_170, hk_400, hk_407, \
                         hk_409, hk_418, hk_420, hk_422, hk_472, hk_479, hk_481, hk_490, \
                         hk_492, hk_494 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_50[k] = -f_265 * hk_148[k]
                  - f_370 * hk_155[k]
                  + f_371 * hk_157[k]
                  - f_265 * hk_166[k]
                  + f_371 * hk_168[k]
                  - f_372 * hk_170[k]
                  - f_265 * hk_400[k]
                  - f_370 * hk_407[k]
                  + f_371 * hk_409[k]
                  - f_265 * hk_418[k]
                  + f_371 * hk_420[k]
                  - f_372 * hk_422[k]
                  + f_370 * hk_472[k]
                  + f_273 * hk_479[k]
                  - f_373 * hk_481[k]
                  + f_370 * hk_490[k]
                  - f_373 * hk_492[k]
                  + f_374 * hk_494[k];
    }

#pragma omp simd aligned(hk_145, hk_150, hk_152, hk_159, hk_161, hk_163, hk_172, hk_174, \
                         hk_176, hk_178, hk_397, hk_402, hk_404, hk_411, hk_413, hk_415, \
                         hk_424, hk_426, hk_428, hk_430, hk_469, hk_474, hk_476, hk_483, \
                         hk_485, hk_487, hk_496, hk_498, hk_500, \
                         hk_502 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_51[k] = f_375 * hk_145[k]
                  + f_250 * hk_150[k]
                  - f_248 * hk_152[k]
                  + f_250 * hk_159[k]
                  - f_137 * hk_161[k]
                  + f_137 * hk_163[k]
                  + f_375 * hk_172[k]
                  - f_248 * hk_174[k]
                  + f_137 * hk_176[k]
                  - f_376 * hk_178[k]
                  + f_375 * hk_397[k]
                  + f_250 * hk_402[k]
                  - f_248 * hk_404[k]
                  + f_250 * hk_411[k]
                  - f_137 * hk_413[k]
                  + f_137 * hk_415[k]
                  + f_375 * hk_424[k]
                  - f_248 * hk_426[k]
                  + f_137 * hk_428[k]
                  - f_376 * hk_430[k]
                  - f_377 * hk_469[k]
                  - f_251 * hk_474[k]
                  + f_137 * hk_476[k]
                  - f_251 * hk_483[k]
                  + f_310 * hk_485[k]
                  - f_310 * hk_487[k]
                  - f_377 * hk_496[k]
                  + f_137 * hk_498[k]
                  - f_310 * hk_500[k]
                  + f_378 * hk_502[k];
    }

#pragma omp simd aligned(hk_146, hk_151, hk_153, hk_160, hk_162, hk_164, hk_173, hk_175, \
                         hk_177, hk_179, hk_398, hk_403, hk_405, hk_412, hk_414, hk_416, \
                         hk_425, hk_427, hk_429, hk_431, hk_470, hk_475, hk_477, hk_484, \
                         hk_486, hk_488, hk_497, hk_499, hk_501, \
                         hk_503 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_52[k] = f_379 * hk_146[k]
                  + f_380 * hk_151[k]
                  - f_381 * hk_153[k]
                  + f_380 * hk_160[k]
                  - f_382 * hk_162[k]
                  + f_383 * hk_164[k]
                  + f_379 * hk_173[k]
                  - f_381 * hk_175[k]
                  + f_383 * hk_177[k]
                  - f_384 * hk_179[k]
                  + f_379 * hk_398[k]
                  + f_380 * hk_403[k]
                  - f_381 * hk_405[k]
                  + f_380 * hk_412[k]
                  - f_382 * hk_414[k]
                  + f_383 * hk_416[k]
                  + f_379 * hk_425[k]
                  - f_381 * hk_427[k]
                  + f_383 * hk_429[k]
                  - f_384 * hk_431[k]
                  - f_385 * hk_470[k]
                  - f_381 * hk_475[k]
                  + f_382 * hk_477[k]
                  - f_381 * hk_484[k]
                  + f_386 * hk_486[k]
                  - f_387 * hk_488[k]
                  - f_385 * hk_497[k]
                  + f_382 * hk_499[k]
                  - f_387 * hk_501[k]
                  + f_388 * hk_503[k];
    }

#pragma omp simd aligned(hk_144, hk_147, hk_149, hk_154, hk_156, hk_158, hk_165, hk_167, \
                         hk_169, hk_171, hk_396, hk_399, hk_401, hk_406, hk_408, hk_410, \
                         hk_417, hk_419, hk_421, hk_423, hk_468, hk_471, hk_473, hk_478, \
                         hk_480, hk_482, hk_489, hk_491, hk_493, \
                         hk_495 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_53[k] = f_375 * hk_144[k]
                  + f_250 * hk_147[k]
                  - f_248 * hk_149[k]
                  + f_250 * hk_154[k]
                  - f_137 * hk_156[k]
                  + f_137 * hk_158[k]
                  + f_375 * hk_165[k]
                  - f_248 * hk_167[k]
                  + f_137 * hk_169[k]
                  - f_376 * hk_171[k]
                  + f_375 * hk_396[k]
                  + f_250 * hk_399[k]
                  - f_248 * hk_401[k]
                  + f_250 * hk_406[k]
                  - f_137 * hk_408[k]
                  + f_137 * hk_410[k]
                  + f_375 * hk_417[k]
                  - f_248 * hk_419[k]
                  + f_137 * hk_421[k]
                  - f_376 * hk_423[k]
                  - f_377 * hk_468[k]
                  - f_251 * hk_471[k]
                  + f_137 * hk_473[k]
                  - f_251 * hk_478[k]
                  + f_310 * hk_480[k]
                  - f_310 * hk_482[k]
                  - f_377 * hk_489[k]
                  + f_137 * hk_491[k]
                  - f_310 * hk_493[k]
                  + f_378 * hk_495[k];
    }

#pragma omp simd aligned(hk_146, hk_151, hk_153, hk_160, hk_164, hk_173, hk_175, hk_177, \
                         hk_398, hk_403, hk_405, hk_412, hk_416, hk_425, hk_427, hk_429, \
                         hk_470, hk_475, hk_477, hk_484, hk_488, hk_497, hk_499, \
                         hk_501 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_54[k] = -f_264 * hk_146[k]
                  - f_264 * hk_151[k]
                  + f_280 * hk_153[k]
                  + f_264 * hk_160[k]
                  - f_389 * hk_164[k]
                  + f_264 * hk_173[k]
                  - f_280 * hk_175[k]
                  + f_389 * hk_177[k]
                  - f_264 * hk_398[k]
                  - f_264 * hk_403[k]
                  + f_280 * hk_405[k]
                  + f_264 * hk_412[k]
                  - f_389 * hk_416[k]
                  + f_264 * hk_425[k]
                  - f_280 * hk_427[k]
                  + f_389 * hk_429[k]
                  + f_265 * hk_470[k]
                  + f_265 * hk_475[k]
                  - f_371 * hk_477[k]
                  - f_265 * hk_484[k]
                  + f_372 * hk_488[k]
                  - f_265 * hk_497[k]
                  + f_371 * hk_499[k]
                  - f_372 * hk_501[k];
    }

#pragma omp simd aligned(hk_144, hk_147, hk_149, hk_154, hk_156, hk_158, hk_165, hk_167, \
                         hk_169, hk_396, hk_399, hk_401, hk_406, hk_408, hk_410, hk_417, \
                         hk_419, hk_421, hk_468, hk_471, hk_473, hk_478, hk_480, hk_482, \
                         hk_489, hk_491, hk_493 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_55[k] = -f_361 * hk_144[k]
                  + f_361 * hk_147[k]
                  + f_362 * hk_149[k]
                  + f_143 * hk_154[k]
                  - f_144 * hk_156[k]
                  - f_363 * hk_158[k]
                  + f_359 * hk_165[k]
                  - f_360 * hk_167[k]
                  + f_145 * hk_169[k]
                  - f_361 * hk_396[k]
                  + f_361 * hk_399[k]
                  + f_362 * hk_401[k]
                  + f_143 * hk_406[k]
                  - f_144 * hk_408[k]
                  - f_363 * hk_410[k]
                  + f_359 * hk_417[k]
                  - f_360 * hk_419[k]
                  + f_145 * hk_421[k]
                  + f_367 * hk_468[k]
                  - f_367 * hk_471[k]
                  - f_144 * hk_473[k]
                  - f_365 * hk_478[k]
                  + f_145 * hk_480[k]
                  + f_369 * hk_482[k]
                  - f_364 * hk_489[k]
                  + f_366 * hk_491[k]
                  - f_368 * hk_493[k];
    }

#pragma omp simd aligned(hk_146, hk_151, hk_153, hk_160, hk_162, hk_173, hk_175, hk_398, \
                         hk_403, hk_405, hk_412, hk_414, hk_425, hk_427, hk_470, hk_475, \
                         hk_477, hk_484, hk_486, hk_497, hk_499 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_56[k] = f_390 * hk_146[k]
                  - f_391 * hk_151[k]
                  - f_392 * hk_153[k]
                  - f_391 * hk_160[k]
                  + f_348 * hk_162[k]
                  + f_390 * hk_173[k]
                  - f_392 * hk_175[k]
                  + f_390 * hk_398[k]
                  - f_391 * hk_403[k]
                  - f_392 * hk_405[k]
                  - f_391 * hk_412[k]
                  + f_348 * hk_414[k]
                  + f_390 * hk_425[k]
                  - f_392 * hk_427[k]
                  - f_350 * hk_470[k]
                  + f_346 * hk_475[k]
                  + f_393 * hk_477[k]
                  + f_346 * hk_484[k]
                  - f_353 * hk_486[k]
                  - f_350 * hk_497[k]
                  + f_393 * hk_499[k];
    }

#pragma omp simd aligned(hk_144, hk_147, hk_149, hk_154, hk_156, hk_165, hk_167, hk_396, \
                         hk_399, hk_401, hk_406, hk_408, hk_417, hk_419, hk_468, hk_471, \
                         hk_473, hk_478, hk_480, hk_489, hk_491 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_57[k] = f_349 * hk_144[k]
                  - f_347 * hk_147[k]
                  - f_350 * hk_149[k]
                  - f_345 * hk_154[k]
                  + f_348 * hk_156[k]
                  + f_345 * hk_165[k]
                  - f_346 * hk_167[k]
                  + f_349 * hk_396[k]
                  - f_347 * hk_399[k]
                  - f_350 * hk_401[k]
                  - f_345 * hk_406[k]
                  + f_348 * hk_408[k]
                  + f_345 * hk_417[k]
                  - f_346 * hk_419[k]
                  - f_354 * hk_468[k]
                  + f_352 * hk_471[k]
                  + f_355 * hk_473[k]
                  + f_351 * hk_478[k]
                  - f_353 * hk_480[k]
                  - f_351 * hk_489[k]
                  + f_348 * hk_491[k];
    }

#pragma omp simd aligned(hk_146, hk_151, hk_160, hk_173, hk_398, hk_403, hk_412, hk_425, \
                         hk_470, hk_475, hk_484, hk_497 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_58[k] = -f_394 * hk_146[k]
                  + f_395 * hk_151[k]
                  - f_395 * hk_160[k]
                  + f_394 * hk_173[k]
                  - f_394 * hk_398[k]
                  + f_395 * hk_403[k]
                  - f_395 * hk_412[k]
                  + f_394 * hk_425[k]
                  + f_396 * hk_470[k]
                  - f_397 * hk_475[k]
                  + f_397 * hk_484[k]
                  - f_396 * hk_497[k];
    }

#pragma omp simd aligned(hk_144, hk_147, hk_154, hk_165, hk_396, hk_399, hk_406, hk_417, \
                         hk_468, hk_471, hk_478, hk_489 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_59[k] = -f_336 * hk_144[k]
                  + f_335 * hk_147[k]
                  - f_334 * hk_154[k]
                  + f_333 * hk_165[k]
                  - f_336 * hk_396[k]
                  + f_335 * hk_399[k]
                  - f_334 * hk_406[k]
                  + f_333 * hk_417[k]
                  + f_340 * hk_468[k]
                  - f_339 * hk_471[k]
                  + f_338 * hk_478[k]
                  - f_337 * hk_489[k];
    }

#pragma omp simd aligned(hk_37, hk_42, hk_51, hk_64, hk_217, hk_222, hk_231, hk_244, hk_289, \
                         hk_294, hk_303, hk_316, hk_541, hk_546, hk_555, hk_568, hk_613, \
                         hk_618, hk_627, hk_640, hk_685, hk_690, hk_699, \
                         hk_712 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_60[k] = f_398 * hk_37[k]
                  - f_399 * hk_42[k]
                  + f_400 * hk_51[k]
                  - f_401 * hk_64[k]
                  + f_402 * hk_217[k]
                  - f_403 * hk_222[k]
                  + f_404 * hk_231[k]
                  - f_405 * hk_244[k]
                  - f_406 * hk_289[k]
                  + f_407 * hk_294[k]
                  - f_408 * hk_303[k]
                  + f_409 * hk_316[k]
                  + f_398 * hk_541[k]
                  - f_399 * hk_546[k]
                  + f_400 * hk_555[k]
                  - f_401 * hk_568[k]
                  - f_406 * hk_613[k]
                  + f_407 * hk_618[k]
                  - f_408 * hk_627[k]
                  + f_409 * hk_640[k]
                  + f_410 * hk_685[k]
                  - f_411 * hk_690[k]
                  + f_412 * hk_699[k]
                  - f_413 * hk_712[k];
    }

#pragma omp simd aligned(hk_40, hk_47, hk_58, hk_220, hk_227, hk_238, hk_292, hk_299, hk_310, \
                         hk_544, hk_551, hk_562, hk_616, hk_623, hk_634, hk_688, hk_695, \
                         hk_706 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_61[k] = f_414 * hk_40[k]
                  - f_415 * hk_47[k]
                  + f_414 * hk_58[k]
                  + f_416 * hk_220[k]
                  - f_417 * hk_227[k]
                  + f_416 * hk_238[k]
                  - f_418 * hk_292[k]
                  + f_419 * hk_299[k]
                  - f_418 * hk_310[k]
                  + f_414 * hk_544[k]
                  - f_415 * hk_551[k]
                  + f_414 * hk_562[k]
                  - f_418 * hk_616[k]
                  + f_419 * hk_623[k]
                  - f_418 * hk_634[k]
                  + f_420 * hk_688[k]
                  - f_421 * hk_695[k]
                  + f_420 * hk_706[k];
    }

#pragma omp simd aligned(hk_37, hk_42, hk_44, hk_51, hk_53, hk_64, hk_66, hk_217, hk_222, \
                         hk_224, hk_231, hk_233, hk_244, hk_246, hk_289, hk_294, hk_296, \
                         hk_303, hk_305, hk_316, hk_318, hk_541, hk_546, hk_548, hk_555, \
                         hk_557, hk_568, hk_570, hk_613, hk_618, hk_620, hk_627, hk_629, \
                         hk_640, hk_642, hk_685, hk_690, hk_692, hk_699, hk_701, hk_712, \
                         hk_714 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_62[k] = -f_422 * hk_37[k]
                  + f_422 * hk_42[k]
                  + f_423 * hk_44[k]
                  + f_424 * hk_51[k]
                  - f_425 * hk_53[k]
                  - f_426 * hk_64[k]
                  + f_427 * hk_66[k]
                  - f_428 * hk_217[k]
                  + f_428 * hk_222[k]
                  + f_425 * hk_224[k]
                  + f_429 * hk_231[k]
                  - f_430 * hk_233[k]
                  - f_431 * hk_244[k]
                  + f_432 * hk_246[k]
                  + f_423 * hk_289[k]
                  - f_423 * hk_294[k]
                  - f_433 * hk_296[k]
                  - f_434 * hk_303[k]
                  + f_435 * hk_305[k]
                  + f_427 * hk_316[k]
                  - f_436 * hk_318[k]
                  - f_422 * hk_541[k]
                  + f_422 * hk_546[k]
                  + f_423 * hk_548[k]
                  + f_424 * hk_555[k]
                  - f_425 * hk_557[k]
                  - f_426 * hk_568[k]
                  + f_427 * hk_570[k]
                  + f_423 * hk_613[k]
                  - f_423 * hk_618[k]
                  - f_433 * hk_620[k]
                  - f_434 * hk_627[k]
                  + f_435 * hk_629[k]
                  + f_427 * hk_640[k]
                  - f_436 * hk_642[k]
                  - f_437 * hk_685[k]
                  + f_437 * hk_690[k]
                  + f_438 * hk_692[k]
                  + f_439 * hk_699[k]
                  - f_440 * hk_701[k]
                  - f_441 * hk_712[k]
                  + f_442 * hk_714[k];
    }

#pragma omp simd aligned(hk_40, hk_49, hk_58, hk_60, hk_220, hk_229, hk_238, hk_240, hk_292, \
                         hk_301, hk_310, hk_312, hk_544, hk_553, hk_562, hk_564, hk_616, \
                         hk_625, hk_634, hk_636, hk_688, hk_697, hk_706, \
                         hk_708 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_63[k] = -f_432 * hk_40[k]
                  + f_443 * hk_49[k]
                  + f_432 * hk_58[k]
                  - f_443 * hk_60[k]
                  - f_444 * hk_220[k]
                  + f_445 * hk_229[k]
                  + f_444 * hk_238[k]
                  - f_445 * hk_240[k]
                  + f_446 * hk_292[k]
                  - f_440 * hk_301[k]
                  - f_446 * hk_310[k]
                  + f_440 * hk_312[k]
                  - f_432 * hk_544[k]
                  + f_443 * hk_553[k]
                  + f_432 * hk_562[k]
                  - f_443 * hk_564[k]
                  + f_446 * hk_616[k]
                  - f_440 * hk_625[k]
                  - f_446 * hk_634[k]
                  + f_440 * hk_636[k]
                  - f_447 * hk_688[k]
                  + f_448 * hk_697[k]
                  + f_447 * hk_706[k]
                  - f_448 * hk_708[k];
    }

#pragma omp simd aligned(hk_37, hk_42, hk_44, hk_51, hk_53, hk_55, hk_64, hk_66, hk_68, \
                         hk_217, hk_222, hk_224, hk_231, hk_233, hk_235, hk_244, hk_246, \
                         hk_248, hk_289, hk_294, hk_296, hk_303, hk_305, hk_307, hk_316, \
                         hk_318, hk_320, hk_541, hk_546, hk_548, hk_555, hk_557, hk_559, \
                         hk_568, hk_570, hk_572, hk_613, hk_618, hk_620, hk_627, hk_629, \
                         hk_631, hk_640, hk_642, hk_644, hk_685, hk_690, hk_692, hk_699, \
                         hk_701, hk_703, hk_712, hk_714, hk_716 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_64[k] = f_449 * hk_37[k]
                  + f_450 * hk_42[k]
                  - f_451 * hk_44[k]
                  + f_452 * hk_51[k]
                  - f_453 * hk_53[k]
                  + f_454 * hk_55[k]
                  - f_452 * hk_64[k]
                  + f_455 * hk_66[k]
                  - f_456 * hk_68[k]
                  + f_457 * hk_217[k]
                  + f_458 * hk_222[k]
                  - f_459 * hk_224[k]
                  + f_460 * hk_231[k]
                  - f_454 * hk_233[k]
                  + f_461 * hk_235[k]
                  - f_460 * hk_244[k]
                  + f_453 * hk_246[k]
                  - f_462 * hk_248[k]
                  - f_463 * hk_289[k]
                  - f_451 * hk_294[k]
                  + f_464 * hk_296[k]
                  - f_465 * hk_303[k]
                  + f_466 * hk_305[k]
                  - f_467 * hk_307[k]
                  + f_465 * hk_316[k]
                  - f_468 * hk_318[k]
                  + f_469 * hk_320[k]
                  + f_449 * hk_541[k]
                  + f_450 * hk_546[k]
                  - f_451 * hk_548[k]
                  + f_452 * hk_555[k]
                  - f_453 * hk_557[k]
                  + f_454 * hk_559[k]
                  - f_452 * hk_568[k]
                  + f_455 * hk_570[k]
                  - f_456 * hk_572[k]
                  - f_463 * hk_613[k]
                  - f_451 * hk_618[k]
                  + f_464 * hk_620[k]
                  - f_465 * hk_627[k]
                  + f_466 * hk_629[k]
                  - f_467 * hk_631[k]
                  + f_465 * hk_640[k]
                  - f_468 * hk_642[k]
                  + f_469 * hk_644[k]
                  + f_470 * hk_685[k]
                  + f_453 * hk_690[k]
                  - f_466 * hk_692[k]
                  + f_471 * hk_699[k]
                  - f_469 * hk_701[k]
                  + f_472 * hk_703[k]
                  - f_471 * hk_712[k]
                  + f_461 * hk_714[k]
                  - f_473 * hk_716[k];
    }

#pragma omp simd aligned(hk_40, hk_47, hk_49, hk_58, hk_60, hk_62, hk_220, hk_227, hk_229, \
                         hk_238, hk_240, hk_242, hk_292, hk_299, hk_301, hk_310, hk_312, \
                         hk_314, hk_544, hk_551, hk_553, hk_562, hk_564, hk_566, hk_616, \
                         hk_623, hk_625, hk_634, hk_636, hk_638, hk_688, hk_695, hk_697, \
                         hk_706, hk_708, hk_710 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_65[k] = f_474 * hk_40[k]
                  + f_475 * hk_47[k]
                  - f_476 * hk_49[k]
                  + f_474 * hk_58[k]
                  - f_476 * hk_60[k]
                  + f_477 * hk_62[k]
                  + f_475 * hk_220[k]
                  + f_478 * hk_227[k]
                  - f_479 * hk_229[k]
                  + f_475 * hk_238[k]
                  - f_479 * hk_240[k]
                  + f_480 * hk_242[k]
                  - f_481 * hk_292[k]
                  - f_482 * hk_299[k]
                  + f_483 * hk_301[k]
                  - f_481 * hk_310[k]
                  + f_483 * hk_312[k]
                  - f_484 * hk_314[k]
                  + f_474 * hk_544[k]
                  + f_475 * hk_551[k]
                  - f_476 * hk_553[k]
                  + f_474 * hk_562[k]
                  - f_476 * hk_564[k]
                  + f_477 * hk_566[k]
                  - f_481 * hk_616[k]
                  - f_482 * hk_623[k]
                  + f_483 * hk_625[k]
                  - f_481 * hk_634[k]
                  + f_483 * hk_636[k]
                  - f_484 * hk_638[k]
                  + f_485 * hk_688[k]
                  + f_486 * hk_695[k]
                  - f_487 * hk_697[k]
                  + f_485 * hk_706[k]
                  - f_487 * hk_708[k]
                  + f_488 * hk_710[k];
    }

#pragma omp simd aligned(hk_37, hk_42, hk_44, hk_51, hk_53, hk_55, hk_64, hk_66, hk_68, hk_70, \
                         hk_217, hk_222, hk_224, hk_231, hk_233, hk_235, hk_244, hk_246, \
                         hk_248, hk_250, hk_289, hk_294, hk_296, hk_303, hk_305, hk_307, \
                         hk_316, hk_318, hk_320, hk_322, hk_541, hk_546, hk_548, hk_555, \
                         hk_557, hk_559, hk_568, hk_570, hk_572, hk_574, hk_613, hk_618, \
                         hk_620, hk_627, hk_629, hk_631, hk_640, hk_642, hk_644, hk_646, \
                         hk_685, hk_690, hk_692, hk_699, hk_701, hk_703, hk_712, hk_714, \
                         hk_716, hk_718 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_66[k] = -f_489 * hk_37[k]
                  - f_490 * hk_42[k]
                  + f_491 * hk_44[k]
                  - f_490 * hk_51[k]
                  + f_492 * hk_53[k]
                  - f_492 * hk_55[k]
                  - f_489 * hk_64[k]
                  + f_491 * hk_66[k]
                  - f_492 * hk_68[k]
                  + f_493 * hk_70[k]
                  - f_494 * hk_217[k]
                  - f_495 * hk_222[k]
                  + f_492 * hk_224[k]
                  - f_495 * hk_231[k]
                  + f_496 * hk_233[k]
                  - f_496 * hk_235[k]
                  - f_494 * hk_244[k]
                  + f_492 * hk_246[k]
                  - f_496 * hk_248[k]
                  + f_384 * hk_250[k]
                  + f_497 * hk_289[k]
                  + f_498 * hk_294[k]
                  - f_499 * hk_296[k]
                  + f_498 * hk_303[k]
                  - f_500 * hk_305[k]
                  + f_500 * hk_307[k]
                  + f_497 * hk_316[k]
                  - f_499 * hk_318[k]
                  + f_500 * hk_320[k]
                  - f_501 * hk_322[k]
                  - f_489 * hk_541[k]
                  - f_490 * hk_546[k]
                  + f_491 * hk_548[k]
                  - f_490 * hk_555[k]
                  + f_492 * hk_557[k]
                  - f_492 * hk_559[k]
                  - f_489 * hk_568[k]
                  + f_491 * hk_570[k]
                  - f_492 * hk_572[k]
                  + f_493 * hk_574[k]
                  + f_497 * hk_613[k]
                  + f_498 * hk_618[k]
                  - f_499 * hk_620[k]
                  + f_498 * hk_627[k]
                  - f_500 * hk_629[k]
                  + f_500 * hk_631[k]
                  + f_497 * hk_640[k]
                  - f_499 * hk_642[k]
                  + f_500 * hk_644[k]
                  - f_501 * hk_646[k]
                  - f_502 * hk_685[k]
                  - f_491 * hk_690[k]
                  + f_503 * hk_692[k]
                  - f_491 * hk_699[k]
                  + f_504 * hk_701[k]
                  - f_504 * hk_703[k]
                  - f_502 * hk_712[k]
                  + f_503 * hk_714[k]
                  - f_504 * hk_716[k]
                  + f_505 * hk_718[k];
    }

#pragma omp simd aligned(hk_38, hk_43, hk_45, hk_52, hk_54, hk_56, hk_65, hk_67, hk_69, hk_71, \
                         hk_218, hk_223, hk_225, hk_232, hk_234, hk_236, hk_245, hk_247, \
                         hk_249, hk_251, hk_290, hk_295, hk_297, hk_304, hk_306, hk_308, \
                         hk_317, hk_319, hk_321, hk_323, hk_542, hk_547, hk_549, hk_556, \
                         hk_558, hk_560, hk_569, hk_571, hk_573, hk_575, hk_614, hk_619, \
                         hk_621, hk_628, hk_630, hk_632, hk_641, hk_643, hk_645, hk_647, \
                         hk_686, hk_691, hk_693, hk_700, hk_702, hk_704, hk_713, hk_715, \
                         hk_717, hk_719 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_67[k] = -f_506 * hk_38[k]
                  - f_256 * hk_43[k]
                  + f_250 * hk_45[k]
                  - f_256 * hk_52[k]
                  + f_251 * hk_54[k]
                  - f_314 * hk_56[k]
                  - f_506 * hk_65[k]
                  + f_250 * hk_67[k]
                  - f_314 * hk_69[k]
                  + f_507 * hk_71[k]
                  - f_375 * hk_218[k]
                  - f_250 * hk_223[k]
                  + f_251 * hk_225[k]
                  - f_250 * hk_232[k]
                  + f_259 * hk_234[k]
                  - f_258 * hk_236[k]
                  - f_375 * hk_245[k]
                  + f_251 * hk_247[k]
                  - f_258 * hk_249[k]
                  + f_508 * hk_251[k]
                  + f_251 * hk_290[k]
                  + f_309 * hk_295[k]
                  - f_136 * hk_297[k]
                  + f_309 * hk_304[k]
                  - f_134 * hk_306[k]
                  + f_509 * hk_308[k]
                  + f_251 * hk_317[k]
                  - f_136 * hk_319[k]
                  + f_509 * hk_321[k]
                  - f_510 * hk_323[k]
                  - f_506 * hk_542[k]
                  - f_256 * hk_547[k]
                  + f_250 * hk_549[k]
                  - f_256 * hk_556[k]
                  + f_251 * hk_558[k]
                  - f_314 * hk_560[k]
                  - f_506 * hk_569[k]
                  + f_250 * hk_571[k]
                  - f_314 * hk_573[k]
                  + f_507 * hk_575[k]
                  + f_251 * hk_614[k]
                  + f_309 * hk_619[k]
                  - f_136 * hk_621[k]
                  + f_309 * hk_628[k]
                  - f_134 * hk_630[k]
                  + f_509 * hk_632[k]
                  + f_251 * hk_641[k]
                  - f_136 * hk_643[k]
                  + f_509 * hk_645[k]
                  - f_510 * hk_647[k]
                  - f_313 * hk_686[k]
                  - f_259 * hk_691[k]
                  + f_248 * hk_693[k]
                  - f_259 * hk_700[k]
                  + f_137 * hk_702[k]
                  - f_316 * hk_704[k]
                  - f_313 * hk_713[k]
                  + f_248 * hk_715[k]
                  - f_316 * hk_717[k]
                  + f_511 * hk_719[k];
    }

#pragma omp simd aligned(hk_36, hk_39, hk_41, hk_46, hk_48, hk_50, hk_57, hk_59, hk_61, hk_63, \
                         hk_216, hk_219, hk_221, hk_226, hk_228, hk_230, hk_237, hk_239, \
                         hk_241, hk_243, hk_288, hk_291, hk_293, hk_298, hk_300, hk_302, \
                         hk_309, hk_311, hk_313, hk_315, hk_540, hk_543, hk_545, hk_550, \
                         hk_552, hk_554, hk_561, hk_563, hk_565, hk_567, hk_612, hk_615, \
                         hk_617, hk_622, hk_624, hk_626, hk_633, hk_635, hk_637, hk_639, \
                         hk_684, hk_687, hk_689, hk_694, hk_696, hk_698, hk_705, hk_707, \
                         hk_709, hk_711 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_68[k] = -f_489 * hk_36[k]
                  - f_490 * hk_39[k]
                  + f_491 * hk_41[k]
                  - f_490 * hk_46[k]
                  + f_492 * hk_48[k]
                  - f_492 * hk_50[k]
                  - f_489 * hk_57[k]
                  + f_491 * hk_59[k]
                  - f_492 * hk_61[k]
                  + f_493 * hk_63[k]
                  - f_494 * hk_216[k]
                  - f_495 * hk_219[k]
                  + f_492 * hk_221[k]
                  - f_495 * hk_226[k]
                  + f_496 * hk_228[k]
                  - f_496 * hk_230[k]
                  - f_494 * hk_237[k]
                  + f_492 * hk_239[k]
                  - f_496 * hk_241[k]
                  + f_384 * hk_243[k]
                  + f_497 * hk_288[k]
                  + f_498 * hk_291[k]
                  - f_499 * hk_293[k]
                  + f_498 * hk_298[k]
                  - f_500 * hk_300[k]
                  + f_500 * hk_302[k]
                  + f_497 * hk_309[k]
                  - f_499 * hk_311[k]
                  + f_500 * hk_313[k]
                  - f_501 * hk_315[k]
                  - f_489 * hk_540[k]
                  - f_490 * hk_543[k]
                  + f_491 * hk_545[k]
                  - f_490 * hk_550[k]
                  + f_492 * hk_552[k]
                  - f_492 * hk_554[k]
                  - f_489 * hk_561[k]
                  + f_491 * hk_563[k]
                  - f_492 * hk_565[k]
                  + f_493 * hk_567[k]
                  + f_497 * hk_612[k]
                  + f_498 * hk_615[k]
                  - f_499 * hk_617[k]
                  + f_498 * hk_622[k]
                  - f_500 * hk_624[k]
                  + f_500 * hk_626[k]
                  + f_497 * hk_633[k]
                  - f_499 * hk_635[k]
                  + f_500 * hk_637[k]
                  - f_501 * hk_639[k]
                  - f_502 * hk_684[k]
                  - f_491 * hk_687[k]
                  + f_503 * hk_689[k]
                  - f_491 * hk_694[k]
                  + f_504 * hk_696[k]
                  - f_504 * hk_698[k]
                  - f_502 * hk_705[k]
                  + f_503 * hk_707[k]
                  - f_504 * hk_709[k]
                  + f_505 * hk_711[k];
    }

#pragma omp simd aligned(hk_38, hk_43, hk_45, hk_52, hk_56, hk_65, hk_67, hk_69, hk_218, \
                         hk_223, hk_225, hk_232, hk_236, hk_245, hk_247, hk_249, hk_290, \
                         hk_295, hk_297, hk_304, hk_308, hk_317, hk_319, hk_321, hk_542, \
                         hk_547, hk_549, hk_556, hk_560, hk_569, hk_571, hk_573, hk_614, \
                         hk_619, hk_621, hk_628, hk_632, hk_641, hk_643, hk_645, hk_686, \
                         hk_691, hk_693, hk_700, hk_704, hk_713, hk_715, \
                         hk_717 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_69[k] = f_512 * hk_38[k]
                  + f_512 * hk_43[k]
                  - f_513 * hk_45[k]
                  - f_512 * hk_52[k]
                  + f_514 * hk_56[k]
                  - f_512 * hk_65[k]
                  + f_513 * hk_67[k]
                  - f_514 * hk_69[k]
                  + f_474 * hk_218[k]
                  + f_474 * hk_223[k]
                  - f_476 * hk_225[k]
                  - f_474 * hk_232[k]
                  + f_477 * hk_236[k]
                  - f_474 * hk_245[k]
                  + f_476 * hk_247[k]
                  - f_477 * hk_249[k]
                  - f_515 * hk_290[k]
                  - f_515 * hk_295[k]
                  + f_516 * hk_297[k]
                  + f_515 * hk_304[k]
                  - f_517 * hk_308[k]
                  + f_515 * hk_317[k]
                  - f_516 * hk_319[k]
                  + f_517 * hk_321[k]
                  + f_512 * hk_542[k]
                  + f_512 * hk_547[k]
                  - f_513 * hk_549[k]
                  - f_512 * hk_556[k]
                  + f_514 * hk_560[k]
                  - f_512 * hk_569[k]
                  + f_513 * hk_571[k]
                  - f_514 * hk_573[k]
                  - f_515 * hk_614[k]
                  - f_515 * hk_619[k]
                  + f_516 * hk_621[k]
                  + f_515 * hk_628[k]
                  - f_517 * hk_632[k]
                  + f_515 * hk_641[k]
                  - f_516 * hk_643[k]
                  + f_517 * hk_645[k]
                  + f_478 * hk_686[k]
                  + f_478 * hk_691[k]
                  - f_518 * hk_693[k]
                  - f_478 * hk_700[k]
                  + f_519 * hk_704[k]
                  - f_478 * hk_713[k]
                  + f_518 * hk_715[k]
                  - f_519 * hk_717[k];
    }

#pragma omp simd aligned(hk_36, hk_39, hk_41, hk_46, hk_48, hk_50, hk_57, hk_59, hk_61, \
                         hk_216, hk_219, hk_221, hk_226, hk_228, hk_230, hk_237, hk_239, \
                         hk_241, hk_288, hk_291, hk_293, hk_298, hk_300, hk_302, hk_309, \
                         hk_311, hk_313, hk_540, hk_543, hk_545, hk_550, hk_552, hk_554, \
                         hk_561, hk_563, hk_565, hk_612, hk_615, hk_617, hk_622, hk_624, \
                         hk_626, hk_633, hk_635, hk_637, hk_684, hk_687, hk_689, hk_694, \
                         hk_696, hk_698, hk_705, hk_707, hk_709 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_70[k] = f_452 * hk_36[k]
                  - f_452 * hk_39[k]
                  - f_455 * hk_41[k]
                  - f_450 * hk_46[k]
                  + f_453 * hk_48[k]
                  + f_456 * hk_50[k]
                  - f_449 * hk_57[k]
                  + f_451 * hk_59[k]
                  - f_454 * hk_61[k]
                  + f_460 * hk_216[k]
                  - f_460 * hk_219[k]
                  - f_453 * hk_221[k]
                  - f_458 * hk_226[k]
                  + f_454 * hk_228[k]
                  + f_462 * hk_230[k]
                  - f_457 * hk_237[k]
                  + f_459 * hk_239[k]
                  - f_461 * hk_241[k]
                  - f_465 * hk_288[k]
                  + f_465 * hk_291[k]
                  + f_468 * hk_293[k]
                  + f_451 * hk_298[k]
                  - f_466 * hk_300[k]
                  - f_469 * hk_302[k]
                  + f_463 * hk_309[k]
                  - f_464 * hk_311[k]
                  + f_467 * hk_313[k]
                  + f_452 * hk_540[k]
                  - f_452 * hk_543[k]
                  - f_455 * hk_545[k]
                  - f_450 * hk_550[k]
                  + f_453 * hk_552[k]
                  + f_456 * hk_554[k]
                  - f_449 * hk_561[k]
                  + f_451 * hk_563[k]
                  - f_454 * hk_565[k]
                  - f_465 * hk_612[k]
                  + f_465 * hk_615[k]
                  + f_468 * hk_617[k]
                  + f_451 * hk_622[k]
                  - f_466 * hk_624[k]
                  - f_469 * hk_626[k]
                  + f_463 * hk_633[k]
                  - f_464 * hk_635[k]
                  + f_467 * hk_637[k]
                  + f_471 * hk_684[k]
                  - f_471 * hk_687[k]
                  - f_461 * hk_689[k]
                  - f_453 * hk_694[k]
                  + f_469 * hk_696[k]
                  + f_473 * hk_698[k]
                  - f_470 * hk_705[k]
                  + f_466 * hk_707[k]
                  - f_472 * hk_709[k];
    }

#pragma omp simd aligned(hk_38, hk_43, hk_45, hk_52, hk_54, hk_65, hk_67, hk_218, hk_223, \
                         hk_225, hk_232, hk_234, hk_245, hk_247, hk_290, hk_295, hk_297, \
                         hk_304, hk_306, hk_317, hk_319, hk_542, hk_547, hk_549, hk_556, \
                         hk_558, hk_569, hk_571, hk_614, hk_619, hk_621, hk_628, hk_630, \
                         hk_641, hk_643, hk_686, hk_691, hk_693, hk_700, hk_702, hk_713, \
                         hk_715 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_71[k] = -f_520 * hk_38[k]
                  + f_521 * hk_43[k]
                  + f_522 * hk_45[k]
                  + f_521 * hk_52[k]
                  - f_425 * hk_54[k]
                  - f_520 * hk_65[k]
                  + f_522 * hk_67[k]
                  - f_427 * hk_218[k]
                  + f_423 * hk_223[k]
                  + f_437 * hk_225[k]
                  + f_423 * hk_232[k]
                  - f_430 * hk_234[k]
                  - f_427 * hk_245[k]
                  + f_437 * hk_247[k]
                  + f_439 * hk_290[k]
                  - f_523 * hk_295[k]
                  - f_430 * hk_297[k]
                  - f_523 * hk_304[k]
                  + f_435 * hk_306[k]
                  + f_439 * hk_317[k]
                  - f_430 * hk_319[k]
                  - f_520 * hk_542[k]
                  + f_521 * hk_547[k]
                  + f_522 * hk_549[k]
                  + f_521 * hk_556[k]
                  - f_425 * hk_558[k]
                  - f_520 * hk_569[k]
                  + f_522 * hk_571[k]
                  + f_439 * hk_614[k]
                  - f_523 * hk_619[k]
                  - f_430 * hk_621[k]
                  - f_523 * hk_628[k]
                  + f_435 * hk_630[k]
                  + f_439 * hk_641[k]
                  - f_430 * hk_643[k]
                  - f_444 * hk_686[k]
                  + f_430 * hk_691[k]
                  + f_445 * hk_693[k]
                  + f_430 * hk_700[k]
                  - f_440 * hk_702[k]
                  - f_444 * hk_713[k]
                  + f_445 * hk_715[k];
    }

#pragma omp simd aligned(hk_36, hk_39, hk_41, hk_46, hk_48, hk_57, hk_59, hk_216, hk_219, \
                         hk_221, hk_226, hk_228, hk_237, hk_239, hk_288, hk_291, hk_293, \
                         hk_298, hk_300, hk_309, hk_311, hk_540, hk_543, hk_545, hk_550, \
                         hk_552, hk_561, hk_563, hk_612, hk_615, hk_617, hk_622, hk_624, \
                         hk_633, hk_635, hk_684, hk_687, hk_689, hk_694, hk_696, hk_705, \
                         hk_707 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_72[k] = -f_426 * hk_36[k]
                  + f_424 * hk_39[k]
                  + f_427 * hk_41[k]
                  + f_422 * hk_46[k]
                  - f_425 * hk_48[k]
                  - f_422 * hk_57[k]
                  + f_423 * hk_59[k]
                  - f_431 * hk_216[k]
                  + f_429 * hk_219[k]
                  + f_432 * hk_221[k]
                  + f_428 * hk_226[k]
                  - f_430 * hk_228[k]
                  - f_428 * hk_237[k]
                  + f_425 * hk_239[k]
                  + f_427 * hk_288[k]
                  - f_434 * hk_291[k]
                  - f_436 * hk_293[k]
                  - f_423 * hk_298[k]
                  + f_435 * hk_300[k]
                  + f_423 * hk_309[k]
                  - f_433 * hk_311[k]
                  - f_426 * hk_540[k]
                  + f_424 * hk_543[k]
                  + f_427 * hk_545[k]
                  + f_422 * hk_550[k]
                  - f_425 * hk_552[k]
                  - f_422 * hk_561[k]
                  + f_423 * hk_563[k]
                  + f_427 * hk_612[k]
                  - f_434 * hk_615[k]
                  - f_436 * hk_617[k]
                  - f_423 * hk_622[k]
                  + f_435 * hk_624[k]
                  + f_423 * hk_633[k]
                  - f_433 * hk_635[k]
                  - f_441 * hk_684[k]
                  + f_439 * hk_687[k]
                  + f_442 * hk_689[k]
                  + f_437 * hk_694[k]
                  - f_440 * hk_696[k]
                  - f_437 * hk_705[k]
                  + f_438 * hk_707[k];
    }

#pragma omp simd aligned(hk_38, hk_43, hk_52, hk_65, hk_218, hk_223, hk_232, hk_245, hk_290, \
                         hk_295, hk_304, hk_317, hk_542, hk_547, hk_556, hk_569, hk_614, \
                         hk_619, hk_628, hk_641, hk_686, hk_691, hk_700, \
                         hk_713 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_73[k] = f_524 * hk_38[k]
                  - f_525 * hk_43[k]
                  + f_525 * hk_52[k]
                  - f_524 * hk_65[k]
                  + f_526 * hk_218[k]
                  - f_527 * hk_223[k]
                  + f_527 * hk_232[k]
                  - f_526 * hk_245[k]
                  - f_416 * hk_290[k]
                  + f_528 * hk_295[k]
                  - f_528 * hk_304[k]
                  + f_416 * hk_317[k]
                  + f_524 * hk_542[k]
                  - f_525 * hk_547[k]
                  + f_525 * hk_556[k]
                  - f_524 * hk_569[k]
                  - f_416 * hk_614[k]
                  + f_528 * hk_619[k]
                  - f_528 * hk_628[k]
                  + f_416 * hk_641[k]
                  + f_529 * hk_686[k]
                  - f_530 * hk_691[k]
                  + f_530 * hk_700[k]
                  - f_529 * hk_713[k];
    }

#pragma omp simd aligned(hk_36, hk_39, hk_46, hk_57, hk_216, hk_219, hk_226, hk_237, hk_288, \
                         hk_291, hk_298, hk_309, hk_540, hk_543, hk_550, hk_561, hk_612, \
                         hk_615, hk_622, hk_633, hk_684, hk_687, hk_694, \
                         hk_705 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_74[k] = f_401 * hk_36[k]
                  - f_400 * hk_39[k]
                  + f_399 * hk_46[k]
                  - f_398 * hk_57[k]
                  + f_405 * hk_216[k]
                  - f_404 * hk_219[k]
                  + f_403 * hk_226[k]
                  - f_402 * hk_237[k]
                  - f_409 * hk_288[k]
                  + f_408 * hk_291[k]
                  - f_407 * hk_298[k]
                  + f_406 * hk_309[k]
                  + f_401 * hk_540[k]
                  - f_400 * hk_543[k]
                  + f_399 * hk_550[k]
                  - f_398 * hk_561[k]
                  - f_409 * hk_612[k]
                  + f_408 * hk_615[k]
                  - f_407 * hk_622[k]
                  + f_406 * hk_633[k]
                  + f_413 * hk_684[k]
                  - f_412 * hk_687[k]
                  + f_411 * hk_694[k]
                  - f_410 * hk_705[k];
    }

#pragma omp simd aligned(hk_73, hk_78, hk_87, hk_100, hk_253, hk_258, hk_267, hk_280, hk_325, \
                         hk_330, hk_339, hk_352, hk_577, hk_582, hk_591, hk_604, hk_649, \
                         hk_654, hk_663, hk_676, hk_721, hk_726, hk_735, \
                         hk_748 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_75[k] = f_110 * hk_73[k]
                  - f_531 * hk_78[k]
                  + f_115 * hk_87[k]
                  - f_532 * hk_100[k]
                  + f_112 * hk_253[k]
                  - f_533 * hk_258[k]
                  + f_11 * hk_267[k]
                  - f_534 * hk_280[k]
                  - f_535 * hk_325[k]
                  + f_536 * hk_330[k]
                  - f_537 * hk_339[k]
                  + f_538 * hk_352[k]
                  + f_110 * hk_577[k]
                  - f_531 * hk_582[k]
                  + f_115 * hk_591[k]
                  - f_532 * hk_604[k]
                  - f_535 * hk_649[k]
                  + f_536 * hk_654[k]
                  - f_537 * hk_663[k]
                  + f_538 * hk_676[k]
                  + f_539 * hk_721[k]
                  - f_535 * hk_726[k]
                  + f_540 * hk_735[k]
                  - f_541 * hk_748[k];
    }

#pragma omp simd aligned(hk_76, hk_83, hk_94, hk_256, hk_263, hk_274, hk_328, hk_335, hk_346, \
                         hk_580, hk_587, hk_598, hk_652, hk_659, hk_670, hk_724, hk_731, \
                         hk_742 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_76[k] = f_542 * hk_76[k]
                  - f_543 * hk_83[k]
                  + f_542 * hk_94[k]
                  + f_544 * hk_256[k]
                  - f_545 * hk_263[k]
                  + f_544 * hk_274[k]
                  - f_546 * hk_328[k]
                  + f_547 * hk_335[k]
                  - f_546 * hk_346[k]
                  + f_542 * hk_580[k]
                  - f_543 * hk_587[k]
                  + f_542 * hk_598[k]
                  - f_546 * hk_652[k]
                  + f_547 * hk_659[k]
                  - f_546 * hk_670[k]
                  + f_548 * hk_724[k]
                  - f_549 * hk_731[k]
                  + f_548 * hk_742[k];
    }

#pragma omp simd aligned(hk_73, hk_78, hk_80, hk_87, hk_89, hk_100, hk_102, hk_253, hk_258, \
                         hk_260, hk_267, hk_269, hk_280, hk_282, hk_325, hk_330, hk_332, \
                         hk_339, hk_341, hk_352, hk_354, hk_577, hk_582, hk_584, hk_591, \
                         hk_593, hk_604, hk_606, hk_649, hk_654, hk_656, hk_663, hk_665, \
                         hk_676, hk_678, hk_721, hk_726, hk_728, hk_735, hk_737, hk_748, \
                         hk_750 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_77[k] = -f_550 * hk_73[k]
                  + f_550 * hk_78[k]
                  + f_551 * hk_80[k]
                  + f_552 * hk_87[k]
                  - f_553 * hk_89[k]
                  - f_554 * hk_100[k]
                  + f_555 * hk_102[k]
                  - f_556 * hk_253[k]
                  + f_556 * hk_258[k]
                  + f_553 * hk_260[k]
                  + f_557 * hk_267[k]
                  - f_558 * hk_269[k]
                  - f_559 * hk_280[k]
                  + f_560 * hk_282[k]
                  + f_561 * hk_325[k]
                  - f_561 * hk_330[k]
                  - f_562 * hk_332[k]
                  - f_560 * hk_339[k]
                  + f_563 * hk_341[k]
                  + f_564 * hk_352[k]
                  - f_565 * hk_354[k]
                  - f_550 * hk_577[k]
                  + f_550 * hk_582[k]
                  + f_551 * hk_584[k]
                  + f_552 * hk_591[k]
                  - f_553 * hk_593[k]
                  - f_554 * hk_604[k]
                  + f_555 * hk_606[k]
                  + f_561 * hk_649[k]
                  - f_561 * hk_654[k]
                  - f_562 * hk_656[k]
                  - f_560 * hk_663[k]
                  + f_563 * hk_665[k]
                  + f_564 * hk_676[k]
                  - f_565 * hk_678[k]
                  - f_564 * hk_721[k]
                  + f_564 * hk_726[k]
                  + f_565 * hk_728[k]
                  + f_566 * hk_735[k]
                  - f_567 * hk_737[k]
                  - f_568 * hk_748[k]
                  + f_569 * hk_750[k];
    }

#pragma omp simd aligned(hk_76, hk_85, hk_94, hk_96, hk_256, hk_265, hk_274, hk_276, hk_328, \
                         hk_337, hk_346, hk_348, hk_580, hk_589, hk_598, hk_600, hk_652, \
                         hk_661, hk_670, hk_672, hk_724, hk_733, hk_742, \
                         hk_744 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_78[k] = -f_560 * hk_76[k]
                  + f_570 * hk_85[k]
                  + f_560 * hk_94[k]
                  - f_570 * hk_96[k]
                  - f_571 * hk_256[k]
                  + f_562 * hk_265[k]
                  + f_571 * hk_274[k]
                  - f_562 * hk_276[k]
                  + f_567 * hk_328[k]
                  - f_572 * hk_337[k]
                  - f_567 * hk_346[k]
                  + f_572 * hk_348[k]
                  - f_560 * hk_580[k]
                  + f_570 * hk_589[k]
                  + f_560 * hk_598[k]
                  - f_570 * hk_600[k]
                  + f_567 * hk_652[k]
                  - f_572 * hk_661[k]
                  - f_567 * hk_670[k]
                  + f_572 * hk_672[k]
                  - f_573 * hk_724[k]
                  + f_574 * hk_733[k]
                  + f_573 * hk_742[k]
                  - f_574 * hk_744[k];
    }

#pragma omp simd aligned(hk_73, hk_78, hk_80, hk_87, hk_89, hk_91, hk_100, hk_102, hk_104, \
                         hk_253, hk_258, hk_260, hk_267, hk_269, hk_271, hk_280, hk_282, \
                         hk_284, hk_325, hk_330, hk_332, hk_339, hk_341, hk_343, hk_352, \
                         hk_354, hk_356, hk_577, hk_582, hk_584, hk_591, hk_593, hk_595, \
                         hk_604, hk_606, hk_608, hk_649, hk_654, hk_656, hk_663, hk_665, \
                         hk_667, hk_676, hk_678, hk_680, hk_721, hk_726, hk_728, hk_735, \
                         hk_737, hk_739, hk_748, hk_750, hk_752 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_79[k] = f_575 * hk_73[k]
                  + f_576 * hk_78[k]
                  - f_577 * hk_80[k]
                  + f_578 * hk_87[k]
                  - f_579 * hk_89[k]
                  + f_580 * hk_91[k]
                  - f_578 * hk_100[k]
                  + f_581 * hk_102[k]
                  - f_582 * hk_104[k]
                  + f_583 * hk_253[k]
                  + f_584 * hk_258[k]
                  - f_585 * hk_260[k]
                  + f_586 * hk_267[k]
                  - f_580 * hk_269[k]
                  + f_587 * hk_271[k]
                  - f_586 * hk_280[k]
                  + f_579 * hk_282[k]
                  - f_588 * hk_284[k]
                  - f_589 * hk_325[k]
                  - f_590 * hk_330[k]
                  + f_587 * hk_332[k]
                  - f_591 * hk_339[k]
                  + f_592 * hk_341[k]
                  - f_593 * hk_343[k]
                  + f_591 * hk_352[k]
                  - f_588 * hk_354[k]
                  + f_594 * hk_356[k]
                  + f_575 * hk_577[k]
                  + f_576 * hk_582[k]
                  - f_577 * hk_584[k]
                  + f_578 * hk_591[k]
                  - f_579 * hk_593[k]
                  + f_580 * hk_595[k]
                  - f_578 * hk_604[k]
                  + f_581 * hk_606[k]
                  - f_582 * hk_608[k]
                  - f_589 * hk_649[k]
                  - f_590 * hk_654[k]
                  + f_587 * hk_656[k]
                  - f_591 * hk_663[k]
                  + f_592 * hk_665[k]
                  - f_593 * hk_667[k]
                  + f_591 * hk_676[k]
                  - f_588 * hk_678[k]
                  + f_594 * hk_680[k]
                  + f_595 * hk_721[k]
                  + f_591 * hk_726[k]
                  - f_596 * hk_728[k]
                  + f_597 * hk_735[k]
                  - f_598 * hk_737[k]
                  + f_599 * hk_739[k]
                  - f_597 * hk_748[k]
                  + f_600 * hk_750[k]
                  - f_601 * hk_752[k];
    }

#pragma omp simd aligned(hk_76, hk_83, hk_85, hk_94, hk_96, hk_98, hk_256, hk_263, hk_265, \
                         hk_274, hk_276, hk_278, hk_328, hk_335, hk_337, hk_346, hk_348, \
                         hk_350, hk_580, hk_587, hk_589, hk_598, hk_600, hk_602, hk_652, \
                         hk_659, hk_661, hk_670, hk_672, hk_674, hk_724, hk_731, hk_733, \
                         hk_742, hk_744, hk_746 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_80[k] = f_602 * hk_76[k]
                  + f_603 * hk_83[k]
                  - f_604 * hk_85[k]
                  + f_602 * hk_94[k]
                  - f_604 * hk_96[k]
                  + f_605 * hk_98[k]
                  + f_603 * hk_256[k]
                  + f_606 * hk_263[k]
                  - f_607 * hk_265[k]
                  + f_603 * hk_274[k]
                  - f_607 * hk_276[k]
                  + f_608 * hk_278[k]
                  - f_609 * hk_328[k]
                  - f_604 * hk_335[k]
                  + f_610 * hk_337[k]
                  - f_609 * hk_346[k]
                  + f_610 * hk_348[k]
                  - f_611 * hk_350[k]
                  + f_602 * hk_580[k]
                  + f_603 * hk_587[k]
                  - f_604 * hk_589[k]
                  + f_602 * hk_598[k]
                  - f_604 * hk_600[k]
                  + f_605 * hk_602[k]
                  - f_609 * hk_652[k]
                  - f_604 * hk_659[k]
                  + f_610 * hk_661[k]
                  - f_609 * hk_670[k]
                  + f_610 * hk_672[k]
                  - f_611 * hk_674[k]
                  + f_612 * hk_724[k]
                  + f_613 * hk_731[k]
                  - f_614 * hk_733[k]
                  + f_612 * hk_742[k]
                  - f_614 * hk_744[k]
                  + f_615 * hk_746[k];
    }

#pragma omp simd aligned(hk_73, hk_78, hk_80, hk_87, hk_89, hk_91, hk_100, hk_102, hk_104, \
                         hk_106, hk_253, hk_258, hk_260, hk_267, hk_269, hk_271, hk_280, \
                         hk_282, hk_284, hk_286, hk_325, hk_330, hk_332, hk_339, hk_341, \
                         hk_343, hk_352, hk_354, hk_356, hk_358, hk_577, hk_582, hk_584, \
                         hk_591, hk_593, hk_595, hk_604, hk_606, hk_608, hk_610, hk_649, \
                         hk_654, hk_656, hk_663, hk_665, hk_667, hk_676, hk_678, hk_680, \
                         hk_682, hk_721, hk_726, hk_728, hk_735, hk_737, hk_739, hk_748, \
                         hk_750, hk_752, hk_754 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_81[k] = -f_616 * hk_73[k]
                  - f_617 * hk_78[k]
                  + f_618 * hk_80[k]
                  - f_617 * hk_87[k]
                  + f_619 * hk_89[k]
                  - f_619 * hk_91[k]
                  - f_616 * hk_100[k]
                  + f_618 * hk_102[k]
                  - f_619 * hk_104[k]
                  + f_620 * hk_106[k]
                  - f_621 * hk_253[k]
                  - f_622 * hk_258[k]
                  + f_619 * hk_260[k]
                  - f_622 * hk_267[k]
                  + f_623 * hk_269[k]
                  - f_623 * hk_271[k]
                  - f_621 * hk_280[k]
                  + f_619 * hk_282[k]
                  - f_623 * hk_284[k]
                  + f_624 * hk_286[k]
                  + f_625 * hk_325[k]
                  + f_626 * hk_330[k]
                  - f_627 * hk_332[k]
                  + f_626 * hk_339[k]
                  - f_628 * hk_341[k]
                  + f_628 * hk_343[k]
                  + f_625 * hk_352[k]
                  - f_627 * hk_354[k]
                  + f_628 * hk_356[k]
                  - f_629 * hk_358[k]
                  - f_616 * hk_577[k]
                  - f_617 * hk_582[k]
                  + f_618 * hk_584[k]
                  - f_617 * hk_591[k]
                  + f_619 * hk_593[k]
                  - f_619 * hk_595[k]
                  - f_616 * hk_604[k]
                  + f_618 * hk_606[k]
                  - f_619 * hk_608[k]
                  + f_620 * hk_610[k]
                  + f_625 * hk_649[k]
                  + f_626 * hk_654[k]
                  - f_627 * hk_656[k]
                  + f_626 * hk_663[k]
                  - f_628 * hk_665[k]
                  + f_628 * hk_667[k]
                  + f_625 * hk_676[k]
                  - f_627 * hk_678[k]
                  + f_628 * hk_680[k]
                  - f_629 * hk_682[k]
                  - f_630 * hk_721[k]
                  - f_631 * hk_726[k]
                  + f_620 * hk_728[k]
                  - f_631 * hk_735[k]
                  + f_624 * hk_737[k]
                  - f_624 * hk_739[k]
                  - f_630 * hk_748[k]
                  + f_620 * hk_750[k]
                  - f_624 * hk_752[k]
                  + f_632 * hk_754[k];
    }

#pragma omp simd aligned(hk_74, hk_79, hk_81, hk_88, hk_90, hk_92, hk_101, hk_103, hk_105, \
                         hk_107, hk_254, hk_259, hk_261, hk_268, hk_270, hk_272, hk_281, \
                         hk_283, hk_285, hk_287, hk_326, hk_331, hk_333, hk_340, hk_342, \
                         hk_344, hk_353, hk_355, hk_357, hk_359, hk_578, hk_583, hk_585, \
                         hk_592, hk_594, hk_596, hk_605, hk_607, hk_609, hk_611, hk_650, \
                         hk_655, hk_657, hk_664, hk_666, hk_668, hk_677, hk_679, hk_681, \
                         hk_683, hk_722, hk_727, hk_729, hk_736, hk_738, hk_740, hk_749, \
                         hk_751, hk_753, hk_755 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_82[k] = -4.1015625 * hk_74[k]
                  - 12.3046875 * hk_79[k]
                  + 24.609375 * hk_81[k]
                  - 12.3046875 * hk_88[k]
                  + 49.21875 * hk_90[k]
                  - 19.6875 * hk_92[k]
                  - 4.1015625 * hk_101[k]
                  + 24.609375 * hk_103[k]
                  - 19.6875 * hk_105[k]
                  + 1.875 * hk_107[k]
                  - 8.203125 * hk_254[k]
                  - 24.609375 * hk_259[k]
                  + 49.21875 * hk_261[k]
                  - 24.609375 * hk_268[k]
                  + 98.4375 * hk_270[k]
                  - 39.375 * hk_272[k]
                  - 8.203125 * hk_281[k]
                  + 49.21875 * hk_283[k]
                  - 39.375 * hk_285[k]
                  + 3.75 * hk_287[k]
                  + 10.9375 * hk_326[k]
                  + 32.8125 * hk_331[k]
                  - 65.625 * hk_333[k]
                  + 32.8125 * hk_340[k]
                  - 131.25 * hk_342[k]
                  + 52.5 * hk_344[k]
                  + 10.9375 * hk_353[k]
                  - 65.625 * hk_355[k]
                  + 52.5 * hk_357[k]
                  - 5.0 * hk_359[k]
                  - 4.1015625 * hk_578[k]
                  - 12.3046875 * hk_583[k]
                  + 24.609375 * hk_585[k]
                  - 12.3046875 * hk_592[k]
                  + 49.21875 * hk_594[k]
                  - 19.6875 * hk_596[k]
                  - 4.1015625 * hk_605[k]
                  + 24.609375 * hk_607[k]
                  - 19.6875 * hk_609[k]
                  + 1.875 * hk_611[k]
                  + 10.9375 * hk_650[k]
                  + 32.8125 * hk_655[k]
                  - 65.625 * hk_657[k]
                  + 32.8125 * hk_664[k]
                  - 131.25 * hk_666[k]
                  + 52.5 * hk_668[k]
                  + 10.9375 * hk_677[k]
                  - 65.625 * hk_679[k]
                  + 52.5 * hk_681[k]
                  - 5.0 * hk_683[k]
                  - 2.1875 * hk_722[k]
                  - 6.5625 * hk_727[k]
                  + 13.125 * hk_729[k]
                  - 6.5625 * hk_736[k]
                  + 26.25 * hk_738[k]
                  - 10.5 * hk_740[k]
                  - 2.1875 * hk_749[k]
                  + 13.125 * hk_751[k]
                  - 10.5 * hk_753[k]
                  + hk_755[k];
    }

#pragma omp simd aligned(hk_72, hk_75, hk_77, hk_82, hk_84, hk_86, hk_93, hk_95, hk_97, hk_99, \
                         hk_252, hk_255, hk_257, hk_262, hk_264, hk_266, hk_273, hk_275, \
                         hk_277, hk_279, hk_324, hk_327, hk_329, hk_334, hk_336, hk_338, \
                         hk_345, hk_347, hk_349, hk_351, hk_576, hk_579, hk_581, hk_586, \
                         hk_588, hk_590, hk_597, hk_599, hk_601, hk_603, hk_648, hk_651, \
                         hk_653, hk_658, hk_660, hk_662, hk_669, hk_671, hk_673, hk_675, \
                         hk_720, hk_723, hk_725, hk_730, hk_732, hk_734, hk_741, hk_743, \
                         hk_745, hk_747 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_83[k] = -f_616 * hk_72[k]
                  - f_617 * hk_75[k]
                  + f_618 * hk_77[k]
                  - f_617 * hk_82[k]
                  + f_619 * hk_84[k]
                  - f_619 * hk_86[k]
                  - f_616 * hk_93[k]
                  + f_618 * hk_95[k]
                  - f_619 * hk_97[k]
                  + f_620 * hk_99[k]
                  - f_621 * hk_252[k]
                  - f_622 * hk_255[k]
                  + f_619 * hk_257[k]
                  - f_622 * hk_262[k]
                  + f_623 * hk_264[k]
                  - f_623 * hk_266[k]
                  - f_621 * hk_273[k]
                  + f_619 * hk_275[k]
                  - f_623 * hk_277[k]
                  + f_624 * hk_279[k]
                  + f_625 * hk_324[k]
                  + f_626 * hk_327[k]
                  - f_627 * hk_329[k]
                  + f_626 * hk_334[k]
                  - f_628 * hk_336[k]
                  + f_628 * hk_338[k]
                  + f_625 * hk_345[k]
                  - f_627 * hk_347[k]
                  + f_628 * hk_349[k]
                  - f_629 * hk_351[k]
                  - f_616 * hk_576[k]
                  - f_617 * hk_579[k]
                  + f_618 * hk_581[k]
                  - f_617 * hk_586[k]
                  + f_619 * hk_588[k]
                  - f_619 * hk_590[k]
                  - f_616 * hk_597[k]
                  + f_618 * hk_599[k]
                  - f_619 * hk_601[k]
                  + f_620 * hk_603[k]
                  + f_625 * hk_648[k]
                  + f_626 * hk_651[k]
                  - f_627 * hk_653[k]
                  + f_626 * hk_658[k]
                  - f_628 * hk_660[k]
                  + f_628 * hk_662[k]
                  + f_625 * hk_669[k]
                  - f_627 * hk_671[k]
                  + f_628 * hk_673[k]
                  - f_629 * hk_675[k]
                  - f_630 * hk_720[k]
                  - f_631 * hk_723[k]
                  + f_620 * hk_725[k]
                  - f_631 * hk_730[k]
                  + f_624 * hk_732[k]
                  - f_624 * hk_734[k]
                  - f_630 * hk_741[k]
                  + f_620 * hk_743[k]
                  - f_624 * hk_745[k]
                  + f_632 * hk_747[k];
    }

#pragma omp simd aligned(hk_74, hk_79, hk_81, hk_88, hk_92, hk_101, hk_103, hk_105, hk_254, \
                         hk_259, hk_261, hk_268, hk_272, hk_281, hk_283, hk_285, hk_326, \
                         hk_331, hk_333, hk_340, hk_344, hk_353, hk_355, hk_357, hk_578, \
                         hk_583, hk_585, hk_592, hk_596, hk_605, hk_607, hk_609, hk_650, \
                         hk_655, hk_657, hk_664, hk_668, hk_677, hk_679, hk_681, hk_722, \
                         hk_727, hk_729, hk_736, hk_740, hk_749, hk_751, \
                         hk_753 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_84[k] = f_633 * hk_74[k]
                  + f_633 * hk_79[k]
                  - f_609 * hk_81[k]
                  - f_633 * hk_88[k]
                  + f_634 * hk_92[k]
                  - f_633 * hk_101[k]
                  + f_609 * hk_103[k]
                  - f_634 * hk_105[k]
                  + f_602 * hk_254[k]
                  + f_602 * hk_259[k]
                  - f_604 * hk_261[k]
                  - f_602 * hk_268[k]
                  + f_605 * hk_272[k]
                  - f_602 * hk_281[k]
                  + f_604 * hk_283[k]
                  - f_605 * hk_285[k]
                  - f_635 * hk_326[k]
                  - f_635 * hk_331[k]
                  + f_636 * hk_333[k]
                  + f_635 * hk_340[k]
                  - f_637 * hk_344[k]
                  + f_635 * hk_353[k]
                  - f_636 * hk_355[k]
                  + f_637 * hk_357[k]
                  + f_633 * hk_578[k]
                  + f_633 * hk_583[k]
                  - f_609 * hk_585[k]
                  - f_633 * hk_592[k]
                  + f_634 * hk_596[k]
                  - f_633 * hk_605[k]
                  + f_609 * hk_607[k]
                  - f_634 * hk_609[k]
                  - f_635 * hk_650[k]
                  - f_635 * hk_655[k]
                  + f_636 * hk_657[k]
                  + f_635 * hk_664[k]
                  - f_637 * hk_668[k]
                  + f_635 * hk_677[k]
                  - f_636 * hk_679[k]
                  + f_637 * hk_681[k]
                  + f_638 * hk_722[k]
                  + f_638 * hk_727[k]
                  - f_639 * hk_729[k]
                  - f_638 * hk_736[k]
                  + f_640 * hk_740[k]
                  - f_638 * hk_749[k]
                  + f_639 * hk_751[k]
                  - f_640 * hk_753[k];
    }

#pragma omp simd aligned(hk_72, hk_75, hk_77, hk_82, hk_84, hk_86, hk_93, hk_95, hk_97, \
                         hk_252, hk_255, hk_257, hk_262, hk_264, hk_266, hk_273, hk_275, \
                         hk_277, hk_324, hk_327, hk_329, hk_334, hk_336, hk_338, hk_345, \
                         hk_347, hk_349, hk_576, hk_579, hk_581, hk_586, hk_588, hk_590, \
                         hk_597, hk_599, hk_601, hk_648, hk_651, hk_653, hk_658, hk_660, \
                         hk_662, hk_669, hk_671, hk_673, hk_720, hk_723, hk_725, hk_730, \
                         hk_732, hk_734, hk_741, hk_743, hk_745 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_85[k] = f_578 * hk_72[k]
                  - f_578 * hk_75[k]
                  - f_581 * hk_77[k]
                  - f_576 * hk_82[k]
                  + f_579 * hk_84[k]
                  + f_582 * hk_86[k]
                  - f_575 * hk_93[k]
                  + f_577 * hk_95[k]
                  - f_580 * hk_97[k]
                  + f_586 * hk_252[k]
                  - f_586 * hk_255[k]
                  - f_579 * hk_257[k]
                  - f_584 * hk_262[k]
                  + f_580 * hk_264[k]
                  + f_588 * hk_266[k]
                  - f_583 * hk_273[k]
                  + f_585 * hk_275[k]
                  - f_587 * hk_277[k]
                  - f_591 * hk_324[k]
                  + f_591 * hk_327[k]
                  + f_588 * hk_329[k]
                  + f_590 * hk_334[k]
                  - f_592 * hk_336[k]
                  - f_594 * hk_338[k]
                  + f_589 * hk_345[k]
                  - f_587 * hk_347[k]
                  + f_593 * hk_349[k]
                  + f_578 * hk_576[k]
                  - f_578 * hk_579[k]
                  - f_581 * hk_581[k]
                  - f_576 * hk_586[k]
                  + f_579 * hk_588[k]
                  + f_582 * hk_590[k]
                  - f_575 * hk_597[k]
                  + f_577 * hk_599[k]
                  - f_580 * hk_601[k]
                  - f_591 * hk_648[k]
                  + f_591 * hk_651[k]
                  + f_588 * hk_653[k]
                  + f_590 * hk_658[k]
                  - f_592 * hk_660[k]
                  - f_594 * hk_662[k]
                  + f_589 * hk_669[k]
                  - f_587 * hk_671[k]
                  + f_593 * hk_673[k]
                  + f_597 * hk_720[k]
                  - f_597 * hk_723[k]
                  - f_600 * hk_725[k]
                  - f_591 * hk_730[k]
                  + f_598 * hk_732[k]
                  + f_601 * hk_734[k]
                  - f_595 * hk_741[k]
                  + f_596 * hk_743[k]
                  - f_599 * hk_745[k];
    }

#pragma omp simd aligned(hk_74, hk_79, hk_81, hk_88, hk_90, hk_101, hk_103, hk_254, hk_259, \
                         hk_261, hk_268, hk_270, hk_281, hk_283, hk_326, hk_331, hk_333, \
                         hk_340, hk_342, hk_353, hk_355, hk_578, hk_583, hk_585, hk_592, \
                         hk_594, hk_605, hk_607, hk_650, hk_655, hk_657, hk_664, hk_666, \
                         hk_677, hk_679, hk_722, hk_727, hk_729, hk_736, hk_738, hk_749, \
                         hk_751 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_86[k] = -f_641 * hk_74[k]
                  + f_642 * hk_79[k]
                  + f_643 * hk_81[k]
                  + f_642 * hk_88[k]
                  - f_553 * hk_90[k]
                  - f_641 * hk_101[k]
                  + f_643 * hk_103[k]
                  - f_555 * hk_254[k]
                  + f_551 * hk_259[k]
                  + f_644 * hk_261[k]
                  + f_551 * hk_268[k]
                  - f_558 * hk_270[k]
                  - f_555 * hk_281[k]
                  + f_644 * hk_283[k]
                  + f_645 * hk_326[k]
                  - f_570 * hk_331[k]
                  - f_646 * hk_333[k]
                  - f_570 * hk_340[k]
                  + f_563 * hk_342[k]
                  + f_645 * hk_353[k]
                  - f_646 * hk_355[k]
                  - f_641 * hk_578[k]
                  + f_642 * hk_583[k]
                  + f_643 * hk_585[k]
                  + f_642 * hk_592[k]
                  - f_553 * hk_594[k]
                  - f_641 * hk_605[k]
                  + f_643 * hk_607[k]
                  + f_645 * hk_650[k]
                  - f_570 * hk_655[k]
                  - f_646 * hk_657[k]
                  - f_570 * hk_664[k]
                  + f_563 * hk_666[k]
                  + f_645 * hk_677[k]
                  - f_646 * hk_679[k]
                  - f_647 * hk_722[k]
                  + f_645 * hk_727[k]
                  + f_648 * hk_729[k]
                  + f_645 * hk_736[k]
                  - f_567 * hk_738[k]
                  - f_647 * hk_749[k]
                  + f_648 * hk_751[k];
    }

#pragma omp simd aligned(hk_72, hk_75, hk_77, hk_82, hk_84, hk_93, hk_95, hk_252, hk_255, \
                         hk_257, hk_262, hk_264, hk_273, hk_275, hk_324, hk_327, hk_329, \
                         hk_334, hk_336, hk_345, hk_347, hk_576, hk_579, hk_581, hk_586, \
                         hk_588, hk_597, hk_599, hk_648, hk_651, hk_653, hk_658, hk_660, \
                         hk_669, hk_671, hk_720, hk_723, hk_725, hk_730, hk_732, hk_741, \
                         hk_743 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_87[k] = -f_554 * hk_72[k]
                  + f_552 * hk_75[k]
                  + f_555 * hk_77[k]
                  + f_550 * hk_82[k]
                  - f_553 * hk_84[k]
                  - f_550 * hk_93[k]
                  + f_551 * hk_95[k]
                  - f_559 * hk_252[k]
                  + f_557 * hk_255[k]
                  + f_560 * hk_257[k]
                  + f_556 * hk_262[k]
                  - f_558 * hk_264[k]
                  - f_556 * hk_273[k]
                  + f_553 * hk_275[k]
                  + f_564 * hk_324[k]
                  - f_560 * hk_327[k]
                  - f_565 * hk_329[k]
                  - f_561 * hk_334[k]
                  + f_563 * hk_336[k]
                  + f_561 * hk_345[k]
                  - f_562 * hk_347[k]
                  - f_554 * hk_576[k]
                  + f_552 * hk_579[k]
                  + f_555 * hk_581[k]
                  + f_550 * hk_586[k]
                  - f_553 * hk_588[k]
                  - f_550 * hk_597[k]
                  + f_551 * hk_599[k]
                  + f_564 * hk_648[k]
                  - f_560 * hk_651[k]
                  - f_565 * hk_653[k]
                  - f_561 * hk_658[k]
                  + f_563 * hk_660[k]
                  + f_561 * hk_669[k]
                  - f_562 * hk_671[k]
                  - f_568 * hk_720[k]
                  + f_566 * hk_723[k]
                  + f_569 * hk_725[k]
                  + f_564 * hk_730[k]
                  - f_567 * hk_732[k]
                  - f_564 * hk_741[k]
                  + f_565 * hk_743[k];
    }

#pragma omp simd aligned(hk_74, hk_79, hk_88, hk_101, hk_254, hk_259, hk_268, hk_281, hk_326, \
                         hk_331, hk_340, hk_353, hk_578, hk_583, hk_592, hk_605, hk_650, \
                         hk_655, hk_664, hk_677, hk_722, hk_727, hk_736, \
                         hk_749 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_88[k] = f_7 * hk_74[k]
                  - f_649 * hk_79[k]
                  + f_649 * hk_88[k]
                  - f_7 * hk_101[k]
                  + f_650 * hk_254[k]
                  - f_651 * hk_259[k]
                  + f_651 * hk_268[k]
                  - f_650 * hk_281[k]
                  - f_652 * hk_326[k]
                  + f_545 * hk_331[k]
                  - f_545 * hk_340[k]
                  + f_652 * hk_353[k]
                  + f_7 * hk_578[k]
                  - f_649 * hk_583[k]
                  + f_649 * hk_592[k]
                  - f_7 * hk_605[k]
                  - f_652 * hk_650[k]
                  + f_545 * hk_655[k]
                  - f_545 * hk_664[k]
                  + f_652 * hk_677[k]
                  + f_653 * hk_722[k]
                  - f_654 * hk_727[k]
                  + f_654 * hk_736[k]
                  - f_653 * hk_749[k];
    }

#pragma omp simd aligned(hk_72, hk_75, hk_82, hk_93, hk_252, hk_255, hk_262, hk_273, hk_324, \
                         hk_327, hk_334, hk_345, hk_576, hk_579, hk_586, hk_597, hk_648, \
                         hk_651, hk_658, hk_669, hk_720, hk_723, hk_730, \
                         hk_741 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_89[k] = f_532 * hk_72[k]
                  - f_115 * hk_75[k]
                  + f_531 * hk_82[k]
                  - f_110 * hk_93[k]
                  + f_534 * hk_252[k]
                  - f_11 * hk_255[k]
                  + f_533 * hk_262[k]
                  - f_112 * hk_273[k]
                  - f_538 * hk_324[k]
                  + f_537 * hk_327[k]
                  - f_536 * hk_334[k]
                  + f_535 * hk_345[k]
                  + f_532 * hk_576[k]
                  - f_115 * hk_579[k]
                  + f_531 * hk_586[k]
                  - f_110 * hk_597[k]
                  - f_538 * hk_648[k]
                  + f_537 * hk_651[k]
                  - f_536 * hk_658[k]
                  + f_535 * hk_669[k]
                  + f_541 * hk_720[k]
                  - f_540 * hk_723[k]
                  + f_535 * hk_730[k]
                  - f_539 * hk_741[k];
    }

#pragma omp simd aligned(hk_1, hk_6, hk_15, hk_28, hk_109, hk_114, hk_123, hk_136, hk_181, \
                         hk_186, hk_195, hk_208, hk_361, hk_366, hk_375, hk_388, hk_433, \
                         hk_438, hk_447, hk_460, hk_505, hk_510, hk_519, \
                         hk_532 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_90[k] = f_398 * hk_1[k]
                  - f_399 * hk_6[k]
                  + f_400 * hk_15[k]
                  - f_401 * hk_28[k]
                  + f_402 * hk_109[k]
                  - f_403 * hk_114[k]
                  + f_404 * hk_123[k]
                  - f_405 * hk_136[k]
                  - f_406 * hk_181[k]
                  + f_407 * hk_186[k]
                  - f_408 * hk_195[k]
                  + f_409 * hk_208[k]
                  + f_398 * hk_361[k]
                  - f_399 * hk_366[k]
                  + f_400 * hk_375[k]
                  - f_401 * hk_388[k]
                  - f_406 * hk_433[k]
                  + f_407 * hk_438[k]
                  - f_408 * hk_447[k]
                  + f_409 * hk_460[k]
                  + f_410 * hk_505[k]
                  - f_411 * hk_510[k]
                  + f_412 * hk_519[k]
                  - f_413 * hk_532[k];
    }

#pragma omp simd aligned(hk_4, hk_11, hk_22, hk_112, hk_119, hk_130, hk_184, hk_191, hk_202, \
                         hk_364, hk_371, hk_382, hk_436, hk_443, hk_454, hk_508, hk_515, \
                         hk_526 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_91[k] = f_414 * hk_4[k]
                  - f_415 * hk_11[k]
                  + f_414 * hk_22[k]
                  + f_416 * hk_112[k]
                  - f_417 * hk_119[k]
                  + f_416 * hk_130[k]
                  - f_418 * hk_184[k]
                  + f_419 * hk_191[k]
                  - f_418 * hk_202[k]
                  + f_414 * hk_364[k]
                  - f_415 * hk_371[k]
                  + f_414 * hk_382[k]
                  - f_418 * hk_436[k]
                  + f_419 * hk_443[k]
                  - f_418 * hk_454[k]
                  + f_420 * hk_508[k]
                  - f_421 * hk_515[k]
                  + f_420 * hk_526[k];
    }

#pragma omp simd aligned(hk_1, hk_6, hk_8, hk_15, hk_17, hk_28, hk_30, hk_109, hk_114, hk_116, \
                         hk_123, hk_125, hk_136, hk_138, hk_181, hk_186, hk_188, hk_195, \
                         hk_197, hk_208, hk_210, hk_361, hk_366, hk_368, hk_375, hk_377, \
                         hk_388, hk_390, hk_433, hk_438, hk_440, hk_447, hk_449, hk_460, \
                         hk_462, hk_505, hk_510, hk_512, hk_519, hk_521, hk_532, \
                         hk_534 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_92[k] = -f_422 * hk_1[k]
                  + f_422 * hk_6[k]
                  + f_423 * hk_8[k]
                  + f_424 * hk_15[k]
                  - f_425 * hk_17[k]
                  - f_426 * hk_28[k]
                  + f_427 * hk_30[k]
                  - f_428 * hk_109[k]
                  + f_428 * hk_114[k]
                  + f_425 * hk_116[k]
                  + f_429 * hk_123[k]
                  - f_430 * hk_125[k]
                  - f_431 * hk_136[k]
                  + f_432 * hk_138[k]
                  + f_423 * hk_181[k]
                  - f_423 * hk_186[k]
                  - f_433 * hk_188[k]
                  - f_434 * hk_195[k]
                  + f_435 * hk_197[k]
                  + f_427 * hk_208[k]
                  - f_436 * hk_210[k]
                  - f_422 * hk_361[k]
                  + f_422 * hk_366[k]
                  + f_423 * hk_368[k]
                  + f_424 * hk_375[k]
                  - f_425 * hk_377[k]
                  - f_426 * hk_388[k]
                  + f_427 * hk_390[k]
                  + f_423 * hk_433[k]
                  - f_423 * hk_438[k]
                  - f_433 * hk_440[k]
                  - f_434 * hk_447[k]
                  + f_435 * hk_449[k]
                  + f_427 * hk_460[k]
                  - f_436 * hk_462[k]
                  - f_437 * hk_505[k]
                  + f_437 * hk_510[k]
                  + f_438 * hk_512[k]
                  + f_439 * hk_519[k]
                  - f_440 * hk_521[k]
                  - f_441 * hk_532[k]
                  + f_442 * hk_534[k];
    }

#pragma omp simd aligned(hk_4, hk_13, hk_22, hk_24, hk_112, hk_121, hk_130, hk_132, hk_184, \
                         hk_193, hk_202, hk_204, hk_364, hk_373, hk_382, hk_384, hk_436, \
                         hk_445, hk_454, hk_456, hk_508, hk_517, hk_526, \
                         hk_528 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_93[k] = -f_432 * hk_4[k]
                  + f_443 * hk_13[k]
                  + f_432 * hk_22[k]
                  - f_443 * hk_24[k]
                  - f_444 * hk_112[k]
                  + f_445 * hk_121[k]
                  + f_444 * hk_130[k]
                  - f_445 * hk_132[k]
                  + f_446 * hk_184[k]
                  - f_440 * hk_193[k]
                  - f_446 * hk_202[k]
                  + f_440 * hk_204[k]
                  - f_432 * hk_364[k]
                  + f_443 * hk_373[k]
                  + f_432 * hk_382[k]
                  - f_443 * hk_384[k]
                  + f_446 * hk_436[k]
                  - f_440 * hk_445[k]
                  - f_446 * hk_454[k]
                  + f_440 * hk_456[k]
                  - f_447 * hk_508[k]
                  + f_448 * hk_517[k]
                  + f_447 * hk_526[k]
                  - f_448 * hk_528[k];
    }

#pragma omp simd aligned(hk_1, hk_6, hk_8, hk_15, hk_17, hk_19, hk_28, hk_30, hk_32, hk_109, \
                         hk_114, hk_116, hk_123, hk_125, hk_127, hk_136, hk_138, hk_140, \
                         hk_181, hk_186, hk_188, hk_195, hk_197, hk_199, hk_208, hk_210, \
                         hk_212, hk_361, hk_366, hk_368, hk_375, hk_377, hk_379, hk_388, \
                         hk_390, hk_392, hk_433, hk_438, hk_440, hk_447, hk_449, hk_451, \
                         hk_460, hk_462, hk_464, hk_505, hk_510, hk_512, hk_519, hk_521, \
                         hk_523, hk_532, hk_534, hk_536 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_94[k] = f_449 * hk_1[k]
                  + f_450 * hk_6[k]
                  - f_451 * hk_8[k]
                  + f_452 * hk_15[k]
                  - f_453 * hk_17[k]
                  + f_454 * hk_19[k]
                  - f_452 * hk_28[k]
                  + f_455 * hk_30[k]
                  - f_456 * hk_32[k]
                  + f_457 * hk_109[k]
                  + f_458 * hk_114[k]
                  - f_459 * hk_116[k]
                  + f_460 * hk_123[k]
                  - f_454 * hk_125[k]
                  + f_461 * hk_127[k]
                  - f_460 * hk_136[k]
                  + f_453 * hk_138[k]
                  - f_462 * hk_140[k]
                  - f_463 * hk_181[k]
                  - f_451 * hk_186[k]
                  + f_464 * hk_188[k]
                  - f_465 * hk_195[k]
                  + f_466 * hk_197[k]
                  - f_467 * hk_199[k]
                  + f_465 * hk_208[k]
                  - f_468 * hk_210[k]
                  + f_469 * hk_212[k]
                  + f_449 * hk_361[k]
                  + f_450 * hk_366[k]
                  - f_451 * hk_368[k]
                  + f_452 * hk_375[k]
                  - f_453 * hk_377[k]
                  + f_454 * hk_379[k]
                  - f_452 * hk_388[k]
                  + f_455 * hk_390[k]
                  - f_456 * hk_392[k]
                  - f_463 * hk_433[k]
                  - f_451 * hk_438[k]
                  + f_464 * hk_440[k]
                  - f_465 * hk_447[k]
                  + f_466 * hk_449[k]
                  - f_467 * hk_451[k]
                  + f_465 * hk_460[k]
                  - f_468 * hk_462[k]
                  + f_469 * hk_464[k]
                  + f_470 * hk_505[k]
                  + f_453 * hk_510[k]
                  - f_466 * hk_512[k]
                  + f_471 * hk_519[k]
                  - f_469 * hk_521[k]
                  + f_472 * hk_523[k]
                  - f_471 * hk_532[k]
                  + f_461 * hk_534[k]
                  - f_473 * hk_536[k];
    }

#pragma omp simd aligned(hk_4, hk_11, hk_13, hk_22, hk_24, hk_26, hk_112, hk_119, hk_121, \
                         hk_130, hk_132, hk_134, hk_184, hk_191, hk_193, hk_202, hk_204, \
                         hk_206, hk_364, hk_371, hk_373, hk_382, hk_384, hk_386, hk_436, \
                         hk_443, hk_445, hk_454, hk_456, hk_458, hk_508, hk_515, hk_517, \
                         hk_526, hk_528, hk_530 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_95[k] = f_474 * hk_4[k]
                  + f_475 * hk_11[k]
                  - f_476 * hk_13[k]
                  + f_474 * hk_22[k]
                  - f_476 * hk_24[k]
                  + f_477 * hk_26[k]
                  + f_475 * hk_112[k]
                  + f_478 * hk_119[k]
                  - f_479 * hk_121[k]
                  + f_475 * hk_130[k]
                  - f_479 * hk_132[k]
                  + f_480 * hk_134[k]
                  - f_481 * hk_184[k]
                  - f_482 * hk_191[k]
                  + f_483 * hk_193[k]
                  - f_481 * hk_202[k]
                  + f_483 * hk_204[k]
                  - f_484 * hk_206[k]
                  + f_474 * hk_364[k]
                  + f_475 * hk_371[k]
                  - f_476 * hk_373[k]
                  + f_474 * hk_382[k]
                  - f_476 * hk_384[k]
                  + f_477 * hk_386[k]
                  - f_481 * hk_436[k]
                  - f_482 * hk_443[k]
                  + f_483 * hk_445[k]
                  - f_481 * hk_454[k]
                  + f_483 * hk_456[k]
                  - f_484 * hk_458[k]
                  + f_485 * hk_508[k]
                  + f_486 * hk_515[k]
                  - f_487 * hk_517[k]
                  + f_485 * hk_526[k]
                  - f_487 * hk_528[k]
                  + f_488 * hk_530[k];
    }

#pragma omp simd aligned(hk_1, hk_6, hk_8, hk_15, hk_17, hk_19, hk_28, hk_30, hk_32, hk_34, \
                         hk_109, hk_114, hk_116, hk_123, hk_125, hk_127, hk_136, hk_138, \
                         hk_140, hk_142, hk_181, hk_186, hk_188, hk_195, hk_197, hk_199, \
                         hk_208, hk_210, hk_212, hk_214, hk_361, hk_366, hk_368, hk_375, \
                         hk_377, hk_379, hk_388, hk_390, hk_392, hk_394, hk_433, hk_438, \
                         hk_440, hk_447, hk_449, hk_451, hk_460, hk_462, hk_464, hk_466, \
                         hk_505, hk_510, hk_512, hk_519, hk_521, hk_523, hk_532, hk_534, \
                         hk_536, hk_538 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_96[k] = -f_489 * hk_1[k]
                  - f_490 * hk_6[k]
                  + f_491 * hk_8[k]
                  - f_490 * hk_15[k]
                  + f_492 * hk_17[k]
                  - f_492 * hk_19[k]
                  - f_489 * hk_28[k]
                  + f_491 * hk_30[k]
                  - f_492 * hk_32[k]
                  + f_493 * hk_34[k]
                  - f_494 * hk_109[k]
                  - f_495 * hk_114[k]
                  + f_492 * hk_116[k]
                  - f_495 * hk_123[k]
                  + f_496 * hk_125[k]
                  - f_496 * hk_127[k]
                  - f_494 * hk_136[k]
                  + f_492 * hk_138[k]
                  - f_496 * hk_140[k]
                  + f_384 * hk_142[k]
                  + f_497 * hk_181[k]
                  + f_498 * hk_186[k]
                  - f_499 * hk_188[k]
                  + f_498 * hk_195[k]
                  - f_500 * hk_197[k]
                  + f_500 * hk_199[k]
                  + f_497 * hk_208[k]
                  - f_499 * hk_210[k]
                  + f_500 * hk_212[k]
                  - f_501 * hk_214[k]
                  - f_489 * hk_361[k]
                  - f_490 * hk_366[k]
                  + f_491 * hk_368[k]
                  - f_490 * hk_375[k]
                  + f_492 * hk_377[k]
                  - f_492 * hk_379[k]
                  - f_489 * hk_388[k]
                  + f_491 * hk_390[k]
                  - f_492 * hk_392[k]
                  + f_493 * hk_394[k]
                  + f_497 * hk_433[k]
                  + f_498 * hk_438[k]
                  - f_499 * hk_440[k]
                  + f_498 * hk_447[k]
                  - f_500 * hk_449[k]
                  + f_500 * hk_451[k]
                  + f_497 * hk_460[k]
                  - f_499 * hk_462[k]
                  + f_500 * hk_464[k]
                  - f_501 * hk_466[k]
                  - f_502 * hk_505[k]
                  - f_491 * hk_510[k]
                  + f_503 * hk_512[k]
                  - f_491 * hk_519[k]
                  + f_504 * hk_521[k]
                  - f_504 * hk_523[k]
                  - f_502 * hk_532[k]
                  + f_503 * hk_534[k]
                  - f_504 * hk_536[k]
                  + f_505 * hk_538[k];
    }

#pragma omp simd aligned(hk_2, hk_7, hk_9, hk_16, hk_18, hk_20, hk_29, hk_31, hk_33, hk_35, \
                         hk_110, hk_115, hk_117, hk_124, hk_126, hk_128, hk_137, hk_139, \
                         hk_141, hk_143, hk_182, hk_187, hk_189, hk_196, hk_198, hk_200, \
                         hk_209, hk_211, hk_213, hk_215, hk_362, hk_367, hk_369, hk_376, \
                         hk_378, hk_380, hk_389, hk_391, hk_393, hk_395, hk_434, hk_439, \
                         hk_441, hk_448, hk_450, hk_452, hk_461, hk_463, hk_465, hk_467, \
                         hk_506, hk_511, hk_513, hk_520, hk_522, hk_524, hk_533, hk_535, \
                         hk_537, hk_539 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_97[k] = -f_506 * hk_2[k]
                  - f_256 * hk_7[k]
                  + f_250 * hk_9[k]
                  - f_256 * hk_16[k]
                  + f_251 * hk_18[k]
                  - f_314 * hk_20[k]
                  - f_506 * hk_29[k]
                  + f_250 * hk_31[k]
                  - f_314 * hk_33[k]
                  + f_507 * hk_35[k]
                  - f_375 * hk_110[k]
                  - f_250 * hk_115[k]
                  + f_251 * hk_117[k]
                  - f_250 * hk_124[k]
                  + f_259 * hk_126[k]
                  - f_258 * hk_128[k]
                  - f_375 * hk_137[k]
                  + f_251 * hk_139[k]
                  - f_258 * hk_141[k]
                  + f_508 * hk_143[k]
                  + f_251 * hk_182[k]
                  + f_309 * hk_187[k]
                  - f_136 * hk_189[k]
                  + f_309 * hk_196[k]
                  - f_134 * hk_198[k]
                  + f_509 * hk_200[k]
                  + f_251 * hk_209[k]
                  - f_136 * hk_211[k]
                  + f_509 * hk_213[k]
                  - f_510 * hk_215[k]
                  - f_506 * hk_362[k]
                  - f_256 * hk_367[k]
                  + f_250 * hk_369[k]
                  - f_256 * hk_376[k]
                  + f_251 * hk_378[k]
                  - f_314 * hk_380[k]
                  - f_506 * hk_389[k]
                  + f_250 * hk_391[k]
                  - f_314 * hk_393[k]
                  + f_507 * hk_395[k]
                  + f_251 * hk_434[k]
                  + f_309 * hk_439[k]
                  - f_136 * hk_441[k]
                  + f_309 * hk_448[k]
                  - f_134 * hk_450[k]
                  + f_509 * hk_452[k]
                  + f_251 * hk_461[k]
                  - f_136 * hk_463[k]
                  + f_509 * hk_465[k]
                  - f_510 * hk_467[k]
                  - f_313 * hk_506[k]
                  - f_259 * hk_511[k]
                  + f_248 * hk_513[k]
                  - f_259 * hk_520[k]
                  + f_137 * hk_522[k]
                  - f_316 * hk_524[k]
                  - f_313 * hk_533[k]
                  + f_248 * hk_535[k]
                  - f_316 * hk_537[k]
                  + f_511 * hk_539[k];
    }

#pragma omp simd aligned(hk_0, hk_3, hk_5, hk_10, hk_12, hk_14, hk_21, hk_23, hk_25, hk_27, \
                         hk_108, hk_111, hk_113, hk_118, hk_120, hk_122, hk_129, hk_131, \
                         hk_133, hk_135, hk_180, hk_183, hk_185, hk_190, hk_192, hk_194, \
                         hk_201, hk_203, hk_205, hk_207, hk_360, hk_363, hk_365, hk_370, \
                         hk_372, hk_374, hk_381, hk_383, hk_385, hk_387, hk_432, hk_435, \
                         hk_437, hk_442, hk_444, hk_446, hk_453, hk_455, hk_457, hk_459, \
                         hk_504, hk_507, hk_509, hk_514, hk_516, hk_518, hk_525, hk_527, \
                         hk_529, hk_531 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_98[k] = -f_489 * hk_0[k]
                  - f_490 * hk_3[k]
                  + f_491 * hk_5[k]
                  - f_490 * hk_10[k]
                  + f_492 * hk_12[k]
                  - f_492 * hk_14[k]
                  - f_489 * hk_21[k]
                  + f_491 * hk_23[k]
                  - f_492 * hk_25[k]
                  + f_493 * hk_27[k]
                  - f_494 * hk_108[k]
                  - f_495 * hk_111[k]
                  + f_492 * hk_113[k]
                  - f_495 * hk_118[k]
                  + f_496 * hk_120[k]
                  - f_496 * hk_122[k]
                  - f_494 * hk_129[k]
                  + f_492 * hk_131[k]
                  - f_496 * hk_133[k]
                  + f_384 * hk_135[k]
                  + f_497 * hk_180[k]
                  + f_498 * hk_183[k]
                  - f_499 * hk_185[k]
                  + f_498 * hk_190[k]
                  - f_500 * hk_192[k]
                  + f_500 * hk_194[k]
                  + f_497 * hk_201[k]
                  - f_499 * hk_203[k]
                  + f_500 * hk_205[k]
                  - f_501 * hk_207[k]
                  - f_489 * hk_360[k]
                  - f_490 * hk_363[k]
                  + f_491 * hk_365[k]
                  - f_490 * hk_370[k]
                  + f_492 * hk_372[k]
                  - f_492 * hk_374[k]
                  - f_489 * hk_381[k]
                  + f_491 * hk_383[k]
                  - f_492 * hk_385[k]
                  + f_493 * hk_387[k]
                  + f_497 * hk_432[k]
                  + f_498 * hk_435[k]
                  - f_499 * hk_437[k]
                  + f_498 * hk_442[k]
                  - f_500 * hk_444[k]
                  + f_500 * hk_446[k]
                  + f_497 * hk_453[k]
                  - f_499 * hk_455[k]
                  + f_500 * hk_457[k]
                  - f_501 * hk_459[k]
                  - f_502 * hk_504[k]
                  - f_491 * hk_507[k]
                  + f_503 * hk_509[k]
                  - f_491 * hk_514[k]
                  + f_504 * hk_516[k]
                  - f_504 * hk_518[k]
                  - f_502 * hk_525[k]
                  + f_503 * hk_527[k]
                  - f_504 * hk_529[k]
                  + f_505 * hk_531[k];
    }

#pragma omp simd aligned(hk_2, hk_7, hk_9, hk_16, hk_20, hk_29, hk_31, hk_33, hk_110, hk_115, \
                         hk_117, hk_124, hk_128, hk_137, hk_139, hk_141, hk_182, hk_187, \
                         hk_189, hk_196, hk_200, hk_209, hk_211, hk_213, hk_362, hk_367, \
                         hk_369, hk_376, hk_380, hk_389, hk_391, hk_393, hk_434, hk_439, \
                         hk_441, hk_448, hk_452, hk_461, hk_463, hk_465, hk_506, hk_511, \
                         hk_513, hk_520, hk_524, hk_533, hk_535, \
                         hk_537 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_99[k] = f_512 * hk_2[k]
                  + f_512 * hk_7[k]
                  - f_513 * hk_9[k]
                  - f_512 * hk_16[k]
                  + f_514 * hk_20[k]
                  - f_512 * hk_29[k]
                  + f_513 * hk_31[k]
                  - f_514 * hk_33[k]
                  + f_474 * hk_110[k]
                  + f_474 * hk_115[k]
                  - f_476 * hk_117[k]
                  - f_474 * hk_124[k]
                  + f_477 * hk_128[k]
                  - f_474 * hk_137[k]
                  + f_476 * hk_139[k]
                  - f_477 * hk_141[k]
                  - f_515 * hk_182[k]
                  - f_515 * hk_187[k]
                  + f_516 * hk_189[k]
                  + f_515 * hk_196[k]
                  - f_517 * hk_200[k]
                  + f_515 * hk_209[k]
                  - f_516 * hk_211[k]
                  + f_517 * hk_213[k]
                  + f_512 * hk_362[k]
                  + f_512 * hk_367[k]
                  - f_513 * hk_369[k]
                  - f_512 * hk_376[k]
                  + f_514 * hk_380[k]
                  - f_512 * hk_389[k]
                  + f_513 * hk_391[k]
                  - f_514 * hk_393[k]
                  - f_515 * hk_434[k]
                  - f_515 * hk_439[k]
                  + f_516 * hk_441[k]
                  + f_515 * hk_448[k]
                  - f_517 * hk_452[k]
                  + f_515 * hk_461[k]
                  - f_516 * hk_463[k]
                  + f_517 * hk_465[k]
                  + f_478 * hk_506[k]
                  + f_478 * hk_511[k]
                  - f_518 * hk_513[k]
                  - f_478 * hk_520[k]
                  + f_519 * hk_524[k]
                  - f_478 * hk_533[k]
                  + f_518 * hk_535[k]
                  - f_519 * hk_537[k];
    }

#pragma omp simd aligned(hk_0, hk_3, hk_5, hk_10, hk_12, hk_14, hk_21, hk_23, hk_25, hk_108, \
                         hk_111, hk_113, hk_118, hk_120, hk_122, hk_129, hk_131, hk_133, \
                         hk_180, hk_183, hk_185, hk_190, hk_192, hk_194, hk_201, hk_203, \
                         hk_205, hk_360, hk_363, hk_365, hk_370, hk_372, hk_374, hk_381, \
                         hk_383, hk_385, hk_432, hk_435, hk_437, hk_442, hk_444, hk_446, \
                         hk_453, hk_455, hk_457, hk_504, hk_507, hk_509, hk_514, hk_516, \
                         hk_518, hk_525, hk_527, hk_529 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_100[k] = f_452 * hk_0[k]
                   - f_452 * hk_3[k]
                   - f_455 * hk_5[k]
                   - f_450 * hk_10[k]
                   + f_453 * hk_12[k]
                   + f_456 * hk_14[k]
                   - f_449 * hk_21[k]
                   + f_451 * hk_23[k]
                   - f_454 * hk_25[k]
                   + f_460 * hk_108[k]
                   - f_460 * hk_111[k]
                   - f_453 * hk_113[k]
                   - f_458 * hk_118[k]
                   + f_454 * hk_120[k]
                   + f_462 * hk_122[k]
                   - f_457 * hk_129[k]
                   + f_459 * hk_131[k]
                   - f_461 * hk_133[k]
                   - f_465 * hk_180[k]
                   + f_465 * hk_183[k]
                   + f_468 * hk_185[k]
                   + f_451 * hk_190[k]
                   - f_466 * hk_192[k]
                   - f_469 * hk_194[k]
                   + f_463 * hk_201[k]
                   - f_464 * hk_203[k]
                   + f_467 * hk_205[k]
                   + f_452 * hk_360[k]
                   - f_452 * hk_363[k]
                   - f_455 * hk_365[k]
                   - f_450 * hk_370[k]
                   + f_453 * hk_372[k]
                   + f_456 * hk_374[k]
                   - f_449 * hk_381[k]
                   + f_451 * hk_383[k]
                   - f_454 * hk_385[k]
                   - f_465 * hk_432[k]
                   + f_465 * hk_435[k]
                   + f_468 * hk_437[k]
                   + f_451 * hk_442[k]
                   - f_466 * hk_444[k]
                   - f_469 * hk_446[k]
                   + f_463 * hk_453[k]
                   - f_464 * hk_455[k]
                   + f_467 * hk_457[k]
                   + f_471 * hk_504[k]
                   - f_471 * hk_507[k]
                   - f_461 * hk_509[k]
                   - f_453 * hk_514[k]
                   + f_469 * hk_516[k]
                   + f_473 * hk_518[k]
                   - f_470 * hk_525[k]
                   + f_466 * hk_527[k]
                   - f_472 * hk_529[k];
    }

#pragma omp simd aligned(hk_2, hk_7, hk_9, hk_16, hk_18, hk_29, hk_31, hk_110, hk_115, hk_117, \
                         hk_124, hk_126, hk_137, hk_139, hk_182, hk_187, hk_189, hk_196, \
                         hk_198, hk_209, hk_211, hk_362, hk_367, hk_369, hk_376, hk_378, \
                         hk_389, hk_391, hk_434, hk_439, hk_441, hk_448, hk_450, hk_461, \
                         hk_463, hk_506, hk_511, hk_513, hk_520, hk_522, hk_533, \
                         hk_535 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_101[k] = -f_520 * hk_2[k]
                   + f_521 * hk_7[k]
                   + f_522 * hk_9[k]
                   + f_521 * hk_16[k]
                   - f_425 * hk_18[k]
                   - f_520 * hk_29[k]
                   + f_522 * hk_31[k]
                   - f_427 * hk_110[k]
                   + f_423 * hk_115[k]
                   + f_437 * hk_117[k]
                   + f_423 * hk_124[k]
                   - f_430 * hk_126[k]
                   - f_427 * hk_137[k]
                   + f_437 * hk_139[k]
                   + f_439 * hk_182[k]
                   - f_523 * hk_187[k]
                   - f_430 * hk_189[k]
                   - f_523 * hk_196[k]
                   + f_435 * hk_198[k]
                   + f_439 * hk_209[k]
                   - f_430 * hk_211[k]
                   - f_520 * hk_362[k]
                   + f_521 * hk_367[k]
                   + f_522 * hk_369[k]
                   + f_521 * hk_376[k]
                   - f_425 * hk_378[k]
                   - f_520 * hk_389[k]
                   + f_522 * hk_391[k]
                   + f_439 * hk_434[k]
                   - f_523 * hk_439[k]
                   - f_430 * hk_441[k]
                   - f_523 * hk_448[k]
                   + f_435 * hk_450[k]
                   + f_439 * hk_461[k]
                   - f_430 * hk_463[k]
                   - f_444 * hk_506[k]
                   + f_430 * hk_511[k]
                   + f_445 * hk_513[k]
                   + f_430 * hk_520[k]
                   - f_440 * hk_522[k]
                   - f_444 * hk_533[k]
                   + f_445 * hk_535[k];
    }

#pragma omp simd aligned(hk_0, hk_3, hk_5, hk_10, hk_12, hk_21, hk_23, hk_108, hk_111, hk_113, \
                         hk_118, hk_120, hk_129, hk_131, hk_180, hk_183, hk_185, hk_190, \
                         hk_192, hk_201, hk_203, hk_360, hk_363, hk_365, hk_370, hk_372, \
                         hk_381, hk_383, hk_432, hk_435, hk_437, hk_442, hk_444, hk_453, \
                         hk_455, hk_504, hk_507, hk_509, hk_514, hk_516, hk_525, \
                         hk_527 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_102[k] = -f_426 * hk_0[k]
                   + f_424 * hk_3[k]
                   + f_427 * hk_5[k]
                   + f_422 * hk_10[k]
                   - f_425 * hk_12[k]
                   - f_422 * hk_21[k]
                   + f_423 * hk_23[k]
                   - f_431 * hk_108[k]
                   + f_429 * hk_111[k]
                   + f_432 * hk_113[k]
                   + f_428 * hk_118[k]
                   - f_430 * hk_120[k]
                   - f_428 * hk_129[k]
                   + f_425 * hk_131[k]
                   + f_427 * hk_180[k]
                   - f_434 * hk_183[k]
                   - f_436 * hk_185[k]
                   - f_423 * hk_190[k]
                   + f_435 * hk_192[k]
                   + f_423 * hk_201[k]
                   - f_433 * hk_203[k]
                   - f_426 * hk_360[k]
                   + f_424 * hk_363[k]
                   + f_427 * hk_365[k]
                   + f_422 * hk_370[k]
                   - f_425 * hk_372[k]
                   - f_422 * hk_381[k]
                   + f_423 * hk_383[k]
                   + f_427 * hk_432[k]
                   - f_434 * hk_435[k]
                   - f_436 * hk_437[k]
                   - f_423 * hk_442[k]
                   + f_435 * hk_444[k]
                   + f_423 * hk_453[k]
                   - f_433 * hk_455[k]
                   - f_441 * hk_504[k]
                   + f_439 * hk_507[k]
                   + f_442 * hk_509[k]
                   + f_437 * hk_514[k]
                   - f_440 * hk_516[k]
                   - f_437 * hk_525[k]
                   + f_438 * hk_527[k];
    }

#pragma omp simd aligned(hk_2, hk_7, hk_16, hk_29, hk_110, hk_115, hk_124, hk_137, hk_182, \
                         hk_187, hk_196, hk_209, hk_362, hk_367, hk_376, hk_389, hk_434, \
                         hk_439, hk_448, hk_461, hk_506, hk_511, hk_520, \
                         hk_533 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_103[k] = f_524 * hk_2[k]
                   - f_525 * hk_7[k]
                   + f_525 * hk_16[k]
                   - f_524 * hk_29[k]
                   + f_526 * hk_110[k]
                   - f_527 * hk_115[k]
                   + f_527 * hk_124[k]
                   - f_526 * hk_137[k]
                   - f_416 * hk_182[k]
                   + f_528 * hk_187[k]
                   - f_528 * hk_196[k]
                   + f_416 * hk_209[k]
                   + f_524 * hk_362[k]
                   - f_525 * hk_367[k]
                   + f_525 * hk_376[k]
                   - f_524 * hk_389[k]
                   - f_416 * hk_434[k]
                   + f_528 * hk_439[k]
                   - f_528 * hk_448[k]
                   + f_416 * hk_461[k]
                   + f_529 * hk_506[k]
                   - f_530 * hk_511[k]
                   + f_530 * hk_520[k]
                   - f_529 * hk_533[k];
    }

#pragma omp simd aligned(hk_0, hk_3, hk_10, hk_21, hk_108, hk_111, hk_118, hk_129, hk_180, \
                         hk_183, hk_190, hk_201, hk_360, hk_363, hk_370, hk_381, hk_432, \
                         hk_435, hk_442, hk_453, hk_504, hk_507, hk_514, \
                         hk_525 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_104[k] = f_401 * hk_0[k]
                   - f_400 * hk_3[k]
                   + f_399 * hk_10[k]
                   - f_398 * hk_21[k]
                   + f_405 * hk_108[k]
                   - f_404 * hk_111[k]
                   + f_403 * hk_118[k]
                   - f_402 * hk_129[k]
                   - f_409 * hk_180[k]
                   + f_408 * hk_183[k]
                   - f_407 * hk_190[k]
                   + f_406 * hk_201[k]
                   + f_401 * hk_360[k]
                   - f_400 * hk_363[k]
                   + f_399 * hk_370[k]
                   - f_398 * hk_381[k]
                   - f_409 * hk_432[k]
                   + f_408 * hk_435[k]
                   - f_407 * hk_442[k]
                   + f_406 * hk_453[k]
                   + f_413 * hk_504[k]
                   - f_412 * hk_507[k]
                   + f_411 * hk_514[k]
                   - f_410 * hk_525[k];
    }

#pragma omp simd aligned(hk_73, hk_78, hk_87, hk_100, hk_325, hk_330, hk_339, hk_352, hk_577, \
                         hk_582, hk_591, hk_604, hk_649, hk_654, hk_663, \
                         hk_676 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_105[k] = -f_655 * hk_73[k]
                   + f_656 * hk_78[k]
                   - f_657 * hk_87[k]
                   + f_658 * hk_100[k]
                   + f_333 * hk_325[k]
                   - f_334 * hk_330[k]
                   + f_335 * hk_339[k]
                   - f_336 * hk_352[k]
                   + f_655 * hk_577[k]
                   - f_656 * hk_582[k]
                   + f_657 * hk_591[k]
                   - f_658 * hk_604[k]
                   - f_333 * hk_649[k]
                   + f_334 * hk_654[k]
                   - f_335 * hk_663[k]
                   + f_336 * hk_676[k];
    }

#pragma omp simd aligned(hk_76, hk_83, hk_94, hk_328, hk_335, hk_346, hk_580, hk_587, hk_598, \
                         hk_652, hk_659, hk_670 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_106[k] = -f_659 * hk_76[k]
                   + f_660 * hk_83[k]
                   - f_659 * hk_94[k]
                   + f_341 * hk_328[k]
                   - f_342 * hk_335[k]
                   + f_341 * hk_346[k]
                   + f_659 * hk_580[k]
                   - f_660 * hk_587[k]
                   + f_659 * hk_598[k]
                   - f_341 * hk_652[k]
                   + f_342 * hk_659[k]
                   - f_341 * hk_670[k];
    }

#pragma omp simd aligned(hk_73, hk_78, hk_80, hk_87, hk_89, hk_100, hk_102, hk_325, hk_330, \
                         hk_332, hk_339, hk_341, hk_352, hk_354, hk_577, hk_582, hk_584, \
                         hk_591, hk_593, hk_604, hk_606, hk_649, hk_654, hk_656, hk_663, \
                         hk_665, hk_676, hk_678 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_107[k] = f_661 * hk_73[k]
                   - f_661 * hk_78[k]
                   - f_391 * hk_80[k]
                   - f_662 * hk_87[k]
                   + f_346 * hk_89[k]
                   + f_663 * hk_100[k]
                   - f_390 * hk_102[k]
                   - f_345 * hk_325[k]
                   + f_345 * hk_330[k]
                   + f_346 * hk_332[k]
                   + f_347 * hk_339[k]
                   - f_348 * hk_341[k]
                   - f_349 * hk_352[k]
                   + f_350 * hk_354[k]
                   - f_661 * hk_577[k]
                   + f_661 * hk_582[k]
                   + f_391 * hk_584[k]
                   + f_662 * hk_591[k]
                   - f_346 * hk_593[k]
                   - f_663 * hk_604[k]
                   + f_390 * hk_606[k]
                   + f_345 * hk_649[k]
                   - f_345 * hk_654[k]
                   - f_346 * hk_656[k]
                   - f_347 * hk_663[k]
                   + f_348 * hk_665[k]
                   + f_349 * hk_676[k]
                   - f_350 * hk_678[k];
    }

#pragma omp simd aligned(hk_76, hk_85, hk_94, hk_96, hk_328, hk_337, hk_346, hk_348, hk_580, \
                         hk_589, hk_598, hk_600, hk_652, hk_661, hk_670, \
                         hk_672 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_108[k] = f_350 * hk_76[k]
                   - f_393 * hk_85[k]
                   - f_350 * hk_94[k]
                   + f_393 * hk_96[k]
                   - f_355 * hk_328[k]
                   + f_356 * hk_337[k]
                   + f_355 * hk_346[k]
                   - f_356 * hk_348[k]
                   - f_350 * hk_580[k]
                   + f_393 * hk_589[k]
                   + f_350 * hk_598[k]
                   - f_393 * hk_600[k]
                   + f_355 * hk_652[k]
                   - f_356 * hk_661[k]
                   - f_355 * hk_670[k]
                   + f_356 * hk_672[k];
    }

#pragma omp simd aligned(hk_73, hk_78, hk_80, hk_87, hk_89, hk_91, hk_100, hk_102, hk_104, \
                         hk_325, hk_330, hk_332, hk_339, hk_341, hk_343, hk_352, hk_354, \
                         hk_356, hk_577, hk_582, hk_584, hk_591, hk_593, hk_595, hk_604, \
                         hk_606, hk_608, hk_649, hk_654, hk_656, hk_663, hk_665, hk_667, \
                         hk_676, hk_678, hk_680 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_109[k] = -f_664 * hk_73[k]
                   - f_665 * hk_78[k]
                   + f_666 * hk_80[k]
                   - f_667 * hk_87[k]
                   + f_362 * hk_89[k]
                   - f_144 * hk_91[k]
                   + f_667 * hk_100[k]
                   - f_365 * hk_102[k]
                   + f_668 * hk_104[k]
                   + f_359 * hk_325[k]
                   + f_143 * hk_330[k]
                   - f_360 * hk_332[k]
                   + f_361 * hk_339[k]
                   - f_144 * hk_341[k]
                   + f_145 * hk_343[k]
                   - f_361 * hk_352[k]
                   + f_362 * hk_354[k]
                   - f_363 * hk_356[k]
                   + f_664 * hk_577[k]
                   + f_665 * hk_582[k]
                   - f_666 * hk_584[k]
                   + f_667 * hk_591[k]
                   - f_362 * hk_593[k]
                   + f_144 * hk_595[k]
                   - f_667 * hk_604[k]
                   + f_365 * hk_606[k]
                   - f_668 * hk_608[k]
                   - f_359 * hk_649[k]
                   - f_143 * hk_654[k]
                   + f_360 * hk_656[k]
                   - f_361 * hk_663[k]
                   + f_144 * hk_665[k]
                   - f_145 * hk_667[k]
                   + f_361 * hk_676[k]
                   - f_362 * hk_678[k]
                   + f_363 * hk_680[k];
    }

#pragma omp simd aligned(hk_76, hk_83, hk_85, hk_94, hk_96, hk_98, hk_328, hk_335, hk_337, \
                         hk_346, hk_348, hk_350, hk_580, hk_587, hk_589, hk_598, hk_600, \
                         hk_602, hk_652, hk_659, hk_661, hk_670, hk_672, \
                         hk_674 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_110[k] = -f_264 * hk_76[k]
                   - f_265 * hk_83[k]
                   + f_280 * hk_85[k]
                   - f_264 * hk_94[k]
                   + f_280 * hk_96[k]
                   - f_389 * hk_98[k]
                   + f_265 * hk_328[k]
                   + f_370 * hk_335[k]
                   - f_371 * hk_337[k]
                   + f_265 * hk_346[k]
                   - f_371 * hk_348[k]
                   + f_372 * hk_350[k]
                   + f_264 * hk_580[k]
                   + f_265 * hk_587[k]
                   - f_280 * hk_589[k]
                   + f_264 * hk_598[k]
                   - f_280 * hk_600[k]
                   + f_389 * hk_602[k]
                   - f_265 * hk_652[k]
                   - f_370 * hk_659[k]
                   + f_371 * hk_661[k]
                   - f_265 * hk_670[k]
                   + f_371 * hk_672[k]
                   - f_372 * hk_674[k];
    }

#pragma omp simd aligned(hk_73, hk_78, hk_80, hk_87, hk_89, hk_91, hk_100, hk_102, hk_104, \
                         hk_106, hk_325, hk_330, hk_332, hk_339, hk_341, hk_343, hk_352, \
                         hk_354, hk_356, hk_358, hk_577, hk_582, hk_584, hk_591, hk_593, \
                         hk_595, hk_604, hk_606, hk_608, hk_610, hk_649, hk_654, hk_656, \
                         hk_663, hk_665, hk_667, hk_676, hk_678, hk_680, \
                         hk_682 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_111[k] = f_506 * hk_73[k]
                   + f_256 * hk_78[k]
                   - f_259 * hk_80[k]
                   + f_256 * hk_87[k]
                   - f_248 * hk_89[k]
                   + f_248 * hk_91[k]
                   + f_506 * hk_100[k]
                   - f_259 * hk_102[k]
                   + f_248 * hk_104[k]
                   - f_669 * hk_106[k]
                   - f_375 * hk_325[k]
                   - f_250 * hk_330[k]
                   + f_248 * hk_332[k]
                   - f_250 * hk_339[k]
                   + f_137 * hk_341[k]
                   - f_137 * hk_343[k]
                   - f_375 * hk_352[k]
                   + f_248 * hk_354[k]
                   - f_137 * hk_356[k]
                   + f_376 * hk_358[k]
                   - f_506 * hk_577[k]
                   - f_256 * hk_582[k]
                   + f_259 * hk_584[k]
                   - f_256 * hk_591[k]
                   + f_248 * hk_593[k]
                   - f_248 * hk_595[k]
                   - f_506 * hk_604[k]
                   + f_259 * hk_606[k]
                   - f_248 * hk_608[k]
                   + f_669 * hk_610[k]
                   + f_375 * hk_649[k]
                   + f_250 * hk_654[k]
                   - f_248 * hk_656[k]
                   + f_250 * hk_663[k]
                   - f_137 * hk_665[k]
                   + f_137 * hk_667[k]
                   + f_375 * hk_676[k]
                   - f_248 * hk_678[k]
                   + f_137 * hk_680[k]
                   - f_376 * hk_682[k];
    }

#pragma omp simd aligned(hk_74, hk_79, hk_81, hk_88, hk_90, hk_92, hk_101, hk_103, hk_105, \
                         hk_107, hk_326, hk_331, hk_333, hk_340, hk_342, hk_344, hk_353, \
                         hk_355, hk_357, hk_359, hk_578, hk_583, hk_585, hk_592, hk_594, \
                         hk_596, hk_605, hk_607, hk_609, hk_611, hk_650, hk_655, hk_657, \
                         hk_664, hk_666, hk_668, hk_677, hk_679, hk_681, \
                         hk_683 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_112[k] = f_670 * hk_74[k]
                   + f_671 * hk_79[k]
                   - f_380 * hk_81[k]
                   + f_671 * hk_88[k]
                   - f_381 * hk_90[k]
                   + f_672 * hk_92[k]
                   + f_670 * hk_101[k]
                   - f_380 * hk_103[k]
                   + f_672 * hk_105[k]
                   - f_493 * hk_107[k]
                   - f_379 * hk_326[k]
                   - f_380 * hk_331[k]
                   + f_381 * hk_333[k]
                   - f_380 * hk_340[k]
                   + f_382 * hk_342[k]
                   - f_383 * hk_344[k]
                   - f_379 * hk_353[k]
                   + f_381 * hk_355[k]
                   - f_383 * hk_357[k]
                   + f_384 * hk_359[k]
                   - f_670 * hk_578[k]
                   - f_671 * hk_583[k]
                   + f_380 * hk_585[k]
                   - f_671 * hk_592[k]
                   + f_381 * hk_594[k]
                   - f_672 * hk_596[k]
                   - f_670 * hk_605[k]
                   + f_380 * hk_607[k]
                   - f_672 * hk_609[k]
                   + f_493 * hk_611[k]
                   + f_379 * hk_650[k]
                   + f_380 * hk_655[k]
                   - f_381 * hk_657[k]
                   + f_380 * hk_664[k]
                   - f_382 * hk_666[k]
                   + f_383 * hk_668[k]
                   + f_379 * hk_677[k]
                   - f_381 * hk_679[k]
                   + f_383 * hk_681[k]
                   - f_384 * hk_683[k];
    }

#pragma omp simd aligned(hk_72, hk_75, hk_77, hk_82, hk_84, hk_86, hk_93, hk_95, hk_97, hk_99, \
                         hk_324, hk_327, hk_329, hk_334, hk_336, hk_338, hk_345, hk_347, \
                         hk_349, hk_351, hk_576, hk_579, hk_581, hk_586, hk_588, hk_590, \
                         hk_597, hk_599, hk_601, hk_603, hk_648, hk_651, hk_653, hk_658, \
                         hk_660, hk_662, hk_669, hk_671, hk_673, \
                         hk_675 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_113[k] = f_506 * hk_72[k]
                   + f_256 * hk_75[k]
                   - f_259 * hk_77[k]
                   + f_256 * hk_82[k]
                   - f_248 * hk_84[k]
                   + f_248 * hk_86[k]
                   + f_506 * hk_93[k]
                   - f_259 * hk_95[k]
                   + f_248 * hk_97[k]
                   - f_669 * hk_99[k]
                   - f_375 * hk_324[k]
                   - f_250 * hk_327[k]
                   + f_248 * hk_329[k]
                   - f_250 * hk_334[k]
                   + f_137 * hk_336[k]
                   - f_137 * hk_338[k]
                   - f_375 * hk_345[k]
                   + f_248 * hk_347[k]
                   - f_137 * hk_349[k]
                   + f_376 * hk_351[k]
                   - f_506 * hk_576[k]
                   - f_256 * hk_579[k]
                   + f_259 * hk_581[k]
                   - f_256 * hk_586[k]
                   + f_248 * hk_588[k]
                   - f_248 * hk_590[k]
                   - f_506 * hk_597[k]
                   + f_259 * hk_599[k]
                   - f_248 * hk_601[k]
                   + f_669 * hk_603[k]
                   + f_375 * hk_648[k]
                   + f_250 * hk_651[k]
                   - f_248 * hk_653[k]
                   + f_250 * hk_658[k]
                   - f_137 * hk_660[k]
                   + f_137 * hk_662[k]
                   + f_375 * hk_669[k]
                   - f_248 * hk_671[k]
                   + f_137 * hk_673[k]
                   - f_376 * hk_675[k];
    }

#pragma omp simd aligned(hk_74, hk_79, hk_81, hk_88, hk_92, hk_101, hk_103, hk_105, hk_326, \
                         hk_331, hk_333, hk_340, hk_344, hk_353, hk_355, hk_357, hk_578, \
                         hk_583, hk_585, hk_592, hk_596, hk_605, hk_607, hk_609, hk_650, \
                         hk_655, hk_657, hk_664, hk_668, hk_677, hk_679, \
                         hk_681 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_114[k] = -f_673 * hk_74[k]
                   - f_673 * hk_79[k]
                   + f_279 * hk_81[k]
                   + f_673 * hk_88[k]
                   - f_674 * hk_92[k]
                   + f_673 * hk_101[k]
                   - f_279 * hk_103[k]
                   + f_674 * hk_105[k]
                   + f_264 * hk_326[k]
                   + f_264 * hk_331[k]
                   - f_280 * hk_333[k]
                   - f_264 * hk_340[k]
                   + f_389 * hk_344[k]
                   - f_264 * hk_353[k]
                   + f_280 * hk_355[k]
                   - f_389 * hk_357[k]
                   + f_673 * hk_578[k]
                   + f_673 * hk_583[k]
                   - f_279 * hk_585[k]
                   - f_673 * hk_592[k]
                   + f_674 * hk_596[k]
                   - f_673 * hk_605[k]
                   + f_279 * hk_607[k]
                   - f_674 * hk_609[k]
                   - f_264 * hk_650[k]
                   - f_264 * hk_655[k]
                   + f_280 * hk_657[k]
                   + f_264 * hk_664[k]
                   - f_389 * hk_668[k]
                   + f_264 * hk_677[k]
                   - f_280 * hk_679[k]
                   + f_389 * hk_681[k];
    }

#pragma omp simd aligned(hk_72, hk_75, hk_77, hk_82, hk_84, hk_86, hk_93, hk_95, hk_97, \
                         hk_324, hk_327, hk_329, hk_334, hk_336, hk_338, hk_345, hk_347, \
                         hk_349, hk_576, hk_579, hk_581, hk_586, hk_588, hk_590, hk_597, \
                         hk_599, hk_601, hk_648, hk_651, hk_653, hk_658, hk_660, hk_662, \
                         hk_669, hk_671, hk_673 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_115[k] = -f_667 * hk_72[k]
                   + f_667 * hk_75[k]
                   + f_365 * hk_77[k]
                   + f_665 * hk_82[k]
                   - f_362 * hk_84[k]
                   - f_668 * hk_86[k]
                   + f_664 * hk_93[k]
                   - f_666 * hk_95[k]
                   + f_144 * hk_97[k]
                   + f_361 * hk_324[k]
                   - f_361 * hk_327[k]
                   - f_362 * hk_329[k]
                   - f_143 * hk_334[k]
                   + f_144 * hk_336[k]
                   + f_363 * hk_338[k]
                   - f_359 * hk_345[k]
                   + f_360 * hk_347[k]
                   - f_145 * hk_349[k]
                   + f_667 * hk_576[k]
                   - f_667 * hk_579[k]
                   - f_365 * hk_581[k]
                   - f_665 * hk_586[k]
                   + f_362 * hk_588[k]
                   + f_668 * hk_590[k]
                   - f_664 * hk_597[k]
                   + f_666 * hk_599[k]
                   - f_144 * hk_601[k]
                   - f_361 * hk_648[k]
                   + f_361 * hk_651[k]
                   + f_362 * hk_653[k]
                   + f_143 * hk_658[k]
                   - f_144 * hk_660[k]
                   - f_363 * hk_662[k]
                   + f_359 * hk_669[k]
                   - f_360 * hk_671[k]
                   + f_145 * hk_673[k];
    }

#pragma omp simd aligned(hk_74, hk_79, hk_81, hk_88, hk_90, hk_101, hk_103, hk_326, hk_331, \
                         hk_333, hk_340, hk_342, hk_353, hk_355, hk_578, hk_583, hk_585, \
                         hk_592, hk_594, hk_605, hk_607, hk_650, hk_655, hk_657, hk_664, \
                         hk_666, hk_677, hk_679 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_116[k] = f_675 * hk_74[k]
                   - f_676 * hk_79[k]
                   - f_351 * hk_81[k]
                   - f_676 * hk_88[k]
                   + f_346 * hk_90[k]
                   + f_675 * hk_101[k]
                   - f_351 * hk_103[k]
                   - f_390 * hk_326[k]
                   + f_391 * hk_331[k]
                   + f_392 * hk_333[k]
                   + f_391 * hk_340[k]
                   - f_348 * hk_342[k]
                   - f_390 * hk_353[k]
                   + f_392 * hk_355[k]
                   - f_675 * hk_578[k]
                   + f_676 * hk_583[k]
                   + f_351 * hk_585[k]
                   + f_676 * hk_592[k]
                   - f_346 * hk_594[k]
                   - f_675 * hk_605[k]
                   + f_351 * hk_607[k]
                   + f_390 * hk_650[k]
                   - f_391 * hk_655[k]
                   - f_392 * hk_657[k]
                   - f_391 * hk_664[k]
                   + f_348 * hk_666[k]
                   + f_390 * hk_677[k]
                   - f_392 * hk_679[k];
    }

#pragma omp simd aligned(hk_72, hk_75, hk_77, hk_82, hk_84, hk_93, hk_95, hk_324, hk_327, \
                         hk_329, hk_334, hk_336, hk_345, hk_347, hk_576, hk_579, hk_581, \
                         hk_586, hk_588, hk_597, hk_599, hk_648, hk_651, hk_653, hk_658, \
                         hk_660, hk_669, hk_671 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_117[k] = f_663 * hk_72[k]
                   - f_662 * hk_75[k]
                   - f_390 * hk_77[k]
                   - f_661 * hk_82[k]
                   + f_346 * hk_84[k]
                   + f_661 * hk_93[k]
                   - f_391 * hk_95[k]
                   - f_349 * hk_324[k]
                   + f_347 * hk_327[k]
                   + f_350 * hk_329[k]
                   + f_345 * hk_334[k]
                   - f_348 * hk_336[k]
                   - f_345 * hk_345[k]
                   + f_346 * hk_347[k]
                   - f_663 * hk_576[k]
                   + f_662 * hk_579[k]
                   + f_390 * hk_581[k]
                   + f_661 * hk_586[k]
                   - f_346 * hk_588[k]
                   - f_661 * hk_597[k]
                   + f_391 * hk_599[k]
                   + f_349 * hk_648[k]
                   - f_347 * hk_651[k]
                   - f_350 * hk_653[k]
                   - f_345 * hk_658[k]
                   + f_348 * hk_660[k]
                   + f_345 * hk_669[k]
                   - f_346 * hk_671[k];
    }

#pragma omp simd aligned(hk_74, hk_79, hk_88, hk_101, hk_326, hk_331, hk_340, hk_353, hk_578, \
                         hk_583, hk_592, hk_605, hk_650, hk_655, hk_664, \
                         hk_677 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_118[k] = -f_677 * hk_74[k]
                   + f_678 * hk_79[k]
                   - f_678 * hk_88[k]
                   + f_677 * hk_101[k]
                   + f_394 * hk_326[k]
                   - f_395 * hk_331[k]
                   + f_395 * hk_340[k]
                   - f_394 * hk_353[k]
                   + f_677 * hk_578[k]
                   - f_678 * hk_583[k]
                   + f_678 * hk_592[k]
                   - f_677 * hk_605[k]
                   - f_394 * hk_650[k]
                   + f_395 * hk_655[k]
                   - f_395 * hk_664[k]
                   + f_394 * hk_677[k];
    }

#pragma omp simd aligned(hk_72, hk_75, hk_82, hk_93, hk_324, hk_327, hk_334, hk_345, hk_576, \
                         hk_579, hk_586, hk_597, hk_648, hk_651, hk_658, \
                         hk_669 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_119[k] = -f_658 * hk_72[k]
                   + f_657 * hk_75[k]
                   - f_656 * hk_82[k]
                   + f_655 * hk_93[k]
                   + f_336 * hk_324[k]
                   - f_335 * hk_327[k]
                   + f_334 * hk_334[k]
                   - f_333 * hk_345[k]
                   + f_658 * hk_576[k]
                   - f_657 * hk_579[k]
                   + f_656 * hk_586[k]
                   - f_655 * hk_597[k]
                   - f_336 * hk_648[k]
                   + f_335 * hk_651[k]
                   - f_334 * hk_658[k]
                   + f_333 * hk_669[k];
    }

#pragma omp simd aligned(hk_1, hk_6, hk_15, hk_28, hk_109, hk_114, hk_123, hk_136, hk_181, \
                         hk_186, hk_195, hk_208, hk_361, hk_366, hk_375, hk_388, hk_433, \
                         hk_438, hk_447, hk_460 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_120[k] = -f_173 * hk_1[k]
                   + f_174 * hk_6[k]
                   - f_161 * hk_15[k]
                   + f_175 * hk_28[k]
                   + f_165 * hk_109[k]
                   - f_166 * hk_114[k]
                   + f_167 * hk_123[k]
                   - f_168 * hk_136[k]
                   + f_176 * hk_181[k]
                   - f_177 * hk_186[k]
                   + f_169 * hk_195[k]
                   - f_178 * hk_208[k]
                   + f_161 * hk_361[k]
                   - f_162 * hk_366[k]
                   + f_163 * hk_375[k]
                   - f_164 * hk_388[k]
                   - f_169 * hk_433[k]
                   + f_170 * hk_438[k]
                   - f_171 * hk_447[k]
                   + f_172 * hk_460[k];
    }

#pragma omp simd aligned(hk_4, hk_11, hk_22, hk_112, hk_119, hk_130, hk_184, hk_191, hk_202, \
                         hk_364, hk_371, hk_382, hk_436, hk_443, \
                         hk_454 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_121[k] = -f_185 * hk_4[k]
                   + f_186 * hk_11[k]
                   - f_185 * hk_22[k]
                   + f_181 * hk_112[k]
                   - f_182 * hk_119[k]
                   + f_181 * hk_130[k]
                   + f_187 * hk_184[k]
                   - f_188 * hk_191[k]
                   + f_187 * hk_202[k]
                   + f_179 * hk_364[k]
                   - f_180 * hk_371[k]
                   + f_179 * hk_382[k]
                   - f_183 * hk_436[k]
                   + f_184 * hk_443[k]
                   - f_183 * hk_454[k];
    }

#pragma omp simd aligned(hk_1, hk_6, hk_8, hk_15, hk_17, hk_28, hk_30, hk_109, hk_114, hk_116, \
                         hk_123, hk_125, hk_136, hk_138, hk_181, hk_186, hk_188, hk_195, \
                         hk_197, hk_208, hk_210, hk_361, hk_366, hk_368, hk_375, hk_377, \
                         hk_388, hk_390, hk_433, hk_438, hk_440, hk_447, hk_449, hk_460, \
                         hk_462 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_122[k] = f_205 * hk_1[k]
                   - f_205 * hk_6[k]
                   - f_206 * hk_8[k]
                   - f_207 * hk_15[k]
                   + f_196 * hk_17[k]
                   + f_208 * hk_28[k]
                   - f_209 * hk_30[k]
                   - f_195 * hk_109[k]
                   + f_195 * hk_114[k]
                   + f_196 * hk_116[k]
                   + f_197 * hk_123[k]
                   - f_198 * hk_125[k]
                   - f_199 * hk_136[k]
                   + f_200 * hk_138[k]
                   - f_210 * hk_181[k]
                   + f_210 * hk_186[k]
                   + f_211 * hk_188[k]
                   + f_212 * hk_195[k]
                   - f_213 * hk_197[k]
                   - f_214 * hk_208[k]
                   + f_215 * hk_210[k]
                   - f_189 * hk_361[k]
                   + f_189 * hk_366[k]
                   + f_190 * hk_368[k]
                   + f_191 * hk_375[k]
                   - f_192 * hk_377[k]
                   - f_193 * hk_388[k]
                   + f_194 * hk_390[k]
                   + f_196 * hk_433[k]
                   - f_196 * hk_438[k]
                   - f_201 * hk_440[k]
                   - f_202 * hk_447[k]
                   + f_203 * hk_449[k]
                   + f_200 * hk_460[k]
                   - f_204 * hk_462[k];
    }

#pragma omp simd aligned(hk_4, hk_13, hk_22, hk_24, hk_112, hk_121, hk_130, hk_132, hk_184, \
                         hk_193, hk_202, hk_204, hk_364, hk_373, hk_382, hk_384, hk_436, \
                         hk_445, hk_454, hk_456 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_123[k] = f_200 * hk_4[k]
                   - f_220 * hk_13[k]
                   - f_200 * hk_22[k]
                   + f_220 * hk_24[k]
                   - f_216 * hk_112[k]
                   + f_217 * hk_121[k]
                   + f_216 * hk_130[k]
                   - f_217 * hk_132[k]
                   - f_221 * hk_184[k]
                   + f_222 * hk_193[k]
                   + f_221 * hk_202[k]
                   - f_222 * hk_204[k]
                   - f_212 * hk_364[k]
                   + f_198 * hk_373[k]
                   + f_212 * hk_382[k]
                   - f_198 * hk_384[k]
                   + f_218 * hk_436[k]
                   - f_219 * hk_445[k]
                   - f_218 * hk_454[k]
                   + f_219 * hk_456[k];
    }

#pragma omp simd aligned(hk_1, hk_6, hk_8, hk_15, hk_17, hk_19, hk_28, hk_30, hk_32, hk_109, \
                         hk_114, hk_116, hk_123, hk_125, hk_127, hk_136, hk_138, hk_140, \
                         hk_181, hk_186, hk_188, hk_195, hk_197, hk_199, hk_208, hk_210, \
                         hk_212, hk_361, hk_366, hk_368, hk_375, hk_377, hk_379, hk_388, \
                         hk_390, hk_392, hk_433, hk_438, hk_440, hk_447, hk_449, hk_451, \
                         hk_460, hk_462, hk_464 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_124[k] = -f_226 * hk_1[k]
                   - f_240 * hk_6[k]
                   + f_227 * hk_8[k]
                   - f_241 * hk_15[k]
                   + f_233 * hk_17[k]
                   - f_228 * hk_19[k]
                   + f_241 * hk_28[k]
                   - f_242 * hk_30[k]
                   + f_243 * hk_32[k]
                   + f_229 * hk_109[k]
                   + f_230 * hk_114[k]
                   - f_153 * hk_116[k]
                   + f_231 * hk_123[k]
                   - f_228 * hk_125[k]
                   + f_232 * hk_127[k]
                   - f_231 * hk_136[k]
                   + f_233 * hk_138[k]
                   - f_234 * hk_140[k]
                   + f_237 * hk_181[k]
                   + f_233 * hk_186[k]
                   - f_139 * hk_188[k]
                   + f_244 * hk_195[k]
                   - f_245 * hk_197[k]
                   + f_154 * hk_199[k]
                   - f_244 * hk_208[k]
                   + f_232 * hk_210[k]
                   - f_246 * hk_212[k]
                   + f_223 * hk_361[k]
                   + f_224 * hk_366[k]
                   - f_225 * hk_368[k]
                   + f_226 * hk_375[k]
                   - f_153 * hk_377[k]
                   + f_138 * hk_379[k]
                   - f_226 * hk_388[k]
                   + f_227 * hk_390[k]
                   - f_228 * hk_392[k]
                   - f_235 * hk_433[k]
                   - f_153 * hk_438[k]
                   + f_236 * hk_440[k]
                   - f_237 * hk_447[k]
                   + f_238 * hk_449[k]
                   - f_239 * hk_451[k]
                   + f_237 * hk_460[k]
                   - f_139 * hk_462[k]
                   + f_154 * hk_464[k];
    }

#pragma omp simd aligned(hk_4, hk_11, hk_13, hk_22, hk_24, hk_26, hk_112, hk_119, hk_121, \
                         hk_130, hk_132, hk_134, hk_184, hk_191, hk_193, hk_202, hk_204, \
                         hk_206, hk_364, hk_371, hk_373, hk_382, hk_384, hk_386, hk_436, \
                         hk_443, hk_445, hk_454, hk_456, hk_458 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_125[k] = -f_256 * hk_4[k]
                   - f_250 * hk_11[k]
                   + f_257 * hk_13[k]
                   - f_256 * hk_22[k]
                   + f_257 * hk_24[k]
                   - f_258 * hk_26[k]
                   + f_250 * hk_112[k]
                   + f_251 * hk_119[k]
                   - f_252 * hk_121[k]
                   + f_250 * hk_130[k]
                   - f_252 * hk_132[k]
                   + f_253 * hk_134[k]
                   + f_259 * hk_184[k]
                   + f_248 * hk_191[k]
                   - f_260 * hk_193[k]
                   + f_259 * hk_202[k]
                   - f_260 * hk_204[k]
                   + f_261 * hk_206[k]
                   + f_247 * hk_364[k]
                   + f_131 * hk_371[k]
                   - f_248 * hk_373[k]
                   + f_247 * hk_382[k]
                   - f_248 * hk_384[k]
                   + f_249 * hk_386[k]
                   - f_136 * hk_436[k]
                   - f_134 * hk_443[k]
                   + f_254 * hk_445[k]
                   - f_136 * hk_454[k]
                   + f_254 * hk_456[k]
                   - f_255 * hk_458[k];
    }

#pragma omp simd aligned(hk_1, hk_6, hk_8, hk_15, hk_17, hk_19, hk_28, hk_30, hk_32, hk_34, \
                         hk_109, hk_114, hk_116, hk_123, hk_125, hk_127, hk_136, hk_138, \
                         hk_140, hk_142, hk_181, hk_186, hk_188, hk_195, hk_197, hk_199, \
                         hk_208, hk_210, hk_212, hk_214, hk_361, hk_366, hk_368, hk_375, \
                         hk_377, hk_379, hk_388, hk_390, hk_392, hk_394, hk_433, hk_438, \
                         hk_440, hk_447, hk_449, hk_451, hk_460, hk_462, hk_464, \
                         hk_466 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_126[k] = f_276 * hk_1[k]
                   + f_262 * hk_6[k]
                   - f_272 * hk_8[k]
                   + f_262 * hk_15[k]
                   - f_269 * hk_17[k]
                   + f_269 * hk_19[k]
                   + f_276 * hk_28[k]
                   - f_272 * hk_30[k]
                   + f_269 * hk_32[k]
                   - f_277 * hk_34[k]
                   - f_267 * hk_109[k]
                   - f_268 * hk_114[k]
                   + f_269 * hk_116[k]
                   - f_268 * hk_123[k]
                   + f_270 * hk_125[k]
                   - f_270 * hk_127[k]
                   - f_267 * hk_136[k]
                   + f_269 * hk_138[k]
                   - f_270 * hk_140[k]
                   + f_271 * hk_142[k]
                   - f_278 * hk_181[k]
                   - f_272 * hk_186[k]
                   + f_279 * hk_188[k]
                   - f_272 * hk_195[k]
                   + f_280 * hk_197[k]
                   - f_280 * hk_199[k]
                   - f_278 * hk_208[k]
                   + f_279 * hk_210[k]
                   - f_280 * hk_212[k]
                   + f_281 * hk_214[k]
                   - f_262 * hk_361[k]
                   - f_263 * hk_366[k]
                   + f_264 * hk_368[k]
                   - f_263 * hk_375[k]
                   + f_265 * hk_377[k]
                   - f_265 * hk_379[k]
                   - f_262 * hk_388[k]
                   + f_264 * hk_390[k]
                   - f_265 * hk_392[k]
                   + f_266 * hk_394[k]
                   + f_272 * hk_433[k]
                   + f_264 * hk_438[k]
                   - f_273 * hk_440[k]
                   + f_264 * hk_447[k]
                   - f_274 * hk_449[k]
                   + f_274 * hk_451[k]
                   + f_272 * hk_460[k]
                   - f_273 * hk_462[k]
                   + f_274 * hk_464[k]
                   - f_275 * hk_466[k];
    }

#pragma omp simd aligned(hk_2, hk_7, hk_9, hk_16, hk_18, hk_20, hk_29, hk_31, hk_33, hk_35, \
                         hk_110, hk_115, hk_117, hk_124, hk_126, hk_128, hk_137, hk_139, \
                         hk_141, hk_143, hk_182, hk_187, hk_189, hk_196, hk_198, hk_200, \
                         hk_209, hk_211, hk_213, hk_215, hk_362, hk_367, hk_369, hk_376, \
                         hk_378, hk_380, hk_389, hk_391, hk_393, hk_395, hk_434, hk_439, \
                         hk_441, hk_448, hk_450, hk_452, hk_461, hk_463, hk_465, \
                         hk_467 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_127[k] = f_299 * hk_2[k]
                   + f_282 * hk_7[k]
                   - f_289 * hk_9[k]
                   + f_282 * hk_16[k]
                   - f_290 * hk_18[k]
                   + f_300 * hk_20[k]
                   + f_299 * hk_29[k]
                   - f_289 * hk_31[k]
                   + f_300 * hk_33[k]
                   - f_301 * hk_35[k]
                   - f_288 * hk_110[k]
                   - f_289 * hk_115[k]
                   + f_290 * hk_117[k]
                   - f_289 * hk_124[k]
                   + f_291 * hk_126[k]
                   - f_292 * hk_128[k]
                   - f_288 * hk_137[k]
                   + f_290 * hk_139[k]
                   - f_292 * hk_141[k]
                   + f_293 * hk_143[k]
                   - f_302 * hk_182[k]
                   - f_291 * hk_187[k]
                   + f_303 * hk_189[k]
                   - f_291 * hk_196[k]
                   + f_304 * hk_198[k]
                   - f_305 * hk_200[k]
                   - f_302 * hk_209[k]
                   + f_303 * hk_211[k]
                   - f_305 * hk_213[k]
                   + f_306 * hk_215[k]
                   - f_282 * hk_362[k]
                   - f_283 * hk_367[k]
                   + f_284 * hk_369[k]
                   - f_283 * hk_376[k]
                   + f_285 * hk_378[k]
                   - f_286 * hk_380[k]
                   - f_282 * hk_389[k]
                   + f_284 * hk_391[k]
                   - f_286 * hk_393[k]
                   + f_287 * hk_395[k]
                   + f_291 * hk_434[k]
                   + f_294 * hk_439[k]
                   - f_295 * hk_441[k]
                   + f_294 * hk_448[k]
                   - f_296 * hk_450[k]
                   + f_297 * hk_452[k]
                   + f_291 * hk_461[k]
                   - f_295 * hk_463[k]
                   + f_297 * hk_465[k]
                   - f_298 * hk_467[k];
    }

#pragma omp simd aligned(hk_0, hk_3, hk_5, hk_10, hk_12, hk_14, hk_21, hk_23, hk_25, hk_27, \
                         hk_108, hk_111, hk_113, hk_118, hk_120, hk_122, hk_129, hk_131, \
                         hk_133, hk_135, hk_180, hk_183, hk_185, hk_190, hk_192, hk_194, \
                         hk_201, hk_203, hk_205, hk_207, hk_360, hk_363, hk_365, hk_370, \
                         hk_372, hk_374, hk_381, hk_383, hk_385, hk_387, hk_432, hk_435, \
                         hk_437, hk_442, hk_444, hk_446, hk_453, hk_455, hk_457, \
                         hk_459 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_128[k] = f_276 * hk_0[k]
                   + f_262 * hk_3[k]
                   - f_272 * hk_5[k]
                   + f_262 * hk_10[k]
                   - f_269 * hk_12[k]
                   + f_269 * hk_14[k]
                   + f_276 * hk_21[k]
                   - f_272 * hk_23[k]
                   + f_269 * hk_25[k]
                   - f_277 * hk_27[k]
                   - f_267 * hk_108[k]
                   - f_268 * hk_111[k]
                   + f_269 * hk_113[k]
                   - f_268 * hk_118[k]
                   + f_270 * hk_120[k]
                   - f_270 * hk_122[k]
                   - f_267 * hk_129[k]
                   + f_269 * hk_131[k]
                   - f_270 * hk_133[k]
                   + f_271 * hk_135[k]
                   - f_278 * hk_180[k]
                   - f_272 * hk_183[k]
                   + f_279 * hk_185[k]
                   - f_272 * hk_190[k]
                   + f_280 * hk_192[k]
                   - f_280 * hk_194[k]
                   - f_278 * hk_201[k]
                   + f_279 * hk_203[k]
                   - f_280 * hk_205[k]
                   + f_281 * hk_207[k]
                   - f_262 * hk_360[k]
                   - f_263 * hk_363[k]
                   + f_264 * hk_365[k]
                   - f_263 * hk_370[k]
                   + f_265 * hk_372[k]
                   - f_265 * hk_374[k]
                   - f_262 * hk_381[k]
                   + f_264 * hk_383[k]
                   - f_265 * hk_385[k]
                   + f_266 * hk_387[k]
                   + f_272 * hk_432[k]
                   + f_264 * hk_435[k]
                   - f_273 * hk_437[k]
                   + f_264 * hk_442[k]
                   - f_274 * hk_444[k]
                   + f_274 * hk_446[k]
                   + f_272 * hk_453[k]
                   - f_273 * hk_455[k]
                   + f_274 * hk_457[k]
                   - f_275 * hk_459[k];
    }

#pragma omp simd aligned(hk_2, hk_7, hk_9, hk_16, hk_20, hk_29, hk_31, hk_33, hk_110, hk_115, \
                         hk_117, hk_124, hk_128, hk_137, hk_139, hk_141, hk_182, hk_187, \
                         hk_189, hk_196, hk_200, hk_209, hk_211, hk_213, hk_362, hk_367, \
                         hk_369, hk_376, hk_380, hk_389, hk_391, hk_393, hk_434, hk_439, \
                         hk_441, hk_448, hk_452, hk_461, hk_463, \
                         hk_465 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_129[k] = -f_312 * hk_2[k]
                   - f_312 * hk_7[k]
                   + f_313 * hk_9[k]
                   + f_312 * hk_16[k]
                   - f_314 * hk_20[k]
                   + f_312 * hk_29[k]
                   - f_313 * hk_31[k]
                   + f_314 * hk_33[k]
                   + f_256 * hk_110[k]
                   + f_256 * hk_115[k]
                   - f_257 * hk_117[k]
                   - f_256 * hk_124[k]
                   + f_258 * hk_128[k]
                   - f_256 * hk_137[k]
                   + f_257 * hk_139[k]
                   - f_258 * hk_141[k]
                   + f_251 * hk_182[k]
                   + f_251 * hk_187[k]
                   - f_315 * hk_189[k]
                   - f_251 * hk_196[k]
                   + f_316 * hk_200[k]
                   - f_251 * hk_209[k]
                   + f_315 * hk_211[k]
                   - f_316 * hk_213[k]
                   + f_307 * hk_362[k]
                   + f_307 * hk_367[k]
                   - f_259 * hk_369[k]
                   - f_307 * hk_376[k]
                   + f_308 * hk_380[k]
                   - f_307 * hk_389[k]
                   + f_259 * hk_391[k]
                   - f_308 * hk_393[k]
                   - f_309 * hk_434[k]
                   - f_309 * hk_439[k]
                   + f_310 * hk_441[k]
                   + f_309 * hk_448[k]
                   - f_311 * hk_452[k]
                   + f_309 * hk_461[k]
                   - f_310 * hk_463[k]
                   + f_311 * hk_465[k];
    }

#pragma omp simd aligned(hk_0, hk_3, hk_5, hk_10, hk_12, hk_14, hk_21, hk_23, hk_25, hk_108, \
                         hk_111, hk_113, hk_118, hk_120, hk_122, hk_129, hk_131, hk_133, \
                         hk_180, hk_183, hk_185, hk_190, hk_192, hk_194, hk_201, hk_203, \
                         hk_205, hk_360, hk_363, hk_365, hk_370, hk_372, hk_374, hk_381, \
                         hk_383, hk_385, hk_432, hk_435, hk_437, hk_442, hk_444, hk_446, \
                         hk_453, hk_455, hk_457 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_130[k] = -f_241 * hk_0[k]
                   + f_241 * hk_3[k]
                   + f_242 * hk_5[k]
                   + f_240 * hk_10[k]
                   - f_233 * hk_12[k]
                   - f_243 * hk_14[k]
                   + f_226 * hk_21[k]
                   - f_227 * hk_23[k]
                   + f_228 * hk_25[k]
                   + f_231 * hk_108[k]
                   - f_231 * hk_111[k]
                   - f_233 * hk_113[k]
                   - f_230 * hk_118[k]
                   + f_228 * hk_120[k]
                   + f_234 * hk_122[k]
                   - f_229 * hk_129[k]
                   + f_153 * hk_131[k]
                   - f_232 * hk_133[k]
                   + f_244 * hk_180[k]
                   - f_244 * hk_183[k]
                   - f_232 * hk_185[k]
                   - f_233 * hk_190[k]
                   + f_245 * hk_192[k]
                   + f_246 * hk_194[k]
                   - f_237 * hk_201[k]
                   + f_139 * hk_203[k]
                   - f_154 * hk_205[k]
                   + f_226 * hk_360[k]
                   - f_226 * hk_363[k]
                   - f_227 * hk_365[k]
                   - f_224 * hk_370[k]
                   + f_153 * hk_372[k]
                   + f_228 * hk_374[k]
                   - f_223 * hk_381[k]
                   + f_225 * hk_383[k]
                   - f_138 * hk_385[k]
                   - f_237 * hk_432[k]
                   + f_237 * hk_435[k]
                   + f_139 * hk_437[k]
                   + f_153 * hk_442[k]
                   - f_238 * hk_444[k]
                   - f_154 * hk_446[k]
                   + f_235 * hk_453[k]
                   - f_236 * hk_455[k]
                   + f_239 * hk_457[k];
    }

#pragma omp simd aligned(hk_2, hk_7, hk_9, hk_16, hk_18, hk_29, hk_31, hk_110, hk_115, hk_117, \
                         hk_124, hk_126, hk_137, hk_139, hk_182, hk_187, hk_189, hk_196, \
                         hk_198, hk_209, hk_211, hk_362, hk_367, hk_369, hk_376, hk_378, \
                         hk_389, hk_391, hk_434, hk_439, hk_441, hk_448, hk_450, hk_461, \
                         hk_463 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_131[k] = f_320 * hk_2[k]
                   - f_321 * hk_7[k]
                   - f_322 * hk_9[k]
                   - f_321 * hk_16[k]
                   + f_196 * hk_18[k]
                   + f_320 * hk_29[k]
                   - f_322 * hk_31[k]
                   - f_209 * hk_110[k]
                   + f_206 * hk_115[k]
                   + f_210 * hk_117[k]
                   + f_206 * hk_124[k]
                   - f_198 * hk_126[k]
                   - f_209 * hk_137[k]
                   + f_210 * hk_139[k]
                   - f_216 * hk_182[k]
                   + f_198 * hk_187[k]
                   + f_217 * hk_189[k]
                   + f_198 * hk_196[k]
                   - f_213 * hk_198[k]
                   - f_216 * hk_209[k]
                   + f_217 * hk_211[k]
                   - f_197 * hk_362[k]
                   + f_317 * hk_367[k]
                   + f_206 * hk_369[k]
                   + f_317 * hk_376[k]
                   - f_192 * hk_378[k]
                   - f_197 * hk_389[k]
                   + f_206 * hk_391[k]
                   + f_318 * hk_434[k]
                   - f_319 * hk_439[k]
                   - f_211 * hk_441[k]
                   - f_319 * hk_448[k]
                   + f_203 * hk_450[k]
                   + f_318 * hk_461[k]
                   - f_211 * hk_463[k];
    }

#pragma omp simd aligned(hk_0, hk_3, hk_5, hk_10, hk_12, hk_21, hk_23, hk_108, hk_111, hk_113, \
                         hk_118, hk_120, hk_129, hk_131, hk_180, hk_183, hk_185, hk_190, \
                         hk_192, hk_201, hk_203, hk_360, hk_363, hk_365, hk_370, hk_372, \
                         hk_381, hk_383, hk_432, hk_435, hk_437, hk_442, hk_444, hk_453, \
                         hk_455 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_132[k] = f_208 * hk_0[k]
                   - f_207 * hk_3[k]
                   - f_209 * hk_5[k]
                   - f_205 * hk_10[k]
                   + f_196 * hk_12[k]
                   + f_205 * hk_21[k]
                   - f_206 * hk_23[k]
                   - f_199 * hk_108[k]
                   + f_197 * hk_111[k]
                   + f_200 * hk_113[k]
                   + f_195 * hk_118[k]
                   - f_198 * hk_120[k]
                   - f_195 * hk_129[k]
                   + f_196 * hk_131[k]
                   - f_214 * hk_180[k]
                   + f_212 * hk_183[k]
                   + f_215 * hk_185[k]
                   + f_210 * hk_190[k]
                   - f_213 * hk_192[k]
                   - f_210 * hk_201[k]
                   + f_211 * hk_203[k]
                   - f_193 * hk_360[k]
                   + f_191 * hk_363[k]
                   + f_194 * hk_365[k]
                   + f_189 * hk_370[k]
                   - f_192 * hk_372[k]
                   - f_189 * hk_381[k]
                   + f_190 * hk_383[k]
                   + f_200 * hk_432[k]
                   - f_202 * hk_435[k]
                   - f_204 * hk_437[k]
                   - f_196 * hk_442[k]
                   + f_203 * hk_444[k]
                   + f_196 * hk_453[k]
                   - f_201 * hk_455[k];
    }

#pragma omp simd aligned(hk_2, hk_7, hk_16, hk_29, hk_110, hk_115, hk_124, hk_137, hk_182, \
                         hk_187, hk_196, hk_209, hk_362, hk_367, hk_376, hk_389, hk_434, \
                         hk_439, hk_448, hk_461 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_133[k] = -f_329 * hk_2[k]
                   + f_330 * hk_7[k]
                   - f_330 * hk_16[k]
                   + f_329 * hk_29[k]
                   + f_325 * hk_110[k]
                   - f_326 * hk_115[k]
                   + f_326 * hk_124[k]
                   - f_325 * hk_137[k]
                   + f_331 * hk_182[k]
                   - f_332 * hk_187[k]
                   + f_332 * hk_196[k]
                   - f_331 * hk_209[k]
                   + f_323 * hk_362[k]
                   - f_324 * hk_367[k]
                   + f_324 * hk_376[k]
                   - f_323 * hk_389[k]
                   - f_327 * hk_434[k]
                   + f_328 * hk_439[k]
                   - f_328 * hk_448[k]
                   + f_327 * hk_461[k];
    }

#pragma omp simd aligned(hk_0, hk_3, hk_10, hk_21, hk_108, hk_111, hk_118, hk_129, hk_180, \
                         hk_183, hk_190, hk_201, hk_360, hk_363, hk_370, hk_381, hk_432, \
                         hk_435, hk_442, hk_453 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_134[k] = -f_175 * hk_0[k]
                   + f_161 * hk_3[k]
                   - f_174 * hk_10[k]
                   + f_173 * hk_21[k]
                   + f_168 * hk_108[k]
                   - f_167 * hk_111[k]
                   + f_166 * hk_118[k]
                   - f_165 * hk_129[k]
                   + f_178 * hk_180[k]
                   - f_169 * hk_183[k]
                   + f_177 * hk_190[k]
                   - f_176 * hk_201[k]
                   + f_164 * hk_360[k]
                   - f_163 * hk_363[k]
                   + f_162 * hk_370[k]
                   - f_161 * hk_381[k]
                   - f_172 * hk_432[k]
                   + f_171 * hk_435[k]
                   - f_170 * hk_442[k]
                   + f_169 * hk_453[k];
    }

#pragma omp simd aligned(hk_73, hk_78, hk_87, hk_100, hk_253, hk_258, hk_267, hk_280, hk_577, \
                         hk_582, hk_591, hk_604 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_135[k] = f_679 * hk_73[k]
                   - f_680 * hk_78[k]
                   + f_681 * hk_87[k]
                   - f_682 * hk_100[k]
                   - f_683 * hk_253[k]
                   + f_684 * hk_258[k]
                   - f_685 * hk_267[k]
                   + f_686 * hk_280[k]
                   + f_679 * hk_577[k]
                   - f_680 * hk_582[k]
                   + f_681 * hk_591[k]
                   - f_682 * hk_604[k];
    }

#pragma omp simd aligned(hk_76, hk_83, hk_94, hk_256, hk_263, hk_274, hk_580, hk_587, \
                         hk_598 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_136[k] = f_687 * hk_76[k]
                   - f_688 * hk_83[k]
                   + f_687 * hk_94[k]
                   - f_689 * hk_256[k]
                   + f_690 * hk_263[k]
                   - f_689 * hk_274[k]
                   + f_687 * hk_580[k]
                   - f_688 * hk_587[k]
                   + f_687 * hk_598[k];
    }

#pragma omp simd aligned(hk_73, hk_78, hk_80, hk_87, hk_89, hk_100, hk_102, hk_253, hk_258, \
                         hk_260, hk_267, hk_269, hk_280, hk_282, hk_577, hk_582, hk_584, \
                         hk_591, hk_593, hk_604, hk_606 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_137[k] = -f_691 * hk_73[k]
                   + f_691 * hk_78[k]
                   + f_692 * hk_80[k]
                   + f_693 * hk_87[k]
                   - f_157 * hk_89[k]
                   - f_694 * hk_100[k]
                   + f_695 * hk_102[k]
                   + f_696 * hk_253[k]
                   - f_696 * hk_258[k]
                   - f_697 * hk_260[k]
                   - f_698 * hk_267[k]
                   + f_699 * hk_269[k]
                   + f_700 * hk_280[k]
                   - f_701 * hk_282[k]
                   - f_691 * hk_577[k]
                   + f_691 * hk_582[k]
                   + f_692 * hk_584[k]
                   + f_693 * hk_591[k]
                   - f_157 * hk_593[k]
                   - f_694 * hk_604[k]
                   + f_695 * hk_606[k];
    }

#pragma omp simd aligned(hk_76, hk_85, hk_94, hk_96, hk_256, hk_265, hk_274, hk_276, hk_580, \
                         hk_589, hk_598, hk_600 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_138[k] = -f_156 * hk_76[k]
                   + f_158 * hk_85[k]
                   + f_156 * hk_94[k]
                   - f_158 * hk_96[k]
                   + f_702 * hk_256[k]
                   - f_125 * hk_265[k]
                   - f_702 * hk_274[k]
                   + f_125 * hk_276[k]
                   - f_156 * hk_580[k]
                   + f_158 * hk_589[k]
                   + f_156 * hk_598[k]
                   - f_158 * hk_600[k];
    }

#pragma omp simd aligned(hk_73, hk_78, hk_80, hk_87, hk_89, hk_91, hk_100, hk_102, hk_104, \
                         hk_253, hk_258, hk_260, hk_267, hk_269, hk_271, hk_280, hk_282, \
                         hk_284, hk_577, hk_582, hk_584, hk_591, hk_593, hk_595, hk_604, \
                         hk_606, hk_608 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_139[k] = f_703 * hk_73[k]
                   + f_307 * hk_78[k]
                   - f_704 * hk_80[k]
                   + f_705 * hk_87[k]
                   - f_309 * hk_89[k]
                   + f_136 * hk_91[k]
                   - f_705 * hk_100[k]
                   + f_131 * hk_102[k]
                   - f_259 * hk_104[k]
                   - f_706 * hk_253[k]
                   - f_707 * hk_258[k]
                   + f_708 * hk_260[k]
                   - f_709 * hk_267[k]
                   + f_132 * hk_269[k]
                   - f_710 * hk_271[k]
                   + f_709 * hk_280[k]
                   - f_711 * hk_282[k]
                   + f_134 * hk_284[k]
                   + f_703 * hk_577[k]
                   + f_307 * hk_582[k]
                   - f_704 * hk_584[k]
                   + f_705 * hk_591[k]
                   - f_309 * hk_593[k]
                   + f_136 * hk_595[k]
                   - f_705 * hk_604[k]
                   + f_131 * hk_606[k]
                   - f_259 * hk_608[k];
    }

#pragma omp simd aligned(hk_76, hk_83, hk_85, hk_94, hk_96, hk_98, hk_256, hk_263, hk_265, \
                         hk_274, hk_276, hk_278, hk_580, hk_587, hk_589, hk_598, hk_600, \
                         hk_602 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_140[k] = f_227 * hk_76[k]
                   + f_153 * hk_83[k]
                   - f_245 * hk_85[k]
                   + f_227 * hk_94[k]
                   - f_245 * hk_96[k]
                   + f_712 * hk_98[k]
                   - f_713 * hk_256[k]
                   - f_714 * hk_263[k]
                   + f_239 * hk_265[k]
                   - f_713 * hk_274[k]
                   + f_239 * hk_276[k]
                   - f_715 * hk_278[k]
                   + f_227 * hk_580[k]
                   + f_153 * hk_587[k]
                   - f_245 * hk_589[k]
                   + f_227 * hk_598[k]
                   - f_245 * hk_600[k]
                   + f_712 * hk_602[k];
    }

#pragma omp simd aligned(hk_73, hk_78, hk_80, hk_87, hk_89, hk_91, hk_100, hk_102, hk_104, \
                         hk_106, hk_253, hk_258, hk_260, hk_267, hk_269, hk_271, hk_280, \
                         hk_282, hk_284, hk_286, hk_577, hk_582, hk_584, hk_591, hk_593, \
                         hk_595, hk_604, hk_606, hk_608, hk_610 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_141[k] = -f_716 * hk_73[k]
                   - f_717 * hk_78[k]
                   + f_365 * hk_80[k]
                   - f_717 * hk_87[k]
                   + f_362 * hk_89[k]
                   - f_362 * hk_91[k]
                   - f_716 * hk_100[k]
                   + f_365 * hk_102[k]
                   - f_362 * hk_104[k]
                   + f_718 * hk_106[k]
                   + f_665 * hk_253[k]
                   + f_719 * hk_258[k]
                   - f_360 * hk_260[k]
                   + f_719 * hk_267[k]
                   - f_366 * hk_269[k]
                   + f_366 * hk_271[k]
                   + f_665 * hk_280[k]
                   - f_360 * hk_282[k]
                   + f_366 * hk_284[k]
                   - f_720 * hk_286[k]
                   - f_716 * hk_577[k]
                   - f_717 * hk_582[k]
                   + f_365 * hk_584[k]
                   - f_717 * hk_591[k]
                   + f_362 * hk_593[k]
                   - f_362 * hk_595[k]
                   - f_716 * hk_604[k]
                   + f_365 * hk_606[k]
                   - f_362 * hk_608[k]
                   + f_718 * hk_610[k];
    }

#pragma omp simd aligned(hk_74, hk_79, hk_81, hk_88, hk_90, hk_92, hk_101, hk_103, hk_105, \
                         hk_107, hk_254, hk_259, hk_261, hk_268, hk_270, hk_272, hk_281, \
                         hk_283, hk_285, hk_287, hk_578, hk_583, hk_585, hk_592, hk_594, \
                         hk_596, hk_605, hk_607, hk_609, hk_611 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_142[k] = -f_721 * hk_74[k]
                   - f_722 * hk_79[k]
                   + f_723 * hk_81[k]
                   - f_722 * hk_88[k]
                   + f_148 * hk_90[k]
                   - f_724 * hk_92[k]
                   - f_721 * hk_101[k]
                   + f_723 * hk_103[k]
                   - f_724 * hk_105[k]
                   + f_725 * hk_107[k]
                   + f_723 * hk_254[k]
                   + f_726 * hk_259[k]
                   - f_727 * hk_261[k]
                   + f_726 * hk_268[k]
                   - f_728 * hk_270[k]
                   + f_729 * hk_272[k]
                   + f_723 * hk_281[k]
                   - f_727 * hk_283[k]
                   + f_729 * hk_285[k]
                   - f_730 * hk_287[k]
                   - f_721 * hk_578[k]
                   - f_722 * hk_583[k]
                   + f_723 * hk_585[k]
                   - f_722 * hk_592[k]
                   + f_148 * hk_594[k]
                   - f_724 * hk_596[k]
                   - f_721 * hk_605[k]
                   + f_723 * hk_607[k]
                   - f_724 * hk_609[k]
                   + f_725 * hk_611[k];
    }

#pragma omp simd aligned(hk_72, hk_75, hk_77, hk_82, hk_84, hk_86, hk_93, hk_95, hk_97, hk_99, \
                         hk_252, hk_255, hk_257, hk_262, hk_264, hk_266, hk_273, hk_275, \
                         hk_277, hk_279, hk_576, hk_579, hk_581, hk_586, hk_588, hk_590, \
                         hk_597, hk_599, hk_601, hk_603 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_143[k] = -f_716 * hk_72[k]
                   - f_717 * hk_75[k]
                   + f_365 * hk_77[k]
                   - f_717 * hk_82[k]
                   + f_362 * hk_84[k]
                   - f_362 * hk_86[k]
                   - f_716 * hk_93[k]
                   + f_365 * hk_95[k]
                   - f_362 * hk_97[k]
                   + f_718 * hk_99[k]
                   + f_665 * hk_252[k]
                   + f_719 * hk_255[k]
                   - f_360 * hk_257[k]
                   + f_719 * hk_262[k]
                   - f_366 * hk_264[k]
                   + f_366 * hk_266[k]
                   + f_665 * hk_273[k]
                   - f_360 * hk_275[k]
                   + f_366 * hk_277[k]
                   - f_720 * hk_279[k]
                   - f_716 * hk_576[k]
                   - f_717 * hk_579[k]
                   + f_365 * hk_581[k]
                   - f_717 * hk_586[k]
                   + f_362 * hk_588[k]
                   - f_362 * hk_590[k]
                   - f_716 * hk_597[k]
                   + f_365 * hk_599[k]
                   - f_362 * hk_601[k]
                   + f_718 * hk_603[k];
    }

#pragma omp simd aligned(hk_74, hk_79, hk_81, hk_88, hk_92, hk_101, hk_103, hk_105, hk_254, \
                         hk_259, hk_261, hk_268, hk_272, hk_281, hk_283, hk_285, hk_578, \
                         hk_583, hk_585, hk_592, hk_596, hk_605, hk_607, \
                         hk_609 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_144[k] = f_731 * hk_74[k]
                   + f_731 * hk_79[k]
                   - f_232 * hk_81[k]
                   - f_731 * hk_88[k]
                   + f_732 * hk_92[k]
                   - f_731 * hk_101[k]
                   + f_232 * hk_103[k]
                   - f_732 * hk_105[k]
                   - f_225 * hk_254[k]
                   - f_225 * hk_259[k]
                   + f_238 * hk_261[k]
                   + f_225 * hk_268[k]
                   - f_733 * hk_272[k]
                   + f_225 * hk_281[k]
                   - f_238 * hk_283[k]
                   + f_733 * hk_285[k]
                   + f_731 * hk_578[k]
                   + f_731 * hk_583[k]
                   - f_232 * hk_585[k]
                   - f_731 * hk_592[k]
                   + f_732 * hk_596[k]
                   - f_731 * hk_605[k]
                   + f_232 * hk_607[k]
                   - f_732 * hk_609[k];
    }

#pragma omp simd aligned(hk_72, hk_75, hk_77, hk_82, hk_84, hk_86, hk_93, hk_95, hk_97, \
                         hk_252, hk_255, hk_257, hk_262, hk_264, hk_266, hk_273, hk_275, \
                         hk_277, hk_576, hk_579, hk_581, hk_586, hk_588, hk_590, hk_597, \
                         hk_599, hk_601 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_145[k] = f_705 * hk_72[k]
                   - f_705 * hk_75[k]
                   - f_131 * hk_77[k]
                   - f_307 * hk_82[k]
                   + f_309 * hk_84[k]
                   + f_259 * hk_86[k]
                   - f_703 * hk_93[k]
                   + f_704 * hk_95[k]
                   - f_136 * hk_97[k]
                   - f_709 * hk_252[k]
                   + f_709 * hk_255[k]
                   + f_711 * hk_257[k]
                   + f_707 * hk_262[k]
                   - f_132 * hk_264[k]
                   - f_134 * hk_266[k]
                   + f_706 * hk_273[k]
                   - f_708 * hk_275[k]
                   + f_710 * hk_277[k]
                   + f_705 * hk_576[k]
                   - f_705 * hk_579[k]
                   - f_131 * hk_581[k]
                   - f_307 * hk_586[k]
                   + f_309 * hk_588[k]
                   + f_259 * hk_590[k]
                   - f_703 * hk_597[k]
                   + f_704 * hk_599[k]
                   - f_136 * hk_601[k];
    }

#pragma omp simd aligned(hk_74, hk_79, hk_81, hk_88, hk_90, hk_101, hk_103, hk_254, hk_259, \
                         hk_261, hk_268, hk_270, hk_281, hk_283, hk_578, hk_583, hk_585, \
                         hk_592, hk_594, hk_605, hk_607 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_146[k] = -f_700 * hk_74[k]
                   + f_696 * hk_79[k]
                   + f_122 * hk_81[k]
                   + f_696 * hk_88[k]
                   - f_157 * hk_90[k]
                   - f_700 * hk_101[k]
                   + f_122 * hk_103[k]
                   + f_124 * hk_254[k]
                   - f_734 * hk_259[k]
                   - f_157 * hk_261[k]
                   - f_734 * hk_268[k]
                   + f_699 * hk_270[k]
                   + f_124 * hk_281[k]
                   - f_157 * hk_283[k]
                   - f_700 * hk_578[k]
                   + f_696 * hk_583[k]
                   + f_122 * hk_585[k]
                   + f_696 * hk_592[k]
                   - f_157 * hk_594[k]
                   - f_700 * hk_605[k]
                   + f_122 * hk_607[k];
    }

#pragma omp simd aligned(hk_72, hk_75, hk_77, hk_82, hk_84, hk_93, hk_95, hk_252, hk_255, \
                         hk_257, hk_262, hk_264, hk_273, hk_275, hk_576, hk_579, hk_581, \
                         hk_586, hk_588, hk_597, hk_599 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_147[k] = -f_694 * hk_72[k]
                   + f_693 * hk_75[k]
                   + f_695 * hk_77[k]
                   + f_691 * hk_82[k]
                   - f_157 * hk_84[k]
                   - f_691 * hk_93[k]
                   + f_692 * hk_95[k]
                   + f_700 * hk_252[k]
                   - f_698 * hk_255[k]
                   - f_701 * hk_257[k]
                   - f_696 * hk_262[k]
                   + f_699 * hk_264[k]
                   + f_696 * hk_273[k]
                   - f_697 * hk_275[k]
                   - f_694 * hk_576[k]
                   + f_693 * hk_579[k]
                   + f_695 * hk_581[k]
                   + f_691 * hk_586[k]
                   - f_157 * hk_588[k]
                   - f_691 * hk_597[k]
                   + f_692 * hk_599[k];
    }

#pragma omp simd aligned(hk_74, hk_79, hk_88, hk_101, hk_254, hk_259, hk_268, hk_281, hk_578, \
                         hk_583, hk_592, hk_605 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_148[k] = f_735 * hk_74[k]
                   - f_736 * hk_79[k]
                   + f_736 * hk_88[k]
                   - f_735 * hk_101[k]
                   - f_687 * hk_254[k]
                   + f_737 * hk_259[k]
                   - f_737 * hk_268[k]
                   + f_687 * hk_281[k]
                   + f_735 * hk_578[k]
                   - f_736 * hk_583[k]
                   + f_736 * hk_592[k]
                   - f_735 * hk_605[k];
    }

#pragma omp simd aligned(hk_72, hk_75, hk_82, hk_93, hk_252, hk_255, hk_262, hk_273, hk_576, \
                         hk_579, hk_586, hk_597 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_149[k] = f_682 * hk_72[k]
                   - f_681 * hk_75[k]
                   + f_680 * hk_82[k]
                   - f_679 * hk_93[k]
                   - f_686 * hk_252[k]
                   + f_685 * hk_255[k]
                   - f_684 * hk_262[k]
                   + f_683 * hk_273[k]
                   + f_682 * hk_576[k]
                   - f_681 * hk_579[k]
                   + f_680 * hk_586[k]
                   - f_679 * hk_597[k];
    }

#pragma omp simd aligned(hk_1, hk_6, hk_15, hk_28, hk_109, hk_114, hk_123, hk_136, hk_361, \
                         hk_366, hk_375, hk_388 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_150[k] = f_8 * hk_1[k]
                   - f_0 * hk_6[k]
                   + f_9 * hk_15[k]
                   - f_10 * hk_28[k]
                   - f_4 * hk_109[k]
                   + f_5 * hk_114[k]
                   - f_6 * hk_123[k]
                   + f_7 * hk_136[k]
                   + f_0 * hk_361[k]
                   - f_1 * hk_366[k]
                   + f_2 * hk_375[k]
                   - f_3 * hk_388[k];
    }

#pragma omp simd aligned(hk_4, hk_11, hk_22, hk_112, hk_119, hk_130, hk_364, hk_371, \
                         hk_382 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_151[k] = f_15 * hk_4[k]
                   - f_16 * hk_11[k]
                   + f_15 * hk_22[k]
                   - f_13 * hk_112[k]
                   + f_14 * hk_119[k]
                   - f_13 * hk_130[k]
                   + f_11 * hk_364[k]
                   - f_12 * hk_371[k]
                   + f_11 * hk_382[k];
    }

#pragma omp simd aligned(hk_1, hk_6, hk_8, hk_15, hk_17, hk_28, hk_30, hk_109, hk_114, hk_116, \
                         hk_123, hk_125, hk_136, hk_138, hk_361, hk_366, hk_368, hk_375, \
                         hk_377, hk_388, hk_390 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_152[k] = -f_21 * hk_1[k]
                   + f_21 * hk_6[k]
                   + f_22 * hk_8[k]
                   + f_28 * hk_15[k]
                   - f_27 * hk_17[k]
                   - f_29 * hk_28[k]
                   + f_30 * hk_30[k]
                   + f_23 * hk_109[k]
                   - f_23 * hk_114[k]
                   - f_20 * hk_116[k]
                   - f_24 * hk_123[k]
                   + f_25 * hk_125[k]
                   + f_26 * hk_136[k]
                   - f_27 * hk_138[k]
                   - f_17 * hk_361[k]
                   + f_17 * hk_366[k]
                   + f_18 * hk_368[k]
                   + f_19 * hk_375[k]
                   - f_20 * hk_377[k]
                   - f_21 * hk_388[k]
                   + f_22 * hk_390[k];
    }

#pragma omp simd aligned(hk_4, hk_13, hk_22, hk_24, hk_112, hk_121, hk_130, hk_132, hk_364, \
                         hk_373, hk_382, hk_384 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_153[k] = -f_34 * hk_4[k]
                   + f_35 * hk_13[k]
                   + f_34 * hk_22[k]
                   - f_35 * hk_24[k]
                   + f_32 * hk_112[k]
                   - f_33 * hk_121[k]
                   - f_32 * hk_130[k]
                   + f_33 * hk_132[k]
                   - f_27 * hk_364[k]
                   + f_31 * hk_373[k]
                   + f_27 * hk_382[k]
                   - f_31 * hk_384[k];
    }

#pragma omp simd aligned(hk_1, hk_6, hk_8, hk_15, hk_17, hk_19, hk_28, hk_30, hk_32, hk_109, \
                         hk_114, hk_116, hk_123, hk_125, hk_127, hk_136, hk_138, hk_140, \
                         hk_361, hk_366, hk_368, hk_375, hk_377, hk_379, hk_388, hk_390, \
                         hk_392 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_154[k] = f_50 * hk_1[k]
                   + f_39 * hk_6[k]
                   - f_51 * hk_8[k]
                   + f_52 * hk_15[k]
                   - f_53 * hk_17[k]
                   + f_54 * hk_19[k]
                   - f_52 * hk_28[k]
                   + f_55 * hk_30[k]
                   - f_56 * hk_32[k]
                   - f_44 * hk_109[k]
                   - f_45 * hk_114[k]
                   + f_46 * hk_116[k]
                   - f_47 * hk_123[k]
                   + f_41 * hk_125[k]
                   - f_48 * hk_127[k]
                   + f_47 * hk_136[k]
                   - f_40 * hk_138[k]
                   + f_49 * hk_140[k]
                   + f_36 * hk_361[k]
                   + f_37 * hk_366[k]
                   - f_38 * hk_368[k]
                   + f_39 * hk_375[k]
                   - f_40 * hk_377[k]
                   + f_41 * hk_379[k]
                   - f_39 * hk_388[k]
                   + f_42 * hk_390[k]
                   - f_43 * hk_392[k];
    }

#pragma omp simd aligned(hk_4, hk_11, hk_13, hk_22, hk_24, hk_26, hk_112, hk_119, hk_121, \
                         hk_130, hk_132, hk_134, hk_364, hk_371, hk_373, hk_382, hk_384, \
                         hk_386 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_155[k] = f_64 * hk_4[k]
                   + f_65 * hk_11[k]
                   - f_66 * hk_13[k]
                   + f_64 * hk_22[k]
                   - f_66 * hk_24[k]
                   + f_67 * hk_26[k]
                   - f_58 * hk_112[k]
                   - f_61 * hk_119[k]
                   + f_62 * hk_121[k]
                   - f_58 * hk_130[k]
                   + f_62 * hk_132[k]
                   - f_63 * hk_134[k]
                   + f_57 * hk_364[k]
                   + f_58 * hk_371[k]
                   - f_59 * hk_373[k]
                   + f_57 * hk_382[k]
                   - f_59 * hk_384[k]
                   + f_60 * hk_386[k];
    }

#pragma omp simd aligned(hk_1, hk_6, hk_8, hk_15, hk_17, hk_19, hk_28, hk_30, hk_32, hk_34, \
                         hk_109, hk_114, hk_116, hk_123, hk_125, hk_127, hk_136, hk_138, \
                         hk_140, hk_142, hk_361, hk_366, hk_368, hk_375, hk_377, hk_379, \
                         hk_388, hk_390, hk_392, hk_394 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_156[k] = -f_77 * hk_1[k]
                   - f_78 * hk_6[k]
                   + f_79 * hk_8[k]
                   - f_78 * hk_15[k]
                   + f_80 * hk_17[k]
                   - f_80 * hk_19[k]
                   - f_77 * hk_28[k]
                   + f_79 * hk_30[k]
                   - f_80 * hk_32[k]
                   + f_81 * hk_34[k]
                   + f_73 * hk_109[k]
                   + f_74 * hk_114[k]
                   - f_71 * hk_116[k]
                   + f_74 * hk_123[k]
                   - f_75 * hk_125[k]
                   + f_75 * hk_127[k]
                   + f_73 * hk_136[k]
                   - f_71 * hk_138[k]
                   + f_75 * hk_140[k]
                   - f_76 * hk_142[k]
                   - f_68 * hk_361[k]
                   - f_69 * hk_366[k]
                   + f_70 * hk_368[k]
                   - f_69 * hk_375[k]
                   + f_71 * hk_377[k]
                   - f_71 * hk_379[k]
                   - f_68 * hk_388[k]
                   + f_70 * hk_390[k]
                   - f_71 * hk_392[k]
                   + f_72 * hk_394[k];
    }

#pragma omp simd aligned(hk_2, hk_7, hk_9, hk_16, hk_18, hk_20, hk_29, hk_31, hk_33, hk_35, \
                         hk_110, hk_115, hk_117, hk_124, hk_126, hk_128, hk_137, hk_139, \
                         hk_141, hk_143, hk_362, hk_367, hk_369, hk_376, hk_378, hk_380, \
                         hk_389, hk_391, hk_393, hk_395 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_157[k] = -f_92 * hk_2[k]
                   - f_93 * hk_7[k]
                   + f_94 * hk_9[k]
                   - f_93 * hk_16[k]
                   + f_95 * hk_18[k]
                   - f_96 * hk_20[k]
                   - f_92 * hk_29[k]
                   + f_94 * hk_31[k]
                   - f_96 * hk_33[k]
                   + f_97 * hk_35[k]
                   + f_88 * hk_110[k]
                   + f_84 * hk_115[k]
                   - f_85 * hk_117[k]
                   + f_84 * hk_124[k]
                   - f_89 * hk_126[k]
                   + f_90 * hk_128[k]
                   + f_88 * hk_137[k]
                   - f_85 * hk_139[k]
                   + f_90 * hk_141[k]
                   - f_91 * hk_143[k]
                   - f_82 * hk_362[k]
                   - f_83 * hk_367[k]
                   + f_84 * hk_369[k]
                   - f_83 * hk_376[k]
                   + f_85 * hk_378[k]
                   - f_86 * hk_380[k]
                   - f_82 * hk_389[k]
                   + f_84 * hk_391[k]
                   - f_86 * hk_393[k]
                   + f_87 * hk_395[k];
    }

#pragma omp simd aligned(hk_0, hk_3, hk_5, hk_10, hk_12, hk_14, hk_21, hk_23, hk_25, hk_27, \
                         hk_108, hk_111, hk_113, hk_118, hk_120, hk_122, hk_129, hk_131, \
                         hk_133, hk_135, hk_360, hk_363, hk_365, hk_370, hk_372, hk_374, \
                         hk_381, hk_383, hk_385, hk_387 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_158[k] = -f_77 * hk_0[k]
                   - f_78 * hk_3[k]
                   + f_79 * hk_5[k]
                   - f_78 * hk_10[k]
                   + f_80 * hk_12[k]
                   - f_80 * hk_14[k]
                   - f_77 * hk_21[k]
                   + f_79 * hk_23[k]
                   - f_80 * hk_25[k]
                   + f_81 * hk_27[k]
                   + f_73 * hk_108[k]
                   + f_74 * hk_111[k]
                   - f_71 * hk_113[k]
                   + f_74 * hk_118[k]
                   - f_75 * hk_120[k]
                   + f_75 * hk_122[k]
                   + f_73 * hk_129[k]
                   - f_71 * hk_131[k]
                   + f_75 * hk_133[k]
                   - f_76 * hk_135[k]
                   - f_68 * hk_360[k]
                   - f_69 * hk_363[k]
                   + f_70 * hk_365[k]
                   - f_69 * hk_370[k]
                   + f_71 * hk_372[k]
                   - f_71 * hk_374[k]
                   - f_68 * hk_381[k]
                   + f_70 * hk_383[k]
                   - f_71 * hk_385[k]
                   + f_72 * hk_387[k];
    }

#pragma omp simd aligned(hk_2, hk_7, hk_9, hk_16, hk_20, hk_29, hk_31, hk_33, hk_110, hk_115, \
                         hk_117, hk_124, hk_128, hk_137, hk_139, hk_141, hk_362, hk_367, \
                         hk_369, hk_376, hk_380, hk_389, hk_391, \
                         hk_393 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_159[k] = f_101 * hk_2[k]
                   + f_101 * hk_7[k]
                   - f_102 * hk_9[k]
                   - f_101 * hk_16[k]
                   + f_103 * hk_20[k]
                   - f_101 * hk_29[k]
                   + f_102 * hk_31[k]
                   - f_103 * hk_33[k]
                   - f_57 * hk_110[k]
                   - f_57 * hk_115[k]
                   + f_59 * hk_117[k]
                   + f_57 * hk_124[k]
                   - f_60 * hk_128[k]
                   + f_57 * hk_137[k]
                   - f_59 * hk_139[k]
                   + f_60 * hk_141[k]
                   + f_98 * hk_362[k]
                   + f_98 * hk_367[k]
                   - f_99 * hk_369[k]
                   - f_98 * hk_376[k]
                   + f_100 * hk_380[k]
                   - f_98 * hk_389[k]
                   + f_99 * hk_391[k]
                   - f_100 * hk_393[k];
    }

#pragma omp simd aligned(hk_0, hk_3, hk_5, hk_10, hk_12, hk_14, hk_21, hk_23, hk_25, hk_108, \
                         hk_111, hk_113, hk_118, hk_120, hk_122, hk_129, hk_131, hk_133, \
                         hk_360, hk_363, hk_365, hk_370, hk_372, hk_374, hk_381, hk_383, \
                         hk_385 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_160[k] = f_52 * hk_0[k]
                   - f_52 * hk_3[k]
                   - f_55 * hk_5[k]
                   - f_39 * hk_10[k]
                   + f_53 * hk_12[k]
                   + f_56 * hk_14[k]
                   - f_50 * hk_21[k]
                   + f_51 * hk_23[k]
                   - f_54 * hk_25[k]
                   - f_47 * hk_108[k]
                   + f_47 * hk_111[k]
                   + f_40 * hk_113[k]
                   + f_45 * hk_118[k]
                   - f_41 * hk_120[k]
                   - f_49 * hk_122[k]
                   + f_44 * hk_129[k]
                   - f_46 * hk_131[k]
                   + f_48 * hk_133[k]
                   + f_39 * hk_360[k]
                   - f_39 * hk_363[k]
                   - f_42 * hk_365[k]
                   - f_37 * hk_370[k]
                   + f_40 * hk_372[k]
                   + f_43 * hk_374[k]
                   - f_36 * hk_381[k]
                   + f_38 * hk_383[k]
                   - f_41 * hk_385[k];
    }

#pragma omp simd aligned(hk_2, hk_7, hk_9, hk_16, hk_18, hk_29, hk_31, hk_110, hk_115, hk_117, \
                         hk_124, hk_126, hk_137, hk_139, hk_362, hk_367, hk_369, hk_376, \
                         hk_378, hk_389, hk_391 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_161[k] = -f_108 * hk_2[k]
                   + f_104 * hk_7[k]
                   + f_109 * hk_9[k]
                   + f_104 * hk_16[k]
                   - f_27 * hk_18[k]
                   - f_108 * hk_29[k]
                   + f_109 * hk_31[k]
                   + f_22 * hk_110[k]
                   - f_18 * hk_115[k]
                   - f_107 * hk_117[k]
                   - f_18 * hk_124[k]
                   + f_25 * hk_126[k]
                   + f_22 * hk_137[k]
                   - f_107 * hk_139[k]
                   - f_104 * hk_362[k]
                   + f_105 * hk_367[k]
                   + f_106 * hk_369[k]
                   + f_105 * hk_376[k]
                   - f_20 * hk_378[k]
                   - f_104 * hk_389[k]
                   + f_106 * hk_391[k];
    }

#pragma omp simd aligned(hk_0, hk_3, hk_5, hk_10, hk_12, hk_21, hk_23, hk_108, hk_111, hk_113, \
                         hk_118, hk_120, hk_129, hk_131, hk_360, hk_363, hk_365, hk_370, \
                         hk_372, hk_381, hk_383 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_162[k] = -f_29 * hk_0[k]
                   + f_28 * hk_3[k]
                   + f_30 * hk_5[k]
                   + f_21 * hk_10[k]
                   - f_27 * hk_12[k]
                   - f_21 * hk_21[k]
                   + f_22 * hk_23[k]
                   + f_26 * hk_108[k]
                   - f_24 * hk_111[k]
                   - f_27 * hk_113[k]
                   - f_23 * hk_118[k]
                   + f_25 * hk_120[k]
                   + f_23 * hk_129[k]
                   - f_20 * hk_131[k]
                   - f_21 * hk_360[k]
                   + f_19 * hk_363[k]
                   + f_22 * hk_365[k]
                   + f_17 * hk_370[k]
                   - f_20 * hk_372[k]
                   - f_17 * hk_381[k]
                   + f_18 * hk_383[k];
    }

#pragma omp simd aligned(hk_2, hk_7, hk_16, hk_29, hk_110, hk_115, hk_124, hk_137, hk_362, \
                         hk_367, hk_376, hk_389 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_163[k] = f_114 * hk_2[k]
                   - f_115 * hk_7[k]
                   + f_115 * hk_16[k]
                   - f_114 * hk_29[k]
                   - f_112 * hk_110[k]
                   + f_113 * hk_115[k]
                   - f_113 * hk_124[k]
                   + f_112 * hk_137[k]
                   + f_110 * hk_362[k]
                   - f_111 * hk_367[k]
                   + f_111 * hk_376[k]
                   - f_110 * hk_389[k];
    }

#pragma omp simd aligned(hk_0, hk_3, hk_10, hk_21, hk_108, hk_111, hk_118, hk_129, hk_360, \
                         hk_363, hk_370, hk_381 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_164[k] = f_10 * hk_0[k]
                   - f_9 * hk_3[k]
                   + f_0 * hk_10[k]
                   - f_8 * hk_21[k]
                   - f_7 * hk_108[k]
                   + f_6 * hk_111[k]
                   - f_5 * hk_118[k]
                   + f_4 * hk_129[k]
                   + f_3 * hk_360[k]
                   - f_2 * hk_363[k]
                   + f_1 * hk_370[k]
                   - f_0 * hk_381[k];
    }
}

}  // namespace simdtrf
