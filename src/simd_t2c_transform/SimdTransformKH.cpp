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


#include "SimdTransformKH.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_kh(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t kh,
             const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 0.205078125 * std::sqrt(6006.0);
    const auto f_1 = 0.41015625 * std::sqrt(6006.0);
    const auto f_2 = 0.041015625 * std::sqrt(6006.0);
    const auto f_3 = 1.025390625 * std::sqrt(6006.0);
    const auto f_4 = 2.05078125 * std::sqrt(6006.0);
    const auto f_5 = 0.615234375 * std::sqrt(6006.0);
    const auto f_6 = 1.23046875 * std::sqrt(6006.0);
    const auto f_7 = 0.123046875 * std::sqrt(6006.0);
    const auto f_8 = 0.029296875 * std::sqrt(6006.0);
    const auto f_9 = 0.05859375 * std::sqrt(6006.0);
    const auto f_10 = 0.005859375 * std::sqrt(6006.0);
    const auto f_11 = 0.328125 * std::sqrt(15015.0);
    const auto f_12 = 1.640625 * std::sqrt(15015.0);
    const auto f_13 = 0.984375 * std::sqrt(15015.0);
    const auto f_14 = 0.046875 * std::sqrt(15015.0);
    const auto f_15 = 0.041015625 * std::sqrt(30030.0);
    const auto f_16 = 0.02734375 * std::sqrt(30030.0);
    const auto f_17 = 0.328125 * std::sqrt(30030.0);
    const auto f_18 = 0.013671875 * std::sqrt(30030.0);
    const auto f_19 = 0.109375 * std::sqrt(30030.0);
    const auto f_20 = 0.205078125 * std::sqrt(30030.0);
    const auto f_21 = 0.13671875 * std::sqrt(30030.0);
    const auto f_22 = 1.640625 * std::sqrt(30030.0);
    const auto f_23 = 0.068359375 * std::sqrt(30030.0);
    const auto f_24 = 0.546875 * std::sqrt(30030.0);
    const auto f_25 = 0.123046875 * std::sqrt(30030.0);
    const auto f_26 = 0.08203125 * std::sqrt(30030.0);
    const auto f_27 = 0.984375 * std::sqrt(30030.0);
    const auto f_28 = 0.005859375 * std::sqrt(30030.0);
    const auto f_29 = 0.00390625 * std::sqrt(30030.0);
    const auto f_30 = 0.046875 * std::sqrt(30030.0);
    const auto f_31 = 0.001953125 * std::sqrt(30030.0);
    const auto f_32 = 0.015625 * std::sqrt(30030.0);
    const auto f_33 = 0.328125 * std::sqrt(5005.0);
    const auto f_34 = 0.65625 * std::sqrt(5005.0);
    const auto f_35 = 1.640625 * std::sqrt(5005.0);
    const auto f_36 = 3.28125 * std::sqrt(5005.0);
    const auto f_37 = 0.984375 * std::sqrt(5005.0);
    const auto f_38 = 1.96875 * std::sqrt(5005.0);
    const auto f_39 = 0.046875 * std::sqrt(5005.0);
    const auto f_40 = 0.09375 * std::sqrt(5005.0);
    const auto f_41 = 0.08203125 * std::sqrt(715.0);
    const auto f_42 = 0.1640625 * std::sqrt(715.0);
    const auto f_43 = 0.984375 * std::sqrt(715.0);
    const auto f_44 = 0.65625 * std::sqrt(715.0);
    const auto f_45 = 0.41015625 * std::sqrt(715.0);
    const auto f_46 = 0.8203125 * std::sqrt(715.0);
    const auto f_47 = 4.921875 * std::sqrt(715.0);
    const auto f_48 = 3.28125 * std::sqrt(715.0);
    const auto f_49 = 0.24609375 * std::sqrt(715.0);
    const auto f_50 = 0.4921875 * std::sqrt(715.0);
    const auto f_51 = 2.953125 * std::sqrt(715.0);
    const auto f_52 = 1.96875 * std::sqrt(715.0);
    const auto f_53 = 0.01171875 * std::sqrt(715.0);
    const auto f_54 = 0.0234375 * std::sqrt(715.0);
    const auto f_55 = 0.140625 * std::sqrt(715.0);
    const auto f_56 = 0.09375 * std::sqrt(715.0);
    const auto f_57 = 0.41015625 * std::sqrt(429.0);
    const auto f_58 = 0.8203125 * std::sqrt(429.0);
    const auto f_59 = 1.09375 * std::sqrt(429.0);
    const auto f_60 = 0.21875 * std::sqrt(429.0);
    const auto f_61 = 2.05078125 * std::sqrt(429.0);
    const auto f_62 = 4.1015625 * std::sqrt(429.0);
    const auto f_63 = 5.46875 * std::sqrt(429.0);
    const auto f_64 = 1.23046875 * std::sqrt(429.0);
    const auto f_65 = 2.4609375 * std::sqrt(429.0);
    const auto f_66 = 3.28125 * std::sqrt(429.0);
    const auto f_67 = 0.65625 * std::sqrt(429.0);
    const auto f_68 = 0.05859375 * std::sqrt(429.0);
    const auto f_69 = 0.1171875 * std::sqrt(429.0);
    const auto f_70 = 0.15625 * std::sqrt(429.0);
    const auto f_71 = 0.03125 * std::sqrt(429.0);
    const auto f_72 = 0.1640625 * std::sqrt(5005.0);
    const auto f_73 = 0.8203125 * std::sqrt(5005.0);
    const auto f_74 = 0.4921875 * std::sqrt(5005.0);
    const auto f_75 = 0.0234375 * std::sqrt(5005.0);
    const auto f_76 = 0.08203125 * std::sqrt(15015.0);
    const auto f_77 = 0.4921875 * std::sqrt(15015.0);
    const auto f_78 = 0.41015625 * std::sqrt(15015.0);
    const auto f_79 = 2.4609375 * std::sqrt(15015.0);
    const auto f_80 = 0.24609375 * std::sqrt(15015.0);
    const auto f_81 = 1.4765625 * std::sqrt(15015.0);
    const auto f_82 = 0.01171875 * std::sqrt(15015.0);
    const auto f_83 = 0.0703125 * std::sqrt(15015.0);
    const auto f_84 = 4.921875 * std::sqrt(429.0);
    const auto f_85 = 0.4921875 * std::sqrt(429.0);
    const auto f_86 = 8.203125 * std::sqrt(429.0);
    const auto f_87 = 16.40625 * std::sqrt(429.0);
    const auto f_88 = 1.640625 * std::sqrt(429.0);
    const auto f_89 = 1.96875 * std::sqrt(4290.0);
    const auto f_90 = 6.5625 * std::sqrt(4290.0);
    const auto f_91 = 0.4921875 * std::sqrt(2145.0);
    const auto f_92 = 0.328125 * std::sqrt(2145.0);
    const auto f_93 = 3.9375 * std::sqrt(2145.0);
    const auto f_94 = 0.1640625 * std::sqrt(2145.0);
    const auto f_95 = 1.3125 * std::sqrt(2145.0);
    const auto f_96 = 1.640625 * std::sqrt(2145.0);
    const auto f_97 = 1.09375 * std::sqrt(2145.0);
    const auto f_98 = 13.125 * std::sqrt(2145.0);
    const auto f_99 = 0.546875 * std::sqrt(2145.0);
    const auto f_100 = 4.375 * std::sqrt(2145.0);
    const auto f_101 = 1.96875 * std::sqrt(1430.0);
    const auto f_102 = 3.9375 * std::sqrt(1430.0);
    const auto f_103 = 6.5625 * std::sqrt(1430.0);
    const auto f_104 = 13.125 * std::sqrt(1430.0);
    const auto f_105 = 0.0703125 * std::sqrt(10010.0);
    const auto f_106 = 0.140625 * std::sqrt(10010.0);
    const auto f_107 = 0.84375 * std::sqrt(10010.0);
    const auto f_108 = 0.5625 * std::sqrt(10010.0);
    const auto f_109 = 0.234375 * std::sqrt(10010.0);
    const auto f_110 = 0.46875 * std::sqrt(10010.0);
    const auto f_111 = 2.8125 * std::sqrt(10010.0);
    const auto f_112 = 1.875 * std::sqrt(10010.0);
    const auto f_113 = 0.3515625 * std::sqrt(6006.0);
    const auto f_114 = 0.703125 * std::sqrt(6006.0);
    const auto f_115 = 0.9375 * std::sqrt(6006.0);
    const auto f_116 = 0.1875 * std::sqrt(6006.0);
    const auto f_117 = 1.171875 * std::sqrt(6006.0);
    const auto f_118 = 2.34375 * std::sqrt(6006.0);
    const auto f_119 = 3.125 * std::sqrt(6006.0);
    const auto f_120 = 0.625 * std::sqrt(6006.0);
    const auto f_121 = 0.984375 * std::sqrt(1430.0);
    const auto f_122 = 3.28125 * std::sqrt(1430.0);
    const auto f_123 = 0.4921875 * std::sqrt(4290.0);
    const auto f_124 = 2.953125 * std::sqrt(4290.0);
    const auto f_125 = 1.640625 * std::sqrt(4290.0);
    const auto f_126 = 9.84375 * std::sqrt(4290.0);
    const auto f_127 = 1.025390625 * std::sqrt(66.0);
    const auto f_128 = 2.05078125 * std::sqrt(66.0);
    const auto f_129 = 0.205078125 * std::sqrt(66.0);
    const auto f_130 = 12.3046875 * std::sqrt(66.0);
    const auto f_131 = 24.609375 * std::sqrt(66.0);
    const auto f_132 = 2.4609375 * std::sqrt(66.0);
    const auto f_133 = 1.845703125 * std::sqrt(66.0);
    const auto f_134 = 3.69140625 * std::sqrt(66.0);
    const auto f_135 = 0.369140625 * std::sqrt(66.0);
    const auto f_136 = 49.21875 * std::sqrt(66.0);
    const auto f_137 = 4.921875 * std::sqrt(66.0);
    const auto f_138 = 0.41015625 * std::sqrt(66.0);
    const auto f_139 = 0.041015625 * std::sqrt(66.0);
    const auto f_140 = 0.4921875 * std::sqrt(66.0);
    const auto f_141 = 1.640625 * std::sqrt(165.0);
    const auto f_142 = 19.6875 * std::sqrt(165.0);
    const auto f_143 = 2.953125 * std::sqrt(165.0);
    const auto f_144 = 39.375 * std::sqrt(165.0);
    const auto f_145 = 0.328125 * std::sqrt(165.0);
    const auto f_146 = 3.9375 * std::sqrt(165.0);
    const auto f_147 = 0.205078125 * std::sqrt(330.0);
    const auto f_148 = 0.13671875 * std::sqrt(330.0);
    const auto f_149 = 1.640625 * std::sqrt(330.0);
    const auto f_150 = 0.068359375 * std::sqrt(330.0);
    const auto f_151 = 0.546875 * std::sqrt(330.0);
    const auto f_152 = 2.4609375 * std::sqrt(330.0);
    const auto f_153 = 19.6875 * std::sqrt(330.0);
    const auto f_154 = 0.8203125 * std::sqrt(330.0);
    const auto f_155 = 6.5625 * std::sqrt(330.0);
    const auto f_156 = 0.369140625 * std::sqrt(330.0);
    const auto f_157 = 0.24609375 * std::sqrt(330.0);
    const auto f_158 = 2.953125 * std::sqrt(330.0);
    const auto f_159 = 0.123046875 * std::sqrt(330.0);
    const auto f_160 = 0.984375 * std::sqrt(330.0);
    const auto f_161 = 4.921875 * std::sqrt(330.0);
    const auto f_162 = 3.28125 * std::sqrt(330.0);
    const auto f_163 = 39.375 * std::sqrt(330.0);
    const auto f_164 = 13.125 * std::sqrt(330.0);
    const auto f_165 = 0.041015625 * std::sqrt(330.0);
    const auto f_166 = 0.02734375 * std::sqrt(330.0);
    const auto f_167 = 0.328125 * std::sqrt(330.0);
    const auto f_168 = 0.013671875 * std::sqrt(330.0);
    const auto f_169 = 0.109375 * std::sqrt(330.0);
    const auto f_170 = 0.4921875 * std::sqrt(330.0);
    const auto f_171 = 3.9375 * std::sqrt(330.0);
    const auto f_172 = 0.1640625 * std::sqrt(330.0);
    const auto f_173 = 1.3125 * std::sqrt(330.0);
    const auto f_174 = 1.640625 * std::sqrt(55.0);
    const auto f_175 = 3.28125 * std::sqrt(55.0);
    const auto f_176 = 19.6875 * std::sqrt(55.0);
    const auto f_177 = 39.375 * std::sqrt(55.0);
    const auto f_178 = 2.953125 * std::sqrt(55.0);
    const auto f_179 = 5.90625 * std::sqrt(55.0);
    const auto f_180 = 78.75 * std::sqrt(55.0);
    const auto f_181 = 0.328125 * std::sqrt(55.0);
    const auto f_182 = 0.65625 * std::sqrt(55.0);
    const auto f_183 = 3.9375 * std::sqrt(55.0);
    const auto f_184 = 7.875 * std::sqrt(55.0);
    const auto f_185 = 0.05859375 * std::sqrt(385.0);
    const auto f_186 = 0.1171875 * std::sqrt(385.0);
    const auto f_187 = 0.703125 * std::sqrt(385.0);
    const auto f_188 = 0.46875 * std::sqrt(385.0);
    const auto f_189 = 1.40625 * std::sqrt(385.0);
    const auto f_190 = 8.4375 * std::sqrt(385.0);
    const auto f_191 = 5.625 * std::sqrt(385.0);
    const auto f_192 = 0.10546875 * std::sqrt(385.0);
    const auto f_193 = 0.2109375 * std::sqrt(385.0);
    const auto f_194 = 1.265625 * std::sqrt(385.0);
    const auto f_195 = 0.84375 * std::sqrt(385.0);
    const auto f_196 = 2.8125 * std::sqrt(385.0);
    const auto f_197 = 16.875 * std::sqrt(385.0);
    const auto f_198 = 11.25 * std::sqrt(385.0);
    const auto f_199 = 0.01171875 * std::sqrt(385.0);
    const auto f_200 = 0.0234375 * std::sqrt(385.0);
    const auto f_201 = 0.140625 * std::sqrt(385.0);
    const auto f_202 = 0.09375 * std::sqrt(385.0);
    const auto f_203 = 0.28125 * std::sqrt(385.0);
    const auto f_204 = 1.6875 * std::sqrt(385.0);
    const auto f_205 = 1.125 * std::sqrt(385.0);
    const auto f_206 = 0.29296875 * std::sqrt(231.0);
    const auto f_207 = 0.5859375 * std::sqrt(231.0);
    const auto f_208 = 0.78125 * std::sqrt(231.0);
    const auto f_209 = 0.15625 * std::sqrt(231.0);
    const auto f_210 = 3.515625 * std::sqrt(231.0);
    const auto f_211 = 7.03125 * std::sqrt(231.0);
    const auto f_212 = 9.375 * std::sqrt(231.0);
    const auto f_213 = 1.875 * std::sqrt(231.0);
    const auto f_214 = 0.52734375 * std::sqrt(231.0);
    const auto f_215 = 1.0546875 * std::sqrt(231.0);
    const auto f_216 = 1.40625 * std::sqrt(231.0);
    const auto f_217 = 0.28125 * std::sqrt(231.0);
    const auto f_218 = 14.0625 * std::sqrt(231.0);
    const auto f_219 = 18.75 * std::sqrt(231.0);
    const auto f_220 = 3.75 * std::sqrt(231.0);
    const auto f_221 = 0.05859375 * std::sqrt(231.0);
    const auto f_222 = 0.1171875 * std::sqrt(231.0);
    const auto f_223 = 0.03125 * std::sqrt(231.0);
    const auto f_224 = 0.703125 * std::sqrt(231.0);
    const auto f_225 = 0.375 * std::sqrt(231.0);
    const auto f_226 = 0.8203125 * std::sqrt(55.0);
    const auto f_227 = 9.84375 * std::sqrt(55.0);
    const auto f_228 = 1.4765625 * std::sqrt(55.0);
    const auto f_229 = 0.1640625 * std::sqrt(55.0);
    const auto f_230 = 1.96875 * std::sqrt(55.0);
    const auto f_231 = 0.41015625 * std::sqrt(165.0);
    const auto f_232 = 2.4609375 * std::sqrt(165.0);
    const auto f_233 = 4.921875 * std::sqrt(165.0);
    const auto f_234 = 29.53125 * std::sqrt(165.0);
    const auto f_235 = 0.73828125 * std::sqrt(165.0);
    const auto f_236 = 4.4296875 * std::sqrt(165.0);
    const auto f_237 = 9.84375 * std::sqrt(165.0);
    const auto f_238 = 59.0625 * std::sqrt(165.0);
    const auto f_239 = 0.08203125 * std::sqrt(165.0);
    const auto f_240 = 0.4921875 * std::sqrt(165.0);
    const auto f_241 = 0.984375 * std::sqrt(165.0);
    const auto f_242 = 5.90625 * std::sqrt(165.0);
    const auto f_243 = 9.84375 * std::sqrt(66.0);
    const auto f_244 = 0.984375 * std::sqrt(66.0);
    const auto f_245 = 16.40625 * std::sqrt(66.0);
    const auto f_246 = 32.8125 * std::sqrt(66.0);
    const auto f_247 = 3.28125 * std::sqrt(66.0);
    const auto f_248 = 7.875 * std::sqrt(165.0);
    const auto f_249 = 26.25 * std::sqrt(165.0);
    const auto f_250 = 0.65625 * std::sqrt(330.0);
    const auto f_251 = 7.875 * std::sqrt(330.0);
    const auto f_252 = 2.625 * std::sqrt(330.0);
    const auto f_253 = 2.1875 * std::sqrt(330.0);
    const auto f_254 = 26.25 * std::sqrt(330.0);
    const auto f_255 = 1.09375 * std::sqrt(330.0);
    const auto f_256 = 8.75 * std::sqrt(330.0);
    const auto f_257 = 15.75 * std::sqrt(55.0);
    const auto f_258 = 26.25 * std::sqrt(55.0);
    const auto f_259 = 52.5 * std::sqrt(55.0);
    const auto f_260 = 0.5625 * std::sqrt(385.0);
    const auto f_261 = 3.375 * std::sqrt(385.0);
    const auto f_262 = 2.25 * std::sqrt(385.0);
    const auto f_263 = 0.9375 * std::sqrt(385.0);
    const auto f_264 = 1.875 * std::sqrt(385.0);
    const auto f_265 = 7.5 * std::sqrt(385.0);
    const auto f_266 = 2.8125 * std::sqrt(231.0);
    const auto f_267 = 0.75 * std::sqrt(231.0);
    const auto f_268 = 4.6875 * std::sqrt(231.0);
    const auto f_269 = 12.5 * std::sqrt(231.0);
    const auto f_270 = 2.5 * std::sqrt(231.0);
    const auto f_271 = 13.125 * std::sqrt(55.0);
    const auto f_272 = 1.96875 * std::sqrt(165.0);
    const auto f_273 = 11.8125 * std::sqrt(165.0);
    const auto f_274 = 6.5625 * std::sqrt(165.0);
    const auto f_275 = 1.845703125 * std::sqrt(6.0);
    const auto f_276 = 3.69140625 * std::sqrt(6.0);
    const auto f_277 = 0.369140625 * std::sqrt(6.0);
    const auto f_278 = 3.076171875 * std::sqrt(6.0);
    const auto f_279 = 6.15234375 * std::sqrt(6.0);
    const auto f_280 = 0.615234375 * std::sqrt(6.0);
    const auto f_281 = 36.9140625 * std::sqrt(6.0);
    const auto f_282 = 73.828125 * std::sqrt(6.0);
    const auto f_283 = 7.3828125 * std::sqrt(6.0);
    const auto f_284 = 1.23046875 * std::sqrt(6.0);
    const auto f_285 = 0.123046875 * std::sqrt(6.0);
    const auto f_286 = 24.609375 * std::sqrt(6.0);
    const auto f_287 = 49.21875 * std::sqrt(6.0);
    const auto f_288 = 4.921875 * std::sqrt(6.0);
    const auto f_289 = 98.4375 * std::sqrt(6.0);
    const auto f_290 = 9.84375 * std::sqrt(6.0);
    const auto f_291 = 12.3046875 * std::sqrt(6.0);
    const auto f_292 = 2.4609375 * std::sqrt(6.0);
    const auto f_293 = 16.40625 * std::sqrt(6.0);
    const auto f_294 = 32.8125 * std::sqrt(6.0);
    const auto f_295 = 3.28125 * std::sqrt(6.0);
    const auto f_296 = 2.953125 * std::sqrt(15.0);
    const auto f_297 = 4.921875 * std::sqrt(15.0);
    const auto f_298 = 59.0625 * std::sqrt(15.0);
    const auto f_299 = 0.984375 * std::sqrt(15.0);
    const auto f_300 = 39.375 * std::sqrt(15.0);
    const auto f_301 = 78.75 * std::sqrt(15.0);
    const auto f_302 = 19.6875 * std::sqrt(15.0);
    const auto f_303 = 26.25 * std::sqrt(15.0);
    const auto f_304 = 0.369140625 * std::sqrt(30.0);
    const auto f_305 = 0.24609375 * std::sqrt(30.0);
    const auto f_306 = 2.953125 * std::sqrt(30.0);
    const auto f_307 = 0.123046875 * std::sqrt(30.0);
    const auto f_308 = 0.984375 * std::sqrt(30.0);
    const auto f_309 = 0.615234375 * std::sqrt(30.0);
    const auto f_310 = 0.41015625 * std::sqrt(30.0);
    const auto f_311 = 4.921875 * std::sqrt(30.0);
    const auto f_312 = 0.205078125 * std::sqrt(30.0);
    const auto f_313 = 1.640625 * std::sqrt(30.0);
    const auto f_314 = 7.3828125 * std::sqrt(30.0);
    const auto f_315 = 59.0625 * std::sqrt(30.0);
    const auto f_316 = 2.4609375 * std::sqrt(30.0);
    const auto f_317 = 19.6875 * std::sqrt(30.0);
    const auto f_318 = 0.08203125 * std::sqrt(30.0);
    const auto f_319 = 0.041015625 * std::sqrt(30.0);
    const auto f_320 = 0.328125 * std::sqrt(30.0);
    const auto f_321 = 3.28125 * std::sqrt(30.0);
    const auto f_322 = 39.375 * std::sqrt(30.0);
    const auto f_323 = 13.125 * std::sqrt(30.0);
    const auto f_324 = 9.84375 * std::sqrt(30.0);
    const auto f_325 = 6.5625 * std::sqrt(30.0);
    const auto f_326 = 78.75 * std::sqrt(30.0);
    const auto f_327 = 26.25 * std::sqrt(30.0);
    const auto f_328 = 0.8203125 * std::sqrt(30.0);
    const auto f_329 = 2.1875 * std::sqrt(30.0);
    const auto f_330 = 1.09375 * std::sqrt(30.0);
    const auto f_331 = 8.75 * std::sqrt(30.0);
    const auto f_332 = 2.953125 * std::sqrt(5.0);
    const auto f_333 = 5.90625 * std::sqrt(5.0);
    const auto f_334 = 4.921875 * std::sqrt(5.0);
    const auto f_335 = 9.84375 * std::sqrt(5.0);
    const auto f_336 = 59.0625 * std::sqrt(5.0);
    const auto f_337 = 118.125 * std::sqrt(5.0);
    const auto f_338 = 0.984375 * std::sqrt(5.0);
    const auto f_339 = 1.96875 * std::sqrt(5.0);
    const auto f_340 = 39.375 * std::sqrt(5.0);
    const auto f_341 = 78.75 * std::sqrt(5.0);
    const auto f_342 = 157.5 * std::sqrt(5.0);
    const auto f_343 = 19.6875 * std::sqrt(5.0);
    const auto f_344 = 26.25 * std::sqrt(5.0);
    const auto f_345 = 52.5 * std::sqrt(5.0);
    const auto f_346 = 0.10546875 * std::sqrt(35.0);
    const auto f_347 = 0.2109375 * std::sqrt(35.0);
    const auto f_348 = 1.265625 * std::sqrt(35.0);
    const auto f_349 = 0.84375 * std::sqrt(35.0);
    const auto f_350 = 0.17578125 * std::sqrt(35.0);
    const auto f_351 = 0.3515625 * std::sqrt(35.0);
    const auto f_352 = 2.109375 * std::sqrt(35.0);
    const auto f_353 = 1.40625 * std::sqrt(35.0);
    const auto f_354 = 4.21875 * std::sqrt(35.0);
    const auto f_355 = 25.3125 * std::sqrt(35.0);
    const auto f_356 = 16.875 * std::sqrt(35.0);
    const auto f_357 = 0.03515625 * std::sqrt(35.0);
    const auto f_358 = 0.0703125 * std::sqrt(35.0);
    const auto f_359 = 0.421875 * std::sqrt(35.0);
    const auto f_360 = 0.28125 * std::sqrt(35.0);
    const auto f_361 = 2.8125 * std::sqrt(35.0);
    const auto f_362 = 11.25 * std::sqrt(35.0);
    const auto f_363 = 5.625 * std::sqrt(35.0);
    const auto f_364 = 33.75 * std::sqrt(35.0);
    const auto f_365 = 22.5 * std::sqrt(35.0);
    const auto f_366 = 0.703125 * std::sqrt(35.0);
    const auto f_367 = 8.4375 * std::sqrt(35.0);
    const auto f_368 = 0.9375 * std::sqrt(35.0);
    const auto f_369 = 1.875 * std::sqrt(35.0);
    const auto f_370 = 7.5 * std::sqrt(35.0);
    const auto f_371 = 0.52734375 * std::sqrt(21.0);
    const auto f_372 = 1.0546875 * std::sqrt(21.0);
    const auto f_373 = 1.40625 * std::sqrt(21.0);
    const auto f_374 = 0.28125 * std::sqrt(21.0);
    const auto f_375 = 0.87890625 * std::sqrt(21.0);
    const auto f_376 = 1.7578125 * std::sqrt(21.0);
    const auto f_377 = 2.34375 * std::sqrt(21.0);
    const auto f_378 = 0.46875 * std::sqrt(21.0);
    const auto f_379 = 10.546875 * std::sqrt(21.0);
    const auto f_380 = 21.09375 * std::sqrt(21.0);
    const auto f_381 = 28.125 * std::sqrt(21.0);
    const auto f_382 = 5.625 * std::sqrt(21.0);
    const auto f_383 = 0.17578125 * std::sqrt(21.0);
    const auto f_384 = 0.3515625 * std::sqrt(21.0);
    const auto f_385 = 0.09375 * std::sqrt(21.0);
    const auto f_386 = 7.03125 * std::sqrt(21.0);
    const auto f_387 = 14.0625 * std::sqrt(21.0);
    const auto f_388 = 18.75 * std::sqrt(21.0);
    const auto f_389 = 3.75 * std::sqrt(21.0);
    const auto f_390 = 37.5 * std::sqrt(21.0);
    const auto f_391 = 7.5 * std::sqrt(21.0);
    const auto f_392 = 3.515625 * std::sqrt(21.0);
    const auto f_393 = 9.375 * std::sqrt(21.0);
    const auto f_394 = 1.875 * std::sqrt(21.0);
    const auto f_395 = 4.6875 * std::sqrt(21.0);
    const auto f_396 = 12.5 * std::sqrt(21.0);
    const auto f_397 = 2.5 * std::sqrt(21.0);
    const auto f_398 = 1.4765625 * std::sqrt(5.0);
    const auto f_399 = 2.4609375 * std::sqrt(5.0);
    const auto f_400 = 29.53125 * std::sqrt(5.0);
    const auto f_401 = 0.4921875 * std::sqrt(5.0);
    const auto f_402 = 13.125 * std::sqrt(5.0);
    const auto f_403 = 0.73828125 * std::sqrt(15.0);
    const auto f_404 = 4.4296875 * std::sqrt(15.0);
    const auto f_405 = 1.23046875 * std::sqrt(15.0);
    const auto f_406 = 7.3828125 * std::sqrt(15.0);
    const auto f_407 = 14.765625 * std::sqrt(15.0);
    const auto f_408 = 88.59375 * std::sqrt(15.0);
    const auto f_409 = 0.24609375 * std::sqrt(15.0);
    const auto f_410 = 1.4765625 * std::sqrt(15.0);
    const auto f_411 = 9.84375 * std::sqrt(15.0);
    const auto f_412 = 118.125 * std::sqrt(15.0);
    const auto f_413 = 29.53125 * std::sqrt(15.0);
    const auto f_414 = 6.5625 * std::sqrt(15.0);
    const auto f_415 = 12.3046875 * std::sqrt(3.0);
    const auto f_416 = 24.609375 * std::sqrt(3.0);
    const auto f_417 = 2.4609375 * std::sqrt(3.0);
    const auto f_418 = 49.21875 * std::sqrt(3.0);
    const auto f_419 = 4.921875 * std::sqrt(3.0);
    const auto f_420 = 65.625 * std::sqrt(3.0);
    const auto f_421 = 131.25 * std::sqrt(3.0);
    const auto f_422 = 13.125 * std::sqrt(3.0);
    const auto f_423 = 39.375 * std::sqrt(3.0);
    const auto f_424 = 78.75 * std::sqrt(3.0);
    const auto f_425 = 7.875 * std::sqrt(3.0);
    const auto f_426 = 52.5 * std::sqrt(30.0);
    const auto f_427 = 31.5 * std::sqrt(30.0);
    const auto f_428 = 2.4609375 * std::sqrt(15.0);
    const auto f_429 = 1.640625 * std::sqrt(15.0);
    const auto f_430 = 0.8203125 * std::sqrt(15.0);
    const auto f_431 = 3.28125 * std::sqrt(15.0);
    const auto f_432 = 13.125 * std::sqrt(15.0);
    const auto f_433 = 8.75 * std::sqrt(15.0);
    const auto f_434 = 105.0 * std::sqrt(15.0);
    const auto f_435 = 4.375 * std::sqrt(15.0);
    const auto f_436 = 35.0 * std::sqrt(15.0);
    const auto f_437 = 7.875 * std::sqrt(15.0);
    const auto f_438 = 5.25 * std::sqrt(15.0);
    const auto f_439 = 63.0 * std::sqrt(15.0);
    const auto f_440 = 2.625 * std::sqrt(15.0);
    const auto f_441 = 21.0 * std::sqrt(15.0);
    const auto f_442 = 9.84375 * std::sqrt(10.0);
    const auto f_443 = 19.6875 * std::sqrt(10.0);
    const auto f_444 = 39.375 * std::sqrt(10.0);
    const auto f_445 = 52.5 * std::sqrt(10.0);
    const auto f_446 = 105.0 * std::sqrt(10.0);
    const auto f_447 = 31.5 * std::sqrt(10.0);
    const auto f_448 = 63.0 * std::sqrt(10.0);
    const auto f_449 = 0.3515625 * std::sqrt(70.0);
    const auto f_450 = 0.703125 * std::sqrt(70.0);
    const auto f_451 = 4.21875 * std::sqrt(70.0);
    const auto f_452 = 2.8125 * std::sqrt(70.0);
    const auto f_453 = 1.40625 * std::sqrt(70.0);
    const auto f_454 = 8.4375 * std::sqrt(70.0);
    const auto f_455 = 5.625 * std::sqrt(70.0);
    const auto f_456 = 1.875 * std::sqrt(70.0);
    const auto f_457 = 3.75 * std::sqrt(70.0);
    const auto f_458 = 22.5 * std::sqrt(70.0);
    const auto f_459 = 15.0 * std::sqrt(70.0);
    const auto f_460 = 1.125 * std::sqrt(70.0);
    const auto f_461 = 2.25 * std::sqrt(70.0);
    const auto f_462 = 13.5 * std::sqrt(70.0);
    const auto f_463 = 9.0 * std::sqrt(70.0);
    const auto f_464 = 1.7578125 * std::sqrt(42.0);
    const auto f_465 = 3.515625 * std::sqrt(42.0);
    const auto f_466 = 4.6875 * std::sqrt(42.0);
    const auto f_467 = 0.9375 * std::sqrt(42.0);
    const auto f_468 = 7.03125 * std::sqrt(42.0);
    const auto f_469 = 9.375 * std::sqrt(42.0);
    const auto f_470 = 1.875 * std::sqrt(42.0);
    const auto f_471 = 18.75 * std::sqrt(42.0);
    const auto f_472 = 25.0 * std::sqrt(42.0);
    const auto f_473 = 5.0 * std::sqrt(42.0);
    const auto f_474 = 5.625 * std::sqrt(42.0);
    const auto f_475 = 11.25 * std::sqrt(42.0);
    const auto f_476 = 15.0 * std::sqrt(42.0);
    const auto f_477 = 3.0 * std::sqrt(42.0);
    const auto f_478 = 4.921875 * std::sqrt(10.0);
    const auto f_479 = 26.25 * std::sqrt(10.0);
    const auto f_480 = 15.75 * std::sqrt(10.0);
    const auto f_481 = 14.765625 * std::sqrt(30.0);
    const auto f_482 = 29.53125 * std::sqrt(30.0);
    const auto f_483 = 7.875 * std::sqrt(30.0);
    const auto f_484 = 47.25 * std::sqrt(30.0);
    const auto f_485 = 1.025390625 * std::sqrt(2.0);
    const auto f_486 = 2.05078125 * std::sqrt(2.0);
    const auto f_487 = 0.205078125 * std::sqrt(2.0);
    const auto f_488 = 3.076171875 * std::sqrt(2.0);
    const auto f_489 = 6.15234375 * std::sqrt(2.0);
    const auto f_490 = 0.615234375 * std::sqrt(2.0);
    const auto f_491 = 24.609375 * std::sqrt(2.0);
    const auto f_492 = 49.21875 * std::sqrt(2.0);
    const auto f_493 = 4.921875 * std::sqrt(2.0);
    const auto f_494 = 98.4375 * std::sqrt(2.0);
    const auto f_495 = 9.84375 * std::sqrt(2.0);
    const auto f_496 = 13.125 * std::sqrt(2.0);
    const auto f_497 = 26.25 * std::sqrt(2.0);
    const auto f_498 = 2.625 * std::sqrt(2.0);
    const auto f_499 = 1.640625 * std::sqrt(5.0);
    const auto f_500 = 21.0 * std::sqrt(5.0);
    const auto f_501 = 0.205078125 * std::sqrt(10.0);
    const auto f_502 = 0.13671875 * std::sqrt(10.0);
    const auto f_503 = 1.640625 * std::sqrt(10.0);
    const auto f_504 = 0.068359375 * std::sqrt(10.0);
    const auto f_505 = 0.546875 * std::sqrt(10.0);
    const auto f_506 = 0.615234375 * std::sqrt(10.0);
    const auto f_507 = 0.41015625 * std::sqrt(10.0);
    const auto f_508 = 3.28125 * std::sqrt(10.0);
    const auto f_509 = 13.125 * std::sqrt(10.0);
    const auto f_510 = 6.5625 * std::sqrt(10.0);
    const auto f_511 = 78.75 * std::sqrt(10.0);
    const auto f_512 = 2.625 * std::sqrt(10.0);
    const auto f_513 = 1.75 * std::sqrt(10.0);
    const auto f_514 = 21.0 * std::sqrt(10.0);
    const auto f_515 = 0.875 * std::sqrt(10.0);
    const auto f_516 = 7.0 * std::sqrt(10.0);
    const auto f_517 = 0.546875 * std::sqrt(15.0);
    const auto f_518 = 1.09375 * std::sqrt(15.0);
    const auto f_519 = 52.5 * std::sqrt(15.0);
    const auto f_520 = 7.0 * std::sqrt(15.0);
    const auto f_521 = 14.0 * std::sqrt(15.0);
    const auto f_522 = 0.01953125 * std::sqrt(105.0);
    const auto f_523 = 0.0390625 * std::sqrt(105.0);
    const auto f_524 = 0.234375 * std::sqrt(105.0);
    const auto f_525 = 0.15625 * std::sqrt(105.0);
    const auto f_526 = 0.05859375 * std::sqrt(105.0);
    const auto f_527 = 0.1171875 * std::sqrt(105.0);
    const auto f_528 = 0.703125 * std::sqrt(105.0);
    const auto f_529 = 0.46875 * std::sqrt(105.0);
    const auto f_530 = 0.9375 * std::sqrt(105.0);
    const auto f_531 = 5.625 * std::sqrt(105.0);
    const auto f_532 = 3.75 * std::sqrt(105.0);
    const auto f_533 = 1.875 * std::sqrt(105.0);
    const auto f_534 = 11.25 * std::sqrt(105.0);
    const auto f_535 = 7.5 * std::sqrt(105.0);
    const auto f_536 = 0.25 * std::sqrt(105.0);
    const auto f_537 = 0.5 * std::sqrt(105.0);
    const auto f_538 = 3.0 * std::sqrt(105.0);
    const auto f_539 = 2.0 * std::sqrt(105.0);
    const auto f_540 = 0.29296875 * std::sqrt(7.0);
    const auto f_541 = 0.5859375 * std::sqrt(7.0);
    const auto f_542 = 0.78125 * std::sqrt(7.0);
    const auto f_543 = 0.15625 * std::sqrt(7.0);
    const auto f_544 = 0.87890625 * std::sqrt(7.0);
    const auto f_545 = 1.7578125 * std::sqrt(7.0);
    const auto f_546 = 2.34375 * std::sqrt(7.0);
    const auto f_547 = 0.46875 * std::sqrt(7.0);
    const auto f_548 = 7.03125 * std::sqrt(7.0);
    const auto f_549 = 14.0625 * std::sqrt(7.0);
    const auto f_550 = 18.75 * std::sqrt(7.0);
    const auto f_551 = 3.75 * std::sqrt(7.0);
    const auto f_552 = 28.125 * std::sqrt(7.0);
    const auto f_553 = 37.5 * std::sqrt(7.0);
    const auto f_554 = 7.5 * std::sqrt(7.0);
    const auto f_555 = 10.0 * std::sqrt(7.0);
    const auto f_556 = 2.0 * std::sqrt(7.0);
    const auto f_557 = 0.2734375 * std::sqrt(15.0);
    const auto f_558 = 3.5 * std::sqrt(15.0);
    const auto f_559 = 0.41015625 * std::sqrt(5.0);
    const auto f_560 = 1.23046875 * std::sqrt(5.0);
    const auto f_561 = 7.3828125 * std::sqrt(5.0);
    const auto f_562 = 5.25 * std::sqrt(5.0);
    const auto f_563 = 31.5 * std::sqrt(5.0);
    const auto f_564 = 2.05078125 * std::sqrt(14.0);
    const auto f_565 = 4.1015625 * std::sqrt(14.0);
    const auto f_566 = 0.41015625 * std::sqrt(14.0);
    const auto f_567 = 6.15234375 * std::sqrt(14.0);
    const auto f_568 = 12.3046875 * std::sqrt(14.0);
    const auto f_569 = 1.23046875 * std::sqrt(14.0);
    const auto f_570 = 24.609375 * std::sqrt(14.0);
    const auto f_571 = 2.4609375 * std::sqrt(14.0);
    const auto f_572 = 49.21875 * std::sqrt(14.0);
    const auto f_573 = 4.921875 * std::sqrt(14.0);
    const auto f_574 = 9.84375 * std::sqrt(14.0);
    const auto f_575 = 19.6875 * std::sqrt(14.0);
    const auto f_576 = 1.96875 * std::sqrt(14.0);
    const auto f_577 = 0.9375 * std::sqrt(14.0);
    const auto f_578 = 1.875 * std::sqrt(14.0);
    const auto f_579 = 0.1875 * std::sqrt(14.0);
    const auto f_580 = 3.28125 * std::sqrt(35.0);
    const auto f_581 = 9.84375 * std::sqrt(35.0);
    const auto f_582 = 19.6875 * std::sqrt(35.0);
    const auto f_583 = 39.375 * std::sqrt(35.0);
    const auto f_584 = 15.75 * std::sqrt(35.0);
    const auto f_585 = 1.5 * std::sqrt(35.0);
    const auto f_586 = 0.41015625 * std::sqrt(70.0);
    const auto f_587 = 0.2734375 * std::sqrt(70.0);
    const auto f_588 = 3.28125 * std::sqrt(70.0);
    const auto f_589 = 0.13671875 * std::sqrt(70.0);
    const auto f_590 = 1.09375 * std::sqrt(70.0);
    const auto f_591 = 1.23046875 * std::sqrt(70.0);
    const auto f_592 = 0.8203125 * std::sqrt(70.0);
    const auto f_593 = 9.84375 * std::sqrt(70.0);
    const auto f_594 = 2.4609375 * std::sqrt(70.0);
    const auto f_595 = 1.640625 * std::sqrt(70.0);
    const auto f_596 = 19.6875 * std::sqrt(70.0);
    const auto f_597 = 6.5625 * std::sqrt(70.0);
    const auto f_598 = 4.921875 * std::sqrt(70.0);
    const auto f_599 = 39.375 * std::sqrt(70.0);
    const auto f_600 = 13.125 * std::sqrt(70.0);
    const auto f_601 = 1.96875 * std::sqrt(70.0);
    const auto f_602 = 1.3125 * std::sqrt(70.0);
    const auto f_603 = 15.75 * std::sqrt(70.0);
    const auto f_604 = 0.65625 * std::sqrt(70.0);
    const auto f_605 = 5.25 * std::sqrt(70.0);
    const auto f_606 = 0.1875 * std::sqrt(70.0);
    const auto f_607 = 0.125 * std::sqrt(70.0);
    const auto f_608 = 1.5 * std::sqrt(70.0);
    const auto f_609 = 0.0625 * std::sqrt(70.0);
    const auto f_610 = 0.5 * std::sqrt(70.0);
    const auto f_611 = 1.09375 * std::sqrt(105.0);
    const auto f_612 = 2.1875 * std::sqrt(105.0);
    const auto f_613 = 3.28125 * std::sqrt(105.0);
    const auto f_614 = 6.5625 * std::sqrt(105.0);
    const auto f_615 = 13.125 * std::sqrt(105.0);
    const auto f_616 = 26.25 * std::sqrt(105.0);
    const auto f_617 = 5.25 * std::sqrt(105.0);
    const auto f_618 = 10.5 * std::sqrt(105.0);
    const auto f_619 = std::sqrt(105.0);
    const auto f_620 = 2.1875 * std::sqrt(15.0);
    const auto f_621 = 1.3125 * std::sqrt(15.0);
    const auto f_622 = 15.75 * std::sqrt(15.0);
    const auto f_623 = 10.5 * std::sqrt(15.0);
    const auto f_624 = 0.125 * std::sqrt(15.0);
    const auto f_625 = 0.25 * std::sqrt(15.0);
    const auto f_626 = 1.5 * std::sqrt(15.0);
    const auto f_627 = std::sqrt(15.0);
    const auto f_628 = 0.546875 * std::sqrt(105.0);
    const auto f_629 = 1.640625 * std::sqrt(105.0);
    const auto f_630 = 2.625 * std::sqrt(105.0);
    const auto f_631 = 0.8203125 * std::sqrt(35.0);
    const auto f_632 = 4.921875 * std::sqrt(35.0);
    const auto f_633 = 2.4609375 * std::sqrt(35.0);
    const auto f_634 = 14.765625 * std::sqrt(35.0);
    const auto f_635 = 29.53125 * std::sqrt(35.0);
    const auto f_636 = 59.0625 * std::sqrt(35.0);
    const auto f_637 = 3.9375 * std::sqrt(35.0);
    const auto f_638 = 23.625 * std::sqrt(35.0);
    const auto f_639 = 0.375 * std::sqrt(35.0);
    const auto f_640 = 2.25 * std::sqrt(35.0);
    const auto f_641 = 6.15234375 * std::sqrt(3.0);
    const auto f_642 = 1.23046875 * std::sqrt(3.0);
    const auto f_643 = 32.8125 * std::sqrt(3.0);
    const auto f_644 = 6.5625 * std::sqrt(3.0);
    const auto f_645 = 19.6875 * std::sqrt(3.0);
    const auto f_646 = 3.9375 * std::sqrt(3.0);
    const auto f_647 = 15.75 * std::sqrt(30.0);
    const auto f_648 = 0.41015625 * std::sqrt(15.0);
    const auto f_649 = 17.5 * std::sqrt(15.0);
    const auto f_650 = 3.9375 * std::sqrt(15.0);
    const auto f_651 = 31.5 * std::sqrt(15.0);
    const auto f_652 = 0.17578125 * std::sqrt(70.0);
    const auto f_653 = 2.109375 * std::sqrt(70.0);
    const auto f_654 = 0.9375 * std::sqrt(70.0);
    const auto f_655 = 11.25 * std::sqrt(70.0);
    const auto f_656 = 7.5 * std::sqrt(70.0);
    const auto f_657 = 0.5625 * std::sqrt(70.0);
    const auto f_658 = 6.75 * std::sqrt(70.0);
    const auto f_659 = 4.5 * std::sqrt(70.0);
    const auto f_660 = 0.87890625 * std::sqrt(42.0);
    const auto f_661 = 2.34375 * std::sqrt(42.0);
    const auto f_662 = 0.46875 * std::sqrt(42.0);
    const auto f_663 = 12.5 * std::sqrt(42.0);
    const auto f_664 = 2.5 * std::sqrt(42.0);
    const auto f_665 = 2.8125 * std::sqrt(42.0);
    const auto f_666 = 7.5 * std::sqrt(42.0);
    const auto f_667 = 1.5 * std::sqrt(42.0);
    const auto f_668 = 2.4609375 * std::sqrt(10.0);
    const auto f_669 = 7.875 * std::sqrt(10.0);
    const auto f_670 = 1.23046875 * std::sqrt(30.0);
    const auto f_671 = 3.9375 * std::sqrt(30.0);
    const auto f_672 = 23.625 * std::sqrt(30.0);
    const auto f_673 = 1.23046875 * std::sqrt(66.0);
    const auto f_674 = 0.24609375 * std::sqrt(66.0);
    const auto f_675 = 6.15234375 * std::sqrt(66.0);
    const auto f_676 = 4.1015625 * std::sqrt(66.0);
    const auto f_677 = 8.203125 * std::sqrt(66.0);
    const auto f_678 = 0.8203125 * std::sqrt(66.0);
    const auto f_679 = 1.96875 * std::sqrt(330.0);
    const auto f_680 = 0.08203125 * std::sqrt(330.0);
    const auto f_681 = 1.23046875 * std::sqrt(330.0);
    const auto f_682 = 9.84375 * std::sqrt(330.0);
    const auto f_683 = 0.41015625 * std::sqrt(330.0);
    const auto f_684 = 0.2734375 * std::sqrt(330.0);
    const auto f_685 = 6.5625 * std::sqrt(55.0);
    const auto f_686 = 0.0703125 * std::sqrt(385.0);
    const auto f_687 = 0.3515625 * std::sqrt(385.0);
    const auto f_688 = 4.21875 * std::sqrt(385.0);
    const auto f_689 = 0.234375 * std::sqrt(385.0);
    const auto f_690 = 0.3515625 * std::sqrt(231.0);
    const auto f_691 = 0.9375 * std::sqrt(231.0);
    const auto f_692 = 0.1875 * std::sqrt(231.0);
    const auto f_693 = 1.7578125 * std::sqrt(231.0);
    const auto f_694 = 1.171875 * std::sqrt(231.0);
    const auto f_695 = 2.34375 * std::sqrt(231.0);
    const auto f_696 = 3.125 * std::sqrt(231.0);
    const auto f_697 = 0.625 * std::sqrt(231.0);
    const auto f_698 = 0.984375 * std::sqrt(55.0);
    const auto f_699 = 4.921875 * std::sqrt(55.0);
    const auto f_700 = 14.765625 * std::sqrt(165.0);
    const auto f_701 = 0.08203125 * std::sqrt(429.0);
    const auto f_702 = 6.15234375 * std::sqrt(429.0);
    const auto f_703 = 12.3046875 * std::sqrt(429.0);
    const auto f_704 = 0.328125 * std::sqrt(4290.0);
    const auto f_705 = 4.921875 * std::sqrt(4290.0);
    const auto f_706 = 0.08203125 * std::sqrt(2145.0);
    const auto f_707 = 0.0546875 * std::sqrt(2145.0);
    const auto f_708 = 0.65625 * std::sqrt(2145.0);
    const auto f_709 = 0.02734375 * std::sqrt(2145.0);
    const auto f_710 = 0.21875 * std::sqrt(2145.0);
    const auto f_711 = 1.23046875 * std::sqrt(2145.0);
    const auto f_712 = 0.8203125 * std::sqrt(2145.0);
    const auto f_713 = 9.84375 * std::sqrt(2145.0);
    const auto f_714 = 0.41015625 * std::sqrt(2145.0);
    const auto f_715 = 3.28125 * std::sqrt(2145.0);
    const auto f_716 = 0.328125 * std::sqrt(1430.0);
    const auto f_717 = 0.65625 * std::sqrt(1430.0);
    const auto f_718 = 4.921875 * std::sqrt(1430.0);
    const auto f_719 = 9.84375 * std::sqrt(1430.0);
    const auto f_720 = 0.01171875 * std::sqrt(10010.0);
    const auto f_721 = 0.0234375 * std::sqrt(10010.0);
    const auto f_722 = 0.09375 * std::sqrt(10010.0);
    const auto f_723 = 0.17578125 * std::sqrt(10010.0);
    const auto f_724 = 0.3515625 * std::sqrt(10010.0);
    const auto f_725 = 2.109375 * std::sqrt(10010.0);
    const auto f_726 = 1.40625 * std::sqrt(10010.0);
    const auto f_727 = 0.1171875 * std::sqrt(6006.0);
    const auto f_728 = 0.15625 * std::sqrt(6006.0);
    const auto f_729 = 0.03125 * std::sqrt(6006.0);
    const auto f_730 = 0.87890625 * std::sqrt(6006.0);
    const auto f_731 = 1.7578125 * std::sqrt(6006.0);
    const auto f_732 = 0.46875 * std::sqrt(6006.0);
    const auto f_733 = 0.1640625 * std::sqrt(1430.0);
    const auto f_734 = 2.4609375 * std::sqrt(1430.0);
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

    const auto *kh_0 = buffer.data(kh + 0);
    const auto *kh_1 = buffer.data(kh + 1);
    const auto *kh_2 = buffer.data(kh + 2);
    const auto *kh_3 = buffer.data(kh + 3);
    const auto *kh_4 = buffer.data(kh + 4);
    const auto *kh_5 = buffer.data(kh + 5);
    const auto *kh_6 = buffer.data(kh + 6);
    const auto *kh_7 = buffer.data(kh + 7);
    const auto *kh_8 = buffer.data(kh + 8);
    const auto *kh_9 = buffer.data(kh + 9);
    const auto *kh_10 = buffer.data(kh + 10);
    const auto *kh_11 = buffer.data(kh + 11);
    const auto *kh_12 = buffer.data(kh + 12);
    const auto *kh_13 = buffer.data(kh + 13);
    const auto *kh_14 = buffer.data(kh + 14);
    const auto *kh_15 = buffer.data(kh + 15);
    const auto *kh_16 = buffer.data(kh + 16);
    const auto *kh_17 = buffer.data(kh + 17);
    const auto *kh_18 = buffer.data(kh + 18);
    const auto *kh_19 = buffer.data(kh + 19);
    const auto *kh_20 = buffer.data(kh + 20);
    const auto *kh_21 = buffer.data(kh + 21);
    const auto *kh_22 = buffer.data(kh + 22);
    const auto *kh_23 = buffer.data(kh + 23);
    const auto *kh_24 = buffer.data(kh + 24);
    const auto *kh_25 = buffer.data(kh + 25);
    const auto *kh_26 = buffer.data(kh + 26);
    const auto *kh_27 = buffer.data(kh + 27);
    const auto *kh_28 = buffer.data(kh + 28);
    const auto *kh_29 = buffer.data(kh + 29);
    const auto *kh_30 = buffer.data(kh + 30);
    const auto *kh_31 = buffer.data(kh + 31);
    const auto *kh_32 = buffer.data(kh + 32);
    const auto *kh_33 = buffer.data(kh + 33);
    const auto *kh_34 = buffer.data(kh + 34);
    const auto *kh_35 = buffer.data(kh + 35);
    const auto *kh_36 = buffer.data(kh + 36);
    const auto *kh_37 = buffer.data(kh + 37);
    const auto *kh_38 = buffer.data(kh + 38);
    const auto *kh_39 = buffer.data(kh + 39);
    const auto *kh_40 = buffer.data(kh + 40);
    const auto *kh_41 = buffer.data(kh + 41);
    const auto *kh_42 = buffer.data(kh + 42);
    const auto *kh_43 = buffer.data(kh + 43);
    const auto *kh_44 = buffer.data(kh + 44);
    const auto *kh_45 = buffer.data(kh + 45);
    const auto *kh_46 = buffer.data(kh + 46);
    const auto *kh_47 = buffer.data(kh + 47);
    const auto *kh_48 = buffer.data(kh + 48);
    const auto *kh_49 = buffer.data(kh + 49);
    const auto *kh_50 = buffer.data(kh + 50);
    const auto *kh_51 = buffer.data(kh + 51);
    const auto *kh_52 = buffer.data(kh + 52);
    const auto *kh_53 = buffer.data(kh + 53);
    const auto *kh_54 = buffer.data(kh + 54);
    const auto *kh_55 = buffer.data(kh + 55);
    const auto *kh_56 = buffer.data(kh + 56);
    const auto *kh_57 = buffer.data(kh + 57);
    const auto *kh_58 = buffer.data(kh + 58);
    const auto *kh_59 = buffer.data(kh + 59);
    const auto *kh_60 = buffer.data(kh + 60);
    const auto *kh_61 = buffer.data(kh + 61);
    const auto *kh_62 = buffer.data(kh + 62);
    const auto *kh_63 = buffer.data(kh + 63);
    const auto *kh_64 = buffer.data(kh + 64);
    const auto *kh_65 = buffer.data(kh + 65);
    const auto *kh_66 = buffer.data(kh + 66);
    const auto *kh_67 = buffer.data(kh + 67);
    const auto *kh_68 = buffer.data(kh + 68);
    const auto *kh_69 = buffer.data(kh + 69);
    const auto *kh_70 = buffer.data(kh + 70);
    const auto *kh_71 = buffer.data(kh + 71);
    const auto *kh_72 = buffer.data(kh + 72);
    const auto *kh_73 = buffer.data(kh + 73);
    const auto *kh_74 = buffer.data(kh + 74);
    const auto *kh_75 = buffer.data(kh + 75);
    const auto *kh_76 = buffer.data(kh + 76);
    const auto *kh_77 = buffer.data(kh + 77);
    const auto *kh_78 = buffer.data(kh + 78);
    const auto *kh_79 = buffer.data(kh + 79);
    const auto *kh_80 = buffer.data(kh + 80);
    const auto *kh_81 = buffer.data(kh + 81);
    const auto *kh_82 = buffer.data(kh + 82);
    const auto *kh_83 = buffer.data(kh + 83);
    const auto *kh_84 = buffer.data(kh + 84);
    const auto *kh_85 = buffer.data(kh + 85);
    const auto *kh_86 = buffer.data(kh + 86);
    const auto *kh_87 = buffer.data(kh + 87);
    const auto *kh_88 = buffer.data(kh + 88);
    const auto *kh_89 = buffer.data(kh + 89);
    const auto *kh_90 = buffer.data(kh + 90);
    const auto *kh_91 = buffer.data(kh + 91);
    const auto *kh_92 = buffer.data(kh + 92);
    const auto *kh_93 = buffer.data(kh + 93);
    const auto *kh_94 = buffer.data(kh + 94);
    const auto *kh_95 = buffer.data(kh + 95);
    const auto *kh_96 = buffer.data(kh + 96);
    const auto *kh_97 = buffer.data(kh + 97);
    const auto *kh_98 = buffer.data(kh + 98);
    const auto *kh_99 = buffer.data(kh + 99);
    const auto *kh_100 = buffer.data(kh + 100);
    const auto *kh_101 = buffer.data(kh + 101);
    const auto *kh_102 = buffer.data(kh + 102);
    const auto *kh_103 = buffer.data(kh + 103);
    const auto *kh_104 = buffer.data(kh + 104);
    const auto *kh_105 = buffer.data(kh + 105);
    const auto *kh_106 = buffer.data(kh + 106);
    const auto *kh_107 = buffer.data(kh + 107);
    const auto *kh_108 = buffer.data(kh + 108);
    const auto *kh_109 = buffer.data(kh + 109);
    const auto *kh_110 = buffer.data(kh + 110);
    const auto *kh_111 = buffer.data(kh + 111);
    const auto *kh_112 = buffer.data(kh + 112);
    const auto *kh_113 = buffer.data(kh + 113);
    const auto *kh_114 = buffer.data(kh + 114);
    const auto *kh_115 = buffer.data(kh + 115);
    const auto *kh_116 = buffer.data(kh + 116);
    const auto *kh_117 = buffer.data(kh + 117);
    const auto *kh_118 = buffer.data(kh + 118);
    const auto *kh_119 = buffer.data(kh + 119);
    const auto *kh_120 = buffer.data(kh + 120);
    const auto *kh_121 = buffer.data(kh + 121);
    const auto *kh_122 = buffer.data(kh + 122);
    const auto *kh_123 = buffer.data(kh + 123);
    const auto *kh_124 = buffer.data(kh + 124);
    const auto *kh_125 = buffer.data(kh + 125);
    const auto *kh_126 = buffer.data(kh + 126);
    const auto *kh_127 = buffer.data(kh + 127);
    const auto *kh_128 = buffer.data(kh + 128);
    const auto *kh_129 = buffer.data(kh + 129);
    const auto *kh_130 = buffer.data(kh + 130);
    const auto *kh_131 = buffer.data(kh + 131);
    const auto *kh_132 = buffer.data(kh + 132);
    const auto *kh_133 = buffer.data(kh + 133);
    const auto *kh_134 = buffer.data(kh + 134);
    const auto *kh_135 = buffer.data(kh + 135);
    const auto *kh_136 = buffer.data(kh + 136);
    const auto *kh_137 = buffer.data(kh + 137);
    const auto *kh_138 = buffer.data(kh + 138);
    const auto *kh_139 = buffer.data(kh + 139);
    const auto *kh_140 = buffer.data(kh + 140);
    const auto *kh_141 = buffer.data(kh + 141);
    const auto *kh_142 = buffer.data(kh + 142);
    const auto *kh_143 = buffer.data(kh + 143);
    const auto *kh_144 = buffer.data(kh + 144);
    const auto *kh_145 = buffer.data(kh + 145);
    const auto *kh_146 = buffer.data(kh + 146);
    const auto *kh_147 = buffer.data(kh + 147);
    const auto *kh_148 = buffer.data(kh + 148);
    const auto *kh_149 = buffer.data(kh + 149);
    const auto *kh_150 = buffer.data(kh + 150);
    const auto *kh_151 = buffer.data(kh + 151);
    const auto *kh_152 = buffer.data(kh + 152);
    const auto *kh_153 = buffer.data(kh + 153);
    const auto *kh_154 = buffer.data(kh + 154);
    const auto *kh_155 = buffer.data(kh + 155);
    const auto *kh_156 = buffer.data(kh + 156);
    const auto *kh_157 = buffer.data(kh + 157);
    const auto *kh_158 = buffer.data(kh + 158);
    const auto *kh_159 = buffer.data(kh + 159);
    const auto *kh_160 = buffer.data(kh + 160);
    const auto *kh_161 = buffer.data(kh + 161);
    const auto *kh_162 = buffer.data(kh + 162);
    const auto *kh_163 = buffer.data(kh + 163);
    const auto *kh_164 = buffer.data(kh + 164);
    const auto *kh_165 = buffer.data(kh + 165);
    const auto *kh_166 = buffer.data(kh + 166);
    const auto *kh_167 = buffer.data(kh + 167);
    const auto *kh_168 = buffer.data(kh + 168);
    const auto *kh_169 = buffer.data(kh + 169);
    const auto *kh_170 = buffer.data(kh + 170);
    const auto *kh_171 = buffer.data(kh + 171);
    const auto *kh_172 = buffer.data(kh + 172);
    const auto *kh_173 = buffer.data(kh + 173);
    const auto *kh_174 = buffer.data(kh + 174);
    const auto *kh_175 = buffer.data(kh + 175);
    const auto *kh_176 = buffer.data(kh + 176);
    const auto *kh_177 = buffer.data(kh + 177);
    const auto *kh_178 = buffer.data(kh + 178);
    const auto *kh_179 = buffer.data(kh + 179);
    const auto *kh_180 = buffer.data(kh + 180);
    const auto *kh_181 = buffer.data(kh + 181);
    const auto *kh_182 = buffer.data(kh + 182);
    const auto *kh_183 = buffer.data(kh + 183);
    const auto *kh_184 = buffer.data(kh + 184);
    const auto *kh_185 = buffer.data(kh + 185);
    const auto *kh_186 = buffer.data(kh + 186);
    const auto *kh_187 = buffer.data(kh + 187);
    const auto *kh_188 = buffer.data(kh + 188);
    const auto *kh_189 = buffer.data(kh + 189);
    const auto *kh_190 = buffer.data(kh + 190);
    const auto *kh_191 = buffer.data(kh + 191);
    const auto *kh_192 = buffer.data(kh + 192);
    const auto *kh_193 = buffer.data(kh + 193);
    const auto *kh_194 = buffer.data(kh + 194);
    const auto *kh_195 = buffer.data(kh + 195);
    const auto *kh_196 = buffer.data(kh + 196);
    const auto *kh_197 = buffer.data(kh + 197);
    const auto *kh_198 = buffer.data(kh + 198);
    const auto *kh_199 = buffer.data(kh + 199);
    const auto *kh_200 = buffer.data(kh + 200);
    const auto *kh_201 = buffer.data(kh + 201);
    const auto *kh_202 = buffer.data(kh + 202);
    const auto *kh_203 = buffer.data(kh + 203);
    const auto *kh_204 = buffer.data(kh + 204);
    const auto *kh_205 = buffer.data(kh + 205);
    const auto *kh_206 = buffer.data(kh + 206);
    const auto *kh_207 = buffer.data(kh + 207);
    const auto *kh_208 = buffer.data(kh + 208);
    const auto *kh_209 = buffer.data(kh + 209);
    const auto *kh_210 = buffer.data(kh + 210);
    const auto *kh_211 = buffer.data(kh + 211);
    const auto *kh_212 = buffer.data(kh + 212);
    const auto *kh_213 = buffer.data(kh + 213);
    const auto *kh_214 = buffer.data(kh + 214);
    const auto *kh_215 = buffer.data(kh + 215);
    const auto *kh_216 = buffer.data(kh + 216);
    const auto *kh_217 = buffer.data(kh + 217);
    const auto *kh_218 = buffer.data(kh + 218);
    const auto *kh_219 = buffer.data(kh + 219);
    const auto *kh_220 = buffer.data(kh + 220);
    const auto *kh_221 = buffer.data(kh + 221);
    const auto *kh_222 = buffer.data(kh + 222);
    const auto *kh_223 = buffer.data(kh + 223);
    const auto *kh_224 = buffer.data(kh + 224);
    const auto *kh_225 = buffer.data(kh + 225);
    const auto *kh_226 = buffer.data(kh + 226);
    const auto *kh_227 = buffer.data(kh + 227);
    const auto *kh_228 = buffer.data(kh + 228);
    const auto *kh_229 = buffer.data(kh + 229);
    const auto *kh_230 = buffer.data(kh + 230);
    const auto *kh_231 = buffer.data(kh + 231);
    const auto *kh_232 = buffer.data(kh + 232);
    const auto *kh_233 = buffer.data(kh + 233);
    const auto *kh_234 = buffer.data(kh + 234);
    const auto *kh_235 = buffer.data(kh + 235);
    const auto *kh_236 = buffer.data(kh + 236);
    const auto *kh_237 = buffer.data(kh + 237);
    const auto *kh_238 = buffer.data(kh + 238);
    const auto *kh_239 = buffer.data(kh + 239);
    const auto *kh_240 = buffer.data(kh + 240);
    const auto *kh_241 = buffer.data(kh + 241);
    const auto *kh_242 = buffer.data(kh + 242);
    const auto *kh_243 = buffer.data(kh + 243);
    const auto *kh_244 = buffer.data(kh + 244);
    const auto *kh_245 = buffer.data(kh + 245);
    const auto *kh_246 = buffer.data(kh + 246);
    const auto *kh_247 = buffer.data(kh + 247);
    const auto *kh_248 = buffer.data(kh + 248);
    const auto *kh_249 = buffer.data(kh + 249);
    const auto *kh_250 = buffer.data(kh + 250);
    const auto *kh_251 = buffer.data(kh + 251);
    const auto *kh_252 = buffer.data(kh + 252);
    const auto *kh_253 = buffer.data(kh + 253);
    const auto *kh_254 = buffer.data(kh + 254);
    const auto *kh_255 = buffer.data(kh + 255);
    const auto *kh_256 = buffer.data(kh + 256);
    const auto *kh_257 = buffer.data(kh + 257);
    const auto *kh_258 = buffer.data(kh + 258);
    const auto *kh_259 = buffer.data(kh + 259);
    const auto *kh_260 = buffer.data(kh + 260);
    const auto *kh_261 = buffer.data(kh + 261);
    const auto *kh_262 = buffer.data(kh + 262);
    const auto *kh_263 = buffer.data(kh + 263);
    const auto *kh_264 = buffer.data(kh + 264);
    const auto *kh_265 = buffer.data(kh + 265);
    const auto *kh_266 = buffer.data(kh + 266);
    const auto *kh_267 = buffer.data(kh + 267);
    const auto *kh_268 = buffer.data(kh + 268);
    const auto *kh_269 = buffer.data(kh + 269);
    const auto *kh_270 = buffer.data(kh + 270);
    const auto *kh_271 = buffer.data(kh + 271);
    const auto *kh_272 = buffer.data(kh + 272);
    const auto *kh_273 = buffer.data(kh + 273);
    const auto *kh_274 = buffer.data(kh + 274);
    const auto *kh_275 = buffer.data(kh + 275);
    const auto *kh_276 = buffer.data(kh + 276);
    const auto *kh_277 = buffer.data(kh + 277);
    const auto *kh_278 = buffer.data(kh + 278);
    const auto *kh_279 = buffer.data(kh + 279);
    const auto *kh_280 = buffer.data(kh + 280);
    const auto *kh_281 = buffer.data(kh + 281);
    const auto *kh_282 = buffer.data(kh + 282);
    const auto *kh_283 = buffer.data(kh + 283);
    const auto *kh_284 = buffer.data(kh + 284);
    const auto *kh_285 = buffer.data(kh + 285);
    const auto *kh_286 = buffer.data(kh + 286);
    const auto *kh_287 = buffer.data(kh + 287);
    const auto *kh_288 = buffer.data(kh + 288);
    const auto *kh_289 = buffer.data(kh + 289);
    const auto *kh_290 = buffer.data(kh + 290);
    const auto *kh_291 = buffer.data(kh + 291);
    const auto *kh_292 = buffer.data(kh + 292);
    const auto *kh_293 = buffer.data(kh + 293);
    const auto *kh_294 = buffer.data(kh + 294);
    const auto *kh_295 = buffer.data(kh + 295);
    const auto *kh_296 = buffer.data(kh + 296);
    const auto *kh_297 = buffer.data(kh + 297);
    const auto *kh_298 = buffer.data(kh + 298);
    const auto *kh_299 = buffer.data(kh + 299);
    const auto *kh_300 = buffer.data(kh + 300);
    const auto *kh_301 = buffer.data(kh + 301);
    const auto *kh_302 = buffer.data(kh + 302);
    const auto *kh_303 = buffer.data(kh + 303);
    const auto *kh_304 = buffer.data(kh + 304);
    const auto *kh_305 = buffer.data(kh + 305);
    const auto *kh_306 = buffer.data(kh + 306);
    const auto *kh_307 = buffer.data(kh + 307);
    const auto *kh_308 = buffer.data(kh + 308);
    const auto *kh_309 = buffer.data(kh + 309);
    const auto *kh_310 = buffer.data(kh + 310);
    const auto *kh_311 = buffer.data(kh + 311);
    const auto *kh_312 = buffer.data(kh + 312);
    const auto *kh_313 = buffer.data(kh + 313);
    const auto *kh_314 = buffer.data(kh + 314);
    const auto *kh_315 = buffer.data(kh + 315);
    const auto *kh_316 = buffer.data(kh + 316);
    const auto *kh_317 = buffer.data(kh + 317);
    const auto *kh_318 = buffer.data(kh + 318);
    const auto *kh_319 = buffer.data(kh + 319);
    const auto *kh_320 = buffer.data(kh + 320);
    const auto *kh_321 = buffer.data(kh + 321);
    const auto *kh_322 = buffer.data(kh + 322);
    const auto *kh_323 = buffer.data(kh + 323);
    const auto *kh_324 = buffer.data(kh + 324);
    const auto *kh_325 = buffer.data(kh + 325);
    const auto *kh_326 = buffer.data(kh + 326);
    const auto *kh_327 = buffer.data(kh + 327);
    const auto *kh_328 = buffer.data(kh + 328);
    const auto *kh_329 = buffer.data(kh + 329);
    const auto *kh_330 = buffer.data(kh + 330);
    const auto *kh_331 = buffer.data(kh + 331);
    const auto *kh_332 = buffer.data(kh + 332);
    const auto *kh_333 = buffer.data(kh + 333);
    const auto *kh_334 = buffer.data(kh + 334);
    const auto *kh_335 = buffer.data(kh + 335);
    const auto *kh_336 = buffer.data(kh + 336);
    const auto *kh_337 = buffer.data(kh + 337);
    const auto *kh_338 = buffer.data(kh + 338);
    const auto *kh_339 = buffer.data(kh + 339);
    const auto *kh_340 = buffer.data(kh + 340);
    const auto *kh_341 = buffer.data(kh + 341);
    const auto *kh_342 = buffer.data(kh + 342);
    const auto *kh_343 = buffer.data(kh + 343);
    const auto *kh_344 = buffer.data(kh + 344);
    const auto *kh_345 = buffer.data(kh + 345);
    const auto *kh_346 = buffer.data(kh + 346);
    const auto *kh_347 = buffer.data(kh + 347);
    const auto *kh_348 = buffer.data(kh + 348);
    const auto *kh_349 = buffer.data(kh + 349);
    const auto *kh_350 = buffer.data(kh + 350);
    const auto *kh_351 = buffer.data(kh + 351);
    const auto *kh_352 = buffer.data(kh + 352);
    const auto *kh_353 = buffer.data(kh + 353);
    const auto *kh_354 = buffer.data(kh + 354);
    const auto *kh_355 = buffer.data(kh + 355);
    const auto *kh_356 = buffer.data(kh + 356);
    const auto *kh_357 = buffer.data(kh + 357);
    const auto *kh_358 = buffer.data(kh + 358);
    const auto *kh_359 = buffer.data(kh + 359);
    const auto *kh_360 = buffer.data(kh + 360);
    const auto *kh_361 = buffer.data(kh + 361);
    const auto *kh_362 = buffer.data(kh + 362);
    const auto *kh_363 = buffer.data(kh + 363);
    const auto *kh_364 = buffer.data(kh + 364);
    const auto *kh_365 = buffer.data(kh + 365);
    const auto *kh_366 = buffer.data(kh + 366);
    const auto *kh_367 = buffer.data(kh + 367);
    const auto *kh_368 = buffer.data(kh + 368);
    const auto *kh_369 = buffer.data(kh + 369);
    const auto *kh_370 = buffer.data(kh + 370);
    const auto *kh_371 = buffer.data(kh + 371);
    const auto *kh_372 = buffer.data(kh + 372);
    const auto *kh_373 = buffer.data(kh + 373);
    const auto *kh_374 = buffer.data(kh + 374);
    const auto *kh_375 = buffer.data(kh + 375);
    const auto *kh_376 = buffer.data(kh + 376);
    const auto *kh_377 = buffer.data(kh + 377);
    const auto *kh_378 = buffer.data(kh + 378);
    const auto *kh_379 = buffer.data(kh + 379);
    const auto *kh_380 = buffer.data(kh + 380);
    const auto *kh_381 = buffer.data(kh + 381);
    const auto *kh_382 = buffer.data(kh + 382);
    const auto *kh_383 = buffer.data(kh + 383);
    const auto *kh_384 = buffer.data(kh + 384);
    const auto *kh_385 = buffer.data(kh + 385);
    const auto *kh_386 = buffer.data(kh + 386);
    const auto *kh_387 = buffer.data(kh + 387);
    const auto *kh_388 = buffer.data(kh + 388);
    const auto *kh_389 = buffer.data(kh + 389);
    const auto *kh_390 = buffer.data(kh + 390);
    const auto *kh_391 = buffer.data(kh + 391);
    const auto *kh_392 = buffer.data(kh + 392);
    const auto *kh_393 = buffer.data(kh + 393);
    const auto *kh_394 = buffer.data(kh + 394);
    const auto *kh_395 = buffer.data(kh + 395);
    const auto *kh_396 = buffer.data(kh + 396);
    const auto *kh_397 = buffer.data(kh + 397);
    const auto *kh_398 = buffer.data(kh + 398);
    const auto *kh_399 = buffer.data(kh + 399);
    const auto *kh_400 = buffer.data(kh + 400);
    const auto *kh_401 = buffer.data(kh + 401);
    const auto *kh_402 = buffer.data(kh + 402);
    const auto *kh_403 = buffer.data(kh + 403);
    const auto *kh_404 = buffer.data(kh + 404);
    const auto *kh_405 = buffer.data(kh + 405);
    const auto *kh_406 = buffer.data(kh + 406);
    const auto *kh_407 = buffer.data(kh + 407);
    const auto *kh_408 = buffer.data(kh + 408);
    const auto *kh_409 = buffer.data(kh + 409);
    const auto *kh_410 = buffer.data(kh + 410);
    const auto *kh_411 = buffer.data(kh + 411);
    const auto *kh_412 = buffer.data(kh + 412);
    const auto *kh_413 = buffer.data(kh + 413);
    const auto *kh_414 = buffer.data(kh + 414);
    const auto *kh_415 = buffer.data(kh + 415);
    const auto *kh_416 = buffer.data(kh + 416);
    const auto *kh_417 = buffer.data(kh + 417);
    const auto *kh_418 = buffer.data(kh + 418);
    const auto *kh_419 = buffer.data(kh + 419);
    const auto *kh_420 = buffer.data(kh + 420);
    const auto *kh_421 = buffer.data(kh + 421);
    const auto *kh_422 = buffer.data(kh + 422);
    const auto *kh_423 = buffer.data(kh + 423);
    const auto *kh_424 = buffer.data(kh + 424);
    const auto *kh_425 = buffer.data(kh + 425);
    const auto *kh_426 = buffer.data(kh + 426);
    const auto *kh_427 = buffer.data(kh + 427);
    const auto *kh_428 = buffer.data(kh + 428);
    const auto *kh_429 = buffer.data(kh + 429);
    const auto *kh_430 = buffer.data(kh + 430);
    const auto *kh_431 = buffer.data(kh + 431);
    const auto *kh_432 = buffer.data(kh + 432);
    const auto *kh_433 = buffer.data(kh + 433);
    const auto *kh_434 = buffer.data(kh + 434);
    const auto *kh_435 = buffer.data(kh + 435);
    const auto *kh_436 = buffer.data(kh + 436);
    const auto *kh_437 = buffer.data(kh + 437);
    const auto *kh_438 = buffer.data(kh + 438);
    const auto *kh_439 = buffer.data(kh + 439);
    const auto *kh_440 = buffer.data(kh + 440);
    const auto *kh_441 = buffer.data(kh + 441);
    const auto *kh_442 = buffer.data(kh + 442);
    const auto *kh_443 = buffer.data(kh + 443);
    const auto *kh_444 = buffer.data(kh + 444);
    const auto *kh_445 = buffer.data(kh + 445);
    const auto *kh_446 = buffer.data(kh + 446);
    const auto *kh_447 = buffer.data(kh + 447);
    const auto *kh_448 = buffer.data(kh + 448);
    const auto *kh_449 = buffer.data(kh + 449);
    const auto *kh_450 = buffer.data(kh + 450);
    const auto *kh_451 = buffer.data(kh + 451);
    const auto *kh_452 = buffer.data(kh + 452);
    const auto *kh_453 = buffer.data(kh + 453);
    const auto *kh_454 = buffer.data(kh + 454);
    const auto *kh_455 = buffer.data(kh + 455);
    const auto *kh_456 = buffer.data(kh + 456);
    const auto *kh_457 = buffer.data(kh + 457);
    const auto *kh_458 = buffer.data(kh + 458);
    const auto *kh_459 = buffer.data(kh + 459);
    const auto *kh_460 = buffer.data(kh + 460);
    const auto *kh_461 = buffer.data(kh + 461);
    const auto *kh_462 = buffer.data(kh + 462);
    const auto *kh_463 = buffer.data(kh + 463);
    const auto *kh_464 = buffer.data(kh + 464);
    const auto *kh_465 = buffer.data(kh + 465);
    const auto *kh_466 = buffer.data(kh + 466);
    const auto *kh_467 = buffer.data(kh + 467);
    const auto *kh_468 = buffer.data(kh + 468);
    const auto *kh_469 = buffer.data(kh + 469);
    const auto *kh_470 = buffer.data(kh + 470);
    const auto *kh_471 = buffer.data(kh + 471);
    const auto *kh_472 = buffer.data(kh + 472);
    const auto *kh_473 = buffer.data(kh + 473);
    const auto *kh_474 = buffer.data(kh + 474);
    const auto *kh_475 = buffer.data(kh + 475);
    const auto *kh_476 = buffer.data(kh + 476);
    const auto *kh_477 = buffer.data(kh + 477);
    const auto *kh_478 = buffer.data(kh + 478);
    const auto *kh_479 = buffer.data(kh + 479);
    const auto *kh_480 = buffer.data(kh + 480);
    const auto *kh_481 = buffer.data(kh + 481);
    const auto *kh_482 = buffer.data(kh + 482);
    const auto *kh_483 = buffer.data(kh + 483);
    const auto *kh_484 = buffer.data(kh + 484);
    const auto *kh_485 = buffer.data(kh + 485);
    const auto *kh_486 = buffer.data(kh + 486);
    const auto *kh_487 = buffer.data(kh + 487);
    const auto *kh_488 = buffer.data(kh + 488);
    const auto *kh_489 = buffer.data(kh + 489);
    const auto *kh_490 = buffer.data(kh + 490);
    const auto *kh_491 = buffer.data(kh + 491);
    const auto *kh_492 = buffer.data(kh + 492);
    const auto *kh_493 = buffer.data(kh + 493);
    const auto *kh_494 = buffer.data(kh + 494);
    const auto *kh_495 = buffer.data(kh + 495);
    const auto *kh_496 = buffer.data(kh + 496);
    const auto *kh_497 = buffer.data(kh + 497);
    const auto *kh_498 = buffer.data(kh + 498);
    const auto *kh_499 = buffer.data(kh + 499);
    const auto *kh_500 = buffer.data(kh + 500);
    const auto *kh_501 = buffer.data(kh + 501);
    const auto *kh_502 = buffer.data(kh + 502);
    const auto *kh_503 = buffer.data(kh + 503);
    const auto *kh_504 = buffer.data(kh + 504);
    const auto *kh_505 = buffer.data(kh + 505);
    const auto *kh_506 = buffer.data(kh + 506);
    const auto *kh_507 = buffer.data(kh + 507);
    const auto *kh_508 = buffer.data(kh + 508);
    const auto *kh_509 = buffer.data(kh + 509);
    const auto *kh_510 = buffer.data(kh + 510);
    const auto *kh_511 = buffer.data(kh + 511);
    const auto *kh_512 = buffer.data(kh + 512);
    const auto *kh_513 = buffer.data(kh + 513);
    const auto *kh_514 = buffer.data(kh + 514);
    const auto *kh_515 = buffer.data(kh + 515);
    const auto *kh_516 = buffer.data(kh + 516);
    const auto *kh_517 = buffer.data(kh + 517);
    const auto *kh_518 = buffer.data(kh + 518);
    const auto *kh_519 = buffer.data(kh + 519);
    const auto *kh_520 = buffer.data(kh + 520);
    const auto *kh_521 = buffer.data(kh + 521);
    const auto *kh_522 = buffer.data(kh + 522);
    const auto *kh_523 = buffer.data(kh + 523);
    const auto *kh_524 = buffer.data(kh + 524);
    const auto *kh_525 = buffer.data(kh + 525);
    const auto *kh_526 = buffer.data(kh + 526);
    const auto *kh_527 = buffer.data(kh + 527);
    const auto *kh_528 = buffer.data(kh + 528);
    const auto *kh_529 = buffer.data(kh + 529);
    const auto *kh_530 = buffer.data(kh + 530);
    const auto *kh_531 = buffer.data(kh + 531);
    const auto *kh_532 = buffer.data(kh + 532);
    const auto *kh_533 = buffer.data(kh + 533);
    const auto *kh_534 = buffer.data(kh + 534);
    const auto *kh_535 = buffer.data(kh + 535);
    const auto *kh_536 = buffer.data(kh + 536);
    const auto *kh_537 = buffer.data(kh + 537);
    const auto *kh_538 = buffer.data(kh + 538);
    const auto *kh_539 = buffer.data(kh + 539);
    const auto *kh_540 = buffer.data(kh + 540);
    const auto *kh_541 = buffer.data(kh + 541);
    const auto *kh_542 = buffer.data(kh + 542);
    const auto *kh_543 = buffer.data(kh + 543);
    const auto *kh_544 = buffer.data(kh + 544);
    const auto *kh_545 = buffer.data(kh + 545);
    const auto *kh_546 = buffer.data(kh + 546);
    const auto *kh_547 = buffer.data(kh + 547);
    const auto *kh_548 = buffer.data(kh + 548);
    const auto *kh_549 = buffer.data(kh + 549);
    const auto *kh_550 = buffer.data(kh + 550);
    const auto *kh_551 = buffer.data(kh + 551);
    const auto *kh_552 = buffer.data(kh + 552);
    const auto *kh_553 = buffer.data(kh + 553);
    const auto *kh_554 = buffer.data(kh + 554);
    const auto *kh_555 = buffer.data(kh + 555);
    const auto *kh_556 = buffer.data(kh + 556);
    const auto *kh_557 = buffer.data(kh + 557);
    const auto *kh_558 = buffer.data(kh + 558);
    const auto *kh_559 = buffer.data(kh + 559);
    const auto *kh_560 = buffer.data(kh + 560);
    const auto *kh_561 = buffer.data(kh + 561);
    const auto *kh_562 = buffer.data(kh + 562);
    const auto *kh_563 = buffer.data(kh + 563);
    const auto *kh_564 = buffer.data(kh + 564);
    const auto *kh_565 = buffer.data(kh + 565);
    const auto *kh_566 = buffer.data(kh + 566);
    const auto *kh_567 = buffer.data(kh + 567);
    const auto *kh_568 = buffer.data(kh + 568);
    const auto *kh_569 = buffer.data(kh + 569);
    const auto *kh_570 = buffer.data(kh + 570);
    const auto *kh_571 = buffer.data(kh + 571);
    const auto *kh_572 = buffer.data(kh + 572);
    const auto *kh_573 = buffer.data(kh + 573);
    const auto *kh_574 = buffer.data(kh + 574);
    const auto *kh_575 = buffer.data(kh + 575);
    const auto *kh_576 = buffer.data(kh + 576);
    const auto *kh_577 = buffer.data(kh + 577);
    const auto *kh_578 = buffer.data(kh + 578);
    const auto *kh_579 = buffer.data(kh + 579);
    const auto *kh_580 = buffer.data(kh + 580);
    const auto *kh_581 = buffer.data(kh + 581);
    const auto *kh_582 = buffer.data(kh + 582);
    const auto *kh_583 = buffer.data(kh + 583);
    const auto *kh_584 = buffer.data(kh + 584);
    const auto *kh_585 = buffer.data(kh + 585);
    const auto *kh_586 = buffer.data(kh + 586);
    const auto *kh_587 = buffer.data(kh + 587);
    const auto *kh_588 = buffer.data(kh + 588);
    const auto *kh_589 = buffer.data(kh + 589);
    const auto *kh_590 = buffer.data(kh + 590);
    const auto *kh_591 = buffer.data(kh + 591);
    const auto *kh_592 = buffer.data(kh + 592);
    const auto *kh_593 = buffer.data(kh + 593);
    const auto *kh_594 = buffer.data(kh + 594);
    const auto *kh_595 = buffer.data(kh + 595);
    const auto *kh_596 = buffer.data(kh + 596);
    const auto *kh_597 = buffer.data(kh + 597);
    const auto *kh_598 = buffer.data(kh + 598);
    const auto *kh_599 = buffer.data(kh + 599);
    const auto *kh_600 = buffer.data(kh + 600);
    const auto *kh_601 = buffer.data(kh + 601);
    const auto *kh_602 = buffer.data(kh + 602);
    const auto *kh_603 = buffer.data(kh + 603);
    const auto *kh_604 = buffer.data(kh + 604);
    const auto *kh_605 = buffer.data(kh + 605);
    const auto *kh_606 = buffer.data(kh + 606);
    const auto *kh_607 = buffer.data(kh + 607);
    const auto *kh_608 = buffer.data(kh + 608);
    const auto *kh_609 = buffer.data(kh + 609);
    const auto *kh_610 = buffer.data(kh + 610);
    const auto *kh_611 = buffer.data(kh + 611);
    const auto *kh_612 = buffer.data(kh + 612);
    const auto *kh_613 = buffer.data(kh + 613);
    const auto *kh_614 = buffer.data(kh + 614);
    const auto *kh_615 = buffer.data(kh + 615);
    const auto *kh_616 = buffer.data(kh + 616);
    const auto *kh_617 = buffer.data(kh + 617);
    const auto *kh_618 = buffer.data(kh + 618);
    const auto *kh_619 = buffer.data(kh + 619);
    const auto *kh_620 = buffer.data(kh + 620);
    const auto *kh_621 = buffer.data(kh + 621);
    const auto *kh_622 = buffer.data(kh + 622);
    const auto *kh_623 = buffer.data(kh + 623);
    const auto *kh_624 = buffer.data(kh + 624);
    const auto *kh_625 = buffer.data(kh + 625);
    const auto *kh_626 = buffer.data(kh + 626);
    const auto *kh_627 = buffer.data(kh + 627);
    const auto *kh_628 = buffer.data(kh + 628);
    const auto *kh_629 = buffer.data(kh + 629);
    const auto *kh_630 = buffer.data(kh + 630);
    const auto *kh_631 = buffer.data(kh + 631);
    const auto *kh_632 = buffer.data(kh + 632);
    const auto *kh_633 = buffer.data(kh + 633);
    const auto *kh_634 = buffer.data(kh + 634);
    const auto *kh_635 = buffer.data(kh + 635);
    const auto *kh_636 = buffer.data(kh + 636);
    const auto *kh_637 = buffer.data(kh + 637);
    const auto *kh_638 = buffer.data(kh + 638);
    const auto *kh_639 = buffer.data(kh + 639);
    const auto *kh_640 = buffer.data(kh + 640);
    const auto *kh_641 = buffer.data(kh + 641);
    const auto *kh_642 = buffer.data(kh + 642);
    const auto *kh_643 = buffer.data(kh + 643);
    const auto *kh_644 = buffer.data(kh + 644);
    const auto *kh_645 = buffer.data(kh + 645);
    const auto *kh_646 = buffer.data(kh + 646);
    const auto *kh_647 = buffer.data(kh + 647);
    const auto *kh_648 = buffer.data(kh + 648);
    const auto *kh_649 = buffer.data(kh + 649);
    const auto *kh_650 = buffer.data(kh + 650);
    const auto *kh_651 = buffer.data(kh + 651);
    const auto *kh_652 = buffer.data(kh + 652);
    const auto *kh_653 = buffer.data(kh + 653);
    const auto *kh_654 = buffer.data(kh + 654);
    const auto *kh_655 = buffer.data(kh + 655);
    const auto *kh_656 = buffer.data(kh + 656);
    const auto *kh_657 = buffer.data(kh + 657);
    const auto *kh_658 = buffer.data(kh + 658);
    const auto *kh_659 = buffer.data(kh + 659);
    const auto *kh_660 = buffer.data(kh + 660);
    const auto *kh_661 = buffer.data(kh + 661);
    const auto *kh_662 = buffer.data(kh + 662);
    const auto *kh_663 = buffer.data(kh + 663);
    const auto *kh_664 = buffer.data(kh + 664);
    const auto *kh_665 = buffer.data(kh + 665);
    const auto *kh_666 = buffer.data(kh + 666);
    const auto *kh_667 = buffer.data(kh + 667);
    const auto *kh_668 = buffer.data(kh + 668);
    const auto *kh_669 = buffer.data(kh + 669);
    const auto *kh_670 = buffer.data(kh + 670);
    const auto *kh_671 = buffer.data(kh + 671);
    const auto *kh_672 = buffer.data(kh + 672);
    const auto *kh_673 = buffer.data(kh + 673);
    const auto *kh_674 = buffer.data(kh + 674);
    const auto *kh_675 = buffer.data(kh + 675);
    const auto *kh_676 = buffer.data(kh + 676);
    const auto *kh_677 = buffer.data(kh + 677);
    const auto *kh_678 = buffer.data(kh + 678);
    const auto *kh_679 = buffer.data(kh + 679);
    const auto *kh_680 = buffer.data(kh + 680);
    const auto *kh_681 = buffer.data(kh + 681);
    const auto *kh_682 = buffer.data(kh + 682);
    const auto *kh_683 = buffer.data(kh + 683);
    const auto *kh_684 = buffer.data(kh + 684);
    const auto *kh_685 = buffer.data(kh + 685);
    const auto *kh_686 = buffer.data(kh + 686);
    const auto *kh_687 = buffer.data(kh + 687);
    const auto *kh_688 = buffer.data(kh + 688);
    const auto *kh_689 = buffer.data(kh + 689);
    const auto *kh_690 = buffer.data(kh + 690);
    const auto *kh_691 = buffer.data(kh + 691);
    const auto *kh_692 = buffer.data(kh + 692);
    const auto *kh_693 = buffer.data(kh + 693);
    const auto *kh_694 = buffer.data(kh + 694);
    const auto *kh_695 = buffer.data(kh + 695);
    const auto *kh_696 = buffer.data(kh + 696);
    const auto *kh_697 = buffer.data(kh + 697);
    const auto *kh_698 = buffer.data(kh + 698);
    const auto *kh_699 = buffer.data(kh + 699);
    const auto *kh_700 = buffer.data(kh + 700);
    const auto *kh_701 = buffer.data(kh + 701);
    const auto *kh_702 = buffer.data(kh + 702);
    const auto *kh_703 = buffer.data(kh + 703);
    const auto *kh_704 = buffer.data(kh + 704);
    const auto *kh_705 = buffer.data(kh + 705);
    const auto *kh_706 = buffer.data(kh + 706);
    const auto *kh_707 = buffer.data(kh + 707);
    const auto *kh_708 = buffer.data(kh + 708);
    const auto *kh_709 = buffer.data(kh + 709);
    const auto *kh_710 = buffer.data(kh + 710);
    const auto *kh_711 = buffer.data(kh + 711);
    const auto *kh_712 = buffer.data(kh + 712);
    const auto *kh_713 = buffer.data(kh + 713);
    const auto *kh_714 = buffer.data(kh + 714);
    const auto *kh_715 = buffer.data(kh + 715);
    const auto *kh_716 = buffer.data(kh + 716);
    const auto *kh_717 = buffer.data(kh + 717);
    const auto *kh_718 = buffer.data(kh + 718);
    const auto *kh_719 = buffer.data(kh + 719);
    const auto *kh_720 = buffer.data(kh + 720);
    const auto *kh_721 = buffer.data(kh + 721);
    const auto *kh_722 = buffer.data(kh + 722);
    const auto *kh_723 = buffer.data(kh + 723);
    const auto *kh_724 = buffer.data(kh + 724);
    const auto *kh_725 = buffer.data(kh + 725);
    const auto *kh_726 = buffer.data(kh + 726);
    const auto *kh_727 = buffer.data(kh + 727);
    const auto *kh_728 = buffer.data(kh + 728);
    const auto *kh_729 = buffer.data(kh + 729);
    const auto *kh_730 = buffer.data(kh + 730);
    const auto *kh_731 = buffer.data(kh + 731);
    const auto *kh_732 = buffer.data(kh + 732);
    const auto *kh_733 = buffer.data(kh + 733);
    const auto *kh_734 = buffer.data(kh + 734);
    const auto *kh_735 = buffer.data(kh + 735);
    const auto *kh_736 = buffer.data(kh + 736);
    const auto *kh_737 = buffer.data(kh + 737);
    const auto *kh_738 = buffer.data(kh + 738);
    const auto *kh_739 = buffer.data(kh + 739);
    const auto *kh_740 = buffer.data(kh + 740);
    const auto *kh_741 = buffer.data(kh + 741);
    const auto *kh_742 = buffer.data(kh + 742);
    const auto *kh_743 = buffer.data(kh + 743);
    const auto *kh_744 = buffer.data(kh + 744);
    const auto *kh_745 = buffer.data(kh + 745);
    const auto *kh_746 = buffer.data(kh + 746);
    const auto *kh_747 = buffer.data(kh + 747);
    const auto *kh_748 = buffer.data(kh + 748);
    const auto *kh_749 = buffer.data(kh + 749);
    const auto *kh_750 = buffer.data(kh + 750);
    const auto *kh_751 = buffer.data(kh + 751);
    const auto *kh_752 = buffer.data(kh + 752);
    const auto *kh_753 = buffer.data(kh + 753);
    const auto *kh_754 = buffer.data(kh + 754);
    const auto *kh_755 = buffer.data(kh + 755);

#pragma omp simd aligned(kh_22, kh_27, kh_36, kh_127, kh_132, kh_141, kh_316, kh_321, kh_330, \
                         kh_589, kh_594, kh_603 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * kh_22[k]
                 - f_1 * kh_27[k]
                 + f_2 * kh_36[k]
                 - f_3 * kh_127[k]
                 + f_4 * kh_132[k]
                 - f_0 * kh_141[k]
                 + f_5 * kh_316[k]
                 - f_6 * kh_321[k]
                 + f_7 * kh_330[k]
                 - f_8 * kh_589[k]
                 + f_9 * kh_594[k]
                 - f_10 * kh_603[k];
    }

#pragma omp simd aligned(kh_25, kh_32, kh_130, kh_137, kh_319, kh_326, kh_592, \
                         kh_599 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_1[k] = f_11 * kh_25[k]
                 - f_11 * kh_32[k]
                 - f_12 * kh_130[k]
                 + f_12 * kh_137[k]
                 + f_13 * kh_319[k]
                 - f_13 * kh_326[k]
                 - f_14 * kh_592[k]
                 + f_14 * kh_599[k];
    }

#pragma omp simd aligned(kh_22, kh_27, kh_29, kh_36, kh_38, kh_127, kh_132, kh_134, kh_141, \
                         kh_143, kh_316, kh_321, kh_323, kh_330, kh_332, kh_589, kh_594, \
                         kh_596, kh_603, kh_605 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_2[k] = -f_15 * kh_22[k]
                 - f_16 * kh_27[k]
                 + f_17 * kh_29[k]
                 + f_18 * kh_36[k]
                 - f_19 * kh_38[k]
                 + f_20 * kh_127[k]
                 + f_21 * kh_132[k]
                 - f_22 * kh_134[k]
                 - f_23 * kh_141[k]
                 + f_24 * kh_143[k]
                 - f_25 * kh_316[k]
                 - f_26 * kh_321[k]
                 + f_27 * kh_323[k]
                 + f_15 * kh_330[k]
                 - f_17 * kh_332[k]
                 + f_28 * kh_589[k]
                 + f_29 * kh_594[k]
                 - f_30 * kh_596[k]
                 - f_31 * kh_603[k]
                 + f_32 * kh_605[k];
    }

#pragma omp simd aligned(kh_25, kh_32, kh_34, kh_130, kh_137, kh_139, kh_319, kh_326, kh_328, \
                         kh_592, kh_599, kh_601 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = -f_33 * kh_25[k]
                 - f_33 * kh_32[k]
                 + f_34 * kh_34[k]
                 + f_35 * kh_130[k]
                 + f_35 * kh_137[k]
                 - f_36 * kh_139[k]
                 - f_37 * kh_319[k]
                 - f_37 * kh_326[k]
                 + f_38 * kh_328[k]
                 + f_39 * kh_592[k]
                 + f_39 * kh_599[k]
                 - f_40 * kh_601[k];
    }

#pragma omp simd aligned(kh_22, kh_27, kh_29, kh_36, kh_38, kh_40, kh_127, kh_132, kh_134, \
                         kh_141, kh_143, kh_145, kh_316, kh_321, kh_323, kh_330, kh_332, \
                         kh_334, kh_589, kh_594, kh_596, kh_603, kh_605, \
                         kh_607 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = f_41 * kh_22[k]
                 + f_42 * kh_27[k]
                 - f_43 * kh_29[k]
                 + f_41 * kh_36[k]
                 - f_43 * kh_38[k]
                 + f_44 * kh_40[k]
                 - f_45 * kh_127[k]
                 - f_46 * kh_132[k]
                 + f_47 * kh_134[k]
                 - f_45 * kh_141[k]
                 + f_47 * kh_143[k]
                 - f_48 * kh_145[k]
                 + f_49 * kh_316[k]
                 + f_50 * kh_321[k]
                 - f_51 * kh_323[k]
                 + f_49 * kh_330[k]
                 - f_51 * kh_332[k]
                 + f_52 * kh_334[k]
                 - f_53 * kh_589[k]
                 - f_54 * kh_594[k]
                 + f_55 * kh_596[k]
                 - f_53 * kh_603[k]
                 + f_55 * kh_605[k]
                 - f_56 * kh_607[k];
    }

#pragma omp simd aligned(kh_23, kh_28, kh_30, kh_37, kh_39, kh_41, kh_128, kh_133, kh_135, \
                         kh_142, kh_144, kh_146, kh_317, kh_322, kh_324, kh_331, kh_333, \
                         kh_335, kh_590, kh_595, kh_597, kh_604, kh_606, \
                         kh_608 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = f_57 * kh_23[k]
                 + f_58 * kh_28[k]
                 - f_59 * kh_30[k]
                 + f_57 * kh_37[k]
                 - f_59 * kh_39[k]
                 + f_60 * kh_41[k]
                 - f_61 * kh_128[k]
                 - f_62 * kh_133[k]
                 + f_63 * kh_135[k]
                 - f_61 * kh_142[k]
                 + f_63 * kh_144[k]
                 - f_59 * kh_146[k]
                 + f_64 * kh_317[k]
                 + f_65 * kh_322[k]
                 - f_66 * kh_324[k]
                 + f_64 * kh_331[k]
                 - f_66 * kh_333[k]
                 + f_67 * kh_335[k]
                 - f_68 * kh_590[k]
                 - f_69 * kh_595[k]
                 + f_70 * kh_597[k]
                 - f_68 * kh_604[k]
                 + f_70 * kh_606[k]
                 - f_71 * kh_608[k];
    }

#pragma omp simd aligned(kh_21, kh_24, kh_26, kh_31, kh_33, kh_35, kh_126, kh_129, kh_131, \
                         kh_136, kh_138, kh_140, kh_315, kh_318, kh_320, kh_325, kh_327, \
                         kh_329, kh_588, kh_591, kh_593, kh_598, kh_600, \
                         kh_602 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = f_41 * kh_21[k]
                 + f_42 * kh_24[k]
                 - f_43 * kh_26[k]
                 + f_41 * kh_31[k]
                 - f_43 * kh_33[k]
                 + f_44 * kh_35[k]
                 - f_45 * kh_126[k]
                 - f_46 * kh_129[k]
                 + f_47 * kh_131[k]
                 - f_45 * kh_136[k]
                 + f_47 * kh_138[k]
                 - f_48 * kh_140[k]
                 + f_49 * kh_315[k]
                 + f_50 * kh_318[k]
                 - f_51 * kh_320[k]
                 + f_49 * kh_325[k]
                 - f_51 * kh_327[k]
                 + f_52 * kh_329[k]
                 - f_53 * kh_588[k]
                 - f_54 * kh_591[k]
                 + f_55 * kh_593[k]
                 - f_53 * kh_598[k]
                 + f_55 * kh_600[k]
                 - f_56 * kh_602[k];
    }

#pragma omp simd aligned(kh_23, kh_30, kh_37, kh_39, kh_128, kh_135, kh_142, kh_144, kh_317, \
                         kh_324, kh_331, kh_333, kh_590, kh_597, kh_604, \
                         kh_606 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = -f_72 * kh_23[k]
                 + f_33 * kh_30[k]
                 + f_72 * kh_37[k]
                 - f_33 * kh_39[k]
                 + f_73 * kh_128[k]
                 - f_35 * kh_135[k]
                 - f_73 * kh_142[k]
                 + f_35 * kh_144[k]
                 - f_74 * kh_317[k]
                 + f_37 * kh_324[k]
                 + f_74 * kh_331[k]
                 - f_37 * kh_333[k]
                 + f_75 * kh_590[k]
                 - f_39 * kh_597[k]
                 - f_75 * kh_604[k]
                 + f_39 * kh_606[k];
    }

#pragma omp simd aligned(kh_21, kh_24, kh_26, kh_31, kh_33, kh_126, kh_129, kh_131, kh_136, \
                         kh_138, kh_315, kh_318, kh_320, kh_325, kh_327, kh_588, kh_591, \
                         kh_593, kh_598, kh_600 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_8[k] = -f_18 * kh_21[k]
                 + f_16 * kh_24[k]
                 + f_19 * kh_26[k]
                 + f_15 * kh_31[k]
                 - f_17 * kh_33[k]
                 + f_23 * kh_126[k]
                 - f_21 * kh_129[k]
                 - f_24 * kh_131[k]
                 - f_20 * kh_136[k]
                 + f_22 * kh_138[k]
                 - f_15 * kh_315[k]
                 + f_26 * kh_318[k]
                 + f_17 * kh_320[k]
                 + f_25 * kh_325[k]
                 - f_27 * kh_327[k]
                 + f_31 * kh_588[k]
                 - f_29 * kh_591[k]
                 - f_32 * kh_593[k]
                 - f_28 * kh_598[k]
                 + f_30 * kh_600[k];
    }

#pragma omp simd aligned(kh_23, kh_28, kh_37, kh_128, kh_133, kh_142, kh_317, kh_322, kh_331, \
                         kh_590, kh_595, kh_604 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_9[k] = f_76 * kh_23[k]
                 - f_77 * kh_28[k]
                 + f_76 * kh_37[k]
                 - f_78 * kh_128[k]
                 + f_79 * kh_133[k]
                 - f_78 * kh_142[k]
                 + f_80 * kh_317[k]
                 - f_81 * kh_322[k]
                 + f_80 * kh_331[k]
                 - f_82 * kh_590[k]
                 + f_83 * kh_595[k]
                 - f_82 * kh_604[k];
    }

#pragma omp simd aligned(kh_21, kh_24, kh_31, kh_126, kh_129, kh_136, kh_315, kh_318, kh_325, \
                         kh_588, kh_591, kh_598 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = f_2 * kh_21[k]
                  - f_1 * kh_24[k]
                  + f_0 * kh_31[k]
                  - f_0 * kh_126[k]
                  + f_4 * kh_129[k]
                  - f_3 * kh_136[k]
                  + f_7 * kh_315[k]
                  - f_6 * kh_318[k]
                  + f_5 * kh_325[k]
                  - f_10 * kh_588[k]
                  + f_9 * kh_591[k]
                  - f_8 * kh_598[k];
    }

#pragma omp simd aligned(kh_85, kh_90, kh_99, kh_232, kh_237, kh_246, kh_463, kh_468, \
                         kh_477 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_11[k] = f_65 * kh_85[k]
                  - f_84 * kh_90[k]
                  + f_85 * kh_99[k]
                  - f_86 * kh_232[k]
                  + f_87 * kh_237[k]
                  - f_88 * kh_246[k]
                  + f_65 * kh_463[k]
                  - f_84 * kh_468[k]
                  + f_85 * kh_477[k];
    }

#pragma omp simd aligned(kh_88, kh_95, kh_235, kh_242, kh_466, kh_473 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_12[k] = f_89 * kh_88[k]
                  - f_89 * kh_95[k]
                  - f_90 * kh_235[k]
                  + f_90 * kh_242[k]
                  + f_89 * kh_466[k]
                  - f_89 * kh_473[k];
    }

#pragma omp simd aligned(kh_85, kh_90, kh_92, kh_99, kh_101, kh_232, kh_237, kh_239, kh_246, \
                         kh_248, kh_463, kh_468, kh_470, kh_477, \
                         kh_479 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = -f_91 * kh_85[k]
                  - f_92 * kh_90[k]
                  + f_93 * kh_92[k]
                  + f_94 * kh_99[k]
                  - f_95 * kh_101[k]
                  + f_96 * kh_232[k]
                  + f_97 * kh_237[k]
                  - f_98 * kh_239[k]
                  - f_99 * kh_246[k]
                  + f_100 * kh_248[k]
                  - f_91 * kh_463[k]
                  - f_92 * kh_468[k]
                  + f_93 * kh_470[k]
                  + f_94 * kh_477[k]
                  - f_95 * kh_479[k];
    }

#pragma omp simd aligned(kh_88, kh_95, kh_97, kh_235, kh_242, kh_244, kh_466, kh_473, \
                         kh_475 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = -f_101 * kh_88[k]
                  - f_101 * kh_95[k]
                  + f_102 * kh_97[k]
                  + f_103 * kh_235[k]
                  + f_103 * kh_242[k]
                  - f_104 * kh_244[k]
                  - f_101 * kh_466[k]
                  - f_101 * kh_473[k]
                  + f_102 * kh_475[k];
    }

#pragma omp simd aligned(kh_85, kh_90, kh_92, kh_99, kh_101, kh_103, kh_232, kh_237, kh_239, \
                         kh_246, kh_248, kh_250, kh_463, kh_468, kh_470, kh_477, kh_479, \
                         kh_481 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_15[k] = f_105 * kh_85[k]
                  + f_106 * kh_90[k]
                  - f_107 * kh_92[k]
                  + f_105 * kh_99[k]
                  - f_107 * kh_101[k]
                  + f_108 * kh_103[k]
                  - f_109 * kh_232[k]
                  - f_110 * kh_237[k]
                  + f_111 * kh_239[k]
                  - f_109 * kh_246[k]
                  + f_111 * kh_248[k]
                  - f_112 * kh_250[k]
                  + f_105 * kh_463[k]
                  + f_106 * kh_468[k]
                  - f_107 * kh_470[k]
                  + f_105 * kh_477[k]
                  - f_107 * kh_479[k]
                  + f_108 * kh_481[k];
    }

#pragma omp simd aligned(kh_86, kh_91, kh_93, kh_100, kh_102, kh_104, kh_233, kh_238, kh_240, \
                         kh_247, kh_249, kh_251, kh_464, kh_469, kh_471, kh_478, kh_480, \
                         kh_482 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_16[k] = f_113 * kh_86[k]
                  + f_114 * kh_91[k]
                  - f_115 * kh_93[k]
                  + f_113 * kh_100[k]
                  - f_115 * kh_102[k]
                  + f_116 * kh_104[k]
                  - f_117 * kh_233[k]
                  - f_118 * kh_238[k]
                  + f_119 * kh_240[k]
                  - f_117 * kh_247[k]
                  + f_119 * kh_249[k]
                  - f_120 * kh_251[k]
                  + f_113 * kh_464[k]
                  + f_114 * kh_469[k]
                  - f_115 * kh_471[k]
                  + f_113 * kh_478[k]
                  - f_115 * kh_480[k]
                  + f_116 * kh_482[k];
    }

#pragma omp simd aligned(kh_84, kh_87, kh_89, kh_94, kh_96, kh_98, kh_231, kh_234, kh_236, \
                         kh_241, kh_243, kh_245, kh_462, kh_465, kh_467, kh_472, kh_474, \
                         kh_476 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = f_105 * kh_84[k]
                  + f_106 * kh_87[k]
                  - f_107 * kh_89[k]
                  + f_105 * kh_94[k]
                  - f_107 * kh_96[k]
                  + f_108 * kh_98[k]
                  - f_109 * kh_231[k]
                  - f_110 * kh_234[k]
                  + f_111 * kh_236[k]
                  - f_109 * kh_241[k]
                  + f_111 * kh_243[k]
                  - f_112 * kh_245[k]
                  + f_105 * kh_462[k]
                  + f_106 * kh_465[k]
                  - f_107 * kh_467[k]
                  + f_105 * kh_472[k]
                  - f_107 * kh_474[k]
                  + f_108 * kh_476[k];
    }

#pragma omp simd aligned(kh_86, kh_93, kh_100, kh_102, kh_233, kh_240, kh_247, kh_249, kh_464, \
                         kh_471, kh_478, kh_480 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = -f_121 * kh_86[k]
                  + f_101 * kh_93[k]
                  + f_121 * kh_100[k]
                  - f_101 * kh_102[k]
                  + f_122 * kh_233[k]
                  - f_103 * kh_240[k]
                  - f_122 * kh_247[k]
                  + f_103 * kh_249[k]
                  - f_121 * kh_464[k]
                  + f_101 * kh_471[k]
                  + f_121 * kh_478[k]
                  - f_101 * kh_480[k];
    }

#pragma omp simd aligned(kh_84, kh_87, kh_89, kh_94, kh_96, kh_231, kh_234, kh_236, kh_241, \
                         kh_243, kh_462, kh_465, kh_467, kh_472, \
                         kh_474 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = -f_94 * kh_84[k]
                  + f_92 * kh_87[k]
                  + f_95 * kh_89[k]
                  + f_91 * kh_94[k]
                  - f_93 * kh_96[k]
                  + f_99 * kh_231[k]
                  - f_97 * kh_234[k]
                  - f_100 * kh_236[k]
                  - f_96 * kh_241[k]
                  + f_98 * kh_243[k]
                  - f_94 * kh_462[k]
                  + f_92 * kh_465[k]
                  + f_95 * kh_467[k]
                  + f_91 * kh_472[k]
                  - f_93 * kh_474[k];
    }

#pragma omp simd aligned(kh_86, kh_91, kh_100, kh_233, kh_238, kh_247, kh_464, kh_469, \
                         kh_478 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = f_123 * kh_86[k]
                  - f_124 * kh_91[k]
                  + f_123 * kh_100[k]
                  - f_125 * kh_233[k]
                  + f_126 * kh_238[k]
                  - f_125 * kh_247[k]
                  + f_123 * kh_464[k]
                  - f_124 * kh_469[k]
                  + f_123 * kh_478[k];
    }

#pragma omp simd aligned(kh_84, kh_87, kh_94, kh_231, kh_234, kh_241, kh_462, kh_465, \
                         kh_472 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = f_85 * kh_84[k]
                  - f_84 * kh_87[k]
                  + f_65 * kh_94[k]
                  - f_88 * kh_231[k]
                  + f_87 * kh_234[k]
                  - f_86 * kh_241[k]
                  + f_85 * kh_462[k]
                  - f_84 * kh_465[k]
                  + f_65 * kh_472[k];
    }

#pragma omp simd aligned(kh_22, kh_27, kh_36, kh_127, kh_132, kh_141, kh_169, kh_174, kh_183, \
                         kh_316, kh_321, kh_330, kh_358, kh_363, kh_372, kh_589, kh_594, \
                         kh_603, kh_631, kh_636, kh_645 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = -f_127 * kh_22[k]
                  + f_128 * kh_27[k]
                  - f_129 * kh_36[k]
                  + f_127 * kh_127[k]
                  - f_128 * kh_132[k]
                  + f_129 * kh_141[k]
                  + f_130 * kh_169[k]
                  - f_131 * kh_174[k]
                  + f_132 * kh_183[k]
                  + f_133 * kh_316[k]
                  - f_134 * kh_321[k]
                  + f_135 * kh_330[k]
                  - f_131 * kh_358[k]
                  + f_136 * kh_363[k]
                  - f_137 * kh_372[k]
                  - f_129 * kh_589[k]
                  + f_138 * kh_594[k]
                  - f_139 * kh_603[k]
                  + f_132 * kh_631[k]
                  - f_137 * kh_636[k]
                  + f_140 * kh_645[k];
    }

#pragma omp simd aligned(kh_25, kh_32, kh_130, kh_137, kh_172, kh_179, kh_319, kh_326, kh_361, \
                         kh_368, kh_592, kh_599, kh_634, kh_641 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = -f_141 * kh_25[k]
                  + f_141 * kh_32[k]
                  + f_141 * kh_130[k]
                  - f_141 * kh_137[k]
                  + f_142 * kh_172[k]
                  - f_142 * kh_179[k]
                  + f_143 * kh_319[k]
                  - f_143 * kh_326[k]
                  - f_144 * kh_361[k]
                  + f_144 * kh_368[k]
                  - f_145 * kh_592[k]
                  + f_145 * kh_599[k]
                  + f_146 * kh_634[k]
                  - f_146 * kh_641[k];
    }

#pragma omp simd aligned(kh_22, kh_27, kh_29, kh_36, kh_38, kh_127, kh_132, kh_134, kh_141, \
                         kh_143, kh_169, kh_174, kh_176, kh_183, kh_185, kh_316, kh_321, \
                         kh_323, kh_330, kh_332, kh_358, kh_363, kh_365, kh_372, kh_374, \
                         kh_589, kh_594, kh_596, kh_603, kh_605, kh_631, kh_636, kh_638, \
                         kh_645, kh_647 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = f_147 * kh_22[k]
                  + f_148 * kh_27[k]
                  - f_149 * kh_29[k]
                  - f_150 * kh_36[k]
                  + f_151 * kh_38[k]
                  - f_147 * kh_127[k]
                  - f_148 * kh_132[k]
                  + f_149 * kh_134[k]
                  + f_150 * kh_141[k]
                  - f_151 * kh_143[k]
                  - f_152 * kh_169[k]
                  - f_149 * kh_174[k]
                  + f_153 * kh_176[k]
                  + f_154 * kh_183[k]
                  - f_155 * kh_185[k]
                  - f_156 * kh_316[k]
                  - f_157 * kh_321[k]
                  + f_158 * kh_323[k]
                  + f_159 * kh_330[k]
                  - f_160 * kh_332[k]
                  + f_161 * kh_358[k]
                  + f_162 * kh_363[k]
                  - f_163 * kh_365[k]
                  - f_149 * kh_372[k]
                  + f_164 * kh_374[k]
                  + f_165 * kh_589[k]
                  + f_166 * kh_594[k]
                  - f_167 * kh_596[k]
                  - f_168 * kh_603[k]
                  + f_169 * kh_605[k]
                  - f_170 * kh_631[k]
                  - f_167 * kh_636[k]
                  + f_171 * kh_638[k]
                  + f_172 * kh_645[k]
                  - f_173 * kh_647[k];
    }

#pragma omp simd aligned(kh_25, kh_32, kh_34, kh_130, kh_137, kh_139, kh_172, kh_179, kh_181, \
                         kh_319, kh_326, kh_328, kh_361, kh_368, kh_370, kh_592, kh_599, \
                         kh_601, kh_634, kh_641, kh_643 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = f_174 * kh_25[k]
                  + f_174 * kh_32[k]
                  - f_175 * kh_34[k]
                  - f_174 * kh_130[k]
                  - f_174 * kh_137[k]
                  + f_175 * kh_139[k]
                  - f_176 * kh_172[k]
                  - f_176 * kh_179[k]
                  + f_177 * kh_181[k]
                  - f_178 * kh_319[k]
                  - f_178 * kh_326[k]
                  + f_179 * kh_328[k]
                  + f_177 * kh_361[k]
                  + f_177 * kh_368[k]
                  - f_180 * kh_370[k]
                  + f_181 * kh_592[k]
                  + f_181 * kh_599[k]
                  - f_182 * kh_601[k]
                  - f_183 * kh_634[k]
                  - f_183 * kh_641[k]
                  + f_184 * kh_643[k];
    }

#pragma omp simd aligned(kh_22, kh_27, kh_29, kh_36, kh_38, kh_40, kh_127, kh_132, kh_134, \
                         kh_141, kh_143, kh_145, kh_169, kh_174, kh_176, kh_183, kh_185, \
                         kh_187, kh_316, kh_321, kh_323, kh_330, kh_332, kh_334, kh_358, \
                         kh_363, kh_365, kh_372, kh_374, kh_376, kh_589, kh_594, kh_596, \
                         kh_603, kh_605, kh_607, kh_631, kh_636, kh_638, kh_645, kh_647, \
                         kh_649 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = -f_185 * kh_22[k]
                  - f_186 * kh_27[k]
                  + f_187 * kh_29[k]
                  - f_185 * kh_36[k]
                  + f_187 * kh_38[k]
                  - f_188 * kh_40[k]
                  + f_185 * kh_127[k]
                  + f_186 * kh_132[k]
                  - f_187 * kh_134[k]
                  + f_185 * kh_141[k]
                  - f_187 * kh_143[k]
                  + f_188 * kh_145[k]
                  + f_187 * kh_169[k]
                  + f_189 * kh_174[k]
                  - f_190 * kh_176[k]
                  + f_187 * kh_183[k]
                  - f_190 * kh_185[k]
                  + f_191 * kh_187[k]
                  + f_192 * kh_316[k]
                  + f_193 * kh_321[k]
                  - f_194 * kh_323[k]
                  + f_192 * kh_330[k]
                  - f_194 * kh_332[k]
                  + f_195 * kh_334[k]
                  - f_189 * kh_358[k]
                  - f_196 * kh_363[k]
                  + f_197 * kh_365[k]
                  - f_189 * kh_372[k]
                  + f_197 * kh_374[k]
                  - f_198 * kh_376[k]
                  - f_199 * kh_589[k]
                  - f_200 * kh_594[k]
                  + f_201 * kh_596[k]
                  - f_199 * kh_603[k]
                  + f_201 * kh_605[k]
                  - f_202 * kh_607[k]
                  + f_201 * kh_631[k]
                  + f_203 * kh_636[k]
                  - f_204 * kh_638[k]
                  + f_201 * kh_645[k]
                  - f_204 * kh_647[k]
                  + f_205 * kh_649[k];
    }

#pragma omp simd aligned(kh_23, kh_28, kh_30, kh_37, kh_39, kh_41, kh_128, kh_133, kh_135, \
                         kh_142, kh_144, kh_146, kh_170, kh_175, kh_177, kh_184, kh_186, \
                         kh_188, kh_317, kh_322, kh_324, kh_331, kh_333, kh_335, kh_359, \
                         kh_364, kh_366, kh_373, kh_375, kh_377, kh_590, kh_595, kh_597, \
                         kh_604, kh_606, kh_608, kh_632, kh_637, kh_639, kh_646, kh_648, \
                         kh_650 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = -f_206 * kh_23[k]
                  - f_207 * kh_28[k]
                  + f_208 * kh_30[k]
                  - f_206 * kh_37[k]
                  + f_208 * kh_39[k]
                  - f_209 * kh_41[k]
                  + f_206 * kh_128[k]
                  + f_207 * kh_133[k]
                  - f_208 * kh_135[k]
                  + f_206 * kh_142[k]
                  - f_208 * kh_144[k]
                  + f_209 * kh_146[k]
                  + f_210 * kh_170[k]
                  + f_211 * kh_175[k]
                  - f_212 * kh_177[k]
                  + f_210 * kh_184[k]
                  - f_212 * kh_186[k]
                  + f_213 * kh_188[k]
                  + f_214 * kh_317[k]
                  + f_215 * kh_322[k]
                  - f_216 * kh_324[k]
                  + f_214 * kh_331[k]
                  - f_216 * kh_333[k]
                  + f_217 * kh_335[k]
                  - f_211 * kh_359[k]
                  - f_218 * kh_364[k]
                  + f_219 * kh_366[k]
                  - f_211 * kh_373[k]
                  + f_219 * kh_375[k]
                  - f_220 * kh_377[k]
                  - f_221 * kh_590[k]
                  - f_222 * kh_595[k]
                  + f_209 * kh_597[k]
                  - f_221 * kh_604[k]
                  + f_209 * kh_606[k]
                  - f_223 * kh_608[k]
                  + f_224 * kh_632[k]
                  + f_216 * kh_637[k]
                  - f_213 * kh_639[k]
                  + f_224 * kh_646[k]
                  - f_213 * kh_648[k]
                  + f_225 * kh_650[k];
    }

#pragma omp simd aligned(kh_21, kh_24, kh_26, kh_31, kh_33, kh_35, kh_126, kh_129, kh_131, \
                         kh_136, kh_138, kh_140, kh_168, kh_171, kh_173, kh_178, kh_180, \
                         kh_182, kh_315, kh_318, kh_320, kh_325, kh_327, kh_329, kh_357, \
                         kh_360, kh_362, kh_367, kh_369, kh_371, kh_588, kh_591, kh_593, \
                         kh_598, kh_600, kh_602, kh_630, kh_633, kh_635, kh_640, kh_642, \
                         kh_644 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = -f_185 * kh_21[k]
                  - f_186 * kh_24[k]
                  + f_187 * kh_26[k]
                  - f_185 * kh_31[k]
                  + f_187 * kh_33[k]
                  - f_188 * kh_35[k]
                  + f_185 * kh_126[k]
                  + f_186 * kh_129[k]
                  - f_187 * kh_131[k]
                  + f_185 * kh_136[k]
                  - f_187 * kh_138[k]
                  + f_188 * kh_140[k]
                  + f_187 * kh_168[k]
                  + f_189 * kh_171[k]
                  - f_190 * kh_173[k]
                  + f_187 * kh_178[k]
                  - f_190 * kh_180[k]
                  + f_191 * kh_182[k]
                  + f_192 * kh_315[k]
                  + f_193 * kh_318[k]
                  - f_194 * kh_320[k]
                  + f_192 * kh_325[k]
                  - f_194 * kh_327[k]
                  + f_195 * kh_329[k]
                  - f_189 * kh_357[k]
                  - f_196 * kh_360[k]
                  + f_197 * kh_362[k]
                  - f_189 * kh_367[k]
                  + f_197 * kh_369[k]
                  - f_198 * kh_371[k]
                  - f_199 * kh_588[k]
                  - f_200 * kh_591[k]
                  + f_201 * kh_593[k]
                  - f_199 * kh_598[k]
                  + f_201 * kh_600[k]
                  - f_202 * kh_602[k]
                  + f_201 * kh_630[k]
                  + f_203 * kh_633[k]
                  - f_204 * kh_635[k]
                  + f_201 * kh_640[k]
                  - f_204 * kh_642[k]
                  + f_205 * kh_644[k];
    }

#pragma omp simd aligned(kh_23, kh_30, kh_37, kh_39, kh_128, kh_135, kh_142, kh_144, kh_170, \
                         kh_177, kh_184, kh_186, kh_317, kh_324, kh_331, kh_333, kh_359, \
                         kh_366, kh_373, kh_375, kh_590, kh_597, kh_604, kh_606, kh_632, \
                         kh_639, kh_646, kh_648 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = f_226 * kh_23[k]
                  - f_174 * kh_30[k]
                  - f_226 * kh_37[k]
                  + f_174 * kh_39[k]
                  - f_226 * kh_128[k]
                  + f_174 * kh_135[k]
                  + f_226 * kh_142[k]
                  - f_174 * kh_144[k]
                  - f_227 * kh_170[k]
                  + f_176 * kh_177[k]
                  + f_227 * kh_184[k]
                  - f_176 * kh_186[k]
                  - f_228 * kh_317[k]
                  + f_178 * kh_324[k]
                  + f_228 * kh_331[k]
                  - f_178 * kh_333[k]
                  + f_176 * kh_359[k]
                  - f_177 * kh_366[k]
                  - f_176 * kh_373[k]
                  + f_177 * kh_375[k]
                  + f_229 * kh_590[k]
                  - f_181 * kh_597[k]
                  - f_229 * kh_604[k]
                  + f_181 * kh_606[k]
                  - f_230 * kh_632[k]
                  + f_183 * kh_639[k]
                  + f_230 * kh_646[k]
                  - f_183 * kh_648[k];
    }

#pragma omp simd aligned(kh_21, kh_24, kh_26, kh_31, kh_33, kh_126, kh_129, kh_131, kh_136, \
                         kh_138, kh_168, kh_171, kh_173, kh_178, kh_180, kh_315, kh_318, \
                         kh_320, kh_325, kh_327, kh_357, kh_360, kh_362, kh_367, kh_369, \
                         kh_588, kh_591, kh_593, kh_598, kh_600, kh_630, kh_633, kh_635, \
                         kh_640, kh_642 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = f_150 * kh_21[k]
                  - f_148 * kh_24[k]
                  - f_151 * kh_26[k]
                  - f_147 * kh_31[k]
                  + f_149 * kh_33[k]
                  - f_150 * kh_126[k]
                  + f_148 * kh_129[k]
                  + f_151 * kh_131[k]
                  + f_147 * kh_136[k]
                  - f_149 * kh_138[k]
                  - f_154 * kh_168[k]
                  + f_149 * kh_171[k]
                  + f_155 * kh_173[k]
                  + f_152 * kh_178[k]
                  - f_153 * kh_180[k]
                  - f_159 * kh_315[k]
                  + f_157 * kh_318[k]
                  + f_160 * kh_320[k]
                  + f_156 * kh_325[k]
                  - f_158 * kh_327[k]
                  + f_149 * kh_357[k]
                  - f_162 * kh_360[k]
                  - f_164 * kh_362[k]
                  - f_161 * kh_367[k]
                  + f_163 * kh_369[k]
                  + f_168 * kh_588[k]
                  - f_166 * kh_591[k]
                  - f_169 * kh_593[k]
                  - f_165 * kh_598[k]
                  + f_167 * kh_600[k]
                  - f_172 * kh_630[k]
                  + f_167 * kh_633[k]
                  + f_173 * kh_635[k]
                  + f_170 * kh_640[k]
                  - f_171 * kh_642[k];
    }

#pragma omp simd aligned(kh_23, kh_28, kh_37, kh_128, kh_133, kh_142, kh_170, kh_175, kh_184, \
                         kh_317, kh_322, kh_331, kh_359, kh_364, kh_373, kh_590, kh_595, \
                         kh_604, kh_632, kh_637, kh_646 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = -f_231 * kh_23[k]
                  + f_232 * kh_28[k]
                  - f_231 * kh_37[k]
                  + f_231 * kh_128[k]
                  - f_232 * kh_133[k]
                  + f_231 * kh_142[k]
                  + f_233 * kh_170[k]
                  - f_234 * kh_175[k]
                  + f_233 * kh_184[k]
                  + f_235 * kh_317[k]
                  - f_236 * kh_322[k]
                  + f_235 * kh_331[k]
                  - f_237 * kh_359[k]
                  + f_238 * kh_364[k]
                  - f_237 * kh_373[k]
                  - f_239 * kh_590[k]
                  + f_240 * kh_595[k]
                  - f_239 * kh_604[k]
                  + f_241 * kh_632[k]
                  - f_242 * kh_637[k]
                  + f_241 * kh_646[k];
    }

#pragma omp simd aligned(kh_21, kh_24, kh_31, kh_126, kh_129, kh_136, kh_168, kh_171, kh_178, \
                         kh_315, kh_318, kh_325, kh_357, kh_360, kh_367, kh_588, kh_591, \
                         kh_598, kh_630, kh_633, kh_640 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = -f_129 * kh_21[k]
                  + f_128 * kh_24[k]
                  - f_127 * kh_31[k]
                  + f_129 * kh_126[k]
                  - f_128 * kh_129[k]
                  + f_127 * kh_136[k]
                  + f_132 * kh_168[k]
                  - f_131 * kh_171[k]
                  + f_130 * kh_178[k]
                  + f_135 * kh_315[k]
                  - f_134 * kh_318[k]
                  + f_133 * kh_325[k]
                  - f_137 * kh_357[k]
                  + f_136 * kh_360[k]
                  - f_131 * kh_367[k]
                  - f_139 * kh_588[k]
                  + f_138 * kh_591[k]
                  - f_129 * kh_598[k]
                  + f_140 * kh_630[k]
                  - f_137 * kh_633[k]
                  + f_132 * kh_640[k];
    }

#pragma omp simd aligned(kh_85, kh_90, kh_99, kh_274, kh_279, kh_288, kh_463, kh_468, kh_477, \
                         kh_505, kh_510, kh_519 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = -f_137 * kh_85[k]
                  + f_243 * kh_90[k]
                  - f_244 * kh_99[k]
                  + f_245 * kh_274[k]
                  - f_246 * kh_279[k]
                  + f_247 * kh_288[k]
                  + f_137 * kh_463[k]
                  - f_243 * kh_468[k]
                  + f_244 * kh_477[k]
                  - f_245 * kh_505[k]
                  + f_246 * kh_510[k]
                  - f_247 * kh_519[k];
    }

#pragma omp simd aligned(kh_88, kh_95, kh_277, kh_284, kh_466, kh_473, kh_508, \
                         kh_515 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = -f_248 * kh_88[k]
                  + f_248 * kh_95[k]
                  + f_249 * kh_277[k]
                  - f_249 * kh_284[k]
                  + f_248 * kh_466[k]
                  - f_248 * kh_473[k]
                  - f_249 * kh_508[k]
                  + f_249 * kh_515[k];
    }

#pragma omp simd aligned(kh_85, kh_90, kh_92, kh_99, kh_101, kh_274, kh_279, kh_281, kh_288, \
                         kh_290, kh_463, kh_468, kh_470, kh_477, kh_479, kh_505, kh_510, \
                         kh_512, kh_519, kh_521 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = f_160 * kh_85[k]
                  + f_250 * kh_90[k]
                  - f_251 * kh_92[k]
                  - f_167 * kh_99[k]
                  + f_252 * kh_101[k]
                  - f_162 * kh_274[k]
                  - f_253 * kh_279[k]
                  + f_254 * kh_281[k]
                  + f_255 * kh_288[k]
                  - f_256 * kh_290[k]
                  - f_160 * kh_463[k]
                  - f_250 * kh_468[k]
                  + f_251 * kh_470[k]
                  + f_167 * kh_477[k]
                  - f_252 * kh_479[k]
                  + f_162 * kh_505[k]
                  + f_253 * kh_510[k]
                  - f_254 * kh_512[k]
                  - f_255 * kh_519[k]
                  + f_256 * kh_521[k];
    }

#pragma omp simd aligned(kh_88, kh_95, kh_97, kh_277, kh_284, kh_286, kh_466, kh_473, kh_475, \
                         kh_508, kh_515, kh_517 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_36[k] = f_184 * kh_88[k]
                  + f_184 * kh_95[k]
                  - f_257 * kh_97[k]
                  - f_258 * kh_277[k]
                  - f_258 * kh_284[k]
                  + f_259 * kh_286[k]
                  - f_184 * kh_466[k]
                  - f_184 * kh_473[k]
                  + f_257 * kh_475[k]
                  + f_258 * kh_508[k]
                  + f_258 * kh_515[k]
                  - f_259 * kh_517[k];
    }

#pragma omp simd aligned(kh_85, kh_90, kh_92, kh_99, kh_101, kh_103, kh_274, kh_279, kh_281, \
                         kh_288, kh_290, kh_292, kh_463, kh_468, kh_470, kh_477, kh_479, \
                         kh_481, kh_505, kh_510, kh_512, kh_519, kh_521, \
                         kh_523 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_37[k] = -f_203 * kh_85[k]
                  - f_260 * kh_90[k]
                  + f_261 * kh_92[k]
                  - f_203 * kh_99[k]
                  + f_261 * kh_101[k]
                  - f_262 * kh_103[k]
                  + f_263 * kh_274[k]
                  + f_264 * kh_279[k]
                  - f_198 * kh_281[k]
                  + f_263 * kh_288[k]
                  - f_198 * kh_290[k]
                  + f_265 * kh_292[k]
                  + f_203 * kh_463[k]
                  + f_260 * kh_468[k]
                  - f_261 * kh_470[k]
                  + f_203 * kh_477[k]
                  - f_261 * kh_479[k]
                  + f_262 * kh_481[k]
                  - f_263 * kh_505[k]
                  - f_264 * kh_510[k]
                  + f_198 * kh_512[k]
                  - f_263 * kh_519[k]
                  + f_198 * kh_521[k]
                  - f_265 * kh_523[k];
    }

#pragma omp simd aligned(kh_86, kh_91, kh_93, kh_100, kh_102, kh_104, kh_275, kh_280, kh_282, \
                         kh_289, kh_291, kh_293, kh_464, kh_469, kh_471, kh_478, kh_480, \
                         kh_482, kh_506, kh_511, kh_513, kh_520, kh_522, \
                         kh_524 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = -f_216 * kh_86[k]
                  - f_266 * kh_91[k]
                  + f_220 * kh_93[k]
                  - f_216 * kh_100[k]
                  + f_220 * kh_102[k]
                  - f_267 * kh_104[k]
                  + f_268 * kh_275[k]
                  + f_212 * kh_280[k]
                  - f_269 * kh_282[k]
                  + f_268 * kh_289[k]
                  - f_269 * kh_291[k]
                  + f_270 * kh_293[k]
                  + f_216 * kh_464[k]
                  + f_266 * kh_469[k]
                  - f_220 * kh_471[k]
                  + f_216 * kh_478[k]
                  - f_220 * kh_480[k]
                  + f_267 * kh_482[k]
                  - f_268 * kh_506[k]
                  - f_212 * kh_511[k]
                  + f_269 * kh_513[k]
                  - f_268 * kh_520[k]
                  + f_269 * kh_522[k]
                  - f_270 * kh_524[k];
    }

#pragma omp simd aligned(kh_84, kh_87, kh_89, kh_94, kh_96, kh_98, kh_273, kh_276, kh_278, \
                         kh_283, kh_285, kh_287, kh_462, kh_465, kh_467, kh_472, kh_474, \
                         kh_476, kh_504, kh_507, kh_509, kh_514, kh_516, \
                         kh_518 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_39[k] = -f_203 * kh_84[k]
                  - f_260 * kh_87[k]
                  + f_261 * kh_89[k]
                  - f_203 * kh_94[k]
                  + f_261 * kh_96[k]
                  - f_262 * kh_98[k]
                  + f_263 * kh_273[k]
                  + f_264 * kh_276[k]
                  - f_198 * kh_278[k]
                  + f_263 * kh_283[k]
                  - f_198 * kh_285[k]
                  + f_265 * kh_287[k]
                  + f_203 * kh_462[k]
                  + f_260 * kh_465[k]
                  - f_261 * kh_467[k]
                  + f_203 * kh_472[k]
                  - f_261 * kh_474[k]
                  + f_262 * kh_476[k]
                  - f_263 * kh_504[k]
                  - f_264 * kh_507[k]
                  + f_198 * kh_509[k]
                  - f_263 * kh_514[k]
                  + f_198 * kh_516[k]
                  - f_265 * kh_518[k];
    }

#pragma omp simd aligned(kh_86, kh_93, kh_100, kh_102, kh_275, kh_282, kh_289, kh_291, kh_464, \
                         kh_471, kh_478, kh_480, kh_506, kh_513, kh_520, \
                         kh_522 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_40[k] = f_183 * kh_86[k]
                  - f_184 * kh_93[k]
                  - f_183 * kh_100[k]
                  + f_184 * kh_102[k]
                  - f_271 * kh_275[k]
                  + f_258 * kh_282[k]
                  + f_271 * kh_289[k]
                  - f_258 * kh_291[k]
                  - f_183 * kh_464[k]
                  + f_184 * kh_471[k]
                  + f_183 * kh_478[k]
                  - f_184 * kh_480[k]
                  + f_271 * kh_506[k]
                  - f_258 * kh_513[k]
                  - f_271 * kh_520[k]
                  + f_258 * kh_522[k];
    }

#pragma omp simd aligned(kh_84, kh_87, kh_89, kh_94, kh_96, kh_273, kh_276, kh_278, kh_283, \
                         kh_285, kh_462, kh_465, kh_467, kh_472, kh_474, kh_504, kh_507, \
                         kh_509, kh_514, kh_516 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_41[k] = f_167 * kh_84[k]
                  - f_250 * kh_87[k]
                  - f_252 * kh_89[k]
                  - f_160 * kh_94[k]
                  + f_251 * kh_96[k]
                  - f_255 * kh_273[k]
                  + f_253 * kh_276[k]
                  + f_256 * kh_278[k]
                  + f_162 * kh_283[k]
                  - f_254 * kh_285[k]
                  - f_167 * kh_462[k]
                  + f_250 * kh_465[k]
                  + f_252 * kh_467[k]
                  + f_160 * kh_472[k]
                  - f_251 * kh_474[k]
                  + f_255 * kh_504[k]
                  - f_253 * kh_507[k]
                  - f_256 * kh_509[k]
                  - f_162 * kh_514[k]
                  + f_254 * kh_516[k];
    }

#pragma omp simd aligned(kh_86, kh_91, kh_100, kh_275, kh_280, kh_289, kh_464, kh_469, kh_478, \
                         kh_506, kh_511, kh_520 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = -f_272 * kh_86[k]
                  + f_273 * kh_91[k]
                  - f_272 * kh_100[k]
                  + f_274 * kh_275[k]
                  - f_144 * kh_280[k]
                  + f_274 * kh_289[k]
                  + f_272 * kh_464[k]
                  - f_273 * kh_469[k]
                  + f_272 * kh_478[k]
                  - f_274 * kh_506[k]
                  + f_144 * kh_511[k]
                  - f_274 * kh_520[k];
    }

#pragma omp simd aligned(kh_84, kh_87, kh_94, kh_273, kh_276, kh_283, kh_462, kh_465, kh_472, \
                         kh_504, kh_507, kh_514 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_43[k] = -f_244 * kh_84[k]
                  + f_243 * kh_87[k]
                  - f_137 * kh_94[k]
                  + f_247 * kh_273[k]
                  - f_246 * kh_276[k]
                  + f_245 * kh_283[k]
                  + f_244 * kh_462[k]
                  - f_243 * kh_465[k]
                  + f_137 * kh_472[k]
                  - f_247 * kh_504[k]
                  + f_246 * kh_507[k]
                  - f_245 * kh_514[k];
    }

#pragma omp simd aligned(kh_22, kh_27, kh_36, kh_127, kh_132, kh_141, kh_169, kh_174, kh_183, \
                         kh_316, kh_321, kh_330, kh_358, kh_363, kh_372, kh_400, kh_405, \
                         kh_414, kh_589, kh_594, kh_603, kh_631, kh_636, kh_645, kh_673, \
                         kh_678, kh_687 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_44[k] = f_275 * kh_22[k]
                  - f_276 * kh_27[k]
                  + f_277 * kh_36[k]
                  + f_278 * kh_127[k]
                  - f_279 * kh_132[k]
                  + f_280 * kh_141[k]
                  - f_281 * kh_169[k]
                  + f_282 * kh_174[k]
                  - f_283 * kh_183[k]
                  + f_280 * kh_316[k]
                  - f_284 * kh_321[k]
                  + f_285 * kh_330[k]
                  - f_286 * kh_358[k]
                  + f_287 * kh_363[k]
                  - f_288 * kh_372[k]
                  + f_287 * kh_400[k]
                  - f_289 * kh_405[k]
                  + f_290 * kh_414[k]
                  - f_280 * kh_589[k]
                  + f_284 * kh_594[k]
                  - f_285 * kh_603[k]
                  + f_291 * kh_631[k]
                  - f_286 * kh_636[k]
                  + f_292 * kh_645[k]
                  - f_293 * kh_673[k]
                  + f_294 * kh_678[k]
                  - f_295 * kh_687[k];
    }

#pragma omp simd aligned(kh_25, kh_32, kh_130, kh_137, kh_172, kh_179, kh_319, kh_326, kh_361, \
                         kh_368, kh_403, kh_410, kh_592, kh_599, kh_634, kh_641, kh_676, \
                         kh_683 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_45[k] = f_296 * kh_25[k]
                  - f_296 * kh_32[k]
                  + f_297 * kh_130[k]
                  - f_297 * kh_137[k]
                  - f_298 * kh_172[k]
                  + f_298 * kh_179[k]
                  + f_299 * kh_319[k]
                  - f_299 * kh_326[k]
                  - f_300 * kh_361[k]
                  + f_300 * kh_368[k]
                  + f_301 * kh_403[k]
                  - f_301 * kh_410[k]
                  - f_299 * kh_592[k]
                  + f_299 * kh_599[k]
                  + f_302 * kh_634[k]
                  - f_302 * kh_641[k]
                  - f_303 * kh_676[k]
                  + f_303 * kh_683[k];
    }

#pragma omp simd aligned(kh_22, kh_27, kh_29, kh_36, kh_38, kh_127, kh_132, kh_134, kh_141, \
                         kh_143, kh_169, kh_174, kh_176, kh_183, kh_185, kh_316, kh_321, \
                         kh_323, kh_330, kh_332, kh_358, kh_363, kh_365, kh_372, kh_374, \
                         kh_400, kh_405, kh_407, kh_414, kh_416, kh_589, kh_594, kh_596, \
                         kh_603, kh_605, kh_631, kh_636, kh_638, kh_645, kh_647, kh_673, \
                         kh_678, kh_680, kh_687, kh_689 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_46[k] = -f_304 * kh_22[k]
                  - f_305 * kh_27[k]
                  + f_306 * kh_29[k]
                  + f_307 * kh_36[k]
                  - f_308 * kh_38[k]
                  - f_309 * kh_127[k]
                  - f_310 * kh_132[k]
                  + f_311 * kh_134[k]
                  + f_312 * kh_141[k]
                  - f_313 * kh_143[k]
                  + f_314 * kh_169[k]
                  + f_311 * kh_174[k]
                  - f_315 * kh_176[k]
                  - f_316 * kh_183[k]
                  + f_317 * kh_185[k]
                  - f_307 * kh_316[k]
                  - f_318 * kh_321[k]
                  + f_308 * kh_323[k]
                  + f_319 * kh_330[k]
                  - f_320 * kh_332[k]
                  + f_311 * kh_358[k]
                  + f_321 * kh_363[k]
                  - f_322 * kh_365[k]
                  - f_313 * kh_372[k]
                  + f_323 * kh_374[k]
                  - f_324 * kh_400[k]
                  - f_325 * kh_405[k]
                  + f_326 * kh_407[k]
                  + f_321 * kh_414[k]
                  - f_327 * kh_416[k]
                  + f_307 * kh_589[k]
                  + f_318 * kh_594[k]
                  - f_308 * kh_596[k]
                  - f_319 * kh_603[k]
                  + f_320 * kh_605[k]
                  - f_316 * kh_631[k]
                  - f_313 * kh_636[k]
                  + f_317 * kh_638[k]
                  + f_328 * kh_645[k]
                  - f_325 * kh_647[k]
                  + f_321 * kh_673[k]
                  + f_329 * kh_678[k]
                  - f_327 * kh_680[k]
                  - f_330 * kh_687[k]
                  + f_331 * kh_689[k];
    }

#pragma omp simd aligned(kh_25, kh_32, kh_34, kh_130, kh_137, kh_139, kh_172, kh_179, kh_181, \
                         kh_319, kh_326, kh_328, kh_361, kh_368, kh_370, kh_403, kh_410, \
                         kh_412, kh_592, kh_599, kh_601, kh_634, kh_641, kh_643, kh_676, \
                         kh_683, kh_685 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_47[k] = -f_332 * kh_25[k]
                  - f_332 * kh_32[k]
                  + f_333 * kh_34[k]
                  - f_334 * kh_130[k]
                  - f_334 * kh_137[k]
                  + f_335 * kh_139[k]
                  + f_336 * kh_172[k]
                  + f_336 * kh_179[k]
                  - f_337 * kh_181[k]
                  - f_338 * kh_319[k]
                  - f_338 * kh_326[k]
                  + f_339 * kh_328[k]
                  + f_340 * kh_361[k]
                  + f_340 * kh_368[k]
                  - f_341 * kh_370[k]
                  - f_341 * kh_403[k]
                  - f_341 * kh_410[k]
                  + f_342 * kh_412[k]
                  + f_338 * kh_592[k]
                  + f_338 * kh_599[k]
                  - f_339 * kh_601[k]
                  - f_343 * kh_634[k]
                  - f_343 * kh_641[k]
                  + f_340 * kh_643[k]
                  + f_344 * kh_676[k]
                  + f_344 * kh_683[k]
                  - f_345 * kh_685[k];
    }

#pragma omp simd aligned(kh_22, kh_27, kh_29, kh_36, kh_38, kh_40, kh_127, kh_132, kh_134, \
                         kh_141, kh_143, kh_145, kh_169, kh_174, kh_176, kh_183, kh_185, \
                         kh_187, kh_316, kh_321, kh_323, kh_330, kh_332, kh_334, kh_358, \
                         kh_363, kh_365, kh_372, kh_374, kh_376, kh_400, kh_405, kh_407, \
                         kh_414, kh_416, kh_418, kh_589, kh_594, kh_596, kh_603, kh_605, \
                         kh_607, kh_631, kh_636, kh_638, kh_645, kh_647, kh_649, kh_673, \
                         kh_678, kh_680, kh_687, kh_689, kh_691 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_48[k] = f_346 * kh_22[k]
                  + f_347 * kh_27[k]
                  - f_348 * kh_29[k]
                  + f_346 * kh_36[k]
                  - f_348 * kh_38[k]
                  + f_349 * kh_40[k]
                  + f_350 * kh_127[k]
                  + f_351 * kh_132[k]
                  - f_352 * kh_134[k]
                  + f_350 * kh_141[k]
                  - f_352 * kh_143[k]
                  + f_353 * kh_145[k]
                  - f_352 * kh_169[k]
                  - f_354 * kh_174[k]
                  + f_355 * kh_176[k]
                  - f_352 * kh_183[k]
                  + f_355 * kh_185[k]
                  - f_356 * kh_187[k]
                  + f_357 * kh_316[k]
                  + f_358 * kh_321[k]
                  - f_359 * kh_323[k]
                  + f_357 * kh_330[k]
                  - f_359 * kh_332[k]
                  + f_360 * kh_334[k]
                  - f_353 * kh_358[k]
                  - f_361 * kh_363[k]
                  + f_356 * kh_365[k]
                  - f_353 * kh_372[k]
                  + f_356 * kh_374[k]
                  - f_362 * kh_376[k]
                  + f_361 * kh_400[k]
                  + f_363 * kh_405[k]
                  - f_364 * kh_407[k]
                  + f_361 * kh_414[k]
                  - f_364 * kh_416[k]
                  + f_365 * kh_418[k]
                  - f_357 * kh_589[k]
                  - f_358 * kh_594[k]
                  + f_359 * kh_596[k]
                  - f_357 * kh_603[k]
                  + f_359 * kh_605[k]
                  - f_360 * kh_607[k]
                  + f_366 * kh_631[k]
                  + f_353 * kh_636[k]
                  - f_367 * kh_638[k]
                  + f_366 * kh_645[k]
                  - f_367 * kh_647[k]
                  + f_363 * kh_649[k]
                  - f_368 * kh_673[k]
                  - f_369 * kh_678[k]
                  + f_362 * kh_680[k]
                  - f_368 * kh_687[k]
                  + f_362 * kh_689[k]
                  - f_370 * kh_691[k];
    }

#pragma omp simd aligned(kh_23, kh_28, kh_30, kh_37, kh_39, kh_41, kh_128, kh_133, kh_135, \
                         kh_142, kh_144, kh_146, kh_170, kh_175, kh_177, kh_184, kh_186, \
                         kh_188, kh_317, kh_322, kh_324, kh_331, kh_333, kh_335, kh_359, \
                         kh_364, kh_366, kh_373, kh_375, kh_377, kh_401, kh_406, kh_408, \
                         kh_415, kh_417, kh_419, kh_590, kh_595, kh_597, kh_604, kh_606, \
                         kh_608, kh_632, kh_637, kh_639, kh_646, kh_648, kh_650, kh_674, \
                         kh_679, kh_681, kh_688, kh_690, kh_692 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_49[k] = f_371 * kh_23[k]
                  + f_372 * kh_28[k]
                  - f_373 * kh_30[k]
                  + f_371 * kh_37[k]
                  - f_373 * kh_39[k]
                  + f_374 * kh_41[k]
                  + f_375 * kh_128[k]
                  + f_376 * kh_133[k]
                  - f_377 * kh_135[k]
                  + f_375 * kh_142[k]
                  - f_377 * kh_144[k]
                  + f_378 * kh_146[k]
                  - f_379 * kh_170[k]
                  - f_380 * kh_175[k]
                  + f_381 * kh_177[k]
                  - f_379 * kh_184[k]
                  + f_381 * kh_186[k]
                  - f_382 * kh_188[k]
                  + f_383 * kh_317[k]
                  + f_384 * kh_322[k]
                  - f_378 * kh_324[k]
                  + f_383 * kh_331[k]
                  - f_378 * kh_333[k]
                  + f_385 * kh_335[k]
                  - f_386 * kh_359[k]
                  - f_387 * kh_364[k]
                  + f_388 * kh_366[k]
                  - f_386 * kh_373[k]
                  + f_388 * kh_375[k]
                  - f_389 * kh_377[k]
                  + f_387 * kh_401[k]
                  + f_381 * kh_406[k]
                  - f_390 * kh_408[k]
                  + f_387 * kh_415[k]
                  - f_390 * kh_417[k]
                  + f_391 * kh_419[k]
                  - f_383 * kh_590[k]
                  - f_384 * kh_595[k]
                  + f_378 * kh_597[k]
                  - f_383 * kh_604[k]
                  + f_378 * kh_606[k]
                  - f_385 * kh_608[k]
                  + f_392 * kh_632[k]
                  + f_386 * kh_637[k]
                  - f_393 * kh_639[k]
                  + f_392 * kh_646[k]
                  - f_393 * kh_648[k]
                  + f_394 * kh_650[k]
                  - f_395 * kh_674[k]
                  - f_393 * kh_679[k]
                  + f_396 * kh_681[k]
                  - f_395 * kh_688[k]
                  + f_396 * kh_690[k]
                  - f_397 * kh_692[k];
    }

#pragma omp simd aligned(kh_21, kh_24, kh_26, kh_31, kh_33, kh_35, kh_126, kh_129, kh_131, \
                         kh_136, kh_138, kh_140, kh_168, kh_171, kh_173, kh_178, kh_180, \
                         kh_182, kh_315, kh_318, kh_320, kh_325, kh_327, kh_329, kh_357, \
                         kh_360, kh_362, kh_367, kh_369, kh_371, kh_399, kh_402, kh_404, \
                         kh_409, kh_411, kh_413, kh_588, kh_591, kh_593, kh_598, kh_600, \
                         kh_602, kh_630, kh_633, kh_635, kh_640, kh_642, kh_644, kh_672, \
                         kh_675, kh_677, kh_682, kh_684, kh_686 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_50[k] = f_346 * kh_21[k]
                  + f_347 * kh_24[k]
                  - f_348 * kh_26[k]
                  + f_346 * kh_31[k]
                  - f_348 * kh_33[k]
                  + f_349 * kh_35[k]
                  + f_350 * kh_126[k]
                  + f_351 * kh_129[k]
                  - f_352 * kh_131[k]
                  + f_350 * kh_136[k]
                  - f_352 * kh_138[k]
                  + f_353 * kh_140[k]
                  - f_352 * kh_168[k]
                  - f_354 * kh_171[k]
                  + f_355 * kh_173[k]
                  - f_352 * kh_178[k]
                  + f_355 * kh_180[k]
                  - f_356 * kh_182[k]
                  + f_357 * kh_315[k]
                  + f_358 * kh_318[k]
                  - f_359 * kh_320[k]
                  + f_357 * kh_325[k]
                  - f_359 * kh_327[k]
                  + f_360 * kh_329[k]
                  - f_353 * kh_357[k]
                  - f_361 * kh_360[k]
                  + f_356 * kh_362[k]
                  - f_353 * kh_367[k]
                  + f_356 * kh_369[k]
                  - f_362 * kh_371[k]
                  + f_361 * kh_399[k]
                  + f_363 * kh_402[k]
                  - f_364 * kh_404[k]
                  + f_361 * kh_409[k]
                  - f_364 * kh_411[k]
                  + f_365 * kh_413[k]
                  - f_357 * kh_588[k]
                  - f_358 * kh_591[k]
                  + f_359 * kh_593[k]
                  - f_357 * kh_598[k]
                  + f_359 * kh_600[k]
                  - f_360 * kh_602[k]
                  + f_366 * kh_630[k]
                  + f_353 * kh_633[k]
                  - f_367 * kh_635[k]
                  + f_366 * kh_640[k]
                  - f_367 * kh_642[k]
                  + f_363 * kh_644[k]
                  - f_368 * kh_672[k]
                  - f_369 * kh_675[k]
                  + f_362 * kh_677[k]
                  - f_368 * kh_682[k]
                  + f_362 * kh_684[k]
                  - f_370 * kh_686[k];
    }

#pragma omp simd aligned(kh_23, kh_30, kh_37, kh_39, kh_128, kh_135, kh_142, kh_144, kh_170, \
                         kh_177, kh_184, kh_186, kh_317, kh_324, kh_331, kh_333, kh_359, \
                         kh_366, kh_373, kh_375, kh_401, kh_408, kh_415, kh_417, kh_590, \
                         kh_597, kh_604, kh_606, kh_632, kh_639, kh_646, kh_648, kh_674, \
                         kh_681, kh_688, kh_690 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_51[k] = -f_398 * kh_23[k]
                  + f_332 * kh_30[k]
                  + f_398 * kh_37[k]
                  - f_332 * kh_39[k]
                  - f_399 * kh_128[k]
                  + f_334 * kh_135[k]
                  + f_399 * kh_142[k]
                  - f_334 * kh_144[k]
                  + f_400 * kh_170[k]
                  - f_336 * kh_177[k]
                  - f_400 * kh_184[k]
                  + f_336 * kh_186[k]
                  - f_401 * kh_317[k]
                  + f_338 * kh_324[k]
                  + f_401 * kh_331[k]
                  - f_338 * kh_333[k]
                  + f_343 * kh_359[k]
                  - f_340 * kh_366[k]
                  - f_343 * kh_373[k]
                  + f_340 * kh_375[k]
                  - f_340 * kh_401[k]
                  + f_341 * kh_408[k]
                  + f_340 * kh_415[k]
                  - f_341 * kh_417[k]
                  + f_401 * kh_590[k]
                  - f_338 * kh_597[k]
                  - f_401 * kh_604[k]
                  + f_338 * kh_606[k]
                  - f_335 * kh_632[k]
                  + f_343 * kh_639[k]
                  + f_335 * kh_646[k]
                  - f_343 * kh_648[k]
                  + f_402 * kh_674[k]
                  - f_344 * kh_681[k]
                  - f_402 * kh_688[k]
                  + f_344 * kh_690[k];
    }

#pragma omp simd aligned(kh_21, kh_24, kh_26, kh_31, kh_33, kh_126, kh_129, kh_131, kh_136, \
                         kh_138, kh_168, kh_171, kh_173, kh_178, kh_180, kh_315, kh_318, \
                         kh_320, kh_325, kh_327, kh_357, kh_360, kh_362, kh_367, kh_369, \
                         kh_399, kh_402, kh_404, kh_409, kh_411, kh_588, kh_591, kh_593, \
                         kh_598, kh_600, kh_630, kh_633, kh_635, kh_640, kh_642, kh_672, \
                         kh_675, kh_677, kh_682, kh_684 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_52[k] = -f_307 * kh_21[k]
                  + f_305 * kh_24[k]
                  + f_308 * kh_26[k]
                  + f_304 * kh_31[k]
                  - f_306 * kh_33[k]
                  - f_312 * kh_126[k]
                  + f_310 * kh_129[k]
                  + f_313 * kh_131[k]
                  + f_309 * kh_136[k]
                  - f_311 * kh_138[k]
                  + f_316 * kh_168[k]
                  - f_311 * kh_171[k]
                  - f_317 * kh_173[k]
                  - f_314 * kh_178[k]
                  + f_315 * kh_180[k]
                  - f_319 * kh_315[k]
                  + f_318 * kh_318[k]
                  + f_320 * kh_320[k]
                  + f_307 * kh_325[k]
                  - f_308 * kh_327[k]
                  + f_313 * kh_357[k]
                  - f_321 * kh_360[k]
                  - f_323 * kh_362[k]
                  - f_311 * kh_367[k]
                  + f_322 * kh_369[k]
                  - f_321 * kh_399[k]
                  + f_325 * kh_402[k]
                  + f_327 * kh_404[k]
                  + f_324 * kh_409[k]
                  - f_326 * kh_411[k]
                  + f_319 * kh_588[k]
                  - f_318 * kh_591[k]
                  - f_320 * kh_593[k]
                  - f_307 * kh_598[k]
                  + f_308 * kh_600[k]
                  - f_328 * kh_630[k]
                  + f_313 * kh_633[k]
                  + f_325 * kh_635[k]
                  + f_316 * kh_640[k]
                  - f_317 * kh_642[k]
                  + f_330 * kh_672[k]
                  - f_329 * kh_675[k]
                  - f_331 * kh_677[k]
                  - f_321 * kh_682[k]
                  + f_327 * kh_684[k];
    }

#pragma omp simd aligned(kh_23, kh_28, kh_37, kh_128, kh_133, kh_142, kh_170, kh_175, kh_184, \
                         kh_317, kh_322, kh_331, kh_359, kh_364, kh_373, kh_401, kh_406, \
                         kh_415, kh_590, kh_595, kh_604, kh_632, kh_637, kh_646, kh_674, \
                         kh_679, kh_688 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_53[k] = f_403 * kh_23[k]
                  - f_404 * kh_28[k]
                  + f_403 * kh_37[k]
                  + f_405 * kh_128[k]
                  - f_406 * kh_133[k]
                  + f_405 * kh_142[k]
                  - f_407 * kh_170[k]
                  + f_408 * kh_175[k]
                  - f_407 * kh_184[k]
                  + f_409 * kh_317[k]
                  - f_410 * kh_322[k]
                  + f_409 * kh_331[k]
                  - f_411 * kh_359[k]
                  + f_298 * kh_364[k]
                  - f_411 * kh_373[k]
                  + f_302 * kh_401[k]
                  - f_412 * kh_406[k]
                  + f_302 * kh_415[k]
                  - f_409 * kh_590[k]
                  + f_410 * kh_595[k]
                  - f_409 * kh_604[k]
                  + f_297 * kh_632[k]
                  - f_413 * kh_637[k]
                  + f_297 * kh_646[k]
                  - f_414 * kh_674[k]
                  + f_300 * kh_679[k]
                  - f_414 * kh_688[k];
    }

#pragma omp simd aligned(kh_21, kh_24, kh_31, kh_126, kh_129, kh_136, kh_168, kh_171, kh_178, \
                         kh_315, kh_318, kh_325, kh_357, kh_360, kh_367, kh_399, kh_402, \
                         kh_409, kh_588, kh_591, kh_598, kh_630, kh_633, kh_640, kh_672, \
                         kh_675, kh_682 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_54[k] = f_277 * kh_21[k]
                  - f_276 * kh_24[k]
                  + f_275 * kh_31[k]
                  + f_280 * kh_126[k]
                  - f_279 * kh_129[k]
                  + f_278 * kh_136[k]
                  - f_283 * kh_168[k]
                  + f_282 * kh_171[k]
                  - f_281 * kh_178[k]
                  + f_285 * kh_315[k]
                  - f_284 * kh_318[k]
                  + f_280 * kh_325[k]
                  - f_288 * kh_357[k]
                  + f_287 * kh_360[k]
                  - f_286 * kh_367[k]
                  + f_290 * kh_399[k]
                  - f_289 * kh_402[k]
                  + f_287 * kh_409[k]
                  - f_285 * kh_588[k]
                  + f_284 * kh_591[k]
                  - f_280 * kh_598[k]
                  + f_292 * kh_630[k]
                  - f_286 * kh_633[k]
                  + f_291 * kh_640[k]
                  - f_295 * kh_672[k]
                  + f_294 * kh_675[k]
                  - f_293 * kh_682[k];
    }

#pragma omp simd aligned(kh_85, kh_90, kh_99, kh_232, kh_237, kh_246, kh_274, kh_279, kh_288, \
                         kh_463, kh_468, kh_477, kh_505, kh_510, kh_519, kh_547, kh_552, \
                         kh_561 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_55[k] = f_415 * kh_85[k]
                  - f_416 * kh_90[k]
                  + f_417 * kh_99[k]
                  + f_416 * kh_232[k]
                  - f_418 * kh_237[k]
                  + f_419 * kh_246[k]
                  - f_420 * kh_274[k]
                  + f_421 * kh_279[k]
                  - f_422 * kh_288[k]
                  + f_415 * kh_463[k]
                  - f_416 * kh_468[k]
                  + f_417 * kh_477[k]
                  - f_420 * kh_505[k]
                  + f_421 * kh_510[k]
                  - f_422 * kh_519[k]
                  + f_423 * kh_547[k]
                  - f_424 * kh_552[k]
                  + f_425 * kh_561[k];
    }

#pragma omp simd aligned(kh_88, kh_95, kh_235, kh_242, kh_277, kh_284, kh_466, kh_473, kh_508, \
                         kh_515, kh_550, kh_557 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_56[k] = f_324 * kh_88[k]
                  - f_324 * kh_95[k]
                  + f_317 * kh_235[k]
                  - f_317 * kh_242[k]
                  - f_426 * kh_277[k]
                  + f_426 * kh_284[k]
                  + f_324 * kh_466[k]
                  - f_324 * kh_473[k]
                  - f_426 * kh_508[k]
                  + f_426 * kh_515[k]
                  + f_427 * kh_550[k]
                  - f_427 * kh_557[k];
    }

#pragma omp simd aligned(kh_85, kh_90, kh_92, kh_99, kh_101, kh_232, kh_237, kh_239, kh_246, \
                         kh_248, kh_274, kh_279, kh_281, kh_288, kh_290, kh_463, kh_468, \
                         kh_470, kh_477, kh_479, kh_505, kh_510, kh_512, kh_519, kh_521, \
                         kh_547, kh_552, kh_554, kh_561, kh_563 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_57[k] = -f_428 * kh_85[k]
                  - f_429 * kh_90[k]
                  + f_302 * kh_92[k]
                  + f_430 * kh_99[k]
                  - f_414 * kh_101[k]
                  - f_297 * kh_232[k]
                  - f_431 * kh_237[k]
                  + f_300 * kh_239[k]
                  + f_429 * kh_246[k]
                  - f_432 * kh_248[k]
                  + f_432 * kh_274[k]
                  + f_433 * kh_279[k]
                  - f_434 * kh_281[k]
                  - f_435 * kh_288[k]
                  + f_436 * kh_290[k]
                  - f_428 * kh_463[k]
                  - f_429 * kh_468[k]
                  + f_302 * kh_470[k]
                  + f_430 * kh_477[k]
                  - f_414 * kh_479[k]
                  + f_432 * kh_505[k]
                  + f_433 * kh_510[k]
                  - f_434 * kh_512[k]
                  - f_435 * kh_519[k]
                  + f_436 * kh_521[k]
                  - f_437 * kh_547[k]
                  - f_438 * kh_552[k]
                  + f_439 * kh_554[k]
                  + f_440 * kh_561[k]
                  - f_441 * kh_563[k];
    }

#pragma omp simd aligned(kh_88, kh_95, kh_97, kh_235, kh_242, kh_244, kh_277, kh_284, kh_286, \
                         kh_466, kh_473, kh_475, kh_508, kh_515, kh_517, kh_550, kh_557, \
                         kh_559 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_58[k] = -f_442 * kh_88[k]
                  - f_442 * kh_95[k]
                  + f_443 * kh_97[k]
                  - f_443 * kh_235[k]
                  - f_443 * kh_242[k]
                  + f_444 * kh_244[k]
                  + f_445 * kh_277[k]
                  + f_445 * kh_284[k]
                  - f_446 * kh_286[k]
                  - f_442 * kh_466[k]
                  - f_442 * kh_473[k]
                  + f_443 * kh_475[k]
                  + f_445 * kh_508[k]
                  + f_445 * kh_515[k]
                  - f_446 * kh_517[k]
                  - f_447 * kh_550[k]
                  - f_447 * kh_557[k]
                  + f_448 * kh_559[k];
    }

#pragma omp simd aligned(kh_85, kh_90, kh_92, kh_99, kh_101, kh_103, kh_232, kh_237, kh_239, \
                         kh_246, kh_248, kh_250, kh_274, kh_279, kh_281, kh_288, kh_290, \
                         kh_292, kh_463, kh_468, kh_470, kh_477, kh_479, kh_481, kh_505, \
                         kh_510, kh_512, kh_519, kh_521, kh_523, kh_547, kh_552, kh_554, \
                         kh_561, kh_563, kh_565 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_59[k] = f_449 * kh_85[k]
                  + f_450 * kh_90[k]
                  - f_451 * kh_92[k]
                  + f_449 * kh_99[k]
                  - f_451 * kh_101[k]
                  + f_452 * kh_103[k]
                  + f_450 * kh_232[k]
                  + f_453 * kh_237[k]
                  - f_454 * kh_239[k]
                  + f_450 * kh_246[k]
                  - f_454 * kh_248[k]
                  + f_455 * kh_250[k]
                  - f_456 * kh_274[k]
                  - f_457 * kh_279[k]
                  + f_458 * kh_281[k]
                  - f_456 * kh_288[k]
                  + f_458 * kh_290[k]
                  - f_459 * kh_292[k]
                  + f_449 * kh_463[k]
                  + f_450 * kh_468[k]
                  - f_451 * kh_470[k]
                  + f_449 * kh_477[k]
                  - f_451 * kh_479[k]
                  + f_452 * kh_481[k]
                  - f_456 * kh_505[k]
                  - f_457 * kh_510[k]
                  + f_458 * kh_512[k]
                  - f_456 * kh_519[k]
                  + f_458 * kh_521[k]
                  - f_459 * kh_523[k]
                  + f_460 * kh_547[k]
                  + f_461 * kh_552[k]
                  - f_462 * kh_554[k]
                  + f_460 * kh_561[k]
                  - f_462 * kh_563[k]
                  + f_463 * kh_565[k];
    }

#pragma omp simd aligned(kh_86, kh_91, kh_93, kh_100, kh_102, kh_104, kh_233, kh_238, kh_240, \
                         kh_247, kh_249, kh_251, kh_275, kh_280, kh_282, kh_289, kh_291, \
                         kh_293, kh_464, kh_469, kh_471, kh_478, kh_480, kh_482, kh_506, \
                         kh_511, kh_513, kh_520, kh_522, kh_524, kh_548, kh_553, kh_555, \
                         kh_562, kh_564, kh_566 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_60[k] = f_464 * kh_86[k]
                  + f_465 * kh_91[k]
                  - f_466 * kh_93[k]
                  + f_464 * kh_100[k]
                  - f_466 * kh_102[k]
                  + f_467 * kh_104[k]
                  + f_465 * kh_233[k]
                  + f_468 * kh_238[k]
                  - f_469 * kh_240[k]
                  + f_465 * kh_247[k]
                  - f_469 * kh_249[k]
                  + f_470 * kh_251[k]
                  - f_469 * kh_275[k]
                  - f_471 * kh_280[k]
                  + f_472 * kh_282[k]
                  - f_469 * kh_289[k]
                  + f_472 * kh_291[k]
                  - f_473 * kh_293[k]
                  + f_464 * kh_464[k]
                  + f_465 * kh_469[k]
                  - f_466 * kh_471[k]
                  + f_464 * kh_478[k]
                  - f_466 * kh_480[k]
                  + f_467 * kh_482[k]
                  - f_469 * kh_506[k]
                  - f_471 * kh_511[k]
                  + f_472 * kh_513[k]
                  - f_469 * kh_520[k]
                  + f_472 * kh_522[k]
                  - f_473 * kh_524[k]
                  + f_474 * kh_548[k]
                  + f_475 * kh_553[k]
                  - f_476 * kh_555[k]
                  + f_474 * kh_562[k]
                  - f_476 * kh_564[k]
                  + f_477 * kh_566[k];
    }

#pragma omp simd aligned(kh_84, kh_87, kh_89, kh_94, kh_96, kh_98, kh_231, kh_234, kh_236, \
                         kh_241, kh_243, kh_245, kh_273, kh_276, kh_278, kh_283, kh_285, \
                         kh_287, kh_462, kh_465, kh_467, kh_472, kh_474, kh_476, kh_504, \
                         kh_507, kh_509, kh_514, kh_516, kh_518, kh_546, kh_549, kh_551, \
                         kh_556, kh_558, kh_560 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_61[k] = f_449 * kh_84[k]
                  + f_450 * kh_87[k]
                  - f_451 * kh_89[k]
                  + f_449 * kh_94[k]
                  - f_451 * kh_96[k]
                  + f_452 * kh_98[k]
                  + f_450 * kh_231[k]
                  + f_453 * kh_234[k]
                  - f_454 * kh_236[k]
                  + f_450 * kh_241[k]
                  - f_454 * kh_243[k]
                  + f_455 * kh_245[k]
                  - f_456 * kh_273[k]
                  - f_457 * kh_276[k]
                  + f_458 * kh_278[k]
                  - f_456 * kh_283[k]
                  + f_458 * kh_285[k]
                  - f_459 * kh_287[k]
                  + f_449 * kh_462[k]
                  + f_450 * kh_465[k]
                  - f_451 * kh_467[k]
                  + f_449 * kh_472[k]
                  - f_451 * kh_474[k]
                  + f_452 * kh_476[k]
                  - f_456 * kh_504[k]
                  - f_457 * kh_507[k]
                  + f_458 * kh_509[k]
                  - f_456 * kh_514[k]
                  + f_458 * kh_516[k]
                  - f_459 * kh_518[k]
                  + f_460 * kh_546[k]
                  + f_461 * kh_549[k]
                  - f_462 * kh_551[k]
                  + f_460 * kh_556[k]
                  - f_462 * kh_558[k]
                  + f_463 * kh_560[k];
    }

#pragma omp simd aligned(kh_86, kh_93, kh_100, kh_102, kh_233, kh_240, kh_247, kh_249, kh_275, \
                         kh_282, kh_289, kh_291, kh_464, kh_471, kh_478, kh_480, kh_506, \
                         kh_513, kh_520, kh_522, kh_548, kh_555, kh_562, \
                         kh_564 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_62[k] = -f_478 * kh_86[k]
                  + f_442 * kh_93[k]
                  + f_478 * kh_100[k]
                  - f_442 * kh_102[k]
                  - f_442 * kh_233[k]
                  + f_443 * kh_240[k]
                  + f_442 * kh_247[k]
                  - f_443 * kh_249[k]
                  + f_479 * kh_275[k]
                  - f_445 * kh_282[k]
                  - f_479 * kh_289[k]
                  + f_445 * kh_291[k]
                  - f_478 * kh_464[k]
                  + f_442 * kh_471[k]
                  + f_478 * kh_478[k]
                  - f_442 * kh_480[k]
                  + f_479 * kh_506[k]
                  - f_445 * kh_513[k]
                  - f_479 * kh_520[k]
                  + f_445 * kh_522[k]
                  - f_480 * kh_548[k]
                  + f_447 * kh_555[k]
                  + f_480 * kh_562[k]
                  - f_447 * kh_564[k];
    }

#pragma omp simd aligned(kh_84, kh_87, kh_89, kh_94, kh_96, kh_231, kh_234, kh_236, kh_241, \
                         kh_243, kh_273, kh_276, kh_278, kh_283, kh_285, kh_462, kh_465, \
                         kh_467, kh_472, kh_474, kh_504, kh_507, kh_509, kh_514, kh_516, \
                         kh_546, kh_549, kh_551, kh_556, kh_558 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_63[k] = -f_430 * kh_84[k]
                  + f_429 * kh_87[k]
                  + f_414 * kh_89[k]
                  + f_428 * kh_94[k]
                  - f_302 * kh_96[k]
                  - f_429 * kh_231[k]
                  + f_431 * kh_234[k]
                  + f_432 * kh_236[k]
                  + f_297 * kh_241[k]
                  - f_300 * kh_243[k]
                  + f_435 * kh_273[k]
                  - f_433 * kh_276[k]
                  - f_436 * kh_278[k]
                  - f_432 * kh_283[k]
                  + f_434 * kh_285[k]
                  - f_430 * kh_462[k]
                  + f_429 * kh_465[k]
                  + f_414 * kh_467[k]
                  + f_428 * kh_472[k]
                  - f_302 * kh_474[k]
                  + f_435 * kh_504[k]
                  - f_433 * kh_507[k]
                  - f_436 * kh_509[k]
                  - f_432 * kh_514[k]
                  + f_434 * kh_516[k]
                  - f_440 * kh_546[k]
                  + f_438 * kh_549[k]
                  + f_441 * kh_551[k]
                  + f_437 * kh_556[k]
                  - f_439 * kh_558[k];
    }

#pragma omp simd aligned(kh_86, kh_91, kh_100, kh_233, kh_238, kh_247, kh_275, kh_280, kh_289, \
                         kh_464, kh_469, kh_478, kh_506, kh_511, kh_520, kh_548, kh_553, \
                         kh_562 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_64[k] = f_316 * kh_86[k]
                  - f_481 * kh_91[k]
                  + f_316 * kh_100[k]
                  + f_311 * kh_233[k]
                  - f_482 * kh_238[k]
                  + f_311 * kh_247[k]
                  - f_323 * kh_275[k]
                  + f_326 * kh_280[k]
                  - f_323 * kh_289[k]
                  + f_316 * kh_464[k]
                  - f_481 * kh_469[k]
                  + f_316 * kh_478[k]
                  - f_323 * kh_506[k]
                  + f_326 * kh_511[k]
                  - f_323 * kh_520[k]
                  + f_483 * kh_548[k]
                  - f_484 * kh_553[k]
                  + f_483 * kh_562[k];
    }

#pragma omp simd aligned(kh_84, kh_87, kh_94, kh_231, kh_234, kh_241, kh_273, kh_276, kh_283, \
                         kh_462, kh_465, kh_472, kh_504, kh_507, kh_514, kh_546, kh_549, \
                         kh_556 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_65[k] = f_417 * kh_84[k]
                  - f_416 * kh_87[k]
                  + f_415 * kh_94[k]
                  + f_419 * kh_231[k]
                  - f_418 * kh_234[k]
                  + f_416 * kh_241[k]
                  - f_422 * kh_273[k]
                  + f_421 * kh_276[k]
                  - f_420 * kh_283[k]
                  + f_417 * kh_462[k]
                  - f_416 * kh_465[k]
                  + f_415 * kh_472[k]
                  - f_422 * kh_504[k]
                  + f_421 * kh_507[k]
                  - f_420 * kh_514[k]
                  + f_425 * kh_546[k]
                  - f_424 * kh_549[k]
                  + f_423 * kh_556[k];
    }

#pragma omp simd aligned(kh_22, kh_27, kh_36, kh_127, kh_132, kh_141, kh_169, kh_174, kh_183, \
                         kh_316, kh_321, kh_330, kh_358, kh_363, kh_372, kh_400, kh_405, \
                         kh_414, kh_589, kh_594, kh_603, kh_631, kh_636, kh_645, kh_673, \
                         kh_678, kh_687, kh_715, kh_720, kh_729 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_66[k] = -f_485 * kh_22[k]
                  + f_486 * kh_27[k]
                  - f_487 * kh_36[k]
                  - f_488 * kh_127[k]
                  + f_489 * kh_132[k]
                  - f_490 * kh_141[k]
                  + f_491 * kh_169[k]
                  - f_492 * kh_174[k]
                  + f_493 * kh_183[k]
                  - f_488 * kh_316[k]
                  + f_489 * kh_321[k]
                  - f_490 * kh_330[k]
                  + f_492 * kh_358[k]
                  - f_494 * kh_363[k]
                  + f_495 * kh_372[k]
                  - f_492 * kh_400[k]
                  + f_494 * kh_405[k]
                  - f_495 * kh_414[k]
                  - f_485 * kh_589[k]
                  + f_486 * kh_594[k]
                  - f_487 * kh_603[k]
                  + f_491 * kh_631[k]
                  - f_492 * kh_636[k]
                  + f_493 * kh_645[k]
                  - f_492 * kh_673[k]
                  + f_494 * kh_678[k]
                  - f_495 * kh_687[k]
                  + f_496 * kh_715[k]
                  - f_497 * kh_720[k]
                  + f_498 * kh_729[k];
    }

#pragma omp simd aligned(kh_25, kh_32, kh_130, kh_137, kh_172, kh_179, kh_319, kh_326, kh_361, \
                         kh_368, kh_403, kh_410, kh_592, kh_599, kh_634, kh_641, kh_676, \
                         kh_683, kh_718, kh_725 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_67[k] = -f_499 * kh_25[k]
                  + f_499 * kh_32[k]
                  - f_334 * kh_130[k]
                  + f_334 * kh_137[k]
                  + f_340 * kh_172[k]
                  - f_340 * kh_179[k]
                  - f_334 * kh_319[k]
                  + f_334 * kh_326[k]
                  + f_341 * kh_361[k]
                  - f_341 * kh_368[k]
                  - f_341 * kh_403[k]
                  + f_341 * kh_410[k]
                  - f_499 * kh_592[k]
                  + f_499 * kh_599[k]
                  + f_340 * kh_634[k]
                  - f_340 * kh_641[k]
                  - f_341 * kh_676[k]
                  + f_341 * kh_683[k]
                  + f_500 * kh_718[k]
                  - f_500 * kh_725[k];
    }

#pragma omp simd aligned(kh_22, kh_27, kh_29, kh_36, kh_38, kh_127, kh_132, kh_134, kh_141, \
                         kh_143, kh_169, kh_174, kh_176, kh_183, kh_185, kh_316, kh_321, \
                         kh_323, kh_330, kh_332, kh_358, kh_363, kh_365, kh_372, kh_374, \
                         kh_400, kh_405, kh_407, kh_414, kh_416, kh_589, kh_594, kh_596, \
                         kh_603, kh_605, kh_631, kh_636, kh_638, kh_645, kh_647, kh_673, \
                         kh_678, kh_680, kh_687, kh_689, kh_715, kh_720, kh_722, kh_729, \
                         kh_731 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_68[k] = f_501 * kh_22[k]
                  + f_502 * kh_27[k]
                  - f_503 * kh_29[k]
                  - f_504 * kh_36[k]
                  + f_505 * kh_38[k]
                  + f_506 * kh_127[k]
                  + f_507 * kh_132[k]
                  - f_478 * kh_134[k]
                  - f_501 * kh_141[k]
                  + f_503 * kh_143[k]
                  - f_478 * kh_169[k]
                  - f_508 * kh_174[k]
                  + f_444 * kh_176[k]
                  + f_503 * kh_183[k]
                  - f_509 * kh_185[k]
                  + f_506 * kh_316[k]
                  + f_507 * kh_321[k]
                  - f_478 * kh_323[k]
                  - f_501 * kh_330[k]
                  + f_503 * kh_332[k]
                  - f_442 * kh_358[k]
                  - f_510 * kh_363[k]
                  + f_511 * kh_365[k]
                  + f_508 * kh_372[k]
                  - f_479 * kh_374[k]
                  + f_442 * kh_400[k]
                  + f_510 * kh_405[k]
                  - f_511 * kh_407[k]
                  - f_508 * kh_414[k]
                  + f_479 * kh_416[k]
                  + f_501 * kh_589[k]
                  + f_502 * kh_594[k]
                  - f_503 * kh_596[k]
                  - f_504 * kh_603[k]
                  + f_505 * kh_605[k]
                  - f_478 * kh_631[k]
                  - f_508 * kh_636[k]
                  + f_444 * kh_638[k]
                  + f_503 * kh_645[k]
                  - f_509 * kh_647[k]
                  + f_442 * kh_673[k]
                  + f_510 * kh_678[k]
                  - f_511 * kh_680[k]
                  - f_508 * kh_687[k]
                  + f_479 * kh_689[k]
                  - f_512 * kh_715[k]
                  - f_513 * kh_720[k]
                  + f_514 * kh_722[k]
                  + f_515 * kh_729[k]
                  - f_516 * kh_731[k];
    }

#pragma omp simd aligned(kh_25, kh_32, kh_34, kh_130, kh_137, kh_139, kh_172, kh_179, kh_181, \
                         kh_319, kh_326, kh_328, kh_361, kh_368, kh_370, kh_403, kh_410, \
                         kh_412, kh_592, kh_599, kh_601, kh_634, kh_641, kh_643, kh_676, \
                         kh_683, kh_685, kh_718, kh_725, kh_727 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_69[k] = f_517 * kh_25[k]
                  + f_517 * kh_32[k]
                  - f_518 * kh_34[k]
                  + f_429 * kh_130[k]
                  + f_429 * kh_137[k]
                  - f_431 * kh_139[k]
                  - f_432 * kh_172[k]
                  - f_432 * kh_179[k]
                  + f_303 * kh_181[k]
                  + f_429 * kh_319[k]
                  + f_429 * kh_326[k]
                  - f_431 * kh_328[k]
                  - f_303 * kh_361[k]
                  - f_303 * kh_368[k]
                  + f_519 * kh_370[k]
                  + f_303 * kh_403[k]
                  + f_303 * kh_410[k]
                  - f_519 * kh_412[k]
                  + f_517 * kh_592[k]
                  + f_517 * kh_599[k]
                  - f_518 * kh_601[k]
                  - f_432 * kh_634[k]
                  - f_432 * kh_641[k]
                  + f_303 * kh_643[k]
                  + f_303 * kh_676[k]
                  + f_303 * kh_683[k]
                  - f_519 * kh_685[k]
                  - f_520 * kh_718[k]
                  - f_520 * kh_725[k]
                  + f_521 * kh_727[k];
    }

#pragma omp simd aligned(kh_22, kh_27, kh_29, kh_36, kh_38, kh_40, kh_127, kh_132, kh_134, \
                         kh_141, kh_143, kh_145, kh_169, kh_174, kh_176, kh_183, kh_185, \
                         kh_187, kh_316, kh_321, kh_323, kh_330, kh_332, kh_334, kh_358, \
                         kh_363, kh_365, kh_372, kh_374, kh_376, kh_400, kh_405, kh_407, \
                         kh_414, kh_416, kh_418, kh_589, kh_594, kh_596, kh_603, kh_605, \
                         kh_607, kh_631, kh_636, kh_638, kh_645, kh_647, kh_649, kh_673, \
                         kh_678, kh_680, kh_687, kh_689, kh_691, kh_715, kh_720, kh_722, \
                         kh_729, kh_731, kh_733 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_70[k] = -f_522 * kh_22[k]
                  - f_523 * kh_27[k]
                  + f_524 * kh_29[k]
                  - f_522 * kh_36[k]
                  + f_524 * kh_38[k]
                  - f_525 * kh_40[k]
                  - f_526 * kh_127[k]
                  - f_527 * kh_132[k]
                  + f_528 * kh_134[k]
                  - f_526 * kh_141[k]
                  + f_528 * kh_143[k]
                  - f_529 * kh_145[k]
                  + f_529 * kh_169[k]
                  + f_530 * kh_174[k]
                  - f_531 * kh_176[k]
                  + f_529 * kh_183[k]
                  - f_531 * kh_185[k]
                  + f_532 * kh_187[k]
                  - f_526 * kh_316[k]
                  - f_527 * kh_321[k]
                  + f_528 * kh_323[k]
                  - f_526 * kh_330[k]
                  + f_528 * kh_332[k]
                  - f_529 * kh_334[k]
                  + f_530 * kh_358[k]
                  + f_533 * kh_363[k]
                  - f_534 * kh_365[k]
                  + f_530 * kh_372[k]
                  - f_534 * kh_374[k]
                  + f_535 * kh_376[k]
                  - f_530 * kh_400[k]
                  - f_533 * kh_405[k]
                  + f_534 * kh_407[k]
                  - f_530 * kh_414[k]
                  + f_534 * kh_416[k]
                  - f_535 * kh_418[k]
                  - f_522 * kh_589[k]
                  - f_523 * kh_594[k]
                  + f_524 * kh_596[k]
                  - f_522 * kh_603[k]
                  + f_524 * kh_605[k]
                  - f_525 * kh_607[k]
                  + f_529 * kh_631[k]
                  + f_530 * kh_636[k]
                  - f_531 * kh_638[k]
                  + f_529 * kh_645[k]
                  - f_531 * kh_647[k]
                  + f_532 * kh_649[k]
                  - f_530 * kh_673[k]
                  - f_533 * kh_678[k]
                  + f_534 * kh_680[k]
                  - f_530 * kh_687[k]
                  + f_534 * kh_689[k]
                  - f_535 * kh_691[k]
                  + f_536 * kh_715[k]
                  + f_537 * kh_720[k]
                  - f_538 * kh_722[k]
                  + f_536 * kh_729[k]
                  - f_538 * kh_731[k]
                  + f_539 * kh_733[k];
    }

#pragma omp simd aligned(kh_23, kh_28, kh_30, kh_37, kh_39, kh_41, kh_128, kh_133, kh_135, \
                         kh_142, kh_144, kh_146, kh_170, kh_175, kh_177, kh_184, kh_186, \
                         kh_188, kh_317, kh_322, kh_324, kh_331, kh_333, kh_335, kh_359, \
                         kh_364, kh_366, kh_373, kh_375, kh_377, kh_401, kh_406, kh_408, \
                         kh_415, kh_417, kh_419, kh_590, kh_595, kh_597, kh_604, kh_606, \
                         kh_608, kh_632, kh_637, kh_639, kh_646, kh_648, kh_650, kh_674, \
                         kh_679, kh_681, kh_688, kh_690, kh_692, kh_716, kh_721, kh_723, \
                         kh_730, kh_732, kh_734 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_71[k] = -f_540 * kh_23[k]
                  - f_541 * kh_28[k]
                  + f_542 * kh_30[k]
                  - f_540 * kh_37[k]
                  + f_542 * kh_39[k]
                  - f_543 * kh_41[k]
                  - f_544 * kh_128[k]
                  - f_545 * kh_133[k]
                  + f_546 * kh_135[k]
                  - f_544 * kh_142[k]
                  + f_546 * kh_144[k]
                  - f_547 * kh_146[k]
                  + f_548 * kh_170[k]
                  + f_549 * kh_175[k]
                  - f_550 * kh_177[k]
                  + f_548 * kh_184[k]
                  - f_550 * kh_186[k]
                  + f_551 * kh_188[k]
                  - f_544 * kh_317[k]
                  - f_545 * kh_322[k]
                  + f_546 * kh_324[k]
                  - f_544 * kh_331[k]
                  + f_546 * kh_333[k]
                  - f_547 * kh_335[k]
                  + f_549 * kh_359[k]
                  + f_552 * kh_364[k]
                  - f_553 * kh_366[k]
                  + f_549 * kh_373[k]
                  - f_553 * kh_375[k]
                  + f_554 * kh_377[k]
                  - f_549 * kh_401[k]
                  - f_552 * kh_406[k]
                  + f_553 * kh_408[k]
                  - f_549 * kh_415[k]
                  + f_553 * kh_417[k]
                  - f_554 * kh_419[k]
                  - f_540 * kh_590[k]
                  - f_541 * kh_595[k]
                  + f_542 * kh_597[k]
                  - f_540 * kh_604[k]
                  + f_542 * kh_606[k]
                  - f_543 * kh_608[k]
                  + f_548 * kh_632[k]
                  + f_549 * kh_637[k]
                  - f_550 * kh_639[k]
                  + f_548 * kh_646[k]
                  - f_550 * kh_648[k]
                  + f_551 * kh_650[k]
                  - f_549 * kh_674[k]
                  - f_552 * kh_679[k]
                  + f_553 * kh_681[k]
                  - f_549 * kh_688[k]
                  + f_553 * kh_690[k]
                  - f_554 * kh_692[k]
                  + f_551 * kh_716[k]
                  + f_554 * kh_721[k]
                  - f_555 * kh_723[k]
                  + f_551 * kh_730[k]
                  - f_555 * kh_732[k]
                  + f_556 * kh_734[k];
    }

#pragma omp simd aligned(kh_21, kh_24, kh_26, kh_31, kh_33, kh_35, kh_126, kh_129, kh_131, \
                         kh_136, kh_138, kh_140, kh_168, kh_171, kh_173, kh_178, kh_180, \
                         kh_182, kh_315, kh_318, kh_320, kh_325, kh_327, kh_329, kh_357, \
                         kh_360, kh_362, kh_367, kh_369, kh_371, kh_399, kh_402, kh_404, \
                         kh_409, kh_411, kh_413, kh_588, kh_591, kh_593, kh_598, kh_600, \
                         kh_602, kh_630, kh_633, kh_635, kh_640, kh_642, kh_644, kh_672, \
                         kh_675, kh_677, kh_682, kh_684, kh_686, kh_714, kh_717, kh_719, \
                         kh_724, kh_726, kh_728 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_72[k] = -f_522 * kh_21[k]
                  - f_523 * kh_24[k]
                  + f_524 * kh_26[k]
                  - f_522 * kh_31[k]
                  + f_524 * kh_33[k]
                  - f_525 * kh_35[k]
                  - f_526 * kh_126[k]
                  - f_527 * kh_129[k]
                  + f_528 * kh_131[k]
                  - f_526 * kh_136[k]
                  + f_528 * kh_138[k]
                  - f_529 * kh_140[k]
                  + f_529 * kh_168[k]
                  + f_530 * kh_171[k]
                  - f_531 * kh_173[k]
                  + f_529 * kh_178[k]
                  - f_531 * kh_180[k]
                  + f_532 * kh_182[k]
                  - f_526 * kh_315[k]
                  - f_527 * kh_318[k]
                  + f_528 * kh_320[k]
                  - f_526 * kh_325[k]
                  + f_528 * kh_327[k]
                  - f_529 * kh_329[k]
                  + f_530 * kh_357[k]
                  + f_533 * kh_360[k]
                  - f_534 * kh_362[k]
                  + f_530 * kh_367[k]
                  - f_534 * kh_369[k]
                  + f_535 * kh_371[k]
                  - f_530 * kh_399[k]
                  - f_533 * kh_402[k]
                  + f_534 * kh_404[k]
                  - f_530 * kh_409[k]
                  + f_534 * kh_411[k]
                  - f_535 * kh_413[k]
                  - f_522 * kh_588[k]
                  - f_523 * kh_591[k]
                  + f_524 * kh_593[k]
                  - f_522 * kh_598[k]
                  + f_524 * kh_600[k]
                  - f_525 * kh_602[k]
                  + f_529 * kh_630[k]
                  + f_530 * kh_633[k]
                  - f_531 * kh_635[k]
                  + f_529 * kh_640[k]
                  - f_531 * kh_642[k]
                  + f_532 * kh_644[k]
                  - f_530 * kh_672[k]
                  - f_533 * kh_675[k]
                  + f_534 * kh_677[k]
                  - f_530 * kh_682[k]
                  + f_534 * kh_684[k]
                  - f_535 * kh_686[k]
                  + f_536 * kh_714[k]
                  + f_537 * kh_717[k]
                  - f_538 * kh_719[k]
                  + f_536 * kh_724[k]
                  - f_538 * kh_726[k]
                  + f_539 * kh_728[k];
    }

#pragma omp simd aligned(kh_23, kh_30, kh_37, kh_39, kh_128, kh_135, kh_142, kh_144, kh_170, \
                         kh_177, kh_184, kh_186, kh_317, kh_324, kh_331, kh_333, kh_359, \
                         kh_366, kh_373, kh_375, kh_401, kh_408, kh_415, kh_417, kh_590, \
                         kh_597, kh_604, kh_606, kh_632, kh_639, kh_646, kh_648, kh_674, \
                         kh_681, kh_688, kh_690, kh_716, kh_723, kh_730, \
                         kh_732 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_73[k] = f_557 * kh_23[k]
                  - f_517 * kh_30[k]
                  - f_557 * kh_37[k]
                  + f_517 * kh_39[k]
                  + f_430 * kh_128[k]
                  - f_429 * kh_135[k]
                  - f_430 * kh_142[k]
                  + f_429 * kh_144[k]
                  - f_414 * kh_170[k]
                  + f_432 * kh_177[k]
                  + f_414 * kh_184[k]
                  - f_432 * kh_186[k]
                  + f_430 * kh_317[k]
                  - f_429 * kh_324[k]
                  - f_430 * kh_331[k]
                  + f_429 * kh_333[k]
                  - f_432 * kh_359[k]
                  + f_303 * kh_366[k]
                  + f_432 * kh_373[k]
                  - f_303 * kh_375[k]
                  + f_432 * kh_401[k]
                  - f_303 * kh_408[k]
                  - f_432 * kh_415[k]
                  + f_303 * kh_417[k]
                  + f_557 * kh_590[k]
                  - f_517 * kh_597[k]
                  - f_557 * kh_604[k]
                  + f_517 * kh_606[k]
                  - f_414 * kh_632[k]
                  + f_432 * kh_639[k]
                  + f_414 * kh_646[k]
                  - f_432 * kh_648[k]
                  + f_432 * kh_674[k]
                  - f_303 * kh_681[k]
                  - f_432 * kh_688[k]
                  + f_303 * kh_690[k]
                  - f_558 * kh_716[k]
                  + f_520 * kh_723[k]
                  + f_558 * kh_730[k]
                  - f_520 * kh_732[k];
    }

#pragma omp simd aligned(kh_21, kh_24, kh_26, kh_31, kh_33, kh_126, kh_129, kh_131, kh_136, \
                         kh_138, kh_168, kh_171, kh_173, kh_178, kh_180, kh_315, kh_318, \
                         kh_320, kh_325, kh_327, kh_357, kh_360, kh_362, kh_367, kh_369, \
                         kh_399, kh_402, kh_404, kh_409, kh_411, kh_588, kh_591, kh_593, \
                         kh_598, kh_600, kh_630, kh_633, kh_635, kh_640, kh_642, kh_672, \
                         kh_675, kh_677, kh_682, kh_684, kh_714, kh_717, kh_719, kh_724, \
                         kh_726 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_74[k] = f_504 * kh_21[k]
                  - f_502 * kh_24[k]
                  - f_505 * kh_26[k]
                  - f_501 * kh_31[k]
                  + f_503 * kh_33[k]
                  + f_501 * kh_126[k]
                  - f_507 * kh_129[k]
                  - f_503 * kh_131[k]
                  - f_506 * kh_136[k]
                  + f_478 * kh_138[k]
                  - f_503 * kh_168[k]
                  + f_508 * kh_171[k]
                  + f_509 * kh_173[k]
                  + f_478 * kh_178[k]
                  - f_444 * kh_180[k]
                  + f_501 * kh_315[k]
                  - f_507 * kh_318[k]
                  - f_503 * kh_320[k]
                  - f_506 * kh_325[k]
                  + f_478 * kh_327[k]
                  - f_508 * kh_357[k]
                  + f_510 * kh_360[k]
                  + f_479 * kh_362[k]
                  + f_442 * kh_367[k]
                  - f_511 * kh_369[k]
                  + f_508 * kh_399[k]
                  - f_510 * kh_402[k]
                  - f_479 * kh_404[k]
                  - f_442 * kh_409[k]
                  + f_511 * kh_411[k]
                  + f_504 * kh_588[k]
                  - f_502 * kh_591[k]
                  - f_505 * kh_593[k]
                  - f_501 * kh_598[k]
                  + f_503 * kh_600[k]
                  - f_503 * kh_630[k]
                  + f_508 * kh_633[k]
                  + f_509 * kh_635[k]
                  + f_478 * kh_640[k]
                  - f_444 * kh_642[k]
                  + f_508 * kh_672[k]
                  - f_510 * kh_675[k]
                  - f_479 * kh_677[k]
                  - f_442 * kh_682[k]
                  + f_511 * kh_684[k]
                  - f_515 * kh_714[k]
                  + f_513 * kh_717[k]
                  + f_516 * kh_719[k]
                  + f_512 * kh_724[k]
                  - f_514 * kh_726[k];
    }

#pragma omp simd aligned(kh_23, kh_28, kh_37, kh_128, kh_133, kh_142, kh_170, kh_175, kh_184, \
                         kh_317, kh_322, kh_331, kh_359, kh_364, kh_373, kh_401, kh_406, \
                         kh_415, kh_590, kh_595, kh_604, kh_632, kh_637, kh_646, kh_674, \
                         kh_679, kh_688, kh_716, kh_721, kh_730 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_75[k] = -f_559 * kh_23[k]
                  + f_399 * kh_28[k]
                  - f_559 * kh_37[k]
                  - f_560 * kh_128[k]
                  + f_561 * kh_133[k]
                  - f_560 * kh_142[k]
                  + f_335 * kh_170[k]
                  - f_336 * kh_175[k]
                  + f_335 * kh_184[k]
                  - f_560 * kh_317[k]
                  + f_561 * kh_322[k]
                  - f_560 * kh_331[k]
                  + f_343 * kh_359[k]
                  - f_337 * kh_364[k]
                  + f_343 * kh_373[k]
                  - f_343 * kh_401[k]
                  + f_337 * kh_406[k]
                  - f_343 * kh_415[k]
                  - f_559 * kh_590[k]
                  + f_399 * kh_595[k]
                  - f_559 * kh_604[k]
                  + f_335 * kh_632[k]
                  - f_336 * kh_637[k]
                  + f_335 * kh_646[k]
                  - f_343 * kh_674[k]
                  + f_337 * kh_679[k]
                  - f_343 * kh_688[k]
                  + f_562 * kh_716[k]
                  - f_563 * kh_721[k]
                  + f_562 * kh_730[k];
    }

#pragma omp simd aligned(kh_21, kh_24, kh_31, kh_126, kh_129, kh_136, kh_168, kh_171, kh_178, \
                         kh_315, kh_318, kh_325, kh_357, kh_360, kh_367, kh_399, kh_402, \
                         kh_409, kh_588, kh_591, kh_598, kh_630, kh_633, kh_640, kh_672, \
                         kh_675, kh_682, kh_714, kh_717, kh_724 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_76[k] = -f_487 * kh_21[k]
                  + f_486 * kh_24[k]
                  - f_485 * kh_31[k]
                  - f_490 * kh_126[k]
                  + f_489 * kh_129[k]
                  - f_488 * kh_136[k]
                  + f_493 * kh_168[k]
                  - f_492 * kh_171[k]
                  + f_491 * kh_178[k]
                  - f_490 * kh_315[k]
                  + f_489 * kh_318[k]
                  - f_488 * kh_325[k]
                  + f_495 * kh_357[k]
                  - f_494 * kh_360[k]
                  + f_492 * kh_367[k]
                  - f_495 * kh_399[k]
                  + f_494 * kh_402[k]
                  - f_492 * kh_409[k]
                  - f_487 * kh_588[k]
                  + f_486 * kh_591[k]
                  - f_485 * kh_598[k]
                  + f_493 * kh_630[k]
                  - f_492 * kh_633[k]
                  + f_491 * kh_640[k]
                  - f_495 * kh_672[k]
                  + f_494 * kh_675[k]
                  - f_492 * kh_682[k]
                  + f_498 * kh_714[k]
                  - f_497 * kh_717[k]
                  + f_496 * kh_724[k];
    }

#pragma omp simd aligned(kh_43, kh_48, kh_57, kh_148, kh_153, kh_162, kh_190, kh_195, kh_204, \
                         kh_337, kh_342, kh_351, kh_379, kh_384, kh_393, kh_421, kh_426, \
                         kh_435, kh_610, kh_615, kh_624, kh_652, kh_657, kh_666, kh_694, \
                         kh_699, kh_708, kh_736, kh_741, kh_750 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_77[k] = -f_564 * kh_43[k]
                  + f_565 * kh_48[k]
                  - f_566 * kh_57[k]
                  - f_567 * kh_148[k]
                  + f_568 * kh_153[k]
                  - f_569 * kh_162[k]
                  + f_568 * kh_190[k]
                  - f_570 * kh_195[k]
                  + f_571 * kh_204[k]
                  - f_567 * kh_337[k]
                  + f_568 * kh_342[k]
                  - f_569 * kh_351[k]
                  + f_570 * kh_379[k]
                  - f_572 * kh_384[k]
                  + f_573 * kh_393[k]
                  - f_574 * kh_421[k]
                  + f_575 * kh_426[k]
                  - f_576 * kh_435[k]
                  - f_564 * kh_610[k]
                  + f_565 * kh_615[k]
                  - f_566 * kh_624[k]
                  + f_568 * kh_652[k]
                  - f_570 * kh_657[k]
                  + f_571 * kh_666[k]
                  - f_574 * kh_694[k]
                  + f_575 * kh_699[k]
                  - f_576 * kh_708[k]
                  + f_577 * kh_736[k]
                  - f_578 * kh_741[k]
                  + f_579 * kh_750[k];
    }

#pragma omp simd aligned(kh_46, kh_53, kh_151, kh_158, kh_193, kh_200, kh_340, kh_347, kh_382, \
                         kh_389, kh_424, kh_431, kh_613, kh_620, kh_655, kh_662, kh_697, \
                         kh_704, kh_739, kh_746 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_78[k] = -f_580 * kh_46[k]
                  + f_580 * kh_53[k]
                  - f_581 * kh_151[k]
                  + f_581 * kh_158[k]
                  + f_582 * kh_193[k]
                  - f_582 * kh_200[k]
                  - f_581 * kh_340[k]
                  + f_581 * kh_347[k]
                  + f_583 * kh_382[k]
                  - f_583 * kh_389[k]
                  - f_584 * kh_424[k]
                  + f_584 * kh_431[k]
                  - f_580 * kh_613[k]
                  + f_580 * kh_620[k]
                  + f_582 * kh_655[k]
                  - f_582 * kh_662[k]
                  - f_584 * kh_697[k]
                  + f_584 * kh_704[k]
                  + f_585 * kh_739[k]
                  - f_585 * kh_746[k];
    }

#pragma omp simd aligned(kh_43, kh_48, kh_50, kh_57, kh_59, kh_148, kh_153, kh_155, kh_162, \
                         kh_164, kh_190, kh_195, kh_197, kh_204, kh_206, kh_337, kh_342, \
                         kh_344, kh_351, kh_353, kh_379, kh_384, kh_386, kh_393, kh_395, \
                         kh_421, kh_426, kh_428, kh_435, kh_437, kh_610, kh_615, kh_617, \
                         kh_624, kh_626, kh_652, kh_657, kh_659, kh_666, kh_668, kh_694, \
                         kh_699, kh_701, kh_708, kh_710, kh_736, kh_741, kh_743, kh_750, \
                         kh_752 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_79[k] = f_586 * kh_43[k]
                  + f_587 * kh_48[k]
                  - f_588 * kh_50[k]
                  - f_589 * kh_57[k]
                  + f_590 * kh_59[k]
                  + f_591 * kh_148[k]
                  + f_592 * kh_153[k]
                  - f_593 * kh_155[k]
                  - f_586 * kh_162[k]
                  + f_588 * kh_164[k]
                  - f_594 * kh_190[k]
                  - f_595 * kh_195[k]
                  + f_596 * kh_197[k]
                  + f_592 * kh_204[k]
                  - f_597 * kh_206[k]
                  + f_591 * kh_337[k]
                  + f_592 * kh_342[k]
                  - f_593 * kh_344[k]
                  - f_586 * kh_351[k]
                  + f_588 * kh_353[k]
                  - f_598 * kh_379[k]
                  - f_588 * kh_384[k]
                  + f_599 * kh_386[k]
                  + f_595 * kh_393[k]
                  - f_600 * kh_395[k]
                  + f_601 * kh_421[k]
                  + f_602 * kh_426[k]
                  - f_603 * kh_428[k]
                  - f_604 * kh_435[k]
                  + f_605 * kh_437[k]
                  + f_586 * kh_610[k]
                  + f_587 * kh_615[k]
                  - f_588 * kh_617[k]
                  - f_589 * kh_624[k]
                  + f_590 * kh_626[k]
                  - f_594 * kh_652[k]
                  - f_595 * kh_657[k]
                  + f_596 * kh_659[k]
                  + f_592 * kh_666[k]
                  - f_597 * kh_668[k]
                  + f_601 * kh_694[k]
                  + f_602 * kh_699[k]
                  - f_603 * kh_701[k]
                  - f_604 * kh_708[k]
                  + f_605 * kh_710[k]
                  - f_606 * kh_736[k]
                  - f_607 * kh_741[k]
                  + f_608 * kh_743[k]
                  + f_609 * kh_750[k]
                  - f_610 * kh_752[k];
    }

#pragma omp simd aligned(kh_46, kh_53, kh_55, kh_151, kh_158, kh_160, kh_193, kh_200, kh_202, \
                         kh_340, kh_347, kh_349, kh_382, kh_389, kh_391, kh_424, kh_431, \
                         kh_433, kh_613, kh_620, kh_622, kh_655, kh_662, kh_664, kh_697, \
                         kh_704, kh_706, kh_739, kh_746, kh_748 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_80[k] = f_611 * kh_46[k]
                  + f_611 * kh_53[k]
                  - f_612 * kh_55[k]
                  + f_613 * kh_151[k]
                  + f_613 * kh_158[k]
                  - f_614 * kh_160[k]
                  - f_614 * kh_193[k]
                  - f_614 * kh_200[k]
                  + f_615 * kh_202[k]
                  + f_613 * kh_340[k]
                  + f_613 * kh_347[k]
                  - f_614 * kh_349[k]
                  - f_615 * kh_382[k]
                  - f_615 * kh_389[k]
                  + f_616 * kh_391[k]
                  + f_617 * kh_424[k]
                  + f_617 * kh_431[k]
                  - f_618 * kh_433[k]
                  + f_611 * kh_613[k]
                  + f_611 * kh_620[k]
                  - f_612 * kh_622[k]
                  - f_614 * kh_655[k]
                  - f_614 * kh_662[k]
                  + f_615 * kh_664[k]
                  + f_617 * kh_697[k]
                  + f_617 * kh_704[k]
                  - f_618 * kh_706[k]
                  - f_537 * kh_739[k]
                  - f_537 * kh_746[k]
                  + f_619 * kh_748[k];
    }

#pragma omp simd aligned(kh_43, kh_48, kh_50, kh_57, kh_59, kh_61, kh_148, kh_153, kh_155, \
                         kh_162, kh_164, kh_166, kh_190, kh_195, kh_197, kh_204, kh_206, \
                         kh_208, kh_337, kh_342, kh_344, kh_351, kh_353, kh_355, kh_379, \
                         kh_384, kh_386, kh_393, kh_395, kh_397, kh_421, kh_426, kh_428, \
                         kh_435, kh_437, kh_439, kh_610, kh_615, kh_617, kh_624, kh_626, \
                         kh_628, kh_652, kh_657, kh_659, kh_666, kh_668, kh_670, kh_694, \
                         kh_699, kh_701, kh_708, kh_710, kh_712, kh_736, kh_741, kh_743, \
                         kh_750, kh_752, kh_754 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_81[k] = -f_557 * kh_43[k]
                  - f_517 * kh_48[k]
                  + f_431 * kh_50[k]
                  - f_557 * kh_57[k]
                  + f_431 * kh_59[k]
                  - f_620 * kh_61[k]
                  - f_430 * kh_148[k]
                  - f_429 * kh_153[k]
                  + f_411 * kh_155[k]
                  - f_430 * kh_162[k]
                  + f_411 * kh_164[k]
                  - f_414 * kh_166[k]
                  + f_429 * kh_190[k]
                  + f_431 * kh_195[k]
                  - f_302 * kh_197[k]
                  + f_429 * kh_204[k]
                  - f_302 * kh_206[k]
                  + f_432 * kh_208[k]
                  - f_430 * kh_337[k]
                  - f_429 * kh_342[k]
                  + f_411 * kh_344[k]
                  - f_430 * kh_351[k]
                  + f_411 * kh_353[k]
                  - f_414 * kh_355[k]
                  + f_431 * kh_379[k]
                  + f_414 * kh_384[k]
                  - f_300 * kh_386[k]
                  + f_431 * kh_393[k]
                  - f_300 * kh_395[k]
                  + f_303 * kh_397[k]
                  - f_621 * kh_421[k]
                  - f_440 * kh_426[k]
                  + f_622 * kh_428[k]
                  - f_621 * kh_435[k]
                  + f_622 * kh_437[k]
                  - f_623 * kh_439[k]
                  - f_557 * kh_610[k]
                  - f_517 * kh_615[k]
                  + f_431 * kh_617[k]
                  - f_557 * kh_624[k]
                  + f_431 * kh_626[k]
                  - f_620 * kh_628[k]
                  + f_429 * kh_652[k]
                  + f_431 * kh_657[k]
                  - f_302 * kh_659[k]
                  + f_429 * kh_666[k]
                  - f_302 * kh_668[k]
                  + f_432 * kh_670[k]
                  - f_621 * kh_694[k]
                  - f_440 * kh_699[k]
                  + f_622 * kh_701[k]
                  - f_621 * kh_708[k]
                  + f_622 * kh_710[k]
                  - f_623 * kh_712[k]
                  + f_624 * kh_736[k]
                  + f_625 * kh_741[k]
                  - f_626 * kh_743[k]
                  + f_624 * kh_750[k]
                  - f_626 * kh_752[k]
                  + f_627 * kh_754[k];
    }

#pragma omp simd aligned(kh_44, kh_49, kh_51, kh_58, kh_60, kh_62, kh_149, kh_154, kh_156, \
                         kh_163, kh_165, kh_167, kh_191, kh_196, kh_198, kh_205, kh_207, \
                         kh_209, kh_338, kh_343, kh_345, kh_352, kh_354, kh_356, kh_380, \
                         kh_385, kh_387, kh_394, kh_396, kh_398, kh_422, kh_427, kh_429, \
                         kh_436, kh_438, kh_440, kh_611, kh_616, kh_618, kh_625, kh_627, \
                         kh_629, kh_653, kh_658, kh_660, kh_667, kh_669, kh_671, kh_695, \
                         kh_700, kh_702, kh_709, kh_711, kh_713, kh_737, kh_742, kh_744, \
                         kh_751, kh_753, kh_755 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_82[k] = -4.1015625 * kh_44[k]
                  - 8.203125 * kh_49[k]
                  + 10.9375 * kh_51[k]
                  - 4.1015625 * kh_58[k]
                  + 10.9375 * kh_60[k]
                  - 2.1875 * kh_62[k]
                  - 12.3046875 * kh_149[k]
                  - 24.609375 * kh_154[k]
                  + 32.8125 * kh_156[k]
                  - 12.3046875 * kh_163[k]
                  + 32.8125 * kh_165[k]
                  - 6.5625 * kh_167[k]
                  + 24.609375 * kh_191[k]
                  + 49.21875 * kh_196[k]
                  - 65.625 * kh_198[k]
                  + 24.609375 * kh_205[k]
                  - 65.625 * kh_207[k]
                  + 13.125 * kh_209[k]
                  - 12.3046875 * kh_338[k]
                  - 24.609375 * kh_343[k]
                  + 32.8125 * kh_345[k]
                  - 12.3046875 * kh_352[k]
                  + 32.8125 * kh_354[k]
                  - 6.5625 * kh_356[k]
                  + 49.21875 * kh_380[k]
                  + 98.4375 * kh_385[k]
                  - 131.25 * kh_387[k]
                  + 49.21875 * kh_394[k]
                  - 131.25 * kh_396[k]
                  + 26.25 * kh_398[k]
                  - 19.6875 * kh_422[k]
                  - 39.375 * kh_427[k]
                  + 52.5 * kh_429[k]
                  - 19.6875 * kh_436[k]
                  + 52.5 * kh_438[k]
                  - 10.5 * kh_440[k]
                  - 4.1015625 * kh_611[k]
                  - 8.203125 * kh_616[k]
                  + 10.9375 * kh_618[k]
                  - 4.1015625 * kh_625[k]
                  + 10.9375 * kh_627[k]
                  - 2.1875 * kh_629[k]
                  + 24.609375 * kh_653[k]
                  + 49.21875 * kh_658[k]
                  - 65.625 * kh_660[k]
                  + 24.609375 * kh_667[k]
                  - 65.625 * kh_669[k]
                  + 13.125 * kh_671[k]
                  - 19.6875 * kh_695[k]
                  - 39.375 * kh_700[k]
                  + 52.5 * kh_702[k]
                  - 19.6875 * kh_709[k]
                  + 52.5 * kh_711[k]
                  - 10.5 * kh_713[k]
                  + 1.875 * kh_737[k]
                  + 3.75 * kh_742[k]
                  - 5.0 * kh_744[k]
                  + 1.875 * kh_751[k]
                  - 5.0 * kh_753[k]
                  + kh_755[k];
    }

#pragma omp simd aligned(kh_42, kh_45, kh_47, kh_52, kh_54, kh_56, kh_147, kh_150, kh_152, \
                         kh_157, kh_159, kh_161, kh_189, kh_192, kh_194, kh_199, kh_201, \
                         kh_203, kh_336, kh_339, kh_341, kh_346, kh_348, kh_350, kh_378, \
                         kh_381, kh_383, kh_388, kh_390, kh_392, kh_420, kh_423, kh_425, \
                         kh_430, kh_432, kh_434, kh_609, kh_612, kh_614, kh_619, kh_621, \
                         kh_623, kh_651, kh_654, kh_656, kh_661, kh_663, kh_665, kh_693, \
                         kh_696, kh_698, kh_703, kh_705, kh_707, kh_735, kh_738, kh_740, \
                         kh_745, kh_747, kh_749 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_83[k] = -f_557 * kh_42[k]
                  - f_517 * kh_45[k]
                  + f_431 * kh_47[k]
                  - f_557 * kh_52[k]
                  + f_431 * kh_54[k]
                  - f_620 * kh_56[k]
                  - f_430 * kh_147[k]
                  - f_429 * kh_150[k]
                  + f_411 * kh_152[k]
                  - f_430 * kh_157[k]
                  + f_411 * kh_159[k]
                  - f_414 * kh_161[k]
                  + f_429 * kh_189[k]
                  + f_431 * kh_192[k]
                  - f_302 * kh_194[k]
                  + f_429 * kh_199[k]
                  - f_302 * kh_201[k]
                  + f_432 * kh_203[k]
                  - f_430 * kh_336[k]
                  - f_429 * kh_339[k]
                  + f_411 * kh_341[k]
                  - f_430 * kh_346[k]
                  + f_411 * kh_348[k]
                  - f_414 * kh_350[k]
                  + f_431 * kh_378[k]
                  + f_414 * kh_381[k]
                  - f_300 * kh_383[k]
                  + f_431 * kh_388[k]
                  - f_300 * kh_390[k]
                  + f_303 * kh_392[k]
                  - f_621 * kh_420[k]
                  - f_440 * kh_423[k]
                  + f_622 * kh_425[k]
                  - f_621 * kh_430[k]
                  + f_622 * kh_432[k]
                  - f_623 * kh_434[k]
                  - f_557 * kh_609[k]
                  - f_517 * kh_612[k]
                  + f_431 * kh_614[k]
                  - f_557 * kh_619[k]
                  + f_431 * kh_621[k]
                  - f_620 * kh_623[k]
                  + f_429 * kh_651[k]
                  + f_431 * kh_654[k]
                  - f_302 * kh_656[k]
                  + f_429 * kh_661[k]
                  - f_302 * kh_663[k]
                  + f_432 * kh_665[k]
                  - f_621 * kh_693[k]
                  - f_440 * kh_696[k]
                  + f_622 * kh_698[k]
                  - f_621 * kh_703[k]
                  + f_622 * kh_705[k]
                  - f_623 * kh_707[k]
                  + f_624 * kh_735[k]
                  + f_625 * kh_738[k]
                  - f_626 * kh_740[k]
                  + f_624 * kh_745[k]
                  - f_626 * kh_747[k]
                  + f_627 * kh_749[k];
    }

#pragma omp simd aligned(kh_44, kh_51, kh_58, kh_60, kh_149, kh_156, kh_163, kh_165, kh_191, \
                         kh_198, kh_205, kh_207, kh_338, kh_345, kh_352, kh_354, kh_380, \
                         kh_387, kh_394, kh_396, kh_422, kh_429, kh_436, kh_438, kh_611, \
                         kh_618, kh_625, kh_627, kh_653, kh_660, kh_667, kh_669, kh_695, \
                         kh_702, kh_709, kh_711, kh_737, kh_744, kh_751, \
                         kh_753 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_84[k] = f_628 * kh_44[k]
                  - f_611 * kh_51[k]
                  - f_628 * kh_58[k]
                  + f_611 * kh_60[k]
                  + f_629 * kh_149[k]
                  - f_613 * kh_156[k]
                  - f_629 * kh_163[k]
                  + f_613 * kh_165[k]
                  - f_613 * kh_191[k]
                  + f_614 * kh_198[k]
                  + f_613 * kh_205[k]
                  - f_614 * kh_207[k]
                  + f_629 * kh_338[k]
                  - f_613 * kh_345[k]
                  - f_629 * kh_352[k]
                  + f_613 * kh_354[k]
                  - f_614 * kh_380[k]
                  + f_615 * kh_387[k]
                  + f_614 * kh_394[k]
                  - f_615 * kh_396[k]
                  + f_630 * kh_422[k]
                  - f_617 * kh_429[k]
                  - f_630 * kh_436[k]
                  + f_617 * kh_438[k]
                  + f_628 * kh_611[k]
                  - f_611 * kh_618[k]
                  - f_628 * kh_625[k]
                  + f_611 * kh_627[k]
                  - f_613 * kh_653[k]
                  + f_614 * kh_660[k]
                  + f_613 * kh_667[k]
                  - f_614 * kh_669[k]
                  + f_630 * kh_695[k]
                  - f_617 * kh_702[k]
                  - f_630 * kh_709[k]
                  + f_617 * kh_711[k]
                  - f_536 * kh_737[k]
                  + f_537 * kh_744[k]
                  + f_536 * kh_751[k]
                  - f_537 * kh_753[k];
    }

#pragma omp simd aligned(kh_42, kh_45, kh_47, kh_52, kh_54, kh_147, kh_150, kh_152, kh_157, \
                         kh_159, kh_189, kh_192, kh_194, kh_199, kh_201, kh_336, kh_339, \
                         kh_341, kh_346, kh_348, kh_378, kh_381, kh_383, kh_388, kh_390, \
                         kh_420, kh_423, kh_425, kh_430, kh_432, kh_609, kh_612, kh_614, \
                         kh_619, kh_621, kh_651, kh_654, kh_656, kh_661, kh_663, kh_693, \
                         kh_696, kh_698, kh_703, kh_705, kh_735, kh_738, kh_740, kh_745, \
                         kh_747 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_85[k] = f_589 * kh_42[k]
                  - f_587 * kh_45[k]
                  - f_590 * kh_47[k]
                  - f_586 * kh_52[k]
                  + f_588 * kh_54[k]
                  + f_586 * kh_147[k]
                  - f_592 * kh_150[k]
                  - f_588 * kh_152[k]
                  - f_591 * kh_157[k]
                  + f_593 * kh_159[k]
                  - f_592 * kh_189[k]
                  + f_595 * kh_192[k]
                  + f_597 * kh_194[k]
                  + f_594 * kh_199[k]
                  - f_596 * kh_201[k]
                  + f_586 * kh_336[k]
                  - f_592 * kh_339[k]
                  - f_588 * kh_341[k]
                  - f_591 * kh_346[k]
                  + f_593 * kh_348[k]
                  - f_595 * kh_378[k]
                  + f_588 * kh_381[k]
                  + f_600 * kh_383[k]
                  + f_598 * kh_388[k]
                  - f_599 * kh_390[k]
                  + f_604 * kh_420[k]
                  - f_602 * kh_423[k]
                  - f_605 * kh_425[k]
                  - f_601 * kh_430[k]
                  + f_603 * kh_432[k]
                  + f_589 * kh_609[k]
                  - f_587 * kh_612[k]
                  - f_590 * kh_614[k]
                  - f_586 * kh_619[k]
                  + f_588 * kh_621[k]
                  - f_592 * kh_651[k]
                  + f_595 * kh_654[k]
                  + f_597 * kh_656[k]
                  + f_594 * kh_661[k]
                  - f_596 * kh_663[k]
                  + f_604 * kh_693[k]
                  - f_602 * kh_696[k]
                  - f_605 * kh_698[k]
                  - f_601 * kh_703[k]
                  + f_603 * kh_705[k]
                  - f_609 * kh_735[k]
                  + f_607 * kh_738[k]
                  + f_610 * kh_740[k]
                  + f_606 * kh_745[k]
                  - f_608 * kh_747[k];
    }

#pragma omp simd aligned(kh_44, kh_49, kh_58, kh_149, kh_154, kh_163, kh_191, kh_196, kh_205, \
                         kh_338, kh_343, kh_352, kh_380, kh_385, kh_394, kh_422, kh_427, \
                         kh_436, kh_611, kh_616, kh_625, kh_653, kh_658, kh_667, kh_695, \
                         kh_700, kh_709, kh_737, kh_742, kh_751 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_86[k] = -f_631 * kh_44[k]
                  + f_632 * kh_49[k]
                  - f_631 * kh_58[k]
                  - f_633 * kh_149[k]
                  + f_634 * kh_154[k]
                  - f_633 * kh_163[k]
                  + f_632 * kh_191[k]
                  - f_635 * kh_196[k]
                  + f_632 * kh_205[k]
                  - f_633 * kh_338[k]
                  + f_634 * kh_343[k]
                  - f_633 * kh_352[k]
                  + f_581 * kh_380[k]
                  - f_636 * kh_385[k]
                  + f_581 * kh_394[k]
                  - f_637 * kh_422[k]
                  + f_638 * kh_427[k]
                  - f_637 * kh_436[k]
                  - f_631 * kh_611[k]
                  + f_632 * kh_616[k]
                  - f_631 * kh_625[k]
                  + f_632 * kh_653[k]
                  - f_635 * kh_658[k]
                  + f_632 * kh_667[k]
                  - f_637 * kh_695[k]
                  + f_638 * kh_700[k]
                  - f_637 * kh_709[k]
                  + f_639 * kh_737[k]
                  - f_640 * kh_742[k]
                  + f_639 * kh_751[k];
    }

#pragma omp simd aligned(kh_42, kh_45, kh_52, kh_147, kh_150, kh_157, kh_189, kh_192, kh_199, \
                         kh_336, kh_339, kh_346, kh_378, kh_381, kh_388, kh_420, kh_423, \
                         kh_430, kh_609, kh_612, kh_619, kh_651, kh_654, kh_661, kh_693, \
                         kh_696, kh_703, kh_735, kh_738, kh_745 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_87[k] = -f_566 * kh_42[k]
                  + f_565 * kh_45[k]
                  - f_564 * kh_52[k]
                  - f_569 * kh_147[k]
                  + f_568 * kh_150[k]
                  - f_567 * kh_157[k]
                  + f_571 * kh_189[k]
                  - f_570 * kh_192[k]
                  + f_568 * kh_199[k]
                  - f_569 * kh_336[k]
                  + f_568 * kh_339[k]
                  - f_567 * kh_346[k]
                  + f_573 * kh_378[k]
                  - f_572 * kh_381[k]
                  + f_570 * kh_388[k]
                  - f_576 * kh_420[k]
                  + f_575 * kh_423[k]
                  - f_574 * kh_430[k]
                  - f_566 * kh_609[k]
                  + f_565 * kh_612[k]
                  - f_564 * kh_619[k]
                  + f_571 * kh_651[k]
                  - f_570 * kh_654[k]
                  + f_568 * kh_661[k]
                  - f_576 * kh_693[k]
                  + f_575 * kh_696[k]
                  - f_574 * kh_703[k]
                  + f_579 * kh_735[k]
                  - f_578 * kh_738[k]
                  + f_577 * kh_745[k];
    }

#pragma omp simd aligned(kh_1, kh_6, kh_15, kh_64, kh_69, kh_78, kh_106, kh_111, kh_120, \
                         kh_211, kh_216, kh_225, kh_253, kh_258, kh_267, kh_295, kh_300, \
                         kh_309, kh_442, kh_447, kh_456, kh_484, kh_489, kh_498, kh_526, \
                         kh_531, kh_540, kh_568, kh_573, kh_582 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_88[k] = -f_485 * kh_1[k]
                  + f_486 * kh_6[k]
                  - f_487 * kh_15[k]
                  - f_488 * kh_64[k]
                  + f_489 * kh_69[k]
                  - f_490 * kh_78[k]
                  + f_491 * kh_106[k]
                  - f_492 * kh_111[k]
                  + f_493 * kh_120[k]
                  - f_488 * kh_211[k]
                  + f_489 * kh_216[k]
                  - f_490 * kh_225[k]
                  + f_492 * kh_253[k]
                  - f_494 * kh_258[k]
                  + f_495 * kh_267[k]
                  - f_492 * kh_295[k]
                  + f_494 * kh_300[k]
                  - f_495 * kh_309[k]
                  - f_485 * kh_442[k]
                  + f_486 * kh_447[k]
                  - f_487 * kh_456[k]
                  + f_491 * kh_484[k]
                  - f_492 * kh_489[k]
                  + f_493 * kh_498[k]
                  - f_492 * kh_526[k]
                  + f_494 * kh_531[k]
                  - f_495 * kh_540[k]
                  + f_496 * kh_568[k]
                  - f_497 * kh_573[k]
                  + f_498 * kh_582[k];
    }

#pragma omp simd aligned(kh_4, kh_11, kh_67, kh_74, kh_109, kh_116, kh_214, kh_221, kh_256, \
                         kh_263, kh_298, kh_305, kh_445, kh_452, kh_487, kh_494, kh_529, \
                         kh_536, kh_571, kh_578 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_89[k] = -f_499 * kh_4[k]
                  + f_499 * kh_11[k]
                  - f_334 * kh_67[k]
                  + f_334 * kh_74[k]
                  + f_340 * kh_109[k]
                  - f_340 * kh_116[k]
                  - f_334 * kh_214[k]
                  + f_334 * kh_221[k]
                  + f_341 * kh_256[k]
                  - f_341 * kh_263[k]
                  - f_341 * kh_298[k]
                  + f_341 * kh_305[k]
                  - f_499 * kh_445[k]
                  + f_499 * kh_452[k]
                  + f_340 * kh_487[k]
                  - f_340 * kh_494[k]
                  - f_341 * kh_529[k]
                  + f_341 * kh_536[k]
                  + f_500 * kh_571[k]
                  - f_500 * kh_578[k];
    }

#pragma omp simd aligned(kh_1, kh_6, kh_8, kh_15, kh_17, kh_64, kh_69, kh_71, kh_78, kh_80, \
                         kh_106, kh_111, kh_113, kh_120, kh_122, kh_211, kh_216, kh_218, \
                         kh_225, kh_227, kh_253, kh_258, kh_260, kh_267, kh_269, kh_295, \
                         kh_300, kh_302, kh_309, kh_311, kh_442, kh_447, kh_449, kh_456, \
                         kh_458, kh_484, kh_489, kh_491, kh_498, kh_500, kh_526, kh_531, \
                         kh_533, kh_540, kh_542, kh_568, kh_573, kh_575, kh_582, \
                         kh_584 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_90[k] = f_501 * kh_1[k]
                  + f_502 * kh_6[k]
                  - f_503 * kh_8[k]
                  - f_504 * kh_15[k]
                  + f_505 * kh_17[k]
                  + f_506 * kh_64[k]
                  + f_507 * kh_69[k]
                  - f_478 * kh_71[k]
                  - f_501 * kh_78[k]
                  + f_503 * kh_80[k]
                  - f_478 * kh_106[k]
                  - f_508 * kh_111[k]
                  + f_444 * kh_113[k]
                  + f_503 * kh_120[k]
                  - f_509 * kh_122[k]
                  + f_506 * kh_211[k]
                  + f_507 * kh_216[k]
                  - f_478 * kh_218[k]
                  - f_501 * kh_225[k]
                  + f_503 * kh_227[k]
                  - f_442 * kh_253[k]
                  - f_510 * kh_258[k]
                  + f_511 * kh_260[k]
                  + f_508 * kh_267[k]
                  - f_479 * kh_269[k]
                  + f_442 * kh_295[k]
                  + f_510 * kh_300[k]
                  - f_511 * kh_302[k]
                  - f_508 * kh_309[k]
                  + f_479 * kh_311[k]
                  + f_501 * kh_442[k]
                  + f_502 * kh_447[k]
                  - f_503 * kh_449[k]
                  - f_504 * kh_456[k]
                  + f_505 * kh_458[k]
                  - f_478 * kh_484[k]
                  - f_508 * kh_489[k]
                  + f_444 * kh_491[k]
                  + f_503 * kh_498[k]
                  - f_509 * kh_500[k]
                  + f_442 * kh_526[k]
                  + f_510 * kh_531[k]
                  - f_511 * kh_533[k]
                  - f_508 * kh_540[k]
                  + f_479 * kh_542[k]
                  - f_512 * kh_568[k]
                  - f_513 * kh_573[k]
                  + f_514 * kh_575[k]
                  + f_515 * kh_582[k]
                  - f_516 * kh_584[k];
    }

#pragma omp simd aligned(kh_4, kh_11, kh_13, kh_67, kh_74, kh_76, kh_109, kh_116, kh_118, \
                         kh_214, kh_221, kh_223, kh_256, kh_263, kh_265, kh_298, kh_305, \
                         kh_307, kh_445, kh_452, kh_454, kh_487, kh_494, kh_496, kh_529, \
                         kh_536, kh_538, kh_571, kh_578, kh_580 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_91[k] = f_517 * kh_4[k]
                  + f_517 * kh_11[k]
                  - f_518 * kh_13[k]
                  + f_429 * kh_67[k]
                  + f_429 * kh_74[k]
                  - f_431 * kh_76[k]
                  - f_432 * kh_109[k]
                  - f_432 * kh_116[k]
                  + f_303 * kh_118[k]
                  + f_429 * kh_214[k]
                  + f_429 * kh_221[k]
                  - f_431 * kh_223[k]
                  - f_303 * kh_256[k]
                  - f_303 * kh_263[k]
                  + f_519 * kh_265[k]
                  + f_303 * kh_298[k]
                  + f_303 * kh_305[k]
                  - f_519 * kh_307[k]
                  + f_517 * kh_445[k]
                  + f_517 * kh_452[k]
                  - f_518 * kh_454[k]
                  - f_432 * kh_487[k]
                  - f_432 * kh_494[k]
                  + f_303 * kh_496[k]
                  + f_303 * kh_529[k]
                  + f_303 * kh_536[k]
                  - f_519 * kh_538[k]
                  - f_520 * kh_571[k]
                  - f_520 * kh_578[k]
                  + f_521 * kh_580[k];
    }

#pragma omp simd aligned(kh_1, kh_6, kh_8, kh_15, kh_17, kh_19, kh_64, kh_69, kh_71, kh_78, \
                         kh_80, kh_82, kh_106, kh_111, kh_113, kh_120, kh_122, kh_124, kh_211, \
                         kh_216, kh_218, kh_225, kh_227, kh_229, kh_253, kh_258, kh_260, \
                         kh_267, kh_269, kh_271, kh_295, kh_300, kh_302, kh_309, kh_311, \
                         kh_313, kh_442, kh_447, kh_449, kh_456, kh_458, kh_460, kh_484, \
                         kh_489, kh_491, kh_498, kh_500, kh_502, kh_526, kh_531, kh_533, \
                         kh_540, kh_542, kh_544, kh_568, kh_573, kh_575, kh_582, kh_584, \
                         kh_586 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_92[k] = -f_522 * kh_1[k]
                  - f_523 * kh_6[k]
                  + f_524 * kh_8[k]
                  - f_522 * kh_15[k]
                  + f_524 * kh_17[k]
                  - f_525 * kh_19[k]
                  - f_526 * kh_64[k]
                  - f_527 * kh_69[k]
                  + f_528 * kh_71[k]
                  - f_526 * kh_78[k]
                  + f_528 * kh_80[k]
                  - f_529 * kh_82[k]
                  + f_529 * kh_106[k]
                  + f_530 * kh_111[k]
                  - f_531 * kh_113[k]
                  + f_529 * kh_120[k]
                  - f_531 * kh_122[k]
                  + f_532 * kh_124[k]
                  - f_526 * kh_211[k]
                  - f_527 * kh_216[k]
                  + f_528 * kh_218[k]
                  - f_526 * kh_225[k]
                  + f_528 * kh_227[k]
                  - f_529 * kh_229[k]
                  + f_530 * kh_253[k]
                  + f_533 * kh_258[k]
                  - f_534 * kh_260[k]
                  + f_530 * kh_267[k]
                  - f_534 * kh_269[k]
                  + f_535 * kh_271[k]
                  - f_530 * kh_295[k]
                  - f_533 * kh_300[k]
                  + f_534 * kh_302[k]
                  - f_530 * kh_309[k]
                  + f_534 * kh_311[k]
                  - f_535 * kh_313[k]
                  - f_522 * kh_442[k]
                  - f_523 * kh_447[k]
                  + f_524 * kh_449[k]
                  - f_522 * kh_456[k]
                  + f_524 * kh_458[k]
                  - f_525 * kh_460[k]
                  + f_529 * kh_484[k]
                  + f_530 * kh_489[k]
                  - f_531 * kh_491[k]
                  + f_529 * kh_498[k]
                  - f_531 * kh_500[k]
                  + f_532 * kh_502[k]
                  - f_530 * kh_526[k]
                  - f_533 * kh_531[k]
                  + f_534 * kh_533[k]
                  - f_530 * kh_540[k]
                  + f_534 * kh_542[k]
                  - f_535 * kh_544[k]
                  + f_536 * kh_568[k]
                  + f_537 * kh_573[k]
                  - f_538 * kh_575[k]
                  + f_536 * kh_582[k]
                  - f_538 * kh_584[k]
                  + f_539 * kh_586[k];
    }

#pragma omp simd aligned(kh_2, kh_7, kh_9, kh_16, kh_18, kh_20, kh_65, kh_70, kh_72, kh_79, \
                         kh_81, kh_83, kh_107, kh_112, kh_114, kh_121, kh_123, kh_125, kh_212, \
                         kh_217, kh_219, kh_226, kh_228, kh_230, kh_254, kh_259, kh_261, \
                         kh_268, kh_270, kh_272, kh_296, kh_301, kh_303, kh_310, kh_312, \
                         kh_314, kh_443, kh_448, kh_450, kh_457, kh_459, kh_461, kh_485, \
                         kh_490, kh_492, kh_499, kh_501, kh_503, kh_527, kh_532, kh_534, \
                         kh_541, kh_543, kh_545, kh_569, kh_574, kh_576, kh_583, kh_585, \
                         kh_587 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_93[k] = -f_540 * kh_2[k]
                  - f_541 * kh_7[k]
                  + f_542 * kh_9[k]
                  - f_540 * kh_16[k]
                  + f_542 * kh_18[k]
                  - f_543 * kh_20[k]
                  - f_544 * kh_65[k]
                  - f_545 * kh_70[k]
                  + f_546 * kh_72[k]
                  - f_544 * kh_79[k]
                  + f_546 * kh_81[k]
                  - f_547 * kh_83[k]
                  + f_548 * kh_107[k]
                  + f_549 * kh_112[k]
                  - f_550 * kh_114[k]
                  + f_548 * kh_121[k]
                  - f_550 * kh_123[k]
                  + f_551 * kh_125[k]
                  - f_544 * kh_212[k]
                  - f_545 * kh_217[k]
                  + f_546 * kh_219[k]
                  - f_544 * kh_226[k]
                  + f_546 * kh_228[k]
                  - f_547 * kh_230[k]
                  + f_549 * kh_254[k]
                  + f_552 * kh_259[k]
                  - f_553 * kh_261[k]
                  + f_549 * kh_268[k]
                  - f_553 * kh_270[k]
                  + f_554 * kh_272[k]
                  - f_549 * kh_296[k]
                  - f_552 * kh_301[k]
                  + f_553 * kh_303[k]
                  - f_549 * kh_310[k]
                  + f_553 * kh_312[k]
                  - f_554 * kh_314[k]
                  - f_540 * kh_443[k]
                  - f_541 * kh_448[k]
                  + f_542 * kh_450[k]
                  - f_540 * kh_457[k]
                  + f_542 * kh_459[k]
                  - f_543 * kh_461[k]
                  + f_548 * kh_485[k]
                  + f_549 * kh_490[k]
                  - f_550 * kh_492[k]
                  + f_548 * kh_499[k]
                  - f_550 * kh_501[k]
                  + f_551 * kh_503[k]
                  - f_549 * kh_527[k]
                  - f_552 * kh_532[k]
                  + f_553 * kh_534[k]
                  - f_549 * kh_541[k]
                  + f_553 * kh_543[k]
                  - f_554 * kh_545[k]
                  + f_551 * kh_569[k]
                  + f_554 * kh_574[k]
                  - f_555 * kh_576[k]
                  + f_551 * kh_583[k]
                  - f_555 * kh_585[k]
                  + f_556 * kh_587[k];
    }

#pragma omp simd aligned(kh_0, kh_3, kh_5, kh_10, kh_12, kh_14, kh_63, kh_66, kh_68, kh_73, \
                         kh_75, kh_77, kh_105, kh_108, kh_110, kh_115, kh_117, kh_119, kh_210, \
                         kh_213, kh_215, kh_220, kh_222, kh_224, kh_252, kh_255, kh_257, \
                         kh_262, kh_264, kh_266, kh_294, kh_297, kh_299, kh_304, kh_306, \
                         kh_308, kh_441, kh_444, kh_446, kh_451, kh_453, kh_455, kh_483, \
                         kh_486, kh_488, kh_493, kh_495, kh_497, kh_525, kh_528, kh_530, \
                         kh_535, kh_537, kh_539, kh_567, kh_570, kh_572, kh_577, kh_579, \
                         kh_581 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_94[k] = -f_522 * kh_0[k]
                  - f_523 * kh_3[k]
                  + f_524 * kh_5[k]
                  - f_522 * kh_10[k]
                  + f_524 * kh_12[k]
                  - f_525 * kh_14[k]
                  - f_526 * kh_63[k]
                  - f_527 * kh_66[k]
                  + f_528 * kh_68[k]
                  - f_526 * kh_73[k]
                  + f_528 * kh_75[k]
                  - f_529 * kh_77[k]
                  + f_529 * kh_105[k]
                  + f_530 * kh_108[k]
                  - f_531 * kh_110[k]
                  + f_529 * kh_115[k]
                  - f_531 * kh_117[k]
                  + f_532 * kh_119[k]
                  - f_526 * kh_210[k]
                  - f_527 * kh_213[k]
                  + f_528 * kh_215[k]
                  - f_526 * kh_220[k]
                  + f_528 * kh_222[k]
                  - f_529 * kh_224[k]
                  + f_530 * kh_252[k]
                  + f_533 * kh_255[k]
                  - f_534 * kh_257[k]
                  + f_530 * kh_262[k]
                  - f_534 * kh_264[k]
                  + f_535 * kh_266[k]
                  - f_530 * kh_294[k]
                  - f_533 * kh_297[k]
                  + f_534 * kh_299[k]
                  - f_530 * kh_304[k]
                  + f_534 * kh_306[k]
                  - f_535 * kh_308[k]
                  - f_522 * kh_441[k]
                  - f_523 * kh_444[k]
                  + f_524 * kh_446[k]
                  - f_522 * kh_451[k]
                  + f_524 * kh_453[k]
                  - f_525 * kh_455[k]
                  + f_529 * kh_483[k]
                  + f_530 * kh_486[k]
                  - f_531 * kh_488[k]
                  + f_529 * kh_493[k]
                  - f_531 * kh_495[k]
                  + f_532 * kh_497[k]
                  - f_530 * kh_525[k]
                  - f_533 * kh_528[k]
                  + f_534 * kh_530[k]
                  - f_530 * kh_535[k]
                  + f_534 * kh_537[k]
                  - f_535 * kh_539[k]
                  + f_536 * kh_567[k]
                  + f_537 * kh_570[k]
                  - f_538 * kh_572[k]
                  + f_536 * kh_577[k]
                  - f_538 * kh_579[k]
                  + f_539 * kh_581[k];
    }

#pragma omp simd aligned(kh_2, kh_9, kh_16, kh_18, kh_65, kh_72, kh_79, kh_81, kh_107, kh_114, \
                         kh_121, kh_123, kh_212, kh_219, kh_226, kh_228, kh_254, kh_261, \
                         kh_268, kh_270, kh_296, kh_303, kh_310, kh_312, kh_443, kh_450, \
                         kh_457, kh_459, kh_485, kh_492, kh_499, kh_501, kh_527, kh_534, \
                         kh_541, kh_543, kh_569, kh_576, kh_583, \
                         kh_585 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_95[k] = f_557 * kh_2[k]
                  - f_517 * kh_9[k]
                  - f_557 * kh_16[k]
                  + f_517 * kh_18[k]
                  + f_430 * kh_65[k]
                  - f_429 * kh_72[k]
                  - f_430 * kh_79[k]
                  + f_429 * kh_81[k]
                  - f_414 * kh_107[k]
                  + f_432 * kh_114[k]
                  + f_414 * kh_121[k]
                  - f_432 * kh_123[k]
                  + f_430 * kh_212[k]
                  - f_429 * kh_219[k]
                  - f_430 * kh_226[k]
                  + f_429 * kh_228[k]
                  - f_432 * kh_254[k]
                  + f_303 * kh_261[k]
                  + f_432 * kh_268[k]
                  - f_303 * kh_270[k]
                  + f_432 * kh_296[k]
                  - f_303 * kh_303[k]
                  - f_432 * kh_310[k]
                  + f_303 * kh_312[k]
                  + f_557 * kh_443[k]
                  - f_517 * kh_450[k]
                  - f_557 * kh_457[k]
                  + f_517 * kh_459[k]
                  - f_414 * kh_485[k]
                  + f_432 * kh_492[k]
                  + f_414 * kh_499[k]
                  - f_432 * kh_501[k]
                  + f_432 * kh_527[k]
                  - f_303 * kh_534[k]
                  - f_432 * kh_541[k]
                  + f_303 * kh_543[k]
                  - f_558 * kh_569[k]
                  + f_520 * kh_576[k]
                  + f_558 * kh_583[k]
                  - f_520 * kh_585[k];
    }

#pragma omp simd aligned(kh_0, kh_3, kh_5, kh_10, kh_12, kh_63, kh_66, kh_68, kh_73, kh_75, \
                         kh_105, kh_108, kh_110, kh_115, kh_117, kh_210, kh_213, kh_215, \
                         kh_220, kh_222, kh_252, kh_255, kh_257, kh_262, kh_264, kh_294, \
                         kh_297, kh_299, kh_304, kh_306, kh_441, kh_444, kh_446, kh_451, \
                         kh_453, kh_483, kh_486, kh_488, kh_493, kh_495, kh_525, kh_528, \
                         kh_530, kh_535, kh_537, kh_567, kh_570, kh_572, kh_577, \
                         kh_579 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_96[k] = f_504 * kh_0[k]
                  - f_502 * kh_3[k]
                  - f_505 * kh_5[k]
                  - f_501 * kh_10[k]
                  + f_503 * kh_12[k]
                  + f_501 * kh_63[k]
                  - f_507 * kh_66[k]
                  - f_503 * kh_68[k]
                  - f_506 * kh_73[k]
                  + f_478 * kh_75[k]
                  - f_503 * kh_105[k]
                  + f_508 * kh_108[k]
                  + f_509 * kh_110[k]
                  + f_478 * kh_115[k]
                  - f_444 * kh_117[k]
                  + f_501 * kh_210[k]
                  - f_507 * kh_213[k]
                  - f_503 * kh_215[k]
                  - f_506 * kh_220[k]
                  + f_478 * kh_222[k]
                  - f_508 * kh_252[k]
                  + f_510 * kh_255[k]
                  + f_479 * kh_257[k]
                  + f_442 * kh_262[k]
                  - f_511 * kh_264[k]
                  + f_508 * kh_294[k]
                  - f_510 * kh_297[k]
                  - f_479 * kh_299[k]
                  - f_442 * kh_304[k]
                  + f_511 * kh_306[k]
                  + f_504 * kh_441[k]
                  - f_502 * kh_444[k]
                  - f_505 * kh_446[k]
                  - f_501 * kh_451[k]
                  + f_503 * kh_453[k]
                  - f_503 * kh_483[k]
                  + f_508 * kh_486[k]
                  + f_509 * kh_488[k]
                  + f_478 * kh_493[k]
                  - f_444 * kh_495[k]
                  + f_508 * kh_525[k]
                  - f_510 * kh_528[k]
                  - f_479 * kh_530[k]
                  - f_442 * kh_535[k]
                  + f_511 * kh_537[k]
                  - f_515 * kh_567[k]
                  + f_513 * kh_570[k]
                  + f_516 * kh_572[k]
                  + f_512 * kh_577[k]
                  - f_514 * kh_579[k];
    }

#pragma omp simd aligned(kh_2, kh_7, kh_16, kh_65, kh_70, kh_79, kh_107, kh_112, kh_121, \
                         kh_212, kh_217, kh_226, kh_254, kh_259, kh_268, kh_296, kh_301, \
                         kh_310, kh_443, kh_448, kh_457, kh_485, kh_490, kh_499, kh_527, \
                         kh_532, kh_541, kh_569, kh_574, kh_583 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_97[k] = -f_559 * kh_2[k]
                  + f_399 * kh_7[k]
                  - f_559 * kh_16[k]
                  - f_560 * kh_65[k]
                  + f_561 * kh_70[k]
                  - f_560 * kh_79[k]
                  + f_335 * kh_107[k]
                  - f_336 * kh_112[k]
                  + f_335 * kh_121[k]
                  - f_560 * kh_212[k]
                  + f_561 * kh_217[k]
                  - f_560 * kh_226[k]
                  + f_343 * kh_254[k]
                  - f_337 * kh_259[k]
                  + f_343 * kh_268[k]
                  - f_343 * kh_296[k]
                  + f_337 * kh_301[k]
                  - f_343 * kh_310[k]
                  - f_559 * kh_443[k]
                  + f_399 * kh_448[k]
                  - f_559 * kh_457[k]
                  + f_335 * kh_485[k]
                  - f_336 * kh_490[k]
                  + f_335 * kh_499[k]
                  - f_343 * kh_527[k]
                  + f_337 * kh_532[k]
                  - f_343 * kh_541[k]
                  + f_562 * kh_569[k]
                  - f_563 * kh_574[k]
                  + f_562 * kh_583[k];
    }

#pragma omp simd aligned(kh_0, kh_3, kh_10, kh_63, kh_66, kh_73, kh_105, kh_108, kh_115, \
                         kh_210, kh_213, kh_220, kh_252, kh_255, kh_262, kh_294, kh_297, \
                         kh_304, kh_441, kh_444, kh_451, kh_483, kh_486, kh_493, kh_525, \
                         kh_528, kh_535, kh_567, kh_570, kh_577 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_98[k] = -f_487 * kh_0[k]
                  + f_486 * kh_3[k]
                  - f_485 * kh_10[k]
                  - f_490 * kh_63[k]
                  + f_489 * kh_66[k]
                  - f_488 * kh_73[k]
                  + f_493 * kh_105[k]
                  - f_492 * kh_108[k]
                  + f_491 * kh_115[k]
                  - f_490 * kh_210[k]
                  + f_489 * kh_213[k]
                  - f_488 * kh_220[k]
                  + f_495 * kh_252[k]
                  - f_494 * kh_255[k]
                  + f_492 * kh_262[k]
                  - f_495 * kh_294[k]
                  + f_494 * kh_297[k]
                  - f_492 * kh_304[k]
                  - f_487 * kh_441[k]
                  + f_486 * kh_444[k]
                  - f_485 * kh_451[k]
                  + f_493 * kh_483[k]
                  - f_492 * kh_486[k]
                  + f_491 * kh_493[k]
                  - f_495 * kh_525[k]
                  + f_494 * kh_528[k]
                  - f_492 * kh_535[k]
                  + f_498 * kh_567[k]
                  - f_497 * kh_570[k]
                  + f_496 * kh_577[k];
    }

#pragma omp simd aligned(kh_43, kh_48, kh_57, kh_148, kh_153, kh_162, kh_190, kh_195, kh_204, \
                         kh_337, kh_342, kh_351, kh_421, kh_426, kh_435, kh_610, kh_615, \
                         kh_624, kh_652, kh_657, kh_666, kh_694, kh_699, \
                         kh_708 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_99[k] = f_641 * kh_43[k]
                  - f_415 * kh_48[k]
                  + f_642 * kh_57[k]
                  + f_641 * kh_148[k]
                  - f_415 * kh_153[k]
                  + f_642 * kh_162[k]
                  - f_643 * kh_190[k]
                  + f_420 * kh_195[k]
                  - f_644 * kh_204[k]
                  - f_641 * kh_337[k]
                  + f_415 * kh_342[k]
                  - f_642 * kh_351[k]
                  + f_645 * kh_421[k]
                  - f_423 * kh_426[k]
                  + f_646 * kh_435[k]
                  - f_641 * kh_610[k]
                  + f_415 * kh_615[k]
                  - f_642 * kh_624[k]
                  + f_643 * kh_652[k]
                  - f_420 * kh_657[k]
                  + f_644 * kh_666[k]
                  - f_645 * kh_694[k]
                  + f_423 * kh_699[k]
                  - f_646 * kh_708[k];
    }

#pragma omp simd aligned(kh_46, kh_53, kh_151, kh_158, kh_193, kh_200, kh_340, kh_347, kh_424, \
                         kh_431, kh_613, kh_620, kh_655, kh_662, kh_697, \
                         kh_704 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_100[k] = f_311 * kh_46[k]
                   - f_311 * kh_53[k]
                   + f_311 * kh_151[k]
                   - f_311 * kh_158[k]
                   - f_327 * kh_193[k]
                   + f_327 * kh_200[k]
                   - f_311 * kh_340[k]
                   + f_311 * kh_347[k]
                   + f_647 * kh_424[k]
                   - f_647 * kh_431[k]
                   - f_311 * kh_613[k]
                   + f_311 * kh_620[k]
                   + f_327 * kh_655[k]
                   - f_327 * kh_662[k]
                   - f_647 * kh_697[k]
                   + f_647 * kh_704[k];
    }

#pragma omp simd aligned(kh_43, kh_48, kh_50, kh_57, kh_59, kh_148, kh_153, kh_155, kh_162, \
                         kh_164, kh_190, kh_195, kh_197, kh_204, kh_206, kh_337, kh_342, \
                         kh_344, kh_351, kh_353, kh_421, kh_426, kh_428, kh_435, kh_437, \
                         kh_610, kh_615, kh_617, kh_624, kh_626, kh_652, kh_657, kh_659, \
                         kh_666, kh_668, kh_694, kh_699, kh_701, kh_708, \
                         kh_710 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_101[k] = -f_405 * kh_43[k]
                   - f_430 * kh_48[k]
                   + f_411 * kh_50[k]
                   + f_648 * kh_57[k]
                   - f_431 * kh_59[k]
                   - f_405 * kh_148[k]
                   - f_430 * kh_153[k]
                   + f_411 * kh_155[k]
                   + f_648 * kh_162[k]
                   - f_431 * kh_164[k]
                   + f_414 * kh_190[k]
                   + f_435 * kh_195[k]
                   - f_519 * kh_197[k]
                   - f_620 * kh_204[k]
                   + f_649 * kh_206[k]
                   + f_405 * kh_337[k]
                   + f_430 * kh_342[k]
                   - f_411 * kh_344[k]
                   - f_648 * kh_351[k]
                   + f_431 * kh_353[k]
                   - f_650 * kh_421[k]
                   - f_440 * kh_426[k]
                   + f_651 * kh_428[k]
                   + f_621 * kh_435[k]
                   - f_623 * kh_437[k]
                   + f_405 * kh_610[k]
                   + f_430 * kh_615[k]
                   - f_411 * kh_617[k]
                   - f_648 * kh_624[k]
                   + f_431 * kh_626[k]
                   - f_414 * kh_652[k]
                   - f_435 * kh_657[k]
                   + f_519 * kh_659[k]
                   + f_620 * kh_666[k]
                   - f_649 * kh_668[k]
                   + f_650 * kh_694[k]
                   + f_440 * kh_699[k]
                   - f_651 * kh_701[k]
                   - f_621 * kh_708[k]
                   + f_623 * kh_710[k];
    }

#pragma omp simd aligned(kh_46, kh_53, kh_55, kh_151, kh_158, kh_160, kh_193, kh_200, kh_202, \
                         kh_340, kh_347, kh_349, kh_424, kh_431, kh_433, kh_613, kh_620, \
                         kh_622, kh_655, kh_662, kh_664, kh_697, kh_704, \
                         kh_706 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_102[k] = -f_478 * kh_46[k]
                   - f_478 * kh_53[k]
                   + f_442 * kh_55[k]
                   - f_478 * kh_151[k]
                   - f_478 * kh_158[k]
                   + f_442 * kh_160[k]
                   + f_479 * kh_193[k]
                   + f_479 * kh_200[k]
                   - f_445 * kh_202[k]
                   + f_478 * kh_340[k]
                   + f_478 * kh_347[k]
                   - f_442 * kh_349[k]
                   - f_480 * kh_424[k]
                   - f_480 * kh_431[k]
                   + f_447 * kh_433[k]
                   + f_478 * kh_613[k]
                   + f_478 * kh_620[k]
                   - f_442 * kh_622[k]
                   - f_479 * kh_655[k]
                   - f_479 * kh_662[k]
                   + f_445 * kh_664[k]
                   + f_480 * kh_697[k]
                   + f_480 * kh_704[k]
                   - f_447 * kh_706[k];
    }

#pragma omp simd aligned(kh_43, kh_48, kh_50, kh_57, kh_59, kh_61, kh_148, kh_153, kh_155, \
                         kh_162, kh_164, kh_166, kh_190, kh_195, kh_197, kh_204, kh_206, \
                         kh_208, kh_337, kh_342, kh_344, kh_351, kh_353, kh_355, kh_421, \
                         kh_426, kh_428, kh_435, kh_437, kh_439, kh_610, kh_615, kh_617, \
                         kh_624, kh_626, kh_628, kh_652, kh_657, kh_659, kh_666, kh_668, \
                         kh_670, kh_694, kh_699, kh_701, kh_708, kh_710, \
                         kh_712 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_103[k] = f_652 * kh_43[k]
                   + f_449 * kh_48[k]
                   - f_653 * kh_50[k]
                   + f_652 * kh_57[k]
                   - f_653 * kh_59[k]
                   + f_453 * kh_61[k]
                   + f_652 * kh_148[k]
                   + f_449 * kh_153[k]
                   - f_653 * kh_155[k]
                   + f_652 * kh_162[k]
                   - f_653 * kh_164[k]
                   + f_453 * kh_166[k]
                   - f_654 * kh_190[k]
                   - f_456 * kh_195[k]
                   + f_655 * kh_197[k]
                   - f_654 * kh_204[k]
                   + f_655 * kh_206[k]
                   - f_656 * kh_208[k]
                   - f_652 * kh_337[k]
                   - f_449 * kh_342[k]
                   + f_653 * kh_344[k]
                   - f_652 * kh_351[k]
                   + f_653 * kh_353[k]
                   - f_453 * kh_355[k]
                   + f_657 * kh_421[k]
                   + f_460 * kh_426[k]
                   - f_658 * kh_428[k]
                   + f_657 * kh_435[k]
                   - f_658 * kh_437[k]
                   + f_659 * kh_439[k]
                   - f_652 * kh_610[k]
                   - f_449 * kh_615[k]
                   + f_653 * kh_617[k]
                   - f_652 * kh_624[k]
                   + f_653 * kh_626[k]
                   - f_453 * kh_628[k]
                   + f_654 * kh_652[k]
                   + f_456 * kh_657[k]
                   - f_655 * kh_659[k]
                   + f_654 * kh_666[k]
                   - f_655 * kh_668[k]
                   + f_656 * kh_670[k]
                   - f_657 * kh_694[k]
                   - f_460 * kh_699[k]
                   + f_658 * kh_701[k]
                   - f_657 * kh_708[k]
                   + f_658 * kh_710[k]
                   - f_659 * kh_712[k];
    }

#pragma omp simd aligned(kh_44, kh_49, kh_51, kh_58, kh_60, kh_62, kh_149, kh_154, kh_156, \
                         kh_163, kh_165, kh_167, kh_191, kh_196, kh_198, kh_205, kh_207, \
                         kh_209, kh_338, kh_343, kh_345, kh_352, kh_354, kh_356, kh_422, \
                         kh_427, kh_429, kh_436, kh_438, kh_440, kh_611, kh_616, kh_618, \
                         kh_625, kh_627, kh_629, kh_653, kh_658, kh_660, kh_667, kh_669, \
                         kh_671, kh_695, kh_700, kh_702, kh_709, kh_711, \
                         kh_713 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_104[k] = f_660 * kh_44[k]
                   + f_464 * kh_49[k]
                   - f_661 * kh_51[k]
                   + f_660 * kh_58[k]
                   - f_661 * kh_60[k]
                   + f_662 * kh_62[k]
                   + f_660 * kh_149[k]
                   + f_464 * kh_154[k]
                   - f_661 * kh_156[k]
                   + f_660 * kh_163[k]
                   - f_661 * kh_165[k]
                   + f_662 * kh_167[k]
                   - f_466 * kh_191[k]
                   - f_469 * kh_196[k]
                   + f_663 * kh_198[k]
                   - f_466 * kh_205[k]
                   + f_663 * kh_207[k]
                   - f_664 * kh_209[k]
                   - f_660 * kh_338[k]
                   - f_464 * kh_343[k]
                   + f_661 * kh_345[k]
                   - f_660 * kh_352[k]
                   + f_661 * kh_354[k]
                   - f_662 * kh_356[k]
                   + f_665 * kh_422[k]
                   + f_474 * kh_427[k]
                   - f_666 * kh_429[k]
                   + f_665 * kh_436[k]
                   - f_666 * kh_438[k]
                   + f_667 * kh_440[k]
                   - f_660 * kh_611[k]
                   - f_464 * kh_616[k]
                   + f_661 * kh_618[k]
                   - f_660 * kh_625[k]
                   + f_661 * kh_627[k]
                   - f_662 * kh_629[k]
                   + f_466 * kh_653[k]
                   + f_469 * kh_658[k]
                   - f_663 * kh_660[k]
                   + f_466 * kh_667[k]
                   - f_663 * kh_669[k]
                   + f_664 * kh_671[k]
                   - f_665 * kh_695[k]
                   - f_474 * kh_700[k]
                   + f_666 * kh_702[k]
                   - f_665 * kh_709[k]
                   + f_666 * kh_711[k]
                   - f_667 * kh_713[k];
    }

#pragma omp simd aligned(kh_42, kh_45, kh_47, kh_52, kh_54, kh_56, kh_147, kh_150, kh_152, \
                         kh_157, kh_159, kh_161, kh_189, kh_192, kh_194, kh_199, kh_201, \
                         kh_203, kh_336, kh_339, kh_341, kh_346, kh_348, kh_350, kh_420, \
                         kh_423, kh_425, kh_430, kh_432, kh_434, kh_609, kh_612, kh_614, \
                         kh_619, kh_621, kh_623, kh_651, kh_654, kh_656, kh_661, kh_663, \
                         kh_665, kh_693, kh_696, kh_698, kh_703, kh_705, \
                         kh_707 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_105[k] = f_652 * kh_42[k]
                   + f_449 * kh_45[k]
                   - f_653 * kh_47[k]
                   + f_652 * kh_52[k]
                   - f_653 * kh_54[k]
                   + f_453 * kh_56[k]
                   + f_652 * kh_147[k]
                   + f_449 * kh_150[k]
                   - f_653 * kh_152[k]
                   + f_652 * kh_157[k]
                   - f_653 * kh_159[k]
                   + f_453 * kh_161[k]
                   - f_654 * kh_189[k]
                   - f_456 * kh_192[k]
                   + f_655 * kh_194[k]
                   - f_654 * kh_199[k]
                   + f_655 * kh_201[k]
                   - f_656 * kh_203[k]
                   - f_652 * kh_336[k]
                   - f_449 * kh_339[k]
                   + f_653 * kh_341[k]
                   - f_652 * kh_346[k]
                   + f_653 * kh_348[k]
                   - f_453 * kh_350[k]
                   + f_657 * kh_420[k]
                   + f_460 * kh_423[k]
                   - f_658 * kh_425[k]
                   + f_657 * kh_430[k]
                   - f_658 * kh_432[k]
                   + f_659 * kh_434[k]
                   - f_652 * kh_609[k]
                   - f_449 * kh_612[k]
                   + f_653 * kh_614[k]
                   - f_652 * kh_619[k]
                   + f_653 * kh_621[k]
                   - f_453 * kh_623[k]
                   + f_654 * kh_651[k]
                   + f_456 * kh_654[k]
                   - f_655 * kh_656[k]
                   + f_654 * kh_661[k]
                   - f_655 * kh_663[k]
                   + f_656 * kh_665[k]
                   - f_657 * kh_693[k]
                   - f_460 * kh_696[k]
                   + f_658 * kh_698[k]
                   - f_657 * kh_703[k]
                   + f_658 * kh_705[k]
                   - f_659 * kh_707[k];
    }

#pragma omp simd aligned(kh_44, kh_51, kh_58, kh_60, kh_149, kh_156, kh_163, kh_165, kh_191, \
                         kh_198, kh_205, kh_207, kh_338, kh_345, kh_352, kh_354, kh_422, \
                         kh_429, kh_436, kh_438, kh_611, kh_618, kh_625, kh_627, kh_653, \
                         kh_660, kh_667, kh_669, kh_695, kh_702, kh_709, \
                         kh_711 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_106[k] = -f_668 * kh_44[k]
                   + f_478 * kh_51[k]
                   + f_668 * kh_58[k]
                   - f_478 * kh_60[k]
                   - f_668 * kh_149[k]
                   + f_478 * kh_156[k]
                   + f_668 * kh_163[k]
                   - f_478 * kh_165[k]
                   + f_509 * kh_191[k]
                   - f_479 * kh_198[k]
                   - f_509 * kh_205[k]
                   + f_479 * kh_207[k]
                   + f_668 * kh_338[k]
                   - f_478 * kh_345[k]
                   - f_668 * kh_352[k]
                   + f_478 * kh_354[k]
                   - f_669 * kh_422[k]
                   + f_480 * kh_429[k]
                   + f_669 * kh_436[k]
                   - f_480 * kh_438[k]
                   + f_668 * kh_611[k]
                   - f_478 * kh_618[k]
                   - f_668 * kh_625[k]
                   + f_478 * kh_627[k]
                   - f_509 * kh_653[k]
                   + f_479 * kh_660[k]
                   + f_509 * kh_667[k]
                   - f_479 * kh_669[k]
                   + f_669 * kh_695[k]
                   - f_480 * kh_702[k]
                   - f_669 * kh_709[k]
                   + f_480 * kh_711[k];
    }

#pragma omp simd aligned(kh_42, kh_45, kh_47, kh_52, kh_54, kh_147, kh_150, kh_152, kh_157, \
                         kh_159, kh_189, kh_192, kh_194, kh_199, kh_201, kh_336, kh_339, \
                         kh_341, kh_346, kh_348, kh_420, kh_423, kh_425, kh_430, kh_432, \
                         kh_609, kh_612, kh_614, kh_619, kh_621, kh_651, kh_654, kh_656, \
                         kh_661, kh_663, kh_693, kh_696, kh_698, kh_703, \
                         kh_705 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_107[k] = -f_648 * kh_42[k]
                   + f_430 * kh_45[k]
                   + f_431 * kh_47[k]
                   + f_405 * kh_52[k]
                   - f_411 * kh_54[k]
                   - f_648 * kh_147[k]
                   + f_430 * kh_150[k]
                   + f_431 * kh_152[k]
                   + f_405 * kh_157[k]
                   - f_411 * kh_159[k]
                   + f_620 * kh_189[k]
                   - f_435 * kh_192[k]
                   - f_649 * kh_194[k]
                   - f_414 * kh_199[k]
                   + f_519 * kh_201[k]
                   + f_648 * kh_336[k]
                   - f_430 * kh_339[k]
                   - f_431 * kh_341[k]
                   - f_405 * kh_346[k]
                   + f_411 * kh_348[k]
                   - f_621 * kh_420[k]
                   + f_440 * kh_423[k]
                   + f_623 * kh_425[k]
                   + f_650 * kh_430[k]
                   - f_651 * kh_432[k]
                   + f_648 * kh_609[k]
                   - f_430 * kh_612[k]
                   - f_431 * kh_614[k]
                   - f_405 * kh_619[k]
                   + f_411 * kh_621[k]
                   - f_620 * kh_651[k]
                   + f_435 * kh_654[k]
                   + f_649 * kh_656[k]
                   + f_414 * kh_661[k]
                   - f_519 * kh_663[k]
                   + f_621 * kh_693[k]
                   - f_440 * kh_696[k]
                   - f_623 * kh_698[k]
                   - f_650 * kh_703[k]
                   + f_651 * kh_705[k];
    }

#pragma omp simd aligned(kh_44, kh_49, kh_58, kh_149, kh_154, kh_163, kh_191, kh_196, kh_205, \
                         kh_338, kh_343, kh_352, kh_422, kh_427, kh_436, kh_611, kh_616, \
                         kh_625, kh_653, kh_658, kh_667, kh_695, kh_700, \
                         kh_709 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_108[k] = f_670 * kh_44[k]
                   - f_314 * kh_49[k]
                   + f_670 * kh_58[k]
                   + f_670 * kh_149[k]
                   - f_314 * kh_154[k]
                   + f_670 * kh_163[k]
                   - f_325 * kh_191[k]
                   + f_322 * kh_196[k]
                   - f_325 * kh_205[k]
                   - f_670 * kh_338[k]
                   + f_314 * kh_343[k]
                   - f_670 * kh_352[k]
                   + f_671 * kh_422[k]
                   - f_672 * kh_427[k]
                   + f_671 * kh_436[k]
                   - f_670 * kh_611[k]
                   + f_314 * kh_616[k]
                   - f_670 * kh_625[k]
                   + f_325 * kh_653[k]
                   - f_322 * kh_658[k]
                   + f_325 * kh_667[k]
                   - f_671 * kh_695[k]
                   + f_672 * kh_700[k]
                   - f_671 * kh_709[k];
    }

#pragma omp simd aligned(kh_42, kh_45, kh_52, kh_147, kh_150, kh_157, kh_189, kh_192, kh_199, \
                         kh_336, kh_339, kh_346, kh_420, kh_423, kh_430, kh_609, kh_612, \
                         kh_619, kh_651, kh_654, kh_661, kh_693, kh_696, \
                         kh_703 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_109[k] = f_642 * kh_42[k]
                   - f_415 * kh_45[k]
                   + f_641 * kh_52[k]
                   + f_642 * kh_147[k]
                   - f_415 * kh_150[k]
                   + f_641 * kh_157[k]
                   - f_644 * kh_189[k]
                   + f_420 * kh_192[k]
                   - f_643 * kh_199[k]
                   - f_642 * kh_336[k]
                   + f_415 * kh_339[k]
                   - f_641 * kh_346[k]
                   + f_646 * kh_420[k]
                   - f_423 * kh_423[k]
                   + f_645 * kh_430[k]
                   - f_642 * kh_609[k]
                   + f_415 * kh_612[k]
                   - f_641 * kh_619[k]
                   + f_644 * kh_651[k]
                   - f_420 * kh_654[k]
                   + f_643 * kh_661[k]
                   - f_646 * kh_693[k]
                   + f_423 * kh_696[k]
                   - f_645 * kh_703[k];
    }

#pragma omp simd aligned(kh_1, kh_6, kh_15, kh_64, kh_69, kh_78, kh_106, kh_111, kh_120, \
                         kh_211, kh_216, kh_225, kh_253, kh_258, kh_267, kh_295, kh_300, \
                         kh_309, kh_442, kh_447, kh_456, kh_484, kh_489, kh_498, kh_526, \
                         kh_531, kh_540 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_110[k] = f_280 * kh_1[k]
                   - f_284 * kh_6[k]
                   + f_285 * kh_15[k]
                   - f_280 * kh_64[k]
                   + f_284 * kh_69[k]
                   - f_285 * kh_78[k]
                   - f_291 * kh_106[k]
                   + f_286 * kh_111[k]
                   - f_292 * kh_120[k]
                   - f_278 * kh_211[k]
                   + f_279 * kh_216[k]
                   - f_280 * kh_225[k]
                   + f_286 * kh_253[k]
                   - f_287 * kh_258[k]
                   + f_288 * kh_267[k]
                   + f_293 * kh_295[k]
                   - f_294 * kh_300[k]
                   + f_295 * kh_309[k]
                   - f_275 * kh_442[k]
                   + f_276 * kh_447[k]
                   - f_277 * kh_456[k]
                   + f_281 * kh_484[k]
                   - f_282 * kh_489[k]
                   + f_283 * kh_498[k]
                   - f_287 * kh_526[k]
                   + f_289 * kh_531[k]
                   - f_290 * kh_540[k];
    }

#pragma omp simd aligned(kh_4, kh_11, kh_67, kh_74, kh_109, kh_116, kh_214, kh_221, kh_256, \
                         kh_263, kh_298, kh_305, kh_445, kh_452, kh_487, kh_494, kh_529, \
                         kh_536 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_111[k] = f_299 * kh_4[k]
                   - f_299 * kh_11[k]
                   - f_299 * kh_67[k]
                   + f_299 * kh_74[k]
                   - f_302 * kh_109[k]
                   + f_302 * kh_116[k]
                   - f_297 * kh_214[k]
                   + f_297 * kh_221[k]
                   + f_300 * kh_256[k]
                   - f_300 * kh_263[k]
                   + f_303 * kh_298[k]
                   - f_303 * kh_305[k]
                   - f_296 * kh_445[k]
                   + f_296 * kh_452[k]
                   + f_298 * kh_487[k]
                   - f_298 * kh_494[k]
                   - f_301 * kh_529[k]
                   + f_301 * kh_536[k];
    }

#pragma omp simd aligned(kh_1, kh_6, kh_8, kh_15, kh_17, kh_64, kh_69, kh_71, kh_78, kh_80, \
                         kh_106, kh_111, kh_113, kh_120, kh_122, kh_211, kh_216, kh_218, \
                         kh_225, kh_227, kh_253, kh_258, kh_260, kh_267, kh_269, kh_295, \
                         kh_300, kh_302, kh_309, kh_311, kh_442, kh_447, kh_449, kh_456, \
                         kh_458, kh_484, kh_489, kh_491, kh_498, kh_500, kh_526, kh_531, \
                         kh_533, kh_540, kh_542 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_112[k] = -f_307 * kh_1[k]
                   - f_318 * kh_6[k]
                   + f_308 * kh_8[k]
                   + f_319 * kh_15[k]
                   - f_320 * kh_17[k]
                   + f_307 * kh_64[k]
                   + f_318 * kh_69[k]
                   - f_308 * kh_71[k]
                   - f_319 * kh_78[k]
                   + f_320 * kh_80[k]
                   + f_316 * kh_106[k]
                   + f_313 * kh_111[k]
                   - f_317 * kh_113[k]
                   - f_328 * kh_120[k]
                   + f_325 * kh_122[k]
                   + f_309 * kh_211[k]
                   + f_310 * kh_216[k]
                   - f_311 * kh_218[k]
                   - f_312 * kh_225[k]
                   + f_313 * kh_227[k]
                   - f_311 * kh_253[k]
                   - f_321 * kh_258[k]
                   + f_322 * kh_260[k]
                   + f_313 * kh_267[k]
                   - f_323 * kh_269[k]
                   - f_321 * kh_295[k]
                   - f_329 * kh_300[k]
                   + f_327 * kh_302[k]
                   + f_330 * kh_309[k]
                   - f_331 * kh_311[k]
                   + f_304 * kh_442[k]
                   + f_305 * kh_447[k]
                   - f_306 * kh_449[k]
                   - f_307 * kh_456[k]
                   + f_308 * kh_458[k]
                   - f_314 * kh_484[k]
                   - f_311 * kh_489[k]
                   + f_315 * kh_491[k]
                   + f_316 * kh_498[k]
                   - f_317 * kh_500[k]
                   + f_324 * kh_526[k]
                   + f_325 * kh_531[k]
                   - f_326 * kh_533[k]
                   - f_321 * kh_540[k]
                   + f_327 * kh_542[k];
    }

#pragma omp simd aligned(kh_4, kh_11, kh_13, kh_67, kh_74, kh_76, kh_109, kh_116, kh_118, \
                         kh_214, kh_221, kh_223, kh_256, kh_263, kh_265, kh_298, kh_305, \
                         kh_307, kh_445, kh_452, kh_454, kh_487, kh_494, kh_496, kh_529, \
                         kh_536, kh_538 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_113[k] = -f_338 * kh_4[k]
                   - f_338 * kh_11[k]
                   + f_339 * kh_13[k]
                   + f_338 * kh_67[k]
                   + f_338 * kh_74[k]
                   - f_339 * kh_76[k]
                   + f_343 * kh_109[k]
                   + f_343 * kh_116[k]
                   - f_340 * kh_118[k]
                   + f_334 * kh_214[k]
                   + f_334 * kh_221[k]
                   - f_335 * kh_223[k]
                   - f_340 * kh_256[k]
                   - f_340 * kh_263[k]
                   + f_341 * kh_265[k]
                   - f_344 * kh_298[k]
                   - f_344 * kh_305[k]
                   + f_345 * kh_307[k]
                   + f_332 * kh_445[k]
                   + f_332 * kh_452[k]
                   - f_333 * kh_454[k]
                   - f_336 * kh_487[k]
                   - f_336 * kh_494[k]
                   + f_337 * kh_496[k]
                   + f_341 * kh_529[k]
                   + f_341 * kh_536[k]
                   - f_342 * kh_538[k];
    }

#pragma omp simd aligned(kh_1, kh_6, kh_8, kh_15, kh_17, kh_19, kh_64, kh_69, kh_71, kh_78, \
                         kh_80, kh_82, kh_106, kh_111, kh_113, kh_120, kh_122, kh_124, kh_211, \
                         kh_216, kh_218, kh_225, kh_227, kh_229, kh_253, kh_258, kh_260, \
                         kh_267, kh_269, kh_271, kh_295, kh_300, kh_302, kh_309, kh_311, \
                         kh_313, kh_442, kh_447, kh_449, kh_456, kh_458, kh_460, kh_484, \
                         kh_489, kh_491, kh_498, kh_500, kh_502, kh_526, kh_531, kh_533, \
                         kh_540, kh_542, kh_544 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_114[k] = f_357 * kh_1[k]
                   + f_358 * kh_6[k]
                   - f_359 * kh_8[k]
                   + f_357 * kh_15[k]
                   - f_359 * kh_17[k]
                   + f_360 * kh_19[k]
                   - f_357 * kh_64[k]
                   - f_358 * kh_69[k]
                   + f_359 * kh_71[k]
                   - f_357 * kh_78[k]
                   + f_359 * kh_80[k]
                   - f_360 * kh_82[k]
                   - f_366 * kh_106[k]
                   - f_353 * kh_111[k]
                   + f_367 * kh_113[k]
                   - f_366 * kh_120[k]
                   + f_367 * kh_122[k]
                   - f_363 * kh_124[k]
                   - f_350 * kh_211[k]
                   - f_351 * kh_216[k]
                   + f_352 * kh_218[k]
                   - f_350 * kh_225[k]
                   + f_352 * kh_227[k]
                   - f_353 * kh_229[k]
                   + f_353 * kh_253[k]
                   + f_361 * kh_258[k]
                   - f_356 * kh_260[k]
                   + f_353 * kh_267[k]
                   - f_356 * kh_269[k]
                   + f_362 * kh_271[k]
                   + f_368 * kh_295[k]
                   + f_369 * kh_300[k]
                   - f_362 * kh_302[k]
                   + f_368 * kh_309[k]
                   - f_362 * kh_311[k]
                   + f_370 * kh_313[k]
                   - f_346 * kh_442[k]
                   - f_347 * kh_447[k]
                   + f_348 * kh_449[k]
                   - f_346 * kh_456[k]
                   + f_348 * kh_458[k]
                   - f_349 * kh_460[k]
                   + f_352 * kh_484[k]
                   + f_354 * kh_489[k]
                   - f_355 * kh_491[k]
                   + f_352 * kh_498[k]
                   - f_355 * kh_500[k]
                   + f_356 * kh_502[k]
                   - f_361 * kh_526[k]
                   - f_363 * kh_531[k]
                   + f_364 * kh_533[k]
                   - f_361 * kh_540[k]
                   + f_364 * kh_542[k]
                   - f_365 * kh_544[k];
    }

#pragma omp simd aligned(kh_2, kh_7, kh_9, kh_16, kh_18, kh_20, kh_65, kh_70, kh_72, kh_79, \
                         kh_81, kh_83, kh_107, kh_112, kh_114, kh_121, kh_123, kh_125, kh_212, \
                         kh_217, kh_219, kh_226, kh_228, kh_230, kh_254, kh_259, kh_261, \
                         kh_268, kh_270, kh_272, kh_296, kh_301, kh_303, kh_310, kh_312, \
                         kh_314, kh_443, kh_448, kh_450, kh_457, kh_459, kh_461, kh_485, \
                         kh_490, kh_492, kh_499, kh_501, kh_503, kh_527, kh_532, kh_534, \
                         kh_541, kh_543, kh_545 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_115[k] = f_383 * kh_2[k]
                   + f_384 * kh_7[k]
                   - f_378 * kh_9[k]
                   + f_383 * kh_16[k]
                   - f_378 * kh_18[k]
                   + f_385 * kh_20[k]
                   - f_383 * kh_65[k]
                   - f_384 * kh_70[k]
                   + f_378 * kh_72[k]
                   - f_383 * kh_79[k]
                   + f_378 * kh_81[k]
                   - f_385 * kh_83[k]
                   - f_392 * kh_107[k]
                   - f_386 * kh_112[k]
                   + f_393 * kh_114[k]
                   - f_392 * kh_121[k]
                   + f_393 * kh_123[k]
                   - f_394 * kh_125[k]
                   - f_375 * kh_212[k]
                   - f_376 * kh_217[k]
                   + f_377 * kh_219[k]
                   - f_375 * kh_226[k]
                   + f_377 * kh_228[k]
                   - f_378 * kh_230[k]
                   + f_386 * kh_254[k]
                   + f_387 * kh_259[k]
                   - f_388 * kh_261[k]
                   + f_386 * kh_268[k]
                   - f_388 * kh_270[k]
                   + f_389 * kh_272[k]
                   + f_395 * kh_296[k]
                   + f_393 * kh_301[k]
                   - f_396 * kh_303[k]
                   + f_395 * kh_310[k]
                   - f_396 * kh_312[k]
                   + f_397 * kh_314[k]
                   - f_371 * kh_443[k]
                   - f_372 * kh_448[k]
                   + f_373 * kh_450[k]
                   - f_371 * kh_457[k]
                   + f_373 * kh_459[k]
                   - f_374 * kh_461[k]
                   + f_379 * kh_485[k]
                   + f_380 * kh_490[k]
                   - f_381 * kh_492[k]
                   + f_379 * kh_499[k]
                   - f_381 * kh_501[k]
                   + f_382 * kh_503[k]
                   - f_387 * kh_527[k]
                   - f_381 * kh_532[k]
                   + f_390 * kh_534[k]
                   - f_387 * kh_541[k]
                   + f_390 * kh_543[k]
                   - f_391 * kh_545[k];
    }

#pragma omp simd aligned(kh_0, kh_3, kh_5, kh_10, kh_12, kh_14, kh_63, kh_66, kh_68, kh_73, \
                         kh_75, kh_77, kh_105, kh_108, kh_110, kh_115, kh_117, kh_119, kh_210, \
                         kh_213, kh_215, kh_220, kh_222, kh_224, kh_252, kh_255, kh_257, \
                         kh_262, kh_264, kh_266, kh_294, kh_297, kh_299, kh_304, kh_306, \
                         kh_308, kh_441, kh_444, kh_446, kh_451, kh_453, kh_455, kh_483, \
                         kh_486, kh_488, kh_493, kh_495, kh_497, kh_525, kh_528, kh_530, \
                         kh_535, kh_537, kh_539 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_116[k] = f_357 * kh_0[k]
                   + f_358 * kh_3[k]
                   - f_359 * kh_5[k]
                   + f_357 * kh_10[k]
                   - f_359 * kh_12[k]
                   + f_360 * kh_14[k]
                   - f_357 * kh_63[k]
                   - f_358 * kh_66[k]
                   + f_359 * kh_68[k]
                   - f_357 * kh_73[k]
                   + f_359 * kh_75[k]
                   - f_360 * kh_77[k]
                   - f_366 * kh_105[k]
                   - f_353 * kh_108[k]
                   + f_367 * kh_110[k]
                   - f_366 * kh_115[k]
                   + f_367 * kh_117[k]
                   - f_363 * kh_119[k]
                   - f_350 * kh_210[k]
                   - f_351 * kh_213[k]
                   + f_352 * kh_215[k]
                   - f_350 * kh_220[k]
                   + f_352 * kh_222[k]
                   - f_353 * kh_224[k]
                   + f_353 * kh_252[k]
                   + f_361 * kh_255[k]
                   - f_356 * kh_257[k]
                   + f_353 * kh_262[k]
                   - f_356 * kh_264[k]
                   + f_362 * kh_266[k]
                   + f_368 * kh_294[k]
                   + f_369 * kh_297[k]
                   - f_362 * kh_299[k]
                   + f_368 * kh_304[k]
                   - f_362 * kh_306[k]
                   + f_370 * kh_308[k]
                   - f_346 * kh_441[k]
                   - f_347 * kh_444[k]
                   + f_348 * kh_446[k]
                   - f_346 * kh_451[k]
                   + f_348 * kh_453[k]
                   - f_349 * kh_455[k]
                   + f_352 * kh_483[k]
                   + f_354 * kh_486[k]
                   - f_355 * kh_488[k]
                   + f_352 * kh_493[k]
                   - f_355 * kh_495[k]
                   + f_356 * kh_497[k]
                   - f_361 * kh_525[k]
                   - f_363 * kh_528[k]
                   + f_364 * kh_530[k]
                   - f_361 * kh_535[k]
                   + f_364 * kh_537[k]
                   - f_365 * kh_539[k];
    }

#pragma omp simd aligned(kh_2, kh_9, kh_16, kh_18, kh_65, kh_72, kh_79, kh_81, kh_107, kh_114, \
                         kh_121, kh_123, kh_212, kh_219, kh_226, kh_228, kh_254, kh_261, \
                         kh_268, kh_270, kh_296, kh_303, kh_310, kh_312, kh_443, kh_450, \
                         kh_457, kh_459, kh_485, kh_492, kh_499, kh_501, kh_527, kh_534, \
                         kh_541, kh_543 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_117[k] = -f_401 * kh_2[k]
                   + f_338 * kh_9[k]
                   + f_401 * kh_16[k]
                   - f_338 * kh_18[k]
                   + f_401 * kh_65[k]
                   - f_338 * kh_72[k]
                   - f_401 * kh_79[k]
                   + f_338 * kh_81[k]
                   + f_335 * kh_107[k]
                   - f_343 * kh_114[k]
                   - f_335 * kh_121[k]
                   + f_343 * kh_123[k]
                   + f_399 * kh_212[k]
                   - f_334 * kh_219[k]
                   - f_399 * kh_226[k]
                   + f_334 * kh_228[k]
                   - f_343 * kh_254[k]
                   + f_340 * kh_261[k]
                   + f_343 * kh_268[k]
                   - f_340 * kh_270[k]
                   - f_402 * kh_296[k]
                   + f_344 * kh_303[k]
                   + f_402 * kh_310[k]
                   - f_344 * kh_312[k]
                   + f_398 * kh_443[k]
                   - f_332 * kh_450[k]
                   - f_398 * kh_457[k]
                   + f_332 * kh_459[k]
                   - f_400 * kh_485[k]
                   + f_336 * kh_492[k]
                   + f_400 * kh_499[k]
                   - f_336 * kh_501[k]
                   + f_340 * kh_527[k]
                   - f_341 * kh_534[k]
                   - f_340 * kh_541[k]
                   + f_341 * kh_543[k];
    }

#pragma omp simd aligned(kh_0, kh_3, kh_5, kh_10, kh_12, kh_63, kh_66, kh_68, kh_73, kh_75, \
                         kh_105, kh_108, kh_110, kh_115, kh_117, kh_210, kh_213, kh_215, \
                         kh_220, kh_222, kh_252, kh_255, kh_257, kh_262, kh_264, kh_294, \
                         kh_297, kh_299, kh_304, kh_306, kh_441, kh_444, kh_446, kh_451, \
                         kh_453, kh_483, kh_486, kh_488, kh_493, kh_495, kh_525, kh_528, \
                         kh_530, kh_535, kh_537 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_118[k] = -f_319 * kh_0[k]
                   + f_318 * kh_3[k]
                   + f_320 * kh_5[k]
                   + f_307 * kh_10[k]
                   - f_308 * kh_12[k]
                   + f_319 * kh_63[k]
                   - f_318 * kh_66[k]
                   - f_320 * kh_68[k]
                   - f_307 * kh_73[k]
                   + f_308 * kh_75[k]
                   + f_328 * kh_105[k]
                   - f_313 * kh_108[k]
                   - f_325 * kh_110[k]
                   - f_316 * kh_115[k]
                   + f_317 * kh_117[k]
                   + f_312 * kh_210[k]
                   - f_310 * kh_213[k]
                   - f_313 * kh_215[k]
                   - f_309 * kh_220[k]
                   + f_311 * kh_222[k]
                   - f_313 * kh_252[k]
                   + f_321 * kh_255[k]
                   + f_323 * kh_257[k]
                   + f_311 * kh_262[k]
                   - f_322 * kh_264[k]
                   - f_330 * kh_294[k]
                   + f_329 * kh_297[k]
                   + f_331 * kh_299[k]
                   + f_321 * kh_304[k]
                   - f_327 * kh_306[k]
                   + f_307 * kh_441[k]
                   - f_305 * kh_444[k]
                   - f_308 * kh_446[k]
                   - f_304 * kh_451[k]
                   + f_306 * kh_453[k]
                   - f_316 * kh_483[k]
                   + f_311 * kh_486[k]
                   + f_317 * kh_488[k]
                   + f_314 * kh_493[k]
                   - f_315 * kh_495[k]
                   + f_321 * kh_525[k]
                   - f_325 * kh_528[k]
                   - f_327 * kh_530[k]
                   - f_324 * kh_535[k]
                   + f_326 * kh_537[k];
    }

#pragma omp simd aligned(kh_2, kh_7, kh_16, kh_65, kh_70, kh_79, kh_107, kh_112, kh_121, \
                         kh_212, kh_217, kh_226, kh_254, kh_259, kh_268, kh_296, kh_301, \
                         kh_310, kh_443, kh_448, kh_457, kh_485, kh_490, kh_499, kh_527, \
                         kh_532, kh_541 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_119[k] = f_409 * kh_2[k]
                   - f_410 * kh_7[k]
                   + f_409 * kh_16[k]
                   - f_409 * kh_65[k]
                   + f_410 * kh_70[k]
                   - f_409 * kh_79[k]
                   - f_297 * kh_107[k]
                   + f_413 * kh_112[k]
                   - f_297 * kh_121[k]
                   - f_405 * kh_212[k]
                   + f_406 * kh_217[k]
                   - f_405 * kh_226[k]
                   + f_411 * kh_254[k]
                   - f_298 * kh_259[k]
                   + f_411 * kh_268[k]
                   + f_414 * kh_296[k]
                   - f_300 * kh_301[k]
                   + f_414 * kh_310[k]
                   - f_403 * kh_443[k]
                   + f_404 * kh_448[k]
                   - f_403 * kh_457[k]
                   + f_407 * kh_485[k]
                   - f_408 * kh_490[k]
                   + f_407 * kh_499[k]
                   - f_302 * kh_527[k]
                   + f_412 * kh_532[k]
                   - f_302 * kh_541[k];
    }

#pragma omp simd aligned(kh_0, kh_3, kh_10, kh_63, kh_66, kh_73, kh_105, kh_108, kh_115, \
                         kh_210, kh_213, kh_220, kh_252, kh_255, kh_262, kh_294, kh_297, \
                         kh_304, kh_441, kh_444, kh_451, kh_483, kh_486, kh_493, kh_525, \
                         kh_528, kh_535 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_120[k] = f_285 * kh_0[k]
                   - f_284 * kh_3[k]
                   + f_280 * kh_10[k]
                   - f_285 * kh_63[k]
                   + f_284 * kh_66[k]
                   - f_280 * kh_73[k]
                   - f_292 * kh_105[k]
                   + f_286 * kh_108[k]
                   - f_291 * kh_115[k]
                   - f_280 * kh_210[k]
                   + f_279 * kh_213[k]
                   - f_278 * kh_220[k]
                   + f_288 * kh_252[k]
                   - f_287 * kh_255[k]
                   + f_286 * kh_262[k]
                   + f_295 * kh_294[k]
                   - f_294 * kh_297[k]
                   + f_293 * kh_304[k]
                   - f_277 * kh_441[k]
                   + f_276 * kh_444[k]
                   - f_275 * kh_451[k]
                   + f_283 * kh_483[k]
                   - f_282 * kh_486[k]
                   + f_281 * kh_493[k]
                   - f_290 * kh_525[k]
                   + f_289 * kh_528[k]
                   - f_287 * kh_535[k];
    }

#pragma omp simd aligned(kh_43, kh_48, kh_57, kh_148, kh_153, kh_162, kh_190, kh_195, kh_204, \
                         kh_337, kh_342, kh_351, kh_379, kh_384, kh_393, kh_610, kh_615, \
                         kh_624, kh_652, kh_657, kh_666 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_121[k] = -f_673 * kh_43[k]
                   + f_132 * kh_48[k]
                   - f_674 * kh_57[k]
                   + f_675 * kh_148[k]
                   - f_130 * kh_153[k]
                   + f_673 * kh_162[k]
                   + f_676 * kh_190[k]
                   - f_677 * kh_195[k]
                   + f_678 * kh_204[k]
                   + f_675 * kh_337[k]
                   - f_130 * kh_342[k]
                   + f_673 * kh_351[k]
                   - f_131 * kh_379[k]
                   + f_136 * kh_384[k]
                   - f_137 * kh_393[k]
                   - f_673 * kh_610[k]
                   + f_132 * kh_615[k]
                   - f_674 * kh_624[k]
                   + f_676 * kh_652[k]
                   - f_677 * kh_657[k]
                   + f_678 * kh_666[k];
    }

#pragma omp simd aligned(kh_46, kh_53, kh_151, kh_158, kh_193, kh_200, kh_340, kh_347, kh_382, \
                         kh_389, kh_613, kh_620, kh_655, kh_662 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_122[k] = -f_272 * kh_46[k]
                   + f_272 * kh_53[k]
                   + f_237 * kh_151[k]
                   - f_237 * kh_158[k]
                   + f_274 * kh_193[k]
                   - f_274 * kh_200[k]
                   + f_237 * kh_340[k]
                   - f_237 * kh_347[k]
                   - f_144 * kh_382[k]
                   + f_144 * kh_389[k]
                   - f_272 * kh_613[k]
                   + f_272 * kh_620[k]
                   + f_274 * kh_655[k]
                   - f_274 * kh_662[k];
    }

#pragma omp simd aligned(kh_43, kh_48, kh_50, kh_57, kh_59, kh_148, kh_153, kh_155, kh_162, \
                         kh_164, kh_190, kh_195, kh_197, kh_204, kh_206, kh_337, kh_342, \
                         kh_344, kh_351, kh_353, kh_379, kh_384, kh_386, kh_393, kh_395, \
                         kh_610, kh_615, kh_617, kh_624, kh_626, kh_652, kh_657, kh_659, \
                         kh_666, kh_668 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_123[k] = f_157 * kh_43[k]
                   + f_172 * kh_48[k]
                   - f_679 * kh_50[k]
                   - f_680 * kh_57[k]
                   + f_250 * kh_59[k]
                   - f_681 * kh_148[k]
                   - f_154 * kh_153[k]
                   + f_682 * kh_155[k]
                   + f_683 * kh_162[k]
                   - f_162 * kh_164[k]
                   - f_154 * kh_190[k]
                   - f_151 * kh_195[k]
                   + f_155 * kh_197[k]
                   + f_684 * kh_204[k]
                   - f_253 * kh_206[k]
                   - f_681 * kh_337[k]
                   - f_154 * kh_342[k]
                   + f_682 * kh_344[k]
                   + f_683 * kh_351[k]
                   - f_162 * kh_353[k]
                   + f_161 * kh_379[k]
                   + f_162 * kh_384[k]
                   - f_163 * kh_386[k]
                   - f_149 * kh_393[k]
                   + f_164 * kh_395[k]
                   + f_157 * kh_610[k]
                   + f_172 * kh_615[k]
                   - f_679 * kh_617[k]
                   - f_680 * kh_624[k]
                   + f_250 * kh_626[k]
                   - f_154 * kh_652[k]
                   - f_151 * kh_657[k]
                   + f_155 * kh_659[k]
                   + f_684 * kh_666[k]
                   - f_253 * kh_668[k];
    }

#pragma omp simd aligned(kh_46, kh_53, kh_55, kh_151, kh_158, kh_160, kh_193, kh_200, kh_202, \
                         kh_340, kh_347, kh_349, kh_382, kh_389, kh_391, kh_613, kh_620, \
                         kh_622, kh_655, kh_662, kh_664 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_124[k] = f_230 * kh_46[k]
                   + f_230 * kh_53[k]
                   - f_183 * kh_55[k]
                   - f_227 * kh_151[k]
                   - f_227 * kh_158[k]
                   + f_176 * kh_160[k]
                   - f_685 * kh_193[k]
                   - f_685 * kh_200[k]
                   + f_271 * kh_202[k]
                   - f_227 * kh_340[k]
                   - f_227 * kh_347[k]
                   + f_176 * kh_349[k]
                   + f_177 * kh_382[k]
                   + f_177 * kh_389[k]
                   - f_180 * kh_391[k]
                   + f_230 * kh_613[k]
                   + f_230 * kh_620[k]
                   - f_183 * kh_622[k]
                   - f_685 * kh_655[k]
                   - f_685 * kh_662[k]
                   + f_271 * kh_664[k];
    }

#pragma omp simd aligned(kh_43, kh_48, kh_50, kh_57, kh_59, kh_61, kh_148, kh_153, kh_155, \
                         kh_162, kh_164, kh_166, kh_190, kh_195, kh_197, kh_204, kh_206, \
                         kh_208, kh_337, kh_342, kh_344, kh_351, kh_353, kh_355, kh_379, \
                         kh_384, kh_386, kh_393, kh_395, kh_397, kh_610, kh_615, kh_617, \
                         kh_624, kh_626, kh_628, kh_652, kh_657, kh_659, kh_666, kh_668, \
                         kh_670 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_125[k] = -f_686 * kh_43[k]
                   - f_201 * kh_48[k]
                   + f_195 * kh_50[k]
                   - f_686 * kh_57[k]
                   + f_195 * kh_59[k]
                   - f_260 * kh_61[k]
                   + f_687 * kh_148[k]
                   + f_187 * kh_153[k]
                   - f_688 * kh_155[k]
                   + f_687 * kh_162[k]
                   - f_688 * kh_164[k]
                   + f_196 * kh_166[k]
                   + f_689 * kh_190[k]
                   + f_188 * kh_195[k]
                   - f_196 * kh_197[k]
                   + f_689 * kh_204[k]
                   - f_196 * kh_206[k]
                   + f_264 * kh_208[k]
                   + f_687 * kh_337[k]
                   + f_187 * kh_342[k]
                   - f_688 * kh_344[k]
                   + f_687 * kh_351[k]
                   - f_688 * kh_353[k]
                   + f_196 * kh_355[k]
                   - f_189 * kh_379[k]
                   - f_196 * kh_384[k]
                   + f_197 * kh_386[k]
                   - f_189 * kh_393[k]
                   + f_197 * kh_395[k]
                   - f_198 * kh_397[k]
                   - f_686 * kh_610[k]
                   - f_201 * kh_615[k]
                   + f_195 * kh_617[k]
                   - f_686 * kh_624[k]
                   + f_195 * kh_626[k]
                   - f_260 * kh_628[k]
                   + f_689 * kh_652[k]
                   + f_188 * kh_657[k]
                   - f_196 * kh_659[k]
                   + f_689 * kh_666[k]
                   - f_196 * kh_668[k]
                   + f_264 * kh_670[k];
    }

#pragma omp simd aligned(kh_44, kh_49, kh_51, kh_58, kh_60, kh_62, kh_149, kh_154, kh_156, \
                         kh_163, kh_165, kh_167, kh_191, kh_196, kh_198, kh_205, kh_207, \
                         kh_209, kh_338, kh_343, kh_345, kh_352, kh_354, kh_356, kh_380, \
                         kh_385, kh_387, kh_394, kh_396, kh_398, kh_611, kh_616, kh_618, \
                         kh_625, kh_627, kh_629, kh_653, kh_658, kh_660, kh_667, kh_669, \
                         kh_671 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_126[k] = -f_690 * kh_44[k]
                   - f_224 * kh_49[k]
                   + f_691 * kh_51[k]
                   - f_690 * kh_58[k]
                   + f_691 * kh_60[k]
                   - f_692 * kh_62[k]
                   + f_693 * kh_149[k]
                   + f_210 * kh_154[k]
                   - f_268 * kh_156[k]
                   + f_693 * kh_163[k]
                   - f_268 * kh_165[k]
                   + f_691 * kh_167[k]
                   + f_694 * kh_191[k]
                   + f_695 * kh_196[k]
                   - f_696 * kh_198[k]
                   + f_694 * kh_205[k]
                   - f_696 * kh_207[k]
                   + f_697 * kh_209[k]
                   + f_693 * kh_338[k]
                   + f_210 * kh_343[k]
                   - f_268 * kh_345[k]
                   + f_693 * kh_352[k]
                   - f_268 * kh_354[k]
                   + f_691 * kh_356[k]
                   - f_211 * kh_380[k]
                   - f_218 * kh_385[k]
                   + f_219 * kh_387[k]
                   - f_211 * kh_394[k]
                   + f_219 * kh_396[k]
                   - f_220 * kh_398[k]
                   - f_690 * kh_611[k]
                   - f_224 * kh_616[k]
                   + f_691 * kh_618[k]
                   - f_690 * kh_625[k]
                   + f_691 * kh_627[k]
                   - f_692 * kh_629[k]
                   + f_694 * kh_653[k]
                   + f_695 * kh_658[k]
                   - f_696 * kh_660[k]
                   + f_694 * kh_667[k]
                   - f_696 * kh_669[k]
                   + f_697 * kh_671[k];
    }

#pragma omp simd aligned(kh_42, kh_45, kh_47, kh_52, kh_54, kh_56, kh_147, kh_150, kh_152, \
                         kh_157, kh_159, kh_161, kh_189, kh_192, kh_194, kh_199, kh_201, \
                         kh_203, kh_336, kh_339, kh_341, kh_346, kh_348, kh_350, kh_378, \
                         kh_381, kh_383, kh_388, kh_390, kh_392, kh_609, kh_612, kh_614, \
                         kh_619, kh_621, kh_623, kh_651, kh_654, kh_656, kh_661, kh_663, \
                         kh_665 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_127[k] = -f_686 * kh_42[k]
                   - f_201 * kh_45[k]
                   + f_195 * kh_47[k]
                   - f_686 * kh_52[k]
                   + f_195 * kh_54[k]
                   - f_260 * kh_56[k]
                   + f_687 * kh_147[k]
                   + f_187 * kh_150[k]
                   - f_688 * kh_152[k]
                   + f_687 * kh_157[k]
                   - f_688 * kh_159[k]
                   + f_196 * kh_161[k]
                   + f_689 * kh_189[k]
                   + f_188 * kh_192[k]
                   - f_196 * kh_194[k]
                   + f_689 * kh_199[k]
                   - f_196 * kh_201[k]
                   + f_264 * kh_203[k]
                   + f_687 * kh_336[k]
                   + f_187 * kh_339[k]
                   - f_688 * kh_341[k]
                   + f_687 * kh_346[k]
                   - f_688 * kh_348[k]
                   + f_196 * kh_350[k]
                   - f_189 * kh_378[k]
                   - f_196 * kh_381[k]
                   + f_197 * kh_383[k]
                   - f_189 * kh_388[k]
                   + f_197 * kh_390[k]
                   - f_198 * kh_392[k]
                   - f_686 * kh_609[k]
                   - f_201 * kh_612[k]
                   + f_195 * kh_614[k]
                   - f_686 * kh_619[k]
                   + f_195 * kh_621[k]
                   - f_260 * kh_623[k]
                   + f_689 * kh_651[k]
                   + f_188 * kh_654[k]
                   - f_196 * kh_656[k]
                   + f_689 * kh_661[k]
                   - f_196 * kh_663[k]
                   + f_264 * kh_665[k];
    }

#pragma omp simd aligned(kh_44, kh_51, kh_58, kh_60, kh_149, kh_156, kh_163, kh_165, kh_191, \
                         kh_198, kh_205, kh_207, kh_338, kh_345, kh_352, kh_354, kh_380, \
                         kh_387, kh_394, kh_396, kh_611, kh_618, kh_625, kh_627, kh_653, \
                         kh_660, kh_667, kh_669 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_128[k] = f_698 * kh_44[k]
                   - f_230 * kh_51[k]
                   - f_698 * kh_58[k]
                   + f_230 * kh_60[k]
                   - f_699 * kh_149[k]
                   + f_227 * kh_156[k]
                   + f_699 * kh_163[k]
                   - f_227 * kh_165[k]
                   - f_175 * kh_191[k]
                   + f_685 * kh_198[k]
                   + f_175 * kh_205[k]
                   - f_685 * kh_207[k]
                   - f_699 * kh_338[k]
                   + f_227 * kh_345[k]
                   + f_699 * kh_352[k]
                   - f_227 * kh_354[k]
                   + f_176 * kh_380[k]
                   - f_177 * kh_387[k]
                   - f_176 * kh_394[k]
                   + f_177 * kh_396[k]
                   + f_698 * kh_611[k]
                   - f_230 * kh_618[k]
                   - f_698 * kh_625[k]
                   + f_230 * kh_627[k]
                   - f_175 * kh_653[k]
                   + f_685 * kh_660[k]
                   + f_175 * kh_667[k]
                   - f_685 * kh_669[k];
    }

#pragma omp simd aligned(kh_42, kh_45, kh_47, kh_52, kh_54, kh_147, kh_150, kh_152, kh_157, \
                         kh_159, kh_189, kh_192, kh_194, kh_199, kh_201, kh_336, kh_339, \
                         kh_341, kh_346, kh_348, kh_378, kh_381, kh_383, kh_388, kh_390, \
                         kh_609, kh_612, kh_614, kh_619, kh_621, kh_651, kh_654, kh_656, \
                         kh_661, kh_663 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_129[k] = f_680 * kh_42[k]
                   - f_172 * kh_45[k]
                   - f_250 * kh_47[k]
                   - f_157 * kh_52[k]
                   + f_679 * kh_54[k]
                   - f_683 * kh_147[k]
                   + f_154 * kh_150[k]
                   + f_162 * kh_152[k]
                   + f_681 * kh_157[k]
                   - f_682 * kh_159[k]
                   - f_684 * kh_189[k]
                   + f_151 * kh_192[k]
                   + f_253 * kh_194[k]
                   + f_154 * kh_199[k]
                   - f_155 * kh_201[k]
                   - f_683 * kh_336[k]
                   + f_154 * kh_339[k]
                   + f_162 * kh_341[k]
                   + f_681 * kh_346[k]
                   - f_682 * kh_348[k]
                   + f_149 * kh_378[k]
                   - f_162 * kh_381[k]
                   - f_164 * kh_383[k]
                   - f_161 * kh_388[k]
                   + f_163 * kh_390[k]
                   + f_680 * kh_609[k]
                   - f_172 * kh_612[k]
                   - f_250 * kh_614[k]
                   - f_157 * kh_619[k]
                   + f_679 * kh_621[k]
                   - f_684 * kh_651[k]
                   + f_151 * kh_654[k]
                   + f_253 * kh_656[k]
                   + f_154 * kh_661[k]
                   - f_155 * kh_663[k];
    }

#pragma omp simd aligned(kh_44, kh_49, kh_58, kh_149, kh_154, kh_163, kh_191, kh_196, kh_205, \
                         kh_338, kh_343, kh_352, kh_380, kh_385, kh_394, kh_611, kh_616, \
                         kh_625, kh_653, kh_658, kh_667 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_130[k] = -f_240 * kh_44[k]
                   + f_143 * kh_49[k]
                   - f_240 * kh_58[k]
                   + f_232 * kh_149[k]
                   - f_700 * kh_154[k]
                   + f_232 * kh_163[k]
                   + f_141 * kh_191[k]
                   - f_237 * kh_196[k]
                   + f_141 * kh_205[k]
                   + f_232 * kh_338[k]
                   - f_700 * kh_343[k]
                   + f_232 * kh_352[k]
                   - f_237 * kh_380[k]
                   + f_238 * kh_385[k]
                   - f_237 * kh_394[k]
                   - f_240 * kh_611[k]
                   + f_143 * kh_616[k]
                   - f_240 * kh_625[k]
                   + f_141 * kh_653[k]
                   - f_237 * kh_658[k]
                   + f_141 * kh_667[k];
    }

#pragma omp simd aligned(kh_42, kh_45, kh_52, kh_147, kh_150, kh_157, kh_189, kh_192, kh_199, \
                         kh_336, kh_339, kh_346, kh_378, kh_381, kh_388, kh_609, kh_612, \
                         kh_619, kh_651, kh_654, kh_661 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_131[k] = -f_674 * kh_42[k]
                   + f_132 * kh_45[k]
                   - f_673 * kh_52[k]
                   + f_673 * kh_147[k]
                   - f_130 * kh_150[k]
                   + f_675 * kh_157[k]
                   + f_678 * kh_189[k]
                   - f_677 * kh_192[k]
                   + f_676 * kh_199[k]
                   + f_673 * kh_336[k]
                   - f_130 * kh_339[k]
                   + f_675 * kh_346[k]
                   - f_137 * kh_378[k]
                   + f_136 * kh_381[k]
                   - f_131 * kh_388[k]
                   - f_674 * kh_609[k]
                   + f_132 * kh_612[k]
                   - f_673 * kh_619[k]
                   + f_678 * kh_651[k]
                   - f_677 * kh_654[k]
                   + f_676 * kh_661[k];
    }

#pragma omp simd aligned(kh_1, kh_6, kh_15, kh_64, kh_69, kh_78, kh_106, kh_111, kh_120, \
                         kh_211, kh_216, kh_225, kh_253, kh_258, kh_267, kh_442, kh_447, \
                         kh_456, kh_484, kh_489, kh_498 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_132[k] = -f_129 * kh_1[k]
                   + f_138 * kh_6[k]
                   - f_139 * kh_15[k]
                   + f_133 * kh_64[k]
                   - f_134 * kh_69[k]
                   + f_135 * kh_78[k]
                   + f_132 * kh_106[k]
                   - f_137 * kh_111[k]
                   + f_140 * kh_120[k]
                   + f_127 * kh_211[k]
                   - f_128 * kh_216[k]
                   + f_129 * kh_225[k]
                   - f_131 * kh_253[k]
                   + f_136 * kh_258[k]
                   - f_137 * kh_267[k]
                   - f_127 * kh_442[k]
                   + f_128 * kh_447[k]
                   - f_129 * kh_456[k]
                   + f_130 * kh_484[k]
                   - f_131 * kh_489[k]
                   + f_132 * kh_498[k];
    }

#pragma omp simd aligned(kh_4, kh_11, kh_67, kh_74, kh_109, kh_116, kh_214, kh_221, kh_256, \
                         kh_263, kh_445, kh_452, kh_487, kh_494 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_133[k] = -f_145 * kh_4[k]
                   + f_145 * kh_11[k]
                   + f_143 * kh_67[k]
                   - f_143 * kh_74[k]
                   + f_146 * kh_109[k]
                   - f_146 * kh_116[k]
                   + f_141 * kh_214[k]
                   - f_141 * kh_221[k]
                   - f_144 * kh_256[k]
                   + f_144 * kh_263[k]
                   - f_141 * kh_445[k]
                   + f_141 * kh_452[k]
                   + f_142 * kh_487[k]
                   - f_142 * kh_494[k];
    }

#pragma omp simd aligned(kh_1, kh_6, kh_8, kh_15, kh_17, kh_64, kh_69, kh_71, kh_78, kh_80, \
                         kh_106, kh_111, kh_113, kh_120, kh_122, kh_211, kh_216, kh_218, \
                         kh_225, kh_227, kh_253, kh_258, kh_260, kh_267, kh_269, kh_442, \
                         kh_447, kh_449, kh_456, kh_458, kh_484, kh_489, kh_491, kh_498, \
                         kh_500 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_134[k] = f_165 * kh_1[k]
                   + f_166 * kh_6[k]
                   - f_167 * kh_8[k]
                   - f_168 * kh_15[k]
                   + f_169 * kh_17[k]
                   - f_156 * kh_64[k]
                   - f_157 * kh_69[k]
                   + f_158 * kh_71[k]
                   + f_159 * kh_78[k]
                   - f_160 * kh_80[k]
                   - f_170 * kh_106[k]
                   - f_167 * kh_111[k]
                   + f_171 * kh_113[k]
                   + f_172 * kh_120[k]
                   - f_173 * kh_122[k]
                   - f_147 * kh_211[k]
                   - f_148 * kh_216[k]
                   + f_149 * kh_218[k]
                   + f_150 * kh_225[k]
                   - f_151 * kh_227[k]
                   + f_161 * kh_253[k]
                   + f_162 * kh_258[k]
                   - f_163 * kh_260[k]
                   - f_149 * kh_267[k]
                   + f_164 * kh_269[k]
                   + f_147 * kh_442[k]
                   + f_148 * kh_447[k]
                   - f_149 * kh_449[k]
                   - f_150 * kh_456[k]
                   + f_151 * kh_458[k]
                   - f_152 * kh_484[k]
                   - f_149 * kh_489[k]
                   + f_153 * kh_491[k]
                   + f_154 * kh_498[k]
                   - f_155 * kh_500[k];
    }

#pragma omp simd aligned(kh_4, kh_11, kh_13, kh_67, kh_74, kh_76, kh_109, kh_116, kh_118, \
                         kh_214, kh_221, kh_223, kh_256, kh_263, kh_265, kh_445, kh_452, \
                         kh_454, kh_487, kh_494, kh_496 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_135[k] = f_181 * kh_4[k]
                   + f_181 * kh_11[k]
                   - f_182 * kh_13[k]
                   - f_178 * kh_67[k]
                   - f_178 * kh_74[k]
                   + f_179 * kh_76[k]
                   - f_183 * kh_109[k]
                   - f_183 * kh_116[k]
                   + f_184 * kh_118[k]
                   - f_174 * kh_214[k]
                   - f_174 * kh_221[k]
                   + f_175 * kh_223[k]
                   + f_177 * kh_256[k]
                   + f_177 * kh_263[k]
                   - f_180 * kh_265[k]
                   + f_174 * kh_445[k]
                   + f_174 * kh_452[k]
                   - f_175 * kh_454[k]
                   - f_176 * kh_487[k]
                   - f_176 * kh_494[k]
                   + f_177 * kh_496[k];
    }

#pragma omp simd aligned(kh_1, kh_6, kh_8, kh_15, kh_17, kh_19, kh_64, kh_69, kh_71, kh_78, \
                         kh_80, kh_82, kh_106, kh_111, kh_113, kh_120, kh_122, kh_124, kh_211, \
                         kh_216, kh_218, kh_225, kh_227, kh_229, kh_253, kh_258, kh_260, \
                         kh_267, kh_269, kh_271, kh_442, kh_447, kh_449, kh_456, kh_458, \
                         kh_460, kh_484, kh_489, kh_491, kh_498, kh_500, \
                         kh_502 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_136[k] = -f_199 * kh_1[k]
                   - f_200 * kh_6[k]
                   + f_201 * kh_8[k]
                   - f_199 * kh_15[k]
                   + f_201 * kh_17[k]
                   - f_202 * kh_19[k]
                   + f_192 * kh_64[k]
                   + f_193 * kh_69[k]
                   - f_194 * kh_71[k]
                   + f_192 * kh_78[k]
                   - f_194 * kh_80[k]
                   + f_195 * kh_82[k]
                   + f_201 * kh_106[k]
                   + f_203 * kh_111[k]
                   - f_204 * kh_113[k]
                   + f_201 * kh_120[k]
                   - f_204 * kh_122[k]
                   + f_205 * kh_124[k]
                   + f_185 * kh_211[k]
                   + f_186 * kh_216[k]
                   - f_187 * kh_218[k]
                   + f_185 * kh_225[k]
                   - f_187 * kh_227[k]
                   + f_188 * kh_229[k]
                   - f_189 * kh_253[k]
                   - f_196 * kh_258[k]
                   + f_197 * kh_260[k]
                   - f_189 * kh_267[k]
                   + f_197 * kh_269[k]
                   - f_198 * kh_271[k]
                   - f_185 * kh_442[k]
                   - f_186 * kh_447[k]
                   + f_187 * kh_449[k]
                   - f_185 * kh_456[k]
                   + f_187 * kh_458[k]
                   - f_188 * kh_460[k]
                   + f_187 * kh_484[k]
                   + f_189 * kh_489[k]
                   - f_190 * kh_491[k]
                   + f_187 * kh_498[k]
                   - f_190 * kh_500[k]
                   + f_191 * kh_502[k];
    }

#pragma omp simd aligned(kh_2, kh_7, kh_9, kh_16, kh_18, kh_20, kh_65, kh_70, kh_72, kh_79, \
                         kh_81, kh_83, kh_107, kh_112, kh_114, kh_121, kh_123, kh_125, kh_212, \
                         kh_217, kh_219, kh_226, kh_228, kh_230, kh_254, kh_259, kh_261, \
                         kh_268, kh_270, kh_272, kh_443, kh_448, kh_450, kh_457, kh_459, \
                         kh_461, kh_485, kh_490, kh_492, kh_499, kh_501, \
                         kh_503 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_137[k] = -f_221 * kh_2[k]
                   - f_222 * kh_7[k]
                   + f_209 * kh_9[k]
                   - f_221 * kh_16[k]
                   + f_209 * kh_18[k]
                   - f_223 * kh_20[k]
                   + f_214 * kh_65[k]
                   + f_215 * kh_70[k]
                   - f_216 * kh_72[k]
                   + f_214 * kh_79[k]
                   - f_216 * kh_81[k]
                   + f_217 * kh_83[k]
                   + f_224 * kh_107[k]
                   + f_216 * kh_112[k]
                   - f_213 * kh_114[k]
                   + f_224 * kh_121[k]
                   - f_213 * kh_123[k]
                   + f_225 * kh_125[k]
                   + f_206 * kh_212[k]
                   + f_207 * kh_217[k]
                   - f_208 * kh_219[k]
                   + f_206 * kh_226[k]
                   - f_208 * kh_228[k]
                   + f_209 * kh_230[k]
                   - f_211 * kh_254[k]
                   - f_218 * kh_259[k]
                   + f_219 * kh_261[k]
                   - f_211 * kh_268[k]
                   + f_219 * kh_270[k]
                   - f_220 * kh_272[k]
                   - f_206 * kh_443[k]
                   - f_207 * kh_448[k]
                   + f_208 * kh_450[k]
                   - f_206 * kh_457[k]
                   + f_208 * kh_459[k]
                   - f_209 * kh_461[k]
                   + f_210 * kh_485[k]
                   + f_211 * kh_490[k]
                   - f_212 * kh_492[k]
                   + f_210 * kh_499[k]
                   - f_212 * kh_501[k]
                   + f_213 * kh_503[k];
    }

#pragma omp simd aligned(kh_0, kh_3, kh_5, kh_10, kh_12, kh_14, kh_63, kh_66, kh_68, kh_73, \
                         kh_75, kh_77, kh_105, kh_108, kh_110, kh_115, kh_117, kh_119, kh_210, \
                         kh_213, kh_215, kh_220, kh_222, kh_224, kh_252, kh_255, kh_257, \
                         kh_262, kh_264, kh_266, kh_441, kh_444, kh_446, kh_451, kh_453, \
                         kh_455, kh_483, kh_486, kh_488, kh_493, kh_495, \
                         kh_497 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_138[k] = -f_199 * kh_0[k]
                   - f_200 * kh_3[k]
                   + f_201 * kh_5[k]
                   - f_199 * kh_10[k]
                   + f_201 * kh_12[k]
                   - f_202 * kh_14[k]
                   + f_192 * kh_63[k]
                   + f_193 * kh_66[k]
                   - f_194 * kh_68[k]
                   + f_192 * kh_73[k]
                   - f_194 * kh_75[k]
                   + f_195 * kh_77[k]
                   + f_201 * kh_105[k]
                   + f_203 * kh_108[k]
                   - f_204 * kh_110[k]
                   + f_201 * kh_115[k]
                   - f_204 * kh_117[k]
                   + f_205 * kh_119[k]
                   + f_185 * kh_210[k]
                   + f_186 * kh_213[k]
                   - f_187 * kh_215[k]
                   + f_185 * kh_220[k]
                   - f_187 * kh_222[k]
                   + f_188 * kh_224[k]
                   - f_189 * kh_252[k]
                   - f_196 * kh_255[k]
                   + f_197 * kh_257[k]
                   - f_189 * kh_262[k]
                   + f_197 * kh_264[k]
                   - f_198 * kh_266[k]
                   - f_185 * kh_441[k]
                   - f_186 * kh_444[k]
                   + f_187 * kh_446[k]
                   - f_185 * kh_451[k]
                   + f_187 * kh_453[k]
                   - f_188 * kh_455[k]
                   + f_187 * kh_483[k]
                   + f_189 * kh_486[k]
                   - f_190 * kh_488[k]
                   + f_187 * kh_493[k]
                   - f_190 * kh_495[k]
                   + f_191 * kh_497[k];
    }

#pragma omp simd aligned(kh_2, kh_9, kh_16, kh_18, kh_65, kh_72, kh_79, kh_81, kh_107, kh_114, \
                         kh_121, kh_123, kh_212, kh_219, kh_226, kh_228, kh_254, kh_261, \
                         kh_268, kh_270, kh_443, kh_450, kh_457, kh_459, kh_485, kh_492, \
                         kh_499, kh_501 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_139[k] = f_229 * kh_2[k]
                   - f_181 * kh_9[k]
                   - f_229 * kh_16[k]
                   + f_181 * kh_18[k]
                   - f_228 * kh_65[k]
                   + f_178 * kh_72[k]
                   + f_228 * kh_79[k]
                   - f_178 * kh_81[k]
                   - f_230 * kh_107[k]
                   + f_183 * kh_114[k]
                   + f_230 * kh_121[k]
                   - f_183 * kh_123[k]
                   - f_226 * kh_212[k]
                   + f_174 * kh_219[k]
                   + f_226 * kh_226[k]
                   - f_174 * kh_228[k]
                   + f_176 * kh_254[k]
                   - f_177 * kh_261[k]
                   - f_176 * kh_268[k]
                   + f_177 * kh_270[k]
                   + f_226 * kh_443[k]
                   - f_174 * kh_450[k]
                   - f_226 * kh_457[k]
                   + f_174 * kh_459[k]
                   - f_227 * kh_485[k]
                   + f_176 * kh_492[k]
                   + f_227 * kh_499[k]
                   - f_176 * kh_501[k];
    }

#pragma omp simd aligned(kh_0, kh_3, kh_5, kh_10, kh_12, kh_63, kh_66, kh_68, kh_73, kh_75, \
                         kh_105, kh_108, kh_110, kh_115, kh_117, kh_210, kh_213, kh_215, \
                         kh_220, kh_222, kh_252, kh_255, kh_257, kh_262, kh_264, kh_441, \
                         kh_444, kh_446, kh_451, kh_453, kh_483, kh_486, kh_488, kh_493, \
                         kh_495 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_140[k] = f_168 * kh_0[k]
                   - f_166 * kh_3[k]
                   - f_169 * kh_5[k]
                   - f_165 * kh_10[k]
                   + f_167 * kh_12[k]
                   - f_159 * kh_63[k]
                   + f_157 * kh_66[k]
                   + f_160 * kh_68[k]
                   + f_156 * kh_73[k]
                   - f_158 * kh_75[k]
                   - f_172 * kh_105[k]
                   + f_167 * kh_108[k]
                   + f_173 * kh_110[k]
                   + f_170 * kh_115[k]
                   - f_171 * kh_117[k]
                   - f_150 * kh_210[k]
                   + f_148 * kh_213[k]
                   + f_151 * kh_215[k]
                   + f_147 * kh_220[k]
                   - f_149 * kh_222[k]
                   + f_149 * kh_252[k]
                   - f_162 * kh_255[k]
                   - f_164 * kh_257[k]
                   - f_161 * kh_262[k]
                   + f_163 * kh_264[k]
                   + f_150 * kh_441[k]
                   - f_148 * kh_444[k]
                   - f_151 * kh_446[k]
                   - f_147 * kh_451[k]
                   + f_149 * kh_453[k]
                   - f_154 * kh_483[k]
                   + f_149 * kh_486[k]
                   + f_155 * kh_488[k]
                   + f_152 * kh_493[k]
                   - f_153 * kh_495[k];
    }

#pragma omp simd aligned(kh_2, kh_7, kh_16, kh_65, kh_70, kh_79, kh_107, kh_112, kh_121, \
                         kh_212, kh_217, kh_226, kh_254, kh_259, kh_268, kh_443, kh_448, \
                         kh_457, kh_485, kh_490, kh_499 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_141[k] = -f_239 * kh_2[k]
                   + f_240 * kh_7[k]
                   - f_239 * kh_16[k]
                   + f_235 * kh_65[k]
                   - f_236 * kh_70[k]
                   + f_235 * kh_79[k]
                   + f_241 * kh_107[k]
                   - f_242 * kh_112[k]
                   + f_241 * kh_121[k]
                   + f_231 * kh_212[k]
                   - f_232 * kh_217[k]
                   + f_231 * kh_226[k]
                   - f_237 * kh_254[k]
                   + f_238 * kh_259[k]
                   - f_237 * kh_268[k]
                   - f_231 * kh_443[k]
                   + f_232 * kh_448[k]
                   - f_231 * kh_457[k]
                   + f_233 * kh_485[k]
                   - f_234 * kh_490[k]
                   + f_233 * kh_499[k];
    }

#pragma omp simd aligned(kh_0, kh_3, kh_10, kh_63, kh_66, kh_73, kh_105, kh_108, kh_115, \
                         kh_210, kh_213, kh_220, kh_252, kh_255, kh_262, kh_441, kh_444, \
                         kh_451, kh_483, kh_486, kh_493 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_142[k] = -f_139 * kh_0[k]
                   + f_138 * kh_3[k]
                   - f_129 * kh_10[k]
                   + f_135 * kh_63[k]
                   - f_134 * kh_66[k]
                   + f_133 * kh_73[k]
                   + f_140 * kh_105[k]
                   - f_137 * kh_108[k]
                   + f_132 * kh_115[k]
                   + f_129 * kh_210[k]
                   - f_128 * kh_213[k]
                   + f_127 * kh_220[k]
                   - f_137 * kh_252[k]
                   + f_136 * kh_255[k]
                   - f_131 * kh_262[k]
                   - f_129 * kh_441[k]
                   + f_128 * kh_444[k]
                   - f_127 * kh_451[k]
                   + f_132 * kh_483[k]
                   - f_131 * kh_486[k]
                   + f_130 * kh_493[k];
    }

#pragma omp simd aligned(kh_43, kh_48, kh_57, kh_148, kh_153, kh_162, kh_337, kh_342, kh_351, \
                         kh_610, kh_615, kh_624 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_143[k] = f_57 * kh_43[k]
                   - f_58 * kh_48[k]
                   + f_701 * kh_57[k]
                   - f_702 * kh_148[k]
                   + f_703 * kh_153[k]
                   - f_64 * kh_162[k]
                   + f_702 * kh_337[k]
                   - f_703 * kh_342[k]
                   + f_64 * kh_351[k]
                   - f_57 * kh_610[k]
                   + f_58 * kh_615[k]
                   - f_701 * kh_624[k];
    }

#pragma omp simd aligned(kh_46, kh_53, kh_151, kh_158, kh_340, kh_347, kh_613, \
                         kh_620 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_144[k] = f_704 * kh_46[k]
                   - f_704 * kh_53[k]
                   - f_705 * kh_151[k]
                   + f_705 * kh_158[k]
                   + f_705 * kh_340[k]
                   - f_705 * kh_347[k]
                   - f_704 * kh_613[k]
                   + f_704 * kh_620[k];
    }

#pragma omp simd aligned(kh_43, kh_48, kh_50, kh_57, kh_59, kh_148, kh_153, kh_155, kh_162, \
                         kh_164, kh_337, kh_342, kh_344, kh_351, kh_353, kh_610, kh_615, \
                         kh_617, kh_624, kh_626 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_145[k] = -f_706 * kh_43[k]
                   - f_707 * kh_48[k]
                   + f_708 * kh_50[k]
                   + f_709 * kh_57[k]
                   - f_710 * kh_59[k]
                   + f_711 * kh_148[k]
                   + f_712 * kh_153[k]
                   - f_713 * kh_155[k]
                   - f_714 * kh_162[k]
                   + f_715 * kh_164[k]
                   - f_711 * kh_337[k]
                   - f_712 * kh_342[k]
                   + f_713 * kh_344[k]
                   + f_714 * kh_351[k]
                   - f_715 * kh_353[k]
                   + f_706 * kh_610[k]
                   + f_707 * kh_615[k]
                   - f_708 * kh_617[k]
                   - f_709 * kh_624[k]
                   + f_710 * kh_626[k];
    }

#pragma omp simd aligned(kh_46, kh_53, kh_55, kh_151, kh_158, kh_160, kh_340, kh_347, kh_349, \
                         kh_613, kh_620, kh_622 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_146[k] = -f_716 * kh_46[k]
                   - f_716 * kh_53[k]
                   + f_717 * kh_55[k]
                   + f_718 * kh_151[k]
                   + f_718 * kh_158[k]
                   - f_719 * kh_160[k]
                   - f_718 * kh_340[k]
                   - f_718 * kh_347[k]
                   + f_719 * kh_349[k]
                   + f_716 * kh_613[k]
                   + f_716 * kh_620[k]
                   - f_717 * kh_622[k];
    }

#pragma omp simd aligned(kh_43, kh_48, kh_50, kh_57, kh_59, kh_61, kh_148, kh_153, kh_155, \
                         kh_162, kh_164, kh_166, kh_337, kh_342, kh_344, kh_351, kh_353, \
                         kh_355, kh_610, kh_615, kh_617, kh_624, kh_626, \
                         kh_628 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_147[k] = f_720 * kh_43[k]
                   + f_721 * kh_48[k]
                   - f_106 * kh_50[k]
                   + f_720 * kh_57[k]
                   - f_106 * kh_59[k]
                   + f_722 * kh_61[k]
                   - f_723 * kh_148[k]
                   - f_724 * kh_153[k]
                   + f_725 * kh_155[k]
                   - f_723 * kh_162[k]
                   + f_725 * kh_164[k]
                   - f_726 * kh_166[k]
                   + f_723 * kh_337[k]
                   + f_724 * kh_342[k]
                   - f_725 * kh_344[k]
                   + f_723 * kh_351[k]
                   - f_725 * kh_353[k]
                   + f_726 * kh_355[k]
                   - f_720 * kh_610[k]
                   - f_721 * kh_615[k]
                   + f_106 * kh_617[k]
                   - f_720 * kh_624[k]
                   + f_106 * kh_626[k]
                   - f_722 * kh_628[k];
    }

#pragma omp simd aligned(kh_44, kh_49, kh_51, kh_58, kh_60, kh_62, kh_149, kh_154, kh_156, \
                         kh_163, kh_165, kh_167, kh_338, kh_343, kh_345, kh_352, kh_354, \
                         kh_356, kh_611, kh_616, kh_618, kh_625, kh_627, \
                         kh_629 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_148[k] = f_9 * kh_44[k]
                   + f_727 * kh_49[k]
                   - f_728 * kh_51[k]
                   + f_9 * kh_58[k]
                   - f_728 * kh_60[k]
                   + f_729 * kh_62[k]
                   - f_730 * kh_149[k]
                   - f_731 * kh_154[k]
                   + f_118 * kh_156[k]
                   - f_730 * kh_163[k]
                   + f_118 * kh_165[k]
                   - f_732 * kh_167[k]
                   + f_730 * kh_338[k]
                   + f_731 * kh_343[k]
                   - f_118 * kh_345[k]
                   + f_730 * kh_352[k]
                   - f_118 * kh_354[k]
                   + f_732 * kh_356[k]
                   - f_9 * kh_611[k]
                   - f_727 * kh_616[k]
                   + f_728 * kh_618[k]
                   - f_9 * kh_625[k]
                   + f_728 * kh_627[k]
                   - f_729 * kh_629[k];
    }

#pragma omp simd aligned(kh_42, kh_45, kh_47, kh_52, kh_54, kh_56, kh_147, kh_150, kh_152, \
                         kh_157, kh_159, kh_161, kh_336, kh_339, kh_341, kh_346, kh_348, \
                         kh_350, kh_609, kh_612, kh_614, kh_619, kh_621, \
                         kh_623 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_149[k] = f_720 * kh_42[k]
                   + f_721 * kh_45[k]
                   - f_106 * kh_47[k]
                   + f_720 * kh_52[k]
                   - f_106 * kh_54[k]
                   + f_722 * kh_56[k]
                   - f_723 * kh_147[k]
                   - f_724 * kh_150[k]
                   + f_725 * kh_152[k]
                   - f_723 * kh_157[k]
                   + f_725 * kh_159[k]
                   - f_726 * kh_161[k]
                   + f_723 * kh_336[k]
                   + f_724 * kh_339[k]
                   - f_725 * kh_341[k]
                   + f_723 * kh_346[k]
                   - f_725 * kh_348[k]
                   + f_726 * kh_350[k]
                   - f_720 * kh_609[k]
                   - f_721 * kh_612[k]
                   + f_106 * kh_614[k]
                   - f_720 * kh_619[k]
                   + f_106 * kh_621[k]
                   - f_722 * kh_623[k];
    }

#pragma omp simd aligned(kh_44, kh_51, kh_58, kh_60, kh_149, kh_156, kh_163, kh_165, kh_338, \
                         kh_345, kh_352, kh_354, kh_611, kh_618, kh_625, \
                         kh_627 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_150[k] = -f_733 * kh_44[k]
                   + f_716 * kh_51[k]
                   + f_733 * kh_58[k]
                   - f_716 * kh_60[k]
                   + f_734 * kh_149[k]
                   - f_718 * kh_156[k]
                   - f_734 * kh_163[k]
                   + f_718 * kh_165[k]
                   - f_734 * kh_338[k]
                   + f_718 * kh_345[k]
                   + f_734 * kh_352[k]
                   - f_718 * kh_354[k]
                   + f_733 * kh_611[k]
                   - f_716 * kh_618[k]
                   - f_733 * kh_625[k]
                   + f_716 * kh_627[k];
    }

#pragma omp simd aligned(kh_42, kh_45, kh_47, kh_52, kh_54, kh_147, kh_150, kh_152, kh_157, \
                         kh_159, kh_336, kh_339, kh_341, kh_346, kh_348, kh_609, kh_612, \
                         kh_614, kh_619, kh_621 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_151[k] = -f_709 * kh_42[k]
                   + f_707 * kh_45[k]
                   + f_710 * kh_47[k]
                   + f_706 * kh_52[k]
                   - f_708 * kh_54[k]
                   + f_714 * kh_147[k]
                   - f_712 * kh_150[k]
                   - f_715 * kh_152[k]
                   - f_711 * kh_157[k]
                   + f_713 * kh_159[k]
                   - f_714 * kh_336[k]
                   + f_712 * kh_339[k]
                   + f_715 * kh_341[k]
                   + f_711 * kh_346[k]
                   - f_713 * kh_348[k]
                   + f_709 * kh_609[k]
                   - f_707 * kh_612[k]
                   - f_710 * kh_614[k]
                   - f_706 * kh_619[k]
                   + f_708 * kh_621[k];
    }

#pragma omp simd aligned(kh_44, kh_49, kh_58, kh_149, kh_154, kh_163, kh_338, kh_343, kh_352, \
                         kh_611, kh_616, kh_625 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_152[k] = f_735 * kh_44[k]
                   - f_123 * kh_49[k]
                   + f_735 * kh_58[k]
                   - f_736 * kh_149[k]
                   + f_737 * kh_154[k]
                   - f_736 * kh_163[k]
                   + f_736 * kh_338[k]
                   - f_737 * kh_343[k]
                   + f_736 * kh_352[k]
                   - f_735 * kh_611[k]
                   + f_123 * kh_616[k]
                   - f_735 * kh_625[k];
    }

#pragma omp simd aligned(kh_42, kh_45, kh_52, kh_147, kh_150, kh_157, kh_336, kh_339, kh_346, \
                         kh_609, kh_612, kh_619 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_153[k] = f_701 * kh_42[k]
                   - f_58 * kh_45[k]
                   + f_57 * kh_52[k]
                   - f_64 * kh_147[k]
                   + f_703 * kh_150[k]
                   - f_702 * kh_157[k]
                   + f_64 * kh_336[k]
                   - f_703 * kh_339[k]
                   + f_702 * kh_346[k]
                   - f_701 * kh_609[k]
                   + f_58 * kh_612[k]
                   - f_57 * kh_619[k];
    }

#pragma omp simd aligned(kh_1, kh_6, kh_15, kh_64, kh_69, kh_78, kh_211, kh_216, kh_225, \
                         kh_442, kh_447, kh_456 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_154[k] = f_8 * kh_1[k]
                   - f_9 * kh_6[k]
                   + f_10 * kh_15[k]
                   - f_5 * kh_64[k]
                   + f_6 * kh_69[k]
                   - f_7 * kh_78[k]
                   + f_3 * kh_211[k]
                   - f_4 * kh_216[k]
                   + f_0 * kh_225[k]
                   - f_0 * kh_442[k]
                   + f_1 * kh_447[k]
                   - f_2 * kh_456[k];
    }

#pragma omp simd aligned(kh_4, kh_11, kh_67, kh_74, kh_214, kh_221, kh_445, \
                         kh_452 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_155[k] = f_14 * kh_4[k]
                   - f_14 * kh_11[k]
                   - f_13 * kh_67[k]
                   + f_13 * kh_74[k]
                   + f_12 * kh_214[k]
                   - f_12 * kh_221[k]
                   - f_11 * kh_445[k]
                   + f_11 * kh_452[k];
    }

#pragma omp simd aligned(kh_1, kh_6, kh_8, kh_15, kh_17, kh_64, kh_69, kh_71, kh_78, kh_80, \
                         kh_211, kh_216, kh_218, kh_225, kh_227, kh_442, kh_447, kh_449, \
                         kh_456, kh_458 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_156[k] = -f_28 * kh_1[k]
                   - f_29 * kh_6[k]
                   + f_30 * kh_8[k]
                   + f_31 * kh_15[k]
                   - f_32 * kh_17[k]
                   + f_25 * kh_64[k]
                   + f_26 * kh_69[k]
                   - f_27 * kh_71[k]
                   - f_15 * kh_78[k]
                   + f_17 * kh_80[k]
                   - f_20 * kh_211[k]
                   - f_21 * kh_216[k]
                   + f_22 * kh_218[k]
                   + f_23 * kh_225[k]
                   - f_24 * kh_227[k]
                   + f_15 * kh_442[k]
                   + f_16 * kh_447[k]
                   - f_17 * kh_449[k]
                   - f_18 * kh_456[k]
                   + f_19 * kh_458[k];
    }

#pragma omp simd aligned(kh_4, kh_11, kh_13, kh_67, kh_74, kh_76, kh_214, kh_221, kh_223, \
                         kh_445, kh_452, kh_454 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_157[k] = -f_39 * kh_4[k]
                   - f_39 * kh_11[k]
                   + f_40 * kh_13[k]
                   + f_37 * kh_67[k]
                   + f_37 * kh_74[k]
                   - f_38 * kh_76[k]
                   - f_35 * kh_214[k]
                   - f_35 * kh_221[k]
                   + f_36 * kh_223[k]
                   + f_33 * kh_445[k]
                   + f_33 * kh_452[k]
                   - f_34 * kh_454[k];
    }

#pragma omp simd aligned(kh_1, kh_6, kh_8, kh_15, kh_17, kh_19, kh_64, kh_69, kh_71, kh_78, \
                         kh_80, kh_82, kh_211, kh_216, kh_218, kh_225, kh_227, kh_229, kh_442, \
                         kh_447, kh_449, kh_456, kh_458, kh_460 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_158[k] = f_53 * kh_1[k]
                   + f_54 * kh_6[k]
                   - f_55 * kh_8[k]
                   + f_53 * kh_15[k]
                   - f_55 * kh_17[k]
                   + f_56 * kh_19[k]
                   - f_49 * kh_64[k]
                   - f_50 * kh_69[k]
                   + f_51 * kh_71[k]
                   - f_49 * kh_78[k]
                   + f_51 * kh_80[k]
                   - f_52 * kh_82[k]
                   + f_45 * kh_211[k]
                   + f_46 * kh_216[k]
                   - f_47 * kh_218[k]
                   + f_45 * kh_225[k]
                   - f_47 * kh_227[k]
                   + f_48 * kh_229[k]
                   - f_41 * kh_442[k]
                   - f_42 * kh_447[k]
                   + f_43 * kh_449[k]
                   - f_41 * kh_456[k]
                   + f_43 * kh_458[k]
                   - f_44 * kh_460[k];
    }

#pragma omp simd aligned(kh_2, kh_7, kh_9, kh_16, kh_18, kh_20, kh_65, kh_70, kh_72, kh_79, \
                         kh_81, kh_83, kh_212, kh_217, kh_219, kh_226, kh_228, kh_230, kh_443, \
                         kh_448, kh_450, kh_457, kh_459, kh_461 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_159[k] = f_68 * kh_2[k]
                   + f_69 * kh_7[k]
                   - f_70 * kh_9[k]
                   + f_68 * kh_16[k]
                   - f_70 * kh_18[k]
                   + f_71 * kh_20[k]
                   - f_64 * kh_65[k]
                   - f_65 * kh_70[k]
                   + f_66 * kh_72[k]
                   - f_64 * kh_79[k]
                   + f_66 * kh_81[k]
                   - f_67 * kh_83[k]
                   + f_61 * kh_212[k]
                   + f_62 * kh_217[k]
                   - f_63 * kh_219[k]
                   + f_61 * kh_226[k]
                   - f_63 * kh_228[k]
                   + f_59 * kh_230[k]
                   - f_57 * kh_443[k]
                   - f_58 * kh_448[k]
                   + f_59 * kh_450[k]
                   - f_57 * kh_457[k]
                   + f_59 * kh_459[k]
                   - f_60 * kh_461[k];
    }

#pragma omp simd aligned(kh_0, kh_3, kh_5, kh_10, kh_12, kh_14, kh_63, kh_66, kh_68, kh_73, \
                         kh_75, kh_77, kh_210, kh_213, kh_215, kh_220, kh_222, kh_224, kh_441, \
                         kh_444, kh_446, kh_451, kh_453, kh_455 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_160[k] = f_53 * kh_0[k]
                   + f_54 * kh_3[k]
                   - f_55 * kh_5[k]
                   + f_53 * kh_10[k]
                   - f_55 * kh_12[k]
                   + f_56 * kh_14[k]
                   - f_49 * kh_63[k]
                   - f_50 * kh_66[k]
                   + f_51 * kh_68[k]
                   - f_49 * kh_73[k]
                   + f_51 * kh_75[k]
                   - f_52 * kh_77[k]
                   + f_45 * kh_210[k]
                   + f_46 * kh_213[k]
                   - f_47 * kh_215[k]
                   + f_45 * kh_220[k]
                   - f_47 * kh_222[k]
                   + f_48 * kh_224[k]
                   - f_41 * kh_441[k]
                   - f_42 * kh_444[k]
                   + f_43 * kh_446[k]
                   - f_41 * kh_451[k]
                   + f_43 * kh_453[k]
                   - f_44 * kh_455[k];
    }

#pragma omp simd aligned(kh_2, kh_9, kh_16, kh_18, kh_65, kh_72, kh_79, kh_81, kh_212, kh_219, \
                         kh_226, kh_228, kh_443, kh_450, kh_457, \
                         kh_459 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_161[k] = -f_75 * kh_2[k]
                   + f_39 * kh_9[k]
                   + f_75 * kh_16[k]
                   - f_39 * kh_18[k]
                   + f_74 * kh_65[k]
                   - f_37 * kh_72[k]
                   - f_74 * kh_79[k]
                   + f_37 * kh_81[k]
                   - f_73 * kh_212[k]
                   + f_35 * kh_219[k]
                   + f_73 * kh_226[k]
                   - f_35 * kh_228[k]
                   + f_72 * kh_443[k]
                   - f_33 * kh_450[k]
                   - f_72 * kh_457[k]
                   + f_33 * kh_459[k];
    }

#pragma omp simd aligned(kh_0, kh_3, kh_5, kh_10, kh_12, kh_63, kh_66, kh_68, kh_73, kh_75, \
                         kh_210, kh_213, kh_215, kh_220, kh_222, kh_441, kh_444, kh_446, \
                         kh_451, kh_453 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_162[k] = -f_31 * kh_0[k]
                   + f_29 * kh_3[k]
                   + f_32 * kh_5[k]
                   + f_28 * kh_10[k]
                   - f_30 * kh_12[k]
                   + f_15 * kh_63[k]
                   - f_26 * kh_66[k]
                   - f_17 * kh_68[k]
                   - f_25 * kh_73[k]
                   + f_27 * kh_75[k]
                   - f_23 * kh_210[k]
                   + f_21 * kh_213[k]
                   + f_24 * kh_215[k]
                   + f_20 * kh_220[k]
                   - f_22 * kh_222[k]
                   + f_18 * kh_441[k]
                   - f_16 * kh_444[k]
                   - f_19 * kh_446[k]
                   - f_15 * kh_451[k]
                   + f_17 * kh_453[k];
    }

#pragma omp simd aligned(kh_2, kh_7, kh_16, kh_65, kh_70, kh_79, kh_212, kh_217, kh_226, \
                         kh_443, kh_448, kh_457 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_163[k] = f_82 * kh_2[k]
                   - f_83 * kh_7[k]
                   + f_82 * kh_16[k]
                   - f_80 * kh_65[k]
                   + f_81 * kh_70[k]
                   - f_80 * kh_79[k]
                   + f_78 * kh_212[k]
                   - f_79 * kh_217[k]
                   + f_78 * kh_226[k]
                   - f_76 * kh_443[k]
                   + f_77 * kh_448[k]
                   - f_76 * kh_457[k];
    }

#pragma omp simd aligned(kh_0, kh_3, kh_10, kh_63, kh_66, kh_73, kh_210, kh_213, kh_220, \
                         kh_441, kh_444, kh_451 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_164[k] = f_10 * kh_0[k]
                   - f_9 * kh_3[k]
                   + f_8 * kh_10[k]
                   - f_7 * kh_63[k]
                   + f_6 * kh_66[k]
                   - f_5 * kh_73[k]
                   + f_0 * kh_210[k]
                   - f_4 * kh_213[k]
                   + f_3 * kh_220[k]
                   - f_2 * kh_441[k]
                   + f_1 * kh_444[k]
                   - f_0 * kh_451[k];
    }
}

}  // namespace simdtrf
